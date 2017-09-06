/*! \mainpage Optimal Control based Trajectory Generation
 *
 * \section intro_sec Introduction
 *
 * Player class is the orchestrator of ball hitting movements.
 * These are 3rd order striking and returning polynomials computed
 * differently for each method.
 *
 *
 */


/**
 * @file player.cpp
 *
 * @brief Table Tennis player class and the functions it uses are stored here.
 *
 * Player class launches 3 different optimization algorithms
 * for striking trajectory generation.
 *
 *  Created on: Feb 8, 2017
 *  Author: okoc
 */

#include <armadillo>
#include <thread>
#include <string>
#include "stdlib.h"
#include "player.hpp"
#include "constants.h"
#include "kalman.h"
#include "optim.h"

using namespace arma;

/**
 * @brief Initialize Table Tennis Player.
 *
 * Table Tennis Player class that can run 3 different trajectory
 * generation algorithms.
 * VHP and FP try to return the ball to the centre of the opponents court,
 * LP tries to just return the ball to the opponents court.
 *
 * @param q0 Initial joint positions.
 * @param filter_ Reference to an input filter (must be initialized in a separate line).
 * @param flags Flags/options for player class, initialized with c++11 (see player.hpp)
 *              or through player.cfg file (see sl_interface)
 */
Player::Player(const vec7 & q0, EKF & filter_, player_flags & flags)
               : filter(filter_), pflags(flags) {

	q_rest_des = q0;
	observations = zeros<mat>(3,pflags.min_obs); // for initializing filter
	times = zeros<vec>(pflags.min_obs); // for initializing filter
	//load_lookup_table(lookup_table);

	double lb[2*NDOF_ACTIVE+1];
	double ub[2*NDOF_ACTIVE+1];
	double SLACK = 0.02;
	double Tmax = 1.0;
	set_bounds(lb,ub,SLACK,Tmax);

	opt = new Optim(q0,lb,ub);
	opt->set_return_time(pflags.time2return);
	opt->set_verbose(pflags.verbosity > 1);
	opt->set_detach(pflags.detach);
}

/*
 *
 * Deconstructor for the Player class.
 * Frees the pointer to Optimization classes
 *
 */
Player::~Player() {

	delete opt;
}

/**
 * If filter is initialized returns true
 * @return
 */
bool Player::filter_is_initialized() const {
	return init_ball_state;
}

/*
 * Filter the blob information with a Kalman Filter.
 * (Extended) KF is used both in simulation mode and for real robot.
 *
 * Checking for new ball that is at least 1 mm away from last observation
 * Checking for also outliers.
 * Resets if the ball suddenly appears on the opponent's court.
 *
 * Ball is valid if ball is a new ball and (in real robot mode)
 * it is not an outlier!
 *
 * Note: we're assuming that time elasped dt = DT = 0.002 seconds every time!
 *
 */
void Player::estimate_ball_state(const vec3 & obs) {

	int verb = pflags.verbosity;
	bool newball = check_new_obs(obs,1e-3);
	valid_obs = false;

	if (check_reset_filter(newball,verb,pflags.t_reset_thresh)) {
		filter = init_filter(pflags.var_model,pflags.var_noise);
		num_obs = 0;
		init_ball_state = false;
		game_state = AWAITING;
		//cout << obs << endl;
		t_obs = 0.0; // t_cumulative
	}

	if (num_obs < pflags.min_obs && newball) {
		times(num_obs) = t_obs;
		observations.col(num_obs) = obs;
		num_obs++;
		if (num_obs == pflags.min_obs) {
			if (verb >= 1)
				cout << "Estimating initial ball state\n";
			vec6 x;
			mat P;
			P.eye(6,6);
			estimate_ball_linear(observations,times,pflags.verbosity > 2,x);
			filter.set_prior(x,P);
			//cout << OBS << TIMES << filter.get_mean() << endl;
		}

	}
	else if (init_ball_state) { // comes here if there are enough balls to start filter
		filter.predict(DT,true); //true);
		if (newball) {
			valid_obs = true;
			if (pflags.outlier_detection)
				valid_obs = !filter.check_outlier(obs,verb > 2);
		}
		if (valid_obs) {
			filter.update(obs);
			//vec x = filter.get_mean();
			//mat P = (filter.get_covar());
			//cout << "OBS:" << obs.t() << "STATE:" << x.t() << "VAR:" << P.diag().t() << endl;
		}

	}
	t_obs += DT;
}

/**
 *
 * @brief Public interface for estimating ball state.
 *
 * This interface allows us to test/debug ball state estimation
 * (which is private).
 *
 * @param obs Ball position observation as a 3-vector.
 * @return Ball state as a 6-vector, if filter is not initialized,
 * returns the observation as positions and zeroes as velocities.
 */
vec6 Player::filt_ball_state(const vec3 & obs) {

	estimate_ball_state(obs);
	try {
		return filter.get_mean();
	}
	catch (const std::exception & exception) {
		return join_vert(obs,zeros<vec>(3));
	}
}

/**
 * @brief Play Table Tennis.
 *
 * Main function for playing Table Tennis. Calls one of three different
 * trajectory generation algorithms (depending on initialization) and
 * updates the desired joint states when the optimization threads have finished.
 *
 * @param qact Actual joint positions, velocities, accelerations.
 * @param ball_obs Ball observations (positions as 3-vector).
 * @param qdes Desired joint positions, velocities, accelerations.
 */
void Player::play(const joint & qact,const vec3 & ball_obs, joint & qdes) {

	estimate_ball_state(ball_obs);

	optim_param(qact);

	// generate movement or calculate next desired step
	calc_next_state(qact, qdes);

}


/**
 * @brief Cheat Table Tennis by getting the exact ball state in simulation.
 *
 * Similar to play() method, this method receives from the simulator the
 * exact ball states, which bypasses then the ball estimation method.
 * Useful for debugging the internal filter used.
 *
 * @param qact Actual joint positions, velocities, accelerations.
 * @param ballstate Ball state (positions AND velocities).
 * @param qdes Desired joint positions, velocities, accelerations.
 */
void Player::cheat(const joint & qact, const vec6 & ballstate, joint & qdes) {

	// resetting legal ball detecting to AWAITING state
	if (ballstate(Y) > 2.0)
		game_state = AWAITING;
	filter.set_prior(ballstate,0.01*eye<mat>(6,6));

	optim_param(qact);

	// generate movement or calculate next desired step
	calc_next_state(qact, qdes);
}

/*
 *
 * Calculate the optimization parameters using an NLOPT nonlinear optimization algorithm
 * in another thread
 *
 * The optimized parameters are: qf, qf_dot, T
 * assuming T_return and q0 are fixed
 *
 */
void Player::optim_param(const joint & qact) {

	mat balls_pred;

	// if ball is fast enough and robot is not moving consider optimization
	if (check_update(qact)) {
		predict_ball(2.0,balls_pred,filter);
		opt->set_des_params(&pred_params);
		opt->update_init_state(qact);
		opt->run();
		}
		else {
			//cout << "Ball is not legal!\n";
		}
	}
}

/*
 * Check MPC flag and update if possible
 *
 * IF MPC IS TURNED OFF
 * if ball is incoming and robot is not moving consider optimization
 *
 * IF MPC IS TURNED ON
 * then additionally consider (after running initial optimization)
 * relaunching optimization if ball is valid (new ball and not an outlier)
 * the frequency of updates is respected, and the ball has not passed the y-limit
 *
 */
bool Player::check_update(const joint & qact) const {

	static int firsttime = true;
	static int counter;
	static vec6 state_last = zeros<vec>(6);
	static wall_clock timer;
	vec6 state_est;
	bool update = false;

	bool activate, passed_lim, incoming, feasible = false;

	if (firsttime) {
		timer.tic();
		firsttime = false;
	}

	try {
		state_est = filter.get_mean();
		counter++;
		feasible = (state_est(DY) > 0.5);
		update = !opt->check_update() && !opt->check_running();
		// ball is incoming
		if (pflags.mpc) {// && t_poly > 0.0) {
			activate = (!pflags.detach) ? counter % 5 == 0 :
					                        timer.toc() > (1.0/pflags.freq_mpc);
			passed_lim = state_est(Y) > pos(Y);
			incoming = state_est(Y) > state_last(Y);
			update = update && valid_obs && activate && feasible && !passed_lim && incoming;
		}
		else {
			update = update && (t_poly == 0.0) && feasible; // only once
		}
		state_last = state_est;
		if (update) {
			//cout << "Consider reacting to ball: " << state_est.t() << endl;
			//cout << num_updates++ << endl;
			timer.tic();
		}
	}
	catch (const std::exception & not_init_error) {
		update = false;
	}

	return update;
}

/*
 * Unfold the next desired state of the 3rd order polynomials in joint space
 * If movement finishes then the desired state velocities and accelerations are zeroed.
 *
 * Multithreading : if after initial lookup, the computation in the other thread terminates, then we
 * synchonize the values between the threads and compute the next desired states given these new polynomial
 * parameters (qf, qf_dot and T_hit)
 *
 */
void Player::calc_next_state(const joint & qact, joint & qdes) {

	// this should be only for MPC?
	if (opt->get_params(qact,poly)) {
		if (pflags.verbosity) {
			std::cout << "Launching/updating strike" << std::endl;
		}
		t_poly = DT;
		opt->set_moving(true);
	}

	// make sure we update after optim finished
	if (t_poly > 0.0) {
		if (!update_next_state(poly,q_rest_des,pflags.time2return,t_poly,qdes)) {
			opt->set_moving(false);
		}
	}

}

/**
 * @brief Method useful for testing performance of different players.
 *
 * Using many balls in simulation requires fast resetting
 * Setting a time threshold as a resetting condition won't work in this case.
 *
 */
void Player::reset_filter(double var_model, double var_noise) {

	filter = init_filter(var_model,var_noise);
	init_ball_state = false;
	num_obs = 0;
	game_state = AWAITING;
	t_obs = 0.0;
}

/**
 * @brief Predict ball with the models fed into the filter
 *
 * Number of prediction steps is given by Nmax in racket
 * parameters
 *
 */
void predict_ball(const double & time_pred, mat & balls_pred, EKF & filter) {

	//static wall_clock timer;
	//timer.tic();
	int N = (int)(time_pred/DT);
	balls_pred = filter.predict_path(DT,N);
	//cout << "Pred. ball time: " << 1000 * timer.toc() << " ms." << endl;
}

/**
 * @brief Least squares to estimate prior given
 * matrix of observations Y of column length N, each row
 * corresponding to one
 *
 * If observations arrive as column vectors then we take
 * transpose of it.
 *
 */
void estimate_ball_linear(const mat & observations,
		                  const vec & times,
					      const bool verbose,
					      vec6 & init_est) {

	int num_samples = times.n_elem;
	mat M = zeros<mat>(num_samples,3);

	// and create the data matrix
	for (int i = 0; i < num_samples; i++) {
		M(i,0) = 1.0;
		M(i,1) = times(i);
		M(i,2) = times(i) * times(i);
	}
	// solving for the parameters
	mat Beta = solve(M,observations.t());
	//mat Beta = pinv(M,0.01) * observations.t();
	init_est = join_horiz(Beta.row(0),Beta.row(1)).t();

	if (verbose) {
		cout << "Times:" << times.t() << endl;
		cout << "Data:\n" << observations.t() << endl;
		cout << "Initial est:" << init_est.t() << endl;
	}
}
