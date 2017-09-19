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
 * @brief Player class and the functions it uses are stored here.
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
#include "kinematics.h"
#include "math.h"

using namespace arma;

/**
 * @brief Initialize Player.
 *
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

	opt = new Optim(q0, flags.active_dofs);
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
			init_ball_state = true;
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
			//cout << "OBS:" << obs.t() << "STATE:" << x.t();
			//<< "VAR:" << P.diag().t() << endl;
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
 * @brief Play the game.
 *
 * Main function for playing a game. Calls one of the optimization
 * algorithms (depending on initialization) and
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
 * @brief Cheat by getting the exact ball state in simulation.
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

	const double time_pred = 2.0;
	mat balls_pred;

	// if ball is fast enough and robot is not moving consider optimization
	if (check_update(qact)) {
		predict_ball(time_pred,balls_pred,filter);
		pred_params.Nmax = (int)(time_pred/DT);
		pred_params.ball_pos = balls_pred.rows(X,Z);
		pred_params.ball_vel = balls_pred.rows(DX,DZ);
		opt->set_des_params(&pred_params);
		opt->update_init_state(qact);
		opt->set_verbose(true);
		opt->run();
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
		feasible = (state_est(DY) >= 0.0);
		update = !opt->check_update() && !opt->check_running();
		// ball is incoming
		if (pflags.mpc) {// && t_poly > 0.0) {
			activate = (!pflags.detach) ? counter % 5 == 0 :
					                        timer.toc() > (1.0/pflags.freq_mpc);
			incoming = state_est(Y) > state_last(Y);
			update = update && valid_obs && activate && feasible && incoming;
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

/**
 * @brief Generate strike and return traj. incrementally
 *
 * Given polynomial parameters saved in poly,
 * move on to the NEXT desired state only (joint pos,vel,acc).
 * @param poly Polynomial parameters updated in OPTIM classes
 * @param q_rest_des FIXED desired resting posture
 * @param time2return FIXED time to return to rest posture after hit
 * @param t The time passed already following trajectory
 * @param qdes Update pos,vel,acc values of this desired joints structure
 * @return
 */
bool update_next_state(const spline_params & poly,
		           const vec7 & q_rest_des,
				   const double time2return,
				   double & t,
				   joint & qdes) {
	mat a,b;
	double tbar;
	bool flag = true;

	if (t <= poly.time2hit) {
		a = poly.a;
		qdes.q = a.col(0)*t*t*t + a.col(1)*t*t + a.col(2)*t + a.col(3);
		qdes.qd = 3*a.col(0)*t*t + 2*a.col(1)*t + a.col(2);
		qdes.qdd = 6*a.col(0)*t + 2*a.col(1);
		t += DT;
		//cout << qdes.q << qdes.qd << qdes.qdd << endl;
	}
	else if (t <= poly.time2hit + time2return) {
		b = poly.b;
		tbar = t - poly.time2hit;
		qdes.q = b.col(0)*tbar*tbar*tbar + b.col(1)*tbar*tbar + b.col(2)*tbar + b.col(3);
		qdes.qd = 3*b.col(0)*tbar*tbar + 2*b.col(1)*tbar + b.col(2);
		qdes.qdd = 6*b.col(0)*tbar + 2*b.col(1);
		t += DT;
	}
	else {
		//printf("Hitting finished!\n");
		t = 0.0;
		flag = false;
		qdes.q = q_rest_des;
		qdes.qd = zeros<vec>(NDOF_ACTIVE);
		qdes.qdd = zeros<vec>(NDOF_ACTIVE);
	}
	return flag;
}

/**
 * @brief Initialize an EKF
 *
 * Called generally when ball state needs to be reset
 * Useful for passing to Player constructor.
 * @param var_model Process noise multiplier for identity matrix.
 * @param var_noise Observation noise mult. for identity matrix.
 * @return EKF Extended Kalman Filter (state uninitialized!)
 */
EKF init_filter(const double var_model, const double var_noise) {

	mat C = eye<mat>(3,6);
	mat66 Q = var_model * eye<mat>(6,6);
	mat33 R = var_noise * eye<mat>(3,3);
	EKF filter = EKF(calc_next_ball,C,Q,R);
	return filter;
}

/**
 * @brief Function that integrates a basketball for an outside filter.
 *
 * Function exposes the integration to filters, e.g. an EKF.
 * They can use then to apply predict() using the this function pointer.
 *
 * The basketball is attached to a string so first finding the angle responsible
 * for the observation, integrating it and then finding the next ball state
 *
 */
vec calc_next_ball(const vec & xnow, const double dt, const void *fp) {

	const vec3 base_pendulum = {0.0, 0.9, 1.0};
	const double string_len = 1.0;
	const double basketball_radius = 0.1213;
	const double gravity = -9.8;
	const double friction = 0.03;

	vec3 ball_pos = xnow.head(3);
	vec3 ball_vel = xnow.tail(3);

	ball_pos -= base_pendulum;

	// to integrate the ball first extract the theta and theta dot
	// integrate them and then revert back to current ball pos
	double theta = atan(ball_pos(Y) / ball_pos(Z));
	double theta_dot = -ball_vel(Y) / ((string_len + basketball_radius) * cos(theta));

	theta_dot += dt * (gravity * sin(theta) - friction * theta_dot);
	theta += dt * theta_dot;

	// we're in the third quadrant
	ball_pos(Y) = base_pendulum(Y) - (string_len + basketball_radius) * sin(theta);
	ball_pos(Z) = base_pendulum(Z) - (string_len + basketball_radius) * cos(theta);
	ball_vel(Y) = -(string_len + basketball_radius) * theta_dot * cos(theta);
	ball_vel(Z) = (string_len + basketball_radius) * theta_dot * sin(theta);
	return join_vert(ball_pos,ball_vel);
}


/**
 * @brief Checks to see if the observation is new (updated)
 *
 * The blobs need to be at least tol apart from each other in distance
 *
 */
bool check_new_obs(const vec3 & obs, double tol) {

	static vec3 last_obs = zeros<vec>(3);

	if (norm(obs - last_obs) > tol) {
		last_obs = obs;
		return true;
	}
	return false;
}

/**
 * @brief Check to see if we want to reset the filter.
 *
 * Basically if a new ball appears 300 ms later than the last new ball
 * we reset the filter.
 *
 */
bool check_reset_filter(const bool newball, const int verbose, const double threshold) {

	bool reset = false;
	static int reset_cnt = 0;
	static bool firsttime = true;
	static wall_clock timer;

	if (firsttime) {
		firsttime = false;
		timer.tic();
	}

	if (newball) {
		if (timer.toc() > threshold) {
			reset = true;
			if (verbose > 0) {
				std::cout << "Resetting filter! Count: " << ++reset_cnt << std::endl;
			}
		}
		timer.tic();
	}
	return reset;
}
