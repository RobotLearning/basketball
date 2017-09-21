/**
 * @file player.hpp
 *
 * @brief Main class for playing Table Tennis declared here.
 *
 *  Created on: Feb 9, 2017
 *      Author: okoc
 */

#ifndef PLAYER_HPP_
#define PLAYER_HPP_

#include "kalman.h"
#include "optim.h"
#include "kinematics.h"

using namespace arma;

/**
 * Finite State machine for Table Tennis
 */
enum game { //trial state
	AWAITING,//!< AWAITING
	HIT,     //!< HIT
};

/**
 * @brief Options passed to Player class (algorithm, saving, corrections, etc.).
 */
struct player_flags {
	bool detach = false; //!< detach optimizations in another thread
	bool outlier_detection = false; //!< outlier detection for real robot
	bool mpc = false; //!< turn on/off corrections
	bool reset = true; //!< reinitializing player class
	bool save = false; //!< saving ball/robot data
	int verbosity = 0; //!< OFF, LOW, HIGH, ALL
	int freq_mpc = 1; //!< frequency of mpc updates if turned on
	int min_obs = 5; //!< number of observations to initialize filter
	double time2return = 0.5; //!< time to return to starting posture after hit
	double var_noise = 0.001; //!< variance of noise process (R)
	double var_model = 0.001; //!< variance of process noise (Q)
	double t_reset_thresh = 0.3; //!< resetting Kalman filter after this many seconds pass without getting valid obs.
	ivec active_dofs = join_vert(LEFT_ARM,RIGHT_ARM);
};

/**
 *
 * @brief Table Tennis Player class for playing Table Tennis.
 *
 * The methods play() or cheat() must be called every DT milliseconds.
 */
class Player {

private:

	// data fields
	bool init_ball_state = false;
	EKF & filter; // filter for the ball estimation
	vec7 q_rest_des; // desired resting joint state
	double t_obs = 0.0; // counting time stamps for resetting filter
	double t_poly = 0.0; // time passed on the hitting spline
	bool valid_obs = true; // ball observed is valid (new ball and not an outlier)
	int num_obs = 0; // number of observations received
	game game_state = AWAITING;
	player_flags pflags;
	optim_des pred_params;
	mat observations; // for initializing filter
	mat times; // for initializing filter
	spline_params poly_left;
	spline_params poly_right;
	Optim *opt_left; // optimizers
	Optim *opt_right;

	// ball estimation
	void estimate_ball_state(const vec3 & obs);

	// optimization for different players
	void optim_param(const joint & qact); // run optimizer for Focused player

	void calc_opt_params(const joint & qact);
	bool check_update(const joint & qact) const; // flag for (re)running optimization
	void calc_next_state(const joint & qact, joint & qdes);

public:

	Player(const vec & q0, EKF & filter, player_flags & flags);
	~Player();

	// auxiliary function, public interface for filter test performance
	vec6 filt_ball_state(const vec3 & obs);
	bool filter_is_initialized() const ;
	void reset_filter(double std_model, double std_noise);

	// main function
	void play(const joint & qact, const vec3 & ball_obs, joint & qdes);

	// cheat function for simulation (with exact knowledge of ball state)
	void cheat(const joint & qact, const vec6 & ballstate, joint & qdes);

};

// ball estimation and filter constructor/state initialization
EKF init_filter(const double var_model = 0.001, const double var_noise = 0.001);
void estimate_ball_linear(const mat & observations,
		                  const vec & times,
					      const bool verbose,
					      vec6 & init_est);
bool check_new_obs(const vec3 & obs, double tol);
void predict_ball(const double & time_pred, mat & balls_pred, EKF & filter);
vec calc_next_ball(const vec & xnow, const double dt, const void *fp);
bool check_reset_filter(const bool newball, const int verbose, const double threshold);
bool check_new_obs(const vec3 & obs, double tol);

// movement generation
bool update_next_state(const vec & q_rest_des,
				   const double time2return,
				   const bool right_arm,
				   spline_params & poly,
				   joint & qdes);


#endif /* PLAYER_HPP_ */
