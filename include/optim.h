/**
 * @file optim.h
 *
 * @brief Interface and declarations for the 3 optimization approaches.
 *
 * Exposes VHP, LP and FP optimizations to outside (e.g. Player class
 * can call them).
 *
 *  Created on: Jun 7, 2016
 *      Author: okoc
 */

#ifndef OPTIMPOLY_H_
#define OPTIMPOLY_H_

// optimization and math libraries
#include <thread>
#include <math.h>
#include <nlopt.h>
#include "string.h" //for bzero
#include "constants.h"

// defines
const int EQ_CONSTR_DIM = 3*NCART;
const int INEQ_CONSTR_DIM = 2*NDOF_ACTIVE + 2*NDOF; // both strike and returning trajectories, min and max
const double MAX_VEL = 10;
const double MAX_ACC = 200;

using namespace arma;

/**
 * @brief Desired racket/ball positions, vels and racket normals for dt*Nmax seconds.
 *
 * This is the structure passed to Optimization class Optim
 * and its descendants.
 * For FP, we compute desired racket positions, velocities and normals
 * based on predicted ball path inside robot workspace.
 * For DP, we only use the ball predicted positions and velocities.
 */
struct optim_des {
	mat racket_pos = zeros<mat>(NCART,1); //!< racket desired pos
	mat racket_vel = zeros<mat>(NCART,1); //!< racket desired vel
	mat racket_normal = zeros<mat>(NCART,1); //!< racket desired normals
	mat ball_pos = zeros<mat>(NCART,1); //!< incoming ball predicted pos.
	mat ball_vel = zeros<mat>(NCART,1); //!< incoming ball predicted vels.
	double dt = DT; //!< time step between each prediction
	int Nmax = 1; //!< max number of time steps for prediction
};

/**
 * @brief 3rd order spline parameters returned from optimizer.
 *
 * If optimization is still running, player does not update trajectories.
 * IF optimization was successful (update is TRUE) and
 * it has terminated (running is FALSE), then player class can update trajectories.
 */
struct spline_params {
	double time2hit = 1.0; //!< free-final time (for hitting ball)
	mat a = zeros<mat>(NDOF_ACTIVE,4); //!< strike poly params of 3rd order
	mat b = zeros<mat>(NDOF_ACTIVE,4); //!< return poly params of 3rd order
};

/**
 * @brief Desired/actual joint positions, velocities, accelerations.
 *
 * Main Player function play() updates desired joint values
 * every 2ms for Barrett WAM.
 */
struct joint {
	vec7 q = zeros<vec>(NDOF_ACTIVE); //!< desired joint pos
	vec7 qd = zeros<vec>(NDOF_ACTIVE); //!< desired joint vels
	vec7 qdd = zeros<vec>(NDOF_ACTIVE); //!< desired joint accs
};

/**
 * @brief Weights for optimization penalties used by DP only.
 *
 * These weights are NOT used in any interesting way so far.
 * It would be interesting to learn/update these weights based on
 * successful returns in real robot experiments.
 */
struct weights {
	double R_strike[NDOF_ACTIVE] = {0.0}; //!< acceleration weights for running cost
	double R_hit = 0.0; //!< weight of dist from racket centre to hit location
};


/**
 * @brief Base class for all optimizers.
 *
 * Containts base class methods, members and virtual methods
 * to be implemented.
 */
class Optim {

protected:
	static const int OPTIM_DIM = 2*NDOF_ACTIVE + 1; //!< dim. of optim problem
	bool lookup = false; //!< use lookup table methods to init. optim params.
	bool verbose = true; //!< verbose output printed
	bool moving = false; //!< robot is already moving so use last computed values to init.
	bool update = false; //!< optim finished and soln. seems valid
	bool running = false; //!< optim still RUNNING
	bool detach = false; //!< detach optim in another thread
	mat lookup_table; //!< lookup table used to init. optim values
	nlopt_opt opt; //!< optimizer from NLOPT library

	double qf[NDOF_ACTIVE] = {0.0}; //!< saved joint positions after optim
	double qfdot[NDOF_ACTIVE] = {0.0}; //!< saved joint velocities after optim
	double T = 1.0; //!< saved hitting time after optim terminates

	void init_last_soln(double *x) const;
	void init_rest_soln(double *x) const;
	double test_soln(const double *x) const;
	void finalize_soln(const double *x, const double dt);
	void optim_rest_posture(vec7 & q_rest_des);
	void update_rest_state(const vec7 & q_rest_new);
	void optim();
public:
	optim_des *param_des; //!< Desired racket and/or ball predicted vals.
	double lb[2*NDOF_ACTIVE+1]; //!< Joint lower limits, joint vel. lower limit and min. hitting time
	double ub[2*NDOF_ACTIVE+1]; //!< Joint upper limits, joint vel. upper limit and max. hitting time
	double qrest[NDOF_ACTIVE] = {0.0}; //!< FIXED Resting posture for optimizers to compute return traj.
	double q0[NDOF_ACTIVE] = {0.0}; //!< Initial joint state needed to compute traj. acc.
	double q0dot[NDOF_ACTIVE] = {0.0}; //!< Initial joint velocities needed to compute traj. acc.
	double time2return = 1.0; //!< FIXED Time to return
	Optim();
	Optim(const vec7 & qrest_, double lb_[], double ub_[]);
	bool check_update();
	bool check_running();
	void set_moving(bool flag);
	void set_detach(bool flag);
	void set_verbose(bool flag);
	bool get_params(const joint & qact, spline_params & p);
	void run_qrest_optim(vec7 & q_rest_des);
	void update_init_state(const joint & qact);
	void set_return_time(const double & time);
	void set_des_params(optim_des *params);
	void run();
};


// functions that all players use
void joint_limits_ineq_constr(unsigned m, double *result,
		                      unsigned n, const double *x, double *grad, void *data);

void calc_strike_poly_coeff(const double *q0, const double *q0dot, const double *x,
		                    double *a1, double *a2);
void calc_return_poly_coeff(const double *q0, const double *q0dot,
		                    const double *x, const double time2return,
		                    double *a1, double *a2);
void calc_strike_extrema_cand(const double *a1, const double *a2, const double T,
		                      const double *q0, const double *q0dot,
							  double *joint_max_cand, double *joint_min_cand);
void calc_return_extrema_cand(const double *a1, const double *a2,
		                      const double *x, const double time2return,
							  double *joint_max_cand, double *joint_min_cand);

#endif /* OPTIMPOLY_H_ */
