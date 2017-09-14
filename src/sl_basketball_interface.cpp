
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <armadillo>
#include <cmath>
#include <sys/time.h>
#include "constants.h"
#include "kalman.h"
#include "player.hpp"
#include "kinematics.h"

using namespace arma;


/* The data structures from SL */
/**
 * @brief (actual) joint space state for each DOF
 */
struct SL_Jstate {
	double   th;   /*!< theta */
	double   thd;  /*!< theta-dot */
	double   thdd; /*!< theta-dot-dot */
	double   ufb;  /*!< feedback portion of command */
	double   u;    /*!< torque command */
	double   load; /*!< sensed torque */
};

/**
 * @brief (desired) joint space state commands for each DOF
 */
struct SL_DJstate { /*!< desired values for controller */
	double   th;   /*!< theta */
	double   thd;  /*!< theta-dot */
	double   thdd; /*!< theta-dot-dot */
	double   uff;  /*!< feedforward torque command */
	double   uex;  /*!< externally imposed torque */
};

/**
 * @brief (actual) Cartesian state
 */
struct SL_Cstate {
	double   x[NCART+1];    /*!< Position [x,y,z] */
	double   xd[NCART+1];   /*!< Velocity */
	double   xdd[NCART+1];  /*!< Acceleration */
};

/**
 * @brief Vision blob info coming from SL (after calibration).
 *
 */
struct blob_state {
	int status; //!< was ball detected reliably in cameras
	double pos[NCART]; //!< ball center cartesian positions from cameras 1 and 2(after calibration)
};

player_flags opt; //!< global structure for setting Player options

#include "sl_basketball_interface.h"

void load_options() {
	//TODO:
}

/**
 * @brief Interface to the PLAYER class that generates desired hitting trajectories.
 *
 * First initializes the player according to the pre-set options
 * and then starts calling play() interface function. Must be called every DT ms.
 *
 *
 * @param joint_state Actual joint positions, velocities, accelerations.
 * @param blobs 3d-positions from cameras stored in blobs[1..NBLOBS]
 * @param joint_des_state Desired joint position, velocity and acceleration commands.
 */
void play(const SL_Jstate joint_state[NDOF+1],
		  const blob_state *blobs,
		  SL_DJstate joint_des_state[NDOF+1]) {

	static vec7 q0;
	static vec3 ball_obs;
	static joint qact;
	static joint qdes;
	static Player *robot = nullptr; // pointer to player
	static EKF filter = init_filter(0.001,0.001);

	if (opt.reset) {
		for (int i = 0; i < NDOF_ACTIVE; i++) {
			qdes.q(i) = q0(i) = joint_state[opt.active_dofs(i)+1].th;
			qdes.qd(i) = 0.0;
			qdes.qdd(i) = 0.0;
		}
		filter = init_filter(0.3,0.001);
		delete robot;
		opt.detach = true;
		robot = new Player(q0,filter,opt);
		opt.reset = false;
	}
	else {
		for (int i = 0; i < NDOF_ACTIVE; i++) {
			qact.q(i) = joint_state[opt.active_dofs(i)+1].th;
			qact.qd(i) = joint_state[opt.active_dofs(i)+1].thd;
			qact.qdd(i) = joint_state[opt.active_dofs(i)+1].thdd;
		}
		fuse_blobs(blobs,ball_obs);
		robot->play(qact,ball_obs,qdes);
	}

	// update desired joint state
	for (int i = 0; i < NDOF_ACTIVE; i++) {
		joint_des_state[opt.active_dofs(i)+1].th = qdes.q(i);
		joint_des_state[opt.active_dofs(i)+1].thd = qdes.qd(i);
		joint_des_state[opt.active_dofs(i)+1].thdd = qdes.qdd(i);
	}
}

/**
 * @brief  CHEAT with exact knowledge of ball state.
 *
 * Interface to the PLAYER class that generates desired hitting trajectories.
 * First initializes the player and then starts calling cheat() interface function.
 *
 * @param joint_state Actual joint positions, velocities, accelerations.
 * @param sim_ball_state Exact simulated ball state (positions and velocities).
 * @param joint_des_state Desired joint position, velocity and acceleration commands.
 */
void cheat(const SL_Jstate joint_state[NDOF+1],
		  const SL_Cstate *sim_ball_state,
		  SL_DJstate joint_des_state[NDOF+1]) {

	static vec7 q0;
	static vec6 ball_state;
	static joint qact;
	static joint qdes;
	static Player *robot = nullptr; // centered player
	static EKF filter = init_filter(0.001,0.001);

	if (opt.reset) {
		for (int i = 0; i < NDOF_ACTIVE; i++) {
			qact.q(i) = joint_state[opt.active_dofs(i)+1].th;
			qact.qd(i) = joint_state[opt.active_dofs(i)+1].thd;
			qact.qdd(i) = joint_state[opt.active_dofs(i)+1].thdd;
		}
		delete robot;
		opt.detach = true;
		robot = new Player(q0,filter,opt);
		opt.reset = false;
	}
	else {
		for (int i = 0; i < NDOF_ACTIVE; i++) {
			qact.q(i) = joint_state[opt.active_dofs(i)+1].th;
			qact.qd(i) = joint_state[opt.active_dofs(i)+1].thd;
			qact.qdd(i) = joint_state[opt.active_dofs(i)+1].thdd;
		}
		for (int i = 0; i < NCART; i++) {
			ball_state(i) = sim_ball_state->x[i+1];
			ball_state(i+NCART) = sim_ball_state->xd[i+1];
		}
		robot->cheat(qact,ball_state,qdes);
	}

	// update desired joint state
	for (int i = 0; i < NDOF_ACTIVE; i++) {
		joint_des_state[opt.active_dofs(i)+1].th = qdes.q(i);
		joint_des_state[opt.active_dofs(i)+1].thd = qdes.qd(i);
		joint_des_state[opt.active_dofs(i)+1].thdd = qdes.qdd(i);
	}
}

/*
 *
 * Fusing multiple blobs
 * TODO:
 *
 */
static void fuse_blobs(const blob_state * blobs, vec3 & obs) {

	if (blobs->status) {
		obs(0) = blobs->pos[0];
		obs(1) = blobs->pos[1];
		obs(2) = blobs->pos[2];
	}
}