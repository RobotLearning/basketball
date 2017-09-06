
#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <armadillo>
#include <cmath>
#include <sys/time.h>
#include "kalman.h"
#include "player.hpp"

using namespace arma;

/*! links of the robot */
enum RobotDOFs {
  BASE=0,

  L_SFE, // 01
  L_SAA, // 02
  L_HR,  // 03
  L_EB,  // 04
  L_WR,  // 05
  L_WFE, // 06
  L_WAA, // 07

  R_SFE, // 08
  R_SAA, // 09
  R_HR,  // 10
  R_EB,  // 11
  R_WR,  // 12
  R_WFE, // 13
  R_WAA, // 14

  L_HFE, // 15
  L_HAA, // 16
  L_HFR, // 17
  L_KFE, // 18
  L_AR,  // 19
  L_AFE, // 20
  L_AAA, // 21

  R_HFE, // 22
  R_HAA, // 23
  R_HFR, // 24
  R_KFE, // 25
  R_AR,  // 26
  R_AFE, // 27
  R_AAA, // 28

  B_TR,  // 29
  B_TAA, // 30
  B_TFE, // 31

  B_HN,  // 32
  B_HT,  // 33
  B_HR,  // 34

  R_EP,  // 35
  R_ET,  // 36
  L_EP,  // 37
  L_ET,  // 38

  N_ROBOT_DOFS
};


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

struct player_flags {
	vec active_dofs = {R_SFE, R_SAA, R_HR, R_EB, R_WR, R_WFE, R_WAA};
	bool reset = false;
};

player_flags opt; //!< global structure for setting Player options

#include "sl_interface.h"

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
		  const blob_state blobs[NBLOBS],
		  SL_DJstate joint_des_state[NDOF+1]) {

	static vec7 q0;
	static vec3 ball_obs;
	static joint qact;
	static joint qdes;
	static Player *robot = nullptr; // centered player
	static EKF filter = init_filter(0.3,0.001);

	if (opt.reset) {
		for (int i = 0; i < NDOF_ACTIVE; i++) {
			qdes.q(i) = q0(i) = joint_state[opt.active_dofs(i)].th;
			qdes.qd(i) = 0.0;
			qdes.qdd(i) = 0.0;
		}
		filter = init_filter(0.3,0.001);
		delete robot;
		robot = new Player(q0,filter,opt);
		opt.reset = false;
	}
	else {
		for (int i = 0; i < NDOF_ACTIVE; i++) {
			qact.q(i) = joint_state[opt.active_dofs(i)].th;
			qact.qd(i) = joint_state[opt.active_dofs(i)].thd;
			qact.qdd(i) = joint_state[opt.active_dofs(i)].thdd;
		}
		fuse_blobs(blobs,ball_obs);
		robot->play(qact,ball_obs,qdes);
	}

	// update desired joint state
	for (int i = 0; i < NDOF; i++) {
		joint_des_state[i+1].th = qdes.q(i);
		joint_des_state[i+1].thd = qdes.qd(i);
		joint_des_state[i+1].thdd = qdes.qdd(i);
	}

}
