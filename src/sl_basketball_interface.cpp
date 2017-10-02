
#include <boost/program_options.hpp>
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

player_flags options; //!< global structure for setting Player options

#include "sl_basketball_interface.h"

/**
 * @brief Set optimization type for BASKETBALL movements.
 *
 * @param alg_num Select between three optimization types: LEFT_HAND, RIGHT_HAND, and BOTH_HANDS
 */
static void set_optim_type(const int opt_num) {

	switch (opt_num) {
		case 0:
			std::cout << "Optimizing only LEFT HAND..." << std::endl;
			options.optim_type = LEFT_HAND_OPT;
			break;
		case 1:
			std::cout << "Optimizing only RIGHT HAND..." << std::endl;
			options.optim_type = RIGHT_HAND_OPT;
			break;
		case 2:
			std::cout << "Optimizing BOTH HANDS..." << std::endl;
			options.optim_type = BOTH_HAND_OPT;
			break;
		default:
			options.optim_type = RIGHT_HAND_OPT;
	}
}

/**
 * @brief Set algorithm and options to initialize Player with.
 *
 * The global variable flags is set here and
 * the play() function will use it to initialize the Player class.
 *
 */
void load_options() {

	namespace po = boost::program_options;
	using namespace std;

	options.reset = true;
	string home = std::getenv("HOME");
	string config_file = home + "/basketball/" + "player.cfg";
	int optim_type;

    try {
		// Declare a group of options that will be
		// allowed in config file
		po::options_description config("Configuration");
		config.add_options()
		    ("outlier_detection", po::value<bool>(&options.outlier_detection)->default_value(true),
			      "OUTLIER DETECTION FOR REAL ROBOT!")
			("hand", po::value<int>(&optim_type)->default_value(0),
				  "optimization type: LEFT_HAND = 0, RIGHT_HAND = 1, BOTH = 2")
			("touch", po::value<bool>(&options.touch)->default_value(true),
				  "only touch the ball if TRUE, hit the ball with des. velocity if FALSE")
			("mpc", po::value<bool>(&options.mpc)->default_value(false),
				 "corrections (MPC)")
			("verbose", po::value<int>(&options.verbosity)->default_value(1),
		         "verbosity level")
		    ("save_data", po::value<bool>(&options.save)->default_value(false),
		         "saving robot/ball data")
			("start_optim_offset", po::value<double>(&options.offset),
				 "start optim offset in y-direction")
			("time2return", po::value<double>(&options.time2return),
						 "time to return to start posture")
			("freq_mpc", po::value<int>(&options.freq_mpc), "frequency of updates")
		    ("min_obs", po::value<int>(&options.min_obs), "minimum obs to start filter")
		    ("var_noise", po::value<double>(&options.var_noise), "std of filter obs noise")
		    ("var_model", po::value<double>(&options.var_model), "std of filter process noise")
		    ("t_reset_threshold", po::value<double>(&options.t_reset_thresh), "filter reset threshold time");
        po::variables_map vm;
        ifstream ifs(config_file.c_str());
        if (!ifs) {
            cout << "can not open config file: " << config_file << "\n";
        }
        else {
            po::store(parse_config_file(ifs, config), vm);
            notify(vm);
        }
    }
    catch(exception& e) {
        cout << e.what() << "\n";
    }
    set_optim_type(optim_type);
    options.detach = true; // always detached in SL/REAL ROBOT!
}

/**
 * @brief Integrate ball state in simulation using symplectic Euler.
 *
 * Using an internal pendulum model to predict ball state.
 * Derives ball state as a function of pendulum angle state.
 * Expose environment variables and constants to SL as an array.
 *
 */
void integrate_ball_state(const double dt, const SL_Cstate robot_state[NENDEFF+1], double ball_state[6], double env_vars[6]) {

	static Ball ball = Ball();
	robot_state_hands hands(robot_state[LEFT_HAND].x,robot_state[RIGHT_HAND].x,
			                robot_state[LEFT_HAND].xd,robot_state[RIGHT_HAND].xd,1);

	ball.integrate_ball_state(hands,dt);
	ball.get_state(ball_state);
	ball.get_env_params(env_vars);
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

	static vec3 ball_obs;
	static vec q0 = zeros<vec>(NDOF_ACTIVE);
	static joint qact;
	static joint qdes;
	static Player *robot = nullptr; // pointer to player
	static EKF filter = init_filter(0.001,0.001);

	if (options.reset) {
		for (int i = 0; i < NDOF_ACTIVE; i++) {
			qdes.q(i) = q0(i) = joint_state[active_dofs(i)+1].th;
			qdes.qd(i) = 0.0;
			qdes.qdd(i) = 0.0;
		}
		filter = init_filter(0.3,0.001);
		delete robot;
		robot = new Player(q0,filter,options);
		options.reset = false;
	}
	else {
		for (int i = 0; i < NDOF_ACTIVE; i++) {
			qact.q(i) = joint_state[active_dofs(i)+1].th;
			qact.qd(i) = joint_state[active_dofs(i)+1].thd;
			qact.qdd(i) = joint_state[active_dofs(i)+1].thdd;
		}
		fuse_blobs(blobs,ball_obs);
		robot->play(qact,ball_obs,qdes);
	}

	// update desired joint state
	for (int i = 0; i < NDOF_ACTIVE; i++) {
		joint_des_state[active_dofs(i)+1].th = qdes.q(i);
		joint_des_state[active_dofs(i)+1].thd = qdes.qd(i);
		joint_des_state[active_dofs(i)+1].thdd = qdes.qdd(i);
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

	static vec q0 = zeros<vec>(NDOF_ACTIVE);
	static vec6 ball_state;
	static joint qact;
	static joint qdes;
	static Player *robot = nullptr; // centered player
	static EKF filter = init_filter(0.001,0.001);

	if (options.reset) {
		for (int i = 0; i < NDOF_ACTIVE; i++) {
			qdes.q(i) = q0(i) = joint_state[active_dofs(i)+1].th;
			qdes.qd(i) = 0.0;
			qdes.qdd(i) = 0.0;
		}
		delete robot;
		robot = new Player(q0,filter,options);
		options.reset = false;
	}
	else {
		for (int i = 0; i < NDOF_ACTIVE; i++) {
			qact.q(i) = joint_state[active_dofs(i)+1].th;
			qact.qd(i) = joint_state[active_dofs(i)+1].thd;
			qact.qdd(i) = joint_state[active_dofs(i)+1].thdd;
		}

		for (int i = 0; i < NCART; i++) {
			ball_state(i) = sim_ball_state->x[i+1];
			ball_state(i+NCART) = sim_ball_state->xd[i+1];
		}
		robot->cheat(qact,ball_state,qdes);
	}

	// update desired joint state
	for (int i = 0; i < NDOF_ACTIVE; i++) {
		joint_des_state[active_dofs(i)+1].th = qdes.q(i);
		joint_des_state[active_dofs(i)+1].thd = qdes.qd(i);
		joint_des_state[active_dofs(i)+1].thdd = qdes.qdd(i);
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
