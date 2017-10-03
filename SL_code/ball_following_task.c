/*
 * ball_following_task.c
 *
 *  Created on: Sep 5, 2017
 *      Author: okan.koc
 */

/* Global includes */
#include "SL_system_headers.h"
#include "SL.h"
#include "SL_user.h"
#include "SL_task_servo.h"
#include "SL_collect_data.h"
#include "SL_dynamics.h"
#include "SL_tasks.h"

#define endEff_idx RIGHT_HAND

/* Internal variables */
static double start_time;
static double ball_speed = 1.0;
static double learning_rate = 0.1;

static Vector x_err;
static Vector q_des;
static Vector qd_des;
static Vector qdd_des;
static Vector gradient;
static Matrix Jt;
static Matrix subJ;

typedef struct {
	int status; //!< was ball detected reliably in cameras
	double pos[N_CART]; //!< ball center cartesian positions from cameras 1 and 2(after calibration)
} blob_state;

blob_state ball_obs;
SL_Cstate sim_ball_state;

#include "sl_basketball_interface.h"

/* Global functions */
void add_ball_following_task(void);

/* Local functions */
static int init_ball_following_task(void);
static int run_ball_following_task(void);
static int change_ball_following_task(void);

static void jacobian_transpose();

const int active_dofs[] = { 7, // number of active degrees of freedom
		              R_SFE, // 08
		  	  	  	  R_SAA, // 09
					  R_HR,  // 10
					  R_EB,  // 11
					  R_WR,  // 12
					  R_WFE, // 13
					  R_WAA, // 14
};

const int N_ACTIVE_DOFS = 7;


static int simulate_ball(void) {
	int j, i;
	static int last_frame_counter = -999;
	static int sendflag = FALSE;
	static double pos[N_CART+1];
	//static double start_pos[N_CART] = {0.272, 0.336, 0.180};
	static double start_pos[N_CART] = {0.372, 0.336, 0.280};

	// Ball simulator
	pos[0] = 0.1213;
	pos[_X_] = start_pos[0] + 0.1*sin(2.*PI*0.1*ball_speed*(task_servo_time-start_time));
	pos[_Y_] = start_pos[1];
	pos[_Z_] = start_pos[2];

	sendUserGraphics("basketball",&(pos[0]),4*sizeof(double));
	raw_blobs[BLOB1].x[1] = pos[_X_];
	raw_blobs[BLOB1].x[2] = pos[_Y_];
	raw_blobs[BLOB1].x[3] = pos[_Z_];

	sim_ball_state.x[1] = pos[_X_];
	sim_ball_state.x[2] = pos[_Y_];
	sim_ball_state.x[3] = pos[_Z_];
	sim_ball_state.xd[1] = 2.*PI*0.1*ball_speed*0.1*cos(2.*PI*0.1*ball_speed*(task_servo_time-start_time));
	for (i = 2; i <= 3; i++) {
		sim_ball_state.xd[i] = 0.0;
	}

	return TRUE;
}


static int init_ball_following_task(void) {

	bzero((char *)&(ball_obs), sizeof(ball_obs));
	bzero((char *)&(sim_ball_state), sizeof(sim_ball_state));

	// Anything done once on static variables and all memory allocations go here

	x_err    = my_vector(1,N_CART);
	q_des    = my_vector(1,N_ACTIVE_DOFS);
	qd_des   = my_vector(1,N_ACTIVE_DOFS);
	qdd_des  = my_vector(1,N_ACTIVE_DOFS);
	gradient = my_vector(1,N_ACTIVE_DOFS);
	Jt       = my_matrix(1,N_ACTIVE_DOFS,1,N_CART);
	subJ     = my_matrix(1,N_CART,1,N_ACTIVE_DOFS);

	int i;
	for (i = 1; i < N_ACTIVE_DOFS; i++) {
		q_des[i] = joint_default_state[active_dofs[i]].th;
	}

	// Add variables for MRDPLOT
	char string[100];
	for (i = _X_; i <= _Z_; ++i) {
		sprintf(string,"Error_%d",i);
		addVarToCollect(&(x_err[i]), string, "m", DOUBLE, FALSE);
	}
	updateDataCollectScript();

	start_time = servo_time;

	// Set the ball speed
	get_double("Ball speed",ball_speed,&ball_speed);

	// Start collecting data
	scd();

	return TRUE;
}


static int run_ball_following_task(void) {
	double t = servo_time - start_time; // Current time

	simulate_ball();

	//jacobian_transpose();
	cheat(joint_state, &sim_ball_state, joint_des_state);

	check_range(joint_des_state); // Check if the trajectory is safe
	SL_InverseDynamics(NULL, joint_des_state, endeff); // Feedforward control
	//SL_InverseDynamics(joint_state, joint_des_state, endeff);

	return TRUE;
}

static void jacobian_transpose() {

	int i, j;

	// We use the Jacobian transpose method: we want to minimize the squared position error via gradient descent
	//  1) q(t+1) = q(t) - learning_rate * gradient
	//     gradient = d(0.5 * x_err' * x_err) / d(q) = - J' * x_err
	//     q(t+1) = q(t) + learning_rate * J' * x_err
	//  2) qd(t+1) = learning_rate * gradient

	// Extract only the needed rows from the Jacobian
	for(i = 1; i <= N_CART; i++) {
		for(j = 1; j <= N_ACTIVE_DOFS; j++) {
			subJ[i][j] = J[i][active_dofs[j]];
		}
	}

	// Compute the Jacobian transpose
	mat_trans(subJ, Jt);

	// Compute the error in the task space
	double max_error = 0.2;
	for(i = 1; i <= N_CART; i++) {
		x_err[i] = (raw_blobs[BLOB1].x[i] - cart_state[RIGHT_HAND].x[i]);
		x_err[i] = x_err[i] > max_error ? max_error : x_err[i] < -max_error ? -max_error : x_err[i]; // Clip error if it is too high to avoid crazy movements
	}

	//printf("Rawblobs: [%f,%f,%f]\n", raw_blobs[BLOB1].x[1], raw_blobs[BLOB1].x[2], raw_blobs[BLOB1].x[3]);

	// Compute the gradient
	mat_vec_mult(Jt,x_err,gradient);

	// Perform a step towards the desired position
	for (i = 1; i <= N_ACTIVE_DOFS; ++i) {
		q_des[i]  = joint_state[active_dofs[i]].th + learning_rate * gradient[i];
		qd_des[i] = learning_rate * gradient[i];
		qdd_des[i] = 0.0;
	}

	/* ------------------------ Solution 1 ------------------------ */
	// Directly set joint_des_state and let SL do everything for you.
	// It works if there is a good model of the robot.
	for (i = 1; i <= N_ACTIVE_DOFS; ++i) {
	    joint_des_state[active_dofs[i]].th   = q_des[i];
	    joint_des_state[active_dofs[i]].thd  = qd_des[i];
	    joint_des_state[active_dofs[i]].thdd = qdd_des[i];
	    joint_des_state[active_dofs[i]].uff  = 0.0;
	}
	/* ------------------------------------------------------------ */

	/* ------------------------ Solution 2 ------------------------ */
	// Override SL PD controller and use your own.
	// Manually set the desired acceleration using your PD and keep
	// the current position and velocity.
	/*float K_p[N_DOFS] = {200., 300., 100., 50., 10., 10., 2.5}; // Gains
	float K_d[N_DOFS] = {7., 15., 5., 2.5, 0.3, 0.3, 0.075};

	for (i = R_EB; i <= R_WAA; ++i) {
		joint_des_state[i].th   = joint_state[i].th;
		joint_des_state[i].thd  = joint_state[i].thd;
		joint_des_state[i].thdd = K_d[i-1] * (qd_des[i] - joint_state[i].thd) + K_p[i-1] * (q_des[i] - joint_state[i].th);
	}*/
	/* ------------------------------------------------------------ */

}


static int change_ball_following_task(void) {
	return TRUE;
}


void add_ball_following_task( void ) {
	addTask("Ball Following Task", init_ball_following_task,
			run_ball_following_task, change_ball_following_task);
}

