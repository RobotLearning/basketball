/*
 * basketball_task.c
 *
 *  Created on: Sep 6, 2017
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
#include "SL_integrate.h"

/* Internal variables */
static double start_time;

typedef struct {
	int status; //!< was ball detected reliably in cameras
	double pos[N_CART]; //!< ball center cartesian positions from cameras 1 and 2(after calibration)
} blob_state;

blob_state ball_obs;
SL_Cstate sim_ball_state;
int simulation;

#include "sl_basketball_interface.h"

/* Global functions */
void add_basketball_task(void);

/* Local functions */
static int init_basketball_task(void);
static int run_basketball_task(void);
static int change_basketball_task(void);


/**
 * Simulates a ball attached to a string in SIMULATION mode
 * Otherwise get the ball observations from SL's blobs[1] structure
 */
static int update_ball_obs(void) {

	int j, i;
	static int firsttime = TRUE;
	static const double dt = 0.002;
	static double theta_dot = 0.0;
	static double theta = 45.0 * (PI/180.0);
	static double ball_state[6];
	static double env_vars[6];

	if (simulation) {
		update_base(&base_state,&base_orient);
		integrate_ball_state(dt, cart_state, ball_state, env_vars);
		sendUserGraphics("basketball_pendulum",env_vars,6*sizeof(double));

		// fill the ball_obs structure for communicating to basketball lib
		ball_obs.status = TRUE;
		for (i = 0; i < 3; i++) {
			sim_ball_state.x[i+1] = ball_state[i];
			sim_ball_state.xd[i+1] = ball_state[i+3];
			ball_obs.pos[i] = ball_state[i];
		}
	}
	else {
		// TODO: Is it blobs[1] that gets the visual ball info?
		ball_obs.status = blobs[1].status;
		for (i = 0; i < 3; i++) {
			ball_obs.pos[i] = blobs[1].blob.x[i+1];
		}
	}

	return TRUE;
}


static int init_basketball_task(void) {

	int i, ready;
	//static double zeroGain[N_DOFS];

	bzero((char *)&(ball_obs), sizeof(ball_obs));
	bzero((char *)&(sim_ball_state), sizeof(sim_ball_state));

	// Start collecting data
	scd();

	load_options();

	for (i = 1; i <= N_DOFS; i++) {
		joint_des_state[i].th = joint_default_state[i].th;
		//printf("default_state[%d] = %f\n", i, joint_default_state[i].th);
	}

	get_int("Use simulated ball?", FALSE, &simulation);
	if (simulation) {
		printf("Testing BASKETBALL on the simulator...\n");
	}

	/* ready to go */
	ready = 999;
	while (ready == 999) {
		if (!get_int("Enter 1 to start or anything else to abort ...",ready, &ready))
			return FALSE;
	}
	if (ready != 1)
		return FALSE;

	start_time = servo_time;
	freezeBase(1);
	//changePIDGains(zeroGain, zeroGain, zeroGain);
	return TRUE;
}


static int run_basketball_task(void) {

	double t = servo_time - start_time; // Current time
	int i, j;

	//for (i = 1; i <= N_DOFS; i++) {
	//	joint_des_state[i].th = joint_default_state[i].th;
	//	//printf("default_state[%d] = %f\n", i, joint_default_state[i].th);
	//}

	//printf("t = %f\n",t);
	//printf("basec = [%f,%f,%f]\n",base_state.x[1],base_state.x[2],base_state.x[3]);
	//printf("baseo = [%f,%f,%f,%f]\n",base_orient.q[1],base_orient.q[2],base_orient.q[3],base_orient.q[4]);
	//printf("effx = [%f,%f,%f]\n",endeff[RIGHT_HAND].x[1],endeff[RIGHT_HAND].x[2],endeff[RIGHT_HAND].x[3]);
	//printf("effo = [%f,%f,%f]\n",endeff[RIGHT_HAND].a[1],endeff[RIGHT_HAND].a[2],endeff[RIGHT_HAND].a[3]);
	//printf("x = [%f,%f,%f]\n",cart_state[RIGHT_HAND].x[1],cart_state[RIGHT_HAND].x[2],cart_state[RIGHT_HAND].x[3]);
	//printf("LEFT HAND = [%f,%f,%f]\n",cart_des_state[LEFT_HAND].x[1],cart_des_state[LEFT_HAND].x[2],cart_des_state[LEFT_HAND].x[3]);
	//printf("RIGHT HAND = [%f,%f,%f]\n",cart_des_state[RIGHT_HAND].x[1],cart_des_state[RIGHT_HAND].x[2],cart_des_state[RIGHT_HAND].x[3]);

	update_ball_obs();
	//cheat(joint_state, &sim_ball_state, joint_des_state);
	play(joint_state, &ball_obs, joint_des_state); // basketball lib generates trajectories for touching the ball

	check_range(joint_des_state); // Check if the trajectory is safe
	SL_InverseDynamics(NULL, joint_des_state, endeff); // Feedforward control
	//SL_InverseDynamics(joint_state, joint_des_state, endeff);

	return TRUE;
}


static int change_basketball_task(void) {
	return TRUE;
}


void add_basketball_task( void ) {
	addTask("Basketball Task", init_basketball_task, run_basketball_task, change_basketball_task);
}


