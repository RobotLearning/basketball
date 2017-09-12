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

/* Internal variables */
static double start_time;
static double ball_speed = 1.0;

typedef struct {
	int status; //!< was ball detected reliably in cameras
	double pos[N_CART]; //!< ball center cartesian positions from cameras 1 and 2(after calibration)
} blob_state;

blob_state ball_obs;

#include "sl_basketball_interface.h"

/* Global functions */
void add_basketball_task(void);

/* Local functions */
static int init_basketball_task(void);
static int run_basketball_task(void);
static int change_basketball_task(void);


/**
 * Simulates a ball attached to a string
 */
static int simulate_ball(void) {
	int j, i;
	static double last_time;
	static int last_frame_counter = -999;
	static int sendflag = FALSE;
	static double ball_speed = 1.0;
	static double gravity = -9.8;
	static double friction = 0.03;
	static int firsttime = TRUE;
	// Angle
	static double theta_dot = 0.0;
	static double theta = 45.0 * (PI/180.0);

	if (firsttime) {
		last_time = servo_time;
		firsttime = FALSE;
	}

	theta_dot += (servo_time - last_time) * (gravity * sin(theta) - friction * theta_dot);
	theta += (servo_time - last_time) * theta_dot;

	sendUserGraphics("basketball_pendulum",&theta,sizeof(double));

	// fill the ball_obs structure for communicating to basketball lib
	double ball_centre[N_CART+1];
	const double base_pendulum[N_CART+1] = {0.0, 0.0, 0.9, 0.9};
	const double basketball_radius = 0.1213;
	const double string_len = 1.0;
	ball_centre[1] = base_pendulum[1]; //x is fixed
	ball_centre[2] = base_pendulum[2] - (string_len + basketball_radius) * sin(theta);
	ball_centre[3] = base_pendulum[3] - (string_len + basketball_radius) * cos(theta);

	ball_obs.status = TRUE;
	for (i = 0; i < 3; i++) {
		ball_obs.pos[i] = ball_centre[i+1];
	}

	last_time = servo_time;
	return TRUE;
}


static int init_basketball_task(void) {

	bzero((char *)&(ball_obs), sizeof(ball_obs));
	start_time = servo_time;

	// Set the ball speed
	get_double("Ball speed",ball_speed,&ball_speed);

	// Start collecting data
	scd();

	return TRUE;
}


static int run_basketball_task(void) {

	double t = servo_time - start_time; // Current time
	int i, j;

	simulate_ball();
	play(joint_state, &ball_obs, joint_des_state); // basketball lib generates trajectories for touching the ball

	check_range(joint_des_state); // Check if the trajectory is safe

	return TRUE;
}


static int change_basketball_task(void) {
	return TRUE;
}


void add_basketball_task( void ) {
	addTask("Basketball Task", init_basketball_task, run_basketball_task, change_basketball_task);
}


