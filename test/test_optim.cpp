/*
 * test_optim.cpp
 *
 * Unit Tests for polynomial optimization
 *
 *  Created on: Feb 17, 2017
 *      Author: okoc
 */

#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE test_basketball
#endif

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <armadillo>
#include <thread>
#include "kinematics.h"
#include "optim.h"
#include "player.hpp"

const double PI = 3.1416;

using namespace arma;

/*
 * Initialize default posture for CB only for the active joints
 */
inline void init_default_posture(const bool right, vec & q0, ivec & active_dofs) {

	q0(0) = 0.0;
	q0(1) = -0.2;
	q0(2) = 0.0;
	q0(3) = 1.57;
	q0(4) = 0.0;
	q0(5) = 0.0;
	q0(6) = 0.0;

	if (right) {
		active_dofs = RIGHT_ARM;
	}
	else {
		active_dofs = LEFT_ARM;
	}
}

/*
 * Initialize default state for incoming basketball
 */
inline void init_default_basketball(vec6 & ball_state) {
	ball_state.zeros();
	const vec3 base_pendulum = {0.00, 0.9, 1.0};
	const double string_len = 1.0;
	const double basketball_radius = 0.1213;
	double theta_init = -45.0 * (PI/180.0);
	ball_state(X) = base_pendulum(X);
	ball_state(Y) = base_pendulum(Y) - (string_len + basketball_radius) * sin(theta_init);
	ball_state(Z) = base_pendulum(Z) - (string_len + basketball_radius) * cos(theta_init);
}

/*
 * Testing the optimization of a Basketball player touching a ball
 * attached to a string (moving in 2d: y and z axis)
 */
BOOST_AUTO_TEST_CASE(test_optim) {

	BOOST_TEST_MESSAGE("Testing the optimization...");
	int N = 1000;
	bool right_side = true;
	vec6 ball_state;
	arma_rng::set_seed(1);
	//arma_rng::set_seed_random();
	vec lb = zeros<vec>(2*NDOF_ACTIVE+1);
	vec ub = zeros<vec>(2*NDOF_ACTIVE+1);
	vec q0 = zeros<vec>(NDOF_ACTIVE);
	joint qact;
	spline_params poly;
	ivec active_dofs;
	init_default_basketball(ball_state);
	init_default_posture(right_side,q0,active_dofs);
	qact.q = q0;
	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);
	mat balls_pred = filter.predict_path(DT,N);
	optim_des params;
	params.Nmax = 1000;
	params.ball_pos = balls_pred.rows(X,Z);
	params.ball_vel = balls_pred.rows(DX,DZ);

	// right side
	Optim opt = Optim(q0,active_dofs,right_side);
	opt.set_des_params(&params);
	opt.update_init_state(qact);
	opt.run();
	bool update_right_side = opt.get_params(qact,poly);

	// left side
	init_default_posture(!right_side,q0,active_dofs);
	opt.update_rest_state(q0);
	opt.set_active_dofs(active_dofs);
	opt.right_arm = false;
	opt.run();
	bool update_left_side = opt.get_params(qact,poly);
	BOOST_TEST(update_right_side);
	BOOST_TEST(update_left_side);
}

/*
 * Testing whether the ball can be touched
 */
BOOST_AUTO_TEST_CASE(test_touch) {

	BOOST_TEST_MESSAGE("Testing if the ball can be touched...");

	arma_rng::set_seed_random();
	//arma_rng::set_seed(5);

	int N = 1000;
	joint qact;
	joint qdes;
	vec6 ball_state;
	vec3 ball_obs;
	double pos_left[NCART], pos_right[NCART];
	ivec active_dofs;
	init_default_basketball(ball_state);
	init_default_posture(true,qact.q,active_dofs);
	EKF filter = init_filter();
	player_flags flags;
	flags.active_dofs = active_dofs;
	flags.detach = false;
	Player robot = Player(qact.q,filter,flags);
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);
	mat balls_pred = filter.predict_path(DT,N);
	qdes.q = qact.q;
	mat xdes = zeros<mat>(2*NCART,N);

	for (int i = 0; i < N; i++) {

		// move the ball
		ball_obs = balls_pred.col(i).head(3);

		// play
		robot.play(qact, ball_obs, qdes);
		//robot.cheat(qact, balls_pred.col(i), qdes);

		// get cartesian state
		get_position(active_dofs,qdes.q.memptr(),pos_left,pos_right);
		xdes(X,i) = pos_right[X];
		xdes(Y,i) = pos_right[Y];
		xdes(Z,i) = pos_right[Z];
		xdes(X+NCART,i) = pos_left[X];
		xdes(Y+NCART,i) = pos_left[Y];
		xdes(Z+NCART,i) = pos_left[Z];

		//usleep(DT*1e6);
		qact.q = qdes.q;
		qact.qd = qdes.qd;
	}

	// test for intersection on Cartesian space
	//balls_pred.save("balls_pred.txt",csv_ascii);
	//xdes.save("robot_cart.txt",csv_ascii);

	// find the closest point between two curves
	mat err = xdes.rows(0,2) - balls_pred.rows(0,2);
	rowvec errnorms = sqrt(sum(err % err,0));
	uword idx = index_min(errnorms);
	BOOST_TEST_MESSAGE("Minimum dist between ball and robot: \n" << errnorms(idx));
	BOOST_TEST(errnorms(idx) < 0.1); // distance should be less than 10 cm
}

/*
 * Testing if the kinematics was copied correctly from SL
 */
BOOST_AUTO_TEST_CASE(test_kinematics) {

	// TODO: Could it be that the difference in default postures causes not-so-small deviations?

	BOOST_TEST_MESSAGE("Testing kinematics close to default posture...");

	ivec active_dofs = RIGHT_ARM;
	double q_active[NDOF_ACTIVE] = {-0.005,-0.186,-0.009,1.521,0.001,-0.001,-0.004,};
	double pos_left[NCART];
	double pos_right[NCART];
	double pos_des_right[NCART] = {0.272,0.336,0.179};

	get_position(active_dofs,q_active,pos_left,pos_right);
	BOOST_TEST(pos_right[0] == pos_des_right[0],boost::test_tools::tolerance(0.01));
	BOOST_TEST(pos_right[1] == pos_des_right[1],boost::test_tools::tolerance(0.01));
	BOOST_TEST(pos_right[2] == pos_des_right[2],boost::test_tools::tolerance(0.01));
}
