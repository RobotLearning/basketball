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

using namespace arma;

/*
 * Initialize default posture for CB only for the active joints
 */
inline void init_default_posture(const bool right, vec & q0) {

	q0(0) = 0.0;
	q0(1) = -0.2;
	q0(2) = 0.0;
	q0(3) = 1.57;
	q0(4) = 0.0;
	q0(5) = 0.0;
	q0(6) = 0.0;
}

/*
 * Testing if the kinematics was copied correctly from SL
 */
BOOST_AUTO_TEST_CASE(test_kinematics) {

	// TODO: Could it be that the difference in default postures causes not-so-small deviations?

	BOOST_TEST_MESSAGE("Testing kinematics close to default posture...");

	ivec active_dofs = join_vert(LEFT_ARM,RIGHT_ARM);
	//double q_active[NDOF_OPT] = {-0.005,-0.186,-0.009,1.521,0.001,-0.001,-0.004,};
	double q_active[NDOF_ACTIVE] = {0.0, -0.2, 0.0, 1.57, 0.0, 0.0, 0.0,
								 0.0, -0.2, 0.0, 1.57, 0.0, 0.0, 0.0};
	vec3 pos_left;
	vec3 pos_right;
	vec3 pos_des_left = {-0.27,0.34,0.21};
	vec3 pos_des_right = {0.27,0.34,0.21};

	calc_cart_pos(active_dofs,q_active,pos_left,pos_right);
	BOOST_TEST(approx_equal(pos_left,pos_des_left,"absdiff", 0.002)); //boost::test_tools::tolerance(0.01)
	BOOST_TEST(approx_equal(pos_right,pos_des_right,"absdiff", 0.002)); //boost::test_tools::tolerance(0.01)
}

/*
 * Testing the optimization of a Basketball player touching a ball
 * attached to a string (moving in 2d: y and z axis)
 */
BOOST_AUTO_TEST_CASE(test_optim) {

	BOOST_TEST_MESSAGE("Testing the optimization for touch/hit...");
	int N = 1000;
	arma_rng::set_seed(1);
	//arma_rng::set_seed_random();
	vec q0 = zeros<vec>(NDOF_OPT);
	joint qact;
	spline_params poly;
	Ball ball = Ball();
	vec6 ball_state = ball.get_state();
	init_default_posture(true,q0);
	qact.q = join_vert(q0,q0);
	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);
	mat balls_pred = filter.predict_path(DT,N);
	optim_des params;
	params.Nmax = 1000;
	params.ball_pos = balls_pred.rows(X,Z);
	params.ball_vel = balls_pred.rows(DX,DZ);

	BOOST_TEST_MESSAGE("TOUCHING from the left side...");
	Optim *opt = new Optim(q0,false);
	opt->set_des_params(&params);
	opt->update_init_state(qact);
	opt->run();
	bool update_left_side = opt->get_params(qact,poly);
	delete opt;

	BOOST_TEST_MESSAGE("TOUCHING from the right side...");
	opt = new Optim(q0,true);
	opt->set_des_params(&params);
	opt->update_init_state(qact);
	opt->run();
	bool update_right_side = opt->get_params(qact,poly);
	delete opt;

	BOOST_TEST(update_right_side);
	BOOST_TEST(update_left_side);

	BOOST_TEST_MESSAGE("HITTING from the left side...");
	opt = new Optim(q0,false,false);
	opt->set_des_params(&params);
	opt->update_init_state(qact);
	opt->run();
	update_right_side = opt->get_params(qact,poly);
	delete opt;

	BOOST_TEST_MESSAGE("HITTING from the right side...");
	opt = new Optim(q0,true,false);
	opt->set_des_params(&params);
	opt->update_init_state(qact);
	opt->run();
	update_left_side = opt->get_params(qact,poly);
	delete opt;

	BOOST_TEST(update_right_side);
	BOOST_TEST(update_left_side);
}

/*
 * Testing whether the ball can be touched
 */
BOOST_AUTO_TEST_CASE(test_touch_player) {

	BOOST_TEST_MESSAGE("Testing if the ball can be TOUCHED with RIGHT HAND...");

	arma_rng::set_seed_random();
	//arma_rng::set_seed(5);
	const double basketball_radius = 0.1213;
	int N = 500;
	joint qact;
	joint qdes;
	vec7 q0;
	vec6 ball_state;
	vec3 ball_obs;
	vec3 pos_left, pos_right;
	ivec active_dofs = join_vert(LEFT_ARM,RIGHT_ARM);
	init_default_posture(true,q0);
	qdes.q = join_vert(q0,q0);
	qact.q = qdes.q;
	EKF filter = init_filter();
	player_flags flags;
	flags.detach = false;
	flags.verbosity = 3;
	flags.optim_type = RIGHT_HAND_OPT;
	flags.touch = true;
	Ball ball = Ball();
	Player robot = Player(qact.q,filter,flags);
	mat66 P; P.eye();
	mat xdes = zeros<mat>(2*NCART,N);
	mat balls = zeros<mat>(NCART,N);

	for (int i = 0; i < N; i++) {

		// move the ball
		ball.integrate_ball_state(DT);
		ball_state = ball.get_state();
		ball_obs = ball_state.head(NCART);
		balls.col(i) = ball_obs;

		// play
		//robot.play(qact, ball_obs, qdes);
		robot.cheat(qact, ball_state, qdes);

		// get cartesian state
		calc_cart_pos(active_dofs,qdes.q.memptr(),pos_left,pos_right);
		xdes.col(i) = join_vert(pos_left,pos_right);

		usleep(DT*1e6);
		qact.q = qdes.q;
		qact.qd = qdes.qd;
	}

	// test for intersection on Cartesian space
	balls.save("balls_pred.txt",csv_ascii);
	xdes.save("robot_cart.txt",csv_ascii);

	// find the closest point between two curves
	mat diff = xdes.rows(0,2) - balls;
	rowvec diff_norms = sqrt(sum(diff % diff,0));
	uword idx = index_min(diff_norms);
	BOOST_TEST_MESSAGE("Minimum dist between ball and robot: \n" << diff_norms(idx));
	BOOST_TEST(diff_norms(idx) <= basketball_radius, boost::test_tools::tolerance(0.01)); // distance should be less than 10 cm
}
