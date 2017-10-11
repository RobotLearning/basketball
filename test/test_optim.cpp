/*
 * test_optim.cpp
 *
 * Unit Tests for polynomial optimization applied to Basketball
 *
 *  Created on: September, 2017
 *      Author: okoc
 */

#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE test_basketball
#endif

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <iostream>
#include <armadillo>
#include <thread>
#include "kinematics.h"
#include "optim.h"
#include "player.hpp"

using namespace arma;
namespace data = boost::unit_test::data;

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

/**
 * @brief Calculate minimum distance between robot hand and ball
 */
inline double calc_min_distance(const alg & optim_type, const mat & xdes, const mat & balls) {

	// first get left hand min distance
	mat diff = xdes.rows(0,2) - balls;
	rowvec diff_norms = sqrt(sum(diff % diff,0));
	uword idx = index_min(diff_norms);
	double min_dist_left = diff_norms(idx);
	// then get right hand min distance
	diff = xdes.rows(3,5) - balls;
	diff_norms = sqrt(sum(diff % diff,0));
	idx = index_min(diff_norms);
	double min_dist_right = diff_norms(idx);

	if (optim_type == LEFT_HAND_OPT) {
		return min_dist_left;
	}
	else if (optim_type == RIGHT_HAND_OPT) {
		return min_dist_right;
	}
	else { // both hands were optimized
		return fmax(min_dist_left,min_dist_right);
	}
}

/*
 * Testing if the kinematics was copied correctly from SL
 */
BOOST_AUTO_TEST_CASE(test_kinematics) {

	// TODO: Could it be that the difference in default postures causes not-so-small deviations?

	BOOST_TEST_MESSAGE("Testing kinematics close to default posture...");

	//double q_active[NDOF_OPT] = {-0.005,-0.186,-0.009,1.521,0.001,-0.001,-0.004,};
	double q_active[NDOF_ACTIVE] = {0.0, -0.2, 0.0, 1.57, 0.0, 0.0, 0.0,
								 0.0, -0.2, 0.0, 1.57, 0.0, 0.0, 0.0};
	vec3 pos_left;
	vec3 pos_right;
	vec3 pos_des_left = {-0.274,0.33,0.213};
	vec3 pos_des_right = {0.274,0.33,0.213};
    //RIGHT HAND = [0.273850,0.329690,0.213394]
	//LEFT HAND = [-0.274531,0.329705,0.213198]

	calc_cart_pos(active_dofs,q_active,pos_left,pos_right);
	//cout << endl << "POS_DES_LEFT: " << pos_des_left.t() << " vs. POS_LEFT:" << pos_left.t();
	//cout << endl << "POS_DES_RIGHT: " << pos_des_right.t() << " vs. POS_RIGHT:" << pos_right.t();
	BOOST_TEST(approx_equal(pos_left,pos_des_left,"absdiff", 0.002)); //boost::test_tools::tolerance(0.01)
	BOOST_TEST(approx_equal(pos_right,pos_des_right,"absdiff", 0.002)); //boost::test_tools::tolerance(0.01)
}

/*
 * Compare jacobian multiplied velocities to numerical differentiated cartesian positions
 */
BOOST_AUTO_TEST_CASE(test_jacobian) {

	BOOST_TEST_MESSAGE("\nComparing exact geometric jacobian to numerical diff. of kinematics ...");

	const double dt = 1e-5;
	//double q_active[NDOF_OPT] = {-0.005,-0.186,-0.009,1.521,0.001,-0.001,-0.004,};
	double q_active[NDOF_ACTIVE] = {0.0, -0.2, 0.0, 1.57, 0.0, 0.0, 0.0,
								    0.0, -0.2, 0.0, 1.57, 0.0, 0.0, 0.0};
	double qdot_active[NDOF_ACTIVE];
	double q_perturb[NDOF_ACTIVE];

	vec3 pos_left, pos_right, vel_left, vel_right;
	vec3 pos_diff_left, pos_diff_right, vel_diff_left, vel_diff_right;

	calc_cart_pos(active_dofs,q_active,pos_left,pos_right);

	for (int i = 0; i < NDOF_ACTIVE; i++) {

		// perturb q
		for (int j = 0; j < NDOF_ACTIVE; j++) {
			qdot_active[j] = 0.0;
			q_perturb[j] = q_active[j];
		}
		qdot_active[i] = 1.0;
		q_perturb[i] += qdot_active[i] * dt;

		// capture new cart. pos
		calc_cart_pos(active_dofs,q_perturb,pos_diff_left,pos_diff_right);

		// and num. diff. positions
		vel_diff_left = (pos_diff_left - pos_left) / dt;
		vel_diff_right = (pos_diff_right - pos_right) / dt;

		// compare with jacobian calculated positions
		calc_cart_pos_and_vel(active_dofs,q_active,qdot_active,pos_left,pos_right,vel_left,vel_right);

		/*cout << endl;
		cout << "VEL_LEFT: " << vel_left.t() << " vs. " << vel_diff_left.t();
		cout << "VEL_RIGHT: " << vel_right.t() << " vs. " << vel_diff_right.t();
		cout << endl;*/
		BOOST_TEST(approx_equal(vel_left,vel_diff_left,"abs_diff",0.002));
		BOOST_TEST(approx_equal(vel_right,vel_diff_right,"abs_diff",0.002));
	}
}

/*
 * Testing the optimization of a Basketball player touching a ball
 * attached to a string (moving in 2d: y and z axis)
 */
BOOST_DATA_TEST_CASE(test_optim, data::xrange(2) * data::xrange(2), touch, hand) {

	std::string touch_str[2] = {"TOUCH", "HIT"};
	std::string hand_str[2] = {"LEFT HAND", "RIGHT HAND"};
	cout << "\nTesting the optimization for " << touch_str[touch]
		 << " with " << hand_str[hand] << endl;
	int N = 1000;
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

	Optim opt = Optim(q0,hand,touch);
	opt.set_des_params(&params);
	opt.update_init_state(qact);
	opt.run();
	BOOST_TEST(opt.get_params(qact,poly));
}

/*
 * Testing whether the ball can be touched.
 * The combinations are : LEFT HAND/RIGHT HAND/BOTH HANDS, HIT/TOUCH, PLAY/CHEAT
 */
BOOST_AUTO_TEST_CASE(test_player) {


	BOOST_TEST_MESSAGE("\nTesting the player for HITTING with BOTH HANDS");

	const double basketball_radius = 0.1213;
	int N = 1000;
	joint qact;
	joint qdes;
	vec7 q0;
	vec6 ball_state;
	vec3 ball_obs;
	init_default_posture(true,q0);
	qdes.q = join_vert(q0,q0);
	qact.q = qdes.q;
	EKF filter = init_filter();
	player_flags flags;
	flags.detach = true;
	flags.verbosity = 3;
	flags.optim_type = BOTH_HAND_OPT;
	flags.touch = true;
	Ball ball = Ball();
	Player robot = Player(qact.q,filter,flags);
	mat66 P; P.eye();
	mat xdes = zeros<mat>(2*NCART,N);
	mat balls_pos = zeros<mat>(NCART,N);
	robot_hands hands;

	for (int i = 0; i < N; i++) {

		// move the ball
		ball.integrate_ball_state(hands,DT);
		ball_state = ball.get_state();
		ball_obs = ball_state.head(NCART);
		balls_pos.col(i) = ball_obs;

		robot.play(qact, ball_obs, qdes);
		//robot.cheat(qact, ball_state, qdes);

		// get cartesian state
		//calc_cart_pos(active_dofs,qdes.q.memptr(),pos_left,pos_right);
		calc_cart_pos_and_vel(active_dofs, qdes, hands);

		xdes.col(i) = join_vert(hands.left_pos, hands.right_pos);

		usleep(DT*1e6);
		qact.q = qdes.q;
		qact.qd = qdes.qd;
	}

	// test for intersection on Cartesian space
	balls_pos.save("balls_pos.txt",csv_ascii);
	xdes.save("robot_cart.txt",csv_ascii);

	// find the closest point between two curves
	double min_dist = calc_min_distance(flags.optim_type,xdes,balls_pos);
	cout << "Minimum dist between ball and robot: \n" << min_dist << endl;
	BOOST_TEST(min_dist <= basketball_radius + 0.01); // distance should be less than radius + 1cm tol.
}
