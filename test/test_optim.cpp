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
inline void init_default_posture(vec & q0) {

	q0(0) = 0.0;
	q0(1) = -0.2;
	q0(2) = 0.0;
	q0(3) = 1.57;
	q0(4) = 0.0;
	q0(5) = 0.0;
	q0(6) = 0.0;
}

/*
 * Initialize default state for incoming basketball
 */
inline void init_default_basketball(vec6 & ball_state) {
	ball_state.zeros();
	const vec3 base_pendulum = {0.0, 1.0, 1.0};
	const double string_len = 1.0;
	const double basketball_radius = 0.1213;
	double theta_init = -45.0 * (PI/180.0);
	ball_state(X) = 0.0;
	ball_state(Y) = base_pendulum(Y) - (string_len + basketball_radius) * sin(theta_init);
	ball_state(Z) = base_pendulum(Z) - (string_len + basketball_radius) * cos(theta_init);
}

/*
 * Testing Basketball player touching a ball
 * attached to a string (moving in 2d: y and z axis)
 */
BOOST_AUTO_TEST_CASE(test_touch_basketball) {

	BOOST_TEST_MESSAGE("Testing if robot can touch the basketball...");
	int N = 1000;
	vec6 ball_state;
	arma_rng::set_seed(1);
	//arma_rng::set_seed_random();
	double lb[2*NDOF_ACTIVE+1], ub[2*NDOF_ACTIVE+1];
	joint qact;
	spline_params poly;
	ivec active_dofs = {R_SFE, R_SAA, R_HR, R_EB, R_WR, R_WFE, R_WAA};
	set_bounds(active_dofs,0.01,1.0,lb,ub);
	init_default_basketball(ball_state);
	init_default_posture(qact.q);
	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);
	mat balls_pred = filter.predict_path(DT,N);

	optim_des params;
	params.Nmax = 1000;
	params.ball_pos = balls_pred.rows(X,Z);
	params.ball_vel = balls_pred.rows(DX,DZ);

	for (int i = 0; i < 2*NDOF_ACTIVE+1; i++) {
		printf("lb[%d] = %f, ub[%d] = %f\n", i,lb[i],i,ub[i]);
	}

	Optim opt = Optim(qact.q.memptr(),lb,ub);
	opt.set_des_params(&params);
	opt.set_active_dofs(active_dofs);
	opt.update_init_state(qact);
	opt.run();
	bool update = opt.get_params(qact,poly);

	BOOST_TEST(update);
}
