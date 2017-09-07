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
#include "utils.h"
#include "optim.h"
#include "player.hpp"

using namespace arma;
int randval;

/*
 * Initialize default posture for CB
 */
inline void init_default_posture(vec7 & q0) {

	q0(0) = 1.0;
	q0(1) = -0.2;
	q0(2) = -0.1;
	q0(3) = 1.8;
	q0(4) = -1.57;
	q0(5) = 0.1;
	q0(6) = 0.3;
}

/*
 * Testing Fixed Player (or Focused Player)
 */
BOOST_AUTO_TEST_CASE(test_optim) {

	cout << "Testing FP Trajectory Optimizer...\n";
	double lb[2*NDOF_ACTIVE+1], ub[2*NDOF_ACTIVE+1];
	double SLACK = 0.01;
	double Tmax = 1.0;
	joint qact;
	spline_params poly;

	// update initial parameters from lookup table
	std::cout << "Looking up a random ball entry..." << std::endl;
	arma_rng::set_seed(randval);
	//arma_rng::set_seed_random();
	vec::fixed<15> strike_params;
	vec6 ball_state;
	init_default_posture(qact.q);
	set_bounds(lb,ub,SLACK,Tmax);
	optim_des params;
	int N = 1000;
	params.Nmax = 1000;

	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);

	vec6 ball_pred;
	double time_land_des = 0.8;
	mat balls_pred = filter.predict_path(DT,N);

	Optim opt = Optim(qact.q.memptr(),lb,ub);
	opt.set_des_params(&params);
	opt.update_init_state(qact);
	opt.run();
	bool update = opt.get_params(qact,poly);

	BOOST_TEST(update);
}
