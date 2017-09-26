/*
 * ball.cpp
 *
 *  Created on: Sep 26, 2017
 *      Author: okan.koc
 */

#include <armadillo>
#include "constants.h"
#include "ball.h"

using namespace arma;

/**
 * @brief Function that integrates a basketball for an outside filter.
 *
 * Function exposes the integration to filters, e.g. an EKF.
 * They can use then to apply predict() using the this function pointer.
 *
 * The basketball is attached to a string so first finding the angle responsible
 * for the observation, integrating it and then finding the next ball state
 *
 */
vec calc_next_ball(const vec & xnow, const double dt, const void *fp) {

	vec3 ball_pos = xnow.head(3);
	vec3 ball_vel = xnow.tail(3);
	double theta, theta_dot;

	// to integrate the ball first extract the theta and theta dot
	// integrate them and then revert back to current ball pos
	calc_angle_from_ball(ball_pos, ball_vel, theta, theta_dot);
	ball_pendulum_model(dt, theta, theta_dot);
	calc_ball_from_angle(theta, theta_dot, ball_pos, ball_vel);
	return join_vert(ball_pos,ball_vel);
}


/**
 * @brief Integrate the angle state of the pendulum string using Symplectic Euler.
 */
void ball_pendulum_model(const double dt, double & theta, double & theta_dot) {

	theta_dot += dt * (gravity * sin(theta) - friction * theta_dot);
	theta += dt * theta_dot;
}

/**
 * Extract angle state from ball positions and velocities (only Y-velocity is used!)
 */
void calc_angle_from_ball(const vec3 & ball_pos, const vec3 & ball_vel, double & theta, double & theta_dot) {

	vec3 pos = ball_pos - base_pendulum;
	theta = atan(pos(Y) / pos(Z));
	theta_dot = -ball_vel(Y) / ((string_len + basketball_radius) * cos(theta));
}

/**
 * @brief Construct ball state purely from pendulum angle pos and velocity.
 *
 * This method disregards previous computations on the ball state (if any) and constructs the whole
 * ball state from angle pos and velocities.
 */
void calc_ball_from_angle(const double theta, const double theta_dot, vec3 & ball_pos, vec3 & ball_vel) {

	// we're in the third quadrant
	ball_pos(X) = base_pendulum(X);
	ball_pos(Y) = base_pendulum(Y) - (string_len + basketball_radius) * sin(theta);
	ball_pos(Z) = base_pendulum(Z) - (string_len + basketball_radius) * cos(theta);
	ball_vel(X) = 0.0;
	ball_vel(Y) = -(string_len + basketball_radius) * theta_dot * cos(theta);
	ball_vel(Z) = (string_len + basketball_radius) * theta_dot * sin(theta);
}

/**
 * @brief Modify the incoming basketball velocity and string angle velocity in case of contact.
 *
 * Contact occurs if the norm of the difference between the arm and the ball is less than
 * or equal to the ball radius.
 * The contact model is an IDEAL MOMENTUM EXCHANGE assuming mass_robot >> mass_ball!
 */
void check_for_contact(const vec3 & robot_pos, const vec3 & robot_vel, const vec3 & ball_pos,
					   const double theta, vec3 & ball_vel, double & theta_dot) {

	if (norm(robot_pos - ball_pos) <= basketball_radius) { // contact occurs

		// change velocities
		ball_vel = -ball_vel + 2*robot_vel;
		// change string angle velocity
		theta_dot = -ball_vel(Y) / ((string_len + basketball_radius) * cos(theta));
	}
}


