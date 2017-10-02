/*
 * ball.cpp
 *
 *  Created on: Sep 26, 2017
 *      Author: okan.koc
 */

#include <boost/program_options.hpp>
#include <armadillo>
#include "constants.h"
#include "ball.h"

using namespace arma;

/**
 * Initialize ball state by guessing it from the string angle
 */
Ball::Ball() {

	load_params("ball.cfg");
	calc_ball_from_angle(param,pos,vel);
}

/**
 * @brief Load ball prediction and other SIM parameters from a CONFIG file
 * @param file_name_relative Relative file name (base is basketball)
 *
 */
void Ball::load_params(const std::string & file_name_relative) {

	const double deg2rad = PI/180;
	namespace po = boost::program_options;
	using namespace std;
	string home = std::getenv("HOME");
	string config_file = home + "/basketball/" + file_name_relative;

	try {
		// Declare a group of options that will be
		// allowed in config file
		po::options_description config("Configuration");
		config.add_options()
			("verbose", po::value<bool>(&verbose)->default_value(false),
					"VERBOSITY OF CONTACT")
			("ball_radius", po::value<double>(&param.radius)->default_value(0.1213),
					"RADIUS OF BALL")
			//("gravity", po::value<double>(&param.gravity)->default_value(-9.80),
			//		"GRAVITY")
			("friction", po::value<double>(&param.friction)->default_value(0.0),
					"FRICTION OF PENDULUM STRING")
			("string_length", po::value<double>(&param.string_len)->default_value(1.0),
					"LENGTH OF PENDULUM STRING")
			("base_pendulum", po::value<vector<double>>(&param.base_pendulum)->multitoken(),
					"BASE OF PENDULUM AS A VECTOR")
			("theta_init", po::value<double>(&param.theta)->default_value(-45.0),
					"INITIAL PENDULUM ANGLE (DEGREES)")
			("theta_dot_init", po::value<double>(&param.theta_dot)->default_value(0.0),
										"INITIAL PENDULUM ANGULAR VELOCITY (DEG/SEC)")
					;
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
	param.theta *= deg2rad;
	param.theta_dot *= deg2rad;
}

/**
 * @brief Integrate ball positions and velocities with Symplectic Euler
 *
 * To integrate the ball first extract the theta and theta dot
 * integrate them and then revert back to current ball pos
 */
void Ball::integrate_ball_state(double dt) {

	ball_pendulum_model(dt, param);
	calc_ball_from_angle(param, pos, vel);
}

/**
 * @brief Integrate ball state and check for contact with robot.
 *
 * To integrate the ball first extract the theta and theta dot
 * integrate them and then revert back to current ball pos
 */
void Ball::integrate_ball_state(const robot_state_hands & robot, double dt) {

	ball_pendulum_model(dt, param);
	calc_ball_from_angle(param, pos, vel);

	if (CHECK_CONTACTS) {
		check_for_contact(robot.left_pos, robot.left_vel, pos, verbose, param, vel);
		check_for_contact(robot.right_pos, robot.right_vel, pos, verbose, param, vel);
	}
}

void Ball::set_state(const vec6 & ball_state) {

	pos = ball_state.head(NCART);
	vel = ball_state.tail(NCART);
}

void Ball::set_state(const double ball_state[2*NCART]) {

	for (int i = 0; i < NCART; i++) {
		pos(i) = ball_state[i];
		vel(i) = ball_state[i+NCART];
	}
}

vec6 Ball::get_state() const {

	return join_vert(pos,vel);
}

void Ball::get_state(double ball_state[2*NCART]) const {

	for (int i = 0; i < NCART; i++) {
		ball_state[i] = pos(i);
		ball_state[i+NCART] = vel(i);
	}
}

void Ball::get_env_params(ball_params & par) const {
	par = param;
}

/**
 * @brief Return environmental constants (loaded from a file by constructor)
 */
void Ball::get_env_params(double env_vars[6]) const {

	env_vars[0] = param.string_len;
	env_vars[1] = param.radius;
	env_vars[2] = param.base_pendulum[X];
	env_vars[3] = param.base_pendulum[Y];
	env_vars[4] = param.base_pendulum[Z];
	env_vars[5] = param.theta;
}

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

	static Ball ball = Ball();
	ball.set_state(xnow);
	ball.integrate_ball_state(dt);
	return ball.get_state();
}


/**
 * @brief Integrate the angle state of the pendulum string using Symplectic Euler.
 */
void ball_pendulum_model(const double dt, ball_params & param) {

	param.theta_dot += dt * (param.gravity * sin(param.theta) - param.friction * param.theta_dot);
	param.theta += dt * param.theta_dot;
}

/**
 * Extract angle state from ball positions and velocities (only Y-velocity is used!)
 */
void calc_angle_from_ball(const vec3 & ball_pos, const vec3 & ball_vel, ball_params & param) {

	vec3 pos = ball_pos;
	for (int i = 0; i < NCART; i++) {
		pos -= param.base_pendulum[i];
	}
	param.theta = atan(pos(Y) / pos(Z));
	param.theta_dot = -ball_vel(Y) / ((param.string_len + param.radius) * cos(param.theta));
}

/**
 * @brief Construct ball state purely from pendulum angle pos and velocity.
 *
 * This method disregards previous computations on the ball state (if any) and constructs the whole
 * ball state from angle pos and velocities.
 */
void calc_ball_from_angle(const ball_params & param, vec3 & ball_pos, vec3 & ball_vel) {

	// we're in the third quadrant
	ball_pos(X) = param.base_pendulum[X];
	ball_pos(Y) = param.base_pendulum[Y] - (param.string_len + param.radius) * sin(param.theta);
	ball_pos(Z) = param.base_pendulum[Z] - (param.string_len + param.radius) * cos(param.theta);
	ball_vel(X) = 0.0;
	ball_vel(Y) = -(param.string_len + param.radius) * param.theta_dot * cos(param.theta);
	ball_vel(Z) = (param.string_len + param.radius) * param.theta_dot * sin(param.theta);
}

/**
 * @brief Modify the incoming basketball velocity and string angle velocity in case of contact.
 *
 * Contact occurs if the norm of the difference between the arm and the ball is less than
 * or equal to the ball radius.
 * The contact model is an IDEAL MOMENTUM EXCHANGE assuming mass_robot >> mass_ball!
 */
void check_for_contact(const vec3 & robot_pos, const vec3 & robot_vel, const vec3 & ball_pos, const bool verbose,
					   ball_params & param, vec3 & ball_vel) {

	if (norm(robot_pos - ball_pos) <= param.radius) { // contact occurs
		if (verbose)
			cout << "Contact between ball and robot!\n";
		// change velocities
		ball_vel = -ball_vel + 2*robot_vel;
		// change string angle velocity
		param.theta_dot = -ball_vel(Y) / ((param.string_len + param.radius) * cos(param.theta));
	}
}


