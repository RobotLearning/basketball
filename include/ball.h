/*
 * ball.h
 *
 *  Created on: Sep 26, 2017
 *      Author: okan.koc
 */

#ifndef INCLUDE_BALL_H_
#define INCLUDE_BALL_H_

using namespace arma;

static const bool CHECK_CONTACTS = true; // turn off for simplified debugging

/**
 * @brief Cartesian state of robot hands (LEFT and RIGHT)
 */
struct robot_hands {

	vec3 left_pos = zeros<vec>(NCART);
	vec3 right_pos = zeros<vec>(NCART);
	vec3 left_vel = zeros<vec>(NCART);
	vec3 right_vel = zeros<vec>(NCART);

	robot_hands() {};

	robot_hands(const double lpos[], const double rpos[],
			          const double lvel[], const double rvel[], const int START_IDX) {
		for (int i = 0; i < NCART; i++) {
			left_pos(i) = lpos[i+START_IDX];
			right_pos(i) = rpos[i+START_IDX];
			left_vel(i) = lvel[i+START_IDX];
			right_vel(i) = rvel[i+START_IDX];
		}
	}

	robot_hands(const vec3 & lpos, const vec3 & rpos, const vec3 & lvel, const vec3 & rvel) {
		left_pos = lpos;
		right_pos = rpos;
		left_vel = lvel;
		right_vel = rvel;
	}

	vec print() {
		return join_vert(join_vert(left_pos,right_pos),join_vert(left_vel,right_vel));
	}
};

/**
 * @brief Ball parameters used to simulate & predict future ball path
 * and to calculate desired parameters.
 */
struct ball_params {

	double radius  = 0.1213; //!< standard basketball radius
	double gravity = -9.8;
	double friction = 0.0;
	double restitution = 0.88; //!< coeff. of restitution between ball and robot
	double string_len = 1.0;
	double theta = -PI/4;
	double theta_dot = 0.0;
	double threshold = 0.5; //!< threshold (in seconds) to consider contact

	// desired parameters
	double theta_dot_des = -1.0;
	std::vector<double> base_pendulum = {0.0, 0.9, 1.0};
};

class Ball {

private:
	vec3 pos;
	vec3 vel;
	bool hit = false;
	bool verbose = false;
	ball_params param;
public:

	Ball();
	void load_params(const std::string & file_name_relative);
	void integrate_ball_state(double dt);
	void integrate_ball_state(const robot_hands & robot, double dt);
	vec6 get_state() const;
	void get_state(double ball_state[]) const;
	void set_state(const vec6 & ball_state);
	void set_state(const double ball_state[]);
	void get_env_params(ball_params & par) const;
	void get_env_params(double env_vars[5]) const;
	bool check_for_hit() const;
	void check_for_contact(const vec3 & robot_pos, const vec3 & robot_vel);
};

// ball models
vec calc_next_ball(const vec & xnow, const double dt, const void *fp);
void ball_pendulum_model(const double dt, ball_params & param);
void calc_angle_from_ball(const vec3 & ball_pos, const vec3 & ball_vel, ball_params & param);
void calc_ball_from_angle(const ball_params & param, vec3 & ball_pos, vec3 & ball_vel);

#endif /* INCLUDE_BALL_H_ */
