/**
 * @file sl_interface.h
 * @brief Interface to SL. Exposes play() and cheat() modes.
 *
 * Before starting play() mode, relevant options for the Player class
 * must be set by calling set_algorithm().
 *
 */
#ifndef SL_INTERF_H
#define SL_INTERF_H

#ifdef __cplusplus
extern "C" {
#endif

// Interface for Player
extern void play(const SL_Jstate joint_state[], const blob_state *blobs, SL_DJstate joint_des_state[]);
extern void cheat(const SL_Jstate joint_state[], const SL_Cstate *sim_ball_state, SL_DJstate joint_des_state[]);
extern void integrate_ball_state(const double dt, double ball_state[6], double env_vars[6]);
extern void load_options();

#ifdef __cplusplus
} // extern "C"
#endif

#ifdef __cplusplus
// internal c++ functions
static void fuse_blobs(const blob_state * blobs, vec3 & obs);
static void set_optim_type(const int opt_num);
static void check_for_contact(const SL_Cstate robot_state[NENDEFF+1], const double theta,
		                      double ball_state[6], double & theta_dot);

#endif

#endif /* SL_INTERF_H */
