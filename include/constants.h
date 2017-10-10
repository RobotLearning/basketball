/**
 * @file constants.h
 *
 * @brief Constants are located here
 *
 * Constants are indices, number of degrees of freedom,
 * table tennis parameters (gravity, coefficients, radii, ...).
 *
 *  Created on: Feb 2, 2017
 *      Author: okoc
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <string>
#include <vector>
#include <armadillo>

using namespace arma;

// indices for simplifying code comprehension
#define X 0
#define Y 1
#define Z 2
#define W 3 //!< for quaternion, not used
#define DX 3
#define DY 4
#define DZ 5

/*! define the DOFs of this robot */
enum RobotDOFs {

	L_SFE,
	L_SAA,
	L_HR,
	L_EB,
	L_WR,
	L_WFE,
	L_WAA,

	R_SFE,
	R_SAA,
	R_HR,
	R_EB,
	R_WR,
	R_WFE,
	R_WAA,

	L_HFE,
	L_HAA,
	L_HFR,
	L_KFE,
	L_AR,
	L_AFE,
	L_AAA,

	R_HFE,
	R_HAA,
	R_HFR,
	R_KFE,
	R_AR,
	R_AFE,
	R_AAA,

	B_TR,
	B_TAA,
	B_TFE,

	B_HN,
	B_HT,
	B_HR,

	R_EP,
	R_ET,
	L_EP,
	L_ET,

	N_MAX_DOFS
};

/*! links of the robot */
enum RobotLinks {
  B_SACRAL=1,

  L_SHOULDER,
  L_SHOULDER_AA,
  L_SHOULDER_HR,
  L_ELBOW,
  L_WRIST,
  L_HAND,

  R_SHOULDER,
  R_SHOULDER_AA,
  R_SHOULDER_HR,
  R_ELBOW,
  R_WRIST,
  R_HAND,

  B_CERVICAL1,
  B_CERVICAL2,
  R_EYE_AXIS,
  R_EYE,
  L_EYE_AXIS,
  L_EYE,
  B_HEAD,

  R_HIP,
  R_HIP_R,
  R_KNEE,
  R_ANKLE,
  R_OUT_TOE,
  R_IN_TOE,
  R_OUT_HEEL,
  R_IN_HEEL,
  R_FOOT,

  L_HIP,
  L_HIP_R,
  L_KNEE,
  L_ANKLE,
  L_OUT_TOE,
  L_IN_TOE,
  L_OUT_HEEL,
  L_IN_HEEL,
  L_FOOT,

  N_ROBOT_LINKS
};

/** endeffector indices of cart_des_state from SL */
enum RobotEndeffectors {

	RIGHT_HAND,
	LEFT_HAND,
	RIGHT_FOOT,
	LEFT_FOOT,

	N_MAX_ENDEFFECTORS
};

const uvec RIGHT_ARM = {R_SFE, R_SAA, R_HR, R_EB, R_WR, R_WFE, R_WAA};
const uvec LEFT_ARM = {L_SFE, L_SAA, L_HR, L_EB, L_WR, L_WFE, L_WAA};
const uvec active_dofs = join_vert(LEFT_ARM,RIGHT_ARM);
const uvec link2endeffmap = {L_HAND,R_HAND,L_FOOT,R_FOOT};

const int NCART = 3;
const int NQUAT = 4;
const int NARMS_ACTIVE = 2;
const int NDOF = N_MAX_DOFS;
const int NDOF_OPT = 7;
const int NDOF_ACTIVE = NDOF_OPT * NARMS_ACTIVE;
const int NBLOBS = 6;
const int NLINK = N_MAX_DOFS;
const int NENDEFF = N_MAX_ENDEFFECTORS;
const double PI = 3.14159265358979323846;
const double DT = 0.002; //!< 500 Hz robot operation

#endif /* INCLUDE_CONSTANTS_H_ */
