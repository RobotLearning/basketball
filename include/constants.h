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

/*! endeffector information */
enum RobotEndeffectors {

	RIGHT_HAND,
	LEFT_HAND,
	RIGHT_FOOT,
	LEFT_FOOT,

	N_MAX_ENDEFFECTORS
};

const ivec RIGHT_ARM = {R_SFE, R_SAA, R_HR, R_EB, R_WR, R_WFE, R_WAA};
const ivec LEFT_ARM = {L_SFE, L_SAA, L_HR, L_EB, L_WR, L_WFE, L_WAA};

const int NCART = 3;
const int NQUAT = 4;
const int NARMS_ACTIVE = 2;
const int NDOF = N_MAX_DOFS;
const int NDOF_OPT = 7;
const int NDOF_ACTIVE = NDOF_OPT * NARMS_ACTIVE;
const int NBLOBS = 6;
const int NLINK = N_MAX_DOFS;
const int NENDEFF = N_MAX_ENDEFFECTORS;
const double DT = 0.002; //!< 500 Hz robot operation

/* Ball variables */
const double ball_radius  = 0.1213; //!< standard basketball radius



#endif /* INCLUDE_CONSTANTS_H_ */
