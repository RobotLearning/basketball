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

// indices for simplifying code comprehension
#define X 0
#define Y 1
#define Z 2
#define W 3 //!< for quaternion, not used
#define DX 3
#define DY 4
#define DZ 5

const std::vector<std::string> joint_names = {
			{"L_SFE"}, // 01
			{"L_SAA"}, // 02
			{"L_HR"},  // 03
			{"L_EB"},  // 04
			{"L_WR"},  // 05
			{"L_WFE"}, // 06
			{"L_WAA"}, // 07

			{"R_SFE"}, // 08
			{"R_SAA"}, // 09
			{"R_HR"},  // 10
			{"R_EB"},  // 11
			{"R_WR"},  // 12
			{"R_WFE"}, // 13
			{"R_WAA"}, // 14

			{"L_HFE"}, // 15
			{"L_HAA"}, // 16
			{"L_HFR"}, // 17
			{"L_KFE"}, // 18
			{"L_AR"},  // 19
			{"L_AFE"}, // 20
			{"L_AAA"}, // 21

			{"R_HFE"}, // 22
			{"R_HAA"}, // 23
			{"R_HFR"}, // 24
			{"R_KFE"}, // 25
			{"R_AR"},  // 26
			{"R_AFE"}, // 27
			{"R_AAA"}, // 28

			{"B_TR"},  // 29
			{"B_TAA"}, // 30
			{"B_TFE"}, // 31

			{"B_HN"},  // 32
			{"B_HT"},  // 33
			{"B_HR"},  // 34

			{"R_EP"},  // 35
			{"R_ET"},  // 36
			{"L_EP"},  // 37
			{"L_ET"},  // 38
};

/*! define the DOFs of this robot */
enum RobotDOFs {

	BASE = 0,

	L_SFE, // 01
	L_SAA, // 02
	L_HR,  // 03
	L_EB,  // 04
	L_WR,  // 05
	L_WFE, // 06
	L_WAA, // 07

	R_SFE, // 08
	R_SAA, // 09
	R_HR,  // 10
	R_EB,  // 11
	R_WR,  // 12
	R_WFE, // 13
	R_WAA, // 14

	L_HFE, // 15
	L_HAA, // 16
	L_HFR, // 17
	L_KFE, // 18
	L_AR,  // 19
	L_AFE, // 20
	L_AAA, // 21

	R_HFE, // 22
	R_HAA, // 23
	R_HFR, // 24
	R_KFE, // 25
	R_AR,  // 26
	R_AFE, // 27
	R_AAA, // 28

	B_TR,  // 29
	B_TAA, // 30
	B_TFE, // 31

	B_HN,  // 32
	B_HT,  // 33
	B_HR,  // 34

	R_EP,  // 35
	R_ET,  // 36
	L_EP,  // 37
	L_ET,  // 38

	N_MAX_DOFS
};

/*! endeffector information */
enum RobotEndeffectors {

	RIGHT_HAND = 1,
	LEFT_HAND,
	RIGHT_FOOT,
	LEFT_FOOT,

	N_MAX_ENDEFFECTORS
};

const int NCART = 3;
const int NQUAT = 4;
const int NDOF = N_MAX_DOFS - 1;
const int NDOF_ACTIVE = 7;
const int NBLOBS = 6;
const int NLINK = N_MAX_DOFS-1;
const int NENDEFF = N_MAX_ENDEFFECTORS-1;
const double DT = 0.002; //!< 500 Hz robot operation

/* Ball variables */
const double ball_radius  = 0.1213; //!< standard basketball radius



#endif /* INCLUDE_CONSTANTS_H_ */
