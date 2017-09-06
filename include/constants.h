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

// indices for simplifying code comprehension
#define X 0
#define Y 1
#define Z 2
#define W 3 //!< for quaternion, not used
#define DX 3
#define DY 4
#define DZ 5

const int NCART = 3;
const int NDOF = 38;
const int NDOF_ACTIVE = 7;
const int NBLOBS = 1;
const double DT = 0.002; //!< 500 Hz robot operation

/* Ball variables */
const double ball_radius  = 0.1213; //!< standard basketball radius



#endif /* INCLUDE_CONSTANTS_H_ */
