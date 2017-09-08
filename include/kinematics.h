/**
 * @file kinematics.h
 *
 * @brief Here we include the kinematics related functions taken from SL
 *
 *  Created on: Jun 22, 2016
 *      Author: okoc
 */

#ifndef KINEMATICS_H_
#define KINEMATICS_H_

#include "math.h"
#include "constants.h"

#define MAG 1

/*! defines that are used to parse the config and prefs files */
#define MIN_THETA    1
#define MAX_THETA    2
#define THETA_OFFSET 3

/* dimensions of the robot */
#define PELVISOFFSET  0.0042 //!< 0.165in
#define PELVIS2THORAX 0.1761 //!< 6.933in
#define THORAX2HEAD   0.3
#define XHIP          0.0889 //!< 3.5in
#define YHIP          0.0400 //!< 1.574in
#define YKNEE         0.0647 //!< 2.548in
#define UPPERLEG      0.3817 //!< 15.029in
#define LOWERLEG      0.3810 //!< 14.999in
#define THORAX2SHOULDER 0.111 //!< 4.38in
#define SHOULDERX     0.2768 //!< 10.899in
#define SHOULDERY     0.024  //!< 0.95 in
#define UPPERARM      0.2577 //!< 10.146 in
#define LOWERARM      0.2408 //!< 9.481 in
#define WRISTY        0.00376 //!< 0.148 in

#define THORAX2NECK   0.3068
#define CERVICAL      0.0000
#define HEAD          0.1356
#define EYEXOFF       0.0366
#define EYEYOFF       0.0674
#define EYE           0.03
#define TOPofHEAD     0.28

#define FOOT          0.22
#define XTOE          (0.05*MAG)
#define YTOE          (0.2*MAG)
#define ZTOE          0.0512
#define XHEEL         (0.05*MAG)
#define YHEEL         (0.10*MAG)
#define ZHEEL         0.0512

// hands length
#define XHAND		0.136

/* special dimensions of neck arrangement */
#define NECK_A 0.0508 //Ludo 0.0606   /*! the distance between the base hinges of the linear neck actuators */
#define NECK_B 0.0314 //Ludo 0.0329   /*! the maximal moment arm for the head-nod movement */
#define NECK_C 0.1589 //Ludo 0.4649   /*! vertical length from base hinges to universal joint */

// FROM MDEFS.H FILE
#define Power(x, y)	(pow((double)(x), (double)(y)))
#define Sqrt(x)		(sqrt((double)(x)))

#define Abs(x)		(fabs((double)(x)))

#define Exp(x)		(exp((double)(x)))
#define Log(x)		(log((double)(x)))

#define Sin(x)		(sin((double)(x)))
#define Cos(x)		(cos((double)(x)))
#define Tan(x)		(tan((double)(x)))

#define ArcSin(x)       (asin((double)(x)))
#define ArcCos(x)       (acos((double)(x)))
#define ArcTan(x)       (atan((double)(x)))

#define Sinh(x)          (sinh((double)(x)))
#define Cosh(x)          (cosh((double)(x)))
#define Tanh(x)          (tanh((double)(x)))


#ifndef E
#define E		2.71828182845904523536029
#endif
#ifndef Pi
#define Pi		3.14159265358979323846264
#endif
#define Degree		0.01745329251994329576924

// kinematics functions from SL
void get_position(const double q[NDOF], double pos[NCART]);
bool read_joint_limits(double *lb, double *ub);

#endif /* KINEMATICS_H_ */
