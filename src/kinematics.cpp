/**
 * @file kinematics.c
 *
 * @brief Here we include the kinematics related functions taken from SL
 *
 * Optimizators call these functions to calculate the Cartesian
 * racket positions, velocities and normals
 *
 *  Created on: Jun 22, 2016
 *      Author: okoc
 */

#include <armadillo>
#include "math.h"
#include "stdlib.h"
#include "constants.h"
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include "kinematics.h"

using namespace arma;

// internal functions used in calculating cartesian quantities
static void kinematics(const double state[NDOF],
		        	   double Xlink[NLINK+1][4],
					   double Xorigin[NDOF+1][4],
					   double Xaxis[NDOF+1][4],
					   double Ahmat[NDOF+1][5][5]);
static void jacobian(const double link[NLINK+1][4],
					 const double origin[NDOF+1][4],
					 const double axis[NDOF+1][4],
					 mat & jac);
static void read_default_state(vec & q_default);


/**
 * @brief Returns the cartesian endeffector positions
 */
void calc_cart_pos(const ivec & active_dofs, const double q_active[], vec & pos_left, vec & pos_right) {

	static double link[NLINK+1][3+1];
	static double origin[NDOF+1][3+1];
	static double axis[NDOF+1][3+1];
	static double amats[NDOF+1][4+1][4+1];
	static vec q_default = zeros<vec>(NDOF);
	static bool firsttime = true;
	thread_local vec q = zeros<vec>(NDOF);

	if (firsttime) {
		read_default_state(q_default);
		firsttime = false;
	}

	q = q_default;
	for (int i = 0; i < active_dofs.n_elem; i++) {
		q(active_dofs[i]) = q_active[i];
	}

	kinematics(q.memptr(),link,origin,axis,amats);
	for (int i = 0; i < NCART; i++) {
		pos_right(i) = link[R_HAND][i+1];
		pos_left(i) = link[L_HAND][i+1];
		//normal[i] = amats[PALM][i+1][2];
	}
}

/**
 * @brief Returns the cartesian positions and velocities of LEFT HAND and RIGHT HAND
 */
void calc_cart_pos_and_vel(const ivec & active_dofs, const double q_active[], const double qdot_active[],
		                   vec3 & pos_left, vec3 & pos_right, vec3 & vel_left, vec3 & vel_right) {

	static double link[NLINK+1][3+1];
	static double origin[NDOF+1][3+1];
	static double axis[NDOF+1][3+1];
	static double amats[NDOF+1][4+1][4+1];
	static vec q_default = zeros<vec>(NDOF);
	static bool firsttime = true;
	thread_local vec q = zeros<vec>(NDOF);
	thread_local mat jac = zeros<mat>(2*NENDEFF*NCART,NDOF);
	thread_local vec qdot = zeros<vec>(NDOF);

	if (firsttime) {
		read_default_state(q_default);
		firsttime = false;
	}

	q = q_default;
	qdot = zeros<vec>(NDOF);
	for (int i = 0; i < active_dofs.n_elem; i++) {
		q(active_dofs[i]) = q_active[i];
		qdot(active_dofs[i]) = qdot_active[i];
	}

	kinematics(q.memptr(),link,origin,axis,amats);

	for (int i = 0; i < NCART; i++) {
		pos_right(i) = link[R_HAND][i+1];
		pos_left(i) = link[L_HAND][i+1];
		//normal[i] = amats[PALM][i+1][2];
	}

	jacobian(link,origin,axis,jac);
	vec vel = jac * qdot;
	//cout << jac;
	//cout << "CART VEL: " << vel.t();
	vel_left = vel(span(0,2));
	vel_right = vel(span(6,8));
}

/*
 *
 * Kinematics from SL
 * TODO: Get rid of 1-based indexing for the 2-5th arguments!
 *
 */
static void kinematics(const double state[NDOF],
		        double Xlink[NLINK+1][4],
				double Xorigin[NDOF+1][4],
				double Xaxis[NDOF+1][4],
		        double Ahmat[NDOF+1][5][5]) {

	static bool firsttime = true;
	static double basec[3+1] = {0.0};
	static double baseo[4+1] = {0.0};
	static double eff_a[NENDEFF+1][NCART+1];
	static double eff_x[NENDEFF+1][NCART+1];

	if (firsttime) {
		firsttime = false;
		// special parameters
		eff_a[RIGHT_HAND+1][1]  = -PI/2.0;
		eff_a[LEFT_HAND+1][1]   = -PI/2.0;
		eff_a[RIGHT_HAND+1][3]  = -PI/2.0;
		eff_a[LEFT_HAND+1][3]   = -PI/2.0;
		eff_x[RIGHT_HAND+1][1]  = XHAND;
		eff_x[LEFT_HAND+1][1]   = XHAND;
		baseo[1] = 1.0;
	}

	static double  sstate29th;
	static double  cstate29th;
	static double  sstate30th;
	static double  cstate30th;
	static double  sstate31th;
	static double  cstate31th;
	static double  sstate1th;
	static double  cstate1th;
	static double  sstate2th;
	static double  cstate2th;
	static double  sstate3th;
	static double  cstate3th;
	static double  sstate4th;
	static double  cstate4th;
	static double  sstate5th;
	static double  cstate5th;
	static double  sstate6th;
	static double  cstate6th;
	static double  sstate7th;
	static double  cstate7th;
	static double  sstate8th;
	static double  cstate8th;
	static double  sstate9th;
	static double  cstate9th;
	static double  sstate10th;
	static double  cstate10th;
	static double  sstate11th;
	static double  cstate11th;
	static double  sstate12th;
	static double  cstate12th;
	static double  sstate13th;
	static double  cstate13th;
	static double  sstate14th;
	static double  cstate14th;
	static double  sstate32th;
	static double  cstate32th;
	static double  sstate33th;
	static double  cstate33th;
	static double  sstate34th;
	static double  cstate34th;
	static double  sstate35th;
	static double  cstate35th;
	static double  sstate36th;
	static double  cstate36th;
	static double  sstate37th;
	static double  cstate37th;
	static double  sstate38th;
	static double  cstate38th;
	static double  sstate23th;
	static double  cstate23th;
	static double  sstate22th;
	static double  cstate22th;
	static double  sstate24th;
	static double  cstate24th;
	static double  sstate25th;
	static double  cstate25th;
	static double  sstate26th;
	static double  cstate26th;
	static double  sstate27th;
	static double  cstate27th;
	static double  sstate28th;
	static double  cstate28th;
	static double  sstate16th;
	static double  cstate16th;
	static double  sstate15th;
	static double  cstate15th;
	static double  sstate17th;
	static double  cstate17th;
	static double  sstate18th;
	static double  cstate18th;
	static double  sstate19th;
	static double  cstate19th;
	static double  sstate20th;
	static double  cstate20th;
	static double  sstate21th;
	static double  cstate21th;

	static double  rseff2a1;
	static double  rceff2a1;
	static double  rseff2a2;
	static double  rceff2a2;
	static double  rseff2a3;
	static double  rceff2a3;
	static double  rseff1a1;
	static double  rceff1a1;
	static double  rseff1a2;
	static double  rceff1a2;
	static double  rseff1a3;
	static double  rceff1a3;
	static double  rseff3a1;
	static double  rceff3a1;
	static double  rseff3a2;
	static double  rceff3a2;
	static double  rseff3a3;
	static double  rceff3a3;
	static double  rseff4a1;
	static double  rceff4a1;
	static double  rseff4a2;
	static double  rceff4a2;
	static double  rseff4a3;
	static double  rceff4a3;

	static double  Hi00[4+1][4+1];
	static double  Hi01[4+1][4+1];
	static double  Hi12[4+1][4+1];
	static double  Hi23[4+1][4+1];
	static double  Hi34[4+1][4+1];
	static double  Hi45[4+1][4+1];
	static double  Hi56[4+1][4+1];
	static double  Hi67[4+1][4+1];
	static double  Hi78[4+1][4+1];
	static double  Hi89[4+1][4+1];
	static double  Hi910[4+1][4+1];
	static double  Hi1011[4+1][4+1];
	static double  Hi312[4+1][4+1];
	static double  Hi1213[4+1][4+1];
	static double  Hi1314[4+1][4+1];
	static double  Hi1415[4+1][4+1];
	static double  Hi1516[4+1][4+1];
	static double  Hi1617[4+1][4+1];
	static double  Hi1718[4+1][4+1];
	static double  Hi1819[4+1][4+1];
	static double  Hi320[4+1][4+1];
	static double  Hi2021[4+1][4+1];
	static double  Hi2122[4+1][4+1];
	static double  Hi2223[4+1][4+1];
	static double  Hi2324[4+1][4+1];
	static double  Hi2425[4+1][4+1];
	static double  Hi2226[4+1][4+1];
	static double  Hi2627[4+1][4+1];
	static double  Hi2728[4+1][4+1];
	static double  Hi2229[4+1][4+1];
	static double  Hi030[4+1][4+1];
	static double  Hi3031[4+1][4+1];
	static double  Hi3132[4+1][4+1];
	static double  Hi3233[4+1][4+1];
	static double  Hi3334[4+1][4+1];
	static double  Hi3435[4+1][4+1];
	static double  Hi3536[4+1][4+1];
	static double  Hi3637[4+1][4+1];
	static double  Hi3638[4+1][4+1];
	static double  Hi3639[4+1][4+1];
	static double  Hi3640[4+1][4+1];
	static double  Hi3641[4+1][4+1];
	static double  Hi042[4+1][4+1];
	static double  Hi4243[4+1][4+1];
	static double  Hi4344[4+1][4+1];
	static double  Hi4445[4+1][4+1];
	static double  Hi4546[4+1][4+1];
	static double  Hi4647[4+1][4+1];
	static double  Hi4748[4+1][4+1];
	static double  Hi4849[4+1][4+1];
	static double  Hi4850[4+1][4+1];
	static double  Hi4851[4+1][4+1];
	static double  Hi4852[4+1][4+1];
	static double  Hi4853[4+1][4+1];

	static double  Ai01[4+1][4+1];
	static double  Ai02[4+1][4+1];
	static double  Ai03[4+1][4+1];
	static double  Ai04[4+1][4+1];
	static double  Ai05[4+1][4+1];
	static double  Ai06[4+1][4+1];
	static double  Ai07[4+1][4+1];
	static double  Ai08[4+1][4+1];
	static double  Ai09[4+1][4+1];
	static double  Ai010[4+1][4+1];
	static double  Ai011[4+1][4+1];
	static double  Ai012[4+1][4+1];
	static double  Ai013[4+1][4+1];
	static double  Ai014[4+1][4+1];
	static double  Ai015[4+1][4+1];
	static double  Ai016[4+1][4+1];
	static double  Ai017[4+1][4+1];
	static double  Ai018[4+1][4+1];
	static double  Ai019[4+1][4+1];
	static double  Ai020[4+1][4+1];
	static double  Ai021[4+1][4+1];
	static double  Ai022[4+1][4+1];
	static double  Ai023[4+1][4+1];
	static double  Ai024[4+1][4+1];
	static double  Ai025[4+1][4+1];
	static double  Ai026[4+1][4+1];
	static double  Ai027[4+1][4+1];
	static double  Ai028[4+1][4+1];
	static double  Ai029[4+1][4+1];
	static double  Ai030[4+1][4+1];
	static double  Ai031[4+1][4+1];
	static double  Ai032[4+1][4+1];
	static double  Ai033[4+1][4+1];
	static double  Ai034[4+1][4+1];
	static double  Ai035[4+1][4+1];
	static double  Ai036[4+1][4+1];
	static double  Ai037[4+1][4+1];
	static double  Ai038[4+1][4+1];
	static double  Ai039[4+1][4+1];
	static double  Ai040[4+1][4+1];
	static double  Ai041[4+1][4+1];
	static double  Ai042[4+1][4+1];
	static double  Ai043[4+1][4+1];
	static double  Ai044[4+1][4+1];
	static double  Ai045[4+1][4+1];
	static double  Ai046[4+1][4+1];
	static double  Ai047[4+1][4+1];
	static double  Ai048[4+1][4+1];
	static double  Ai049[4+1][4+1];
	static double  Ai050[4+1][4+1];
	static double  Ai051[4+1][4+1];
	static double  Ai052[4+1][4+1];
	static double  Ai053[4+1][4+1];

	/* Need [n_joints+1]x[3+1] matrices: Xorigin,Xmcog,Xaxis, and Xlink[nLinks+1][3+1] */

	/* sine and cosine precomputation */

	sstate1th=Sin(state[0]);
	cstate1th=Cos(state[0]);

	sstate2th=Sin(state[1]);
	cstate2th=Cos(state[1]);

	sstate3th=Sin(state[2]);
	cstate3th=Cos(state[2]);

	sstate4th=Sin(state[3]);
	cstate4th=Cos(state[3]);

	sstate5th=Sin(state[4]);
	cstate5th=Cos(state[4]);

	sstate6th=Sin(state[5]);
	cstate6th=Cos(state[5]);

	sstate7th=Sin(state[6]);
	cstate7th=Cos(state[6]);

	sstate8th=Sin(state[7]);
	cstate8th=Cos(state[7]);

	sstate9th=Sin(state[8]);
	cstate9th=Cos(state[8]);

	sstate10th=Sin(state[9]);
	cstate10th=Cos(state[9]);

	sstate11th=Sin(state[10]);
	cstate11th=Cos(state[10]);

	sstate12th=Sin(state[11]);
	cstate12th=Cos(state[11]);

	sstate13th=Sin(state[12]);
	cstate13th=Cos(state[12]);

	sstate14th=Sin(state[13]);
	cstate14th=Cos(state[13]);

	sstate15th=Sin(state[14]);
	cstate15th=Cos(state[14]);

	sstate16th=Sin(state[15]);
	cstate16th=Cos(state[15]);

	sstate17th=Sin(state[16]);
	cstate17th=Cos(state[16]);

	sstate18th=Sin(state[17]);
	cstate18th=Cos(state[17]);

	sstate19th=Sin(state[18]);
	cstate19th=Cos(state[18]);

	sstate20th=Sin(state[19]);
	cstate20th=Cos(state[19]);

	sstate21th=Sin(state[20]);
	cstate21th=Cos(state[20]);

	sstate32th=Sin(state[31]);
	cstate32th=Cos(state[31]);

	sstate33th=Sin(state[32]);
	cstate33th=Cos(state[32]);

	sstate34th=Sin(state[33]);
	cstate34th=Cos(state[33]);

	sstate35th=Sin(state[34]);
	cstate35th=Cos(state[34]);

	sstate36th=Sin(state[35]);
	cstate36th=Cos(state[35]);

	sstate37th=Sin(state[36]);
	cstate37th=Cos(state[36]);

	sstate38th=Sin(state[37]);
	cstate38th=Cos(state[37]);

	sstate23th=Sin(state[22]);
	cstate23th=Cos(state[22]);

	sstate22th=Sin(state[21]);
	cstate22th=Cos(state[21]);

	sstate24th=Sin(state[23]);
	cstate24th=Cos(state[23]);

	sstate25th=Sin(state[24]);
	cstate25th=Cos(state[24]);

	sstate26th=Sin(state[25]);
	cstate26th=Cos(state[25]);

	sstate27th=Sin(state[26]);
	cstate27th=Cos(state[26]);

	sstate28th=Sin(state[27]);
	cstate28th=Cos(state[27]);

	sstate29th=Sin(state[28]);
	cstate29th=Cos(state[28]);

	sstate30th=Sin(state[29]);
	cstate30th=Cos(state[29]);

	sstate31th=Sin(state[30]);
	cstate31th=Cos(state[30]);


	/* rotation matrix sine and cosine precomputation */

	rseff2a1=Sin(eff_a[2][1]);
	rceff2a1=Cos(eff_a[2][1]);

	rseff2a2=Sin(eff_a[2][2]);
	rceff2a2=Cos(eff_a[2][2]);

	rseff2a3=Sin(eff_a[2][3]);
	rceff2a3=Cos(eff_a[2][3]);


	rseff1a1=Sin(eff_a[1][1]);
	rceff1a1=Cos(eff_a[1][1]);

	rseff1a2=Sin(eff_a[1][2]);
	rceff1a2=Cos(eff_a[1][2]);

	rseff1a3=Sin(eff_a[1][3]);
	rceff1a3=Cos(eff_a[1][3]);


	rseff3a1=Sin(eff_a[3][1]);
	rceff3a1=Cos(eff_a[3][1]);

	rseff3a2=Sin(eff_a[3][2]);
	rceff3a2=Cos(eff_a[3][2]);

	rseff3a3=Sin(eff_a[3][3]);
	rceff3a3=Cos(eff_a[3][3]);


	rseff4a1=Sin(eff_a[4][1]);
	rceff4a1=Cos(eff_a[4][1]);

	rseff4a2=Sin(eff_a[4][2]);
	rceff4a2=Cos(eff_a[4][2]);

	rseff4a3=Sin(eff_a[4][3]);
	rceff4a3=Cos(eff_a[4][3]);



	/* inverse homogeneous rotation matrices */
	Hi00[1][1]=-1 + 2*Power(baseo[1],2) + 2*Power(baseo[2],2);
	Hi00[1][2]=2*(baseo[2]*baseo[3] - baseo[1]*baseo[4]);
	Hi00[1][3]=2*(baseo[1]*baseo[3] + baseo[2]*baseo[4]);
	Hi00[1][4]=basec[1];

	Hi00[2][1]=2*(baseo[2]*baseo[3] + baseo[1]*baseo[4]);
	Hi00[2][2]=-1 + 2*Power(baseo[1],2) + 2*Power(baseo[3],2);
	Hi00[2][3]=2*(-(baseo[1]*baseo[2]) + baseo[3]*baseo[4]);
	Hi00[2][4]=basec[2];

	Hi00[3][1]=2*(-(baseo[1]*baseo[3]) + baseo[2]*baseo[4]);
	Hi00[3][2]=2*(baseo[1]*baseo[2] + baseo[3]*baseo[4]);
	Hi00[3][3]=-1 + 2*Power(baseo[1],2) + 2*Power(baseo[4],2);
	Hi00[3][4]=basec[3];


	Hi01[1][1]=cstate29th;
	Hi01[1][2]=-sstate29th;

	Hi01[2][1]=-sstate29th;
	Hi01[2][2]=-cstate29th;
	Hi01[2][4]=-PELVISOFFSET;

	Hi01[3][4]=PELVIS2THORAX;


	Hi12[1][1]=-sstate30th;
	Hi12[1][2]=-cstate30th;

	Hi12[3][1]=cstate30th;
	Hi12[3][2]=-sstate30th;


	Hi23[1][1]=cstate31th;
	Hi23[1][2]=-sstate31th;

	Hi23[3][1]=-sstate31th;
	Hi23[3][2]=-cstate31th;


	Hi34[1][1]=0.7071067811865475*cstate1th;
	Hi34[1][2]=-0.7071067811865475*sstate1th;
	Hi34[1][4]=-THORAX2SHOULDER;

	Hi34[2][1]=-sstate1th;
	Hi34[2][2]=-cstate1th;

	Hi34[3][1]=0.7071067811865475*cstate1th;
	Hi34[3][2]=-0.7071067811865475*sstate1th;


	Hi45[1][1]=-0.7071067811865475*cstate2th - 0.7071067811865475*sstate2th;
	Hi45[1][2]=-0.7071067811865475*cstate2th + 0.7071067811865475*sstate2th;

	Hi45[3][1]=0.7071067811865475*cstate2th - 0.7071067811865475*sstate2th;
	Hi45[3][2]=-0.7071067811865475*cstate2th - 0.7071067811865475*sstate2th;
	Hi45[3][4]=-SHOULDERX;


	Hi56[1][1]=cstate3th;
	Hi56[1][2]=-sstate3th;
	Hi56[1][4]=-SHOULDERY;

	Hi56[3][1]=-sstate3th;
	Hi56[3][2]=-cstate3th;


	Hi67[2][1]=cstate4th;
	Hi67[2][2]=-sstate4th;

	Hi67[3][1]=sstate4th;
	Hi67[3][2]=cstate4th;
	Hi67[3][4]=-UPPERARM;


	Hi78[1][1]=cstate5th;
	Hi78[1][2]=-sstate5th;

	Hi78[3][1]=-sstate5th;
	Hi78[3][2]=-cstate5th;


	Hi89[2][1]=sstate6th;
	Hi89[2][2]=cstate6th;
	Hi89[2][4]=WRISTY;

	Hi89[3][1]=-cstate6th;
	Hi89[3][2]=sstate6th;
	Hi89[3][4]=-LOWERARM;


	Hi910[1][1]=cstate7th;
	Hi910[1][2]=-sstate7th;

	Hi910[3][1]=sstate7th;
	Hi910[3][2]=cstate7th;


	Hi1011[1][1]=rceff2a2*rceff2a3;
	Hi1011[1][2]=-(rceff2a2*rseff2a3);
	Hi1011[1][3]=rseff2a2;
	Hi1011[1][4]=eff_x[2][1];

	Hi1011[2][1]=rceff2a3*rseff2a1*rseff2a2 + rceff2a1*rseff2a3;
	Hi1011[2][2]=rceff2a1*rceff2a3 - rseff2a1*rseff2a2*rseff2a3;
	Hi1011[2][3]=-(rceff2a2*rseff2a1);
	Hi1011[2][4]=eff_x[2][2];

	Hi1011[3][1]=-(rceff2a1*rceff2a3*rseff2a2) + rseff2a1*rseff2a3;
	Hi1011[3][2]=rceff2a3*rseff2a1 + rceff2a1*rseff2a2*rseff2a3;
	Hi1011[3][3]=rceff2a1*rceff2a2;
	Hi1011[3][4]=eff_x[2][3];


	Hi312[1][1]=0.7071067811865475*cstate8th;
	Hi312[1][2]=-0.7071067811865475*sstate8th;
	Hi312[1][4]=-THORAX2SHOULDER;

	Hi312[2][1]=-sstate8th;
	Hi312[2][2]=-cstate8th;

	Hi312[3][1]=-0.7071067811865475*cstate8th;
	Hi312[3][2]=0.7071067811865475*sstate8th;


	Hi1213[1][1]=-0.7071067811865475*cstate9th - 0.7071067811865475*sstate9th;
	Hi1213[1][2]=-0.7071067811865475*cstate9th + 0.7071067811865475*sstate9th;

	Hi1213[3][1]=-0.7071067811865475*cstate9th + 0.7071067811865475*sstate9th;
	Hi1213[3][2]=0.7071067811865475*cstate9th + 0.7071067811865475*sstate9th;
	Hi1213[3][4]=SHOULDERX;


	Hi1314[1][1]=cstate10th;
	Hi1314[1][2]=-sstate10th;
	Hi1314[1][4]=-SHOULDERY;

	Hi1314[3][1]=sstate10th;
	Hi1314[3][2]=cstate10th;


	Hi1415[2][1]=cstate11th;
	Hi1415[2][2]=-sstate11th;

	Hi1415[3][1]=-sstate11th;
	Hi1415[3][2]=-cstate11th;
	Hi1415[3][4]=UPPERARM;


	Hi1516[1][1]=cstate12th;
	Hi1516[1][2]=-sstate12th;

	Hi1516[3][1]=sstate12th;
	Hi1516[3][2]=cstate12th;


	Hi1617[2][1]=sstate13th;
	Hi1617[2][2]=cstate13th;
	Hi1617[2][4]=WRISTY;

	Hi1617[3][1]=cstate13th;
	Hi1617[3][2]=-sstate13th;
	Hi1617[3][4]=LOWERARM;


	Hi1718[1][1]=cstate14th;
	Hi1718[1][2]=-sstate14th;

	Hi1718[3][1]=-sstate14th;
	Hi1718[3][2]=-cstate14th;


	Hi1819[1][1]=rceff1a2*rceff1a3;
	Hi1819[1][2]=-(rceff1a2*rseff1a3);
	Hi1819[1][3]=rseff1a2;
	Hi1819[1][4]=eff_x[1][1];

	Hi1819[2][1]=rceff1a3*rseff1a1*rseff1a2 + rceff1a1*rseff1a3;
	Hi1819[2][2]=rceff1a1*rceff1a3 - rseff1a1*rseff1a2*rseff1a3;
	Hi1819[2][3]=-(rceff1a2*rseff1a1);
	Hi1819[2][4]=eff_x[1][2];

	Hi1819[3][1]=-(rceff1a1*rceff1a3*rseff1a2) + rseff1a1*rseff1a3;
	Hi1819[3][2]=rceff1a3*rseff1a1 + rceff1a1*rseff1a2*rseff1a3;
	Hi1819[3][3]=rceff1a1*rceff1a2;
	Hi1819[3][4]=eff_x[1][3];


	Hi320[1][1]=sstate32th;
	Hi320[1][2]=cstate32th;
	Hi320[1][4]=-THORAX2NECK;

	Hi320[2][1]=-cstate32th;
	Hi320[2][2]=sstate32th;


	Hi2021[2][1]=sstate33th;
	Hi2021[2][2]=cstate33th;
	Hi2021[2][4]=-CERVICAL;

	Hi2021[3][1]=-cstate33th;
	Hi2021[3][2]=sstate33th;


	Hi2122[1][1]=cstate34th;
	Hi2122[1][2]=-sstate34th;

	Hi2122[3][1]=-sstate34th;
	Hi2122[3][2]=-cstate34th;


	Hi2223[1][1]=cstate35th;
	Hi2223[1][2]=-sstate35th;
	Hi2223[1][4]=EYEXOFF;

	Hi2223[2][1]=sstate35th;
	Hi2223[2][2]=cstate35th;
	Hi2223[2][4]=-EYEYOFF;

	Hi2223[3][4]=-HEAD;


	Hi2324[2][1]=sstate36th;
	Hi2324[2][2]=cstate36th;

	Hi2324[3][1]=cstate36th;
	Hi2324[3][2]=-sstate36th;


	Hi2425[2][4]=-EYE;


	Hi2226[1][1]=cstate37th;
	Hi2226[1][2]=-sstate37th;
	Hi2226[1][4]=-EYEXOFF;

	Hi2226[2][1]=sstate37th;
	Hi2226[2][2]=cstate37th;
	Hi2226[2][4]=-EYEYOFF;

	Hi2226[3][4]=-HEAD;


	Hi2627[2][1]=sstate38th;
	Hi2627[2][2]=cstate38th;

	Hi2627[3][1]=cstate38th;
	Hi2627[3][2]=-sstate38th;


	Hi2728[2][4]=-EYE;


	Hi2229[3][4]=-TOPofHEAD;


	Hi030[1][1]=-sstate23th;
	Hi030[1][2]=-cstate23th;
	Hi030[1][4]=XHIP;

	Hi030[3][1]=-cstate23th;
	Hi030[3][2]=sstate23th;


	Hi3031[1][1]=0.26720974913585105*cstate22th + 0.9636383917044586*sstate22th;
	Hi3031[1][2]=0.9636383917044586*cstate22th - 0.26720974913585105*sstate22th;

	Hi3031[3][1]=-0.9636383917044586*cstate22th + 0.26720974913585105*sstate22th;
	Hi3031[3][2]=0.26720974913585105*cstate22th + 0.9636383917044586*sstate22th;


	Hi3132[1][1]=cstate24th;
	Hi3132[1][2]=-sstate24th;
	Hi3132[1][4]=YHIP;

	Hi3132[3][1]=-sstate24th;
	Hi3132[3][2]=-cstate24th;


	Hi3233[1][1]=0.9636374101817434*cstate25th - 0.26721328877550665*sstate25th;
	Hi3233[1][2]=-0.26721328877550665*cstate25th - 0.9636374101817434*sstate25th;
	Hi3233[1][4]=YKNEE;

	Hi3233[3][1]=-0.26721328877550665*cstate25th - 0.9636374101817434*sstate25th;
	Hi3233[3][2]=-0.9636374101817434*cstate25th + 0.26721328877550665*sstate25th;
	Hi3233[3][4]=UPPERLEG;


	Hi3334[1][1]=cstate26th;
	Hi3334[1][2]=-sstate26th;

	Hi3334[3][1]=sstate26th;
	Hi3334[3][2]=cstate26th;


	Hi3435[1][1]=-sstate27th;
	Hi3435[1][2]=-cstate27th;

	Hi3435[3][1]=cstate27th;
	Hi3435[3][2]=-sstate27th;
	Hi3435[3][4]=LOWERLEG;


	Hi3536[1][1]=cstate28th;
	Hi3536[1][2]=-sstate28th;

	Hi3536[3][1]=-sstate28th;
	Hi3536[3][2]=-cstate28th;


	Hi3637[1][4]=ZTOE;

	Hi3637[2][4]=-XTOE;

	Hi3637[3][4]=YTOE;


	Hi3638[1][4]=ZTOE;

	Hi3638[2][4]=XTOE;

	Hi3638[3][4]=YTOE;


	Hi3639[1][4]=ZHEEL;

	Hi3639[2][4]=-XHEEL;

	Hi3639[3][4]=-YHEEL;


	Hi3640[1][4]=ZHEEL;

	Hi3640[2][4]=XHEEL;

	Hi3640[3][4]=-YHEEL;


	Hi3641[1][1]=rceff3a2*rceff3a3;
	Hi3641[1][2]=-(rceff3a2*rseff3a3);
	Hi3641[1][3]=rseff3a2;
	Hi3641[1][4]=eff_x[3][1];

	Hi3641[2][1]=rceff3a3*rseff3a1*rseff3a2 + rceff3a1*rseff3a3;
	Hi3641[2][2]=rceff3a1*rceff3a3 - rseff3a1*rseff3a2*rseff3a3;
	Hi3641[2][3]=-(rceff3a2*rseff3a1);
	Hi3641[2][4]=eff_x[3][2];

	Hi3641[3][1]=-(rceff3a1*rceff3a3*rseff3a2) + rseff3a1*rseff3a3;
	Hi3641[3][2]=rceff3a3*rseff3a1 + rceff3a1*rseff3a2*rseff3a3;
	Hi3641[3][3]=rceff3a1*rceff3a2;
	Hi3641[3][4]=eff_x[3][3];


	Hi042[1][1]=sstate16th;
	Hi042[1][2]=cstate16th;
	Hi042[1][4]=-XHIP;

	Hi042[3][1]=-cstate16th;
	Hi042[3][2]=sstate16th;


	Hi4243[1][1]=0.26720974913585105*cstate15th + 0.9636383917044586*sstate15th;
	Hi4243[1][2]=0.9636383917044586*cstate15th - 0.26720974913585105*sstate15th;

	Hi4243[3][1]=0.9636383917044586*cstate15th - 0.26720974913585105*sstate15th;
	Hi4243[3][2]=-0.26720974913585105*cstate15th - 0.9636383917044586*sstate15th;


	Hi4344[1][1]=cstate17th;
	Hi4344[1][2]=-sstate17th;
	Hi4344[1][4]=YHIP;

	Hi4344[3][1]=sstate17th;
	Hi4344[3][2]=cstate17th;


	Hi4445[1][1]=0.9636374101817434*cstate18th - 0.26721328877550665*sstate18th;
	Hi4445[1][2]=-0.26721328877550665*cstate18th - 0.9636374101817434*sstate18th;
	Hi4445[1][4]=YKNEE;

	Hi4445[3][1]=0.26721328877550665*cstate18th + 0.9636374101817434*sstate18th;
	Hi4445[3][2]=0.9636374101817434*cstate18th - 0.26721328877550665*sstate18th;
	Hi4445[3][4]=-UPPERLEG;


	Hi4546[1][1]=cstate19th;
	Hi4546[1][2]=-sstate19th;

	Hi4546[3][1]=-sstate19th;
	Hi4546[3][2]=-cstate19th;


	Hi4647[1][1]=-sstate20th;
	Hi4647[1][2]=-cstate20th;

	Hi4647[3][1]=-cstate20th;
	Hi4647[3][2]=sstate20th;
	Hi4647[3][4]=-LOWERLEG;


	Hi4748[1][1]=cstate21th;
	Hi4748[1][2]=-sstate21th;

	Hi4748[3][1]=sstate21th;
	Hi4748[3][2]=cstate21th;


	Hi4849[1][4]=ZTOE;

	Hi4849[2][4]=XTOE;

	Hi4849[3][4]=-YTOE;


	Hi4850[1][4]=ZTOE;

	Hi4850[2][4]=-XTOE;

	Hi4850[3][4]=-YTOE;


	Hi4851[1][4]=ZHEEL;

	Hi4851[2][4]=XHEEL;

	Hi4851[3][4]=YHEEL;


	Hi4852[1][4]=ZHEEL;

	Hi4852[2][4]=-XHEEL;

	Hi4852[3][4]=YHEEL;


	Hi4853[1][1]=rceff4a2*rceff4a3;
	Hi4853[1][2]=-(rceff4a2*rseff4a3);
	Hi4853[1][3]=rseff4a2;
	Hi4853[1][4]=eff_x[4][1];

	Hi4853[2][1]=rceff4a3*rseff4a1*rseff4a2 + rceff4a1*rseff4a3;
	Hi4853[2][2]=rceff4a1*rceff4a3 - rseff4a1*rseff4a2*rseff4a3;
	Hi4853[2][3]=-(rceff4a2*rseff4a1);
	Hi4853[2][4]=eff_x[4][2];

	Hi4853[3][1]=-(rceff4a1*rceff4a3*rseff4a2) + rseff4a1*rseff4a3;
	Hi4853[3][2]=rceff4a3*rseff4a1 + rceff4a1*rseff4a2*rseff4a3;
	Hi4853[3][3]=rceff4a1*rceff4a2;
	Hi4853[3][4]=eff_x[4][3];



	/* per link inverse homogeneous rotation matrices */
	Ai01[1][1]=Hi00[1][1]*Hi01[1][1] + Hi00[1][2]*Hi01[2][1];
	Ai01[1][2]=Hi00[1][1]*Hi01[1][2] + Hi00[1][2]*Hi01[2][2];
	Ai01[1][3]=-Hi00[1][3];
	Ai01[1][4]=Hi00[1][4] + Hi00[1][2]*Hi01[2][4] + Hi00[1][3]*Hi01[3][4];

	Ai01[2][1]=Hi00[2][1]*Hi01[1][1] + Hi00[2][2]*Hi01[2][1];
	Ai01[2][2]=Hi00[2][1]*Hi01[1][2] + Hi00[2][2]*Hi01[2][2];
	Ai01[2][3]=-Hi00[2][3];
	Ai01[2][4]=Hi00[2][4] + Hi00[2][2]*Hi01[2][4] + Hi00[2][3]*Hi01[3][4];

	Ai01[3][1]=Hi00[3][1]*Hi01[1][1] + Hi00[3][2]*Hi01[2][1];
	Ai01[3][2]=Hi00[3][1]*Hi01[1][2] + Hi00[3][2]*Hi01[2][2];
	Ai01[3][3]=-Hi00[3][3];
	Ai01[3][4]=Hi00[3][4] + Hi00[3][2]*Hi01[2][4] + Hi00[3][3]*Hi01[3][4];


	Ai02[1][1]=Ai01[1][1]*Hi12[1][1] + Ai01[1][3]*Hi12[3][1];
	Ai02[1][2]=Ai01[1][1]*Hi12[1][2] + Ai01[1][3]*Hi12[3][2];
	Ai02[1][3]=-Ai01[1][2];
	Ai02[1][4]=Ai01[1][4];

	Ai02[2][1]=Ai01[2][1]*Hi12[1][1] + Ai01[2][3]*Hi12[3][1];
	Ai02[2][2]=Ai01[2][1]*Hi12[1][2] + Ai01[2][3]*Hi12[3][2];
	Ai02[2][3]=-Ai01[2][2];
	Ai02[2][4]=Ai01[2][4];

	Ai02[3][1]=Ai01[3][1]*Hi12[1][1] + Ai01[3][3]*Hi12[3][1];
	Ai02[3][2]=Ai01[3][1]*Hi12[1][2] + Ai01[3][3]*Hi12[3][2];
	Ai02[3][3]=-Ai01[3][2];
	Ai02[3][4]=Ai01[3][4];


	Ai03[1][1]=Ai02[1][1]*Hi23[1][1] + Ai02[1][3]*Hi23[3][1];
	Ai03[1][2]=Ai02[1][1]*Hi23[1][2] + Ai02[1][3]*Hi23[3][2];
	Ai03[1][3]=Ai02[1][2];
	Ai03[1][4]=Ai02[1][4];

	Ai03[2][1]=Ai02[2][1]*Hi23[1][1] + Ai02[2][3]*Hi23[3][1];
	Ai03[2][2]=Ai02[2][1]*Hi23[1][2] + Ai02[2][3]*Hi23[3][2];
	Ai03[2][3]=Ai02[2][2];
	Ai03[2][4]=Ai02[2][4];

	Ai03[3][1]=Ai02[3][1]*Hi23[1][1] + Ai02[3][3]*Hi23[3][1];
	Ai03[3][2]=Ai02[3][1]*Hi23[1][2] + Ai02[3][3]*Hi23[3][2];
	Ai03[3][3]=Ai02[3][2];
	Ai03[3][4]=Ai02[3][4];


	Ai04[1][1]=Ai03[1][1]*Hi34[1][1] + Ai03[1][2]*Hi34[2][1] + Ai03[1][3]*Hi34[3][1];
	Ai04[1][2]=Ai03[1][1]*Hi34[1][2] + Ai03[1][2]*Hi34[2][2] + Ai03[1][3]*Hi34[3][2];
	Ai04[1][3]=0.7071067811865475*Ai03[1][1] - 0.7071067811865475*Ai03[1][3];
	Ai04[1][4]=Ai03[1][4] + Ai03[1][1]*Hi34[1][4];

	Ai04[2][1]=Ai03[2][1]*Hi34[1][1] + Ai03[2][2]*Hi34[2][1] + Ai03[2][3]*Hi34[3][1];
	Ai04[2][2]=Ai03[2][1]*Hi34[1][2] + Ai03[2][2]*Hi34[2][2] + Ai03[2][3]*Hi34[3][2];
	Ai04[2][3]=0.7071067811865475*Ai03[2][1] - 0.7071067811865475*Ai03[2][3];
	Ai04[2][4]=Ai03[2][4] + Ai03[2][1]*Hi34[1][4];

	Ai04[3][1]=Ai03[3][1]*Hi34[1][1] + Ai03[3][2]*Hi34[2][1] + Ai03[3][3]*Hi34[3][1];
	Ai04[3][2]=Ai03[3][1]*Hi34[1][2] + Ai03[3][2]*Hi34[2][2] + Ai03[3][3]*Hi34[3][2];
	Ai04[3][3]=0.7071067811865475*Ai03[3][1] - 0.7071067811865475*Ai03[3][3];
	Ai04[3][4]=Ai03[3][4] + Ai03[3][1]*Hi34[1][4];


	Ai05[1][1]=Ai04[1][1]*Hi45[1][1] + Ai04[1][3]*Hi45[3][1];
	Ai05[1][2]=Ai04[1][1]*Hi45[1][2] + Ai04[1][3]*Hi45[3][2];
	Ai05[1][3]=-Ai04[1][2];
	Ai05[1][4]=Ai04[1][4] + Ai04[1][3]*Hi45[3][4];

	Ai05[2][1]=Ai04[2][1]*Hi45[1][1] + Ai04[2][3]*Hi45[3][1];
	Ai05[2][2]=Ai04[2][1]*Hi45[1][2] + Ai04[2][3]*Hi45[3][2];
	Ai05[2][3]=-Ai04[2][2];
	Ai05[2][4]=Ai04[2][4] + Ai04[2][3]*Hi45[3][4];

	Ai05[3][1]=Ai04[3][1]*Hi45[1][1] + Ai04[3][3]*Hi45[3][1];
	Ai05[3][2]=Ai04[3][1]*Hi45[1][2] + Ai04[3][3]*Hi45[3][2];
	Ai05[3][3]=-Ai04[3][2];
	Ai05[3][4]=Ai04[3][4] + Ai04[3][3]*Hi45[3][4];


	Ai06[1][1]=Ai05[1][1]*Hi56[1][1] + Ai05[1][3]*Hi56[3][1];
	Ai06[1][2]=Ai05[1][1]*Hi56[1][2] + Ai05[1][3]*Hi56[3][2];
	Ai06[1][3]=Ai05[1][2];
	Ai06[1][4]=Ai05[1][4] + Ai05[1][1]*Hi56[1][4];

	Ai06[2][1]=Ai05[2][1]*Hi56[1][1] + Ai05[2][3]*Hi56[3][1];
	Ai06[2][2]=Ai05[2][1]*Hi56[1][2] + Ai05[2][3]*Hi56[3][2];
	Ai06[2][3]=Ai05[2][2];
	Ai06[2][4]=Ai05[2][4] + Ai05[2][1]*Hi56[1][4];

	Ai06[3][1]=Ai05[3][1]*Hi56[1][1] + Ai05[3][3]*Hi56[3][1];
	Ai06[3][2]=Ai05[3][1]*Hi56[1][2] + Ai05[3][3]*Hi56[3][2];
	Ai06[3][3]=Ai05[3][2];
	Ai06[3][4]=Ai05[3][4] + Ai05[3][1]*Hi56[1][4];


	Ai07[1][1]=Ai06[1][2]*Hi67[2][1] + Ai06[1][3]*Hi67[3][1];
	Ai07[1][2]=Ai06[1][2]*Hi67[2][2] + Ai06[1][3]*Hi67[3][2];
	Ai07[1][3]=Ai06[1][1];
	Ai07[1][4]=Ai06[1][4] + Ai06[1][3]*Hi67[3][4];

	Ai07[2][1]=Ai06[2][2]*Hi67[2][1] + Ai06[2][3]*Hi67[3][1];
	Ai07[2][2]=Ai06[2][2]*Hi67[2][2] + Ai06[2][3]*Hi67[3][2];
	Ai07[2][3]=Ai06[2][1];
	Ai07[2][4]=Ai06[2][4] + Ai06[2][3]*Hi67[3][4];

	Ai07[3][1]=Ai06[3][2]*Hi67[2][1] + Ai06[3][3]*Hi67[3][1];
	Ai07[3][2]=Ai06[3][2]*Hi67[2][2] + Ai06[3][3]*Hi67[3][2];
	Ai07[3][3]=Ai06[3][1];
	Ai07[3][4]=Ai06[3][4] + Ai06[3][3]*Hi67[3][4];


	Ai08[1][1]=Ai07[1][1]*Hi78[1][1] + Ai07[1][3]*Hi78[3][1];
	Ai08[1][2]=Ai07[1][1]*Hi78[1][2] + Ai07[1][3]*Hi78[3][2];
	Ai08[1][3]=Ai07[1][2];
	Ai08[1][4]=Ai07[1][4];

	Ai08[2][1]=Ai07[2][1]*Hi78[1][1] + Ai07[2][3]*Hi78[3][1];
	Ai08[2][2]=Ai07[2][1]*Hi78[1][2] + Ai07[2][3]*Hi78[3][2];
	Ai08[2][3]=Ai07[2][2];
	Ai08[2][4]=Ai07[2][4];

	Ai08[3][1]=Ai07[3][1]*Hi78[1][1] + Ai07[3][3]*Hi78[3][1];
	Ai08[3][2]=Ai07[3][1]*Hi78[1][2] + Ai07[3][3]*Hi78[3][2];
	Ai08[3][3]=Ai07[3][2];
	Ai08[3][4]=Ai07[3][4];


	Ai09[1][1]=Ai08[1][2]*Hi89[2][1] + Ai08[1][3]*Hi89[3][1];
	Ai09[1][2]=Ai08[1][2]*Hi89[2][2] + Ai08[1][3]*Hi89[3][2];
	Ai09[1][3]=Ai08[1][1];
	Ai09[1][4]=Ai08[1][4] + Ai08[1][2]*Hi89[2][4] + Ai08[1][3]*Hi89[3][4];

	Ai09[2][1]=Ai08[2][2]*Hi89[2][1] + Ai08[2][3]*Hi89[3][1];
	Ai09[2][2]=Ai08[2][2]*Hi89[2][2] + Ai08[2][3]*Hi89[3][2];
	Ai09[2][3]=Ai08[2][1];
	Ai09[2][4]=Ai08[2][4] + Ai08[2][2]*Hi89[2][4] + Ai08[2][3]*Hi89[3][4];

	Ai09[3][1]=Ai08[3][2]*Hi89[2][1] + Ai08[3][3]*Hi89[3][1];
	Ai09[3][2]=Ai08[3][2]*Hi89[2][2] + Ai08[3][3]*Hi89[3][2];
	Ai09[3][3]=Ai08[3][1];
	Ai09[3][4]=Ai08[3][4] + Ai08[3][2]*Hi89[2][4] + Ai08[3][3]*Hi89[3][4];


	Ai010[1][1]=Ai09[1][1]*Hi910[1][1] + Ai09[1][3]*Hi910[3][1];
	Ai010[1][2]=Ai09[1][1]*Hi910[1][2] + Ai09[1][3]*Hi910[3][2];
	Ai010[1][3]=-Ai09[1][2];
	Ai010[1][4]=Ai09[1][4];

	Ai010[2][1]=Ai09[2][1]*Hi910[1][1] + Ai09[2][3]*Hi910[3][1];
	Ai010[2][2]=Ai09[2][1]*Hi910[1][2] + Ai09[2][3]*Hi910[3][2];
	Ai010[2][3]=-Ai09[2][2];
	Ai010[2][4]=Ai09[2][4];

	Ai010[3][1]=Ai09[3][1]*Hi910[1][1] + Ai09[3][3]*Hi910[3][1];
	Ai010[3][2]=Ai09[3][1]*Hi910[1][2] + Ai09[3][3]*Hi910[3][2];
	Ai010[3][3]=-Ai09[3][2];
	Ai010[3][4]=Ai09[3][4];


	Ai011[1][1]=Ai010[1][1]*Hi1011[1][1] + Ai010[1][2]*Hi1011[2][1] + Ai010[1][3]*Hi1011[3][1];
	Ai011[1][2]=Ai010[1][1]*Hi1011[1][2] + Ai010[1][2]*Hi1011[2][2] + Ai010[1][3]*Hi1011[3][2];
	Ai011[1][3]=Ai010[1][1]*Hi1011[1][3] + Ai010[1][2]*Hi1011[2][3] + Ai010[1][3]*Hi1011[3][3];
	Ai011[1][4]=Ai010[1][4] + Ai010[1][1]*Hi1011[1][4] + Ai010[1][2]*Hi1011[2][4] + Ai010[1][3]*Hi1011[3][4];

	Ai011[2][1]=Ai010[2][1]*Hi1011[1][1] + Ai010[2][2]*Hi1011[2][1] + Ai010[2][3]*Hi1011[3][1];
	Ai011[2][2]=Ai010[2][1]*Hi1011[1][2] + Ai010[2][2]*Hi1011[2][2] + Ai010[2][3]*Hi1011[3][2];
	Ai011[2][3]=Ai010[2][1]*Hi1011[1][3] + Ai010[2][2]*Hi1011[2][3] + Ai010[2][3]*Hi1011[3][3];
	Ai011[2][4]=Ai010[2][4] + Ai010[2][1]*Hi1011[1][4] + Ai010[2][2]*Hi1011[2][4] + Ai010[2][3]*Hi1011[3][4];

	Ai011[3][1]=Ai010[3][1]*Hi1011[1][1] + Ai010[3][2]*Hi1011[2][1] + Ai010[3][3]*Hi1011[3][1];
	Ai011[3][2]=Ai010[3][1]*Hi1011[1][2] + Ai010[3][2]*Hi1011[2][2] + Ai010[3][3]*Hi1011[3][2];
	Ai011[3][3]=Ai010[3][1]*Hi1011[1][3] + Ai010[3][2]*Hi1011[2][3] + Ai010[3][3]*Hi1011[3][3];
	Ai011[3][4]=Ai010[3][4] + Ai010[3][1]*Hi1011[1][4] + Ai010[3][2]*Hi1011[2][4] + Ai010[3][3]*Hi1011[3][4];


	Ai012[1][1]=Ai03[1][1]*Hi312[1][1] + Ai03[1][2]*Hi312[2][1] + Ai03[1][3]*Hi312[3][1];
	Ai012[1][2]=Ai03[1][1]*Hi312[1][2] + Ai03[1][2]*Hi312[2][2] + Ai03[1][3]*Hi312[3][2];
	Ai012[1][3]=-0.7071067811865475*Ai03[1][1] - 0.7071067811865475*Ai03[1][3];
	Ai012[1][4]=Ai03[1][4] + Ai03[1][1]*Hi312[1][4];

	Ai012[2][1]=Ai03[2][1]*Hi312[1][1] + Ai03[2][2]*Hi312[2][1] + Ai03[2][3]*Hi312[3][1];
	Ai012[2][2]=Ai03[2][1]*Hi312[1][2] + Ai03[2][2]*Hi312[2][2] + Ai03[2][3]*Hi312[3][2];
	Ai012[2][3]=-0.7071067811865475*Ai03[2][1] - 0.7071067811865475*Ai03[2][3];
	Ai012[2][4]=Ai03[2][4] + Ai03[2][1]*Hi312[1][4];

	Ai012[3][1]=Ai03[3][1]*Hi312[1][1] + Ai03[3][2]*Hi312[2][1] + Ai03[3][3]*Hi312[3][1];
	Ai012[3][2]=Ai03[3][1]*Hi312[1][2] + Ai03[3][2]*Hi312[2][2] + Ai03[3][3]*Hi312[3][2];
	Ai012[3][3]=-0.7071067811865475*Ai03[3][1] - 0.7071067811865475*Ai03[3][3];
	Ai012[3][4]=Ai03[3][4] + Ai03[3][1]*Hi312[1][4];


	Ai013[1][1]=Ai012[1][1]*Hi1213[1][1] + Ai012[1][3]*Hi1213[3][1];
	Ai013[1][2]=Ai012[1][1]*Hi1213[1][2] + Ai012[1][3]*Hi1213[3][2];
	Ai013[1][3]=Ai012[1][2];
	Ai013[1][4]=Ai012[1][4] + Ai012[1][3]*Hi1213[3][4];

	Ai013[2][1]=Ai012[2][1]*Hi1213[1][1] + Ai012[2][3]*Hi1213[3][1];
	Ai013[2][2]=Ai012[2][1]*Hi1213[1][2] + Ai012[2][3]*Hi1213[3][2];
	Ai013[2][3]=Ai012[2][2];
	Ai013[2][4]=Ai012[2][4] + Ai012[2][3]*Hi1213[3][4];

	Ai013[3][1]=Ai012[3][1]*Hi1213[1][1] + Ai012[3][3]*Hi1213[3][1];
	Ai013[3][2]=Ai012[3][1]*Hi1213[1][2] + Ai012[3][3]*Hi1213[3][2];
	Ai013[3][3]=Ai012[3][2];
	Ai013[3][4]=Ai012[3][4] + Ai012[3][3]*Hi1213[3][4];


	Ai014[1][1]=Ai013[1][1]*Hi1314[1][1] + Ai013[1][3]*Hi1314[3][1];
	Ai014[1][2]=Ai013[1][1]*Hi1314[1][2] + Ai013[1][3]*Hi1314[3][2];
	Ai014[1][3]=-Ai013[1][2];
	Ai014[1][4]=Ai013[1][4] + Ai013[1][1]*Hi1314[1][4];

	Ai014[2][1]=Ai013[2][1]*Hi1314[1][1] + Ai013[2][3]*Hi1314[3][1];
	Ai014[2][2]=Ai013[2][1]*Hi1314[1][2] + Ai013[2][3]*Hi1314[3][2];
	Ai014[2][3]=-Ai013[2][2];
	Ai014[2][4]=Ai013[2][4] + Ai013[2][1]*Hi1314[1][4];

	Ai014[3][1]=Ai013[3][1]*Hi1314[1][1] + Ai013[3][3]*Hi1314[3][1];
	Ai014[3][2]=Ai013[3][1]*Hi1314[1][2] + Ai013[3][3]*Hi1314[3][2];
	Ai014[3][3]=-Ai013[3][2];
	Ai014[3][4]=Ai013[3][4] + Ai013[3][1]*Hi1314[1][4];


	Ai015[1][1]=Ai014[1][2]*Hi1415[2][1] + Ai014[1][3]*Hi1415[3][1];
	Ai015[1][2]=Ai014[1][2]*Hi1415[2][2] + Ai014[1][3]*Hi1415[3][2];
	Ai015[1][3]=-Ai014[1][1];
	Ai015[1][4]=Ai014[1][4] + Ai014[1][3]*Hi1415[3][4];

	Ai015[2][1]=Ai014[2][2]*Hi1415[2][1] + Ai014[2][3]*Hi1415[3][1];
	Ai015[2][2]=Ai014[2][2]*Hi1415[2][2] + Ai014[2][3]*Hi1415[3][2];
	Ai015[2][3]=-Ai014[2][1];
	Ai015[2][4]=Ai014[2][4] + Ai014[2][3]*Hi1415[3][4];

	Ai015[3][1]=Ai014[3][2]*Hi1415[2][1] + Ai014[3][3]*Hi1415[3][1];
	Ai015[3][2]=Ai014[3][2]*Hi1415[2][2] + Ai014[3][3]*Hi1415[3][2];
	Ai015[3][3]=-Ai014[3][1];
	Ai015[3][4]=Ai014[3][4] + Ai014[3][3]*Hi1415[3][4];


	Ai016[1][1]=Ai015[1][1]*Hi1516[1][1] + Ai015[1][3]*Hi1516[3][1];
	Ai016[1][2]=Ai015[1][1]*Hi1516[1][2] + Ai015[1][3]*Hi1516[3][2];
	Ai016[1][3]=-Ai015[1][2];
	Ai016[1][4]=Ai015[1][4];

	Ai016[2][1]=Ai015[2][1]*Hi1516[1][1] + Ai015[2][3]*Hi1516[3][1];
	Ai016[2][2]=Ai015[2][1]*Hi1516[1][2] + Ai015[2][3]*Hi1516[3][2];
	Ai016[2][3]=-Ai015[2][2];
	Ai016[2][4]=Ai015[2][4];

	Ai016[3][1]=Ai015[3][1]*Hi1516[1][1] + Ai015[3][3]*Hi1516[3][1];
	Ai016[3][2]=Ai015[3][1]*Hi1516[1][2] + Ai015[3][3]*Hi1516[3][2];
	Ai016[3][3]=-Ai015[3][2];
	Ai016[3][4]=Ai015[3][4];


	Ai017[1][1]=Ai016[1][2]*Hi1617[2][1] + Ai016[1][3]*Hi1617[3][1];
	Ai017[1][2]=Ai016[1][2]*Hi1617[2][2] + Ai016[1][3]*Hi1617[3][2];
	Ai017[1][3]=-Ai016[1][1];
	Ai017[1][4]=Ai016[1][4] + Ai016[1][2]*Hi1617[2][4] + Ai016[1][3]*Hi1617[3][4];

	Ai017[2][1]=Ai016[2][2]*Hi1617[2][1] + Ai016[2][3]*Hi1617[3][1];
	Ai017[2][2]=Ai016[2][2]*Hi1617[2][2] + Ai016[2][3]*Hi1617[3][2];
	Ai017[2][3]=-Ai016[2][1];
	Ai017[2][4]=Ai016[2][4] + Ai016[2][2]*Hi1617[2][4] + Ai016[2][3]*Hi1617[3][4];

	Ai017[3][1]=Ai016[3][2]*Hi1617[2][1] + Ai016[3][3]*Hi1617[3][1];
	Ai017[3][2]=Ai016[3][2]*Hi1617[2][2] + Ai016[3][3]*Hi1617[3][2];
	Ai017[3][3]=-Ai016[3][1];
	Ai017[3][4]=Ai016[3][4] + Ai016[3][2]*Hi1617[2][4] + Ai016[3][3]*Hi1617[3][4];


	Ai018[1][1]=Ai017[1][1]*Hi1718[1][1] + Ai017[1][3]*Hi1718[3][1];
	Ai018[1][2]=Ai017[1][1]*Hi1718[1][2] + Ai017[1][3]*Hi1718[3][2];
	Ai018[1][3]=Ai017[1][2];
	Ai018[1][4]=Ai017[1][4];

	Ai018[2][1]=Ai017[2][1]*Hi1718[1][1] + Ai017[2][3]*Hi1718[3][1];
	Ai018[2][2]=Ai017[2][1]*Hi1718[1][2] + Ai017[2][3]*Hi1718[3][2];
	Ai018[2][3]=Ai017[2][2];
	Ai018[2][4]=Ai017[2][4];

	Ai018[3][1]=Ai017[3][1]*Hi1718[1][1] + Ai017[3][3]*Hi1718[3][1];
	Ai018[3][2]=Ai017[3][1]*Hi1718[1][2] + Ai017[3][3]*Hi1718[3][2];
	Ai018[3][3]=Ai017[3][2];
	Ai018[3][4]=Ai017[3][4];


	Ai019[1][1]=Ai018[1][1]*Hi1819[1][1] + Ai018[1][2]*Hi1819[2][1] + Ai018[1][3]*Hi1819[3][1];
	Ai019[1][2]=Ai018[1][1]*Hi1819[1][2] + Ai018[1][2]*Hi1819[2][2] + Ai018[1][3]*Hi1819[3][2];
	Ai019[1][3]=Ai018[1][1]*Hi1819[1][3] + Ai018[1][2]*Hi1819[2][3] + Ai018[1][3]*Hi1819[3][3];
	Ai019[1][4]=Ai018[1][4] + Ai018[1][1]*Hi1819[1][4] + Ai018[1][2]*Hi1819[2][4] + Ai018[1][3]*Hi1819[3][4];

	Ai019[2][1]=Ai018[2][1]*Hi1819[1][1] + Ai018[2][2]*Hi1819[2][1] + Ai018[2][3]*Hi1819[3][1];
	Ai019[2][2]=Ai018[2][1]*Hi1819[1][2] + Ai018[2][2]*Hi1819[2][2] + Ai018[2][3]*Hi1819[3][2];
	Ai019[2][3]=Ai018[2][1]*Hi1819[1][3] + Ai018[2][2]*Hi1819[2][3] + Ai018[2][3]*Hi1819[3][3];
	Ai019[2][4]=Ai018[2][4] + Ai018[2][1]*Hi1819[1][4] + Ai018[2][2]*Hi1819[2][4] + Ai018[2][3]*Hi1819[3][4];

	Ai019[3][1]=Ai018[3][1]*Hi1819[1][1] + Ai018[3][2]*Hi1819[2][1] + Ai018[3][3]*Hi1819[3][1];
	Ai019[3][2]=Ai018[3][1]*Hi1819[1][2] + Ai018[3][2]*Hi1819[2][2] + Ai018[3][3]*Hi1819[3][2];
	Ai019[3][3]=Ai018[3][1]*Hi1819[1][3] + Ai018[3][2]*Hi1819[2][3] + Ai018[3][3]*Hi1819[3][3];
	Ai019[3][4]=Ai018[3][4] + Ai018[3][1]*Hi1819[1][4] + Ai018[3][2]*Hi1819[2][4] + Ai018[3][3]*Hi1819[3][4];


	Ai020[1][1]=Ai03[1][1]*Hi320[1][1] + Ai03[1][2]*Hi320[2][1];
	Ai020[1][2]=Ai03[1][1]*Hi320[1][2] + Ai03[1][2]*Hi320[2][2];
	Ai020[1][3]=Ai03[1][3];
	Ai020[1][4]=Ai03[1][4] + Ai03[1][1]*Hi320[1][4];

	Ai020[2][1]=Ai03[2][1]*Hi320[1][1] + Ai03[2][2]*Hi320[2][1];
	Ai020[2][2]=Ai03[2][1]*Hi320[1][2] + Ai03[2][2]*Hi320[2][2];
	Ai020[2][3]=Ai03[2][3];
	Ai020[2][4]=Ai03[2][4] + Ai03[2][1]*Hi320[1][4];

	Ai020[3][1]=Ai03[3][1]*Hi320[1][1] + Ai03[3][2]*Hi320[2][1];
	Ai020[3][2]=Ai03[3][1]*Hi320[1][2] + Ai03[3][2]*Hi320[2][2];
	Ai020[3][3]=Ai03[3][3];
	Ai020[3][4]=Ai03[3][4] + Ai03[3][1]*Hi320[1][4];


	Ai021[1][1]=Ai020[1][2]*Hi2021[2][1] + Ai020[1][3]*Hi2021[3][1];
	Ai021[1][2]=Ai020[1][2]*Hi2021[2][2] + Ai020[1][3]*Hi2021[3][2];
	Ai021[1][3]=Ai020[1][1];
	Ai021[1][4]=Ai020[1][4] + Ai020[1][2]*Hi2021[2][4];

	Ai021[2][1]=Ai020[2][2]*Hi2021[2][1] + Ai020[2][3]*Hi2021[3][1];
	Ai021[2][2]=Ai020[2][2]*Hi2021[2][2] + Ai020[2][3]*Hi2021[3][2];
	Ai021[2][3]=Ai020[2][1];
	Ai021[2][4]=Ai020[2][4] + Ai020[2][2]*Hi2021[2][4];

	Ai021[3][1]=Ai020[3][2]*Hi2021[2][1] + Ai020[3][3]*Hi2021[3][1];
	Ai021[3][2]=Ai020[3][2]*Hi2021[2][2] + Ai020[3][3]*Hi2021[3][2];
	Ai021[3][3]=Ai020[3][1];
	Ai021[3][4]=Ai020[3][4] + Ai020[3][2]*Hi2021[2][4];


	Ai022[1][1]=Ai021[1][1]*Hi2122[1][1] + Ai021[1][3]*Hi2122[3][1];
	Ai022[1][2]=Ai021[1][1]*Hi2122[1][2] + Ai021[1][3]*Hi2122[3][2];
	Ai022[1][3]=Ai021[1][2];
	Ai022[1][4]=Ai021[1][4];

	Ai022[2][1]=Ai021[2][1]*Hi2122[1][1] + Ai021[2][3]*Hi2122[3][1];
	Ai022[2][2]=Ai021[2][1]*Hi2122[1][2] + Ai021[2][3]*Hi2122[3][2];
	Ai022[2][3]=Ai021[2][2];
	Ai022[2][4]=Ai021[2][4];

	Ai022[3][1]=Ai021[3][1]*Hi2122[1][1] + Ai021[3][3]*Hi2122[3][1];
	Ai022[3][2]=Ai021[3][1]*Hi2122[1][2] + Ai021[3][3]*Hi2122[3][2];
	Ai022[3][3]=Ai021[3][2];
	Ai022[3][4]=Ai021[3][4];


	Ai023[1][1]=Ai022[1][1]*Hi2223[1][1] + Ai022[1][2]*Hi2223[2][1];
	Ai023[1][2]=Ai022[1][1]*Hi2223[1][2] + Ai022[1][2]*Hi2223[2][2];
	Ai023[1][3]=Ai022[1][3];
	Ai023[1][4]=Ai022[1][4] + Ai022[1][1]*Hi2223[1][4] + Ai022[1][2]*Hi2223[2][4] + Ai022[1][3]*Hi2223[3][4];

	Ai023[2][1]=Ai022[2][1]*Hi2223[1][1] + Ai022[2][2]*Hi2223[2][1];
	Ai023[2][2]=Ai022[2][1]*Hi2223[1][2] + Ai022[2][2]*Hi2223[2][2];
	Ai023[2][3]=Ai022[2][3];
	Ai023[2][4]=Ai022[2][4] + Ai022[2][1]*Hi2223[1][4] + Ai022[2][2]*Hi2223[2][4] + Ai022[2][3]*Hi2223[3][4];

	Ai023[3][1]=Ai022[3][1]*Hi2223[1][1] + Ai022[3][2]*Hi2223[2][1];
	Ai023[3][2]=Ai022[3][1]*Hi2223[1][2] + Ai022[3][2]*Hi2223[2][2];
	Ai023[3][3]=Ai022[3][3];
	Ai023[3][4]=Ai022[3][4] + Ai022[3][1]*Hi2223[1][4] + Ai022[3][2]*Hi2223[2][4] + Ai022[3][3]*Hi2223[3][4];


	Ai024[1][1]=Ai023[1][2]*Hi2324[2][1] + Ai023[1][3]*Hi2324[3][1];
	Ai024[1][2]=Ai023[1][2]*Hi2324[2][2] + Ai023[1][3]*Hi2324[3][2];
	Ai024[1][3]=-Ai023[1][1];
	Ai024[1][4]=Ai023[1][4];

	Ai024[2][1]=Ai023[2][2]*Hi2324[2][1] + Ai023[2][3]*Hi2324[3][1];
	Ai024[2][2]=Ai023[2][2]*Hi2324[2][2] + Ai023[2][3]*Hi2324[3][2];
	Ai024[2][3]=-Ai023[2][1];
	Ai024[2][4]=Ai023[2][4];

	Ai024[3][1]=Ai023[3][2]*Hi2324[2][1] + Ai023[3][3]*Hi2324[3][1];
	Ai024[3][2]=Ai023[3][2]*Hi2324[2][2] + Ai023[3][3]*Hi2324[3][2];
	Ai024[3][3]=-Ai023[3][1];
	Ai024[3][4]=Ai023[3][4];


	Ai025[1][1]=-Ai024[1][3];
	Ai025[1][2]=-Ai024[1][2];
	Ai025[1][3]=-Ai024[1][1];
	Ai025[1][4]=Ai024[1][4] + Ai024[1][2]*Hi2425[2][4];

	Ai025[2][1]=-Ai024[2][3];
	Ai025[2][2]=-Ai024[2][2];
	Ai025[2][3]=-Ai024[2][1];
	Ai025[2][4]=Ai024[2][4] + Ai024[2][2]*Hi2425[2][4];

	Ai025[3][1]=-Ai024[3][3];
	Ai025[3][2]=-Ai024[3][2];
	Ai025[3][3]=-Ai024[3][1];
	Ai025[3][4]=Ai024[3][4] + Ai024[3][2]*Hi2425[2][4];


	Ai026[1][1]=Ai022[1][1]*Hi2226[1][1] + Ai022[1][2]*Hi2226[2][1];
	Ai026[1][2]=Ai022[1][1]*Hi2226[1][2] + Ai022[1][2]*Hi2226[2][2];
	Ai026[1][3]=Ai022[1][3];
	Ai026[1][4]=Ai022[1][4] + Ai022[1][1]*Hi2226[1][4] + Ai022[1][2]*Hi2226[2][4] + Ai022[1][3]*Hi2226[3][4];

	Ai026[2][1]=Ai022[2][1]*Hi2226[1][1] + Ai022[2][2]*Hi2226[2][1];
	Ai026[2][2]=Ai022[2][1]*Hi2226[1][2] + Ai022[2][2]*Hi2226[2][2];
	Ai026[2][3]=Ai022[2][3];
	Ai026[2][4]=Ai022[2][4] + Ai022[2][1]*Hi2226[1][4] + Ai022[2][2]*Hi2226[2][4] + Ai022[2][3]*Hi2226[3][4];

	Ai026[3][1]=Ai022[3][1]*Hi2226[1][1] + Ai022[3][2]*Hi2226[2][1];
	Ai026[3][2]=Ai022[3][1]*Hi2226[1][2] + Ai022[3][2]*Hi2226[2][2];
	Ai026[3][3]=Ai022[3][3];
	Ai026[3][4]=Ai022[3][4] + Ai022[3][1]*Hi2226[1][4] + Ai022[3][2]*Hi2226[2][4] + Ai022[3][3]*Hi2226[3][4];


	Ai027[1][1]=Ai026[1][2]*Hi2627[2][1] + Ai026[1][3]*Hi2627[3][1];
	Ai027[1][2]=Ai026[1][2]*Hi2627[2][2] + Ai026[1][3]*Hi2627[3][2];
	Ai027[1][3]=-Ai026[1][1];
	Ai027[1][4]=Ai026[1][4];

	Ai027[2][1]=Ai026[2][2]*Hi2627[2][1] + Ai026[2][3]*Hi2627[3][1];
	Ai027[2][2]=Ai026[2][2]*Hi2627[2][2] + Ai026[2][3]*Hi2627[3][2];
	Ai027[2][3]=-Ai026[2][1];
	Ai027[2][4]=Ai026[2][4];

	Ai027[3][1]=Ai026[3][2]*Hi2627[2][1] + Ai026[3][3]*Hi2627[3][1];
	Ai027[3][2]=Ai026[3][2]*Hi2627[2][2] + Ai026[3][3]*Hi2627[3][2];
	Ai027[3][3]=-Ai026[3][1];
	Ai027[3][4]=Ai026[3][4];


	Ai028[1][1]=-Ai027[1][3];
	Ai028[1][2]=-Ai027[1][2];
	Ai028[1][3]=-Ai027[1][1];
	Ai028[1][4]=Ai027[1][4] + Ai027[1][2]*Hi2728[2][4];

	Ai028[2][1]=-Ai027[2][3];
	Ai028[2][2]=-Ai027[2][2];
	Ai028[2][3]=-Ai027[2][1];
	Ai028[2][4]=Ai027[2][4] + Ai027[2][2]*Hi2728[2][4];

	Ai028[3][1]=-Ai027[3][3];
	Ai028[3][2]=-Ai027[3][2];
	Ai028[3][3]=-Ai027[3][1];
	Ai028[3][4]=Ai027[3][4] + Ai027[3][2]*Hi2728[2][4];


	Ai029[1][1]=Ai022[1][1];
	Ai029[1][2]=-Ai022[1][2];
	Ai029[1][3]=-Ai022[1][3];
	Ai029[1][4]=Ai022[1][4] + Ai022[1][3]*Hi2229[3][4];

	Ai029[2][1]=Ai022[2][1];
	Ai029[2][2]=-Ai022[2][2];
	Ai029[2][3]=-Ai022[2][3];
	Ai029[2][4]=Ai022[2][4] + Ai022[2][3]*Hi2229[3][4];

	Ai029[3][1]=Ai022[3][1];
	Ai029[3][2]=-Ai022[3][2];
	Ai029[3][3]=-Ai022[3][3];
	Ai029[3][4]=Ai022[3][4] + Ai022[3][3]*Hi2229[3][4];


	Ai030[1][1]=Hi00[1][1]*Hi030[1][1] + Hi00[1][3]*Hi030[3][1];
	Ai030[1][2]=Hi00[1][1]*Hi030[1][2] + Hi00[1][3]*Hi030[3][2];
	Ai030[1][3]=Hi00[1][2];
	Ai030[1][4]=Hi00[1][4] + Hi00[1][1]*Hi030[1][4];

	Ai030[2][1]=Hi00[2][1]*Hi030[1][1] + Hi00[2][3]*Hi030[3][1];
	Ai030[2][2]=Hi00[2][1]*Hi030[1][2] + Hi00[2][3]*Hi030[3][2];
	Ai030[2][3]=Hi00[2][2];
	Ai030[2][4]=Hi00[2][4] + Hi00[2][1]*Hi030[1][4];

	Ai030[3][1]=Hi00[3][1]*Hi030[1][1] + Hi00[3][3]*Hi030[3][1];
	Ai030[3][2]=Hi00[3][1]*Hi030[1][2] + Hi00[3][3]*Hi030[3][2];
	Ai030[3][3]=Hi00[3][2];
	Ai030[3][4]=Hi00[3][4] + Hi00[3][1]*Hi030[1][4];


	Ai031[1][1]=Ai030[1][1]*Hi3031[1][1] + Ai030[1][3]*Hi3031[3][1];
	Ai031[1][2]=Ai030[1][1]*Hi3031[1][2] + Ai030[1][3]*Hi3031[3][2];
	Ai031[1][3]=-Ai030[1][2];
	Ai031[1][4]=Ai030[1][4];

	Ai031[2][1]=Ai030[2][1]*Hi3031[1][1] + Ai030[2][3]*Hi3031[3][1];
	Ai031[2][2]=Ai030[2][1]*Hi3031[1][2] + Ai030[2][3]*Hi3031[3][2];
	Ai031[2][3]=-Ai030[2][2];
	Ai031[2][4]=Ai030[2][4];

	Ai031[3][1]=Ai030[3][1]*Hi3031[1][1] + Ai030[3][3]*Hi3031[3][1];
	Ai031[3][2]=Ai030[3][1]*Hi3031[1][2] + Ai030[3][3]*Hi3031[3][2];
	Ai031[3][3]=-Ai030[3][2];
	Ai031[3][4]=Ai030[3][4];


	Ai032[1][1]=Ai031[1][1]*Hi3132[1][1] + Ai031[1][3]*Hi3132[3][1];
	Ai032[1][2]=Ai031[1][1]*Hi3132[1][2] + Ai031[1][3]*Hi3132[3][2];
	Ai032[1][3]=Ai031[1][2];
	Ai032[1][4]=Ai031[1][4] + Ai031[1][1]*Hi3132[1][4];

	Ai032[2][1]=Ai031[2][1]*Hi3132[1][1] + Ai031[2][3]*Hi3132[3][1];
	Ai032[2][2]=Ai031[2][1]*Hi3132[1][2] + Ai031[2][3]*Hi3132[3][2];
	Ai032[2][3]=Ai031[2][2];
	Ai032[2][4]=Ai031[2][4] + Ai031[2][1]*Hi3132[1][4];

	Ai032[3][1]=Ai031[3][1]*Hi3132[1][1] + Ai031[3][3]*Hi3132[3][1];
	Ai032[3][2]=Ai031[3][1]*Hi3132[1][2] + Ai031[3][3]*Hi3132[3][2];
	Ai032[3][3]=Ai031[3][2];
	Ai032[3][4]=Ai031[3][4] + Ai031[3][1]*Hi3132[1][4];


	Ai033[1][1]=Ai032[1][1]*Hi3233[1][1] + Ai032[1][3]*Hi3233[3][1];
	Ai033[1][2]=Ai032[1][1]*Hi3233[1][2] + Ai032[1][3]*Hi3233[3][2];
	Ai033[1][3]=Ai032[1][2];
	Ai033[1][4]=Ai032[1][4] + Ai032[1][1]*Hi3233[1][4] + Ai032[1][3]*Hi3233[3][4];

	Ai033[2][1]=Ai032[2][1]*Hi3233[1][1] + Ai032[2][3]*Hi3233[3][1];
	Ai033[2][2]=Ai032[2][1]*Hi3233[1][2] + Ai032[2][3]*Hi3233[3][2];
	Ai033[2][3]=Ai032[2][2];
	Ai033[2][4]=Ai032[2][4] + Ai032[2][1]*Hi3233[1][4] + Ai032[2][3]*Hi3233[3][4];

	Ai033[3][1]=Ai032[3][1]*Hi3233[1][1] + Ai032[3][3]*Hi3233[3][1];
	Ai033[3][2]=Ai032[3][1]*Hi3233[1][2] + Ai032[3][3]*Hi3233[3][2];
	Ai033[3][3]=Ai032[3][2];
	Ai033[3][4]=Ai032[3][4] + Ai032[3][1]*Hi3233[1][4] + Ai032[3][3]*Hi3233[3][4];


	Ai034[1][1]=Ai033[1][1]*Hi3334[1][1] + Ai033[1][3]*Hi3334[3][1];
	Ai034[1][2]=Ai033[1][1]*Hi3334[1][2] + Ai033[1][3]*Hi3334[3][2];
	Ai034[1][3]=-Ai033[1][2];
	Ai034[1][4]=Ai033[1][4];

	Ai034[2][1]=Ai033[2][1]*Hi3334[1][1] + Ai033[2][3]*Hi3334[3][1];
	Ai034[2][2]=Ai033[2][1]*Hi3334[1][2] + Ai033[2][3]*Hi3334[3][2];
	Ai034[2][3]=-Ai033[2][2];
	Ai034[2][4]=Ai033[2][4];

	Ai034[3][1]=Ai033[3][1]*Hi3334[1][1] + Ai033[3][3]*Hi3334[3][1];
	Ai034[3][2]=Ai033[3][1]*Hi3334[1][2] + Ai033[3][3]*Hi3334[3][2];
	Ai034[3][3]=-Ai033[3][2];
	Ai034[3][4]=Ai033[3][4];


	Ai035[1][1]=Ai034[1][1]*Hi3435[1][1] + Ai034[1][3]*Hi3435[3][1];
	Ai035[1][2]=Ai034[1][1]*Hi3435[1][2] + Ai034[1][3]*Hi3435[3][2];
	Ai035[1][3]=-Ai034[1][2];
	Ai035[1][4]=Ai034[1][4] + Ai034[1][3]*Hi3435[3][4];

	Ai035[2][1]=Ai034[2][1]*Hi3435[1][1] + Ai034[2][3]*Hi3435[3][1];
	Ai035[2][2]=Ai034[2][1]*Hi3435[1][2] + Ai034[2][3]*Hi3435[3][2];
	Ai035[2][3]=-Ai034[2][2];
	Ai035[2][4]=Ai034[2][4] + Ai034[2][3]*Hi3435[3][4];

	Ai035[3][1]=Ai034[3][1]*Hi3435[1][1] + Ai034[3][3]*Hi3435[3][1];
	Ai035[3][2]=Ai034[3][1]*Hi3435[1][2] + Ai034[3][3]*Hi3435[3][2];
	Ai035[3][3]=-Ai034[3][2];
	Ai035[3][4]=Ai034[3][4] + Ai034[3][3]*Hi3435[3][4];


	Ai036[1][1]=Ai035[1][1]*Hi3536[1][1] + Ai035[1][3]*Hi3536[3][1];
	Ai036[1][2]=Ai035[1][1]*Hi3536[1][2] + Ai035[1][3]*Hi3536[3][2];
	Ai036[1][3]=Ai035[1][2];
	Ai036[1][4]=Ai035[1][4];

	Ai036[2][1]=Ai035[2][1]*Hi3536[1][1] + Ai035[2][3]*Hi3536[3][1];
	Ai036[2][2]=Ai035[2][1]*Hi3536[1][2] + Ai035[2][3]*Hi3536[3][2];
	Ai036[2][3]=Ai035[2][2];
	Ai036[2][4]=Ai035[2][4];

	Ai036[3][1]=Ai035[3][1]*Hi3536[1][1] + Ai035[3][3]*Hi3536[3][1];
	Ai036[3][2]=Ai035[3][1]*Hi3536[1][2] + Ai035[3][3]*Hi3536[3][2];
	Ai036[3][3]=Ai035[3][2];
	Ai036[3][4]=Ai035[3][4];


	Ai037[1][1]=Ai036[1][1];
	Ai037[1][2]=Ai036[1][2];
	Ai037[1][3]=Ai036[1][3];
	Ai037[1][4]=Ai036[1][4] + Ai036[1][1]*Hi3637[1][4] + Ai036[1][2]*Hi3637[2][4] + Ai036[1][3]*Hi3637[3][4];

	Ai037[2][1]=Ai036[2][1];
	Ai037[2][2]=Ai036[2][2];
	Ai037[2][3]=Ai036[2][3];
	Ai037[2][4]=Ai036[2][4] + Ai036[2][1]*Hi3637[1][4] + Ai036[2][2]*Hi3637[2][4] + Ai036[2][3]*Hi3637[3][4];

	Ai037[3][1]=Ai036[3][1];
	Ai037[3][2]=Ai036[3][2];
	Ai037[3][3]=Ai036[3][3];
	Ai037[3][4]=Ai036[3][4] + Ai036[3][1]*Hi3637[1][4] + Ai036[3][2]*Hi3637[2][4] + Ai036[3][3]*Hi3637[3][4];


	Ai038[1][1]=Ai036[1][1];
	Ai038[1][2]=Ai036[1][2];
	Ai038[1][3]=Ai036[1][3];
	Ai038[1][4]=Ai036[1][4] + Ai036[1][1]*Hi3638[1][4] + Ai036[1][2]*Hi3638[2][4] + Ai036[1][3]*Hi3638[3][4];

	Ai038[2][1]=Ai036[2][1];
	Ai038[2][2]=Ai036[2][2];
	Ai038[2][3]=Ai036[2][3];
	Ai038[2][4]=Ai036[2][4] + Ai036[2][1]*Hi3638[1][4] + Ai036[2][2]*Hi3638[2][4] + Ai036[2][3]*Hi3638[3][4];

	Ai038[3][1]=Ai036[3][1];
	Ai038[3][2]=Ai036[3][2];
	Ai038[3][3]=Ai036[3][3];
	Ai038[3][4]=Ai036[3][4] + Ai036[3][1]*Hi3638[1][4] + Ai036[3][2]*Hi3638[2][4] + Ai036[3][3]*Hi3638[3][4];


	Ai039[1][1]=Ai036[1][1];
	Ai039[1][2]=Ai036[1][2];
	Ai039[1][3]=Ai036[1][3];
	Ai039[1][4]=Ai036[1][4] + Ai036[1][1]*Hi3639[1][4] + Ai036[1][2]*Hi3639[2][4] + Ai036[1][3]*Hi3639[3][4];

	Ai039[2][1]=Ai036[2][1];
	Ai039[2][2]=Ai036[2][2];
	Ai039[2][3]=Ai036[2][3];
	Ai039[2][4]=Ai036[2][4] + Ai036[2][1]*Hi3639[1][4] + Ai036[2][2]*Hi3639[2][4] + Ai036[2][3]*Hi3639[3][4];

	Ai039[3][1]=Ai036[3][1];
	Ai039[3][2]=Ai036[3][2];
	Ai039[3][3]=Ai036[3][3];
	Ai039[3][4]=Ai036[3][4] + Ai036[3][1]*Hi3639[1][4] + Ai036[3][2]*Hi3639[2][4] + Ai036[3][3]*Hi3639[3][4];


	Ai040[1][1]=Ai036[1][1];
	Ai040[1][2]=Ai036[1][2];
	Ai040[1][3]=Ai036[1][3];
	Ai040[1][4]=Ai036[1][4] + Ai036[1][1]*Hi3640[1][4] + Ai036[1][2]*Hi3640[2][4] + Ai036[1][3]*Hi3640[3][4];

	Ai040[2][1]=Ai036[2][1];
	Ai040[2][2]=Ai036[2][2];
	Ai040[2][3]=Ai036[2][3];
	Ai040[2][4]=Ai036[2][4] + Ai036[2][1]*Hi3640[1][4] + Ai036[2][2]*Hi3640[2][4] + Ai036[2][3]*Hi3640[3][4];

	Ai040[3][1]=Ai036[3][1];
	Ai040[3][2]=Ai036[3][2];
	Ai040[3][3]=Ai036[3][3];
	Ai040[3][4]=Ai036[3][4] + Ai036[3][1]*Hi3640[1][4] + Ai036[3][2]*Hi3640[2][4] + Ai036[3][3]*Hi3640[3][4];


	Ai041[1][1]=Ai036[1][1]*Hi3641[1][1] + Ai036[1][2]*Hi3641[2][1] + Ai036[1][3]*Hi3641[3][1];
	Ai041[1][2]=Ai036[1][1]*Hi3641[1][2] + Ai036[1][2]*Hi3641[2][2] + Ai036[1][3]*Hi3641[3][2];
	Ai041[1][3]=Ai036[1][1]*Hi3641[1][3] + Ai036[1][2]*Hi3641[2][3] + Ai036[1][3]*Hi3641[3][3];
	Ai041[1][4]=Ai036[1][4] + Ai036[1][1]*Hi3641[1][4] + Ai036[1][2]*Hi3641[2][4] + Ai036[1][3]*Hi3641[3][4];

	Ai041[2][1]=Ai036[2][1]*Hi3641[1][1] + Ai036[2][2]*Hi3641[2][1] + Ai036[2][3]*Hi3641[3][1];
	Ai041[2][2]=Ai036[2][1]*Hi3641[1][2] + Ai036[2][2]*Hi3641[2][2] + Ai036[2][3]*Hi3641[3][2];
	Ai041[2][3]=Ai036[2][1]*Hi3641[1][3] + Ai036[2][2]*Hi3641[2][3] + Ai036[2][3]*Hi3641[3][3];
	Ai041[2][4]=Ai036[2][4] + Ai036[2][1]*Hi3641[1][4] + Ai036[2][2]*Hi3641[2][4] + Ai036[2][3]*Hi3641[3][4];

	Ai041[3][1]=Ai036[3][1]*Hi3641[1][1] + Ai036[3][2]*Hi3641[2][1] + Ai036[3][3]*Hi3641[3][1];
	Ai041[3][2]=Ai036[3][1]*Hi3641[1][2] + Ai036[3][2]*Hi3641[2][2] + Ai036[3][3]*Hi3641[3][2];
	Ai041[3][3]=Ai036[3][1]*Hi3641[1][3] + Ai036[3][2]*Hi3641[2][3] + Ai036[3][3]*Hi3641[3][3];
	Ai041[3][4]=Ai036[3][4] + Ai036[3][1]*Hi3641[1][4] + Ai036[3][2]*Hi3641[2][4] + Ai036[3][3]*Hi3641[3][4];


	Ai042[1][1]=Hi00[1][1]*Hi042[1][1] + Hi00[1][3]*Hi042[3][1];
	Ai042[1][2]=Hi00[1][1]*Hi042[1][2] + Hi00[1][3]*Hi042[3][2];
	Ai042[1][3]=-Hi00[1][2];
	Ai042[1][4]=Hi00[1][4] + Hi00[1][1]*Hi042[1][4];

	Ai042[2][1]=Hi00[2][1]*Hi042[1][1] + Hi00[2][3]*Hi042[3][1];
	Ai042[2][2]=Hi00[2][1]*Hi042[1][2] + Hi00[2][3]*Hi042[3][2];
	Ai042[2][3]=-Hi00[2][2];
	Ai042[2][4]=Hi00[2][4] + Hi00[2][1]*Hi042[1][4];

	Ai042[3][1]=Hi00[3][1]*Hi042[1][1] + Hi00[3][3]*Hi042[3][1];
	Ai042[3][2]=Hi00[3][1]*Hi042[1][2] + Hi00[3][3]*Hi042[3][2];
	Ai042[3][3]=-Hi00[3][2];
	Ai042[3][4]=Hi00[3][4] + Hi00[3][1]*Hi042[1][4];


	Ai043[1][1]=Ai042[1][1]*Hi4243[1][1] + Ai042[1][3]*Hi4243[3][1];
	Ai043[1][2]=Ai042[1][1]*Hi4243[1][2] + Ai042[1][3]*Hi4243[3][2];
	Ai043[1][3]=Ai042[1][2];
	Ai043[1][4]=Ai042[1][4];

	Ai043[2][1]=Ai042[2][1]*Hi4243[1][1] + Ai042[2][3]*Hi4243[3][1];
	Ai043[2][2]=Ai042[2][1]*Hi4243[1][2] + Ai042[2][3]*Hi4243[3][2];
	Ai043[2][3]=Ai042[2][2];
	Ai043[2][4]=Ai042[2][4];

	Ai043[3][1]=Ai042[3][1]*Hi4243[1][1] + Ai042[3][3]*Hi4243[3][1];
	Ai043[3][2]=Ai042[3][1]*Hi4243[1][2] + Ai042[3][3]*Hi4243[3][2];
	Ai043[3][3]=Ai042[3][2];
	Ai043[3][4]=Ai042[3][4];


	Ai044[1][1]=Ai043[1][1]*Hi4344[1][1] + Ai043[1][3]*Hi4344[3][1];
	Ai044[1][2]=Ai043[1][1]*Hi4344[1][2] + Ai043[1][3]*Hi4344[3][2];
	Ai044[1][3]=-Ai043[1][2];
	Ai044[1][4]=Ai043[1][4] + Ai043[1][1]*Hi4344[1][4];

	Ai044[2][1]=Ai043[2][1]*Hi4344[1][1] + Ai043[2][3]*Hi4344[3][1];
	Ai044[2][2]=Ai043[2][1]*Hi4344[1][2] + Ai043[2][3]*Hi4344[3][2];
	Ai044[2][3]=-Ai043[2][2];
	Ai044[2][4]=Ai043[2][4] + Ai043[2][1]*Hi4344[1][4];

	Ai044[3][1]=Ai043[3][1]*Hi4344[1][1] + Ai043[3][3]*Hi4344[3][1];
	Ai044[3][2]=Ai043[3][1]*Hi4344[1][2] + Ai043[3][3]*Hi4344[3][2];
	Ai044[3][3]=-Ai043[3][2];
	Ai044[3][4]=Ai043[3][4] + Ai043[3][1]*Hi4344[1][4];


	Ai045[1][1]=Ai044[1][1]*Hi4445[1][1] + Ai044[1][3]*Hi4445[3][1];
	Ai045[1][2]=Ai044[1][1]*Hi4445[1][2] + Ai044[1][3]*Hi4445[3][2];
	Ai045[1][3]=-Ai044[1][2];
	Ai045[1][4]=Ai044[1][4] + Ai044[1][1]*Hi4445[1][4] + Ai044[1][3]*Hi4445[3][4];

	Ai045[2][1]=Ai044[2][1]*Hi4445[1][1] + Ai044[2][3]*Hi4445[3][1];
	Ai045[2][2]=Ai044[2][1]*Hi4445[1][2] + Ai044[2][3]*Hi4445[3][2];
	Ai045[2][3]=-Ai044[2][2];
	Ai045[2][4]=Ai044[2][4] + Ai044[2][1]*Hi4445[1][4] + Ai044[2][3]*Hi4445[3][4];

	Ai045[3][1]=Ai044[3][1]*Hi4445[1][1] + Ai044[3][3]*Hi4445[3][1];
	Ai045[3][2]=Ai044[3][1]*Hi4445[1][2] + Ai044[3][3]*Hi4445[3][2];
	Ai045[3][3]=-Ai044[3][2];
	Ai045[3][4]=Ai044[3][4] + Ai044[3][1]*Hi4445[1][4] + Ai044[3][3]*Hi4445[3][4];


	Ai046[1][1]=Ai045[1][1]*Hi4546[1][1] + Ai045[1][3]*Hi4546[3][1];
	Ai046[1][2]=Ai045[1][1]*Hi4546[1][2] + Ai045[1][3]*Hi4546[3][2];
	Ai046[1][3]=Ai045[1][2];
	Ai046[1][4]=Ai045[1][4];

	Ai046[2][1]=Ai045[2][1]*Hi4546[1][1] + Ai045[2][3]*Hi4546[3][1];
	Ai046[2][2]=Ai045[2][1]*Hi4546[1][2] + Ai045[2][3]*Hi4546[3][2];
	Ai046[2][3]=Ai045[2][2];
	Ai046[2][4]=Ai045[2][4];

	Ai046[3][1]=Ai045[3][1]*Hi4546[1][1] + Ai045[3][3]*Hi4546[3][1];
	Ai046[3][2]=Ai045[3][1]*Hi4546[1][2] + Ai045[3][3]*Hi4546[3][2];
	Ai046[3][3]=Ai045[3][2];
	Ai046[3][4]=Ai045[3][4];


	Ai047[1][1]=Ai046[1][1]*Hi4647[1][1] + Ai046[1][3]*Hi4647[3][1];
	Ai047[1][2]=Ai046[1][1]*Hi4647[1][2] + Ai046[1][3]*Hi4647[3][2];
	Ai047[1][3]=Ai046[1][2];
	Ai047[1][4]=Ai046[1][4] + Ai046[1][3]*Hi4647[3][4];

	Ai047[2][1]=Ai046[2][1]*Hi4647[1][1] + Ai046[2][3]*Hi4647[3][1];
	Ai047[2][2]=Ai046[2][1]*Hi4647[1][2] + Ai046[2][3]*Hi4647[3][2];
	Ai047[2][3]=Ai046[2][2];
	Ai047[2][4]=Ai046[2][4] + Ai046[2][3]*Hi4647[3][4];

	Ai047[3][1]=Ai046[3][1]*Hi4647[1][1] + Ai046[3][3]*Hi4647[3][1];
	Ai047[3][2]=Ai046[3][1]*Hi4647[1][2] + Ai046[3][3]*Hi4647[3][2];
	Ai047[3][3]=Ai046[3][2];
	Ai047[3][4]=Ai046[3][4] + Ai046[3][3]*Hi4647[3][4];


	Ai048[1][1]=Ai047[1][1]*Hi4748[1][1] + Ai047[1][3]*Hi4748[3][1];
	Ai048[1][2]=Ai047[1][1]*Hi4748[1][2] + Ai047[1][3]*Hi4748[3][2];
	Ai048[1][3]=-Ai047[1][2];
	Ai048[1][4]=Ai047[1][4];

	Ai048[2][1]=Ai047[2][1]*Hi4748[1][1] + Ai047[2][3]*Hi4748[3][1];
	Ai048[2][2]=Ai047[2][1]*Hi4748[1][2] + Ai047[2][3]*Hi4748[3][2];
	Ai048[2][3]=-Ai047[2][2];
	Ai048[2][4]=Ai047[2][4];

	Ai048[3][1]=Ai047[3][1]*Hi4748[1][1] + Ai047[3][3]*Hi4748[3][1];
	Ai048[3][2]=Ai047[3][1]*Hi4748[1][2] + Ai047[3][3]*Hi4748[3][2];
	Ai048[3][3]=-Ai047[3][2];
	Ai048[3][4]=Ai047[3][4];


	Ai049[1][1]=Ai048[1][1];
	Ai049[1][2]=Ai048[1][2];
	Ai049[1][3]=Ai048[1][3];
	Ai049[1][4]=Ai048[1][4] + Ai048[1][1]*Hi4849[1][4] + Ai048[1][2]*Hi4849[2][4] + Ai048[1][3]*Hi4849[3][4];

	Ai049[2][1]=Ai048[2][1];
	Ai049[2][2]=Ai048[2][2];
	Ai049[2][3]=Ai048[2][3];
	Ai049[2][4]=Ai048[2][4] + Ai048[2][1]*Hi4849[1][4] + Ai048[2][2]*Hi4849[2][4] + Ai048[2][3]*Hi4849[3][4];

	Ai049[3][1]=Ai048[3][1];
	Ai049[3][2]=Ai048[3][2];
	Ai049[3][3]=Ai048[3][3];
	Ai049[3][4]=Ai048[3][4] + Ai048[3][1]*Hi4849[1][4] + Ai048[3][2]*Hi4849[2][4] + Ai048[3][3]*Hi4849[3][4];


	Ai050[1][1]=Ai048[1][1];
	Ai050[1][2]=Ai048[1][2];
	Ai050[1][3]=Ai048[1][3];
	Ai050[1][4]=Ai048[1][4] + Ai048[1][1]*Hi4850[1][4] + Ai048[1][2]*Hi4850[2][4] + Ai048[1][3]*Hi4850[3][4];

	Ai050[2][1]=Ai048[2][1];
	Ai050[2][2]=Ai048[2][2];
	Ai050[2][3]=Ai048[2][3];
	Ai050[2][4]=Ai048[2][4] + Ai048[2][1]*Hi4850[1][4] + Ai048[2][2]*Hi4850[2][4] + Ai048[2][3]*Hi4850[3][4];

	Ai050[3][1]=Ai048[3][1];
	Ai050[3][2]=Ai048[3][2];
	Ai050[3][3]=Ai048[3][3];
	Ai050[3][4]=Ai048[3][4] + Ai048[3][1]*Hi4850[1][4] + Ai048[3][2]*Hi4850[2][4] + Ai048[3][3]*Hi4850[3][4];


	Ai051[1][1]=Ai048[1][1];
	Ai051[1][2]=Ai048[1][2];
	Ai051[1][3]=Ai048[1][3];
	Ai051[1][4]=Ai048[1][4] + Ai048[1][1]*Hi4851[1][4] + Ai048[1][2]*Hi4851[2][4] + Ai048[1][3]*Hi4851[3][4];

	Ai051[2][1]=Ai048[2][1];
	Ai051[2][2]=Ai048[2][2];
	Ai051[2][3]=Ai048[2][3];
	Ai051[2][4]=Ai048[2][4] + Ai048[2][1]*Hi4851[1][4] + Ai048[2][2]*Hi4851[2][4] + Ai048[2][3]*Hi4851[3][4];

	Ai051[3][1]=Ai048[3][1];
	Ai051[3][2]=Ai048[3][2];
	Ai051[3][3]=Ai048[3][3];
	Ai051[3][4]=Ai048[3][4] + Ai048[3][1]*Hi4851[1][4] + Ai048[3][2]*Hi4851[2][4] + Ai048[3][3]*Hi4851[3][4];


	Ai052[1][1]=Ai048[1][1];
	Ai052[1][2]=Ai048[1][2];
	Ai052[1][3]=Ai048[1][3];
	Ai052[1][4]=Ai048[1][4] + Ai048[1][1]*Hi4852[1][4] + Ai048[1][2]*Hi4852[2][4] + Ai048[1][3]*Hi4852[3][4];

	Ai052[2][1]=Ai048[2][1];
	Ai052[2][2]=Ai048[2][2];
	Ai052[2][3]=Ai048[2][3];
	Ai052[2][4]=Ai048[2][4] + Ai048[2][1]*Hi4852[1][4] + Ai048[2][2]*Hi4852[2][4] + Ai048[2][3]*Hi4852[3][4];

	Ai052[3][1]=Ai048[3][1];
	Ai052[3][2]=Ai048[3][2];
	Ai052[3][3]=Ai048[3][3];
	Ai052[3][4]=Ai048[3][4] + Ai048[3][1]*Hi4852[1][4] + Ai048[3][2]*Hi4852[2][4] + Ai048[3][3]*Hi4852[3][4];


	Ai053[1][1]=Ai048[1][1]*Hi4853[1][1] + Ai048[1][2]*Hi4853[2][1] + Ai048[1][3]*Hi4853[3][1];
	Ai053[1][2]=Ai048[1][1]*Hi4853[1][2] + Ai048[1][2]*Hi4853[2][2] + Ai048[1][3]*Hi4853[3][2];
	Ai053[1][3]=Ai048[1][1]*Hi4853[1][3] + Ai048[1][2]*Hi4853[2][3] + Ai048[1][3]*Hi4853[3][3];
	Ai053[1][4]=Ai048[1][4] + Ai048[1][1]*Hi4853[1][4] + Ai048[1][2]*Hi4853[2][4] + Ai048[1][3]*Hi4853[3][4];

	Ai053[2][1]=Ai048[2][1]*Hi4853[1][1] + Ai048[2][2]*Hi4853[2][1] + Ai048[2][3]*Hi4853[3][1];
	Ai053[2][2]=Ai048[2][1]*Hi4853[1][2] + Ai048[2][2]*Hi4853[2][2] + Ai048[2][3]*Hi4853[3][2];
	Ai053[2][3]=Ai048[2][1]*Hi4853[1][3] + Ai048[2][2]*Hi4853[2][3] + Ai048[2][3]*Hi4853[3][3];
	Ai053[2][4]=Ai048[2][4] + Ai048[2][1]*Hi4853[1][4] + Ai048[2][2]*Hi4853[2][4] + Ai048[2][3]*Hi4853[3][4];

	Ai053[3][1]=Ai048[3][1]*Hi4853[1][1] + Ai048[3][2]*Hi4853[2][1] + Ai048[3][3]*Hi4853[3][1];
	Ai053[3][2]=Ai048[3][1]*Hi4853[1][2] + Ai048[3][2]*Hi4853[2][2] + Ai048[3][3]*Hi4853[3][2];
	Ai053[3][3]=Ai048[3][1]*Hi4853[1][3] + Ai048[3][2]*Hi4853[2][3] + Ai048[3][3]*Hi4853[3][3];
	Ai053[3][4]=Ai048[3][4] + Ai048[3][1]*Hi4853[1][4] + Ai048[3][2]*Hi4853[2][4] + Ai048[3][3]*Hi4853[3][4];



	/* joint ID: 0 */
	Xorigin[0][1]=Hi00[1][4];
	Xorigin[0][2]=Hi00[2][4];
	Xorigin[0][3]=Hi00[3][4];

	/* link: {basec$0$$x[[1]], basec$0$$x[[2]], basec$0$$x[[3]]} */
	Xlink[0][1]=Hi00[1][4];
	Xlink[0][2]=Hi00[2][4];
	Xlink[0][3]=Hi00[3][4];

	Ahmat[0][1][1]=Hi00[1][1];
	Ahmat[0][1][2]=Hi00[1][2];
	Ahmat[0][1][3]=Hi00[1][3];
	Ahmat[0][1][4]=Hi00[1][4];

	Ahmat[0][2][1]=Hi00[2][1];
	Ahmat[0][2][2]=Hi00[2][2];
	Ahmat[0][2][3]=Hi00[2][3];
	Ahmat[0][2][4]=Hi00[2][4];

	Ahmat[0][3][1]=Hi00[3][1];
	Ahmat[0][3][2]=Hi00[3][2];
	Ahmat[0][3][3]=Hi00[3][3];
	Ahmat[0][3][4]=Hi00[3][4];

	Ahmat[0][4][4]=1;


	/* joint ID: 29 */
	Xorigin[29][1]=Ai01[1][4];
	Xorigin[29][2]=Ai01[2][4];
	Xorigin[29][3]=Ai01[3][4];

	Xaxis[29][1]=Ai01[1][3];
	Xaxis[29][2]=Ai01[2][3];
	Xaxis[29][3]=Ai01[3][3];

	/* link: {0, -PELVISOFFSET, PELVIS2THORAX} */
	Xlink[1][1]=Ai01[1][4];
	Xlink[1][2]=Ai01[2][4];
	Xlink[1][3]=Ai01[3][4];

	Ahmat[1][1][1]=Ai03[1][1];
	Ahmat[1][1][2]=Ai03[1][2];
	Ahmat[1][1][3]=Ai03[1][3];
	Ahmat[1][1][4]=Ai03[1][4];

	Ahmat[1][2][1]=Ai03[2][1];
	Ahmat[1][2][2]=Ai03[2][2];
	Ahmat[1][2][3]=Ai03[2][3];
	Ahmat[1][2][4]=Ai03[2][4];

	Ahmat[1][3][1]=Ai03[3][1];
	Ahmat[1][3][2]=Ai03[3][2];
	Ahmat[1][3][3]=Ai03[3][3];
	Ahmat[1][3][4]=Ai03[3][4];

	Ahmat[1][4][4]=1;


	/* joint ID: 30 */
	Xorigin[30][1]=Ai02[1][4];
	Xorigin[30][2]=Ai02[2][4];
	Xorigin[30][3]=Ai02[3][4];

	Xaxis[30][1]=Ai02[1][3];
	Xaxis[30][2]=Ai02[2][3];
	Xaxis[30][3]=Ai02[3][3];

	/* joint ID: 31 */
	Xorigin[31][1]=Ai03[1][4];
	Xorigin[31][2]=Ai03[2][4];
	Xorigin[31][3]=Ai03[3][4];

	Xaxis[31][1]=Ai03[1][3];
	Xaxis[31][2]=Ai03[2][3];
	Xaxis[31][3]=Ai03[3][3];

	/* joint ID: 1 */
	Xorigin[1][1]=Ai04[1][4];
	Xorigin[1][2]=Ai04[2][4];
	Xorigin[1][3]=Ai04[3][4];

	Xaxis[1][1]=Ai04[1][3];
	Xaxis[1][2]=Ai04[2][3];
	Xaxis[1][3]=Ai04[3][3];

	/* link: {-THORAX2SHOULDER, 0, 0} */
	Xlink[2][1]=Ai04[1][4];
	Xlink[2][2]=Ai04[2][4];
	Xlink[2][3]=Ai04[3][4];

	Ahmat[2][1][1]=Ai04[1][1];
	Ahmat[2][1][2]=Ai04[1][2];
	Ahmat[2][1][3]=Ai04[1][3];
	Ahmat[2][1][4]=Ai04[1][4];

	Ahmat[2][2][1]=Ai04[2][1];
	Ahmat[2][2][2]=Ai04[2][2];
	Ahmat[2][2][3]=Ai04[2][3];
	Ahmat[2][2][4]=Ai04[2][4];

	Ahmat[2][3][1]=Ai04[3][1];
	Ahmat[2][3][2]=Ai04[3][2];
	Ahmat[2][3][3]=Ai04[3][3];
	Ahmat[2][3][4]=Ai04[3][4];

	Ahmat[2][4][4]=1;


	/* joint ID: 2 */
	Xorigin[2][1]=Ai05[1][4];
	Xorigin[2][2]=Ai05[2][4];
	Xorigin[2][3]=Ai05[3][4];

	Xaxis[2][1]=Ai05[1][3];
	Xaxis[2][2]=Ai05[2][3];
	Xaxis[2][3]=Ai05[3][3];

	/* link: {0, 0, -SHOULDERX} */
	Xlink[3][1]=Ai05[1][4];
	Xlink[3][2]=Ai05[2][4];
	Xlink[3][3]=Ai05[3][4];

	Ahmat[3][1][1]=Ai05[1][1];
	Ahmat[3][1][2]=Ai05[1][2];
	Ahmat[3][1][3]=Ai05[1][3];
	Ahmat[3][1][4]=Ai05[1][4];

	Ahmat[3][2][1]=Ai05[2][1];
	Ahmat[3][2][2]=Ai05[2][2];
	Ahmat[3][2][3]=Ai05[2][3];
	Ahmat[3][2][4]=Ai05[2][4];

	Ahmat[3][3][1]=Ai05[3][1];
	Ahmat[3][3][2]=Ai05[3][2];
	Ahmat[3][3][3]=Ai05[3][3];
	Ahmat[3][3][4]=Ai05[3][4];

	Ahmat[3][4][4]=1;


	/* joint ID: 3 */
	Xorigin[3][1]=Ai06[1][4];
	Xorigin[3][2]=Ai06[2][4];
	Xorigin[3][3]=Ai06[3][4];

	Xaxis[3][1]=Ai06[1][3];
	Xaxis[3][2]=Ai06[2][3];
	Xaxis[3][3]=Ai06[3][3];

	/* link: {-SHOULDERY, 0, 0} */
	Xlink[4][1]=Ai06[1][4];
	Xlink[4][2]=Ai06[2][4];
	Xlink[4][3]=Ai06[3][4];

	Ahmat[4][1][1]=Ai06[1][1];
	Ahmat[4][1][2]=Ai06[1][2];
	Ahmat[4][1][3]=Ai06[1][3];
	Ahmat[4][1][4]=Ai06[1][4];

	Ahmat[4][2][1]=Ai06[2][1];
	Ahmat[4][2][2]=Ai06[2][2];
	Ahmat[4][2][3]=Ai06[2][3];
	Ahmat[4][2][4]=Ai06[2][4];

	Ahmat[4][3][1]=Ai06[3][1];
	Ahmat[4][3][2]=Ai06[3][2];
	Ahmat[4][3][3]=Ai06[3][3];
	Ahmat[4][3][4]=Ai06[3][4];

	Ahmat[4][4][4]=1;


	/* joint ID: 4 */
	Xorigin[4][1]=Ai07[1][4];
	Xorigin[4][2]=Ai07[2][4];
	Xorigin[4][3]=Ai07[3][4];

	Xaxis[4][1]=Ai07[1][3];
	Xaxis[4][2]=Ai07[2][3];
	Xaxis[4][3]=Ai07[3][3];

	/* link: {0, 0, -UPPERARM} */
	Xlink[5][1]=Ai07[1][4];
	Xlink[5][2]=Ai07[2][4];
	Xlink[5][3]=Ai07[3][4];

	Ahmat[5][1][1]=Ai08[1][1];
	Ahmat[5][1][2]=Ai08[1][2];
	Ahmat[5][1][3]=Ai08[1][3];
	Ahmat[5][1][4]=Ai08[1][4];

	Ahmat[5][2][1]=Ai08[2][1];
	Ahmat[5][2][2]=Ai08[2][2];
	Ahmat[5][2][3]=Ai08[2][3];
	Ahmat[5][2][4]=Ai08[2][4];

	Ahmat[5][3][1]=Ai08[3][1];
	Ahmat[5][3][2]=Ai08[3][2];
	Ahmat[5][3][3]=Ai08[3][3];
	Ahmat[5][3][4]=Ai08[3][4];

	Ahmat[5][4][4]=1;


	/* joint ID: 5 */
	Xorigin[5][1]=Ai08[1][4];
	Xorigin[5][2]=Ai08[2][4];
	Xorigin[5][3]=Ai08[3][4];

	Xaxis[5][1]=Ai08[1][3];
	Xaxis[5][2]=Ai08[2][3];
	Xaxis[5][3]=Ai08[3][3];

	/* joint ID: 6 */
	Xorigin[6][1]=Ai09[1][4];
	Xorigin[6][2]=Ai09[2][4];
	Xorigin[6][3]=Ai09[3][4];

	Xaxis[6][1]=Ai09[1][3];
	Xaxis[6][2]=Ai09[2][3];
	Xaxis[6][3]=Ai09[3][3];

	/* link: {0, WRISTY, -LOWERARM} */
	Xlink[6][1]=Ai09[1][4];
	Xlink[6][2]=Ai09[2][4];
	Xlink[6][3]=Ai09[3][4];

	Ahmat[6][1][1]=Ai010[1][1];
	Ahmat[6][1][2]=Ai010[1][2];
	Ahmat[6][1][3]=Ai010[1][3];
	Ahmat[6][1][4]=Ai010[1][4];

	Ahmat[6][2][1]=Ai010[2][1];
	Ahmat[6][2][2]=Ai010[2][2];
	Ahmat[6][2][3]=Ai010[2][3];
	Ahmat[6][2][4]=Ai010[2][4];

	Ahmat[6][3][1]=Ai010[3][1];
	Ahmat[6][3][2]=Ai010[3][2];
	Ahmat[6][3][3]=Ai010[3][3];
	Ahmat[6][3][4]=Ai010[3][4];

	Ahmat[6][4][4]=1;


	/* joint ID: 7 */
	Xorigin[7][1]=Ai010[1][4];
	Xorigin[7][2]=Ai010[2][4];
	Xorigin[7][3]=Ai010[3][4];

	Xaxis[7][1]=Ai010[1][3];
	Xaxis[7][2]=Ai010[2][3];
	Xaxis[7][3]=Ai010[3][3];

	/* link: {eff$2$$x[[1]], eff$2$$x[[2]], eff$2$$x[[3]]} */
	Xlink[7][1]=Ai011[1][4];
	Xlink[7][2]=Ai011[2][4];
	Xlink[7][3]=Ai011[3][4];

	Ahmat[7][1][1]=Ai011[1][1];
	Ahmat[7][1][2]=Ai011[1][2];
	Ahmat[7][1][3]=Ai011[1][3];
	Ahmat[7][1][4]=Ai011[1][4];

	Ahmat[7][2][1]=Ai011[2][1];
	Ahmat[7][2][2]=Ai011[2][2];
	Ahmat[7][2][3]=Ai011[2][3];
	Ahmat[7][2][4]=Ai011[2][4];

	Ahmat[7][3][1]=Ai011[3][1];
	Ahmat[7][3][2]=Ai011[3][2];
	Ahmat[7][3][3]=Ai011[3][3];
	Ahmat[7][3][4]=Ai011[3][4];

	Ahmat[7][4][4]=1;


	/* joint ID: 8 */
	Xorigin[8][1]=Ai012[1][4];
	Xorigin[8][2]=Ai012[2][4];
	Xorigin[8][3]=Ai012[3][4];

	Xaxis[8][1]=Ai012[1][3];
	Xaxis[8][2]=Ai012[2][3];
	Xaxis[8][3]=Ai012[3][3];

	/* link: {-THORAX2SHOULDER, 0, 0} */
	Xlink[8][1]=Ai012[1][4];
	Xlink[8][2]=Ai012[2][4];
	Xlink[8][3]=Ai012[3][4];

	Ahmat[8][1][1]=Ai012[1][1];
	Ahmat[8][1][2]=Ai012[1][2];
	Ahmat[8][1][3]=Ai012[1][3];
	Ahmat[8][1][4]=Ai012[1][4];

	Ahmat[8][2][1]=Ai012[2][1];
	Ahmat[8][2][2]=Ai012[2][2];
	Ahmat[8][2][3]=Ai012[2][3];
	Ahmat[8][2][4]=Ai012[2][4];

	Ahmat[8][3][1]=Ai012[3][1];
	Ahmat[8][3][2]=Ai012[3][2];
	Ahmat[8][3][3]=Ai012[3][3];
	Ahmat[8][3][4]=Ai012[3][4];

	Ahmat[8][4][4]=1;


	/* joint ID: 9 */
	Xorigin[9][1]=Ai013[1][4];
	Xorigin[9][2]=Ai013[2][4];
	Xorigin[9][3]=Ai013[3][4];

	Xaxis[9][1]=Ai013[1][3];
	Xaxis[9][2]=Ai013[2][3];
	Xaxis[9][3]=Ai013[3][3];

	/* link: {0, 0, SHOULDERX} */
	Xlink[9][1]=Ai013[1][4];
	Xlink[9][2]=Ai013[2][4];
	Xlink[9][3]=Ai013[3][4];

	Ahmat[9][1][1]=Ai013[1][1];
	Ahmat[9][1][2]=Ai013[1][2];
	Ahmat[9][1][3]=Ai013[1][3];
	Ahmat[9][1][4]=Ai013[1][4];

	Ahmat[9][2][1]=Ai013[2][1];
	Ahmat[9][2][2]=Ai013[2][2];
	Ahmat[9][2][3]=Ai013[2][3];
	Ahmat[9][2][4]=Ai013[2][4];

	Ahmat[9][3][1]=Ai013[3][1];
	Ahmat[9][3][2]=Ai013[3][2];
	Ahmat[9][3][3]=Ai013[3][3];
	Ahmat[9][3][4]=Ai013[3][4];

	Ahmat[9][4][4]=1;


	/* joint ID: 10 */
	Xorigin[10][1]=Ai014[1][4];
	Xorigin[10][2]=Ai014[2][4];
	Xorigin[10][3]=Ai014[3][4];

	Xaxis[10][1]=Ai014[1][3];
	Xaxis[10][2]=Ai014[2][3];
	Xaxis[10][3]=Ai014[3][3];

	/* link: {-SHOULDERY, 0, 0} */
	Xlink[10][1]=Ai014[1][4];
	Xlink[10][2]=Ai014[2][4];
	Xlink[10][3]=Ai014[3][4];

	Ahmat[10][1][1]=Ai014[1][1];
	Ahmat[10][1][2]=Ai014[1][2];
	Ahmat[10][1][3]=Ai014[1][3];
	Ahmat[10][1][4]=Ai014[1][4];

	Ahmat[10][2][1]=Ai014[2][1];
	Ahmat[10][2][2]=Ai014[2][2];
	Ahmat[10][2][3]=Ai014[2][3];
	Ahmat[10][2][4]=Ai014[2][4];

	Ahmat[10][3][1]=Ai014[3][1];
	Ahmat[10][3][2]=Ai014[3][2];
	Ahmat[10][3][3]=Ai014[3][3];
	Ahmat[10][3][4]=Ai014[3][4];

	Ahmat[10][4][4]=1;


	/* joint ID: 11 */
	Xorigin[11][1]=Ai015[1][4];
	Xorigin[11][2]=Ai015[2][4];
	Xorigin[11][3]=Ai015[3][4];

	Xaxis[11][1]=Ai015[1][3];
	Xaxis[11][2]=Ai015[2][3];
	Xaxis[11][3]=Ai015[3][3];

	/* link: {0, 0, UPPERARM} */
	Xlink[11][1]=Ai015[1][4];
	Xlink[11][2]=Ai015[2][4];
	Xlink[11][3]=Ai015[3][4];

	Ahmat[11][1][1]=Ai016[1][1];
	Ahmat[11][1][2]=Ai016[1][2];
	Ahmat[11][1][3]=Ai016[1][3];
	Ahmat[11][1][4]=Ai016[1][4];

	Ahmat[11][2][1]=Ai016[2][1];
	Ahmat[11][2][2]=Ai016[2][2];
	Ahmat[11][2][3]=Ai016[2][3];
	Ahmat[11][2][4]=Ai016[2][4];

	Ahmat[11][3][1]=Ai016[3][1];
	Ahmat[11][3][2]=Ai016[3][2];
	Ahmat[11][3][3]=Ai016[3][3];
	Ahmat[11][3][4]=Ai016[3][4];

	Ahmat[11][4][4]=1;


	/* joint ID: 12 */
	Xorigin[12][1]=Ai016[1][4];
	Xorigin[12][2]=Ai016[2][4];
	Xorigin[12][3]=Ai016[3][4];

	Xaxis[12][1]=Ai016[1][3];
	Xaxis[12][2]=Ai016[2][3];
	Xaxis[12][3]=Ai016[3][3];

	/* joint ID: 13 */
	Xorigin[13][1]=Ai017[1][4];
	Xorigin[13][2]=Ai017[2][4];
	Xorigin[13][3]=Ai017[3][4];

	Xaxis[13][1]=Ai017[1][3];
	Xaxis[13][2]=Ai017[2][3];
	Xaxis[13][3]=Ai017[3][3];

	/* link: {0, WRISTY, LOWERARM} */
	Xlink[12][1]=Ai017[1][4];
	Xlink[12][2]=Ai017[2][4];
	Xlink[12][3]=Ai017[3][4];

	Ahmat[12][1][1]=Ai018[1][1];
	Ahmat[12][1][2]=Ai018[1][2];
	Ahmat[12][1][3]=Ai018[1][3];
	Ahmat[12][1][4]=Ai018[1][4];

	Ahmat[12][2][1]=Ai018[2][1];
	Ahmat[12][2][2]=Ai018[2][2];
	Ahmat[12][2][3]=Ai018[2][3];
	Ahmat[12][2][4]=Ai018[2][4];

	Ahmat[12][3][1]=Ai018[3][1];
	Ahmat[12][3][2]=Ai018[3][2];
	Ahmat[12][3][3]=Ai018[3][3];
	Ahmat[12][3][4]=Ai018[3][4];

	Ahmat[12][4][4]=1;


	/* joint ID: 14 */
	Xorigin[14][1]=Ai018[1][4];
	Xorigin[14][2]=Ai018[2][4];
	Xorigin[14][3]=Ai018[3][4];

	Xaxis[14][1]=Ai018[1][3];
	Xaxis[14][2]=Ai018[2][3];
	Xaxis[14][3]=Ai018[3][3];

	/* link: {eff$1$$x[[1]], eff$1$$x[[2]], eff$1$$x[[3]]} */
	Xlink[13][1]=Ai019[1][4];
	Xlink[13][2]=Ai019[2][4];
	Xlink[13][3]=Ai019[3][4];

	Ahmat[13][1][1]=Ai019[1][1];
	Ahmat[13][1][2]=Ai019[1][2];
	Ahmat[13][1][3]=Ai019[1][3];
	Ahmat[13][1][4]=Ai019[1][4];

	Ahmat[13][2][1]=Ai019[2][1];
	Ahmat[13][2][2]=Ai019[2][2];
	Ahmat[13][2][3]=Ai019[2][3];
	Ahmat[13][2][4]=Ai019[2][4];

	Ahmat[13][3][1]=Ai019[3][1];
	Ahmat[13][3][2]=Ai019[3][2];
	Ahmat[13][3][3]=Ai019[3][3];
	Ahmat[13][3][4]=Ai019[3][4];

	Ahmat[13][4][4]=1;


	/* joint ID: 32 */
	Xorigin[32][1]=Ai020[1][4];
	Xorigin[32][2]=Ai020[2][4];
	Xorigin[32][3]=Ai020[3][4];

	Xaxis[32][1]=Ai020[1][3];
	Xaxis[32][2]=Ai020[2][3];
	Xaxis[32][3]=Ai020[3][3];

	/* link: {-THORAX2NECK, 0, 0} */
	Xlink[14][1]=Ai020[1][4];
	Xlink[14][2]=Ai020[2][4];
	Xlink[14][3]=Ai020[3][4];

	Ahmat[14][1][1]=Ai020[1][1];
	Ahmat[14][1][2]=Ai020[1][2];
	Ahmat[14][1][3]=Ai020[1][3];
	Ahmat[14][1][4]=Ai020[1][4];

	Ahmat[14][2][1]=Ai020[2][1];
	Ahmat[14][2][2]=Ai020[2][2];
	Ahmat[14][2][3]=Ai020[2][3];
	Ahmat[14][2][4]=Ai020[2][4];

	Ahmat[14][3][1]=Ai020[3][1];
	Ahmat[14][3][2]=Ai020[3][2];
	Ahmat[14][3][3]=Ai020[3][3];
	Ahmat[14][3][4]=Ai020[3][4];

	Ahmat[14][4][4]=1;


	/* joint ID: 33 */
	Xorigin[33][1]=Ai021[1][4];
	Xorigin[33][2]=Ai021[2][4];
	Xorigin[33][3]=Ai021[3][4];

	Xaxis[33][1]=Ai021[1][3];
	Xaxis[33][2]=Ai021[2][3];
	Xaxis[33][3]=Ai021[3][3];

	/* link: {0, -CERVICAL, 0} */
	Xlink[15][1]=Ai021[1][4];
	Xlink[15][2]=Ai021[2][4];
	Xlink[15][3]=Ai021[3][4];

	Ahmat[15][1][1]=Ai022[1][1];
	Ahmat[15][1][2]=Ai022[1][2];
	Ahmat[15][1][3]=Ai022[1][3];
	Ahmat[15][1][4]=Ai022[1][4];

	Ahmat[15][2][1]=Ai022[2][1];
	Ahmat[15][2][2]=Ai022[2][2];
	Ahmat[15][2][3]=Ai022[2][3];
	Ahmat[15][2][4]=Ai022[2][4];

	Ahmat[15][3][1]=Ai022[3][1];
	Ahmat[15][3][2]=Ai022[3][2];
	Ahmat[15][3][3]=Ai022[3][3];
	Ahmat[15][3][4]=Ai022[3][4];

	Ahmat[15][4][4]=1;


	/* joint ID: 34 */
	Xorigin[34][1]=Ai022[1][4];
	Xorigin[34][2]=Ai022[2][4];
	Xorigin[34][3]=Ai022[3][4];

	Xaxis[34][1]=Ai022[1][3];
	Xaxis[34][2]=Ai022[2][3];
	Xaxis[34][3]=Ai022[3][3];

	/* joint ID: 35 */
	Xorigin[35][1]=Ai023[1][4];
	Xorigin[35][2]=Ai023[2][4];
	Xorigin[35][3]=Ai023[3][4];

	Xaxis[35][1]=Ai023[1][3];
	Xaxis[35][2]=Ai023[2][3];
	Xaxis[35][3]=Ai023[3][3];

	/* link: {EYEXOFF, -EYEYOFF, -HEAD} */
	Xlink[16][1]=Ai023[1][4];
	Xlink[16][2]=Ai023[2][4];
	Xlink[16][3]=Ai023[3][4];

	Ahmat[16][1][1]=Ai024[1][1];
	Ahmat[16][1][2]=Ai024[1][2];
	Ahmat[16][1][3]=Ai024[1][3];
	Ahmat[16][1][4]=Ai024[1][4];

	Ahmat[16][2][1]=Ai024[2][1];
	Ahmat[16][2][2]=Ai024[2][2];
	Ahmat[16][2][3]=Ai024[2][3];
	Ahmat[16][2][4]=Ai024[2][4];

	Ahmat[16][3][1]=Ai024[3][1];
	Ahmat[16][3][2]=Ai024[3][2];
	Ahmat[16][3][3]=Ai024[3][3];
	Ahmat[16][3][4]=Ai024[3][4];

	Ahmat[16][4][4]=1;


	/* joint ID: 36 */
	Xorigin[36][1]=Ai024[1][4];
	Xorigin[36][2]=Ai024[2][4];
	Xorigin[36][3]=Ai024[3][4];

	Xaxis[36][1]=Ai024[1][3];
	Xaxis[36][2]=Ai024[2][3];
	Xaxis[36][3]=Ai024[3][3];

	/* link: {0, -EYE, 0} */
	Xlink[17][1]=Ai025[1][4];
	Xlink[17][2]=Ai025[2][4];
	Xlink[17][3]=Ai025[3][4];

	Ahmat[17][1][1]=Ai025[1][1];
	Ahmat[17][1][2]=Ai025[1][2];
	Ahmat[17][1][3]=Ai025[1][3];
	Ahmat[17][1][4]=Ai025[1][4];

	Ahmat[17][2][1]=Ai025[2][1];
	Ahmat[17][2][2]=Ai025[2][2];
	Ahmat[17][2][3]=Ai025[2][3];
	Ahmat[17][2][4]=Ai025[2][4];

	Ahmat[17][3][1]=Ai025[3][1];
	Ahmat[17][3][2]=Ai025[3][2];
	Ahmat[17][3][3]=Ai025[3][3];
	Ahmat[17][3][4]=Ai025[3][4];

	Ahmat[17][4][4]=1;


	/* joint ID: 37 */
	Xorigin[37][1]=Ai026[1][4];
	Xorigin[37][2]=Ai026[2][4];
	Xorigin[37][3]=Ai026[3][4];

	Xaxis[37][1]=Ai026[1][3];
	Xaxis[37][2]=Ai026[2][3];
	Xaxis[37][3]=Ai026[3][3];

	/* link: {-EYEXOFF, -EYEYOFF, -HEAD} */
	Xlink[18][1]=Ai026[1][4];
	Xlink[18][2]=Ai026[2][4];
	Xlink[18][3]=Ai026[3][4];

	Ahmat[18][1][1]=Ai027[1][1];
	Ahmat[18][1][2]=Ai027[1][2];
	Ahmat[18][1][3]=Ai027[1][3];
	Ahmat[18][1][4]=Ai027[1][4];

	Ahmat[18][2][1]=Ai027[2][1];
	Ahmat[18][2][2]=Ai027[2][2];
	Ahmat[18][2][3]=Ai027[2][3];
	Ahmat[18][2][4]=Ai027[2][4];

	Ahmat[18][3][1]=Ai027[3][1];
	Ahmat[18][3][2]=Ai027[3][2];
	Ahmat[18][3][3]=Ai027[3][3];
	Ahmat[18][3][4]=Ai027[3][4];

	Ahmat[18][4][4]=1;


	/* joint ID: 38 */
	Xorigin[38][1]=Ai027[1][4];
	Xorigin[38][2]=Ai027[2][4];
	Xorigin[38][3]=Ai027[3][4];

	Xaxis[38][1]=Ai027[1][3];
	Xaxis[38][2]=Ai027[2][3];
	Xaxis[38][3]=Ai027[3][3];

	/* link: {0, -EYE, 0} */
	Xlink[19][1]=Ai028[1][4];
	Xlink[19][2]=Ai028[2][4];
	Xlink[19][3]=Ai028[3][4];

	Ahmat[19][1][1]=Ai028[1][1];
	Ahmat[19][1][2]=Ai028[1][2];
	Ahmat[19][1][3]=Ai028[1][3];
	Ahmat[19][1][4]=Ai028[1][4];

	Ahmat[19][2][1]=Ai028[2][1];
	Ahmat[19][2][2]=Ai028[2][2];
	Ahmat[19][2][3]=Ai028[2][3];
	Ahmat[19][2][4]=Ai028[2][4];

	Ahmat[19][3][1]=Ai028[3][1];
	Ahmat[19][3][2]=Ai028[3][2];
	Ahmat[19][3][3]=Ai028[3][3];
	Ahmat[19][3][4]=Ai028[3][4];

	Ahmat[19][4][4]=1;


	/* link: {0, 0, -TOPofHEAD} */
	Xlink[20][1]=Ai029[1][4];
	Xlink[20][2]=Ai029[2][4];
	Xlink[20][3]=Ai029[3][4];

	Ahmat[20][1][1]=Ai029[1][1];
	Ahmat[20][1][2]=Ai029[1][2];
	Ahmat[20][1][3]=Ai029[1][3];
	Ahmat[20][1][4]=Ai029[1][4];

	Ahmat[20][2][1]=Ai029[2][1];
	Ahmat[20][2][2]=Ai029[2][2];
	Ahmat[20][2][3]=Ai029[2][3];
	Ahmat[20][2][4]=Ai029[2][4];

	Ahmat[20][3][1]=Ai029[3][1];
	Ahmat[20][3][2]=Ai029[3][2];
	Ahmat[20][3][3]=Ai029[3][3];
	Ahmat[20][3][4]=Ai029[3][4];

	Ahmat[20][4][4]=1;


	/* joint ID: 23 */
	Xorigin[23][1]=Ai030[1][4];
	Xorigin[23][2]=Ai030[2][4];
	Xorigin[23][3]=Ai030[3][4];

	Xaxis[23][1]=Ai030[1][3];
	Xaxis[23][2]=Ai030[2][3];
	Xaxis[23][3]=Ai030[3][3];

	/* link: {XHIP, 0, 0} */
	Xlink[21][1]=Ai030[1][4];
	Xlink[21][2]=Ai030[2][4];
	Xlink[21][3]=Ai030[3][4];

	Ahmat[21][1][1]=Ai031[1][1];
	Ahmat[21][1][2]=Ai031[1][2];
	Ahmat[21][1][3]=Ai031[1][3];
	Ahmat[21][1][4]=Ai031[1][4];

	Ahmat[21][2][1]=Ai031[2][1];
	Ahmat[21][2][2]=Ai031[2][2];
	Ahmat[21][2][3]=Ai031[2][3];
	Ahmat[21][2][4]=Ai031[2][4];

	Ahmat[21][3][1]=Ai031[3][1];
	Ahmat[21][3][2]=Ai031[3][2];
	Ahmat[21][3][3]=Ai031[3][3];
	Ahmat[21][3][4]=Ai031[3][4];

	Ahmat[21][4][4]=1;


	/* joint ID: 22 */
	Xorigin[22][1]=Ai031[1][4];
	Xorigin[22][2]=Ai031[2][4];
	Xorigin[22][3]=Ai031[3][4];

	Xaxis[22][1]=Ai031[1][3];
	Xaxis[22][2]=Ai031[2][3];
	Xaxis[22][3]=Ai031[3][3];

	/* joint ID: 24 */
	Xorigin[24][1]=Ai032[1][4];
	Xorigin[24][2]=Ai032[2][4];
	Xorigin[24][3]=Ai032[3][4];

	Xaxis[24][1]=Ai032[1][3];
	Xaxis[24][2]=Ai032[2][3];
	Xaxis[24][3]=Ai032[3][3];

	/* link: {YHIP, 0, 0} */
	Xlink[22][1]=Ai032[1][4];
	Xlink[22][2]=Ai032[2][4];
	Xlink[22][3]=Ai032[3][4];

	Ahmat[22][1][1]=Ai032[1][1];
	Ahmat[22][1][2]=Ai032[1][2];
	Ahmat[22][1][3]=Ai032[1][3];
	Ahmat[22][1][4]=Ai032[1][4];

	Ahmat[22][2][1]=Ai032[2][1];
	Ahmat[22][2][2]=Ai032[2][2];
	Ahmat[22][2][3]=Ai032[2][3];
	Ahmat[22][2][4]=Ai032[2][4];

	Ahmat[22][3][1]=Ai032[3][1];
	Ahmat[22][3][2]=Ai032[3][2];
	Ahmat[22][3][3]=Ai032[3][3];
	Ahmat[22][3][4]=Ai032[3][4];

	Ahmat[22][4][4]=1;


	/* joint ID: 25 */
	Xorigin[25][1]=Ai033[1][4];
	Xorigin[25][2]=Ai033[2][4];
	Xorigin[25][3]=Ai033[3][4];

	Xaxis[25][1]=Ai033[1][3];
	Xaxis[25][2]=Ai033[2][3];
	Xaxis[25][3]=Ai033[3][3];

	/* link: {YKNEE, 0, UPPERLEG} */
	Xlink[23][1]=Ai033[1][4];
	Xlink[23][2]=Ai033[2][4];
	Xlink[23][3]=Ai033[3][4];

	Ahmat[23][1][1]=Ai034[1][1];
	Ahmat[23][1][2]=Ai034[1][2];
	Ahmat[23][1][3]=Ai034[1][3];
	Ahmat[23][1][4]=Ai034[1][4];

	Ahmat[23][2][1]=Ai034[2][1];
	Ahmat[23][2][2]=Ai034[2][2];
	Ahmat[23][2][3]=Ai034[2][3];
	Ahmat[23][2][4]=Ai034[2][4];

	Ahmat[23][3][1]=Ai034[3][1];
	Ahmat[23][3][2]=Ai034[3][2];
	Ahmat[23][3][3]=Ai034[3][3];
	Ahmat[23][3][4]=Ai034[3][4];

	Ahmat[23][4][4]=1;


	/* joint ID: 26 */
	Xorigin[26][1]=Ai034[1][4];
	Xorigin[26][2]=Ai034[2][4];
	Xorigin[26][3]=Ai034[3][4];

	Xaxis[26][1]=Ai034[1][3];
	Xaxis[26][2]=Ai034[2][3];
	Xaxis[26][3]=Ai034[3][3];

	/* joint ID: 27 */
	Xorigin[27][1]=Ai035[1][4];
	Xorigin[27][2]=Ai035[2][4];
	Xorigin[27][3]=Ai035[3][4];

	Xaxis[27][1]=Ai035[1][3];
	Xaxis[27][2]=Ai035[2][3];
	Xaxis[27][3]=Ai035[3][3];

	/* link: {0, 0, LOWERLEG} */
	Xlink[24][1]=Ai035[1][4];
	Xlink[24][2]=Ai035[2][4];
	Xlink[24][3]=Ai035[3][4];

	Ahmat[24][1][1]=Ai036[1][1];
	Ahmat[24][1][2]=Ai036[1][2];
	Ahmat[24][1][3]=Ai036[1][3];
	Ahmat[24][1][4]=Ai036[1][4];

	Ahmat[24][2][1]=Ai036[2][1];
	Ahmat[24][2][2]=Ai036[2][2];
	Ahmat[24][2][3]=Ai036[2][3];
	Ahmat[24][2][4]=Ai036[2][4];

	Ahmat[24][3][1]=Ai036[3][1];
	Ahmat[24][3][2]=Ai036[3][2];
	Ahmat[24][3][3]=Ai036[3][3];
	Ahmat[24][3][4]=Ai036[3][4];

	Ahmat[24][4][4]=1;


	/* joint ID: 28 */
	Xorigin[28][1]=Ai036[1][4];
	Xorigin[28][2]=Ai036[2][4];
	Xorigin[28][3]=Ai036[3][4];

	Xaxis[28][1]=Ai036[1][3];
	Xaxis[28][2]=Ai036[2][3];
	Xaxis[28][3]=Ai036[3][3];

	/* link: {ZTOE, -XTOE, YTOE} */
	Xlink[25][1]=Ai037[1][4];
	Xlink[25][2]=Ai037[2][4];
	Xlink[25][3]=Ai037[3][4];

	Ahmat[25][1][1]=Ai037[1][1];
	Ahmat[25][1][2]=Ai037[1][2];
	Ahmat[25][1][3]=Ai037[1][3];
	Ahmat[25][1][4]=Ai037[1][4];

	Ahmat[25][2][1]=Ai037[2][1];
	Ahmat[25][2][2]=Ai037[2][2];
	Ahmat[25][2][3]=Ai037[2][3];
	Ahmat[25][2][4]=Ai037[2][4];

	Ahmat[25][3][1]=Ai037[3][1];
	Ahmat[25][3][2]=Ai037[3][2];
	Ahmat[25][3][3]=Ai037[3][3];
	Ahmat[25][3][4]=Ai037[3][4];

	Ahmat[25][4][4]=1;


	/* link: {ZTOE, XTOE, YTOE} */
	Xlink[26][1]=Ai038[1][4];
	Xlink[26][2]=Ai038[2][4];
	Xlink[26][3]=Ai038[3][4];

	Ahmat[26][1][1]=Ai038[1][1];
	Ahmat[26][1][2]=Ai038[1][2];
	Ahmat[26][1][3]=Ai038[1][3];
	Ahmat[26][1][4]=Ai038[1][4];

	Ahmat[26][2][1]=Ai038[2][1];
	Ahmat[26][2][2]=Ai038[2][2];
	Ahmat[26][2][3]=Ai038[2][3];
	Ahmat[26][2][4]=Ai038[2][4];

	Ahmat[26][3][1]=Ai038[3][1];
	Ahmat[26][3][2]=Ai038[3][2];
	Ahmat[26][3][3]=Ai038[3][3];
	Ahmat[26][3][4]=Ai038[3][4];

	Ahmat[26][4][4]=1;


	/* link: {ZHEEL, -XHEEL, -YHEEL} */
	Xlink[27][1]=Ai039[1][4];
	Xlink[27][2]=Ai039[2][4];
	Xlink[27][3]=Ai039[3][4];

	Ahmat[27][1][1]=Ai039[1][1];
	Ahmat[27][1][2]=Ai039[1][2];
	Ahmat[27][1][3]=Ai039[1][3];
	Ahmat[27][1][4]=Ai039[1][4];

	Ahmat[27][2][1]=Ai039[2][1];
	Ahmat[27][2][2]=Ai039[2][2];
	Ahmat[27][2][3]=Ai039[2][3];
	Ahmat[27][2][4]=Ai039[2][4];

	Ahmat[27][3][1]=Ai039[3][1];
	Ahmat[27][3][2]=Ai039[3][2];
	Ahmat[27][3][3]=Ai039[3][3];
	Ahmat[27][3][4]=Ai039[3][4];

	Ahmat[27][4][4]=1;


	/* link: {ZHEEL, XHEEL, -YHEEL} */
	Xlink[28][1]=Ai040[1][4];
	Xlink[28][2]=Ai040[2][4];
	Xlink[28][3]=Ai040[3][4];

	Ahmat[28][1][1]=Ai040[1][1];
	Ahmat[28][1][2]=Ai040[1][2];
	Ahmat[28][1][3]=Ai040[1][3];
	Ahmat[28][1][4]=Ai040[1][4];

	Ahmat[28][2][1]=Ai040[2][1];
	Ahmat[28][2][2]=Ai040[2][2];
	Ahmat[28][2][3]=Ai040[2][3];
	Ahmat[28][2][4]=Ai040[2][4];

	Ahmat[28][3][1]=Ai040[3][1];
	Ahmat[28][3][2]=Ai040[3][2];
	Ahmat[28][3][3]=Ai040[3][3];
	Ahmat[28][3][4]=Ai040[3][4];

	Ahmat[28][4][4]=1;


	/* link: {eff$3$$x[[1]], eff$3$$x[[2]], eff$3$$x[[3]]} */
	Xlink[29][1]=Ai041[1][4];
	Xlink[29][2]=Ai041[2][4];
	Xlink[29][3]=Ai041[3][4];

	Ahmat[29][1][1]=Ai041[1][1];
	Ahmat[29][1][2]=Ai041[1][2];
	Ahmat[29][1][3]=Ai041[1][3];
	Ahmat[29][1][4]=Ai041[1][4];

	Ahmat[29][2][1]=Ai041[2][1];
	Ahmat[29][2][2]=Ai041[2][2];
	Ahmat[29][2][3]=Ai041[2][3];
	Ahmat[29][2][4]=Ai041[2][4];

	Ahmat[29][3][1]=Ai041[3][1];
	Ahmat[29][3][2]=Ai041[3][2];
	Ahmat[29][3][3]=Ai041[3][3];
	Ahmat[29][3][4]=Ai041[3][4];

	Ahmat[29][4][4]=1;


	/* joint ID: 16 */
	Xorigin[16][1]=Ai042[1][4];
	Xorigin[16][2]=Ai042[2][4];
	Xorigin[16][3]=Ai042[3][4];

	Xaxis[16][1]=Ai042[1][3];
	Xaxis[16][2]=Ai042[2][3];
	Xaxis[16][3]=Ai042[3][3];

	/* link: {-XHIP, 0, 0} */
	Xlink[30][1]=Ai042[1][4];
	Xlink[30][2]=Ai042[2][4];
	Xlink[30][3]=Ai042[3][4];

	Ahmat[30][1][1]=Ai043[1][1];
	Ahmat[30][1][2]=Ai043[1][2];
	Ahmat[30][1][3]=Ai043[1][3];
	Ahmat[30][1][4]=Ai043[1][4];

	Ahmat[30][2][1]=Ai043[2][1];
	Ahmat[30][2][2]=Ai043[2][2];
	Ahmat[30][2][3]=Ai043[2][3];
	Ahmat[30][2][4]=Ai043[2][4];

	Ahmat[30][3][1]=Ai043[3][1];
	Ahmat[30][3][2]=Ai043[3][2];
	Ahmat[30][3][3]=Ai043[3][3];
	Ahmat[30][3][4]=Ai043[3][4];

	Ahmat[30][4][4]=1;


	/* joint ID: 15 */
	Xorigin[15][1]=Ai043[1][4];
	Xorigin[15][2]=Ai043[2][4];
	Xorigin[15][3]=Ai043[3][4];

	Xaxis[15][1]=Ai043[1][3];
	Xaxis[15][2]=Ai043[2][3];
	Xaxis[15][3]=Ai043[3][3];

	/* joint ID: 17 */
	Xorigin[17][1]=Ai044[1][4];
	Xorigin[17][2]=Ai044[2][4];
	Xorigin[17][3]=Ai044[3][4];

	Xaxis[17][1]=Ai044[1][3];
	Xaxis[17][2]=Ai044[2][3];
	Xaxis[17][3]=Ai044[3][3];

	/* link: {YHIP, 0, 0} */
	Xlink[31][1]=Ai044[1][4];
	Xlink[31][2]=Ai044[2][4];
	Xlink[31][3]=Ai044[3][4];

	Ahmat[31][1][1]=Ai044[1][1];
	Ahmat[31][1][2]=Ai044[1][2];
	Ahmat[31][1][3]=Ai044[1][3];
	Ahmat[31][1][4]=Ai044[1][4];

	Ahmat[31][2][1]=Ai044[2][1];
	Ahmat[31][2][2]=Ai044[2][2];
	Ahmat[31][2][3]=Ai044[2][3];
	Ahmat[31][2][4]=Ai044[2][4];

	Ahmat[31][3][1]=Ai044[3][1];
	Ahmat[31][3][2]=Ai044[3][2];
	Ahmat[31][3][3]=Ai044[3][3];
	Ahmat[31][3][4]=Ai044[3][4];

	Ahmat[31][4][4]=1;


	/* joint ID: 18 */
	Xorigin[18][1]=Ai045[1][4];
	Xorigin[18][2]=Ai045[2][4];
	Xorigin[18][3]=Ai045[3][4];

	Xaxis[18][1]=Ai045[1][3];
	Xaxis[18][2]=Ai045[2][3];
	Xaxis[18][3]=Ai045[3][3];

	/* link: {YKNEE, 0, -UPPERLEG} */
	Xlink[32][1]=Ai045[1][4];
	Xlink[32][2]=Ai045[2][4];
	Xlink[32][3]=Ai045[3][4];

	Ahmat[32][1][1]=Ai046[1][1];
	Ahmat[32][1][2]=Ai046[1][2];
	Ahmat[32][1][3]=Ai046[1][3];
	Ahmat[32][1][4]=Ai046[1][4];

	Ahmat[32][2][1]=Ai046[2][1];
	Ahmat[32][2][2]=Ai046[2][2];
	Ahmat[32][2][3]=Ai046[2][3];
	Ahmat[32][2][4]=Ai046[2][4];

	Ahmat[32][3][1]=Ai046[3][1];
	Ahmat[32][3][2]=Ai046[3][2];
	Ahmat[32][3][3]=Ai046[3][3];
	Ahmat[32][3][4]=Ai046[3][4];

	Ahmat[32][4][4]=1;


	/* joint ID: 19 */
	Xorigin[19][1]=Ai046[1][4];
	Xorigin[19][2]=Ai046[2][4];
	Xorigin[19][3]=Ai046[3][4];

	Xaxis[19][1]=Ai046[1][3];
	Xaxis[19][2]=Ai046[2][3];
	Xaxis[19][3]=Ai046[3][3];

	/* joint ID: 20 */
	Xorigin[20][1]=Ai047[1][4];
	Xorigin[20][2]=Ai047[2][4];
	Xorigin[20][3]=Ai047[3][4];

	Xaxis[20][1]=Ai047[1][3];
	Xaxis[20][2]=Ai047[2][3];
	Xaxis[20][3]=Ai047[3][3];

	/* link: {0, 0, -LOWERLEG} */
	Xlink[33][1]=Ai047[1][4];
	Xlink[33][2]=Ai047[2][4];
	Xlink[33][3]=Ai047[3][4];

	Ahmat[33][1][1]=Ai048[1][1];
	Ahmat[33][1][2]=Ai048[1][2];
	Ahmat[33][1][3]=Ai048[1][3];
	Ahmat[33][1][4]=Ai048[1][4];

	Ahmat[33][2][1]=Ai048[2][1];
	Ahmat[33][2][2]=Ai048[2][2];
	Ahmat[33][2][3]=Ai048[2][3];
	Ahmat[33][2][4]=Ai048[2][4];

	Ahmat[33][3][1]=Ai048[3][1];
	Ahmat[33][3][2]=Ai048[3][2];
	Ahmat[33][3][3]=Ai048[3][3];
	Ahmat[33][3][4]=Ai048[3][4];

	Ahmat[33][4][4]=1;


	/* joint ID: 21 */
	Xorigin[21][1]=Ai048[1][4];
	Xorigin[21][2]=Ai048[2][4];
	Xorigin[21][3]=Ai048[3][4];

	Xaxis[21][1]=Ai048[1][3];
	Xaxis[21][2]=Ai048[2][3];
	Xaxis[21][3]=Ai048[3][3];

	/* link: {ZTOE, XTOE, -YTOE} */
	Xlink[34][1]=Ai049[1][4];
	Xlink[34][2]=Ai049[2][4];
	Xlink[34][3]=Ai049[3][4];

	Ahmat[34][1][1]=Ai049[1][1];
	Ahmat[34][1][2]=Ai049[1][2];
	Ahmat[34][1][3]=Ai049[1][3];
	Ahmat[34][1][4]=Ai049[1][4];

	Ahmat[34][2][1]=Ai049[2][1];
	Ahmat[34][2][2]=Ai049[2][2];
	Ahmat[34][2][3]=Ai049[2][3];
	Ahmat[34][2][4]=Ai049[2][4];

	Ahmat[34][3][1]=Ai049[3][1];
	Ahmat[34][3][2]=Ai049[3][2];
	Ahmat[34][3][3]=Ai049[3][3];
	Ahmat[34][3][4]=Ai049[3][4];

	Ahmat[34][4][4]=1;


	/* link: {ZTOE, -XTOE, -YTOE} */
	Xlink[35][1]=Ai050[1][4];
	Xlink[35][2]=Ai050[2][4];
	Xlink[35][3]=Ai050[3][4];

	Ahmat[35][1][1]=Ai050[1][1];
	Ahmat[35][1][2]=Ai050[1][2];
	Ahmat[35][1][3]=Ai050[1][3];
	Ahmat[35][1][4]=Ai050[1][4];

	Ahmat[35][2][1]=Ai050[2][1];
	Ahmat[35][2][2]=Ai050[2][2];
	Ahmat[35][2][3]=Ai050[2][3];
	Ahmat[35][2][4]=Ai050[2][4];

	Ahmat[35][3][1]=Ai050[3][1];
	Ahmat[35][3][2]=Ai050[3][2];
	Ahmat[35][3][3]=Ai050[3][3];
	Ahmat[35][3][4]=Ai050[3][4];

	Ahmat[35][4][4]=1;


	/* link: {ZHEEL, XHEEL, YHEEL} */
	Xlink[36][1]=Ai051[1][4];
	Xlink[36][2]=Ai051[2][4];
	Xlink[36][3]=Ai051[3][4];

	Ahmat[36][1][1]=Ai051[1][1];
	Ahmat[36][1][2]=Ai051[1][2];
	Ahmat[36][1][3]=Ai051[1][3];
	Ahmat[36][1][4]=Ai051[1][4];

	Ahmat[36][2][1]=Ai051[2][1];
	Ahmat[36][2][2]=Ai051[2][2];
	Ahmat[36][2][3]=Ai051[2][3];
	Ahmat[36][2][4]=Ai051[2][4];

	Ahmat[36][3][1]=Ai051[3][1];
	Ahmat[36][3][2]=Ai051[3][2];
	Ahmat[36][3][3]=Ai051[3][3];
	Ahmat[36][3][4]=Ai051[3][4];

	Ahmat[36][4][4]=1;


	/* link: {ZHEEL, -XHEEL, YHEEL} */
	Xlink[37][1]=Ai052[1][4];
	Xlink[37][2]=Ai052[2][4];
	Xlink[37][3]=Ai052[3][4];

	Ahmat[37][1][1]=Ai052[1][1];
	Ahmat[37][1][2]=Ai052[1][2];
	Ahmat[37][1][3]=Ai052[1][3];
	Ahmat[37][1][4]=Ai052[1][4];

	Ahmat[37][2][1]=Ai052[2][1];
	Ahmat[37][2][2]=Ai052[2][2];
	Ahmat[37][2][3]=Ai052[2][3];
	Ahmat[37][2][4]=Ai052[2][4];

	Ahmat[37][3][1]=Ai052[3][1];
	Ahmat[37][3][2]=Ai052[3][2];
	Ahmat[37][3][3]=Ai052[3][3];
	Ahmat[37][3][4]=Ai052[3][4];

	Ahmat[37][4][4]=1;


	/* link: {eff$4$$x[[1]], eff$4$$x[[2]], eff$4$$x[[3]]} */
	Xlink[38][1]=Ai053[1][4];
	Xlink[38][2]=Ai053[2][4];
	Xlink[38][3]=Ai053[3][4];

	Ahmat[38][1][1]=Ai053[1][1];
	Ahmat[38][1][2]=Ai053[1][2];
	Ahmat[38][1][3]=Ai053[1][3];
	Ahmat[38][1][4]=Ai053[1][4];

	Ahmat[38][2][1]=Ai053[2][1];
	Ahmat[38][2][2]=Ai053[2][2];
	Ahmat[38][2][3]=Ai053[2][3];
	Ahmat[38][2][4]=Ai053[2][4];

	Ahmat[38][3][1]=Ai053[3][1];
	Ahmat[38][3][2]=Ai053[3][2];
	Ahmat[38][3][3]=Ai053[3][3];
	Ahmat[38][3][4]=Ai053[3][4];

	Ahmat[38][4][4]=1;


}

static void jacobian(const double link[NLINK+1][4],
					 const double origin[NDOF+1][4],
					 const double axis[NDOF+1][4],
					 mat & jac) {

	static int firsttime = true;
	static mat Jlist = zeros<mat>(NENDEFF,NDOF);

	if (firsttime) {
		firsttime = false;
		Jlist(1,span(7,13)) = ones<rowvec>(7); // right arm
		Jlist(1,span(28,30)) = ones<rowvec>(3);
		Jlist(0,span(0,6)) = ones<rowvec>(7); // left arm
		Jlist(0,span(28,30)) = ones<rowvec>(3);
		Jlist(2,span(21,27)) = ones<rowvec>(7);
		Jlist(3,span(14,20)) = ones<rowvec>(7);
	}

	static vec6 c;
	for (int i = 0; i < NENDEFF; i++) {
		for (int j = 0; j < NDOF; ++j) {
			if (Jlist(i,j) != 0) {
				c(X) = axis[j+1][2] * (link[link2endeffmap(i)][3] - origin[j+1][3]) - axis[j+1][3] * (link[link2endeffmap(i)][2]-origin[j+1][2]);
				c(Y) = axis[j+1][3] * (link[link2endeffmap(i)][1] - origin[j+1][1]) - axis[j+1][1] * (link[link2endeffmap(i)][3]-origin[j+1][3]);
				c(Z) = axis[j+1][1] * (link[link2endeffmap(i)][2] - origin[j+1][2]) - axis[j+1][2] * (link[link2endeffmap(i)][1]-origin[j+1][1]);
				c(DX) = axis[j+1][1];
				c(DY) = axis[j+1][2];
				c(DZ) = axis[j+1][3];
				jac(span(i*2*NCART,(i+1)*2*NCART-1),j) = c;
			}
			else {
				jac(span(i*2*NCART,(i+1)*2*NCART-1),j) = zeros<vec>(2*NCART);
			}
		}
	}
}

/**
 *
 * @brief Reads joint limits from file.
 *
 * @param lb Array of lower bound values to be loaded.
 * @param ub Array of upper bound values to be loaded.
 * @return If can load the joint limits successfully then returns 1.
 *
 */
void read_joint_limits(vec & lb, vec & ub) {

	using namespace std;
	mat limits;
	string homename = getenv("HOME");
	string fullname = homename + "/basketball/config/joint_limits.cfg";
	limits.load(fullname);
	lb = limits.col(0);
	ub = limits.col(1);
}

/**
 * @brief Reads the default joint states (initial posture) from file.
 *
 *
 */
static void read_default_state(vec & q_default) {

	using namespace std;
	mat limits;
	string homename = getenv("HOME");
	string fullname = homename + "/basketball/config/joint_limits.cfg";
	limits.load(fullname);
	q_default = limits.col(2);
}
