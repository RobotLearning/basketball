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
		double jac[2*NCART][NDOF]);
static void read_default_state(vec & q_default);


/**
 * @brief Returns the cartesian endeffector positions
 */
void get_position(const ivec & active_dofs, const double q_active[NDOF_ACTIVE],
		          double pos[NCART]) {

	const int RIGHT_PALM = R_WAA;
	static double link[NLINK+1][3+1];
	static double origin[NDOF+1][3+1];
	static double axis[NDOF+1][3+1];
	static double amats[NDOF+1][4+1][4+1];
	static vec q = zeros<vec>(NDOF);
	static bool firsttime = true;

	if (firsttime) {
		read_default_state(q);
		firsttime = false;
	}
	for (int i = 0; i < NDOF_ACTIVE; i++) {
		q(active_dofs[i]) = q_active[i];
	}

	kinematics(q.memptr(),link,origin,axis,amats);
	for (int i = 0; i < NCART; i++) {
		pos[i] = link[RIGHT_PALM][i+1];
		//normal[i] = amats[PALM][i+1][2];
	}
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
		eff_a[RIGHT_HAND][1]  = -M_PI/2.;
		eff_a[LEFT_HAND][1]   = -M_PI/2.;
		eff_a[RIGHT_HAND][3]  = -M_PI/2.;
		eff_a[LEFT_HAND][3]   = -M_PI/2.;
		eff_x[RIGHT_HAND][1]  = XHAND;
		eff_x[LEFT_HAND][1]   = XHAND;
		baseo[2] = 1.0;
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

	static double  Si00[3+1][3+1];
	static double  Si01[3+1][3+1];
	static double  Si12[3+1][3+1];
	static double  Si23[3+1][3+1];
	static double  Si34[3+1][3+1];
	static double  Si45[3+1][3+1];
	static double  Si56[3+1][3+1];
	static double  Si67[3+1][3+1];
	static double  Si78[3+1][3+1];
	static double  Si89[3+1][3+1];
	static double  Si910[3+1][3+1];
	static double  Si1011[3+1][3+1];
	static double  Si312[3+1][3+1];
	static double  Si1213[3+1][3+1];
	static double  Si1314[3+1][3+1];
	static double  Si1415[3+1][3+1];
	static double  Si1516[3+1][3+1];
	static double  Si1617[3+1][3+1];
	static double  Si1718[3+1][3+1];
	static double  Si1819[3+1][3+1];
	static double  Si320[3+1][3+1];
	static double  Si2021[3+1][3+1];
	static double  Si2122[3+1][3+1];
	static double  Si2223[3+1][3+1];
	static double  Si2324[3+1][3+1];
	static double  Si2226[3+1][3+1];
	static double  Si2627[3+1][3+1];
	static double  Si030[3+1][3+1];
	static double  Si3031[3+1][3+1];
	static double  Si3132[3+1][3+1];
	static double  Si3233[3+1][3+1];
	static double  Si3334[3+1][3+1];
	static double  Si3435[3+1][3+1];
	static double  Si3536[3+1][3+1];
	static double  Si3641[3+1][3+1];
	static double  Si042[3+1][3+1];
	static double  Si4243[3+1][3+1];
	static double  Si4344[3+1][3+1];
	static double  Si4445[3+1][3+1];
	static double  Si4546[3+1][3+1];
	static double  Si4647[3+1][3+1];
	static double  Si4748[3+1][3+1];
	static double  Si4853[3+1][3+1];

	static double  v[3+1];

	static double  vv[3+1];

	/* this function assumes that the array Xlink[nLinks+1][3+1]
	 is passed as an argument. Only the real number of links are computed */

	/* sine and cosine precomputation */
	sstate29th=Sin(state[29]);
	cstate29th=Cos(state[29]);

	sstate30th=Sin(state[30]);
	cstate30th=Cos(state[30]);

	sstate31th=Sin(state[31]);
	cstate31th=Cos(state[31]);

	sstate1th=Sin(state[1]);
	cstate1th=Cos(state[1]);

	sstate2th=Sin(state[2]);
	cstate2th=Cos(state[2]);

	sstate3th=Sin(state[3]);
	cstate3th=Cos(state[3]);

	sstate4th=Sin(state[4]);
	cstate4th=Cos(state[4]);

	sstate5th=Sin(state[5]);
	cstate5th=Cos(state[5]);

	sstate6th=Sin(state[6]);
	cstate6th=Cos(state[6]);

	sstate7th=Sin(state[7]);
	cstate7th=Cos(state[7]);

	sstate8th=Sin(state[8]);
	cstate8th=Cos(state[8]);

	sstate9th=Sin(state[9]);
	cstate9th=Cos(state[9]);

	sstate10th=Sin(state[10]);
	cstate10th=Cos(state[10]);

	sstate11th=Sin(state[11]);
	cstate11th=Cos(state[11]);

	sstate12th=Sin(state[12]);
	cstate12th=Cos(state[12]);

	sstate13th=Sin(state[13]);
	cstate13th=Cos(state[13]);

	sstate14th=Sin(state[14]);
	cstate14th=Cos(state[14]);

	sstate32th=Sin(state[32]);
	cstate32th=Cos(state[32]);

	sstate33th=Sin(state[33]);
	cstate33th=Cos(state[33]);

	sstate34th=Sin(state[34]);
	cstate34th=Cos(state[34]);

	sstate35th=Sin(state[35]);
	cstate35th=Cos(state[35]);

	sstate36th=Sin(state[36]);
	cstate36th=Cos(state[36]);

	sstate37th=Sin(state[37]);
	cstate37th=Cos(state[37]);

	sstate38th=Sin(state[38]);
	cstate38th=Cos(state[38]);

	sstate23th=Sin(state[23]);
	cstate23th=Cos(state[23]);

	sstate22th=Sin(state[22]);
	cstate22th=Cos(state[22]);

	sstate24th=Sin(state[24]);
	cstate24th=Cos(state[24]);

	sstate25th=Sin(state[25]);
	cstate25th=Cos(state[25]);

	sstate26th=Sin(state[26]);
	cstate26th=Cos(state[26]);

	sstate27th=Sin(state[27]);
	cstate27th=Cos(state[27]);

	sstate28th=Sin(state[28]);
	cstate28th=Cos(state[28]);

	sstate16th=Sin(state[16]);
	cstate16th=Cos(state[16]);

	sstate15th=Sin(state[15]);
	cstate15th=Cos(state[15]);

	sstate17th=Sin(state[17]);
	cstate17th=Cos(state[17]);

	sstate18th=Sin(state[18]);
	cstate18th=Cos(state[18]);

	sstate19th=Sin(state[19]);
	cstate19th=Cos(state[19]);

	sstate20th=Sin(state[20]);
	cstate20th=Cos(state[20]);

	sstate21th=Sin(state[21]);
	cstate21th=Cos(state[21]);


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



	/* inverse rotation matrices */
	Si00[1][1]=-1 + 2*Power(baseo[1],2) + 2*Power(baseo[2],2);
	Si00[1][2]=2*(baseo[2]*baseo[3] - baseo[1]*baseo[4]);
	Si00[1][3]=2*(baseo[1]*baseo[3] + baseo[2]*baseo[4]);

	Si00[2][1]=2*(baseo[2]*baseo[3] + baseo[1]*baseo[4]);
	Si00[2][2]=-1 + 2*Power(baseo[1],2) + 2*Power(baseo[3],2);
	Si00[2][3]=2*(-(baseo[1]*baseo[2]) + baseo[3]*baseo[4]);

	Si00[3][1]=2*(-(baseo[1]*baseo[3]) + baseo[2]*baseo[4]);
	Si00[3][2]=2*(baseo[1]*baseo[2] + baseo[3]*baseo[4]);
	Si00[3][3]=-1 + 2*Power(baseo[1],2) + 2*Power(baseo[4],2);


	Si01[1][1]=cstate29th;
	Si01[1][2]=-sstate29th;

	Si01[2][1]=-sstate29th;
	Si01[2][2]=-cstate29th;


	Si12[1][1]=-sstate30th;
	Si12[1][2]=-cstate30th;

	Si12[3][1]=cstate30th;
	Si12[3][2]=-sstate30th;


	Si23[1][1]=cstate31th;
	Si23[1][2]=-sstate31th;

	Si23[3][1]=-sstate31th;
	Si23[3][2]=-cstate31th;


	Si34[1][1]=0.7071067811865475*cstate1th;
	Si34[1][2]=-0.7071067811865475*sstate1th;

	Si34[2][1]=-sstate1th;
	Si34[2][2]=-cstate1th;

	Si34[3][1]=0.7071067811865475*cstate1th;
	Si34[3][2]=-0.7071067811865475*sstate1th;


	Si45[1][1]=-0.7071067811865475*cstate2th - 0.7071067811865475*sstate2th;
	Si45[1][2]=-0.7071067811865475*cstate2th + 0.7071067811865475*sstate2th;

	Si45[3][1]=0.7071067811865475*cstate2th - 0.7071067811865475*sstate2th;
	Si45[3][2]=-0.7071067811865475*cstate2th - 0.7071067811865475*sstate2th;


	Si56[1][1]=cstate3th;
	Si56[1][2]=-sstate3th;

	Si56[3][1]=-sstate3th;
	Si56[3][2]=-cstate3th;


	Si67[2][1]=cstate4th;
	Si67[2][2]=-sstate4th;

	Si67[3][1]=sstate4th;
	Si67[3][2]=cstate4th;


	Si78[1][1]=cstate5th;
	Si78[1][2]=-sstate5th;

	Si78[3][1]=-sstate5th;
	Si78[3][2]=-cstate5th;


	Si89[2][1]=sstate6th;
	Si89[2][2]=cstate6th;

	Si89[3][1]=-cstate6th;
	Si89[3][2]=sstate6th;


	Si910[1][1]=cstate7th;
	Si910[1][2]=-sstate7th;

	Si910[3][1]=sstate7th;
	Si910[3][2]=cstate7th;


	Si1011[1][1]=rceff2a2*rceff2a3;
	Si1011[1][2]=-(rceff2a2*rseff2a3);
	Si1011[1][3]=rseff2a2;

	Si1011[2][1]=rceff2a3*rseff2a1*rseff2a2 + rceff2a1*rseff2a3;
	Si1011[2][2]=rceff2a1*rceff2a3 - rseff2a1*rseff2a2*rseff2a3;
	Si1011[2][3]=-(rceff2a2*rseff2a1);

	Si1011[3][1]=-(rceff2a1*rceff2a3*rseff2a2) + rseff2a1*rseff2a3;
	Si1011[3][2]=rceff2a3*rseff2a1 + rceff2a1*rseff2a2*rseff2a3;
	Si1011[3][3]=rceff2a1*rceff2a2;


	Si312[1][1]=0.7071067811865475*cstate8th;
	Si312[1][2]=-0.7071067811865475*sstate8th;

	Si312[2][1]=-sstate8th;
	Si312[2][2]=-cstate8th;

	Si312[3][1]=-0.7071067811865475*cstate8th;
	Si312[3][2]=0.7071067811865475*sstate8th;


	Si1213[1][1]=-0.7071067811865475*cstate9th - 0.7071067811865475*sstate9th;
	Si1213[1][2]=-0.7071067811865475*cstate9th + 0.7071067811865475*sstate9th;

	Si1213[3][1]=-0.7071067811865475*cstate9th + 0.7071067811865475*sstate9th;
	Si1213[3][2]=0.7071067811865475*cstate9th + 0.7071067811865475*sstate9th;


	Si1314[1][1]=cstate10th;
	Si1314[1][2]=-sstate10th;

	Si1314[3][1]=sstate10th;
	Si1314[3][2]=cstate10th;


	Si1415[2][1]=cstate11th;
	Si1415[2][2]=-sstate11th;

	Si1415[3][1]=-sstate11th;
	Si1415[3][2]=-cstate11th;


	Si1516[1][1]=cstate12th;
	Si1516[1][2]=-sstate12th;

	Si1516[3][1]=sstate12th;
	Si1516[3][2]=cstate12th;


	Si1617[2][1]=sstate13th;
	Si1617[2][2]=cstate13th;

	Si1617[3][1]=cstate13th;
	Si1617[3][2]=-sstate13th;


	Si1718[1][1]=cstate14th;
	Si1718[1][2]=-sstate14th;

	Si1718[3][1]=-sstate14th;
	Si1718[3][2]=-cstate14th;


	Si1819[1][1]=rceff1a2*rceff1a3;
	Si1819[1][2]=-(rceff1a2*rseff1a3);
	Si1819[1][3]=rseff1a2;

	Si1819[2][1]=rceff1a3*rseff1a1*rseff1a2 + rceff1a1*rseff1a3;
	Si1819[2][2]=rceff1a1*rceff1a3 - rseff1a1*rseff1a2*rseff1a3;
	Si1819[2][3]=-(rceff1a2*rseff1a1);

	Si1819[3][1]=-(rceff1a1*rceff1a3*rseff1a2) + rseff1a1*rseff1a3;
	Si1819[3][2]=rceff1a3*rseff1a1 + rceff1a1*rseff1a2*rseff1a3;
	Si1819[3][3]=rceff1a1*rceff1a2;


	Si320[1][1]=sstate32th;
	Si320[1][2]=cstate32th;

	Si320[2][1]=-cstate32th;
	Si320[2][2]=sstate32th;


	Si2021[2][1]=sstate33th;
	Si2021[2][2]=cstate33th;

	Si2021[3][1]=-cstate33th;
	Si2021[3][2]=sstate33th;


	Si2122[1][1]=cstate34th;
	Si2122[1][2]=-sstate34th;

	Si2122[3][1]=-sstate34th;
	Si2122[3][2]=-cstate34th;


	Si2223[1][1]=cstate35th;
	Si2223[1][2]=-sstate35th;

	Si2223[2][1]=sstate35th;
	Si2223[2][2]=cstate35th;


	Si2324[2][1]=sstate36th;
	Si2324[2][2]=cstate36th;

	Si2324[3][1]=cstate36th;
	Si2324[3][2]=-sstate36th;



	Si2226[1][1]=cstate37th;
	Si2226[1][2]=-sstate37th;

	Si2226[2][1]=sstate37th;
	Si2226[2][2]=cstate37th;


	Si2627[2][1]=sstate38th;
	Si2627[2][2]=cstate38th;

	Si2627[3][1]=cstate38th;
	Si2627[3][2]=-sstate38th;




	Si030[1][1]=-sstate23th;
	Si030[1][2]=-cstate23th;

	Si030[3][1]=-cstate23th;
	Si030[3][2]=sstate23th;


	Si3031[1][1]=0.26720974913585105*cstate22th + 0.9636383917044586*sstate22th;
	Si3031[1][2]=0.9636383917044586*cstate22th - 0.26720974913585105*sstate22th;

	Si3031[3][1]=-0.9636383917044586*cstate22th + 0.26720974913585105*sstate22th;
	Si3031[3][2]=0.26720974913585105*cstate22th + 0.9636383917044586*sstate22th;


	Si3132[1][1]=cstate24th;
	Si3132[1][2]=-sstate24th;

	Si3132[3][1]=-sstate24th;
	Si3132[3][2]=-cstate24th;


	Si3233[1][1]=0.9636374101817434*cstate25th - 0.26721328877550665*sstate25th;
	Si3233[1][2]=-0.26721328877550665*cstate25th - 0.9636374101817434*sstate25th;

	Si3233[3][1]=-0.26721328877550665*cstate25th - 0.9636374101817434*sstate25th;
	Si3233[3][2]=-0.9636374101817434*cstate25th + 0.26721328877550665*sstate25th;


	Si3334[1][1]=cstate26th;
	Si3334[1][2]=-sstate26th;

	Si3334[3][1]=sstate26th;
	Si3334[3][2]=cstate26th;


	Si3435[1][1]=-sstate27th;
	Si3435[1][2]=-cstate27th;

	Si3435[3][1]=cstate27th;
	Si3435[3][2]=-sstate27th;


	Si3536[1][1]=cstate28th;
	Si3536[1][2]=-sstate28th;

	Si3536[3][1]=-sstate28th;
	Si3536[3][2]=-cstate28th;






	Si3641[1][1]=rceff3a2*rceff3a3;
	Si3641[1][2]=-(rceff3a2*rseff3a3);
	Si3641[1][3]=rseff3a2;

	Si3641[2][1]=rceff3a3*rseff3a1*rseff3a2 + rceff3a1*rseff3a3;
	Si3641[2][2]=rceff3a1*rceff3a3 - rseff3a1*rseff3a2*rseff3a3;
	Si3641[2][3]=-(rceff3a2*rseff3a1);

	Si3641[3][1]=-(rceff3a1*rceff3a3*rseff3a2) + rseff3a1*rseff3a3;
	Si3641[3][2]=rceff3a3*rseff3a1 + rceff3a1*rseff3a2*rseff3a3;
	Si3641[3][3]=rceff3a1*rceff3a2;


	Si042[1][1]=sstate16th;
	Si042[1][2]=cstate16th;

	Si042[3][1]=-cstate16th;
	Si042[3][2]=sstate16th;


	Si4243[1][1]=0.26720974913585105*cstate15th + 0.9636383917044586*sstate15th;
	Si4243[1][2]=0.9636383917044586*cstate15th - 0.26720974913585105*sstate15th;

	Si4243[3][1]=0.9636383917044586*cstate15th - 0.26720974913585105*sstate15th;
	Si4243[3][2]=-0.26720974913585105*cstate15th - 0.9636383917044586*sstate15th;


	Si4344[1][1]=cstate17th;
	Si4344[1][2]=-sstate17th;

	Si4344[3][1]=sstate17th;
	Si4344[3][2]=cstate17th;


	Si4445[1][1]=0.9636374101817434*cstate18th - 0.26721328877550665*sstate18th;
	Si4445[1][2]=-0.26721328877550665*cstate18th - 0.9636374101817434*sstate18th;

	Si4445[3][1]=0.26721328877550665*cstate18th + 0.9636374101817434*sstate18th;
	Si4445[3][2]=0.9636374101817434*cstate18th - 0.26721328877550665*sstate18th;


	Si4546[1][1]=cstate19th;
	Si4546[1][2]=-sstate19th;

	Si4546[3][1]=-sstate19th;
	Si4546[3][2]=-cstate19th;


	Si4647[1][1]=-sstate20th;
	Si4647[1][2]=-cstate20th;

	Si4647[3][1]=-cstate20th;
	Si4647[3][2]=sstate20th;


	Si4748[1][1]=cstate21th;
	Si4748[1][2]=-sstate21th;

	Si4748[3][1]=sstate21th;
	Si4748[3][2]=cstate21th;






	Si4853[1][1]=rceff4a2*rceff4a3;
	Si4853[1][2]=-(rceff4a2*rseff4a3);
	Si4853[1][3]=rseff4a2;

	Si4853[2][1]=rceff4a3*rseff4a1*rseff4a2 + rceff4a1*rseff4a3;
	Si4853[2][2]=rceff4a1*rceff4a3 - rseff4a1*rseff4a2*rseff4a3;
	Si4853[2][3]=-(rceff4a2*rseff4a1);

	Si4853[3][1]=-(rceff4a1*rceff4a3*rseff4a2) + rseff4a1*rseff4a3;
	Si4853[3][2]=rceff4a3*rseff4a1 + rceff4a1*rseff4a2*rseff4a3;
	Si4853[3][3]=rceff4a1*rceff4a2;



	/* calculation of link coordinates */
	/* link: {basec$0$$x[[1]], basec$0$$x[[2]], basec$0$$x[[3]]} */
	Xlink[0][1]=basec[1];
	Xlink[0][2]=basec[2];
	Xlink[0][3]=basec[3];

	v[1]=basec[1] - PELVISOFFSET*Si00[1][2] + PELVIS2THORAX*Si00[1][3];
	v[2]=basec[2] - PELVISOFFSET*Si00[2][2] + PELVIS2THORAX*Si00[2][3];
	v[3]=basec[3] - PELVISOFFSET*Si00[3][2] + PELVIS2THORAX*Si00[3][3];

	/* link: {0, -PELVISOFFSET, PELVIS2THORAX} */
	Xlink[1][1]=v[1];
	Xlink[1][2]=v[2];
	Xlink[1][3]=v[3];

	v[1]=-(THORAX2SHOULDER*Si23[1][1]);
	v[3]=-(THORAX2SHOULDER*Si23[3][1]);

	vv[1]=v[1]*Si12[1][1];
	vv[2]=-v[3];
	vv[3]=v[1]*Si12[3][1];

	v[1]=vv[1]*Si01[1][1] + vv[2]*Si01[1][2];
	v[2]=-PELVISOFFSET + vv[1]*Si01[2][1] + vv[2]*Si01[2][2];
	v[3]=PELVIS2THORAX - vv[3];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {-THORAX2SHOULDER, 0, 0} */
	Xlink[2][1]=vv[1];
	Xlink[2][2]=vv[2];
	Xlink[2][3]=vv[3];

	v[1]=-0.7071067811865475*SHOULDERX - THORAX2SHOULDER;
	v[3]=0.7071067811865475*SHOULDERX;

	vv[1]=v[1]*Si23[1][1];
	vv[2]=v[3];
	vv[3]=v[1]*Si23[3][1];

	v[1]=vv[1]*Si12[1][1] + vv[2]*Si12[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si12[3][1] + vv[2]*Si12[3][2];

	vv[1]=v[1]*Si01[1][1] + v[2]*Si01[1][2];
	vv[2]=-PELVISOFFSET + v[1]*Si01[2][1] + v[2]*Si01[2][2];
	vv[3]=PELVIS2THORAX - v[3];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {0, 0, -SHOULDERX} */
	Xlink[3][1]=v[1];
	Xlink[3][2]=v[2];
	Xlink[3][3]=v[3];

	v[1]=-(SHOULDERY*Si45[1][1]);
	v[3]=-SHOULDERX - SHOULDERY*Si45[3][1];

	vv[1]=-THORAX2SHOULDER + 0.7071067811865475*v[3] + v[1]*Si34[1][1];
	vv[2]=v[1]*Si34[2][1];
	vv[3]=-0.7071067811865475*v[3] + v[1]*Si34[3][1];

	v[1]=vv[1]*Si23[1][1] + vv[2]*Si23[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si23[3][1] + vv[2]*Si23[3][2];

	vv[1]=v[1]*Si12[1][1] + v[2]*Si12[1][2];
	vv[2]=-v[3];
	vv[3]=v[1]*Si12[3][1] + v[2]*Si12[3][2];

	v[1]=vv[1]*Si01[1][1] + vv[2]*Si01[1][2];
	v[2]=-PELVISOFFSET + vv[1]*Si01[2][1] + vv[2]*Si01[2][2];
	v[3]=PELVIS2THORAX - vv[3];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {-SHOULDERY, 0, 0} */
	Xlink[4][1]=vv[1];
	Xlink[4][2]=vv[2];
	Xlink[4][3]=vv[3];

	v[1]=-SHOULDERY;
	v[2]=-UPPERARM;

	vv[1]=v[1]*Si45[1][1] + v[2]*Si45[1][2];
	vv[3]=-SHOULDERX + v[1]*Si45[3][1] + v[2]*Si45[3][2];

	v[1]=-THORAX2SHOULDER + 0.7071067811865475*vv[3] + vv[1]*Si34[1][1];
	v[2]=vv[1]*Si34[2][1];
	v[3]=-0.7071067811865475*vv[3] + vv[1]*Si34[3][1];

	vv[1]=v[1]*Si23[1][1] + v[2]*Si23[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si23[3][1] + v[2]*Si23[3][2];

	v[1]=vv[1]*Si12[1][1] + vv[2]*Si12[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si12[3][1] + vv[2]*Si12[3][2];

	vv[1]=v[1]*Si01[1][1] + v[2]*Si01[1][2];
	vv[2]=-PELVISOFFSET + v[1]*Si01[2][1] + v[2]*Si01[2][2];
	vv[3]=PELVIS2THORAX - v[3];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {0, 0, -UPPERARM} */
	Xlink[5][1]=v[1];
	Xlink[5][2]=v[2];
	Xlink[5][3]=v[3];

	v[1]=WRISTY*Si78[1][2];
	v[2]=-LOWERARM;
	v[3]=WRISTY*Si78[3][2];

	vv[1]=v[3];
	vv[2]=v[1]*Si67[2][1] + v[2]*Si67[2][2];
	vv[3]=-UPPERARM + v[1]*Si67[3][1] + v[2]*Si67[3][2];

	v[1]=-SHOULDERY + vv[1]*Si56[1][1] + vv[2]*Si56[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si56[3][1] + vv[2]*Si56[3][2];

	vv[1]=v[1]*Si45[1][1] + v[2]*Si45[1][2];
	vv[2]=-v[3];
	vv[3]=-SHOULDERX + v[1]*Si45[3][1] + v[2]*Si45[3][2];

	v[1]=-THORAX2SHOULDER + 0.7071067811865475*vv[3] + vv[1]*Si34[1][1] + vv[2]*Si34[1][2];
	v[2]=vv[1]*Si34[2][1] + vv[2]*Si34[2][2];
	v[3]=-0.7071067811865475*vv[3] + vv[1]*Si34[3][1] + vv[2]*Si34[3][2];

	vv[1]=v[1]*Si23[1][1] + v[2]*Si23[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si23[3][1] + v[2]*Si23[3][2];

	v[1]=vv[1]*Si12[1][1] + vv[2]*Si12[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si12[3][1] + vv[2]*Si12[3][2];

	vv[1]=v[1]*Si01[1][1] + v[2]*Si01[1][2];
	vv[2]=-PELVISOFFSET + v[1]*Si01[2][1] + v[2]*Si01[2][2];
	vv[3]=PELVIS2THORAX - v[3];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {0, WRISTY, -LOWERARM} */
	Xlink[6][1]=v[1];
	Xlink[6][2]=v[2];
	Xlink[6][3]=v[3];

	v[1]=eff_x[2][1]*Si910[1][1] + eff_x[2][2]*Si910[1][2];
	v[2]=-eff_x[2][3];
	v[3]=eff_x[2][1]*Si910[3][1] + eff_x[2][2]*Si910[3][2];

	vv[1]=v[3];
	vv[2]=WRISTY + v[1]*Si89[2][1] + v[2]*Si89[2][2];
	vv[3]=-LOWERARM + v[1]*Si89[3][1] + v[2]*Si89[3][2];

	v[1]=vv[1]*Si78[1][1] + vv[2]*Si78[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si78[3][1] + vv[2]*Si78[3][2];

	vv[1]=v[3];
	vv[2]=v[1]*Si67[2][1] + v[2]*Si67[2][2];
	vv[3]=-UPPERARM + v[1]*Si67[3][1] + v[2]*Si67[3][2];

	v[1]=-SHOULDERY + vv[1]*Si56[1][1] + vv[2]*Si56[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si56[3][1] + vv[2]*Si56[3][2];

	vv[1]=v[1]*Si45[1][1] + v[2]*Si45[1][2];
	vv[2]=-v[3];
	vv[3]=-SHOULDERX + v[1]*Si45[3][1] + v[2]*Si45[3][2];

	v[1]=-THORAX2SHOULDER + 0.7071067811865475*vv[3] + vv[1]*Si34[1][1] + vv[2]*Si34[1][2];
	v[2]=vv[1]*Si34[2][1] + vv[2]*Si34[2][2];
	v[3]=-0.7071067811865475*vv[3] + vv[1]*Si34[3][1] + vv[2]*Si34[3][2];

	vv[1]=v[1]*Si23[1][1] + v[2]*Si23[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si23[3][1] + v[2]*Si23[3][2];

	v[1]=vv[1]*Si12[1][1] + vv[2]*Si12[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si12[3][1] + vv[2]*Si12[3][2];

	vv[1]=v[1]*Si01[1][1] + v[2]*Si01[1][2];
	vv[2]=-PELVISOFFSET + v[1]*Si01[2][1] + v[2]*Si01[2][2];
	vv[3]=PELVIS2THORAX - v[3];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {eff$2$$x[[1]], eff$2$$x[[2]], eff$2$$x[[3]]} */
	Xlink[7][1]=v[1];
	Xlink[7][2]=v[2];
	Xlink[7][3]=v[3];

	v[1]=-(THORAX2SHOULDER*Si23[1][1]);
	v[3]=-(THORAX2SHOULDER*Si23[3][1]);

	vv[1]=v[1]*Si12[1][1];
	vv[2]=-v[3];
	vv[3]=v[1]*Si12[3][1];

	v[1]=vv[1]*Si01[1][1] + vv[2]*Si01[1][2];
	v[2]=-PELVISOFFSET + vv[1]*Si01[2][1] + vv[2]*Si01[2][2];
	v[3]=PELVIS2THORAX - vv[3];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {-THORAX2SHOULDER, 0, 0} */
	Xlink[8][1]=vv[1];
	Xlink[8][2]=vv[2];
	Xlink[8][3]=vv[3];

	v[1]=-0.7071067811865475*SHOULDERX - THORAX2SHOULDER;
	v[3]=-0.7071067811865475*SHOULDERX;

	vv[1]=v[1]*Si23[1][1];
	vv[2]=v[3];
	vv[3]=v[1]*Si23[3][1];

	v[1]=vv[1]*Si12[1][1] + vv[2]*Si12[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si12[3][1] + vv[2]*Si12[3][2];

	vv[1]=v[1]*Si01[1][1] + v[2]*Si01[1][2];
	vv[2]=-PELVISOFFSET + v[1]*Si01[2][1] + v[2]*Si01[2][2];
	vv[3]=PELVIS2THORAX - v[3];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {0, 0, SHOULDERX} */
	Xlink[9][1]=v[1];
	Xlink[9][2]=v[2];
	Xlink[9][3]=v[3];

	v[1]=-(SHOULDERY*Si1213[1][1]);
	v[3]=SHOULDERX - SHOULDERY*Si1213[3][1];

	vv[1]=-THORAX2SHOULDER - 0.7071067811865475*v[3] + v[1]*Si312[1][1];
	vv[2]=v[1]*Si312[2][1];
	vv[3]=-0.7071067811865475*v[3] + v[1]*Si312[3][1];

	v[1]=vv[1]*Si23[1][1] + vv[2]*Si23[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si23[3][1] + vv[2]*Si23[3][2];

	vv[1]=v[1]*Si12[1][1] + v[2]*Si12[1][2];
	vv[2]=-v[3];
	vv[3]=v[1]*Si12[3][1] + v[2]*Si12[3][2];

	v[1]=vv[1]*Si01[1][1] + vv[2]*Si01[1][2];
	v[2]=-PELVISOFFSET + vv[1]*Si01[2][1] + vv[2]*Si01[2][2];
	v[3]=PELVIS2THORAX - vv[3];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {-SHOULDERY, 0, 0} */
	Xlink[10][1]=vv[1];
	Xlink[10][2]=vv[2];
	Xlink[10][3]=vv[3];

	v[1]=-SHOULDERY;
	v[2]=-UPPERARM;

	vv[1]=v[1]*Si1213[1][1] + v[2]*Si1213[1][2];
	vv[3]=SHOULDERX + v[1]*Si1213[3][1] + v[2]*Si1213[3][2];

	v[1]=-THORAX2SHOULDER - 0.7071067811865475*vv[3] + vv[1]*Si312[1][1];
	v[2]=vv[1]*Si312[2][1];
	v[3]=-0.7071067811865475*vv[3] + vv[1]*Si312[3][1];

	vv[1]=v[1]*Si23[1][1] + v[2]*Si23[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si23[3][1] + v[2]*Si23[3][2];

	v[1]=vv[1]*Si12[1][1] + vv[2]*Si12[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si12[3][1] + vv[2]*Si12[3][2];

	vv[1]=v[1]*Si01[1][1] + v[2]*Si01[1][2];
	vv[2]=-PELVISOFFSET + v[1]*Si01[2][1] + v[2]*Si01[2][2];
	vv[3]=PELVIS2THORAX - v[3];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {0, 0, UPPERARM} */
	Xlink[11][1]=v[1];
	Xlink[11][2]=v[2];
	Xlink[11][3]=v[3];

	v[1]=WRISTY*Si1516[1][2];
	v[2]=-LOWERARM;
	v[3]=WRISTY*Si1516[3][2];

	vv[1]=-v[3];
	vv[2]=v[1]*Si1415[2][1] + v[2]*Si1415[2][2];
	vv[3]=UPPERARM + v[1]*Si1415[3][1] + v[2]*Si1415[3][2];

	v[1]=-SHOULDERY + vv[1]*Si1314[1][1] + vv[2]*Si1314[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si1314[3][1] + vv[2]*Si1314[3][2];

	vv[1]=v[1]*Si1213[1][1] + v[2]*Si1213[1][2];
	vv[2]=v[3];
	vv[3]=SHOULDERX + v[1]*Si1213[3][1] + v[2]*Si1213[3][2];

	v[1]=-THORAX2SHOULDER - 0.7071067811865475*vv[3] + vv[1]*Si312[1][1] + vv[2]*Si312[1][2];
	v[2]=vv[1]*Si312[2][1] + vv[2]*Si312[2][2];
	v[3]=-0.7071067811865475*vv[3] + vv[1]*Si312[3][1] + vv[2]*Si312[3][2];

	vv[1]=v[1]*Si23[1][1] + v[2]*Si23[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si23[3][1] + v[2]*Si23[3][2];

	v[1]=vv[1]*Si12[1][1] + vv[2]*Si12[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si12[3][1] + vv[2]*Si12[3][2];

	vv[1]=v[1]*Si01[1][1] + v[2]*Si01[1][2];
	vv[2]=-PELVISOFFSET + v[1]*Si01[2][1] + v[2]*Si01[2][2];
	vv[3]=PELVIS2THORAX - v[3];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {0, WRISTY, LOWERARM} */
	Xlink[12][1]=v[1];
	Xlink[12][2]=v[2];
	Xlink[12][3]=v[3];

	v[1]=eff_x[1][1]*Si1718[1][1] + eff_x[1][2]*Si1718[1][2];
	v[2]=eff_x[1][3];
	v[3]=eff_x[1][1]*Si1718[3][1] + eff_x[1][2]*Si1718[3][2];

	vv[1]=-v[3];
	vv[2]=WRISTY + v[1]*Si1617[2][1] + v[2]*Si1617[2][2];
	vv[3]=LOWERARM + v[1]*Si1617[3][1] + v[2]*Si1617[3][2];

	v[1]=vv[1]*Si1516[1][1] + vv[2]*Si1516[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si1516[3][1] + vv[2]*Si1516[3][2];

	vv[1]=-v[3];
	vv[2]=v[1]*Si1415[2][1] + v[2]*Si1415[2][2];
	vv[3]=UPPERARM + v[1]*Si1415[3][1] + v[2]*Si1415[3][2];

	v[1]=-SHOULDERY + vv[1]*Si1314[1][1] + vv[2]*Si1314[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si1314[3][1] + vv[2]*Si1314[3][2];

	vv[1]=v[1]*Si1213[1][1] + v[2]*Si1213[1][2];
	vv[2]=v[3];
	vv[3]=SHOULDERX + v[1]*Si1213[3][1] + v[2]*Si1213[3][2];

	v[1]=-THORAX2SHOULDER - 0.7071067811865475*vv[3] + vv[1]*Si312[1][1] + vv[2]*Si312[1][2];
	v[2]=vv[1]*Si312[2][1] + vv[2]*Si312[2][2];
	v[3]=-0.7071067811865475*vv[3] + vv[1]*Si312[3][1] + vv[2]*Si312[3][2];

	vv[1]=v[1]*Si23[1][1] + v[2]*Si23[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si23[3][1] + v[2]*Si23[3][2];

	v[1]=vv[1]*Si12[1][1] + vv[2]*Si12[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si12[3][1] + vv[2]*Si12[3][2];

	vv[1]=v[1]*Si01[1][1] + v[2]*Si01[1][2];
	vv[2]=-PELVISOFFSET + v[1]*Si01[2][1] + v[2]*Si01[2][2];
	vv[3]=PELVIS2THORAX - v[3];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {eff$1$$x[[1]], eff$1$$x[[2]], eff$1$$x[[3]]} */
	Xlink[13][1]=v[1];
	Xlink[13][2]=v[2];
	Xlink[13][3]=v[3];

	v[1]=-(THORAX2NECK*Si23[1][1]);
	v[3]=-(THORAX2NECK*Si23[3][1]);

	vv[1]=v[1]*Si12[1][1];
	vv[2]=-v[3];
	vv[3]=v[1]*Si12[3][1];

	v[1]=vv[1]*Si01[1][1] + vv[2]*Si01[1][2];
	v[2]=-PELVISOFFSET + vv[1]*Si01[2][1] + vv[2]*Si01[2][2];
	v[3]=PELVIS2THORAX - vv[3];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {-THORAX2NECK, 0, 0} */
	Xlink[14][1]=vv[1];
	Xlink[14][2]=vv[2];
	Xlink[14][3]=vv[3];

	v[1]=-THORAX2NECK - CERVICAL*Si320[1][2];
	v[2]=-(CERVICAL*Si320[2][2]);

	vv[1]=v[1]*Si23[1][1] + v[2]*Si23[1][2];
	vv[3]=v[1]*Si23[3][1] + v[2]*Si23[3][2];

	v[1]=vv[1]*Si12[1][1];
	v[2]=-vv[3];
	v[3]=vv[1]*Si12[3][1];

	vv[1]=v[1]*Si01[1][1] + v[2]*Si01[1][2];
	vv[2]=-PELVISOFFSET + v[1]*Si01[2][1] + v[2]*Si01[2][2];
	vv[3]=PELVIS2THORAX - v[3];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {0, -CERVICAL, 0} */
	Xlink[15][1]=v[1];
	Xlink[15][2]=v[2];
	Xlink[15][3]=v[3];

	v[1]=EYEXOFF*Si2122[1][1] - EYEYOFF*Si2122[1][2];
	v[2]=-HEAD;
	v[3]=EYEXOFF*Si2122[3][1] - EYEYOFF*Si2122[3][2];

	vv[1]=v[3];
	vv[2]=-CERVICAL + v[1]*Si2021[2][1] + v[2]*Si2021[2][2];
	vv[3]=v[1]*Si2021[3][1] + v[2]*Si2021[3][2];

	v[1]=-THORAX2NECK + vv[1]*Si320[1][1] + vv[2]*Si320[1][2];
	v[2]=vv[1]*Si320[2][1] + vv[2]*Si320[2][2];
	v[3]=vv[3];

	vv[1]=v[1]*Si23[1][1] + v[2]*Si23[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si23[3][1] + v[2]*Si23[3][2];

	v[1]=vv[1]*Si12[1][1] + vv[2]*Si12[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si12[3][1] + vv[2]*Si12[3][2];

	vv[1]=v[1]*Si01[1][1] + v[2]*Si01[1][2];
	vv[2]=-PELVISOFFSET + v[1]*Si01[2][1] + v[2]*Si01[2][2];
	vv[3]=PELVIS2THORAX - v[3];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {EYEXOFF, -EYEYOFF, -HEAD} */
	Xlink[16][1]=v[1];
	Xlink[16][2]=v[2];
	Xlink[16][3]=v[3];

	v[2]=-(EYE*Si2324[2][2]);
	v[3]=-(EYE*Si2324[3][2]);

	vv[1]=EYEXOFF + v[2]*Si2223[1][2];
	vv[2]=-EYEYOFF + v[2]*Si2223[2][2];
	vv[3]=-HEAD + v[3];

	v[1]=vv[1]*Si2122[1][1] + vv[2]*Si2122[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si2122[3][1] + vv[2]*Si2122[3][2];

	vv[1]=v[3];
	vv[2]=-CERVICAL + v[1]*Si2021[2][1] + v[2]*Si2021[2][2];
	vv[3]=v[1]*Si2021[3][1] + v[2]*Si2021[3][2];

	v[1]=-THORAX2NECK + vv[1]*Si320[1][1] + vv[2]*Si320[1][2];
	v[2]=vv[1]*Si320[2][1] + vv[2]*Si320[2][2];
	v[3]=vv[3];

	vv[1]=v[1]*Si23[1][1] + v[2]*Si23[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si23[3][1] + v[2]*Si23[3][2];

	v[1]=vv[1]*Si12[1][1] + vv[2]*Si12[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si12[3][1] + vv[2]*Si12[3][2];

	vv[1]=v[1]*Si01[1][1] + v[2]*Si01[1][2];
	vv[2]=-PELVISOFFSET + v[1]*Si01[2][1] + v[2]*Si01[2][2];
	vv[3]=PELVIS2THORAX - v[3];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {0, -EYE, 0} */
	Xlink[17][1]=v[1];
	Xlink[17][2]=v[2];
	Xlink[17][3]=v[3];

	v[1]=-(EYEXOFF*Si2122[1][1]) - EYEYOFF*Si2122[1][2];
	v[2]=-HEAD;
	v[3]=-(EYEXOFF*Si2122[3][1]) - EYEYOFF*Si2122[3][2];

	vv[1]=v[3];
	vv[2]=-CERVICAL + v[1]*Si2021[2][1] + v[2]*Si2021[2][2];
	vv[3]=v[1]*Si2021[3][1] + v[2]*Si2021[3][2];

	v[1]=-THORAX2NECK + vv[1]*Si320[1][1] + vv[2]*Si320[1][2];
	v[2]=vv[1]*Si320[2][1] + vv[2]*Si320[2][2];
	v[3]=vv[3];

	vv[1]=v[1]*Si23[1][1] + v[2]*Si23[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si23[3][1] + v[2]*Si23[3][2];

	v[1]=vv[1]*Si12[1][1] + vv[2]*Si12[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si12[3][1] + vv[2]*Si12[3][2];

	vv[1]=v[1]*Si01[1][1] + v[2]*Si01[1][2];
	vv[2]=-PELVISOFFSET + v[1]*Si01[2][1] + v[2]*Si01[2][2];
	vv[3]=PELVIS2THORAX - v[3];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {-EYEXOFF, -EYEYOFF, -HEAD} */
	Xlink[18][1]=v[1];
	Xlink[18][2]=v[2];
	Xlink[18][3]=v[3];

	v[2]=-(EYE*Si2627[2][2]);
	v[3]=-(EYE*Si2627[3][2]);

	vv[1]=-EYEXOFF + v[2]*Si2226[1][2];
	vv[2]=-EYEYOFF + v[2]*Si2226[2][2];
	vv[3]=-HEAD + v[3];

	v[1]=vv[1]*Si2122[1][1] + vv[2]*Si2122[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si2122[3][1] + vv[2]*Si2122[3][2];

	vv[1]=v[3];
	vv[2]=-CERVICAL + v[1]*Si2021[2][1] + v[2]*Si2021[2][2];
	vv[3]=v[1]*Si2021[3][1] + v[2]*Si2021[3][2];

	v[1]=-THORAX2NECK + vv[1]*Si320[1][1] + vv[2]*Si320[1][2];
	v[2]=vv[1]*Si320[2][1] + vv[2]*Si320[2][2];
	v[3]=vv[3];

	vv[1]=v[1]*Si23[1][1] + v[2]*Si23[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si23[3][1] + v[2]*Si23[3][2];

	v[1]=vv[1]*Si12[1][1] + vv[2]*Si12[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si12[3][1] + vv[2]*Si12[3][2];

	vv[1]=v[1]*Si01[1][1] + v[2]*Si01[1][2];
	vv[2]=-PELVISOFFSET + v[1]*Si01[2][1] + v[2]*Si01[2][2];
	vv[3]=PELVIS2THORAX - v[3];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {0, -EYE, 0} */
	Xlink[19][1]=v[1];
	Xlink[19][2]=v[2];
	Xlink[19][3]=v[3];

	v[2]=-TOPofHEAD;

	vv[2]=-CERVICAL + v[2]*Si2021[2][2];
	vv[3]=v[2]*Si2021[3][2];

	v[1]=-THORAX2NECK + vv[2]*Si320[1][2];
	v[2]=vv[2]*Si320[2][2];
	v[3]=vv[3];

	vv[1]=v[1]*Si23[1][1] + v[2]*Si23[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si23[3][1] + v[2]*Si23[3][2];

	v[1]=vv[1]*Si12[1][1] + vv[2]*Si12[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si12[3][1] + vv[2]*Si12[3][2];

	vv[1]=v[1]*Si01[1][1] + v[2]*Si01[1][2];
	vv[2]=-PELVISOFFSET + v[1]*Si01[2][1] + v[2]*Si01[2][2];
	vv[3]=PELVIS2THORAX - v[3];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {0, 0, -TOPofHEAD} */
	Xlink[20][1]=v[1];
	Xlink[20][2]=v[2];
	Xlink[20][3]=v[3];

	v[1]=basec[1] + XHIP*Si00[1][1];
	v[2]=basec[2] + XHIP*Si00[2][1];
	v[3]=basec[3] + XHIP*Si00[3][1];

	/* link: {XHIP, 0, 0} */
	Xlink[21][1]=v[1];
	Xlink[21][2]=v[2];
	Xlink[21][3]=v[3];

	v[1]=YHIP*Si3031[1][1];
	v[3]=YHIP*Si3031[3][1];

	vv[1]=XHIP + v[1]*Si030[1][1];
	vv[2]=v[3];
	vv[3]=v[1]*Si030[3][1];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {YHIP, 0, 0} */
	Xlink[22][1]=v[1];
	Xlink[22][2]=v[2];
	Xlink[22][3]=v[3];

	v[1]=YHIP + YKNEE*Si3132[1][1];
	v[2]=UPPERLEG;
	v[3]=YKNEE*Si3132[3][1];

	vv[1]=v[1]*Si3031[1][1] + v[2]*Si3031[1][2];
	vv[2]=-v[3];
	vv[3]=v[1]*Si3031[3][1] + v[2]*Si3031[3][2];

	v[1]=XHIP + vv[1]*Si030[1][1] + vv[2]*Si030[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si030[3][1] + vv[2]*Si030[3][2];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {YKNEE, 0, UPPERLEG} */
	Xlink[23][1]=vv[1];
	Xlink[23][2]=vv[2];
	Xlink[23][3]=vv[3];

	v[2]=-LOWERLEG;

	vv[1]=YKNEE + v[2]*Si3233[1][2];
	vv[3]=UPPERLEG + v[2]*Si3233[3][2];

	v[1]=YHIP + vv[1]*Si3132[1][1];
	v[2]=vv[3];
	v[3]=vv[1]*Si3132[3][1];

	vv[1]=v[1]*Si3031[1][1] + v[2]*Si3031[1][2];
	vv[2]=-v[3];
	vv[3]=v[1]*Si3031[3][1] + v[2]*Si3031[3][2];

	v[1]=XHIP + vv[1]*Si030[1][1] + vv[2]*Si030[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si030[3][1] + vv[2]*Si030[3][2];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {0, 0, LOWERLEG} */
	Xlink[24][1]=vv[1];
	Xlink[24][2]=vv[2];
	Xlink[24][3]=vv[3];

	v[1]=ZTOE*Si3536[1][1] - XTOE*Si3536[1][2];
	v[2]=YTOE;
	v[3]=ZTOE*Si3536[3][1] - XTOE*Si3536[3][2];

	vv[1]=v[1]*Si3435[1][1] + v[2]*Si3435[1][2];
	vv[2]=-v[3];
	vv[3]=LOWERLEG + v[1]*Si3435[3][1] + v[2]*Si3435[3][2];

	v[1]=vv[1]*Si3334[1][1] + vv[2]*Si3334[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si3334[3][1] + vv[2]*Si3334[3][2];

	vv[1]=YKNEE + v[1]*Si3233[1][1] + v[2]*Si3233[1][2];
	vv[2]=v[3];
	vv[3]=UPPERLEG + v[1]*Si3233[3][1] + v[2]*Si3233[3][2];

	v[1]=YHIP + vv[1]*Si3132[1][1] + vv[2]*Si3132[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si3132[3][1] + vv[2]*Si3132[3][2];

	vv[1]=v[1]*Si3031[1][1] + v[2]*Si3031[1][2];
	vv[2]=-v[3];
	vv[3]=v[1]*Si3031[3][1] + v[2]*Si3031[3][2];

	v[1]=XHIP + vv[1]*Si030[1][1] + vv[2]*Si030[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si030[3][1] + vv[2]*Si030[3][2];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {ZTOE, -XTOE, YTOE} */
	Xlink[25][1]=vv[1];
	Xlink[25][2]=vv[2];
	Xlink[25][3]=vv[3];

	v[1]=ZTOE*Si3536[1][1] + XTOE*Si3536[1][2];
	v[2]=YTOE;
	v[3]=ZTOE*Si3536[3][1] + XTOE*Si3536[3][2];

	vv[1]=v[1]*Si3435[1][1] + v[2]*Si3435[1][2];
	vv[2]=-v[3];
	vv[3]=LOWERLEG + v[1]*Si3435[3][1] + v[2]*Si3435[3][2];

	v[1]=vv[1]*Si3334[1][1] + vv[2]*Si3334[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si3334[3][1] + vv[2]*Si3334[3][2];

	vv[1]=YKNEE + v[1]*Si3233[1][1] + v[2]*Si3233[1][2];
	vv[2]=v[3];
	vv[3]=UPPERLEG + v[1]*Si3233[3][1] + v[2]*Si3233[3][2];

	v[1]=YHIP + vv[1]*Si3132[1][1] + vv[2]*Si3132[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si3132[3][1] + vv[2]*Si3132[3][2];

	vv[1]=v[1]*Si3031[1][1] + v[2]*Si3031[1][2];
	vv[2]=-v[3];
	vv[3]=v[1]*Si3031[3][1] + v[2]*Si3031[3][2];

	v[1]=XHIP + vv[1]*Si030[1][1] + vv[2]*Si030[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si030[3][1] + vv[2]*Si030[3][2];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {ZTOE, XTOE, YTOE} */
	Xlink[26][1]=vv[1];
	Xlink[26][2]=vv[2];
	Xlink[26][3]=vv[3];

	v[1]=ZHEEL*Si3536[1][1] - XHEEL*Si3536[1][2];
	v[2]=-YHEEL;
	v[3]=ZHEEL*Si3536[3][1] - XHEEL*Si3536[3][2];

	vv[1]=v[1]*Si3435[1][1] + v[2]*Si3435[1][2];
	vv[2]=-v[3];
	vv[3]=LOWERLEG + v[1]*Si3435[3][1] + v[2]*Si3435[3][2];

	v[1]=vv[1]*Si3334[1][1] + vv[2]*Si3334[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si3334[3][1] + vv[2]*Si3334[3][2];

	vv[1]=YKNEE + v[1]*Si3233[1][1] + v[2]*Si3233[1][2];
	vv[2]=v[3];
	vv[3]=UPPERLEG + v[1]*Si3233[3][1] + v[2]*Si3233[3][2];

	v[1]=YHIP + vv[1]*Si3132[1][1] + vv[2]*Si3132[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si3132[3][1] + vv[2]*Si3132[3][2];

	vv[1]=v[1]*Si3031[1][1] + v[2]*Si3031[1][2];
	vv[2]=-v[3];
	vv[3]=v[1]*Si3031[3][1] + v[2]*Si3031[3][2];

	v[1]=XHIP + vv[1]*Si030[1][1] + vv[2]*Si030[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si030[3][1] + vv[2]*Si030[3][2];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {ZHEEL, -XHEEL, -YHEEL} */
	Xlink[27][1]=vv[1];
	Xlink[27][2]=vv[2];
	Xlink[27][3]=vv[3];

	v[1]=ZHEEL*Si3536[1][1] + XHEEL*Si3536[1][2];
	v[2]=-YHEEL;
	v[3]=ZHEEL*Si3536[3][1] + XHEEL*Si3536[3][2];

	vv[1]=v[1]*Si3435[1][1] + v[2]*Si3435[1][2];
	vv[2]=-v[3];
	vv[3]=LOWERLEG + v[1]*Si3435[3][1] + v[2]*Si3435[3][2];

	v[1]=vv[1]*Si3334[1][1] + vv[2]*Si3334[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si3334[3][1] + vv[2]*Si3334[3][2];

	vv[1]=YKNEE + v[1]*Si3233[1][1] + v[2]*Si3233[1][2];
	vv[2]=v[3];
	vv[3]=UPPERLEG + v[1]*Si3233[3][1] + v[2]*Si3233[3][2];

	v[1]=YHIP + vv[1]*Si3132[1][1] + vv[2]*Si3132[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si3132[3][1] + vv[2]*Si3132[3][2];

	vv[1]=v[1]*Si3031[1][1] + v[2]*Si3031[1][2];
	vv[2]=-v[3];
	vv[3]=v[1]*Si3031[3][1] + v[2]*Si3031[3][2];

	v[1]=XHIP + vv[1]*Si030[1][1] + vv[2]*Si030[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si030[3][1] + vv[2]*Si030[3][2];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {ZHEEL, XHEEL, -YHEEL} */
	Xlink[28][1]=vv[1];
	Xlink[28][2]=vv[2];
	Xlink[28][3]=vv[3];

	v[1]=eff_x[3][1]*Si3536[1][1] + eff_x[3][2]*Si3536[1][2];
	v[2]=eff_x[3][3];
	v[3]=eff_x[3][1]*Si3536[3][1] + eff_x[3][2]*Si3536[3][2];

	vv[1]=v[1]*Si3435[1][1] + v[2]*Si3435[1][2];
	vv[2]=-v[3];
	vv[3]=LOWERLEG + v[1]*Si3435[3][1] + v[2]*Si3435[3][2];

	v[1]=vv[1]*Si3334[1][1] + vv[2]*Si3334[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si3334[3][1] + vv[2]*Si3334[3][2];

	vv[1]=YKNEE + v[1]*Si3233[1][1] + v[2]*Si3233[1][2];
	vv[2]=v[3];
	vv[3]=UPPERLEG + v[1]*Si3233[3][1] + v[2]*Si3233[3][2];

	v[1]=YHIP + vv[1]*Si3132[1][1] + vv[2]*Si3132[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si3132[3][1] + vv[2]*Si3132[3][2];

	vv[1]=v[1]*Si3031[1][1] + v[2]*Si3031[1][2];
	vv[2]=-v[3];
	vv[3]=v[1]*Si3031[3][1] + v[2]*Si3031[3][2];

	v[1]=XHIP + vv[1]*Si030[1][1] + vv[2]*Si030[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si030[3][1] + vv[2]*Si030[3][2];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {eff$3$$x[[1]], eff$3$$x[[2]], eff$3$$x[[3]]} */
	Xlink[29][1]=vv[1];
	Xlink[29][2]=vv[2];
	Xlink[29][3]=vv[3];

	v[1]=basec[1] - XHIP*Si00[1][1];
	v[2]=basec[2] - XHIP*Si00[2][1];
	v[3]=basec[3] - XHIP*Si00[3][1];

	/* link: {-XHIP, 0, 0} */
	Xlink[30][1]=v[1];
	Xlink[30][2]=v[2];
	Xlink[30][3]=v[3];

	v[1]=YHIP*Si4243[1][1];
	v[3]=YHIP*Si4243[3][1];

	vv[1]=-XHIP + v[1]*Si042[1][1];
	vv[2]=-v[3];
	vv[3]=v[1]*Si042[3][1];

	v[1]=basec[1] + vv[1]*Si00[1][1] + vv[2]*Si00[1][2] + vv[3]*Si00[1][3];
	v[2]=basec[2] + vv[1]*Si00[2][1] + vv[2]*Si00[2][2] + vv[3]*Si00[2][3];
	v[3]=basec[3] + vv[1]*Si00[3][1] + vv[2]*Si00[3][2] + vv[3]*Si00[3][3];

	/* link: {YHIP, 0, 0} */
	Xlink[31][1]=v[1];
	Xlink[31][2]=v[2];
	Xlink[31][3]=v[3];

	v[1]=YHIP + YKNEE*Si4344[1][1];
	v[2]=UPPERLEG;
	v[3]=YKNEE*Si4344[3][1];

	vv[1]=v[1]*Si4243[1][1] + v[2]*Si4243[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si4243[3][1] + v[2]*Si4243[3][2];

	v[1]=-XHIP + vv[1]*Si042[1][1] + vv[2]*Si042[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si042[3][1] + vv[2]*Si042[3][2];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {YKNEE, 0, -UPPERLEG} */
	Xlink[32][1]=vv[1];
	Xlink[32][2]=vv[2];
	Xlink[32][3]=vv[3];

	v[2]=-LOWERLEG;

	vv[1]=YKNEE + v[2]*Si4445[1][2];
	vv[3]=-UPPERLEG + v[2]*Si4445[3][2];

	v[1]=YHIP + vv[1]*Si4344[1][1];
	v[2]=-vv[3];
	v[3]=vv[1]*Si4344[3][1];

	vv[1]=v[1]*Si4243[1][1] + v[2]*Si4243[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si4243[3][1] + v[2]*Si4243[3][2];

	v[1]=-XHIP + vv[1]*Si042[1][1] + vv[2]*Si042[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si042[3][1] + vv[2]*Si042[3][2];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {0, 0, -LOWERLEG} */
	Xlink[33][1]=vv[1];
	Xlink[33][2]=vv[2];
	Xlink[33][3]=vv[3];

	v[1]=ZTOE*Si4748[1][1] + XTOE*Si4748[1][2];
	v[2]=YTOE;
	v[3]=ZTOE*Si4748[3][1] + XTOE*Si4748[3][2];

	vv[1]=v[1]*Si4647[1][1] + v[2]*Si4647[1][2];
	vv[2]=v[3];
	vv[3]=-LOWERLEG + v[1]*Si4647[3][1] + v[2]*Si4647[3][2];

	v[1]=vv[1]*Si4546[1][1] + vv[2]*Si4546[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si4546[3][1] + vv[2]*Si4546[3][2];

	vv[1]=YKNEE + v[1]*Si4445[1][1] + v[2]*Si4445[1][2];
	vv[2]=-v[3];
	vv[3]=-UPPERLEG + v[1]*Si4445[3][1] + v[2]*Si4445[3][2];

	v[1]=YHIP + vv[1]*Si4344[1][1] + vv[2]*Si4344[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si4344[3][1] + vv[2]*Si4344[3][2];

	vv[1]=v[1]*Si4243[1][1] + v[2]*Si4243[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si4243[3][1] + v[2]*Si4243[3][2];

	v[1]=-XHIP + vv[1]*Si042[1][1] + vv[2]*Si042[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si042[3][1] + vv[2]*Si042[3][2];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {ZTOE, XTOE, -YTOE} */
	Xlink[34][1]=vv[1];
	Xlink[34][2]=vv[2];
	Xlink[34][3]=vv[3];

	v[1]=ZTOE*Si4748[1][1] - XTOE*Si4748[1][2];
	v[2]=YTOE;
	v[3]=ZTOE*Si4748[3][1] - XTOE*Si4748[3][2];

	vv[1]=v[1]*Si4647[1][1] + v[2]*Si4647[1][2];
	vv[2]=v[3];
	vv[3]=-LOWERLEG + v[1]*Si4647[3][1] + v[2]*Si4647[3][2];

	v[1]=vv[1]*Si4546[1][1] + vv[2]*Si4546[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si4546[3][1] + vv[2]*Si4546[3][2];

	vv[1]=YKNEE + v[1]*Si4445[1][1] + v[2]*Si4445[1][2];
	vv[2]=-v[3];
	vv[3]=-UPPERLEG + v[1]*Si4445[3][1] + v[2]*Si4445[3][2];

	v[1]=YHIP + vv[1]*Si4344[1][1] + vv[2]*Si4344[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si4344[3][1] + vv[2]*Si4344[3][2];

	vv[1]=v[1]*Si4243[1][1] + v[2]*Si4243[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si4243[3][1] + v[2]*Si4243[3][2];

	v[1]=-XHIP + vv[1]*Si042[1][1] + vv[2]*Si042[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si042[3][1] + vv[2]*Si042[3][2];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {ZTOE, -XTOE, -YTOE} */
	Xlink[35][1]=vv[1];
	Xlink[35][2]=vv[2];
	Xlink[35][3]=vv[3];

	v[1]=ZHEEL*Si4748[1][1] + XHEEL*Si4748[1][2];
	v[2]=-YHEEL;
	v[3]=ZHEEL*Si4748[3][1] + XHEEL*Si4748[3][2];

	vv[1]=v[1]*Si4647[1][1] + v[2]*Si4647[1][2];
	vv[2]=v[3];
	vv[3]=-LOWERLEG + v[1]*Si4647[3][1] + v[2]*Si4647[3][2];

	v[1]=vv[1]*Si4546[1][1] + vv[2]*Si4546[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si4546[3][1] + vv[2]*Si4546[3][2];

	vv[1]=YKNEE + v[1]*Si4445[1][1] + v[2]*Si4445[1][2];
	vv[2]=-v[3];
	vv[3]=-UPPERLEG + v[1]*Si4445[3][1] + v[2]*Si4445[3][2];

	v[1]=YHIP + vv[1]*Si4344[1][1] + vv[2]*Si4344[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si4344[3][1] + vv[2]*Si4344[3][2];

	vv[1]=v[1]*Si4243[1][1] + v[2]*Si4243[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si4243[3][1] + v[2]*Si4243[3][2];

	v[1]=-XHIP + vv[1]*Si042[1][1] + vv[2]*Si042[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si042[3][1] + vv[2]*Si042[3][2];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {ZHEEL, XHEEL, YHEEL} */
	Xlink[36][1]=vv[1];
	Xlink[36][2]=vv[2];
	Xlink[36][3]=vv[3];

	v[1]=ZHEEL*Si4748[1][1] - XHEEL*Si4748[1][2];
	v[2]=-YHEEL;
	v[3]=ZHEEL*Si4748[3][1] - XHEEL*Si4748[3][2];

	vv[1]=v[1]*Si4647[1][1] + v[2]*Si4647[1][2];
	vv[2]=v[3];
	vv[3]=-LOWERLEG + v[1]*Si4647[3][1] + v[2]*Si4647[3][2];

	v[1]=vv[1]*Si4546[1][1] + vv[2]*Si4546[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si4546[3][1] + vv[2]*Si4546[3][2];

	vv[1]=YKNEE + v[1]*Si4445[1][1] + v[2]*Si4445[1][2];
	vv[2]=-v[3];
	vv[3]=-UPPERLEG + v[1]*Si4445[3][1] + v[2]*Si4445[3][2];

	v[1]=YHIP + vv[1]*Si4344[1][1] + vv[2]*Si4344[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si4344[3][1] + vv[2]*Si4344[3][2];

	vv[1]=v[1]*Si4243[1][1] + v[2]*Si4243[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si4243[3][1] + v[2]*Si4243[3][2];

	v[1]=-XHIP + vv[1]*Si042[1][1] + vv[2]*Si042[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si042[3][1] + vv[2]*Si042[3][2];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {ZHEEL, -XHEEL, YHEEL} */
	Xlink[37][1]=vv[1];
	Xlink[37][2]=vv[2];
	Xlink[37][3]=vv[3];

	v[1]=eff_x[4][1]*Si4748[1][1] + eff_x[4][2]*Si4748[1][2];
	v[2]=-eff_x[4][3];
	v[3]=eff_x[4][1]*Si4748[3][1] + eff_x[4][2]*Si4748[3][2];

	vv[1]=v[1]*Si4647[1][1] + v[2]*Si4647[1][2];
	vv[2]=v[3];
	vv[3]=-LOWERLEG + v[1]*Si4647[3][1] + v[2]*Si4647[3][2];

	v[1]=vv[1]*Si4546[1][1] + vv[2]*Si4546[1][2];
	v[2]=vv[3];
	v[3]=vv[1]*Si4546[3][1] + vv[2]*Si4546[3][2];

	vv[1]=YKNEE + v[1]*Si4445[1][1] + v[2]*Si4445[1][2];
	vv[2]=-v[3];
	vv[3]=-UPPERLEG + v[1]*Si4445[3][1] + v[2]*Si4445[3][2];

	v[1]=YHIP + vv[1]*Si4344[1][1] + vv[2]*Si4344[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si4344[3][1] + vv[2]*Si4344[3][2];

	vv[1]=v[1]*Si4243[1][1] + v[2]*Si4243[1][2];
	vv[2]=v[3];
	vv[3]=v[1]*Si4243[3][1] + v[2]*Si4243[3][2];

	v[1]=-XHIP + vv[1]*Si042[1][1] + vv[2]*Si042[1][2];
	v[2]=-vv[3];
	v[3]=vv[1]*Si042[3][1] + vv[2]*Si042[3][2];

	vv[1]=basec[1] + v[1]*Si00[1][1] + v[2]*Si00[1][2] + v[3]*Si00[1][3];
	vv[2]=basec[2] + v[1]*Si00[2][1] + v[2]*Si00[2][2] + v[3]*Si00[2][3];
	vv[3]=basec[3] + v[1]*Si00[3][1] + v[2]*Si00[3][2] + v[3]*Si00[3][3];

	/* link: {eff$4$$x[[1]], eff$4$$x[[2]], eff$4$$x[[3]]} */
	Xlink[38][1]=vv[1];
	Xlink[38][2]=vv[2];
	Xlink[38][3]=vv[3];

}

static void jacobian(const double link[NLINK+1][4],
		const double origin[NDOF+1][4],
		const double axis[NDOF+1][4],
		double jac[2*NCART][NDOF]) {

	static const int PALM = 6;
	static double c[2*NCART];
	for (int j = 1; j <= NDOF; ++j) {
		c[0] = axis[j][2] * (link[PALM][3] - origin[j][3]) - axis[j][3] * (link[PALM][2]-origin[j][2]);
		c[1] = axis[j][3] * (link[PALM][1] - origin[j][1]) - axis[j][1] * (link[PALM][3]-origin[j][3]);
		c[2] = axis[j][1] * (link[PALM][2] - origin[j][2]) - axis[j][2] * (link[PALM][1]-origin[j][1]);
		c[3] = axis[j][1];
		c[4] = axis[j][2];
		c[5] = axis[j][3];
		//rev_geo_jac_col(lp[PALM], jop[j], jap[j], c);
		for (int i = 0; i < 2*NCART; i++)
			jac[i][j-1] = c[i];
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

	mat limits;
	limits.load("config/joint_limits.cfg");
	lb = limits.col(0);
	ub = limits.col(1);
}

/**
 * @brief Reads the default joint states (initial posture) from file.
 *
 *
 */
static void read_default_state(vec & q_default) {

	mat limits;
	limits.load("config/joint_limits.cfg");
	q_default = limits.col(2);
}
