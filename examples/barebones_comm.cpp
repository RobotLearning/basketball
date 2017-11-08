/*
 * barebones_comm.cpp
 *
 *  Created on: Sep 12, 2017
 *      Author: okan.koc
 */

#include <stdio.h>
#include <stdlib.h>
#include "sys/time.h"
#include <math.h>
#include <nlopt.h>
#include <iostream>
#include <armadillo>
#include "barebones_interface.h"

using namespace std;
using namespace arma;

#define DIM 4 // dimension of the problem
#define CONSTR_DIM 2 // number of constraints

/**
 * Here we're testing the communication with SL
 * by linking the simplest function possible
 *
 */
void simplest_func_ever() {

	//cout << "Hello world!";
	mat A = zeros<mat>(4,5);
	mat B = zeros<mat>(4,5);

	cout << A*B.t() << endl;
}


