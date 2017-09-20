/**
 * @file optimpoly.cpp
 * @brief Nonlinear optimization in C using the NLOPT library
 * @author Okan
 * @date 30/05/2016
 *
 */

#include <armadillo>
#include <thread>
#include "constants.h"
#include "utils.h"
#include "stdlib.h"
#include "math.h"
#include "kinematics.h"
#include "optim.h"

// termination
static bool check_optim_result(const int res);

// optimization related methods
static double costfunc(unsigned n, const double *x, double *grad, void *my_func_data);
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *f_data);
static void first_order_hold(const optim_des* data, const double T, double ball_pos[NCART], double ball_vel[NCART]);

/**
 * @brief Initialize the NLOPT optimization procedure here for FP
 * @param qrest_ Fixed resting posture
 * @param lb_ Fixed joint lower limits
 * @param ub_ Fixed joint upper limits
 */
Optim::Optim(const vec & qrest_, const ivec & active_dofs_, const bool right) : qrest(qrest_), right_arm(right) {

	double tol_eq[EQ_CONSTR_DIM];
	double tol_ineq[INEQ_CONSTR_DIM];
	const_vec(EQ_CONSTR_DIM,1e-2,tol_eq);
	const_vec(INEQ_CONSTR_DIM,1e-3,tol_ineq);

	double SLACK = 0.02;
	double Tmax = 2.0;
	active_dofs = active_dofs_;
	set_bounds(active_dofs,SLACK,Tmax,lb,ub);

	opt = nlopt_create(NLOPT_LN_COBYLA, OPTIM_DIM);
	nlopt_set_xtol_rel(opt, 1e-2);
	nlopt_set_lower_bounds(opt, lb.memptr());
	nlopt_set_upper_bounds(opt, ub.memptr());
	nlopt_set_min_objective(opt, costfunc, this);
	nlopt_add_inequality_mconstraint(opt, INEQ_CONSTR_DIM, joint_limits_ineq_constr, this, tol_ineq);
	nlopt_add_equality_mconstraint(opt, EQ_CONSTR_DIM, kinematics_eq_constr, this, tol_eq);
}

/**
 * @brief Update the initial state of optimization to PLAYER's current joint states.
 * @param qact Initial joint states acquired from sensors
 */
void Optim::update_init_state(const joint & qact) {
	for (int i = 0; i < NDOF_ACTIVE; i++) {
		q0[i] = qact.q(i);
		q0dot[i] = qact.qd(i);
	}
}

void Optim::set_active_dofs(const ivec & act_dofs) {
	active_dofs = act_dofs;
}

/**
 * @brief Tells the player optimization thread is still BUSY.
 *
 * If the (detached) thread is still running then table tennis player does not
 * update/launch new trajectories.
 * @return running
 */
bool Optim::check_running() {
	return running;
}

/**
 * @brief If the optimization was successful notify the Player class
 *
 * If the optimization was successful, update is turned ON and the table tennis
 * player can launch/update the polynomial trajectories.
 * @return update
 */
bool Optim::check_update() {
	return update;
}

/**
 * @brief If the robot starts moving the optimization is notified via this function
 *
 * If the robot is moving, this means last optimization was feasible, hence
 * we can initialize the new optimization from the last optimized parameter values.
 * @param flag_move
 */
void Optim::set_moving(bool flag_move) {
	moving = flag_move;
}

/**
 * @brief Detach the optimization thread.
 *
 * If the optimization is performed on SL or REAL_ROBOT then the optimization
 * thread should be detached.
 * @param flag_detach
 */
void Optim::set_detach(bool flag_detach) {
	detach = flag_detach;
}

/**
 * Set the final time for the returning trajectory
 * @param ret_time
 */
void Optim::set_return_time(const double & ret_time) {
	time2return = ret_time;
}

/**
 * @brief Print verbose optimization output (detailed optimization results are printed)
 * @param flag_verbose
 */
void Optim::set_verbose(bool flag_verbose) {
	verbose = flag_verbose;
}

/**
 * @brief If optimization succeeded, update polynomial parameters p
 *
 * If the optimizers finished running and terminated successfully,
 * then generates the striking and returning polynomial parameters
 * from qf, qfdot and T, and given the actual joint states qact, qactdot
 *
 */
bool Optim::get_params(const joint & qact, spline_params & p) {

	bool flag = false;
	if (update && !running) {
		vec qf_ = zeros<vec>(NDOF_ACTIVE);
		vec qfdot_ = zeros<vec>(NDOF_ACTIVE);
		for (int i = 0; i < NDOF_ACTIVE; i++) {
			qf_(i) = qf[i];
			qfdot_(i) = qfdot[i];
		}
		vec qnow = qact.q;
		vec qdnow = qact.qd;
		p.a.col(0) = 2.0 * (qnow - qf_) / pow(T,3) + (qfdot_ + qdnow) / pow(T,2);
		p.a.col(1) = 3.0 * (qf_ - qnow) / pow(T,2) - (qfdot_ + 2.0*qdnow) / T;
		p.a.col(2) = qdnow;
		p.a.col(3) = qnow;
		//cout << "A = \n" << p.a << endl;
		p.b.col(0) = 2.0 * (qf_ - qrest) / pow(time2return,3) + (qfdot_) / pow(time2return,2);
		p.b.col(1) = 3.0 * (qrest - qf_) / pow(time2return,2) - (2.0*qfdot_) / time2return;
		p.b.col(2) = qfdot_;
		p.b.col(3) = qf_;
		p.time2hit = T;
		//cout << "B = \n" << p.b << endl;
		flag = true;
		update = false;
	}
	return flag;
}

/**
 * @brief Update the rest state from outside
 *
 * If there is an additional optimization somewhere else
 * that optimizes for the resting posture, notify the optim classes
 * @param q_rest_new
 */
void Optim::update_rest_state(const vec & q_rest_new) {
	qrest = q_rest_new;
}

/**
 * @brief Set desired optimization parameters before running optim.
 *
 * @param params_ Desired optimization parameters are racket and/or ball values
 * predicted or computed by player class.
 */
void Optim::set_des_params(optim_des *params_) {
	param_des = params_;
}


/**
 * @brief Runs the optimization.
 *
 * Detaches the optimization if detach is set to TRUE. The
 * optimization method is shared by all Optim class descendants (VHP,FP,DP).
 */
void Optim::run() {

	std::thread t = std::thread(&Optim::optim,this);
	if (detach) {
		t.detach();
	}
	else {
		t.join();
	}
}

/**
 * @brief NLOPT optimization happens here.
 */
void Optim::optim() {

	update = false;
	running = true;
	double x[OPTIM_DIM];

	if (moving) {
		init_last_soln(x);
	}
	else {
		init_rest_soln(x);
	}

	double init_time = get_time();
	double past_time = 0.0;
	double minf; // the minimum objective value, upon return //
	int res; // error code

	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
		past_time = (get_time() - init_time)/1e3;
		if (verbose) {
			printf("NLOPT failed with exit code %d!\n", res);
		    printf("NLOPT took %f ms\n", past_time);
		}
	}
	else {
		past_time = (get_time() - init_time)/1e3;
		if (verbose) {
			printf("NLOPT success with exit code %d!\n", res);
			printf("NLOPT took %f ms\n", past_time);
			printf("Found minimum at f = %0.10g\n", minf);
		}
	    if (test_soln(x) < 1e-2)
	    	finalize_soln(x,past_time);
	}
	if (verbose)
		check_optim_result(res);
	running = false;
}

/**
 * @brief Initialize the optimization parameters the last optimized solution values.
 * @param x Optim params
 */
void Optim::init_last_soln(double x[]) const {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF_ACTIVE; i++) {
		x[i] = qf[i];
		x[i+NDOF_ACTIVE] = qfdot[i];
	}
	x[2*NDOF_ACTIVE] = T;
	//cout << "Initialization from T = " << T << endl;

}

/**
 * @brief Initialize the optim params to fixed resting posture.
 *
 * Initializes the optim params to qf fixed to q_rest, zero velocites,
 * and 0.5 hitting time.
 * @param x Optim params
 */
void Optim::init_rest_soln(double x[]) const {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF_ACTIVE; i++) {
		x[i] = qrest(i);
		x[i+NDOF_ACTIVE] = 0.0;
	}
	x[2*NDOF_ACTIVE] = 0.5;
}

/**
 * @brief Finalize solution if more than 50 ms is available for hitting.
 * @param x Optim params
 * @param time_elapsed Time elapsed during optimization
 */
void Optim::finalize_soln(const double x[], double time_elapsed) {

	if (x[2*NDOF_ACTIVE] > fmax(time_elapsed/1e3,0.05)) {
		// initialize first dof entries to q0
		for (int i = 0; i < NDOF_ACTIVE; i++) {
			qf[i] = x[i];
			qfdot[i] = x[i+NDOF_ACTIVE];
		}
		T = x[2*NDOF_ACTIVE];
		if (detach)
			T -= (time_elapsed/1e3);
		update = true;
	}
}

/**
 * @brief Test solution with hard kinematics constraints
 *
 * If constraints are violated then do not update/init. trajectories!
 *
 * @param x Optim params
 * @return Maximum value of constraint violations.
 */
double Optim::test_soln(const double x[]) const {

	// give info on constraint violation
	double *grad = 0;
	static double max_acc_violation; // at hitting time
	static double kin_violation[EQ_CONSTR_DIM];
	static double lim_violation[INEQ_CONSTR_DIM]; // joint limit violations on strike and return
	kinematics_eq_constr(EQ_CONSTR_DIM, kin_violation,
			             OPTIM_DIM, x, grad, (void*)this);
	joint_limits_ineq_constr(INEQ_CONSTR_DIM, lim_violation,
			                 OPTIM_DIM, x, grad, (void*)this);
	double cost = costfunc(OPTIM_DIM, x, grad, (void*)this);

	if (verbose) {
		// give info on solution vector
		print_optim_vec(x);
		printf("f = %.2f\n",cost);
		printf("Position constraint violation: [%.2f %.2f %.2f]\n",
				kin_violation[0],kin_violation[1],kin_violation[2]);
		for (int i = 0; i < INEQ_CONSTR_DIM; i++) {
			if (lim_violation[i] > 0.0)
				printf("Joint limit violated by %.2f on joint %d\n", lim_violation[i], i % NDOF_ACTIVE + 1);
		}
	}
	max_acc_violation = calc_max_acc_violation(x,q0,q0dot);

	return fmax(fmax(max_abs_array(kin_violation,EQ_CONSTR_DIM),
			    max_array(lim_violation,INEQ_CONSTR_DIM)),
			    max_acc_violation);
}

/*
 * Calculates the cost function for table tennis trajectory generation optimization
 * to find spline (3rd order strike+return) polynomials
 */
static double costfunc(unsigned n, const double *x, double *grad, void *my_func_params) {

	double a1[NDOF_ACTIVE];
	double a2[NDOF_ACTIVE];
	double T = x[2*NDOF_ACTIVE];

	if (grad) {
		static double h = 1e-6;
		static double val_plus, val_minus;
		static double xx[2*NDOF_ACTIVE+1];
		for (unsigned i = 0; i < n; i++)
			xx[i] = x[i];
		for (unsigned i = 0; i < n; i++) {
			xx[i] += h;
			val_plus = costfunc(n, xx, NULL, my_func_params);
			xx[i] -= 2*h;
			val_minus = costfunc(n, xx, NULL, my_func_params);
			grad[i] = (val_plus - val_minus) / (2*h);
			xx[i] += h;
		}
	}

	Optim *opt = (Optim*) my_func_params;
	double *q0 = opt->q0;
	double *q0dot = opt->q0dot;

	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);

	return T * (3*T*T*inner_prod(NDOF_ACTIVE,a1,a1) +
			3*T*inner_prod(NDOF_ACTIVE,a1,a2) + inner_prod(NDOF_ACTIVE,a2,a2));
}

/*
 * This is the constraint that makes sure we hit/touch the ball
 */
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *my_function_data) {

	static double des_pos[NCART];
	static double des_vel[NCART];
	static double pos_right[NCART];
	static double pos_left[NCART];
	static double qfdot[NDOF_ACTIVE];
	static double vel[NCART];
	static double qf[NDOF_ACTIVE];
	double T = x[2*NDOF_ACTIVE];

	Optim *opt = (Optim*) my_function_data;
	optim_des* data = opt->param_des;

	if (grad) {
		static double h = 1e-6;
		static double res_plus[EQ_CONSTR_DIM], res_minus[EQ_CONSTR_DIM];
		static double xx[2*NDOF_ACTIVE+1];
		for (unsigned i = 0; i < n; i++)
			xx[i] = x[i];
		for (unsigned i = 0; i < n; i++) {
			xx[i] += h;
			kinematics_eq_constr(m, res_plus, n, xx, NULL, my_function_data);
			xx[i] -= 2*h;
			kinematics_eq_constr(m, res_minus, n, xx, NULL, my_function_data);
			xx[i] += h;
			for (unsigned j = 0; j < m; j++)
				grad[j*n + i] = (res_plus[j] - res_minus[j]) / (2*h);
		}
	}

	// interpolate at time T to get the desired racket parameters
	first_order_hold(data,T,des_pos,des_vel);

	// extract state information from optimization variables
	for (int i = 0; i < NDOF_ACTIVE; i++) {
		qf[i] = x[i];
		qfdot[i] = x[i+NDOF_ACTIVE];
	}

	// compute the actual racket pos,vel and normal
	get_position(opt->active_dofs,qf,pos_left,pos_right);

	// deviations from the desired racket frame
	for (int i = 0; i < NCART; i++) {
		if (opt->right_arm)
			result[i] = pos_right[i] - des_pos[i];
		else
			result[i] = pos_left[i] - des_pos[i];
	}
}

/*
 * First order hold to interpolate linearly at time T
 * between ball pos,vel
 *
 * IF T is nan, ball variables are assigned to zero-element of
 * relevant ball entries
 *
 */
static void first_order_hold(const optim_des* data, const double T,
		                     double ball_pos[NCART], double ball_vel[NCART]) {

	double deltat = data->dt;
	if (std::isnan(T)) {
		printf("Warning: T value is nan!\n");

		for(int i = 0; i < NCART; i++) {
			ball_pos[i] = data->ball_pos(i,0);
			ball_vel[i] = data->ball_vel(i,0);
		}
	}
	else {
		int N = (int) (T/deltat);
		double Tdiff = T - N*deltat;
		int Nmax = data->Nmax;

		for (int i = 0; i < NCART; i++) {
			if (N < Nmax - 1) {
				ball_pos[i] = data->ball_pos(i,N) +
						(Tdiff/deltat) * (data->ball_pos(i,N+1) - data->ball_pos(i,N));
				ball_vel[i] = data->ball_vel(i,N) +
						(Tdiff/deltat) * (data->ball_vel(i,N+1) - data->ball_vel(i,N));
			}
			else {
				ball_pos[i] = data->ball_pos(i,Nmax-1);
				ball_vel[i] = data->ball_vel(i,Nmax-1);
			}
		}
	}
}

/*
 * This is the inequality constraint that makes sure we never exceed the
 * joint limits during the striking and returning motion
 *
 */
void joint_limits_ineq_constr(unsigned m, double *result,
		unsigned n, const double *x, double *grad, void *my_func_params) {

	static double a1[NDOF_ACTIVE];
	static double a2[NDOF_ACTIVE];
	static double a1ret[NDOF_ACTIVE]; // coefficients for the returning polynomials
	static double a2ret[NDOF_ACTIVE];
	static double joint_strike_max_cand[NDOF_ACTIVE];
	static double joint_strike_min_cand[NDOF_ACTIVE];
	static double joint_return_max_cand[NDOF_ACTIVE];
	static double joint_return_min_cand[NDOF_ACTIVE];
	static vec qdot_rest = zeros<vec>(NDOF_ACTIVE);

	Optim *opt = (Optim*) my_func_params;
	double *q0 = opt->q0;
	double *q0dot = opt->q0dot;
	double Tret = opt->time2return;

	if (grad) {
		static double h = 1e-6;
		static double res_plus[INEQ_CONSTR_DIM], res_minus[INEQ_CONSTR_DIM];
		static double xx[2*NDOF_ACTIVE+1];
		for (unsigned i = 0; i < n; i++)
			xx[i] = x[i];
		for (unsigned i = 0; i < n; i++) {
			xx[i] += h;
			joint_limits_ineq_constr(m, res_plus, n, xx, NULL, my_func_params);
			xx[i] -= 2*h;
			joint_limits_ineq_constr(m, res_minus, n, xx, NULL, my_func_params);
			xx[i] += h;
			for (unsigned j = 0; j < m; j++)
				grad[j*n + i] = (res_plus[j] - res_minus[j]) / (2*h);
		}
	}

	// calculate the polynomial coeffs which are used for checking joint limits
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);
	calc_return_poly_coeff(opt->qrest,qdot_rest,x,Tret,a1ret,a2ret);
	// calculate the candidate extrema both for strike and return
	calc_strike_extrema_cand(a1,a2,x[2*NDOF_ACTIVE],q0,q0dot,
			joint_strike_max_cand,joint_strike_min_cand);
	calc_return_extrema_cand(a1ret,a2ret,x,Tret,joint_return_max_cand,joint_return_min_cand);

	/* deviations from joint min and max */
	for (int i = 0; i < NDOF_ACTIVE; i++) {
		result[i] = joint_strike_max_cand[i] - opt->ub(i);
		result[i+NDOF_ACTIVE] = opt->lb(i) - joint_strike_min_cand[i];
		result[i+2*NDOF_ACTIVE] = joint_return_max_cand[i] - opt->ub(i);
		result[i+3*NDOF_ACTIVE] = opt->lb(i) - joint_return_min_cand[i];
		//printf("%f %f %f %f\n", result[i],result[i+DOF],result[i+2*DOF],result[i+3*DOF]);
	}

}

/*
 * Calculate the polynomial coefficients from the optimized variables qf,qfdot,T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
void calc_strike_poly_coeff(const double *q0, const double *q0dot, const double *x,
		                    double *a1, double *a2) {

	double T = x[2*NDOF_ACTIVE];

	for (int i = 0; i < NDOF_ACTIVE; i++) {
		a1[i] = (2/pow(T,3))*(q0[i]-x[i]) + (1/(T*T))*(q0dot[i] + x[i+NDOF_ACTIVE]);
		a2[i] = (3/(T*T))*(x[i]-q0[i]) - (1/T)*(x[i+NDOF_ACTIVE] + 2*q0dot[i]);
	}

	return;
}

/*
 * Calculate the returning polynomial coefficients from the optimized variables qf,qfdot
 * and time to return constant T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
void calc_return_poly_coeff(const vec & qrest, const vec & q0dot,
		                    const double *x, const double T,
		                    double *a1, double *a2) {

	for (int i = 0; i < NDOF_ACTIVE; i++) {
		a1[i] = (2/pow(T,3))*(x[i]-qrest(i)) + (1/(T*T))*(q0dot(i) + x[i+NDOF_ACTIVE]);
		a2[i] = (3/(T*T))*(qrest(i)-x[i]) - (1/T)*(2*x[i+NDOF_ACTIVE] + q0dot(i));
	}

}

/*
 * Calculate the extrema candidates for each joint (2*7 candidates in total)
 * For the striking polynomial
 * Clamp to [0,T]
 *
 */
void calc_strike_extrema_cand(const double *a1, const double *a2, const double T,
		                      const double *q0, const double *q0dot,
		                      double *joint_max_cand, double *joint_min_cand) {

	static double cand1, cand2;

	for (int i = 0; i < NDOF_ACTIVE; i++) {
		cand1 = fmin(T,fmax(0,(-a2[i] + sqrt(a2[i]*a2[i] - 3*a1[i]*q0dot[i]))/(3*a1[i])));
		cand2 =  fmin(T,fmax(0,(-a2[i] - sqrt(a2[i]*a2[i] - 3*a1[i]*q0dot[i]))/(3*a1[i])));
		cand1 = a1[i]*pow(cand1,3) + a2[i]*pow(cand1,2) + q0dot[i]*cand1 + q0[i];
		cand2 = a1[i]*pow(cand2,3) + a2[i]*pow(cand2,2) + q0dot[i]*cand2 + q0[i];
		joint_max_cand[i] = fmax(cand1,cand2);
		joint_min_cand[i] = fmin(cand1,cand2);
	}
}

/*
 * Calculate the extrema candidates for each joint (2*7 candidates in total)
 * For the return polynomial
 * Clamp to [0,TIME2RETURN]
 *
 */
void calc_return_extrema_cand(const double *a1, const double *a2,
		                      const double *x, const double Tret,
		                      double *joint_max_cand, double *joint_min_cand) {

	static double cand1, cand2;

	for (int i = 0; i < NDOF_ACTIVE; i++) {
		cand1 = fmin(Tret, fmax(0,(-a2[i] + sqrt(a2[i]*a2[i] - 3*a1[i]*x[i+NDOF_ACTIVE]))/(3*a1[i])));
		cand2 =  fmin(Tret, fmax(0,(-a2[i] - sqrt(a2[i]*a2[i] - 3*a1[i]*x[i+NDOF_ACTIVE]))/(3*a1[i])));
		cand1 = a1[i]*pow(cand1,3) + a2[i]*pow(cand1,2) + x[i+NDOF_ACTIVE]*cand1 + x[i];
		cand2 = a1[i]*pow(cand2,3) + a2[i]*pow(cand2,2) + x[i+NDOF_ACTIVE]*cand2 + x[i];
		joint_max_cand[i] = fmax(cand1,cand2);
		joint_min_cand[i] = fmin(cand1,cand2);
	}
}

/*
 * Since the accelerations of third order polynomials
 * are linear functions of time we check the values
 * at start of traj, t = 0 and end of traj, t = T_hit
 * which are given by 6*a1*T + a2 and a2 respectively
 *
 * Returns the max acc
 */
double calc_max_acc_violation(const double x[2*NDOF_ACTIVE+1],
		const double q0[NDOF_ACTIVE],
		const double q0dot[NDOF_ACTIVE]) {

	double T = x[2*NDOF_ACTIVE];
	double a1[NDOF_ACTIVE], a2[NDOF_ACTIVE];
	double acc_abs_max = 0.0;
	double acc_max_cand;

	calc_strike_poly_coeff(q0,q0dot,x,a1,a2); // get a1,a2 out

	for (int i = 0; i < NDOF_ACTIVE; i++) {
		acc_max_cand = fmax(fabs(6*a1[i]*T + a2[i]),fabs(a2[i]));
		//printf("qdd_max[%d] = %f\n", i, acc_max_cand);
		if (acc_max_cand > MAX_ACC && acc_max_cand > acc_abs_max) {
			acc_abs_max = acc_max_cand;
		}
	}
	return acc_abs_max;
}

/*
 * Give info about the optimization after termination
 *
 */
static bool check_optim_result(const int res) {

	bool flag = false;
	switch (res) {
	case NLOPT_SUCCESS:
		printf("Success!\n");
		flag = true;
		break;
	case NLOPT_STOPVAL_REACHED:
		printf("Optimization stopped because stopval (above) was reached.\n");
		flag = true;
		break;
	case NLOPT_FTOL_REACHED:
		printf("Optimization stopped because ftol_rel "
				"or ftol_abs (above) was reached.\n");
		flag = true;
		break;
	case NLOPT_XTOL_REACHED:
		flag = true;
		printf("Optimization stopped because xtol_rel or xtol_abs (above) was reached.\n");
		break;
	case NLOPT_MAXEVAL_REACHED:
		flag = true;
		printf("Optimization stopped because maxeval (above) was reached.\n");
		break;
	case NLOPT_MAXTIME_REACHED:
		flag = true;
		printf("Optimization stopped because maxtime (above) was reached.\n");
		break;
	case NLOPT_FAILURE:
		printf("Epic fail!\n");
		break;
	case NLOPT_INVALID_ARGS:
		printf("Invalid arguments (e.g. lower bounds are bigger than "
				"upper bounds, an unknown algorithm was specified, etcetera).\n");
		break;
	case NLOPT_OUT_OF_MEMORY:
		printf("Ran out of memory!\n");
		break;
	case NLOPT_ROUNDOFF_LIMITED:
		printf("Halted because roundoff errors limited progress."
			"(In this case, the optimization still typically returns a useful result.\n");
		break;
	case NLOPT_FORCED_STOP:
		printf("Halted because of a forced termination: "
				"the user called nlopt_force_stop(opt)"
				"on the optimization’s nlopt_opt object "
				"opt from the user’s objective function or constraints.\n");
		break;

	}
	return flag;
}

/*
 * Set upper and lower bounds on the optimization.
 * First loads the joint limits and then puts some slack
 */
void set_bounds(const ivec & active_dofs, const double SLACK, const double Tmax, vec & lb, vec & ub) {

	vec lb_full = zeros<vec>(NDOF);
	vec ub_full = zeros<vec>(NDOF);
	read_joint_limits(lb_full,ub_full);
	// lower bounds and upper bounds for qf are the joint limits
	for (int i = 0; i < NDOF_ACTIVE; i++) {
		ub(i) = ub_full(active_dofs(i)) - SLACK;
		lb(i) = lb_full(active_dofs(i)) + SLACK;
		ub(NDOF_ACTIVE + i) = MAX_VEL;
		lb(NDOF_ACTIVE + i) = -MAX_VEL;
	}
	// constraints on final time
	ub(2*NDOF_ACTIVE) = Tmax;
	lb(2*NDOF_ACTIVE) = 0.01;
}
