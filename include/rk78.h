#ifndef RK78_H
#define RK78_H

// ----------------------------------------------
// rk78
// ----------------------------------------------
// Purpose:
//    Runge-Kutta-Fehlberg 7(8) integrator
//
// Parameters:
//    t: pointer to the current time (initially pointing to the time of the initial conditions and at the end of the function pointing to the time at t + h)
//    x: pointer to the current state vector (initially pointing to the initial conditions and at the end of the function pointing to the state vector at t + h)
//    h: pointer to the wished step size (initially the step at which we want to integrate the system and at the end of the function pointing to the real step size of the next step)
//    hmin: minimum step size (if the step size is < 0, then hmin < 0)
//    hmax: maximum step size (if the step size is < 0, then hmax < 0)
//    tol: tolerance, which corresponds to the (theoretical) error between the real solution and the next step of the solution
//    n: dimension of the field
//    field: pointer to the function that calculates the field
//    param: pointer to the parameters of the field function
//
// Return value:
//    0: success
//    otherwise: error
// ----------------------------------------------
int rk78(double *t, double x[], double *h, double hmin, double hmax, double tol, int n, int (*field)(int n, double t, double x[], double f[], void *prm), void *prm);

#endif  // RK78_H
