#ifndef FLOW_H
#define FLOW_H
// ----------------------------------------------
// flow
// ----------------------------------------------
// Purpose:
//    Integrates the field from t to t + T
//
// Parameters:
//    t: pointer to the current time (initially pointing to the time of the initial conditions and at the end of the function pointing to the time at t + T)
//    x: pointer to the current state vector (initially pointing to the initial conditions and at the end of the function pointing to the state vector at t + T)
//    h: pointer to the wished step size (initially the solution time of the next step will as as close as possible to *t + h, where *t is the inital time)
//    T: integration time
//    hmin: minimum step size
//    hmax: maximum step size
//    tol: tolerance, which corresponds to the (theoretical) error between the real solution and the next step of the solution
//    maxNumSteps: maximum number of steps
//    n: dimension of the field
//    field: pointer to the function that calculates the field
//    param: pointer to the parameters of the field function
//
// Return value:
//    0: success
//    otherwise: error
// ----------------------------------------------
int flow(double *t, double x[], double *h, double T, double hmin, double hmax, double tol, int maxNumSteps, int n, int (*field)(int n, double t, double x[], double f[], void *prm), void *prm);

#endif  // FLOW_H
