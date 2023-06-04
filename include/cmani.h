#ifndef CMANI_H
#define CMANI_H

#include <stdlib.h>
#include <string.h>

// ---------------------------------------------
// cmani_gdg
// ---------------------------------------------
// Purpose:
//   Compute the value of the function g and its differential dg
//
// Arguments:
//   m: dimension of the position space (or velocity space)
//   x0: initial condition
//   xf: final condition where we want to arrive
//   dt: time step for the maneuver
//   dv: initial guess for the increment of velocity for the maneuver
//   g: contains the value of the function g at the end of the function
//   dg: contains the value of the differential dg at the end of the function
//   h: initial step size for the integration of the variational equations
//   hmin: minimum step size for the integration of the variational equations
//   hmax: maximum step size for the integration of the variational equations
//   tolfl: tolerance for the integration of the variational equations
//   maxNumSteps: maximum number of steps for the integration of the variational equations
//   field: function that computes the field
//   prm: parameters for the field function
//
// Return value:
//   0: success
//   1: error
// ---------------------------------------------
int cmani_gdg(int m, double x0[], double xf[], double dt, double dv[], double g[], double dg[], double h, double hmin, double hmax, double tolfl, int maxNumSteps, int (*field)(int n, double t, double x[], double f[], void *prm), void *prm);

// ---------------------------------------------
// cmani
// ---------------------------------------------
// Purpose:
//   Compute the value increment dv necessary for the maneuver
//
// Arguments:
//   m: dimension of the position space (or velocity space)
//   x0: initial condition
//   xf: final condition where we want to arrive
//   dt: time step for the maneuver
//   dv: initial guess for the increment of velocity for the maneuver. At the end of the function, it contains the value of the increment of velocity for the maneuver
//   tol: tolerance for the Newton method
//   maxit: maximum number of iterations for the Newton method
//   h: initial step size for the integration of the variational equations
//   hmin: minimum step size for the integration of the variational equations
//   hmax: maximum step size for the integration of the variational equations
//   tolfl: tolerance for the integration of the variational equations
//   maxNumSteps: maximum number of steps for the integration of the variational equations
//   field: function that computes the field
//   prm: parameters for the field function
//
// Return value:
//   0: success
//   1: error
// ---------------------------------------------
int cmani(int m, double x0[], double xf[], double dt, double dv[], double tol, int maxit, double h, double hmin, double hmax, double tolfl, int maxNumSteps, int (*field)(int n, double t, double x[], double f[], void *prm), void *prm);

#endif  // CMANI_H
