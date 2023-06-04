#include <math.h>
#include <stdio.h>

#include "../include/cmani.h"
#include "../include/fields.h"

#define M_PI_2 1.57079632679489661923

int main(int argc, char const *argv[]) {
  int m = 1;
  int n = 2 * m;
  double *x0 = (double *)calloc(n, sizeof(double));
  double *xf = (double *)calloc(n, sizeof(double));
  double *dv = (double *)calloc(n, sizeof(double));  // calloc initializes to 0
  x0[0] = 1.;
  xf[1] = -sqrt(2. * (1. - cos(1.)));

  double dt = M_PI_2, tolnwt = 1.e-12, h = 0.01, hmin = 0.01, hmax = 0.05, tolfl = 1.e-12;
  int maxitnwt = 100, maxNumSteps = 100;
  args_pendulum prm;
  prm.w = 0.;
  prm.mu = 0.;

  if (cmani(m, x0, xf, dt, dv, tolnwt, maxitnwt, h, hmin, hmax, tolfl, maxNumSteps, pendulum, &prm)) {
    printf("Error in cmani\n");
    return 1;
  }
  printf("dv = [%.16f %.16f]\n", dv[0], dv[1]);
  return 0;
}
