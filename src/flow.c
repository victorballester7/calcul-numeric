#include "../include/flow.h"

#include <math.h>
#include <stdio.h>

#include "../include/rk78.h"

int flow(double *t, double x[], double *h, double T, double hmin, double hmax, double tol, int maxNumSteps, int n, int (*field)(int n, double t, double x[], double f[], void *prm), void *prm) {
  double t0 = *t;
  int count = 0;
  while (*t < t0 + T && count < maxNumSteps) {
    if (*t + *h > t0 + T)
      *h = t0 + T - *t;
    if (rk78(t, x, h, hmin, hmax, tol, n, field, prm)) return 1;
    count++;
  }
  if (count == maxNumSteps)
    return 1;  // this means count = maxNumSteps
  else {
    *t = t0 + T;
    return 0;
  }
}
