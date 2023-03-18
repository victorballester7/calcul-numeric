#include <math.h>
#include <stdio.h>

#include "../include/fields.h"
#include "../include/flow.h"
#include "../include/rk78.h"

int main(int argc, char const *argv[]) {
  double hmin = 0.01, hmax = 0.05, tol = 0.000001;
  int np = 100, n = 1;
  double x[n], t, t0, h = 0.01, T = 1;

  t = t0 = 1;
  x[0] = 1;
  if (!flow(&t, x, &h, T, hmin, hmax, tol, np, n, ivp2, NULL))
    printf("tolerance: %lf\n(t + T) ^ 2:\n\tReal:  %.10lf\n\tAprox: %.10lf\n\n", tol, (t0 + T) * (t0 + T), x[0]);
  else
    printf("Error in the numerical computation.\n");

  n = 2;
  t = t0 = 0;
  x[0] = 1;
  x[1] = 1;
  if (!flow(&t, x, &h, T, hmin, hmax, tol, np, n, harmonicOscillator, NULL))
    printf("tolerance: % lf\ncos(t + T) + sin(t + T):\n\tReal:  %.10lf\n\tAprox: %.10lf\n", tol, cos(t0 + T) + sin(t0 + T), x[0]);
  else
    printf("Error in the numerical computation.\n");
  return 0;
}
