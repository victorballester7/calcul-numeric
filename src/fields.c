#include "../include/fields.h"

#include <math.h>
#include <stdio.h>

int pendulum(int n, double t, double x[], double f[], void *prm) {
  double r = x[0], v = x[1];
  args_pendulum *prm2 = prm;

  f[0] = v;
  f[1] = -sin(r) - prm2->mu * v + sin(prm2->w * t);

  return 0;
}

int ivp2(int n, double t, double x[], double f[], void *prm) {
  f[0] = 2 * x[0] / t;
  return 0;
}

int harmonicOscillator(int n, double t, double x[], double f[], void *prm) {
  f[0] = x[1];
  f[1] = -x[0];

  return 0;
}

int lorenz(int n, double t, double x[], double f[], void *prm) {
  double X = x[0], Y = x[1], Z = x[2];
  args_lorenz *prm2 = prm;

  f[0] = prm2->sigma * (Y - X);
  f[1] = -X * Z + prm2->rho * X - Y;
  f[2] = X * Y - prm2->beta * Z;

  return 0;
}

int rtbp(int n, double t, double x[], double f[], void *prm) {
  double X = x[0], Y = x[1], Z = x[2];
  double u = x[3], v = x[4], w = x[5];
  args_rtbps *prm2 = prm;
  double mu = prm2->mu;
  double rho1 = sqrt((X - mu) * (X - mu) + Y * Y + Z * Z);
  double rho2 = sqrt((X - mu + 1) * (X - mu) + Y * Y + Z * Z);

  f[0] = u;
  f[1] = v;
  f[2] = w;
  f[3] = 2 * v + X - (X - mu) * (1 - mu) / (rho1 * rho1 * rho1) - (X - mu + 1) * mu / (rho2 * rho2 * rho2);
  f[4] = -2 * u + Y - Y * (1 - mu) / (rho1 * rho1 * rho1) - Y * mu / (rho2 * rho2 * rho2);
  f[5] = -Z * (1 - mu) / (rho1 * rho1 * rho1) - Z * mu / (rho2 * rho2 * rho2);

  return 0;
}
