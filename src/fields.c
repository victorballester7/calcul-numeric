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

int exemple2(int n, double t, double x[], double f[], void *prm) {
#define N 2
  args_2 *prm2 = prm;
  double r2 = x[0] * x[0] + x[1] * x[1];
  f[0] = prm2->alpha * x[0] * (1 - r2) - x[1];
  f[1] = prm2->alpha * x[1] * (1 - r2) + x[0];

  if (n > N) {
#define A(i, j) x[N + (i) + (j)*N]
#define F(i, j) f[N + (i) + (j)*N]
    double f11, f12, f21, f22;  // fij = dfi/dxj
    f11 = prm2->alpha * (1 - r2) - 2 * prm2->alpha * x[0] * x[0];
    f12 = -2 * prm2->alpha * x[1] * x[0] - 1;
    f21 = -2 * prm2->alpha * x[1] * x[0] + 1;
    f22 = prm2->alpha * (1 - r2) - 2 * prm2->alpha * x[1] * x[1];
    F(0, 0) = f11 * A(0, 0) + f12 * A(1, 0);
    F(0, 1) = f11 * A(0, 1) + f12 * A(1, 1);
    F(1, 0) = f21 * A(0, 0) + f22 * A(1, 0);
    F(1, 1) = f21 * A(0, 1) + f22 * A(1, 1);
#undef A
#undef F
  }
#undef N
  return 0;
}

// already done
// int rtbp(int n, double t, double x[], double f[], void *prm) {
//   double X = x[0], Y = x[1], Z = x[2];
//   double u = x[3], v = x[4], w = x[5];
//   args_rtbps *prm2 = prm;
//   double mu = prm2->mu;
//   double rho1 = sqrt((X - mu) * (X - mu) + Y * Y + Z * Z);
//   double rho2 = sqrt((X - mu + 1) * (X - mu) + Y * Y + Z * Z);

//   f[0] = u;
//   f[1] = v;
//   f[2] = w;
//   f[3] = 2 * v + X - (X - mu) * (1 - mu) / (rho1 * rho1 * rho1) - (X - mu + 1) * mu / (rho2 * rho2 * rho2);
//   f[4] = -2 * u + Y - Y * (1 - mu) / (rho1 * rho1 * rho1) - Y * mu / (rho2 * rho2 * rho2);
//   f[5] = -Z * (1 - mu) / (rho1 * rho1 * rho1) - Z * mu / (rho2 * rho2 * rho2);

//   return 0;
// }
