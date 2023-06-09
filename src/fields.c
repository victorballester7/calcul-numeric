#include "../include/fields.h"

#include <math.h>
#include <stdio.h>

int pendulum(int n, double t, double x[], double f[], void *prm) {
#define N 2
  double r = x[0], v = x[1];
  args_pendulum *prm2 = prm;

  f[0] = v;
  f[1] = -sin(r) - prm2->mu * v + sin(prm2->w * t);

  if (n > N) {  // variational equation
#define A(i, j) x[N + (i)*N + (j)]
#define F(i, j) f[N + (i)*N + (j)]
    double f11, f12, f21, f22;  // fij = dfi/dxj
    f11 = 0;
    f12 = 1;
    f21 = -cos(r);
    f22 = -prm2->mu;
    F(0, 0) = f11 * A(0, 0) + f12 * A(1, 0);
    F(0, 1) = f11 * A(0, 1) + f12 * A(1, 1);
    F(1, 0) = f21 * A(0, 0) + f22 * A(1, 0);
    F(1, 1) = f21 * A(0, 1) + f22 * A(1, 1);
  }
#undef A
#undef F
#undef N
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
#define A(i, j) x[N + (i)*N + (j)]
#define F(i, j) f[N + (i)*N + (j)]
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

#define N 6
#define MU mu
#define M1 (1 - MU)
#define M2 MU
#define U x[3]
#define V x[4]
#define W x[5]

int rtbps(int n, double t, double *x, double *f, void *prm) {
  double mu, x_mu, x_mu_1, y, z, y2, z2, temp1, temp2;
  double rho1_2, rho2_2, rho1, rho2, rho1_3, rho2_3, pot3, potx3;
  mu = *((double *)prm);
  x_mu = x[0] - MU;
  x_mu_1 = x_mu + 1;
  y = x[1];
  y2 = y * y;
  z = x[2];
  z2 = z * z;
  temp1 = y2 + z2;
  rho1_2 = x_mu * x_mu + temp1;
  rho2_2 = x_mu_1 * x_mu_1 + temp1;
  rho1 = sqrt(rho1_2);
  rho2 = sqrt(rho2_2);
  rho1_3 = rho1_2 * rho1;
  rho2_3 = rho2_2 * rho2;
  temp1 = -M1 / rho1_3;
  temp2 = -M2 / rho2_3;
  pot3 = temp1 + temp2;
  potx3 = temp1 * x_mu + temp2 * x_mu_1;
  // Equations of the system
  f[0] = U;
  f[1] = V;
  f[2] = W;
  f[3] = 2 * V + x[0] + potx3;
  f[4] = -2 * U + y * (1 + pot3);
  f[5] = z * pot3;
  // Variational equations
  if (n > N) {
#define A0(i, j) x[N + (i)*N + (j)]  // initial matrix for the variational equations
#define Af(i, j) f[N + (i)*N + (j)]  // final matrix for the variational equations = Df * A
    double sum_pow5, sum_x_pow5, Omega_xx, Omega_xy, Omega_xz, Omega_yy, Omega_yz, Omega_zz;
    temp1 = M1 / (rho1_3 * rho1_2);
    temp2 = M2 / (rho2_3 * rho2_2);
    sum_pow5 = temp1 + temp2;
    temp1 *= x_mu;
    temp2 *= x_mu_1;
    sum_x_pow5 = temp1 + temp2;
    Omega_xx = 1 + pot3 + 3 * (temp1 * x_mu + temp2 * x_mu_1);
    Omega_xy = 3 * y * sum_x_pow5;
    Omega_xz = 3 * z * sum_x_pow5;
    Omega_yy = 1 + pot3 + 3 * y2 * sum_pow5;
    Omega_yz = 3 * y * z * sum_pow5;
    Omega_zz = pot3 + 3 * z2 * sum_pow5;
    for (int j = 0; j < N; j++) {  // we fill the matrix Af by columns
      Af(0, j) = A0(3, j);
      Af(1, j) = A0(4, j);
      Af(2, j) = A0(5, j);
      Af(3, j) = Omega_xx * A0(0, j) + Omega_xy * A0(1, j) + Omega_xz * A0(2, j) + 2 * A0(4, j);
      Af(4, j) = Omega_xy * A0(0, j) + Omega_yy * A0(1, j) + Omega_yz * A0(2, j) - 2 * A0(3, j);
      Af(5, j) = Omega_xz * A0(0, j) + Omega_yz * A0(1, j) + Omega_zz * A0(2, j);
    }
#undef A0
#undef Af
  }
  return 0;
}
#undef W
#undef V
#undef U
