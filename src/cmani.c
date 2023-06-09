#include "../include/cmani.h"

#include <stdio.h>

#include "../include/flow.h"
#include "../include/gauss.h"
int cmani_gdg(int m, double x0[], double xf[], double dt, double dv[], double g[], double dg[], double h, double hmin, double hmax, double tolfl, int maxNumSteps, int (*field)(int n, double t, double x[], double f[], void *prm), void *prm) {
  double t = 0;
  int n = 2 * m * (1 + 2 * m);
  double x[n];  // x is the vector with r0 and v0 + dv_0
  for (int i = 0; i < m; i++)
    x[i] = x0[i];
  for (int i = m; i < 2 * m; i++)
    x[i] = x0[i] + dv[i - m];
  for (int i = 0; i < 2 * m; i++) {
    for (int j = 0; j < 2 * m; j++)
      x[2 * m + i * 2 * m + j] = (i == j) ? 1 : 0;  // identity matrix (initial condition for the variational equations)
  }
  if (flow(&t, x, &h, dt / 2, hmin, hmax, tolfl, maxNumSteps, n, field, prm)) exit(1);

  double DvPhi05[(2 * m) * m];  // DvPhi(x_0.5)

  for (int i = 0; i < 2 * m; i++) {
    for (int j = 0; j < m; j++)
      DvPhi05[i * m + j] = x[2 * m + i * 2 * m + j + m];
  }

  for (int i = m; i < 2 * m; i++)
    x[i] += dv[i];
  // we reset the initial condition for the variational equations to the identity matrix
  for (int i = 0; i < 2 * m; i++) {
    for (int j = 0; j < 2 * m; j++)
      x[2 * m + i * 2 * m + j] = (i == j) ? 1 : 0;
  }
  if (flow(&t, x, &h, dt / 2, hmin, hmax, tolfl, maxNumSteps, n, field, prm)) exit(1);

  for (int i = 0; i < 2 * m; i++)
    g[i] = x[i] - xf[i];
  double DPhi15[(2 * m) * (2 * m)];  // DPhi(x_1.5)
  double DvPhi15[(2 * m) * m];       // DvPhi(x_1.5)

  memcpy(DPhi15, x + 2 * m, 2 * m * 2 * m * sizeof(double));

  for (int i = 0; i < 2 * m; i++) {
    for (int j = 0; j < m; j++)
      DvPhi15[i * m + j] = x[2 * m + i * 2 * m + j + m];
  }

  for (int i = 0; i < 2 * m; i++) {
    for (int j = 0; j < m; j++) {
      // first matrix
      dg[i * 2 * m + j] = 0;
      for (int k = 0; k < 2 * m; k++)  // DPhi15 * DvPhi05
        dg[i * 2 * m + j] += DPhi15[i * 2 * m + k] * DvPhi05[k * m + j];
      // second matrix
      dg[i * 2 * m + j + m] = DvPhi15[i * m + j];
    }
  }
  return 0;
}

int cmani(int m, double x0[], double xf[], double dt, double dv[], double tol, int maxit, double h, double hmin, double hmax, double tolfl, int maxNumSteps, int (*field)(int n, double t, double x[], double f[], void *prm), void *prm) {
  int iter = 0;
  int n = 2 * m;
  double *g = (double *)malloc(n * sizeof(double));
  double *dg = (double *)malloc(n * n * sizeof(double));
  if (cmani_gdg(m, x0, xf, dt, dv, g, dg, h, hmin, hmax, tolfl, maxNumSteps, field, prm)) exit(1);

  double error = tol + 1;

  if (cmani_gdg(m, x0, xf, dt, dv, g, dg, h, hmin, hmax, tolfl, maxNumSteps, field, prm)) exit(1);
  while (error > tol && iter < maxit) {
    // Newton's method
    gauss(n, dg, g, tol);
    for (int i = 0; i < n; i++)
      dv[i] -= g[i];  // new iterate
    // we evaluate G(dv) = g and its Jacobian dG(dv) = dg
    if (cmani_gdg(m, x0, xf, dt, dv, g, dg, h, hmin, hmax, tolfl, maxNumSteps, field, prm)) exit(1);
    error = 0;
    for (int i = 0; i < n; i++)
      error += fabs(g[i]);  // error of g = G(dv)
    iter++;
  }
  free(g);
  free(dg);
  if (iter == maxit) return 1;
  return 0;
}
