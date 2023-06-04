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
  // memcpy(x, x0, m * sizeof(double));
  for (int i = m; i < 2 * m; i++)
    x[i] = x0[i] + dv[i - m];
  for (int i = 0; i < 2 * m; i++) {
    for (int j = 0; j < 2 * m; j++)
      x[2 * m + i * 2 * m + j] = (i == j) ? 1 : 0;  // identity matrix (initial condition for the variational equations)
  }
  printf("abans-1r-flow\n");
  printf("t = %lf, h = %lf, dt = %lf, tolfl = %g, maxNumSteps = %d, n = %d\n", t, h, dt, tolfl, maxNumSteps, n);
  if (flow(&t, x, &h, dt / 2, hmin, hmax, tolfl, maxNumSteps, n, field, prm))
    return 1;
  for (int i = m; i < 2 * m; i++)
    x[i] += dv[i];

  double DvPhi05[(2 * m) * m];  // DvPhi(x_0.5)

  for (int i = 0; i < 2 * m; i++) {
    for (int j = 0; j < m; j++)
      DvPhi05[i * m + j] = x[2 * m + i * 2 * m + j + m];
  }

  // we reset the initial condition for the variational equations to the identity matrix
  for (int i = 0; i < 2 * m; i++) {
    for (int j = 0; j < 2 * m; j++)
      x[2 * m + i * 2 * m + j] = (i == j) ? 1 : 0;
  }
  printf("abans-2n-flow\n");
  if (flow(&t, x, &h, dt / 2, hmin, hmax, tolfl, maxNumSteps, n, field, prm))
    return 1;

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
  double *x = (double *)malloc(n * sizeof(double));
  double *dg = (double *)malloc(n * n * sizeof(double));
  printf("hola0\n");
  if (cmani_gdg(m, x0, xf, dt, dv, g, dg, h, hmin, hmax, tolfl, maxNumSteps, field, prm)) {
    free(g);
    free(dg);
    return 1;
  }
  printf("hola1\n");
  double error = tol + 1;
  while (error > tol && iter < maxit) {
    // Newton's method
    resgauss(n, dg, g, x);  // x = -dg^-1 * g
    memcpy(g, x, n * sizeof(double));
    printf("hola---\n");
    // gauss(n, dg, g, tol);
    for (int i = 0; i < n; i++)
      dv[i] -= g[i];  // new iterate
    // -----------------
    // we evaluate G(dv) = g and its Jacobian dG(dv) = dg
    if (cmani_gdg(m, x0, xf, dt, dv, g, dg, h, hmin, hmax, tolfl, maxNumSteps, field, prm)) {
      free(g);
      free(dg);
      return 1;
    }
    error = 0;
    printf("dv = [");
    for (int i = 0; i < n; i++) {
      printf("%g ", dv[i]);
    }
    printf("]\n");
    for (int i = 0; i < n; i++)
      error += fabs(g[i]);  // error of g = G(dv)
    printf("error = %g\n", error);
    iter++;
  }
  printf("hola2, iter = %d\n", iter);
  free(g);
  free(dg);
  if (iter == maxit) return 1;
  return 0;
}
