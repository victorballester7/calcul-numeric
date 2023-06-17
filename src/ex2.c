#include <math.h>
#include <stdio.h>

#include "../include/fields.h"
#include "../include/flow.h"
#include "../include/rk78.h"

int main(int argc, char const *argv[]) {
  double hmin = 0.01, hmax = 0.05, tol = 1.e-6;
  int numMax = 20, n;
  double t, h, T = 0.5;
  args_2 prm;
  double x0 = 1, y0 = 1;

  prm.alpha = 0.4;
  t = 0;
  h = 0.01;
  n = 6;
  double x[6] = {x0, y0, 1, 0, 0, 1};
  flow(&t, x, &h, T, hmin, hmax, tol, numMax, n, exemple2, &prm);
  double y[4] = {x[2], x[3], x[4], x[5]};
  printf("Derivades del flux integrant les equacions variacionals en el punt (x, y) = (%lf, %lf):\n%lf %lf %lf %lf\n\n", x[0], x[1], y[0], y[1], y[2], y[3]);

  n = 2;
  double delta = 0.0001;
  double df11, df12, df21, df22;

  x[0] = x0 + delta;
  x[1] = y0;
  t = 0;
  h = 0.01;
  flow(&t, x, &h, T, hmin, hmax, tol, numMax, n, exemple2, &prm);
  df11 = x[0];
  df21 = x[1];

  x[0] = x0 - delta;
  x[1] = y0;
  t = 0;
  h = 0.01;
  flow(&t, x, &h, T, hmin, hmax, tol, numMax, n, exemple2, &prm);
  df11 = (df11 - x[0]) / (2 * delta);  // (fx(x0 + delta, y0) - fx(x0 - delta, y0)) / (2 * delta)
  df21 = (df21 - x[1]) / (2 * delta);  // (fy(x0 + delta, y0) - fy(x0 - delta, y0)) / (2 * delta)

  x[0] = x0;
  x[1] = y0 + delta;
  t = 0;
  h = 0.01;
  flow(&t, x, &h, T, hmin, hmax, tol, numMax, n, exemple2, &prm);
  df12 = x[0];
  df22 = x[1];

  x[0] = x0;
  x[1] = y0 - delta;
  t = 0;
  h = 0.01;
  flow(&t, x, &h, T, hmin, hmax, tol, numMax, n, exemple2, &prm);
  df12 = (df12 - x[0]) / (2 * delta);  // (fx(x0, y0 + delta) - fx(x0, y0 - delta)) / (2 * delta)
  df22 = (df22 - x[1]) / (2 * delta);  // (fy(x0, y0 + delta) - fy(x0, y0 - delta)) / (2 * delta)

  printf("Derivades del flux integrant per diferències finites el flux al mateix punt:\n%lf %lf %lf %lf\n\n", df11, df12, df21, df22);

  printf("Errors:\n%g %g %g %g\n", fabs(df11 - y[0]), fabs(df12 - y[1]), fabs(df21 - y[2]), fabs(df22 - y[3]));
  return 0;
}
