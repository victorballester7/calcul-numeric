#include <stdio.h>

#include "../include/fields.h"
#include "../include/flow.h"
#include "../include/rk78.h"

int main(int argc, char const *argv[]) {
  double hmin = 0.01, hmax = 0.05, tol = 0.000001;
  int np, n = 3;
  double x[n], t, h, T;
  args_lorenz prm;

  // check correct execution
  if (argc != 6 || sscanf(argv[1], "%lf", &prm.sigma) != 1 || sscanf(argv[2], "%lf", &prm.rho) != 1 || sscanf(argv[3], "%lf", &prm.beta) != 1 || sscanf(argv[4], "%lf", &T) != 1 || sscanf(argv[5], "%d", &np) != 1) {
    printf("Execute as ./bin/lorenz_int sigma rho beta tf nt\n");
    return -1;
  }

  // file handling
  FILE *input, *output;
  input = fopen("data/initial_conditions_lorenz.txt", "r");
  output = fopen("data/output_lorenz.txt", "w");
  if (input == NULL || output == NULL) {
    printf("Error opening the files.\n");
    return -1;
  }

  while (fscanf(input, "%lf %lf %lf", &x[0], &x[1], &x[2]) != EOF) {
    t = 0;
    h = 0.01;
    fprintf(output, "\n\"Initial conditions: (x, y, z, \u03C3, \u03C1, \u03B2) = (%g, %g, %g, %g, %g, %g)\"\n", x[0], x[1], x[2], prm.sigma, prm.rho, prm.beta);
    fprintf(output, "%lf %lf %lf %lf\n", t, x[0], x[1], x[2]);
    for (int i = 0; i < np; i++) {
      flow(&t, x, &h, T * 1. / np, hmin, hmax, tol, np, n, lorenz, &prm);
      fprintf(output, "%lf %lf %lf %lf\n", t, x[0], x[1], x[2]);
    }
    fprintf(output, "\n");
  }

  fclose(input);
  fclose(output);
  return 0;
}
