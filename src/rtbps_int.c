#include <stdio.h>

#include "../include/fields.h"
#include "../include/flow.h"
#include "../include/rk78.h"
#include "../include/rtbps.h"

int main(int argc, char const *argv[]) {
  double hmin = 0.01, hmax = 0.05, tol = 0.000001;
  int np, n = 6;
  double x[n], t, h, T = 6;
  args_rtbps prm;

  // check correct execution
  if (argc != 2 || sscanf(argv[1], "%lf", &prm.mu) != 1) {
    printf("Execute as ./bin/rtbps_int mu\n");
    return -1;
  }

  // file handling
  FILE *input, *output;
  input = fopen("data/initial_conditions_rtbps.txt", "r");
  output = fopen("data/output_rtbps.txt", "w");
  if (input == NULL || output == NULL) {
    printf("Error opening the files.\n");
    return -1;
  }

  while (fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %d", &x[0], &x[1], &x[2], &x[3], &x[4], &x[5], &h, &np) != EOF) {
    printf("numero np=%d\n", np);
    t = 0;
    fprintf(output, "\n\"Initial conditions: (x, y, z, vx, vy, vz, \u03BC) = (%g, %g, %g, %g, %g, %g, %g)\"\n", x[0], x[1], x[2], x[3], x[4], x[5], prm.mu);
    fprintf(output, "%lf %lf %lf %lf\n", t, x[0], x[1], x[2]);
    for (int i = 0; i < np; i++) {
      flow(&t, x, &h, T * 1. / np, hmin, hmax, tol, np, n, rtbps, &prm);
      fprintf(output, "%lf %lf %lf %lf\n", t, x[0], x[1], x[2]);
    }
    fprintf(output, "\n");
  }

  fclose(input);
  fclose(output);
  return 0;
}
