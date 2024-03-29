#include <stdio.h>
#include <stdlib.h>

#include "../include/cmani.h"
#include "../include/fields.h"
#include "../include/rk78.h"

int main(int argc, char *argv[]) {
  double tolnwt;
  int maxitnwt;
  double h = 0.001, hmin = 0.001, hmax = 0.05, tolf = 1.e-12, maxNumSteps = 100;
  int m = 3, n = 2 * m;
  args_rtbps prm;

  printf("tolf = %e\n", tolf);

  // check correct execution
  if (argc != 4 || sscanf(argv[1], "%lf", &prm.mu) != 1 || sscanf(argv[2], "%lf", &tolnwt) != 1 || sscanf(argv[3], "%d", &maxitnwt) != 1) {
    printf("Execute as ./bin/rtbps_int mu tolnwt maxitnwt\n");
    return -1;
  }

  // file handling
  FILE *input, *output;
  input = fopen("data/input_cmani_rtbp.txt", "r");
  output = fopen("data/output_cmani_rtbp.txt", "w");
  if (input == NULL || output == NULL) {
    printf("Error opening the files.\n");
    return -1;
  }
  double dt;
  double *x0 = malloc(n * sizeof(double));
  double *xf = malloc(n * sizeof(double));
  double *dv = calloc(n, sizeof(double));  // calloc initializes to 0

  while (fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &dt, &x0[0], &x0[1], &x0[2], &x0[3], &x0[4], &x0[5], &xf[0], &xf[1], &xf[2], &xf[3], &xf[4], &xf[5]) != EOF) {
    for (int i = 0; i < n; i++) dv[i] = 0;
    if (cmani(m, x0, xf, dt, dv, tolnwt, maxitnwt, h, hmin, hmax, tolf, maxNumSteps, rtbps, &prm)) {
      printf("Error in cmani\n");
      return 1;
    }
    // print solution
    fprintf(output, "%.16g %.16g %.16g %.16g %.16g %.16g\n", dv[0], dv[1], dv[2], dv[3], dv[4], dv[5]);
  }

  fclose(input);
  fclose(output);
  return 0;
}
