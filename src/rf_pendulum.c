#include <stdio.h>

#include "../include/fields.h"
#include "../include/rk78.h"

int main(int argc, char const *argv[]) {
  double hmin, hmax, tol;
  int np, n = 2;
  double x[n], t, h;
  args_pendulum prm;
  prm.mu = 0;
  prm.w = 0;

  // check correct execution
  if (argc != 4 || sscanf(argv[1], "%lf", &hmin) != 1 || sscanf(argv[2], "%lf", &hmax) != 1 || sscanf(argv[3], "%lf", &tol) != 1) {
    printf("Execute as ./bin/rf_pendol hmin hmax tol\n");
    return -1;
  }

  // file handling
  FILE *input, *output;
  input = fopen("data/initial_conditions_pendulum.txt", "r");
  output = fopen("data/output_pendulum.txt", "w");
  if (input == NULL || output == NULL) {
    printf("Error opening the files.\n");
    return -1;
  }

  while (fscanf(input, "%lf %lf %lf %d", &x[0], &x[1], &hmin, &np) != EOF) {
    t = 0;
    h = 0.001;
    fprintf(output, "\n\"Initial conditions: (x, v) = (%g, %g)\"\n", x[0], x[1]);
    fprintf(output, "%lf %lf %lf %lf\n", t, x[0], x[1], h);
    for (int i = 0; i < np; i++) {
      rk78(&t, x, &h, hmin, hmax, tol, n, pendulum, &prm);
      fprintf(output, "%lf %lf %lf %lf\n", t, x[0], x[1], h);
    }
    fprintf(output, "\n");
  }

  fclose(input);
  fclose(output);
  return 0;
}
