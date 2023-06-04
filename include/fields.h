#ifndef FIELDS_H
#define FIELDS_H

typedef struct args_pendulum {
  double w;
  double mu;
} args_pendulum;

typedef struct args_lorenz {
  double sigma;
  double rho;
  double beta;
} args_lorenz;

typedef struct args_rtbps {
  double mu;
} args_rtbps;

typedef struct args_2 {
  double alpha;
} args_2;

#define RTBPS_N 6

int pendulum(int n, double t, double x[], double f[], void *prm);
int ivp2(int n, double t, double x[], double f[], void *prm);
int harmonicOscillator(int n, double t, double x[], double f[], void *prm);
int lorenz(int n, double t, double x[], double f[], void *prm);
int rtbp(int n, double t, double x[], double f[], void *prm);
int exemple2(int n, double t, double x[], double f[], void *prm);
int rtbps(int n, double t, double *x, double *f, void *prm);

#endif  // FIELDS_H