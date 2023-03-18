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

int pendulum(int n, double t, double x[], double f[], void *prm);
int ivp2(int n, double t, double x[], double f[], void *prm);
int harmonicOscillator(int n, double t, double x[], double f[], void *prm);
int lorenz(int n, double t, double x[], double f[], void *prm);
int rtbp(int n, double t, double x[], double f[], void *prm);
