int flow(double *t, double x[], double *h, double T, double hmin, double hmax, double tol, int np, int n, int (*field)(int n, double t, double x[], double f[], void *prm), void *prm);