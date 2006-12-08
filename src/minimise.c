#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>
#include <math.h>

void            dist(double *, int, int, double *);
double          fnorm(double *, int);
void            grad(double *, double *, int, int, double *);
void            get_M(double *, double *, int, double, double, double, double, double *);
double          Energy(double *, double *, int, double, double, double, double);
void            XsubsY(double *, double *, int, double *);
void            aX(double, double *, int, double *);
void            equal(double *, int, double *);

const double    myInf = 1e36;

void
doubledist(double *x, int *n, int *d, double *ddist)
{
  int             i;
  double         *dis = (double *) malloc(*n * *n * sizeof(double));
  dist(x, *n, *d, dis);
  for (i = 0; i < *n * *n; i++)
    ddist[i] = 2 * dis[i];
}

void
boxcox(double *X1, double *D0, int *n, int *d,
       double *lam, double *mu, double *nu, double *c, double *energy,
       int *niter)
{
  double          minsize = 1e-4;
  double          tmp;
  double          stepsize = 0.1, s0 = 10, s1 = s1 - 1;
  int             i = 0, j, count = 0, count1 = 0, count2 = 0;
  double         *X0 = (double *) malloc(*n * *d * sizeof(double));
  double         *tmpX = (double *) malloc(*n * *d * sizeof(double));
  double         *D1 = (double *) malloc(*n * *n * sizeof(double));
  double         *M = (double *) malloc(*n * *n * sizeof(double));
  double         *Grad = (double *) malloc(*n * *d * sizeof(double));
  double         *normGrad = (double *) malloc(*n * *d * sizeof(double));

  dist(X1, *n, *d, D1);

  while (stepsize > minsize && count < *niter) {
    if (s1 > s0 && count > 1) {
      stepsize = 0.5 * stepsize;
      aX(stepsize, normGrad, *n * *d, tmpX);
      XsubsY(X0, tmpX, *n * *d, X1);
    } else {
      stepsize = 1.05 * stepsize;
      equal(X1, *n * *d, X0);
      get_M(D0, D1, *n, *lam, *mu, *nu, *c, M);
      grad(X0, M, *n, *d, Grad);
      aX(fnorm(X0, *n * *d) / fnorm(Grad, *n * *d), Grad, *n * *d, normGrad);
      aX(stepsize, normGrad, *n * *d, tmpX);
      XsubsY(X0, tmpX, *n * *d, X1);
    }
    count++;
    s0 = s1;
    dist(X1, *n, *d, D1);
    s1 = Energy(D0, D1, *n, *lam, *mu, *nu, *c);
//    if (count % 50 == 0) {
//      printf("Inter: %d\n", count);
//      printf("Energy: %f\n", s1);
//    }
  }
  *energy = s1;
}

/* ------------------ Functions ------------------- */


void
dist(double *x, int n, int d, double *dis)
{
  int             i, j, k;
  double          s, r;

  for (i = 0; i < n; i++)
    for (j = 0; j < i; j++) {
      s = 0.0;
      for (k = 0; k < d; k++) {
        r = x[k + i * d] - x[k + j * d];
        s += r * r;
      }
      dis[i + j * n] = dis[j + i * n] = sqrt(s);
    }
}

double
fnorm(double *x, int len)
{
  int             i;
  double          s = 0;

  for (i = 0; i < len; i++)
    s += pow(x[i], 2);
  return sqrt(s);
}

void
grad(double *X0, double *M, int n, int d, double *Grad)
{
  int             i, j, k;
  double          tmp;

  for (k = 0; k < d; k++)
    for (i = 0; i < n; i++) {
      tmp = 0.0;
      for (j = 0; j < n; j++)
        tmp += (X0[k + i * d] - X0[k + j * d]) * M[i + j * n];
      Grad[k + i * d] = tmp;
    }
}

void
get_M(double *D0, double *D1,
      int n, double lam, double mu, double nu, double c,
      double *M)
{
  int             i, j;
  double          s, d0, d1;
  /* Diagnal element is set to be 0 */
  for (i = 0; i < n; i++)
    M[i + i * n] = 0;

  for (i = 0; i < n; i++)
    for (j = 0; j < i; j++) {
      d0 = D0[i + j * n];
      d1 = D1[i + j * n];

      if (d0 < myInf)
        s = pow(d0, nu) * pow(d1, mu + 1 / lam - 2) - pow(d0, nu + 1 / lam) * pow(d1, mu - 2);
      else
        s = -c * pow(d1, mu - 2);

      M[i + j * n] = M[j + i * n] = s;
    }
}

double
Energy(double *D0, double *D1, int n, double lam, double mu, double nu, double c)
{
  int             i, j, count;
  double          d0, d1, bc1, bc2, s, en = 0;

  for (i = 0; i < n; i++)
    for (j = 0; j < i; j++) {
      d0 = D0[i + j * n];
      d1 = D1[i + j * n];

      if (mu + 1 / lam == 0)
        bc1 = log(d1);
      else
        bc1 = (pow(d1, mu + 1 / lam) - 1) / (mu + 1 / lam);

      if (mu == 0)
        bc2 = log(d1);
      else
        bc2 = (pow(d1, mu) - 1) / mu;

      if (d0 < myInf)
        s = pow(d0, nu) * bc1 - pow(d0, nu + 1 / lam) * bc2;
      else
        s = -c * bc2;

      en = en + 2 * s;
    }
  return en;
}

void
XsubsY(double *X, double *Y, int len, double *Z)
{
  int             i;
  for (i = 0; i < len; i++)
    Z[i] = X[i] - Y[i];
}

void
aX(double a, double *X, int len, double *Z)
{
  int             i;
  for (i = 0; i < len; i++)
    Z[i] = a * X[i];
}

void
equal(double *X, int len, double *Y)
{
  int             i = 0;
  for (i = 0; i < len; i++)
    Y[i] = X[i];
}
