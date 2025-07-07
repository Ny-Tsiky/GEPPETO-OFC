#define NRANSI
#include "nrutil.h"
#define EPS 1.0e-6

extern int nvar;
extern double x1,x2;

// JNI 16/Jun/2000 Added "extern" (already defined in SHOOTF.c)
extern int kmax,kount;
extern double *xp,**yp,dxsav;

void shoot(int n, double v[], double f[])
{

void derivs(double x, double y[], double dydx[]);

void load(double x1, double v[], double y[]);

void odeint(double ystart[], int nvar, double x1, double x2,


double eps, double h1, double hmin, int *nok, int *nbad,


void (*derivs)(double, double [], double []),


void (*rkqs)(double [], double [], int, double *, double, double,


double [], double *, double *, void (*)(double, double [], double [])));

void rkqs(double y[], double dydx[], int n, double *x,


double htry, double eps, double yscal[], double *hdid, double *hnext,


void (*derivs)(double, double [], double []));

void score(double xf, double y[], double f[]);

int nbad,nok;

double h1,hmin=0.0,*y;


y=dvector(1,nvar);

kmax=0;

h1=(x2-x1)/100.0;

load(x1,v,y);

odeint(y,nvar,x1,x2,EPS,h1,hmin,&nok,&nbad,derivs,rkqs);

score(x2,y,f);

free_dvector(y,1,nvar);
}
#undef EPS
#undef NRANSI
