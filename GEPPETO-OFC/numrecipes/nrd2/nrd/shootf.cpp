#define NRANSI
#include "nrutil.h"
#define EPS 1.0e-6

//int nn2,nvar;
//double x1,x2,xf;
// PB 03/Oct/2003 These extern can't be found anywhere
extern int nn2,nvar;
extern double x1,x2,xf;

// JNI 16/Jun/2000 These are already defined in SHOOT.c & HYPGEO.c
// LIB complains about this.
int kmax,kount;
double *xp,**yp,dxsav;

void shootf(int n, double v[], double f[])
{

void derivs(double x, double y[], double dydx[]);

void load1(double x1, double v1[], double y[]);

void load2(double x2, double v2[], double y[]);

void odeint(double ystart[], int nvar, double x1, double x2,


double eps, double h1, double hmin, int *nok, int *nbad,


void (*derivs)(double, double [], double []),


void (*rkqs)(double [], double [], int, double *, double, double,


double [], double *, double *, void (*)(double, double [], double [])));

void rkqs(double y[], double dydx[], int n, double *x,


double htry, double eps, double yscal[], double *hdid, double *hnext,


void (*derivs)(double, double [], double []));

void score(double xf, double y[], double f[]);

int i,nbad,nok;

double h1,hmin=0.0,*f1,*f2,*y;


f1=dvector(1,nvar);

f2=dvector(1,nvar);

y=dvector(1,nvar);

kmax=0;

h1=(x2-x1)/100.0;

load1(x1,v,y);

odeint(y,nvar,x1,xf,EPS,h1,hmin,&nok,&nbad,derivs,rkqs);

score(xf,y,f1);

load2(x2,&v[nn2],y);

odeint(y,nvar,x2,xf,EPS,h1,hmin,&nok,&nbad,derivs,rkqs);

score(xf,y,f2);

for (i=1;i<=n;i++) f[i]=f1[i]-f2[i];

free_dvector(y,1,nvar);

free_dvector(f2,1,nvar);

free_dvector(f1,1,nvar);
}
#undef EPS
#undef NRANSI
