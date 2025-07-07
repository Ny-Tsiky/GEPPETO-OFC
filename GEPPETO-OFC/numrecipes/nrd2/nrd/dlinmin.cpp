#define NRANSI
#include "nrutil.h"
#define TOL 6.0e-8

// JNI 16/Jun/2000 Added "extern" (already defined in module: LINMIN.c)
extern int ncom;
extern double *pcom,*xicom,(*nrfunc)(double []);
void (*nrdfun)(double [], double []);

void dlinmin(double p[], double xi[], int n, double *fret, double (*func)(double []),

void (*dfunc)(double [], double []))
{

double dbrent(double ax, double bx, double cx,


double (*f)(double), double (*df)(double), double tol, double *xmin);

double f1dim(double x);

double df1dim(double x);

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,


double *fc, double (*func)(double));

int j;

double xx,xmin,fx,fb,fa,bx,ax;


ncom=n;

pcom=dvector(1,n);

xicom=dvector(1,n);

nrfunc=func;

nrdfun=dfunc;

for (j=1;j<=n;j++) {


pcom[j]=p[j];


xicom[j]=xi[j];

}

ax=0.0;

xx=1.0;

mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);

*fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,&xmin);

for (j=1;j<=n;j++) {


xi[j] *= xmin;


p[j] += xi[j];

}

free_dvector(xicom,1,n);

free_dvector(pcom,1,n);
}
#undef TOL
#undef NRANSI
