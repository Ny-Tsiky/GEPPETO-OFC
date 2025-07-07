#include <math.h>
#define NRANSI
#include "nrutil.h"

void ks2d1s(double x1[], double y1[], unsigned long n1,

void (*quadvl)(double, double, double *, double *, double *, double *),

double *d1, double *prob)
{

void pearsn(double x[], double y[], unsigned long n, double *r, double *prob,


double *z);

double probks(double alam);

void quadct(double x, double y, double xx[], double yy[], unsigned long nn,


double *fa, double *fb, double *fc, double *fd);

unsigned long j;

double dum,dumm,fa,fb,fc,fd,ga,gb,gc,gd,r1,rr,sqen;


*d1=0.0;

for (j=1;j<=n1;j++) {


quadct(x1[j],y1[j],x1,y1,n1,&fa,&fb,&fc,&fd);


(*quadvl)(x1[j],y1[j],&ga,&gb,&gc,&gd);


*d1=DMAX(*d1,fabs(fa-ga));


*d1=DMAX(*d1,fabs(fb-gb));


*d1=DMAX(*d1,fabs(fc-gc));


*d1=DMAX(*d1,fabs(fd-gd));

}

pearsn(x1,y1,n1,&r1,&dum,&dumm);

sqen=sqrt((double)n1);

rr=sqrt(1.0-r1*r1);

*prob=probks(*d1*sqen/(1.0+rr*(0.25-0.75/sqen)));
}
#undef NRANSI
