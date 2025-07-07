#include <math.h>
#define NRANSI
#include "nrutil.h"

extern long idum;
extern double tt;

double amotsa(double **p, double y[], double psum[], int ndim, double pb[],

double *yb, double (*funk)(double []), int ihi, double *yhi, double fac)
{

double ran1(long *idum);

int j;

double fac1,fac2,yflu,ytry,*ptry;


ptry=dvector(1,ndim);

fac1=(1.0-fac)/ndim;

fac2=fac1-fac;

for (j=1;j<=ndim;j++)


ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;

ytry=(*funk)(ptry);

if (ytry <= *yb) {


for (j=1;j<=ndim;j++) pb[j]=ptry[j];


*yb=ytry;

}

yflu=ytry-tt*log(ran1(&idum));

if (yflu < *yhi) {


y[ihi]=ytry;


*yhi=yflu;


for (j=1;j<=ndim;j++) {



psum[j] += ptry[j]-p[ihi][j];



p[ihi][j]=ptry[j];


}

}

free_dvector(ptry,1,ndim);

return yflu;
}
#undef NRANSI
