#include <math.h>
#define NRANSI
#include "nrutil.h"

void usrfun(double *x,int n,double *fvec,double **fjac);
#define FREETHINGS free_dmatrix(fjac,1,n,1,n),free_dvector(fvec,1,n),free_dvector(p,1,n),free_ivector(indx,1,n);

void mnewt(int ntrial, double x[], int n, double tolx, double tolf)
{

void lubksb(double **a, int n, int *indx, double b[]);

void ludcmp(double **a, int n, int *indx, double *d);

int k,i,*indx;

double errx,errf,d,*fvec,**fjac,*p;


indx=ivector(1,n);

p=dvector(1,n);

fvec=dvector(1,n);

fjac=dmatrix(1,n,1,n);

for (k=1;k<=ntrial;k++) {


usrfun(x,n,fvec,fjac);


errf=0.0;


for (i=1;i<=n;i++) errf += fabs(fvec[i]);


if (errf <= tolf) { FREETHINGS return; }


for (i=1;i<=n;i++) p[i] = -fvec[i];


ludcmp(fjac,n,indx,&d);


lubksb(fjac,n,indx,p);


errx=0.0;


for (i=1;i<=n;i++) {



errx += fabs(p[i]);



x[i] += p[i];


}


if (errx <= tolx) { FREETHINGS return; }

}

FREETHINGS
}
#undef FREETHINGS
#undef NRANSI
