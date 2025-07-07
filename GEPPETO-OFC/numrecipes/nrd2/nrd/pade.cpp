#include <math.h>
#define NRANSI
#include "nrutil.h"
#define BIG 1.0e30

void pade(double cof[], int n, double *resid)
{

void lubksb(double **a, int n, int *indx, double b[]);

void ludcmp(double **a, int n, int *indx, double *d);

void mprove(double **a, double **alud, int n, int indx[], double b[],


double x[]);

int j,k,*indx;

double d,rr,rrold,sum,**q,**qlu,*x,*y,*z;


indx=ivector(1,n);

q=dmatrix(1,n,1,n);

qlu=dmatrix(1,n,1,n);

x=dvector(1,n);

y=dvector(1,n);

z=dvector(1,n);

for (j=1;j<=n;j++) {


y[j]=x[j]=cof[n+j];


for (k=1;k<=n;k++) {



q[j][k]=cof[j-k+n];



qlu[j][k]=q[j][k];


}

}

ludcmp(qlu,n,indx,&d);

lubksb(qlu,n,indx,x);

rr=BIG;

for (;;) {


rrold=rr;


for (j=1;j<=n;j++) z[j]=x[j];


mprove(q,qlu,n,indx,y,x);


for (rr=0.0,j=1;j<=n;j++)



rr += DSQR(z[j]-x[j]);


if (rr >= rrold) break;

}

*resid=sqrt(rr);

for (k=1;k<=n;k++) {


for (sum=cof[k],j=1;j<=k;j++) sum -= x[j]*cof[k-j];


y[k]=sum;

}

for (j=1;j<=n;j++) {


cof[j]=y[j];


cof[j+n] = -x[j];

}

free_dvector(z,1,n);

free_dvector(y,1,n);

free_dvector(x,1,n);

free_dmatrix(qlu,1,n,1,n);

free_dmatrix(q,1,n,1,n);

free_ivector(indx,1,n);
}
#undef BIG
#undef NRANSI
