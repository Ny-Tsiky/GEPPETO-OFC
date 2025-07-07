#define NRANSI
#include "nrutil.h"

void fred2(int n, double a, double b, double t[], double f[], double w[],

double (*g)(double), double (*ak)(double, double))
{

void gauleg(double x1, double x2, double x[], double w[], int n);

void lubksb(double **a, int n, int *indx, double b[]);

void ludcmp(double **a, int n, int *indx, double *d);

int i,j,*indx;

double d,**omk;


indx=ivector(1,n);

omk=dmatrix(1,n,1,n);

gauleg(a,b,t,w,n);

for (i=1;i<=n;i++) {


for (j=1;j<=n;j++)



omk[i][j]=(double)(i == j)-(*ak)(t[i],t[j])*w[j];


f[i]=(*g)(t[i]);

}

ludcmp(omk,n,indx,&d);

lubksb(omk,n,indx,f);

free_dmatrix(omk,1,n,1,n);

free_ivector(indx,1,n);
}
#undef NRANSI
