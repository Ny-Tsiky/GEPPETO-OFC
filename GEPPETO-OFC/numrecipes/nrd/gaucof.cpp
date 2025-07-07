#include <math.h>
#define NRANSI
#include "nrutil.h"

void gaucof(int n, double a[], double b[], double amu0, double x[], double w[])
{

void eigsrt(double d[], double **v, int n);

void tqli(double d[], double e[], int n, double **z);

int i,j;

double **z;


z=dmatrix(1,n,1,n);

for (i=1;i<=n;i++) {


if (i != 1) b[i]=sqrt(b[i]);


for (j=1;j<=n;j++) z[i][j]=(double)(i == j);

}

tqli(a,b,n,z);

eigsrt(a,z,n);

for (i=1;i<=n;i++) {


x[i]=a[i];


w[i]=amu0*z[1][i]*z[1][i];

}

free_dmatrix(z,1,n,1,n);
}
#undef NRANSI
