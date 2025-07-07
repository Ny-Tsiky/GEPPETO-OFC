#include <math.h>
#define NRANSI
#include "nrutil.h"
#define PI 3.14159265

double x;

void quadmx(double **a, int n)
{

void kermom(double w[], double y, int m);

void wwghts(double wghts[], int n, double h,


void (*kermom)(double [], double ,int));

int j,k;

double h,*wt,xx,cx;


wt=dvector(1,n);

h=PI/(n-1);

for (j=1;j<=n;j++) {


x=xx=(j-1)*h;


wwghts(wt,n,h,kermom);


cx=cos(xx);


for (k=1;k<=n;k++) a[j][k]=wt[k]*cx*cos((k-1)*h);


++a[j][j];

}

free_dvector(wt,1,n);
}
#undef PI
#undef NRANSI
