#include <math.h>
#define NRANSI
#include "nrutil.h"
#define PI 3.141592653589793

void chebft(double a, double b, double c[], int n, double (*func)(double))
{

int k,j;

double fac,bpa,bma,*f;


f=dvector(0,n-1);

bma=0.5*(b-a);

bpa=0.5*(b+a);

for (k=0;k<n;k++) {


double y=cos(PI*(k+0.5)/n);


f[k]=(*func)(y*bma+bpa);

}

fac=2.0/n;

for (j=0;j<n;j++) {


double sum=0.0;


for (k=0;k<n;k++)



sum += f[k]*cos(PI*j*(k+0.5)/n);


c[j]=fac*sum;

}

free_dvector(f,0,n-1);
}
#undef PI
#undef NRANSI
