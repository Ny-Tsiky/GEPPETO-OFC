#define NRANSI
#include "nrutil.h"

void cyclic(double a[], double b[], double c[], double alpha, double beta,

double r[], double x[], unsigned long n)
{

void tridag(double a[], double b[], double c[], double r[], double u[],


unsigned long n);

unsigned long i;

double fact,gamma,*bb,*u,*z;


if (n <= 2) nrerror("n too small in cyclic");

bb=dvector(1,n);

u=dvector(1,n);

z=dvector(1,n);

gamma = -b[1];

bb[1]=b[1]-gamma;

bb[n]=b[n]-alpha*beta/gamma;

for (i=2;i<n;i++) bb[i]=b[i];

tridag(a,bb,c,r,x,n);

u[1]=gamma;

u[n]=alpha;

for (i=2;i<n;i++) u[i]=0.0;

tridag(a,bb,c,u,z,n);

fact=(x[1]+beta*x[n]/gamma)/


(1.0+z[1]+beta*z[n]/gamma);

for (i=1;i<=n;i++) x[i] -= fact*z[i];

free_dvector(z,1,n);

free_dvector(u,1,n);

free_dvector(bb,1,n);
}
#undef NRANSI
