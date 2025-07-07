#define NRANSI
#include "nrutil.h"

void stoerm(double y[], double d2y[], int nv, double xs, double htot, int nstep,

double yout[], void (*derivs)(double, double [], double []))
{

int i,n,neqns,nn;

double h,h2,halfh,x,*ytemp;


ytemp=dvector(1,nv);

h=htot/nstep;

halfh=0.5*h;

neqns=nv/2;

for (i=1;i<=neqns;i++) {


n=neqns+i;


ytemp[i]=y[i]+(ytemp[n]=h*(y[n]+halfh*d2y[i]));

}

x=xs+h;

(*derivs)(x,ytemp,yout);

h2=h*h;

for (nn=2;nn<=nstep;nn++) {


for (i=1;i<=neqns;i++)



ytemp[i] += (ytemp[(n=neqns+i)] += h2*yout[i]);


x += h;


(*derivs)(x,ytemp,yout);

}

for (i=1;i<=neqns;i++) {


n=neqns+i;


yout[n]=ytemp[n]/h+halfh*yout[i];


yout[i]=ytemp[i];

}

free_dvector(ytemp,1,nv);
}
#undef NRANSI
