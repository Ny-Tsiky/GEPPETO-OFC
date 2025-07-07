#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0

#define FREEALL free_dvector(xi,1,n);free_dvector(pnew,1,n);free_dmatrix(hessin,1,n,1,n);free_dvector(hdg,1,n);free_dvector(g,1,n);free_dvector(dg,1,n);

void dfpmin(double p[], int n, double gtol, int *iter, double *fret,

double(*func)(double []), void (*dfunc)(double [], double []))
{

void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],


 double *f, double stpmax, int *check, double (*func)(double []));

int check,i,its,j;

double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;

double *dg,*g,*hdg,**hessin,*pnew,*xi;


dg=dvector(1,n);

g=dvector(1,n);

hdg=dvector(1,n);

hessin=dmatrix(1,n,1,n);

pnew=dvector(1,n);

xi=dvector(1,n);

fp=(*func)(p);

(*dfunc)(p,g);

for (i=1;i<=n;i++) {


for (j=1;j<=n;j++) hessin[i][j]=0.0;


hessin[i][i]=1.0;


xi[i] = -g[i];


sum += p[i]*p[i];

}

stpmax=STPMX*DMAX(sqrt(sum),(double)n);

for (its=1;its<=ITMAX;its++) {


*iter=its;


lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,func);


fp = *fret;


for (i=1;i<=n;i++) {



xi[i]=pnew[i]-p[i];



p[i]=pnew[i];


}


test=0.0;


for (i=1;i<=n;i++) {



temp=fabs(xi[i])/DMAX(fabs(p[i]),1.0);



if (temp > test) test=temp;


}


if (test < TOLX) {



FREEALL



return;


}


for (i=1;i<=n;i++) dg[i]=g[i];


(*dfunc)(p,g);


test=0.0;


den=DMAX(*fret,1.0);


for (i=1;i<=n;i++) {



temp=fabs(g[i])*DMAX(fabs(p[i]),1.0)/den;



if (temp > test) test=temp;


}


if (test < gtol) {



FREEALL



return;


}


for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];


for (i=1;i<=n;i++) {



hdg[i]=0.0;



for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];


}


fac=fae=sumdg=sumxi=0.0;


for (i=1;i<=n;i++) {



fac += dg[i]*xi[i];



fae += dg[i]*hdg[i];



sumdg += DSQR(dg[i]);



sumxi += DSQR(xi[i]);


}


if (fac*fac > EPS*sumdg*sumxi) {



fac=1.0/fac;



fad=1.0/fae;



for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];



for (i=1;i<=n;i++) {




for (j=1;j<=n;j++) {





hessin[i][j] += fac*xi[i]*xi[j]





-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];




}



}


}


for (i=1;i<=n;i++) {



xi[i]=0.0;



for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];


}

}

nrerror("too many iterations in dfpmin");

FREEALL
}
#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX
#undef FREEALL
#undef NRANSI
