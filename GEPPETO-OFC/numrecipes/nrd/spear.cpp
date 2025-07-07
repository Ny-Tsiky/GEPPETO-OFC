#include <math.h>
#define NRANSI
#include "nrutil.h"

void spear(double data1[], double data2[], unsigned long n, double *d, double *zd,

double *probd, double *rs, double *probrs)
{

double betai(double a, double b, double x);

void crank(unsigned long n, double w[], double *s);

double erfcc(double x);

void sort2(unsigned long n, double arr[], double brr[]);

unsigned long j;

double vard,t,sg,sf,fac,en3n,en,df,aved,*wksp1,*wksp2;


wksp1=dvector(1,n);

wksp2=dvector(1,n);

for (j=1;j<=n;j++) {


wksp1[j]=data1[j];


wksp2[j]=data2[j];

}

sort2(n,wksp1,wksp2);

crank(n,wksp1,&sf);

sort2(n,wksp2,wksp1);

crank(n,wksp2,&sg);

*d=0.0;

for (j=1;j<=n;j++)


*d += DSQR(wksp1[j]-wksp2[j]);

en=n;

en3n=en*en*en-en;

aved=en3n/6.0-(sf+sg)/12.0;

fac=(1.0-sf/en3n)*(1.0-sg/en3n);

vard=((en-1.0)*en*en*DSQR(en+1.0)/36.0)*fac;

*zd=(*d-aved)/sqrt(vard);

*probd=erfcc(fabs(*zd)/1.4142136);

*rs=(1.0-(6.0/en3n)*(*d+(sf+sg)/12.0))/sqrt(fac);

fac=(*rs+1.0)*(1.0-(*rs));

if (fac > 0.0) {


t=(*rs)*sqrt((en-2.0)/fac);


df=en-2.0;


*probrs=betai(0.5*df,0.5,df/(df+t*t));

} else


*probrs=0.0;

free_dvector(wksp2,1,n);

free_dvector(wksp1,1,n);
}
#undef NRANSI
