#include <math.h>
#define TINY 1.0e-20

void pearsn(double x[], double y[], unsigned long n, double *r, double *prob,

double *z)
{

double betai(double a, double b, double x);

double erfcc(double x);

unsigned long j;

double yt,xt,t,df;

double syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;


for (j=1;j<=n;j++) {


ax += x[j];


ay += y[j];

}

ax /= n;

ay /= n;

for (j=1;j<=n;j++) {


xt=x[j]-ax;


yt=y[j]-ay;


sxx += xt*xt;


syy += yt*yt;


sxy += xt*yt;

}

*r=sxy/sqrt(sxx*syy);

*z=0.5*log((1.0+(*r)+TINY)/(1.0-(*r)+TINY));

df=n-2;

t=(*r)*sqrt(df/((1.0-(*r)+TINY)*(1.0+(*r)+TINY)));

*prob=betai(0.5*df,0.5,df/(df+t*t));
}
#undef TINY
