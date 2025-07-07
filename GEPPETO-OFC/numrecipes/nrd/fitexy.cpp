#include <math.h>
#define NRANSI
#include "nrutil.h"
#define POTN 1.571000
#define BIG 1.0e30
#define PI 3.14159265
#define ACC 1.0e-3

// JNI 16/Jun/2000 Added "extern" (already defined in module: NEWT.c)
extern int    nn;
// JNI 16/Jun/2000 Added "extern" (already defined in module: RKDUMB.c)
extern double *xx;
// JNI 16/Jun/2000 Added "extern" (already defined in module: MEDFIT.c)
extern double  aa;

double *yy,*sx,*sy,*ww,offs;

void fitexy(double x[], double y[], int ndat, double sigx[], double sigy[],

double *a, double *b, double *siga, double *sigb, double *chi2, double *q)
{

void avevar(double data[], unsigned long n, double *ave, double *var);

double brent(double ax, double bx, double cx,


double (*f)(double), double tol, double *xmin);

double chixy(double bang);

void fit(double x[], double y[], int ndata, double sig[], int mwt,


double *a, double *b, double *siga, double *sigb, double *chi2, double *q);

double gammq(double a, double x);

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,


double *fc, double (*func)(double));

double zbrent(double (*func)(double), double x1, double x2, double tol);

int j;

double swap,amx,amn,varx,vary,ang[7],ch[7],scale,bmn,bmx,d1,d2,r2,


dum1,dum2,dum3,dum4,dum5;


xx=dvector(1,ndat);

yy=dvector(1,ndat);

sx=dvector(1,ndat);

sy=dvector(1,ndat);

ww=dvector(1,ndat);

avevar(x,ndat,&dum1,&varx);

avevar(y,ndat,&dum1,&vary);

scale=sqrt(varx/vary);

nn=ndat;

for (j=1;j<=ndat;j++) {


xx[j]=x[j];


yy[j]=y[j]*scale;


sx[j]=sigx[j];


sy[j]=sigy[j]*scale;


ww[j]=sqrt(DSQR(sx[j])+DSQR(sy[j]));

}

fit(xx,yy,nn,ww,1,&dum1,b,&dum2,&dum3,&dum4,&dum5);

offs=ang[1]=0.0;

ang[2]=atan(*b);

ang[4]=0.0;

ang[5]=ang[2];

ang[6]=POTN;

for (j=4;j<=6;j++) ch[j]=chixy(ang[j]);

mnbrak(&ang[1],&ang[2],&ang[3],&ch[1],&ch[2],&ch[3],chixy);

*chi2=brent(ang[1],ang[2],ang[3],chixy,ACC,b);

*chi2=chixy(*b);

*a=aa;

*q=gammq(0.5*(nn-2),*chi2*0.5);

for (r2=0.0,j=1;j<=nn;j++) r2 += ww[j];

r2=1.0/r2;

bmx=BIG;

bmn=BIG;

offs=(*chi2)+1.0;

for (j=1;j<=6;j++) {


if (ch[j] > offs) {



d1=fabs(ang[j]-(*b));



while (d1 >= PI) d1 -= PI;



d2=PI-d1;



if (ang[j] < *b) {




swap=d1;




d1=d2;




d2=swap;



}



if (d1 < bmx) bmx=d1;



if (d2 < bmn) bmn=d2;


}

}

if (bmx < BIG) {


bmx=zbrent(chixy,*b,*b+bmx,ACC)-(*b);


amx=aa-(*a);


bmn=zbrent(chixy,*b,*b-bmn,ACC)-(*b);


amn=aa-(*a);


*sigb=sqrt(0.5*(bmx*bmx+bmn*bmn))/(scale*DSQR(cos(*b)));


*siga=sqrt(0.5*(amx*amx+amn*amn)+r2)/scale;

} else (*sigb)=(*siga)=BIG;

*a /= scale;

*b=tan(*b)/scale;

free_dvector(ww,1,ndat);

free_dvector(sy,1,ndat);

free_dvector(sx,1,ndat);

free_dvector(yy,1,ndat);

free_dvector(xx,1,ndat);
}
#undef POTN
#undef BIG
#undef PI
#undef ACC
#undef NRANSI
