#include <math.h>
#define RTPIO2 1.2533141

void sphbes(int n, double x, double *sj, double *sy, double *sjp, double *syp)
{

void bessjy(double x, double xnu, double *rj, double *ry, double *rjp,


double *ryp);

void nrerror(char error_text[]);

double factor,order,rj,rjp,ry,ryp;


if (n < 0 || x <= 0.0) nrerror("bad arguments in sphbes");

order=n+0.5;

bessjy(x,order,&rj,&ry,&rjp,&ryp);

factor=RTPIO2/sqrt(x);

*sj=factor*rj;

*sy=factor*ry;

*sjp=factor*rjp-(*sj)/(2.0*x);

*syp=factor*ryp-(*sy)/(2.0*x);
}
#undef RTPIO2
