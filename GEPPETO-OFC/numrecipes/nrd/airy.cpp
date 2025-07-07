#include <math.h>
#define PI 3.1415927
#define THIRD (1.0/3.0)
#define TWOTHR (2.0*THIRD)
#define ONOVRT 0.57735027

void airy(double x, double *ai, double *bi, double *aip, double *bip)
{

void bessik(double x, double xnu, double *ri, double *rk, double *rip,


double *rkp);

void bessjy(double x, double xnu, double *rj, double *ry, double *rjp,


double *ryp);

double absx,ri,rip,rj,rjp,rk,rkp,rootx,ry,ryp,z;


absx=fabs(x);

rootx=sqrt(absx);

z=TWOTHR*absx*rootx;

if (x > 0.0) {


bessik(z,THIRD,&ri,&rk,&rip,&rkp);


*ai=rootx*ONOVRT*rk/PI;


*bi=rootx*(rk/PI+2.0*ONOVRT*ri);


bessik(z,TWOTHR,&ri,&rk,&rip,&rkp);


*aip = -x*ONOVRT*rk/PI;


*bip=x*(rk/PI+2.0*ONOVRT*ri);

} else if (x < 0.0) {


bessjy(z,THIRD,&rj,&ry,&rjp,&ryp);


*ai=0.5*rootx*(rj-ONOVRT*ry);


*bi = -0.5*rootx*(ry+ONOVRT*rj);


bessjy(z,TWOTHR,&rj,&ry,&rjp,&ryp);


*aip=0.5*absx*(ONOVRT*ry+rj);


*bip=0.5*absx*(ONOVRT*rj-ry);

} else {


*ai=0.35502805;


*bi=(*ai)/ONOVRT;


*aip = -0.25881940;


*bip = -(*aip)/ONOVRT;

}
}
#undef PI
#undef THIRD
#undef TWOTHR
#undef ONOVRT
