#include <math.h>
#define NRANSI
#include "nrutil.h"
#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0

double dfridr(double (*func)(double), double x, double h, double *err)
{

int i,j;

double errt,fac,hh,**a,ans;


if (h == 0.0) nrerror("h must be nonzero in dfridr.");

a=dmatrix(1,NTAB,1,NTAB);

hh=h;

a[1][1]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);

*err=BIG;

for (i=2;i<=NTAB;i++) {


hh /= CON;


a[1][i]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);


fac=CON2;


for (j=2;j<=i;j++) {



a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);



fac=CON2*fac;



errt=DMAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));



if (errt <= *err) {




*err=errt;




ans=a[j][i];



}


}


if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) {



free_dmatrix(a,1,NTAB,1,NTAB);



return ans;


}

}

free_dmatrix(a,1,NTAB,1,NTAB);

return ans;
}
#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE
#undef NRANSI
