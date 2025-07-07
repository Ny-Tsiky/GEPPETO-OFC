#include <math.h>
#define NRANSI
#include "nrutil.h"
#define TINY 1.0e-30

void cntab2(int **nn, int ni, int nj, double *h, double *hx, double *hy,

double *hygx, double *hxgy, double *uygx, double *uxgy, double *uxy)
{

int i,j;

double sum=0.0,p,*sumi,*sumj;


sumi=dvector(1,ni);

sumj=dvector(1,nj);

for (i=1;i<=ni;i++) {


sumi[i]=0.0;


for (j=1;j<=nj;j++) {



sumi[i] += nn[i][j];



sum += nn[i][j];


}

}

for (j=1;j<=nj;j++) {


sumj[j]=0.0;


for (i=1;i<=ni;i++)



sumj[j] += nn[i][j];

}

*hx=0.0;

for (i=1;i<=ni;i++)


if (sumi[i]) {



p=sumi[i]/sum;



*hx -= p*log(p);


}

*hy=0.0;

for (j=1;j<=nj;j++)


if (sumj[j]) {



p=sumj[j]/sum;



*hy -= p*log(p);


}

*h=0.0;

for (i=1;i<=ni;i++)


for (j=1;j<=nj;j++)



if (nn[i][j]) {




p=nn[i][j]/sum;




*h -= p*log(p);



}

*hygx=(*h)-(*hx);

*hxgy=(*h)-(*hy);

*uygx=(*hy-*hygx)/(*hy+TINY);

*uxgy=(*hx-*hxgy)/(*hx+TINY);

*uxy=2.0*(*hx+*hy-*h)/(*hx+*hy+TINY);

free_dvector(sumj,1,nj);

free_dvector(sumi,1,ni);
}
#undef TINY
#undef NRANSI
