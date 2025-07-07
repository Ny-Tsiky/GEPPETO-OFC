#include <math.h>
#define NRANSI
#include "nrutil.h"

void memcof(double data[], int n, int m, double *xms, double d[])
{

int k,j,i;

double p=0.0,*wk1,*wk2,*wkm;


wk1=dvector(1,n);

wk2=dvector(1,n);

wkm=dvector(1,m);

for (j=1;j<=n;j++) p += DSQR(data[j]);

*xms=p/n;

wk1[1]=data[1];

wk2[n-1]=data[n];

for (j=2;j<=n-1;j++) {


wk1[j]=data[j];


wk2[j-1]=data[j];

}

for (k=1;k<=m;k++) {


double num=0.0,denom=0.0;


for (j=1;j<=(n-k);j++) {



num += wk1[j]*wk2[j];



denom += DSQR(wk1[j])+DSQR(wk2[j]);


}


d[k]=2.0*num/denom;


*xms *= (1.0-DSQR(d[k]));


for (i=1;i<=(k-1);i++)



d[i]=wkm[i]-d[k]*wkm[k-i];


if (k == m) {



free_dvector(wkm,1,m);



free_dvector(wk2,1,n);



free_dvector(wk1,1,n);



return;


}


for (i=1;i<=k;i++) wkm[i]=d[i];


for (j=1;j<=(n-k-1);j++) {



wk1[j] -= wkm[k]*wk2[j];



wk2[j]=wk2[j+1]-wkm[k]*wk1[j+1];


}

}

nrerror("never get here in memcof.");
}
#undef NRANSI
