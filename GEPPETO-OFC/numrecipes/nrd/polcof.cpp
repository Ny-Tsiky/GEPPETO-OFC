#include <math.h>
#define NRANSI
#include "nrutil.h"

void polcof(double xa[], double ya[], int n, double cof[])
{

void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

int k,j,i;

double xmin,dy,*x,*y;


x=dvector(0,n);

y=dvector(0,n);

for (j=0;j<=n;j++) {


x[j]=xa[j];


y[j]=ya[j];

}

for (j=0;j<=n;j++) {


polint(x-1,y-1,n+1-j,0.0,&cof[j],&dy);


xmin=1.0e38;


k = -1;


for (i=0;i<=n-j;i++) {



if (fabs(x[i]) < xmin) {




xmin=fabs(x[i]);




k=i;



}



if (x[i]) y[i]=(y[i]-cof[j])/x[i];


}


for (i=k+1;i<=n-j;i++) {



y[i-1]=y[i];



x[i-1]=x[i];


}

}

free_dvector(y,0,n);

free_dvector(x,0,n);
}
#undef NRANSI
