#define NRANSI
#include "nrutil.h"

void voltra(int n, int m, double t0, double h, double *t, double **f,

double (*g)(int, double), double (*ak)(int, int, double, double))
{

void lubksb(double **a, int n, int *indx, double b[]);

void ludcmp(double **a, int n, int *indx, double *d);

int i,j,k,l,*indx;

double d,sum,**a,*b;


indx=ivector(1,m);

a=dmatrix(1,m,1,m);

b=dvector(1,m);

t[1]=t0;

for (k=1;k<=m;k++) f[k][1]=(*g)(k,t[1]);

for (i=2;i<=n;i++) {


t[i]=t[i-1]+h;


for (k=1;k<=m;k++) {



sum=(*g)(k,t[i]);



for (l=1;l<=m;l++) {




sum += 0.5*h*(*ak)(k,l,t[i],t[1])*f[l][1];




for (j=2;j<i;j++)





sum += h*(*ak)(k,l,t[i],t[j])*f[l][j];




a[k][l]=(k == l)-0.5*h*(*ak)(k,l,t[i],t[i]);



}



b[k]=sum;


}


ludcmp(a,m,indx,&d);


lubksb(a,m,indx,b);


for (k=1;k<=m;k++) f[k][i]=b[k];

}

free_dvector(b,1,m);

free_dmatrix(a,1,m,1,m);

free_ivector(indx,1,m);
}
#undef NRANSI
