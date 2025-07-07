#define NRANSI
#include "nrutil.h"

void mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],

int ma, double **covar, double **alpha, double *chisq,

void (*funcs)(double, double [], double *, double [], int), double *alamda)
{

void covsrt(double **covar, int ma, int ia[], int mfit);

void gaussj(double **a, int n, double **b, int m);

void mrqcof(double x[], double y[], double sig[], int ndata, double a[],


int ia[], int ma, double **alpha, double beta[], double *chisq,


void (*funcs)(double, double [], double *, double [], int));

int j,k,l,m;

static int mfit;

static double ochisq,*atry,*beta,*da,**oneda;


if (*alamda < 0.0) {


atry=dvector(1,ma);


beta=dvector(1,ma);


da=dvector(1,ma);


for (mfit=0,j=1;j<=ma;j++)



if (ia[j]) mfit++;


oneda=dmatrix(1,mfit,1,1);


*alamda=0.001;


mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);


ochisq=(*chisq);


for (j=1;j<=ma;j++) atry[j]=a[j];

}

for (j=0,l=1;l<=ma;l++) {


if (ia[l]) {



for (j++,k=0,m=1;m<=ma;m++) {




if (ia[m]) {





k++;





covar[j][k]=alpha[j][k];




}



}



covar[j][j]=alpha[j][j]*(1.0+(*alamda));



oneda[j][1]=beta[j];


}

}

gaussj(covar,mfit,oneda,1);

for (j=1;j<=mfit;j++) da[j]=oneda[j][1];

if (*alamda == 0.0) {


covsrt(covar,ma,ia,mfit);


free_dmatrix(oneda,1,mfit,1,1);


free_dvector(da,1,ma);


free_dvector(beta,1,ma);


free_dvector(atry,1,ma);


return;

}

for (j=0,l=1;l<=ma;l++)


if (ia[l]) atry[l]=a[l]+da[++j];

mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);

if (*chisq < ochisq) {


*alamda *= 0.1;


ochisq=(*chisq);


for (j=0,l=1;l<=ma;l++) {



if (ia[l]) {




for (j++,k=0,m=1;m<=ma;m++) {





if (ia[m]) {






k++;






alpha[j][k]=covar[j][k];





}




}




beta[j]=da[j];




a[l]=atry[l];



}


}

} else {


*alamda *= 10.0;


*chisq=ochisq;

}
}
#undef NRANSI
