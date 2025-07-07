#include <math.h>
#include <stdio.h>
#define NRANSI
#include "nrutil.h"
#define WINDOW(j,a,b) (1.0-fabs((((j)-1)-(a))*(b)))       /* Bartlett */

void spctrm(FILE *fp, double p[], int m, int k, int ovrlap)
{

void four1(double data[], unsigned long nn, int isign);

int mm,m44,m43,m4,kk,joffn,joff,j2,j;

double w,facp,facm,*w1,*w2,sumw=0.0,den=0.0;


mm=m+m;

m43=(m4=mm+mm)+3;

m44=m43+1;

w1=dvector(1,m4);

w2=dvector(1,m);

facm=m;

facp=1.0/m;

for (j=1;j<=mm;j++) sumw += DSQR(WINDOW(j,facm,facp));

for (j=1;j<=m;j++) p[j]=0.0;

if (ovrlap)


for (j=1;j<=m;j++) fscanf(fp,"%f",&w2[j]);

for (kk=1;kk<=k;kk++) {


for (joff = -1;joff<=0;joff++) {



if (ovrlap) {




for (j=1;j<=m;j++) w1[joff+j+j]=w2[j];




for (j=1;j<=m;j++) fscanf(fp,"%f",&w2[j]);




joffn=joff+mm;




for (j=1;j<=m;j++) w1[joffn+j+j]=w2[j];



} else {




for (j=joff+2;j<=m4;j+=2)





fscanf(fp,"%f",&w1[j]);



}


}


for (j=1;j<=mm;j++) {



j2=j+j;



w=WINDOW(j,facm,facp);



w1[j2] *= w;



w1[j2-1] *= w;


}


four1(w1,mm,1);


p[1] += (DSQR(w1[1])+DSQR(w1[2]));


for (j=2;j<=m;j++) {



j2=j+j;



p[j] += (DSQR(w1[j2])+DSQR(w1[j2-1])




+DSQR(w1[m44-j2])+DSQR(w1[m43-j2]));


}


den += sumw;

}

den *= m4;

for (j=1;j<=m;j++) p[j] /= den;

free_dvector(w2,1,m);

free_dvector(w1,1,m4);
}
#undef WINDOW
#undef NRANSI
