#include <math.h>
#define NRANSI
#include "nrutil.h"
//#define GET_PSUM for (n=1;n<=ndim;n++) {\
//for (sum=0.0,m=1;m<=mpts;m++) sum += p[m][n];\
//psum[n]=sum;}
extern long idum;
double tt;

void GET_PSUM(double **p, double *psum, int ndim)
{
  int m,n;
  double sum;
  for (n=1;n<=ndim;n++) {
    for (sum=0.0,m=1;m<=ndim+1;m++) sum += p[m][n];
    psum[n]=sum;
  }
}


void amebsa(double **p, double y[], int ndim, double pb[], double *yb, double ftol,double (*funk)(double []), int *iter, double temptr)
{
  double amotsa(double **p, double y[], double psum[], int ndim, double pb[],double *yb, double (*funk)(double []), int ihi, double *yhi, double fac);
  
  double ran1(long *idum);
  int i,ihi,ilo,j,m,n,mpts=ndim+1;
  double rtol,sum,swap,yhi,ylo,ynhi,ysave,yt,ytry,*psum;
  psum=dvector(1,ndim);
  tt = -temptr;
  GET_PSUM(p,psum,ndim);
    for (;;) {
      ilo=1;
      ihi=2;
      ynhi=ylo=y[1]+tt*log(ran1(&idum));
      yhi=y[2]+tt*log(ran1(&idum));
      if (ylo > yhi) {
	ihi=1;
	ilo=2;
	ynhi=yhi;
	yhi=ylo;
	ylo=ynhi;
      }
      for (i=3;i<=mpts;i++) {
	yt=y[i]+tt*log(ran1(&idum));
	if (yt <= ylo) {
	  ilo=i;
	  ylo=yt;
	}
	if (yt > yhi) {
	  ynhi=yhi;
	  ihi=i;
	  yhi=yt;
	} else if (yt > ynhi) {
	  ynhi=yt;
	}
      }
      rtol=2.0*fabs(yhi-ylo)/(fabs(yhi)+fabs(ylo));
      if (rtol < ftol || *iter < 0) {
	swap=y[1];
	y[1]=y[ilo];
	y[ilo]=swap;
	for (n=1;n<=ndim;n++) {
	  swap=p[1][n];
	  p[1][n]=p[ilo][n];
	  p[ilo][n]=swap;
	}
	break;
      }
      *iter -= 2;
      ytry=amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,-1.0);
      if (ytry <= ylo) {
	ytry=amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,2.0);
      } else if (ytry >= ynhi) {
	ysave=yhi;
	ytry=amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,0.5);
	if (ytry >= ysave) {
	  for (i=1;i<=mpts;i++) {
	    if (i != ilo) {
	      for (j=1;j<=ndim;j++) {
		psum[j]=0.5*(p[i][j]+p[ilo][j]);
		p[i][j]=psum[j];
	      }
	      y[i]=(*funk)(psum);
	    }
	  }
	  *iter -= ndim;
	  GET_PSUM(p,psum,ndim);
	    }
      } else ++(*iter);
    }
  free_dvector(psum,1,ndim);
}
#undef GET_PSUM
#undef NRANSI
