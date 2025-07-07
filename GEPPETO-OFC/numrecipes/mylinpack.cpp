#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <conio.h>  // if Windows really needs this, put an #ifdef _WIN32
#include <ctype.h>
//#include <fstream.h> // test
#include "nr.h"
#include "nrutil.h"
#include "mylinpack.h"

#define NO_MATLAB

/*
extern double **yp;
extern double *xp;
extern int kount;
*/

#ifndef NO_MATLAB
#include "mex.h"  // comment out for pure C++
#undef DB
#define DB mexPrintf("Got to line %d in mex file %s.\n",__LINE__,__FILE__);

#endif



static const double rel_error= 1E-12;        
// ERF and ERFC implementation from Steve Strand, Jan 2004
//calculate 12 significant figures
//you can adjust rel_error to trade off between accuracy and speed
//but don't ask for > 15 figures (assuming usual 52 bit mantissa in a double)


double erf(double x)
//erf(x) = 2/sqrt(pi)*integral(exp(-t^2),t,0,x)
//       = 2/sqrt(pi)*[x - x^3/3 + x^5/5*2! - x^7/7*3! + ...]
//       = 1-erfc(x)
{
    static const double two_sqrtpi=  1.128379167095512574;        // 2/sqrt(pi)
    if (fabs(x) > 2.2) {
        return 1.0 - erfc(x);        //use continued fraction when fabs(x) > 2.2
    }
    double sum= x, term= x, xsqr= x*x;
    int j= 1;
    do {
        term*= xsqr/j;
        sum-= term/(2*j+1);
        ++j;
        term*= xsqr/j;
        sum+= term/(2*j+1);
        ++j;
    } while (fabs(term/sum) > rel_error);
    return two_sqrtpi*sum;
}


double erfc(double x)
//erfc(x) = 2/sqrt(pi)*integral(exp(-t^2),t,x,inf)
//        = exp(-x^2)/sqrt(pi) * [1/x+ (1/2)/x+ (2/2)/x+ (3/2)/x+ (4/2)/x+ ...]
//        = 1-erf(x)
//expression inside [] is a continued fraction so '+' means add to denominator only
{
    static const double one_sqrtpi=  0.564189583547756287;        // 1/sqrt(pi)
    if (fabs(x) < 2.2) {
        return 1.0 - erf(x);        //use series when fabs(x) < 2.2
    }
    //    if (signbit(x)) {               //continued fraction only valid for x>0
    if (x < 0) {               //continued fraction only valid for x>0
        return 2.0 - erfc(-x);
    }
    double a=1, b=x;                //last two convergent numerators
    double c=x, d=x*x+0.5;          //last two convergent denominators
    double q1, q2= b/d;             //last two convergents (a/c and b/d)
    double n= 1.0, t;
    do {
        t= a*n+b*x;
        a= b;
        b= t;
        t= c*n+d*x;
        c= d;
        d= t;
        n+= 0.5;
        q1= q2;
        q2= b/d;
      } while (fabs(q1-q2)/q2 > rel_error);
    return one_sqrtpi*exp(-x*x)*q2;
}

/**
double randn2()
{
  return(0);
}//test
*/

double randn2()
// from erf()
{
  double z,zmin,zmax,err,r=rand()/(double)RAND_MAX;
  
  zmin=-100;
  zmax=100;
  err=1;
  while (fabs(err)>1e-10) {
    z=(zmin+zmax)/2;
    err=r - (1+erf(z/sqrt2))/2;
    if (err>0) zmin=z; else zmax=z;
  }
  return(z);
}




#include "rand.cpp"


double randn()
//return a normally distributed random number
//with mean 0 and standard deviation 1
//See Knuth, "Art of Computer Programming", vol. 2, 2nd edition, page 117
{
  double u1, u2, r1, s;
  static double r2;
  static bool r2_is_set= false;
  static const double two31        = 2147483648.0;	// 2^31
  static const double two31lesshalf= 2147483647.5;	// 2^32 - 1/2
  
  if (r2_is_set) {	//return random r2 generated on last call
    r2_is_set= false;
    return r2;
  }
  
  do {	//pick a point inside unit circle, but not (0,0) to avoid log(0)
    u1= (rand2()-two31lesshalf)/two31;	//in range -(1-2^-32) .. (1-2^-32)
    u2= (rand2()-two31lesshalf)/two31;
    s= u1*u1+u2*u2;
  } while (s>1.0);
  
  //generate two randoms r1 and r2, return r1 and save r2
  s= sqrt(-2.0*log(s)/s);
  r1= u1*s;
  r2= u2*s;
  
  r2_is_set= true;
  return r1;
}


void randn(double *v, int n)
{
  int i;

  for ( i=1; i<=n; i++) v[i] = randn();  
}

double *randn(int n)
{
  double *v = dvector(1,n);

  randn(v, n);
  return(v);
}

void randn(double **A, int m, int n)
{
  int i, j;

  for ( i=1; i<=m; i++) 
    for ( j=1; j<=n; j++) 
      A[i][j] = randn();
}


double **randn(int m, int n)
{
  double **A = dmatrix(1,m,1,n);

  randn(A, m, n);
  return(A);
}

void rand(double **A, int m, int n)
{
  int i, j;

  for ( i=1; i<=m; i++) 
    for ( j=1; j<=n; j++) 
      A[i][j] = rand();
}


double **rand(int m, int n)
{
  double **A = dmatrix(1,m,1,n);

  rand(A, m, n);
  return(A);
}


/*  Now as inline function
double sign(double x)
{
  if (x > EPSILON) return(1);
  else if (x < -EPSILON ) return(-1);
  else return(0);
}
*/

#ifdef _WIN32   // crappy Windows needs this

double rint( double x) 
// Copyright (C) 2001 Tor M. Aamodt, University of Toronto 
// Permisssion to use for all purposes commercial and otherwise granted. 
// THIS MATERIAL IS PROVIDED "AS IS" WITHOUT WARRANTY, OR ANY CONDITION OR 
// OTHER TERM OF ANY KIND INCLUDING, WITHOUT LIMITATION, ANY WARRANTY 
// OF MERCHANTABILITY, SATISFACTORY QUALITY, OR FITNESS FOR A PARTICULAR 
// PURPOSE. 
{ 
    if( x > 0 ) { 
        __int64 xint = (__int64) (x+0.5); 
        if( xint % 2 ) { 
            // then we might have an even number... 
            double diff = x - (double)xint; 
            if( diff == -0.5 ) 
                return double(xint-1); 
        } 
        return double(xint); 
    } else { 
        __int64 xint = (__int64) (x-0.5); 
        if( xint % 2 ) { 
            // then we might have an even number... 
            double diff = x - (double)xint; 
            if( diff == 0.5 ) 
                return double(xint+1); 
        } 
        return double(xint); 
    } 
}
#endif

/* MATRIX MULTIPLICATION */

void 
matmult(int m, double **A, int n,  double **B, int p, double**C)
{
  int i,j,k;
  double S;

  for (i=1; i<=m; i++) 
    for (j=1; j<=p; j++) {
      S=0;
      for (k=1; k<=n; k++) S+=A[i][k]*B[k][j];
      C[i][j]=S;
    }
}

void 
vecmatmult(double *u, int n,  double **B, int p, double *v)
{
  int j,k;
  double S;

  for (j=1; j<=p; j++) {
    S=0;
    for (k=1; k<=n; k++) S+=u[k]*B[k][j];
    v[j]=S;
  }
}

void 
matvecmult(int m, double **A, int n,  double *u, double *v)
{
  int i,k;
  double S;

  for (i=1; i<=m; i++) {
    S=0;
    for (k=1; k<=n; k++) S+=A[i][k]*u[k];
    v[i]=S;
  }
}

/* MATRIX TRANSPOSITION */

void 
transpose(double **A, int m, int n, double **B)
{
  int i,j;

  for (i=1; i<=m; i++)
    for (j=1; j<=n; j++) B[j][i]=A[i][j];
}

/* VECTOR SUBTRACTION */

void 
vectsub(double *u, double *v, int n, double *w)
{
  int i;

  for(i=1; i<=n; i++) w[i]=u[i]-v[i];
}

/* DOT PRODUCT */

double 
dotprod(double *u, double *v, int n)
{
  int i;
  double d=0;

  for(i=1; i<=n; i++) d+=u[i]*v[i];
  return(d);
}

/* "TOD" PRODUCT : u*v' */

void todprod(double *u, double *v, int n, double **A)
{
  int i,j;

  for(i=1; i<=n; i++) 
    for(j=1; j<=n; j++) 
      A[i][j] = u[i]*v[j];
}


double **
todprod(double *u, double *v, int n)
{
  double **A = dmatrix(1,n,1,n);

  todprod(u,v,n,A);
  return(A);
}

/* SQUARED EUCLIDIAN NORM */

double 
squared(double *v, int n)
{
  int i;
  double sq=0;

  for (i=1; i<=n; i++) sq+=v[i]*v[i];
  return(sq);
}


/* EUCLIDIAN NORM */

double 
norm2(double *v, int n)
{
  return(sqrt(squared(v,n)));
}

double
norm2(double  **A, int m, int n)
{
  int i;
  double sq=0;

  for (i=1; i<=m; i++) sq+=squared(A[i],n);
  return(sqrt(sq));
}

double 
sqnorm2(double *v, int n)
{
  return(squared(v,n));
}

double
sqnorm2(double  **A, int m, int n)
{
  int i;
  double sq=0;

  for (i=1; i<=m; i++) sq+=squared(A[i],n);
  return(sq);
}


/* GENERIC p-NORM (slower) */

double 
norm(double *v, int n, int p)
{
  int i;
  double norm=0;

  for (i=1; i<=n; i++) norm+=pow(fabs(v[i]),p); 
  return(pow(norm,1/p));
  
}

/* CROSS PRODUCT */

void
crossprod(double* u, double* v, double* w)
{
  int i;

  for (i=1; i<=3; i++) w[i]=u[i%3 + 1]*v[(i+1)%3+1] - v[i%3 + 1]*u[(i+1)%3+1];
}

/* SIGN INVERSION */

void
invsign(double **A, int m, int n, double **B)
{
  int i, j;

  for (i=1; i<=m; i++) for (j=1; j<=n; j++) B[i][j]=-A[i][j];
}

void
invsign(double **A, int m, int n)
{
  int i, j;

  for (i=1; i<=m; i++) for (j=1; j<=n; j++) A[i][j]=-A[i][j];
}

void
invsign(double *u, int m, double *v)
{
  int i;

  for (i=1; i<=m; i++) v[i]=-u[i];
}

void
invsign(double *u, int m)
{
  int i;

  for (i=1; i<=m; i++) u[i]=-u[i];
}

//////////////////////////////////////////////////////////////////////////

void diag(double **A, double *v, int n)
{
  int i;

  zeros(A,n,n);  
  for (i=1; i<=n; i++)
    A[i][i] = v[i];
}

double **diag(double *v, int n)
{
  double **A = dmatrix(1,n,1,n);
  
  diag(A,v,n);
  return(A);
}


void apply(double **A, int m, int n, double (*f)(double), double **B)
{
  int i,j;

  for ( i=1; i<=m; i++) 
    for ( j=1; j<=n; j++) B[i][j] = f(A[i][j]);
}

double **apply(double **A, int m, int n, double (*f)(double))
{
  double **B = dmatrix(1,m,1,n);
  
  apply(A,m,n,f,B);
  return(B);
}

void apply(double *u, int m, double (*f)(double), double *v)
{
  int i;

  for ( i=1; i<=m; i++) 
    v[i] = f(u[i]);
}

double *apply(double *u, int m, double (*f)(double))
{
  double *v = dvector(1,m);
  
  apply(u,m,f,v);
  return(v);
}


void sliding_mean(double **A, int m, int n, double **mA, int N)
{
  int i,j;

  for ( i=1; i<=m; i++) 
    for ( j=1; j<=n; j++) 
      mA[i][j] = (mA[i][j]*(N-1) + A[i][j])/N;
}


void sliding_mean(double *u, int m, double *mu, int N)
{
  int i;

  for ( i=1; i<=m; i++) mu[i] = (mu[i]*(N-1) + u[i])/N;
}




/* GAUSS-JORDAN ELIMINATION                                         */
/* a[1..n][1..n] is the input matrix                                */
/* b[1..n][1..m] is the set of right-hand side vectors              */
/* a is replaced by its inverse ; b by the set of solutions vectors */

void 
gaussjordan(double **a, int n, double **b, int m)
{
  int *indxc, *indxr, *ipiv;
  int i, icol, irow, j, k, l, ll;
  double big, dum, pivinv, temp;
  
  indxc=ivector(1,n);
  indxr=ivector(1,n);
  ipiv=ivector(1,n);
  for (j=1; j<=n; j++) ipiv[j]=0;
  for (i=1; i<=n; i++) {
    big=0.0;
    for (j=1; j<=n; j++)
      if (ipiv[j] != 1)
	for (k=1; k<=n; k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  } else if (ipiv[k] > 1) nrerror("gaussj : Singular Matrix-1");
	}
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=1; l<=n; l++) swap(a[irow][l],a[icol][l]);
      for (l=1; l<=m; l++) swap(b[irow][l],b[icol][l]);
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) nrerror("gaussj : Singular Matrix-2");
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=1; l<=n; l++) a[icol][l]*=pivinv;
    for (l=1; l<=m; l++) b[icol][l]*=pivinv;
    for (ll=1; ll<=n; ll++)
      if (ll != icol) {
	dum=a[ll][icol];
	a[ll][icol]=0.0;
	for (l=1; l<=n; l++) a[ll][l]-=a[icol][l]*dum;
	for (l=1; l<=m; l++) b[ll][l]-=b[icol][l]*dum;
      }
  }
  for (l=n; l>=1; l--) {
    if (indxr[l] != indxc[l])
      for (k=1; k<=n; k++) 
	swap(a[k][indxr[l]], a[k][indxc[l]]);
  }
  free_ivector(ipiv,1,n);
  free_ivector(indxr,1,n);
  free_ivector(indxc,1,n);
}


void inv(double **A, int m, double **B)
{
  double **dum=dmatrix(1,m,1,1);

  if (B!=A) {
    copy(B,A,m);
    gaussjordan(B,m,dum,1);
  }
  else gaussjordan(A,m,dum,1); // gaussj couldn't (???!) compile

  free_dmatrix(dum,1,m,1,1);
}
  
double **inv(double **A, int m)
{
  double **B = dmatrix(1,m,1,m);
  inv(A,m,B);

  return(B);
}

//*************************************************
//*************************************************
//*************************************************

/*
void matadd(double **a, double **b, double **c, int n)
{
  int i,j;

  for (j=1;j<=n;j++)
    for (i=1;i<=n;i++)
     c[i][j]=a[i][j]+b[i][j];
}


void matadd(double **a, double **b, double **c, int m, int n)
{
  int i,j;

  for (j=1;j<=n;j++)
    for (i=1;i<=m;i++)
     c[i][j]=a[i][j]+b[i][j];
}
*/
    
void vecadd(double *a, double *b, double *c, int n)
{
  int i;

  for (i=1;i<=n;i++)
    c[i]=a[i]+b[i];
}     

/*
void matsub(double **a, double **b, double **c, int n)
{
  int i,j;
  
  for (j=1;j<=n;j++)
    for (i=1;i<=n;i++)
      c[i][j]=a[i][j]-b[i][j];
}

void matsub(double **a, double **b, double **c, int m, int n)
{
  int i,j;

  for (j=1;j<=n;j++)
    for (i=1;i<=m;i++)
      c[i][j]=a[i][j]-b[i][j];
}
*/

void vecsub(double *a, double *b, double *c, int n)
{
  int i;

  for (i=1;i<=n;i++)
    c[i]=a[i]-b[i];
}     

/*
void copy(double **aout, double **ain, int n)
{
  int i,j;

  for (i=1;i<=n;i++)
    for (j=1;j<=n;j++)
      aout[j][i]=ain[j][i];
}
*/

void copy(double **aout, double **ain, int m, int n)
{
  int i,j;

  for (i=1;i<=n;i++)
    for (j=1;j<=m;j++)
      aout[j][i]=ain[j][i];
}

void copy(double **aout, double **ain, int m1, int m2, int n1, int n2)
{
  int i,j;

  for (i=n1;i<=n2;i++)
    for (j=m1;j<=m2;j++)
      aout[j][i]=ain[j-m1+1][i-n1+1];
}

void veccopy(double *aout, double *ain, int n)
{
  int i;

  for (i=1;i<=n;i++)
      aout[i]=ain[i];
}

void veccopy(double *aout, double *ain, int n1, int n2)
{
  int i;

  for (i=n1;i<=n2;i++)
      aout[i]=ain[i-n1+1];
}

void veccopy(int *aout, int *ain, int n)
{
  int i;

  for (i=1;i<=n;i++)
      aout[i]=ain[i];
}

void veccopy(int *aout, int *ain, int n1, int n2)
{
  int i;

  for (i=n1;i<=n2;i++)
      aout[i]=ain[i-n1+1];
}


//*************************************************

double **matextract(double **A, int m1, int m2, int n1, int n2)
{
  int i,j;
  double **B=dmatrix(1,m2-m1+1,1,n2-n1+1);

  for (i=m1; i<=m2; i++) 
    for (j=n1; j<=n2; j++) 
      B[i-m1+1][j-n1+1]=A[i][j];
  return(B);
}


double *vecextract(double *u, int m1, int m2)
{
  int i;
  double *v=dvector(1,m2-m1+1);

  for (i=m1; i<=m2; i++) v[i-m1+1]=u[i];
  return(v);
}
  
void matextract(double **A, int m1, int m2, int n1, int n2, double **B)
{
  int i,j;

  for (i=m1; i<=m2; i++) 
    for (j=n1; j<=n2; j++) 
      B[i-m1+1][j-n1+1]=A[i][j];
}


void vecextract(double *u, int m1, int m2, double *v)
{
  int i;

  for (i=m1; i<=m2; i++) v[i-m1+1]=u[i];
}
  

//*************************************************

void matassign(double **A, int m, int n, double **B, int m1, int n1)
{
  int i,j;

  for (i=1; i<=m; i++) 
    for (j=1; j<=n; j++) 
      B[i+m1-1][j+n1-1]=A[i][j];  
}

void vecassign(double *u, int m, double *v, int m1)
{
  int i;

  for (i=1; i<=m; i++) 
      v[i+m1-1]=u[i];  
}


//*************************************************


void scalarmult(double **A, int m, int n, double x)
{
  int i, j;

  for (i=1; i<=m; i++) 
    for (j=1; j<=n; j++) 
      A[i][j]*=x;
}

void scalarmult(double *u, int m, double x)
{
  int i;

  for (i=1; i<=m; i++) 
      u[i]*=x;
}

void scalarmult(double **A, int m, int n, double x, double **B)
{
  int i, j;

  for (i=1; i<=m; i++) 
    for (j=1; j<=n; j++) 
      B[i][j] = x*A[i][j];
}

void scalarmult(double *u, int m, double x, double *v)
{
  int i;

  for (i=1; i<=m; i++) 
      v[i] = x*u[i];
}

void mat2vec(double **A, int m, char flag, double *u)
{
  int i;
  
  switch (flag) {
  case 'c': // input == column
    for (i=1; i<=m; i++) u[i]=A[i][1];
    break;
  case 'r': // row
    for (i=1; i<=m; i++) u[i]=A[1][i];  // usually best replaced by u=A[1]
    break;
  }
}

void vec2mat(double *u, int m, char flag, double **A)
{
  int i;
  
  switch (flag) {
  case 'c': // input == column
    for (i=1; i<=m; i++) A[i][1]=u[i];
    break;
  case 'r': // row
    for (i=1; i<=m; i++) A[1][i]=u[i];  // usually best replaced by A[1]=u
    break;
  }
}


//*************************************************

double **formm(double *u, int m, int n)
{
  int i, j, k=1;
  double **A=dmatrix(1,m,1,n);

  for (j=1; j<=n; j++) 
    for (i=1; i<=m; i++) 
      A[i][j]=u[k++];   // Matlab's matrix representation

  return(A);
}


void formm(double *u, int m, int n, double **A)
{
  int i, j, k=1;

  for (j=1; j<=n; j++) 
    for (i=1; i<=m; i++) 
      A[i][j]=u[k++];   // Matlab's matrix representation
}



double **forms(double *u, int p)
{
  int i,j,k=1;
  int n=(int)rint(sqrt(2*p+0.25) - 0.5);
  double **A=dmatrix(1,n,1,n);

  for (j=1; j<=n ; j++) 
    for (i=j; i<=n; i++) 
      A[i][j]=u[k++];
  for (j=1; j<=n ; j++) 
    for (i=1; i<j; i++) 
      A[i][j]=A[j][i];  // symmetry

  return(A);
}



void forms(double *u, int p, double **A)
{
  int i,j,k=1;
  int n=(int)rint(sqrt(2*p+0.25) - 0.5);

  for (j=1; j<=n ; j++) 
    for (i=j; i<=n; i++) 
      A[i][j]=u[k++];
  for (j=1; j<=n ; j++) 
    for (i=1; i<j; i++) 
      A[i][j]=A[j][i];  // symmetry
}

//*************************************************


double **zeros(int m, int n)
{
  int i,j;
  double **A=dmatrix(1,m,1,n);
  
  for (i=1; i<=m; i++) 
    for (j=1; j<=n; j++) 
      A[i][j]=0;
  
  return(A);
}

void zeros(double **A, int m, int n)
{
  int i,j;
  
  for (i=1; i<=m; i++) 
    for (j=1; j<=n; j++) 
      A[i][j]=0;
}

double *zeros(int m)
{
  int i;
  double *u=dvector(1,m);
  
  for (i=1; i<=m; i++) 
      u[i]=0;

  return(u);
}

void zeros(double *u, int m )
{
  int i;
  
  for (i=1; i<=m; i++) 
      u[i]=0;
}


double **ones(int m, int n)
{
  int i,j;
  double **A=dmatrix(1,m,1,n);

  for (i=1; i<=m; i++) 
    for (j=1; j<=n; j++) 
      A[i][j]=1;
  
  return(A);
}

void ones(double **A, int m, int n)
{
  int i,j;
  
  for (i=1; i<=m; i++) 
    for (j=1; j<=n; j++) 
      A[i][j]=1;
}

double *ones(int m)
{
  int i;
  double *u=dvector(1,m);
  
  for (i=1; i<=m; i++) 
      u[i]=1;

  return(u);
}

void ones(double *u, int m)
{
  int i;
  
  for (i=1; i<=m; i++) 
      u[i]=1;
}

double **eye(int m)
{
  int i;
  double **A=zeros(m,m);
  
  for (i=1; i<=m; i++) 
      A[i][i]=1;

  return(A);
}

void eye(double **A, int m)
{
  int i;
  
  zeros(A,m,m);
  
  for (i=1; i<=m; i++) 
      A[i][i]=1;
}


//*************************************************


double *veclinspace(double a, double b, int N)
{
  int i;
  double *u=NULL;

  if (N<=1) {
    if (fabs(a-b)<EPSILON && N==1) {
      u[1] = a;
    }
    else
      fprintf(stderr,"Veclinspace(a,b,N): N must be >1\n");
  }
  else {
    u=dvector(1,N);
    for (i=1; i<=N; i++) u[i]=a + (b-a)*(i-1)/(N-1);
  }
  return(u);
}

double **matlinspace(double a, double b, int N, char flag)
{
  int i;
  double **A=NULL;

  if (N<=1) {
    if (fabs(a-b)<EPSILON && N==1) {
      A[1][1] = a;
    }
    else
      fprintf(stderr,"Matlinspace(a,b,N,flag): N must be >1\n");
  }
  else 
    switch (flag) {
    case 'c':
      A=dmatrix(1,N,1,1);
      for (i=1; i<=N; i++) A[i][1]=a + (b-a)*(i-1)/(N-1);
      break;
    case 'r':
      A=dmatrix(1,1,1,N);
      for (i=1; i<=N; i++) A[1][i]=a + (b-a)*(i-1)/(N-1);
      break;
    }
  return(A);
}


//*************************************************


void matdisplay(double **M, int m, int n, const char *matname)
{
  int i, j;

  fprintf(stdout,"%s (%d x %d):\n",matname,m,n);
  for (i=1; i<=m; i++) {
    for (j=1; j<=n; j++) fprintf(stdout,"%3.6g\t",M[i][j]);
    fprintf(stdout,"\n");
  }
  fprintf(stdout,"\n");
}

void vecdisplay(double *M, int m, const char *matname)
{
  int i;

  fprintf(stdout,"%s (%d):\n",matname,m);
  for (i=1; i<=m; i++) {
    fprintf(stdout,"%3.6f\t",M[i]);
    fprintf(stdout,"\n");
  }
  fprintf(stdout,"\n");
}



//*************************************************

#ifdef _WIN32 // Crappy Windows needs this
int isnan(double x)
{
  return(x!=x);
}
#endif

//*************************************************



void table1(double **x, double tindex, int m, int n, double *res)
{
  int i,j, ima, imi;
  double tm, xma, xmi;

  imi=1;
  ima=m;
  while (ima-imi>1) {
    i=imi+(ima-imi)/2;
    tm=x[i][1]; //first column of x = time
    if (tm>tindex) ima=i; else imi=i;
  }
  for (j=1; j<=n; j++) {
    xma=x[ima][j+1]; 
    xmi=x[imi][j+1];
    res[j]=xmi + (xma-xmi)*(tindex-x[imi][1])/(x[ima][1]-x[imi][1]);
  }
}


double *table1(double **x, double tindex, int m, int n)
{
  double *res=dvector(1,n);

  table1(x,tindex,m,n,res);

  return(res); // NR-type vector
}


void circtable1(double **x, int index, double tindex, int m, int n, double *res)
{
  int i,j, ima, imi;
  double tm, xma, xmi;
  double *t = dvector(1,m);

  for ( i=index+1; i<=m; i++ ) t[i-index] = x[i][1];   // unwrap time markers
  for ( i=1; i<=index; i++ ) t[m+i-index] = x[i][1];

  //  for ( i=1; i<=m; i++ ) t[i] = x[(i-index+m-1)%m +1][i]; // unwrap time markers (shorter to write/longer to execute!)

  imi=1;
  ima=m;
  while (ima-imi>1) {
    i=imi+(ima-imi)/2;
    tm=t[i];
    if (tm>tindex) ima=i; else imi=i;
  }
  imi = (imi+index-1)%m + 1;
  ima = (ima+index-1)%m + 1;
  for (j=1; j<=n; j++) {
    xma=x[ima][j+1]; 
    xmi=x[imi][j+1];
    res[j]=xmi + (xma-xmi)*(tindex-x[imi][1])/(x[ima][1]-x[imi][1]);
  }

  free_dvector(t,1,m);
}


double *circtable1(double **x, int index, double tindex, int m, int n)
{
  double *res=dvector(1,n);

  circtable1(x,index,tindex,m,n,res);

  return(res); // NR-type vector
}

//*************************************************

void mat2col(double **smat,  int m, int n, double *dmat)
{
  int i, j;

  for (j=1; j<=n ; j++) 
    for (i=1; i<=m ; i++)   
      dmat[i+m*(j-1)]=smat[i][j];  // Fortran representation
}

void col2mat(double *smat, int m, int n, double **dmat)
{
  int i, j;

  for (j=1; j<=n ; j++) 
    for (i=1; i<=m ; i++)         
      dmat[i][j]=smat[i+m*(j-1)];  // Fortran representation
}


// version with useful dimension
void col2symmat2(double *smat, int n, double **dmat)
{
  int i, j, k=1;

  for (j=1; j<=n ; j++) 
    for (i=j; i<=n; i++) 
      dmat[j][i]=smat[k++];
  for (j=1; j<=n ; j++) 
    for (i=1; i<j; i++) 
      dmat[j][i]=dmat[i][j];  // symmetry
}


// OLD VERSION

void col2symmat(double *smat, int m, double **dmat)
{
  int i, j, n, k=1;

  n=(int)rint(sqrt(2*m+0.25) - 0.5);

  for (j=1; j<=n ; j++) 
    for (i=j; i<=n; i++) 
      dmat[j][i]=smat[k++];
  for (j=1; j<=n ; j++) 
    for (i=1; i<j; i++) 
      dmat[j][i]=dmat[i][j];  // symmetry
}


void symmat2col(double **smat, int m, double *dmat)
{
  int i, j, k=1;

  for (j=1; j<=m ; j++) 
    for (i=j; i<=m; i++) 
      dmat[k++]=smat[j][i];  
}


void colconcat(double *a, int n1, double *b, int n2)
{
  int i;

  for (i=1; i<=n2; i++) a[n1+i]=b[i];  // WARNING: proper allocation needed
}


// THIS IS A PROPER IMPLEMENTATION
// version with allocation
// default behavior is horiz concat (cf. mylinpack.h)
double **matconcat(double **a, int m1, int n1, double **b, int m2, int n2, char dir)
{
  int i, j;
  double **c=NULL;

  if ( dir == 'h') { // horiz concat
    if (m1==m2) {
      c=dmatrix(1,m1,1,n1+n2);
      copy(c,a,m1,n1);
      for (i=1; i<=m1; i++)
	for (j=n1+1; j<=n1+n2; j++)
	  c[i][j]=b[i][j-n1];  
    }
    else
      fprintf(stderr,"Error in matconcat: Matrices of improper size\n");    
  }
  else {  // vert concat
    if (n1==n2) {
      c=dmatrix(1,m1+m2,1,n1);
      copy(c,a,m1,n1);
      for (j=1; j<=n1; j++)
	for (i=m1+1; i<=m1+m2; i++)
	  c[i][j]=b[i-m1][j];
    }
    else
      fprintf(stderr,"Error in matconcat: Matrices of improper size\n");
  }
  return(c);
}

// without allocation
// default is horiz concat
void matconcat(double **a, int m1, int n1, double **b, int m2, int n2, double **c, char dir)
{
  int i, j;

  if ( dir == 'h') { // horiz concat
    if (m1==m2) {
      copy(c,a,m1,n1);
      for (i=1; i<=m1; i++)
	for (j=n1+1; j<=n1+n2; j++)
	  c[i][j]=b[i][j-n1];  
    }
    else 
      fprintf(stderr,"Error in matconcat: Matrices of improper size\n");    
  }
  else {  // vert concat
    if (n1==n2) {
      copy(c,a,m1,n1);
      for (j=1; j<=n1; j++)
	for (i=m1+1; i<=m1+m2; i++)
	  c[i][j]=b[i-m1][j];
    }
    else
      fprintf(stderr,"Error in matconcat: Matrices of improper size\n");
  }
}

void diagconcat(double **A, int m1, int n1, double **B, int m2, int n2, double **C)
{
  int i, j;
  double **Z1 = zeros(m1,n2);
  double **Z2 = zeros(m2,n1);
  double **temp1, **temp2;

  temp1 = matconcat(A,m1,n1,Z1,m1,n2);
  temp2 = matconcat(Z2,m2,n1,B,m2,n2);
  matconcat(temp1,m1,n1+n2, temp2,m2,n1+n2, C, 'v');

  free_dmatrix(Z1,1,m1,1,n2);
  free_dmatrix(Z2,1,m2,1,n1);
  free_dmatrix(temp1,1,m1,1,n1+n2);
  free_dmatrix(temp2,1,m2,1,n1+n2);
}


double **diagconcat(double **A, int m1, int n1, double **B, int m2, int n2)
{
  double **C = dmatrix(1,m1+m2,1,n1+n2);

  diagconcat(A,m1,n1,B,m2,n2,C);

  return(C);
}


// without alloc
void vecconcat(double *a, int m1, double *b, int m2, double *c)
{
  int i;

  veccopy(c,a,m1);
  for (i=m1+1; i<=m1+m2; i++)
    c[i]=b[i-m1]; 
}

// with alloc
double *vecconcat(double *a, int m1, double *b, int m2)
{
  double *c=NULL;

  c=dvector(1,m1+m2);
  vecconcat(a,m1,b,m2,c);
  return(c);
}

//*************************************************

double mean(double *u, int n)
{
  int i;
  double mu = 0;
  
  for ( i=1; i<=n; i++ )
    mu += u[i];
  
  return(mu/n);
}

void covmat(double *u, double **A, int n)
{
  int i, j;
  double meanu;
  double *v = dvector(1,n);

  meanu = mean(u,n);

  for ( i=1; i<=n; i++ ) v[i] = u[i] - meanu;

  todprod(v,v,n,A);
}


double **covmat(double *u, int n)
{
  double **A = dmatrix(1,n,1,n);
  
  covmat(u,A,n);

  return(A);
}





//*************************************************





void treuler(double ystart[], int nvar, double x1, double x2,
        double eps, int Nsteps, void (*derivs)(double, double [], double []),int *pkount, double *xp, double **yp)
{
  int i, c;
  double x, xstep;
  double *y=dvector(1,nvar);
  double *dy=dvector(1,nvar);
  double *dydx=dvector(1,nvar);

  *pkount=Nsteps;  // for compatibility reasons with odeint
  xstep=fabs(x1-x2)/(Nsteps-1);
  for (i=1; i<=nvar; i++) yp[i][1]=ystart[i];

  if (x1<x2) {
    for (x=x1, c=1; c<Nsteps; x+=xstep, c++) {
      for (i=1; i<=nvar; i++) y[i]=yp[i][c];
      (*derivs)(x,y,dydx);
      for (i=1; i<=nvar; i++) {
	dy[i]=dydx[i]*xstep;
	yp[i][c+1]=y[i] + dy[i];
	//	rdy[i]=dy[i]/y[i];
      }
      xp[c]=x;
      //      if (norm2(rdy,nvar)>eps) xstep/=2;
    }
    xp[c]=x2;
  }
  else {
    for (x=x1, c=1; c<Nsteps; x-=xstep, c++) {
      for (i=1; i<=nvar; i++) y[i]=yp[i][c];
      (*derivs)(x,y,dydx);
      for (i=1; i<=nvar; i++) {
	dy[i]=-dydx[i]*xstep;
	yp[i][c+1]=y[i] + dy[i];
	//	rdy[i]=dy[i]/y[i];
      }
      xp[c]=x;
      //      if (norm2(rdy,nvar)>eps) xstep/=2;
    }
    xp[c]=x2;
  }
  free_dvector(y,1,nvar);
  free_dvector(dy,1,nvar);
  free_dvector(dydx,1,nvar);

  //  for (c=1; c<=Nsteps; c++) 
  //    printf("%g\t%g\t%g\n",xp[c],yp[2][c],yp[3][c]);
}


void backwards_euler(double ystart[], int nvar, double x1, double x2,
        double eps, int Nsteps, void (*derivs)(double, double [], double []),int *pkount, double *xp, double **yp)
{
  int i, c;
  double x, xstep;
  double *y=dvector(1,nvar);
  double *yguess=dvector(1,nvar);
  double *ydiff=dvector(1,nvar);
  double *dy=dvector(1,nvar);
  double *dydx=dvector(1,nvar);

  *pkount=Nsteps;
  xstep=fabs(x1-x2)/(Nsteps-1);
  veccopy(yp[1],ystart,nvar);

  if (x1<x2) {
    for (x=x1, c=1; c<Nsteps; x+=xstep, c++) {
      (*derivs)(x,yp[c],dydx);
      veccopy(y,yp[c],nvar);
      do {
	for (i=1; i<=nvar; i++) yguess[i] = yp[c][i] + dydx[i]*xstep;
	(*derivs)(x,yguess,dydx);
	vecsub(yguess,y,ydiff,nvar);
	veccopy(y,yguess,nvar);
      }
      while(norm2(ydiff,nvar)>eps);
      for (i=1; i<=nvar; i++) {
	dy[i]=dydx[i]*xstep;
	yp[c+1][i]=yguess[i];
	//	rdy[i]=dy[i]/y[i];
      }
      xp[c]=x;
      //      if (norm2(rdy,nvar)>eps) xstep/=2;
    }
    xp[c]=x2;
  }
  else {
    for (x=x1, c=1; c<Nsteps; x-=xstep, c++) {
      (*derivs)(x,yp[c],dydx);
      veccopy(y,yp[c],nvar);
      do {
	for (i=1; i<=nvar; i++) yguess[i] = yp[c][i] - dydx[i]*xstep;
	(*derivs)(x,yguess,dydx);
	vecsub(yguess,y,ydiff,nvar);
	veccopy(y,yguess,nvar);
      }
      while(norm2(ydiff,nvar)>eps);
      for (i=1; i<=nvar; i++) {
	dy[i]=-dydx[i]*xstep;
	yp[c+1][i]=yguess[i];
	//	rdy[i]=dy[i]/y[i];
      }
      xp[c]=x;
      //      if (norm2(rdy,nvar)>eps) xstep/=2;
    }
    xp[c]=x2;
  }
  free_dvector(y,1,nvar);
  free_dvector(yguess,1,nvar);
  free_dvector(ydiff,1,nvar);
  free_dvector(dy,1,nvar);
  free_dvector(dydx,1,nvar);
}



void euler(double ystart[], int nvar, double x1, double x2,
        double eps, int Nsteps, void (*derivs)(double, double [], double []),int *pkount, double *xp, double **yp)
{
  int i, c;
  double x, xstep;
  double *y=dvector(1,nvar);
  double *dy=dvector(1,nvar);
  double *dydx=dvector(1,nvar);

  *pkount=Nsteps;
  xstep=fabs(x1-x2)/(Nsteps-1);
  veccopy(yp[1],ystart,nvar);

  if (x1<x2) {
    for (x=x1, c=1; c<Nsteps; x+=xstep, c++) {
      (*derivs)(x,yp[c],dydx);
      for (i=1; i<=nvar; i++) {
	dy[i]=dydx[i]*xstep;
	yp[c+1][i]=yp[c][i] + dy[i];
	//	rdy[i]=dy[i]/y[i];
      }
      xp[c]=x;
      //      if (norm2(rdy,nvar)>eps) xstep/=2;
    }
    xp[c]=x2;
  }
  else {
    for (x=x1, c=1; c<Nsteps; x-=xstep, c++) {
      (*derivs)(x,yp[c],dydx);
      for (i=1; i<=nvar; i++) {
	dy[i]=-dydx[i]*xstep;
	yp[c+1][i]=yp[c][i] + dy[i];
	//	rdy[i]=dy[i]/y[i];
      }
      xp[c]=x;
      //      if (norm2(rdy,nvar)>eps) xstep/=2;
    }
    xp[c]=x2;
  }
  free_dvector(y,1,nvar);
  free_dvector(dy,1,nvar);
  free_dvector(dydx,1,nvar);
}


#define MAXSTP 10000
#define TINY 1.0e-30


// personal version of odeint
void rungekutta4(double ystart[], int nvar, double x1, double x2, double eps, double h1,

double hmin, int *nok, int *nbad,

void (*derivs)(double, double [], double []),

void (*rkqs)(double [], double [], int, double *, double, double, double [],

	     double *, double *, void (*)(double, double [], double [])), 
	 int kmax, int *pkount, double *xp, double **yp)
{

  int nstp,i,kount;
  double dxsav,xsav,x,hnext,hdid,h;
  double *yscal,*y,*dydx;
  

  yscal=dvector(1,nvar);
  y=dvector(1,nvar);
  dydx=dvector(1,nvar);
  x=x1;
  h=SIGN(h1,x2-x1);

  *nok = (*nbad) = kount = 0;
  dxsav = h1;  // my choice... good?

  for (i=1;i<=nvar;i++) y[i]=ystart[i];
  if (kmax > 0) xsav=x-dxsav*2.0;
  
  
  for (nstp=1;nstp<=MAXSTP;nstp++) {
    (*derivs)(x,y,dydx);
    for (i=1;i<=nvar;i++)
      yscal[i]=fabs(dydx[i]*h)+TINY;
    //      yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
    if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
      xp[++kount]=x;
      for (i=1;i<=nvar;i++) yp[kount][i]=y[i];
      xsav=x;
    }
    if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
    (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
    if (hdid == h) ++(*nok); else ++(*nbad);
    if ((x-x2)*(x2-x1) >= 0.0) {
      for (i=1;i<=nvar;i++) ystart[i]=y[i];
      if (kmax) {
	xp[++kount]=x;
	for (i=1;i<=nvar;i++) yp[kount][i]=y[i];
      }
      free_dvector(dydx,1,nvar);
      free_dvector(y,1,nvar);
      free_dvector(yscal,1,nvar);
      *pkount = kount;
      return;
    }
    if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
    h=hnext;
  }
  
  nrerror("Too many steps in routine odeint");
}
#undef MAXSTP
#undef TINY



//*************************************************
//*************************************************
//*************************************************
//*************************************************

//  FLOAT VERSIONS 

//*************************************************
//*************************************************
//*************************************************
//*************************************************

/*

// MATRIX MULTIPLICATION 

void 
fmatmult(int m, float **A, int n,  float **B, int p, float**C)
{
  int i,j,k;
  float S;

  for (i=1; i<=m; i++) 
    for (j=1; j<=p; j++) {
      S=0;
      for (k=1; k<=n; k++) S+=A[i][k]*B[k][j];
      C[i][j]=S;
    }
}

// MATRIX TRANSPOSITION 

void 
ftranspose(float **A, int m, int n, float **B)
{
  int i,j;

  for (i=1; i<=m; i++)
    for (j=1; j<=n; j++) B[j][i]=A[i][j];
}

// VECTOR SUBTRACTION 

void 
fvectsub(float *u, float *v, int n, float *w)
{
  int i;

  for(i=1; i<=n; i++) w[i]=u[i]-v[i];
}

// SQUARED EUCLIDIAN NORM 

float 
fsquared(float *v, int n)
{
  int i;
  float sq=0;

  for (i=1; i<=n; i++) sq+=v[i]*v[i];
  return(sq);
}


// EUCLIDIAN NORM 

float 
fnorm2(float *v, int n)
{
  return(sqrt(fsquared(v,n)));
}

float 
fnorm2(float **A, int m, int n)
{
  int i;
  float sq=0;

  for (i=1; i<=m; i++) sq+=fsquared(A[i],n);
  return(sqrt(sq));
}

// GENERIC p-NORM (slower) 

float 
fnorm(float *v, int n, int p)
{
  int i;
  float norm=0;

  for (i=1; i<=n; i++) norm+=pow(fabs(v[i]),p); 
  return(pow(norm,1/p));
  
}

// CROSS PRODUCT 

void
fcrossprod(float* u, float* v, float* w)
{
  int i;

  for (i=1; i<=3; i++) w[i]=u[i%3 + 1]*v[(i+1)%3+1] - v[i%3 + 1]*u[(i+1)%3+1];
}

// SIGN INVERSION 

void
finvsign(float **A, int m, int n, float **B)
{
  int i, j;

  for (i=1; i<=m; i++) for (j=1; j<=n; j++) B[i][j]=-A[i][j];
}

// GAUSS-JORDAN ELIMINATION                                         
// a[1..n][1..n] is the input matrix                                
// b[1..n][1..m] is the set of right-hand side vectors              
// a is replaced by its inverse ; b by the set of solutions vectors 

void 
fgaussjordan(float **a, int n, float **b, int m)
{
  int *indxc, *indxr, *ipiv;
  int i, icol, irow, j, k, l, ll;
  float big, dum, pivinv, temp;
  
  indxc=ivector(1,n);
  indxr=ivector(1,n);
  ipiv=ivector(1,n);
  for (j=1; j<=n; j++) ipiv[j]=0;
  for (i=1; i<=n; i++) {
    big=0.0;
    for (j=1; j<=n; j++)
      if (ipiv[j] != 1)
	for (k=1; k<=n; k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  } else if (ipiv[k] > 1) nrerror("gaussj : Singular Matrix-1");
	}
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=1; l<=n; l++) swap(a[irow][l],a[icol][l]);
      for (l=1; l<=m; l++) swap(b[irow][l],b[icol][l]);
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) nrerror("gaussj : Singular Matrix-2");
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=1; l<=n; l++) a[icol][l]*=pivinv;
    for (l=1; l<=m; l++) b[icol][l]*=pivinv;
    for (ll=1; ll<=n; ll++)
      if (ll != icol) {
	dum=a[ll][icol];
	a[ll][icol]=0.0;
	for (l=1; l<=n; l++) a[ll][l]-=a[icol][l]*dum;
	for (l=1; l<=m; l++) b[ll][l]-=b[icol][l]*dum;
      }
  }
  for (l=n; l>=1; l--) {
    if (indxr[l] != indxc[l])
      for (k=1; k<=n; k++) 
	swap(a[k][indxr[l]], a[k][indxc[l]]);
  }
  free_ivector(ipiv,1,n);
  free_ivector(indxr,1,n);
  free_ivector(indxc,1,n);
}
  
    
void fmatadd(float **a, float **b, float **c, int n)
{
  int i,j;

  for (j=1;j<=n;j++)
    for (i=1;i<=n;i++)
     c[i][j]=a[i][j]+b[i][j];
}

void fmatadd(float **a, float **b, float **c, int m, int n)
{
  int i,j;

  for (j=1;j<=n;j++)
    for (i=1;i<=m;i++)
     c[i][j]=a[i][j]+b[i][j];
}

    
void fvecadd(float *a, float *b, float *c, int n)
{
  int i;

  for (i=1;i<=n;i++)
    c[i]=a[i]+b[i];
}     


void fmatsub(float **a, float **b, float **c, int n)
{
  int i,j;
  
  for (j=1;j<=n;j++)
    for (i=1;i<=n;i++)
      c[i][j]=a[i][j]-b[i][j];
}


void fmatsub(float **a, float **b, float **c, int m, int n)
{
  int i,j;

  for (j=1;j<=n;j++)
    for (i=1;i<=m;i++)
      c[i][j]=a[i][j]-b[i][j];
}

void fcopy(float **aout, float **ain, int n)
{
  int i,j;

  for (i=1;i<=n;i++)
    for (j=1;j<=n;j++)
      aout[j][i]=ain[j][i];
}

void fcopy(float **aout, float **ain, int m, int n)
{
  int i,j;

  for (i=1;i<=n;i++)
    for (j=1;j<=m;j++)
      aout[j][i]=ain[j][i];
}

void fveccopy(float *aout, float *ain, int n)
{
  int i;

  for (i=1;i<=n;i++)
      aout[i]=ain[i];
}




void matrix2dmatrix(float **src, double **dest, int m, int n)
{
  int i,j;

  for (i=1; i<=m; i++)
    for (j=1; j<=n; j++)
      dest[i][j]=src[i][j];
}


void dmatrix2matrix(double **src, float **dest, int m, int n)
{
  int i,j;

  for (i=1; i<=m; i++)
    for (j=1; j<=n; j++)
      dest[i][j]=(float)src[i][j];
}


float *d2f(double *src, int l)
{
  int i;
  float *res=(float*)malloc(l*sizeof(float));

  for (i=0; i<l; i++) res[i]=(float)src[i];
  return(res);
}



//*************************************************

float **fmatextract(float **A, int m1, int m2, int n1, int n2)
{
  int i,j;
  float **B=matrix(1,m2-m1+1,1,n2-n1+1);

  for (i=m1; i<=m2; i++) 
    for (j=n1; j<=n2; j++) 
      B[i-m1+1][j-n1+1]=A[i][j];
  return(B);
}


float *fvecextract(float *u, int m1, int m2)
{
  int i;
  float *v=vector(1,m2-m1+1);

  for (i=m1; i<=m2; i++) v[i-m1+1]=u[i];
  return(v);
}
  
void fmatextract(float **A, int m1, int m2, int n1, int n2, float **B)
{
  int i,j;

  for (i=m1; i<=m2; i++) 
    for (j=n1; j<=n2; j++) 
      B[i-m1+1][j-n1+1]=A[i][j];
}


void fvecextract(float *u, int m1, int m2, float *v)
{
  int i;

  for (i=m1; i<=m2; i++) v[i-m1+1]=u[i];
}
  

//*************************************************


void fscalarmult(float **A, int m, int n, float x)
{
  int i, j;

  for (i=1; i<=m; i++) 
    for (j=1; j<=n; j++) 
      A[i][j]*=x;
}

void fscalarmult(float *u, int m, float x)
{
  int i;

  for (i=1; i<=m; i++) 
      u[i]*=x;
}

void fmat2vec(float **A, int m, char flag, float *u)
{
  int i;
  
  switch (flag) {
  case 'c': // input == column
    for (i=1; i<=m; i++) u[i]=A[i][1];
    break;
  case 'r': // row
    for (i=1; i<=m; i++) u[i]=A[1][i];  // usually best replaced by u=A[1]
    break;
  }
}

void fvec2mat(float *u, int m, char flag, float **A)
{
  int i;
  
  switch (flag) {
  case 'c': // input == column
    for (i=1; i<=m; i++) A[i][1]=u[i];
    break;
  case 'r': // row
    for (i=1; i<=m; i++) A[1][i]=u[i];  // usually best replaced by A[1]=u
    break;
  }
}


//*************************************************

float **fformm(float *u, int m, int n)
{
  int i, j, k=1;
  float **A=matrix(1,m,1,n);

  for (j=1; j<=n; j++) 
    for (i=1; i<=m; i++) 
      A[i][j]=u[k++];   // Matlab's matrix representation

  return(A);
}


void fformm(float *u, int m, int n, float **A)
{
  int i, j, k=1;

  for (j=1; j<=n; j++) 
    for (i=1; i<=m; i++) 
      A[i][j]=u[k++];   // Matlab's matrix representation
}



float **fforms(float *u, int p)
{
  int i,j,k=1;
  int n=(int)rint(sqrt(2*p+0.25) - 0.5);
  float **A=matrix(1,n,1,n);

  for (j=1; j<=n ; j++) 
    for (i=j; i<=n; i++) 
      A[i][j]=u[k++];
  for (j=1; j<=n ; j++) 
    for (i=1; i<j; i++) 
      A[i][j]=A[j][i];  // symmetry

  return(A);
}



void fforms(float *u, int p, float **A)
{
  int i,j,k=1;
  int n=(int)rint(sqrt(2*p+0.25) - 0.5);

  for (j=1; j<=n ; j++) 
    for (i=j; i<=n; i++) 
      A[i][j]=u[k++];
  for (j=1; j<=n ; j++) 
    for (i=1; i<j; i++) 
      A[i][j]=A[j][i];  // symmetry
}

//*************************************************


float **fzeros(int m, int n)
{
  int i,j;
  float **A=matrix(1,m,1,n);
  
  for (i=1; i<=m; i++) 
    for (j=1; j<=n; j++) 
      A[i][j]=0;
  
  return(A);
}

void fzeros(int m, int n, float **A)
{
  int i,j;
  
  for (i=1; i<=m; i++) 
    for (j=1; j<=n; j++) 
      A[i][j]=0;
}

float *fzeros(int m)
{
  int i;
  float *u=vector(1,m);
  
  for (i=1; i<=m; i++) 
      u[i]=0;

  return(u);
}

void fzeros(int m, float *u)
{
  int i;
  
  for (i=1; i<=m; i++) 
      u[i]=0;
}

float **fones(int m, int n)
{
  int i,j;
  float **A=matrix(1,m,1,n);

  for (i=1; i<=m; i++) 
    for (j=1; j<=n; j++) 
      A[i][j]=1;
  
  return(A);
}

void fones(int m, int n, float **A)
{
  int i,j;
  
  for (i=1; i<=m; i++) 
    for (j=1; j<=n; j++) 
      A[i][j]=1;
}

float *fones(int m)
{
  int i;
  float *u=vector(1,m);
  
  for (i=1; i<=m; i++) 
      u[i]=1;

  return(u);
}

void fones(int m, float *u)
{
  int i;
  
  for (i=1; i<=m; i++) 
      u[i]=1;
}


//*************************************************


float *fveclinspace(float a, float b, int N)
{
  int i;
  float *u=NULL;

  if (N<=1) fprintf(stderr,"Veclinspace(a,b,N): N must be >1\n");
  else {
    u=vector(1,N);
    for (i=1; i<=N; i++) u[i]=a + (b-a)*(i-1)/(N-1);
  }
  return(u);
}

float **fmatlinspace(float a, float b, int N, char flag)
{
  int i;
  float **A=NULL;

  if (N<=1) 
    fprintf(stderr,"Matlinspace(a,b,N,flag): N must be >1\n");
  else 
    switch (flag) {
    case 'c':
      A=matrix(1,N,1,1);
      for (i=1; i<=N; i++) A[i][1]=a + (b-a)*(i-1)/(N-1);
      break;
    case 'r':
      A=matrix(1,1,1,N);
      for (i=1; i<=N; i++) A[1][i]=a + (b-a)*(i-1)/(N-1);
      break;
    }
  return(A);
}


//*************************************************

void fmatdisplay(float **M, int m, int n, const char *matname)
{
  int i, j;

  fprintf(stdout,"%s (%d x %d):\n",matname,m,n);
  for (i=1; i<=m; i++) {
    for (j=1; j<=n; j++) fprintf(stdout,"%3.6f\t",M[i][j]);
    fprintf(stdout,"\n");
  }
  fprintf(stdout,"\n");
}

void fvecdisplay(float *M, int m, const char *matname)
{
  int i;

  fprintf(stdout,"%s (%d):\n",matname,m);
  for (i=1; i<=m; i++) {
    fprintf(stdout,"%3.6f\t",M[i]);
    fprintf(stdout,"\n");
  }
  fprintf(stdout,"\n");
}



//*************************************************


int fisnan(float x)
{
  return(x!=x);
}





*/

