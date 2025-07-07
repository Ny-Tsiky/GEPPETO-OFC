/* EMBRYONIC LINEAR ALGEBRA PACKAGE DEFS -- P.B. 27/10/96 */
/* from Numerical Recipes */

//#include <curses.h>

//#define swap(a,b) {temp=(a); (a)=(b); (b)=temp;}
//#define sqr(a) ((a)*(a))
//#define max(a,b) ((a>b)?a:b)
//#define min(a,b) ((a<b)?a:b)



#define EPSILON 1e-10


inline void swap(double a, double b) {double temp=a; a=b; b=temp;}
inline void swap(int a, int b) {int temp=a; a=b; b=temp;}
inline int sqr(int a) {return(a*a);}
inline double sqr(double a) {return(a*a);}
inline double max(double a, double b) {return( (a > b) ? a : b);}
inline double min(double a, double b) {return( (a < b) ? a : b);}
inline double pospart(double x) {return((x>0) ? x : 0);}
inline int sign(double x) {return( (x > EPSILON) ? 1 : (x < -EPSILON) ? -1 : 0);}


#ifndef Pi
#define Pi                   3.14159265358979323846264
#endif
const double sqrt2 = sqrt(2);


#define ESC 27
#define DB {fprintf(stdout,"Got to line %d in file %s\n",__LINE__,__FILE__); fflush(stdout);}
#define DBK {fprintf(stdout,"Got to line %d in file %s (ESC to stop)\n",__LINE__,__FILE__); fflush(stdout); if (getchar()==ESC) exit(0);}


#ifdef WIN32
double rint( double x);
#endif

double erf(double x);
double erfc(double x);
double randn();
double randn2();
void randn(double *v, int n);
double *randn(int n);
void randn(double **A, int m, int n);
double **randn(int m, int n);
void rand(double **A, int m, int n);
double **rand(int m, int n);

//double sign(double x);

void matmult(int m, double **A, int n,  double **B, int p, double**C);
void vecmatmult(double *u, int n,  double **B, int p, double *v);
void matvecmult(int m, double **A, int n,  double *u, double *v);
void transpose(double **A, int m, int n, double **B);
void vectsub(double *u, double *v, int n, double *w);
double dotprod(double *u, double *v, int n);
void todprod(double *u, double *v, int n, double **A);
double **todprod(double *u, double *v, int n);
double norm(double *v, int n, int p);
void crossprod(double *u, double *v, double* w);
double squared(double *v, int n);
double norm2(double *v, int n);     /* faster */
void invsign(double **A, int m, int n, double **B);
void invsign(double **A, int m, int n);
void invsign(double *u, int m, double *v);
void invsign(double *u, int m);
void diag(double **A, double *v, int n);
double **diag(double *v, int n);

void apply(double **A, int m, int n, double (*f)(double), double **B);
double **apply(double **A, int m, int n, double (*f)(double));
void apply(double *u, int m, double (*f)(double), double *v);
double *apply(double *u, int m, double (*f)(double));

void sliding_mean(double **A, int m, int n, double **mA, int N);
void sliding_mean(double *u, int m, double *mu, int N);

void gaussjordan(double **a, int n, double **b, int m);
void inv(double **A, int m, double **B);
double **inv(double **A, int m);

    
void matadd(double **a, double **b, double **c, int n);
void matadd(double **a, double **b, double **c, int m, int n);
void vecadd(double *a, double *b, double *c, int n); 
void matsub(double **a, double **b, double **c, int n);
void matsub(double **a, double **b, double **c, int m, int n);
void vecsub(double *a, double *b, double *c, int n); 
void copy(double **aout, double **ain, int n);
void copy(double **aout, double **ain, int m, int n);
void copy(double **aout, double **ain, int m1, int m2, int n1, int n2);
void veccopy(double *aout, double *ain, int n); 
void veccopy(double *aout, double *ain, int n1, int n2); 

double **matextract(double **A, int m1, int m2, int n1, int n2);
double *vecextract(double *u, int m1, int m2);
void matextract(double **A, int m1, int m2, int n1, int n2, double **B);
void vecextract(double *u, int m1, int m2, double *v);
void matassign(double **A, int m, int n, double **B, int m1, int n1);
void vecassign(double *u, int m, double *v, int m1);


void scalarmult(double **A, int m, int n, double x);
void scalarmult(double *u, int m, double x);
void scalarmult(double **A, int m, int n, double x, double **B);
     void scalarmult(double *u, int m, double x, double *v);
void mat2vec(double **A, int m, char flag, double *u);
void vec2mat(double *u, int m, char flag, double **A);
double **formm(double *u, int m, int n);
void formm(double *u, int m, int n, double **A);
double **forms(double *u, int p);
void forms(double *u, int p, double **A);

double **zeros(int m, int n);
void zeros(double **A, int m, int n);
double *zeros(int m);
void zeros(double *u, int m);
double **ones(int m, int n);
void ones(double **A, int m, int n);
double *ones(int m);
void ones(double *u, int m);
double **eye(int m);
void eye(double **A, int m);

double *veclinspace(double a, double b, int N);
double **matlinspace(double a, double b, int N, char flag);

void matdisplay(double **M, int m, int n, const char *matname);
void matdisplay(double **M, int m, int n, const char *matname);
void vecdisplay(double *M, int m, const char *matname);
void vecdisplay(double *M, int m, const char *matname);

//int isnan(double x);

void table1(double **x, double t, int m, int n, double *res);
double *table1(double **x, double t, int m, int n);
void circtable1(double **x, int circ_index, double t, int m, int n, double *res);
double *circtable1(double **x, int circ_index, double t, int m, int n);
void mat2col(double **smat,  int m, int n, double *dmat);
void col2mat(double *smat, int m, int n, double **dmat);
void col2symmat(double *smat, int m, double **dmat);
void col2symmat2(double *smat, int m, double **dmat);
void symmat2col(double **smat, int m, double *dmat);
void colconcat(double *a, int n1, double *b, int n2);
double **matconcat(double **a, int m1, int n1, double **b, int m2, int n2, char dir = 'h');
void matconcat(double **a, int m1, int n1, double **b, int m2, int n2, double **c, char dir = 'h');
void diagconcat(double **A, int m1, int n1, double **B, int m2, int n2, double **C);
double **diagconcat(double **A, int m1, int n1, double **B, int m2, int n2);
double *vecconcat(double *a, int m1, double *b, int m2);
void vecconcat(double *a, int m1, double *b, int m2, double *c);

double mean(double *u, int n);
void covmat(double *u, double **A, int n);
double **covmat(double *u, int n);

void euler(double ystart[], int nvar, double x1, double x2, double eps, int Nssteps, void (*derivs)(double, double [], double []),int *pkount, double *xp, double **yp);
void backwards_euler(double ystart[], int nvar, double x1, double x2, double eps, int Nssteps, void (*derivs)(double, double [], double []),int *pkount, double *xp, double **yp);
void treuler(double ystart[], int nvar, double x1, double x2, double eps, int Nsteps, void (*derivs)(double, double [], double []),int *pkount, double *xp, double **yp);
void rungekutta4(double ystart[], int nvar, double x1, double x2, double eps, double h1,double hmin, int *nok, int *nbad,
void (*derivs)(double, double [], double []),
void (*rkqs)(double [], double [], int, double *, double, double, double [],
	     double *, double *, void (*)(double, double [], double [])), 
int kmax, int *pkount, double *xp, double **yp);



// FLOAT VERSIONS

void fmatmult(int m, float **A, int n,  float **B, int p, float**C);
void ftranspose(float **A, int m, int n, float **B);
void fvectsub(float *u, float *v, int n, float *w);
float fnorm(float *v, int n, int p);
void fcrossprod(float *u, float *v, float* w);
float fsquared(float *v, int n);
float fnorm2(float *v, int n);     /* faster */
float fnorm2(float **A, int m, int n);     /* matrix version */
void finvsign(float **A, int m, int n, float **B);
void fgaussjordan(float **a, int n, float **b, int m);

    
void fmatadd(float **a, float **b, float **c, int n);
void fmatadd(float **a, float **b, float **c, int m, int n);
void fvecadd(float *a, float *b, float *c, int n); // NO double VERSION
void fmatsub(float **a, float **b, float **c, int n);
void fmatsub(float **a, float **b, float **c, int m, int n);
void fcopy(float **aout, float **ain, int n);
void fcopy(float **aout, float **ain, int m, int n);
void fveccopy(float *aout, float *ain, int n); // NO double VERSION


void matrix2dmatrix(float **src, double **dest, int m, int n);
void dmatrix2matrix(double **src, float **dest, int m, int n);
float *d2f(double *src, int l);


float **fmatextract(float **A, int m1, int n1, int m2, int n2);
float *fvecextract(float *u, int m1, int m2);
void fmatextract(float **A, int m1, int n1, int m2, int n2, float **B);
void fvecextract(float *u, int m1, int m2, float *v);

void fscalarmult(float **A, int m, int n, float x);
void fscalarmult(float *u, int m, float x);
void fmat2vec(float **A, int m, char flag, float *u);
void fvec2mat(float *u, int m, char flag, float **A);
float **fformm(float *u, int m, int n);
void fformm(float *u, int m, int n, float **A);
float **fforms(float *u, int p);
void fforms(float *u, int p, float **A);

float **fzeros(int m, int n);
double **zeros(int m, int n);
void fzeros(int m, int n, float **A);
float *fzeros(int m);
void fzeros(int m, float *u);
float **fones(int m, int n);
double **ones(int m, int n);
void fones(int m, int n, float **A);
float *fones(int m);
void fones(int m, float *u);

float *fveclinspace(float a, float b, int N);
float **fmatlinspace(float a, float b, int N, char flag);

void fmatdisplay(double **M, int m, int n, const char *matname);
void fmatdisplay(float **M, int m, int n, const char *matname);
void fvecdisplay(double *M, int m, const char *matname);
void fvecdisplay(float *M, int m, const char *matname);

int fisnan(float x);
