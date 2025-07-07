#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

#define NR_END 1


#define nrerror2                 NR_nrerror2
#define nrerror                  NR_nrerror
#define vector                   NR_dvector
#define fvector                  NR_fvector
#define ivector                  NR_ivector
#define dvector                  NR_dvector
#define cvector                  NR_cvector
#define lvector                  NR_lvector
#define matrix                   NR_dmatrix
#define fmatrix                  NR_fmatrix
#define dmatrix                  NR_dmatrix
#define imatrix                  NR_imatrix
#define submatrix                NR_subdmatrix
#define convert_matrix           NR_convert_dmatrix
#define subdmatrix               NR_subdmatrix
#define convert_dmatrix          NR_convert_dmatrix
#define d3tensor                 NR_d3tensor
#define free_vector              NR_free_dvector
#define free_ivector             NR_free_ivector
#define free_fvector             NR_free_fvector
#define free_dvector             NR_free_dvector
#define free_cvector             NR_free_cvector
#define free_lvector             NR_free_lvector
#define free_matrix              NR_free_dmatrix
#define free_imatrix             NR_free_imatrix
#define free_fmatrix             NR_free_fmatrix
#define free_dmatrix             NR_free_dmatrix
#define free_submatrix           NR_free_subdmatrix
#define free_convert_matrix      NR_free_convert_dmatrix
#define free_subdmatrix          NR_free_subdmatrix
#define free_convert_dmatrix     NR_free_convert_dmatrix
#define free_d3tensor            NR_free_d3tensor


void nrerror(char error_text[]); 
void nrerror2( char error_text[] );    /* Non-exit version JNI 16-Nov-00 */
float *fvector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
float **fmatrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
double **subdmatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
double **convert_dmatrix(double *a, long nrl, long nrh, long ncl, long nch);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_fvector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_fmatrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_subdmatrix(double **b, long nrl, long nrh, long ncl, long nch);
void free_convert_dmatrix(double **b, long nrl, long nrh, long ncl, long nch);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);




/////////////////////////////////////////////////////////////////////////////
// Comment top and comment out this to use personal version of dmatrix (e.g.)
// -- it links correctly with NRd.lib with NR_* calls
/////////////////////////////////////////////////////////////////////////////

/* 
void NR_nrerror(char error_text[]);
void NR_nrerror2( char error_text[] );    // Non-exit version JNI 16-Nov-00 
float *NR_vector(long nl, long nh);
int *NR_ivector(long nl, long nh);
unsigned char *NR_cvector(long nl, long nh);
unsigned long *NR_lvector(long nl, long nh);
double *NR_dvector(long nl, long nh);
float **NR_matrix(long nrl, long nrh, long ncl, long nch);
double **NR_dmatrix(long nrl, long nrh, long ncl, long nch);
int **NR_imatrix(long nrl, long nrh, long ncl, long nch);
double **NR_submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
double **NR_convert_matrix(double *a, long nrl, long nrh, long ncl, long nch);
double ***NR_f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void NR_free_vector(double *v, long nl, long nh);
void NR_free_ivector(int *v, long nl, long nh);
void NR_free_cvector(unsigned char *v, long nl, long nh);
void NR_free_lvector(unsigned long *v, long nl, long nh);
void NR_free_dvector(double *v, long nl, long nh);
void NR_free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void NR_free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void NR_free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void NR_free_submatrix(double **b, long nrl, long nrh, long ncl, long nch);
void NR_free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch);
void NR_free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);

*/





//static float sqrarg;
//#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

/*
static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))
*/

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#endif /* _NR_UTILS_H_ */




