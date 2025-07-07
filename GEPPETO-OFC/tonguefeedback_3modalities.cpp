#include <math.h>
#include "mylinpack.h"
#include "nr.h"
#include "nrutil.h"
#include "tongue_geom_audio_tactile_autoenc.h"
#include "autoencoder.h"
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"

#ifndef NO_MATLAB
#include "mex.h"  
#undef DB
#define DB mexPrintf("Got to line %d in mex file %s.\n",__LINE__,__FILE__);

#endif

void mexmatdisplay(double **M, int m, int n, const char *matname);
void mexvecdisplay(double *M, int m, const char *matname);

//extern int DOF, ns, nc;
extern int ns, nc;
int *ny_modal;
int Nmodal = 3;

bool feedback_initialized = false;

int ny;
extern const int NF = 3; // formants
extern const int NT;
extern const int Nproprio = 9;
extern const int D = 2;

static int offset;

const double DELTA = 1e-3; // for numerical differentiation 
double DELTA_weights[] = {0,0.6936,0.3982,1.0000,0.4401,0.6407};

extern double *averagetongue, **svdmatrix, **JD;
double **t_svdmatrix;
extern int Nsvd;
static double **deltas;


// not local vars to save time & deallocate cleanly
static double *tong, *tongvel, *tongacc;
static double *tongdeltas;
static double *form1;
static double *form2;
static double *prop1;
static double *prop2;
static double *tact1;
static double *tact2;

static double **Z1, **Z2, **Z3, **Z4, **Z5, **I1, **I2;
static double **temp1;
static double **temp2;
static double **temp3;
static double **temp4;
static double **temp5;
static double **temp6;
static double **temp7;
static double **temp8;
static double **formJacobian;
static double **tactJacobian;
static double **propJacobian;

static double *tongcontour;
static double area[N_VTSECTIONS];

void tongue2formantsAndGuess(double *stong, double *formantsGuess, double *formants)
{
  static double tubelength;
 
  decode(stong,tongcontour);
  vecadd(tongcontour,averagetongue,tongcontour,2*N_NODES);
  tongue2VTarea(tongcontour, area, tubelength);

  VTarea2formantLocationAndGuess(area, tubelength, formantsGuess, formants);
  hz2bark(formants,formants,NF);
}


void tongue2formantsAndGuess(double *stong, double *formantsGuess, double *formants, double *vt_area)
{
  static double tubelength;

  decode(stong,tongcontour);
  vecadd(tongcontour,averagetongue,tongcontour,2*N_NODES);
  tongue2VTarea(tongcontour, area, tubelength);
  veccopy(vt_area,area-1,N_VTSECTIONS);

  VTarea2formantLocationAndGuess(area, tubelength, formantsGuess, formants);
  hz2bark(formants,formants,NF);
}

void tongue2prop(double *stong, double *prop)
{
	double *svd_tong = dvector(1,Nsvd);
	double *svd_tongvel = dvector(1,Nsvd);
  double *svd_tongacc = dvector(1,Nsvd); 
	double *decoded_tong = dvector(1,2*N_NODES);
	double *decoded_tongvel = dvector(1,2*N_NODES);
	double *decoded_tongacc = dvector(1,2*N_NODES);

  veccopy(tong,stong,DOF);
  veccopy(tongvel,stong+DOF,DOF);
  vecsub(stong+DOF,stong+(D+1)*DOF,tongacc,DOF); 
 
	decode(tong,decoded_tong);
  matvecmult(2*N_NODES,JD,DOF,tongvel,decoded_tongvel);
  matvecmult(2*N_NODES,JD,DOF,tongacc,decoded_tongacc);

	matvecmult(Nsvd,t_svdmatrix,2*N_NODES,decoded_tong,svd_tong);
  matvecmult(Nsvd,t_svdmatrix,2*N_NODES,decoded_tongvel,svd_tongvel);
  matvecmult(Nsvd,t_svdmatrix,2*N_NODES,decoded_tongacc,svd_tongacc);

	veccopy(prop,svd_tong,Nproprio/3);
  veccopy(prop+Nproprio/3,svd_tongvel,Nproprio/3); 
	veccopy(prop+2*Nproprio/3,svd_tongacc,Nproprio/3);

	free_dvector(svd_tong,1,Nsvd);
	free_dvector(svd_tongvel,1,Nsvd);
	free_dvector(svd_tongacc,1,Nsvd);
	free_dvector(decoded_tong,1,2*N_NODES);
	free_dvector(decoded_tongvel,1,2*N_NODES);
	free_dvector(decoded_tongacc,1,2*N_NODES);
}





void tongue2formants(double *stong, double *formants)
{
  static double tubelength;
 
  decode(stong,tongcontour);
  vecadd(tongcontour,averagetongue,tongcontour,2*N_NODES);
  tongue2VTarea(tongcontour, area, tubelength);
  VTarea2formantLocation(area, tubelength, formants);
  hz2bark(formants,formants,NF);
}

void tongue2formants(double *stong, double *formants, double *vt_area)
{
  static double tubelength;

  decode(stong,tongcontour);
  vecadd(tongcontour,averagetongue,tongcontour,2*N_NODES);
  tongue2VTarea(tongcontour, area, tubelength);
  veccopy(vt_area,area-1,N_VTSECTIONS);
  VTarea2formantLocation(area, tubelength, formants);
	hz2bark(formants,formants,NF);
}

void tongue2formantsNumericalJacobian(double *stong, double **jacobian)
{
  int i, j;
  double delta;
  bool weird;
  int tried;

  veccopy(tong,stong,DOF);

  for (j = 1; j <= DOF; j++) {
    delta = DELTA*DELTA_weights[j];
    tried = 0;
    do {
      weird = false;
      tong[j] -= delta; 
      tongue2formants(tong, form1);
      tong[j] += 2*delta; 
      tongue2formants(tong, form2);
      tong[j] -= delta;  // back to original data
      for (i = 1; i <= NF; i++) {
        jacobian[i][j] = (form2[i]-form1[i])/2/delta;
        weird |= (fabs(form2[i]-form1[i]) > .1); 
      }
      delta /= 10;
      tried++;
    }
    while (weird && tried < 5);
    
    if (weird && tried == 5) {
      mexPrintf("Coord %d, delta %g:\n",j, delta);
      mexvecdisplay(form1,NF,"form1");
      mexvecdisplay(form2,NF,"form2");
      mexvecdisplay(stong,DOF,"tongpos");
      mexErrMsgTxt("Problem computing the derivatives of audio feedback. Stopping here.");
    }
    
  }
}

// computed once at init
void tongue2propJacobian(double **jacobian)
{
  double **proprioD;
  double **minusproprioD;
 
  proprioD = dmatrix(1,Nproprio/3,1,DOF);
  minusproprioD = dmatrix(1,Nproprio/3,1,DOF); 
  matmult(Nproprio/3,t_svdmatrix,2*N_NODES,JD,DOF,proprioD); 
  scalarmult(proprioD,Nproprio/3,DOF, -1., minusproprioD);

  zeros(jacobian,Nproprio,(D+2)*DOF);
  matassign(proprioD,Nproprio/3,DOF, jacobian,1,1);
  matassign(proprioD,Nproprio/3,DOF, jacobian,Nproprio/3+1,DOF+1);
  matassign(proprioD,Nproprio/3,DOF, jacobian,2*Nproprio/3+1,DOF+1);
  matassign(minusproprioD,Nproprio/3,DOF, jacobian,2*Nproprio/3+1,(D+1)*DOF+1);
  
  free_dmatrix(proprioD,1,Nproprio/3,1,DOF);
  free_dmatrix(minusproprioD,1,Nproprio/3,1,DOF);
}

void tongue2propNumericalJacobian(double *stong, double **jacobian)
{
  int i, j, k;
  
   double *propTong = dvector(1,(D+2)*DOF);

   veccopy(propTong,stong,(D+2)*DOF);

  for (j = 1; j <= (D+2)*DOF; j++) {
    k = (j-1)%DOF+1; 
    propTong[j] -= DELTA*DELTA_weights[k]; 
    tongue2prop(propTong, prop1);
    
    propTong[j] += 2*DELTA*DELTA_weights[k];
    tongue2prop(propTong, prop2);
    
    propTong[j] -= DELTA*DELTA_weights[k];
    for (i = 1; i <= Nproprio; i++) jacobian[i][j] = (prop2[i]-prop1[i])/(2*DELTA*DELTA_weights[k]); 
  }

	free_dvector(propTong,1,(D+2)*DOF);
}



void tongue2formantsNumericalJacobianAndGuess(double *stong, double *Fguess, double **jacobian)
{
  int i, j;
  veccopy(tong,stong,DOF);

  for (j = 1; j <= DOF; j++) {
    tong[j] -= DELTA*DELTA_weights[j]; 
    tongue2formantsAndGuess(tong, Fguess, form1);
    
    tong[j] += 2*DELTA*DELTA_weights[j]; 
    tongue2formantsAndGuess(tong, Fguess, form2);
    
    tong[j] -= DELTA*DELTA_weights[j]; 
    for (i = 1; i <= NF; i++) jacobian[i][j] = (form2[i]-form1[i])/2/(DELTA*DELTA_weights[j]); 
  }
}


void tongue2tactileNumericalJacobian(double *stong, double **jacobian)
{
  int i, j;
  veccopy(tong,stong,DOF);

  for (j = 1; j <= DOF; j++) {
    tong[j] -= DELTA*DELTA_weights[j];
    tongue2tactile(tong, tact1);
    tong[j] += 2*DELTA*DELTA_weights[j];
    tongue2tactile(tong, tact2);
    tong[j] -= DELTA*DELTA_weights[j]; // back to original data
    for (i = 1; i <= NT; i++) jacobian[i][j] = (tact2[i]-tact1[i])/(2*DELTA*DELTA_weights[j]);
  }
}

void feedbackAndGuess(double* s, double *u, double *formantsGuess, double *y)
{
  int i,j,k;
  tongue2formantsAndGuess(s+offset, formantsGuess, y); // order : auditory, proprio, tactile
  tongue2prop(s+offset, y+NF);
  tongue2tactile(s+offset, y+NF+Nproprio);
}

void feedbackAndGuess(double* s, double *u, double *y, double *formantsGuess, double *vt_area)
{
  int i,j,k;
  tongue2formantsAndGuess(s+offset, formantsGuess, y, vt_area); // order : auditory, proprio, tactile
  tongue2prop(s+offset, y+NF);
  tongue2tactile(s+offset, y+NF+Nproprio);
}


void feedback(double* s, double *u, double *y)
{
  int i,j,k;
  tongue2formants(s+offset, y) ; // order : auditory, proprio, tactile
  tongue2prop(s+offset, y+NF);
  tongue2tactile(s+offset, y+NF+Nproprio);
}

void feedback(double* s, double *u, double *y, double *vt_area)
{
  int i,j,k;
  tongue2formants(s+offset, y, vt_area) ; // order : auditory, proprio, tactile
  tongue2prop(s+offset, y+NF);
  tongue2tactile(s+offset, y+NF+Nproprio);
}

#include <unistd.h> // for usleep

void feedback_ds(double* s, double *u, double **Hs)
{
  zeros(temp2,Nproprio,ns-offset);
  zeros(temp3,NT,ns-offset);
 
  tongue2formantsNumericalJacobian(s+offset, formJacobian);
  tongue2tactileNumericalJacobian(s+offset, tactJacobian);
  tongue2propNumericalJacobian(s+offset, propJacobian);

  matconcat(formJacobian,NF,DOF,Z1,NF,ns-offset-DOF,temp1);
  matconcat(propJacobian,Nproprio,(D+2)*DOF,Z2,Nproprio,ns-offset-(D+2)*DOF,temp2);
  matconcat(tactJacobian,NT,DOF,Z3,NT,ns-offset-DOF,temp3);
  
  matconcat(temp1,NF,ns-offset,temp2,Nproprio,ns-offset,temp7,'v');
  matconcat(temp7,NF+Nproprio,ns-offset,temp3,NT,ns-offset,temp8,'v');
  matconcat(Z4,ny,offset,temp8,ny,ns-offset,Hs);
  
  bool weird = false;
  for (int i = 1; i <= ny; i++) for (int j = 1; j <= ns; j++) if (isnan(Hs[i][j])) weird = true;
  if (weird) mexmatdisplay(Hs,ny,ns,"Hoho weird Hs");
}

void feedback_dsAndGuess(double* s, double *Fguess, double *u, double **Hs)
{
  zeros(temp2,Nproprio,ns-offset);
  zeros(temp3,NT,ns-offset);

  tongue2formantsNumericalJacobianAndGuess(s+offset, Fguess, formJacobian);
  tongue2tactileNumericalJacobian(s+offset, tactJacobian);

  matconcat(formJacobian,NF,DOF,Z1,NF,ns-offset-DOF,temp1);
  matconcat(propJacobian,Nproprio,(D+2)*DOF,Z2,Nproprio,ns-offset-(D+2)*DOF,temp2);
  matconcat(tactJacobian,NT,DOF,Z3,NT,ns-offset-DOF,temp3);
  
  matconcat(temp1,NF,ns-offset,temp2,Nproprio,ns-offset,temp7,'v');
  matconcat(temp7,NF+Nproprio,ns-offset,tactJacobian,NT,ns-offset,temp8,'v');
  matconcat(Z4,ny,offset,temp8,ny,ns-offset,Hs); 
}

void init_feedback() 
{
  ny = NF + NT + Nproprio ;   // full proprio feedback + audio feedback + tactile
  ny_modal = ivector(1,Nmodal);
  ny_modal[1] = NF;
  ny_modal[2] = Nproprio;
  ny_modal[3] = NT;
  offset = 1;

	t_svdmatrix = dmatrix(1,Nsvd,1,2*N_NODES);
	transpose(svdmatrix,2*N_NODES,Nsvd,t_svdmatrix); // for proprio

  deltas = eye(DOF);
  scalarmult(deltas, DOF, DOF, DELTA);
 
  tong  = dvector(1,DOF);
  tongvel  = dvector(1,DOF);
	tongacc  = dvector(1,DOF);
  tongdeltas  = dvector(1,DOF);
  form1 = dvector(1,NF);
  form2 = dvector(1,NF);
  prop1 = dvector(1,Nproprio);
  prop2 = dvector(1,Nproprio);
  tact1 = dvector(1,NT);
  tact2 = dvector(1,NT);

  Z1 = dmatrix(1,NF,1,ns-offset-DOF);
  Z2 = dmatrix(1,Nproprio,1,ns-offset-(D+2)*DOF);
  Z3 = dmatrix(1,NT,1,ns-offset-DOF);
  Z4 = dmatrix(1,ny,1,offset);
  temp1 = dmatrix(1,NF,1,ns-offset);
  temp2 = dmatrix(1,Nproprio,1,ns-offset);
  temp3 = dmatrix(1,NT,1,ns-offset);
  temp7 = dmatrix(1,NF+Nproprio,1,ns-offset);
  temp8 = dmatrix(1,ny,1,ns-offset);
  formJacobian = dmatrix(1,NF,1,DOF);
  tactJacobian = dmatrix(1,NT,1,DOF);
  propJacobian = dmatrix(1,Nproprio,1,(D+2)*DOF);
  tongue2propJacobian(propJacobian);

  tongcontour = dvector(1,2*N_NODES);

  zeros(Z1,NF,ns-offset-DOF);
  zeros(Z2,Nproprio,ns-offset-(D+2)*DOF);
  zeros(Z3,NT,ns-offset-DOF);
  zeros(Z4,ny,offset);

  feedback_initialized = true;
}


void free_feedback_variables()
{
  free_ivector(ny_modal,1,3);

  free_dvector(tong,1,DOF);
  free_dvector(tongvel,1,DOF);
	free_dvector(tongacc,1,DOF);
  free_dvector(tongdeltas,1,DOF);
  free_dvector(form1,1,NF);
  free_dvector(form2,1,NF);
  free_dvector(prop1,1,Nproprio);
  free_dvector(prop2,1,Nproprio);
  free_dvector(tact1,1,NT);
  free_dvector(tact2,1,NT);

  free_dmatrix(temp1,1,NF,1,ns-offset);
  free_dmatrix(temp2,1,Nproprio,1,ns-offset);
  free_dmatrix(temp3,1,NT,1,ns-offset);
  free_dmatrix(temp7,1,NF+Nproprio,1,ns-offset);
  free_dmatrix(temp8,1,ny,1,ns-offset);
  free_dmatrix(formJacobian,1,NF,1,DOF);
  free_dmatrix(tactJacobian,1,NT,1,DOF);
  free_dmatrix(propJacobian,1,Nproprio,1,(D+2)*DOF);
  free_dmatrix(Z1,1,NF,1,ns-offset-DOF);
  free_dmatrix(Z2,1,Nproprio,1,ns-offset-(D+2)*DOF);
  free_dmatrix(Z3,1,NT,1,ns-offset-DOF);
  free_dmatrix(Z4,1,ny,1,offset);

  free_dmatrix(t_svdmatrix,1,Nsvd,1,2*N_NODES);

  free_dvector(tongcontour,1,2*N_NODES);

}
