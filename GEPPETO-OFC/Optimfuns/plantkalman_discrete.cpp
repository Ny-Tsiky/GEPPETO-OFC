
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mylinpack.h"
#include "nr.h"
#include "nrutil.h"   
#include "kalman8.h"
#include "plantfuns_via.h"
#include "plantfeedback.h"

#ifndef NO_MATLAB
#include "mex.h"  
#endif

#undef DB
#define DB mexPrintf("Got to line %d in mex file %s.\n",__LINE__,__FILE__);

bool debuggo = false;


extern int ns, nc, DOF, ny, Nmodal, NF;
extern int *ny_modal;

static double deltat;

inline void veccopy(int *aout, int *ain, int n)
{
  int i;

  for (i=1;i<=n;i++)
      aout[i]=ain[i];
}

void Bark2Hz(double *bark_y, double *y, int NF);

void kalman_A(double* s, double t, double** A) 
{
  statetransition_matrix(t,A);
}

void kalman_B(double* s, double* u, double t, double** B)
{
  control_matrix(t,B);
}

void kalman_H(double* s, double t, double** H) 
{
  observation_matrix(t,H);
}


void kalman_f(double* s, double* u_old, double* u, double t, double* snext, double *out_HS, double *out_CS)
{
  discrete_forward(u_old,u,s,t,deltat,snext, 0, out_HS, out_CS, 0);
}

void kalman_df_ds(double* s, double* u, double t, double** Fs)
{
  discrete_adjoint_dx2(u,s,t,deltat,Fs, 0);
}

void kalman_df_du(double* s, double* u_old, double* u, double t, double** Fu)
{
  discrete_adjoint_du2(u_old,u,s,t,deltat,Fu, 0);  
}

void kalman_h(double* s, double *u, double t, double *y)
{
  feedback(s,u,y); 
}

void kalman_dh_ds(double* s, double *u, double t, double** Hs)
{ 
  feedback_ds(s,u,Hs); 
}


void free_kalman_filter(kalman_type kalmtype)
{
  if (kalmtype == DISCRETE_EKF) free_discrete_ekf_variables();
  else if (kalmtype == DISCRETE) free_kalman_variables();
}

void init_kalman_filter(double *s, noise_params *noise, double dt, double *delay, kalman_type kalmtype)  
{
  double *ustart;
  double *ystart;
  init_feedback();
  ustart = vecextract(s,ns-nc+1,ns); 
  ystart = dvector(1,ny);
  
  feedback(s,ustart,ystart);
  
  switch (kalmtype) {
    case DISCRETE_EKF:
    
    static ekf_params ekf;
    ekf.noise_type = MIXED;
    ekf.sigQp = noise->var_u_sdn;   
    ekf.sigRp = noise->var_s_sdn;
    ekf.sigQ = noise->var_u_add;
    ekf.sigR = noise->var_s_add;
    ekf.sigS = noise->var_x_add;
    ekf.Nmodal = Nmodal; 
    ekf.nx = ns;
    ekf.ny = ny;
    ekf.ny_modal = ny_modal;
    ekf.nu = nc;
    ekf.delay = delay;
    ekf.dt = deltat = dt;
    ekf.f = kalman_f;
    ekf.dfdx = kalman_df_ds;
    ekf.dfdu = kalman_df_du;
    ekf.h = kalman_h;
    ekf.dhdx = kalman_dh_ds;
    init_ekf(&ekf);
    init_discrete_ekf_variables_static(s, ustart, ystart);

    break;

    case DISCRETE:
    static kalman_params kalm;
    
    kalm.noise_type = MIXED;
    kalm.sigRp = noise->var_s_sdn;
    kalm.sigQ = noise->var_u_add;
    kalm.sigR = noise->var_s_add;
    kalm.sigS = noise->var_x_add;
    kalm.Nmodal = Nmodal;

    kalm.nx = ns;
    kalm.ny = ny;
    kalm.nu = nc;
    kalm.delay = dvector(1,3);
    veccopy(kalm.delay, delay, 3); 
    kalm.dt = deltat = dt;
    kalm.Afun = kalman_A;
    kalm.Bfun = kalman_B;
    kalm.Hfun = kalman_H;
    init_kalman(&kalm);
    init_kalman_variables_static(s, ustart, ystart);
    break;
          
    case CONTINUOUS:
    init_kalman(&kalm);
    break;
          
    case CONTINUOUS_EKF:
    init_kalman(&kalm);
    break;

  }

  free_dvector(ustart,1,nc);
  free_dvector(ystart,1,ny);
}



void update_kalman_noises(noise_params *noise, kalman_type kalmtype)  
{
    switch (kalmtype) {
    case DISCRETE_EKF:
    
    static ekf_params ekf;
    ekf.noise_type = MIXED;
    ekf.sigQp = noise->var_u_sdn;
    ekf.sigRp = noise->var_s_sdn;
    ekf.sigQ = noise->var_u_add;
    ekf.sigR = noise->var_s_add;
    ekf.sigS = noise->var_x_add;

    update_ekf_noiseparams(&ekf);

    break;

    case DISCRETE:
    
    case CONTINUOUS:
            break;
          
    case CONTINUOUS_EKF:
            break;
  }
}

