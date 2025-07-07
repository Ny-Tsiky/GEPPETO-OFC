#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mylinpack.h"
#include "nr.h"
#include "nrutil.h"   // to shunt the better but longer-named NR_* prefix functions
#include <limits.h>
#include "kalman8.h"

#include <unistd.h> // for usleep
#include <iostream> // for exception handling

#define USE_MATLAB

#include "mex.h" 

#undef DB
#define DB mexPrintf("Got to line %d in mex file %s.\n",__LINE__,__FILE__);

void mexmatdisplay(double **M, int m, int n, const char *matname);
void mexvecdisplay(double *M, int m, const char *matname);
bool alert = false;

void inverse(double **X, double **X_inv, int m)
{
  double x[m*m];
  double res[m*m];
  double smax;
  int b_i;
  int i;
  int j;
  int jA;
  int jp1j;
  int k;
  int kAcol;
  int x_tmp;
  int tmp = 0;
  signed char ipiv[m];
  signed char p[m];
  for (i = 1; i <= m; i++) {
      for (j = 1; j <= m; j++) {
        res[tmp] = 0.0;
        x[tmp] = X[i][j];
          tmp++;
      }
  }
  for (i = 0; i < m; i++) {
    ipiv[i] = (signed char)(i + 1);
  }
  for (j = 0; j < m-1; j++) {
    int b_tmp;
    int mmj_tmp;
    mmj_tmp = (m-2) - j;
    b_tmp = j * (m+1);
    jp1j = b_tmp + 2;
    jA = m - j;
    kAcol = 0;
    smax = fabs(x[b_tmp]);
    for (k = 2; k <= jA; k++) {
      double s;
      s = fabs(x[(b_tmp + k) - 1]);
      if (s > smax) {
        kAcol = k - 1;
        smax = s;
      }
    }
    if (x[b_tmp + kAcol] != 0.0) {
      if (kAcol != 0) {
        jA = j + kAcol;
        ipiv[j] = (signed char)(jA + 1);
        for (k = 0; k < m; k++) {
          kAcol = j + k * m;
          smax = x[kAcol];
          x_tmp = jA + k * m;
          x[kAcol] = x[x_tmp];
          x[x_tmp] = smax;
        }
      }
      i = (b_tmp - j) + m;
      for (b_i = jp1j; b_i <= i; b_i++) {
        x[b_i - 1] /= x[b_tmp];
      }
    }
    jA = b_tmp;
    for (kAcol = 0; kAcol <= mmj_tmp; kAcol++) {
      smax = x[(b_tmp + kAcol * m) + m];
      if (smax != 0.0) {
        i = jA + (m+2);
        jp1j = (jA - j) + m*2;
        for (x_tmp = i; x_tmp <= jp1j; x_tmp++) {
          x[x_tmp - 1] += x[((b_tmp + x_tmp) - jA) - (m+1)] * -smax;
        }
      }
      jA += m;
    }
  }
  for (i = 0; i < m; i++) {
    p[i] = (signed char)(i + 1);
  }
  for (k = 0; k < m-1; k++) {
    signed char i1;
    i1 = ipiv[k];
    if (i1 > k + 1) {
      jA = p[i1 - 1];
      p[i1 - 1] = p[k];
      p[k] = (signed char)jA;
    }
  }
  for (k = 0; k < m; k++) {
    x_tmp = m * (p[k] - 1);
    res[k + x_tmp] = 1.0;
    for (j = k + 1; j < m+1; j++) {
      i = (j + x_tmp) - 1;
      if (res[i] != 0.0) {
        jp1j = j + 1;
        for (b_i = jp1j; b_i < m+1; b_i++) {
          jA = (b_i + x_tmp) - 1;
          res[jA] -= res[i] * x[(b_i + m * (j - 1)) - 1];
        }
      }
    }
  }
  for (j = 0; j < m; j++) {
    jA = m * j;
    for (k = m-1; k >= 0; k--) {
      kAcol = m * k;
      i = k + jA;
      smax = res[i];
      if (smax != 0.0) {
        res[i] = smax / x[k + kAcol];
        for (b_i = 0; b_i < k; b_i++) {
          x_tmp = b_i + jA;
          res[x_tmp] -= res[i] * x[b_i + kAcol];
        }
      }
    }
  };
    tmp = 0;
    for (i = 1; i <= m; i++) {
      for (j = 1; j <= m; j++) {
         X_inv[i][j] = res[tmp];
          tmp++;
      }
  }
}
// ===============================================================================


// Kalman

static double **A, **B, **H, **Q, **R, **S, **Qp, **Rp;
static double **tA, **tB, **tH;
static double **E, **K, **tK;
static double *xhat = NULL, **yqueue, *y_old;
static int m, n, nx, nd, p, delaysteps;
static int first_time = 1;
static double **security;
const double security_shift = 1e-10;

static int Nmodal;
static int *modal_delay_steps, *m_modal_start, *m_modal;

static int delayed = 0; // flag for delayed kalman

// EKF
static void (*f)(double *x, double *u_old, double *u, double t,  double *xdot, double *out_HS, double *out_CS);
static void (*dfdx)(double *x, double *u, double t,  double **Fx);
static void (*dfdu)(double *x, double *u_old, double *u, double t, double **Fu);
static void (*h)(double *x, double *u, double t, double *z);
static void (*dhdx)(double *x, double *u, double t, double **H);
static double **tu, **txhat, **ty;
static double **ku, **kxhat, **ky; // without time column (useless in discrete case)
static double *mean_u, *mean_y, *mean_u_saved, *mean_y_saved;
static double internal_time;
static int Nq, tindex, refresh_queues = 0; 
static double *modal_delay, static_delay, dt, internal_delay;
static double sigQ, *sigR, sigQp, *sigRp, sigS;
//static double sigR2, sigRp2;
static enum noise noise_type;

const double SIG_NO_FEEDBACK = 1e8; // variance of "no feedback" pseudosignal

extern int LSTM_SIZE;

// odeint parameters
const double eps = 1e-3;  
const double hmin = 1e-8;
static int kmax;
static int keuler = 100;

void update_gaussian(int *available_feedback);
void update_sdn(double *u, double *y, int *available_feedback, double *mu_saved, double *my_saved, bool set_mean);
void update_poisson(double *u, double *y, int *available_feedback, double *mu_saved, double *my_saved, bool set_mean);


void calc_gains()
{
  double **temp1, **temp2, **temp3, **temp4, **temp5, **temp6, **temp7, **temp8, **temp9;

  temp1 = dmatrix(1,n,1,m);
  temp2 = dmatrix(1,m,1,m);
  temp8 = dmatrix(1,m,1,m);
  temp3 = dmatrix(1,n,1,n);
  temp4 = dmatrix(1,n,1,n);
  temp5 = dmatrix(1,n,1,n);
  temp6 = dmatrix(1,n,1,p);
  temp7 = dmatrix(1,n,1,n);
  temp9 = dmatrix(1,p,1,p);


  matmult(n,A,n,E,n,temp4);
  matmult(n,temp4,n,tA,n,temp5);  
  matadd(Q,Qp,temp9,p);
  matmult(n,B,p,temp9,p,temp6);
  matmult(n,temp6,p,tB,n,temp7);
  matadd(temp5,temp7,temp7,n);
  matadd(temp7,S,E,n);

  matmult(n,E,n,tH,m,temp1);
  matmult(m,H,n,temp1,m,temp2);  
  matadd(temp2,Rp,temp2,m);
  matadd(temp2,R,temp2,m);
  matadd(temp2,security,temp2,m);  // add some diag component to ensure invertibility
  inv(temp2,m,temp2);
  matmult(n,temp1,m,temp2,m,K);

  eye(temp3,n);
  matmult(n,K,m,H,n,temp4);
  matsub(temp3,temp4,temp4,n);  
  matmult(n,temp4,n,E,n,temp3);
  copy(E,temp3,n);

  free_dmatrix(temp1,1,n,1,m);
  free_dmatrix(temp2,1,m,1,m);
  free_dmatrix(temp8,1,m,1,m);
  free_dmatrix(temp3,1,n,1,n);
  free_dmatrix(temp4,1,n,1,n);
  free_dmatrix(temp5,1,n,1,n);
  free_dmatrix(temp6,1,n,1,p);
  free_dmatrix(temp7,1,n,1,n);
  free_dmatrix(temp9,1,p,1,p);
}

void update_kalman(double *u, double *xest, double t, kalman_params *params)
{
  if (!delayed) {
    params->Afun(xest,t,A);
    params->Bfun(xest,u,t,B);
    params->Hfun(xest,t,H);
  }
  else {
    double **A0, **B0, **H0;
    double **Z1,**Z2, **Z3, **Z4, **Z5, **I, **temp1, **temp2, **temp3;

    A0 = dmatrix(1,nx,1,nx);
    B0 = dmatrix(1,nx,1,p);
    H0 = dmatrix(1,m,1,nx);

    params->Afun(xest,t,A0);
    params->Bfun(xest,u,t,B0);
    params->Hfun(xest,t,H0);

    Z1 = zeros(nx,nd);
    Z2 = zeros(nd,nx);
    I = eye(nd);
    temp1 = matconcat(A0,nx,nx,Z1,nx,nd);
    temp2 = matconcat(I, nd,nd, Z2, nd,nx);
    matconcat(temp1, nx,n, temp2, nd,n, A,'v');

    Z3 = zeros(nd,p);
    matconcat(B0,nx,p,Z3,nd,p, B,'v');  

    Z4 = zeros(m,nd);    
    matconcat(Z4,m,nd, H0,m,nx, H);

    free_dmatrix(temp1,1,nx,1,n);
    free_dmatrix(temp2,1,nd,1,n);
    free_dmatrix(temp3,1,nx,1,nx);
    free_dmatrix(Z1,1,nx,1,nd);
    free_dmatrix(Z2,1,nd,1,nx);
    free_dmatrix(Z3,1,nd,1,p);
    free_dmatrix(Z4,1,m,1,nd);
    free_dmatrix(I,1,nd,1,nd);
    free_dmatrix(A0,1,nx,1,nx);
    free_dmatrix(B0,1,nx,1,p);
    free_dmatrix(H0,1,m,1,nx);
  }

}


void init_kalman_variables_static(double *xest, double *ustart, double *ystart)
{
  int i;

  if (delayed) {
    y_old = dvector(1,m);
    
    yqueue = dmatrix(1,delaysteps,1,m);
    for (i=1; i<=delaysteps; i++)
      veccopy(yqueue[i],ystart,m);
  }

  xhat = dvector(1,n);
  veccopy(xhat,xest,nx);
  
  zeros(E,n,n); 
  
  for (i=1; i<delaysteps; i++)
    vecconcat(xhat,nx*i, xest, nx, xhat);

}


void init_kalman(kalman_params *params)
{
  int i, modal;
  double **Z1, **Z2, **temp1, **temp2;
  
  if (params == NULL) 
    fprintf(stderr,"init_kalman: Error, empty parameter structure\n");

  Nmodal = params->Nmodal;
  modal_delay = dvector(1,Nmodal);
  modal_delay_steps = ivector(1,Nmodal);
  m_modal = ivector(1,Nmodal);
  m_modal_start = ivector(1,Nmodal+1);
  
  dt = params->dt;
  
  static_delay = 0;
  for (modal = 1; modal <= Nmodal; modal++) m_modal[modal] = params->ny_modal[modal]; 
  m_modal_start[1] = 0;
  for (modal = 1; modal <= Nmodal; modal++) {
    if (modal > 1) m_modal_start[modal] = m_modal_start[modal-1] + m_modal[modal-1];
    modal_delay_steps[modal] = params->delay[modal];
    if (modal_delay_steps[modal] > delaysteps) delaysteps = modal_delay_steps[modal]; 
  }
  
  m = m_modal_start[Nmodal+1]  = m_modal_start[Nmodal] + m_modal[Nmodal]; 
  
  n = params->nx;
  p = params->nu;

  nd = n*delaysteps;
  nx = n;
  n += nd;

  delayed = delaysteps>0;

  //  if (delayed) {
  A = dmatrix(1,n,1,n);
  B = dmatrix(1,n,1,p);
  H = dmatrix(1,m,1,n); 
  tA = dmatrix(1,n,1,n);
  tB = dmatrix(1,p,1,n);
  tH = dmatrix(1,n,1,m);
      
    //  }

  E = zeros(n,n);
  K = dmatrix(1,n,1,m);
  
  transpose(A,n,n,tA);
  transpose(B,n,p,tB);  
  transpose(H,m,n,tH);

  noise_type = params->noise_type;


  if (!delayed) {
    S = eye(n);
    scalarmult(S,n,n,params->sigS);
  } 
  else {
    Z1 = zeros(nd,n);
    temp1 = eye(nx);
    scalarmult(temp1,nx,nx,params->sigS);
    Z2 = zeros(nx,nd);
    temp2 = dmatrix(1,nx,1,n);
    matconcat(temp1,nx,nx,Z2,nx,nd,temp2);
    S = matconcat(temp2, nx,n, Z1, nd,n,'v');
  }

  double **Rmodal;
  for (int modal = 1; modal <= Nmodal; modal++) {
    Rmodal = eye(m_modal[modal]);
  	scalarmult(Rmodal, m_modal[modal], m_modal[modal], sigR[modal]);
    matassign(Rmodal, m_modal[modal], m_modal[modal], R, m_modal_start[modal]+1, m_modal_start[modal]+1);
    free_dmatrix(Rmodal, 1,m_modal[modal], 1,m_modal[modal]);
  }
  Q = eye(p);
  scalarmult(Q,p,p,params->sigQ);
  Qp = zeros(p,p);   // goret
  Rp = zeros(m,m);   // goret itou


  security = eye(m);
  scalarmult(security,m,m,security_shift);

  free_dmatrix(Z1,1,nd,1,n);
  free_dmatrix(Z2,1,nx,1,nd);
  free_dmatrix(temp1,1,nx,1,nx);
  free_dmatrix(temp2,1,nx,1,n);

  first_time = 0;
}



void reset_kalman()
{
  first_time = 1;
}

void kalman(double t, double *u, double *y, double *xest, int *consider_feedback, kalman_params *params, int reset)
{
  double *temp0, *temp1, *temp2, *temp3, *temp4, *temp5;
  static int index;
  
  if (first_time || reset) {
    if (!first_time) free_kalman_variables();
    init_kalman(params);
    init_kalman_variables_static(xest,u,y);
    first_time = 0;
    index = 0;
  }

  temp1 = dvector(1,n);
  temp2 = dvector(1,n);
  temp3 = dvector(1,m);
  temp4 = dvector(1,m);
  temp5 = dvector(1,n);          

  if (delayed) {
    //update sensory queue
    index = index % delaysteps + 1;
    veccopy(y_old,yqueue[index],m); 
    veccopy(yqueue[index],y,m); 
  }
  else 
    y_old = y;
  

  switch (noise_type) {      
  case POISSON: update_poisson(u,y,consider_feedback,mean_u,mean_y,false); break;
  case GAUSSIAN: update_gaussian(consider_feedback); break;
  case MIXED:
  case SDN: update_sdn(u,y,consider_feedback,mean_u,mean_y,false); break;
  }    

  update_kalman(u,xest,t, params);

  calc_gains();

  matvecmult(n,A,n,xhat,temp1);
  matvecmult(n,B,p,u,temp2);
  vecadd(temp1,temp2,xhat,n); 

  matvecmult(m,H,n,xhat,temp3); vecdisplay(temp0,m,"Hxhat");
  vecsub(y_old,temp3,temp4,m); vecdisplay(y_old,m,"y_old");
  matvecmult(n,K,m,temp4,temp5); // error term
  if (consider_feedback)
    vecadd(xhat,temp5,xhat,n); 


  veccopy(xest,xhat,nx);  // first nx coords = last estimate
  
  free_dvector(temp1,1,n);
  free_dvector(temp2,1,n);
  free_dvector(temp3,1,m);
  free_dvector(temp4,1,m);
  free_dvector(temp5,1,n);    

}


void free_kalman_variables()
{
  if (modal_delay) { // minimal safety
     free_dvector(modal_delay,1,Nmodal);
     free_ivector(modal_delay_steps,1,Nmodal);
     free_ivector(m_modal,1,Nmodal);
     free_ivector(m_modal_start,1,Nmodal+1);
  }
  
  if (K) {
    free_dmatrix(E,1,n,1,n);
    free_dmatrix(K,1,n,1,m);

    free_dmatrix(A,1,n,1,n);
    free_dmatrix(B,1,n,1,p);
    free_dmatrix(H,1,m,1,n);
    free_dmatrix(tA,1,n,1,n);
    free_dmatrix(tB,1,p,1,n);
    free_dmatrix(tH,1,n,1,m);

    free_dvector(xhat,1,n);

    if (delayed) {
      free_dvector(y_old,1,m);
      free_dmatrix(yqueue,1,delaysteps,1,m);
    }
  }
}


//*********************************************************************************
//*********************************************************************************
//*********************************************************************************


// XXXX TODO/ adapt to modality-specific FB availability
void calc_ekf_gains(double t, double *xhat, double *u)
{
  double **temp1, **temp2, **temp3, **temp4;
  double **tH;


  tH = dmatrix(1,n,1,m);

  temp1 = dmatrix(1,n,1,n);
  temp2 = dmatrix(1,m,1,m);
  temp3 = dmatrix(1,n,1,n);
  temp4 = dmatrix(1,n,1,n);
  
  dhdx(xhat, u, t, H);

  transpose(H,m,n,tH);
  matmult(n,E,n,tH,m,temp1);
  matadd(R,Rp,temp2,m);
  matadd(temp2,security,temp2,m);  // add some diag component to ensure invertibility
  inv(temp2,m,temp2);
  matmult(n,temp1,m,temp2,m,K);

  free_dmatrix(tH,1,n,1,m);

  free_dmatrix(temp1,1,n,1,n);
  free_dmatrix(temp2,1,m,1,m);
  free_dmatrix(temp3,1,n,1,n);
  free_dmatrix(temp4,1,n,1,n);
}



void update_ekf_noiseparams(ekf_params *params)
{
  int modal, i, k;
  double *diR = dvector(1,m);
    

  sigQ = params->sigQ;
  sigQp = params->sigQp;
  for (i = 1; i <= Nmodal; i++) {
    sigR[i] = params->sigR[i];
    sigRp[i] = params->sigRp[i];
  }
  sigS = params->sigS;

  eye(Q,p);
  eye(S,n);
  scalarmult(Q,p,p,sigQ);
  scalarmult(S,n,n,sigS);

  for (k = 1, modal = 1; modal <= Nmodal; modal++) 
    for (i = 0; i < m_modal[modal]; i++, k++) {
      diR[k] = params->sigR[modal];
    }
  diag(R,diR,m);  

  zeros(Qp,p,p);
  zeros(Rp,m,m);  
  
  zeros(E,n,n); 
  
  free_dvector(diR,1,m);
}


void update_ekf(ekf_params *params)
{
    
  f = params->f; 
  dfdx = params->dfdx;
  dfdu = params->dfdu;
  h = params->h;
  dhdx = params->dhdx;

  update_ekf_noiseparams(params);
}



void get_ekf_params(ekf_params *params)
{
	params->Nmodal = Nmodal;
  
  params->ny_modal = ivector(1,Nmodal);
  params->delay = dvector(1,Nmodal);
	for (int i=1; i<= Nmodal; i++) {
    params->ny_modal[i] = m_modal[i];
    params->delay[i] = modal_delay[i];
  }
  
  params->nx = n;
  params->nu = p;

  params->dt = dt;
  params->noise_type = noise_type;

  params->f = f;
  params->dfdx = dfdx;
  params->dfdu = dfdu;
  params->h = h;
  params->dhdx = dhdx;


  params->sigQ = sigQ;
  params->sigQp = sigQp;
  for (int modal = 1; modal <= Nmodal; modal++) {
    params->sigR[modal] = sigR[modal];
    params->sigRp[modal] = sigRp[modal];
  }
  params->sigS = sigS;
}

void init_ekf_variables(double **txest, double **tustart, double **tystart)
{
  int i;
  double **u, **t, **x, **y;

    xhat = dvector(1,n);
    tu = dmatrix(1,Nq,1,p+1);
    txhat = dmatrix(1,Nq,1,n+1);
    ty = dmatrix(1,Nq,1,m+1);

    vecextract(tu[Nq],2,n+1, xhat); // ??? bizarre

  copy(txhat, txest, Nq, n+1);
  copy(tu, tustart, Nq, p+1);
  copy(ty, tystart, Nq, m+1);

  zeros(E,n,n); 
  
  internal_time = tustart[Nq][1];

  tindex = 1;

  free_dmatrix(t,1,Nq,1,1);

  first_time = 0;
}



void init_ekf_variables_static(double *xest, double *ustart, double *ystart, double **Estart)
{
  int i;
  double **u, **t, **x, **y;

    xhat = dvector(1,n);
    tu = dmatrix(1,Nq,1,p+1);
    txhat = dmatrix(1,Nq,1,n+1);
    ty = dmatrix(1,Nq,1,m+1);

  veccopy(xhat,xest,n);

  t = matlinspace(-Nq*dt,-dt,Nq,'c');

  u = dmatrix(1,Nq,1,p);
  for ( i=1; i<=Nq; i++) veccopy(u[i],ustart,p);
  matconcat(t,Nq,1,u,Nq,p, tu); // control value circular queue

  x = dmatrix(1,Nq,1,n);
  for ( i=1; i<=Nq; i++) veccopy(x[i],xest,n);
  matconcat(t,Nq,1,x,Nq,n, txhat); // state estimate circular queue
  
  y = dmatrix(1,Nq,1,m);
  for ( i=1; i<=Nq; i++) veccopy(y[i],ystart,m);
  matconcat(t,Nq,1,y,Nq,m, ty); // sensory feedback circular queue
  
  // reset covariance estimate
  if (Estart) copy(E, Estart, n);
  else zeros(E,n,n); 

  internal_time = 0;

  tindex = 1;

  free_dmatrix(t,1,Nq,1,1);
  free_dmatrix(u,1,Nq,1,p);
  free_dmatrix(x,1,Nq,1,n);
  free_dmatrix(y,1,Nq,1,m);

  first_time = 0;
}


// HYP: t = 0 at start
void init_ekf(ekf_params *params) 
{
  int i, modal;
  if (params == NULL) 
    fprintf(stderr,"init_ekf: Error, empty parameter structure\n");
  Nmodal = params->Nmodal;
  modal_delay = dvector(1,Nmodal);
  modal_delay_steps = ivector(1,Nmodal);
  m_modal = ivector(1,Nmodal);
  m_modal_start = ivector(1,Nmodal+1);
  dt = params->dt;//mexPrintf("dt = %g\n",dt);

  static_delay = 0;
  for (modal = 1; modal <= Nmodal; modal++) {
    m_modal[modal] = params->ny_modal[modal]; 
  }
 
  m_modal_start[1] = 0;
  modal_delay[1] = params->delay[1];
  modal_delay_steps[1] = int(ceil(modal_delay[1]/dt));
  if (modal_delay[1] > static_delay) static_delay = modal_delay[1]; 
  for (modal = 2; modal <= Nmodal; modal++) {
    m_modal_start[modal] = m_modal[modal-1] + m_modal_start[modal-1];
    modal_delay[modal] = params->delay[modal];
    modal_delay_steps[modal] = int(ceil(modal_delay[modal]/dt));
    if (modal_delay[modal] > static_delay) static_delay = modal_delay[modal]; // EKF delay is max feedback delay
  }
 
  m = m_modal_start[Nmodal+1]  = m_modal_start[Nmodal] + m_modal[Nmodal]; // dimensionality of FB
  
  n = params->nx;
  p = params->nu;
  
  E = zeros(n,n);  
  A = dmatrix(1,n,1,n);
  tA = dmatrix(1,n,1,n);
  B = dmatrix(1,n,1,p);
  tB = dmatrix(1,p,1,n);
  K = dmatrix(1,n,1,m);
  tK = dmatrix(1,m,1,n);
  H = dmatrix(1,m,1,n);
  noise_type = params->noise_type; 
  Nq = (int) ceil(static_delay/dt) + 1;
  internal_delay = (Nq-1)*dt;
  Q = dmatrix(1,p,1,p);
  R = dmatrix(1,m,1,m);
  S = dmatrix(1,n,1,n);
  Qp = dmatrix(1,p,1,p);
  Rp = dmatrix(1,m,1,m);
  
  sigR = dvector(1,Nmodal);
  sigRp = dvector(1,Nmodal);
  
  mean_u = dvector(1,p);
  mean_y = dvector(1,m);
  mean_u_saved = dvector(1,p);
  mean_y_saved = dvector(1,m);
  
  update_ekf(params);
  security = eye(m);
  scalarmult(security,m,m,security_shift);
  
  kmax = 500;
}

/**********************************/
/**********************************/
bool is_symmetrical(double **A, int n)
{
  for (int i = 1; i<= n; i++)
    for (int j = 1; j <= i; j++)
      if (fabs(A[i][j] - A[j][i]) > EPSILON) return(false);
  
  return(true);
}


void calc_discrete_ekf_gains(double t, double *xhat, double *u, int *available_feedback)
{
  double **temp1, **temp2, **temp3, **temp4;
  double **itemp2;
  double **tH;
  static int counter = 0;
  double **Ebefore;
  int k, modal;
  
  tH = dmatrix(1,n,1,m);

  temp1 = dmatrix(1,n,1,m);
  temp2 = dmatrix(1,m,1,m);
  temp3 = dmatrix(1,n,1,n);
  temp4 = dmatrix(1,n,1,n);
  itemp2 = dmatrix(1,m,1,m);
  Ebefore = dmatrix(1,n,1,n);

  dhdx(xhat, u, t, H);
  
  copy(Ebefore,E,n);

  transpose(H,m,n,tH);
  matmult(n,E,n,tH,m,temp1);
  matmult(m,H,n,temp1,m,temp2); 
  matadd(temp2,Rp,temp2,m); 
  matadd(temp2,R,temp2,m);
  matadd(temp2,security,temp2,m); 

  inverse(temp2, itemp2, m); 
  
  matmult(n,temp1,m,itemp2,m,K);   
  if (alert) {
  mexmatdisplay(temp2,m,m,"temp2 "); // XXXX
  mexmatdisplay(itemp2,9,9,"itemp2 upper quadrant "); // XXXX
  mexmatdisplay(K,n,m,"K "); // XXXX
  mexmatdisplay(tH,n,m,"tH "); // XXXX

mexErrMsgTxt("Fuckup.");
  }
  eye(temp3, n);
  matmult(n,K,m,H,n,temp4);
  matsub(temp3,temp4,temp4,n);  // temp4 == I - KH
  matmult(n,temp4,n,E,n,temp3);
  copy(E,temp3,n); // E(+) = (I-KH)E(-)
  free_dmatrix(tH,1,n,1,m);
  free_dmatrix(Ebefore,1,n,1,n);
  free_dmatrix(temp1,1,n,1,m);
  free_dmatrix(temp2,1,m,1,m);
  free_dmatrix(temp3,1,n,1,n);
  free_dmatrix(temp4,1,n,1,n);
  free_dmatrix(itemp2,1,m,1,m);
} 




void init_discrete_ekf_variables(double **kxest, double **kustart, double **kystart)
{
    xhat = dvector(1,n);
    ku = dmatrix(1,Nq,1,p);
    kxhat = dmatrix(1,Nq,1,n);
    ky = dmatrix(1,Nq,1,m);


  copy(kxhat, kxest, Nq, n);
  copy(ku, kustart, Nq, p);
  copy(ky, kystart, Nq, m);

  veccopy(xhat, kxhat[Nq], n); 
  tindex = 1;

  first_time = 0;
}



void init_discrete_ekf_variables_static(double *xest, double *ustart, double *ystart, double **Estart)
{
  int i;
    xhat = dvector(1,n);
    ku = dmatrix(1,Nq,1,p);
    kxhat = dmatrix(1,Nq,1,n);
    ky = dmatrix(1,Nq,1,m);
    
  veccopy(xhat,xest,n);
  

  for ( i=1; i<=Nq; i++) veccopy(ku[i],ustart,p);
  for ( i=1; i<=Nq; i++) veccopy(kxhat[i],xhat,n); 
  for ( i=1; i<=Nq; i++) veccopy(ky[i],ystart,m); 
  
  if (Estart) copy(E, Estart, n);
  else zeros(E,n,n); 
  
  tindex = 1;

  first_time = 0;
}


/**********************************/
/**********************************/

void update_gaussian(int *available_feedback)
{
  double *diR = dvector(1,m); 
  int i, k , modal;   
    
  k = 1;
    
  //mexPrintf("sigR[1] = %f\n",sigR[1]);
  //mexPrintf("sigR[2] = %f\n",sigR[2]);
  //mexPrintf("sigR[3] = %f\n",sigR[3]);
  for (modal = 1; modal <= Nmodal; modal++) {
    if (available_feedback[modal])
      for (i = 0; i < m_modal[modal]; i++, k++) diR[k] = sigR[modal];
    else 
      for (i = 0; i < m_modal[modal]; i++, k++) diR[k] = SIG_NO_FEEDBACK;
  }
  diag(R,diR,m); 
  free_dvector(diR,1,m);
}






void update_sdn(double *u, double *y, int *available_feedback, double *mu_saved, double *my_saved, bool set_mean)
{
  static double *mu = zeros(p);
  double *du; 
  static double *my = zeros(m);
  double *dy;
  double *dy_modal;
  
  if (set_mean) {
    veccopy(mu, mu_saved, p); 
    veccopy(my, my_saved, m);
  }
  sliding_mean(u,p,mu,10);
  sliding_mean(y,m,my,10); 
  if (!set_mean) {
    veccopy(mu_saved, mu, p); 
    veccopy(my_saved, my, m);
  }
  
  //  vecdisplay(u,p,"u");

  du = apply(mu,p,sqr); 
  scalarmult(du, p, sigQp);
  diag(Qp,du,p);  

  dy = apply(my,m,sqr);
  for (int modal = 1; modal <= Nmodal; modal++) {
    dy_modal = vecextract(dy,m_modal_start[modal]+1,m_modal_start[modal+1]);
    if (available_feedback[modal]) 
        scalarmult(dy_modal, m_modal[modal], sigRp[modal]);//
    else {
        ones(dy_modal, m_modal[modal]);
        scalarmult(dy_modal, m_modal[modal], SIG_NO_FEEDBACK); 
    }
    vecassign(dy_modal, m_modal[modal], dy, m_modal_start[modal]+1);
    free_dvector(dy_modal,1,m_modal[modal]);
  }
   diag(Rp,dy, m);  
   
  if (isnan(Rp[m][m])) {mexvecdisplay(y,m,"y"); mexvecdisplay(my,m,"my"); mexvecdisplay(sigRp,Nmodal,"sigRp"); mexErrMsgTxt("SDN fuckup.");}

  free_dvector(dy,1,m);
  free_dvector(du,1,p); 
}







//  uodate noise covariance for motor & sensory signals
void update_poisson(double *u, double *y, int *available_feedback, double *mu_saved, double *my_saved, bool set_mean)
{
  static double *mu = zeros(p);
  static double *my = zeros(m);
  double *dy = dvector(1,m);
  double *dy_modal;
  
  if (set_mean) {
    veccopy(mu, mu_saved, p); // set mean from memorized value
    veccopy(my, my_saved, m);
  }
  sliding_mean(u,p,mu,10);
  sliding_mean(y,m,my,10); 
  if (!set_mean) {
    veccopy(mu_saved, mu, p); // save mean
    veccopy(my_saved, my, m);
  }

  diag(Qp, mu, p);  
  scalarmult(Qp, p,p, sigQp);

  for (int modal = 1; modal <= Nmodal; modal++) {
    dy_modal = vecextract(my,m_modal_start[modal]+1,m_modal_start[modal+1]);
    if (available_feedback[modal]) scalarmult(dy_modal, m_modal[modal], sigRp[modal]);
    else scalarmult(dy_modal, m_modal[modal], SIG_NO_FEEDBACK);
    vecassign(dy_modal, m_modal[modal], dy, m_modal_start[modal]+1);
    free_dvector(dy_modal,1,m_modal[modal]);
  }
  diag(Rp,dy, m);  
  
  free_dvector(dy,1,m);
}



void f_derivs(double t, double *x, double *dxdt)
{
  double *u, *u_old;

  u = circtable1(tu,tindex,t,Nq,p);
	
	if (t<dt) u_old = circtable1(tu,tindex,t,Nq,p); 
	else u_old = circtable1(tu,tindex,t-dt,Nq,p);
	

  f(x,u_old,u,t,dxdt,x,x);

  free_dvector(u,1,p);
  free_dvector(u_old,1,p);
}


void cov_derivs(double t, double *flat_P, double *flat_Pdot)
{
  double *u, *x, *u_old;
  double **P = dmatrix(1,n,1,n);
  double **Pdot = dmatrix(1,n,1,n);
  double **temp1 = dmatrix(1,n,1,n);
  double **temp2 = dmatrix(1,n,1,n);
  double **temp3 = dmatrix(1,n,1,p);
  double **temp4 = dmatrix(1,n,1,n);
  double **temp5 = dmatrix(1,p,1,p);
  double **temp6 = dmatrix(1,n,1,n);
  double **temp7 = dmatrix(1,m,1,m);
  double **temp8 = dmatrix(1,n,1,m);

  x = circtable1(txhat,tindex,t,Nq,n);
  u = circtable1(tu,tindex,t,Nq,p);
  
  if (t<dt) u_old = circtable1(tu,tindex,t,Nq,p); 
  else u_old = circtable1(tu,tindex,t-dt,Nq,p);
  dfdx(x,u,t,A);
  dfdu(x,u_old,u,t,B);
  transpose(A,n,n,tA);
  transpose(B,n,p,tB);
  transpose(K,n,m,tK);
  col2symmat2(flat_P,n,P);

  matmult(n,A,n,P,n,temp1);
  matmult(n,P,n,tA,n,temp2);
  matadd(temp1,temp2,temp1,n);

  matadd(Q,Qp,temp5,p);
  matmult(n,B,p,temp5,p,temp3);
  matmult(n,temp3,p,tB,n,temp4);
  
  matadd(R,Rp,temp7,m);
  matmult(n,K,m,temp7,m,temp8);
  matmult(n,temp8,m,tK,n,temp6);

  matadd(temp1,temp4,temp4,n);
  matsub(temp4,temp6,temp4,n);
  matadd(temp4,S,Pdot,n);

  symmat2col(Pdot,n,flat_Pdot);

  free_dvector(x,1,n);
  free_dvector(u,1,p);
  free_dvector(u_old,1,p);
  free_dmatrix(P,1,n,1,n);
  free_dmatrix(Pdot,1,n,1,n);
  free_dmatrix(temp1,1,n,1,n);
  free_dmatrix(temp2,1,n,1,n);
  free_dmatrix(temp3,1,n,1,p);
  free_dmatrix(temp4,1,n,1,n);
  free_dmatrix(temp5,1,p,1,p);
  free_dmatrix(temp6,1,n,1,n);
  free_dmatrix(temp7,1,m,1,m);
  free_dmatrix(temp8,1,n,1,m);
}


void ekf_propagate(double t, double *xhat, double *u, double T, double deltat, double *xhat_pred)
{
  int i, nok, nbad, N, n2=(n*(n+1))/2;
  double *xstart = dvector(1,n);
  double *xp = dvector(1,kmax);
  double **yp = dmatrix(1,kmax,1,n2);
  double *time = dvector(1,1);
  double sample_time;


  if (refresh_queues) {
    double *timeu = dvector(1,p+1);
    double *internal_u;

    time[1] = t;
    vecconcat(time,1,u,p,timeu);
    veccopy(tu[tindex],timeu,p+1);  // temporarily insert u in circular queue

    internal_u = circtable1(tu,tindex,internal_time,Nq,p+1);

    time[1] = internal_time;    
    vecconcat(time,1,internal_u,p,timeu);
    veccopy(tu[tindex],timeu,p+1);  // insert proper u in circular queue    

    free_dvector(timeu,1,p+1);
    free_dvector(internal_u,1,p);

  }

  veccopy(xstart,xhat,n);

  backwards_euler(xstart, n, internal_time-T, t, eps, keuler, f_derivs,&N,xp,yp);
  veccopy(xhat_pred, yp[N],n); // estimate at present time


  if (refresh_queues) {
    double **mxp = dmatrix(1,kmax,1,1);
    double *timex = dvector(1,n+1);
    double *flat_E = dvector(1,n2);
    double **ts = dmatrix(1,kmax,1,n+1);
    double *xhat_queued = dvector(1,n);


    vec2mat(xp,N,'c',mxp);
    matconcat(mxp,N,1,yp,N,n,ts);
    
    for ( i=0; i<Nq; i++ ) {
      sample_time = internal_time-i*deltat;
      table1(ts,sample_time,N,n, xhat_queued); // temporary estimate at next internal time step
      time[1] = sample_time;    
      vecconcat(time,1,xhat_queued,n,timex);
      veccopy(txhat[(tindex-i+Nq-1)%Nq+1],timex,n+1);
    }
    veccopy(xhat,xhat_queued,n);  // temporary estimate at last internal time step

    tindex = tindex%Nq + 1;
    refresh_queues = 0;
    
    symmat2col(E,n,flat_E);
    backwards_euler(flat_E, n2, internal_time-T, internal_time-T+deltat, eps, keuler, cov_derivs,&N,xp,yp);

    col2symmat2(yp[N],n,E);
    
    free_dmatrix(mxp,1,kmax,1,1);
    free_dvector(timex,1,n+1);
    free_dvector(flat_E,1,n2);
    free_dmatrix(ts,1,kmax,1,n+1);
    free_dvector(xhat_queued,1,n);
  }
 

  free_dvector(xstart,1,n);
  free_dvector(xp,1,kmax);
  free_dmatrix(yp,1,kmax,1,n2);
  free_dvector(time,1,1);
}




// extended Kalman filter
void ekf(double t, double *u, double *y, double *xest, int *consider_feedback, ekf_params *params, int reset)
{
  double *xhat_now, *yhat, *temp4, *temp5;
  double **I;
  double *time = dvector(1,1);
  double *timey = dvector(1,m+1);
  double *internal_y;
  double *u_old, *y_old, t_old;
  int circindex;
  int *available_feedback = ivector(1,Nmodal);

  if (first_time || reset) {
    if (!first_time) free_ekf_variables();
    init_ekf(params);
    init_ekf_variables_static(xest,u,y);
  }

  refresh_queues = (t>=internal_time+dt);

  if (refresh_queues) {


    internal_time += dt;
    time[1] = t;
    vecconcat(time,1,y,m,timey);
    veccopy(ty[tindex],timey,m+1);  // insert y in circular queue to get resampled value
    internal_y = circtable1(ty,tindex,internal_time,Nq,m+1);


    time[1] = internal_time;    
    vecconcat(time,1,internal_y,m,timey);
    veccopy(ty[tindex],timey,m+1);  // insert proper y in circular queue    
    free_dvector(internal_y,1,m);
  

    yhat = dvector(1,m);
    temp4 = dvector(1,m);
    temp5 = dvector(1,n);          

    xhat_now = dvector(1,n);

    
    
    
    
    for (int i=0; i < Nq-1; i++) {
      circindex = (tindex+i)%Nq+1;
      
      u_old = vecextract(tu[circindex],2,p+1); // 2nd oldest value
      y_old = vecextract(ty[circindex],2,m+1); // oldest value
      vecextract(txhat[circindex],2,n+1,xhat); // 2nd oldest value
      t_old = internal_time +i*dt - internal_delay;

      for (int modal=1; modal <= Nmodal; modal++) available_feedback[modal] = int(consider_feedback[modal] && (Nq-i > modal_delay_steps[modal]));
        
      switch (noise_type) {      
      case POISSON: update_poisson(u_old, y_old, available_feedback, mean_u, mean_y, false);
      case GAUSSIAN: update_gaussian(available_feedback); break;
      case MIXED:
      case SDN: update_sdn(u_old, y_old, available_feedback, mean_u, mean_y, false);
      }
      
      calc_ekf_gains(t_old,xhat,u_old); 
        
    h(xhat, u_old, t_old, yhat); // h(xhat{t}(-))
    vecsub(y_old,yhat,temp4,m);
    matvecmult(n,K,m,temp4,temp5); // error term
    
      vecadd(xhat,temp5,xhat,n); 

      ekf_propagate(t, xhat, u, internal_delay, dt, xhat_now); 
    } 
        
    free_dvector(yhat,1,m);
    free_dvector(y_old,1,m);
    free_dvector(u_old,1,p);
    free_dvector(temp4,1,m);
    free_dvector(temp5,1,n);    
  }

  veccopy(xest,xhat_now,n);
  
  free_dvector(xhat_now,1,n);
  free_ivector(available_feedback,1,Nmodal);
}


/*******************/

void discrete_ekf(double t, double *u, double *y, double *xest, int *consider_feedback, double *HS, double *CS, ekf_params *params, int reset) 
{
  double *yhat, *sensory_error, *xhat_corr;
  double *u_old, *u_older, *y_old, t_old;
  double **temp1 = dmatrix(1,n,1,n);
  double **temp2 = dmatrix(1,n,1,n);
  double **temp3 = dmatrix(1,n,1,p);
  double **temp4 = dmatrix(1,n,1,n);
  double **temp5 = dmatrix(1,p,1,p);
  int new_tindex, circindex, ind;
  int *available_feedback = ivector(1,Nmodal);
  static int counter = 0;
  double **Esaved; 
  bool DBG = false; 
  
    u_older = dvector(1,p);

  if (first_time || reset) {
   if (!first_time) free_discrete_ekf_variables();
    init_ekf(params);
    init_discrete_ekf_variables_static(xest,u,y);
    set_LSTM_states(t,dt,HS,CS,0);
  }

  veccopy(u_older,ku[tindex],p);
  
  veccopy(ky[tindex],y,m);   
  veccopy(ku[tindex],u,p); 
  
  y_old = dvector(1,m);
  u_old = dvector(1,p);
  yhat = dvector(1,m);
  sensory_error = dvector(1,m);
  xhat_corr = dvector(1,n);
  t_old = t - (Nq-1)*dt;
  new_tindex = tindex%Nq+1; 
  xhat = kxhat[new_tindex]; 

  Esaved = dmatrix(1,n,1,n);
  
  veccopy(mean_u, mean_u_saved, p); 
  veccopy(mean_y, mean_y_saved, m);

  circindex = new_tindex;
  for (int i=0; i < Nq-1; i++) {
    veccopy(u_older,u_old,p);
    veccopy(u_old,ku[circindex],p); 
    veccopy(y_old,ky[circindex],m); 
    
    h(xhat, u_old, t_old, yhat); 
    vecsub(y_old,yhat,sensory_error,m); 
 
    for (int modal=1; modal <= Nmodal; modal++) {
      available_feedback[modal] = int(consider_feedback[modal] && (Nq-i > modal_delay_steps[modal]));
      if (!available_feedback[modal])
        for (int j = m_modal_start[modal]+1; j <= m_modal_start[modal+1]; j++) sensory_error[j] = 0; 
    }
 
    switch (noise_type) {     
      case POISSON: update_poisson(u_old, y_old, available_feedback, mean_u, mean_y, i==0); break;
      case GAUSSIAN: update_gaussian(available_feedback); break;
      case MIXED: 
      case SDN: update_sdn(u_old, y_old, available_feedback, mean_u, mean_y, i==0); 
    }
    
    mexvecdisplay(xhat,6,"xhat"); 
    calc_discrete_ekf_gains(t_old,xhat,u_old, available_feedback); 
    matvecmult(n,K,m,sensory_error,xhat_corr); 
    vecadd(xhat,xhat_corr,xhat,n); 
    dfdx(xhat,u_old,t_old,A);  
    dfdu(xhat,u_older,u_old,t_old,B);
    if (DBG) mexmatdisplay(A,n,n,"A");  
    transpose(A,n,n,tA);
    transpose(B,n,p,tB);
    matmult(n,A,n,E,n,temp1);
    matmult(n,temp1,n,tA,n,temp2);  
    matadd(Q,Qp,temp5,p);
    matmult(n,B,p,temp5,p,temp3);
    matmult(n,temp3,p,tB,n,temp4);   if (DBG) mexmatdisplay(temp4,n,n,"temp4");
    matadd(temp2,temp4,temp4,n);
    matadd(temp4,S,E,n);
    f(xhat, u_older, u_old, t_old, xhat, HS, CS);  
    t_old = t-dt*(Nq-i-2); 
    circindex = circindex%Nq+1;
    
    veccopy(kxhat[circindex],xhat,n); 
    

    if (i == 0) {
      copy(Esaved,E,n); 
     
      veccopy(mean_u_saved, mean_u, p); 
      veccopy(mean_y_saved, mean_y, m);
    }

  } 
  veccopy(xest,xhat,n);  
  tindex = new_tindex;


  for (int cc = 2; cc <= 6; cc++) { 
    if (xest[cc] > 1) mexErrMsgTxt("WTF K+?");
    if (xest[cc] < 0) mexErrMsgTxt("WTF K?");
  }

  copy(E,Esaved,n); 
  free_dmatrix(Esaved,1,n,1,n);
  free_dmatrix(temp1,1,n,1,n);
  free_dmatrix(temp2,1,n,1,n);
  free_dmatrix(temp3,1,n,1,p);
  free_dmatrix(temp4,1,n,1,n);
  free_dmatrix(temp5,1,p,1,p);

  free_ivector(available_feedback,1,Nmodal);
}


void discrete_ekf(double t, double *u, double *y, double *xest, int *consider_feedback, double *yest, double **outE, ekf_params *params, int reset) 
{
  double *yhat, *sensory_error, *xhat_corr, *new_xhat;
  double *u_old, *u_older, *y_old, t_old;
  double **temp1 = dmatrix(1,n,1,n);
  double **temp2 = dmatrix(1,n,1,n);
  double **temp3 = dmatrix(1,n,1,p);
  double **temp4 = dmatrix(1,n,1,n);
  double **temp5 = dmatrix(1,p,1,p);
  int new_tindex, circindex, ind;
  int *available_feedback = ivector(1,Nmodal);
  static int counter = 0;
  double **Esaved;
  double *HS = zeros(LSTM_SIZE);
  double *CS = zeros(LSTM_SIZE);
  bool DBG = false;
  
  u_older = dvector(1,p); 

  if (first_time || reset) {
   if (!first_time) free_discrete_ekf_variables();
    init_ekf(params);
    init_discrete_ekf_variables_static(xest,u,y);
    set_LSTM_states(t,dt,HS,CS,0);
  }

  veccopy(u_older,ku[tindex],p);
  
  veccopy(ky[tindex],y,m);     
  veccopy(ku[tindex],u,p); 
  
  y_old = dvector(1,m);
  u_old = dvector(1,p);
  yhat = dvector(1,m);
  sensory_error = dvector(1,m);
  xhat_corr = dvector(1,n);
  new_xhat = dvector(1,n);
  
  t_old = t - (Nq-1)*dt;
  new_tindex = tindex%Nq+1; 
  xhat = kxhat[new_tindex]; 
  
  Esaved = dmatrix(1,n,1,n);
  
  veccopy(mean_u, mean_u_saved, p); 
  veccopy(mean_y, mean_y_saved, m);

  circindex = new_tindex;
 
  for (int i=0; i < Nq-1; i++) {
    veccopy(u_older,u_old,p);
    veccopy(u_old,ku[circindex],p); 
    veccopy(y_old,ky[circindex],m); 
    h(xhat, u_old, t_old, yhat); 
    vecsub(y_old,yhat,sensory_error,m); 
    
    for (int modal=1; modal <= Nmodal; modal++) {
      available_feedback[modal] = int(consider_feedback[modal] && (Nq-i > modal_delay_steps[modal]));
      if (!available_feedback[modal])
        for (int j = m_modal_start[modal]+1; j <= m_modal_start[modal+1]; j++) sensory_error[j] = 0; 
    }
  
    switch (noise_type) {     
      case POISSON: update_poisson(u_old, y_old, available_feedback, mean_u, mean_y, i==0); 
      case GAUSSIAN: update_gaussian(available_feedback); 
      case MIXED:
      case SDN: update_sdn(u_old, y_old, available_feedback, mean_u, mean_y, i==0); 
    }
   
    calc_discrete_ekf_gains(t_old,xhat,u_old, available_feedback); 
    matvecmult(n,K,m,sensory_error,xhat_corr); 
    vecadd(xhat,xhat_corr,xhat,n); 
    dfdx(xhat,u_old,t_old,A);  
    dfdu(xhat,u_older,u_old,t_old,B);
    if (DBG) mexmatdisplay(A,n,n,"A");  
   
    transpose(A,n,n,tA);
    transpose(B,n,p,tB);
    matmult(n,A,n,E,n,temp1);
    matmult(n,temp1,n,tA,n,temp2);  
    matadd(Q,Qp,temp5,p);
    matmult(n,B,p,temp5,p,temp3);
    matmult(n,temp3,p,tB,n,temp4);   if (DBG) mexmatdisplay(temp4,n,n,"temp4");
    matadd(temp2,temp4,temp4,n);
    matadd(temp4,S,E,n); 
    f(xhat, u_older, u_old, t_old, new_xhat, HS, CS);   
    veccopy(xhat,new_xhat,n); 
    t_old = t-dt*(Nq-i-2);  
    
    circindex = circindex%Nq+1;
    veccopy(kxhat[circindex],xhat,n); 
    

    if (i == 0) {
      copy(Esaved,E,n);
      veccopy(mean_u_saved, mean_u, p);
      veccopy(mean_y_saved, mean_y, m);
    }
  }
	
  veccopy(xest,xhat,n);  
  veccopy(yest,yhat,m);
  copy(outE,Esaved,n); 
  tindex = new_tindex;

  copy(E,Esaved,n); 
    
  free_dmatrix(Esaved,1,n,1,n);
  free_dmatrix(temp1,1,n,1,n);
  free_dmatrix(temp2,1,n,1,n);
  free_dmatrix(temp3,1,n,1,p);
  free_dmatrix(temp4,1,n,1,n);
  free_dmatrix(temp5,1,p,1,p);
  free_ivector(available_feedback,1,Nmodal);
  free_dvector(xhat_corr,1,n);
  free_dvector(new_xhat,1,n);
  free_dvector(HS,1,LSTM_SIZE);
  free_dvector(CS,1,LSTM_SIZE);
}

/*******************/


void free_ekf_variables()
{
  if (modal_delay) { 
     free_dvector(modal_delay,1,Nmodal);
     free_ivector(modal_delay_steps,1,Nmodal);
     free_ivector(m_modal,1,Nmodal);
     free_ivector(m_modal_start,1,Nmodal+1);
  }

  if (A) {
    free_dmatrix(E,1,n,1,n);
    free_dmatrix(A,1,n,1,n); 
    free_dmatrix(tA,1,n,1,n);
    free_dmatrix(B,1,n,1,p); 
    free_dmatrix(tB,1,p,1,n);
    free_dmatrix(K,1,n,1,m);
    free_dmatrix(tK,1,m,1,n);
  
    free_dvector(sigR,1,Nmodal);
    free_dvector(sigRp,1,Nmodal);
  
    free_dmatrix(Q,1,p,1,p);
    free_dmatrix(R,1,m,1,m);
    free_dmatrix(S,1,n,1,n);
    free_dmatrix(Qp,1,p,1,p);
    free_dmatrix(Rp,1,m,1,m);
    free_dmatrix(H,1,m,1,n);
    free_dmatrix(security,1,m,1,m);
  }
  
  if (txhat) {
    free_dmatrix(txhat,1,Nq,1,n+1);
    free_dmatrix(tu,1,Nq,1,p+1);
    free_dmatrix(ty,1,Nq,1,m+1);
    free_dvector(xhat,1,n);
    
    free_dvector(mean_u,1,p);
    free_dvector(mean_y,1,m);
    free_dvector(mean_u_saved,1,p);
    free_dvector(mean_y_saved,1,m);
  }
}


void free_discrete_ekf_variables()
{
  if (modal_delay) { 
     free_dvector(modal_delay,1,Nmodal);
     free_ivector(modal_delay_steps,1,Nmodal);
     free_ivector(m_modal,1,Nmodal);
     free_ivector(m_modal_start,1,Nmodal+1);
  }

  if (A) { // minimal safety
    free_dmatrix(E,1,n,1,n);
    free_dmatrix(A,1,n,1,n); 
    free_dmatrix(tA,1,n,1,n);
    free_dmatrix(B,1,n,1,p); 
    free_dmatrix(tB,1,p,1,n);
    free_dmatrix(K,1,n,1,m);
    free_dmatrix(tK,1,m,1,n);

    free_dvector(sigR,1,Nmodal);
    free_dvector(sigRp,1,Nmodal);

    free_dmatrix(Q,1,p,1,p);
    free_dmatrix(R,1,m,1,m);
    free_dmatrix(S,1,n,1,n);
    free_dmatrix(Qp,1,p,1,p);
    free_dmatrix(Rp,1,m,1,m);
    free_dmatrix(H,1,m,1,n);
    free_dmatrix(security,1,m,1,m);
  }
  
  if (kxhat) {
    free_dmatrix(kxhat,1,Nq,1,n);
    free_dmatrix(ku,1,Nq,1,p);
    free_dmatrix(ky,1,Nq,1,m);
    free_dvector(xhat,1,n);
    
    free_dvector(mean_u,1,p);
    free_dvector(mean_y,1,m);
    free_dvector(mean_u_saved,1,p);
    free_dvector(mean_y_saved,1,m);    
  }
}


