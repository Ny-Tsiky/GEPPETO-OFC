#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mex.h"
#include "mylinpack.h"
#include "nr.h"
#include "nrutil.h"  
#include "kalman8.h"
#include "plantfuns_via.h"
#include "plantfeedback.h"

#undef DB
#define DB mexPrintf("Got to line %d in mex file %s.\n",__LINE__,__FILE__);

extern bool feedback_initialized;

bool resetPlantLSTM;
bool resetPlantLSTMderivatives;

static int ns, np, nc;
extern int ny, nyp, NF, Nproprio, NT, Nmodal, LSTM_SIZE;
extern int *ny_modal;
int *modal_start;
int *consider_feedback;
noise_params noisepar;

const double u0[] = {50,50,44,51,95,18,90,63};

int N;
int N_VTSECTIONS = 44; 

double **u, **s, **sest, **clean_y, **y, **yest, *nu, **E;
double *HS_Kalman, *CS_Kalman, *HS_Plant, *CS_Plant;

FILE *fout;

void mexmatdisplay(double **M, int m, int n, const char *matname);
void mexvecdisplay(double *M, int m, const char *matname);

void alloc_estim_variables()
{
//  y = dvector(1,ny);
}

void free_estim_variables()
{
//  free_dvector(y,1,ny);
}


void copy_columnwise(double **A, double *p, int m, int n)
{
  int i,j,k = 0;
  
  for ( j=1; j<=n; j++ )
    for ( i=1; i<=m; i++ )
      A[i][j] = p[k++];
}


void myExitFcn() {
  mexPrintf("MEX-file is being unloaded, freeing memory\n\n");
  free_plant_variables(); 
  free_feedback_variables();
  free_dvector(noisepar.var_s_sdn,1,Nmodal);
  free_dvector(noisepar.var_s_add,1,Nmodal);
  free_ivector(consider_feedback,1,Nmodal);
  free_ivector(modal_start,1,Nmodal);
	}

 void Bark2Hz(double *bark_y, double *y, int NF){
	 for (int i=1; i<=NF; i++) {
		if (bark_y[i] < 2) bark_y[i] = (bark_y[i]-0.3)/0.85;
		if (bark_y[i] > 20.1) bark_y[i] = (bark_y[i]+4.422)/1.22;
		y[i] = 1960*((bark_y[i]+0.53)/(26.28-bark_y[i]));
	 }
  }

void state_estimate(double t, double *u, double *s, double sqdt, double *delay, double *SDN_Snoisestd, double *Gaussian_Snoisestd, double *sest, double *clean_y, double *y, double *yest, double **E, double *y_pert, double t_pert_on, double t_pert_full)
{
  int k, modal;
  int huhu=0;
  feedback(s, u, clean_y); 
  if (isnan(clean_y[ny])) {mexvecdisplay(y,ny,"state_estimate: y");}
	for ( modal = 1; modal <= Nmodal; modal++ ) {
    for (k = modal_start[modal]+1; k <= modal_start[modal]+ny_modal[modal]; k++) {
       y[k] = clean_y[k]*(1+SDN_Snoisestd[modal]*randn()*sqdt) + Gaussian_Snoisestd[modal]*randn()*sqdt; 
    }
  }
  
  if (t > t_pert_on) vecadd(y, y_pert, y, ny); 
  
   for (k = 1; k<=Nmodal; k++) {
     consider_feedback[k] = int(t>delay[k]);
   }

  	discrete_ekf(t, u, y, sest, consider_feedback, yest, E);
}

void noisy_forward( double **u, double *s0, double *sest0, double **s, double **sest, double **clean_y, double **y, double **yest, double *tspan, double dt, double *delay, noise_params *noise, double **state_cov)
{
  int i,k;
  double **clean_u, **old_u;
  double sqdt = sqrt(dt);
  double SDN_Unoise_std = sqrt(noise->var_u_sdn);
  double Gaussian_Unoise_std = sqrt(noise->var_u_add);
  double *SDN_Snoise_std =  dvector(1,Nmodal);
  double *Gaussian_Snoise_std = dvector(1,Nmodal);
  double Gaussian_Xnoise = sqrt(noise->var_x_add);
  clean_u = dmatrix(1,N,1,nc);
  veccopy(s[1],s0,ns);
  veccopy(sest[1],sest0,ns);
  
  for ( i=1; i<=N; i++ ) {
    veccopy(clean_u[i],u[i],nc);
    for ( k=1; k<=nc; k++ ) {
			u[i][k] += pospart((u0[k]-u[i][k]))*SDN_Unoise_std*randn()*sqdt + Gaussian_Unoise_std*randn()*sqdt;
    }
  }
     
  for (int modal = 1; modal <= Nmodal; modal++) {
    SDN_Snoise_std[modal] = sqrt(noise->var_s_sdn[modal]);
    Gaussian_Snoise_std[modal] = sqrt(noise->var_s_add[modal]);
  }
  for ( i=1; i<=N; i++ ) {
    if (i == 1) discrete_forward(u[1],u[i],s[i],tspan[i],dt,s[i+1],1, HS_Plant, CS_Plant,1); 
    else {
        discrete_forward(u[i-1],u[i],s[i],tspan[i],dt,s[i+1],1, HS_Plant, CS_Plant,1); 
    }
    state_estimate(tspan[i], clean_u[i], s[i], sqdt, delay, SDN_Snoise_std, Gaussian_Snoise_std, sest[i+1], clean_y[i], y[i], yest[i], state_cov, noise->y_perturb, noise->t_onset, noise->t_full);
	}
    free_dmatrix(clean_u,1,N,1,nc);
	free_dvector(SDN_Snoise_std,1,Nmodal);
  	free_dvector(Gaussian_Snoise_std,1,Nmodal);
  }


void alloc_optim_variables()
{
  u = dmatrix(1,N,1,nc);
  nu = dvector(1,np);
  s = dmatrix(1,N+1,1,ns);
  sest = dmatrix(1,N+1,1,ns);
  clean_y = dmatrix(1,N,1,ny);
  y = dmatrix(1,N,1,ny);
  yest = dmatrix(1,N,1,ny);
  E = dmatrix(1,ns,1,ns);
  
  HS_Kalman = dvector(1,LSTM_SIZE);
  CS_Kalman = dvector(1,LSTM_SIZE);
  HS_Plant = dvector(1,LSTM_SIZE);
  CS_Plant = dvector(1,LSTM_SIZE);
}

void free_optim_variables()
{
  free_dmatrix(u,1,N,1,nc);
  free_dvector(nu,1,np);
  free_dmatrix(s,1,N+1,1,ns);
  free_dmatrix(sest,1,N+1,1,ns);
  free_dmatrix(clean_y,1,N,1,ny);
  free_dmatrix(y,1,N,1,ny);
  free_dmatrix(yest,1,N,1,ny);
  free_dmatrix(E,1,ns,1,ns);  
  
  free_dvector(CS_Kalman,1,LSTM_SIZE);
  free_dvector(HS_Kalman,1,LSTM_SIZE);
  free_dvector(CS_Plant,1,LSTM_SIZE);
  free_dvector(HS_Plant,1,LSTM_SIZE);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *s0, *sest0, *tspan, *p, *perturb;
  double *t, *f, *ls, *lsest, *lcy, *ly, *lyest, *lu, *lE;
  double *out_HS_Kalman, *out_CS_Kalman, *out_HS_Plant, *out_CS_Plant, *in_HS_Kalman, *in_CS_Kalman, *in_HS_Plant, *in_CS_Plant, *Eest0;
  double dt;  
  double *delay;
  int sz=8;
  static bool first_time = true;
    
  int i,j,k;
  mxArray *p_goalseq, *p_goalcodeseq;
  int Na, n;
  bool resetKalman;

  mexAtExit(myExitFcn);
    
  if (first_time) {
    mexinit(&ns,&np,&nc);
    if (!feedback_initialized) init_feedback();
    noisepar.var_s_sdn = dvector(1,Nmodal);
    noisepar.var_s_add = dvector(1,Nmodal);
    noisepar.y_perturb = zeros(ny);
    consider_feedback = ivector(1,Nmodal); for (i = 1; i<=Nmodal; i++) consider_feedback[i] = 1; 
    modal_start = ivector(1,Nmodal);
    modal_start[1] = 0;
	  for (i = 2; i <= Nmodal; i++) modal_start[i] = modal_start[i-1] + ny_modal[i-1];
  }  
    
  if (nrhs != 14) mexErrMsgTxt("Wrong number of parameters for simnoisy_discrete8");
  int c = 0;
  p = mxGetPr(prhs[c])-1;
  N = mxGetM(prhs[c++]);
  s0 = mxGetPr(prhs[c++])-1;
  sest0 = mxGetPr(prhs[c++])-1;
  Eest0 = mxGetPr(prhs[c++])-1;
  tspan = mxGetPr(prhs[c++])-1;
  perturb = mxGetPr(prhs[c])-1;
  memcpy(noisepar.y_perturb+1, mxGetData(mxGetField(prhs[c],0,"audiovector")),sz*NF);
  memcpy(&(noisepar.t_onset),  mxGetData(mxGetField(prhs[c],0,"t_onset")),sz);
  memcpy(&(noisepar.t_full), mxGetData(mxGetField(prhs[c++],0,"t_full")),sz);
  delay = mxGetPr(prhs[c++])-1;  
  memcpy(&(noisepar.var_u_sdn),mxGetData(mxGetField(prhs[c],0,"var_u_sdn")),sz);
  memcpy(&(noisepar.var_u_add),mxGetData(mxGetField(prhs[c],0,"var_u_add")),sz);
  memcpy(noisepar.var_s_sdn+1, mxGetData(mxGetField(prhs[c],0,"var_s_sdn")),sz*Nmodal);
  memcpy(noisepar.var_s_add+1, mxGetData(mxGetField(prhs[c],0,"var_s_add")),sz*Nmodal);
  memcpy(&(noisepar.var_x_add),mxGetData(mxGetField(prhs[c++],0,"var_x_add")),sz);
  resetKalman = (bool)*mxGetLogicals(prhs[c++]);
  resetPlantLSTM = resetPlantLSTMderivatives = (bool)*mxGetLogicals(prhs[c++]);
  in_HS_Kalman = mxGetPr(prhs[c++])-1;
  in_CS_Kalman = mxGetPr(prhs[c++])-1;
  in_HS_Plant =  mxGetPr(prhs[c++])-1;
  in_CS_Plant =  mxGetPr(prhs[c])-1;

   p_goalseq = mexGetVariable("global", "goalseq");
   if (!p_goalseq) {mexErrMsgTxt("goalseq pas trouve\n"); return;}
   p_goalcodeseq = mexGetVariable("global", "goalcodeseq");
   if (!p_goalcodeseq) {mexErrMsgTxt("goalcodeseq pas trouve\n"); return;}
   int Ngoals = mxGetN(p_goalseq);
   char *goalcodeseq = (char*)malloc((Ngoals+1)*sizeof(char));
   mxGetString(p_goalcodeseq, goalcodeseq, Ngoals+1);
   init_goals(mxGetPr(p_goalseq), goalcodeseq, mxGetM(p_goalseq), Ngoals); 
   free(goalcodeseq);
    
  n = mxGetN(prhs[1]);
  if ( n != ns ) {
    n = mxGetM(prhs[1]);
    if ( n != ns ) mexErrMsgTxt("Wrong initial state vector size");
  }

  n = mxGetN(prhs[4]);
  if ( n != N ) {
    n = mxGetM(prhs[4]);
    if ( n != N ) mexErrMsgTxt("Wrong time vector size (inconsistent with command matrix)");
  }

  if ( N <= 2 ) mexErrMsgTxt("Too few time steps");
  dt = (tspan[2] - tspan[1]);
    
	if (resetPlantLSTM) set_LSTM_states(tspan[1], dt, in_HS_Plant, in_CS_Plant, 1);
    
  alloc_optim_variables();
    
	 if (first_time) {   
        init_kalman_filter(s0,&noisepar,dt,delay,DISCRETE_EKF);    
		first_time = false;
  }
    
  col2mat(p,N,nc,u);
	if (resetKalman) {
		double *ustart = u[1];
		double *ystart = dvector(1,ny);
		feedback(s0,ustart,ystart);
    col2mat(Eest0,ns,ns,E);
  	init_discrete_ekf_variables_static(s0, ustart, ystart, E);
    
    update_kalman_noises(&noisepar, DISCRETE_EKF); 
    set_LSTM_states(tspan[1], dt, in_HS_Kalman, in_CS_Kalman, 0);
		free_dvector(ystart,1,ny);
	}

  plhs[0] = mxCreateDoubleMatrix(N+1,ns,mxREAL);
  ls = mxGetPr(plhs[0]);

  if (nlhs>1) {
    plhs[1] = mxCreateDoubleMatrix(N+1,ns,mxREAL);
    lsest = mxGetPr(plhs[1]);
  }

  if (nlhs>2) {
    plhs[2] = mxCreateDoubleMatrix(N,ny,mxREAL);
    lcy = mxGetPr(plhs[2]);
  }

  if (nlhs>3) {
    plhs[3] = mxCreateDoubleMatrix(N,ny,mxREAL);
    ly = mxGetPr(plhs[3]);
  }

  if (nlhs>4) {
    plhs[4] = mxCreateDoubleMatrix(N,ny,mxREAL);
    lyest = mxGetPr(plhs[4]);
  }

  if (nlhs>5) {
    plhs[5] = mxCreateDoubleMatrix(N,nc,mxREAL);
    lu = mxGetPr(plhs[5]);
  }

	if (nlhs>6) {
    plhs[6] = mxCreateDoubleMatrix(ns,ns,mxREAL);
    lE = mxGetPr(plhs[6]);
  }

	if (nlhs>7) {
    plhs[7] = mxCreateDoubleMatrix(1,LSTM_SIZE,mxREAL);
    out_HS_Kalman = mxGetPr(plhs[7]);
  }

	if (nlhs>8) {
    plhs[8] = mxCreateDoubleMatrix(1,LSTM_SIZE,mxREAL);
    out_CS_Kalman = mxGetPr(plhs[8]);
  }

  noisy_forward(u,s0,sest0,s,sest,clean_y,y,yest,tspan,dt,delay,&noisepar,E);

  k = 0;
  for ( j=1; j<=ns; j++ )
    for ( i=1; i<=N+1; i++ ) {
      ls[k] = s[i][j];
      if (nlhs > 1) lsest[k] = sest[i][j];
      k++;
    }
  if (nlhs > 2) {
    k = 0;
    for ( j=1; j<=ny; j++ )
      for ( i=1; i<=N; i++ )
        lcy[k++] = clean_y[i][j]; 
  }

  if (nlhs > 3) {
    k = 0;
    for ( j=1; j<=ny; j++ )
      for ( i=1; i<=N; i++ )
        ly[k++] = y[i][j];   
  }

  if (nlhs > 4) {
    k = 0;
    for ( j=1; j<=ny; j++ )
      for ( i=1; i<=N; i++ )
        lyest[k++] = yest[i][j]; 
  }
  
  if (nlhs > 5) {
    k = 0;
    for ( j=1; j<=nc; j++ )
      for ( i=1; i<=N; i++ )
        lu[k++] = u[i][j]; 
  }

	if (nlhs > 6) {
    k = 0;
    for ( j=1; j<=ns; j++ )
      for ( i=1; i<=ns; i++ )
        lE[k++] = E[i][j]; 
  }

	if (nlhs > 7) {
    k = 0;
    for ( j=1; j<=LSTM_SIZE; j++ )
        out_HS_Kalman[k++] = HS_Kalman[j]; 
  }

	if (nlhs > 8) {
    k = 0;
    for ( j=1; j<=LSTM_SIZE; j++ )
        out_CS_Kalman[k++] = CS_Kalman[j]; 
  }

  free_optim_variables();

}
