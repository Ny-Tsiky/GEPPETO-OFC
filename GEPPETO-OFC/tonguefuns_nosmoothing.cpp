#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <cstring>
#include "mylinpack.h"
#include "nr.h"
#include "nrutil.h"   
#include "unistd.h" 
#include "plantfeedback.h" 
#include "tongue_geom_audio_tactile_autoenc.h"
#include "autoencoder.h"
#include "my_lstmForward_c11.h"
#include "goal_geometry.h"


#ifndef M_PI
const double M_PI = 3.14159265358979;
#endif


#ifndef NO_MATLAB
#include "mex.h" 
#undef DB
#define DB mexPrintf("Got to line %d in mex file %s.\n",__LINE__,__FILE__);

#endif

#define SMOOTH_POS_CONSTRAINT     

const double VEL_CONSTR_SCALING  = 1e-9; 
const double DVEL_CONSTR_SCALING = 1e-9; 
const double PRECISION_TRADEOFF_SCALING = 1e-2; 
const double AUDITORY_GOAL_SCALING = 0.98; 
const double TACTILE_GOAL_SCALING = 0.01; 
const double PROPRIO_GOAL_SCALING = 0.01; 
const double MIN_TOUCH_LEVEL = 0.1;
const double MAX_TACT_ERROR = 10; 
const double SOFT_POS_BOUNDARY = 1e-2;

const double DT_EPS = 1e-3;

int  nc, ns, np, nyp, Nsvd;
extern int ny;
extern const int NF;
extern const int NT; 
const int NTe = 3; 
extern const int Nproprio; 
int Ngoals;
int Ngoal_modalities;


static double **mean_goal;
static double **cov_goal;
static double *Tg, **Wg;
static char *Cg;

static double **tu, **ts;
double *Fguess1 = NULL; 

extern bool autoenc_initialized;
extern bool resetLSTM;
extern bool resetLSTMderivatives;
extern bool saved;


int nn, nh1, nh2;
double *MIN_NETINPUT,  *MAX_NETINPUT;
double *net_input; 
double *dmap_minmax;
double **Jnet; 

int LSTM_SIZE = 20; 
double **CS_actual, **HS_actual, **CS_kalman, **HS_kalman;
double *saved_CS_kalman, *saved_HS_kalman, *saved_CS_actual, *saved_HS_actual;
int MAX_STEPS = 50000; 

double **svdmatrix, *averagetongue; 
static double invtau_a, invtau_e;
const double u0[] = {50,50,44,51,95,18,90,63};
static double *v, *a;
const double DELTA = 1e-2; 
static double *sdeltas, *scopy; 
const double Kext = 1; // TBD
const double Bext = 1; // TBD

double T;
const double tau_a=0.005, tau_e=0.005;   // s

static double *th_goal;

const int D = 2; 
const double sa = 10.0;
const double re = 0.0;
const double invr = 1e-6;

extern const double TONGUE_MARGIN; //= 1; // mm
const double TONGUEVEL_MARGIN = 1; // mm/s

//**************************************************
// ODE params

extern int N_ODEsteps;
extern double eps; //not used

//**************************************************

double global_sens_cost, global_tact_cost, global_neuro_cost;

//**************************************************


void mexmatdisplay(double **M, int m, int n, const char *matname)
{
  int i, j;

  mexPrintf("%s (%d x %d):\n",matname,m,n);
  for (i=1; i<=m; i++) {
    for (j=1; j<=n; j++) mexPrintf("%3.6f\t",M[i][j]);
    mexPrintf("\n");
  }
  mexPrintf("\n");
}

void mexvecdisplay(double *M, int m, const char *matname)
{
  int i;

  mexPrintf("%s (%d):\n",matname,m);
  for (i=1; i<=m; i++) {
    mexPrintf("%3.6f\t",M[i]);
    // mexPrintf("\n");
  }
  mexPrintf("\n");
}

void reshape(double *v, double **M, int m, int n)
{
  int i, j, k;

  k = 1;
  for (j=1; j<=n; j++)
    for (i=1; i<=m; i++)
      M[i][j] = v[k++];
}

const double SMOOTH_EPS = 1e-1;
static double SMOOTH_THRESH = log(exp(SMOOTH_EPS)-1);

double gfun(double x, double s=1, double re=0)
{
  double y;

  re -= SMOOTH_THRESH/s;
  y = s*(x-re);
  if (y<700) return(log(1+exp(y))/s);
  else return(x-re);
}

double gfunderiv(double x, double s=1, double re=0)
{
  double y;

  re -= SMOOTH_THRESH/s;
  y = s*(x-re);
  if (y<700) return(exp(y)/(1+exp(y)));
  else return(1);
}

double gfundderiv(double x, double s=1, double re=0)
{
  double y, z;

  re -= SMOOTH_THRESH/s;
  y = s*(x-re);
  z = exp(y);
  if (y<700) return(z/sqr(1+z)*s);
  else return(0);
}


//**************************************************
//**************************************************
double softabs(double x, double zeroshift=.1)
{
  return((sqr(x)+fabs(x)+zeroshift)/(fabs(x)+zeroshift)-1+zeroshift);
}
double softabs_deriv(double x, double zeroshift=.1)
{
  return((x*fabs(x)+2*x*zeroshift)/sqr(fabs(x)+zeroshift));
}


double softpospart(double x, double zeroshift)
{
  double y = (sqr(x)+fabs(x)+zeroshift*2)/(fabs(x)+zeroshift*2)-1+zeroshift*2; // softabs(x,zeroshift*2)
  return((x + y)/2);
}

double softpospart_deriv(double x, double zeroshift)
{
  double y = (x*fabs(x)+4*x*zeroshift)/sqr(fabs(x)+zeroshift*2); 
  return((1 + y)/2);
}

double softstep(double x, double a, double b, double shift)
{
  double zshift = shift*(b-a);
  return((softpospart(x-a, zshift) - softpospart(x-b, zshift))/(b-a));
}

double softstep_deriv(double x, double a, double b, double shift)
{
  double zshift = shift*(b-a);
  return((softpospart_deriv(x-a, zshift) - softpospart_deriv(x-b, zshift))/(b-a));
}


//**************************************************
//**************************************************


double tansig(double x)
{
  return(2/(1+exp(-2*x))-1);
}

double dtansig_dx(double x)
{
  double y = 4/(exp(x)+exp(-x));
  return(sqr(y));
}

double sig(double x)
{
  return(1/(1+exp(-x)));
}

double dsig_dx(double x)
{
  double y = exp(-x)/sqr(1+exp(-x));
  return(y);
}

double my_tanh(double x)
{
  return((exp(2*x)-1)/(1+exp(2*x)));
}

double dtanh_dx(double x)
{
  double y = 4/sqr(exp(x)+exp(-x))   ;
  return(y);
}


//**************************************************
//**************************************************

 void tactile_error(double *tact, double *tact_goal, double *tact_error)
{
  int i, i_goal;
  double tact_size = 0;
  double centroid = 0,  mass = 0;
  double local_tact_error, tact_peakpos = tact_goal[1];
  double list_tact_goal;

  for (i = 1; i <= NT; i++) {
    centroid += i*tact[i];
    mass += tact[i];
    tact_size += softstep(tact[i], MIN_TOUCH_LEVEL, 1, MIN_TOUCH_LEVEL/10);
  }

	if (centroid < EPS) {  
		tact_error[1] = tact_error[2] = tact_error[3] = MAX_TACT_ERROR;
		return;
	}

  centroid /= mass;
  list_tact_goal = tact_goal[1]; 
  tact_error[1] = (centroid - list_tact_goal);
  for (i = 2; i <= 14; i++){ 
    list_tact_goal += 0.5;
    local_tact_error = (centroid - list_tact_goal);
    if (tact_error[1] > local_tact_error) {
      tact_error[1] = local_tact_error;
      tact_peakpos  = list_tact_goal; 
    }
  }

	i_goal = int(tact_peakpos);
	tact_error[2] = (tact[i_goal] - tact_goal[2])*10; 
    if (tact_error[2]<0) tact_error[2] = 0;
    if (tact_error[2]>10) tact_error[2] = 10;

  tact_error[3] = 0;
}


// dtact_error/ds
void tactile_error_dx(double *tact, double *tact_goal, double **Hs, double** Htact)
{
  int i, j, k, i_goal;
  double tact_size = 0;
  double centroid = 0,  mass = 0, local_tact_error, tact_peakpos = tact_goal[1];
  double list_tact_goal, min_tact_error;

  for (i = 1; i <= NTe; i++) for (j = 1; j <= ns; j++) Htact[i][j] = 0;

  for (i = 1; i <= NT; i++) {
    centroid += i*tact[i];
    mass += tact[i];
    tact_size += softstep(tact[i], MIN_TOUCH_LEVEL, 1, MIN_TOUCH_LEVEL/10);
  } 

	if (centroid < EPS) { // no contact here
		for (j = 2; j <= 1+DOF; j++) 
			for (i = 1; i < 3; i++)
				Htact[i][j] = 0;
		return;
	} 
  centroid /= mass;
  list_tact_goal = tact_goal[1];
  min_tact_error = (centroid - list_tact_goal);
  for (i = 2; i <= 14; i++){ 
    list_tact_goal += 0.5;
    local_tact_error = (centroid - list_tact_goal);
    if ( local_tact_error < min_tact_error) {
      min_tact_error = local_tact_error; // find min
      tact_peakpos  = list_tact_goal; 
    }
  }

	i_goal = int(tact_peakpos);

  for (j = 2; j <= 1+DOF; j++) {
    for (i = 1; i <= NT; i++) {
      Htact[1][j] += (i*mass - centroid)/sqr(mass)*Hs[NF+Nproprio+i][j];  
      Htact[3][j] += 0;
    }

    Htact[2][j] += Hs[NF+Nproprio+i_goal][j]*10;
  }
}    

//**************************************************
//**************************************************

void statetransition_matrix(double t, double **A){;}
void control_matrix(double t, double **B){;}
void observation_matrix(double t, double **H){;}

//**************************************************
//**************************************************


void my_predict(const double input[], int step, double output[DOF], double **CS, double **HS, int actual)
{
  float out_layer1[LSTM_SIZE];
  float out_layer2[nh2];
  float float_input[nn];
  float f;
  int i;
  int i1;

  for (i = 0; i < nn; i++) {
    float_input[i] = static_cast<float>(input[i]);
  }
  lstmForwardUsingExplicitLoops(float_input, W, R, gate_bias, CS, HS, step, actual);

  for (i = 0; i < LSTM_SIZE; i++) {
    f = 0.0F;
    for (i1 = 0; i1 < LSTM_SIZE; i1++) {
      f += Wout1[i + LSTM_SIZE * i1] * HS[step+1][i1]; 
    }
    out_layer1[i] = f + bias1[i];
  }
 
  for (i = 0; i < nh2; i++) {
    f = 0.0F;
    for (i1 = 0; i1 < LSTM_SIZE; i1++) {
      f += Wout2[i + nh2 * i1] * out_layer1[i1]; 
    }
    out_layer2[i] = f + bias2[i];
  }
  for (i = 0; i < DOF; i++) {
    f = 0.0F;
    for (i1 = 0; i1 < nh2; i1++) {
      f += Wout3[i + DOF * i1] * out_layer2[i1];
    }
    output[i] = f + bias3[i];
  }
}

void contactpush_and_velocity(double *u, double *s, double t, double dt, double *a, double *netinput, int actual)
{
  static int oldstep = 0;
  int step = int(round(t/dt))+1; 

  vecconcat(s+1,DOF,u,nc,netinput); 
  map_minmax(netinput,nn,MIN_NETINPUT,MAX_NETINPUT);
   if (actual == 1) {
    if (step > 0) my_predict(netinput+1, step, a+1, CS_actual, HS_actual, actual); 
    else zeros(a,DOF);
  }
  if (actual == 0) {
    if (step > 0) my_predict(netinput+1, step, a+1, CS_kalman, HS_kalman, actual); 
    else zeros(a,DOF);
  }

  if (actual){
  if ((step > 0) && (step-oldstep > 1) ) {char ha[256]; snprintf(ha,255,"Skipping a step in simulation! (%d->%d)",oldstep, step); mexErrMsgTxt(ha);}
  else {    oldstep = step;}}
}

void compute_net_input(double *u, double *s, double *netinput)
{
  vecconcat(s+1,DOF,u,nc,netinput); 
  map_minmax(netinput,nn,MIN_NETINPUT,MAX_NETINPUT);
}

void set_LSTM_states(double t, double dt, double *HS, double*CS, int actual)
{
  int step = int(round(t/dt))+1;
  
  if (actual) {
    veccopy(CS_actual[step]-1, CS, LSTM_SIZE);  
    veccopy(HS_actual[step]-1, HS, LSTM_SIZE);
  }
  else {
    veccopy(CS_kalman[step]-1, CS, LSTM_SIZE);
    veccopy(HS_kalman[step]-1, HS, LSTM_SIZE);
  }  
}

// Jacobian 
void dvelocity_dinput(double *net_input, double t, double dt, double **J, int actual)
{
  int i, j, k, p, b_i, i1;
  double f, f1, f2, f3, f4;
  float out_layer1[LSTM_SIZE];
  float out_layer2[nh2];
  int step = int(round(t/dt))+1; 
  double *CS_now, *HS_now, *CS_next, *HS_next; // pointers to relevant state cells
  
  // 1 LSTM layer
  float G[4*LSTM_SIZE];
  float df[nn],df1[nn],df2[nn],df3[nn];// , c[LSTM_SIZE];
  float dG[4*LSTM_SIZE][nn], dHS[LSTM_SIZE][nn];
  float dCS;
  
  if (step < 1) {
    zeros(J,DOF,nn);  
    return;
  }
  
  if (actual) {
    CS_now = CS_actual[step];
    HS_now = HS_actual[step];
    CS_next = CS_actual[step+1];
    HS_next = HS_actual[step+1];
  }
  else {
    CS_now = CS_kalman[step];
    HS_now = HS_kalman[step];
    CS_next = CS_kalman[step+1];
    HS_next = HS_kalman[step+1];
  }
  
  for (i = 0; i < LSTM_SIZE; i++) {
    f  = 0;
    f1 = 0;
    f2 = 0;
    f3 = 0;
    for (p = 0; p < nn; p++) {
      b_i = i + 4 * LSTM_SIZE * p;
      f4 = net_input[p+1]; 
      f  += W[b_i] * f4;
      f1 += W[b_i + LSTM_SIZE]*f4;
      f2 += W[b_i + 2*LSTM_SIZE]*f4;
      f3 += W[b_i + 3*LSTM_SIZE]*f4;
    }
    for (k = 0; k < LSTM_SIZE; k++) { 
      b_i = i + 4 * LSTM_SIZE * k;
      f4 = HS_now[k];
      f  += R[b_i]*f4;
      f1 += R[b_i + LSTM_SIZE]*f4;
      f2 += R[b_i + 2*LSTM_SIZE]*f4;
      f3 += R[b_i + 3*LSTM_SIZE]*f4;
    }
    G[i] = 1.0F / (std::exp(-(f + gate_bias[i]))+ 1.0F);                      // input gate
    G[i + LSTM_SIZE] = 1.0F / (std::exp(-(f1 + gate_bias[i + LSTM_SIZE])) + 1.0F);          // forget gate
    G[i + 2*LSTM_SIZE] = std::tanh(f2 + gate_bias[i + 2*LSTM_SIZE]);                              // cell gate
    G[i + 3*LSTM_SIZE] = 1.0F / (std::exp(-(f3 + gate_bias[i + 3*LSTM_SIZE])) + 1.0F);          // output gate

    CS_next[i] = G[i+2*LSTM_SIZE] * G[i]  + G[i+LSTM_SIZE] * CS_now[i];
    HS_next[i] = std::tanh(CS_next[i]) * G[i + 3*LSTM_SIZE];

    for (p = 0; p < nn; p++) { // df/dx_p
      b_i = i + 4 * LSTM_SIZE * p;
      df[p]  = W[b_i];
      df1[p] = W[b_i + LSTM_SIZE];
      df2[p] = W[b_i + 2*LSTM_SIZE];
      df3[p] = W[b_i + 3*LSTM_SIZE];


      dG[i][p] = dsig_dx(f + gate_bias[i])*df[p];                      // input gate partial derivative
      dG[i + LSTM_SIZE][p] = dsig_dx(f1 + gate_bias[i + LSTM_SIZE])*df1[p];          // forget gate partial derivative
      dG[i + 2*LSTM_SIZE][p] = dtanh_dx(f2 + gate_bias[i + 2*LSTM_SIZE])*df2[p];         // cell gate partial derivative
      dG[i + 3*LSTM_SIZE][p] = dsig_dx(f3 + gate_bias[i + 3*LSTM_SIZE])*df3[p];          // output gate partial derivative

      dCS = dG[i+2*LSTM_SIZE][p] * G[i] + G[i+2*LSTM_SIZE]*dG[i][p] +  dG[i+LSTM_SIZE][p] * float(CS_now[i]); 
      dHS[i][p] = dtanh_dx(float(CS_next[i]))*dCS*G[i + 3*LSTM_SIZE] + tanh(float(CS_next[i]))*dG[i + 3*LSTM_SIZE][p];
    }
  }

  for (p = 0; p < nn; p++) {

    for (i = 0; i < LSTM_SIZE; i++) {
      f = 0.0F;
      for (i1 = 0; i1 < LSTM_SIZE; i1++) {
        f += Wout1[i + LSTM_SIZE * i1] * dHS[i1][p]; 
      }
      out_layer1[i] = f;
    }
    for (i = 0; i < nh2; i++) {
      f = 0.0F;
      for (i1 = 0; i1 < LSTM_SIZE; i1++) {
        f += Wout2[i + nh2 * i1] * out_layer1[i1]; 
      }
      out_layer2[i] = f; 
    }
    for (i = 0; i < DOF; i++) {
      f = 0.0F;
      for (i1 = 0; i1 < nh2; i1++) {
        f += Wout3[i + DOF * i1] * out_layer2[i1];
      }
      J[i+1][p+1] = f;
    }
  }
}

void dmod_velocity_dx(double *u, double *s, double t, double dt, double *f1, int actual)
{
  int i, j;
  double clamp, clamp1, clamp2;

 
  compute_net_input(u, s, net_input);
  dvelocity_dinput(net_input, t, dt, Jnet, actual); 

  int threshold = 1; 
  double norm_gradient = 0;
  
  if (norm_gradient > threshold) {
    for ( j=1; j<=DOF; j++ )
      for ( i=1; i<=DOF; i++ )
	f1[j+1+DOF+i*ns] += (Jnet[j][i]*dmap_minmax[i]*threshold)/norm_gradient;
  }
  else {
    for ( j=1; j<=DOF; j++ )
      for ( i=1; i<=DOF; i++ )
	f1[j+1+DOF+i*ns] += Jnet[j][i]*dmap_minmax[i];//*(1+clamp);
  }
}

void dmod_velocity_dx2(double *u, double *s, double t, double dt, double **f1mat, int actual)
{
  int i, j;
  double clamp, clamp1, clamp2;

  compute_net_input(u, s, net_input);
  dvelocity_dinput(net_input, t, dt, Jnet, actual); 
  int threshold = 1; 
  double norm_gradient = 0;
  
  if (norm_gradient > threshold) {
    for ( j=1; j<=DOF; j++ )
      for ( i=1; i<=DOF; i++ )
	f1mat[j+1+DOF][i+1] += (Jnet[j][i]*dmap_minmax[i]*threshold)/norm_gradient;
  }
  else {
    for ( j=1; j<=DOF; j++ )
      for ( i=1; i<=DOF; i++ )
	f1mat[j+1+DOF][i+1] += Jnet[j][i]*dmap_minmax[i];
  }
}

// df_v/du, matrix version
void dmod_velocity_du2(double *u, double *s, double t, double dt, double **f1mat, int actual)
{
  int i, j;

  compute_net_input(u, s, net_input);
    for ( j=1; j<=DOF; j++ )
      for ( i=1; i<=nc; i++ )
	f1mat[j+1+DOF][i] += Jnet[j][i+DOF]*dmap_minmax[i+DOF];
  }



void discrete_forward(double *old_u, double *u, double *s, double t, double dt, double *f1, int actual, double *out_HS, double *out_CS, int perturb)
{
  double normu2=0, clamp;
  int i;
  bool is_k, is_t;
  if (t<Tg[Ngoals]+0.05) {
      for (i=1; i<=nc; i++)  {
         if (old_u[i]-u[i]>0) normu2 += sqr(old_u[i]-u[i]);
      }
  }
  f1[1] = s[1] + normu2/2*dt;

  double *y = dvector(1,(ny)); 
  double *s_err = dvector(1,np);
  double *tact = dvector(1,NT);
  double *tact_goal = dvector(1,NTe);
  double *tact_error = dvector(1,NTe);
  double *tact_cost_vector = dvector(1,NTe);
  int step;

	if (t<dt) {
		global_tact_cost = 0;
		global_sens_cost = 0;
		global_neuro_cost = 0;	
		Fguess1[1] = 495;
		Fguess1[2] = 1696.9;
		Fguess1[3] = 2634.9;
	}

  feedback(s, NULL, y);    

  for (int g = 1; g <= Ngoals; g++)
    if (fabs(t-Tg[g]+dt) < DT_EPS) { 
      ellipticzone_error(Cg[g], 'a', y, s_err,1);  
      ellipticzone_error(Cg[g], 'p', y+NF, s_err+NF, PROPRIO_GOAL_SCALING);  
      for (i=NF + Nproprio/3+1; i<=NF+2*Nproprio/3; i++) s_err[i] = y[i]*VEL_CONSTR_SCALING; 
        
      for (i=NF+2*Nproprio/3+1; i<=NF + Nproprio  ; i++) s_err[i] = y[i]*DVEL_CONSTR_SCALING; 
      is_k = sqr(mean_goal[g][1]-221)<=1;
      is_t = sqr(mean_goal[g][1]-187)<=1;
      if (is_k || is_t) { 
        if (is_k) {
				  tact_goal[1] =12; // 12 -> 19 defined in tactile_error
				  tact_goal[2] = 1;
        }
        else { // "t"
          tact_goal[1] = 26; // 26 -> 33, defined in tactile_error
				  tact_goal[2] = 1;
        }
				for (i=1; i<=NT; i++) tact[i] = y[NF + Nproprio + i];
				tactile_error(tact, tact_goal, tact_error);
				for (i=1; i<=NF; i++) s_err[i] = s_err[i]*AUDITORY_GOAL_SCALING; // special scaling for consonants
				for (i=1;  i<=NTe; i++) s_err[NF+Nproprio+i] = TACTILE_GOAL_SCALING*tact_error[i]; // special scaling for consonants
				f1[1] += PRECISION_TRADEOFF_SCALING*sqr(norm2(s_err,NF+Nproprio+NTe))/2;
				for (i=1; i<=NTe; i++) tact_cost_vector[i] = s_err[NF+Nproprio+i];
				global_tact_cost+= sqr(norm2(tact_cost_vector,NTe));
      } 
		  else {
			  f1[1] += PRECISION_TRADEOFF_SCALING*sqr(norm2(s_err,NF+Nproprio))/2;
		  }
    }
	
  free_dvector(y,1,ny);
  free_dvector(s_err,1,NF+Nproprio+NTe);
  free_dvector(tact,1,NT);
  free_dvector(tact_goal,1,NTe);
  free_dvector(tact_error,1,NTe);
  free_dvector(tact_cost_vector,1,NTe);
 
  for ( i=2; i<=1+DOF; i++ ) f1[i] = softstep(s[i] + s[i+DOF]*dt, 0, 1, SOFT_POS_BOUNDARY); 
  
  contactpush_and_velocity(u, s, t, dt, a, net_input, actual);
  
  if (out_HS && out_CS) {
    step = int(round(t/dt))+1;
    if (step<1) step = 1;
    if (actual == 1 ) {
      veccopy(out_HS, HS_actual[step+1]-1, LSTM_SIZE); 
      veccopy(out_CS, CS_actual[step+1]-1, LSTM_SIZE);
    }
    else {
      veccopy(out_HS, HS_kalman[step+1]-1, LSTM_SIZE);
      veccopy(out_CS, CS_kalman[step+1]-1, LSTM_SIZE);
    }
  }
  //
  
  for ( i=2+DOF; i<=1+2*DOF; i++ ) f1[i] = a[i-DOF-1]; 

  for ( i=2+2*DOF; i<=1+(D+2)*DOF; i++ ) f1[i] = s[i-DOF]; 
}

void discrete_constraints(double *u, double *s, double *cstr, double t, double dt, double *f1)
{
  double normu2=0, clamp;
  int i;
 
   for (i=1; i<=nc; i++) f1[i]= sqr(u[i] - cstr[ns-2*nc+i]);
}

void discrete_constraints_ds(double *u, double *s, double *cstr, double t, double dt, double *f1)
{  
   int i,j;
  bool is_k, is_t;

  for (i=1; i<=ns; i++)
    for (j=1; j<=ns; j++) f1[i+(j-1)*ns] = 0;  
}


// df/dx
void discrete_adjoint_dx2(double *u, double *s, double t, double dt, double **f1, int actual)
{
  int i,j;
  bool is_k, is_t;

  for (i=1; i<=ns; i++)
    for (j=1; j<=ns; j++) f1[i][j] = 0;

  double **Hss = dmatrix(1,ny,1,ns);
  double **Haudio = dmatrix(1,NF,1,ns);
  double **Htact = dmatrix(1,NTe,1,ns);
  double **Hprop = dmatrix(1,Nproprio,1,ns);
  double *y = dvector(1,ny);
  double *temp = dvector(1,NF);
  double *tact = dvector(1,NT);
  double *tact_goal = dvector(1,NTe);
  double *temp_tact = dvector(1,NTe);
  double *temp_prop = dvector(1,Nproprio);
  
  if (Fguess1 == NULL) mexErrMsgTxt("Hoho 2.\n");

  feedback(s, NULL, y);
  
  feedback_ds(s, NULL, Hss);
  for (int g = 1; g <= Ngoals; g++)
    if (fabs(t-Tg[g]+dt) < DT_EPS) { 
      
      ellipticzone_error(Cg[g], 'a', y, temp); 
      ellipticzone_error(Cg[g], 'p', y+NF, temp_prop, PROPRIO_GOAL_SCALING); 
      ellipticzone_error_ds(Cg[g], 'a', y, Hss, Haudio);
      ellipticzone_error_ds(Cg[g], 'p', y, Hss, Hprop);
      
      is_k = sqr(mean_goal[g][1]-221)<=1;
      is_t = sqr(mean_goal[g][1]-137)<=1;

      if (is_k || is_t) {
        if (is_k) {
          tact_goal[1] = 12; 
          tact_goal[2] = 1; 
        }
        else {
          tact_goal[1] = 26; 
          tact_goal[2] = 1; 
        }
        for (i=1; i<=NF; i++) temp[i] = temp[i]*AUDITORY_GOAL_SCALING;
        for (i=1; i<=NT; i++) tact[i] = y[NF + Nproprio + i];  
        tactile_error(tact,tact_goal,temp_tact);
        tactile_error_dx(tact, tact_goal, Hss, Htact);
        for (j=1; j<=DOF; j++) for (i=1; i<=NF; i++) f1[1][j+1] += PRECISION_TRADEOFF_SCALING*(AUDITORY_GOAL_SCALING*Haudio[i][j+1]*temp[i] + sqr(TACTILE_GOAL_SCALING)*Htact[i][j+1]*temp_tact[i] + PROPRIO_GOAL_SCALING*Hprop[i][j+1]*temp_prop[i]);
      }
      else {
        for (j=1; j<=DOF; j++) for (i=1; i<=NF; i++) f1[1][j+1] += PRECISION_TRADEOFF_SCALING*(Haudio[i][j+1]*temp[i]+Hprop[i][j+1]*temp_prop[i]*PROPRIO_GOAL_SCALING); 
      }
      for (j=1; j<=DOF; j++) for (i=Nproprio/3+1; i<=2*Nproprio/3; i++) f1[1][j+DOF+1] += PRECISION_TRADEOFF_SCALING*sqr(VEL_CONSTR_SCALING)*Hss[i][1+DOF+j]*y[NF+i];
      for (j=1; j<=DOF; j++) for (i=2*Nproprio/3+1; i<=Nproprio;   i++) f1[1][j+DOF+1] += PRECISION_TRADEOFF_SCALING*sqr(DVEL_CONSTR_SCALING)*Hss[i][1+DOF+j]*y[NF+i];
      for (j=1; j<=DOF; j++) for (i=2*Nproprio/3+1; i<=Nproprio;   i++) f1[1][j+(D+1)*DOF+1] += PRECISION_TRADEOFF_SCALING*sqr(DVEL_CONSTR_SCALING)*Hss[i][1+(D+1)*DOF+j]*y[NF+i];  
    }

  free_dvector(y,1,ny);
  free_dvector(tact,1,NT);
  free_dvector(tact_goal,1,NTe);
  free_dvector(temp,1,NF);
  free_dvector(temp_tact,1,NTe);
  free_dvector(temp_prop,1,Nproprio);
  free_dmatrix(Hss,1,ny,1,ns);
  free_dmatrix(Haudio,1,NF,1,ns);
  free_dmatrix(Htact,1,NTe,1,ns);
  free_dmatrix(Hprop,1,Nproprio,1,ns);

  f1[1][1] = 1.0;

  for ( i=2; i<=1+DOF; i++ ) {
        f1[i][i] = softstep_deriv(s[i] + s[i+DOF]*dt, 0, 1, SOFT_POS_BOUNDARY);
        f1[i][i+DOF] = softstep_deriv(s[i] + s[i+DOF]*dt, 0, 1, SOFT_POS_BOUNDARY)*dt;
  }
  
  dmod_velocity_dx2(u, s, t, dt, f1, actual);
  for (i=2+2*DOF; i<=1+(D+2)*DOF; i++) f1[i][i-DOF] = 1;

}



// df/du
void discrete_adjoint_du2(double *old_u ,double *u, double *s, double t, double dt, double **f2, int actual)
{
  int i,j;

  for (i=1; i<=ns; i++) for (j=1; j<=nc; j++) f2[i][j] = 0;
  for (i=1; i<=nc; i++) f2[1][i] = (-gfun(old_u[i]-u[i],sa) - (old_u[i]-u[i])*gfunderiv(old_u[i]-u[i],sa))/2*dt;
  
  dmod_velocity_du2(u, s, t, dt, f2, actual);

}

void discrete_adjoint_dx(double *u, double *s, double t, double dt, double *f1, int actual)
{
 
}



void discrete_adjoint_du(double *old_u, double *u, double *s, double t, double dt, double *f2)
{
  
}





void discrete_cost(double *s, double *cost, double *sens_cost, double *tact_cost, double *neuro_cost)
{
  *cost = s[1];
  *sens_cost = global_sens_cost;
  *tact_cost = global_tact_cost;
  *neuro_cost = global_neuro_cost;
}

void costNumericalJacobian(double *s, double *jacobian)
{
  int i, j;
  double *mys = dvector(1,ns);
  double cost1, cost2, *cost3, *cost4, *cost5;

  veccopy(mys,s,ns);

  for (j = 1; j <= ns; j++) {
    mys[j] -= DELTA;
    discrete_cost(mys, &cost1, cost3, cost4, cost5);
    mys[j] += 2*DELTA;
    discrete_cost(mys, &cost2, cost3, cost4, cost5);
    mys[j] -= DELTA; 
    jacobian[j] = (cost2-cost1)/2/DELTA;
  }

  free_dvector(mys,1,ns);
}



void discrete_cost_dx(double *s, double *f2)
{
  int i;
  f2[1] = 1;
  for (i=2; i<=ns; i++) f2[i] = 0;
}



#ifndef NO_MATLAB
bool mexinit(int *pns, int *pnp, int *pnc)
{
  int i;
  mxArray *psvdMvar, *pmuvar, *pmivar, *pMivar ;
  double *pM;

  psvdMvar= mexGetVariable("global", "svdmatrix");
  if (!psvdMvar) {mexErrMsgTxt("svdmatrix : matrice pas trouvee\n"); usleep(5e6); return(false);}
  Nsvd = mxGetN(psvdMvar);
  pmuvar= mexGetVariable("global", "mu");
  if (!pmuvar) {mexErrMsgTxt("mu : vecteur pas trouve\n"); usleep(5e6); return(false);}
  pmivar= mexGetVariable("global", "minInput");
  if (!pmivar) {mexErrMsgTxt("minInput : vecteur pas trouve\n"); usleep(5e6); return(false);}
  pMivar= mexGetVariable("global", "maxInput");
  if (!pMivar) {mexErrMsgTxt("maxInput : vecteur pas trouve\n"); usleep(5e6); return(false);}


  if (!autoenc_initialized && !autoencoder_init()) return(false);

  init_goal_geometry(); 
  nc = 7;
  *pns = ns = 1 + (D+2)*DOF;
  *pnp = np = NF + Nproprio + NT;
  nyp = 2*(DOF+nc); 
  *pnc = nc;

  nn = DOF + nc;
  nh2 = 12; 

  net_input = dvector(1,nn);

  MIN_NETINPUT = dvector(1,nn);
  MAX_NETINPUT = dvector(1,nn);
  dmap_minmax = dvector(1,nn);

  Jnet = dmatrix(1,DOF,1,nn);

  CS_actual = dmatrix(1,MAX_STEPS, 0,LSTM_SIZE-1); 
  HS_actual = dmatrix(1,MAX_STEPS, 0,LSTM_SIZE-1);
  CS_kalman = dmatrix(1,MAX_STEPS, 0,LSTM_SIZE-1); 
  HS_kalman = dmatrix(1,MAX_STEPS, 0,LSTM_SIZE-1);

  saved_CS_kalman = dvector(1,LSTM_SIZE); 
  saved_HS_kalman = dvector(1,LSTM_SIZE);	
  saved_CS_actual = dvector(1,LSTM_SIZE); 
  saved_HS_actual = dvector(1,LSTM_SIZE);	

  veccopy(MIN_NETINPUT, mxGetPr(pmivar)-1, nn);
  veccopy(MAX_NETINPUT, mxGetPr(pMivar)-1, nn);

  for (i = 1; i<=nn; i++) dmap_minmax[i] = 2/(MAX_NETINPUT[i]-MIN_NETINPUT[i]);

  svdmatrix = dmatrix(1,2*N_NODES,1,DOF);
  averagetongue = dvector(1,2*N_NODES);

  v = dvector(1,DOF);
  a = dvector(1,DOF);
  sdeltas = dvector(1,1+2*DOF);
  scopy = dvector(1,ns);
    
  pM = mxGetPr(psvdMvar);
  matconvert(pM, svdmatrix, 2*N_NODES, DOF);

  veccopy(averagetongue, mxGetPr(pmuvar)-1, 2*N_NODES);


  invtau_e = 1.0/tau_e;
  invtau_a = 1.0/tau_a;

  Fguess1 = dvector(1,NF);

  init_feedback(); 
  return(true);
}



void init_state_constraints(double *thf)
{
  th_goal = dvector(1,NF+DOF+NTe+2);   

  for (int i=1; i<=NF+DOF+NTe+2; i++) th_goal[i] = thf[i];
}
#endif


void  init_goals(double *goalseq, char *goalcodeseq, int Ndims, int Ng)
{
  int i, j, k1, k2, g, m;
  Ngoals = Ng;
  Cg = (char*) malloc((Ngoals+1)*sizeof(char))-1; 
  Tg = dvector(1,Ngoals); 
  Ngoal_modalities = Ndims - (ny-NT-2*Nproprio/3) - 1; 
  Wg = dmatrix(1,Ngoals,1,3);

  mean_goal = dmatrix(1,Ngoals,1,ny-NT-2*Nproprio/3);  
  for (g = 1, k1 = 0, k2 = 0; g <= Ng; g++) {
    Cg[g] = goalcodeseq[g-1];
    for (i=1; i<=ny-NT-2*Nproprio/3; i++) mean_goal[g][i] = goalseq[k1++];
    Tg[g] = goalseq[k1++]; 
    for (m = 1; m <= 3; m++) Wg[g][m] = goalseq[k1++];
  }
  
}

void free_plant_variables()
{
#ifndef NO_MATLAB
  
  if (svdmatrix) free_dmatrix(svdmatrix,1,2*N_NODES,1,DOF);

  if (averagetongue) free_dvector(averagetongue,1,2*N_NODES);

  if (MIN_NETINPUT) free_dvector(MIN_NETINPUT ,1,nn);
  if (MAX_NETINPUT) free_dvector(MAX_NETINPUT ,1,nn);
  if (net_input)  free_dvector(net_input,1,nn);
	
  if (saved_HS_kalman) free_dvector(saved_HS_kalman,1,LSTM_SIZE);
  if (saved_CS_kalman) free_dvector(saved_CS_kalman,1,LSTM_SIZE);
  if (saved_HS_actual) free_dvector(saved_HS_actual,1,LSTM_SIZE);
  if (saved_CS_actual) free_dvector(saved_CS_actual,1,LSTM_SIZE);

  if (dmap_minmax) free_dvector(dmap_minmax,1,nn);
  if (Jnet) free_dmatrix(Jnet,1,DOF,1,nn);

  if (CS_actual) free_dmatrix(CS_actual, 1,MAX_STEPS, 0,LSTM_SIZE-1);
  if (HS_actual) free_dmatrix(HS_actual, 1,MAX_STEPS, 0,LSTM_SIZE-1);
  if (CS_kalman) free_dmatrix(CS_kalman, 1,MAX_STEPS, 0,LSTM_SIZE-1);
  if (HS_kalman) free_dmatrix(HS_kalman, 1,MAX_STEPS, 0,LSTM_SIZE-1);

#endif

  if (v) free_dvector(v,1,DOF);
  if (a) free_dvector(a,1,DOF);
  if (sdeltas) free_dvector(sdeltas,1,1+2*DOF);
  if (scopy) free_dvector(scopy,1,ns);
  if (Fguess1) free_dvector(Fguess1,1,NF);

}

void free_state_constraints_variables()
{
  free_dvector(th_goal,1,NF+DOF+NTe+2);
}


void free_goals_variables()
{
  if (mean_goal) free_dmatrix(mean_goal,1,Ngoals,1,ny);
  if (cov_goal) free_dmatrix(cov_goal,1,Ngoals*NF,1,NF);

  if (Cg) free(Cg+1);
  if (Tg) free_dvector(Tg,1,Ngoals);
  if (Wg) free_dmatrix(Wg,1,Ngoals,1,Ngoal_modalities); 

}
