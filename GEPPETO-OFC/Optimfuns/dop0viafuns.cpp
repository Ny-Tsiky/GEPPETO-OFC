#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "mylinpack.h"
#include "nr.h"
#include "nrutil.h"
#include "kalman8.h" 
#include "plantfuns_via.h"
#include "unistd.h" 

#undef DB
#define DB mexPrintf("Got to line %d in mex file %s.\n",__LINE__,__FILE__);

bool resetLSTM;
bool resetLSTMderivatives;

void mexmatdisplay(double **M, int m, int n, const char *matname);
void mexvecdisplay(double *M, int m, const char *matname);


  
void myExitFcn() {
  mexPrintf("MEX-file is being unloaded, freeing memory\n\n");
  free_plant_variables();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *old_u, *u, *s, *cstr, dt, t, *f1, *f2, *sens_cost, *tact_cost, *neuro_cost, *ls, *HS, *CS, **m, **straj, *tspan;
  double *in_HS, *in_CS, *out_HS, *out_CS;
  int i, j, flag, c, perturb;
  static int ns, np, nc;
  static int first_time=1;
  mxArray *p_goalseq, *p_goalcodeseq;
	
  if (nrhs!=10) { mexErrMsgTxt("dop0viafuns requires 10 args\n");}

  mexAtExit(myExitFcn);
	
  if (first_time) {
    if (!mexinit(&ns,&np,&nc)) return;
    first_time=0;  
  }
	

  c=0;
  old_u=mxGetPr(prhs[c++])-1;
  u=mxGetPr(prhs[c++])-1;
  s=mxGetPr(prhs[c++])-1;  
  cstr=mxGetPr(prhs[c++])-1;  
 
  in_HS=mxGetPr(prhs[c++])-1;
  in_CS=mxGetPr(prhs[c++])-1; 

  dt=*mxGetPr(prhs[c++]);   
  tspan = mxGetPr(prhs[c++]);
  t=*tspan;   

  flag=(int)*mxGetPr(prhs[c++]);
  perturb = (int)*mxGetPr(prhs[c++]);

  switch (flag) {
    case 1:
    {
      double normu2=0;

      if (nlhs>0) 
        mxDestroyArray(plhs[0]);
			plhs[0]=mxCreateDoubleMatrix(ns,1,mxREAL);
      f1=mxGetPr(plhs[0])-1;

		if (nlhs>1) 
        mxDestroyArray(plhs[1]);
        plhs[1]=mxCreateDoubleMatrix(1,20,mxREAL);
        out_HS=mxGetPr(plhs[1])-1;

		if (nlhs>2) 
        mxDestroyArray(plhs[2]);
        plhs[2]=mxCreateDoubleMatrix(1,20,mxREAL);
        out_CS=mxGetPr(plhs[2])-1;
      
      discrete_forward(old_u, u, s, t, dt, f1, 0, out_HS, out_CS, perturb); 
    }
    break;
    
    case 3: 
    {
      if (nlhs>0) 
        mxDestroyArray(plhs[0]);

      plhs[0]=mxCreateDoubleMatrix(ns,ns,mxREAL);
      f1=mxGetPr(plhs[0])-1;
      discrete_adjoint_dx(u, s, t, dt, f1, 0);  
      if (nlhs>1) {
        mxDestroyArray(plhs[1]);
        plhs[1]=mxCreateDoubleMatrix(ns,nc,mxREAL);
        f2=mxGetPr(plhs[1])-1;

      discrete_adjoint_du(old_u, u, s, t, dt, f2);
      }
    }
    break;

    case 2:
    if (nlhs>0)

    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    f1 = mxGetPr(plhs[0]);

    plhs[2]=mxCreateDoubleMatrix(1,1,mxREAL);
    sens_cost = mxGetPr(plhs[2]);
	
    plhs[3]=mxCreateDoubleMatrix(1,1,mxREAL);
    tact_cost = mxGetPr(plhs[3]);
	
    plhs[4]=mxCreateDoubleMatrix(1,1,mxREAL);
    neuro_cost = mxGetPr(plhs[4]);
	
    discrete_cost(s, f1,sens_cost, tact_cost, neuro_cost);
      
    if (nlhs>1) {
      mxDestroyArray(plhs[1]);
      plhs[1]=mxCreateDoubleMatrix(1,ns,mxREAL);
      f2=mxGetPr(plhs[1])-1;
	
      discrete_cost_dx(s, f2);
    }
    break;
      
    case 4:
    p_goalseq = mexGetVariable("global", "goalseq");
    if (!p_goalseq) {mexErrMsgTxt("goalseq pas trouve\n"); return;}
    p_goalcodeseq = mexGetVariable("global", "goalcodeseq");
    if (!p_goalcodeseq) {mexErrMsgTxt("goalcodeseq pas trouve\n"); return;}
    
    {
      int Ngoals = mxGetN(p_goalseq);
      char *goalcodeseq = (char*)malloc((Ngoals+1)*sizeof(char));
      mxGetString(p_goalcodeseq, goalcodeseq, Ngoals+1); 
      init_goals(mxGetPr(p_goalseq), goalcodeseq, mxGetM(p_goalseq), Ngoals); 
      free(goalcodeseq);
    }
    break;
    
    case 5:  
    free_goals_variables();
    break;   
	
    case 6:
    set_LSTM_states(t,dt,in_HS,in_CS,0); 
    break;
      
    case 7:
    if (nlhs>0) 
      mxDestroyArray(plhs[0]);
    plhs[0]=mxCreateDoubleMatrix(ns,1,mxREAL);
    f1=mxGetPr(plhs[0])-1;
    discrete_constraints(old_u, s, cstr, t, dt, f1);
	  
    if (nlhs>1) {
      mxDestroyArray(plhs[1]);
      plhs[1]=mxCreateDoubleMatrix(ns,ns,mxREAL);
      f2=mxGetPr(plhs[1])-1;
	
      discrete_constraints_ds(u, s, cstr, t, dt, f2);
    }
    break;
    
    
    case 8: 

    int N = max(mxGetM(prhs[7]), mxGetN(prhs[7])); 

    plhs[0] = mxCreateDoubleMatrix(ns,N+1,mxREAL);
    ls = mxGetPr(plhs[0]);
    
    if (nlhs>1) {
      plhs[1]=mxCreateDoubleMatrix(1,20,mxREAL);
      out_HS=mxGetPr(plhs[1]);
    }
    if (nlhs>2) {
      plhs[2]=mxCreateDoubleMatrix(1,20,mxREAL);
      out_CS=mxGetPr(plhs[2]);
    }

    const int LSTM_SIZE = 20;
    HS = dvector(1,LSTM_SIZE);
    CS = dvector(1,LSTM_SIZE);
    m  = dmatrix(1,N,1,nc);
    straj  = dmatrix(1,N+1,1,ns);
    
    int k = 1; 
    for ( i=1; i<=N; i++ )
      for ( j=1; j<=nc; j++ )
        m[i][j] = u[k++];  

    tspan--; 

    veccopy(straj[1],s,ns);

    discrete_forward(m[1], m[1], straj[1], tspan[1], dt, straj[2], 0, HS, CS, perturb); 
    for (i = 2; i <= N; i++)
      discrete_forward(m[i-1], m[i], straj[i], tspan[i], dt, straj[i+1], 0, HS, CS, perturb); 

    k = 0;
    for ( i=1; i<=N+1; i++ )
      for ( j=1; j<=ns; j++ )
        ls[k++] = straj[i][j]; 
    if (nlhs>1) {
      k = 0;
      for ( j=1; j<=LSTM_SIZE; j++ )
        out_HS[k++] = HS[j]; // 
    }
    if (nlhs>2) {
      k = 0;
      for ( j=1; j<=LSTM_SIZE; j++ )
        out_CS[k++] = CS[j]; // 
    }
    
      free_dvector(HS,1,LSTM_SIZE);
      free_dvector(CS,1,LSTM_SIZE);
      free_dmatrix(m ,1,N,1,nc);
      free_dmatrix(straj ,1,N+1,1,ns);
      
    break;
      
  }
}

