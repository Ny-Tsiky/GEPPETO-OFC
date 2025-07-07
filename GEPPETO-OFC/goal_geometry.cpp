#include <string.h>
#include <math.h>
#include "mex.h"
#include "mylinpack.h"
#include "nr.h"
#include "nrutil.h"
#include "autoencoder.h" // for hz2bark
#include "goal_geometry.h"

extern int NF, Nproprio, ns;

void mexvecdisplay(double *v, int N, char *s);
void mexmatdisplay(double **M, int m, int n, const char *matname);


double *audio_mean[N_SpeechGoals], *proprio_mean[N_SpeechGoals];
double **audio_invcovmat[N_SpeechGoals], **proprio_invcovmat[N_SpeechGoals];

void ellipticzone_error(char target_code, char target_modality, double *y, double *err, double scaling)
{
	int N, i, k;
	double *s_err = dvector(1,NF);
	double *mean_goal, **invcov_goal;
  double *y_centered, *y_circle, *Py;
	double threshold = 2.8; // from Hotelling T2 asympotic distribution for 3D data
	double norm_y;
  
  k = int(strchr(SpeechGoalName, target_code) - SpeechGoalName); // index of target_code  
  if (k < 0 || k >= N_SpeechGoals) {zeros(err,N); return;} // safety
  // mexPrintf("goal code = %d\n",k); // DBG
  
  switch (target_modality) {
    case 'a':
    N = NF;
    mean_goal = audio_mean[k];
    invcov_goal = audio_invcovmat[k];  
    break;
    case 'p':
    N = Nproprio/3;
    mean_goal = proprio_mean[k];
    invcov_goal = proprio_invcovmat[k];
    break;
  }

  y_centered = dvector(1,N);
  y_circle = dvector(1,N);
  Py = dvector(1,N);
  
	vecsub(y,mean_goal,y_centered, N);
	matvecmult(N,invcov_goal,N, y_centered, y_circle);
	norm_y = norm2(y_circle,N);
  
	if (norm_y < threshold)
			zeros(err,N);
	else { // we are not inside the confidence ellipse
			for (i=1; i<=N; i++) Py[i] = y_centered[i]/norm_y*threshold; // nearest point on ellipse border
			for (i=1; i<=N; i++) err[i] = scaling*(y[i] - (Py[i]+mean_goal[i]));
	}

//mexvecdisplay(err,NF,"ellipse error");

	free_dvector(y_centered,1,N);
	free_dvector(y_circle,1,N);
	free_dvector(Py,1,N);
}



// H = dh/ds
void ellipticzone_error_ds(char target_code, char target_modality, double *y, double **H, double **He, double scaling)
{
	int N, i, k;
	double *s_err = dvector(1,NF);
	double *mean_goal, **invcov_goal;
  double *y_centered, *y_circle;
	double threshold = 2.8; // from Hotelling T2 asympotic distribution for 3D data
	double norm_y;
  
  k = int(strchr(SpeechGoalName, target_code) - SpeechGoalName); // index of target_code  
  if (k < 0 || k >= N_SpeechGoals) {zeros(He,N,N); return;} // safety
  // mexPrintf("goal code = %d\n",k); // DBG
  
  switch (target_modality) {
    case 'a':
    N = NF;
    mean_goal = audio_mean[k];
    invcov_goal = audio_invcovmat[k]; 
    break;
    case 'p':
    N = Nproprio/3;
    mean_goal = proprio_mean[k];
    invcov_goal = proprio_invcovmat[k];
    break;
  }

  y_centered = dvector(1,N);
  y_circle = dvector(1,N);
  
	vecsub(y,mean_goal,y_centered, N);
	matvecmult(N,invcov_goal,N, y_centered, y_circle);
	norm_y = norm2(y_circle,N);
  
	if (norm_y < threshold) {
		zeros(He,N,ns);
	}
	else {
  	double **IN = dmatrix(1,N,1,N);
  	double **temp = dmatrix(1,N,1,N);
  	double **temp2 = dmatrix(1,N,1,N);
  	double **temp3 = dmatrix(1,N,1,ns);
  	double **tinvcov_goal = dmatrix(1,N,1,N);
    
    // derr/dy = scaling.[I - thresh/normy(I - ycent*ycent'*C'*C/normy^2) ]
		eye(IN,N);
		todprod(y_centered,y_centered,N,temp); // outer prod
		transpose(invcov_goal,N,N,tinvcov_goal); // should not be necessary but since Tsiky's cov matrices are not always symmetrical (...) XXXX
		matmult(N,temp,N,tinvcov_goal,N,temp2);
		matmult(N,temp2,N,invcov_goal,N,temp); // (YY')C'C
		scalarmult(temp,N,N,1./sqr(norm_y));
		matsub(IN,temp,temp2,N,N);
		scalarmult(temp2,N,N,threshold/norm_y);
		matsub(IN,temp2,temp,N,N);
		scalarmult(temp,N,N,scaling);
    switch (target_modality) {
      case 'a':
  		matextract(H,1,NF,1,ns,temp3); // dh/dy_audio
      break;
      case 'p':
  		matextract(H,NF+1,NF+Nproprio/3,1,ns,temp3);  // dh/dy_proprio
      break;
    }
		matmult(N,temp,N,temp3,ns,He);
    
    free_dmatrix(IN,1,N,1,N);
  	free_dmatrix(temp, 1,N,1,N);
  	free_dmatrix(temp2,1,N,1,N);
  	free_dmatrix(temp3,1,N,1,ns);
  	free_dmatrix(tinvcov_goal,1,N,1,N);
	}
//mexvecdisplay(normalized_s_err,NF,"ellipse error");

	free_dvector(y_centered,1,N);
	free_dvector(y_circle,1,N);
}




void init_goal_geometry()
{
  mxArray *p_audioInvCov, *p_proprioInvCov, *p_audioMean, *p_proprioMean;
  double *audioMean, *proprioMean, *audioInvCov, *proprioInvCov;
  int i, j, k, l, m, totalnum_goals;

  p_audioMean = mexGetVariable("global", "audioGoalMean");
  if (!p_audioMean) {mexErrMsgTxt("audioGoalMean pas trouve\n"); return;}
  audioMean = mxGetPr(p_audioMean);
  
  p_proprioMean = mexGetVariable("global", "proprioGoalMean");
  if (!p_proprioMean) {mexErrMsgTxt("proprioGoalMean pas trouve\n"); return;}
  proprioMean = mxGetPr(p_proprioMean);

  p_audioInvCov = mexGetVariable("global", "audioGoalInvCov");
  if (!p_audioInvCov) {mexErrMsgTxt("audioGoalInvCov pas trouve\n"); return;}
  audioInvCov = mxGetPr(p_audioInvCov);

  p_proprioInvCov = mexGetVariable("global", "proprioGoalInvCov");
  if (!p_proprioInvCov) {mexErrMsgTxt("proprioGoalInvCov pas trouve\n"); return;}
  proprioInvCov = mxGetPr(p_proprioInvCov);


  totalnum_goals = mxGetM(p_audioMean);
  
  for (k = 0, l = 0, m = 0; k < totalnum_goals; k++) {
    // set ellipse center for each goal
    audio_mean[k] = dvector(1,NF);
    for (i = 1; i <= NF; i++) audio_mean[k][i] = audioMean[k+(i-1)*N_SpeechGoals]; //Hz2bark(audioMean[k+(i-1)*N_SpeechGoals]); // more consistent if done before in Matlab!
    
    proprio_mean[k] = dvector(1,Nproprio/3);
    for (i = 1; i <= NF; i++) proprio_mean[k][i] = proprioMean[k+(i-1)*N_SpeechGoals];
    
    // set inverse cov mat for each goal
    audio_invcovmat[k] = dmatrix(1,NF,1,NF);
    for (j = 1; j <= NF; j++)
      for (i = 1; i <= NF; i++)
        audio_invcovmat[k][i][j] = audioInvCov[l++];

    proprio_invcovmat[k] = dmatrix(1,Nproprio/3,1,Nproprio/3);
    for (j = 1; j <= Nproprio/3; j++)
      for (i = 1; i <= Nproprio/3; i++)
        proprio_invcovmat[k][i][j] = proprioInvCov[m++];
  }
  
  // mexvecdisplay(audio_mean[0], NF, "audio mean target for goal /i/"); // DBG XXXX
  // mexmatdisplay(audio_invcovmat[0], NF, NF, "audio cov mat for goal /i/"); // DBG XXXX
}