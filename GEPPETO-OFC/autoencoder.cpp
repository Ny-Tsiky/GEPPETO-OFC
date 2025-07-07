#include  "autoencoder.h"
#include <math.h>
#include "mylinpack.h"
#include "nr.h"
#include "nrutil.h"
#include <stdlib.h>
#include "mex.h"
#include <string.h>

#ifndef NO_MATLAB
#include "mex.h"  
#undef DB
#define DB mexPrintf("Got to line %d in mex file %s.\n",__LINE__,__FILE__);

#endif

static double *bE;
static double **WE;
static double *m;
static double *M;
static double *bD;
static double **WD;
double **JD;
const int N_NODES = 16;
extern double **svdmatrix;

bool autoenc_initialized = false;


void mexvecdisplay(double *v, int N, char *s)
{
	mexPrintf("Vector %s:\n",s);
	for (int i=1; i<=N; i++) mexPrintf("%g\t",v[i]);
	mexPrintf("\n");
}

void mexmatdisplay(double **A, int m, int n, char *s)
{
	mexPrintf("Matrix %s (%dx%d):\n",s,m,n);
	for (int i=1; i<=m; i++) {
		for (int j=1; j<=n; j++) mexPrintf("%g\t",A[i][j]);
		mexPrintf("\n");
	}
	mexPrintf("\n");
}

void satlin(int n,double *p){
  for(int i = 1; i<=n; i++){
    if(p[i]<=0)
      p[i] = 0;
    if(p[i]>=1)
      p[i] = 1;
  }
}


void map_minmax(double *v, const int n, const double *m, const double *M)
{
  int i;
  for (i=1; i<=n; i++){
		if (m[i] < M[i]) v[i] = 2*(v[i] - m[i])/(M[i]-m[i]) - 1;
	}
}

void unmap_minmax(double *v, const int n, const double *m, const double *M)
{
  int i;

  for (i=1; i<=n; i++) v[i] = (v[i] + 1)*(M[i]-m[i])/2 + m[i];
}

void
matvecmult1(int m, double **A, int n,  double *u, double *v)
{
  int i,k;
  double S;

  for (i=1; i<=m; i++) {
    S=0;
    for (k=1; k<=n; k++)
		{
			S+=A[i][k]*u[k];
		}
    v[i]=S;
  }
}

void encode(double *p, double *result){
  double *temp, *tmp;
  int n1 = 2*N_NODES;
  int n2 = DOF;

 tmp = dvector(1,n1);
 temp = dvector(1,n2);
	if (!autoenc_initialized)  if (!autoencoder_init()) {mexErrMsgTxt("Erreur a l'initialisation autoencodeur\n"); }

	veccopy(tmp,p,n1);
  map_minmax(tmp,n1,m,M);
  matvecmult(n2,WE,n1,tmp,temp);
  vecadd(temp,bE,result,n2);
  satlin(n2,result);


  free_dvector(tmp,1,n1);
  free_dvector(temp,1,n2);
}

void decode(double *p, double *result){
  double *temp;
  int n2 = DOF;
  int n1 = 2*N_NODES;

  temp = dvector(1,n1);
	if (!autoenc_initialized) if (!autoencoder_init()) {mexErrMsgTxt("Erreur a l'initialisation autoencodeur (2)\n"); };
  matvecmult1(n1,WD,n2,p,temp);
  vecadd(temp,bD,result,n1);
  unmap_minmax(result,n1,m,M);
  free_dvector(temp,1,n1);
  
  if (result[N_NODES+1] < 0.1) result[N_NODES+1] = 0.1;
  if (result[N_NODES+2] < 0.2) result[N_NODES+2] = 0.2;
  if (result[N_NODES+3] < 0.3) result[N_NODES+3] = 0.3;
  // **************************************************
}

double Hz2bark(double fHz)
{
	double fbark = ((26.81*fHz)/(1960+fHz))-0.53;
	if (fbark < 2)    fbark = fbark + (0.15)*(2-fbark);
	if (fbark > 20.1) fbark = fbark + (0.22)*(fbark-20.1);
  return(fbark);
}

void hz2bark(double *y, double *bark_y, int NF){
	for (int i=1; i<=NF; i++) {
    bark_y[i] = Hz2bark(y[i]);
  }
}

void bark2hz(double *bark_y, double *y, int NF){
	 for (int i=1; i<=NF; i++) {
		if (bark_y[i] < 2) bark_y[i] = (bark_y[i]-0.3)/0.85;
		if (bark_y[i] > 20.1) bark_y[i] = (bark_y[i]+4.422)/1.22;
		y[i] = 1960*((bark_y[i]+0.53)/(26.28-bark_y[i]));
	 }
}

//**************************************************
//**************************************************

void matconvert(double *pA, double **A, int m, int n)
{
  int i, j, k = 0;

  if (A == NULL) A = dmatrix(1,m,1,n);

  for ( j=1; j<=n; j++ )
    for ( i=1; i<=m; i++ )
      A[i][j] = pA[k++];  //checked
}
//**************************************************
//**************************************************

void DistanceEllipse(int g, int NF, double **mean_goal, double *y, double *normalized_s_err)
{
	double test = 0;
	double *s_err = dvector(1,NF);
	double *bark_mean_goal = dvector(1,NF);
    double *new_y = dvector(1,NF);
	double *y_circle = dvector(1,NF);
	double *Py = dvector(1,NF);
	double threshold = 2.8; // from Hotelling T2 asympotic distribution for 3D data
	double **covX = dmatrix(1,NF,1,NF);
	double norm_y;

	if (sqr(mean_goal[g][1]-297)<=1) {
		covX[1][1] = 45.9692; 
		covX[1][2] = 38.4313; 
		covX[1][3] = 27.6545; 
		covX[2][1] = 38.4313; 
		covX[2][2] = 47.8899;  
		covX[2][3] = 24.0003; 
		covX[3][1] = 27.6545; 
		covX[3][2] = 24.0003; 
		covX[3][3] = 29.0260; 
	}
	if (sqr(mean_goal[g][1]-383)<=1) {
		covX[1][1] = 11.5083;
		covX[1][2] = 1.5335;  
		covX[1][3] = 2.4926; 
		covX[2][1] = 3.0670; 
		covX[2][2] = 6.2217; 
		covX[2][3] = 0.4883; 
		covX[3][1] = 2.4926;
		covX[3][2] = 0.2441; 
		covX[3][3] = 22.1739;                           
	}
	if (sqr(mean_goal[g][1]-453)<=1) {
		covX[1][1] = 30.5021;  
		covX[1][2] = 10.6101; 
		covX[1][3] = -5.2994;
		covX[2][1] = 10.6101; 
		covX[2][2] = 13.2495; 
		covX[2][3] = -2.1955; 
		covX[3][1] = -5.2994; 
		covX[3][2] = -2.1955; 
		covX[3][3] = 15.9632;                   
    }
	if (sqr(mean_goal[g][1]-554)<=1) {
		covX[1][1] = 37.1704; 
		covX[1][2] = 3.1192; 
		covX[1][3] = 1.6680; 
		covX[2][1] = 3.1192; 
		covX[2][2] = 7.3697;
		covX[2][3] = -1.8138;
		covX[3][1] = 1.6680; 
		covX[3][2] = -1.8138; 
		covX[3][3] = 10.6251; 
	}
	if (sqr(mean_goal[g][1]-471)<=1) {
		covX[1][1] = 5.4333; 
		covX[1][2] = -0.6941; 
		covX[1][3] = -0.8982; 
		covX[2][1] = -0.6941;
		covX[2][2] = 8.4192; 
		covX[2][3] = -4.9052;  
		covX[3][1] = -0.8982; 
		covX[3][2] = -4.9052; 
		covX[3][3] = 11.4231; 
	}
	if (sqr(mean_goal[g][1]-449)<=1) {
		covX[1][1] = 5.5042; 
		covX[1][2] = -1.2432; 
		covX[1][3] = 0.2507; 
		covX[2][1] = -1.2432; 
		covX[2][2] = 7.5762;
		covX[2][3] = -2.0822;
		covX[3][1] = 0.2507;
		covX[3][2] = -2.0822; 
		covX[3][3] = 12.0402;  
	}
	if (sqr(mean_goal[g][1]-221)<=1) { //
		covX[1][1] = 4.8385;
		covX[1][2] = -0.2272; 
		covX[1][3] = 0.7820; 
		covX[2][1] = -0.2272;  
		covX[2][2] = 2.2762; 
		covX[2][3] = -0.6481;
		covX[3][1] = 0.7820;
		covX[3][2] = -0.6481; 
		covX[3][3] = 5.0038;
	}
	if (sqr(mean_goal[g][1]-187)<=1) {  
		covX[1][1] = 2.4693; 
		covX[1][2] = 0.1683; 
		covX[1][3] = -1.8044; 
		covX[2][1] = 0.1683; 
		covX[2][2] = 6.6635; 
		covX[2][3] = -2.3169;  
		covX[3][1] = -1.8044; 
		covX[3][2] = -2.3169; 
		covX[3][3] = 2.0422;
	}
	
	for (int i=1; i<=NF; i++) bark_mean_goal[i] = mean_goal[g][i];
	hz2bark(bark_mean_goal,bark_mean_goal,NF);
	for (int i=1; i<=NF; i++) new_y[i] = y[i]-bark_mean_goal[i];
	matvecmult(NF,covX,NF,new_y,y_circle);
	norm_y = norm2(y_circle,NF);
	if (norm_y < threshold)
			for (int i=1; i<=NF; i++) normalized_s_err[i] = 0;
	else { 
			for (int i=1; i<=NF; i++) Py[i] = new_y[i]/norm_y*threshold;
			for (int i=1; i<=NF; i++) normalized_s_err[i] = (y[i] - (Py[i]+bark_mean_goal[i]));
	}

	free_dvector(s_err,1,NF);
	free_dvector(bark_mean_goal,1,NF);
	free_dvector(new_y,1,NF);
	free_dvector(y_circle,1,NF);
	free_dvector(Py,1,NF);
	free_dmatrix(covX,1,NF,1,NF);
}

void DistanceProprio(int g, int Nproprio, int NF, double PROPRIO_GOAL_SCALING, double **mean_goal, double *svd_s, double *normalized_s_err)
{
	double test = 0;
	double *s_err = dvector(1,Nproprio);
    double *new_s = dvector(1,Nproprio);
	double *s_circle = dvector(1,Nproprio);
	double *Ps = dvector(1,Nproprio);
	double threshold = 2.8; 
	double **covX = dmatrix(1,Nproprio,1,Nproprio);
	double norm_s;

	if (sqr(mean_goal[g][1]-297)<=1) {
		covX[1][1] = 0.7314; 
		covX[1][2] = 0.0163; 
		covX[1][3] = -0.0904; 
		covX[2][1] = 0.0163; 
		covX[2][2] = 0.9067; 
		covX[2][3] = -0.0174; 
		covX[3][1] = -0.0904; 
		covX[3][2] = -0.0174; 
		covX[3][3] = 0.5109; 
	}
	if (sqr(mean_goal[g][1]-383)<=1) {
		covX[1][1] = 0.4062; 
		covX[1][2] = 0.0012; 
		covX[1][3] = -0.2300; 
		covX[2][1] = -0.0453; 
		covX[2][2] = 0.6040;
		covX[2][3] = -0.0048; 
		covX[3][1] = -0.2257; 
		covX[3][2] = -0.0444; 
		covX[3][3] = 0.5423; 
	}
	if (sqr(mean_goal[g][1]-453)<=1) {
		covX[1][1] = 1.8290; 
		covX[1][2] = -0.2094;
		covX[1][3] = -0.4870; 
		covX[2][1] = -0.2094; 
		covX[2][2] = 2.5701; 
		covX[2][3] = -0.0478; 
		covX[3][1] = -0.4870; 
		covX[3][2] = -0.0478;
		covX[3][3] = 1.3064; 	                                
	}
	if (sqr(mean_goal[g][1]-554)<=1) {
		covX[1][1] = 0.6551; 
		covX[1][2] = -0.1589; 
		covX[1][3] = -0.3498; 
		covX[2][1] = -0.1589; 
		covX[2][2] = 0.4987; 
		covX[2][3] = -0.0991; 
		covX[3][1] = -0.3498; 
		covX[3][2] = -0.0991; 
		covX[3][3] = 0.9343; 	
	}
	if (sqr(mean_goal[g][1]-471)<=1) {
		covX[1][1] = 0.4786; 
		covX[1][2] = -0.1877;
		covX[1][3] = -0.1802; 
		covX[2][1] = -0.1877; 
		covX[2][2] = 0.5246; 
		covX[2][3] = 0.0044; 
		covX[3][1] = -0.1802; 
		covX[3][2] = 0.0044; 
		covX[3][3] = 0.4894; 
	}
	if (sqr(mean_goal[g][1]-449)<=1) {
		covX[1][1] = 0.2743; 
		covX[1][2] = 0.0235; 
		covX[1][3] = -0.0584; 
		covX[2][1] = 0.0235; 
		covX[2][2] = 0.6008; 
		covX[2][3] = -0.0353; 
		covX[3][1] = -0.0584; 
		covX[3][2] = -0.0353; 
		covX[3][3] = 0.5156; 
	}
	if (sqr(mean_goal[g][1]-221)<=1) { //
		covX[1][1] = 0.1570; 
		covX[1][2] = -0.0258; 
		covX[1][3] = -0.0868; 
		covX[2][1] = -0.0258; 
		covX[2][2] = 0.5789; 
		covX[2][3] = 0.0862; 
		covX[3][1] = -0.0868; 
		covX[3][2] = 0.0862; 
		covX[3][3] = 0.3628; 
	}
	if (sqr(mean_goal[g][1]-187)<=1) { 
		covX[1][1] = 0.7384; 
		covX[1][2] = 0.0192; 
		covX[1][3] = 0.5655; 
		covX[2][1] =  0.0192; 
		covX[2][2] = 0.9382; 
		covX[2][3] =  0.0662; 
		covX[3][1] = 0.5655; 
		covX[3][2] =0.0662; 
		covX[3][2] =0.0662; 
		covX[3][3] = 2.5696; 
	}
	
	for (int i=1; i<=Nproprio; i++) new_s[i] = svd_s[i]-mean_goal[g][i+NF];
	matvecmult(Nproprio,covX,Nproprio,new_s,s_circle);
	norm_s = norm2(s_circle,Nproprio);
	if (norm_s < threshold)
			for (int i=1; i<=Nproprio; i++) normalized_s_err[i] = 0;
	else { 
			for (int i=1; i<=Nproprio; i++) Ps[i] = new_s[i]/norm_s*threshold;
			for (int i=1; i<=Nproprio; i++) normalized_s_err[i+NF] = (PROPRIO_GOAL_SCALING*(svd_s[i] - (Ps[i]+mean_goal[g][i+NF])));
    }
	free_dvector(s_err,1,Nproprio);
	free_dvector(new_s,1,Nproprio);
	free_dvector(s_circle,1,Nproprio);
	free_dvector(Ps,1,Nproprio);
	free_dmatrix(covX,1,Nproprio,1,Nproprio);
}


void DistanceEllipse_ds(int g, int NF, int ns, double **mean_goal, double *y, double **Hs, double **Hes)
{
	double test = 0;
	double *bark_mean_goal = dvector(1,NF);
  double *new_y = dvector(1,NF);
	double *y_circle = dvector(1,NF);
	double norm_y;
	double threshold = 2.8;
	double **covX = dmatrix(1,NF,1,NF);
	double **tcovX = dmatrix(1,NF,1,NF);
	double **I3 = dmatrix(1,NF,1,NF);
	double **temp = dmatrix(1,NF,1,NF);
	double **temp2 = dmatrix(1,NF,1,NF);
	double **temp3 = dmatrix(1,NF,1,ns);

	if (sqr(mean_goal[g][1]-297)<=1) {
		covX[1][1] = 0.6125; 
		covX[1][2] = -0.8821; 
		covX[1][3] = -1.1130; 
		covX[2][1] = -0.8821; 
		covX[2][2] = 2.4709; 
		covX[2][3] = 0.2952; 
		covX[3][1] = -1.1130;  
		covX[3][2] = 0.2952; 
		covX[3][3] = 3.4471;
	}
	if (sqr(mean_goal[g][1]-383)<=1) {
		covX[1][1] = 0.5397; 
		covX[1][2] = -1.3515; 
		covX[1][3] = -0.2185; 
		covX[2][1] = -1.3515; 
		covX[2][2] = 3.4581; 
		covX[2][3] = 0.1097; 
		covX[3][1] = -0.2185; 
		covX[3][2] = 0.1097; 
		covX[3][3] = 2.6847;  		                          
	}
	if (sqr(mean_goal[g][1]-453)<=1) {
		covX[1][1] = 1.0960; 
		covX[1][2] = -0.6802; 
		covX[1][3] = 0.1442;
		covX[2][1] = -0.6802; 
		covX[2][2] = 2.6458; 
		covX[2][3] = 0.1429; 
		covX[3][1] = 0.1442; 
		covX[3][2] = 0.1429;
		covX[3][3] = 3.7438;  	                    
    }
	if (sqr(mean_goal[g][1]-554)<=1) {
		covX[1][1] = 1.2526; 
		covX[1][2] = 0.1423; 
		covX[1][3] = -0.7842;  
		covX[2][1] = 0.1423; 
		covX[2][2] = 1.7301; 
		covX[2][3] = 1.0122; 
		covX[3][1] = -0.7842; 
		covX[3][2] = 1.0122; 
		covX[3][3] = 2.6802; 
	}
	if (sqr(mean_goal[g][1]-471)<=1) {
		covX[1][1] = 2.0332; 
		covX[1][2] = 0.4445; 
		covX[1][3] = 0.7377;  
		covX[2][1] = 0.4445; 
		covX[2][2] = 2.7713; 
		covX[2][3] = 1.9709;
		covX[3][1] = 0.7377; 
		covX[3][2] = 1.9709; 
		covX[3][3] = 3.2749; 
	}
	if (sqr(mean_goal[g][1]-449)<=1) {
		covX[1][1] = 2.0578; 
		covX[1][2] = 0.4529;
		covX[1][3] = 0.0088; 
		covX[2][1] = 0.4529;
		covX[2][2] = 2.4506; 
		covX[2][3] = 0.6801; 
		covX[3][1] = 0.0088;
		covX[3][2] = 0.6801; 
		covX[3][3] = 2.8620; 
	}
	if (sqr(mean_goal[g][1]-221)<=1) { 
		covX[1][1] = 1.4072; 
		covX[1][2] = 0.7147; 
		covX[1][3] = -1.5170;
		covX[2][1] = 0.7147;
		covX[2][2] = 4.2991; 
		covX[2][3] = 1.8277; 
		covX[3][1] = -1.5170; 
		covX[3][2] = 1.8277;
		covX[3][3] = 3.3505; 
	}
	if (sqr(mean_goal[g][1]-187)<=1) { 
		covX[1][1] =  14.8396; 
		covX[1][2] =  -0.1875; 
		covX[1][3] =  14.9051; 
		covX[2][1] = -0.1875; 
		covX[2][2] =  3.2284; 
		covX[2][3] =  -0.1042; 
		covX[3][1] =  14.9051;  
		covX[3][2] =  -0.1042;
		covX[3][3] =  14.9730; 
	}
	
	//
	for (int i=1; i<=NF; i++) bark_mean_goal[i] = mean_goal[g][i];
	hz2bark(bark_mean_goal,bark_mean_goal,NF);
	for (int i=1; i<=NF; i++) new_y[i] = y[i]-bark_mean_goal[i];
	matvecmult(NF,covX,NF,new_y,y_circle);
	norm_y = norm2(y_circle, NF);

	if (norm_y < threshold) {
		zeros(Hes,NF,ns);
	}
	else {
		eye(I3,NF);
		todprod(new_y,new_y,NF,temp); 
		transpose(covX,NF,NF,tcovX);
		matmult(NF,temp,NF,tcovX,NF,temp2);
		matmult(NF,temp2,NF,covX,NF,temp);
		scalarmult(temp,NF,NF,1./sqr(norm_y));
		matsub(I3,temp,temp2,NF,NF);
		scalarmult(temp2,NF,NF,threshold/norm_y);
		matsub(I3,temp2,temp,NF,NF);
		matextract(Hs,1,NF,1,ns,temp3);
		matmult(NF,temp,NF,temp3,ns,Hes);
	}

	free_dvector(bark_mean_goal,1,NF);
	free_dvector(new_y,1,NF);
	free_dvector(y_circle,1,NF);
	free_dmatrix(covX,1,NF,1,NF);
	free_dmatrix(tcovX,1,NF,1,NF);
  free_dmatrix(I3,1,NF,1,NF);
	free_dmatrix(temp,1,NF,1,NF);
	free_dmatrix(temp2,1,NF,1,NF);
	free_dmatrix(temp3,1,NF,1,ns);

}


void DistanceProprio_ds(int g, int Nproprio, int NF, int ns, double **mean_goal,double *svd_s, double **Hs, double **Hprop)
{
	double test = 0;
	double *s_err = dvector(1,Nproprio);
  double *new_s = dvector(1,Nproprio);
	double *s_circle = dvector(1,Nproprio);
	double *Ps = dvector(1,Nproprio);
	const double threshold = 2.8; 
	double **covX = dmatrix(1,Nproprio,1,Nproprio);
	double **tcovX = dmatrix(1,Nproprio,1,Nproprio);
	double **I3 = dmatrix(1,Nproprio,1,Nproprio);
	double **temp = dmatrix(1,Nproprio,1,Nproprio);
	double **temp2 = dmatrix(1,Nproprio,1,Nproprio);
	double **temp3 = dmatrix(1,Nproprio,1,ns);
	double norm_s;

	if (sqr(mean_goal[g][1]-297)<=1) {
		covX[1][1] = 1.3982; 
		covX[1][2] = -0.0205; 
		covX[1][3] = 0.2468; 
		covX[2][1] = -0.0205; 
		covX[2][2] = 1.1039; 
		covX[2][3] = 0.0340; 
		covX[3][1] = 0.2468; 
		covX[3][2] = 0.0340; 
		covX[3][3] = 2.0023; 
	}
	if (sqr(mean_goal[g][1]-383)<=1) {
		covX[1][1] = 3.2356; 
		covX[1][2] = 0.0946; 
		covX[1][3] = 1.3733; 
		covX[2][1] = 0.2535; 
		covX[2][2] = 1.6641; 
		covX[2][3] = 0.1222; 
		covX[3][1] = 1.3675; 
		covX[3][2] = 0.1757; 
		covX[3][3] = 2.4256; 
	}
	if (sqr(mean_goal[g][1]-453)<=1) {
		covX[1][1] = 0.6145; 
		covX[1][2] = 0.0544; 
		covX[1][3] = 0.2310; 
		covX[2][1] = 0.0544; 
		covX[2][2] = 0.3942; 
		covX[2][3] = 0.0347; 
		covX[3][1] = 0.2310; 
		covX[3][2] = 0.0347; 
		covX[3][3] = 0.8528; 		                                
	}
	if (sqr(mean_goal[g][1]-554)<=1) {
		covX[1][1] = 2.2451; 
		covX[1][2] = 0.9015; 
		covX[1][3] = 0.9361; 
		covX[2][1] = 0.9015; 
		covX[2][2] = 2.4103; 
		covX[2][3] = 0.5931; 
		covX[3][1] = 0.9361; 
		covX[3][2] = 0.5931; 
		covX[3][3] = 1.4836; 	                               
	}
	if (sqr(mean_goal[g][1]-471)<=1) {
		covX[1][1] = 2.8879; 
		covX[1][2] = 1.0242;
		covX[1][3] = 1.0539; 
		covX[2][1] = 1.0242; 
		covX[2][2] = 2.2695; 
		covX[2][3] = 0.3566;
		covX[3][1] = 1.0539; 
		covX[3][2] = 0.3566; 
		covX[3][3] = 2.4279; 
	}
	if (sqr(mean_goal[g][1]-449)<=1) {
		covX[1][1] = 3.7451; 
		covX[1][2] = -0.1222;
		covX[1][3] = 0.4158; 
		covX[2][1] = -0.1222; 
		covX[2][2] = 1.6752;
		covX[2][3] = 0.1009;
		covX[3][1] = 0.4158;
		covX[3][2] = 0.1009; 
		covX[3][3] = 1.9934; 
	}
	if (sqr(mean_goal[g][1]-221)<=1) { //
		covX[1][1] = 4.5915; 
		covX[1][2] = 0.4124; 
		covX[1][3] = 0.5140; 
		covX[2][1] = 0.4124; 
		covX[2][2] = 1.6329; 
		covX[2][3] = -0.4595; 
		covX[3][1] = 0.5140; 
		covX[3][2] = -0.4595; 
		covX[3][3] = 2.4342; 
	}
	if (sqr(mean_goal[g][1]-187)<=1) { 
		covX[1][1] = 1.6287; 
		covX[1][2] = -0.0080; 
		covX[1][3] = -0.3582; 
		covX[2][1] = -0.0080; 
		covX[2][2] = 1.0679; 
		covX[2][3] = -0.0258; 
		covX[3][1] = -0.3582; 
		covX[3][2] = -0.0258;
		covX[3][3] = 0.4687;
    }
	
	for (int i=1; i<=Nproprio; i++) new_s[i] = svd_s[i]-mean_goal[g][i+NF];
	matvecmult(Nproprio,covX,Nproprio,new_s,s_circle);
	norm_s = norm2(s_circle,Nproprio);
	if (norm_s < threshold)
			zeros(Hprop,Nproprio,ns);
	else { 
			eye(I3,Nproprio);
			todprod(svd_s,svd_s,Nproprio,temp); 
			transpose(covX,Nproprio,Nproprio,tcovX);
			matmult(Nproprio,temp,Nproprio,tcovX,Nproprio,temp2);
			matmult(Nproprio,temp2,Nproprio,covX,Nproprio,temp);
			scalarmult(temp,Nproprio,Nproprio,1./sqr(norm_s));
			matsub(I3,temp,temp2,Nproprio,Nproprio);
			scalarmult(temp2,Nproprio,Nproprio,threshold/norm_s);
			matsub(I3,temp2,temp,Nproprio,Nproprio);
			matextract(Hs,NF+1,NF+Nproprio,1,ns,temp3); 
			matmult(Nproprio,temp,Nproprio,temp3,ns,Hprop);
	}

	free_dvector(s_err,1,Nproprio);
	free_dvector(new_s,1,Nproprio);
	free_dvector(s_circle,1,Nproprio);
	free_dvector(Ps,1,Nproprio);
	free_dmatrix(covX,1,Nproprio,1,Nproprio);
	free_dmatrix(tcovX,1,Nproprio,1,Nproprio);
    free_dmatrix(I3,1,Nproprio,1,Nproprio);
	free_dmatrix(temp,1,Nproprio,1,Nproprio);
	free_dmatrix(temp2,1,Nproprio,1,Nproprio);
	free_dmatrix(temp3,1,Nproprio,1,ns);

}


bool autoencoder_init()
{
  mxArray *pW1, *pW2, *pb1, *pb2, *pm, *pM;
  double *psvdM1, *psvdM2;
  int sz = sizeof(double);

  pW1= mexGetVariable("global", "WE");
  if (!pW1) {mexErrMsgTxt("WE: matrice pas trouvee\n"); }
  pW2= mexGetVariable("global", "WD");
  if (!pW2) {mexErrMsgTxt("WD : matrice pas trouvee\n"); }
  pb1 = mexGetVariable("global", "bE");
  if (!pb1) {mexErrMsgTxt("bE : vecteur pas trouve\n"); }
  pb2 = mexGetVariable("global", "bD");
  if (!pb2) {mexErrMsgTxt("bD : vecteur pas trouve\n"); }
  pm = mexGetVariable("global", "xmin");
  if (!pm) {mexErrMsgTxt("xmin : vecteur pas trouve\n"); }
  pM = mexGetVariable("global", "xmax");
  if (!pM) {mexErrMsgTxt("xmax : vecteur pas trouve\n"); }

  WE = dmatrix(1,DOF,1,2*N_NODES);
  psvdM1 = mxGetPr(pW1);
  matconvert(psvdM1, WE, DOF, 2*N_NODES);
  WD = dmatrix(1,2*N_NODES,1,DOF);
  psvdM2 = mxGetPr(pW2);
  matconvert(psvdM2, WD,2*N_NODES, DOF);

  bE = dvector(1,DOF);
  memcpy(bE+1, mxGetPr(pb1), DOF*sz);
  bD = dvector(1,2*N_NODES);
  memcpy(bD+1, mxGetPr(pb2), 2*N_NODES*sz);
  m = dvector(1,2*N_NODES);
  memcpy(m+1, mxGetPr(pm), 2*N_NODES*sz);
  M = dvector(1,2*N_NODES);
  memcpy(M+1, mxGetPr(pM), 2*N_NODES*sz);

  // Jacobian of decoding fn
  JD = dmatrix(1,2*N_NODES,1,DOF);
  
  double **JDunscaled, **scalingmat, *valrange;
  JDunscaled = dmatrix(1,2*N_NODES,1,DOF);
  scalingmat = dmatrix(1,2*N_NODES,1,2*N_NODES);
  valrange = dvector(1,2*N_NODES);
  
  copy(JDunscaled,WD,2*N_NODES,DOF);
  vecsub(M,m,valrange, 2*N_NODES);
  scalarmult(valrange, 2*N_NODES, 1./2);
  diag(scalingmat,valrange,2*N_NODES);
  matmult(2*N_NODES,scalingmat,2*N_NODES,JDunscaled,DOF,JD);

  free_dmatrix(JDunscaled,1,2*N_NODES,1,2*N_NODES);
  free_dmatrix(scalingmat,1,2*N_NODES,1,2*N_NODES);
  free_dvector(valrange,1,2*N_NODES);


  autoenc_initialized = true;
  return(autoenc_initialized);
}

void free_variables()
{
  free_dmatrix(WE,1,DOF,1,2*N_NODES);
  free_dmatrix(WD,1,2*N_NODES,1,DOF);
  free_dmatrix(JD,1,2*N_NODES,1,DOF);
  free_dvector(bE,1,DOF);
  free_dvector(bD,1,2*N_NODES);
  free_dvector(m,1,2*N_NODES);
  free_dvector(M,1,2*N_NODES);

}
