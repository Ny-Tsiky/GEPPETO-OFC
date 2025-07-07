#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include "nrutil.h"
#include "mylinpack.h"
#include "security.h"

/* returns the vector of n regression coefs plus (index 0) the determination coefficient R*/


double *LinearRegression(int n_params, int n_samples, double to_fit[],...)
{
  va_list argptr;
  double **params;
  int i,j,k;
  double  **mult, **dum, *dummy, *result, *fit, *mean, variance, S=0;

  n_params++;  /* constant parameter added */

  params=(double **) Calloc(n_params+NR_END,sizeof(double*));
  params[n_params]=(double *) Calloc(n_samples+NR_END,sizeof(double));

  va_start(argptr,to_fit);
  for (i=1; i<n_params; i++) params[i]=va_arg(argptr,double*)-NR_END;
  for (j=1; j<=n_samples; j++) params[n_params][j]=1;

  /* params = X' in Saporta */
  
  to_fit--;  /* so to_fit[1] is the first coord */

  mult=dmatrix(1,n_params,1,n_params);
  dum=dmatrix(1,n_params,1,n_samples);
  dummy=dvector(1,n_samples); 
  mean=dvector(1,n_samples);
  fit=dvector(1,n_samples);
  result=dvector(0,n_params);




  for (i=1; i<=n_params; i++) 
    for (j=1; j<=n_params; j++){
      mult[i][j]=0; 
      for (k=1; k<=n_samples; k++) mult[i][j]+=params[i][k]*params[j][k];
    }

  gaussjordan(mult,n_params,dum,1);

  for (i=1; i<=n_params; i++) {
    result[i]=0;
    for (j=1; j<=n_samples; j++) {
      dum[i][j]=0;
      for (k=1; k<=n_params; k++) dum[i][j]+=mult[i][k]*params[k][j];
      result[i]+=dum[i][j]*to_fit[j];  
    }
  }

  for (i=1; i<=n_samples; i++) {
      fit[i]=0; 
      for (k=1; k<=n_params; k++) fit[i]+=params[k][i]*result[k];
    }

  for (i=1; i<=n_samples; i++) S+=to_fit[i];
  S/=n_samples;
  for (i=1; i<=n_samples; i++) mean[i]=S;

  vectsub(mean,to_fit,n_samples,dummy);
  vectsub(mean,fit,n_samples,mean);
  variance=squared(dummy,n_samples);
  if (variance>EPSILON) result[0]=sqrt(squared(mean,n_samples)/variance);
  else result[0]=1;

  free_dmatrix(mult,1,n_params,1,n_params);
  free_dmatrix(dum,1,n_params,1,n_samples);
  free_dvector(dummy,1,n_samples);
  free_dvector(mean,1,n_samples);
  free_dvector(fit,1,n_samples);
  free(params[n_params]);
  free(params); 

  va_end(argptr);
  return(result);
}


