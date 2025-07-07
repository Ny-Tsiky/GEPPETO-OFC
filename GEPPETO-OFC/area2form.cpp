#include <math.h>
#include "mylinpack.h"
#include "nr.h"
#include "nrutil.h"
#include "tongue_geom_audio_autoenc.h" 

#include "mex.h" 

double *averagetongue; 
int NF = 3;

void area2formants(double *tube_area, double tube_length, double *formants)
{
	VTarea2formantLocation(tube_area, tube_length, formants);
	}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *tube_area;
	double tube_length; 
	double *formants = dvector(1,NF);

  tube_area = mxGetPr(prhs[0]); 
  tube_length = *mxGetPr(prhs[1]);  
	
  plhs[0] = mxCreateDoubleMatrix(NF,1,mxREAL);
  formants = mxGetPr(plhs[0])-1;

  area2formants(tube_area, tube_length, formants); 
}

