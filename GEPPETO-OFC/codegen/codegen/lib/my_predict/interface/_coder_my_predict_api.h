//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: _coder_my_predict_api.h
//
// MATLAB Coder version            : 5.3
// C/C++ source code generated on  : 29-Dec-2021 17:47:39
//

#ifndef _CODER_MY_PREDICT_API_H
#define _CODER_MY_PREDICT_API_H

// Include Files
#include "emlrt.h"
#include "tmwtypes.h"
#include <algorithm>
#include <cstring>

// Variable Declarations
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

// Function Declarations
void my_predict(real_T in[10], real32_T out[4]);

void my_predict_api(const mxArray *prhs, const mxArray **plhs);

void my_predict_atexit();

void my_predict_initialize();

void my_predict_terminate();

void my_predict_xil_shutdown();

void my_predict_xil_terminate();

#endif
//
// File trailer for _coder_my_predict_api.h
//
// [EOF]
//
