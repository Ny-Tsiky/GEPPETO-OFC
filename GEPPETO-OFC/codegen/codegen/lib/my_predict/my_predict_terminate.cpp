//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: my_predict_terminate.cpp
//
// MATLAB Coder version            : 5.3
// C/C++ source code generated on  : 29-Dec-2021 17:47:39
//

// Include Files
#include "my_predict_terminate.h"
#include "my_predict.h"
#include "my_predict_data.h"

// Function Definitions
//
// Arguments    : void
// Return Type  : void
//
void my_predict_terminate()
{
  my_predict_free();
  omp_destroy_nest_lock(&my_predict_nestLockGlobal);
  isInitialized_my_predict = false;
}

//
// File trailer for my_predict_terminate.cpp
//
// [EOF]
//
