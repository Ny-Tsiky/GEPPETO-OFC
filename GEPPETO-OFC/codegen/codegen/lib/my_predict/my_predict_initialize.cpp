//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: my_predict_initialize.cpp
//
// MATLAB Coder version            : 5.3
// C/C++ source code generated on  : 29-Dec-2021 17:47:39
//

// Include Files
#include "my_predict_initialize.h"
#include "my_predict.h"
#include "my_predict_data.h"

// Function Definitions
//
// Arguments    : void
// Return Type  : void
//
void my_predict_initialize()
{
  //omp_init_nest_lock(&my_predict_nestLockGlobal);
  my_predict_init();
  isInitialized_my_predict = true;
}

//
// File trailer for my_predict_initialize.cpp
//
// [EOF]
//
