//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: fullyConnectedForward.h
//
// MATLAB Coder version            : 5.3
// C/C++ source code generated on  : 29-Dec-2021 17:47:39
//

#ifndef FULLYCONNECTEDFORWARD_H
#define FULLYCONNECTEDFORWARD_H

// Include Files
#include "rtwtypes.h"
#include "omp.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace layer {
void b_iAddBiasApplyActivation(float X[8], const float bias[8]);

void c_iAddBiasApplyActivation(float X[4], const float bias[4]);

void iAddBiasApplyActivation(float X[16], const float bias[16]);

} // namespace layer
} // namespace coder

#endif
//
// File trailer for fullyConnectedForward.h
//
// [EOF]
//
