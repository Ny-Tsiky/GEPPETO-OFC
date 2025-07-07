//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: lstmForward.h
//
// MATLAB Coder version            : 5.3
// C/C++ source code generated on  : 29-Dec-2021 17:47:39
//

#ifndef LSTMFORWARD_H
#define LSTMFORWARD_H

// Include Files
#include "rtwtypes.h"
#include "omp.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace layer {
void lstmForwardUsingExplicitLoops(const float X[10], const float W[640],
                                   const float R[1024], const float b[64],
//                                   const float c0[16], const float b_y0[16],
                                   float CS[16], float YT[16]);

}
} // namespace coder

#endif
//
// File trailer for lstmForward.h
//
// [EOF]
//
