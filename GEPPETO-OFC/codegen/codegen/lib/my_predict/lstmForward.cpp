//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: lstmForward.cpp
//
// MATLAB Coder version            : 5.3
// C/C++ source code generated on  : 29-Dec-2021 17:47:39
//

// Include Files
#include "lstmForward.h"
#include <algorithm>
#include <cmath>
#include <cstring>

// Function Definitions
//
// Arguments    : const float X[10]
//                const float W[640]
//                const float R[1024]
//                const float b[64]
//                const float c0[16]
//                const float b_y0[16]
//                float CS[16]
//                float YT[16]
// Return Type  : void
//
namespace coder {
namespace layer {
void lstmForwardUsingExplicitLoops(const float X[10], const float W[640],
                                   const float R[1024], const float b[64],
//                                   const float c0[16], const float b_y0[16],
                                   float CS[16], float YT[16])
{
  float G[64];
  float cellGateOp[16];
  float forgetGateOp[16];
  float ipGateOp[16];
  float outputGateOp[16];
  float f;
  float f1;
  float f2;
  float f3;
  float f4;
  int b_i;
  int i;
  int p;
  std::memset(&G[0], 0, 64U * sizeof(float));
//#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
    p, f, f1, f2, f3, b_i, f4)

  for (i = 0; i < 16; i++) {
    f = G[i];
    f1 = G[i + 16];
    f2 = G[i + 32];
    f3 = G[i + 48];
    for (p = 0; p < 10; p++) {
      b_i = i + (p << 6);
      f4 = X[p];
      f += W[b_i] * f4;
      f1 += W[b_i + 16] * f4;
      f2 += W[b_i + 32] * f4;
      f3 += W[b_i + 48] * f4;
    }
    for (p = 0; p < 16; p++) {
      b_i = i + (p << 6);
      f4 = YT[p];
      f += R[b_i] * f4;
      f1 += R[b_i + 16] * f4;
      f2 += R[b_i + 32] * f4;
      f3 += R[b_i + 48] * f4;
    }
    G[i] = 1.0F / (std::exp(-(f + b[i])) + 1.0F);
    G[i + 16] = 1.0F / (std::exp(-(f1 + b[i + 16])) + 1.0F);
    G[i + 32] = std::tanh(f2 + b[i + 32]);
    G[i + 48] = 1.0F / (std::exp(-(f3 + b[i + 48])) + 1.0F);
  }
  std::copy(&G[0], &G[16], &ipGateOp[0]);
  std::copy(&G[16], &G[32], &forgetGateOp[0]);
  std::copy(&G[32], &G[48], &cellGateOp[0]);
  std::copy(&G[48], &G[64], &outputGateOp[0]);
//#pragma omp parallel for num_threads(omp_get_max_threads()) private(f)

  for (i = 0; i < 16; i++) {
    f = cellGateOp[i] * ipGateOp[i] + forgetGateOp[i] * CS[i];// c0[i];
    CS[i] = f;
    YT[i] = std::tanh(f) * outputGateOp[i];
  }
}

} // namespace layer
} // namespace coder

//
// File trailer for lstmForward.cpp
//
// [EOF]
//
