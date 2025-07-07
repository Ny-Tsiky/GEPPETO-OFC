//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: fullyConnectedForward.cpp
//
// MATLAB Coder version            : 5.3
// C/C++ source code generated on  : 29-Dec-2021 17:47:39
//

// Include Files
#include "fullyConnectedForward.h"

// Function Definitions
//
// Arguments    : float X[8]
//                const float bias[8]
// Return Type  : void
//
namespace coder {
namespace layer {
void b_iAddBiasApplyActivation(float X[8], const float bias[8])
{
  int v1;
//#pragma omp parallel for num_threads(omp_get_max_threads()) private(v1)

  for (int iElem = 0; iElem < 8; iElem++) {
    v1 = iElem - ((iElem / 8) << 3);
    X[v1] += bias[v1];
  }
}

//
// Arguments    : float X[4]
//                const float bias[4]
// Return Type  : void
//
void c_iAddBiasApplyActivation(float X[4], const float bias[4])
{
  X[0] += bias[0];
  X[1] += bias[1];
  X[2] += bias[2];
  X[3] += bias[3];
}

//
// Arguments    : float X[16]
//                const float bias[16]
// Return Type  : void
//
void iAddBiasApplyActivation(float X[16], const float bias[16])
{
  int v1;
//#pragma omp parallel for num_threads(omp_get_max_threads()) private(v1)

  for (int iElem = 0; iElem < 16; iElem++) {
    v1 = iElem - ((iElem / 16) << 4);
    X[v1] += bias[v1];
  }
}

} // namespace layer
} // namespace coder

//
// File trailer for fullyConnectedForward.cpp
//
// [EOF]
//
