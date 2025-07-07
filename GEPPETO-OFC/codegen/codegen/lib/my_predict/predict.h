//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: predict.h
//
// MATLAB Coder version            : 5.3
// C/C++ source code generated on  : 29-Dec-2021 17:47:39
//

#ifndef PREDICT_H
#define PREDICT_H

// Include Files
#include "rtwtypes.h"
#include "omp.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
namespace coder {
namespace ctarget {
typedef struct coder_ctarget_DeepLearningNetwork_tag_0 DeepLearningNetwork;

}
} // namespace coder

// Function Declarations
namespace coder {
namespace ctarget {
void DeepLearningNetwork_predict(DeepLearningNetwork *obj,
                                 const double varargin_1[10],
                                 double varargout_1[4]);

}
} // namespace coder

#endif
//
// File trailer for predict.h
//
// [EOF]
//
