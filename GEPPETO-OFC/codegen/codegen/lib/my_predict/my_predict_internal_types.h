//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: my_predict_internal_types.h
//
// MATLAB Coder version            : 5.3
// C/C++ source code generated on  : 29-Dec-2021 17:47:39
//

#ifndef MY_PREDICT_INTERNAL_TYPES_H
#define MY_PREDICT_INTERNAL_TYPES_H

// Include Files
//#include "my_predict_types.h"
#include "rtwtypes.h"

// Type Definitions
struct cell_wrap_3 {
  float f1[16];
};

struct cell_wrap_4 {
  cell_wrap_3 f1[2];
};

namespace coder {
namespace ctarget {
struct coder_ctarget_DeepLearningNetwork_tag_0 {
  boolean_T matlabCodegenIsDeleted;
  //cell_wrap_4 State[1];
  float HS[16];
  float CS[16];
  boolean_T IsInitialized;
};
typedef coder_ctarget_DeepLearningNetwork_tag_0 DeepLearningNetwork;

} // namespace ctarget
} // namespace coder

#endif
//
// File trailer for my_predict_internal_types.h
//
// [EOF]
//
