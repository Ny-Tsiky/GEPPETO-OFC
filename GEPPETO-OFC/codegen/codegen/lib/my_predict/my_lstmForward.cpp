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
#include "my_lstmForward.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include "mex.h"  // comment out for pure C++

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


extern int MAX_STEPS;
extern int nn;
extern int LSTM_SIZE;
//extern double *saved_CS, *saved_HS;

void lstmForwardUsingExplicitLoops(const float X[], const float W[],
                                   const float R[], const float gate_bias[],
                                   double **CS, double **HS, int step)
{
  float G[4*LSTM_SIZE];
  float cellGateOp[LSTM_SIZE];
  float forgetGateOp[LSTM_SIZE];
  float ipGateOp[LSTM_SIZE];
  float outputGateOp[LSTM_SIZE];
  float f;
  float f1;
  float f2;
  float f3;
  float f4;
  int b_i;
  int i;
  int p;
  std::memset(G, 0, 4*LSTM_SIZE * sizeof(float));
//#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
    p, f, f1, f2, f3, b_i, f4)

if (step == 1) {
	memset(HS[1], 0, LSTM_SIZE * sizeof(double));
  memset(CS[1], 0, LSTM_SIZE * sizeof(double));
}

  char errmsg[128];
  sprintf(errmsg, "LSTM_forward : etat non sauvegarde (pas de temps : %d, max = %d)\n", step, MAX_STEPS);
  if (step > MAX_STEPS) {mexErrMsgTxt(errmsg); return;}

/*
  for (i = 0; i < LSTM_SIZE; i++) {
    f = G[i];
    f1 = G[i + LSTM_SIZE];
    f2 = G[i + 2*LSTM_SIZE];
    f3 = G[i + 3*LSTM_SIZE];
    for (p = 0; p < nn; p++) {
      b_i = i + (p << 6);
      f4 = X[p];
      f += W[b_i] * f4;
      f1 += W[b_i + LSTM_SIZE] * f4;
      f2 += W[b_i + 2*LSTM_SIZE] * f4;
      f3 += W[b_i + 3*LSTM_SIZE] * f4;
    }
    for (p = 0; p < LSTM_SIZE; p++) {
      b_i = i + (p << 6);
      f4 = HS[step][p];
      f += R[b_i] * f4;
      f1 += R[b_i + LSTM_SIZE] * f4;
      f2 += R[b_i + 2*LSTM_SIZE] * f4;
      f3 += R[b_i + 3*LSTM_SIZE] * f4;
    }
    G[i] = 1.0F / (std::exp(-(f + b[i])) + 1.0F);
    G[i + LSTM_SIZE] = 1.0F / (std::exp(-(f1 + b[i + LSTM_SIZE])) + 1.0F);
    G[i + 2*LSTM_SIZE] = std::tanh(f2 + b[i + 2*LSTM_SIZE]);
    G[i + 3*LSTM_SIZE] = 1.0F / (std::exp(-(f3 + b[i + 3*LSTM_SIZE])) + 1.0F);
  }
  std::copy(&G[0], &G[LSTM_SIZE], &ipGateOp[0]);
  std::copy(&G[LSTM_SIZE], &G[2*LSTM_SIZE], &forgetGateOp[0]);
  std::copy(&G[2*LSTM_SIZE], &G[3*LSTM_SIZE], &cellGateOp[0]);
  std::copy(&G[3*LSTM_SIZE], &G[4*LSTM_SIZE], &outputGateOp[0]);
//#pragma omp parallel for num_threads(omp_get_max_threads()) private(f)

  for (i = 0; i < LSTM_SIZE; i++) {
    f = cellGateOp[i] * ipGateOp[i] + forgetGateOp[i] * CS[step][i];// c0[i];
    CS[step+1][i] = f;
    HS[step+1][i] = std::tanh(f) * outputGateOp[i];
  }
  */
  
  
   for (i = 0; i < LSTM_SIZE; i++) {
    f = 0;
    f1 = 0;
    f2 = 0;
    f3 = 0;
    for (p = 0; p < 10; p++) {
      //b_i = i + (p << 6);
	  b_i = i + 4*LSTM_SIZE * p;
      f4 = X[p];
      f += W[b_i] * f4;
      f1 += W[b_i + LSTM_SIZE] * f4;
      f2 += W[b_i + 2*LSTM_SIZE]*f4;
      f3 += W[b_i + 3*LSTM_SIZE] * f4;
    }
    for (p = 0; p < LSTM_SIZE; p++) {
      //b_i = i + (p << 6);
	    b_i = i + 4*LSTM_SIZE * p;
      f4 = HS[step][p];
      f += R[b_i] * f4;
      f1 += R[b_i + LSTM_SIZE] * f4;
      f2 += R[b_i + 2*LSTM_SIZE]*f4;
      f3 += R[b_i + 3*LSTM_SIZE] * f4;
    }
    G[i] = 1.0F / (std::exp(-(f + gate_bias[i])) + 1.0F);
    G[i + LSTM_SIZE] = 1.0F / (std::exp(-(f1 + gate_bias[i + LSTM_SIZE])) + 1.0F);
    G[i + 2*LSTM_SIZE] = std::tanh(f2 + gate_bias[i + 2*LSTM_SIZE]);
    G[i + 3*LSTM_SIZE] = 1.0F / (std::exp(-(f3 + gate_bias[i + 3*LSTM_SIZE])) + 1.0F);
  
    f = G[i+2*LSTM_SIZE] * G[i]  + G[i+LSTM_SIZE] * CS[step][i];// c0[i];
    CS[step+1][i] = f;
    HS[step+1][i] = std::tanh(f) * G[i + 3*LSTM_SIZE];
  }
}


//
// File trailer for lstmForward.cpp
//
// [EOF]
//
