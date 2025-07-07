
enum noise {GAUSSIAN, POISSON, SDN, MIXED};
typedef enum kalmantype {DISCRETE, CONTINUOUS, DISCRETE_EKF, CONTINUOUS_EKF} kalman_type;


typedef struct kalman_par {
  enum noise noise_type;  // type of multiplicative noise (Poisson or SDN)
  double *delay;   //  in time steps
  double dt;
  void (*Afun)(double *x, double t, double **A);
  void (*Bfun)(double *x,double *u,double t, double **B);
  void (*Hfun)(double *x,double t, double **H);
  double sigQ;  // additive motor noise variance scaling factor
  double sigQp; // multiplicative motor noise variance scaling factor
  double *sigR;  // additive sensory ... now a vector, variance per modality
  double *sigRp; // multiplicative sensory ...
  double sigS;  // additive process noise variance scaling factor
  int nx;
  int ny;
  int Nmodal;
  int *ny_modal; // vector sizes of modality-specific FB
  int nu;
} kalman_params;


typedef struct ekf_par {
  enum noise noise_type;
  double *delay; // now a vector of modality-specific delays
  double dt;
  void (*f)(double* x, double* old_u, double* u, double t, double* xdot, double *out_HS, double *out_CS); //old_u to compute u_{i+1}-u_i
  void (*dfdx)(double* x, double* u, double t, double** Fx);
  void (*dfdu)(double* x, double* old_u, double* u, double t, double** Fu);
  void (*h)(double *x, double *u, double t, double *z);
  void (*dhdx)(double *x, double *u, double t, double **H);
  double sigQ;  // additive motor noise variance scaling factor
  double sigQp; // multiplicative motor noise variance scaling factor
  double *sigR;  // additive sensory ... now a vector, variance per modality
  double *sigRp; // multiplicative sensory ...
  //double sigR;  // additive sensory ...
  //double sigRp; // multiplicative sensory ...
  double sigS;  // additive process noise variance scaling factor
  int nx;
  int ny;
  int Nmodal;
  int *ny_modal; // vector sizes of modality-specific FB
  int nu;
} ekf_params;


void init_kalman_variables_static(double *s, double *ustart, double *ystart);
void init_kalman(kalman_params *params);
void update_kalman(kalman_params *params);
void reset_kalman();
void kalman(double *u, double *y, double *xhat, int *consider_feedback, kalman_params *params=NULL, int reset=0);
void free_kalman_variables();

void init_ekf_variables(double **txest, double **tustart, double **tystart);
void init_ekf_variables_static(double *xest, double *ustart, double *ystart, double **Estart=NULL);
void init_ekf(ekf_params *params); //, double *xest, double *ustart, double *ystart);
void update_ekf_noiseparams(ekf_params *params);
void update_ekf(ekf_params *params);
void get_ekf_params(ekf_params *params);

void init_discrete_ekf_variables(double **kxest, double **kustart, double **kystart);
void init_discrete_ekf_variables_static(double *xest, double *ustart, double *ystart, double **Estart=NULL);
void set_LSTM_states(double t, double dt, double *HS, double*CS, int actual);

//void reset_ekf();
void ekf(double t, double *u, double *y, double *xhat, int *consider_feedback, ekf_params *params=NULL,  int reset=0);
void ekf(double t, double *u, double *y, double *xest, double time_ahead, double *xest_ahead, int *consider_feedback, ekf_params *params=NULL, int reset=0);
void discrete_ekf(double t, double *u, double *y, double *xest, int *consider_feedback, double *out_HS, double *out_CS, ekf_params *params=NULL, int reset=0);
void discrete_ekf(double t, double *u, double *y, double *xest, int *consider_feedback, double *yest, double **E, ekf_params *params=NULL, int reset=0);
void free_ekf_variables();
void free_discrete_ekf_variables();




