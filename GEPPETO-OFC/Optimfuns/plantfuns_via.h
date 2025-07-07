typedef struct noise_par {
	double var_u_sdn;
	double var_s1_sdn; 
	double var_s2_sdn; 
	double *var_s_sdn;
	double var_u_add;
	double var_s1_add;  
	double var_s2_add;  
	double *var_s_add;
	double var_x_add;
	int dim_s2;         
	int Nmodal;
  
  double *y_perturb;  
  double t_onset, t_full;
} noise_params;



void calc_dynamics(double *s);
void funback(double *u, double *s, double t, double **f1, double **f2);
void firsthalf_flatfunback(double *u, double *s, double t, double *f1);
void secondhalf_flatfunback(double *u, double *s, double t, double *f2);
void firsthalffunback(double *u, double *s, double t, double **f1);
void halffunback(double *u, double *s, double t, double **f2);
void funforward(double *u, double *s, double t, double *f1);
void funforward2(double *u, double *s, double t, double *f1);

void calc_Phis(double *s, double *f1, double **f2); 
void calc_Phist(double *s, double t, double *f1, double **f2, double *f3);
void firsthalf_calc_flatPhis(double *s, double *f1);
void secondhalf_calc_flatPhis(double *s, double *f2);
void firstpart_calc_flatPhist(double *s, double t, double *f1);
void secondpart_calc_flatPhist(double *s, double t, double *f2);
void thirdpart_calc_flatPhist(double *s, double t, double *f3);

bool mexinit(int *pns, int *pnp, int *pnc);
void init(int *pns, int *pnp, int *pnc, int Nsteps, int kmax, int Nmax, double *thi, double *thf, double *tbds, double ***ptu, double ***pts, double **ps0, double **pxp, double ***pyp);
void init_arm(int *pns, int *pnp, int *pnc, int Nsteps, double *thi, double ***pu, double **ps0);
void init_cost_constraints(double *constr);
void init_state_constraints(double *constr);
void free_state_constraints_variables();
void free_plant_variables();
void free_arm_variables();
void update_kalman_filter(double *u, double *sest, double *y, double dt);
void update_kalman_noises(noise_params *noise, kalman_type kalmtype);
void init_kalman_filter(double *s, noise_params *noise, double dt, double *delay, int *ny_modal, kalman_type kalmtype);
void init_kalman_filter(double *s, int *pny, double var_u, double var_s, double dt, double delay);
void init_kalman_filter(double *s, noise_params *noise, double dt, double *delay, kalman_type kalmtype);
void init_kalman_filter(double *s, int *pny, double var_u, double var_s, double var_u_add, double var_s_add, double dt, double delay);  
void init_kalman_filter(double *s, int *pny, double var_u, double var_s, double var_u_add, double var_s_add, double var_x_add, double dt, double delay);  
void init_kalman_filter(double *s, int *pny, noise_params *noise, double dt, double delay);  
void init_kalman_filter(double *s, noise_params *noise, double dt, double delay, kalman_type kalmtype);  
void init_stupid_filter(double *s, int *pny, double dt, double delay);
void init_kalman_filter_test(int *pny, double var_u, double var_s, double dt, double delay, enum noise noise_type);
void free_kalman_filter(kalman_type kalmtype);
void free_kalman_filter();
void xinit_kalman_filter(double *s, int *pnx, int *pny, double var_u, double var_s, double dt);

void discrete_forward(double *old_u, double *u, double *s, double t, double dt, double *f1, int actual, double *out_HS, double *out_CS);
void discrete_forward(double *old_u, double *u, double *s, double t, double dt, double *f1, int actual, double *out_HS, double *out_CS, int perturb);
void discrete_forward_biomeca(double *old_u, double *u, double *s, double t, double dt, double *f1);
void discrete_constraints(double *u, double *s, double *cstr, double t, double dt, double *f1);
void discrete_constraints_ds(double *u, double *s, double *cstr, double t, double dt, double *f1);
void discrete_adjoint_dx(double *u, double *s, double t, double dt, double *f1, int actual);
void discrete_adjoint_du(double *old_u, double *u, double *s, double t, double dt, double *f2);
void discrete_adjoint_dx2(double *u, double *s, double t, double dt, double **f1, int actual);
void discrete_adjoint_du2(double *u_old,double *u, double *s, double t, double dt, double **f2, int actual=0);
void discrete_constraints(double *s, double *f1);
void discrete_constraints_dx(double *s, double *f2);
void discrete_constraints_dx2(double *s, double **f2);

void statetransition_matrix(double t, double **A);
void control_matrix(double t, double **B);
void observation_matrix(double t, double **H);

void discrete_via_forward(double *u, double *s, double t, double dt, double *via_t, double *f1);
void discrete_via_adjoint_dx(double *u, double *s, double t, double dt, double *via_t, double *f1);
void init_via_state_constraints(double *constr, int n_viapts);
void pass_via_times_ptr(double *viat);

void discrete_cost(double *s, double *f1, double *sens_cost, double *tact_cost, double *neuro_cost); 
void discrete_cost_dx(double *s, double *f2);
void init_goals(double *goalseq, int Ndims, int Ng);
void init_goals(double *goalseq, char *goalcodeseq, int Ndims, int Ng);

void init_lstm_dop0viafuns(double *HS_Kalman, double *CS_Kalman);
void init_lstm_simnoisy(double *HS_Kalman, double *CS_Kalman, double *HS_LSTM, double *CS_LSTM); 
void free_goals_variables();

void dirkinematics(double *theta, double*x);
int invkinematics(double *x, double *theta);
void jacobian(double *theta, double **J);
void calc_feedbackforce(double *s, double *extfce);
void calc_feedbackderivs(double *s, double **derivs);
void init_environment_variables();
void getforcetorques(double *s, double *ft);
void getforcetorques(double *s, double t, double *ft);

