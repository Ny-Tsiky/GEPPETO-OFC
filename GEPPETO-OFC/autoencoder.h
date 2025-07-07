const int DOF = 5;

void satlin(int n,double *p);

void encode(double *p,double *result);

void decode(double *p,double *result);

void map_minmax(double *v, const int n, const double *m, const double *M);

void unmap_minmax(double *v, const int n, const double *m, const double *M);

void matconvert(double *pA, double **A, int m, int n);

void DistanceEllipse(int g, int NF, double **mean_goal, double *y, double *s_err);
void DistanceProprio(int g, int Nproprio, int NF, double PROPRIO_GOAL_SCALING, double **mean_goal, double *svd_s, double *normalized_s_err);
void DistanceEllipse_ds(int g,int NF, int ns, double **mean_goal, double *y, double **Hs, double **Hes);
void DistanceProprio_ds(int g, int Nproprio, int NF, int ns, double **mean_goal,double *svd_s, double **Hs, double **Hprop);

double Hz2bark(double fHz);
void hz2bark(double *y,double *bark_y,int NF);
 
void bark2hz(double *bark_y,double *y,int NF);

bool autoencoder_init();
