const double EPS = 1e-10;

const int N_NODES = 16; 
const int N_VTSECTIONS = 44;

bool tongue2VTarea(double *tong, double *VTarea, double &sectLen);
bool tongue2VTarea(double *tong, double *VTarea, double &sectLen, double*slicing_diam, double *, double *);
void VTarea2formantLocation(double *vtArea, double vtTubeLength, double *F);
void VTarea2formantLocationAndGuess(double *vtArea, double vtTubeLength, double *Fguess, double *F);

void heapsort(int *x, int N);
void heapsort(double *x, int N);
bool segmentIntersect(const double xs1, const double ys1, const double xs2, const double ys2, double xc1, double yc1, double xc2, double yc2, double &xi, double &yi);
bool lineIntersect(const double xs1, const double ys1, const double xs2, const double ys2, double xc1, double yc1, double xc2, double yc2, double &xi, double &yi);
double lineIntersect2(const double xs1, const double ys1, const double xs2, const double ys2, double xc1, double yc1, double xc2, double yc2, double &xi, double &yi);
bool contourIntersect(int sl, const double* xs1, const double* ys1, const double* xs2, const double* ys2, double* xc, double* yc, int N, double* xsamp, double* ysamp);

void interp1bl(double *x, double* y, int nx, double*xi, int m, double * yi);
void VTareaPBS92(double* tubeLen, double* tubeSagit, double* area, double &sectLen);
void elecSpectrum(double * w, double* area, double* rZ_oral, double* iZ_oral, double tubeLen, double* rH, double* iH);
void VTarea2form(double *area, double tubeLen, double* rHeq, double* iHeq);
double cubicImax(double *x, int N);
double cubicImin(double *x, int N);
double** findFormantLocationWithDifferentiation(double *Hamp, double *freqs);
void findFormantLocationWithDifferentiationAndGuess(double *Hamp, double *freqs, double *Fguess, double *F);

void tonguePalateContact(double *s, double dt, double& dx, double& dxp, bool anticipation = false); // tongue contact with palate and pharynx
void tonguePalateContactDBG(double *s, double dt, double& dx, double& dxp); /
