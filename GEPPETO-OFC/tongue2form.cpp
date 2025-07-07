#define _USE_MATH_DEFINES
#include <math.h>
#include "tongue2form.h"
#include "mex.h"

#define sqr(x) ((x)*(x))



const int N_FORMANTS = 4;
extern const int NF;


// Perrier Boë Sock 1992
const double BETA = 1.5;
const double THRESH1 = 6; // mm
const double THRESH2 = 10; // mm
const double DZETAP[] = {1.2288, 8.1271, 8.1271, 8.1271, 4.5175, 4.5175, 4.5175, 4.5175, 4.5175, 4.5175, 4.6841, 4.6841, 4.6841, 4.6841, 4.5853, 4.5853, 4.5853, 4.5853, 4.5853, 5.0596, 5.0596, 5.0596, 5.0596, 5.0596, 5.0596, 5.0596, 3.9528, 3.9528, 3.9528}; 
const double DZETAG[] = {1.2288, 8.6646, 8.6646, 8.6646, 4.5175, 4.5175, 4.5175, 4.5175, 4.5175, 4.5175, 4.6841, 4.6841, 4.6841, 4.6841, 5.2178, 5.2178, 5.2178, 5.2178, 5.2178, 6.1823, 6.1823, 6.1823, 6.1823, 6.1823, 6.1823, 6.1823, 8.2219, 8.2219, 8.2219};
              
const double TONGUE_ANCHOR[] = {89, 50}; // basis of tongue
const double MIN_AREA = 1e-6; // 1 mm^2

// Badin-Fant formant estimation
const double CS = 351; // m/s                                       // sound velocity (at 35 deg, apparently) at 30 deg, 349 and density 1.16
const double RHO = 1.14; // kg/m^3                                  // air density
const double LAMBDA = 0.023; //J.s/(m K)  // 5.5e-5 (cal.s/(cm.K))  // coef of heat transfer for air (beware, not the usual in W/(m^2.K)
const double ETA = 1.4;                                             // adiabatic gas constant (dimensionless)
const double MU = 1.86e-5;  // Pa.s.  // was 1.86e-4 (dyn.s/cm^2)   // air viscosity
const double CP = 1.006e3; // J/(kg.K)  // 0.24; //cal/(g.K)        // specific heat of air (at cst pressure)
const double BP = 1.6e4; // N.s/m //1600; // dynes.s/cm^3           // component Rw of wall acoustic impedance (Badin Fant 1984 p.79), g/cm^2/s (ie, rayl)
const double MP = 14; // kg/m^2 // 1.4; // g/cm²                    // component Lw of wall acoustic impedance  (rayl.s)


// formant computation
const double FS = 16000;
const int N_FREQ = 500;
const double F_MAX = FS/2.2;
const double F_MIN = 150;

// vocal tract slicing
const double DSAGIT_LARYNX = 7.5; // 3.0 % mm (PP/RW value)
const double DSAGIT_LIP = 15; // mm (PP value)
const double LENGTH_LIP = 5; // mm (PP value)
const int NONLIP = 29; // index of last alveolar pt (from there, lip tube)
const int N_SLICES = NONLIP+2;

static bool vt2formlib_loaded = false;

const int N_VTcontInt1Pts = 8;
const int N_VTcontInt2Pts = 16;

const double xVTcontInt1[] = {86.0778,  85.7250,  86.7833,  89.2528,  94.1917,  98.4250,  99.4833,  100.5417};
const double yVTcontInt1[] = {19.9917,  20.3444,  23.1667,  25.9889,  30.5750,  35.5139,  40.1,  43.6278};
const double xVTcontInt2[] = {57.8556,  59.9722,  58.2083,  57.1500,  55.3861,  52.9167,  47.9778,  44.0972,  41.6278,  43.0389,  40.2167,  37.0417,  33.5139,  29.6333,  25.7528,  22.9306};
const double yVTcontInt2[] = {62.6778,  67.9694,  72.9083,  78.9056,  83.8444,  88.4306,  93.3694,  95.8389,  98.3083,  88.7833,  91.9583,  94.4278,  95.1333,  94.4278,  91.2528,  87.0194};
const double xVTptExt[] = {99.0704,  98.8440,  100.1323,  104.9305,  118.0030,  119.1289,  120.1487,  121.3280,  122.0038,  122.7036,  123.6620,  124.3822,  122.3089,  120.2356,  114.1112,  109.0201,  104.3535,  99.6934,  94.5411,  89.1054,  82.8561,  75.7324,  70.3503,  65.8133,  61.3508,  56.7259,  52.0281,  47.1782,  41.7896,  37.1,  32.1};
const double yVTptExt[] = {18.9586,  25.2843,  31.3265,  36.6655,  40.3465,  46.4013,  52.4973,  58.5614,  64.7263,  70.8664,  77.0748,  84.0064,  90.4387,  96.8710,  100.6924,  103.8689,  106.7806,  109.5044,  111.4091,  113.8172,  115.4532,  116.3369,  116.0667,  114.0150,  111.7838,  109.7031,  108.0386,  106.7402,  103.4951,  100.8090,  101.6825};
// const double xVTslices1[] = {79.40,  80.30,  81.30,  82.3,  83.30,  84.20,  85.20,  86.20,  87.20,  88.10,  89.10,  89.10,  89.10,  89.10,  89.10,  89.10,  89.10,  89.10,  89.10,  89.10,  89.10,  89.10,  84.50,  79.90,  75.30,  70.60,  66.0,  61.40,  55.65,  37.10,  32.10,  27.1};
// const double yVTslices1[] = {22.9,  29.0,  35.1,  41.2,  47.3,  53.4,  59.5,  65.6,  71.7,  77.8,  84.0,  84.0,  84.0,  84.0,  84.0,  84.0,  84.0,  84.0,  84.0,  84.0,  84.0,  84.0,  82.0,  80.1,  78.2,  76.3,  74.4,  72.5,  70.125,  70.125,  70.125,  70.125};
// VTslices1 slightly modified to allow more room for tongue mvt:
const double xVTslices1[] = {76.7050,  77.6050,  78.6050,  79.6050,  80.6050,  81.5050,  82.5050,  83.5050,  84.5050,  85.4050,  85.7000,  85.7000,  85.7000,  85.7000,  85.7000,  85.7000,  85.7000,  85.7000,  85.7000,  85.7000,  85.7000,  85.7000,  85.7000,  80.9550,  76.3550,  71.6550,  67.0550,  62.4550,  56.7050,  37.1000,  32.1000,  27.1000};
const double yVTslices1[] = {23.4400,  29.5400,  35.6400,  41.7400,  47.8400,  53.9400,  60.0400,  66.1400,  72.2400,  78.3400,  79.0000,  79.0000,  79.0000,  79.0000,  79.0000,  79.0000,  79.0000,  79.0000,  79.0000,  79.0000,  79.0000,  79.0000,  79.0000,  77.5600,  75.6600,  73.7600,  71.8600,  69.9600,  67.5850,  67.5850,  67.5850,  67.5850};
const double xVTslices2[] = {133.3,  134.2,  135.2,  136.2,  137.2,  138.1,  139.1,  140.1,  141.1,  142.0,  143.0,  144.1,  143.1,  139.9,  134.8,  128.0,  119.7,  110.2,  99.8,  89.11,  78.4,  68.1,  63.4,  58.8,  54.2,  49.5,  44.9,  40.3,  34.55,  37.1,  32.1,  27.1};
const double yVTslices2[] = {12.1,  18.2,  24.3,  30.4,  36.5,  42.6,  48.7,  54.8,  60.9,  67.0,  73.2,  84.01,  94.7,  105.0,  114.5,  122.8,  129.7,  134.8,  137.9,  138.9,  137.9,  134.8,  132.8,  130.9,  129.0,  127.1,  125.2,  123.3,  120.925,  120.925,  120.925,  120.925};









#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
const int M = 7;
const int NSTACK = 50;

// NR's sort function adapted to standard C indexing
void heapsort(int *x, int N)
{
  int i, ir = N-1, j, k, l = 0;
  int jstack = -1, istack[NSTACK];
  int a,temp;


  for (;;) {
    if (ir-l < M) {
      for (j=l+1; j<=ir; j++) {

        a = x[j];
        for (i = j-1; i >= 0; i--) {
          if (x[i] <= a) break;
          x[i+1] = x[i];
        }
        x[i+1] = a;
      }
      if (jstack < 0) break;
      
      ir = istack[jstack--];
      l = istack[jstack--];
    } 
    else {
      k = (l+ir) >> 1;
      SWAP(x[k],x[l+1])
      if (x[l+1] > x[ir]) {
        SWAP(x[l+1],x[ir])
      }
      if (x[l] > x[ir]) {
        SWAP(x[l],x[ir])
      }
      if (x[l+1] > x[l]) {
        SWAP(x[l+1],x[l])
      }
      i = l+1;
      j = ir;
      a = x[l];
      
      for (;;) {
        do i++; while (x[i] < a);
        do j--; while (x[j] > a);
        if (j < i) break;
        SWAP(x[i],x[j]);
      }
      
      x[l] = x[j];
      x[j] = a;
      jstack += 2;
      if (jstack > NSTACK) {mexErrMsgTxt("NSTACK too small in sort."); return;} //fprintf(stderr,"NSTACK too small in sort.");
      if (ir-i+1 >= j-l) {
        istack[jstack] = ir;
        istack[jstack-1]=i;
        ir = j-1;
      } 
      else {
        istack[jstack] = j-1;
        istack[jstack-1]=l;
        l = i;
      }
    }
  }
}

// NR's sort function adapted to standard C indexing
void heapsort(double *x, int N)
{
  int i, ir = N-1, j, k, l = 0;
  int jstack = -1, istack[NSTACK];
  double a,temp;


  for (;;) {
    if (ir-l < M) {
      for (j=l+1; j<=ir; j++) {

        a = x[j];
        for (i = j-1; i >= 0; i--) {
          if (x[i] <= a) break;
          x[i+1] = x[i];
        }
        x[i+1] = a;
      }
      if (jstack < 0) break;
      
      ir = istack[jstack--];
      l = istack[jstack--];
    } 
    else {
      k = (l+ir) >> 1;
      SWAP(x[k],x[l+1])
      if (x[l+1] > x[ir]) {
        SWAP(x[l+1],x[ir])
      }
      if (x[l] > x[ir]) {
        SWAP(x[l],x[ir])
      }
      if (x[l+1] > x[l]) {
        SWAP(x[l+1],x[l])
      }
      i = l+1;
      j = ir;
      a = x[l];
      
      for (;;) {
        do i++; while (x[i] < a);
        do j--; while (x[j] > a);
        if (j < i) break;
        SWAP(x[i],x[j]);
      }
      
      x[l] = x[j];
      x[j] = a;
      jstack += 2;
      if (jstack > NSTACK) {mexErrMsgTxt("NSTACK too small in sort."); return;} 
      if (ir-i+1 >= j-l) {
        istack[jstack] = ir;
        istack[jstack-1]=i;
        ir = j-1;
      } 
      else {
        istack[jstack] = j-1;
        istack[jstack-1]=l;
        l = i;
      }
    }
  }
}


// S1 + alpha.S1S2 = C1 + beta.C1C2, alpha and beta both in [0 1] (that's the idea, the sign of alpha and beta may be difft from that of alph and bet)
bool segmentIntersect(const double xs1, const double ys1, const double xs2, const double ys2, double xc1, double yc1, double xc2, double yc2, double &xi, double &yi)
{
  int i, j, k;
  double xDeltaSlice, yDeltaSlice, xDeltaContour, yDeltaContour;
  double xD1, xD2, yD1, yD2, alph, bet, den;
  bool ok = false;
  
    xDeltaSlice = xs2 - xs1;
    yDeltaSlice = ys2 - ys1;
    xDeltaContour = xc2 - xc1;
    yDeltaContour = yc2 - yc1;
    xD1 = xc1 - xs1;
    yD1 = yc1 - ys1;

    if (fabs(xDeltaSlice) > EPS) {
    den = yDeltaContour-xDeltaContour/xDeltaSlice*yDeltaSlice;
    if (fabs(den) > EPS) {
      bet = (-yD1 + yDeltaSlice/xDeltaSlice*xD1)/den;
      alph = (xD1 + bet*xDeltaContour)/xDeltaSlice;
      
    }
    else return(ok);      
    }
    else if (fabs(xDeltaContour) > EPS)
    {
        bet = -xD1/xDeltaContour;
        alph = (yD1 + bet*yDeltaContour)/yDeltaSlice;
    }
    else return(ok);

  if (alph >= 0 && alph <= 1 && bet >= 0 && bet <= 1) {
    xi = xs1 + alph*xDeltaSlice;
    yi = ys1 + alph*yDeltaSlice;
    ok = true;
  }
  return(ok);
}  

// version with infinite intersecting grid line; not entirely correct (can overestimate vocal tract size in the region where grid lines intersect) but robust to any tongue motion
// S1 + alpha.S1S2 = C1 + beta.C1C2
double lineIntersect2(const double xs1, const double ys1, const double xs2, const double ys2, double xc1, double yc1, double xc2, double yc2, double &xi, double &yi)
{
  int i, j, k;
  double xDeltaSlice, yDeltaSlice, xDeltaContour, yDeltaContour;
  double xD1, xD2, yD1, yD2, alph, retalph = -1, bet, den;

  
  xDeltaSlice = xs2 - xs1;
  yDeltaSlice = ys2 - ys1;
  xDeltaContour = xc2 - xc1;
  yDeltaContour = yc2 - yc1;
  xD1 = xc1 - xs1;
  yD1 = yc1 - ys1;

  if (fabs(xDeltaSlice) > EPS) {
    den = yDeltaContour-xDeltaContour/xDeltaSlice*yDeltaSlice;
    if (fabs(den) > EPS) {
      bet = (-yD1 + yDeltaSlice/xDeltaSlice*xD1)/den;
      alph = (xD1 + bet*xDeltaContour)/xDeltaSlice;
    }
    else return(retalph);
  }
  else if (fabs(xDeltaContour) > EPS)
  {
    bet = -xD1/xDeltaContour;
    alph = (yD1 + bet*yDeltaContour)/yDeltaSlice;
  }
  else return(retalph);

  if (bet >= 0 && bet <= 1) {
    xi = xs1 + alph*xDeltaSlice;
    yi = ys1 + alph*yDeltaSlice;
    retalph = alph;
  }
  return(retalph);
}  

bool contourIntersect(int sl, const double* xs1, const double* ys1, const double* xs2, const double* ys2, double* xc, double* yc, int N, double* xsamp, double* ysamp)
{
  int i;
  double xi, yi; 
  // bool ok;
  double relpos, minrelpos = 0;//
    
  for (i = 0; i < N-1 ; i++) {
    relpos = lineIntersect2(xs1[sl],ys1[sl],xs2[sl],ys2[sl],xc[i],yc[i],xc[i+1],yc[i+1],xi,yi);
    if (relpos > minrelpos) {
      xsamp[sl] = xi; 
      ysamp[sl] = yi;
      minrelpos = relpos;
    }
  } 
  return(minrelpos > 0);
}


double contourIntersectDBG(int sl, const double* xs1, const double* ys1, const double* xs2, const double* ys2, double* xc, double* yc, int N, double* xsamp, double* ysamp)
{
  int i;
  double xi, yi; 
  // bool ok;
  double relpos, minrelpos = 0;
    
  for (i = 0; i < N-1 ; i++) {
    relpos = lineIntersect2(xs1[sl],ys1[sl],xs2[sl],ys2[sl],xc[i],yc[i],xc[i+1],yc[i+1],xi,yi);
    mexPrintf("relpos (i=%d) = %g;  seg : (%g %g)-(%g %g)\n",i,relpos,xc[i],yc[i],xc[i+1],yc[i+1]);
    if (relpos > minrelpos) { // choose intersecting point that is furthest from slice base (point s1), that is, closest to exterior contour
      xsamp[sl] = xi; 
      ysamp[sl] = yi;
      minrelpos = relpos;
    }
  } 
  return(minrelpos);
}

// from a MEX by B Luong - called by my interp1bl Matlab function
// x must be sorted in ascending order; x and y have the same length
void interp1bl(double *x, double* y, int nx, double*xi, int m, double * yi)
{
  int k, i1, i9, imid;
  double xik;
  
  // Loop over the points
  for (k=m; k--;) // Reverse of for (k=0; k<m; k++) {...}
  {
      // Get data value
      xik = xi[k];
      
      i1=0;
      i9=nx-1;
      while (i9>i1+1) // Dichotomy search
      {
          imid = (i1+i9+1)/2;
          if (x[imid]<xik) i1=imid;
          else i9=imid;
      } // of while loop
      if (i1==i9)
          yi[k] = y[i1];
      else
          yi[k] = y[i1] + (y[i9]-y[i1])*(xik-x[i1])/(x[i9]-x[i1]);
  } 
  
}



void VTareaPBS92(double* tubeLen, double* tubeSagit, double* area, double &sectLen)
{
// % area: matrix of vocal tract tube areas (Ntimesteps x NBSECT_FINAL)
// % sectlen: vector of vocal tract tube lengths (Ntimesteps x 1)
// % WARNING: input values are in mm or mm^2 and return values are in SI
  
  int i;
  double dzeta;
  double areaInt[N_SLICES], VTpos[N_SLICES+1], VTsampPos[N_VTSECTIONS+1], intVol[N_SLICES+1], sampVol[N_VTSECTIONS+1], ar;
  
  for (i = 0; i < NONLIP; i++) { 
    if ((tubeLen[i] <= THRESH1) || ((i < NONLIP-1) && (tubeLen[i+1] <= THRESH1)) || ((i > 0) && (tubeLen[i-1] <= THRESH1))) //Pmask
      dzeta = DZETAP[i];
    else if (tubeLen[i] >= THRESH2)
      dzeta = DZETAG[i];
    else
      dzeta = DZETAP[i] + (tubeSagit[i]-THRESH1) * (DZETAG[i]-DZETAP[i]) / (THRESH2-THRESH1);  // this could be C-inf with gfun
    areaInt[i] = dzeta*pow(tubeSagit[i],BETA);
  }

  for (i = NONLIP; i < N_SLICES; i++) areaInt[i] = M_PI*sqr(tubeSagit[i])/4;
  
  VTpos[0] = 0;
  intVol[0] = 0;
  for (i = 1; i < N_SLICES+1; i++) {
    VTpos[i] = tubeLen[i-1] + VTpos[i-1];
    intVol[i] = tubeLen[i-1]*areaInt[i-1] + intVol[i-1]; // integrated volume for interpolation
  }
  sectLen = VTpos[N_SLICES]/N_VTSECTIONS;
  
  
  // numerical integration of volume, resampling, and differentiation  // this is a bit more tricky to C-infinitize
  for ( i = 0; i <= N_VTSECTIONS; i++) VTsampPos[i] = sectLen*i;
  interp1bl(VTpos, intVol, N_SLICES+1, VTsampPos, N_VTSECTIONS+1, sampVol);
  for ( i = 0; i < N_VTSECTIONS; i++) {
    ar = (sampVol[i+1] - sampVol[i])/sectLen/1e6; // in m^2
    area[i] = (ar > MIN_AREA) ? ar : MIN_AREA; // safety
  }
  sectLen = sectLen/1e3; // in m 
  
}




// *********** //
// already defined in tonguefuns as softpospart but this code may be used in other context so must be self-contained
const double SMOOTH_EPS = 1e-1;

double softpositivepart(double x, double zeroshift=SMOOTH_EPS)
{
  double y = (sqr(x)+fabs(x)+zeroshift*2)/(fabs(x)+zeroshift*2)-1+zeroshift*2; // softabs(x,zeroshift*2)
  return((x + y)/2);
}


// *********** //





// special version for NR-style vector tong of length 2*N_NODES
void tongue2VTarea(double *tong, double *VTarea, double &sectLen)
{
  int i, j, k, sl;
  bool ok = true;
  double *xVTcontInt, *yVTcontInt;
  double xVTptInt[N_SLICES], yVTptInt[N_SLICES]; 
  double tubeLen[N_SLICES], tubeSagittal[N_SLICES];
  double S;
  int N_VTcontIntPts;
  double newcoord;

  N_VTcontIntPts = N_VTcontInt1Pts+N_NODES+1+N_VTcontInt2Pts;  
  xVTcontInt = new double[N_VTcontIntPts];
  yVTcontInt = new double[N_VTcontIntPts];
  
  k = 0;
  for ( i = 0; i < N_VTcontInt1Pts; i++, k++) {
    xVTcontInt[k] = xVTcontInt1[i];
    yVTcontInt[k] = yVTcontInt1[i];
  }
  xVTcontInt[k] = TONGUE_ANCHOR[0];
  yVTcontInt[k++] = TONGUE_ANCHOR[1];
  for ( i = 1; i <= N_NODES; i++, k++) { 
    xVTcontInt[k] = tong[i];
    yVTcontInt[k] = tong[i+N_NODES];
  }

 // find overlap between internal contour and tongue contour
  for (j = 0; (j < N_NODES) && (xVTcontInt2[j] > tong[N_NODES]); j++); 
  for (i = 0; i < j; i++, k++) {
    xVTcontInt[k] = xVTcontInt2[j]; // remove part of contInt2 that is overlapping tongue
    yVTcontInt[k] = yVTcontInt2[j]; 
  }
  for ( i = j; i < N_VTcontInt2Pts; i++, k++) {
    xVTcontInt[k] = xVTcontInt2[i];
    yVTcontInt[k] = yVTcontInt2[i];
  }
  
  
  // all the following assumes that slice numbering goes in same direction as contourpoint numbering!
  

  for (sl = 0; sl < N_SLICES && ok; sl++) {
      ok = contourIntersect(sl, xVTslices1, yVTslices1, xVTslices2, yVTslices2, xVTcontInt, yVTcontInt, N_VTcontIntPts, xVTptInt, yVTptInt);//, sl == 0);
    
      if (!ok) {
        double relp = contourIntersectDBG(sl, xVTslices1, yVTslices1, xVTslices2, yVTslices2, xVTcontInt, yVTcontInt, N_VTcontIntPts, xVTptInt, yVTptInt);// DBG
        char msg[100]; sprintf(msg,"Vocal tract slicing incorrect (slice %d: (%g, %g)-(%g, %g); relp=%g)",sl,xVTslices1[sl],yVTslices1[sl],xVTslices2[sl],yVTslices2[sl],relp); mexErrMsgTxt(msg); return;
      };
    }

  for ( i = 13; i < N_SLICES; i++) {
    newcoord = yVTptExt[i] - softpositivepart(yVTptExt[i] - yVTptInt[i]); 
    xVTptInt[i] = xVTptInt[i] + (newcoord - yVTptInt[i])/(yVTptExt[i] - yVTptInt[i])*(xVTptExt[i]-xVTptInt[i]); // move point along slice
    yVTptInt[i] = newcoord;
  }
  for ( i = 0; i < 13; i++) {
    newcoord = xVTptExt[i] - softpositivepart(xVTptExt[i] - xVTptInt[i]); // smooth constraint
    yVTptInt[i] = yVTptInt[i] + (newcoord - xVTptInt[i])/(xVTptExt[i] - xVTptInt[i])*(yVTptExt[i]-yVTptInt[i]); // move point along slice
    xVTptInt[i] = newcoord;
  }

  // % Fix lip aperture (needed for last but one tube...)
  yVTptInt[N_SLICES-2] = yVTptExt[N_SLICES-2] - DSAGIT_LIP;
  yVTptInt[N_SLICES-1] = yVTptExt[N_SLICES-1] - DSAGIT_LIP;

    for ( i = 0; i < N_SLICES-1; i++) {
    S = fabs((xVTptInt[i]-xVTptExt[i+1])*(yVTptExt[i]-yVTptInt[i+1]) - (xVTptExt[i]-xVTptInt[i+1])*(yVTptInt[i]-yVTptExt[i+1])) / 2.0;
    tubeLen[i] = sqrt(sqr(xVTptInt[i]+xVTptExt[i]-xVTptInt[i+1]-xVTptExt[i+1]) + sqr(yVTptInt[i]+yVTptExt[i]-yVTptInt[i+1]-yVTptExt[i+1])) / 2.0; //  equivalent tube
    tubeSagittal[i] = S/tubeLen[i]; // equivalent tube
  }

  tubeSagittal[N_SLICES-1] = DSAGIT_LIP; // on impose une dist. sagittale aux lèvres
  tubeLen[N_SLICES-1] = LENGTH_LIP; // idem pour la longueur

  // % fix laryngeal sagittal distance, RW 09/2012
  for (i = 0; i< 3; i++) tubeSagittal[i] = DSAGIT_LARYNX; // Rajouté en avril 04 % corrects formants, from a Fant paper on actual VT estimation

  VTareaPBS92(tubeLen, tubeSagittal, VTarea, sectLen);
  
  delete[] xVTcontInt;
  delete[] yVTcontInt;
} 

//**********************************************************************************************************************************************//
//**********************************************************************************************************************************************//

void elecSpectrum(double * w, double* area, double* rZ_oral, double* iZ_oral, double tubeLen, double* rH, double* iH)
{
  int i, j;
  double S[N_VTSECTIONS], L[N_VTSECTIONS], C[N_VTSECTIONS];
  double R_coef[N_FREQ], G_coef[N_FREQ], YP_coef[N_FREQ];
  double R, G, rYP, iYP, rZ, iZ, rY, iY;
  double ra, rb, rc, rd;
  double ia, ib, ic, id;
  double rA, rB, rC, rD, nrA, nrB, nrC, nrD;
  double iA, iB, iC, iD, niA, niB, niC, niD;
  double sqmodulus;
  double rHi, iHi;
  
  
  for (j = 0; j < N_FREQ; j++) {
    R_coef[j] = sqrt(RHO*MU/2*w[j]) ;
    G_coef[j] = (ETA-1)/(RHO*sqr(CS))*sqrt(LAMBDA*w[j]/(2*CP*RHO));
    YP_coef[j] = 1./(sqr(BP)+sqr(MP)*sqr(w[j]));
  }
  
 
  for ( i = 0; i<N_VTSECTIONS; i++) {
    S[i] = 2*sqrt(area[i]*M_PI) ;
    L[i] = RHO/area[i]*tubeLen;
    C[i] = area[i]*tubeLen/RHO/sqr(CS);
  }
  
  for (j = 0; j < N_FREQ; j++) {
    for ( i = 0; i<N_VTSECTIONS; i++) {

      R = S[i]*tubeLen/sqr(area[i]) * R_coef[j] ;  
      G = S[i]*tubeLen * G_coef[j] ;

      rYP = S[i]*tubeLen * ( BP*YP_coef[j] );
      iYP = S[i]*tubeLen * ( -MP*w[j]*YP_coef[j] );

      rZ = R;
      iZ = L[i]*w[j];
      rY = G + rYP;
      iY = C[i]*w[j] + iYP;
            
      ra = 1 + (rZ*rY - iZ*iY)/2;
      ia = (iZ*rY + rZ*iY)/2;
      rb = - (rZ + ((sqr(rZ) - sqr(iZ))*rY - 2*rZ*iZ*iY)/4);
      ib = - (iZ + (2*rZ*iZ*rY + (sqr(rZ) - sqr(iZ))*iY)/4);
      rc = - rY;
      ic = - iY;
      rd = ra;
      id = ia;


      if (i == 0) {
        rA = ra; iA = ia;
        rB = rb; iB = ib;
        rC = rc; iC = ic;
        rD = rd; iD = id;        
      }
      else {
        nrA = ra*rA + rb*rC - (ia*iA + ib*iC);
        niA = ra*iA + rb*iC + (ia*rA + ib*rC);
        nrB = ra*rB + rb*rD - (ia*iB + ib*iD);
        niB = ra*iB + rb*iD + (ia*rB + ib*rD);
        nrC = rc*rA + rd*rC - (ic*iA + id*iC);
        niC = rc*iA + rd*iC + (ic*rA + id*rC);
        nrD = rc*rB + rd*rD - (ic*iB + id*iD);
        niD = rc*iB + rd*iD + (ic*rB + id*rD);
        rA = nrA; iA = niA;
        rB = nrB; iB = niB;
        rC = nrC; iC = niC;
        rD = nrD; iD = niD;
      }
    } 


  rHi = rA - (rC*rZ_oral[j] - iC*iZ_oral[j]);
  iHi = iA - (iC*rZ_oral[j] + rC*iZ_oral[j]);
  sqmodulus = sqr(rHi)+sqr(iHi);
  rH[j] =   rHi/sqmodulus;
  iH[j] = - iHi/sqmodulus;
  
}
  
}





void   VTarea2form(double *area, double tubeLen, double* rHeq, double* iHeq)
{
  int i;
  double w[N_FREQ], rZ_oral[N_FREQ], iZ_oral[N_FREQ];

  for ( i = 0; i < N_FREQ; i++) {
    w[i] = 2*M_PI*(F_MIN + i*(F_MAX-F_MIN)/(N_FREQ-1));
    rZ_oral[i] = RHO/(2*M_PI*CS)*sqr(w[i]);
    iZ_oral[i] = 8*RHO/(3*M_PI*sqrt(M_PI*area[N_VTSECTIONS-1]))*w[i];
  }
  
  elecSpectrum(w, area, rZ_oral, iZ_oral, tubeLen, rHeq, iHeq);  
}



double cubicImax(double *x, int N)
{
  int i;
  double *y, *I, den;
  double ymax, maxI;
  
  y = new double[N];
  I = new double[N];
  
  y[0] = x[0];
  y[N-1] = x[N-1];
  I[0] = I[N-1] = 0;

  // find extrema of cubic fits (could be pos of min or max, later max() will eliminate irrelevant I values)
  for ( i = 1; i < N-1; i++) {
    den = x[i-1]+x[i+1]-2*x[i];
    if (fabs(den) > EPS) {
      y[i] = -sqr(x[i+1]-x[i-1])/8./den + x[i];  // extremal value
      I[i] = (x[i-1]-x[i+1])/(2*den); // location of extremum, between -1 & 1, 0 corresponds to index value i
    }
    else {
      y[i] = x[i];
      I[i] = 0;
    }
    if ( fabs(I[i]) > 1 ) y[i] = x[i];
  }

  // max()
  ymax = y[0]; maxI = 0;
  for ( i = 1; i < N-1; i++) {
    if (y[i] > ymax) {ymax = y[i]; maxI = I[i]+i;}
  }
  
  delete[] y;
  delete[] I;
  return(maxI);
}

double cubicImin(double *x, int N)
{
  int i;
  double *y, *I, den;
  double ymin, minI;
  
  y = new double[N];
  I = new double[N];
  
  y[0] = x[0];
  y[N-1] = x[N-1];
  I[0] = I[N-1] = 0;
  
  // find extrema of cubic fits (could be pos of min or max, later min() will eliminate irrelevant I values)
  for ( i = 1; i < N-1; i++) {
    den = x[i-1]+x[i+1]-2*x[i];
    if (fabs(den) > EPS) {
      y[i] = -sqr(x[i+1]-x[i-1])/8./den + x[i]; // extremal value
      I[i] = (x[i-1]-x[i+1])/(2*den); // location of extremum, between -1 & 1, 0 corresponds to index value i
    }
    else {
      y[i] = x[i];
      I[i] = 0;
    }
    if ( fabs(I[i]) > 1 ) y[i] = x[i];
  }

  // min()
  ymin = y[0]; minI = 0;
  for ( i = 1; i < N-1; i++) {
    if (y[i] < ymin) {ymin = y[i]; minI = I[i]+i;}
  }

  delete[] y;
  delete[] I;
  return(minI);
}



void findFormantLocationWithDifferentiation(double *Hamp, double *freqs, double *F)
{
  const int N_MAX_PEAKS = 20; 
  const double SIGMOID_THRESHOLD = 1e-3;
  double freqPerBin = freqs[1]-freqs[0];
  double zerocrossThreshold = 0.6/freqPerBin;
  double freqTol = 50/freqPerBin; // 50 Hz
  int i, j, k, ka = 0, kb = 0, kc = 0;
  static double dHamp[N_FREQ], ddHamp[N_FREQ]; // static to save some allocation time
  static int pklocs1a[N_MAX_PEAKS];
  static int pklocs1b[N_MAX_PEAKS];
  static int pklocs1c[N_MAX_PEAKS];
  static double finepklocs1a[N_MAX_PEAKS];
  static double finepklocs1b[N_MAX_PEAKS];
  static double finepklocs1c[N_MAX_PEAKS];
  static double pklocs1[N_MAX_PEAKS];
  int Nform;
  double findx;
  
  
  dHamp[0] = 0;
  ddHamp[0] = 0;
  for ( i = 1; i < N_FREQ; i++) dHamp[i] = Hamp[i] - Hamp[i-1];
  for ( i = 1; i < N_FREQ; i++) ddHamp[i] = dHamp[i] - dHamp[i-1];
  
  for ( i = 0; i < N_FREQ; i++) {    
    if ( (i < N_FREQ-4) && (dHamp[i+2]*dHamp[i+3] < 0) && (dHamp[i+2] > 0) ) {pklocs1a[ka++] = i+2;}

    if ( (i < N_FREQ-4) && (ddHamp[i+2]*ddHamp[i+3] < 0) && (ddHamp[i+2] < 0) && (dHamp[i+2] < zerocrossThreshold) && (dHamp[i+2] > 0)) {pklocs1b[kb++] = i+2;}

    if ( (i < N_FREQ-4) && (ddHamp[i+2]*ddHamp[i+3] < 0) && (ddHamp[i+2] > 0) && (dHamp[i+2] > -zerocrossThreshold) && (dHamp[i+2] < 0)) {pklocs1c[kc++] = i+2;}
  }
  

  for ( k = 0; k < ka; k++) {
    findx = cubicImax(Hamp+pklocs1a[k]-2, 5);
    finepklocs1a[k] = findx + pklocs1a[k]-2;
  }

  for ( k = 0; k < kb; k++) {
    findx = cubicImin(dHamp+pklocs1b[k]-2, 5);
    finepklocs1b[k] = findx + pklocs1b[k]-2;
    for (j = 0; j< ka; j++) if (fabs(finepklocs1b[k] - finepklocs1a[j]) < freqTol) finepklocs1b[k] = 0; // tag for deletion
  }
  
  for ( k = 0; k < kc; k++) {
    findx = cubicImax(dHamp+pklocs1c[k]-2, 5);
    finepklocs1c[k] = findx + pklocs1c[k]-2;
    for (j = 0; j< ka; j++) if (fabs(finepklocs1c[k] - finepklocs1a[j]) < freqTol) finepklocs1c[k] = 0; // tag for deletion
  }

  j = 0;
  for (k = 0; k < ka; k++) pklocs1[j++] = finepklocs1a[k];
  for (k = 0; k < kb; k++) if (finepklocs1b[k] > 0) pklocs1[j++] = finepklocs1b[k];
  for (k = 0; k < kc; k++) if (finepklocs1c[k] > 0) pklocs1[j++] = finepklocs1c[k];
  Nform = j;  
  heapsort(pklocs1, Nform);
  
  
  for (j = 0; j < NF; j++) 
    if (j < Nform) F[j+1] = pklocs1[j]*freqPerBin + freqs[0]; // if() clause for safety but this should never happen in practice
}


 void VTarea2formantLocation(double *vtArea, double vtTubeLength, double *F)
{
  // % Compute location of formants peaks from an area function according to Badin's program
  // % vt_aire2frm modified by Ralf and PP in 2010 to solve situations in which
  // % the original program did not find any formants. This enhanced version
  // % looks for shoulders in the spectral envelope by computing the derivative
  // % of the envelope.

  int i;
  double rHeq[N_FREQ], iHeq[N_FREQ], Hamp[N_FREQ];
  double freqs[N_FREQ];
    
  VTarea2form(vtArea, vtTubeLength, rHeq, iHeq); 

  for ( i = 0; i < N_FREQ; i++ ) {
    Hamp[i] = 20*log10(sqrt(sqr(rHeq[i])+sqr(iHeq[i]))); // complex modulus
    freqs[i] = F_MIN + i*(F_MAX-F_MIN)/(N_FREQ-1);
  }
  
  // % procedure to detect poles based on differentiation
  findFormantLocationWithDifferentiation(Hamp, freqs, F);
} 
