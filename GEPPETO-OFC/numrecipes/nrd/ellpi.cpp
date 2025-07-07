#include <math.h>
#define NRANSI
#include "nrutil.h"

double ellpi(double phi, double en, double ak)
{

double rf(double x, double y, double z);

double rj(double x, double y, double z, double p);

double cc,enss,q,s;


s=sin(phi);

enss=en*s*s;

cc=DSQR(cos(phi));

q=(1.0-s*ak)*(1.0+s*ak);

return s*(rf(cc,q,1.0)-enss*rj(cc,q,1.0,1.0+enss)/3.0);
}
#undef NRANSI
