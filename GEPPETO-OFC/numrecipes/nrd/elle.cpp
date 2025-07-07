#include <math.h>
#define NRANSI
#include "nrutil.h"

double elle(double phi, double ak)
{

double rd(double x, double y, double z);

double rf(double x, double y, double z);

double cc,q,s;


s=sin(phi);

cc=DSQR(cos(phi));

q=(1.0-s*ak)*(1.0+s*ak);

return s*(rf(cc,q,1.0)-(DSQR(s*ak))*rd(cc,q,1.0)/3.0);
}
#undef NRANSI
