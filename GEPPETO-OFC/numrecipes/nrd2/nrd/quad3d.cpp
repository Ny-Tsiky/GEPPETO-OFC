static double xsav,ysav;
static double (*nrfunc)(double,double,double);

double quad3d(double (*func)(double, double, double), double x1, double x2)
{

double qgaus(double (*func)(double), double a, double b);

double f1(double x);


nrfunc=func;

return qgaus(f1,x1,x2);
}

double f1(double x)
{

double qgaus(double (*func)(double), double a, double b);

double f2(double y);

double yy1(double),yy2(double);


xsav=x;

return qgaus(f2,yy1(x),yy2(x));
}

double f2(double y)
{

double qgaus(double (*func)(double), double a, double b);

double f3(double z);

double z1(double,double),z2(double,double);


ysav=y;

return qgaus(f3,z1(xsav,y),z2(xsav,y));
}

double f3(double z)
{

return (*nrfunc)(xsav,ysav,z);
}
