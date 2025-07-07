double fredin(double x, int n, double a, double b, double t[], double f[],

double w[], double (*g)(double), double (*ak)(double, double))
{

int i;

double sum=0.0;


for (i=1;i<=n;i++) sum += (*ak)(x,t[i])*w[i]*f[i];

return (*g)(x)+sum;
}
