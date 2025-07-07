void fpoly(double x, double p[], int np)
{

int j;


p[1]=1.0;

for (j=2;j<=np;j++) p[j]=p[j-1]*x;
}
