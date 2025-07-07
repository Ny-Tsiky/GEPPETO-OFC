void qrsolv(double **a, int n, double c[], double d[], double b[])
{

void rsolv(double **a, int n, double d[], double b[]);

int i,j;

double sum,tau;


for (j=1;j<n;j++) {


for (sum=0.0,i=j;i<=n;i++) sum += a[i][j]*b[i];


tau=sum/c[j];


for (i=j;i<=n;i++) b[i] -= tau*a[i][j];

}

rsolv(a,n,d,b);
}
