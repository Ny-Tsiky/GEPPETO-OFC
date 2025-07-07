void matadd(double **a, double **b, double **c, int n)
{

int i,j;


for (j=1;j<=n;j++)


for (i=1;i<=n;i++)



c[i][j]=a[i][j]+b[i][j];
}



void matadd(double **a, double **b, double **c, int m, int n)
{

int i,j;


for (j=1;j<=n;j++)


for (i=1;i<=m;i++)



c[i][j]=a[i][j]+b[i][j];
}
