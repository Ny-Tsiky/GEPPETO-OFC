void ftest(double data1[], unsigned long n1, double data2[], unsigned long n2,

double *f, double *prob)
{

void avevar(double data[], unsigned long n, double *ave, double *var);

double betai(double a, double b, double x);

double var1,var2,ave1,ave2,df1,df2;


avevar(data1,n1,&ave1,&var1);

avevar(data2,n2,&ave2,&var2);

if (var1 > var2) {


*f=var1/var2;


df1=n1-1;


df2=n2-1;

} else {


*f=var2/var1;


df1=n2-1;


df2=n1-1;

}

*prob = 2.0*betai(0.5*df2,0.5*df1,df2/(df2+df1*(*f)));

if (*prob > 1.0) *prob=2.0-*prob;
}
