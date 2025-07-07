void avevar(double data[], unsigned long n, double *ave, double *var)
{

unsigned long j;

double s,ep;


for (*ave=0.0,j=1;j<=n;j++) *ave += data[j];

*ave /= n;

*var=ep=0.0;

for (j=1;j<=n;j++) {


s=data[j]-(*ave);


ep += s;


*var += s*s;

}

*var=(*var-ep*ep/n)/(n-1);
}
