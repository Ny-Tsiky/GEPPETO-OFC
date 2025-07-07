#define NRANSI
#include "nrutil.h"

void correl(double data1[], double data2[], unsigned long n, double ans[])
{

void realft(double data[], unsigned long n, int isign);

void twofft(double data1[], double data2[], double fft1[], double fft2[],


unsigned long n);

unsigned long no2,i;

double dum,*fft;


fft=dvector(1,n<<1);

twofft(data1,data2,fft,ans,n);

no2=n>>1;

for (i=2;i<=n+2;i+=2) {


ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/no2;


ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/no2;

}

ans[2]=ans[n+1];

realft(ans,n,-1);

free_dvector(fft,1,n<<1);
}
#undef NRANSI
