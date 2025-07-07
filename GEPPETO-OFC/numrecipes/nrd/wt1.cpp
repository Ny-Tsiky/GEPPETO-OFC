void wt1(double a[], unsigned long n, int isign,

void (*wtstep)(double [], unsigned long, int))
{

unsigned long nn;


if (n < 4) return;

if (isign >= 0) {


for (nn=n;nn>=4;nn>>=1) (*wtstep)(a,nn,isign);

} else {


for (nn=4;nn<=n;nn<<=1) (*wtstep)(a,nn,isign);

}
}
