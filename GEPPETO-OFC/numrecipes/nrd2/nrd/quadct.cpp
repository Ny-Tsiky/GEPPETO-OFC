void quadct(double x, double y, double xx[], double yy[], unsigned long nn,

double *fa, double *fb, double *fc, double *fd)
{

unsigned long k,na,nb,nc,nd;

double ff;

na=nb=nc=nd=0;

for (k=1;k<=nn;k++) {


if (yy[k] > y) {



xx[k] > x ? ++na : ++nb;


} else {



xx[k] > x ? ++nd : ++nc;


}

}

ff=1.0/nn;

*fa=ff*na;

*fb=ff*nb;

*fc=ff*nc;

*fd=ff*nd;
}
