void piksr2(int n, double arr[], double brr[])
{

int i,j;

double a,b;


for (j=2;j<=n;j++) {


a=arr[j];


b=brr[j];


i=j-1;


while (i > 0 && arr[i] > a) {



arr[i+1]=arr[i];



brr[i+1]=brr[i];



i--;


}


arr[i+1]=a;


brr[i+1]=b;

}
}
