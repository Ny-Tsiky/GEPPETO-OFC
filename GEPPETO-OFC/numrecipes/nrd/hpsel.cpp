void hpsel(unsigned long m, unsigned long n, double arr[], double heap[])
{

void sort(unsigned long n, double arr[]);

void nrerror(char error_text[]);

unsigned long i,j,k;

double swap;


if (m > n/2 || m < 1) nrerror("probable misuse of hpsel");

for (i=1;i<=m;i++) heap[i]=arr[i];

sort(m,heap);

for (i=m+1;i<=n;i++) {


if (arr[i] > heap[1]) {



heap[1]=arr[i];



for (j=1;;) {




k=j << 1;




if (k > m) break;




if (k != m && heap[k] > heap[k+1]) k++;




if (heap[j] <= heap[k]) break;




swap=heap[k];




heap[k]=heap[j];




heap[j]=swap;




j=k;



}


}

}
}
