#include<stdlib.h>
#include<stdio.h>
#include<math.h>

int tridiag(double diag[], double a[], double c[], double r[], double x[], int n)
{
        int j;
        double bet,*gam;

        gam = malloc(n*sizeof(double));
        if (fabs(diag[0]) == 0.0) return -1;    // Error
        x[0]=r[0]/(bet=diag[0]);
        for (j=1;j<n;j++) {
             gam[j]=c[j-1]/bet;
             bet=diag[j]-a[j]*gam[j];
             if (fabs(bet) == 0.0) return -1;   // Error
             x[j]=(r[j]-a[j]*x[j-1])/bet;
        }
        for (j=n-2;j>=0;j--)
             x[j] -= gam[j+1]*x[j+1];
        free(gam); 
        return 0;
}
