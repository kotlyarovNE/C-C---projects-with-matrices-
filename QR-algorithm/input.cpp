#include "function.h"
#include <stdio.h>
#include<math.h>
#include<stdlib.h>


int input(double* a, int n, int k, FILE *fin){
int i,j;
//Заполняем матрицу по функции
if(k!=0){
    for(i=0; i < n; i++){
            for(j=0; j < n; j++){
                a[i*n+j] = f(k,n,i+1,j+1);
            }
    }
}
//Заполняем матрицу из файлика
else{

    for (i = 0; i < n; ++i)
		{
			for (j = 0; j < n; ++j){
				if (fscanf(fin, "%lf", &a[i * n + j]) != 1)
					return -1;}
        }
fclose(fin);
}

return 0;
}


double DIFF(double *a, double tr, int n){

double tmp = 0;

    for(int i = 0; i < n; i++){
        tmp += a[i];}

return fabs(tr - tmp);
}

double trace(double *a, int n){
double tmp = 0;

    for(int i = 0; i < n; i++){
        tmp += a[i * n + i];}

return tmp;
}

double second_INVAR(double *a, double *b, int n){
int i;
double tmp = 0, pmt = 0;

    for(i = 0; i < n * n; i++){
        tmp += a[i] * a[i];
    }

    for(i = 0; i < n; i++){
        pmt += b[i] * b[i];
    }

return fabs(tmp - pmt);
}

double NORM_MATRIX(double *a, int n, int k){
int i,j;
double norm = 0, pmt = 0;

    for(i = 0; i < k; i++){
        norm += fabs(a[i]);
    }

    for(i = 1; i < k; i++){
        for(j = 0; j < k; j++){
            pmt += fabs(a[i * n + j]);
        }
        if(pmt > norm){
            norm = pmt;
        }
        pmt = 0;
    }
return norm;
}
