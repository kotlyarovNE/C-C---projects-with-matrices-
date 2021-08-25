#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "function.h"
#include<iostream>
#include<time.h>
#include <locale.h>
using namespace ::std;

int n, m, k, p, rezult = 0;
double *a = NULL, *a1 = NULL, *a2 = NULL, *a3 = NULL, temp = 0, tr = 0;
double *cosPhi = NULL, *sinPhi = NULL, *eigenvalues = NULL, *E = NULL, *tmp = NULL;
double e = 0;

int main(int argc,char *argv[]){
    FILE *file = NULL;
    clock_t t;

    if ( argc != 6 && argc != 5 && argc != 1 )  {
        printf("incorrect data entered");
        return -1;}

    if ( argc == 1 )
        return -1;

     n = atoi(argv[1]);
     m = atoi(argv[2]);
     e = atof(argv[3]);
     k = atoi(argv[4]);

     if(m > n){
        m = n;
        printf("the output will be made by default (by max matrix size)\n");}

     if(k!=0 && argc != 5){
        printf("incorrect data entered");
        return -1;
     }

    if((k==0) &&(argc ==6)){

        file = fopen(argv[5],"r");
		if (!file)
		{
			printf("Can't open file!\n");

			return -1;
		}

        if(file==NULL)
        {
            return -1;
        }

		while(fscanf(file, "%lf", &temp)==1)
            p++;

        if(p != n*n){
            fclose(file);
            printf("not enough data in the file || too much data in the file\n");
        }

        if(!feof(file)){
            fclose(file);
            return -1;
        }

        fseek(file, 0, 0);
    }


    a = (double*)malloc(n * n * sizeof(double));
    E = (double*)malloc(n * n * sizeof(double));
    a1 = (double*)malloc((n-1)*sizeof(double));
    a3 = (double*)malloc((n-1)*sizeof(double));
    a2 = (double*)malloc((n)*sizeof(double));
    cosPhi = (double*)malloc((n-1)*sizeof(double));
    sinPhi = (double*)malloc((n-1)*sizeof(double));
    eigenvalues = (double*)malloc((n)*sizeof(double));
    tmp = (double*)malloc((n * n)*sizeof(double));

    if (!(a)){
        fclose(file);
        free(a);
        free(a1);
        free(a2);
        free(a3);
        free(cosPhi);
        free(sinPhi);
        free(E);
        free(eigenvalues);
        free(tmp);
        return -1;
        }



    rezult = input(a,n,k,file); //ввод матрицы по формуле или из файла


    if (rezult == -1)
	{
		printf("Error in reading from file!\n");
        free(a);
        free(a1);
        free(a2);
        free(a3);
        free(cosPhi);
        free(sinPhi);
        free(E);
        free(eigenvalues);
        free(tmp);
        return -1;
	}

    COPY(a, tmp, n);
    //Вывод матрицы А
    printf("\nMatrix A:\n");
    output(a, m, m, n);  //out matr
    printf("\n");

    tr = trace(a, n); //след исходной матрицы

    printf("Calculating...\n");
    printf("\n");

    t = clock();
    id(E,n);
    SOLV(a, n, a1, a2, a3, cosPhi, sinPhi, eigenvalues, e, k, m, E);
    //Rotation(n, a);
    //output(a, m, m, n);
    t = clock() - t;



    switch (rezult)
	{
	case -1:

        printf("Error in reading from file!\n");
        break;

	case 0:

        printf("\n");
        printf("Eigenvalues of the matrix:\n");
        for(int i = 0; i < m; i++){
            printf("    %.20lf\n",eigenvalues[i]);}

        printf("Diffrence between TRACE and SUM of eigenvalues : %10.3e\n", DIFF(eigenvalues,tr,n));
        printf("Second invariant: %10.3e \n", second_INVAR(tmp, eigenvalues, n));
		printf("Solution time\t\t= %.2f sec.\n\n", (double)t / CLOCKS_PER_SEC);
		//printf("norm A = %.20lf", NORM_MATRIX(tmp, n));
        break;

	default:
		printf("Unknown error!\n");
        break;
	}

    free(a);
    free(a1);
    free(a2);
    free(a3);
    free(cosPhi);
    free(sinPhi);
    free(E);
    free(eigenvalues);
    free(tmp);
    return 0;
}

