#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "function.h"
#include<iostream>
#include<time.h>
#include <locale.h>
using namespace ::std;

int n,m,k,i,j,p,rezult=0;
double *a = NULL,temp;
double* b = NULL;
double* x = NULL;
double* S = NULL;
int* Index ;

int main(int argc,char *argv[]){
    FILE * file;
    clock_t t;

    if ( argc != 5 && argc != 4 && argc != 1 )  {
        printf("incorrect data entered");
        return -1;}

    if ( argc == 1)
        return -1;

     n = atoi(argv[1]);
     m = atoi(argv[2]);
     k = atoi(argv[3]);

     if(m > n){
        m = n;
        printf("the output will be made by default (by max matrix size)\n");}

     if(k!=0 && argc != 4){
        printf("incorrect data entered");
        return -1;
     }

    if((k==0) &&(argc ==5)){

        file = fopen(argv[4],"r");
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

        if(p!=n*n){
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
    b = (double*)malloc(n * sizeof(double));
	x = (double*)malloc(n * sizeof(double));
	Index = (int*)malloc(n * sizeof(int));
	S = (double*)malloc(n * sizeof(double));

    if (!(Index && a && b && x)){
        fclose(file);
        free(a);
        free(b);
        free(x);
        return -1;
        }



    rezult = input(a,b,n,k,file); //input matr

    if (rezult==-1)
	{
		printf("Error in reading from file!\n");
        free(a);
		free(b);
		free(x);
		free(Index);
        return -1;
	}


    printf("\nMatrix A:\n");
    output(a, m, m);  //out matr
    printf("\n");
    printf("Vector B:\n");
    output(b,1,m);  //out vector b
    printf("\n");

    printf("Calculating...\n");
    t = clock();
    rezult = SolveSystem(n, a, b, x, Index, S);
    t = clock() - t;

    switch (rezult)
	{
	case -1:
		printf("Can't solve - the matrix is degenerate.\n");

		break;
	case 0:
		printf("\nrezulting vector x :\n");
		output(x,1,m);
		printf("\n");

		printf("Solution time\t\t= %.2f sec.\n\n", (double)t / CLOCKS_PER_SEC);

        printf("Norm of residuals: ||Ax - b||/||b||\t= %e\n", Error(n, a, b, S));

		printf("Solution accuracy\t= %e\n", accur(n, x));


        break;
	default:
		printf("Unknown error!\n");

		break;
	}



    free(a);
    free(b);
    free(x);
    free(Index);
    return 0;
}

