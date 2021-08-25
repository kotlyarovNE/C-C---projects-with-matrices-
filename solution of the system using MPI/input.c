#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "functions.h"
#include <mpi.h>
#include<math.h>

void input(int k, int n, double *A, double *B, int my_rank, int p) {
	int i = 0, j = 0, rows;
    double x;

	if (my_rank + 1 > n%p) rows = n/p;
	else rows = n/p + 1;

	for (i = 0; i < rows; i++)
	{
		for (x = 0.0, j = 0; j < n; j++){
            A[i * n + j] = f(k, n, my_rank + p * i, j);
			if ((j % 2) == 0) x += A[i * n + j];
		}
		B[i] = x;
	}
}


void input_from_file(int n, double *A, double *a, double *b, int my_rank, int p){
    int i,j;
    int rows ;
    double x;

    if (my_rank + 1 > n%p) rows = n/p;
	else rows = n/p + 1;

	for (i = 0; i < rows; i++)
	{
		for (x = 0.0, j = 0; j < n; j++)
		{
			a[i * n + j] = A[(my_rank + p * i) * n + j];
			if ((j % 2) == 0) x += a[i * n + j];
		}
		b[i] = x;
	}
}






double f(int k, int n, int i, int j)
{
	if (k == 1) return (n - max(i+1, j+1) +1);
	if (k == 2) return (max(i+1, j+1)) ;
	if (k == 3) return max(i-j, j-i);
	if (k == 4) return 1.0/(i+j+1);
	else return 0;
}




int matr_file(int n, char *filename, double *A){
	FILE *IN ;

	IN = fopen(filename, "r");
	if (IN == NULL)
		return -1;

    for (int i = 0; i < n * n; i++) {
		if (fscanf(IN, "%lf", &A[i]) != 1) {
			fclose(IN);
			return -2;
		}
	}

	fclose(IN);
	return 0;
}



int max(int x, int y) {
	if ( x > y ) return x;
	else return y;
}





void B_Matr(double *a, double*b, int n, int p, int my_rank){
int i, j, first_row, last_row;
double x;

    first_row = n * my_rank;
	first_row /= p;

	last_row = n * (my_rank + 1);
	last_row = last_row / p - 1;

for (i = 0; i < n; i++)
	{
		for (x = 0.0, j = 0; j < n; j++){
			if ((j % 2) == 0) x += a[i * n + j];
		}
		b[i] = x;
	}
}


void inp_A_B(int n, int k, double *A, double *B){
int i,j;
double tmp;

    for(i=0;i < n; i++){
                for(j=0; j < n; j++){
                    A[i*n+j] = f(k,n,i,j);
                }
        }

    for (i = 0; i < n; ++i)
            {
                tmp = 0.0;
                for (j = 0; j < n; ++j)
                {
                    if (j % 2 == 0)
                    tmp += A[i * n + j];
                }
                B[i] = tmp;
            }

}

void input_b_after_solution(double *B, double *b, int n, int my_rank, int p, double *buf){
int i;
double x[1];
MPI_Status status;


    for (i = 0; i < n; i++)
	{
		if (my_rank == 0)
		{
			if (my_rank == i%p){

                    B[i] = b[i/p];
            }
			else
			{
				MPI_Recv(x, 1, MPI_DOUBLE, i%p, 0, MPI_COMM_WORLD, &status);
                B[i] = x[0];
			}
		}
		else if (my_rank == i%p)
		{
			x[0] = b[i/p];
			MPI_Send(x, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}



	for (int k = 1; k < p; k++){
	    if(my_rank == 0){
            MPI_Send(B, n, MPI_DOUBLE, k, 0, MPI_COMM_WORLD);
        }

        if(my_rank == k){
            MPI_Recv(buf, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
            for(i = 0;i < n; i++)
                B[i] = buf[i];
        }
	}
}



