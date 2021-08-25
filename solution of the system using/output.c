#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "functions.h"


void OutputMatrix(int n, int m, double *a, double *b, double *x, int my_rank, int p)
{
	int i, j;
	MPI_Status status;


	for (i = 0; i < m; i++)
	{
		if (my_rank == 0)
		{
			if (my_rank == i%p)
			{
				printf("| ");
				for (j = 0; j < m; j++)
					printf("%10.3e ", a[i/p * n + j]);
				printf("|   %10.3e\n", b[i/p]);
			}
			else
			{
				MPI_Recv(x, m + 1, MPI_DOUBLE, i%p, 0, MPI_COMM_WORLD, &status);
				printf("| ");
				for (j = 0; j < m; j++)
					printf("%10.3e ", x[j]);
				printf("|   %10.3e\n", x[m]);
			}
		}
		else if (my_rank == i%p)
		{
			for (j = 0; j < m; j++)
				x[j] = a[i/p * n + j];
			x[m] = b[i/p];
			MPI_Send(x, m + 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}
}

void OutputVector(int m, double *b, double *x, int my_rank, int p)
{
	int i;
	MPI_Status status;

	for (i = 0; i < m; i++)
	{
		if (my_rank == 0)
		{
			if (my_rank == i%p)
				printf("%10.3g ", b[i/p]);
			else
			{
				MPI_Recv(x, 1, MPI_DOUBLE, i%p, 0, MPI_COMM_WORLD, &status);
				printf("%10.3g ", x[0]);
			}
		}
		else if (my_rank == i%p)
		{
			x[0] = b[i/p];
			MPI_Send(x, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}
}


void output(double* a, int strok, int stolb, int n){
    int i,j;

    for(i = 0; i < strok; i++){
        for(j = 0; j < stolb; j++){
          printf("%lf ",a[i*n+j]);
        }
        printf("\n");
    }
}



void S(double *b, int *index, double *b_buf, int my_rank, int p, int n, double *b_rez){
int i;
double x[1];
MPI_Status status;

for (i = 0; i < n; i++)
	{
		if (my_rank == 0)
		{
			if (my_rank == i%p){

                    b_buf[i] = b[i/p];
            }
			else
			{
				MPI_Recv(x, 1, MPI_DOUBLE, i%p, 0, MPI_COMM_WORLD, &status);
                b_buf[i] = x[0];
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
            MPI_Send(b_buf, n, MPI_DOUBLE, k, 0, MPI_COMM_WORLD);
        }

        if(my_rank == k){
            MPI_Recv(b_rez, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
            for(i = 0;i < n; i++)
                b_buf[i] = b_rez[i];
        }
	}



	for(i = 0; i < n; i++){
        b_rez[index[i]] = b_buf[i];
	}
}
