#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "functions.h"
#include <math.h>

double Error(int n, double* a, double *B, double *vect, int my_rank, int p) //Норма невязки, запотел...
{
	int i,j;
	double rezult1 = 0; //||Ax - b||
	double rezult2 = 0.0; //||B||
	double tmp = 0, x[1], y[1];
	double Nev_error = 0;
	MPI_Status status;

    for (i = 0; i < n; i++)
	{
		if (my_rank == 0)
		{
			if (my_rank == i%p){
                tmp = 0.0;
                for(j = 0; j < n; j++){
                    tmp += a[i/p * n + j] * B[j];}
                tmp -= vect[i/p];
                rezult1 += tmp * tmp;
            }
			else
			{
				MPI_Recv(x, 1, MPI_DOUBLE, i%p, 0, MPI_COMM_WORLD, &status);
                rezult1 += x[0];
			}
		}
		else if (my_rank == i%p)
		{
		    tmp = 0.0;
		    for(j = 0; j < n; j++){
                    tmp += a[i/p * n + j] * B[j];}
            tmp -= vect[i/p];
            x[0] = tmp * tmp;
			MPI_Send(x, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}


    MPI_Barrier(MPI_COMM_WORLD);
	//Теперь посчитаем норму b:

	for (i = 0; i < n; i++)
	{
		if (my_rank == 0)
		{
			if (my_rank == i%p){
                rezult2 += vect[i/p]*vect[i/p];
            }
			else
			{
				MPI_Recv(y, 1, MPI_DOUBLE, i%p, 0, MPI_COMM_WORLD, &status);
                rezult2 += y[0];
			}
		}
		else if (my_rank == i%p)
		{
			y[0] = vect[i/p]*vect[i/p];
			MPI_Send(y, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}



	if(my_rank == 0)
    {
        Nev_error = sqrt(rezult1/rezult2);
        printf("%10.3e\n ", Nev_error);
	}


return Nev_error;
}


double accur(int n, double* x, int my_rank)
{
	double tmp;
	int i;
	double rez;

	rez = 0.0;
	for (i = 0; i < n; i++)
	{
		tmp = x[i];

		if (i % 2 == 0)
			tmp -= 1.0;

		rez +=tmp*tmp;
	}

	if(my_rank == 0)
        printf("%lf ", sqrt(rez));

	return sqrt(rez);
}


void NORM(int n, double *a, double *b, int p, int my_rank){
int i, j, rows;
double NORM = 0, tmp = 0, final_norm = 0;

    if (my_rank + 1 > n%p) rows = n/p;
	else rows = n/p + 1;

	for(i = 0; i < n; i++)
        NORM += fabs(a[i]);


	for(i = 1; i < rows; i++){
        tmp = 0;
        for(j = 0; j < n; j++){
            tmp += fabs(a[i * rows + j]);
        }
        if(tmp > NORM)
            NORM = tmp;
	}

	MPI_Allreduce(&NORM, &final_norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if(final_norm < 1e-13){
        for(i = 0; i < n; i++){
            if(my_rank == i%p){
                for(j = 0; j < n; j++){
                    a[i/p * n + j] /= final_norm;
                }
                b[i/p] /= final_norm;
            }
        }
    }

}
