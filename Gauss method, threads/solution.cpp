#include <pthread.h>
#include <math.h>
#include "function.h"
#include "struct.hpp"

void synchronize(int total_threads)
{
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
	static int threads_in = 0;
	static int threads_out = 0;

	pthread_mutex_lock(&mutex);

	threads_in++;
	if (threads_in >= total_threads)
	{
		threads_out = 0;
		pthread_cond_broadcast(&condvar_in);
	} else
		while (threads_in < total_threads)
			pthread_cond_wait(&condvar_in,&mutex);

	threads_out++;
	if (threads_out >= total_threads)
	{
		threads_in = 0;
		pthread_cond_broadcast(&condvar_out);
	} else
		while (threads_out < total_threads)
			pthread_cond_wait(&condvar_out,&mutex);

	pthread_mutex_unlock(&mutex);
}

int SolveSystem(int n, double *a, double *b, double *x, double *y, int *index, int my_rank, int total_threads, maxElem *maxx, int *flag)
{
	int i, j, k;
	int  k1, k2;
	int first_row;
	int last_row;
	double tmp = 0, maxElem;
	
	if(my_rank == 0){
		for (i = 0; i < n; i++)
			index[i] = i;}

	for (i = 0; i < n; ++i)
	{
	//synchronize(total_threads);
        first_row = (n - i) * my_rank;
        first_row = first_row/total_threads + i;
        last_row = (n - i) * (my_rank + 1);
        last_row = last_row/total_threads + i;
        //k1 = i;
        //k2 = i;
        maxx[my_rank].elem = fabs(a[first_row * n + i]);
        maxx[my_rank].rowIndex = first_row;
        maxx[my_rank].colIndex = i;
        for (j = first_row; j < last_row; ++j){
            for (k = i; k < n; ++k){
                if (fabs( a[j * n + k] ) > maxx[my_rank].elem)
                {
                    maxx[my_rank].elem = fabs( a[j * n + k] );
                    maxx[my_rank].rowIndex = j;
                    maxx[my_rank].colIndex = k;
                }
            }
        }

        //mass_strok[2 * my_rank] = k1;
        //mass_strok[2 * my_rank+1] = k2;
        //mass_stolb[my_rank] = k2;

        synchronize(total_threads);

        if(my_rank == 0)
        {

            /*k1 = mass_strok[0];
            k2 = mass_strok[1];
            for(k = 0; k < total_threads - 2; k++){
                if( a[k1 * n + k2] < a[mass_strok[2*k+2]*n + mass_strok[2*k+3]] ){
                    k1 = mass_strok[2*k+2];
                    k2 = mass_strok[2*k+3];
                }
            }

			j = index[i];
			index[i] = index[k2];
			index[k2] = j;

            for(j = 0; j < total_threads * 2; j++){
                mass_strok[j] = 0;
            }*/
            maxElem = maxx[0].elem;
            k1 = maxx[0].rowIndex;
            k2 = maxx[0].colIndex;
            for(k = 1; k < total_threads; ++k)
            {
                if(maxElem < maxx[k].elem)
                {
                    maxElem = maxx[k].elem;
                    k1 = maxx[k].rowIndex;
                    k2 = maxx[k].colIndex;
                }
            }

            j = index[i];
	    index[i] = index[k2];
	    index[k2] = j;

            if(fabs(maxElem) < 1e-13){
                *flag = 1;} //det = 0;

            else{
			for (j = 0; j < n; j++)
			{
				tmp = a[i * n + j];
				a[i * n + j] = a[k1 * n + j];
				a[k1 * n + j] = tmp;
			}

			for (j = 0; j < n; j++)
			{
				tmp = a[j * n + i];
				a[j * n + i] = a[j * n + k2];
				a[j * n + k2] = tmp;
			}

			tmp = b[i];
			b[i] = b[k1];
			b[k1] = tmp;

			//if (fabs(a[i * n + i]) < 1e-13)
                	//*flag = 1;

            		tmp = a[i * n + i];

			tmp = 1.0/tmp;
			for (j = i; j < n; j++)
				a[i * n + j] *= tmp;
			b[i] *= tmp;
            }
	}


        synchronize(total_threads);
        if (*flag != 0) return -1;

        first_row = (n - 1 - i) * my_rank;
        first_row = first_row/total_threads + i + 1;
        last_row = (n - 1 - i) * (my_rank + 1);
        last_row = last_row/total_threads + i + 1;
        for (j = first_row; j < last_row; j++)
        {
            tmp = a[j * n + i];
            for (k = i; k < n; k++)
                a[j * n + k] -= tmp * a[i * n + k];
            b[j] -= tmp * b[i];
        }
        synchronize(total_threads);
    }
    
	//synchronize(total_threads);


    if (my_rank == 0)
        {
		for (i = n - 1; i >= 0; --i)
		{
			tmp = b[i];
			for (j = i + 1; j < n; ++j)
				tmp -= a[i * n + j] * x[index[j]];
			x[index[i]] = tmp;
		}

		for( i = 0; i < n; ++i){
            		y[i] = x[index[i]];}
        }

return 0;
}
