#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include "functions.h"
#include<string.h>


struct Local
{
    double max_loc;
    int thread;
};


int solve(int n, double *a, int my_rank, int p, double *b, double *vect, int *buffer_Stolb, double *BUFFER){
int i , j, k, rows, first, det[1], loc_i;
int indMax1, indMax2;
double tmp;
MPI_Status status;
struct Local loc, glob;


    det[0] = 0;
    loc.thread = my_rank;

    if (my_rank + 1 > n%p) rows = n/p;
	else rows = n/p + 1;

	for(i = 0; i < n; i++)
        buffer_Stolb[i] = i;



        for(i = 0; i < n; i++){
            indMax1 = i; indMax2 = i;
            if(my_rank == i%p){
                loc.max_loc = fabs(a[i/p * n + i]);
                for(j = i; j < n; j++){
                    for(k = i; k < n; k++){
                        if(fabs(a[j/p * n + k]) > fabs(loc.max_loc)){
                            indMax1 = j;
                            indMax2 = k;
                            loc.max_loc = fabs(a[j/p * n + k]);
                        }
                    }
                }
            }

            else{
                if (my_rank > i%p) first = i;
                else first = i + p;
                loc.max_loc = fabs(a[first/p * n + i]);
                for(j = first; j < n; j++){
                    for(k = i; k < n; k++){
                        if(fabs(a[j/p * n + k]) > fabs(loc.max_loc)){
                            indMax1 = j;
                            indMax2 = k;
                            loc.max_loc = fabs(a[j/p * n + k]);
                        }
                    }
                }
            }


        MPI_Allreduce(&loc, &glob, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);//Нащли максимум

        MPI_Bcast(&indMax1, 1, MPI_INT, glob.thread, MPI_COMM_WORLD);

        MPI_Bcast(&indMax2, 1, MPI_INT, glob.thread, MPI_COMM_WORLD);

        MPI_Bcast(&buffer_Stolb[i], 1, MPI_INT, glob.thread, MPI_COMM_WORLD);

        MPI_Bcast(&buffer_Stolb[indMax2], 1, MPI_INT, glob.thread, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        j = buffer_Stolb[i];

		buffer_Stolb[i] = buffer_Stolb[indMax2];

		buffer_Stolb[indMax2] = j;

		if(my_rank == 0){
            if(fabs(glob.max_loc) < 1e-13){
                det[0] = 1;
            }
		}

		MPI_Bcast(&det[0], 1, MPI_INT, i%p, MPI_COMM_WORLD);
		if(det[0] == 1)
            return -1;

        //Меняем местами столбцы
        for (j = 0; j < rows; j++)
            {
                tmp = a[j * n + i];
                a[j * n + i] = a[j * n + indMax2];
                a[j * n + indMax2] = tmp;
            }





    if(indMax1!=i){

        if ((glob.thread == (i % p)) && (glob.thread == my_rank)) // Если строчки лежат в одном процессе
        {
            for (j = 0; j < n; j++)
            {
                tmp = a[i / p * n + j];
                a[i / p * n + j] = a[indMax1 / p * n + j];
                a[indMax1 / p * n + j] = tmp;
            }

            tmp = b[i / p];
            b[i / p] = b[indMax1 / p];
            b[indMax1 / p] = tmp;

        }
        else if (glob.thread == my_rank)
        {
            loc_i = indMax1 / p;
            memcpy(BUFFER, a + loc_i * n, n * sizeof(double));
            memcpy(BUFFER + n, b + loc_i, 1 * sizeof(double));

            MPI_Sendrecv_replace(BUFFER, n + 1, MPI_DOUBLE, i % p, 0, i % p, 0, MPI_COMM_WORLD, &status);

            memcpy(a + loc_i * n, BUFFER, n * sizeof(double));
            memcpy(b + loc_i , BUFFER + n, 1 * sizeof(double));
        }
        else if ((i % p) == my_rank)
        {
            loc_i = i / p;
            memcpy(BUFFER, a + loc_i * n, n * sizeof(double));
            memcpy(BUFFER + n, b + loc_i, 1 * sizeof(double));

            MPI_Sendrecv_replace(BUFFER, n + 1, MPI_DOUBLE, glob.thread, 0, glob.thread, 0, MPI_COMM_WORLD, &status);

            memcpy(a + loc_i * n, BUFFER, n * sizeof(double));
            memcpy(b + loc_i, BUFFER + n, 1 * sizeof(double));
        }


    }









        MPI_Barrier(MPI_COMM_WORLD);


        //ТЕПЕРЬ САМО РЕШЕНИЕ ПОЙДЕТ
        if (my_rank == i%p) //если у нас сейчас (i mod p) процесс
		{
            tmp = 1.0/a[i/p * n + i];
			for (j = i; j < n; j++)
				a[i/p * n + j] *= tmp;
			b[i/p] *= tmp;
			vect[i/p] = b[i/p];

			if (i == n - 1) continue; //подошли к концу и дальше не надо пересылать

            //надо записать в буфер терерь
			for (j = i; j < n; j++)
				BUFFER[j] = a[i/p * n + j]; //записываем в буфер строчку матрицы размера (n-i)

			BUFFER[n] = b[i/p]; //помимо матрицы и сам вектор b надо закинуть в буфер

			MPI_Bcast(BUFFER, n + 1, MPI_DOUBLE, i%p, MPI_COMM_WORLD); //идёт перессылка всем процессам
			for (j = i/p + 1; j < rows; j++) //обнуляем подстолбец размера rows
			{
				tmp = a[j * n + i];
				for (k = i; k < n; k++)
					a[j * n + k] -= tmp * a[i/p * n + k]; //сдандартно вычитаем строки
				b[j] -= tmp * b[i/p];//вычитаем и вектор тоже
				vect[j] = b[j];
			}
		}
		else // если у нас не (i mod p) процесс
		{
			if (i == n - 1) continue; //последний щаг не нужно дальше пересылать ничего это конец


			//первые строчки процесса, которые работали на прошлых шагах цикла в if сдвигаем на 1 позицию вперед
			//ну потому что размер матрицы уменьшился сам по себе на 1 по сравнению с предыдущим шагом
			//а мы матрицу храним по блокам через p-1 строк на каждом узле
			//таким образом, мы постепенно движемся к низу

			if (my_rank > i%p) first = i/p;
			else first = i/p + 1;

			MPI_Bcast(BUFFER, n + 1, MPI_DOUBLE, i%p, MPI_COMM_WORLD);//эта функция должна быть вызвана во всех процесса группы com с одинаковым значенем root
			for (j = first; j < rows; j++)
			{
				tmp = a[j * n + i];
				for (k = i; k < n; k++)
					a[j * n + k] -= tmp * BUFFER[k];//вычитаем то первую строчку
				b[j] -= tmp * BUFFER[n]; //аналогично работаем и с вектором
                vect[j] = b[j];
			}
		}



    MPI_Barrier(MPI_COMM_WORLD);
    }//ТУТ ЗАКАНЧИВАЕТСЯ ОСНОВНОЙ ЦИКЛ

    //Обратный ход пошел

    for (i = n - 1; i >= 1; i--)
	{
		if (my_rank == i%p)
		{
			MPI_Bcast(&b[i/p], 1, MPI_DOUBLE, i%p, MPI_COMM_WORLD);
			for (j = i/p - 1; j >= 0 ; j--)
				b[j] -= b[i/p] * a[j * n + i];
		}
		else
		{
			MPI_Bcast(BUFFER, 1, MPI_DOUBLE, i%p, MPI_COMM_WORLD);

			if (my_rank < i%p) first = i/p;
			else if (i/p - 1 >= 0) first = i/p - 1;
			else continue;

			for (j = first; j >= 0; j--)
				b[j] -= BUFFER[0] * a[j * n + i];
		}
	}


return 0;
}

