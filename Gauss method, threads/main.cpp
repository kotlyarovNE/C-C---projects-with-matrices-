#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "function.h"
#include "struct.hpp"
#include<iostream>
#include<time.h>
#include <pthread.h>
using namespace ::std;

// Аргyменты для потока

typedef struct
{
	int n;
	double *a;
	double *b;
	double *x;
	double *y;
	int *index;
	int my_rank;
	int total_threads;
	bool DEG_ORNOT;
	long int time;
	//int *mass_strok;
	//int *mass_stolb;
	//int I;
	maxElem *maxx;
	int *flag;
} ARGS;



long int thread_time = 0; 	//сyммарное время работы всех задач
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER; //объект типа мютекс для синра к общ. времени

void *Solution(void *p_arg) //решение для одной задачи
{
	ARGS *arg = (ARGS*)p_arg;
	long int t1;

	synchronize(arg->total_threads);
	t1 = get_full_time();
	if (SolveSystem(arg->n, arg->a, arg->b, arg->x, arg->y, arg->index, arg->my_rank, arg->total_threads, arg->maxx, arg->flag) == -1 ){
		arg->DEG_ORNOT = true;
		}
    synchronize(arg->total_threads);

	t1 = get_full_time() - t1;	//выдает время затраченное на тред
	arg->time = t1;
	thread_time+=t1;
	synchronize(arg->total_threads);

	return NULL;
}

int main(int argc,char *argv[]){

    FILE *file = NULL;
    int n, m, k, i, p = 0, rezult=0, zero = 0;
    double*a = NULL, temp = 0;
    double* b = NULL;
    double* x = NULL;
    double* y = NULL;
    int* index = NULL;
    ARGS *args;         //массив аргyментов для созданных задач
    int total_threads;
    pthread_t *threads; // массив идентефикаторов созданных задач
    maxElem *maxx;
    //int *mass_strok = NULL;
	//int *mass_stolb = NULL;

    if ( argc != 6 && argc != 5 && argc != 1 )  {
        printf("incorrect data entered\n");
        return -1;}

    if ( argc == 1)
        return -1;

     n = atoi(argv[1]);
     m = atoi(argv[2]);
     k = atoi(argv[3]);
     total_threads = atoi(argv[4]);

     if(total_threads < 1){
		printf("incorrect numbers of threads\n");
		return -1;}

     if(m > n){
        m = n;
        printf("the output will be made by default (by max matrix size)\n");}

     if(k!=0 && argc != 5){
        printf("incorrect data entered\n");
        return -1;
     }

    if((k==0) && (argc!=6)){
	printf("file name not entered!\n");
	return -1;}

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
    y = (double*)malloc(n * sizeof(double));
    index = (int*)malloc(n * sizeof(int));
    threads = (pthread_t*)malloc(total_threads * sizeof(pthread_t));
    //maxx=new maxElem[threadsCount];
    maxx = (maxElem*)malloc(total_threads * sizeof(maxElem));;
    args = (ARGS*)malloc(total_threads * sizeof(ARGS));

    //mass_strok = (int*)malloc(2*total_threads * sizeof(int));
    //mass_stolb = (int*)malloc(total_threads * sizeof(int));

    if (!(index && a && b && x && threads && args && y && maxx))
    {
        printf("Not enough memory!\n");
        fclose(file);
        free(a);
        free(b);
        free(x);
        free(y);
        free(threads);
        free(args);
        free(index);
        //free(mass_strok);
        //free(mass_stolb);
        free(maxx);
        return -1;
    }



    rezult = input(a,b,n,k,file); //input matr

    if (rezult == -1){
		printf("Error in reading from file!\n");
        free(a);
		free(b);
		free(x);
		free(y);
		free(index);
		free(threads);
		free(args);
		free(maxx);
		//free(mass_strok);
        //free(mass_stolb);
        return -1;
	}


    printf("\nMatrix A:\n");
    output(a, m, m);
    printf("\n");
    printf("Vector B:\n");
    output(b,1,m);
    printf("\n");
    Norm(a,b,n);
    printf("Calculating...\n");

    //инициализация аргyментов задач
    for (i = 0; i < total_threads; i++)
	{
		args[i].n = n;
		args[i].a = a;
		args[i].b = b;
		args[i].x = x;
		args[i].y = y;
		args[i].index = index;
		args[i].my_rank = i;
		args[i].total_threads = total_threads;
		args[i].DEG_ORNOT = false;
		args[i].time = 0;
		//args[i].mass_strok = mass_strok;
		//args[i].mass_stolb = mass_stolb;
		args[i].maxx = maxx;
		args[i].flag = &zero;
	}

    for (i = 0; i < total_threads; i++) //запyскаем задачи
		if (pthread_create(threads + i, 0, Solution, args + i))
		{
			printf("Cannot create thread %d!\n", i);

			if (a) free(a);
			if (b) free(b);
			if (x) free(x);
			if (y) free(y);
			if (index) free(index);
			if (threads) free(threads);
			if (args) free(args);
			//if (mass_strok) free(mass_strok);
			//if (mass_stolb) free(mass_stolb);
			if (maxx) free(maxx);
			return -1;
		}

	for (i = 0; i < total_threads; i++)
		if (pthread_join(threads[i], 0))    //ожидаем окончания задач
		{
			printf("Cannot wait thread %d!\n", i);

			if (a) free(a);
			if (b) free(b);
			if (x) free(x);
			if (y) free(y);
			if (index) free(index);
			if (threads) free(threads);
			if (args) free(args);
			//if (mass_strok) free(mass_strok);
			//if (mass_stolb) free(mass_stolb);
			if (maxx) free(maxx);
            return -1;
		}

	//Проверка на вырожденность
	p = 0;
	for(i = 0; i < total_threads; i++){
		if(args[i].DEG_ORNOT){
			p = -1;
			break;}
	}

  switch (p)
	{
	case -1:
		printf("Can't solve - the matrix is degenerate.\n");

		break;
	case 0:
		printf("\nrezulting vector x :\n");
		output(x,1,m);
		printf("\n");
        printf("Norm of residuals: ||Ax - b||/||b||\t= %e\n", Error(n, a, b, y));
        printf("Solution accuracy\t= %e\n", accur(n, x));

		printf("Time = %ld\n",args[total_threads-1].time );
		printf("All total tred time / total threads = %ld\n", thread_time/total_threads);
		break;
	default:
		printf("Unknown error!\n");
		break;
	}

    free(a);
    free(b);
    free(x);
    free(index);
    free(threads);
    free(args);
    free(y);
    //free(mass_strok);
    //free(mass_stolb);
    free(maxx);
    return 0;
}

