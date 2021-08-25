#include "function.h"
#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include <sys/resource.h>
#include <sys/time.h>

int input(double* a,double* b,int n, int k, FILE *fin){
int i,j;
double tmp;
//Заполняем матрицу по функции
if(k!=0){
    for(i=0;i < n; i++){
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


//Заполняем вектор b по странной формуле
for (i = 0; i < n; ++i)
		{
			tmp = 0.0;
			for (j = 0; j < n; ++j)
			{
                if (j % 2 == 0)
                tmp += a[i * n + j];
			}
			b[i] = tmp;
		}

return 0;
}



double accur(int n, double* x)
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

	return sqrt(rez);
}




double Error(int n, double* a, double* b, double* x)
{
	int i;
	int j;
	double rezult1 = 0.0;
	double rezult2 = 0.0;
	double tmp;

	for (i = 0; i < n; i++)
	{
		tmp = 0.0;
		for (j = 0; j < n; ++j)
			tmp += a[i * n + j] * x[j];
		tmp -= b[i];

		rezult1 +=tmp*tmp;
	}

	for(i=0; i<n; i++){
        rezult2+=b[i]*b[i];
	}
	return sqrt(rezult1/rezult2);
}

long int get_time(void)
{
	struct rusage buf;

	getrusage(RUSAGE_SELF, &buf);

	return buf.ru_utime.tv_sec * 100 + buf.ru_utime.tv_usec/10000;
}

long int get_full_time(void)
{
	struct timeval buf;

	gettimeofday(&buf, 0);

	return buf.tv_sec * 100 + buf.tv_usec/10000;
}


void Norm(double* a, double *b, int n){
double NormA = 0, tmp = 0;
int i,j;	
for(j = 0; j < n; j++)
	NormA += fabs(a[j]);
	
for(i = 1; i < n; i++){
	for(j = 0; j < n; j++){
		tmp+=fabs(a[i*n+j]);}
	if(NormA <= tmp) 
		NormA = tmp;
	tmp = 0;
}
		
	
if(fabs(NormA) <= 1e-13){
	for(i = 0; i < n*n; i++)
		a[i] = a[i]/NormA;

	for(i = 0; i < n; i++)
		b[i] = b[i]/NormA;
		}	
}



