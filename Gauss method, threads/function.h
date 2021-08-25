#include <stdio.h>
#include<math.h>
#include <stdlib.h>


int input(double* a,double* b, int n, int k, FILE *fin);
double Error(int n, double* a, double* b, double* x);
double accur(int n, double* x);
void output(double* a, int m, int n);
double f(int k, int n, int i, int j);
//int SolveSystem(int n, double *a, double *b, double *x, double *y, int *index, int my_rank, int total_threads, maxElem *maxx, int *flag);
long int get_time(void);
long int get_full_time(void);
void Norm(double* a, double *b, int n);
//void synchronize(int total_threads);

