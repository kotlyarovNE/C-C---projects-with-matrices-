#ifndef struct_hpp
#define struct_hpp

#include <stdio.h>

typedef struct
{
    double elem;
    int rowIndex;
    int colIndex;
} maxElem;

int SolveSystem(int n, double *a, double *b, double *x, double *y, int *index, int my_rank, int total_threads, maxElem *maxx, int *flag);
void synchronize(int total_threads);

#endif /* struct_hpp */
