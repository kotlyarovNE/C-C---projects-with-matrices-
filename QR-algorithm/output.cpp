#include "function.h"
#include <stdio.h>
#include<math.h>
#include <stdlib.h>

void output(double* a, int STROK, int STOLB,  int n){
    int i, j;

    for(i = 0; i < STROK; i++){
        for(j = 0; j < STOLB; j++){
          printf("%10.3e ",a[i * n + j]);
        }
        printf("\n");
    }
}

void COPY(double *a, double *E, int n){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            E[i * n + j] = a[i * n + j];
        }
    }
}

void id(double *E, int n){
//Заведу единичную матрицу, она нигде не используется в алгоритме, просто для удобства вывода матрицы Q в QR разложении
	for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i == j){ E[i * n + j ] = 1; }
            else{ E[i * n + j] = 0; }
        }
	}
}

void WATCH_PROCESS(double *a, int n, int m){
    printf("\nMatrix A_k = R_k*Q_k :\n");
    output(a, m, m, n);  //out matr*/
}

void WATCH_QR(double *a, double *cosPhi, double *sinPhi, double *E, int n, int m, int i){
    printf("\nStep %d:", i + 1);
    id(E,n);
    //Вывод матрицы R:
    printf("\nMatrix R:\n");
    output(a, m, m, n);  //out matr

    //Вывожу матрицу Q
    printf("\n");
    printf("\nMatrix Q:\n");
    MULT1(E, cosPhi, sinPhi, n);
    output(E, m, m, n);

    COPY(a,E,n);

    //Проверка QR алгоритма
    MULT1(E, cosPhi, sinPhi, n);
    printf("\nMatrix A = Q*R:\n");
    output(a, m, m, n);  //out matr
}
