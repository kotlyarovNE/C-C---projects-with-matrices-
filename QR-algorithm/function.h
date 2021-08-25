#include <stdio.h>
#include<math.h>
#include <stdlib.h>

int input(double* a, int n, int k, FILE *fin);

void output(double* a, int STROK, int STOLB,  int n);

double f(int k, int n, int i, int j);



double DIFF(double *a, double tr, int n);

double second_INVAR(double *a, double *b, int n);

double trace(double *a, int n);

double NORM_MATRIX(double *a, int n, int k);

void SHIFT(double *a, double s, int n, int k);

double quadratic_equation(int k, double b, double c);

void MINUS_SIMMETRIX(double *a, int n);


void QR(double *R, double *a1, double *a2, double *a3, double *cosPhi, double *sinPhi, int n, int k, double e);

void Rotation(int n, double* a);

void SOLV(double *A, int n, double *a1, double *a2, double *a3, double *cosPhi, double *sinPhi, double *eigenvalues, double e, int k, int m, double *E);

void MULT(double *a, double *cosPhi, double *sinPhi, int n, int j);

void MULT1(double *a, double *cosPhi, double *sinPhi, int n);



void id(double *E, int n);

void WATCH_PROCESS(double *a, int n, int m);

void WATCH_QR(double *a, double *cosPhi, double *sinPhi, double *E, int n, int m, int i);

void COPY(double *a, double *E, int n);

