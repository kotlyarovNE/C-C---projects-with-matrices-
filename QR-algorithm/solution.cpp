#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "function.h"

void MINUS_SIMMETRIX(double *a, int n){ //��� ���� ������� ������������ ���, �� ��� �������� �� ��������� ��������� �����.
int i;
    for(i = 0; i < n; i++){
        a[i] *= -1;
    }
}

void SHIFT(double *a, double s, int n, int k){
int i;
    for(i = 0; i < k; i++){
        a[i * n + i] -= s;
    }
}

double quadratic_equation(int k, double b, double c){
double MAX = 0, MIN = 0, tmp;
    MAX = ((-b+sqrt(b*b-4*c))/2);
    MIN = ((-b-sqrt(b*b-4*c))/2);

    if(MAX < MIN){
        tmp = MAX;
        MAX = MIN;
        MIN = tmp;
    }

    if(k == 1){
        return MAX;
    }
    if(k == 2){
        return MIN;
    }
return 0;
}

void SOLV(double *A, int n, double *a1, double *a2, double *a3, double *cosPhi, double *sinPhi, double *eigenvalues, double e, int k, int m, double *E){
    double Norm_matr, s;
    //Norm_matr = NORM_MATRIX(A, n, n);
    m += 0;
    k += 0;
    //Rotation(n, A);
    //Norm_matr = NORM_MATRIX(A, n, n);//����� �������
    //������� ������� � ���������������� ���� �������
    //if(k != 2){Rotation(n, A);}
    Rotation(n, A);
    //else{MINUS_SIMMETRIX(A, n*n);}
    //Norm_matr = NORM_MATRIX(A,n);
    //����� �������� �����:
    //���� QR ���������� �������, Q_k*R_k, ����� ����� � ������������: R_k*Q_k = R_k * (�(T_k,k+1)^T), k = 1...n-1
    //� ����� QR ���������� ��� ���� ����� �������

    //����������� ����� ���:
    //����� ���������� ����������� �������� �� i-�� ����, ����� n-2 �����, ����� ����� ���������� ���������
    //����� �� i-�� ���� ����� a_nn, ��������� QR ���������� � ���� �� ��������, � ,����� ������� a_n,n-1 = 0, ���������� �.�
    //� ��������� � ������� �� ������ ������
    id(E,n);
    for(int i = n; i > 2; --i){
        Norm_matr = NORM_MATRIX(A, n, n);//����� �������

        while(fabs(A[(i - 1) * n + i - 2]) > e * Norm_matr){ //��� �������� ������ �.� (a_n,n-1 = 0)

        //s = A[(i - 1) * n + i - 1];

        if(k!=2) {s = A[(i - 1) * n + i - 1];}//����� s_k

        else {s = 0;}

        SHIFT(A, s, n, i);//������� ����� �� i ����

        QR(A, a1, a2, a3, cosPhi, sinPhi, n, i, e);//������ QR ���������� ������� A_k-s_k*I

        //WATCH_QR(A, cosPhi, sinPhi, E, n, m, (n - i));//������ �� QR ������������

        MULT(A, cosPhi, sinPhi, n, i); //R_k*Q_k = R_k * (�(T_k,k+1)^T), k = 1...n-1

        SHIFT(A, -s, n, i); //������� ����� � �������� �������

        //WATCH_PROCESS(A, n, m);//�_� �������
        }

        eigenvalues[i - 1] = A[(i - 1) * n + i - 1];; //����������� �������� �� ��������� �������
    }

    //�������� ��� 2 �.�, ��� ������������ ��� ������� ����������� ���������

    eigenvalues[0] = quadratic_equation(1, -(A[0]+A[n+1]),(A[0]*A[n+1]-A[1]*A[n]) );
    eigenvalues[1] = quadratic_equation(2, -(A[0]+A[n+1]),(A[0]*A[n+1]-A[1]*A[n]) );

    //if(k == 2){MINUS_SIMMETRIX(eigenvalues, n);}
}



void QR(double *R, double *a1, double *a2, double *a3, double *cosPhi, double *sinPhi, int n, int k, double e){
    int i;
	double x;
	double y;
	double r;
	double cos;
	double sin;

	for (i = 0; i < k - 1; ++i){
        a3[i] = R[i * n + i + 1];
        a2[i] = R[i* n + i];
        a1[i] = R[(i + 1) * n + i];
	}
	a2[k - 1] = R[n * (k - 1) + k - 1];
	//� �������� ��������� ��������� ������ a2, a1 - ������ ������������, a3 -  �������

	for(i = 1; i < k; ++i){ // k = 1...n-1 steps,
        x = a2[i-1];
        y = a1[i-1];
        r = sqrt(x * x + y * y);

        if (fabs(y) < e)
            continue;

        r = sqrt(x * x + y * y);

        if (r < e){
            if(a2[i] < 0){
                cosPhi[i - 1] = 1;}
            else{
                cosPhi[i - 1] =  -1;}
			sinPhi[i - 1] = 0.0;
        }

        else{
            cos = x / r;
            sin = -y / r;
            cosPhi[i - 1] = cos; //���������� ������ ��� ������� Q - �������������
            sinPhi[i - 1] = sin;
        }

        R[i * n + i - 1] = 0; //�������� ������������
        R[(i - 1) * n + (i - 1)] = r; //������� ������� � ������� ����� ������

        //����� �������� �� T_ij �����, i = 1...n-1, j = i+1
        //�������� ������ i � i+1 ������
        //�� ���� ������ 4 �������� �������� �� i-�� ����, ��� ��� 2 ��� �� ��������
        R[(i - 1) * n + i] = cos * a3[i - 1] - sin * a2[i];
        R[i * n + i] = sin * a3[i - 1] + cos * a2[i];
        R[(i - 1) * n + i + 1] = -sin * a3[i];
        //�� ��������� ���� 2 �������� �� ����������
        if( i != (k - 1)){
            R[i * n + i + 1] = cos * a3[i];
            a3[i] = R[i * n + i + 1];
            }
        //�������� �������� ������� A, ������� �������� � ��� ��������
        a2[i] = R[i * n + i];
        }

}


void Rotation(int n, double* a)
{
	int i, j, k;
	double x, y, r, t, A_ii, A_ij, A_ji, A_jj, cosPhi, sinPhi;

	for (i = 1; i < n - 1; ++i) // n-2 steps
	{
		for (j = i + 1; j < n; ++j)
		{
			x = a[i * n + i - 1]; //������� �������
			y = a[j * n + i - 1]; //����� �������� ���
            r = sqrt(x * x + y * y);

			if (r < 1e-12) continue;
            cosPhi = x/r;
			sinPhi = -y/r;
            a[j * n + i - 1] = a[(i - 1) * n + j] = 0.0;
			a[i * n + i - 1] = a[(i - 1) * n + i] = r;

            for (k = i + 1; k < n; k++){ //�������� �� ������� ������������ ��������

				if (k == j) {continue;}
                x = a[i * n + k];
				y = a[j * n + k];
				a[k * n + i] = a[i * n + k] = x * cosPhi - y * sinPhi; //������� ������������
				a[k * n + j] = a[j * n + k] = x * sinPhi + y * cosPhi;}

			//4 �������� �� ii ij ji jj ������ �������� �� ��������
			x = a[i * n + i];
			y = a[j * n + j];
			r = a[i * n + j];
			t = a[j * n + i];

            //�����*
			A_ii = x * cosPhi - t * sinPhi;
			A_ji = x * sinPhi + t * cosPhi;
			A_ij = r * cosPhi - y * sinPhi;
			A_jj = r * sinPhi + y * cosPhi;
            //������*
			a[i * n + i] = A_ii * cosPhi - A_ij * sinPhi;
			a[j * n + i] = A_ii * sinPhi + A_ij * cosPhi;
			a[i * n + j] = a[j * n + i];
			a[j * n + j] = A_ji * sinPhi + A_jj * cosPhi;
		}
	}
}

//������������ ������� R_k �� (T_k,k+1)^T ������, ��� ����� ��� ��������� ���� � ��������� �� ����� �.�.
//��� ��������� ������, �������� ������ k � k+1 ������� ������� R_k, k = 1...n-1
void MULT(double *a, double *cosPhi, double *sinPhi, int n, int j){
    int i, k;
    double x, y;
    for(k = 1; k < j; ++k){
        for(i = 0; i < k + 1; i++){
            x = a[i * n + k - 1];
            y = a[i * n + k];
            a[i * n + k - 1] = x * cosPhi[k - 1] - y * sinPhi[k - 1];
            a[i * n + k] = x * sinPhi[k - 1] + y * cosPhi[k - 1];
        }
    }
}
//��� ��������� �� ������������ ������ (T_k,k+1)^T ������, ����� ��� ��������, �� � ��� ����������� ������ ������� Q_k �� k ����
void MULT1(double *a, double *cosPhi, double *sinPhi, int n){
    int i,k;
    double x,y;
    for(k = n-1; k > 0; k--){
        for(i = 0; i < n; i++){
            x = a[(k - 1) * n + i];
            y = a[k * n + i];
            a[(k-1) * n + i] = cosPhi[k-1] * x + sinPhi[k-1] * y;
            a[k * n + i] = -sinPhi[k-1] * x + cosPhi[k-1] * y;
        }
    }
}

