#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include "functions.h"

int main(int argc, char **argv)
{
	double t; // B - копи€ ј, чтобы в конце посчитать норму нев€зки, t - астро. врем€
	int rezult, n, m, k, p, error = 0, err;		// p - общее число процессов, curr - номер текущего процесса; пор€док задани€ аргументов в командной строке такой же, но p - уже не параметр командной строки
	int max_rows, my_rank, *Index, determ[1];
	double *a, *b, *x, *A, *B, *vect, *buf;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	if ( (argc != 5) && (argc != 4) )  {
		printf("Error! Incorrecr number of parameters\n");
		error = 1;
	}

	MPI_Allreduce(&error, &err, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if (err == 1) {
		if (my_rank == 0) printf("Error! Incorrect parameters!\n");
		MPI_Finalize();
		return -1;
	}

	n = atoi(argv[1]);
	m = atoi(argv[2]);
	k = atoi(argv[3]);


	if ( ( (k == 0) && (argc != 5) ) || ( (k != 0) && (argc != 4) ) || (p < 0) || (p > n)) {
		error = 1;
	}

	MPI_Allreduce(&error, &err, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if (err == 1) {
		if (my_rank == 0) printf("Error! Incorrect parameters!\n");
		MPI_Finalize();
		return -1;
	}

	if (p > n) {p = n;}

	max_rows = n/p + n%p;   //максимум на 1 процесс строк


	a = (double*)malloc(max_rows*n*sizeof(double)); // т.к будем пересылать между процессами, выдел€ем максимум, как в Ѕогачеве
	if (a==NULL) error = 1;

    b = (double*)malloc(max_rows*sizeof(double));
	if (b==NULL) error = 1;

    x = (double*)malloc((n + 1)*sizeof(double));
	if (x==NULL) error = 1;

    A = (double*)malloc(n*n*sizeof(double)); //ƒл€ чтени€ из файла
	if (A==NULL) error = 1;

	B = (double*)malloc(n*sizeof(double)); //дл€ подстчета нормы нев€зки, нам нужен старый вектор
	if (B==NULL) error = 1;

	vect = (double*)malloc(max_rows*sizeof(double));
	if (vect==NULL) error = 1;

	buf = (double*)malloc(n*sizeof(double));
	if (buf==NULL) error = 1;

	Index = (int *)malloc(n * sizeof(int));
	if (Index==NULL)
		error = 1;

    MPI_Allreduce(&error, &err, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); // собираем ошибки
	if (err == 1) {
        free(a);
		free(b);
		free(x);
		free(A);
        free(B);
        free(vect);
        free(buf);
        free(Index);
		if (my_rank == 0) printf("Error! Not enough memory! \n");
		MPI_Finalize();
		return -1;
	}

	//“еперь заполнение матрицы:
	if (k!=0) {
		input(k, n, a, b, my_rank, p);
		if(my_rank == 0){
            inp_A_B(n, k, A, B);}
	}
	else {
		rezult = matr_file(n, argv[4], A);
		if (rezult == -1) {
			if (my_rank == 0) printf("Error! The file cannot be opened! \n");
			free(a);
			free(b);
			free(x);
			free(A);
			free(B);
			free(vect);
            free(buf);
            free(Index);
			MPI_Finalize();
			return -1;
		}
		if (rezult == -2) {
			if (my_rank == 0) printf("Error! The file doesn't meet the requirements! \n");
			free(a);
			free(b);
			free(x);
			free(A);
			free(B);
			free(vect);
            free(buf);
            free(Index);
			MPI_Finalize();
			return -1;
		}

		input_from_file(n, A, a, b, my_rank, p);
	}

    OutputMatrix(n, m, a, b, x, my_rank, p);
    NORM(n, a, b, p, my_rank);
	if (my_rank == 0) printf("Matrix A:\n\n");
	OutputMatrix(n, m, a, b, x, my_rank, p);

	MPI_Barrier(MPI_COMM_WORLD);
	t = MPI_Wtime();

    //NORM(n, a, b, p, my_rank);

	determ[0] = solve(n, a, my_rank, p, b, vect, Index, x);


	MPI_Barrier(MPI_COMM_WORLD);
	t = MPI_Wtime() - t;


    if(determ [0] == -1){
        if (my_rank == 0) printf("  CAN'T SOLVE - MATRIX DEGENERATE!\n");
    }

    else{

        if (my_rank == 0) printf("\nA:\n");
        OutputMatrix(n, m, a, vect, x, my_rank, p);

        //¬ывод решени€
        if (my_rank == 0) printf("\nSolution:\n");

        S(b, Index, buf, my_rank, p, n, B);//ѕерестановка обратна€
        MPI_Barrier(MPI_COMM_WORLD);

        if (my_rank == 0){output(B, 1, m, n);}

        if (my_rank == 0) printf("\n\nSolution time = %e\n", t);


        if(my_rank == 0){
            printf("\nSolution Accuracy: ");}

        MPI_Barrier(MPI_COMM_WORLD);
        accur(n, B, my_rank); //разница между точным и прибл. решением


        if (my_rank == 0) {printf("\n||Ax-b||*||b||^-1: ");}
        input_b_after_solution(B, b, n, my_rank, p, buf); //заполнение матрицы B дл€ нормы нев€зки

        MPI_Barrier(MPI_COMM_WORLD);
        Error(n, a, B, vect, my_rank, p); //норма нев€зки

    }


	free(a);
	free(b);
	free(x);
	free(A);
    free(B);
    free(vect);
    free(buf);
    free(Index);
	MPI_Finalize();

	return 0;
}
