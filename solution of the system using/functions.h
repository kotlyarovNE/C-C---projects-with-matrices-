
int matr_file(int n, char *filename, double *A); //чтения матр из файла

void input_from_file(int n, double *A, double *a, double *b, int my_rank, int p);//заполн. матр. из файла

double f(int k, int n, int i, int j);//функия заполнения матр

void input(int k, int n, double *A, double *B, int my_rank, int p);//ввод матр по формуле

int max(int x, int y);



void inp_A_B(int n, int k, double *A, double *B);//Копия матрицы A, B, A - нужна для чтения из файла, B - для невязки надо

void input_b_after_solution(double *B, double *b, int n, int my_rank, int p, double *buf);//Тут копия вектора измененного для невязки

void output(double* a, int strok,int stolb, int n);//обычный самый вывод матр A




void OutputMatrix(int n, int m, double *a, double *b, double *x, int my_rank, int p);//Нормальный вывод матрицы A

void OutputVector(int m, double *b, double *x, int my_rank, int p);//Нормальный вывод вектора B




double Error(int n, double* a, double *B, double *vect, int my_rank, int p);//Невязка

double accur(int n, double* x, int my_rank);//разность между точным и прибл. решением

void NORM(int n, double *a, double *b, int p, int my_rank);




int solve(int n, double *a, int my_rank, int p, double *b, double *vect, int *buffer_Stolb, double *BUFFER);

void S(double *b, int *Index, double *buf, int my_rank, int p, int n, double *B);



