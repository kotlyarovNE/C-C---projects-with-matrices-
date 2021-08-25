#include "function.h"
#include <stdio.h>
#include<math.h>
#include<stdlib.h>

int SolveSystem(int n, double* a, double* b, double* x, int* Index, double *S)
{
	int i;
	int j;
	int k;
	int indMax1; // макс по строке
	int indMax2; //макс по столбцу
	double tmp;
	double max; //макс в подматрице

	Norm(a, b, n);

	for (i = 0; i < n; ++i)
		Index[i] = i;

	for (i = 0; i < n; ++i)
	{
		max = fabs(a[i * n + i]);
		indMax1 = i;
		indMax2 = i;

		for (j = i; j < n; ++j)
			for (k = i; k < n; ++k)
				if (max < fabs(a[j * n + k]))
				{
					max = fabs(a[j * n + k]);
					indMax1 = j;
					indMax2 = k;
				}

		k = Index[i];
		Index[i] = Index[indMax2];  // индекс i теперь номер столбца с макс элеменом
		Index[indMax2] = k;   //номер столбца с макс элементом теперь i

		for (j = 0; j < n; ++j)   //мен€ю местами строчки i и с макс элементом
		{
			tmp = a[i * n + j];
			a[i * n + j] = a[indMax1 * n + j];
			a[indMax1 * n + j] = tmp;
		}

		tmp = b[i];
		b[i] = b[indMax1];  //мен€ю и сам вектор b в соотвествии со строкой с макс элементом
		b[indMax1] = tmp;

		for (j = 0; j < n; ++j)    //мен€ю местами столбец i и с макс элементом
		{
			tmp = a[j * n + i];
			a[j * n + i] = a[j * n + indMax2];
			a[j * n + indMax2] = tmp;
		}

		if (fabs(a[i * n + i]) < 1e-13)  // если на i-ом шаге наш макс элемент в подматрице очень мал, то матрица вырождена и алгорит останавливаетс€
			return -1;

		tmp = 1.0 / a[i * n + i];
		for (j = i; j < n; ++j)           //делим на макс элемент строку ну и вектор тоже
			a[i * n + j] *= tmp;
		b[i] *= tmp;

		for (j = i + 1; j < n; ++j)     //из строк, лежащих ниже главной нашей в подматрице, мы вычитаем главную, умноженную на первый элемент этих строк (обнул€ем подстолец)
		{
			tmp = a[j * n + i];
			for (k = i; k < n; ++k)
				a[j * n + k] -= a[i * n + k] * tmp;
			b[j] -= b[i] * tmp;         //и вектор не забываем преобразовать также
		}
	}

	for (i = n - 1; i >= 0; --i)        //решать систему уже легко, решаем с конца, мы мен€ли местами столбцы, это нехорошо, но мы запомнили как мен€ли их в index[i]
	{
		tmp = b[i];
		for (j = i + 1; j < n; ++j)
			tmp -= a[i * n + j] * x[Index[j]];
		x[Index[i]] = tmp;
	}

	for( i = 0; i < n; i++)
        S[i] = x[Index[i]];

	return 0;
}
