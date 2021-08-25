#include "function.h"
#include <stdio.h>
#include<math.h>
#include <stdlib.h>

void output(double* a, int m, int n){
    int i,j;

    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
          printf("%10.3e ",a[i*n+j]);
        }
        printf("\n");
    }
}
