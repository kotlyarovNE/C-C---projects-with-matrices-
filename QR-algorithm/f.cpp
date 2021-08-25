#include "function.h"
#include<math.h>
#include <stdlib.h>
#include <stdio.h>
#include<iostream>
using namespace ::std;

double f(int k,int n, int i, int j)
{
    if(k==1) return (double(n-max(i,j)+1));

    if(k==2) {
        if(i == j){return 2;}
        if(abs(i - j) == 1){return -1;}
        else{return 0;}
    }

    if(k==3){
        if((i == j) && (i < n)){return 1;}
        if(j == n){return i;}
        if(i == n){return j;}
        else{return 0;}
    }

    if(k==4) return (1/double((i+j-1)));

    else return 0;
return 0;
}
