#include "function.h"
#include<math.h>
#include <stdlib.h>
#include <stdio.h>
#include<iostream>
using namespace ::std;

double f(int k,int n, int i, int j)
{
    if(k==1) return (double(n-max(i,j)+1));
    if(k==2) return (double(max(i,j)));
    if(k==3) return (double(abs(i-j)));
    if(k==4) return (1/double((i+j-1)));
    else return 0;
return 0;
}
