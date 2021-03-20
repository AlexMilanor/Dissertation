#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

double 
cosine(double x)
{
    return cos(x);
}

double 
squared(double x)
{
    return pow(x,2.0); 
}

double (*function)(double x);


typedef double (*def_potential)(double);
// note that the typedef name is indeed myFuncDef

def_potential assign_potential(char *function_name) {
    
    def_potential functionPtr;

    if (strcmp(function_name, "cosine")==0){
        functionPtr = &cosine;
    }
    else if (strcmp(function_name, "squared")==0){
        functionPtr = &squared;
    }

    return functionPtr;
}


double
calculate(double x)
{
    return function(x);
}


void main(int argc, char *argv[]){

    function = assign_potential(argv[1]);
    

    double a = 5.0;
    double y;

    y = calculate(a);

    printf("Valor de y: %.15lf\n",y);

}