/*
*********************************************
*                                           *
*       Testing Welford West algorithm      *
*                                           *
*********************************************

Author: Alexandre A. A. Almeida
First Created: 20/03/2021

*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "./test_include_welford_west.h"


int main(int argc, char *argv[])
{

    double mean;
    double variance;
    double stddev;

    double value;
    int n=1;

    FILE * fr = fopen("input_welford_west_test.txt", "rt");
     
    char tempBuffer[100]; // temporary buffer with the line values
    
    // While not end of file
    while(!feof(fr))
    {
  
      // reads a line from the stream, save the number of characters in the line
        fgets(tempBuffer,100,fr);
        sscanf(tempBuffer, "%lf", &value);
        if (feof(fr)) { break; }
//        printf("Entrada \n");
//        printf("Valor: %lf \n", value);

        if (n==0)
        {
            mean = initialize_mean(value);
            variance = initialize_variance();
        }
        else
        {
            variance = update_variance(variance, mean, value, (double)n);
            mean = update_mean(mean, value, (double)n);
        }
        
//        printf("Mn = %lf. \n", mean);
//        printf("Sn = %lf. \n", variance);
        
        n++;
    }
    
    fclose(fr);

    variance = variance/((double)n-1.0);
    stddev = sqrt(variance);

    printf("A média final é de %.12lf. \n", mean);
    printf("A variância final é de %.12lf. \n", variance);
    printf("O desvio padrão final é de %.12lf. \n", stddev);

    return 0;
}
