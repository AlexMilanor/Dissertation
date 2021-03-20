
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Teste
double 
initialize_variance()
{
    /* 
    Welford-West numerical schema for the mean and variance computation
    "Algorithms for Computing the Sample Variance: Analysis and Recommendations", Chan et al., 1983
    */
    return 0.0;
}

double 
update_variance(double prev_var, double prev_mean, double point, double n_th)
{
    /* 
    Welford-West numerical schema for the mean and variance computation
    "Algorithms for Computing the Sample Variance: Analysis and Recommendations", Chan et al., 1983
    */
    double new_var;
    new_var = prev_var + (n_th-1.0)*(point - prev_mean)*((point-prev_mean)/n_th);
    return new_var;
}


double 
initialize_mean(double first_point)
{
    /* 
    Welford-West numerical schema for the mean and variance computation
    "Algorithms for Computing the Sample Variance: Analysis and Recommendations", Chan et al., 1983
    */
    return first_point;
}


double 
update_mean(double prev_mean, double point, double n_th)
{
    /* 
    Welford-West numerical schema for the mean and variance computation
    "Algorithms for Computing the Sample Variance: Analysis and Recommendations", Chan et al., 1983
    */
    double new_mean;
    new_mean = prev_mean + (1.0/n_th)*(point-prev_mean);
    return new_mean;
}


int 
take_sqrt_of_vector(double *vector, int length)
{
    /*
     Takes the square root of each element in a (double) array 
     */
    for(int i=0;i<length;i++){
        vector[i] = sqrt(vector[i]);
    }
}


int 
divide_vector_by_value(double *vector, int length, double value)
{
    /*
    Divides each element in a (double) array by a given value 
    */
    for(int i=0;i<length;i++){
        vector[i] = vector[i]/value;
    }
}


int 
transform_sqrSum_in_variance(double *sqrSum_vector, int length, double number_of_points)
{
    /*
    Transform the sum of squares from the Welford-West algorithm into the variance 
    */
    divide_vector_by_value(sqrSum_vector, length, number_of_points-1.0);
}


int 
transform_variance_in_std_dev(double *var_vector, int length)
{
    /*
    Transform the variance into the standar deviation 
    */
    take_sqrt_of_vector(var_vector, length);
}


