/*

Code to simulate Stochastic Runge Kutta of fourth order.
Used specifically for Harmonic Chain
(maybe with anharmonic potential).

The baths are of Langevin type.

Author: Alexandre A. A. Almeida
Date: 05/06/2019

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <math.h>

// Pi decimal expansion taken from the Online Encyclopedia of Integer Sequences
// Source: http://oeis.org/A000796
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279
#endif /* M_PI */

// 1/Pi decimal expansion taken from the Online Encyclopedia of Integer Sequences
// Source: http://oeis.org/A049541
#ifndef M_1_PI
#define M_1_PI 0.318309886183790671537767526745
#endif /* M_1_PI */

double
no_potential(double x, double A)
{
    /* Force assuming no external potential */
    return 0.0;
}


double 
cosine_potential(double x, double V)
{
    /* Force due to a cosine potential */
    double a0=1.0; // external potential length
    return -V*sin(2.0*M_PI*x);
}


double 
phi_4_potential(double x, double k)
{
    /* Force due to a quartic (phi4) potential */
    return -k*pow(x,3);
}


double (*potential)(double x, double A);


// https://stackoverflow.com/questions/840501/how-do-function-pointers-in-c-work
typedef double (*def_potential)(double,double);

def_potential 
assign_potential(char *potential_name) 
{    
    def_potential functionPtr;

    if (strcmp(potential_name,"none")==0){
        functionPtr = no_potential;
    }

    else if (strcmp(potential_name,"cosine")==0){
        functionPtr = cosine_potential;
    }

    else if (strcmp(potential_name,"phi4")==0){
        functionPtr = phi_4_potential;
    }

    else {
        functionPtr = no_potential;
    }

    return functionPtr;
}


double 
force_bath(double gamma, double k, double x0, double x1, double v0, double V)
{
    /* Defining the deterministic forces on the Langevin Bath */
    return -gamma*v0 - k*(x0 - x1) - k*x0 + potential(x0, V);
}


double 
spring_force(double k, double x1, double x2)
{
    /* Force of interaction between neighboring particles */
    return -k*(x1 - x2);
}


double 
force(double k1, double k2, double x0, double x1, double x2, double V)
{
    /* Defining the deterministic forces on x1 */
    return spring_force(k1, x1, x0) + spring_force(k2, x1, x2) + potential(x1, V);
}

double
sign(double x)
{
    if (x >= 0.0)
    {
        return 1.0;
    }
    else 
    {
        return -1.0;
    }
}

double 
inter_force(double k, double mu, double x1, double x2)
{
    /* Force in the interphase of the left and right chain */
    double force_interf;
    force_interf = -k*sign(x1-x2)*pow(fabs(x1 - x2),mu-1.0);
    return force_interf;
}


double 
dxdt(double v)
{
    /* Derivative of x for each particle */
    return v;
}


double 
x_m(double delta, double x0, double K0, int k)
{
    /* Defining the variables for the next step on SRK4 routine */

    double x;

    if(k < 3){
        x = x0 + 0.5*delta*K0;
    }

    else{
        x = x0 + delta*K0;
    }

    return x;
}


double 
step_next(double delta, double x0, double K0, double K1, double K2, double K3)
{
    /* Calculating x(t+delta(t)) for SRK4 */
    return x0 + (K0+2.0*K1+2.0*K2+K3)*delta/6.0;
}


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
    // new_mean = prev_mean + (1.0/n_th)*(point-prev_mean);
    new_mean = prev_mean + point;
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


int 
write_vector_to_row(char *name, double *vector, int length, FILE *write_file)
{
    /*
    Write all the elements of an array to the entire row of a file. 
    */

    fprintf(write_file, "%s",name);

    for(int i=0;i<length;i++){
        fprintf(write_file, ",%.12g",vector[i]);
    }

    fprintf(write_file, "\n");
}



int 
diode_srk4(char *filename, char *filename_vars, char *potential_name, 
           double *x0, double *v0, int N_C, double gama, 
           double *V, double mu, double tau, double tempo, double TransientTime, 
           double *A, double *k, unsigned int *semente)
{

    /*
    Main code for running the 4th order Stochastic Runge Kutta (SRK4).
    ----------------------------------------
    char *filename: filename of the file
    char *filename_vars: filename of the file with dynamic variables
    double *x0: initial conditions for position
    double *v0: initial conditions for velocity
    int N_C: number of particles in the chain
    double gama: drag coefficient
    double *V: array of amplitudes of the external potential (left chain and right chain)
    double mu: exponent of the potential in the middle
    double tau: timestep
    double tempo: time interval
    double TransientTime: time of the transient
    double *A: Array of amplitudes of the heat baths (left heat bath and right heat bath)
    double *k: Array of spring constants (left chain and right chain)
    unsigned int *semente: random number generator seeds (left heat bath and right heat bath)
    */

    /* Main code for the SRK4 in the harmonic chain */

    /* Defining output files */
    FILE *file; // Declaring file for summaries
    FILE *file_vars; // Declaring file for variables

    file = fopen(filename, "w"); // Creating csv files
    file_vars = fopen(filename_vars, "w"); // Creating csv files 

    /* Creating column headers */
    fprintf(file, "index");
    fprintf(file_vars, "t");
    
    for(int i=0;i<N_C;i++)
    {
        fprintf(file, ",%d",i+1);
        fprintf(file_vars, ",x%d,v%d", i+1, i+1);
    }
    
    fprintf(file, "\n");
    fprintf(file_vars, "\n");


    /* Defining the external potential */
    potential = assign_potential(potential_name);


    /* Defining the Random Number Generators (RNGs) */
    gsl_rng *re;
    gsl_rng *rd;

    // Using Mersenne Twister generator
    re = gsl_rng_alloc(gsl_rng_mt19937);
    rd = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(re, semente[0]);
    gsl_rng_set(rd, semente[1]);

    /* Defining main variables */
    double epse, epsd; //Noise of heat beaths
    double x[N_C][5]; //Position of particles
    double v[N_C][5]; //Velocity of particles
    double Kx[N_C][4]; //Matrix K (of Stochastic Runge Kutta)
    double Kv[N_C][4]; //Matrix K (of Stochastic Runge Kutta)
    double Temp[N_C]; //Kinetic temperature of the particles
    double J_flux[N_C]; //Heat flux between particles



    /* Saving initial Conditions*/
    fprintf(file, "x_t_0");

    for(int i=0;i<N_C;i++)
    {
        x[i][0]=x0[i];
        fprintf(file, ",%.12g", x[i][0]);
    }
    
    fprintf(file, "\n");

    fprintf(file, "v_t_0");
    
    for(int i=0;i<N_C;i++)
    {
        v[i][0]=v0[i];
        fprintf(file, ",%.12g", v[i][0]);
    }
    
    fprintf(file, "\n");

    fprintf(file_vars, "0.0");
    for(int i=0;i<N_C;i++)
    {
        fprintf(file_vars, ",%.8g,%.8g", x[i][0], v[i][0]);
    }
    fprintf(file_vars, "\n");


    /* Vector of outputs */
    double x_mean[N_C]; // time mean of position
    double v_mean[N_C]; // time mean of velocity
    double T_mean[N_C]; // time mean of kinetical Temperature
    double J_mean[N_C]; // time mean of heat flux

    double x_std[N_C]; // time standard deviation of position
    double v_std[N_C]; // time standard deviation of velocity
    double T_std[N_C]; // time standard deviation of kinetical Temperature
    double J_std[N_C]; // time standard deviation of heat flux

    /* Initializing Outputs*/
    int initialized_outputs[N_C];

    for (int i=0;i<N_C;i++)
    {
        initialized_outputs[i]=0;
    }


    /* Running the program */
    // Iteration for the time interval
    long int time_points; // Number of points of time
    time_points = (long int)(tempo/tau); // time_points = time_range/tau
    long int t; // Time
    long int transient_point = (long int)(TransientTime/tau + 0.5);


    for(t=1;t<time_points+1;t++)
    {

        //Noise of the Heat Baths
        epse = gsl_ran_ugaussian(re);
        epsd = gsl_ran_ugaussian(rd);

        /* ---------------------- Three first steps of SRK4 ---------------------*/
        for(int m=1;m<4;m++)
        {

            //Left Heat Bath
            Kv[0][m-1] = force_bath(gama, k[0], x[0][m-1], x[1][m-1], v[0][m-1], V[0]);
            Kx[0][m-1] = dxdt(v[0][m-1]);
            v[0][m] = x_m(tau, v[0][0], Kv[0][m-1], m) + 0.5*sqrt(A[0]*tau)*epse;
            x[0][m] = x_m(tau, x[0][0], Kx[0][m-1], m);

            //Middle Particles
            for(int i=1;i<N_C-1;i++)
            {
                if(i < ((N_C/2)-1))
                {
                    Kv[i][m-1] = force(k[0],k[0], x[i-1][m-1],x[i][m-1],x[i+1][m-1], V[0]);
                    Kx[i][m-1] = dxdt(v[i][m-1]);
                    v[i][m] = x_m(tau, v[i][0], Kv[i][m-1], m);
                    x[i][m] = x_m(tau, x[i][0], Kx[i][m-1], m);
                }
                else if(i == N_C/2 - 1)
                {
                    Kv[i][m-1] = spring_force(k[0], x[i][m-1],x[i-1][m-1]) + inter_force(k[1], mu, x[i][m-1], x[i+1][m-1]) + potential(x[i][m-1], V[0]);
                    Kx[i][m-1] = dxdt(v[i][m-1]);
                    v[i][m] = x_m(tau, v[i][0], Kv[i][m-1], m);
                    x[i][m] = x_m(tau, x[i][0], Kx[i][m-1], m);
                }
                else if(i == N_C/2)
                {
                    Kv[i][m-1] = spring_force(k[2], x[i][m-1],x[i+1][m-1]) + inter_force(k[1], mu, x[i][m-1], x[i-1][m-1]) + potential(x[i][m-1], V[1]);
                    Kx[i][m-1] = dxdt(v[i][m-1]);
                    v[i][m] = x_m(tau, v[i][0], Kv[i][m-1], m);
                    x[i][m] = x_m(tau, x[i][0], Kx[i][m-1], m);
                }
                else if(i > N_C/2)
                {
                    Kv[i][m-1] = force(k[2],k[2], x[i-1][m-1],x[i][m-1],x[i+1][m-1], V[1]);
                    Kx[i][m-1] = dxdt(v[i][m-1]);
                    v[i][m] = x_m(tau, v[i][0], Kv[i][m-1], m);
                    x[i][m] = x_m(tau, x[i][0], Kx[i][m-1], m);
                }
            }

            //Right Heat Bath
            Kv[N_C-1][m-1] = force_bath(gama, k[2], x[N_C-1][m-1], x[N_C-2][m-1], v[N_C-1][m-1], V[1]);
            Kx[N_C-1][m-1] = dxdt(v[N_C-1][m-1]);
            v[N_C-1][m] = x_m(tau, v[N_C-1][0], Kv[N_C-1][m-1], m) + 0.5*sqrt(A[1]*tau)*epsd;
            x[N_C-1][m] = x_m(tau, x[N_C-1][0], Kx[N_C-1][m-1], m);
        }
        /* -------------------------------------------------------------------- */

        //x(t+d) = x(t) + (1/3)*(0.5*K_0 + K_1 + K_2 + 0.5*K_3)*d + (1/3)*(0.5*M_0z + M_1z + M_2z + 0.5*M_3z) 


        /* ------------ Fourth Step (Getting x(t+delta)) ------------------------*/

        //Left Heat Bath
        Kv[0][3] = force_bath(gama, k[0], x[0][3], x[1][3], v[0][3], V[0]);
        Kx[0][3] = dxdt(v[0][3]);
        v[0][4] = step_next(tau, v[0][0], Kv[0][0], Kv[0][1], Kv[0][2], Kv[0][3]) + sqrt(A[0]*tau)*epse;
        x[0][4] = step_next(tau, x[0][0], Kx[0][0], Kx[0][1], Kx[0][2], Kx[0][3]);

        //Middle Particles
        for(int i=1;i<N_C-1;i++)
        {
            if(i < N_C/2 - 1)
            {
                Kv[i][3] = force(k[0], k[0], x[i-1][3], x[i][3], x[i+1][3], V[0]);
                Kx[i][3] = dxdt(v[i][3]);
                v[i][4] = step_next(tau, v[i][0], Kv[i][0], Kv[i][1], Kv[i][2], Kv[i][3]);
                x[i][4] = step_next(tau, x[i][0], Kx[i][0], Kx[i][1], Kx[i][2], Kx[i][3]);
            }
            else if(i == N_C/2 - 1)
            {
                Kv[i][3] = spring_force(k[0], x[i][3],x[i-1][3]) + inter_force(k[1], mu, x[i][3], x[i+1][3]) + potential(x[i][3], V[0]);
                Kx[i][3] = dxdt(v[i][3]);
                v[i][4] = step_next(tau, v[i][0], Kv[i][0], Kv[i][1], Kv[i][2], Kv[i][3]);
                x[i][4] = step_next(tau, x[i][0], Kx[i][0], Kx[i][1], Kx[i][2], Kx[i][3]);
            }
            else if(i == N_C/2)
            {
                Kv[i][3] = spring_force(k[2], x[i][3],x[i+1][3]) + inter_force(k[1], mu, x[i][3], x[i-1][3]) + potential(x[i][3], V[1]);
                Kx[i][3] = dxdt(v[i][3]);
                v[i][4] = step_next(tau, v[i][0], Kv[i][0], Kv[i][1], Kv[i][2], Kv[i][3]);
                x[i][4] = step_next(tau, x[i][0], Kx[i][0], Kx[i][1], Kx[i][2], Kx[i][3]);
            }
            else if(i > N_C/2)
            {
                Kv[i][3] = force(k[2], k[2], x[i-1][3], x[i][3], x[i+1][3], V[1]);
                Kx[i][3] = dxdt(v[i][3]);
                v[i][4] = step_next(tau, v[i][0], Kv[i][0], Kv[i][1], Kv[i][2], Kv[i][3]);
                x[i][4] = step_next(tau, x[i][0], Kx[i][0], Kx[i][1], Kx[i][2], Kx[i][3]);
            }

        }

        //Right Heat Bath
        Kv[N_C-1][3] = force_bath(gama, k[2], x[N_C-1][3], x[N_C-2][3], v[N_C-1][3], V[1]);
        Kx[N_C-1][3] = dxdt(v[N_C-1][3]);
        v[N_C-1][4] = step_next(tau, v[N_C-1][0], Kv[N_C-1][0], Kv[N_C-1][1], Kv[N_C-1][2], Kv[N_C-1][3]) + sqrt(A[1]*tau)*epsd;
        x[N_C-1][4] = step_next(tau, x[N_C-1][0], Kx[N_C-1][0], Kx[N_C-1][1], Kx[N_C-1][2], Kx[N_C-1][3]);


        /*---------------------------------------------------------------------------*/

        //Next iteration
        for(int i=0;i<N_C;i++)
        {
            x[i][0]=x[i][4];
            v[i][0]=v[i][4];
        }

        // We save the results after the transient time
        if (t>transient_point-1)
        {

            // Saving time to file_vars
            fprintf(file_vars, "%.8g", (double)(t*tau));


            for (int i=0;i<N_C;i++)
            {

                // Saving variables to file_vars
                fprintf(file_vars, ",%.8g,%.8g", x[i][4], v[i][4]);

                // Calculating thermodynamic variables
                Temp[i] = pow(v[i][0],2.0);

                if (i == 0)
                {
                    J_flux[i]=0.0;
                }
                else if (i <= N_C/2 - 1)
                {
                    J_flux[i] = 0.5*(v[i][4]+v[i-1][4])*spring_force(k[0], x[i-1][4], x[i][4]);
                }
                else if(i == N_C/2)
                {
                    J_flux[i] = 0.5*(v[i][4]+v[i-1][4])*inter_force(k[1], mu, x[i-1][4], x[i][4]);
                }
                else if(i > N_C/2)
                {
                    J_flux[i] = 0.5*(v[i][4]+v[i-1][4])*spring_force(k[2], x[i-1][4], x[i][4]);
                }

                if (!initialized_outputs[i])
                {
                    x_std[i]=initialize_variance();
                    v_std[i]=initialize_variance();
                    T_std[i]=initialize_variance();
                    J_std[i]=initialize_variance();

                    x_mean[i]=initialize_mean(x[i][4]);
                    v_mean[i]=initialize_mean(v[i][4]);
                    T_mean[i]=initialize_mean(Temp[i]);
                    J_mean[i]=initialize_mean(J_flux[i]);

                    initialized_outputs[i]=1;
                }
                else 
                {
                    x_std[i] = update_variance(x_std[i], x_mean[i], x[i][4], (double)t);
                    v_std[i] = update_variance(v_std[i], v_mean[i], v[i][4], (double)t);
                    T_std[i] = update_variance(T_std[i], T_mean[i], Temp[i], (double)t);
                    J_std[i] = update_variance(J_std[i], J_mean[i], J_flux[i], (double)t);

                    x_mean[i] = update_mean(x_mean[i],x[i][4],(double)t);
                    v_mean[i] = update_mean(v_mean[i],v[i][4],(double)t);
                    T_mean[i] = update_mean(T_mean[i],Temp[i],(double)t);
                    J_mean[i] = update_mean(J_mean[i],J_flux[i],(double)t);
                }

            }

        }
    }

    // Now we save the simulation results
    divide_vector_by_value(x_mean, N_C, time_points - transient_point);
    divide_vector_by_value(v_mean, N_C, time_points - transient_point);
    divide_vector_by_value(T_mean, N_C, time_points - transient_point);
    divide_vector_by_value(J_mean, N_C, time_points - transient_point);

    write_vector_to_row("x_mean", x_mean, N_C, file);
    write_vector_to_row("v_mean", v_mean, N_C, file);
    write_vector_to_row("T_mean", T_mean, N_C, file);
    write_vector_to_row("J_mean", J_mean, N_C, file);

    transform_sqrSum_in_variance(x_std, N_C, time_points - transient_point);
    transform_sqrSum_in_variance(v_std, N_C, time_points - transient_point);
    transform_sqrSum_in_variance(T_std, N_C, time_points - transient_point);
    transform_sqrSum_in_variance(J_std, N_C, time_points - transient_point);

    transform_variance_in_std_dev(x_std, N_C);
    transform_variance_in_std_dev(v_std, N_C);
    transform_variance_in_std_dev(T_std, N_C);
    transform_variance_in_std_dev(J_std, N_C);

    write_vector_to_row("x_std", x_std, N_C, file);
    write_vector_to_row("v_std", v_std, N_C, file);
    write_vector_to_row("T_std", T_std, N_C, file);
    write_vector_to_row("J_std", J_std, N_C, file);


    /* Cleaning memory */
    fclose(file); // Closing files
    fclose(file_vars); //Closing files
    gsl_rng_free(re); // Closing left RNG
    gsl_rng_free(rd); // Closing right RNG



    return 0;
}
