/*
*********************************************
*                                           *
*       Simulation of Langevin Equation     *
*             for one particle              *
*                                           *
*********************************************

Author: Alexandre A. A. Almeida
Date: 07/03/2019

*/

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <math.h>



int stochastic_euler(char *name1, double x0, double v0, double A, double gamma, double tau, double tempo, unsigned int semente){

    /* Defining simulation parameters
     Idea 1: Defining mass(m), coeff (gamma), Temperature(T), Boltzmann(kB)
     Idea 2: Defining coeff (gamma) and amplitude (A)
     Relation: A = 2*gamma*kB*T/m

     Defining timestep (tau) in seconds
     Defining number of timesteps(N) for 1 minute (60 s)
    */
    //Declaring parameters
    long int N;
    N = (int)tempo/tau;


    /* Defining the Random Number Generator (RNG) */
    // Declaring RNG
    gsl_rng *r;
    // Using Mersenne Twister generator
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, semente);


    /* Defining output files */
    // Declaring files
    FILE *file;
    // Creating csv files
    file = fopen(name1, "w");
    // Creating column headers
    fprintf(file, "time,x,v\n");

    /*Defining time, space, velocity and random variables*/
    // Declaring time t, space x, velocity v, random eps
    long int t;
    double x, v;
    double eps;
    // Declarint initial space and velocity
    x = x0;
    v = v0;
    t = 0;

    /* Running the program */
    // Initial point
    fprintf(file, "%.8g,%.8g,%.8g\n",t*tau,x,v);
    eps = gsl_ran_ugaussian(r);

    // Next points
    for(t=1;t<N+1;t++){


        x = x + tau*v;
        v = v - gamma*tau*v + sqrt(tau*A)*eps;
        eps = gsl_ran_ugaussian(r);

        fprintf(file, "%.8g,%.8g,%.8g\n",t*tau,x,v);

    }

    // Closing files
    fclose(file);
    // Cleaning memory
    gsl_rng_free(r);

    return 0;
}

int stochastic_rk4(char *name1, double x0, double v0, double A, double gamma, double tau, double tempo, unsigned int semente){

    //Declaring parameters
    long int N;
    N = (int)tempo/tau;

    /* Defining the Random Number Generator (RNG) */
    // Declaring RNG
    gsl_rng *r;
    // Using Mersenne Twister generator
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, semente);


    /* Defining output files */
    // Declaring files
    FILE *file;
    // Creating csv files
    file = fopen(name1, "w");
    // Creating column headers
    fprintf(file, "time,x,v\n");

    /*Defining time, space, velocity and random variables*/
    // Declaring time t, space x, velocity v, random eps
    long int t;
    double v[5];
    double Kv[4];
    double x[5];
    double Kx[4];
    double eps;
    // Declarint initial space and velocity
    x[0] = x0;
    v[0] = v0;
    t = 0;

    /* Running the program */
    // Initial point
    fprintf(file, "%.8g,%.8g,%.8g\n",t*tau,x[0],v[0]);

    // Next points
    for(t=1;t<N+1;t++){

        eps = gsl_ran_ugaussian(r);

        Kv[0] = - gamma*v[0];
        Kx[0] = v[0];

        v[1] = v[0] + 0.5*tau*Kv[0] + 0.5*sqrt(A*tau)*eps;
        x[1] = x[0] + 0.5*tau*Kx[0];
        Kv[1] = - gamma*v[1];
        Kx[1] = v[1];

        v[2] = v[0] + 0.5*tau*Kv[1] + 0.5*sqrt(A*tau)*eps;
        x[2] = x[0] + 0.5*tau*Kx[1];
        Kv[2] = - gamma*v[2];
        Kx[2] = v[2];

        v[3] = v[0] + tau*Kv[2] + sqrt(A*tau)*eps;
        x[3] = x[0] + tau*Kx[2];
        Kv[3] = - gamma*v[3];
        Kx[3] = v[3];

        v[4] = v[0] + (Kv[0]+2.0*Kv[1]+2.0*Kv[2]+Kv[3])*tau/6.0 + sqrt(A*tau)*eps;
        x[4] = x[0] + (Kx[0]+2.0*Kx[1]+2.0*Kx[2]+Kx[3])*tau/6.0;

        fprintf(file, "%.8g,%.8g,%.8g\n",t*tau,x[4],v[4]);

        x[0] = x[4];
        v[0] = v[4];

    }

    // Closing files
    fclose(file);
    // Cleaning memory
    gsl_rng_free(r);

    return 0;
}

int stochastic_rk2(char *name1, double x0, double v0, double A, double gamma, double tau, double tempo, unsigned int semente){

    //Declaring parameters
    long int N;
    N = (int)tempo/tau;

    /* Defining the Random Number Generator (RNG) */
    // Declaring RNG
    gsl_rng *r;
    // Using Mersenne Twister generator
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, semente);


    /* Defining output files */
    // Declaring files
    FILE *file;
    // Creating csv files
    file = fopen(name1, "w");
    // Creating column headers
    fprintf(file, "time,x,v\n");

    /*Defining time, space, velocity and random variables*/
    // Declaring time t, space x, velocity v, random eps
    long int t;
    double v[3];
    double Kv[2];
    double x[3];
    double Kx[2];
    double eps;
    // Declarint initial space and velocity
    x[0] = x0;
    v[0] = v0;
    t = 0;

    /* Running the program */
    // Initial point
    fprintf(file, "%.8g,%.8g,%.8g\n",t*tau,x[0],v[0]);

    // Next points
    for(t=1;t<N+1;t++){

        eps = gsl_ran_ugaussian(r);

        Kv[0] = - gamma*v[0];
        Kx[0] = v[0];

        v[1] = v[0] + tau*Kv[0] + sqrt(A*tau)*eps;
        x[1] = x[0] + tau*Kx[0];
        Kv[1] = - gamma*v[1];
        Kx[1] = v[1];

        v[2] = v[0] + 0.5*tau*(Kv[0]+Kv[1]) + sqrt(A*tau)*eps;
        x[2] = x[0] + 0.5*tau*(Kx[0]+Kx[1]);

        fprintf(file, "%.8g,%.8g,%.8g\n",t*tau,x[2],v[2]);

        x[0] = x[2];
        v[0] = v[2];

    }

    // Closing files
    fclose(file);
    // Cleaning memory
    gsl_rng_free(r);

    return 0;
}
