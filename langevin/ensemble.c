/*
*********************************************
*                                           *
*       Simulation of Langevin Equation     *
*             for ensemble of               *
*                particles                  *
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
#include <string.h>
#include <omp.h>
#include "particle.h"
#include "rand_seed.h"

int main(void){

    /* Defining simulation parameters
     Number of particles (N_p);
    */
    //Declaring parameters
    const int N_p = 1000;

    /* Defining output files names
     Name of file
     datatype (dtpe)
     number of particle (number)(n)
     */
    // Declaring strings
    char s_euler[80];
    char s_rk2[80];
    char s_rk4[80];
    char dtpe[10];
    char number1[10];
    char number2[10];
    char number3[10];
    int n1, n2, n3;

    // Defining plot type
    strcpy(dtpe, ".csv");

    // Simulation parameters
    double A = 0.1; //Noise amplitude
    double gamma = 0.1; //Drag Coefficient
    double tau = 0.01; //Timestep
    double tempo = 60.0; //Simulation Time in seconds

    // Defining seeds for the simulations
    unsigned int semente1, semente2, semente3;

    //Initial values
    double x0=0.0;
    double v0=1.0;

//    #pragma omp parallel
//    {
//    #pragma omp sections
//    {
//    #pragma omp section
    for(n1=1;n1<=N_p;n1++){

        strcpy(s_euler, "data/euler_particle");
        sprintf(number1, "%d", n1);
        strcat(s_euler, number1);
        strcat(s_euler, dtpe);

        /* Using urandom as seed
        semente1 = rand_seed();
        */
        semente1 = n1;

        stochastic_euler(s_euler, x0, v0, A, gamma, tau, tempo, semente1);

    }

//    #pragma omp section
    for(n2=1;n2<=N_p;n2++){

        strcpy(s_rk2, "data/srk2_particle");
        sprintf(number2, "%d", n2);
        strcat(s_rk2, number2);
        strcat(s_rk2, dtpe);

        /* Using urandom as seed
        semente2 = rand_seed();
        */
        semente2 = n2;

        stochastic_rk2(s_rk2, x0, v0, A, gamma, tau, tempo, semente2);

    }

//    #pragma omp section
    for(n3=1;n3<=N_p;n3++){

        strcpy(s_rk4, "data/srk4_particle");
        sprintf(number3, "%d", n3);
        strcat(s_rk4, number3);
        strcat(s_rk4, dtpe);

        /* Using urandom as seed
        semente3 = rand_seed();
        */
        semente3 = n3;

        stochastic_rk4(s_rk4, x0, v0, A, gamma, tau, tempo, semente3);

    }
//    }
//    }

    return 0;
}
