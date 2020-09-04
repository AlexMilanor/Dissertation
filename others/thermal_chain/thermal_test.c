/*
*********************************************
*                                           *
*       Simulation of Fourier Chain         *
*             for sample                    *
*                                           *
*********************************************

Author: Alexandre A. A. Almeida
Date: 01/04/2019

*/

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <math.h>
#include <string.h>
#include "diode.h"

int main(int argc, char **argv){

/* =========================  Define control variables =========================== */
    /* Auxiliary Variables*/
    int n, l; // Control variables

    /*Testing if number of input arguments is correct*/
    if (argc < 11){
        printf("Too few arguments.\n");
        return 0;
    }
    if (argc>11){
        printf("Too many arguments.\n");
        return 0;
    }

/* ========================= Inicializando Parâmetros ============================= */
    /*Testing the number of particles inputted*/
    char *stop1;
    int arg1 = strtol(argv[1],&stop1,10);
    if (*stop1!='\0'||stop1==argv[1]){
        printf("Number of particles must be an integer.\n");
        return 0;
    }

    /*Testing the chain parameters inputted*/
    char *stop2_9[8];
    double args2_9[8];
    for(n=2;n<=9;n++){
        args2_9[n-2] = strtod(argv[n],&stop2_9[n-2]);
	if (*stop2_9[n-2]!='\0'||stop2_9[n-2]==argv[n]){
        	printf("Chain parameters must be double.\n");
        	return 0;
		}
	}

	/* Saving the filename inputted*/
    char *data_dir=argv[10];

    /* Defining output files names */
    char plot_srk4[3][80]; // Name of file
    char dtpe[3][10]; // Data Type
    char number3[3][10]; // Number of simulation

    /* Defining plot type */
    strcpy(dtpe[0], ".csv"); // Thread_ID=0
    strcpy(dtpe[1], ".csv"); // Thread_ID=1
    strcpy(dtpe[2], ".csv"); // Thread_ID=2

    /* Defining simulation parameters */
    const int number_of_samples=120; 
    const int number_of_particles=arg1; // (argv[1]) (number_of_particles=50)
    const double timestep=0.1; 
    const double tempo=1000000.0;
    const int threads=1; // paralell threads
    int ID[3]; // ID of each thread

    /* Defining physical parameters of the chain*/
    /* Variables */
    const double temp_diff=args2_9[0]; // (argv[2]) (temp_diff=0.5)
    const double chains_ratio=args2_9[1]; // (argv[3]) (chains_ratio=1.0)
    const double mid_spring_const=args2_9[2]; // (argv[4]) (mid_spring_const=1.0)
    const double temp_mean=args2_9[3]; // (argv[5]) (temp_mean=0.07)
    const double poly_power=args2_9[4]; // (argv[6]) (poly_power=2.0)

    /* Fixed or Calculated */
    const double left_ext_potential=args2_9[5]/(4.0*M_PI*M_PI);_diff  // (argv[7]?) (VL=2.0/(4.0*M_PI*M_PI))
    const double right_ext_potential = chains_ratio*left_ext_potential; 
    /*Uniting everything in an array*/
    double Ext_Potential[2]={left_ext_potential, right_ext_potential}; 

    const double left_spring_const=args2_9[6]; // (argv[8]) (left_spring_const=1.0)
    const double right_spring_const = chains_ratio*left_spring_const; 
    /*Uniting everything in an array*/
    double Spring_Const[3]={left_spring_const, mid_spring_const, right_spring_const};

    const double drag_coefficient=args2_9[7]; // (argv[9]) (drag_coefficient=1.0)

    /*Amplitude of heat baths*/
    const double left_temp  = temp_mean*(1.0+temp_diff); // Reduced temperature in the hot bath (1/(m*beta))
    const double right_temp = temp_mean*(1.0-temp_diff); // Reduced temperature in the cold bath (1/(m*beta))
    const double left_amp  = 2.0*drag_coefficient*left_temp; 
    const double right_amp = 2.0*drag_coefficient*right_temp; 
    double Baths_Amp[2] = {left_amp, right_amp};  // Uniting in an array

    /* Declaring RNG seeds */
    unsigned int semente[2];

    /* Initial Values */
    double x0[number_of_particles];
    double v0[number_of_particles];

    for(n=1;n<number_of_particles-1;n++){
        x0[n] = 0.0;
        v0[n] = 0.25;
    }
    x0[0] = 0.0;
    x0[number_of_particles-1] = 0.0;
    v0[0] = 1.0;
    v0[number_of_particles-1] = 0.0;







/* ====================== Run Simulation =========================== */

/* Parallel programming */
// omp_set_num_threads(threads); // Number of threads
// #pragma omp parallel
// {

// #pragma omp for
    for(l=1;l<=number_of_samples;l++){

        //printf("(%d, ID %d)\n",l%(number_of_samples/threads), ID[omp_get_thread_num()]); // Output which simulation we are at.

        /* Simulating the system */
        //ID[omp_get_thread_num()] = omp_get_thread_num();
        // strcpy(plot_srk4[ID[omp_get_thread_num()]], "data/thermal_diode_");
        // sprintf(number3[ID[omp_get_thread_num()]], "%d", l);
        // strcat(plot_srk4[ID[omp_get_thread_num()]], number3[ID[omp_get_thread_num()]]);
        // strcat(plot_srk4[ID[omp_get_thread_num()]], dtpe[ID[omp_get_thread_num()]]);

        ID[0]=0;
        strcpy(plot_srk4[ID[0]], data_dir); // (argv[10]?) "data/poly_power_2.0_t106/thermal_code_"
        sprintf(number3[ID[0]], "%d", l);
        strcat(plot_srk4[ID[0]], number3[ID[0]]);
        strcat(plot_srk4[ID[0]], dtpe[ID[0]]);
	

        semente[0] = (2*l) + 1;
        semente[1] = (2*l);
        diode_srk4(plot_srk4[ID[0]], x0, v0, number_of_particles, drag_coefficient, 
                   Ext_Potential, poly_power, timestep, tempo, Baths_Amp, Spring_Const, semente);

    }

//}
    return 1;
}
