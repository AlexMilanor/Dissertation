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
#include "initial_conditions.h"



int main(int argc, char **argv){

    /* ========================= Simulation Parameters ========================*/

    /* system parameters */
    const int number_of_samples = 1; 
    const int number_of_particles = 50;
    

    /* simulation time */
    const double timeStep = 0.01; 
    const double timeRange = 125000.0;
    const double timeTransient = 100000.0;


    /* Defining physical parameters of the chain*/
    /* Independent Variables */
    const double temp_mean = 0.07; // Average Temperature between heat baths
    const double temp_diff = 1.0; // Temperature difference between heat baths (in % of temp_mean)
    const double chains_ratio = 1.0; // ratio = k_Right / k_Left = V_Right / V_Left
    const double mid_spring_const = 1.0; // Spring constant between left and right sides of the chain
    const double poly_power = 2.0; // Exponent of the potential on center (chains interphase)


    /* External Potential */
    const double left_ext_potential_amp = 5.0/(4.0*M_PI*M_PI);
    const double right_ext_potential_amp = chains_ratio*left_ext_potential_amp; 
    double Ext_Potential_Amp[2] = {left_ext_potential_amp, right_ext_potential_amp};


    /* Interaction Potential */
    const double left_spring_const = 1.0;
    const double right_spring_const = chains_ratio*left_spring_const; 
    double Spring_Const[3]={left_spring_const, mid_spring_const, right_spring_const};


    /* Heat Baths */
    const double drag_coefficient=1.0;
    const double left_temp  = temp_mean*(1.0+0.5*temp_diff); // Reduced temperature in the left bath (1/(m*beta))
    const double right_temp = temp_mean*(1.0-0.5*temp_diff); // Reduced temperature in the right bath (1/(m*beta))
    const double left_amp  = 2.0*drag_coefficient*left_temp; // Noise amplitude in the left bath
    const double right_amp = 2.0*drag_coefficient*right_temp; // Noise amplitude in the right bath
    double Baths_Amp[2] = {left_amp, right_amp};

    /* Initial Values */
    double x_0[number_of_particles];
    double v_0[number_of_particles];

    initialize_positions(x_0, number_of_particles);
    initialize_velocities(v_0, number_of_particles, temp_mean);


/* ====================== Simulations =========================== */

    for (int n=1 ; n <= number_of_samples; n++) {

        /* ==================== Simulation control variables ========================= */
        /* Defining output filename */
        /* Output filename : [data_dir]/thermal_diode_[n].[extension] */
        char archname[80];

        char data_dir[] = "data/diode_20200716/"; // Name of directory
        char plot_filename[] = "thermal_diode_"; // Name of file
        char sim_number[10]; // Number of simulation
        sprintf(sim_number, "%d", n);
        char output_extension[] = ".csv"; // File extension

        strcpy(archname, data_dir);
        strcat(archname, plot_filename);
        strcat(archname, sim_number);
        strcat(archname, output_extension);

        /* ====================== Running Simulations =========================== */
        /* Declaring RNG seeds */
        unsigned int semente[2];
        semente[0] = (2*n) + 1;
        semente[1] = (2*n);
        
        diode_srk4(archname, 
                   x_0, v_0, number_of_particles, 
                   drag_coefficient, 
                   Ext_Potential_Amp, 
                   poly_power, 
                   timeStep, 
                   timeRange,
                   timeTransient,
                   Baths_Amp, 
                   Spring_Const, 
                   semente);

        printf("Fim da Simulação");

    }

    return 0;
}
