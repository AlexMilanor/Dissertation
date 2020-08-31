/*
*********************************************
*                                           *
*       Simulation of Fourier Chain         *
*             for sample                    *
*                                           *
*********************************************

Author: Alexandre A. A. Almeida
First Created: 01/04/2019

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
#include "read_params.h"

int main(int argc, char *argv[]){

    /* 
    ========================= Simulation Parameters ========================
    */

    /*
    ======================== Parsing the input file ========================
    */

    if (valid_input(argc, argv)==0){
        printf("Invalid input file!\n");
        return 0;
    }

    char input_filename[100];
    strcpy(input_filename, argv[1]);

    int values_system[2]
    double values_simulations[3]
    double values_physics[8]
    char values_files[3][100];

    get_input_params(input_filename, values_system, values_simulations, values_physics, values_files)

    /* Position of each variable in the arrays:

    PARAMS_SYSTEM = {0:"Number_of_Samples", 1:"Number_of_Particles_in_the_chain"};

    PARAMS_SIMULATION = {0:"Time_step", 1:"Total_Simulation_Time", 2:"Transient_Time"};

    PARAMS_PHYSICS = {0:"Mean_Temperature", 1:"Temperature_Difference_(%)", 2:"Left_to_Right_Ratio",
                      3:"Interphase_Spring_Constant", 4:"Interphase_Polynomial_Power",
                      5:"Left_Potential_Amplitude", 6:"Left_Spring_Constant", 7:"Drag_Coefficient"};

    PARAMS_FILES = {0:"Data_Directory", 1:"Output_filename", 2:"Output_extension"};
    */


    /*
    =================== Read the params from the input file ========================
    */

    /* system parameters */
    const int number_of_samples = values_system[0]; 
    const int number_of_particles = values_system[1];
    

    /* simulation time */
    const double timeStep = values_simulations[0]; 
    const double timeRange = values_simulations[1];
    const double timeTransient = values_simulations[2];


    /* Defining physical parameters of the chain*/
    /* Independent Variables */
    const double temp_mean = values_physics[0]; // Average Temperature between heat baths
    const double temp_diff = values_physics[1]; // Temperature difference between heat baths (in % of temp_mean)
    const double chains_ratio = values_physics[2]; // ratio = k_Right / k_Left = V_Right / V_Left
    const double mid_spring_const = values_physics[3]; // Spring constant between left and right sides of the chain
    const double poly_power = values_physics[4]; // Exponent of the potential on center (chains interphase)


    /* External Potential */
    const double left_ext_potential_amp = values_physics[5]*0.25*M_1_PI*M_1_PI; // Normalize dividing by (2\pi)^2
    const double right_ext_potential_amp = chains_ratio*left_ext_potential_amp; 
    double Ext_Potential_Amp[2] = {left_ext_potential_amp, right_ext_potential_amp};


    /* Interaction Potential */
    const double left_spring_const = values_physics[6];
    const double right_spring_const = chains_ratio*left_spring_const; 
    double Spring_Const[3]={left_spring_const, mid_spring_const, right_spring_const};


    /* Heat Baths */
    const double drag_coefficient= values_physics[7];
    const double left_temp  = temp_mean*(1.0+0.5*temp_diff); // Reduced temperature in the left bath (1/(m*beta))
    const double right_temp = temp_mean*(1.0-0.5*temp_diff); // Reduced temperature in the right bath (1/(m*beta))
    const double left_amp  = 2.0*drag_coefficient*left_temp; // Noise amplitude in the left bath
    const double right_amp = 2.0*drag_coefficient*right_temp; // Noise amplitude in the right bath
    double Baths_Amp[2] = {left_amp, right_amp};

    /*
    =================== Initialize velocity and position ==================
    */

    /* Initial Values */
    double x_0[number_of_particles];
    double v_0[number_of_particles];

    initialize_positions(x_0, number_of_particles);
    initialize_velocities(v_0, number_of_particles, temp_mean);


    /* 
    ===================== Simulations ===========================
    */

    for (int n=1 ; n <= number_of_samples; n++) {

        /* ==================== Simulation control variables ========================= */
        /* Defining output filename */        
        char data_dir[100]; // Name of directory
        strcpy(data_dir, values_files[0]);

        char plot_filename[100]; // Name of file
        strcpy(plot_filename, values_files[1]);

        char sim_number[10]; // Number of simulation
        sprintf(sim_number, "%d", n);

        char output_extension[100]; // File extension
        strcpy(output_extension,values_files[2]);


        char archname[100]; // Final name

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

        printf("Fim da Simulação\n");

    }

    return 0;
}
