
#include <stdlib.h>
#include <stdio.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <math.h>

void initialize_positions (double *x_0, int number_of_particles) {

    for (int n = 0; n < number_of_particles; n++) {
        x_0[n] = 0.0;
    }

}

void initialize_velocities (double *v_0, int number_of_particles, double Temperature) {
    /* Declaring Random Number Generators (RNGs) : Mersenne Twister */
    // gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);

    /* Initializing velocities */
    // for (int n = 0; n < number_of_particles; n++) {
        
    //     double std_deviation = sqrt(Temperature);
    //     double initial_velocity = gsl_ran_gaussian(r, std_deviation);
    //     v_0[n] = initial_velocity;
    // }

    /* Getting the Center of mass velocity */
    // double center_of_mass_velocity = 0.0;
    // for (int n = 0; n < number_of_particles; n++) {

    //     center_of_mass_velocity += v_0[n];

    // }
    // center_of_mass_velocity = center_of_mass_velocity/number_of_particles;

    /* Zeroing the center of mass velocity */
    // for (int n = 0; n < number_of_particles; n++) {
    //     v_0[n] = v_0[n] - center_of_mass_velocity;
    // }

    for(int n=0; n<number_of_particles;n++){
        v_0[n] = 0.0;
    }


}