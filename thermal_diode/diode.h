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
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <math.h>

double potential(double x, double V){
    /* Force due to the external potential */
    double a0=1; // external potential length
    return -(V*2.0*M_PI/a0)*sin((2.0*M_PI/a0)*x);
}

double force_bath(double gamma, double k, double x0, double x1, double v0, double V){
/* Defining the deterministic forces on the Langevin Bath */
    return -gamma*v0 - k*(x0-x1)+potential(x0, V);
}

double spring_force(double k, double x1, double x2){
/* Force of interaction between neighboring particles */
    return -k*(x1 - x2);
}

double force(double k1, double k2, double x0, double x1, double x2, double V){
/* Defining the deterministic forces on x1 */
    return spring_force(k1, x1, x0) + spring_force(k2, x1, x2) + potential(x1, V);
}

double inter_force(double k, double mu, double x1, double x2){
/* Force in the interphase of the left and right chain */
    double force_interf;
    force_interf = -k*pow(fabs(x1 - x2),mu-2.0)*(x1-x2);
    return force_interf;
}

double dxdt(double v){
/* Deterministic velocity for any particle */
    return v;
}

double x_m(double delta, double x0, double K0, int k){
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

double step_next(double delta, double x0, double K0, double K1, double K2, double K3){
/* Calculating x(t+delta(t)) for SRK4 */
    return x0 + (K0+2.0*K1+2.0*K2+K3)*delta/6.0;
}

int diode_srk4(char *name, double *x0, double *v0, int N_C, double gama, 
               double *V, double mu, double tau, double tempo, double TransientTime, 
               double *A, double *k, unsigned int *semente){

/* Main code for the SRK4 in the harmonic chain */

    int i, m;  // iteration variables
    long int N; // Number of points
    N = (long int)tempo/tau; // N = time_range/tau

    /* Defining the Random Number Generators (RNGs) */
    gsl_rng *re;
    gsl_rng *rd;

    // Using Mersenne Twister generator
    re = gsl_rng_alloc(gsl_rng_mt19937);
    rd = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(re, semente[0]);
    gsl_rng_set(rd, semente[1]);

    /* Defining output files */
    FILE *file; // Declaring file
    file = fopen(name, "w"); // Creating csv files

    /* Creating column headers */
    fprintf(file, "time");
    for(i=0;i<N_C;i++){
        fprintf(file, ",x%d,v%d",i+1,i+1);
    }
    fprintf(file, "\n");

    long int t; // Time
    double epse, epsd; //Noise of heat beaths
    double x[N_C][5]; //Position of particles
    double v[N_C][5]; //Velocity of particles
    double Kx[N_C][4]; //Matrix K (of Stochastic Runge Kutta)
    double Kv[N_C][4]; //Matrix K (of Stochastic Runge Kutta)

    /* Initial Conditions*/
    t = 0;
    fprintf(file, "%.8g",(double)t*tau);
    for(i=0;i<N_C;i++){
        x[i][0]=x0[i];
        v[i][0]=v0[i];
        fprintf(file, ",%.8g,%.8g",x[i][0], v[i][0]);
    }
    fprintf(file, "\n");


    /* Running the program */
    // Iteration for the time interval
    for(t=1;t<N+1;t++){

        // We save the results after the transient time
        if((double)t*tau>TransientTime){
            fprintf(file, "%.8g",(double)t*tau);
        }

        //Noise of the Heat Baths
        epse = gsl_ran_ugaussian(re);
        epsd = gsl_ran_ugaussian(rd);
        /*----------------------- First N/2 particles ---------------------------*/

        /* ---------------------- Three first steps of SRK4 ---------------------*/
        for(m=1;m<4;m++){

            //Left Heat Bath
            Kv[0][m-1] = force_bath(gama, k[0], x[0][m-1], x[1][m-1], v[0][m-1], V[0]);
            Kx[0][m-1] = dxdt(v[0][m-1]);
            v[0][m] = x_m(tau, v[0][0], Kv[0][m-1], m) + 0.5*sqrt(A[0]*tau)*epse;
            x[0][m] = x_m(tau, x[0][0], Kx[0][m-1], m);

            //Right Heat Bath
            Kv[N_C-1][m-1] = force_bath(gama, k[2], x[N_C-1][m-1], x[N_C-2][m-1], v[N_C-1][m-1], V[1]);
            Kx[N_C-1][m-1] = dxdt(v[N_C-1][m-1]);
            v[N_C-1][m] = x_m(tau, v[N_C-1][0], Kv[N_C-1][m-1], m) + 0.5*sqrt(A[1]*tau)*epsd;
            x[N_C-1][m] = x_m(tau, x[N_C-1][0], Kx[N_C-1][m-1], m);

            //Middle Particles
            for(i=1;i<N_C-1;i++){
                if(i < ((N_C/2)-1)){
                    Kv[i][m-1] = force(k[0],k[0], x[i-1][m-1],x[i][m-1],x[i+1][m-1], V[0]);
                    Kx[i][m-1] = dxdt(v[i][m-1]);
                    v[i][m] = x_m(tau, v[i][0], Kv[i][m-1], m);
                    x[i][m] = x_m(tau, x[i][0], Kx[i][m-1], m);
                }
                else if(i > N_C/2){
                    Kv[i][m-1] = force(k[2],k[2], x[i-1][m-1],x[i][m-1],x[i+1][m-1], V[1]);
                    Kx[i][m-1] = dxdt(v[i][m-1]);
                    v[i][m] = x_m(tau, v[i][0], Kv[i][m-1], m);
                    x[i][m] = x_m(tau, x[i][0], Kx[i][m-1], m);
                }
                else if(i == N_C/2 - 1){
                    Kv[i][m-1] = spring_force(k[0], x[i][m-1],x[i-1][m-1]) + inter_force(k[1], mu, x[i][m-1], x[i+1][m-1]) + potential(x[i][m-1], V[0]);
                    Kx[i][m-1] = dxdt(v[i][m-1]);
                    v[i][m] = x_m(tau, v[i][0], Kv[i][m-1], m);
                    x[i][m] = x_m(tau, x[i][0], Kx[i][m-1], m);
                }
                else if(i == N_C/2){
                    Kv[i][m-1] = spring_force(k[2], x[i][m-1],x[i+1][m-1]) + inter_force(k[1], mu, x[i][m-1], x[i-1][m-1]) + potential(x[i][m-1], V[1]);
                    Kx[i][m-1] = dxdt(v[i][m-1]);
                    v[i][m] = x_m(tau, v[i][0], Kv[i][m-1], m);
                    x[i][m] = x_m(tau, x[i][0], Kx[i][m-1], m);
                }
            }
        }
        /* -------------------------------------------------------------------- */


        /* ------------ Fourth Step (Getting x(t+delta)) ------------------------*/
        //Left Heat Bath
        Kv[0][3] = force_bath(gama, k[0], x[0][3], x[1][3], v[0][3], V[0]);
        Kx[0][3] = dxdt(v[0][3]);
        v[0][4] = step_next(tau, v[0][0], Kv[0][0], Kv[0][1], Kv[0][2], Kv[0][3]) + sqrt(A[0]*tau)*epse;
        x[0][4] = step_next(tau, x[0][0], Kx[0][0], Kx[0][1], Kx[0][2], Kx[0][3]);
        
        // We save the results after the transient time
        if((double)t*tau>TransientTime){
            fprintf(file, ",%.8g,%.8g",x[0][4],v[0][4]);
        }

        //Middle Particles
        for(i=1;i<N_C-1;i++){
            if(i < N_C/2 - 1){
                Kv[i][3] = force(k[0], k[0], x[i-1][3], x[i][3], x[i+1][3], V[0]);
                Kx[i][3] = dxdt(v[i][3]);
                v[i][4] = step_next(tau, v[i][0], Kv[i][0], Kv[i][1], Kv[i][2], Kv[i][3]);
                x[i][4] = step_next(tau, x[i][0], Kx[i][0], Kx[i][1], Kx[i][2], Kx[i][3]);
            }
            else if(i == N_C/2 - 1){
                Kv[i][3] = spring_force(k[0], x[i][3],x[i-1][3]) + inter_force(k[1], mu, x[i][3], x[i+1][3]) + potential(x[i][3], V[0]);
                Kx[i][3] = dxdt(v[i][3]);
                v[i][4] = step_next(tau, v[i][0], Kv[i][0], Kv[i][1], Kv[i][2], Kv[i][3]);
                x[i][4] = step_next(tau, x[i][0], Kx[i][0], Kx[i][1], Kx[i][2], Kx[i][3]);
            }
            else if(i == N_C/2){
                Kv[i][3] = spring_force(k[2], x[i][3],x[i+1][3]) + inter_force(k[1], mu, x[i][3], x[i-1][3]) + potential(x[i][3], V[1]);
                Kx[i][3] = dxdt(v[i][3]);
                v[i][4] = step_next(tau, v[i][0], Kv[i][0], Kv[i][1], Kv[i][2], Kv[i][3]);
                x[i][4] = step_next(tau, x[i][0], Kx[i][0], Kx[i][1], Kx[i][2], Kx[i][3]);
            }
            else if(i > N_C/2){
                Kv[i][3] = force(k[2], k[2], x[i-1][3], x[i][3], x[i+1][3], V[1]);
                Kx[i][3] = dxdt(v[i][3]);
                v[i][4] = step_next(tau, v[i][0], Kv[i][0], Kv[i][1], Kv[i][2], Kv[i][3]);
                x[i][4] = step_next(tau, x[i][0], Kx[i][0], Kx[i][1], Kx[i][2], Kx[i][3]);
            }

        // We save the results after the transient time
        if((double)t*tau>TransientTime){
            fprintf(file, ",%.8g,%.8g",x[i][4],v[i][4]);
        }

    }
        //Right Heat Bath
        Kv[N_C-1][3] = force_bath(gama, k[2], x[N_C-1][3], x[N_C-2][3], v[N_C-1][3], V[1]);
        Kx[N_C-1][3] = dxdt(v[N_C-1][3]);
        v[N_C-1][4] = step_next(tau, v[N_C-1][0], Kv[N_C-1][0], Kv[N_C-1][1], Kv[N_C-1][2], Kv[N_C-1][3]) + sqrt(A[1]*tau)*epsd;
        x[N_C-1][4] = step_next(tau, x[N_C-1][0], Kx[N_C-1][0], Kx[N_C-1][1], Kx[N_C-1][2], Kx[N_C-1][3]);


        // We save the results after the transient time
        if((double)t*tau>TransientTime){
            fprintf(file, ",%.8g,%.8g",x[N_C-1][4],v[N_C-1][4]);
        }
        /*---------------------------------------------------------------------------*/

        //Next iteration
        for(i=0;i<N_C;i++){
            x[i][0]=x[i][4];
            v[i][0]=v[i][4];
        }

        // We save the results after the transient time
        if((double)t*tau>TransientTime){
            fprintf(file, "\n");
        }
    }

    /* Cleaning memory */
    fclose(file); // Closing files
    gsl_rng_free(re); // Closing left RNG
    gsl_rng_free(rd); // Closing right RNG

    return 0;
}
