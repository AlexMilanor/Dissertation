#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

int grab_index(const char **values_array, int array_size, char *string_to_find){

    for (int n=0 ; n <= array_size; n++) {

        if (strcmp(string_to_find,values_array[n])==0){

            return n;

        }

    }

    return -1;

}

int main(int argc, char **argv){

    // Arrays that gives the position of the variables

    const char *PARAMS_SYSTEM[2] = {"Number_of_Samples", "Number_of_Particles_in_the_chain"};
    int values_system[2];


    const char *PARAMS_SIMULATION[3] = {"Time_step", "Total_Simulation_Time", "Transient_Time"};
    double values_simulations[3];


    const char *PARAMS_PHYSICS[8] = {"Mean_Temperature", "Temperature_Difference_(%)", "Left_to_Right_Ratio",
                                    "Interphase_Spring_Constant", "Interphase_Polynomial_Power",
                                    "Left_Potential_Amplitude", "Left_Spring_Constant", "Drag_Coefficient"};
    double values_physics[8];


    const char *PARAMS_FILES[3] = {"Data_Directory", "Output_filename", "Output_extension"};
    char values_files[3][50];




    FILE * fr = fopen(argv[1], "rt");

    // Checking if the file exists
    if (fr == NULL) {
        printf("file %s not found\n", argv[1]);
    }

    char section[50] = "---"; // name of the input_params section
    char tempBuffer[50]; // temporary buffer with the line values
    char parameter[50]; // name of the parameter to be defined
    double value; // value of the parameter defined in input_params
    char value_string[50]; // value of parameters for the output files
    int value_int; // value of parameters for chain and ensemble size

    // While not end of file
    while(!feof(fr)){
    
        // reads a line from the stream, save the number of characters in the line
        fgets(tempBuffer,50,fr);
        int length_tempBuffer = strlen(tempBuffer); 

        // if the line is not empty, we read it
        if (length_tempBuffer > 1) {

            // check for change of input_params section
            if (tempBuffer[0] == '-') {
                sscanf(tempBuffer, "--- %s ---",section);
                printf("%s\n",section);
            }
            else {        

                // for the System section, get what is the physical system (chain and ensemble size)
                if (strcmp(section,"System")==0)
                {
                    sscanf(tempBuffer, "%s : %d", parameter, &value_int);

                    int array_size = sizeof(PARAMS_SYSTEM)/sizeof(PARAMS_SYSTEM[0]);
                    int index = grab_index(PARAMS_SYSTEM, array_size, parameter);
                    values_system[index] = value_int;

                    printf("System: %s = %d \n",parameter,values_system[index]);
                }

                // for the Simulation section, get the numerical method parameters
                else if (strcmp(section,"Simulation")==0)
                {
                    sscanf(tempBuffer, "%s : %lf", parameter, &value);

                    int array_size = sizeof(PARAMS_SIMULATION)/sizeof(PARAMS_SIMULATION[0]);
                    int index = grab_index(PARAMS_SIMULATION, array_size, parameter);
                    values_simulations[index] = value;

                    printf("Simulation: %s = %lf \n",parameter,values_simulations[index]);
                }

                // for the Physics section, get the physical parameters of the system
                else if (strcmp(section,"Physics")==0)
                {
                    sscanf(tempBuffer, "%s : %lf", parameter, &value);

                    int array_size = sizeof(PARAMS_PHYSICS)/sizeof(PARAMS_PHYSICS[0]);
                    int index = grab_index(PARAMS_PHYSICS, array_size, parameter);

                    if (strcmp(parameter, "Left_Potential_Amplitude")){

                        // we normalize the amplitude by dividing it by (4*M_PI*M_PI)
                        value = value*M_1_PI*M_1_PI*0.25
                        // for numerical reasons, we instead multiplied by (1/4)*(1/M_PI)*(1/M_PI)

                    }

                    values_physics[index] = value;

                    printf("Physics: %s = %lf \n",parameter,values_physics[index]);
                }

                // for the Files section, get the output file parameters
                else if (strcmp(section,"Files")==0)
                {
                    sscanf(tempBuffer, "%s : %s", parameter,value_string);

                    int array_size = sizeof(PARAMS_FILES)/sizeof(PARAMS_FILES[0]);
                    int index = grab_index(PARAMS_FILES, array_size, parameter);
                    strcpy(values_files[index], value_string); 

                    printf("Files: %s = %s \n",parameter,values_files[index]);
                }

                // Nonexistant section
                else {
                    printf("ERROR: Section %s; Field %s\n",section,parameter);
                    break;
                }
            }

        }   
        
    }

    fclose(fr);

    return 0;
}