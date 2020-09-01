#include <stdlib.h>
#include <stdio.h>
#include <string.h>


char *get_filename_extension(char *filename){
    /*
    Function that gets the extension written in the filename. Taken from stackoverflow.
    Source:
    https://stackoverflow.com/questions/5309471/getting-file-extension-in-c
    */
    char *dot = strrchr(filename, '.');
    if (!dot || dot == filename) return "";
    return dot+1;

}

int file_exists(char *filename){
    /* 
    Function that checks if there is a file with a given filename and if we can read it.
    Source:
    https://stackoverflow.com/questions/230062/whats-the-best-way-to-check-if-a-file-exists-in-c
    */
    FILE *file;
    if ((file = fopen(filename, "r"))){
        fclose(file);
        return 1;
    }
    return 0;
}


int valid_input(int arg_count, char *arg_values[]){
    if (arg_count > 2){
        printf("Only one input file accepted.\n");
        return 0;
    }
    else if (arg_count != 2){
        printf("No input file defined.\n");
        return 0;
    }
    else {
        char input_filename[100];
        strcpy(input_filename, arg_values[1]);
        
        char *extension = get_filename_extension(input_filename);

        if (strcmp(extension,"txt")!=0){
            printf("Extension not supported\n");
            return 0;
        }
        else {
            if (file_exists(input_filename)){
                return 1;
            }
            else {
                return 0;
            }
        }
    }
}


int grab_index(const char **values_array, int array_size, char *string_to_find){

    for (int n=0 ; n <= array_size; n++) {

        if (strcmp(string_to_find,values_array[n])==0){

            return n;

        }

    }

    return -1;

}

int get_input_params(char *input_filename, 
                      int *values_system, 
                      double *values_simulations, 
                      double *values_physics, 
                      char (*values_files)[100]){

    // Arrays that gives the position of the variables

    const char *PARAMS_SYSTEM[2] = {"Number_of_Samples", "Number_of_Particles_in_the_chain"};


    const char *PARAMS_SIMULATION[3] = {"Time_step", "Total_Simulation_Time", "Transient_Time"};


    const char *PARAMS_PHYSICS[8] = {"Mean_Temperature", "Temperature_Difference_(%)", "Left_to_Right_Ratio",
                                    "Interphase_Spring_Constant", "Interphase_Polynomial_Power",
                                    "Left_Potential_Amplitude", "Left_Spring_Constant", "Drag_Coefficient"};


    const char *PARAMS_FILES[3] = {"Data_Directory", "Output_filename", "Output_extension"};



    FILE * fr = fopen(input_filename, "rt");

    // Checking if the file exists
    if (fr == NULL) {
        printf("File %s not found.\n", input_filename);
    }

    char section[100] = "---"; // name of the input_params section
    char tempBuffer[100]; // temporary buffer with the line values
    char parameter[100]; // name of the parameter to be defined
    double value; // value of the parameter defined in input_params
    char value_string[100]; // value of parameters for the output files
    int value_int; // value of parameters for chain and ensemble size

    // While not end of file
    while(!feof(fr)){
    
        // reads a line from the stream, save the number of characters in the line
        fgets(tempBuffer,100,fr);
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