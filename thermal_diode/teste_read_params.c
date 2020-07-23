#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char **argv){
    FILE * fr = fopen(argv[1], "rt");

    if (fr == NULL) {
        printf("file %s not found\n", argv[1]);
    }

    char section[50] = "---";
    char tempbuff[50];
    char field[50];
    double value;
    char value_string[50];
    int value_int;
    while(!feof(fr)){
    
        fgets(tempbuff,50,fr);
        int length_tempbuff = strlen(tempbuff); 

        if (length_tempbuff > 1) {

            if (tempbuff[0] == '-') {
                sscanf(tempbuff, "--- %s ---",section);
                printf("%s\n",section);
            }
            else {                

                if (strcmp(section,"Simulation")==0)
                {
                    sscanf(tempbuff, "%s : %lf", field, &value);
                    printf("Simulation: %s = %lf \n",field,value);
                }
                else if (strcmp(section,"System")==0)
                {
                    sscanf(tempbuff, "%s : %d", field, &value_int);
                    printf("System: %s = %d \n",field,value_int);
                }
                else if (strcmp(section,"Physics")==0)
                {
                    sscanf(tempbuff, "%s : %lf", field, &value);
                    printf("Physics: %s = %lf \n",field,value);
                }
                else if (strcmp(section,"Files")==0)
                {
                    sscanf(tempbuff, "%s : %s", field,value_string);
                    printf("Files: %s = %s \n",field,value_string);
                }
                else {
                    printf("ERROR: Section %s; Field %s\n",section,field);
                    break;
                }
            }

        }   
        
    }

    fclose(fr);

    return 0;
}