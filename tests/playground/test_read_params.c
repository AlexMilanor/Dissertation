#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "read_params.h"

int main(int argc, char *argv[]){

    int is_input_valid = valid_input(argc,argv);

    if (is_input_valid==1){
        printf("Valid Input! \n");
    }
    else {
        printf("Invalid Input! \n");
    }
    return 0;
}