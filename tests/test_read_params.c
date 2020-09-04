#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../include/read_params.h"

int main(int argc, char *argv[]){

    char *arg_values[2] = {"./tests/test_read_params", "./tests/test_input_params.txt"};

    int is_input_valid = valid_input(2,arg_values);

    if (is_input_valid==1){
        printf("CORRECT! \n");
    }
    else {
        printf("FALSE! \n");
    }
    return 0;
}