#include <stdlib.h>
#include <stdio.h>
#include <string.h>

const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

int main(int argc, char **argv){

    printf("%d\n",argc);
    printf("%s\n",argv[0]);

    char *dot = strrchr(argv[1], '.');

    if (strcmp(dot+1, "txt")==0){
        printf("%s\n",dot+1);        
    }
    else {printf("Format not supported.\n");}

    return 0;
}