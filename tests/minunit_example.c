/* file minunit_example.c */
 
 #include <stdio.h>
 #include "minunit.h"
 
 int tests_run = 0;
 
 int foo = 7;
 int bar = 4;
 
 static char * test_foo() {
     mu_assert("error, foo != 7", foo == 7);
     return 0;
 }
 
 static char * test_bar() {
     mu_assert("error, bar != 5", bar == 3);
     return 0;
 }

  static char * all_tests() {
     mu_run_test(test_bar);
     mu_run_test(test_foo);
     return 0;
 }

  
 int main(int argc, char **argv) {
     mu_run_test(test_bar);
     mu_run_test(test_foo);
     printf("Tests run: %d\n", tests_run);
 
    int result = 0;
     return result != 0;
 }