#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Pi decimal expansion taken from the Online Encyclopedia of Integer Sequences
// Source: http://oeis.org/A000796
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279
#endif /* M_PI */

// 1/Pi decimal expansion taken from the Online Encyclopedia of Integer Sequences
// Source: http://oeis.org/A049541
#ifndef M_1_PI
#define M_1_PI 0.318309886183790671537767526745
#endif /* M_1_PI */

void main(){

    double a = 10.0*4.0*M_PI*M_PI;
    double b = 15.0*4.0*M_PI*M_PI;
    double c = 20.0*4.0*M_PI*M_PI;
    double d = 25.0*4.0*M_PI*M_PI;
    double e = 30.0*4.0*M_PI*M_PI;
    double f = 35.0*4.0*M_PI*M_PI;
    double g = 40.0*4.0*M_PI*M_PI;


    printf("Valor de 10.0: %.15lf\n",a);
    printf("Valor de 15.0: %.15lf\n",b);
    printf("Valor de 20.0: %.15lf\n",c);
    printf("Valor de 25.0: %.15lf\n",d);
    printf("Valor de 30.0: %.15lf\n",e);
    printf("Valor de 40.0: %.15lf\n",f);
    printf("Valor de 45.0: %.15lf\n",g);

}
