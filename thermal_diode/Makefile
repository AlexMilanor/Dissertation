thermal_diode: thermal_diode.o
	gcc  -o ./thermal_diode ./thermal_diode.o -lm -L/usr/lib/x86_64-linux-gnu/ -lgsl -lgslcblas

thermal_diode.o: thermal_diode.c
	gcc -I/usr/include/gsl -c ./thermal_diode.c -o ./thermal_diode.o

clean:
	rm thermal_diode.o thermal_diode