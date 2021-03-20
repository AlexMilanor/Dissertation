.PHONY: clean


thermal_diode_2: src/thermal_diode_2.o
	gcc -static -o ./thermal_diode_2 ./src/thermal_diode_2.o -L/usr/lib/x86_64-linux-gnu/ -lgsl -lgslcblas -lm

src/thermal_diode_2.o: src/thermal_diode_2.c
	gcc -I/usr/include/gsl -c ./src/thermal_diode_2.c -o ./src/thermal_diode_2.o

thermal_diode: src/thermal_diode.o
	gcc -static -o ./thermal_diode ./src/thermal_diode.o -L/usr/lib/x86_64-linux-gnu/ -lgsl -lgslcblas -lm

src/thermal_diode.o: src/thermal_diode.c
	gcc -I/usr/include/gsl -c ./src/thermal_diode.c -o ./src/thermal_diode.o

tests/test_read_params: tests/test_read_params.o
	gcc  -o ./tests/test_read_params ./tests/test_read_params.o 

tests/test_read_params.o: tests/test_read_params.c
	gcc -c ./tests/test_read_params.c -o ./tests/test_read_params.o

clean:
	rm src/thermal_diode.o tests/test_read_params.o src/thermal_diode_2.o
