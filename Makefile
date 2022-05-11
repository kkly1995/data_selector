selector:
	gfortran -fopenmp -O3 -o bin/select_data src/misc.f90 src/tables.f90 src/selector.f90

adder:
	gfortran -fopenmp -O3 -o bin/add_data src/misc.f90 src/tables.f90 src/adder.f90

clean:
	rm -f *mod
