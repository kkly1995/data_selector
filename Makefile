selector:
	gfortran -Ofast -o bin/select_data src/misc.f90 src/tables.f90 src/selector.f90
	gfortran -Ofast -o bin/select_short src/misc.f90 src/tables.f90 src/short_selector.f90

adder:
	gfortran -Ofast -o bin/add_data src/misc.f90 src/tables.f90 src/adder.f90
	gfortran -Ofast -o bin/add_short src/misc.f90 src/tables.f90 src/short_adder.f90

clean:
	rm -f *mod
