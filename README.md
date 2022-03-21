# installation

tested with GNU Fortran 11.2.0, modify `Makefile` to mix it up. this is actually two separate programs, `select_data` and `add_data`. for a description of what they do, see `doc/main.pdf`.

# usage

`select_data [data set file] [number of points to select] [initial temperature]`

`add_data [fixed data set file] [candidate data set file] [number of points to select] [initial temperature]`

the data set files should have as the first line `# N M`, where N indicates how many data points there are, and M indicates how many components each data point has. this first line is then followed by the data itself, which should just be a N by M table of numbers. for example, in python, `numpy.loadtxt()` should be able to parse the data and result in an array with shape `(N, M)`.

# output

the first few lines of output will echo how the program parsed your inputs. check these lines to ensure that you are getting what you really want. once annealing begins, the total energy and acceptance rate at each temperature is printed out.

when annealing finishes, the program will perform a check on the total energy. the reason for this is that at each iteration it is the change in energy, not the total energy, which is calculated. this calculation is independent of the total energy calculation, and obviously much cheaper. whenever an exchange is accepted, this change in energy is added to the previous total energy to estimate the new total energy. this will NOT be identical to the final calculation of the total energy, because of significant rounding. so long as these numbers match in the first 6 or so digits, there is no cause for alarm, because the total energy is never relevant to the annealing procedure.

after this check the indices of the selected points are printed, with the indices starting from 0 (as opposed to 1, which is how it is in Fortran). after the indices are printed, the respective data points are also printed, in table format. the print formatting is left up to the compiler via `print *`. if this proves to be a hinderance, the formatting may be cleaned up in the future.
