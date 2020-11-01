# Fortran Air Pollution Lattice Boltzmann
## Compilation
All code must be compiled using f95. A standard compilation would be
```bash
gfortran -std=f95 -c LB_D3Q19.f90
gfortran -std=f95 -c Constants.f90
gfortran -std=f95 main.f90 LB_D3Q19.f90 Constants.f90
```
However, makefiles are included to save time and effort in this matter for some files, 
and diferent debugging, profiling and optimization checks are also automatized. Please do 
check the following table:

### Implemented Automatization
| Shell Command   | What it does oversimplified                                               |
|-----------------|---------------------------------------------------------------------------|
| make            | Compiles using f95 standard and runs program using time                   |
| make debug      | Compiles with all warning flags, saves debug info using -g flag and runs  |
| make valgrind   | Compiles with -g flag and runs using valgrind                             |
| make cachegrind | Same as valgrind but using runs using cachegrind                          |
| make gprof      | Compiles with -pg flag, runs and display gprof info                       |
| make perf       | Same as gprof, but runs and displays perf info                            |