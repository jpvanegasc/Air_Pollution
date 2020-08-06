# Air Pollution Lattice Boltzmann
Computational Fluid Dynamics model of difusing pollutants over a city block.

## Author
* Juan Pablo Vanegas Correa. Contact: jpvanegasc@unal.edu.co.

## Acknowlegements
* Jose Daniel Muñoz
* Pierre Sagaut

## References
* Wind comfort assessment by means of large eddy simulation with lattice Boltzmann method in full scale city area

## Notes on Programming Languages
This project was built primary on C++ and CUDA. C++ was used primarily to make small developments and 
test its functionality. CUDA is intended for simulating bigger systems, and the main results from this work are intended to come from the CUDA module.
Just for fun I've included a Lattice-Boltzmann developed in FORTRAN, but It's not my priority to 
reproduce results using it.

Programming languages are divided by folders, and each one of them has 
a different README with its respective compilation guide.