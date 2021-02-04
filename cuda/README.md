# CUDA Lattice Boltzmann
## Compilation
Code must be compiled using the nvcc compiler. A standard compilation would be 
```bash
nvcc main.cu
```
## Warnings
The f array is converted to a 1D array, so:
* Not sure if the use of dim3 objects is useful
* The way I'm accessing positions using threadId may not work.