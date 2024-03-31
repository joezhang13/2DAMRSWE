# 2DAMRSWE
This is a finite volume solver for 2D nonlinear shallow-water equations. 
A two-level adaptivel mesh refinement (AMR) is used with a physics-based indicator. 
The program runs in parallel under the Message Passing Interface (MPI).

The report of this project can be found at:\
https://drive.google.com/file/d/18O7br3vqqR6iqZmpEcgfYDgKxafIp_a9/view?usp=sharing

## Dependencies
gcc/8.2.0\
openmpi/4.0.6

## Usage
The codes in `amr.cpp` simulate the evolution of the water surface in a square domain with a Gaussian-shaped profile in the initial condition. The number of grids $N$ and the number of processors $P$ are specified as the input parameters. In the timing folder, the codes are modified to record the running time of each part of the program with different $P$, which measures the performance of the parallel solver.

## Contact
Zhou Zhang\
Email: joezhang@umich.edu
