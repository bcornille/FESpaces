# FESpaces

This repository contains programs for solving the Poisson equation using Mixed Finite Element Methods.

## Compiling

### Poisson1D

`g++ -std=c++11 -march=native -O3 -o Poisson1D Poisson1D.cpp`

### Poisson2D

`g++ -std=c++14 -march=native -O3 -o Poisson2D Poisson2D.cpp`

## Running

Each program requires a filename for input and output.
For example,

### Poisson1D

`./Poisson1D in.json out.json`

### Poisson2D

`./Poisson2D in2d.json out2d.json`

## Input

Input files are specifed and read using the JSON standard.
Example files `in.json` and `in2d.json` show what keywords must be specified for valid input.
In general, the options should be self-explanatory.

### Poisson1D

#### Formulation
 - Standard
 - Mixed
 - Mimetic

#### Function
 - Two
 - ExpX3pX

### Poisson2D

#### Formulation
 - Standard
 - Mixed
 - Dual-Mixed

#### Function
 - SinXSinY

## Output

Both programs output the error in the solution to the terminal.
If plotting is turned on then the program will call the python plotting utilities
and plot the calculated and analytic solutions.
