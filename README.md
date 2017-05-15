# FESpaces

## Compiling

### Poisson1D

`g++ -std=c++11 -march=native -O3 -o Poisson1D Poisson1D.cpp`

### Poisson2D

`g++ -std=c++14 -march=native -O3 -o Poisson2D Poisson2D.cpp`

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
