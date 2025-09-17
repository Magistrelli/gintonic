# GINTONIC - Gross-pitaevskii INtegrator on a TOrus with Non-zero vortICity

Main developer: [Fabio Magistrelli](https://orcid.org/0009-0005-0976-7851) - fabio.magistrelli@uni-jena.de  
Co-developer: [Marco Antonelli](https://orcid.org/0000-0002-5470-4308) - antonelli@lpccaen.in2p3.fr

## Description

The code integrates the Gross-Pitaevskii equation (GPE) to simulate the dynamics of a superfluid made of spinless particles and immersed in an external background potential.
A set of Quasi-Periodic Boundary Conditions (QPBC) is employed to study the bulk properties of the superfluid, allowing for systems with non-zero net vorticity to be represented on a torus.

Reference paper (ref [1] in the code): [***TODO*** link to the paper](link)



## Scripts description

First, run
`make`
to compile.

You can run `./run_tests.sh` to check the installation.



[`run.sh`](run.sh): Automatically runs `GINTONIC` and optionally some post-process scripts (`plots.py`) after creating the directories and copied the files needed to run and store the results (help: `./run.sh -h`; automatic compilation with flag `-c`)  
[`plots.py`](plots.py): Plots all the results (initial condition, time evolution of the phase, density and velocity of the superfluid and of the QPBC integration constants, conserved quantities, final condition)  
[`run_tests.sh`](run_tests.sh): Run the standard test routine comparing against the results of the last commit of the branch `main`  

[`src`](src): source code; the main scripts and the libraries and their headers have the `.C`,`.cpp` and `.h` extensions, respectively  

- [`src/gin.C`](gin.C): Simulates the dynamic evolution of a superfluid immersed in a certain pinning landscape and initialized with either a certain lattice of vortices, random or uniform density and phase or a couple of opposite soliton waves. It can accept `<vBgT>` and `<vBgEv>` parameters from terminal to evolve with a time-dependent background velocity (see usage in header)  

- [`Vector2.*`](src/Vector2.h): Defines the classes Vector2 and ComplexFunc, which respectively represent a 2D vector (or complex number) and a complex (or 2D) function  
- [`functions.*`](src/functions.h): Defines some useful general functions that are used by the classes defined in Jacobi.h, Pinning.h, Superfluid.h and PDE.h  
- [`Jacobi.*`](src/Jacobi.h): Implements the Jacobi theta function and the initial superfluid's phase condition described in ref. [1], Ch. 7  
- [`Random.*`](src/Random.h): Random numbers generator; needs `Primes` and `seed.in` to initial the random number generator  
- [`Pinning.*`](src/Pinning.h): Implements the class Pinning, which defines the pinning potential as a scalar field  
- [`Superfluid.*`](src/Superfluid.h): Implements the class Superfluid, which describes a superfluid configuration  
- [`PDE.*`](src/PDE.h): Implements the classes PDE to solve a generic partial differential equation and GPE (which inherits from PDE and Superfluid) to solve the GPE and then simulate the superfluid's dynamics

[`parameters`](parameters): Contains ready-to-use example parameter files;  

- [`parameters_example.dat`](parameters/parameters_example.dat): Sets the physical parameters of the superfluid and the computational parameters of the space time grid  
- [`inVort_example.dat`](parameters/inVort_example.dat): Define the initial vortex lattice  
- [`prepare_inVort_set.py`](parameters/prepare_inVort_set.py): Generates a series of inVort.dat parameter files for a specified array of x and y position (single charged vortices)  

[`tests`](tests): Contains a series of folders with parameter files with standard tests (run with `run_tests.sh`);  


### Input files

For the appropriate formatting for global parameter file check [`parameters_example.dat`](parameters/parameters_example.dat).  
Please always use `if_active_particles=0`, as the implementation of the backreaction of the superfluid on the external potenrial is **not** tested.

Appropriate formatting for the vortex input file (e.g. [`inVort_example.dat`](parameters/inVort_example.dat)):
 <br>
 three different columns <br>
 the first line is the legend (also reported below) and the second must be empty <br>
 then, in each row the following parameters must be specified: <br>
```
  xv0[csi_units]  yv0[csi_units]  vortex_charge
```

