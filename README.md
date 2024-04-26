<div style="width: 65%;">
<p align="right" width="100%">
   <img src="https://github.com/Eines-Informatiques-Avancades/awesome-readme/blob/master/icon.png" alt="Logo" width="75" style="position: absolute; top: 10px; right: 10px;">
</p>  
   
## MD code 	
This repository hosts Fortran90 programs tailored for serial molecular dynamics simulations.  
Data analysis of the generated data is conducted utilizing gnuplot.


## Table of Contents
1. [Install](#install)
   1. [Requirements](#requirements)
2. [Serie code](#serie-code)
	1. [Usage](#usage)
	2. [Program structure](#program-structure)
	   1. [Main Program](#main-program)
	   2. [Modules](#modules)
3. [Parallel code](#parallel-code)
	1. [Usage](#usage)
4. [Expected Results](#expected-results)
   1. [Justification of the results](#justification-of-the-results)
5. [Contributors](#contributors)


## Installation
In order to run the program, the user must download the prerequisites listed below.

### Requirements
- gfortran: The GNU Fortran compiler is required to compile and execute Fortran 90 code.
  - Install gfortran in Ubuntu with: 

    ```bash
    sudo apt-get install gfortran
    ```
- gnuplot: Gnuplot is a command-line driven plotting tool. It's used for generating plots.
  - Install gnuplot in Ubuntu with:

    ```bash
    sudo apt-get install gnuplot
    ```
## Serie Code

### Usage
To run the program, type:  

```
make main
make run
```
Compiles the modules and generates the main program. Then executes the main program.  

```
make stats
make plot
```
Generates statiscs file with info on mean and standard deviation of observables. Runs plots from data in datafile and moves them to a new folder.  

For more commands type make help

###  Program structure

####  Main program
Main.f90 contains the main simulation loop. It requires the use of the modules detailed in the subsequent section.  
Includes the system initialization, thermalization, and production run processes.  
Data files 'thermodynamics.dat' and 'resultsrdflong_def.dat' will be generated containing the simulation results.

####  Modules
The modules required to compile the main program are:

- forces.f90: Allows for the calculation of forces acting on the particles.
- initialize.f90: Contains subroutines for system initialization and applying periodic boundary conditions.
- integrate.f90: Contains subroutines for the Velocity Verlet integrator, Andersen thermostat, and Gaussian distribution.
- thermodynamics.f90: Contains subroutines for calculating different observables of the simulation.  

Additionally, the following files are used for proper functioning:
- input.txt: Includes input simulation parameters such as the number of particles, timestep, cutoff, temperature, among others.
- data2plots.gnu: Generates figures from the data generated in the simulation.  
- data2plots_adSI.gnu: Generates figures in SI units.  
- binning.f90: Performs block average method and calculates the mean and standard deviation of observables.

## Parallel Code

### Usage
To run the parallel program, type:  


```
make enter_interactive
module load intel_compiler_suite/2021.3
module load openmpi/4.1.4_ics-2021.3
```
Enters interactive environment and loads required modules.

```
make main
make runX
```
Compile the modules and generates the main program.  
Then execute the main program. Replace X with the number of cores you want to use to run the program.

```
make stats
make plot
```
Generates statiscs file with info on mean and standard deviation of observables. Runs plots from data in datafile and moves them to a new folder.  
For more commands type make help

##  Expected Results
The generated plots should look like:  

![Energy components evolution over time](https://github.com/Eines-Informatiques-Avancades/Project-I/blob/master/serie_code/Figs_SIunits/Esvst_SI.png)  
Figure 1: Energy components evolution during the simulation.

![Pressure vs Temperature in reduced units](https://github.com/Eines-Informatiques-Avancades/Project-I/blob/master/serie_code/Figs_SIunits/PvsT_SI.png)  
Figure 2: Pressure of the system respect to its temperature.

![Radial distribution function](https://github.com/Eines-Informatiques-Avancades/Project-I/blob/master/serie_code/Figs_redunits/RDF.png)  
Figure 3: Radial distribution function of the simulation.

###  Justification of the results
In this simulation, the system was modeled at a temperature T = 300K and a density of 3.749 g/L, mirroring the conditions typical of Kr gas.
The interaction between particles was defined using the Lennard-Jones potential, with parameters ε/kB = 162.58 K and σ = 3.6274 angstrom.  
The parameters were extracted from Rutkai et al. (Rutkai G.; Thol M.; Span R.; Vrabec J. How Well Does the Lennard-Jones Potential Represent the Thermodynamic Properties of Noble Gases?. Mol. Phys. 2017, 115, 1104–1121.)  

## Contributors
[Pau Franquesa](https://github.com/PFranqV) (Initial conditions and pbc)  
[Pau Eritja](https://github.com/PauEritja) (Forces)  
[Safae el Jelloui](https://github.com/SafaJeBo) (Statistics and visualization of results)  
[Marcos Carrera](https://github.com/Marcos-C-A) (Integration)  
[Berta Bori](https://github.com/bbobru) (Coordination and integration of the code)

</div>
