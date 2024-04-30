
# MAKEFILE FOR PARALEL
# This Makefile is prepared to send jobs to cluster CERQT2. However, it supports interactive node usage.  
## For interactive node calculations load openmpi/4.1.4_ics-2021.3

COMPILADOR=mpif90
#OPTIMITZADOR=-O2
numproc=2
dirs2tar := Figs_SIunits Figs_redunits paralel_data VMD paralel_plot simulation_data statistics trajectories
files2tar := $(wildcard *.f90) $(wildcard *.gnu) $(wildcard *.sub) $(wildcard *.txt) $(wildcard *.sh) main
name4tar := EIA_paralel.tar

#datafile=thermodynamics.dat
compressed_file=P1_serie.tar.gz
EPSS=$(wildcard *.png)
datafile=$(wildcard thermo_*.dat)
## main: generates the main program
main : MOD_INIT.o forces.o integrate.o thermodynamics.o binning_gestor.o main.o
	$(COMPILADOR) -o $@ $^

%.o : %.f90
	$(COMPILADOR) -c $^

## cluster_ALLRUN: runs code with all number of processors in cluster
.PHONY : cluster_ALLRUN
cluster_ALLRUN:
	rm -f 1procs.out   # Avoid overlaping
	qsub -pe smp 1 -o "1procs.out" submit07_runner.sub
	rm -f 2procs.out   # Avoid overlaping
	qsub -pe smp 2 -o "2procs.out" submit07_runner.sub
	rm -f 4procs.out   # Avoid overlaping
	qsub -pe smp 4 -o "4procs.out" submit07_runner.sub
	rm -f 8procs.out   # Avoid overlaping
	qsub -pe smp 8 -o "8procs.out" submit07_runner.sub
	rm -f 16procs.out   # Avoid overlaping
	qsub -pe smp 16 -o "16procs.out" submit07_runner.sub
	rm -f 32procs.out   # Avoid overlaping
	qsub -pe smp 32 -o "32procs.out" submit07_runner.sub
	rm -f 40procs.out   # Avoid overlaping
	qsub -pe smp 40 -o "40procs.out" submit07_runner.sub

## cluster_compile: compile modules in cluster
.PHONY: cluster_compile
cluster_compile:
	qsub submit07_compiler.sub

.PHONY : run1
run1: main
	mpirun -np 1 main<input.txt

.PHONY : run2
run2: main
	mpirun -np 2 main<input.txt

.PHONY : run4
run4: main
	mpirun -np 4 main<input.txt

.PHONY : run8
run8: main
	mpirun -np 8 main<input.txt

.PHONY : run16
run16: main
	mpirun -np 16 main<input.txt

.PHONY : run32
run32: main
	mpirun -np 32 main<input.txt

.PHONY : run40
run40: main
	mpirun -np 40 main<input.txt

## enter_interactive: enters interactive node in iqtc07
.PHONY : enter_interactive
enter_interactive: 
	qrsh -q iqtc07.q -pe smp 1

## plot: runs plots from data in datafile and move them
.PHONY: plot
plot: $(datafile) 
	mkdir -p statistics
	mkdir -p statistics/Plots
	gnuplot statistics2plots.gnu
	mv *_mean.dat statistics
	mv *_block.dat statistics
	mv *.png statistics/Plots
	mv fit.log statistics
	gnuplot data2plots.gnu  
	mkdir -p Figs_redunits
	mv *.png Figs_redunits 
	gnuplot data2plots_adSI.gnu 
	mkdir -p Figs_SIunits
	mv *.png Figs_SIunits
	gnuplot VMD/VMDplots.gnu
	mv *.png VMD

## paralel_plot: analyzes data in paralel, clears directories
.PHONY: paralel_plot
paralel_plot:
	chmod u+x ALLRUNS.sh
	./ALLRUNS.sh
	gnuplot analysis_parallel.gnu
	mkdir -p paralel_plot
	mv *.png paralel_plot
	mkdir -p paralel_data
	rm md_compile.out
	mv *.out paralel_data
	mv PARALELIZATION.dat paralel_data
	mkdir -p trajectories
	mv *.xyz trajectories
	mkdir -p simulation_data
	mv results_rdf.dat simulation_data
	mv thermo_kin+pot.dat simulation_data
	mv thermo_temp+press.dat simulation_data
	mv thermo_tot+msd.dat simulation_data
	make super_clean

## clean: removes all .o and .mod files
.PHONY : clean
clean:
	rm -f *.o
	rm -f *.mod

## super_clean: removes all .o, .mod, .dat, .tar.gz and .eps files
.PHONY : super_clean
super_clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(EPSS)
	rm -f *.dat
	rm -f *.tar.gz
	rm -f *.err
## compress: compresses all the .f90 files and the Makefile to a .tar.gz file
.PHONY : compress
compress:
	tar -cvf $(name4tar)  $(dirs2tar) $(files2tar)
	rm -f -rR $(dirs2tar) $(files2tar)

## decompress: decompresses .tar.gz in current directory
.PHONY: decompress
decompress: $(name4tar)
	tar -xvf $(name4tar)
	rm $(name4tar)

## variables: prints the variables used: the compiler and the optimitzator
.PHONY : variables
variables:
	@echo COMPILADOR: $(COMPILADOR)
	#@echo OPTIMITZADOR: $(OPTIMITZADOR)
	@echo numproc: $(nuproc)
	@echo datafile: $(datafile)
	@echo compressed_file: $(compressed_file)
	@echo EPSS: $(EPSS)
	
## help: displays this help message
.PHONY : help
help:
	@grep '^##' Makefile







