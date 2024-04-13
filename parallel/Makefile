
# MAKEFILE FOR PARALEL
# !!! COMPILING AND RUNNING IS DONE EXPECTING YOU ARE IN INTERACTIVE MODE !!! 
# TO ACCESS IT DO:
#			qrsh -q iqtc07.q -pe smp 1
#			module load openmpi/2.0.1_ics-2015.0

COMPILADOR=mpif90
#OPTIMITZADOR=-O2
numproc=2

datafile=thermodynamics.dat
compressed_file=P1_serie.tar.gz
EPSS=$(wildcard *.png)

## main: generates the main program
main : MOD_INIT.o main.o
	$(COMPILADOR) -o $@ $^

%.o : %.f90
	$(COMPILADOR) -c $^

binning : binning.o
	$(COMPILADOR) -o $@ $^

## run: executes the main program
.PHONY : run
run: main
	mpirun -np $(numproc) main<input.txt

## enter_interactive: enters interactive node in iqtc07
.PHONY : enter_interactive
enter_interactive: 
	qrsh -q iqtc07.q -pe smp 1

## load_mpi: load mpi module when in interactive mode
.PHONY :load_mpi
load_mpi:
	module load openmpi/2.0.1_ics-2015.0

## plot: runs plots from data in datafile and move them (DOES NOT WORK YET!!!)
.PHONY: plot
plot: $(datafile)
	gnuplot data2plots.gnu  
	mkdir -p Figs_redunits
	mv *.png Figs_redunits 
	gnuplot data2plots_adSI.gnu 
	mkdir -p Figs_SIunits
	mv *.png Figs_SIunits

## stats: generates statistics (DOES NOT WORK YET!!!)
.PHONY: stats
stats: $(datafile) binning
	./binning 
	mkdir -p statistics
	mv *mean.dat statistics

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

## compress: compresses all the .f90 files and the Makefile to a .tar.gz file
.PHONY : compress
compress:
	tar -czvf $(compressed_file) * ../README.md 

## decompress: decompresses .tar.gz in current directory 
.PHONY: decompress
decompress: $(compressed_file)
	tar -xvzf $(compressed_file)

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


