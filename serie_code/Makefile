COMPILADOR=gfortran
OPTIMITZADOR=-O2

datafile=thermodynamics.dat
compressed_file=P1_serie.tar.gz
EPSS=$(wildcard *.png)
dirs2tar="Figs_SIunits" "Figs_redunits" "VMD" "simulation_data" "statistics" "trajectories"
files2tar = $(wildcard *.f90) $(wildcard *.gnu) $(wildcard *.txt) main
name4tar="EIA_serie.tar"

## main: generates the main program
main : initialize.o forces.o integrate.o thermodynamics.o binning_gestor.o main.o
	$(COMPILADOR) $(OPTIMITZADOR) -o $@ $^

%.o : %.f90
	$(COMPILADOR) -c $^

# binning : binning.o
#	$(COMPILADOR) $(OPTIMITZADOR) -o $@ $^

## run: executes the main program
.PHONY : run 
run: main
	mkdir -p statistics
	./main < input.txt > SERIE.dat
	gnuplot statistics2plots.gnu
	mv fit.log statistics
	mv *_mean.dat statistics
	mv *_block.dat statistics
	mkdir -p statistics/Plots
	mv *.png statistics/Plots

## plot: runs plots from data in datafile and move them
.PHONY: plot
plot: $(datafile)
	gnuplot data2plots.gnu  
	mkdir -p Figs_redunits
	mv *.png Figs_redunits 
	gnuplot data2plots_adSI.gnu 
	mkdir -p Figs_SIunits
	mv *.png Figs_SIunits
	mkdir -p trajectories
	mv *.xyz trajectories
	mkdir -p simulation_data
	mv *.dat simulation_data
	gnuplot VMD/VMDplots.gnu
	mv *.png VMD
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
	rm -f *.tar

## compress: compresses all the .f90 files and the Makefile to a .tar.gz file
.PHONY : compress
compress: 
	tar -cvf $(name4tar) $(files2tar) $(dirs2tar)
	rm -f -rR $(dirs2tar) $(files2tar)

## decompress: decompresses .tar.gz in current directory
.PHONY: decompress
decompress: 
	tar -xvf $(name4tar)
	rm $(name4tar)


## variables: prints the variables used: the compiler and the optimitzator
.PHONY : variables
variables:
	@echo COMPILADOR: $(COMPILADOR)
	@echo OPTIMITZADOR: $(OPTIMITZADOR)
	@echo datafile: $(datafile)
	@echo compressed_file: $(compressed_file)
	@echo EPSS: $(EPSS)
	
## help: displays this help message
.PHONY : help
help:
	@grep '^##' Makefile
