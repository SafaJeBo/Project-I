COMPILADOR=gfortran
OPTIMITZADOR=-O2

datafile=thermodynamics.dat
compressed_file=P1_serie.tar.gz
EPSS=$(wildcard *.png)

## main: generates the main program
main : initialize.o forces.o integrate.o thermodynamics.o main.o
	$(COMPILADOR) $(OPTIMITZADOR) -o $@ $^

%.o : %.f90
	$(COMPILADOR) -c $^

binning : binning.o
	$(COMPILADOR) $(OPTIMITZADOR) -o $@ $^

## run: executes the main program
.PHONY : run
run: main
	./main<input.txt

## plot: runs plots from data in datafile and move them
.PHONY: plot
plot: $(datafile)
	gnuplot data2plots.gnu  
	mkdir -p Figs_redunits
	mv *.png Figs_redunits 
	gnuplot data2plots_adSI.gnu 
	mkdir -p Figs_SIunits
	mv *.png Figs_SIunits

## stats: generates statistics
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
	@echo OPTIMITZADOR: $(OPTIMITZADOR)
	@echo datafile: $(datafile)
	@echo compressed_file: $(compressed_file)
	@echo EPSS: $(EPSS)
	
## help: displays this help message
.PHONY : help
help:
	@grep '^##' Makefile


