COMPILADOR=gfortran
OPTIMITZADOR=-o

datafile=thermodynamics.dat
zip_init=mainmods.zip
EPSS=$(wildcard *.png)

## main: generates the main program
main : initialize.o forces.o integrate.o thermodynamics.o main.o
	$(COMPILADOR) $(OPTIMITZADOR) $@ $^

%.o : %.f90
	$(COMPILADOR) -c $^

binning : binning.o
	$(COMPILADOR) $(OPTIMITZADOR) $@ $^

## run: executes the main program
.PHONY : run
run:
	./main<input.txt

## plot: runs plots from data in datafile and move them
.PHONY: plot
plot: $(datafile)
	gnuplot data2plots.gnu  
	mkdir Figs_redunits
	mv *.png Figs_redunits 
	gnuplot data2plots_adSI.gnu 
	mkdir Figs_SIunits
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

## compress: compresses all the .f90 files and the Makefile to a ex1.tar.gz file
.PHONY : compress
compress:
	tar -czvf P1_serie.tar.gz * ../README.md 

## decompress: decompresses zip 
.PHONY: decompress
decompress: $(zip_init)
	unzip $(zip_init)

## variables: prints the variables used: the compiler and the optimitzator
.PHONY : variables
variables:
	@echo COMPILADOR: $(COMPILADOR)
	@echo OPTIMITZADOR: $(OPTIMITZADOR)
	@echo datafile: $(datafile)
	@echo zip_init: $(zip_init)
	@echo EPSS: $(EPSS)
	
## help: displays this help message
.PHONY : help
help:
	@grep '^##' Makefile


