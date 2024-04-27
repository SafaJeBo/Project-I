## all runs are performed consecutively in parallel. From nproc large to small. 
##Every simulation produces a (#nprocs)procs.out file containing info. on parallelization and time

make run40 > 40procs.out
make run32 > 32procs.out
make run16 > 16procs.out
make run8 > 8procs.out
make run4 > 4procs.out
make run2 > 2procs.out

## numbers on last line are extracted and saved ito a file as nprocs, time columns

tac 2procs.out | nl -ba | sed -n '3p' > procsVStime.dat
tac 4procs.out | nl -ba | sed -n '2p' >> procsVStime.dat
tac 8procs.out | nl -ba | sed -n '2p' >> procsVStime.dat
tac 16procs.out | nl -ba | sed -n '2p' >> procsVStime.dat
tac 32procs.out | nl -ba | sed -n '2p' >> procsVStime.dat
tac 40procs.out | nl -ba | sed -n '2p' >> procsVStime.dat

# last sed prints the -(-#lastline) and i dont want it bc not cosmetic
awk '{$1=""; sub(/^       /, ""); print}' procsVStime.dat > PARALELIZATION.dat
rm -f procsVStime.dat

