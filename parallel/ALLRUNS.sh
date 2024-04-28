## all runs are performed consecutively in parallel. From nproc large to small. 
##Every simulation produces a (#nprocs)procs.out file containing info. on parallelization and time

make run40 > 40procs.out
echo 40
make run32 > 32procs.out
echo 32
make run16 > 16procs.out
echo 16
make run8 > 8procs.out
echo 8
make run4 > 4procs.out
echo 4
make run2 > 2procs.out
echo 2
make run1 > 1procs.out
echo 1

## numbers on last line are extracted and saved ito a file as nprocs, time columns

# tac 2procs.out | nl -ba | sed -n '3p' > procsVStime.dat
# tac 4procs.out | nl -ba | sed -n '2p' >> procsVStime.dat
# tac 8procs.out | nl -ba | sed -n '2p' >> procsVStime.dat
# tac 16procs.out | nl -ba | sed -n '2p' >> procsVStime.dat
# tac 32procs.out | nl -ba | sed -n '2p' >> procsVStime.dat
# tac 40procs.out | nl -ba | sed -n '2p' >> procsVStime.dat

# # last sed prints the -(-#lastline) and i dont want it bc not cosmetic
# awk '{$1=""; sub(/^       /, ""); print}' procsVStime.dat > PARALELIZATION6.dat
# rm -f procsVStime.dat

time1=$(tail -n 2 1procs.out | head -n 1 | sed "s/ //g" | cut -d "=" -f2)
for n in 1 2 4 8 16 32 40
do
	file=$n"procs.out"
	time=$(tail -n 2 $file | head -n 1 | sed "s/ //g" | cut -d "=" -f2)
	echo $n $time $time1
done > PARALELIZATION.dat
