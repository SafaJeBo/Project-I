## all runs are performed consecutively in parallel. From nproc large to small. 
##Every simulation produces a (#nprocs)procs.out file containing info. on parallelization and time

time1=$(tail -n 1 1procs.out | head -n 1 | sed "s/ //g" | cut -d "=" -f2)
for n in 1 2 4 8 16 32 40
do
	file=$n"procs.out"
	time=$(tail -n 1 $file | head -n 1 | sed "s/ //g" | cut -d "=" -f2)
	echo $n $time $time1
done > PARALELIZATION.dat
