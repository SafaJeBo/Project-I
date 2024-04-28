# Analysis of parallelization scalability

set term png

### TIME VS NPROCS ###
set output "scalability_timevsP.png"

set autoscale
set xlabel "N. processors"
set ylabel "Time (s)"
set title "Time vs Procs"

plot "PARALELIZATION.dat" u 1:2 t "Time" with lp lt 2 ps 1 

replot
unset output

### SPEEDUP ###
set output "scalability_speedupvsP.png"

set autoscale
set xlabel "N. processors"
set ylabel "Speedup=f(P)"
set title "Speedup"

plot "PARALELIZATION.dat" u 1:($3/$2) t "Time" with lp lt 1 ps 1

replot
unset output

### EFFICIENCY ###
set output "scalability_EffvsP.png"

set autoscale
set xlabel "N. processors"
set ylabel "Efficiency=f(P)"
set title "Efficiency per processor"

plot "PARALELIZATION.dat" u 1:($3/$2/$1) t "Time" with lp lt 1 ps 1

replot
unset output