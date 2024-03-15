#gnuplot script for plotting E, U, and Kin as a function of iteration
set term postscript eps

set output 'Esvst.eps'

set autoscale
set xlabel "Time (simulation)"
set ylabel "Energy"
set title "Energy components evolution during simulation"
set xr [0:100]
set yr [-750:350]
set xtics 5000
set key vertical
set key outside

# Set line types and colors
set linetype 1 lc rgb "red" ps 1
set linetype 2 lc rgb "blue" lw 0.5
set linetype 3 lc rgb "dark-green" lw 1

# Plotting
set size ratio 1/10
plot 'thermodynamics.dat' u 1:2 t 'Ekin' with lp lt 2 ps 1, \
'thermodynamics.dat' u 1:4 t 'Etot' with lp lt 3 ps 1, \
'thermodynamics.dat' u 1:3 t 'ULJ' with points lt 1 ps 1 

replot
pause 1
unset output

#Eergy oscilation
set output 'Evst_osc.eps'

set autoscale
set xlabel "Time (simulation)"
set ylabel "Energy"
set title "Energy oscilation during simulation"
set xr [0:200]
set yr [-450:-250]

set linetype 4 lc rgb "dark-orange" 
plot 'thermodynamics.dat' u 1:4 t 'Etot' with lp lt 4 ps 0.5, \

rep
pause 10
unset output

#Tinst evolution
set output 'Tvst.eps'

set autoscale
set xlabel "Time (simulation)"
set ylabel "Temperature"
set title "Instant T evolution during simulation"
set xr [0:200]
set yr [0:2]

set linetype 5 lc rgb "dark-cyan"
plot 'thermodynamics.dat' u 1:5 t 'Tinst' with lp lt 5 ps 0.5, \

rep
pause 5
unset output

#Pressure evolution
set output 'Pvst.eps'

set autoscale
set xlabel "Time (simulation)"
set ylabel "Energy"
set title "Pressure evolution during simulation"
set xr [0:200]
set yr [0:10]

set linetype 6 lc rgb "gold"
plot 'thermodynamics.dat' u 1:7 t 'Pressure' with lp lt 7 ps 0.5, \

rep
pause 5
unset output

#TvsP
set output 'TvsP.eps'

set autoscale
set xlabel "Temprature"
set ylabel "Pressure"
set title "Temp. vs Press. in reduced units"
set xr [0:200]
set yr [0:10]

set linetype 7 lc rgb "violet" lw 1 dashtype 1 
plot 'thermodynamics.dat' u 5:7 t 'T=f(P)' with lp lt 7, \

rep
pause 5
unset output

#Radial distribution function/ pair distribution function

set output 'RDF.eps'

set autoscale
set title 'RDF for two fixed particles in a LJ fluid'
set xlabel 'r'
set ylabel'RDF'

set linetype 8 lc rgb "dark-blue" lw 1 pt 7 dt 2

plot 'resultsrdflong_def.dat' u 1:2 t'RDF' with lp lt 8 ps 0.75

rep
pause 5
unset output


