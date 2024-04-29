# gnuplot script for plotting E, U, and Kin as a function of iteration
set term png  

################################### RDF  ###############################
set output 'gdr_serie.png'

# Set plot  parameters
set autoscale
set xlabel "r"
set ylabel "g(r)"
set title "Radial Distribution Function"
set key inside right top
set border lw 2
set key font ",12"

# Plotting
plot 'gdr_serie.dat' u 1:2  with lp pt 7 lt 10 lc rgb "red"

# Reset plot
#unset plot

# Save the plot
set output

################################### RMSD  ###############################
set output 'rmsd_serie.png'

# Set plot  parameters
set autoscale
set xlabel "Frames"
set ylabel "RMSD"
set title "Root-Mean-Square Deviation"
set key inside bottom right
set border lw 2
set key font ",12"

# Plotting
plot 'rmsd_serie.dat' u 1:2  with lp pt 20 lt 1 lc rgb "red"

# Reset plot
#unset plot

# Save the plot
set output
