# gnuplot script for plotting E, U, and Kin as a function of iteration
set term png size 1200,600 

################################### TOTAL ENERGY VS BLOCK SIZE ###############################
set output 'Etot_BA.png'


set multiplot layout 1,2

# Set plot 1 parameters
set autoscale
set xlabel "Block size"
set ylabel "Total Energy Standard Deviation"
set title "Total Energy Standard Deviation vs. Block size"
set key inside bottom right 
set border lw 2
set key font ",12"
set fit quiet

# Define fitting function
f(x)= a - b * exp(-c * x) 

# Fit the function to the data
fit f(x) 'Etot_block.dat' u 1:3 via a, b, c 

# Plotting
plot 'Etot_block.dat' u 1:3 t 'Data' with p pt 14 lc rgb "red", \
         f(x) t 'Fit' w l lw 1 lc rgb "blue" lt 30

# Set plot 2 parameters
set autoscale
set xlabel "Block size"
set ylabel "Total Energy"
set title "Total Energy vs. Block size"
set key inside bottom right
set border lw 2
set key font ",12"

# Plot data with error bars
plot 'Etot_block.dat' using 1:2:3 with yerrorbars title 'Data with Error Bars'

# Reset multiplot
unset multiplot

# Save the plot
set output


###################################### POTENTIAL ENERGY VS BLOCK SIZE ##########################
set output 'Epot_BA.png'
set multiplot layout 1,2

# Set plot 1 parameters
set autoscale
set xlabel "Block size"
set ylabel "Potential Energy Standard Deviation"
set title "Potential Energy Standard Deviation vs. Block size"
set key inside bottom right
set border lw 2
set key font ",12"

# Define fitting function
f(x)= a - b * exp(-c * x)

# Fit the function to the datfit f(x) 'Epot_block.dat' u 1:3 via a, b, c
fit f(x) 'Epot_block.dat' u 1:3 via a, b, c

# Plotting
plot 'Epot_block.dat' u 1:3 t 'Data' with p pt 14 lc rgb "red", \
         f(x) t 'Fit' w l lw 1 lc rgb "blue" lt 30

# Set plot 2 parameters
set autoscale
set xlabel "Block size"
set ylabel "Potential Energy"
set title "Potential Energy vs. Block size"
set key inside bottom right
set border lw 2
set key font ",12"

# Plot data with error bars
plot 'Epot_block.dat' using 1:2:3 with yerrorbars title 'Data with Error Bars'

# Reset multiplot
unset multiplot

# Save the plot
set output


########################################## KINETIC ENERGY VS BLOCK SIZE ############################
set output 'Ekin_BA.png'
set multiplot layout 1,2

# Set plot 1 parameters
set autoscale
set xlabel "Block size"
set ylabel "Kinetic Energy Standard Deviation"
set title "Kinetic Energy Standard Deviation vs. Block size"
set key inside bottom right
set border lw 2
set key font ",12"

# Define fitting function
f(x)= a - b * exp(-c * x)

# Fit the function to the datfit f(x) 'Ekin_block.dat' u 1:3 via a, b, c
fit f(x) 'Ekin_block.dat' u 1:3 via a, b, c

# Plotting
plot 'Ekin_block.dat' u 1:3 t 'Data' with p pt 14 lc rgb "red", \
         f(x) t 'Fit' w l lw 1 lc rgb "blue" lt 30

# Set plot 2 parameters
set autoscale
set xlabel "Block size"
set ylabel "Kinetic Energy"
set title "Kinetic Energy vs. Block size"
set key inside bottom right
set border lw 2
set key font ",12"

# Plot data with error bars
plot 'Ekin_block.dat' using 1:2:3 with yerrorbars title 'Data with Error Bars'

# Reset multiplot
unset multiplot

# Save the plot
set output

################################################## MSD VS BLOCK SIZE ##############################
set output 'MSD_BA.png'
set multiplot layout 1,2

# Set plot 1 parameters
set autoscale
set xlabel "Block size"
set ylabel "MSD Standard Deviation"
set title "MSD Standard Deviation vs. Block size"
set key inside bottom right
set border lw 2
set key font ",12"

# Define fitting function
f(x)= a - b * exp(-c * x)

# Fit the function to the datfit f(x) 'Epot_block.dat' u 1:3 via a, b, c
fit f(x) 'MSD_block.dat' u 1:3 via a, b, c

# Plotting
plot 'MSD_block.dat' u 1:3 t 'Data' with p pt 14 lc rgb "red", \
         f(x) t 'Fit' w l lw 1 lc rgb "blue" lt 30

# Set plot 2 parameters
set autoscale
set xlabel "Block size"
set ylabel "MSD"
set title "MSD vs. Block size"
set key inside bottom right
set border lw 2
set key font ",12"

# Plot data with error bars
plot 'MSD_block.dat' using 1:2:3 with yerrorbars title 'Data with Error Bars'

# Reset multiplot
unset multiplot

# Save the plot
set output


############################################### PRESSURE VS BLOCK SIZE ############################
set output 'Press_BA.png'
set multiplot layout 1,2

# Set plot 1 parameters
set autoscale
set xlabel "Block size"
set ylabel "Pressure Energy Standard Deviation"
set title "Pressure Standard Deviation vs. Block size"
set key inside bottom right
set border lw 2
set key font ",12"

# Define fitting function
f(x)= a - b * exp(-c * x)

# Fit the function to the datfit f(x) 'Press_block.dat' u 1:3 via a, b, c
fit f(x) 'Press_block.dat' u 1:3 via a, b, c

# Plotting
plot 'Press_block.dat' u 1:3 t 'Data' with p pt 14 lc rgb "red", \
         f(x) t 'Fit' w l lw 1 lc rgb "blue" lt 30

# Set plot 2 parameters
set autoscale
set xlabel "Block size"
set ylabel "Pressure"
set title "Pressure vs. Block size"
set key inside bottom right
set border lw 2
set key font ",12"

# Plot data with error bars
plot 'Press_block.dat' using 1:2:3 with yerrorbars title 'Data with Error Bars'

# Reset multiplot
unset multiplot

# Save the plot
set output


######################################### TEMPERATURE VS BLOCK SIZE ###################################
set output 'Temp_BA.png'
set multiplot layout 1,2

# Set plot 1 parameters
set autoscale
set xlabel "Block size"
set ylabel "Temperature Standard Deviation"
set title "Temperature Standard Deviation vs. Block size"
set key inside bottom right
set border lw 2
set key font ",12"

# Define fitting function
f(x)= a - b * exp(-c * x)

# Fit the function to the datfit f(x) 'Temp_block.dat' u 1:3 via a, b, c
fit f(x) 'Temp_block.dat' u 1:3 via a, b, c

# Plotting
plot 'Temp_block.dat' u 1:3 t 'Data' with p pt 14 lc rgb "red", \
         f(x) t 'Fit' w l lw 1 lc rgb "blue" lt 30

# Set plot 2 parameters
set autoscale
set xlabel "Block size"
set ylabel "Temperature"
set title "Temperature vs. Block size"
set key inside bottom right
set border lw 2
set key font ",12"

# Plot data with error bars
plot 'Temp_block.dat' using 1:2:3 with yerrorbars title 'Data with Error Bars'

# Reset multiplot
unset multiplot

# Save the plot
set output

