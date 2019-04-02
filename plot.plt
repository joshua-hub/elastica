!/usr/bin/gnuplot -persist

#  TITLES
set title "as10 melsom 42593249"
set xlabel "x"
set ylabel "y"

# AXES SETTINGS
#set yrange [-0.2:1]
#set xrange [-0.2:1]
#set size square
#set xtics ('-2pi' -2*pi, '-pi' -pi, 0, 'pi' pi, '2pi' 2*pi) 
#set ytics 1
#set tics scale 0.75

# LEGEND
#set key default
#set key top left
#set key outside
#set key top left

# FILE OUTPUT
set term postscript color
set output "as10-plot-melsom-42593249.ps"

plot "test" using 3:4 title "Pressure = 40" w l, \
     "test1" using 3:4 title "Pressure = 60" w l, \
     "test2" using 3:4 title "Pressure = 100" w l
