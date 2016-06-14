#!/usr/bin/gnuplot -persist
set terminal jpeg font arial 12 size 800,600
set output "graph/ENERGY_L=4.jpg"
set grid x y
set xlabel "T"
set ylabel "E"
plot "result/ENERGY_L=4" using 2:1 title "ENERGY_L-4" with lines lt rgb "red"