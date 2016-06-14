#!/usr/bin/gnuplot -persist
set terminal jpeg font arial 12 size 800,600
set output "result/PP_POREZI_L=8.jpg"
set grid x y
set xlabel "T"
set ylabel "PP"
plot "results/PP_POREZI_L=8.dat" using 1:2 title "spin-glass-8" with lines lt rgb "red"
