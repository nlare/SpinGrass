#!/usr/bin/gnuplot -persist
set terminal jpeg font arial 12 size 800,600
set output "graph/C_v_L=4.jpg"
set grid x y
set xlabel "T"
set ylabel "C_v"
plot "result/C_v_L=4" using 2:1 title "C_v_L-4" with lines lt rgb "red"