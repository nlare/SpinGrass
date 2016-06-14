#!/usr/bin/gnuplot -persist
set terminal jpeg font arial 12 size 800,600
set output "graph/PP_PARISI_L=8.jpg"
set grid x y
set xlabel "T"
set ylabel "PP"
plot "result/PP_PARISI_L=8" using 2:1 title "PP_PARISI_L-8" with lines lt rgb "red"