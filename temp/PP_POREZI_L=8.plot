#!/usr/bin/gnuplot -persist
set terminal jpeg font arial 12 size 800,600
set output "graph/PP_POREZI_L=8.jpg"
set grid x y
set xlabel "T"
set ylabel "PP"
plot "result/PP_POREZI_L=8" using 2:1 title "landau-wang-8" with lines lt rgb "red"