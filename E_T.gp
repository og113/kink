#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
#set terminal postscript eps color enhanced size 10,5
#set output 'figure.eps';

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "E vs T"
set xlabel "T"
set ylabel "E"
plot "results/13.02.15_pi_assembled.dat" using 4:($5/18.9) with points

pause -1
