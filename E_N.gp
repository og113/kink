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
set title "E vs N"
set xlabel "N"
set ylabel "E"
plot "./data/mainAction.dat" using 9:8 with points

pause -1
