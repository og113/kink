#gnuplot program to plot theta against T from mainAction.dat

#if you want to save directly to a file, use the following two lines of code
#set terminal postscript eps color enhanced size 10,5
#set output 'figure.eps';

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "T vs theta"
set xlabel "T"
set ylabel "theta"
set grid
plot "results/main_data.dat" using 5:7 with points

pause -1
