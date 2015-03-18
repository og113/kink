#gnuplot program to plot W against E and N from mainAction.dat

#if you want to save directly to a file, use the following two lines of code
#set terminal postscript eps color enhanced size 10,5
#set output 'figure.eps';

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "W(N,E)"
set xlabel "N"
set ylabel "E"
set zlabel "W"
set grid
splot "results/main_data.dat" using 9:8:10 with points

pause -1
