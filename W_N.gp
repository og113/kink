#gnuplot program to plot W against N from mainAction.dat

#if you want to save directly to a file, use the following two lines of code
#set terminal postscript eps color enhanced size 10,5
#set output 'figure.eps';

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "W(N)"
set xlabel "N"
set ylabel "W"
plot f using 9:10 with points

pause -1
