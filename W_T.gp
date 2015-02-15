#gnuplot program to plot W against T from mainAction.dat

#if you want to save directly to a file, use the following two lines of code
#set terminal postscript eps color enhanced size 10,5
#set output 'figure.eps';

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "W(T,N)"
set xlabel "T"
set ylabel "W"
set grid
#plot "./data/mainAction.dat" using 8:10 with points
plot "results/13.02.15_pi_assembled.dat" using 4:(2.0*($7-$5*$4)/26.31894507) with points

pause -1
