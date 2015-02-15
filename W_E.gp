#gnuplot program to plot W against E from mainAction.dat

#if you want to save directly to a file, use the following two lines of code
#set terminal postscript eps color enhanced size 10,5
#set output 'figure.eps';

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "W(E,N)"
set xlabel "E"
set ylabel "W"
set grid
#plot "./data/mainAction.dat" using 8:10 with points
plot "results/15.02.15_pi_assembled.dat" using ($5/18.9):(2.0*($7-$5*$4)/26.31894507) with points, \
	"results/13.02.15_pi_assembled.dat" using ($5/18.9):(2.0*($7-$5*$4)/26.31894507) with points

pause -1
