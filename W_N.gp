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
plot "results/main_data.dat" using 9:11 with points
#plot "results/15.02.15_pi_assembled.dat" using 6:(2.0*($7-$5*$4)/26.31894507) with points, \
#	"results/13.02.15_pi_assembled.dat" using 6:(2.0*($7-$5*$4)/26.31894507) with points, \
#	"results/16.02.15_main_Tb_0.8_0.796.dat" using 9:($7/26.31894507) with points
pause -1
