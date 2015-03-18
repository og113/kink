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
plot "results/main_data.dat" using 5:11 with points
#plot "results/15.02.15_pi_assembled.dat" using 4:(2.0*($7-$5*$4)/26.31894507) with points, \
#	"results/13.02.15_pi_assembled.dat" using 4:(2.0*($7-$5*$4)/26.31894507) with points, \
#	"results/16.02.15_main_Tb_0.8_0.796.dat" using 5:($10/26.31894507) with points

pause -1
