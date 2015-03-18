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
plot "results/main_data.dat" using ($8/18.9):11 with points
#plot "data/mainAction_Tb_0.8_0.79.dat" using ($8/18.9):11 with points
#plot "results/15.02.15_pi_assembled.dat" using ($5/18.9):(2.0*($7-$5*$4)/26.31894507) with points, \
#	"results/13.02.15_pi_assembled.dat" using ($5/18.9):(2.0*($7-$5*$4)/26.31894507) with points, \
#	"results/16.02.15_main_Tb_0.8_0.796.dat" using ($8/18.9):($10/26.31894507) with points

pause -1
