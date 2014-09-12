#gnuplot program to plot vector output from pi.cc
#to plot from command line type gnuplot -e "f='data/.....dat'" repi.gp
#where the .... denotes the relevant filename
#plots real part only

#if you want to save directly to a file, use the following two lines of code
#set terminal postscript eps color enhanced size 10,5
#set output 'figure.eps';

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "re(phi)"
set xlabel "x"
set ylabel "t"
set zlabel "re(phi)"
set grid
splot f using 3:($2-$1):5 with points

pause -1
