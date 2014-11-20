#gnuplot program to plot vector output from pi.cc
#to plot from command line type gnuplot -e "f='data/.....dat'" repi.gp
#where the .... denotes the relevant filename
#plots real part only

#if you want to save directly to a file, use the following two lines of code
<<<<<<< HEAD
#set terminal postscript eps color enhanced size 10,5
#set output './pics/fig2.eps';
=======
#set term png size 1600,800
#set output './pics/fig2.png';
>>>>>>> e14eddf6d352e5753edd56b6669c6bc48c3fb11a

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "re(phi)"
set xlabel "x"
set ylabel "re(t)-im(t)"
set zlabel "re(phi)"
set grid
set hidden3d
splot f using 3:($1-$2):4 with lines

pause -1
