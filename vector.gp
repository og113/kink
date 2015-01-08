#gnuplot program to plot vector output from sphaleronY0.cc
#to plot from command line type gnuplot -e "f='data/.....dat'" repi.gp
#where the .... denotes the relevant filename
#plots real part only

#if you want to save directly to a file, use the following two lines of code
set term png size 1600,800
set output './pics/instantonGuess.png'

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "phi"
set xlabel "t"
set ylabel "x"
set zlabel "phi"
set grid
set hidden3d
splot f using 1:2:3 with points

pause -1
