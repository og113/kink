#gnuplot program to plot vector output from pi.cc
#to plot from command line type gnuplot -e "f='data/.....dat'" pi.gp
#where the .... denotes the relevant filename
#plots real next to imaginary parts

#if you want to save directly to a file, use the following two lines of code
if (o eq 'pic') {
set term png size 1600,800
set output './pics/pi.png'
}

set multiplot layout 1, 2 title "phi(x,t)";
set xtics rotate
set bmargin 5

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
splot f using 3:($2-$1):4 with points

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "im(phi)"
set xlabel "x"
set ylabel "t"
set zlabel "im(phi)"
set grid
splot f using 3:($2-$1):5 with points

unset multiplot

pause -1
