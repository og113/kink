#gnuplot program to plot vector output from pi.cc
#to plot from command line type gnuplot -e "f='data/.....dat'" pic.gp
#where the .... denotes the relevant filename
#plots magnitude versus phase

#if you want to save directly to a file, use the following two lines of code
set terminal postscript eps color enhanced size 10,5
set output './pics/figure.eps';

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "phi"
set xlabel "x"
set ylabel "t"
set zlabel "sign(re(phi))*mag(phi)"
set cblabel "phase(phi)"
set grid
set palette rgb 30,31,32;
sign(x) = x/sqrt(x*x)
mag(x,y) = sqrt(x*x + y*y)
phase(x,y) = atan(y/x)
splot f using 3:($1-$2):(sign($4)*mag($4,$5)):(phase($4,$5)) title 'phi(x,t)' with points palette

pause -1
