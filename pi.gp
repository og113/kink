#gnuplot program to plot vector output from pi.cc
reset
unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "phi(t,x)"
set xlabel "x"
set ylabel "t"
set zlabel "phi"
set grid
#set pm3d at s
set palette rgb 30,31,32;
mag(x,y) = sqrt(x*x + y*y)
phase(x,y) = atan(y/x)
splot "./data/piEarly00.dat" using 3:($1-$2):5:(phase($5,$6)) title 'phi(x,t)' with points palette
