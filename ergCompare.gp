#gnuplot program to plot more than one erg line on the same plot

#if you want to save directly to a file, use the following two lines of code
set terminal postscript eps color enhanced size 10,5
set output './pics/linE_dE.eps';

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "E(t)"
set xlabel "t"
set ylabel "E"
set key

f(x) = a*sin(b*(x-c)) + d
a = 0.0002; b = 2.0; c = 1.5; d = 0.005;           # initial guess
fit f(x) "./data/141110215218mainlinErg_0_0.dat" using ($1*18.0/99.0):2 via a, b, c, d
title_f(a,b,c,d) = sprintf('f(x) = %.4g * sin( %.4g * (x - %.4g ) ) + %.4g ', a, b,c, d)

g(x) = i*sin(j*(x-k)) + l 
i = 0.0002; j = 2.0; k = 1.5; l = 0.005;
fit g(x) "./data/141110210540mainlinErg_1_0.dat" using ($1*18.0/99.0):2 via i, j, k, l
title_g(a,b,c,d) = sprintf('g(x) = %.4g * sin( %.4g * (x - %.4g ) ) + %.4g ', i, j, k, l)

h(x) = m*sin(n*(x-o)) + p 
m = 0.0002; n = 2.0; o = 1.5; p = 0.005;
fit h(x) "./data/141110210450mainlinErg_0_0.dat" using ($1*18.0/99.0):2 via m, n, o, p
title_h(m,n,o,p) = sprintf('h(x) = %.4g * sin( %.4g * (x - %.4g ) ) + %.4g ', m, n, o, p)

plot "./data/141110215218mainlinErg_0_0.dat" using ($1*18.0/99.0):2 title 'dE=0.049' with lines, \
	"./data/141110210540mainlinErg_1_0.dat" using ($1*18.0/99.0):2 title 'dE=0.0495' with lines, \
    "./data/141110210450mainlinErg_0_0.dat" using ($1*18.0/99.0):2 title 'dE=0.05' with lines, \
    f(x) title title_f(a,b,c,d), \
    g(x) title title_g(i,j,k,l), \
    h(x) title title_h(m,n,o,p)
    
