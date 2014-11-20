#gnuplot program to plot more than one erg line on the same plot

#if you want to save directly to a file, use the following two lines of code
set term png size 1600,800
set output './pics/linE_Na300.png';

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

#f(x) = a*sin(b*(x-c)) + d
#a = 0.5; b = 2.0; c = 0.0; d = 0.005;           # initial guess
#fit f(x) "./data/141112082242mainlinErg_0_0.dat" using ($1*17.0/199.0):2 via a, b, c, d
#title_f(a,b,c,d) = sprintf('f(x) = %.4g * sin( %.4g * (x - %.4g ) ) + %.4g ', a, b,c, d)

g(x) = i*sin(j*(x+k))*(l*x+m) + n 
i = 0.004; j = 2.0; k = 0.5; l = 0.001; m = 90; n=0.03;
FIT_LIMIT = 1e-5;
fit g(x) "./data/141113110734mainlinErg_0_0.dat" using ($1*18.0/79.0):2 via i, j, k, l, m, n
title_g(i,j,k,l,m,n) = sprintf('g(x) = %.4g * sin( %.4g * (x + %.4g ) )*( %.4g * x - %.4g ) + %.4g ', i, j, k, l, m, n)

#h(x) = m*sin(n*(x-o)) + p 
#m = 0.0002; n = 2.0; o = 1.5; p = 0.005;
#fit h(x) "./data/141111183022mainlinErg_0_7.dat" using ($1*18.0/119.0):2 via m, n, o, p
#title_h(m,n,o,p) = sprintf('h(x) = %.4g * sin( %.4g * (x - %.4g ) ) + %.4g ', m, n, o, p)

plot "./data/141113110734mainlinErg_0_0.dat" using ($1*18.0/79.0):2 title 'Na = 300' with lines, \
	g(x) title title_g(i,j,k,l,m,n)#, \
	#"./data/141112082242mainlinErg_0_0.dat" using ($1*17.0/119.0):2 title 'Na = 120' with lines, \
    #f(x) title title_f(a,b,c,d) , \
    #"./data/141110210450mainlinErg_0_0.dat" using ($1*18.0/119.0):2 title 'dE=0.05' with lines, \
    #h(x) title title_h(m,n,o,p)
