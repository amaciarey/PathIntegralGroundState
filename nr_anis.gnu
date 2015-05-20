input='nr_vpi'

data_file=sprintf("%s.out",input)
plot_file=sprintf("%s.eps",input)

# Fitting function

f(x)=a*((1-c)*exp(-b*x*x)+c)

fit f(x) data_file u 1:2 via a,b,c

set term postscript eps color enhanced 24
set output plot_file

set logscale y
set xlabel 'r'
set ylabel 'OBDM'
set yrange [:1.1]

plot data_file u 1:(($2+$4+$6+$8+$10+$12)/a):(($3+$5+$7+$9+$11+$13)/a) w e pt 7 lc 1 t '{/Symbol r}(r,0)',\
     data_file u 1:(($2-$4+$6-$8+$10-$12)/a):(($3+$5+$7+$9+$11+$13)/a) w e pt 7 lc 3 t '{/Symbol r}(0,r)'

set term wxt
