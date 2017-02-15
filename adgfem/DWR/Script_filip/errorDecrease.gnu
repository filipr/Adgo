set terminal postscript eps color enhanced latex

set size 0.5, 0.7
set size ratio -1
#unset key
#unset border

set output "errorDecrease.eps"
set logscale
set format y "%.2t*10^%+03T";
p 'reconstr.dat' u 4:44 title 'DWR' w lp, 'res.dat' u 4:44 title 'RES' w lp, 'plus.dat' u 4:44 title 'DWR plus' w lp 

