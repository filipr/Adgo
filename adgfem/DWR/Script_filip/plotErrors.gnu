set terminal postscript eps color enhanced
set size 0.7, 0.7
set xlabel "DOF"

stats 'tableRezL2_P1' using 1:5 
# STATS variable get the range
set xrange[STATS_min_x*0.95:STATS_max_x*1.05]

set logscale y
set logscale x
set xtics 1,10,1e8

set format y "10^{%02T}"
set format x "10^{%01T}" 

set grid

set output 'errsResL2.eps'
p [][] 'tableRezL2_P1' u 1:3 t'{/Symbol h}_{E}' lw 2.0 w lp,  'tableRezL2_P1' u 1:5 t'{/Symbol h}_{S}' lw 2.0 w lp, 'tableRezL2_P1' u 1:6 t'{/Symbol h}_{A}' lw 2.0 w lp

