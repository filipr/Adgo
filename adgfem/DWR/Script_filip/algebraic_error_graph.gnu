set terminal postscript eps color enhanced

set size 0.7, 0.7
#set size ratio -1
#unset key
#unset border

#stats 'aDWR_errors' using 4:44 
# STATS variable get the range
#set xrange[STATS_min_x*0.95:STATS_max_x*1.05]

set output "algebraicErrors1.eps"
set logscale y
#set xtics 1,10,1e8
#set xtics add (STATS_min_x, STATS_max_x)

set format y "10^{%02T}"
#set format x "10^{%01T}"  

p[][] 'aDWR_errors' u 1:4 title 'Exact' w lp, 'aDWR_errors' u 1:6 title 'S' w lp, 'aDWR_errors' u 1:8 title 'S-D' w lp, 'aDWR_errors' u 1:10 title 'A' w lp, 'aDWR_errors' u 1:11 title 'A-D' w lp

#p[][] 'aDWR_errors' u 2:4 title 'Exact' w lp, 'aDWR_errors' u 2:6 title 'S' w lp, 'aDWR_errors' u 3:8 title 'S-D' w lp, 'aDWR_errors' u 2:10 title 'A' w lp, 'aDWR_errors' u 3:11 title 'A-D' w lp


#p[][] 'aDWR_errors' u 1:4 title 'Exact' w lp, 'aDWR_errors' u 1:6 title 'S' w lp, 'aDWR_errors' u 1:7 title 'S-abs' w lp



set output "primal_algebraicErrors.eps"
set logscale y
set format y "10^{%02T}"


p[][] 'aDWR_errors' u 1:4 title 'Exact' w lp, 'aDWR_errors' u 1:6 title 'S' w lp, 'aDWR_errors' u 1:7 title 'S-Abs' w lp, 'aDWR_errors' u 1:10 title 'A' w lp, 'aDWR_errors' u 1:11 title 'A-D' w lp

#p[][] 'aDWR_errors' u 2:4 title 'Exact' w lp, 'aDWR_errors' u 2:6 title 'S' w lp, 'aDWR_errors' u 2:7 title 'S-Abs' w lp, 'aDWR_errors' u 2:10 title 'A' w lp, 'aDWR_errors' u 3:11 title 'A-D' w lp
