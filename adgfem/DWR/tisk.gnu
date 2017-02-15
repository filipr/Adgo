#set terminal postscript portrait monochrome dashed "Helvetica" 14
#set terminal postscript eps 
set terminal postscript eps color enhanced
set size 0.5, 0.7

# [Houston, Suli, Wihler] ityp = 14
set output 'U_CUT-00005.eps'
plot [0:1] [:] 'U_CUT-00005' lw 1.0 w l,x*(1-x)*x*(1-x)*(1-2*x)*exp(-20*(2*x-1)**2) t'exact'lw 0.5 w l

set output 'U_CUTa-00005.eps'
plot [0.25:0.75] [:] 'U_CUT-00005' lw 1.0 w l,x*(1-x)*x*(1-x)*(1-2*x)*exp(-20*(2*x-1)**2) t'exact'lw 0.5 w l

set size 0.75, 1.0
set size ratio -1
set output 'U_ISO-00005.eps'
plot 'U_ISO-00005' lw 0.25 w l     
set output 'mesh-00005.eps'
plot 'mesh-00005' with lines 



### method RES, residual based error estimates
set size ratio 0

set size 0.7, 0.7
set xlabel "DOF"
# HP adaptation - use order1.dat
#set xlabel "DOF^{1/3}"
# h adaptation - use order.dat

stats 'order.dat' using 4:44 
# STATS variable get the range
set xrange[STATS_min_x*0.95:STATS_max_x*1.05]

set logscale y
set logscale x
set xtics 1,10,1e8
#set xtics add (STATS_min_x, STATS_max_x)

set format y "10^{%02T}"
set format x "10^{%01T}" 
#set format y "%.1t*10^{%02T}"
#set format x "%.1t*10^{%01T}"

set grid


#set ylabel "errors and estimates L2, H1, eta_S"
#set output 'errs1.eps'
#p 'order1.dat' u 4:11 t'L_2' lw 2.0 w lp, 'order1.dat' u 4:12 t'H_1' lw 2.0 w lp, 'order1.dat' u 4:42 t'{/Symbol h}_{S}' lw 2.0 w lp 

#set ylabel "abs. estimates"
#set output 'errs1.eps'
#p [][] 'order1.dat' u 4:44 t'{/Symbol h}_{E}' lw 2.0 w lp,'order1.dat' u 4:43 t'{/Symbol h}_{S_{abs}}' lw 2.0 w lp, 'order1.dat' u 4:49 t'{/Symbol h}_{S_{dual-abs}}' lw 2.0 w lp


#set ylabel "abs. estimates"
set output 'errs1.eps'
p [][] 'order.dat' u 4:44 t'{/Symbol h}_{E}' lw 2.0 w lp,'order.dat' u 4:42 t'{/Symbol h}_{S_{primal}}' lw 2.0 w lp, 'order.dat' u 4:51 t'{/Symbol h}_{S_{aver}}' lw 2.0 w lp

#set ylabel "error and estimate"
#set output 'errs2.eps'
#p  'order1.dat' u 4:27 t'X' lw 2.0 w lp,'order1.dat' u 4:42 t'{/Symbol h}_{S}' lw 2.0 w lp

#set ylabel "estimates etaS,etaSabs"
set output 'errs2.eps'
p [][] 'order.dat' u 4:44 t'{/Symbol h}_{E}' lw 2.0 w lp,'order.dat' u 4:43 t'{/Symbol h}_{S_{abs}}' lw 2.0 w lp , 'order.dat' u 4:42 t'{/Symbol h}_{S}' lw 2.0 w lp

#set ylabel "dual estimates etaS,etaSabs, etaA"
set output 'errs3.eps'
#p [][] 'order.dat' u 4:44 t'{/Symbol h}_{E}' lw 2.0 w lp,'order.dat' u 4:42 t'{/Symbol h}_{S}' lw 2.0 w lp , 'order.dat' u 4:43 t'{/Symbol h}_{S_{abs}}' lw 2.0 w lp, 'order.dat' u 4:48 t'{/Symbol h}_{S_{dual}}' lw 2.0 w lp, 'order.dat' u 4:49 t'{/Symbol h}_{S_{dual-abs}}' lw 2.0 w lp

# true + average ABS + etaA
p [][] 'order.dat' u 4:44 t'{/Symbol h}_{E}' lw 2.0 w lp,  'order.dat' u 4:51 t'{/Symbol h}_{S}' lw 2.0 w lp, 'order.dat' u 4:41 t'{/Symbol h}_{A}' lw 2.0 w lp
#p [][] 'order.dat' u 4:41 t'{/Symbol h}_{A}' lw 2.0 w lp

#set ylabel "etaA"
#set output 'errs3.eps'
#p [][] 'order1.dat' u 4:41 t'{/Symbol h}_{A}' lw 2.0 w lp


unset logscale

#set size 0.7,1.0
#set size 0.5,0.85

#set logscale x
set xlabel " "
set ylabel " "
