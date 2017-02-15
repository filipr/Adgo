#!/bin/csh

if (  $#argv != 6) then
  echo ' Syntax: gengnu.sh  <ifig>   <problem>  <method>  <icase> <par> <Re>'
  echo ' Preparing of postscripts of figure number  <ifig> '
  echo ' problem:  gamm, cyl, shock, naca, step, sod, pulse JK, scalar11, scalarBL ARC ARC2'
  echo ' method:  RES, HO_res, pNeu, inter'
  echo ' icase, par, Re :: data from *.ini file'
else

#@ i=$argv[1]

set i = $argv[1]

set problem = $argv[2]
set method = $argv[3]

set icase = $argv[4]
set par = $argv[5]
set Re = $argv[6]


if($method == meshes) then

cat << EOF  > tisk.gnu
#set terminal postscript portrait monochrome dashed "Helvetica" 14
#set terminal postscript eps 
set terminal postscript eps color enhanced

EOF


   if($problem == scalarL ) then
 
cat << EOF  >> tisk.gnu
set size 0.6, 0.9
set size ratio -1
set output 'meshA-$i.eps'
plot[] []  'mesh-$i' with lines ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'

set output 'meshB-$i.eps'
plot[-1E-1:1E-1] [-1E-1:1E-1]  'mesh-$i' with lines,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'

set output 'meshC-$i.eps'
plot[-1E-2:1E-2] [-1E-2:1E-2]  'mesh-$i' with lines ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'

set output 'meshD-$i.eps'
plot[-1E-3:1E-3] [-1E-3:1E-3]  'mesh-$i' with lines ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'

set output 'meshE-$i.eps'
plot[-1E-4:1E-4] [-1E-4:1E-4]  'mesh-$i' with lines ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'

set output 'meshF-$i.eps'
plot[-1E-5:1E-5] [-1E-5:1E-5]  'mesh-$i' with lines ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'


### method RES, residual based error estimates
set size ratio 0

EOF

cat << EOF  >> eps_to_pdf.sh
rm -f meshA-$i.pdf meshB-$i.pdf meshC-$i.pdf meshD-$i.pdf meshE-$i.pdf
epstopdf meshA-$i.eps
epstopdf meshB-$i.eps
epstopdf meshC-$i.eps
epstopdf meshD-$i.eps
epstopdf meshE-$i.eps
epstopdf meshF-$i.eps
EOF

 
  else if($problem == scalar   ) then
cat << EOF  >> tisk.gnu
set size 0.75, 1.0
set size ratio -1
set output 'meshA-$i.eps'
plot[] []  'mesh-$i' with lines  ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'

set output 'meshB-$i.eps'
plot[0:0.2] [0.:0.2]  'mesh-$i' with lines  ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'

set output 'meshC-$i.eps'
plot[0:0.2] [0.2:0.4]  'mesh-$i' with lines  ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'

set output 'meshD-$i.eps'
plot[0.:0.2] [0.8:1.0]  'mesh-$i' with lines  ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'

set output 'meshE-$i.eps'
plot[0.5:0.7] [0.:0.2]  'mesh-$i' with lines  ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'

set output 'meshF-$i.eps'
plot[-0.1:0.1] [-1.0:-0.8]  'mesh-$i' with lines  ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'


### method RES, residual based error estimates
set size ratio 0

EOF

cat << EOF  >> eps_to_pdf.sh
rm -f meshA-$i.pdf meshB-$i.pdf meshC-$i.pdf meshD-$i.pdf meshE-$i.pdf
epstopdf meshA-$i.eps
epstopdf meshB-$i.eps
epstopdf meshC-$i.eps
epstopdf meshD-$i.eps
epstopdf meshE-$i.eps
epstopdf meshF-$i.eps
EOF

  else if($problem == scalar11   ) then
cat << EOF  >> tisk.gnu
set size 0.75, 1.0
set size ratio -1
set output 'meshA-$i.eps'
plot[] []  'mesh-$i' with lines  ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'

set output 'meshB-$i.eps'
plot[-0.1:0.5] [-1.:-0.4]  'mesh-$i' with lines  ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'

set output 'meshC-$i.eps'
plot[-0.1:0.2] [-1.:-0.7]  'mesh-$i' with lines  ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'

set output 'meshD-$i.eps'
plot[0.3:0.9] [-0.3:0.3]  'mesh-$i' with lines  ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'

set output 'meshE-$i.eps'
plot[0.6:1] [0.4:0.8]  'mesh-$i' with lines  ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'

#set output 'meshF-$i.eps'
#plot[-0.1:0.1] [-1.0:-0.8]  'mesh-$i' with lines  ,'h_split-$i' t 'h-split','p_split-$i' t 'p-split'


### method RES, residual based error estimates
set size ratio 0

EOF

cat << EOF  >> eps_to_pdf.sh
rm -f meshA-$i.pdf meshB-$i.pdf meshC-$i.pdf meshD-$i.pdf meshE-$i.pdf
epstopdf meshA-$i.eps
epstopdf meshB-$i.eps
epstopdf meshC-$i.eps
epstopdf meshD-$i.eps
epstopdf meshE-$i.eps
#epstopdf meshF-$i.eps
EOF

 else


    echo 'NO other method implemented in gengnu.sh'


    endif

else


# ======================================================================
# Generace souboru pro gnuplot
# ======================================================================

cat << EOF  > tisk.gnu
#set terminal postscript portrait monochrome dashed "Helvetica" 14
#set terminal postscript eps 
set terminal postscript eps color enhanced
set size 0.5, 0.7

EOF

if($problem == gamm || $problem == cyl ) then
#  Not for step
cat << EOF  >> tisk.gnu

set size 0.8, 1.2
set output 'P_WALLS-$i.eps'
plot 'P_WALLS-$i' with lines
#set output 'M_WALLS-$i.eps'
#plot  'M_WALLS-$i' with lines
EOF

else if($problem == sod) then
#  Not for step
cat << EOF  >> tisk.gnu

set output 'P_WALLS-$i.eps'
plot 'P_WALLS-$i' with lines
set output 'R_WALLS-$i.eps'
plot  'RO_WALLS-$i' with lines
set output 'V_WALLS-$i.eps'
plot  'V_WALLS-$i' with lines
EOF

else if($problem == vocal) then
#  Not for step
cat << EOF  >> tisk.gnu

  set size 1.6, 0.4
  set size ratio -1
  set output 'p_iso-$i.eps'
  p 'P_ISO-$i'w l

  set size 0.9, 1.
  set output 'p_isoA-$i.eps'
  p [0.2:0.8] [-0.2:0.2] 'P_ISO-$i'w l

  set grid
  set size ratio 0
  set size 0.6, 1.0
  set output 'M_BC-$i.eps'
  p 'Min-$i' u 3:2  t'M inlet' lw 2.5 w l, 'Mout-$i' u 3:2 t'M outlet' lw 2.5 w l

  set output 'p_BC-$i.eps'
  p 'Pin-$i' u 3:2 t'p inlet' lw 2.5 w l, 'Pout-$i' u 3:2  t'p outlet' lw 2.5 w l
 
  set output 'V_BC-$i.eps'
  p 'Vin-$i' u 3:2 t'|v| inlet' lw 2.5 w l, 'Vout-$i' u 3:2  t'|v| outlet' lw 2.5 w l
 
EOF

else if($problem == gamm) then
### GAMM channel
cat << EOF  >> tisk.gnu
set size 1., 0.5
set view 86, 32
#set nokey
#set output 'M_3D-$i.eps'
#splot 'M_3D-$i' with lines 

set nokey
set noxtics
set noytics
set noborder



set size 1., 0.72
set size ratio -1
set output 'M_ISO-$i.eps'
plot [-1:1] [0:1] 'M_ISO-$i' with lines 
 ###  set output 'P_tr-$i.eps'
#set output 'LM-P_iso-$i.eps'
##set output 'P_ISO-$i.eps'
#plot [-1:1] [0:1] 'P_ISO-$i' with lines 

##set output 'M_ISO-det-$i.eps'
##plot [-4:4] [-4:4] 'M_ISO-$i' with lines 
##set output 'P_ISO-det-$i.eps'
##plot [-4:4] [-4:4] 'P_ISO-$i' with lines 

#set output 'M_ISO-$i.eps'
#plot [-0.5:1.5] [-1:1] 'M_ISO-$i' with lines 
#set output 'P_ISO-$i.eps'
#plot [-0.5:1.5] [-1:1] 'P_ISO-$i' with lines 

EOF

else if($problem == cyl) then
### GAMM channel
cat << EOF  >> tisk.gnu
### cylinder
set size 1., 0.7
set size ratio -1
set output 'M_ISO-$i.eps'
plot [-4:4] [-4:4] 'M_ISO-$i' with lines 
set output 'P_ISO-$i.eps'
plot [-4:4] [-4:4] 'P_ISO-$i' with lines 

#set output 'Vor-mesh.eps'
#plot [0:10] [0:10] 'gnu.1' with lines 

EOF

else if($problem == naca) then
### GAMM channel
cat << EOF  >> tisk.gnu
### cylinder
set size 0.7, 0.7
set size ratio -1
set output 'M_ISO-$i.eps'
plot [-0.3:1.5] [-0.4:0.8] 'M_ISO-$i' with lines 

set size 0.5, 0.5
set size ratio 0
set output 'P_WALLS-$i.eps'
plot [-0.05:1.05] []  'P_WALLS-$i' with lines 

set output 'PC_WALLS-$i.eps'
plot [-0.05:1.05] []  'PC_WALLS-$i' with lines 

#set output 'Vor-mesh.eps'
#plot [0:10] [0:10] 'gnu.1' with lines 

EOF

else if($problem == shock) then
### GAMM channel
cat << EOF  >> tisk.gnu
### shock
set size 1., 1.0

set output 'P_CUT-$i.eps'
plot [0:2] [0.8:1.6] 'P_CUT-$i' with lines 

set size 0.75, 1.0

set size ratio -1
set output 'P_ISO-$i.eps'
plot [0:2] [0:2] 'P_ISO-$i' with lines 


EOF

else if($problem == vortex) then
### GAMM channel
cat << EOF  >> tisk.gnu
### shock
set size 1., 1.0

set output 'M_CUT-$i.eps'
plot  'M_CUT-$i' with lines 

set size 0.75, 1.0

set size ratio -1
set output 'M_ISO-$i.eps'
plot [0:10] [0:10] 'M_ISO-$i' with lines 
set output 'mesh-$i.eps'
plot [0:10] [0:10] 'mesh-$i' with lines 


EOF

else if($problem == se1050) then
### GAMM channel
cat << EOF  >> tisk.gnu
### shock
set size 1., 1.0

set output 'P_CUT-$i.eps'
plot  'P_CUT-$i' with lines 

set size 0.75, 1.0

set size ratio -1
set output 'M_ISO-$i.eps'
plot  'M_ISO-$i' with lines 

set output 'M_ISOa-$i.eps'
plot [0.02:0.15] [-0.22:0] 'M_ISO-$i' with lines 

set output 'mesh-$i.eps'
plot  'mesh-$i' with lines 


EOF
endif


if($icase == 3) then
cat << EOF  >> tisk.gnu
# boundary layer
set output 'U_CUT-$i.eps'
 plot [0:1] [0:1.25] 'U_CUT-$i' lw 4.0 w l,(-exp(-$Re) +(-1+exp(-$Re))*(1-x) + exp(-x*$Re) )**2 t'exact'lw 4.0 w l

set output 'U_CUTa-$i.eps'
plot [0.:0.1] [0:1.25] 'U_CUT-$i' lw 4.0 w l,(-exp(-$Re) +(-1+exp(-$Re))*(1-x) + exp(-x*$Re) )**2 t'exact'lw 4.0 w l

EOF

else if($icase == 39) then
cat << EOF  >> tisk.gnu
# boundary layer
set output 'U_CUT-$i.eps'
 plot [0:1] [0:1.25] 'U_CUT-$i' lw 4.0 w l,(-exp(-$par) +(-1+exp(-$par))*(1-x) + exp(-x*$par) )**2 t'exact'lw 4.0 w l

set output 'U_CUTa-$i.eps'
plot [0.:0.1] [0:1.25] 'U_CUT-$i' lw 4.0 w l,(-exp(-$par) +(-1+exp(-$par))*(1-x) + exp(-x*$par) )**2 t'exact'lw 4.0 w l

EOF

else if($icase == 4 || $icase == 10) then
cat << EOF  >> tisk.gnu
#singul corner
set output 'U_CUT-$i.eps'
plot [0:1] [:] 'U_CUT-$i' lw 1.0 w l,2*x*x*(1-x)*(1-x)*(2*x*x)**($par/2) t'exact'lw 0.5 w l

set output 'U_CUTa-$i.eps'
plot [0:0.1] [:] 'U_CUT-$i' lw 1.0 w l,2*x*x*(1-x)*(1-x)*(2*x*x)**($par/2) t'exact'lw 0.5 w l

EOF
else if($icase == 14) then
cat << EOF  >> tisk.gnu
# [Houston, Suli, Wihler] ityp = 14
set output 'U_CUT-$i.eps'
plot [0:1] [:] 'U_CUT-$i' lw 1.0 w l,x*(1-x)*x*(1-x)*(1-2*x)*exp(-20*(2*x-1)**2) t'exact'lw 0.5 w l

set output 'U_CUTa-$i.eps'
plot [0.25:0.75] [:] 'U_CUT-$i' lw 1.0 w l,x*(1-x)*x*(1-x)*(1-2*x)*exp(-20*(2*x-1)**2) t'exact'lw 0.5 w l

EOF

else if($icase == 15) then
cat << EOF  >> tisk.gnu
## [Houston et all], ityp = 15 Lshape singular

set output 'U_CUT-$i.eps'
plot [-1:1] [:] 'U_CUT-$i' lw 1.0 w l, (abs(x*2**0.5))**(2./3) * sin(2./3*pi/4) t 'exact'

set output 'U_CUTa-$i.eps'
plot [-0.1:0.1] [:] 'U_CUT-$i' lw 1.0 w l, (abs(x*2**0.5))**(2./3) * sin(2./3*pi/4) t 'exact'

EOF

else if($icase == 17) then
cat << EOF  >> tisk.gnu
# Boundary layer
set output 'U_CUT-$i.eps'
plot [0:1] [0:1.25] 'U_CUT-$i' lw 4.0 w l, x - (exp((x-1)*$Re) - exp(-$Re))/(1. - exp(-$Re) ) t'exact'lw 4.0 w l

set output 'U_CUTa-$i.eps'
plot [0:0.999] [0:1.25] 'U_CUT-$i' lw 4.0 w l, x - (exp((x-1)*$Re) - exp(-$Re))/(1. - exp(-$Re) ) t'exact'lw 4.0 w l

EOF

else if($icase == 40) then
cat << EOF  >> tisk.gnu
# LL-domain
set output 'U_CUT-$i.eps'
plot [-1:0]  'U_CUT-$i' lw 2.0 w l, (2*x**2)**(1./3)* sin (pi/4 * 2/3) t'exact'lw 2.0 w l

set output 'U_CUTa-$i.eps'
plot [-0.1:0.] 'U_CUT-$i' lw 2.0 w l, (2*x**2)**(1./3)* sin (pi/4 * 2/3) t'exact'lw 2.0 w l

EOF

else if($icase == 40) then
cat << EOF  >> tisk.gnu
# LL-domain
set output 'U_CUT-$i.eps'
plot [-1:0]  'U_CUT-$i' lw 2.0 w l, (2*x**2)**(1./3)* sin (pi/4 * 2/3) t'exact'lw 2.0 w l

set output 'U_CUTa-$i.eps'
plot [-0.1:0.] 'U_CUT-$i' lw 2.0 w l, (2*x**2)**(1./3)* sin (pi/4 * 2/3) t'exact'lw 2.0 w l

EOF

else if($icase == 24) then
cat << EOF  >> tisk.gnu
# Barenblatt, porus media'
set output 'U_CUT-$i.eps'
plot [-6:6]  [-0.1:1.2]'U_CUT-$i' lw 2.0 w l,'BBexact0.0' lw 2.0 w l, 'BBexact0.5' lw 2.0 w l, 'BBexact1.0' lw 2.0 w l

EOF

else if($icase == 2) then

else if($icase == 64) then

else if($icase > 0) then

echo 'UNKNOWN ICASE IN GENGNU '

cat << EOF  >> tisk.gnu
# boundary layer
set output 'U_CUT-$i.eps'
 plot  'U_CUT-$i' lw 4.0 w l

set output 'U_CUTa-$i.eps'
plot  'U_CUT-$i' lw 4.0 w l
EOF

endif

 if($icase > 0 && $icase <= 62) then
cat << EOF  >> tisk.gnu
set size 0.75, 1.0
set size ratio -1
set output 'U_ISO-$i.eps'
plot 'U_ISO-$i' lw 0.25 w l     
set output 'mesh-$i.eps'
plot 'mesh-$i' with lines 


EOF

 else if ($icase == 63) then
# battery, high domain
cat << EOF  >> tisk.gnu
set size 0.35, 1.0
set size ratio -1

set output 'U_ISO-$i.eps'
plot 'U_ISO-$i' lw 0.25 w l

set output 'mesh-$i.eps'
plot 'mesh-$i' with lines 
EOF

endif


# method based graphs

if($problem == scalar || $problem == scalarL ||  $problem == scalarBL  ||  $problem == BL ||$problem == ARC || $problem == ARC2  || $problem == vortex|| $problem == scalar11 ||  $problem == scalar66 || $problem == battery || $problem == battery_sim || $problem == naca ) then
### scalar, square domain 
~/adgfem/Scripts/sqrt3_order order.dat order1.dat



if($method == pNeu ) then

cat << EOF  >> tisk.gnu
set size ratio 0
set size 0.7, 0.7
set logscale y
set grid
set xlabel "DOF^{1/3}"

set ylabel "errors and estimates L2, H1"

set output 'errs1.eps'
p 'order1.dat' u 4:11 t'L_2' lw 2.0 w lp, 'order1.dat' u 4:12 t'H_1' lw 2.0 w lp, 'order1.dat' u 4:46 t'{/Symbol h}_{tot}' lw 2.0 w lp

set ylabel "error and estimate"
set output 'errs2.eps'
p  'order1.dat' u 4:12 t'H_1' lw 2.0 w lp,'order1.dat' u 4:46 t'{/Symbol h}_{tot}' lw 2.0 w lp

set ylabel "error and estimate"
set output 'errs3.eps'
p [][1E-7:] 'order1.dat' u 4:46 t'{/Symbol h}_{tot}' lw 2.0 w lp, 'order1.dat' u 4:41 t'{/Symbol h}_{F}' lw 2.0 w lp,'order1.dat' u 4:42 t'{/Symbol h}_{rez}' lw 2.0 w lp,'order1.dat' u 4:44 t'{/Symbol h}_{pot}' lw 2.0 w lp, 'order1.dat' u 4:45 t'{/Symbol h}_{BC}' lw 2.0 w lp, 'order1.dat' u 4:47 t'{/Symbol h}_{HG}' lw 2.0 w lp

#p [][1E-7:] 'order1.dat' u 4:22 t'{/Symbol h}_{L^2}' lw 2.0 w lp, 'order1.dat' u 4:23 t'{/Symbol h}_{L^{oo}}' lw 2.0 w lp,'order1.dat' u 4:24 t'{/Symbol h}_{H^1}' lw 2.0 w lp
unset logscale

#set size 0.7,1.0
set size 0.5,0.85

#set logscale x
set xlabel " "
set ylabel " "
EOF

else if($method == RES ) then
cat << EOF  >> tisk.gnu

### method RES, residual based error estimates
set size ratio 0

set size 0.7, 0.7
set logscale y
set grid
set xlabel "DOF^{1/3}"

set ylabel "errors and estimates L2, H1"

set output 'errs1.eps'
p 'order1.dat' u 4:11 t'L_2' lw 2.0 w lp, 'order1.dat' u 4:12 t'H_1' lw 2.0 w lp, 'order1.dat' u 4:42 t'{/Symbol h}_{S}' lw 2.0 w lp

set ylabel "error and estimate"
set output 'errs2.eps'
p  'order1.dat' u 4:27 t'X' lw 2.0 w lp,'order1.dat' u 4:42 t'{/Symbol h}_{S}' lw 2.0 w lp, 'order1.dat' u 4:45 t'{/Symbol h}_{Sp}' lw 2.0 w lp

set ylabel "error and estimate"
set output 'errs3.eps'
p [][] 'order1.dat' u 4:41 t'{/Symbol h}_{A}' lw 2.0 w lp, 'order1.dat' u 4:42 t'{/Symbol h}_{S}' lw 2.0 w lp,'order1.dat' u 4:43 t'{/Symbol h}_{T}' lw 2.0 w lp,'order1.dat' u 4:44 t'{/Symbol h}_{ST}' lw 2.0 w lp

unset logscale

#set size 0.7,1.0
#set size 0.5,0.85

#set logscale x
set xlabel " "
set ylabel " "
EOF

else if($method == DWR ) then
cat << EOF  >> tisk.gnu

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
EOF

    endif

endif

endif
