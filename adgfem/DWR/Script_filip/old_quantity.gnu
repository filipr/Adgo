set terminal postscript eps color enhanced

set size 0.5, 0.7
set size ratio -1

unset key
unset border
unset tics

#set ylabel "homogeneous Neumann BC" font ",18"
#set label 1 at 0.05,0.05 "u_D = 0" font ",24"
#set label 2 at 0.4, 0.05 "u_D = 1" font ",24"
#set label 3 at 0.7, 0.05 "u_D = 0" font ",24"

set output "quantity.eps"
p 'gnu.00'lw 2.0 w l, 'insideSubmesh.gnu' lt 1 lc 3 lw 2.0 w l

set output "init_mesh.eps"
p 'mesh-00000' w l, 'insideSubmesh.gnu' lt 1 lc 3 lw 2.0 w l



#p 'Gnu_dwr.00'lw 1.0 w l, 'insideSubmesh.gnu' lt 1 lc 3 lw 2.0 w lgnup
