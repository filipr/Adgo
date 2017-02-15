# DRAW meshes with the submesh of the dual solution and ISOlines of the solution 
# the dual_submesh.gnu file with the subgrid is needed
#
# number of the images
l = 25
n = 25

#set terminal postscript eps color enhanced
set terminal postscript eps monochrome dashed enhanced

set size 0.5, 0.5
set size ratio -1

unset key
unset border
unset tics

#set ylabel "homogeneous Neumann BC" font ",18"
#set label 1 at 0.05,0.05 "u_D = 0" font ",24"
#set label 2 at 0.4, 0.05 "u_D = 1" font ",24"
#set label 3 at 0.7, 0.05 "u_D = 0" font ",24"

outfile = 'mesh_000.eps'
meshfile = 'mesh-00000'
set output outfile

set label 1 at 0.05,0.25 "T_f" font ",18"
set label 2 at 0.85,0.75 "T_g" font ",18"
p meshfile lw 0.5 w l, 'dual_submesh.gnu' lt 1 lw 5.0 w l, 'primal_submesh.gnu' lt 1 lw 5.0 w l

unset label 

# meshes with subdomain
do for [t=l:n] {
  outfile = sprintf('mesh_%03.0f.eps',t)
  meshfile = sprintf('mesh-%05.0f', t)
  set output outfile
  name = sprintf('plus_{%01.0f}',t)
  #set title name
  
  p meshfile lw 0.5 w l, 'dual_submesh.gnu' lt 1 lw 5.0 w l, 'primal_submesh.gnu' lt 1 lw 5.0 w l
  #p meshfile lw 0.5 w l, 'dual_submesh.gnu' lt 1 lc 3 lw 2.0 w l, 'primal_submesh.gnu' lt 1 lc 3 lw 2.0 w l
}

# isolines with subdomain
do for [t=l:n] {
  outfile = sprintf('iso_%03.0f.eps',t)
  meshfile = sprintf('U_ISO-%05.0f', t)
  set output outfile
  set title 'Isolines'
  p meshfile lw 0.2 w l, 'dual_submesh.gnu' lt 1 lc 3 lw 2.0 w l, 'primal_submesh.gnu' lt 1 lc 3 lw 2.0 w l
}

# 3d solution with subdomain
#do for [t=l:n] {
#  outfile = sprintf('u3D_%03.0f.eps',t)
#  meshfile = sprintf('U_3D-%05.0f', t)
#  set output outfile
#  set title 'Isolines'
#   sp meshfile w l
#}