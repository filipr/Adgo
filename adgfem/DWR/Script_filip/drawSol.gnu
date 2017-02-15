set terminal postscript eps color enhanced


unset ztics
unset key

#set ylabel "homogeneous Neumann BC" font ",18"
#set label 1 at 0.05,0.05 "u_D = 0" font ",24"
#set label 2 at 0.4, 0.05 "u_D = 1" font ",24"
#set label 3 at 0.7, 0.05 "u_D = 0" font ",24"

set output "primal_sol.eps"
sp 'plot_sol_primal.gnu' w l


set output "dual_sol.eps"
sp 'plot_sol.gnu' w l
