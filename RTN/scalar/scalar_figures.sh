#!/bin/csh

cat <<EOF >> ExeE
~/adgfem/Scripts/make_fig.sh Fig_${file_name}_gK.eps $dir/tri-00000 $dir/est-00000 S3  pNeu sol_VS color Ex 
EOF


if( $mi == $argv[1]) then
 cat << EOF >> tisk-FPg.tex

\vspace{2mm}

% mesh N=$mi, 
\$P_$pi\$ %, $ti-STDG
EOF
endif
cat << EOF >> tisk-FPg.tex
  \includegraphics[height=30mm]{Fig_${file_name}_gK.eps}
EOF



cat <<EOF >> ExeD
~/adgfem/Scripts/make_fig.sh Fig_${file_name}_F-1.eps $dir/tri-00000 $dir/est-00000 S9  pNeu sol_VS color Ex 
~/adgfem/Scripts/make_fig.sh Fig_${file_name}_F-2.eps $dir/tri-00000 $dir/est-00000 S10 pNeu sol_VS color Ex 
~/adgfem/Scripts/make_fig.sh Fig_${file_name}_F-3.eps $dir/tri-00000 $dir/est-00000 S11 pNeu sol_VS color Ex 
~/adgfem/Scripts/make_fig.sh Fig_${file_name}_F-4.eps $dir/tri-00000 $dir/est-00000 S12 pNeu sol_VS color Ex 

~/adgfem/Scripts/make_fig.sh Fig_${file_name}_P-1.eps $dir/tri-00000 $dir/est-00000 S13 pNeu sol_VS color Ex 
~/adgfem/Scripts/make_fig.sh Fig_${file_name}_P-2.eps $dir/tri-00000 $dir/est-00000 S14 pNeu sol_VS color Ex 
~/adgfem/Scripts/make_fig.sh Fig_${file_name}_P-3.eps $dir/tri-00000 $dir/est-00000 S15 pNeu sol_VS color Ex 
~/adgfem/Scripts/make_fig.sh Fig_${file_name}_P-4.eps $dir/tri-00000 $dir/est-00000 S16 pNeu sol_VS color Ex 

~/adgfem/Scripts/make_fig.sh Fig_${file_name}_Fo-1.eps $dir/tri-00000 $dir/est-00000 S17 pNeu sol_VS color Ex 
~/adgfem/Scripts/make_fig.sh Fig_${file_name}_Fo-2.eps $dir/tri-00000 $dir/est-00000 S18 pNeu sol_VS color Ex 
~/adgfem/Scripts/make_fig.sh Fig_${file_name}_Fo-3.eps $dir/tri-00000 $dir/est-00000 S19 pNeu sol_VS color Ex 
~/adgfem/Scripts/make_fig.sh Fig_${file_name}_Po-1.eps $dir/tri-00000 $dir/est-00000 S20 pNeu sol_VS color Ex 
~/adgfem/Scripts/make_fig.sh Fig_${file_name}_Po-2.eps $dir/tri-00000 $dir/est-00000 S21 pNeu sol_VS color Ex 
~/adgfem/Scripts/make_fig.sh Fig_${file_name}_Po-3.eps $dir/tri-00000 $dir/est-00000 S22 pNeu sol_VS color Ex 

EOF


 cat << EOF >> tisk-FP.tex
\vspace{2mm}

 mesh N=$mi, \$P_$pi\$ %, $ti-BDF

  \includegraphics[height=30mm]{Fig_${file_name}_F-1.eps}
  \includegraphics[height=30mm]{Fig_${file_name}_F-2.eps}
  \includegraphics[height=30mm]{Fig_${file_name}_F-3.eps}
  \includegraphics[height=30mm]{Fig_${file_name}_F-4.eps}

  EOC\qquad
  \includegraphics[height=30mm]{Fig_${file_name}_Fo-1.eps}
  \includegraphics[height=30mm]{Fig_${file_name}_Fo-2.eps}
  \includegraphics[height=30mm]{Fig_${file_name}_Fo-3.eps}

\vspace{2mm}

  \includegraphics[height=30mm]{Fig_${file_name}_P-1.eps}
  \includegraphics[height=30mm]{Fig_${file_name}_P-2.eps}
  \includegraphics[height=30mm]{Fig_${file_name}_P-3.eps}
  \includegraphics[height=30mm]{Fig_${file_name}_P-4.eps}

   EOC\qquad
  \includegraphics[height=30mm]{Fig_${file_name}_Po-1.eps}
  \includegraphics[height=30mm]{Fig_${file_name}_Po-2.eps}
  \includegraphics[height=30mm]{Fig_${file_name}_Po-3.eps}

\vspace{2mm}

  \hrule
EOF

 cat << EOF >> tisk-FPo.tex
\includegraphics[height=50mm]{u-cut-$Cname.eps}
EOF



#   rm -f figA-r$mi.P$pi.T$ti-$tau.gnu
#   touch figA-r$mi.P$pi.T$ti-$tau.gnu

# cat << EOF >> figA-r$mi.P$pi.T$ti-$tau.gnu
# set terminal postscript eps color
# set logscale y
# set size 0.6, 0.7

# #set xtics 2000

# #set xlabel "CPU(s)"
# set xlabel "time"
# #set xlabel "iter"
# #set ylabel "steady state reziduum"
# set ylabel "\| u - u_h \|_X"
# #set ylabel "tau"
# #set ylabel "Delta c_L"
# EOF

# cp  figA-r$mi.P$pi.T$ti-$tau.gnu  figB-r$mi.P$pi.T$ti-$tau.gnu

# cat << EOF >> figA-r$mi.P$pi.T$ti-$tau.gnu
# set output 'figA-r$mi.P$pi.T$ti-$tau.eps'
# EOF

# cat << EOF >> figB-r$mi.P$pi.T$ti-$tau.gnu
# set ylabel "tau"

# set output 'figB-r$mi.P$pi.T$ti-$tau.eps'
# EOF



# cat <<EOF > smazA
# Sca.r$mi.P$pi.T$ti.t$tau/scalar.conv 3:38  algeb
# Sca.r$mi.P$pi.T$ti.t$tau/scalar.conv 3:39  eta_S
# Sca.r$mi.P$pi.T$ti.t$tau/scalar.conv 3:40  eta_T
# Sca.r$mi.P$pi.T$ti.t$tau/scalar.conv 3:41  eta_ST
# EOF


# cat <<EOF > smazB
# Sca.r$mi.P$pi.T$ti.t$tau/scalar.conv 3:4  tau
# Sca.r$mi.P$pi.T$ti.t$tau/scalar.conv 3:28  err_L2
# Sca.r$mi.P$pi.T$ti.t$tau/scalar.conv 3:46  eta_A/eta_ST
# Sca.r$mi.P$pi.T$ti.t$tau/scalar.conv 3:47  eta_S/eta_T
# EOF


# cat << EOF >> figures-A.tex
# mesh $mi, \$P_$pi\$, $ti-BDF, \$\tau =$tau  \$

#  \includegraphics[width=80mm]{figA-r$mi.P$pi.T$ti-$tau.eps} %60mm (3 figure),45 mm (4 figs)
#  \includegraphics[width=80mm]{figB-r$mi.P$pi.T$ti-$tau.eps} %60mm (3 figure),45 mm (4 figs)

# EOF

# cat << EOF >> Exe
# #gnuplot figA-r$mi.P$pi.T$ti-$tau.gnu
# #gnuplot figB-r$mi.P$pi.T$ti-$tau.gnu

# EOF


