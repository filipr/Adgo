#!/bin/csh

# take file order.dat for AMAtdp (with unsucessfull steps)

rm -f tab.tex 

# ../../Scripts/remake_ord   tab_init 

# ../../Scripts/remake_ord DDD/Sca_m250_P1_T1_t1.00E-03/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P2_T1_t1.00E-03/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P3_T1_t1.00E-03/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P5_T1_t1.00E-03/order.dat order.dat remesh

# ../../Scripts/remake_ord   hline

# ../../Scripts/remake_ord DDD/Sca_m250_P1_T1_t1.00E-04/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P2_T1_t1.00E-04/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P3_T1_t1.00E-04/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P5_T1_t1.00E-04/order.dat order.dat remesh
# #../../Scripts/remake_ord DDD/Sca_m250_P10_T1_t1.00E-04/order.dat order.dat remesh

# ../../Scripts/remake_ord   hline
# ../../Scripts/remake_ord DDD/Sca_m250_P1_T1_t1.00E-05/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P2_T1_t1.00E-05/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P3_T1_t1.00E-05/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P5_T1_t1.00E-05/order.dat order.dat remesh
# #../../Scripts/remake_ord DDD/Sca_m250_P10_T1_t1.00E-05/order.dat order.dat remesh

# ../../Scripts/remake_ord   hline

# ../../Scripts/remake_ord DDD/Sca_m250_P1_T2_t1.00E-03/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P2_T2_t1.00E-03/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P3_T2_t1.00E-03/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P5_T2_t1.00E-03/order.dat order.dat remesh
# #../../Scripts/remake_ord DDD/Sca_m250_P10_T2_t1.00E-03/order.dat order.dat remesh

# ../../Scripts/remake_ord   hline

# ../../Scripts/remake_ord DDD/Sca_m250_P1_T2_t1.00E-04/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P2_T2_t1.00E-04/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P3_T2_t1.00E-04/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P5_T2_t1.00E-04/order.dat order.dat remesh
# #../../Scripts/remake_ord DDD/Sca_m250_P10_T2_t1.00E-04/order.dat order.dat remesh

# ../../Scripts/remake_ord   hline
# ../../Scripts/remake_ord DDD/Sca_m250_P1_T2_t1.00E-05/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P2_T2_t1.00E-05/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P3_T2_t1.00E-05/order.dat order.dat remesh
# ../../Scripts/remake_ord DDD/Sca_m250_P5_T2_t1.00E-05/order.dat order.dat remesh
# #../../Scripts/remake_ord DDD/Sca_m250_P10_T2_t1.00E-05/order.dat order.dat remesh


# ../../Scripts/remake_ord   tab_end 


# mv tab.tex tab1.tex


../../Scripts/remake_ord   tab_init 

../../Scripts/remake_ord EEE/Sca_m250_P1_T1_t1.00E-03/order.dat order.dat remesh
../../Scripts/remake_ord EEE/Sca_m250_P1_T1_t1.00E-04/order.dat order.dat remesh
../../Scripts/remake_ord EEE/Sca_m250_P1_T1_t1.00E-05/order.dat order.dat remesh

../../Scripts/remake_ord   hline

../../Scripts/remake_ord EEE/Sca_m250_P1_T2_t1.00E-03/order.dat order.dat remesh
../../Scripts/remake_ord EEE/Sca_m250_P1_T2_t1.00E-04/order.dat order.dat remesh
../../Scripts/remake_ord EEE/Sca_m250_P1_T2_t1.00E-05/order.dat order.dat remesh

../../Scripts/remake_ord   tab_end 



cat  << EOF > tiskT.tex

 \documentclass[10pt]{paper}
 \usepackage{a4wide}
 \usepackage{epsfig}
 \usepackage{rotating}
 \usepackage{amssymb}
 \usepackage{amsmath}
\def\Thpm {{{\cal T}}_{hp}^{(m)}}
\def\EIpm {{E_I }}
 \begin{document}
 \thispagestyle{empty}
  
\begin{sidewaystable}
%KNM26

%\input{tab1.tex}

\vspace{10mm}

NBK

\input{tab.tex}
\\end{sidewaystable}

\end{document}

EOF

cp tab.tex mF_tab.tex


pdflatex tiskT.tex


set Adir = EEE/Sca_m250_P1_T2_t1.00E-03
../../Scripts/remake_ord $Adir/order.dat order.dat remesh

cat  << EOF > tiskT.gnu
set terminal postscript color eps enhanced
set size 1.2, 0.65

set key spacing 1.5
set grid
set logscale y

set output "mF_L2.eps"
p [-0.02:1.02] [1E-4:1E-2] \
           'order.dat' u 10:11 t "||e_h(t_m)||_{L^2({/Symbol W})}" w lp, \
           'order.dat' u 10:31  t "||e_h ||_{L^{oo}(0,t_m;L^2({/Symbol W}))}" lw 2.5 lt 3 w l

set output "mF_Ei.eps"
p [-0.02:1.02] [1E-4:1E-2] 'order.dat' u 10:23 t "||E_I^{(m)}||_{L^2({/Symbol W})}" w lp,'remesh' u 10:23 t "remeshing" lt 3 w p


EOF
gnuplot tiskT.gnu


# ~/adgfem/Scripts/make_fig.sh mF_E-3_hp0.eps $Adir/triA00000 $Adir/solA00000 RO scalar gridHP color Ex  -1.0  1.0  -1.0  1.0 1 4
# ~/adgfem/Scripts/make_fig.sh mF_E-3_hp1.eps $Adir/triA00089 $Adir/solA00089 RO scalar gridHP color Ex  -1.0  1.0  -1.0  1.0 1 4
# ~/adgfem/Scripts/make_fig.sh mF_E-3_hp2.eps $Adir/triA00152 $Adir/solA00152 RO scalar gridHP_VS color Ex  -1.0  1.0  -1.0  1.0 1 4

set Adir = EEE/Sca_m250_P1_T2_t1.00E-05

# ~/adgfem/Scripts/make_fig.sh mF_E-5_hp0.eps $Adir/triA00010 $Adir/solA00010 RO scalar gridHP color Ex  -1.0  1.0  -1.0  1.0 1 7
# ~/adgfem/Scripts/make_fig.sh mF_E-5_hp1.eps $Adir/triA00221 $Adir/solA00221 RO scalar gridHP color Ex  -1.0  1.0  -1.0  1.0 1 7
# ~/adgfem/Scripts/make_fig.sh mF_E-5_hp2.eps $Adir/triA00429 $Adir/solA00429 RO scalar gridHP_VS color Ex  -1.0  1.0  -1.0  1.0 1 7


echo '   '
ls -l *.eps


mv mF*.eps mF*.tex ~/papers/ESCO14/

echo '  '
echo ' files moved:   mv mF*.eps mF*.tex ~/papers/ESCO14/'
echo '  '
