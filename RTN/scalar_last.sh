#!/bin/csh

cat << EOF > num
0 0 
STDGM
smaz_ini.ini 
pNeu
EOF


cat << EOF >> Exe

#~/adgfem/Scripts/gen_tisk_texAD
~/adgfem/SRC_O/gen_tisk_texAD

#./ExeD

#./ExeE

#gnuplot tisk.gnu

gnuplot tisk1.gnu
#Filip - latex2pdf does not work
#latex2pdf template
latex template
dvipdf template

EOF

cat << EOF >> template.tex

\vspace{5mm}

%\input{tabALL.tex}

\end{sidewaystable}

%\newpage

%\input{tisk-FP.tex}

%\newpage

%\input{tisk-FPg.tex}

%\newpage

\input{tisk-FPo.tex}

%\input{figures-res.tex}

%\newpage

%\input{figures-A.tex}

\end{document}

EOF

cat mmake1 mmake2 mmake3 mmake4 > mmake

