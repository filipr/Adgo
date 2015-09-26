#!/bin/csh

if (  $#argv != 9) then
   echo 'Syntax:  scalar_prepare.sh  <mesh_min> <mesh_max> <p_min> <p_max> <T_min> <T_max> <tau_max> <tau_split_num> <CASE>'

else 


rm -f mmake
touch mmake
chmod u+x mmake

rm -f mmake1
touch mmake1
chmod u+x mmake1

rm -f mmake2
touch mmake2
chmod u+x mmake2

rm -f mmake3
touch mmake3
chmod u+x mmake3

rm -f mmake4
touch mmake4
chmod u+x mmake4


#rm -f prepare
#touch prepare
#chmod u+x prepare

rm -f Exe
touch Exe
chmod u+x Exe

rm -f ExeD
touch ExeD
chmod u+x ExeD

rm -f ExeE
touch ExeE
chmod u+x ExeE


cat << EOF > ExeD
#!/bin/csh
EOF

cat << EOF > ExeE
#!/bin/csh
EOF

cat << EOF > Exe
#!/bin/csh

rm -f tabALL.tex
touch tabALL.tex

rm -f convfile
touch convfile

rm -f orderP*.dat

EOF

rm -f template.tex
touch template.tex
cat << EOF > template.tex
\documentclass{article}
\usepackage{amssymb}
\usepackage{amsmath}
\usepackage{epsfig}
\usepackage{rotating}
%\usepackage{a4wide}
\usepackage{GGGraphics}

\hoffset=-40mm
\voffset=-10mm
\topmargin -2cm
\textheight 280mm
\textwidth 170mm

\def\Tr{{T_h}}
\def\ehX{{\|e_h\|_X}}
\def\ehJ{{\|e_h\|_J}}
\def\etaX{{\eta}}
\def\dof{{\sf dof}}
\def\Nhp{{N_{hp}}}
 \def\ehp{{e_{hp}}}
 \def\etaA{{\eta_\mathrm{A}}}
 \def\etaS{{\eta_\mathrm{S}}}
 \def\etaT{{\eta_\mathrm{T}}}
 \def\etaST{{\eta_\mathrm{ST}}}
 \def\errX{{\|e_h\|_X}}
 \def\errY{{\|e_h\|_Y}}

\begin{document}
%\thispagestyle{empty}
%\today

%\leftmargin=+0.5in
%two figures one next to one
%\hspace*{56 mm} %insert the text

\begin{sidewaystable}[h]
command: scalar.sh $argv[1] $argv[2] $argv[3] $argv[4] $argv[5] $argv[6] $argv[7] $argv[8] $argv[9]
EOF


rm -f figures-res.tex
touch figures-res.tex

rm -f figures-A.tex
touch figures-A.tex

rm -f tisk.tex
touch tisk.tex

rm -f tisk-FP.tex
touch tisk-FP.tex

rm -f tisk-FPg.tex
touch tisk-FPg.tex

rm -f tisk-FPo.tex
touch tisk-FPo.tex

rm -f tisk.gnu
touch tisk.gnu

rm -f tisk1.gnu
touch tisk1.gnu

# cat << EOF >> tisk.gnu
# set terminal postscript eps color

# EOF

 cat << EOF >> tisk1.gnu
   set terminal postscript eps color enhanced

EOF

endif
