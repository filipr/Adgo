#!/bin/csh

if (  $#argv != 6) then
   echo '    '
   echo 'Syntax:  make_movie  <lev_min>   <lev_max>  <problem> <adapt> <color/BW> <tisk-tdp ready?>'
   echo '   '
   echo ' problem:  gamm, cyl, shock, naca, step, sod, DMR, pulse scalar scalarL scalar11 scalar66 JK scalaBL ARC ARC2 ALG2 pNeu'
   echo ' adapt:  A for multilevel adaptation '
   echo ' color/BW: Color or BlackWhite'  
   echo ' tisk-tdp ready Y/N (Y = tisk-tdp was already performed)'
   echo ' '
   echo '   '
else


    if($argv[4] == A) then
      set tri = "triA"
      set sol = "solA"
      set exa = "exaA"
      set err = "errA"
      set est = "estA"
    else
      set tri = "tri-"
      set sol = "sol-"
      set exa = "exa-"
      set err = "err-"
      set est = "est-"
    endif

#if($argv[6] == 'N' || $argv[6] == 'n') then
#  ../Scripts/vis-tdp $argv[1] $argv[2] $argv[3] $argv[4] $argv[5]
#endif


  @ j=$argv[1]

  while( $j <= $argv[2])
    set num_file = `echo $j | ~/adgfem/Scripts/setnum `

    #echo ' Plotting gnu.* files No. ' $num_file

   if($argv[4] == A || $argv[4] == a ) then
       set tri_file = $tri$num_file
    else
       set tri_file = $tri'00000'
    endif

    set sol_file = $sol$num_file
    set exa_file = $exa$num_file
    set err_file = $err$num_file
    set est_file = $est$num_file

    #echo  $tri_file $sol_file

    if($argv[3] == 'shock') then
       set p_max = 7
       set x_max = 2

    else if($argv[3] == 'mF') then
       set p_max = 6
       set x_max = 1

    else
       echo 'other case not implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    endif

     if($argv[6] == 'N' || $argv[6] == 'n') then

      ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color Ex 1 $p_max 
     endif

     ~/adgfem/SRC_O/plotdgm  $tri_file $sol_file  > smaz

     cat << EOF  > tisk.gnu
       set terminal postscript eps color enhanced
       set size 0.8, 0.8
       set size ratio -1
    
       unset key
       set output "P_ISO-$num_file.eps"
      p  [0:$x_max][0:$x_max] 'gnu.01' w l
EOF

    gnuplot tisk.gnu

   #pstopnm -portrait -xmax 1200  P_ISO-$num_file.eps
   #cjpeg -quality 90 P_ISO-$num_file.eps001.ppm > P_ISO-$num_file.jpg

   #pstopnm -portrait -xmax 1200  hpmesh-$num_file.eps
   #cjpeg -quality 90 hpmesh-$num_file.eps001.ppm > hpmesh-$num_file.jpg
    
#    cp hpmesh-$num_file.eps Fig1.eps
#    cp hpmeshS-$num_file.eps Fig2.eps


  cat << EOF  > resA.tex
 \documentclass[10pt]{article}
 \usepackage{epsfig}
 \usepackage{rotating}
 \usepackage{amssymb}
 \usepackage{amsmath}
% \usepackage{GGGraphics}
\hoffset=-41mm
 \voffset=-8mm
 \topmargin -2cm
 %\leftmargin -25mm
 \textheight 265mm
 \textwidth 195mm

\begin{document}
 \thispagestyle{empty}
  
 \begin{center}

\includegraphics[height=     5.00cm] {hpmesh-$num_file.eps}
%\includegraphics[height=     5.00cm] {hpmeshS-$num_file.eps}
\includegraphics[height=     5.00cm] {P_ISO-$num_file.eps}

 \end{center}
 \end{document}
EOF

 #   cat << EOF  > resA.gnu
 # set terminal postscript eps color enhanced
 # set size 0.7, 0.7
 # set size ratio -1

 # set output "P_ISO-$num_file.eps"
 # p [0:2] [0:2]'P_ISO-$num_file' w l
 # EOF
 # gnuplot resA.gnu


   latex2eps resA > smaz
   if ( -r resA.epsi  ) then
     mv resA.epsi resA.eps
   endif
   #mv resA.eps fig$num_file.eps

   pstopnm -portrait -xmax 1200  resA.eps 

   #mv resA.eps001.ppm  Fig-$num_file.ppm
   #cp resA.eps001.ppm  Fig-$num_file.ppm

   cjpeg -quality 90 resA.eps001.ppm > Fig-$num_file.jpg

   @ j++
end

mencoder "mf://*.jpg" -mf fps=5 -o valel.avi -vf scale -ovc lavc -lavcopts vcodec=msmpeg4v2
#mencoder "mf://*.jpg" -mf fps=10 -o valel.avi -vf scale -ovc lavc -lavcopts vcodec=msmpeg4v2

endif
echo 'Movie was created, for replay use, e.g.:'
echo ' '
echo 'mplayer valel.avi -loop 0'
echo ' '

endif
