#!/bin/csh

if (  $#argv != 9) then
  echo ' Syntax: tisk_ps.sh <ifig_min>  <ifig_max>   <problem>  <method> <ini_file> <icase> <par> <Re>  <type_of_table>'
  echo ' Preparing of postscripts of figure series <ifig_min> - <ifig_max> '
  echo ' problem:  gamm, cyl, shock, naca, step, sod, pulse JK, scalar11, scalarBL'
  echo ' method:  RES, HO_res, pNeu, inter'
  echo ' ini_file:  G*.ini file used for computation'
  echo ' type of table: RES, RES-ST-book, AMA, AMApap, .... (see .f90)'
else

  echo 'Preparing of postscripts of figure series' $argv[1]' -  '$argv[2]

@ i=$argv[1]
@ j=$argv[2]
set method = $argv[4]
   
# ./vis-tdp $argv[1]  $argv[2] > smaz
if($method == meshes) then
 cat << EOF  > eps_to_pdf.sh
#!/bin/csh
EOF
endif

while( $i <= $argv[2])
  set num_file = `echo $i | ~/adgfem/Scripts/setnum `


  ## PS postscript figures
  #  echo 'gengnu' $num_file $argv[3] $argv[5] $argv[6] $argv[7]
  ~/adgfem/Scripts/gengnu.sh $num_file $argv[3] $argv[4] $argv[6] $argv[7] $argv[8]

  ### JPEG figures  
  #~/adgfem/Scripts/gengnu-jpg $num_file $argv[3]

  gnuplot tisk.gnu 

    if($argv[3] == shock || $argv[3] == vocal) then
    ### JPEG figures  
    ~/adgfem/Scripts/gengnu-jpg $num_file $argv[3]
     gnuplot tisk.gnu 
     #mv   P_ISO-$num_file.ps P_ISO-$num_file.ps

     #pstopnm -portrait -xmax 1200  P_ISOvalel.ps 
     #cjpeg -quality 80 valel001.ppm > $val_file
    endif

  @ i+=1
end

echo 'Postscripts of figures have been plotted'

 if($argv[4] == meshes) then
   chmod u+x ./eps_to_pdf.sh
   ./eps_to_pdf.sh

    echo 'Pdf files of meshe of figures have been plotted'
 endif




  rm -f num
cat << EOF  > num
    $argv[1]  $argv[2]
    $argv[3]
    $argv[4]
    $argv[5]
    $argv[9]
EOF


# if($argv[3] == 'scalar' || $argv[3] == 'scalar11' || $argv[3] == 'scalarL' || $argv[3] == 'scalarBL'|| $argv[3] == 'BL'  || $argv[3] == 'vortex'|| $argv[3] == 'JK' || $argv[3] == 'ARC' || $argv[3] == 'ARC2'  ) then
#   ~/adgfem/Scripts/Setorder.sh order.dat  4 $argv[8] tabA.tex
# #  ~/adgfem/Scripts/Setorderx  4
# #cp tab.tex blE-3-tab.tex
# #cat tab_beg.tex tab.tex tab_end.tex > tabA.tex
# endif

if($argv[9] == "N" || $argv[9] == "n" || $argv[9] == "none") then
    rm -f tabA.tex
    touch tabA.tex
else
  if($argv[3] == scalarST) then
 #    ~/adgfem/SRC/Setorder.sh order.dat  3 $argv[8] tabA.tex  
 #    ~/adgfem/Scripts/Setorder.sh order.dat  3 STDGM tabA.tex
 echo 'scalarST???'
 else
     ~/adgfem/Scripts/Setorder.sh order.dat  4 $argv[9] tabA.tex
  endif
endif

#  ~/adgfem/Scripts/Setorderx  4
#cp tab.tex blE-3-tab.tex
#cat tab_beg.tex tab.tex tab_end.tex > tabA.tex


#endif

#~/adgfem/Scripts/gen_tisk_tex

#~/adgfem/Scripts/gen_tisk_texAD
~/adgfem/SRC_O/gen_tisk_texAD


#latex2pdf tisk

#Filip
latex tisk  
dvipdf tisk.dvi

#latex tisk  > smaz

#dvips -q -o tisk.ps tisk 
#echo 'Final file "tisk.ps" has been created'

#Filip uncomment


echo 'Final file "tisk.pdf" has been created'

## movie

 if($argv[3] == shock || $argv[3] == gamm || $argv[3] == vocal) then
    mencoder "mf://*.jpg" -mf fps=5 -o P_ISO.avi -ovc lavc -lavcopts vcodec=mpeg4    
 ##mencoder "mf://*.jpg" -mf fps=5 -o P_ISO.avi -ovc copy -lavcopts vcodec=mpeg4    

   echo ' '
   echo ' mplayer P_ISO.avi -loop 0'
   echo ' '

 endif

endif
