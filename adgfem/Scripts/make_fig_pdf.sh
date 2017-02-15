#!/bin/csh

# + 3
if (  $#argv != 8 && $#argv != 10 && $#argv != 12 && $#argv != 14) then
   echo 'Syntax:  make_fig.sh  <fig_name>  <tri> <sol> <quantity>   <type> <figType> <color/BW> <axes>  {<xmin> <xmax> <ymin> <ymax> } { <rmin> <rmax>}'
   echo ' fig_name: output pdf file'
   echo 'tri: file with triangulation' 
   echo 'sol: file solution on tri with triangulation' 
   echo 'quantity: quantity to draw, RO, V1, M, S2, S234, ...'
   echo 'type:  type of problem scalar, scalar11 ... '
   echo 'figType:  gridHP, gridHP_VS, sol, solM, solM_VS, ...'
   echo ' color/BW:  color fig or BW figure' 
   echo 'axes:  format E = esx.y  or F= fx.y '
   echo '{<xmin> <xmax> <ymin> <ymax>}  frame of figure, obligatory, else automatic'
   echo '{<rmin> <rmax> }  range of coloring figure, obligatory, else automatic'

   
else

if($#argv == 8) then
    ~/adgfem/Scripts/make_fig.sh smaz_smaz_smaz.eps $argv[2] $argv[3] $argv[4] $argv[5] $argv[6] $argv[7] $argv[8]

else if($#argv == 10) then
    ~/adgfem/Scripts/make_fig.sh smaz_smaz_smaz.eps $argv[2] $argv[3] $argv[4] $argv[5] $argv[6] $argv[7] $argv[8] $argv[9] $argv[10]

else if($#argv == 12) then
    ~/adgfem/Scripts/make_fig.sh smaz_smaz_smaz.eps $argv[2] $argv[3] $argv[4] $argv[5] $argv[6] $argv[7] $argv[8] $argv[9] $argv[10] $argv[11] $argv[12]

else if($#argv == 14) then
    ~/adgfem/Scripts/make_fig.sh smaz_smaz_smaz.eps $argv[2] $argv[3] $argv[4] $argv[5] $argv[6] $argv[7] $argv[8] $argv[9] $argv[10] $argv[11] $argv[12] $argv[13] $argv[14] 

endif

  epstopdf smaz_smaz_smaz.eps
  mv smaz_smaz_smaz.pdf $argv[1]
  rm smaz_smaz_smaz.eps 
endif
