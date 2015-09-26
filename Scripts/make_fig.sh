#!/bin/csh

# + 3
if (  $#argv != 8 && $#argv != 10 && $#argv != 12 && $#argv != 14) then
   echo 'Syntax:  make_fig.sh  <fig_name>  <tri> <sol> <quantity>   <type> <figType> <color/BW> <axes>  {<xmin> <xmax> <ymin> <ymax> } { <rmin> <rmax>}'
   echo ' fig_name: output eps file'
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

cp $argv[2] tri
cp $argv[3] solx

# manual setting of colors?
@ mancol = 0
  if (  $#argv == 10 || $#argv == 14) then
    @ mancol = 1
endif

#manual settings of frame
@ manbnd = 0
  if (  $#argv == 12 || $#argv == 14) then
    @ manbnd = 1
endif

cat << EOF  > movie.cfg
  $argv[4]
  $mancol
EOF

# automatic settings of colors
if($mancol == 0) then

cat << EOF  >> movie.cfg
   0.
   1. 
EOF

endif

# manual settings of colors
 if($mancol != 0) then

 if (  $#argv == 10 ) then
cat << EOF  >> movie.cfg
   $argv[9]
   $argv[10] 
EOF

else
cat << EOF  >> movie.cfg
   $argv[13]
   $argv[14] 
EOF
endif
endif

cat << EOF  >> movie.cfg
  $manbnd
EOF

if($manbnd == 0) then
cat << EOF  >> movie.cfg
   0.  1.
   0.  1.
EOF
else 
# manual settings of frame
cat << EOF  >> movie.cfg
   $argv[9]  $argv[10] 
   $argv[11]  $argv[12] 
EOF
endif

cat << EOF  >> movie.cfg
   $argv[8]
   $argv[7]
   $argv[5]
   1
   $argv[6]

   #sol 
   #sol_VS
   #solM
   #solM_VS 
   #gridHP
   #gridHP_VS
EOF


~/adgfem/SRC/dgfem2fvmx

~/adgfem/SRC/cfigx

# pstopnm -portrait -xmax 1200  GGGfig.eps 
# mv GGGfig.eps001.ppm  $argv[1].ppm

mv GGGfig.eps $argv[1]



### #cjpeg -quality 80 GGGfig.eps001.ppm > $argv[1].jpg

#rm GGGfig.eps

#

endif

