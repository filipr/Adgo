#!/bin/csh

if (  $#argv != 5 && $#argv != 7 && $#argv != 9  && $#argv != 11) then
   echo 'Syntax:  gen_movie.cfg.sh  <quantity>  <axes>  <color/BW> <type> <figType>'
   echo 'or'
   echo 'Syntax:  gen_movie.cfg.sh  <quantity>  <axes>  <color/BW> <type> <figType> <rmin> <rmax> ' 
   echo 'or'
   echo 'Syntax:  gen_movie.cfg.sh  <quantity>  <axes>  <color/BW> <type> <figType> <x_min>  <x_max>  <y_min>  <y_max>  '
   echo 'or'
   echo 'Syntax:  gen_movie.cfg.sh  <quantity>  <axes>  <color/BW> <type> <figType> <x_min>  <x_max>  <y_min>  <y_max>  <rmin> <rmax> ' 
   echo 'axes:  format E = esx.y  or F= fx.y '
   echo 'type:  scalar, scalar11 ... '
   echo 'figType:  gridHP, gridHP_VS'
else

  echo  ' FORMAT of gen_movie.cf.sh WAS CHANGED ! ! ! '

  cat << EOF  > movie.cfg
  $argv[1]
EOF

if (  $#argv == 5 || $#argv == 9) then
  cat << EOF  >> movie.cfg
   0    # automatic colors
   0.
   1. 
EOF
endif

if( $#argv == 7) then
  cat << EOF  >> movie.cfg
   1    # given colors
   $argv[6]
   $argv[7]
EOF
endif

if( $#argv == 11) then
  cat << EOF  >> movie.cfg
   1    # gievn colors
   $argv[10]
   $argv[11]
EOF
endif

# automatic size of the figure
if (  $#argv == 5 || $#argv == 7) then
  cat << EOF  >> movie.cfg
   0            # automatic frame
   0.   1.
   0.    1.
EOF

else
# a priori given size of figure
  cat << EOF  >> movie.cfg
   1             # manual frame
   $argv[6] $argv[7]
   $argv[8] $argv[9]
EOF

endif

 cat << EOF  >> movie.cfg
   $argv[2]
   $argv[3]
   $argv[4]
   1
   $argv[5]
   gridHP  # NOT USED
   gridHP_VS  # NOT USED
EOF


endif

