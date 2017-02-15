#!/bin/csh

if (  $#argv != 8 && $#argv != 9 ) then
  echo ' Syntax: vis-tisk.sh <ifig_min>  <ifig_max>   <problem>  <method> <ini_file> <adapt>  <color/BW> <type_of_table>  (<SKIP vis-tdp.sh>)'
  echo ' Creation of figure series <ifig_min> - <ifig_max> '
  echo ' problem:  gamm, cyl, shock, naca, naca-FF, step, sod, blasius scalar scalarS scalarL pulse JK ARC ARC2'
  echo ' method:  RES, HO_res, pNeu, inter'
  echo ' ini_file:  G*.ini file used for computation'
  echo ' adapt:  A for multilevel adaptation '
  echo ' color/BW: Color or BlackWhite'
  echo ' type of table: N(one)  RES  RESoT  RESo4  RESoSpap  RES-ST-book  AMA  AMApap .... (see .f90)'
  echo ' SKIP vis-tdp.sh: if Y, we do not create colour figures again'
else

set skip = no

if($#argv == 9 ) then
    if( $argv[9] == "Y") then
      set skip = yes
   endif
endif

if($skip == "no") then
 ~/adgfem/Scripts/vis-tdp.sh $argv[1]  $argv[2] $argv[3] $argv[4] $argv[6] $argv[7]
endif

cp $argv[5] smaz.ini
 set case = `echo  1 | ~/adgfem/Scripts/read_icaseAD `
 set par  = `echo  2 | ~/adgfem/Scripts/read_icaseAD `
 set Re   = `echo  3 | ~/adgfem/Scripts/read_icaseAD `
echo $case $par $Re

if($argv[3] == "ALG2") then
  ~/adgfem/Scripts/tisk_ALG.sh  $argv[1]  $argv[2]  $argv[3]  $argv[4]  $argv[5] $case $par $Re

else

 ~/adgfem/Scripts/tisk_ps.sh $argv[1]  $argv[2]  $argv[3]  $argv[4]  $argv[5] $case $par $Re $argv[8]

##### ~/zdroj/adgfem/Scripts/tisk_jpg $argv[1]  $argv[2]  $argv[3]  $argv[4] $argv[5]
endif
