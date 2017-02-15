#!/bin/csh

#if (  $#argv != 2) then
#  echo '1st argument is the name of the .ini file'
#  echo ' 2nd argument is a number of mesh adaptation to draw'
#  echo ' Syntax: name of the new folder'
#else

# Du/Dx
#./rm-files 
#../SRC_O/Adgfem subdomainDer.ini 
#./../Scripts/vis-tisk.sh 1 8 scalar DWR subdomainDer.ini A color DWR 
#./copy_files.sh subdomainDer-P2-250_mte5_res

# Peak 
#./rm-files 
#../SRC_O/Adgfem peak.ini 
#./../Scripts/vis-tisk.sh 0 8 scalar DWR peak.ini A color DWR
#./copy_files.sh peak_P2_primal_plus

# nonlinear 
./rm-files 
../SRC_O/Adgfem nonlinear14.ini 
./../Scripts/vis-tisk.sh 0 8 scalar DWR nonlinear14.ini A color DWR
./copy_files.sh Smaller-nonlinear_P2_primal_plus

# Convection ###########################################################x
#./rm-files 
#../SRC_O/Adgfem convection.ini 
#./../Scripts/vis-tisk.sh 0 15 scalar DWR convection.ini A color DWR 
#./copy_files.sh convection_P2_primal_5percent 

# CROSS domain problem 
#./rm-files 
#../SRC_O/Adgfem cross.ini 
#./../Scripts/vis-tisk.sh 0 3 scalar DWR cross.ini A color DWR
#./copy_files.sh cross_P1_primal_plus



#endif
