#!/bin/csh

#if (  $#argv != 2) then
#  echo '1st argument is the name of the .ini file'
#  echo ' 2nd argument is a number of mesh adaptation to draw'
#  echo ' Syntax: name of the new folder'
#else

# Du/Dx
#./rm-files 
#../SRC_O/Adgfem subdomainDer.ini 
#./../Scripts/vis-tisk.sh 0 5 scalar DWR subdomainDer.ini A color DWR 
#./copy_files.sh subdomainDer_P1_DWRplus_PU

# Peak 
#./rm-files 
#../SRC_O/Adgfem peak.ini 
#./../Scripts/vis-tisk.sh 0 5 scalar DWR peak.ini A color DWR
#./copy_files.sh peak_Ritz_P1

##./copy_files.sh peak_ALG_P1_primal-2

# nonlinear 
./rm-files 
../SRC_O/Adgfem nonlinear14.ini 
./../Scripts/vis-tisk.sh 0 5 scalar DWR nonlinear14.ini A color DWR
./copy_files.sh SNA_P2plus_adapt_rezL2

# Convection ###########################################################x
#./rm-files 
#../SRC_O/Adgfem convection.ini 
#./../Scripts/vis-tisk.sh 0 15 scalar DWR convection.ini A color DWR 
#./copy_files.sh convection_P2_primal_5percent


# CROSS domain problem 
#./rm-files 
#../SRC_O/Adgfem cross.ini 
#./../Scripts/vis-tisk.sh 0 2 scalar DWR cross.ini A BW DWR 
#./copy_files.sh cross_P1plus_ref

#./rm-files 
#../SRC_O/Adgfem nonlinear_58.ini 
#./../Scripts/vis-tisk.sh 0 8 scalar DWR nonlinear_58.ini A BW DWR 
#./copy_files.sh nonlinSingLineP1


#endif
