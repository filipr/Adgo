#!/bin/csh

# STDGM test case = 42  !!!!!!!!!!!!!!!
# STDGM test case = 45  for stiff ODE !!!!!!!!!!!!!!!

if (  $#argv != 9) then
  echo ' Syntax: runSTDGM <p_space> <tau_min>  <tau_num> <q_min> <q_max> <problem> <icase> <mesh> <figures>'
  echo ' problem:  gamm, cyl, shock, naca, naca-FF, step, sod, blasius scalar scalarL pulse JK ARC ARC2'
  echo 'figures = T/F'
else
echo 'Run STDGM'  

rm order.dat 
touch order.dat

if ( $argv[6] == scalar) then
set problem = 'sca'
else 
echo 'N-S equation are not implemented yet.'
endif

set Re = 1. #E+3
set par = 0. #parameter in 1st row in .ini file
set icase = $argv[7] 

# POLYNOMIAL APROXIMATION IN SPACE 
set p = $argv[1]
if ($icase == 1 || $icase == 42) then
set mesh = square01.str.$argv[8]

else if ($icase == 2) then
set mesh = square-11.str.$argv[8]
#else if ($icase == 45) then 
#set mesh = grid.1e
#echo "Grid - one triangle"
else 
echo "Mesh set as square01.str.$argv[8], is it the right one for this problem?" 
set mesh = square01.str.$argv[8]
endif

set ttime = 0.999999
# Final time


set steps = 110
#number of time steps per time level

#set tau = $argv[2] 

set j = 0



foreach q (`seq $argv[4] 1 $argv[5]`) 

set tau = $argv[2] 

foreach i (`seq 1 1 $argv[3]`) 

#set file_name = "ADtest_$tau" 
set file_name = "test" 
echo $file_name

## 			Syntax: AD_gen_ini_file i_name model Re mod par Ttim  ssTOL ctol     grid        Qk prf IPG cW Pk Tdis Tk tau  Ttol  S_est tol_mx tol_mi Lq madt alev maxI stopNS tol MG TOL outputtime isol nBC M  alp)
	~/adgfem/Scripts/AD_gen_ini_file_ST $file_name $problem $Re $icase $par $ttime 1E-6 1e+5 $mesh   1   0  I  20  $p STDG $q fixed $tau   RES  1.  1.  0.0  HGhp   0  $steps  rezL2 1.E-6  no 1.E-6 1. $j 8 1.E-06

echo 'tau = ' $tau 'P^q = ' $q 

	./../SRC/AAdgfem "$file_name.ini" 

#	rm "$file_name.ini"

#echo 'tau = ' $tau 

set tau = `echo $tau 0.5 | ~/adgfem/Scripts/multi `

#echo $tau

@ j = $j + 1
end
#end of foreach tau

#tau = ` echo "2 * $tau" | bc -l ` 
#echo $tau

#echo "tau = $tau"

end
#end of foreach q

@ ic = (1 + $argv[5] - $argv[4])
echo "ic $ic"
@ ic = $ic * $argv[3]
#echo "$argv[3], $argv[5], $argv[4] ic = $ic"
# ./../Scripts/vis-tisk 0 $ic scalarST $file_name.ini n color N

if ( $argv[9] == T) then
 ./../Scripts/vis-tisk 0 $ic scalarST $file_name.ini n color STDGM
else 
./../Scripts/vis-tisk 0 0 scalarST $file_name.ini n color STDGM
endif

endif
 
#../Scripts/AD_gen_ini_file ADtestST sca 1E+3 1  0.0 1.0   1E-6 1e+5 square01.str.8   1   0  I  20  2 STDG 1 fixed 0.02   RES  1E-4  1E-3  0.0  HGhp   0  10  rezL2 1E-6  no 8 1.E-06
