#!/bin/csh

if (  $#argv != 25 && $#argv != 27) then
   echo 'Syntax:  gen_ini_file.sh ini_name  grid  Qk prf dim mod par Pk IES Tk BDF iter adt alev  TOL   tol Ttim  tau  CFL  BDFt  Re tol_mx tol_mi Lq  nBC  (M  alp) '
else

  set ndim = $argv[5]
  set tol = $argv[16]
  set tol2 = `echo $tol 10. | ~/adgfem/Scripts/multi `

cat << EOF  > $argv[1].ini
2                              ** space dimension
'../Grids/$argv[2].grid'       ** file with a triangulation
$argv[3]                              ** k, Q_k approximation of curved boundary
'../Grids/profiles.$argv[4]'          ** file with curved boundary
test.sol                       ** solution file
$ndim  $argv[6]  $argv[7]                      ** number of PDEs, 4 for NS eqs,  model problem
$argv[8]                              ** k, P_k approximation of the solution
$argv[9]    $argv[10]    $argv[11]                  ** I/E/S; t-order; BDF-STDG
30  1E-9  1.00E-01  30         ** Newton: max_iter, tol1,  tol2, minimal update of C
0.0   0.0   0.0   0            ** RTN error estimates: gamma_rem, gamma_alg, gamma_lin, nu
0   20                         ** NIPG = 1, IIPG = 0, SIPG = -1; sigma
1                              ** IC: 0=G0.rsol, -1=resultsx, -2=dgm.sol
G0.rsol                        ** IC file read when type of IC = 0A
test.conv                      ** history of convergence
$argv[12]                             ** maximum number of additional iterations
$argv[14]   $argv[13]                      ** maximum level of adaptations
$argv[15]                          ** stopping criterium for steady state
$tol                           ** cD
$tol2                       ** cL
$tol                           ** cM
$argv[17]                           **  Final Time
$argv[18]                           ** if nonzero, explicit time step 
0.5                            ** if nonzero, explicit tolerance for GMRES
0.                             ** output time
0                              ** isol_start, initial value of sol* file
$argv[19]                           ** CFL number
$argv[20]                            ** tolerance for ABDF
$argv[21]                           ** Reynolds number (2-D) / viscosity (1-D)
$argv[25]                              ** number of prescribed Inlet/Outlet BC
EOF

if($ndim == 4 && $#argv == 27) then
set rho = 1.

set v1 = `echo $argv[27]  | ~/adgfem/Scripts/vel_1 `
set v2 = `echo $argv[27]  | ~/adgfem/Scripts/vel_2 `
set p  = `echo $argv[26]  | ~/adgfem/Scripts/mach `

cat << EOF  >>  $argv[1].ini
1 0 1.  $v1  $v2  $p ** Inlet,  M = $argv[26], alpha = $argv[27] 
2 0 1.  $v1  $v2  $p ** Inlet,  M = $argv[26], alpha = $argv[27] 
EOF

else if($ndim == 1 && $#argv == 25) then
  @ j=1
  while( $j <= $argv[25])
cat << EOF  >>  $argv[1].ini
$j 0  1.  ** 
EOF
    @ j++
  end

else
    echo 'BAD PARAMETERS'
    stop

endif

endif

cat << EOF  >>  $argv[1].ini
2.0   0.   0.5    0.           ** Stabilization parameters
$argv[22]    $argv[23]    $argv[24]            ** state%tol_max, state%tol_min, state%Lq 



EOF


endif

