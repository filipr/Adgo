#!/bin/csh



if (  $#argv != 9) then
   echo 'Syntax:  scalar.sh  <mesh_min> <mesh_max> <p_min> <p_max> <T_min> <T_max> <tau_max> <tau_split_num> <CASE>'

   echo "<CASE> = A,  case 44 - growing top (P_4 in space, exp in time)"
   echo "<CASE> = B,  case 45 - moving peak ( exp in space, passive in time)"
   echo "<CASE> = STDG case 42 only for the table"
   echo "<CASE> = ODE case 45 only for the table"
   
   echo "<CASE> = RTN case 7 only for the table"
   echo "<CASE> = RTNtime case 42 only for the table"
   echo "<CASE> = RTN2 case 2 only for the table"
   echo "<CASE> = RTN23 case 2 only for the table"
   echo "<CASE> = DUA case 7 only for the table"

   echo "<CASE> = C,  test of RES for STDG in TIME: case 49 - growing top (P_4 in space, +exp in time)"
   echo "<CASE> = D,  test of RES for STDG in SPACE: case 48 - growing top (P_4 in space, -exp in time)"
   echo "<CASE> = F,  test of RES for STDG in SPACE-TIME: case 7 - growing exp in space and time"
   echo "<CASE> = G,  test of RES for STDG in ALGEB: case 51 - growing exp in space and time"
   echo "<CASE> = H,  test of RES for STDG in singul: case 52 - growing exp in space and time"
   echo "<CASE> = I,  test of RES for STDG in singul: case 2 - moving front"

else 

#Preparation of the tex and run files / uniform for all variants
./scalar/scalar_prepare.sh $argv[1] $argv[2] $argv[3] $argv[4] $argv[5] $argv[6] $argv[7] $argv[8] $argv[9]




# Set the variable dependent on the method (mesh, icase, Re etc.)
source scalar/scalar_setVariables.sh

cat << EOF >> template.tex
method: {$IPG}IPG

\input{des.tex}
EOF

# Main cycle making the ini files and preparing their related tex files 
@ iini = 0

@ icn = 0

@ pi=$argv[3]

while( $pi <= $argv[4])

   @ ti=$argv[5]

   while( $ti <= $argv[6])

cat << EOF >> Exe
touch orderP$pi.T$ti.dat
EOF

cat << EOF >> template.tex
\vspace{5mm}

\input{tabP$pi.T$ti.tex}

EOF

      @ mi=$argv[1]

      while( $mi <= $argv[2])

         #@ cfl=$argv[7]
         set tau=$argv[7]
         #@ ratio = 10


	 ## ALGEB
         @ ratio = 2
         #@ ratio = 5

	 set tau = `echo $tau 1. | ~/adgfem/Scripts/div `
	 

#         if($ipg ==  -1) then
if($mi ==  2) then
      cat << EOF >> figures-res.tex
   EOF

   cat << EOF >> tisk.tex
   EOF
endif

cat << EOF >> figures-res.tex
mesh $mi, \$P_$pi\$, $ti-STDG

 \includegraphics[width=80mm]{fig-r$mi.P$pi.T$ti.eps}   %60mm (3 figure), 45 mm (4 figures
 \includegraphics[width=80mm]{fig1-r$mi.P$pi.T$ti.eps}   %60mm (3 figure), 45 mm (4 figures

EOF


rm -f smaz smaz1
touch smaz
touch smaz1

	 @ icfl = 1
	 @ iter = 10 

         while( $icfl <= $argv[8])



	   set icn1 = `echo $icn 3 | ~/adgfem/Scripts/moddiv `
	   @ icn1++
	   #set icn1 = 1

	   ### SIPG
	   #set sigma = `echo 20 $pi | ~/adgfem/Scripts/power `

	   #set file_name = m$mi.P$pi
	   set file_name = m{$mi}_P{$pi}
	   set mesh = $mesh_name.$mi
	   

	   if($argv[9] == 'E') then
	      set mesh = $mesh_name$mi
      endif


	   set Cname = m${mi}_P${pi}_T${ti}_t${tau}
	   set dir = Sca_${Cname}

	   # Num proc
	    cp mmake$icn1 mmake


	   echo $mi $pi $ti  $IPG $sigma $tau '#p'$icn1  '#case' $icn 'dir:' $dir

	   
	   
	~/adgfem/Scripts/AD_gen_ini_file_ST 1 $file_name $problem $Re $icase $paramet $ttime 1E-8 1e+5 $mesh   1   0  $IPG  $sigma  $pi STDG $ti $tau_type $tau $estim_space  1.  1.  0.0  RGhp   0  $Miter  $non_solver $algeb  no $lin_tol $nBC 
	
	if($iini == 0) then
	   ~/adgfem/Scripts/AD_gen_ini_file_ST 0 smaz_ini $problem $Re $icase $paramet $ttime 1E-8 1e+5 $mesh   1   0  $IPG  $sigma  $pi STDG $ti $tau_type $tau $estim_space  1.  1.  0.0  RGhp   0  $Miter  $non_solver $algeb no $lin_tol $nBC 

	   @ iini++
   endif


# preparing the script which runs Adgfem 
    cat << EOF >> mmake
    cd $dir
    ../../SRC_O/Adgfem $file_name.ini
    ../../SRC_O/plotdgm tri-00000 sol-00000
    cd ../

EOF

	   # Num proc
       cp mmake mmake$icn1 

        rm -r -f $dir/
	mkdir $dir 
	mv $file_name.ini $dir/
	cp plot.dgm $dir/
	
#preparing the script which extracts the information from the order.dat files
cat <<EOF >> Exe
  mv orderP$pi.T$ti.dat 2.dat
  tail -n 1 $dir/order.dat > 1.dat

  cat 2.dat 1.dat >  orderP$pi.T$ti.dat

    rm 2.dat 1.dat
EOF

# make_figures /make_fig.sh + input the figures to the tex files tisk***.tex
source scalar/scalar_figures.sh



           #./scalar.INI $mi $pi $ti $IPG $sigma $tau $iter

	    #cp mmake mmake$icn1 

	   @ icn++

cat <<EOF >> smaz
Sca.r$mi.P$pi.T$ti.t$tau/scalar.conv 3:30  tau=$tau
EOF
##Sca.r$mi.P$pi.T$ti.t$tau/scalar.conv r$mi.P$pi.T$ti

cat <<EOF >> smaz1
Sca.r$mi.P$pi.T$ti.t$tau/scalar.conv 3:28  tau=$tau
EOF




# plot the numerical and exact solution
source scalar/scalar_plotSolution.sh





         #./append_gnu figA-r$mi.P$pi.T$ti-$tau.gnu smazA
         #./append_gnu figB-r$mi.P$pi.T$ti-$tau.gnu smazB

           set tau = `echo $tau $ratio | ~/adgfem/Scripts/div `

	   @ icfl++
	   @ iter *= 2

         end

cat << EOF >> Exe
#gnuplot fig-r$mi.P$pi.T$ti.gnu

#gnuplot fig1-r$mi.P$pi.T$ti.gnu
EOF

if($argv[9] == 'A') then
 cat << EOF >> tisk-FPo.tex


EOF
endif

         #./append_gnu fig-r$mi.P$pi.T$ti.gnu smaz
         #./append_gnu fig1-r$mi.P$pi.T$ti.gnu smaz1


	if($argv[9] == 'E') then
	    @ mi *= 2
	else
	    @ mi *= 4
	endif
      end

      
   #mv order.dat order-P$pi.dat


echo 'Type_tab' $type_tab

# Prepare the file for computing orders of the computation
cat << EOF >> Exe
#cp orderP$pi.T$ti.dat orderP$pi.xlt

../Scripts/o_Setorder.sh orderP$pi.T$ti.dat $iord $type_tab tabP$pi.T$ti.tex

EOF

if($pi == $argv[3]) then
cat << EOF >> Exe
  cat tab_beg.tex >> tabALL.tex
EOF
endif
cat << EOF >> Exe
  cat tab.tex >> tabALL.tex
EOF
if($pi == $argv[4]) then
cat << EOF >> Exe
  cat tab_end.tex >> tabALL.tex
EOF
endif

      @ ti++
   end

if($argv[9] == 'B') then
 cat << EOF >> tisk-FPo.tex


EOF
endif

   @ pi++
   
end
#end of $pi cyklus
@ pi--

# connecting files together, calling gnuplot and latex, uniform for all variants
source scalar/scalar_last.sh

endif
