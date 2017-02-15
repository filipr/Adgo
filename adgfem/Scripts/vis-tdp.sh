#!/bin/csh

if (  $#argv != 6) then
   echo 'Syntax:  vis-tdp.sh   <lev_min>   <lev_max>  <problem>  <method>  <adapt> <color/BW>'
   echo ' problem:  gamm, cyl, shock, naca, step, sod, DMR, pulse scalar scalarL scalar11 scalar66 JK scalaBL ARC ARC2 ALG2 pNeu se1050'
  echo ' method:  RES, HO_res, pNeu, inter DWR'
   echo ' adapt:  A for multilevel adaptation '
   echo ' color/BW: Color or BlackWhite'  
   echo ' '
else
  rm -f res.tex
  touch res.tex

  @ j=$argv[1]
#  @ k=$argv[1]

  set problem = $argv[3]
  set method = $argv[4]
  set color = $argv[6]

    if($argv[5] == A) then
     echo 'Plotting of figures from files "triA'$argv[1]'" - triA'$argv[2]'" and "solA'$argv[1]'" - solA'$argv[2]'"'
    else
     echo 'Plotting of figures from files "tri" and  "sol'$argv[1]'" - sol'$argv[2]'"'
    endif
    
    if($argv[5] == A) then
      set tri = "triA"
      set sol = "solA"
      set exa = "exaA"
      set err = "errA"
      set est = "estA"
    else
      set tri = "tri-"
      set sol = "sol-"
      set exa = "exa-"
      set err = "err-"
      set est = "est-"
    endif

  while( $j <= $argv[2])
    set num_file = `echo $j | ~/adgfem/Scripts/setnum `

    echo ' Plotting gnu.* files No. ' $num_file

    if($argv[5] == A || $argv[5] == a ) then
       set tri_file = $tri$num_file
    else
       set tri_file = $tri'00000'
    endif

    set sol_file = $sol$num_file
    set exa_file = $exa$num_file
    set err_file = $err$num_file
    set est_file = $est$num_file

    echo  $tri_file $sol_file

   #if($problem != scalar11) then
    ~/adgfem/SRC_O/plotdgm $tri_file $sol_file > smaz
   # endif

    #mv gnu.00 hpmesh$num_file

#  Euler equations: cylinder,
    if($problem == DMR  ) then
      mv gnu.01 RO_ISO-$num_file
      #mv gnu.02 RO_cut0.25-$num_file
      #mv gnu.03 RO_cut0.50-$num_file
      #mv gnu.04 RO_cut0.75-$num_file
    endif


    if($problem == scalarBL) then
     cp gnu.00 mesh-$num_file
      cp gnu.01 U_ISO-$num_file
      cp gnu.02 U_3D-$num_file
      cp gnu.03 U_CUT-$num_file

       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.0  1.0  0.0  1.0 

       ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.9 1.0  0.9  1.0 

       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.5  0.6  0.9  1.0 

#       ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.9  1.0  0.5  0.6

        ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $est_file E6 scalar sol_VS color E 

   ~/adgfem/SRC_O/plotdgm $tri_file $est_file > smaz
    mv gnu.06 errL8-$num_file
    endif


 
    if($problem == scalar) then
      cp gnu.00 mesh-$num_file

      if($method == meshes) then


     else

      cp gnu.01 U_ISO-$num_file
      cp gnu.02 U_3D-$num_file
      cp gnu.03 U_CUT-$num_file

      # new for Filip DWR uOverSubdomain or pointValue
      # Ei ~ q(i) in OutputElementEstims in io_sub.f90 
      # change which figures (hpmeshD, U_CUT, and more) in gen_tisk_texAD.f90
       if($method == DWR) then
         cp gnu.01 U_ISO-$num_file
         cp gnu.02 U_3D-$num_file
         cp gnu.04 U_CUT-$num_file # gnu.04 for icase = 14
         
        #SQUARE 01
         #whole domain
       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0.0  1.0  0.0  1.0  0  9
         # peak
       ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0.24 0.26   0.24 0.26  0  9
         # supp of target functional ../Grids/square01.4uns.250 ../Subgrids/point0-8125in250
       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0.74  0.76  0.74 0.76  0 9  
       
       
         #SQUARE CROSS -2 to 2.0
         #whole domain
#       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  -2.0  2.0  -2.0  2.0  0  9
#         # peak
#       ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  1.2 1.4   0.2 0.4  0  9
#         # supp of target functional ../Grids/square01.4uns.250 ../Subgrids/point0-8125in250
#       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  1.0  2.0  -1.0 1.0  0 9  

        #SQUARE 01 case 14 - nonlinear, square01-point012.grid, point012 ~ 180 elems 
         #whole domain
#       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0.0  1.0  0.0  1.0  0  9
#         # bottom
#       ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0.25 0.75   0.0 0.5  0  9
#         # supp of target functional ../Grids/square01.4uns.250 ../Subgrids/point0-8125in250
#       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0.49  0.51  0.11 0.13  0 9 
       
       # LL-shaped 
#       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  -1.0  1.0  -1.0  1.0  0  9
#       # peak
#       ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  -0.6 -0.4   0.4 0.6  0  9
#         # supp of target functional ../Grids/square01.4uns.250 ../Subgrids/LL-hexagon ?
#       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0.4  0.6  -0.6 -0.4  0 9  

        #SQUARE 01 CONVECTION ICASE=16
         #whole domain
#       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0.0  1.0  0.0  1.0  0  9
#         # bottom left quarter
#       ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0.0 0.5   0.0 0.5  0  9
#         # supp of target functional ../Grids/square01.4uns.250 ../Subgrids/square0306.subgrid
#       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0.3 0.35 0.6 0.65  0 9  
       #OLD
                # supp of target functional ../Grids/square01.4uns.250 ../Subgrids/point0-8125in250
#       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0.25 0.375 0.5625 0.6875  0 9  
       
       # L2  Ei ~ q(i) in OutputElementEstims in io_sub.f90
       ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $est_file E14 $method sol_VS color E 
       # space error indicators |etaS| 
       ~/adgfem/Scripts/make_fig.sh hpmeshG-$num_file.eps $tri_file $est_file E2 $method sol_VS color E 
#       # algebraic indicators 
#       ~/adgfem/Scripts/make_fig.sh hpmeshH-$num_file.eps $tri_file $est_file E1 $method sol_VS color E 

       # space error DUAL indicators |eta_Sdual| 
       ~/adgfem/Scripts/make_fig.sh hpmeshH-$num_file.eps $tri_file $est_file E4 $method sol_VS color E 

       # space error DUAL indicators |eta_Sdual| 
#       ~/adgfem/Scripts/make_fig.sh hpmeshI-$num_file.eps $tri_file $est_file E4 $method sol_VS color E
       
       # algebraic error DUAL indicators |eta_Adual| 
#       ~/adgfem/Scripts/make_fig.sh hpmeshJ-$num_file.eps $tri_file $est_file E3 $method sol_VS color E
#       # space AVERAGE indicators |eta_aver| 
#       ~/adgfem/Scripts/make_fig.sh hpmeshJ-$num_file.eps $tri_file $est_file E5 $method sol_VS color E

       #sign indicators eta_sign
       ~/adgfem/Scripts/make_fig.sh hpmeshI-$num_file.eps $tri_file $est_file E6 $method sol_VS color E
       
       #dual sign indicators eta_sign
       ~/adgfem/Scripts/make_fig.sh hpmeshJ-$num_file.eps $tri_file $est_file E7 $method sol_VS color E
       
          #~/adgfem/Scripts/make_fig.sh hpmeshI-$num_file.eps $tri_file $est_file E17 HO_rec sol_VS color E 


      else  # scalar HERE

       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color Ex  0.0  1.0  0.0  1.0  0 9

       ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0. 0.025  0.  0.025  0  9

#       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0. 0.001  0.  0.001  0  9

       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0.5  0.51  0.  0.01  0 9 

#       ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0.  0.1  0.5  0.6

       
          # AMA -estimates
           #~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $est_file E8 $method sol_VS color E 
           #~/adgfem/Scripts/make_fig.sh hpmeshG-$num_file.eps $tri_file $est_file E9 $method sol_VS color E 

	      # pNeu error estimates
          ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $est_file E14 $method sol_VS color E 
          ~/adgfem/Scripts/make_fig.sh hpmeshG-$num_file.eps $tri_file $est_file E10 $method sol_VS color E 
          ~/adgfem/Scripts/make_fig.sh hpmeshH-$num_file.eps $tri_file $est_file E16 $method sol_VS color E 
          #~/adgfem/Scripts/make_fig.sh hpmeshI-$num_file.eps $tri_file $est_file E17 HO_rec sol_VS color E 

       endif



      ~/adgfem/SRC_O/plotdgm $tri_file $est_file > smaz
       mv gnu.06 errL8-$num_file

       endif
    endif 

   #  if($problem == pNeu) then
   #    cp gnu.00 mesh-$num_file
   #    cp gnu.01 U_ISO-$num_file
   #    cp gnu.02 U_3D-$num_file
   #    cp gnu.03 U_CUT-$num_file

   #    ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0  9
   #    ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E -0.025 0.025  -0.025  0.025  0  9

   #    ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $err_file S1 pNeu sol_VS color E 
   #    ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $est_file S5 pNeu sol_VS color E 
   #    ~/adgfem/Scripts/make_fig.sh hpmeshG-$num_file.eps $tri_file $est_file S17 pNeu sol_VS color E 


   # ~/adgfem/SRC/plotdgm $tri_file $est_file > smaz
   #  mv gnu.06 errL8-$num_file
   #  endif 
    
    if($problem == scalarST) then
      cp gnu.00 mesh-$num_file
      cp gnu.01 U_ISO-$num_file
      cp gnu.02 U_3D-$num_file
      cp gnu.03 U_CUT-$num_file

#       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.0  1.0  0.0  1.0  0  9

    endif


    if($problem == scalar11) then


      if($method == meshes) then


     else
      cp gnu.00 mesh-$num_file
      cp gnu.01 U_ISO-$num_file
      cp gnu.05 U_CUT-$num_file

    ~/adgfem/Scripts/make_fig.sh hpmesh-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color Ex  -1.2  1.0  -1.1  1.1  1  9 

#    ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color Ex  -0.1  0.4  -1.0  -0.5  1  9 

#    ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color Ex  -0.1  0.1  -1.0  -0.8  1  9 

    ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO scalar sol_VS color Ex  

    #~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $est_file E14 $method sol_VS color E 

    ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $est_file E2 $method sol_VS color Ex  -1.2  -0.8 -0.2 0.2

    ~/adgfem/Scripts/make_fig.sh hpmeshS-$num_file.eps $tri_file $est_file E3 $method sol_VS color Ex  -1.2  -0.8 -0.2 0.2

    ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $est_file E4 $method sol_VS color Ex 

    ~/adgfem/Scripts/make_fig.sh hpmeshG-$num_file.eps $tri_file $est_file E19 $method sol_VS color Ex -1.2  -0.8 -0.2 0.2

    ~/adgfem/Scripts/make_fig.sh hpmeshH-$num_file.eps $tri_file $est_file E20 $method sol_VS color Ex -1.2  -0.8 -0.2 0.2 

#    ~/adgfem/Scripts/make_fig.sh hpmeshI-$num_file.eps $tri_file $est_file S18 $method sol_VS color Ex   -1.0  1.0  -1.0  1.0 
#    ~/adgfem/Scripts/make_fig.sh hpmeshJ-$num_file.eps $tri_file $est_file S18 $method sol_VS color Ex   -0.1  0.4  -1.0  -0.5  

#   ~/adgfem/SRC/plotdgm $tri_file $est_file > smaz
#    mv gnu.06 errL8-$num_file

    endif
    endif 
    
    

    if($problem == scalar66) then
      cp gnu.00 mesh-$num_file
      cp gnu.01 U_ISO-$num_file
      cp gnu.03 U_CUT-$num_file

       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  -6.0  6.0  -6.0  6.0  

       ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  1. 5.0  1.  5.0  

#       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0. 0.001  0.  0.001  1 8

#       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.5  0.51  0.  0.01 

#       ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.  0.1  0.5  0.6

#        ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $est_file E6 scalar sol_VS color E 
#        ~/adgfem/Scripts/make_fig.sh hpmeshG-$num_file.eps $tri_file $est_file E7 scalar sol_VS color E 

#   ~/adgfem/SRC/plotdgm $tri_file $est_file > smaz
#    mv gnu.06 errL8-$num_file
    endif
    
    
    if($problem == ARC) then
      cp gnu.00 mesh-$num_file
      cp gnu.01 U_ISO-$num_file
      cp gnu.02 U_3D-$num_file
      cp gnu.03 U_CUT-$num_file

       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.0  1.0  0.0  1.0 1 7

       ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.22 0.26  0.22  0.26 1 7 

       #~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.46 0.49  0.46  0.49 1 7 

       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.  0.1  0.3  0.4 1 7

       ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.  0.1  0.62  0.72 1 7 

        ~/adgfem/Scripts/make_fig.sh hpmeshG-$num_file.eps $tri_file $est_file E8 scalar sol_VS color E 
        ~/adgfem/Scripts/make_fig.sh hpmeshH-$num_file.eps $tri_file $est_file E10 scalar sol_VS color E 

   ~/adgfem/SRC_O/plotdgm $tri_file $est_file > smaz
    mv gnu.06 errL8-$num_file
    endif


    if($problem == ARC2) then
      cp gnu.00 mesh-$num_file
      cp gnu.01 U_ISO-$num_file
      cp gnu.02 U_3D-$num_file
      cp gnu.04 U_CUT-$num_file

       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.0  2.0  0.0  1.0 

       ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.8 1.2  0.  0.2  1 9

       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  1.6 2.0  0.8  1.0 1 9 

       ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.5  0.51  0.  0.01  1 9

       #~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.  0.1  0.5  0.6 1 9

        ~/adgfem/Scripts/make_fig.sh hpmeshG-$num_file.eps $tri_file $est_file E8 scalar sol_VS color E 0. 2. 0. 1. 
        ~/adgfem/Scripts/make_fig.sh hpmeshH-$num_file.eps $tri_file $est_file E9 scalar sol_VS color E  0. 2. 0. 1.

   ~/adgfem/SRC_O/plotdgm $tri_file $est_file > smaz
    mv gnu.06 errL8-$num_file
    endif


    if($problem == BL) then
     cp gnu.00 mesh-$num_file
      cp gnu.01 U_ISO-$num_file
      cp gnu.02 U_3D-$num_file
      cp gnu.03 U_CUT-$num_file

       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.0  1.0  0.0  1.0 

       ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.95 1.0  0.4999999 0.55  

       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color F  0.995 1.0  0.4999999 0.505  

       ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color F  0.9995 1.0  0.4999999 0.5005  


#        ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $est_file E6 scalar sol_VS color# E 

   ~/adgfem/SRC_O/plotdgm $tri_file $est_file > smaz
    mv gnu.06 errL8-$num_file
    endif



    if( $problem == scalarL ) then
      cp gnu.00 mesh-$num_file
      cp gnu.01 U_ISO-$num_file
      #cp gnu.02 U_3D-$num_file
      cp gnu.03 U_CUT-$num_file

      if($method == meshes) then


    else
       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color Ex  -1.0  1.0  -1.0  1.0 1 9

       ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color Ex  -0.05  0.05 -0.05  0.05 1 9

       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color Ex  -0.001  0.001 -0.001  0.001 1 9

       ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  -0.0001  0.0001 -0.0001  0.0001 1 9 

    ~/adgfem/Scripts/make_fig.sh hpmeshG-$num_file.eps $tri_file $est_file E10 $method sol_VS color Ex 
#-0.05  0.05 -0.05  0.05

    ~/adgfem/Scripts/make_fig.sh hpmeshH-$num_file.eps $tri_file $est_file E14 $method sol_VS color Ex 
#-0.05  0.05 -0.05  0.05

    #~/adgfem/Scripts/make_fig.sh hpmeshI-$num_file.eps $tri_file $est_file S12 $method sol_VS color Ex 0. 1. 
#-0.05  0.05 -0.05  0.05

    ~/adgfem/Scripts/make_fig.sh hpmeshI-$num_file.eps $tri_file $est_file S18 $method sol_VS color Ex -1. 1. -1. 1.

    ~/adgfem/Scripts/make_fig.sh hpmeshJ-$num_file.eps $tri_file $est_file S18 $method sol_VS color Ex -0.05  0.05 -0.05  0.05 

   #~/adgfem/SRC/plotdgm $tri_file $est_file > smaz
   # mv gnu.06 errL8-$num_file

    endif

endif


    if( $problem == pNeu ) then
      cp gnu.00 mesh-$num_file
      cp gnu.01 U_ISO-$num_file
      cp gnu.02 U_3D-$num_file
      cp gnu.03 U_CUT-$num_file

       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  

       ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0. 0.1 0. 0.1 

       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0 0.01 0. 0.01

       ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.25 0.75 0.25 0.75

    ~/adgfem/Scripts/make_fig.sh hpmeshG-$num_file.eps $tri_file $est_file S7 pNeuA sol_VS color E

    ~/adgfem/Scripts/make_fig.sh hpmeshH-$num_file.eps $tri_file $est_file S6 pNeuA sol_VS color E  

    ~/adgfem/Scripts/make_fig.sh hpmeshI-$num_file.eps $tri_file $est_file E9 pNeuA sol_VS color E 

    ~/adgfem/Scripts/make_fig.sh hpmeshJ-$num_file.eps $tri_file $est_file E10 pNeuA sol_VS color E 

   #~/adgfem/SRC/plotdgm $tri_file $est_file > smaz
   # mv gnu.06 errL8-$num_file

endif


    if( $problem == JK ) then
      cp gnu.00 mesh-$num_file
      cp gnu.01 U_ISO-$num_file
      cp gnu.02 U_3D-$num_file
      cp gnu.03 U_CUT-$num_file

       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.0  1.0  0.0  1.0 

       ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.2  0.3 0.2  0.3 

       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  0.255  0.265 0.24  0.25 

 #      ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color E  -0.0005  0.0005 -0.0005  0.0005 

    ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $est_file E6 scalar sol_VS color E 
   #~/adgfem/SRC/plotdgm $tri_file $est_file > smaz
   # mv gnu.06 errL8-$num_file

endif



    if($problem == shishkin  ) then
      mv gnu.01 u-FVM-$num_file

    endif

    if($problem == vocal  ) then
      mv gnu.00 mesh-$num_file
      mv gnu.02 P_ISO-$num_file
      mv gnu.05 Min-$num_file
      mv gnu.06 Mout-$num_file
      mv gnu.07 Pin-$num_file
      mv gnu.08 Pout-$num_file
      mv gnu.09 Vin-$num_file
      mv gnu.10 Vout-$num_file
    endif

    if($problem == shock  ) then
      mv gnu.00 mesh-$num_file
      mv gnu.01 P_ISO-$num_file
      mv gnu.02 P_CUT-$num_file

       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO scalar gridHP_VS color Ex  0.  2.0  0.0  2.0  1 7

       ~/adgfem/Scripts/make_fig.sh  hpmeshS-$num_file.eps $tri_file $sol_file P scalar sol_VS color Ex  0.0  2.0  0.0  2.0 0.9 1.5

    endif


    if($problem == gamm || $problem == cyl) then
       mv gnu.01 M_ISO-$num_file
    #  Not for step
      #mv gnu.02 M_WALLS-$num_file
      #mv gnu.03 P_ISO-$num_file
      mv gnu.02 P_WALLS-$num_file

      #mv gnu.06 RO_ISO-$num_file

       ~/adgfem/Scripts/make_fig.sh hpmesh-$num_file.eps $tri_file $sol_file M naca gridHP_VS color E -1. 1. 0. 1. 1 7

    endif

    if($problem == battery) then
      mv gnu.00 mesh-$num_file
      mv gnu.01 U_ISO-$num_file
      mv gnu.03 U_CUT-$num_file

	~/adgfem/Scripts/make_fig.sh hpmesh-$num_file.eps $tri_file $sol_file RO battery gridHP_VS color N 0. 8.4 0. 24.   0 9
	~/adgfem/Scripts/make_fig.sh U_sol-$num_file.eps $tri_file $sol_file RO battery sol_VS color N  0. 8.4 0. 24. 

        if($method == RES) then
   	~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $est_file E10  battery sol_VS color N  0. 8.4 0. 24. 
        else if($method == pNeu) then
   	~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $est_file E1  battery sol_VS color N  0. 8.4 0. 24. 
        endif
#-6 +3  # 1E-6  1E+3

    endif

    if($problem == battery_sim) then
      cp gnu.00 mesh-$num_file

      if($method == meshes) then


     else

      cp gnu.01 U_ISO-$num_file
      cp gnu.02 U_3D-$num_file
      cp gnu.03 U_CUT-$num_file


       ~/adgfem/Scripts/make_fig.sh  hpmesh-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0.0  1.0  0.0  1.0  0  9

       ~/adgfem/Scripts/make_fig.sh  hpmeshD-$num_file.eps $tri_file $est_file S18 $method sol_VS color E  0.0  1.0  0.0  1.0  

       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file RO $method sol_VS color E  

#       ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $sol_file RO $method gridHP_VS color E  0.  0.1  0.5  0.6

       
          # AMA -estimates
           #~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $est_file E8 $method sol_VS color E 
           #~/adgfem/Scripts/make_fig.sh hpmeshG-$num_file.eps $tri_file $est_file E9 $method sol_VS color E 

	      # pNeu error estimates
          ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $est_file E10 $method sol_VS color E 
          #~/adgfem/Scripts/make_fig.sh hpmeshG-$num_file.eps $tri_file $est_file E11 $method sol_VS color E 
#          ~/adgfem/Scripts/make_fig.sh hpmeshG-$num_file.eps $tri_file $est_file E10 $method sol_VS color E 
#          ~/adgfem/Scripts/make_fig.sh hpmeshH-$num_file.eps $tri_file $est_file E16 $method sol_VS color E 
          #~/adgfem/Scripts/make_fig.sh hpmeshI-$num_file.eps $tri_file $est_file E17 HO_rec sol_VS color E 


      ~/adgfem/SRC_O/plotdgm $tri_file $est_file > smaz
       mv gnu.06 errL8-$num_file

       endif
    endif 


    if($problem == naca) then
       mv gnu.00 mesh-$num_file
       mv gnu.01 M_ISO-$num_file
       mv gnu.02 P_WALLS-$num_file
       mv gnu.03 PC_WALLS-$num_file

       ~/adgfem/Scripts/make_fig.sh hpmesh-$num_file.eps $tri_file $sol_file M naca gridHP_VS color E -0.3   1.5 -0.4   0.8 1 8

       ~/adgfem/Scripts/make_fig.sh hpmeshD-$num_file.eps $tri_file $sol_file M naca gridHP_VS color E -0.07  0.07 -0.07   0.07  1 8  

       #~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file M naca gridHP_VS color E 0.85   1.15 -0.15   0.15 1 8

       #~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file M naca gridHP_VS color E -0.5 2.5 -1.5 1.5  1 8

       ~/adgfem/Scripts/make_fig.sh hpmeshE-$num_file.eps $tri_file $sol_file M naca gridHP_VS color E  0.93 1.07 -0.07 0.07  1 8

       ~/adgfem/Scripts/make_fig.sh hpmeshT-$num_file.eps $tri_file $sol_file M naca gridHP_VS color E 1 8

        ~/adgfem/Scripts/make_fig.sh hpmeshF-$num_file.eps $tri_file $est_file E10 naca sol_VS color E -4   6 -5.   5.

        ~/adgfem/Scripts/make_fig.sh hpmeshG-$num_file.eps $tri_file $est_file E10 naca sol_VS color E -0.2   1. -0.6   0.6

    endif

    if($problem == naca-FF) then
       mv gnu.01 M_ISO-$num_file
       mv gnu.02 P_WALLS-$num_file
    endif


    if($problem == sod) then
    #  Not for step
      mv gnu.01 M_ISO-$num_file
      mv gnu.04 P_WALLS-$num_file
      mv gnu.07 V_WALLS-$num_file
      mv gnu.09 RO_WALLS-$num_file
    endif

    if($problem == vortex) then
    #  Not for step
      mv gnu.00 mesh-$num_file
      mv gnu.01 M_ISO-$num_file
      mv gnu.02 M_CUT-$num_file
       ~/adgfem/Scripts/make_fig.sh hpmesh-$num_file.eps $tri_file $sol_file M naca gridHP_VS color E
    endif

    if($problem == se1050) then
    #  Not for step
      mv gnu.00 mesh-$num_file
      mv gnu.01 M_ISO-$num_file
      mv gnu.04 P_CUT-$num_file
    endif

    @ j++
  end

endif

echo 'Figures have been plotted'
