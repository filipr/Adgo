program SetOrder
  use  paramets

  integer:: maxmesh, ival, ivar
  real, dimension (:,:), allocatable :: x, y,  z, ord
  real, dimension (:,:), allocatable :: xlog, ylog
  integer, dimension(:), allocatable :: idof, ix
  integer, dimension(:,:), allocatable :: ipar
  logical :: logEOC, ordA
  integer :: ife, ife_b, ife_e, ife_txt, iEOC, iline, iform
  character*150 text(1:30)
  character*7  EOC
  character*15 typ
  character*50 input_file
  real :: told

  !RTNst variables
  real :: d_k, dualError, cTn

  ! setting of indexes of appropriate errors or estims using paramet.f90
  i_errL2 = 0 + L2
  i_errH1 = 0 + H1
  i_errH1_discrete = 0 + H1_discrete
  i_errL8 = 0 + L8
  i_interLq = 0 + interLq
  i_interL8 = 0 + interL8
  i_interH1 = 0 + interH1
  i_XX = 0 + XX
  i_Ejumps  = 0 + Ejumps
  i_EjumpsG = 0 + EjumpsG
  i_EjumpsPi = 0 + EjumpsPi

  i_errL8L2 = max_errS + L8L2
  i_errL2L2 = max_errS + L2L2
  i_errL2H1 = max_errS + L2H1
  i_errL2H10 = max_errS + L2H10
  i_Snorm1 = max_errS + Snorm1
  i_Snorm2 = max_errS + Snorm2
  i_Snorm3 = max_errS + Snorm3

  i_resA = max_errS +  max_errSTnorm + resA
  i_resS = max_errS +  max_errSTnorm + resS
  i_resT = max_errS +  max_errSTnorm + resT
  i_resST = max_errS +  max_errSTnorm + resST

  ! HO_rec
  !i_HO_estim_L2_p0 = max_errS +  max_errSTnorm + HO_estim_L2_p0
  !i_HO_estim_H1_p0 = max_errS +  max_errSTnorm + HO_estim_H1_p0
  !i_HO_trunc_L2_p0 = max_errS +  max_errSTnorm + HO_trunc_L2_p0
  !i_HO_trunc_H1_p0 = max_errS +  max_errSTnorm + HO_trunc_H1_p0

  i_HO_estim_L2_p1 = max_errS +  max_errSTnorm + HO_estim_L2_p1
  i_HO_estim_H1_p1 = max_errS +  max_errSTnorm + HO_estim_H1_p1
  i_HO_trunc_L2_p1 = max_errS +  max_errSTnorm + HO_trunc_L2_p1
  i_HO_trunc_H1_p1 = max_errS +  max_errSTnorm + HO_trunc_H1_p1

  i_HO_estim_L2_p2 = max_errS +  max_errSTnorm + HO_estim_L2_p2
  i_HO_estim_H1_p2 = max_errS +  max_errSTnorm + HO_estim_H1_p2
  i_HO_trunc_L2_p2 = max_errS +  max_errSTnorm + HO_trunc_L2_p2
  i_HO_trunc_H1_p2 = max_errS +  max_errSTnorm + HO_trunc_H1_p2
  i_HO_recovery    = max_errS +  max_errSTnorm + HO_recovery
  i_HO_rec_estim    = max_errS +  max_errSTnorm + HO_rec_estim

  ! p_robust
  i_P_flux = max_errS +  max_errSTnorm + P_flux
  i_P_rez = max_errS +  max_errSTnorm + P_rez
  i_P_FR = max_errS +  max_errSTnorm + P_FR
  i_P_pot = max_errS +  max_errSTnorm + P_pot
  i_P_BC = max_errS +  max_errSTnorm + P_BC
  i_P_tot = max_errS +  max_errSTnorm + P_tot
  i_P_sF = max_errS +  max_errSTnorm + P_sF
  i_P_su = max_errS +  max_errSTnorm + P_su
  i_P_FDu = max_errS +  max_errSTnorm + P_FDu
  i_P_potP = max_errS +  max_errSTnorm + P_potP
  i_P_F_p1 = max_errS +  max_errSTnorm + P_F_p1
  i_P_F_p2 = max_errS +  max_errSTnorm + P_F_p2
  i_P_F_p3 = max_errS +  max_errSTnorm + P_F_p3
  i_P_F_p4 = max_errS +  max_errSTnorm + P_F_p4
  i_P_s_p1 = max_errS +  max_errSTnorm + P_s_p1
  i_P_s_p2 = max_errS +  max_errSTnorm + P_s_p2
  i_P_s_p3 = max_errS +  max_errSTnorm + P_s_p3
  i_P_s_p4 = max_errS +  max_errSTnorm + P_s_p4

  ! RTNst
  i_RTN_all = max_errS +  max_errSTnorm + RTNall
  i_RTN_eta = max_errS +  max_errSTnorm + RTNeta
  i_RTN_rez = max_errS +  max_errSTnorm + RTNrez
  i_RTN_flux = max_errS +  max_errSTnorm + RTNflux
  i_RTN_radau = max_errS +  max_errSTnorm + RTNradau
  i_RTN_radau2 = max_errS +  max_errSTnorm + RTNradau2
  i_RTN_jump = max_errS +  max_errSTnorm + RTNjump
  i_RTNfluxnorm = max_errS +  max_errSTnorm + RTNfluxnorm


   ! DUA
  i_DUA_total = max_errS +  max_errSTnorm + total
  i_DUA_rez = max_errS +  max_errSTnorm + Rn
  i_DUA_flux = max_errS +  max_errSTnorm + DFn
  i_DUA_NC = max_errS +  max_errSTnorm + NC1n

  !DWR
  i_DWR_A = max_errS +  max_errSTnorm  + dwrA
  i_DWR_S = max_errS +  max_errSTnorm  + dwrS
  i_DWR_S_abs = max_errS +  max_errSTnorm  + dwrS_abs
  i_DWR_exa =   max_errS +  max_errSTnorm  + dwrE
  i_DWR_dualS = max_errS +  max_errSTnorm  + dwr_dualS
  i_DWR_dualS_abs = max_errS +  max_errSTnorm  + dwr_dualS_abs
  i_DWR_dualA = max_errS +  max_errSTnorm  + dwr_dualA
  i_DWR_Juh = max_errS +  max_errSTnorm  + dwr_Juh
  i_DWR_aver = max_errS +  max_errSTnorm  + dwr_aver
  i_DWR_aver_abs = max_errS +  max_errSTnorm  + dwr_aver_abs

  ! CFD
  i_CFD_cD = CFD_cD
  i_CFD_cL = CFD_cL
  i_CFD_cM = CFD_cM

  i_CPU  =  1  ! array z
  i_area =  2
  i_estim  =  6  ! array z

  ! ival = number of real values in one row in order.dat at the input
  ! ncounts = number of real columns in tables at the output
  logEOC  = .false.

  iEOC = 3

  !iline = 30  ! after "iline" lines the  \hline is plotted
  !iline = 21  ! after "iline" lines the  \hline is plotted
  iline = 16  ! after "iline" lines the  \hline is plotted
  !iline = 5  ! after "iline" lines the  \hline is plotted

  if (command_argument_count() == 3) then
    call get_command_argument(1,input_file)

    call get_command_argument(2,typ)
    open(12,file='smaz', status='unknown')
    write(12,*) typ
    close(12)
    open(12,file='smaz', status='OLD')
    read(12,*) indexi
    close(12)

    call get_command_argument(3,typ)

    !print*,'%Insert dependence variable: '
    !print*,'%   1 (h_max), 2 (h_aver), 3 (tau),  4 (DOF), 5(ST DOF), 6(p)'
    !!read(*,*) indexi
    !write(*,*) ' Dependence variable = ', indexi
  else
    print*,'% Syntax:  Setorderx  <order.dat> <indexi> <type_of_table>'
    print*,'% indexi:  1 (h_max), 2 (h_aver), 3 (tau),  4 (DOF), 5(ST DOF), 6(p)'
    print*,'% type of table: IPG ERR RES STDGM  AMAsteady AMAtdp AMAst  AMAstO RES_STDGM', &
      ' RES_STDGMp pNeu pNeuP pNeuAM pNeuHP HO_rec RTNst DWR CFD_RES'
    stop
 endif

 !typ = 'STDGM'     ! new ouput for ST DGM
 !typ = 'AMAst'
 !typ = 'AMAstO'
 !typ = RES_STDGM
 !!typ = 'CFD-AMA'



  maxmesh = 500
  ivar = 4
  Nx = 10
  Ny0 =  max_errS + max_errSTnorm + max_eta  ! = 50?
  Ny = Ny0 + 3 !FR changed from 2
  Nz = 11 !12
  i_tol = 9 ! state%space%adapt%tol_max
  i_CW  = 10
  i_LAiter = 11
  i_IPG = 12

  !allocate( x(1:maxmesh,ival+ivar), y(maxmesh,ival+ivar),r(ival,3), ix(maxmesh) )
  allocate( x(1:maxmesh,1:Nx), y(1:maxmesh,1:Ny), ord(0:maxmesh, 1:Ny), z(1:maxmesh, 1:Nz) )
  allocate( xlog(1:maxmesh, 1:2), ylog(1:maxmesh,1:Ny)  )
  allocate( idof(maxmesh), ipar(maxmesh, 1:5), ix(1:maxmesh) )

  EOC = 'EOC'

  imesh0 =1

  ! reading of the data
  ires = 11
  open(ires, file=input_file,status='OLD')

  do k = 1,maxmesh
     read(ires, *, end=110) ix(k), x(k,2:3), idof(k), ipar(k, 1:5), x(k,10), y(k, 1:Ny0), z(k, 1:Nz)

     !write(*,*) Nx, Ny, Nz, x(k, 1:4),y(k, 1:4)
  enddo
110 imesh = k-1

  print*,'Number of records (lines) =',imesh
  close(ires)

  ! special values
  i_DG = Ny0 + 1
  y(1:imesh,i_DG)  = (y(1:imesh,i_errH1_discrete)**2 + y(1:imesh,i_EjumpsPi)**2)**0.5
  !write(*,'(a8, 30es12.4)') 'H1', y(1:imesh,i_errH1_discrete)
  !write(*,'(a8, 30es12.4)') 'J', y(1:imesh,i_EjumpsPi)
  !write(*,'(a8, i5, 30es12.4)') 'DG', i_DG, y(1:imesh,i_DG)

  i_etaDG = Ny0 + 2
  y(1:imesh,i_etaDG) = (y(1:imesh,i_P_tot)**2 + y(1:imesh,i_EjumpsPi)**2)**0.5
  !write(*,'(a8, 30es12.4)') 'eta', y(1:imesh,i_P_tot)
  !write(*,'(a8, 30es12.4)') 'J', y(1:imesh,i_EjumpsPi)
  !write(*,'(a8, i5, 30es12.4)') 'eDG', i_etaDG, y(1:imesh,i_etaDG)


  !special weighted L2 - H1 norm for RTNst estims
  !works only if h_k = h and \tau_m = \tau !!!!!!!
  i_RTNstNorm = Ny0 + 3
  i_RTNstNormOld = Ny0 + 4

  do i = imesh0, imesh
   d_k = sqrt( x(i,2)**2. + x(i,3)**2. )
   !d_K = 1.0
   !cTn = ( 1.0 / x(i,3)**2. ) + ( 1.0 / x(i,2)**2 )

   !Mila
   ! FR commented because caused problems in DWR
!   y(i, i_RTNstNormOld) = sqrt( (d_k**2.0 / x(i,3)**2.0 ) * y(i,i_errL2L2)**2. + &
!      (d_k**2.0 / x(i,2)**2.0 ) * y(i,i_errL2H1)**2.0 )

      !!! new version - temporal L2 + \sigma(u) norm

    y(i, i_RTNstNorm) = sqrt( (d_k**2.0 / x(i,3)**2.0 ) * y(i,i_errL2L2)**2. + &
      (d_k**2.0 / x(i,2)**2.0 ) * y(i,i_RTNfluxnorm)**2.0 )


! Vohralik
!      y(i, i_RTNstNorm) = sqrt( ( cTn**(-1.0) / x(i,3)**2.0 ) * y(i,i_errL2L2)**2. + &
!      ( cTn**(-1.0) / x(i,2)**2.0 ) * y(i,i_errL2H1)**2.0 )

  enddo

 if( typ == 'HO_rec' .or.typ == 'HO_recC' .or. typ == 'HO_recA'  ) then
     y(1:imesh,i_DG)  = (y(1:imesh,i_errH1)**2 + y(1:imesh,i_EjumpsG)**2)**0.5
     y(1:imesh,i_etaDG) = (y(1:imesh,i_HO_estim_H1_p2)**2 + y(1:imesh,i_EjumpsG)**2)**0.5
  endif


  ! end of special values

  if(indexi == 1) then ! h_max
     xlog(1:imesh, 1) =  x(1:imesh,2)

  elseif(indexi == 2) then ! h_aver
     xlog(1:imesh, 1) = z(1:imesh,i_area) /((ix(1:imesh)/2)**0.5)

  elseif(indexi == 3) then ! tau
     xlog(1:imesh, 1) =  x(1:imesh,3)

  elseif(indexi == 4) then ! 1/sqrt(DOF)
     !xlog(1:imesh,1) = 1./((idof(1:imesh))**0.5)
     xlog(1:imesh,1) = 1./((idof(1:imesh)) )  !**(1./3))

  elseif(indexi == 5) then ! 1/(ST DOF)
     xlog(1:imesh,4) = (x(1:imesh, 3) / idof(1:imesh) ) **(1./3)
     !write(*,'(12es12.4)') x(1:imesh, 3)
     !write(*,'(12es12.4)') 1.*idof(1:imesh)
     !write(*,'(12es12.4)') (x(1:imesh, 3) / idof(1:imesh))**(-1)
     !write(*,'(12es12.4)') x(1:imesh, 4)
     indexi = 4

  elseif(indexi == 6) then ! p - degree of polynomial approximation
     xlog(1:imesh,1) = ipar(1:imesh,1)
  endif

  ! computing of the logarithm
  xlog(:,2) = 0.
  ylog(:,:) = 0.
  do k=1,imesh
     if( xlog(k, 1) > 0.)  xlog(k, 2)  = log(xlog( k, 1) )
     do i=1,Ny
        if( y(k, i) > 1E-20) ylog(k,i ) = log( y(k, i ) )
     enddo
  enddo

  !do k=1,imesh
  !   print*,'!!!',k,ylog(k,:)
  !enddo

  !Ny =3

  ! computation of the regresion
  do l=1,Ny
     ! global
     call Regres(imesh-imesh0 +1, xlog(imesh0:imesh,2), ylog(imesh0:imesh,l), ord(0, l) )

    ! locals
     do k=imesh0,imesh-1
        call Regres(2, xlog(k:k+1,2), ylog(k:k+1, l), ord(k+1, l) )
     enddo
  end do

  !do k=1, 2
  !   write(*,'(a6,60es12.4)')'dat',x(k,2),xlog(k, 1:2),ylog(k,1:Ny)
  !enddo

  !do k=2, 2
  !   write(*,'(a6,60es12.4)')'ord', ord(k, 1:Ny)
  !enddo

  !x(imesh0:imesh, 2) = ix(imesh0:imesh)

  ! SQUARE ROOT of number of elements is N
  !ix(imesh0:imesh) = (ix(imesh0:imesh)/2)**0.5
  !idof(imesh0:imesh) = (idof(imesh0:imesh))**0.5

  ife_b = 25
  open (ife_b,file='tab_beg.tex', status='UNKNOWN')

  ife_e = 26
  open (ife_e,file='tab_end.tex', status='UNKNOWN')

  ife = 22
  open (ife,file='tab.tex', status='UNKNOWN')

  ife_txt = 23
  open (ife_txt,file='tab.txt', status='UNKNOWN')




  if(typ == 'IPG') then

     ncounts = 3 ! 8   ! number of columns in tab.tex
     logEOC  = .false.
     iEOC = 2
     iline = 1000000

     !text(1) = "$|E_I(u_h)|_{L^2}$"
     text(1) = "$\|e_h \|_{L^2}$" !(\Omega)}$"
     text(2) = "$|e_h|_{H^1}$" !(\Omega)}$"
     text(3) = " iter " !(\Omega)}$"

     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"

     write(ife_b, *) "\begin{tabular}{crrr|cc|r}"

     write(ife_b, *) "\hline"
     write(ife_b, *) " method & $p$ &  $h$  & $C_W$ "  !!! &

     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh
        ! errors

        if(  int(0.5+z(i, i_IPG)) == 1) write(ife, *) 'NIPG'
        if(  abs( z(i, i_IPG)) < 0.1) write(ife, *) 'IIPG'
        if(  (0.5+z(i, i_IPG)) < 0) write(ife, *) 'SIPG'

        write(ife, *) ' &  ',ipar(k,1)
        write(ife,'(a8,i6,a6)') '& $ 1/ ',int( (x(i,2)/2**0.5 )**(-1)+0.5),'$ '
        write(ife,'(a2, f6.0)') '&', z(i, i_CW)

        !if(i > 1) then

         !!     write(ife, 952) y(i,i_errL2), ord(i,i_errL2)
         !     write(ife, 952) y(i,i_errH1), ord(i,i_errH1)
        write(ife, 201) y(i,i_errL2)
        write(ife, 201) y(i,i_errH1)
        write(ife,'(a2, i8)') '&', int(0.5+z(i, i_LAiter))

        !write(ife, 301) z(i,i_CPU)  ! CPU time
        write(ife, 300)
     end do
     write(ife, *) "\hline"


  elseif(typ == 'ERR' ) then
     ncounts = 4   ! number of columns in tab.tex
     !logEOC  = .false.
     logEOC  = .true.
     iEOC = 3
     iline = 1000000
     !if(typ == 'AMAstO') logEOC  = .true.

     !text(1) = "$|E_I(u_h)|_{L^2}$"
     text(1) = "$\| e_h\|_{L^2}$" !(\Omega)}$"
     text(2) = "$\| e_h\|_{H^1}$" !(\Omega)}$"
     !text(3) = "$\| e_h\|_{X}$" !(\Omega)}$"
     text(3) = "$ J_h(e_h)$" !(\Omega)}$"
     !text(5) = "$\eta_A$" !(\Omega)}$"
     !text(6) = "$\eta_S$"  !!!$|E_I(u_h)|_\infty$"
     !text(7) = "$i_{\mathrm{eff}}$"  !!!$|E_I(u_h)|_\infty$"
     text(4) = "CPU(s)"

     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"

     write(ife_b, *) "\begin{tabular}{crr||cc|cc|cc|r}"

     write(ife_b, *) "\hline"
     write(ife_b, *) " $p$ &  $h$  & DoF "  !!! &

     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh
        write(ife, *) ipar(k,1)
        write(ife,'(a8,i6,a6)') '& $ 1/ ',int( (x(i,2)/2**0.5 )**(-1)+0.5),'$ '
        write(ife, '(a3, i 8)')  '&', idof(i)

        if(i > 1) then
           write(ife, 952) y(i,i_errL2), ord(i,i_errL2)
           write(ife, 952) y(i,i_errH1), ord(i,i_errH1)
           write(ife, 952) y(i,i_Ejumps), ord(i,i_Ejumps)
           !write(ife, 952) y(i,i_errL8L2), ord(i,i_errL8L2)
           !write(ife, 952) y(i,i_errL2L2), ord(i,i_errL2L2)
           !write(ife, 952) y(i,i_errL2H1), ord(i,i_errL2H1)
           !write(ife, 952) y(i,i_interLq), ord(i,i_interLq)
           !write(ife, 952) y(i,i_interL8), ord(i,i_interL8)
        else

           write(ife, 992) y(i,i_errL2)
           write(ife, 992) y(i,i_errH1)
           !write(ife, 201) y(i,i_XX)
           write(ife, 992) y(i,i_Ejumps)
           !write(ife, 201) y(i,i_resA)
           !write(ife, 201) y(i,i_resS)
           !write(ife, 701) y(i,i_resS) / max(1E-15, y(i,i_XX))    ! i_eff
        endif

        write(ife, 301) z(i,i_CPU)  ! CPU time
           !write(ife, 300)

        !if(i > 1e+6) then
        !if(i > 1) then
        !   write(ife,*) '  \multicolumn{3}{c||}{{\scriptsize (EOC)}} '
        !
        !   write(ife, 621) ord(i,i_errL2)
        !   write(ife, 621) ord(i,i_errH1)
        !   write(ife, 621) ord(i,i_XX)
        !   write(ife, 621) ord(i,i_Ejumps)
        !   write(ife, 621) ord(i,i_resA)
        !   write(ife, 621) ord(i,i_resS)
        !   write(ife, *) ' & '
        !   write(ife, 300)
        !
        !endif

        !write(ife, 600) x(i, 31) / x(i,27) ! i_{eff}
        !write(ife, 600) x(i, 16) / x(i,7) ! i_{eff}
        !write(ife, 600) (x(i, 16)**2 + x(i,8)**2)**0.5 / (x(i,7)**2 + x(i,8)**2)**0.5 ! i_{eff}

     end do
     write(ife, *) "\hline"


  elseif(typ == 'RES' ) then
     ncounts = 8   ! number of columns in tab.tex
     logEOC  = .false.
     !logEOC  = .true.
     iEOC = 7
     iline = 1000000
     !if(typ == 'AMAstO') logEOC  = .true.

     !text(1) = "$|E_I(u_h)|_{L^2}$"
     text(1) = "$\| e_h\|_{L^2}$" !(\Omega)}$"
     text(2) = "$\| e_h\|_{H^1}$" !(\Omega)}$"
     text(3) = "$\| e_h\|_{X}$" !(\Omega)}$"
     text(4) = "$ J_h(e_h)$" !(\Omega)}$"
     text(5) = "$\eta_A$" !(\Omega)}$"
     text(6) = "$\eta_S$"  !!!$|E_I(u_h)|_\infty$"
     text(7) = "$i_{\mathrm{eff}}$"  !!!$|E_I(u_h)|_\infty$"
     text(8) = "CPU(s)"

     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"

     write(ife_b, *) "\begin{tabular}{crr||cccc|cc|c|r}"

     write(ife_b, *) "\hline"
     write(ife_b, *) " $\ell$ &  $N_h$  & $\Nhp $ "  !!! &

     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh
         write(ife, 803) i-1,  ix(i), idof(i)

        write(ife, 201) y(i,i_errL2)
        write(ife, 201) y(i,i_errH1)
        write(ife, 201) y(i,i_XX)
        write(ife, 201) y(i,i_Ejumps)
        write(ife, 201) y(i,i_resA)
        write(ife, 201) y(i,i_resS)
        write(ife, 701) y(i,i_resS) / max(1E-15, y(i,i_XX))    ! i_eff
        write(ife, 301) z(i,i_CPU)  ! CPU time
        !write(ife, 300)

        !if(i > 1e+6) then
        if(i > 1) then
           write(ife,*) '  \multicolumn{3}{c||}{{\scriptsize (EOC)}} '

           write(ife, 621) ord(i,i_errL2)
           write(ife, 621) ord(i,i_errH1)
           write(ife, 621) ord(i,i_XX)
           write(ife, 621) ord(i,i_Ejumps)
           write(ife, 621) ord(i,i_resA)
           write(ife, 621) ord(i,i_resS)
           write(ife, *) ' & '
           write(ife, 300)

        endif

        !write(ife, 600) x(i, 31) / x(i,27) ! i_{eff}
        !write(ife, 600) x(i, 16) / x(i,7) ! i_{eff}
        !write(ife, 600) (x(i, 16)**2 + x(i,8)**2)**0.5 / (x(i,7)**2 + x(i,8)**2)**0.5 ! i_{eff}

     end do
     write(ife, *) "\hline"


  elseif(typ == 'CFD_RES' ) then
     ncounts = 7   ! number of columns in tab.tex
     logEOC  = .false.
     !logEOC  = .true.
     iEOC = 7
     iline = 1000000
     !if(typ == 'AMAstO') logEOC  = .true.

     !text(1) = "$|E_I(u_h)|_{L^2}$"
     text(1) = "$ c_D$" !(\Omega)}$"
     text(2) = "$ c_L$" !(\Omega)}$"
     text(3) = "$ c_M$" !(\Omega)}$"
     text(4) = "$ J_h(e_h)$" !(\Omega)}$"
     text(5) = "$\eta_A$" !(\Omega)}$"
     text(6) = "$\eta_S$"  !!!$|E_I(u_h)|_\infty$"
     !text(7) = "$i_{\mathrm{eff}}$"  !!!$|E_I(u_h)|_\infty$"
     text(7) = "CPU(s)"

     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"

     write(ife_b, *) "\begin{tabular}{crr||cccc|cc|r}"

     write(ife_b, *) "\hline"
     write(ife_b, *) " $\ell$ &  $N_h$  & $\Nhp $ "  !!! &

     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh
         write(ife, 803) i-1,  ix(i), idof(i)

        write(ife, 291) y(i,i_CFD_cD)
        write(ife, 291) y(i,i_CFD_cL)
        write(ife, 291) y(i,i_CFD_cM)
        write(ife, 201) y(i,i_Ejumps)
        write(ife, 201) y(i,i_resA)
        write(ife, 201) y(i,i_resS)
        !write(ife, 701) y(i,i_resS) ! / max(1E-6, y(i,i_XX))    ! i_eff
        write(ife, 301) z(i,i_CPU)  ! CPU time
        !write(ife, 300)

        if(i > 1e+6) then
        !if(i > 1) then
           write(ife,*) '  \multicolumn{3}{c||}{{\scriptsize (EOC)}} '

           write(ife, 621) ord(i,i_errL2)
           write(ife, 621) ord(i,i_errH1)
           write(ife, 621) ord(i,i_XX)
           write(ife, 621) ord(i,i_Ejumps)
           write(ife, 621) ord(i,i_resA)
           write(ife, 621) ord(i,i_resS)
           write(ife, *) ' & '
           write(ife, 300)

        endif

        !write(ife, 600) x(i, 31) / x(i,27) ! i_{eff}
        !write(ife, 600) x(i, 16) / x(i,7) ! i_{eff}
        !write(ife, 600) (x(i, 16)**2 + x(i,8)**2)**0.5 / (x(i,7)**2 + x(i,8)**2)**0.5 ! i_{eff}

     end do
     write(ife, *) "\hline"


  elseif(typ == 'DWR' ) then
     ncounts = 10   ! number of columns in tab.tex
     logEOC  = .false.
     iEOC = 7
     iline = 1000000
     !if(typ == 'AMAstO') logEOC  = .true.

     !text(1) = "$|E_I(u_h)|_{L^2}$"
!     text(1) = "$\| e_h\|_{L^2}$" !(\Omega)}$"
!     text(2) = "$\| e_h\|_{H^1}$" !(\Omega)}$"
     !text(3) = "$\| e_h\|_{X}$" !(\Omega)}$"
     text(1) = "$ J_h(u_h)$" !(\Omega)}$"
     text(2) = "$ J_h(e_h)$" !(\Omega)}$"
     text(3) = "$\eta_A$" !(\Omega)}$"
     text(4) = "$\eta_S$"  !!!$|E_I(u_h)|_\infty$"
     text(5) = "$\eta_{S_{dual}}$"  !!!$|E_I(u_h)|_\infty$"
     text(6) = "$\eta_{S_{aver}}$"  !!!$|E_I(u_h)|_\infty$"
     text(7) = "$i_{\mathrm{primal}}$"  !!!$|E_I(u_h)|_\infty$"
     text(8) = "$i_{\mathrm{dual}}$"  !!!$|E_I(u_h)|_\infty$"
     text(9) = "$i_{\mathrm{aver}}$"  !!!$|E_I(u_h)|_\infty$"
     text(10) = "CPU(s)"

     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"

!     write(ife_b, *) "\begin{tabular}{crr||cccc|ccc|cc|r}"
     write(ife_b, *) "\begin{tabular}{crr||cc|cccc|ccc|r}"

     write(ife_b, *) "\hline"
     write(ife_b, *) " $\ell$ &  $N_h$  & $\Nhp $ "  !!! &

     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh
         write(ife, 803) i-1,  ix(i), idof(i)

!        write(ife, 201) y(i,i_errL2)
!        write(ife, 201) y(i,i_errH1)

!        write(ife, 201) y(i,i_XX)
        write(ife, 201) y(i,i_DWR_Juh)
        write(ife, 201) y(i,i_DWR_exa)
        write(ife, 201) y(i,i_DWR_A)
        write(ife, 201) y(i,i_DWR_S)
        write(ife, 201) y(i,i_DWR_dualS)
        write(ife, 201) y(i,i_DWR_aver)
!        write(ife, 201) y(i,i_DWR_S_abs)
        write(ife, 701) y(i,i_DWR_S) / y(i,i_DWR_exa)    ! i_eff
        write(ife, 701) y(i,i_DWR_dualS) / y(i,i_DWR_exa)    ! i_eff
!        write(ife, 701) y(i,i_DWR_S_abs) / y(i,i_DWR_exa)    ! i_eff
        write(ife, 701) y(i,i_DWR_aver) / y(i,i_DWR_exa)    ! i_eff

        write(ife, 301) z(i,i_CPU)  ! CPU time

!        write(ife, 201) y(i,i_errL2)
!        write(ife, 201) y(i,i_errH1)
!!        write(ife, 201) y(i,i_XX)
!        write(ife, 201) y(i,i_DWR_exa)
!        write(ife, 201) y(i,i_DWR_A)
!        write(ife, 201) y(i,i_DWR_S)
!        write(ife, 201) y(i,i_DWR_S_abs)
!        write(ife, 701) y(i,i_DWR_S) / y(i,i_DWR_exa)    ! i_eff
!        write(ife, 301) z(i,i_CPU)  ! CPU time
!        !write(ife, 300)

       ! if(i > 1e+6) then
           write(ife,*) '  \multicolumn{3}{c||}{{\scriptsize (EOC)}} '

           !write(ife, 621) ord(i,i_errL2)
           !write(ife, 621) ord(i,i_errH1)
           write(ife, *) ' & ' !ord(i,i_DWR_Juh)
           write(ife, 621) ord(i,i_DWR_exa)
           write(ife, 621) ord(i,i_DWR_A)
           write(ife, 621) ord(i,i_DWR_S)
           write(ife, 621) ord(i,i_DWR_dualS)
           write(ife, 621) ord(i,i_DWR_aver)
!           write(ife, 621) ord(i,i_DWR_S_abs)
           write(ife, *) ' & '
           write(ife, *) ' & '
           write(ife, *) ' & '
           write(ife, 300)

        !endif

        !write(ife, 600) x(i, 31) / x(i,27) ! i_{eff}
        !write(ife, 600) x(i, 16) / x(i,7) ! i_{eff}
        !write(ife, 600) (x(i, 16)**2 + x(i,8)**2)**0.5 / (x(i,7)**2 + x(i,8)**2)**0.5 ! i_{eff}

     end do
     write(ife, *) "\hline"

  else if(typ == 'AMAst' .or. typ == 'AMAstO'  ) then

     ncounts = 8   ! number of columns in tab.tex
     logEOC  = .false.
     iEOC = 7
     iline = 1000000
     if(typ == 'AMAstO') logEOC  = .true.

     !text(1) = "$|E_I(u_h)|_{L^2}$"
     text(1) = "$\|\ehp\|_{L^2}$" !(\Omega)}$"
     text(2) = "$\|\ehp\|_{H^1}$" !(\Omega)}$"
     text(3) = "$\|\ehp\|_{L^{\infty}(L^2)}$" !(\Omega)}$"
     text(4) = "$\|\ehp\|_{L^2(L^2)}$" !(\Omega)}$"
     text(5) = "$\|\ehp\|_{L^2(H^1)}$" !(\Omega)}$"
     text(6) = "$\|E_I\|_{L^2}$"  !!!$|E_I(u_h)|_\infty$"
     text(7) = "$\|E_I\|_{L^\infty}$"  !!!$|E_I(u_h)|_\infty$"
     text(8) = "CPU(s)"

     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"

     if(typ == 'AMAstO')  then
        write(ife_b, *) "\begin{tabular}{crrr||cc|cc||cc|cc|cc||cc|cc|r}"
     else
        write(ife_b, *) "\begin{tabular}{crrr||cc|ccc|cc|r}"
     endif

     write(ife_b, *) "\hline"
     write(ife_b, *) " $\ell$ &  $N_h$  & $\Nhp $ & $t$"  !!! &

     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh
        if( abs(told - x(i,10))> 1E-5 )  write(ife,*) "\hline"
        write(ife, 803) i-1,  ix(i), idof(i)
        write(ife, '(a1, f10.2)') '&', x(i,10)
        told = x(i,10)


        ! errors
        if( typ == 'AMAstO'  ) then
           if(i > 1) then

              write(ife, 952) y(i,i_errL2), ord(i,i_errL2)
              write(ife, 952) y(i,i_errH1), ord(i,i_errH1)
              write(ife, 952) y(i,i_errL8L2), ord(i,i_errL8L2)
              write(ife, 952) y(i,i_errL2L2), ord(i,i_errL2L2)
              write(ife, 952) y(i,i_errL2H1), ord(i,i_errL2H1)
              write(ife, 952) y(i,i_interLq), ord(i,i_interLq)
              write(ife, 952) y(i,i_interL8), ord(i,i_interL8)
           else
              write(ife, 992) y(i,i_errL2)
              write(ife, 992) y(i,i_errH1)
              write(ife, 992) y(i,i_errL8L2)
              write(ife, 992) y(i,i_errL2L2)
              write(ife, 992) y(i,i_errL2H1)
              write(ife, 992) y(i,i_interLq)
              write(ife, 992) y(i,i_interL8)
           endif
        else
           write(ife, 201) y(i,i_errL2)
           write(ife, 201) y(i,i_errH1)
           write(ife, 201) y(i,i_errL8L2)
           write(ife, 201) y(i,i_errL2L2)
           write(ife, 201) y(i,i_errL2H1)
           write(ife, 201) y(i,i_interLq)
           write(ife, 201) y(i,i_interL8)
        endif

        !write(ife, 600) x(i, 31) / x(i,27) ! i_{eff}
        !write(ife, 600) x(i, 16) / x(i,7) ! i_{eff}
        !write(ife, 600) (x(i, 16)**2 + x(i,8)**2)**0.5 / (x(i,7)**2 + x(i,8)**2)**0.5 ! i_{eff}

        write(ife, 301) z(i,i_CPU)  ! CPU time
        !write(ife, 300)
     end do
     write(ife, *) "\hline"


  else if(typ == 'STDGM'  ) then
     ncounts = 7   ! number of columns in tab.tex
     logEOC  = .true.
     iEOC = 7
     iline = 1000000

     text(1) = "$\|e_h\|_{L^\infty(L^2)}$"
     text(2) = "$\|e_h\|_{L^2(L^2)}$"
     text(3) = "$\|e_h\|_{L^2(H^1)}$"
     text(4) = "$\|e_h\|_{L^2(H^1_0)}$"
     text(5) = "$\|e_h(t_m)\|_{L^2(\Omega)}$"
     text(6) = "$\|e_h(t_r)\|_{L^2(\Omega)}$"
     text(7) = "$\|e_h(t_i)\|_{L^2(\Omega)}$"

     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"

     write(ife_b, *) "\begin{tabular}{cccc|cc|cc|cc|cc|cc|cc|cc}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " $h$ & $\tau $ &  $P_k$ & $T_m$ "

     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh

        write(ife, 99) x(i,2), x(i,3)  ! h_aver, tau
        write(ife, 100) ipar(i, 1:2)   ! P_k, T_m, ..

        if(i > 1) then
           write(ife, 952) y(i,i_errL8L2), ord(i,i_errL8L2)
           write(ife, 952) y(i,i_errL2L2), ord(i,i_errL2L2)
           write(ife, 952) y(i,i_errL2H1), ord(i,i_errL2H1)
           write(ife, 952) y(i,i_errL2H10), ord(i,i_errL2H10)
           write(ife, 952) y(i,i_Snorm1) , ord(i,i_Snorm1)
           write(ife, 952) y(i,i_Snorm2) , ord(i,i_Snorm2)
           write(ife, 952) y(i,i_Snorm3) , ord(i,i_Snorm3)

        else
           write(ife, 992) y(i,i_errL8L2)
           write(ife, 992) y(i,i_errL2L2)
           write(ife, 992) y(i,i_errL2H1)
           write(ife, 992) y(i,i_errL2H10)
           write(ife, 992) y(i,i_Snorm1)
           write(ife, 992) y(i,i_Snorm2)
           write(ife, 992) y(i,i_Snorm3)
        endif

        !write(ife, 600) x(i, 31) / x(i,27) ! i_{eff}
        !write(ife, 600) x(i, 16) / x(i,7) ! i_{eff}
        !write(ife, 600) (x(i, 16)**2 + x(i,8)**2)**0.5 / (x(i,7)**2 + x(i,8)**2)**0.5 ! i_{eff}

        !write(ife, 301) z(i,i_CPU)  ! CPU time
        write(ife, 300)
     end do
     write(ife, *) "\hline"

else if(typ == 'STDGM2'  ) then
     ncounts = 4  ! number of columns in tab.tex
     logEOC  = .true.
     iEOC = 7
     iline = 1000000

  else if(typ == 'RES_STDGM' ) then
     ncounts = 11   ! number of columns in tab.tex
     !logEOC  = .true.
     logEOC  = .false.
     iEOC = 5
     iline = 4

     text(1) = "$\|e_h\|_{L^\infty(L^2)}$"
     text(2) = "$\|e_h\|_{L^2(L^2)}$"
     text(3) = "$\|e_{h\tau}\|_{L^2(H^1)}$"
     text(4) = "$ \eta_A$"
     text(5) = "$\eta_S $"
     text(6) = "$\eta_T $"
     text(7) = "$\eta_{ST}$"
     ! adaptST
     text(8) = "$\eta_T / \eta_S$"
     text(9) = "$\frac{\eta}{e_{L^2H^1}}$"
     text(10) = "$\frac{\eta}{e_{L^2L^2}}$"
     text(11) = "CPU(s)"

     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"

     if(logEOC) then
        write(ife_b, *) "\begin{tabular}{cccc|c|cc|cc|cc|cc|cc|c|c|r}"
        write(ife_b, *) "\hline"
        write(ife_b, *) " $h$ & $\tau $ &  $p$ & $q$ "
     else
        write(ife_b, *) "\begin{tabular}{cccc|ccc|cccc|c|cc|r}"
        write(ife_b, *) "\hline"
        ! adaptST
        write(ife_b, *) "   $p$ & $q$ & $N$ & $\tau $ "
        !write(ife_b, *) "  $p$ & $q$ & $\omega $ "
     endif



     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh

         if(logEOC) then
            write(ife, 89) x(i,2), x(i,3)  ! h_aver, tau
            write(ife, 100) ipar(i, 1:2)   ! P_k, T_m, ..
         else
            write(ife,*) '%', x(i,2), ix(i)  ! h, N
            write(ife, 109) ipar(i, 1:2)   ! P_k, T_m, ..
            write(ife, 95)  ix(i)  !  N
            write(ife, 96) x(i,3)  ! tau
            ! adaptST
            !write(ife, *) ' & adapt'
            !write(ife, '(a3,es8.1)' ) ' & ', 0.4 / 4**(-1)
         endif

        !write(ife, 198) ipar(i, 1), idof(i)   ! P_k, DOF, ..

        if(.not. logEOC) then
           write(ife,287)  y(i,i_errL8L2), y(i,i_errL2L2), &
                y(i,i_errL2H1), &
                y(i,i_resA), y(i,i_resS), y(i,i_resT), y(i,i_resST)


        else
           if(i > 1) then
              write(ife, 952) y(i,i_errL8L2), ord(i,i_errL8L2)
              write(ife, 852) y(i,i_errL2L2), ord(i,i_errL2L2)
              write(ife, 852) y(i,i_errL2H1), ord(i,i_errL2H1)
              write(ife, 852) y(i,i_resA) , ord(i,i_resA)
              write(ife, 852) y(i,i_resS) , ord(i,i_resS)
              write(ife, 852) y(i,i_resT) , ord(i,i_resT)
              write(ife, 852) y(i,i_resST) , ord(i,i_resST)

           else
              write(ife, 992) y(i,i_errL8L2)
              write(ife, 892) y(i,i_errL2L2)
              write(ife, 892) y(i,i_errL2H1)

              write(ife, 892) y(i,i_resA)
              write(ife, 892) y(i,i_resS)
              write(ife, 892) y(i,i_resT)
              write(ife, 892) y(i,i_resST)
           endif
        endif


        ! adaptST
        if(y(i,i_resT) / y(i,i_resS) > 1E-2) then
           write(ife, 201) y(i,i_resT) / y(i,i_resS)
        else
           write(ife, 1201) y(i,i_resT) / y(i,i_resS)
        endif


        write(ife, 521) y(i, i_resST) / y(i,i_errL2H1) ! i_{eff}
        write(ife, 521) y(i, i_resST) / y(i,i_errL2L2) ! i_{eff}
        !write(ife, 600) (x(i, 16)**2 + x(i,8)**2)**0.5 / (x(i,7)**2 + x(i,8)**2)**0.5 ! i_{eff}

        write(ife, 301) z(i,i_CPU)  ! CPU time
        !write(ife, 300)

        !if(.not. logEOC .and. i > 1) then
        !   write(ife,*) '  \multicolumn{4}{c|}{{\scriptsize (EOC)}} '
        !
        !   write(ife,627)  ord(i,i_errL8L2), ord(i,i_errL2L2), &
        !        ord(i,i_errL2H1), &
        !        ord(i,i_resA), ord(i,i_resS), ord(i,i_resT), ord(i,i_resST)
        !    write(ife, *) ' & &  & '
        !   write(ife, 300)
        !endif


        if( mod(i, iline) == 0 )  write(ife,*) "\hline"

     end do
     write(ife, *) "\hline"

  else if( typ == 'RES_STDGMp'  ) then
     ncounts = 2   ! number of columns in tab.tex
     !logEOC  = .false.
     logEOC  = .true.
     iEOC = 1
     iline = 5
     ordA = .false.
     !ordA = .true.
     if(ordA) ncounts = 7
     if(ordA) iline = 100

     !text(1) = "$\|e_h\|_{L^\infty(L^2)}$"
     text(1) = "$\|e_h\|_{L^2(L^2)}$"
     !text(1) = "$\|e_{h\tau}\|_{L^2(H^1)}$"
     !text(2) = "$ \eta_A$"
     !text(3) = "$\eta_S $"
     !text(4) = "$\eta_T $"
     !text(5) = "$\eta_{ST}$"
     !text(6) = "$i_X$"
     text(2) = "CPU(s)"
     if(ordA) then
        text(6) = "$i_A$"
        text(7) = "CPU(s)"
     endif

     write(ife_b, *) "{\footnotesize           %  vetsi"
     write(ife_b, *) "%{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"

     if(.not. ordA) then
        if(logEOC) then
           !!write(ife_b, *) "\begin{tabular}{cccc|cc}"
           write(ife_b, *) "\begin{tabular}{cccc|cc|r}"
        else
           write(ife_b, *) "\begin{tabular}{cccc|c}"
        endif

        write(ife_b, *) "\hline"
        write(ife_b, *) " $h$ &   $p$ & $\tau $ & $q$ "
     else
        write(ife_b, *) "\begin{tabular}{c|ccccc|c|c}"
        write(ife_b, *) "\hline"
        write(ife_b, *) " $c_A$ "
     endif




     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh

        if(.not. ordA) then
           !write(ife, 99) x(i,2), x(i,3)  ! h_aver, tau

           !!write(ife,'(a6,i6,a6)') '$ 1/ ',int( (x(i,2)/2**0.5 )**(-1)+0.5),'$ '
           !write(ife, 98) x(i,2)   ! h_aver
           write(ife, 98) x(i,2)   ! h_max

           write(ife, 108) ipar(i, 1)   ! P_k

           !!write(ife,'(a6,i6,a6)') '&$ 1/ ',int(  x(i,3)**(-1) +0.5),'$ '
           write(ife, 97) x(i,3)   ! tau
           !write(ife,*) '& adapt'

           write(ife, 108) ipar(i, 2)   ! T_m


           if( logEOC) then
              if(i > 1) then
                  write(ife, 852) y(i,i_errL2L2), ord(i,i_errL2L2)
               else
                  write(ife, 892) y(i,i_errL2L2)
               endif
            else

               !write(ife, 198) ipar(i, 1), idof(i)   ! P_k, DOF, ..
               write(ife, 286) & ! y(i,i_errL2L2),
                    y(i,i_errL2H1), &
                    y(i,i_resA), y(i,i_resS), y(i,i_resT), y(i,i_resST),&
                    y(i,i_resST) /  y(i,i_errL2H1)
            endif

            !write(ife, 300)
            write(ife, 301) z(i,i_CPU)  ! CPU time

            !write(ife,'(a3, es12.2,a3,f9.1,a3)') ' % ', y(i,i_resA)/ y(i,i_resS), ' & ',z(i,i_CPU) , '\\'

        else
           write(ife, *) '%  ',x(i,2), x(i,3)  ! h_aver, tau
           write(ife, *)'%   ', ipar(i, 1:2)   ! P_k, T_m, ..
           !write(ife, '( es12.3)') 1./ ( 2.**(i-1))
           write(ife,'(a6,i6,a6)') ' $ 1/ ', 2**(i-1),'$ '

           write(ife, 286) & ! y(i,i_errL2L2),
                y(i,i_errL2H1), &
                y(i,i_resA), y(i,i_resS), y(i,i_resT), y(i,i_resST),&
                y(i,i_resA) /  y(i,i_resS)

           write(ife, 301) z(i,i_CPU)  ! CPU time
        endif

        if(.not. logEOC .and. (i >1 .and. .not. ordA) ) then
           write(ife,*) ' \multicolumn{4}{c|}{{\scriptsize (EOC)}} '
           write(ife, 625) &!ord(i,i_errL2L2),
                ord(i,i_errL2H1), &
                ord(i,i_resA), ord(i,i_resS), ord(i,i_resT), ord(i,i_resST)
           write(ife,*) ' &  \\'
        endif


        if( mod(i, iline) == 0 )  write(ife,*) "\hline"

     end do
     write(ife, *) "\hline"


  else if(typ == 'STDGM_book' ) then
     ncounts = 3   ! number of columns in tab.tex
     logEOC  = .true.
     !logEOC  = .false.
     iEOC = 5
     iline = 5

     text(1) = "$\|e_h\|_{L^\infty(L^2)}$"
     text(2) = "$\|e_h\|_{L^2(L^2)}$"
     text(3) = "$\|e_{h\tau}\|_{L^2(H^1)}$"

     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"

     if(logEOC) then
        write(ife_b, *) "\begin{tabular}{cc|cc|cc|cc}"
        write(ife_b, *) "\hline"
        !write(ife_b, *) " $h$ & $\tau $ &  $p$ & $q$ "
        write(ife_b, *) " $\tau $ &  $q$ "
     else
        write(ife_b, *) "\begin{tabular}{ccc|c|cccc|c|c|r}"
        write(ife_b, *) "\hline"
        ! adaptST
        write(ife_b, *) "  $p$ & $q$ & $\tau $ "
        !write(ife_b, *) "  $p$ & $q$ & $\omega $ "
     endif



     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh

         if(logEOC) then
            !write(ife, 89) x(i,2), x(i,3)  ! h_aver, tau
            !write(ife, 100) ipar(i, 1:2)   ! P_k, T_m, ..
            write(ife, 98) x(i,3)  ! tau
            write(ife, 90) ipar(i, 2)   ! P_k, T_m, ..
         else
            write(ife,*) '%', x(i,2), ix(i)  ! h, N
            write(ife, 109) ipar(i, 1:2)   ! P_k, T_m, ..
            write(ife, 96) x(i,3)  ! tau
            ! adaptST
            !write(ife, *) ' & adapt'
            !write(ife, '(a3,es8.1)' ) ' & ', 0.4 / 4**(-1)
         endif

        !write(ife, 198) ipar(i, 1), idof(i)   ! P_k, DOF, ..

        if(.not. logEOC) then
           write(ife,205) & !y(i,i_errL2L2),
                y(i,i_errL2H1), &
                y(i,i_resA), y(i,i_resS), y(i,i_resT), y(i,i_resST)
        else
           if(i > 1) then
              write(ife, 952) y(i,i_errL8L2), ord(i,i_errL8L2)
              write(ife, 852) y(i,i_errL2L2), ord(i,i_errL2L2)
              write(ife, 852) y(i,i_errL2H1), ord(i,i_errL2H1)
              !write(ife, 852) y(i,i_resA) , ord(i,i_resA)
              !write(ife, 852) y(i,i_resS) , ord(i,i_resS)
              !write(ife, 852) y(i,i_resT) , ord(i,i_resT)
              !write(ife, 852) y(i,i_resST) , ord(i,i_resST)

           else
              write(ife, 992) y(i,i_errL8L2)
              write(ife, 892) y(i,i_errL2L2)
              write(ife, 892) y(i,i_errL2H1)

              !write(ife, 892) y(i,i_resA)
              !write(ife, 892) y(i,i_resS)
              !write(ife, 892) y(i,i_resT)
              !write(ife, 892) y(i,i_resST)
           endif
        endif

        !write(ife, 301) z(i,i_CPU)  ! CPU time
        write(ife, 300)


        if( mod(i, iline) == 0 )  write(ife,*) "\hline"

     end do
     write(ife, *) "\hline"



  else if(typ == 'AMAtdp'  ) then
     ncounts = 7   ! number of columns in tab.tex
     logEOC  = .true.
     iEOC = 5
     iline = 4

     text(1) = "$\|e_h\|_{L^\infty(L^2)}$"
     text(2) = "$\|e_h\|_{L^2(L^2)}$"
     text(3) = "$\|e_h\|_{L^2(H^1)}$"
     text(4) = "$\| E_I\|_{L^q}$"
     text(5) = "$\| E_I\|_{L^\infty}$"
     text(6) = "$T$"
     text(7) = "CPU(s)"

     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"

     !write(ife_b, *) "\begin{tabular}{cccc|c|cc|cc|cc|cc|cc|cr}"
     !write(ife_b, *) "\hline"
     !write(ife_b, *) " $h$ & $\tau $ &  $P_k$ & $T_m$  & DOF "

     write(ife_b, *) "\begin{tabular}{cccc|c|c|cc|cc|cc|cc|cc|cr}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " $h$ & $\tau $ &  $P_k$ & $T_m$ &ndimL & DOF "


     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh

        write(ife, 99) x(i,2), x(i,3)  ! h_aver, tau
        write(ife, 100) ipar(i, 1:2)   ! P_k, T_m, ..
        write(ife, 95) ipar(i,3)   !  Q_k
        write(ife, 95) idof(i)   !  DOF,

        if(i > 1) then
           write(ife, 952) y(i,i_errL8L2), ord(i,i_errL8L2)
           write(ife, 852) y(i,i_errL2L2), ord(i,i_errL2L2)
           write(ife, 852) y(i,i_errL2H1), ord(i,i_errL2H1)
           write(ife, 852) y(i,i_interLq) , ord(i,i_interLq)
           write(ife, 852) y(i,i_interL8) , ord(i,i_interL8)

        else
           write(ife, 992) y(i,i_errL8L2)
           write(ife, 892) y(i,i_errL2L2)
           write(ife, 892) y(i,i_errL2H1)

           write(ife, 892) y(i,i_interLq)
           write(ife, 892) y(i,i_interL8)

        endif


        write(ife, 600) x(i,10)  ! state%time%ttime
        write(ife, 301) z(i,i_CPU)  ! CPU time
        !write(ife, 300)
        if( mod(i, iline) == 0 )  write(ife,*) "\hline"

     end do
     write(ife, *) "\hline"

  else if(typ == 'AMAsteady'  ) then
     ncounts = 11  ! number of columns in tab.tex
     !logEOC  = .true.
     logEOC  = .false.
     iEOC = 1
     iline = 12

     text(1) = "$\|u - u_h \|_{L^2}$"
     text(2) = "$|u - u_h|_{H^1}$"
     text(3) = "$\|u_h - \tilde{u}_h\|_{L^2}$"
     text(4) = "$|u_h - \tilde{u}_h|_{H^1}$"
     text(5) = "$\| E_I\|_{L^2}$"
     !text(3) = "$\|u - \Pi_h u\|_{L^\infty}$"
     text(6) = "$\| E_I\|_{L^\infty}$"
     !text(3) = "$\frac{\| E_I\|_{L^2}}{\|u - \Pi_h u\|_{L^2}}  $"
     !text(6) = "$\frac{\| E_I\|_{L^\infty}}{\|u - \Pi_h u\|_{L^\infty}}  $"
     text(7) = "$| E_I|_{H^1}$"

     text(8) = "$i^{\mathrm{ho}}_{L^2}$"
     text(9) = "$i^{\mathrm{ho}}_{H^1}$"
     text(10) = "$i^{\mathrm{int}}_{L^2}$"
     text(11) = "$i^{\mathrm{int}}_{H^1}$"
     !text(6) = "$T$"
     !text(12) = "CPU(s)"

     write(ife_b, *) "\tabcolsep 2.0pt"

     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"


     write(ife_b, *) "\begin{tabular}{rr|cc|cc|ccc|cccc}" !c|c}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " $\#\Thp $ &  $\Nhp$ "


     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh

        !write(ife, '(es9.1)')  z(i,i_tol)   ! TOLERANCE

        !write(ife, 99) x(i,2), x(i,3)  ! h_aver, tau
        !write(ife, 100) ipar(i, 1:2)   ! P_k, T_m, ..
        !write(ife, 95) ipar(i,3)   !  Q_k
        write(ife, 78) ix(i)   !  nelem
        write(ife, 95) idof(i)   !  DOF,

        !if(i > 1) then
        !   write(ife, 852) y(i,i_errL2), ord(i,i_errL2)*2 !!!!!!!!!!!!!!!!!!!!!!!!!
        !   !write(ife, 852) y(i,i_interLq) , ord(i,i_interLq)
        !   !write(ife, 852) y(i,i_errL8), ord(i,i_errL8)
        !   !write(ife, 852) y(i,i_errL2L2), ord(i,i_errL2L2)
        !   !write(ife, 852) y(i,i_errL2H1), ord(i,i_errL2H1)
        !   !write(ife, 852) y(i,i_interL8) , ord(i,i_interL8)
        !
        !else
        !   write(ife, 892) y(i,i_errL2)
           !write(ife, 892) y(i,i_interLq)

           !write(ife, 892) y(i,i_errL8)
           !write(ife, 892) y(i,i_errL2L2)
           !write(ife, 892) y(i,i_errL2H1)

           !write(ife, 892) y(i,i_interL8)

        !endif
        !write(ife, 281) y(i,i_interLq)/  y(i,i_errL2)
        !write(ife, 281)   y(i,i_errL2) /z(i,i_tol)

        write(ife, 281)   y(i,i_errL2)
        write(ife, 281)   y(i,i_errH1)

        write(ife, 281)   y(i,i_HO_estim_L2_p2)
        write(ife, 281)   y(i,i_HO_estim_H1_p2)

        write(ife, 281) y(i,i_interLq)
        write(ife, 281) y(i,i_interL8)
        write(ife, 281) y(i,i_interH1)

        write(ife, 501) y(i,i_HO_estim_L2_p2)/  y(i,i_errL2)
        write(ife, 501) y(i,i_HO_estim_H1_p2)/  y(i,i_errH1)

        write(ife, 501) y(i,i_interLq)/  y(i,i_errL2)
        write(ife, 501) y(i,i_interH1)/  y(i,i_errL2)


        !write(ife, 600) x(i,10)  ! state%time%ttime
        !write(ife, 301) z(i,i_CPU)  ! CPU time
        write(ife, 300)
        if( mod(i, iline) == 0 )  write(ife,*) "\hline"

     end do
     write(ife, *) "\hline"

  else if(typ == 'pNeu' .or. typ == 'pNeuP' .or. typ == 'pNeuAM' .or. typ == 'pNeuAMa' &
       .or. typ == 'pNeuHP'  ) then
     !ncounts = 18   ! number of columns in tab.tex
     ncounts = 16   ! number of columns in tab.tex
     logEOC  = .false.
     iEOC = 10
     !iline = 100000
     iline = 6
      if(typ == 'pNeuAM' ) then
         ncounts = 11

         text(1) = "$\|G(u-u_h)\|$"
         text(2) = "$\|u-u_h\|_{\mathrm{J}}$"
         text(3) = "$|u-u_h|_{\mathrm{DG}}$"
         text(4) = "$\|\nabla u_h + \sigma\|$"
         text(5) = "$\frac{h_K}{\pi}\|f-\nabla\cdot\sigma\|$"
         text(6) = "$\|G(u_h-s_h)\|$"
         text(7) = "$\eta_{\mathrm{BC}}$"
         text(8) = "$\eta_{tot}$"
         text(9) = "$\eta_{\mathrm{DG}}$"
         text(10) = "$I^{\mathrm{eff}}$"
         text(11) = "$I^{\mathrm{eff}}_{\mathrm{DG}}$"

      elseif(typ == 'pNeuAMa' ) then
         ncounts = 11

         text(1) = "$\norm{\GGG(\uu \! - \! \uh)}$"
         text(2) = "$\norm{\uu \! - \! \uh}_{\mathrm{J}}$"
         text(3) = "$\norm{\uu \! - \! \uh}_{\mathrm{DG}}$"
         !text(4) = "$\norm{\GGG \uh \! + \! \frh}$"
         !text(5) = " $\norm{\GGG (\uh \! - \! \prh)}$"
         text(4) = " $\eta_{\mathrm{flux}}$"
         text(5) = " $\eta_{\mathrm{osc}}$"
         text(6) = "$\eta_{\mathrm{pot}}$"
         text(7) = "$\eta_{\mathrm{BC}}$"
         text(8) = " $\eta$"
         text(9) = " $\eta_{\mathrm{DG}}$"
         text(10) = " $I^{\mathrm{eff}}$"
         text(11) = " $I^{\mathrm{eff}}_{\mathrm {DG}}$"
      elseif(typ == 'pNeuP' ) then
         ncounts = 7

         !text(1) = "$\|e_h\|_0$"
         text(1) = "$|e_h|_1$"
         text(2) = "$\eta_{\mathrm{F}}$"
         text(3) = "$\eta_{\mathrm{rez}}$"
         text(4) = "$\eta_{\mathrm{pot}}$"
         text(5) = "$\eta_{\mathrm{BC}}$"
         text(6) = "$\eta_{\mathrm{tot}}$"
         text(7) = "$i_{\mathrm{eff}}$"

      elseif(typ == 'pNeuHP' ) then
         ncounts = 7
         iline = 1000
         !text(1) = "$\|e_h\|_0$"
         text(1) = "$\|\nabla(u-u_h)\|$"
         text(2) = "$\eta_{\mathrm{F}}$"
         text(3) = "$\eta_{\mathrm{rez}}$"
         text(4) = "$\eta_{\mathrm{pot}}$"
         text(5) = "$\eta_{\mathrm{BC}}$"
         text(6) = "$\eta_{\mathrm{tot}}$"
         text(7) = "$I^{\mathrm{eff}}$"
         !text(8) = "$\eta_{\mathrm{tot}}^{p-1}$"

      else

         text(1) = "$\|e_h\|_0$"
         text(2) = "$|e_h|_1$"
         text(3) = "$\|\nabla u_h + \sigma\|$"
         text(4) = "$\frac{h_K}{\pi}\|f-\nabla\cdot\sigma\|$"
         text(5) = "$\|\nabla(u_h-s_h)\|$"
         text(6) = "$\eta_{BC}$"
         text(7) = "$\eta_{tot}$"
         text(8) = "$i_{\mathrm{eff}}$"

         text(9) = "$\eta_{tot}^{p-1} $"

         text(10) = "$\|\nabla s_h + \sigma\|$"
         text(11) = "$\|\nabla (s_h -u) \|$"
         text(12) = "$\|\nabla u + \sigma\|$"

         text(13) = "$\| \Pi_{-1}\sigma\|$"
         text(14) = "$\| \Pi_{-2}\sigma\|$"
         text(15) = "$\| \Pi_{-3}\sigma\|$"
         text(16) = "$\| \Pi_{-4}\sigma\|$"
         ! text(12) = "$\| \Pi_{-1}\nabla s\|$"
         ! text(13) = "$\| \Pi_{-2}\nabla s\|$"
         ! text(14) = "$\| \Pi_{-3}\nabla s\|$"
         ! text(15) = "$\| \Pi_{-4}\nabla s\|$"
         text(17) = "$\| \Pi_{-1} s\|$"
         text(18) = "$\| \Pi_{-2} s\|$"
         text(19) = "$\| \Pi_{-3} s\|$"
         text(20) = "$\| \Pi_{-4} s\|$"
      endif

      if(typ == 'pNeuHP' .or. typ == 'pNeuP') then
         write(ife_b, *) "{\footnotesize           %  vetsi"
         write(ife_b, *) "%{\scriptsize             %  mensi"
         write(ife_b, *) "%{\tiny                   %  nejmensi"
         write(ife_b, *) "%{                   %  normal"
      else
         write(ife_b, *) "%{\footnotesize           %  vetsi"
         write(ife_b, *) "%{\scriptsize             %  mensi"
         write(ife_b, *) "{\tiny                   %  nejmensi"
      endif

     if(typ == 'pNeuP') then
        write(ife_b, *) "\begin{tabular}{|ccr|c|ccccc|r|}"

     elseif(typ == 'pNeuHP') then
        !write(ife_b, *) "\begin{tabular}{|rrr|c|ccccc|r|c|}"
        write(ife_b, *) "\begin{tabular}{|rrr|c|ccccc|r|}"

     elseif(typ == 'pNeuAM') then
        write(ife_b, *) "\begin{tabular}{|cc|ccc|cccc|cc|rr|}"

     elseif(typ == 'pNeuAMa') then
        write(ife_b, *) " \begin{tabular}{|cc|ccc|cccc|cc|rr|} "

     else
        !write(ife_b, *) "\begin{tabular}{|ccr|cc|ccccc|r|cccc|cccc|}"
        write(ife_b, *) "\begin{tabular}{|ccr|cc|ccccc|r|c|ccc|cccc|cccc|}"
     endif

     write(ife_b, *) "\hline"

     if(typ == 'pNeuHP') then
        write(ife_b, *) "lev &  $\#\Th$ &  DoF "

     elseif(typ == 'pNeuAMa') then
        write(ife_b, *) " $h$ &  $p$   "

     elseif(typ /= 'pNeuAM') then
        write(ife_b, *) " $N$ &  $P_k$  & dof "

     elseif(typ == 'pNeuAM') then
        write(ife_b, *) " $N$ &  $P_k$   "

     else
        write(ife_b, *) " $h$ &  $P_k$  "
     endif

     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh

        !write(ife, 99) x(i,2), x(i,3)  ! h_aver, tau
        !write(ife, 100) ipar(i, 1:2)   ! P_k, T_m, ..

        !write(ife, 88) xi(k)/2**0.5   ! h ???
        !write(ife, 98) x(i,2)   ! h

        if(typ == 'pNeuHP') then
           write(ife, '(i4, a4)') i-1,' & '
           write(ife, 78) ix(i)   ! nelem
           write(ife, 95) idof(i)

        elseif(typ == 'pNeuP') then
           write(ife, 78) ix(i)   ! nelem
           write(ife, 90) ipar(i, 1:1)           ! P_k, ..
           write(ife, 95) idof(i)           ! dof

        else

           if(typ == 'pNeuAMa') then
              ii = 2**(i-1)
              if(i==1)  write(ife, *) '$h_0$ '
              if(i>1)  write(ife, '(a25,i3,a4 )') '$\approx{h_0/}{',ii,'}$ '
              !write(ife, 98) x(i,2)
              write(ife, 90) ipar(i, 1:1)
           else

              write(ife, 78) ix(i)   ! nelem
              if(typ /= 'pNeuHP') write(ife, 90) ipar(i, 1:1)           ! P_k, ..
              if(typ /= 'pNeuAM') write(ife, 95) idof(i)           ! dof

              write(ife, 78) ix(i)   ! nelem
              if(typ /= 'pNeuHP') write(ife, 90) ipar(i, 1:1)           ! P_k, ..
              if(typ /= 'pNeuAM') write(ife, 95) idof(i)           ! dof

           endif
        endif

        !  write(*,'(40i5)') i_errL2, i_errH1, i_P_flux, i_P_rez, i_P_pot, i_P_tot, -999, &
        !       i_P_sF, i_P_su, i_P_FDu, i_P_F_p1, i_P_F_p2, i_P_F_p3, i_P_F_p4, i_P_s_p1, i_P_s_p2, i_P_s_p3, i_P_s_p4

        if(typ == 'pNeu') write(ife, 201) y(i,i_errL2)
        !write(ife, 201) y(i,i_errH1_discrete)
        write(ife, 201) y(i,i_errH1)

        if(typ == 'pNeuAM' .or. typ == 'pNeuAMa') then
           write(ife, 201) y(i,i_EjumpsPi)
           write(ife, 201) y(i,i_DG)
        endif
        write(ife, 201) y(i,i_P_flux)
        write(ife, 201) y(i,i_P_rez)
        write(ife, 201) y(i,i_P_pot)
        !if(typ /= 'pNeuAM')
        write(ife, 201) y(i,i_P_BC)
        write(ife, 201) y(i,i_P_tot)
        if(typ == 'pNeuAM'.or. typ == 'pNeuAMa')  write(ife, 201) y(i,i_etaDG)

        !write(ife, 701) y(i,i_P_tot) / max(1E-15, y(i,i_errH1_discrete))    ! i_eff
        write(ife, 701) y(i,i_P_tot) / max(1E-15, y(i,i_errH1))    ! i_eff

        if(typ == 'pNeuAM' .or. typ == 'pNeuAMa') &
             write(ife, 701) y(i,i_etaDG)/ y(i,i_DG)

        !if(typ == 'pNeuHP') write(ife, 201) y(i,i_P_potP)

        if(typ == 'pNeu') then
           write(ife, 201) y(i,i_P_sF)
           write(ife, 201) y(i,i_P_su)
           write(ife, 201) y(i,i_P_FDu)
           write(ife, 201) y(i,i_P_F_p1)
           write(ife, 201) y(i,i_P_F_p2)
           write(ife, 201) y(i,i_P_F_p3)
           write(ife, 201) y(i,i_P_F_p4)
           write(ife, 201) y(i,i_P_s_p1)
           write(ife, 201) y(i,i_P_s_p2)
           write(ife, 201) y(i,i_P_s_p3)
           write(ife, 201) y(i,i_P_s_p4)
        endif

        write(ife, 300)

        if(i > 1 .and. typ /= 'pNeuAMa' .and. typ /= 'pNeuHP') then   ! FOR SHORTER TABLE fo SISC
           if( typ == 'pNeuHP') then
              write(ife,*) '  \multicolumn{2}{|c|}{{\footnotesize (EOC)}} '

           elseif(typ /= 'pNeuAM') then

              write(ife,*) '  \multicolumn{3}{|c|}{{\scriptsize (EOC)}} '
              if(typ == 'pNeu') write(ife, 621) ord(i,i_errL2)

           elseif(typ == 'pNeuAM') then

              write(ife,*) '  \multicolumn{2}{|c|}{{\scriptsize (EOC)}} '
           else

              write(ife,*) '  \multicolumn{3}{|c|}{{\scriptsize (EOC)}} '
           endif

           iform = 621
           if( typ == 'pNeuP' .or. typ == 'pNeuHP') then
              !iform = 671
              write(ife, 621) ord(i,i_errH1_discrete)
              write(ife, 621) ord(i,i_P_flux)
              write(ife, 621) ord(i,i_P_rez)
              write(ife, 621) ord(i,i_P_pot)
              write(ife, 621) ord(i,i_P_BC)
              write(ife, 621) ord(i,i_P_tot)
              write(ife, *) ' &  '
              if( typ == 'pNeuHP' ) write(ife, 671) ord(i,i_P_potP)

           else
              write(ife, 621) ord(i,i_errH1_discrete)
              if(typ == 'pNeuAM') then
                 write(ife, 621) ord(i,i_EjumpsPi)
                 write(ife, 621) ord(i,i_DG)
              endif

              write(ife, 621) ord(i,i_P_flux)
              write(ife, 621) ord(i,i_P_rez)
              write(ife, 621) ord(i,i_P_pot)
              !if(typ /= 'pNeuAM')
              write(ife, 621) ord(i,i_P_BC)
              write(ife, 621) ord(i,i_P_tot)
              if(typ == 'pNeuAM')  write(ife, 621) ord(i,i_etaDG)
              write(ife, *) ' & '
              if(typ == 'pNeuAM' )  write(ife, *) ' & '

              if( typ /= 'pNeuHP' .and. typ /= 'pNeuAM') then
                 write(ife, 621) ord(i,i_P_potP)
                 write(ife, 621) ord(i,i_P_sF)
                 write(ife, 621) ord(i,i_P_su)
                 write(ife, 621) ord(i,i_P_FDu)
              endif

              if(typ == 'pNeu') then
                 write(ife, 621) ord(i,i_P_F_p1)
                 write(ife, 621) ord(i,i_P_F_p2)
                 write(ife, 621) ord(i,i_P_F_p3)
                 write(ife, 621) ord(i,i_P_F_p4)
                 write(ife, 621) ord(i,i_P_s_p1)
                 write(ife, 621) ord(i,i_P_s_p2)
                 write(ife, 621) ord(i,i_P_s_p3)
                 write(ife, 621) ord(i,i_P_s_p4)
              endif
           endif


        write(ife, 300)

        endif

        !write(ife, 600) x(i, 31) / x(i,27) ! i_{eff}
        !write(ife, 600) x(i, 16) / x(i,7) ! i_{eff}
        !write(ife, 600) (x(i, 16)**2 + x(i,8)**2)**0.5 / (x(i,7)**2 + x(i,8)**2)**0.5 ! i_{eff}

        !write(ife, 301) z(i,i_CPU)  ! CPU time

        if( mod(i, iline) == 0 )  write(ife,*) "\hline"

     end do
     write(ife, *) "\hline"

  else if( typ == 'HO_rec'  ) then
     ncounts = 6   ! number of columns in tab.tex
     logEOC  = .false.
     iEOC = 12
     iline = 500
     ordA = .false.
     !ordA = .true.
     !if(ordA) ncounts = 16
     if(ordA) iline = 100

     text(1) = "$\|e_h\|_{L^2}$"
     text(2) = "$\eta_{L^2}(u_h)$"
     !text(3) = "$\eta_{L^2}(\Pi_h u)$"
     text(3) = "$i_{L^2}$"
     !text(5) = "$\eta(u_h)/\eta(\Pi u_h)$"

     text(4) = "$\|e_h\|_{H^1}$"
     text(5) = "$\eta_{H^1}(u_h)$"
     !text(8) = "$\eta_{H^1}(\Pi_h u)$"
     text(6) = "$i_{H^1}$"
     !text(10) = "$\eta(u_h)/\eta(\Pi u_h)$"

     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"

     write(ife_b, *) "\begin{tabular}{|cc|ccc|ccc|}"
     write(ife_b, *) "\hline"
     if( typ == 'HO_rec') write(ife_b, *) " $h$  &  $p$  "

     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh

        !write(ife, 99) x(i,2), x(i,3)  ! h_aver, tau
        !write(ife, 98) x(i,2)  ! h_aver, tau
        write(ife,'(a6,i6,a6)') '$ 1/ ',int( (x(i,2)/2**0.5 )**(-1)+0.5),'$ '
        write(ife, 90) ipar(i, 1)   ! P_k, T_m, ..

        write(ife, 282) y(i,i_errL2), y(i,i_HO_estim_L2_p2) !,y(i,i_HO_estim_L2_p2)
        write(ife, 521)     y(i,i_HO_estim_L2_p2) / y(i,i_errL2)
        !write(ife, 521)     y(i,i_HO_estim_L2_p2) / y(i,i_errL2)
        !write(ife, 281)  y(i,i_HO_estim_L2_p1)/y(i,i_HO_estim_L2_p2)
        !write(ife, 282) y(i,i_HO_trunc_L2_p1),y(i,i_HO_trunc_L2_p2)

        write(ife, 282) y(i,i_errH1), y(i,i_HO_estim_H1_p2) !,y(i,i_HO_estim_H1_p2)
        write(ife, 521)     y(i,i_HO_estim_H1_p2) / y(i,i_errH1)
        !write(ife, 521)     y(i,i_HO_estim_H1_p2) / y(i,i_errH1)
        !write(ife, 281)  y(i,i_HO_estim_H1_p1)/y(i,i_HO_estim_H1_p2)
        !write(ife, 282) y(i,i_HO_trunc_H1_p1),y(i,i_HO_trunc_H1_p2)

        write(ife, 300)

        if(i >1 ) then
           write(ife,*) ' \multicolumn{2}{|c|}{{\scriptsize (EOC)}} '
           write(ife, 622) ord(i,i_errL2), ord(i,i_HO_estim_L2_p2) !,ord(i,i_HO_estim_L2_p2)
           write(ife,*) ' &  '
           !write(ife, 622) ord(i,i_HO_trunc_L2_p1),ord(i,i_HO_trunc_L2_p2)

           write(ife, 622) ord(i,i_errH1), ord(i,i_HO_estim_H1_p2)
           write(ife,*) ' &'
           !write(ife, 622) ord(i,i_HO_trunc_H1_p1),ord(i,i_HO_trunc_H1_p2)

           write(ife,*) '  \\'
        endif


        if( mod(i, iline) == 0 )  write(ife,*) "\hline"

     end do
     write(ife, *) "\hline"


  else if( typ == 'HO_recC' .or. typ == 'HO_recA'  ) then
     if(typ == 'HO_recA') ncounts = 4   ! number of columns in tab.tex
     if(typ == 'HO_recC') ncounts = 4   ! number of columns in tab.tex
     !logEOC  = .false.
     logEOC  = .true.
     iEOC = 2
     iline = 500
     ordA = .false.
     !ordA = .true.
     !if(ordA) ncounts = 16
     if(ordA) iline = 100

     !text(1) = "$\|e_h\|_{L^2(\Omega)}$"
     !text(2) = "$\|\Eh \|_{L^2(\Omega)}$"
     !text(3) = "$i_{L^2}$"

     text(1) = "$|e_h|_{H^1(\Omega,\Th)}$"
     text(2) = "$|\Eh|_{H^1(\Omega,\Th)}$"
     !text(3) = "$\|u_h\|_{J}$"
     text(3) = "$\ieff$"
     text(4) = "{\scriptsize\%cpu}"
     !text(5) = "$\ieffDG$"
     !text(3) = "$i_{H^1}$"
     !text(4) = "$|e_h|_{H^1(\Omega,\Th)}$"
     !text(5) = "$|\Eh|_{H^1(\Omega,\Th)}$"
     !text(6) = "$i_{H^1}$"

     !text(10) = "$\eta(u_h)/\eta(\Pi u_h)$"
     !text(11) = "$|R_h u_h-u_h|$"
     !text(9) = "$|R_h u_h - \Pi^{p-1} u_h|$"
     !text(10) = "$|R_h u_h - \Pi^{p-2} u_h|$"
     !text(14) = "$i'_{H^1}$"

     write(ife_b, *) "%{                   %  normal"
     write(ife_b, *) "{\footnotesize           %  vetsi"
     write(ife_b, *) "%{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"

     if( typ == 'HO_recC') then
        if(.not. logEOC )   write(ife_b, *) "\begin{tabular}{|cc|ccc|ccc|}"
        if(      logEOC )   write(ife_b, *) "\begin{tabular}{|cc|cr|cr|c|r|}"
        write(ife_b, *) "\hline"
        write(ife_b, *) " $h$  &  $p$  "

     else if( typ == 'HO_recA') then

        if(.not. logEOC )  write(ife_b, *) "\begin{tabular}{|ccr|ccc|ccc|}"
        if(      logEOC )  write(ife_b, *) "\begin{tabular}{|ccr|cr|cr|c|r|}"
        write(ife_b, *) "\hline"
        write(ife_b, *) " lev &  $\#\Th$  &  DOF  "

     endif


     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh

        if( typ == 'HO_recC') then
           !write(ife, 99) x(i,2), x(i,3)  ! h_aver, tau
           !write(ife, 98) x(i,2)  ! h_aver, tau
           write(ife,'(a6,i6,a6)') '$ 1/ ',int( (x(i,2)/2**0.5 )**(-1)+0.5),'$ '
           write(ife, 90) ipar(i, 1)   ! P_k
        else if( typ == 'HO_recA') then
           write(ife, '(i4, a4)') i-1,' & '
           write(ife, 78) ix(i)   ! nelem
           write(ife, 95)  idof(i)  ! DOF
        endif

        if(logEOC) then

            if(i == 1) write(ife, 992) y(i,i_errH1)
            if(i >  1) write(ife, 902) y(i,i_errH1) ,  ord(i,i_errH1)

            if(i == 1) write(ife, 992) y(i,i_HO_estim_H1_p2)
            if(i >  1) write(ife, 902) y(i,i_HO_estim_H1_p2) ,  ord(i, i_HO_estim_H1_p2)

            !if(i == 1) write(ife, 992) y(i,i_EjumpsG)
            !if(i >  1) write(ife, 902) y(i,i_EjumpsG) ,  ord(i, i_EjumpsG)

            write(ife, 501)     y(i,i_HO_estim_H1_p2) / y(i,i_errH1)
            write(ife, 571)     z(i,i_estim) / z(i,i_CPU) *100
            !write(ife, 521)     y(i,i_etaDG) / y(i,i_DG)


           !if(i == 1) write(ife, 992) y(i,i_DG)
           !if(i >  1) write(ife, 902) y(i,i_DG) ,  ord(i,i_DG)

           !if(i == 1) write(ife, 992) y(i,i_etaDG)
           !if(i >  1) write(ife, 902) y(i,i_etaDG) ,  ord(i, i_etaDG)

           !write(ife, 521)     y(i,i_etaDG) / y(i,i_DG)

           write(ife, 300)

        else
           write(ife, 282) y(i,i_errL2), y(i,i_HO_estim_L2_p2) !,y(i,i_HO_estim_L2_p2)
           write(ife, 521)     y(i,i_HO_estim_L2_p2) / y(i,i_errL2)
           !!write(ife, 521)     y(i,i_HO_estim_L2_p2) / y(i,i_errL2)
           !!write(ife, 281)  y(i,i_HO_estim_L2_p1)
           !!write(ife, 281)  y(i,i_HO_estim_L2_p0)
           !!write(ife, 283) y(i,i_HO_trunc_L2_p2),y(i,i_HO_trunc_L2_p1),y(i,i_HO_trunc_L2_p0)
           !!write(ife, 521)     y(i,i_HO_trunc_L2_p2) / y(i,i_errL2)

           write(ife, 282) y(i,i_errH1), y(i,i_HO_estim_H1_p2) !,y(i,i_HO_estim_H1_p2)
           write(ife, 521)     y(i,i_HO_estim_H1_p2) / y(i,i_errH1)
           !write(ife, 521)     y(i,i_HO_estim_H1_p2) / y(i,i_errH1)
           !write(ife, 281)  y(i,i_HO_estim_H1_p1)
           !write(ife, 281)  y(i,i_HO_estim_H1_p0)
           !write(ife, 281)  y(i,i_HO_estim_H1_p1)/y(i,i_HO_estim_H1_p2)
           !write(ife, 283) y(i,i_HO_trunc_H1_p2),y(i,i_HO_trunc_H1_p1),y(i,i_HO_trunc_H1_p0)
           !write(ife, 521)     y(i,i_HO_trunc_H1_p2) / y(i,i_errH1)

           write(ife, 300)

           if(i >1 ) then
              if( typ == 'HO_recC') then
                 write(ife,*) ' \multicolumn{2}{|c|}{{\tiny (EOC)}} '
              elseif( typ == 'HO_recA') then
                 write(ife,*) ' \multicolumn{3}{|c|}{{\tiny (EOC)}} '
              endif
              write(ife, 622) ord(i,i_errL2), ord(i,i_HO_estim_L2_p2) !,ord(i,i_HO_estim_L2_p2)
              write(ife,*) ' &  '
              !!write(ife, 622) ord(i,i_HO_estim_L2_p1), ord(i,i_HO_estim_L2_p0)
              !!write(ife, 623) ord(i,i_HO_trunc_L2_p2),ord(i,i_HO_trunc_L2_p1), ord(i,i_HO_trunc_L2_p0)
              !!write(ife,*) ' &  '

              write(ife, 622) ord(i,i_errH1), ord(i,i_HO_estim_H1_p2)
              write(ife,*) ' &'
              !write(ife, 623) ord(i,i_HO_trunc_H1_p2), ord(i,i_HO_trunc_H1_p1), ord(i,i_HO_trunc_H1_p0)
              !write(ife, 622) ord(i,i_HO_estim_H1_p1), ord(i,i_HO_estim_H1_p0)
              !write(ife,*) ' &  '

              write(ife,*) '  \\'
           endif

        endif

        if( mod(i, iline) == 0 )  write(ife,*) "\hline"

     end do
     write(ife, *) "\hline"





  else if( typ == 'HO_recB'  ) then
     ncounts = 9  ! number of columns in tab.tex
     logEOC  = .false.
     iEOC = 12
     iline = 500
     ordA = .false.
     !ordA = .true.
     !if(ordA) ncounts = 16
     if(ordA) iline = 100

     !text(2) = "$\eta_{L^2}^p$"
     !text(1) = "$\|e_h\|_{L^2}$"
     !text(2) = "$\eta_{L^2}(u_h)$"
     !text(3) = "$\theta_{L^2}^p$"
     !text(4) = "$i_{L^2}^p$"

     text(1) = "$\|e_h\|_{H^1}$"
     text(2) = "$\eta_{H^1}(u_h)$"
     text(3) = "$\theta_{H^1}^p$"
     text(4) = "$i_{H^1}^p$"
     text(5) = "$\ll$"
     text(6) = "$|u - Ru_h|_{H^1}$"
     text(7) = "$|u_h - Ru_h|_{H^1}$"
     text(8) = "$i_{R}$"
     text(9) = "$\ll$"



     write(ife_b, *) "{\footnotesize           %  vetsi"
     write(ife_b, *) "%{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"

     if(.not. ordA) then
        write(ife_b, *) "\begin{tabular}{|cc||ccc|cc|cc|c|c|}"
        write(ife_b, *) "\hline"
        write(ife_b, *) " $h$  &  $p$  "
     else
        write(ife_b, *) "\begin{tabular}{c|cc|cc|cc|cc|cc}"
        write(ife_b, *) "\hline"
        write(ife_b, *) " $c_A$ "
     endif

     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh

        if(.not. ordA) then
           !write(ife, 99) x(i,2), x(i,3)  ! h_aver, tau
           !write(ife, 98) x(i,2)  ! h_aver, tau
           !if( typ == 'HO_rec') then
              write(ife,'(a6,i6,a6)') '$ 1/ ',int( (x(i,2)/2**0.5 )**(-1)+0.5),'$ '
              write(ife, 90) ipar(i, 1)   ! P_k, T_m, ..
           !endif
           !if( typ == 'HO_recA')  write(ife, '(i5,a3,i8)' ) i, '&',  idof(i)
           !write(ife,'(a6,i6,a6)') '&$ 1/ ',int(  x(i,3)**(-1) +0.5),'$ '
           !!write(ife,*) '& adapt'

           !!write(ife, 97) x(i,3)   ! tau
           !write(ife, 100) ipar(i, 1:2)   ! P_k, T_m, ..
           !write(ife, 198) ipar(i, 1), idof(i)   ! P_k, DOF, ..


           !write(ife, 283) y(i,i_errL2), y(i,i_HO_estim_L2_p2),y(i,i_HO_trunc_L2_p2)
           !write(ife, 521) y(i,i_HO_estim_L2_p2) /  y(i,i_errL2)
           write(ife, 283) y(i,i_errH1), y(i,i_HO_estim_H1_p2),y(i,i_HO_trunc_H1_p2)
           write(ife, 521) y(i,i_HO_estim_H1_p2) /  y(i,i_errH1)
           write(ife, 281) y(i,i_HO_trunc_H1_p2) /  y(i,i_HO_estim_H1_p2)

           write(ife, 282) y(i,i_HO_recovery), y(i,i_HO_rec_estim)
           write(ife, 521) y(i,i_HO_rec_estim) /  y(i,i_errH1)
           write(ife, 281) y(i,i_HO_recovery)/y(i,i_HO_rec_estim)


           write(ife, 300)

        else
        endif

        if(i >1 .and. .not. ordA) then
           write(ife,*) ' \multicolumn{2}{|c||}{{\scriptsize (EOC)}} '
           !write(ife, 624) ord(i,i_errL2), ord(i,i_errH1), ord(i,i_HO_L2), ord(i,i_HO_H1)
           !write(ife, 623) ord(i,i_errL2), ord(i,i_HO_estim_L2_p2),ord(i,i_HO_trunc_L2_p2)
           !write(ife,*) ' & '
           write(ife, 623) ord(i,i_errH1), ord(i,i_HO_estim_H1_p2),ord(i,i_HO_trunc_H1_p2)
           write(ife,*) ' & &'
           write(ife, 622) ord(i,i_HO_recovery), ord(i,i_HO_rec_estim)
           write(ife,*) ' & '
           write(ife,*) ' & '
           !write(ife, 622) ord(i,i_HO_trunc_L2_p1),ord(i,i_HO_trunc_L2_p2)

           write(ife,*) '  \\'
        endif


        if( mod(i, iline) == 0 )  write(ife,*) "\hline"

     end do
     write(ife, *) "\hline"



  else if(typ == 'RTNst'  ) then
     ncounts = 8   ! number of columns in tab.tex
     logEOC  = .true.
     iEOC = 8
     iline = 1000000

     text(1) = "$ \|e_h\|_{\sigma} $"
   !  text(1) = "$\|e_h\|_{L^2(L^2)}$"
     text(2) = "$\|e_h\|_{L^2(H^1)}$"
     !text(3) = "$ \sqrt(\eta) + \sqrt(J(u))$"
     text(3) = "$ \eta $"
     text(4) = "$ \eta_{rez} $"
     text(5) = "$ \eta_{flux} $"
     text(6) = "$ \eta_{rad} $"
     text(7) = "$ \eta_{rad2} $"
     text(8) = "$ J(u) $"


     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"


     write(ife_b, *) "\begin{tabular}{cccc|cc|cc|cc|cc|cc|cc|cc|cc}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " $h$ & $\tau $ &  $P_k$ & $T_m$ "

     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh


        write(ife, 99) x(i,2), x(i,3)  ! h_aver, tau
        write(ife, 100) ipar(i, 1:2)   ! P_k, T_m, ..

        if(i > 1) then

           write(ife, 952) y(i,i_RTNstNorm), ord(i,i_RTNstNorm)
           !write(ife, 952) y(i,i_errL2L2), ord(i,i_errL2L2)
           write(ife, 952) y(i,i_errL2H1), ord(i,i_errL2H1)

           write(ife, 952) y(i,i_RTN_eta), ord(i,i_RTN_eta)
           write(ife, 952) y(i,i_RTN_rez), ord(i,i_RTN_rez)
           write(ife, 952) y(i,i_RTN_flux), ord(i,i_RTN_flux)
           write(ife, 952) y(i,i_RTN_radau), ord(i,i_RTN_radau) ! write(ife, 952) y(i,i_RTN_radau), ord(i,i_RTN_radau)
           write(ife, 952) y(i,i_RTN_radau2), ord(i,i_RTN_radau2)
           write(ife, 952) y(i,i_RTN_jump) , ord(i,i_RTN_jump)

        else
           write(ife, 992) y(i,i_RTNstNorm)
           !write(ife, 992) y(i,i_errL2L2)
           write(ife, 992) y(i,i_errL2H1)

           write(ife, 992) y(i,i_RTN_eta)
           write(ife, 992) y(i,i_RTN_rez)
           write(ife, 992) y(i,i_RTN_flux)
           write(ife, 992) y(i,i_RTN_radau) ! write(ife, 992) y(i,i_RTN_radau)
           write(ife, 992) y(i,i_RTN_radau2)
           write(ife, 992) y(i,i_RTN_jump)

        endif

        !write(ife, 600) x(i, 31) / x(i,27) ! i_{eff}
        !write(ife, 600) x(i, 16) / x(i,7) ! i_{eff}
        !write(ife, 600) (x(i, 16)**2 + x(i,8)**2)**0.5 / (x(i,7)**2 + x(i,8)**2)**0.5 ! i_{eff}

        !write(ife, 301) z(i,i_CPU)  ! CPU time
        write(ife, 300)
     end do
     write(ife, *) "\hline"


  else if(typ == 'DUA'  ) then
     ncounts = 6   ! number of columns in tab.tex
     logEOC  = .true.
     iEOC = 6
     iline = 1000000

     text(1) = "$\|e_h\|_{L^2(L^2)}$"
     text(2) = "$\|e_h\|_{L^2(H^1)}$"
     text(3) = "$ \eta $"
     text(4) = "$ \eta_{rez} $"
     text(5) = "$ \eta_{flux} $"
     text(6) = "$ \eta_{NC} $"



     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"

     write(ife_b, *) "\begin{tabular}{cccc|cc|cc|cc|cc|cc|cc|cc}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " $h$ & $\tau $ &  $P_k$ & $T_m$ "

     do i = 1, ncounts
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i)
        endif
     enddo
     write(ife_b, *) "\\   \hline"

     do i=imesh0,imesh

        write(ife, 99) x(i,2), x(i,3)  ! h_aver, tau
        write(ife, 100) ipar(i, 1:2)   ! P_k, T_m, ..

        if(i > 1) then

           write(ife, 952) y(i,i_errL2L2), ord(i,i_errL2L2)
           write(ife, 952) y(i,i_errL2H1), ord(i,i_errL2H1)

!           write(ife, 952) y(i,i_DUAeta) + y(i,i_DUAjump), ord(i,i_DUA_all)
           write(ife, 952) y(i,i_DUA_total), ord(i,i_DUA_total)
           write(ife, 952) y(i,i_DUA_rez), ord(i,i_DUA_rez)
           write(ife, 952) y(i,i_DUA_flux), ord(i,i_DUA_flux)
           write(ife, 952) y(i,i_DUA_NC), ord(i,i_DUA_NC)

        else
           write(ife, 992) y(i,i_errL2L2)
           write(ife, 992) y(i,i_errL2H1)

!           write(ife, 992) y(i,i_DUAeta) + y(i,i_DUAjump)
           write(ife, 992) y(i,i_DUA_total)
           write(ife, 992) y(i,i_DUA_rez)
           write(ife, 992) y(i,i_DUA_flux)
           write(ife, 992) y(i,i_DUA_NC)
        endif

        !write(ife, 600) x(i, 31) / x(i,27) ! i_{eff}
        !write(ife, 600) x(i, 16) / x(i,7) ! i_{eff}
        !write(ife, 600) (x(i, 16)**2 + x(i,8)**2)**0.5 / (x(i,7)**2 + x(i,8)**2)**0.5 ! i_{eff}

        !write(ife, 301) z(i,i_CPU)  ! CPU time
        write(ife, 300)
     end do
     write(ife, *) "\hline"

  else
     print*
     print*
     print*,'Unknown type of Table in ../SRC_O/Setorderx.f90'
     print*,'% type of table: RES STDGM  AMAst  AMAstO AMAtdp RES_STDGM pNeu pNeuP'
     print*
     print*
     stop

  endif


  write(ife_e, *) "\end{tabular}"
  write(ife_e, *) " } %%% end of font size  "

88 format(es9.1 )
97 format(' & 'es10.2 )
96 format(' & 'es9.1 )
98 format(es10.2 )
99 format(es10.2, ' & 'es10.2 )
89 format(es9.1, ' & 'es9.1 )

101 format(i8, '\hspace{2mm} & 'es12.3, '\hspace{2mm}' )
103 format(i8, ' & ' i5 )
111 format(i8, ' & 'es11.3 )
102 format(i8, ' & ', i8, ' & ',es9.1 )
151 format(i8, 2(' & 'es10.2) )
108 format(1(' & ', i3) )
100 format(2(' & ', i3) )
109 format(i3, ' & ', i3 )
80  format( i3 )
78  format( i8 )
90  format(1(' & ', i3) )
95  format(1(' & ', i8) )
198  format(2(' & ', i8) )

201 format(  1(' & ', es10.2) )
202 format(  2(' & ', es10.2) )
203 format(  3(' & ', es10.2) )
204 format(  4(' & ', es10.2) )
205 format(  5(' & ', es10.2) )
206 format(  6(' & ', es10.2) )
207 format(  7(' & ', es10.2) )
208 format(  8(' & ', es10.2) )
210 format( 10(' & ', es10.2) )
211 format( 11(' & ', es10.2) )
212 format( 12(' & ', es10.2) )
213 format( 13(' & ', es10.2) )

281 format(  1(' & ', es11.3) )
282 format(  2(' & ', es11.3) )
283 format(  3(' & ', es11.3) )
284 format(  4(' & ', es11.3) )
286 format(  6(' & ', es11.3) )
287 format(  7(' & ', es11.3) )

291 format(  1(' & ', es12.4) )
292 format(  2(' & ', es12.4) )


1282 format(  2(' & {\bf ', es11.3,'}') )

1201 format(  1(' & {\bf ', es10.2,'}') )
1202 format(  2(' & {\bf ', es10.2,'}') )
1203 format(  3(' & {\bf ', es10.2,'}') )
1204 format(  4(' & {\bf ', es10.2,'}') )

251 format(  1(' & ', es9.1) )
252 format(  2(' & ', es9.1) )
253 format(  3(' & ', es9.1) )
254 format(  4(' & ', es9.1) )
1251 format(  1(' & {\bf ', es9.1,'}') )
1252 format(  2(' & {\bf ', es9.1,'}') )
1254 format(  4(' & {\bf ', es9.1,'}') )

266 format( 3(' &', es12.4, '&',es9.1) )

1602 format(  2(' & {\bf ', f6.2,'}') )


301 format( '&', f8.1, '\\' )
309 format( '&', i8, '\\' )
300 format('\\')

571 format( 1(' & ', f7.1) )
501 format( 1(' & ', f5.2) )
504 format( 4(' & ', f5.1) )

521 format( 1(' & ', f8.3) )
522 format( 2(' & ', f8.3) )

600 format( 1(' & ',f6.2) )
601 format( 1(' & (', f6.2,')') )
602 format( 2(' & (', f6.2,')') )
603 format( 3(' & (', f6.2,')') )
604 format( 4(' & (', f6.2,')') )
605 format( 5(' & (', f6.2,')') )
606 format( 6(' & (', f6.2,')') )
607 format( 7(' & (', f6.2,')') )
608 format( 8(' & (', f6.2,')') )
610 format(10(' & (', f6.2,')') )
612 format(12(' & (', f6.2,')') )

621 format(  (' & {\tiny(', f7.2,')}') )
622 format( 2(' & {\tiny(', f7.2,')}') )
623 format( 3(' & {\tiny(', f7.2,')}') )
624 format( 4(' & {\tiny(', f7.2,')}') )
625 format( 5(' & {\tiny(', f7.2,')}') )
626 format( 6(' & {\tiny(', f7.2,')}') )
627 format( 7(' & {\tiny(', f7.2,')}') )

671 format(  (' & {\scriptsize(', f7.2,')}') )
672 format( 2(' & {\scriptsize(', f7.2,')}') )
673 format( 3(' & {\scriptsize(', f7.2,')}') )
674 format( 4(' & {\scriptsize(', f7.2,')}') )
675 format( 5(' & {\scriptsize(', f7.2,')}') )
676 format( 6(' & {\scriptsize(', f7.2,')}') )
677 format( 7(' & {\scriptsize(', f7.2,')}') )

701 format( 1(' & ', f8.2) )
704 format( 4(' & ', f8.2) )

721 format( 1(' & ', f10.4) )
741 format( 1(' & ', f12.6) )

752 format( ' &\hspace{2mm} ', es12.4, ' & \hspace{2mm}', f8.2,'\hspace{2mm}' )
792 format( ' &\hspace{2mm} ', es12.4, ' & \hspace{2mm}', ' -- \hspace{2mm}' )

742 format( ' & ', es11.1, ' & ', f8.2 )
782 format( ' & ', es11.1, ' & ', ' -- ' )


802 format(i8, ' & ', i8 )
803 format(i3, '&' , i8, ' & ', i8 )
804 format(i3, '&' , i8, ' & ', i8, '\ &', i8,'\quad ' )
8041 format(i3, '&' , '--  & ', i8, '\ &', i8,'\quad ')

813 format(i8, '&' , i8, ' & ', i8 )

852 format( ' & ', es11.3, ' & ', f7.1 )
892 format( ' & ', es11.3, ' & ', ' -- ' )

902 format( ' & ', es10.2, ' & ', f6.2 )
952 format( ' & ', es10.2, ' & ', f7.1 )
992 format( ' & ', es10.2, ' & ', ' -- ' )


999 format( '& & \hspace{0mm}', f8.2, '& & \hspace{0mm}', f8.2, '& & \hspace{2mm}', f8.2, '\hspace{2mm} \\ ' )


end program SetOrder

subroutine Regres(n, x, y, order )
  integer, intent(in) :: n     ! number of datas
  real, dimension(1:n), intent(in) :: x, y   ! arrays of data
  real, intent(out) :: order
  real :: sumx, sumy, sumx2, sumy2, sumxy
  real :: rb, ra, rr
  integer :: i

  sumx = sum( x(1:n) )
  sumy = sum( y(1:n) )
  sumx2 = dot_product( x(1:n), x(1:n) )
  sumy2 = dot_product( y(1:n), y(1:n) )
  sumxy = dot_product( x(1:n), y(1:n) )

  if( n*sumx2-sumx*sumx /= 0.) then
     rb = (n*sumxy -sumx*sumy)/(n*sumx2-sumx*sumx)
  else
     rb = -99999.
  endif

  !ra = (sumy - rb*sumx)/n
  !rr = (n*sumxy-sumx*sumy)/sqrt((n*sumx2-sumx*sumx)*(n*sumy2-sumy*sumy))
  !ra = exp(ra)

  order = rb

end subroutine Regres
