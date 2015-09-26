  program AD_Gen_test_f90
  character(len=128):: ini_file
  character(len=64):: chA, chB, chC, chD, chE, ch_model, time_disc
  character(len=64):: mesh
  integer::arg, ini,dat, len, imod
  integer, dimension(:), allocatable :: ibc
  real :: Ttime, Re, par, shock_capture
  integer :: idir, nibc
  
  real :: pi 

  pi = 2 * asin(1.) 

  dat =11
  open(dat,file='smaz',status='UNKNOWN')

  do arg=1, command_argument_count() 
     call get_command_argument(arg, chA)
     write(dat,*) chA
  enddo
  close(dat)

  open(dat,file='smaz',status='OLD')

  read(dat,*) idir    ! level of directories

  !name of the file
  read(dat,*) ini_file
  i = len_trim(ini_file)
  ini_file(i+1:i+4) ='.ini'

  ini = 10
  open(ini,file=ini_file,status='UNKNOWN')

  ! line model
  read(dat,*) ch_model
  read(dat,*) Re
  read(dat,*) imod
  read(dat,*) par

  IC = 1
  if(ch_model == 'NSe') then
     IC = imod
     imod = 0
  endif

  !write(ini,*) '123456789*123456789*123456789*123456789|123456789*123456789*123456789*123456789*123456789'

  write(ini, '(a7, es12.4,i5, es12.4, a56)') &  ! ** after 40 chacaters'
       ch_model, Re, imod, par, &
       "** model: 'scalar/NSe/turb/wet_steam',Re, isca, par"
  
  ! line final time - steady state
  read(dat,*) Ttime 
  write(ini,'(es10.4,  a75)') Ttime,  &
       "** final Time, =0 or >1E+10 for steady state"

  ! stopping criteria for steady state (rez and coeff_diffs)
  read(dat,*) stopR
  read(dat,*) stopC

  !print*,stopR

  write(ini,'(es8.2, 3es10.2, a51)') stopR, stopC, 10*stopC, stopC, &
       "** stopping criteria for steady state: Rez, cDLM"


  ! grid file
  read(dat,*) chA
  i = len_trim(chA)
  mesh(1:10) = chA(1:10)

  if(idir == 0) then
     chB(1:10) ="'../Grids/"
     chB(10+1:10+i) = chA(1:i)
     chB(11+i:16+i) = ".grid'"
     do j=17+i, 37
        chB(j:j) = " "
     enddo
  else
     chB(1:13) ="'../../Grids/"
     chB(13+1:13+i) = chA(1:i)
     chB(14+i:19+i) = ".grid'"
     do j=20+i, 37
        chB(j:j) = " "
     enddo
  endif

  write(ini,'(i1, a2, a37,  a35)') 2, ' ', chB(1:37),  &
       "** space dim (2,3), file with grid"

  ! prof file
  read(dat,*) iP
  read(dat,*) chA
  i = len_trim(chA)
 
  if(idir == 0) then
     chB(1:19) = "'../Grids/profiles."
     chB(19+1:19+i) = chA(1:i)
     chB(20+i:20+i) = "'"
     do j=21+i, 37
        chB(j:j) = " "
     enddo
  else
     chB(1:22) = "'../../Grids/profiles."
     chB(22+1:22+i) = chA(1:i)
     chB(23+i:23+i) = "'"
     do j=24+i, 37
        chB(j:j) = " "
     enddo
  endif


  write(ini,'(i1, a2, a37,  a44)') iP, ' ', chB(1:37),  &
       "** Q_k degree of curved boundary, prof_file"


  write(ini, '(i2,a98)' ) IC,'   G0.rsol                             ** IC: 1(BC) 0(G0.rsol) (-1=resultsx) (-2=dgm.sol) >0(code)'

!  write(ini, *)'G0.rsol           ** IC file read when type of IC = 0A'  

  ! name of output files 
  write(ini,'(a4,a76)') "test","       ** name of output files (*.sol, *.conv)"

  ! space discretization
  read(dat,*) chA
  read(dat,*) cW
  read(dat,*) iP
  write(ini,'(a1,a3, es14.4, i3,  a56)') chA(1:1),"IPG",  cW, iP,   &
       "** {N}{I}{S}IPG, C_w , p: P_p approx"

  ! time discretization
  read(dat,*) time_disc
  read(dat,*) iP
  write(ini,'(a6, i3,  a82)') time_disc(1:6), iP,   &
       "** time discret: 'RK/SEMI/BDF/STDG', q: P_q approx"

  ! choice of the time step
  read(dat,*) chA
  read(dat,*) val

  if(chA == 'fixed' ) then
     write(ini,'(a6, es12.4,  a6, a76)') chA(1:6), val, '---',  &
          "** tau: fixed (tau), cfl (cfl), exp (cfl), adapt (tol)|type"
  else
     write(ini,'(a6, es12.4,  a6, a76)') chA(1:6), val, 'tRES',  &
          "** tau: fixed (tau), cfl (cfl), exp (cfl), adapt (tol)|type"
  endif

  read(dat,*) chA
  read(dat,*) tol_max
  read(dat,*) tol_min
  read(dat,*) tol_Lq

  write(ini,'(a6, 2es12.4, f6.2, a60)') chA, tol_max,tol_min, tol_Lq,  &
     '** %estim_space, state%tol_max, state%tol_min, state%Lq'


  read(dat,*) chA
  read(dat,*) i_adapt_iter

  write(ini,'(a6, i5, a6, a78)') chA(1:6), i_adapt_iter, 'y',   &
       "** mesh adapt 'RGh/RGp/RGhp/AMAh/AMA/hp', adapt levels"

  read(dat,*) max_iter

  if(max_iter == -1) then
     ! spatial case, only one time step, one Newton, tol GMRES small
     write(ini, '(i6, a81)') 1,'  ** maximum number of time steps for mesh level'

     read(dat,*) chA, tol
     
     write(ini, '(a6,a8,es12.4,a15,a51)' ) 'newton  ', chA(1:6), tol ,'   1    1   ', &
          '** non_solver, stop_crit, tol, max_iter, min C update'
     
     read(dat,*) chA
     write(ini, *) 'none      GMRES_ILU   1E-20               ', &
                '** linear MultiGrid?,  lin_solver   tol ' 

  elseif(max_iter == -2) then
     ! spatial case, only one time step, one Newton, tol GMRES small
     write(ini, '(i6, a81)') 1,'  ** maximum number of time steps for mesh level'

     read(dat,*) chA, tol
     
     write(ini, '(a6,a8,es12.4,a15,a51)' ) 'newton  ', chA(1:6), tol ,'   1    1   ', &
          '** non_solver, stop_crit, tol, max_iter, min C update'
     
     read(dat,*) chA
     write(ini, *) 'none      GMRES_ILU   1E-10               ', &
                '** linear MultiGrid?,  lin_solver   tol ' 

  else

     write(ini, '(i6, a81)') max_iter,'  ** maximum number of time steps for mesh level'

     read(dat,*) chA, tol
     
     write(ini, '(a6,a8,es12.4,a15,a51)' ) 'newton  ', chA(1:6), tol ,'   30    30   ', &
          '** non_solver, stop_crit, tol, max_iter, min C update'
          
     
     read(dat,*) chA, tol
     print*, toler
     if(chA /= "lin") then
        if(time_disc == 'STDG') then
           write(ini, '(a21, es8.2 , a52)') 'none      GMRES_ILU  F', tol , &
          '** linear MultiGrid?,  lin_solver   tol' 
        else
           write(ini, '(a21, es8.2 , a52)') 'none      GMRES_ILU  ', tol , &
          '** linear MultiGrid?,  lin_solver   tol' 
        endif
        
     else
        write(ini, '(a21, es8.2 , a52)') 'LINMG      GMRES_ILU  ', tol , &
          '** linear MultiGrid?,  lin_solver   tol' 
     endif
  endif

  if(TTime > 1E+10 .or. ch_model /= 'NSe') then
     write(ini, *)'0.     0                                ** output time, isol_start, initial value of sol* file'
  else
     write(ini, '(es12.4,a60)') Ttime/20,'     0                       ** output time, isol_start, initial value of sol* file'
  endif



  read(dat,*) nBC
  allocate(ibc(1: abs(nbc) ) )
  if(nbc > 0) then
     do i=1,nBC
        ibc(i) = i
     enddo
  else
     do i=1, abs(nBC)
        read(dat, *)  ibc(i) 
     enddo
     
  endif

  !write(*,'(a4,8i5)') '@@@',nbc, ibc(:)
  nbc = abs(nbc)


!  write(ini, *)' 0   20            ** NIPG = 1, IIPG = 0, SIPG = -1; sigma'

!write(ini, *) max_iter,' argv[13]   ** maximum number of additional iterations'
!  write(ini, *) i_adapt_iter, chb(1:5),'$argv[15]   $argv[14]                      ** maximum level of adaptations'
  ! write(ini, *) stopR, '$argv[16]                          ** stopping criterium for steady state'
  ! write(ini, *)stopC,'$tol                      ** cD'
  ! write(ini, *)stopC*10, '$tol2                  ** cL'
  ! write(ini, *)stopC,'$tol                           ** cM'
  !write(ini, *) fix_tau, '$argv[18]                           ** if nonzero, explicit time step '
!  write(ini, *)'0.5                            ** if nonzero, explicit tolerance for GMRES'
!  write(ini, *) CFL,'$argv[19]                           ** CFL number'
!  write(ini, *)BDFt,'$argv[20]                            ** tolerance for ABDF'
!  write(ini, *)nBC, '$argv[24]                              ** number of prescribed Inlet/Outlet BC'

!EOF

  write(ini, '(i3,a77)') nBC, '                              ** number of prescribed Inlet/Outlet BC'

  shock_capture = 0.
  !print*,'###',command_argument_count()
  if(ch_model == 'NSe' .and. command_argument_count() == 33) then

     read(dat,*) rmach
     read(dat,*) alpha
     read(dat,*) shock_capture
     
     rho = 1.
     
     v1 = cos(alpha/180* pi)
     v2 = sin(alpha/180* pi)
     p = 1./(1.4 * rmach * rmach)

     niBC = 1

     if(IC == 7 )then
        rho = 1.
        v1 = 1.
        v2 = 1.
        p = 1.
        niBC = 6
     endif

     
     write(ini,'(2i3,4es14.6,a15, f5.2,a8, f5.2)') niBC, 0, 1., v1,  v2,  p, &
          ' ** Inlet,  M = ', rmach, ', alpha =', alpha
     if(nbc > 1) &
          write(ini,'(2i3,4es14.6,a15, f5.2,a8, f5.2)') niBC+1, 1, 1., v1,  v2,  p, &
          ' ** Outlet,  M = ', rmach, ', alpha =', alpha
    
  elseif(ch_model == 'sca' .or. ch_model == 'scalar' ) then
    do i= 1, nBC
       write(ini,'(2i5, es12.4, a71)') ibc(i), 0, 1.,'   ** Dirichlet BC: wall index, Dirichlet, val (unused)'
    enddo
 else
    print*,'NON implemented !!!!!!!!!!!!!!!!!!!!11'
    print*,ch_model, command_argument_count()
 endif
 
! write(ini,* ) '2.0   0.   0.5    0.                    ** Stabilization parameters'
 write(ini,'(4f7.2,a40)' ) 2., shock_capture, 0.5, shock_capture,'    ** Stabilization parameters'
 write(ini, *) '  0.0    0.0    0.0    0                ** RTN error estims: gamma_rem, gamma_alg, gamma_lin, nu'

  close(ini)


  close(dat)


end program AD_Gen_test_f90

function iTrans_integer(string)
  integer:: iTrans_integer
  character(len=128), intent(in):: string 

  open(99,file='smaz',status='unknown')
  write(99,*) string(1: len_trim(string) )
  close(99)
  open(99,file='smaz',status='OLD')
  read(99,*) iTrans_integer
  close(99)

end function iTrans_integer

function Trans_real(string)
  real:: Trans_real
  character(len=128), intent(in):: string 

  open(99,file='smaz',status='unknown')
  write(99,*) string(1: len_trim(string) )
  close(99)
  open(99,file='smaz',status='OLD')
  read(99,*) Trans_real
  close(99)

end function Trans_real
