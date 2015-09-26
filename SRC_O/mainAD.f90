!> main program
program adgfem
  use mesh_oper
  use mesh_oper3D
  use problem_oper
  use euler_problem
  use solve_problem
  !FEM  use conforming_FEM
  use marking
  use mesh_adaptation
  use hp_adaptation
  !use angen
  use main_data
  use compute_oper
  use computeAD_oper
  !!!use color_figures
  use  helpers

  implicit none

  ! name of files
!  character(len=50) :: gridfile     ! file with mesh
!  character(len=50) :: solfile      ! solution with basis coefficients
!  character(len=50) :: rsolfile     ! solution in Lagrange nodes
!  character(len=50) :: convfile     ! histry of computation
!  character(len=50) :: prof_file    ! nodes lying on curved parts of boundary


  character(len=50) :: command_name, tri_name, sol_name, exa_name, err_name, est_name
  character(len=1) :: ch1
  character(len=3) :: ch3
!  character(len=5) :: ch5
  character(len=10) :: ch10
  character(len=20) :: chA
  integer :: i, i0, isca
  !> real time instants (CPU)
  real :: t2, val
  integer:: input
  character(len=128):: input_file

  call cpu_time(state%start_time)

  allocate(grid)

  !call TEST_ILU()
  !stop

  if (command_argument_count() > 0) then
    call get_command_argument(1,input_file)
    input = 10
    open(input,file=input_file,action='read')
  else
    input = 5
  endif

  ! READING OF PROBLEM DATA
  print*,' #            A D G F E M'
  write(*,*)' # ---------------------------------------------------------'

  state%SP = .false.

  ! model problem
  read(input,*) state%modelName, state%model%Re, isca, t2
  if(state%modelName == 'sca') state%modelName = 'scalar'

  if(state%modelName == 'NSe') then ! Euler or Navier-Stokes equations
     ndim = 4
     state%model%kappa = 1.4
     state%model%kappa1 = state%model%kappa - 1.
     state%model%Pr = 0.72   ! Prandtl number (almost all gases)

     if(state%model%Re == 0.) then
        print*,' # Compressible Euler equations'
        state%model%Re1 = 0.
     elseif(state%model%Re >0.) then
        state%model%Re1 = 1./state%model%Re
        write(*,'(a46, es12.4)')' # Compressible Navier-Stokes equations, Re = ',state%model%Re

     else
        print*,'# Reynolds number is negative',state%model%Re,' STABILIZATION !!!'
        !stop
     endif

  elseif (state%modelName == 'scalar') then  ! scalar case
     ndim = 1
     if(state%model%Re >0.) then
        state%model%Re1 = 1./state%model%Re
     else
        state%model%Re1 = 0.
     endif

     if(isca <= 0) then
        print*,' isca for scalar problem must be positive !'
        stop
     else
        allocate(state%scalar)
        state%model%icase = isca
        state%model%param1 = t2
        call InitScalarCase( )
        write(*,'(a25, es10.4, a10, i3, a8,es11.4 )') &
             '  # Scalar equation: ve= ',state%model%Re1, &
             ', icase = ',state%model%icase, &
             ', par = ',state%model%param1
     endif

  elseif (state%modelName == '2eqs') then  ! two equations
     ndim = 2
     if(state%model%Re >0.) then
        state%model%Re1 = 1./state%model%Re
     else
        state%model%Re1 = 0.
     endif

     if(isca <= 0) then
        print*,' isca for scalar problem must be positive !'
        stop
     else
        allocate(state%scalar)
        state%model%icase = isca
        state%model%param1 = t2
        call InitScalarCase( )
        write(*,'(a25, es10.4, a10, i3, a8,es11.4 )') &
             '  # Two equations: ve= ',state%model%Re1, &
             ', icase = ',state%model%icase, &
             ', par = ',state%model%param1
     endif

  elseif(state%modelName == 'wet_steam') then  ! scalar case
     ndim = 8
     state%model%kappa = 1.4
     state%model%kappa1 = state%model%kappa - 1.
     state%model%Pr = 0.72   ! Prandtl number (almost all gases)

     if(state%model%Re == 0.) then
        print*,' BAD setting in wet_steam'
        stop
     elseif(state%model%Re >0.) then
        state%model%Re1 = 1./state%model%Re
        print*,' # Simulation of wet steam flow, Re=',state%model%Re
     endif

  elseif(state%modelName == 'incNS') then ! incompressible  Navier-Stokes equations
     ndim = 3
     state%SP = .true.
     Bdim = 3
     state%model%kappa = 1.4
     state%model%kappa1 = state%model%kappa - 1.
     state%model%Pr = 0.72   ! Prandtl number (almost all gases)

  elseif(state%modelName == 'visNS') then ! incompressible  Navier-Stokes equations with memory
     ndim = 7
     state%SP = .true.

  else
     print*,' Model ',state%modelName,' is not implemented'
     stop
  endif

  ! time-dependent?, stopping final time
  read(input,*) state%time%FinTime
  if(state%time%FinTime < 1E+10 .or. state%time%FinTime == 0. )  then
     ! time dependent problem
     state%time_dependent = .true.
     write(*,'(a41,es12.4)') '  # Time dependent problem, final time = ', state%time%FinTime
  else
     ! stationary problem, we seek steady-state solution
     state%time_dependent = .false.
     state%time%FinTime = 1E+30
     write(*,'(a49,es12.4)') &
          ' # Steady-state problem, fictitious final time =', state%time%FinTime
  endif

  ! stopping tolerance for steady-state reziduum and cD, cL and cM
  read(input,*) state%conv_rez, state%EcD_tol, state%EcL_tol, state%EcM_tol


  ! space dimension (d=2,3), grid
  read(input,*) nbDim,  gridfile

  ! degree of approximation of curved parts of boundaries, name of file
  read(input,*) grid%curved_deg,  prof_file

  write(*,'(a4,i1,a8,a30, a7,i1, a20,a25)') &
       ' # ',nbDim,'D mesh: ', gridfile,'   Q_',grid%curved_deg, &
       ' boundary in file: ',prof_file


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! type of initial condition, =0 from file, =1 from BC
  read(input,*) state%type_IC, rsolfile
  !read(input,*) rsolfile
  if(state%type_IC == 1) then
     if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
        write(*,*) ' # IC taken from the exact solution'
     else
        write(*,*) ' # IC taken from BC'
     endif

  elseif(state%type_IC == 0) then
      write(*,*) ' # IC taken from file',rsolfile(1:len_trim(rsolfile))

  elseif(state%type_IC == -1) then
      write(*,*) ' # IC taken from file "resultsx" (piecewise constant)'

  elseif(state%type_IC == -2) then
      write(*,*) ' # IC taken from file "dgm.sol" (DG for visualization)'

  elseif(state%type_IC >0) then
      write(*,*) ' # IC given explicitly in subroutine SetOneElementIC(0) '

  else
     print*,'UNKNOWN type of IC'
     stop
  endif

  ! name of output file
  read(input,*) solfile
  i = len_trim(solfile)
  convfile = solfile
  solfile(i+1:i+4) = '.sol'
  convfile(i+1:i+5) = '.conv'

  print*,' # Output in files "',solfile(1:i+4),'" and "',convfile(1:i+5),'"'
  print*

  !! ###############################################################

  ! space approximation of the solution
  if(.not. state%SP) then
     read(input,*) state%space%disc_space, state%space%sigma, state%space%deg
     state%space%degP = -1
  elseif(state%modelName == 'incNS') then
      read(input,*) state%space%disc_space, state%space%sigma, state%space%deg, state%space%degP
  else
     stop 'UNKNOWN TYPE in mainAD.f90 (34)'
  endif


  if(state%space%disc_space == 'NIPG') then
     state%space%m_IPG = 1
  elseif(state%space%disc_space == 'IIPG') then
     state%space%m_IPG = 0
  elseif(state%space%disc_space == 'SIPG') then
     state%space%m_IPG = -1
  elseif(state%space%disc_space == 'FEM') then
     state%space%m_IPG = 5
     print*,'FEM NOT TESTED !!!'
  else
     print*,'Bad space discretization, only SIPG, NIPG, IIPG are implemented'
     stop
  endif
  if(state%space%deg < 0) then
     print*,'Default approximation degree has to be nonnegative ',state%space%deg
     stop
  endif
  if(state%space%deg < 10)  then
     write(*,'(a4,a24, a4,a7,es12.4, a4,i1)') &
          '  # ','Space discretization by ', &
          state%space%disc_space, ', c_W =',state%space%sigma, ", P_",state%space%deg
  else

  endif


  !! ###############################################################

  ! time approximation of the solution
  read(input,*) state%time%disc_time, state%time%deg

  state%time%cn = .false.  ! OLD techniques
  state%time%tdg = .false. ! OLD techniques
  state%time%stdgm  = .false.

  if(state%time%disc_time == 'RK') then
     if(state%time%deg /= 1 ) then
        print*,' Runge-Kutta method only with Tdeg  = 1 '
        stop
     endif
     state%time%time_method = 'E'    ! explicit
     if(state%time%deg <= 0) then
        print*,'Degree of time degree of BDF has to be >0 ',state%time%deg
        stop
     endif
     write(*,'(a45,i1)') '# Time discretization by RK method order = ',state%time%deg

  elseif (state%time%disc_time == 'SEMI') then
     state%time%time_method = 'S'    ! semi-implicit
     write(*,'(a53,i1)') &
          '# Semi-implicit time discretization by BDF, order = ',state%time%deg

  elseif (state%time%disc_time == 'BDF') then
     state%time%time_method = 'I'    ! implicit
     if(state%time%deg <= 0) then
        print*,'Degree of time degree of BDF has to be >0 ',state%time%deg
        stop
     endif
     write(*,'(a55,i1)') &
          ' # Fully implicit time discretization by BDF, order = ',state%time%deg


  elseif (state%time%disc_time == 'STDG') then
     state%time%stdgm  = .true.
     state%time%time_method = 'I'
     if(state%time%deg < 0) then
        print*,'Degree of time degree of ST DG  has to be >=0 ',state%time%deg
        stop
     endif

     write(*,'(a57,i1)') &
          ' # Fully implicit time discretization by ST DG, order = ',state%time%deg

  elseif(state%time%disc_time == 'CN') then ! Crank - Nicolson
     state%time%cn = .true.
     print*,'NON TESTED'

  elseif (state%time%disc_time == 'TDG') then ! extrapolated STDGM
     state%time%tdg = .true.
     print*,'NON TESTED'

  else
     print*,'Unknown time integration method: ',state%time%disc_time,' ??'
     stop
  endif

  !! ###############################################################


  ! choice of the time step
  read(input, *) state%time%tau_choice, val, state%time%estim_time
  if(state%time%tau_choice == 'fixed') then
      state%time%tau_fixed = .true.
      state%time%tau_new = val
      state%time%tau_fixed_size = val
      !if (state%time%tau_new <= 0.D+00) state%time%tau_fixed = .false.
      !state%time%tau_fixed_size = state%time%tau_new

    write(*,'(a25,es12.4)') &
          '  # Fixed time step tau = ',state%time%tau_new

   elseif(state%time%tau_choice == 'cfl') then
      state%time%CFL = val
    write(*,'(a20,es12.4)') &
          '  # Fixed CFL number = ',state%time%CFL

   elseif(state%time%tau_choice == 'exp') then
      state%time%CFL = val

    write(*,'(a40,es12.4)') &
          '  # Exponentially increasing  CFL number = ',state%time%CFL

   elseif(state%time%tau_choice == 'adapt') then
      state%time%BDFtol = val

      if(state%time%estim_time == 'loc') then
         write(*,'(a57,es12.4)') &
              ' # Adaptively chosen time step based on LOC estim: tol =',state%time%BDFtol
      elseif(state%time%estim_time == 'tRES') then
         write(*,'(a57,es12.4)') &
              ' # Adaptively chosen time step based on RES estim: tol =',state%time%BDFtol
      else
         print*,' Unknown type of the adaptive choice of the time step:  ',state%time%estim_time
         print*, ' Only techniques "loc" and "tRES" implemented !'
         stop
      endif
  else
     print*,' Unknown type of the choice of the time step:',state%time%tau_choice
     print*, 'Possibilities are: fixed, cfl, exp, adapt'
  endif

  ! (space) error estimation method
  read(input,*) state%space%estim_space, state%space%adapt%tol_max, state%space%adapt%tol_min, state%space%adapt%Lq
  write(*,'(a30,a6, a14, 3es9.2)') '  # Type of error estimation: ', state%space%estim_space, &
       ', tolerances:', state%space%adapt%tol_max, state%space%adapt%tol_min, state%space%adapt%Lq

  ! mesh adaptation method
  read(input,*) state%space%adapt%adapt_space, state%space%adapt%max_adapt_level, ch1

  if(state%space%adapt%adapt_space == 'HGhp' .or. state%space%adapt%adapt_space == 'HGh' &
       .or. state%space%adapt%adapt_space == 'HGp') then
     state%space%adapt%adapt_method = 'RES'
     state%space%adapt%adapt_type = 'HG'

  elseif(state%space%adapt%adapt_space == 'RGhp' .or. state%space%adapt%adapt_space == 'RGh' &
       .or. state%space%adapt%adapt_space == 'RGp') then
     state%space%adapt%adapt_method = 'RES'
     state%space%adapt%adapt_type = 'RG'

  elseif(state%space%adapt%adapt_space == 'AMAhp' .or. state%space%adapt%adapt_space == 'AMAh') then
     state%space%adapt%adapt_method = 'Ahp'
     state%space%adapt%adapt_type = 'Ahp'

  elseif(state%space%adapt%adapt_space == 'ANIhp' .or. state%space%adapt%adapt_space == 'ANIh') then
     state%space%adapt%adapt_method = 'ANI'
     state%space%adapt%adapt_type = 'Ahp'

  elseif(state%space%adapt%adapt_space == 'IMAhp' .or. state%space%adapt%adapt_space == 'IMAh') then
     state%space%adapt%adapt_method = 'Ihp'
     state%space%adapt%adapt_type = 'Ihp'

  elseif(state%space%adapt%adapt_space == 'none' .or. state%space%adapt%adapt_space == '-') then
     state%space%adapt%max_adapt_level = 0

  else
     print*,'Unknown mesh adaptation method, possibilities: '
     print*,'          HGh/HGp/HGhp / RGh/RGp/RGhp / AMAh/AMAhp  IMAhp/IMAh'
     stop
  endif

  if( state%space%adapt%adapt_type == 'Ahp' .or.  state%space%adapt%adapt_type == 'Ihp')  allocate(AMA)

  write(*,'(a32,a6,a15, i5)') &
       ' # Mesh adaptation technique:  ',state%space%adapt%adapt_space,&
       ', max levels = ',state%space%adapt%max_adapt_level
  print*

  state%tri_solA = .false.
  if(ch1 == 'Y' .or. ch1 == 'y' .or. ch1 == 'A' .or. ch1 == 'a') state%tri_solA = .true.

  ! maximal number of time steps for each mesh level
  read(input,*) state%time%maxiter
  write(*,'(a54,i6)') '  # Maximal number of time steps for each mesh level =', state%time%maxiter

  ! nonlinear solver

  allocate(state%nlSolver)
  !read(input,*) state%nlSolver%max_iter, state%nlSolver%tol, state%nlSolver%tol2, &
  !     state%nlSolver%min_update

  read(input,*) chA, state%nlSolver%non_alg_stop, state%nlSolver%tol, state%nlSolver%max_iter,&
       state%nlSolver%min_update

  if(chA == 'pseud') then
     print*,' Pseudo-time stepping'
     state%time%time_method = 'P'    ! pseudo-time stepping implicit

  elseif(chA /= 'newton' .and. chA /= 'Newton') then
     print*,' Nonlinear algebraic solver, only Newton is implemented'
     stop
  else
     if(state%nlSolver%non_alg_stop == 'aRES' .or. state%nlSolver%non_alg_stop == 'rezL2' ) then
        write(*,'(a45,a6,a8, es9.2, a11, i3, a13, i3)') &
             '  # Newton-like solver: stopping criterion = ', &
             state%nlSolver%non_alg_stop,', tol = ', state%nlSolver%tol, &
             ', max_iter=',state%nlSolver%max_iter, ', min_update=', state%nlSolver%min_update
        if(state%nlSolver%non_alg_stop == 'aRES') then
           state%nlSolver%tol2 = state%nlSolver%tol
           state%nlSolver%tol = 1E-15
        else
           state%nlSolver%tol2 = -1.  ! NOT USED
        endif

     else
        print*,' Unknown stopping criterion for Newton method,',  &
             ' only "aRES" and "rezL2" are implemented'
        stop
     endif
  endif
  !FR does aRES a tRES work with DWR? it does not right now
   if(state%space%estim_space /= 'DWR') then

      if( (state%nlSolver%non_alg_stop == 'aRES' .or.state%time%estim_time == 'tRES') .and. &
         state%space%estim_space /= 'RES') then
         print*,'TROUBLES in input datas'
         stop
      endif

   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! linear algebraic solver
  read(input,*) ch10, state%linSolver%name, state%linSolver%tol

  state%MGsolver = .false.
  !
  !print*,ch10,len(ch10),ch10=='none',state%linSolver%name, state%linSolver%tol
  !
  if( ch10 == 'UMFPACK' ) then
     state%MGsolver = .true.
     state%linSolver%name = ch10
     write(*,*) ' # direct solver is used -- UMFPACK'
     !
  elseif( ch10 == 'AGMG' ) then
     state%MGsolver = .true.
     state%linSolver%name = ch10
     write(*,*) ' # algebraic multigrid approach -- AGMG by Yvan Notay'
     !
  elseif( ch10 == 'JACOBI' ) then
     state%MGsolver = .true.
     state%linSolver%name = ch10
     write(*,*), ' # iterative solver is used -- block Jacobi'
     !
  elseif( ch10 == 'GS' ) then
     state%MGsolver = .true.
     state%linSolver%name = ch10
     write(*,*), ' # iterative solver is used -- block Gaus-Seidel'
     !
  elseif( ch10 == 'MG1JACJAC' ) then
     state%MGsolver = .true.
     state%linSolver%name = ch10
     write(*,*), ' # linear multigrid approach -- MG with Jacobi as 1x smoother and exact solver'
     !
  elseif( ch10 == 'MG2JACJAC' ) then
     state%MGsolver = .true.
     state%linSolver%name = ch10
     write(*,*), ' # linear multigrid approach -- MG with Jacobi as 2x smoother and exact solver'
     !
  elseif( ch10 == 'MG3JACJAC' ) then
     state%MGsolver = .true.
     state%linSolver%name = ch10
     write(*,*), ' # linear multigrid approach -- MG with Jacobi as 3x smoother and exact solver'
     !
  elseif( ch10 == 'MG1GSGS' ) then
     state%MGsolver = .true.
     state%linSolver%name = ch10
     write(*,*), ' # linear multigrid approach -- MG with G-S as 1x smoother and exact solver'
     !
  elseif( ch10 == 'MG2GSGS' ) then
     state%MGsolver = .true.
     state%linSolver%name = ch10
     write(*,*), ' # linear multigrid approach -- MG with G-S as 2x smoother and exact solver'
     !
  elseif( ch10 == 'MG3GSGS' ) then
     state%MGsolver = .true.
     state%linSolver%name = ch10
     write(*,*), ' # linear multigrid approach -- MG with G-S as 3x smoother and exact solver'
     !
  elseif( ch10 == 'MG1JACGS' ) then
     state%MGsolver = .true.
     state%linSolver%name = ch10
     write(*,*), ' # linear multigrid approach -- with Jacobi as 1x smoother and G-S as exact solver'
     !
  elseif( ch10 == 'MG2JACGS' ) then
     state%MGsolver = .true.
     state%linSolver%name = ch10
     write(*,*), ' # linear multigrid approach -- with Jacobi as 2x smoother and G-S as exact solver'
     !
  elseif( ch10 == 'MG3JACGS' ) then
     state%MGsolver = .true.
     state%linSolver%name = ch10
     write(*,*), ' # linear multigrid approach -- with Jacobi as 3x smoother and G-S as exact solver'
     !
  elseif( ch10 == 'MG_Jacobi1' ) then
     state%MGsolver = .true.
     state%linSolver%name = 'MG_Jacobi1'
     write(*,*) ' # linear multigrid approach -- bJacobi smoother and 1xJacobi iter. as exact solution'
     !
  elseif( ch10 == 'MG_JacobiX' ) then
     state%MGsolver = .true.
     state%linSolver%name = 'MG_JacobiX'
     write(*,*) ' # linear multigrid approach -- bJacobi smoother and 10xJacobi iteration as exact solution'
     !
  elseif( ch10 == 'MGxGMRES' ) then
     state%MGsolver = .true.
     state%linSolver%name = 'MGxGMRES'
     write(*,*) ' # linear multigrid approach -- bJacobi smoother and GMRES as exact solution'
     !
  elseif( ch10 == 'MG_ILU') then
     state%MGsolver = .true.
     state%linSolver%name = 'MG_ILU'
     write(*,*) ' # linear multigrid approach -- ILU smoother'
     !
  else

     if(state%linSolver%name == "GMRES") then
        write(*,'(a54,es10.4)') &
             '  # GMRES linear iterative solver without prec, tol = ',state%linSolver%tol

     elseif(state%linSolver%name == "GMRES_D") then
        write(*,'(a62,es10.4)') &
             '  # GMRES linear iterative solver with block Diag prec, tol = ',state%linSolver%tol

     elseif(state%linSolver%name == "GMRES_ILU") then
        write(*,'(a58,es10.4)') &
             '  # GMRES linear iterative solver with ILU(0) prec, tol = ',state%linSolver%tol

     elseif(state%linSolver%name == "ILU") then
        write(*,'(a58,es10.4)') &
             '  # SMOOTHING ILU(0) prec, tol = ',state%linSolver%tol
     else
        print*,'UNKNOWN linear iterative solver -- STOP'
        STOP
     endif
  endif

  state%linSolver%tol_fixed = .true.
  if (state%linSolver%tol <= 0.D+00) state%linSolver%tol_fixed = .false.

  ! output time and  init number of sol* files
  read(input,*) state%time%OutTime, state%isol
  if(state%time%OutTime > 0. ) &
       write(*,'(a43,es10.4,a15,i3)') &
       '  # Output files with time equidistance  = ', state%time%OutTime, &
       ' numbered from ',state%isol

  ! reading of boundary conditions
  read(input,*) state%numBC
  allocate (state%BC(1:state%numBC) )
  do i=1,state%numBC
     allocate(state%BC(i)%w(1:ndim))
     allocate(state%BC(i)%ww(1:ndim))

     read(input,*) state%BC(i)%ibc, state%BC(i)%inout, state%BC(i)%w(1:ndim)

     if(state%SP) then
        state%BC(i)%ww(1:ndim) = state%BC(i)%w(1:ndim)

     else
        state%BC(i)%ww(1) = state%BC(i)%w(1)

        if(ndim >= 4) then
           state%BC(i)%ww(2:nbDim+1) = state%BC(i)%w(1) * state%BC(i)%w(2:nbDim+1)
           state%BC(i)%ww(nbDim+2) = state%BC(i)%w(nbDim+2)/state%model%kappa1  &
                + dot_product(state%BC(i)%w(2:nbDim+1), state%BC(i)%w(2:nbDim+1)) &
                * state%BC(i)%w(1)/2

           state%BC(i)%press_extrap  = 0.
           if(state%BC(i)%inout .eq. 0) then
              state%rho_infty  = state%BC(i)%w(1)
              state%v_infty  = dot_product(state%BC(i)%w(2:nbDim+1),&
                   state%BC(i)%w(2:nbDim+1))**0.5
              state%p_infty  = state%BC(i)%w(nbDim+2)
              state%alpha_infty = acos (state%BC(i)%w(2)/ state%v_infty)
              state%theta_infty  = state%p_infty/state%rho_infty/state%model%kappa1
           endif

           ! turbulence models
           !if(ndim > 2+nbDim) then
           !   state%BC(i)%ww(nbDim+3:ndim) = state%BC(i)%w(nbDim+3:ndim)
           !endif
           if(state%modelName == 'wet_steam') then
              state%BC(i)%ww(5:8) = state%BC(i)%w(1)*state%BC(i)%w(5:8)   ! wet steam components of state vector: \rho*\omega; \rho*Q2; \rho*Q1; \rho*Q0;
           endif

           !write(*,*) 'EDE w',i , state%BC(i)%w(:)
           !write(*,*) 'EDE ww', state%BC(i)%ww(:)
           !print*,'------------------------------------'
        endif
     endif  ! end of SP
  enddo  ! i=1,state%numBC


  if(ndim >= 4 ) then
     write(*,'(a37,6es9.2)')'  # Far field: rho, v, alpha, p, M, E: ', &
          state%rho_infty, state%v_infty, state%alpha_infty/3.1415926*180, &
          state%p_infty, state%v_infty/(state%model%kappa*state%p_infty/state%rho_infty)**0.5, &
          state%p_infty/state%model%kappa1-0.5*state%rho_infty*state%v_infty**2.

!     write(*,'(a37,6es9.2)')'  # Far field: rho, v, alpha, p, M, E, rho*(omega, Q2, Q1, Q0): ', &
!          state%rho_infty, state%v_infty, state%alpha_infty/3.1415926*180, &
!          state%p_infty, state%v_infty/(state%model%kappa*state%p_infty/state%rho_infty)**0.5, &
!          state%p_infty/state%model%kappa1-0.5*state%rho_infty*state%v_infty**2., &
!          state%BC(1)%ww(5), state%BC(1)%ww(6), state%BC(1)%ww(7), state%BC(1)%ww(8)

  endif

  read(input,*) state%ST_Vp, state%ST_Vc, state%ST_Ep, state%ST_Ec
  write(*,'(a30,4es9.2)') '  # Stabilization parameters : ',  &
       state%ST_Vp, state%ST_Vc, state%ST_Ep, state%ST_Ec


  !user-specified parameters for RTN reconstruction-based adaptive solution algorithm
  read(input,*) gamma_rem, gamma_alg, gamma_lin, nu, stop_crit
  write(*,'(a45, 3es12.4, i4, a2)') &
       '  # gamma_rem, gamma_alg, gamma_lin, nu, stop_crit = ', &
       gamma_rem, gamma_alg, gamma_lin, nu, stop_crit


  write(*,*)' # ---------------------------------------------------------'
  if (input > 5) close(input)

  !call Write_state_settings( )
 ! stop


  !! PREPROCESSING
  ! general association, fixed during all computations and adaptations
  call InitProblem()

  call ReadMesh(gridfile, grid)

  call PlotMesh(grid, 'mesh0')
  if (nbDim == 3) then
     call PlotMesh3D(grid)
  endif
  !print*,'# Mesh plotted'

  call SeekNeighbours(grid)
  !print*,'Neigbours found'



  if (nbDim==2) then
     call ReadCurvedBoundary(prof_file)
     call SeekCurvedBoundary(grid)
     !  call SplineCurvedBoundary( )
     !print*,'Curved boundary found'
  endif

  call ComputeGeometry(grid)
  !print*,'Geometry computed'


   !for DWR - false- primal problem, true - dual
  state%dual = .false.
  call PrepareProblem(grid, rsolfile)

!  if (state%space%estim_space = 'DWR') then
!   call InitDualProblem()
!  endif

!  print*, 'stopping after PrepareProblem'
!   stop
  !call TEST_ILU1()
  !stop

  call cpu_time(t2)
  print*,'# Problem prepared for computation within ',t2-state%start_time, ' s'
  print*,'# ------------------------------------------------------'


  call InitSolveProblem(convfile)
  !print*,'InitSolveProblem -- done'

  if(state%space%m_IPG == 5 ) then   ! conforming FEM
     call SetFileNames(.false., command_name, tri_name, sol_name, exa_name, err_name, &
          est_name,  .false.)

     call OutputDGFEMtri(tri_name)
     if(state%time%OutTime > 0.) call OutputDGFEMsol(sol_name)


     !FEM      call ConformingFEM(convfile)

     if(state%time%OutTime <= 0.) call OutputDGFEMsol(sol_name)

     call cpu_time(t2)
     print*,'# ADGFEM finished after ',t2-state%start_time, ' s'
     stop
  endif

  state%Set_R_s = 'Set_R_s_scalar'

  ! original variant
  !call Compute(convfile, solfile )

  ! NEW variant
  call ComputeAD(convfile, solfile )
  !call ComputeAD_OLD(convfile, solfile )

  call cpu_time(t2)

  print*,'# ADGFEM finished after ',t2-state%start_time, ' s'
  !stop

  call CleanMemory ( )

end program


