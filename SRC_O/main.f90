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
  !> file containg grid in generalised ANGENER format
!  character(len=50) :: gridfile
!  !> file containg output solution, i.e.. basis coefficients
!  character(len=50) :: solfile
!  !> file containg input solution if any
!  character(len=50) :: rsolfile
!  !> file containing history of computation,
!  !>
!  !> e.g., time derivative residuum,
!  !> physical time, size of time step, local CFL number, number of iteration
!  !> of linear algebra problems, coefficients \f$ c_D,\ c_L,\ c_M\f$, etc.
!  character(len=50) :: convfile
!  !> profile file, a dense list of nodes lying of nonpolygonal parts of boundary
!  character(len=50) :: prof_file

  character(len=50) :: command_name, tri_name, sol_name, exa_name, err_name, est_name
  character(len=3) :: ch3
  character(len=5) :: ch5
  !> number of equations
  integer :: i, i0, isca
  !> real time instants (CPU)
  real t2
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
  ! number of space dimensions of the computational domain
  read(input,*) nbDim
  write(*,*) ' # Dimenzion of problem :',nbDim,'D'

  ! grid for computation
  read(input,*) gridfile

  ! degree of approximation of curved parts of boundaries
  read(input,*) grid%curved_deg

  write(*,'(a12,a30,a18,i5)') '# Grid file:',gridfile,', bound degree = ',grid%curved_deg

  read(input,*) prof_file
  write(*,*) '# File with curved bondary:',prof_file

  read(input,*) solfile
  print*,'# Output of the solution in file ',solfile

  ! number of equations
  read(input,*) ndim, isca, t2
  if(ndim < 1) then
     print*,'Number of PDEs has to be positive',ndim
     stop
  endif

  if(nbDim == 2 .and. ndim == 1) then        ! 2D scalar equation
     state%modelName = 'scalar'

  elseif(nbDim == 2 .and. ndim == 4) then
     state%modelName = 'NSe'
  else
      print*,' Model is not implemented'
   end if

   if(ndim >= 4 .and. ndim <= 6)  then
     !  Euler/N.S. equations
     state%model%kappa = 1.4
     state%model%kappa1 = state%model%kappa - 1.
     state%model%Pr = 0.72   ! Prandtl number (almost all gases)
  endif

  ! scalar case
  if(isca > 0) then
     allocate(state%scalar)
     state%model%icase = isca
     state%model%param1 = t2
     call InitScalarCase( )
  endif


  ! default degree of approximation of the solution
  read(input,*) state%space%deg
  if(state%space%deg < 0) then
     print*,'Default approximation degree has to be nonnegative ',state%space%deg
     stop
  endif
  print*,'# Default degree of polynomial approximation = ',state%space%deg

  !state%time%time_method = 'E'    ! explicit
  !state%time%time_method = 'I'    ! fully implicit
  state%time%time_method = 'S'    ! semi-implicit

  ! order of BDF time discretization
  read(input,*) state%time%time_method, state%time%deg, ch3

  if(state%time%time_method == 'E') then     ! explicit
     print*,'# Explicit time discretization'
  elseif(state%time%time_method == 'I' ) then   ! fully implicit
     print*,'# Fully implicit time discretization'
  elseif(state%time%time_method == 'S') then    ! semi-implicit
     print*,'$ Semi-implicit time discretization'
  else
     print*,'Unknown type of time discretization'
     stop
  endif

  state%time%cn = .false.
  state%time%stdgm  = .false.
  state%time%tdg = .false.

  if(ch3 == 'BDF') then

  elseif(ch3 == 'CN') then ! Crank - Nicolson
     state%time%cn = .true.
  elseif (ch3 == 'TDG') then
     state%time%tdg = .true.
  elseif (ch3 == 'STG') then
     state%time%stdgm  = .true.
     state%time%disc_time = 'STDG'
  else
     print*,'Unknown time integration method'
     stop
  endif

  if( state%time%stdgm ) then
     if(state%time%deg < 0) then
        print*,'Degree of time approximation  of ST DG  has to be positive ',state%time%deg
        stop
     endif
  else
     if(state%time%deg <= 0) then
        print*,'Order of BDF time discretization has to be positive ',state%time%deg
        stop
     endif
  endif

  print*,'# Order of time discretization = ',state%time%deg
  if(state%time%tdg) print*,'# Time discontinuous Galerkin'
  if(state%time%stdgm) print*, '# Space-time discontinuous Galerkin'
  if(.not. state%time%tdg .and. .not. state%time%cn .and. .not. state%time%stdgm) print*,'# Backward difference formulae'
  if(state%time%cn) print*,'# Crank - Nicolson'


  allocate(state%nlSolver)
  read(input,*) state%nlSolver%max_iter, state%nlSolver%tol, state%nlSolver%tol2, &
       state%nlSolver%min_update

  write(*,'(a20,i3,a9,2es12.4,a8,i2)')'# Newton: max_iter=', state%nlSolver%max_iter, &
       ',  tols: ',state%nlSolver%tol, state%nlSolver%tol2, 'update=', state%nlSolver%min_update

  if(state%nlSolver%tol < 0.) then
     state%nlSolver%non_alg_stop = 'aRES'
  else
     state%nlSolver%non_alg_stop = 'precL2'
  endif

  !user-specified parameters for RTN reconstruction-based adaptive solution algorithm
  read(input,*) gamma_rem, gamma_alg, gamma_lin, nu, stop_crit
  write(*,'(a50, 3es12.4, i4, a2)') &
       ' # gamma_rem, gamma_alg, gamma_lin, nu, stop_crit = ', &
       gamma_rem, gamma_alg, gamma_lin, nu, stop_crit

  !type of stabilization: for IPG:  NIPG = 1, IIPG = 0, SIPG = -1
  ! penalty parameter
  read(input,*) state%space%m_IPG, state%space%sigma   !!!!!, state%space%pen_deg

  print*,'# '
  if( state%space%m_IPG == 0) then
     print*,'# IIPG stabilization,  sigma =',state%space%sigma
  elseif( state%space%m_IPG == 1) then
     print*,'# NIPG stabilization,  sigma =',state%space%sigma
  elseif( state%space%m_IPG == -1) then
     print*,'# SIPG stabilization,  sigma =',state%space%sigma
  elseif( state%space%m_IPG == 5) then
     print*,'# conforming finite elements'
  else
     print*,'BAD CHOICE OF STABILIZATION'
     print*,'It should be 1, 0, -1'
  endif


  ! type of initial condition, =0 from file, =1 from BC
  read(input,*) state%type_IC
 ! if(state%type_IC == 0)
  read(input,*) rsolfile


  read(input,*) convfile
  write(*,*) '# History of convergence in file ',convfile

  ! maximum number of time steps
  read(input,*) state%time%maxiter

  print*, 'Max number of time steps', state%time%maxiter

  ! maximum number of levels of adaptation
  read(input,*) state%space%adapt%max_adapt_level, state%space%adapt%adapt_method

  print*, 'Adapt method ' , state%space%adapt%adapt_method
  ! stopping tolerance for steady-state reziduum
  read(input,*) state%conv_rez

  read(input,*) state%EcD_tol !
  read(input,*) state%EcL_tol ! stopping criterium for cD, cL and cM
  read(input,*) state%EcM_tol !

  ! stopping final time
  read(input,*) state%time%FinTime
  if(state%time%FinTime < 100. )  then
     state%time_dependent = .true.
     write(*,*) '# Time dependent problem, final time = ', state%time%FinTime
  else
     state%time_dependent = .false.
     write(*,*) '# Steady-state problem, fictitious final time = ', state%time%FinTime
  endif


  read(input,*) state%time%tau_new
  state%time%tau_fixed = .true.
  if (state%time%tau_new <= 0.D+00) state%time%tau_fixed = .false.
  state%time%tau_fixed_size = state%time%tau_new

  read(input,*) state%linSolver%tol, ch5
  state%linSolver%tol_fixed = .true.
  if (state%linSolver%tol <= 0.D+00) state%linSolver%tol_fixed = .false.

  state%MGsolver = .false.
  if(ch5 == 'LINMG') state%MGsolver = .true.

  state%linSolver%name = "GMRES_ILU"   ! compatiblity of old-Adgfem witn new-AAdgfem
  write(*,'(a31,i7,a7,e11.4,a11,e11.4)') &
       ' # Stopping criteria: max_iter=',state%time%maxiter, &
       ', time=', state%time%FinTime,', reziduum=',state%conv_rez

  ! output time
  read(input,*) state%time%OutTime
  if(state%time%OutTime > 100 ) &
       write(*,*) '# Output files with time equidistance  = ', state%time%OutTime

  ! reading if init number of sol* files
  read(input,*) state%isol
  !if(state%type_IC == 1)  state%isol = 0
  write(*,*) '# Initial values of sol* file = ',state%isol


  ! maximal CFL number
  read(input,*) state%time%CFL
  write(*,*) '# Maximum CFL number = ',abs(state%time%CFL)

  ! tolerance for BDF
  read(input,*) state%time%BDFtol
  if(state%time%BDFtol > 0.) then
     write(*,*) '# Tolerance for adaptive BDF = ',state%time%BDFtol
       if(state%space%adapt%adapt_method == 'RES') then
          state%time%estim_time = 'tRES'
       else
          state%time%estim_time = 'loc'
       endif
  else
     write(*,*) '# No adaptive time step'
  endif

  ! setting of the  time step choice type
  if( state%time%tau_fixed ) then
     state%time%tau_choice = 'fixed'  ! fixed time step with tau=state%time%tau_fixed_size
  else
     if( state%time%BDFtol > 0.) then
        state%time%tau_choice = 'adapt'  ! adaptively chosen time step
     else
        if(state%time%CFL < 0.) then
           state%time%CFL = abs(state%time%CFL)
           state%time%tau_choice = 'cfl'    ! fixed CFL number
        else
           state%time%tau_choice = 'exp'    ! a priori given increase of the time step
        endif
     endif
  endif

  read(input,*) state%model%Re   ! Reynolds number

  if(state%model%Re < 0.) then
     print*,'# Reynolds number/viscosity has to be nonnegative',state%model%Re
     print*,' STABILIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!'
     !stop
  elseif(state%model%Re == 0.) then
     print*,'# Inviscid case: Reynolds number/viscosity is zero',state%model%Re
     state%model%Re1 = 0.
  else
     state%model%Re1 = 1./state%model%Re
     print*,'# Reynolds number/viscosity = ',state%model%Re,'/',state%model%Re1
  endif

  ! reading of boundary conditions
  read(input,*) state%numBC
  allocate (state%BC(1:state%numBC) )
  do i=1,state%numBC
     allocate(state%BC(i)%w(1:ndim))
     allocate(state%BC(i)%ww(1:ndim))

     read(input,*) state%BC(i)%ibc, state%BC(i)%inout, state%BC(i)%w(1:ndim)

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
        if(ndim > 2+nbDim) then
           state%BC(i)%ww(nbDim+3:ndim) = state%BC(i)%w(nbDim+3:ndim)
        endif
     endif
  enddo

  if(ndim >= 4 .and. ndim <= nbDim + 4) then
     write(*,'(a38,6es9.2)')' # Far field: rho, v, alpha, p, M, E: ', &
          state%rho_infty, state%v_infty, state%alpha_infty/3.1415926*180, &
          state%p_infty, state%v_infty/(state%model%kappa*state%p_infty/state%rho_infty)**0.5, &
          state%p_infty/state%model%kappa1-0.5*state%rho_infty*state%v_infty**2.
     write(*,*)
  endif

  read(input,*) state%ST_Vp, state%ST_Vc, state%ST_Ep, state%ST_Ec
  write(*,'(a30,4es9.2)') '# Stabilization parameters  :',  &
       state%ST_Vp, state%ST_Vc, state%ST_Ep, state%ST_Ec

  read(input,*) state%space%adapt%tol_max, state%space%adapt%tol_min, state%space%adapt%Lq

  write(*,*)'# ---------------------------------------------------------'
  if (input > 5) close(input)

  !call Write_state_settings( )
  !stop


  !! PREPROCESSING
  ! general association, fixed during all computations and adaptations
  call InitProblem()

  call ReadMesh(gridfile, grid)

  call PlotMesh(grid, 'mesh0')
  if (nbDim == 3) then
    call PlotMesh3D(grid)
  endif
  print*,'# Mesh plotted'

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


  call PrepareProblem(grid, rsolfile)


  !call TEST_ILU1()
  !stop

  call cpu_time(t2)
  print*,'# Computation prepared within ',t2-state%start_time, ' s'
  print*,'# ------------------------------------------------------'

  call InitSolveProblem(convfile)
  !print*,'InitSolveProblem -- done'

  if(state%space%m_IPG == 5 ) then   ! conforming FEM
     call SetFileNames(.false., command_name, tri_name, sol_name, exa_name, err_name, &
          est_name, .false.)

     call OutputDGFEMtri(tri_name)
     if(state%time%OutTime > 0.) call OutputDGFEMsol(sol_name)


     !FEM   call ConformingFEM(convfile)

     if(state%time%OutTime <= 0.) call OutputDGFEMsol(sol_name)

     call cpu_time(t2)
     print*,'# ADGFEM finished after ',t2-state%start_time, ' s'
     stop
  endif

  state%Set_R_s = 'Set_R_s_scalar'

  ! oroginal variant
  call Compute(convfile, solfile )

  ! NEW variant
  !call ComputeAD(convfile, solfile )

  call cpu_time(t2)

  print*,'# ADGFEM finished after ',t2-state%start_time, ' s'
  !stop

end program adgfem


