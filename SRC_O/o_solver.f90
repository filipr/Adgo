module solver_mod
   use AMAdata
   use data_mod
   use define_state
   use main_data
   use mesh_mod
   use model_mod
   use modelNS_mod
   use scalar_mod



 implicit none

 private

 public :: allocMesh
! public :: initProblem_o
 public :: readMainData


! public :: test_alloc

 contains


 !> read the basic data from .ini file and set state variables
 subroutine readMainData()
   integer:: input
   character(len=128):: input_file
   integer :: isca, i
   real :: t2, Reynolds, val
   character(len=1) :: ch1
!   character(len=3) :: ch3
!!  character(len=5) :: ch5
   character(len=10) :: ch10 !used for linsolver
   character(len=20) :: chA, non_solver, lin_solver
   character(len=20) ::  non_alg_stop  ! type of stopping criterion for nonlinear solver
   integer :: max_iter        ! maximal number of the Newton iterations
   integer :: min_update      ! number of minimal update of the Flux matrix nonlinear method
   real :: tol        ! tolerances for the nonlinear solver
   real :: FinTime
   real :: tol_min, tol_max, Lq

   integer :: curved_deg ! temporal variable for grid%curved_deg
   integer :: tdeg ! temporal for time disc degree

   character(len=20) :: disc_space
   real :: sigma
   integer :: deg, degP
   character(len=20) :: disc_time

   character(len=20) :: adapt_space
   integer :: max_adapt_level


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

   !FERROR
   !state%SP = .false.

   ! init model - PROBLEM

   read(input,*) state%modelName, Reynolds, isca, t2

   call state%readModelData( Reynolds, isca, t2)

   ! time-dependent?, stopping final time
   read(input,*) FinTime
   if( FinTime < 1E+10 .or. FinTime == 0. )  then
      ! time dependent problem
      state%time_dependent = .true.
      write(*,'(a41,es12.4)') '  # Time dependent problem, final time = ', FinTime
   else
      ! stationary problem, we seek steady-state solution
      state%time_dependent = .false.
      FinTime = 1E+30
      write(*,'(a49,es12.4)') &
          ' # Steady-state problem, fictitious final time =', FinTime
   endif

   !FERROR Where should be this?
   ! stopping tolerance for steady-state reziduum and cD, cL and cM
   read(input,*) state%conv_rez, state%EcD_tol, state%EcL_tol, state%EcM_tol

   ! space dimension (d=2,3), grid
   read(input,*) nbDim,  gridfile

   ! degree of approximation of curved parts of boundaries, name of file
   read(input,*) curved_deg,  prof_file

   write(*,'(a4,i1,a8,a30, a7,i1, a20,a25)') &
       ' # ',nbDim,'D mesh: ', gridfile,'   Q_',curved_deg, &
       ' boundary in file: ',prof_file

   ! Initial conditions
   read(input,*) state%type_IC, rsolfile
   call state%printInitialConditions( rsolfile )

   ! name of output file
   read(input,*) solfile
   i = len_trim(solfile)
   convfile = solfile
   solfile(i+1:i+4) = '.sol'
   convfile(i+1:i+5) = '.conv'

   print*,' # Output in files "',solfile(1:i+4),'" and "',convfile(1:i+5),'"'
   print*

   ! Space discretization method
   allocate( Space_t :: state%space )

   if(.not. state%SP) then
     read(input,*) disc_space, sigma, deg
     !FERROR degP -> deg
     state%space%degP = -1
   elseif(state%modelName == 'incNS') then
      read(input,*) state%space%disc_space, state%space%sigma, state%space%deg, state%space%degP
   else
     stop 'UNKNOWN TYPE in o_main.f90 (34)'
   endif

   call state%space%initDGdata( disc_space, sigma, deg )
   !call state%initSpaceDiscMethod()

   !Time discretization method
   read(input,*) disc_time, tdeg
   !FR temporarily
   !state%time%deg = -50
   !state%time%deg = tdeg

   ! allocate state%time (BDF or STDG)
   call state%allocTimeDiscMethod( disc_time, tdeg, FinTime )

!!!!!!!!!!!!!!!! Choice of the time step !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   read(input, *) state%time%tau_choice, val, state%time%estim_time

   call state%time%initTimeStepAdapt( val )

   ! (space) error estimation method
   read(input,*) state%space%estim_space, tol_max, tol_min, Lq

   call state%space%allocEstims( state%space%estim_space, Lq )

!!!!!!!!!!!!!!!! mesh adaptation method !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   read(input,*) adapt_space, max_adapt_level, ch1
   state%tri_solA = .false.
   if(ch1 == 'Y' .or. ch1 == 'y' .or. ch1 == 'A' .or. ch1 == 'a') state%tri_solA = .true.


   call state%space%allocAdaptation( tol_max, tol_min, adapt_space , max_adapt_level, Lq )
   if( state%space%adapt%adapt_type == 'Ahp' .or.  &
       state%space%adapt%adapt_type == 'Ihp') &
         allocate(AMA)

   call state%printSpaceErrorEstims()

  ! allocate mesh - depends on the adapt method
  call allocMesh( adapt_space, ch1, state%space%adapt%adapt_type  )
  grid%curved_deg = curved_deg

  ! maximal number of time steps for each mesh level
  read(input,*) state%time%maxiter
  write(*,'(a54,i6)') '  # Maximal number of time steps for each mesh level =', state%time%maxiter

!!!!!!!!!!!!!!!! NONLINEAR SOLVER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   read(input,*) non_solver, non_alg_stop, tol, max_iter, min_update

   if( non_solver == 'pseud' ) then
     print*,' Pseudo-time stepping'
     state%time%time_method = 'P'    ! pseudo-time stepping implicit
     stop 'Pseudo-time stepping not implement in ADGo'
   endif

   call state%readNonlinearSolData( non_solver, non_alg_stop, tol, max_iter, min_update )

!!!!!!!!!!!!!!!! LINEAR SOLVER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  read(input,*) ch10, lin_solver, tol

  call state%readLinearSolData( lin_solver, tol, ch10 )

!!!!!!!!!!!!!!!! output time and  init number of sol* files !!!!!!!!!!!!!!!!!!!!
  read(input,*) state%time%OutTime, state%isol
  if(state%time%OutTime > 0. ) &
       write(*,'(a43,es10.4,a15,i3)') &
       '  # Output files with time equidistance  = ', state%time%OutTime, &
       ' numbered from ',state%isol

!!!!!!!!!!!!!!!! reading boundary conditions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  read(input,*) state%numBC

  allocate ( state%BC(1:state%numBC) )

  do i=1,state%numBC
     allocate ( state%BC(i)%w(1:ndim) )
     allocate ( state%BC(i)%ww(1:ndim) )

     !print*, 'ndim=' , ndim

     call state%BC(i)%init( input, state%model )

!FERROR - temporarily ?where to put this?
      write(debug,*) 'reading boundary conditions  where to put v_infty, p_infty... (o_solver)?'
     if( ndim  >= 4 .and. state%BC(i)%inout .eq. 0) then
           state%rho_infty  = state%BC(i)%w(1)
           state%v_infty  = dot_product(state%BC(i)%w(2:nbDim+1),&
                state%BC(i)%w(2:nbDim+1))**0.5
           state%p_infty  = state%BC(i)%w(nbDim+2)
           state%alpha_infty = acos (state%BC(i)%w(2)/ state%v_infty)
           state%theta_infty  = state%p_infty/state%rho_infty/state%model%kappa1

           state%BC(:)%press_extrap =  state%p_infty
           state%BC(:)%rho_extrap = state%rho_infty

      endif
   enddo

  if(ndim >= 4 ) then
     write(*,'(a37,6es9.2)')'  # Far field: rho, v, alpha, p, M, E: ', &
          state%rho_infty, state%v_infty, state%alpha_infty/3.1415926*180, &
          state%p_infty, state%v_infty/(state%model%kappa*state%p_infty/state%rho_infty)**0.5, &
          state%p_infty/state%model%kappa1-0.5*state%rho_infty*state%v_infty**2.
  endif

  !FERROR : WHERE TO PUT THE FOLLOWING?

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


 end subroutine readMainData



 subroutine allocMesh( adapt_space, ch1, adapt_type )
   character(len=20), intent( in ) :: adapt_space
   character(len=1), intent ( in ) :: ch1
   character(len=8), intent ( in ) :: adapt_type      !  Ahp, HG, RG, derived from adapt_space


   select case ( adapt_space )

   case( 'HGhp ')
      allocate( MeshHG_t :: grid )
   case( 'HGh' )
      allocate( MeshHG_t :: grid )
   case( 'HGp' )
      allocate( MeshHG_t :: grid )

   case( 'RGhp' )
      allocate( MeshRG_t :: grid )
   case( 'RGh' )
      allocate( MeshRG_t :: grid )
   case( 'RGp' )
      allocate( MeshRG_t :: grid )

   case( 'AMAhp' )
      allocate( MeshAMA_t :: grid )
      select type ( grid )
      type is ( MeshAMA_t )
!         if( adapt_type == 'Ahp' .or. adapt_type == 'Ihp')  allocate(grid%AMA)
      end select
   case( 'AMAh' )
      allocate( MeshAMA_t :: grid )
      select type ( grid )
      type is ( MeshAMA_t )
 !        if( adapt_type == 'Ahp' .or. adapt_type == 'Ihp')  allocate(grid%AMA)
      end select
   case( 'AMAp' )
      allocate( MeshAMA_t :: grid )
      select type ( grid )
      type is ( MeshAMA_t )
  !       if( adapt_type == 'Ahp' .or. adapt_type == 'Ihp')  allocate(grid%AMA)
      end select

   case( 'ANIhp' )
      stop 'not the case 3eset34dss in o_solver.f90'
   case( 'ANIh' )
      stop 'not the case 3eset34dss in o_solver.f90'
   case( 'ANIp' )
      stop 'not the case 3eset34dss in o_solver.f90'

   case( 'IMAhp' )
      allocate( MeshAMA_t :: grid )
   case( 'IMAh' )
      allocate( MeshAMA_t :: grid )
   case( 'IMAp' )
      allocate( MeshAMA_t :: grid )


   case default
      print*, 'No adaptation. '
      allocate( mesh :: grid )

   end select


   write(*,'(a32,a6,a15, i5)') &
       ' # Mesh adaptation technique:  ',adapt_space,&
       ', max levels = ',state%space%adapt%max_adapt_level


 end subroutine allocMesh



! subroutine test_alloc(model)
!  class(Model_t) :: model
!  class(Model_t), allocatable :: what
!
! select type(model)
! type is (Scalar_t)
!   allocate (Scalar_t::what )
!   print*, 'Scalar'
!
! type is (NavierStokes_t)
!   allocate( NavierStokes_t :: what )
!
!   select type (what)
!   type is (NavierStokes_t)
!  !   what = model ! copy sh
!     print*, 'model%kappa:', model%kappa
!     what%kappa = 4545
!     print*, 'what%kappa:', what%kappa
!
!   end select
!   print*, 'N-S'
!!   what%convection = 4
!!   what%kappa = 5
!!   print*, 'KAPPA =', what%kappa
! end select
!
!
! end subroutine test_alloc



end module solver_mod
