module nonlinear_mod

   implicit none
   !> from Newton_type
   !> structure of the nonlinear solver, only Newton-like method implemented
   type, public :: NonlinearSol_t
      character(len=20) :: name
      character(len=20) ::  non_alg_stop  ! type of stopping criterion
      integer :: max_iter        ! maximal number of the Newton iterations
      integer :: iter            ! number of actual Newton iterations
      integer :: Aiter           ! number of actual Newton iterations in one adaptive cycle
      real :: tol, tol2          ! tolerances for the Newton
      integer :: min_update      ! number of minimal update of the Flux matrix Newton method
      integer :: updates         ! number of performed updates of the Flux matrix Newton method
      real :: norm_res
      real :: lambda, lambda1, lambda_old
      real :: theta
      real :: res, res0, res1, res_ref    ! residuum \f$ \| f(w)\|_{\ell^2} \f$
      real, dimension(:), allocatable ::  x, b, b1, rr  ! arrays for computation

      logical :: implicitly                ! .true. = implicit performace !from state

      contains

      procedure :: init => initNonlinear
      procedure :: InitNLSolverSettings

   end type NonlinearSol_t

   type, public, extends ( NonlinearSol_t ) :: Newton_t

      contains

      procedure :: init => initNewton

   end type Newton_t

   contains

   subroutine initNonlinear( this, name, non_alg_stop, tol, max_iter, min_update )
      class( NonlinearSol_t ), intent ( inout ) :: this
      character(len = 20), intent( in ) :: name
      character(len = 20), intent( in ) :: non_alg_stop
      real, intent( in ) :: tol
      integer, intent( in ) :: max_iter
      integer, intent( in ) :: min_update

      print*, 'Only abstract type Nonlinear_t. Should be changed to ABSTRACT!'


   end subroutine initNonlinear

        !> initialization of for the Newton methods, max_iter, tolerance ...
  subroutine InitNLSolverSettings( this, time_method )
    class( NonlinearSol_t ), intent ( inout ) :: this
    character :: time_method

    this%Aiter = 0
    this%updates = 0

    if( time_method == 'E') then       ! explicit time discretization
       this%implicitly = .false.
       this%max_iter = 1

    elseif( time_method == 'I') then   ! fully implicit time discretization
       this%implicitly = .true.

    elseif( time_method == 'S') then   ! semi-implicit time discretization
       this%implicitly = .true.
       this%max_iter = 1

   elseif( time_method == 'P') then   ! implicit method with the pseudo-time stepping
       this%implicitly = .false.
       !this%Newton%max_iter = 30


    endif

  end subroutine InitNLSolverSettings


   subroutine initNewton( this, name, non_alg_stop, tol, max_iter, min_update )
      class( Newton_t ), intent ( inout ) :: this
      character(len = 20), intent( in ) :: name
      character(len = 20), intent( in ) :: non_alg_stop
      real, intent( in ) :: tol
      integer, intent( in ) :: max_iter
      integer, intent( in ) :: min_update

      this%name = name
      this%non_alg_stop = non_alg_stop
      this%tol = tol
      this%max_iter = max_iter
      this%min_update = min_update

       if(this%non_alg_stop == 'aRES' .or. this%non_alg_stop == 'rezL2' ) then
        write(*,'(a45,a6,a8, es9.2, a11, i3, a13, i3)') &
             '  # Newton-like solver: stopping criterion = ', &
             this%non_alg_stop,', tol = ', this%tol, &
             ', max_iter=', this%max_iter, ', min_update=', this%min_update
        if( this%non_alg_stop == 'aRES' ) then
           this%tol2 = this%tol
           this%tol = 1E-15
        else
           this%tol2 = -1.  ! NOT USED
        endif

     else
        print*,' Unknown stopping criterion for Newton method,',  &
             ' only "aRES" and "rezL2" are implemented'
        stop
     endif

   end subroutine initNewton





end module nonlinear_mod
