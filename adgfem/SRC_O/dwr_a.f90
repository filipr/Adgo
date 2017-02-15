module dwr_alg_mod
!  use dual_element_mod
!  use matrix_oper
!  use lapack_oper
!  use main_data
   use elemental_mod
   use target_functional_mod

   implicit none

   type :: DWR_alg_t
      integer :: iter! number of nonlinear Newton steps
      integer :: restart_primal ! set number of iterations of GMRES primal
      integer :: restart_dual ! set number of iterations of GMRES primal
      integer :: iter_lin_primal
      integer :: iter_lin_dual
      real :: C_A    ! how many time should be alg. error smaller than the space estimate 0.1 - 0.001
      real :: C_L    ! safety parameter ~ 0.1
      logical :: linPdone, linDdone, nlDone ! linPrimalProblem is computed exactly enough, similarly for dual and nonlinear criterion
      logical :: optimal ! TRUE - algorithmic - compute only the estimates, which are truly needed , FALSE - compute all estimates - to compare the decrease
      real :: linTol , nlTol ! actuall stopping criteria for linear and nonlinear problem
      character(len=20) :: file_error_name
      character(len=20) :: file_error_new

   contains
      procedure :: init => initDWR_alg
      procedure :: update => updateDWR_alg
      procedure :: writeErrorFile
!      procedure :: updateDWR
!    procedure :: allocDualElems

   end type DWR_alg_t

   contains

    subroutine initDWR_alg( this )
      class(DWR_alg_t), intent(inout) :: this

      integer :: i
      integer :: iFile

      iFile = 47

      ! set here !
      this%file_error_name = "aDWR_errors"
      this%file_error_new = "aDWR_nlErrors"
      this%restart_primal = 25
      this%restart_dual = 25
      this%C_A = 0.05       ! ratio between lin, nonlin and discr estimates
      this%C_L = 0.5        ! safety parameter
      this%optimal = .true.

      this%nlTol = this%C_A
      this%linTol = this%C_A * this%C_A

      this%iter_lin_primal = 0
      this%iter_lin_dual = 0
      this%iter = 0

      write(*,*) '# DWR: algebraic DWR stopping criterion is used. Parameters set in initDWR_alg.'
      !print*, 'restart_primal = ', this%restart_primal
      !print*, 'restart_dual = ', this%restart_dual
      write(*,'(a9,f5.2)') ' # C_A = ' , this%C_A
      ! clear the file_error_name
      open( iFile, file = this%file_error_name, action="write", status="replace" )
      close( iFile )

   end subroutine initDWR_alg

   !> called from DWR%update, after each adaptation
   subroutine updateDWR_alg( this )
      class(DWR_alg_t), intent(inout) :: this

      ! not used ?
      this%linPdone = .false.
      this%linDdone = .false.
      this%nlDone = .false.

      this%iter = 0 ! count
      this%iter_lin_primal = 0 ! count
      this%iter_lin_dual = 0 ! count
      ! after adaptation linTol cannot be connected to estimNL from last mesh
      this%linTol = this%nlTol * this%C_A

!      print*, 'linTol = ', this%linTol
!      print*, 'nlTol = ', this%nlTol

      this%iter_lin_primal = 0
      this%iter_lin_dual = 0
      this%iter = 0

   end subroutine


   ! write one line intp the file for errors from state%estim()
   ! now it is cleaned after every adaptation
   subroutine writeErrorFile( this )
      class( DWR_alg_t ) :: this
!      integer, intent(in) :: iter ! number of the line
!      integer, intent(in) :: iter_lin_primal !
!      integer, intent(in) :: iter_lin_dual ! number of the lin. alg. iterations
      integer :: iFile = 42

      ! ALL NUMBERS NEED to be squarerooted !!!!!!!!!

      open( iFile, file = this%file_error_name , action="write", position="append" )
        write(iFile,'(i6, i6, i6, 8es14.6)') &
!         state%space%adapt%adapt_level, &
         this%iter , this%iter_lin_primal, this%iter_lin_dual, &
         sqrt(state%estim( dwrE, 1)) , sqrt(state%estim( dwr_aver, 1)) ,  &
         sqrt(state%estim( dwrS, 1)) , sqrt(state%estim( dwrS_abs, 1)) , &
         sqrt(state%estim( dwr_dualS, 1)) , sqrt(state%estim( dwr_dualS_abs, 1)) , &
         sqrt(state%estim( dwrA, 1)) , sqrt(state%estim( dwr_dualA, 1))
      close(iFile)

  end subroutine writeErrorFile




end module dwr_alg_mod
