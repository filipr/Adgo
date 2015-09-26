module boundary_mod
   use model_mod
   use modelNS_mod
   use paramets
   use wetsteam_mod
   use scalar_mod

   implicit none

   !> boundary conditions
   type, public :: typeBC
     integer :: ibc                             !!!! index of boundary component
     integer :: inout                           ! 0=inlet, 1=outlet
     real, dimension(:), allocatable :: w       ! prescribed constant vector (physical)
     real, dimension(:), allocatable :: ww      ! prescribed constant vector (conservative)
     real :: press_extrap                       ! extrapolated pressure
     real :: rho_extrap

   contains

     procedure :: init => initBoundaryCond


   end type typeBC

   contains

   subroutine initBoundaryCond( this, input, problem)
      class( typeBC ), intent( inout ) :: this
      integer, intent ( in ) :: input
      class( Model_t ), intent( in ) :: problem
!      character(len=20), intent( in ) :: model
!      real, intent ( in ), optional :: kappa1

      read(input,*) this%ibc, this%inout, this%w(1:problem%ndim)

      this%ww(1) = this%w(1)

      select type ( problem )
      type is ( Scalar_t )

      type is (TwoEqs_t)

      type is ( NavierStokes_t )   ! (ndim >= 4)
         this%ww(2:nbDim+1) = this%w(1) * this%w(2:nbDim+1)
         this%ww(nbDim+2) = this%w(nbDim+2)/problem%kappa1  &
             + dot_product(this%w(2:nbDim+1), this%w(2:nbDim+1)) &
             * this%w(1)/2

         this%press_extrap  = 0.

      type is ( WetSteam_t )
         this%ww(2:nbDim+1) = this%w(1) * this%w(2:nbDim+1)
         this%ww(nbDim+2) = this%w(nbDim+2)/problem%kappa1  &
             + dot_product(this%w(2:nbDim+1), this%w(2:nbDim+1)) &
             * this%w(1)/2

         this%press_extrap  = 0.

         this%ww(5:8) = this%w(1)*this%w(5:8)   ! wet steam components of state vector: \rho*\omega; \rho*Q2; \rho*Q1; \rho*Q0;


      class default
         stop 'Unknown model in initBoundaryCond'

      end select
!FERROR when is ndim= 3 ??
!      if ( ndim==3 ) then
!	      this%ww(2:ndim) = this%w(2:ndim)
   end subroutine initBoundaryCond


end module boundary_mod
