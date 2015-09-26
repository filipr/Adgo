module model_mod

 implicit none

 private

 type, public :: Model_t
   integer :: ndim
   real :: Re
   real :: Re1

   real :: kappa, kappa1                    ! problem parameters for NS only
   real :: Pr                           ! problem parameters for NS only


   !SCALAR DATA
   integer :: icase                        ! index of the scalar problem
   integer :: idiff                        ! index of diffusion
   integer :: iconv                        ! index of convection
   integer :: ireac                        ! index of reaction
   integer :: iexact                       ! index of the exact solution + RHS
   real :: conv_rat                        ! parameter for the extrapolation of Jacobi
   real :: param1

   !NS DATA


!FERROR put the following here from state ???
  ! RHS_presented = .true.
  ! homogenDirichlet = .false.

 contains
   procedure :: initModel
  ! procedure :: printModel

   procedure :: init => initModel
 end type Model_t


! type, EXTENDS( Model_t ), public :: Scalar_t
!  ! integer :: isca
!
!     integer :: icase                        ! index of the scalar problem
!     integer :: idiff                        ! index of diffusion
!     integer :: iconv                        ! index of convection
!     integer :: ireac                        ! index of reaction
!     integer :: iexact                       ! index of the exact solution + RHS
!     real :: conv_rat                        ! parameter for the extrapolation of Jacobi
!     real :: param1
!
! contains
!   procedure :: initScalarCase
!
!   procedure :: init => initScalarCase
!
! end type Scalar_t


! type, EXTENDS( Model_t ), public :: NavierStokes_t
!   real :: kappa, kappa1
!   real :: Pr
!
!
! contains
!   procedure :: init => initNavierStokes
!
! end type NavierStokes_t

! interface initScalarCase
!   module procedure initScalarCase
! end interface
!
! interface
!   subroutine initScalarCase( this, Re, isca, t2 )
!    import :: Scalar_t
!    class (Scalar_t), intent(inout) :: this
!    real, intent(in) :: Re
!    integer, intent(in), optional :: isca
!    real, intent(in), optional :: t2
!   end subroutine initScalarCase
! end interface
!
!!!
!
! interface
!   subroutine initNavierStokes( this, Re, isca, t2)
!      import :: NavierStokes_t
!      class (NavierStokes_t), intent(inout) :: this
!      real, intent(in) ::Re
!      integer, intent(in), optional :: isca
!      real, intent(in), optional :: t2
!   end subroutine initNavierStokes
! end interface

 contains

 subroutine initModel ( this, Re, isca, t2)
    class (Model_t), intent(inout) :: this
    real, intent(in) ::Re
    integer, intent(in), optional :: isca
    real, intent(in), optional :: t2

   print*, 'Wrong call! Model is just an abstract type.'
   stop
   this%Re = Re


 end subroutine initModel



! subroutine initNavierStokes( this, Re, isca, t2)
!    class (NavierStokes_t), intent(inout) :: this
!    real, intent(in) ::Re
!    integer, intent(in), optional :: isca
!    real, intent(in), optional :: t2
!
!    this%Re = Re
!    this%kappa = 1.4
!    this%kappa1 = this%kappa - 1.
!    this%Pr = 0.72   ! Prandtl number (almost all gases)
!
!
!    if ( Re == 0.) then
!        print*,' # Compressible Euler equations'
!        this%Re1 = 0.
!     elseif  ( Re > 0.) then
!        this%Re1 = 1./this%Re
!        print*,' # Compressible Navier-Stokes equations, Re=',this%Re
!
!     else
!        print*,'# Reynolds number is negative',this%Re,' STABILIZATION !!!'
!        !stop
!     endif
!
!  end subroutine initNavierStokes

! subroutine printModel(this)
!   class( Model_t ), intent(in) :: this
!
!   print*, 'Model Re:', this%Re
!   print*, 'Model Re1:', this%Re1
!
!   select type ( this )
!
!   class is (Scalar_t)
!      print*, 'Model isca:'!, this%icase
!
!   class is (NavierStokes_t)
!      print*, 'Model kappa:'!, this%kappa
!
!   end select
!
! end subroutine printModel

!   subroutine initModel ( model, model_kind )
!      class( Model_t ), intent(inout) :: model
!      character( len = 20 ), intent(in) :: model_kind
!
!      select case(model_kind)
!         case ( 'sca')
!            allocate ( Scalar_t::model )
!
!            select type (model)
!               type is (Scalar_t)
!                  call model%initScalar( 1,2,3)
!            end select
!
!         case ( 'NSe' )
!            allocate( NavierStokes_t::model)
!
!         case default
!            stop 'Unknown model in initModel'
!
!      end select
!
!   end subroutine initModel
!

end module model_mod
