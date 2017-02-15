module model_mod

 implicit none

 private

 type, public :: Model_t
   integer :: ndim
   real :: Re
   real :: Re1

   real :: kappa, kappa1                    ! problem parameters for NS only
   real :: Pr                           ! problem parameters for NS only

   real :: p0                           ! problem parameters for Pedestrian (pedes) only
   real :: rho_char

   !real:: v_max                       ! given in 'eikonal.f90', module PedestrianSettings
   !real:: rho_max
   !real :: alpha_ped
   real :: tau

   !SCALAR DATA
   integer :: icase                        ! index of the scalar problem
   integer :: idiff                        ! index of diffusion
   integer :: iconv                        ! index of convection
   integer :: ireac                        ! index of reaction
   integer :: iexact                       ! index of the exact solution + RHS
   real :: conv_rat                        ! parameter for the extrapolation of Jacobi
   real :: param1

   logical :: known_sol                    ! is the exact solution a priori known
   logical :: subdomainRHS              ! .true. = rhs is computed differently from the Set_Model_Data, only on a subdomain and possibly multiplied by dPhi/dx(?)
   integer :: rhsTestFunDerivative        ! RHS 0 = f*phi (implicitly), 1 = f * dPhi/dx, 2 = f * dPhi/dy
   logical :: convective                ! .true. = convective terms are presented
   logical :: varying_time_term         ! .true. = a function in front of the time derivative term

   !NS DATA


!FERROR put the following here from state ???
  ! RHS_presented = .true.
  ! homogenDirichlet = .false.

 contains
   procedure :: initModel
  ! procedure :: printModel

   procedure :: init => initModel
 end type Model_t

 contains

 subroutine initModel ( this, Re, isca, t2)
    class (Model_t), intent(inout) :: this
    real, intent(in) ::Re
    integer, intent(in), optional :: isca
    real, intent(in), optional :: t2

   print*, 'Wrong call! Model is just an abstract type.'
   stop
   this%Re = Re
   this%convective = .true.
   this%varying_time_term = .false.     ! .true. = a function in front of the time derivative term

 end subroutine initModel

end module model_mod
