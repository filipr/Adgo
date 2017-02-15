!> definition of model which are simulated: Euler, Navier-Stokes, convection-diffusion
module wetsteam_mod

 ! use main_data
  use model_mod
!  use f_mapping
!  use mesh_oper
!  use define_state
!  use blocks_integ
!  use model3DNS
!  use model2DNS
!  use modelTurb2e
!  use modelFE
!  use modelLaplace

  implicit none

 type, EXTENDS( Model_t ), public :: WetSteam_t
!   real :: kappa, kappa1
!   real :: Pr


 contains
   procedure :: init => initWetSteam

 end type WetSteam_t


 contains
 !> initialization of the NS case
 !> isca, t2 not used
  subroutine initWetSteam( this, Re, isca, t2)
    class (WetSteam_t), intent(inout) :: this
    real, intent(in) ::Re
    integer, intent(in), optional :: isca
    real, intent(in), optional :: t2

    stop 'initWetSteam not implemented yet'
    !??? ndim?
    this%ndim = 8
    this%convective = .true.

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

  end subroutine initWetSteam

end module wetsteam_mod

