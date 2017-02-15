
!> definition of models which are simulated: Euler, Navier-Stokes, convection-diffusion
module modelNS_mod

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

  !> Euler or the Navier-Stokes equations
  type, EXTENDS( Model_t ), public :: NavierStokes_t
     !   real :: kappa, kappa1
     !   real :: Pr
     
     
   contains
     procedure :: init => initNavierStokes
     
  end type NavierStokes_t
  

  !> pedestrian flow 
  type, EXTENDS( Model_t ), public :: Pedestrian_t
     !real :: kappa, kappa1
     !real :: Pr
     !real :: p0
     
   contains
     procedure :: init => initPedestrian
     
  end type Pedestrian_t

   !> pedestrian flow 
  type, EXTENDS( Model_t ), public :: ShallowWater_t
     !real :: kappa, kappa1
     !real :: Pr
     !real :: p0
     
   contains
     procedure :: init => initShallowWater
     
  end type ShallowWater_t

  public :: initNavierStokes
  public :: initPedestrian
  public :: initShallowWater


 contains

 !> initialization of the NS case
 !> isca, t2 not used
  subroutine initNavierStokes( this, Re, isca, t2)
    class (NavierStokes_t), intent(inout) :: this
    real, intent(in) ::Re
    integer, intent(in), optional :: isca
    real, intent(in), optional :: t2

    this%ndim = 4

    this%Re = Re
    this%kappa = 1.4
    this%kappa1 = this%kappa - 1.
    this%Pr = 0.72   ! Prandtl number (almost all gases)
    this%convective = .true.


    if ( Re == 0.) then
        print*,' # Compressible Euler equations'
        this%Re1 = 0.
     elseif  ( Re > 0.) then
        this%Re1 = 1./this%Re
        print*,' # Compressible Navier-Stokes equations, Re=',this%Re

     else
        print*,'# Reynolds number is negative',this%Re,' STABILIZATION !!!'
        !stop
     endif

  end subroutine initNavierStokes


 !> initialization of the pedestrian flow problem
  subroutine initPedestrian( this, Re, isca, t2)
    class (Pedestrian_t), intent(inout) :: this
    real, intent(in) ::Re
    integer, intent(in), optional :: isca
    real, intent(in), optional :: t2

    this%ndim = 3

    this%Re = 0.  ! Re is p0
    this%Re1 = 0.
    !!this%Pr = 0.005  ! used as the minimal allowed density
    !this%Pr = 1E-3  ! used as the minimal allowed density  WORKS for P_0 approximation
    this%Pr = 5E-3  ! used as the minimal allowed density
    !this%Pr = 1E-2  ! used as the minimal allowed density

    
    this%kappa = t2   
    this%kappa1 = this%kappa - 1.

    ! chacateristic density
    this%rho_char = 1.
    !this%rho_char = 0.1

    !this%p0 = Re
    this%p0 = Re * this%rho_char**this%kappa

    !this%v_max = 2       ! given in 'eikonal.f90', module PedestrianSettings
    !this%rho_max = 9
    !this%alpha_ped = 7.5
    this%tau = 0.61

    write(*,'(a34, f8.2,a8,es12.4)')' # Pedestrian flow (SWE): kappa =', this%kappa,',   p0 =', this%p0
    !  elseif  ( Re > 0.) then
    !     this%Re1 = 1./this%Re
    !     print*,' # Compressible Navier-Stokes equations, Re=',this%Re

    !  else
    !     print*,'# Reynolds number is negative',this%Re,' STABILIZATION !!!'
    !     !stop
    !  endif

  end subroutine initPedestrian



 !> initialization of the shallow water problem
 !> isca, t2 not used
  subroutine initShallowWater( this, Re, isca, t2)
    class (ShallowWater_t), intent(inout) :: this
    real, intent(in) ::Re
    integer, intent(in), optional :: isca
    real, intent(in), optional :: t2

    this%ndim = 3

    this%Re = 0.  ! Re is p0
    this%Re1 = 0.
    !!this%Pr = 0.005  ! used as the minimal allowed density
    !this%Pr = 1E-3  ! used as the minimal allowed density  WORKS for P_0 approximation
    this%Pr = 5E-3  ! used as the minimal allowed density
    !this%Pr = 1E-2  ! used as the minimal allowed density

    
    this%kappa = t2   
    this%kappa1 = this%kappa - 1.

    ! chacateristic density
    this%rho_char = 1.
    !this%rho_char = 0.1

    !this%p0 = Re
    this%p0 = Re * this%rho_char**this%kappa

    !this%v_max = 2       ! given in 'eikonal.f90', module PedestrianSettings
    !this%rho_max = 9
    !this%alpha_ped = 7.5
    this%tau = 0.61

    write(*,'(a34, f8.2,a8,es12.4)')' # Pedestrian flow (SWE): kappa =', this%kappa,',   p0 =', this%p0
    !  elseif  ( Re > 0.) then
    !     this%Re1 = 1./this%Re
    !     print*,' # Compressible Navier-Stokes equations, Re=',this%Re

    !  else
    !     print*,'# Reynolds number is negative',this%Re,' STABILIZATION !!!'
    !     !stop
    !  endif

  end subroutine initShallowWater

end module modelNS_mod
