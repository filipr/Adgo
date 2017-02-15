module scalar_mod

 ! use main_data
  use model_mod
  use paramets
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

 type, EXTENDS( Model_t ), public :: Scalar_t
  ! integer :: isca

!     integer :: icase                        ! index of the scalar problem
!     integer :: idiff                        ! index of diffusion
!     integer :: iconv                        ! index of convection
!     integer :: ireac                        ! index of reaction
!     integer :: iexact                       ! index of the exact solution + RHS
!     real :: conv_rat                        ! parameter for the extrapolation of Jacobi
!     real :: param1

 contains
   !procedure :: initScalarCase

   procedure :: init => initScaCase

 end type Scalar_t


 type, EXTENDS( Model_t ), public :: TwoEqs_t
  ! integer :: isca

!     integer :: icase                        ! index of the scalar problem
!     integer :: idiff                        ! index of diffusion
!     integer :: iconv                        ! index of convection
!     integer :: ireac                        ! index of reaction
!     integer :: iexact                       ! index of the exact solution + RHS
!     real :: conv_rat                        ! parameter for the extrapolation of Jacobi
!     real :: param1

 contains
   !procedure :: initScalarCase

   procedure :: init => initTwoEqs

 end type TwoEqs_t




!  !> scalar equation
!  public:: InitScalarCase
!  public:: Eval_Convective_Coeffs
!  public:: Eval_Diffusion_Coeffs
!  public:: Set_Model_Data
!
!  public:: Set_f_s_scalar
!  public:: Der_f_s_scalar
!  public:: Set_NumFlux_scalar
!
!  public:: Set_A_s_scalar
!  public:: Set_Ppm_scalar
!  public:: Set_K_sk_scalar
!  public:: Set_R_s_scalar
!
!  public:: Exact_Scalar
!  public:: Der_Exact_Scalar
!  public:: RHS_Scalar
!  public:: Der_TimeDer_Exact_Scalar
!  public:: TimeDer_Exact_Scalar
!
!  public:: Exact_Sol

! interface
!   subroutine initScalarCase( this, Re, isca, t2 )
!    import :: Scalar_t
!    class (Scalar_t), intent(inout) :: this
!    real, intent(in) :: Re
!    integer, intent(in), optional :: isca
!    real, intent(in), optional :: t2
!   end subroutine initScalarCase
! end interface


contains




  !> Initialization of the scalar case
  subroutine initScaCase( this, Re, isca, t2 )
    class (Scalar_t), intent(inout) :: this
    real, intent(in) :: Re
    integer, intent(in), optional :: isca
    real, intent(in), optional :: t2

   ! scalar => state%scalar
   this%ndim = 1

   this%subdomainRHS = .false.
   this%rhsTestFunDerivative = 0 ! rhs: \int f*phi dx
   ! is the exact solution a priori known?
   this%known_sol = .false.

   this%Re = Re
   this%convective = .true.

   if( Re > 0.) then
      this%Re1  = 1./ this%Re
   else
      this%Re1 = 0.
   endif

   if ( present( isca ) .and. present( t2 ) ) then

     if(isca <= 0) then
        stop ' isca for scalar problem must be positive !'
     else
        !allocate(state%scalar)
        this%icase = isca
        this%param1 = t2

        write(*,'(a25, es10.4, a10, i3, a8,es11.4 )') &
             '  # Scalar equation: ve= ',this%Re1, &
             ', icase = ',this%icase, &
             ', par = ',this%param1
     endif

    this%ireac = 0

    select case (this%icase)
    case(:0)
       print*,' UNKNOWN TYPE of scalar%icase !!!'
       stop

    case(1)     ! technical test case
       this%idiff = 1
       this%iconv = 0
       this%iexact = 1
       this%conv_rat = 1.0
       this%ireac = 0

       this%known_sol = .true.
    ! (2): moving front [Vlasak, Dolejsi, Hajek], set epsilon in *.ini, O=[-1,1]^2
    case(2)
       this%idiff = 1
       this%iconv = 3
       this%iexact = 4
       this%conv_rat = 2.0

    case(3)     ! boundary layer, [Roos, Dolejsi], set epsilon in *.ini
       this%idiff = 1
       this%iconv = 2
       !!!!this%iconv = 11
       this%iexact = 6
       this%conv_rat = 1.0

    case(4)     ! singul corner, anisotropic diff, Burgers
       this%idiff = 4
       this%iconv = 4
       this%iexact = 3
       this%conv_rat = 2.0
       this%known_sol = .true.
       !this%param1 = given in *.ini file

    case(5)     ! heat conduction, Nicaise, time (exp(-t))
       this%idiff = 1
       this%iconv = 0
       this%iexact = 17
       this%conv_rat = 1.0
    case(6)     ! singularly pertubed problem, [Raalte 2004]
       this%idiff = 1
       this%iconv = 10
       this%iexact = 27
       this%conv_rat = 1.0
    case(7)     ! heat equation with u = exp(2t + x + y)
       this%idiff = 1
       this%iconv = 0
       this%iexact = 19
       this%conv_rat = 1.0
    case(8)     ! heat equation Nicaise 2006  ?, u = exp(-t)*x*y*(1-x)*(1-y)
       this%idiff = 1
       this%iconv = 0
       this%iexact = 17
       this%conv_rat = 1.0
    case(9)     ! test case
       this%idiff = 1
       this%iconv = 0
       this%iexact = 12
       this%conv_rat = 1.0
    case(10)     ! singul corner, laplace, steady
       this%idiff = 1
       this%iconv = 0
       this%iexact = 3
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(11)     ! singul corner, laplace, steady
       this%idiff = 1
       this%iconv = 0
       this%iexact = 3
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(12)     ! singul corner, laplace, steady
       this%idiff = 1
       this%iconv = 0
       this%iexact = 3
       this%conv_rat = 1.0
       !this%param1 = 0.25
       !this%param1 = given in *.ini file
    case(13)     ! ????
       this%idiff = 1
       this%iconv = 4   !4
       this%iexact = 10   ! 6;  10 = circular layer
       this%conv_rat = 2.0
       !this%param1 = given in *.ini file

    case(14)     ! NONLINEAR elliptic [Houston, Sulli, Robson 2007]
       this%idiff = 7
       this%iconv = 0
       this%iexact = 28
       this%conv_rat = 1.0

       this%known_sol = .true.

    case(15)     ! NONLINEAR elliptic [Houston, Sulli, Robson 2007] 2nd
       this%idiff = 8  ! do not changes, rhs is computaed directly
       this%iconv = 0
       this%iexact = 29
       this%conv_rat = 1.0
       !this%param1 = 2.
       !this%param1 = 2. / 3.
       !this%param1 = given in *.ini file


    case(16)     !  John, Knobloch, stabilization necessary
       this%idiff = 1
       this%iconv = 6
       this%iexact = 21
       this%conv_rat = 1.0


    case(17)     !  linear convection-diffusion, one boundary layer
       this%idiff = 1
       this%iconv = 10
       this%iexact = 30
       this%conv_rat = 1.0

    case(18)     ! LINEAR elliptic. L shape domain [Eibner, Melenk 2007]
       this%idiff = 1  ! do not changes, rhs is computaed directly
       this%iconv = 0
       this%iexact = 31
       this%conv_rat = 1.0

    case(19)     ! singul corner, Laplace
       this%idiff = 1
       this%iconv = 0
       this%iexact = 32
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(20)     ! singul line, Laplace
       this%idiff = 1
       this%iconv = 0
       this%iexact = 33
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(21)     ! Laplace problem for anisotropic mesh refinement
       this%idiff = 1
       this%iconv = 0
       this%iexact = 35  !34
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(22)     ! interior layer
       this%idiff = 1
       this%iconv = 0
       this%iexact = 35
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(23)     ! Kacur: degenerate parabolic problem  (Eymard, Hilhorst, Vohralik 2006)
       this%idiff = 5
       this%iconv = 12
       this%iexact = 15
       this%conv_rat = 2.0
       !this%param1 = given in *.ini file

    case(24)     ! Barenblatt, porus media flow, Radu et all 2008
       this%idiff = 6
       this%iconv = 0
       this%iexact = 26
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(25)   !sin
       this%idiff = 1
       this%iconv = 0
       this%iexact = 11
       this%conv_rat = 2.0

    case(26)   !sin
       this%idiff = 1
       this%iconv = 0
       this%iexact = 12
       this%conv_rat = 2.0

    case(27)   !linear convection diffusion, exponential and parabolic BL
       this%idiff = 1
       this%iconv = 10
       this%iexact = 36
       this%conv_rat = 1.0

    case(28)   !poisson problem with two parabolic boundary layers [AMA-00]
       this%idiff = 1
       this%iconv = 0
       this%iexact = 37
       this%conv_rat = 1.0

    case(29)   !convection dominated flow problem [Knopp, Lube,  Rapin CMAME 02]
       this%idiff = 1
       this%iconv = 13
       this%iexact = 38
       this%conv_rat = 1.0

       this%known_sol = .true.

    case(30)   !
       this%idiff = 1
       this%iconv = 14
       this%iexact = 39
       this%conv_rat = 1.0

    case(31)  !test for ALG2, computing of jumps of normal component of flux reconstruction on a face with HG
       this%idiff = 1
       this%iconv = 0
       this%iexact = 40
       this%conv_rat = 2.0

    case(32)  ! laplace for Strakos
       this%idiff = 1
       this%iconv = 0
       this%iexact = 13
       this%conv_rat = 2.0

    case(33)  !test for ALG2, FNC estimator when HG nodes are present, quadratic func.
       this%idiff = 1
       this%iconv = 0
       this%iexact = 41
       this%conv_rat = 2.0

    case(34)  !test for ALG2, FNC estimator when HG nodes are present, linear func.
       this%idiff = 1
       this%iconv = 0
       this%iexact = 42
       this%conv_rat = 2.0

    case(35)     ! singul corner, laplace, steady
       this%idiff = 1
       this%iconv = 0
       this%iexact = 3
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(36)
       ! Laplace on L-shaped domain [Ainsworth 2005, Robust AEE for NFE approx]
       ! DIFFERENT WITH LL-SHAPED DOMAIN
       this%idiff = 1
       this%iconv = 0
       this%iexact = 43
       this%conv_rat = 1.0

    case(37)   !sin   for book
       this%idiff = 1
       this%iconv = 0
       this%iexact = 44
       this%conv_rat = 2.0

    case(38)   !poisson problem with two identical parabolic boundary layers
       this%idiff = 1
       this%iconv = 0
       this%iexact = 45
       this%conv_rat = 1.0

    case(39)  ! exponential boundary layers, [Roos, Dolejsi], set epsilon in *.ini
       this%idiff = 1
       this%iconv = 0
       !!!this%iexact = 6
       this%iexact = 14
       this%conv_rat = 1.0


    case(40)
       ! Laplace on LL-shaped domain [Vohralik ESAIM]
       ! DIFFERENT WITH L-SHAPED DOMAIN !!!!!
       this%idiff = 1
       this%iconv = 0
       this%iexact = 46
       this%conv_rat = 1.0

    case(41)     !  linear convection-diffusion, parabolic and exponential BLs
       this%idiff = 1
       this%iconv = 10
       this%iexact = 47
       this%conv_rat = 1.0

    case(42)     ! test case for STDGM
       this%idiff = 1
       this%iconv = 1
       !this%ireac = 1
       this%iexact = 48 !51 sin in space, 48 - polynomial in space
       this%conv_rat = 1.0

    case(43)
       print*,'### ATTENTION with alg_estim.f90:'

    case(44)     ! linear convection-reaction equation, increase in time
       this%idiff = 0
       this%iconv = 1
       this%ireac = 1
       this%iexact = 18
       this%conv_rat = 1.0

    case(45)     ! moving peak
       this%idiff = 0
       this%iconv = 1
       this%ireac = 1
       this%iexact = 23
       this%conv_rat = 1.0

    case(46) ! STDGM ODE
       this%idiff = 0
       this%iconv = 0
       this%ireac = 1
       this%iexact = 49
       this%conv_rat = 1.0

    case(47) ! linear convection: nonconstant for LevelSet
       this%idiff = 0
       this%iconv = 15
       this%ireac = 0
       this%iexact = 50
       this%conv_rat = 1.0

    case(48)     ! linear convection-reaction equation, decrease in time
       this%idiff = 1
       this%iconv = 1
       this%ireac = 0
       this%iexact = 5
       this%conv_rat = 1.0

    case(49)     ! linear convection-reaction equation, increase in time
       this%idiff = 1
       this%iconv = 1
       this%ireac = 0
       this%iexact = 18
       this%conv_rat = 1.0

    case(50)   !sin   for AM
       this%idiff = 1
       this%iconv = 0
       this%iexact = 12
       this%conv_rat = 2.0

    case(51)     ! non-linear convection-reaction equation, increase in time
       this%idiff = 4
       this%iconv = 4
       this%ireac = 0
       this%iexact = 18
       this%conv_rat = 2.0

    case(52)     ! non-linear convection-reaction equation, increase in time
       this%idiff = 4
       this%iconv = 4
       this%ireac = 0
       this%iexact = 2
       this%conv_rat = 2.0  ! verify !!!

    case(53)   !poisson problem with two identical parabolic boundary layers
       this%idiff = 1
       this%iconv = 0
       this%iexact = 52
       this%conv_rat = 1.0

    case(54)   !boundary line singularity
       this%idiff = 1
       this%iconv = 0
       this%iexact = 53
       this%conv_rat = 1.0

    case(55)   ! wave front
       this%idiff = 1
       this%iconv = 0
       this%iexact = 54
       this%conv_rat = 1.0

    case(56)   ! interior line singularity
       this%idiff = 1
       this%iconv = 0
       this%iexact = 55
       this%conv_rat = 1.0

    case(57)   ! wave front, convection-diffusion, with "outside" convection
       this%idiff = 2 !7  !1
       this%iconv = 16
       this%iexact = 54
       this%conv_rat = 1.0

    case(58)   ! interior line singularity, nonlinear diffusion , in Dolejsi, May, Roskovec,Solin 2016 for ESCO
       this%idiff = 4
       this%iconv = 0
       this%iexact = 55
       this%conv_rat = 1.0
       this%known_sol = .true. !use to compute J(u) in DWR


    case(59)   ! nonlinear hyperbolic problem
       this%idiff = 0
       this%iconv = 3
       this%iexact = 56
       this%conv_rat = 2.0


    case(60)   ! linear hyperbolic problem
       this%idiff = 0
       this%iconv = 1
       this%iexact = 56
       this%conv_rat = 1.0


    case(61)   ! linear hyperbolic problem, rotating peak
       this%idiff = 0
       this%iconv = 13
       this%iexact = 57
       this%conv_rat = 1.0

    case(62)   ! multiple difficulties
       this%idiff = 1
       this%iconv = 0
       this%iexact = 58
       this%conv_rat = 1.0

    case(63)   ! battery
       this%idiff = 9
       this%iconv = 0
       this%iexact = 59
       this%conv_rat = 1.0

    case(64)  ! linear test case for STATIONARY DWR method
       this%idiff = 1
       this%iconv = 0
       this%ireac = 0
!       this%iexact =  60 ! PEAK in (0.25,0.25) !11
!       this%known_sol = .true.
       this%iexact =  62 ! PEAK in (0.25,0.25) Unknown sol
       this%known_sol = .false.
       this%subdomainRHS = .true.
       this%rhsTestFunDerivative = 0

       this%conv_rat = 1.0
    case(65)
       ! ! LL-shaped test case for STATIONARY DWR
       ! DIFFERENT WITH L-SHAPED DOMAIN !!!!!
       this%idiff = 1
       this%iconv = 0
       this%iexact = 60   ! 46 Vohralik , 60 - DWR peak, check the location of the PEAK
       this%known_sol = .true.

!       this%iexact =  62 ! PEAK in (-0.5,0.5) Unknown sol
!       this%known_sol = .false.

       this%conv_rat = 1.0
    case(66)
      ! DWR nonlinear testcase
      ! NONLINEAR elliptic
      this%idiff = 8
      this%iconv = 0
      this%iexact =  60 ! PEAK in (0.25,0.25) !28
      this%conv_rat = 1.0

      this%known_sol = .true.

    case(67)   ! simplified battery
       this%idiff = 10
       this%iconv = 0
       this%iexact = 61
       this%conv_rat = 1.0

       ! test case for dudx DWR
    case(68)
       this%idiff = 1
       this%iconv = 0
       this%ireac = 0
       this%iexact =  63 !
       this%known_sol = .false.
       this%subdomainRHS = .true.
       this%rhsTestFunDerivative = 1 ! f *dphi/dx

       this%conv_rat = 1.0

    case(69)   ! porous media flow,  Forchheimer 2-term law
       this%idiff = 11
       this%iconv = 0
       this%iexact = 64
       this%conv_rat = 1.0

    case(70)  ! linear test case for STATIONARY DWR method, domain cross, -lapl=1, Ainsworth,Rankin 2012 - Guaranteed computable bounds..., Example 2
       this%idiff = 1
       this%iconv = 0
       this%ireac = 0
       this%iexact =  63 ! unknown solution
       this%known_sol = .false.
       this%conv_rat = 1.0

    case(71)   ! porous media flow,  Darcy
       this%idiff = 12
       this%iconv = 0
       this%iexact = 65
       this%conv_rat = 1.0

    case(72)   ! porous media flow,  Forchheimer 2-term law
       this%idiff = 13
       this%iconv = 0
       this%iexact = 65
       this%conv_rat = 1.0

    case(73:)
       print*,' UNKNOWN TYPE of this%isca !!!'
       stop

    end select
   else
      print*, 'Problem: initScalarCase should by always called with parameters isca and t2!'

   end if !present isca

   ! NO convection
   if(this%iconv == 0) this%convective = .false.

  end subroutine initScaCase

 !> Initialization of the scalar case
  !>
  !> (2): moving front [Vlasak, Dolejsi, Hajek], set epsilon in *.ini, O=[-1,1]^2
  subroutine initTwoEqs( this, Re, isca, t2 )
    class (TwoEqs_t), intent(inout) :: this
    real, intent(in) :: Re
    integer, intent(in), optional :: isca
    real, intent(in), optional :: t2


   this%ndim = 2

   this%Re = Re

   if( Re > 0.) then
        this%Re1  = 1./ this%Re
     else
        this%Re1 = 0.
     endif

   if ( present( isca ) .and. present( t2 ) ) then


     if(isca <= 0) then
        stop ' isca for scalar problem must be positive !'
     else
        !allocate(state%scalar)
        this%icase = isca
        this%param1 = t2

        write(*,'(a25, es10.4, a10, i3, a8,es11.4 )') &
             '  # 2 same equations: ve= ',this%Re1, &
             ', icase = ',this%icase, &
             ', par = ',this%param1
     endif

    this%ireac = 0

    select case (this%icase)
    case(:0)
       print*,' UNKNOWN TYPE of scalar%icase !!!'
       stop

    case(1)     ! technical test case
       this%idiff = 1
       this%iconv = 0
       this%iexact = 1
       this%conv_rat = 1.0

    case(2)     ! moving front [Vlasak, Dolejsi, Hajek], set epsilon in *.ini
       this%idiff = 1
       this%iconv = 3
       this%iexact = 4
       this%conv_rat = 2.0

    case(3)     ! boundary layer, [Roos, Dolejsi], set epsilon in *.ini
       this%idiff = 1
       this%iconv = 2
       !!!!this%iconv = 11
       this%iexact = 6
       this%conv_rat = 1.0

    case(4)     ! singul corner, anisotropic diff, Burgers
       this%idiff = 4
       this%iconv = 4
       this%iexact = 3
       this%conv_rat = 2.0
       !this%param1 = given in *.ini file

    case(5)     ! heat conduction, Nicaise, time (exp(-t))
       this%idiff = 1
       this%iconv = 0
       this%iexact = 17
       this%conv_rat = 1.0
    case(6)     ! singularly pertubed problem, [Raalte 2004]
       this%idiff = 1
       this%iconv = 10
       this%iexact = 27
       this%conv_rat = 1.0
    case(7)     ! heat equation with u = exp(2t + x + y)
       this%idiff = 1
       this%iconv = 0
       this%iexact = 19
       this%conv_rat = 1.0
    case(8)     ! heat equation Nicaise 2006  ?, u = exp(-t)*x*y*(1-x)*(1-y)
       this%idiff = 1
       this%iconv = 0
       this%iexact = 17
       this%conv_rat = 1.0
    case(9)     ! test case
       this%idiff = 1
       this%iconv = 0
       this%iexact = 12
       this%conv_rat = 1.0
    case(10)     ! singul corner, laplace, steady
       this%idiff = 1
       this%iconv = 0
       this%iexact = 3
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(11)     ! singul corner, laplace, steady
       this%idiff = 1
       this%iconv = 0
       this%iexact = 3
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(12)     ! singul corner, laplace, steady
       this%idiff = 1
       this%iconv = 0
       this%iexact = 3
       this%conv_rat = 1.0
       !this%param1 = 0.25
       !this%param1 = given in *.ini file
    case(13)     ! ????
       this%idiff = 1
       this%iconv = 4   !4
       this%iexact = 10   ! 6;  10 = circular layer
       this%conv_rat = 2.0
       !this%param1 = given in *.ini file

    case(14)     ! NONLINEAR elliptic [Houston, Sulli, Robson 2007]
       this%idiff = 7
       this%iconv = 0
       this%iexact = 28
       this%conv_rat = 1.0

    case(15)     ! NONLINEAR elliptic [Houston, Sulli, Robson 2007] 2nd
       this%idiff = 8  ! do not changes, rhs is computaed directly
       this%iconv = 0
       this%iexact = 29
       this%conv_rat = 1.0
       !this%param1 = 2.
       !this%param1 = 2. / 3.
       !this%param1 = given in *.ini file


    case(16)     !  John, Knobloch, stabilization necessary
       this%idiff = 1
       this%iconv = 6
       this%iexact = 21
       this%conv_rat = 1.0


    case(17)     !  linear convection-diffusion, one boundary layer
       this%idiff = 1
       this%iconv = 10
       this%iexact = 30
       this%conv_rat = 1.0

    case(18)     ! LINEAR elliptic. L shape domain [Eibner, Melenk 2007]
       this%idiff = 1  ! do not changes, rhs is computaed directly
       this%iconv = 0
       this%iexact = 31
       this%conv_rat = 1.0

    case(19)     ! singul corner, Laplace
       this%idiff = 1
       this%iconv = 0
       this%iexact = 32
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(20)     ! singul line, Laplace
       this%idiff = 1
       this%iconv = 0
       this%iexact = 33
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(21)     ! Laplace problem for anisotropic mesh refinement
       this%idiff = 1
       this%iconv = 0
       this%iexact = 35  !34
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(22)     ! interior layer
       this%idiff = 1
       this%iconv = 0
       this%iexact = 35
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(23)     ! Kacur: degenerate parabolic problem  (Eymard, Hilhorst, Vohralik 2006)
       this%idiff = 5
       this%iconv = 12
       this%iexact = 15
       this%conv_rat = 2.0
       !this%param1 = given in *.ini file

    case(24)     ! Barenblatt, porus media flow, Radu et all 2008
       this%idiff = 6
       this%iconv = 0
       this%iexact = 26
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(25)   !sin
       this%idiff = 1
       this%iconv = 0
       this%iexact = 11
       this%conv_rat = 2.0

    case(26)   !sin
       this%idiff = 1
       this%iconv = 0
       this%iexact = 12
       this%conv_rat = 2.0

    case(27)   !linear convection diffusion, exponential and parabolic BL
       this%idiff = 1
       this%iconv = 10
       this%iexact = 36
       this%conv_rat = 1.0

    case(28)   !poisson problem with two parabolic boundary layers [AMA-00]
       this%idiff = 1
       this%iconv = 0
       this%iexact = 37
       this%conv_rat = 1.0

    case(29)   !convection dominated flow problem [Knopp, Lube,  Rapin CMAME 02]
       this%idiff = 1
       this%iconv = 13
       this%iexact = 38
       this%conv_rat = 1.0

    case(30)   ! [Hall-PhD]
       this%idiff = 1
       this%iconv = 14
       this%iexact = 39
       this%conv_rat = 1.0

    case(31)  !test for ALG2, computing of jumps of normal component of flux reconstruction on a face with HG
       this%idiff = 1
       this%iconv = 0
       this%iexact = 40
       this%conv_rat = 2.0

    case(32)  ! laplace for Strakos
       this%idiff = 1
       this%iconv = 0
       this%iexact = 13
       this%conv_rat = 2.0

    case(33)  !test for ALG2, FNC estimator when HG nodes are present, quadratic func.
       this%idiff = 1
       this%iconv = 0
       this%iexact = 41
       this%conv_rat = 2.0

    case(34)  !test for ALG2, FNC estimator when HG nodes are present, linear func.
       this%idiff = 1
       this%iconv = 0
       this%iexact = 42
       this%conv_rat = 2.0

    case(35)     ! singul corner, laplace, steady
       this%idiff = 1
       this%iconv = 0
       this%iexact = 3
       this%conv_rat = 1.0
       !this%param1 = given in *.ini file

    case(36)
       ! Laplace on L-shaped domain [Ainsworth 2005, Robust AEE for NFE approx]
       ! DIFFERENT WITH LL-SHAPED DOMAIN
       this%idiff = 1
       this%iconv = 0
       this%iexact = 43
       this%conv_rat = 1.0

    case(37)   !sin   for book
       this%idiff = 1
       this%iconv = 0
       this%iexact = 44
       this%conv_rat = 2.0

    case(38)   !poisson problem with two identical parabolic boundary layers
       this%idiff = 1
       this%iconv = 0
       this%iexact = 45
       this%conv_rat = 1.0

    case(39)  ! exponential boundary layers, [Roos, Dolejsi], set epsilon in *.ini
       this%idiff = 1
       this%iconv = 0
       !!!this%iexact = 6
       this%iexact = 14
       this%conv_rat = 1.0


    case(40)
       ! Laplace on LL-shaped domain [Vohralik ESAIM]
       ! DIFFERENT WITH L-SHAPED DOMAIN !!!!!
       this%idiff = 1
       this%iconv = 0
       this%iexact = 46
       this%conv_rat = 1.0

    case(41)     !  linear convection-diffusion, parabolic and exponential BLs
       this%idiff = 1
       this%iconv = 10
       this%iexact = 47
       this%conv_rat = 1.0

    case(42)     ! test case for STDGM
       this%idiff = 0 !1
       this%iconv = 1 !1
       this%ireac = 1 !1
       this%iexact = 48 !51 sin in space, 48 - polynomial in space
       this%conv_rat = 1.0

    case(43)
       print*,'### ATTENTION with alg_estim.f90:'

    case(44)     ! linear convection-reaction equation, increase in time
       this%idiff = 0
       this%iconv = 1
       this%ireac = 1
       this%iexact = 18
       this%conv_rat = 1.0

    case(45)     ! moving peak
       this%idiff = 0
       this%iconv = 1
       this%ireac = 1
       this%iexact = 23
       this%conv_rat = 1.0

    case(46) ! STDGM ODE
       this%idiff = 0
       this%iconv = 0
       this%ireac = 1
       this%iexact = 49
       this%conv_rat = 1.0

    case(47) ! linear convection: nonconstant for LevelSet
       this%idiff = 0
       this%iconv = 15
       this%ireac = 0
       this%iexact = 50
       this%conv_rat = 1.0

    case(48)     ! linear convection-reaction equation, decrease in time
       this%idiff = 1
       this%iconv = 1
       this%ireac = 0
       this%iexact = 5
       this%conv_rat = 1.0

    case(49)     ! linear convection-reaction equation, increase in time
       this%idiff = 1
       this%iconv = 1
       this%ireac = 0
       this%iexact = 18
       this%conv_rat = 1.0

    case(50)   !sin   for AM
       this%idiff = 1
       this%iconv = 0
       this%iexact = 12
       this%conv_rat = 2.0

    case(51)     ! non-linear convection-reaction equation, increase in time
       this%idiff = 4
       this%iconv = 4
       this%ireac = 0
       this%iexact = 18
       this%conv_rat = 1.0

    case(52)     ! non-linear convection-reaction equation, increase in time
       this%idiff = 4
       this%iconv = 4
       this%ireac = 0
       this%iexact = 2
       this%conv_rat = 1.0

    case(53:)
       print*,' UNKNOWN TYPE of this%isca !!!'
       stop

    end select
   else
      print*, 'Problem: initTwoEqs should by always called with parameters isca and t2!'

   end if !present isca

  end subroutine initTwoEqs

end module scalar_mod
