!> definition of models which are simulated: Euler, Navier-Stokes, convection-diffusion
module model_oper
  use main_data
  use f_mapping
  use mesh_oper
  use define_state
  use blocks_integ
!  use model3DNS
  use model2DNS
  use modelPedes
  use modelSWE
  use modelTurb2e
  use modelFE
  use modelLaplace
  use modelPorous

  implicit none

  !> scalar equation
!  public:: InitScalarCase moved to o_scalar.f90
  public:: Eval_Diffusion_Coeffs
  public:: Eval_Convective_Coeffs
  public:: Set_Model_Data
  public :: Set_Battery
  public:: Set_f_s_scalar
  public:: Der_f_s_scalar
  public:: Set_NumFlux_scalar

  public:: Set_A_s_scalar
  public:: Set_Ppm_scalar
  public:: Set_K_sk_scalar
  public:: Set_R_s_scalar

  public:: Exact_Scalar
  public:: Der_Exact_Scalar
  public:: RHS_Scalar
  public:: Der_TimeDer_Exact_Scalar
  public:: TimeDer_Exact_Scalar

  public:: Detect_apriori_known_singularity

contains

  !> evaluation of diffusion coefficients and their derivatives
  !> \f$ K_{s,k}^{i,j},\ s,k=1,2 (\mbox{space dimension}),\ i,j=1,\dots, ndim\f$,
  !> \f$ ider =0 \Rightarrow K(u),\f$ or \f$ ider =0 \Rightarrow K(|\nabla u|),\f$
  !> \f$ ider =1 \Rightarrow \frac{\rm d}{{\rm d} u} K(u) \f$ or
  !> \f$ ider =1 \Rightarrow \frac{\rm d}{{\rm d} |\nabla u|} K(|\nabla u|) \f$
  function Eval_Diffusion_Coeffs(u, Du, s, k, i, j, Re_1, ider, xi)
    real :: Eval_Diffusion_Coeffs
    real, intent(in) :: u            ! solution
    real, dimension(1:nbDim), intent(in) :: Du ! derivative of the solution
    real, dimension(1:nbDim), intent(in) :: xi ! physical coordinate
    integer, intent(in) :: s,k,i,j   ! indexes of coefficients
    real, dimension(1:iRe), intent(in) :: Re_1         ! viscosity
    integer, intent(in) :: ider      ! =0 => K(u), =1 => d K(u)/ d u
    integer :: imod     ! IMOD
    real :: m, uu, rK, rKp, val1, val2, val3
    real :: viscos, compress, permeab, a0, a1

    imod = state%model%idiff
    !imod = 1    ! Laplace, linear diffusion
    !imod = 2    ! linear diffusion with different coeficients
    !imod = 3    ! nonlinear diffusion, atan
    !imod = 4    ! nonlinear diffusion, atan, anisotrop
    !imod = 5    ! Kacur: degenerate parabolic problem (Eymard, Hilhorst, Vohralik 2006)
    !imod = 6    ! Barenblatt, porus media flow, Radu et all 2008
    !imod = 7    ! NONLINEAR elliptic [Houston, Sulli, Robson 2007]
    !imod = 8    ! NONLINEAR elliptic [Houston, Sulli, Robson 2007] second

    select case (imod)
    case(0)   ! no diffusion
          Eval_Diffusion_Coeffs = 0.

    case(1)     ! linear diffusion
       if(ider == 0) then
          ! functions
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = Re_1(1)
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = Re_1(1)

       else
          ! derivatives
          Eval_Diffusion_Coeffs = 0.
       endif


    case(2)     ! linear diffusion with different coefficients
       if(ider == 0) then
          ! functions
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = 2*Re_1(1)
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.5*Re_1(1)
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = Re_1(1)

       else
          ! derivatives
          Eval_Diffusion_Coeffs = 0.
       endif

    case(3)     ! nonlinear diffusion, atan
       if(ider == 0) then
          ! functions
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = Re_1(1) *( 2 + atan(u) )
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = Re_1(1) *( 2 + atan(u) )

       else
          ! derivatives
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = Re_1(1) /(1+ u*u)
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = Re_1(1) /(1+ u*u)


       endif


    case(4)     ! nonlinear diffusion, atan, anisotrop
       if(ider == 0) then
          ! functions
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = Re_1(1) *( 2 + atan(u) )
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.25* Re_1(1) * ( 2 - atan(u) )
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.

          ! FR TODO ???? in articles there is 0.5*Re_1(1) *( 4.0 + atan(u) )
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = 0.5*Re_1(1) *( 2 + atan(u) )

       else
          ! derivatives
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = Re_1(1) /(1+ u*u)
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = -0.25* Re_1(1) / ( 1 + u*u )
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = 0.5*Re_1(1) /(1+ u*u)

       endif


    case(5)     ! Kacur: degenerate parabolic problem (Eymard, Hilhorst, Vohralik 2006)
       if(ider == 0) then
          ! functions
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = Re_1(1) * 2 * u
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = Re_1(1) * 2 * u

       else
          ! derivatives
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = Re_1(1) * 2.
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = Re_1(1) * 2.

       endif

    case(6)  ! Barenblatt, porus media flow, Radu et all 2008
       m = state%model%param1       !parameter from *.ini file
       uu = u
       if(u < 0.) uu = 0.
       if(ider == 0) then
          ! functions
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = m * uu**(m-1.)
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = m * uu**(m-1.)
       else
          ! derivatives
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = m * (m-1) * uu**(m-2.)
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = m * (m-1) * uu**(m-2.)
       endif

    case(7)     !NONLINEAR elliptic [Houston, Sulli, Robson 2007]
       ! ider
       uu = ( Du(1)*Du(1) + Du(2)*Du(2) )**0.5
       if(ider == 0) then
          ! functions
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = Re_1(1) * (2 + 1./(1+ uu) )
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = Re_1(1) * (2 + 1./(1+ uu) )

       else
          ! derivatives
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = -Re_1(1) / (1 + uu )**2
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = -Re_1(1) / (1 + uu )**2

       endif

    case(8)     !NONLINEAR elliptic [Houston, Sulli, Robson 2007] second
       ! ider
       uu = ( Du(1)*Du(1) + Du(2)*Du(2) )**0.5
       if(ider == 0) then
          ! functions
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = Re_1(1) * (1 + exp(-uu*uu) )
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = Re_1(1) * (1 + exp(-uu*uu) )

       else
          ! derivatives
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = Re_1(1) * exp(-uu*uu)*(-2)*uu
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = Re_1(1) * exp(-uu*uu)*(-2)*uu

       endif

    case(9)     ! battery

       call Set_Battery(xi(1),xi(2), val1, val2, val3)

       if(ider == 0) then
          ! functions
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = val1
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = val2

       else
          ! derivatives
          Eval_Diffusion_Coeffs = 0.
       endif


    case(10)     ! battery

       call Set_Battery_Simplified(xi(1),xi(2), val1, val2, val3)

       if(ider == 0) then
          ! functions
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = val1
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = val2

       else
          ! derivatives
          Eval_Diffusion_Coeffs = 0.
       endif

    case(11)     !porous media flow,  Forchheimer 2-term law
       ! ider
       uu = sqrt( Du(1)*Du(1) + Du(2)*Du(2) )

       rK = Re_1(1) * 2./ ( 1+ sqrt(1 +4 * uu))

       rKp = - Re_1(1) * 4./ ( 1+ sqrt(1 +4 * uu))**2  / sqrt(1+ 4. * uu)

       if(ider == 0) then
          ! functions
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = rK
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = rK

       else
          ! derivatives
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = rKp
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = rKp

       endif

    case(12)     ! porous media flow, case A1

       !call  Set_porous_media_A1(xi(1),xi(2), val1, val2, val3)
       call  Set_porous_media_A1(xi(2), -xi(1), val1, val2, val3)

       ! val1 = permeability
       ! coef = permeability / viscosity / compressibility
       val1 = val1 / 1.3E-03  / 5E-10

       if(ider == 0) then
          ! functions
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = val1 * Re_1(1)
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = val1 * Re_1(1)

          !print*,'#DE#DE#',  val1 * Re_1(1)
       else
          ! derivatives
          Eval_Diffusion_Coeffs = 0.
       endif

    case(13)     !porous media flow,  Forchheimer 2-term law
       ! ider
       !call  Set_porous_media_A2(xi(2), -xi(1), val1, val2, val3)

       !if( val1/Re_1(2) > 1E+3 .or. val1/Re_1(2) < 1E-3) then
       !   write(*,'(a8,8es12.4)') 'diffgfee:', val1, Re_1(2), val1/Re_1(2), xi(1:2)
       !endif

       !permeab = val1
       permeab = Re_1(2)  ! uses the precomputed value

       viscos = 1.3E-03
       compress = 5E-10

       a0 = viscos / permeab
       a1 = 550 /sqrt(permeab )
       !a1 = 0.

       if(( a0 + sqrt(a0*a0 + 4 * a1 * uu)) <=  0.) &
            write(*,'(a8,4es12.4)') 'permeaB:', permeab, a0, a1

       uu = sqrt( Du(1)*Du(1) + Du(2)*Du(2) )

       rK = Re_1(1) * 2./ ( a0 + sqrt(a0*a0 + 4 * a1 * uu))

       rKp = - Re_1(1) * 1./ ( a0 + sqrt(a0*a0 +4 * a1 * uu))**2  / sqrt(a0*a0+ 4. * a1*  uu) * 4 *a1

       write(20, '(6es12.4)') xi(1:2), rK / compress

       if(ider == 0) then
          ! functions
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = rK / compress
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = rK  / compress

          !print*,'#DE#DE#',  rK / compress
       else
          ! derivatives
          if( s==1 .and. k==1)  Eval_Diffusion_Coeffs = rKp  / compress
          if( s==1 .and. k==2)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==1)  Eval_Diffusion_Coeffs = 0.
          if( s==2 .and. k==2)  Eval_Diffusion_Coeffs = rKp  / compress

       endif


    case(14:)
       stop 'UNKNOWN TYPE in Eval_Diffusion_Coeffs'

    end select

  end function Eval_Diffusion_Coeffs



  !> setting of the battery from Mitchel
  subroutine Set_porous_media_A1(xii, yii, val1, val2, f)
    real, intent(in) :: xii, yii
    real, intent(inout) :: val1, val2, f
    real, dimension(:), allocatable :: x, y
    real, dimension(:,:), allocatable :: pq
    real :: xi, yi, qq
    integer :: ityp
    integer :: i
    !real:: eps = 1E-5
    real:: eps = 1E-5


    allocate(pq(1:2, 1:3) )

    ! coefficient pq(:, 1 ) = a0 =   K / mu * (1/kappa)  Darcy
    ! coefficient pq(:, 2 ) = a1 =   K / mu * (1/kappa)  Darcy
    ! K = permeability of the medium = 1E-12 ..  1E-15
    ! mu = viscosity of the fluid   = 1.3E-3  (for water)
    ! 1/kappa = compressibility 5E-10 for water

    pq(1, 1:3) = (/ 1E-12, 1E-12, 0. /)
    pq(2, 1:3) = (/ 1E-15, 1E-15, 0. /)


    allocate(x(1:4), y(1:4) )

    x(1:4) = (/ -2.0 , -0.25,  0.25,  2.5 /)
    y(1:4) = (/ -2.0 , -0.10,  0.10,  2.  /)

    xi = xii
    yi = yii


    ! NEW variant with boundary smoothing
    val1 = 0.;    val2 = 0;   f = 0.;

    call Set_Battery_Cells(x(1), x(2), y(1), y(4), eps, xi, yi, qq)
    ityp = 1;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)

    call Set_Battery_Cells(x(3), x(4), y(1), y(4), eps, xi, yi, qq)
    ityp = 1;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)

    call Set_Battery_Cells(x(2), x(3), y(2), y(3), eps, xi, yi, qq)
    ityp = 1;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)

    call Set_Battery_Cells(x(2), x(3), y(1), y(2), eps, xi, yi, qq)
    ityp = 2;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)

    call Set_Battery_Cells(x(2), x(3), y(3), y(4), eps, xi, yi, qq)
    ityp = 2;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)


    !stop "373ud39dj3w"
    !val1 = 1.
    !val2 = 1.
    !f = 0.
    deallocate(x, y, pq)
  end subroutine Set_porous_media_A1


  !> setting of the battery from Mitchel
  subroutine Set_porous_media_A2(xii, yii, val1, val2, f)
    real, intent(in) :: xii, yii
    real, intent(inout) :: val1, val2, f
    real, dimension(:,:), allocatable :: x
    real, dimension(:,:), allocatable :: pq
    real, dimension (:), allocatable :: xi
    real ::  qq1, qq2, qq
    integer :: ityp
    integer :: i
    !real:: eps = 1E-5
    real:: eps = 1E-5


    allocate(pq(1:2, 1:3) )

    ! coefficient pq(:, 1 ) = a0 =   K / mu * (1/kappa)  Darcy
    ! coefficient pq(:, 2 ) = a1 =   K / mu * (1/kappa)  Darcy
    ! K = permeability of the medium = 1E-12 ..  1E-15
    ! mu = viscosity of the fluid   = 1.3E-3  (for water)
    ! 1/kappa = compressibility 5E-10 for water

    pq(1, 1:3) = (/ 1E-12, 1E-12, 0. /)
    pq(2, 1:3) = (/ 1E-15, 1E-15, 0. /)


    allocate(x(1:8, 1:2) )

    x(1, 1:2) = (/  0.35 ,  1.00 /)
    x(2, 1:2) = (/ -0.15 ,  1.00 /)
    x(3, 1:2) = (/ -0.25 ,  0.10 /)
    x(4, 1:2) = (/  0.25 ,  0.10 /)

    x(5, 1:2) = (/ -0.30 , -1.00 /)
    x(6, 1:2) = (/  0.30 , -1.00 /)
    x(7, 1:2) = (/  0.25 , -0.10 /)
    x(8, 1:2) = (/ -0.25 , -0.10 /)

    allocate(xi(1:2) )
    xi(1) = xii
    xi(2) = yii


    ! NEW variant with boundary smoothing
    val1 = 0.;    val2 = 0;   f = 0.;

    qq1 = InsideQuadrilaterall(xi, x(2,:), x(3,:), x(4,:), x(1,:), eps, .false.)

    qq2 = InsideQuadrilaterall(xi, x(6,:), x(7,:), x(8,:), x(5,:), eps, .false.)


     ! two multicomponents
    ityp = 2;  val1 = val1 + qq1*pq(ityp, 1);   val2 = val2 + qq1*pq(ityp, 2); f = f + qq1*pq(ityp, 3)
    ityp = 2;  val1 = val1 + qq2*pq(ityp, 1);   val2 = val2 + qq2*pq(ityp, 2); f = f + qq2*pq(ityp, 3)

    ! the rest
    qq = 1.- qq1 - qq2
    ityp = 1;  val1 = val1 + qq*pq(ityp, 1);  val2 = val2 + qq*pq(ityp, 2); f = f + qq*pq(ityp, 3)


    ! if(qq1 <= 0. .and. qq2 <= 0. ) then
    ! !    write(21,*) -xi(2), xi(1),  qq1, qq2

    !  elseif(qq1 >= 1. .or. qq2 >= 1. ) then

    ! !    write(22,*) -xi(2), xi(1),  qq1, qq2

    !  else
    !     write(*,*) -xi(2), xi(1), qq1,qq2, qq
    !  endif

    deallocate(x,  xi, pq)

  end subroutine Set_porous_media_A2


  !> setting of the battery from Mitchel
  subroutine Set_Battery(xii, yii, val1, val2, f)
    real, intent(in) :: xii, yii
    real, intent(inout) :: val1, val2, f
    real, dimension(:), allocatable :: x, y
    real, dimension(:,:), allocatable :: pq
    real :: xi, yi, qq
    integer :: ityp
    integer :: i
    !real:: eps = 1E-5
    real:: eps = 1E-5



    allocate(pq(1:5, 1:3) )

    pq(1, 1:3) = (/ 25., 25., 0. /)
    pq(2, 1:3) = (/  7., 0.8, 1. /)
    pq(3, 1:3) = (/  5., 0.0001, 1. /)
    !!!pq(3, 1:3) = (/  5., 0.1, 1. /)
    pq(4, 1:3) = (/ 0.2, 0.2, 0. /)
    pq(5, 1:3) = (/ 0.05, 0.05, 0. /)

    allocate(x(1:5), y(1:8) )
    x(1) = 0.
    x(2) = 6.1
    x(3) = 6.5
    x(4) = 8.0
    x(5) = 8.4

    y(1) = 0.
    y(2) = 0.8
    y(3) = 1.6
    y(4) = 3.6
    y(5) = 18.8
    y(6) = 21.2
    y(7) = 23.2
    y(8) = 24.0

    xi = xii
    yi = yii

    ! NEW variant with boundary smoothing
    val1 = 0.;    val2 = 0;   f = 0.;

    call Set_Battery_Cells(x(1) -2*eps , x(2), y(2), y(3), eps, xi, yi, qq)
    ityp = 5;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)

    call Set_Battery_Cells(x(1) -2*eps , x(2), y(3), y(4), eps, xi, yi, qq)
    ityp = 2;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)

    call Set_Battery_Cells(x(1) -2*eps , x(2), y(4), y(5), eps, xi, yi, qq)
    ityp = 3;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)

    call Set_Battery_Cells(x(1) -2*eps , x(2), y(5), y(6), eps, xi, yi, qq)
    ityp = 2;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)



    call Set_Battery_Cells(x(2) , x(3), y(2), y(6), eps, xi, yi, qq)
    ityp = 4;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)



    call Set_Battery_Cells(x(1) -2*eps , x(4), y(1) -2*eps, y(2), eps, xi, yi, qq)
    ityp = 1;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)

    call Set_Battery_Cells(x(1) -2*eps , x(4), y(7), y(8) +2*eps, eps, xi, yi, qq)
    ityp = 1;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)

    call Set_Battery_Cells(x(4), x(5)+2*eps, y(1)-2*eps, y(8) +2*eps, eps, xi, yi, qq)
    ityp = 1;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)


    call Set_Battery_Cells(x(3), x(4), y(2), y(7), eps, xi, yi, qq)
    ityp = 5;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)

    call Set_Battery_Cells(x(1)-2*eps, x(4), y(6), y(7), eps, xi, yi, qq)
    ityp = 5;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)


    !write(101,*) xi, yi, val1, val2, f


    !stop "373ud39dj3w"
    !val1 = 1.
    !val2 = 1.
    !f = 0.
    deallocate(x, y, pq)
  end subroutine Set_Battery


  !> setting of the battery from Mitchel
  subroutine Set_Battery_Simplified(xii, yii, val1, val2, f)
    real, intent(in) :: xii, yii
    real, intent(inout) :: val1, val2, f
    real, dimension(:), allocatable :: x, y
    real, dimension(:,:), allocatable :: pq
    real :: xi, yi, qq
    integer :: ityp
    integer :: i
    !real:: eps = 1E-5
    real:: eps = 1E-5



    allocate(pq(1:2, 1:3) )

    pq(1, 1:3) = (/ 1E+1, 1E+1, 0.1 /)
    pq(2, 1:3) = (/ 1E0, 1E0, 10. /)

    allocate(x(1:3), y(1:4) )
    x(1) = 0.
    x(2) = 0.5
    x(3) = 1.0

    y(1) = 0.
    y(2) = 0.25
    y(3) = 0.75
    y(4) = 1.


    xi = xii
    yi = yii


    ! NEW variant with boundary smoothing
    val1 = 0.;    val2 = 0;   f = 0.;

    call Set_Battery_Cells(x(1) -2*eps , x(2), y(2), y(3), eps, xi, yi, qq)
    ityp = 1;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)

    call Set_Battery_Cells(x(1) -2*eps , x(2), y(1)-2*eps, y(2), eps, xi, yi, qq)
    ityp = 2;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)

    call Set_Battery_Cells(x(1) -2*eps , x(2), y(3), y(4)+2*eps, eps, xi, yi, qq)
    ityp = 2;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)

    call Set_Battery_Cells(x(2), x(3)+2*eps , y(1) -2*eps, y(4)+2*eps, eps, xi, yi, qq)
    ityp = 2;   val1 = val1 +  qq*pq(ityp, 1);   val2 = val2 +  qq*pq(ityp, 2); f = f +  qq*pq(ityp, 3)


    !write(101,*) xi, yi, val1, val2, f


    !stop "373ud39dj3w"
    !val1 = 1.
    !val2 = 1.
    !f = 0.
    deallocate(x, y, pq)
  end subroutine Set_Battery_Simplified

  !> add one cell of the battery given by frame [x1, x2] x [y1, y2]
  !> tolerance eps, xi, yi the particular node
  subroutine Set_Battery_Cells(x1, x2, y1, y2, eps, xi, yi, qq)
    real, intent(in) :: x1, x2, y1, y2, eps, xi, yi
    real, intent(inout) :: qq
    real :: sinx, siny
    integer :: it

    ! x_1 coordinate
    if(xi > x2 + eps .or. xi < x1 - eps) then
       sinx = 0.

    elseif( xi < x2 - eps .and. xi > x1 + eps) then
       sinx = 1.

    elseif( xi >= x2 - eps .and. xi <= x2 + eps) then
       sinx =   (sin( - (xi - x2 )/ (2*eps)*pi  ) + 1 ) / 2

    elseif( xi >= x1 - eps .and. xi <= x1 + eps) then
       sinx =   (sin( (xi - x1 )/ (2*eps)*pi  ) + 1 ) / 2

    else
       write(*,'(a30, 10es12.4)')'Trouble in Set_Battery_Cells X:', x1, x2, y1, y2, eps, xi, yi
    endif


    ! x_2 coordinate
    if(yi > y2 + eps .or. yi < y1 - eps) then
       siny = 0.

    elseif( yi < y2 - eps .and. yi > y1 + eps) then
       siny = 1.

    elseif( yi >= y2 - eps .and. yi <= y2 + eps) then
       siny =   (sin( - (yi - y2 )/ (2*eps)*pi  ) + 1 ) / 2

    elseif( yi >= y1 - eps .and. yi <= y1 + eps) then
       siny =   (sin( (yi - y1 )/ (2*eps)*pi  ) + 1 ) / 2

    else
       write(*,'(a30, 10es12.4)')'Trouble in Set_Battery_Cells Y:', x1, x2, y1, y2, eps, xi, yi
    endif


    qq = sinx * siny

    !if(qq < 0.) &
    !write(*,'(a8,i5, 40es12.4)') 'S B C:',it, x1, x2,  eps, x1-eps, x1+eps, xi, yi, sinx,qq




  end subroutine Set_Battery_Cells

  !> evaluation of convective coefficients and their derivatives
  !> \f$ f_{s}^{i},\ s=1,2 (\mbox{space dimension}),\ i,j=1,\dots, ndim\f$,
  !> \f$ ider =0 \Rightarrow f(u),\f$
  !> \f$ ider =1 \Rightarrow \frac{\rm d}{{\rm d} u} f(u) \f$
  function Eval_Convective_Coeffs(u, s,  i,  ider, x)
    real :: Eval_Convective_Coeffs
    real, intent(in) :: u            ! solution
    integer, intent(in) :: s, i      ! indexes of coefficients
    integer, intent(in) :: ider      ! =0 => f(u), =1 => d f(u)/ d u
    integer :: imod                ! IMOD
    real, dimension(1:nbDim), intent(in) :: x
    real :: r

    imod = state%model%iconv
    !imod = 0    !  NO convection
    !imod = 1    !  linear convection - (u, u)
    !imod = 2    !  linear convection - (-u, -u)
    !imod = 3    !  Burgers equation  - (u^2/2,  u^2/2)
    !imod = 4    !  Burgers equation  - (-u^2/2, -u^2/2
    !imod = 5    !  1D (x) Burgers equation
    !imod = 6    ! John, Knobloch
    !imod = 7    !  linear convection - (u, u) - perturbed
    !imod = 8    !  stronger "Burgers" equation  - (u^3/3,  u^3/3)
    !imod = 9    ! linear convection (0.6, 0.3)  [ Hilhorst, Vohralik ]
    !imod = 10   ! linear convection (1.0, 0.)   [ Raalte]
    !imod = 12   ! 1D Burgers equation
    !imod = 13   ! linear convection (-y, x)   [Knopp, Lube,  Rapin CMAME 02]
    !imod = 14   !  triple layer  [Hall PhD]?
    !imod = 15   ! non-constant  velocity, given by file
    select case (imod)
    case(0)     ! NO convection
       if(ider == 0) then       ! functions
          Eval_Convective_Coeffs = 0.
       else  ! derivatives
          Eval_Convective_Coeffs = 0.
       endif

    case(1)     ! linear convection - positive
       if(ider == 0) then       ! functions
          Eval_Convective_Coeffs = u
       else  ! derivatives
          Eval_Convective_Coeffs = 1.
       endif

    case(2)     ! linear convection - negative
       if(ider == 0) then       ! functions
          Eval_Convective_Coeffs = -u
       else  ! derivatives
          Eval_Convective_Coeffs = -1.
       endif

    case(3)     ! Burgers equation  - positive
       if(ider == 0) then       ! functions
          Eval_Convective_Coeffs =  u * u /2
       else  ! derivatives
          Eval_Convective_Coeffs =  u
       endif

    case(4)     ! Burgers equation  - positive
       if(ider == 0) then       ! functions
          Eval_Convective_Coeffs =  - u * u /2
       else  ! derivatives
          Eval_Convective_Coeffs =  - u
       endif

    case(5)     ! 1D (x) Burgers equation
       if(s == 1) then
          if(ider == 0) then       ! functions
             Eval_Convective_Coeffs =   0.5 * u * u
          else  ! derivatives
             Eval_Convective_Coeffs =  0.5* 2* u
          endif
       else
          Eval_Convective_Coeffs = 0.
       endif

    case(6)     ! John Knobloch
       if(ider == 0) then       ! functions
          if(s == 1) then       ! functions
             Eval_Convective_Coeffs = cos(-pi/3) * u
          else
             Eval_Convective_Coeffs = sin(-pi/3) * u
          endif
       else  ! derivatives
          if(s == 1) then       ! functions
             Eval_Convective_Coeffs = cos(-pi/3)
          else
             Eval_Convective_Coeffs = sin(-pi/3)
          endif

       endif


    case(7)     ! perturbed simple advection
       if(ider == 0) then       ! functions
          if(s == 1) then       ! functions
             Eval_Convective_Coeffs = 1.1 * u
          else
             Eval_Convective_Coeffs = 1. * u
          endif
       else  ! derivatives
          if(s == 1) then       ! functions
             Eval_Convective_Coeffs = 1.1
          else
             Eval_Convective_Coeffs = 1.
          endif

       endif

    case(8)     ! stronger "Burgers" equation
       if(ider == 0) then       ! functions
          Eval_Convective_Coeffs =  u * u * u / 3
       else  ! derivatives
          Eval_Convective_Coeffs =  u * u
       endif


    case(9)    ! linear convection (0.6, 0.3)  [ Hilhorst, Vohralik ]
       if(ider == 0) then       ! functions
          if(s == 1) then       ! functions
             Eval_Convective_Coeffs = 0.6 * u
          else
             Eval_Convective_Coeffs = 0.3 * u
          endif
       else  ! derivatives
          if(s == 1) then       ! functions
             Eval_Convective_Coeffs = 0.6
          else
             Eval_Convective_Coeffs = 0.3
          endif

       endif

    case(10)    ! linear convection (1.0, 0.)  [ Raaltek ]
       if(ider == 0) then       ! functions
          if(s == 1) then       ! functions
             Eval_Convective_Coeffs = u
          else
             Eval_Convective_Coeffs = 0.
          endif
       else  ! derivatives
          if(s == 1) then       ! functions
             Eval_Convective_Coeffs = 1.0
          else
             Eval_Convective_Coeffs = 0.
          endif

       endif

    case(11)     ! perturbed simple advection
       if(ider == 0) then       ! functions
          if(s == 1) then       ! functions
             Eval_Convective_Coeffs = -u
          else
             Eval_Convective_Coeffs = -u   * 0.1
          endif
       else  ! derivatives
          if(s == 1) then       ! functions
             Eval_Convective_Coeffs = -1.
          else
             Eval_Convective_Coeffs = -1. * 0.1
          endif

       endif

    case(12)     ! 1D Burgers equation
       if(s == 1) then       ! first component
          if(ider == 0) then       ! functions
             Eval_Convective_Coeffs =  u * u /2
          else  ! derivatives
             Eval_Convective_Coeffs =  u
          endif
       else
          Eval_Convective_Coeffs = 0.
       endif


    case(13)     ! ! [Knopp, Lube ..]
       if(ider == 0) then       ! functions
          if(s == 1) then       ! first component
             Eval_Convective_Coeffs =  -x(2) * u
          else
             Eval_Convective_Coeffs =  x(1) * u
          endif
       else ! derivatives
          if(s == 1) then       ! first component
             Eval_Convective_Coeffs =  -x(2)
          else
             Eval_Convective_Coeffs =  x(1)
          endif
       endif
    case(14)     !  triple layer   [Hall PhD]?
       if(ider == 0) then       ! functions
          if(s == 1) then       ! first component
             !if(x(1) >= 1.) Eval_Convective_Coeffs =  u
             if(x(1) >=  1.) Eval_Convective_Coeffs =  x(2) * u
             if(x(1) <  1.) Eval_Convective_Coeffs =  x(2) * u
          else
             !if(x(1) >= 1.) Eval_Convective_Coeffs =  u / 10.
             if(x(1) >=  1.) Eval_Convective_Coeffs =  (1- x(1) )**2 * u
             if(x(1) <  1.) Eval_Convective_Coeffs =  (1- x(1) )**2 * u
          endif
       else ! derivatives
          if(s == 1) then       ! first component
             !if(x(1) >= 1.) Eval_Convective_Coeffs =  1.
             if(x(1) >=  1.)  Eval_Convective_Coeffs =  x(2)
             if(x(1) <  1.)  Eval_Convective_Coeffs =  x(2)
          else
             !if(x(1) >= 1.) Eval_Convective_Coeffs =  1./10
             if(x(1) >=  1.)  Eval_Convective_Coeffs =  (1 - x(1) )**2
             if(x(1) <  1.)  Eval_Convective_Coeffs =  (1 - x(1) )**2
          endif
       endif

    case(15)     !  non-constant velocity for LevelSet
       !x(1:2) contains the prescribed velocity field (overwitten in  SetElementsIC_LevelSet())
       if(ider == 0) then       ! functions
          if(s == 1) then       ! first component
             Eval_Convective_Coeffs =  x(1) * u
          else
             Eval_Convective_Coeffs =  x(2) * u
          endif
       else ! derivatives
          if(s == 1) then       ! first component
             Eval_Convective_Coeffs =  x(1)
          else
             Eval_Convective_Coeffs =  x(2)
          endif
       endif
    case(16)     ! ! central velocity
       if(ider == 0) then       ! functions
          if(s == 1) then       ! first component
             Eval_Convective_Coeffs =  -x(2) * u
          else
             Eval_Convective_Coeffs =  x(1) * u
          endif
       else ! derivatives
          if(s == 1) then       ! first component
             Eval_Convective_Coeffs =  -x(2)
          else
             Eval_Convective_Coeffs =  x(1)
          endif
       endif

       case(17:)
          stop 'UNKNOWN Eval_Convective_Coeffs'
    end select

  end function Eval_Convective_Coeffs

  !> evaluation of reaction coefficients
  function Eval_Reaction_Coeffs(u, Du, x, ider)
    real :: Eval_Reaction_Coeffs
    real, intent(in) :: u            ! solution
    real, dimension(1:nbDim), intent(in) :: Du  ! gradient of the solution
    real, dimension(1:nbDim), intent(in) :: x
    integer, intent(in) :: ider  ! ider = 0 => reaction, ider = 1 -> its derivative
    integer :: imod                ! IMOD
    real :: r

    imod = state%model%ireac
    !imod = 0    !  NO reaction
    !imod = 1    !  linear reaction - (u, u)

    select case (imod)
    case(0)     ! NO reaction
       Eval_Reaction_Coeffs = 0.

    case(1)     ! linear reaction
       if(ider == 0) then       ! functions
          Eval_Reaction_Coeffs = state%model%param1 * u
       else
          Eval_Reaction_Coeffs = state%model%param1  ! its derivative
       endif
    case(2:)
       stop 'UNKNOWN Type in Eval_Reaction_Coeffs'
    end select

  end function Eval_Reaction_Coeffs


  !> setting of the problem data, RHS, exact solution, ...
  !> ityp = 1 ==> w is exact solution,
  !> ityp = 2 ==> w is the right-hand-side
  !> ityp = 3 ==> w is the derivative of the right-hand-side with respect to x
  !> ityp = 4 ==> w is the derivative of the right-hand-side with respect to y
  !> ityp = 5 ==> w is the derivative of the time der. of exact solution with respect to x
  !> ityp = 6 ==> w is the derivative of the time der. of exact solution with respect to y
  !> ityp = 7 ==> w is the time derivative of the exact solution
  !> ityp = 8 ==> w is the derivative of the exact solution with respect x
  !> ityp = 9 ==> w is the derivative of the exact solution with respect y
  subroutine Set_Model_Data(x, t, w, ityp)
    real, dimension(1:nbDim), intent(in) :: x
    real, intent(in) :: t
    integer, intent(in) :: ityp
    real, dimension(1:ndim), intent(out) :: w
    real, dimension(:), allocatable :: du, Re_1
    real, dimension(:,:), allocatable :: d2u

    integer:: imod   ! considered model problem
    real :: r, alpha, beta, delta, v, fac
    real :: u, ut  ! solution, 1st and 2nd order derivatives
    real :: ev, evt, evt2, evtt,ev1,ev2,ev3,ev4,ev0, evx, evy, a1, a2,a3,a4,a5,a6,a7,a8
    real :: a0, c1, c2, evxd, evyd, evxd2, evyd2, evtd
    real :: ev1x, ev1y, ev2x, ev2y, ev1xx, ev1yy, ev2xx, ev2xy, ev2yy
    real :: utx, uty, deltaux, deltauy, xx
    real :: conv(1:nbDim), reac, phi, Dphi
    real :: m, r0, r2, x0, y0, v1, v2, rm, rt, rtt, rf, rft,xc, yc
    real :: z, z1, z2, z11, z12, z22, f, ff, f1, f2, f11, f12, f22, uz, uf, uzz, uzf, uff, rhs
    integer :: i,j
    real, dimension(1:nbDim) :: x_peak ! DWR
    real :: p_L, p_w, p_p, p_BL

    allocate( Re_1(1:iRe) )
    Re_1(:) = 1.
    allocate( du(1:nbDim) , d2u(1:nbDim,1:nbDim)  )

    imod = state%model%iexact
!    print*,'###########', imod


    !imod = 1   ! simple test case   IMOD
    !imod = 2   ! singular corner [0,1] x [ 0,1], time convergence
    !imod = 3   ! singular corner [0,1] x [ 0,1], space convergence
    !imod = 4   ! moving front with zero RHS [-1,1][ x [-1,1]
    !imod = 14   ! fixed front with zero RHS [-1,1][ x [-1,1]
    !imod = 5   ! top on [0,1] x [ 0,1], time convergence
    !imod = 6    ! H-G Roos 1 exp(-t)
    !imod = 7    ! H-G Roos 1a exp(t)
    !imod = 8    ! simple passive
    !imod = 9     ! periodic BC, sin(2 pi (x+y - 2*t) )
    !imod = 10   ! generalization of Houston, Sulli, circle interior layer
    !imod = 11    ! smooth bounded solution (sin)
    !imod = 12    ! smooth bounded solution (sin)
    !imod = 13    ! laplace for Strakos
    !imod = 15    ! Kacur: degenerate parabolic problem,  (Eymard, Hilhorst, Vohralik 2006)
    !imod = 17    ! Nicaise
    !imod = 18    ! Nicaise modified (u = exp(t) )
    !imod= 19     ! Vohralik  u= exp(t)*exp(x)*exp(y)
    !imod = 21    ! John Knobloch
    !imod = 22    ! simple linear advection
    !imod = 23    ! moving peak
    !imod = 24    ! moving peak [ Hilhorst, Vohralik], linear convection (0.6, 0.3)
    !imod = 25    ! circular arc
    !imod = 26    ! Barenblatt, porus media flow, Radu et all 2008
    !imod = 27    ! [Raalte 2004]
    !imod = 28    !  NONLINEAR elliptic [Houston, Sulli, Robson 2007]
    !imod = 29    !  NONLINEAR elliptic [Houston, Sulli, Robson 2007] second
    !imod = 30    !  linear convection-diffusion, one boundary layer

    !imod = 38   !convection dominated flow problem [Knopp, Lube,  Rapin CMAME 02]
    !imod = 39   ! [Hall -Phd]

    ! imod = 60, ! 62 DWR peak
    ! imod = 63 ! DWR du/dx in a subDomain, used also in Ainsworth,Rankin CROSS domain problem


    select case (imod)
    case(1)

       u = 0.0
       ut = 0.0
       du(1:2) = 0.
       d2u(:,:) = 0.

       ! P_0 solution
       !u  = 1.

       ! P_1 solution
       !u = x(1)
       !du(1) = 1.

       ! P_2 solution
!       u     = x(1)*x(1) + t
!       du(1) =    2*x(1)
!       d2u(1,1) = 2.
!       ut = 1.0

       ! P_3 solution
       u     = x(1)*x(1)*x(1) * t
       du(1) =    3*x(1)*x(1) * t
       d2u(1,1) = 6*x(1) * t
       ut = x(1)*x(1)*x(1)

       !u     = x(1)**3 + 2*x(2)**3
       !du(1) =    3*x(1)**2
       !du(2) =    6*x(2)**2
       !d2u(1,1) = 6*x(1)
       !d2u(2,2) = 12*x(2)

       !u     = (1 - x(1))**3
       !du(1) =    -3*(1- x(1))**2
       !d2u(1,1) = 6*(1 - x(1))

       ! P_4 solution
       !u     =       x(1)**4
       !du(1) =     4*x(1)**3
       !d2u(1,1) = 12*x(1)**2

       ! P_7 solution
       !u     =       x(1)**7
       !du(1) =     7*x(1)**6
       !d2u(1,1) = 42*x(1)**5

!       u     =       x(1)**3  !+ 30 *   x(2)**3
!       du(1) =    3.*x(1)**2
!       !du(2) =    90*x(2)**2
!       d2u(1,1) = 6 * x(1)
!       d2u(2,1) = 0.
!       d2u(1,2) = 0.
       !d2u(2,2) = 180 * x(2)

!       u =  t*( 5.0*x(1)**2.0 + 10.0*x(2) )
!       ut = ( 5.0*x(1)**2.0 + 10.0*x(2) )
!       du(1) = 5.0*2.0*x(1)*t
!       du(2) = 10.0*t
!       d2u(:,:) = 0.
!       d2u(1,1) = 10.0 * t

!        u = ( 1 + t ) * ( exp( ) )
!        ut =  (x(1) + 2*x(2))
!        du(1) = 1. * t
!        du(2) = 2.0 * t
!        d2u(:,:) = 0.

!       u =  t*( 2.0*x(1)**(5.0) + 3.0*x(2) )
!       ut = ( 2.0*x(1)**(5.0) + 3.0*x(2) )
!       du(1) = t*2.0*5.0*x(1)**(4.0)
!       du(2) = t*3.0
!       d2u(:,:) = 0.
!       d2u(1,1) = t*2.0*20.0*x(1)**(3.0)


       !u = t**1.5 + 1 + x(1) + x(2)
       !u = t**4.5 + 1 + x(1)**2 + x(2)**3
       !u = t**4.5 + 1 + x(1) + 2*x(2)
       !u = exp(t) + x(1) + 2*x(2)
       !u = exp(t)
       !u = t*t + 1
       !u = t + x(1)**4
       !u = x(1)*x(2)*(1-x(1) ) *(1 - x(2) ) + 0.1  + 0.05*x(1)**4 + 0.1*x(2)**4
       !u =  0.1  + 0.05*x(1)**1 + 0.1*x(2)**1
       !u =  x(1)**4 + x(2)**3.
       !u =  1. - x(1)**2 - x(2)**2.
       !u = t**1.5 + x(1)**4 + x(2)**3.
       !u = x(1)
       !u = 1. + t**2 + 2. * x(1) * x(2)**2.
		 !u = 2.*t + 3. * x(1) * x(2)

       !ut = 2.
       !ut = 1.
       !ut = 1.5 * t**0.5
       !ut = 4.5*t**3.5
       !ut = exp(t)
       !!ut = 2*t
       !ut = x(1)
       !ut = x(1) + x(2)
       !ut = 10* exp(10*t)/exp(10.) *(x(1) + x(2))
       !du(1:2) = exp(10*t)/exp(10.)

       !ut = 10* exp(10*t)/exp(10.)
       !du(1:2) = 0.



       !FRF!d2u(1,1) = ( m*t*x(2) )**2. * (-1.)*sin( m*t*x(1)*x(2) )
       !FRF!d2u(2,2) = ( m*t*x(1) )**2. * (-1.)*sin( m*t*x(1)*x(2) )
       !FRF!d2u(1,2) = (m*t)**2. * x(1)*(2)  * (-1.)*sin(m*t*x(1)*x(2)) + 20.*Pi*t*cos(20.*Pi*t*x(1)*x(2))
       !FRF!d2u(2,1) = d2u(1,2)
       !d2u(1,2) = 3.
       !d2u(2,1) = d2u(1,2)
       !d2u(2,2) = 4. * x(1)

       ! m = 1.0
       !  u = m* x(1)*(1-x(1)) * x(2)*(1-x(2))
       !  ut = 0.

       !  du(1) =m* (1 - 2*x(1)) * x(2)*(1-x(2))
       !  du(2) =m* (1 - 2*x(2)) * x(1)*(1-x(1))

       !  d2u(1,1) = m* (- 2)* x(2)*(1-x(2))
       !  d2u(1,2) = m*(1 - 2*x(1))* (1 - 2*x(2))
       !  d2u(2,1) = d2u(1,2)
       !  d2u(2,2) = m*(- 2)* x(1)*(1-x(1))


    case(2:3)   ! singular corner [0,1] x [ 0,1]
       if(imod == 2) then
          !delta = 10.
          delta = 1.
       else
          !delta = -10.
          delta = -1.
       endif

       !alpha = 7.00
       !alpha = 3.50
       !alpha = 2.00
       !alpha = -0.50
       !alpha = -1.00
       !alpha = -1.50
       !alpha = -2.00
       !alpha = -2.50

       alpha = state%model%param1  ! defined in InitScalarCase( )

       a1 = alpha/2
       a2 = alpha/2-1.
       a3 = alpha/2-2.

       if(imod == 2) then
          !evt = (exp(delta * t) - 1.)/ (exp(delta * 0.5 )  - 1.)     ! unsteady increasing
          !evtt =  delta* exp(delta * t ) / (exp(delta)  - 1.)

          evt = 1. - exp(delta * (t+0.))     ! unsteady descresing
          evtt =  -delta* exp(delta * (t+0.) )

          !if(dot_product(x, x) < 1E-3) write(*,'(a3,5es12.4)') '!!!',t,evt,evtt,delta,&
          !     dot_product(x, x)
       else

          evt = 1.                           ! steady
          evtt = 0.
          if(t == 0.) evt = 0.9

       endif


       r = x(1)*x(1)+x(2)*x(2)

       if( r > 0.) then
          ev0 = 2*x(1)*x(2)*(1. -x(1))*(1.-x(2))
          ev1 = r**a1
          ev2 = r**a2
          ev3 = r**a3

          if(alpha == 0.) then
             ev1 = 1.
             ev2 = 0.
             ev3 = 0.
          endif

          u = evt* ev1 * ev0
          ut = evtt * ev1 * ev0


          du(1) = evt*(a1* ev2 * 2* x(1)* ev0 + 2* ev1 * x(2) *(1-x(2))*(1-2*x(1)))
          du(2) = evt*(a1* ev2 * 2* x(2)* ev0 + 2* ev1 * x(1) *(1-x(1))*(1-2*x(2)))


          d2u(1,1) = evt* (a1*a2*ev3 * 4*x(1)*x(1) * ev0  &
               + 4*a1* ev2 *x(2)*(1-x(2))*(2*x(1)-3*x(1)*x(1)) &
               + 4*a1* ev2 *x(1)*x(2)*(1-x(2))*(1-2*x(1))  - 4*ev1*x(2)*(1-x(2)))


          d2u(2,2) = evt* (a1*a2*ev3 * 4* x(2)*x(2) * ev0  &
               + 4*a1* ev2 *x(1)*(1-x(1))*(2*x(2)-3*x(2)*x(2)) &
               + 4*a1* ev2 *x(1)*x(2)*(1-x(1))*(1-2*x(2))  - 4*ev1*x(1)*(1-x(1)))

          d2u(1,2) = evt*(a1* a2 * ev3* 4*x(1)*x(2)*ev0  &
               + a1 *ev2 *2*x(1)*x(1)*(1-x(1))*(1-2*x(2))*2. &
               + a1*ev2*2*x(2)*x(2)*(1-x(2))*(1-2*x(1))*2. + ev1*(1-2*x(1))*(1-2*x(2))*2.)

          d2u(2,1)  = d2u(1,2)

       else
          if(alpha >= 2.) then
             u = 0.
             ut = 0.
             du(:) = 0.
             d2u(:,:) = 0.
          else
             !print*,'@@@ echo zero R in case (2)(3) in model.f90'
             !print*,'r=',r,', xi=',x(1:2)
             u = 0.
             w(1) = 1E+5
             if(ityp == 2) return
          endif

       endif

    case(4)   ! moving front with zero RHS [-1,1][ x [-1,1]
       if ( abs( (x(1)+ x(2) +1 - t)/2.*state%model%Re ) < 250.0 ) then
         ev1 = exp( (x(1)+ x(2) +1 - t)/2.*state%model%Re )
       else if ( (x(1)+ x(2) +1 - t)/2.*state%model%Re > 0.0 ) then
         ev1 = 1E+90
       else
         ev1 = 1E-90
       endif

       u =  1./(1.+ ev1 )

       du(1:2) = -u*u * ev1 * state%model%Re/2
       ut      = -u*u * ev1 * state%model%Re/2

       if(u <= 1.D-80) u = 0.

       !write(50+state%time%iter,*) x(1),x(2),u, t
       !write(*,'(15es12.4)') x(1),x(2), t, u, du(1:2)


    case(5)   ! top [0,1] x [ 0,1], time dependent
       delta = 1.  ! shiff
       alpha = -10.   ! growth

       ! exponential growth
       !evt = (delta + exp(alpha * t)) !/ (delta + exp(alpha ))
       !evtt = alpha* exp(alpha * t) !/ (delta + exp(alpha ))

       ! linear growth
       alpha = 1.
       evt = delta + alpha * t
       evtt = alpha


       evx = x(1) * (x(1)-1)
       evy = x(2) * (x(2)-1)
       evxd = 2*x(1) - 1
       evyd = 2*x(2) - 1

       u = evt * evx * evy
       ut = evtt * evx * evy

       du(1) = evxd * evy  * evt
       du(2) = evx  * evyd * evt

       d2u(1,1) = 2*evt * evy
       d2u(2,2) = 2*evt * evx

       d2u(1,2) =  evt * evxd * evyd
       d2u(2,1) = d2u(1,2)

       utx = evtt * evxd * evy
       uty = evtt * evx * evyd

       deltaux = 2* evt * evxd   ! ONLY FOR heat equation !!!!!
       deltauy = 2* evt * evyd

       !if(state%time%iter <= 3) write(100,*) state%time%iter, x(1), x(2), u


    case(6)   ! boundary layer, time independent, exp(-t)
       !if(imod == 2) then
       delta = 10.
       c1 = -exp(-state%model%Re)
       c2 = -1 - c1

       !evt = 1- exp(-delta * t)
       !evtt =  delta*exp(-delta * t)

       evt = 1.
       evtt = 0.

       evx = c1 + c2 * (1 - x(1))  + exp((-x(1)) *state%model%Re)
       evy = c1 + c2 * (1 - x(2))  + exp((-x(2)) *state%model%Re)

       evxd = -c2 - state%model%Re * exp((-x(1)) *state%model%Re)
       evyd = -c2 - state%model%Re * exp((-x(2)) *state%model%Re)

       evxd2 = state%model%Re * state%model%Re * exp((-x(1)) *state%model%Re)
       evyd2 = state%model%Re * state%model%Re * exp((-x(2)) *state%model%Re)

       u = evt * evx * evy
       ut = evtt * evx * evy

       du(1) = evxd * evy  * evt
       du(2) = evx  * evyd * evt

       d2u(1,1) = evxd2 * evy   * evt
       d2u(2,2) = evx   * evyd2 * evt

       d2u(1,2) =  evxd * evyd * evt
       d2u(2,1) = d2u(1,2)

    case(7)   ! boundary layer, time dependent, exp(+t)
       !if(imod == 2) then
       delta = 10.
       c1 = -exp(-state%model%Re)
       c2 = -1 - c1

       evt = (exp(delta * t)  - 1) / (exp(delta)  - 1)
       evtt =  delta * exp(delta * t) / (exp(delta)  - 1)

       evt = 1.
       evtt = 0.

       evx = c1 + c2 * (1 - x(1))  + exp((-x(1)) *state%model%Re)
       evy = c1 + c2 * (1 - x(2))  + exp((-x(2)) *state%model%Re)

       evxd = -c2 - state%model%Re * exp((-x(1)) *state%model%Re)
       evyd = -c2 - state%model%Re * exp((-x(2)) *state%model%Re)

       evxd2 = state%model%Re * state%model%Re * exp((-x(1)) *state%model%Re)
       evyd2 = state%model%Re * state%model%Re * exp((-x(2)) *state%model%Re)

       u = evt * evx * evy
       ut = evtt * evx * evy

       du(1) = evxd * evy  * evt
       du(2) = evx  * evyd * evt

       d2u(1,1) = evxd2 * evy   * evt
       d2u(2,2) = evx   * evyd2 * evt


       d2u(1,2) =  evxd * evyd * evt
       d2u(2,1) = d2u(1,2)

       utx = evtt * evxd * evy
       uty = evtt * evx * evyd

       deltaux = - state%model%Re * d2u(1,1) + evxd * evyd2 * evt
       deltauy = evxd2 * evyd   * evt - state%model%Re * d2u(2,2)


    case(8)   ! simple passing of the vertex
       u = exp(-200*((x(1)-0.25)**2 + (x(2) - 0.25)**2) )

    case(9)  ! periodic BC, sin(2 pi (x+y - 2*t) )

       !u = sin(2* pi*  (x(1)+ x(2)  - 2*t) )

       u = (sin(2* pi *(x(1) - x(2)) ) + 5)**0.5

    case(10)   ! generalization of Houston, Sulli, circular layer [0,1] x [ 0,1]
       r0 = 0.5

       !alpha = 7.00
       alpha = 3.00
       !alpha = 2.00
       !alpha = 0.50
       !alpha = -1.50
       !alpha = -2.00

       a1 = alpha
       a2 = alpha-1.
       a3 = alpha-2.

       evt = 1.                           ! steady
       evtt = 0.
       if(t == 0.) evt = 0.9

       r2 = x(1)*x(1)+x(2)*x(2)
       r = r2 ** 0.5

       if( r < r0) then
          u = 0.
          ut = 0.
          du(:) = 0.
          d2u(:,:) = 0.
       else
          ev1 = r**a1
          ev2 = r**a2
          ev3 = r**a3

          if(alpha == 0.) then
             ev1 = 1.
             ev2 = 0.
             ev3 = 0.
          endif

          u = evt* (r - r0)**a1
          ut = evtt * ev0


          du(1) = evt* a1* (r - r0)** a2 * (0.5)* r2**(-0.5) * 2* x(1)
          du(1) = evt* a1* (r - r0)** a2 * (0.5)* r2**(-0.5) * 2* x(2)


          d2u(1,1) = evt* ( a1* (r - r0)** a2 * (0.5)* r2**(-0.5) * 2  &
               +  a1* (r - r0)** a2 * (0.5)* (-0.5)* r2**(-1.5) * 4 * x(1) * x(1) &
               + a1*a2*(r-r0)**a3 * (0.25)* r2**(-1.0) * 4* x(1) * x(1) )

          d2u(2,2) = evt* ( a1* (r - r0)** a2 * (0.5)* r2**(-0.5) * 2  &
               +  a1* (r - r0)** a2 * (0.5)* (-0.5)* r2**(-1.5) * 4 * x(2) * x(2) &
               + a1*a2*(r-r0)**a3 * (0.25)* r2**(-1.0) * 4* x(2) * x(2) )


          d2u(1,2) = evt*(  a1* (r - r0)** a2 * (0.5)* (-0.5)* r2**(-1.5)*4 * x(1) * x(2) &
               + a1*a2*(r-r0)**a3 * (0.25)* r2**(-1.0) * 4* x(1) * x(2) )

          d2u(2,1)  = d2u(1,2)

       endif


    case(11)  !smooth bounded solution (sin)

       u = sin(2* pi *(x(1) + x(2)) )
       ut = 0.

       du(1) = 2* pi * cos(2* pi *(x(1) + x(2)) )
       du(2) = 2* pi * cos(2* pi *(x(1) + x(2)) )

       d2u(1,1) = -4 *pi * pi * u
       d2u(2,2) = -4 *pi * pi * u
       d2u(1,2) = -4 *pi * pi * u
       d2u(2,1)  = d2u(1,2)

    case(12)  !smooth bounded solution (sin)

       u = sin( pi *x(1) )  * sin( pi *x(2) )
       ut = 0.

       du(1) =  pi * cos( pi *x(1) )  * sin( pi *x(2) )
       du(2) =  pi * sin( pi *x(1) )  * cos( pi *x(2) )


       d2u(1,1) = - pi * pi * u
       d2u(2,2) = - pi * pi * u
       d2u(1,2) =   pi * pi * cos(2* pi *x(1) )  * cos(2* pi *x(2) )
       d2u(2,1)  = d2u(1,2)

    case(13)  !laplace for strakos
       evx = x(1)* (1-x(1) )
       evy = x(2)* (1-x(2) )

       u = 16* evx * evy
       if(t <=0. ) u = 0.

       ut = 0.

       du(1) = 16* (1- 2*x(1)) * evy
       du(2) = 16 * (1- 2*x(2)) * evx

       d2u(1,1) = -32 * evy
       d2u(2,2) = -32 * evx
       d2u(1,2) = 16 *(1- 2*x(1))*(1- 2*x(2))
       d2u(2,1)  = d2u(1,2)

    case(14)   ! exponential boundary layer for Laplace: state%model%Re -> scalar%param1
       alpha = state%model%param1  ! defined in InitScalarCase( )

       c1 = -exp(-alpha)
       c2 = -1 - c1

       !delta = 10.
       !evt = 1- exp(-delta * t)
       !evtt =  delta*exp(-delta * t)

       evt = 1.
       evtt = 0.

       evx = c1 + c2 * (1 - x(1))  + exp((-x(1)) *alpha)
       evy = c1 + c2 * (1 - x(2))  + exp((-x(2)) *alpha)

       evxd = -c2 - alpha * exp((-x(1)) *alpha)
       evyd = -c2 - alpha * exp((-x(2)) *alpha)

       evxd2 = alpha * alpha * exp((-x(1)) *alpha)
       evyd2 = alpha * alpha * exp((-x(2)) *alpha)

       u = evt * evx * evy
       ut = evtt * evx * evy

       du(1) = evxd * evy  * evt
       du(2) = evx  * evyd * evt

       d2u(1,1) = evxd2 * evy   * evt
       d2u(2,2) = evx   * evyd2 * evt

       d2u(1,2) =  evxd * evyd * evt
       d2u(2,1) = d2u(1,2)

       !write(99,*) x(:), u, alpha

    case(16)   ! top [0,1] x [ 0,1], time dependent
       !if(imod == 2) then
       delta = 1.
       c1 = -exp(-state%model%Re)
       c2 = -1 - c1

       evt = 1- exp(-delta * t)
       evtt =  delta*exp(-delta * t)

       evx = c1 + c2*x(1)  + exp((x(1) - 1) *state%model%Re)
       evy = c1 + c2*x(2)  + exp((x(2) - 1) *state%model%Re)

       evxd = c2 + state%model%Re * exp((x(1) - 1) *state%model%Re)
       evyd = c2 + state%model%Re * exp((x(2) - 1) *state%model%Re)

       evxd2 = state%model%Re * state%model%Re * exp((x(1) - 1) *state%model%Re)
       evyd2 = state%model%Re * state%model%Re * exp((x(2) - 1) *state%model%Re)

       u = evt * evx * evy
       ut = evtt * evx * evy

       du(1) = evxd * evy  * evt
       du(2) = evx  * evyd * evt

       d2u(1,1) = evxd2 * evy   * evt
       d2u(2,2) = evx   * evyd2 * evt


       d2u(1,2) =  evxd * evyd * evt
       d2u(2,1) = d2u(1,2)

    case(15)   ! Kacur: degenerate parabolic problem (Eymard, Hilhorst, Vohralik 2006)
       r = 0.25  ! time shifting
       v = 0.5  ! convective speed, if /= 0.5 then CHANGES IN CONVECTIVE TERMS NECESSARY !!!

       if(x(1) <= v*t + r) then
          evt = exp(state%model%Re * v*(x(1) - v * t - r)/2  )

          u = 1. - evt

          ut = state%model%Re * v *v / 2 * evt

          du(1) = -state%model%Re * v / 2 * evt
          du(2) = 0.

          d2u(:,:) = 0.

          d2u(1,1) = -state%model%Re * state%model%Re * v *v / 4 * evt
       else
          u = 0.
          ut = 0.
          du(:) = 0.
          d2u(:,:) = 0.
       endif


    case(17)   ! Nicaise
       evt = exp(-t)
       evtt = - exp(-t)
       evx = x(1) * (x(1)-1)
       evy = x(2) * (x(2)-1)
       evxd = 2*x(1) - 1
       evyd = 2*x(2) - 1

!!!
       !evt = 1.
       !evtt = 0.

       u = evt * evx * evy
       ut = evtt * evx * evy

       du(1) = evxd * evy  * evt
       du(2) = evx  * evyd * evt

       d2u(1,1) = 2*evt * evy
       d2u(2,2) = 2*evt * evx

       d2u(1,2) =  evt * evxd * evyd
       d2u(2,1) = d2u(1,2)

       utx = evtt * evxd * evy
       uty = evtt * evx * evyd

       deltaux = 2* evt * evxd
       deltauy = 2* evt * evyd

       !if(state%time%iter <= 3) write(100,*) state%time%iter, x(1), x(2), u

    case(18)   ! Nicaise
       delta = 0.1  ! shiff
       alpha = 10.   ! growth

       ! exponential growth
       evt = (delta + exp(alpha * t)) !/ (delta + exp(alpha ))
       evtt = alpha* exp(alpha * t) !/ (delta + exp(alpha ))
       !evt = (delta + exp(alpha * t)) / (delta + exp(alpha ))
       !evtt = alpha* exp(alpha * t) / (delta + exp(alpha ))

       ! linear growth
       !evt = delta * t
       !evtt = delta

       ! time independent
       !evt = 1.
       !evtt = 0.


       evx = x(1) * (x(1)-1)
       evy = x(2) * (x(2)-1)
       evxd = 2*x(1) - 1
       evyd = 2*x(2) - 1

       u = evt * evx * evy + delta
       ut = evtt * evx * evy

       du(1) = evxd * evy  * evt
       du(2) = evx  * evyd * evt

       d2u(1,1) = 2*evt * evy
       d2u(2,2) = 2*evt * evx

       d2u(1,2) =  evt * evxd * evyd
       d2u(2,1) = d2u(1,2)

       utx = evtt * evxd * evy
       uty = evtt * evx * evyd

       deltaux = 2* evt * evxd   ! ONLY FOR heat equation !!!!!
       deltauy = 2* evt * evyd

       !if(state%time%iter <= 3) write(100,*) state%time%iter, x(1), x(2), u

    case(19)   ! vohralik

       evt = exp(2.*t)
       evtt = 2.*exp(2.*t)

       !evt = 1.  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !evtt = 0. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       evx = exp(x(1))
       evy = exp(x(2))

       evxd = exp(x(1) )
       evyd = exp(x(2) )

       u = evt * evx * evy
       ut = evtt * evx * evy

       du(1) = evxd * evy  * evt
       du(2) = evx  * evyd * evt

       d2u(1,1) = evt * evy * evx
       d2u(2,2) = evt * evx * evy

       d2u(1,2) =  evt * evxd * evyd
       d2u(2,1) = d2u(1,2)

    case(21)   ! John Knobloch - analytical solution unknown
       if(x(2) < 0.7 - 3**0.5 * x(1) .or. (x(1) >=0.999999) .or. &
            (x(1) > 0.7/(3**0.5) .and. x(2) < 0.0000001) ) then
          u = 0.
       else
          u = 1.
       endif
       !!write(99,'(12es12.4)') x(:), u

       !if(x(1) >= 0.9999999 .or. x(2) <= 0.7) then
       !   u = 0.
       !   !write(*,'(a6,12es12.4)') '!!!!',x(:), u
       !else
       !   u = 1.
       !endif

       ut = 0.
       du(:) = 0.
       d2u(:,:) = 0.


    case(22)   ! simple linear advection
       if(x(2) >= x(1)/ 1.1 ) then
          u = 0.
       else
          u = 1.
       endif
       !!write(99,'(12es12.4)') x(:), u

       ut = 0.
       du(:) = 0.
       d2u(:,:) = 0.

    case(23)   ! top [0,1] x [ 0,1] moving peak

       delta = -100.

       evx = (x(1) - t - 0.25)
       evy = (x(2) - t - 0.25)

       u = exp( delta * (evx**2  + evy**2) )
       ut = u * delta *(- 2*evx - 2*evy)

       du(1) = u * delta* 2 * evx
       du(2) = u * delta* 2 * evy

       d2u(1,1) = delta*( u * delta * 2* evx * 2 * evx + u * 2 )
       d2u(2,2) = delta*( u * delta * 2* evy * 2 * evy + u * 2 )

       d2u(1,2) = delta* u * delta * 2* evx * 2 * evy
       d2u(2,1) = d2u(1,2)

    case(24)   ! top [0,1] x [ 0,1] moving peak [ Hilhorst, Vohralik]

       !v1 = 0.8   ! must be same as imod == 9 in convective terms !!!!!
       !v2 = 0.4
       v1 = 0.6   ! must be same as imod == 9 in convective terms !!!!!
       v2 = 0.3

       x0 = 0.35
       y0 = 0.45

       delta = 200.
       fac = 50.

       ev1 = delta * state%model%Re1 * t + 1

       evt = 1./ev1
       evt2 = evt * evt
       evtt = -1./ev1/ev1 * delta * state%model%Re1

       evx = x(1) - x0 -  v1 * t
       evy = x(2) - y0 -  v2 * t

       ev0 = exp( -fac * (evx**2  + evy**2)/ ev1 )

       u = evt * ev0

       ut =  evtt * ev0 + evt * ev0 * (-fac) &
            * ( ( 2* evx *(-v1) + 2*evy *(-v2) )* evt + (evx**2  + evy**2)*evtt)

       du(1) = evt * ev0 * (-fac*2*evx * evt)
       du(2) = evt * ev0 * (-fac*2*evy * evt)

       d2u(1,1) = evt * ( ev0 *(-fac *2)*evt + ev0 * (-fac*2*evx*evt)**2  )
       d2u(2,2) = evt * ( ev0 *(-fac *2)*evt + ev0 * (-fac*2*evy*evt)**2 )

       d2u(1,2) = evt * (                    ev0 * (-fac*2*evx*evt) * (-fac*2*evy*evt) )

       d2u(2,1) = d2u(1,2)


    case(25)   ! circular arc

       u = x(2) / x(1)

       ut =  0.

       du(1) = - x(2) / x(1) / x(1)
       du(2) = 1./ x(1)

       d2u(1,1) = 2. * x(2) / x(1) / x(1) / x(1)
       d2u(2,2) =  0.

       d2u(1,2) = -1. / x(1) / x(1)

       d2u(2,1) = d2u(1,2)

       !write(98,'(20es12.4)') x(1:2), u, ut, d2u(1,1),d2u(2,2)

    case(26)   ! Barenblatt, porus media flow, Radu et all 2008
       m = state%model%param1       !parameter from *.ini file

       rm = (m-1.)/(4*m*m)
       r2 = x(1)*x(1) + x(2)*x(2)
       rt = (t+1)**(-1./m)
       rtt = (-1./m) * (t+1)**(-1./m - 1)


       ev = 1. - rm * r2 * rt

       if(ev < 0) then
          u = 0.
          ut = 0.
          du(:) = 0.
          d2u(:,:) = 0.
       else
          rf = 1./(t+1)
          rft = -1./(t+1)**2

          a0 = m/(m-1)
          a1 = a0 - 1.
          !a2 = a1 - 1.

          ev0 = ev**a0
          ev1 = a0 * ev**a1
          !ev2 = a0 *a1 * ev**a2

          u = ( rf * ev0)**(1./m)

          !write(*,'(20es12.4)') x(1),x(2),rf, ev, 1./m, u
          ut = 1./m * ( rf * ev0)**(1./m-1) &
               *( rft * ev0 &
               + rf * ev1 *(- rm * r2 * rtt ) )

          du(1) = 1./m * ( rf * ev0)**(1./m-1) *rf * ev1 * (-rm * rt * 2*x(1) )
          du(2) = 1./m * ( rf * ev0)**(1./m-1) *rf * ev1 * (-rm * rt * 2*x(2) )
          d2u(:,:) = 0.
          d2u(1,1) = 1E+50  ! it is not important, zero RHS
          d2u(2,2) = 1E+50  ! it is not important, zero RHS
       endif

       !!(1./(1+1)*(1- (2-1)/(4*2*2)*(x*x + y*y)/(1+1)**(1./2))**(1./(2-1)) )**(1./2)

       !write(98,'(20es12.4)') x(1:2), u, ut, d2u(1,1),d2u(2,2)

    case(27)   ! [Raalte 2004]
       if(x(1) <= 0.01 .or. x(2) >= 0.99 ) then
          u = 0.
       else
          u = 1.
       endif
       ! unknown exact solution
       ut = 0.
       du(:) = 0.
       d2u(:,:) = 0.

    case(28) !  NONLINEAR elliptic [Houston, Sulli, Robson 2007]
       r = 20.

       ut = 0.

       evx = x(1)*(1- x(1) )
       evy = x(2)*(1- x(2) )*(1-2*x(2))

       evxd = 1. - 2*x(1)
       evxd2 = -2.

       evyd = 1. - 6*x(2) + 6*x(2)**2
       evyd2 = -6 + 12*x(2)

       ev0 = exp(-r*(2*x(1) -1)**2)
       ev1 = -r*2*(2*x(1) -1) * 2 * ev0
       ev2 = -r*8 * ev0 - r*2*(2*x(1) -1) * 2 * ev1

       u = evx* evy * ev0
       if(t == 0.) u = u*0.9

       du(1) = evy*(evxd * ev0 + evx * ev1)
       du(2) = evx * evyd * ev0

       d2u(1,1) = evy*(evxd2*ev0 + 2* evxd*ev1 + evx * ev2)
       d2u(2,2) = evx * evyd2 * ev0
       d2u(1,2) = evyd*(evxd * ev0 + evx * ev1)
       d2u(2,1) = d2u(1,2)

    case(29) !  NONLINEAR elliptic [Houston, Sulli, Robson 2007] second
       !print*,'....A',x(1:2), t
       if(x(1) == 0.) then
          if(x(2) >  0.) f = pi /2
          if(x(2) <  0.) f = 3.*pi/2
          if(x(2) == 0.) f = 0.
       else
          f = atan(x(2) / x(1) )
          if(x(1) < 0 ) f = f + pi
          if(x(1) >0. .and. x(2) < 0.) then  ! bad quadrant, a correction necessary
             if(abs(x(1)) > abs(x(2) )  ) then
                f = 0.
             else
                f = 3.*pi/2
             end if
          endif
       endif


       ! NEW
       alpha = state%model%param1
       ut = 0.
       r2 = x(1)*x(1) + x(2)*x(2)
       r = r2**0.5

       z = r**alpha

       u = z * sin (alpha * f)
       !if(t == 0.) u = u*0.9

       if( u < 0.) then
          write(*,'(a6,12es14.6)') '???>>',u,z, alpha, f, x(:)
          stop
       endif

       if(r > 0) then
          !rhs = -16./ 81 * sin(2./3*f) / r2 * exp(-4./9 / z)
          rhs = 2.*r**(3*alpha - 4)*alpha**3*(alpha -1) &
               * exp(-r**(2*alpha - 2) * alpha**2) * sin(alpha*f)

          uz = alpha * r**(alpha -1) * sin(alpha * f)
          uf = alpha * r**alpha * cos(alpha * f)

          du(1) = cos( f )* uz  - sin(f) / r * uf
          du(2) = sin( f )* uz  + cos(f) / r * uf
       else
          r = 1E-8
          rhs = 2.*r**(3*alpha - 4)*alpha**3*(alpha -1) &
               * exp(-r**(2*alpha - 2) * alpha**2) * sin(alpha*f)

          uz = alpha * r**(alpha -1) * sin(alpha * f)
          uf = alpha * r**alpha * cos(alpha * f)

          du(1) = cos( f )* uz  - sin(f) / r * uf
          du(2) = sin( f )* uz  + cos(f) / r * uf
       endif


       !if(abs (x(1) ) < 1E-8) then
       !write(*,'(a6,12es12.4)') '## !',x(1:2), rhs,f12, rhs - f12
       !endif

    case(30)   ! linear convection-diffusion, one boundary layer

       c1 = 1./(1. - exp(-state%model%Re) )
       ev = exp( (x(1)-1.) * state%model%Re) - exp(-state%model%Re)
       evx = state%model%Re * exp( (x(1)-1.) * state%model%Re)
       evxd = state%model%Re * evx

       u = x(1) - ev  * c1

       ut =  0.

       du(1) = 1. -evx * c1
       du(2) = 0.

       d2u(1,1) = -evxd * c1
       d2u(2,2) = 0.
       d2u(1,2) = 0.
       d2u(2,1) = d2u(1,2)

       !!write(*,'(20es12.4)') x(1:2), u, du(1), exp(state%model%Re), c1, ev, &
       ! d2u(1,1), -state%model%Re1 * d2u(1,1) + du(1) - 1.

    case(31)    ! LINEAR elliptic. L shape domain [Eibner, Melenk 2007]
       if(x(1) == 0.) then
          if(x(2) >  0.) f = pi /2
          if(x(2) <  0.) f = 3.*pi/2
          if(x(2) == 0.) f = 0.
       else
          f = atan(x(2) / x(1) )
          if(x(1) < 0 ) f = f + pi
          if(x(1) >0. .and. x(2) < 0.) then
             f = f + 2*pi
             ! bad quadrant, a correction necessary
             !if(abs(x(1)) > abs(x(2) )  ) then
             !   f = 0.
             !else
             !   f = 3.*pi/2
             !end if
          endif

          !!r = dot_product(x(1:2), x(1:2))**0.5
          !!write(99,'(8es12.4)') x(1:2), f, r, r*cos(f), r*sin(f), &
          !!     r*cos(f)-x(1), r*sin(f)-x(2)

       endif

       ! NEW
       alpha = 2./3
       r2 = x(1)*x(1) + x(2)*x(2)
       r = r2**0.5

       u = r**(2./3)*sin(alpha*f)*(1-r2*cos(f)**2) *(1-r2*sin(f)**2)
       !z = r**alpha*f

       ut = 0.

       rhs = 2./3*r**(2./3)*(-8*sin(alpha*f)*r**2*cos(f)**2 &
            + 8*r**2*sin(alpha*f)*cos(f)**4 &
            +10*sin(alpha*f)-3*sin(alpha*f)*r**2 &
            +4*r**2*cos(alpha*f)*cos(f)*sin(f) &
            -8*r**2*cos(alpha*f)*cos(f)**3*sin(f))

       ev1 = r**(alpha)
       ev2 = sin(alpha*f)
       ev3 = (1-r2*cos(f)**2)
       ev4 = (1-r2*sin(f)**2)

       !write(*,'(a6,20es14.6)') '???',u, ev1*ev2*ev3*ev4, u - ev1*ev2*ev3*ev4

       uz = ev2*(alpha * r**(alpha - 1)*ev3*ev4 + ev1* (-2.)*r*cos(f)**2 * ev4 &
            + ev1*ev3 * (-2.*r*sin(f)**2 ) )

       uf = ev1*(alpha*cos(alpha * f) * ev3 * ev4 + ev2* ev4 * 2*r2*cos(f)*sin(f) &
            + ev2 * ev3 * (-r2)*2.*sin(f)*cos(f) )

       ! comparison with the output of maple
       !write(*,*), uz, &
       !     -2./3*sin(2./3*f)*(-1+4*r**2-7*r**4*cos(f)**2+7*r**4*cos(f)**4)/(r**(1./3)), &
       !     uz +2./3*sin(2./3*f)*(-1+4*r**2-7*r**4*cos(f)**2+7*r**4*cos(f)**4)/(r**(1./3))
       !write(*,*) uf, &
       !     -2./3*r**(2./3)*(-cos(2./3*f)+r**2*cos(2./3*f)-r**4*cos(2./3*f)*cos(f)**2 &
       !     +r**4*cos(2./3*f)*cos(f)**4+3*r**4*sin(2./3*f)*cos(f)*sin(f) &
       !     -6*r**4*sin(2./3*f)*cos(f)**3*sin(f)), &
       !     uf + 2./3*r**(2./3)*(-cos(2./3*f)+r**2*cos(2./3*f)-r**4*cos(2./3*f)*cos(f)**2 &
       !     +r**4*cos(2./3*f)*cos(f)**4+3*r**4*sin(2./3*f)*cos(f)*sin(f) &
       !     -6*r**4*sin(2./3*f)*cos(f)**3*sin(f))


       du(1) = cos( f )* uz  - sin(f) / r * uf
       du(2) = sin( f )* uz  + cos(f) / r * uf


    case(32)   ! singular corner [0,1] x [ 0,1]
         print*, 'mopdel 32'
       alpha = state%model%param1  ! defined in InitScalarCase( )

       a1 = alpha/2
       a2 = alpha/2-1.
       a3 = alpha/2-2.

       evt = 1.                           ! steady
       evtt = 0.
       if(t == 0.) evt = 0.9

       r = x(1)*x(1)+x(2)*x(2)

       if( r > 0.) then
          ev1 = r**a1
          ev2 = r**a2
          ev3 = r**a3

          if(alpha == 0.) then
             ev1 = 1.
             ev2 = 0.
             ev3 = 0.
          endif

          u = ev1
          ut = 0.


          du(1) = a1* ev2 * 2* x(1)
          du(2) = a1* ev2 * 2* x(2)


          d2u(1,1) = a1*a2*ev3 * 4 *x(1)*x(1) + a1* ev2 * 2

          d2u(2,2) = a1*a2*ev3 * 4* x(2)*x(2) + a1* ev2 * 2

          d2u(1,2) = a1* a2 * ev3* 4*x(1)*x(2)
          d2u(2,1)  = d2u(1,2)

       else
          if(alpha >= 2.) then
             u = 0.
             ut = 0.
             du(:) = 0.
             d2u(:,:) = 0.
          else
             print*,'@@@ echo zero R in (32)in model.f90'
             u = 0.
             w(1) = 1E+5
             if(ityp == 2) return
          endif

       endif

    case(33)   ! singular line [0,1] x [ 0,1]

       alpha = state%model%param1  ! defined in InitScalarCase( )

       a1 = alpha
       a2 = alpha-1.
       a3 = alpha-2.

       evt = 1.                           ! steady
       evtt = 0.
       if(t == 0.) evt = 0.9

       r = x(1)

       if( r > 0.) then
          ev1 = r**a1
          ev2 = r**a2
          ev3 = r**a3

          if(alpha == 0.) then
             ev1 = 1.
             ev2 = 0.
             ev3 = 0.
          endif

          u = ev1
          ut = 0.


          du(1) = a1* ev2
          du(2) = 0.


          d2u(1,1) = a1*a2*ev3

          d2u(2,2) = 0.

          d2u(1,2) = 0.
          d2u(2,1)  = 0.

       else
          if(alpha >= 2.) then
             u = 0.
             ut = 0.
             du(:) = 0.
             d2u(:,:) = 0.
          else
             !print*,'@@@ echo zero R in (33)in model.f90', r, x(1)
             u = 0.
             w(1) = 1E+5
             if(ityp == 2) return
          endif

       endif


    case(34)   ! ! Laplace problem for anisotropic mesh refinement

       u = sin(5*x(1) )

       ut =  0.

       du(1) = 5*cos(5*x(1) )
       du(2) = 0.

       d2u(1,1) = -25*sin(5*x(1) )
       d2u(2,2) =  0.
       d2u(1,2) = 0.
       d2u(2,1) = d2u(1,2)

       !write(98,'(20es12.4)') x(1:2), u, ut, d2u(1,1),d2u(2,2)



    case(35)   ! ! interior layer

       alpha = state%model%param1

       u = atan(alpha*x(1) )

       ut =  0.

       du(1) = alpha/(1 + (alpha * x(1))**2)
       du(2) = 0.

       d2u(1,1) = -2*alpha**3*x(1)/((1+alpha**2 * x(1)**2)**2)
       d2u(2,2) =  0.
       d2u(1,2) = 0.
       d2u(2,1) = d2u(1,2)

       !write(98,'(20es12.4)') x(1:2), u, ut, d2u(1,1),d2u(2,2)

    case(36)   !linear convection diffusion, exponential and parabolic BL, EXACT SOLUTION UNKNOWN
       u = 100* x(1) * x(2) * (1-x(1) ) * (1.-x(2) )   ! homogeneous boundary conditions
       !u = 1.
       du(1) = 1.  ! approximation of the derivative
       du(2) = 0.  ! approximation of the derivative
       d2u(:,:) = 0.

       rhs = 1.

    case(37)   !linear convection diffusion, exponential and parabolic BL, EXACT SOLUTION UNKNOWN
       evx = 1 - x(1)**10
       evy = 1 - x(2)**20

       evxd =  -10. * x(1)**9
       evyd =  -20. * x(2)**19

       evxd2 =   -90. * x(1)**8
       evyd2 =  -380. * x(2)**18


       u = evx * evy
       ut = 0.

       du(1) = evxd * evy
       du(2) = evx  * evyd

       d2u(1,1) = evxd2 * evy
       d2u(1,2) = evxd  * evyd
       d2u(2,1) = evxd  * evyd
       d2u(2,2) = evx   * evyd2

    case(38)
       r = (x(1)*x(1) + x(2) * x(2) )**0.5

       u = 0.
       if(r >= 1./3 .and. r <= 2./3) u = 1.

       ut = 0.
       du(:) = 0.
       d2u(:, :) = 0.

    case(39) !!! triple layer
       u = 0.

       !if(x(2) <= x(1) - 1./8 .and. x(2) >= x(1) - 3./4 .and. x(2) < 0.1) then
       !   u = 1.0
       !endif
       if(x(2) <= x(1) - 1./8 .and. x(2) >= x(1) - 1./2 .and. x(2) < 0.1) then
          !u = max(0., (abs (1. - 2* x(2)) )**0.2 )
          u = 1.0
       endif
       if(x(2) <= x(1) - 1./2 .and. x(2) >= x(1) - 3./4 .and. x(2) < 0.1) then
          !u = max(0., (abs (1. - 2* x(2)) )**0.2 )
          u = 2.0
       endif

       ut = 0.
       du(:) = 0.
       d2u(:, :) = 0.

     case(40) !test for ALG2, jump of normal component of flux reconstruction should be on a face x(1) = 0.5 and x(2) \in [0, 0.5]
       if(x(1) < 0.5) then
        if (x(2) < 0.25) then
          u = 1+ x(1) !1+ x(1)
          du(1) = 1.
          du(2) = 0.
          d2u(:,:) = 0.
        else
          u = 1- x(1) + 2*x(2) !1 + x(1)
          du(1) = - 1. !- 1
          du(2) = 2. !2
          d2u(:,:) = 0.
        endif
       else
        u = 2 + 2*x(1) + 2*x(2) !2 + 2*x(1)
        du(1) = 2. !2
        du(2) = 2. !2
        d2u(:,:) = 0.
       endif

     case(41)  !test for ALG2, FNC estimator when HG nodes are present, quadratic func.

       u = x(1)*(1- x(2))

       ut = 0.

       du(1) = 1- x(2)
       du(2) = - x(1)

       d2u(1,1) = 0
       d2u(2,2) = 0
       d2u(1,2) = -1
       d2u(2,1)  = -1

     case(42)  !test for ALG2, FNC estimator when HG nodes are present, linear func.

       if(x(1) < 0.5) then
        if (x(2) < 0.25) then
          u = 5.25 + x(1) - 16*x(2)
          du(1) = 1.
          du(2) = -16.
          d2u(:,:) = 0.
        else
          u = 1 + x(1) + x(2)
          du(1) = 1.
          du(2) = 1.
          d2u(:,:) = 0.
        endif
       else
        u = 1 + x(1) + x(2)
        du(1) = 1.
        du(2) = 1.
        d2u(:,:) = 0.
       endif

     case(43)  ! Laplace on L-shaped domain [Ainsworth 2005, Robust AEE for NFE approx]

       ! if(x(1) == 0.) then
       !    if(x(2) >  0.) f = pi /2
       !    if(x(2) <  0.) f = 3.*pi/2
       !    if(abs(x(2)) <= 1E-15) f = 0.

       ! elseif(abs(x(2)) <= 1E-15) then
       !    if(x(1) >=  0.) f = 0.
       !    if(x(1) <  0.)  f = pi

       ! else
       !    f = atan(x(2) / x(1) )
       !    if(x(1) < 0 ) f = f + pi
       !    if(x(1) >0. .and. x(2) < 0.) then

       !       !write(*,'(a6,8es12.4)') 'EDW@#@W',x(:),f, f+2*pi, f/2/pi*360, (f+2*pi)/2/pi*360
       !       f = f + 2*pi


       !       !! bad quadrant, a correction necessary
       !       !if(abs(x(1)) > abs(x(2) )  ) then
       !       !   f = 0.
       !       !else
       !       !   f = 3.*pi/2
       !       !endif
       !    endif

       ! endif

        if(r > 0) then
          f = atan2(x(2), x(1) )
          if(f <= 0. ) f = f + 2 * pi
          !!!!if(f >= 2*pi - 1E-6 ) f = f - 2 * pi

          if(abs(x(1)) <= 1E-14 .and. x(2) >  0. )then
             f = pi /2
          elseif(abs(x(1)) <= 1E-12 .and. x(2) <  0. )then
             f = 3.*pi/2
          endif

          if( f > 2*pi - 1E-1  ) f = f- 2* pi ! the bad quandrant

          !!if(abs( f - f1) > 1E-8) write(*,'(a6,5es14.6)') '$$$',x(:), f,f1, f - f1
          !write(*,*) 'model.f90:::::',x(:), f,f/2/pi*360
       else
          f = 0.
       endif


       alpha = 2./3
       r2 = x(1)*x(1) + x(2)*x(2)
       r = r2**0.5

       u = r**(alpha)*sin(alpha*f)

       ut = 0.

       ev1 = r**(alpha)
       ev2 = sin(alpha*f)

       if(r > 0) then
          uz = ev2*(alpha * r**(alpha - 1) )

          uf = ev1*(alpha*cos(alpha * f) )

          du(1) = cos( f )* uz  - sin(f) / r * uf
          du(2) = sin( f )* uz  + cos(f) / r * uf

       else
          du(1) =1E+10
          du(2) =1E+10
       endif


       !if(x(1) >0. .and. abs(x(2)) < 1.E-4) &
       !     write(*,'(a6,8es12.4)') 'EDW@#@W',x(:),f,f/2/pi*360,u


       ! RHS = 0, evaluation is not necessary

    case(44)  !smooth bounded solution (sin)
       evx = sin(2 * pi * x(1) )
       evy = sin(2 * pi * x(2) )

       evxd = 2 * pi * cos(2 * pi * x(1))
       evyd = 2 * pi * cos(2 * pi * x(2))

       evxd2 = -4 * pi * pi * sin(2 * pi * x(1))
       evyd2 = -4 * pi * pi * sin(2 * pi * x(2))

       u = evx * evy
       ut = 0.

       du(1) = evxd * evy
       du(2) = evx  * evyd

       d2u(1,1) = evxd2 * evy
       d2u(2,2) = evx   * evyd2
       d2u(1,2) = evxd  * evyd
       d2u(2,1)  = d2u(1,2)

    case(45)   !poisson problem with two identical parabolic boundary layers
       evx = 1 - x(1)**10
       evy = 1 - x(2)**10

       evxd =  -10. * x(1)**9
       evyd =  -10. * x(2)**9

       evxd2 =  -90. * x(1)**8
       evyd2 =  -90. * x(2)**8


       u = evx * evy
       ut = 0.

       du(1) = evxd * evy
       du(2) = evx  * evyd

       d2u(1,1) = evxd2 * evy
       d2u(1,2) = evxd  * evyd
       d2u(2,1) = evxd  * evyd
       d2u(2,2) = evx   * evyd2

    case(46)  !! Laplace on LL-shaped domain [Vohralik ESAIM]

       if(abs(x(1)) <= 1E-14 .and. abs(x(2)) <= 1E-14  )then
          f = 0.
       elseif(abs(x(1)) <= 1E-14 .and. x(2) >  0. )then
          f = pi /2

       elseif(abs(x(1)) <= 1E-12 .and. x(2) <  0. )then
          f = 3.*pi/2

       elseif( abs(x(2)) <= 1.0E-12 .and. x(1) >=  0 ) then
          f = 2 * pi

       elseif( abs(x(2)) <= 1.0E-14 .and. x(1) <  0. ) then
          f = pi

       else
          f = atan(x(2) / x(1) )
          if(x(1) < 0 ) f = f + pi
          if(x(1) >0. .and. x(2) < 0.) then
             f = f + 2*pi
          endif

       endif

       ! !if( x(1) > 0.55 .and. x(1) < 0.7 .and. abs(x(2)) < 1E-3) &
       ! !     print*,'$$$',x(:),f,  atan(x(2) / x(1) )
       f1 = f

       alpha = 2./3
       r2 = x(1)*x(1) + x(2)*x(2)
       r = r2**0.5

       if(r > 0) then
          f = atan2(x(2), x(1) )
          if(f <= 1E-6 ) f = f + 2 * pi

          if(abs( f - f1) > 1E-8) write(*,'(a6,5es14.6)') '$$$',x(:), f,f1, f - f1
          !write(20,*) x(:), f,f1
       else
          f = 0.
       endif

       u = r**(alpha)*sin(alpha*f)

       ev1 = (x(1) - r*cos(f))**2 + (x(2) - r*sin(f))**2
       if(ev1 > 1e-18) print*,'###### echo in model.f90  ###########',x(:),r,f,ev1

       ut = 0.

       ev1 = r**(alpha)
       ev2 = sin(alpha*f)

       if(r > 0) then
          uz = ev2*(alpha * r**(alpha - 1) )

          uf = ev1*(alpha*cos(alpha * f) )

          du(1) = cos( f )* uz  - sin(f) / r * uf
          du(2) = sin( f )* uz  + cos(f) / r * uf

          ! RHS = 0, evaluation is not necessary
       else
          ! used only for visualization, (I hope :-) )
          !du(1:2) = 1E+40
          du(1:2) = 0.
       endif

    case(47)   ! linear convection-diffusion, , parabolic and exponential BLs

       c1 = 1./(1. - exp(-state%model%Re) )
       ev = exp( (x(1)-1.) * state%model%Re) - exp(-state%model%Re)
       evx = state%model%Re * exp( (x(1)-1.) * state%model%Re)
       evxd = state%model%Re * evx

       ! EXACT SOLUTION UNKNOWN
       u = (x(1) - ev  * c1) * x(1)*x(2)*(1-x(1)) *(1-x(2))

       ut =  0.

       du(1) = 1. -evx * c1
       du(2) = 0.

       d2u(1,1) = -evxd * c1
       d2u(2,2) = 0.
       d2u(1,2) = 0.
       d2u(2,1) = d2u(1,2)

       !!write(*,'(20es12.4)') x(1:2), u, du(1), exp(state%model%Re), c1, ev, &
       ! d2u(1,1), -state%model%Re1 * d2u(1,1) + du(1) - 1.

    case(48) ! linear test case for STDGM exp in time, polynomial in space
       m = 10.
       !m = 20. * Pi
       evx = x(1) - x(1)**2.
       evy = x(2) - x(2)**2.
       evt = 16. * (exp(m*t) - 1. ) / (exp(m) - 1.)
       evtd =16. * ( m * exp(m*t) ) / (exp(m) - 1.)
      ! evt = 16. * cos(m *t)
      ! evtd =16. * m * (-1.) * sin(m * t)
       u = evx * evy * evt
       ut = evx * evy * evtd

       du(1:2) = 0.
       evxd = (1. - 2.*x(1))
       evyd = (1. - 2.*x(2))
       du(1) = evt * evxd * evy
       du(2) = evt * evyd * evx

       d2u(:,:) = 0.
       d2u(1,1) = (-2.) * evt * evy
       d2u(2,2) = (-2.) * evt * evx
       d2u(1,2) = evt * evxd * evyd
       d2u(2,1) = d2u(1,2)
!       ! END
!       m=12.
!       evt = m * t**7
!       evtd = 7 * m * t**6
!
!       evx = sin(x(1)*x(2))
!
!       u = evt * evx
!       ut = evtd * evx
!       du(1:2) = 0.
!
!       du(1) = evt * cos(x(1)*x(2)) * x(2)
!       du(2) = evt * cos(x(1)*x(2)) * x(1)
!
!       d2u(:,:) = 0.
!       d2u(1,1) = evt * (-1.) * sin(x(1)*x(2)) * x(2)**2
!       d2u(2,2) = evt * (-1.) * sin(x(1)*x(2)) * x(1)**2
!       d2u(1,2) = evt * (-1.) * sin(x(1)*x(2)) * x(1) * x(2)
!       d2u(2,1) = d2u(1,2)


    case(49) ! STIFF EQUATION for STDGM
    !CASE SUPER STIFF should be used with paramet > 1000
       z = 2.0 * pi
       phi = cos((pi * 0.25) + z * t)
       Dphi = (-1) * sin((pi * 0.25) + z * t ) * z

       u = phi
       ut = Dphi

       m = state%model%param1

       rhs = m * phi + Dphi


       du(1) = 0.
       du(2) = 0.
       d2u(:,:) = 0.
!       END OF SUPER STIFF

!!!        MODERATE STIFF
!       m = state%model%param1
!!       phi = sin((pi * 0.25) + t)
!!     !  Dphi = cos((pi * 0.25) + t )

!       u = ( (-1.)/(1. - m) ) * exp( (-1.)*t )
!       ut = (-1.) * u

!       rhs = exp(-t)
!!
!!
!       du(1) = 0.
!       du(2) = 0.
!       d2u(:,:) = 0.
!!!        END



    case(50) ! nonconstant linear convection

       c1 = 1./(1. - exp(-1.0D+03) )
       ev = exp( (x(1)-1.) * 1.0D+03) - exp(-1.0D+03)

       u = (x(1) - ev  * c1) * x(1)*x(2)*(1-x(1)) *(1-x(2))
       !print*,'@@@@',x(:), c1, ev, u

       du(1) = 0.
       du(2) = 0.
       d2u(:,:) = 0.

       rhs = 0.

    case(51)        ! STDGM STIFF SYSTEM WITH LAPLACE
    !     ! u' = \lapl u - \lapl f + f'

    !NEPOCITA ???
!
!    m = 2. !2.* Pi
!	 u = sin(m*t*x(1)*x(2))
!
!    ut = cos( m*t*x(1)*x(2) ) *x(1)*x(2)*m
!
!    du(1:2) = 0.
!    du(1) = cos( m*t*x(1)*x(2) ) *m*t*x(2)
!    du(2) = cos( m*t*x(1)*x(2) ) *m*t*x(1)
!
!
!    d2u(:,:) = 0.
!    d2u(1,1) = ( m*t*x(2) )**2. * (-1.)*sin( m*t*x(1)*x(2) )
!    d2u(2,2) = ( m*t*x(1) )**2. * (-1.)*sin( m*t*x(1)*x(2) )
!    d2u(1,2) = (m*t)**2. * x(1)*x(2)  * (-1.)*sin(m*t*x(1)*x(2)) + m * t * cos(m*t*x(1)*x(2))
!    d2u(2,1) = d2u(1,2)
!
!    f = ut !f'
!    ff = (d2u(1,1) + d2u(2,2) ) / state%model%Re ! lapl f
!    rhs = f - ff

    ! m = 2. * Pi
     m = 8.
     ev = 2. * Pi * ( x(1) + x(2) )
     evx = 2. * Pi !* x(2)
     evy = 2. * Pi !* x(1)
     !evt = exp( t )
    ! evt = cos(m * t) !(exp(m*t) - 1. ) / (exp(m) - 1.)
    ! evtd =(-1.) * m * sin(m*t) ! (  m * exp(m*t) ) / (exp(m) - 1.) !
     evt = (exp(m*t) - 1. ) / (exp(m) - 1.)
     evtd = (  m * exp(m*t) ) / (exp(m) - 1.)
!
     u =  evt * sin(ev)
     ut = evtd * sin(ev)
!
!   !  f =  2. * Pi * (x(1) + x(2)) * cos(ev)    !f'
!   !  ff = (-2.) * (2 * Pi * t)**(2.) * sin(ev)  !ff
!
!    ! ut = 2. * Pi * (x(1) + x(2)) * cos(ev)

!    ! m = state%model%param1

!     !rhs = (-1.) * ff + f
!
     du(1) = evx * cos(ev) * evt
     du(2) = evy * cos(ev) * evt
     d2u(:,:) = 0.

     d2u(1,1) = evx * evx * (-1.) * sin(ev) * evt
     d2u(2,2) = evy * evy * (-1.) * sin(ev) * evt
     d2u(1,2) = evx * evy * (-1.) * sin(ev) * evt !+  2. * Pi * cos(ev) * evt
     d2u(2,1) = d2u(1,2)

    case(52) ! nonconstant linear convection

       u = exp(x(1) * x(2) )
       ut = 0.

       du(1) = x(2) * u
       du(2) = x(1) * u

       d2u(1,1) = x(2) * x(2) * u
       d2u(1,2) = u + x(1) * x(2) * u
       d2u(2,1) = d2u(1, 2)
       d2u(2,2) = x(1) * x(1) * u
       !print*,'###',x(1:2), u, du(:), d2u(:,:)

    case(53)  !boundary line singularity
       m = state%model%param1
       xx = max(x(1), 0. ) ! 1E-10)
       u =  xx**m
       ut = 0.

       du(1) = m * xx**(m-1)
       du(2) = 0.

       d2u(1,1) = m* (m-1) * xx**(m-2)
       d2u(1,2) = 0.
       d2u(2,1) = 0.
       d2u(2,2) = 0.

       if(xx == 0.) then
          if(m > 1) then
             du(1) = 0.
          else
             du(1) = 1E+20
          endif

          if(m > 2) then
             d2u(1,1) = 0.
          else
             d2u(1,1) = 1E+20
          endif

       endif


       !write(*,'(a4,18es12.4)')'###',x(1:2), u, du(:), d2u(:,:)


    case(54)  !wave front
       r0 = 0.25
       m = 50.
       xc = 0.5
       yc = 0.5

       ev1 = (x(1) - xc)**2 + (x(2) - yc)**2

       ev1x = 2*(x(1) - xc)
       ev1y = 2*(x(2) - yc)

       ev1xx = 2.
       ev1yy = 2.

       ev2  = m*( ev1**0.5 - r0 )
       ev2x = m * 0.5* ev1**(-0.5) * ev1x
       ev2y = m * 0.5* ev1**(-0.5) * ev1y

       ev2xx = m * 0.5* ( -0.5*ev1**(-1.5) * ev1x *ev1x +  ev1**(-0.5) * ev1xx)

       ev2xy = m * 0.5* ( -0.5*ev1**(-1.5) * ev1y *ev1x  )

       ev2yy = m * 0.5* ( -0.5*ev1**(-1.5) * ev1y *ev1y +  ev1**(-0.5) * ev1yy)


       u =  atan(ev2 )
       ut = 0.

       du(1) = 1./ ( 1 + ev2**2) * ev2x
       du(2) = 1./ ( 1 + ev2**2) * ev2y



       d2u(1,1) = -1./(1. + ev2**2)**2. * 2.* ev2 * ev2x * ev2x + 1./ ( 1 + ev2**2) * ev2xx
       d2u(1,2) = -1./(1. + ev2**2)**2. * 2.* ev2 * ev2y * ev2x + 1./ ( 1 + ev2**2) * ev2xy
       d2u(2,1) = d2u(1,2)
       d2u(2,2) = -1./(1. + ev2**2)**2. * 2.* ev2 * ev2y * ev2y + 1./ ( 1 + ev2**2) * ev2yy


       !write(*,'(a4,18es12.4)')'###',x(1:2), u, du(:), d2u(:,:)

    case(55)  !interior line singularity
       alpha = state%model%param1
       beta = 0.6

       ev1 = cos(pi * x(2) / 2)
       ev1y  = -sin(pi * x(2) / 2) * pi / 2
       ev1yy = -cos(pi * x(2) / 2) * pi *pi / 4

       xx = x(1) - beta*(x(2) + 1)
       if(xx > 0) then
          ev2  = (xx)**alpha
          ev2x = alpha * xx**(alpha - 1)
          ev2y = alpha * xx**(alpha - 1) *(-beta)

          ev2xx = alpha * (alpha - 1.) * xx**(alpha - 2)
          ev2xy = alpha * (alpha - 1.) * xx**(alpha - 2) *(-beta)
          ev2yy = alpha * (alpha - 1.) * xx**(alpha - 2) *beta**2
       else
          ev2 = 0.
          ev2x = 0.
          ev2y = 0.

          ev2xx = 0.
          ev2xy = 0.
          ev2yy = 0.
       endif

       u =  ev1 + ev2
       ut = 0.

       du(1) =       ev2x
       du(2) = ev1y + ev2y



       d2u(1,1) =        ev2xx
       d2u(2,1) =        ev2xy
       d2u(1,2) =        ev2xy
       d2u(2,2) = ev1yy + ev2yy


       !write(*,'(a4,18es12.4)')'###',x(1:2), u, du(:), d2u(:,:)

    case(56)  ! nonlinear hyperbolic problem

       ev = sin(2* pi *(x(1) + x(2) - 2 *t) )

       u = ev  !+ 1.5

       ev0 = 2 *pi* cos(2* pi *(x(1) + x(2) - 2 *t) )

       ut = -2. * ev0

       du(1) =  ev0
       du(2) =   ev0

       d2u(1,1) = -2 * pi * u
       d2u(2,2) = -2 * pi * u
       d2u(1,2) = -2 * pi * u
       d2u(2,1)  = d2u(1,2)


    case(57)   ! rotating peak

       delta = -100.

       r0 = 0.5

       !xc = r0 * cos(2*pi*t)
       !yc = r0 * sin(2*pi*t)

       xc = r0 * cos(t)
       yc = r0 * sin(t)

       evx = (x(1) - xc)
       evy = (x(2) - yc)

       u = exp( delta * (evx**2  + evy**2) )

       !ut = u * delta *( 2*evx*r0*sin(2*pi*t) *2*pi- 2*evy*r0*cos(2*pi*t)*2*pi )
       ut = u * delta *( 2*evx*r0*sin(t) - 2*evy*r0*cos(t) )

       du(1) = u * delta* 2 * evx
       du(2) = u * delta* 2 * evy

       d2u(1,1) = delta*( u * delta * 2* evx * 2 * evx + u * 2 )
       d2u(2,2) = delta*( u * delta * 2* evy * 2 * evy + u * 2 )

       d2u(1,2) = delta* u * delta * 2* evx * 2 * evy
       d2u(2,1) = d2u(1,2)

    case(58)  !  multiple difficulties : L-shape + interior wave + peak + boundary layer
       p_L  = 1.
       p_w  = 1.
       p_p  = 1.
       p_BL = 1.

       ut = 0.
       u = 0.
       du(:) = 0.
       d2u(:,:) = 0.

       !! L-shape singularity (as case 43)
       alpha = 2./3
       r2 = x(1)*x(1) + x(2)*x(2)
       r = sqrt(r2)

       if(r > 0) then
          f = atan2(x(2), x(1) )
          if(f <= 0. ) f = f + 2 * pi

          if(abs(x(1)) <= 1E-14 .and. x(2) >  0. )then
             f = pi /2
          elseif(abs(x(1)) <= 1E-14 .and. x(2) <  0. )then
             f = 3.*pi/2
          endif

          if( f > 2*pi - 1E-1  ) f = f- 2* pi ! the bad quandrant

          !!if(abs( f - f1) > 1E-8) write(*,'(a6,5es14.6)') '$$$',x(:), f,f1, f - f1
          !write(*,*) 'model.f90:::::',x(:), f,f/2/pi*360
       else
          f = 0.
       endif

       u = u + p_L * r**(alpha)*sin(alpha*f)

       ev1 = r**(alpha)
       ev2 = sin(alpha*f)

       if(r > 0) then
          uz = ev2*(alpha * r**(alpha - 1) )

          uf = ev1*(alpha*cos(alpha * f) )

          du(1) = du(1) + p_L * (cos( f )* uz  - sin(f) / r * uf)
          du(2) = du(2) + p_L * (sin( f )* uz  + cos(f) / r * uf)

          !d2u(:,:)= 0. ! the 2nd derivatives are non-zero but their do not bring a contribution to RHS

          ! RHS = 0, evaluation is not necessary
       else
          ! used only for visualization, (I hope :-) )
          !du(1:2) = 1E+40
          !du(1:2) = 0.
       endif

       ! interior wave
       r0 = 0.75
       m = 200.
       xc = 0.
       yc = -0.75

       ev1 = (x(1) - xc)**2 + (x(2) - yc)**2

       ev1x = 2*(x(1) - xc)
       ev1y = 2*(x(2) - yc)

       ev1xx = 2.
       ev1yy = 2.

       if(ev1 <= 0) ev1 = 1E-14
       ev2  = m*( ev1**0.5 - r0 )
       ev2x = m * 0.5* ev1**(-0.5) * ev1x
       ev2y = m * 0.5* ev1**(-0.5) * ev1y

       ev2xx = m * 0.5* ( -0.5*ev1**(-1.5) * ev1x *ev1x +  ev1**(-0.5) * ev1xx)

       ev2xy = m * 0.5* ( -0.5*ev1**(-1.5) * ev1y *ev1x  )

       ev2yy = m * 0.5* ( -0.5*ev1**(-1.5) * ev1y *ev1y +  ev1**(-0.5) * ev1yy)


       u =  u + p_w * atan(ev2 )

       du(1) = du(1) + p_w/ ( 1 + ev2**2) * ev2x
       du(2) = du(2) + p_w/ ( 1 + ev2**2) * ev2y

       d2u(1,1) = d2u(1,1) +p_w*(-1./(1. + ev2**2)**2.* 2.*ev2 *ev2x *ev2x + 1./ ( 1 + ev2**2) * ev2xx)
       d2u(1,2) = d2u(1,2) +p_w*(-1./(1. + ev2**2)**2.* 2.*ev2 *ev2y *ev2x + 1./ ( 1 + ev2**2) * ev2xy)
       d2u(2,1) = d2u(1,2)
       d2u(2,2) = d2u(2,2) +p_w*(-1./(1. + ev2**2)**2.*2.* ev2 *ev2y *ev2y + 1./ ( 1 + ev2**2) * ev2yy)


       ! PEAK
       m = -1000.
       xc = -sqrt(5.)/4
       yc = -1./4

       ev1 = (x(1) - xc)**2 + (x(2) - yc)**2
       ev0 = exp( m * ev1)

       u = u + p_p * ev0
       du(1) = du(1) + p_p * ev0 * m * 2*(x(1) - xc)
       du(2) = du(2) + p_p * ev0 * m * 2*(x(2) - yc)

       d2u(1,1) = d2u(1,1) + p_p *( 2. * m * ev0 + (2.* m *(x(1) - xc))**2 * ev0)
       d2u(1,2) = d2u(1,2) + p_p *( (2.* m *(x(1) - xc))* (2.* m *(x(2) - yc)) * ev0 )
       d2u(2,1) = d2u(1,2)
       d2u(2,2) = d2u(2,2) + p_p * (2. * m * ev0 + (2.* m *(x(2) - yc))**2 * ev0)


       !print*,'de3de3d3',u,ut, du(:), d2u(:)
       ! exponential boundary layer
       m = - 100.
       ev0 = exp(m *(1+ x(2) ) )
       u = u + p_BL * ev0

       du(2) = du(2) + p_BL* m * ev0

       d2u(2,2) = d2u(2,2) + p_BL *  m*m * ev0


       !write(*,*) 'de3de3',x(1:2), u

    case(59)  ! battery

       !u =    x(2) *(24 - x(2) ) * (1. - x(1)/8.4 )  + 1
       u = 1.
       !if(x(1) < 1E-6) u = 1

       ut = 0.

       du(1) =  0.  !0.01 * x(2) *(24 - x(2) )/ (-8.4)
       du(2) =  0.  !0.01 * (24 - 2*x(2))  * (1. - x(1)/8.4 )

       d2u(1,1) = 0.
       d2u(2,2) = 0.
       d2u(1,2) = 0.
       d2u(2,1)  = d2u(1,2)

       call Set_Battery(x(1),x(2), ev1, ev2, rhs)

       !if(x(1) < 1. .and. x(2) < 1. )   print*,'####', x(:), u, rhs

    case(60) ! DWR TEST CASE - peak of the solution at (0.25, 0.25)
      ! u = log(r(x)), where r(x) = sqrt( |x-x_peak|^2 + eps^2 )
      ! we use log(sqrt(x)) = 1/2log(x)

!       x_peak = (/ 0.25, 0.25 /) !square (0,1)
       c1 = 0.0625
       x_peak = (/ 0.25 - c1, 0.25 - c1  /)
!       x_peak = (/ -0.5, 0.5 /)  ! LL-shaped

         ! should be dependent on the mesh; if to small - peak will not be visible
!       ev = 0.03 !1 ! epsilon
       ev = c1  / 2.0 !1 ! epsilon

       r0 = Pi*ev*ev
       rt = distance( x , x_peak )
       z1 = x(1) - x_peak(1)
       z2 = x(2) - x_peak(2)
       c1 = rt*rt + ev*ev ! |x-a|^2 + eps^2


       u = (-1)*0.5* log( c1  )

       ut = 0.0

       du(1) = (-1)*( z1 ) / ( c1 ) ! ( rt*rt + ev*ev)**(-1.0) * 2.0 * (x(1) - x_peak(1))
       du(2) = (-1)*( z2 ) / ( c1 ) ! ( rt*rt + ev*ev)**(-1.0) * 2.0 * (x(2) - x_peak(2))
!
       d2u(1,1) =  (-1)*(c1 - 2.0*z1*z1 ) / (c1*c1)
       d2u(2,2) =  (-1)*(c1 - 2.0*z2*z2) / (c1*c1)
       d2u(1,2) = (-1)*( -2.0*z1*z2 ) / (c1*c1)
       d2u(2,1) = d2u(1,2)

    case(61)  ! battery

       !u =    x(2) *(24 - x(2) ) * (1. - x(1)/8.4 )  + 1
       u = 1.
       !if(x(1) < 1E-6) u = 1

       ut = 0.

       du(1) =  0.  !0.01 * x(2) *(24 - x(2) )/ (-8.4)
       du(2) =  0.  !0.01 * (24 - 2*x(2))  * (1. - x(1)/8.4 )

       d2u(1,1) = 0.
       d2u(2,2) = 0.
       d2u(1,2) = 0.
       d2u(2,1)  = d2u(1,2)

       call Set_Battery_Simplified(x(1),x(2), ev1, ev2, rhs)

       !if(x(1) < 1. .and. x(2) < 1. )   print*,'####', x(:), u, rhs

   case(62) ! DWR TEST CASE - peak of the solution at (0.25, 0.25),
      ! connected with subdomain in twoSquaresIn250.subgrid
      ! UNKNOWN EXACT SOLUTION
      ! unknown solution
      ! rhs
      u = 0.0
      ut = 0.0
      du(:) = 0.0
      d2u(:,:) = 0.0
!      ! the volume of the subdomain in twoSquaresIn250.subgrid
      r0 = 0.125*0.125
     ! the volume of the subdomain in twoSquaresIn164.subgrid
!      r0 = 0.01*0.01
      rhs = 1.0  / r0


      ! we use log(sqrt(x)) = 1/2log(x)

       ! should be dependent on the mesh; if to small - peak will not be visible
      ! for finer mesh e.g. square01-025075.grid
!       ev = 0.01 !1 ! epsilon
!       x_peak = (/ 0.25, 0.25 /) !square (0,1)

      ! for square01.4uns.250
!       ev = 0.1 !1 ! epsilon
!      ! x_peak = (/ 0.25, 0.25 /) - (/ ev*0.5, ev*0.5 /)
!
!      ! for LL-shape.str.N8.grid
!
!       x_peak = (/ -0.5, 0.5 /)  ! LL-shaped
!
!
!       r0 = Pi*ev*ev
!       rt = distance( x , x_peak )
!       z1 = x(1) - x_peak(1)
!       z2 = x(2) - x_peak(2)
!       c1 = rt*rt + ev*ev ! |x-a|^2 + eps^2
!
!
!        u = 1.0
!        ut = 0.0
!        du(:) = 0.0
!        d2u(:,:) = 0.0
!
!
!       if ( rt < ev  ) then
!!         print*, 'true True'
!          rhs = 1.0 / r0
!!          rhs = 1.0 / ( ev*ev )
!!          print*, 'rhs:', rhs
!       else
!          rhs = 0.0
!       endif

!       if ( ( abs( x_peak(1) - x(1) ) < ev*0.5  ) .and. &
!            ( abs( x_peak(2) - x(2) ) < ev*0.5  ) ) then
!         rhs = 1.0 / ( ev*ev )
!!          print*, 'rhs:', rhs
!       else
!          rhs = 0.0
!       endif
    ! test case for DWR target functional dudx
    ! used also for computation of the problem on the domain CROSS
    case(63)
      ! unknown solution
      ! rhs - dw/dx
      u = 0.0
      ut = 0.0
      du(:) = 0.0
      d2u(:,:) = 0.0

      rhs = 1.0 !/ r0


    case(64) !  ! porous media flow,  Forchheimer 2-term law

       evt = exp(-2.* t)

       evx = x(1)*(1- x(1) )
       evy = x(2)*(1- x(2) )

       evxd = 1. - 2*x(1)
       evxd2 = -2.

       evyd = 1. - 2*x(2)
       evyd2 = -2

       u =  evt * evx* evy
       ut = -2.* u

       du(1) =  evt * evxd * evy
       du(2) =  evt * evx  * evyd

       d2u(1,1) =  evt * evxd2 * evy
       d2u(1,2) =  evt * evxd  * evyd
       d2u(2,1) =  d2u(1,2)
       d2u(2,2) =  evt * evx   * evyd2


    case(65) !  ! porous media flow,  Darcy /Forchheimer 2-term law

       u =  0.  !0.1 ! + 1E+00 * sin(5*pi * x(1)) * sin(5 * pi * x(2))
       !if( abs(x(1)) <1E-6) u = 1.E+3
       !if(  x(1)  < -0.999 .and. x(2) >= -0.1000001  .and. x(2) <= 0.100001 ) u = 1.E+3
       if(  x(2)  > 0.999 .and. x(1) >= -0.3000001  .and. x(1) <= 0.000001) then
          !r = 0.05
          !fac = 1.

          !! snmoothing of IC
          !if(x(2) -1. < r * (x(2) - 1.) ) then
          !  fac = (sin((x(2)-(1+fac/2))/fac*pi )+1 ) / 2
          !endif
          !!write(*,*) x(1:2), t, fac

          fac = 1.
          if(t <= 1) fac = 2*log(1+t)
          if(t >= 2)  fac = exp(-2*(t-2))
          fac = min(1., fac)

          !if(x(1) > -0.01) print*,'###',t,fac
          !!!if(t > 1) u = u * exp(-5*(t-1))

          u = 1.E+3 * fac


       endif


       rhs = 0.
       du(:) = 0.
       d2u(:,:) = 0.

    case(66:)
       stop 'UNKNOWN type in set_Model_Data'
    end select

    if(ityp == 1) then
       ! output of the exact solution
       w(1) = u

       !write(98,'(20es12.4)') x(1:2), w(1), ut, d2u(1,1),d2u(2,2) !, &
    elseif(ityp == 2) then        ! output of the right-hand-side

       if(imod == 8 .or. imod == 9 .or. imod == 21 .or. imod == 22 &
            .or. imod == -24  .or. imod == 25 .or. imod == 26      &
            .or. imod == 27 .or. imod == 38 .or. imod == 39 .or. imod == 43 &
            .or. imod == 46) then ! zero RHS
          w(1) = 0.

       elseif(imod == 4 .or. imod == 50  ) then ! zero right hand side
          w(1) = 0.

       elseif(imod == 47 ) then ! zero righ hand side
          w(1) = 1.

       !elseif(imod == 49 .or. imod == 51) then !49,51 - STIFF systems
       elseif(imod == 49) then !49 - STIFF systems
          w(1) = rhs

       elseif(imod == 59 .or. imod == 61 .or. imod == 65)  then ! 59, battery
          w(1) = rhs

       !DWR peak (TRY) UNKNOWN SOL, du/dx unknown sol
       elseif(imod == 62 ) then
          w(1) = rhs

       !DWR du/dx unknown sol
       elseif( imod == 63 ) then
          w(1) = rhs

       elseif(imod == 29 .or. imod == 31 .or. imod == 36 ) then ! direct computation
          w(1) = rhs

       elseif(imod == 28 .or. imod == 64) then ! viscous coefficient depends on | \nabla u |
          if(state%model%Re /= 0.) then
             Re_1(1) = 1./state%model%Re
          else
             Re_1(1) = 0.
          endif

          w(1) =  ut !! time derivative

          !! general convection
          do i=1,nbDim
             w(1) = w(1) + Eval_Convective_Coeffs(u, i, 1, 1, x(1:nbDim) ) * du(i)
          enddo

          !! general (nonlinear) diffusion: $ -\nabla (\mu(|\nabla u|) \nabla u)$
          do i=1,nbDim
             j = mod(i,nbDim) + 1

             w(1) = w(1) - Eval_Diffusion_Coeffs(u, du(1:nbDim), i, i, 1, 1, Re_1, 0,x(1:2))* d2u(i,i) &
                     - Eval_Diffusion_Coeffs(u, du(1:nbDim), i, i, 1, 1, Re_1, 1,x(1:2)) * du(i) &
                     * 2 * (du(i) * d2u(i,i) + du(j)* d2u(i,j) ) / (2*(du(1)**2 + du(2)**2)**0.5)

             !write(*,'(a6,2i5, 8es12.4)')'SWE#:',i,j, w(1), &
             !     Eval_Diffusion_Coeffs(u, du(1:nbDim), i, i, 1, 1, Re_1, 0,x(1:2))* d2u(i,i) &
             !     - Eval_Diffusion_Coeffs(u, du(1:nbDim), i, i, 1, 1, Re_1, 1,x(1:2)) * du(i) &
             !     * 2 * (du(i) * d2u(i,i) + du(j)* d2u(i,j) ) / (2*(du(1)**2 + du(2)**2)**0.5), &
             !     Eval_Diffusion_Coeffs(u, du(1:nbDim), i, i, 1, 1, Re_1, 0,x(1:2)) , &
             !     - Eval_Diffusion_Coeffs(u, du(1:nbDim), i, i, 1, 1, Re_1, 1,x(1:2)),&
             !     du(1:2), d2u(i,i)

          enddo

          !! reaction
          if(state%model%ireac > 0)  w(1) = w(1) + Eval_Reaction_Coeffs(u, du(1:2), x(1:nbDim),0 )


       else

          if(state%model%Re /= 0.) then
             Re_1(1) = 1./state%model%Re
          else
             Re_1(1) = 0.
          endif

          !! time derivative
          w(1) =  ut

          ! FEM: linear convection
          !w(1) = w(1) +Scal_Lin(1, x(1), x(2), t)*du(1)+Scal_Lin(2, x(1), x(2), t)*du(2)
          ! FEM: linear reaction
          !w(1) = w(1) + Scal_Lin(0, x(1), x(2), t)*u

          !! general convection
          do i=1,nbDim
             w(1) = w(1) + Eval_Convective_Coeffs(u, i, 1, 1, x(1:nbDim) ) * du(i)
          enddo

          !! general (nonlinear) diffusion
          if(state%model%idiff > 0) then
             do i=1,nbDim
                do j=1,nbDim
                   w(1) = w(1) - Eval_Diffusion_Coeffs(u, du(1:nbDim), i, j, 1, 1, Re_1, 0, x(1:2)) &
                        * d2u(i,j) &
                        - Eval_Diffusion_Coeffs(u, du(1:nbDim), i, j, 1, 1, Re_1, 1, x(1:2)) &
                        * du(i) * du(j)
                enddo
             enddo
          endif

          !! reaction
          if(state%model%ireac > 0)  w(1) = w(1) + Eval_Reaction_Coeffs(u, du(1:2), x(1:nbDim),0 )


          !write(*,'(20es12.4)') x(1:2), w(1)
          !write(400+state%time%iter,'(20es12.4)') x(1:2), w(1)


          !write(500+state%time%iter,'(20es12.4)') x(1:2), w(1), u, ut, du(1:2), d2u(1:2,1:2)!, &

          !!if(abs(w(1) ) > 1E-8) then
          !write(102,'(20es14.6)') x(1:2), w(1), u, ut, du(1:2), d2u(1:2,1:2), &
          !     Eval_Convective_Coeffs(u, 1, 1, 1) * du(1), &
          !     Eval_Convective_Coeffs(u, 2, 1, 1) * du(2), &
          !     Eval_Diffusion_Coeffs(u, du(1:nbDim), 1, 1, 1, 1, Re_1, 0), &
          !     Eval_Diffusion_Coeffs(u, du(1:nbDim), 1, 1, 1, 1, Re_1, 0)* d2u(1,1)
               !Eval_Diffusion_Coeffs(u, du(1:nbDim), 1, 2, 1, 1, Re_1, 0) * d2u(1,2), &
               !Eval_Diffusion_Coeffs(u, du(1:nbDim), 2, 1, 1, 1, Re_1, 0) * d2u(2,1), &
               !Eval_Diffusion_Coeffs(u, du(1:nbDim), 2, 2, 1, 1, Re_1, 0) * d2u(2,2)
       endif


    elseif(ityp == 3) then
       ! output is the derivative of the right-hand-side with respect to x
       w(1) = utx - deltaux

    elseif(ityp == 4) then
       ! output is the derivative of the right-hand-side with respect to y
       w(1) = uty - deltauy

    elseif(ityp == 5) then
       ! output is the derivative of the time derivative of the exact sol. with respect to x
       w(1) = utx

    elseif(ityp == 6) then
       ! output is the derivative of the time derivative of the exact sol. with respect to y
       w(1) = uty

    elseif(ityp == 7) then
       ! output is the time derivative of the exact solution
       w(1) = ut

    elseif(ityp == 8) then
       ! output is the derivative of the exact solution with respect x
       w(1) = du(1)
    elseif(ityp == 9) then
       ! output is the derivative of the exact solution with respect x
       w(1) = du(2)

    else
       print*,'Unknown ityp = ', ityp,' in Set_Model_Data'
       stop
    endif

    deallocate(d2u, du)

    if(ndim > 1) then
       do i=2, ndim
          w(i) = w(1)
       enddo
    endif

    !write(*,*) 'WEDESWES', ityp, imod, w(:)

  end subroutine Set_Model_Data



  !> compute Vectors \f$ f_s = {\bf f_s}({\bf w}) ,\ s=1,2\f$
  !> for scalar equation
  subroutine Set_f_s_scalar(ndimL, nbDim, Qdof, w, f_s, x, ie )
    integer, intent(in) :: Qdof, ndimL, nbDim
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
    real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
    integer, intent(in) :: ie
    integer, parameter :: ityp = 0   ! function f_s
    integer :: i, k, l

    do l=1,ndimL
       do i=1,Qdof
          do k=1,nbDim
             f_s(i, k, l) = Eval_Convective_Coeffs( w(i,l), k, 1, ityp, x(i, 1:nbDim)  )
          enddo
       enddo
    enddo

  end subroutine Set_f_s_scalar

  !> compute Vectors \f$ {\bf f}'_s = \frac{d}{d {\bf w}}{\bf f_s}({\bf w}) ,\ s=1,2\f$
  !> for scalar equation
  subroutine Der_f_s_scalar(ndimL, Qdof, w, Df_s, x )
    integer, intent(in) :: Qdof, ndimL
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:nbDim, 1:ndimL, 1:ndimL), intent(inout) :: Df_s
    integer, parameter :: ityp = 1   ! derivative f_s'(u)
    real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
    integer :: i,k,l

    Df_s(:,:,:,:) = 0.

    do l=1,ndimL
       do i=1,Qdof
          do k=1,nbDim
             ! diagonal blocks !!!
             Df_s(i, k, l,l) = Eval_Convective_Coeffs( w(i,l), k, 1, ityp, x(i, 1:nbDim) )
          enddo
       enddo
    enddo

  end subroutine Der_f_s_scalar

  !> compute Vectors on numerical flux
  !> \f$ H(w1, w2, n) \approx \sum_{s=1}^2{\bf f_s}({\bf w})n_s\f$ in integ nodes
  !> for scalar equation
  subroutine Set_NumFlux_scalar(ndimL, nbDim, Qdof, wi, wj, nc, xi, H, area_1, ie )
    integer, intent(in) :: Qdof, ndimL, nbDim, ie
    real, dimension(1:Qdof, 1:ndimL), intent(in):: wi, wj ! state  w in integ nodes
    real, dimension(1:Qdof, 1:nbDim), intent(in):: nc        ! outer normal in integ nodes
    real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
    real, dimension(1:Qdof,1:ndimL), intent(inout) :: H   ! numer flux H in  -- " --
    real, intent(in) :: area_1
    real, dimension(:,:), allocatable :: wa
    real, dimension(:,:,:), allocatable ::  f_s
    real, dimension(:,:,:,:), allocatable ::  Df_s
    integer :: l,k
    real :: val

    allocate(wa(1:Qdof,1:ndimL) )
    allocate(f_s(1:Qdof,1:nbDim, 1:ndimL) )
    allocate(Df_s(1:Qdof,1:nbDim, 1:ndimL, 1:ndimL) )

    wa(1:Qdof,1:ndimL) = (wi(1:Qdof,1:ndimL) + wj(1:Qdof,1:ndimL) )/2
    ! derivatives
    call Der_f_s_scalar(ndimL, Qdof, wa(1:Qdof,1:ndimL), &
         Df_s(1:Qdof,1:nbDim, 1:ndimL, 1:ndimL), xi(1:Qdof,1 :nbDim))

    do k=1,ndimL
       do l=1,Qdof
          val = Df_s(l,1,k,k) * nc(l, 1)  + Df_s(l,2,k,k) * nc(l, 2)
          if(val > 0) then
             wa(l,1:ndim)=  wi(l, 1:ndimL)
          else
             wa(l,1:ndim)=  wj(l, 1:ndimL)
          endif
          !!elem%max_eigenvals = max(elem%max_eigenvals, abs(val) * area_1 )
          state%max_eigenvals = max(state%max_eigenvals, abs(val) * area_1 )
       enddo
    enddo

    call Set_f_s_scalar(ndim, nbDim, Qdof, wa, f_s(1:Qdof,1:nbDim, 1:ndimL), xi, ie )

    write(*,'(a6,16es12.4)') 'f 1:',f_s(1:Qdof, 1, :)
    write(*,'(a6,16es12.4)') 'f 2:',f_s(1:Qdof, 2, :)

    do l=1,ndimL
       H(1:Qdof,l) = f_s(1:Qdof, 1, l)*nc(1:Qdof,1) + f_s(1:Qdof, 2, l)*nc(1:Qdof,2)
    enddo

    deallocate(wa, f_s, Df_s)

  end subroutine Set_NumFlux_scalar

  !> compute matrices \f$ A_s = \frac{D{\bf f_s}({\bf w})}{D{\bf w}},\ s=1,2\f$
  !> for scalar equation
  subroutine Set_A_s_scalar(ndimL, nbDim, Qdof, w, A_s, xi, ie)
    integer, intent(in) :: Qdof, nbdim, ndimL
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: A_s
    ! matrices A_s in  -- " --
    real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
    integer, intent(in) :: ie

    call Der_f_s_scalar(ndimL, Qdof, w, A_s(1:Qdof, 1:nbDim, 1:ndimL, 1:ndimL ), &
         xi(1:Qdof, 1:nbDim) )

  end subroutine Set_A_s_scalar

  !> compute matrices
  !> \f$ P^{\pm} = \left(\frac{D({\bf f_1}({\bf w})n_1+{\bf f_2}({\bf w})n_2}{D{\bf w}}
  !>  \right)^{\pm}\f$
  !> for scalar equation
  subroutine Set_Ppm_scalar(ndimL, nbDim, Qdof, w, n, xi, Ppm, one_over_area, elem)
    integer, intent(in) :: Qdof, ndimL, nbDim
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: Ppm
    ! matrices Ppm in  -- " --
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: n   ! outer normal
    real, dimension(1:Qdof, 1:nbDim),intent(in) ::  xi          ! node on the edge?
    real, intent(in), optional :: one_over_area !
    type( element ), intent( inout ), optional :: elem !not used

    real, dimension(:,:,:,:), allocatable ::  Df_s
    real :: val, rat
    integer :: l,k

    rat =  state%model%conv_rat

    !rat = 1.
    !print *,'!!!only for Burgers !!!  IMOD  !!!!!!!!!', state%max_eigenvals
    !rat = 2.
    !rat = 3.

    allocate(Df_s(1:Qdof,1:nbDim, 1:ndimL, 1:ndimL) )
    call Der_f_s_scalar(ndimL, Qdof, w(1:Qdof, 1:ndimL), &
         Df_s(1:Qdof,1:nbDim, 1:ndimL, 1:ndimL), xi(1:Qdof, 1:nbDim) )

    Ppm(1:Qdof, 1:nbDim, 1:ndimL, 1:ndimL) = 0.

    elem%max_eigenvals = 0.
    do k=1,ndimL
       do l=1,Qdof
          val = df_s(l,1,k,k) * n(l, 1)  + df_s(l,2,k,k) * n(l, 2)
          if(val > 0) then

             Ppm(l, 1, k, k) = val/rat
          else
             Ppm(l, 2, k, k) = val/rat
          endif

          !write(*,'(a4,i5,30es10.2)') '####', l, val,rat, Ppm(l,1,1:ndim, 1:ndim), &
          !  Ppm(l,2,1:ndim, 1:ndim), w(l, 1:ndim), df_s(l,1,:), df_s(l,2,:)
          elem%max_eigenvals = max(elem%max_eigenvals, abs(val) * one_over_area )
          state%max_eigenvals = max(state%max_eigenvals, abs(val) * one_over_area )
       enddo
    enddo


    deallocate(Df_s)


  end subroutine Set_Ppm_scalar


  !> compute viscous fluxes R_s, s=1,2 for scalar equation
  !> in integ nodes
  subroutine Set_R_s_scalar(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
    integer, intent(in) :: ndimL, nbDim, iRe, Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi !physical cooedinates
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    !real, intent(in) :: Re_1                     ! inverse of Reynolds number
    real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
    integer :: i,j,k,l

    R_s(:, :, :) = 0.

    do l=1,ndimL
       do i=1,nbDim
          do j=1,nbDim
             do k=1,Qdof

                !print*, iRe, k , size(Re_1(:,1)) , size(Re_1(1,:))
!                print*, 'RE_1:', Re_1
                R_s(k, i, l) = R_s(k, i, l)   &
                     + Dw(k, 1, j)* &
                     Eval_Diffusion_Coeffs(w(k,l), Dw(k,l,1:nbDim),i,j,l,l, &
                     Re_1(1:iRe, k), 0, xi(k, 1:2))
                !if(k==1) &
                     !write(*,'(a8, 4i4, es12.4)')
!                     print*, 'R_S:',k,i,l, j, R_s(k, i, l), Dw(k, 1, j), &
!                     Eval_Diffusion_Coeffs(w(k,l), Dw(k,l,1:nbDim),i,j,l,l, Re_1(k), 0, xi(k, 1:2))
             enddo
          enddo
       enddo
    enddo

  end subroutine Set_R_s_scalar

  !> compute "matrices" 1x1 K_sk, s,k=1,2 for scalar equations in integ nodes
  subroutine Set_K_sk_scalar(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, K_sk, xi)
    integer, intent(in) :: ndimL, nbDim, iRe, Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi !physical cooedinates
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,ndimL), intent(inout) :: K_sk
    integer :: i,j,k,l

    K_sk(:, :, :, : ,:) = 0.

    do l=1,ndimL
       do i=1,nbDim
          do j=1,nbDim
             do k=1,Qdof
                K_sk(k, i, j, l,l) = &
                     Eval_Diffusion_Coeffs(w(k,l), Dw(k,l,1:nbDim),i, j,l,l,Re_1(1:iRe, k),&
                     0, xi(k, 1:nbDim) )
             enddo
          enddo
       enddo
    enddo


  end subroutine Set_K_sk_scalar


  !> compute reactive terms S,  for scalar equation in integ nodes
  subroutine Set_S_scalar(ndimL, nbDim, Qdof, xi, w, Dw, S)
    integer, intent(in) :: ndimL, nbDim, Qdof
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL), intent(inout) :: S
    integer :: k,l

    do l=1,ndimL
       do k=1,Qdof
          S(k, l) = Eval_Reaction_Coeffs(w(k,l), Dw(k,l,1:nbDim), xi(k, 1:nbDim), 0 )
       enddo
    enddo

  end subroutine Set_S_scalar

  !> compute derivative of the reactive terms S,  for scalar equation in integ nodes
  subroutine Set_DS_scalar(ndimL, nbDim, Qdof, xi, w, Dw, DS)
    integer, intent(in) :: ndimL, nbDim, Qdof
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:ndimL), intent(inout) :: DS
    integer :: k, l

    DS(:, :, :) = 0.

    do l=1, ndimL
       do k=1,Qdof
          DS(k, l, l) = Eval_Reaction_Coeffs(w(k,l), Dw(k,l,1:nbDim), xi(k, 1:nbDim), 1 )
       enddo
    enddo

  end subroutine Set_DS_scalar

  !> evaluate exact \f$wi\f$ at node \f$x\f$ for scalar equation
  subroutine Exact_Scalar(x, wi, t)
    real, dimension(1:nbDim), intent(in) :: x
    real, dimension(1:ndim), intent(out) :: wi
    real, intent(in) :: t

    call Set_Model_Data(x, t, wi, 1)

 end subroutine Exact_Scalar

 !> evaluate the space derivation of the exact solution
 !> \f$Dwi\f$ at node \f$x\f$ for scalar equation
 subroutine Der_Exact_Scalar(x, Dwi, t)
   real, dimension(1:nbDim), intent(in) :: x
   real, dimension(1:ndim, 1:2), intent(out) :: Dwi
   real, intent(in) :: t

   call Set_Model_Data(x, t, Dwi(:,1), 8)
   call Set_Model_Data(x, t, Dwi(:,2), 9)


 end subroutine Der_Exact_Scalar

 !> evaluate RHS  \f$ f \f$ at node \f$x\f$ for scalar equation
 subroutine RHS_Scalar(x, f, t)
   real, dimension(1:nbDim), intent(in) :: x
   real, dimension(1:ndim), intent(out) :: f
   real, intent(in) :: t

   call Set_Model_Data(x, t, f, 2)

 end subroutine RHS_Scalar


 !> evaluate the gradient of RHS  \f$ f \f$ at node \f$x\f$ for scalar equation
 subroutine Der_RHS_Scalar(x, Nablaf, t)
   real, dimension(1:nbDim), intent(in) :: x
   real, dimension(1:ndim, 1:nbDim), intent(out) :: Nablaf
   real, intent(in) :: t

   call Set_Model_Data(x, t, Nablaf(1:ndim, 1), 3)
   call Set_Model_Data(x, t, Nablaf(1:ndim, 2), 4)

 end subroutine Der_RHS_Scalar

 !> evaluate the gradient of the time derivative of the difference of \f$ {\bf w}_h\f$
 !> and the exact solution
 subroutine Der_TimeDer_Exact_Scalar(x, utgrad, t)
   real, dimension(1:2), intent(in) :: x
   real, dimension(1:ndim, 1:2), intent(out) :: utgrad
   real, intent(in) :: t

   call Set_Model_Data(x, t, utgrad(1:ndim, 1), 5)
   call Set_Model_Data(x, t, utgrad(1:ndim, 2), 6)

 end subroutine Der_TimeDer_Exact_Scalar

 !> evaluate RHS  \f$ f \f$ at node \f$x\f$ for scalar equation
 subroutine TimeDer_Exact_Scalar( x, f, t)
   real, dimension(1:nbDim), intent(in) :: x
   real, dimension(1:ndim), intent(out) :: f
   real, intent(in) :: t

   call Set_Model_Data(x, t, f, 7)

 end subroutine TimeDer_Exact_Scalar




  subroutine Detect_apriori_known_singularity(elem, singularity)
    type(element), intent(in):: elem      ! elem = element
    logical, intent(inout) :: singularity
    real:: beta, alpha, xx
    integer :: i, j, ip, im

    singularity = .false.

    select case(state%model%icase)

    case(40)  ! L-shape domain, singularity at the origin

       if(dot_product(grid%x(elem%face(idx,1),1:2),grid%x(elem%face(idx,1),1:2)) > 1E-15 .and.&
            dot_product(grid%x(elem%face(idx,2),1:2), grid%x(elem%face(idx,2),1:2))>1E-15.and.&
            dot_product(grid%x(elem%face(idx,3),1:2),grid%x(elem%face(idx,3),1:2))>1E-15 ) then
          singularity = .false.
       else
          singularity = .true.
       endif


    case(56, 58)  !interior line singularity x - beta*(y - 1 )
       beta = 0.6
       ip = 0
       im = 0
       do j=1, elem%flen
           xx = grid%x(elem%face(idx, j), 1) - beta*(grid%x(elem%face(idx, j), 2) + 1)
           if(xx > 0) ip = ip + 1
           if(xx < 0) im = im + 1
        enddo

        if(ip > 0 .and. im > 0) singularity = .true. ! triangle has vertices up and bellow line

     case default
        if(elem%i == 1) print *, 'Singularity not detected'

    end select

    !if(singularity ) then
    !   print*,'Singul :', elem%xc
    !   write(91, *) elem%xc
    !endif


  end subroutine Detect_apriori_known_singularity


end module model_oper
