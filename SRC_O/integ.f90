!> volume (triangular, quadrilateral) and edge quadrature rules
!>
!> contain also the values of the basis functions (phi) and their derivatives
!> (Dphi) in the corresponding integration nodes,

module integration
  use lapack_oper
  use quadrature_mod
  use tquadrature_mod

  implicit none
  public

!  !>  quadratic rule for volume integrals
!  type, public :: volume_rule
!     integer::  Qdeg    !  order of approximation
!     integer::  Qdof    !  number of integration nodes
!     real, pointer, dimension(:) :: weights
!     real, pointer, dimension(:,:) :: lambda
!     real, dimension(:,:), pointer :: phi
!     real, dimension(:,:,:), pointer :: Dphi
!     logical :: def
!
!     contains
!
!     procedure :: createVrule
!
!  end type volume_rule

!  !>  quadratic rule for edge integrals
!  type :: Gauss_rule
!     integer::  Qdeg         !  type of a rule, accurate for 2*Qdeg-1 polynom
!     integer::  Qdof         !  number of integration nodes
!     real, pointer, dimension(:) :: weights
!     real, pointer, dimension(:) :: lambda
!     real, dimension(:,:,:,:,:), pointer :: phi
!     real, dimension(:,:,:,:,:,:), pointer :: Dphi
!     real, dimension(:,:), pointer :: Leg_phi      ! Legendre polynomials
!
!  end type Gauss_rule

!  !> Lagrangian nodes for visualisation
!  type :: Lagrang_rule         !
!     integer::  Qnum    !
!     integer::  Qdeg    !  order of approximation
!     integer::  Qdof    !  number of integration nodes
!     real, pointer, dimension(:,:) :: lambda
!     real, dimension(:,:), pointer :: phi   ! DGFE basis in Lagrangian integ nodes
!     real, dimension(:,:), pointer :: psi   ! coefficients of Lagrangian basis
!     integer, dimension(:,:), pointer :: psi_pow  ! powers of x_1 and x_2 in in Lagrangian basis
!  end type Lagrang_rule

!  !>  rule for the time intervals
!  type :: Time_rule
!     integer::  Qdeg         !  type of a rule, accurate for 2*Qdeg-1 polynom / resp 2*Qdeg - 2 for Radau
!     integer::  Qdof         !  number of integration nodes
!     real, pointer, dimension(:) :: weights
!     real, pointer, dimension(:) :: lambda
!     real, dimension(:,:), pointer :: phi
!     real, dimension(:,:), pointer :: Dphi
!
!  end type Time_rule



  !> Adaptive Backward Differential Formulae (\f$n\f$-steps, \f$ n = \f$ deg ), pamaters
  !> \f$ \alpha_i,\ i=0,\dots, n\f$,   \f$ \gamma\f$, \f$ \delta_i,\ i=0,\dots, n+1\f$,
  !>
  !> \f$  y' = F(y),  y^k \approx y(t_k)\f$,
  !> \f$ \sum_{i=0}^{n} \alpha_i y^{k-i} = \tau_k F(y^k) \f$,
  !> local discretization error \f$ L_k \f$:
  !> \f$ \sum_{i=0}^{n} \alpha_i y^{k-i} = \tau_k F(y^k) + L_k\f$,
  !> \f$ L_k \approx \gamma \tau_k^{n+1} y^{(n+1)} \f$,
  !> for the aproximation of the  \f$ y^{(n+1)} \f$ we use backward differences, i.e.
  !> \f$ y^{(n+1)} \approx \frac{1}{\tau_k^{n+1}} \sum_{i=0}^{n+1} \delta_i y^{k-i} \f$
  !> for semi-implicit method, extrapolation:
  !> \f$ y^{k} \approx  \sum_{i=1}^{n} \beta_i y^{k-i}\f$,
  !> \f$ \beta_i \f$ = BDF%extrap(i),
  !> for error estimate via explicit extrapolation:
  !> \f$ y^{k} \approx  \sum_{i=1}^{n+1} \bar{\beta_i} y^{k-i}\f$,
  !> \f$ \bar{\beta}_i \f$ = BDF%Bextrap(i)
  type :: ABDF
     integer :: max_deg      ! maximal implemented degree
     integer :: deg          ! degree of this scheme
     !integer :: deg2         ! degree of the second scheme for the error estim.
     real, pointer, dimension(:) :: alpha
     real, pointer, dimension(:) :: Bextrap
     real, pointer, dimension(:) :: extrap
     real, pointer, dimension(:) :: delta
     real  :: gamm
  end type ABDF

  public :: Create_volume_rule
  public :: Create_Gauss_rule
!  public :: Create_time_rule
!  public :: Create_Radau_rule
  public :: Create_L_rule
  public :: Eval_L_rule
  public :: Eval_L_rule1
  public :: Eval_L_Direct
  public :: Add_4volume_rule
  public :: AllocateABDF
!  public :: SetABDF
!  public :: FakeBDFLagrangeExtrapolation

! moved to paramets
!  public:: QnumOffset
!  !integer,parameter:: QnumOffset = 21
!  integer :: QnumOffset = maxVrule
contains

  !> generates one edge (Gauss) quadrature rule
  subroutine Create_Gauss_rule(G_rule, Qdeg)
    type(Gauss_rule), intent(inout) :: G_rule
    integer, intent(in) :: Qdeg

    G_rule%Qdeg = max(Qdeg, 1)
    G_rule%Qdof = G_rule%Qdeg

    allocate(G_rule%weights(1:G_rule%Qdeg))
    allocate(G_rule%lambda(1:G_rule%Qdeg))

    select case (Qdeg)
    case(:1)
       G_rule%weights(1) = 1.
       G_rule%lambda(1) = 0.5
    case(2)
       G_rule%weights(1:2) = 0.5
       G_rule%lambda(1) = (1.-0.57735026918962576451)/2
       G_rule%lambda(2) = 1.57735026918962576451/2
    case(3)
       G_rule%weights(1) = 5./18
       G_rule%weights(2) = 8./18
       G_rule%weights(3) = 5./18
       G_rule%lambda(1) = (1.-0.77459666924148337704)/2
       G_rule%lambda(2) = 0.5
       G_rule%lambda(3) = 1.77459666924148337704/2
    case(4)
       G_rule%weights(1) = 0.34785484513745385737/2
       G_rule%weights(2) = 0.65214515486254614263/2
       G_rule%weights(3) = G_rule%weights(2)
       G_rule%weights(4) = G_rule%weights(1)
       G_rule%lambda(1) = (1-0.86113631159405257522)/2
       G_rule%lambda(2) = (1-0.33998104358485626480)/2
       G_rule%lambda(3) = 1.33998104358485626480/2
       G_rule%lambda(4) = 1.86113631159405257522/2
    case(5)
       G_rule%weights(1) = 0.23692688505618908751/2.
       G_rule%weights(2) = 0.47862867049936646804/2.
       G_rule%weights(3) = 0.56888888888888888889/2.
       G_rule%weights(4) = G_rule%weights(2)
       G_rule%weights(5) = G_rule%weights(1)

       G_rule%lambda(1) = (1-0.90617984593866399280)/2.
       G_rule%lambda(2) = (1-0.53846931010568309104)/2.
       G_rule%lambda(3) = 0.5
       G_rule%lambda(4) = 1.53846931010568309104/2.
       G_rule%lambda(5) = 1.90617984593866399280/2.
    case(6)
       G_rule%weights(1) = 0.17132449237917034504/2.
       G_rule%weights(2) = 0.36076157304813860757/2.
       G_rule%weights(3) = 0.46791393457269104739/2.
       G_rule%weights(4) = G_rule%weights(3)
       G_rule%weights(5) = G_rule%weights(2)
       G_rule%weights(6) = G_rule%weights(1)

       G_rule%lambda(1) = (1.-0.93246951420315202781)/2.
       G_rule%lambda(2) = (1.-0.66120938646626451366)/2.
       G_rule%lambda(3) = (1.-0.23861918608319690863)/2.
       G_rule%lambda(4) =  1.23861918608319690863/2.
       G_rule%lambda(5) =  1.66120938646626451366/2.
       G_rule%lambda(6) =  1.93246951420315202781/2.
     case(7)
       G_rule%weights(1) = 0.12948496616886969327/2.
       G_rule%weights(2) = 0.27970539148927666790/2.
       G_rule%weights(3) = 0.38183005050511894495/2.
       G_rule%weights(4) = 0.41795918367346938776/2.
       G_rule%weights(5) = G_rule%weights(3)
       G_rule%weights(6) = G_rule%weights(2)
       G_rule%weights(7) = G_rule%weights(1)

       G_rule%lambda(1) = (1.-0.94910791234275852453)/2.
       G_rule%lambda(2) = (1.-0.74153118559939443986)/2.
       G_rule%lambda(3) = (1.-0.40584515137739716691)/2.
       G_rule%lambda(4) =  0.5
       G_rule%lambda(5) =  1.40584515137739716691/2.
       G_rule%lambda(6) =  1.74153118559939443986/2.
       G_rule%lambda(7) =  1.94910791234275852453/2.

      case( 8)
       G_rule%weights( 1) =  0.5061426814518818E-01
       G_rule%weights( 2) =  0.1111905172266872E+00
       G_rule%weights( 3) =  0.1568533229389436E+00
       G_rule%weights( 4) =  0.1813418916891809E+00
       G_rule%weights( 5) =  0.1813418916891809E+00
       G_rule%weights( 6) =  0.1568533229389436E+00
       G_rule%weights( 7) =  0.1111905172266872E+00
       G_rule%weights( 8) =  0.5061426814518818E-01

       G_rule%lambda( 1) =  0.1985507175123191E-01
       G_rule%lambda( 2) =  0.1016667612931866E+00
       G_rule%lambda( 3) =  0.2372337950418355E+00
       G_rule%lambda( 4) =  0.4082826787521751E+00
       G_rule%lambda( 5) =  0.5917173212478249E+00
       G_rule%lambda( 6) =  0.7627662049581645E+00
       G_rule%lambda( 7) =  0.8983332387068134E+00
       G_rule%lambda( 8) =  0.9801449282487681E+00

      case( 9)
       G_rule%weights( 1) =  0.4063719418078721E-01
       G_rule%weights( 2) =  0.9032408034742863E-01
       G_rule%weights( 3) =  0.1303053482014677E+00
       G_rule%weights( 4) =  0.1561735385200015E+00
       G_rule%weights( 5) =  0.1651196775006299E+00
       G_rule%weights( 6) =  0.1561735385200015E+00
       G_rule%weights( 7) =  0.1303053482014677E+00
       G_rule%weights( 8) =  0.9032408034742863E-01
       G_rule%weights( 9) =  0.4063719418078721E-01

       G_rule%lambda( 1) =  0.1591988024618696E-01
       G_rule%lambda( 2) =  0.8198444633668212E-01
       G_rule%lambda( 3) =  0.1933142836497048E+00
       G_rule%lambda( 4) =  0.3378732882980955E+00
       G_rule%lambda( 5) =  0.5000000000000000E+00
       G_rule%lambda( 6) =  0.6621267117019045E+00
       G_rule%lambda( 7) =  0.8066857163502952E+00
       G_rule%lambda( 8) =  0.9180155536633179E+00
       G_rule%lambda( 9) =  0.9840801197538130E+00

      case(10)
       G_rule%weights( 1) =  0.3333567215434402E-01
       G_rule%weights( 2) =  0.7472567457529031E-01
       G_rule%weights( 3) =  0.1095431812579911E+00
       G_rule%weights( 4) =  0.1346333596549982E+00
       G_rule%weights( 5) =  0.1477621123573765E+00
       G_rule%weights( 6) =  0.1477621123573765E+00
       G_rule%weights( 7) =  0.1346333596549982E+00
       G_rule%weights( 8) =  0.1095431812579911E+00
       G_rule%weights( 9) =  0.7472567457529031E-01
       G_rule%weights(10) =  0.3333567215434402E-01

       G_rule%lambda( 1) =  0.1304673574141413E-01
       G_rule%lambda( 2) =  0.6746831665550773E-01
       G_rule%lambda( 3) =  0.1602952158504878E+00
       G_rule%lambda( 4) =  0.2833023029353764E+00
       G_rule%lambda( 5) =  0.4255628305091844E+00
       G_rule%lambda( 6) =  0.5744371694908156E+00
       G_rule%lambda( 7) =  0.7166976970646236E+00
       G_rule%lambda( 8) =  0.8397047841495122E+00
       G_rule%lambda( 9) =  0.9325316833444923E+00
       G_rule%lambda(10) =  0.9869532642585859E+00

      case(11)
       G_rule%weights( 1) =  0.2783428355808687E-01
       G_rule%weights( 2) =  0.6279018473245233E-01
       G_rule%weights( 3) =  0.9314510546386705E-01
       G_rule%weights( 4) =  0.1165968822959953E+00
       G_rule%weights( 5) =  0.1314022722551234E+00
       G_rule%weights( 6) =  0.1364625433889503E+00
       G_rule%weights( 7) =  0.1314022722551234E+00
       G_rule%weights( 8) =  0.1165968822959953E+00
       G_rule%weights( 9) =  0.9314510546386705E-01
       G_rule%weights(10) =  0.6279018473245233E-01
       G_rule%weights(11) =  0.2783428355808687E-01

       G_rule%lambda( 1) =  0.1088567092697151E-01
       G_rule%lambda( 2) =  0.5646870011595234E-01
       G_rule%lambda( 3) =  0.1349239972129753E+00
       G_rule%lambda( 4) =  0.2404519353965941E+00
       G_rule%lambda( 5) =  0.3652284220238275E+00
       G_rule%lambda( 6) =  0.5000000000000000E+00
       G_rule%lambda( 7) =  0.6347715779761725E+00
       G_rule%lambda( 8) =  0.7595480646034058E+00
       G_rule%lambda( 9) =  0.8650760027870247E+00
       G_rule%lambda(10) =  0.9435312998840477E+00
       G_rule%lambda(11) =  0.9891143290730284E+00

      case(12)
       G_rule%weights( 1) =  0.2358766819325595E-01
       G_rule%weights( 2) =  0.5346966299765912E-01
       G_rule%weights( 3) =  0.8003916427167312E-01
       G_rule%weights( 4) =  0.1015837133615329E+00
       G_rule%weights( 5) =  0.1167462682691774E+00
       G_rule%weights( 6) =  0.1245735229067013E+00
       G_rule%weights( 7) =  0.1245735229067013E+00
       G_rule%weights( 8) =  0.1167462682691774E+00
       G_rule%weights( 9) =  0.1015837133615329E+00
       G_rule%weights(10) =  0.8003916427167312E-01
       G_rule%weights(11) =  0.5346966299765912E-01
       G_rule%weights(12) =  0.2358766819325595E-01

       G_rule%lambda( 1) =  0.9219682876640378E-02
       G_rule%lambda( 2) =  0.4794137181476255E-01
       G_rule%lambda( 3) =  0.1150486629028477E+00
       G_rule%lambda( 4) =  0.2063410228566913E+00
       G_rule%lambda( 5) =  0.3160842505009099E+00
       G_rule%lambda( 6) =  0.4373832957442655E+00
       G_rule%lambda( 7) =  0.5626167042557344E+00
       G_rule%lambda( 8) =  0.6839157494990901E+00
       G_rule%lambda( 9) =  0.7936589771433087E+00
       G_rule%lambda(10) =  0.8849513370971523E+00
       G_rule%lambda(11) =  0.9520586281852375E+00
       G_rule%lambda(12) =  0.9907803171233596E+00

      case(13)
       G_rule%weights( 1) =  0.2024200238265791E-01
       G_rule%weights( 2) =  0.4606074991886418E-01
       G_rule%weights( 3) =  0.6943675510989367E-01
       G_rule%weights( 4) =  0.8907299038097291E-01
       G_rule%weights( 5) =  0.1039080237684442E+00
       G_rule%weights( 6) =  0.1131415901314486E+00
       G_rule%weights( 7) =  0.1162757766154369E+00
       G_rule%weights( 8) =  0.1131415901314486E+00
       G_rule%weights( 9) =  0.1039080237684442E+00
       G_rule%weights(10) =  0.8907299038097291E-01
       G_rule%weights(11) =  0.6943675510989367E-01
       G_rule%weights(12) =  0.4606074991886418E-01
       G_rule%weights(13) =  0.2024200238265791E-01

       G_rule%lambda( 1) =  0.7908472640705932E-02
       G_rule%lambda( 2) =  0.4120080038851104E-01
       G_rule%lambda( 3) =  0.9921095463334506E-01
       G_rule%lambda( 4) =  0.1788253302798299E+00
       G_rule%lambda( 5) =  0.2757536244817765E+00
       G_rule%lambda( 6) =  0.3847708420224326E+00
       G_rule%lambda( 7) =  0.5000000000000000E+00
       G_rule%lambda( 8) =  0.6152291579775674E+00
       G_rule%lambda( 9) =  0.7242463755182235E+00
       G_rule%lambda(10) =  0.8211746697201701E+00
       G_rule%lambda(11) =  0.9007890453666549E+00
       G_rule%lambda(12) =  0.9587991996114890E+00
       G_rule%lambda(13) =  0.9920915273592941E+00

      case(14)
       G_rule%weights( 1) =  0.1755973016587599E-01
       G_rule%weights( 2) =  0.4007904357988007E-01
       G_rule%weights( 3) =  0.6075928534395157E-01
       G_rule%weights( 4) =  0.7860158357909681E-01
       G_rule%weights( 5) =  0.9276919873896890E-01
       G_rule%weights( 6) =  0.1025992318606479E+00
       G_rule%weights( 7) =  0.1076319267315789E+00
       G_rule%weights( 8) =  0.1076319267315789E+00
       G_rule%weights( 9) =  0.1025992318606479E+00
       G_rule%weights(10) =  0.9276919873896890E-01
       G_rule%weights(11) =  0.7860158357909681E-01
       G_rule%weights(12) =  0.6075928534395157E-01
       G_rule%weights(13) =  0.4007904357988007E-01
       G_rule%weights(14) =  0.1755973016587599E-01

       G_rule%lambda( 1) =  0.6858095651593843E-02
       G_rule%lambda( 2) =  0.3578255816821324E-01
       G_rule%lambda( 3) =  0.8639934246511749E-01
       G_rule%lambda( 4) =  0.1563535475941573E+00
       G_rule%lambda( 5) =  0.2423756818209230E+00
       G_rule%lambda( 6) =  0.3404438155360551E+00
       G_rule%lambda( 7) =  0.4459725256463282E+00
       G_rule%lambda( 8) =  0.5540274743536718E+00
       G_rule%lambda( 9) =  0.6595561844639448E+00
       G_rule%lambda(10) =  0.7576243181790770E+00
       G_rule%lambda(11) =  0.8436464524058427E+00
       G_rule%lambda(12) =  0.9136006575348825E+00
       G_rule%lambda(13) =  0.9642174418317868E+00
       G_rule%lambda(14) =  0.9931419043484062E+00

     case(15:)
       G_rule%weights( 1) =  0.1537662099805857E-01
       G_rule%weights( 2) =  0.3518302374405406E-01
       G_rule%weights( 3) =  0.5357961023358598E-01
       G_rule%weights( 4) =  0.6978533896307712E-01
       G_rule%weights( 5) =  0.8313460290849704E-01
       G_rule%weights( 6) =  0.9308050000778109E-01
       G_rule%weights( 7) =  0.9921574266355583E-01
       G_rule%weights( 8) =  0.1012891209627806E+00
       G_rule%weights( 9) =  0.9921574266355583E-01
       G_rule%weights(10) =  0.9308050000778109E-01
       G_rule%weights(11) =  0.8313460290849704E-01
       G_rule%weights(12) =  0.6978533896307712E-01
       G_rule%weights(13) =  0.5357961023358598E-01
       G_rule%weights(14) =  0.3518302374405406E-01
       G_rule%weights(15) =  0.1537662099805857E-01

       G_rule%lambda( 1) =  0.6003740989757256E-02
       G_rule%lambda( 2) =  0.3136330379964708E-01
       G_rule%lambda( 3) =  0.7589670829478640E-01
       G_rule%lambda( 4) =  0.1377911343199150E+00
       G_rule%lambda( 5) =  0.2145139136957306E+00
       G_rule%lambda( 6) =  0.3029243264612183E+00
       G_rule%lambda( 7) =  0.3994029530012828E+00
       G_rule%lambda( 8) =  0.5000000000000000E+00
       G_rule%lambda( 9) =  0.6005970469987173E+00
       G_rule%lambda(10) =  0.6970756735387817E+00
       G_rule%lambda(11) =  0.7854860863042694E+00
       G_rule%lambda(12) =  0.8622088656800850E+00
       G_rule%lambda(13) =  0.9241032917052137E+00
       G_rule%lambda(14) =  0.9686366962003530E+00
       G_rule%lambda(15) =  0.9939962590102427E+00

    end select

    if(Qdeg > 15) then
      print*,'Warning !! Edge quadrature rule of degree ', Qdeg, &
      ' is not implemented'
      print*,'           the maximal default quadrature rule of degree', &
      ' 15 is used'
    endif
  end subroutine Create_Gauss_rule

  !> generates one volume (Dunavant triangular) quadrature rule
  subroutine Create_volume_rule(V_rule, Qdeg)
    type(volume_rule), intent(inout) :: V_rule
    integer, intent(in) :: Qdeg

    V_rule%def = .true.

    V_rule%Qdeg = Qdeg
    select case (Qdeg)
    case( :1  )
       V_rule%Qdof =   1
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  1.0000000000000000E+00

       V_rule%lambda(  1,1:3) = (/ 3.3333333333333298E-01,  3.3333333333333298E-01,  3.3333333333333404E-01/)


    case(  2  )
       V_rule%Qdof =   3
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  3.3333333333333298E-01
       V_rule%weights(  2) =  3.3333333333333298E-01
       V_rule%weights(  3) =  3.3333333333333298E-01

       V_rule%lambda(  1,1:3) = (/ 6.6666666666666696E-01,  1.6666666666666699E-01,  1.6666666666666605E-01/)
       V_rule%lambda(  2,1:3) = (/ 1.6666666666666699E-01,  1.6666666666666699E-01,  6.6666666666666607E-01/)
       V_rule%lambda(  3,1:3) = (/ 1.6666666666666699E-01,  6.6666666666666696E-01,  1.6666666666666605E-01/)


    case(  3  )
       V_rule%Qdof =   4
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) = -5.6250000000000000E-01
       V_rule%weights(  2) =  5.2083333333333304E-01
       V_rule%weights(  3) =  5.2083333333333304E-01
       V_rule%weights(  4) =  5.2083333333333304E-01

       V_rule%lambda(  1,1:3) = (/ 3.3333333333333298E-01,  3.3333333333333298E-01,  3.3333333333333404E-01/)
       V_rule%lambda(  2,1:3) = (/ 5.9999999999999998E-01,  2.0000000000000001E-01,  2.0000000000000001E-01/)
       V_rule%lambda(  3,1:3) = (/ 2.0000000000000001E-01,  2.0000000000000001E-01,  5.9999999999999998E-01/)
       V_rule%lambda(  4,1:3) = (/ 2.0000000000000001E-01,  5.9999999999999998E-01,  2.0000000000000001E-01/)


    case(  4  )
       V_rule%Qdof =   6
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  2.2338158967801100E-01
       V_rule%weights(  2) =  2.2338158967801100E-01
       V_rule%weights(  3) =  2.2338158967801100E-01
       V_rule%weights(  4) =  1.0995174365532200E-01
       V_rule%weights(  5) =  1.0995174365532200E-01
       V_rule%weights(  6) =  1.0995174365532200E-01

       V_rule%lambda(  1,1:3) = (/ 1.0810301816807000E-01,  4.4594849091596500E-01,  4.4594849091596500E-01/)
       V_rule%lambda(  2,1:3) = (/ 4.4594849091596500E-01,  4.4594849091596500E-01,  1.0810301816807000E-01/)
       V_rule%lambda(  3,1:3) = (/ 4.4594849091596500E-01,  1.0810301816807000E-01,  4.4594849091596500E-01/)
       V_rule%lambda(  4,1:3) = (/ 8.1684757298045896E-01,  9.1576213509771007E-02,  9.1576213509770035E-02/)
       V_rule%lambda(  5,1:3) = (/ 9.1576213509771007E-02,  9.1576213509771007E-02,  8.1684757298045796E-01/)
       V_rule%lambda(  6,1:3) = (/ 9.1576213509771007E-02,  8.1684757298045896E-01,  9.1576213509770035E-02/)


    case(  5  )
       V_rule%Qdof =   7
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  2.2500000000000001E-01
       V_rule%weights(  2) =  1.3239415278850600E-01
       V_rule%weights(  3) =  1.3239415278850600E-01
       V_rule%weights(  4) =  1.3239415278850600E-01
       V_rule%weights(  5) =  1.2593918054482700E-01
       V_rule%weights(  6) =  1.2593918054482700E-01
       V_rule%weights(  7) =  1.2593918054482700E-01

       V_rule%lambda(  1,1:3) = (/ 3.3333333333333298E-01,  3.3333333333333298E-01,  3.3333333333333404E-01/)
       V_rule%lambda(  2,1:3) = (/ 5.9715871789770003E-02,  4.7014206410511500E-01,  4.7014206410511500E-01/)
       V_rule%lambda(  3,1:3) = (/ 4.7014206410511500E-01,  4.7014206410511500E-01,  5.9715871789770003E-02/)
       V_rule%lambda(  4,1:3) = (/ 4.7014206410511500E-01,  5.9715871789770003E-02,  4.7014206410511500E-01/)
       V_rule%lambda(  5,1:3) = (/ 7.9742698535308698E-01,  1.0128650732345600E-01,  1.0128650732345702E-01/)
       V_rule%lambda(  6,1:3) = (/ 1.0128650732345600E-01,  1.0128650732345600E-01,  7.9742698535308798E-01/)
       V_rule%lambda(  7,1:3) = (/ 1.0128650732345600E-01,  7.9742698535308698E-01,  1.0128650732345702E-01/)


    case(  6  )
       V_rule%Qdof =  12
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  1.1678627572637899E-01
       V_rule%weights(  2) =  1.1678627572637899E-01
       V_rule%weights(  3) =  1.1678627572637899E-01
       V_rule%weights(  4) =  5.0844906370206999E-02
       V_rule%weights(  5) =  5.0844906370206999E-02
       V_rule%weights(  6) =  5.0844906370206999E-02
       V_rule%weights(  7) =  8.2851075618374001E-02
       V_rule%weights(  8) =  8.2851075618374001E-02
       V_rule%weights(  9) =  8.2851075618374001E-02
       V_rule%weights( 10) =  8.2851075618374001E-02
       V_rule%weights( 11) =  8.2851075618374001E-02
       V_rule%weights( 12) =  8.2851075618374001E-02

       V_rule%lambda(  1,1:3) = (/ 5.0142650965817903E-01,  2.4928674517091001E-01,  2.4928674517091096E-01/)
       V_rule%lambda(  2,1:3) = (/ 2.4928674517091001E-01,  2.4928674517091001E-01,  5.0142650965818003E-01/)
       V_rule%lambda(  3,1:3) = (/ 2.4928674517091001E-01,  5.0142650965817903E-01,  2.4928674517091096E-01/)
       V_rule%lambda(  4,1:3) = (/ 8.7382197101699599E-01,  6.3089014491502005E-02,  6.3089014491502005E-02/)
       V_rule%lambda(  5,1:3) = (/ 6.3089014491502005E-02,  6.3089014491502005E-02,  8.7382197101699599E-01/)
       V_rule%lambda(  6,1:3) = (/ 6.3089014491502005E-02,  8.7382197101699599E-01,  6.3089014491502005E-02/)
       V_rule%lambda(  7,1:3) = (/ 5.3145049844817001E-02,  3.1035245103378400E-01,  6.3650249912139900E-01/)
       V_rule%lambda(  8,1:3) = (/ 3.1035245103378400E-01,  6.3650249912139900E-01,  5.3145049844816994E-02/)
       V_rule%lambda(  9,1:3) = (/ 6.3650249912139900E-01,  5.3145049844817001E-02,  3.1035245103378400E-01/)
       V_rule%lambda( 10,1:3) = (/ 3.1035245103378400E-01,  5.3145049844817001E-02,  6.3650249912139900E-01/)
       V_rule%lambda( 11,1:3) = (/ 6.3650249912139900E-01,  3.1035245103378400E-01,  5.3145049844816994E-02/)
       V_rule%lambda( 12,1:3) = (/ 5.3145049844817001E-02,  6.3650249912139900E-01,  3.1035245103378400E-01/)


    case(  7  )
       V_rule%Qdof =  13
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) = -1.4957004446768199E-01
       V_rule%weights(  2) =  1.7561525743320799E-01
       V_rule%weights(  3) =  1.7561525743320799E-01
       V_rule%weights(  4) =  1.7561525743320799E-01
       V_rule%weights(  5) =  5.3347235608838001E-02
       V_rule%weights(  6) =  5.3347235608838001E-02
       V_rule%weights(  7) =  5.3347235608838001E-02
       V_rule%weights(  8) =  7.7113760890256997E-02
       V_rule%weights(  9) =  7.7113760890256997E-02
       V_rule%weights( 10) =  7.7113760890256997E-02
       V_rule%weights( 11) =  7.7113760890256997E-02
       V_rule%weights( 12) =  7.7113760890256997E-02
       V_rule%weights( 13) =  7.7113760890256997E-02

       V_rule%lambda(  1,1:3) = (/ 3.3333333333333298E-01,  3.3333333333333298E-01,  3.3333333333333404E-01/)
       V_rule%lambda(  2,1:3) = (/ 4.7930806784191998E-01,  2.6034596607903998E-01,  2.6034596607904004E-01/)
       V_rule%lambda(  3,1:3) = (/ 2.6034596607903998E-01,  2.6034596607903998E-01,  4.7930806784192004E-01/)
       V_rule%lambda(  4,1:3) = (/ 2.6034596607903998E-01,  4.7930806784191998E-01,  2.6034596607904004E-01/)
       V_rule%lambda(  5,1:3) = (/ 8.6973979419556802E-01,  6.5130102902216006E-02,  6.5130102902215978E-02/)
       V_rule%lambda(  6,1:3) = (/ 6.5130102902216006E-02,  6.5130102902216006E-02,  8.6973979419556802E-01/)
       V_rule%lambda(  7,1:3) = (/ 6.5130102902216006E-02,  8.6973979419556802E-01,  6.5130102902215978E-02/)
       V_rule%lambda(  8,1:3) = (/ 4.8690315425316003E-02,  3.1286549600487401E-01,  6.3844418856981000E-01/)
       V_rule%lambda(  9,1:3) = (/ 3.1286549600487401E-01,  6.3844418856981000E-01,  4.8690315425315989E-02/)
       V_rule%lambda( 10,1:3) = (/ 6.3844418856981000E-01,  4.8690315425316003E-02,  3.1286549600487401E-01/)
       V_rule%lambda( 11,1:3) = (/ 3.1286549600487401E-01,  4.8690315425316003E-02,  6.3844418856981000E-01/)
       V_rule%lambda( 12,1:3) = (/ 6.3844418856981000E-01,  3.1286549600487401E-01,  4.8690315425315989E-02/)
       V_rule%lambda( 13,1:3) = (/ 4.8690315425316003E-02,  6.3844418856981000E-01,  3.1286549600487401E-01/)


    case(  8  )
       V_rule%Qdof =  16
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  1.4431560767778701E-01
       V_rule%weights(  2) =  9.5091634267284994E-02
       V_rule%weights(  3) =  9.5091634267284994E-02
       V_rule%weights(  4) =  9.5091634267284994E-02
       V_rule%weights(  5) =  1.0321737053471799E-01
       V_rule%weights(  6) =  1.0321737053471799E-01
       V_rule%weights(  7) =  1.0321737053471799E-01
       V_rule%weights(  8) =  3.2458497623198003E-02
       V_rule%weights(  9) =  3.2458497623198003E-02
       V_rule%weights( 10) =  3.2458497623198003E-02
       V_rule%weights( 11) =  2.7230314174435000E-02
       V_rule%weights( 12) =  2.7230314174435000E-02
       V_rule%weights( 13) =  2.7230314174435000E-02
       V_rule%weights( 14) =  2.7230314174435000E-02
       V_rule%weights( 15) =  2.7230314174435000E-02
       V_rule%weights( 16) =  2.7230314174435000E-02

       V_rule%lambda(  1,1:3) = (/ 3.3333333333333298E-01,  3.3333333333333298E-01,  3.3333333333333404E-01/)
       V_rule%lambda(  2,1:3) = (/ 8.1414823414554000E-02,  4.5929258829272301E-01,  4.5929258829272301E-01/)
       V_rule%lambda(  3,1:3) = (/ 4.5929258829272301E-01,  4.5929258829272301E-01,  8.1414823414553972E-02/)
       V_rule%lambda(  4,1:3) = (/ 4.5929258829272301E-01,  8.1414823414554000E-02,  4.5929258829272301E-01/)
       V_rule%lambda(  5,1:3) = (/ 6.5886138449648002E-01,  1.7056930775175999E-01,  1.7056930775175999E-01/)
       V_rule%lambda(  6,1:3) = (/ 1.7056930775175999E-01,  1.7056930775175999E-01,  6.5886138449648002E-01/)
       V_rule%lambda(  7,1:3) = (/ 1.7056930775175999E-01,  6.5886138449648002E-01,  1.7056930775175999E-01/)
       V_rule%lambda(  8,1:3) = (/ 8.9890554336593798E-01,  5.0547228317030998E-02,  5.0547228317031026E-02/)
       V_rule%lambda(  9,1:3) = (/ 5.0547228317030998E-02,  5.0547228317030998E-02,  8.9890554336593798E-01/)
       V_rule%lambda( 10,1:3) = (/ 5.0547228317030998E-02,  8.9890554336593798E-01,  5.0547228317031026E-02/)
       V_rule%lambda( 11,1:3) = (/ 8.3947774099580007E-03,  2.6311282963463800E-01,  7.2849239295540402E-01/)
       V_rule%lambda( 12,1:3) = (/ 2.6311282963463800E-01,  7.2849239295540402E-01,  8.3947774099579764E-03/)
       V_rule%lambda( 13,1:3) = (/ 7.2849239295540402E-01,  8.3947774099580007E-03,  2.6311282963463800E-01/)
       V_rule%lambda( 14,1:3) = (/ 2.6311282963463800E-01,  8.3947774099580007E-03,  7.2849239295540402E-01/)
       V_rule%lambda( 15,1:3) = (/ 7.2849239295540402E-01,  2.6311282963463800E-01,  8.3947774099579764E-03/)
       V_rule%lambda( 16,1:3) = (/ 8.3947774099580007E-03,  7.2849239295540402E-01,  2.6311282963463800E-01/)


    case(  9  )
       V_rule%Qdof =  19
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  9.7135796282799003E-02
       V_rule%weights(  2) =  3.1334700227139002E-02
       V_rule%weights(  3) =  3.1334700227139002E-02
       V_rule%weights(  4) =  3.1334700227139002E-02
       V_rule%weights(  5) =  7.7827541004774001E-02
       V_rule%weights(  6) =  7.7827541004774001E-02
       V_rule%weights(  7) =  7.7827541004774001E-02
       V_rule%weights(  8) =  7.9647738927209999E-02
       V_rule%weights(  9) =  7.9647738927209999E-02
       V_rule%weights( 10) =  7.9647738927209999E-02
       V_rule%weights( 11) =  2.5577675658698000E-02
       V_rule%weights( 12) =  2.5577675658698000E-02
       V_rule%weights( 13) =  2.5577675658698000E-02
       V_rule%weights( 14) =  4.3283539377289001E-02
       V_rule%weights( 15) =  4.3283539377289001E-02
       V_rule%weights( 16) =  4.3283539377289001E-02
       V_rule%weights( 17) =  4.3283539377289001E-02
       V_rule%weights( 18) =  4.3283539377289001E-02
       V_rule%weights( 19) =  4.3283539377289001E-02

       V_rule%lambda(  1,1:3) = (/ 3.3333333333333298E-01,  3.3333333333333298E-01,  3.3333333333333404E-01/)
       V_rule%lambda(  2,1:3) = (/ 2.0634961602524999E-02,  4.8968251919873801E-01,  4.8968251919873701E-01/)
       V_rule%lambda(  3,1:3) = (/ 4.8968251919873801E-01,  4.8968251919873801E-01,  2.0634961602523982E-02/)
       V_rule%lambda(  4,1:3) = (/ 4.8968251919873801E-01,  2.0634961602524999E-02,  4.8968251919873701E-01/)
       V_rule%lambda(  5,1:3) = (/ 1.2582081701412701E-01,  4.3708959149293702E-01,  4.3708959149293597E-01/)
       V_rule%lambda(  6,1:3) = (/ 4.3708959149293702E-01,  4.3708959149293702E-01,  1.2582081701412595E-01/)
       V_rule%lambda(  7,1:3) = (/ 4.3708959149293702E-01,  1.2582081701412701E-01,  4.3708959149293597E-01/)
       V_rule%lambda(  8,1:3) = (/ 6.2359292876193495E-01,  1.8820353561903300E-01,  1.8820353561903205E-01/)
       V_rule%lambda(  9,1:3) = (/ 1.8820353561903300E-01,  1.8820353561903300E-01,  6.2359292876193395E-01/)
       V_rule%lambda( 10,1:3) = (/ 1.8820353561903300E-01,  6.2359292876193495E-01,  1.8820353561903205E-01/)
       V_rule%lambda( 11,1:3) = (/ 9.1054097321109495E-01,  4.4729513394452997E-02,  4.4729513394452053E-02/)
       V_rule%lambda( 12,1:3) = (/ 4.4729513394452997E-02,  4.4729513394452997E-02,  9.1054097321109406E-01/)
       V_rule%lambda( 13,1:3) = (/ 4.4729513394452997E-02,  9.1054097321109495E-01,  4.4729513394452053E-02/)
       V_rule%lambda( 14,1:3) = (/ 3.6838412054736001E-02,  2.2196298916076601E-01,  7.4119859878449801E-01/)
       V_rule%lambda( 15,1:3) = (/ 2.2196298916076601E-01,  7.4119859878449801E-01,  3.6838412054735981E-02/)
       V_rule%lambda( 16,1:3) = (/ 7.4119859878449801E-01,  3.6838412054736001E-02,  2.2196298916076598E-01/)
       V_rule%lambda( 17,1:3) = (/ 2.2196298916076601E-01,  3.6838412054736001E-02,  7.4119859878449801E-01/)
       V_rule%lambda( 18,1:3) = (/ 7.4119859878449801E-01,  2.2196298916076601E-01,  3.6838412054735981E-02/)
       V_rule%lambda( 19,1:3) = (/ 3.6838412054736001E-02,  7.4119859878449801E-01,  2.2196298916076598E-01/)


    case( 10  )
       V_rule%Qdof =  25
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  9.0817990382753996E-02
       V_rule%weights(  2) =  3.6725957756466997E-02
       V_rule%weights(  3) =  3.6725957756466997E-02
       V_rule%weights(  4) =  3.6725957756466997E-02
       V_rule%weights(  5) =  4.5321059435527999E-02
       V_rule%weights(  6) =  4.5321059435527999E-02
       V_rule%weights(  7) =  4.5321059435527999E-02
       V_rule%weights(  8) =  7.2757916845420004E-02
       V_rule%weights(  9) =  7.2757916845420004E-02
       V_rule%weights( 10) =  7.2757916845420004E-02
       V_rule%weights( 11) =  7.2757916845420004E-02
       V_rule%weights( 12) =  7.2757916845420004E-02
       V_rule%weights( 13) =  7.2757916845420004E-02
       V_rule%weights( 14) =  2.8327242531056999E-02
       V_rule%weights( 15) =  2.8327242531056999E-02
       V_rule%weights( 16) =  2.8327242531056999E-02
       V_rule%weights( 17) =  2.8327242531056999E-02
       V_rule%weights( 18) =  2.8327242531056999E-02
       V_rule%weights( 19) =  2.8327242531056999E-02
       V_rule%weights( 20) =  9.4216669637330001E-03
       V_rule%weights( 21) =  9.4216669637330001E-03
       V_rule%weights( 22) =  9.4216669637330001E-03
       V_rule%weights( 23) =  9.4216669637330001E-03
       V_rule%weights( 24) =  9.4216669637330001E-03
       V_rule%weights( 25) =  9.4216669637330001E-03

       V_rule%lambda(  1,1:3) = (/ 3.3333333333333298E-01,  3.3333333333333298E-01,  3.3333333333333404E-01/)
       V_rule%lambda(  2,1:3) = (/ 2.8844733232685001E-02,  4.8557763338365700E-01,  4.8557763338365800E-01/)
       V_rule%lambda(  3,1:3) = (/ 4.8557763338365700E-01,  4.8557763338365700E-01,  2.8844733232685993E-02/)
       V_rule%lambda(  4,1:3) = (/ 4.8557763338365700E-01,  2.8844733232685001E-02,  4.8557763338365800E-01/)
       V_rule%lambda(  5,1:3) = (/ 7.8103684902992598E-01,  1.0948157548503699E-01,  1.0948157548503702E-01/)
       V_rule%lambda(  6,1:3) = (/ 1.0948157548503699E-01,  1.0948157548503699E-01,  7.8103684902992598E-01/)
       V_rule%lambda(  7,1:3) = (/ 1.0948157548503699E-01,  7.8103684902992598E-01,  1.0948157548503702E-01/)
       V_rule%lambda(  8,1:3) = (/ 1.4170721941487999E-01,  3.0793983876412101E-01,  5.5035294182099903E-01/)
       V_rule%lambda(  9,1:3) = (/ 3.0793983876412101E-01,  5.5035294182099903E-01,  1.4170721941487996E-01/)
       V_rule%lambda( 10,1:3) = (/ 5.5035294182099903E-01,  1.4170721941487999E-01,  3.0793983876412101E-01/)
       V_rule%lambda( 11,1:3) = (/ 3.0793983876412101E-01,  1.4170721941487999E-01,  5.5035294182099903E-01/)
       V_rule%lambda( 12,1:3) = (/ 5.5035294182099903E-01,  3.0793983876412101E-01,  1.4170721941487996E-01/)
       V_rule%lambda( 13,1:3) = (/ 1.4170721941487999E-01,  5.5035294182099903E-01,  3.0793983876412101E-01/)
       V_rule%lambda( 14,1:3) = (/ 2.5003534762685999E-02,  2.4667256063990300E-01,  7.2832390459741103E-01/)
       V_rule%lambda( 15,1:3) = (/ 2.4667256063990300E-01,  7.2832390459741103E-01,  2.5003534762685964E-02/)
       V_rule%lambda( 16,1:3) = (/ 7.2832390459741103E-01,  2.5003534762685999E-02,  2.4667256063990298E-01/)
       V_rule%lambda( 17,1:3) = (/ 2.4667256063990300E-01,  2.5003534762685999E-02,  7.2832390459741103E-01/)
       V_rule%lambda( 18,1:3) = (/ 7.2832390459741103E-01,  2.4667256063990300E-01,  2.5003534762685964E-02/)
       V_rule%lambda( 19,1:3) = (/ 2.5003534762685999E-02,  7.2832390459741103E-01,  2.4667256063990298E-01/)
       V_rule%lambda( 20,1:3) = (/ 9.5408154002989999E-03,  6.6803251012199999E-02,  9.2365593358750098E-01/)
       V_rule%lambda( 21,1:3) = (/ 6.6803251012199999E-02,  9.2365593358749998E-01,  9.5408154003000234E-03/)
       V_rule%lambda( 22,1:3) = (/ 9.2365593358749998E-01,  9.5408154002989999E-03,  6.6803251012201026E-02/)
       V_rule%lambda( 23,1:3) = (/ 6.6803251012199999E-02,  9.5408154002989999E-03,  9.2365593358750098E-01/)
       V_rule%lambda( 24,1:3) = (/ 9.2365593358749998E-01,  6.6803251012199999E-02,  9.5408154003000234E-03/)
       V_rule%lambda( 25,1:3) = (/ 9.5408154002989999E-03,  9.2365593358749998E-01,  6.6803251012201026E-02/)


    case( 11  )
       V_rule%Qdof =  27
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  9.2700632896100001E-04
       V_rule%weights(  2) =  9.2700632896100001E-04
       V_rule%weights(  3) =  9.2700632896100001E-04
       V_rule%weights(  4) =  7.7149534914812995E-02
       V_rule%weights(  5) =  7.7149534914812995E-02
       V_rule%weights(  6) =  7.7149534914812995E-02
       V_rule%weights(  7) =  5.9322977380774002E-02
       V_rule%weights(  8) =  5.9322977380774002E-02
       V_rule%weights(  9) =  5.9322977380774002E-02
       V_rule%weights( 10) =  3.6184540503418003E-02
       V_rule%weights( 11) =  3.6184540503418003E-02
       V_rule%weights( 12) =  3.6184540503418003E-02
       V_rule%weights( 13) =  1.3659731002678000E-02
       V_rule%weights( 14) =  1.3659731002678000E-02
       V_rule%weights( 15) =  1.3659731002678000E-02
       V_rule%weights( 16) =  5.2337111962204003E-02
       V_rule%weights( 17) =  5.2337111962204003E-02
       V_rule%weights( 18) =  5.2337111962204003E-02
       V_rule%weights( 19) =  5.2337111962204003E-02
       V_rule%weights( 20) =  5.2337111962204003E-02
       V_rule%weights( 21) =  5.2337111962204003E-02
       V_rule%weights( 22) =  2.0707659639141000E-02
       V_rule%weights( 23) =  2.0707659639141000E-02
       V_rule%weights( 24) =  2.0707659639141000E-02
       V_rule%weights( 25) =  2.0707659639141000E-02
       V_rule%weights( 26) =  2.0707659639141000E-02
       V_rule%weights( 27) =  2.0707659639141000E-02

       V_rule%lambda(  1,1:3) = (/-6.9222096541516995E-02,  5.3461104827075800E-01,  5.3461104827075900E-01/)
       V_rule%lambda(  2,1:3) = (/ 5.3461104827075800E-01,  5.3461104827075800E-01, -6.9222096541516009E-02/)
       V_rule%lambda(  3,1:3) = (/ 5.3461104827075800E-01, -6.9222096541516995E-02,  5.3461104827075900E-01/)
       V_rule%lambda(  4,1:3) = (/ 2.0206139406828999E-01,  3.9896930296585498E-01,  3.9896930296585503E-01/)
       V_rule%lambda(  5,1:3) = (/ 3.9896930296585498E-01,  3.9896930296585498E-01,  2.0206139406829005E-01/)
       V_rule%lambda(  6,1:3) = (/ 3.9896930296585498E-01,  2.0206139406828999E-01,  3.9896930296585503E-01/)
       V_rule%lambda(  7,1:3) = (/ 5.9338019913743500E-01,  2.0330990043128200E-01,  2.0330990043128300E-01/)
       V_rule%lambda(  8,1:3) = (/ 2.0330990043128200E-01,  2.0330990043128200E-01,  5.9338019913743600E-01/)
       V_rule%lambda(  9,1:3) = (/ 2.0330990043128200E-01,  5.9338019913743500E-01,  2.0330990043128300E-01/)
       V_rule%lambda( 10,1:3) = (/ 7.6129817543483702E-01,  1.1935091228258100E-01,  1.1935091228258198E-01/)
       V_rule%lambda( 11,1:3) = (/ 1.1935091228258100E-01,  1.1935091228258100E-01,  7.6129817543483802E-01/)
       V_rule%lambda( 12,1:3) = (/ 1.1935091228258100E-01,  7.6129817543483702E-01,  1.1935091228258198E-01/)
       V_rule%lambda( 13,1:3) = (/ 9.3527010377744801E-01,  3.2364948111276000E-02,  3.2364948111275986E-02/)
       V_rule%lambda( 14,1:3) = (/ 3.2364948111276000E-02,  3.2364948111276000E-02,  9.3527010377744801E-01/)
       V_rule%lambda( 15,1:3) = (/ 3.2364948111276000E-02,  9.3527010377744801E-01,  3.2364948111275986E-02/)
       V_rule%lambda( 16,1:3) = (/ 5.0178138310495002E-02,  3.5662064826129303E-01,  5.9320121342821197E-01/)
       V_rule%lambda( 17,1:3) = (/ 3.5662064826129303E-01,  5.9320121342821297E-01,  5.0178138310494003E-02/)
       V_rule%lambda( 18,1:3) = (/ 5.9320121342821297E-01,  5.0178138310495002E-02,  3.5662064826129203E-01/)
       V_rule%lambda( 19,1:3) = (/ 3.5662064826129303E-01,  5.0178138310495002E-02,  5.9320121342821197E-01/)
       V_rule%lambda( 20,1:3) = (/ 5.9320121342821297E-01,  3.5662064826129303E-01,  5.0178138310494003E-02/)
       V_rule%lambda( 21,1:3) = (/ 5.0178138310495002E-02,  5.9320121342821297E-01,  3.5662064826129203E-01/)
       V_rule%lambda( 22,1:3) = (/ 2.1022016536166001E-02,  1.7148898030404200E-01,  8.0748900315979200E-01/)
       V_rule%lambda( 23,1:3) = (/ 1.7148898030404200E-01,  8.0748900315979200E-01,  2.1022016536166005E-02/)
       V_rule%lambda( 24,1:3) = (/ 8.0748900315979200E-01,  2.1022016536166001E-02,  1.7148898030404200E-01/)
       V_rule%lambda( 25,1:3) = (/ 1.7148898030404200E-01,  2.1022016536166001E-02,  8.0748900315979200E-01/)
       V_rule%lambda( 26,1:3) = (/ 8.0748900315979200E-01,  1.7148898030404200E-01,  2.1022016536166005E-02/)
       V_rule%lambda( 27,1:3) = (/ 2.1022016536166001E-02,  8.0748900315979200E-01,  1.7148898030404200E-01/)


    case( 12  )
       V_rule%Qdof =  33
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  2.5731066440454999E-02
       V_rule%weights(  2) =  2.5731066440454999E-02
       V_rule%weights(  3) =  2.5731066440454999E-02
       V_rule%weights(  4) =  4.3692544538038003E-02
       V_rule%weights(  5) =  4.3692544538038003E-02
       V_rule%weights(  6) =  4.3692544538038003E-02
       V_rule%weights(  7) =  6.2858224217885006E-02
       V_rule%weights(  8) =  6.2858224217885006E-02
       V_rule%weights(  9) =  6.2858224217885006E-02
       V_rule%weights( 10) =  3.4796112930709000E-02
       V_rule%weights( 11) =  3.4796112930709000E-02
       V_rule%weights( 12) =  3.4796112930709000E-02
       V_rule%weights( 13) =  6.1662610515590003E-03
       V_rule%weights( 14) =  6.1662610515590003E-03
       V_rule%weights( 15) =  6.1662610515590003E-03
       V_rule%weights( 16) =  4.0371557766381003E-02
       V_rule%weights( 17) =  4.0371557766381003E-02
       V_rule%weights( 18) =  4.0371557766381003E-02
       V_rule%weights( 19) =  4.0371557766381003E-02
       V_rule%weights( 20) =  4.0371557766381003E-02
       V_rule%weights( 21) =  4.0371557766381003E-02
       V_rule%weights( 22) =  2.2356773202303001E-02
       V_rule%weights( 23) =  2.2356773202303001E-02
       V_rule%weights( 24) =  2.2356773202303001E-02
       V_rule%weights( 25) =  2.2356773202303001E-02
       V_rule%weights( 26) =  2.2356773202303001E-02
       V_rule%weights( 27) =  2.2356773202303001E-02
       V_rule%weights( 28) =  1.7316231108659000E-02
       V_rule%weights( 29) =  1.7316231108659000E-02
       V_rule%weights( 30) =  1.7316231108659000E-02
       V_rule%weights( 31) =  1.7316231108659000E-02
       V_rule%weights( 32) =  1.7316231108659000E-02
       V_rule%weights( 33) =  1.7316231108659000E-02

       V_rule%lambda(  1,1:3) = (/ 2.3565220452389998E-02,  4.8821738977380502E-01,  4.8821738977380497E-01/)
       V_rule%lambda(  2,1:3) = (/ 4.8821738977380502E-01,  4.8821738977380502E-01,  2.3565220452389957E-02/)
       V_rule%lambda(  3,1:3) = (/ 4.8821738977380502E-01,  2.3565220452389998E-02,  4.8821738977380497E-01/)
       V_rule%lambda(  4,1:3) = (/ 1.2055121541107899E-01,  4.3972439229445998E-01,  4.3972439229446103E-01/)
       V_rule%lambda(  5,1:3) = (/ 4.3972439229445998E-01,  4.3972439229445998E-01,  1.2055121541108005E-01/)
       V_rule%lambda(  6,1:3) = (/ 4.3972439229445998E-01,  1.2055121541107899E-01,  4.3972439229446103E-01/)
       V_rule%lambda(  7,1:3) = (/ 4.5757922997576800E-01,  2.7121038501211597E-01,  2.7121038501211603E-01/)
       V_rule%lambda(  8,1:3) = (/ 2.7121038501211597E-01,  2.7121038501211597E-01,  4.5757922997576805E-01/)
       V_rule%lambda(  9,1:3) = (/ 2.7121038501211597E-01,  4.5757922997576800E-01,  2.7121038501211603E-01/)
       V_rule%lambda( 10,1:3) = (/ 7.4484770891682806E-01,  1.2757614554158600E-01,  1.2757614554158594E-01/)
       V_rule%lambda( 11,1:3) = (/ 1.2757614554158600E-01,  1.2757614554158600E-01,  7.4484770891682794E-01/)
       V_rule%lambda( 12,1:3) = (/ 1.2757614554158600E-01,  7.4484770891682806E-01,  1.2757614554158594E-01/)
       V_rule%lambda( 13,1:3) = (/ 9.5736529909357904E-01,  2.1317350453210000E-02,  2.1317350453210964E-02/)
       V_rule%lambda( 14,1:3) = (/ 2.1317350453210000E-02,  2.1317350453210000E-02,  9.5736529909358004E-01/)
       V_rule%lambda( 15,1:3) = (/ 2.1317350453210000E-02,  9.5736529909357904E-01,  2.1317350453210964E-02/)
       V_rule%lambda( 16,1:3) = (/ 1.1534349453469800E-01,  2.7571326968551402E-01,  6.0894323577978793E-01/)
       V_rule%lambda( 17,1:3) = (/ 2.7571326968551402E-01,  6.0894323577978804E-01,  1.1534349453469794E-01/)
       V_rule%lambda( 18,1:3) = (/ 6.0894323577978804E-01,  1.1534349453469800E-01,  2.7571326968551396E-01/)
       V_rule%lambda( 19,1:3) = (/ 2.7571326968551402E-01,  1.1534349453469800E-01,  6.0894323577978793E-01/)
       V_rule%lambda( 20,1:3) = (/ 6.0894323577978804E-01,  2.7571326968551402E-01,  1.1534349453469794E-01/)
       V_rule%lambda( 21,1:3) = (/ 1.1534349453469800E-01,  6.0894323577978804E-01,  2.7571326968551396E-01/)
       V_rule%lambda( 22,1:3) = (/ 2.2838332222257000E-02,  2.8132558098993998E-01,  6.9583608678780307E-01/)
       V_rule%lambda( 23,1:3) = (/ 2.8132558098993998E-01,  6.9583608678780295E-01,  2.2838332222257063E-02/)
       V_rule%lambda( 24,1:3) = (/ 6.9583608678780295E-01,  2.2838332222257000E-02,  2.8132558098994004E-01/)
       V_rule%lambda( 25,1:3) = (/ 2.8132558098993998E-01,  2.2838332222257000E-02,  6.9583608678780307E-01/)
       V_rule%lambda( 26,1:3) = (/ 6.9583608678780295E-01,  2.8132558098993998E-01,  2.2838332222257063E-02/)
       V_rule%lambda( 27,1:3) = (/ 2.2838332222257000E-02,  6.9583608678780295E-01,  2.8132558098994004E-01/)
       V_rule%lambda( 28,1:3) = (/ 2.5734050548330001E-02,  1.1625191590759699E-01,  8.5801403354407302E-01/)
       V_rule%lambda( 29,1:3) = (/ 1.1625191590759699E-01,  8.5801403354407302E-01,  2.5734050548329987E-02/)
       V_rule%lambda( 30,1:3) = (/ 8.5801403354407302E-01,  2.5734050548330001E-02,  1.1625191590759698E-01/)
       V_rule%lambda( 31,1:3) = (/ 1.1625191590759699E-01,  2.5734050548330001E-02,  8.5801403354407302E-01/)
       V_rule%lambda( 32,1:3) = (/ 8.5801403354407302E-01,  1.1625191590759699E-01,  2.5734050548329987E-02/)
       V_rule%lambda( 33,1:3) = (/ 2.5734050548330001E-02,  8.5801403354407302E-01,  1.1625191590759698E-01/)


    case( 13  )
       V_rule%Qdof =  37
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  5.2520923400802000E-02
       V_rule%weights(  2) =  1.1280145209330000E-02
       V_rule%weights(  3) =  1.1280145209330000E-02
       V_rule%weights(  4) =  1.1280145209330000E-02
       V_rule%weights(  5) =  3.1423518362454002E-02
       V_rule%weights(  6) =  3.1423518362454002E-02
       V_rule%weights(  7) =  3.1423518362454002E-02
       V_rule%weights(  8) =  4.7072502504194001E-02
       V_rule%weights(  9) =  4.7072502504194001E-02
       V_rule%weights( 10) =  4.7072502504194001E-02
       V_rule%weights( 11) =  4.7363586536355001E-02
       V_rule%weights( 12) =  4.7363586536355001E-02
       V_rule%weights( 13) =  4.7363586536355001E-02
       V_rule%weights( 14) =  3.1167529045794000E-02
       V_rule%weights( 15) =  3.1167529045794000E-02
       V_rule%weights( 16) =  3.1167529045794000E-02
       V_rule%weights( 17) =  7.9757714650739997E-03
       V_rule%weights( 18) =  7.9757714650739997E-03
       V_rule%weights( 19) =  7.9757714650739997E-03
       V_rule%weights( 20) =  3.6848402728732001E-02
       V_rule%weights( 21) =  3.6848402728732001E-02
       V_rule%weights( 22) =  3.6848402728732001E-02
       V_rule%weights( 23) =  3.6848402728732001E-02
       V_rule%weights( 24) =  3.6848402728732001E-02
       V_rule%weights( 25) =  3.6848402728732001E-02
       V_rule%weights( 26) =  1.7401463303822001E-02
       V_rule%weights( 27) =  1.7401463303822001E-02
       V_rule%weights( 28) =  1.7401463303822001E-02
       V_rule%weights( 29) =  1.7401463303822001E-02
       V_rule%weights( 30) =  1.7401463303822001E-02
       V_rule%weights( 31) =  1.7401463303822001E-02
       V_rule%weights( 32) =  1.5521786839044999E-02
       V_rule%weights( 33) =  1.5521786839044999E-02
       V_rule%weights( 34) =  1.5521786839044999E-02
       V_rule%weights( 35) =  1.5521786839044999E-02
       V_rule%weights( 36) =  1.5521786839044999E-02
       V_rule%weights( 37) =  1.5521786839044999E-02

       V_rule%lambda(  1,1:3) = (/ 3.3333333333333298E-01,  3.3333333333333298E-01,  3.3333333333333404E-01/)
       V_rule%lambda(  2,1:3) = (/ 9.9036301205910008E-03,  4.9504818493970498E-01,  4.9504818493970404E-01/)
       V_rule%lambda(  3,1:3) = (/ 4.9504818493970498E-01,  4.9504818493970498E-01,  9.9036301205900346E-03/)
       V_rule%lambda(  4,1:3) = (/ 4.9504818493970498E-01,  9.9036301205910008E-03,  4.9504818493970404E-01/)
       V_rule%lambda(  5,1:3) = (/ 6.2566729780851996E-02,  4.6871663510957401E-01,  4.6871663510957401E-01/)
       V_rule%lambda(  6,1:3) = (/ 4.6871663510957401E-01,  4.6871663510957401E-01,  6.2566729780851982E-02/)
       V_rule%lambda(  7,1:3) = (/ 4.6871663510957401E-01,  6.2566729780851996E-02,  4.6871663510957401E-01/)
       V_rule%lambda(  8,1:3) = (/ 1.7095732639744701E-01,  4.1452133680127701E-01,  4.1452133680127601E-01/)
       V_rule%lambda(  9,1:3) = (/ 4.1452133680127701E-01,  4.1452133680127701E-01,  1.7095732639744599E-01/)
       V_rule%lambda( 10,1:3) = (/ 4.1452133680127701E-01,  1.7095732639744701E-01,  4.1452133680127601E-01/)
       V_rule%lambda( 11,1:3) = (/ 5.4120085591433698E-01,  2.2939957204283101E-01,  2.2939957204283201E-01/)
       V_rule%lambda( 12,1:3) = (/ 2.2939957204283101E-01,  2.2939957204283101E-01,  5.4120085591433797E-01/)
       V_rule%lambda( 13,1:3) = (/ 2.2939957204283101E-01,  5.4120085591433698E-01,  2.2939957204283201E-01/)
       V_rule%lambda( 14,1:3) = (/ 7.7115100960733995E-01,  1.1442449519632999E-01,  1.1442449519633005E-01/)
       V_rule%lambda( 15,1:3) = (/ 1.1442449519632999E-01,  1.1442449519632999E-01,  7.7115100960733995E-01/)
       V_rule%lambda( 16,1:3) = (/ 1.1442449519632999E-01,  7.7115100960733995E-01,  1.1442449519633005E-01/)
       V_rule%lambda( 17,1:3) = (/ 9.5037721727308200E-01,  2.4811391363459001E-02,  2.4811391363459001E-02/)
       V_rule%lambda( 18,1:3) = (/ 2.4811391363459001E-02,  2.4811391363459001E-02,  9.5037721727308200E-01/)
       V_rule%lambda( 19,1:3) = (/ 2.4811391363459001E-02,  9.5037721727308200E-01,  2.4811391363459001E-02/)
       V_rule%lambda( 20,1:3) = (/ 9.4853828379578994E-02,  2.6879499705876098E-01,  6.3635117456166002E-01/)
       V_rule%lambda( 21,1:3) = (/ 2.6879499705876098E-01,  6.3635117456166002E-01,  9.4853828379579008E-02/)
       V_rule%lambda( 22,1:3) = (/ 6.3635117456166002E-01,  9.4853828379578994E-02,  2.6879499705876098E-01/)
       V_rule%lambda( 23,1:3) = (/ 2.6879499705876098E-01,  9.4853828379578994E-02,  6.3635117456166002E-01/)
       V_rule%lambda( 24,1:3) = (/ 6.3635117456166002E-01,  2.6879499705876098E-01,  9.4853828379579008E-02/)
       V_rule%lambda( 25,1:3) = (/ 9.4853828379578994E-02,  6.3635117456166002E-01,  2.6879499705876098E-01/)
       V_rule%lambda( 26,1:3) = (/ 1.8100773278806999E-02,  2.9173006673428797E-01,  6.9016915998690498E-01/)
       V_rule%lambda( 27,1:3) = (/ 2.9173006673428797E-01,  6.9016915998690498E-01,  1.8100773278807047E-02/)
       V_rule%lambda( 28,1:3) = (/ 6.9016915998690498E-01,  1.8100773278806999E-02,  2.9173006673428803E-01/)
       V_rule%lambda( 29,1:3) = (/ 2.9173006673428797E-01,  1.8100773278806999E-02,  6.9016915998690498E-01/)
       V_rule%lambda( 30,1:3) = (/ 6.9016915998690498E-01,  2.9173006673428797E-01,  1.8100773278807047E-02/)
       V_rule%lambda( 31,1:3) = (/ 1.8100773278806999E-02,  6.9016915998690498E-01,  2.9173006673428803E-01/)
       V_rule%lambda( 32,1:3) = (/ 2.2233076674090000E-02,  1.2635738549166900E-01,  8.5140953783424100E-01/)
       V_rule%lambda( 33,1:3) = (/ 1.2635738549166900E-01,  8.5140953783424100E-01,  2.2233076674089997E-02/)
       V_rule%lambda( 34,1:3) = (/ 8.5140953783424100E-01,  2.2233076674090000E-02,  1.2635738549166900E-01/)
       V_rule%lambda( 35,1:3) = (/ 1.2635738549166900E-01,  2.2233076674090000E-02,  8.5140953783424100E-01/)
       V_rule%lambda( 36,1:3) = (/ 8.5140953783424100E-01,  1.2635738549166900E-01,  2.2233076674089997E-02/)
       V_rule%lambda( 37,1:3) = (/ 2.2233076674090000E-02,  8.5140953783424100E-01,  1.2635738549166900E-01/)


    case( 14  )
       V_rule%Qdof =  42
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  2.1883581369429000E-02
       V_rule%weights(  2) =  2.1883581369429000E-02
       V_rule%weights(  3) =  2.1883581369429000E-02
       V_rule%weights(  4) =  3.2788353544125001E-02
       V_rule%weights(  5) =  3.2788353544125001E-02
       V_rule%weights(  6) =  3.2788353544125001E-02
       V_rule%weights(  7) =  5.1774104507292001E-02
       V_rule%weights(  8) =  5.1774104507292001E-02
       V_rule%weights(  9) =  5.1774104507292001E-02
       V_rule%weights( 10) =  4.2162588736993002E-02
       V_rule%weights( 11) =  4.2162588736993002E-02
       V_rule%weights( 12) =  4.2162588736993002E-02
       V_rule%weights( 13) =  1.4433699669777001E-02
       V_rule%weights( 14) =  1.4433699669777001E-02
       V_rule%weights( 15) =  1.4433699669777001E-02
       V_rule%weights( 16) =  4.9234036024000003E-03
       V_rule%weights( 17) =  4.9234036024000003E-03
       V_rule%weights( 18) =  4.9234036024000003E-03
       V_rule%weights( 19) =  2.4665753212564000E-02
       V_rule%weights( 20) =  2.4665753212564000E-02
       V_rule%weights( 21) =  2.4665753212564000E-02
       V_rule%weights( 22) =  2.4665753212564000E-02
       V_rule%weights( 23) =  2.4665753212564000E-02
       V_rule%weights( 24) =  2.4665753212564000E-02
       V_rule%weights( 25) =  3.8571510787061003E-02
       V_rule%weights( 26) =  3.8571510787061003E-02
       V_rule%weights( 27) =  3.8571510787061003E-02
       V_rule%weights( 28) =  3.8571510787061003E-02
       V_rule%weights( 29) =  3.8571510787061003E-02
       V_rule%weights( 30) =  3.8571510787061003E-02
       V_rule%weights( 31) =  1.4436308113534000E-02
       V_rule%weights( 32) =  1.4436308113534000E-02
       V_rule%weights( 33) =  1.4436308113534000E-02
       V_rule%weights( 34) =  1.4436308113534000E-02
       V_rule%weights( 35) =  1.4436308113534000E-02
       V_rule%weights( 36) =  1.4436308113534000E-02
       V_rule%weights( 37) =  5.0102288385009998E-03
       V_rule%weights( 38) =  5.0102288385009998E-03
       V_rule%weights( 39) =  5.0102288385009998E-03
       V_rule%weights( 40) =  5.0102288385009998E-03
       V_rule%weights( 41) =  5.0102288385009998E-03
       V_rule%weights( 42) =  5.0102288385009998E-03

       V_rule%lambda(  1,1:3) = (/ 2.2072179275642999E-02,  4.8896391036217901E-01,  4.8896391036217801E-01/)
       V_rule%lambda(  2,1:3) = (/ 4.8896391036217901E-01,  4.8896391036217901E-01,  2.2072179275641979E-02/)
       V_rule%lambda(  3,1:3) = (/ 4.8896391036217901E-01,  2.2072179275642999E-02,  4.8896391036217801E-01/)
       V_rule%lambda(  4,1:3) = (/ 1.6471056131909201E-01,  4.1764471934045400E-01,  4.1764471934045400E-01/)
       V_rule%lambda(  5,1:3) = (/ 4.1764471934045400E-01,  4.1764471934045400E-01,  1.6471056131909201E-01/)
       V_rule%lambda(  6,1:3) = (/ 4.1764471934045400E-01,  1.6471056131909201E-01,  4.1764471934045400E-01/)
       V_rule%lambda(  7,1:3) = (/ 4.5304494338232298E-01,  2.7347752830883898E-01,  2.7347752830883804E-01/)
       V_rule%lambda(  8,1:3) = (/ 2.7347752830883898E-01,  2.7347752830883898E-01,  4.5304494338232204E-01/)
       V_rule%lambda(  9,1:3) = (/ 2.7347752830883898E-01,  4.5304494338232298E-01,  2.7347752830883804E-01/)
       V_rule%lambda( 10,1:3) = (/ 6.4558893517491300E-01,  1.7720553241254300E-01,  1.7720553241254400E-01/)
       V_rule%lambda( 11,1:3) = (/ 1.7720553241254300E-01,  1.7720553241254300E-01,  6.4558893517491400E-01/)
       V_rule%lambda( 12,1:3) = (/ 1.7720553241254300E-01,  6.4558893517491300E-01,  1.7720553241254400E-01/)
       V_rule%lambda( 13,1:3) = (/ 8.7640023381825505E-01,  6.1799883090873003E-02,  6.1799883090871949E-02/)
       V_rule%lambda( 14,1:3) = (/ 6.1799883090873003E-02,  6.1799883090873003E-02,  8.7640023381825394E-01/)
       V_rule%lambda( 15,1:3) = (/ 6.1799883090873003E-02,  8.7640023381825505E-01,  6.1799883090871949E-02/)
       V_rule%lambda( 16,1:3) = (/ 9.6121807750259802E-01,  1.9390961248700999E-02,  1.9390961248700978E-02/)
       V_rule%lambda( 17,1:3) = (/ 1.9390961248700999E-02,  1.9390961248700999E-02,  9.6121807750259802E-01/)
       V_rule%lambda( 18,1:3) = (/ 1.9390961248700999E-02,  9.6121807750259802E-01,  1.9390961248700978E-02/)
       V_rule%lambda( 19,1:3) = (/ 5.7124757403648002E-02,  1.7226668782135601E-01,  7.7060855477499601E-01/)
       V_rule%lambda( 20,1:3) = (/ 1.7226668782135601E-01,  7.7060855477499601E-01,  5.7124757403647974E-02/)
       V_rule%lambda( 21,1:3) = (/ 7.7060855477499601E-01,  5.7124757403648002E-02,  1.7226668782135598E-01/)
       V_rule%lambda( 22,1:3) = (/ 1.7226668782135601E-01,  5.7124757403648002E-02,  7.7060855477499601E-01/)
       V_rule%lambda( 23,1:3) = (/ 7.7060855477499601E-01,  1.7226668782135601E-01,  5.7124757403647974E-02/)
       V_rule%lambda( 24,1:3) = (/ 5.7124757403648002E-02,  7.7060855477499601E-01,  1.7226668782135598E-01/)
       V_rule%lambda( 25,1:3) = (/ 9.2916249356972000E-02,  3.3686145979634502E-01,  5.7022229084668297E-01/)
       V_rule%lambda( 26,1:3) = (/ 3.3686145979634502E-01,  5.7022229084668297E-01,  9.2916249356972014E-02/)
       V_rule%lambda( 27,1:3) = (/ 5.7022229084668297E-01,  9.2916249356972000E-02,  3.3686145979634502E-01/)
       V_rule%lambda( 28,1:3) = (/ 3.3686145979634502E-01,  9.2916249356972000E-02,  5.7022229084668297E-01/)
       V_rule%lambda( 29,1:3) = (/ 5.7022229084668297E-01,  3.3686145979634502E-01,  9.2916249356972014E-02/)
       V_rule%lambda( 30,1:3) = (/ 9.2916249356972000E-02,  5.7022229084668297E-01,  3.3686145979634502E-01/)
       V_rule%lambda( 31,1:3) = (/ 1.4646950055653999E-02,  2.9837288213625801E-01,  6.8698016780808802E-01/)
       V_rule%lambda( 32,1:3) = (/ 2.9837288213625801E-01,  6.8698016780808802E-01,  1.4646950055653973E-02/)
       V_rule%lambda( 33,1:3) = (/ 6.8698016780808802E-01,  1.4646950055653999E-02,  2.9837288213625801E-01/)
       V_rule%lambda( 34,1:3) = (/ 2.9837288213625801E-01,  1.4646950055653999E-02,  6.8698016780808802E-01/)
       V_rule%lambda( 35,1:3) = (/ 6.8698016780808802E-01,  2.9837288213625801E-01,  1.4646950055653973E-02/)
       V_rule%lambda( 36,1:3) = (/ 1.4646950055653999E-02,  6.8698016780808802E-01,  2.9837288213625801E-01/)
       V_rule%lambda( 37,1:3) = (/ 1.2683309328719999E-03,  1.1897449769695700E-01,  8.7975717137017095E-01/)
       V_rule%lambda( 38,1:3) = (/ 1.1897449769695700E-01,  8.7975717137017095E-01,  1.2683309328720416E-03/)
       V_rule%lambda( 39,1:3) = (/ 8.7975717137017095E-01,  1.2683309328719999E-03,  1.1897449769695705E-01/)
       V_rule%lambda( 40,1:3) = (/ 1.1897449769695700E-01,  1.2683309328719999E-03,  8.7975717137017095E-01/)
       V_rule%lambda( 41,1:3) = (/ 8.7975717137017095E-01,  1.1897449769695700E-01,  1.2683309328720416E-03/)
       V_rule%lambda( 42,1:3) = (/ 1.2683309328719999E-03,  8.7975717137017095E-01,  1.1897449769695705E-01/)


    case( 15  )
       V_rule%Qdof =  48
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  1.9168756428489999E-03
       V_rule%weights(  2) =  1.9168756428489999E-03
       V_rule%weights(  3) =  1.9168756428489999E-03
       V_rule%weights(  4) =  4.4249027271144999E-02
       V_rule%weights(  5) =  4.4249027271144999E-02
       V_rule%weights(  6) =  4.4249027271144999E-02
       V_rule%weights(  7) =  5.1186548718851997E-02
       V_rule%weights(  8) =  5.1186548718851997E-02
       V_rule%weights(  9) =  5.1186548718851997E-02
       V_rule%weights( 10) =  2.3687735870687999E-02
       V_rule%weights( 11) =  2.3687735870687999E-02
       V_rule%weights( 12) =  2.3687735870687999E-02
       V_rule%weights( 13) =  1.3289775690021001E-02
       V_rule%weights( 14) =  1.3289775690021001E-02
       V_rule%weights( 15) =  1.3289775690021001E-02
       V_rule%weights( 16) =  4.7489166081920000E-03
       V_rule%weights( 17) =  4.7489166081920000E-03
       V_rule%weights( 18) =  4.7489166081920000E-03
       V_rule%weights( 19) =  3.8550072599592998E-02
       V_rule%weights( 20) =  3.8550072599592998E-02
       V_rule%weights( 21) =  3.8550072599592998E-02
       V_rule%weights( 22) =  3.8550072599592998E-02
       V_rule%weights( 23) =  3.8550072599592998E-02
       V_rule%weights( 24) =  3.8550072599592998E-02
       V_rule%weights( 25) =  2.7215814320624001E-02
       V_rule%weights( 26) =  2.7215814320624001E-02
       V_rule%weights( 27) =  2.7215814320624001E-02
       V_rule%weights( 28) =  2.7215814320624001E-02
       V_rule%weights( 29) =  2.7215814320624001E-02
       V_rule%weights( 30) =  2.7215814320624001E-02
       V_rule%weights( 31) =  2.1820773667970000E-03
       V_rule%weights( 32) =  2.1820773667970000E-03
       V_rule%weights( 33) =  2.1820773667970000E-03
       V_rule%weights( 34) =  2.1820773667970000E-03
       V_rule%weights( 35) =  2.1820773667970000E-03
       V_rule%weights( 36) =  2.1820773667970000E-03
       V_rule%weights( 37) =  2.1505319847730998E-02
       V_rule%weights( 38) =  2.1505319847730998E-02
       V_rule%weights( 39) =  2.1505319847730998E-02
       V_rule%weights( 40) =  2.1505319847730998E-02
       V_rule%weights( 41) =  2.1505319847730998E-02
       V_rule%weights( 42) =  2.1505319847730998E-02
       V_rule%weights( 43) =  7.6739426310490000E-03
       V_rule%weights( 44) =  7.6739426310490000E-03
       V_rule%weights( 45) =  7.6739426310490000E-03
       V_rule%weights( 46) =  7.6739426310490000E-03
       V_rule%weights( 47) =  7.6739426310490000E-03
       V_rule%weights( 48) =  7.6739426310490000E-03

       V_rule%lambda(  1,1:3) = (/-1.3945833716486000E-02,  5.0697291685824297E-01,  5.0697291685824308E-01/)
       V_rule%lambda(  2,1:3) = (/ 5.0697291685824297E-01,  5.0697291685824297E-01, -1.3945833716485945E-02/)
       V_rule%lambda(  3,1:3) = (/ 5.0697291685824297E-01, -1.3945833716486000E-02,  5.0697291685824308E-01/)
       V_rule%lambda(  4,1:3) = (/ 1.3718729143395500E-01,  4.3140635428302299E-01,  4.3140635428302199E-01/)
       V_rule%lambda(  5,1:3) = (/ 4.3140635428302299E-01,  4.3140635428302299E-01,  1.3718729143395403E-01/)
       V_rule%lambda(  6,1:3) = (/ 4.3140635428302299E-01,  1.3718729143395500E-01,  4.3140635428302199E-01/)
       V_rule%lambda(  7,1:3) = (/ 4.4461271030571098E-01,  2.7769364484714398E-01,  2.7769364484714504E-01/)
       V_rule%lambda(  8,1:3) = (/ 2.7769364484714398E-01,  2.7769364484714398E-01,  4.4461271030571203E-01/)
       V_rule%lambda(  9,1:3) = (/ 2.7769364484714398E-01,  4.4461271030571098E-01,  2.7769364484714504E-01/)
       V_rule%lambda( 10,1:3) = (/ 7.4707021791749195E-01,  1.2646489104125400E-01,  1.2646489104125405E-01/)
       V_rule%lambda( 11,1:3) = (/ 1.2646489104125400E-01,  1.2646489104125400E-01,  7.4707021791749195E-01/)
       V_rule%lambda( 12,1:3) = (/ 1.2646489104125400E-01,  7.4707021791749195E-01,  1.2646489104125405E-01/)
       V_rule%lambda( 13,1:3) = (/ 8.5838322805062806E-01,  7.0808385974686000E-02,  7.0808385974685945E-02/)
       V_rule%lambda( 14,1:3) = (/ 7.0808385974686000E-02,  7.0808385974686000E-02,  8.5838322805062806E-01/)
       V_rule%lambda( 15,1:3) = (/ 7.0808385974686000E-02,  8.5838322805062806E-01,  7.0808385974685945E-02/)
       V_rule%lambda( 16,1:3) = (/ 9.6206965951785295E-01,  1.8965170241072998E-02,  1.8965170241074053E-02/)
       V_rule%lambda( 17,1:3) = (/ 1.8965170241072998E-02,  1.8965170241072998E-02,  9.6206965951785395E-01/)
       V_rule%lambda( 18,1:3) = (/ 1.8965170241072998E-02,  9.6206965951785295E-01,  1.8965170241074053E-02/)
       V_rule%lambda( 19,1:3) = (/ 1.3373416196662100E-01,  2.6131137114008701E-01,  6.0495446689329202E-01/)
       V_rule%lambda( 20,1:3) = (/ 2.6131137114008701E-01,  6.0495446689329102E-01,  1.3373416196662197E-01/)
       V_rule%lambda( 21,1:3) = (/ 6.0495446689329102E-01,  1.3373416196662100E-01,  2.6131137114008796E-01/)
       V_rule%lambda( 22,1:3) = (/ 2.6131137114008701E-01,  1.3373416196662100E-01,  6.0495446689329202E-01/)
       V_rule%lambda( 23,1:3) = (/ 6.0495446689329102E-01,  2.6131137114008701E-01,  1.3373416196662197E-01/)
       V_rule%lambda( 24,1:3) = (/ 1.3373416196662100E-01,  6.0495446689329102E-01,  2.6131137114008796E-01/)
       V_rule%lambda( 25,1:3) = (/ 3.6366677396916999E-02,  3.8804676709026897E-01,  5.7558655551281401E-01/)
       V_rule%lambda( 26,1:3) = (/ 3.8804676709026897E-01,  5.7558655551281401E-01,  3.6366677396917013E-02/)
       V_rule%lambda( 27,1:3) = (/ 5.7558655551281401E-01,  3.6366677396916999E-02,  3.8804676709026897E-01/)
       V_rule%lambda( 28,1:3) = (/ 3.8804676709026897E-01,  3.6366677396916999E-02,  5.7558655551281401E-01/)
       V_rule%lambda( 29,1:3) = (/ 5.7558655551281401E-01,  3.8804676709026897E-01,  3.6366677396917013E-02/)
       V_rule%lambda( 30,1:3) = (/ 3.6366677396916999E-02,  5.7558655551281401E-01,  3.8804676709026897E-01/)
       V_rule%lambda( 31,1:3) = (/-1.0174883126570999E-02,  2.8571222004991598E-01,  7.2446266307665508E-01/)
       V_rule%lambda( 32,1:3) = (/ 2.8571222004991598E-01,  7.2446266307665497E-01, -1.0174883126570944E-02/)
       V_rule%lambda( 33,1:3) = (/ 7.2446266307665497E-01, -1.0174883126570999E-02,  2.8571222004991603E-01/)
       V_rule%lambda( 34,1:3) = (/ 2.8571222004991598E-01, -1.0174883126570999E-02,  7.2446266307665508E-01/)
       V_rule%lambda( 35,1:3) = (/ 7.2446266307665497E-01,  2.8571222004991598E-01, -1.0174883126570944E-02/)
       V_rule%lambda( 36,1:3) = (/-1.0174883126570999E-02,  7.2446266307665497E-01,  2.8571222004991603E-01/)
       V_rule%lambda( 37,1:3) = (/ 3.6843869875878003E-02,  2.1559966407228401E-01,  7.4755646605183801E-01/)
       V_rule%lambda( 38,1:3) = (/ 2.1559966407228401E-01,  7.4755646605183801E-01,  3.6843869875877983E-02/)
       V_rule%lambda( 39,1:3) = (/ 7.4755646605183801E-01,  3.6843869875878003E-02,  2.1559966407228398E-01/)
       V_rule%lambda( 40,1:3) = (/ 2.1559966407228401E-01,  3.6843869875878003E-02,  7.4755646605183801E-01/)
       V_rule%lambda( 41,1:3) = (/ 7.4755646605183801E-01,  2.1559966407228401E-01,  3.6843869875877983E-02/)
       V_rule%lambda( 42,1:3) = (/ 3.6843869875878003E-02,  7.4755646605183801E-01,  2.1559966407228398E-01/)
       V_rule%lambda( 43,1:3) = (/ 1.2459809331198999E-02,  1.0357561657638600E-01,  8.8396457409241502E-01/)
       V_rule%lambda( 44,1:3) = (/ 1.0357561657638600E-01,  8.8396457409241602E-01,  1.2459809331197974E-02/)
       V_rule%lambda( 45,1:3) = (/ 8.8396457409241602E-01,  1.2459809331198999E-02,  1.0357561657638498E-01/)
       V_rule%lambda( 46,1:3) = (/ 1.0357561657638600E-01,  1.2459809331198999E-02,  8.8396457409241502E-01/)
       V_rule%lambda( 47,1:3) = (/ 8.8396457409241602E-01,  1.0357561657638600E-01,  1.2459809331197974E-02/)
       V_rule%lambda( 48,1:3) = (/ 1.2459809331198999E-02,  8.8396457409241602E-01,  1.0357561657638498E-01/)


    case( 16  )
       V_rule%Qdof =  52
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  4.6875697427641999E-02
       V_rule%weights(  2) =  6.4058785785849996E-03
       V_rule%weights(  3) =  6.4058785785849996E-03
       V_rule%weights(  4) =  6.4058785785849996E-03
       V_rule%weights(  5) =  4.1710296739386997E-02
       V_rule%weights(  6) =  4.1710296739386997E-02
       V_rule%weights(  7) =  4.1710296739386997E-02
       V_rule%weights(  8) =  2.6891484250064001E-02
       V_rule%weights(  9) =  2.6891484250064001E-02
       V_rule%weights( 10) =  2.6891484250064001E-02
       V_rule%weights( 11) =  4.2132522761649999E-02
       V_rule%weights( 12) =  4.2132522761649999E-02
       V_rule%weights( 13) =  4.2132522761649999E-02
       V_rule%weights( 14) =  3.0000266842772998E-02
       V_rule%weights( 15) =  3.0000266842772998E-02
       V_rule%weights( 16) =  3.0000266842772998E-02
       V_rule%weights( 17) =  1.4200098925024000E-02
       V_rule%weights( 18) =  1.4200098925024000E-02
       V_rule%weights( 19) =  1.4200098925024000E-02
       V_rule%weights( 20) =  3.5824623512730001E-03
       V_rule%weights( 21) =  3.5824623512730001E-03
       V_rule%weights( 22) =  3.5824623512730001E-03
       V_rule%weights( 23) =  3.2773147460626997E-02
       V_rule%weights( 24) =  3.2773147460626997E-02
       V_rule%weights( 25) =  3.2773147460626997E-02
       V_rule%weights( 26) =  3.2773147460626997E-02
       V_rule%weights( 27) =  3.2773147460626997E-02
       V_rule%weights( 28) =  3.2773147460626997E-02
       V_rule%weights( 29) =  1.5298306248441000E-02
       V_rule%weights( 30) =  1.5298306248441000E-02
       V_rule%weights( 31) =  1.5298306248441000E-02
       V_rule%weights( 32) =  1.5298306248441000E-02
       V_rule%weights( 33) =  1.5298306248441000E-02
       V_rule%weights( 34) =  1.5298306248441000E-02
       V_rule%weights( 35) =  2.3862441928390000E-03
       V_rule%weights( 36) =  2.3862441928390000E-03
       V_rule%weights( 37) =  2.3862441928390000E-03
       V_rule%weights( 38) =  2.3862441928390000E-03
       V_rule%weights( 39) =  2.3862441928390000E-03
       V_rule%weights( 40) =  2.3862441928390000E-03
       V_rule%weights( 41) =  1.9084792755898999E-02
       V_rule%weights( 42) =  1.9084792755898999E-02
       V_rule%weights( 43) =  1.9084792755898999E-02
       V_rule%weights( 44) =  1.9084792755898999E-02
       V_rule%weights( 45) =  1.9084792755898999E-02
       V_rule%weights( 46) =  1.9084792755898999E-02
       V_rule%weights( 47) =  6.8500545465420004E-03
       V_rule%weights( 48) =  6.8500545465420004E-03
       V_rule%weights( 49) =  6.8500545465420004E-03
       V_rule%weights( 50) =  6.8500545465420004E-03
       V_rule%weights( 51) =  6.8500545465420004E-03
       V_rule%weights( 52) =  6.8500545465420004E-03

       V_rule%lambda(  1,1:3) = (/ 3.3333333333333298E-01,  3.3333333333333298E-01,  3.3333333333333404E-01/)
       V_rule%lambda(  2,1:3) = (/ 5.2389161031230003E-03,  4.9738054194843800E-01,  4.9738054194843900E-01/)
       V_rule%lambda(  3,1:3) = (/ 4.9738054194843800E-01,  4.9738054194843800E-01,  5.2389161031239917E-03/)
       V_rule%lambda(  4,1:3) = (/ 4.9738054194843800E-01,  5.2389161031230003E-03,  4.9738054194843900E-01/)
       V_rule%lambda(  5,1:3) = (/ 1.7306112290129499E-01,  4.1346943854935198E-01,  4.1346943854935303E-01/)
       V_rule%lambda(  6,1:3) = (/ 4.1346943854935198E-01,  4.1346943854935198E-01,  1.7306112290129605E-01/)
       V_rule%lambda(  7,1:3) = (/ 4.1346943854935198E-01,  1.7306112290129499E-01,  4.1346943854935303E-01/)
       V_rule%lambda(  8,1:3) = (/ 5.9082801866016998E-02,  4.7045859906699100E-01,  4.7045859906699200E-01/)
       V_rule%lambda(  9,1:3) = (/ 4.7045859906699100E-01,  4.7045859906699100E-01,  5.9082801866017998E-02/)
       V_rule%lambda( 10,1:3) = (/ 4.7045859906699100E-01,  5.9082801866016998E-02,  4.7045859906699200E-01/)
       V_rule%lambda( 11,1:3) = (/ 5.1889250006095800E-01,  2.4055374996952100E-01,  2.4055374996952100E-01/)
       V_rule%lambda( 12,1:3) = (/ 2.4055374996952100E-01,  2.4055374996952100E-01,  5.1889250006095800E-01/)
       V_rule%lambda( 13,1:3) = (/ 2.4055374996952100E-01,  5.1889250006095800E-01,  2.4055374996952100E-01/)
       V_rule%lambda( 14,1:3) = (/ 7.0406841155485400E-01,  1.4796579422257300E-01,  1.4796579422257300E-01/)
       V_rule%lambda( 15,1:3) = (/ 1.4796579422257300E-01,  1.4796579422257300E-01,  7.0406841155485400E-01/)
       V_rule%lambda( 16,1:3) = (/ 1.4796579422257300E-01,  7.0406841155485400E-01,  1.4796579422257300E-01/)
       V_rule%lambda( 17,1:3) = (/ 8.4906962468505198E-01,  7.5465187657473995E-02,  7.5465187657474023E-02/)
       V_rule%lambda( 18,1:3) = (/ 7.5465187657473995E-02,  7.5465187657473995E-02,  8.4906962468505198E-01/)
       V_rule%lambda( 19,1:3) = (/ 7.5465187657473995E-02,  8.4906962468505198E-01,  7.5465187657474023E-02/)
       V_rule%lambda( 20,1:3) = (/ 9.6680719475395005E-01,  1.6596402623024999E-02,  1.6596402623024951E-02/)
       V_rule%lambda( 21,1:3) = (/ 1.6596402623024999E-02,  1.6596402623024999E-02,  9.6680719475395005E-01/)
       V_rule%lambda( 22,1:3) = (/ 1.6596402623024999E-02,  9.6680719475395005E-01,  1.6596402623024951E-02/)
       V_rule%lambda( 23,1:3) = (/ 1.0357569224525200E-01,  2.9655559657988700E-01,  5.9986871117486096E-01/)
       V_rule%lambda( 24,1:3) = (/ 2.9655559657988700E-01,  5.9986871117486096E-01,  1.0357569224525204E-01/)
       V_rule%lambda( 25,1:3) = (/ 5.9986871117486096E-01,  1.0357569224525200E-01,  2.9655559657988706E-01/)
       V_rule%lambda( 26,1:3) = (/ 2.9655559657988700E-01,  1.0357569224525200E-01,  5.9986871117486096E-01/)
       V_rule%lambda( 27,1:3) = (/ 5.9986871117486096E-01,  2.9655559657988700E-01,  1.0357569224525204E-01/)
       V_rule%lambda( 28,1:3) = (/ 1.0357569224525200E-01,  5.9986871117486096E-01,  2.9655559657988706E-01/)
       V_rule%lambda( 29,1:3) = (/ 2.0083411655416002E-02,  3.3772306340307900E-01,  6.4219352494150495E-01/)
       V_rule%lambda( 30,1:3) = (/ 3.3772306340307900E-01,  6.4219352494150495E-01,  2.0083411655416050E-02/)
       V_rule%lambda( 31,1:3) = (/ 6.4219352494150495E-01,  2.0083411655416002E-02,  3.3772306340307906E-01/)
       V_rule%lambda( 32,1:3) = (/ 3.3772306340307900E-01,  2.0083411655416002E-02,  6.4219352494150495E-01/)
       V_rule%lambda( 33,1:3) = (/ 6.4219352494150495E-01,  3.3772306340307900E-01,  2.0083411655416050E-02/)
       V_rule%lambda( 34,1:3) = (/ 2.0083411655416002E-02,  6.4219352494150495E-01,  3.3772306340307906E-01/)
       V_rule%lambda( 35,1:3) = (/-4.3410026141390001E-03,  2.0474828164281200E-01,  7.9959272097132694E-01/)
       V_rule%lambda( 36,1:3) = (/ 2.0474828164281200E-01,  7.9959272097132705E-01, -4.3410026141390556E-03/)
       V_rule%lambda( 37,1:3) = (/ 7.9959272097132705E-01, -4.3410026141390001E-03,  2.0474828164281195E-01/)
       V_rule%lambda( 38,1:3) = (/ 2.0474828164281200E-01, -4.3410026141390001E-03,  7.9959272097132694E-01/)
       V_rule%lambda( 39,1:3) = (/ 7.9959272097132705E-01,  2.0474828164281200E-01, -4.3410026141390556E-03/)
       V_rule%lambda( 40,1:3) = (/-4.3410026141390001E-03,  7.9959272097132705E-01,  2.0474828164281195E-01/)
       V_rule%lambda( 41,1:3) = (/ 4.1941786468010001E-02,  1.8935849213062300E-01,  7.6869972140136700E-01/)
       V_rule%lambda( 42,1:3) = (/ 1.8935849213062300E-01,  7.6869972140136800E-01,  4.1941786468009001E-02/)
       V_rule%lambda( 43,1:3) = (/ 7.6869972140136800E-01,  4.1941786468010001E-02,  1.8935849213062200E-01/)
       V_rule%lambda( 44,1:3) = (/ 1.8935849213062300E-01,  4.1941786468010001E-02,  7.6869972140136700E-01/)
       V_rule%lambda( 45,1:3) = (/ 7.6869972140136800E-01,  1.8935849213062300E-01,  4.1941786468009001E-02/)
       V_rule%lambda( 46,1:3) = (/ 4.1941786468010001E-02,  7.6869972140136800E-01,  1.8935849213062200E-01/)
       V_rule%lambda( 47,1:3) = (/ 1.4317320230681000E-02,  8.5283615682656994E-02,  9.0039906408666204E-01/)
       V_rule%lambda( 48,1:3) = (/ 8.5283615682656994E-02,  9.0039906408666104E-01,  1.4317320230681968E-02/)
       V_rule%lambda( 49,1:3) = (/ 9.0039906408666104E-01,  1.4317320230681000E-02,  8.5283615682657965E-02/)
       V_rule%lambda( 50,1:3) = (/ 8.5283615682656994E-02,  1.4317320230681000E-02,  9.0039906408666204E-01/)
       V_rule%lambda( 51,1:3) = (/ 9.0039906408666104E-01,  8.5283615682656994E-02,  1.4317320230681968E-02/)
       V_rule%lambda( 52,1:3) = (/ 1.4317320230681000E-02,  9.0039906408666104E-01,  8.5283615682657965E-02/)


    case( 17  )
       V_rule%Qdof =  61
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  3.3437199290803001E-02
       V_rule%weights(  2) =  5.0934154405070002E-03
       V_rule%weights(  3) =  5.0934154405070002E-03
       V_rule%weights(  4) =  5.0934154405070002E-03
       V_rule%weights(  5) =  1.4670864527638000E-02
       V_rule%weights(  6) =  1.4670864527638000E-02
       V_rule%weights(  7) =  1.4670864527638000E-02
       V_rule%weights(  8) =  2.4350878353672001E-02
       V_rule%weights(  9) =  2.4350878353672001E-02
       V_rule%weights( 10) =  2.4350878353672001E-02
       V_rule%weights( 11) =  3.1107550868968999E-02
       V_rule%weights( 12) =  3.1107550868968999E-02
       V_rule%weights( 13) =  3.1107550868968999E-02
       V_rule%weights( 14) =  3.1257111218620001E-02
       V_rule%weights( 15) =  3.1257111218620001E-02
       V_rule%weights( 16) =  3.1257111218620001E-02
       V_rule%weights( 17) =  2.4815654339665000E-02
       V_rule%weights( 18) =  2.4815654339665000E-02
       V_rule%weights( 19) =  2.4815654339665000E-02
       V_rule%weights( 20) =  1.4056073070557000E-02
       V_rule%weights( 21) =  1.4056073070557000E-02
       V_rule%weights( 22) =  1.4056073070557000E-02
       V_rule%weights( 23) =  3.1946761737789999E-03
       V_rule%weights( 24) =  3.1946761737789999E-03
       V_rule%weights( 25) =  3.1946761737789999E-03
       V_rule%weights( 26) =  8.1196553189930003E-03
       V_rule%weights( 27) =  8.1196553189930003E-03
       V_rule%weights( 28) =  8.1196553189930003E-03
       V_rule%weights( 29) =  8.1196553189930003E-03
       V_rule%weights( 30) =  8.1196553189930003E-03
       V_rule%weights( 31) =  8.1196553189930003E-03
       V_rule%weights( 32) =  2.6805742283163000E-02
       V_rule%weights( 33) =  2.6805742283163000E-02
       V_rule%weights( 34) =  2.6805742283163000E-02
       V_rule%weights( 35) =  2.6805742283163000E-02
       V_rule%weights( 36) =  2.6805742283163000E-02
       V_rule%weights( 37) =  2.6805742283163000E-02
       V_rule%weights( 38) =  1.8459993210822000E-02
       V_rule%weights( 39) =  1.8459993210822000E-02
       V_rule%weights( 40) =  1.8459993210822000E-02
       V_rule%weights( 41) =  1.8459993210822000E-02
       V_rule%weights( 42) =  1.8459993210822000E-02
       V_rule%weights( 43) =  1.8459993210822000E-02
       V_rule%weights( 44) =  8.4768685343280005E-03
       V_rule%weights( 45) =  8.4768685343280005E-03
       V_rule%weights( 46) =  8.4768685343280005E-03
       V_rule%weights( 47) =  8.4768685343280005E-03
       V_rule%weights( 48) =  8.4768685343280005E-03
       V_rule%weights( 49) =  8.4768685343280005E-03
       V_rule%weights( 50) =  1.8292796770024999E-02
       V_rule%weights( 51) =  1.8292796770024999E-02
       V_rule%weights( 52) =  1.8292796770024999E-02
       V_rule%weights( 53) =  1.8292796770024999E-02
       V_rule%weights( 54) =  1.8292796770024999E-02
       V_rule%weights( 55) =  1.8292796770024999E-02
       V_rule%weights( 56) =  6.6656320041650003E-03
       V_rule%weights( 57) =  6.6656320041650003E-03
       V_rule%weights( 58) =  6.6656320041650003E-03
       V_rule%weights( 59) =  6.6656320041650003E-03
       V_rule%weights( 60) =  6.6656320041650003E-03
       V_rule%weights( 61) =  6.6656320041650003E-03

       V_rule%lambda(  1,1:3) = (/ 3.3333333333333298E-01,  3.3333333333333298E-01,  3.3333333333333404E-01/)
       V_rule%lambda(  2,1:3) = (/ 5.6589188864520001E-03,  4.9717054055677401E-01,  4.9717054055677401E-01/)
       V_rule%lambda(  3,1:3) = (/ 4.9717054055677401E-01,  4.9717054055677401E-01,  5.6589188864519802E-03/)
       V_rule%lambda(  4,1:3) = (/ 4.9717054055677401E-01,  5.6589188864520001E-03,  4.9717054055677401E-01/)
       V_rule%lambda(  5,1:3) = (/ 3.5647354750751002E-02,  4.8217632262462501E-01,  4.8217632262462401E-01/)
       V_rule%lambda(  6,1:3) = (/ 4.8217632262462501E-01,  4.8217632262462501E-01,  3.5647354750749982E-02/)
       V_rule%lambda(  7,1:3) = (/ 4.8217632262462501E-01,  3.5647354750751002E-02,  4.8217632262462401E-01/)
       V_rule%lambda(  8,1:3) = (/ 9.9520061958437003E-02,  4.5023996902078200E-01,  4.5023996902078101E-01/)
       V_rule%lambda(  9,1:3) = (/ 4.5023996902078200E-01,  4.5023996902078200E-01,  9.9520061958435990E-02/)
       V_rule%lambda( 10,1:3) = (/ 4.5023996902078200E-01,  9.9520061958437003E-02,  4.5023996902078101E-01/)
       V_rule%lambda( 11,1:3) = (/ 1.9946752124520600E-01,  4.0026623937739703E-01,  4.0026623937739697E-01/)
       V_rule%lambda( 12,1:3) = (/ 4.0026623937739703E-01,  4.0026623937739703E-01,  1.9946752124520595E-01/)
       V_rule%lambda( 13,1:3) = (/ 4.0026623937739703E-01,  1.9946752124520600E-01,  4.0026623937739697E-01/)
       V_rule%lambda( 14,1:3) = (/ 4.9571746405809503E-01,  2.5214126797095299E-01,  2.5214126797095199E-01/)
       V_rule%lambda( 15,1:3) = (/ 2.5214126797095299E-01,  2.5214126797095299E-01,  4.9571746405809403E-01/)
       V_rule%lambda( 16,1:3) = (/ 2.5214126797095299E-01,  4.9571746405809503E-01,  2.5214126797095199E-01/)
       V_rule%lambda( 17,1:3) = (/ 6.7590599068307700E-01,  1.6204700465846100E-01,  1.6204700465846200E-01/)
       V_rule%lambda( 18,1:3) = (/ 1.6204700465846100E-01,  1.6204700465846100E-01,  6.7590599068307800E-01/)
       V_rule%lambda( 19,1:3) = (/ 1.6204700465846100E-01,  6.7590599068307700E-01,  1.6204700465846200E-01/)
       V_rule%lambda( 20,1:3) = (/ 8.4824823547850803E-01,  7.5875882260746000E-02,  7.5875882260745972E-02/)
       V_rule%lambda( 21,1:3) = (/ 7.5875882260746000E-02,  7.5875882260746000E-02,  8.4824823547850803E-01/)
       V_rule%lambda( 22,1:3) = (/ 7.5875882260746000E-02,  8.4824823547850803E-01,  7.5875882260745972E-02/)
       V_rule%lambda( 23,1:3) = (/ 9.6869054606435601E-01,  1.5654726967822000E-02,  1.5654726967821993E-02/)
       V_rule%lambda( 24,1:3) = (/ 1.5654726967822000E-02,  1.5654726967822000E-02,  9.6869054606435601E-01/)
       V_rule%lambda( 25,1:3) = (/ 1.5654726967822000E-02,  9.6869054606435601E-01,  1.5654726967821993E-02/)
       V_rule%lambda( 26,1:3) = (/ 1.0186928826919000E-02,  3.3431986736365799E-01,  6.5549320380942300E-01/)
       V_rule%lambda( 27,1:3) = (/ 3.3431986736365799E-01,  6.5549320380942300E-01,  1.0186928826919017E-02/)
       V_rule%lambda( 28,1:3) = (/ 6.5549320380942300E-01,  1.0186928826919000E-02,  3.3431986736365799E-01/)
       V_rule%lambda( 29,1:3) = (/ 3.3431986736365799E-01,  1.0186928826919000E-02,  6.5549320380942300E-01/)
       V_rule%lambda( 30,1:3) = (/ 6.5549320380942300E-01,  3.3431986736365799E-01,  1.0186928826919017E-02/)
       V_rule%lambda( 31,1:3) = (/ 1.0186928826919000E-02,  6.5549320380942300E-01,  3.3431986736365799E-01/)
       V_rule%lambda( 32,1:3) = (/ 1.3544087167103600E-01,  2.9222153779694399E-01,  5.7233759053202005E-01/)
       V_rule%lambda( 33,1:3) = (/ 2.9222153779694399E-01,  5.7233759053202005E-01,  1.3544087167103597E-01/)
       V_rule%lambda( 34,1:3) = (/ 5.7233759053202005E-01,  1.3544087167103600E-01,  2.9222153779694393E-01/)
       V_rule%lambda( 35,1:3) = (/ 2.9222153779694399E-01,  1.3544087167103600E-01,  5.7233759053202005E-01/)
       V_rule%lambda( 36,1:3) = (/ 5.7233759053202005E-01,  2.9222153779694399E-01,  1.3544087167103597E-01/)
       V_rule%lambda( 37,1:3) = (/ 1.3544087167103600E-01,  5.7233759053202005E-01,  2.9222153779694393E-01/)
       V_rule%lambda( 38,1:3) = (/ 5.4423924290583001E-02,  3.1957488542319001E-01,  6.2600119028622703E-01/)
       V_rule%lambda( 39,1:3) = (/ 3.1957488542319001E-01,  6.2600119028622803E-01,  5.4423924290581960E-02/)
       V_rule%lambda( 40,1:3) = (/ 6.2600119028622803E-01,  5.4423924290583001E-02,  3.1957488542318896E-01/)
       V_rule%lambda( 41,1:3) = (/ 3.1957488542319001E-01,  5.4423924290583001E-02,  6.2600119028622703E-01/)
       V_rule%lambda( 42,1:3) = (/ 6.2600119028622803E-01,  3.1957488542319001E-01,  5.4423924290581960E-02/)
       V_rule%lambda( 43,1:3) = (/ 5.4423924290583001E-02,  6.2600119028622803E-01,  3.1957488542318896E-01/)
       V_rule%lambda( 44,1:3) = (/ 1.2868560833637001E-02,  1.9070422419229199E-01,  7.9642721497407098E-01/)
       V_rule%lambda( 45,1:3) = (/ 1.9070422419229199E-01,  7.9642721497407098E-01,  1.2868560833637022E-02/)
       V_rule%lambda( 46,1:3) = (/ 7.9642721497407098E-01,  1.2868560833637001E-02,  1.9070422419229202E-01/)
       V_rule%lambda( 47,1:3) = (/ 1.9070422419229199E-01,  1.2868560833637001E-02,  7.9642721497407098E-01/)
       V_rule%lambda( 48,1:3) = (/ 7.9642721497407098E-01,  1.9070422419229199E-01,  1.2868560833637022E-02/)
       V_rule%lambda( 49,1:3) = (/ 1.2868560833637001E-02,  7.9642721497407098E-01,  1.9070422419229202E-01/)
       V_rule%lambda( 50,1:3) = (/ 6.7165782413524000E-02,  1.8048321164874601E-01,  7.5235100593772997E-01/)
       V_rule%lambda( 51,1:3) = (/ 1.8048321164874601E-01,  7.5235100593772897E-01,  6.7165782413525027E-02/)
       V_rule%lambda( 52,1:3) = (/ 7.5235100593772897E-01,  6.7165782413524000E-02,  1.8048321164874703E-01/)
       V_rule%lambda( 53,1:3) = (/ 1.8048321164874601E-01,  6.7165782413524000E-02,  7.5235100593772997E-01/)
       V_rule%lambda( 54,1:3) = (/ 7.5235100593772897E-01,  1.8048321164874601E-01,  6.7165782413525027E-02/)
       V_rule%lambda( 55,1:3) = (/ 6.7165782413524000E-02,  7.5235100593772897E-01,  1.8048321164874703E-01/)
       V_rule%lambda( 56,1:3) = (/ 1.4663182224828000E-02,  8.0711313679563995E-02,  9.0462550409560805E-01/)
       V_rule%lambda( 57,1:3) = (/ 8.0711313679563995E-02,  9.0462550409560805E-01,  1.4663182224827959E-02/)
       V_rule%lambda( 58,1:3) = (/ 9.0462550409560805E-01,  1.4663182224828000E-02,  8.0711313679563954E-02/)
       V_rule%lambda( 59,1:3) = (/ 8.0711313679563995E-02,  1.4663182224828000E-02,  9.0462550409560805E-01/)
       V_rule%lambda( 60,1:3) = (/ 9.0462550409560805E-01,  8.0711313679563995E-02,  1.4663182224827959E-02/)
       V_rule%lambda( 61,1:3) = (/ 1.4663182224828000E-02,  9.0462550409560805E-01,  8.0711313679563954E-02/)


    case( 18  )
       V_rule%Qdof =  70
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  3.0809939937647000E-02
       V_rule%weights(  2) =  9.0724366794039998E-03
       V_rule%weights(  3) =  9.0724366794039998E-03
       V_rule%weights(  4) =  9.0724366794039998E-03
       V_rule%weights(  5) =  1.8761316939594000E-02
       V_rule%weights(  6) =  1.8761316939594000E-02
       V_rule%weights(  7) =  1.8761316939594000E-02
       V_rule%weights(  8) =  1.9441097985477000E-02
       V_rule%weights(  9) =  1.9441097985477000E-02
       V_rule%weights( 10) =  1.9441097985477000E-02
       V_rule%weights( 11) =  2.7753948610810000E-02
       V_rule%weights( 12) =  2.7753948610810000E-02
       V_rule%weights( 13) =  2.7753948610810000E-02
       V_rule%weights( 14) =  3.2256225351457002E-02
       V_rule%weights( 15) =  3.2256225351457002E-02
       V_rule%weights( 16) =  3.2256225351457002E-02
       V_rule%weights( 17) =  2.5074032616921999E-02
       V_rule%weights( 18) =  2.5074032616921999E-02
       V_rule%weights( 19) =  2.5074032616921999E-02
       V_rule%weights( 20) =  1.5271927971832000E-02
       V_rule%weights( 21) =  1.5271927971832000E-02
       V_rule%weights( 22) =  1.5271927971832000E-02
       V_rule%weights( 23) =  6.7939220229630004E-03
       V_rule%weights( 24) =  6.7939220229630004E-03
       V_rule%weights( 25) =  6.7939220229630004E-03
       V_rule%weights( 26) = -2.2230987299200001E-03
       V_rule%weights( 27) = -2.2230987299200001E-03
       V_rule%weights( 28) = -2.2230987299200001E-03
       V_rule%weights( 29) =  6.3319140764059997E-03
       V_rule%weights( 30) =  6.3319140764059997E-03
       V_rule%weights( 31) =  6.3319140764059997E-03
       V_rule%weights( 32) =  6.3319140764059997E-03
       V_rule%weights( 33) =  6.3319140764059997E-03
       V_rule%weights( 34) =  6.3319140764059997E-03
       V_rule%weights( 35) =  2.7257538049138001E-02
       V_rule%weights( 36) =  2.7257538049138001E-02
       V_rule%weights( 37) =  2.7257538049138001E-02
       V_rule%weights( 38) =  2.7257538049138001E-02
       V_rule%weights( 39) =  2.7257538049138001E-02
       V_rule%weights( 40) =  2.7257538049138001E-02
       V_rule%weights( 41) =  1.7676785649465000E-02
       V_rule%weights( 42) =  1.7676785649465000E-02
       V_rule%weights( 43) =  1.7676785649465000E-02
       V_rule%weights( 44) =  1.7676785649465000E-02
       V_rule%weights( 45) =  1.7676785649465000E-02
       V_rule%weights( 46) =  1.7676785649465000E-02
       V_rule%weights( 47) =  1.8379484638070001E-02
       V_rule%weights( 48) =  1.8379484638070001E-02
       V_rule%weights( 49) =  1.8379484638070001E-02
       V_rule%weights( 50) =  1.8379484638070001E-02
       V_rule%weights( 51) =  1.8379484638070001E-02
       V_rule%weights( 52) =  1.8379484638070001E-02
       V_rule%weights( 53) =  8.1047328081919993E-03
       V_rule%weights( 54) =  8.1047328081919993E-03
       V_rule%weights( 55) =  8.1047328081919993E-03
       V_rule%weights( 56) =  8.1047328081919993E-03
       V_rule%weights( 57) =  8.1047328081919993E-03
       V_rule%weights( 58) =  8.1047328081919993E-03
       V_rule%weights( 59) =  7.6341290707250004E-03
       V_rule%weights( 60) =  7.6341290707250004E-03
       V_rule%weights( 61) =  7.6341290707250004E-03
       V_rule%weights( 62) =  7.6341290707250004E-03
       V_rule%weights( 63) =  7.6341290707250004E-03
       V_rule%weights( 64) =  7.6341290707250004E-03
       V_rule%weights( 65) =  4.6187660794000000E-05
       V_rule%weights( 66) =  4.6187660794000000E-05
       V_rule%weights( 67) =  4.6187660794000000E-05
       V_rule%weights( 68) =  4.6187660794000000E-05
       V_rule%weights( 69) =  4.6187660794000000E-05
       V_rule%weights( 70) =  4.6187660794000000E-05

       V_rule%lambda(  1,1:3) = (/ 3.3333333333333298E-01,  3.3333333333333298E-01,  3.3333333333333404E-01/)
       V_rule%lambda(  2,1:3) = (/ 1.3310382738157000E-02,  4.9334480863092101E-01,  4.9334480863092200E-01/)
       V_rule%lambda(  3,1:3) = (/ 4.9334480863092101E-01,  4.9334480863092101E-01,  1.3310382738157989E-02/)
       V_rule%lambda(  4,1:3) = (/ 4.9334480863092101E-01,  1.3310382738157000E-02,  4.9334480863092200E-01/)
       V_rule%lambda(  5,1:3) = (/ 6.1578811516085998E-02,  4.6921059424195699E-01,  4.6921059424195699E-01/)
       V_rule%lambda(  6,1:3) = (/ 4.6921059424195699E-01,  4.6921059424195699E-01,  6.1578811516086018E-02/)
       V_rule%lambda(  7,1:3) = (/ 4.6921059424195699E-01,  6.1578811516085998E-02,  4.6921059424195699E-01/)
       V_rule%lambda(  8,1:3) = (/ 1.2743720822598900E-01,  4.3628139588700598E-01,  4.3628139588700499E-01/)
       V_rule%lambda(  9,1:3) = (/ 4.3628139588700598E-01,  4.3628139588700598E-01,  1.2743720822598803E-01/)
       V_rule%lambda( 10,1:3) = (/ 4.3628139588700598E-01,  1.2743720822598900E-01,  4.3628139588700499E-01/)
       V_rule%lambda( 11,1:3) = (/ 2.1030765865316800E-01,  3.9484617067341599E-01,  3.9484617067341599E-01/)
       V_rule%lambda( 12,1:3) = (/ 3.9484617067341599E-01,  3.9484617067341599E-01,  2.1030765865316803E-01/)
       V_rule%lambda( 13,1:3) = (/ 3.9484617067341599E-01,  2.1030765865316800E-01,  3.9484617067341599E-01/)
       V_rule%lambda( 14,1:3) = (/ 5.0041086239368604E-01,  2.4979456880315701E-01,  2.4979456880315695E-01/)
       V_rule%lambda( 15,1:3) = (/ 2.4979456880315701E-01,  2.4979456880315701E-01,  5.0041086239368604E-01/)
       V_rule%lambda( 16,1:3) = (/ 2.4979456880315701E-01,  5.0041086239368604E-01,  2.4979456880315695E-01/)
       V_rule%lambda( 17,1:3) = (/ 6.7713561251231502E-01,  1.6143219374384299E-01,  1.6143219374384199E-01/)
       V_rule%lambda( 18,1:3) = (/ 1.6143219374384299E-01,  1.6143219374384299E-01,  6.7713561251231402E-01/)
       V_rule%lambda( 19,1:3) = (/ 1.6143219374384299E-01,  6.7713561251231502E-01,  1.6143219374384199E-01/)
       V_rule%lambda( 20,1:3) = (/ 8.4680354502925703E-01,  7.6598227485370998E-02,  7.6598227485371970E-02/)
       V_rule%lambda( 21,1:3) = (/ 7.6598227485370998E-02,  7.6598227485370998E-02,  8.4680354502925803E-01/)
       V_rule%lambda( 22,1:3) = (/ 7.6598227485370998E-02,  8.4680354502925703E-01,  7.6598227485371970E-02/)
       V_rule%lambda( 23,1:3) = (/ 9.5149512129309999E-01,  2.4252439353449999E-02,  2.4252439353450013E-02/)
       V_rule%lambda( 24,1:3) = (/ 2.4252439353449999E-02,  2.4252439353449999E-02,  9.5149512129309999E-01/)
       V_rule%lambda( 25,1:3) = (/ 2.4252439353449999E-02,  9.5149512129309999E-01,  2.4252439353450013E-02/)
       V_rule%lambda( 26,1:3) = (/ 9.1370726556607096E-01,  4.3146367216964999E-02,  4.3146367216964042E-02/)
       V_rule%lambda( 27,1:3) = (/ 4.3146367216964999E-02,  4.3146367216964999E-02,  9.1370726556606996E-01/)
       V_rule%lambda( 28,1:3) = (/ 4.3146367216964999E-02,  9.1370726556607096E-01,  4.3146367216964042E-02/)
       V_rule%lambda( 29,1:3) = (/ 8.4305362024199998E-03,  3.5891149494094399E-01,  6.3265796885663605E-01/)
       V_rule%lambda( 30,1:3) = (/ 3.5891149494094399E-01,  6.3265796885663605E-01,  8.4305362024199582E-03/)
       V_rule%lambda( 31,1:3) = (/ 6.3265796885663605E-01,  8.4305362024199998E-03,  3.5891149494094393E-01/)
       V_rule%lambda( 32,1:3) = (/ 3.5891149494094399E-01,  8.4305362024199998E-03,  6.3265796885663605E-01/)
       V_rule%lambda( 33,1:3) = (/ 6.3265796885663605E-01,  3.5891149494094399E-01,  8.4305362024199582E-03/)
       V_rule%lambda( 34,1:3) = (/ 8.4305362024199998E-03,  6.3265796885663605E-01,  3.5891149494094393E-01/)
       V_rule%lambda( 35,1:3) = (/ 1.3118655173718799E-01,  2.9440247675195702E-01,  5.7441097151085496E-01/)
       V_rule%lambda( 36,1:3) = (/ 2.9440247675195702E-01,  5.7441097151085496E-01,  1.3118655173718802E-01/)
       V_rule%lambda( 37,1:3) = (/ 5.7441097151085496E-01,  1.3118655173718799E-01,  2.9440247675195708E-01/)
       V_rule%lambda( 38,1:3) = (/ 2.9440247675195702E-01,  1.3118655173718799E-01,  5.7441097151085496E-01/)
       V_rule%lambda( 39,1:3) = (/ 5.7441097151085496E-01,  2.9440247675195702E-01,  1.3118655173718802E-01/)
       V_rule%lambda( 40,1:3) = (/ 1.3118655173718799E-01,  5.7441097151085496E-01,  2.9440247675195708E-01/)
       V_rule%lambda( 41,1:3) = (/ 5.0203151565675001E-02,  3.2501780164181399E-01,  6.2477904679251106E-01/)
       V_rule%lambda( 42,1:3) = (/ 3.2501780164181399E-01,  6.2477904679251195E-01,  5.0203151565674065E-02/)
       V_rule%lambda( 43,1:3) = (/ 6.2477904679251195E-01,  5.0203151565675001E-02,  3.2501780164181304E-01/)
       V_rule%lambda( 44,1:3) = (/ 3.2501780164181399E-01,  5.0203151565675001E-02,  6.2477904679251106E-01/)
       V_rule%lambda( 45,1:3) = (/ 6.2477904679251195E-01,  3.2501780164181399E-01,  5.0203151565674065E-02/)
       V_rule%lambda( 46,1:3) = (/ 5.0203151565675001E-02,  6.2477904679251195E-01,  3.2501780164181304E-01/)
       V_rule%lambda( 47,1:3) = (/ 6.6329263810916000E-02,  1.8473755966604599E-01,  7.4893317652303804E-01/)
       V_rule%lambda( 48,1:3) = (/ 1.8473755966604599E-01,  7.4893317652303704E-01,  6.6329263810916972E-02/)
       V_rule%lambda( 49,1:3) = (/ 7.4893317652303704E-01,  6.6329263810916000E-02,  1.8473755966604696E-01/)
       V_rule%lambda( 50,1:3) = (/ 1.8473755966604599E-01,  6.6329263810916000E-02,  7.4893317652303804E-01/)
       V_rule%lambda( 51,1:3) = (/ 7.4893317652303704E-01,  1.8473755966604599E-01,  6.6329263810916972E-02/)
       V_rule%lambda( 52,1:3) = (/ 6.6329263810916000E-02,  7.4893317652303704E-01,  1.8473755966604696E-01/)
       V_rule%lambda( 53,1:3) = (/ 1.1996194566236001E-02,  2.1879680001332100E-01,  7.6920700542044296E-01/)
       V_rule%lambda( 54,1:3) = (/ 2.1879680001332100E-01,  7.6920700542044296E-01,  1.1996194566236046E-02/)
       V_rule%lambda( 55,1:3) = (/ 7.6920700542044296E-01,  1.1996194566236001E-02,  2.1879680001332105E-01/)
       V_rule%lambda( 56,1:3) = (/ 2.1879680001332100E-01,  1.1996194566236001E-02,  7.6920700542044296E-01/)
       V_rule%lambda( 57,1:3) = (/ 7.6920700542044296E-01,  2.1879680001332100E-01,  1.1996194566236046E-02/)
       V_rule%lambda( 58,1:3) = (/ 1.1996194566236001E-02,  7.6920700542044296E-01,  2.1879680001332105E-01/)
       V_rule%lambda( 59,1:3) = (/ 1.4858100590125000E-02,  1.0117959713640801E-01,  8.8396230227346695E-01/)
       V_rule%lambda( 60,1:3) = (/ 1.0117959713640801E-01,  8.8396230227346695E-01,  1.4858100590125045E-02/)
       V_rule%lambda( 61,1:3) = (/ 8.8396230227346695E-01,  1.4858100590125000E-02,  1.0117959713640805E-01/)
       V_rule%lambda( 62,1:3) = (/ 1.0117959713640801E-01,  1.4858100590125000E-02,  8.8396230227346695E-01/)
       V_rule%lambda( 63,1:3) = (/ 8.8396230227346695E-01,  1.0117959713640801E-01,  1.4858100590125045E-02/)
       V_rule%lambda( 64,1:3) = (/ 1.4858100590125000E-02,  8.8396230227346695E-01,  1.0117959713640805E-01/)
       V_rule%lambda( 65,1:3) = (/-3.5222015287949000E-02,  2.0874755282586002E-02,  1.0143472600053629E+00/)
       V_rule%lambda( 66,1:3) = (/ 2.0874755282586002E-02,  1.0143472600053629E+00, -3.5222015287948910E-02/)
       V_rule%lambda( 67,1:3) = (/ 1.0143472600053629E+00, -3.5222015287949000E-02,  2.0874755282586095E-02/)
       V_rule%lambda( 68,1:3) = (/ 2.0874755282586002E-02, -3.5222015287949000E-02,  1.0143472600053629E+00/)
       V_rule%lambda( 69,1:3) = (/ 1.0143472600053629E+00,  2.0874755282586002E-02, -3.5222015287948910E-02/)
       V_rule%lambda( 70,1:3) = (/-3.5222015287949000E-02,  1.0143472600053629E+00,  2.0874755282586095E-02/)


    case( 19  )
       V_rule%Qdof =  73
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  3.2906331388918998E-02
       V_rule%weights(  2) =  1.0330731891272000E-02
       V_rule%weights(  3) =  1.0330731891272000E-02
       V_rule%weights(  4) =  1.0330731891272000E-02
       V_rule%weights(  5) =  2.2387247263016000E-02
       V_rule%weights(  6) =  2.2387247263016000E-02
       V_rule%weights(  7) =  2.2387247263016000E-02
       V_rule%weights(  8) =  3.0266125869467999E-02
       V_rule%weights(  9) =  3.0266125869467999E-02
       V_rule%weights( 10) =  3.0266125869467999E-02
       V_rule%weights( 11) =  3.0490967802198000E-02
       V_rule%weights( 12) =  3.0490967802198000E-02
       V_rule%weights( 13) =  3.0490967802198000E-02
       V_rule%weights( 14) =  2.4159212741641001E-02
       V_rule%weights( 15) =  2.4159212741641001E-02
       V_rule%weights( 16) =  2.4159212741641001E-02
       V_rule%weights( 17) =  1.6050803586800999E-02
       V_rule%weights( 18) =  1.6050803586800999E-02
       V_rule%weights( 19) =  1.6050803586800999E-02
       V_rule%weights( 20) =  8.0845802617839999E-03
       V_rule%weights( 21) =  8.0845802617839999E-03
       V_rule%weights( 22) =  8.0845802617839999E-03
       V_rule%weights( 23) =  2.0793620274849999E-03
       V_rule%weights( 24) =  2.0793620274849999E-03
       V_rule%weights( 25) =  2.0793620274849999E-03
       V_rule%weights( 26) =  3.8848769049810001E-03
       V_rule%weights( 27) =  3.8848769049810001E-03
       V_rule%weights( 28) =  3.8848769049810001E-03
       V_rule%weights( 29) =  3.8848769049810001E-03
       V_rule%weights( 30) =  3.8848769049810001E-03
       V_rule%weights( 31) =  3.8848769049810001E-03
       V_rule%weights( 32) =  2.5574160612022001E-02
       V_rule%weights( 33) =  2.5574160612022001E-02
       V_rule%weights( 34) =  2.5574160612022001E-02
       V_rule%weights( 35) =  2.5574160612022001E-02
       V_rule%weights( 36) =  2.5574160612022001E-02
       V_rule%weights( 37) =  2.5574160612022001E-02
       V_rule%weights( 38) =  8.8809035733379994E-03
       V_rule%weights( 39) =  8.8809035733379994E-03
       V_rule%weights( 40) =  8.8809035733379994E-03
       V_rule%weights( 41) =  8.8809035733379994E-03
       V_rule%weights( 42) =  8.8809035733379994E-03
       V_rule%weights( 43) =  8.8809035733379994E-03
       V_rule%weights( 44) =  1.6124546761731001E-02
       V_rule%weights( 45) =  1.6124546761731001E-02
       V_rule%weights( 46) =  1.6124546761731001E-02
       V_rule%weights( 47) =  1.6124546761731001E-02
       V_rule%weights( 48) =  1.6124546761731001E-02
       V_rule%weights( 49) =  1.6124546761731001E-02
       V_rule%weights( 50) =  2.4919418174909999E-03
       V_rule%weights( 51) =  2.4919418174909999E-03
       V_rule%weights( 52) =  2.4919418174909999E-03
       V_rule%weights( 53) =  2.4919418174909999E-03
       V_rule%weights( 54) =  2.4919418174909999E-03
       V_rule%weights( 55) =  2.4919418174909999E-03
       V_rule%weights( 56) =  1.8242840118951002E-02
       V_rule%weights( 57) =  1.8242840118951002E-02
       V_rule%weights( 58) =  1.8242840118951002E-02
       V_rule%weights( 59) =  1.8242840118951002E-02
       V_rule%weights( 60) =  1.8242840118951002E-02
       V_rule%weights( 61) =  1.8242840118951002E-02
       V_rule%weights( 62) =  1.0258563736198999E-02
       V_rule%weights( 63) =  1.0258563736198999E-02
       V_rule%weights( 64) =  1.0258563736198999E-02
       V_rule%weights( 65) =  1.0258563736198999E-02
       V_rule%weights( 66) =  1.0258563736198999E-02
       V_rule%weights( 67) =  1.0258563736198999E-02
       V_rule%weights( 68) =  3.7999288553020000E-03
       V_rule%weights( 69) =  3.7999288553020000E-03
       V_rule%weights( 70) =  3.7999288553020000E-03
       V_rule%weights( 71) =  3.7999288553020000E-03
       V_rule%weights( 72) =  3.7999288553020000E-03
       V_rule%weights( 73) =  3.7999288553020000E-03

       V_rule%lambda(  1,1:3) = (/ 3.3333333333333298E-01,  3.3333333333333298E-01,  3.3333333333333404E-01/)
       V_rule%lambda(  2,1:3) = (/ 2.0780025853987000E-02,  4.8960998707300601E-01,  4.8960998707300696E-01/)
       V_rule%lambda(  3,1:3) = (/ 4.8960998707300601E-01,  4.8960998707300601E-01,  2.0780025853987971E-02/)
       V_rule%lambda(  4,1:3) = (/ 4.8960998707300601E-01,  2.0780025853987000E-02,  4.8960998707300696E-01/)
       V_rule%lambda(  5,1:3) = (/ 9.0926214604215003E-02,  4.5453689269789299E-01,  4.5453689269789199E-01/)
       V_rule%lambda(  6,1:3) = (/ 4.5453689269789299E-01,  4.5453689269789299E-01,  9.0926214604214017E-02/)
       V_rule%lambda(  7,1:3) = (/ 4.5453689269789299E-01,  9.0926214604215003E-02,  4.5453689269789199E-01/)
       V_rule%lambda(  8,1:3) = (/ 1.9716663870113799E-01,  4.0141668064943098E-01,  4.0141668064943103E-01/)
       V_rule%lambda(  9,1:3) = (/ 4.0141668064943098E-01,  4.0141668064943098E-01,  1.9716663870113804E-01/)
       V_rule%lambda( 10,1:3) = (/ 4.0141668064943098E-01,  1.9716663870113799E-01,  4.0141668064943103E-01/)
       V_rule%lambda( 11,1:3) = (/ 4.8889669119380502E-01,  2.5555165440309802E-01,  2.5555165440309696E-01/)
       V_rule%lambda( 12,1:3) = (/ 2.5555165440309802E-01,  2.5555165440309802E-01,  4.8889669119380397E-01/)
       V_rule%lambda( 13,1:3) = (/ 2.5555165440309802E-01,  4.8889669119380502E-01,  2.5555165440309696E-01/)
       V_rule%lambda( 14,1:3) = (/ 6.4584411569574096E-01,  1.7707794215212999E-01,  1.7707794215212905E-01/)
       V_rule%lambda( 15,1:3) = (/ 1.7707794215212999E-01,  1.7707794215212999E-01,  6.4584411569573996E-01/)
       V_rule%lambda( 16,1:3) = (/ 1.7707794215212999E-01,  6.4584411569574096E-01,  1.7707794215212905E-01/)
       V_rule%lambda( 17,1:3) = (/ 7.7987789354409598E-01,  1.1006105322795200E-01,  1.1006105322795202E-01/)
       V_rule%lambda( 18,1:3) = (/ 1.1006105322795200E-01,  1.1006105322795200E-01,  7.7987789354409598E-01/)
       V_rule%lambda( 19,1:3) = (/ 1.1006105322795200E-01,  7.7987789354409598E-01,  1.1006105322795202E-01/)
       V_rule%lambda( 20,1:3) = (/ 8.8894275149632096E-01,  5.5528624251839999E-02,  5.5528624251839041E-02/)
       V_rule%lambda( 21,1:3) = (/ 5.5528624251839999E-02,  5.5528624251839999E-02,  8.8894275149631996E-01/)
       V_rule%lambda( 22,1:3) = (/ 5.5528624251839999E-02,  8.8894275149632096E-01,  5.5528624251839041E-02/)
       V_rule%lambda( 23,1:3) = (/ 9.7475627244554297E-01,  1.2621863777229000E-02,  1.2621863777228029E-02/)
       V_rule%lambda( 24,1:3) = (/ 1.2621863777229000E-02,  1.2621863777229000E-02,  9.7475627244554197E-01/)
       V_rule%lambda( 25,1:3) = (/ 1.2621863777229000E-02,  9.7475627244554297E-01,  1.2621863777228029E-02/)
       V_rule%lambda( 26,1:3) = (/ 3.6114178484120000E-03,  3.9575478735694303E-01,  6.0063379479464496E-01/)
       V_rule%lambda( 27,1:3) = (/ 3.9575478735694303E-01,  6.0063379479464496E-01,  3.6114178484120130E-03/)
       V_rule%lambda( 28,1:3) = (/ 6.0063379479464496E-01,  3.6114178484120000E-03,  3.9575478735694303E-01/)
       V_rule%lambda( 29,1:3) = (/ 3.9575478735694303E-01,  3.6114178484120000E-03,  6.0063379479464496E-01/)
       V_rule%lambda( 30,1:3) = (/ 6.0063379479464496E-01,  3.9575478735694303E-01,  3.6114178484120130E-03/)
       V_rule%lambda( 31,1:3) = (/ 3.6114178484120000E-03,  6.0063379479464496E-01,  3.9575478735694303E-01/)
       V_rule%lambda( 32,1:3) = (/ 1.3446675453078000E-01,  3.0792998388043602E-01,  5.5760326158878404E-01/)
       V_rule%lambda( 33,1:3) = (/ 3.0792998388043602E-01,  5.5760326158878404E-01,  1.3446675453077994E-01/)
       V_rule%lambda( 34,1:3) = (/ 5.5760326158878404E-01,  1.3446675453078000E-01,  3.0792998388043596E-01/)
       V_rule%lambda( 35,1:3) = (/ 3.0792998388043602E-01,  1.3446675453078000E-01,  5.5760326158878404E-01/)
       V_rule%lambda( 36,1:3) = (/ 5.5760326158878404E-01,  3.0792998388043602E-01,  1.3446675453077994E-01/)
       V_rule%lambda( 37,1:3) = (/ 1.3446675453078000E-01,  5.5760326158878404E-01,  3.0792998388043596E-01/)
       V_rule%lambda( 38,1:3) = (/ 1.4446025776115000E-02,  2.6456694840652001E-01,  7.2098702581736496E-01/)
       V_rule%lambda( 39,1:3) = (/ 2.6456694840652001E-01,  7.2098702581736496E-01,  1.4446025776115035E-02/)
       V_rule%lambda( 40,1:3) = (/ 7.2098702581736496E-01,  1.4446025776115000E-02,  2.6456694840652006E-01/)
       V_rule%lambda( 41,1:3) = (/ 2.6456694840652001E-01,  1.4446025776115000E-02,  7.2098702581736496E-01/)
       V_rule%lambda( 42,1:3) = (/ 7.2098702581736496E-01,  2.6456694840652001E-01,  1.4446025776115035E-02/)
       V_rule%lambda( 43,1:3) = (/ 1.4446025776115000E-02,  7.2098702581736496E-01,  2.6456694840652006E-01/)
       V_rule%lambda( 44,1:3) = (/ 4.6933578838178003E-02,  3.5853935220595101E-01,  5.9452706895587104E-01/)
       V_rule%lambda( 45,1:3) = (/ 3.5853935220595101E-01,  5.9452706895587104E-01,  4.6933578838177947E-02/)
       V_rule%lambda( 46,1:3) = (/ 5.9452706895587104E-01,  4.6933578838178003E-02,  3.5853935220595096E-01/)
       V_rule%lambda( 47,1:3) = (/ 3.5853935220595101E-01,  4.6933578838178003E-02,  5.9452706895587104E-01/)
       V_rule%lambda( 48,1:3) = (/ 5.9452706895587104E-01,  3.5853935220595101E-01,  4.6933578838177947E-02/)
       V_rule%lambda( 49,1:3) = (/ 4.6933578838178003E-02,  5.9452706895587104E-01,  3.5853935220595096E-01/)
       V_rule%lambda( 50,1:3) = (/ 2.8611203505669999E-03,  1.5780740596859499E-01,  8.3933147368083805E-01/)
       V_rule%lambda( 51,1:3) = (/ 1.5780740596859499E-01,  8.3933147368083905E-01,  2.8611203505659599E-03/)
       V_rule%lambda( 52,1:3) = (/ 8.3933147368083905E-01,  2.8611203505669999E-03,  1.5780740596859397E-01/)
       V_rule%lambda( 53,1:3) = (/ 1.5780740596859499E-01,  2.8611203505669999E-03,  8.3933147368083805E-01/)
       V_rule%lambda( 54,1:3) = (/ 8.3933147368083905E-01,  1.5780740596859499E-01,  2.8611203505659599E-03/)
       V_rule%lambda( 55,1:3) = (/ 2.8611203505669999E-03,  8.3933147368083905E-01,  1.5780740596859397E-01/)
       V_rule%lambda( 56,1:3) = (/ 2.2386142409791601E-01,  7.5050596975911002E-02,  7.0108797892617303E-01/)
       V_rule%lambda( 57,1:3) = (/ 7.5050596975911002E-02,  7.0108797892617303E-01,  2.2386142409791598E-01/)
       V_rule%lambda( 58,1:3) = (/ 7.0108797892617303E-01,  2.2386142409791601E-01,  7.5050596975910960E-02/)
       V_rule%lambda( 59,1:3) = (/ 7.5050596975911002E-02,  2.2386142409791601E-01,  7.0108797892617303E-01/)
       V_rule%lambda( 60,1:3) = (/ 7.0108797892617303E-01,  7.5050596975911002E-02,  2.2386142409791598E-01/)
       V_rule%lambda( 61,1:3) = (/ 2.2386142409791601E-01,  7.0108797892617303E-01,  7.5050596975910960E-02/)
       V_rule%lambda( 62,1:3) = (/ 3.4647074816760000E-02,  1.4242160111338301E-01,  8.2293132406985703E-01/)
       V_rule%lambda( 63,1:3) = (/ 1.4242160111338301E-01,  8.2293132406985703E-01,  3.4647074816759965E-02/)
       V_rule%lambda( 64,1:3) = (/ 8.2293132406985703E-01,  3.4647074816760000E-02,  1.4242160111338298E-01/)
       V_rule%lambda( 65,1:3) = (/ 1.4242160111338301E-01,  3.4647074816760000E-02,  8.2293132406985703E-01/)
       V_rule%lambda( 66,1:3) = (/ 8.2293132406985703E-01,  1.4242160111338301E-01,  3.4647074816759965E-02/)
       V_rule%lambda( 67,1:3) = (/ 3.4647074816760000E-02,  8.2293132406985703E-01,  1.4242160111338298E-01/)
       V_rule%lambda( 68,1:3) = (/ 1.0161119296278000E-02,  6.5494628082938003E-02,  9.2434425262078401E-01/)
       V_rule%lambda( 69,1:3) = (/ 6.5494628082938003E-02,  9.2434425262078401E-01,  1.0161119296277984E-02/)
       V_rule%lambda( 70,1:3) = (/ 9.2434425262078401E-01,  1.0161119296278000E-02,  6.5494628082937989E-02/)
       V_rule%lambda( 71,1:3) = (/ 6.5494628082938003E-02,  1.0161119296278000E-02,  9.2434425262078401E-01/)
       V_rule%lambda( 72,1:3) = (/ 9.2434425262078401E-01,  6.5494628082938003E-02,  1.0161119296277984E-02/)
       V_rule%lambda( 73,1:3) = (/ 1.0161119296278000E-02,  9.2434425262078401E-01,  6.5494628082937989E-02/)


    case( 20  )
       V_rule%Qdof =  79
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))

       V_rule%weights(  1) =  3.3057055541623998E-02
       V_rule%weights(  2) =  8.6701918566299996E-04
       V_rule%weights(  3) =  8.6701918566299996E-04
       V_rule%weights(  4) =  8.6701918566299996E-04
       V_rule%weights(  5) =  1.1660052716448000E-02
       V_rule%weights(  6) =  1.1660052716448000E-02
       V_rule%weights(  7) =  1.1660052716448000E-02
       V_rule%weights(  8) =  2.2876936356421001E-02
       V_rule%weights(  9) =  2.2876936356421001E-02
       V_rule%weights( 10) =  2.2876936356421001E-02
       V_rule%weights( 11) =  3.0448982673937999E-02
       V_rule%weights( 12) =  3.0448982673937999E-02
       V_rule%weights( 13) =  3.0448982673937999E-02
       V_rule%weights( 14) =  3.0624891725355000E-02
       V_rule%weights( 15) =  3.0624891725355000E-02
       V_rule%weights( 16) =  3.0624891725355000E-02
       V_rule%weights( 17) =  2.4368057676800000E-02
       V_rule%weights( 18) =  2.4368057676800000E-02
       V_rule%weights( 19) =  2.4368057676800000E-02
       V_rule%weights( 20) =  1.5997432032024000E-02
       V_rule%weights( 21) =  1.5997432032024000E-02
       V_rule%weights( 22) =  1.5997432032024000E-02
       V_rule%weights( 23) =  7.6983018156020003E-03
       V_rule%weights( 24) =  7.6983018156020003E-03
       V_rule%weights( 25) =  7.6983018156020003E-03
       V_rule%weights( 26) = -6.3206049748799995E-04
       V_rule%weights( 27) = -6.3206049748799995E-04
       V_rule%weights( 28) = -6.3206049748799995E-04
       V_rule%weights( 29) =  1.7511343011929999E-03
       V_rule%weights( 30) =  1.7511343011929999E-03
       V_rule%weights( 31) =  1.7511343011929999E-03
       V_rule%weights( 32) =  1.6465839189576000E-02
       V_rule%weights( 33) =  1.6465839189576000E-02
       V_rule%weights( 34) =  1.6465839189576000E-02
       V_rule%weights( 35) =  1.6465839189576000E-02
       V_rule%weights( 36) =  1.6465839189576000E-02
       V_rule%weights( 37) =  1.6465839189576000E-02
       V_rule%weights( 38) =  4.8390335404850000E-03
       V_rule%weights( 39) =  4.8390335404850000E-03
       V_rule%weights( 40) =  4.8390335404850000E-03
       V_rule%weights( 41) =  4.8390335404850000E-03
       V_rule%weights( 42) =  4.8390335404850000E-03
       V_rule%weights( 43) =  4.8390335404850000E-03
       V_rule%weights( 44) =  2.5804906534650000E-02
       V_rule%weights( 45) =  2.5804906534650000E-02
       V_rule%weights( 46) =  2.5804906534650000E-02
       V_rule%weights( 47) =  2.5804906534650000E-02
       V_rule%weights( 48) =  2.5804906534650000E-02
       V_rule%weights( 49) =  2.5804906534650000E-02
       V_rule%weights( 50) =  8.4710910544410004E-03
       V_rule%weights( 51) =  8.4710910544410004E-03
       V_rule%weights( 52) =  8.4710910544410004E-03
       V_rule%weights( 53) =  8.4710910544410004E-03
       V_rule%weights( 54) =  8.4710910544410004E-03
       V_rule%weights( 55) =  8.4710910544410004E-03
       V_rule%weights( 56) =  1.8354914106279999E-02
       V_rule%weights( 57) =  1.8354914106279999E-02
       V_rule%weights( 58) =  1.8354914106279999E-02
       V_rule%weights( 59) =  1.8354914106279999E-02
       V_rule%weights( 60) =  1.8354914106279999E-02
       V_rule%weights( 61) =  1.8354914106279999E-02
       V_rule%weights( 62) =  7.0440467790799997E-04
       V_rule%weights( 63) =  7.0440467790799997E-04
       V_rule%weights( 64) =  7.0440467790799997E-04
       V_rule%weights( 65) =  7.0440467790799997E-04
       V_rule%weights( 66) =  7.0440467790799997E-04
       V_rule%weights( 67) =  7.0440467790799997E-04
       V_rule%weights( 68) =  1.0112684927462000E-02
       V_rule%weights( 69) =  1.0112684927462000E-02
       V_rule%weights( 70) =  1.0112684927462000E-02
       V_rule%weights( 71) =  1.0112684927462000E-02
       V_rule%weights( 72) =  1.0112684927462000E-02
       V_rule%weights( 73) =  1.0112684927462000E-02
       V_rule%weights( 74) =  3.5739093859499999E-03
       V_rule%weights( 75) =  3.5739093859499999E-03
       V_rule%weights( 76) =  3.5739093859499999E-03
       V_rule%weights( 77) =  3.5739093859499999E-03
       V_rule%weights( 78) =  3.5739093859499999E-03
       V_rule%weights( 79) =  3.5739093859499999E-03

       V_rule%lambda(  1,1:3) = (/ 3.3333333333333298E-01,  3.3333333333333298E-01,  3.3333333333333404E-01/)
       V_rule%lambda(  2,1:3) = (/-1.9009287044000000E-03,  5.0095046435220003E-01,  5.0095046435219992E-01/)
       V_rule%lambda(  3,1:3) = (/ 5.0095046435220003E-01,  5.0095046435220003E-01, -1.9009287044000622E-03/)
       V_rule%lambda(  4,1:3) = (/ 5.0095046435220003E-01, -1.9009287044000000E-03,  5.0095046435219992E-01/)
       V_rule%lambda(  5,1:3) = (/ 2.3574084130542999E-02,  4.8821295793472902E-01,  4.8821295793472796E-01/)
       V_rule%lambda(  6,1:3) = (/ 4.8821295793472902E-01,  4.8821295793472902E-01,  2.3574084130541961E-02/)
       V_rule%lambda(  7,1:3) = (/ 4.8821295793472902E-01,  2.3574084130542999E-02,  4.8821295793472796E-01/)
       V_rule%lambda(  8,1:3) = (/ 8.9726636099435000E-02,  4.5513668195028301E-01,  4.5513668195028201E-01/)
       V_rule%lambda(  9,1:3) = (/ 4.5513668195028301E-01,  4.5513668195028301E-01,  8.9726636099433987E-02/)
       V_rule%lambda( 10,1:3) = (/ 4.5513668195028301E-01,  8.9726636099435000E-02,  4.5513668195028201E-01/)
       V_rule%lambda( 11,1:3) = (/ 1.9600748136342100E-01,  4.0199625931828897E-01,  4.0199625931829003E-01/)
       V_rule%lambda( 12,1:3) = (/ 4.0199625931828897E-01,  4.0199625931828897E-01,  1.9600748136342205E-01/)
       V_rule%lambda( 13,1:3) = (/ 4.0199625931828897E-01,  1.9600748136342100E-01,  4.0199625931829003E-01/)
       V_rule%lambda( 14,1:3) = (/ 4.8821418048115700E-01,  2.5589290975942097E-01,  2.5589290975942203E-01/)
       V_rule%lambda( 15,1:3) = (/ 2.5589290975942097E-01,  2.5589290975942097E-01,  4.8821418048115806E-01/)
       V_rule%lambda( 16,1:3) = (/ 2.5589290975942097E-01,  4.8821418048115700E-01,  2.5589290975942203E-01/)
       V_rule%lambda( 17,1:3) = (/ 6.4702348800978804E-01,  1.7648825599510601E-01,  1.7648825599510595E-01/)
       V_rule%lambda( 18,1:3) = (/ 1.7648825599510601E-01,  1.7648825599510601E-01,  6.4702348800978804E-01/)
       V_rule%lambda( 19,1:3) = (/ 1.7648825599510601E-01,  6.4702348800978804E-01,  1.7648825599510595E-01/)
       V_rule%lambda( 20,1:3) = (/ 7.9165828932648297E-01,  1.0417085533675800E-01,  1.0417085533675903E-01/)
       V_rule%lambda( 21,1:3) = (/ 1.0417085533675800E-01,  1.0417085533675800E-01,  7.9165828932648397E-01/)
       V_rule%lambda( 22,1:3) = (/ 1.0417085533675800E-01,  7.9165828932648297E-01,  1.0417085533675903E-01/)
       V_rule%lambda( 23,1:3) = (/ 8.9386207231813997E-01,  5.3068963840930003E-02,  5.3068963840930031E-02/)
       V_rule%lambda( 24,1:3) = (/ 5.3068963840930003E-02,  5.3068963840930003E-02,  8.9386207231813997E-01/)
       V_rule%lambda( 25,1:3) = (/ 5.3068963840930003E-02,  8.9386207231813997E-01,  5.3068963840930031E-02/)
       V_rule%lambda( 26,1:3) = (/ 9.1676256960794200E-01,  4.1618715196028999E-02,  4.1618715196028999E-02/)
       V_rule%lambda( 27,1:3) = (/ 4.1618715196028999E-02,  4.1618715196028999E-02,  9.1676256960794200E-01/)
       V_rule%lambda( 28,1:3) = (/ 4.1618715196028999E-02,  9.1676256960794200E-01,  4.1618715196028999E-02/)
       V_rule%lambda( 29,1:3) = (/ 9.7683615718635597E-01,  1.1581921406822000E-02,  1.1581921406822031E-02/)
       V_rule%lambda( 30,1:3) = (/ 1.1581921406822000E-02,  1.1581921406822000E-02,  9.7683615718635597E-01/)
       V_rule%lambda( 31,1:3) = (/ 1.1581921406822000E-02,  9.7683615718635597E-01,  1.1581921406822031E-02/)
       V_rule%lambda( 32,1:3) = (/ 4.8741583664839001E-02,  3.4485577022900099E-01,  6.0640264610616001E-01/)
       V_rule%lambda( 33,1:3) = (/ 3.4485577022900099E-01,  6.0640264610616001E-01,  4.8741583664838994E-02/)
       V_rule%lambda( 34,1:3) = (/ 6.0640264610616001E-01,  4.8741583664839001E-02,  3.4485577022900099E-01/)
       V_rule%lambda( 35,1:3) = (/ 3.4485577022900099E-01,  4.8741583664839001E-02,  6.0640264610616001E-01/)
       V_rule%lambda( 36,1:3) = (/ 6.0640264610616001E-01,  3.4485577022900099E-01,  4.8741583664838994E-02/)
       V_rule%lambda( 37,1:3) = (/ 4.8741583664839001E-02,  6.0640264610616001E-01,  3.4485577022900099E-01/)
       V_rule%lambda( 38,1:3) = (/ 6.3141159486049996E-03,  3.7784326959485398E-01,  6.1584261445654098E-01/)
       V_rule%lambda( 39,1:3) = (/ 3.7784326959485398E-01,  6.1584261445654098E-01,  6.3141159486050369E-03/)
       V_rule%lambda( 40,1:3) = (/ 6.1584261445654098E-01,  6.3141159486049996E-03,  3.7784326959485404E-01/)
       V_rule%lambda( 41,1:3) = (/ 3.7784326959485398E-01,  6.3141159486049996E-03,  6.1584261445654098E-01/)
       V_rule%lambda( 42,1:3) = (/ 6.1584261445654098E-01,  3.7784326959485398E-01,  6.3141159486050369E-03/)
       V_rule%lambda( 43,1:3) = (/ 6.3141159486049996E-03,  6.1584261445654098E-01,  3.7784326959485404E-01/)
       V_rule%lambda( 44,1:3) = (/ 1.3431652054734800E-01,  3.0663547906235700E-01,  5.5904800039029501E-01/)
       V_rule%lambda( 45,1:3) = (/ 3.0663547906235700E-01,  5.5904800039029501E-01,  1.3431652054734800E-01/)
       V_rule%lambda( 46,1:3) = (/ 5.5904800039029501E-01,  1.3431652054734800E-01,  3.0663547906235700E-01/)
       V_rule%lambda( 47,1:3) = (/ 3.0663547906235700E-01,  1.3431652054734800E-01,  5.5904800039029501E-01/)
       V_rule%lambda( 48,1:3) = (/ 5.5904800039029501E-01,  3.0663547906235700E-01,  1.3431652054734800E-01/)
       V_rule%lambda( 49,1:3) = (/ 1.3431652054734800E-01,  5.5904800039029501E-01,  3.0663547906235700E-01/)
       V_rule%lambda( 50,1:3) = (/ 1.3973893962392001E-02,  2.4941936277474200E-01,  7.3660674326286601E-01/)
       V_rule%lambda( 51,1:3) = (/ 2.4941936277474200E-01,  7.3660674326286601E-01,  1.3973893962391987E-02/)
       V_rule%lambda( 52,1:3) = (/ 7.3660674326286601E-01,  1.3973893962392001E-02,  2.4941936277474197E-01/)
       V_rule%lambda( 53,1:3) = (/ 2.4941936277474200E-01,  1.3973893962392001E-02,  7.3660674326286601E-01/)
       V_rule%lambda( 54,1:3) = (/ 7.3660674326286601E-01,  2.4941936277474200E-01,  1.3973893962391987E-02/)
       V_rule%lambda( 55,1:3) = (/ 1.3973893962392001E-02,  7.3660674326286601E-01,  2.4941936277474197E-01/)
       V_rule%lambda( 56,1:3) = (/ 7.5549132909764005E-02,  2.1277572480280199E-01,  7.1167514228743400E-01/)
       V_rule%lambda( 57,1:3) = (/ 2.1277572480280199E-01,  7.1167514228743400E-01,  7.5549132909764005E-02/)
       V_rule%lambda( 58,1:3) = (/ 7.1167514228743400E-01,  7.5549132909764005E-02,  2.1277572480280199E-01/)
       V_rule%lambda( 59,1:3) = (/ 2.1277572480280199E-01,  7.5549132909764005E-02,  7.1167514228743400E-01/)
       V_rule%lambda( 60,1:3) = (/ 7.1167514228743400E-01,  2.1277572480280199E-01,  7.5549132909764005E-02/)
       V_rule%lambda( 61,1:3) = (/ 7.5549132909764005E-02,  7.1167514228743400E-01,  2.1277572480280199E-01/)
       V_rule%lambda( 62,1:3) = (/-8.3681532082269996E-03,  1.4696543605323900E-01,  8.6140271715498795E-01/)
       V_rule%lambda( 63,1:3) = (/ 1.4696543605323900E-01,  8.6140271715498695E-01, -8.3681532082259535E-03/)
       V_rule%lambda( 64,1:3) = (/ 8.6140271715498695E-01, -8.3681532082269996E-03,  1.4696543605324006E-01/)
       V_rule%lambda( 65,1:3) = (/ 1.4696543605323900E-01, -8.3681532082269996E-03,  8.6140271715498795E-01/)
       V_rule%lambda( 66,1:3) = (/ 8.6140271715498695E-01,  1.4696543605323900E-01, -8.3681532082259535E-03/)
       V_rule%lambda( 67,1:3) = (/-8.3681532082269996E-03,  8.6140271715498695E-01,  1.4696543605324006E-01/)
       V_rule%lambda( 68,1:3) = (/ 2.6686063258714001E-02,  1.3772697882892301E-01,  8.3558695791236304E-01/)
       V_rule%lambda( 69,1:3) = (/ 1.3772697882892301E-01,  8.3558695791236304E-01,  2.6686063258713949E-02/)
       V_rule%lambda( 70,1:3) = (/ 8.3558695791236304E-01,  2.6686063258714001E-02,  1.3772697882892296E-01/)
       V_rule%lambda( 71,1:3) = (/ 1.3772697882892301E-01,  2.6686063258714001E-02,  8.3558695791236304E-01/)
       V_rule%lambda( 72,1:3) = (/ 8.3558695791236304E-01,  1.3772697882892301E-01,  2.6686063258713949E-02/)
       V_rule%lambda( 73,1:3) = (/ 2.6686063258714001E-02,  8.3558695791236304E-01,  1.3772697882892296E-01/)
       V_rule%lambda( 74,1:3) = (/ 1.0547719294141000E-02,  5.9696109149006998E-02,  9.2975617155685197E-01/)
       V_rule%lambda( 75,1:3) = (/ 5.9696109149006998E-02,  9.2975617155685297E-01,  1.0547719294140029E-02/)
       V_rule%lambda( 76,1:3) = (/ 9.2975617155685297E-01,  1.0547719294141000E-02,  5.9696109149006027E-02/)
       V_rule%lambda( 77,1:3) = (/ 5.9696109149006998E-02,  1.0547719294141000E-02,  9.2975617155685197E-01/)
       V_rule%lambda( 78,1:3) = (/ 9.2975617155685297E-01,  5.9696109149006998E-02,  1.0547719294140029E-02/)
       V_rule%lambda( 79,1:3) = (/ 1.0547719294141000E-02,  9.2975617155685297E-01,  5.9696109149006027E-02/)

     case(  21 )
        !Wandzura's work extends the work of Dunavant by providing degree
        ! 5,10,15,20,25, and 30 rules with positive weights for the triangle.
        ! Copied on 3rd July 2008 from:
        ! http://people.scs.fsu.edu/~burkardt/f_src/wandzura/wandzura.f90
        !Stephen Wandzura, Hong Xiao, Symmetric Quadrature Rules on a Triangle,
        !Computers and Mathematics with Applications, Volume 45, Number 12, June 2003, pages 1829-1840.

        V_rule%Qdof =    85
        allocate(V_rule%weights(1:V_rule%Qdof))
        allocate(V_rule%lambda(1:V_rule%Qdof,3))

        V_rule%weights(  1) =   2.7610426997699521E-02
        V_rule%weights(  2) =   1.7790295473267400E-03
        V_rule%weights(  3) =   1.7790295473267400E-03
        V_rule%weights(  4) =   1.7790295473267400E-03
        V_rule%weights(  5) =   2.0112398113961170E-02
        V_rule%weights(  6) =   2.0112398113961170E-02
        V_rule%weights(  7) =   2.0112398113961170E-02
        V_rule%weights(  8) =   2.6817847259331569E-02
        V_rule%weights(  9) =   2.6817847259331569E-02
        V_rule%weights( 10) =   2.6817847259331569E-02
        V_rule%weights( 11) =   2.4523133801502010E-02
        V_rule%weights( 12) =   2.4523133801502010E-02
        V_rule%weights( 13) =   2.4523133801502010E-02
        V_rule%weights( 14) =   1.6394578410695391E-02
        V_rule%weights( 15) =   1.6394578410695391E-02
        V_rule%weights( 16) =   1.6394578410695391E-02
        V_rule%weights( 17) =   1.4795907398649600E-02
        V_rule%weights( 18) =   1.4795907398649600E-02
        V_rule%weights( 19) =   1.4795907398649600E-02
        V_rule%weights( 20) =   4.5792822777042507E-03
        V_rule%weights( 21) =   4.5792822777042507E-03
        V_rule%weights( 22) =   4.5792822777042507E-03
        V_rule%weights( 23) =   1.6518265155762170E-03
        V_rule%weights( 24) =   1.6518265155762170E-03
        V_rule%weights( 25) =   1.6518265155762170E-03
        V_rule%weights( 26) =   2.3491709085755839E-03
        V_rule%weights( 27) =   2.3491709085755839E-03
        V_rule%weights( 28) =   2.3491709085755839E-03
        V_rule%weights( 29) =   2.3491709085755839E-03
        V_rule%weights( 30) =   2.3491709085755839E-03
        V_rule%weights( 31) =   2.3491709085755839E-03
        V_rule%weights( 32) =   4.4659257541817933E-03
        V_rule%weights( 33) =   4.4659257541817933E-03
        V_rule%weights( 34) =   4.4659257541817933E-03
        V_rule%weights( 35) =   4.4659257541817933E-03
        V_rule%weights( 36) =   4.4659257541817933E-03
        V_rule%weights( 37) =   4.4659257541817933E-03
        V_rule%weights( 38) =   6.0995668079079721E-03
        V_rule%weights( 39) =   6.0995668079079721E-03
        V_rule%weights( 40) =   6.0995668079079721E-03
        V_rule%weights( 41) =   6.0995668079079721E-03
        V_rule%weights( 42) =   6.0995668079079721E-03
        V_rule%weights( 43) =   6.0995668079079721E-03
        V_rule%weights( 44) =   6.8910813271882030E-03
        V_rule%weights( 45) =   6.8910813271882030E-03
        V_rule%weights( 46) =   6.8910813271882030E-03
        V_rule%weights( 47) =   6.8910813271882030E-03
        V_rule%weights( 48) =   6.8910813271882030E-03
        V_rule%weights( 49) =   6.8910813271882030E-03
        V_rule%weights( 50) =   7.9974750724781628E-03
        V_rule%weights( 51) =   7.9974750724781628E-03
        V_rule%weights( 52) =   7.9974750724781628E-03
        V_rule%weights( 53) =   7.9974750724781628E-03
        V_rule%weights( 54) =   7.9974750724781628E-03
        V_rule%weights( 55) =   7.9974750724781628E-03
        V_rule%weights( 56) =   7.3861342853360238E-03
        V_rule%weights( 57) =   7.3861342853360238E-03
        V_rule%weights( 58) =   7.3861342853360238E-03
        V_rule%weights( 59) =   7.3861342853360238E-03
        V_rule%weights( 60) =   7.3861342853360238E-03
        V_rule%weights( 61) =   7.3861342853360238E-03
        V_rule%weights( 62) =   1.2799331878648261E-02
        V_rule%weights( 63) =   1.2799331878648261E-02
        V_rule%weights( 64) =   1.2799331878648261E-02
        V_rule%weights( 65) =   1.2799331878648261E-02
        V_rule%weights( 66) =   1.2799331878648261E-02
        V_rule%weights( 67) =   1.2799331878648261E-02
        V_rule%weights( 68) =   1.7258071175696551E-02
        V_rule%weights( 69) =   1.7258071175696551E-02
        V_rule%weights( 70) =   1.7258071175696551E-02
        V_rule%weights( 71) =   1.7258071175696551E-02
        V_rule%weights( 72) =   1.7258071175696551E-02
        V_rule%weights( 73) =   1.7258071175696551E-02
        V_rule%weights( 74) =   1.8672945902935469E-02
        V_rule%weights( 75) =   1.8672945902935469E-02
        V_rule%weights( 76) =   1.8672945902935469E-02
        V_rule%weights( 77) =   1.8672945902935469E-02
        V_rule%weights( 78) =   1.8672945902935469E-02
        V_rule%weights( 79) =   1.8672945902935469E-02
        V_rule%weights( 80) =   2.2818224058395259E-02
        V_rule%weights( 81) =   2.2818224058395259E-02
        V_rule%weights( 82) =   2.2818224058395259E-02
        V_rule%weights( 83) =   2.2818224058395259E-02
        V_rule%weights( 84) =   2.2818224058395259E-02
        V_rule%weights( 85) =   2.2818224058395259E-02

       V_rule%lambda(  1,1:3) = (/  3.3333333333332998E-01 ,  3.3333333333332998E-01 ,  3.3333333333334003E-01/)
       V_rule%lambda(  2,1:3) = (/  1.5006493244300000E-03 ,  4.9924967533779002E-01 ,  4.9924967533777997E-01/)
       V_rule%lambda(  3,1:3) = (/  4.9924967533779002E-01 ,  4.9924967533779002E-01 ,  1.5006493244199559E-03/)
       V_rule%lambda(  4,1:3) = (/  4.9924967533779002E-01 ,  1.5006493244300000E-03 ,  4.9924967533777997E-01/)
       V_rule%lambda(  5,1:3) = (/  9.4139751938949995E-02 ,  4.5293012403052002E-01 ,  4.5293012403053001E-01/)
       V_rule%lambda(  6,1:3) = (/  4.5293012403052002E-01 ,  4.5293012403052002E-01 ,  9.4139751938959959E-02/)
       V_rule%lambda(  7,1:3) = (/  4.5293012403052002E-01 ,  9.4139751938949995E-02 ,  4.5293012403053001E-01/)
       V_rule%lambda(  8,1:3) = (/  2.0447212408953000E-01 ,  3.9776393795524001E-01 ,  3.9776393795523002E-01/)
       V_rule%lambda(  9,1:3) = (/  3.9776393795524001E-01 ,  3.9776393795524001E-01 ,  2.0447212408951998E-01/)
       V_rule%lambda( 10,1:3) = (/  3.9776393795524001E-01 ,  2.0447212408953000E-01 ,  3.9776393795523002E-01/)
       V_rule%lambda( 11,1:3) = (/  4.7099959493443000E-01 ,  2.6450020253279000E-01 ,  2.6450020253278000E-01/)
       V_rule%lambda( 12,1:3) = (/  2.6450020253279000E-01 ,  2.6450020253279000E-01 ,  4.7099959493442001E-01/)
       V_rule%lambda( 13,1:3) = (/  2.6450020253279000E-01 ,  4.7099959493443000E-01 ,  2.6450020253278000E-01/)
       V_rule%lambda( 14,1:3) = (/  5.7796207181585002E-01 ,  2.1101896409208001E-01 ,  2.1101896409206997E-01/)
       V_rule%lambda( 15,1:3) = (/  2.1101896409208001E-01 ,  2.1101896409208001E-01 ,  5.7796207181584003E-01/)
       V_rule%lambda( 16,1:3) = (/  2.1101896409208001E-01 ,  5.7796207181585002E-01 ,  2.1101896409206997E-01/)
       V_rule%lambda( 17,1:3) = (/  7.8452878565745998E-01 ,  1.0773560717127000E-01 ,  1.0773560717127002E-01/)
       V_rule%lambda( 18,1:3) = (/  1.0773560717127000E-01 ,  1.0773560717127000E-01 ,  7.8452878565745998E-01/)
       V_rule%lambda( 19,1:3) = (/  1.0773560717127000E-01 ,  7.8452878565745998E-01 ,  1.0773560717127002E-01/)
       V_rule%lambda( 20,1:3) = (/  9.2186182432439001E-01 ,  3.9069087837799998E-02 ,  3.9069087837809990E-02/)
       V_rule%lambda( 21,1:3) = (/  3.9069087837799998E-02 ,  3.9069087837799998E-02 ,  9.2186182432440000E-01/)
       V_rule%lambda( 22,1:3) = (/  3.9069087837799998E-02 ,  9.2186182432439001E-01 ,  3.9069087837809990E-02/)
       V_rule%lambda( 23,1:3) = (/  9.7765124054134001E-01 ,  1.1174379729330001E-02 ,  1.1174379729329994E-02/)
       V_rule%lambda( 24,1:3) = (/  1.1174379729330001E-02 ,  1.1174379729330001E-02 ,  9.7765124054134001E-01/)
       V_rule%lambda( 25,1:3) = (/  1.1174379729330001E-02 ,  9.7765124054134001E-01 ,  1.1174379729329994E-02/)
       V_rule%lambda( 26,1:3) = (/  5.3496181873399998E-03 ,  6.3549665908349998E-02 ,  9.3110071590430998E-01/)
       V_rule%lambda( 27,1:3) = (/  6.3549665908349998E-02 ,  9.3110071590430998E-01 ,  5.3496181873400189E-03/)
       V_rule%lambda( 28,1:3) = (/  9.3110071590430998E-01 ,  5.3496181873399998E-03 ,  6.3549665908350012E-02/)
       V_rule%lambda( 29,1:3) = (/  6.3549665908349998E-02 ,  5.3496181873399998E-03 ,  9.3110071590430998E-01/)
       V_rule%lambda( 30,1:3) = (/  9.3110071590430998E-01 ,  6.3549665908349998E-02 ,  5.3496181873400189E-03/)
       V_rule%lambda( 31,1:3) = (/  5.3496181873399998E-03 ,  9.3110071590430998E-01 ,  6.3549665908350012E-02/)
       V_rule%lambda( 32,1:3) = (/  7.9548170661999998E-03 ,  1.5710691894070999E-01 ,  8.3493826399309001E-01/)
       V_rule%lambda( 33,1:3) = (/  1.5710691894070999E-01 ,  8.3493826399309001E-01 ,  7.9548170661999928E-03/)
       V_rule%lambda( 34,1:3) = (/  8.3493826399309001E-01 ,  7.9548170661999998E-03 ,  1.5710691894070999E-01/)
       V_rule%lambda( 35,1:3) = (/  1.5710691894070999E-01 ,  7.9548170661999998E-03 ,  8.3493826399309001E-01/)
       V_rule%lambda( 36,1:3) = (/  8.3493826399309001E-01 ,  1.5710691894070999E-01 ,  7.9548170661999928E-03/)
       V_rule%lambda( 37,1:3) = (/  7.9548170661999998E-03 ,  8.3493826399309001E-01 ,  1.5710691894070999E-01/)
       V_rule%lambda( 38,1:3) = (/  1.0422398281260001E-02 ,  3.9564211436437002E-01 ,  5.9393548735436996E-01/)
       V_rule%lambda( 39,1:3) = (/  3.9564211436437002E-01 ,  5.9393548735435997E-01 ,  1.0422398281270007E-02/)
       V_rule%lambda( 40,1:3) = (/  5.9393548735435997E-01 ,  1.0422398281260001E-02 ,  3.9564211436438002E-01/)
       V_rule%lambda( 41,1:3) = (/  3.9564211436437002E-01 ,  1.0422398281260001E-02 ,  5.9393548735436996E-01/)
       V_rule%lambda( 42,1:3) = (/  5.9393548735435997E-01 ,  3.9564211436437002E-01 ,  1.0422398281270007E-02/)
       V_rule%lambda( 43,1:3) = (/  1.0422398281260001E-02 ,  5.9393548735435997E-01 ,  3.9564211436438002E-01/)
       V_rule%lambda( 44,1:3) = (/  1.0964414796120000E-02 ,  2.7316757071290998E-01 ,  7.1586801449097004E-01/)
       V_rule%lambda( 45,1:3) = (/  2.7316757071290998E-01 ,  7.1586801449097004E-01 ,  1.0964414796119981E-02/)
       V_rule%lambda( 46,1:3) = (/  7.1586801449097004E-01 ,  1.0964414796120000E-02 ,  2.7316757071290998E-01/)
       V_rule%lambda( 47,1:3) = (/  2.7316757071290998E-01 ,  1.0964414796120000E-02 ,  7.1586801449097004E-01/)
       V_rule%lambda( 48,1:3) = (/  7.1586801449097004E-01 ,  2.7316757071290998E-01 ,  1.0964414796119981E-02/)
       V_rule%lambda( 49,1:3) = (/  1.0964414796120000E-02 ,  7.1586801449097004E-01 ,  2.7316757071290998E-01/)
       V_rule%lambda( 50,1:3) = (/  3.8566712085460003E-02 ,  1.0178538248502000E-01 ,  8.5964790542951997E-01/)
       V_rule%lambda( 51,1:3) = (/  1.0178538248502000E-01 ,  8.5964790542951997E-01 ,  3.8566712085460031E-02/)
       V_rule%lambda( 52,1:3) = (/  8.5964790542951997E-01 ,  3.8566712085460003E-02 ,  1.0178538248502003E-01/)
       V_rule%lambda( 53,1:3) = (/  1.0178538248502000E-01 ,  3.8566712085460003E-02 ,  8.5964790542951997E-01/)
       V_rule%lambda( 54,1:3) = (/  8.5964790542951997E-01 ,  1.0178538248502000E-01 ,  3.8566712085460031E-02/)
       V_rule%lambda( 55,1:3) = (/  3.8566712085460003E-02 ,  8.5964790542951997E-01 ,  1.0178538248502003E-01/)
       V_rule%lambda( 56,1:3) = (/  3.5580507817219997E-02 ,  4.4665854917641001E-01 ,  5.1776094300636999E-01/)
       V_rule%lambda( 57,1:3) = (/  4.4665854917641001E-01 ,  5.1776094300636999E-01 ,  3.5580507817219997E-02/)
       V_rule%lambda( 58,1:3) = (/  5.1776094300636999E-01 ,  3.5580507817219997E-02 ,  4.4665854917641001E-01/)
       V_rule%lambda( 59,1:3) = (/  4.4665854917641001E-01 ,  3.5580507817219997E-02 ,  5.1776094300636999E-01/)
       V_rule%lambda( 60,1:3) = (/  5.1776094300636999E-01 ,  4.4665854917641001E-01 ,  3.5580507817219997E-02/)
       V_rule%lambda( 61,1:3) = (/  3.5580507817219997E-02 ,  5.1776094300636999E-01 ,  4.4665854917641001E-01/)
       V_rule%lambda( 62,1:3) = (/  4.9670816362760002E-02 ,  1.9901079414949999E-01 ,  7.5131838948773999E-01/)
       V_rule%lambda( 63,1:3) = (/  1.9901079414949999E-01 ,  7.5131838948773000E-01 ,  4.9670816362770015E-02/)
       V_rule%lambda( 64,1:3) = (/  7.5131838948773000E-01 ,  4.9670816362760002E-02 ,  1.9901079414951001E-01/)
       V_rule%lambda( 65,1:3) = (/  1.9901079414949999E-01 ,  4.9670816362760002E-02 ,  7.5131838948773999E-01/)
       V_rule%lambda( 66,1:3) = (/  7.5131838948773000E-01 ,  1.9901079414949999E-01 ,  4.9670816362770015E-02/)
       V_rule%lambda( 67,1:3) = (/  4.9670816362760002E-02 ,  7.5131838948773000E-01 ,  1.9901079414951001E-01/)
       V_rule%lambda( 68,1:3) = (/  5.8519725084330003E-02 ,  3.2426118369228002E-01 ,  6.1721909122339003E-01/)
       V_rule%lambda( 69,1:3) = (/  3.2426118369228002E-01 ,  6.1721909122339003E-01 ,  5.8519725084329954E-02/)
       V_rule%lambda( 70,1:3) = (/  6.1721909122339003E-01 ,  5.8519725084330003E-02 ,  3.2426118369227996E-01/)
       V_rule%lambda( 71,1:3) = (/  3.2426118369228002E-01 ,  5.8519725084330003E-02 ,  6.1721909122339003E-01/)
       V_rule%lambda( 72,1:3) = (/  6.1721909122339003E-01 ,  3.2426118369228002E-01 ,  5.8519725084329954E-02/)
       V_rule%lambda( 73,1:3) = (/  5.8519725084330003E-02 ,  6.1721909122339003E-01 ,  3.2426118369227996E-01/)
       V_rule%lambda( 74,1:3) = (/  1.2149778700439000E-01 ,  2.0853136321012999E-01 ,  6.6997084978548005E-01/)
       V_rule%lambda( 75,1:3) = (/  2.0853136321012999E-01 ,  6.6997084978546995E-01 ,  1.2149778700440006E-01/)
       V_rule%lambda( 76,1:3) = (/  6.6997084978546995E-01 ,  1.2149778700439000E-01 ,  2.0853136321014004E-01/)
       V_rule%lambda( 77,1:3) = (/  2.0853136321012999E-01 ,  1.2149778700439000E-01 ,  6.6997084978548005E-01/)
       V_rule%lambda( 78,1:3) = (/  6.6997084978546995E-01 ,  2.0853136321012999E-01 ,  1.2149778700440006E-01/)
       V_rule%lambda( 79,1:3) = (/  1.2149778700439000E-01 ,  6.6997084978546995E-01 ,  2.0853136321014004E-01/)
       V_rule%lambda( 80,1:3) = (/  1.4071084494394001E-01 ,  3.2317056653625997E-01 ,  5.3611858851979999E-01/)
       V_rule%lambda( 81,1:3) = (/  3.2317056653625997E-01 ,  5.3611858851979999E-01 ,  1.4071084494394004E-01/)
       V_rule%lambda( 82,1:3) = (/  5.3611858851979999E-01 ,  1.4071084494394001E-01 ,  3.2317056653625997E-01/)
       V_rule%lambda( 83,1:3) = (/  3.2317056653625997E-01 ,  1.4071084494394001E-01 ,  5.3611858851979999E-01/)
       V_rule%lambda( 84,1:3) = (/  5.3611858851979999E-01 ,  3.2317056653625997E-01 ,  1.4071084494394004E-01/)
       V_rule%lambda( 85,1:3) = (/  1.4071084494394001E-01 ,  5.3611858851979999E-01 ,  3.2317056653625997E-01/)
     case( 22 )
        V_rule%Qdof =   126
        allocate(V_rule%weights(1:V_rule%Qdof))
        allocate(V_rule%lambda(1:V_rule%Qdof,3))

        V_rule%weights(  1) =   8.0055818800204171E-03
        V_rule%weights(  2) =   8.0055818800204171E-03
        V_rule%weights(  3) =   8.0055818800204171E-03
        V_rule%weights(  4) =   1.5947076832390501E-02
        V_rule%weights(  5) =   1.5947076832390501E-02
        V_rule%weights(  6) =   1.5947076832390501E-02
        V_rule%weights(  7) =   1.3109141230795530E-02
        V_rule%weights(  8) =   1.3109141230795530E-02
        V_rule%weights(  9) =   1.3109141230795530E-02
        V_rule%weights( 10) =   1.9583000965635620E-02
        V_rule%weights( 11) =   1.9583000965635620E-02
        V_rule%weights( 12) =   1.9583000965635620E-02
        V_rule%weights( 13) =   1.6470885441537270E-02
        V_rule%weights( 14) =   1.6470885441537270E-02
        V_rule%weights( 15) =   1.6470885441537270E-02
        V_rule%weights( 16) =   8.5472790740921002E-03
        V_rule%weights( 17) =   8.5472790740921002E-03
        V_rule%weights( 18) =   8.5472790740921002E-03
        V_rule%weights( 19) =   8.1618858572264918E-03
        V_rule%weights( 20) =   8.1618858572264918E-03
        V_rule%weights( 21) =   8.1618858572264918E-03
        V_rule%weights( 22) =   6.1211465399837791E-03
        V_rule%weights( 23) =   6.1211465399837791E-03
        V_rule%weights( 24) =   6.1211465399837791E-03
        V_rule%weights( 25) =   2.9084982649366649E-03
        V_rule%weights( 26) =   2.9084982649366649E-03
        V_rule%weights( 27) =   2.9084982649366649E-03
        V_rule%weights( 28) =   6.9227524566199629E-04
        V_rule%weights( 29) =   6.9227524566199629E-04
        V_rule%weights( 30) =   6.9227524566199629E-04
        V_rule%weights( 31) =   1.2482891992773970E-03
        V_rule%weights( 32) =   1.2482891992773970E-03
        V_rule%weights( 33) =   1.2482891992773970E-03
        V_rule%weights( 34) =   1.2482891992773970E-03
        V_rule%weights( 35) =   1.2482891992773970E-03
        V_rule%weights( 36) =   1.2482891992773970E-03
        V_rule%weights( 37) =   3.4047529088030220E-03
        V_rule%weights( 38) =   3.4047529088030220E-03
        V_rule%weights( 39) =   3.4047529088030220E-03
        V_rule%weights( 40) =   3.4047529088030220E-03
        V_rule%weights( 41) =   3.4047529088030220E-03
        V_rule%weights( 42) =   3.4047529088030220E-03
        V_rule%weights( 43) =   3.3596543260640509E-03
        V_rule%weights( 44) =   3.3596543260640509E-03
        V_rule%weights( 45) =   3.3596543260640509E-03
        V_rule%weights( 46) =   3.3596543260640509E-03
        V_rule%weights( 47) =   3.3596543260640509E-03
        V_rule%weights( 48) =   3.3596543260640509E-03
        V_rule%weights( 49) =   1.7161565394967541E-03
        V_rule%weights( 50) =   1.7161565394967541E-03
        V_rule%weights( 51) =   1.7161565394967541E-03
        V_rule%weights( 52) =   1.7161565394967541E-03
        V_rule%weights( 53) =   1.7161565394967541E-03
        V_rule%weights( 54) =   1.7161565394967541E-03
        V_rule%weights( 55) =   1.4808563167156060E-03
        V_rule%weights( 56) =   1.4808563167156060E-03
        V_rule%weights( 57) =   1.4808563167156060E-03
        V_rule%weights( 58) =   1.4808563167156060E-03
        V_rule%weights( 59) =   1.4808563167156060E-03
        V_rule%weights( 60) =   1.4808563167156060E-03
        V_rule%weights( 61) =   3.5113126107286850E-03
        V_rule%weights( 62) =   3.5113126107286850E-03
        V_rule%weights( 63) =   3.5113126107286850E-03
        V_rule%weights( 64) =   3.5113126107286850E-03
        V_rule%weights( 65) =   3.5113126107286850E-03
        V_rule%weights( 66) =   3.5113126107286850E-03
        V_rule%weights( 67) =   7.3935501497064838E-03
        V_rule%weights( 68) =   7.3935501497064838E-03
        V_rule%weights( 69) =   7.3935501497064838E-03
        V_rule%weights( 70) =   7.3935501497064838E-03
        V_rule%weights( 71) =   7.3935501497064838E-03
        V_rule%weights( 72) =   7.3935501497064838E-03
        V_rule%weights( 73) =   7.9830874773765582E-03
        V_rule%weights( 74) =   7.9830874773765582E-03
        V_rule%weights( 75) =   7.9830874773765582E-03
        V_rule%weights( 76) =   7.9830874773765582E-03
        V_rule%weights( 77) =   7.9830874773765582E-03
        V_rule%weights( 78) =   7.9830874773765582E-03
        V_rule%weights( 79) =   4.3559626131580414E-03
        V_rule%weights( 80) =   4.3559626131580414E-03
        V_rule%weights( 81) =   4.3559626131580414E-03
        V_rule%weights( 82) =   4.3559626131580414E-03
        V_rule%weights( 83) =   4.3559626131580414E-03
        V_rule%weights( 84) =   4.3559626131580414E-03
        V_rule%weights( 85) =   7.3650567014178318E-03
        V_rule%weights( 86) =   7.3650567014178318E-03
        V_rule%weights( 87) =   7.3650567014178318E-03
        V_rule%weights( 88) =   7.3650567014178318E-03
        V_rule%weights( 89) =   7.3650567014178318E-03
        V_rule%weights( 90) =   7.3650567014178318E-03
        V_rule%weights( 91) =   1.0963572846419550E-02
        V_rule%weights( 92) =   1.0963572846419550E-02
        V_rule%weights( 93) =   1.0963572846419550E-02
        V_rule%weights( 94) =   1.0963572846419550E-02
        V_rule%weights( 95) =   1.0963572846419550E-02
        V_rule%weights( 96) =   1.0963572846419550E-02
        V_rule%weights( 97) =   1.1749961743541121E-02
        V_rule%weights( 98) =   1.1749961743541121E-02
        V_rule%weights( 99) =   1.1749961743541121E-02
        V_rule%weights(100) =   1.1749961743541121E-02
        V_rule%weights(101) =   1.1749961743541121E-02
        V_rule%weights(102) =   1.1749961743541121E-02
        V_rule%weights(103) =   1.0015600713798570E-02
        V_rule%weights(104) =   1.0015600713798570E-02
        V_rule%weights(105) =   1.0015600713798570E-02
        V_rule%weights(106) =   1.0015600713798570E-02
        V_rule%weights(107) =   1.0015600713798570E-02
        V_rule%weights(108) =   1.0015600713798570E-02
        V_rule%weights(109) =   1.3309640787628680E-02
        V_rule%weights(110) =   1.3309640787628680E-02
        V_rule%weights(111) =   1.3309640787628680E-02
        V_rule%weights(112) =   1.3309640787628680E-02
        V_rule%weights(113) =   1.3309640787628680E-02
        V_rule%weights(114) =   1.3309640787628680E-02
        V_rule%weights(115) =   1.4154446505226140E-02
        V_rule%weights(116) =   1.4154446505226140E-02
        V_rule%weights(117) =   1.4154446505226140E-02
        V_rule%weights(118) =   1.4154446505226140E-02
        V_rule%weights(119) =   1.4154446505226140E-02
        V_rule%weights(120) =   1.4154446505226140E-02
        V_rule%weights(121) =   1.4881379561168010E-02
        V_rule%weights(122) =   1.4881379561168010E-02
        V_rule%weights(123) =   1.4881379561168010E-02
        V_rule%weights(124) =   1.4881379561168010E-02
        V_rule%weights(125) =   1.4881379561168010E-02
        V_rule%weights(126) =   1.4881379561168010E-02

       V_rule%lambda(  1,1:3) = (/  2.7946483073169999E-02 ,  4.8602675846340998E-01 ,  4.8602675846342003E-01/)
       V_rule%lambda(  2,1:3) = (/  4.8602675846340998E-01 ,  4.8602675846340998E-01 ,  2.7946483073180040E-02/)
       V_rule%lambda(  3,1:3) = (/  4.8602675846340998E-01 ,  2.7946483073169999E-02 ,  4.8602675846342003E-01/)
       V_rule%lambda(  4,1:3) = (/  1.3117860132765000E-01 ,  4.3441069933616999E-01 ,  4.3441069933618004E-01/)
       V_rule%lambda(  5,1:3) = (/  4.3441069933616999E-01 ,  4.3441069933616999E-01 ,  1.3117860132766002E-01/)
       V_rule%lambda(  6,1:3) = (/  4.3441069933616999E-01 ,  1.3117860132765000E-01 ,  4.3441069933618004E-01/)
       V_rule%lambda(  7,1:3) = (/  2.2022172951207000E-01 ,  3.8988913524396002E-01 ,  3.8988913524396995E-01/)
       V_rule%lambda(  8,1:3) = (/  3.8988913524396002E-01 ,  3.8988913524396002E-01 ,  2.2022172951207997E-01/)
       V_rule%lambda(  9,1:3) = (/  3.8988913524396002E-01 ,  2.2022172951207000E-01 ,  3.8988913524396995E-01/)
       V_rule%lambda( 10,1:3) = (/  4.0311353196039001E-01 ,  2.9844323401980000E-01 ,  2.9844323401980999E-01/)
       V_rule%lambda( 11,1:3) = (/  2.9844323401980000E-01 ,  2.9844323401980000E-01 ,  4.0311353196040001E-01/)
       V_rule%lambda( 12,1:3) = (/  2.9844323401980000E-01 ,  4.0311353196039001E-01 ,  2.9844323401980999E-01/)
       V_rule%lambda( 13,1:3) = (/  5.3191165532525997E-01 ,  2.3404417233736999E-01 ,  2.3404417233737004E-01/)
       V_rule%lambda( 14,1:3) = (/  2.3404417233736999E-01 ,  2.3404417233736999E-01 ,  5.3191165532526008E-01/)
       V_rule%lambda( 15,1:3) = (/  2.3404417233736999E-01 ,  5.3191165532525997E-01 ,  2.3404417233737004E-01/)
       V_rule%lambda( 16,1:3) = (/  6.9706333078196003E-01 ,  1.5146833460902001E-01 ,  1.5146833460901996E-01/)
       V_rule%lambda( 17,1:3) = (/  1.5146833460902001E-01 ,  1.5146833460902001E-01 ,  6.9706333078196003E-01/)
       V_rule%lambda( 18,1:3) = (/  1.5146833460902001E-01 ,  6.9706333078196003E-01 ,  1.5146833460901996E-01/)
       V_rule%lambda( 19,1:3) = (/  7.7453221290801000E-01 ,  1.1273389354599000E-01 ,  1.1273389354600000E-01/)
       V_rule%lambda( 20,1:3) = (/  1.1273389354599000E-01 ,  1.1273389354599000E-01 ,  7.7453221290801999E-01/)
       V_rule%lambda( 21,1:3) = (/  1.1273389354599000E-01 ,  7.7453221290801000E-01 ,  1.1273389354600000E-01/)
       V_rule%lambda( 22,1:3) = (/  8.4456861581694997E-01 ,  7.7715692091529995E-02 ,  7.7715692091520031E-02/)
       V_rule%lambda( 23,1:3) = (/  7.7715692091529995E-02 ,  7.7715692091529995E-02 ,  8.4456861581693998E-01/)
       V_rule%lambda( 24,1:3) = (/  7.7715692091529995E-02 ,  8.4456861581694997E-01 ,  7.7715692091520031E-02/)
       V_rule%lambda( 25,1:3) = (/  9.3021381277141002E-01 ,  3.4893093614300000E-02 ,  3.4893093614289980E-02/)
       V_rule%lambda( 26,1:3) = (/  3.4893093614300000E-02 ,  3.4893093614300000E-02 ,  9.3021381277140003E-01/)
       V_rule%lambda( 27,1:3) = (/  3.4893093614300000E-02 ,  9.3021381277141002E-01 ,  3.4893093614289980E-02/)
       V_rule%lambda( 28,1:3) = (/  9.8548363075812995E-01 ,  7.2581846209300001E-03 ,  7.2581846209400528E-03/)
       V_rule%lambda( 29,1:3) = (/  7.2581846209300001E-03 ,  7.2581846209300001E-03 ,  9.8548363075814005E-01/)
       V_rule%lambda( 30,1:3) = (/  7.2581846209300001E-03 ,  9.8548363075812995E-01 ,  7.2581846209400528E-03/)
       V_rule%lambda( 31,1:3) = (/  1.2923527044399999E-03 ,  2.2721445215336000E-01 ,  7.7149319514220005E-01/)
       V_rule%lambda( 32,1:3) = (/  2.2721445215336000E-01 ,  7.7149319514218995E-01 ,  1.2923527044500505E-03/)
       V_rule%lambda( 33,1:3) = (/  7.7149319514218995E-01 ,  1.2923527044399999E-03 ,  2.2721445215337005E-01/)
       V_rule%lambda( 34,1:3) = (/  2.2721445215336000E-01 ,  1.2923527044399999E-03 ,  7.7149319514220005E-01/)
       V_rule%lambda( 35,1:3) = (/  7.7149319514218995E-01 ,  2.2721445215336000E-01 ,  1.2923527044500505E-03/)
       V_rule%lambda( 36,1:3) = (/  1.2923527044399999E-03 ,  7.7149319514218995E-01 ,  2.2721445215337005E-01/)
       V_rule%lambda( 37,1:3) = (/  5.3997012721200000E-03 ,  4.3501055485356999E-01 ,  5.5958974387431004E-01/)
       V_rule%lambda( 38,1:3) = (/  4.3501055485356999E-01 ,  5.5958974387431004E-01 ,  5.3997012721199722E-03/)
       V_rule%lambda( 39,1:3) = (/  5.5958974387431004E-01 ,  5.3997012721200000E-03 ,  4.3501055485356999E-01/)
       V_rule%lambda( 40,1:3) = (/  4.3501055485356999E-01 ,  5.3997012721200000E-03 ,  5.5958974387431004E-01/)
       V_rule%lambda( 41,1:3) = (/  5.5958974387431004E-01 ,  4.3501055485356999E-01 ,  5.3997012721199722E-03/)
       V_rule%lambda( 42,1:3) = (/  5.3997012721200000E-03 ,  5.5958974387431004E-01 ,  4.3501055485356999E-01/)
       V_rule%lambda( 43,1:3) = (/  6.3840030339800003E-03 ,  3.2030959927219999E-01 ,  6.7330639769382006E-01/)
       V_rule%lambda( 44,1:3) = (/  3.2030959927219999E-01 ,  6.7330639769381995E-01 ,  6.3840030339800680E-03/)
       V_rule%lambda( 45,1:3) = (/  6.7330639769381995E-01 ,  6.3840030339800003E-03 ,  3.2030959927220004E-01/)
       V_rule%lambda( 46,1:3) = (/  3.2030959927219999E-01 ,  6.3840030339800003E-03 ,  6.7330639769382006E-01/)
       V_rule%lambda( 47,1:3) = (/  6.7330639769381995E-01 ,  3.2030959927219999E-01 ,  6.3840030339800680E-03/)
       V_rule%lambda( 48,1:3) = (/  6.3840030339800003E-03 ,  6.7330639769381995E-01 ,  3.2030959927220004E-01/)
       V_rule%lambda( 49,1:3) = (/  5.0282115019900002E-03 ,  9.1750322280009997E-02 ,  9.0322146621800004E-01/)
       V_rule%lambda( 50,1:3) = (/  9.1750322280009997E-02 ,  9.0322146621800004E-01 ,  5.0282115019899681E-03/)
       V_rule%lambda( 51,1:3) = (/  9.0322146621800004E-01 ,  5.0282115019900002E-03 ,  9.1750322280009969E-02/)
       V_rule%lambda( 52,1:3) = (/  9.1750322280009997E-02 ,  5.0282115019900002E-03 ,  9.0322146621800004E-01/)
       V_rule%lambda( 53,1:3) = (/  9.0322146621800004E-01 ,  9.1750322280009997E-02 ,  5.0282115019899681E-03/)
       V_rule%lambda( 54,1:3) = (/  5.0282115019900002E-03 ,  9.0322146621800004E-01 ,  9.1750322280009969E-02/)
       V_rule%lambda( 55,1:3) = (/  6.8267586217800004E-03 ,  3.8010835858719998E-02 ,  9.5516240551950005E-01/)
       V_rule%lambda( 56,1:3) = (/  3.8010835858719998E-02 ,  9.5516240551949005E-01 ,  6.8267586217899481E-03/)
       V_rule%lambda( 57,1:3) = (/  9.5516240551949005E-01 ,  6.8267586217800004E-03 ,  3.8010835858729948E-02/)
       V_rule%lambda( 58,1:3) = (/  3.8010835858719998E-02 ,  6.8267586217800004E-03 ,  9.5516240551950005E-01/)
       V_rule%lambda( 59,1:3) = (/  9.5516240551949005E-01 ,  3.8010835858719998E-02 ,  6.8267586217899481E-03/)
       V_rule%lambda( 60,1:3) = (/  6.8267586217800004E-03 ,  9.5516240551949005E-01 ,  3.8010835858729948E-02/)
       V_rule%lambda( 61,1:3) = (/  1.0016199639930000E-02 ,  1.5742521848530999E-01 ,  8.3255858187476006E-01/)
       V_rule%lambda( 62,1:3) = (/  1.5742521848530999E-01 ,  8.3255858187475995E-01 ,  1.0016199639930057E-02/)
       V_rule%lambda( 63,1:3) = (/  8.3255858187475995E-01 ,  1.0016199639930000E-02 ,  1.5742521848531005E-01/)
       V_rule%lambda( 64,1:3) = (/  1.5742521848530999E-01 ,  1.0016199639930000E-02 ,  8.3255858187476006E-01/)
       V_rule%lambda( 65,1:3) = (/  8.3255858187475995E-01 ,  1.5742521848530999E-01 ,  1.0016199639930057E-02/)
       V_rule%lambda( 66,1:3) = (/  1.0016199639930000E-02 ,  8.3255858187475995E-01 ,  1.5742521848531005E-01/)
       V_rule%lambda( 67,1:3) = (/  2.5757813173390001E-02 ,  2.3988965977853000E-01 ,  7.3435252704807996E-01/)
       V_rule%lambda( 68,1:3) = (/  2.3988965977853000E-01 ,  7.3435252704807996E-01 ,  2.5757813173390043E-02/)
       V_rule%lambda( 69,1:3) = (/  7.3435252704807996E-01 ,  2.5757813173390001E-02 ,  2.3988965977853005E-01/)
       V_rule%lambda( 70,1:3) = (/  2.3988965977853000E-01 ,  2.5757813173390001E-02 ,  7.3435252704807996E-01/)
       V_rule%lambda( 71,1:3) = (/  7.3435252704807996E-01 ,  2.3988965977853000E-01 ,  2.5757813173390043E-02/)
       V_rule%lambda( 72,1:3) = (/  2.5757813173390001E-02 ,  7.3435252704807996E-01 ,  2.3988965977853005E-01/)
       V_rule%lambda( 73,1:3) = (/  3.0227898119920001E-02 ,  3.6194311812606000E-01 ,  6.0782898375401995E-01/)
       V_rule%lambda( 74,1:3) = (/  3.6194311812606000E-01 ,  6.0782898375401995E-01 ,  3.0227898119920049E-02/)
       V_rule%lambda( 75,1:3) = (/  6.0782898375401995E-01 ,  3.0227898119920001E-02 ,  3.6194311812606006E-01/)
       V_rule%lambda( 76,1:3) = (/  3.6194311812606000E-01 ,  3.0227898119920001E-02 ,  6.0782898375401995E-01/)
       V_rule%lambda( 77,1:3) = (/  6.0782898375401995E-01 ,  3.6194311812606000E-01 ,  3.0227898119920049E-02/)
       V_rule%lambda( 78,1:3) = (/  3.0227898119920001E-02 ,  6.0782898375401995E-01 ,  3.6194311812606006E-01/)
       V_rule%lambda( 79,1:3) = (/  3.0504990107159999E-02 ,  8.3551960954830001E-02 ,  8.8594304893801001E-01/)
       V_rule%lambda( 80,1:3) = (/  8.3551960954830001E-02 ,  8.8594304893801001E-01 ,  3.0504990107159985E-02/)
       V_rule%lambda( 81,1:3) = (/  8.8594304893801001E-01 ,  3.0504990107159999E-02 ,  8.3551960954829987E-02/)
       V_rule%lambda( 82,1:3) = (/  8.3551960954830001E-02 ,  3.0504990107159999E-02 ,  8.8594304893801001E-01/)
       V_rule%lambda( 83,1:3) = (/  8.8594304893801001E-01 ,  8.3551960954830001E-02 ,  3.0504990107159985E-02/)
       V_rule%lambda( 84,1:3) = (/  3.0504990107159999E-02 ,  8.8594304893801001E-01 ,  8.3551960954829987E-02/)
       V_rule%lambda( 85,1:3) = (/  4.5956547362570002E-02 ,  1.4844322073242000E-01 ,  8.0560023190500996E-01/)
       V_rule%lambda( 86,1:3) = (/  1.4844322073242000E-01 ,  8.0560023190500996E-01 ,  4.5956547362570044E-02/)
       V_rule%lambda( 87,1:3) = (/  8.0560023190500996E-01 ,  4.5956547362570002E-02 ,  1.4844322073242006E-01/)
       V_rule%lambda( 88,1:3) = (/  1.4844322073242000E-01 ,  4.5956547362570002E-02 ,  8.0560023190500996E-01/)
       V_rule%lambda( 89,1:3) = (/  8.0560023190500996E-01 ,  1.4844322073242000E-01 ,  4.5956547362570044E-02/)
       V_rule%lambda( 90,1:3) = (/  4.5956547362570002E-02 ,  8.0560023190500996E-01 ,  1.4844322073242006E-01/)
       V_rule%lambda( 91,1:3) = (/  6.7442800540279998E-02 ,  2.8373970872753002E-01 ,  6.4881749073218997E-01/)
       V_rule%lambda( 92,1:3) = (/  2.8373970872753002E-01 ,  6.4881749073218997E-01 ,  6.7442800540280012E-02/)
       V_rule%lambda( 93,1:3) = (/  6.4881749073218997E-01 ,  6.7442800540279998E-02 ,  2.8373970872753002E-01/)
       V_rule%lambda( 94,1:3) = (/  2.8373970872753002E-01 ,  6.7442800540279998E-02 ,  6.4881749073218997E-01/)
       V_rule%lambda( 95,1:3) = (/  6.4881749073218997E-01 ,  2.8373970872753002E-01 ,  6.7442800540280012E-02/)
       V_rule%lambda( 96,1:3) = (/  6.7442800540279998E-02 ,  6.4881749073218997E-01 ,  2.8373970872753002E-01/)
       V_rule%lambda( 97,1:3) = (/  7.0045091415910005E-02 ,  4.0689937511878999E-01 ,  5.2305553346529998E-01/)
       V_rule%lambda( 98,1:3) = (/  4.0689937511878999E-01 ,  5.2305553346529998E-01 ,  7.0045091415910032E-02/)
       V_rule%lambda( 99,1:3) = (/  5.2305553346529998E-01 ,  7.0045091415910005E-02 ,  4.0689937511879004E-01/)
       V_rule%lambda(100,1:3) = (/  4.0689937511878999E-01 ,  7.0045091415910005E-02 ,  5.2305553346529998E-01/)
       V_rule%lambda(101,1:3) = (/  5.2305553346529998E-01 ,  4.0689937511878999E-01 ,  7.0045091415910032E-02/)
       V_rule%lambda(102,1:3) = (/  7.0045091415910005E-02 ,  5.2305553346529998E-01 ,  4.0689937511879004E-01/)
       V_rule%lambda(103,1:3) = (/  8.3911524640120000E-02 ,  1.9411398702488999E-01 ,  7.2197448833499001E-01/)
       V_rule%lambda(104,1:3) = (/  1.9411398702488999E-01 ,  7.2197448833499001E-01 ,  8.3911524640120000E-02/)
       V_rule%lambda(105,1:3) = (/  7.2197448833499001E-01 ,  8.3911524640120000E-02 ,  1.9411398702488999E-01/)
       V_rule%lambda(106,1:3) = (/  1.9411398702488999E-01 ,  8.3911524640120000E-02 ,  7.2197448833499001E-01/)
       V_rule%lambda(107,1:3) = (/  7.2197448833499001E-01 ,  1.9411398702488999E-01 ,  8.3911524640120000E-02/)
       V_rule%lambda(108,1:3) = (/  8.3911524640120000E-02 ,  7.2197448833499001E-01 ,  1.9411398702488999E-01/)
       V_rule%lambda(109,1:3) = (/  1.2037553567714999E-01 ,  3.2413434700069998E-01 ,  5.5549011732215003E-01/)
       V_rule%lambda(110,1:3) = (/  3.2413434700069998E-01 ,  5.5549011732214004E-01 ,  1.2037553567715997E-01/)
       V_rule%lambda(111,1:3) = (/  5.5549011732214004E-01 ,  1.2037553567714999E-01 ,  3.2413434700070998E-01/)
       V_rule%lambda(112,1:3) = (/  3.2413434700069998E-01 ,  1.2037553567714999E-01 ,  5.5549011732215003E-01/)
       V_rule%lambda(113,1:3) = (/  5.5549011732214004E-01 ,  3.2413434700069998E-01 ,  1.2037553567715997E-01/)
       V_rule%lambda(114,1:3) = (/  1.2037553567714999E-01 ,  5.5549011732214004E-01 ,  3.2413434700070998E-01/)
       V_rule%lambda(115,1:3) = (/  1.4806689915737001E-01 ,  2.2927748355597999E-01 ,  6.2265561728664998E-01/)
       V_rule%lambda(116,1:3) = (/  2.2927748355597999E-01 ,  6.2265561728664998E-01 ,  1.4806689915737004E-01/)
       V_rule%lambda(117,1:3) = (/  6.2265561728664998E-01 ,  1.4806689915737001E-01 ,  2.2927748355598002E-01/)
       V_rule%lambda(118,1:3) = (/  2.2927748355597999E-01 ,  1.4806689915737001E-01 ,  6.2265561728664998E-01/)
       V_rule%lambda(119,1:3) = (/  6.2265561728664998E-01 ,  2.2927748355597999E-01 ,  1.4806689915737004E-01/)
       V_rule%lambda(120,1:3) = (/  1.4806689915737001E-01 ,  6.2265561728664998E-01 ,  2.2927748355598002E-01/)
       V_rule%lambda(121,1:3) = (/  1.9177186586733000E-01 ,  3.2561812259598000E-01 ,  4.8261001153669003E-01/)
       V_rule%lambda(122,1:3) = (/  3.2561812259598000E-01 ,  4.8261001153668998E-01 ,  1.9177186586733003E-01/)
       V_rule%lambda(123,1:3) = (/  4.8261001153668998E-01 ,  1.9177186586733000E-01 ,  3.2561812259598000E-01/)
       V_rule%lambda(124,1:3) = (/  3.2561812259598000E-01 ,  1.9177186586733000E-01 ,  4.8261001153669003E-01/)
       V_rule%lambda(125,1:3) = (/  4.8261001153668998E-01 ,  3.2561812259598000E-01 ,  1.9177186586733003E-01/)
       V_rule%lambda(126,1:3) = (/  1.9177186586733000E-01 ,  4.8261001153668998E-01 ,  3.2561812259598000E-01/)
     case( 23 )
        V_rule%Qdof =   175
        allocate(V_rule%weights(1:V_rule%Qdof))
        allocate(V_rule%lambda(1:V_rule%Qdof,3))

        V_rule%weights(  1) =   1.5579960202899200E-02
        V_rule%weights(  2) =   3.1772337005341340E-03
        V_rule%weights(  3) =   3.1772337005341340E-03
        V_rule%weights(  4) =   3.1772337005341340E-03
        V_rule%weights(  5) =   1.0483426635730771E-02
        V_rule%weights(  6) =   1.0483426635730771E-02
        V_rule%weights(  7) =   1.0483426635730771E-02
        V_rule%weights(  8) =   1.3209459577743630E-02
        V_rule%weights(  9) =   1.3209459577743630E-02
        V_rule%weights( 10) =   1.3209459577743630E-02
        V_rule%weights( 11) =   1.4975006966271499E-02
        V_rule%weights( 12) =   1.4975006966271499E-02
        V_rule%weights( 13) =   1.4975006966271499E-02
        V_rule%weights( 14) =   1.4987904443384190E-02
        V_rule%weights( 15) =   1.4987904443384190E-02
        V_rule%weights( 16) =   1.4987904443384190E-02
        V_rule%weights( 17) =   1.3338864741021660E-02
        V_rule%weights( 18) =   1.3338864741021660E-02
        V_rule%weights( 19) =   1.3338864741021660E-02
        V_rule%weights( 20) =   1.0889171113902011E-02
        V_rule%weights( 21) =   1.0889171113902011E-02
        V_rule%weights( 22) =   1.0889171113902011E-02
        V_rule%weights( 23) =   8.1894406608934607E-03
        V_rule%weights( 24) =   8.1894406608934607E-03
        V_rule%weights( 25) =   8.1894406608934607E-03
        V_rule%weights( 26) =   5.5753875886077851E-03
        V_rule%weights( 27) =   5.5753875886077851E-03
        V_rule%weights( 28) =   5.5753875886077851E-03
        V_rule%weights( 29) =   3.1912164734119760E-03
        V_rule%weights( 30) =   3.1912164734119760E-03
        V_rule%weights( 31) =   3.1912164734119760E-03
        V_rule%weights( 32) =   1.2967151443270450E-03
        V_rule%weights( 33) =   1.2967151443270450E-03
        V_rule%weights( 34) =   1.2967151443270450E-03
        V_rule%weights( 35) =   2.9826282613491719E-04
        V_rule%weights( 36) =   2.9826282613491719E-04
        V_rule%weights( 37) =   2.9826282613491719E-04
        V_rule%weights( 38) =   9.9890568507889641E-04
        V_rule%weights( 39) =   9.9890568507889641E-04
        V_rule%weights( 40) =   9.9890568507889641E-04
        V_rule%weights( 41) =   9.9890568507889641E-04
        V_rule%weights( 42) =   9.9890568507889641E-04
        V_rule%weights( 43) =   9.9890568507889641E-04
        V_rule%weights( 44) =   4.6285084917325331E-04
        V_rule%weights( 45) =   4.6285084917325331E-04
        V_rule%weights( 46) =   4.6285084917325331E-04
        V_rule%weights( 47) =   4.6285084917325331E-04
        V_rule%weights( 48) =   4.6285084917325331E-04
        V_rule%weights( 49) =   4.6285084917325331E-04
        V_rule%weights( 50) =   1.2344513363824129E-03
        V_rule%weights( 51) =   1.2344513363824129E-03
        V_rule%weights( 52) =   1.2344513363824129E-03
        V_rule%weights( 53) =   1.2344513363824129E-03
        V_rule%weights( 54) =   1.2344513363824129E-03
        V_rule%weights( 55) =   1.2344513363824129E-03
        V_rule%weights( 56) =   5.7071985224320615E-04
        V_rule%weights( 57) =   5.7071985224320615E-04
        V_rule%weights( 58) =   5.7071985224320615E-04
        V_rule%weights( 59) =   5.7071985224320615E-04
        V_rule%weights( 60) =   5.7071985224320615E-04
        V_rule%weights( 61) =   5.7071985224320615E-04
        V_rule%weights( 62) =   1.1269461258776241E-03
        V_rule%weights( 63) =   1.1269461258776241E-03
        V_rule%weights( 64) =   1.1269461258776241E-03
        V_rule%weights( 65) =   1.1269461258776241E-03
        V_rule%weights( 66) =   1.1269461258776241E-03
        V_rule%weights( 67) =   1.1269461258776241E-03
        V_rule%weights( 68) =   1.7478669494073371E-03
        V_rule%weights( 69) =   1.7478669494073371E-03
        V_rule%weights( 70) =   1.7478669494073371E-03
        V_rule%weights( 71) =   1.7478669494073371E-03
        V_rule%weights( 72) =   1.7478669494073371E-03
        V_rule%weights( 73) =   1.7478669494073371E-03
        V_rule%weights( 74) =   1.1828188150316569E-03
        V_rule%weights( 75) =   1.1828188150316569E-03
        V_rule%weights( 76) =   1.1828188150316569E-03
        V_rule%weights( 77) =   1.1828188150316569E-03
        V_rule%weights( 78) =   1.1828188150316569E-03
        V_rule%weights( 79) =   1.1828188150316569E-03
        V_rule%weights( 80) =   1.9908392946750338E-03
        V_rule%weights( 81) =   1.9908392946750338E-03
        V_rule%weights( 82) =   1.9908392946750338E-03
        V_rule%weights( 83) =   1.9908392946750338E-03
        V_rule%weights( 84) =   1.9908392946750338E-03
        V_rule%weights( 85) =   1.9908392946750338E-03
        V_rule%weights( 86) =   1.9004127950359800E-03
        V_rule%weights( 87) =   1.9004127950359800E-03
        V_rule%weights( 88) =   1.9004127950359800E-03
        V_rule%weights( 89) =   1.9004127950359800E-03
        V_rule%weights( 90) =   1.9004127950359800E-03
        V_rule%weights( 91) =   1.9004127950359800E-03
        V_rule%weights( 92) =   4.4983658088174512E-03
        V_rule%weights( 93) =   4.4983658088174512E-03
        V_rule%weights( 94) =   4.4983658088174512E-03
        V_rule%weights( 95) =   4.4983658088174512E-03
        V_rule%weights( 96) =   4.4983658088174512E-03
        V_rule%weights( 97) =   4.4983658088174512E-03
        V_rule%weights( 98) =   3.4787194602747189E-03
        V_rule%weights( 99) =   3.4787194602747189E-03
        V_rule%weights(100) =   3.4787194602747189E-03
        V_rule%weights(101) =   3.4787194602747189E-03
        V_rule%weights(102) =   3.4787194602747189E-03
        V_rule%weights(103) =   3.4787194602747189E-03
        V_rule%weights(104) =   4.1023990367239534E-03
        V_rule%weights(105) =   4.1023990367239534E-03
        V_rule%weights(106) =   4.1023990367239534E-03
        V_rule%weights(107) =   4.1023990367239534E-03
        V_rule%weights(108) =   4.1023990367239534E-03
        V_rule%weights(109) =   4.1023990367239534E-03
        V_rule%weights(110) =   4.0217615497441621E-03
        V_rule%weights(111) =   4.0217615497441621E-03
        V_rule%weights(112) =   4.0217615497441621E-03
        V_rule%weights(113) =   4.0217615497441621E-03
        V_rule%weights(114) =   4.0217615497441621E-03
        V_rule%weights(115) =   4.0217615497441621E-03
        V_rule%weights(116) =   6.0331646607950659E-03
        V_rule%weights(117) =   6.0331646607950659E-03
        V_rule%weights(118) =   6.0331646607950659E-03
        V_rule%weights(119) =   6.0331646607950659E-03
        V_rule%weights(120) =   6.0331646607950659E-03
        V_rule%weights(121) =   6.0331646607950659E-03
        V_rule%weights(122) =   3.9462903021295981E-03
        V_rule%weights(123) =   3.9462903021295981E-03
        V_rule%weights(124) =   3.9462903021295981E-03
        V_rule%weights(125) =   3.9462903021295981E-03
        V_rule%weights(126) =   3.9462903021295981E-03
        V_rule%weights(127) =   3.9462903021295981E-03
        V_rule%weights(128) =   6.6440445376802684E-03
        V_rule%weights(129) =   6.6440445376802684E-03
        V_rule%weights(130) =   6.6440445376802684E-03
        V_rule%weights(131) =   6.6440445376802684E-03
        V_rule%weights(132) =   6.6440445376802684E-03
        V_rule%weights(133) =   6.6440445376802684E-03
        V_rule%weights(134) =   8.2543058560784581E-03
        V_rule%weights(135) =   8.2543058560784581E-03
        V_rule%weights(136) =   8.2543058560784581E-03
        V_rule%weights(137) =   8.2543058560784581E-03
        V_rule%weights(138) =   8.2543058560784581E-03
        V_rule%weights(139) =   8.2543058560784581E-03
        V_rule%weights(140) =   6.4960566334064107E-03
        V_rule%weights(141) =   6.4960566334064107E-03
        V_rule%weights(142) =   6.4960566334064107E-03
        V_rule%weights(143) =   6.4960566334064107E-03
        V_rule%weights(144) =   6.4960566334064107E-03
        V_rule%weights(145) =   6.4960566334064107E-03
        V_rule%weights(146) =   9.2527781441466023E-03
        V_rule%weights(147) =   9.2527781441466023E-03
        V_rule%weights(148) =   9.2527781441466023E-03
        V_rule%weights(149) =   9.2527781441466023E-03
        V_rule%weights(150) =   9.2527781441466023E-03
        V_rule%weights(151) =   9.2527781441466023E-03
        V_rule%weights(152) =   9.1649207262942799E-03
        V_rule%weights(153) =   9.1649207262942799E-03
        V_rule%weights(154) =   9.1649207262942799E-03
        V_rule%weights(155) =   9.1649207262942799E-03
        V_rule%weights(156) =   9.1649207262942799E-03
        V_rule%weights(157) =   9.1649207262942799E-03
        V_rule%weights(158) =   1.1569524628097671E-02
        V_rule%weights(159) =   1.1569524628097671E-02
        V_rule%weights(160) =   1.1569524628097671E-02
        V_rule%weights(161) =   1.1569524628097671E-02
        V_rule%weights(162) =   1.1569524628097671E-02
        V_rule%weights(163) =   1.1569524628097671E-02
        V_rule%weights(164) =   1.1761116467609170E-02
        V_rule%weights(165) =   1.1761116467609170E-02
        V_rule%weights(166) =   1.1761116467609170E-02
        V_rule%weights(167) =   1.1761116467609170E-02
        V_rule%weights(168) =   1.1761116467609170E-02
        V_rule%weights(169) =   1.1761116467609170E-02
        V_rule%weights(170) =   1.3824702182165400E-02
        V_rule%weights(171) =   1.3824702182165400E-02
        V_rule%weights(172) =   1.3824702182165400E-02
        V_rule%weights(173) =   1.3824702182165400E-02
        V_rule%weights(174) =   1.3824702182165400E-02
        V_rule%weights(175) =   1.3824702182165400E-02

       V_rule%lambda(  1,1:3) = (/  3.3333333333332998E-01 ,  3.3333333333332998E-01 ,  3.3333333333334003E-01/)
       V_rule%lambda(  2,1:3) = (/  7.3301164327699998E-03 ,  4.9633494178361998E-01 ,  4.9633494178361004E-01/)
       V_rule%lambda(  3,1:3) = (/  4.9633494178361998E-01 ,  4.9633494178361998E-01 ,  7.3301164327600477E-03/)
       V_rule%lambda(  4,1:3) = (/  4.9633494178361998E-01 ,  7.3301164327699998E-03 ,  4.9633494178361004E-01/)
       V_rule%lambda(  5,1:3) = (/  8.2995675802959995E-02 ,  4.5850216209852002E-01 ,  4.5850216209852002E-01/)
       V_rule%lambda(  6,1:3) = (/  4.5850216209852002E-01 ,  4.5850216209852002E-01 ,  8.2995675802959967E-02/)
       V_rule%lambda(  7,1:3) = (/  4.5850216209852002E-01 ,  8.2995675802959995E-02 ,  4.5850216209852002E-01/)
       V_rule%lambda(  8,1:3) = (/  1.5098095612540999E-01 ,  4.2450952193729002E-01 ,  4.2450952193729996E-01/)
       V_rule%lambda(  9,1:3) = (/  4.2450952193729002E-01 ,  4.2450952193729002E-01 ,  1.5098095612541995E-01/)
       V_rule%lambda( 10,1:3) = (/  4.2450952193729002E-01 ,  1.5098095612540999E-01 ,  4.2450952193729996E-01/)
       V_rule%lambda( 11,1:3) = (/  2.3590585989217000E-01 ,  3.8204707005392002E-01 ,  3.8204707005390998E-01/)
       V_rule%lambda( 12,1:3) = (/  3.8204707005392002E-01 ,  3.8204707005392002E-01 ,  2.3590585989215995E-01/)
       V_rule%lambda( 13,1:3) = (/  3.8204707005392002E-01 ,  2.3590585989217000E-01 ,  3.8204707005390998E-01/)
       V_rule%lambda( 14,1:3) = (/  4.3802430840785000E-01 ,  2.8098784579607999E-01 ,  2.8098784579607000E-01/)
       V_rule%lambda( 15,1:3) = (/  2.8098784579607999E-01 ,  2.8098784579607999E-01 ,  4.3802430840784001E-01/)
       V_rule%lambda( 16,1:3) = (/  2.8098784579607999E-01 ,  4.3802430840785000E-01 ,  2.8098784579607000E-01/)
       V_rule%lambda( 17,1:3) = (/  5.4530204829192996E-01 ,  2.2734897585402999E-01 ,  2.2734897585404004E-01/)
       V_rule%lambda( 18,1:3) = (/  2.2734897585402999E-01 ,  2.2734897585402999E-01 ,  5.4530204829193996E-01/)
       V_rule%lambda( 19,1:3) = (/  2.2734897585402999E-01 ,  5.4530204829192996E-01 ,  2.2734897585404004E-01/)
       V_rule%lambda( 20,1:3) = (/  6.5088177698254002E-01 ,  1.7455911150872999E-01 ,  1.7455911150872999E-01/)
       V_rule%lambda( 21,1:3) = (/  1.7455911150872999E-01 ,  1.7455911150872999E-01 ,  6.5088177698254002E-01/)
       V_rule%lambda( 22,1:3) = (/  1.7455911150872999E-01 ,  6.5088177698254002E-01 ,  1.7455911150872999E-01/)
       V_rule%lambda( 23,1:3) = (/  7.5348314559713003E-01 ,  1.2325842720143999E-01 ,  1.2325842720142997E-01/)
       V_rule%lambda( 24,1:3) = (/  1.2325842720143999E-01 ,  1.2325842720143999E-01 ,  7.5348314559712004E-01/)
       V_rule%lambda( 25,1:3) = (/  1.2325842720143999E-01 ,  7.5348314559713003E-01 ,  1.2325842720142997E-01/)
       V_rule%lambda( 26,1:3) = (/  8.3983154221560996E-01 ,  8.0084228892200002E-02 ,  8.0084228892190037E-02/)
       V_rule%lambda( 27,1:3) = (/  8.0084228892200002E-02 ,  8.0084228892200002E-02 ,  8.3983154221559997E-01/)
       V_rule%lambda( 28,1:3) = (/  8.0084228892200002E-02 ,  8.3983154221560996E-01 ,  8.0084228892190037E-02/)
       V_rule%lambda( 29,1:3) = (/  9.0445106518420004E-01 ,  4.7774467407900000E-02 ,  4.7774467407899958E-02/)
       V_rule%lambda( 30,1:3) = (/  4.7774467407900000E-02 ,  4.7774467407900000E-02 ,  9.0445106518420004E-01/)
       V_rule%lambda( 31,1:3) = (/  4.7774467407900000E-02 ,  9.0445106518420004E-01 ,  4.7774467407899958E-02/)
       V_rule%lambda( 32,1:3) = (/  9.5655897063971995E-01 ,  2.1720514680140000E-02 ,  2.1720514680140048E-02/)
       V_rule%lambda( 33,1:3) = (/  2.1720514680140000E-02 ,  2.1720514680140000E-02 ,  9.5655897063971995E-01/)
       V_rule%lambda( 34,1:3) = (/  2.1720514680140000E-02 ,  9.5655897063971995E-01 ,  2.1720514680140048E-02/)
       V_rule%lambda( 35,1:3) = (/  9.9047064476913005E-01 ,  4.7646776154400003E-03 ,  4.7646776154299511E-03/)
       V_rule%lambda( 36,1:3) = (/  4.7646776154400003E-03 ,  4.7646776154400003E-03 ,  9.9047064476911995E-01/)
       V_rule%lambda( 37,1:3) = (/  4.7646776154400003E-03 ,  9.9047064476913005E-01 ,  4.7646776154299511E-03/)
       V_rule%lambda( 38,1:3) = (/  9.2537119334999999E-04 ,  4.1529527091330998E-01 ,  5.8377935789334001E-01/)
       V_rule%lambda( 39,1:3) = (/  4.1529527091330998E-01 ,  5.8377935789334001E-01 ,  9.2537119335001083E-04/)
       V_rule%lambda( 40,1:3) = (/  5.8377935789334001E-01 ,  9.2537119334999999E-04 ,  4.1529527091330998E-01/)
       V_rule%lambda( 41,1:3) = (/  4.1529527091330998E-01 ,  9.2537119334999999E-04 ,  5.8377935789334001E-01/)
       V_rule%lambda( 42,1:3) = (/  5.8377935789334001E-01 ,  4.1529527091330998E-01 ,  9.2537119335001083E-04/)
       V_rule%lambda( 43,1:3) = (/  9.2537119334999999E-04 ,  5.8377935789334001E-01 ,  4.1529527091330998E-01/)
       V_rule%lambda( 44,1:3) = (/  1.3859258555600001E-03 ,  6.1189909785350001E-02 ,  9.3742416435909004E-01/)
       V_rule%lambda( 45,1:3) = (/  6.1189909785350001E-02 ,  9.3742416435909004E-01 ,  1.3859258555599593E-03/)
       V_rule%lambda( 46,1:3) = (/  9.3742416435909004E-01 ,  1.3859258555600001E-03 ,  6.1189909785349959E-02/)
       V_rule%lambda( 47,1:3) = (/  6.1189909785350001E-02 ,  1.3859258555600001E-03 ,  9.3742416435909004E-01/)
       V_rule%lambda( 48,1:3) = (/  9.3742416435909004E-01 ,  6.1189909785350001E-02 ,  1.3859258555599593E-03/)
       V_rule%lambda( 49,1:3) = (/  1.3859258555600001E-03 ,  9.3742416435909004E-01 ,  6.1189909785349959E-02/)
       V_rule%lambda( 50,1:3) = (/  3.6824154559099999E-03 ,  1.6490869013691001E-01 ,  8.3140889440718002E-01/)
       V_rule%lambda( 51,1:3) = (/  1.6490869013691001E-01 ,  8.3140889440718002E-01 ,  3.6824154559099709E-03/)
       V_rule%lambda( 52,1:3) = (/  8.3140889440718002E-01 ,  3.6824154559099999E-03 ,  1.6490869013690998E-01/)
       V_rule%lambda( 53,1:3) = (/  1.6490869013691001E-01 ,  3.6824154559099999E-03 ,  8.3140889440718002E-01/)
       V_rule%lambda( 54,1:3) = (/  8.3140889440718002E-01 ,  1.6490869013691001E-01 ,  3.6824154559099709E-03/)
       V_rule%lambda( 55,1:3) = (/  3.6824154559099999E-03 ,  8.3140889440718002E-01 ,  1.6490869013690998E-01/)
       V_rule%lambda( 56,1:3) = (/  3.9032234241600000E-03 ,  2.5035062231999999E-02 ,  9.7106171434384003E-01/)
       V_rule%lambda( 57,1:3) = (/  2.5035062231999999E-02 ,  9.7106171434384003E-01 ,  3.9032234241599684E-03/)
       V_rule%lambda( 58,1:3) = (/  9.7106171434384003E-01 ,  3.9032234241600000E-03 ,  2.5035062231999968E-02/)
       V_rule%lambda( 59,1:3) = (/  2.5035062231999999E-02 ,  3.9032234241600000E-03 ,  9.7106171434384003E-01/)
       V_rule%lambda( 60,1:3) = (/  9.7106171434384003E-01 ,  2.5035062231999999E-02 ,  3.9032234241599684E-03/)
       V_rule%lambda( 61,1:3) = (/  3.9032234241600000E-03 ,  9.7106171434384003E-01 ,  2.5035062231999968E-02/)
       V_rule%lambda( 62,1:3) = (/  3.2332481550099998E-03 ,  3.0606446515109997E-01 ,  6.9070228669389000E-01/)
       V_rule%lambda( 63,1:3) = (/  3.0606446515109997E-01 ,  6.9070228669389000E-01 ,  3.2332481550100267E-03/)
       V_rule%lambda( 64,1:3) = (/  6.9070228669389000E-01 ,  3.2332481550099998E-03 ,  3.0606446515109997E-01/)
       V_rule%lambda( 65,1:3) = (/  3.0606446515109997E-01 ,  3.2332481550099998E-03 ,  6.9070228669389000E-01/)
       V_rule%lambda( 66,1:3) = (/  6.9070228669389000E-01 ,  3.0606446515109997E-01 ,  3.2332481550100267E-03/)
       V_rule%lambda( 67,1:3) = (/  3.2332481550099998E-03 ,  6.9070228669389000E-01 ,  3.0606446515109997E-01/)
       V_rule%lambda( 68,1:3) = (/  6.4674321122399998E-03 ,  1.0707328373022000E-01 ,  8.8645928415754005E-01/)
       V_rule%lambda( 69,1:3) = (/  1.0707328373022000E-01 ,  8.8645928415754005E-01 ,  6.4674321122399486E-03/)
       V_rule%lambda( 70,1:3) = (/  8.8645928415754005E-01 ,  6.4674321122399998E-03 ,  1.0707328373021995E-01/)
       V_rule%lambda( 71,1:3) = (/  1.0707328373022000E-01 ,  6.4674321122399998E-03 ,  8.8645928415754005E-01/)
       V_rule%lambda( 72,1:3) = (/  8.8645928415754005E-01 ,  1.0707328373022000E-01 ,  6.4674321122399486E-03/)
       V_rule%lambda( 73,1:3) = (/  6.4674321122399998E-03 ,  8.8645928415754005E-01 ,  1.0707328373021995E-01/)
       V_rule%lambda( 74,1:3) = (/  3.2474754913299998E-03 ,  2.2995754934557999E-01 ,  7.6679497516309003E-01/)
       V_rule%lambda( 75,1:3) = (/  2.2995754934557999E-01 ,  7.6679497516308004E-01 ,  3.2474754913399684E-03/)
       V_rule%lambda( 76,1:3) = (/  7.6679497516308004E-01 ,  3.2474754913299998E-03 ,  2.2995754934558996E-01/)
       V_rule%lambda( 77,1:3) = (/  2.2995754934557999E-01 ,  3.2474754913299998E-03 ,  7.6679497516309003E-01/)
       V_rule%lambda( 78,1:3) = (/  7.6679497516308004E-01 ,  2.2995754934557999E-01 ,  3.2474754913399684E-03/)
       V_rule%lambda( 79,1:3) = (/  3.2474754913299998E-03 ,  7.6679497516308004E-01 ,  2.2995754934558996E-01/)
       V_rule%lambda( 80,1:3) = (/  8.6750908067500000E-03 ,  3.3703663330577999E-01 ,  6.5428827588746996E-01/)
       V_rule%lambda( 81,1:3) = (/  3.3703663330577999E-01 ,  6.5428827588745997E-01 ,  8.6750908067600441E-03/)
       V_rule%lambda( 82,1:3) = (/  6.5428827588745997E-01 ,  8.6750908067500000E-03 ,  3.3703663330579003E-01/)
       V_rule%lambda( 83,1:3) = (/  3.3703663330577999E-01 ,  8.6750908067500000E-03 ,  6.5428827588746996E-01/)
       V_rule%lambda( 84,1:3) = (/  6.5428827588745997E-01 ,  3.3703663330577999E-01 ,  8.6750908067600441E-03/)
       V_rule%lambda( 85,1:3) = (/  8.6750908067500000E-03 ,  6.5428827588745997E-01 ,  3.3703663330579003E-01/)
       V_rule%lambda( 86,1:3) = (/  1.5597026467310000E-02 ,  5.6256576182060002E-02 ,  9.2814639735062998E-01/)
       V_rule%lambda( 87,1:3) = (/  5.6256576182060002E-02 ,  9.2814639735062998E-01 ,  1.5597026467310017E-02/)
       V_rule%lambda( 88,1:3) = (/  9.2814639735062998E-01 ,  1.5597026467310000E-02 ,  5.6256576182060022E-02/)
       V_rule%lambda( 89,1:3) = (/  5.6256576182060002E-02 ,  1.5597026467310000E-02 ,  9.2814639735062998E-01/)
       V_rule%lambda( 90,1:3) = (/  9.2814639735062998E-01 ,  5.6256576182060002E-02 ,  1.5597026467310017E-02/)
       V_rule%lambda( 91,1:3) = (/  1.5597026467310000E-02 ,  9.2814639735062998E-01 ,  5.6256576182060022E-02/)
       V_rule%lambda( 92,1:3) = (/  1.7976721253690001E-02 ,  4.0245137521239999E-01 ,  5.7957190353390997E-01/)
       V_rule%lambda( 93,1:3) = (/  4.0245137521239999E-01 ,  5.7957190353390997E-01 ,  1.7976721253690042E-02/)
       V_rule%lambda( 94,1:3) = (/  5.7957190353390997E-01 ,  1.7976721253690001E-02 ,  4.0245137521240004E-01/)
       V_rule%lambda( 95,1:3) = (/  4.0245137521239999E-01 ,  1.7976721253690001E-02 ,  5.7957190353390997E-01/)
       V_rule%lambda( 96,1:3) = (/  5.7957190353390997E-01 ,  4.0245137521239999E-01 ,  1.7976721253690042E-02/)
       V_rule%lambda( 97,1:3) = (/  1.7976721253690001E-02 ,  5.7957190353390997E-01 ,  4.0245137521240004E-01/)
       V_rule%lambda( 98,1:3) = (/  1.7124245353890000E-02 ,  2.4365470201083000E-01 ,  7.3922105263528004E-01/)
       V_rule%lambda( 99,1:3) = (/  2.4365470201083000E-01 ,  7.3922105263528004E-01 ,  1.7124245353889955E-02/)
       V_rule%lambda(100,1:3) = (/  7.3922105263528004E-01 ,  1.7124245353890000E-02 ,  2.4365470201082995E-01/)
       V_rule%lambda(101,1:3) = (/  2.4365470201083000E-01 ,  1.7124245353890000E-02 ,  7.3922105263528004E-01/)
       V_rule%lambda(102,1:3) = (/  7.3922105263528004E-01 ,  2.4365470201083000E-01 ,  1.7124245353889955E-02/)
       V_rule%lambda(103,1:3) = (/  1.7124245353890000E-02 ,  7.3922105263528004E-01 ,  2.4365470201082995E-01/)
       V_rule%lambda(104,1:3) = (/  2.2883405346580000E-02 ,  1.6538958561452999E-01 ,  8.1172700903889006E-01/)
       V_rule%lambda(105,1:3) = (/  1.6538958561452999E-01 ,  8.1172700903887995E-01 ,  2.2883405346590058E-02/)
       V_rule%lambda(106,1:3) = (/  8.1172700903887995E-01 ,  2.2883405346580000E-02 ,  1.6538958561454004E-01/)
       V_rule%lambda(107,1:3) = (/  1.6538958561452999E-01 ,  2.2883405346580000E-02 ,  8.1172700903889006E-01/)
       V_rule%lambda(108,1:3) = (/  8.1172700903887995E-01 ,  1.6538958561452999E-01 ,  2.2883405346590058E-02/)
       V_rule%lambda(109,1:3) = (/  2.2883405346580000E-02 ,  8.1172700903887995E-01 ,  1.6538958561454004E-01/)
       V_rule%lambda(110,1:3) = (/  3.2737597287770002E-02 ,  9.9301874495849998E-02 ,  8.6796052821638003E-01/)
       V_rule%lambda(111,1:3) = (/  9.9301874495849998E-02 ,  8.6796052821639003E-01 ,  3.2737597287759976E-02/)
       V_rule%lambda(112,1:3) = (/  8.6796052821639003E-01 ,  3.2737597287770002E-02 ,  9.9301874495839965E-02/)
       V_rule%lambda(113,1:3) = (/  9.9301874495849998E-02 ,  3.2737597287770002E-02 ,  8.6796052821638003E-01/)
       V_rule%lambda(114,1:3) = (/  8.6796052821639003E-01 ,  9.9301874495849998E-02 ,  3.2737597287759976E-02/)
       V_rule%lambda(115,1:3) = (/  3.2737597287770002E-02 ,  8.6796052821639003E-01 ,  9.9301874495839965E-02/)
       V_rule%lambda(116,1:3) = (/  3.3821012342340001E-02 ,  3.0847833306904998E-01 ,  6.5770065458861005E-01/)
       V_rule%lambda(117,1:3) = (/  3.0847833306904998E-01 ,  6.5770065458860005E-01 ,  3.3821012342349965E-02/)
       V_rule%lambda(118,1:3) = (/  6.5770065458860005E-01 ,  3.3821012342340001E-02 ,  3.0847833306905992E-01/)
       V_rule%lambda(119,1:3) = (/  3.0847833306904998E-01 ,  3.3821012342340001E-02 ,  6.5770065458861005E-01/)
       V_rule%lambda(120,1:3) = (/  6.5770065458860005E-01 ,  3.0847833306904998E-01 ,  3.3821012342349965E-02/)
       V_rule%lambda(121,1:3) = (/  3.3821012342340001E-02 ,  6.5770065458860005E-01 ,  3.0847833306905992E-01/)
       V_rule%lambda(122,1:3) = (/  3.5547614460019999E-02 ,  4.6066831859210999E-01 ,  5.0378406694787004E-01/)
       V_rule%lambda(123,1:3) = (/  4.6066831859210999E-01 ,  5.0378406694787004E-01 ,  3.5547614460019972E-02/)
       V_rule%lambda(124,1:3) = (/  5.0378406694787004E-01 ,  3.5547614460019999E-02 ,  4.6066831859210999E-01/)
       V_rule%lambda(125,1:3) = (/  4.6066831859210999E-01 ,  3.5547614460019999E-02 ,  5.0378406694787004E-01/)
       V_rule%lambda(126,1:3) = (/  5.0378406694787004E-01 ,  4.6066831859210999E-01 ,  3.5547614460019972E-02/)
       V_rule%lambda(127,1:3) = (/  3.5547614460019999E-02 ,  5.0378406694787004E-01 ,  4.6066831859210999E-01/)
       V_rule%lambda(128,1:3) = (/  5.0539790306870003E-02 ,  2.1881529945393000E-01 ,  7.3064491023919997E-01/)
       V_rule%lambda(129,1:3) = (/  2.1881529945393000E-01 ,  7.3064491023919997E-01 ,  5.0539790306870030E-02/)
       V_rule%lambda(130,1:3) = (/  7.3064491023919997E-01 ,  5.0539790306870003E-02 ,  2.1881529945393002E-01/)
       V_rule%lambda(131,1:3) = (/  2.1881529945393000E-01 ,  5.0539790306870003E-02 ,  7.3064491023919997E-01/)
       V_rule%lambda(132,1:3) = (/  7.3064491023919997E-01 ,  2.1881529945393000E-01 ,  5.0539790306870030E-02/)
       V_rule%lambda(133,1:3) = (/  5.0539790306870003E-02 ,  7.3064491023919997E-01 ,  2.1881529945393002E-01/)
       V_rule%lambda(134,1:3) = (/  5.7014714915730000E-02 ,  3.7920955156026998E-01 ,  5.6377573352400001E-01/)
       V_rule%lambda(135,1:3) = (/  3.7920955156026998E-01 ,  5.6377573352399002E-01 ,  5.7014714915740006E-02/)
       V_rule%lambda(136,1:3) = (/  5.6377573352399002E-01 ,  5.7014714915730000E-02 ,  3.7920955156027997E-01/)
       V_rule%lambda(137,1:3) = (/  3.7920955156026998E-01 ,  5.7014714915730000E-02 ,  5.6377573352400001E-01/)
       V_rule%lambda(138,1:3) = (/  5.6377573352399002E-01 ,  3.7920955156026998E-01 ,  5.7014714915740006E-02/)
       V_rule%lambda(139,1:3) = (/  5.7014714915730000E-02 ,  5.6377573352399002E-01 ,  3.7920955156027997E-01/)
       V_rule%lambda(140,1:3) = (/  6.4152806421199998E-02 ,  1.4296081941819000E-01 ,  7.9288637416061003E-01/)
       V_rule%lambda(141,1:3) = (/  1.4296081941819000E-01 ,  7.9288637416061003E-01 ,  6.4152806421199970E-02/)
       V_rule%lambda(142,1:3) = (/  7.9288637416061003E-01 ,  6.4152806421199998E-02 ,  1.4296081941818997E-01/)
       V_rule%lambda(143,1:3) = (/  1.4296081941819000E-01 ,  6.4152806421199998E-02 ,  7.9288637416061003E-01/)
       V_rule%lambda(144,1:3) = (/  7.9288637416061003E-01 ,  1.4296081941819000E-01 ,  6.4152806421199970E-02/)
       V_rule%lambda(145,1:3) = (/  6.4152806421199998E-02 ,  7.9288637416061003E-01 ,  1.4296081941818997E-01/)
       V_rule%lambda(146,1:3) = (/  8.0501148287629998E-02 ,  2.8373128210592002E-01 ,  6.3576756960644998E-01/)
       V_rule%lambda(147,1:3) = (/  2.8373128210592002E-01 ,  6.3576756960644998E-01 ,  8.0501148287629998E-02/)
       V_rule%lambda(148,1:3) = (/  6.3576756960644998E-01 ,  8.0501148287629998E-02 ,  2.8373128210592002E-01/)
       V_rule%lambda(149,1:3) = (/  2.8373128210592002E-01 ,  8.0501148287629998E-02 ,  6.3576756960644998E-01/)
       V_rule%lambda(150,1:3) = (/  6.3576756960644998E-01 ,  2.8373128210592002E-01 ,  8.0501148287629998E-02/)
       V_rule%lambda(151,1:3) = (/  8.0501148287629998E-02 ,  6.3576756960644998E-01 ,  2.8373128210592002E-01/)
       V_rule%lambda(152,1:3) = (/  1.0436706813453001E-01 ,  1.9673744100443999E-01 ,  6.9889549086102998E-01/)
       V_rule%lambda(153,1:3) = (/  1.9673744100443999E-01 ,  6.9889549086102998E-01 ,  1.0436706813453003E-01/)
       V_rule%lambda(154,1:3) = (/  6.9889549086102998E-01 ,  1.0436706813453001E-01 ,  1.9673744100444002E-01/)
       V_rule%lambda(155,1:3) = (/  1.9673744100443999E-01 ,  1.0436706813453001E-01 ,  6.9889549086102998E-01/)
       V_rule%lambda(156,1:3) = (/  6.9889549086102998E-01 ,  1.9673744100443999E-01 ,  1.0436706813453003E-01/)
       V_rule%lambda(157,1:3) = (/  1.0436706813453001E-01 ,  6.9889549086102998E-01 ,  1.9673744100444002E-01/)
       V_rule%lambda(158,1:3) = (/  1.1384489442875000E-01 ,  3.5588914121165999E-01 ,  5.3026596435959006E-01/)
       V_rule%lambda(159,1:3) = (/  3.5588914121165999E-01 ,  5.3026596435958995E-01 ,  1.1384489442875007E-01/)
       V_rule%lambda(160,1:3) = (/  5.3026596435958995E-01 ,  1.1384489442875000E-01 ,  3.5588914121166004E-01/)
       V_rule%lambda(161,1:3) = (/  3.5588914121165999E-01 ,  1.1384489442875000E-01 ,  5.3026596435959006E-01/)
       V_rule%lambda(162,1:3) = (/  5.3026596435958995E-01 ,  3.5588914121165999E-01 ,  1.1384489442875007E-01/)
       V_rule%lambda(163,1:3) = (/  1.1384489442875000E-01 ,  5.3026596435958995E-01 ,  3.5588914121166004E-01/)
       V_rule%lambda(164,1:3) = (/  1.4536348771551999E-01 ,  2.5981868535190999E-01 ,  5.9481782693257002E-01/)
       V_rule%lambda(165,1:3) = (/  2.5981868535190999E-01 ,  5.9481782693256002E-01 ,  1.4536348771552998E-01/)
       V_rule%lambda(166,1:3) = (/  5.9481782693256002E-01 ,  1.4536348771551999E-01 ,  2.5981868535191999E-01/)
       V_rule%lambda(167,1:3) = (/  2.5981868535190999E-01 ,  1.4536348771551999E-01 ,  5.9481782693257002E-01/)
       V_rule%lambda(168,1:3) = (/  5.9481782693256002E-01 ,  2.5981868535190999E-01 ,  1.4536348771552998E-01/)
       V_rule%lambda(169,1:3) = (/  1.4536348771551999E-01 ,  5.9481782693256002E-01 ,  2.5981868535191999E-01/)
       V_rule%lambda(170,1:3) = (/  1.8994565282198000E-01 ,  3.2192318123129998E-01 ,  4.8813116594672001E-01/)
       V_rule%lambda(171,1:3) = (/  3.2192318123129998E-01 ,  4.8813116594672001E-01 ,  1.8994565282198000E-01/)
       V_rule%lambda(172,1:3) = (/  4.8813116594672001E-01 ,  1.8994565282198000E-01 ,  3.2192318123129998E-01/)
       V_rule%lambda(173,1:3) = (/  3.2192318123129998E-01 ,  1.8994565282198000E-01 ,  4.8813116594672001E-01/)
       V_rule%lambda(174,1:3) = (/  4.8813116594672001E-01 ,  3.2192318123129998E-01 ,  1.8994565282198000E-01/)
       V_rule%lambda(175,1:3) = (/  1.8994565282198000E-01 ,  4.8813116594672001E-01 ,  3.2192318123129998E-01/)

    case(24:)
       print*,'Warning !! Volume quadrature rule of degree ', Qdeg, &
            ' is not implemented'
       print*,'           the maximal default quadrature rule of degree', &
            ' 23 is used'
    end select


  end subroutine Create_volume_rule

  !> generates one volume (Dunavant triangular) quadrature rule
  subroutine Create_volume_rule_OLD(V_rule, Qdeg)
    type(volume_rule), intent(inout) :: V_rule
    integer, intent(in) :: Qdeg

    V_rule%Qdeg = Qdeg
    select case (Qdeg)
    case(:1)
       V_rule%Qdof = 1
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))
       V_rule%weights(:) = 1.
       V_rule%lambda(1,1:3) = 1./3
    case(2)
       V_rule%Qdof = 3
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))
       V_rule%weights(1:3)  = 1./3
       V_rule%lambda(1,1:3) = (/4./6, 1./6, 1./6 /)
       V_rule%lambda(2,1:3) = (/1./6, 4./6, 1./6 /)
       V_rule%lambda(3,1:3) = (/1./6, 1./6, 4./6 /)

      case( 3)
       V_rule%Qdof =  4
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))
       V_rule%weights( 1) = -0.5625000000000000E+00
       V_rule%weights( 2: 4) =  0.5208333333333330E+00

       V_rule%lambda( 1,1:3) = (/ 0.3333333333333330E+00, 0.3333333333333330E+00, 0.3333333333333340E+00/)
       V_rule%lambda( 2,1:3) = (/ 0.6000000000000000E+00, 0.2000000000000000E+00, 0.2000000000000000E+00/)
       V_rule%lambda( 3,1:3) = (/ 0.2000000000000000E+00, 0.2000000000000000E+00, 0.6000000000000000E+00/)
       V_rule%lambda( 4,1:3) = (/ 0.2000000000000000E+00, 0.6000000000000000E+00, 0.2000000000000000E+00/)


    case(4)
       V_rule%Qdof = 6
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))
       V_rule%weights(1:3) = 0.223381589678011
       V_rule%weights(4:6) = 0.109951743655322

       V_rule%lambda(1,1:3) = (/0.108103018168070, 0.445948490915965, 0.445948490915965 /)
       V_rule%lambda(2,1:3) = (/0.445948490915965, 0.108103018168070, 0.445948490915965 /)
       V_rule%lambda(3,1:3) = (/0.445948490915965, 0.445948490915965, 0.108103018168070 /)
       V_rule%lambda(4,1:3) = (/ 0.816847572980459, 0.091576213509771, 0.091576213509771/)
       V_rule%lambda(5,1:3) = (/ 0.091576213509771, 0.816847572980459, 0.091576213509771/)
       V_rule%lambda(6,1:3) = (/ 0.091576213509771, 0.091576213509771, 0.816847572980459/)

    case(5)
       V_rule%Qdof = 7
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))
       V_rule%weights(1) =  0.225
       V_rule%weights(2:4) =  0.132394152788506
       V_rule%weights(5:7) =  0.125939180544827

       V_rule%lambda(1,1:3) = (/0.333333333333333, 0.333333333333333, 0.333333333333333 /)
       V_rule%lambda(2,1:3) = (/0.059715871789770, 0.470142064105115, 0.470142064105115 /)
       V_rule%lambda(3,1:3) = (/0.470142064105115, 0.059715871789770, 0.470142064105115 /)
       V_rule%lambda(4,1:3) = (/0.470142064105115, 0.470142064105115, 0.059715871789770 /)
       V_rule%lambda(5,1:3) = (/0.797426985353087, 0.101286507323456, 0.101286507323456 /)
       V_rule%lambda(6,1:3) = (/0.101286507323456, 0.797426985353087, 0.101286507323456 /)
       V_rule%lambda(7,1:3) = (/0.101286507323456, 0.101286507323456, 0.797426985353087 /)
    case(6)
       V_rule%Qdof = 12
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))
       V_rule%weights(1:3) =  0.116786275726379
       V_rule%weights(4:6) =  0.050844906370207
       V_rule%weights(7:12) =  0.082851075618374

       V_rule%lambda(1,1:3) = (/0.501426509658179, 0.249286745170910, 0.249286745170910/)
       V_rule%lambda(2,1:3) = (/0.249286745170910, 0.501426509658179, 0.249286745170910/)
       V_rule%lambda(3,1:3) = (/0.249286745170910, 0.249286745170910, 0.501426509658179/)
       V_rule%lambda(4,1:3) = (/0.873821971016996, 0.063089014491502, 0.063089014491502/)
       V_rule%lambda(5,1:3) = (/0.063089014491502, 0.873821971016996, 0.063089014491502/)
       V_rule%lambda(6,1:3) = (/0.063089014491502, 0.063089014491502, 0.873821971016996/)
       V_rule%lambda(7,1:3) = (/0.053145049844817, 0.310352451033784, 0.636502499121399/)
       V_rule%lambda(8,1:3) = (/0.053145049844817, 0.636502499121399, 0.310352451033784/)
       V_rule%lambda(9,1:3) = (/0.310352451033784, 0.053145049844817, 0.636502499121399/)
       V_rule%lambda(10,1:3)= (/0.310352451033784, 0.636502499121399, 0.053145049844817/)
       V_rule%lambda(11,1:3)= (/0.636502499121399, 0.053145049844817, 0.310352451033784/)
       V_rule%lambda(12,1:3)= (/0.636502499121399, 0.310352451033784, 0.053145049844817/)
    case(7)
       V_rule%Qdof = 13
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))
       V_rule%weights(1) =  -0.149570044467682
       V_rule%weights(2:4) =  0.175615257433208
       V_rule%weights(5:7) =  0.053347235608838
       V_rule%weights(8:13) =  0.077113760890257

       V_rule%lambda(1,1:3) = (/0.333333333333333, 0.333333333333333, 0.333333333333333/)
       V_rule%lambda(2,1:3) = (/0.479308067841920, 0.260345966079040, 0.260345966079040/)
       V_rule%lambda(3,1:3) = (/0.260345966079040, 0.479308067841920, 0.260345966079040/)
       V_rule%lambda(4,1:3) = (/0.260345966079040, 0.260345966079040, 0.479308067841920/)
       V_rule%lambda(5,1:3) = (/0.869739794195568, 0.065130102902216, 0.065130102902216/)
       V_rule%lambda(6,1:3) = (/0.065130102902216, 0.869739794195568, 0.065130102902216/)
       V_rule%lambda(7,1:3) = (/0.065130102902216, 0.065130102902216, 0.869739794195568/)
       V_rule%lambda(8,1:3) = (/0.048690315425316, 0.312865496004874, 0.638444188569810/)
       V_rule%lambda(9,1:3) = (/0.048690315425316, 0.638444188569810, 0.312865496004874/)
       V_rule%lambda(10,1:3)= (/0.312865496004874, 0.048690315425316, 0.638444188569810/)
       V_rule%lambda(11,1:3)= (/0.312865496004874, 0.638444188569810, 0.048690315425316/)
       V_rule%lambda(12,1:3)= (/0.638444188569810, 0.048690315425316, 0.312865496004874/)
       V_rule%lambda(13,1:3)= (/0.638444188569810, 0.312865496004874, 0.048690315425316/)

    case(8)
       V_rule%Qdof = 16
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))
       V_rule%weights(1) =  0.144315607677787
       V_rule%weights(2:4) =  0.095091634267285
       V_rule%weights(5:7) =  0.103217370534718
       V_rule%weights(8:10) =  0.032458497623198
       V_rule%weights(11:16) =  0.027230314174435

       V_rule%lambda(1,1:3) = (/0.333333333333333, 0.333333333333333, 0.333333333333333/)
       V_rule%lambda(2,1:3) = (/0.081414823414554, 0.459292588292723, 0.459292588292723/)
       V_rule%lambda(3,1:3) = (/0.459292588292723, 0.081414823414554, 0.459292588292723/)
       V_rule%lambda(4,1:3) = (/0.459292588292723, 0.459292588292723, 0.081414823414554/)
       V_rule%lambda(5,1:3) = (/0.658861384496480, 0.170569307751760, 0.170569307751760/)
       V_rule%lambda(6,1:3) = (/0.170569307751760, 0.658861384496480, 0.170569307751760/)
       V_rule%lambda(7,1:3) = (/0.170569307751760, 0.170569307751760, 0.658861384496480/)
       V_rule%lambda(8,1:3) = (/0.898905543365938, 0.050547228317031, 0.050547228317031/)
       V_rule%lambda(9,1:3) = (/0.050547228317031, 0.898905543365938, 0.050547228317031/)
       V_rule%lambda(10,1:3) = (/0.050547228317031, 0.050547228317031, 0.898905543365938/)
       V_rule%lambda(11,1:3) = (/0.008394777409958, 0.263112829634638, 0.728492392955404/)
       V_rule%lambda(12,1:3) = (/0.008394777409958, 0.728492392955404, 0.263112829634638/)
       V_rule%lambda(13,1:3) = (/0.263112829634638, 0.008394777409958, 0.728492392955404/)
       V_rule%lambda(14,1:3) = (/0.728492392955404, 0.008394777409958, 0.263112829634638/)
       V_rule%lambda(15,1:3) = (/0.263112829634638, 0.728492392955404, 0.008394777409958/)
       V_rule%lambda(16,1:3) = (/0.728492392955404, 0.263112829634638, 0.008394777409958/)

      case( 9)
       V_rule%Qdof = 19
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))
       V_rule%weights( 1) =  0.9713579628279900E-01
       V_rule%weights( 2: 4) =  0.3133470022713900E-01
       V_rule%weights( 5: 7) =  0.7782754100477400E-01
       V_rule%weights( 8:10) =  0.7964773892721000E-01
       V_rule%weights(11:13) =  0.2557767565869800E-01
       V_rule%weights(14:19) =  0.4328353937728900E-01

       V_rule%lambda( 1,1:3) = (/ 0.3333333333333330E+00, 0.3333333333333330E+00, 0.3333333333333340E+00/)
       V_rule%lambda( 2,1:3) = (/ 0.2063496160252500E-01, 0.4896825191987380E+00, 0.4896825191987370E+00/)
       V_rule%lambda( 3,1:3) = (/ 0.4896825191987380E+00, 0.4896825191987380E+00, 0.2063496160252398E-01/)
       V_rule%lambda( 4,1:3) = (/ 0.4896825191987380E+00, 0.2063496160252500E-01, 0.4896825191987370E+00/)
       V_rule%lambda( 5,1:3) = (/ 0.1258208170141270E+00, 0.4370895914929370E+00, 0.4370895914929360E+00/)
       V_rule%lambda( 6,1:3) = (/ 0.4370895914929370E+00, 0.4370895914929370E+00, 0.1258208170141260E+00/)
       V_rule%lambda( 7,1:3) = (/ 0.4370895914929370E+00, 0.1258208170141270E+00, 0.4370895914929360E+00/)
       V_rule%lambda( 8,1:3) = (/ 0.6235929287619350E+00, 0.1882035356190330E+00, 0.1882035356190320E+00/)
       V_rule%lambda( 9,1:3) = (/ 0.1882035356190330E+00, 0.1882035356190330E+00, 0.6235929287619340E+00/)
       V_rule%lambda(10,1:3) = (/ 0.1882035356190330E+00, 0.6235929287619350E+00, 0.1882035356190320E+00/)
       V_rule%lambda(11,1:3) = (/ 0.9105409732110950E+00, 0.4472951339445300E-01, 0.4472951339445208E-01/)
       V_rule%lambda(12,1:3) = (/ 0.4472951339445300E-01, 0.4472951339445300E-01, 0.9105409732110941E+00/)
       V_rule%lambda(13,1:3) = (/ 0.4472951339445300E-01, 0.9105409732110950E+00, 0.4472951339445208E-01/)
       V_rule%lambda(14,1:3) = (/ 0.3683841205473600E-01, 0.2219629891607660E+00, 0.7411985987844980E+00/)
       V_rule%lambda(15,1:3) = (/ 0.2219629891607660E+00, 0.7411985987844980E+00, 0.3683841205473604E-01/)
       V_rule%lambda(16,1:3) = (/ 0.7411985987844980E+00, 0.3683841205473600E-01, 0.2219629891607660E+00/)
       V_rule%lambda(17,1:3) = (/ 0.2219629891607660E+00, 0.3683841205473600E-01, 0.7411985987844980E+00/)
       V_rule%lambda(18,1:3) = (/ 0.7411985987844980E+00, 0.2219629891607660E+00, 0.3683841205473604E-01/)
       V_rule%lambda(19,1:3) = (/ 0.3683841205473600E-01, 0.7411985987844980E+00, 0.2219629891607660E+00/)

      case(10)
       V_rule%Qdof = 25
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))
       V_rule%weights( 1) =  0.9081799038275400E-01
       V_rule%weights( 2: 4) =  0.3672595775646700E-01
       V_rule%weights( 5: 7) =  0.4532105943552800E-01
       V_rule%weights( 8:13) =  0.7275791684542000E-01
       V_rule%weights(14:19) =  0.2832724253105700E-01
       V_rule%weights(20:25) =  0.9421666963733000E-02

       V_rule%lambda( 1,1:3) = (/ 0.3333333333333330E+00, 0.3333333333333330E+00, 0.3333333333333340E+00/)
       V_rule%lambda( 2,1:3) = (/ 0.2884473323268500E-01, 0.4855776333836570E+00, 0.4855776333836580E+00/)
       V_rule%lambda( 3,1:3) = (/ 0.4855776333836570E+00, 0.4855776333836570E+00, 0.2884473323268599E-01/)
       V_rule%lambda( 4,1:3) = (/ 0.4855776333836570E+00, 0.2884473323268500E-01, 0.4855776333836580E+00/)
       V_rule%lambda( 5,1:3) = (/ 0.7810368490299260E+00, 0.1094815754850370E+00, 0.1094815754850370E+00/)
       V_rule%lambda( 6,1:3) = (/ 0.1094815754850370E+00, 0.1094815754850370E+00, 0.7810368490299260E+00/)
       V_rule%lambda( 7,1:3) = (/ 0.1094815754850370E+00, 0.7810368490299260E+00, 0.1094815754850370E+00/)
       V_rule%lambda( 8,1:3) = (/ 0.1417072194148800E+00, 0.3079398387641210E+00, 0.5503529418209990E+00/)
       V_rule%lambda( 9,1:3) = (/ 0.3079398387641210E+00, 0.5503529418209990E+00, 0.1417072194148800E+00/)
       V_rule%lambda(10,1:3) = (/ 0.5503529418209990E+00, 0.1417072194148800E+00, 0.3079398387641210E+00/)
       V_rule%lambda(11,1:3) = (/ 0.3079398387641210E+00, 0.1417072194148800E+00, 0.5503529418209990E+00/)
       V_rule%lambda(12,1:3) = (/ 0.5503529418209990E+00, 0.3079398387641210E+00, 0.1417072194148800E+00/)
       V_rule%lambda(13,1:3) = (/ 0.1417072194148800E+00, 0.5503529418209990E+00, 0.3079398387641210E+00/)
       V_rule%lambda(14,1:3) = (/ 0.2500353476268600E-01, 0.2466725606399030E+00, 0.7283239045974110E+00/)
       V_rule%lambda(15,1:3) = (/ 0.2466725606399030E+00, 0.7283239045974110E+00, 0.2500353476268602E-01/)
       V_rule%lambda(16,1:3) = (/ 0.7283239045974110E+00, 0.2500353476268600E-01, 0.2466725606399029E+00/)
       V_rule%lambda(17,1:3) = (/ 0.2466725606399030E+00, 0.2500353476268600E-01, 0.7283239045974110E+00/)
       V_rule%lambda(18,1:3) = (/ 0.7283239045974110E+00, 0.2466725606399030E+00, 0.2500353476268602E-01/)
       V_rule%lambda(19,1:3) = (/ 0.2500353476268600E-01, 0.7283239045974110E+00, 0.2466725606399029E+00/)
       V_rule%lambda(20,1:3) = (/ 0.9540815400299000E-02, 0.6680325101220000E-01, 0.9236559335875010E+00/)
       V_rule%lambda(21,1:3) = (/ 0.6680325101220000E-01, 0.9236559335875000E+00, 0.9540815400300051E-02/)
       V_rule%lambda(22,1:3) = (/ 0.9236559335875000E+00, 0.9540815400299000E-02, 0.6680325101220097E-01/)
       V_rule%lambda(23,1:3) = (/ 0.6680325101220000E-01, 0.9540815400299000E-02, 0.9236559335875010E+00/)
       V_rule%lambda(24,1:3) = (/ 0.9236559335875000E+00, 0.6680325101220000E-01, 0.9540815400300051E-02/)
       V_rule%lambda(25,1:3) = (/ 0.9540815400299000E-02, 0.9236559335875000E+00, 0.6680325101220097E-01/)

      case(11)
       V_rule%Qdof = 27
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))
       V_rule%weights( 1: 3) =  0.9270063289610000E-03
       V_rule%weights( 4: 6) =  0.7714953491481299E-01
       V_rule%weights( 7: 9) =  0.5932297738077400E-01
       V_rule%weights(10:12) =  0.3618454050341800E-01
       V_rule%weights(13:15) =  0.1365973100267800E-01
       V_rule%weights(16:21) =  0.5233711196220400E-01
       V_rule%weights(22:27) =  0.2070765963914100E-01

       V_rule%lambda( 1,1:3) = (/-0.6922209654151699E-01, 0.5346110482707580E+00, 0.5346110482707590E+00/)
       V_rule%lambda( 2,1:3) = (/ 0.5346110482707580E+00, 0.5346110482707580E+00,-0.6922209654151601E-01/)
       V_rule%lambda( 3,1:3) = (/ 0.5346110482707580E+00,-0.6922209654151699E-01, 0.5346110482707590E+00/)
       V_rule%lambda( 4,1:3) = (/ 0.2020613940682900E+00, 0.3989693029658550E+00, 0.3989693029658550E+00/)
       V_rule%lambda( 5,1:3) = (/ 0.3989693029658550E+00, 0.3989693029658550E+00, 0.2020613940682900E+00/)
       V_rule%lambda( 6,1:3) = (/ 0.3989693029658550E+00, 0.2020613940682900E+00, 0.3989693029658550E+00/)
       V_rule%lambda( 7,1:3) = (/ 0.5933801991374350E+00, 0.2033099004312820E+00, 0.2033099004312831E+00/)
       V_rule%lambda( 8,1:3) = (/ 0.2033099004312820E+00, 0.2033099004312820E+00, 0.5933801991374360E+00/)
       V_rule%lambda( 9,1:3) = (/ 0.2033099004312820E+00, 0.5933801991374350E+00, 0.2033099004312831E+00/)
       V_rule%lambda(10,1:3) = (/ 0.7612981754348370E+00, 0.1193509122825810E+00, 0.1193509122825820E+00/)
       V_rule%lambda(11,1:3) = (/ 0.1193509122825810E+00, 0.1193509122825810E+00, 0.7612981754348380E+00/)
       V_rule%lambda(12,1:3) = (/ 0.1193509122825810E+00, 0.7612981754348370E+00, 0.1193509122825820E+00/)
       V_rule%lambda(13,1:3) = (/ 0.9352701037774480E+00, 0.3236494811127600E-01, 0.3236494811127599E-01/)
       V_rule%lambda(14,1:3) = (/ 0.3236494811127600E-01, 0.3236494811127600E-01, 0.9352701037774480E+00/)
       V_rule%lambda(15,1:3) = (/ 0.3236494811127600E-01, 0.9352701037774480E+00, 0.3236494811127599E-01/)
       V_rule%lambda(16,1:3) = (/ 0.5017813831049500E-01, 0.3566206482612930E+00, 0.5932012134282120E+00/)
       V_rule%lambda(17,1:3) = (/ 0.3566206482612930E+00, 0.5932012134282130E+00, 0.5017813831049400E-01/)
       V_rule%lambda(18,1:3) = (/ 0.5932012134282130E+00, 0.5017813831049500E-01, 0.3566206482612920E+00/)
       V_rule%lambda(19,1:3) = (/ 0.3566206482612930E+00, 0.5017813831049500E-01, 0.5932012134282120E+00/)
       V_rule%lambda(20,1:3) = (/ 0.5932012134282130E+00, 0.3566206482612930E+00, 0.5017813831049400E-01/)
       V_rule%lambda(21,1:3) = (/ 0.5017813831049500E-01, 0.5932012134282130E+00, 0.3566206482612920E+00/)
       V_rule%lambda(22,1:3) = (/ 0.2102201653616600E-01, 0.1714889803040420E+00, 0.8074890031597920E+00/)
       V_rule%lambda(23,1:3) = (/ 0.1714889803040420E+00, 0.8074890031597920E+00, 0.2102201653616598E-01/)
       V_rule%lambda(24,1:3) = (/ 0.8074890031597920E+00, 0.2102201653616600E-01, 0.1714889803040420E+00/)
       V_rule%lambda(25,1:3) = (/ 0.1714889803040420E+00, 0.2102201653616600E-01, 0.8074890031597920E+00/)
       V_rule%lambda(26,1:3) = (/ 0.8074890031597920E+00, 0.1714889803040420E+00, 0.2102201653616598E-01/)
       V_rule%lambda(27,1:3) = (/ 0.2102201653616600E-01, 0.8074890031597920E+00, 0.1714889803040420E+00/)

      case(12)
       V_rule%Qdof = 33
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))
       V_rule%weights( 1: 3) =  0.2573106644045500E-01
       V_rule%weights( 4: 6) =  0.4369254453803800E-01
       V_rule%weights( 7: 9) =  0.6285822421788501E-01
       V_rule%weights(10:12) =  0.3479611293070900E-01
       V_rule%weights(13:15) =  0.6166261051559000E-02
       V_rule%weights(16:21) =  0.4037155776638100E-01
       V_rule%weights(22:27) =  0.2235677320230300E-01
       V_rule%weights(28:33) =  0.1731623110865900E-01

       V_rule%lambda( 1,1:3) = (/ 0.2356522045239000E-01, 0.4882173897738050E+00, 0.4882173897738050E+00/)
       V_rule%lambda( 2,1:3) = (/ 0.4882173897738050E+00, 0.4882173897738050E+00, 0.2356522045238996E-01/)
       V_rule%lambda( 3,1:3) = (/ 0.4882173897738050E+00, 0.2356522045239000E-01, 0.4882173897738050E+00/)
       V_rule%lambda( 4,1:3) = (/ 0.1205512154110790E+00, 0.4397243922944600E+00, 0.4397243922944610E+00/)
       V_rule%lambda( 5,1:3) = (/ 0.4397243922944600E+00, 0.4397243922944600E+00, 0.1205512154110800E+00/)
       V_rule%lambda( 6,1:3) = (/ 0.4397243922944600E+00, 0.1205512154110790E+00, 0.4397243922944610E+00/)
       V_rule%lambda( 7,1:3) = (/ 0.4575792299757680E+00, 0.2712103850121160E+00, 0.2712103850121160E+00/)
       V_rule%lambda( 8,1:3) = (/ 0.2712103850121160E+00, 0.2712103850121160E+00, 0.4575792299757681E+00/)
       V_rule%lambda( 9,1:3) = (/ 0.2712103850121160E+00, 0.4575792299757680E+00, 0.2712103850121160E+00/)
       V_rule%lambda(10,1:3) = (/ 0.7448477089168281E+00, 0.1275761455415860E+00, 0.1275761455415859E+00/)
       V_rule%lambda(11,1:3) = (/ 0.1275761455415860E+00, 0.1275761455415860E+00, 0.7448477089168279E+00/)
       V_rule%lambda(12,1:3) = (/ 0.1275761455415860E+00, 0.7448477089168281E+00, 0.1275761455415859E+00/)
       V_rule%lambda(13,1:3) = (/ 0.9573652990935790E+00, 0.2131735045321000E-01, 0.2131735045321093E-01/)
       V_rule%lambda(14,1:3) = (/ 0.2131735045321000E-01, 0.2131735045321000E-01, 0.9573652990935800E+00/)
       V_rule%lambda(15,1:3) = (/ 0.2131735045321000E-01, 0.9573652990935790E+00, 0.2131735045321093E-01/)
       V_rule%lambda(16,1:3) = (/ 0.1153434945346980E+00, 0.2757132696855140E+00, 0.6089432357797879E+00/)
       V_rule%lambda(17,1:3) = (/ 0.2757132696855140E+00, 0.6089432357797880E+00, 0.1153434945346979E+00/)
       V_rule%lambda(18,1:3) = (/ 0.6089432357797880E+00, 0.1153434945346980E+00, 0.2757132696855140E+00/)
       V_rule%lambda(19,1:3) = (/ 0.2757132696855140E+00, 0.1153434945346980E+00, 0.6089432357797879E+00/)
       V_rule%lambda(20,1:3) = (/ 0.6089432357797880E+00, 0.2757132696855140E+00, 0.1153434945346979E+00/)
       V_rule%lambda(21,1:3) = (/ 0.1153434945346980E+00, 0.6089432357797880E+00, 0.2757132696855140E+00/)
       V_rule%lambda(22,1:3) = (/ 0.2283833222225700E-01, 0.2813255809899400E+00, 0.6958360867878031E+00/)
       V_rule%lambda(23,1:3) = (/ 0.2813255809899400E+00, 0.6958360867878030E+00, 0.2283833222225706E-01/)
       V_rule%lambda(24,1:3) = (/ 0.6958360867878030E+00, 0.2283833222225700E-01, 0.2813255809899401E+00/)
       V_rule%lambda(25,1:3) = (/ 0.2813255809899400E+00, 0.2283833222225700E-01, 0.6958360867878031E+00/)
       V_rule%lambda(26,1:3) = (/ 0.6958360867878030E+00, 0.2813255809899400E+00, 0.2283833222225706E-01/)
       V_rule%lambda(27,1:3) = (/ 0.2283833222225700E-01, 0.6958360867878030E+00, 0.2813255809899401E+00/)
       V_rule%lambda(28,1:3) = (/ 0.2573405054833000E-01, 0.1162519159075970E+00, 0.8580140335440730E+00/)
       V_rule%lambda(29,1:3) = (/ 0.1162519159075970E+00, 0.8580140335440730E+00, 0.2573405054833000E-01/)
       V_rule%lambda(30,1:3) = (/ 0.8580140335440730E+00, 0.2573405054833000E-01, 0.1162519159075970E+00/)
       V_rule%lambda(31,1:3) = (/ 0.1162519159075970E+00, 0.2573405054833000E-01, 0.8580140335440730E+00/)
       V_rule%lambda(32,1:3) = (/ 0.8580140335440730E+00, 0.1162519159075970E+00, 0.2573405054833000E-01/)
       V_rule%lambda(33,1:3) = (/ 0.2573405054833000E-01, 0.8580140335440730E+00, 0.1162519159075970E+00/)

      case(13)
       V_rule%Qdof = 37
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))
       V_rule%weights( 1) =  0.5252092340080200E-01
       V_rule%weights( 2: 4) =  0.1128014520933000E-01
       V_rule%weights( 5: 7) =  0.3142351836245400E-01
       V_rule%weights( 8:10) =  0.4707250250419400E-01
       V_rule%weights(11:13) =  0.4736358653635500E-01
       V_rule%weights(14:16) =  0.3116752904579400E-01
       V_rule%weights(17:19) =  0.7975771465074000E-02
       V_rule%weights(20:25) =  0.3684840272873200E-01
       V_rule%weights(26:31) =  0.1740146330382200E-01
       V_rule%weights(32:37) =  0.1552178683904500E-01

       V_rule%lambda( 1,1:3) = (/ 0.3333333333333330E+00, 0.3333333333333330E+00, 0.3333333333333340E+00/)
       V_rule%lambda( 2,1:3) = (/ 0.9903630120591001E-02, 0.4950481849397050E+00, 0.4950481849397040E+00/)
       V_rule%lambda( 3,1:3) = (/ 0.4950481849397050E+00, 0.4950481849397050E+00, 0.9903630120590035E-02/)
       V_rule%lambda( 4,1:3) = (/ 0.4950481849397050E+00, 0.9903630120591001E-02, 0.4950481849397040E+00/)
       V_rule%lambda( 5,1:3) = (/ 0.6256672978085200E-01, 0.4687166351095740E+00, 0.4687166351095740E+00/)
       V_rule%lambda( 6,1:3) = (/ 0.4687166351095740E+00, 0.4687166351095740E+00, 0.6256672978085198E-01/)
       V_rule%lambda( 7,1:3) = (/ 0.4687166351095740E+00, 0.6256672978085200E-01, 0.4687166351095740E+00/)
       V_rule%lambda( 8,1:3) = (/ 0.1709573263974470E+00, 0.4145213368012770E+00, 0.4145213368012760E+00/)
       V_rule%lambda( 9,1:3) = (/ 0.4145213368012770E+00, 0.4145213368012770E+00, 0.1709573263974460E+00/)
       V_rule%lambda(10,1:3) = (/ 0.4145213368012770E+00, 0.1709573263974470E+00, 0.4145213368012760E+00/)
       V_rule%lambda(11,1:3) = (/ 0.5412008559143370E+00, 0.2293995720428310E+00, 0.2293995720428320E+00/)
       V_rule%lambda(12,1:3) = (/ 0.2293995720428310E+00, 0.2293995720428310E+00, 0.5412008559143380E+00/)
       V_rule%lambda(13,1:3) = (/ 0.2293995720428310E+00, 0.5412008559143370E+00, 0.2293995720428320E+00/)
       V_rule%lambda(14,1:3) = (/ 0.7711510096073400E+00, 0.1144244951963300E+00, 0.1144244951963300E+00/)
       V_rule%lambda(15,1:3) = (/ 0.1144244951963300E+00, 0.1144244951963300E+00, 0.7711510096073400E+00/)
       V_rule%lambda(16,1:3) = (/ 0.1144244951963300E+00, 0.7711510096073400E+00, 0.1144244951963300E+00/)
       V_rule%lambda(17,1:3) = (/ 0.9503772172730820E+00, 0.2481139136345900E-01, 0.2481139136345900E-01/)
       V_rule%lambda(18,1:3) = (/ 0.2481139136345900E-01, 0.2481139136345900E-01, 0.9503772172730820E+00/)
       V_rule%lambda(19,1:3) = (/ 0.2481139136345900E-01, 0.9503772172730820E+00, 0.2481139136345900E-01/)
       V_rule%lambda(20,1:3) = (/ 0.9485382837957899E-01, 0.2687949970587610E+00, 0.6363511745616600E+00/)
       V_rule%lambda(21,1:3) = (/ 0.2687949970587610E+00, 0.6363511745616600E+00, 0.9485382837957901E-01/)
       V_rule%lambda(22,1:3) = (/ 0.6363511745616600E+00, 0.9485382837957899E-01, 0.2687949970587610E+00/)
       V_rule%lambda(23,1:3) = (/ 0.2687949970587610E+00, 0.9485382837957899E-01, 0.6363511745616600E+00/)
       V_rule%lambda(24,1:3) = (/ 0.6363511745616600E+00, 0.2687949970587610E+00, 0.9485382837957901E-01/)
       V_rule%lambda(25,1:3) = (/ 0.9485382837957899E-01, 0.6363511745616600E+00, 0.2687949970587610E+00/)
       V_rule%lambda(26,1:3) = (/ 0.1810077327880700E-01, 0.2917300667342880E+00, 0.6901691599869051E+00/)
       V_rule%lambda(27,1:3) = (/ 0.2917300667342880E+00, 0.6901691599869050E+00, 0.1810077327880699E-01/)
       V_rule%lambda(28,1:3) = (/ 0.6901691599869050E+00, 0.1810077327880700E-01, 0.2917300667342880E+00/)
       V_rule%lambda(29,1:3) = (/ 0.2917300667342880E+00, 0.1810077327880700E-01, 0.6901691599869051E+00/)
       V_rule%lambda(30,1:3) = (/ 0.6901691599869050E+00, 0.2917300667342880E+00, 0.1810077327880699E-01/)
       V_rule%lambda(31,1:3) = (/ 0.1810077327880700E-01, 0.6901691599869050E+00, 0.2917300667342880E+00/)
       V_rule%lambda(32,1:3) = (/ 0.2223307667409000E-01, 0.1263573854916690E+00, 0.8514095378342410E+00/)
       V_rule%lambda(33,1:3) = (/ 0.1263573854916690E+00, 0.8514095378342410E+00, 0.2223307667409002E-01/)
       V_rule%lambda(34,1:3) = (/ 0.8514095378342410E+00, 0.2223307667409000E-01, 0.1263573854916690E+00/)
       V_rule%lambda(35,1:3) = (/ 0.1263573854916690E+00, 0.2223307667409000E-01, 0.8514095378342410E+00/)
       V_rule%lambda(36,1:3) = (/ 0.8514095378342410E+00, 0.1263573854916690E+00, 0.2223307667409002E-01/)
       V_rule%lambda(37,1:3) = (/ 0.2223307667409000E-01, 0.8514095378342410E+00, 0.1263573854916690E+00/)

      case(14)
       V_rule%Qdof = 42
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))
       V_rule%weights( 1: 3) =  0.2188358136942900E-01
       V_rule%weights( 4: 6) =  0.3278835354412500E-01
       V_rule%weights( 7: 9) =  0.5177410450729200E-01
       V_rule%weights(10:12) =  0.4216258873699300E-01
       V_rule%weights(13:15) =  0.1443369966977700E-01
       V_rule%weights(16:18) =  0.4923403602400000E-02
       V_rule%weights(19:24) =  0.2466575321256400E-01
       V_rule%weights(25:30) =  0.3857151078706100E-01
       V_rule%weights(31:36) =  0.1443630811353400E-01
       V_rule%weights(37:42) =  0.5010228838501000E-02

       V_rule%lambda( 1,1:3) = (/ 0.2207217927564300E-01, 0.4889639103621790E+00, 0.4889639103621780E+00/)
       V_rule%lambda( 2,1:3) = (/ 0.4889639103621790E+00, 0.4889639103621790E+00, 0.2207217927564198E-01/)
       V_rule%lambda( 3,1:3) = (/ 0.4889639103621790E+00, 0.2207217927564300E-01, 0.4889639103621780E+00/)
       V_rule%lambda( 4,1:3) = (/ 0.1647105613190920E+00, 0.4176447193404540E+00, 0.4176447193404540E+00/)
       V_rule%lambda( 5,1:3) = (/ 0.4176447193404540E+00, 0.4176447193404540E+00, 0.1647105613190920E+00/)
       V_rule%lambda( 6,1:3) = (/ 0.4176447193404540E+00, 0.1647105613190920E+00, 0.4176447193404540E+00/)
       V_rule%lambda( 7,1:3) = (/ 0.4530449433823230E+00, 0.2734775283088390E+00, 0.2734775283088380E+00/)
       V_rule%lambda( 8,1:3) = (/ 0.2734775283088390E+00, 0.2734775283088390E+00, 0.4530449433823220E+00/)
       V_rule%lambda( 9,1:3) = (/ 0.2734775283088390E+00, 0.4530449433823230E+00, 0.2734775283088380E+00/)
       V_rule%lambda(10,1:3) = (/ 0.6455889351749130E+00, 0.1772055324125430E+00, 0.1772055324125440E+00/)
       V_rule%lambda(11,1:3) = (/ 0.1772055324125430E+00, 0.1772055324125430E+00, 0.6455889351749140E+00/)
       V_rule%lambda(12,1:3) = (/ 0.1772055324125430E+00, 0.6455889351749130E+00, 0.1772055324125440E+00/)
       V_rule%lambda(13,1:3) = (/ 0.8764002338182550E+00, 0.6179988309087300E-01, 0.6179988309087192E-01/)
       V_rule%lambda(14,1:3) = (/ 0.6179988309087300E-01, 0.6179988309087300E-01, 0.8764002338182539E+00/)
       V_rule%lambda(15,1:3) = (/ 0.6179988309087300E-01, 0.8764002338182550E+00, 0.6179988309087192E-01/)
       V_rule%lambda(16,1:3) = (/ 0.9612180775025980E+00, 0.1939096124870100E-01, 0.1939096124870099E-01/)
       V_rule%lambda(17,1:3) = (/ 0.1939096124870100E-01, 0.1939096124870100E-01, 0.9612180775025980E+00/)
       V_rule%lambda(18,1:3) = (/ 0.1939096124870100E-01, 0.9612180775025980E+00, 0.1939096124870099E-01/)
       V_rule%lambda(19,1:3) = (/ 0.5712475740364800E-01, 0.1722666878213560E+00, 0.7706085547749960E+00/)
       V_rule%lambda(20,1:3) = (/ 0.1722666878213560E+00, 0.7706085547749960E+00, 0.5712475740364797E-01/)
       V_rule%lambda(21,1:3) = (/ 0.7706085547749960E+00, 0.5712475740364800E-01, 0.1722666878213560E+00/)
       V_rule%lambda(22,1:3) = (/ 0.1722666878213560E+00, 0.5712475740364800E-01, 0.7706085547749960E+00/)
       V_rule%lambda(23,1:3) = (/ 0.7706085547749960E+00, 0.1722666878213560E+00, 0.5712475740364797E-01/)
       V_rule%lambda(24,1:3) = (/ 0.5712475740364800E-01, 0.7706085547749960E+00, 0.1722666878213560E+00/)
       V_rule%lambda(25,1:3) = (/ 0.9291624935697200E-01, 0.3368614597963450E+00, 0.5702222908466830E+00/)
       V_rule%lambda(26,1:3) = (/ 0.3368614597963450E+00, 0.5702222908466830E+00, 0.9291624935697196E-01/)
       V_rule%lambda(27,1:3) = (/ 0.5702222908466830E+00, 0.9291624935697200E-01, 0.3368614597963451E+00/)
       V_rule%lambda(28,1:3) = (/ 0.3368614597963450E+00, 0.9291624935697200E-01, 0.5702222908466830E+00/)
       V_rule%lambda(29,1:3) = (/ 0.5702222908466830E+00, 0.3368614597963450E+00, 0.9291624935697196E-01/)
       V_rule%lambda(30,1:3) = (/ 0.9291624935697200E-01, 0.5702222908466830E+00, 0.3368614597963451E+00/)
       V_rule%lambda(31,1:3) = (/ 0.1464695005565400E-01, 0.2983728821362580E+00, 0.6869801678080880E+00/)
       V_rule%lambda(32,1:3) = (/ 0.2983728821362580E+00, 0.6869801678080880E+00, 0.1464695005565397E-01/)
       V_rule%lambda(33,1:3) = (/ 0.6869801678080880E+00, 0.1464695005565400E-01, 0.2983728821362580E+00/)
       V_rule%lambda(34,1:3) = (/ 0.2983728821362580E+00, 0.1464695005565400E-01, 0.6869801678080880E+00/)
       V_rule%lambda(35,1:3) = (/ 0.6869801678080880E+00, 0.2983728821362580E+00, 0.1464695005565397E-01/)
       V_rule%lambda(36,1:3) = (/ 0.1464695005565400E-01, 0.6869801678080880E+00, 0.2983728821362580E+00/)
       V_rule%lambda(37,1:3) = (/ 0.1268330932872000E-02, 0.1189744976969570E+00, 0.8797571713701710E+00/)
       V_rule%lambda(38,1:3) = (/ 0.1189744976969570E+00, 0.8797571713701710E+00, 0.1268330932872042E-02/)
       V_rule%lambda(39,1:3) = (/ 0.8797571713701710E+00, 0.1268330932872000E-02, 0.1189744976969570E+00/)
       V_rule%lambda(40,1:3) = (/ 0.1189744976969570E+00, 0.1268330932872000E-02, 0.8797571713701710E+00/)
       V_rule%lambda(41,1:3) = (/ 0.8797571713701710E+00, 0.1189744976969570E+00, 0.1268330932872042E-02/)
       V_rule%lambda(42,1:3) = (/ 0.1268330932872000E-02, 0.8797571713701710E+00, 0.1189744976969570E+00/)

      case(15)
       V_rule%Qdof = 48
       allocate(V_rule%weights(1:V_rule%Qdof))
       allocate(V_rule%lambda(1:V_rule%Qdof,3))
       V_rule%weights( 1: 3) =  0.1916875642849000E-02
       V_rule%weights( 4: 6) =  0.4424902727114500E-01
       V_rule%weights( 7: 9) =  0.5118654871885200E-01
       V_rule%weights(10:12) =  0.2368773587068800E-01
       V_rule%weights(13:15) =  0.1328977569002100E-01
       V_rule%weights(16:18) =  0.4748916608192000E-02
       V_rule%weights(19:24) =  0.3855007259959300E-01
       V_rule%weights(25:30) =  0.2721581432062400E-01
       V_rule%weights(31:36) =  0.2182077366797000E-02
       V_rule%weights(37:42) =  0.2150531984773100E-01
       V_rule%weights(43:48) =  0.7673942631049000E-02

       V_rule%lambda( 1,1:3) = (/-0.1394583371648600E-01, 0.5069729168582430E+00, 0.5069729168582431E+00/)
       V_rule%lambda( 2,1:3) = (/ 0.5069729168582430E+00, 0.5069729168582430E+00,-0.1394583371648594E-01/)
       V_rule%lambda( 3,1:3) = (/ 0.5069729168582430E+00,-0.1394583371648600E-01, 0.5069729168582431E+00/)
       V_rule%lambda( 4,1:3) = (/ 0.1371872914339550E+00, 0.4314063542830230E+00, 0.4314063542830220E+00/)
       V_rule%lambda( 5,1:3) = (/ 0.4314063542830230E+00, 0.4314063542830230E+00, 0.1371872914339540E+00/)
       V_rule%lambda( 6,1:3) = (/ 0.4314063542830230E+00, 0.1371872914339550E+00, 0.4314063542830220E+00/)
       V_rule%lambda( 7,1:3) = (/ 0.4446127103057110E+00, 0.2776936448471440E+00, 0.2776936448471450E+00/)
       V_rule%lambda( 8,1:3) = (/ 0.2776936448471440E+00, 0.2776936448471440E+00, 0.4446127103057120E+00/)
       V_rule%lambda( 9,1:3) = (/ 0.2776936448471440E+00, 0.4446127103057110E+00, 0.2776936448471450E+00/)
       V_rule%lambda(10,1:3) = (/ 0.7470702179174920E+00, 0.1264648910412540E+00, 0.1264648910412540E+00/)
       V_rule%lambda(11,1:3) = (/ 0.1264648910412540E+00, 0.1264648910412540E+00, 0.7470702179174920E+00/)
       V_rule%lambda(12,1:3) = (/ 0.1264648910412540E+00, 0.7470702179174920E+00, 0.1264648910412540E+00/)
       V_rule%lambda(13,1:3) = (/ 0.8583832280506281E+00, 0.7080838597468600E-01, 0.7080838597468597E-01/)
       V_rule%lambda(14,1:3) = (/ 0.7080838597468600E-01, 0.7080838597468600E-01, 0.8583832280506281E+00/)
       V_rule%lambda(15,1:3) = (/ 0.7080838597468600E-01, 0.8583832280506281E+00, 0.7080838597468597E-01/)
       V_rule%lambda(16,1:3) = (/ 0.9620696595178529E+00, 0.1896517024107300E-01, 0.1896517024107403E-01/)
       V_rule%lambda(17,1:3) = (/ 0.1896517024107300E-01, 0.1896517024107300E-01, 0.9620696595178539E+00/)
       V_rule%lambda(18,1:3) = (/ 0.1896517024107300E-01, 0.9620696595178529E+00, 0.1896517024107403E-01/)
       V_rule%lambda(19,1:3) = (/ 0.1337341619666210E+00, 0.2613113711400870E+00, 0.6049544668932920E+00/)
       V_rule%lambda(20,1:3) = (/ 0.2613113711400870E+00, 0.6049544668932910E+00, 0.1337341619666219E+00/)
       V_rule%lambda(21,1:3) = (/ 0.6049544668932910E+00, 0.1337341619666210E+00, 0.2613113711400880E+00/)
       V_rule%lambda(22,1:3) = (/ 0.2613113711400870E+00, 0.1337341619666210E+00, 0.6049544668932920E+00/)
       V_rule%lambda(23,1:3) = (/ 0.6049544668932910E+00, 0.2613113711400870E+00, 0.1337341619666219E+00/)
       V_rule%lambda(24,1:3) = (/ 0.1337341619666210E+00, 0.6049544668932910E+00, 0.2613113711400880E+00/)
       V_rule%lambda(25,1:3) = (/ 0.3636667739691700E-01, 0.3880467670902690E+00, 0.5755865555128140E+00/)
       V_rule%lambda(26,1:3) = (/ 0.3880467670902690E+00, 0.5755865555128140E+00, 0.3636667739691701E-01/)
       V_rule%lambda(27,1:3) = (/ 0.5755865555128140E+00, 0.3636667739691700E-01, 0.3880467670902690E+00/)
       V_rule%lambda(28,1:3) = (/ 0.3880467670902690E+00, 0.3636667739691700E-01, 0.5755865555128140E+00/)
       V_rule%lambda(29,1:3) = (/ 0.5755865555128140E+00, 0.3880467670902690E+00, 0.3636667739691701E-01/)
       V_rule%lambda(30,1:3) = (/ 0.3636667739691700E-01, 0.5755865555128140E+00, 0.3880467670902690E+00/)
       V_rule%lambda(31,1:3) = (/-0.1017488312657100E-01, 0.2857122200499160E+00, 0.7244626630766551E+00/)
       V_rule%lambda(32,1:3) = (/ 0.2857122200499160E+00, 0.7244626630766550E+00,-0.1017488312657089E-01/)
       V_rule%lambda(33,1:3) = (/ 0.7244626630766550E+00,-0.1017488312657100E-01, 0.2857122200499160E+00/)
       V_rule%lambda(34,1:3) = (/ 0.2857122200499160E+00,-0.1017488312657100E-01, 0.7244626630766551E+00/)
       V_rule%lambda(35,1:3) = (/ 0.7244626630766550E+00, 0.2857122200499160E+00,-0.1017488312657089E-01/)
       V_rule%lambda(36,1:3) = (/-0.1017488312657100E-01, 0.7244626630766550E+00, 0.2857122200499160E+00/)
       V_rule%lambda(37,1:3) = (/ 0.3684386987587800E-01, 0.2155996640722840E+00, 0.7475564660518380E+00/)
       V_rule%lambda(38,1:3) = (/ 0.2155996640722840E+00, 0.7475564660518380E+00, 0.3684386987587795E-01/)
       V_rule%lambda(39,1:3) = (/ 0.7475564660518380E+00, 0.3684386987587800E-01, 0.2155996640722840E+00/)
       V_rule%lambda(40,1:3) = (/ 0.2155996640722840E+00, 0.3684386987587800E-01, 0.7475564660518380E+00/)
       V_rule%lambda(41,1:3) = (/ 0.7475564660518380E+00, 0.2155996640722840E+00, 0.3684386987587795E-01/)
       V_rule%lambda(42,1:3) = (/ 0.3684386987587800E-01, 0.7475564660518380E+00, 0.2155996640722840E+00/)
       V_rule%lambda(43,1:3) = (/ 0.1245980933119900E-01, 0.1035756165763860E+00, 0.8839645740924150E+00/)
       V_rule%lambda(44,1:3) = (/ 0.1035756165763860E+00, 0.8839645740924160E+00, 0.1245980933119795E-01/)
       V_rule%lambda(45,1:3) = (/ 0.8839645740924160E+00, 0.1245980933119900E-01, 0.1035756165763850E+00/)
       V_rule%lambda(46,1:3) = (/ 0.1035756165763860E+00, 0.1245980933119900E-01, 0.8839645740924150E+00/)
       V_rule%lambda(47,1:3) = (/ 0.8839645740924160E+00, 0.1035756165763860E+00, 0.1245980933119795E-01/)
       V_rule%lambda(48,1:3) = (/ 0.1245980933119900E-01, 0.8839645740924160E+00, 0.1035756165763850E+00/)

      end select

      if(Qdeg > 15) then
         print*,'Warning !! Volume quadrature rule of degree ', Qdeg, &
              ' is not implemented'
         print*,'           the maximal default quadrature rule of degree', &
              ' 15 is used'
      endif
    end subroutine Create_volume_rule_OLD


  !> generates one Langangian quadrature rule
  subroutine Create_L_rule(L_rule, Qnum)
    type(Lagrang_rule), intent(inout) :: L_rule
    integer, intent(in) :: Qnum
    real, dimension(:,:), allocatable :: M
    real, dimension(:), allocatable :: b
    integer :: i,j, k, l, n, Qdiff, Qdeg, Qdof

    Qdiff = 0
    if(Qnum < QnumOffset) then     !triangle
       Qdeg = Qnum
       L_rule%Qdeg = Qdeg
       L_rule%Qdof = (Qnum+1)*(Qnum+2)/2
    else
       Qdeg = Qnum - QnumOffset
       L_rule%Qdeg = Qdeg
       L_rule%Qdof = (Qnum-QnumOffset +1)*(Qnum-QnumOffset +2)
       Qdiff = (Qnum-QnumOffset +1)*(Qnum-QnumOffset +2)/2
    endif


    allocate(L_rule%lambda(1:L_rule%Qdof,1:2))

    select case (L_rule%Qdeg)
    case(:0)
       L_rule%lambda(1,1:2) = (/ 1./3, 1./3 /)
       if(Qdiff > 0) L_rule%lambda(2,1:2) = (/ 2./3, 2./3 /)

    case(1:)
       k = 0
       do i=0, L_rule%Qdeg
          do j=0, L_rule%Qdeg - i
             k = k + 1
             L_rule%lambda(k,1:2) = (/ 1.*j/Qdeg, 1.*i/Qdeg /)
             ! additional for quadrilaterals
             if(Qdiff > 0) &
                  L_rule%lambda(Qdiff+k,1:2) = (/ 1.*(Qdeg-j)/Qdeg, 1.*(Qdeg-i)/Qdeg /)
          enddo
       enddo
    end select

    ! only for triangles
    if(Qdiff == 0) then     !triangle
       Qdof = L_rule%Qdof
       allocate(L_rule%psi(1:Qdof,1:Qdof))
       allocate(L_rule%psi_pow(1:Qdof,1:2))

       ! settings of powers of x_1 and x_2 in Lagr. basis
       k = 0
       do i=0, Qdeg
          do j=0, i
             k = k + 1
             L_rule%psi_pow(k, 1) = i - j
             L_rule%psi_pow(k, 2) = j
          enddo
       enddo
       !write(*,'(a6,20i5)') '@@eds@',k,    L_rule%psi_pow(:, 1)
       !write(*,'(a6,20i5)') '@@wqa@',Qdof, L_rule%psi_pow(:, 2)
       !print*,'-------------------------'

       allocate( M(1:Qdof, 1:Qdof),  b(1:Qdof) )

       ! matrix with monopolynomial terms x^i y^j at the Lagrange node x_l
       do l=1,Qdof   ! index of Lagrange node
          do n=1, Qdof ! index of each monopolynomial terms x^i y^j
             i = L_rule%psi_pow(n, 1)
             j = L_rule%psi_pow(n, 2)

             M(l, n) = L_rule%lambda(l,1)**i * L_rule%lambda(l,2)**j

          enddo
       enddo



       ! inversion of M
       call MblockInverse(Qdof, M(1:Qdof, 1:Qdof) )

       ! evaluation of basis coeficients of Lagr. polynom
       do k=1, Qdof ! index of test function
          b(:) = 0.
          b(k) = 1.   ! RHS, psi_k (x_m) = \delta_{k m}

          L_rule%psi(k, 1:Qdof) = matmul( M(1:Qdof, 1:Qdof), b(1:Qdof) )

          !write(*,'(a6,i5, 20es12.4)') 'M',k,  M(k, 1:Qdof)
          !write(*,'(a6,i5, 200es12.4)') 'psi',k,  L_rule%psi(k, 1:Qdof)
       enddo
       !print*,'****************************************'

       deallocate(M, b)
    endif
    !if(Qdeg == 2) stop
  end subroutine Create_L_rule

  !> generates one Langangian quadrature rule
  subroutine Create_L_rule1(L_rule, Qnum)
    type(Lagrang_rule), intent(inout) :: L_rule
    integer, intent(in) :: Qnum

    if(Qnum < QnumOffset) then     !triangle
       L_rule%Qdeg = Qnum
       L_rule%Qdof = (Qnum+1)*(Qnum+2)/2
    else
       L_rule%Qdeg = Qnum - QnumOffset
       L_rule%Qdof = (Qnum-QnumOffset +1)*(Qnum-QnumOffset +2)
    endif

    allocate(L_rule%lambda(1:L_rule%Qdof,1:2))

    select case (L_rule%Qdeg)
    case(:0)
       L_rule%lambda(1,1:2) = (/ 1./3, 1./3 /)
    case(1)
       L_rule%lambda(1,1:2) = (/ 0., 0. /)
       L_rule%lambda(2,1:2) = (/ 1., 0. /)
       L_rule%lambda(3,1:2) = (/ 0., 1. /)
    case(2)
       L_rule%lambda(1,1:2) = (/ 0., 0. /)
       L_rule%lambda(2,1:2) = (/ 1., 0. /)
       L_rule%lambda(3,1:2) = (/ 0., 1. /)
       L_rule%lambda(4,1:2) = (/ 1./2, 0./2 /)
       L_rule%lambda(5,1:2) = (/ 1./2, 1./2 /)
       L_rule%lambda(6,1:2) = (/ 0./2, 1./2 /)
    case(3)
       L_rule%lambda(1,1:2) = (/ 0., 0. /)
       L_rule%lambda(2,1:2) = (/ 1., 0. /)
       L_rule%lambda(3,1:2) = (/ 0., 1. /)
       L_rule%lambda(4,1:2) = (/ 1./3, 0./3 /)
       L_rule%lambda(5,1:2) = (/ 2./3, 0./3 /)
       L_rule%lambda(6,1:2) = (/ 2./3, 1./3 /)
       L_rule%lambda(7,1:2) = (/ 1./3, 2./3 /)
       L_rule%lambda(8,1:2) = (/ 0./3, 2./3 /)
       L_rule%lambda(9,1:2) = (/ 0./3, 1./3 /)
       L_rule%lambda(10,1:2)= (/ 1./3, 1./3 /)
    end select

    if(L_rule%Qdeg > 3) then
       print*,'Warning !! Volume quadrature rule of degree ', L_rule%Qdeg, &
            ' is not implemented'
    endif


    if(Qnum >= QnumOffset) then  ! aditional nodes for  quadrilaterals
       select case (L_rule%Qdeg)
       case(:0)
          L_rule%lambda(2,1:2) = (/ 2./3, 2./3 /)
       case(1)
          L_rule%lambda(4,1:2) = (/ 1., 0. /)
          L_rule%lambda(5,1:2) = (/ 1., 1. /)
          L_rule%lambda(6,1:2) = (/ 0., 1. /)
       case(2)
          L_rule%lambda( 7,1:2) = (/ 1., 0. /)
          L_rule%lambda( 8,1:2) = (/ 1., 1. /)
          L_rule%lambda( 9,1:2) = (/ 0., 1. /)
          L_rule%lambda(10,1:2) = (/ 1., 1./2 /)
          L_rule%lambda(11,1:2) = (/ 1./2, 1. /)
          L_rule%lambda(12,1:2) = (/ 1./2, 1./2 /)
       case(3)
          L_rule%lambda(11,1:2) = (/ 1., 0. /)
          L_rule%lambda(12,1:2) = (/ 1., 1. /)
          L_rule%lambda(13,1:2) = (/ 0., 1. /)
          L_rule%lambda(14,1:2) = (/ 1. , 1./3 /)
          L_rule%lambda(15,1:2) = (/ 1. , 2./3 /)
          L_rule%lambda(16,1:2) = (/ 2./3, 1.  /)
          L_rule%lambda(17,1:2) = (/ 1./3, 1. /)
          L_rule%lambda(18,1:2) = (/ 1./3, 2./3 /)
          L_rule%lambda(19,1:2) = (/ 2./3, 1./3 /)
          L_rule%lambda(20,1:2)=  (/ 2./3, 2./3 /)
       end select

    endif

  end subroutine Create_L_rule1

  !> evaluate the Langrangian test function in a node given by
  !> barycentric coordinates
  subroutine Eval_L_rule1(L_rule, num, lambda, L_psi)
    type(Lagrang_rule), intent(inout) :: L_rule
    integer,  intent(in) :: num
    real, dimension(1:num,1:3), intent(in) :: lambda
    real, dimension(1:L_rule%Qdof, 1:num ), intent(out) :: L_psi

    select case (L_rule%Qdeg)
       case(:1)
          L_psi(1, 1:num) = lambda(1:num,3)
          L_psi(2, 1:num) = lambda(1:num,1)
          L_psi(3, 1:num) = lambda(1:num,2)
       case(2)
          L_psi(1, 1:num) = lambda(1:num,3) * (2*lambda(1:num,3) -1)
          L_psi(2, 1:num) = lambda(1:num,1) * (2*lambda(1:num,1) -1)
          L_psi(3, 1:num) = lambda(1:num,2) * (2*lambda(1:num,2) -1)
          L_psi(4, 1:num) = 4 * lambda(1:num,3) * lambda(1:num,1)
          L_psi(5, 1:num) = 4 * lambda(1:num,1) * lambda(1:num,2)
          L_psi(6, 1:num) = 4 * lambda(1:num,2) * lambda(1:num,3)
       case(3)
          L_psi(1, 1:num) = lambda(1:num,3)*(3*lambda(1:num,3)-1)*(3*lambda(1:num,3)-2)/2
          L_psi(2, 1:num) = lambda(1:num,1)*(3*lambda(1:num,1)-1)*(3*lambda(1:num,1)-2)/2
          L_psi(3, 1:num) = lambda(1:num,2)*(3*lambda(1:num,2)-1)*(3*lambda(1:num,2)-2)/2
          L_psi(4, 1:num) = 4.5 * lambda(1:num,3)*lambda(1:num,1)* (3*lambda(1:num,3)-1)
          L_psi(5, 1:num) = 4.5 * lambda(1:num,3)*lambda(1:num,1)* (3*lambda(1:num,1)-1)
          L_psi(6, 1:num) = 4.5 * lambda(1:num,1)*lambda(1:num,2)* (3*lambda(1:num,1)-1)
          L_psi(7, 1:num) = 4.5 * lambda(1:num,1)*lambda(1:num,2)* (3*lambda(1:num,2)-1)
          L_psi(8, 1:num) = 4.5 * lambda(1:num,2)*lambda(1:num,3)* (3*lambda(1:num,2)-1)
          L_psi(9, 1:num) = 4.5 * lambda(1:num,2)*lambda(1:num,3)* (3*lambda(1:num,3)-1)
          L_psi(10, 1:num) = 27 *  lambda(1:num,1)*lambda(1:num,2)* lambda(1:num,3)
       case(4)

          L_psi(1, 1:num)= lambda(1:num,1)*(4*lambda(1:num,1)-1) &
               *(2*lambda(1:num,1)-1)*(4*lambda(1:num,1)-3)/3
          L_psi(2, 1:num)= lambda(1:num,2)*(4*lambda(1:num,2)-1)&
               *(2*lambda(1:num,2)-1)*(4*lambda(1:num,2)-3)/3
          L_psi(3, 1:num)= lambda(1:num,3)*(4*lambda(1:num,3)-1)&
               *(2*lambda(1:num,3)-1)*(4*lambda(1:num,3)-3)/3

          L_psi(4, 1:num)= 16*lambda(1:num,1)*lambda(1:num,2)&
               *(4*lambda(1:num,1)-1)*(2*lambda(1:num,1)-1) /3
          L_psi(5, 1:num)= 4* lambda(1:num,1)*lambda(1:num,2)&
               *(4*lambda(1:num,1)-1)*(4*lambda(1:num,2)-1)
          L_psi(6, 1:num)= 16*lambda(1:num,1)*lambda(1:num,2)&
               *(4*lambda(1:num,2)-1)*(2*lambda(1:num,2)-1) /3

          L_psi(7, 1:num)= 16*lambda(1:num,2)*lambda(1:num,3)&
               *(4*lambda(1:num,2)-1)*(2*lambda(1:num,2)-1) /3
          L_psi(8, 1:num)= 4* lambda(1:num,2)*lambda(1:num,3)&
               *(4*lambda(1:num,2)-1)*(4*lambda(1:num,3)-1)
          L_psi(9, 1:num)= 16*lambda(1:num,2)*lambda(1:num,3)&
               *(4*lambda(1:num,3)-1)*(2*lambda(1:num,3)-1) /3

          L_psi(10, 1:num)= 16*lambda(1:num,3)*lambda(1:num,1)&
               *(4*lambda(1:num,3)-1)*(2*lambda(1:num,3)-1) /3
          L_psi(11, 1:num)= 4* lambda(1:num,3)*lambda(1:num,1)&
               *(4*lambda(1:num,3)-1)*(4*lambda(1:num,1)-1)
          L_psi(12, 1:num)= 16*lambda(1:num,3)*lambda(1:num,1)&
               *(4*lambda(1:num,1)-1)*(2*lambda(1:num,1)-1) /3

          L_psi(13, 1:num)= 32*lambda(1:num,1)*lambda(1:num,2)&
               *lambda(1:num,3)*(4*lambda(1:num,1)-1)
          L_psi(14, 1:num)= 32*lambda(1:num,1)*lambda(1:num,2)&
               *lambda(1:num,3)*(4*lambda(1:num,2)-1)
          L_psi(15, 1:num)= 32*lambda(1:num,1)*lambda(1:num,2)&
               *lambda(1:num,3)*(4*lambda(1:num,3)-1)

       case(5:)
          print*,'NOT IMPLEMENTED in Eval_L_rule1'
          stop
       end select
     end subroutine Eval_L_rule1


  !> evaluate the Langrangian test function  in a node given by
  !> barycentric coordinates directly from L_rule%psi
  subroutine Eval_L_Direct(L_rule, num, lambda, L_psi)
    type(Lagrang_rule), intent(inout) :: L_rule
    integer,  intent(in) :: num
    real, dimension(1:num,1:2), intent(in) :: lambda
    real, dimension(1:L_rule%Qdof, 1:num ), intent(out) :: L_psi
    integer :: Qdof, i, j, k, l

    Qdof = L_rule%Qdof

    do k=1, Qdof  ! index of the Lagr. basis function
       L_psi(k, 1:num) = 0.

       do l=1,Qdof ! index of each particular term x^i y^j
          i = L_rule%psi_pow(l,1)
          j = L_rule%psi_pow(l,2)

          L_psi(k, 1:num)  = L_psi(k, 1:num) &
               + L_rule%psi(k, l) * lambda(1:num, 1)**i * lambda(1:num, 2)**j
          !write(*,'(a6, 5i5, 30es12.4)') 'L_dir', L_rule%Qdeg, k,l,i,j, L_psi(k, 1:num)
       enddo ! l

    enddo ! k

   ! do k=1,Qdof
   !    write(*,'(a12, 3i5,300es12.4)') 'L_rule%psi',k,L_rule%psi_pow(k,1:2), L_rule%psi(k, :)
   ! enddo


  end subroutine Eval_L_Direct


  !> evaluate the Langrangian test function in a node given by
  !> barycentric coordinates
  subroutine Eval_L_rule(L_rule, num, lambda, L_psi)
    type(Lagrang_rule), intent(inout) :: L_rule
    integer,  intent(in) :: num
    real, dimension(1:num,1:3), intent(in) :: lambda
    real, dimension(1:L_rule%Qdof, 1:num ), intent(out) :: L_psi

    select case (L_rule%Qdeg)
    case(:1)
       L_psi(1, 1:num) = lambda(1:num,3)
       L_psi(2, 1:num) = lambda(1:num,1)
       L_psi(3, 1:num) = lambda(1:num,2)
    case(2)
       L_psi(1, 1:num) = lambda(1:num,3) * (2*lambda(1:num,3) -1)
       L_psi(3, 1:num) = lambda(1:num,1) * (2*lambda(1:num,1) -1)
       L_psi(6, 1:num) = lambda(1:num,2) * (2*lambda(1:num,2) -1)
       L_psi(2, 1:num) = 4 * lambda(1:num,3) * lambda(1:num,1)
       L_psi(5, 1:num) = 4 * lambda(1:num,1) * lambda(1:num,2)
       L_psi(4, 1:num) = 4 * lambda(1:num,2) * lambda(1:num,3)
    case(3)
       L_psi(1, 1:num) = lambda(1:num,3)*(3*lambda(1:num,3)-1)*(3*lambda(1:num,3)-2)/2
       L_psi(4, 1:num) = lambda(1:num,1)*(3*lambda(1:num,1)-1)*(3*lambda(1:num,1)-2)/2
       L_psi(10, 1:num) = lambda(1:num,2)*(3*lambda(1:num,2)-1)*(3*lambda(1:num,2)-2)/2
       L_psi(2, 1:num) = 4.5 * lambda(1:num,3)*lambda(1:num,1)* (3*lambda(1:num,3)-1)
       L_psi(3, 1:num) = 4.5 * lambda(1:num,3)*lambda(1:num,1)* (3*lambda(1:num,1)-1)
       L_psi(7, 1:num) = 4.5 * lambda(1:num,1)*lambda(1:num,2)* (3*lambda(1:num,1)-1)
       L_psi(9, 1:num) = 4.5 * lambda(1:num,1)*lambda(1:num,2)* (3*lambda(1:num,2)-1)
       L_psi(8, 1:num) = 4.5 * lambda(1:num,2)*lambda(1:num,3)* (3*lambda(1:num,2)-1)
       L_psi(5, 1:num) = 4.5 * lambda(1:num,2)*lambda(1:num,3)* (3*lambda(1:num,3)-1)
       L_psi(6, 1:num) = 27 *  lambda(1:num,1)*lambda(1:num,2)* lambda(1:num,3)
    case(4)

       L_psi(1, 1:num)= lambda(1:num,3)*(4*lambda(1:num,3)-1) &
            *(2*lambda(1:num,3)-1)*(4*lambda(1:num,3)-3)/3
       L_psi(5, 1:num)= lambda(1:num,1)*(4*lambda(1:num,1)-1)&
            *(2*lambda(1:num,1)-1)*(4*lambda(1:num,1)-3)/3
       L_psi(15, 1:num)= lambda(1:num,2)*(4*lambda(1:num,2)-1)&
            *(2*lambda(1:num,2)-1)*(4*lambda(1:num,2)-3)/3

       L_psi(2, 1:num)= 16*lambda(1:num,3)*lambda(1:num,1)&
            *(4*lambda(1:num,3)-1)*(2*lambda(1:num,3)-1) /3
       L_psi(3, 1:num)= 4* lambda(1:num,3)*lambda(1:num,1)&
            *(4*lambda(1:num,3)-1)*(4*lambda(1:num,1)-1)
       L_psi(4, 1:num)= 16*lambda(1:num,3)*lambda(1:num,1)&
            *(4*lambda(1:num,1)-1)*(2*lambda(1:num,1)-1) /3

       L_psi(9, 1:num)= 16*lambda(1:num,1)*lambda(1:num,2)&
            *(4*lambda(1:num,1)-1)*(2*lambda(1:num,1)-1) /3
       L_psi(12, 1:num)= 4* lambda(1:num,1)*lambda(1:num,2)&
            *(4*lambda(1:num,1)-1)*(4*lambda(1:num,2)-1)
       L_psi(14, 1:num)= 16*lambda(1:num,1)*lambda(1:num,2)&
            *(4*lambda(1:num,2)-1)*(2*lambda(1:num,2)-1) /3

       L_psi(13, 1:num)= 16*lambda(1:num,2)*lambda(1:num,3)&
            *(4*lambda(1:num,2)-1)*(2*lambda(1:num,2)-1) /3
       L_psi(10, 1:num)= 4* lambda(1:num,2)*lambda(1:num,3)&
            *(4*lambda(1:num,2)-1)*(4*lambda(1:num,3)-1)
       L_psi( 6, 1:num)= 16*lambda(1:num,2)*lambda(1:num,3)&
            *(4*lambda(1:num,3)-1)*(2*lambda(1:num,3)-1) /3

       L_psi( 7, 1:num)= 32*lambda(1:num,3)*lambda(1:num,1)&
            *lambda(1:num,2)*(4*lambda(1:num,3)-1)
       L_psi( 8, 1:num)= 32*lambda(1:num,3)*lambda(1:num,1)&
            *lambda(1:num,2)*(4*lambda(1:num,1)-1)
       L_psi(11, 1:num)= 32*lambda(1:num,3)*lambda(1:num,1)&
            *lambda(1:num,2)*(4*lambda(1:num,2)-1)

    case(5:)
       print*,'NOT IMPLEMENTED in Eval_L_rule'
       stop
    end select
  end subroutine Eval_L_rule


  !> generates one volume (quadrilateral "bi-Gauss" ) quadrature rule
  subroutine Add_4volume_rule(Qdeg, G_rule, V_rule)
    type(Gauss_rule), intent(in) :: G_rule
    type(volume_rule), intent(inout) :: V_rule
    integer, intent(in) :: Qdeg
    integer :: Qdof, i, j

    V_rule%def = .true.

    Qdof = Qdeg*Qdeg
    V_rule%Qdeg = Qdeg
    V_rule%Qdof = Qdof

    allocate(V_rule%weights(1:V_rule%Qdof))
    allocate(V_rule%lambda(1:V_rule%Qdof,2))

    do i=1,Qdeg
       do j=1,Qdeg
          V_rule%weights((i-1)*Qdeg + j) = G_rule%weights(i)*G_rule%weights(j)
          V_rule%lambda((i-1)*Qdeg + j,1) = G_rule%lambda(i)
          V_rule%lambda((i-1)*Qdeg + j,2) = G_rule%lambda(j)
       enddo
    enddo
  end subroutine Add_4volume_rule


!  !> generates integ rule for time intervals
!  subroutine Create_time_rule(T_rule, Qdeg)
!    type(time_rule), intent(inout) :: T_rule
!    integer, intent(in) :: Qdeg
!
!    T_rule%Qdeg = max(Qdeg,1)
!!    T_rule%Qdof = T_rule%Qdeg
!
!
!    allocate(T_rule%weights(1:T_rule%Qdeg))
!    allocate(T_rule%lambda(1:T_rule%Qdeg))
!
!    select case (Qdeg)
!    case(:1)
!       T_rule%weights(1) = 1.
!       T_rule%lambda(1) = 0.5
!    case(2)
!       T_rule%weights(1:2) = 0.5
!       T_rule%lambda(1) = (1.-0.57735026918962576451)/2
!       T_rule%lambda(2) = 1.57735026918962576451/2
!    case(3)
!       T_rule%weights(1) = 5./18
!       T_rule%weights(2) = 8./18
!       T_rule%weights(3) = 5./18
!       T_rule%lambda(1) = (1.-0.77459666924148337704)/2
!       T_rule%lambda(2) = 0.5
!       T_rule%lambda(3) = 1.77459666924148337704/2
!    case(4)
!       T_rule%weights(1) = 0.34785484513745385737/2
!       T_rule%weights(2) = 0.65214515486254614263/2
!       T_rule%weights(3) = T_rule%weights(2)
!       T_rule%weights(4) = T_rule%weights(1)
!       T_rule%lambda(1) = (1-0.86113631159405257522)/2
!       T_rule%lambda(2) = (1-0.33998104358485626480)/2
!       T_rule%lambda(3) = 1.33998104358485626480/2
!       T_rule%lambda(4) = 1.86113631159405257522/2
!    case(5)
!       T_rule%weights(1) = 0.23692688505618908751/2.
!       T_rule%weights(2) = 0.47862867049936646804/2.
!       T_rule%weights(3) = 0.56888888888888888889/2.
!       T_rule%weights(4) = T_rule%weights(2)
!       T_rule%weights(5) = T_rule%weights(1)
!
!       T_rule%lambda(1) = (1-0.90617984593866399280)/2.
!       T_rule%lambda(2) = (1-0.53846931010568309104)/2.
!       T_rule%lambda(3) = 0.5
!       T_rule%lambda(4) = 1.53846931010568309104/2.
!       T_rule%lambda(5) = 1.90617984593866399280/2.
!    case(6)
!       T_rule%weights(1) = 0.17132449237917034504/2.
!       T_rule%weights(2) = 0.36076157304813860757/2.
!       T_rule%weights(3) = 0.46791393457269104739/2.
!       T_rule%weights(4) = T_rule%weights(3)
!       T_rule%weights(5) = T_rule%weights(2)
!       T_rule%weights(6) = T_rule%weights(1)
!
!       T_rule%lambda(1) = (1.-0.93246951420315202781)/2.
!       T_rule%lambda(2) = (1.-0.66120938646626451366)/2.
!       T_rule%lambda(3) = (1.-0.23861918608319690863)/2.
!       T_rule%lambda(4) =  1.23861918608319690863/2.
!       T_rule%lambda(5) =  1.66120938646626451366/2.
!       T_rule%lambda(6) =  1.93246951420315202781/2.
!     case(7)
!       T_rule%weights(1) = 0.12948496616886969327/2.
!       T_rule%weights(2) = 0.27970539148927666790/2.
!       T_rule%weights(3) = 0.38183005050511894495/2.
!       T_rule%weights(4) = 0.41795918367346938776/2.
!       T_rule%weights(5) = T_rule%weights(3)
!       T_rule%weights(6) = T_rule%weights(2)
!       T_rule%weights(7) = T_rule%weights(1)
!
!       T_rule%lambda(1) = (1.-0.94910791234275852453)/2.
!       T_rule%lambda(2) = (1.-0.74153118559939443986)/2.
!       T_rule%lambda(3) = (1.-0.40584515137739716691)/2.
!       T_rule%lambda(4) =  0.5
!       T_rule%lambda(5) =  1.40584515137739716691/2.
!       T_rule%lambda(6) =  1.74153118559939443986/2.
!       T_rule%lambda(7) =  1.94910791234275852453/2.
!
!      case( 8)
!       T_rule%weights( 1) =  0.5061426814518818E-01
!       T_rule%weights( 2) =  0.1111905172266872E+00
!       T_rule%weights( 3) =  0.1568533229389436E+00
!       T_rule%weights( 4) =  0.1813418916891809E+00
!       T_rule%weights( 5) =  0.1813418916891809E+00
!       T_rule%weights( 6) =  0.1568533229389436E+00
!       T_rule%weights( 7) =  0.1111905172266872E+00
!       T_rule%weights( 8) =  0.5061426814518818E-01
!
!       T_rule%lambda( 1) =  0.1985507175123191E-01
!       T_rule%lambda( 2) =  0.1016667612931866E+00
!       T_rule%lambda( 3) =  0.2372337950418355E+00
!       T_rule%lambda( 4) =  0.4082826787521751E+00
!       T_rule%lambda( 5) =  0.5917173212478249E+00
!       T_rule%lambda( 6) =  0.7627662049581645E+00
!       T_rule%lambda( 7) =  0.8983332387068134E+00
!       T_rule%lambda( 8) =  0.9801449282487681E+00
!
!      case( 9)
!       T_rule%weights( 1) =  0.4063719418078721E-01
!       T_rule%weights( 2) =  0.9032408034742863E-01
!       T_rule%weights( 3) =  0.1303053482014677E+00
!       T_rule%weights( 4) =  0.1561735385200015E+00
!       T_rule%weights( 5) =  0.1651196775006299E+00
!       T_rule%weights( 6) =  0.1561735385200015E+00
!       T_rule%weights( 7) =  0.1303053482014677E+00
!       T_rule%weights( 8) =  0.9032408034742863E-01
!       T_rule%weights( 9) =  0.4063719418078721E-01
!
!       T_rule%lambda( 1) =  0.1591988024618696E-01
!       T_rule%lambda( 2) =  0.8198444633668212E-01
!       T_rule%lambda( 3) =  0.1933142836497048E+00
!       T_rule%lambda( 4) =  0.3378732882980955E+00
!       T_rule%lambda( 5) =  0.5000000000000000E+00
!       T_rule%lambda( 6) =  0.6621267117019045E+00
!       T_rule%lambda( 7) =  0.8066857163502952E+00
!       T_rule%lambda( 8) =  0.9180155536633179E+00
!       T_rule%lambda( 9) =  0.9840801197538130E+00
!
!      case(10)
!       T_rule%weights( 1) =  0.3333567215434402E-01
!       T_rule%weights( 2) =  0.7472567457529031E-01
!       T_rule%weights( 3) =  0.1095431812579911E+00
!       T_rule%weights( 4) =  0.1346333596549982E+00
!       T_rule%weights( 5) =  0.1477621123573765E+00
!       T_rule%weights( 6) =  0.1477621123573765E+00
!       T_rule%weights( 7) =  0.1346333596549982E+00
!       T_rule%weights( 8) =  0.1095431812579911E+00
!       T_rule%weights( 9) =  0.7472567457529031E-01
!       T_rule%weights(10) =  0.3333567215434402E-01
!
!       T_rule%lambda( 1) =  0.1304673574141413E-01
!       T_rule%lambda( 2) =  0.6746831665550773E-01
!       T_rule%lambda( 3) =  0.1602952158504878E+00
!       T_rule%lambda( 4) =  0.2833023029353764E+00
!       T_rule%lambda( 5) =  0.4255628305091844E+00
!       T_rule%lambda( 6) =  0.5744371694908156E+00
!       T_rule%lambda( 7) =  0.7166976970646236E+00
!       T_rule%lambda( 8) =  0.8397047841495122E+00
!       T_rule%lambda( 9) =  0.9325316833444923E+00
!       T_rule%lambda(10) =  0.9869532642585859E+00
!
!      case(11)
!       T_rule%weights( 1) =  0.2783428355808687E-01
!       T_rule%weights( 2) =  0.6279018473245233E-01
!       T_rule%weights( 3) =  0.9314510546386705E-01
!       T_rule%weights( 4) =  0.1165968822959953E+00
!       T_rule%weights( 5) =  0.1314022722551234E+00
!       T_rule%weights( 6) =  0.1364625433889503E+00
!       T_rule%weights( 7) =  0.1314022722551234E+00
!       T_rule%weights( 8) =  0.1165968822959953E+00
!       T_rule%weights( 9) =  0.9314510546386705E-01
!       T_rule%weights(10) =  0.6279018473245233E-01
!       T_rule%weights(11) =  0.2783428355808687E-01
!
!       T_rule%lambda( 1) =  0.1088567092697151E-01
!       T_rule%lambda( 2) =  0.5646870011595234E-01
!       T_rule%lambda( 3) =  0.1349239972129753E+00
!       T_rule%lambda( 4) =  0.2404519353965941E+00
!       T_rule%lambda( 5) =  0.3652284220238275E+00
!       T_rule%lambda( 6) =  0.5000000000000000E+00
!       T_rule%lambda( 7) =  0.6347715779761725E+00
!       T_rule%lambda( 8) =  0.7595480646034058E+00
!       T_rule%lambda( 9) =  0.8650760027870247E+00
!       T_rule%lambda(10) =  0.9435312998840477E+00
!       T_rule%lambda(11) =  0.9891143290730284E+00
!
!      case(12)
!       T_rule%weights( 1) =  0.2358766819325595E-01
!       T_rule%weights( 2) =  0.5346966299765912E-01
!       T_rule%weights( 3) =  0.8003916427167312E-01
!       T_rule%weights( 4) =  0.1015837133615329E+00
!       T_rule%weights( 5) =  0.1167462682691774E+00
!       T_rule%weights( 6) =  0.1245735229067013E+00
!       T_rule%weights( 7) =  0.1245735229067013E+00
!       T_rule%weights( 8) =  0.1167462682691774E+00
!       T_rule%weights( 9) =  0.1015837133615329E+00
!       T_rule%weights(10) =  0.8003916427167312E-01
!       T_rule%weights(11) =  0.5346966299765912E-01
!       T_rule%weights(12) =  0.2358766819325595E-01
!
!       T_rule%lambda( 1) =  0.9219682876640378E-02
!       T_rule%lambda( 2) =  0.4794137181476255E-01
!       T_rule%lambda( 3) =  0.1150486629028477E+00
!       T_rule%lambda( 4) =  0.2063410228566913E+00
!       T_rule%lambda( 5) =  0.3160842505009099E+00
!       T_rule%lambda( 6) =  0.4373832957442655E+00
!       T_rule%lambda( 7) =  0.5626167042557344E+00
!       T_rule%lambda( 8) =  0.6839157494990901E+00
!       T_rule%lambda( 9) =  0.7936589771433087E+00
!       T_rule%lambda(10) =  0.8849513370971523E+00
!       T_rule%lambda(11) =  0.9520586281852375E+00
!       T_rule%lambda(12) =  0.9907803171233596E+00
!
!      case(13)
!       T_rule%weights( 1) =  0.2024200238265791E-01
!       T_rule%weights( 2) =  0.4606074991886418E-01
!       T_rule%weights( 3) =  0.6943675510989367E-01
!       T_rule%weights( 4) =  0.8907299038097291E-01
!       T_rule%weights( 5) =  0.1039080237684442E+00
!       T_rule%weights( 6) =  0.1131415901314486E+00
!       T_rule%weights( 7) =  0.1162757766154369E+00
!       T_rule%weights( 8) =  0.1131415901314486E+00
!       T_rule%weights( 9) =  0.1039080237684442E+00
!       T_rule%weights(10) =  0.8907299038097291E-01
!       T_rule%weights(11) =  0.6943675510989367E-01
!       T_rule%weights(12) =  0.4606074991886418E-01
!       T_rule%weights(13) =  0.2024200238265791E-01
!
!       T_rule%lambda( 1) =  0.7908472640705932E-02
!       T_rule%lambda( 2) =  0.4120080038851104E-01
!       T_rule%lambda( 3) =  0.9921095463334506E-01
!       T_rule%lambda( 4) =  0.1788253302798299E+00
!       T_rule%lambda( 5) =  0.2757536244817765E+00
!       T_rule%lambda( 6) =  0.3847708420224326E+00
!       T_rule%lambda( 7) =  0.5000000000000000E+00
!       T_rule%lambda( 8) =  0.6152291579775674E+00
!       T_rule%lambda( 9) =  0.7242463755182235E+00
!       T_rule%lambda(10) =  0.8211746697201701E+00
!       T_rule%lambda(11) =  0.9007890453666549E+00
!       T_rule%lambda(12) =  0.9587991996114890E+00
!       T_rule%lambda(13) =  0.9920915273592941E+00
!
!      case(14)
!       T_rule%weights( 1) =  0.1755973016587599E-01
!       T_rule%weights( 2) =  0.4007904357988007E-01
!       T_rule%weights( 3) =  0.6075928534395157E-01
!       T_rule%weights( 4) =  0.7860158357909681E-01
!       T_rule%weights( 5) =  0.9276919873896890E-01
!       T_rule%weights( 6) =  0.1025992318606479E+00
!       T_rule%weights( 7) =  0.1076319267315789E+00
!       T_rule%weights( 8) =  0.1076319267315789E+00
!       T_rule%weights( 9) =  0.1025992318606479E+00
!       T_rule%weights(10) =  0.9276919873896890E-01
!       T_rule%weights(11) =  0.7860158357909681E-01
!       T_rule%weights(12) =  0.6075928534395157E-01
!       T_rule%weights(13) =  0.4007904357988007E-01
!       T_rule%weights(14) =  0.1755973016587599E-01
!
!       T_rule%lambda( 1) =  0.6858095651593843E-02
!       T_rule%lambda( 2) =  0.3578255816821324E-01
!       T_rule%lambda( 3) =  0.8639934246511749E-01
!       T_rule%lambda( 4) =  0.1563535475941573E+00
!       T_rule%lambda( 5) =  0.2423756818209230E+00
!       T_rule%lambda( 6) =  0.3404438155360551E+00
!       T_rule%lambda( 7) =  0.4459725256463282E+00
!       T_rule%lambda( 8) =  0.5540274743536718E+00
!       T_rule%lambda( 9) =  0.6595561844639448E+00
!       T_rule%lambda(10) =  0.7576243181790770E+00
!       T_rule%lambda(11) =  0.8436464524058427E+00
!       T_rule%lambda(12) =  0.9136006575348825E+00
!       T_rule%lambda(13) =  0.9642174418317868E+00
!       T_rule%lambda(14) =  0.9931419043484062E+00
!
!     case(15:)
!       T_rule%weights( 1) =  0.1537662099805857E-01
!       T_rule%weights( 2) =  0.3518302374405406E-01
!       T_rule%weights( 3) =  0.5357961023358598E-01
!       T_rule%weights( 4) =  0.6978533896307712E-01
!       T_rule%weights( 5) =  0.8313460290849704E-01
!       T_rule%weights( 6) =  0.9308050000778109E-01
!       T_rule%weights( 7) =  0.9921574266355583E-01
!       T_rule%weights( 8) =  0.1012891209627806E+00
!       T_rule%weights( 9) =  0.9921574266355583E-01
!       T_rule%weights(10) =  0.9308050000778109E-01
!       T_rule%weights(11) =  0.8313460290849704E-01
!       T_rule%weights(12) =  0.6978533896307712E-01
!       T_rule%weights(13) =  0.5357961023358598E-01
!       T_rule%weights(14) =  0.3518302374405406E-01
!       T_rule%weights(15) =  0.1537662099805857E-01
!
!       T_rule%lambda( 1) =  0.6003740989757256E-02
!       T_rule%lambda( 2) =  0.3136330379964708E-01
!       T_rule%lambda( 3) =  0.7589670829478640E-01
!       T_rule%lambda( 4) =  0.1377911343199150E+00
!       T_rule%lambda( 5) =  0.2145139136957306E+00
!       T_rule%lambda( 6) =  0.3029243264612183E+00
!       T_rule%lambda( 7) =  0.3994029530012828E+00
!       T_rule%lambda( 8) =  0.5000000000000000E+00
!       T_rule%lambda( 9) =  0.6005970469987173E+00
!       T_rule%lambda(10) =  0.6970756735387817E+00
!       T_rule%lambda(11) =  0.7854860863042694E+00
!       T_rule%lambda(12) =  0.8622088656800850E+00
!       T_rule%lambda(13) =  0.9241032917052137E+00
!       T_rule%lambda(14) =  0.9686366962003530E+00
!       T_rule%lambda(15) =  0.9939962590102427E+00
!
!    end select
!
!    if(Qdeg > 15) then
!      print*,'Warning !! Time quadrature rule of degree ', Qdeg, &
!      ' is not implemented'
!      print*,'           the maximal default quadrature rule of degree', &
!      ' 15 is used'
!    endif
!
!
!
!
!    !select case (Qdeg)
!    !case( :1  )  ! AT THIS MOMENT ASSOCIATED WITH GAUSS, NOT NECESSARY TO DEFINE
!       !V_rule%Qdof =   1
!       !allocate(V_rule%weights(1:V_rule%Qdof))
!       !allocate(V_rule%lambda(1:V_rule%Qdof,3))
!
!    !end select
!  end subroutine Create_time_rule

!  !> generates one edge (Radau) quadrature rule
!  subroutine Create_Radau_rule(R_rule, Qdeg)
!    type(Time_rule), intent(inout) :: R_rule
!    integer, intent(in) :: Qdeg
!
!    R_rule%Qdeg = max(Qdeg, 1)
!!    R_rule%Qdof = R_rule%Qdeg
!
!    allocate(R_rule%weights(1:R_rule%Qdeg))
!    allocate(R_rule%lambda(1:R_rule%Qdeg))
!
!    select case (Qdeg)
!    case(:1)
!       R_rule%weights(1) = 1.
!       R_rule%lambda(1) = 1.
!    case(2)
!       R_rule%weights(1) =   0.7500000000000000
!       R_rule%weights(2) =   0.2500000000000000
!       R_rule%lambda(1) =   0.3333333333333333
!       R_rule%lambda(2) =   1.0000000000000000
!    case(3)
!       R_rule%weights(1) =   0.3764030627004669
!       R_rule%weights(2) =   0.5124858261884215
!       R_rule%weights(3) =   0.1111111111111111
!       R_rule%lambda(1) =   0.1550510257216822
!       R_rule%lambda(2) =   0.6449489742783179
!       R_rule%lambda(3) =   1.0000000000000000
!    case(4)
!       R_rule%weights(1) =   0.2204622111767679
!       R_rule%weights(2) =   0.3881934688431719
!       R_rule%weights(3) =   0.3288443199800598
!       R_rule%weights(4) =   0.0625000000000000
!       R_rule%lambda(1) =   0.0885879595127039
!       R_rule%lambda(2) =   0.4094668644407347
!       R_rule%lambda(3) =   0.7876594617608470
!       R_rule%lambda(4) =   1.0000000000000000
!    case(5)
!       R_rule%weights(1) =   0.1437135607912257
!       R_rule%weights(2) =   0.2813560151494621
!       R_rule%weights(3) =   0.3118265229757412
!       R_rule%weights(4) =   0.2231039010835709
!       R_rule%weights(5) =   0.0400000000000000
!       R_rule%lambda(1) =   0.0571041961145177
!       R_rule%lambda(2) =   0.2768430136381238
!       R_rule%lambda(3) =   0.5835904323689168
!       R_rule%lambda(4) =   0.8602401356562195
!       R_rule%lambda(5) =   1.0000000000000000
!    case(6)
!       R_rule%weights(1) =   0.1007941926267402
!       R_rule%weights(2) =   0.2084506671559546
!       R_rule%weights(3) =   0.2604633915947875
!       R_rule%weights(4) =   0.2426935942344849
!       R_rule%weights(5) =   0.1598203766102554
!       R_rule%weights(6) =   0.0277777777777778
!       R_rule%lambda(1) =   0.0398098570514687
!       R_rule%lambda(2) =   0.1980134178736082
!       R_rule%lambda(3) =   0.4379748102473862
!       R_rule%lambda(4) =   0.6954642733536361
!       R_rule%lambda(5) =   0.9014649142011736
!       R_rule%lambda(6) =   1.0000000000000000
!    case(7)
!       R_rule%weights(1) =   0.0744942355560100
!       R_rule%weights(2) =   0.1591021157336509
!       R_rule%weights(3) =   0.2123518895029778
!       R_rule%weights(4) =   0.2235549145072832
!       R_rule%weights(5) =   0.1904749368221156
!       R_rule%weights(6) =   0.1196137446126562
!       R_rule%weights(7) =   0.0204081632653061
!       R_rule%lambda(1) =   0.0293164271597849
!       R_rule%lambda(2) =   0.1480785996684843
!       R_rule%lambda(3) =   0.3369846902811543
!       R_rule%lambda(4) =   0.5586715187715501
!       R_rule%lambda(5) =   0.7692338620300545
!       R_rule%lambda(6) =   0.9269456713197410
!       R_rule%lambda(7) =   1.0000000000000000
!    case(8)
!       R_rule%weights(1) =   0.0572544073721283
!       R_rule%weights(2) =   0.1248239506649321
!       R_rule%weights(3) =   0.1735073978172506
!       R_rule%weights(4) =   0.1957860837262468
!       R_rule%weights(5) =   0.1882587726945593
!       R_rule%weights(6) =   0.1520653103233925
!       R_rule%weights(7) =   0.0926790774014897
!       R_rule%weights(8) =   0.0156250000000000
!       R_rule%lambda(1) =   0.0224793864387125
!       R_rule%lambda(2) =   0.1146790531609042
!       R_rule%lambda(3) =   0.2657898227845895
!       R_rule%lambda(4) =   0.4528463736694446
!       R_rule%lambda(5) =   0.6473752828868303
!       R_rule%lambda(6) =   0.8197593082631076
!       R_rule%lambda(7) =   0.9437374394630779
!       R_rule%lambda(8) =   1.0000000000000000
!    case(9)
!       R_rule%weights(1) =   0.0453572524616435
!       R_rule%weights(2) =   0.1002766490122754
!       R_rule%weights(3) =   0.1431933481786156
!       R_rule%weights(4) =   0.1688469834879647
!       R_rule%weights(5) =   0.1741365013864834
!       R_rule%weights(6) =   0.1584218878352190
!       R_rule%weights(7) =   0.1235946891022966
!       R_rule%weights(8) =   0.0738270095231576
!       R_rule%weights(9) =   0.0123456790123457
!       R_rule%lambda(1) =   0.0177799151473634
!       R_rule%lambda(2) =   0.0913236078997939
!       R_rule%lambda(3) =   0.2143084793956307
!       R_rule%lambda(4) =   0.3719321645832723
!       R_rule%lambda(5) =   0.5451866848034267
!       R_rule%lambda(6) =   0.7131752428555695
!       R_rule%lambda(7) =   0.8556337429578544
!       R_rule%lambda(8) =   0.9553660447100302
!       R_rule%lambda(9) =   1.0000000000000000
!    case(10)
!       R_rule%weights(1) =   0.0368085027433803
!       R_rule%weights(2) =   0.0821880063684609
!       R_rule%weights(3) =   0.1195967158571898
!       R_rule%weights(4) =   0.1453050824164592
!       R_rule%weights(5) =   0.1567912286134692
!       R_rule%weights(6) =   0.1529296438622114
!       R_rule%weights(7) =   0.1340974189205893
!       R_rule%weights(8) =   0.1021350659395004
!       R_rule%weights(9) =   0.0601483352787409
!       R_rule%weights(10) =   0.0100000000000000
!       R_rule%lambda(1) =   0.0144124096488766
!       R_rule%lambda(2) =   0.0743873897091961
!       R_rule%lambda(3) =   0.1761166561629953
!       R_rule%lambda(4) =   0.3096675799276378
!       R_rule%lambda(5) =   0.4619704010810110
!       R_rule%lambda(6) =   0.6181172346952940
!       R_rule%lambda(7) =   0.7628230151850396
!       R_rule%lambda(8) =   0.8819210212100013
!       R_rule%lambda(9) =   0.9637421871167906
!       R_rule%lambda(10) =   1.0000000000000000
!    case(11)
!       R_rule%weights(1) =   0.0304625489060658
!       R_rule%weights(2) =   0.0685168410666014
!       R_rule%weights(3) =   0.1010815542700119
!       R_rule%weights(4) =   0.1254626888485642
!       R_rule%weights(5) =   0.1396806665516916
!       R_rule%weights(6) =   0.1425827819705036
!       R_rule%weights(7) =   0.1339335430948421
!       R_rule%weights(8) =   0.1144330619244883
!       R_rule%weights(9) =   0.0856588096033299
!       R_rule%weights(10) =   0.0499230409539840
!       R_rule%weights(11) =   0.0082644628099174
!       R_rule%lambda(1) =   0.0119176134324156
!       R_rule%lambda(2) =   0.0617320718771481
!       R_rule%lambda(3) =   0.1471114496430702
!       R_rule%lambda(4) =   0.2611596760084562
!       R_rule%lambda(5) =   0.3946398468857868
!       R_rule%lambda(6) =   0.5367387657156606
!       R_rule%lambda(7) =   0.6759444616766651
!       R_rule%lambda(8) =   0.8009789210368988
!       R_rule%lambda(9) =   0.9017109877901468
!       R_rule%lambda(10) =   0.9699709678385136
!       R_rule%lambda(11) =   1.0000000000000000
!    case(12:)
!       R_rule%weights(1) =   0.0256240496036348
!       R_rule%weights(2) =   0.0579537401458690
!       R_rule%weights(3) =   0.0863853196566543
!       R_rule%weights(4) =   0.1089344395130960
!       R_rule%weights(5) =   0.1240607804020050
!       R_rule%weights(6) =   0.1307328302760666
!       R_rule%weights(7) =   0.1284956690763539
!       R_rule%weights(8) =   0.1175015575724929
!       R_rule%weights(9) =   0.0984992674130449
!       R_rule%weights(10) =   0.0727818344269976
!       R_rule%weights(11) =   0.0420860674693405
!       R_rule%weights(12) =   0.0069444444444444
!       R_rule%lambda(1) =   0.0100182804616804
!       R_rule%lambda(2) =   0.0520354511271806
!       R_rule%lambda(3) =   0.1246192251444431
!       R_rule%lambda(4) =   0.2228406070438378
!       R_rule%lambda(5) =   0.3400081579146652
!       R_rule%lambda(6) =   0.4681376130895840
!       R_rule%lambda(7) =   0.5984972797671392
!       R_rule%lambda(8) =   0.7222032848909679
!       R_rule%lambda(9) =   0.8308248996228186
!       R_rule%lambda(10) =   0.9169583865525948
!       R_rule%lambda(11) =   0.9747263796024797
!       R_rule%lambda(12) =   1.0000000000000000
!
!    end select
!
!    if(Qdeg > 12) then
!      print*,'Warning !! Radau rule of degree ', Qdeg, &
!      ' is not implemented'
!      print*,'           the maximal default quadrature rule of degree', &
!      ' 12 is used'
!    endif
!  end subroutine Create_Radau_rule

  !> allocate ABDF method
  subroutine AllocateABDF(BDF, max_deg)
    type(ABDF), intent(inout) :: BDF
    integer, intent(in) :: max_deg

    BDF%max_deg = max_deg
    allocate(BDF%alpha(0:max_deg+1) )
    allocate(BDF%delta(0:max_deg+1) )   ! for estimation of the local error

    allocate(BDF%extrap(1:max_deg) )
    !allocate(BDF%Bextrap(1:max_deg+1) )
  end subroutine AllocateABDF


!  !> set coefficients \f$ \alpha_k,\ \beta_k,\ \dots \f$ of ABDF method
!  subroutine SetABDF(BDF, deg, tau)
!    type(ABDF), intent(inout) :: BDF
!    integer, intent(in) :: deg
!    real, dimension (1:deg+1), intent(in) :: tau
!    real :: theta0, theta1, theta2, A0, A1, A2, B0, B1, C0, C_BDF, C_EXT
!    real, dimension(:,:), allocatable :: lag
!    real :: factorial, prod
!    integer :: i,j,k, version
!
!    if(deg > BDF%max_deg) then
!       print*,'Error in SetABDF, try to set degree ',deg,&
!            ', but there is allocated the degree ',BDF%max_deg
!       stop
!    endif
!
!    version = 1
!    if(version == 1) then ! NEW version
!       BDF%deg = deg
!
!       ! Lagrangian interpolation \sum_{i=0}^{deg(+1)} w(i,:) L_i(t)
!
!       allocate(lag(0:deg+1, 0:deg+1) )
!
!       factorial = 1.
!       do i=0, deg+1
!          lag(i,i) = 0.
!          do j= i+1, deg + 1
!             lag(i, j) = sum(tau(i+1:j) ) ! lag(i,j) = t_{k-i} - t_{k-j}
!             lag(j, i) = - lag(i, j)
!          enddo
!          if(i> 0) factorial = factorial * i   ! factorial = (deg+1)!
!       enddo
!
!       do i= 0, deg
!          ! alpha(i) = (L_i)' (t_k)
!          BDF%alpha(i) = 0.
!
!          do j= 0, deg
!             if(j /= i) then
!                prod = 1./lag(i,j)
!
!                do k=0, deg
!                   if(k /= j .and. k/=i )  prod = prod * lag(0,k)/(lag(i,k))
!                enddo
!
!                BDF%alpha(i) = BDF%alpha(i) + prod
!             endif
!          enddo
!
!          BDF%alpha(i) = BDF%alpha(i) * tau(1)
!       enddo
!
!
!       ! coefficients delta for local error estimate
!       ! delta(i) = (deg+1)-th derivative of L_i, = const
!
!       lag(0:deg+1, 0:deg+1) = lag(0:deg+1, 0:deg+1) / tau(1)
!
!       do i=0, deg+1
!          lag(i,i) = 1.
!       enddo
!
!       do i=0, deg + 1
!          BDF%delta(i) =  factorial / product(lag(i, 0:deg+1) )
!       enddo
!
!       ! coefficient gamma
!       BDF%gamm = -1.
!       do i=1,deg
!          BDF%gamm = BDF%gamm * sum(tau(1:i)) /(tau(1) * (i+1))
!       enddo
!
!
!       !write(*,'(a7,i5,es12.4,a1,8es12.4)') 'alpha',deg,sum(BDF%alpha(0:deg)), '|', &
!       !     BDF%alpha(0:deg)
!       !write(*,'(a7,i5,es12.4,a1,8es12.4)') 'delta',deg,sum(BDF%delta(0:deg+1)),'|', &
!       !     BDF%delta(0:deg+1)
!       !write(*,'(a7,i5,2es12.4)') 'gamma',deg, BDF%gamm, 1./(deg+1)
!       !print*,'----------------------'
!
!
!       BDF%deg = deg
!       select case (deg)
!       case(:1)   ! EBDF extrapolated BDF
!          BDF%extrap(1) = 1.
!
!       case(2)
!          theta0 = tau(1)/tau(2)
!          theta1 = tau(2)/tau(3)
!
!          BDF%extrap(1) = theta0 + 1.
!          BDF%extrap(2) = -theta0
!
!       case(3)
!          theta0 = tau(1)/tau(2)
!          theta1 = tau(2)/tau(3)
!          theta2 = tau(3)/tau(4)
!
!          A0 = theta0 + 1
!          A1 = theta1*A0 + 1
!          A2 = theta2*A1 + 1
!
!          B0 = theta1 + 1
!          B1 = theta2*B0 + 1
!
!          BDF%extrap(1) = A0*A1/B0
!          BDF%extrap(2) = -theta0*A1
!          BDF%extrap(3) = theta0*theta1**2 *A0/B0
!
!       case(4:)
!
!          !if(state%time%iter <= 5) &
!          !     print*,'m-step extrapolated BDF scheme for m>3 is not implemented'
!          !stop
!
!       end select
!
!
!    else
!
!       BDF%deg = deg
!       select case (deg)
!       case(:1)   ! EBDF extrapolated BDF
!          theta0 = tau(1)/tau(2)
!
!          BDF%alpha(0:1) = (/ 1., -1./)
!
!          BDF%extrap(1) = 1.
!
!          BDF%Bextrap(1) = theta0 + 1.
!          BDF%Bextrap(2) = -theta0
!
!          ! 0^th order extrapolation
!          !BDF%Bextrap(1) = 1.
!          !BDF%Bextrap(2) = 0.
!
!          ! for estimation of the local discretization error via finite differences
!          BDF%delta(0) =  2.* theta0 /( 1. +  theta0)
!          BDF%delta(1) = -2.* theta0
!          BDF%delta(2) =  2.* theta0  * theta0 /( 1. +  theta0)
!
!
!          !C_BDF = -0.5
!          !C_EXT = 0.5*(theta0+1)/theta0
!          !BDF%gamm = C_BDF/( BDF%alpha(0)* C_EXT - C_BDF)
!
!          BDF%gamm = -theta0/(2*theta0+1)
!       case(2)
!          theta0 = tau(1)/tau(2)
!          theta1 = tau(2)/tau(3)
!
!          A0 = theta0 + 1
!          A1 = theta1*A0 + 1
!
!          B0 = theta1 + 1
!
!          BDF%alpha(1) = -A0
!          BDF%alpha(2) = theta0**2/A0
!
!          BDF%alpha(0) = -( BDF%alpha(1) + BDF%alpha(2))
!
!
!          BDF%extrap(1) = theta0 + 1.
!          BDF%extrap(2) = -theta0
!
!          BDF%Bextrap(1) = A0*A1/B0
!          BDF%Bextrap(2) = -theta0*A1
!          BDF%Bextrap(3) = theta0*theta1**2 *A0/B0
!
!
!          ! for estimation of the local discretization error via finite differences
!          BDF%delta(0) =  6.* theta0 * theta0 * theta1 /( 1. +  theta0) &
!               /(1+theta1+theta0*theta1)
!          BDF%delta(1) = -6.* theta0 * theta0 * theta1 /( 1. +  theta1)
!          BDF%delta(2) =  6.* theta0 * theta0 * theta0 * theta1 /( 1. +  theta0)
!          BDF%delta(3) = - 6.* (theta0 * theta1)**3 /( 1. +  theta1) &
!               /(1+theta1+theta0*theta1)
!
!
!
!          C_BDF = -A0/theta0
!          C_EXT = A0*A1/theta0**2/theta1
!
!          BDF%gamm = C_BDF/( BDF%alpha(0)* C_EXT - C_BDF)
!
!       case(3)
!          theta0 = tau(1)/tau(2)
!          theta1 = tau(2)/tau(3)
!          theta2 = tau(3)/tau(4)
!
!          A0 = theta0 + 1
!          A1 = theta1*A0 + 1
!          A2 = theta2*A1 + 1
!
!          B0 = theta1 + 1
!          B1 = theta2*B0 + 1
!
!          C0 = theta2 + 1
!
!          BDF%alpha(1) = -A0*A1/B0
!          BDF%alpha(2) = theta0**2*A1/A0
!          BDF%alpha(3) = -theta0**2 *theta1**3 *A0/B0/A1
!          BDF%alpha(0) = -( BDF%alpha(1) + BDF%alpha(2) + BDF%alpha(3))
!
!
!          BDF%extrap(1) = A0*A1/B0
!          BDF%extrap(2) = -theta0*A1
!          BDF%extrap(3) = theta0*theta1**2 *A0/B0
!
!          BDF%Bextrap(1) = A0/ B0 * A1 /B1 * A2
!          BDF%Bextrap(2) = -theta0/ C0 * A1 * A2
!          BDF%Bextrap(3) = theta0*theta1**2/B0 *A0*A2
!          BDF%Bextrap(4) = -theta0* theta1**2 * theta2**3* A0*A1 /C0 /B1
!
!          C_BDF = -A0*A1/theta0**2 /theta1
!
!          C_EXT = A0*A1*A2/theta0**3/theta1**2/theta2
!
!!!!print*,'C_EXT = ',C_EXT, C_EXT*tau(1)**4
!
!          BDF%gamm = C_BDF/( BDF%alpha(0)* C_EXT - C_BDF)
!
!
!
!       end select
!
!       if(deg > 3) then
!          print*,'m-step BDF scheme for m>3 is not implemented'
!          stop
!       endif
!
!       !write(*,'(a7,i5,8es12.4)') 'ALPHA',deg,sum(BDF%alpha(0:deg)), BDF%alpha(0:deg)
!       !write(*,'(a7,i5,8es12.4)') 'DELTA',deg,BDF%delta(0:deg+1)
!
!    endif
!
!
!  end subroutine SetABDF

!  !> set up a fake BDF formula to simply perform Langrange extrapolation
!  subroutine FakeBDFLagrangeExtrapolation(BDF, deg, tau)
!    type(ABDF), intent(inout) :: BDF
!    integer, intent(in) :: deg
!    real, dimension (1:deg+1), intent(in) :: tau
!    integer:: n,i
!
!    BDF%deg = deg + 1
!    n = size(tau)
!    ! FIXME: This is the classical formula. Go via DD eventually!
!    do i = 1,n
!       ! the leading minus is due to minus used in BDF extrapolation
!       BDF%alpha(i) = -product(-tau/(tau(i) - tau),mask=tau/=tau(i))
!    end do
!  end subroutine FakeBDFLagrangeExtrapolation


  !> generates Legenre polynomials in the Gauss quadrature nodes
  subroutine Init_Legenre_pols(G_rule, max_deg)
    type(Gauss_rule), target, intent(inout) :: G_rule
    integer, intent(in) :: max_deg ! = MaxDegreeImplemented
    real, dimension(:), pointer :: xi ! Gauss integ nodes
    real, dimension(:), allocatable :: ti ! Gauss integ nodes transformed to (-1,1)
    integer::  i, Qdof

    Qdof = G_rule%Qdeg

    allocate(G_rule%Leg_phi(0:max_deg,0:Qdof), ti(1:Qdof) )  ! ti(0,*) = L^2 norm of Leg_phi

    xi => G_rule%lambda(1:Qdof)
    ti(1:Qdof) = 2*xi(1:Qdof) - 1.


    G_rule%Leg_phi(0, 1:Qdof) = 1.
    G_rule%Leg_phi(1, 1:Qdof) = ti(1:Qdof)

    !write(*,'(a6,60es12.4)') 'ti: ',ti(1:Qdof)

    do i=1, max_deg -1
       G_rule%Leg_phi(i+1, 1:Qdof) = (2.*i + 1) * ti(1:Qdof) * G_rule%Leg_phi(i, 1:Qdof) &
            - 1. * i *  G_rule%Leg_phi(i-1, 1:Qdof)

       G_rule%Leg_phi(i+1, 1:Qdof) =  G_rule%Leg_phi(i+1, 1:Qdof) / (i+1) !/ 2**0.5

    enddo

    ! setting of the normalization factors
    do i=0, max_deg
       G_rule%Leg_phi(i,0) = dot_product(G_rule%weights(1:Qdof), &
            G_rule%Leg_phi(i, 1:Qdof) * G_rule%Leg_phi(i, 1:Qdof) )
    enddo

    deallocate(ti)

  end subroutine Init_Legenre_pols

! FRERROR problem with Qdof - it is an argument but defined as local too
  !> evaluation of the Legenre polynomials in the given set of nodes
  subroutine Eval_Legenre_pols(G_rule, dof, Qdof, xi, phi)
    type(Gauss_rule), intent(inout) :: G_rule
    integer, intent(in) :: dof ! = number of legendre functions
    integer, intent(in) :: Qdof
    real, dimension(1:Qdof), intent(in) :: xi ! sought nodes
    real, dimension(0:dof, 1:Qdof), intent(inout) :: phi ! Gauss integ nodes
    real, dimension(:), allocatable :: ti ! Gauss integ nodes transformed to (-1,1)
    integer:: i!, Qdof

    allocate( ti(1:Qdof) )

    ti(1:Qdof) = 2*xi(1:Qdof) - 1.

    !print*,'###',Qdof, size(phi , 1), dof
    phi(0, 1:Qdof) = 1.
    if(dof > 0) &
         phi(1, 1:Qdof) = ti(1:Qdof)

    !write(*,'(a6,60es12.4)') 'ti: ',ti(1:Qdof)

    do i=1, dof -1
       phi(i+1, 1:Qdof) = (2.*i + 1) * ti(1:Qdof) * phi(i, 1:Qdof) &
            - 1. * i * phi(i-1, 1:Qdof)

       phi(i+1, 1:Qdof) =  phi(i+1, 1:Qdof) / (i+1) !/ 2**0.5

    enddo

    deallocate(ti)

  end subroutine Eval_Legenre_pols

end module integration


