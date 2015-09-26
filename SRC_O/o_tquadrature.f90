
module tquadrature_mod

  use quadrature_mod
  use paramets
  use lapack_oper

  implicit none

     !>  rule for the time intervals
   type :: Time_rule
      integer::  Qdeg         !  type of a rule, accurate for 2*Qdeg-1 polynom / resp 2*Qdeg - 2 for Radau
    !  integer::  Qdof         !  number of integration nodes NOT USED same as Qdeg
      real, allocatable, dimension(:) :: weights
      real, allocatable, dimension(:) :: lambda
      real, dimension(:,:), allocatable :: phi
      real, dimension(:,:), allocatable :: Dphi

      real, dimension(:), allocatable :: RadauPol !only for Radau Quadrature

   contains

      procedure :: createTimeRule => createAbstractTimeRule
!      procedure :: EvalTrulePhi
      procedure, non_overridable :: SetTruleBasis
      procedure, non_overridable :: WriteTruleBasis

  end type Time_rule

  type, extends( Time_rule ) :: RadauTime_rule

   contains
      procedure :: createTimeRule => createRadauRule
      procedure :: evalRadauPolynomial
      procedure :: testRadauPol
  end type

  type, extends( Time_rule ) :: GaussTime_rule

   contains
      procedure :: createTimeRule => createTimeGaussRule
  end type

  public :: evalTphi
  public :: evalTdphi

  public :: SetLegendrePolynomials
  public :: Eval_LegendrePolynomialsOneNode

  public :: evalTimeFunctionValue
  public :: evalTimeFunctionDerivativeValue

  contains

  !> abstract - shouldnt be called
   subroutine createAbstractTimeRule( this, Qdeg)
      class( Time_rule ), intent( inout ) :: this
      integer, intent( in ) :: Qdeg

      print*, 'this is only abstract procedure - should not be used, use Radau or Gauss T_rule', &
         'Works only as a template for child types'
      stop

   end subroutine createAbstractTimeRule

  !TIME RULE
   !> generates one Time Gauss quadrature rule (same as Gauss rule)
   subroutine createTimeGaussRule( this, Qdeg)
      class( GaussTime_rule ), intent( inout ) :: this
      integer, intent( in ) :: Qdeg

      this%Qdeg = max(Qdeg,1)
   !   this%Qdof = this%Qdeg


    allocate(this%weights(1:this%Qdeg))
    allocate(this%lambda(1:this%Qdeg))

    select case (Qdeg)
    case(:1)
       this%weights(1) = 1.
       this%lambda(1) = 0.5
    case(2)
       this%weights(1:2) = 0.5
       this%lambda(1) = (1.-0.57735026918962576451)/2
       this%lambda(2) = 1.57735026918962576451/2
    case(3)
       this%weights(1) = 5./18
       this%weights(2) = 8./18
       this%weights(3) = 5./18
       this%lambda(1) = (1.-0.77459666924148337704)/2
       this%lambda(2) = 0.5
       this%lambda(3) = 1.77459666924148337704/2
    case(4)
       this%weights(1) = 0.34785484513745385737/2
       this%weights(2) = 0.65214515486254614263/2
       this%weights(3) = this%weights(2)
       this%weights(4) = this%weights(1)
       this%lambda(1) = (1-0.86113631159405257522)/2
       this%lambda(2) = (1-0.33998104358485626480)/2
       this%lambda(3) = 1.33998104358485626480/2
       this%lambda(4) = 1.86113631159405257522/2
    case(5)
       this%weights(1) = 0.23692688505618908751/2.
       this%weights(2) = 0.47862867049936646804/2.
       this%weights(3) = 0.56888888888888888889/2.
       this%weights(4) = this%weights(2)
       this%weights(5) = this%weights(1)

       this%lambda(1) = (1-0.90617984593866399280)/2.
       this%lambda(2) = (1-0.53846931010568309104)/2.
       this%lambda(3) = 0.5
       this%lambda(4) = 1.53846931010568309104/2.
       this%lambda(5) = 1.90617984593866399280/2.
    case(6)
       this%weights(1) = 0.17132449237917034504/2.
       this%weights(2) = 0.36076157304813860757/2.
       this%weights(3) = 0.46791393457269104739/2.
       this%weights(4) = this%weights(3)
       this%weights(5) = this%weights(2)
       this%weights(6) = this%weights(1)

       this%lambda(1) = (1.-0.93246951420315202781)/2.
       this%lambda(2) = (1.-0.66120938646626451366)/2.
       this%lambda(3) = (1.-0.23861918608319690863)/2.
       this%lambda(4) =  1.23861918608319690863/2.
       this%lambda(5) =  1.66120938646626451366/2.
       this%lambda(6) =  1.93246951420315202781/2.
     case(7)
       this%weights(1) = 0.12948496616886969327/2.
       this%weights(2) = 0.27970539148927666790/2.
       this%weights(3) = 0.38183005050511894495/2.
       this%weights(4) = 0.41795918367346938776/2.
       this%weights(5) = this%weights(3)
       this%weights(6) = this%weights(2)
       this%weights(7) = this%weights(1)

       this%lambda(1) = (1.-0.94910791234275852453)/2.
       this%lambda(2) = (1.-0.74153118559939443986)/2.
       this%lambda(3) = (1.-0.40584515137739716691)/2.
       this%lambda(4) =  0.5
       this%lambda(5) =  1.40584515137739716691/2.
       this%lambda(6) =  1.74153118559939443986/2.
       this%lambda(7) =  1.94910791234275852453/2.

      case( 8)
       this%weights( 1) =  0.5061426814518818E-01
       this%weights( 2) =  0.1111905172266872E+00
       this%weights( 3) =  0.1568533229389436E+00
       this%weights( 4) =  0.1813418916891809E+00
       this%weights( 5) =  0.1813418916891809E+00
       this%weights( 6) =  0.1568533229389436E+00
       this%weights( 7) =  0.1111905172266872E+00
       this%weights( 8) =  0.5061426814518818E-01

       this%lambda( 1) =  0.1985507175123191E-01
       this%lambda( 2) =  0.1016667612931866E+00
       this%lambda( 3) =  0.2372337950418355E+00
       this%lambda( 4) =  0.4082826787521751E+00
       this%lambda( 5) =  0.5917173212478249E+00
       this%lambda( 6) =  0.7627662049581645E+00
       this%lambda( 7) =  0.8983332387068134E+00
       this%lambda( 8) =  0.9801449282487681E+00

      case( 9)
       this%weights( 1) =  0.4063719418078721E-01
       this%weights( 2) =  0.9032408034742863E-01
       this%weights( 3) =  0.1303053482014677E+00
       this%weights( 4) =  0.1561735385200015E+00
       this%weights( 5) =  0.1651196775006299E+00
       this%weights( 6) =  0.1561735385200015E+00
       this%weights( 7) =  0.1303053482014677E+00
       this%weights( 8) =  0.9032408034742863E-01
       this%weights( 9) =  0.4063719418078721E-01

       this%lambda( 1) =  0.1591988024618696E-01
       this%lambda( 2) =  0.8198444633668212E-01
       this%lambda( 3) =  0.1933142836497048E+00
       this%lambda( 4) =  0.3378732882980955E+00
       this%lambda( 5) =  0.5000000000000000E+00
       this%lambda( 6) =  0.6621267117019045E+00
       this%lambda( 7) =  0.8066857163502952E+00
       this%lambda( 8) =  0.9180155536633179E+00
       this%lambda( 9) =  0.9840801197538130E+00

      case(10)
       this%weights( 1) =  0.3333567215434402E-01
       this%weights( 2) =  0.7472567457529031E-01
       this%weights( 3) =  0.1095431812579911E+00
       this%weights( 4) =  0.1346333596549982E+00
       this%weights( 5) =  0.1477621123573765E+00
       this%weights( 6) =  0.1477621123573765E+00
       this%weights( 7) =  0.1346333596549982E+00
       this%weights( 8) =  0.1095431812579911E+00
       this%weights( 9) =  0.7472567457529031E-01
       this%weights(10) =  0.3333567215434402E-01

       this%lambda( 1) =  0.1304673574141413E-01
       this%lambda( 2) =  0.6746831665550773E-01
       this%lambda( 3) =  0.1602952158504878E+00
       this%lambda( 4) =  0.2833023029353764E+00
       this%lambda( 5) =  0.4255628305091844E+00
       this%lambda( 6) =  0.5744371694908156E+00
       this%lambda( 7) =  0.7166976970646236E+00
       this%lambda( 8) =  0.8397047841495122E+00
       this%lambda( 9) =  0.9325316833444923E+00
       this%lambda(10) =  0.9869532642585859E+00

      case(11)
       this%weights( 1) =  0.2783428355808687E-01
       this%weights( 2) =  0.6279018473245233E-01
       this%weights( 3) =  0.9314510546386705E-01
       this%weights( 4) =  0.1165968822959953E+00
       this%weights( 5) =  0.1314022722551234E+00
       this%weights( 6) =  0.1364625433889503E+00
       this%weights( 7) =  0.1314022722551234E+00
       this%weights( 8) =  0.1165968822959953E+00
       this%weights( 9) =  0.9314510546386705E-01
       this%weights(10) =  0.6279018473245233E-01
       this%weights(11) =  0.2783428355808687E-01

       this%lambda( 1) =  0.1088567092697151E-01
       this%lambda( 2) =  0.5646870011595234E-01
       this%lambda( 3) =  0.1349239972129753E+00
       this%lambda( 4) =  0.2404519353965941E+00
       this%lambda( 5) =  0.3652284220238275E+00
       this%lambda( 6) =  0.5000000000000000E+00
       this%lambda( 7) =  0.6347715779761725E+00
       this%lambda( 8) =  0.7595480646034058E+00
       this%lambda( 9) =  0.8650760027870247E+00
       this%lambda(10) =  0.9435312998840477E+00
       this%lambda(11) =  0.9891143290730284E+00

      case(12)
       this%weights( 1) =  0.2358766819325595E-01
       this%weights( 2) =  0.5346966299765912E-01
       this%weights( 3) =  0.8003916427167312E-01
       this%weights( 4) =  0.1015837133615329E+00
       this%weights( 5) =  0.1167462682691774E+00
       this%weights( 6) =  0.1245735229067013E+00
       this%weights( 7) =  0.1245735229067013E+00
       this%weights( 8) =  0.1167462682691774E+00
       this%weights( 9) =  0.1015837133615329E+00
       this%weights(10) =  0.8003916427167312E-01
       this%weights(11) =  0.5346966299765912E-01
       this%weights(12) =  0.2358766819325595E-01

       this%lambda( 1) =  0.9219682876640378E-02
       this%lambda( 2) =  0.4794137181476255E-01
       this%lambda( 3) =  0.1150486629028477E+00
       this%lambda( 4) =  0.2063410228566913E+00
       this%lambda( 5) =  0.3160842505009099E+00
       this%lambda( 6) =  0.4373832957442655E+00
       this%lambda( 7) =  0.5626167042557344E+00
       this%lambda( 8) =  0.6839157494990901E+00
       this%lambda( 9) =  0.7936589771433087E+00
       this%lambda(10) =  0.8849513370971523E+00
       this%lambda(11) =  0.9520586281852375E+00
       this%lambda(12) =  0.9907803171233596E+00

      case(13)
       this%weights( 1) =  0.2024200238265791E-01
       this%weights( 2) =  0.4606074991886418E-01
       this%weights( 3) =  0.6943675510989367E-01
       this%weights( 4) =  0.8907299038097291E-01
       this%weights( 5) =  0.1039080237684442E+00
       this%weights( 6) =  0.1131415901314486E+00
       this%weights( 7) =  0.1162757766154369E+00
       this%weights( 8) =  0.1131415901314486E+00
       this%weights( 9) =  0.1039080237684442E+00
       this%weights(10) =  0.8907299038097291E-01
       this%weights(11) =  0.6943675510989367E-01
       this%weights(12) =  0.4606074991886418E-01
       this%weights(13) =  0.2024200238265791E-01

       this%lambda( 1) =  0.7908472640705932E-02
       this%lambda( 2) =  0.4120080038851104E-01
       this%lambda( 3) =  0.9921095463334506E-01
       this%lambda( 4) =  0.1788253302798299E+00
       this%lambda( 5) =  0.2757536244817765E+00
       this%lambda( 6) =  0.3847708420224326E+00
       this%lambda( 7) =  0.5000000000000000E+00
       this%lambda( 8) =  0.6152291579775674E+00
       this%lambda( 9) =  0.7242463755182235E+00
       this%lambda(10) =  0.8211746697201701E+00
       this%lambda(11) =  0.9007890453666549E+00
       this%lambda(12) =  0.9587991996114890E+00
       this%lambda(13) =  0.9920915273592941E+00

      case(14)
       this%weights( 1) =  0.1755973016587599E-01
       this%weights( 2) =  0.4007904357988007E-01
       this%weights( 3) =  0.6075928534395157E-01
       this%weights( 4) =  0.7860158357909681E-01
       this%weights( 5) =  0.9276919873896890E-01
       this%weights( 6) =  0.1025992318606479E+00
       this%weights( 7) =  0.1076319267315789E+00
       this%weights( 8) =  0.1076319267315789E+00
       this%weights( 9) =  0.1025992318606479E+00
       this%weights(10) =  0.9276919873896890E-01
       this%weights(11) =  0.7860158357909681E-01
       this%weights(12) =  0.6075928534395157E-01
       this%weights(13) =  0.4007904357988007E-01
       this%weights(14) =  0.1755973016587599E-01

       this%lambda( 1) =  0.6858095651593843E-02
       this%lambda( 2) =  0.3578255816821324E-01
       this%lambda( 3) =  0.8639934246511749E-01
       this%lambda( 4) =  0.1563535475941573E+00
       this%lambda( 5) =  0.2423756818209230E+00
       this%lambda( 6) =  0.3404438155360551E+00
       this%lambda( 7) =  0.4459725256463282E+00
       this%lambda( 8) =  0.5540274743536718E+00
       this%lambda( 9) =  0.6595561844639448E+00
       this%lambda(10) =  0.7576243181790770E+00
       this%lambda(11) =  0.8436464524058427E+00
       this%lambda(12) =  0.9136006575348825E+00
       this%lambda(13) =  0.9642174418317868E+00
       this%lambda(14) =  0.9931419043484062E+00

     case(15:)
       this%weights( 1) =  0.1537662099805857E-01
       this%weights( 2) =  0.3518302374405406E-01
       this%weights( 3) =  0.5357961023358598E-01
       this%weights( 4) =  0.6978533896307712E-01
       this%weights( 5) =  0.8313460290849704E-01
       this%weights( 6) =  0.9308050000778109E-01
       this%weights( 7) =  0.9921574266355583E-01
       this%weights( 8) =  0.1012891209627806E+00
       this%weights( 9) =  0.9921574266355583E-01
       this%weights(10) =  0.9308050000778109E-01
       this%weights(11) =  0.8313460290849704E-01
       this%weights(12) =  0.6978533896307712E-01
       this%weights(13) =  0.5357961023358598E-01
       this%weights(14) =  0.3518302374405406E-01
       this%weights(15) =  0.1537662099805857E-01

       this%lambda( 1) =  0.6003740989757256E-02
       this%lambda( 2) =  0.3136330379964708E-01
       this%lambda( 3) =  0.7589670829478640E-01
       this%lambda( 4) =  0.1377911343199150E+00
       this%lambda( 5) =  0.2145139136957306E+00
       this%lambda( 6) =  0.3029243264612183E+00
       this%lambda( 7) =  0.3994029530012828E+00
       this%lambda( 8) =  0.5000000000000000E+00
       this%lambda( 9) =  0.6005970469987173E+00
       this%lambda(10) =  0.6970756735387817E+00
       this%lambda(11) =  0.7854860863042694E+00
       this%lambda(12) =  0.8622088656800850E+00
       this%lambda(13) =  0.9241032917052137E+00
       this%lambda(14) =  0.9686366962003530E+00
       this%lambda(15) =  0.9939962590102427E+00

    end select

    if(Qdeg > 15) then
      print*,'Warning !! Time quadrature rule of degree ', Qdeg, &
      ' is not implemented'
      print*,'           the maximal default quadrature rule of degree', &
      ' 15 is used'
    endif

   end subroutine createTimeGaussRule

 !> generates one Time Radau quadrature rule
   subroutine createRadauRule( this, Qdeg)
      class( RadauTime_rule ), intent( inout ) :: this
      integer, intent( in ) :: Qdeg

    this%Qdeg = max(Qdeg, 1)
  !  this%Qdof = this%Qdeg

    allocate( this%weights( 1:this%Qdeg ) )
    allocate( this%lambda( 1:this%Qdeg ) )

    select case (Qdeg)
    case(:1)
       this%weights(1) = 1.
       this%lambda(1) = 1.
    case(2)
       this%weights(1) =   0.7500000000000000
       this%weights(2) =   0.2500000000000000
       this%lambda(1) =   0.3333333333333333
       this%lambda(2) =   1.0000000000000000
    case(3)
       this%weights(1) =   0.3764030627004669
       this%weights(2) =   0.5124858261884215
       this%weights(3) =   0.1111111111111111
       this%lambda(1) =   0.1550510257216822
       this%lambda(2) =   0.6449489742783179
       this%lambda(3) =   1.0000000000000000
    case(4)
       this%weights(1) =   0.2204622111767679
       this%weights(2) =   0.3881934688431719
       this%weights(3) =   0.3288443199800598
       this%weights(4) =   0.0625000000000000
       this%lambda(1) =   0.0885879595127039
       this%lambda(2) =   0.4094668644407347
       this%lambda(3) =   0.7876594617608470
       this%lambda(4) =   1.0000000000000000
    case(5)
       this%weights(1) =   0.1437135607912257
       this%weights(2) =   0.2813560151494621
       this%weights(3) =   0.3118265229757412
       this%weights(4) =   0.2231039010835709
       this%weights(5) =   0.0400000000000000
       this%lambda(1) =   0.0571041961145177
       this%lambda(2) =   0.2768430136381238
       this%lambda(3) =   0.5835904323689168
       this%lambda(4) =   0.8602401356562195
       this%lambda(5) =   1.0000000000000000
    case(6)
       this%weights(1) =   0.1007941926267402
       this%weights(2) =   0.2084506671559546
       this%weights(3) =   0.2604633915947875
       this%weights(4) =   0.2426935942344849
       this%weights(5) =   0.1598203766102554
       this%weights(6) =   0.0277777777777778
       this%lambda(1) =   0.0398098570514687
       this%lambda(2) =   0.1980134178736082
       this%lambda(3) =   0.4379748102473862
       this%lambda(4) =   0.6954642733536361
       this%lambda(5) =   0.9014649142011736
       this%lambda(6) =   1.0000000000000000
    case(7)
       this%weights(1) =   0.0744942355560100
       this%weights(2) =   0.1591021157336509
       this%weights(3) =   0.2123518895029778
       this%weights(4) =   0.2235549145072832
       this%weights(5) =   0.1904749368221156
       this%weights(6) =   0.1196137446126562
       this%weights(7) =   0.0204081632653061
       this%lambda(1) =   0.0293164271597849
       this%lambda(2) =   0.1480785996684843
       this%lambda(3) =   0.3369846902811543
       this%lambda(4) =   0.5586715187715501
       this%lambda(5) =   0.7692338620300545
       this%lambda(6) =   0.9269456713197410
       this%lambda(7) =   1.0000000000000000
    case(8)
       this%weights(1) =   0.0572544073721283
       this%weights(2) =   0.1248239506649321
       this%weights(3) =   0.1735073978172506
       this%weights(4) =   0.1957860837262468
       this%weights(5) =   0.1882587726945593
       this%weights(6) =   0.1520653103233925
       this%weights(7) =   0.0926790774014897
       this%weights(8) =   0.0156250000000000
       this%lambda(1) =   0.0224793864387125
       this%lambda(2) =   0.1146790531609042
       this%lambda(3) =   0.2657898227845895
       this%lambda(4) =   0.4528463736694446
       this%lambda(5) =   0.6473752828868303
       this%lambda(6) =   0.8197593082631076
       this%lambda(7) =   0.9437374394630779
       this%lambda(8) =   1.0000000000000000
    case(9)
       this%weights(1) =   0.0453572524616435
       this%weights(2) =   0.1002766490122754
       this%weights(3) =   0.1431933481786156
       this%weights(4) =   0.1688469834879647
       this%weights(5) =   0.1741365013864834
       this%weights(6) =   0.1584218878352190
       this%weights(7) =   0.1235946891022966
       this%weights(8) =   0.0738270095231576
       this%weights(9) =   0.0123456790123457
       this%lambda(1) =   0.0177799151473634
       this%lambda(2) =   0.0913236078997939
       this%lambda(3) =   0.2143084793956307
       this%lambda(4) =   0.3719321645832723
       this%lambda(5) =   0.5451866848034267
       this%lambda(6) =   0.7131752428555695
       this%lambda(7) =   0.8556337429578544
       this%lambda(8) =   0.9553660447100302
       this%lambda(9) =   1.0000000000000000
    case(10)
       this%weights(1) =   0.0368085027433803
       this%weights(2) =   0.0821880063684609
       this%weights(3) =   0.1195967158571898
       this%weights(4) =   0.1453050824164592
       this%weights(5) =   0.1567912286134692
       this%weights(6) =   0.1529296438622114
       this%weights(7) =   0.1340974189205893
       this%weights(8) =   0.1021350659395004
       this%weights(9) =   0.0601483352787409
       this%weights(10) =   0.0100000000000000
       this%lambda(1) =   0.0144124096488766
       this%lambda(2) =   0.0743873897091961
       this%lambda(3) =   0.1761166561629953
       this%lambda(4) =   0.3096675799276378
       this%lambda(5) =   0.4619704010810110
       this%lambda(6) =   0.6181172346952940
       this%lambda(7) =   0.7628230151850396
       this%lambda(8) =   0.8819210212100013
       this%lambda(9) =   0.9637421871167906
       this%lambda(10) =   1.0000000000000000
    case(11)
       this%weights(1) =   0.0304625489060658
       this%weights(2) =   0.0685168410666014
       this%weights(3) =   0.1010815542700119
       this%weights(4) =   0.1254626888485642
       this%weights(5) =   0.1396806665516916
       this%weights(6) =   0.1425827819705036
       this%weights(7) =   0.1339335430948421
       this%weights(8) =   0.1144330619244883
       this%weights(9) =   0.0856588096033299
       this%weights(10) =   0.0499230409539840
       this%weights(11) =   0.0082644628099174
       this%lambda(1) =   0.0119176134324156
       this%lambda(2) =   0.0617320718771481
       this%lambda(3) =   0.1471114496430702
       this%lambda(4) =   0.2611596760084562
       this%lambda(5) =   0.3946398468857868
       this%lambda(6) =   0.5367387657156606
       this%lambda(7) =   0.6759444616766651
       this%lambda(8) =   0.8009789210368988
       this%lambda(9) =   0.9017109877901468
       this%lambda(10) =   0.9699709678385136
       this%lambda(11) =   1.0000000000000000
    case(12:)
       this%weights(1) =   0.0256240496036348
       this%weights(2) =   0.0579537401458690
       this%weights(3) =   0.0863853196566543
       this%weights(4) =   0.1089344395130960
       this%weights(5) =   0.1240607804020050
       this%weights(6) =   0.1307328302760666
       this%weights(7) =   0.1284956690763539
       this%weights(8) =   0.1175015575724929
       this%weights(9) =   0.0984992674130449
       this%weights(10) =   0.0727818344269976
       this%weights(11) =   0.0420860674693405
       this%weights(12) =   0.0069444444444444
       this%lambda(1) =   0.0100182804616804
       this%lambda(2) =   0.0520354511271806
       this%lambda(3) =   0.1246192251444431
       this%lambda(4) =   0.2228406070438378
       this%lambda(5) =   0.3400081579146652
       this%lambda(6) =   0.4681376130895840
       this%lambda(7) =   0.5984972797671392
       this%lambda(8) =   0.7222032848909679
       this%lambda(9) =   0.8308248996228186
       this%lambda(10) =   0.9169583865525948
       this%lambda(11) =   0.9747263796024797
       this%lambda(12) =   1.0000000000000000

    end select

    if(Qdeg > 12) then
      print*,'Warning !! Radau rule of degree ', Qdeg, &
      ' is not implemented'
      print*,'           the maximal default quadrature rule of degree', &
      ' 12 is used'
    endif


   end subroutine createRadauRule


  ! setting of time test functions \f$\phi(t)\f$ in time integ nodes
  subroutine SetTruleBasis ( this )
    class( Time_rule ), intent(inout) :: this
    integer:: Qdeg, dof, j, i
    real :: t
    real, dimension(:), allocatable :: xi
!    real, dimension(:,:), allocatable :: phi
!    real, dimension(:,:), allocatable :: dphi

!old
!    this%phi(1:dof, -2) = 0.
!    this%Dphi(1:dof, -2) = 0.
!
!   do j = -1, Qdeg ! in all quad nodes + starting and ending point of the interval
!      if(j == -1) then
!          t = 0.
!      elseif(j == 0) then
!          t =1.
!      else
!          t = this%lambda(j)
!      endif
!
!
!      this%phi( 1:dof , j ) = evalTphi( (/ (i, i = 1, dof) /) , t )
!      this%Dphi( 1:dof,j ) = evalTdphi( (/ (i, i=1, dof ) /) , t )
!
!
!   end do ! j


   !New version using LegendrePolynomials
   Qdeg = this%Qdeg ! number of nodes of the quadrature
   dof =  MaxTimeDegree + 2 ! 2nd +1 for Tdof_plus = number of basis functions

   ! time basis functions
   allocate( this%phi(1:dof, -2:Qdeg) ) ! at Gauss/Radau integ nodes, -1,0  =  phi at end points, -2 special point if needed
   allocate( this%Dphi(1:dof, -2:Qdeg) )
   allocate( xi(-1:Qdeg) )
!   allocate( phi(1:dof, -1:Qdeg) )
!   allocate( dphi(1:dof, -1:Qdeg) )

   do j = -1, Qdeg ! in all quad nodes + starting and ending point of the interval
      if(j == -1) then
          xi(j) = 0.
      elseif(j == 0) then
          xi(j) =1.
      else
          xi(j) = this%lambda(j)
      endif
   end do ! j

   call SetLegendrePolynomials( dof, Qdeg+2, xi(-1:Qdeg), this%phi(1:dof, -1:Qdeg),  this%Dphi(1:dof, -1:Qdeg) )


      !TEST against the older version
     ! phi(1:dof,-1:Qdeg) = phi(1:dof,-1:Qdeg) - this%phi( 1:dof , -1:Qdeg )
     ! dphi(1:dof,-1:Qdeg) =  this%dphi( 1:dof , -1:Qdeg ) - dphi(1:dof,-1:Qdeg)

      ! values of phi in integ nodes
!    print*,'&&&& dof =',dof,'Qdeg = ',Qdeg
!    print*, 'lambda bounds: ', lbound(this%lambda(:)), ubound(this%lambda(:))
   !  write(*,'(a6,30es13.5)') 'nodes:',this%lambda(1:Qdeg)
!     do i=1,dof
!      do j = -1,Qdeg
!         if ( abs(phi(i,j)) > 1.0e-12 ) &
!          write(*,*) 'phi_',i,phi(i,j)
!         if ( abs(dphi(i,j)) > 1.0e-12 ) &
!          write(*,*) 'Dphi_',i,Dphi(i,j)
!      enddo
!
!     enddo
!     print*,'.....................'
!
!      !test of the L^2-orthogonality
!     do i=1,min(dof, Qdeg)
!        do j=1,min(dof, Qdeg)
!           t = dot_product( this%weights(1:Qdeg), this%phi(i,1:Qdeg)*this%phi(j,1:Qdeg) )
!           !write(*,*)  this%weights(1:Qdeg), phi(i,1:Qdeg)
!           write(*,'(a6,2i5,es16.8)') 'L^2 ort',i,j,t
!           enddo
!        print*,'--------------------------------------------'
!     enddo


     deallocate ( xi )

  end subroutine SetTruleBasis

  subroutine WriteTruleBasis( this )
      class( Time_rule ), intent(in) :: this
      integer :: i,j , Qdeg

      print*, 'Printing T_rule basis to file T_rulePhi.'

      open (58, file="T_rulePhi", action ="write", status="replace")

     ! do j = 1, maxTrule
        ! T_rule => state%time%T_rule(j)
         Qdeg = this%Qdeg
       !  Qdof = T_rule%Qdof

      write(58,*) , 'Trule:' , Qdeg
      do i = 1, MaxTimeDegree + 2
         write(58, * ) , this%phi(i,:)
      enddo
      write(58,*) '---------------'

     ! enddo !j
      close(58)

  end subroutine WriteTruleBasis


  subroutine evalRadauPolynomial( this )
   class( RadauTime_rule ), intent(inout) :: this

   real,dimension(:,:), allocatable :: phi
   integer :: i, n

   n = this%Qdeg + 1 ! dim of the system, # of basis coeficients

   allocate ( this%RadauPol(1:n) )
   allocate ( phi(1:n, 1:n) ) ! needs to be transposed this%phi

   !solve the moment problem
   phi(1,1:n) = this%phi(1:n, -1 ) ! value in t=0

   do i =2,n
      phi(i,1:n) = this%phi(1:n, i-1 ) !other values \phi( t_i )
   enddo

   !rhs
   this%RadauPol(1) = 1.0
   this%RadauPol(2:n) = 0.0

   ! in RadauPol are saved the coeffs of RadauPol against the standard Legendre basis functions
   call SolveLocalMatrixProblem(n, phi(1:n,1:n), 1, this%RadauPol(1:n) )


   deallocate (phi)





  end subroutine evalRadauPolynomial

  subroutine testRadauPol( this )
   class( RadauTime_rule ), intent(inout) :: this

   integer :: i,n
   real :: f

   n = this%Qdeg+1

   f = evalTimeFunctionValue( n, this%RadauPol(1:n), 0.0 )
   print*, 'radau(0) =' ,f

   do i = 1,this%Qdeg
      f = evalTimeFunctionValue( n, this%RadauPol(1:n), this%lambda(i) )
      print*, 'radau(t_',i,') =' ,f
   end do



  end subroutine testRadauPol


!!!! public functions !!!!!!!!!!


!not used
  elemental function evalTphi( j, t ) result ( val )
   real :: val
   integer, intent( in ) :: j ! j-th basis function
   real, intent ( in ) :: t ! time node in interval [0,1]

   select case ( j )
      case(1)
         val = 1.
      case(2)
         val = 12.**0.5 *(t - 0.5)
      case(3)
         !val = 30.*2.**0.5 *(t*t - t  + 1./6.)
         val = 3.*20.**0.5 *(t*t - t  + 1./6.)
      case(4)
         val = 7.**0.5 * (20. * t**(3.) - 30. * t**(2.) + 12. * t - 1. )
      case(5)
         val = 3. * ( 70. * t**(4.) - 140. * t**(3.) + 90. * t**(2.) - 20. * t + 1. )
      case(6)
         val = 11.**0.5 * ( 252. * t**(5.) - 630. * t**(4.) + 560. * t**(3.) - 210. * t**(2.) + 30. * t - 1. )
      case(7:)
!         print*,' order >= 7 in evalTphi in tquadrature.f90 not implemented'
!         stop ' order >= 7 in EvalTphi in o_tquadrature.f90 not implemented'
      end select

  end function evalTphi

!not used
   elemental function evalTdphi( j, t ) result ( dval )
      real :: dval
      integer, intent( in ) :: j ! j-th basis function
      real, intent ( in ) :: t ! time node in interval [0,1]


      select case ( j )
      case(1)
         dval = 0.
      case(2)
         dval = 12.**0.5
      case(3)
         dval = 3.*20.**0.5 * (2.*t - 1.)
      case(4)
         dval = 7.**0.5 * (60. * t*t - 60. * t + 12.)
      case(5)
         dval = 60.*( 14. * t**(3.) - 21.  * t**(2.) + 9.  * t - 1. )
      case(6)
         dval = 11.**0.5 *(1260. * t**(4.) - 2520. *t**(3.) + 1680. *t**(2.) - 420. * t + 30. )
      case(7:)
 !        print*,' order >= 7 in SetTrulePhi in tquadrature.f90 not implemented'
!         stop ' order >= 7 in EvalTdPhi in o_tquadrature.f90 not implemented'
      end select


   end function evalTdphi


     !> evaluation of the Legenre polynomials and its derivatives in the given set of nodes
  !> eval values of ON Legendre polynomials on (0,1) in the given set of nodes in [0,1]
  pure subroutine Eval_LegendrePolynomialsOneNode( dof, xi, phi, dphi)
    !class(Gauss_rule), intent(inout) :: this
    integer, intent(in) :: dof ! = number of legendre functions
    real, intent(in) :: xi ! sought node
    real, dimension(0:dof-1), intent(out) :: phi ! Gauss integ nodes
    real, dimension(0:dof-1), intent(out), optional :: dphi ! Gauss integ nodes

    real :: ti ! Gauss integ nodes transformed to (-1,1)
    integer:: i,j
    real :: eps

   !print*, 'Qdof = ', Qdof

   ti = 2.0*xi - 1. ! from 0,1 to -1,1

   phi(0) = 1.
   if ( dof > 1) &
      phi(1) = ti

   do i=1,dof-2

    phi(i+1) = (2.0*i + 1.0) * ti * phi(i) &
         - ( 1. * i * phi(i-1) )

    phi(i+1) =  phi(i+1) / (i+1)
   enddo



   if ( present( dphi ) ) then

      eps = 1.0e-9

      dphi(0) = 0.0

      do i = 1, dof-1


            ! the relation does not hold for xi = 0,1 (resp ti=-1,1)
            if ( abs( xi ) <  eps ) then !xi == 0
               dphi(i) = (-1.0)**(i+1) * i * (i + 1.0)
            else if ( abs(xi - 1.0 ) < eps ) then !xi == 1
               dphi(i) = i * (i + 1.0)
            else !recursive relation for derivative of the Legendre polynomials
               dphi(i) = ( 2*i * ( ti*ti - 1.0 )**(-1.0) ) &
                  * ( ti*phi(i) - phi(i-1)  )
            end if

      end do

   end if


   !normalize on (0,1)
   do i = 0, dof-1
      phi(i) = phi(i) * sqrt(2.0*i + 1.0)
   end do

   if ( present( dphi ) ) then
      do i = 0, dof-1
         dphi(i) = dphi(i) * sqrt(2.0*i + 1.0)
      end do
   end if

   ! deallocate(ti)

  end subroutine Eval_LegendrePolynomialsOneNode

     !> evaluation of the Legenre polynomials and its derivatives in the given set of nodes
  !> eval values of ON Legendre polynomials on (0,1) in the given set of nodes in [0,1]
  subroutine SetLegendrePolynomials( dof, Qdof, xi, phi, dphi)
    !class(Gauss_rule), intent(inout) :: this
    integer, intent(in) :: dof ! = number of legendre functions
    integer, intent(in) :: Qdof
    real, dimension(1:Qdof), intent(in) :: xi ! sought nodes
    real, dimension(0:dof-1, 1:Qdof), intent(out) :: phi ! Gauss integ nodes
    real, dimension(0:dof-1, 1:Qdof), intent(out), optional :: dphi ! Gauss integ nodes

    real, dimension(:), allocatable :: ti ! Gauss integ nodes transformed to (-1,1)
    integer:: i,j
    real :: eps


   !print*, 'Qdof = ', Qdof

   allocate( ti(1:Qdof) )

   ti(1:Qdof) = 2*xi(1:Qdof) - 1. ! from 0,1 to -1,1

   phi(0, 1:Qdof) = 1.
   phi(1, 1:Qdof) = ti(1:Qdof)

   !write(*,'(a6,60es12.4)') 'ti: ',ti(1:Qdof)

   do i=1,dof-2

    phi(i+1, 1:Qdof) = (2.0*i + 1.0) * ti(1:Qdof) * phi(i, 1:Qdof) &
         - ( 1. * i * phi(i-1, 1:Qdof) )

    phi(i+1, 1:Qdof) =  phi(i+1, 1:Qdof) / (i+1)
   enddo



   if ( present( dphi ) ) then

      eps = 1.0e-9

      dphi(0, 1:Qdof) = 0.0
   !      dphi(1, 1:Qdof) = 2.0
      do i = 1, dof-1
         do j = 1,Qdof

            ! the relation does not hold for xi = 0,1 (resp ti=-1,1)
            if ( abs(xi(j)) <  eps ) then !xi(j) == 0
               dphi(i,j) = (-1.0)**(i+1) * i * (i + 1.0)
            else if ( abs(xi(j) - 1.0 ) < eps ) then !xi(j) == 1
               dphi(i,j) = i * (i + 1.0)
            else !recursive relation for derivative of the Legendre polynomials
               dphi(i,j) = ( 2*i * ( ti(j)*ti(j) - 1.0 )**(-1.0) ) &
                  * ( ti(j)*phi(i,j) - phi(i-1,j)  )
            endif

         end do

      end do

   end if


   !normalize on (0,1)
   do i = 0, dof-1
      phi(i, 1:Qdof) = phi(i, 1:Qdof) * sqrt(2.0*i + 1.0)
   end do

   if ( present( dphi ) ) then
      do i = 0, dof-1
         dphi(i, 1:Qdof) = dphi(i, 1:Qdof) * sqrt(2.0*i + 1.0)
      end do
   end if

   ! deallocate(ti)

  end subroutine SetLegendrePolynomials

   ! computed function is given by coefs with respect to the Legendre time basis functions
   ! computes function value in time=ti \in [0,1]
  pure function evalTimeFunctionValue( dof, f, ti) result ( funValue )
   integer, intent(in) :: dof ! # of basis functions
   real, intent(in) :: ti ! time in the refference interval
   real, dimension(1:dof), intent(in) :: f !coefs of the function
   real :: funValue
   real, dimension(1:dof) :: phi

   call Eval_LegendrePolynomialsOneNode( dof, ti, phi(1:dof) )

   funValue = dot_product( f(1:dof), phi(1:dof) )

  end function evalTimeFunctionValue

   ! computed derivative of a function given by coefs with respect to the Legendre time basis functions
   ! computes function derivative value in time=ti \in [0,1]
  pure function evalTimeFunctionDerivativeValue( dof, f, ti) result ( funValue )
   integer, intent(in) :: dof ! # of basis functions
   real, intent(in) :: ti ! time in the refference interval
   real, dimension(1:dof), intent(in) :: f !coefs of the function
   real :: funValue
   real, dimension(1:dof) :: phi, dphi

   call Eval_LegendrePolynomialsOneNode( dof, ti, phi(1:dof), dphi(1:dof) )

   funValue = dot_product( f(1:dof), dphi(1:dof) )

  end function evalTimeFunctionDerivativeValue




end module tquadrature_mod
