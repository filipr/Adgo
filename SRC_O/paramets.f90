!> global parameters
module paramets
  implicit none
  public


  integer :: nbDim  ! number of the problem space dimension, read from ini file
  integer :: ndim   ! number of the equations, the same state%ndim
  integer :: Bdim   ! number of block types for saddle-point (SP) problems

  !> saddle point problems, indexes of blocks
  integer, parameter :: bVV = 1      ! array STblock(1, : ) velocity-velocity
  integer, parameter :: bVP = 2      ! array STblock(2, : ) velocity-pressure
  integer, parameter :: bPV = 3      ! array STblock(3, : ) pressure-velocity

!> saddle point problems, indexes of unknowns type
  integer, parameter :: wV1 = 1      ! array elem*wSP(1, : , :) - 1st component of velocity
  integer, parameter :: wV2 = 2      ! array elem*wSP(2, : , :) - 2st component of velocity
  integer, parameter :: wP = 3       ! array elem*wSP(3, : , :) - pressure

  !> global variables for arrays of elem%face
  integer, parameter :: idx = 1      ! 1st index of elem%face(*, : ) = vertices
  integer, parameter :: neigh = 2    ! 2nd index of elem%face(*, : ) = neighbours
  integer, parameter :: nei_i = 3    ! 3rd index of elem%face(*, : ) = index of neigh
  integer, parameter :: fGnum =  4   ! 4th index of elem%face(*, : ) = index of Gauss quadr.
  integer, parameter :: fGdof =  5   ! 5th index of elem%face(*, : ) = Gdof of Gauss quadr.
  integer, parameter :: fdeg  = 6    ! 6th index of elem%face(*, : ) = max(elem%deg, elem1%deg)
  integer, parameter :: fdof  = 7    ! 7th index of elem%face(*, : ) = elem1%dof
  integer, parameter :: rot  = 8     ! 8th index of elem%face(*, : ) = rotation of face for 3D comp.
  integer, parameter :: nbnode  = 9  ! 9th index of elem%face(*, : ) = number of vertices of face
  integer, parameter :: fTdeg  = 10  !10th index of elem%face(*, : ) = time degree on neighbouring element
  integer, parameter :: fTdof  = 11  !11th index of elem%face(*, : ) = DOF on  neighbouring element
  integer, parameter :: fdegP  = 12    ! 6th index of elem%face(*, : ) = max(elem%deg, elem1%degP)
  integer, parameter :: fdofP  = 13    ! 7th index of elem%face(*, : ) = elem1%dofP

  integer, parameter :: max_face = 13  ! maximum of the previous indexes (only for allocation)

  !> global variables for arrays of elem%vec
  integer, parameter :: rhs = 1      ! 1st index of array elem%vec(*, : )  = RHS
  integer, parameter :: rhsM = 2     ! 2nd index of array elem%vec(*, : )  = mass RHS
  integer, parameter :: res_vec = 3  ! 3rd index of array elem%vec(*, : )  = residual vector
  integer, parameter :: aRe  = 4     ! 4th index of array elem%vec(*, : )  = artificial viscosity
  integer, parameter :: rhsT = 5     ! 5th index of array elem%vec(*, : )  = mass RHS for time estim
  integer, parameter :: rhsOLD = 6   ! 6th index of array elem%vec(*, : )  = RHS at k-1 level
  integer, parameter :: res_func = 7 ! 7th index of array elem%vec(*, : )  = coefficients of residual function (made up from residual vector) in basis

  integer, parameter :: max_vecs = 7  ! maximum of the previous indexes (only for allocation)

  !> global variables for arrays of elem%eta for aposteriori RTN error estimates
  integer, parameter :: DFn = 1
  integer, parameter :: DFnS = 2
  integer, parameter :: DFnT = 3
  integer, parameter :: Rn = 4
  integer, parameter :: DFRn = 5
  integer, parameter :: NC1n = 6
  integer, parameter :: NC2n = 7
  integer, parameter :: Osc = 8
  integer, parameter :: total = 9
  integer :: IC = 7  ! only for DUA !!FR problem - it is used in ama-angener as variable
  integer, parameter :: eN1 = 10
  integer, parameter :: eN2 = 11
  integer, parameter :: eN1p = 12
  integer, parameter :: eN2p = 13
  integer, parameter :: eN3p = 14
  integer, parameter :: quadra = 15

  !> global variables for arrays of elem%eta for a posteriori RTNst error estimation
  integer, parameter :: RTNall = 1 !TODO comment FR ! not used
  integer, parameter :: RTNeta = 2
  integer, parameter :: RTNrez = 3
  integer, parameter :: RTNflux = 4
  integer, parameter :: RTNradau = 5
  integer, parameter :: RTNradau2 = 6
  integer, parameter :: RTNjump = 7
  integer, parameter :: RTNfluxnorm = 8

  !> global variables for arrays of elem%eta for aposteriori HELmholtz error estimates
  integer, parameter :: Hrez = 1
  integer, parameter :: Hjump = 2
  integer, parameter :: HjumpD = 3
  integer, parameter :: Hjump_1 = 4
  integer, parameter :: Heta1 = 5
  integer, parameter :: Heta = 6

  !> global variables for arrays of elem%eta for RTN-based a posteriori ALGebraic error estimates
  integer, parameter :: FNCD = 1  !Discretization component of Flux NonConformity estimator
  integer, parameter :: FNCA = 2  !Algebraic component of Flux NonConformity estimator
  integer, parameter :: FD = 3    !Flux Discretization estimator: \|\nabla u_h + d_h\|_{L^2}
  integer, parameter :: FA = 4    !Flux Algebraic estimator: \|a_h\|_{L^2}
  integer, parameter :: Resid = 5 !Residual estimator
  integer, parameter :: PNC = 6   !Potential NonConformity estimator
  integer, parameter :: Disc = 7  !Discretization error estimator
  integer, parameter :: Alg = 8   !Algebraic error estimator
  integer, parameter :: Rem = 9   !algebraic Reminder estimator
  integer, parameter :: Tot = 10  !guaranteed upper bound
  integer, parameter :: FT = 11   !Total flux estimator: \|\nabla u_h + t_h\|_{L^2}
  integer, parameter :: FNCT = 12 !Total Flux NonConformity estimator

  !> global variables for arrays of elem%eta for aposteriori RESidual error estimates
  integer, parameter :: resA = 1
  integer, parameter :: resS = 2
  integer, parameter :: resT = 3
  integer, parameter :: resST = 4
  integer, parameter :: min_resT_S = 5
  integer, parameter :: max_resT_S = 6
  integer, parameter :: resA_ST = 7
  integer, parameter :: min_resT_S_loc = 8
  integer, parameter :: resA_ST_loc = 9

  !> global variables for arrays of elem%eta for aposteriori HO_rec
  integer, parameter :: HO_estim_L2_p0 = 1  !
  integer, parameter :: HO_estim_H1_p0 = 2  !
  integer, parameter :: HO_trunc_L2_p0 = 3  !
  integer, parameter :: HO_trunc_H1_p0 = 4  !

  integer, parameter :: HO_estim_L2_p1 = 5
  integer, parameter :: HO_estim_H1_p1 = 6
  integer, parameter :: HO_trunc_L2_p1 = 7
  integer, parameter :: HO_trunc_H1_p1 = 8

  integer, parameter :: HO_estim_L2_p2 = 9
  integer, parameter :: HO_estim_H1_p2 = 10
  integer, parameter :: HO_trunc_L2_p2 = 11
  integer, parameter :: HO_trunc_H1_p2 = 12
  integer, parameter :: HO_recovery    = 13
  integer, parameter :: HO_rec_estim   = 14

  !integer, parameter :: HO_hp_H0_Pm = 13
  !integer, parameter :: HO_hp_H0_P0 = 14
  !integer, parameter :: HO_hp_H0_Pp = 15

  !integer, parameter :: HO_hp_Hp_Pm = 16
  !integer, parameter :: HO_hp_Hp_P0 = 17
  !integer, parameter :: HO_hp_Hp_Pp = 18


  !> global variables for arrays of elem%eta for aposteriori pNeu
  integer, parameter :: P_flux = 1
  integer, parameter :: P_rez = 2
  integer, parameter :: P_FR = 3
  integer, parameter :: P_pot = 4
  integer, parameter :: P_BC = 5
  integer, parameter :: P_tot = 6
  integer, parameter :: P_sF = 7
  integer, parameter :: P_su = 8
  integer, parameter :: P_FDu = 9
  integer, parameter :: P_F_p1 = 10
  integer, parameter :: P_F_p2 = 11
  integer, parameter :: P_F_p3 = 12
  integer, parameter :: P_F_p4 = 13
  integer, parameter :: P_s_p1 = 14
  integer, parameter :: P_s_p2 = 15
  integer, parameter :: P_s_p3 = 16
  integer, parameter :: P_s_p4 = 17

  !FILIP: is it OK ?
  integer, parameter :: P_potP = 18


  integer, parameter :: max_eta = 20  ! maximum of the RTN and HEL indexes (only for allocation)


  !> global parameters for DWR target functionals / Quantity of Interest
  integer, parameter :: DWR_pointVal = 1


  !> global variables for arrays of elem%MGw
  integer, parameter :: MGv = 1
  integer, parameter :: MGr = 2
  integer, parameter :: max_MG = 2  ! maximum of the previous indexes (only for allocation)

  !> global variables for RTN reconstruction-based adaptive solution algorithm
  real :: gamma_rem
  real :: gamma_alg
  real :: gamma_lin
  integer :: nu
  character :: stop_crit

  !> global variables for errSTnorm estimates
  integer, parameter :: L8L2 = 1 ! L^\infty(0,T,L^2(\Omega)) norm (computed in Gauss/Radau integration nodes)
  integer, parameter :: L2L2 = 2 ! L^2(0,T,L^2(\Omega)) norm
  integer, parameter :: L2H1 = 3 ! L^2(0,T,H^2(\Omega)) norm
  integer, parameter :: L2H10 = 4 !L^2(0,T,H^2(\Omega)) semi-norm
 ! integer, parameter :: L2DG = 5 ! L^2(0,T,DG) norm
  integer, parameter :: L2L2eH1 = 5 ! L^2(0,T,L^2 + e H^1) norm
  !integer, parameter :: L2X_sc = 6 ! used for algeb. stop criteria AEE
  integer, parameter :: L2F = 6 ! L^2(0,T, F) norm type of dual norm
  integer, parameter :: NC = 7 !
  integer, parameter :: Snorm1 = 8
  integer, parameter :: Snorm2 = 9
  integer, parameter :: Snorm3 = 10
  integer, parameter :: max_errSTnorm = 10 ! maximum of the previous indices

  !> global variables for state%err
  integer, parameter :: L2 = 1 ! L^\infty(0,T,L^2(\Omega)) norm (computed in Gauss/Radau integration nodes)
  integer, parameter :: H1 = 2 ! L^2(0,T,L^2(\Omega)) norm
  integer, parameter :: L2_old = 3 ! L^2(0,T,H^2(\Omega)) norm
  integer, parameter :: H1_old = 4 !L^2(0,T,H^2(\Omega)) semi-norm
  integer, parameter :: H1_discrete = 5 !  H^1 discrete in pNeu for SIPG and NIPG [Ern, Voh, SINUM 15]
  !!integer, parameter :: SSL2 = 5 ! L^2(0,T,L^2 + e H^1) norm
  integer, parameter :: SSL8 = 6 ! used for algeb. stop criteria AEE
  integer, parameter :: SSnew = 7
  !integer, parameter :: SSnew8 = 8
  integer, parameter :: IC_L2 = 8 ! initial condition
  integer, parameter :: IC_H1 = 9 ! initial condition
  integer, parameter :: Terr = 10
  integer, parameter :: Terr_loc = 11
  integer, parameter :: interLq = 12
  integer, parameter :: interL8 = 13
  integer, parameter :: interH1 = 14  ! interH1 and algeb ARE USED SIMULTANEOUSLY !!!!!!
  integer, parameter :: algeb = 14
  integer, parameter :: err_0 = 15
  integer, parameter :: L8 = 16       ! L^\infty(\Omega))
  integer, parameter :: XX = 17       ! X norm
  integer, parameter :: Ejumps = 18       ! sum_\Gamma [ u_h]
  integer, parameter :: EjumpsG = 19       ! sum_\Gamma [ u_h]/ h_G
  integer, parameter :: EjumpsPi = 20       ! sum_\Gamma [\Pi^0 u_h]/ h_G
  integer, parameter :: max_errS = 20 ! maximum of the previous indices

  !parameters moved from state -> what is implemented

   ! minimal implemented degree of SPACE polynomial approximation
   integer, parameter :: MinDegreeImplemented = 1
   ! maximal implemented degree of SPACE polynomial approximation - connected with quadratures and basis functions
   integer, parameter :: MaxDegreeImplemented = 10

   !  time basis functions computed recurrently from Legendre pol -> only restr is the maxRadau degree = 12
   integer, parameter :: MaxTimeDegree = 10
   !integer, parameter :: MaxTimeDegree = 4            !  maximal implemented time-degree ???FERROR is it connected with BDF or STDG
   integer, parameter :: MaxRTNImplemented = MaxDegreeImplemented       ! maximal implemented RTN

   !maximal degrees of implemented quadratures
   integer, parameter :: maxVrule = 23 ! max volume_rule
   integer, parameter :: maxGrule = 15
   integer :: maxTrule = 12 ! maxRrule or maxGrule - depends on which quad is chosen in initTime
   integer, parameter :: maxRrule = 12 !Radau quadrature rule
   integer, parameter :: maxLrule = MaxDegreeImplemented

   !Multigrid parameters
   integer, parameter :: MaxMGlevel = 9

   ! parameter for quadrilaterals quadrature indexes
   integer :: QnumOffset = maxVrule

   integer :: max_RGlevel = 8              ! maximal implemented level of RD

   !DEBUG file paramets
   character(len=50) :: debug_file = 'debug'
   integer :: debug = 68

end module paramets
