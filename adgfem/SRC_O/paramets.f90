!> global parameters
module paramets
  implicit none
  public

  !  used in  ElementViscous..., number of arrays, 1= Re, 2,.. coeffs, iRe = 1 + #coeffs !!!!!
  integer, parameter  :: iRe = 4

  integer :: nbDim  ! number of the problem space dimension, read from ini file
  integer :: ndim   ! number of the equations, the same state%ndim
  integer :: Bdim   ! number of block types for saddle-point (SP) problems

  !> saddle point problems, indexes of blocks
  integer, parameter :: bVV = 1      ! array STblock(1, : ) velocity-velocity
  integer, parameter :: bVP = 2      ! array STblock(2, : ) velocity-pressure
  integer, parameter :: bPV = 3      ! array STblock(3, : ) pressure-velocity

!> saddle point problems, indexes of unknowns type
!  integer, parameter :: wV1 = 1      ! array elem*wSP(1, : , :) - 1st component of velocity
!  integer, parameter :: wV2 = 2      ! array elem*wSP(2, : , :) - 2st component of velocity
!  integer, parameter :: wP = 3       ! array elem*wSP(3, : , :) - pressure

  !> global variables for arrays of elem%face
  integer, parameter :: idx = 1      ! 1st index of elem%face(*, : ) = vertices
  integer, parameter :: neigh = 2    ! 2nd index of elem%face(*, : ) = neighbours
  integer, parameter :: nei_i = 3    ! 3rd index of elem%face(*, : ) = index of neigh, i.e. if j = elem(i)%face(neigh,1) and k = elem(i)%face(nei_i,1) then i = elem(j)%face(neigh, k) - (at which position is elem(i) in local structure of elem(j)
  integer, parameter :: fGnum =  4   ! 4th index of elem%face(*, : ) = index of Gauss quadr.
  integer, parameter :: fGdof =  5   ! 5th index of elem%face(*, : ) = Gdof of Gauss quadr.
  integer, parameter :: fdeg  = 6    ! 6th index of elem%face(*, : ) = max(elem%deg, elem1%deg)
  integer, parameter :: fdof  = 7    ! 7th index of elem%face(*, : ) = elem1%dof
  integer, parameter :: rot  = 8     ! 8th index of elem%face(*, : ) = rotation of face for 3D comp.
  integer, parameter :: nbnode  = 9  ! 9th index of elem%face(*, : ) = number of vertices of face
  integer, parameter :: fTdeg  = 10  !10th index of elem%face(*, : ) = time degree on neighbouring element
  integer, parameter :: fTdof  = 11  !11th index of elem%face(*, : ) = DOF on  neighbouring element
  integer, parameter :: fdof_plus = 12    ! 12th index of elem%face(*, : ) = elem1%dof_plus
!  integer, parameter :: fdegP  = 12    ! 6th index of elem%face(*, : ) = max(elem%deg, elem1%degP)
!  integer, parameter :: fdofP  = 13    ! 7th index of elem%face(*, : ) = elem1%dofP

  integer, parameter :: max_face = 12  ! maximum of the previous indexes (only for allocation)

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
  integer, parameter :: resSr = 5
  integer, parameter :: min_resT_S = 6
  integer, parameter :: max_resT_S = 7
  integer, parameter :: resA_ST = 8
  integer, parameter :: min_resT_S_loc = 9
  integer, parameter :: resA_ST_loc = 10
  integer, parameter :: res_HO_p0 = 11
  integer, parameter :: res_HO_p1 = 12
  integer, parameter :: res_HO_p2 = 13
  !!integer, parameter :: estim_loc = 20

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
  integer, parameter :: P_HG = 7
  integer, parameter :: P_sF = 8
  integer, parameter :: P_su = 9
  integer, parameter :: P_FDu = 10
  integer, parameter :: P_F_p1 = 11
  integer, parameter :: P_F_p2 = 12
  integer, parameter :: P_F_p3 = 13
  integer, parameter :: P_F_p4 = 14
  integer, parameter :: P_s_p1 = 15
  integer, parameter :: P_s_p2 = 16
  integer, parameter :: P_s_p3 = 17
  integer, parameter :: P_s_p4 = 18
  integer, parameter :: P_potP = 19
  integer, parameter :: P_potPP = 20

  !> global parameters for arrays of elem%eta for aposteriori DWR
  integer, parameter :: dwrA = 1
  integer, parameter :: dwrS = 2
  integer, parameter :: dwrS_abs = 3
  integer, parameter :: dwrE = 4  ! Exact error
  integer, parameter :: dwrST = 5 ! not used
  integer, parameter :: dwrT = 6 !not used
  integer, parameter :: dwr_dualA = 7
  integer, parameter :: dwr_dualS = 8
  integer, parameter :: dwr_dualS_abs = 9
  integer, parameter :: dwr_Juh = 10
  integer, parameter :: dwr_aver = 11
  integer, parameter :: dwr_aver_abs = 12
  integer, parameter :: dwr_etaS = 13
  integer, parameter :: dwr_sign = 14
  integer, parameter :: dwr_dual_sign = 15
  integer, parameter :: dwrLinP = 16
  integer, parameter :: dwrLinD = 17

  integer, parameter :: dwr_max_eta = 17


  integer, parameter :: max_eta = 20  ! maximum of the RTN and HEL indexes (only for allocation)


  !> global parameters for DWR target functionals / Quantity of Interest
!  integer, parameter :: DWR_pointVal = 1


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
  integer, parameter :: Ejumps = 18       ! sum_\Gamma [ u_h]^2
  integer, parameter :: EjumpsG = 19       ! sum_\Gamma [ u_h]^2/ h_G
  integer, parameter :: EjumpsPi = 20       ! sum_\Gamma [\Pi^0 u_h]^2/ h_G
  integer, parameter :: CFD_cD = 15       ! c_D
  integer, parameter :: CFD_cL = 16       ! c_L
  integer, parameter :: CFD_cM = 17       ! c_M

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
   integer, parameter :: maxLrule =  22 ! MaxDegreeImplemented ! FR ???

   !Multigrid parameters
   integer, parameter :: MaxMGlevel = 9

   ! parameter for quadrilaterals quadrature indexes
   integer :: QnumOffset = maxVrule

   integer :: max_RGlevel = 30       ! maximal implemented level of RD

   !DEBUG file paramets
   character(len=50) :: debug_file = 'debug' ! '../DWR/debug'
   integer :: debug = 68


 contains

   !> compute the order of the derivative \f$ \partial_x^i \partial_y^j \f$ in the sequence:
   !> \f$ \cdot, \partial_x, \partial_y, \partial_{xx}, \partial_{xy}, \partial_{yy}, \dots \f$
   function Deriv_idx(i, j)
     integer :: Deriv_idx
     integer :: i, j

     Deriv_idx = (i+j + 1) * (i+j+2) /2 - i

   end function Deriv_idx

   !> compute the indexes i and j from the given order of the derivative
   !> \f$ \partial_x^i \partial_y^i \f$ in the sequence:
   !> \f$ \cdot, \partial_x, \partial_y, \partial_{xx}, \partial_{xy}, \partial_{yy}, \dots \f$
   subroutine Deriv_idx_reverse(iphi, i, j)
     integer, intent(in) :: iphi
     integer, intent(inout) :: i, j
!     integer, dimension (:,:), allocatable :: loc_idx
     integer ::  l
     real :: r

     ! allocate(loc_idx(1:45, 1:2) )

     ! loc_idx( 1, 1:2)  = (/ 0, 0 /)
     ! loc_idx( 2, 1:2)  = (/ 1, 0 /)
     ! loc_idx( 3, 1:2)  = (/ 0, 1 /)

     ! loc_idx( 4, 1:2)  = (/ 2, 0 /)
     ! loc_idx( 5, 1:2)  = (/ 1, 1 /)
     ! loc_idx( 6, 1:2)  = (/ 0, 2 /)

     ! loc_idx( 7, 1:2)  = (/ 3, 0 /)
     ! loc_idx( 8, 1:2)  = (/ 2, 1 /)
     ! loc_idx( 9, 1:2)  = (/ 1, 2 /)
     ! loc_idx(10, 1:2)  = (/ 0, 3 /)

     ! loc_idx(11, 1:2)  = (/ 4, 0 /)
     ! loc_idx(12, 1:2)  = (/ 3, 1 /)
     ! loc_idx(13, 1:2)  = (/ 2, 2 /)
     ! loc_idx(14, 1:2)  = (/ 1, 3 /)
     ! loc_idx(15, 1:2)  = (/ 0, 4 /)

     ! loc_idx(16, 1:2)  = (/ 5, 0 /)
     ! loc_idx(17, 1:2)  = (/ 3, 1 /)
     ! loc_idx(18, 1:2)  = (/ 3, 2 /)
     ! loc_idx(19, 1:2)  = (/ 2, 3 /)
     ! loc_idx(20, 1:2)  = (/ 1, 4 /)
     ! loc_idx(21, 1:2)  = (/ 0, 5 /)

     ! if(iphi > 21) then
     !    stop "insufficient array in paramets.f90 94ur943ju"
     ! endif
     ! !i1 = loc_idx(iphi, 1)
     ! !j1 = loc_idx(iphi, 2)

     r = ( (sqrt(1. + 8*iphi) - 1) / 2 + 0.00001) ! in order to avoid rounding errors

     ! rounding of r to the higher integer
     l = int(r)
     if(l <= r - 0.0001) l = l + 1

     i = l * (l+1) / 2 - iphi
     j = l  - i - 1

     !print*, '73d3hd', iphi,  loc_idx(iphi, 1:2), '||',i , j , '|', r, l
     !print*, '73d3hd', iphi, '||',i , j , '|', r, l

     !deallocate(loc_idx)

   end subroutine Deriv_idx_reverse

end module paramets
