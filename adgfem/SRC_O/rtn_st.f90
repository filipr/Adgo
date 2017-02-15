module rtn_st_mod

   use blocks_integ
   use data_mod
   use define_state
   use dual_estim
   use element_mod
   use eval_rav_tho_ned
   use eval_rtn_st_mod
   use eval_sol
   use loc_rav_tho_ned
   use main_data
   use paramets
   use stdgm_mod




implicit none

   public :: ComputeRTNstEstim
   public :: ComputeElemRTNstEstim
   public :: ComputeRTNMomenta
   public :: exportRTNst
   public :: writeInitDualProblemFile
   public :: setDomainCorners
   public :: SetRTNstRhs
   public :: setWeightingConstants




contains

!> estimates the computational error using the approach based on the flux reconstruction
   subroutine ComputeRTNstEstim( estim , tQnum )
      real, dimension(1:max_eta,1:ndim), intent(out) :: estim
      class(element), pointer :: elem
      integer, intent(in) :: tQnum
      integer :: i !, spaceLevels, timeLevels
      integer :: Fdeg, tdeg
      class( Time_rule), pointer :: T_rule
      logical :: impl

      ! for nonlinear problems !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! we need to call ComputeSTDGM_Term with implicitly = true
      ! to compute the constants CK, Cb used for weighting the estimates
      impl = state%nlSolver%implicitly
      state%nlSolver%implicitly  = .true.
      call ComputeSTDGM_Terms( .false. )
      state%nlSolver%implicitly  = impl


      write( debug,* ) 'ComputeRTNstEstim - We expect that the degree',&
         'of time polynomial approximation is the same for all triangles! '

      if(state%time%iter == 1) &
         print*, 'IC estimator should be implemented!'

      estim(1:max_eta,1:ndim) = 0.0

      !initialize loc_RTN for further use in ComputeElemRTN_STEstim
      do Fdeg = 1, state%space%deg
         if ( .not. state%loc_RTN(Fdeg)%defined ) &
            call Init_Loc_RTN( state%loc_RTN(Fdeg) , Fdeg )
      enddo

      !!! Compute Radau polynomial !!!

      ! degree of the polynomial Radau reconstruction
      tdeg =  state%time%deg + 1

      T_rule => state%time%T_rule( tdeg )
      select type ( T_rule )
         type is ( RadauTime_rule )
            !print*, 'time%deg=', state%time%deg, ' # nodes: ',T_rule%Qdeg

            if ( .not. allocated( T_rule%RadauPol ) ) &
               call T_rule%evalRadauPolynomial()
            !call T_rule%testRadauPol()

         type is ( GaussTime_rule )
            stop 'Radau polynomial cannot be computed with Gauss quadrature'

         class default
            stop ' Abstract Time_rule - should be specified (Radau or Gauss) '
      end select


      do i = 1, grid%nelem
         elem => grid%elem(i)
         call ComputeElemRTNstEstim( elem, tQnum )

         estim( RTNrez, 1:ndim ) = &
           estim( RTNrez, 1:ndim ) + ( elem%eta( RTNrez, 1:ndim ) )**(2.0)
         estim( RTNflux, 1:ndim ) = &
           estim( RTNflux, 1:ndim ) + ( elem%eta( RTNflux, 1:ndim ) )**(2.0)
         ! |(R-u)'|
         estim( RTNradau, 1:ndim ) = &
            estim( RTNradau, 1:ndim ) + ( elem%eta( RTNradau, 1:ndim ) )**(2.0)
         !New version |R-u| ... space time L2 norm
         estim( RTNradau2, 1:ndim ) = &
            estim( RTNradau2, 1:ndim ) + ( elem%eta( RTNradau2, 1:ndim ) )**(2.0)

         estim( RTNjump, 1:ndim ) = &
           estim( RTNjump, 1:ndim ) + elem%eta( RTNjump, 1:ndim )

         !OLD
!         estim( RTNeta, 1:ndim ) = estim( RTNeta, 1:ndim ) + &
!            ( elem%eta(RTNrez, 1:ndim) + elem%eta(RTNflux, 1:ndim) + elem%eta(RTNradau, 1:ndim) )**(2.0)
         !NEW Radau2 - ST L2 norm
         estim( RTNeta, 1:ndim ) = estim( RTNeta, 1:ndim ) + &
            ( elem%eta(RTNrez, 1:ndim) + elem%eta(RTNflux, 1:ndim) + elem%eta(RTNradau2, 1:ndim) )**(2.0)

          ! new flux error
          estim( RTNfluxnorm, 1:ndim ) = &
            estim( RTNfluxnorm, 1:ndim ) + ( elem%eta( RTNfluxnorm, 1:ndim ) )**(2.0)

      enddo ! i

      print*, 'Eta_rez in the ', state%time%iter,'-th time step:', sqrt( estim( RTNrez, 1:ndim ) )
      print*, 'Eta_flux in the ', state%time%iter,'-th time step:', sqrt( estim( RTNflux, 1:ndim ) )
      print*, 'Eta_radau in the ', state%time%iter,'-th time step:', sqrt( estim( RTNradau, 1:ndim ) )
      print*, 'Eta_radau2 in the ', state%time%iter,'-th time step:', sqrt( estim( RTNradau2, 1:ndim ) )
      print*, 'Eta_jump in the ', state%time%iter,'-th time step:', sqrt( estim( RTNjump, 1:ndim ) ) !, '^2' ,estim( RTNjump, 1:ndim )
      print*, 'Flux norm in the ', state%time%iter,'-th time step:', sqrt( estim( RTNfluxnorm, 1:ndim ) )


     !print*, 'Max of d_K:' , maxval( grid%elem(:)%dK )

!      print*, 'Estimator in the ', state%time%iter,'-th time step:', sqrt( estim( RTNall, 1:ndim ) ) , &
!         '+', sqrt( estim( RTNjump, 1:ndim ) )

!      stop 'End of ComputeRTNstEstim!'

   end subroutine ComputeRTNstEstim


   !> compute error estimates via the dual (residual) form using RTN flux reconstruction for one element
   subroutine ComputeElemRTNstEstim( elem , tQnum )
      class(element), intent(inout) :: elem
      integer, intent(in) :: tQnum !number of time moments
      real, dimension(:,:), allocatable :: Mspace ! inverse of Momentum matrix on elem
      real, dimension(:,:,:), allocatable :: flux   ! basis coefficients of flux (1:ndim - future), 1:tdeg, 1:Fdof
      real, dimension(:,:,:), allocatable :: RadauReconstruct
      real, dimension(:,:,:), allocatable :: rhs !, space_rhs
      integer :: Fdeg, Fdof
      integer :: Qdof, Tdof
      real :: dK, cTn, cKn, Cnc, Ckb, h
      real, dimension(:,:,:), allocatable :: fun
      character(len=20) :: outputfile
      real :: errL2, errH1
      real, dimension(:), allocatable :: errL8

      allocate( errL8(1:elem%TQnum ) )

      Tdof = elem%Tdof

      if (.not. allocated(elem%eta) ) stop 'elem%eta is not allocated in ComputeElemRTNEstim'
      elem%eta(1:max_eta, 1:ndim) = 0.0

      Fdeg = elem%deg
      Fdof = SetRTNdof(Fdeg)
      Qdof = elem%Qdof

      !loc_RTN => state%loc_RTN(Fdeg)
      !TODO it should be defined!
      if(.not. state%loc_RTN(Fdeg)%defined ) &
         stop 'Problem in ComputeElemRTNEstim: loc%RTN is not defined!'

      allocate(Mspace(1:Fdof, 1:Fdof) )


      ! Same as ComputeRTNMomenta for HG nodes too
      ! call ComputeLocRTNMomentumsElem( elem, Fdeg, Fdof, Mspace )

      ! evaluate the momentums of the local RTN basis on K
      Mspace(1:Fdof, 1:Fdof) = ComputeRTNMomenta( elem, Fdeg, Fdof )
!
!      write(outputfile, "(A6,I0)") 'Mspace' ,1 !, tQnum
!
!      call WriteLocalMatrix( Mspace, outputfile)


!!!!!!! SETTING SPACE-TIME RHS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !SIPG,NIPG - additional boundary member of RHS for the volume part

      allocate( rhs( 1:ndim, 1:Fdof, 1:Tdof ) )
      ! rhs(m,i) = coefficient of t with respect to basis $\phi_i(t) \varphi(x)_j$
      rhs(1:ndim,1:Fdof, 1:Tdof) = SetRTNstRhs( elem, Fdeg, Fdof, tQNum)
!
!      write(outputfile, "(A6,I0)") 'rhs' ,1 !, tQnum
!
      !call WriteLocalMatrix( rhs(1,:,:) , 'rhs')



!!!!!!COMPUTE THE SPACE-TIME FLUX RECONSTRUCTION $t_{h,\tau}$
      allocate(flux(1:ndim, 1:Fdof,1:Tdof) )
      ! Set the basis coefficients of the RTN flux reconstruction with respect to the basis \hat\phi (made from standard basis functions)
      flux(1:ndim,1:Fdof, 1:Tdof) = SetRTNstFluxReconstr( Tdof, Fdof, &
         Mspace(1:Fdof, 1:Fdof), rhs(1:ndim, 1:Fdof, 1:Tdof) )

!!!!!!! compute the estimates !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elem%eta( RTNrez, 1:ndim) = 0.0
      elem%eta( RTNflux, 1:ndim) = 0.0
      elem%eta( RTNradau, 1:ndim) = 0.0 !evalRadauEta( elem, Rad
      elem%eta( RTNradau2, 1:ndim) = 0.0
      elem%eta(RTNfluxnorm, 1:ndim) = 0.0

!      if(state%time%iter == 1) then
!
!         print*, 'IC estimator should be implemented!'
!         !from aposter.f90
!         !call IC_Estimator(elem, Rdof, potential(0, 1:ndim, 1:Rdof), etaIC)
!         !state%err(IC_L2) = state%err(IC_L2) + sum(etaIC(:) )
!      endif


   !!!!!!!EVAL eta Flux = \f$ \norm{\sigma(u,\nabla u ) + \bkt} \f$ !!!!!!!!!!!!!!!
      elem%eta( RTNflux, 1:ndim) = RTNstFluxEstimator( elem, Fdeg, Fdof, Tdof, tQnum, &
          flux(1:ndim,1:Fdof, 1:Tdof) )

   !!!!!!EVAL eta RTNradau!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! set the Radau space-time reconstruction
      allocate ( RadauReconstruct( 1:ndim, 1:elem%dof, 1:elem%Tdof+1 ) )
      RadauReconstruct( 1:ndim, 1:elem%dof, 1:elem%Tdof+1 ) = evalRadauReconstruct( elem )
      allocate( fun(1:ndim,1:elem%dof,1:Tdof+1) , source = 0.0 )

      fun(1:ndim,1:elem%dof,1:Tdof) = RadauReconstruct( 1:ndim, 1:elem%dof, 1:Tdof ) - &
         elem%wST(1:ndim,1:elem%dof,1:Tdof)
      fun(1:ndim,1:elem%dof,Tdof+1) = RadauReconstruct( 1:ndim, 1:elem%dof, Tdof+1)

      elem%eta( RTNradau, 1:ndim ) = &
         evalH1L2STNorm_Elem(elem, state%time%tau(1), ndim, elem%dof, Tdof+1, fun )

      elem%eta( RTNradau2, 1:ndim ) = &
         evalL2STNorm_Elem( elem, state%time%tau(1), ndim, elem%dof, Tdof+1, fun )


      deallocate(fun)

   !!!!!!! EVAL eta RTNrez !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elem%eta( RTNrez, 1:ndim ) = RTNstRezEstimator( elem, elem%dof, Fdeg, Fdof, Tdof, tQnum, &
         flux(1:ndim, 1:Fdof, 1:Tdof), RadauReconstruct(1:ndim, 1:elem%dof, 1:Tdof+1) )

   !!!!!!! EVAL eta Jump   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elem%eta( RTNjump, 1:ndim ) = RTNstJumpEstimator( elem, tQnum )

   !!!!!!! EVAL eta RTNfluxnorm
      !e_FR = tau^-2 ||u-uh||^2 + h^-2 || \sigma(u) - \sigma(uh) ||^2
      elem%eta( RTNfluxnorm, 1:ndim ) = RTNstComputeFluxNorm( Set_R_s_scalar, Set_f_s_scalar, elem )


      !!!!!!! Weight the estimates !!!!!!!!!!!!!!!!!!!

      call setWeightingConstants( elem, dK, Cnc )
      elem%dK = dK

      elem%eta(RTNrez,1:ndim) = elem%CP * dK * sqrt(elem%eta(RTNrez,1:ndim))
      elem%eta(RTNflux,1:ndim) = ( dK / elem%diam ) * sqrt(elem%eta(RTNflux,1:ndim))
      !H1L2 norm
      elem%eta(RTNradau,1:ndim) = dK * sqrt(elem%eta(RTNradau,1:ndim))
      !ST L2 norm
      elem%eta(RTNradau2, 1:ndim) = ( dK / state%time%tau(1) ) * sqrt(elem%eta(RTNradau2, 1:ndim))

      !RTNjump is SQUARED !!!
      ! new C-D constant Cnc
      elem%eta(RTNjump, 1:ndim) = Cnc * elem%eta(RTNjump, 1:ndim)

      elem%eta(RTNfluxnorm, 1:ndim) = dK * sqrt( elem%eta(RTNfluxnorm, 1:ndim) )
      !print*, 'RTNfluxnorm is:' , elem%eta(RTNfluxnorm, 1:ndim)

      !      print*, 'ETA Rez:', elem%eta(RTNrez,:)
      !      print*, 'ETA Flux:', elem%eta(RTNflux,:)
      !      print*, 'ETA Radau:', elem%eta(RTNradau,:)
      !      print*, 'ETA Jump:', elem%eta(RTNjump,:)


      deallocate(Mspace, rhs, flux, RadauReconstruct)

   end subroutine ComputeElemRTNstEstim


  !> evaluate the momentums of the basis RTN function of degree Fdeg on elem
  !> and the resulting matrix is stored in MMRE (for nbdim=2 only)
  !> vector RTN basis functions made of scalar space basis functions on elem, (1,0)*phi,(0,1)*phi, x*phi
   function ComputeRTNMomenta(elem, Fdeg, Fdof) result ( MMRE )
      class(element), intent(in) :: elem
      integer, intent(in) :: Fdeg, Fdof

      real, dimension(1:Fdof, 1:Fdof) :: MMRE

!    type(basis_rtn_fe), pointer :: loc_RTN
      real, dimension(:,:,:), allocatable :: psi ! basis of RTN_K
      real, dimension(:,:), allocatable :: qi
      integer :: i, iedge, it, j
      real, dimension(:,:), allocatable :: x
      integer :: Qdof, Gnum, Gdof, Mdof
      real :: edge_momentum

      !real, dimension(:,:), allocatable :: help


      Mdof = maxval(elem%face(fGdof,:) )
      Mdof = max(Mdof, elem%Qdof)

      allocate( psi(1:Fdof, 1:nbDim, 1: Mdof)  )
      allocate(  qi(1:Mdof, 1:nbDim  ) )

      ! indices of vector RTN basis functions made of scalar space basis functions on elem, (1,0)*phi,(0,1)*phi, x*phi
      associate ( loc_RTN => state%loc_RTN(Fdeg) )


      i = 0   ! index of momentum

      MMRE(1:Fdof, 1:Fdof) = 0.0


      if(elem%HGnode) then
         stop 'ComputeRTNMomenta not implemented for HG nodes'
!       if(elem%type /= 3) print*,' TROUBLES in dual_estim.f90 with HG'
!
!       allocate( HGvertex(4) )
!       HGvertex(1:3)   = elem%HGvertex(1:3 )
!       HGvertex(4) = elem%flen + 1
      endif

      ! face momenta
      do iedge=1,3 ! loop over triagle edges
         ! index of G_rule quadrature for this edge
         Gnum = elem%face(fGnum,iedge)
         ! number of Quadrature nodes
         Gdof = elem%face(fGdof,iedge)
         associate( G_rule => state%space%G_rule(Gnum) )


         ! eval RTN basis functions on the edge iedge in integration nodes
         call Eval_Loc_RTN_Edge(elem, loc_RTN, iedge, psi(1:Fdof, 1:nbDim, 1:Gdof ) )

         do j=1, Fdeg+1  ! Fdeg+1 degrees of freedom over face
            i = i + 1

            call EvalRTNMomentumEdgeTestFunction( iedge, j, G_rule, qi(1:Gdof, 1) )

            !compute i-th row of the Momentum matrix MMRe
            do it = 1, Fdof

               call IntegrateFunctionNormalEdge(elem, iedge, &
                     psi(it, 1:nbDim, 1:Gdof),  qi(1:Gdof,1), edge_momentum )

                     MMRE(i, it) =  edge_momentum

            enddo !it
         enddo !j
         end associate !G_rule
      enddo !iedge

      write(debug,*) 'Control in ComputeRTNMomenta'
      !TODO remove CONTROL
      if ( i /= 3*(Fdeg + 1) ) then
         print*, 'Wrong index of i:', i, 'Fdeg=', Fdeg
      endif

      ! VOLUME MOMENTUMS
      Qdof = elem%Qdof

      ! volume integ nodes
      allocate(x(1:Qdof, 1:nbDim))
      associate ( V_rule => state%space%V_rule(elem%Qnum) )

      x(1:Qdof, 1:nbDim) = V_rule%lambda(1:Qdof,1:nbDim)
      ! eval RTN basis functions in volume integration nodes
      call Eval_Loc_RTN_Elem(elem, loc_RTN, psi(1:Fdof, 1:nbDim, 1:Qdof) )

      do j=1, Fdeg*(Fdeg+1)  ! Fdeg*(Fdeg+1) degrees of freedom over element
         i = i+1
         call EvalMomentumVolumeD( Fdeg, i, Qdof, x(1:Qdof, 1:nbDim), &
         qi(1:Qdof, 1:nbDim) )

         !compute i-th row of the Momentum matrix MMRe
         do it = 1, Fdof

            call IntegrateFunctionsVec(elem,  psi(it, 1:nbDim, 1:Qdof), &
            qi(1:Qdof, 1:nbDim),  MMRE(i, it) )

         enddo !it
      end do !j

      end associate ! loc_RTN

      deallocate(x, qi, psi)
      end associate ! loc_RTN


   end function ComputeRTNMomenta



   function SetRTNstRhs( elem, Fdeg, Fdof, tQnum) result(rhs)
      class(element), intent(inout) :: elem
      integer, intent(in) :: Fdeg
      integer, intent(in) :: Fdof
      integer, intent(in) :: tQnum !number of quadrature nodes

      real, dimension(1:ndim, 1:Fdof, 1:elem%Tdof) :: rhs

      real, dimension(1:ndim, 1:Fdof, 1:tQnum) :: space_rhs ! space part of the rhs in time integration nodes

      real, dimension(:,:,:,:,:), allocatable :: sigma   ! 0:numberOfEdges(=3),1:nbDim, 1:tdeg, 1:G(Q)dof,  1:ndim
      real, dimension(:,:,:,:,:), allocatable :: sigmaIPG   ! 1:numberOfEdges(=3),1:nbDim, 1:tdeg, 1:MGdof,  1:ndim

      integer :: i, m,d
      integer :: Qdof, Mdof, MGnum, TDof

      Qdof = elem%Qdof
      Tdof = elem%Tdof

      MGnum = maxval( elem%face(fGnum,:) )
      Mdof = max( elem%Qdof , state%space%G_rule(MGnum)%Qdof )


      ! go through time integration nodes
      ! compute flux at k-th time level - Set the RHS for the local Moment Problem

      !allocate( sigma(0:3, 1:nbDim, 1:tQnum, 1:Mdof, 1:ndim) )
      allocate( sigma(0:3, 1:ndim,1:nbDim, 1:Mdof, 1:tQnum), source = 0.0 )


      !NIPG/SIPG need a boundary member in volume sigma
      if (state%space%m_IPG /= 0 ) then

         allocate( sigmaIPG(1:3, 1:ndim,1:nbDim, 1:Mdof, 1:tQnum), source = 0.0 )

         !!! Set the flux sigma on the edge and volume also
         call SetRTNstFluxesScalar( elem, tQnum, Mdof, &
            sigma(0:3, 1:ndim, 1:nbDim, 1:Mdof, 1:tQnum), &
            sigmaIPG(1:3, 1:ndim, 1:nbDim, 1:Mdof, 1:tQnum)  )

!         print*, 'IPG:', sigmaIPG(1,1,1,:,:)
        ! stop 'in SetRTNstRhs'

         ! compute the space part of the RHS in time integration nodes 1:tdeg
         space_rhs(1:ndim,1:Fdof,1:tQnum) =  SetRTNstSpaceRhs( elem, tQnum, Fdeg, Fdof, Mdof, &
            sigma(0:3, 1:ndim,1:nbDim, 1:Mdof, 1:tQnum), &
            sigmaIPG(1:3, 1:ndim,1:nbDim, 1:Mdof, 1:tQnum) )

      !IIPG
      else
         !!! Set the flux sigma on the edge and volume also
         call SetRTNstFluxesScalar( elem, tQnum, Mdof, sigma(0:3, 1:ndim,1:nbDim, 1:Mdof, 1:tQnum) )
         ! compute the space part of the RHS in time integration nodes 1:tdeg

         space_rhs(1:ndim,1:Fdof,1:tQnum) =  SetRTNstSpaceRhs( elem, tQnum, Fdeg, Fdof, Mdof, &
            sigma(0:3, 1:ndim,1:nbDim, 1:Mdof, 1:tQnum) )
      endif

      !time integration of the RHS -> rhs for all test function
      associate ( Tphi => state%time%T_rule(tQnum)%phi )

      do m = 1, Tdof
         do i = 1, Fdof
            do d = 1,ndim
            ! in fact the time basis func. satisfy \phi_*\phi_j = \tau \delta_{ij}
            ! hence we need to divide the RHS by \tau, which is done by integrating over 0,1 instead od 0,\tau in EvalTimeScalarProduct
!            rhs(d,i,m) = EvalTimeScalarProduct( &
!               tQnum, state%time%tau(1), Tphi(m, 1:tQnum), space_rhs(d,i,1:tQnum) )
            rhs(d,i,m) = EvalTimeScalarProduct( &
               tQnum, 1.0 , Tphi(m, 1:tQnum), space_rhs(d,i,1:tQnum) )

            end do !d
         end do !i
      end do !m
      end associate !Tphi

      !print*, 'Warning in SetRTNstRhs (EvalTimeScalarProduct)'

   end function SetRTNstRhs

   !> export rtnst solution to file rtn_m.txt
   !> this is used to compute the dual solution which enables to compute the exact error EST(u_h)
   !> works only if Omega = (0,1)^2  and the triangles are half squares, (LB corner is a node)!!!!
   subroutine exportRTNst( grid , spaceLevels, timeLevels )
      class( mesh ) , intent(in) :: grid
      ! how many inner values in each elem is computed (degree of Lagrange nodes )
      integer, intent(in) :: spaceLevels
      ! number of equidistantly distributed time levels + 1 , i.e. # of nodes
      integer, intent(in) :: timeLevels
      class(element), pointer :: elem
      integer :: i, j, k, ifile, N, i1, j1, l, gi, kk
      character(len=90) :: fileName
      character(len=90) :: folderName
      real, allocatable, dimension(:,:) :: Fxi ! solution in Lagr. nodes
      integer :: nodes, nelem, nCoord, dof, Tdof
      real :: timeNode, realTime, hx
      real :: max_dK, min_dK
      type(Lagrang_rule), pointer :: L_rule
      real, allocatable, dimension(:,:,:) :: sol , F
      real, allocatable, dimension(:,:) :: w_jump, Re_1
      real, allocatable, dimension(:,:) :: w, wDt, rhs, divF
      real, allocatable, dimension(:,:,:) :: wDx, dWdx, R_s, f_s
      real, dimension(2) :: lb_corner, rt_corner
      integer :: steps

      call setDomainCorners( state%model%icase, lb_corner , rt_corner )

      ! used for control only ! may be commented
      call PlotSolPrimal( state%time%iter )

      nelem = grid%nelem

      ! Parameters of the computation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! number of element in x-direction
      N = nint( sqrt( nelem / 2. ) )

      !  x-size of the little triangle, =  x-sizeofTriangle / spaceLevels
      !  y-size = x-size since Omega is supposed to be square domain
      hx = ( rt_corner(1) - lb_corner(1) ) / ( N * spaceLevels  )
      nCoord = N * spaceLevels + 1! number of node in x-direction
!      print*, 'hx = ', hx, nCoord

      ! maxval of the weighting parameters - it should be the same for all elements
      max_dK = maxval( grid%elem(:)%dK )
      min_dK = minval( grid%elem(:)%dK )
      if ( abs( max_dk / min_dK ) > 5 ) then
         print*, 'max DK:', max_dK, 'min dK:', min_dK
         stop 'max DK is more than 10x min DK! Could be problem when export'
      endif

      ! max_degree implemented
      if ((spaceLevels > maxLrule) ) then !.or. (timeLevels > maxTrule)) then
         print*, 'The needed degree of quadrature is larger than implemented!'
         print*, 'spaceLevels:', spaceLevels, 'maxLrule:', maxLrule !,  &
                 !'timeLevels:', timeLevels, 'maxTrule:', maxTrule
         stop
      endif

      L_rule => state%space%L_rule( spaceLevels )
      ! physical coordinates of the Lagrange nodes
      allocate( Fxi(1:L_rule%Qdof, 1:2 ) )

      ! number of nodes in each element
      nodes = (spaceLevels+1)*(spaceLevels+2)/2
      if ( nodes /= L_rule%Qdof ) then
         print*, 'nodes:=' , nodes, L_rule%Qdof
         stop 'wrong number of nodes vs Lagrange Qdof'
      endif


      ! needed in Set_R_s
      allocate(Re_1(1:iRe, 1:nodes) )
      Re_1(1:iRe, 1:nodes) = state%model%Re1

      ifile = 29
      ! write the number of levels into a file
      ! folder
      steps = nint(state%time%FinTime / state%time%tau(1) )
      !print*, 'step:', steps, 'tau' , state%time%tau(1)

      write(folderName,'(A11,I2.2,A2,I2.2,A2,I2.2,A6,I6.6,A6,I5.5,A1)') &
                        'export/case', state%model%icase, &
                        '_p', state%space%deg, &
                        '_q' , state%time%deg, '_nelem', grid%nelem, &
                        '_steps', steps, '/'
      call system('mkdir -p ' // adjustl(trim( folderName ) ) )

      ! write the ini file initDualProblem.txt
      call writeInitDualProblemFile( folderName, grid, spaceLevels, timeLevels )

      write(fileName, '(A,A4,I4.4,A4)')  trim(folderName) , 'rtn_', state%time%iter, '.txt'
      fileName = adjustl(trim( fileName ) )

      open(ifile, file = fileName, action="write", status="UNKNOWN" )
      ! number of noded in x-direction, number of elems in row nelem = ((n_coord-1)*2)**2, # time levels, time-begin, time-end
      write(ifile,*) lb_corner(1:2), rt_corner(1:2)
      write(ifile,*) nCoord, grid%h, timeLevels, &
                     state%time%ctime - state%time%tau(1) , state%time%ctime, max_dK

      ! solution at each element in Lagrange nodes
      ! 1 - x coord
      ! 2 - y coord
      ! 3 - time coord
      ! 4 - value
      allocate( sol(1:grid%nelem, 1:nodes,1:6), source = 0.0 )
      allocate( F(1:nCoord, 1:nCoord,1:8), source = 0.0 )
      ! jump of the solution in Lagrange nodes
      allocate( w_jump(1:ndim, 1:nodes), source =0.0 )

      allocate( w(1:ndim,1:nodes), source = 0.0 )
      allocate( wDt(1:ndim,1:nodes), source = 0.0 )
      allocate( wDx(1:ndim,1:nodes, 1:nbDim), source = 0.0 )
      allocate( dWdx(1:nodes,1:ndim,1:nbDim), source = 0.0 )
      allocate( rhs(1:ndim, 1:nodes), source = 0.0 )
      !allocate( divF(1:nodes, 1:ndim), source = 0.0 )

      allocate( R_s(1:nodes, 1:nbDim, 1:ndim), source = 0.0 )
      allocate( f_s(1:nodes, 1:nbDim, 1:ndim), source = 0.0 )


      do k = 1, timeLevels

         ! zero arrays
         F(1:nCoord, 1:nCoord,1:8) = 0.0
         ! real physical time
         realTime = ( state%time%ctime - state%time%tau(1) ) + ((k-1)* state%time%tau(1) / (timeLevels-1))
         ! relative time in [0,1] used for computing w in time moments
         timeNode = (k-1) / (timeLevels-1.0)
!         print*, 'time node = ', timeNode, realTime

         ! compute the solution (and other functions) for each element in all Lagrange nodes
         do i = 1,grid%nelem
            elem => grid%elem(i)
            dof = elem%dof
            Tdof = elem%Tdof
            !coordinates
            call ComputeF( elem , nodes, L_rule%lambda( 1:nodes, 1:2 ) , &
                              Fxi( 1:nodes, 1:nbDim) )

            ! compute the residual !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! add the jump in the solution {u}_m-1
            if ( k == 1) then
               !                                          lagrDeg
               w_jump(1:ndim,1:nodes) = evalwSTjumpInLagrangeNodes( elem , spaceLevels )
            else
               w_jump(1:ndim,1:nodes) = 0.0
            end if

            ! compute w, dw/dt, grad(w) in Lagrange nodes
            call evalSTfunInRealTimeLagrange_dxdt( elem, ndim, dof, Tdof, &
                                          elem%wST(1:ndim, 1:dof, 1:Tdof), &
                                          timeNode, spaceLevels, nodes, &
                                          w(1:ndim,1:nodes), &
                                          wDt(1:ndim, 1:nodes), &
                                          wDx(1:ndim, 1:nodes, 1:nbDim) )
            ! transpose - needed to Set_R_s
            do kk = 1,ndim
               dWdx(1:nodes, kk, 1:nbDim ) = wDx(kk, 1:nodes, 1:nbDim)
            end do ! kk

            ! diffusive terms
            call Set_R_s_scalar( ndim, nbDim, iRe, nodes, transpose(w(1:ndim, 1:nodes)), &
                          dWdx(1:nodes,1:ndim, 1:nbDim), Re_1, &
                          R_s(1:nodes, 1:nbDim, 1:ndim), Fxi( 1:nodes, 1:nbDim) )

            ! convective terms
            call Set_f_s_scalar( ndim, nbDim, nodes, transpose(w(1:ndim, 1:nodes)), &
                          f_s(1:nodes, 1:nbDim, 1:ndim), Fxi( 1:nodes, 1:nbDim), &
                          elem%i )
            ! total flux
            f_s(1:nodes, 1:nbDim, 1:ndim)  = R_s(1:nodes, 1:nbDim, 1:ndim) &
               - f_s(1:nodes, 1:nbDim, 1:ndim)

            ! we do NOT use divF, but directly the fluxes
            ! try
            ! call EvalDiv_F(elem, nodes-2, nodes, Fxi(1:nodes, 1:nbDim), &
            !      f_s(1:nodes, 1:nbDim, 1:ndim), DivF(1:nodes, 1:ndim) )


            ! eval rhs
            do gi=1,nodes   ! integ nodes
               call RHS_Scalar( Fxi(gi, 1:nbDim), rhs(1:ndim,gi), realTime )
            end do ! gi

            sol(i, 1:nodes,1:2) = Fxi( 1:nodes, 1:2 )
            !sol(i, 1:nodes,3) = realTime ! not needed

            sol(i, 1:nodes,3) = rhs(1, 1:nodes) - wDt(1,1:nodes)
            ! flux = - \sigma
!            sol(i, 1:nodes,4:5) = (-1.) * wDx(1,1:nodes,1:nbDim)
            sol(i, 1:nodes,4:5) = (-1.) * f_s(1:nodes,1:nbDim,1)
            ! -jump of the solution
            sol(i, 1:nodes, 6) = (-1.) * w_jump(1,1:nodes)

         end do


         ! Fill F( 1:N*spaceLevels + 1, 1:N*spaceLevels + 1, 1: 5)
         ! move the functions saved element-wise to global array
         ! 1-3 - coords
         ! 4 - how many times a value was added - we need to do an average in the end
         ! 5 = f - u'
         ! 6-7  = -\sigma
         ! 8 = - {u_h}_m-1
         ! since one value may filled from multiple elements

         do i = 1,nelem
            do l = 1, nodes
               ! x-index: i = ( x_i- lbCorner_x )/(little_h) + 1 ( x_i = lb_corner_x +  (i-1)*h_l
               i1 = nint( ( ( sol(i,l,1) - lb_corner(1) ) / hx ) + 1 )
               if ( abs(i1 -  ((( sol(i,l,1) - lb_corner(1) ) / hx ) + 1) ) > 1E-5 ) then
                  print*, 'i1 and sol:', i, l, i1, ((sol(i,l,1) - lb_corner(1))/hx)+1
                  stop 'rounding error in exportRTNst'
               endif
               !y-index ! nint -closest int
               j1 = nint(( (sol(i,l,2)-lb_corner(2)) / hx ) + 1 )
               ! fill the coords
               F(i1,j1, 1:3) = (/ sol(i,l,1) , sol(i,l,2), realTime /)
               ! fill the counter needed for the average
               F(i1,j1, 4) = F(i1,j1, 4) + 1
               ! fill the values
               ! f - u'
               F(i1,j1, 5) = F(i1,j1, 5) + sol(i,l,3)
               !-\sigma
               F(i1,j1, 6:7) = F(i1,j1, 6:7) + sol(i,l,4:5)
               ! - jump
               F(i1,j1, 8) = F(i1,j1, 8) + sol(i,l,6)
               ! TRY - per partes without the edges
               !F(i1,j1, 9) = F(i1,j1, 9) + sol(i,l,7)

            end do ! l
         end do ! i

         ! write into the file
         do j = 1, nCoord ! N*spaceLevels + 1
            do i = 1, nCoord ! N*spaceLevels + 1
               write(ifile, * ) F(i,j,1:3),  F(i,j,5) / F(i,j,4) , &
                                F(i,j,6) / F(i,j,4), F(i,j,7) / F(i,j,4), &
                                F(i,j,8) / F(i,j,4) , &
                                F(i,j,4)
                                !F(i,j,9) / F(i,j,4)  ! TRY divF
            end do ! i
         end do ! j

      end do ! k

      deallocate( sol, F, w_jump, w, wDt, wDx, dWdx )
      deallocate( R_s, f_s )

      deallocate( Fxi, Re_1 )
      nullify( L_rule )

      close(ifile)

   end subroutine exportRTNst

   !> set  the corners of the domain, used dividing the mesh in Lagrange nodes
   !> in order to export the problem to Fenics
   subroutine setDomainCorners( icase, lb_corner, rt_corner )
      integer, intent(in) :: icase ! isca
      real, dimension(1:nbDim), intent(out) :: lb_corner, rt_corner

       !CASE 2
      select case ( icase )
         case ( 2 )
            print*, ' SQUARE DOMAIN (-1,1)^2 ! '
            lb_corner = (/ -1. , -1. /)
            rt_corner = (/ 1. , 1. /)
         ! CASE 24 BARRENBATT
         case ( 24 )
            print*, ' SQUARE DOMAIN (-6,6)^2 ! '
            lb_corner = (/ -6. , -6. /)
            rt_corner = (/ 6. , 6. /)

         case ( 42 )
            print*, ' SQUARE DOMAIN (0,1)^2 ! '
            lb_corner = (/ 0. , 0. /)
            rt_corner = (/ 1. , 1. /)
         case default
            stop 'Select the size of the domain for export in exportRTNst!'
      end select

   end subroutine setDomainCorners

   subroutine writeInitDualProblemFile( folderName, grid, spaceLevels, timeLevels )
      character(len=90) :: folderName
      class( mesh ), intent(in) :: grid
      integer, intent(in) :: spaceLevels, timeLevels
      integer :: ifile
      character(len=90) :: fileName

      ifile = 29

      write(fileName, *)  trim(folderName) , 'initDualProblem.txt'
      fileName = adjustl(trim( fileName ) )

      ! it is overwritten in each time step - so in the end the global values should be there
      open(ifile, file = fileName, action="write", status="UNKNOWN")
         write(ifile,*)  nint(state%time%FinTime / state%time%tau(1) ), &
                         spaceLevels, timeLevels
         write(ifile,*)  grid%nelem, grid%h, state%time%tau(1), state%space%deg, state%time%deg
         write(ifile,*)  state%model%icase, state%model%Re1, state%space%m_IPG
         write(ifile,*)  sqrt( state%estim( RTNrez, 1 ) ), &
                         sqrt( state%estim( RTNflux, 1 ) ), &
                         sqrt( state%estim( RTNradau2, 1 ) ), &
                         sqrt( state%estim( RTNjump, 1 ) )
         write(ifile,*)  sqrt( state%estim( RTNeta, 1) ), &
                         sqrt( state%estim( RTNfluxnorm, 1 ) )
      close(ifile)

   end subroutine writeInitDualProblemFile

   !> set constants use for weighting the estimates in RNTstElem
   subroutine setWeightingConstants( elem, dK, Cnc )
      class(element), intent(in) :: elem
      real, intent(out) :: dK , Cnc
      real :: h , Ckb

      ! ERN weighting d_K^2 = CTn^-1
      dK = ( ( elem%Cb + elem%CK ) / elem%diam**(2.0) &
                  + state%time%FinTime / state%time%tau(1)**(2.0) )**(-0.5)
      !dK = sqrt( elem%diam**(2.0) + state%time%tau(1)**(2.0) )
      !dK = sqrt( elem%diam**(2.0) / ( elem%Cb + elem%CK  ) &
       !               + state%time%tau(1)**(2.0)/state%time%FinTime )

      ! Constant in front of the NC term
      ! instead h there should be h_Gamma !
      h = elem%diam
      ! may be quite LARGE ???
      Ckb = (elem%CKo**2)/h + h*elem%Cbo**2 + (state%model%Re1*state%space%sigma)**2 / h
      Cnc = dK**(2.0) * Ckb / h**2.
      ! OLD version
      !Cnc = dK**(2.) * h**(-3.)

      if (elem%i == 1) then
         print*, 'Actual dK:', dK
         print*, 'CONSTANTS h_2:', h**(-2.), 'tau-2:', state%time%tau(1)**(-2.)
         print*, 'CKo:' , elem%CKo, 'Cbo:' , elem%Cbo , 'eps:' , state%model%Re1
         print*, 'CK:' , elem%CK, 'Cb:' , elem%Cb, 'Ckb:', Ckb
         print*, 'Cnc old:' , dK**(2.) * h**(-3.) , &
                  'new:' , (dK**(2.0) * Ckb) / (h**2.)
         print*, 'DK old:', sqrt( elem%diam**(2.0) / ( elem%Cb + elem%CK  ) &
                      + state%time%tau(1)**(2.0)/state%time%FinTime ), &
                  'new:' , ( ( elem%Cb + elem%CK ) / elem%diam**(2.0) &
                  + state%time%FinTime / state%time%tau(1)**(2.0) )**(-0.5)
      endif
      !print*, 'dk, cp, h:', dk, elem%CP, elem%diam, state%time%tau(1)

   end subroutine setWeightingConstants

end module rtn_st_mod
