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
   public :: SetRTNstRhs



contains

!> estimates the computational error using the approach based on the flux reconstruction
   subroutine ComputeRTNstEstim( estim , tQnum )
      real, dimension(1:max_eta,1:ndim), intent(out) :: estim
      class(element), pointer :: elem
      integer, intent(in) :: tQnum
      integer :: i
      integer :: Fdeg, tdeg
      class( Time_rule), pointer :: T_rule

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
         estim( RTNradau, 1:ndim ) = &
            estim( RTNradau, 1:ndim ) + ( elem%eta( RTNradau, 1:ndim ) )**(2.0)
         !New version
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



!         print*, 'Eta_rez in the ', state%time%iter,'-th time step:', sqrt( estim( RTNrez, 1:ndim ) )
!         print*, 'Eta_flux in the ', state%time%iter,'-th time step:', sqrt( estim( RTNflux, 1:ndim ) )
!         print*, 'Eta_radau in the ', state%time%iter,'-th time step:', sqrt( estim( RTNradau, 1:ndim ) )
!         print*, 'Eta_jump in the ', state%time%iter,'-th time step:', sqrt( estim( RTNjump, 1:ndim ) )
!         stop 'after first element RTNst estims!'

      enddo

      print*, 'Eta_rez in the ', state%time%iter,'-th time step:', sqrt( estim( RTNrez, 1:ndim ) )
      print*, 'Eta_flux in the ', state%time%iter,'-th time step:', sqrt( estim( RTNflux, 1:ndim ) )
      print*, 'Eta_radau in the ', state%time%iter,'-th time step:', sqrt( estim( RTNradau, 1:ndim ) )
      print*, 'Eta_radau2 in the ', state%time%iter,'-th time step:', sqrt( estim( RTNradau2, 1:ndim ) )
      print*, 'Eta_jump in the ', state%time%iter,'-th time step:', sqrt( estim( RTNjump, 1:ndim ) ), '^2' ,estim( RTNjump, 1:ndim )
      print*, 'Flux norm in the ', state%time%iter,'-th time step:', sqrt( estim( RTNfluxnorm, 1:ndim ) )


     ! print*, 'Warning dK set to 1.0!!!'


!      print*, 'Estimator in the ', state%time%iter,'-th time step:', sqrt( estim( RTNall, 1:ndim ) ) , &
!         '+', sqrt( estim( RTNjump, 1:ndim ) )

!      stop 'End of ComputeRTNstEstim!'

   end subroutine ComputeRTNstEstim


   !> compute error estimates via the dual (residual) form using RTN flux reconstruction for one element
   subroutine ComputeElemRTNstEstim( elem , tQnum )
      class(element), intent(inout) :: elem
      integer, intent(in) :: tQnum !number of time moments

!      real, dimension(1:max_eta, 1:ndim), intent(inout) :: estim
      !type(basis_rtn_fe), pointer:: loc_RTN
      real, dimension(:,:), allocatable :: Mspace ! inverse of Momentum matrix on elem
      real, dimension(:,:,:), allocatable :: flux   ! basis coefficients of flux (1:ndim - future), 1:tdeg, 1:Fdof
      real, dimension(:,:,:), allocatable :: RadauReconstruct
      real, dimension(:,:,:), allocatable :: rhs !, space_rhs
      integer :: Fdeg, Fdof
      integer :: Qdof, Tdof
      real :: dK, cTn, cKn
      real, dimension(:,:,:), allocatable :: fun
      character(len=20) :: outputfile


      !print*, 'ComputeElemRTNstEstim for elem(', elem%i, '):'

!    ! scaling factor computed in errorFlux.f90
!    !elem%CTn = state%time%FinTime/state%time%tau(1)**2 + (elem%CK + elem%Cb)/elem%diam**2

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
      elem%eta( RTNfluxnorm, 1:ndim ) = RTNstComputeFluxNorm( Set_R_s_scalar, Set_f_s_scalar, elem )
!      print*, 'fluxnorm:' , elem%eta( RTNflux, 1:ndim )


   !!!!!!! Weight the estimates !!!!!!!!!!!!!!!!!!!
      !OLD VERSION
      dK = sqrt( elem%diam**(2.0) + state%time%tau(1)**(2.0) )

      ! Temporarily
      !dK = 1.0


      ! print*, 'dk, cp, h:', dk, elem%CP, elem%diam, state%time%tau(1)

      elem%eta(RTNrez,1:ndim) = elem%CP * dK * sqrt(elem%eta(RTNrez,1:ndim))
      elem%eta(RTNflux,1:ndim) = ( dK / elem%diam ) * sqrt(elem%eta(RTNflux,1:ndim))
      !H1L2 norm
      elem%eta(RTNradau,1:ndim) = dK * sqrt(elem%eta(RTNradau,1:ndim))
      !ST L2 norm
      elem%eta(RTNradau2, 1:ndim) = ( dk / state%time%tau(1) ) * sqrt(elem%eta(RTNradau2, 1:ndim))
      !RTNjump is SQUARED
      elem%eta(RTNjump, 1:ndim) = dK**(2.) * elem%diam**(-3.) * elem%eta(RTNjump, 1:ndim)

      elem%eta(RTNfluxnorm, 1:ndim) = sqrt( elem%eta(RTNfluxnorm, 1:ndim) )


!   !!!!!! TRY second version from Vohralik !!!!!!!!!!!!!!!!
!      !                                                        CK ???
!      cTn = ( state%time%FinTime / state%time%tau(1)**2.0 ) + ( 1.0 / elem%diam**(2.0) )
!      cKn = (state%space%sigma**(2.0) + 1.0 ) / elem%diam
!
!      elem%eta(RTNrez,1:ndim) = elem%CP * cTn**(-0.5) * sqrt(elem%eta(RTNrez,1:ndim))
!      elem%eta(RTNflux,1:ndim) = cTn**(-0.5) * elem%diam**(-1.0) * sqrt(elem%eta(RTNflux,1:ndim))
!      elem%eta(RTNradau,1:ndim) = cTn**(-0.5) * sqrt(elem%eta(RTNradau,1:ndim))
!
!      !RTNjump is SQUARED
!      elem%eta(RTNjump, 1:ndim) = cTn**(-1.0) * elem%diam**(-2.) *cKn* elem%eta(RTNjump, 1:ndim)
!!END Vohralik !!!!!!


!      print*, 'ETA Rez:', elem%eta(RTNrez,:)
!      print*, 'ETA Flux:', elem%eta(RTNflux,:)
!      print*, 'ETA Radau:', elem%eta(RTNradau,:)
!      print*, 'ETA Jump:', elem%eta(RTNjump,:)

      deallocate(Mspace, rhs, flux, RadauReconstruct)

      !stop 'END of ComputeElemRTNstEstim'

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


end module rtn_st_mod
