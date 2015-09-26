module eval_rtn_st_mod
      use blocks_integ
   use data_mod
   use define_state
   use element_mod
   use euler_problem
   use eval_rav_tho_ned
   use eval_sol
   use geometry, only: WriteLocalMatrix
   use loc_rav_tho_ned
   use main_data
   use stdgm_mod


implicit none

   public :: evalRadauReconstruct
   public :: EvalRTNFluxEdgeIPG
   public :: EvalRTNFluxEdgeNormal
   public :: EvalRTNFluxVolume
   public :: evalRTNstFunDiv_Q
   public :: EvalRTNstFunInIntNodes
   public :: EvalSTsigma ! compute the advective-diffusive flux from the equation

   public :: RTNstFluxEstimator
   public :: RTNstRezEstimator

   public :: Set_RTN_integ_nodeNew
   public :: SetRTNstFluxesScalar
   public :: SetRTNstFluxReconstr
   public :: SetRTNstSpaceRhs


contains

    !> compute the jump estimator \sum_{\Gamma_T} || \jump{v}||^2_{T \times I_m}
   function RTNstJumpEstimator( elem, tQnum ) result (eta)
      class(element), intent(inout) :: elem
      integer, intent(in) :: tQnum
      real, dimension(1:ndim) :: eta

      real, dimension(:,:), allocatable :: jump_wi
      real, dimension(1:ndim) :: temp
      real, dimension(1:ndim,1:tQnum) ::jump
      real :: cTime, tau, dK
      integer :: i, iEdge, m, Gnum, Gdof
      type( Gauss_rule ), pointer :: G_rule

      cTime = state%time%ctime
      tau = state%time%tau(1)

      jump(:,:) = 0.0

      do m = 1, tQnum

         !set the right %ctime - boundary elements need to take value from the exact sol
         state%time%ctime =  (state%time%ttime - tau ) + &
               tau * state%time%T_rule(tQnum)%lambda( m )

         !to computation of fluxes we need the solution w in time integ nodes for all neighboring elements
         call Transfer_wST_to_w_Patch(elem, m, tQnum)

         !edges
         do iEdge = 1, elem%flen

            Gnum = elem%face(fGnum,iEdge)
            G_rule => state%space%G_rule(Gnum)
                  !eval \hat\sigma on the edge (multiplied by the outer normal)
            Gdof = G_rule%Qdof

            allocate( jump_wi(1:Gdof,1:ndim) )
            call ElementEdgeJump( elem, iEdge, jump_wi(1:Gdof,1:ndim) )

            call IntegrateFunctionVecEdge( elem, iEdge, ndim, jump_wi(1:Gdof, 1:ndim)**(2.0) , temp(1:ndim) )
            deallocate( jump_wi )

            jump (1:ndim, m) = jump(1:ndim, m) + temp(1:ndim)

         end do !iEdge

      end do !m

      !dK = sqrt( elem%diam**(2.) + tau**(2.) )

      eta(1:ndim) = IntegrateTimeFunctionVec( tQnum, tau, jump(1:ndim, 1:tQnum) )


      state%time%ctime = cTime

   end function RTNstJumpEstimator

   !> compute the Rezidual RTNst estimator \f$ f - R' - \nabla \cdot {\bf t}_{h_\tau} \f$
   function RTNstRezEstimator( elem, dof, Fdeg, Fdof, Tdof, tQnum, flux, Radau ) result (eta)
      class(element) , intent(inout) :: elem
      integer, intent(in) :: dof, Fdeg, Fdof
      integer, intent(in) :: Tdof, tQnum
      real, dimension(1:ndim,1:Fdof, 1:Tdof) :: flux
      real, dimension(1:ndim,1:dof, 1:Tdof+1) :: Radau
      real, dimension(1:ndim) :: eta

      real, dimension(1:elem%Qdof, 1:nbDim):: x, Fx
      real, dimension(1:ndim, 1:elem%Qdof, 1:tQnum) :: div_flux
      real, dimension(1:ndim, 1:elem%Qdof, 1:tQnum) :: f
      real, dimension(1:ndim, 1:elem%Qdof, 1:tQnum) :: RadauReconstructDer_Q
      integer :: Qdof, iNode, l
      real :: cTime, tau

      Qdof = elem%Qdof
      tau = state%time%tau(1)

      ! eval divergence of the flux reconstruction in ST integ nodes
      div_flux(1:ndim, 1:Qdof, 1:tQnum) = evalRTNstFunDiv_Q( elem, Fdeg, Fdof, Tdof, tQnum, &
         flux(1:ndim, 1:Fdof, 1:Tdof) )

      ! compute RHS f in ST integ nodes
      x(1:Qdof, 1:nbDim) = state%space%V_rule(elem%Qnum)%lambda(1:Qdof,1:nbDim)
      call ComputeF(elem, Qdof, x, Fx)

      do iNode = 1, tQnum
         !Warning cTime!!!!!!!!!!!!!!
         cTime = (state%time%ttime - tau) + tau*state%time%T_rule(tQnum)%lambda(iNode)
         do l=1, Qdof
            call RHS_Scalar( Fx(l, 1:nbDim), f(1:ndim, l, iNode), cTime )
         enddo !l
      end do !iNode

!      print*, 'RHS = ', f(1,1,1)


      ! eval RadauReconstruct in ST integ nodes
      do iNode = 1, tQnum
         RadauReconstructDer_Q(1:ndim, 1:Qdof, iNode) = evalSTfunInIntTime_Der( elem, &
            tau, ndim, elem%dof, Tdof+1, Radau(1:ndim, 1:elem%dof, 1:Tdof + 1 ), iNode, tQnum )
      end do !iNode

!      print*, 'Radau: ' , RadauReconstructDer_Q(1:ndim, :, :)
!      print*, 'div flux:', div_flux(1:ndim, 1:Qdof, 1:tQnum)
!      print*,
!      print*, 'minus:', RadauReconstructDer_Q(1:ndim, :, :) - div_flux(1:ndim, 1:Qdof, 1:tQnum)
!      print*, 'plus:', RadauReconstructDer_Q(1:ndim, :, :) + div_flux(1:ndim, 1:Qdof, 1:tQnum)

      f(1:ndim, 1:Qdof, 1:tQnum) = f(1:ndim, 1:Qdof, 1:tQnum) &
         - RadauReconstructDer_Q(1:ndim, 1:Qdof, 1: tQnum) - div_flux(1:ndim, 1:Qdof, 1:tQnum)

!      print*, 'f = ' , f

      eta(1:ndim) = evalL2STNormQ_Elem( elem, tau, tQnum, f(1:ndim, 1:elem%Qdof, 1:tQnum) )

!      if ( eta(1) > 1.E-20) then
!       if (elem%i == 8) then
!         print*, 'Elem(', elem%i ,')ETA Rez:', eta(:)
!
!         print*, 'flux basis coefs' , flux(1,1,:)

!         print*, 'Radau basis coefs:', Radau(1,1,:)
!
!         do iNode = 1, tQnum
!
!            u(1:ndim,1:Qdof,iNode) = evalSTfunInIntTime(elem, ndim, elem%dof, &
!               elem%Tdof, elem%wST, iNode, tQnum)
!         end do
!
!         print*, 'Radau: ' , RadauReconstructDer_Q(1:ndim, 1, :)
!         !print*, 'U = ', u(1,1,1)
!         print*, 'div flux:', div_flux(1:ndim, 1, :)
!         print*,
!         print*, 'plus:', RadauReconstructDer_Q(1:ndim, 1,:) + div_flux(1:ndim, 1, : )
!         print*,
!         print*, 'f: ', f(1,1,:)

!      endif
!      stop

   end function RTNstRezEstimator

   function RTNstFluxEstimator( elem, Fdeg, Fdof, Tdof, tQnum, flux ) result (eta)
      class(element) , intent(inout) :: elem
      integer, intent(in) :: Fdeg, Fdof
      integer, intent(in) :: Tdof, tQnum
      real, dimension(1:ndim,1:Fdof, 1:Tdof) :: flux

      real, dimension(1:ndim) :: eta

      real, dimension(:,:,:,:), allocatable :: fluxQ, fluxCD
      integer :: Qdof

      Qdof = elem%Qdof

      ! flux in space and time integ. nodes
      allocate( fluxQ(1:ndim, 1:nbDim, 1:Qdof, 1:tQnum) )


      fluxQ(1:ndim, 1:nbDim, 1:Qdof, 1:tQnum) = &
         EvalRTNstFunInIntNodes( elem, Tdof, tQNum, Fdeg, Fdof, flux(1:ndim,1:Fdof, 1:Tdof) )


      allocate( fluxCD(1:ndim, 1:nbDim, 1:Qdof, 1:tQnum) )
      fluxCD(1:ndim, 1:nbDim, 1:Qdof, 1:tQnum) = EvalSTsigma( elem,  tQnum)

      ! integrate flux and fluxCD
      eta(1:ndim) = evalL2STNormVecQ_Elem( elem, state%time%tau(1), tQnum, &
         fluxCD(1:ndim, 1:nbDim, 1:Qdof, 1:tQnum) - fluxQ(1:ndim, 1:nbDim, 1:Qdof, 1:tQnum) )

      deallocate( fluxQ, fluxCD )

   end function RTNstFluxEstimator


 !> compute the flux error for element STDG
  !> \f$ c_f \|u_h -u\| + c_K \|R_s(u_h) - R_s(u) \| + c_f \|f_s(u_h) - f_s\| \f$
   function RTNstComputeFluxNorm( Set_R_s, Set_f_s, elem ) result(eta)
   interface
      subroutine Set_R_s(ndimL, nbDim, Qdof, w, Dw, Re_1, R_s)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
      end subroutine Set_R_s
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
      end subroutine Set_f_s
    end interface
      class(element), intent(inout) :: elem

      real, dimension(1:ndim) :: eta

      real, dimension(:,:,:,:), allocatable :: fluxUh
      real, dimension(:,:), allocatable :: wE
      real, dimension(:,:,:), allocatable :: R_sE, f_sE
      real, dimension(:,:,:), allocatable :: DwE
      real, dimension(:,:,:,:), allocatable :: func
      real, dimension(:), allocatable :: Re1
      real ::  cTime
      integer :: Qdof, alpha, i, tQnum


      cTime = state%time%ctime
      Qdof = elem%Qdof
      allocate( Re1(1:Qdof) )
      Re1(1:Qdof) = state%model%Re1
      tQnum = elem%tQnum


!      print*, 'Re1:' , Re1(1)

      allocate ( fluxUh( 1:ndim, 1:nbDim, 1:elem%Qdof, elem%tQnum ) )

      allocate ( wE(1:Qdof, 1:ndim), DwE(1:Qdof, 1:ndim, 1:nbDim), source = 0.0 )

      allocate( func(1:ndim,1:nbDim,1:elem%Qdof,1:tQnum), source = 0.0 )

      allocate ( f_sE(1:Qdof, 1:nbDim, 1:ndim), R_sE(1:Qdof, 1:nbDim, 1:ndim), source = 0.0 )

      do alpha = 1, elem%tQnum

            !set the right %ctime - boundary elements need to take value from the exact sol
            state%time%ctime =  (state%time%ttime - state%time%tau(1) ) + &
               state%time%tau(1) * state%time%T_rule(tQnum)%lambda(alpha)


            ! setting of the exact solution in integ nodes
            call SetExactSolutionQnodes(elem, state%space%V_rule(elem%Qnum), &
                wE(1:Qdof, 1:ndim), DwE(1:Qdof, 1:ndim, 1:nbDim) )

            call Transfer_wST_to_w_Elem(elem , alpha, tQnum)

            ! DIFFUSIVE TERMS
            ! exact value
            call Set_R_s(ndim, nbDim, Qdof, wE(1:Qdof,1:ndim), DwE(1:Qdof, 1:ndim, 1:nbDim), Re1(1:Qdof),&
                R_sE(1:Qdof, 1:nbDim, 1:ndim) )


                ! CONVECTIVE TERMS
            ! exact value
            call Set_f_s(ndim, nbDim, Qdof, wE(1:Qdof,1:ndim), f_sE(1:Qdof, 1:nbDim, 1:ndim), &
                elem%xi(0,1:Qdof, 1:nbDim) )


            !compute the volume flux \sigma in current time
            fluxUh(1:ndim, 1:nbDim, 1:elem%Qdof, alpha) = EvalRTNFluxVolume( &
                Set_f_s_scalar, Set_R_s_scalar, elem, state%space%V_rule( elem%Qnum ) )

       do i=1,nbDim
          func(1:ndim, i, 1:elem%Qdof, alpha) =  fluxUh(1:ndim, i, 1:elem%Qdof, alpha) &
               -  transpose(f_sE(1:Qdof, i, 1:ndim)) + transpose(R_sE(1:Qdof,i,1:ndim))
       enddo
      end do !alpha

      eta( 1:ndim) =  evalL2STNormVecQ_Elem(elem, state%time%tau(1), elem%tQnum, &
        func(1:ndim,1:nbDim,1:elem%Qdof,1:tQnum) )

      deallocate( wE, DwE, func, fluxUh, Re1 )


   end function RTNstComputeFluxNorm

   !> eval of RTNst function in space and time integration nodes
   !> TODO: psi should be probably allocatable part of elem (is computed twice while computing RTNst estims)
   function EvalRTNstFunInIntNodes(elem, Tdof, tQnum, Fdeg, Fdof, fun) result(funQ)
      class(element), intent(in) :: elem
      integer, intent(in) :: Tdof, tQnum
      integer, intent(in) :: Fdeg, Fdof
      real, dimension(1:ndim,1:Fdof, 1:Tdof), intent(in) :: fun

      real,dimension(1:ndim,1:nbDim,1:elem%Qdof,1:tQnum) :: funQ !function in integ nodes

      class(Time_rule), pointer :: T_rule
      type(volume_rule), pointer :: V_rule
      integer ::  iNode, d, Qdof, nb
      real, dimension(1:Fdof, 1:nbDim, 1:elem%Qdof) :: psi
      real, dimension(1:ndim,1:Fdof,1:tQnum) :: wi

      Qdof = elem%Qdof

      T_rule => state%time%T_rule(tQnum)
      V_rule => state%space%V_rule( elem%Qnum )

      associate ( loc_RTN => state%loc_RTN(Fdeg) )

      call Eval_Loc_RTN_Elem(elem, loc_RTN, psi(1:Fdof, 1:nbDim, 1:Qdof) )

      do iNode = 1,tQnum
        wi(1:ndim, 1:Fdof, iNode) = evalSTfunInIntTimeDof( elem, ndim, Fdof, Tdof, &
            fun(1:ndim,1:Fdof,1:Tdof), iNode, tQnum )

         do d = 1, ndim
            do nb = 1, nbDim
               funQ(d, nb, 1:Qdof, iNode) = matmul( wi(d, 1:Fdof, iNode), psi(1:Fdof, nb, 1:Qdof) )
!                dot_product( psi(1:Fdof, nb, i) , wi(d, 1:Fdof, iNode) )
            end do ! nb
         end do !d

      end do ! iNode

      end associate ! loc_RTN

   end function EvalRTNstFunInIntNodes

   function EvalRTNFluxEdgeNormal( Set_Ppm, Set_R_s, elem, iedge, Gdof) result (flux)
      class(element), intent(inout) :: elem
      integer, intent(in) :: iedge
      integer, intent(in) :: Gdof

      real, dimension(1:ndim,1:Gdof) :: flux
      interface
         subroutine Set_Ppm( ndimL, nbDim, Qdof, w, n, xi, Ppm, one_over_area, elem)
            import :: element
            integer, intent(in) :: Qdof, ndimL, nbDim
            real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
            real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: Ppm
                                                  ! matrices Ppm in  -- " --
            real, dimension(1:Qdof, 1:nbDim), intent(in) :: n   ! outer normal
            real, dimension(1:Qdof, 1:nbDim),intent(in) ::  xi                    ! node on the edge?
            real, intent(in), optional :: one_over_area
            type(element), intent(inout), optional :: elem
         end subroutine
         subroutine Set_R_s(ndimL, nbDim, Qdof, w, Dw, Re_1, R_s)
            integer, intent(in) :: ndimL, nbDim, Qdof
            real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
            real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
            real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
            !real, intent(in) :: Re_1                     ! inverse of Reynolds number
            real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
         end subroutine Set_R_s
      end interface

      real, dimension(:,:), allocatable :: Rflux, Cflux, jump_wi
      integer :: Gnum, Qdof


      Gnum = elem%face(fGnum, iedge)
      associate ( G_rule => state%space%G_rule(Gnum) )
         Qdof = G_rule%Qdof

         !control - may not work for different Gdof - elem%nc in special nodes
         if ( Gdof /= Qdof) &
            stop 'Problem in EvalRTNFluxEdgeNormal: Gdof/=Qdof'

         ! diffusive fluxes
         allocate(Rflux(1:Gdof, 1:ndim) , source = 0.0 )
         call Eval_aver_R_s_Edge(Set_R_s, elem, iedge, Rflux(1:Gdof, 1:ndim) )

         !print*, 'R_flux',iedge,':' , Rflux(1, 1)

         ! convective fluxes
         allocate(Cflux(1:Gdof, 1:ndim), source = 0.0 )
         call Eval_NFlux_Edge(Set_Ppm, elem, iedge, Cflux(1:Gdof, 1:ndim) )

         !print*, 'C_flux:' , Cflux(1, 1)

         allocate( jump_wi(1:Gdof,1:ndim) )
         call ElementEdgeJump(elem, iedge, jump_wi(1:Gdof,1:ndim) )

         !print*, 'Jump:' , jump_wi(1, 1)
         !print*, 'Sigma and RE1:' , state%space%sigma, state%model%Re1

         write(debug,*) ' control the size of the normal and also h^-1 and other coeffs before jump_wi '
         ! total flux, already premultiplied by dn

         !curved edge
         if (elem%ibcur > 0 .and. elem%jcur == iedge ) then
            stop 'Curved edge not implemented in EvalRTNFluxEdgeNormal'
!            flux(1:ndim, 1:Gdof) = transpose( ( -Rflux(1:Gdof, 1:ndim) + Cflux(1:Gdof, 1:ndim)  &
!               + state%space%sigma * state%model%Re1 * jump_wi(1:Gdof, 1:ndim) ) / elem%dnc(1:Gdof) )
         else
            flux(1:ndim, 1:Gdof) = transpose( ( -Rflux(1:Gdof, 1:ndim) + Cflux(1:Gdof, 1:ndim)  &
            + state%space%sigma * state%model%Re1 * jump_wi(1:Gdof, 1:ndim) ) / elem%dn(iedge ) )  ! /dn(Gdof) * dn(Gdof) = 1
         endif

      end associate

   end function EvalRTNFluxEdgeNormal

   function EvalRTNFluxVolume( Set_f_s, Set_R_s, elem, V_rule ) result (flux)
      class(element), intent(in) :: elem
      type(volume_rule), intent(in) :: V_rule

      real, dimension(1:ndim, 1:nbDim, 1:V_rule%Qdof ) :: flux
      interface
         subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x )
            integer, intent(in) :: Qdof, ndimL, nbDim
            real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
            real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
            real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
         end subroutine Set_f_s
         subroutine Set_R_s(ndimL, nbDim, Qdof, w, Dw, Re_1, R_s)
            integer, intent(in) :: ndimL, nbDim, Qdof
            real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
            real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
            real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
            !real, intent(in) :: Re_1                     ! inverse of Reynolds number
            real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
         end subroutine Set_R_s
      end interface

      real, dimension(:,:,:), allocatable :: R_s
      real, dimension(:,:,:), allocatable :: f_s
      integer :: Qdof,d

      Qdof = V_rule%Qdof                   !number of volume integration nodes

      !TODO - shouldnt be %Re before diffusive flux?
      ! TODO FERROR- it call Eval_R_s_Elem which works with the solution w(0, -- does it work right with STDG ???

      allocate(R_s(1:Qdof, 1:nbDim, 1:ndim) )
      call Eval_R_s_Elem(Set_R_s, elem, R_s(1:Qdof, 1:nbDim, 1:ndim) )

      allocate(f_s(1:Qdof, 1:nbDim, 1:ndim) )
      call Eval_f_s_Elem(Set_f_s, elem, f_s(1:Qdof, 1:nbDim, 1:ndim) )

!      print*, 'R_S: ' , R_s(1, 1:nbDim, 1:ndim)
!      print*, 'F_S: ' , f_s(1, 1:nbDim, 1:ndim)

      do d =1,nbDim
         flux(1:ndim, d, 1:Qdof) = &
            transpose( f_s(1:Qdof, d, 1:ndim) - R_s(1:Qdof, d, 1:ndim) )
      end do !d

   end function EvalRTNFluxVolume

   ! work only for scalar eq now!!!
   ! WARNING %time%ctime must be set before and wST -> w has to be done also in the right time integ node
   function EvalRTNFluxEdgeIPG( Set_K_sk, elem, iedge, Qdof ) result (flux)
      class(element), intent(in) :: elem
      integer, intent(in) :: iedge
      integer, intent(in) :: Qdof

      real, dimension(1:ndim, 1:nbDim, 1:Qdof) :: flux
      interface
      subroutine Set_K_sk(ndimL, nbDim,Qdof, w, Dw, Re_1, K_sk)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,ndimL), intent(inout) :: K_sk
       end subroutine Set_K_sk
      end interface
      real, dimension(:), allocatable :: Re1
      real, dimension(:,:), allocatable :: jump_wi, wi
      real, dimension(:,:,:), allocatable :: Dwi
      real, dimension(:,:,:,:,:), allocatable :: Ksk
      integer :: Gnum, Gdof, j
      real :: omega_edge

      Gnum = elem%face(fGnum, iedge)

      associate ( G_rule => state%space%G_rule( Gnum ) ) !choice of edge quadrature
         Gdof = G_rule%Qdof !number of edge integ nodes

         !control - on curved edges the normal is saved only in the specific integ nodes
         if (Qdof /= Gdof) &
            stop'Problem in EvalRTNFluxEdgeIPG Qdof /= %face(fGnum,iedge) '

         if(elem%face(neigh, iedge) > 0) then
            omega_edge  = 0.5
         else
            omega_edge  = 1.0
         endif

         ! jump on edge
         allocate( jump_wi(1:Gdof,1:ndim) )
         call ElementEdgeJump(elem, iedge, jump_wi(1:Gdof,1:ndim) )

         ! matrix K on the edge
         allocate( wi(1:Gdof,1:ndim) )
         call Eval_w_Edge(elem, iedge, wi, .false.) ! false - from the inside

         allocate( Dwi(1:Gdof, 1:ndim, 1:nbDim) )
         call Eval_Dw_Edge(elem, iedge, Dwi(1:Gdof, 1:ndim, 1:nbDim), .false.) ! false - from the inside

         allocate(Re1(1:Gdof))
         Re1(1:Gdof) = state%model%Re1

         allocate(Ksk(1:Gdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim) )
         call Set_K_sk( ndim, nbDim, Gdof, wi(1:Gdof,1:ndim), &
            Dwi(1:Gdof, 1:ndim, 1:nbDim), Re1(1:Gdof) , &
            Ksk(1:Gdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim) )


         write(debug,*) ' EvalRTNFluxEdgeIPG - what to do with Ksk for ndim > 1 '
         write(debug,*) 'EvalRTNFluxEdgeIPG - Watch if the transposition is right?'

         if (ndim > 1) then
            stop 'Problem in EvalRTNFluxEdgeIPG, implemented for ndim==1 only'
         else
            ! curved edges :
            if(elem%face(neigh,iedge) > 0 .or. elem%F%iFlin) then ! inner or straight boundary  edge
               do j = 1, Gdof
                  ! Warning the size of elem%n equals elem%dn and NOT 1
                  flux(1,1:nbDim, j) =  ( 1 / elem%dn(iedge) ) * omega_edge * jump_wi(j, 1) * &
                     matmul( elem%n(iedge, 1:nbDim) , Ksk(j, 1:nbDim, 1:nbDim, 1, 1) )
               end do !j
            else !curved edge
               do j = 1, Gdof
                  flux(1,1:nbDim, j) =  ( 1 / elem%dnc(j) ) * omega_edge * jump_wi(j, 1) * &
                     matmul( elem%nc(j, 1:nbDim) , Ksk(j, 1:nbDim, 1:nbDim, 1, 1) )
                  stop 'EvalRTNFluxEdgeIPG not tested for curved edges'
               end do !j
            end if
         endif

      end associate

   end function EvalRTNFluxEdgeIPG

   ! compute RTN fluxes sigma in the time integ nodes
   subroutine SetRTNstFluxesScalar( elem, tQnum, Qdof, sigma, sigmaIPG )
      class(element), intent(inout) :: elem
      integer, intent(in) :: tQnum
      integer, intent(in) :: Qdof ! should be max of elem%Qdof, G_rule(elem%face(fGnum,:))%Qdof
      real , dimension(0:3, 1:ndim, 1:nbDim, 1:Qdof, 1:tQnum), intent(out) :: sigma ! 0:numberOfEdges(=3),0 for volume

      real , dimension(1:3, 1:ndim, 1:nbDim, 1:Qdof, 1:tQnum), intent(out), optional :: sigmaIPG ! for SIPG/NIPG only
      integer :: alpha, Gnum, Gdof, iedge
      real :: cTime
      !real, dimension( 1:Qdof, 1:nbDim, 1:ndim) :: help

      cTime = state%time%ctime

     ! print*, 'First CTIME:' , cTime
     ! state%time%ctime = state%time%ctime - state%time%tau(1)
     ! print*, '%cTime after - tau:', state%time%ctime


      !SIPG/NIPG
      if ( present( sigmaIPG ) ) then

         do alpha = 1, tQnum
            !set the right %ctime - boundary elements need to take value from the exact sol
            state%time%ctime =  (state%time%ttime - state%time%tau(1) ) + &
               state%time%tau(1) * state%time%T_rule(tQnum)%lambda(alpha)

            !to computation of fluxes we need the solution w in time integ nodes for all neighboring elements
            ! TODO: we also need to switch %ctime correctly -> boundary elements values from the outside are taken from the exact solution
            call Transfer_wST_to_w_Patch(elem, alpha, tQnum)

            !compute the flux \sigma in current time
            !edges
            do iedge= 1,3
               Gnum = elem%face(fGnum,iedge)
               associate( G_rule => state%space%G_rule(Gnum) )
                  !eval \hat\sigma on the edge (multiplied by the outer normal)
                  Gdof = G_rule%Qdof
                   sigma( iedge, 1:ndim, 1, 1:Gdof, alpha) = EvalRTNFluxEdgeNormal( &
                     Set_Ppm_scalar, Set_R_s_scalar, elem, iedge, Gdof)

            !!!!!!! another part to the volume flux sigma for NIPG/SIPG!!!!!!!!!!!!!!
                  sigmaIPG(iedge, 1:ndim, 1:nbDim, 1:Gdof, alpha ) = &
                     EvalRTNFluxEdgeIPG( Set_K_sk_scalar, elem, iedge, Gdof)
               end associate

            end do !iedge

            !volume
            sigma(0, 1:ndim, 1:nbDim, 1:Qdof, alpha) = EvalRTNFluxVolume( &
                Set_f_s_scalar, Set_R_s_scalar, elem, state%space%V_rule( elem%Qnum ) )

         end do !alpha
      !IIPG - no sigmaIPG
      else
         !print*, ' SetRTNstFluxesScalar ' , tQNum

         do alpha = 1, tQnum
            !set the right %ctime - boundary elements need to take value from the exact sol
            state%time%ctime =  (state%time%ttime - state%time%tau(1) ) + &
               state%time%tau(1) * state%time%T_rule(tQnum)%lambda(alpha)

            !print*, 'cTime - ' , state%time%ctime
            !to computation of fluxes we need the solution w in time integ nodes for all neighboring elements
            ! TODO: we also need to switch %ctime correctly -> boundary elements values from the outside are taken from the exact solution
            call Transfer_wST_to_w_Patch(elem, alpha, tQnum)



            !print*, 'elem%W', elem%w(0,:)
            !print*, 'WST:', elem%wST(:,:,:)

            !compute the flux \sigma in current time
            !edges
            do iedge= 1,3
               Gnum = elem%face(fGnum,iedge)
               associate( G_rule => state%space%G_rule(Gnum) )
                  !eval \hat\sigma on the edge (multiplied by the outer normal)
                  Gdof = G_rule%Qdof
                  sigma( iedge, 1:ndim, 1, 1:Gdof, alpha ) = EvalRTNFluxEdgeNormal( &
                     Set_Ppm_scalar, Set_R_s_scalar, elem, iedge, Gdof)
               end associate
            end do !iedge

            !volume fluxes
            sigma(0, 1:ndim, 1:nbDim, 1:Qdof, alpha ) = EvalRTNFluxVolume( &
                Set_f_s_scalar, Set_R_s_scalar, elem, state%space%V_rule( elem%Qnum ) )



!            help(1:Qdof, 1:nbDim, 1:ndim) = EvalRTNFluxVolume( &
!                Set_f_s_scalar, Set_R_s_scalar, elem, state%space%V_rule( elem%Qnum ) )
!
!                        print*, nbDim, Qdof, ndim
!            print*, 'HELP: ', help(1:Qdof,1,1)
!            print*, '-------------'
!            print*, 'HELP: ', help(1:Qdof,2,1)
!            print*, '-------------'
!
!                           print*, 'HELP: ', help(1:Qdof,1,1)
!            print*, '-------------'
!            print*, 'HELP: ', help(1:Qdof,2,1)
!            print*, '-------------'
!


         end do !alpha


      endif

!      print*,'SetRTNstFluxesScalar:', nbDim, Qdof, ndim
!            print*, 'SIGMA1: ', sigmaIPG(1,1,1,:,1)
!            print*, '-------------'
!            print*, 'SIGMA2: ', sigmaIPG(1,1,2,:,1)
!            print*, '-------------'
!            print*, 'SIGMA2: ', sigma(0:3,1,2,1,:)
!            print*, '-------------'

      state%time%ctime = cTime
     !stop 'END SetRTNstFluxesScalar'

   end subroutine SetRTNstFluxesScalar

   !> eval the advection-diffusion flux from the equation
   !> \f$ \sigma(u,\nabla u) = -{\bf K}(u)\nabla u + \vec \phi(u) \f$
   function EvalSTsigma( elem, tQnum) result (sigma)
      class(element), intent(inout) :: elem
      integer, intent(in) :: tQnum ! #number of T_Rule nodes
      real, dimension(1:ndim,1:nbDim, 1:elem%Qdof, 1:tQnum) :: sigma

      integer :: alpha
      real :: cTime

      cTime = state%time%ctime

      do alpha = 1, tQnum

            !set the right %ctime - boundary elements need to take value from the exact sol
            state%time%ctime =  (state%time%ttime - state%time%tau(1) ) + &
               state%time%tau(1) * state%time%T_rule(tQnum)%lambda(alpha)

            call Transfer_wST_to_w_Elem(elem , alpha, tQnum)

            !compute the volume flux \sigma in current time
            sigma(1:ndim, 1:nbDim, 1:elem%Qdof, alpha) = EvalRTNFluxVolume( &
                Set_f_s_scalar, Set_R_s_scalar, elem, state%space%V_rule( elem%Qnum ) )

      end do !alpha

      state%time%ctime = cTime

   end function EvalSTsigma

   function SetRTNstSpaceRhs( elem, tQnum, Fdeg, Fdof, Mdof, sigma, sigmaIPG) result (space_rhs)
      class(element), intent(in) :: elem
      integer, intent(in) :: tQnum
      integer, intent(in) :: Fdeg, Fdof
      integer, intent(in) :: Mdof !  max ( elem%Qdof , state%space%G_rule(MGnum)%Qdof )
      real, dimension( 0:3, 1:ndim, 1:nbDim, 1:Mdof, 1:tQnum ), intent(in) :: sigma
      real, dimension( 1:3, 1:ndim, 1:nbDim, 1:Mdof, 1:tQnum ), intent(in), optional :: sigmaIPG   ! for NIPG/SIPG only

      real, dimension(1:ndim,1:Fdof, 1:tQnum) :: space_rhs

      real, dimension(1:Mdof, 1:nbDim) :: qi
      real, dimension(:,:,:), allocatable :: x ! for sigmaIPG
      integer :: i, j, d
      integer :: iTime, iedge, Gdof, Gnum, Qdof
      real :: temp
      !real, dimension(:,:,:), allocatable :: help

      Qdof = elem%Qdof

      i = 0
!!!!!!!eval edge test functions and space part of RHS
      do iedge= 1, 3
         Gnum = elem%face(fGnum,iedge)
         ! number of Quadrature nodes
         Gdof = elem%face(fGdof,iedge)
         associate( G_rule => state%space%G_rule(Gnum) )

         do j = 1,Fdeg + 1 ! number of edge test functions
            i = i+1
            call EvalRTNMomentumEdgeTestFunction( iedge, j, G_rule, qi( 1:Gdof, 1) )

            do iTime = 1, tQnum
               !compute space_rhs(iTime, i)

               do d = 1,ndim
                  space_rhs(d, i, iTime) = EvalScalarProdEdge(elem, iedge,  &
                     sigma(iedge, d, 1, 1:Gdof, iTime), qi(1:Gdof, 1) )
               end do !d
            end do !iTime
         end do !j

         end associate
      enddo

      allocate(x(0:3, 1:Mdof, 1:nbDim))
         !--------- ^^^  0= volume, 1,2,3= 1st, 2nd, 3rd face ! setting integ nodes on edges
      call Set_RTN_integ_nodeNew( elem, x(0:3, 1:Mdof, 1:nbDim) )




!!!!!!!eval volume test functions and space part of RHS
      do j=1, Fdeg*(Fdeg+1)  ! Fdeg*(Fdeg+1) degrees of freedom over element
         i = i+1
         ! eval test function volume
         call EvalMomentumVolumeD( Fdeg, i, Qdof, x(0,1:Qdof, 1:nbDim), &
         qi(1:Qdof, 1:nbDim) )

         !compute i-th of the space part of RHS
         do iTime = 1, tQnum
            ! eval space L2 scalar product
            do d=1, ndim
               call IntegrateFunctionsVec( elem, sigma(0, d, 1:nbDim, 1:Qdof, iTime ) , &
                  qi(1:Qdof, 1:nbDim), space_rhs(d, i, iTime) )

               !print*, 'SIGMA:', sigma(0, d, 1:nbDim, 1:Qdof, iTime )
            end do !d
         end do !iTime
      end do !j

!!!!!!! in SIPG/NIPG case we need also test functions qi on edges
      if ( present(sigmaIPG) ) then
         !TODO: SMAZ
         !allocate( help(1:ndim,1:3*Fdeg*(Fdeg+1),1:tQnum) )
         !control
         if ( state%space%m_IPG == 0 ) &
            stop 'SetRTNstSpaceRhs: There is no need to compute sigmaIPG for IIPG'

         do iedge = 1,3
            i = i - Fdeg*(Fdeg+1)
            Gnum = elem%face(fGnum,iedge)
            ! number of Quadrature nodes
            Gdof = elem%face(fGdof,iedge)

            associate( G_rule => state%space%G_rule(Gnum) )
               do j = 1, Fdeg*(Fdeg+1)
                  i = i + 1

                  ! eval RTN test function in edge integ nodes
                  call EvalMomentumVolumeD(Fdeg, i, Gdof, x(iedge, 1:Gdof, 1:nbDim), &
                  qi(1:Gdof, 1:nbDim))

                  do iTime = 1, tQnum
                     do d = 1, ndim
                        call IntegrateFunctionsEdgeVec(elem, iedge, Gnum, Gdof, &
                        sigmaIPG(iedge, d, 1:nbDim, 1:Gdof, iTime), qi(1:Gdof, 1:nbDim), temp )
                        !THETA - SIPG=-1,NIPG=1,IIPG=0 (reversely than in Vohralik scripts)
                        !help(d, Fdeg*(Fdeg+1)*(iedge-1) +1, iTime) = temp
                        space_rhs(d, i, iTime) = space_rhs(d, i, iTime) - (state%space%m_IPG*temp)
                     end do !d
                  end do !iTime
               end do !j
            end associate

         end do !iedge

      !deallocate(help)

      else
         if (state%space%m_IPG /= 0 ) &
            stop 'SetRTNstSpaceRhs: sigmaIPG should be computed for S/NIPG'
      endif

      deallocate(x)

      !TODO: Remove control
      if ( i /= Fdof ) stop 'SetRTNstSpaceRhs: i /= Fdof'

   end function SetRTNstSpaceRhs

   ! computes divergence in integration nodes of the RTNst function given by the coeffs of st basis functions
   function evalRTNstFunDiv_Q( elem, Fdeg, Fdof, Tdof, tQnum, fun ) result( fun_div )
      class(element), intent(in) :: elem
      integer, intent(in) :: Fdeg, Fdof
      integer, intent(in) :: Tdof, tQnum
      real, dimension(1:ndim, 1:Fdof, 1:Tdof), intent(in) :: fun
      real, dimension(1:ndim, 1:elem%Qdof, 1:tQnum) :: fun_div

      integer :: iNode,d, Qdof
      real, dimension(1:ndim, 1:Fdof, 1:tQnum) :: ff
      real, dimension(1:Fdof, 1:elem%Qdof) :: divPsi

      Qdof = elem%Qdof

      associate ( loc_RTN => state%loc_RTN(Fdeg) )
         call Eval_Loc_RTN_Div_Elem( elem, loc_RTN, divPsi(1:Fdof, 1:Qdof) )
      end associate !loc_RTN

      do iNode = 1,tQnum

        ff(1:ndim, 1:Fdof, iNode) = evalSTfunInIntTimeDof( elem, ndim, Fdof, Tdof, &
            fun(1:ndim,1:Fdof,1:Tdof), iNode, tQnum )

         do d = 1, ndim
               fun_div(d, 1:Qdof, iNode) = matmul( ff(d, 1:Fdof, iNode), divPsi(1:Fdof, 1:Qdof) )
         end do !d

      end do ! iNode

   end function evalRTNstFunDiv_Q

  !> from the computed RTN momentums compute
  !> flux (from RTN) as  basis coeficients of LocRTN basis
  !> MM stores RTN momentums of LocRTN basis functions
  function SetRTNstFluxReconstr( Tdof, Fdof, MM, rhs) result ( flux )
    integer, intent(in) :: Tdof
    integer, intent(in) :: Fdof
    real, dimension(1:Fdof, 1:Fdof), intent(in) :: MM
    real, dimension(1:ndim, 1:Fdof, 1:Tdof), intent(in) :: rhs
    real, dimension(1:ndim, 1:Fdof, 1:Tdof) :: flux

    real, dimension(:,:), allocatable :: x

    integer :: d

    write(debug, *) 'Control SetRTNstFluxReconstr for ndim>1 !!!'


    flux( 1:ndim, 1:Fdof, 1:Tdof ) = rhs( 1:ndim, 1:Fdof, 1:Tdof )

    do d = 1,ndim
      call SolveLocalMatrixProblem(Fdof, MM(1:Fdof,1:Fdof), Tdof, flux( d, 1:Fdof, 1:Tdof) )
    end do !d

  end function SetRTNstFluxReconstr

  !> setting the volume and edge integ. nodes (barycentric)
  subroutine Set_RTN_integ_nodeNew( elem, xi)
    class(element), intent(in) :: elem
    real, dimension(0:3 ,1:max( elem%Qdof, maxval( elem%face(fGdof,:) ) ) ,1:nbDim ), intent(out) :: xi
    real, dimension(1:nbDim) :: x0, a0
    real :: t
    integer :: ie, l, Qdof

    !type(volume_rule), pointer :: V_rule
    !type(Gauss_rule), pointer :: G_rule

    associate( V_rule => state%space%V_rule(elem%Qnum) )
    ! volume integ. nodes
    Qdof = V_rule%Qdof

    !MGdof = max(G_rule(:)%Qdof)
!    if ( Mdof /= Qdof .and. Mdof /= MGdof ) &
!      stop 'Problem in Set_RTN_integ_nodeNew - # of nodes is different from Mdof'

    xi(:,:,:) = 0.0
    xi(0,1:Qdof, 1:2) = V_rule%lambda(1:Qdof,1:2)

    end associate

    ! face integ. nodes
    do ie = 1, 3   ! loop through edges

       if(ie == 1 ) then
          x0(1:2) = (/0., 0./)
          a0(1:2) = (/1., 0./)
       elseif(ie == 2 ) then
          x0(1:2) = (/1., 0./)
          a0(1:2) = (/-1., 1./)
       elseif(ie == 3) then
          x0(1:2) = (/0., 1./)
          a0(1:2) = (/0., -1./)
       endif

       associate ( G_rule => state%space%G_rule( elem%face(fGnum, ie) ) )

       do l=1, G_rule%Qdof  !loop though 1D Gauss quadrature rule
          t = G_rule%lambda(l)
          xi(ie,l,1:2) = x0(1:2) +  a0(1:2) * t
       enddo ! l

       end associate
    enddo ! ie

  end subroutine Set_RTN_integ_nodeNew

   !> computes the Radau reconstruction $ R = u_{h, \tau} - \jump{u_{h,\tau} r_h} $, where r_h is the Radau polynomial of degree elem%tdeg+1
   function evalRadauReconstruct( elem ) result ( reconstr )
      class( element ), intent(in) :: elem
      real, dimension(1:ndim, 1:elem%dof, 1:elem%Tdof+1) :: reconstr ! one higher degree than wST

      integer :: i, dof, Tdof
      class( Time_rule ), pointer :: T_rule
      real, dimension( 1:ndim, 1:elem%dof ) :: wi

      dof = elem%dof
      Tdof = elem%Tdof
      T_rule => state%time%T_rule( Tdof )

      if ( .not. allocated(T_rule%RadauPol) ) &
         stop 'Radau pol is not allocated in evalRadauReconstruct'

      do i = 1, ndim
         ! eval the solution in the time moment t_{m-1}^+, phi(i,-1) = phi_i(0)
         wi(i, 1:dof) = matmul( elem%wST( i, 1:elem%dof, 1:elem%Tdof) , T_rule%phi(1:elem%Tdof,-1) )
      end do

      !eval the jump $\{ u_h(t_{m-1})\}$
      wi(1:ndim, 1:dof) = wi(1:ndim, 1:dof) - elem%wSTfinAD( 1:ndim, 1:dof )

      ! solution wST \in P^q => the last coefficient is zero

      ! the jump of wST multiplied by the Radau polynomial
      do i =1,  Tdof + 1
         reconstr( 1:ndim, 1:dof,i ) =  T_rule%RadauPol(i) * wi(1:ndim,1:dof)
      end do

      reconstr( 1:ndim,1:dof,1:Tdof ) = elem%wST(1:ndim,1:dof,1:Tdof) - reconstr( 1:ndim,1:dof,1:Tdof)
      reconstr( 1:ndim,1:dof, Tdof+1 ) = - reconstr( 1:ndim,1:dof, Tdof+1 )

   end function evalRadauReconstruct


end module eval_rtn_st_mod
