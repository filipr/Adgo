!> aposteriori error estimation subroutines for Vohralik approach for nonlinear problem
!> estimates in the dual (residual) form using RTN flux reconstruction
module dual_estim
  use rav_tho_nedelec
  use eval_sol
  use apost_estimation
  use euler_problem
  use loc_rav_tho_ned
  use eval_rav_tho_ned
  use model_oper
   use eval_jumps

  implicit none

  public:: ComputeDualEstim         !
  public:: ComputeElemDualEstim     !
  public:: ConstructFluxCD          !?
  public:: ResidualElemEstimCD      !
  public:: FluxElemEstimCD          !
  public:: ComputeLocRTNMomentumsElem
  public:: FluxMM2RTN

  !public:: CheckRTNMomentums   !! should be created based on ConstructFluxCD

contains
  !> estimates the computational error using the approach based on the flux reconstruction
  subroutine ComputeDualEstim( )
    class(element), pointer :: elem
    integer :: dof, Fdeg, Fdof
    integer :: i
    real, dimension(:,:), allocatable :: estim, sum_estim

    real, dimension(:,:), pointer :: phi, phiT   !for TEST
    real, dimension(:,:), allocatable :: wi, q !for TEST
    real :: val
    integer :: k, Qdof, ist, l !for TEST

    !print*,'##############  dual_estim.f90 ',state%time%iter

    allocate(estim(1:max_eta, 1:ndim), sum_estim(1:max_eta, 1:ndim) )
    sum_estim(:,:) = 0.

    !if(state%time%iter == 1) state%err(IC_L2) = 0.

    do i = 1, grid%nelem
       elem => grid%elem(i)
       call ComputeElemDualEstim(elem, estim(1:max_eta, 1:ndim) )

       if(state%time%iter > 0) then
          sum_estim( DFn, 1:ndim) = sum_estim( DFn, 1:ndim) + estim( DFn, 1:ndim)
          sum_estim(DFnS, 1:ndim) = sum_estim(DFnS, 1:ndim) + estim(DFnS, 1:ndim)*2
          sum_estim(DFnT, 1:ndim) = sum_estim(DFnT, 1:ndim) + estim(DFnT, 1:ndim)
          sum_estim(  Rn, 1:ndim) = sum_estim(  Rn, 1:ndim) + estim(  Rn, 1:ndim)
          sum_estim(NC1n, 1:ndim) = sum_estim(NC1n, 1:ndim) + estim(NC1n, 1:ndim)

          sum_estim(DFRn, 1:ndim) = sum_estim(DFRn, 1:ndim) &
               + ( estim( DFn, 1:ndim)**0.5 + estim( Rn, 1:ndim)**0.5 )**2

          sum_estim(IC, 1:ndim) = sum_estim(IC, 1:ndim) + estim(IC, 1:ndim)
       endif
    enddo  ! end of i=1,grid%nelem

    !if(state%time%iter == 1) state%err(IC_L2) = state%err(IC_L2)**0.5

    state%estim(1:IC, 1:ndim) = state%estim(1:IC, 1:ndim) + sum_estim(1:IC, 1:ndim)

    state%estim(total, 1:ndim) = state%estim(DFRn, 1:ndim)**0.5 &
         + state%estim(NC1n, 1:ndim)**0.5 + state%estim(IC, 1:ndim)**0.5


    !write(100,'(a6,7es12.4)') 'sum?',state%estim( 1:IC, 1)**0.5, state%estim(total, 1:ndim)


    deallocate(estim, sum_estim)

  end subroutine ComputeDualEstim

  !> compute error estimates via the dual (residual) form using RTN flux reconstruction
  !> for one elements
  subroutine ComputeElemDualEstim(elem, estim )
    type(element), intent(inout) :: elem
    real, dimension(1:max_eta, 1:ndim), intent(inout) :: estim
    type(basis_rtn_fe), pointer:: loc_RTN
    real, dimension(:,:), allocatable :: MMinvRE ! inverse of Momentum matrix on elem
    !!real, dimension(:,:,:), allocatable :: flux  ! flux in integ nodes
    real, dimension(:,:), allocatable :: fluxi   ! basis coefficients of flux
    real, dimension(:), allocatable :: fluxmomentum  !for testing
    !real, dimension(:,:,:), allocatable :: CDfluxes   ! physical fluxes
    integer :: Fdeg, Fdof, ib



    ! scaling factor computed in errorFlux.f90
    !elem%CTn = state%time%FinTime/state%time%tau(1)**2 + (elem%CK + elem%Cb)/elem%diam**2

    !Fdeg = state%space%deg-1
    Fdeg = state%space%deg

    !Fdof = state%RTN(Fdeg)%dof
    Fdof = SetRTNdof(Fdeg)

    loc_RTN => state%loc_RTN(Fdeg)
    if(.not. loc_RTN%defined ) call Init_Loc_RTN(state%loc_RTN(Fdeg), Fdeg )

    allocate(MMinvRE(1:Fdof, 1:Fdof) )
    allocate(fluxi(1:ndim, 1:Fdof) )
!    allocate(CDfluxes(1:elem%Qdof, 1:nbDim, 1:ndim)  )

    ! evaluate the momentums of the local RTN basis on K
    !call ComputeRTNMomentumsRealElem(elem, Fdeg, Fdof, MMinvRE)  ! RE = real element
    call ComputeLocRTNMomentumsElem(elem, Fdeg, Fdof, MMinvRE)  ! RE = real element

    ! flux at k-th time level
    if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then

       ! allocation of arays for storing of RTN fluxes, RTNflux(l,*,*),at t_{k-l} level
       if(state%time%iter == 0) then
          allocate(elem%RTNflux(0:1, 1:ndim, 1:Fdof) )

          state%time%ctime = 0.
          call ConstructFluxCD(Set_f_s_scalar, Set_Ppm_scalar, Set_R_s_scalar,Set_K_sk_scalar,&
               elem, Fdeg, Fdof, MMinvRE, elem%RTNflux(0, 1:ndim, 1:Fdof) )

          deallocate(MMinvRE, fluxi )
          return

       else
          ! storing of the older values
          elem%RTNflux(1, 1:ndim, 1:Fdof) = elem%RTNflux(0, 1:ndim, 1:Fdof)

          state%time%ctime = state%time%ttime
          ! old version
          !call ConstructFluxCD(Set_f_s_scalar, Set_Ppm_scalar, Set_R_s_scalar,Set_K_sk_scalar,&
          !     elem, Fdeg, Fdof, MMinvRE, fluxi(1:ndim, 1:Fdof))

          ! new version
          call ConstructFluxCD(Set_f_s_scalar, Set_Ppm_scalar, Set_R_s_scalar,Set_K_sk_scalar,&
               elem, Fdeg, Fdof, MMinvRE, elem%RTNflux(0, 1:ndim, 1:Fdof) )

          ! only for TEST !!!  OLD VERSION, PROBABLY BAD
          !allocate(fluxmomentum(1:Fdof) )
          !call CheckRTNMomentums(elem, Fdeg, Fdof, fluxi, fluxmomentum)
          !deallocate(fluxmomentum )

          call ResidualElemEstimCD(elem, Fdeg, Fdof, fluxi, estim(Rn, 1:ndim))

          call FluxElemEstimCD(Set_f_s_scalar, Set_R_s_scalar, &
               elem, Fdeg, Fdof, fluxi,  estim(DFn:DfnT, 1:ndim))
       endif

    else
       print*,'Case nbDim ==',nbDim,' & ndim ==',ndim,' not implemented in dual_estim.f90'
       stop
    endif

    ! computed in errorFlux.f90
    !call NonConfEstimCD(elem, estim(NC1n,1:ndim) )
    !estim(NC1n,1:ndim) = elem%eta(NC1n, 1:ndim)
    estim(NC1n,1:ndim) = 0.   ! completly done in errorF.f90

    ! we sum the same values for each time step
    estim(IC,1:ndim) = elem%eta(IC,1)**2 / (elem%CTn*state%time%tau(1) )

    ! \eta_space include {\not nonconformity term} and redidual term!!!
    !estim(DFnS,1:ndim) = ( estim(DFnS,1:ndim)**0.5 + estim(Rn, 1:ndim)**0.5 &
    !     + elem%eta(NC1n, 1:ndim)**0.5 )**2
    estim(DFnS,1:ndim) = ( estim(DFnS,1:ndim)**0.5 + estim(Rn, 1:ndim)**0.5 )**2

    ! for adaptation
    elem%eta(DFn, 1:ndim)  = estim(DFn, 1:ndim)
    elem%eta(DFnS, 1:ndim)  = estim(DFnS, 1:ndim)
    elem%eta(DFnT, 1:ndim)  = estim(DFnT, 1:ndim)
    elem%eta(Rn, 1:ndim)  = estim(Rn, 1:ndim)

    elem%eta(resST, 1) = elem%eta(DFn, 1)**0.5 + elem%eta(Rn, 1)**0.5  + elem%eta(NC1n, 1)**0.5

    write(*,'(a6,i5,5es12.4)') &
        '######',elem%i, elem%eta(resST, 1), elem%eta(DFn, 1), elem%eta(Rn, 1), elem%eta(NC1n, 1)

    !write(1000+state%time%iter,*) elem%xc(:),estim(DFn, 1:ndim)**0.5, estim(Rn, 1:ndim)**0.5
    !write(*,'(a12,i5,6es12.4)') 'eta DF, R :', elem%i, &
    !     estim(DFn, 1:ndim)**0.5, estim(Rn, 1:ndim)**0.5, estim(NC1n, 1:ndim)**0.5

    deallocate(MMinvRE, fluxi ) !,CDfluxes )

  end subroutine ComputeElemDualEstim

  !> evaluate noncoformity estimator \f$ \int_{\partial K} [u_h]^s\, d s\f$
  subroutine NonConfEstimCD(elem, etaNC1n )
    type(element), intent(in) :: elem
    real, dimension(1:ndim), intent(inout) :: etaNC1n
    real, dimension(1:ndim) :: jumps
    integer:: ie

    etaNC1n(1:ndim) = 0.
    do ie = 1, elem%flen

       ! jumps of the solution
       call IntegElementEdgeJump(elem, ie, jumps(1:ndim) )
       etaNC1n(1:ndim) = etaNC1n(1:ndim) + jumps(1:ndim)
    enddo

  end subroutine NonConfEstimCD



  !> recontruction of the flux at element using the RTN momentums for
  !> general nonlinear convection-diffusion equation
  subroutine ConstructFluxCD(Set_f_s, Set_Ppm, Set_R_s, Set_K_sk, elem, Fdeg, Fdof, MMinvRE,&
       fluxi) !, CDfluxes)
    interface
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x, ie )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
        integer, intent(in) :: ie
       end subroutine Set_f_s
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
      subroutine Set_R_s(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
      end subroutine Set_R_s
      subroutine Set_K_sk(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, K_sk, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,ndimL), intent(inout) :: K_sk
       end subroutine Set_K_sk
    end interface

    type(element), intent(inout) :: elem
    integer, intent(in) :: Fdeg, Fdof               ! degree and DOF of the RTN element
    real, dimension(1:Fdof, 1:Fdof), intent(in) :: MMinvRE ! (Momentum matrix)^-1 on elem
    real, dimension(1:ndim, 1:Fdof), intent(inout) ::  fluxi
    real, dimension(:, :), allocatable :: flux_momentums
    real, dimension(:,:,:), allocatable :: x, phiRE, CDfluxes
    real, dimension(:,:), allocatable :: Rflux, Cflux
    real, dimension(:,:,:), allocatable :: R_s, f_s
    real, dimension(:), allocatable :: func
    real, dimension(:,:,:), allocatable :: func3
    real, dimension(:,:), allocatable :: jump_wi, wi, qi
    real, dimension(:,:,:,:,:), allocatable :: Ksk
    class(element), pointer :: elem1

    type(volume_rule), pointer :: V_rule
    type(Gauss_rule), pointer :: G_rule
    integer ::  Qdof, Gdof, Gdeg
    integer :: ie, i, k, j1, l, k1,  ib, ibnd
    integer :: HGidx, Tface, Tface1, Tface2
    integer, dimension(:), allocatable :: HGvertex
    real :: val, omega_edge, edge_momentum
    real, dimension(:,:), allocatable :: xi, Fx
    real, dimension(:,:), allocatable :: Re1



    write(debug,*) 'Warning in dual_estim 3 dimensional array func was previously 1 dim'
    write(debug,*)  'Many problems with calling Set_K_sk(ndim, Gdof, wi(1:Gdof,1:ndim), func3, Re1 , Ksk) - wrong arguments'


    ib = 0                               !index of momentum

    allocate(flux_momentums(1:ndim, 1:Fdof) )
    allocate( CDfluxes(1:elem%Qdof, 1:nbDim, 1:ndim) )

    flux_momentums(:, :) = 0

    if(elem%HGnode) then
       if(elem%type /= 3) print*,' TROUBLES in dual_estim.f90 with HG'

       allocate( HGvertex(4) )
       HGvertex(1:3)   = elem%HGvertex(1:3 )
       HGvertex(4) = elem%flen + 1

       !write(*,'(i5,a2,10i5)') elem%i,'|',elem%HGvertex(:)
       !write(*,'(i5,a2,10i5)') elem%i,'|',HGvertex(:)
       !write(*,'(i5,a2,10i5)') elem%i,'|',elem%HGface(1,:)
       !write(*,'(i5,a2,10i5)') elem%i,'|',elem%HGface(2,:)
       !print*,'-------------------'
    endif


    ! face momentums of flux
    do ie = 1, 3!loop through edges of triangle  (meshes with HG nodes  included)

       if(.not. elem%HGnode) then
          Tface1 = ie   ! only one element edge on the triangle edge
          Tface2 = ie
          HGidx = 1
       else
          Tface1 = HGvertex(ie)
          Tface2 = HGvertex(ie+1) - 1
       endif

       !flux_momentums(:, :) = 0

       ib = (ie - 1) * (Fdeg + 1) ! index of the momentums for HG nodes
       do i=1, Fdeg+1    ! Fdeg+1 degrees of freedom over face
          ib = ib + 1    ! ib is index of the momentum

          !flux_momentums(:, ib) = 0

       do Tface = Tface1, Tface2
          if( elem%HGnode) HGidx = elem%HGface(2,Tface)

          edge_momentum = 0

          !if ( elem%HGnode) then
             !write(*,'(a6,i5)') 'Fdeg:', Fdeg
             !write(*,'(a4,7i5)') 'HG:',elem%i, elem%flen, ie, Tface, Tface1, Tface2, HGidx
             !write(*,'(a10,6i5)') 'HGvertex:', HGvertex(:)
             !pause
          !endif

          Gdeg = elem%face(fGnum, Tface)
          G_rule => state%space%G_rule(Gdeg) !choice of edge quadrature
          Gdof = G_rule%Qdof !number of edge integration nodes

          ! diffusive fluxes
          allocate(Rflux(1:Gdof, 1:ndim) )
          call Eval_aver_R_s_Edge(Set_R_s, elem, Tface, Rflux(1:Gdof, 1:ndim) )

          if (ib == 0) then
          if(elem%i == 1 .and. Tface == 2) &
               write(*,'(a6,2i5,20es14.6)') &
               'Rflux=',elem%i,elem%face(neigh,Tface), -Rflux(1:Gdof, 1:ndim)
          if(elem%i == 1 .and. Tface == 2) write(*,'(a6,i5)') 'ib', ib

          if(elem%i == 6 .and. Tface == 3) &
               write(*,'(a6,2i5,20es14.6)') &
               'Rflux=',elem%i,elem%face(neigh,Tface), -Rflux(1:Gdof, 1:ndim)
          if(elem%i == 6 .and. Tface == 3) write(*,'(a6,i5)') 'ib', ib

          endif


          ! convective fluxes
          allocate(Cflux(1:Gdof, 1:ndim) )
          Cflux(1:Gdof, 1:ndim) = 0.
          call Eval_NFlux_Edge(Set_Ppm, elem, Tface, Cflux(1:Gdof, 1:ndim) )



          ! jumps of the solutions
          allocate(jump_wi(1:Gdof,1:ndim))
          call ElementEdgeJump(elem, Tface, jump_wi(1:Gdof,1:ndim) )

          if (ib == 0) then
          if(elem%i == 1 .and. Tface == 2) &
               write(*,'(a6,2i5,20es14.6)') 'jump=',elem%i,elem%face(neigh,Tface), &
               jump_wi(1:Gdof, 1:ndim), state%space%sigma * state%model%Re1
          if(elem%i == 6 .and. Tface == 3) &
               write(*,'(a6,2i5,20es14.6)') 'jump=',elem%i,elem%face(neigh,Tface), &
               jump_wi(1:Gdof, 1:ndim), state%space%sigma * state%model%Re1

          write(*,*) 'attention in dual_estims with C_W !!!!!!!!!!!!!!!!!!!!!'

          endif

          ! total flux, already premultiplied by dn
          Rflux(1:Gdof, 1:ndim) = -Rflux(1:Gdof, 1:ndim) + Cflux(1:Gdof, 1:ndim) &
               + state%space%sigma * state%model%Re1 * jump_wi(1:Gdof, 1:ndim)  ! /dn(Gdof) * dn(Gdof) = 1

          ! test function for the evaluation of the momentum
          allocate( qi(1:Gdof, 1:nbDim  ) )

          !ib = (ie - 1) * (Fdeg + 1) ! index of the momentums for HG nodes
          !do i=1, Fdeg+1    ! Fdeg+1 degrees of freedom over face
           !  ib = ib + 1    ! ib is index of the momentum

             !flux_momentums(:, ib) = 0   !have to be before the Tface cycle!
             ! evaluate test function q in integ nodes

             !if(.not. elem%HGnode) then
             !   call EvalMomentumFace(ie, Fdeg, ib, Gdof, G_rule, qi(1:Gdof, 1))
                !write(*,'(a4,i5,20es12.4)') 'qiA:',HGidx, qi(1:Gdof, 1)
             !endif
             call EvalMomentumFaceHG(ie, HGidx, Fdeg, ib, Gdof, G_rule, qi(1:Gdof, 1))

             !if (ib == 6) then
             !if(elem%i == 1 .and. Tface == 2) &
             !write(*,'(a4,i5,20es12.4)') 'qiB:',HGidx, qi(1:Gdof, 1)
             !if(elem%i == 6 .and. Tface == 3) &
             !write(*,'(a4,i5,20es12.4)') 'qiB:',HGidx, qi(1:Gdof, 1)
             !endif

             qi(1:Gdof,1) = G_rule%weights(1:Gdof) * qi(1:Gdof,1)

             if (ib == 0) then
             if(elem%i == 1 .and. Tface == 2 .and. ib == 4) &
             write(*,'(a4,i5,20es12.4)') 'qiB:',HGidx, qi(1:Gdof, 1)
             if(elem%i == 6 .and. Tface == 3 .and. ib == 7) &
             write(*,'(a4,i5,20es12.4)') 'qiB:',HGidx, qi(1:Gdof, 1)


             if(elem%i == 1 .and. Tface == 2 .and. ib == 4) &
               write(*,'(a20,20es14.6)') 'Rflux(1:Gdof, 1:ndim)=', Rflux(1:Gdof, 1:ndim)
             if(elem%i == 6 .and. Tface == 3 .and. ib == 7) &
               write(*,'(a20,20es14.6)') 'Rflux(1:Gdof, 1:ndim)=', Rflux(1:Gdof, 1:ndim)
             endif



             do k=1, ndim
                edge_momentum = dot_product( qi(1:Gdof,1),  Rflux(1:Gdof, k) )

                if (ib == 0) then
                if(elem%i == 1 .and. Tface == 2 .and. ib == 4) &
                  write(*,'(a25, i2, es14.6)') 'elem%i, edge_momentum', elem%i, edge_momentum
                if(elem%i == 6 .and. Tface == 3 .and. ib == 7) &
                  write(*,'(a25, i2, es14.6)') 'elem%i, edge_momentum', elem%i, edge_momentum
                if(elem%i == 1 .and. Tface == 2 .and. ib == 4) pause
                if(elem%i == 6 .and. Tface == 3 .and. ib == 7) pause
                endif

                flux_momentums(k, ib) = flux_momentums(k, ib) + edge_momentum
                !write(*,'(a18,20es12.4)') 'flux_mome(1,ib) pred:', flux_momentums(k, ib)
                !flux_momentums(k, ib) = dot_product( qi(1:Gdof,1),  Rflux(1:Gdof, k) )
                !write(*,'(a18,20es12.4)') 'flux_mome(1,ib) po:', flux_momentums(k, ib)
             enddo

             !if(elem%i == 1 .or. elem%i == 6)
             !write(*,'(a28, i2, 20es12.4)') 'elem%i, flux_mome(1,:):', elem%i, flux_momentums(1, :)
             !if(elem%i == 1 .or. elem%i == 6) pause

          !enddo  !i

          deallocate(jump_wi)
          deallocate(Rflux)
          deallocate(Cflux)

          deallocate(qi)
       enddo !Tface
      enddo  !i
    enddo  ! ie

    !if(elem%i == 1 .or. elem%i == 6) &
    !write(*,'(a28, i2, 20es12.4)') 'elem%i, flux_mome(1,:):', elem%i, flux_momentums(1, :)
    !if(elem%i == 1 .or. elem%i == 6) pause


    !volumes momentums of flux
    V_rule => state%space%V_rule(elem%Qnum)    !choice of volume quadrature
    Qdof = V_rule%Qdof                   !number of volume integration nodes

    do i=1, Fdeg*(Fdeg+1)       ! deg*(deg+1) degrees of freedom over element
       ib = ib + 1
       flux_momentums(:, ib) = 0

       if(state%space%m_IPG /= 0.) then  ! SIPG or NIPG variants
          ! face integrals
          do Tface=1, 3 !loop through edges of triangle  (meshes with HG nodes - not included)

             G_rule => state%space%G_rule(elem%face(fGnum, Tface)) !choice of edge quadrature
             Gdof = state%space%G_rule(elem%face(fGnum, Tface))%Qdof !number of edge integ nodes

             if(elem%face(neigh, Tface) > 0) then
                omega_edge  = 0.5
             else
                omega_edge  = 1.0
             endif

             ! jump on edge
             allocate(jump_wi(1:Gdof,1:ndim))
             call ElementEdgeJump(elem, Tface, jump_wi)

             ! matrix K on the edge
             allocate(wi(1:Gdof,1:ndim))
             call Eval_w_Edge(elem, Tface, wi, .false.)
!             Eval_Dw_Edge v block

             allocate(func3(1:Gdof, 1:ndim, 1:nbDim) )
             ! TODO: FR func3 has to be a 3 dimension array! ,isnt func3 unused???

             ! should be Dw due to call of Set_K_sk
             func3(:,:,:) = state%model%Re1
             allocate(Re1(1:iRe, 1:Gdof))
             Re1(1:iRe, 1:Gdof) = state%model%Re1
             allocate(Ksk(1:Gdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim) )
             call Set_K_sk(ndim, nbDim, iRe, Gdof, wi(1:Gdof,1:ndim), func3, Re1 , &
                  Ksk(1:Gdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim), elem%xi(0, 1:Qdof, 1:nbDim) )
             print*,'Here in dual.estim.f90 it is bad, average on edge should by given',&
                  'probably ??'


             ! value of the test function on the edge
             allocate(x(0:3, 1:max(Qdof, Gdof), 1:nbDim))
             !--------- ^^^  0= volume, 1,2,3= 1st, 2nd, 3rd face ! ...will be used later

             call Set_RTN_integ_node(V_rule, G_rule, Qdof, Gdof, x)

             allocate( qi(1:Gdof, 1:nbDim  ) )
             call EvalMomentumVolumeD(Fdeg, ib, Gdof, x(Tface, 1:Gdof, 1:nbDim), &
               qi(1:Gdof, 1:nbDim))

             if (elem%face(neigh, Tface) > 0) then  !interior face
                 omega_edge = 0.5
             else                   !boundary face
                 omega_edge = 1.
             endif

             do k=1, ndim
                val = 0
                do l=1, Gdof
                   val = val + G_rule%weights(l) * &
                        dot_product(elem%n(Tface, 1:nbDim), qi(l, 1:nbDim)) * &
                        jump_wi(l, k)

                   !write(*,'(a6,i2,20es12.4)') 'Mvol1:',l,qi(l,1:2), val , &
                   !     G_rule%weights(l) * &
                   !    dot_product(elem%n(Tface, 1:nbDim), qi(l, 1:nbDim)) * &
                   !    jump_wi(l, k), &
                   !    dot_product(elem%n(Tface, 1:nbDim), qi(l, 1:nbDim)), &
                   !    jump_wi(l, k)

                enddo
                val = val * omega_edge
                flux_momentums(k, ib) = flux_momentums(k, ib) + val

                !write(203,'(2i3,a12,20es12.4)') Tface,ib, 'edges', val,flux_momentums(k, :)

             enddo !k

             !write(*,'(a6,20es12.4)') 'M***1:',flux_momentums(1, ib), val
             !do l=1,Gdof
             !   write(*,'(a6,i2,20es12.4)') 'Mvol1:',l,qi(l,1:2),flux_momentums(1,ib), val
             !enddo

             deallocate(jump_wi, wi, Ksk, func3, qi)

             if (Tface /= 3) deallocate(x)
             ! there are saved values for volume integ. rule in x(0, ...)
             !deallocate(func )
          enddo  !Tface

          flux_momentums(:, ib) = -state%space%m_IPG * flux_momentums(:, ib)
       else

          ! value of the test function on the edge
          allocate(x(0:3, 1:max(Qdof, Gdof), 1:nbDim))
          !--------- ^^^  0= volume, 1,2,3= 1st, 2nd, 3rd face ! ...will be used later

          !G_rule => state%space%G_rule(elem%face(fGnum, 1))  ! not used in Set_RTN_integ_node
          !call Set_RTN_integ_node(V_rule, G_rule, Qdof, Gdof, x)

          ! volume integ. nodes for volume integrals
          x(0,1:Qdof, 1:2) = V_rule%lambda(1:Qdof,1:2)

       endif

       !write(203,'(a12,8es12.4)') 'flux A', flux_momentums(1, :)

       ! volume integrals

       allocate(R_s(1:Qdof, 1:nbDim, 1:ndim) )
       call Eval_R_s_Elem(Set_R_s, elem, R_s(1:Qdof, 1:nbDim, 1:ndim) )
       !write(*,'(a6,2i5,20es14.6)') 'Rs=',elem%i,1,R_s(1:Qdof, 1, 1:ndim)
       !write(*,'(a6,2i5,20es14.6)') 'Rs=',elem%i,2,R_s(1:Qdof, 2, 1:ndim)

       allocate(f_s(1:Qdof, 1:nbDim, 1:ndim) )

!stop 'all the following'
!print*, 'HERE was stop all the following'

       call Eval_f_s_Elem(Set_f_s, elem, f_s(1:Qdof, 1:nbDim, 1:ndim) )

       !write(*,'(a6,2i5,20es14.6)') 'fs=',elem%i,1,f_s(1:Qdof, 1, 1:ndim)
       !write(*,'(a6,2i5,20es14.6)') 'fs=',elem%i,2,f_s(1:Qdof, 2, 1:ndim)

       ! used in flux estimator
       CDfluxes(1:Qdof, 1:nbDim, 1:ndim) = f_s(1:Qdof, 1:nbDim, 1:ndim) &
            - R_s(1:Qdof, 1:nbDim, 1:ndim)

       allocate( qi(1:Qdof, 1:nbDim  ) )
       call EvalMomentumVolumeD(Fdeg, ib, Qdof, x(0, 1:Qdof, 1:nbDim), qi(1:Qdof, 1:nbDim))

       allocate(func(1:Qdof) )

       do k=1,ndim
          func(1:Qdof) = 0.
          do l=1,nbDim
             func(1:Qdof) = func(1:Qdof) + qi(1:Qdof, l) * CDfluxes(1:Qdof, l, k)
          enddo

          !write(*,'(a6,2i5,20es14.6)') 'func=1',elem%i,k, &
          !     (f_s(1:Qdof, 1, 1) - R_s(1:Qdof, 1, 1) )
          !
          !write(*,'(a6,2i5,20es14.6)') 'func=2',elem%i,k, &
          !     (f_s(1:Qdof, 2, 1) - R_s(1:Qdof, 2, 1) )
          !
          !write(*,'(a6,2i5,20es14.6)') 'funcS',elem%i,k,func(1:Qdof)

          call IntegrateFunction(elem, func(1:Qdof), val)

          flux_momentums(k, ib) = flux_momentums(k, ib) + val
       enddo

       deallocate(qi, x, R_s, f_s, func)

    enddo  !i

    if(ib /= Fdof) then
       print*,'ERROR: Number of momentums and DOF of the RTN element disagree!!', ib, Fdof
       stop
    endif

    !write(*,*) 'elem%i, elem%HGnode =', elem%i, elem%HGnode


    !if (elem%i == 6) then
    !  print*, 'v------------In ConstructFluxCD sub.----------------v'
    !  write(*,'(a15)') 'MMinvRE(ib, it-vse):'
    !  do ib= 1, Fdof
    !    write(*,'(8es14.6)') MMinvRE(ib, 1:Fdof)
    !  enddo
    !endif
    !pause

    call FluxMM2RTN(Fdof, MMinvRE(1:Fdof, 1:Fdof), flux_momentums(1:ndim, 1:Fdof), &
         fluxi(1:ndim, 1:Fdof) )

    !if (elem%i > 5 .and. elem%i < 9) then
       !write(*,'(a12,20es12.4)') 'fluxi(1,:):', fluxi(1,:)
       !pause
    !endif

    !pause

    deallocate(flux_momentums, CDfluxes )

    if(elem%HGnode) deallocate( HGvertex )

    !print*,'  #  end subroutine ConstructFluxCD'
    !stop

  end subroutine ConstructFluxCD


  !> evaluate resiual estimator
  !> \f$ \eta_{R,T}^n := (c_P h_T m_T \f$
  !> \f$ \|f -\partial_t u_{h\tau} -\nabla\cdot {\bf t}_{h\tau}\|_{T\times I_m} )^2\f$
  subroutine ResidualElemEstimCD(elem, Fdeg, Fdof, fluxi, etaRnT)
    type(element), intent(in) :: elem
    integer, intent(in) :: Fdeg, Fdof
    real, dimension(1:ndim, 1:Fdof), intent(in) :: fluxi
    real, dimension(1:ndim), intent(inout) :: etaRnT
    type(basis_rtn_fe), pointer :: loc_RTN
    type(Gauss_rule), pointer :: G_rule
    real, dimension(:,:), allocatable :: DivPsi ! divergence of local RTN basis functions
    real, dimension(:,:), allocatable :: wt     ! time derivative of the approximate sol
    real, dimension(:,:), allocatable :: RTNflux_loc  ! flux reconstruction pw linear in time
    real, dimension(:,:), allocatable :: temp, temp1, x, Fx, f
    real, dimension(1:ndim) :: val
    real :: t, mT
    integer :: k, l, it, dof, Qdof, Qnum, Gnum, Gdof

    dof = elem%dof
    Qdof = elem%Qdof
    Qnum = elem%Qnum

    loc_RTN => state%loc_RTN(Fdeg)

    ! righ hand side (source terms)
    allocate(x(1:Qdof, 1:nbDim))
    x(1:Qdof, 1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)

    allocate(Fx(1:Qdof, 1:nbDim))
    allocate(f(1:ndim, 1:Qdof) )
    allocate(RTNflux_loc(1:ndim, 1:Fdof) )

    ! Fx contains physical coordinates of integ nodes
    call ComputeF(elem, Qdof, x, Fx)

    ! time derivative of the approximate solution
    allocate( wt(1:Qdof, 1:ndim) )
    call Eval_w_t_Elem(elem, wt)


    ! divergence of the RTN basis functions on elem in integ nodes
    allocate(DivPsi(1:Fdof, 1:Qdof) )
    call Eval_Loc_RTN_Div_Elem(elem, loc_RTN, DivPsi )

    allocate( temp(1:ndim, 1:Qdof), temp1(1:ndim, 1:Qdof) )
    ! time independent
!    do k=1, ndim
!       temp(k, 1:Qdof) = - matmul(fluxi(k, 1:Fdof), DivPsi(1:Fdof, 1:Qdof) )
!
!       temp(k, 1:Qdof) = temp(k, 1:Qdof)  - wt(1:Qdof, k)
!
!       !if(elem%i == 2) &
!       !call PlotElemFunctionQ(740+state%time%iter, elem, 'V', Qnum, Qdof,  temp(k,1:Qdof) )
!
!    enddo


    !SS Gnum = 1   ! 15 is the maximal one !!!!!!!!!
    Gnum = 2   ! 15 is the maximal one !!!!!!!!!
    G_rule => state%space%G_rule(Gnum)
    Gdof = G_rule%Qdof

    etaRnT(1:ndim) = 0.
    ! integration over Gauss integ nodes in time
    do it=1, Gdof
       t = G_rule%lambda(it)
       !SS t = 0.
       state%time%ctime = state%time%ttime - state%time%tau(1) * t

       ! time derivative and the divergence of the flux
       !if(elem%i < 2) then
       !   write(*,'(a8,i5,30es12.4)') '..old',elem%i, fluxi(1, 1:Fdof)
       !   write(*,'(a8,i5,30es12.4)') '..new',elem%i, elem%RTNflux(0, 1, 1:Fdof)
       !   write(*,'(a8,i5,30es12.4)') '..new',elem%i, elem%RTNflux(1, 1, 1:Fdof)
       !   write(*,'(a8,i5,30es12.4)') 'diff ',elem%i, fluxi(1, 1:Fdof)-elem%RTNflux(0, 1, 1:Fdof)
        !  write(*,'(a8,i5,30es12.4)') 'diff2 ',elem%i,elem%RTNflux(1, 1, 1:Fdof) -elem%RTNflux(0, 1, 1:Fdof)
        !  write(*,*)'---------------------------------------------------'
       !endif

       ! linear interpolation
       !! t = 0.  ! PIECEWISE CONSTANT RECONSTRUCTION
       RTNflux_loc(1:ndim, 1:Fdof) = t * elem%RTNflux(1, 1:ndim, 1:Fdof) &
            + ( 1. - t)  *  elem%RTNflux(0, 1:ndim, 1:Fdof)

       do k=1, ndim
          temp(k, 1:Qdof) = - matmul(RTNflux_loc(k, 1:Fdof), DivPsi(1:Fdof, 1:Qdof) )
          !temp(k, 1:Qdof) = - matmul(fluxi(k, 1:Fdof), DivPsi(1:Fdof, 1:Qdof) )
          !SS
          temp(k, 1:Qdof) = temp(k, 1:Qdof)  - wt(1:Qdof, k)
       enddo

       ! source (RHS) term
       do l=1, Qdof
          call RHS_Scalar(Fx(l, 1:nbDim), f(1:ndim, l), state%time%ctime)
       enddo

       temp1(1:ndim, 1:Qdof) = temp(1:ndim, 1:Qdof) + f(1:ndim, 1:Qdof)

       !call PlotElemFunctionQ(840+state%time%iter, elem, 'V', Qnum, Qdof,  temp(1,1:Qdof) )


       call IntegrateSquareVectorFunction(elem, temp1(1:ndim, 1:Qdof), val(1:ndim) )

       etaRnT(1:ndim) = etaRnT(1:ndim) + G_rule%weights(it) * val(1:ndim)

       !if(elem%i == grid%nelem) write(100,'(a8,i5, 6es12.4)') 'val :',it, val(:)
    enddo ! it =1,Gdof ! time integration


    etaRnT(1:ndim) = etaRnT(1:ndim) * state%time%tau(1)    ! length of the time interval

    ! old version May 2011
    !mT = state%time%tau(1)/(state%time%FinTime**0.5 * elem%diam)
    !if(elem%CK >0) mT = min(mT, 1./elem%CK**0.5)
    !if(elem%Cb >0) mT = min(mT, 1./elem%Cb**0.5)
    !mT = (mT* elem%diam/ Pi)**2
    !etaRnT(1:ndim) = etaRnT(1:ndim) * mT


    ! new version June 2011
    !etaRnT(1:ndim) = etaRnT(1:ndim) / (Pi * Pi * elem%CTn)  ! etaRnt stores \eta_R^2 !!!
    ! new version October 2011 Poincare -> Friedrichs
    etaRnT(1:ndim) = etaRnT(1:ndim) / (2. * elem%CTn)  ! etaRnt stores \eta_R^2 !!!

    !write(*,'(a6,i5,12es12.4)') '#Rn#',elem%i, &
    !     1./(Pi * elem%CTn**0.5), mT,1./(Pi * elem%CTn**0.5)-  mT, &
    !     state%time%tau(1)/(state%time%FinTime**0.5 * elem%diam), 1./elem%CK**0.5,1./elem%Cb**0.5

    !if(elem%i == grid%nelem) write(100,'(a8,i5, 6es12.4)') 'etaRn:',it, etaRnT(:), &
    !     state%time%tau(1), mT, (state%time%tau(1)/Pi)**2,elem%CK, elem%Cb


    deallocate(DivPsi, temp, temp1, x, Fx, f, RTNflux_loc)

    !write(*,'(a6,i5,12es14.6)') 'etaRn:',elem%i, etaRnT(1:ndim), mT

  end subroutine ResidualElemEstimCD


  !> evaluate flux convective diffusive estimator
  !> \f$ \eta_{DCF,T}^n := ( \f$
  !> \f$ \|{\bf K}(u_h)\nabla u_h - \vec{f}(u_h) +  {\bf t}_{h\tau}\|_{T\times I_m} )^2\f$
  subroutine FluxElemEstimCD(Set_f_s, Set_R_s,  elem, Fdeg, Fdof, fluxi, &
       !!CDfluxes,
       etaDFn)
    interface
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x, ie )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
         integer, intent(in) :: ie
       end subroutine Set_f_s
      subroutine Set_R_s(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
      end subroutine Set_R_s
    end interface

    type(element), intent(in) :: elem
    integer, intent(in) :: Fdeg, Fdof
    real, dimension(1:ndim, 1:Fdof), intent(in) :: fluxi
    !!real, dimension(1:elem%Qdof, 1:nbDim, 1:ndim), intent(in) :: CDfluxes
    real, dimension(:, :, :), allocatable :: fluxes, fluxesF
    real, dimension(1:3, 1:ndim), intent(inout) :: etaDFn ! 1 = total, 2 = space, 3 = time
    type(basis_rtn_fe), pointer :: loc_RTN
    type(Gauss_rule), pointer :: G_rule
    real, dimension(:,:,:), allocatable :: psi  ! local RTN basis functions
    real, dimension(:,:), allocatable :: wt     ! time derivative of the approximate sol
    real, dimension(:,:), allocatable :: RTNflux_loc  ! flux reconstruction pw linear in time
    real, dimension(:,:), allocatable :: x, Fx, f
    real, dimension(:,:,:), allocatable :: f_s, R_s, temp
    real, dimension(1:ndim) :: val
    real :: t, mT
    integer :: it, k, l, dof, Qdof, Qnum, Gnum, Gdof

    dof = elem%dof
    Qdof = elem%Qdof
    Qnum = elem%Qnum

    loc_RTN => state%loc_RTN(Fdeg)

    ! divergence of the RTN basis functions on elem in integ nodes
    allocate(psi(1:Fdof, 1:nbDim, 1:Qdof) )
    call Eval_Loc_RTN_Elem(elem, loc_RTN, psi )

    allocate( temp(1:3, 1:ndim, 1:Qdof) )  ! first index: 1 = total, 2 = space, 3 = time
    allocate(fluxes(1:Qdof, 1:nbDim, 1:ndim), fluxesF(1:Qdof, 1:nbDim, 1:ndim) )
    allocate(R_s(1:Qdof, 1:nbDim, 1:ndim) )
    allocate(f_s(1:Qdof, 1:nbDim, 1:ndim) )

    allocate(RTNflux_loc(1:ndim, 1:Fdof) )

    ! fluxes at time level t_k
    state%time%ctime = state%time%ttime
    call Eval_R_s_Elem(Set_R_s, elem, R_s(1:Qdof, 1:nbDim, 1:ndim) )

    call Eval_f_s_Elem(Set_f_s, elem, f_s(1:Qdof, 1:nbDim, 1:ndim) )

    fluxesF(1:Qdof, 1:nbDim, 1:ndim) = f_s(1:Qdof, 1:nbDim, 1:ndim) &
            - R_s(1:Qdof, 1:nbDim, 1:ndim)

    ! setting of integ  rule
    Gnum = 2   ! 15 is the maximal one !!!!!!!!!
    G_rule => state%space%G_rule(Gnum)
    Gdof = G_rule%Qdof

    etaDFn(1:3,1:ndim) = 0.

    ! integration over Gauss integ nodes in time
    do it=1, Gdof
       t = G_rule%lambda(it)
       state%time%ctime = state%time%ttime - state%time%tau(1) * t

       call Eval_R_s_Elem_at_time(Set_R_s, elem, R_s(1:Qdof, 1:nbDim, 1:ndim), 1.-t )
       !write(*,'(a6,2i5,20es14.6)') 'Rs=',elem%i,1,R_s(1:Qdof, 1, 1:ndim)
       !write(*,'(a6,2i5,20es14.6)') 'Rs=',elem%i,2,R_s(1:Qdof, 2, 1:ndim)

       call Eval_f_s_Elem_at_time(Set_f_s, elem, f_s(1:Qdof, 1:nbDim, 1:ndim), 1.-t )
       !write(*,'(a6,2i5,20es14.6)') 'fs=',elem%i,1,f_s(1:Qdof, 1, 1:ndim)
       !write(*,'(a6,2i5,20es14.6)') 'fs=',elem%i,2,f_s(1:Qdof, 2, 1:ndim)


       ! used in flux estimator
       fluxes(1:Qdof, 1:nbDim, 1:ndim) = f_s(1:Qdof, 1:nbDim, 1:ndim) &
            - R_s(1:Qdof, 1:nbDim, 1:ndim)

       ! linear interpolation
       !! t = 0.  ! PIECEWISE CONSTANT RECONSTRUCTION
       RTNflux_loc(1:ndim, 1:Fdof) = t * elem%RTNflux(1, 1:ndim, 1:Fdof) &
            + ( 1. - t)  *  elem%RTNflux(0, 1:ndim, 1:Fdof)

       do k=1, ndim
          temp(1:3, k, 1:Qdof) = 0.
          do l=1,nbDim
             ! total estimate
             temp(1, k, 1:Qdof) = temp(1, k, 1:Qdof) + &
                  (matmul(RTNflux_loc(k, 1:Fdof), psi(1:Fdof, l, 1:Qdof) ) &
                  - fluxes(1:Qdof, l, k) )**2

             ! space estimate (can be outside of the cycles)
             temp(2, k, 1:Qdof) = temp(1, k, 1:Qdof) + &
                  (matmul(RTNflux_loc(k, 1:Fdof), psi(1:Fdof, l, 1:Qdof) ) &
                  - fluxesF(1:Qdof, l, k) )**2

             ! time estimate
             temp(3, k, 1:Qdof) = temp(1, k, 1:Qdof) + &
                  (fluxesF(1:Qdof, l, k) - fluxes(1:Qdof, l, k) )**2

             !write(*,'(a8,2i5,20es12.4)') 'theta',elem%i,l, &
             !     matmul(fluxi(k, 1:Fdof), psi(1:Fdof, l, 1:Qdof) )
             !
             !write(*,'(a8,2i5,20es12.4)') 'flux   ',elem%i,l, fluxes(1:Qdof, l, k)
             !write(*,'(a8,2i5,20es12.4)') 'fluxCD ',elem%i,l, CDfluxes(1:Qdof, l, k)


          enddo
          !call PlotElemFunctionQ(840+state%time%iter, elem, 'V', Qnum, Qdof,  temp(k,1:Qdof) )

       enddo

       do l=1,3
          call IntegrateVectorFunction(elem, temp(l, 1:ndim, 1:Qdof), val(1:ndim))
          etaDFn(l, 1:ndim) = etaDFn(l, 1:ndim) + G_rule%weights(it) * val(1:ndim)
       enddo

    enddo

    etaDFn(1:3, 1:ndim) = etaDFn(1:3, 1:ndim) * state%time%tau(1)    ! length of the time interval

    ! etaDFn stores \eta_DF^2 !!
    etaDFn(1:3, 1:ndim) = etaDFn(1:3, 1:ndim) / (elem%diam**2 * elem%CTn)

    !write(*,'(a6,i5,12es12.4)') '#DFn',elem%i, &
    !     1./(elem%diam * elem%CTn**0.5), mT,1./(elem%diam * elem%CTn**0.5)-  mT

    deallocate(psi, temp, R_s, f_s, fluxes, fluxesF, RTNflux_loc )

    !write(*,'(a6,i5,12es14.6)') 'etaDF:',elem%i, etaDFn(1:ndim)

  end subroutine FluxElemEstimCD

  !> evaluate the momentums of the basis RTN function of degree Fdeg on elem
  !> and the resulting matrix is stored in MMRE
  subroutine ComputeLocRTNMomentumsElem(elem, Fdeg, Fdof, MMRE)
    type(element), intent(in) :: elem
    integer, intent(in) :: Fdeg, Fdof
    real, dimension(1:Fdof, 1:Fdof), intent(inout) :: MMRE
    type(basis_rtn_fe), pointer :: loc_RTN
    real, dimension(:,:,:), allocatable :: psi
    real, dimension(:,:), allocatable :: x, qi
    integer :: Qdof, Gdof,  Mdof
    integer :: ie, j, ib, it, i
    integer :: HGidx, Tface, Tface1, Tface2
    integer, dimension(:), allocatable :: HGvertex
    real :: edge_momentum

    Mdof = maxval(elem%face(fGdof,:) )
    Mdof = max(Mdof, elem%Qdof)

    allocate( psi(1:Fdof, 1:nbDim, 1: Mdof)  )
    allocate(  qi(1:Mdof, 1:nbDim  ) )

    loc_RTN => state%loc_RTN(Fdeg)

    ib = 0   ! index of momentum

    MMRE(:, :) = 0


    if(elem%HGnode) then
       if(elem%type /= 3) print*,' TROUBLES in dual_estim.f90 with HG'

       allocate( HGvertex(4) )
       HGvertex(1:3)   = elem%HGvertex(1:3 )
       HGvertex(4) = elem%flen + 1
    endif

    ! face momentums
    do ie=1,3 ! loop over triagle edges (including HG nodes )

       if(.not. elem%HGnode) then
          Tface1 = ie   ! only one element edge on the triangle edge
          Tface2 = ie
          HGidx = 1
       else
          Tface1 = HGvertex(ie)
          Tface2 = HGvertex(ie+1) - 1
       endif

       !Gdof = elem%face(fGdof,ie)
       !call Eval_Loc_RTN_Edge(elem, loc_RTN, ie, psi(1:Fdof, 1:nbDim, 1:Gdof) )

       ib = (ie - 1) * (Fdeg + 1) ! index of the momentums for HG nodes

       do j=1, Fdeg+1  ! Fdeg+1 degrees of freedom over face
          ib = ib + 1

          do Tface = Tface1, Tface2
             if( elem%HGnode) HGidx = elem%HGface(2,Tface)

             Gdof = elem%face(fGdof,Tface)

             edge_momentum = 0

             call Eval_Loc_RTN_Edge_HG(elem, loc_RTN, Tface, HGidx, &
                  psi(1:Fdof, 1:nbDim, 1:Gdof), .false. )

             call EvalMomentumFaceHG(ie, HGidx, Fdeg, ib, Gdof, &
                  state%space%G_rule(elem%face(fGnum,Tface)), qi(1:Gdof, 1))

             !call EvalMomentumFace(ie, Fdeg, ib, Gdof, &
             ! state%space%G_rule(elem%face(fGnum,ie) ),  qi(1:Gdof, 1))

             do it = 1, Fdof
                call IntegrateFunctionNormalEdge(elem, Tface, &
                     psi(it, 1:nbDim, 1:Gdof),  qi(1:Gdof,1), edge_momentum )

                !call IntegrateFunctionNormalEdge(elem,ie,psi(it,1:nbDim,1:Gdof),&
                !     qi(1:Gdof,1), MMRE(ib, it) )

                MMRE(ib, it) = MMRE(ib, it) + edge_momentum
             enddo
          enddo !Tface
       enddo !j
    enddo !ie

    ! VOLUME MOMENTUMS
    Qdof = elem%Qdof

    ! volume integ nodes
    allocate(x(1:Qdof, 1:nbDim))
    x(1:Qdof, 1:2) = state%space%V_rule(elem%Qnum)%lambda(1:Qdof,1:2)

    call Eval_Loc_RTN_Elem(elem, loc_RTN, psi(1:Fdof, 1:nbDim, 1:Qdof) )

    do j=1, Fdeg*(Fdeg+1)  ! Fdeg*(Fdeg+1) degrees of freedom over element
       ib = ib + 1

       call EvalMomentumVolumeD(Fdeg, ib, Qdof, x(1:Qdof, 1:nbDim), &
            qi(1:Qdof, 1:nbDim))

       do it = 1, Fdof
          call IntegrateFunctionsVec(elem,  psi(it, 1:nbDim, 1:Qdof), &
               qi(1:Qdof, 1:nbDim),  MMRE(ib, it))

       enddo !it

    enddo !j

    !!! MMRE is not inverse !!!!
    !call MblockInverse(Fdof, MMRE )

    deallocate(x, qi, psi)

    !print*,'  end subroutine ComputeRTNMomentumsRealElem'

  end subroutine ComputeLocRTNMomentumsElem


  !> from the computed RTN momentums (stored in momentums(:,:)) compute
  !> flux (from RTN) as  basis coeficients of LocRTN basis
  !> MMRE stores RTN momentums of LocRTN basis functions
  subroutine FluxMM2RTN(Fdof, MMRE, momentums, flux )
    integer, intent(in) :: Fdof
    real, dimension (1:Fdof, 1:Fdof), intent(in) :: MMRE
    real, dimension(1:ndim, 1:Fdof), intent(in) :: momentums
    real, dimension(1:ndim, 1:Fdof), intent(out) :: flux
    real, dimension (:,:), allocatable :: x
    real, dimension (:, :), allocatable :: MMinvRE
    integer, dimension(:), allocatable :: ipiv
    external :: DGESV
    integer :: k, info


    allocate(MMinvRE(1:Fdof, 1:Fdof) )
    MMinvRE(1:Fdof, 1:Fdof) = MMRE(1:Fdof, 1:Fdof)

    ! NEW VARIANT using LAPACK
    allocate(ipiv(1:Fdof) )
    do k=1, Fdof
       ipiv(k) = k
    enddo

    allocate(x(1:Fdof, 1:ndim) )

    do k=1,ndim
       x(1:Fdof, k) = momentums(k, 1:Fdof)
    enddo

   !print*, 'v------------In FluxMM2RTN sub. before DGESV----------------v'
    !  do k = 1, Fdof
     !    write(*,'(8es12.4)') MMinvRE(k, :)
     ! enddo
     ! pause

    call DGESV(Fdof, ndim, MMinvRE, Fdof, ipiv, x, Fdof, info)
    if(INFO /= 0) print*,'Warning in dual_estim.f90 LAPACK: DGESV,  INFO = ',info

    !print*, 'v------------In FluxMM2RTN sub. after DGESV----------------v'
    !if(INFO /= 0) then
      !do k = 1, Fdof
      !   write(*,'(8es12.4)') MMinvRE(k, :)
      !enddo
      !pause
    !endif

    do k=1,ndim
       flux(k, 1:Fdof) = x(1:Fdof, k)
    enddo



    deallocate(ipiv, x)

    ! OLD VARIANT
    !MMinvRE(1:Fdof, 1:Fdof) = MMRE(1:Fdof, 1:Fdof)
    !call MblockInverse(Fdof, MMinvRE )
    !do k=1, ndim
    !   flux(k, 1:Fdof) = matmul(MMinvRE(1:Fdof, 1:Fdof), momentums(k, 1:Fdof))
    !!!   write(*,'(a6,20es12.4)') 'Old', flux(k,:)
    !enddo

    deallocate(MMinvRE)

  end subroutine FluxMM2RTN


end module dual_estim
