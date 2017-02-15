!> aposteriori error estimation subroutines for Vohralik approach
module apost_estimation
  use main_data  ! contains type(mesh) :: grid for computation
  use geometry
  use problem_oper
  use euler_problem
  use rav_tho_nedelec
  use inviscid_fluxes
  use f_mapping
  use blocks_integ
  use basis
  use eval_sol
  use errorDual


  implicit none

  public:: ComputeApostEstim
  public:: ConstructOswald
  public:: ConstructPotential
  public:: ComputeRTNMomentumsRealElem
  public:: ConstructFlux
  public:: ComputeElemApostEstim
  public:: DiffusiveFluxElemEstimator
  public:: DiffusiveFluxEstimator
  public:: DiffusiveFluxElemEstimator2
  public:: Eval_Dpsi
  public:: NonconformityElemEstimator1
  public:: NonconformityElemEstimator2
  public:: PiolasTransformation
  public:: IC_Estimator
  public:: DataOscillation
  public:: ResidualElemEstimator
  public:: ComputefluxMomentumsRealElem !for testing
contains
  !> estimates the computational error using the approach based on the flux reconstruction
  subroutine ComputeApostEstim( )
    class(element), pointer :: elem
    real, dimension(:,:,:,:), allocatable :: Oswald  ! Oswald inter in Lagr. nodes
    integer :: dof, Rdeg, Rdof, Fdeg, Fdof
    integer :: i
    real, dimension(:,:), allocatable :: estim

    real, dimension(:,:), pointer :: phi, phiT   !for TEST
    real, dimension(:,:), allocatable :: wi, q !for TEST
    real :: val
    integer :: k, Qdof, ist, l !for TEST

    !open(50, file='solution', status='UNKNOWN')
    !do i = 1, grid%nelem
    !   elem => grid%elem(i)
    !   call PlotElemSolution(50, elem)
    !enddo
    !close(50)

    allocate(estim(1:max_eta, 1:ndim) )
    estim(:,:) = 0.

    ! degree of polynomial reconstruction
    !Rdeg = min(maxval(grid%elem(1:grid%nelem)%deg),  MaxRTNImplemented )
    !Rdeg = min(maxval(grid%elem(1:grid%nelem)%deg),  MaxRTNImplemented+1 )
    !Rdeg = 3  ! for bubble functions
    !Fdeg = 3  !min(Rdeg,  MaxRTNImplemented )

    Rdeg = max(maxval(grid%elem(:)%deg), 3)

    !Fdeg = Rdeg - 1
    !Fdeg = Rdeg
    Fdeg = state%space%deg

    !print*,'@@@',Fdeg, Rdeg, maxval(grid%elem(:)%deg )

    Fdof = state%RTN(Fdeg)%dof

    Rdof = (Rdeg+1)*(Rdeg+2)/2


    !write(*,'(a3,5(a8,i5))') '##','Fdeg=',Fdeg,', Fdof=',Fdof, &
    !     ', Rdeg=',Rdeg, ', Rdof=',Rdof, ', Qdof=',grid%elem(1)%Qdof

    ! construction of the Oswald, contains values in Lagr. nodes
    allocate(Oswald(0:1,1:grid%nelem,1:ndim, 1:Rdof) )

    ! Oswald at k-th time level
    call ConstructOswald(0, Rdeg, Rdof, Oswald(0, 1:grid%nelem, 1:ndim, 1:Rdof))

    ! Oswald at (k-1)-th time level
    call ConstructOswald(1, Rdeg, Rdof, Oswald(1, 1:grid%nelem, 1:ndim, 1:Rdof))

    if(state%time%iter == 1) state%err(IC_L2) = 0.

    do i = 1, grid%nelem
       elem => grid%elem(i)
       call ComputeElemApostEstim(elem, Oswald(0:1, i, 1:ndim, 1:Rdof), &
            estim(1:max_eta, 1:ndim), Rdeg, Rdof, Fdeg, Fdof )

    enddo  ! end of i=1,grid%nelem

    if(state%time%iter == 1) state%err(IC_L2) = state%err(IC_L2)**0.5

    ! NEW
    estim(DFn, 1:ndim)  = estim(DFn, 1:ndim) * state%time%tau(1) ! (eta_{DF}^n)^2
    estim(Rn, 1:ndim)   = estim(Rn, 1:ndim)  * 2 * state%time%tau(1) ! (eta_{R}^n)^2 = 2*tau_n * sum_T of (eta_{R,T}^n)^2
    estim(NC1n, 1:ndim) = estim(NC1n, 1:ndim)* state%time%tau(1) ! (eta_{NC,1}^n)^2
    estim(NC2n, 1:ndim) = estim(NC2n, 1:ndim)* state%time%tau(1) ! (eta_{NC,2}^n)^2
    !estim(Osc, 1:ndim)  = estim(Osc, 1:ndim) * 0.5 * state%time%tau(1) ! 0.5*tau_n *sum_T ||\nabla(f(t_{n-1})-f^{\tilde}n)||_T^2
    estim(DFRn, 1:ndim) = estim(DFRn, 1:ndim) * state%time%tau(1) /2

    !val = -1000.
    call ComputeH1ErrorTimeDerDual(2, val )
    !call ComputeH1ErrorTimeDerDualMore(2, val )
    estim(Osc, 1:ndim) = val


    ! \sum_{j=1}^{n-1} (eta_{X}^j)^2 + (eta_{X}^n)^2,   X \in {DF; R; NC,1; NC,2; Osc}
    state%estim(1:Osc, 1:ndim) = state%estim(1:Osc, 1:ndim) + estim(1:Osc, 1:ndim)


    !state%estim(total,1:ndim) = (state%estim(NC1n,1:ndim))**0.5 &
    !     + (state%estim(NC2n,1:ndim))**0.5 &
    !     + 3 * (state%estim(DFn,1:ndim) + state%estim(Rn,1:ndim))**0.5 &
    !     + state%err(IC_L2)*2**0.5 + 3* (state%estim(Osc,1:ndim))**0.5

    state%estim(total,1:ndim) = (state%estim(NC1n,1:ndim))**0.5 &
         + (state%estim(NC2n,1:ndim))**0.5 &
         + 3 * state%estim(DFRn,1:ndim)**0.5 &
         + state%err(IC_L2)* 2**0.5 + 3* (state%estim(Osc,1:ndim))**0.5

    !write(*, '(a22,6es12.4)') 'estim(1:Osc, 1:ndim):', estim(1:Osc, 1)


    !state%errSTnorm(L2H1) =  sqrt(state%errSTnorm(L2H1))
    !write(*, '(a16,2es12.4)') 'state%errSTnorm(L2H1):', state%errSTnorm(L2H1)**0.5
    !write(*, '(a16,es12.4)') 'state%eta(1):', state%eta(1)

    deallocate(estim)
    deallocate(Oswald)

    !print*,'  ???  end subroutine ComputeApostEstim'
    !stop

  end subroutine ComputeApostEstim

  !> recontruction of the flux as an RTN element using the RTN momentums
  !>
  subroutine ConstructFlux(elem, Fdeg, Fdof, flux, MMinvRE, fluxi)
    type(element), intent(inout) :: elem
    integer, intent(in) :: Fdeg, Fdof               ! degree and DOF of the RTN element
    real, dimension(1:ndim, 1:elem%Qdof, 1:nbDim), intent(inout) :: flux !flux in integ nds
    real, dimension(1:Fdof, 1:Fdof), intent(in) :: MMinvRE ! (Momentum matrix)^-1 on elem
    real, dimension(1:ndim, 1:Fdof), intent(inout) ::  fluxi
    !real, dimension(1:ndim, 1:Fdof) :: flux_momentums, fluxi
    real, dimension(:, :), allocatable :: flux_momentums
    real, dimension(:,:,:), allocatable :: Dwi, Dwi1, x, phiRE
    real, dimension(:,:), allocatable :: jump_wi, wi, wi1, qi, clen_aver
    class(element), pointer :: elem1

    type(volume_rule), pointer :: V_rule
    type(Gauss_rule), pointer :: G_rule
    integer ::  Qdof, Gdof
    integer :: i, j, k, j1, l, k1,  ib, ibnd
    real :: val, omega_edge
    real, dimension(:,:), allocatable :: wExact     ! exact solution  in integ nodes
    real, dimension(:,:), allocatable :: xi, Fx

    V_rule => state%space%V_rule(elem%Qnum) !choice of volume quadrature
    Qdof = V_rule%Qdof    !number of volume integration nodes
    ib = 0              !index of momentum

    allocate(flux_momentums(1:ndim, 1:Fdof) )

    ! face momentums of flux
    do j = 1, 3 !loop through edges of triangle  (meshes with HG nodes - not included)
       G_rule => state%space%G_rule(elem%face(fGnum, j)) !choice of edge quadrature
       Gdof = state%space%G_rule(elem%face(fGnum, j))%Qdof !number of edge integration nodes

       allocate( qi(1:Gdof, 1:nbDim  ) )
       allocate(x(0:3, 1:max(Qdof, Gdof), 1:nbDim))   !information for - Qdof, V_rule is useless
       !--------- ^^^  0= volume, 1,2,3= 1st, 2nd, 3rd face

       call Set_RTN_integ_node(V_rule, G_rule, Qdof, Gdof, x)

       allocate(Dwi(1:Gdof,1:ndim,1:nbDim))
       call Eval_Dw_Edge(elem, j, Dwi, .false.)    !Dwi are multiplied by the Gauss weights


       !write(*,'(4es12.4)') Dwi(1:Gdof, ndim, 1)  !TEST
       !write(*,'(4es12.4)') Dwi(1:Gdof, ndim, 2)  !TEST

       allocate(wi(1:Gdof,1:ndim))
       call Eval_w_Edge(elem, j, wi, .false.)
       !write(206,'(a2,i2,a3,4es12.4)') 'j', j, 'wi', wi(1:Gdof, ndim)  !TEST

       allocate(clen_aver(1:Gdof,1:ndim))
       allocate(jump_wi(1:Gdof,1:ndim))

       k1 = 0
       if (elem%face(neigh, j) > 0) then        !interior edge
          k = elem%face(neigh, j)               !index of element sharing edge j with elem
          j1 = elem%face(nei_i, j)

          k1 = k

          elem1 => grid%elem(k)

          allocate(Dwi1(1:elem1%face(fGdof,j1),1:ndim,1:nbDim))
          call Eval_Dw_Edge(elem1, j1, Dwi1, .true.)

          allocate(wi1(1:elem1%face(fGdof,j1),1:ndim))
          call Eval_w_Edge(elem1, j1, wi1, .true.)
          !write(*,'(a2,i2,a4,4es12.4)') 'j', j, 'wi1', wi1(1:elem%face(fGdof,j1), ndim)  !TEST

          do l=1, Gdof
             do k=1, ndim
                clen_aver(l, k) = dot_product( -elem%n(j, 1:nbDim), &
                     (Dwi(l, k, 1:nbDim) + Dwi1(l, k, 1:nbDim)) /2. )
                ! C_W*(h_edge)^(-1)...
                !jump_wi(l, k) =  state%space%sigma * (1./elem%dn(j)) * (wi(l, k) - wi1(l, k))

                !jump_wi(l, k) =  state%space%sigma * state%model%Re1 &
                !     * (wi(l, k) - wi1(l, k)) / elem%dn(j)
                ! VD sigma 
                jump_wi(l, k) =   max(elem%d_gamma , elem1%d_gamma) * state%model%Re1 &
                     * (wi(l, k) - wi1(l, k)) / elem%dn(j)
             enddo !k
          enddo  !l

          !write(*,'(a6,2i5, 6es12.4)') '####',elem%i,k1,sign, jump_wi(:,1)

          deallocate(wi1, Dwi1)

       else                       !boundary edge

          allocate(wExact(1:Gdof, 1:ndim) )
          allocate(Fx(1:Gdof, 1:nbDim))

          call ComputeF(elem, Gdof, x(j, 1:Gdof, 1:nbDim), &
               Fx(1:Gdof, 1:nbDim) )

          do l=1, Gdof
             if(state%modelName == 'scalar' .or.state%modelName == '2eqs' &
                  .or.state%modelName == 'porous') then

                call Exact_Scalar(Fx(l,1:nbDim), wExact(l, 1:ndim), state%time%ctime)
                 !write(100+state%time%iter, *) Fx(j,1:nbDim), wExact(j, 1:ndim)
             else
                print*,'Non implemented ^%$^ in aposter.f90'
                !call Exact_Ringleb(Fx(j,1:nbDim), wExact(j, 1:ndim) )
             endif
          enddo

          !write(*,'(i3,a7,4es12.4)') j, 'wExact:', wExact(1:Gdof, 1:ndim)

          !write(204, '(a6,10es12.4)') 'bedge',-elem%n(j, 1:nbDim), Dwi(1, 1, 1:nbDim)
          do l=1, Gdof
             do k=1, ndim
                clen_aver(l, k) = dot_product( -elem%n(j, 1:nbDim), Dwi(l, k, 1:nbDim) )
                ! C_W*(h_edge)^(-1)...
                !jump_wi(l, k) =  state%space%sigma * state%model%Re1 &
                !     * (wi(l, k) - wExact(l, k)) / elem%dn(j)

                ! VD sigma
                jump_wi(l, k) =  elem%d_gamma * state%model%Re1 &
                     * (wi(l, k) - wExact(l, k)) / elem%dn(j)

                !write(*,'(a6,2i3,6es12.4)') 'JUMP:', &
                !     k,l,wi(l, k), wExact(l,k), wi(l, k) - wExact(l, k)
             enddo !k
             !write(*,'(a3,i3)') 'k:', k
             !write(*,'(a3)') 'n('
             !write(*,'(i3)') j
             !write(*,'(a3)') '):'
             !write(*,'(2F6.2)') elem%n(j, 1:nbDim)
             !write(*,'(F6.2)') elem%dn(j)
             !write(*,'(F6.2)') state%space%C_W
             !write(*,'(4es12.4)') wi(1:Gdof, ndim)
             !write(*,'(i3)') ndim
             !write(*,'(i3,i3,es12.4)') l, k, wi(l, k)
             !write(*,'(i3,es12.4)') l, Dwi(l, 1, 1)
             !write(*,'(i3,es12.4)') l, Dwi(l, 1, 2)
             !write(*,'(a5,i5,a9,2es12.4)') 'Dwi(', l, '1, 1:2):', Dwi(l, 1, 1:2)
             !write(*,'(a25,i3,i3)') 'index of Gauss quadr.', j, elem%face(fGnum,j)
             !write(*,'(a12,i3,a2,i3)') 'dof of face', j, ':', Gdof
             !write(*,'(a13,i3,a8,es12.4)') 'n*nabla w on', j, 'th edge:', clen_aver(l, ndim)
             !write(*,'(a24,i3,a8,es12.4)') 'C_W*(h_edge)^(-1)*w on', j, 'th edge:', jump_wi(l, ndim)
          enddo  !l

          deallocate(wExact)
          deallocate(Fx)

       endif

       !!!!!!!!! SMAZ
       !allocate(Fx(1:Gdof, 1:nbDim))
       !call ComputeF(elem, Gdof, x(j, 1:Gdof, 1:nbDim), Fx(1:Gdof, 1:nbDim) )
       !do l=1,Gdof
       !   select case(elem%i)
       !   case(1)
       !      Dwi(l, 1, 1) = 2.+ 3.*Fx(l,1)**2
       !      Dwi(l, 1, 2) = 1.
       !   case(2)
       !      Dwi(l, 1, 1) = 1.+ 3.*Fx(l,1)**2
       !      Dwi(l, 1, 2) = 2.
       !   end select

       !!   write(44,*) Fx(l , 1:nbDim),  jump_wi(l, 1)*elem%dn(j), clen_aver(l, 1),&
       !!        wi(l, 1), wExact(1, l), wi(l, 1) - wExact(1, l)
       !!   !G_rule%weights(l), , qi(l, 1), elem%dn(j),  &
       !enddo
       !!
       !deallocate(Fx)

       !!!!!!!!! SMAZ
       !deallocate(wExact)

       do i=1, Fdeg+1    ! Fdeg+1 degrees of freedom over face
          ib = ib + 1
          flux_momentums(:, ib) = 0

          call EvalMomentumFace(j, Fdeg, ib, Gdof, G_rule, qi(1:Gdof, 1))

          !if (ib == 3) then
              !write(*,'(a2,i2)') 'ib', ib
              !write(*,'(a2,i2a16,2F6.2)') 'j', j, 'elem%n(j, 1:nbDim)', elem%n(j, 1:nbDim)
              !write(*,'(a12,i3,a2,i3)') 'dof of face', j, ':', Gdof
              !write(*,'(a14,2es12.4)') 'qi(1:Gdof, 1)', qi(1:Gdof, 1)  !TEST
              !write(*,'(a23,2es12.4)') 'G_rule%weights(1:Gdof)', G_rule%weights(1:Gdof)
              !write(*,'(a23,2es12.4)') 'jump_wi(1:Gdof, ndim)', jump_wi(1:Gdof, ndim)
              !write(*,'(a23,2es12.4)') 'clen_aver(1:Gdof, ndim)', clen_aver(1:Gdof, ndim)
          !endif!TEST

          do k=1, ndim
             val = 0
             do l=1, Gdof
                !!! NO -- Dwi are multiplied by the weights the Gauss quadrature already.
                !Dwi are multiplied by the size of the edge through elem%n
                !val = val + G_rule%weights(l) *(jump_wi(l, k) * qi(l, 1) * elem%dn(j)  &
                !     + clen_aver(l, k) * qi(l, 1) )

                val = val + G_rule%weights(l) * qi(l, 1) &
                     *(jump_wi(l, k)  * elem%dn(j) + clen_aver(l, k)  )

                ! TEST
                !val = val + G_rule%weights(l) * qi(l, 1) &
                !     * dot_product(Dwi(l,1,1:2) , elem%n(j, 1:2) )


                !write(*,'(a12,3i5,8es12.4)') 'flux A', i, ib, l, val, &
                !     G_rule%weights(l), jump_wi(l, k), qi(l, 1), elem%dn(j),  &
                !     clen_aver(l, k)

             enddo !l
             flux_momentums(k, ib) = val
             !flux_momentums(k, ib) = val *sign  ! NEW TEST HERE
          enddo  !k

          !if(elem%i == 1 .or. elem%i == 2 ) then
          !   write(*,'(a6,4i2, es12.4,a1,4es12.4a1,4es12.4)') &
          !        'Moment',elem%i,j,k1, ib, flux_momentums(1, ib), ':', &
          !        !jump_wi(:, 1)  * elem%dn(j) + clen_aver(:, 1), '|', &
          !        dot_product(Dwi(1,1,1:2) , elem%n(j, 1:2) ), &
          !        dot_product(Dwi(2,1,1:2) , elem%n(j, 1:2) ),&
          !        dot_product(Dwi(3,1,1:2) , elem%n(j, 1:2) ),&
          !        dot_product(Dwi(4,1,1:2) , elem%n(j, 1:2) ),'|',&
          !        qi(:, 1)
          !   if( i == Fdeg + 1) print*
          !   !write(*,'(a6,20es12.4)') '>',G_rule%lambda(1:Gdof)**2
          !   !write(*,'(a6,20es12.4)') '<',(1.-G_rule%lambda(1:Gdof))**2
          !endif

       enddo  !i

       deallocate(jump_wi)
       deallocate(wi)
       deallocate(qi)
       deallocate(clen_aver)
       deallocate(Dwi)
       deallocate(x)

    enddo  !j


    !volumes momentums of flux
    do i=1, Fdeg*(Fdeg+1)       ! deg*(deg+1) degrees of freedom over element
       ib = ib + 1
       flux_momentums(:, ib) = 0

       do j=1, 3 !loop through edges of triangle  (meshes with HG nodes - not included)

          G_rule => state%space%G_rule(elem%face(fGnum, j)) !choice of edge quadrature
          Gdof = state%space%G_rule(elem%face(fGnum, j))%Qdof !number of edge integ nodes
          allocate(x(0:3, 1:max(Qdof, Gdof), 1:nbDim)) !information for - Qdof, V_rule is redundant now
          !--------- ^^^  0= volume, 1,2,3= 1st, 2nd, 3rd face ! ...will be used later

          call Set_RTN_integ_node(V_rule, G_rule, Qdof, Gdof, x)

          allocate( qi(1:Gdof, 1:nbDim  ) )

          allocate(jump_wi(1:Gdof,1:ndim))

          allocate(wi(1:Gdof,1:ndim))
          call Eval_w_Edge(elem, j, wi, .false.)

          if (elem%face(neigh, j) > 0) then  !interior edge
             k = elem%face(neigh, j)         !index of element sharing edge j with elem
             j1 = elem%face(nei_i, j)

             elem1 => grid%elem(k)

             allocate(wi1(1:elem1%face(fGdof,j1),1:ndim))
             call Eval_w_Edge(elem1, j1, wi1, .true.)

             do l=1, Gdof
                do k=1, ndim
                   jump_wi(l, k) =  wi(l, k) - wi1(l, k)
                enddo !k
             enddo   !l

             omega_edge = 0.5    !parameter in definition of momentums of flux

             deallocate(wi1)

          else              !boundary edge

             ibnd = -elem%face(neigh, j)

             allocate(wExact(1:Gdof, 1:ndim) )
             allocate(Fx(1:Gdof, 1:nbDim))

             call ComputeF(elem, Gdof, x(j, 1:Gdof, 1:nbDim), &
                  Fx(1:Gdof, 1:nbDim) )


             do l=1, Gdof
                if(state%modelName == 'scalar' .or. state%modelName == '2eqs' &
                     .or. state%modelName == 'porous') then
                   call Exact_Scalar(Fx(l,1:nbDim), wExact(l,1:ndim), state%time%ctime)

                 !write(100+state%time%iter, *) Fx(j,1:nbDim), wExact(j, 1:ndim)
                !else
                   !call Exact_Ringleb(Fx(j,1:nbDim), wExact(j, 1:ndim) )
                endif
             enddo

             !print*,'### ctime =', state%time%ctime,',  elem%i =',elem%i,',  j= ',j, &
             !elem%face(neigh, j)
             !do l=1,Gdof
             !   write(*,'(a6,20es14.6)' ) '!@@',x(j,l,1:2), Fx(l,1:2), wExact(l,1:ndim)
             !enddo

             !write(*,'(a6,20es14.6)' ) '@@@',Fx(:,1)
             !write(*,'(a6,20es14.6)' ) '@@@',Fx(:,2)
             !write(*,*)
             !write(*,'(a6,20es14.6)' ) '@@@',grid%b_edge(ibnd)%x_div(:,1)
             !write(*,'(a6,20es14.6)' ) '@@@',grid%b_edge(ibnd)%x_div(:,2)
             !write(*,*) '---------------------------------------------------------'


             do l=1, Gdof
                do k=1, ndim
                   jump_wi(l, k) =  (wi(l, k) - wExact(l, k))
                enddo !k

                !write(400+state%time%iter, '(20es14.6)') &
                !     Fx(l,1:nbDim), wi(l, 1) - wExact(l,1), wi(l,1), wExact(l, 1)

             enddo   !l

             omega_edge = 1.    !parameter in definition of momentums of flux

             deallocate(wExact)
             deallocate(Fx)

          endif

          call EvalMomentumVolume(Fdeg, ib, Gdof, j, x(0:3, 1:Gdof, 1:nbDim), &
               qi(1:Gdof, 1:nbDim))

          !if(elem%i == 2 .or. elem%i == 2 ) then
          !   write(*,'(a6,4i5, es12.4,a1,4es12.4a1,4es12.4)') &
          !        'Momen.',elem%i,i, ib, j, flux_momentums(1, ib), ':', &
          !        qi(:, 1),'|',qi(:,2)
          !   !if( i == Fdeg*(Fdeg + 1) .and. j==3 ) print*
          !endif

          do k=1, ndim
             val = 0
             do l=1, Gdof
                val = val + G_rule%weights(l) * &
                    dot_product(elem%n(j, 1:nbDim), qi(l, 1:nbDim)) * &
                    jump_wi(l, k)

                !write(*,'(a6,i2,20es12.4)') 'Mvol1:',l,qi(l,1:2), val , &
                !     G_rule%weights(l) * &
                !    dot_product(elem%n(j, 1:nbDim), qi(l, 1:nbDim)) * &
                !    jump_wi(l, k), &
                !    dot_product(elem%n(j, 1:nbDim), qi(l, 1:nbDim)), &
                !    jump_wi(l, k)

             enddo
             val = val * omega_edge
             flux_momentums(k, ib) = flux_momentums(k, ib) + val

             !write(203,'(2i3,a12,20es12.4)') j,ib, 'edges', val,flux_momentums(k, :)

          enddo !k

          !write(*,'(a6,20es12.4)') 'M***1:',flux_momentums(1, ib), val
          !do l=1,Gdof
          !   write(*,'(a6,i2,20es12.4)') 'Mvol1:',l,qi(l,1:2),flux_momentums(1,ib), val
          !enddo

          deallocate(jump_wi, wi, qi)

          if (j /= 3) deallocate(x) !there are saved values for volume integ. rule in x(0, ...)

       enddo  !j


       !write(203,'(a12,8es12.4)') 'flux A', flux_momentums(1, :)
       ! CHECK HERE
       !flux_momentums(:, ib) = state%space%m_IPG * flux_momentums(:, ib)
       flux_momentums(:, ib) = -state%space%m_IPG * flux_momentums(:, ib)

       !write(203,'(a12,8es12.4)') 'flux A', flux_momentums(1, :)

       allocate( qi(1:Qdof, 1:nbDim  ) )

       call EvalMomentumVolume(Fdeg, ib, Qdof, 0, x(0:3, 1:Qdof, 1:nbDim), qi(1:Qdof, 1:nbDim))

       allocate(Dwi(1:Qdof, 1:ndim, 1:nbDim))
       call Eval_Dw_Elem(elem, Dwi)



       !!!!!!!!! SMAZ
       !allocate(Fx(1:Qdof, 1:nbDim))
       !call ComputeF(elem, Qdof, x(0, 1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim) )
       !do l=1,Qdof
       !   select case(elem%i)
       !   case(1)
       !      Dwi(l, 1, 1) = 2. + 3.*Fx(l,1)**2
       !      Dwi(l, 1, 2) = 1.
       !   case(2)
       !      Dwi(l, 1, 1) = 1. + 3.*Fx(l,1)**2
       !      Dwi(l, 1, 2) = 2.
       !   end select
       !
       !!   write(44,*) Fx(l , 1:nbDim),  jump_wi(l, 1)*elem%dn(j), clen_aver(l, 1),&
       !!        wi(l, 1), wExact(l,1), wi(l, 1) - wExact(l, 1)
       !!   !G_rule%weights(l), , qi(l, 1), elem%dn(j),  &
       !enddo
       !
       !deallocate(Fx)
       !
       !!!!!!!!! SMAZ



       do k=1, ndim
          val = 0.
          do l=1, Qdof

             !if(elem%i==1 .and. l== 1 .and.  j==1) print*,'TEST 2345a'
             !flux_momentums(k, ib) = 0.
             val = val + 0.5 * V_rule%weights(l) * abs(elem%F%JF0) * &
                 dot_product( -Dwi(l, k, 1:nbDim), qi(l, 1:nbDim) )

             !write(*,'(a6,i2,20es12.4)') 'Mvol2:',l,qi(l,1:2), val , &
             !     0.5 * V_rule%weights(l) * abs(elem%F%JF0) * &
             !    dot_product( -Dwi(l, k, 1:nbDim), qi(l, 1:nbDim) ), &
             !    -Dwi(l, k, 1:nbDim), V_rule%weights(l), abs(elem%F%JF0)

             ! TEST
             !val = val + 0.5 * V_rule%weights(l) * abs(elem%F%JF0) * &
             !    dot_product( Dwi(l, k, 1:nbDim), qi(l, 1:nbDim) )
          enddo
          !write(201,'(2i3,a12,2es12.4)') i,ib, 'nablaw*r_h:', val,flux_momentums(k, ib)
          flux_momentums(k, ib) = flux_momentums(k, ib) + val

          !if(elem%i == 1 .or. elem%i == 2 ) then
          !   write(*,'(a6,4i4, es12.4,a1,4es12.4a1,4es12.4)') &
          !        'Moment',elem%i,i, ib, -1, flux_momentums(1, ib), ':', &
          !        qi(1:4, 1),'|',Dwi(1:4, 1, 1)
          !   write(*,'(a6,4i4, es12.4,a1,4es12.4a1,4es12.4)') &
          !        'Moment',elem%i,i, ib, -2, flux_momentums(1, ib), ':', &
          !        qi(1:4, 2),'|',Dwi(1:4, 1, 2)
          !   !if( i == Fdeg*(Fdeg + 1)  )
          !   print*
          !endif

       enddo !k

       !write(*,'(a6,20es12.4)') 'M***2:',flux_momentums(1, ib), val
       !do l=1,Qdof
       !   write(*,'(a6,i2,20es12.4)') 'Mvol2:',l,qi(l,1:2),flux_momentums(1, ib), val
       !enddo


       !!write(*,'(a6,i2, 20es12.4)') 'Moment',ib, flux_momentums(1, ib)

       deallocate(qi, Dwi, x)

    enddo  !i

    flux_momentums(1:ndim, 1:Fdof) = flux_momentums(1:ndim, 1:Fdof) * state%model%Re1
    !write(*,'(a8, 20es12.4)') 'Moments',flux_momentums(1, :)


    if(ib /= Fdof) then
       print*,'ERROR: Number of momentums of flux and DOF of the RTN element disagree!!', ib, Fdof
       stop
    endif

    !do l=1,Fdof
    !   !write(*,'(10es12.4)') state%RTN(Fdeg)%MM(l, 1:Fdof)
    !enddo

    !Is RTN(Fdeg) evaluated in the same volume integ nodes as the basis functions of S_hp ?
    if(.not. state%RTN(Fdeg)%Vnode(V_rule%Qdeg)%def) &
       call Init_RTN_FE_Vnode(state%RTN(Fdeg), state%space%V_rule(elem%Qnum) )

    allocate(phiRE(1:Fdof, 1:nbDim, 1:Qdof))

    allocate(x(0:3, 1:max(Qdof, Gdof), 1:nbDim))    !information for - Gdof, G_rule is redundant
    !--------- ^^^  0= volume, 1,2,3= 1st, 2nd, 3rd face

    !print*,'????',V_rule%Qdeg, G_rule%Qdeg, Qdof, Gdof
    call Set_RTN_integ_node(V_rule, G_rule, Qdof, Gdof, x)

    ! SMAZ E1
    !allocate(Fx(1:Qdof, 1:nbDim))
    !call ComputeF(elem, Qdof, x(0, 1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim) )

    do l=1, Fdof
       call PiolasTransformation(elem, &
            state%RTN(Fdeg)%Vnode(V_rule%Qdeg)%phi(l, 1:nbDim, 1:Qdof), &
            phiRE(l, 1:nbDim, 1:Qdof),  x(0, 1:Qdof, 1:nbDim) )

       ! SMAZ E1
       !do i=1,Qdof
       !   write(1000+100*elem%i+l,*) Fx(i, 1:2), phiRE(l, 1:2, i),x(0,i,1:2), &
       !        state%RTN(Fdeg)%Vnode(V_rule%Qdeg)%phi(l, 1:2, i)
       !enddo

    enddo !l

    !SMAZ E1
    !deallocate(Fx)



    !write(201,'(A4,100es12.4)') 'Mat', flux_momentums(1, 1:Fdof)
    !print*,'############################   MMinvRE'
    !do i=1,Fdof
    !   write(*,'(i2,a2,30es12.4)') i,': ',MMinvRE(i,:)
    !enddo

    do k=1, ndim
       fluxi(k, 1:Fdof) = MATMUL(MMinvRE, flux_momentums(k, 1:Fdof))
       do i=1, Qdof
          !if (i==4) then         !TEST
              !write(*,'(a6,i5)') 'Fdeg:', Fdeg
              !write(*,'(a6,i5)') 'Qnum:', state%RTN(Fdeg)%Qnum
              !write(*,'(a6,i5)') 'Qdof:', state%RTN(Fdeg)%Qdof
              !print*,'Vypis phi(1:Fdof, 1, 0, 4).'
              !print*,'Vypis Vnodephi(1:Fdof, 1, 4).'
              !do l=1,Fdof
                 !write(*,'(es12.4)') state%RTN(Fdeg)%phi(l, 1, 0, i)
                 !write(*,'(es12.4)') state%RTN(Fdeg)%Vnode(V_rule%Qdeg)%phi(l, 1, i)
              !enddo
          !endif           !TEST
          flux(k, i, 1) = dot_product(fluxi(k, 1:Fdof), phiRE(1:Fdof, 1, i))  !state%RTN(Fdeg)%phi(1:Fdof, 1, 0, i)
          flux(k, i, 2) = dot_product(fluxi(k, 1:Fdof), phiRE(1:Fdof, 2, i))  !state%RTN(Fdeg)%phi(1:Fdof, 2, 0, i)

          !write(840+state%time%iter,*) x(0, i, 1:nbDim),  flux(k, i,1:2)

       enddo
    enddo

    !write(*,'(a17,100es12.4)') 'flux_momentums:', flux_momentums(1, 1:Fdof)
    !write(40+state%time%iter,'(a17,100es12.4)') 'flux_momentums:', flux_momentums(1, 1:Fdof)

    deallocate(x)
    deallocate(phiRE)
    deallocate(flux_momentums )

    !TEST----------------------
    !write(*,'(a9,i5)') 'elem%Qnum:', elem%Qnum
    !write(*,'(a9,i5)') 'Qdeg:', V_rule%Qdeg
    !write(*,'(a6,i5)') 'Qdof:', Qdof
    !write(*,'(a6,i5)') 'Fdof:', Fdof
    !write(*,'(a12,8es12.4)') 'fluxmoment:',flux_momentums(1, 1:Fdof)
    !write(*,'(a6,8es12.4)') 'fluxi:',fluxi(1, 1:Fdof)
    !write(*,'(a6,10es12.4)') 'flux1:',flux(1, 1:Qdof, 1)
    !write(*,'(a6,10es12.4)') 'flux2:',flux(1, 1:Qdof, 2)
    !TEST----------------------


    !print*,'  ???  end subroutine ConstructFlux'
    !stop

  end subroutine ConstructFlux

  !> Piola's transformation of RTN basis from reference to real element
  !> evaluated in integ nodes of elem%Qnum rule
  subroutine PiolasTransformation(elem, phihat, phiRE, x)  ! RE = real element
    type(element), intent(in) :: elem
    real, dimension(1:nbDim, 1:elem%Qdof), intent(in) :: phihat
    real, dimension(1:nbDim, 1:elem%Qdof), intent(out) :: phiRE
    real, dimension(1:elem%Qdof, 1:nbDim), intent(in) :: x
    real, dimension(:,:,:), allocatable :: DF
    integer :: Qdof
    integer :: l

    Qdof = elem%Qdof

    allocate(DF(1:Qdof, 1:nbDim, 1:nbDim) )
    call ComputeDF(elem, Qdof, x(1:Qdof, 1:nbDim), DF(1:Qdof, 1:nbDim, 1:nbDim) )

    do l=1, Qdof
       phiRE(1:nbDim, l) = MATMUL(DF(l, 1:nbDim, 1:nbDim), phihat(1:nbDim, l) )/abs(elem%F%JF0)  ! only for linear F!!
    enddo !l

    deallocate(DF)

  end subroutine PiolasTransformation

  !> computation of momentums of RTN basis on the real element
  subroutine ComputeRTNMomentumsRealElem(elem, Fdeg, Fdof, MMinvRE)
    type(element), intent(in) :: elem
    integer, intent(in) :: Fdeg, Fdof
    real, dimension(1:Fdof, 1:Fdof), intent(inout) :: MMinvRE !(Momentum matrix)^-1 on elem
    real, dimension(:,:,:,:), pointer :: phi  ! arguments: (1:Fdof, 1:nbDim, 0:3, 1:Q(G)dof)
    real, dimension(:), pointer :: weights
    real, dimension(:,:), allocatable :: qi
    real, dimension(:,:,:), allocatable :: x, DF
    integer :: Qdof, Gdof, Qdeg
    integer :: i,j, ib, it, l
    !real, dimension(1:nbDim) :: vektor

    Qdof = state%space%V_rule(state%RTN(Fdeg)%Qnum)%Qdof
    Gdof = state%space%G_rule(state%RTN(Fdeg)%Gnum)%Qdof


    allocate( qi(1:max(Qdof, Gdof), 1:nbDim  ) )
    phi => state%RTN(Fdeg)%phi(1:Fdof, 1:nbDim, 0:3, 1:max(Qdof, Gdof))


    !write(*,'(a7,i5)') 'Fdeg:', Fdeg
    !write(*,'(a7,i5)') 'Qdof:', Qdof
    !write(*,'(a7,i5)') 'Gdof:', Gdof
    !write(*,'(a7,i5)') 'Fdof:', Fdof
    !do l=1, Fdof
    !  write(*,'(a7,i5)') 'kolbaz', l
    !  do i=1, Qdof
    !    write(*,'(es12.4)') phi(l, 1, 1, i)       !pole je implicitne cislovano od 1!!!
    !    write(*,'(es12.4)') phi(l, 2, 1, i)
    !    print*, ' '
    !    write(*,'(es12.4)') state%RTN(Fdeg)%phi(l, 1, 0, i)
    !    write(*,'(es12.4)') state%RTN(Fdeg)%phi(l, 2, 0, i)
    !    print*, '--- '
    !  enddo !i
    !enddo !l

  !  weights => state%space%V_rule(state%RTN(Fdeg)%Qnum)%weights(1:Qdof)
   ! do i=1, Qdof
   !     write(*,'(es12.4)') weights(i)
   !     write(*,'(es12.4)') weights(i)
   !     print*, ' '
   !     write(*,'(es12.4)') state%space%V_rule(state%RTN(Fdeg)%Qnum)%weights(i)
   !     write(*,'(es12.4)') state%space%V_rule(state%RTN(Fdeg)%Qnum)%weights(i)
   !     print*, '--- '
   ! enddo !i

    weights => state%space%G_rule(state%RTN(Fdeg)%Gnum)%weights(1:Gdof)

    allocate(x(0:3, 1:max(Qdof, Gdof), 1:nbDim))
    !--------- ^^^  0= volume, 1,2,3= 1st, 2nd, 3rd face

    call Set_RTN_integ_node(state%space%V_rule(state%RTN(Fdeg)%Qnum), &
         state%space%G_rule(state%RTN(Fdeg)%Gnum), Qdof, Gdof, x)

    allocate(DF(1:max(Qdof, Gdof), 1:nbDim, 1:nbDim) )

    ib = 0   ! index of momentum

    ! face momentums
    do i=1,3 ! loop over triagle edges

       call ComputeDF(elem, Gdof, x(i, 1:Gdof, 1:nbDim), DF(1:Gdof, 1:nbDim, 1:nbDim) )
       !write(*,'(a22)') '*****VYPIS DF:*****'
       !write(*,'(2es12.4)') DF(2, 1, 1), DF(2, 1, 2)
       !write(*,'(2es12.4)') DF(2, 2, 1), DF(2, 2, 2)
       !allocate(Fx(1:Gdof, 1:nbDim))
       !call ComputeF(elem, Gdof, x(i, 1:Gdof, 1:nbDim), Fx(1:Gdof, 1:nbDim) )

       do j=1, Fdeg+1  ! Fdeg+1 degrees of freedom over face
          ib = ib + 1

          MMinvRE(ib, :) = 0

          call EvalMomentumFace(i, Fdeg, ib, Gdof, state%space%G_rule(state%RTN(Fdeg)%Gnum), &
               qi(1:Gdof, 1))

          !write(*,'(a8, 3i5,6es12.4)') '%%', elem%i,i,ib, qi(1:Gdof, 1)
          do it = 1, Fdof

            !TEST-------------------------
            !if (ib == 4) then
              !if (it == 3) then
                 !write(*,'(a7,i5)') 'Gdof:', Gdof
                 !write(*,'(i5,a9)') i, 'ta hrana'
                 !write(*,'(a7,2es12.4)') 'vahy:', weights(1:Gdof)
                 !write(*,'(a7,es12.4)') '1./detDF:', 1./abs(elem%F%JF0)
                 !write(*,'(a15,2es12.4)') 'testovaci fce:', qi(1:Gdof, 1)
                 !write(*,'(a7,2es12.4)') 'psi_3^:', state%RTN(Fdeg)%phi(it, 1:nbDim, i, 1)
                 !write(*,'(a7,2es12.4)') 'psi_3^:', state%RTN(Fdeg)%phi(it, 1:nbDim, i, Gdof)
                 !write(*,'(a7,2es12.4)') 'psi_3:', MATMUL(DF(1, 1:nbDim, 1:nbDim), state%RTN(Fdeg)%phi(it, 1:nbDim, i, 1 ))
                 !write(*,'(a7,2es12.4)') 'psi_3:', MATMUL(DF(Gdof, 1:nbDim, 1:nbDim), &
                                          !state%RTN(Fdeg)%phi(it, 1:nbDim, i, Gdof ))
                 !write(*,'(a9,2es12.4)') 'normala:', elem%n(i, 1:nbDim)

              !endif
            !endif
            !TEST-------------------------

            do l = 1, Gdof
               MMinvRE(ib, it) = MMinvRE(ib, it) + weights(l)/abs(elem%F%JF0) * qi(l, 1) * &
                    !dot_product(MATMUL(DF(l, 1:nbDim, 1:nbDim), state%RTN(Fdeg)%phi(it, 1:nbDim, i, l)),&
                    dot_product(MATMUL(DF(l, 1:nbDim, 1:nbDim), phi(it, 1:nbDim, i+1, l)),&
                    elem%n(i, 1:nbDim) )

               !write(*,'(a6,5i5, 12es16.8)') '???/', i,ib,it, l, Fdeg,&
               !     state%RTN(Fdeg)%phi(it, 1, i, l), &
               !     phi(it, 1, i+1, l), state%RTN(Fdeg)%phi(it, 1, i, l) - &
               !     phi(it, 1, i+1, l)
            enddo !l
          enddo !it

       enddo !j

    enddo !i


    !allocate(DF(1:Qdof, 1:nbDim, 1:nbDim) )

    call ComputeDF(elem, Qdof, x(0, 1:Qdof, 1:nbDim), DF(1:Qdof, 1:nbDim, 1:nbDim) )

    weights => state%space%V_rule(state%RTN(Fdeg)%Qnum)%weights(1:Qdof)

    ! volumes momentums
    do j=1, Fdeg*(Fdeg+1)  ! Fdeg*(Fdeg+1) degrees of freedom over element
       ib = ib + 1

       MMinvRE(ib, :) = 0

       call EvalMomentumVolume(Fdeg, ib, Qdof, 0, x(0:3, 1:Qdof, 1:nbDim), &
            qi(1:Qdof, 1:nbDim))
       do it = 1, Fdof

          do l = 1, Qdof
             MMinvRE(ib, it) = MMinvRE(ib, it) + 0.5* weights(l) &
                  * dot_product(qi(l, 1:nbDim) ,   &
                  !MATMUL(DF(l, 1:nbDim, 1:nbDim), state%RTN(Fdeg)%phi(it, 1:nbDim, 0, l)) )
                  MATMUL(DF(l, 1:nbDim, 1:nbDim), phi(it, 1:nbDim, 1, l)) )
            !vektor = MATMUL(DF(l, 1:nbDim, 1:nbDim), state%RTN(Fdeg)%phi(it, 1:nbDim, 0, l))
            !MMinvRE(ib, it) = vektor(1)

          enddo !l
       enddo !it

    enddo !j

    !write(*,'(a36)') 'matrix momentums na real elem:'
    !j = 0
    !do ib = 1, Fdof
    !   j = j + 1
    !   write(*,'(i2, a10, 50es12.4)') j, 'th moment', MMinvRE(ib, 1:Fdof)
    !enddo
    !print*, ' '

    call MblockInverse(Fdof, MMinvRE )

    !write(*,'(a48)') 'inverse of matrix momentums na real elem:'
    j = 0
    do ib = 1, Fdof
       j = j + 1
       !write(*,'(10es12.4)') MMinvRE(ib, 1:Fdof)
    enddo

    deallocate(qi)
    deallocate(x)
    deallocate(DF)

    !print*,'  ???  end subroutine ComputeRTNMomentumsRealElem'
    !stop

  end subroutine ComputeRTNMomentumsRealElem


  !> TEST: computation of momentums of attained flux on the real element
  subroutine ComputefluxMomentumsRealElem(elem, Fdeg, Fdof, fluxi, fluxmomentum)
    type(element), intent(in) :: elem
    integer, intent(in) :: Fdeg, Fdof
    real, dimension(1:ndim, 1:Fdof), intent(in) ::  fluxi
    real, dimension(1:Fdof), intent(out) :: fluxmomentum
    real, dimension(:,:,:,:), pointer :: phi  ! arguments: (1:Fdof, 1:nbDim, 0:3, 1:Q(G)dof)
    real, dimension(:,:,:,:), allocatable :: phiRE
    real, dimension(:), pointer :: weights
    real, dimension(:,:), allocatable :: qi, Fx
    real, dimension(:,:,:), allocatable :: x, DF, fluxRE
    real:: sign
    integer :: Qdof, Gdof, Qdeg
    integer :: i,j, ib, it, l, k

    Qdof = state%space%V_rule(state%RTN(Fdeg)%Qnum)%Qdof
    Gdof = state%space%G_rule(state%RTN(Fdeg)%Gnum)%Qdof


    allocate( qi(1:max(Qdof, Gdof), 1:nbDim  ) )
    phi => state%RTN(Fdeg)%phi(1:Fdof, 1:nbDim, 0:3, 1:max(Qdof, Gdof))
    allocate( phiRE(1:Fdof, 1:nbDim, 1:4, 1:max(Qdof, Gdof)) )
    allocate( fluxRE(1:nbDim, 1:4, 1:max(Qdof, Gdof)) )

    weights => state%space%G_rule(state%RTN(Fdeg)%Gnum)%weights(1:Gdof)

    allocate(x(0:3, 1:max(Qdof, Gdof), 1:nbDim))
    !--------- ^^^  0= volume, 1,2,3= 1st, 2nd, 3rd face

    call Set_RTN_integ_node(state%space%V_rule(state%RTN(Fdeg)%Qnum), &
         state%space%G_rule(state%RTN(Fdeg)%Gnum), Qdof, Gdof, x)

    allocate(DF(1:max(Qdof, Gdof), 1:nbDim, 1:nbDim) )

    ib = 0   ! index of momentum


    ! face momentums
    if(elem%flen > 3) print*,'TRouble in ComputefluxMomentumsRealElem at aposter.f90'
    do i=1,3 ! loop over triagle edges


       call ComputeDF(elem, Gdof, x(i, 1:Gdof, 1:nbDim), DF(1:Gdof, 1:nbDim, 1:nbDim) )

       do k=1, Fdof
         do l=1, Gdof
          phiRE(k, 1:nbDim, i+1, l) = &
               MATMUL(DF(l, 1:nbDim, 1:nbDim), phi(k, 1:nbDim, i+1, l) )/abs(elem%F%JF0)
          ! only for linear F!!
         enddo !l
       enddo !k

       do  l=1, Gdof
         do k=1, nbDim
           fluxRE(k, i+1, l) = dot_product(fluxi(ndim,1:Fdof), phiRE(1:Fdof, k, i+1, l) )
         enddo !k
       enddo !l


       ! SMAZ E2
       !sign = 1.
       !if(elem%i < elem%face(neigh,i) ) sign = -1.
       !allocate(Fx(1:Gdof, 1:2) )
       !call ComputeF(elem, Gdof, x(i, 1:Gdof, 1:nbDim), Fx(1:Gdof, 1:nbDim) )
       !do  l=1, Gdof
       !   write(100+elem%i, *) Fx(l,1:2), fluxRE(1:2,i+1,l)
       !   write(200+elem%i, *) Fx(l,1:2), dot_product(fluxRE(1:2,i+1,l),elem%n(i,1:2) ),&
       !        sign*dot_product(fluxRE(1:2,i+1,l),elem%n(i,1:2) )
       !enddo
       !deallocate(Fx)
       ! SMAZ E2

       do j=1, Fdeg+1  ! Fdeg+1 degrees of freedom over face
          ib = ib + 1

          fluxmomentum(ib) = 0

          call EvalMomentumFace(i, Fdeg, ib, Gdof, state%space%G_rule(state%RTN(Fdeg)%Gnum), &
               qi(1:Gdof, 1))

          do l = 1, Gdof
               fluxmomentum(ib) = fluxmomentum(ib) + weights(l) * qi(l, 1) * &
                    dot_product(fluxRE(1:nbDim, i+1, l), elem%n(i, 1:nbDim))
          enddo !l


       enddo !j


    enddo !i

    call ComputeDF(elem, Qdof, x(0, 1:Qdof, 1:nbDim), DF(1:Qdof, 1:nbDim, 1:nbDim) )

    weights => state%space%V_rule(state%RTN(Fdeg)%Qnum)%weights(1:Qdof)

    do k=1, Fdof
         do l=1, Qdof
          phiRE(k, 1:nbDim, 1, l) = &
               MATMUL(DF(l, 1:nbDim, 1:nbDim), phi(k, 1:nbDim, 1, l) )/abs(elem%F%JF0)
          ! only for linear F!!
         enddo !l
    enddo !k

    do  l=1, Qdof
      do k=1, nbDim
        fluxRE(k, 1, l) = dot_product(fluxi(ndim,1:Fdof), phiRE(1:Fdof, k, 1, l) )
      enddo !k
    enddo !l

    ! SMAZ E2
    !allocate(Fx(1:Qdof, 1:2) )
    !call ComputeF(elem, Qdof, x(i, 1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim) )
    !do  l=1, Gdof
    !   write(100+elem%i, *) Fx(l,1:2), fluxRE(1:2,1,l)
    !enddo
    !deallocate(Fx)
    ! SMAZ E2

    ! volumes momentums
    do j=1, Fdeg*(Fdeg+1)  ! Fdeg*(Fdeg+1) degrees of freedom over element
       ib = ib + 1

       fluxmomentum(ib) = 0

       call EvalMomentumVolume(Fdeg, ib, Qdof, 0, x(0:3, 1:Qdof, 1:nbDim), &
            qi(1:Qdof, 1:nbDim))

       do l = 1, Qdof
             fluxmomentum(ib) = fluxmomentum(ib) + 0.5* weights(l)*abs(elem%F%JF0)  &
                  * dot_product(qi(l, 1:nbDim),  fluxRE(1:nbDim, 1, l))
       enddo !l


    enddo !j

    !write(*,'(a15,100es12.4)') 'fluxmomentum:', fluxmomentum(1:Fdof)
    !write(50+state%time%iter,'(a17,100es12.4)') 'flux_momentums:', fluxmomentum(1:Fdof)

    deallocate(qi)
    deallocate(x)
    deallocate(DF)
    deallocate(phiRE)
    deallocate(fluxRE)

    !print*,'  ???  end subroutine ComputefluxMomentumsRealElem'
    !stop


  end subroutine ComputefluxMomentumsRealElem

  !> recontruction of the potential from the Oswald interpolation satifying
  !> the mean value conservation
  subroutine ConstructPotential(elem, Rdeg, Oswald, Qnum, Rdof, potential, time_lvl )
    type(element), intent(inout) :: elem
    integer, intent(in) :: Rdeg !  degree of reconstruction
    integer, intent(in)  :: Rdof !  DOF of reconstruction
    integer, intent(in)  :: Qnum !  index of V_rule
    real, dimension(1:ndim, 1:Rdof), intent(inout) :: Oswald  ! Oswald inter in Lagr. nodes
    real, dimension(1:ndim, 1:Rdof), intent(inout) :: potential ! potential recons. in basis
    integer, intent(in) ::  time_lvl

    real, dimension(:), pointer :: weights
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:), allocatable :: bubble, Rbubble
    real, dimension(:,:), allocatable :: wi, q
    integer :: Bdeg =3, Bdof = 10      ! degree of bubble function (cubic) !!!
    integer :: dof, Qdof
    integer :: k, l, i, j
    real :: coeff, val, valOS, valW, valB

    ! recomputation of the Oswald interpolation from Lagr. nodes into basis
    call Lagr2BasisDiff(elem, Rdeg,  Oswald(1:ndim, 1:Rdof), Qnum, Rdof, &
         potential(1:ndim, 1:Rdof) )

    !call PlotElemFunctionQ(500+state%time%iter, elem, 'L', Rdeg, Rdof, Oswald(1,1:Rdof) )
    !call PlotElemFunctionB(600+state%time%iter, elem, 'V', Rdof, potential(1,1:Rdof) )

    !!write(*,'(a20,20es12.4)') 'POTENTIAL=',potential(1:ndim, 1:Rdof)

    ! bubble function, cubic
    allocate(bubble(1:ndim,1:Bdof), Rbubble(1:ndim, 1:Bdof))
    ! bubble function in Lagr. node
    Rbubble(1:ndim,1:Bdof) =  0.
    Rbubble(1:ndim,6) = 1.   ! center node

    ! recomputation of the bubble function into basis
    call Lagr2BasisDiff(elem, 3,  Rbubble(1:ndim, 1:Bdof), Qnum, Bdof, &
         bubble(1:ndim, 1:Bdof) )

    !call PlotElemFunctionB(800+state%time%iter, elem, 'V', Bdof, bubble(1,1:Bdof) )

    dof = elem%dof

    Qdof = state%space%V_rule(Qnum)%Qdof
    phi => state%space%V_rule(Qnum)%phi(1:max(Rdof, dof, Bdof), 1:Qdof)
    weights => state%space%V_rule(Qnum)%weights(1:Qdof)

    ! potential = Oswald + coef * bubble   ,  hence compute coeff
    do k=1,ndim
       valB = 0.  ! integral of bubble function
       valW = 0.  ! integral of the solution
       valOS = 0. ! intergal of the Oswald interpolation

       do l=1,max(Rdof, dof, Bdof)    ! val = \int_K \phi_l
          val = dot_product(phi(l, 1:Qdof), weights(1:Qdof)) * elem%F%JF0 * 0.5

          !if(l <= dof)  valW = valW + val* elem%w(0,(k-1)*dof + l) ! ERROR H3S2
          if(l <= dof)  valW = valW + val* elem%w(time_lvl,(k-1)*dof + l)
          if(l <= Rdof) valOS = valOS + val*potential(k, l)
          if(l <= Bdof)   valB = valB + val*bubble(k, l)
       enddo

       coeff = (valW - valOS)/ valB

       !write(*,'(a4,F6.2)') 'valW', valW
       !write(*,'(a6,F6.2)') 'valOS', valOS
       !write(*,'(a24,10F6.2)') 'potential(k, 1:Bdof)', potential(k, 1:Bdof)

       potential(k, 1:Bdof) = potential(k, 1:Bdof) + coeff*bubble(k, 1:Bdof)
    enddo

    deallocate(bubble, Rbubble)

    !call PlotElemFunctionB(900+state%time%iter, elem, 'V', Rdof, potential(1,1:Rdof) )
    !call PlotElemFunctionB(950+state%time%iter, elem, 'V', Bdof, potential(1,1:Bdof) )

  end subroutine ConstructPotential

  !> reconstruction of the Oswald interpolation, ONLY for triangles
  !> at time level \f$ k - time_level \f$
  subroutine ConstructOswald(time_level, Rdeg, Rdof, Oswald)
    integer, intent(in) :: time_level
    integer, intent(in) :: Rdeg  ! degree of Oswald reconstruction
    integer, intent(in) :: Rdof  ! DOF of Oswald reconstruction
    real, dimension(1:grid%nelem,1:ndim, 1:Rdof), intent(inout) :: Oswald
    class(element), pointer :: elem, elem1
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:), allocatable :: wi, q
    integer, dimension(:), allocatable :: weight

    integer :: nelem, dof, deg, Qdof, Qnum
    integer :: i, j, k, ist, l, l1, l2, ie, ii, jj, iie

    nelem = grid%nelem

    ! degree of Oswald reconstruction
    !Qnum = state%space%Qdeg(Rdeg,1)
    Qnum = Rdeg

    Qdof = state%space%L_rule(Qnum)%Qdof
    if(Qdof /= Rdof) then
       print*,'Trouble in aposter.f90', Qnum, Qdof, Rdof
       stop
    endif

    !print*,'Qdof, Rdof = ', Qdof, Rdof

    allocate(wi(1:ndim, 1:state%space%max_dof), q(1:ndim, 1:Qdof) )

    phi => state%space%L_rule(Qnum)%phi(1:state%space%max_dof, 1:Qdof)

    ! setting of discontinuous piecewise polynomial solution into Lagr. integ nodes
    ! temporay in array Oswald
    do i=1,nelem
       elem => grid%elem(i)
       if(elem%type /= 3) then
          print*,'Elem ',i,' in ConstructOswald is not a triangle'
          stop
       endif

       dof = elem%dof
       do k=1,ndim
          ist = (k-1)*dof + 1
          wi(k, 1:dof) = elem%w(time_level,ist:ist+dof-1)

          ! evaluation of w in the Langrangian nodes
          do l=1, Qdof
             q(k, l) = dot_product( wi(k, 1:dof), phi(1:dof,l) )
          enddo
       enddo
       Oswald(i,1:ndim, 1:Rdof) = q(1:ndim, 1:Qdof)
       !!write(*,'(a5,i5,20es12.4)') 'pot',i,q(1:ndim, 1:Qdof)
    enddo


    deallocate(q)

    ! avereging of the solution in the nodes of the mesh
    allocate(q(1:grid%npoin, 1:ndim))
    allocate(weight(1:grid%npoin))

    q(:,:) = 0
    weight(:) = 0

    do i=1,nelem
       elem => grid%elem(i)
       do j=1, elem%type

          l = VertexIndex(Rdeg, j)
          k = elem%face(idx,j)
          if(elem%HGnode) k = elem%face(idx,elem%HGvertex(j) ) ! for hanging node

          q(k, 1:ndim) =  q(k, 1:ndim) + Oswald(i, 1:ndim, l)
          weight(k) = weight(k) + 1
       enddo

    enddo

    do k=1,grid%npoin
       q(k, 1:ndim) = q(k, 1:ndim) /weight(k)

       !!write(*,'(a5,i5,20es12.4)') 'nod',k,q(k,1:ndim)
    enddo


    ! inserting new values in nodes
    do i=1,nelem
       elem => grid%elem(i)
       do j=1, elem%type

          l = VertexIndex(Rdeg, j)
          k = elem%face(idx,j)
          if(elem%HGnode) k = elem%face(idx,elem%HGvertex(j) ) ! for hanging node

          Oswald(i, 1:ndim, l) = q(k, 1:ndim)
       enddo
       !!write(*,'(a5,i5,20es12.4)') 'POT',i,Oswald(i, 1:ndim, 1:Qdof)
    enddo

    ! averaging of the solution on edges (without HG nodes)
    do i=1,nelem
       elem => grid%elem(i)
       do j=1, elem%flen
          ie = j
          if(elem%HGnode ) then
             if(elem%HGface(2, j) > 1 ) goto 100 ! no average
             ie = elem%HGface(1,j)
          endif

          if(elem%face(neigh, j) > 0 ) then
             ii = elem%face(neigh, j)
             jj = elem%face(nei_i, j)

             iie = jj

             elem1 => grid%elem(ii)
             if(elem1%HGnode)then
                if(elem1%HGface(2, jj) > 1) goto 100 ! no average
                iie = elem1%HGface(1, jj)
             endif

             do l= 1, Rdeg - 1
                l1 = EdgeIndex(Rdeg, ie, l)
                l2 = EdgeIndex(Rdeg, iie, Rdeg - l)

                Oswald(i,1:ndim, l1) = (Oswald(i,1:ndim, l1) + Oswald(ii,1:ndim, l2) )/ 2.
                Oswald(ii,1:ndim, l2) = Oswald(i,1:ndim, l1)


             enddo
          endif
100       continue
       enddo
    enddo
    deallocate( wi, q)


    ! exact solution in integ nodes
    !do i= 1, nelem
    !   elem => grid%elem(i)
    !   call PlotElemFunctionQ(300+state%time%iter+10*time_level, elem, 'L', Qnum, Qdof, &
    !        Oswald(i, 1, 1:Qdof) )
    !enddo

    !   do j=1, 3
    !      do l= 1, i - 1
    !         print*,'???', i,j,l, EdgeIndex(i, j, l)
    !      enddo
    !      print*
    !   enddo
    !   print*,'------------------------'
    !enddo
    !stop

  end subroutine ConstructOswald



    !> evaluation of the derivatives of the RTN basis functions on element
  subroutine Eval_Dpsi(elem, Fdeg, Der)
    type(element), intent(in):: elem
    integer, intent(in) :: Fdeg
    real, dimension(1:state%RTN(Fdeg)%dof,1:nbDim,1:elem%Qdof,1:nbDim), intent(inout):: Der
    real, dimension(:,:,:,:), pointer:: Dpsi ! pointers to test functions
    real, dimension(:,:), allocatable :: x, Fx
    real, dimension(:,:,:), allocatable :: DF, B1
    integer :: Fdof, Qnum, Qdof, j, l, i

    Qnum = elem%Qnum
    Qdof = elem%Qdof
    Fdof = state%RTN(Fdeg)%dof

    allocate(x(1:Qdof, 1:nbDim), DF(1:Qdof, 1:nbDim, 1:nbDim) )
    allocate(B1(1:Qdof, 1:nbDim, 1:nbDim) )
    x(1:Qdof, 1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)
    call ComputeDF(elem, Qdof, x(1:Qdof, 1:nbDim), DF(1:Qdof, 1:nbDim, 1:nbDim) )

    !write(*, '(a12,i2)') 'RTN-Vrule-Qdof:', state%space%V_rule(state%RTN(Fdeg)%Qnum)%Qdof
    !write(*, '(a12,i2)') 'elem%Qdof:', elem%Qdof

    ! Are derivatives of RTN(Fdeg) bas. f. evaluated in the same volume integ nodes
    ! as the basis functions of S_hp ?
    if(.not. state%RTN(Fdeg)%Vnode(state%space%V_rule(Qnum)%Qdeg)%defder) &
       call Init_RTN_FE_Vnode_Der(state%RTN(Fdeg), state%space%V_rule(Qnum) )


    !Dpsi => state%RTN(Fdeg)%Dphi(1:Fdof, 1:nbDim, 0:3, 1:Qdof, 1:nbDim)
    Dpsi => state%RTN(Fdeg)%Vnode(state%space%V_rule(Qnum)%Qdeg)%Dphi(1:Fdof, 1:nbDim, 1:Qdof, 1:nbDim)

    ! SMAZ E3  -  new test is Der of RTN basis functions are computed correctly
    do l=1,Qdof
       B1(l,1,1) = DF(l,2,2)
       B1(l,1,2) = -DF(l,1,2)
       B1(l,2,1) = -DF(l,2,1)
       B1(l,2,2) = DF(l,1,1)

       B1(l,1:2,1:2) = B1(l,1:2,1:2) /abs(elem%F%JF0)
    enddo

    do j=1, Fdof
       do l=1, Qdof
          Der(j, 1:2, l, 1:2) = &
                matmul(DF(l,1:2,1:2), MATMUL(Dpsi(j, 1:2, l, 1:2), B1(l,1:2,1:2) ) )
          Der(j, 1:2, l, 1:2) = Der(j, 1:2, l, 1:2) /abs(elem%F%JF0)
       enddo
    enddo

    !do j=1, Fdof
    !   do l=1, Qdof
    !      write(33,'(3i5,20es14.6)') elem%i,j,l,Der(j, 1, l, 1:2),Der(j, 2, l, 1:2)
    !   enddo
    !enddo
    deallocate(B1)

    return

    Der(:,:,:,:) = 0.
    ! SMAZ E3


    if(elem%F%iFlin) then        ! linear element
       do j=1, Fdof
          Der(j, 1:nbDim, 1:Qdof, 1) = &
               elem%F%D1F0(1,1) * Dpsi(j, 1:nbDim, 1:Qdof, 1) &
               +elem%F%D1F0(1,2)* Dpsi(j, 1:nbDim, 1:Qdof, 2)

          do l=1, Qdof
             Der(j,1:nbDim, l, 1) = MATMUL(DF(l,1:nbDim, 1:nbDim), Der(j, 1:nbDim, l, 1))&
                  /abs(elem%F%JF0)
          enddo !l

          Der(j, 1:nbDim, 1:Qdof, 2) = &
               elem%F%D1F0(2,1) * Dpsi(j, 1:nbDim, 1:Qdof, 1) &
               +elem%F%D1F0(2,2)* Dpsi(j, 1:nbDim, 1:Qdof, 2)

          do l=1, Qdof
             Der(j, 1:nbDim, l, 2) = MATMUL(DF(l, 1:nbDim, 1:nbDim), Der(j, 1:nbDim, l, 2))&
                  /abs(elem%F%JF0)
          enddo !l

       enddo
    else                         ! curved element
       do j=1, Fdof
          do i=1, 2
            Der(j, i, 1:Qdof, 1) =  &
                 elem%F%V%D1F(1:Qdof, 1,1) * Dpsi(j, i, 1:Qdof, 1) &
                 +elem%F%V%D1F(1:Qdof, 1,2)* Dpsi(j, i, 1:Qdof, 2)

            Der(j, i, 1:Qdof, 2) = &
                 elem%F%V%D1F(1:Qdof, 2,1) * Dpsi(j, i, 1:Qdof, 1) &
                 +elem%F%V%D1F(1:Qdof, 2,2)* Dpsi(j, i, 1:Qdof, 2)

          enddo !i

          do l=1, Qdof
          Der(j, 1:nbDim, l, 1) = 1./abs(elem%F%V%JF(l)) * MATMUL(DF(l, 1:nbDim, 1:nbDim), Der(j, 1:nbDim, l, 1) )
          Der(j, 1:nbDim, l, 2) = 1./abs(elem%F%V%JF(l)) * MATMUL(DF(l, 1:nbDim, 1:nbDim), Der(j, 1:nbDim, l, 2) )
          enddo !l

       enddo
    endif


    !do j=1, Fdof
    !   do l=1, Qdof
    !      write(34,'(3i5,20es14.6)') elem%i,j,l,Der(j, 1, l, 1:2),Der(j, 2, l, 1:2)
    !   enddo
    !enddo
    ! SMAZ E3

    ! SMAZ E1
    !allocate(Fx(1:Qdof, 1:nbDim))
    !call ComputeF(elem, Qdof, x(1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim) )
    !do l=1, Fdof
    !   do i=1,Qdof
    !      write(2000+100*elem%i+l,*) Fx(i, 1:2), Der(l, 1, i, 1:2),Der(l, 2, i, 1:2), &
    !           x(i,1:2), Dpsi(l, 1, i, 1:2), Dpsi(l, 2, i, 1:2)
    !   enddo
    !
    !enddo !l
    !
    !deallocate(Fx)




    deallocate(x, DF)
    !print*, 'End of subroutine Eval_Dpsi'

  end subroutine Eval_Dpsi

  subroutine ResidualElemEstimator(elem, Rdof, potential, Fdeg, fluxi, etaRnT)
    type(element), intent(inout) :: elem
    integer, intent(in) :: Fdeg, Rdof
    real, dimension(0:1, 1:ndim, 1:Rdof), intent(in) :: potential  ! potential reconstr. in basis
    real, dimension(1:ndim, 1:state%RTN(Fdeg)%dof), intent(in) :: fluxi
    real, dimension(1:ndim), intent(inout) :: etaRnT
    !real, dimension(1:state%RTN(Fdeg)%dof,1:nbDim,1:elem%Qdof,1:nbDim) :: Der
    real, dimension(:,:,:,:), allocatable :: Der
    real, dimension(:), pointer :: weights
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:), allocatable :: x, Fx, divpsi, smaz
    real, dimension(:), allocatable :: intercalc1, intercalc2
    !real, dimension(1:ndim, 1:elem%Qdof) :: f
    real, dimension(:, :), allocatable :: f
    integer :: i, k, l, Qdof, Qnum, Fdof, j

    Qdof = elem%Qdof
    Qnum = elem%Qnum
    Fdof = state%RTN(Fdeg)%dof

    !if(elem%i ==1)  print*,'Qnum in ap1.f90', Qnum, Qdof

    weights => state%space%V_rule(Qnum)%weights(1:Qdof)
    phi => state%space%V_rule(Qnum)%phi(1:Rdof, 1:Qdof)


    allocate(x(1:Qdof, 1:nbDim))
    x(1:Qdof, 1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)

    allocate(Fx(1:Qdof, 1:nbDim))
    allocate(f(1:ndim, 1:Qdof) )

    call ComputeF(elem, Qdof, x, Fx)

    !print*,'ctime in ap1.f90', state%time%ctime
    do l=1, Qdof
       call RHS_Scalar(Fx(l, 1:nbDim), f(1:ndim, l), state%time%ctime)
    enddo

    !write(*, '(a12,12es12.4)') 'f(1:Qdof):', f(ndim, 1:Qdof)

    !stop
    allocate(Der(1:Fdof, 1:nbDim, 1:Qdof, 1:nbDim) )
    Der(1:Fdof, 1:nbDim, 1:Qdof, 1:nbDim) = 0.

    call Eval_Dpsi(elem, Fdeg, Der)

    !TEST
    !write(*, '(a6,i2)') 'Qdof:', Qdof
    !write(*, '(a6,i2)') 'Fdof:', Fdof
    !write(*, '(a7,15es12.4)') 'fluxi:', fluxi(1, 1:Fdof)
    !write(*,*) 'Elem = ',elem%i
    !do j=1, Fdof
    !write(*, '(a25,16es12.4)') 'Dpsi(j, 1, 1:Qdof, 1):', Der(j, 1, 1:Qdof, 1)
    !write(*, '(a25,16es12.4)') 'Dpsi(j, 1, 1:Qdof, 2):', Der(j, 1, 1:Qdof, 2)
    !write(*, '(a2)') ' '
    !write(*, '(a25,16es12.4)') 'Dpsi(j, 2, 1:Qdof, 1):', Der(j, 2, 1:Qdof, 1)
    !write(*, '(a25,16es12.4)') 'Dpsi(j, 2, 1:Qdof, 2):', Der(j, 2, 1:Qdof, 2)
    !print*,'-------------------------------------------------'
    !enddo !j
    !TEST

    !print*, 'Subroutine Eval_Dpsi went well.'

    !write(*,'(a6,i2)') 'Rdof:', Rdof
    !write(*,'(a6,i2)') 'Fdof:', Fdof
    !k = 1
    !write(*,'(a24,10F6.2)') 'potential(0,k,1:Rdof)', potential(0, k, 1:Rdof)
    !write(*,'(a24,10F6.2)') 'potential(1,k,1:Rdof)', potential(1, k, 1:Rdof)
    !write(*,'(a20,12es12.4)') 'fluxi(1,1:Fdof):', fluxi(1, 1:Fdof)
    !write(*,'(a12,8es12.4)') 'etaRnT(k):', etaRnT(k)
    !write(*,'(a6,i2)') 'Fdeg:', Fdeg


    allocate(divpsi(1:Fdof, 1:Qdof) )
    divpsi(1:Fdof, 1:Qdof) = 0
    do j=1, Fdof
       divpsi(j, 1:Qdof) = Der(j, 1, 1:Qdof, 1) + Der(j, 2, 1:Qdof, 2)
       !write(*,'(i2,a12,50es12.4)') j, 'divpsi(j,:)', divpsi(j, 1:Qdof)
    enddo !j

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !allocate(smaz(1:Qdof, 1:2))
    !do l=1, Qdof
    !   smaz(l,1) = dot_product(fluxi(1, 1:Fdof), divpsi(1:Fdof, l) )
    !   smaz(l,2) = dot_product(elem%w(0,1:elem%dof) - elem%w(1,1:elem%dof) , &
    !        state%space%V_rule(Qnum)%phi(1:elem%dof, l) )/ state%time%tau(1)
    !enddo

    !!print*,'elem=',elem%i
    !!if(elem%i ==2 .or. elem%i == 3) &
    !call PlotElemFunctionQ(740+state%time%iter, elem, 'V', Qnum, Qdof,  smaz(1:Qdof,1) )
    !call PlotElemFunctionQ(840+state%time%iter, elem, 'V', Qnum, Qdof,  smaz(1:Qdof,2) )
    !call PlotElemFunctionQ(940+state%time%iter, elem, 'V', Qnum, Qdof,  &
    !     smaz(1:Qdof,1)+smaz(1:Qdof,2) )
    !deallocate(smaz)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    deallocate(x)
    deallocate(Fx)

    etaRnT(1:ndim) = 0

    allocate(intercalc1(1:ndim), intercalc2(1:ndim))

    do k=1, ndim
       do l=1, Qdof
          intercalc1(k) = dot_product(potential(0, k, 1:Rdof) - potential(1, k, 1:Rdof), &
                          phi(1:Rdof, l) ) / state%time%tau(1)
          intercalc2(k) = dot_product(fluxi(k, 1:Fdof), divpsi(1:Fdof, l) )

          etaRnT(k) = etaRnT(k) + (f(k, l) - intercalc1(k) - state%model%Re1 * intercalc2(k) )**2 &
                      * elem%F%JF0 * 0.5 * weights(l)

          !write(60+state%time%iter,*) elem%xc(:),f(k, l) - intercalc1(k) - intercalc2(k), &
          !     f(k, l), intercalc1(k), intercalc2(k)
       enddo !l
       !write(*,'(a12,es12.4)') 'etaRnT(k):', etaRnT(k)
    enddo !k

    do k=1, ndim
       etaRnT(k) = etaRnT(k) * (1./(Pi**2)) * (elem%diam)**2
       !write(*,'(a12,es12.4)') 'etaRnT(k)^2:', etaRnT(k)
    enddo !k

    deallocate(divpsi, intercalc1, intercalc2)
    deallocate(Der)

  end subroutine ResidualElemEstimator




  subroutine DiffusiveFluxElemEstimator(elem, Rdof, potential, flux, etaDFnT)
    type(element), intent(inout) :: elem
    integer, intent(in) :: Rdof
    real, dimension(1:ndim, 1:Rdof), intent(in) :: potential  ! potential  in basis
    real, dimension(1:ndim, 1:elem%Qdof, 1:nbDim), intent(in) :: flux !flux in integ n.
    real, dimension(1:ndim),intent(inout) :: etaDFnT
    real, dimension(:,:), allocatable :: func
    real, dimension(:,:,:), allocatable :: Dpotent
    integer :: i, k, l, Qdof, Qnum

    Qdof = elem%Qdof
    Qnum = elem%Qnum

    allocate(Dpotent(1:ndim, 1:Qdof, 1:nbDim ),  func(1:ndim, 1:Qdof) )

    ! evaluation of the gradient of potential s in integ. nodes
    call Eval_DVec_Elem(elem, Rdof, potential(1:ndim, 1:Rdof), &
         Dpotent(1:ndim, 1:Qdof, 1:nbDim) )

    ! \nabla s := \nabla s + \theta (graident of potential + flux)
    Dpotent(1:ndim, 1:Qdof, 1:nbDim) = Dpotent(1:ndim, 1:Qdof, 1:nbDim) &
         + flux(1:ndim, 1:Qdof, 1:nbDim)

    ! \nabla s \cdot \nabla s
    do k=1,ndim
       do l=1,Qdof
          func(k,l) = dot_product(Dpotent(k,l,1:nbDim), Dpotent(k,l,1:nbDim) )
       enddo
    enddo

    call IntegrateVectorFunction(elem, func(1:ndim, 1:Qdof), etaDFnT(1:ndim))

    deallocate(Dpotent,  func )

  end subroutine DiffusiveFluxElemEstimator



  subroutine DiffusiveFluxEstimator(etaDFnT, etaDFn)
    real, dimension(0:1, 1:ndim) :: etaDFnT
    real, dimension(1:ndim), intent(inout) :: etaDFn
    integer :: k

    do k=1, ndim
       etaDFn(k) = etaDFn(k) + state%time%tau(1) * (etaDFnT(0, k) + etaDFnT(1, k))
    enddo   !k

  end subroutine DiffusiveFluxEstimator


  !etaDFnT2 = ||\nabla s^n - \nabla s^(n-1)||_{T}^2
  subroutine DiffusiveFluxElemEstimator2(elem, Rdof, potential, etaDFnT2)
    type(element), intent(inout) :: elem
    integer, intent(in) :: Rdof
    real, dimension(0:1, 1:ndim, 1:Rdof), intent(in) :: potential  ! potential reconstr. in basis
    real, dimension(1:ndim), intent(inout) :: etaDFnT2
    real, dimension(:,:), allocatable :: func
    real, dimension(:,:,:,:), allocatable :: Dpotent
    integer :: i, k, l, Qdof, Qnum

    Qdof = elem%Qdof
    Qnum = elem%Qnum

    allocate(Dpotent(0:1, 1:ndim, 1:Qdof, 1:nbDim ),  func(1:ndim, 1:Qdof) )

    ! evaluation of the gradient of potential s^n in integ. nodes
    call Eval_DVec_Elem(elem, Rdof, potential(0, 1:ndim, 1:Rdof), &
         Dpotent(0, 1:ndim, 1:Qdof, 1:nbDim) )

    ! evaluation of the gradient of potential s^(n-1) in integ. nodes
    call Eval_DVec_Elem(elem, Rdof, potential(1, 1:ndim, 1:Rdof), &
         Dpotent(1, 1:ndim, 1:Qdof, 1:nbDim) )

    ! \nabla (s^n - s^(n-1)) \cdot \nabla (s^n - s^(n-1))
    do k=1,ndim
       do l=1,Qdof
          func(k,l) = dot_product(Dpotent(0,k,l,1:nbDim)-Dpotent(1,k,l,1:nbDim), Dpotent(0,k,l,1:nbDim)-Dpotent(1,k,l,1:nbDim) )
       enddo
    enddo

    call IntegrateVectorFunction(elem, func(1:ndim, 1:Qdof), etaDFnT2(1:ndim))

    deallocate(Dpotent,  func )

  end subroutine DiffusiveFluxElemEstimator2


  subroutine NonconformityElemEstimator1(elem, Rdof, potential, etaNC1nT, time_lvl)
    type(element), intent(inout) :: elem
    real, dimension(1:ndim, 1:Rdof), intent(in) :: potential  ! potential reconstr. in basis
    real, dimension(1:ndim), intent(inout) :: etaNC1nT
    real, dimension(1:ndim, 1:nbDim) :: intercalc
    !real, dimension(1:Rdof,1:nbDim,1:elem%Qdof) :: Der
    real, dimension(:,:,:), allocatable :: Der
    real, dimension(:), pointer :: weights
    real, dimension(:,:,:), pointer:: Dphi ! pointers to test functions
    real, dimension(:,:), allocatable :: wi
    !real, dimension(:,:), allocatable ::x, Fx    !!  SMAZ E4
    integer :: i, k, l, Qdof, Qnum, Rdof, ist, j, time_lvl

    Qdof = elem%Qdof
    Qnum = elem%Qnum
    allocate( wi(1:ndim, 1:state%space%max_dof) )
    weights => state%space%V_rule(Qnum)%weights(1:Qdof)


    !print*, '@@@@',Rdof, elem%dof

    allocate (Der(1:max(elem%dof, Rdof),1:nbDim,1:elem%Qdof))

    ! derivatives of the reference test functions
    Dphi => state%space%V_rule(Qnum)%Dphi(1:max(Rdof, elem%dof), 1:nbDim, 1:Qdof)
    do j=1, max(elem%dof, Rdof)
       Der(j, 1, 1:Qdof) = &
            elem%F%D1F0(1,1) * Dphi(j, 1, 1:Qdof) &
            +elem%F%D1F0(1,2)* Dphi(j, 2, 1:Qdof)

       Der(j, 2, 1:Qdof) = &
            elem%F%D1F0(2,1) * Dphi(j, 1, 1:Qdof) &
            +elem%F%D1F0(2,2)* Dphi(j, 2, 1:Qdof)
    enddo !j

    etaNC1nT(:) = 0

    ! SMAZ E4
    !allocate(Fx(1:Qdof, 1:nbDim))
    !call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:nbDim) )

    do k=1, ndim
       ist = (k-1)*elem%dof + 1
       wi(k, 1:elem%dof) = elem%w(time_lvl,ist:ist+elem%dof-1)
       do l=1, Qdof
          do i=1, 2
             intercalc(k, i) = dot_product(potential(k, 1:Rdof), Der(1:Rdof, i, l) ) &
                  - dot_product(wi(k, 1:elem%dof), Der(1:elem%dof, i, l) )
          enddo !i

          ! SMAZ E4
          !write(50+time_lvl ,*), Fx(l,1:2), &
          !     intercalc(k, 1:2), &
          !     dot_product(potential(k, 1:Rdof), Der(1:Rdof, 1, l) ), &
          !     dot_product(potential(k, 1:Rdof), Der(1:Rdof, 2, l) ), &
          !     dot_product(wi(k, 1:elem%dof), Der(1:elem%dof, 1, l) ), &
          !     dot_product(wi(k, 1:elem%dof), Der(1:elem%dof, 2, l) )

          etaNC1nT(k) = etaNC1nT(k) &
               + weights(l)* dot_product(intercalc(k, 1:nbDim), intercalc(k, 1:nbDim)) * elem%F%JF0 * 0.5
          !write(*,'(a12,es12.4)') 'etaNC1nT(k):', etaNC1nT(k)
       enddo !l
    enddo   !k

    deallocate(Der)

    ! SMAZ E4
    !deallocate(Fx)


  end subroutine NonconformityElemEstimator1


  subroutine NonconformityElemEstimator2(elem, Rdof, potential, etaNC2nT)
    type(element), intent(inout) :: elem
    real, dimension(0:1, 1:ndim, 1:Rdof), intent(in) :: potential  ! potential reconstr. in basis
    real, dimension(1:ndim), intent(inout) :: etaNC2nT
    real, dimension(:), pointer :: weights
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:,:), allocatable :: wi
    integer :: k, Qdof, Qnum, ist, Rdof, l


    Qdof = elem%Qdof
    Qnum = elem%Qnum
    weights => state%space%V_rule(Qnum)%weights(1:Qdof)
    phi => state%space%V_rule(Qnum)%phi(1:max(Rdof, elem%dof), 1:Qdof)
    allocate( wi(0:1, 1:ndim, 1:state%space%max_dof) )
    etaNC2nT(1:ndim) = 0

    do k=1, ndim
       ist = (k-1)*elem%dof + 1
          wi(0, k, 1:elem%dof) = elem%w(0,ist:ist+elem%dof-1)
          wi(1, k, 1:elem%dof) = elem%w(1,ist:ist+elem%dof-1)
          !write(*, '(a3,F12.8)') 'pi:', Pi
          !write(*, '(a5,i2)') 'Rdof:', Rdof
          !write(*, '(a10,i2)') 'elem%dof:', elem%dof
       do l=1, Qdof
          etaNC2nT(k) = etaNC2nT(k) + weights(l)*(dot_product(potential(0, k, 1:Rdof), phi(1:Rdof, l)) &
           - dot_product(wi(0, k, 1:elem%dof), phi(1:elem%dof, l)) &
           - dot_product(potential(1, k, 1:Rdof), phi(1:Rdof, l)) &
           + dot_product(wi(1, k, 1:elem%dof), phi(1:elem%dof, l)) )**2 * elem%F%JF0 * 0.5
          !write(*, '(a14,es12.4)') 'eta_{NC2,T}^n^2:', etaNC2nT(k)
       enddo !l
       etaNC2nT(k) = etaNC2nT(k)* (1./(state%time%tau(1)**2)) * (1./(Pi**2))*(elem%diam)**2
       !write(*, '(a14,es12.4)') 'eta_{NC2,T}^n^2:', etaNC2nT(k)
    enddo   !k


    deallocate(wi)

  end subroutine NonconformityElemEstimator2


  subroutine IC_Estimator(elem, Rdof, potential, etaIC)
    type(element), intent(inout) :: elem
    real, dimension( 1:ndim, 1:Rdof), intent(in) :: potential  ! potential reconstr. in basis
    real, dimension(1:ndim), intent(inout) :: etaIC
    real, dimension(:,:), pointer :: phi
    real, dimension(:), allocatable :: potent
    real, dimension(:,:), allocatable :: xi, Fx, wExact
    integer :: k, Qdof, Qnum, ist, Rdof, j, l


    Qdof = elem%Qdof
    Qnum = elem%Qnum
    phi => state%space%V_rule(Qnum)%phi(1:max(Rdof, elem%dof), 1:Qdof)
    allocate( potent(1:Qdof) )

    allocate(xi(1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim))
    allocate(wExact(1:ndim,1: Qdof))

    ! exact solution in integ nodes
   xi(1:Qdof, 1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)   !
   call ComputeF(elem, Qdof, xi(1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim) )
   do j=1,Qdof
      call Exact_Scalar(Fx(j,1:nbDim), wExact(1:ndim, j), state%time%ctime)
      !write(100+state%time%iter, *) Fx(j,1:nbDim), wExact(1:ndim, j)
   enddo

   !call PlotElemFunctionB(950+state%time%iter, elem, 'V', Rdof, potential(1,1:Rdof) )

   do k=1, ndim
   ! potential in integ nodes
      do j=1, Qdof
         potent(j) = dot_product(potential(k, 1:Rdof), phi(1:Rdof, j) )
         !write(100+state%time%iter, *) Fx(j,1:nbDim), wExact(1:ndim, j),potent(j)
      enddo
      potent(1:Qdof) = (potent(1:Qdof) - wExact(k, 1:Qdof) )**2

      call IntegrateFunction(elem, potent(1:Qdof), etaIC(k)  )
    enddo   !k

    deallocate(potent, xi, Fx, wExact)

  end subroutine IC_Estimator

  subroutine DataOscillation(elem, etaOscT)
    type(element), intent(inout) :: elem
    real, dimension(1:ndim), intent(inout) :: etaOscT
    real, dimension(:), pointer :: weights
    real, dimension(1:ndim, 1:nbDim) :: intercalc
    real, dimension(:,:), allocatable :: x, Fx
    real, dimension(1:ndim, 1:nbDim, 1:elem%Qdof, 0:1) :: Nablaf  ! 0...t_n, 1...t_{n-1}
    integer :: k, Qdof, Qnum, l

    Qdof = elem%Qdof
    Qnum = elem%Qnum
    weights => state%space%V_rule(Qnum)%weights(1:Qdof)

    allocate(x(1:Qdof, 1:nbDim))
    x(1:Qdof, 1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)
    allocate(Fx(1:Qdof, 1:nbDim))

    call ComputeF(elem, Qdof, x, Fx)

    do l=1, Qdof
       call Der_RHS_Scalar(Fx(l, 1:nbDim), Nablaf(1:ndim, 1:nbDim, l, 0), state%time%ttime)
       call Der_RHS_Scalar(Fx(l, 1:nbDim), Nablaf(1:ndim, 1:nbDim, l, 1), &
            state%time%ttime - state%time%tau(1))
    enddo

    etaOscT(1:ndim) = 0

    do k=1, ndim
       do l=1, Qdof
          intercalc(k, 1:nbDim) = Nablaf(k, 1:nbDim, l, 1) - Nablaf(k, 1:nbDim, l, 0)
          etaOscT(k) = etaOscT(k) + weights(l) * elem%F%JF0 * 0.5 &
               * dot_product(intercalc(k, 1:nbDim), intercalc(k, 1:nbDim))
          !write(*,'(a12,es12.4)') 'etaOscT(k):', etaOscT(k)
       enddo !l
    enddo !k

    deallocate(x)
    deallocate(Fx)

  end subroutine DataOscillation

  subroutine ComputeElemApostEstim(elem, Oswald, estim, Rdeg, Rdof, Fdeg, Fdof)
    type(element), intent(inout) :: elem
    real, dimension(1:max_eta, 1:ndim), intent(inout) :: estim
    real, dimension(0:1, 1:ndim, 1:Rdof) :: Oswald
    real, dimension(:,:,:), allocatable :: potential  ! potential reconstr. in basis
    real, dimension(:,:), allocatable :: MMinvRE ! inverse of Momentum matrix on real element
    real, dimension(:,:,:), allocatable :: flux  ! flux reconstr. in basis
    real, dimension(:,:), allocatable :: fluxi  ! coefficients of flux reconstr.
    real, dimension(:,:), allocatable :: etaDFnT, etaNC1nT
    !real, dimension(1:ndim) :: etaNC2nT, etaRnT, etaOscT, etaIC
    real, dimension(:), allocatable :: etaNC2nT, etaRnT, etaOscT, etaIC, etaDFnT2
    !integer :: Qdof, dof, RTNdeg, Qdeg
    integer ::  Rdeg, Rdof, Qnum, Fdeg, Fdof
    real, dimension(1:Fdof) :: fluxmomentum  !for testing

    Qnum = state%space%Qdeg(2*Rdeg,1)
    !Qnum = 2*Rdeg

    allocate(etaNC2nT(1:ndim), etaRnT(1:ndim), etaOscT(1:ndim), etaIC(1:ndim))


    allocate(potential(0:1, 1:ndim, 1:Rdof) )
    allocate(MMinvRE(1:Fdof, 1:Fdof) )
    allocate(flux(1:ndim, 1:state%space%max_Qdof, 1:nbDim) )
    allocate(fluxi(1:ndim, 1:Fdof) )
    allocate(etaDFnT(0:1, 1:ndim), etaNC1nT(0:1, 1:ndim))
    allocate(etaDFnT2(1:ndim))

    call ConstructPotential(elem, Rdeg, Oswald(0, 1:ndim, 1:Rdof), Qnum, Rdof, &
         potential(0, 1:ndim, 1:Rdof), 0 )

    call ConstructPotential(elem, Rdeg, Oswald(1, 1:ndim, 1:Rdof), Qnum, Rdof, &
         potential(1, 1:ndim, 1:Rdof), 1 )

    !write(*,'(a20,20es12.4)')  'Oswald=',Oswald(1, 1:ndim, 1:Rdof)
    !write(*,'(a20,20es12.4)')  'Potential=',potential(1, 1:ndim, 1:Rdof)
    !write(*,*)

    if(state%time%iter == 1) then
       call IC_Estimator(elem, Rdof, potential(0, 1:ndim, 1:Rdof), etaIC)
       state%err(IC_L2) = state%err(IC_L2) + sum(etaIC(:) )
    endif


    call ComputeRTNMomentumsRealElem(elem, Fdeg, Fdof, MMinvRE)  ! RE = real element

    ! flux at k-th time level
    call ConstructFlux(elem, Fdeg, Fdof, flux(1:ndim, 1:elem%Qdof, 1:nbDim), &
         MMinvRE, fluxi(1:ndim, 1:Fdof))

    !the following suborutine only for testing
    !call ComputefluxMomentumsRealElem(elem, Fdeg, Fdof, fluxi, fluxmomentum)

    !call PlotElemSolution(140+state%time%iter, elem)
    !call PlotElemFunctionQ(540+state%time%iter, elem, 'V', elem%Qnum, elem%Qdof,  &
    !     flux(1, 1:elem%Qdof, 1) )
    !call PlotElemFunctionQ(640+state%time%iter, elem, 'V', elem%Qnum, elem%Qdof,  &
    !     flux(1, 1:elem%Qdof, 2) )


    !if (elem%F%iFlin)  !write(*,'(a4,i2,a6)') 'elem',i, 'true'
    !  true for linear mapping F (triangles or paralelogram)
    !write(*,'(a14,i2,e12.4)')  'area of elem',i, elem%area

    call DiffusiveFluxElemEstimator(elem, Rdof, potential(0, 1:ndim, 1:Rdof), &
         flux(1:ndim, 1:elem%Qdof, 1:nbDim), etaDFnT(0, 1:ndim))

    call DiffusiveFluxElemEstimator(elem, Rdof, potential(1, 1:ndim, 1:Rdof), &
         flux(1:ndim, 1:elem%Qdof, 1:nbDim), etaDFnT(1, 1:ndim))


    call DiffusiveFluxElemEstimator2(elem, Rdof, potential(0:1, 1:ndim, 1:Rdof), &
         etaDFnT2(1:ndim) )

    !write(*,'(a20,12es12.4)') 'fluxi(1,1:Fdof):',fluxi(1, 1:Fdof)

    call ResidualElemEstimator(elem, Rdof, potential(0:1, 1:ndim, 1:Rdof), Fdeg,&
         fluxi(1:ndim, 1:Fdof), etaRnT)

    !write(*,'(a20,12es12.4)') 'fluxi(1,1:Fdof):',fluxi(1, 1:Fdof)


    call NonconformityElemEstimator1(elem, Rdof, potential(0, 1:ndim, 1:Rdof), &
         etaNC1nT(0, 1:ndim), 0)
    call NonconformityElemEstimator1(elem, Rdof, potential(1, 1:ndim, 1:Rdof), &
         etaNC1nT(1, 1:ndim), 1)
    !write(*,*)  potential(1, 1:ndim, 1:Rdof)

    !write(*,'(a7,i2,a3)') '***elem', i, '***'

    call NonconformityElemEstimator2(elem, Rdof, potential(0:1, 1:ndim, 1:Rdof), &
         etaNC2nT(1:ndim))

    ! OLD
    !call DataOscillation(elem, etaOscT(1:ndim))


    estim(DFn, 1:ndim)  = estim(DFn, 1:ndim) + etaDFnT(0, 1:ndim) + etaDFnT(1, 1:ndim) ! sum_T of (||\nabla s^(n-1) + \theta^n||_{T}^2 + ||\nabla s^n + \theta^n||_{T}^2)
    estim(Rn, 1:ndim)   = estim(Rn, 1:ndim)  + etaRnT(1:ndim)   ! sum_T of (eta_{R,T}^n)^2

    ! OLD less accurate
    !estim(NC1n, 1:ndim) = estim(NC1n, 1:ndim)+ etaNC1nT(1, 1:ndim) + etaNC1nT(0, 1:ndim)

    ! NEW we integrate square of the linear function over I_m: (g_a^2 + g_b^2 + g_a g_b)/3
    estim(NC1n, 1:ndim) = estim(NC1n, 1:ndim)+ &
         (etaNC1nT(1, 1:ndim) + etaNC1nT(0, 1:ndim)  &
         + (etaNC1nT(1, 1:ndim)* etaNC1nT(0, 1:ndim))**0.5) /3

    ! sum_T of (||\nabla s^(n-1) - \nabla u_{htau}^n||_{T}^2 + ||\nabla s^n - \nabla u_{htau}^n||_{T}^2)

    estim(NC2n, 1:ndim) = estim(NC2n, 1:ndim)+ etaNC2nT(1:ndim) ! sum_T of (eta_{NC2,T}^n)^2
    ! OLD
    !estim(Osc, 1:ndim)  = estim(Osc, 1:ndim) + etaOscT(1:ndim)  !sum_T of || \nabla (f(t_{n-1}) - f^{\tilde}n) ||_{T}^2

    estim(DFRn, 1:ndim)  = estim(DFRn, 1:ndim) &
         + (etaRnT(1:ndim)**0.5 + etaDFnT(0,1:ndim)**0.5)**2 &
         + (etaRnT(1:ndim)**0.5 + etaDFnT(1,1:ndim)**0.5)**2


    elem%eta(Rn,   1:ndim) = etaRnT(1:ndim)**0.5  ! etaRnT
    elem%eta(DFn , 1:ndim) = etaDFnT(0, 1:ndim)**0.5 ! ||\nabla s^n + \theta^n||_T
    elem%eta(NC1n, 1:ndim) = etaNC1nT(0, 1:ndim)**0.5 ! ||\nabla s^n - \nabla u_{htau}^n||_T
    elem%eta(NC2n, 1:ndim) = etaNC2nT(1:ndim)**0.5 ! etaNC2nT

    !OLD
    !elem%eta(Osc,  1:ndim) = etaOscT(1:ndim)**0.5 ! || \nabla (f(t_{n-1}) - f^{\tilde}n) ||_{T}

    !write(*,'(a6,3es12.4,a2,12es12.4)') '???', &
    !     estim(DFn, 1), estim(Rn, 1),estim(DFRn, 1),'|', &
    !     etaDFnT(0, 1:ndim) , etaDFnT(1, 1:ndim), etaRnT(1:ndim), &
    !     (etaRnT(1:ndim)**0.5 + (etaDFnT(0,1:ndim)**0.5 +etaDFnT(1,1:ndim)**0.5)/2)**2


    deallocate(potential, MMinvRE, flux, fluxi, etaDFnT, etaNC1nT)
    deallocate(etaNC2nT, etaRnT, etaOscT, etaIC )
    deallocate(etaDFnT2)

    elem%rezid = 0.

    ! doesn't  work for quadrilaterals
    !if(.not. state%RTN(RTNdeg)%Vnode(Qdeg)%def) &
    !     call Init_RTN_FE_Vnode(state%RTN(RTNdeg), state%space%V_rule(elem%Qnum) )

    ! zde uz jsou k dispozici   state%RTN(RTNdeg)%Vnode(Qdeg)%phi(:, : , :)
    ! ktere predstavuji vujadreni RTN bazovych funkci ve stejnych uzlech jako
    ! bazove funkce S_{hp}

    !stop

  end subroutine ComputeElemApostEstim


end module apost_estimation



