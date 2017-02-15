!> aposteriori error estimation based on Helmholtz decomposition
module helmholtz_estim
  use paramets
  use main_data  ! contains type(mesh) :: grid for computation
  use problem_oper
  use euler_problem
  use rav_tho_nedelec
  use inviscid_fluxes
  use f_mapping
  use blocks_integ
  use basis
  use eval_sol
  use errorDual
  use eval_jumps

  implicit none

  public:: ComputeHelmholtzEstim
  public:: HelmholtzEstimElem
  public:: ElementJumps

contains
  !> evaluate the error estimates based on the Helmholtz decomposition
  subroutine  ComputeHelmholtzEstim( )
    class(element), pointer :: elem
    real, dimension(:,:), allocatable :: estimL
    integer :: i

    !do i=1,grid%nelem
    !   elem => grid%elem(i)
    !   elem%w(0:1,:) = 0.
    !   elem%w(0:1,1) = 0.
    !enddo

    allocate(estimL(1:max_eta, 1:ndim) )


    estimL(:,:) = 0.

    do i=1,grid%nelem
       elem => grid%elem(i)
       !if(i == 1) elem%w(0:1,1) = 1.
       !if(i == 1) elem%w(0:1,2) = 1.
       !if(i == 1) elem%w(0:1,3) = 1.

       call HelmholtzEstimElem(elem)
       elem%eta(Heta1, 1:ndim) = sum(elem%eta(Hrez, 1:ndim) + elem%eta(Hjump_1, 1:ndim) &
            + elem%eta(HjumpD, 1:ndim) )

       !!elem%eta(Heta2, 1:ndim) =   elem%eta(Hjump, 1:ndim)  !! Heta2 is not defined !!

       estimL(Hrez,   1:ndim) = estimL(Hrez,   1:ndim) + elem%eta(Hrez,   1:ndim)**2
       estimL(Hjump,  1:ndim) = estimL(Hjump,  1:ndim) + elem%eta(Hjump,  1:ndim)**2
       estimL(Hjump_1,1:ndim) = estimL(Hjump_1,1:ndim) + elem%eta(Hjump_1,1:ndim)**2
       estimL(HjumpD, 1:ndim) = estimL(HjumpD, 1:ndim) + elem%eta(HjumpD, 1:ndim)**2
       estimL(Heta1,  1:ndim) = estimL(Heta1,  1:ndim) + elem%eta(Heta1,  1:ndim)**2
    enddo
    estimL(Hrez,   1:ndim) = state%time%tau(1) * estimL(Hrez,   1:ndim)
    estimL(Hjump_1,1:ndim) = state%time%tau(1) * estimL(Hjump_1,1:ndim)
    estimL(HjumpD, 1:ndim) = state%time%tau(1) * estimL(HjumpD, 1:ndim)
    estimL(Heta1,  1:ndim) = state%time%tau(1) * estimL(Heta1,  1:ndim)

    estimL(Heta,  1:ndim) = estimL(Heta1,  1:ndim) + estimL(Hjump,  1:ndim)

100 format(a6,2es12.4,'|',14es12.4)
!    write(22,100) &
!         '$$$',state%space%h, state%time%tau(1), estimL(Hrez:Heta1, 1:ndim)

    ! accumulation of the estimates
    state%estim(Hrez:Heta, 1:ndim) = state%estim(Hrez:Heta, 1:ndim) + estimL(Hrez:Heta, 1:ndim)

    deallocate(estimL)

  end subroutine ComputeHelmholtzEstim


  !> evaluate the error estimates based on the Helmholtz decomposition for element
  subroutine HelmholtzEstimElem(elem)
    type(element), intent(inout) :: elem


    call ElementResiduum(elem, elem%eta(Hrez,1:ndim) )

    call ElementJumps(elem, elem%eta(Hjump,1:ndim), elem%eta(HjumpD,1:ndim), elem%eta(Hjump_1,1:ndim) )

    !write(*,'(a6,i5,12es14.6)') 'estim',elem%i, elem%eta(Hrez:HjumpD,1:ndim)

  end subroutine HelmholtzEstimElem

  !> evaluate the residuum
  !> \f$ \|f(t_n) + \Delta u_h^n - \frac{1}{\tau_n}(u_h^n - u_h^{n-1}) \|_K\f$,
  !> elem = \f$K\in T_h\f$
  subroutine ElementResiduum(elem, rez )
    type(element), intent(inout) :: elem
    real, dimension(1:ndim), intent(out) ::  rez
    real, dimension(:,:,:), allocatable :: Dwi  ! derivative of the solution in integ nodes
    real, dimension(:,:,:), allocatable :: D2wi  ! second derivative of the solution
    real, dimension(:,:,:), allocatable :: Dphi ! derivative of the test functions
    real, dimension(:), allocatable :: vector, vector2
    real, dimension(:,:), allocatable :: Dwt
    real, dimension(:,:), allocatable ::  Fx, f
    integer :: dof, Qdof, i, j, k, k1, k2, l

    dof = elem%dof
    Qdof = elem%Qdof

    allocate(Dwi(1:Qdof, 1:ndim, 1:nbDim) )
    call Eval_Dw_Elem(elem, Dwi)

    allocate(D2wi(1:Qdof, 1:ndim, 1:nbDim) )

    allocate(Dphi(1:dof, 1:nbDim, 1:Qdof ) )
    call  Eval_Dphi(elem, dof,  Dphi)

    allocate(vector(1: dof*ndim), vector2(1: dof*ndim) )

    allocate(Fx(1:Qdof, 1:nbDim), f(1:ndim, 1:Qdof) )
    call ComputeF(elem, Qdof, state%space%V_rule(elem%Qnum)%lambda(1:Qdof, 1:nbDim), &
         Fx(1:Qdof, 1:nbDim) )

    do l=1, Qdof
       call RHS_Scalar(Fx(l, 1:nbDim), f(1:ndim, l), state%time%ctime)
    enddo


    do j=1,nbDim
       vector(:) = 0.
       call EvalVectorB(elem, Dwi(1:Qdof,1:ndim, j), dof, vector)

       do k=1,ndim
          k1 = (k-1)*dof + 1
          k2 = k*dof
          ! vector2 contain basis coefficients of Dw in basi phi ...
          vector2(k1:k2) = matmul(elem%MassInv%Mb(1:dof,1:dof), vector(k1:k2) )

          ! evaluation of the second order derivatives in integ nodes
          do l=1,Qdof
             D2wi(l, k, j) = dot_product(vector2(k1:k2), Dphi(1:dof, j, l) )
             !write(20+j,'(20es12.4)') Fx(l,1:2), D2wi(l,k,j),Dwi(l,1,1:2), &
             !     dot_product(vector2(k1:k2), state%space%V_rule(elem%Qnum)%phi(1:dof, l) )
          enddo
       enddo
    enddo
    deallocate(vector, vector2 )



    ! setting of the residuum in integ nodes
    allocate(Dwt(1:Qdof, 1:ndim) )
    call Eval_w_t_Elem(elem, Dwt(1:Qdof, 1:ndim) )  ! time derivative of w_h in integ

    do k=1,ndim
       k1 = (k-1)*dof + 1
       k2 = k*dof

       ! -(w^k - w^{k-1})/tau + RHS + \epsilon*\Delta u_h
       f(k,1:Qdof) = Dwt(1:Qdof, k) -  f(k, 1:Qdof) &
            - state%model%Re1*( D2wi(1:Qdof, k, 1) + D2wi(1:Qdof, k, 2) )

       f(k,1:Qdof) = f(k,1:Qdof)**2   ! square for the L2 norm
    enddo

    ! integration of  over K, array rez(:) is overwritten in IntegrateVectorFunction
    call IntegrateVectorFunction(elem, f(1:ndim, 1:Qdof), rez(1:ndim) )
    rez(1:ndim) = elem%diam * rez(1:ndim)**0.5

    deallocate(Dwt, f, Fx)

  end subroutine ElementResiduum

  !> evaluate \f$ \int_{\partial K} [{\bf w}]^2\, dS \f$,
  subroutine ElementJumps(elem, Hjump, HjumpD, Hjump_1)
    type(element), intent(inout) :: elem
    real, dimension(1:ndim), intent(out) ::  Hjump, HjumpD, Hjump_1
    real, dimension(:), allocatable :: jumps
    real :: hG
    integer :: ie

    allocate(jumps(1:ndim) )

    Hjump(:) = 0.
    HjumpD(:) = 0.
    Hjump_1(:) = 0.

    do ie = 1, elem%flen
       hG = elem%dn(ie)

       ! jumps of the solution
       call IntegElementEdgeJump(elem, ie, jumps(1:ndim) )
       Hjump(1:ndim)   = Hjump(1:ndim)   + (jumps(1:ndim) * hG )**0.5
       Hjump_1(1:ndim) = Hjump_1(1:ndim) + (jumps(1:ndim) / hG )**0.5
       !write(*,'(a6,2i5,12es14.6)') 'jumpSS',elem%i,ie,jumps(:)


       call IntegElementEdgeJumpDer(elem, ie, jumps(1:ndim) )
       HjumpD(1:ndim) = HjumpD(1:ndim) + (jumps(1:ndim) * hG )**0.5
       !write(*,'(a6,2i5,12es14.6)') 'jumpDD',elem%i,ie,jumps(:),&
       !     ((3.4641*4/2**0.5)**2)*(2**0.5/2)

    enddo

    deallocate(jumps)
  end subroutine ElementJumps

end module helmholtz_estim

