!> setting of matrix elements per blocks
!>
!> i.e., evaluating of various volume and/or edge integrals
!> of test functions
module blocks_integ

  use define_state ! not used
  use tquadrature_mod  ! not used
  use time_mod  ! not used
  use main_data
  use f_mapping

  implicit none

  public:: Eval_V_Weights
  public:: Eval_V_Weights_plus
  public:: Eval_Dphi
  public:: Eval_Dphi_plus
  public:: Eval_Dphi_Edge
  public:: Eval_Phi_Edge
  public:: EvalScalarProdEdge

  public:: IntegrateFunction
  public:: IntegrateVectorFunction
  public:: IntegrateVectorFunction2
  public:: IntegrateFunctionsVec
  public:: IntegrateSquareVectorFunction
  public:: IntegrateSquareVectorFunction2

  public:: IntegrateFunctionEdge
  public:: IntegrateFunctionVecEdge
  public:: IntegrateFunctionsEdgeVec

  public:: IntegrateFunctionNormalEdge
  public:: IntegrateSquareFunctionNormalEdge
  public:: IntegrateVectorB
  public:: IntegrateVectorBplus
  public:: IntegrateVectorDplus
  public:: IntegrateBlockBB
  public:: IntegrateBlockBBplus
  public:: IntegrateBlockBBmass
  public:: IntegrateBlockD2
  public:: IntegrateBlockD2plus
  public:: IntegrateBlockJumps
  public :: IntegrateTimeFunction
  public :: IntegrateTimeFunctionVec

  public:: IntegrateEdgeBlockBB

  public:: EvalBlockBB
  public:: EvalBlockBD
  public:: EvalBlockDB
  public:: EvalBlockDD
  public:: EvalVectorB
  public:: EvalVectorD
  public:: EvalEdgeBlockDiagBB
  public:: EvalEdgeBlockBB
  public:: EvalEdgeBlockBD
  public:: EvalEdgeBlockDB
  public:: EvalEdgeVectorB
  public:: EvalEdgeVectorD

  public:: ExplEvalD
  public:: ExplEdgeB
contains

  !> evaluation of weights for volume quadrature multiplied by Jacobian,
  !>
  !>  (5 -elem%type)  = \f$ |\hat{K}|\f$
  pure subroutine Eval_V_Weights(elem, weights)
    type(element), intent(in):: elem        ! elem = element
    real, dimension(1:elem%Qdof), intent(out):: weights

    if(elem%F%iFlin) then ! linear element, constant Jacobian
       weights(1:elem%Qdof)  = state%space%V_rule(elem%Qnum)%weights(1:elem%Qdof) &
            * elem%F%JF0 / (5 -elem%type)
    else
       weights(1:elem%Qdof)  = state%space%V_rule(elem%Qnum)%weights(1:elem%Qdof) &
            * elem%F%V%JF(1:elem%Qdof)/ (5 -elem%type)
    endif
  end subroutine Eval_V_Weights

  !> evaluation of weights for volume quadrature multiplied by Jacobian,
  !> with given quadrature rule
  !>  (5 -elem%type)  = \f$ |\hat{K}|\f$
  subroutine Eval_V_Weights_plus(elem, V_rule, weights)
    type(element), intent(in):: elem        ! elem = element
    type(volume_rule), intent(in) :: V_rule
    real, dimension(1: V_rule%Qdof), intent(inout):: weights
    real, dimension(:,:,:), allocatable :: DF !, D1F
    real, dimension(:), allocatable :: JF
    integer :: Qdof, Qnum

    Qdof = V_rule%Qdof  ! /= elem%Qdof !!!! in GENERAL

    if(elem%F%iFlin) then ! linear element, constant Jacobian
       weights(1:Qdof)  = V_rule%weights(1:Qdof) * elem%F%JF0 / (5 -elem%type)
    else
       ! OLD variant
       !weights(1:Qdof)  = state%space%V_rule(Qnum)%weights(1:Qdof) &
       !     * elem%F%JF0 / (5 -elem%type)

       allocate( DF(1:Qdof, 1:nbDim, 1:nbDim), JF(1:Qdof) )

       call ComputeDF(elem, Qdof, V_rule%lambda(1:Qdof, 1:nbDim), &
            DF(1:Qdof, 1:nbDim, 1:nbDim) )

       JF(1:Qdof) = DF(1:Qdof,1,1)*DF(1:Qdof,2,2) - DF(1:Qdof,1,2)*DF(1:Qdof,2,1)

       weights(1:Qdof)  = V_rule%weights(1:Qdof) * JF(1:Qdof)/(5-elem%type)

       deallocate(DF, JF)

    endif
  end subroutine Eval_V_Weights_plus


  !> evaluation of the basis functions on edge of element
  subroutine Eval_Phi_Edge(elem, dofA, ie, phi, opposite)
    type(element), intent(in):: elem        ! elem = element
    integer, intent(in) :: dofA           ! actual dof, usually dofA=dof
    integer, intent(in) :: ie               ! index of the edge
    ! derivetives of the test function in integ. nodes on edge ie
    real, dimension(1:dofA,1:elem%face(fGdof,ie)), intent(inout):: phi
    logical, intent(in) :: opposite   ! opposite orientation of integ. nodes on  edge
    real, dimension(:,:), pointer:: Rphi ! pointers to test functions
    integer :: Qnum, Qdof, l

    Qnum = elem%face(fGnum,ie)
    Qdof = state%space%G_rule(Qnum)%Qdof


    ! the reference test functions
    if(elem%HGnode) then
       Rphi => state%space%G_rule(Qnum)%phi(elem%type, elem%HGface(1, ie), elem%HGface(2, ie), &
            1:dofA, 1:Qdof)
    else
       Rphi => state%space%G_rule(Qnum)%phi(elem%type, ie, 1, 1:dofA, 1:Qdof)
    endif

    if(state%print) then
       if(elem%HGnode) then
          do l=1,dofA
             write(*,'(a5,i5,3i3,20es10.2)')'PhiHG',elem%i, elem%HGface(1:nbDim, ie),l,&
                  Rphi(l,1:Qdof)
          enddo
       else
          do l=1,dofA
             write(*,'(a5,i5,3i3,20es10.2)')'Phi  ',elem%i, ie, 1,l,&
                  Rphi(l,1:Qdof)
          enddo
       endif
       print*
    endif

    ! if opposite then reorder
    if(opposite) then
       do l=1, Qdof
          phi(1:dofA,l) = Rphi(1:dofA, Qdof -l + 1)
       enddo
    else
       phi(1:dofA, 1:Qdof) = Rphi(1:dofA, 1:Qdof)
    endif
  end subroutine Eval_Phi_Edge


  !> evaluation of the derivation of the basis functions on element at elem integ nodes
  subroutine Eval_Dphi(elem, dofA,  Der)
    type(element), intent(in):: elem        ! elem = element
    integer, intent(in) :: dofA           ! actual dof, usually dofA=dof
    ! derivetives of the test function in integ. nodes on elem
    real, dimension(1:dofA,1:nbDim,1:elem%Qdof), intent(inout):: Der
    real, dimension(:,:,:), pointer:: Dphi ! pointers to test functions
    !real, dimension(:), allocatable:: temp
    integer :: Qnum, Qdof, j

    Qnum = elem%Qnum
    Qdof = elem%Qdof

    ! derivatives of the reference test functions
    Dphi => state%space%V_rule(Qnum)%Dphi(1:dofA, 1:nbDim, 1:Qdof)

    if(elem%F%iFlin) then        ! linear element
       do j=1, dofA
          Der(j, 1, 1:Qdof) = &
               elem%F%D1F0(1,1) * Dphi(j, 1, 1:Qdof) &
               +elem%F%D1F0(1,2)* Dphi(j, 2, 1:Qdof)

          Der(j, 2, 1:Qdof) = &
               elem%F%D1F0(2,1) * Dphi(j, 1, 1:Qdof) &
               +elem%F%D1F0(2,2)* Dphi(j, 2, 1:Qdof)
       enddo
    else                         ! curved element
       do j=1, dofA
          Der(j, 1, 1:Qdof) =  &
               elem%F%V%D1F(1:Qdof, 1,1) * Dphi(j, 1, 1:Qdof) &
               +elem%F%V%D1F(1:Qdof, 1,2)* Dphi(j, 2, 1:Qdof)

          Der(j, 2, 1:Qdof) = &
               elem%F%V%D1F(1:Qdof, 2,1) * Dphi(j, 1, 1:Qdof) &
               +elem%F%V%D1F(1:Qdof, 2,2)* Dphi(j, 2, 1:Qdof)
       enddo
    endif

  end subroutine Eval_Dphi


  !> evaluation of the derivation of the basis functions on element at integ nodes
  !> of the given quadrature
  subroutine Eval_Dphi_plus(elem, V_rule, dofA,  Der)
    type(element), intent(in):: elem        ! elem = element
    type(volume_rule), target, intent(in) :: V_rule
    integer, intent(in) :: dofA           ! actual dof, usually dofA=dof
    ! derivetives of the test function in integ. nodes on elem
    real, dimension(1:dofA,1:nbDim,1:V_rule%Qdof), intent(inout):: Der
    real, dimension(:,:,:), pointer:: Dphi ! pointers to test functions
    real, dimension(:,:,:), allocatable :: DF , D1F
    real, dimension(:), allocatable :: JF
    integer :: Qdof, j

    Qdof = V_rule%Qdof  ! /= elem%Qdof !!!! in GENERAL


    ! derivatives of the reference test functions
    Dphi => V_rule%Dphi(1:dofA, 1:nbDim, 1:Qdof)

    if(elem%F%iFlin) then        ! linear element
       do j=1, dofA
          Der(j, 1, 1:Qdof) = &
               elem%F%D1F0(1,1) * Dphi(j, 1, 1:Qdof) &
               +elem%F%D1F0(1,2)* Dphi(j, 2, 1:Qdof)

          Der(j, 2, 1:Qdof) = &
               elem%F%D1F0(2,1) * Dphi(j, 1, 1:Qdof) &
               +elem%F%D1F0(2,2)* Dphi(j, 2, 1:Qdof)
       enddo
    else                         ! curved element

       allocate( DF(1:Qdof, 1:nbDim, 1:nbDim),  D1F(1:Qdof, 1:nbDim, 1:nbDim), &
            JF(1:Qdof) )

       call ComputeDF(elem, Qdof, V_rule%lambda(1:Qdof, 1:nbDim), &
            DF(1:Qdof, 1:nbDim, 1:nbDim) )

       JF(1:Qdof) = DF(1:Qdof,1,1)*DF(1:Qdof,2,2) - DF(1:Qdof,1,2)*DF(1:Qdof,2,1)

       !transpose and inverse of DF/Dx
       D1F(1:Qdof,1,1) =  DF(1:Qdof,2,2) / JF(1:Qdof)
       D1F(1:Qdof,2,1) = -DF(1:Qdof,1,2) / JF(1:Qdof)
       D1F(1:Qdof,1,2) = -DF(1:Qdof,2,1) / JF(1:Qdof)
       D1F(1:Qdof,2,2) =  DF(1:Qdof,1,1) / JF(1:Qdof)

       do j=1, dofA
         Der(j, 1, 1:Qdof) =  &
              D1F(1:Qdof, 1,1) * Dphi(j, 1, 1:Qdof) &
              +D1F(1:Qdof, 1,2)* Dphi(j, 2, 1:Qdof)

         Der(j, 2, 1:Qdof) = &
              D1F(1:Qdof, 2,1) * Dphi(j, 1, 1:Qdof) &
              +D1F(1:Qdof, 2,2)* Dphi(j, 2, 1:Qdof)
       enddo
       deallocate (DF, D1F, JF)

       !write(*,'(a6,i5,12es12.4)') 'NEW',Qdof, Der(1:2, 1:2, 6:8)
       !print*
       !stop

    endif

  end subroutine Eval_Dphi_plus

  !> evaluation of the derivatives of the basis functions on edge of element
  !> multiplied by the weights of the Gauss quadrature
  subroutine Eval_Dphi_Edge(elem, dofA, ie, Der, opposite)
    type(element), intent(in):: elem        ! elem = element
    integer, intent(in) :: dofA           ! actual dof, usually dofA=dof
    integer, intent(in) :: ie               ! index of the edge
    ! derivetives of the test function in integ. nodes on edge ie
    real, dimension(1:dofA, 1:nbDim, 1:elem%face(fGdof,ie)), intent(inout):: Der
    logical, intent(in) :: opposite   ! opposite orientation of integ. nodes on edge
    real, dimension(:,:,:), pointer:: Dphi ! pointers to test functions
    real, dimension(:), allocatable:: temp
    integer :: Qnum, Qdof, j, l

    Qnum = elem%face(fGnum,ie)
    Qdof = state%space%G_rule(Qnum)%Qdof

    ! derivatives of the reference test functions
    if(elem%HGnode) then
       Dphi => state%space%G_rule(Qnum)%Dphi(elem%type, elem%HGface(1, ie), elem%HGface(2, ie), &
            1:dofA,1:nbDim, 1:Qdof)
    else
       Dphi => state%space%G_rule(Qnum)%Dphi(elem%type, ie, 1, 1:dofA, 1:nbDim, 1:Qdof)
    endif

    if(elem%F%iFlin) then        ! linear element
       do j=1, dofA
          Der(j, 1, 1:Qdof) = state%space%G_rule(Qnum)%weights(1:Qdof) &
               *(elem%F%D1F0(1,1) * Dphi(j, 1, 1:Qdof) &
               + elem%F%D1F0(1,2) * Dphi(j, 2, 1:Qdof))

          Der(j, 2, 1:Qdof) = state%space%G_rule(Qnum)%weights(1:Qdof) &
               *(elem%F%D1F0(2,1) * Dphi(j, 1, 1:Qdof) &
               + elem%F%D1F0(2,2) * Dphi(j, 2, 1:Qdof))
       enddo
    else                         ! curved element
       do j=1, dofA
          Der(j, 1, 1:Qdof) = state%space%G_rule(Qnum)%weights(1:Qdof) &
               *(elem%F%E(ie)%D1F(1:Qdof, 1,1) * Dphi(j, 1, 1:Qdof) &
               + elem%F%E(ie)%D1F(1:Qdof, 1,2) * Dphi(j, 2, 1:Qdof) )

          Der(j, 2, 1:Qdof) = state%space%G_rule(Qnum)%weights(1:Qdof) &
               *(elem%F%E(ie)%D1F(1:Qdof, 2,1) * Dphi(j, 1, 1:Qdof) &
               + elem%F%E(ie)%D1F(1:Qdof, 2,2) * Dphi(j, 2, 1:Qdof) )

       enddo
    endif

    ! if opposite then reorder
    if(opposite) then
       allocate(temp(1:Qdof) )

       do j=1,dofA
          ! first component
          do l=1, Qdof
             temp(l) = Der(j, 1, Qdof -l + 1)
          enddo
          Der(j, 1, 1:Qdof) = temp(1:Qdof)
          ! second component
          do l=1, Qdof
             temp(l) = Der(j, 2, Qdof -l + 1)
          enddo
          Der(j, 2, 1:Qdof) = temp(1:Qdof)
       enddo

       deallocate(temp)
    endif
  end subroutine Eval_Dphi_Edge



  !> integrate \f$ val =  \int_{K_{elem}} f \ dx \f$
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  subroutine IntegrateFunction(elem, func, val)
    type(element), intent (in) :: elem
    real, dimension(1:elem%Qdof), intent(in) :: func
    real, intent(inout) :: val
    real, dimension(:), allocatable :: weights
    integer :: Qdof

    Qdof = elem%Qdof

    allocate( weights(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof) )

    val = dot_product(weights(1:Qdof), func(1:Qdof) )

    deallocate(weights)

  end subroutine IntegrateFunction

  pure function IntegrateTimeFunction(Tdeg, tau, func) result (integral)
    integer, intent(in) :: Tdeg ! number of int. nodes
    real, intent(in) :: tau !length of the time interval
    real, dimension(1:Tdeg), intent(in) :: func !integrated func
    !class(Time_rule), pointer :: T_rule
    real :: integral


    !T_rule => state%time%T_rule(Tdeg)

    integral = tau * dot_product( state%time%T_rule(Tdeg)%weights(1:Tdeg), func(1:Tdeg) )

  end function IntegrateTimeFunction

     !> integral = \int_{I_m} f \ dt \f$
  pure function IntegrateTimeFunctionVec(Tdeg, tau, func) result (integral)
    integer, intent(in) :: Tdeg ! number of int. nodes
    real, intent(in) :: tau !length of the time interval
    real, dimension(1:ndim, 1:Tdeg), intent(in) :: func !integrated func
    !class(Time_rule), pointer :: T_rule
    real, dimension(1:ndim) :: integral


    !T_rule => state%time%T_rule(Tdeg)

    integral(1:ndim) = tau * &
       matmul( func(1:ndim, 1:Tdeg), state%time%T_rule(Tdeg)%weights(1:Tdeg) )

  end function IntegrateTimeFunctionVec

   !> sca_prod = \int_{I_m} f q \ dt \f$
  function EvalTimeScalarProduct(Tdeg, tau, func1, func2) result (sca_prod)
    integer, intent(in) :: Tdeg ! number of int. nodes
    real, intent(in) :: tau !length of the time interval
    real, dimension(1:Tdeg), intent(in) :: func1, func2 !integrated functions
    real, dimension(1:Tdeg) :: temp
    class(Time_rule), pointer :: T_rule
    real :: sca_prod
    real :: Tdof

    T_rule => state%time%T_rule(Tdeg)

    temp(1:tdeg) = func1(1:Tdeg) * func2(1:Tdeg)

    sca_prod = tau * dot_product( T_rule%weights(1:Tdeg), temp(1:Tdeg) )

  end function EvalTimeScalarProduct



  !> integrate \f$ val =  \int_{K_{elem}} {\bf f} \cdot {\bf q} \ dx \f$
  !>
  !> \f$ {\bf f}, {\bf g} \in R^d \f$ are  given by its values in integration nodes 1..Qdof
  subroutine IntegrateFunctionsVec(elem, f, g, val)
    type(element), intent (in) :: elem
    real, dimension(1:nbDim, 1:elem%Qdof), intent(in) :: f
    real, dimension(1:elem%Qdof, 1:nbDim), intent(in) :: g
    real, intent(out) :: val
    real, dimension(:,:), allocatable :: weights
    integer :: Qdof, i

    Qdof = elem%Qdof

    allocate( weights(1:2, 1:Qdof) )
    call Eval_V_Weights(elem, weights(1, 1:Qdof) )

    weights(2, 1:Qdof) = 0.
    do i=1,nbDim
       weights(2, 1:Qdof) = weights(2, 1:Qdof) + f(i, 1:Qdof) * g(1:Qdof, i)
    enddo

    val = dot_product(weights(1, 1:Qdof), weights(2, 1:Qdof) )

    deallocate(weights)

  end subroutine IntegrateFunctionsVec


  !> integrate \f$ val =  \int_{K_{elem}} {\bf f} \cdot {\bf q} \ dx \f$
  !>
  !> \f$ {\bf f}, {\bf g} \in R^d \f$ are  given by its values in integration nodes 1..Qdof
  subroutine IntegrateFunctionsVec_plus(elem, V_rule, f, g, val)
    type(element), intent (in) :: elem
    type(volume_rule), intent(in) :: V_rule
    real, dimension(1:nbDim, 1:V_rule%Qdof), intent(in) :: f
    real, dimension(1:V_rule%Qdof, 1:nbDim), intent(in) :: g
    real, intent(inout) :: val
    real, dimension(:,:), allocatable :: weights
    integer :: Qdof, i

    Qdof = V_rule%Qdof

    allocate( weights(1:2, 1:Qdof) )
    call Eval_V_Weights_plus(elem, V_rule, weights(1, 1:Qdof) )

    weights(2, 1:Qdof) = 0.
    do i=1,nbDim
       weights(2, 1:Qdof) = weights(2, 1:Qdof) + f(i, 1:Qdof) * g(1:Qdof, i)
    enddo

    val = dot_product(weights(1, 1:Qdof), weights(2, 1:Qdof) )

    deallocate(weights)

  end subroutine IntegrateFunctionsVec_plus

   !> integrate \f$ val =  \int_{\Gamma} {\bf f} \cdot {\bf q} \ dx \f$
  !>
  !> \f$ {\bf f}, {\bf g} \in R^d \f$ are  given by its values in integration nodes 1..Gdof
  subroutine IntegrateFunctionsEdgeVec(elem, iedge, Gnum, Gdof, f, g, integral)
    type(element), intent (in) :: elem
    integer, intent(in) :: iedge
    integer, intent(in) :: Gnum
    integer, intent(in) :: Gdof
    real, dimension(1:nbDim, 1:Gdof), intent(in) :: f
    real, dimension(1:Gdof, 1:nbDim), intent(in) :: g
    real, intent(out) :: integral
    real, dimension(:), allocatable :: weights
    real, dimension(:), allocatable :: temp
    integer :: i

    if (Gdof /= state%space%G_rule(Gnum)%Qdof ) &
      stop 'IntegrateFunctionsEdgeVec wrong Gnum and Gdof?'

    allocate( weights(1:Gdof) )
    weights(1:Gdof)  = state%space%G_rule(Gnum)%weights(1:Gdof)

   !TODO:
   ! write(debug, *) 'FR - isnt wrong? what if this edge is not the curved one?'
!    if(elem%ibcur > 0) then
    if  ( elem%ibcur > 0 .and. elem%jcur == iedge ) then
     print*, 'FR - IntegrateFunctionsEdgeVec - is it right? '
       weights(1:Gdof) = weights(1:Gdof) * elem%dnc(1:Gdof)
    else
       weights(1:Gdof) = weights(1:Gdof) * elem%dn(iedge)
    endif


    allocate( temp(1:Gdof) , source = 0.0 )
    do i=1,nbDim
       temp(1:Gdof) = temp(1:Gdof) + ( f(i, 1:Gdof) * g(1:Gdof, i) )
    enddo

    integral = dot_product( weights(1:Gdof), temp(1:Gdof) )

    deallocate(weights, temp)

  end subroutine IntegrateFunctionsEdgeVec


  !> integrate \f$ val =  \int_{K_{elem}} f \ dx \f$
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  subroutine IntegrateVectorFunction(elem, func, val)
    type(element), intent (in) :: elem
    real, dimension(1:ndim, 1:elem%Qdof), intent(in) :: func
    real, dimension(1:ndim), intent(out) :: val
    real, dimension(:), allocatable :: weights
    integer :: Qdof

    Qdof = elem%Qdof

    allocate( weights(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof) )

    val(1:ndim) = matmul(func(1:ndim,1:Qdof), weights(1:Qdof)  )

    deallocate(weights)

  end subroutine IntegrateVectorFunction


  !> integrate \f$ val =  \int_{K_{elem}} f \ dx \f$
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  subroutine IntegrateVectorFunction2(elem, func, val)
    type(element), intent (in) :: elem
    real, dimension(1:elem%Qdof, 1:ndim), intent(in) :: func
    real, dimension(1:ndim), intent(out) :: val
    real, dimension(:), allocatable :: weights
    integer :: Qdof

    Qdof = elem%Qdof

    allocate( weights(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof) )

    val(1:ndim) = matmul(weights(1:Qdof), func(1:Qdof, 1:ndim)   )

    deallocate(weights)

  end subroutine IntegrateVectorFunction2

  !> integrate \f$ val =  \int_{K_{elem}} f^2 \ dx \f$
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  subroutine IntegrateSquareVectorFunction(elem, func, val)
    type(element), intent (in) :: elem
    real, dimension(1:ndim, 1:elem%Qdof), intent(in) :: func
    real, dimension(1:ndim), intent(out) :: val
    real, dimension(:), allocatable :: weights
    integer :: Qdof

    Qdof = elem%Qdof

    allocate( weights(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof) )

    val(1:ndim) = matmul(func(1:ndim,1:Qdof)*func(1:ndim,1:Qdof), weights(1:Qdof)  )

    deallocate(weights)

  end subroutine IntegrateSquareVectorFunction


  !> integrate \f$ val =  \int_{K_{elem}} f^2 \ dx \f$
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  pure subroutine IntegrateSquareVectorFunction2(elem, func, val)
    type(element), intent (in) :: elem
    real, dimension(1:elem%Qdof, 1:ndim), intent(in) :: func
    real, dimension(1:ndim), intent(out) :: val
    real, dimension(:), allocatable :: weights
    integer :: Qdof

    Qdof = elem%Qdof

    allocate( weights(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof) )

    val(1:ndim) = matmul(weights(1:Qdof), func(1:Qdof, 1:ndim) * func(1:Qdof,1:ndim)  )

    deallocate(weights)

  end subroutine IntegrateSquareVectorFunction2

  !> integrate \f$ val =  \int_{K_{elem}} f \ dx \f$
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  pure subroutine IntegrateFunctionEdge(elem, ie, func, val)
    type(element), intent (inout) :: elem
    integer, intent(in) :: ie
    !real, dimension(1:elem%Qdof), intent(in) :: func  ! HERE WAS ERROR !!!!
    real, dimension(1:elem%face(fGdof,ie)), intent(in) :: func
    real, intent(inout) :: val
    real, dimension(:), allocatable :: weights
    integer :: Qdof, Qnum

    Qnum = elem%face(fGnum,ie)
    Qdof = elem%face(fGdof,ie)

    allocate( weights(1:Qdof) )


    if( elem%ibcur > 0 .and. elem%jcur == ie ) then
      !print*, 'Control curved edge in EvalScalarProdEdge - is it right?'
       !print*, 'FR IntegrateFunctionEdge - isnt wrong? what if this edge is not the curved one?'
       weights(1:Qdof) = state%space%G_rule(Qnum)%weights(1:Qdof) * elem%dnc(1:Qdof)
    else
       weights(1:Qdof) = state%space%G_rule(Qnum)%weights(1:Qdof) * elem%dn(ie)
    endif

    val = dot_product(weights(1:Qdof) , func(1:Qdof) )

    deallocate(weights)

  end subroutine IntegrateFunctionEdge

    !> integrate \f$ val =  \int_{K_{elem}} f \ dx \f$
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  pure subroutine IntegrateFunctionVecEdge(elem, ie, d, func, val)
    type(element), intent (inout) :: elem
    integer, intent(in) :: ie
    integer, intent(in) :: d !dimension
    real, dimension(1:d, 1:elem%face(fGdof,ie)), intent(in) :: func
    real, dimension(1:d), intent(inout) :: val

    real, dimension(:), allocatable :: weights
    integer :: Qdof, Qnum

    Qnum = elem%face(fGnum,ie)
    Qdof = elem%face(fGdof,ie)

    allocate( weights(1:Qdof) )


    if( elem%ibcur > 0 .and. elem%jcur == ie ) then
      !print*, 'FR IntegrateFunctionEdge - isnt wrong? what if this edge is not the curved one?'
       weights(1:Qdof) = state%space%G_rule(Qnum)%weights(1:Qdof) * elem%dnc(1:Qdof)
    else
       weights(1:Qdof) = state%space%G_rule(Qnum)%weights(1:Qdof) * elem%dn(ie)
    endif

    val(1:d) = matmul( func(1:d, 1:Qdof), weights(1:Qdof) )

    deallocate(weights)

  end subroutine IntegrateFunctionVecEdge




    !> integrate \f$ val =  \int_{K_{elem}} f g \ dx \f$
  !>
  !> \f$ f,g \f$ = func are given by their values in integration nodes 1..Qdof
  pure function EvalScalarProdEdge(elem, ie, func1, func2) result (sca_prod)
    type(element), intent (in) :: elem
    integer, intent(in) :: ie
    real, dimension(1:elem%face(fGdof,ie)), intent(in) :: func1, func2
    real :: sca_prod

    real, dimension(:), allocatable :: weights
    integer :: Qdof, Qnum

    Qnum = elem%face(fGnum,ie)
    Qdof = elem%face(fGdof,ie)

    allocate( weights(1:Qdof) )

    if( elem%ibcur > 0 .and. elem%jcur == ie ) then
!      print*, 'Control curved edge in EvalScalarProdEdge - is it right?'
      weights(1:Qdof) = state%space%G_rule(Qnum)%weights(1:Qdof) * elem%dnc(1:Qdof)
    else
      weights(1:Qdof) = state%space%G_rule(Qnum)%weights(1:Qdof) * elem%dn(ie)
    endif

      sca_prod = dot_product( weights(1:Qdof) , func1(1:Qdof)*func2(1:Qdof) )

   deallocate( weights )

  end function EvalScalarProdEdge




  !> integrate \f$ val =  \int_{K_{elem}} {\bf f}(x)\cdot{\bf n}(x) q(x) \ dx \f$
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  subroutine IntegrateFunctionNormalEdge(elem, ie, vf, q, val)
    type(element), intent (in) :: elem
    integer, intent(in) :: ie
    real, dimension(1:nbDim, 1:elem%face(fGdof,ie)), intent(in) :: vf
    real, dimension(1:elem%face(fGdof,ie)), intent(in) :: q
    real, intent(out) :: val
    real, dimension(:,:), allocatable :: temp
    integer :: Qdof, Qnum, i

    Qnum = elem%face(fGnum,ie)
    Qdof = elem%face(fGdof,ie)

    if (elem%i == 6) then
        !print*, 'v------------In IntegrateFunctionNormalEdge sub.----------------v'
        !write(*,'(a20,2i3)') 'ie, Qnum:', ie, Qnum
        !write(*,'(a20,2i3)') 'ie, Qdof:', ie, Qdof
        !print*, ' '
        !write(*,'(a12, 8es14.6)') 'vf(1, 1:Qdof):', vf(1, 1:Qdof)
        !write(*,'(a12, 8es14.6)') 'vf(2, 1:Qdof):', vf(2, 1:Qdof)
        !print*, ' '
        !write(*,'(a12, 8es14.6)') 'elem%n(ie, 1):', elem%n(ie, 1)
        !write(*,'(a12, 8es14.6)') 'elem%n(ie, 2):', elem%n(ie, 2)
    endif

    allocate( temp(1:2, 1:Qdof) )
    temp(1, 1:Qdof)  = state%space%G_rule(Qnum)%weights(1:Qdof) * q(1:Qdof)
    temp(2, 1:Qdof)  = 0.

    if( elem%ibcur > 0 .and. elem%jcur == ie ) then
       do i=1,nbDim
          temp(2, 1:Qdof)  = temp(2, 1:Qdof) + vf(i, 1:Qdof) * elem%nc(1:Qdof, i)
       enddo
    else
       do i=1,nbDim
          temp(2, 1:Qdof)  = temp(2, 1:Qdof) + vf(i, 1:Qdof) * elem%n(ie, i)
       enddo
    endif

!    if ((elem%i == 1 .or. elem%i == 3) .and. ie == 1) then
!        print*, 'v------------In IntegrateFunctionNormalEdge sub.----------------v'
!        print*, 'QI:' , q(1:Qdof)
!        print*, 'weight:',state%space%G_rule(Qnum)%weights(1:Qdof)
!        write(*,'(a15,8es14.6)') 'temp(1, 1:Qdof):', temp(1, 1:Qdof)
!        !write(*,'(a15,8es14.6)') 'temp(2, 1:Qdof):', temp(2, 1:Qdof)
!        !print*, ' '
!    endif

    val = dot_product(temp(1, 1:Qdof), temp(2, 1:Qdof) )

    deallocate(temp)

  end subroutine IntegrateFunctionNormalEdge


  !> integrate \f$ val =  \int_{ie} ({\bf f}(x)\cdot{\bf n}(x))^2 \ dS \f$
  !> \f$ ie \f$ = (sub)edge of a simplex elem
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  subroutine IntegrateSquareFunctionNormalEdge(elem, ie, vf, val)
    type(element), intent (in) :: elem
    integer, intent(in) :: ie
    real, dimension(1:nbDim, 1:elem%face(fGdof,ie)), intent(in) :: vf
    real, intent(inout) :: val
    real, dimension(:,:), allocatable :: temp
    integer :: Qdof, Qnum, i, l

    Qnum = elem%face(fGnum,ie)
    Qdof = elem%face(fGdof,ie)

    allocate( temp(1:2, 1:Qdof) )
    temp(1, 1:Qdof)  = state%space%G_rule(Qnum)%weights(1:Qdof)
    temp(2, 1:Qdof)  = 0.

    if(elem%ibcur > 0) then
       do i=1,nbDim
          temp(2, 1:Qdof)  = temp(2, 1:Qdof) + vf(i, 1:Qdof) * elem%nc(1:Qdof, i)
       enddo
    else
       do i=1,nbDim
          temp(2, 1:Qdof)  = temp(2, 1:Qdof) + vf(i, 1:Qdof) * elem%n(ie, i)
       enddo
    endif

    do l=1, Qdof
       temp(2, l) = (temp(2, l))**2
    enddo

    val = dot_product(temp(1, 1:Qdof), temp(2, 1:Qdof) )

    deallocate(temp)

  end subroutine IntegrateSquareFunctionNormalEdge

  !> integrate \f$ Mblock(R,C) =  \sum_{\Gamma \in F_h} \frac{1}{|\Gamma|}
  !> \int_{\Gamma}  \nabla\phi_{R}\cdot\nabla\phi_C\ dx \f$
  !> \f$ \forall R, C \f$
  subroutine IntegrateBlockJumps(elem, dofA, Mblock)
    type(element), intent (inout) :: elem
    integer, intent(in) :: dofA           ! actual dof, usually dofA=dof
!    real, dimension(1:elem%Qdof), intent(in) :: func
    real, dimension(1:dofA, 1:dofA), intent(inout) :: Mblock
    real, dimension(:), allocatable :: temp  !> weights in integ. nodes
    real, dimension(:, : ), allocatable :: phi
    integer :: ie, i, j, Qdof, Qnum


    do ie=1,elem%flen
       Qnum = elem%face(fGnum,ie)
       Qdof = state%space%G_rule(Qnum)%Qdof

       allocate(phi(1:dofA, 1:Qdof) )
       call  Eval_Phi_Edge(elem, dofA, ie, phi, .false. )

       allocate(temp(1:Qdof) )
       do i=1, dofA
          temp(1:Qdof) = state%space%G_rule(Qnum)%weights(1:Qdof) * phi(i, 1:Qdof)

          do j=1, dofA
             Mblock(i, j) = Mblock(i, j)  &
                  + dot_product(temp(1:Qdof), phi(j, 1:Qdof) )
          enddo
       enddo

       deallocate(temp, phi)

    enddo
  end subroutine IntegrateBlockJumps

  !> integrate \f$ Vector(R) =  \int_{K_{elem}} f \phi_{R} \ dx \f$
  !> \f$ \forall R \f$
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  subroutine IntegrateVectorB(elem, dofA, func, Vector)
    type(element), intent (in) :: elem
    integer, intent(in) :: dofA
    real, dimension(1:elem%Qdof), intent(in) :: func
    real, dimension(1:dofA), intent(inout) :: Vector
    real, dimension(:), allocatable :: weights
    real, dimension(:, :), pointer :: phi
    integer :: Qdof

    Qdof = elem%Qdof

    phi => state%space%V_rule(elem%Qnum)%phi(1:dofA, 1:Qdof)

    allocate( weights(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof) )

    weights(1:Qdof) = weights(1:Qdof) * func(1:Qdof)

    Vector(1:dofA) = matmul(phi(1:dofA, 1:Qdof), weights(1:Qdof) )

    deallocate(weights)

  end subroutine IntegrateVectorB


  !> integrate \f$ Vector(R) =  \int_{K_{elem}} f \phi_{R} \ dx \f$
  !> \f$ \forall R \f$ with a given quadrature
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  subroutine IntegrateVectorBplus(elem, Qnum, dofA, func, Vector)
    type(element), intent (in) :: elem
    integer, intent(in) :: Qnum, dofA
    real, dimension(1:state%space%V_rule(Qnum)%Qdof), intent(in) :: func
    real, dimension(1:dofA), intent(inout) :: Vector
    real, dimension(:), allocatable :: weights
    real, dimension(:, :), pointer :: phi
    integer ::  Qdof

    Qdof = state%space%V_rule(Qnum)%Qdof

    phi => state%space%V_rule(Qnum)%phi(1:dofA, 1:Qdof)

    allocate( weights(1:Qdof) )
    call Eval_V_Weights_plus(elem, state%space%V_rule(Qnum), weights(1:Qdof) )

    weights(1:Qdof) = weights(1:Qdof) * func(1:Qdof)
    Vector(1:dofA) = matmul(phi(1:dofA, 1:Qdof), weights(1:Qdof) )

    deallocate(weights)

  end subroutine IntegrateVectorBplus


  !> integrate \f$ Vector(R) =  \int_{K_{elem}} \nabla f \nabla \phi_{R} \ dx \f$
  !> \f$ \forall R \f$ with a given quadrature
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  subroutine IntegrateVectorDplus(elem, Qnum, dofA, func, Vector)
    type(element), intent (in) :: elem
    integer, intent(in) :: Qnum, dofA
    real, dimension(1:state%space%V_rule(Qnum)%Qdof, 1:nbDim), intent(in) :: func
    real, dimension(1:dofA), intent(inout) :: Vector
    real, dimension(:), allocatable :: weights
    real, dimension(:, :, :), pointer :: Der
    type(volume_rule), pointer :: V_rule
    integer ::  i, n, Qdof

    V_rule => state%space%V_rule(Qnum)
    Qdof = V_rule%Qdof

    allocate( Der(1:dofA, 1:nbDim, 1:Qdof) )
    call  Eval_Dphi_plus(elem, V_rule, dofA, Der(1:dofA, 1:nbDim, 1:Qdof) )

    allocate( weights(1:Qdof) )
    call Eval_V_Weights_plus(elem, V_rule, weights(1:Qdof) )

    do i=1, dofA
       do n=1,nbDim
          Vector(i) = Vector(i) + &
               dot_product(weights(1:Qdof), Der(i,n,1:Qdof) * func(1:Qdof, n) )
       enddo
    enddo ! loop i

    deallocate(weights, Der)

  end subroutine IntegrateVectorDplus


  !> integrate \f$ Mblock(R,C) =  \int_{K_{elem}} f \phi_{R} \phi_C\ dx \f$
  !> \f$ \forall R, C \f$
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..elem%Qdof
  subroutine IntegrateBlockBB(elem, dofA, func, Mblock)
    type(element), intent (inout) :: elem
    integer, intent(in) :: dofA
    real, dimension(1:elem%Qdof), intent(in) :: func
    real, dimension(1:dofA, 1:dofA), intent(inout) :: Mblock
    real, dimension(:), allocatable :: weights, temp
    real, dimension(:, :), pointer :: phi
    integer :: i, Qdof

    Qdof = elem%Qdof

    phi => state%space%V_rule(elem%Qnum)%phi(1:dofA, 1:Qdof)

    allocate(weights(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof) )

    weights(1:Qdof) = weights(1:Qdof) * func(1:Qdof)

    allocate(temp(1:Qdof) )
    do i=1, dofA
       temp(1:Qdof) = weights(1:Qdof) * phi(i, 1:Qdof)
       Mblock(i, 1:dofA) = Mblock(i, 1:dofA) + matmul(phi(1:dofA, 1:Qdof), temp(1:Qdof) )
    enddo

    deallocate(weights, temp)

  end subroutine IntegrateBlockBB


  !> integrate \f$ Mblock(R,C) =  \int_{K_{elem}} f \phi_{R} \phi_C\ dx \f$
  !> \f$ \forall R, C \f$ with the given integ rule
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..elem%Qdof
  subroutine IntegrateBlockBBplus(elem, Qnum, dofA, func, Mblock)
    type(element), intent (inout) :: elem
    integer, intent(in) :: Qnum, dofA
    real, dimension(1:state%space%V_rule(Qnum)%Qdof), intent(in) :: func
    real, dimension(1:dofA, 1:dofA), intent(inout) :: Mblock
    real, dimension(:), allocatable :: weights, temp
    real, dimension(:, :), pointer :: phi
    integer :: i, Qdof

    Qdof = state%space%V_rule(Qnum)%Qdof

    phi => state%space%V_rule(Qnum)%phi(1:dofA, 1:Qdof)

    allocate(weights(1:Qdof) )
    call Eval_V_Weights_plus(elem, state%space%V_rule(Qnum), weights(1:Qdof) )

    weights(1:Qdof) = weights(1:Qdof) * func(1:Qdof)

    allocate(temp(1:Qdof) )
    do i=1, dofA
       temp(1:Qdof) = weights(1:Qdof) * phi(i, 1:Qdof)
       Mblock(i, 1:dofA) = Mblock(i, 1:dofA) + matmul(phi(1:dofA, 1:Qdof), temp(1:Qdof) )
    enddo

    deallocate(weights, temp)

  end subroutine IntegrateBlockBBplus

  !calculates rectangular matrix dofA x dofB
  !> integrate \f$ Mblock(R,C) =  \int_{K_{elem}} f \phi_{R} \phi_C\ dx \f$
  !> \f$ \forall R, C \f$ with the given integ rule
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..elem%Qdof
  subroutine IntegrateBlockBBmass(elem, Qnum, dofA, dofB, Mblock)
    type(element), intent (inout) :: elem
    integer, intent(in) :: Qnum
    integer, intent(in) :: dofA , dofB
  !  real, dimension(1:state%space%V_rule(Qnum)%Qdof), intent(in) :: func
    real, dimension(1:dofA, 1:dofA), intent(inout) :: Mblock
    real, dimension(:), allocatable :: weights, temp
    real, dimension(:, :), pointer :: phi
    integer :: i, Qdof



    Qdof = state%space%V_rule(Qnum)%Qdof

    phi => state%space%V_rule(Qnum)%phi(1:dofA, 1:Qdof)

    allocate(weights(1:Qdof) )
    call Eval_V_Weights_plus(elem, state%space%V_rule(Qnum), weights(1:Qdof) )

  !  weights(1:Qdof) = weights(1:Qdof) * func(1:Qdof)

    allocate(temp(1:Qdof) )
    do i=1, dofA
       temp(1:Qdof) = weights(1:Qdof) * phi(i, 1:Qdof)
       Mblock(i, 1:dofB) = Mblock(i, 1:dofB) + matmul(phi(1:dofB, 1:Qdof), temp(1:Qdof) )
    enddo

    deallocate(weights, temp)

  end subroutine IntegrateBlockBBmass

  !> integrate \f$ Mblock(R,C) =
  !> \int_{K_{elem}} f \nabla\phi_{R}\cdot\nabla\phi_C\ dx \f$
  !> \f$ \forall R, C \f$
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  subroutine IntegrateBlockD2(elem, dofA, func, Mblock)
    type(element), intent (inout) :: elem
    integer, intent(in) :: dofA           ! actual dof, usually dofA=dof
    real, dimension(1:elem%Qdof), intent(in) :: func
    real, dimension(1:dofA, 1:dofA), intent(inout) :: Mblock
    real, dimension(:), allocatable :: weights, temp1, temp2  !> weights in integ. nodes
    real, dimension(:, :, :), allocatable :: Der
    integer :: i, j, Qdof

    Qdof = elem%Qdof

    allocate(Der(1:dofA, 1:nbDim, 1:Qdof) )
    call  Eval_Dphi(elem, dofA, Der)

    allocate(weights(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof) )

    weights(1:Qdof) = weights(1:Qdof) * func(1:Qdof)

    allocate(temp1(1:Qdof), temp2(1:Qdof) )
    do i=1, dofA
       temp1(1:Qdof) = weights(1:Qdof)* Der(i, 1, 1:Qdof)
       temp2(1:Qdof) = weights(1:Qdof)* Der(i, 2, 1:Qdof)

       do j=1, dofA
          Mblock(i, j) = Mblock(i, j)  &
               + dot_product(temp1(1:Qdof), Der(j, 1, 1:Qdof) ) &
               + dot_product(temp2(1:Qdof), Der(j, 2, 1:Qdof) )
       enddo
    enddo

    deallocate(weights, Der, temp1, temp2)

  end subroutine IntegrateBlockD2


  !> integrate \f$ Mblock(R,C) =
  !> \int_{K_{elem}} f \nabla\phi_{R}\cdot\nabla\phi_C\ dx \f$
  !> \f$ \forall R, C \f$
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  subroutine IntegrateBlockD2plus(elem, Qnum, dofA, func, Mblock)
    type(element), intent (inout) :: elem
    integer, intent(in) :: Qnum, dofA           ! actual dof, usually dofA=dof
    real, dimension(1:state%space%V_rule(Qnum)%Qdof), intent(in) :: func
    real, dimension(1:dofA, 1:dofA), intent(inout) :: Mblock
    real, dimension(:), allocatable :: weights, temp1, temp2  !> weights in integ. nodes
    real, dimension(:, :, :), allocatable :: Der
    type(volume_rule), pointer :: V_rule
    integer :: i, j, Qdof

    V_rule => state%space%V_rule(Qnum)
    Qdof = V_rule%Qdof

    allocate(Der(1:dofA, 1:nbDim, 1:Qdof) )
    call  Eval_Dphi_plus(elem, V_rule, dofA, Der(1:dofA, 1:nbDim, 1:Qdof) )

    allocate(weights(1:Qdof) )
    call Eval_V_Weights_plus(elem, state%space%V_rule(Qnum), weights(1:Qdof) )

    weights(1:Qdof) = weights(1:Qdof) * func(1:Qdof)

    allocate(temp1(1:Qdof), temp2(1:Qdof) )
    do i=1, dofA
       temp1(1:Qdof) = weights(1:Qdof)* Der(i, 1, 1:Qdof)
       temp2(1:Qdof) = weights(1:Qdof)* Der(i, 2, 1:Qdof)

       do j=1, dofA
          Mblock(i, j) = Mblock(i, j)  &
               + dot_product(temp1(1:Qdof), Der(j, 1, 1:Qdof) ) &
               + dot_product(temp2(1:Qdof), Der(j, 2, 1:Qdof) )
       enddo
    enddo

    deallocate(weights, Der, temp1, temp2)

  end subroutine IntegrateBlockD2plus


  !> integrate \f$Mblock(R,C)
  !> =  \int_{\Gamma} f\ \phi_R|_{elemR}\  \phi_C|_{elemC}\ dS \qquad \forall R, C\f$
  !>
  !> \f$f\f$ = func  given by its values in integ nodes 1..Qdof,
  !> \f$f\f$ IS ALREADY PRE-MULTIPLIED BY A SEGMENT INCREASE (=size of the face)
  !> \f$\Gamma \subset K_{elemR}\cap K_{elemC}\f$,
  !> R = "row" test function, "C" = "column" test function
  subroutine IntegrateEdgeBlockBB(elemR, ieR, elemC, ieC, func, Mblock)
    type(element), intent (inout) :: elemR, elemC
    integer, intent (in) :: ieR, ieC   ! inedex of the edge
    real, dimension(1:elemR%face(fGdof,ieR)), intent(in) :: func
    real, dimension(1:elemR%dof, 1:elemC%dof), intent(inout) :: Mblock
    real, dimension(:), allocatable :: temp, temp1
    !real, dimension(:, :), pointer :: phiR, phi_temp
    real, dimension(:, :), allocatable :: phiR, phiC
    integer :: i, j, Qnum, Qdof
    integer :: elemRdof, elemCdof

    !! seting of degree of the Gauss quadrature
    Qnum =elemR%face(fGnum,ieR)
    Qdof =elemR%face(fGdof,ieR)

    elemRdof = elemR%dof
    elemCdof = elemC%dof

    !elemRdof = max(1, (elemR%deg + 0)*(elemR%deg + 1) /2)
    !elemCdof = max(1, (elemC%deg + 0)*(elemC%deg + 1) /2)
    !print*,'!!!',elemRdof, elemR%dof

    ! test functions on "rows" element
    allocate(phiR( 1:elemR%dof, 1:Qdof) )
    call Eval_Phi_Edge(elemR, elemR%dof, ieR, phiR, .false.)

    ! test functions on "columns" element
    allocate(phiC( 1:elemC%dof, 1:Qdof) )
    if(elemR%i == elemC%i) then  ! identical elements
       call Eval_Phi_Edge(elemC, elemC%dof, ieC, phiC, .false.)
    else  ! different elements, oposite orientation of interation nodes
       call Eval_Phi_Edge(elemC, elemC%dof, ieC, phiC, .true.)
    endif


    allocate(temp(1:Qdof), temp1(1:Qdof) )
    temp(1:Qdof) = func(1:Qdof) * state%space%G_rule(Qnum)%weights(1:Qdof)

    do i=1, elemRdof
       temp1(1:Qdof) = temp(1:Qdof) * phiR(i, 1:Qdof)
       do j=1, elemC%dof
          Mblock(i, j) = Mblock(i, j) + dot_product(temp1(1:Qdof), phiC(j, 1:Qdof) )
       enddo
    enddo

    deallocate(temp, temp1)
    deallocate(phiC, phiR)

  end subroutine IntegrateEdgeBlockBB





!!! EVAL SUBROUTINES  -- big blocks

  !> global variant:
  !> integrate \f$ Vector(R) =  \int_{K_{elem}} f \phi_{R} \ dx \f$
  !> \f$ \forall R \f$
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  subroutine EvalVectorB(elem, func, dofA, Vector)
    type(element), intent (inout) :: elem
    real, dimension(1:elem%Qdof,1:ndim), intent(in) :: func
    integer, intent(in) :: dofA     ! dimension of Vector
    real, dimension(1:dofA *ndim), intent(inout) :: Vector
    real, dimension(:), allocatable :: temp, weights
    real, dimension(:, :), pointer :: phi
    integer :: Qdof, Qnum, i, k, kst

    Qdof = elem%Qdof
    Qnum = elem%Qnum

    phi => state%space%V_rule(Qnum)%phi(1:dofA, 1:Qdof)

    allocate( weights(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof) )

    allocate( temp(1:Qdof) )

    do k=1,ndim
       temp(1:Qdof) = weights(1:Qdof) * func(1:Qdof,k)
       kst = dofA*(k-1)
       do i=1, dofA
          Vector(kst+i) = Vector(kst+i) + dot_product(temp(1:Qdof), phi(i, 1:Qdof) )
          !if(i==1) then
          !   write(*,'(a6,i5,26es12.4)') 'func:',i, func(:,1)
          !   write(*,'(a6,i5,26es12.4)') 'temp:',i, temp(:)
          !   write(*,'(a6,i5,26es12.4)') 'phi: ',i, phi(i,:)
          !   print*, 'ddot:', dot_product(temp(1:Qdof), phi(i, 1:Qdof) )
          !   print*,'_____________________________________________'
          !endif

       enddo ! loop i
    enddo ! loop k

    deallocate(temp, weights)

  end subroutine EvalVectorB


  !> global variant:
  !> integrate \f$ Vector(R) =
  !> \int_{K_{elem}} (f_1 \partial_1 \phi_{R} + f_2 \partial_2 \phi_{R} )\ dx \f$
  !> \f$ \forall R \f$
  !>
  !> \f$ f \f$ = func is given by its values in integration nodes 1..Qdof
  subroutine EvalVectorD(elem, func, dofA, Vector)
    type(element), intent (inout) :: elem
    real, dimension(1:elem%Qdof, 1:nbDim, 1:ndim), intent(in) :: func
    integer, intent(in) :: dofA     ! dimension of Vector
    real, dimension(1:dofA*ndim), intent(inout) :: Vector
    real, dimension(:), allocatable :: weights
    real, dimension(:, :, :), allocatable :: Der
    !real, dimension(:, :), pointer :: phi
    integer :: Qdof, Qnum, i, k, kst

    Qdof = elem%Qdof
    Qnum = elem%Qnum

    allocate( weights(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof) )

    allocate(Der(1:dofA, 1:nbDim, 1:Qdof))
    call Eval_Dphi(elem, dofA, Der)

    do k=1,ndim
       kst = dofA*(k-1)
       do i=1, dofA
          Vector(kst+i) = Vector(kst+i) + &
               dot_product(weights(1:Qdof), Der(i,1,1:Qdof) * func(1:Qdof, 1, k) &
                  + Der(i,2,1:Qdof) *  func(1:Qdof, 2, k) )
       enddo ! loop i
    enddo ! loop k

    deallocate(Der, weights)

  end subroutine EvalVectorD


  !> global variant:
  !> integrate \f$ Mblock(R,C) =  \int_{K_{elem}}  f
  !> \phi_{C} \phi_R\ dx \f$
  !> \f$ \forall R, C \f$
  !>
  !> \f$ f = \f$ = func is given by its values
  !> in integration nodes 1..Qdof
  subroutine EvalBlockBB(elem, func, Mblock)
    type(element), intent (inout) :: elem
    real, dimension(1:elem%Qdof, ndim, ndim), intent(in) :: func
    real, dimension(1:ndim*elem%dof, 1:ndim*elem%dof), intent(inout) :: Mblock
    real, dimension(:), allocatable :: weights
    real, dimension(:, :), pointer :: phi
    integer :: Qdof, i, j, k, k1, col, row

    Qdof = elem%Qdof

    phi => state%space%V_rule(elem%Qnum)%phi(1:elem%dof, 1:Qdof)

    allocate(weights(1:Qdof))
    call Eval_V_Weights(elem, weights(1:Qdof) )


    do i=1, elem%dof
       do k=1,ndim                 ! k = index of component of w
          row = (k-1)*elem%dof
          do k1=1,ndim          ! k1 = index of component of w
             col = (k1-1)*elem%dof
             do j=1, elem%dof
                Mblock(row+i, col+j) = Mblock(row+i, col+j)  &
                     + dot_product( weights(1:Qdof)* phi(i, 1:Qdof) , &
                      func(1:Qdof, k, k1) *  phi(j, 1:Qdof)  )

             enddo ! loop j
          enddo !loop k1
       enddo ! loop k
    enddo  ! loop i

    deallocate(weights)

  end subroutine EvalBlockBB


  !> global variant:
  !> integrate \f$ Mblock(R,C) =  \int_{K_{elem}} \sum_{s=1}^2 f_{s}
  !> \phi_{C} \partial_s \phi_R\ dx \f$
  !> \f$ \forall R, C \f$
  !>
  !> \f$ f = \{f_{s}\}_{s=1}^2 \f$ = func is given by its values
  !> in integration nodes 1..Qdof
  subroutine EvalBlockBD(elem, func, Mblock)
    type(element), intent (inout) :: elem
    real, dimension(1:elem%Qdof, 1:nbDim, ndim, ndim), intent(in) :: func
    real, dimension(1:ndim*elem%dof, 1:ndim*elem%dof), intent(inout) :: Mblock
    real, dimension(:), allocatable :: weights, temp
    real, dimension(:, :, :), allocatable :: Der
    real, dimension(:, :), pointer :: phi
    integer :: Qdof, i, j, k, k1, col, row

    Qdof = elem%Qdof

    phi => state%space%V_rule(elem%Qnum)%phi(1:elem%dof, 1:Qdof)

    allocate(weights(1:Qdof), temp(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof) )

    allocate(Der(1:elem%dof, 1:nbDim, 1:Qdof))
    call Eval_Dphi(elem, elem%dof, Der)

    do i=1, elem%dof
       do k=1,ndim                 ! k = index of component of w
          row = (k-1)*elem%dof
          do k1=1,ndim          ! k1 = index of component of w
             col = (k1-1)*elem%dof

             temp(1:Qdof) = weights(1:Qdof) *(Der(i,1,1:Qdof) * func(1:Qdof, 1, k, k1) &
                  + Der(i,2,1:Qdof) *  func(1:Qdof, 2, k, k1) )
             !write(*,*) '%%%$$$', func(1:Qdof, 1, k, k1)  ! wet steam
             do j=1, elem%dof
                Mblock(row+i, col+j) = Mblock(row+i, col+j) &
                + dot_product(temp(1:Qdof), phi(j, 1:Qdof) )
             enddo ! loop j
          enddo !loop k1
       enddo ! loop k
    enddo  ! loop i

    deallocate(temp, weights )
    deallocate(Der )

  end subroutine EvalBlockBD


  !> global variant:
  !> integrate \f$ Mblock(R,C) =  \int_{K_{elem}} \sum_{s=1}^2 f_{s}
  !> \partial_s \phi_{C} \phi_R\ dx \f$
  !> \f$ \forall R, C \f$
  !>
  !> \f$ f = \{f_{s}\}_{s=1}^2 \f$ = func is given by its values
  !> in integration nodes 1..Qdof
  subroutine EvalBlockDB(elem, func, Mblock)
    type(element), intent (inout) :: elem
    real, dimension(1:elem%Qdof, 1:nbDim, ndim, ndim), intent(in) :: func
    real, dimension(1:ndim*elem%dof, 1:ndim*elem%dof), intent(inout) :: Mblock
    real, dimension(:), allocatable :: weights, temp
    real, dimension(:, :, :), allocatable :: Der
    real, dimension(:, :), pointer :: phi
    integer :: Qdof, i, j, k, k1, col, row

    Qdof = elem%Qdof

    phi => state%space%V_rule(elem%Qnum)%phi(1:elem%dof, 1:Qdof)

    allocate(weights(1:Qdof), temp(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof) )

    allocate(Der(1:elem%dof, 1:nbDim, 1:Qdof))
    call Eval_Dphi(elem, elem%dof, Der)

    do i=1, elem%dof
       do k=1,ndim                 ! k = index of component of w
          row = (k-1)*elem%dof
          do k1=1,ndim          ! k1 = index of component of w
             col = (k1-1)*elem%dof

!             temp(1:Qdof) = weights(1:Qdof) *(Der(i,1,1:Qdof) * func(1:Qdof, 1, k, k1) &
!                  + Der(i,2,1:Qdof) *  func(1:Qdof, 2, k, k1) )

             do j=1, elem%dof

                Mblock(row+i, col+j) = Mblock(row+i, col+j)  &
                     + dot_product( weights(1:Qdof)* phi(i, 1:Qdof) , &
                     ( func(1:Qdof,1,k,k1) *  Der(j,1,1:Qdof)  &
                     + func(1:Qdof,2,k,k1) * Der(j,2,1:Qdof) ) )

!                Mblock(row+i, col+j) = Mblock(row+i, col+j) &
!                + dot_product(temp(1:Qdof), phi(j, 1:Qdof) )


             enddo ! loop j
          enddo !loop k1
       enddo ! loop k
    enddo  ! loop i

    deallocate(temp, weights )
    deallocate(Der )

  end subroutine EvalBlockDB


  !> global variant:
  !> integrate \f$ Mblock(R,C) =  \int_{K_{elem}} \sum_{k,s=1}^2 f_{s,k}
  !> \partial_k\phi_{C} \partial_s \phi_R\ dx \f$
  !> \f$ \forall R, C \f$
  !>
  !> \f$ f = \{f_{s,k}\}_{s,k=1}^2 \f$ = func is given by its values
  !> in integ nodes 1..Qdof
  subroutine EvalBlockDD(elem, func, Mblock)
    type(element), intent (inout) :: elem
    real, dimension(1:elem%Qdof,1:nbDim,1:nbDim,1:ndim, 1:ndim), intent(in) :: func
    real, dimension(1:ndim*elem%dof, 1:ndim*elem%dof), intent(inout) :: Mblock
    real, dimension(:), allocatable     :: weights     !> weights in integ nodes
    real, dimension(:, :, :), allocatable :: Der      !> derivatives of test functions
    integer :: Qdof, i, j,  k, k1, row, col

    Qdof = elem%Qdof

    ! setting of integ. weights
    allocate(weights(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof) )

    ! derivatives of test functions
    allocate(Der(1:elem%dof, 1:nbDim, 1:Qdof))
    call Eval_Dphi(elem, elem%dof, Der)

    do i=1, elem%dof
       ! evaluation of Block elements
       do j=1, elem%dof
          do k=1,ndim                 ! k = index of component of w
             row = (k-1)*elem%dof
             do k1=1,ndim          ! k1 = index of component of w
                col = (k1-1)*elem%dof


                Mblock(row+i, col+j) = Mblock(row+i, col+j)  &
                     + dot_product( weights(1:Qdof), &
                     func(1:Qdof,1,1,k,k1) * Der(i,1,1:Qdof) * Der(j,1,1:Qdof)  &
                     + func(1:Qdof,1,2,k,k1) * Der(i,1,1:Qdof) * Der(j,2,1:Qdof)  &
                     + func(1:Qdof,2,1,k,k1) * Der(i,2,1:Qdof) * Der(j,1,1:Qdof)  &
                     + func(1:Qdof,2,2,k,k1) * Der(i,2,1:Qdof) * Der(j,2,1:Qdof) )
             enddo ! loop k1
          enddo ! loop k
       enddo ! loop j
    enddo ! loop i

    deallocate(weights, Der )

  end subroutine EvalBlockDD

  !> global variant:
  !> integrate \f$Mblock(R,C)
  !> =  \int_{\Gamma} f\ \phi_R|_{elemR}\  \phi_C|_{elemC}\ dS \qquad \forall R, C\f$
  !>
  !> \f$f\f$ = func  given by its values in integ nodes 1..Qdof,
  !> \f$f\f$ IS ALREADY PRE-MULTIPLIED BY A SEGMENT INCREASE (=size of the face)
  !> \f$\Gamma \subset K_{elemR}\cap K_{elemC}\f$,
  !> R = "row" test function, "C" = "column" test function,
  !> ONLY DIAGONAL BLOCKS
  subroutine EvalEdgeBlockDiagBB(elemR, ieR, elemC, ieC, func, Mblock)
    type(element), intent (inout) :: elemR, elemC
    integer, intent (in) :: ieR, ieC   ! inedex of the edge
    real, dimension(1:elemR%face(fGdof,ieR) ), intent(in) :: func
    real, dimension(1:elemR%dof*ndim, 1:elemC%dof*ndim),intent(inout) :: Mblock
    real, dimension(:), allocatable :: temp, temp1
    !real, dimension(:, :), pointer :: phiR, phi_temp
    real, dimension(:, :), allocatable :: phiC, phiR
    integer :: Qdof, Qnum, i, j, k, k1, row, col

    !! seting of degree of the Gauss quadrature
    Qnum =elemR%face(fGnum,ieR)
    Qdof =elemR%face(fGdof,ieR)

    ! test functions on "rows" element
    allocate(phiR( 1:elemR%dof, 1:Qdof) )
    call Eval_Phi_Edge(elemR, elemR%dof, ieR, phiR, .false.)


    ! test functions on "columns" element
    allocate(phiC( 1:elemC%dof, 1:Qdof) )
    if(elemR%i == elemC%i) then  ! identical elements
       call Eval_Phi_Edge(elemC, elemC%dof, ieC, phiC, .false.)
    else  ! different elements, oposite orientation of interation nodes
       call Eval_Phi_Edge(elemC, elemC%dof, ieC, phiC, .true.)
    endif

    allocate(temp(1:Qdof), temp1(1:Qdof) )

    temp(1:Qdof) = func(1:Qdof) * state%space%G_rule(Qnum)%weights(1:Qdof)

    do k=1,ndim                 ! k = index of component of w
       row = (k-1)*elemR%dof

       !do k1=1,ndim          ! k1 = index of component of w
       k1 = k ! only diagonal blocks, penalization
       col = (k1-1)*elemC%dof

       do i=1, elemR%dof
          temp1(1:Qdof) = temp(1:Qdof) * phiR(i, 1:Qdof)
          do j=1, elemC%dof
             Mblock(row+i, col+j) = Mblock(row+i, col+j) &
                  + dot_product(temp1(1:Qdof), phiC(j, 1:Qdof) )
          enddo  ! loop j
       enddo   !loop i
    !enddo !loop k1
    enddo !loop k

    deallocate(temp, temp1)
    deallocate(phiR, phiC)

  end subroutine EvalEdgeBlockDiagBB

  !> similar as EvalEdgeBlockDiagBB, only for 2nd and 33rd component
  subroutine EvalEdgeBlockDiagBB23(elemR, ieR, elemC, ieC, func, Mblock)
    type(element), intent (inout) :: elemR, elemC
    integer, intent (in) :: ieR, ieC   ! inedex of the edge
    real, dimension(1:elemR%face(fGdof,ieR) ), intent(in) :: func
    real, dimension(1:elemR%dof*ndim, 1:elemC%dof*ndim),intent(inout) :: Mblock
    real, dimension(:), allocatable :: temp, temp1
    !real, dimension(:, :), pointer :: phiR, phi_temp
    real, dimension(:, :), allocatable :: phiC, phiR
    integer :: Qdof, Qnum, i, j, k, k1, row, col

    !! seting of degree of the Gauss quadrature
    Qnum =elemR%face(fGnum,ieR)
    Qdof =elemR%face(fGdof,ieR)

    ! test functions on "rows" element
    allocate(phiR( 1:elemR%dof, 1:Qdof) )
    call Eval_Phi_Edge(elemR, elemR%dof, ieR, phiR, .false.)


    ! test functions on "columns" element
    allocate(phiC( 1:elemC%dof, 1:Qdof) )
    if(elemR%i == elemC%i) then  ! identical elements
       call Eval_Phi_Edge(elemC, elemC%dof, ieC, phiC, .false.)
    else  ! different elements, oposite orientation of interation nodes
       call Eval_Phi_Edge(elemC, elemC%dof, ieC, phiC, .true.)
    endif

    allocate(temp(1:Qdof), temp1(1:Qdof) )

    temp(1:Qdof) = func(1:Qdof) * state%space%G_rule(Qnum)%weights(1:Qdof)

    !do k=1,ndim                 ! k = index of component of w
    do k=2,3
       row = (k-1)*elemR%dof

       !do k1=1,ndim          ! k1 = index of component of w
       k1 = k ! only diagonal blocks, penalization
       col = (k1-1)*elemC%dof

       do i=1, elemR%dof
          temp1(1:Qdof) = temp(1:Qdof) * phiR(i, 1:Qdof)
          do j=1, elemC%dof
             Mblock(row+i, col+j) = Mblock(row+i, col+j) &
                  + dot_product(temp1(1:Qdof), phiC(j, 1:Qdof) )
          enddo  ! loop j
       enddo   !loop i
    !enddo !loop k1
    enddo !loop k

    deallocate(temp, temp1)
    deallocate(phiR, phiC)

  end subroutine EvalEdgeBlockDiagBB23



  !> global variant:
  !> integrate \f$Mblock(R,C)
  !> =  \int_{\Gamma} f\ \phi_R|_{elemR}\  \phi_C|_{elemC}\ dS \qquad \forall R, C\f$
  !>
  !> \f$f\f$ = func  given by its values in integ nodes 1..Qdof,
  !> \f$f\f$ IS ALREADY PRE-MULTIPLIED BY A SEGMENT INCREASE (=size of the face)
  !> \f$\Gamma \subset K_{elemR}\cap K_{elemC}\f$,
  !> R = "row" test function, "C" = "column" test function
  subroutine EvalEdgeBlockBB(elemR, ieR, elemC, ieC, func, Mblock)
    type(element), intent (inout) :: elemR, elemC
    integer, intent (in) :: ieR, ieC   ! inedex of the edge
    real, dimension(1:elemR%face(fGdof,ieR), 1:ndim, 1:ndim), intent(in) :: func
    real, dimension(1:elemR%dof*ndim, 1:elemC%dof*ndim),intent(inout) :: Mblock
    real, dimension(:), allocatable :: temp, temp1
    !real, dimension(:, :), pointer :: phiR, phi_temp
    real, dimension(:, :), allocatable :: phiC, phiR
    integer :: Qdof, Qnum, i, j, k, k1, row, col

    !! seting of degree of the Gauss quadrature
    Qnum =elemR%face(fGnum,ieR)
    Qdof =elemR%face(fGdof,ieR)

    ! test functions on "rows" element
    allocate(phiR( 1:elemR%dof, 1:Qdof) )
    call Eval_Phi_Edge(elemR, elemR%dof, ieR, phiR, .false.)


    ! test functions on "columns" element
    allocate(phiC( 1:elemC%dof, 1:Qdof) )
    if(elemR%i == elemC%i) then  ! identical elements
       call Eval_Phi_Edge(elemC, elemC%dof, ieC, phiC, .false.)
    else  ! different elements, oposite orientation of interation nodes
       call Eval_Phi_Edge(elemC, elemC%dof, ieC, phiC, .true.)
    endif


    !if(elemC%i == 455 .or. elemC%i == 2564  .or. elemC%i == 2638  .or. elemC%i == 2639 &
    !     .or. elemR%i==455 .or. elemR%i==2564 .or. elemR%i==2638 .or. elemR%i == 2639)&
    !     then
    !   do l=1,elemC%dof
    !      write(*,'(a5,5i5,20es12.4)')'Dphi',elemC%i, elemR%i, ieC,ieR,l,phiC(l,1:Qdof)
    !   enddo
    !   print*,
    !   do l=1,elemR%dof
    !      write(*,'(a5,5i5,20es12.4)')'Dphi',elemC%i, elemR%i, ieC,ieR,l,phiR(l,1:Qdof)
    !   enddo
    !   print*,'_________________________________________'
    !endif

    allocate(temp(1:Qdof), temp1(1:Qdof) )

    do k=1,ndim                 ! k = index of component of w
       row = (k-1)*elemR%dof
       do k1=1,ndim          ! k1 = index of component of w
          col = (k1-1)*elemC%dof

          temp(1:Qdof) = func(1:Qdof,k,k1) * state%space%G_rule(Qnum)%weights(1:Qdof)

          do i=1, elemR%dof
             temp1(1:Qdof) = temp(1:Qdof) * phiR(i, 1:Qdof)
             do j=1, elemC%dof
                Mblock(row+i, col+j) = Mblock(row+i, col+j) &
                     + dot_product(temp1(1:Qdof), phiC(j, 1:Qdof) )
             enddo  ! loop j
          enddo   !loop i
       enddo !loop k1
    enddo !loop k

    deallocate(temp, temp1)
    deallocate(phiR, phiC)

  end subroutine EvalEdgeBlockBB



  !> global variant:
  !> integrate \f$Mblock(R,C)
  !> =  \int_{\Gamma}  \sum_{s=1}^2\left(\sum_{k=1}^2 f_{k,s}
  !> \ \partial_k\phi_R|_{elemR}\ \right)\, n_s\
  !> \phi_C|_{elemC}\ dS \qquad \forall R, C\f$
  !>
  !> \f$ f = \{f_{s,k}\}_{s,k=1}^2 \f$ = func is given by its values
  !> in integ nodes 1..Qdof,
  !> \f$\Gamma \subset K_{elemR}\cap K_{elemC}\f$,
  !> \f$ n = (n_1, n_2) =\f$ unit normal to \f$ \Gamma \f$,
  !> R = "row" test function, "C" = "column" test function
  subroutine EvalEdgeBlockBD(elemR, ieR, elemC, ieC, func, Mblock)
    type(element), intent (inout) :: elemR, elemC
    integer, intent (in) :: ieR, ieC   ! inedex of the edge
    real, dimension(1:elemR%face(fGdof,ieR),1:nbDim,1:nbDim,1:ndim,1:ndim),intent(in) :: func
    real, dimension(1:ndim*elemR%dof, 1:ndim*elemC%dof), intent(inout) :: Mblock
    real, dimension(:), allocatable ::  temp1
    real, dimension(:, :, :), allocatable :: DerR        ! test functions on elem
    !real, dimension(:, :), pointer :: phi_temp ! deriv. of test functions on elem
    real, dimension(:, :), allocatable :: phiC
    integer :: Qdof, i, j, k, k1, row, col

    !! seting of degree of the Gauss quadrature
    Qdof =elemR%face(fGdof,ieR)

    ! test functions on "column" element
    allocate(phiC( 1:elemC%dof, 1:Qdof) )
    if(elemR%i == elemC%i) then  ! identical elements
       call Eval_Phi_Edge(elemC, elemC%dof, ieC, phiC, .false.)
    else  ! different elements, oposite orientation of interation nodes
       call Eval_Phi_Edge(elemC, elemC%dof, ieC, phiC, .true.)
    endif


    ! derivatives of test functions on "row" element
    allocate(DerR(1:elemR%dof, 1:nbDim, 1:Qdof))
    call Eval_Dphi_Edge(elemR, elemR%dof, ieR, DerR, .false.)


    allocate( temp1(1:Qdof) )

    do j=1, elemR%dof
       do k=1,ndim                 ! k = index of component of w
          row = (k-1)*elemR%dof

          do k1=1,ndim          ! k1 = index of component of w
             col = (k1-1)*elemC%dof

             ! evaluation of < * >n terms
             if(elemR%face(neigh,ieR) > 0 .or. elemR%F%iFlin) then
                ! inner edge or boundary edge of linear element
                temp1(1:Qdof) = &
                     func(1:Qdof,1, 1, k, k1) * DerR(j,1,1:Qdof) * elemC%n(ieC, 1) &
                     + func(1:Qdof,1, 2, k, k1) * DerR(j,2,1:Qdof) * elemC%n(ieC, 1) &
                     + func(1:Qdof,2, 1, k, k1) * DerR(j,1,1:Qdof) * elemC%n(ieC, 2) &
!transp              + func(1:Qdof,2, 1, k, k1) * DerR(j,2,1:Qdof) * elemC%n(ieC, 1) &
!transp              + func(1:Qdof,1, 2, k, k1) * DerR(j,1,1:Qdof) * elemC%n(ieC, 2) &
                     + func(1:Qdof,2, 2, k, k1) * DerR(j,2,1:Qdof) * elemC%n(ieC, 2)
             else
                ! curved edge
                temp1(1:Qdof) = &
                     func(1:Qdof,1, 1, k, k1) * DerR(j,1,1:Qdof) * elemC%nc(1:Qdof, 1) &
                     + func(1:Qdof,1, 2, k, k1) * DerR(j,2,1:Qdof) * elemC%nc(1:Qdof, 1) &
                     + func(1:Qdof,2, 1, k, k1) * DerR(j,1,1:Qdof) * elemC%nc(1:Qdof, 2) &
! transp             + func(1:Qdof,2, 1, k, k1) * DerR(j,2,1:Qdof) * elemC%nc(1:Qdof, 1) &
!transp              + func(1:Qdof,1, 2, k, k1) * DerR(j,1,1:Qdof) * elemC%nc(1:Qdof, 2) &
                     + func(1:Qdof,2, 2, k, k1) * DerR(j,2,1:Qdof) * elemC%nc(1:Qdof, 2)
             endif

             do i=1,elemC%dof  ! multiplication by [*] term
                Mblock(row+j, col+i) = Mblock(row+ j, col+i) &
                     + dot_product(temp1(1:Qdof), phiC(i, 1:Qdof) )
             enddo

             ! nasledujici je pomalejsi nez-li predchozi
             !Mblock(row+j, col+1:col+elemC%dof) = Mblock(row+ j, col+1:col+elemC%dof) &
             !     + matmul(phiC(1:elemC%dof, 1:Qdof), temp1(1:Qdof) )


          enddo ! k1 loop
       enddo ! k loop

    enddo  ! j loop

    deallocate(temp1)
    deallocate(phiC)
    deallocate(DerR)

  end subroutine EvalEdgeBlockBD


  !> global variant:
  !> integrate \f$Mblock(R,C)
  !> =  \int_{\Gamma}  \sum_{s=1}^2\left(\sum_{k=1}^2 f_{s,k}
  !> \ \partial_k\phi_C|_{elemC}\ \right)\, n_s\
  !> \phi_R|_{elemR}\ dS \qquad \forall R, C\f$
  !>
  !> \f$ f = \{f_{s,k}\}_{s,k=1}^2 \f$ = func is given by its values
  !> in integ nodes 1..Qdof,
  !> \f$\Gamma \subset K_{elemR}\cap K_{elemC}\f$,
  !> \f$ n = (n_1, n_2) =\f$ unit normal to \f$ \Gamma \f$,
  !> R = "row" test function, "C" = "column" test function
  subroutine EvalEdgeBlockDB(elemR, ieR, elemC, ieC, func, Mblock)
    type(element), intent (inout) :: elemR, elemC
    integer, intent (in) :: ieR, ieC   ! inedex of the edge
    real, dimension(1:elemR%face(fGdof,ieR),1:nbDim,1:nbDim,1:ndim,1:ndim),intent(in) :: func
    real, dimension(1:elemR%dof*ndim, 1:elemC%dof*ndim), intent(inout) :: Mblock
    real, dimension(:), allocatable :: temp1
    real, dimension(:, :), allocatable :: phiR        ! test functions on elem
    real, dimension(:, :, :), allocatable :: DerC ! derivative of test functions
    integer :: Qdof, i, j, k, k1, row, col

    !! seting of degree of the Gauss quadrature
    Qdof =elemR%face(fGdof,ieR)

    ! test functions on "rows" element
    allocate(phiR( 1:elemR%dof, 1:Qdof) )
    call Eval_Phi_Edge(elemR, elemR%dof, ieR, phiR, .false.)

    ! derivatives of test functions on "columns" element
    allocate(DerC( 1:elemC%dof, 1:nbDim, 1:Qdof) )

    if(elemR%i == elemC%i) then  ! identical elements
       call Eval_Dphi_Edge(elemC, elemC%dof, ieC, DerC, .false.)
    else
       call Eval_Dphi_Edge(elemC, elemC%dof, ieC, DerC, .true.)
    endif


    allocate(temp1(1:Qdof) )

    do j=1, elemC%dof
       do k=1,ndim                 ! k = index of component of w
          row = (k-1)*elemR%dof
          do k1=1,ndim          ! k1 = index of component of w
             col = (k1-1)*elemC%dof

             if(elemR%face(neigh,ieR) > 0 .or. elemR%F%iFlin) then
                ! inner edge or boundary edge of linear element
                temp1(1:Qdof) = &
                     func(1:Qdof, 1, 1, k, k1) * DerC(j,1,1:Qdof) * elemR%n(ieR, 1) &
                     + func(1:Qdof, 1, 2, k, k1) * DerC(j,2,1:Qdof) * elemR%n(ieR, 1) &
                     + func(1:Qdof, 2, 1, k, k1) * DerC(j,1,1:Qdof) * elemR%n(ieR, 2) &
                     + func(1:Qdof, 2, 2, k, k1) * DerC(j,2,1:Qdof) * elemR%n(ieR, 2)
             else
                ! curved edge
                temp1(1:Qdof) = &
                     func(1:Qdof, 1, 1, k, k1) * DerC(j,1,1:Qdof) * elemR%nc(1:Qdof, 1) &
                     + func(1:Qdof, 1, 2, k, k1) * DerC(j,2,1:Qdof) * elemR%nc(1:Qdof, 1) &
                     + func(1:Qdof, 2, 1, k, k1) * DerC(j,1,1:Qdof) * elemR%nc(1:Qdof, 2) &
                     + func(1:Qdof, 2, 2, k, k1) * DerC(j,2,1:Qdof) * elemR%nc(1:Qdof, 2)
             endif

             do i=1,elemR%dof
                Mblock(row+i, col+j) = Mblock(row+i, col+j) &
                     + dot_product(temp1(1:Qdof), phiR(i, 1:Qdof) )
             enddo

             ! nasledujici je pomalejsi nez-li to predchozi !!! ???
             !Mblock(row+1:row+elemR%dof, col+j) = Mblock(row+1:row+elemR%dof, col+j) &
             !        + matmul(phiR(1:elemR%dof, 1:Qdof), temp1(1:Qdof) )


          enddo ! loop k1
       enddo ! loop k
    enddo ! loop j

    deallocate(DerC, phiR,  temp1)

  end subroutine EvalEdgeBlockDB


  !> global variant:
  !> integrate \f$ Vector(R)
  !> =  \int_{\Gamma}  f \phi_R|_{elemR} w_B\ dS \qquad \forall R\f$
  !>
  !> \f$ f = \f$  func is given by its values
  !> in integ nodes 1..Qdof,
  !> \f$f\f$ IS ALREADY PRE-MULTIPLIED BY A SEGMENT INCREASE (=size of the face)
  !> \f$ w_B = \f$ = wB is given by its values
  !> in integ nodes 1..Qdof,
  !> \f$\Gamma \subset K_{elemR}\cap K_{elemC}\f$,
  !> R = "row" test function, "C" = "column" test function
  subroutine EvalEdgeVectorB(elemR, ieR, func, wB, dofA, Vector)
    type(element), intent (inout) :: elemR
    integer, intent (in) :: ieR   ! inedex of the edge
    real, dimension(1:elemR%face(fGdof,ieR)), intent(in) :: func
    real, dimension(1:elemR%face(fGdof,ieR), 1:ndim), intent(in) :: wB
    integer, intent(in) :: dofA     ! dimension of Vector
    real, dimension(1:dofA*ndim), intent(inout) :: Vector
    real, dimension(:), allocatable :: temp
    real, dimension(:, :), allocatable :: phiR        ! test functions on elem

    integer :: Qdof, i, Qnum, k, row, elemRdof

    ! "size" of the edge
    ! func IS ALREADY PRE-MULTIPLIED BY A SEGMENT INCREASE (=size of the face)

    !! seting of degree of the Gauss quadrature
    Qnum = elemR%face(fGnum,ieR)
    Qdof = elemR%face(fGdof,ieR)

    elemRdof = elemR%dof
    !elemRdof = max(1, (elemR%deg + 0)*(elemR%deg + 1) /2)
    !print*,'!!!',elemRdof, elemR%dof

    ! test functions on "row" element
    allocate(phiR( 1:dofA, 1:Qdof) )
    call Eval_Phi_Edge(elemR, dofA, ieR, phiR, .false.)

    !phiR => state%space%G_rule(Qnum)%phi(elemR%type, ieR, 1, 1:elemRdof, 1:Qdof)

    allocate(temp(1:Qdof) )

    temp(1:Qdof) = state%space%G_rule(Qnum)%weights(1:Qdof) *func(1:Qdof)

    do k=1,ndim
       row = (k-1)*dofA

          do i=1,dofA
             Vector(row+i) = Vector(row+i) &
                  + dot_product(phiR(i,1:Qdof),  temp(1:Qdof)* wB(1:Qdof, k) )
          enddo

       ! nasledujici je pomalejsi nez-li to predchozi !!! ???
       !Vector(row+1:row+elemRdof) = Vector(row+1:row+elemRdof) &
       !     + matmul(phiR(1:elemRdof,1:Qdof),  temp(1:Qdof)* wB(1:Qdof, k) )
    enddo ! loop k

    deallocate(temp, phiR)

  end subroutine EvalEdgeVectorB


  !> global variant:
  !> integrate \f$ Vector(R)
  !> =  \int_{\Gamma}  \sum_{s=1}^2\left(\sum_{k=1}^2 f_{k,s}
  !> \ \partial_k\phi_R|_{elemR}\ \right)\, n_s\
  !> w_B\ dS \qquad \forall R\f$
  !>
  !> \f$ f = \{f_{s,k}\}_{s,k=1}^2 \f$ = func is given by its values
  !> in integ nodes 1..Qdof,
  !> \f$ w_B = \f$ = wB is given by its values
  !> in integ nodes 1..Qdof,
  !> \f$\Gamma \subset K_{elemR}\cap K_{elemC}\f$,
  !> \f$ n = (n_1, n_2) =\f$ unit normal to \f$ \Gamma \f$,
  !> R = "row" test function, "C" = "column" test function
  subroutine EvalEdgeVectorD(elemR, ieR, func, wB, dofA, Vector)
    type(element), intent (inout) :: elemR
    integer, intent (in) :: ieR   ! inedex of the edge
    real, dimension(1:elemR%face(fGdof,ieR),1:nbDim,1:nbDim,1:ndim,1:ndim),intent(in):: func
    real, dimension(1:elemR%face(fGdof,ieR), 1:ndim), intent(in) :: wB
    integer, intent(in) :: dofA     ! dimension of Vector
    real, dimension(1:dofA*ndim), intent(inout) :: Vector
    real, dimension(:), allocatable :: temp1
    real, dimension(:, :, :), allocatable :: DerR        ! test functions on elem
    integer ::  Qdof, j, k, k1, row

    !! seting of degree of the Gauss quadrature
    Qdof = elemR%face(fGdof,ieR)

    ! derivatives of test functions on "row" element
    allocate(DerR(1:dofA, 1:nbDim, 1:Qdof))

    call Eval_Dphi_Edge(elemR, dofA, ieR, DerR, .false.)

    allocate( temp1(1:Qdof) )

    do j=1, dofA
       do k=1,ndim                 ! k = index of component of w
          row = (k-1)*dofA
          do k1=1,ndim                 ! k = index of component of w
             if(elemR%F%iFlin .or. elemR%face(neigh,ieR)>0) then
                temp1(1:Qdof) = &
                     func(1:Qdof, 1, 1,k,k1) * DerR(j,1,1:Qdof) * elemR%n(ier, 1) &
                     + func(1:Qdof, 1, 2,k,k1) * DerR(j,2,1:Qdof) * elemR%n(ieR, 1) &
                     + func(1:Qdof, 2, 1, k, k1) * DerR(j,1,1:Qdof) * elemR%n(ieR, 2) &
!transp               + func(1:Qdof, 2, 1, k, k1) * DerR(j,2,1:Qdof) * elemR%n(ieR, 1) &
!transp               + func(1:Qdof, 1, 2, k, k1) * DerR(j,1,1:Qdof) * elemR%n(ieR, 2) &
                     + func(1:Qdof, 2, 2, k, k1) * DerR(j,2,1:Qdof) * elemR%n(ieR, 2)
             else
                temp1(1:Qdof) = &
                     func(1:Qdof, 1, 1, k, k1) * DerR(j,1,1:Qdof) * elemR%nc(1:Qdof, 1) &
                     + func(1:Qdof, 1, 2, k, k1) * DerR(j,2,1:Qdof) * elemR%nc(1:Qdof, 1) &
                     + func(1:Qdof, 2, 1, k, k1) * DerR(j,1,1:Qdof) * elemR%nc(1:Qdof, 2) &
!transp               + func(1:Qdof, 2, 1, k, k1) * DerR(j,2,1:Qdof) * elemR%nc(1:Qdof, 1)&
!transp               + func(1:Qdof, 1, 2, k, k1) * DerR(j,1,1:Qdof) * elemR%nc(1:Qdof, 2)&
                     + func(1:Qdof, 2, 2, k, k1) * DerR(j,2,1:Qdof) * elemR%nc(1:Qdof, 2)
             endif

             Vector(row+j) = Vector(row+j) + dot_product(temp1(1:Qdof), wB( 1:Qdof, k1) )

          enddo ! loop k1
       enddo ! loop k

    enddo ! loop j

    deallocate(temp1, DerR)

  end subroutine EvalEdgeVectorD


  !> global variant:
  !> integrate \f$ Vector(R) =  \int_{K_{elem}} \sum_{s=1}^2 f_{s}
  !> \partial_s \phi_R\ dx \f$
  !> \f$ \forall R \f$
  !>
  !> \f$ f = \{f_{s}\}_{s=1}^2 \f$ = func is given by its values
  !> in integration nodes 1..Qdof
  subroutine ExplEvalD(elem, func, dofA, Vector)
    type(element), intent (inout) :: elem
    real, dimension(1:elem%Qdof, 1:nbDim, ndim), intent(in) :: func
    integer, intent(in) :: dofA     ! dimension of Vector
    real, dimension(1: dofA *ndim), intent(inout) :: Vector
    real, dimension(:), allocatable :: weights, temp
    real, dimension(:, :, :), allocatable :: Der
    integer :: Qdof, i, k,  row

    Qdof = elem%Qdof

    allocate(weights(1:Qdof), temp(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof) )

    allocate(Der(1:dofA, 1:nbDim, 1:Qdof))
    call Eval_Dphi(elem,dofA, Der)

    do k=1,ndim                 ! k = index of component of w
       row = (k-1)*dofA
       do i=1, dofA

          Vector(row+i) = Vector(row+i) + dot_product(weights(1:Qdof), &
               Der(i,1,1:Qdof) * func(1:Qdof, 1, k) &
               + Der(i,2,1:Qdof) *  func(1:Qdof, 2, k) )
       enddo ! loop i
    enddo  ! loop k

    deallocate(temp, weights )
    deallocate(Der )

  end subroutine ExplEvalD


  !> global variant:
  !> integrate \f$Vector(R)
  !> =  \int_{\Gamma} f\ \phi_R|_{elem} \ dS \qquad \forall R\f$
  !>
  !> \f$f\f$ = func  given by its values in integ nodes 1..Qdof,
  !> \f$f\f$ IS ALREADY PRE-MULTIPLIED BY A SEGMENT INCREASE (=size of the face)
  !> \f$\Gamma \subset K_{elem}\f$,
  !> R = "row" test function,
  subroutine ExplEdgeB(elem, ie, func, dofA, Vector)
    type(element), intent (inout) :: elem
    integer, intent (in) :: ie  ! inedex of the edge
    real, dimension(1:elem%face(fGdof,ie), 1:ndim), intent(in) :: func
    integer, intent(in) :: dofA     ! dimension of Vector
    real, dimension(1:dofA*ndim),intent(inout) :: Vector
    real, dimension(:), allocatable :: temp
    real, dimension(:, :), allocatable :: phi
    integer :: Qdof, Qnum, i,  k,  row

    !! seting of degree of the Gauss quadrature
    Qnum =elem%face(fGnum,ie)
    Qdof =elem%face(fGdof,ie)

    ! test functions on "rows" element
    allocate(phi( 1:dofA, 1:Qdof) )
    call Eval_Phi_Edge(elem, dofA, ie, phi, .false.)
    !phi => state%space%G_rule(Qnum)%phi(elem%type, ie, 1, 1:elem%dof, 1:Qdof)

    allocate(temp(1:Qdof) )

    do k=1,ndim                 ! k = index of component of w
       row = (k-1)*dofA

       temp(1:Qdof) = func(1:Qdof,k) * state%space%G_rule(Qnum)%weights(1:Qdof)

       do i=1, dofA
          Vector(row+i) = Vector(row+i) + dot_product(temp(1:Qdof), phi(i, 1:Qdof) )
       enddo

    enddo !loop k

    deallocate(temp, phi)

  end subroutine ExplEdgeB


end module blocks_integ
