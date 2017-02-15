!> evaluate the solution on element and faces
module eval_sol
  use mesh_oper
  use main_data
  use data_mod
  use integration
  use blocks_integ
  use model_oper
  use basis

  implicit none

  public:: Eval_func_Elem
  public:: EvalFuncInIntNodes
  public:: Eval_vec_func_Elem
  public:: Eval_res_func_Elem
  public:: Eval_w_Elem
  public:: Eval_w_Elem_plus
  public:: Eval_w_Elem_time
  public:: Eval_w_Elem_at_time
  public:: Eval_w_t_Elem
  public:: Eval_w_w_Elem
  public:: Eval_Dw_Elem
  public:: Eval_Dw_Elem_plus
  public:: Eval_Dw_Elem_time
  public:: Eval_Dw_Elem_at_time
  public:: Eval_aver_w_Elem
  public:: Eval_aver_Dw_Elem
  public:: Eval_wht_Elem

!  public:: Eval_wAEE_Elem !ALG
!  public:: Eval_DwAEE_Elem !ALG


  public:: Eval_func_Edge
  public:: Eval_w_Edge
  public:: Eval_w_EdgeProj
  public:: Eval_w_Edge_time
  public:: Eval_aver_w_Edge
  public:: Eval_Dw_Edge

  public:: Eval_f_s_Elem
  public:: Eval_f_s_Elem_at_time
  public:: Eval_NFlux_Edge
  public:: Eval_R_s_Elem
  public:: Eval_R_s_Elem_at_time
  public:: Eval_R_s_Edge
  public:: Eval_aver_R_s_Edge

  public:: Eval_phi_Qnode
  public:: SetLagrCoeffs

  public:: PlotElemSolution
  public:: PlotElemSolution3D
  public:: PlotElemFunctionB
  public:: PlotElemFunctionQ

  public:: PlotElemFunction3D
  public:: PlotElem_D_Function3D

  public:: BarycCoord
  public:: Trans_Integ_Nodes_to_Basis

  public:: Energy_Elem_projection
  public:: Energy_Elem_deg_projection

  public:: Compute_means_values_DW
  public:: Eval_Canonical_Functions
  public:: Canonical_Functions
  public:: Set_canonical_basis_index_deriv

contains
  !> evaluation of the scalar function \f$ r \f$ on element \f$ elem \f$  in
  !> integ nodes, i.e. recomputation of elem%w into wi in volume integ. nodes
  subroutine Eval_func_Elem(elem, r, wi)
    type(element), intent(in):: elem      ! elem = element
    real, dimension(1:elem%dof), intent(in):: r ! r in integ. nodes
    real, dimension(1:elem%Qdof), intent(inout):: wi ! w in integ. nodes
    real, dimension(:,:), pointer:: phi ! local store arrays
    integer :: dof, Qdof, i

    dof = elem%dof
    Qdof = elem%Qdof

    phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)

    do i=1,Qdof
       wi(i) = dot_product(r(1:dof), phi(1:dof,i) )
    enddo

  end subroutine Eval_func_Elem

  !> evaluation of the scalar function \f$ r \f$ on element \f$ elem \f$  in
  !> integ nodes, i.e. recomputation function given as basis coefficients into wi in volume integ. nodes
  !> based on the fact that basis functions phi are the same for all degrees, only save in different nodes
  subroutine EvalFuncInIntNodes( V_rule, dof, func, wi )
    type( volume_rule ), intent(in) :: V_rule
    integer, intent(in) :: dof
    real, dimension(1:dof), intent(in):: func ! function with respect to the basis functions phi
    real, dimension(1:V_rule%Qdof), intent(out):: wi ! w in integ. nodes
!    real, dimension(:,:), pointer:: phi ! local store arrays
    integer :: i, Qdof

    Qdof = V_rule%Qdof
!    phi => V_rule%phi(1:dof,1:Qdof)

!    do i=1,Qdof
!       wi(i) = dot_product(r(1:dof), phi(1:dof,i) )
!    enddo

    wi(1:Qdof) = matmul( func(1:dof) , V_rule%phi(1:dof,1:Qdof) )

  end subroutine EvalFuncInIntNodes

  !> evaluation of the scalar function \f$ r \f$ on element \f$ elem \f$  in
  !> integ nodes, i.e. recomputation of elem%w into wi in volume integ. nodes
  subroutine Eval_vec_func_Elem(elem, dof, r, wi)
    type(element), intent(in):: elem      ! elem = element
    integer, intent(in) :: dof
    real, dimension(1:ndim, 1:dof), intent(in):: r ! r in integ. nodes
    real, dimension(1:ndim, 1:elem%Qdof), intent(inout):: wi ! w in integ. nodes
    real, dimension(:,:), pointer:: phi ! local store arrays
    integer :: Qdof, i

    Qdof = elem%Qdof

    phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)


    do i=1,Qdof
       wi(1:ndim,i) = matmul(r(1:ndim,1:dof), phi(1:dof,i) )
    enddo
  end subroutine Eval_vec_func_Elem

  !> evaluation of residual function \f$ r_T^i \f$ on element \f$ elem \f$ in
  !> integ nodes
  subroutine Eval_res_func_Elem(elem)
    type(element), intent(inout):: elem      ! elem = element
    real, dimension(:,:), pointer:: phi ! local store arrays
    real, dimension(:), pointer :: weights
    real, dimension(:,:), pointer :: MassInv
    integer :: dof, Qdof, i, k, kst, j, l
    real, dimension(:), allocatable :: test
    real, dimension(:), allocatable :: func
    real :: val

    dof = elem%dof
    Qdof = elem%Qdof

    allocate( func(1:Qdof) )

    MassInv => elem%MassInv%Mb

    do k=1, ndim              ! k = index of component of residual function
       kst = dof*(k-1) + 1
          elem%vec( res_func, kst:kst+dof-1) = matmul( MassInv(1:dof, 1:dof ), elem%vec( res_vec, kst:kst+dof-1) )

    !write(*,'(a6,2i5,50es12.4)') 'EW1',elem%i,dof, elem%vec( res_vec, kst:kst+dof-1)

    enddo

    phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)

    do k=1, ndim              ! k = index of component of residual function
       kst = dof*(k-1) + 1
       do i=1,Qdof
          elem%res_func(k, i) = dot_product(elem%vec( res_func, kst:kst+dof-1), phi(1:dof, i) )
       enddo
    enddo


    !TEST
    !allocate( test(1:dof) )
    !test(:) = 0
    func(:) = 0
    weights => state%space%V_rule(elem%Qnum)%weights(1:Qdof)
    do k=1, ndim
       kst = dof*(k-1) + 1
       l = kst
       do j=1, dof
          !do i=1, Qdof
          !   test(j) = test(j) + ( weights(i)*elem%res_func(k, i)*phi(j, i) )* elem%F%JF0 * 0.5
          !enddo
          !write(*,'(a5, i2, a1, a2, i2, a2, 2es16.8)') 'test', j, ',', 'R', l, ':', test(j), elem%vec( res_vec, l)
          func(1:Qdof) = elem%res_func(k, 1:Qdof)*phi(j, 1:Qdof)
          call IntegrateFunction(elem, func, val)
          !write(*,'(a5, i2, a1, a2, i2, a2, 2es16.8)') 'val', j, ',', 'R', l, ':', val, elem%vec( res_vec, l)
          l = l + 1
       enddo
    enddo

    deallocate( func )

    !deallocate(ipiv, ident, Mb )

  end subroutine Eval_res_func_Elem


  !> evaluation of the state vector \f$w\f$ on element \f$ elem \f$  in
  !> integ nodes, i.e. recomputation of elem%w into wi in volume integ. nodes
  subroutine Eval_w_Elem(elem, wi)
    type(element), intent(in):: elem      ! elem = element
    real, dimension(1:elem%Qdof,1:ndim), intent(inout):: wi ! w in integ. nodes
    real, dimension(:,:), pointer:: phi ! local store arrays
    integer :: dof, Qdof, k, kst, i

    dof = elem%dof
    Qdof = elem%Qdof

    phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)

    do k=1,ndim              ! k = index of component of w
       kst = dof*(k-1) + 1
       do i=1,Qdof
          wi(i,k) = dot_product(elem%w(0,kst:kst+dof-1), phi(1:dof,i) )
       enddo
       ! Faster ?
       !wi(1:Qdof, k) = matmul(elem%w(0,kst:kst+dof-1), phi(1:dof,1:Qdof) )
    enddo

  end subroutine Eval_w_Elem


  !> evaluation of the state vector \f$w\f$ on element \f$ elem \f$  in
  !> integ nodes, i.e. recomputation of elem%w into wi in volume integ. nodes
  subroutine Eval_w_Elem_plus(elem, V_rule, wi)
    type(element), intent(in):: elem      ! elem = element
    type(volume_rule), target, intent(in) :: V_rule
    real, dimension(1:V_rule%Qdof, 1:ndim ), intent(inout):: wi ! w in integ. nodes
    real, dimension(:,:), pointer:: phi ! local store arrays
    integer :: dof, Qdof, k, kst, i

    dof = elem%dof
    Qdof = V_rule%Qdof

    phi => V_rule%phi(1:dof,1:Qdof)

    do k=1,ndim              ! k = index of component of w
       kst = dof*(k-1) + 1

       wi(1:Qdof, k) = matmul(elem%w(0,kst:kst+dof-1),  phi(1:dof,1:Qdof) )
    enddo

  end subroutine Eval_w_Elem_plus

  !> evaluation of the state vector \f$w\f$ on element \f$ elem \f$  in
  !> integ nodes, i.e. recomputation of elem%w into wi in volume integ. nodes
  subroutine Eval_w_array_Elem_plus(elem, V_rule, dof, ww, wi)
    type(element), intent(in):: elem      ! elem = element
    type(volume_rule), target, intent(in) :: V_rule
    integer, intent(in) :: dof ! DOF of the array ww
    real, dimension ( 1:dof), intent(in) :: ww
    real, dimension(1:V_rule%Qdof ), intent(inout):: wi ! w in integ. nodes
    real, dimension(:,:), pointer:: phi ! local store arrays
    integer ::  Qdof,  i

    Qdof = V_rule%Qdof

    phi => V_rule%phi(1:dof,1:Qdof)

    wi(1:Qdof) = matmul(ww( 1:dof),  phi(1:dof,1:Qdof) )

  end subroutine Eval_w_array_Elem_plus

  !> evaluation of the state vector \f$w\f$ on element \f$ elem \f$  in
  !> integ nodes, i.e. recomputation of elem%w into wi in volume integ. nodes
  subroutine Eval_Dw_Elem_plus(elem, V_rule, Dwi)
    type(element), intent(in):: elem      ! elem = element
    type(volume_rule), target, intent(in) :: V_rule
    real, dimension(0:2, 1:ndim, 1:V_rule%Qdof ), intent(inout):: Dwi ! w in integ. nodes
    real, dimension(:,:), pointer:: phi ! local store arrays
    real, dimension(:,:,:), allocatable:: Dphi ! local store arrays
    integer :: dof, Qdof, k, kst, i

    dof = elem%dof
    Qdof = V_rule%Qdof

    allocate( Dphi(1:dof,1:nbDim, 1:Qdof) )

    phi => V_rule%phi(1:dof,1:Qdof)
    call Eval_Dphi_plus(elem, V_rule, dof,  Dphi(1:dof,1:nbDim, 1:Qdof)  )

    do k=1,ndim              ! k = index of component of w
       kst = dof*(k-1) + 1

       Dwi(0, k, 1:Qdof) = matmul(elem%w(0,kst:kst+dof-1),  phi(1:dof,1:Qdof) )
       Dwi(1, k, 1:Qdof) = matmul(elem%w(0,kst:kst+dof-1), Dphi(1:dof, 1, 1:Qdof) )
       Dwi(2, k, 1:Qdof) = matmul(elem%w(0,kst:kst+dof-1), Dphi(1:dof, 2, 1:Qdof) )
    enddo

    deallocate(Dphi)
  end subroutine Eval_Dw_Elem_plus


  !> evaluation of the state vector \f$w\f$ on element \f$ elem \f$  in
  !> integ nodes, i.e. recomputation of elem%w into wi in volume integ. nodes
  !> at  time  t_{k-1}
  subroutine Eval_w_Elem_time(elem, wi)
    type(element), intent(in):: elem      ! elem = element
    real, dimension(1:elem%Qdof,1:ndim), intent(inout):: wi ! w in integ. nodes
    real, dimension(:,:), pointer:: phi ! local store arrays
    integer :: dof, Qdof, k, kst, i

    dof = elem%dof
    Qdof = elem%Qdof

    phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)

    do k=1,ndim              ! k = index of component of w
       kst = dof*(k-1) + 1
       do i=1,Qdof
          wi(i,k) = dot_product(elem%w(1, kst:kst+dof-1), phi(1:dof,i) )
       enddo
    enddo

  end subroutine Eval_w_Elem_time


  !> evaluation of the state vector \f$w\f$ on element \f$ elem \f$  in
  !> integ nodes, i.e. recomputation of elem%w into wi in volume integ. nodes
  !> at  time  t, t=1 => t_k, t=0 => t_{k-1},  BDF ONLY
  subroutine Eval_w_Elem_at_time(elem, wi, t)
    type(element), intent(in):: elem      ! elem = element
    real, dimension(1:elem%Qdof,1:ndim), intent(inout):: wi ! w in integ. nodes
    real, intent(in) :: t
    real, dimension(:,:), pointer:: phi ! local store arrays
    integer :: dof, Qdof, k, kst, i

    dof = elem%dof
    Qdof = elem%Qdof

    phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)

    do k=1,ndim              ! k = index of component of w
       kst = dof*(k-1) + 1
       do i=1,Qdof
          wi(i,k) = dot_product(t*elem%w(0, kst:kst+dof-1) +(1-t)*elem%w(1, kst:kst+dof-1),&
               phi(1:dof,i) )
       enddo
    enddo

  end subroutine Eval_w_Elem_at_time

!  !> evaluation of the computational solution based on AEE ALG (state vector) \f$wc\f$ on element \f$ elem \f$  in
!  !> integ nodes, i.e. recomputation of elem%wc into wi in volume integ. nodes
!  subroutine Eval_wAEE_Elem(elem, wi)
!    type(element), intent(in):: elem      ! elem = element
!    real, dimension(1:elem%Qdof,1:ndim), intent(inout):: wi ! w in integ. nodes
!    real, dimension(:,:), pointer:: phi ! local store arrays
!    integer :: dof, Qdof, k, kst, i
!
!    dof = elem%dof
!    Qdof = elem%Qdof
!
!    phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)
!
!    do k=1,ndim              ! k = index of component of w
!       kst = dof*(k-1) + 1
!       do i=1,Qdof
!          wi(i,k) = dot_product(elem%wc(0, kst:kst+dof-1), phi(1:dof,i) )
!       enddo
!    enddo
!
!  end subroutine Eval_wAEE_Elem
!
!
!  !> evaluation of the derivatives of the computational solution based on AEE ALG (state vector) \f$wc\f$
!  !> on element \f$ elem \f$  in integ nodes
!  subroutine Eval_DwAEE_Elem(elem, Dwi)
!    type(element), intent(in) :: elem      ! elem = element
!    real, dimension(1:elem%Qdof,1:ndim,1:nbDim), intent(inout):: Dwi ! Dw in integ. nodes
!    real, dimension(:,:,:), allocatable :: Dphi
!    integer :: dof, Qdof, k, kst, i
!
!    dof = elem%dof
!    Qdof = elem%Qdof
!
!    ! derivatives on real element
!    allocate(Dphi(1:dof, 1:nbDim, 1:Qdof) )
!
!    call Eval_Dphi(elem, dof, Dphi)
!
!    do k=1,ndim              ! k = index of component of Dw
!       kst = dof*(k-1) + 1
!       do i=1,Qdof
!          Dwi(i, k, 1) = dot_product(elem%wc(0,kst:kst+dof-1), Dphi(1:dof, 1, i) )
!          Dwi(i, k, 2) = dot_product(elem%wc(0,kst:kst+dof-1), Dphi(1:dof, 2, i) )
!       enddo
!    enddo
!    deallocate(Dphi)
!  end subroutine Eval_DwAEE_Elem


  !> evaluation of the state vector \f${\bf w}_{h\tau}\f$ with its gradient
  !> on element \f$ elem \f$  in
  !> integ nodes, i.e. recomputation of elem%w into wi in volume integ. nodes
  !> at  time  t, given by the precomputed Lagrangain coefficients,
  !> \f${\bf w}_{h\tau}\f$ is piecewise polynomial with respect to time from Lagr. reconstr.
  subroutine Eval_wht_Elem(Tdeg, Lag_coef, elem, wi, Dwi)
    integer, intent(in) :: Tdeg    ! degree of Lagrangian reconstruction
    real, dimension(0:Tdeg), intent(in) :: Lag_coef
    type(element), intent(in):: elem      ! elem = element
    real, dimension(1:elem%Qdof,1:ndim), intent(inout):: wi ! w in integ. nodes
    real, dimension(1:elem%Qdof,1:ndim, 1:nbDim), intent(inout):: Dwi ! Dw in integ. nodes
    real, dimension(:,:), pointer:: phi ! local store arrays
    real, dimension(:), allocatable:: ws ! local store arrays
    real, dimension(:,:,:), allocatable :: Dphi
    integer :: dof, ndof, Qdof, k, kst, i

    dof = elem%dof
    ndof = dof * ndim
    Qdof = elem%Qdof
    allocate(ws (1:ndof) )

    phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)

    ! derivatives on real element
    allocate(Dphi(1:dof, 1:nbDim, 1:Qdof) )
    call Eval_Dphi(elem, dof, Dphi)

    ! setting of the Lagr. reconst. of the state vector at time t
    ws(1:ndof) = elem%w(0 , 1:ndof) * Lag_coef(Tdeg)  ! w
    do i=1, Tdeg
       ws(1:ndof) = ws(1:ndof) + elem%w(i, 1:ndof) * Lag_coef(Tdeg - i)  ! reverse ordering
    enddo

    !do i=0, Tdeg
    !   write(*,'(a6,2i5,20es16.8)') 'ws ',elem%i,i,elem%w(i, :)
    !enddo

    ! evaluation at volume integ. nodes
    do k=1,ndim              ! k = index of component of w
       kst = dof*(k-1) + 1
       do i=1,Qdof
          wi(i,k) = dot_product(ws(kst:kst+dof-1),  phi(1:dof,i) )

          Dwi(i, k, 1) = dot_product(ws(kst:kst+dof-1), Dphi(1:dof, 1, i) )
          Dwi(i, k, 2) = dot_product(ws(kst:kst+dof-1), Dphi(1:dof, 2, i) )

       enddo
    enddo

    deallocate(ws, Dphi )

  end subroutine Eval_wht_Elem




  !> evaluation of the state vector \f$\frac{1}{\tau_k}(w^{k}-w^{k-1}) \f$
  !> on element \f$ elem \f$  in
  !> integ nodes, i.e. recomputation of elem%w into wi in volume integ. nodes
  subroutine Eval_w_t_Elem(elem, wi)
    type(element), intent(in):: elem      ! elem = element
    real, dimension(1:elem%Qdof,1:ndim), intent(inout):: wi ! w in integ. nodes
    real, dimension(:,:), pointer:: phi ! local store arrays
    integer :: dof, Qdof, k, kst, i

    dof = elem%dof
    Qdof = elem%Qdof

    phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)

    !write(*,'(a4,2i5)') 'es:',elem%dof, dof
    do k=1,ndim              ! k = index of component of w
       kst = dof*(k-1) + 1

       do i=1,Qdof
          wi(i,k) = dot_product(elem%w(0,kst:kst+dof-1) - elem%w(1,kst:kst+dof-1), &
               phi(1:dof,i) ) /state%time%tau(1)
       enddo
    enddo

  end subroutine Eval_w_t_Elem

  !> evaluation of the state vector \f$ w^{k}-w^{k-1} \f$
  !> on element \f$ elem \f$  in
  !> integ nodes, i.e. recomputation of elem%w into wi in volume integ. nodes
  subroutine Eval_w_w_Elem(elem, wi)
    type(element), intent(in):: elem      ! elem = element
    real, dimension(1:elem%Qdof,1:ndim), intent(inout):: wi ! w in integ. nodes
    real, dimension(:,:), pointer:: phi ! local store arrays
    integer :: dof, Qdof, k, kst, i

    dof = elem%dof
    Qdof = elem%Qdof

    phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)

    !write(*,'(a4,2i5)') 'es:',elem%dof, dof
    do k=1,ndim              ! k = index of component of w
       kst = dof*(k-1) + 1

       do i=1,Qdof
          wi(i,k) = dot_product(elem%w(0,kst:kst+dof-1) - elem%w(1,kst:kst+dof-1), &
               phi(1:dof,i) )
       enddo
    enddo

  end subroutine Eval_w_w_Elem

  !> evaluation of  \f$\bar{w}=\frac{1}{|K|}\int_K w \, dx \f$ on element \f$ K = elem \f$
  subroutine Eval_aver_w_Elem(elem, w)
    type(element), intent(in):: elem      ! elem = element
    real, dimension(1:ndim), intent(inout):: w ! w in integ. nodes
    real, dimension(:), allocatable :: wi ! w in integ. nodes
    real, dimension(:,:), pointer:: phi ! local store arrays
    real, dimension(:), pointer:: weights ! local store arrays
    integer :: dof, Qdof, k, kst, i

    dof = elem%dof
    Qdof = elem%Qdof

    allocate(wi(1:Qdof) ) ! w in integ. nodes

    phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)
    weights => state%space%V_rule(elem%Qnum)%weights(1:Qdof)

    do k=1,ndim              ! k = index of component of w
       kst = dof*(k-1) + 1
       wi(1:Qdof) = matmul(elem%w(0,kst:kst+dof-1), phi(1:dof, 1:Qdof) )
       w(k) = dot_product(wi(1:Qdof), weights(1:Qdof) )
    enddo

    deallocate(wi)
  end subroutine Eval_aver_w_Elem


  !> evaluation of  \f$\bar{w}=\frac{1}{|K|}\int_K w \, dx \f$ on element \f$ K = elem \f$
  subroutine Eval_min_w_Elem(elem, w)
    type(element), intent(in):: elem      ! elem = element
    real, dimension(1:ndim), intent(inout):: w ! w in integ. nodes
    real, dimension(:), allocatable :: wi ! w in integ. nodes
    real, dimension(:,:), pointer:: phi ! local store arrays
    real, dimension(:), pointer:: weights ! local store arrays
    integer :: dof, Qdof, k, kst, i

    dof = elem%dof

    Qdof = state%space%L_rule(elem%deg+1)%Qdof
    phi => state%space%L_rule(elem%deg+1)%phi(1:dof,1:Qdof)

    !Qdof = elem%Qdof
    !phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)
    !weights => state%space%V_rule(elem%Qnum)%weights(1:Qdof)

    allocate(wi(1:Qdof) ) ! w in integ. nodes


    do k=1,ndim              ! k = index of component of w
       kst = dof*(k-1) + 1
       wi(1:Qdof) = matmul(elem%w(0,kst:kst+dof-1), phi(1:dof, 1:Qdof) )
       w(k) = minval(wi(1:Qdof))
    enddo

    deallocate(wi)
  end subroutine Eval_min_w_Elem

  !> evaluation of average of the the derivatives of the state vector \f$w\f$ on element \f$ elem \f$  in
  !> integ nodes
  subroutine Eval_aver_Dw_Elem(elem, Dw)
    type(element), intent(in) :: elem      ! elem = element
    real, dimension(1:ndim,1:nbDim), intent(inout):: Dw ! Dw in integ. nodes
    real, dimension(:,:,:), allocatable:: Dwi ! Dw in integ. nodes
    real, dimension(:), pointer:: weights ! local store arrays
    integer :: dof, Qdof, k, kst, i

    dof = elem%dof
    Qdof = elem%Qdof

    ! derivatives on real element
    allocate( Dwi(1:Qdof,1:ndim,1:nbDim) )
    call Eval_Dw_Elem(elem, Dwi)

    weights => state%space%V_rule(elem%Qnum)%weights(1:Qdof)

    do k=1,ndim              ! k = index of component of w
       Dw(k, 1:nbDim) = matmul( weights(1:Qdof), Dwi(1:Qdof,k,1:nbDim) )
    enddo

    deallocate( Dwi )
  end subroutine Eval_aver_Dw_Elem


  !> evaluation of the derivatives of the state vector \f$w\f$ on element \f$ elem \f$  in
  !> integ nodes
  subroutine Eval_Dw_Elem(elem, Dwi)
    type(element), intent(in) :: elem      ! elem = element
    real, dimension(1:elem%Qdof,1:ndim,1:nbDim), intent(inout):: Dwi ! Dw in integ. nodes
    real, dimension(:,:,:), allocatable :: Dphi
    integer :: dof, Qdof, k, kst, i

    dof = elem%dof
    Qdof = elem%Qdof

    ! derivatives on real element
    allocate(Dphi(1:dof, 1:nbDim, 1:Qdof) )

    call Eval_Dphi(elem, dof, Dphi)

    do k=1,ndim              ! k = index of component of Dw
       kst = dof*(k-1) + 1
       do i=1,Qdof
          Dwi(i, k, 1) = dot_product(elem%w(0,kst:kst+dof-1), Dphi(1:dof, 1, i) )
          Dwi(i, k, 2) = dot_product(elem%w(0,kst:kst+dof-1), Dphi(1:dof, 2, i) )
       enddo
    enddo
    deallocate(Dphi)
  end subroutine Eval_Dw_Elem


  !> evaluation of the derivatives of the state vector \f$w\f$ on element \f$ elem \f$  in
  !> integ nodes at time t_{k-1}
  subroutine Eval_Dw_Elem_time(elem, Dwi)
    type(element), intent(in) :: elem      ! elem = element
    real, dimension(1:elem%Qdof,1:ndim,1:nbDim), intent(inout):: Dwi ! Dw in integ. nodes
    real, dimension(:,:,:), allocatable :: Dphi
    integer :: dof, Qdof, k, kst, i

    dof = elem%dof
    Qdof = elem%Qdof

    ! derivatives on real element
    allocate(Dphi(1:dof, 1:nbDim, 1:Qdof) )

    call Eval_Dphi(elem, dof, Dphi)

    do k=1,ndim              ! k = index of component of Dw
       kst = dof*(k-1) + 1
       do i=1,Qdof
          Dwi(i, k, 1) = dot_product(elem%w(1,kst:kst+dof-1), Dphi(1:dof, 1, i) )
          Dwi(i, k, 2) = dot_product(elem%w(1,kst:kst+dof-1), Dphi(1:dof, 2, i) )
       enddo
    enddo
    deallocate(Dphi)
  end subroutine Eval_Dw_Elem_time


  !> evaluation of the derivatives of the state vector \f$w\f$ on element \f$ elem \f$  in
  !> integ nodes at time t, t=1 => t_k, t=0 => t_{k-1}
  subroutine Eval_Dw_Elem_at_time(elem, Dwi, t)
    type(element), intent(in) :: elem      ! elem = element
    real, dimension(1:elem%Qdof,1:ndim,1:nbDim), intent(inout):: Dwi ! Dw in integ. nodes
    real, intent(in) :: t
    real, dimension(:,:,:), allocatable :: Dphi
    integer :: dof, Qdof, k, kst, i

    dof = elem%dof
    Qdof = elem%Qdof

    ! derivatives on real element
    allocate(Dphi(1:dof, 1:nbDim, 1:Qdof) )

    call Eval_Dphi(elem, dof, Dphi)

    do k=1,ndim              ! k = index of component of Dw
       kst = dof*(k-1) + 1
       do i=1,Qdof
          Dwi(i, k, 1) = dot_product( &
               t*elem%w(0,kst:kst+dof-1) +(1-t)*elem%w(1,kst:kst+dof-1), Dphi(1:dof, 1, i) )
          Dwi(i, k, 2) = dot_product( &
               t*elem%w(0,kst:kst+dof-1) +(1-t)*elem%w(1,kst:kst+dof-1), Dphi(1:dof, 2, i) )
       enddo
    enddo
    deallocate(Dphi)
  end subroutine Eval_Dw_Elem_at_time


  !> function \f$func\f$  given on element \f$ elem \f$ by its basis coefficients,
  !> subroutines evaluate it in integ nodes
  subroutine Eval_DVec_Elem(elem, dof, w,  Dwi)
    type(element), intent(in) :: elem      ! elem = element
    integer, intent(in) :: dof
    real, dimension(1:ndim, 1:dof), intent(in):: w ! r in integ. nodes
    real, dimension(1:ndim, 1:elem%Qdof,1:nbDim), intent(inout):: Dwi !Dw in integ nodes
    real, dimension(:,:,:), allocatable :: Dphi
    integer ::  Qdof, k

    Qdof = elem%Qdof

    ! derivatives on real element
    allocate(Dphi(1:dof, 1:nbDim, 1:Qdof) )

    call Eval_Dphi(elem, dof, Dphi)

    do k=1,nbDim              ! k = index of component of Dw
       Dwi(1:ndim, 1:Qdof, k) = matmul(w(1:ndim,1:dof), Dphi(1:dof, k, 1:Qdof) )
    enddo

    deallocate(Dphi)
  end subroutine Eval_DVec_Elem


  !> evaluation of the scalar function \f$ r \f$ on the \f${ie}\f$-th
  !> edge of element \f$ {elem} \f$ in
  !> integ nodes, i.e. recomputation of elem%w into wi in edge integ. nodes
  subroutine Eval_func_Edge(elem, ie, r, wi, opposite)
    type(element), intent(in):: elem ! elem = element
    integer, intent(in) :: ie               ! index of the edge
    real, dimension(1:elem%dof), intent(in):: r !w in integ. nodes
    real, dimension(1:elem%face(fGdof,ie)), intent(inout):: wi !w in integ. nodes
    logical, intent(in) :: opposite
    real, dimension(:,:), allocatable:: phi ! test functions
    real, dimension(:,:),   pointer:: phiA ! pointers to test functions
    integer :: dof, Qdof, k, l, kst

    dof = elem%dof
    Qdof = elem%face(fGdof,ie)

    allocate(phi(1:dof, 1:Qdof))
    call Eval_Phi_Edge(elem, dof, ie, phi, opposite)

    do l=1,Qdof
       wi(l) = dot_product(r(1:dof), phi(1:dof,l) )
    enddo
    deallocate(phi)

  end subroutine Eval_func_Edge

  !> evaluation of the state vector \f$w\f$ on the \f${ie}\f$-th
  !> edge of element \f$ {elem} \f$ in
  !> integ nodes, i.e. recomputation of elem%w into wi in edge integ. nodes
  subroutine Eval_w_Edge(elem, ie, wi, opposite)
    type(element), intent(in):: elem ! elem = element
    integer, intent(in) :: ie               ! index of the edge
    real, dimension(1:elem%face(fGdof,ie),1:ndim), intent(inout):: wi !w in integ. nodes
    logical, intent(in) :: opposite
    real, dimension(:,:), allocatable:: phi ! test functions
    real, dimension(:,:),   pointer:: phiA ! pointers to test functions
    integer :: dof, Qdof, k, l, kst

    dof = elem%dof
    Qdof = elem%face(fGdof,ie)

    allocate(phi(1:dof, 1:Qdof))
    call Eval_Phi_Edge(elem, dof, ie, phi, opposite)

    do k=1,ndim              ! k = index of component of w
       kst = dof*(k-1) + 1

       do l=1,Qdof
          wi(l,k) = dot_product(elem%w(0,kst:kst+dof-1), phi(1:dof,l) )
       enddo
    enddo
    deallocate(phi)
  end subroutine Eval_w_Edge


  !> evaluation of the state vector \f$w\f$ and \f$\Pi^{p-1} w\f$  on the \f${ie}\f$-th
  !> edge of element \f$ {elem} \f$ in
  !> integ nodes, i.e. recomputation of elem%w into wi in edge integ. nodes
  subroutine Eval_w_EdgeProj(elem, ie, wi, wi1, opposite)
    type(element), intent(in):: elem ! elem = element
    integer, intent(in) :: ie               ! index of the edge
    real, dimension(1:elem%face(fGdof,ie),1:ndim), intent(inout):: wi, wi1 !w in integ. nodes
    logical, intent(in) :: opposite
    real, dimension(:,:), allocatable:: phi ! test functions
    real, dimension(:,:),   pointer:: phiA ! pointers to test functions
    integer :: dof, dof1, Qdof, k, l, kst

    dof = elem%dof
    dof1 = dof - elem%deg - 1   !  DOF of P-1 projection
    if(elem%deg == 0) dof1 = dof

    Qdof = elem%face(fGdof,ie)

    allocate(phi(1:dof, 1:Qdof))
    call Eval_Phi_Edge(elem, dof, ie, phi, opposite)

    do k=1,ndim              ! k = index of component of w
       kst = dof*(k-1) + 1

       do l=1,Qdof
          wi(l,k) = dot_product(elem%w(0,kst:kst+dof-1), phi(1:dof,l) )
          wi1(l,k) = dot_product(elem%w(0,kst:kst+dof1-1), phi(1:dof1,l) )
       enddo
    enddo
    !if(elem%i == 100) print*,'########',1,dof, dof1

    deallocate(phi)
  end subroutine Eval_w_EdgeProj

  !> evaluation of the state vector \f$w\f$ on the \f${ie}\f$-th
  !> edge of element \f$ {elem} \f$ in
  !> integ nodes, i.e. recomputation of elem%w into wi in edge integ. nodes at t_{k-1}
  subroutine Eval_w_Edge_time(elem, ie, wi, opposite)
    type(element), intent(in):: elem ! elem = element
    integer, intent(in) :: ie               ! index of the edge
    real, dimension(1:elem%face(fGdof,ie),1:ndim), intent(inout):: wi !w in integ. nodes
    logical, intent(in) :: opposite
    real, dimension(:,:), allocatable:: phi ! test functions
    real, dimension(:,:),   pointer:: phiA ! pointers to test functions
    integer :: dof, Qdof, k, l, kst

    dof = elem%dof
    Qdof = elem%face(fGdof,ie)

    allocate(phi(1:dof, 1:Qdof))
    call Eval_Phi_Edge(elem, dof, ie, phi, opposite)

    do k=1,ndim              ! k = index of component of w
       kst = dof*(k-1) + 1

       do l=1,Qdof
          wi(l,k) = dot_product(elem%w(1,kst:kst+dof-1), phi(1:dof,l) )
       enddo
    enddo
    deallocate(phi)
  end subroutine Eval_w_Edge_time

  !> evaluation of the average of state vector \f$w\f$ on the \f$ {ie} \f$ -th
  !> inner edge of element \f$ elem \f$  and its neighbour \f$ elem1\f$  in
  !> integ nodes, i.e. recomputation of elem%w into state%wi in edge integ. nodes
  subroutine Eval_aver_w_Edge(elem, elem1, ie, Qdof, wi)
    use  pedes_averaging
    type(element), intent(in):: elem, elem1 ! elem = element
    integer, intent(in) :: ie               ! index of the edge
    integer, intent(in) :: Qdof             ! number of integ. nodes
    real, dimension(1:Qdof, 1:ndim), intent(inout) :: wi ! state vector in integ. nodes
!    real, dimension(1:Qdof, 1:ndim)  :: wL, wR
    real, dimension(:,:), allocatable :: phi, phi1 ! test functions
    real, dimension(:,:), pointer :: phiA, phiA1 ! pointers to test functions
    real :: wik
    integer :: dof, dof1, Qnum, k, kst
    integer :: l1, kst1, l, ie1


    dof = elem%dof
    dof1 = elem1%dof

    ie1 = elem%face(nei_i,ie)

    !! seting of degree of the Gauss quadrature
    Qnum = elem%face(fGnum,ie)

    if(Qdof .ne. state%space%G_rule(Qnum)%Qdof) then
       print*,'Uncompatible Qdof in Eval_aver_w'
       stop
    endif

    allocate(phi(1:dof, 1:Qdof))
    call Eval_Phi_Edge(elem, dof, ie, phi, .false.)

    allocate(phi1(1:dof1, 1:Qdof))
    call Eval_Phi_Edge(elem1, dof1, ie1, phi1, .true.)

    !!!phiA  => state%space%G_rule(Qnum)%phi(elem%type,  ie,  1, 1:dof, 1:Qdof)
    !!!phiA1 => state%space%G_rule(Qnum)%phi(elem1%type, ie1, 1, 1:dof1, 1:Qdof)

    ! evaluation of (w_ie^+ + w_ie^-)/2 in integ. nodes
    do l=1, Qdof
       !!!l1 = Qdof - l + 1
       do k=1, ndim
          kst = (k-1)*dof + 1
          kst1 = (k-1)*dof1 + 1
          wi(l,k) = (dot_product(phi(1:dof ,l),  elem%w(0,kst: kst+dof-1) ) &
               + dot_product(phi1(1:dof1 ,l), elem1%w(0,kst1: kst1+dof1-1) ) ) /2.

          !wL(l,k) = dot_product(phi(1:dof ,l),  elem%w(0,kst: kst+dof-1) )
          !wR(l,k) = dot_product(phi1(1:dof1 ,l), elem1%w(0,kst1: kst1+dof1-1) )

          !!!wik = (dot_product(phiA(1:dof ,l),  elem%w(0,kst: kst+dof-1) ) &
          !!!     + dot_product(phiA1(1:dof1 ,l1), elem1%w(0,kst1: kst1+dof1-1) ) ) /2.
          !!!if(wi(l,k) - wik .ne. 0.) &
          !!!     write(*,'(a3,4i5,3es14.6)') '?!?',elem%i,elem1%i,l,k,wi(l,k),wik, wi(l,k) - wik
       enddo

       !if(state%print)then
       !if(elem%i ==1 .and. elem1%i==2 .and. ie ==1)then
       !   write(*,'(a5,5i5)')'test',elem%i, elem1%i, ie,ie1
       !   write(*,'(a5,4e12.4)')'wi',wi(l,1:4)
       !   write(*,'(a5,4e12.4)')'wL',wL(l,1:4)
       !   write(*,'(a5,4e12.4)')'wR',wR(l,1:4)
       !   write(*,'(a2,i2,40es12.4)') '~~',l,elem1%w(0,1:7)
       !endif

      ! if ( abs( wi(l,1) - 1.0) > 0.01 ) print*, 'wi in Eval_aver_w_Edge: ' , wi(l,1)

!       if ( abs( wi(l,1) - 1.0) > 0.1 ) then
!         print*, 'elem%i , Qdof, l , Qnum', elem%i, Qdof, l, Qnum
!         print*, 'wi in Eval_aver_w_Edge: ' , wi(l,1)
!         stop
!       endif

       if(wi(l,1) < 0. .and. ndim == 4) then
          open(44, file='bad_node', status = 'UNKNOWN')
          write(44,*) elem%xc(:)
          write(44,*) (elem%xc(:) + elem1%xc(:))/2
          write(44,*) elem1%xc(:)
          close(44)

          ! write(*,'(a5,5i5,4e12.4)')'INNER',elem%i, elem1%i, ie,ie1,l,wi(l,1:4)
          ! write(*,*) 'tau = ',state%time%tau(1), state%max_eigenvals*state%time%tau(1)
          ! write(*,*) 'xc = ',elem%xc(:)
          ! write(*,'(a4,i3,40es12.4)') 'deg=',elem%deg,elem%w(0,:)
          ! write(*,'(a4,i3,40es12.4)') 'deg=',elem1%deg,elem1%w(0,:)

          ! print*,'Elem :................................'
          ! write(*,'(a6,20i5)') 'elem:',elem%i,elem%flen, ie
          ! write(*,'(a6,20i5)') 'face:',elem%face(1,:)
          ! write(*,'(a6,20i5)') 'neig:',elem%face(2,:)
          ! write(*,'(a6,20i5)') 'ne_i:',elem%face(3,:)
          ! if(elem%HGnode) then
          !    write(*,'(a6,20i5)') 'HGfa:',elem%HGface(1,:)
          !    write(*,'(a6,20i5)') 'HGfa:',elem%HGface(2,:)
          !    write(*,'(a6,20i5)') 'HGver:',elem%HGvertex(:)
          ! endif
          ! print*,'Elem1:'
          ! write(*,'(a6,20i5)') 'elem1:',elem1%i,elem%flen, ie
          ! write(*,'(a6,20i5)') 'face:',elem1%face(1,:)
          ! write(*,'(a6,20i5)') 'neig:',elem1%face(2,:)
          ! write(*,'(a6,20i5)') 'ne_i:',elem1%face(3,:)
          ! if(elem1%HGnode) then
          !    write(*,'(a6,20i5)') 'HGfa:',elem1%HGface(1,:)
          !    write(*,'(a6,20i5)') 'HGfa:',elem1%HGface(2,:)
          !    write(*,'(a6,20i5)') 'HGver:',elem1%HGvertex(:)
          ! endif
          ! print*,'................................'

          ! do l1=1,dof
          !    write(*,'(1i5,10es12.4)') l1,phi(l1,:)
          ! enddo
          ! do l1=1,dof1
          !    write(*,'(1i5,10es12.4)') l1,phi1(l1,:)
          ! enddo
          ! print*,'--------------------------'
          ! stop

          write(*,'(a20,2(a4,es10.2),a8,2i5,a9,i2, a17)') &
               'nonpositive rho: ', &
               ' e =', wi(l,4), &
               ' r =', wi(l,1), &
               ', elems=',elem%i,elem1%i,  &
               ', Newton= ', state%nlSolver%iter , &
               ', file "bad_node"'

          wi(l,1) = 0.001

       endif
    enddo !! l

    ! in the case of the limiting, we use the averaged values
    if(state%modelName == 'pedes' ) then
       if( limit_w_bar(elem%i) == 1 .or. limit_w_bar(elem1%i) == 1) then
          !write(*,'(a5, 30es14.6)') 'orig:', wi(1:Qdof,2) / wi(1:Qdof,1)
          !print*,elem%xc, elem%i
          !print*,elem1%xc, elem1%i
          !stop'deud437d'
          ! first element
          if( limit_w_bar(elem%i) == 1) then
             wi(1:Qdof,1) = w_bar(elem%i, 1)
             wi(1:Qdof,2) = w_bar(elem%i, 2)
             wi(1:Qdof,3) = w_bar(elem%i, 3)
          else
             do k=1, ndim
                kst = (k-1)*dof + 1
                !kst1 = (k-1)*dof1 + 1
                wi(1:Qdof,k) = matmul( elem%w(0,kst: kst+dof-1), phi(1:dof ,1:Qdof) )
             enddo
          endif

          ! second element
          if( limit_w_bar(elem1%i) == 1) then
             wi(1:Qdof,1) =  wi(1:Qdof,1) + w_bar(elem1%i, 1)
             wi(1:Qdof,2) =  wi(1:Qdof,2) + w_bar(elem1%i, 2)
             wi(1:Qdof,3) =  wi(1:Qdof,3) + w_bar(elem1%i, 3)
          else
             do k=1, ndim
                !kst = (k-1)*dof + 1
                kst1 = (k-1)*dof1 + 1
                wi(1:Qdof,k) = wi(1:Qdof,k)  &
                     + matmul( elem1%w(0,kst1: kst1+dof1-1), phi1(1:dof1 ,1:Qdof) )
             enddo
          endif
          wi(1:Qdof,1:3)  = wi(1:Qdof,1:3) / 2

          !write(*,'(a5, 30es14.6)') 'new:',  wi(1:Qdof,2) / wi(1:Qdof,1)
          !print*,'--------------------------------'
          !stop'dejd38dyu93iwks'
       endif
    endif

    deallocate(phi, phi1)

  end subroutine Eval_aver_w_Edge


  !> evaluation of the derivatives of the state vector \f$w\f$ on the \f$ {ie}\f$ -th
  !> edge of element \f$ elem \f$ in
  !> integ nodes, i.e. recomputation of elem%w into wi in edge integ. nodes
  !> values are multiplied by the Gauss weights
  subroutine Eval_Dw_Edge(elem, ie, Dwi, opposite)
    type(element), intent(in):: elem ! elem = element
    integer, intent(in) :: ie               ! index of the edge
    real, dimension(1:elem%face(fGdof,ie),1:ndim,1:nbDim), intent(inout):: Dwi !w in integ. nodes
    logical, intent(in) :: opposite
    real, dimension(:,:,:), allocatable :: Dphi
    integer :: dof, Qdof, Qnum, k, l, kst

    dof = elem%dof
    Qdof = elem%face(fGdof,ie)
    Qnum = elem%face(fGnum,ie)

    ! derivatives on real element
    allocate(Dphi(1:dof, 1:nbDim, 1:Qdof) )

    call Eval_Dphi_Edge(elem, dof, ie, Dphi, opposite)


    do k=1,ndim              ! k = index of component of w
       kst = dof*(k-1) + 1

       ! Dphi are premultiplied by the weights !!!!!
       do l=1,Qdof
          Dwi(l,k,1) = dot_product(elem%w(0,kst:kst+dof-1), Dphi(1:dof, 1, l) ) &
               / state%space%G_rule(Qnum)%weights(l)
          Dwi(l,k,2) = dot_product(elem%w(0,kst:kst+dof-1), Dphi(1:dof, 2, l) ) &
               / state%space%G_rule(Qnum)%weights(l)
       enddo
    enddo

    deallocate(Dphi)
  end subroutine Eval_Dw_Edge

  !> evaluate the inviscid fluxes f_s
  !> in integ nodes on  elem,
  subroutine Eval_f_s_Elem(Set_f_s, elem,  f_s)
    interface
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x, ie )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
         integer, intent(in) :: ie
      end subroutine Set_f_s
    end interface
    type(element), intent (in) :: elem   ! element
    real, dimension(1:elem%Qdof, 1:nbDim, 1:ndim), intent(out) :: f_s ! output
    real, dimension(:,:), allocatable :: wi
    integer :: Qdof

    Qdof = elem%Qdof
    allocate( wi(1:Qdof,1:ndim) )

    call Eval_w_Elem(elem,  wi)
    call Set_f_s(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), f_s(1:Qdof, 1:nbDim, 1:ndim), &
         elem%xi(0,1:Qdof, 1:nbDim), elem%i)

    deallocate(wi )
  end subroutine Eval_f_s_Elem

  !> evaluate the inviscid fluxes f_s
  !> in integ nodes on  elem t, t=1 => t_k, t=0 => t_{k-1}
  subroutine Eval_f_s_Elem_at_time(Set_f_s, elem,  f_s, t)
    interface
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x, ie )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
         integer, intent(in) :: ie
      end subroutine Set_f_s
    end interface
    type(element), intent (in) :: elem   ! element
    real, dimension(1:elem%Qdof, 1:nbDim, 1:ndim), intent(out) :: f_s ! output
    real, intent(in) :: t
    real, dimension(:,:), allocatable :: wi
    integer :: Qdof

    Qdof = elem%Qdof
    allocate( wi(1:Qdof,1:ndim) )

    call Eval_w_Elem_at_time(elem,  wi, t)
    call Set_f_s(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), f_s(1:Qdof, 1:nbDim, 1:ndim), &
         elem%xi(0,1:Qdof, 1:nbDim), elem%i)

    deallocate(wi )
  end subroutine Eval_f_s_Elem_at_time


  !> evaluate numerical flux
  !> \f$ H(u, v, n) \f$
  !> in integ nodes on the ie-th edge of elem
  subroutine Eval_NFlux_Edge(Set_Ppm, elem, ie, Rflux)
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
    end interface
    type(element), intent (inout) :: elem   ! element
    integer, intent(in) :: ie            ! index of the edge
    real, dimension(1:elem%face(fGdof,ie), 1:ndim), intent(out) :: Rflux ! output
    class(element), pointer :: elem1
    real, dimension(:,:), allocatable :: wi, wii, nc
    real, dimension(:,:, :, :), allocatable :: Ppm
    real, dimension(1:nbDim) :: xc
    integer :: l, Qdof, ii, k, ie1

    Qdof = elem%face(fGdof,ie)

    allocate(wi(1:Qdof,1:ndim), wii(1:Qdof,1:ndim) )

    ii  = elem%face(neigh, ie)
    ie1 = elem%face(nei_i, ie)

    if(ii > 0) then
       elem1 => grid%elem(ii)
       call Eval_aver_w_Edge(elem, elem1, ie, Qdof, wi(1:Qdof,1:ndim))
       xc(:) = (elem%xc(1:nbDim)+elem1%xc(1:nbDim))/2

    else
       call Eval_w_Edge(elem, ie, wi(1:Qdof,1:ndim), .false.)
       xc(:) = elem%xc(1:nbDim)
    endif

    allocate(nc(1:Qdof, 1:nbDim) )
    if(elem%ibcur > 0) then
       !courved edge
       nc(1:Qdof,1:nbDim) = elem%nc(1:Qdof,1:nbDim)
    else
     ! straight edge
       nc(1:Qdof,1) = elem%n(ie,1)
       nc(1:Qdof,2) = elem%n(ie,2)
    endif


    allocate(Ppm(1:Qdof, 1:2, 1:ndim, 1:ndim) )
    ! Vijajasyndaram
    call Set_Ppm(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), nc(1:Qdof,1:nbDim), &
         xc(1:nbDim), Ppm(1:Qdof,1:2, 1:ndim, 1:ndim), 1./elem%area, elem )


    call Eval_w_Edge(elem,  ie,  wi,  .false.)

    if(ii > 0) then
       call Eval_w_Edge(elem1, ie1, wii, .true.)
    else

       do k=1,Qdof
          xc(1:nbDim) = grid%b_edge(-elem%face(neigh,ie))%x_div(k, 1:nbDim)
          call Exact_Scalar(xc(1:nbDim), wii(k,1:ndim), state%time%ctime )

          !if(elem%i <= 3) &
          ! write(*,'(a5,2i5, 13e12.4,a5)') 'ri:', elem%i, ie, xc(1:nbDim), wi(k,1), wii(k,1:1)
       enddo

       !!wii(:,:) = wi(:,:)
    endif

    do l=1,Qdof
       Rflux(l,1:ndim) = matmul(Ppm(l,1,1:ndim, 1:ndim), wi(l, 1:ndim) ) &
            + matmul(Ppm(l,2,1:ndim, 1:ndim), wii(l, 1:ndim) )
    enddo

    deallocate(wi, wii, nc, Ppm )

  end subroutine Eval_NFlux_Edge


  !> evaluate the viscous fluxes R_s
  !> \f$ \sum_{s=1}^d R_s(w,\nabla w) n_s  \f$
  !> in integ nodes on  elem,
  subroutine Eval_R_s_Elem(Set_R_s, elem,  R_s)
    interface
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

    type(element), intent (in) :: elem   ! element
    real, dimension(1:elem%Qdof, 1:nbDim, 1:ndim), intent(out) :: R_s ! output
    real, dimension(:,:), allocatable :: wi
    real, dimension(:,:,:), allocatable :: Dwi
    real, dimension(:,:), allocatable :: Re_1
    integer :: Qdof, i, j


    Qdof = elem%Qdof

    allocate( wi(1:Qdof,1:ndim),  Dwi(1:Qdof, 1:ndim, 1:nbDim), Re_1(1:iRe, 1:Qdof)  )
    Re_1(1:iRe, 1:Qdof) = state%model%Re1

    !print*, 'ctime in Eval_R_s_Elem ', state%time%ctime

    call Eval_w_Elem(elem,  wi)
    call Eval_Dw_Elem(elem, Dwi(1:Qdof,1:ndim, 1:nbDim) )

!    print*, 'WI:' , wi(1,1:ndim)
!    print*, 'DWI:', Dwi(1,1:ndim, 1:nbDim)

    call Set_R_s(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof, 1:ndim, 1:nbDim), &
         Re_1(1:iRe, 1:Qdof), R_s(1:Qdof, 1:nbDim, 1:ndim), elem%xi(0,1:Qdof, 1:2))

    deallocate(wi, Dwi, Re_1 )

  end subroutine Eval_R_s_Elem


  !> evaluate the viscous fluxes R_s
  !> \f$ \sum_{s=1}^d R_s(w,\nabla w) n_s  \f$
  !> in integ nodes on elem, at time t, t=1 => t_k, t=0 => t_{k-1}
  subroutine Eval_R_s_Elem_at_time(Set_R_s, elem,  R_s, t )
    interface
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
    type(element), intent (in) :: elem   ! element
    real, intent(in) :: t
    real, dimension(1:elem%Qdof, 1:nbDim, 1:ndim), intent(out) :: R_s ! output
    real, dimension(:,:), allocatable :: wi
    real, dimension(:,:,:), allocatable :: Dwi
    real, dimension(:,:), allocatable :: Re_1
    integer :: Qdof, i, j

    Qdof = elem%Qdof

    allocate( wi(1:Qdof,1:ndim),  Dwi(1:Qdof, 1:ndim, 1:nbDim), Re_1(1:iRe, 1:Qdof)  )
    Re_1(1:iRe, 1:Qdof) = state%model%Re1

    call Eval_w_Elem_at_time(elem,  wi, t)
    call Eval_Dw_Elem_at_time(elem, Dwi(1:Qdof,1:ndim, 1:nbDim), t )

    call Set_R_s(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof, 1:ndim, 1:nbDim), &
         Re_1(1:iRe, 1:Qdof), R_s(1:Qdof, 1:nbDim, 1:ndim), elem%xi(0,1:Qdof, 1:2) )

    deallocate(wi, Dwi, Re_1 )

  end subroutine Eval_R_s_Elem_at_time

  !> evaluate the viscous fluxes R_s
  !> \f$ \sum_{s=1}^d R_s(w,\nabla w) n_s  \f$
  !> in integ nodes on the ie-th edge of elem,
  !> is opposite = .true., the orientation \f$ (n_1,\dots,n_d)\f$ is opposite
  subroutine Eval_R_s_Edge(Set_R_s, elem, ie, Rflux, opposite)
    interface
      subroutine Set_R_s(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi )
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
          real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
      end subroutine Set_R_s
    end interface
    type(element), intent (in) :: elem   ! element
    integer, intent(in) :: ie            ! index of the edge
    real, dimension(1:elem%face(fGdof,ie), 1:ndim), intent(out) :: Rflux ! output
    logical, intent(in) :: opposite
    real, dimension(:,:), allocatable :: wi
    real, dimension(:,:,:), allocatable :: Dwi
    real, dimension(:,:,:), allocatable :: R_s
    real, dimension(:,:), allocatable :: Re_1
    integer :: Qdof, i, j

    Qdof = elem%face(fGdof,ie)

!    stop 'Problem: Re_1 should be 2 dimensional array and integer iRe is not defined here'
    ! iRe is defined in paramets
    allocate( wi(1:Qdof,1:ndim),  Dwi(1:Qdof,1:ndim,1:nbDim), Re_1(1:iRe, 1:Qdof)  )
    allocate( R_s(1:Qdof, 1:nbDim, 1:ndim)  )
    Re_1(:,:) = state%model%Re1

    call Eval_w_Edge(elem, ie,  wi, opposite)
    call Eval_Dw_Edge(elem,  ie,  Dwi(1:Qdof,1:ndim,1:nbDim), opposite)

    ! FR added Re_1(1:iRe,1:Qdof), is it correct?
    call Set_R_s(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof, 1:ndim, 1:nbDim), &
         Re_1(1:iRe,1:Qdof), R_s(1:Qdof, 1:nbDim, 1:ndim), elem%xi(ie,1:Qdof, 1:2)  )

!    call Set_R_s(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof, 1:ndim, 1:nbDim), &
!         Re_1(1:Qdof), R_s(1:Qdof, 1:nbDim, 1:ndim), elem%xi(ie,1:Qdof, 1:2)  )

    Rflux(1:Qdof, 1:ndim) = 0.

    !print*, 'R_S = ', R_s(1, 1:nbDim, 1), 'normal:', elem%n(ie, :)

    if(elem%face(neigh,ie) > 0 .or. elem%F%iFlin) then ! inner or straight boundary  edge
       do i=1,nbDim
          Rflux(1:Qdof, 1:ndim) = Rflux(1:Qdof, 1:ndim) &
               + R_s(1:Qdof, i, 1:ndim)* elem%n(ie, i)
       enddo
    else
       do i=1,nbDim
          do j=1,Qdof
             Rflux(j, 1:ndim) = Rflux(j, 1:ndim) +  R_s(j, i, 1:ndim) * elem%nc(j, i)
          enddo
       enddo
    endif

    ! for opposite the orientation of the normal is oposite
    if(opposite) Rflux(1:Qdof, 1:ndim) = - Rflux(1:Qdof, 1:ndim)

    deallocate(wi, Dwi, Re_1, R_s )
  end subroutine Eval_R_s_Edge

  !> evaluate average of the viscous fluxes
  !> \f$ \langle\sum_{s=1}^d R_s(w,\nabla w) n_s \rangle \f$
  !> in integ nodes on the ie-th edge of elem
  subroutine Eval_aver_R_s_Edge(Set_R_s, elem, ie, Rflux)
    interface
      subroutine Set_R_s(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1     ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
      end subroutine Set_R_s
    end interface
    type(element), intent (in) :: elem   ! element
    integer, intent(in) :: ie            ! index of the edge
    real, dimension(1:elem%face(fGdof,ie), 1:ndim), intent(out) :: Rflux ! output
    class(element), pointer :: elem1
    real, dimension(:,:, :), allocatable :: flux_LR
    integer :: Qdof, k, ie1

    Qdof = elem%face(fGdof,ie)
    allocate(flux_LR(1:2, 1:Qdof,  1:ndim) )

    ! fluxes from the left and right
    call Eval_R_s_Edge(Set_R_s, elem, ie, flux_LR(1,1:Qdof,1:ndim), .false.)

    k = elem%face(neigh, ie)
    if(k > 0 ) then ! inner edge
       elem1 => grid%elem(k)
       ie1 = elem%face(nei_i,ie)

       call Eval_R_s_Edge(Set_R_s, elem1, ie1, flux_LR(2,1:Qdof,1:ndim), .true.)
       ! their average
       Rflux(1:Qdof, 1:ndim) = &
            ( flux_LR(1, 1:Qdof, 1:ndim) + flux_LR(2, 1:Qdof, 1:ndim) )/2.

    else ! boundary edge
       Rflux(1:Qdof, 1:ndim) =  flux_LR(1, 1:Qdof, 1:ndim)
    endif

    if (elem%i == 3) then
        !write(*,'(a6,i5)') 'Qdof:', Qdof
        !write(*,'(a30,4es14.6)') 'flux_LR(1,1:Qdof,1:ndim):', flux_LR(1,1:Qdof,1)
        !write(*,'(a30,4es14.6)') 'flux_LR(2,1:Qdof,1:ndim):', flux_LR(2,1:Qdof,1)
    endif

    deallocate(flux_LR )

  end subroutine Eval_aver_R_s_Edge

  !> evaluation of the coefficients of the Lagrangian reconstruction for \f$ n\f$-BDF,
  !> from nodes \f$ t_{k-i},\ i=0,\dots, n \f$ at node t_rel\f$\in (0,1)\f$,
  !> t_rel corresponds to \f$ t^* = t_{k-1} + t_{rel}(t_k - t_{k-1})\f$
  !> i.e. \f$ w(t_{\rm rel}) = \sum_{i=0}^{Tdeg} w(t_{k-i}) Lagcoef_i \f$
  subroutine SetLagrCoeffs(Tdeg, t_rel, Lag_coef )
    integer, intent(in) :: Tdeg  ! degree of the Lagr. inter. pol.
    real, intent(in) :: t_rel     ! relative position of interp. node within the last interval
    real, dimension(0:1, 0:Tdeg), intent(out) :: Lag_coef ! coeff of Lagran. inter and its deriv
    real, dimension(:), allocatable :: t
    !real, dimension(:), pointer :: taus
    real :: tp
    integer :: i,j,k

    associate ( taus => state%time%tau(1:Tdeg) )

       ! setting of times
       allocate( t(0:Tdeg) )
       t(Tdeg) = taus(1)
       t(Tdeg-1) = 0.
       do i=Tdeg-2,0,-1
          j = Tdeg - i
          t(i) = t(i+1) - taus(j)
       enddo

       ! node, where we interpolate
       tp = t_rel * taus(1)

       do i=0, Tdeg
          ! setting of the coefficients of the Lagr. interpol. function
          Lag_coef(0, i) = product((t(:) - tp)/(t(:) - t(i)), mask=t/=t(i) )

          !its derivative
          Lag_coef(1, i) = 0.
          do j = 0, Tdeg
             if(i/= j) Lag_coef(1, i) = Lag_coef(1, i) &
                  + product((t(:) - tp)/(t(:) - t(i)), mask=(t/=t(i) .and. t/=t(j) ) )/(t(i)-t(j) )
          enddo
       enddo

       !write(*,'(a6,1i5, 4es12.4 )') 'PW',Tdeg, T_rel, tp
       !write(*,'(a6,6es12.4)') 'taus',state%time%tau(:)
       !write(*,'(a6,6es12.4)') 't_k',t(:)
       !write(*,'(a6,6es12.4)') 'coefs ',Lag_coef(0,:), sum(Lag_coef(0,:) )
       !write(*,'(a6,6es12.4)') 'Dcoefs',Lag_coef(1,:), sum(Lag_coef(1,:) )
       !write(*,*)'---------------------------------------'

       deallocate( t)
    end associate

  end subroutine SetLagrCoeffs

  ! plot the actual solution and its derivatives in integ nodes
  subroutine PlotElemSolution(ifile, elem)
    integer, intent(in) :: ifile        ! number of the file
    type(element), intent (in) :: elem     ! element where the basis is given
    real, dimension(:,:), allocatable :: wi, xi
    real, dimension(:,:,:), allocatable :: Dwi
    integer :: Qnum, Qdof, j

    Qdof = elem%Qdof
    Qnum = elem%Qnum

    allocate( xi(1:Qdof,1:nbDim), wi(1:Qdof,1:ndim), Dwi(1:Qdof,1:ndim,1:nbDim) )

    call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim), xi(1:Qdof,1:nbDim) )

    call Eval_w_Elem(elem, wi)

    call Eval_Dw_Elem(elem, Dwi)

    do j=1,Qdof
       write(ifile, *) xi(j, 1:nbDim), wi(j, 1:ndim), Dwi(j, 1:ndim, 1), Dwi(j, 1:ndim, 2)
    enddo
    write(ifile,'(x)')

    deallocate(xi, wi, Dwi)
  end subroutine PlotElemSolution



  ! 3D plot the actual solution
  subroutine PlotElemSolution3D(ifile, elem)
    integer, intent(in) :: ifile        ! number of the file
    type(element), intent (in) :: elem     ! element where the basis is given
    real, dimension(:,:), allocatable :: wi, xi, Fxi
    real, dimension(:,:,:), allocatable :: Dwi
    real, dimension(:,:), allocatable :: phiT
    type(Lagrang_rule), pointer :: L_rule
    integer :: Qnum, Qdof, Qdeg, dof2, j, j1,j2, k, nface, it

    nface = elem%deg  !3, 5

    L_rule => state%space%L_rule(nface)

    Qdeg = L_rule%Qdeg
    Qdof = L_rule%Qdof

    allocate(  xi(1:Qdof, 1:2), Fxi(1:Qdof, 1:2), wi(1:Qdof,1:ndim), Dwi(1:Qdof,1:ndim,1:nbDim) )

    xi(1:Qdof, 1:2) = L_rule%lambda(1:Qdof, 1:2)

    call ComputeF(elem, Qdof, xi(1:Qdof,1:nbDim), Fxi(1:Qdof,1:nbDim) )

    dof2 = elem%dof
    ! DG test functions in xi
    allocate( phiT(1:dof2, 1:Qdof) )

    call PHI_orthonormal(Qdof, nbDim, xi(1:Qdof, 1:nbDim), 3, dof2, phiT(1:dof2, 1:Qdof) )

    do j=1,Qdof
       do k=1, ndim
          wi(j, k) = dot_product(phiT(1:dof2, j), elem%w(0,(k-1)*dof2 + 1 : k*dof2) )
          !write(*,'(a8,40es12.4)') '???LK',wi(j,k), phiT(1:dof2, j)
       enddo
    enddo


    ! first orientation
    it= 0
    do j=0, Qdeg
       do k=0, Qdeg - j
          it = it + 1
          write(ifile, *) Fxi(it, 1:nbDim), wi(it, 1:ndim) !, Dwi(it, 1:ndim, 1), Dwi(it, 1:ndim, 2)
       enddo
       write(ifile,'(x)')
    enddo
    write(ifile,'(x)')

   ! second orientation
    do j=0, L_rule%Qdeg
       it = nface + 1 - j
       do k=0, L_rule%Qdeg - j
          write(ifile, *) Fxi(it, 1:nbDim), wi(it, 1:ndim) !, Dwi(it, 1:ndim, 1), Dwi(it, 1:ndim, 2)
          it = it + nface - k
       enddo
       write(ifile,'(x)')
    enddo
    write(ifile,'(x)')

    ! third orientation
    do j=0, L_rule%Qdeg
       it = j + 1
       do k=0, L_rule%Qdeg - j
          write(ifile, *) Fxi(it, 1:nbDim), wi(it, 1:ndim) !, Dwi(it, 1:ndim, 1), Dwi(it, 1:ndim, 2)
          it = it + nface +1 - k
       enddo
       write(ifile,'(x)')
    enddo
    write(ifile,'(x)')


    deallocate(xi, Fxi, wi, Dwi)


  end subroutine PlotElemSolution3D



  !> 3D plot the actual solution
  subroutine PlotElemFunction3D(ifile, elem,  dof2, w)
    integer, intent(in) :: ifile        ! number of the file
    type(element), intent (in) :: elem     ! element where the basis is given
    integer, intent(in) :: dof2
    real, dimension(1:dof2), intent(in) :: w
    real, dimension(:,:), allocatable :: wi, xi, Fxi
    real, dimension(:,:,:), allocatable :: Dwi
    real, dimension(:,:,:), allocatable :: DphiT, Dphi
    type(Lagrang_rule), pointer :: L_rule
    integer :: Qnum, Qdof, Qdeg,  j, j1,j2, k, nface, it


    if (dof2 == 1) then
      nface = 1
    else
       nface = elem%deg + 3 !2  !4 !5
    endif


    L_rule => state%space%L_rule(nface)

    Qdeg = L_rule%Qdeg
    Qdof = L_rule%Qdof

    allocate(  xi(1:Qdof, 1:2), Fxi(1:Qdof, 1:2), wi(1:Qdof,1:ndim), Dwi(1:Qdof,1:ndim,0:nbDim) )

    xi(1:Qdof, 1:2) = L_rule%lambda(1:Qdof, 1:2)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! xi(1:Qdof, 1) = 4* xi(1:Qdof, 1) - 0.5 !
    ! xi(1:Qdof, 2) = 4* xi(1:Qdof, 2) - 0.5 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call ComputeF(elem, Qdof, xi(1:Qdof,1:nbDim), Fxi(1:Qdof,1:nbDim) )

    ! DG test functions in xi
    allocate( DphiT(1:dof2, 0:2, 1:Qdof) )  ! second index= 0 => phi, 1,2 => Dphi
    allocate( Dphi(1:dof2, 0:2, 1:Qdof) )  ! second index= 0 => phi, 1,2 => Dphi

    ! Dphi on reference element
    call PHI_orthonormal(Qdof, nbDim, xi(1:Qdof, 1:nbDim), 3, dof2, DphiT(1:dof2, 0, 1:Qdof), &
         Dphi(1:dof2, 1:2, 1:Qdof) )

    ! Dphi on real element
    do j=1, dof2
       DphiT(j, 1, 1:Qdof) = &
            elem%F%D1F0(1,1) * Dphi(j, 1, 1:Qdof) &
            +elem%F%D1F0(1,2)* Dphi(j, 2, 1:Qdof)

       DphiT(j, 2, 1:Qdof) = &
            elem%F%D1F0(2,1) * Dphi(j, 1, 1:Qdof) &
            +elem%F%D1F0(2,2)* Dphi(j, 2, 1:Qdof)
    enddo

    do j=1,Qdof
       do k=1, 1 !!ndim
          wi(j, k) = dot_product(DphiT(1:dof2, 0, j), w(1 : dof2) )
          Dwi(j, k, 1) = dot_product(DphiT(1:dof2, 1, j), w(1 : dof2) )
          Dwi(j, k, 2) = dot_product(DphiT(1:dof2, 2, j), w(1 : dof2) )
          Dwi(j, k, 0) = sqrt( Dwi(j, k, 1)**2 +  Dwi(j, k, 2)**2)
       enddo
    enddo


    ! first orientation
    it= 0
    do j=0, Qdeg
       do k=0, Qdeg - j
          it = it + 1
          write(ifile, *) Fxi(it, 1:nbDim), wi(it, 1:ndim), &
               Dwi(it, 1:ndim, 0), Dwi(it, 1:ndim, 1), Dwi(it, 1:ndim, 2)
       enddo
       write(ifile,'(x)')
    enddo
    write(ifile,'(x)')

   ! second orientation
    do j=0, L_rule%Qdeg
       it = nface + 1 - j
       do k=0, L_rule%Qdeg - j
          write(ifile, *) Fxi(it, 1:nbDim), wi(it, 1:ndim), &
               Dwi(it, 1:ndim, 0), Dwi(it, 1:ndim, 1), Dwi(it, 1:ndim, 2)
          it = it + nface - k
       enddo
       write(ifile,'(x)')
    enddo
    write(ifile,'(x)')

    ! third orientation
    do j=0, L_rule%Qdeg
       it = j + 1
       do k=0, L_rule%Qdeg - j
          write(ifile, *) Fxi(it, 1:nbDim), wi(it, 1:ndim), &
               Dwi(it, 1:ndim, 0), Dwi(it, 1:ndim, 1), Dwi(it, 1:ndim, 2)
          it = it + nface +1 - k
       enddo
       write(ifile,'(x)')
    enddo
    write(ifile,'(x)')


    deallocate(xi, Fxi, wi, Dwi, Dphi, DphiT)


  end subroutine PlotElemFunction3D


  !> 3D plot the derivative of the actual solution
  subroutine PlotElem_D_Function3D(ifile, elem,  dof2, w)
    integer, intent(in) :: ifile        ! number of the file
    type(element), intent (in) :: elem     ! element where the basis is given
    integer, intent(in) :: dof2
    real, dimension(1:dof2), intent(in) :: w
    real, dimension(:,:), allocatable :: xi, Fxi
    real, dimension(:,:), allocatable :: wi
    real, dimension(:,:,:), allocatable :: Dwi
    real, dimension(:,:), allocatable :: phiT
    real, dimension(:,:,:), allocatable :: DphiT, Dphi
    type(Lagrang_rule), pointer :: L_rule
    integer :: Qnum, Qdof, Qdeg,  j, j1,j2, k, nface, it

    !print*,'elem%i', elem%i
    nface = elem%deg + 4 !5

    L_rule => state%space%L_rule(nface)

    Qdeg = L_rule%Qdeg
    Qdof = L_rule%Qdof

    allocate(  xi(1:Qdof, 1:2), Fxi(1:Qdof, 1:2), wi(1:Qdof,1:ndim), Dwi(1:Qdof,1:ndim,1:nbDim) )

    xi(1:Qdof, 1:2) = L_rule%lambda(1:Qdof, 1:2)

    call ComputeF(elem, Qdof, xi(1:Qdof,1:nbDim), Fxi(1:Qdof,1:nbDim) )

    ! DG test functions in xi
    allocate( phiT(1:dof2, 1:Qdof) )
    allocate( DphiT(1:dof2, 1:nbDim, 1:Qdof), Dphi(1:dof2, 1:nbDim, 1:Qdof) )

    call PHI_orthonormal(Qdof, nbDim, xi(1:Qdof, 1:nbDim), 3, dof2, phiT(1:dof2, 1:Qdof), &
         DphiT(1:dof2, 1:nbDim, 1:Qdof) )

    do j=1, dof2
       if( elem%F%iFlin) then
          Dphi(j, 1, 1:Qdof) = elem%F%D1F0(1,1) * DphiT(j, 1, 1:Qdof) &
               + elem%F%D1F0(1,2)* DphiT(j, 2, 1:Qdof)

          Dphi(j, 2, 1:Qdof) = elem%F%D1F0(2,1) * DphiT(j, 1, 1:Qdof) &
               + elem%F%D1F0(2,2)* DphiT(j, 2, 1:Qdof)
       else
          ! first index  in elem%F%V%D1F( *, , ) is only AN APPROXIMATION !!!!
          Dphi(j, 1, 1:Qdof) = elem%F%V%D1F(  1   ,1,1) * DphiT(j, 1, 1:Qdof) &
               + elem%F%V%D1F( 1  ,1,2)* DphiT(j, 2, 1:Qdof)

          Dphi(j, 2, 1:Qdof) = elem%F%V%D1F( 1 ,2,1) * DphiT(j, 1, 1:Qdof) &
               + elem%F%V%D1F( 1 ,2,2)* DphiT(j, 2, 1:Qdof)

       endif

    enddo

    do j=1,Qdof
       do k=1, 1 !!ndim
          wi(j, k) = dot_product(phiT(1:dof2, j), w(1 : dof2) )
          Dwi(j, k, 1) = dot_product(Dphi(1:dof2, 1, j), w(1 : dof2) )
          Dwi(j, k, 2) = dot_product(Dphi(1:dof2, 2, j), w(1 : dof2) )
       enddo
    enddo


    ! first orientation
    it= 0
    do j=0, Qdeg
       do k=0, Qdeg - j
          it = it + 1
          write(ifile, *) Fxi(it, 1:nbDim), wi(it, 1:1), Dwi(it, 1:1, 1), Dwi(it, 1:1, 2)
       enddo
       write(ifile,'(x)')
    enddo
    write(ifile,'(x)')

   ! second orientation
    do j=0, L_rule%Qdeg
       it = nface + 1 - j
       do k=0, L_rule%Qdeg - j
          write(ifile, *) Fxi(it, 1:nbDim), wi(it, 1:1), Dwi(it, 1:1, 1), Dwi(it, 1:1, 2)
          it = it + nface - k
       enddo
       write(ifile,'(x)')
    enddo
    write(ifile,'(x)')

    ! third orientation
    do j=0, L_rule%Qdeg
       it = j + 1
       do k=0, L_rule%Qdeg - j
          write(ifile, *) Fxi(it, 1:nbDim), wi(it, 1:1), Dwi(it, 1:1, 1), Dwi(it, 1:1, 2)
          it = it + nface +1 - k
       enddo
       write(ifile,'(x)')
    enddo
    write(ifile,'(x)')


    deallocate(xi, Fxi, wi, Dwi, phiT, DphiT, Dphi)


  end subroutine PlotElem_D_Function3D


  ! plot the function given by its Basis coefficients w(i), i=1,...,dof
  subroutine PlotElemFunctionB(ifile, elem, Zrule, dof,  w )
    type(element), intent (in) :: elem     ! element where the basis is given
    integer, intent(in) :: dof, ifile
    character, intent(in) :: Zrule
    real, dimension(1:dof), intent(in) :: w
    real, dimension(:,:), allocatable :: xi
    real, dimension(:,:), pointer :: phi
    integer :: Qnum, Qdof, i

    write(ifile,*) '#  Output file generated by PlotElemFunctionB subroutine'
    write(ifile,*) '# Qnum = ',elem%Qnum,',  Qdof = ',elem%Qdof

    Qdof = elem%Qdof
    Qnum = elem%Qnum

    allocate( xi(1:Qdof,1:nbDim) )

    call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim), &
         xi(1:Qdof,1:nbDim) )

    if(Zrule == 'V') then
       phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)

    elseif(Zrule == 'L') then
       phi => state%space%L_rule(elem%Qnum)%phi(1:dof,1:Qdof)

    else
       print*,'Unknowb Zrule in eval_sol.f90'
       stop
    endif

    do i=1,Qdof
       write(ifile, *) xi(i, 1:nbDim), dot_product(w(1:dof), phi(1:dof,i) )
    enddo
    write(ifile,'(x)')

    deallocate(xi)

  end subroutine PlotElemFunctionB


  !> plot a function given by its values in integ nodes
  subroutine PlotElemFunctionQ(ifile, elem, Zrule, Qnum, Qdof,  wi )
    type(element), intent(in) :: elem
    integer, intent(in) :: Qdof, Qnum, ifile
    character, intent(in) :: Zrule
    real, dimension(1:Qdof), intent(in) :: wi
    real, dimension(:,:), pointer :: lambda
    real, dimension(:,:), allocatable :: xi
    integer :: i

    write(ifile,*) '#  Output file generated by PlotElemFunctionQ subroutine'

    if(Zrule == 'L') then
       lambda => state%space%L_rule(Qnum)%lambda(1:Qdof,1:nbDim)
    elseif(Zrule == 'V') then
       lambda => state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)
    else
       print*,'Unknown Zrule in eval_sol.f90'
       stop
    endif

    allocate(xi(1:Qdof, 1:nbDim))

    call ComputeF(elem, Qdof, lambda(1:Qdof,1:nbDim), xi(1:Qdof,1:nbDim) )

    do i=1,Qdof
       write(ifile, *) xi(i, 1:nbDim), wi(i), &
            2*(xi(i,1)*(1-xi(i,1))+xi(i,2)*(1-xi(i,2))), &
            wi(i) - 2*(xi(i,1)*(1-xi(i,1))+xi(i,2)*(1-xi(i,2)))
    enddo
    write(ifile,'(x)')
    write(ifile,'(x)')
    write(ifile,'(x)')

    deallocate(xi)

  end subroutine PlotElemFunctionQ

  ! plot the RTN function given by its Basis coefficients w(i), i=1,...,dof
  subroutine PlotElemRTNfunctionB(ifile, elem, Fdeg, Fdof,  w )
    type(element), intent (in) :: elem     ! element where the basis is given
    integer, intent(in) :: Fdof, ifile, Fdeg
    real, dimension(1:Fdof), intent(in) :: w
    real, dimension(:,:), allocatable :: xi
    real, dimension(:,:), pointer :: phi, phi1
    integer :: Qnum, Qdof, i

    write(*,*) ' PlotElemRTNfunctionB does not work, use of Piola trnasf is necessary'
    stop

    write(ifile,*) '#  Output file generated by PlotElemFunctionB subroutine'
    !write(ifile,*) '# Qnum = ',elem%Qnum,',  Qdof = ',elem%Qdof

    Qnum = state%RTN(Fdeg)%Qnum
    Qdof = state%RTN(Fdeg)%Qdof

    phi =>  state%RTN(Fdeg)%phi(1:Fdof,1, 0, 1:Qdof)
    phi1 => state%RTN(Fdeg)%phi(1:Fdof,2, 0, 1:Qdof)

    if(Fdof /= state%RTN(Fdeg)%dof) then
       print *,'Bad values in PlotElemRTNfunctionB  in eval_sol.f90'
       stop
    endif

    allocate( xi(1:Qdof,1:nbDim) )

    call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim), &
         xi(1:Qdof,1:nbDim) )

    do i=1,Qdof
       write(ifile, *) xi(i, 1:nbDim), dot_product(w(1:Fdof), phi(1:Fdof,i) ), &
            dot_product(w(1:Fdof), phi1(1:Fdof,i) )
    enddo
    write(ifile,'(x)')

    deallocate(xi)

  end subroutine PlotElemRTNfunctionB

  !> evaluation of test functions and their derivatives
  !> in an arbitrary  given set of integ nodes
  subroutine Eval_phi_Qnode(elem, dofA, Qdof, xi, phi, Dphi)
    type(element), intent(in):: elem      ! elem = element
    integer, intent(in) :: Qdof
    integer, intent(in) :: dofA
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
    real, dimension(1:dofA, 1:Qdof), intent(inout):: phi
    real, dimension(1:dofA, 1:nbDim, 1:Qdof), optional:: Dphi
    real, dimension(:,:,:), allocatable :: DphiK
    real, dimension(:,:,:), allocatable :: DF , D1F
    real, dimension(:), allocatable :: JF
    integer :: i, j, n, n1


    if(.not. present(Dphi)) then

       call PHI_orthonormal(Qdof, nbDim, xi(1:Qdof, 1:nbDim), 3, dofA, phi(1:dofA, 1:Qdof))

    else
       allocate(DphiK(1:dofA, 1:nbDim, 1:Qdof) )

       call PHI_orthonormal(Qdof, nbDim, xi(1:Qdof, 1:nbDim),3, dofA, phi(1:dofA, 1:Qdof),&
            DphiK(1:dofA, 1:nbDim, 1:Qdof)  )

       Dphi(:, :, :) = 0.

       if(elem%F%iFlin) then        ! linear element
          ! evaluation of derivative of test functions, only for linear elements
          do j=1, dofA
             do n=1,nbDim
                do n1 = 1, nbDim
                   Dphi(j, n, 1:Qdof) = Dphi(j, n, 1:Qdof)  &
                        + elem%F%D1F0(n, n1) * DphiK(j, n1, 1:Qdof)
                enddo
             enddo
          enddo
       else                         ! curved element
          allocate( DF(1:Qdof, 1:nbDim, 1:nbDim),  &
               D1F(1:Qdof, 1:nbDim, 1:nbDim),   JF(1:Qdof) )

          call ComputeDF(elem, Qdof, xi(1:Qdof, 1:nbDim), &
               DF(1:Qdof, 1:nbDim, 1:nbDim) )

          JF(1:Qdof) = DF(1:Qdof,1,1)*DF(1:Qdof,2,2) &
               - DF(1:Qdof,1,2)*DF(1:Qdof,2,1)

          !transpose and inverse of DF/Dx
          D1F(1:Qdof,1,1) =  DF(1:Qdof,2,2) / JF(1:Qdof)
          D1F(1:Qdof,2,1) = -DF(1:Qdof,1,2) / JF(1:Qdof)
          D1F(1:Qdof,1,2) = -DF(1:Qdof,2,1) / JF(1:Qdof)
          D1F(1:Qdof,2,2) =  DF(1:Qdof,1,1) / JF(1:Qdof)

          do j=1, dofA
             Dphi(j, 1, 1:Qdof) = &
                  D1F(1:Qdof,1,1) * DphiK(j, 1, 1:Qdof) &
                  + D1F(1:Qdof,1,2)* DphiK(j, 2, 1:Qdof)

             Dphi(j, 2, 1:Qdof) = &
                  D1F(1:Qdof,2,1) * DphiK(j, 1, 1:Qdof) &
                  + D1F(1:Qdof,2,2)* DphiK(j, 2, 1:Qdof)
          enddo

          !!write(*,'(a3,2i5,12es12.4)') 'new',elem%i,Qdof, Dphi(1:2, 1:2, 4:6)

          deallocate(DF, D1F, JF)

          deallocate(DphiK)

       endif !! if(elem%F%iFlin)

    endif  !! if (.not. present(Dphi))

  end subroutine Eval_phi_Qnode


  !> evaluation of the state vector \f$w\f$ on element \f$ elem \f$  in
  !> directly given nodes (by xi)
  subroutine Eval_w_Qnode(elem, Qdof, xi, t, wi, Dwt, Dwx)
    type(element), intent(in):: elem      ! elem = element
    integer, intent(in) :: Qdof
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
    real, intent(in) :: t       ! relative time between w^{k-1} and w^{k}
    real, dimension(1:Qdof,1:ndim), intent(inout):: wi, Dwt ! w and Dwt in integ. nodes
    real, dimension(1:Qdof,1:ndim, 1:nbDim), intent(inout):: Dwx ! Dwx in integ. nodes
    real, dimension(:,:,:), allocatable :: DphiK, DphiL
    real, dimension(:,:), allocatable :: phiL ! local store arrays
    integer :: dof, k, kst, i, j

    dof = elem%dof

    ! shape functions in new coordinates
    allocate(phiL ( 1:dof, 1:Qdof) )

    ! derivatives on real element
    allocate( DphiL(1:dof, 1:nbDim, 1:Qdof) )

    call  Eval_phi_Qnode(elem, dof, Qdof, xi, phiL, DphiL)

    do k=1,ndim              ! k = index of component of w
       kst = dof*(k-1) + 1
       do i=1,Qdof
          wi(i,k) = dot_product(elem%w(0,kst:kst+dof-1)*t + elem%w(1,kst:kst+dof-1)*(1-t),&
               phiL(1:dof,i) )

          Dwt(i,k) = dot_product(elem%w(0,kst:kst+dof-1) - elem%w(1,kst:kst+dof-1), &
               phiL(1:dof,i) ) / state%time%tau(1)

          Dwx(i, k, 1) = dot_product( &
               t*elem%w(0,kst:kst+dof-1) +(1-t)*elem%w(1,kst:kst+dof-1), DphiL(1:dof, 1, i) )

          Dwx(i, k, 2) = dot_product( &
               t*elem%w(0,kst:kst+dof-1) +(1-t)*elem%w(1,kst:kst+dof-1), DphiL(1:dof, 2, i) )

       enddo
    enddo

    deallocate(phiL, DphiL)

  end subroutine Eval_w_Qnode


  !> evaluate the divergence of the flux F (given in integ nodes)  in integ nodes
  subroutine EvalDiv_F(elem, Qnum, Qdof, xi, F, DivF )
    type(element), intent(in):: elem      ! elem = element
    integer, intent(in) :: Qnum, Qdof
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
    real, dimension(1:Qdof, 1:nbDim, 1:ndim), intent(in) :: F
    real, dimension(1:Qdof, 1:ndim), intent(out) :: DivF
    type(volume_rule), pointer :: V_rule
    real, dimension(:,:), allocatable :: temp, temp1 ! local store arrays
    real, dimension(:,:), allocatable :: phiL, Mass ! local store arrays
    real, dimension(:,:,:), allocatable :: DphiK, DphiL
    integer :: dof, k, kst, i,j

    V_rule => state%space%V_rule(Qnum)

    print*, 'Qdof = ', Qdof, 'V_rule%Qdof:' , V_rule%Qdof, 'Qnum:' , Qnum

    dof = elem%dof

    ! shape functions in new coordinates
    allocate(phiL ( 1:dof, 1:Qdof) )

    ! derivatives on real element
    allocate( DphiL(1:dof, 1:nbDim, 1:Qdof) )

    call  Eval_phi_Qnode(elem, dof, Qdof, xi, phiL, DphiL)


    ! mass matrix, we ignore the element area since G*area/ area = G in fact
    allocate(Mass(1:dof, 1:dof))
    do i=1,dof
       do j=1,dof
          Mass(i,j) = dot_product(V_rule%weights(1:Qdof),  phiL(i, 1:Qdof) * phiL(j, 1:Qdof) )
       enddo
    enddo

    call MblockInverse(dof, Mass)

    allocate(temp (1:2, 1:dof) )


    do k=1,ndim

       DivF(1:Qdof, k) = 0.

       do i=1,nbDim
          temp(1, 1:dof) = matmul(phiL(1:dof, 1:Qdof), F(1:Qdof, i, k)*V_rule%weights(1:Qdof) )
          temp(2, 1:dof) = matmul( Mass(1:dof, 1:dof), temp(1, 1:dof) ) ! basis coeffs

          do j=1,dof  ! divergence
             DivF(1:Qdof, k) = DivF(1:Qdof, k) +  DphiL(j, i, 1:Qdof)  * temp(2, j)
          enddo
       enddo
    enddo


    deallocate(phiL, DphiL, Mass, temp)
  end subroutine EvalDiv_F




  !> Fx  real physical coordinates of nodes outside elem
  !> xi barycentric coordinates of Fx  with respect elem
  subroutine BarycCoord(elem, Qdof1, Fx,  xi )
    type(element), intent(in) :: elem
    integer, intent(in) :: Qdof1
    real, dimension(1:Qdof1, 1:2), intent(in) :: Fx
    real, dimension(1:Qdof1, 1:2), intent(inout) :: xi
    real, dimension(:,:), allocatable :: xp
    real :: D, Dx, Dy
    integer :: i

    allocate(xp (1:3, 1:2) ) ! coordinates of vertices of elem
    xp(1, 1:2) = grid%x(elem%face(idx, 1), 1:2)
    xp(2, 1:2) = grid%x(elem%face(idx, 2), 1:2)
    xp(3, 1:2) = grid%x(elem%face(idx, 3), 1:2)

    if(elem%HGnode) then
       xp(1, 1:2) = grid%x(elem%face(idx,elem%HGvertex(1) ), 1:2 )
       xp(2, 1:2) = grid%x(elem%face(idx,elem%HGvertex(2) ), 1:2 )
       xp(3, 1:2) = grid%x(elem%face(idx,elem%HGvertex(3) ), 1:2 )
    endif

    D = (xp(2, 1) - xp(1, 1)) * (xp(3, 2) - xp(1, 2)) - (xp(2, 2) - xp(1, 2)) * (xp(3, 1) - xp(1, 1))
    ! cyclus of nodes
    do i=1, Qdof1
       Dx = (Fx(i, 1) - xp(1, 1))*(xp(3, 2) - xp(1, 2)) - (Fx(i, 2) - xp(1, 2))*(xp(3, 1) - xp(1, 1))
       Dy = (xp(2, 1) - xp(1, 1))*(Fx(i, 2) - xp(1, 2)) - (xp(2, 2) - xp(1, 2))*(Fx(i, 1) - xp(1, 1))

       xi(i, 1) = Dx / D
       xi(i, 2) = Dy / D
    enddo

  end subroutine BarycCoord


  !> recompute the function given in intge nodes to the basis coefficients
  subroutine Trans_Integ_Nodes_to_Basis(elem, dofP, ndimL, V_rule, wi, ww)
    class(element), intent(in), target :: elem
    integer, intent(in) :: dofP, ndimL
    type(volume_rule), intent(in), target :: V_rule
    real, dimension (1:ndimL, 1:V_rule%Qdof), intent(in) :: wi
    real, dimension (1:ndimL, 1:dofP), intent(inout) :: ww
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:), allocatable :: qq
!    real, dimension(:,:), allocatable :: Mass
    type( Mblock), allocatable :: Mass
    integer :: Qdof, k, Qnum

    Qdof = V_rule%Qdof
    phi => V_rule%phi(1:dofP, 1:Qdof)

    allocate(qq(1:dofP, 1:ndimL) ) ! RHS , function integrated by the test functions

    if ( DegFromDofTriang(dofP) > V_rule%Qdeg) &
      stop 'quadrature rule used for integration in Trans_Integ_Nodes_to_Basis is not accurate enough!'

    do k=1,ndimL
       ! QNUM ?
       call IntegrateVectorBplus(elem, V_rule%Qnum, dofP, wi(k, 1:Qdof), qq(1:dofP, k) )
!       write(*,'(a8,40es12.4)') 'qq:', qq(1:dofP, k)
    enddo

    if(dofP == elem%dof) then
       !do k=1,ndimL
       ww(1:ndimL, 1:dofP) = transpose( matmul(elem%MassInv%Mb(1:dofP,1:dofP), qq(1:dofP,1:ndimL) ) )
       !enddo
    elseif(dofP <= size(elem%Mass%Mb, 1) ) then

       allocate( Mass )
       call InitMblock( Mass,dofP,dofP )

       Mass%Mb(1:dofP, 1:dofP) = elem%Mass%Mb(1:dofP, 1:dofP)

       call SolveLocalMatrixProblem(dofP, Mass%Mb(1:dofP, 1:dofP), ndimL, qq(1:dofP,1:ndimL))
       ww(1:ndimL, 1:dofP) =  transpose ( qq(1:dofP,1:ndimL) )

       call DeleteMblock( Mass )
       deallocate(Mass)

    else
       allocate( Mass )
       call InitMblock( Mass,dofP,dofP )
       ! ??? QNUM, only one part of MASS , copy the first block ?.
       Qnum = DegFromDofTriang( dofP )*2
!       print*, 'Qnum=' , Qnum, V_rule%Qnum

       ! compute the bigger mass matrix
       call IntegrateBlockBBmass( elem, Qnum, dofP, dofP , Mass%Mb(1:dofP,1:dofP) )
       !       write(*,'(a8,40es12.2)') 'mass:', Mass%Mb( dofP,:)

       call SolveLocalMatrixProblem( dofP, Mass%Mb(1:dofP, 1:dofP), ndimL, qq(1:dofP,1:ndimL))
       ww(1:ndimL, 1:dofP) =  transpose ( qq(1:dofP,1:ndimL) )

       call DeleteMblock( Mass )
       deallocate( Mass )
    endif

    deallocate(qq)

  end subroutine Trans_Integ_Nodes_to_Basis

  !> projection of the DG-solution from elem%w(0,:) onto lower degree space using the corresponding
  !> operator stored in elem%block(:), Mb For the neighbouring elements we do not degree dof
  subroutine Energy_Elem_projection(elem)
    class(element), target, intent(inout) :: elem ! elem = element
    class(element), pointer  :: elem1 ! elem = element
    real, dimension(:, :), allocatable :: vec
    real, dimension(:, :), allocatable :: A
    integer :: dof, deg, dofM, degM, ndof, mdof, max_dof, dof1

    integer :: ip, i,j,k, ni

    deg = elem%deg
    dof = elem%dof

    ndof = dof*ndim
    allocate(A(1:ndof, 1:ndof) ) ! solution vector on elem and neighbours

    !print*,' max_dof = ', max_dof
    allocate(vec(0:1, 1:ndof) )   ! vec(0, :) solution, vec(1, :) RHS


    do ip = 1, 2
       degM = deg - ip

       if(degM < 0) then
          elem%wSS(1, ip,:) = 0.

       else

          dofM = ( degM + 1 ) * ( degM + 2 ) / 2

          ndof = dof*ndim

          A(1:ndof, 1:ndof) = elem%block(0)%Mb(1:ndof, 1:ndof)

          vec(1, 1:ndof) = matmul(A(1:ndof, 1:ndof), elem%w(0, 1:ndof))


          if(ndim > 1) &
               stop "bad setting of 'vec' in Energy_Elem_projection, skip rows and columns in A!!!"

          !ndof = dof*ndim  ! original solution
          ndof = dofM * ndim ! projection to lower order

          call SolveLocalMatrixProblem(ndof, A(1:ndof, 1:ndof), 1, vec(1, 1:ndof) )

          elem%wSS(1, ip,:) = 0.
          elem%wSS(1, ip, 1:ndof) = vec(1, 1:ndof)
          ! !!!!elem%wS(0, k*dof + 1 : k*dof + dofM) = vec(1, 1:ndof)
          !write(*,'(a8, 2i5, 40es12.4)') 'solution',0, ndof, vec(1,1:ndof)

          !print*,'----------------------------------------'

       endif

    enddo

    deallocate(vec, A)

    !stop" subroutine Energy_Elem_projection"
  end subroutine Energy_Elem_projection

  !> projection of the DG-solution from elem%w(0,:) onto lower degree space DEG using the corresponding
  !> operator stored in elem%block(:), Mb For the neighbouring elements we do not degree dof
  subroutine Energy_Elem_deg_projection(elem, degM, wp)
    class(element), target, intent(in) :: elem ! elem = element
    integer, intent(in):: degM
    real, dimension(elem%dof, 1:ndim), intent(inout):: wp
    real, dimension(:, :), allocatable :: vec
    real, dimension(:, :), allocatable :: A, AP
    integer :: dof, deg, dofM, ndof, mdof, max_dof, dof1

    integer :: ip, i,j,k, ni, nj, mi, mj

    deg = elem%deg
    dof = elem%dof

    ndof = dof*ndim
    allocate(A(1:ndof, 1:ndof) ) ! solution vector on elem and neighbours

    !print*,' max_dof = ', max_dof
    allocate(vec(0:1, 1:ndof) )   ! vec(0, :) solution, vec(1, :) RHS

    wp(:,:) = 0.

    if(degM >= 0) then

       dofM = ( degM + 1 ) * ( degM + 2 ) / 2

       ndof = dof*ndim

       A(1:ndof, 1:ndof) = elem%block(0)%Mb(1:ndof, 1:ndof)

       vec(1, 1:ndof) = matmul(A(1:ndof, 1:ndof), elem%w(0, 1:ndof))


       if(ndim > 1) &
               stop "verify ALGORITHM in in Energy_Elem_deg_projection, skip rows and columns in A!!!"

       !ndof = dof*ndim  ! original solution
       ndof = dofM * ndim ! projection to lower order
       allocate(AP(1:ndof, 1:ndof) )

       do i=1, ndim
          do j=1,ndim
             ni = (i-1)*dof
             nj = (j-1)*dof

             mi = (i-1)*dofM
             mj = (j-1)*dofM

             AP(mi+1: mi+dofM, mj+j : mj+dofM) =  A(ni+1: ni+dofM, nj+j : nj+dofM)
          enddo
       enddo


       call SolveLocalMatrixProblem(ndof, AP(1:ndof, 1:ndof), 1, vec(1, 1:ndof) )

       do k=1, ndim
          wP(1:ndof, k) = vec(1, (k-1)*dofM + 1 : k*dofM)
       enddo

       ! !!!!elem%wS(0, k*dof + 1 : k*dof + dofM) = vec(1, 1:ndof)
          !write(*,'(a8, 2i5, 40es12.4)') 'solution',0, ndof, vec(1,1:ndof)

          !print*,'----------------------------------------'

    endif

    deallocate(vec, A, AP)

    !stop" subroutine Energy_Elem_projection"
  end subroutine Energy_Elem_deg_projection


  subroutine Set_Elem_Coeffs_Decay(elem)
   type(element), target, intent(in):: elem ! elem = element
    real, dimension(:), pointer :: weights
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:,:), pointer :: Dphi
    real, dimension(:,:), allocatable :: wi, wii
    real, dimension(:,:,:), allocatable :: wwP   ! projections
    integer :: i, Qnum, Qdof, dof, dofL, deg, k, kst, l, ifile, ifile2, ifile3, ifile4
    real, dimension (:, :), allocatable :: ss
    real :: val, val1, sumy, sumyi, sumx, sumx2, sumy2, sumxy, order
    integer :: ni, min_deg
    logical :: iprint, singularity

    iprint = .false.
    !iprint = .true.

    call Detect_apriori_known_singularity(elem, singularity)
    if(singularity) iprint = .true.

    if( mod(elem%i, 1) == 0) iprint = .true.

    ! files outputs
    if(singularity) then
       ifile = state%space%adapt%adapt_level + 550
    else
       ifile = state%space%adapt%adapt_level + 500
    endif
    ifile2 = ifile + 100
    ifile3 = ifile + 200
    ifile4 = ifile + 300


    !if(elem%i == 1 .or. elem%i == 101 .or. elem%i == 111) iprint = .true.   ! mesh LL-shape.uns.365.grid
    !if(elem%i == 2 .or. elem%i == 38 ) iprint = .true.  ! mesh LL-shape.uns.100.grid


    dof = elem%dof

    ! projection of the solution of the element
    allocate(wwP(0:elem%deg, 1:dof, 1:ndim) )
    !if(iprint) &
    !   write(*,'(a6,3i5, 50es12.4)') 'proj W:',elem%i, elem%deg, deg, elem%w(0, 1:dof)
    do k=1, ndim
       wwP(elem%deg, 1:dof, k) = elem%w(0, (k-1)*dof + 1 : k*dof)
    enddo


    !write(*,'(a6,3i5, 50es12.4)') 'proj W:',elem%i, elem%deg, deg, wwP(elem%deg, 1:dof, 1)
    do deg = elem%deg, 0, - 1
       call Energy_Elem_deg_projection(elem, deg, wwP(deg, 1:dof, 1:ndim) )
    !!   if(iprint) then
    !!   write(*,'(a6,3i5, 50es12.4)') 'proj W:',elem%i, elem%deg, deg, wwP(deg, 1:dof, 1)
    !!   endif
    enddo

    !if(iprint)  print*

    ! estimate of the regularity

    allocate(ss(0:elem%deg, -1:7) )
    ss(:,:)= 0.


    Qdof = elem%Qdof
    weights => state%space%V_rule(elem%Qnum)%weights(1:Qdof)
    phi => state%space%V_rule(elem%Qnum)%phi(1:dof, 1:Qdof)

    ! derivatives of test functions in integ nodes
    allocate( Dphi(1:dof, 1:nbDim, 1:Qdof) )
    call Eval_Dphi(elem, dof, Dphi(1:dof, 1:nbDim, 1:Qdof) )

    allocate(wi(1:Qdof, 1:ndim) )


    do deg = elem%deg, 0, - 1
       ! "Legendre coefficients" ( for the first component (density) )
       k= deg*(deg+1)/2 + 1
       l = (deg+1)*(deg+2)/ 2

       ! basis coefficients corresponding to degre deg
       ss(deg, 0) =  sqrt( dot_product( wwP(elem%deg, k:l, 1),  wwP(elem%deg, k:l, 1)) )

       ! root test
       if(deg > 1) then
          !ss(deg, 1) = log( (2.*deg +1)/ (2 * ss(deg, 0)**2) ) / 2 / log(1.* deg) - 0.5
          !ss(deg, 1) =  ss(deg, 0)**(1./deg)
          ss(deg, 1) = log(1./ss(deg, 0)**2 ) / 1 / log(1.*deg)
       endif

       ! difference w_h and \Pi^E_{p-l} w_h
       do k=1, ndim

          ! L^2-norm
          wi(1:Qdof, k ) = matmul(wwP(elem%deg, 1:dof, k) - wwP(deg, 1:dof, k), phi(1:dof, 1:Qdof) )
          ss(deg, 2) = ss(deg, 2) + dot_product(weights(1:Qdof), wi(1:Qdof, k)**2 )

          ! H^1-semi-norm
          wi(1:Qdof, k ) =  &
               (matmul(wwP(elem%deg, 1:dof, k) - wwP(deg, 1:dof, k), Dphi(1:dof, 1, 1:Qdof) ))**2 + &
               (matmul(wwP(elem%deg, 1:dof, k) - wwP(deg, 1:dof, k), Dphi(1:dof, 2, 1:Qdof) ))**2

          ss(deg, 5) = ss(deg, 5) + dot_product(weights(1:Qdof), wi(1:Qdof, k)**2 )
       enddo

       ss(deg, 2) = sqrt( ss(deg, 2) )
       ss(deg, 5) = sqrt( ss(deg, 5) )
    enddo

    !ratio test
    do deg = elem%deg, 1, - 1
       val = ( abs(ss(deg, 0)) / max(1E-25, abs(ss(deg-1, 0)) ) )**2 !*(2*deg -1) / (2*deg + 1)

       ss(deg, -1) =  val !deg*(1.-val) / (1.+val)
       !write(*,'(a8, 2i5, 30es14.6)') ';39u39',elem%i, deg, val, ss(deg, 0), ss(deg-1, 0), ss(deg, -1)
    enddo



    ! convergence rates
    !write(*,'(a6, i5, 90es11.3)' )'dx39jde3',elem%deg, ss(0:elem%deg, 2)
    do deg = elem%deg-1, 1, - 1
       ! EOC in L^2-norm
       ss(deg, 3) = (log(ss(deg, 2) / ss(deg-1, 2)) - log(elem%diam) ) / log(1.*(deg ) / (deg+1))
       ss(deg, 4) = (log(ss(deg, 2) / ss(deg-1, 2)) ) / log(1.*(deg ) / (deg+1) )

       ! EOC in H^1-seminorm
       if(deg > 1) then
          ss(deg, 6) = (log(ss(deg, 5) / ss(deg-1, 5)) - log(elem%diam) ) / log(1.*(deg -1 ) / (deg))
          ss(deg, 7) = (log(ss(deg, 5) / ss(deg-1, 5)) ) / log(1.*(deg-1 ) / (deg) )
       endif

    enddo
    !write(*,'(a6, 2i5, 90es11.3)' )'dx39jde3',elem%i,elem%deg, ss(0:elem%deg, 2),ss(0:elem%deg, 3)
    !print*,'---'


    ! least squares for decays of energy projections in H^1-seminorm
    sumx = 0.; sumy = 0.; sumx2 = 0.; sumy2 = 0.; sumxy = 0.; ni = 0

    min_deg = max(2, elem%deg -2)  ! minimal value has to be at least 2
    min_deg = 0
    do deg = elem%deg, min_deg, - 1
       ! least squares:
       if(deg < elem%deg .and. deg >= 2 ) then
          ni = ni + 1
          sumx = sumx + log(1. * deg)
          sumx2 = sumx2 + log(1. * deg)**2
          sumy = sumy + log(ss(deg,5))
          sumy2 = sumy2 + log(ss(deg,5) )**2
          sumxy = sumxy + log(1. * deg) * log(ss(deg,5))
       endif
       !write(*,'(a8, i5, 30es12.4)') 'L R 23:',ni, 1.*deg, ss(deg,5), sumx, sumx2, sumy, sumy2, sumxy
    enddo

    if( ni * sumx2-sumx*sumx /= 0.) then
       order =  - (ni*sumxy -sumx*sumy)/(ni*sumx2-sumx*sumx) + 1.5
    else
       order = 1.
    endif

    if(elem%deg == 3) then
       order = 1.5 + log(ss(1,5) / ss(2, 5)) / log(2.)
    endif

    !write(ifile4, *) elem%i, elem%deg, elem%deg, order, elem%xc, val

    !if(ni >=2) stop "hj9eu-]32[\q"

    ! least squares for coeffs decay
    sumy = 0;  sumyi = 0.
    do i = 0, elem%deg
       sumy = sumy + log( abs(ss(i, 0) ) )
       sumyi = sumyi + i * log( abs(ss(i, 0) ) )
    enddo

    deg = elem%deg
    val = 6 *( 2*sumyi - deg * sumy) / deg/(deg+1)/(deg+2)
    !write(*,'(a8, 2i5, 30es12.4)' )'REG ###', elem%i, deg, elem%xc, val, 10**val

    write(ifile4, *) elem%i, elem%deg, elem%deg, order, elem%xc, 10**val, val


    do deg = elem%deg, 0, - 1
       if(iprint) then
          write(ifile, *) elem%i, elem%deg, elem%i - 0.1 * (elem%deg - deg), elem%xc, &
               ss(deg, -1), ss(deg, 0), 10**val, &
               ss(deg, 1:7)
               !max(1E-25, ss(deg, 1:7) )

          if(deg > 1)  &
          write(ifile2, *) elem%i, elem%deg, elem%i - 0.1 * (elem%deg - deg), elem%xc,&
               ss(deg, -1), ss(deg, 0), 10**val, &
               ss(deg, 1:7)
               !max(1E-25, ss(deg, 1:7) )

          if(deg > 0 ) & !.and. deg < elem%deg)  &
          write(ifile3, *) elem%i, elem%deg, elem%i - 0.1 * (elem%deg - deg), elem%xc,&
               ss(deg, -1), ss(deg, 0), 10**val, &
               ss(deg, 1:7)
               !max(1E-25, ss(deg, 1:7) )
       endif
    enddo


    write(ifile, '(x)')
    write(ifile2, '(x)')
    write(ifile3, '(x)')

    deallocate(ss, wwP, wi, Dphi)
  end subroutine Set_Elem_Coeffs_Decay


  !> compute the values of the all derivatives of elem%wS(1:ndimL, 1:elem%dof) in BARYCENTERS
  !> results stored in elem%wSD(1:ndimL, i), i = 1 ... (deg+1)*(deg+2) / 2
  !> itype = 1, => deg = elem%deg
  !> itype = 2, => deg = elem%HO_deg
  subroutine Compute_means_values_DW(gridN, ndimL, itype)
    type(mesh), intent(inout), target  :: gridN
    integer, intent(in) :: ndimL
    integer, intent(in) :: itype
    class(element), pointer :: elem
    real, dimension(:,:), allocatable :: psi    ! canonical basis functions
    real, dimension(:), allocatable :: weights ! Inversion of local mass matrix of order deg+1!!
    real, dimension(:,:), allocatable :: Mass ! Local mass matrix of order deg+1!!
    real, dimension(:,:), allocatable :: MassInv ! Inversion of local mass matrix of order deg+1!!
    real, dimension(:,:), allocatable ::  vec  ! temporary array
    real, dimension(:),   allocatable ::    wi ! temporary array

    real, dimension(:,:), allocatable ::  Fx  ! temporary array

    real, dimension(:,:), pointer :: phi
    real, dimension(:,:), allocatable :: val
    real :: norm, rphiD
    integer :: ie, i, j,k, iphi, iphiD
    integer ::  deg, Qdof, degP, dofP, degS, ideg
    type(volume_rule), pointer :: V_rule


    do ie = 1, gridN%nelem
       elem => gridN%elem(ie)

       if(itype == 1) then
          degP = elem%deg
          dofP = elem%dof

       elseif(itype == 2) then
          degP =  elem%HO_deg
          dofP = DOFtriang( degP )  ! dof of the HO reconstruction

       endif


       allocate(elem%wSD(1:2, 1:ndimL, 1:dofP) ) ! 1st idx = 1 coeffs, 2nd idx = derivatives
       elem%wSD = 0.

       V_rule => state%space%V_rule(elem%Qnum)
       Qdof = V_rule%Qdof

       ! computation of the HO derivatives by the projection
       allocate( psi(1:dofP, 0:Qdof) )

       call Eval_Canonical_Functions(elem, degP, dofP, V_rule, psi(1:dofP, 0:Qdof) )


       ! mass matrix
       allocate(weights(1:Qdof) )
       call Eval_V_Weights_plus(elem, V_rule, weights(1:Qdof) )

       allocate( Mass(1:dofP, 1:dofP) )
       allocate( vec ( 1:dofP, 1:ndimL ) )
       allocate( wi (1:Qdof) )

       phi => V_rule%phi(1:dofP, 1:Qdof)


       do i=1,dofP
          do j=1,i
             !if(state%space%adapt%adapt_level >= 7) then
             !   write(*,'(a8,i5, 300es12.4)') '### we', i, weights(1:Qdof)
             !   write(*,'(a8,i5, 300es12.4)') '### ps', j, psi(j, 1:Qdof)
             !endif


             Mass(i,j) = dot_product(weights(1:Qdof) *  psi(i, 1:Qdof), psi(j, 1:Qdof) )
             Mass(j,i) = Mass(i,j)
          enddo

          do k = 1, ndimL

             !print*,'DegP:', ie, degP, dofP, Qdof,k,ndimL
             wi(1:Qdof) = matmul( elem%wS(k, 1:dofP), phi(1:dofP, 1:Qdof) )

             vec(i, k) = dot_product(weights(1:Qdof) *  psi(i, 1:Qdof), wi(1:Qdof) )
          enddo
       enddo


       !do l=1, dofP
       !    write(*,'(i5, 30es12.4)') l, Mass(l, 1:dofP)
       ! enddo
       ! print*
       ! do ideg=ideg_min, ideg_max
       !    write(*,'(i5, 30es12.4)') ideg, Dww (ideg, 1, 1:dofP)
       ! enddo
       ! print*
       ! do ideg=ideg_min, ideg_max
       !    write(*,'(i5, 30es12.4)') ideg, vec(ideg, :, 1)
       ! enddo
       ! print*,'------------------------------------------------'

       ! evaluation of the function in canonical basis

       allocate( MassInv(1:dofP, 1:dofP) )

       MassInv(1:dofP, 1:dofP) = Mass(1:dofP, 1:dofP)

       call SolveLocalMatrixProblem(dofP, MassInv(1:dofP, 1:dofP), ndimL, vec(1:dofP,1:ndimL))

       ! the basis coefficients in the canonical basis
       do k=1,ndimL
          elem%wSD(1, k, 1:dofP) = vec(1:dofP, k)
       enddo

       ! derivatives in the barycentre
       do iphi =1, dofP
          call Deriv_idx_reverse(iphi, i, j)
          elem%wSD(2, 1:ndimL, iphi) = vec(iphi, 1:ndimL) * factorial(i) * factorial(j)/elem%diam**(i+j)
       enddo


       !write(*,'(a8, i5, 80es12.4)') 'sols:', degP, vec(1:dofP, 1)
       !write(*,'(a8, i5, 80es12.4)') 'ders:', degP, elem%wSD(2,1, 1:dofP)
       !print*,
       !deallocate(val)


       !allocate(Fx(1:Qdof, 1:2) )
       !call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )
       !
       !do i=1, Qdof
       !   write(22, *) Fx(i, 1:2),  dot_product( vec(1:dofP, 1), psi(1:dofP, i) )
       !enddo
       !deallocate(Fx)

       ! the MEAN VALUES
       if(itype == 2) then
          ! computing of the mean values of the derivatives
          allocate(val (1:ndimL, 1:Qdof))
          do ideg = 0, degP   ! order of the derivative
             do i= ideg, 0, -1  ! d^ideg / d x^i / d y^{ideg-i}
                val(:,:) = 0.

                do iphi = 1, dofP  ! index of the basis function
                   ! computing of the derivatives in ineg nodes
                   call Set_canonical_basis_index_deriv(ideg, i, iphi, iphiD, rphiD)

                   ! the correct normalization since iphi and iphiD has different degress of elem%diam
                   rphiD = rphiD  / elem%diam**ideg

                   if(iphiD > 0) then
                      do k=1,ndimL
                         val(k, 1:Qdof) = val(k, 1:Qdof) + vec(iphi, k) * psi(iphiD, 1:Qdof) * rphiD
                         !    ------------ correct -------     ^^^^  -------  ^^^^^

                         ! verifying !!
                         ! if(ideg == 2 .and. i == 2) then
                         !    do j=1, Qdof
                         !       write(*,'(a6, 5i5, es12.4, a2, 10es12.4) ') &
                         !            'D u_h', ideg, i, ideg - i , iphi, iphiD, rphiD, &
                         !            '|',  vec(iphi, k), psi(iphiD, j), val(k, j)
                         !    enddo
                         !    print*,'---------------____________'
                         ! endif

                      enddo

                   endif
                enddo

                ! computing of the mean values of the derivatives
                do k=1,ndimL
                   j = Deriv_idx(i, ideg - i)
                   elem%wSD(2, k, j) = dot_product(weights(1:Qdof), val(k, 1:Qdof) )
                enddo

             enddo

          end do

          deallocate(val)


          elem%wSD(2, 1:ndimL, 1:dofP) = elem%wSD(2, 1:ndimL, 1:dofP) / elem%area

       endif

       !if( sqrt(dot_product(elem%xc, elem%xc)) < 0.2) then
       !   write(*,'(a8, i5, 30es12.4)') 'Derivs:',elem%i, elem%wSD(2, 1, 1:dofP)
       !   print*,'----------------------------------', elem%xc(:)
       !endif




       deallocate(Mass, vec, wi, psi)

       deallocate(MassInv, weights)

       !stop '73d93djjsdyu923w'

    enddo  ! do i =1 gridN%nelem

    !if(elem%xc(1) > 0.95 .and. elem%xc(2) > 0.95) &
    !     print*,'-----------------------------------'




  end subroutine Compute_means_values_DW


  !> evaluate the canonical basis functions in integ nodes
  !> \f$ \psi_{ij} = \frac{1}{h_K^{i+1}} (x-x_c)^i (y - y_c)^j \f$
  !> \f$ \mbox{psi}(i, 0) = \int_K \psi_i d x\f$ !!!!
  subroutine Eval_Canonical_Functions(elem, degP, dofP, V_rule, psi)
    class(element), intent(inout):: elem
    integer,  intent(in):: degP, dofP
    type(volume_rule), intent(in) :: V_rule
    real, dimension(1:dofP, 0:V_rule%Qdof), intent(inout) :: psi  ! (:, 0) the MEAN VALUE !!!
    real, dimension(:,:), allocatable :: Fx     ! physical coordinates in integ nodes
    integer :: i,j,k, Qdof
    real :: val

    Qdof = V_rule%Qdof

    allocate(Fx(1:Qdof, 1:2) )
    call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )

    Fx(1:Qdof, 1) = (Fx(1:Qdof, 1) - elem%xc(1)) / elem%diam
    Fx(1:Qdof, 2) = (Fx(1:Qdof, 2) - elem%xc(2)) / elem%diam

    ! canonical basis functions (x- xc)**j (y - yc)**(i-j)
    call Canonical_Functions(degP, dofP, Qdof, Fx(1:Qdof, 1:2), psi(1:dofP, 1:Qdof) )

    ! the basis function has to be mean value free, i.e., \int_K \psi(k, .) dx = 0,  k>1
    ! no influence to the derivatives
    do k= 2, dofP
       ! integral over elem
       call IntegrateFunction_plus(elem, V_rule, psi(k, 1:Qdof), val)
       psi(k , 0) = val / elem%area
       psi(k, 1:Qdof) = psi(k, 1:Qdof) -  psi(k , 0)

    end do
    psi(1 , 0) = 0.   ! no shift

    deallocate(Fx)
  end subroutine Eval_Canonical_Functions

  !> canonical basis functions (x- xc)**j (y - yc)**(i-j),see subroutine Eval_Canonical_Functions
  subroutine Canonical_Functions(degP, dofP, Qdof, Fx, psi )
    integer,  intent(in):: degP, dofP, Qdof
    real, dimension(1:Qdof,1:nbDim), intent(in) :: Fx     ! physical coordinates in integ nodes
    real, dimension(1:dofP, 1:Qdof), intent(inout) :: psi
    integer :: k, i, j

    k = 0
    do i=0, degP
       do j=0, i
          k = k + 1
          if( i + j == 0) then
             psi(k, 1:Qdof) = 1.

          elseif( j == 0) then
             psi(k, 1:Qdof) = Fx(1:Qdof, 1)**(i - j)

          elseif( j - i == 0) then
             psi(k, 1:Qdof) = Fx(1:Qdof, 2)**j

          else
             psi(k, 1:Qdof) = Fx(1:Qdof, 1)**(i - j) * Fx(1:Qdof, 2)**j

          endif

          !psi(k, 1:Qdof) =  psi(k, 1:Qdof) !/ elem%diam**i
          !do l=1, Qdof
          !   write(100+k, *) Fx(l, 1:2) + elem%xc(1:2), psi(k,l)
          !enddo

       enddo
    enddo

  end subroutine Canonical_Functions



  !> compute the derivatives of the canonical basis
  !> \f$ \psi_{ij} = \frac{1}{h_K^{i+1}} (x-x_c)^i (y - y_c)^j \f$
  !> obviously \f$ \frac{d^{n}}{d x^i d y^{n - i}} \psi_{iphi} =  rphi \[psi_{iphiD} \f$
  subroutine Set_canonical_basis_index_deriv(n, i, iphi, iphiD, rphiD)
    integer, intent(in) :: n, i, iphi
    integer, intent(inout) :: iphiD
    real, intent(inout) :: rphiD
    integer :: ix, iy     !  i,j, indexes of iphi function (according Eval_Canonical_Functions)
    integer :: iDx, iDy   !  i,j, indexes of Diphi
    integer :: ider_x, ider_y   !  i,j, indexes of the derivative

    if(n == 0) then ! no derivative
       iphiD = iphi
       rphiD = 1

    else
       ! indexes of the canonical basis function

       call Deriv_idx_reverse(iphi, ix, iy)

       ! indexes of the actual derivative
       ider_x = i
       ider_y = n-i

       iDx = ix - ider_x
       iDy = iy - ider_y

       if(iDx >= 0 .and. iDy >= 0) then

          iphiD = Deriv_idx(iDx, iDy)

          rphiD = factorial(ix) / factorial(iDx) * factorial(iy) / factorial(iDy)
       else
          iphiD = 0
          rphiD = 0.
       endif

    endif

    !write(*,'(a6, 4i5, a2, 5i5, es12.4) ') 'idx D', n, ider_x, ider_y, iphi, '|', &
    !     ix, iy, iDx, iDy, &
    !     iphiD, rphiD

  end subroutine Set_canonical_basis_index_deriv

end module eval_sol
