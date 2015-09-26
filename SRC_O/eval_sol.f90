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
  public:: Eval_vec_func_Elem
  public:: Eval_res_func_Elem
  public:: Eval_w_Elem
  public:: Eval_w_Elem_time
  public:: Eval_w_Elem_at_time
  public:: Eval_w_t_Elem
  public:: Eval_w_w_Elem
  public:: Eval_Dw_Elem
  public:: Eval_Dw_Elem_time
  public:: Eval_Dw_Elem_at_time
  public:: Eval_aver_w_Elem
  public:: Eval_aver_Dw_Elem
  public:: Eval_wht_Elem

  public:: Eval_wAEE_Elem
  public:: Eval_DwAEE_Elem


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
  public:: PlotElemFunctionB
  public:: PlotElemFunctionQ
  
  public:: PlotElemFunction3D
  public:: PlotElem_D_Function3D

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
  subroutine Eval_w_Elem_plus(elem, V_rule, Dwi)
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
  end subroutine Eval_w_Elem_plus


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

  !> evaluation of the computational solution based on AEE ALG (state vector) \f$wc\f$ on element \f$ elem \f$  in
  !> integ nodes, i.e. recomputation of elem%wc into wi in volume integ. nodes
  subroutine Eval_wAEE_Elem(elem, wi)
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
          wi(i,k) = dot_product(elem%wc(0, kst:kst+dof-1), phi(1:dof,i) )
       enddo
    enddo

  end subroutine Eval_wAEE_Elem


  !> evaluation of the derivatives of the computational solution based on AEE ALG (state vector) \f$wc\f$
  !> on element \f$ elem \f$  in integ nodes
  subroutine Eval_DwAEE_Elem(elem, Dwi)
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
          Dwi(i, k, 1) = dot_product(elem%wc(0,kst:kst+dof-1), Dphi(1:dof, 1, i) )
          Dwi(i, k, 2) = dot_product(elem%wc(0,kst:kst+dof-1), Dphi(1:dof, 2, i) )
       enddo
    enddo
    deallocate(Dphi)
  end subroutine Eval_DwAEE_Elem


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
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
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
         elem%xi(0,1:Qdof, 1:nbDim))

    deallocate(wi )
  end subroutine Eval_f_s_Elem

  !> evaluate the inviscid fluxes f_s
  !> in integ nodes on  elem t, t=1 => t_k, t=0 => t_{k-1}
  subroutine Eval_f_s_Elem_at_time(Set_f_s, elem,  f_s, t)
    interface
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
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
         elem%xi(0,1:Qdof, 1:nbDim))

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
      subroutine Set_R_s(ndimL, nbDim, Qdof, w, Dw, Re_1, R_s)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
      end subroutine Set_R_s
    end interface

    type(element), intent (in) :: elem   ! element
    real, dimension(1:elem%Qdof, 1:nbDim, 1:ndim), intent(out) :: R_s ! output
    real, dimension(:,:), allocatable :: wi
    real, dimension(:,:,:), allocatable :: Dwi
    real, dimension(:), allocatable :: Re_1
    integer :: Qdof, i, j


    Qdof = elem%Qdof

    allocate( wi(1:Qdof,1:ndim),  Dwi(1:Qdof, 1:ndim, 1:nbDim), Re_1(1:Qdof)  )
    Re_1(:) = state%model%Re1

    !print*, 'ctime in Eval_R_s_Elem ', state%time%ctime

    call Eval_w_Elem(elem,  wi)
    call Eval_Dw_Elem(elem, Dwi(1:Qdof,1:ndim, 1:nbDim) )

!    print*, 'WI:' , wi(1,1:ndim)
!    print*, 'DWI:', Dwi(1,1:ndim, 1:nbDim)

    call Set_R_s(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof, 1:ndim, 1:nbDim), &
         Re_1(1:Qdof), R_s(1:Qdof, 1:nbDim, 1:ndim) )

    deallocate(wi, Dwi, Re_1 )

  end subroutine Eval_R_s_Elem


  !> evaluate the viscous fluxes R_s
  !> \f$ \sum_{s=1}^d R_s(w,\nabla w) n_s  \f$
  !> in integ nodes on elem, at time t, t=1 => t_k, t=0 => t_{k-1}
  subroutine Eval_R_s_Elem_at_time(Set_R_s, elem,  R_s, t )
    interface
      subroutine Set_R_s(ndimL, nbDim, Qdof, w, Dw, Re_1, R_s)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
      end subroutine Set_R_s
    end interface
    type(element), intent (in) :: elem   ! element
    real, intent(in) :: t
    real, dimension(1:elem%Qdof, 1:nbDim, 1:ndim), intent(out) :: R_s ! output
    real, dimension(:,:), allocatable :: wi
    real, dimension(:,:,:), allocatable :: Dwi
    real, dimension(:), allocatable :: Re_1
    integer :: Qdof, i, j

    Qdof = elem%Qdof

    allocate( wi(1:Qdof,1:ndim),  Dwi(1:Qdof, 1:ndim, 1:nbDim), Re_1(1:Qdof)  )
    Re_1(:) = state%model%Re1

    call Eval_w_Elem_at_time(elem,  wi, t)
    call Eval_Dw_Elem_at_time(elem, Dwi(1:Qdof,1:ndim, 1:nbDim), t )

    call Set_R_s(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof, 1:ndim, 1:nbDim), &
         Re_1(1:Qdof), R_s(1:Qdof, 1:nbDim, 1:ndim) )

    deallocate(wi, Dwi, Re_1 )

  end subroutine Eval_R_s_Elem_at_time

  !> evaluate the viscous fluxes R_s
  !> \f$ \sum_{s=1}^d R_s(w,\nabla w) n_s  \f$
  !> in integ nodes on the ie-th edge of elem,
  !> is opposite = .true., the orientation \f$ (n_1,\dots,n_d)\f$ is opposite
  subroutine Eval_R_s_Edge(Set_R_s, elem, ie, Rflux, opposite)
    interface
      subroutine Set_R_s(ndimL, nbDim, Qdof, w, Dw, Re_1, R_s)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
      end subroutine Set_R_s
    end interface
    type(element), intent (in) :: elem   ! element
    integer, intent(in) :: ie            ! index of the edge
    real, dimension(1:elem%face(fGdof,ie), 1:ndim), intent(out) :: Rflux ! output
    logical, intent(in) :: opposite
    real, dimension(:,:), allocatable :: wi
    real, dimension(:,:,:), allocatable :: Dwi
    real, dimension(:,:,:), allocatable :: R_s
    real, dimension(:), allocatable :: Re_1
    integer :: Qdof, i, j

    Qdof = elem%face(fGdof,ie)

    allocate( wi(1:Qdof,1:ndim),  Dwi(1:Qdof,1:ndim,1:nbDim), Re_1(1:Qdof)  )
    allocate( R_s(1:Qdof, 1:nbDim, 1:ndim)  )
    Re_1(:) = state%model%Re1

    call Eval_w_Edge(elem, ie,  wi, opposite)
    call Eval_Dw_Edge(elem,  ie,  Dwi(1:Qdof,1:ndim,1:nbDim), opposite)

    call Set_R_s(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof, 1:ndim, 1:nbDim), &
         Re_1(1:Qdof), R_s(1:Qdof, 1:nbDim, 1:ndim) )

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
      subroutine Set_R_s(ndimL, nbDim, Qdof, w, Dw, Re_1, R_s)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
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

    nface = elem%deg + 3 !5

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
    real, dimension(:,:), allocatable :: phiT
    type(Lagrang_rule), pointer :: L_rule
    integer :: Qnum, Qdof, Qdeg,  j, j1,j2, k, nface, it

    nface = elem%deg + 4 !5

    L_rule => state%space%L_rule(nface)
    
    Qdeg = L_rule%Qdeg
    Qdof = L_rule%Qdof
    
    allocate(  xi(1:Qdof, 1:2), Fxi(1:Qdof, 1:2), wi(1:Qdof,1:ndim), Dwi(1:Qdof,1:ndim,1:nbDim) )

    xi(1:Qdof, 1:2) = L_rule%lambda(1:Qdof, 1:2) 

    call ComputeF(elem, Qdof, xi(1:Qdof,1:nbDim), Fxi(1:Qdof,1:nbDim) )

    ! DG test functions in xi
    allocate( phiT(1:dof2, 1:Qdof) ) 

    call PHI_orthonormal(Qdof, nbDim, xi(1:Qdof, 1:nbDim), 3, dof2, phiT(1:dof2, 1:Qdof) )

    do j=1,Qdof
       do k=1, 1 !!ndim
          wi(j, k) = dot_product(phiT(1:dof2, j), w(1 : dof2) )
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
       Dphi(j, 1, 1:Qdof) = elem%F%D1F0(1,1) * DphiT(j, 1, 1:Qdof) +elem%F%D1F0(1,2)* DphiT(j, 2, 1:Qdof)

       Dphi(j, 2, 1:Qdof) = elem%F%D1F0(2,1) * DphiT(j, 1, 1:Qdof) +elem%F%D1F0(2,2)* DphiT(j, 2, 1:Qdof)
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


end module eval_sol
