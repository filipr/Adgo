!> space high-order  iterpolation on unstructured grids
module ama_L2interpol

  use main_data
  use problem_oper
  use geometry
  use mesh_oper
  use st_interpol

  implicit none

  public:: SimpleInterpolDGsolution
  public:: AdvancedInterpolDGsolution

  public:: LeastSquareInterpolationH1
  public:: LeastSquareInterpolationL2
  public:: LeastSquareInterpolationWeighted
  public:: LeastSquareRecovery_Dw_mean_values
  public:: LeastSquareInterpolationNodes
  public:: LeastSquareInterpolationPP
  public:: P_0_to_P2_Interpolation

contains

 !> a simple reconstruction of a high order approximation from gridO to gridN,
 !> new grid: pointer gridN,
 !> old grid: pointer gridO;
 !> subroutine uses array AMA%iaegr computed by ANGENER
 !> the solution from the old grid is set to integ nodes of the new one
 subroutine SimpleInterpolDGsolution(gridN, gridO )
    type(mesh), intent(inout), target	:: gridN, gridO
    class(element), pointer :: elem, elemO
    integer :: Qdof, dof, max_tri, si, sj, ssi, ssj, ss_plus, dofO
    integer :: ie, ip, i, j, ii, jj, k, kst, l,  nsq
    real, dimension(:,:), allocatable :: xi
    real, dimension(:,:,:), allocatable :: wi
    real, dimension(1:2, 1:2) :: rmx
    integer, dimension(:,:,:), allocatable :: itri
    integer, dimension(:,:), allocatable :: ntri
    real, dimension(:,:), allocatable :: xt !
    real, dimension(:,:), allocatable :: phi ! local store arrays
    real, dimension(:), allocatable :: qi, qj
    real, dimension(:,:), pointer :: phiL ! local store arrays
    real :: x0(1:2)
    integer:: k_low, k_top


    k_low = 0
    k_top = state%time%deg+1
    !k_top = 0

    !integer, dimension(:), allocatable :: itli
    !real, dimension(:,:), allocatable :: tlr

    allocate(xt(1:3, 1:2) ) ! coordinates of vertices of elem

    ! frame of the computational domain
    rmx(1, 1) = minval(gridO%x(:, 1) )
    rmx(1, 2) = maxval(gridO%x(:, 1) )
    rmx(2, 1) = minval(gridO%x(:, 2) )
    rmx(2, 2) = maxval(gridO%x(:, 2) )

    ! splitting of the domain onto nsq x nsq square cells
    nsq = (gridO%nelem/2)**0.5
    !!!!nsq = 15

    !max_tri = 250
    max_tri = gridO%nelem

    allocate( ntri(1:nsq, 1:nsq), itri(1:nsq, 1:nsq, 1: max_tri) )
    ntri(:,:) = 0

    ! seeking of the corresponding triangles from the old grid to the square cells
    do ie=1, gridO%nelem
       elem => gridO%elem(ie)

       !write(99,*) elem%xc(:)

       i = int(nsq * (elem%xc(1) - rmx(1, 1) ) /  (rmx(1,2) - rmx(1, 1) )) + 1
       j = int(nsq * (elem%xc(2) - rmx(2, 1) ) /  (rmx(2,2) - rmx(2, 1) )) + 1
       ntri(i, j) = ntri(i, j) + 1
       k = ntri(i,j)
       if( k <= max_tri) then
          itri(i,j, k) = ie
       else
          print*,' problem in anisot.f90: k > max_tri', i,j,ie, k, max_tri, nsq
          stop
       endif
    enddo

    !do i=1,nsq
    !   do j=1,nsq
    !      k = ntri(i,j)
    !      write(*,'(20i5)') i,j, k, idx(i,j,1:k)
    !   enddo
    !enddo


    do ie=1, gridN%nelem
       elem => gridN%elem(ie)

       if(elem%to_recompute) then
          Qdof = elem%Qdof
          allocate( xi(1:Qdof, 1:nbDim), wi(k_low:k_top, 1:Qdof, 1:ndim) )
          !integration nodes on K
          call ComputeF(elem, Qdof, state%space%V_rule(elem%Qnum)%lambda(1:Qdof,1:nbDim), &
               xi(1:Qdof, 1:nbDim) )

          !allocate(itli(1:Qdof), tlr(Qdof, 1:3) )
          !itli(:) = 0
          !tlr(:,:) = 0.

          !for each integ node, we seek the square element
          do ip=1,Qdof
             !write(98,*) xi(ip, 1:2), elem%i, ip
             ! starting relative coordinates of
             ssi = int(nsq * (elem%xc(1) - rmx(1, 1) ) /(rmx(1,2) -rmx(1, 1) ))+ 1
             ssj = int(nsq * (elem%xc(2) - rmx(2, 1) ) /(rmx(2,2) -rmx(2, 1) ))+ 1


             do ss_plus = 0, nsq
                do ii = -ss_plus, ss_plus
                   do jj = -ss_plus, ss_plus
                      ! checking of the correct square cell
                      if(max(abs(ii), abs(jj)) == ss_plus) then
                         si = ssi + ii
                         sj = ssj + jj

                         if(si >= 1 .and. si <= nsq  &
                              .and. sj >= 1 .and. sj <= nsq) then
                            ! we go through the list of elements

                            do l=1, ntri(si, sj)
                               elemO => gridO%elem( itri(si, sj, l) )

                               xt(1, 1:2) = gridO%x(elemO%face(idx, 1), 1:2)
                               xt(2, 1:2) = gridO%x(elemO%face(idx, 2), 1:2)
                               xt(3, 1:2) = gridO%x(elemO%face(idx, 3), 1:2)

                               call BarycCoordOne(xt(1:3, 1:2), xi(ip, 1:2), x0(1:2))

                               !write(*,'(a4,5i4,5es12.4)')'...',ie, ip, si, sj,  idx(si, sj, l), &
                               !     xi(ip, 1:2),  x0(1:2), 1. - xi(ip,1) - xi(ip, 2)

                               ! the element found ??
                               if(x0(1) >= 0. .and. x0(2) >= 0. .and. x0(1) + x0(2) <= 1.) &
                                    goto 10

                            enddo ! l=1, ntri(si, sj)
                         endif
                      endif
                   enddo ! jj
                enddo ! ii
             enddo ! ss_plus = 0, nsq
10           continue

             !write(*,'(a4,6i4, 5es12.4)') &
             !     '###',ie, ip, ss_plus, si, sj,  elemO%i , &
             !     xi(ip, 1:2),  x0(1:2), 1. - x0(1) - x0(2)

             ! evaluation of the solution at the integ. node
             dofO = elemO%dof
             allocate(phi ( 1:dofO, 1:1) )
             !test function in x0
!!!call PHI_orthonormal(1, x0(1:nbDim), 3, dofO, phi(1:dofO, 1:1) )
             call Eval_phi_Qnode(elemO, dofO, 1, x0(1:nbDim), phi(1:dofO, 1:1) )

             do k=1,ndim              ! k = index of component of w
                kst = dofO*(k-1) + 1
                do i=k_low, k_top
                   wi(i,ip,k) = dot_product(elemO%w(i,kst:kst+dofO-1), phi(1:dofO,1))
                enddo
             enddo
             deallocate (phi)

             !write(97, *) xi(ip, 1:2), wi(ip, 1:ndim)

          enddo !! ip=1,Qdof

          ! recomputation of the solution from the integ nodes to the
          ! basis coefficients
          dof = elem%dof
          allocate( qi(1:dof), qj(1:dof) )

          ! evaluating of the coefficients of basis expansion
          do k=1, ndim
             do i=0,0
                qi(:) = 0.
                call IntegrateVectorB(elem, dof, wi(i, 1:Qdof, k), qi(1:dof) )

                qj(1:dof) = matmul(elem%MassInv%Mb(1:dof,1:dof), qi(1:dof) )


                elem%w(i, (k-1)*dof + 1: k*dof) = qj(1:dof)
             enddo
          enddo

          deallocate( xi, wi, qi, qj )
       end if  ! if(elem%to_recompute)
    enddo  ! ie = 1, gridN%nelem

    deallocate(xt)

  end subroutine SimpleInterpolDGsolution



  !> an advanced reconstruction of a high order approximation from gridO to gridN,
  !> new grid: pointer gridN,
  !> old grid: pointer gridO;
  !> subroutine uses array AMA%iaegr computed by ANGENER
  !> the solution from the old grid is set to integ nodes of the new one
  !> Similarly as in SimpleInterpolDGsolution,
  !> but with an adaptive quadrature in the projection
  !> IF type_of_WW IS PRESENTED THEN :
  !> if type_of_WW = 1 RECOMPUTATION from elem%w TO ARRAYS ELEM%wSS !!!!!!
  !> if type_of_WW = 2 RECOMPUTATION from elem%wSS(1, 0:1, :) TO ARRAYS ELEM%wSS !!!!!!
  subroutine AdvancedInterpolDGsolution(gridN, gridO, type_of_WW )
    type(mesh), intent(inout), target	:: gridN
    type(mesh), intent(in), target	:: gridO
    integer, intent(in), optional :: type_of_WW
    class(element), pointer :: elem, elemO
    type(volume_rule), pointer :: V_rule
    integer :: Qdof, dof, max_tri,  dofO
    integer :: ie, ip, i, j, k, kst,  nsq, n1, n2, nmir
    integer :: il,lev, max_il, il_out, Tdeg1, iTd
    real, dimension(:,:), allocatable :: hxi, xi
    real, dimension(:,:, :), allocatable :: wi, qi, ww
    real, dimension(1:2, 1:2) :: rmx
    integer, dimension(:,:,:), allocatable :: itri
    integer, dimension(:,:), allocatable :: ntri
    real, dimension(:,:), allocatable :: phi ! local store arrays
    real, dimension(:,:), pointer :: phiL ! local store arrays
    real, dimension(:), pointer :: weights
    real, dimension(:,:,:), allocatable :: DF !, D1F
    real, dimension(:), allocatable :: JF
    real :: diff, norm, tol
    real :: x0(1:2), xS(1:2), xM(1:2)
    integer :: itest, type_of_WWWW

    type_of_WWWW = 0
    if( present(type_of_WW)) then
       type_of_WWWW = type_of_WW
    endif


    Tdeg1 = state%time%deg+1
    if( state%time%disc_time == 'STDG') Tdeg1 = 0

    !if(type_of_WWWW == 2) Tdeg1 = 1 ! original for p, p-1
    if(type_of_WWWW == 2) Tdeg1 = 2  ! new for p , p-1, p-2

    itest = -10

    if(state%isol == 5) itest = 321

    ! square of the tolerance for convergence of the adaptive numerical integration
    tol = 1E-12

    ! we create a Cartesian  grid covering the computational domain
    ! for each square cell a list of corresponding triangles is found

    call FindTriangleIndexs(gridO, nsq, rmx, max_tri, ntri, itri)

    max_il = 1
    if(state%time_dependent) max_il = 3

    ! over all new elements
    do ie=1, gridN%nelem
       elem => gridN%elem(ie)

       !if(elem%xc(1) > 0.2 .and. elem%xc(1) < 0.3 .and. &
       !     elem%xc(2) < -0.85) then
       !   print*,elem%xc(:), elem%i,'  $$$'
       !   !itest = ie
       !endif


       dof  = elem%dof
       Qdof = elem%Qdof
       V_rule => state%space%V_rule(elem%Qnum)

       allocate( weights(1:Qdof) )
       allocate( qi(0:Tdeg1, 1:dof, 1:ndim), ww(0:Tdeg1,1:dof*ndim,1:max_il) )
       allocate( wi(0:Tdeg1, 1:Qdof, 1:ndim) )

       ! hxi = barycentric coordinates, xi = physical coordinates
       allocate( hxi(1:Qdof, 1:nbDim), xi(1:Qdof, 1:nbDim))

       ! lev = level of refinement of subtriangles
       do il = 1,max_il
          lev = il**2

          qi(:, :, :) = 0.

          do n2 = 1, lev
             do n1 = 1, lev - n2 +1

                ! xS shift from the small reference subinterval to the actual one
                xS(1:2) = (/ 1. *(n1 -1)/ lev, 1. *(n2 -1)/ lev /)

                ! xM "mirror" point for the uper-right triangles
                xM(1:2) = xS(1:2) + (/ 0.5/lev, 0.5/lev  /)

                !!write(*,'(a6, 3i5,8es12.4)') 'ndx:',lev, n1, n2, xS(1:2), xM(1:2)
                do nmir =1,2
                   !!print*,'#######',n1, n2, lev - n2 +1
                   if(nmir == 1 .or. n1 < lev - n2 +1) then ! mirror of sub-triangle

                      !print*,'@@@@@@ A0b', Qdof, elem%Qnum

                      ! scale
                      hxi(1:Qdof,1:nbDim) = V_rule%lambda(1:Qdof,1:nbDim) / lev

                      ! shift
                      do i=1,Qdof
                         hxi(i,1:nbDim) = hxi(i,1:nbDim) + xS(1:nbDim)
                      enddo

                         ! Mirror
                      if(nmir == 2) then
                        do i=1,Qdof
                            hxi(i,1:nbDim) = - hxi(i,1:nbDim) + 2*xM(1:nbDim)
                         enddo
                      endif

                      !if(elem%i == itest) then
                      !   do i=1,Qdof
                      !      write(20+lev,*) hxi(i, :), lev, ie
                      !   enddo
                      !endif

                      !integration nodes on K
                      call ComputeF(elem, Qdof, hxi(1:Qdof, 1:nbDim), xi(1:Qdof, 1:nbDim) )

                      !for each integ node, we seek the square element
                      do ip=1,Qdof

                         ! we seek the old triangle elem0 containg node xi(xp, 1:2)
                         call FindTriangleCoords(gridO, elemO, xi(ip, 1:2), x0(1:2),&
                              nsq, max_tri, ntri, itri, rmx )

                         ! evaluation of the solution at the integ. node
                         dofO = elemO%dof
                         allocate(phi ( 1:dofO, 1:1) )

                         !test functions in x0
                         call Eval_phi_Qnode(elemO, dofO, 1, x0(1:nbDim), &
                              phi(1:dofO, 1:1) )

                         !if(elem%i == itest) then
                         !   if(lev == 1) print*,'@@@@', ip, elemO%i
                         !   write(100+lev,'(2i5,40es12.4)') ip, elemO%i, hxi(ip, :), xi(ip, 1:2), phi(1:dofO, 1:1)
                         !write(200+lev,'(2i5,40es12.4)') ip, elemO%i, hxi(ip, :), xi(ip, 1:2), elemO%w(0,:)
                         !endif


                         do k=1,ndim              ! k = index of component of w
                            kst = dofO*(k-1) + 1
                            do iTd = 0, Tdeg1
                               if(type_of_WWWW == 2) then
                                  wi(iTd,ip,k) = dot_product(elemO%wSS(1, iTd,kst:kst+dofO-1), &
                                       phi(1:dofO,1) )
                               else
                                  wi(iTd,ip,k) =dot_product(elemO%w(iTd,kst:kst+dofO-1), phi(1:dofO,1) )
                               endif
                            enddo
                         enddo

                         !if(elem%i == itest) write(30+lev,*) xi(ip, 1:2), wi(ip, :)

                         deallocate (phi)

                      enddo !! ip=1,Qdof

                      ! recomputation of the solution from the integ nodes to the
                      ! basis coefficients
                      allocate(phi ( 1:dof, 1:Qdof) )
                      call Eval_phi_Qnode(elem, dof, Qdof,  hxi(1:Qdof,1:nbDim), &
                           phi(1:dof, 1:Qdof) )

                      !!weights(1:Qdof) = V_rule%weights(1:Qdof)

                      if(elem%F%iFlin) then ! linear element
                         weights(1:Qdof) = V_rule%weights(1:Qdof) &
                              * elem%F%JF0 / (5-elem%type) /lev**2
                      else  ! curved element
                         allocate( DF(1:Qdof, 1:nbDim, 1:nbDim), JF(1:Qdof) )


                         call ComputeDF(elem, Qdof, hxi(1:Qdof, 1:nbDim), &
                              DF(1:Qdof, 1:nbDim, 1:nbDim) )

                         JF(1:Qdof) = DF(1:Qdof,1,1)*DF(1:Qdof,2,2) &
                              - DF(1:Qdof,1,2)*DF(1:Qdof,2,1)

                         weights(1:Qdof)  = V_rule%weights(1:Qdof) &
                              * JF(1:Qdof)/(5-elem%type)/lev**2

                         deallocate(DF, JF)

                         !write(*,'(a4,i5,12es12.4)') 'OLD',elem%i, &
                         !     V_rule%weights(1:12) * elem%area /lev**2
                         !
                         !write(*,'(a4,i5,12es12.4)') 'neW',elem%i, &
                         !     weights(1:12)
                         !
                         !print*

                      endif

                      ! evaluating of the coefficients of basis expansion
                      do k=1, ndim
                         ! integration of \int_K \phi_i w dx
                         do i=1,dof
                            do iTd = 0, Tdeg1
                               qi(iTd, i,k) = qi(iTd, i,k) +  & !!elem%F%JF0 /2. /lev**2 &
                                    !dot_product(V_rule%weights(1:Qdof), &
                                    dot_product(weights(1:Qdof), &
                                    phi( i, 1:Qdof) * wi(iTd, 1:Qdof,k) )
                            enddo

                         enddo
                      enddo ! do k

                      deallocate(phi)
                   endif ! mirror or the subelement
                end do  !do nmir = 1,2
             enddo ! do n1
          enddo    ! do n2


          do k=1, ndim
             do iTd = 0, Tdeg1
                ww(iTd, (k-1)*dof +1 : k*dof, il) &
                     = matmul(elem%MassInv%Mb(1:dof,1:dof), qi(iTd,1:dof, k) )
             enddo
          enddo ! do k


          !norm = 0.
          !do k=1,ndim
          if(il >= 2) then
             norm = dot_product( ww(0,:,il), ww(0,:, il) )
             diff = dot_product( ww(0,:,il)- ww(0,:,il -1), ww(0,:, il) -ww(0,:,il-1) )
             !write(*,'(a6,2i5,es18.10,4es12.4)') '@!!',ie, il, norm, diff
             if(norm <= 1E-35) goto 100
             if( diff / norm <= tol) then
                il_out = il
                goto 100
             endif

          endif
          il_out = il
       enddo  ! lev = 1, max_il

100    continue

       ! filling the projected solution
       if( type_of_WWWW == 1  ) then
          ! filling only elem%w(0,:) for AMA
          elem%wSS(1,1, 1:dof*ndim) = ww(0, 1:dof*ndim, il_out)

       elseif( type_of_WWWW == 2  ) then
          do iTd=0,Tdeg1
             elem%wSS(1, iTd, 1:dof*ndim) = ww(iTd, 1:dof*ndim, il_out)
          enddo

       else
          ! Full (standard) filling
          do iTd=0,Tdeg1
             elem%w(iTd, 1:dof*ndim) = ww(iTd, 1:dof*ndim, il_out)
          enddo


          if( state%time%disc_time == 'STDG') then
             do k=1,ndim
                elem%wSTfinAD(k,1:elem%dof) = elem%w(0, (k-1)*elem%dof + 1: k*elem%dof)
                elem%wST(k, 1:dof, 1) = elem%wSTfinAD(k,1:elem%dof)
             enddo
             elem%wST(1:ndim, 1:dof, 2:elem%Tdof) = 0
          endif
       endif

       deallocate( ww, weights )

       deallocate( hxi, xi, wi )
       deallocate( qi )

    enddo  ! ie = 1, gridN%nelem

  end subroutine AdvancedInterpolDGsolution





  !> higher order reconstruction via least square method in the \f$ H^1\f$-norm,
  !> \f$ P^{degP} \f$ approximation
  !> over elem and its neighbours given
  !> by \f$ \int_{D_e} \sum_j w_j \varphi_j \varphi_i dx =  \int_{D_e} w_h \varphi_i dx\f$,
  !> \f$ i=1,\dots dofP \f$
  !> USES elem%wS !!!
  subroutine LeastSquareInterpolationH1(elem, ndimL, MG, Qnum, degP, dofP, Dw, R_type  )
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    logical, intent(in) :: MG    ! .false. => use elem%dof, .true. => use elem%MGdof,
    integer, intent(in) :: Qnum  ! the prescribed quadrature
    integer, intent(in) :: degP, dofP  ! the desired degree of the reconstruction
    real, dimension(1:ndimL, 0:degP, 0:degP, 1:dofP), intent(inout) :: Dw
    integer, intent(in) :: R_type  ! = -1 support, 0 - all neighbours, 1- neigh 1,2, 2- neigh 2,3, ..
    class(element), pointer :: elem1
    type(volume_rule), pointer :: V_rule
    real, dimension(:,:), pointer :: phi0 ! the test functions
    real, dimension(:,:,:), allocatable :: Dphi0 ! the test functions
    real, dimension(:,:), allocatable :: phi ! the test functions
    real, dimension(:,:), allocatable :: phiW ! the test functions multiplied by weights
    real, dimension(:,:,:), allocatable :: Dphi ! Der of test functions
    real, dimension(:,:,:), allocatable :: DphiW ! Der of test functions multiplied by weights
    real, dimension(:,:), allocatable :: xi ! barycentric coordinates
    real, dimension(:,:), allocatable :: Fx ! real physical coordinates
    real, dimension(:,:), allocatable :: Mass ! mass matrix
    real, dimension(:,:), allocatable :: Stiff ! mass matrix
    real, dimension(:,:), allocatable :: rhs ! right-hand-side
    real, dimension(:,:), allocatable :: rhsSff ! right-hand-side
    real, dimension(:), allocatable :: weights, wi
    real, dimension(:,:), allocatable :: Dwi
    integer :: dof, Qdof,  Qnum1, Qdof1, dof1, N_elem
    integer :: i, i1, j, k, l, l1, l2, n, n1

    ! Qdiff is the increase of the order of quadrature rule used for the reconstruction
    dof = elem%dof
    if(MG)  dof = elem%MGdof

    if(R_type == 10 .or. R_type == -10) dof = DOFtriang( max(0, elem%deg - 1) )
    if(R_type == 20 .or. R_type == -20) dof = DOFtriang( max(0, elem%deg - 2) )
    if(R_type == 50 .or. R_type == -50) dof = DOFtriang( max(0, elem%deg + 1) )

    !!!!!!Qnum = elem%Qnum + Qdiff
    Qdof = state%space%V_rule(Qnum)%Qdof


    ! Mass(i,j) = \int_{K and its neighbours} \phi_i \phi_j dx
    ! rhs(i) = \int_{K and its neighbours} \phi_i w_h dx
    allocate(Mass(1:dofP, 1:dofP), Stiff(1:dofP, 1:dofP))
    allocate(rhs(1:ndimL, 1:dofP), rhsSff(1:ndimL, 1:dofP) )

    Mass(:,:) = 0.
    Stiff(:,:) = 0.
    rhs(:,:) = 0.
    rhsSff(:,:) = 0.

    !integerals over elem
    allocate(weights(1:Qdof) )
    weights(:) = 1.
    !call IntegrateBlockBB(elem, dofP, weights(1:Qdof), Mass(1:dofP,1:dofP) )
    call IntegrateBlockBBplus(elem, Qnum, dofP, weights(1:Qdof), Mass(1:dofP,1:dofP))

    call IntegrateBlockD2plus(elem, Qnum, dofP, weights(1:Qdof),Stiff(1:dofP,1:dofP))

    deallocate(weights )

    phi0 => state%space%V_rule(Qnum)%phi(1:dofP, 1:Qdof)

    allocate( Dphi0(1:dofP, 1:nbDim, 1:Qdof) )
    call Eval_Dphi_plus(elem, state%space%V_rule(Qnum), dofP, Dphi0(1:dofP, 1:nbDim, 1:Qdof) )


    ! SMAZ - begin
    !allocate(Fx(1:Qdof, 1:2) )
    !call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )
    !do l=1,dofP
    !   do l1 =1,Qdof
    !      write(100+l,*) Fx(l1, 1:2), phi0(l, l1), Dphi0(l, :, l1), 'C',elem%i,l
    !   enddo
    !enddo
    !deallocate(Fx)
    ! SMAZ -end

    allocate( wi(1:Qdof) , Dwi(1:Qdof, 1:nbDim) )
    do k=1, ndimL
       ! product (wi, phi0)
       wi(1:Qdof) = matmul(elem%wS(k, 1 : dof), phi0(1:dof, 1:Qdof) )
       !!WS!! wi(1:Qdof) = matmul(elem%w(0, (k-1)*dof +1 : k*dof), phi0(1:dof, 1:Qdof) )

       !if(elem%i == 2) then
       !   write(*,'(a6, 2i5,500es12.4)') 'wS aa', elem%i, dof, elem%wS(1, 1 : dof)
       !endif


       !call IntegrateVectorB(elem, dofP, wi(1:Qdof), rhs(k, 1:dofP) )
       call IntegrateVectorBplus(elem, Qnum, dofP, wi(1:Qdof), rhs(k, 1:dofP) )

       ! product (Dwi, Dphi0)
       do n=1,nbDim
          !!WS!! Dwi(1:Qdof, n) = matmul(elem%w(0, (k-1)*dof +1 : k*dof), &
          Dwi(1:Qdof, n) = matmul(elem%wS(k,1 : dof), Dphi0(1:dof, n, 1:Qdof))
       enddo
       call IntegrateVectorDplus(elem, Qnum, dofP, Dwi(1:Qdof, 1:nbDim), &
            rhsSff(k, 1:dofP))

       ! write(*,'(a6,200es12.4)') 'wS:', elem%wS(k, 1 : dof)
       ! write(*,'(a6,200es12.4)') 'wi:', wi(1:Qdof)
       ! write(*,'(a6,200es12.4)') 'Dwi:', Dwi(1:Qdof,1)
       ! write(*,'(a6,200es12.4)') 'Dwi:', Dwi(1:Qdof,2)
       ! write(*,'(a6,200es12.4)') 'rhs:', rhs(k, 1:dofP)
       ! write(*,'(a6,200es12.4)') 'rhsSff:', rhsSff(k, 1:dofP)

    enddo


    !allocate(Fx(1:Qdof, 1:2) )
    !call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )
    !do l=1, Qdof
    !   write(91,'(30es12.4)') Fx(l, 1:2),  state%space%V_rule(Qnum)%lambda(l, 1:2), phi0(1:dofP, l)
    !   write(93,'(30es12.4)') Fx(l, 1:2),  state%space%V_rule(Qnum)%lambda(l, 1:2), wi(l)
    !enddo
    !deallocate(Fx)

    deallocate(wi, Dwi, Dphi0 )

    !!!goto 111

    ! adding of integrals over neighbouring elements
    ! phi0  basis functions from elem1 in elem1's integ nodes
    ! phi   basis functions from elem  in elem1's integ nodes, i.e., outside of elem


    ! type of the patch
    if(     R_type == 0 .or. R_type ==  10 .or. R_type ==  20 .or. R_type ==  50) then
       N_elem = elem%flen   ! all neighbours

    elseif(R_type == -1 .or. R_type == -10 .or. R_type == -20 .or. R_type == -50) then
       N_elem = elem%isupp   ! elements sharing a vertex

    elseif(R_type ==  1) then
       N_elem = 2            ! first and second neighbours
    elseif(R_type ==  2) then
       N_elem = 2            ! second and third neighbours
    elseif(R_type ==  3) then
       N_elem = 2            ! third and first neighbours

    else
       stop 'UNKNOWN R_type in LeastSquareInterpolationH1'
    endif

    do j=1, N_elem
       if(R_type ==  0 .or. R_type ==  10 .or. R_type ==  20 .or. R_type ==  50) i1 = elem%face(neigh,j)
       if(R_type == -1 .or. R_type == -10 .or. R_type == -20 .or. R_type == -50) i1 = elem%supp(j, 1)

       if(R_type ==  1) i1 = elem%face(neigh,j)
       if(R_type ==  2) i1 = elem%face(neigh,j+1)
       if(R_type ==  3) i1 = elem%face(neigh,mod(j+1,3)+1)


       if(i1 > 0) then
          elem1 => grid%elem(i1)
          dof1  = elem1%dof
          if(MG)  dof1 = elem1%MGdof

          if(R_type == 10 .or. R_type == -10) dof1 = DOFtriang( max(0, elem1%deg - 1) )
          if(R_type == 20 .or. R_type == -20) dof1 = DOFtriang( max(0, elem1%deg  -2) )
          if(R_type == 50 .or. R_type == -50) dof1 = DOFtriang( max(0, elem1%deg + 1) )

          Qnum1 = Qnum  !! probably more accurate than the previous one
          Qdof1 = state%space%V_rule(Qnum1)%Qdof

          !Qdof1 = elem1%Qdof

          ! Fx integ nodes on elem1 - real physical cordinates
          ! xi barycentric coordinates of Fx (integ nodes on elem1) with respect elem
          allocate(Fx(1:Qdof1, 1:2), xi(1:Qdof1, 1:2) )
          call ComputeF(elem1, Qdof1, state%space%V_rule(Qnum1)%lambda(1:Qdof1, 1:2), &
               Fx(1:Qdof1, 1:2) )
          call BarycCoord(elem, Qdof1, Fx(1:Qdof1, 1:2),  xi(1:Qdof1, 1:2) )

          allocate(phi(1:dofP, 1:Qdof1) ) ! value of the test function in integ nodes
          allocate(Dphi(1:dofP, 1:nbDim, 1:Qdof1) ) ! value of Deriv test functions


          call Eval_phi_Qnode(elem, dofP, Qdof1, xi(1:Qdof1, 1:nbDim), &
               phi(1:dofP, 1:Qdof1), Dphi(1:dofP, 1:nbDim, 1:Qdof1) )


          allocate(wi(1:Qdof1),  weights(1:Qdof1) )
          allocate(Dwi(1:Qdof1, 1:nbDim)  )
          call Eval_V_Weights_plus(elem1, state%space%V_rule(Qnum1), weights(1:Qdof1))

          allocate(phiW(1:dofP, 1:Qdof1) )           ! phi multiplied by weights
          allocate(DphiW(1:dofP, 1:nbDim, 1:Qdof1) ) ! Dphi multiplied by weights

          do l=1, dofP
             phiW(l, 1:Qdof1) = phi(l, 1:Qdof1) * weights(1:Qdof1)

             do n=1,nbDim
                DphiW(l, n, 1:Qdof1) = Dphi(l, n, 1:Qdof1) * weights(1:Qdof1)
             enddo

             !do l1 =1,Qdof
             !   write(100+l,*) Fx(l1, 1:2), phi(l, l1), Dphi(l, :, l1), xi(l1, :), j, i1
             !write(100+10*j + l,*) &
             !        Fx(l1, 1:2), phi(l, l1), Dphi(l, :, l1), xi(l1, :), j, i1
             !enddo
          enddo

          !if(elem%i == 1) then
             !write(*,'(a6,3i5,50es12.4)') &
             !     'N_a',elem%i, degP, Qdof1, Mass(1,2),  & !  phi(1, 1:40)
             !     Mass(1, 1:15)
             !return
          !endif


          ! adding of Mass and Stiff  matrices
          do l1 = 1, dofP
             do l2 = 1, dofP
                Mass(l1, l2) = Mass(l1, l2) &
                     + dot_product(phiW(l1, 1:Qdof1), phi(l2, 1:Qdof1) )

                do n=1,nbDim
                   Stiff(l1, l2) = Stiff(l1, l2) &
                        + dot_product(DphiW(l1, n, 1:Qdof1), Dphi(l2, n, 1:Qdof1) )
                enddo
             enddo
          enddo

          !if(elem%i == 1) then
          !   write(*,'(a6,3i5,50es12.4)') &
          !        'N_a',elem1%i, degP, Qdof1, Mass(1,2), & !Mass(1, 1:15)
          !        phi(2, 1:40)
          !
          !   return
          !endif


          ! adding of rhs
          phi0 => state%space%V_rule(Qnum1)%phi(1:dof1, 1:Qdof1)

          allocate( Dphi0(1:dof1, 1:nbDim, 1:Qdof1) )
          call Eval_Dphi_plus(elem1, state%space%V_rule(Qnum1), dof1, Dphi0(1:dof1, 1:nbDim, 1:Qdof1) )

          !if(elem%i == 2) then
          !    write(*,'(a6, 2i5,500es12.4)') 'wS aa', i1, dof1, elem1%wS(1, 1 : dof1)
          !endif


          do k=1, ndimL
             ! wi = values of w at integ nodes of elem1
             !!WS!! wi(1:Qdof1) = matmul(elem1%w(0, (k-1)*dof1 +1 : k*dof1), &
             wi(1:Qdof1) = matmul(elem1%wS(k, 1 : dof1),  phi0(1:dof1, 1:Qdof1) )
             ! product (wi, phi0)
             do l=1, dofP
                rhs(k, l) = rhs(k, l)  &
                      + dot_product(wi(1:Qdof1), phiW(l, 1:Qdof1) )
             enddo


             ! Dwi in integ nodes
             do n=1,nbDim
                !!WS!! Dwi(1:Qdof1, n) = matmul(elem1%w(0, (k-1)*dof1 +1 : k*dof1), &
                Dwi(1:Qdof1, n) = matmul(elem1%wS(k, 1 : dof1), Dphi0(1:dof1, n, 1:Qdof1))
             enddo

             ! product (Dwi, Dphi0)
             do l=1, dofP
                do n=1,nbDim
                   rhsSff(k, l) = rhsSff(k, l)  &
                        + dot_product(Dwi(1:Qdof1, n), DphiW(l, n, 1:Qdof1) )
                enddo

             enddo

       ! write(*,'(a6,200es12.4)') 'wS:', elem1%wS(k, 1 : dof1)
       ! write(*,'(a6,200es12.4)') 'wi:', wi(1:Qdof1)
       ! write(*,'(a6,200es12.4)') 'Dwi:', Dwi(1:Qdof1,1)
       ! write(*,'(a6,200es12.4)') 'Dwi:', Dwi(1:Qdof1,2)
       ! write(*,'(a6,200es12.4)') 'rhs:', rhs(k, 1:dofP)
       ! write(*,'(a6,200es12.4)') 'rhsSff:', rhsSff(k, 1:dofP)

          enddo

          !do l=1, Qdof1
          !   write(91,'(30es12.4)') Fx(l, 1:2),  xi(l, 1:2), phi(1:dofP, l)
          !   write(92,'(30es12.4)') Fx(l, 1:2),  xi(l, 1:2), wi(l)
          !enddo

          deallocate(Fx, xi, phi, phiW, wi, weights)
          deallocate(Dphi, DphiW, Dwi, Dphi0)

       endif  !(i1 > 0)

    enddo ! j=1, N_elem

!!!111 continue

    !do i=1, dofP
    !   write(*,'(a6,i5,30es12.4)') 'mass',i, Mass(i, 1:dofP)
    !enddo
    !do i=1, dofP
    !   write(33,'(a6,i5,30es12.4)') 'mass',i, Mass(i, :),  rhs(1, i)
    !enddo
    !write(33,*) '####', state%space%adapt%adapt_level, elem%i, degP,'#####################'

    ! the H1-norm (comment for the L^2-norm)
    Mass(1:dofP, 1:dofP) = Mass(1:dofP, 1:dofP) + Stiff(1:dofP, 1:dofP)
    rhs(1:ndimL, 1:dofP) = rhs(1:ndimL, 1:dofP) + rhsSff(1:ndimL, 1:dofP)


    ! evaluation of the coefficients of the reconstruction
    call MblockInverse(dofP, Mass(1:dofP, 1:dofP) )

    !do i=1, dofP
    !   write(*,'(a6,i5,30es12.4)') 'mass',i, Mass(i, 1:dofP)
    !enddo

    do k=1, ndimL
       ! zero order derivative = a reconstruction of  function w
       Dw(k, 0, 0, 1:dofP)  = matmul(Mass(1:dofP, 1:dofP), rhs(k, 1:dofP) )

       !! only for test values of reconstruction in integ nodes
       !if(elem%i == 289) print*,'#############',elem%i, k
       !if(elem%xc(1) > 0.95  .and. elem%xc(2) < 0.4 .and. elem%xc(2)> 0.2 ) then
       !phi0 => state%space%V_rule(Qnum)%phi(1:dofP, 1:Qdof)
       !allocate(Fx(1:Qdof, 1:2) )
       !call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )
       !do l=1,Qdof
       !   write(80+state%space%adapt%adapt_level,'(30es12.4)') &
       !        Fx(l, 1:2), dot_product(Dw(k,0, 0, 1:dofP), phi0(1:dofP, l) )
       !enddo
       !deallocate(Fx)
       !endif

    enddo

    deallocate( Mass, Stiff, rhs, rhsSff )

    !write(*,'(a6,i5,30es12.4)') '## sol',0,Dw(1, 0, 0, 1:dofP)
    !if(elem%i == 1) stop

  end subroutine LeastSquareInterpolationH1


  !> higher order reconstruction via least square method in the \f$ L^2\f$-norm,
  !> \f$ P^{degP} \f$ approximation
  !> over elem and its neighbours given
  !> by \f$ \int_{D_e} \sum_j w_j \varphi_j \varphi_i dx =  \int_{D_e} w_h \varphi_i dx\f$,
  !> \f$ i=1,\dots dofP \f$
  subroutine LeastSquareInterpolationL2(elem, ndimL, MG, Qnum, degP, dofP, Dw, R_type  )
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    logical, intent(in) :: MG    ! .false. => use elem%dof, .true. => use elem%MGdof,
    integer, intent(in) :: Qnum  ! the prescribed quadrature
    integer, intent(in) :: degP, dofP  ! the desired degree of the reconstruction
    real, dimension(1:ndimL, 0:degP, 0:degP, 1:dofP), intent(inout) :: Dw
    integer, intent(in) :: R_type
    class(element), pointer :: elem1
    type(volume_rule), pointer :: V_rule
    real, dimension(:,:), pointer :: phi0 ! the test functions
    real, dimension(:,:), allocatable :: phi ! the test functions
    real, dimension(:,:), allocatable :: phiW ! the test functions multiplied by weights
    real, dimension(:,:), allocatable :: xi ! barycentric coordinates
    real, dimension(:,:), allocatable :: Fx ! real physical coordinates
    real, dimension(:,:), allocatable :: Mass ! mass matrix
    real, dimension(:,:), allocatable :: rhs ! right-hand-side
    real, dimension(:), allocatable :: weights, wi
    integer :: dof, Qdof, Qnum1, Qdof1, dof1, N_elem
    integer :: i, i1, j, k, l, l1, l2

    ! Qdiff is the increase of the order of quadrature rule used for the reconstruction
    dof = elem%dof
    if(MG)  dof = elem%MGdof

    !!!!!!Qnum = elem%Qnum + Qdiff
    Qdof = state%space%V_rule(Qnum)%Qdof

    ! Mass(i,j) = \int_{K and its neighbours} \phi_i \phi_j dx
    ! rhs(i) = \int_{K and its neighbours} \phi_i w_h dx
    allocate(Mass(1:dofP, 1:dofP), rhs(1:ndimL, 1:dofP) )

    Mass(:,:) = 0.
    rhs(:,:) = 0.

    allocate(weights(1:Qdof) )
    weights(:) = 1.
    !call IntegrateBlockBB(elem, dofP, weights(1:Qdof), Mass(1:dofP,1:dofP) )
    call IntegrateBlockBBplus(elem, Qnum, dofP, weights(1:Qdof), Mass(1:dofP,1:dofP) )

    deallocate(weights )

    phi0 => state%space%V_rule(Qnum)%phi(1:dofP, 1:Qdof)

    allocate( wi(1:Qdof) )
    do k=1, ndimL
       !!WS!! wi(1:Qdof) = matmul(elem%w(0, (k-1)*dof +1 : k*dof), phi0(1:dof, 1:Qdof) )
       wi(1:Qdof) = matmul(elem%wS(k, 1 : dof), phi0(1:dof, 1:Qdof) )
       !call IntegrateVectorB(elem, dofP, wi(1:Qdof), rhs(k, 1:dofP) )
       call IntegrateVectorBplus(elem, Qnum, dofP, wi(1:Qdof), rhs(k, 1:dofP) )
    enddo

    !allocate(Fx(1:Qdof, 1:2) )
    !call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )
    !do l=1, Qdof
    !   write(91,'(30es12.4)') Fx(l, 1:2),  state%space%V_rule(Qnum)%lambda(l, 1:2), phi0(1:dofP, l)
    !   write(93,'(30es12.4)') Fx(l, 1:2),  state%space%V_rule(Qnum)%lambda(l, 1:2), wi(l)
    !enddo
    !deallocate(Fx)

    deallocate(wi)


    ! adding of integrals over neighbouring elements
    ! phi0  basis functions from elem1 in elem1's integ nodes
    ! phi   basis functions from elem  in elem1's integ nodes, i.e., outside of elem

     ! type of the patch
    if(R_type == 0) then
       N_elem = elem%flen   ! all neighbours
    elseif(R_type == -1) then
       N_elem = elem%isupp   ! elements sharing a vertex
    elseif(R_type ==  1) then
       N_elem = 2            ! first and second neighbours
    elseif(R_type ==  2) then
       N_elem = 2            ! second and third neighbours
    elseif(R_type ==  3) then
       N_elem = 2            ! thirs and first neighbours
    else
       stop 'UNKNOWN R_type in LeastSquareInterpolationH1'
    endif

    do j=1, N_elem
       if(R_type == 0)  i1 = elem%face(neigh,j)
       if(R_type == -1) i1 = elem%supp(j, 1)
       if(R_type ==  1) i1 = elem%face(neigh,j)
       if(R_type ==  2) i1 = elem%face(neigh,j+1)
       if(R_type ==  3) i1 = elem%face(neigh,mod(j+1,3)+1)

       if(i1 > 0) then
          elem1 => grid%elem(i1)
          dof1  = elem1%dof
          if(MG)  dof1 = elem1%MGdof

          Qnum1 = Qnum  !! probably more accurate than the previous one
          Qdof1 = state%space%V_rule(Qnum1)%Qdof

          ! Fx integ nodes on elem1 - real physical cordinates
          ! xi barycentric coordinates of Fx (integ nodes on elem1) with respect elem
          allocate(Fx(1:Qdof1, 1:2), xi(1:Qdof1, 1:2) )
          call ComputeF(elem1, Qdof1, state%space%V_rule(Qnum1)%lambda(1:Qdof1, 1:2), Fx(1:Qdof1, 1:2) )
          call BarycCoord(elem, Qdof1, Fx(1:Qdof1, 1:2),  xi(1:Qdof1, 1:2) )

          allocate(phi(1:dofP, 1:Qdof1) ) ! value of the test function in integ nodes
          call PHI_orthonormal(Qdof1, nbDim, xi(1:Qdof1, 1:nbDim), 3, dofP, phi(1:dofP, 1:Qdof1) )

          allocate(wi(1:Qdof1), weights(1:Qdof1) )
          call Eval_V_Weights_plus(elem1, state%space%V_rule(Qnum1), weights)

          allocate(phiW(1:dofP, 1:Qdof1) ) ! phi multiplied by weights
          do l=1, dofP
             phiW(l, 1:Qdof1) = phi(l, 1:Qdof1) * weights(1:Qdof1)
          enddo

          ! adding of Mass matrix
          do l1 = 1, dofP
             do l2 = 1, dofP
                Mass(l1, l2) = Mass(l1, l2) &
                     + dot_product(phiW(l1, 1:Qdof1), phi(l2, 1:Qdof1) )
             enddo
          enddo

          ! adding of rhs
          phi0 => state%space%V_rule(Qnum1)%phi(1:dof1, 1:Qdof1)

          do k=1, ndimL
             ! wi = values of w at integ nodes of elem1
             !!WS!! wi(1:Qdof1) = matmul(elem1%w(0, (k-1)*dof1 +1 : k*dof1), phi0(1:dof1, 1:Qdof1) )
             wi(1:Qdof1) = matmul(elem1%wS(k, 1:dof1), phi0(1:dof1, 1:Qdof1) )

             do l=1, dofP
                rhs(k, l) = rhs(k, l) + dot_product(wi(1:Qdof1), phiW(l, 1:Qdof1) )
             enddo
          enddo


          !do l=1, Qdof1
          !   write(91,'(30es12.4)') Fx(l, 1:2),  xi(l, 1:2), phi(1:dofP, l)
          !   write(92,'(30es12.4)') Fx(l, 1:2),  xi(l, 1:2), wi(l)
          !enddo

          deallocate(Fx, xi, phi, phiW, wi, weights)
       endif  !(i1 > 0)
    enddo ! j=1, N_elem


    !do i=1, dofP
    !   write(*,'(a6,i5,30es12.4)') 'mass',i, Mass(i, 1:dofP)
    !enddo
    !do i=1, dofP
    !   write(34,'(a6,i5,30es12.4)') 'mass',i, Mass(i, :),  rhs(1, i)
    !enddo
    !write(34,*) '####', state%space%adapt%adapt_level, elem%i, degP,'#####################'




    ! evaluation of the coefficients of the reconstruction
    call MblockInverse(dofP, Mass(1:dofP, 1:dofP) )

    !do i=1, dofP
    !   write(*,'(a6,i5,30es12.4)') 'mass',i, Mass(i, 1:dofP)
    !enddo

    do k=1, ndimL
       ! zero order derivative = a reconstruction of  function w
       Dw(k, 0, 0, 1:dofP)  = matmul(Mass(1:dofP, 1:dofP), rhs(k, 1:dofP) )

       !! only for test values of reconstruction in integ nodes
       !if(elem%i == 289) print*,'#############',elem%i, k
       !if(elem%xc(1) > 0.95  .and. elem%xc(2) < 0.4 .and. elem%xc(2)> 0.2 ) then
          phi0 => state%space%V_rule(Qnum)%phi(1:dofP, 1:Qdof)
          allocate(Fx(1:Qdof, 1:2) )
          call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )
          !do l=1,Qdof
          !   write(80+state%space%adapt%adapt_level,'(30es12.4)') &
          !        Fx(l, 1:2), dot_product(Dw(k,0, 0, 1:dofP), phi0(1:dofP, l) )
          !enddo
          deallocate(Fx)
       !endif

    enddo

    !write(*,'(a6,i5,30es12.4)') '## sol',0,Dw(1, 0, 0, 1:dofP)

    !stop
    deallocate( Mass, rhs )

  end subroutine LeastSquareInterpolationL2




  !> higher order reconstruction via least square method in the \f$ H^1\f$-norm,
  !> \f$ P^{degP} \f$ approximation
  !> over elem and its neighbours given
  !> by \f$ \int_{D_e} \sum_j w_j \varphi_j \varphi_i dx =  \int_{D_e} w_h \varphi_i dx\f$,
  !> \f$ i=1,\dots dofP \f$
  !> USES elem%wS !!!
  !> works only if grid%ElemSupports is already done
  subroutine LeastSquareInterpolationWeighted(elem, ndimL, MG, Qnum, degP, dofP, Dw, R_type  )
    class(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    logical, intent(in) :: MG    ! .false. => use elem%dof, .true. => use elem%MGdof,
    integer, intent(in) :: Qnum  ! the prescribed quadrature
    integer, intent(in) :: degP, dofP  ! the desired degree of the reconstruction
    real, dimension(1:ndimL, 0:degP, 0:degP, 1:dofP), intent(inout) :: Dw
    integer, intent(in) :: R_type  ! = -1 H1-norm, 0 - L2-norm
    class(element), pointer :: elem1
    type(volume_rule), pointer :: V_rule
    real, dimension(:,:), pointer :: phi0 ! the test functions
    real, dimension(:,:,:), allocatable :: Dphi0 ! the test functions
    real, dimension(:,:), allocatable :: phi ! the test functions
    real, dimension(:,:), allocatable :: phiW ! the test functions multiplied by weights
    real, dimension(:,:,:), allocatable :: Dphi ! Der of test functions
    real, dimension(:,:,:), allocatable :: DphiW ! Der of test functions multiplied by weights
    real, dimension(:,:), allocatable :: xi ! barycentric coordinates
    real, dimension(:,:), allocatable :: Fx ! real physical coordinates
    real, dimension(:,:), allocatable :: Mass ! mass matrix
    real, dimension(:,:), allocatable :: Stiff ! mass matrix
    real, dimension(:,:), allocatable :: rhs ! right-hand-side
    real, dimension(:,:), allocatable :: rhsSff ! right-hand-side
    real, dimension(:), allocatable :: weights, wi
    real, dimension(:,:), allocatable :: Dwi
    integer :: dof, Qdof,  Qnum1, Qdof1, dof1, N_elem
    integer :: i, i1, j, k, l, l1, l2, n, n1
    integer:: p_min, p_max, p_diff
    real :: weight, L2weight


    if(degP > MaxDegreeImplemented ) stop "Too high degP in LeastSquareInterpolationWeighted"

    ! Qdiff is the increase of the order of quadrature rule used for the reconstruction
    dof = elem%dof
    if(MG)  dof = elem%MGdof

    !if(abs(R_type) == 10 ) dof = DOFtriang( max(0, elem%deg - 1) )
    !if(abs(R_type) == 20 ) dof = DOFtriang( max(0, elem%deg - 2) )
    !if(abs(R_type) == 50 ) dof = DOFtriang( max(0, elem%deg + 1) )

    !!!!!!Qnum = elem%Qnum + Qdiff
    Qdof = state%space%V_rule(Qnum)%Qdof


    ! Mass(i,j) = \int_{K and its neighbours} \phi_i \phi_j dx
    ! rhs(i) = \int_{K and its neighbours} \phi_i w_h dx
    allocate(Mass(1:dofP, 1:dofP), Stiff(1:dofP, 1:dofP))
    allocate(rhs(1:dofP, 1:ndimL), rhsSff(1:dofP, 1:ndimL) )

    Mass(:,:) = 0.
    Stiff(:,:) = 0.
    rhs(:,:) = 0.
    rhsSff(:,:) = 0.

    !integerals over elem
    allocate(weights(1:Qdof) )
    weights(:) = 1.
    call IntegrateBlockBBplus(elem, Qnum, dofP, weights(1:Qdof), Mass(1:dofP,1:dofP))

    call IntegrateBlockD2plus(elem, Qnum, dofP, weights(1:Qdof),Stiff(1:dofP,1:dofP))

    deallocate(weights )

    phi0 => state%space%V_rule(Qnum)%phi(1:dofP, 1:Qdof)

    allocate( Dphi0(1:dofP, 1:nbDim, 1:Qdof) )
    call Eval_Dphi_plus(elem, state%space%V_rule(Qnum), dofP, Dphi0(1:dofP, 1:nbDim, 1:Qdof) )



    allocate( wi(1:Qdof) , Dwi(1:Qdof, 1:nbDim) )
    do k=1, ndimL
       ! product (wi, phi0)
       wi(1:Qdof) = matmul(elem%wS(k, 1 : dof), phi0(1:dof, 1:Qdof) )
       !!WS!! wi(1:Qdof) = matmul(elem%w(0, (k-1)*dof +1 : k*dof), phi0(1:dof, 1:Qdof) )

       !if(elem%i == 1199) write(*,'(a6,2i5, 30es12.4)') 'de45^',elem%i, dof, wi(:)
       !call IntegrateVectorB(elem, dofP, wi(1:Qdof), rhs(1:dofP, k) )
       call IntegrateVectorBplus(elem, Qnum, dofP, wi(1:Qdof), rhs(1:dofP, k) )

       ! product (Dwi, Dphi0)
       do n=1,nbDim
          !!WS!! Dwi(1:Qdof, n) = matmul(elem%w(0, (k-1)*dof +1 : k*dof), &
          Dwi(1:Qdof, n) = matmul(elem%wS(k,1 : dof), Dphi0(1:dof, n, 1:Qdof))
       enddo
       call IntegrateVectorDplus(elem, Qnum, dofP, Dwi(1:Qdof, 1:nbDim), &
            rhsSff(1:dofP, k))

    enddo


    deallocate(wi, Dwi, Dphi0 )

    !!!goto 111

    ! seeking of the variantion of elem%deg for the neighbours
    p_min = elem%deg
    p_max = elem%deg

    !do j=1,elem%flen
    !   i1 = elem%face(neigh, j)

    N_elem = elem%isupp
    do j=1, N_elem
       i1 = elem%supp(j, 1)
       if(i1 > 0) then
          p_min = min(p_min, grid%elem(i1)%deg)
          p_max = max(p_max, grid%elem(i1)%deg)
       endif
    enddo
    p_diff = p_max - p_min

    ! adding of integrals over neighbouring elements
    ! phi0  basis functions from elem1 in elem1's integ nodes
    ! phi   basis functions from elem  in elem1's integ nodes, i.e., outside of elem


    N_elem = elem%isupp   ! elements sharing a vertex

    do j=1, N_elem

       i1 = elem%supp(j, 1)

       if(i1 > 0) then
          if(elem%supp(j, 2) == 1) weight = 1.   ! do dot comment this line

          !if(elem%supp(j, 2) == 2) weight = 0.5   ! interpol estimates
          if(elem%supp(j, 2) == 2) weight = 1.0  !0.5  !0.5   !0.05  ! 0.005

          if( state%space%estim_space == 'HO_rec') then
             if(elem%supp(j, 2) == 2) weight = 0.02 * (p_diff -1)  ! AMC 2015, HO_reconstruction
          endif


          if(elem%deg == 0) weight = 1.

          if(weight > 0.) then
             !if(elem%supp(j,2)==2 .and. p_diff > 0) &
             !     print*,'Weight-',elem%i,j,i1, elem%supp(j,2) ,p_diff , weight

             elem1 => grid%elem(i1)
             dof1  = elem1%dof
             if(MG)  dof1 = elem1%MGdof

             !if( abs(R_type) == 10 ) dof1 = DOFtriang( max(0, elem1%deg - 1) )
             !if( abs(R_type) == 20 ) dof1 = DOFtriang( max(0, elem1%deg  -2) )
             !if( abs(R_type) == 50 ) dof1 = DOFtriang( max(0, elem1%deg + 1) )

             Qnum1 = Qnum  !! probably more accurate than the previous one
             Qdof1 = state%space%V_rule(Qnum1)%Qdof

             !Qdof1 = elem1%Qdof

             ! Fx integ nodes on elem1 - real physical cordinates
             ! xi barycentric coordinates of Fx (integ nodes on elem1) with respect elem
             allocate(Fx(1:Qdof1, 1:2), xi(1:Qdof1, 1:2) )
             call ComputeF(elem1, Qdof1, state%space%V_rule(Qnum1)%lambda(1:Qdof1, 1:2), &
                  Fx(1:Qdof1, 1:2) )
             call BarycCoord(elem, Qdof1, Fx(1:Qdof1, 1:2),  xi(1:Qdof1, 1:2) )

             allocate(phi(1:dofP, 1:Qdof1) ) ! value of the test function in integ nodes
             allocate(Dphi(1:dofP, 1:nbDim, 1:Qdof1) ) ! value of Deriv test functions


             call Eval_phi_Qnode(elem, dofP, Qdof1, xi(1:Qdof1, 1:nbDim), &
                  phi(1:dofP, 1:Qdof1), Dphi(1:dofP, 1:nbDim, 1:Qdof1) )


             allocate(wi(1:Qdof1),  weights(1:Qdof1) )
             allocate(Dwi(1:Qdof1, 1:nbDim)  )
             call Eval_V_Weights_plus(elem1, state%space%V_rule(Qnum1), weights(1:Qdof1))

             allocate(phiW(1:dofP, 1:Qdof1) )           ! phi multiplied by weights
             allocate(DphiW(1:dofP, 1:nbDim, 1:Qdof1) ) ! Dphi multiplied by weights

             do l=1, dofP
                phiW(l, 1:Qdof1) = phi(l, 1:Qdof1) * weights(1:Qdof1)

                do n=1,nbDim
                   DphiW(l, n, 1:Qdof1) = Dphi(l, n, 1:Qdof1) * weights(1:Qdof1)
                enddo
             enddo

             ! adding of Mass and Stiff  matrices
             do l1 = 1, dofP
                do l2 = 1, dofP
                   Mass(l1, l2) = Mass(l1, l2) &
                        + weight * dot_product(phiW(l1, 1:Qdof1), phi(l2, 1:Qdof1) )

                   do n=1,nbDim
                      Stiff(l1, l2) = Stiff(l1, l2) &
                           + weight * dot_product(DphiW(l1, n, 1:Qdof1), Dphi(l2, n, 1:Qdof1) )
                   enddo
                enddo
             enddo

             ! adding of rhs
             phi0 => state%space%V_rule(Qnum1)%phi(1:dof1, 1:Qdof1)

             allocate( Dphi0(1:dof1, 1:nbDim, 1:Qdof1) )
             call Eval_Dphi_plus(elem1, state%space%V_rule(Qnum1), dof1, Dphi0(1:dof1, 1:nbDim, 1:Qdof1) )


             do k=1, ndimL
                ! wi = values of w at integ nodes of elem1
                !!WS!! wi(1:Qdof1) = matmul(elem1%w(0, (k-1)*dof1 +1 : k*dof1), &
                wi(1:Qdof1) = matmul(elem1%wS(k, 1 : dof1),  phi0(1:dof1, 1:Qdof1) )
                !write(*,'(a8, 2i5, 20es12.4)') 'de3:',elem%i,k, wi(1), (elem1%xc(1)/40)**2, elem1%xc
                ! product (wi, phi0)
                do l=1, dofP
                   rhs(l, k) = rhs(l, k)  &
                        + weight * dot_product(wi(1:Qdof1), phiW(l, 1:Qdof1) )
                enddo


                ! Dwi in integ nodes
                do n=1,nbDim
                   !!WS!! Dwi(1:Qdof1, n) = matmul(elem1%w(0, (k-1)*dof1 +1 : k*dof1), &
                   Dwi(1:Qdof1, n) = matmul(elem1%wS(k, 1 : dof1), Dphi0(1:dof1, n, 1:Qdof1))
                enddo

                ! product (Dwi, Dphi0)
                do l=1, dofP
                   do n=1,nbDim
                      rhsSff(l, k) = rhsSff(l, k)  &
                           + weight * dot_product(Dwi(1:Qdof1, n), DphiW(l, n, 1:Qdof1) )
                   enddo

                enddo

             enddo

             deallocate(Fx, xi, phi, phiW, wi, weights)
             deallocate(Dphi, DphiW, Dwi, Dphi0)

          endif   ! if (weight > 0.) then

       endif  !(i1 > 0)

    enddo ! j=1, N_elem

!!!111 continue

    ! the H1-norm (comment for the L^2-norm)
    L2weight = 1.0
    if(R_type < 0) then
       Mass(1:dofP, 1:dofP) = Mass(1:dofP, 1:dofP) +  L2weight * Stiff(1:dofP, 1:dofP)
       rhs(1:dofP, 1:ndimL) = rhs(1:dofP, 1:ndimL) +  L2weight * rhsSff(1:dofP, 1:ndimL)
    endif

    call SolveLocalMatrixProblem(dofP, Mass(1:dofP, 1:dofP), ndimL, rhs(1:dofP, 1:ndimL) )
    Dw(1:ndimL, 0, 0, 1:dofP) = transpose( rhs(1:dofP, 1:ndimL))

    ! evaluation of the coefficients of the reconstruction
    !call MblockInverse(dofP, Mass(1:dofP, 1:dofP) )

    !do k=1, ndimL
    !   ! zero order derivative = a reconstruction of  function w
    !   Dw(k, 0, 0, 1:dofP)  = matmul(Mass(1:dofP, 1:dofP), rhs(1:dofP, k) )
    !
    !   !if(elem%i == 1199) write(*,'(a6,2i5, 30es12.4)') 'DE45^',elem%i, dofP,  Dw(k, 0, 0, 1:dofP)
    !
    !enddo


    deallocate( Mass, Stiff, rhs, rhsSff )

  end subroutine LeastSquareInterpolationWeighted



  !> higher order reconstruction keeping the derivatives up to degree elem%deg on elem
  !> keeping derivatives on neighbouring elements in least square
  !> USES elem%wS and elem%wSD computed in Compute_means_values_DW !!!
  subroutine LeastSquareRecovery_Dw_mean_values(elem, ndimL, MG, Qnum, degP, dofP, Dw, R_type  )
    class(element), intent(inout), target :: elem
    class(element), pointer :: elem1
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    logical, intent(in) :: MG    ! .false. => use elem%dof, .true. => use elem%MGdof,
    integer, intent(in) :: Qnum  ! the prescribed quadrature
    integer, intent(in) :: degP, dofP  ! the desired degree of the reconstruction
    real, dimension(1:ndimL, 0:degP, 0:degP, 1:dofP), intent(inout) :: Dw
    integer, intent(in) :: R_type  ! = -1 H1-norm, 0 - L2-norm
    type(volume_rule), pointer :: V_rule
    real, dimension(:, :) , allocatable :: wR, wRi, ww
    real, dimension(:,:), allocatable :: psi    ! canonical basis functions
    real, dimension(:, :), pointer :: phi       ! DG basis functions
    real, dimension(:,:), allocatable :: Mass ! Local mass matrix of order deg+1!!
    real, dimension(:,:,:), allocatable :: DpsiK ! canonical basis functions in barycentre of neighbours
    real, dimension(:,:), allocatable :: A, M  ! matrices for the least squares
    real, dimension(:,:), allocatable :: b, x  !vectores for the least squares

    real, dimension(:,:), allocatable :: Fx, FxK

    integer :: dof, Qdof, i, j, ie, l, ideg, iphi, iphiD, k, nx, ny, max_deg
    real :: rphiD

    dof = elem%dof

    ! quadrature rule
    V_rule => state%space%V_rule(Qnum)
    Qdof = V_rule%Qdof
    phi => V_rule%phi(1:dofP, 1:Qdof)

    allocate(Fx(1:Qdof, 1:2) )
    call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )

    ! canonical basis functions in integ nodes
    allocate( psi(1:dofP, 0:Qdof) )

    call Eval_Canonical_Functions(elem, degP, dofP, V_rule, psi(1:dofP, 0:Qdof) )

     ! derivatives of the canonical basis functions in barycentres of neighbours
    allocate(DpsiK(1:dofP, 1:dofP,  1:elem%flen) )  ! 1st= function idx, 2nd = derivative
    DpsiK = 0.

    ! relative coordinates of the barycentre of the neighbour
    allocate(FxK(1:elem%flen, 1:nbDim) )

    nx = 0           ! maximal number of constrains
    max_deg = 0      ! maximal degree of neighbours

    do j = 1, elem%flen
       ie = elem%face(neigh, j)
       if(ie > 0) then ! the neighbour
          elem1 => grid%elem(ie)

          nx = nx + elem1%dof
          max_deg = max(max_deg, elem1%deg)

          FxK(j, 1:2) = (elem1%xc(1:2) - elem%xc(1:2)) / elem%diam
          ! Fx - elem%xc
          call Canonical_Functions(degP, dofP, 1, FxK(j, 1:2), DpsiK(1:dofP, 1, j) )

          !write(33, *) FxK(j, 1:2)*elem%diam +  elem%xc(1:2), DpsiK(1:dofP, 1, j)

          DpsiK(2:dofP,1, j) = DpsiK(2:dofP, 1, j) - psi(2:dofP, 0) !normalization to the zero-mean value

          !write(34, *) FxK(j, 1:2)*elem%diam +  elem%xc(1:2), DpsiK(1:dofP, 1, j)

          ! derivatives of the basis functions in the neighbouring barycentres
          do ideg = 0, degP   ! order of the derivative
             do i = ideg, 0, -1  ! d^ideg / d x^i / d y^{ideg-i}

                k = Deriv_idx(i, ideg - i)

                do iphi = 1, dofP  ! index of the basis function
                   ! computing of the derivatives in ineg nodes
                   call Set_canonical_basis_index_deriv(ideg, i, iphi, iphiD, rphiD)

                   ! the correct normalization since iphi and iphiD has different degress of elem%diam
                   rphiD = rphiD  / elem%diam**ideg

                   if(iphiD > 0) then
                      DpsiK( iphi, k,  j) =  ( DpsiK(iphiD, 1, j) + psi(iphiD, 0) ) * rphiD
                      !      ^^^^  --------------  ^^^^^  OK       ^^^^^^ the mean-value-free factor
                   endif
                enddo
             enddo
          enddo

          !do k=1,dofP
          !   write(*,'(a8, 3i5, 30es12.4)') 'ders:',j, ie, k, DpsiK(k, 1:dofP,  j)
          !enddo
          !print*,'---------------------------', elem1%xc, elem%xc

       end if
    end do

    ! setting of the least squares
    ny = dofP - dof  ! number of unknowns
    !nx              ! maximal number of constrains

    if(nx < ny) then
       print*,' nx, ny = ', nx, ny
       stop "insufficient number  of constrains"
    endif

    allocate(A(1:nx, 1:ny), b(1:nx, 1:ndimL) )

    ! number of real constrains = row_index
    nx = 0
    do ideg = 0, max_deg
       do j = 1, elem%flen
          ie = elem%face(neigh, j)
          if(ie > 0) then ! the neighbour
             elem1 => grid%elem(ie)

             if(elem1%deg >= ideg) then
                ! derivatives of degree deg
                do l= ideg * (ideg+1)/2 +1,  ( ideg+1)*(ideg+2)/2
                   ! matrix row
                   nx = nx + 1
                   A(nx, 1 : ny) = DpsiK(dof+1: dofP, l, j)

                   ! rhs row
                   do k=1, ndimL
                      b(nx, k) = elem1%wSD(2, k, l) &
                           - dot_product( elem%wSD(1, k, 1:dof), DpsiK(1: dof, l, j) )
                   enddo
                enddo
             endif
          endif ! (ie > 0) then ! the neighbour
       end do ! j = 1, elem%flen

       !if(nx >= ny ) goto 11 ! we have anough of constrains

    enddo ! deg = 0, max_deg

11  continue
    !do j=1,Qdof
    !   write(33, *) Fx(j, 1:2), psi(1:dofP, j)
    !enddo

    !do j=1, nx
    !   write(*,'(a4,i5,es12.4,a2,40es12.4)') 'Ab:', j, b(j, 1),'|',A(j, :)
    !enddo
    !print*,' nx, ny = ', nx, ny, elem%xc

    ! least square problem
    allocate(M(1:ny, 1:ny), x(1:ny, 1:ndimL) )
    M(1:ny, 1:ny) = matmul(transpose(A(1:nx, 1:ny)), A(1:nx, 1:ny) )
    x(1:ny, 1:ndimL) = matmul(transpose(A(1:nx, 1:ny)), b(1:nx, 1:ndimL) )

    call SolveLocalMatrixProblem(ny, M(1:ny, 1:ny), ndimL, x(1:ny, 1:ndimL) )


    ! recovered functions
    allocate( wR (1:ndimL, 1:dofP) ) ! coeffs of canonical basis
    allocate( wRi(1:ndimL, 1:Qdof) ) ! integ nodes
    allocate( ww (1:ndimL, 1:dofP) ) ! coeffs of DG basis

    ww = 0.
    wR = 0.

    ! functions up to the degree deg are stored
    wR(1:ndimL, 1:dof) = elem%wSD(1, 1:ndimL, 1:dof)

    ! adding of the recovery update
    wR(1:ndimL, dof+1: dofP) = transpose( x(1:ny, 1:ndimL) )

    !write(*,'(a8,i5,30es12.4)') ' wR:',dofP,  wR(1, 1:dofP)

    ! recovery in integ nodes
    wRi(1:ndimL, 1:Qdof) = matmul(  wR(1:ndimL, 1:dofP), psi(1:dofP, 1:Qdof) )

    ! recovery in DG basis coefficients
    call Trans_Integ_Nodes_to_Basis(elem, dofP, ndimL, V_rule, &
         wRi(1:ndimL, 1:Qdof), ww(1:ndimL, 1:dofP) )


    !write(*,'(a8,i5,30es12.4)') ' wSD:',dofP,  elem%wSD(1, 1, 1:dof)
    !write(*,'(a8,i5,30es12.4)') ' wR:',dofP,  ww(1, 1:dofP)
    call PlotElemFunction3D(10+state%space%adapt%adapt_level, elem,  dofP, ww(1, 1:dofP)  )


    !do i=1,Qdof
    !   write(22, *) Fx(i, 1:2),  wRi(1:ndimL, i)
    !enddo

    deallocate (Fx)
    deallocate(FxK)
    deallocate(DpsiK)

    !stop ":uehde3i"
  end subroutine LeastSquareRecovery_Dw_mean_values


  !> gradient recovery procedure by the least square fitting
  !> USES elem%w
  subroutine LeastSquareInterpolationNodes(elem, RDwh  )
    class(element), intent(inout) :: elem
    real, dimension(1:ndim,  1:elem%dof, 1:2 ), intent(inout) :: RDwh
    class(element), pointer :: elem1
    type(volume_rule), pointer :: V_rule
    integer, dimension(:), pointer :: idx_V_rule
    real, dimension(:,:), allocatable :: A, b, x
    real, dimension(:,:), allocatable :: phi ! the test functions
    real, dimension(:,:,:), allocatable :: Dphi ! Der of test functions
    real, dimension(:,:), allocatable :: xi ! barycentric coordinates
    real, dimension(:,:), allocatable :: Fx ! real physical coordinates

    integer :: deg, dof, Qdof,  Qnum1, Qdof1, dof1, N_nodes, N_elem
    integer :: i, i1, j, k, l, l1, l2, n, n1, ip
    real :: weight

    deg = elem%deg
    dof = elem%dof

    allocate(idx_V_rule(0:elem%flen) )

    ! setting of quadratrure nodes for reconstruction
    !idx_V_rule(0) = elem%deg +2
    !idx_V_rule(0) = elem%Qnum
    idx_V_rule(0) = state%space%Qdeg(elem%deg -1, 1)
    do j=1,elem%flen
       if(elem%face(neigh, j) > 0) then

          elem1 => grid%elem(elem%face(neigh, j))
          !idx_V_rule(j) = elem1%deg +2
          !idx_V_rule(j) = elem1%Qnum
          idx_V_rule(j) = state%space%Qdeg(elem1%deg -1, 1)
       endif
    enddo

    !Qdof = state%space%V_rule(Qnum)%Qdof
    !Qdof = elem%Qdof

    !Qdof = state%space%V_rule( state%space%Qdeg(elem%deg -1, 1) )%Qdof  !!!! elem%Qnum-1)%Qdof
    Qdof = state%space%V_rule(idx_V_rule(0) )%Qdof  !!!! elem%Qnum-1)%Qdof

    N_nodes = Qdof

    do j=1,elem%flen
       if(elem%face(neigh, j) > 0) then
          elem1 => grid%elem(elem%face(neigh, j))

          N_nodes = N_nodes + state%space%V_rule( idx_V_rule(j) )%Qdof
       endif
    enddo


    if(N_nodes < dof) then
       print*,'Under-determined system in ama-L2interpol.f90'
       print*,'#DEDE#@',elem%i, Qdof, N_nodes
       stop
    endif

    allocate(A(1:N_nodes, 1: dof), b(1:N_nodes, 1: 2*ndim) )
    A(:,:) = 0.

    ! evaluation of matrix coefficients, test functions in the nodes for recontruction

    ! the central element
    ip = 0
    V_rule => state%space%V_rule( idx_V_rule(0) )
    Qdof = V_rule%Qdof


    allocate( Fx(1:Qdof, 1:2) )
    call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )

    allocate( Dphi(1:dof, 1:nbDim, 1:Qdof) )  ! phi nad Dphi are different test functions!
    call Eval_Dphi_plus(elem, V_rule, dof,  Dphi(1:dof, 1:nbDim, 1:Qdof) )

    do l = 1, Qdof
       ip = ip + 1
       A(ip, 1:dof) = V_rule%phi(1:dof, l)

       do k=1,ndim
          b(ip, k )       = dot_product(elem%w(0, (k-1)*dof + 1: k*dof), Dphi(1:dof, 1, l ) )
          b(ip, ndim + k) = dot_product(elem%w(0, (k-1)*dof + 1: k*dof), Dphi(1:dof, 2, l ) )
       enddo

       !write(20,*) V_rule%lambda(l, 1:2), Fx(l, 1:2), A(ip, 1:dof)
       !write(30,*) V_rule%lambda(l, 1:2), Fx(l, 1:2), b(ip, : )
    enddo

    deallocate(Dphi, Fx)

    ! neighbouring elements
    do j=1, elem%flen
        if(elem%face(neigh, j) > 0) then
          elem1 => grid%elem(elem%face(neigh, j))

          V_rule => state%space%V_rule( idx_V_rule(j) )
          Qdof1 = V_rule%Qdof

          ! Fx integ nodes on elem1 - real physical cordinates
          ! xi barycentric coordinates of Fx (integ nodes on elem1) with respect elem
          allocate(Fx(1:Qdof1, 1:2), xi(1:Qdof1, 1:2) )
          call ComputeF(elem1, Qdof1, V_rule%lambda(1:Qdof1, 1:2), Fx(1:Qdof1, 1:2) )
          call BarycCoord(elem, Qdof1, Fx(1:Qdof1, 1:2),  xi(1:Qdof1, 1:2) )


          allocate(phi(1:dof, 1:Qdof1) ) ! value of the test function in integ nodes

          call Eval_phi_Qnode(elem, dof, Qdof1, xi(1:Qdof1, 1:nbDim), phi(1:dof, 1:Qdof1) )


          ! for RHS
          dof1 = elem1%dof
          allocate( Dphi(1:dof1, 1:nbDim, 1:Qdof1) )  ! phi nad Dphi are different test functions!
          call Eval_Dphi_plus(elem1, V_rule, dof1,  Dphi(1:dof1, 1:nbDim, 1:Qdof1) )

          do l = 1, Qdof1
             ip = ip + 1
             A(ip, 1:dof) = phi(1:dof, l)

             do k=1,ndim
                b(ip, k )       = dot_product(elem1%w(0, (k-1)*dof1 + 1: k*dof1), Dphi(1:dof1, 1, l ))
                b(ip, ndim + k) = dot_product(elem1%w(0, (k-1)*dof1 + 1: k*dof1), Dphi(1:dof1, 2, l ))
             enddo

             !write(20,*) xi(l, 1:2), Fx(l, 1:2), A(ip, 1:dof)
             !write(30,*) xi(l, 1:2), Fx(l, 1:2), b(ip, :)
          enddo


          deallocate (xi, Fx, phi, Dphi)

       endif
    enddo


    ! ! least squares - NEW variant
    allocate( x(1:dof, 2*ndim) )
    call SolveMatrixProblemLeastSquares(N_nodes, dof, A(1:N_nodes, 1:dof), 2*ndim,  &
         b(1:N_nodes, 1: 2*ndim), x(1:dof, 1:2*ndim) )

    ! setting  of the output
    do k=1,ndim
       RDwh(k, 1:dof, 1) = x(1:dof, k)
       RDwh(k, 1:dof, 2) = x(1:dof, k + ndim)
    enddo

    !print*
    !do l=1,dof
    !   write(*,'(a8,i5, 30es12.4)') 'sol:', l, x(l,:)
    !enddo

    ! graphical checking
    ! allocate( Fx(1:Qdof, 1:2) )
    ! call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )

    ! do l = 1, Qdof
    !    do k=1,ndim
    !       write(40,*) V_rule%lambda(l, 1:2), Fx(l, 1:2),  &
    !            dot_product ( x(1:dof, k ),      V_rule%phi(1:dof,  l ) ) ,  &
    !            dot_product ( x(1:dof,ndim+ k ), V_rule%phi(1:dof,  l ) )
    !    enddo

    ! enddo

    ! deallocate(Fx)

    deallocate(idx_V_rule )
    deallocate(A, b,x)

  end subroutine LeastSquareInterpolationNodes


  !> Polynomial Preserving least squares interpolation
  !> USES elem%w
  subroutine LeastSquareInterpolationPP(elem, dofP, RDwh  )
    class(element), intent(inout) :: elem
    integer, intent(in) :: dofP
    real, dimension(1:ndim,  1:dofP, 0:2 ), intent(inout) :: RDwh
    class(element), pointer :: elem1
    type(volume_rule), pointer :: V_rule
    integer, dimension(:), pointer :: idx_V_rule
    real, dimension(:,:), allocatable :: A, b, x
    real, dimension(:,:), allocatable :: phiP ! the test functions
    real, dimension(:,:,:), allocatable :: Dphi, DphiP ! Der of test functions
    real, dimension(:,:), allocatable :: xi ! barycentric coordinates
    real, dimension(:,:), allocatable :: Fx ! real physical coordinates
    real, dimension(:), allocatable :: weights
    real, dimension(:,:), allocatable :: wi
    logical :: iprint = .false.

    integer :: deg, dof,  dofM, Qdof,  Qnum1, Qdof1, dof1, N_eq
    integer :: i, i1, j, k, l,  ip

    !iprint = .true.

    deg = elem%deg
    dof = elem%dof

    allocate(idx_V_rule(0:elem%flen) )

    ! setting of quadratrure nodes for reconstruction
    idx_V_rule(0) = state%space%Qdeg(elem%deg + 1, 1)
    do j=1,elem%flen
       if(elem%face(neigh, j) > 0) then

          elem1 => grid%elem(elem%face(neigh, j))
          idx_V_rule(j) = state%space%Qdeg(elem1%deg + 1, 1)
       endif
    enddo

    !if(dofP /= DOFtriang( deg+1 ) ) then
    !   print*,'## mishmatch DoF of the polynomial reconstruction of degree elem%deg + 1 !!!'
    !endif

    ! number of equations
    N_eq = dof  !  element itself

    ! its neighbours
    do j=1,elem%flen
       if(elem%face(neigh, j) > 0) then
          elem1 => grid%elem(elem%face(neigh, j))

          N_eq = N_eq + elem1%dof
       endif
    enddo


    if(N_eq < dofP) then
       print*,'Under-determined system in ama-L2interpol.f90'
       print*,'#DED23E#@',elem%i, dofP, N_eq
       stop
    endif

    allocate(A(1:N_eq, 1: dofP), b(1:N_eq, 1: 3*ndim) )
    A(:,:) = 0.
    b(:,:) = 0.

    ! evaluation of matrix coefficients

    ! the central element
    ip = 0
    V_rule => state%space%V_rule( idx_V_rule(0) )
    Qdof = V_rule%Qdof


    !if(iprint) then
    !   print*
    !   write(*,'(a8,i5, 300es12.4)') 'orig:', elem%i, elem%w(0, 1:dof )
    !endif

    allocate( Fx(1:Qdof, 1:2) )
    call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )

    allocate(weights(1:Qdof), wi(0:2, 1:Qdof) )
    call Eval_V_Weights_plus(elem, V_rule, weights(1:Qdof))

    allocate( Dphi(1:dofP, 1:nbDim, 1:Qdof) )
    call Eval_Dphi_plus(elem, V_rule, dofP, Dphi(1:dofP, 1:nbDim, 1:Qdof) )

    do l = 1, dof
       ip = ip + 1

       A(ip, 1:dofP) = matmul(V_rule%phi(1:dofP, 1:Qdof), weights(1:Qdof) *  V_rule%phi(l, 1:Qdof) )

       do k=1,ndim
          wi(0, 1:Qdof)  = matmul(elem%w(0, (k-1)*dof + 1: k*dof), V_rule%phi(1:dof, 1:Qdof ) )
          wi(1, 1:Qdof)  = matmul(elem%w(0, (k-1)*dof + 1: k*dof), Dphi(1:dof, 1, 1:Qdof ) )
          wi(2, 1:Qdof)  = matmul(elem%w(0, (k-1)*dof + 1: k*dof), Dphi(1:dof, 2, 1:Qdof ) )

          b(ip, k          ) = dot_product( weights(1:Qdof) *  V_rule%phi(l, 1:Qdof), wi(0,1:Qdof) )
          b(ip, k + ndim   ) = dot_product( weights(1:Qdof) *  V_rule%phi(l, 1:Qdof), wi(1,1:Qdof) )
          b(ip, k + 2*ndim ) = dot_product( weights(1:Qdof) *  V_rule%phi(l, 1:Qdof), wi(2,1:Qdof) )
       enddo
    enddo

    if(iprint) then
       do l=1,Qdof
          write(20, *) Fx(l, :), wi(0:2,l)
       enddo
    endif


    deallocate(weights, wi, Dphi )
    deallocate(Fx)

    ! neighbouring elements
    do j=1,elem%flen
       if(elem%face(neigh, j) > 0) then

          elem1 => grid%elem(elem%face(neigh, j))

          dof1 = elem1%dof
          dofM = max(dof1, dofP)

          V_rule => state%space%V_rule( idx_V_rule(j) )
          Qdof1 = V_rule%Qdof


          allocate(weights(1:Qdof1), wi(0:2, 1:Qdof1) )
          call Eval_V_Weights_plus(elem1, V_rule, weights(1:Qdof1))

          allocate(Fx(1:Qdof1, 1:2), xi(1:Qdof1, 1:2) )
          call ComputeF(elem1, Qdof1, V_rule%lambda(1:Qdof1, 1:2), Fx(1:Qdof1, 1:2) )
          call BarycCoord(elem, Qdof1, Fx(1:Qdof1, 1:2),  xi(1:Qdof1, 1:2) )


          allocate( Dphi(1:dof1, 1:nbDim, 1:Qdof1) )
          call Eval_Dphi_plus(elem1, V_rule, dof1, Dphi(1:dof1, 1:nbDim, 1:Qdof1) )

          !print*,'AAAAAAA  17.3', elem%i, dofM, Qdof1


          !if(elem%i >= 21 .and. elem%i <= 22) print*,'AAed6fAA  17.*5', dofM, Qdof1, nbDim, elem%i

          allocate(phiP(1:dofM, 1:Qdof1) ) ! value of the test function in integ nodes

          !if(elem%i >= 21 .and. elem%i <= 22) print*,'Aed3eAAA  17.35', dofM, Qdof1, nbDim

          !print*,'AAAAAAA  17.x', elem%i, dofM, nbDim, Qdof1

          allocate(DphiP(1:dofM, 1:nbDim, 1:Qdof1) ) ! value of Deriv test functions

          !print*,'AAAAAAA  17.4', elem%i, dofP, j, elem%face(neigh, j)
          call Eval_phi_Qnode(elem, dofM, Qdof1, xi(1:Qdof1, 1:nbDim), &
               phiP(1:dofM, 1:Qdof1), DphiP(1:dofM, 1:nbDim, 1:Qdof1) )

          !print*,'AAAAAAA  17.5', elem%i, dofP, j, elem%face(neigh, j)

          do l = 1, dof1
             ip = ip + 1

             A(ip, 1:dofP) = matmul(phiP(1:dofP, 1:Qdof1), weights(1:Qdof1) *  phiP(l, 1:Qdof1) )

             do k=1,ndim
                wi(0, 1:Qdof1)  = matmul(elem1%w(0, (k-1)*dof1 + 1: k*dof1), V_rule%phi(1:dof1, 1:Qdof1))
                wi(1, 1:Qdof1)  = matmul(elem1%w(0, (k-1)*dof1 + 1: k*dof1), Dphi(1:dof1, 1, 1:Qdof1))
                wi(2, 1:Qdof1)  = matmul(elem1%w(0, (k-1)*dof1 + 1: k*dof1), Dphi(1:dof1, 2, 1:Qdof1))

                b(ip, k        )  = dot_product( weights(1:Qdof1) *  phiP(l, 1:Qdof1), wi(0,1:Qdof1))
                b(ip, k + ndim )  = dot_product( weights(1:Qdof1) *  phiP(l, 1:Qdof1), wi(1,1:Qdof1))
                b(ip, k + 2*ndim) = dot_product( weights(1:Qdof1) *  phiP(l, 1:Qdof1), wi(2,1:Qdof1))
             enddo
          enddo

          if(iprint) then
             do l=1,Qdof1
                write(30+j, *) Fx(l, :), wi(0:2,l)
             enddo
          endif


          !print*,'......F'
          deallocate(weights)
          deallocate(wi)
          deallocate(phiP)
          deallocate(DphiP)
          deallocate(Dphi)
          deallocate(xi)
          deallocate(Fx )

          !print*,'......G'

       endif
    enddo


    !do i=1,N_eq
    !   write(*,'(a5,i5,3es12.4, a2, 30es12.4)') 'b, A:',i, b(i, 1:3), '|',A(i, 1:)
    !enddo

    ! ! least squares - NEW variant
    allocate( x(1:dofP, 3*ndim) )
    call SolveMatrixProblemLeastSquares(N_eq, dofP, A(1:N_eq, 1:dofP), 3*ndim,  &
         b(1:N_eq, 1: 3*ndim), x(1:dofP, 1:3*ndim) )

    ! setting  of the output
    do k=1,ndim
       RDwh(k, 1:dofP, 0) = x(1:dofP, k)
       RDwh(k, 1:dofP, 1) = x(1:dofP, k + ndim)
       RDwh(k, 1:dofP, 2) = x(1:dofP, k + 2*ndim)
    enddo

    !if(iprint) then
    !   !print*
    !   write(*,'(a8,i5, 300es12.4)') 'sol:', elem%i, x(1:dofP, 1)
    !   !write(*,'(a8,i5, 300es12.4)') 'sol:', 1, x(1:dofP, 2)
    !   !write(*,'(a8,i5, 300es12.4)') 'sol:', 2, x(1:dofP, 3)
    !endif

    ! graphical checking
    if(iprint) then
       V_rule => state%space%V_rule( idx_V_rule(0) )
       allocate( Fx(1:Qdof, 1:2) )
       call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )

       do l = 1, Qdof
          do k=1,ndim
             write(40,*) Fx(l, 1:2), &
                  dot_product (RDwh(k, 1:dofP, 0 ),  V_rule%phi(1:dofP,  l ) ), &
                  dot_product (RDwh(k, 1:dofP, 1 ),  V_rule%phi(1:dofP,  l ) ), &
                  dot_product (RDwh(k, 1:dofP, 2 ),  V_rule%phi(1:dofP,  l ) )
          enddo

       enddo

       deallocate(Fx)
    endif


     !stop 'EDT#%WT'

     deallocate(idx_V_rule )
     deallocate(A, b,x)

  end subroutine LeastSquareInterpolationPP



  !> computation of the secondorder derivatives for P_0 approximation
  !> higher order reconstruction via least square method in the \f$ H^1\f$-norm,
  !> \f$ P^{degP} \f$ approximation
  !> over elem and its neighbours given
  !> by \f$ \int_{D_e} \sum_j w_j \varphi_j \varphi_i dx =  \int_{D_e} w_h \varphi_i dx\f$,
  !> \f$ i=1,\dots dofP \f$
  !> USES elem%wS !!!
  subroutine P_0_to_P2_Interpolation(elem, ndimL, MG, Qnum, degP, dofP, Dw, R_type  )
    class(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    logical, intent(in) :: MG    ! .false. => use elem%dof, .true. => use elem%MGdof,
    integer, intent(in) :: Qnum  ! the prescribed quadrature
    integer, intent(in) :: degP, dofP  ! the desired degree of the reconstruction
    real, dimension(1:ndimL, 0:degP, 0:degP, 1:dofP), intent(inout) :: Dw
    integer, intent(in) :: R_type  ! = -1 H1-norm, 0 - L2-norm
    class(element), pointer :: elem1
    type(volume_rule), pointer :: V_rule

  end subroutine P_0_to_P2_Interpolation


end module ama_L2interpol

