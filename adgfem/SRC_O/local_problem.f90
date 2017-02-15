!> solution of the local problem for one element used for the error estimates

module local_problem
  use main_data
  use problem_oper
  use lin_solvers
  use inviscid_fluxes
  use time_sets
  use matrix_oper_int
  use euler_problem
  use ama_L2interpol
  implicit none
  
  public:: SolveLocalProblem
  public:: SolveLocalNewton
  public:: ComputeLocalBlocks
  public:: PrepareLocalElement
  public:: DeleteLocalElement
  public:: ComputeLocalMatrixBlock
  public:: SolvePPlocalProblem

contains

  !> solve the local problem on one element
  subroutine SolveLocalProblem(elem, degP, dofP, Rhw) 
    class(element), intent(in) :: elem
    integer, intent(in) :: degP, dofP
    real, dimension(-2:1, 1:ndim, 1:dofP), intent(inout) :: Rhw
    class(element), pointer :: elemP, elem1, elemP1
    integer :: degP1
    logical :: loc_implicitly
    integer :: i,j,k, l, l1, imp, ideg, flenP, flen

    Rhw(:,:,:) = 0.

    ! compute only diag blocks for local problems for EE
    !state%only_diag = .true.
    state%local_problem = .true.

    !print*,'SolveLocalProblem started'
    loc_implicitly = state%nlSolver%implicitly

    ! allocation of new elements for the reconstruction
    flenP = elem%flen
    grid%flenP = flenP

    allocate(grid%elemP(0: flenP) ) 
    grid%elemP(0: flenP)%deg = -1

    allocate(grid%elemP_idx(0:flenP, 1:flenP) )
    grid%elemP_idx(:,:) = -1000

    !allocate(elemP)
    elemP => grid%elemP(0) 
    call PrepareLocalElement(elemP, degP, elem)
    if( elemP%dof /= dofP) stop'Inconsistency (1) in SolveLocalProblem in local_problem.f90'
    elemP%i = -elem%i

    do j = 1, flenP
       k = elemP%face(neigh,j)
       if(k > 0) then
          grid%elemP_idx(0,j) = -j       !pointers to neighbouring elements


          elem1  => grid%elem(k)
          elemP1 => grid%elemP(j)
          !degP1 = elem1%deg + (degP - elem%deg)
          degP1 = degP

          call PrepareLocalElement(elemP1, degP1, elem1)
          elemP1%i = -elem1%i
       endif
    enddo

    !pointers to neighbouring elements of the neighbouring elements
    do j = 1, flenP
       k = elemP%face(neigh,j)
       if(k > 0) then
          elem1  => grid%elem(k)

          flen = elem1%flen
          if(flen /= flenP) stop 'TRouble 35d4276a at local_problem.f90'
          do l=1,flen
             grid%elemP_idx(j,l) = elem1%face(neigh, l)
             if(elem1%face(neigh, l) <= 0) grid%elemP_idx(j,l) = -1000
          enddo
          
          do l=1, flen
             if(elem1%face(neigh, l) > 0 ) then
                do l1 = 0, flen

                   if( elem1%face(neigh,l) ==  -grid%elemP(l1)%i ) then
                      grid%elemP_idx(j, l) = -l1
                   endif
                enddo
             endif
          enddo

       endif
    enddo


    ! do i=0, flenP
    !    do j=1, flenP
    !       write(*,'(a9, 20i5)') 'elemP:::',i,j,grid%elemP_idx(i,j), grid%elemP(i)%i
    !    enddo
    !    print*
    ! enddo


    ! solution of the local problems
    !do ideg = 1, 1
    ideg = 1

    elemP => grid%elemP(0) 
    !elemP%deg = elem%deg + ideg
    !elemP%dof = DOFtriang( elemP%deg)
    
    !write(*,'(a6,2i5,200es12.4)') 'orig',ideg, elem%dof, elem%w(0,:)
    !write(*,'(a6,2i5,200es12.4)') 'WEDS',ideg, elemP%dof, Rhw(ideg, 1, :)
    !print*
    ! computing of the auxiliarly local problem via Newton method
    !call SolveLocalNewton( )   ! NOT VERIFIED
    call SolvePPlocalProblem( )
    
    ! setting of the reconstruction
    do k=1,ndim
       j = (k-1)*elemP%dof + 1
       i = k*elemP%dof
       Rhw(ideg, k,1: elemP%dof) = elemP%w(0,j:i)
    enddo
    
    !write(*,'(a6,2i5,200es12.4)') 'orig',ideg, elem%dof, elem%w(0,:)
    !write(*,'(a6,2i5,200es12.4)') 'WEDS',ideg, elemP%dof, Rhw(ideg, 1, :)
    !print*
    !print*,'------------------------------------------'
    !!enddo


    ! deallocation of the auxiliarlu element
    call DeleteLocalElement(grid%elemP(0)) 
    !deallocate(elemP)

    do j = 1, flenP
       if(elem%face(neigh,j) > 0) call DeleteLocalElement(grid%elemP(j)) 
    enddo


    ! back the original values
    !state%only_diag = .false.
    state%local_problem = .false.
    state%nlSolver%implicitly = loc_implicitly

    !stop 'SolveLocalProblem, elem%i wsuy563ed'

    deallocate(grid%elemP)
    deallocate(grid%elemP_idx) 

  end subroutine SolveLocalProblem


  subroutine SolveLocalNewton( )
    integer :: flenP
    class(element), pointer :: elemP
    real, dimension(:), allocatable :: x, b, b1  ! vectors for newton
    real, dimension(:,:), allocatable :: mA
    !logical :: vector_update
    real :: res0, res1
    real :: Newton_tol, Newton_max_iter, theta, lambda, lambda_old
    integer :: iter, ndof, l
    integer :: i,j,k, ie

    flenP = grid%flenP

    Newton_tol = 1E-8
    Newton_max_iter = 3

    elemP => grid%elemP(0)

    ! arrays for Newton
    ndof = elemP%dof * ndim
    allocate(x(1:ndof), b(1:ndof), b1(1:ndof), mA(1:ndof, 1:ndof)  )

    state%nlSolver%implicitly = .true.
    ! Newton iterations
    do iter = 1, Newton_max_iter
       !print*,'###############################    Newton_iter', iter,'       elem=',elemP%i

       if( state%nlSolver%implicitly ) then
          call ComputeLocalBlocks( )  ! matrix block
       endif

       state%nlSolver%implicitly = .false.
       !print*,'------------------------**************-----------------------------------------------'

       if(iter == 1) then
          call ComputeLocalBlocks( )  ! vector

          b(1:ndof) = elemP%vec(rhs, 1:ndof)
          res0 = VectorNorm(b(1:ndof))

       else
          b(1:ndof) = b1(1:ndof)
          res0 = res1
       endif

       if(elemP%i == -1) then
          do k=1,ndim
             do i=1,elemP%dof
                j = (k-1)*elemP%dof + i
                write(*,'(a3,3i5,300es12.4)') 'As',k,i,j, &
                     elemP%block(0)%Mb(j, :), &
                     elemP%vec(rhs, j)
             enddo
          enddo
          print*
          do ie=1,elemP%flen
             print*,'----------',ie, elemP%face(neigh, ie)
             if(elemP%face(neigh, ie) > 0) then
                do k=1,ndim
                   do i=1,elemP%dof
                      j = (k-1)*elemP%dof + i
                      write(*,'(a3,3i5,300es12.4)') 'As',k,i,j, &
                           elemP%block(ie)%Mb(j, :), &
                           elemP%vec(rhs, j)
                   enddo
                enddo
             endif
          enddo
       endif

       
       if(iter == 1 ) &
            write(*,'(a6,i5,3es12.4,a3,300es12.4)') 'iters:',0, res0,-2., -2.,'|',elemP%w(0,1:6)



       !if(iter == 1) write(*,'(a6,i5,3es12.4,a3,300es12.4)') 'iters:',0, res0,0.,0.,'|',elemP%w(0,:)


       x(1:ndof) = b(1:ndof)

       mA(1:ndof, 1:ndof) = elemP%block(0)%Mb(1:ndof, 1:ndof)

       ! print*,'--------------------------------------------------'
       ! do j = 1,ndof
       !    write(*,'(a3,i5,300es12.4)') 'A',j, mA(j,:), x(j)
       ! enddo
       !print*,'--------------------------------------------------',iter, res0, res1

       call SolveLocalMatrixProblem(ndof, mA(1:ndof, 1:ndof), 1, x(1:ndof) )

       lambda_old = 0.
       lambda = 1.0   ! initialization of lambda (= damping factor

       do l=1,3    ! iterations, seeking the optimal damping factor

          elemP%w(0,1:ndof) = elemP%w(0,1:ndof) + (lambda - lambda_old) * x(1:ndof)
          lambda_old = lambda

          state%nlSolver%implicitly = .false.
          call ComputeLocalBlocks( )  ! vector

          b1(1:ndof) = elemP%vec(rhs, 1:ndof)
          res1 = VectorNorm(b1(1:ndof))

          theta = res1 / res0
          write(*,'(a6,i5,3es12.4,a3,300es12.4)') 'iters:',iter, res1,theta, lambda,'|',elemP%w(0,1:6)

          if(res1 < Newton_tol) goto 20

          if(theta > 0.95) then
             lambda = lambda / 2
          else
             goto 10
          endif
       enddo
10     continue
       if(lambda < 1.) state%nlSolver%implicitly = .true.
       
    enddo  ! do iter =1, Newton_max_iter
20  continue

    deallocate( x, b, b1, mA)

  end subroutine SolveLocalNewton


  !> compute either a local matrix block or a block vector
  subroutine ComputeLocalBlocks( )
    integer :: flenP
    class(element), pointer :: elemP, elemP0
    integer :: i, j, k

    flenP = grid%flenP
    elemP0 => grid%elemP(0)

    do i=0, flenP
       grid%elemP_idx_actual = i

       !if(i == 0 .or.  elemP0%face(neigh, max(i,1)) > 0 ) then
       elemP => grid%elemP(i)
       if(elemP%deg >= 0) then ! otherwise no element (we are at the boundary)

          !write(*,'(a6,10i5)') 'elemP:',i,elemP%i, &
          !     size(elemP%block(0)%Mb, 1), &
          !     size(elemP%block(0)%Mb, 1), &
          !     size(elemP%block(3)%Mb, 1), &
          !     size(elemP%block(3)%Mb, 2)

          if(state%modelName == 'scalar') then        ! 2D scalar equation
             call ComputeLocalMatrixBlock(Set_f_s_scalar, Set_A_s_scalar, Set_Ppm_scalar, &
                  Set_R_s_scalar, Set_K_sk_scalar, Set_S_scalar, Set_DS_scalar, elemP)
          else
             stop 'Other variants in SolveLocalProblem not YET implemented in local_problem.f90'
          endif

          do k=0, elemP%flen
             if(k==0 .or. elemP%face(neigh, max(k,1)) > 0) then
                !print*,'#####',elemP%i,k
                !do j=1,elemP%dof
                !   write(*,'(a3,2i5,300es12.4)') 'WDKN',i,j,elemP%vec(RHS,j), elemP%block(k)%Mb(j, :)
                !enddo
                !print*
             endif
          enddo
       end if
       !print*,'### ed49 ##################################', i, flenP
    enddo



  end subroutine ComputeLocalBlocks
 


  !> preparation of the element for the local problem
  subroutine PrepareLocalElement(elemP, degP, elem)
    class(element), intent(inout) :: elemP
    class(element), intent(in) :: elem
    integer, intent(in) :: degP
    integer :: Qnum, dof, flen, FFdof, k, j, nsl, nsl1

    elemP%deg = degP
    elemP%dof = DOFtriang(degP)
    elemP%Tdeg = elem%Tdeg

    elemP%deg_plus = .false.

    ! volume quadratures 
    !elemP%Qdeg = state%space%Qdeg(elemP%deg, 1)
    elemP%Qdeg = elem%Qdeg

    elemP%Qnum = elemP%Qdeg
    elemP%Qdof = state%space%V_rule(elemP%Qnum)%Qdof
    
    ! geoemetry
    flen = elem%flen
    elemP%flen = flen
    elemP%type = elem%type

    allocate(elemP%iBC(1:flen))
    elemP%ibc(1:flen) = elem%ibc(1:flen) 

    elemP%HGnode = .false.

    allocate(elemP%F)
    elemP%F%deg = elem%F%deg
    elemP%F%dof = elem%F%dof
    FFdof = elemP%F%dof

    allocate(elemP%F%F(1:FFdof,1:nbDim))
    elemP%F%F(1:FFdof,1:nbDim) = elem%F%F(1:FFdof,1:nbDim)

    elemP%F%iFlin = elem%F%iFlin
    elemP%ibcur = elem%ibcur

    elemP%rezid = elem%rezid

    ! edges
    allocate(elemP%face(1:max_face, 1:elemP%flen))
    elemP%face(1:max_face, 1:elemP%flen) = elem%face(1:max_face, 1:elemP%flen)

    do j=1,elemP%flen
       if(elemP%face(neigh, j) > 0) then
          elemP%face(fdeg,j) = elemP%face(fdeg,j) + 1
          elemP%face(fdof,j) = DOFtriang( elemP%face(fdeg,j))
       endif
    enddo
          

    ! removed to the call PrepareLocalElement
    !elemP%i = -elem%i

    ! if(  elem%Qnum == elemP%Qnum ) then the following can be replaced by a direct copy from elem
    call ComputeIntegNode(elemP)
    
    !do j=1,flen
    !   write(*,'(a8,20i5)') 'EDE@W',elemP%face(:,j)
    !enddo


    ! normals
    allocate( elemP%n(1:flen,1:nbDim))
    allocate( elemP%dn(1:flen))
    elemP%n(1:flen,1:nbDim) = elem%n(1:flen,1:nbDim)
    elemP%dn(1:flen) = elem%dn(1:flen)

    elemP%area = elem%area

    if(elemP%F%iFlin ) then
       ! element with linear F ==> constant Jacobian
       allocate(elemP%F%D1F0(1:nbDim,1:nbDim) )
       elemP%F%D1F0(1:nbDim, 1:nbDim) = elem%F%D1F0(1:nbDim, 1:nbDim)
       elemP%F%JF0 = elem%F%JF0
    else
       stop'NOT YET implemented ETDSETSWUE in local_problem.f90'
       !allocate(elemP%F%V%D1F(1:Qdof,1:nbDim,1:nbDim) )
       !allocate(elemP%F%V%JF(1:Qdof) )
    endif



    ! arrays
    allocate(elemP%w(0:0, 1:ndim*elemP%dof) ) 
    elemP%w(:,:) = 0.
    do k=1, ndim
       j = (k-1)*elemP%dof
       elemP%w(0, j+1: j+elem%dof) = elem%w(0, (k-1)*elem%dof + 1: k*elem%dof) 
    enddo


    ! local mass matrix
    Qnum = state%space%Qdeg(elem%deg+1,1)
    dof = elemP%dof

    allocate(elemP%Mass)
    call InitMblock(elemP%Mass,dof,dof)
    call IntegrateBlockBBmass(elemP, Qnum, dof, dof, elemP%Mass%Mb(1:dof,1:dof) )

    ! matrix block
    allocate( elemP%block( 0:flen) )
    !!allocate( elemP%block( 0:0) )
    nsl = ndim*dof
    call InitMblock(elemP%block(0), nsl, nsl )
 
    ! off diagonal blocks
    do j=1,elemP%flen
       if(elemP%face(neigh,j) > 0) then
          nsl1 = elemP%face(fdof,j) * ndim
          
          call InitMblock(elemP%block(j), nsl, nsl1)
       endif
       
    enddo
    
    ! vector
    allocate( elemP%vec(1:max_vecs, 1:dof * ndim) ) 
    elemP%vec(:,:) = 0.

  end subroutine PrepareLocalElement

  !> deallocation of the element
  subroutine DeleteLocalElement(elemP) 
    class(element), intent(inout) :: elemP
    integer :: j
    
    deallocate(elemP%xi) 
    
    deallocate( elemP%n)
    deallocate( elemP%dn )

    if(elemP%F%iFlin ) then
       ! element with linear F ==> constant Jacobian
       deallocate(elemP%F%D1F0  )
    else
       stop'NOT YET implemented ETDSETSWUE in local_problem.f90'
       !allocate(elemP%F%V%D1F(1:Qdof,1:nbDim,1:nbDim) )
       !allocate(elemP%F%V%JF(1:Qdof) )
    endif

    deallocate(elemP%F%F)
    deallocate(elemP%F)


    ! arrays
    deallocate(elemP%w, elemP%vec ) 

    deallocate(elemP%Mass%Mb)
    deallocate(elemP%Mass)

    ! matrix block
    deallocate(elemP%block(0)%Mb) 
 
    ! off diagonal blocks, used only for the consistency 
    !do j=1,elemP%flen
    !   if(elemP%face(neigh,j) > 0) deallocate (elemP%block(j)%Mb) 
    !enddo

    deallocate( elemP%block )

    deallocate(elemP%face)

    deallocate(elemP%iBC)

    
    !!print*,'NOT COMPLETED ??'

  end subroutine DeleteLocalElement


  !> compute the matrix block and the RHS for the local problem
  subroutine ComputeLocalMatrixBlock(Set_f_s, Set_A_s, Set_Ppm, Set_R_s, Set_K_sk, Set_S, Set_DS, elemP)
    external :: Set_f_s, Set_A_s, Set_Ppm, Set_R_s, Set_K_sk, Set_S, Set_DS
    class(element), intent(inout) :: elemP
    integer :: i, j

    !print*,'ComputeLocalMatrixBlock started', '  dofP = ', elemP%dof, elemP%Qdof

    if(state%nlSolver%implicitly) then
       elemP%block(0)%Mb(:,:) = 0.
       do j=1,elemP%flen
          if(elemP%face(neigh,j) > 0) then
             elemP%block(j)%Mb(:,:) = 0.
          endif
       enddo
    endif

    elemP%vec(:,:) = 0.
   
    call ComputeOneElementTerms(Set_f_s, Set_A_s, Set_Ppm, Set_R_s, Set_K_sk, Set_S, Set_DS, elemP)


    ! do i=1,elemP%dof
    !    write(*,'(a3,i5,300es12.4)') 'WDKN',i,elemP%vec(RHS,i), elemP%block(0)%Mb(i, :)
    ! enddo
    ! print*,'))))))))))))))))))))))))))))))))))))))))))))))))'

  end subroutine ComputeLocalMatrixBlock

  !> solve the Polynomially Preserving local problem
  subroutine SolvePPlocalProblem( )
    integer :: flenP
    class(element), pointer :: elemP, elemP0
    real, dimension(:,:,:), allocatable :: phi  ! patch basis functions in basis coefficients
    real, dimension(:,:),   allocatable :: Qphi ! patch basis functions in integ nodes
    real, dimension(:,:), allocatable :: Fx, xi
    real, dimension(:), allocatable :: weights
    real, dimension(:,:), pointer :: phi0

    !real, dimension(:), allocatable :: x, b, b1  ! vectors for newton
    real, dimension(:,:), allocatable :: mA, mB, x
    type(volume_rule), pointer :: V_rule
    !logical :: vector_update
    real :: res0, res1
    real :: Newton_tol, Newton_max_iter, theta, lambda, lambda_old
    integer :: iter, ndof, l, N_eq, ndofP
    integer :: i,j,k, ie, ip, dofP, dof, dofM, Qdof

    flenP = grid%flenP
    elemP0 => grid%elemP(0)
    dofP = elemP0%dof

    ! the volume quadrature used for the whole stencil 
    V_rule => state%space%V_rule( state%space%Qdeg(elemP0%deg, 1 ) )
    Qdof = V_rule%Qdof

    ! setting of the patch basis functions in basis coefficients and integ nodes
    allocate(  phi(0:flenP, 1:dofP, 1:dofP) )
    allocate( Qphi( 1:dofP, 1:Qdof) )

    phi(:,:,:) = 0.
    
    do j=1,dofP
       phi(0, j, j) = 1.  ! canonical basis function on the central element
    enddo
    

    ! elements from the stencil
    do i=1, flenP
       elemP => grid%elemP(i)
       if(elemP%deg >=0) then

          ! setting of the patch basis functions in integ nodes of the element from the patch 
          allocate(Fx(1:Qdof, 1:2), xi(1:Qdof, 1:2) )
          call ComputeF(elemP, Qdof, V_rule%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )
          call BarycCoord(elemP0, Qdof, Fx(1:Qdof, 1:2),  xi(1:Qdof, 1:2) )

          call Eval_phi_Qnode(elemP, dofP, Qdof, xi(1:Qdof, 1:nbDim), Qphi(1:dofP, 1:Qdof) )


          ! evaluation of the patch basis functions in basis coeffs on the element from the patch 
          allocate(weights(1:Qdof)  )
          call Eval_V_Weights_plus(elemP, V_rule, weights(1:Qdof))
          phi0 => V_rule%phi(1:dofP, 1:Qdof)
          
          allocate(mA(1:dofP, 1:dofP), mB(1:dofP, 1:dofP) )
          do j=1,dofP
             mA(j, 1:dofP) = matmul( phi0(1:dofP, 1:Qdof), weights(1:Qdof) * phi0(j, 1:Qdof) )
             mB(1:dofP, j) = matmul( phi0(1:dofP, 1:Qdof), weights(1:Qdof) * Qphi(j, 1:Qdof) )
          enddo
          
          call SolveLocalMatrixProblem(dofP, mA(1:dofP, 1:dofP), dofP, mB(1:dofP, 1:dofP) )

          phi(i, 1:dofP, 1:dofP) = transpose( mB(1:dofP, 1:dofP) )  

          deallocate(weights, Fx, xi, mA, mB)

       endif
    enddo  ! = 1, flenP


    !print*,'SolvePP:',dofP, Qdof, elemP0%deg
    !stop


    ! counting the number of equations
    N_eq = 0
    do i=0, flenP  ! we go over all elements from the patch
       elemP => grid%elemP(i)
       N_eq = N_eq + DOFtriang(elemP%deg -1 ) * ndim
    enddo
    if(N_eq < dofP) then
       print*,'small number of equations N_eq = ', n_eq,',     Dof = ', dofP
       stop
    endif


    ndofP = ndim * dofP
    allocate(mA(1:N_eq, 1: ndofP), mb(1:N_eq, 1) ) 
    mA(:,:) = 0.
    mb(:,:) = 0.

    ! evaluation of the forms a_h( psi_i, phi_j)
    state%nlSolver%implicitly = .false.
    

    do ie = 0, dofP ! we go over all patch basis functions (ie = 0 is for the RHS)

       ! setting of w := psi_i
       do i=0, flenP  ! we go over all elements from the patch
          elemP => grid%elemP(i)
          dof = elemP%dof

          if(elemP%deg >= 0) then  ! otherwise no element (we are at the boundary)

             if(ie == 0) then
                elemP%w(0,1:dof) = 0. ! setting of the RHS
             else
                ! setting of test function of the patch
                elemP%w(0,1:dof) = phi(i, ie, 1:dof) 
             endif
             
             !! graphical verification of the basis functions
             !call PlotElemSolution3D(20+ie, elemP)
          endif
       enddo
       
       call ComputeLocalBlocks( )  ! matrix block
       
       ip = 0
       do i=0, flenP
          elemP => grid%elemP(i)
          if(elemP%deg >= 0) then  ! otherwise no element (we are at the boundary)
             ndof = elemP%dof * ndim
             dof = elemP%dof 
             dofM = DOFtriang(elemP%deg -1 )
             
             !!write(*,'(a8,3i5,60es12.4)') 'RHS:', ie, elemP%i, elemP%dof, elemP%vec(rhs,1:ndof)
             
             if(ie == 0) then ! RHS

                do k=0, ndim-1
                   mb(ip+ k*dofM + 1 : ip + (k+1)*dofM, 1) = elemP%vec(rhs, k*dof+1 : k*dof + dofM)
                enddo

             else ! the matric

                do k=0, ndim-1
                   mA(ip+ k*dofM + 1 : ip + (k+1)*dofM, ie) =  &
                      - elemP%vec(rhs, k*dof+1 : k*dof + dofM) + mb(ip+ k*dofM + 1 : ip + (k+1)*dofM, 1)
                enddo
             endif
             ip = ip + dofM * ndim

          endif
          
       enddo
       !print*
    enddo  ! do ie = 0, dofP

    !do i=1,N_eq
    !   write(*,'(a5,i5,1es12.4, a2, 30es12.4)') 'b, A:',i, mb(i, 1), '|',mA(i, 1:)
    !enddo
    
    ! least squares - 
    allocate( x(1:ndofP, 1) )
    call SolveMatrixProblemLeastSquares(N_eq, ndofP, mA(1:N_eq, 1:ndofP), 1,  &
         mb(1:N_eq, 1), x(1:ndofP, 1) )

    ! setting  of the output
    elemP0%w(0,1:ndofP) = x(1:ndofP, 1)
    !call PlotElemSolution3D(50, elemP0)
    

    deallocate(mA, mB, x)
    deallocate(phi, Qphi)
  end subroutine SolvePPlocalProblem



end module local_problem

