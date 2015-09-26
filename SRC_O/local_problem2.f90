!> solution of the local problem for one element used for the error estimates

module local_problem2
  use main_data
  use problem_oper
  use lin_solvers
  use inviscid_fluxes
  use time_sets
  use matrix_oper_int
  use euler_problem
  use ama_L2interpol
  use  mesh_mod
  implicit none

  public:: SolveLocalProblem2
  public:: SolveLocalNewton2
  public:: ComputeLocalBlocks2
  public:: PrepareLocalElement2
  public:: DeleteLocalElement2
  public:: ComputeLocalMatrixBlock2
  public:: SolvePPlocalProblem2

contains

  !> solve the local problem on one element
  subroutine SolveLocalProblem2(elem, degP, dofP, Rhw)
    class(element), intent(in) :: elem
    integer, intent(in) :: degP, dofP
    real, dimension(1:ndim, 1:dofP), intent(inout) :: Rhw
    !class (mesh), allocatable :: gridL
    class(element), pointer :: elemP, elem1, elemP1
    integer :: degP1
    logical :: loc_implicitly
    integer :: i,j,k, l, l1, imp, ideg, flenP, flen


    !print*,'grids:', grid%nelem

    !Rhw(:,:) = 0.

    ! compute only diag blocks for local problems for EE
    !state%only_diag = .true.
    state%local_problem = .true.

    !print*,'SolveLocalProblem started'
    loc_implicitly = state%nlSolver%implicitly

    ! setting of the local patch
    allocate( MeshAMA_t :: gridL )

    !call Set_Local_Patch_interior(elem, degP)
    !call Set_Local_Patch_neighbours(elem)
    call Set_Local_Patch_oneK(elem)

    call Prepare_Local_grid( degP )

    !call Prepare_Local_grid_solution( degP )

    ! works only for one element
    call Prepare_Local_grid_solution_Rhw( degP, dofP, Rhw(1:ndim, 1:dofP) )

    ! soultion of the DG problem on the patch
    call SolveLocalNewton2( ) !gridL )

    !do i=1, 1 !gridL%nelem
    !   elemP => gridL%elem(i)
    !   call PlotElemSolution3D( 60, elemP)
    !enddo

    !stop 'edr34ews34'

    ! setting of the reconstruction
    elemP => gridL%elem(1)

    do k=1,ndim
       j = (k-1)*elemP%dof + 1
       i = k*elemP%dof
       Rhw(k,1: elemP%dof) = elemP%w(0,j:i)
    enddo

    !write(*,'(a6,2i5,200es12.4)') 'orig',ideg, elem%dof, elem%w(0,:)
    !write(*,'(a6,2i5,200es12.4)') 'WEDS',ideg, elemP%dof, Rhw( 1, :)
    !print*
    !print*,'------------------------------------------'
    !!enddo


    ! deallocation of the auxiliarlu element
    do i=1, gridL%nelem
       call DeleteLocalElement2(gridL%elem(i))
    enddo

    deallocate( gridL)


    ! back the original values
    !state%only_diag = .false.
    state%local_problem = .false.
    state%nlSolver%implicitly = loc_implicitly

    !stop 'SolveLocalProblem, elem%i wsuy563ed'

    !deallocate(grid%elemP)
    !deallocate(grid%elemP_idx)

  end subroutine SolveLocalProblem2


  !> setting of the local patch consisting of four subelements
  subroutine Set_Local_Patch_interior(elem, degP) !, gridL)
    class(element), intent(in) :: elem
    !class (mesh),  intent(inout), target :: gridL
    class(element), pointer :: elemP, elem1, elemP1
    integer, dimension (:,:), allocatable:: iloc
    real, dimension (:,:), allocatable:: xi, wi, Qphi,w
    integer, intent(in) :: degP
    integer :: i, j, k, l, l1, flenP, flen
    integer :: dofP, dof, Qdof,Fdof


    dof  = elem%dof
    dofP = DOFtriang(degP)

    flenP = elem%flen
    grid%flenP = flenP
    if(flenP /= 3) stop 'non-implemented cas 4dets5etd76e in local_problem2.f90'

    ! allocation of local mesh for the reconstruction
    call gridL%init(4, 6, 6, 4, 2)

    ! creation of the local grid
    ! spltting of elem onto sub-elements
    gridL%x(1, :) = grid%x(elem%face(idx, 1), : )
    gridL%x(2, :) = grid%x(elem%face(idx, 2), : )
    gridL%x(3, :) = grid%x(elem%face(idx, 3), : )
    gridL%x(4, :) = (gridL%x(1, :) + gridL%x(2, :) ) / 2.
    gridL%x(5, :) = (gridL%x(2, :) + gridL%x(3, :) ) / 2.
    gridL%x(6, :) = (gridL%x(3, :) + gridL%x(1, :) ) / 2.

    allocate(iloc(1:4, 1:3) ) ! number of sub-triangles, 3= triangle
    iloc(1, 1:3) = (/ 4, 5, 6 /)
    iloc(2, 1:3) = (/ 4, 6, 1 /)
    iloc(3, 1:3) = (/ 5, 4, 2 /)
    iloc(4, 1:3) = (/ 6, 5, 3 /)

    do i = 1, gridL%nelem
       elemP => gridL%elem(i)
       call elemP%init(i , 3, 3, iloc(i, 1:3) )

       !elemP%deg_cur = 1  !!!!
       !elemP%ibcur = -1  !!!!
       !elemP%jcur = 0  !!!!
       !elemP%iBC(:) = -1  ! = 0 means Neuman BC, > 1 Dirichlet taken from *.ini, -1 Dirichlet(probably?)
    enddo

    ! boundary segments
    gridL%b_edge(1)%lbn(1:2) = (/ 1, 4 /);    gridL%b_edge(1)%ibc = elem%face(neigh, 1)
    gridL%b_edge(2)%lbn(1:2) = (/ 4, 2 /);    gridL%b_edge(2)%ibc = elem%face(neigh, 1)
    gridL%b_edge(3)%lbn(1:2) = (/ 2, 5 /);    gridL%b_edge(3)%ibc = elem%face(neigh, 2)
    gridL%b_edge(4)%lbn(1:2) = (/ 5, 3 /);    gridL%b_edge(4)%ibc = elem%face(neigh, 2)
    gridL%b_edge(5)%lbn(1:2) = (/ 3, 6 /);    gridL%b_edge(5)%ibc = elem%face(neigh, 3)
    gridL%b_edge(6)%lbn(1:2) = (/ 6, 1 /);    gridL%b_edge(6)%ibc = elem%face(neigh, 3)

    gridL%b_edge(1:gridL%nbelm)%icurv = 0

    call gridL%plot('meshL')

    call gridL%seekNeighbours()

    ! setting of further arrays
    do i = 1, gridL%nelem
       elemP => gridL%elem(i)

       elemP%deg_cur = 1  !!!!
       elemP%ibcur = -1  !!!!
       elemP%jcur = 0  !!!!
       elemP%iBC(:) = -1  ! = 0 means Neuman BC, > 1 Dirichlet taken from *.ini, -1 Dirichlet(probably?)
    enddo

    !write(*,'(a6,30i5)') 'elem:',elem%i, elem%face(idx, :), elem%face(neigh, :)

    ! do i=1,gridL%nelem
    !    write(*,'(a6,30i5)') 'elemL:',i, gridL%elem(i)%face(idx, :), gridL%elem(i)%face(neigh, :)
    ! enddo

    ! do i=1,gridL%nbelm
    !    write(*,'(a6,30i5)') 'B_elL:',i, gridL%b_edge(i)%lbn(1:2), gridL%b_edge(i)%ibc
    ! enddo

    call gridL%computeGeometry()

    do i=1,gridL%nelem
       elemP => gridL%elem(i)
       elemP%deg = degP
       elemP%Tdeg = 0
       elemP%degP = 0

       call elemP%initElementDof()


       allocate(elemP%w(0:1, 1:ndim*elemP%dof) )
    enddo

   call SetEdgeDegrees(gridL)

    do i=1,gridL%nelem
       elemP => gridL%elem(i)

       call SetElementQuadraturesDegrees(elemP )

       call PrepareOneElement(elemP )

       call ComputeIntegNode(elemP )

       call ComputeLocalMassMatrix(elemP )

       call InitMatrixElementShape(elemP )

       !write(*,*) 'i=',i,'Qdof = ',elemP%Qdof, ' icur =', elemP%deg_cur

    enddo

    ! assempling matrix block together
    call InitGlobalMatrixShape(gridL)

   ! setting of boundary conditions
   ! call SetConstBC(gridL)

    ! setting of the initial approximation
    !Qdof = gridL%elem(1)%Qdof  ! Qdof is the same for the whole patch
    elemP => gridL%elem(1)
    Fdof = max(elemP%Qdof, maxval( elemP%face(fGdof, :)  ) )

    allocate( xi(1:Fdof, 1:2) )
    allocate( wi(1:ndim, 1:Fdof) )
    allocate( w(1:ndim, 1:dofP) )
    allocate( Qphi(1:dofP, 1:Fdof) )

    ! we go over elements of the local grid
    do i=1, gridL%nelem
       elemP => gridL%elem(i)
       Qdof = elemP%Qdof


      ! setting of the patch basis functions in integ nodes of the element from the patch
       call BarycCoord(elem, Qdof, elemP%xi(0, 1:Qdof, 1:2),  xi(1:Qdof, 1:2) )

       call Eval_phi_Qnode(elem, dofP, Qdof, xi(1:Qdof, 1:nbDim), Qphi(1:dofP, 1:Qdof) )

       ! solution in integ nodes
       do k=1, ndim
          wi(k, 1:Qdof) = matmul( elem%w(0, (k-1)*dof +1: k*dof),  Qphi(1:dof, 1:Qdof) )

          call IntegrateVectorB( elemP, dofP, wi(k, 1:Qdof), w(k, 1:dofP) )

          elemP%w(0, (k-1)*dofP +1: k*dofP) = matmul(elemP%MassInv%Mb(1:dofP, 1:dofP),  w(k, 1:dofP) )
       enddo

       !call PlotElemSolution3D( 40, elemP)

       !do l=1,Qdof
       !   write(26, *) elemP%xi(0,l, 1:2), wi(:, l)
       !enddo

       ! setting of the BC
       allocate(elemP%wSS(1:elemP%flen, 1: elemP%face(fGdof, 1), 1:ndim ) ) ! the same degree

       do j=1, elemP%flen
          if(elemP%face(neigh, j) < 0) then
             Qdof = elemP%face(fGdof, j)


             ! setting of the patch basis functions in integ nodes of the element from the patch
             call BarycCoord(elem, Qdof, elemP%xi(j, 1:Qdof, 1:2),  xi(1:Qdof, 1:2) )

             call Eval_phi_Qnode(elem, dofP, Qdof, xi(1:Qdof, 1:nbDim), Qphi(1:dofP, 1:Qdof) )

             ! solution in integ nodes
             do k=1, ndim
                !wi(k, 1:Qdof) = matmul( elem%w(0, (k-1)*dof +1: k*dof),  Qphi(1:dof, 1:Qdof) )
                elemP%wSS(j, 1:Qdof, k) = matmul( elem%w(0, (k-1)*dof +1: k*dof),  Qphi(1:dof, 1:Qdof) )
             enddo

             !do l=1,Qdof
             !   write(50, *) elemP%xi(j,l, 1:2), elemP%wSS(j, l, 1:ndim)
             !enddo
             !write(50,*) '############'

          endif

       enddo

    enddo

    deallocate(xi, wi, w,Qphi)


  end subroutine Set_Local_Patch_interior



  !> setting of the local patch consisting of one element
  subroutine Set_Local_Patch_oneK(elem)
    class(element), intent(in) :: elem
    class(element), pointer :: elemP, elem1, elemP1
    integer, dimension (:,:), allocatable:: iloc
    integer :: i, j, j1, j2, k, l, l1, ib, flenP, flen
    integer :: j_, j1_, j2_
    integer :: dofP, dof, Qdof,Fdof


    dof  = elem%dof

    flenP = elem%flen
    if(flenP /= 3) stop 'non-implemented cas 4dxxdc5etd76e in local_problem2.f90'


    ! allocation of local mesh for the reconstruction
    call gridL%init(1, 3, 3, 1, 2)

    ! creation of the local grid
    ! nodes of element elem
    gridL%x(1, :) = grid%x(elem%face(idx, 1), : )
    gridL%x(2, :) = grid%x(elem%face(idx, 2), : )
    gridL%x(3, :) = grid%x(elem%face(idx, 3), : )


    !! vertexes of triangles of local mesh, idx = 4 index of the mother element
    allocate(iloc(1:gridL%nelem, 1:4) )
    iloc(1, 1:3) = (/ 1, 2, 3 /)
    iloc(1, 4) = elem%i

    i  = 1  ! index of triagles
    ib = 0  ! index of boundary segment


    do j=1, flenP
       j1 = mod(j , 3) + 1
       j2 = mod(j1, 3) + 1

       k = elem%face(neigh, j)
       if(k > 0) then
          elem1 => grid%elem(k)

          !j_  =  elem%face(nei_i, j)
          !j1_ = mod(j_  , 3) + 1
          !j2_ = mod(j1_ , 3) + 1

          ! two boundary segment
          ib = ib + 1
          gridL%b_edge(ib)%lbn(1) = j
          gridL%b_edge(ib)%lbn(2) = j1
          gridL%b_edge(ib)%ibc = k

       else ! one boundary segment
          ib = ib + 1
          gridL%b_edge(ib)%lbn(1) = j
          gridL%b_edge(ib)%lbn(2) = j1
          gridL%b_edge(ib)%ibc = -1

       endif
    enddo
    gridL%b_edge(1:gridL%nbelm)%icurv = 0


    do i = 1, gridL%nelem
       elemP => gridL%elem(i)
       call elemP%init(i , 3, 3, iloc(i, 1:3) )

       elemP%psplit = iloc(i,4)  ! storing of the index of mother elements
       !elemP%deg_cur = 1  !!!!
       !elemP%ibcur = -1  !!!!
       !elemP%jcur = 0  !!!!
       !elemP%iBC(:) = -1  ! = 0 means Neuman BC, > 1 Dirichlet taken from *.ini, -1 Dirichlet(probably?)
    enddo


    ! do i=1,gridL%npoin
    !    write(*,'(a6, i8, 2es15.6)') 'xp:',i, gridL%x(i,:)
    ! enddo
    ! do i=1,gridL%nelem
    !    write(*,'(a6,30i8)') 'elemL:',i, gridL%elem(i)%face(idx, :), gridL%elem(i)%face(neigh, :)
    ! enddo

    ! do i=1,gridL%nbelm
    !    write(*,'(a6,30i8)') 'B_elL:',i, gridL%b_edge(i)%lbn(1:2), gridL%b_edge(i)%ibc
    ! enddo

    call gridL%plot('meshL')

    call gridL%seekNeighbours()

    ! setting of further arrays
    do i = 1, gridL%nelem
       elemP => gridL%elem(i)

       elemP%deg_cur = 1  !!!!
       elemP%ibcur = -1  !!!!
       elemP%jcur = 0  !!!!
       elemP%iBC(:) = -1  ! = 0 means Neuman BC, > 1 Dirichlet taken from *.ini, -1 Dirichlet(probably?)
    enddo

    !write(*,'(a6,30i5)') 'elem:',elem%i, elem%face(idx, :), elem%face(neigh, :)

    ! do i=1,gridL%nelem
    !    write(*,'(a6,30i5)') 'elemL:',i, gridL%elem(i)%face(idx, :), gridL%elem(i)%face(neigh, :)
    ! enddo

    ! do i=1,gridL%nbelm
    !    write(*,'(a6,30i5)') 'B_elL:',i, gridL%b_edge(i)%lbn(1:2), gridL%b_edge(i)%ibc
    ! enddo

   ! setting of boundary conditions
   ! call SetConstBC(gridL)

    deallocate(iloc)

  end subroutine Set_Local_Patch_oneK

  !> setting of the local patch consisting of four neighbouring elements
  subroutine Set_Local_Patch_neighbours(elem)
    class(element), intent(in) :: elem
    class(element), pointer :: elemP, elem1, elemP1
    integer, dimension (:,:), allocatable:: iloc
    integer :: i, j, j1, j2, k, l, l1, ib, flenP, flen
    integer :: j_, j1_, j2_
    integer :: dofP, dof, Qdof,Fdof


    dof  = elem%dof

    flenP = elem%flen
    if(flenP /= 3) stop 'non-implemented cas 4dxxdc5etd76e in local_problem2.f90'

    ! counting of the number of neighboursm then number of elements is  i + 1,
    i = 0
    do j=1, flenP
       if(elem%face(neigh, j) > 0) i = i + 1
    enddo

    ! number nodes is 3 + i
    ! number boundary segments is 3 + i

    ! allocation of local mesh for the reconstruction
    call gridL%init(i + 1, i + 3, 3 + i, 1, 2)

    ! creation of the local grid
    ! nodes of element elem
    gridL%x(1, :) = grid%x(elem%face(idx, 1), : )
    gridL%x(2, :) = grid%x(elem%face(idx, 2), : )
    gridL%x(3, :) = grid%x(elem%face(idx, 3), : )

    ! the nodes of neighbouring elements
    i = 3
    do j=1, flenP
       k = elem%face(neigh, j)
       if(k > 0) then
          elem1 => grid%elem(k)
          i = i + 1
          j1 = mod(elem%face(nei_i, j) , 3) + 1
          j2 = mod(j1, 3) + 1
          gridL%x(i, :) =  grid%x(elem1%face(idx, j2), : )
       endif

    enddo

    !! vertexes of triangles of local mesh, idx = 4 index of the mother element
    allocate(iloc(1:gridL%nelem, 1:4) )
    iloc(1, 1:3) = (/ 1, 2, 3 /)
    iloc(1, 4) = elem%i

    i  = 1  ! index of triagles
    ib = 0  ! index of boundary segment


    do j=1, flenP
       j1 = mod(j , 3) + 1
       j2 = mod(j1, 3) + 1

       k = elem%face(neigh, j)
       if(k > 0) then
          elem1 => grid%elem(k)
          i = i + 1
          iloc(i, 1) = j1
          iloc(i, 2) = j
          iloc(i, 3) = i + 2  ! index of the node corresponding to the vertex of neighbouring element
          iloc(i, 4) = elem1%i   ! index of the mother element


          j_  =  elem%face(nei_i, j)
          j1_ = mod(j_  , 3) + 1
          j2_ = mod(j1_ , 3) + 1

          ! two boundary segment
          ib = ib + 1
          gridL%b_edge(ib)%lbn(1) = j
          gridL%b_edge(ib)%lbn(2) = i + 2
          gridL%b_edge(ib)%ibc = 1

          gridL%b_edge(ib)%lbn(3) = elem1%i                ! global index of the neighbour
          gridL%b_edge(ib)%lbn(4) = elem1%face(neigh, j1_) ! global index of the neigh of the neigh

          ib = ib + 1
          gridL%b_edge(ib)%lbn(1) = i + 2
          gridL%b_edge(ib)%lbn(2) = j1
          gridL%b_edge(ib)%ibc = 1

          gridL%b_edge(ib)%lbn(3) = elem1%i                ! global index of the neighbour
          gridL%b_edge(ib)%lbn(4) = elem1%face(neigh, j2_) ! global index of the neigh of the neigh


       else ! one boundary segment
          ib = ib + 1
          gridL%b_edge(ib)%lbn(1) = j
          gridL%b_edge(ib)%lbn(2) = j1
          gridL%b_edge(ib)%ibc = -1

          gridL%b_edge(ib)%lbn(3) = elem%i
          gridL%b_edge(ib)%lbn(4) = -1

       endif
    enddo
    gridL%b_edge(1:gridL%nbelm)%icurv = 0


    do i = 1, gridL%nelem
       elemP => gridL%elem(i)
       call elemP%init(i , 3, 3, iloc(i, 1:3) )

       elemP%psplit = iloc(i,4)  ! storing of the index of mother elements
       !elemP%deg_cur = 1  !!!!
       !elemP%ibcur = -1  !!!!
       !elemP%jcur = 0  !!!!
       !elemP%iBC(:) = -1  ! = 0 means Neuman BC, > 1 Dirichlet taken from *.ini, -1 Dirichlet(probably?)
    enddo


    ! do i=1,gridL%npoin
    !    write(*,'(a6, i8, 2es15.6)') 'xp:',i, gridL%x(i,:)
    ! enddo
    ! do i=1,gridL%nelem
    !    write(*,'(a6,30i8)') 'elemL:',i, gridL%elem(i)%face(idx, :), gridL%elem(i)%face(neigh, :)
    ! enddo

    ! do i=1,gridL%nbelm
    !    write(*,'(a6,30i8)') 'B_elL:',i, gridL%b_edge(i)%lbn(1:2), gridL%b_edge(i)%ibc
    ! enddo

    call gridL%plot('meshL')

    call gridL%seekNeighbours()

    ! setting of further arrays
    do i = 1, gridL%nelem
       elemP => gridL%elem(i)

       elemP%deg_cur = 1  !!!!
       elemP%ibcur = -1  !!!!
       elemP%jcur = 0  !!!!
       elemP%iBC(:) = -1  ! = 0 means Neuman BC, > 1 Dirichlet taken from *.ini, -1 Dirichlet(probably?)
    enddo

    !write(*,'(a6,30i5)') 'elem:',elem%i, elem%face(idx, :), elem%face(neigh, :)

    ! do i=1,gridL%nelem
    !    write(*,'(a6,30i5)') 'elemL:',i, gridL%elem(i)%face(idx, :), gridL%elem(i)%face(neigh, :)
    ! enddo

    ! do i=1,gridL%nbelm
    !    write(*,'(a6,30i5)') 'B_elL:',i, gridL%b_edge(i)%lbn(1:2), gridL%b_edge(i)%ibc
    ! enddo

   ! setting of boundary conditions
   ! call SetConstBC(gridL)

    deallocate(iloc)

  end subroutine Set_Local_Patch_neighbours

  subroutine SolveLocalNewton2( ) !gridL )
    !class (mesh),  intent(inout), target :: gridL
    integer :: flenP
    class(element), pointer :: elemP
    real, dimension(:), allocatable :: x, b, b1  ! vectors for newton
    real, dimension(:,:), allocatable :: mA
    !logical :: vector_update
    real :: res0, res1
    real :: Newton_tol, Newton_max_iter, theta, lambda, lambda_old
    integer :: iter, ndof, dof
    integer :: i, j, k, l, is


    Newton_tol = 1E-8
    Newton_max_iter = 10

    ndof = sum(gridL%elem(:)%dof) * ndim
    !print*,'ndim = ', ndof

    elemP => gridL%elem(1)

    ! arrays for Newton
    allocate(x(1:ndof), b(1:ndof), b1(1:ndof), mA(1:ndof, 1:ndof)  )
    mA(:,:) = 0.
    x(:) = 0.

    state%nlSolver%implicitly = .true.
    ! Newton iterations
    do iter = 1, Newton_max_iter
       !print*,'###############################    Newton_iter', iter,'       elem=',elemP%i

       if( state%nlSolver%implicitly ) then
          call ComputeLocalBlocks2( ) ! matrix block
          !print*,'Matrix updated'
       endif

       state%nlSolver%implicitly = .false.
       !print*,'------------------------**************-----------------------------------------------'

       if(iter == 1) then
          call ComputeLocalBlocks2( ) ! vector

          !do i=1,ndof
          !   write(*,'(i5,a2,500es9.1)') i,': ',mA(i, :)
          !enddo

          ! filling of the vector
          call FillLocalVector(ndof, b(1:ndof) )
          !b(1:ndof) = elemP%vec(rhs, 1:ndof)
          res0 = VectorNorm(b(1:ndof))

       else
          b(1:ndof) = b1(1:ndof)
          res0 = res1
       endif


       !if(iter == 1 ) then
       !   print*,'----------------------------------------'
       !   write(*,'(a6,i5,3es12.4,a3,300es12.4)') &
       !        'iters:',0, res0,-2., -2.,'|',gridL%elem(1)%w(0,1:6)
       !endif

       x(1:ndof) = b(1:ndof)


       ! print*,'--------------------------------------------------'
       ! do j = 1,ndof
       !    write(*,'(a3,i5,300es12.4)') 'A',j, mA(j,:), x(j)
       ! enddo
       !print*,'--------------------------------------------------',iter, res0, res1

       call FillLocalMatrix(ndof, mA(1:ndof, 1:ndof) ) ! matrix mA is changed


       call SolveLocalMatrixProblem(ndof, mA(1:ndof, 1:ndof), 1, x(1:ndof) )



       lambda_old = 0.
       lambda = 1.0   ! initialization of lambda (= damping factor

       do l=1,5    ! iterations, seeking the optimal damping factor

          call UpdateLocalSolution(ndof, x(1:ndof), lambda, lambda_old)

          lambda_old = lambda

          state%nlSolver%implicitly = .false.
          call ComputeLocalBlocks2( ) !gridL )  ! vector
          call FillLocalVector(ndof, b1(1:ndof) )

          res1 = VectorNorm(b1(1:ndof))

          theta = res1 / res0
          !write(*,'(a6,i5,3es12.4,a3,300es12.4)') &
          !    'iters:',iter, res1,theta, lambda,'|', gridL%elem(1)%w(0,1:6) !gridL%elem(1)%vec(rhs,1:6)

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

  end subroutine SolveLocalNewton2


  !> compute either a local matrix block or a block vector
  subroutine ComputeLocalBlocks2( ) !gridL )
    !class (mesh),  intent(inout), target :: gridL
    integer :: flenP
    class(element), pointer :: elemP, elemP0
    integer :: i, j, k


    do i=1, gridL%nelem
       elemP => gridL%elem(i)

       if(state%modelName == 'scalar') then        ! 2D scalar equation
          call ComputeLocalMatrixBlock2(Set_f_s_scalar, Set_A_s_scalar, Set_Ppm_scalar, &
               Set_R_s_scalar, Set_K_sk_scalar, Set_S_scalar, Set_DS_scalar, elemP)
       else
          stop 'Other variants in SolveLocalProblem not YET implemented in local_problem.f90'
       endif

       do k=1, elemP%flen
          if( elemP%face(neigh,k) > 0) then
             !print*,'#####',elemP%i,k
             !do j=1,elemP%dof
             !   write(*,'(a3,2i5,300es12.4)') 'WDKN',i,j,elemP%vec(RHS,j), elemP%block(k)%Mb(j, :)
             !enddo
             !print*
          endif
       enddo
       !print*,'### ed49 ##################################', i, flenP
    enddo



  end subroutine ComputeLocalBlocks2



  !> preparation of the element for the local problem
  subroutine PrepareLocalElement2(elemP, degP, elem)
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

  end subroutine PrepareLocalElement2

  !> deallocation of the element
  subroutine DeleteLocalElement2(elemP)
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

  end subroutine DeleteLocalElement2


  !> compute the matrix block and the RHS for the local problem
  subroutine ComputeLocalMatrixBlock2(Set_f_s, Set_A_s, Set_Ppm, Set_R_s, Set_K_sk, Set_S, Set_DS, elemP)
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

    elemP%vec(RHS,:) = 0.

    call ComputeOneElementTerms(Set_f_s, Set_A_s, Set_Ppm, Set_R_s, Set_K_sk, Set_S, Set_DS, elemP)


    if(elemP%i == -1) then
       write(*,*) 'elemP%dof = ', elemP%dof
       do i=1,elemP%dof
          write(*,'(a3,i5,300es12.4)') 'WDKN',i,elemP%vec(RHS,i), elemP%block(0)%Mb(i, :)
          !write(*,'(a3,i5,300es12.4)') 'WDKN',i, elemP%block(0)%Mb(i, :)
       enddo
       print*,'))))))))))))))))))))))))))))))))))))))))))))))))'
    endif

  end subroutine ComputeLocalMatrixBlock2

  !> solve the Polynomially Preserving local problem
  subroutine SolvePPlocalProblem2( )
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

       call ComputeLocalBlocks2( )  ! matrix block

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
  end subroutine SolvePPlocalProblem2

  !> assembling of the matrix of the local  problem on the patch
  subroutine  FillLocalMatrix(nsize, A )
    integer, intent(in) :: nsize
    real, dimension(1:nsize, 1:nsize), intent(inout) :: A
    class(element), pointer :: elem, elem1
    integer :: i, k, in, ir, ir1, ic, ic1, ndof1, ndof

    do i=1,gridL%nelem
       elem => gridL%elem(i)
       ndof = elem%dof * ndim

       ir = elem%ncv
       ir1 = ir + ndof - 1

       A(ir:ir1, ir:ir1 ) = elem%block(0)%Mb(1:ndof,1:ndof)

       do k=1,elem%flen
          in = elem%face(neigh,k)
          if(in >0) then
             elem1 => gridL%elem(in)
             ndof1 = elem1%dof * ndim

             ic = elem1%ncv
             ic1 = ic + ndof1 - 1

             A(ir:ir1, ic:ic1 ) = elem%block(k)%Mb(1:ndof,1:ndof1)
          endif
       enddo
    enddo


  end subroutine FillLocalMatrix


  !> assembling of the RHS of the local  problem on the patch
  subroutine  FillLocalVector(nsize, b )
    integer, intent(in) :: nsize
    real, dimension(1:nsize), intent(inout) :: b
    class(element), pointer :: elem
    integer :: i,  ndof, is

    is = 1
    do i=1,gridL%nelem
       elem => gridL%elem(i)
       ndof = elem%dof * ndim

       !print*,'RFT^%$EDR;', is, elem%ncv

       b(is : is + ndof -1) = elem%vec(rhs, 1:ndof)
       is = is + ndof
    enddo
  end subroutine FillLocalVector


  !> updating of the solution after the Newton iteration
  !> lambda, lambda_old are damping parameters
  subroutine UpdateLocalSolution(nsize, x, lambda, lambda_old)
    integer, intent(in) :: nsize
    real, dimension(1:nsize), intent(inout) :: x
    real, intent(in) :: lambda, lambda_old
    class(element), pointer :: elemP
    integer :: i,  dof, is
    real :: update

    update = lambda - lambda_old

    do i=1,gridL%nelem
       elemP => gridL%elem(i)
       is  = elemP%ncv
       dof = elemP%dof * ndim
       elemP%w(0,1:dof) = elemP%w(0,1:dof) + update * x(is : is + dof - 1)
    enddo

  end subroutine UpdateLocalSolution


  subroutine Prepare_Local_grid( degP )
    integer, intent(in) :: degP
    class(element), pointer :: elemP
    integer :: i

    call gridL%computeGeometry()

    do i=1,gridL%nelem
       elemP => gridL%elem(i)
       elemP%deg = degP
       elemP%Tdeg = 0
       elemP%degP = 0

       call elemP%initElementDof()


       allocate(elemP%w(0:1, 1:ndim*elemP%dof) )
    enddo

   call SetEdgeDegrees(gridL)

    do i=1,gridL%nelem
       elemP => gridL%elem(i)

       call SetElementQuadraturesDegrees(elemP )

       call PrepareOneElement(elemP )

       call ComputeIntegNode(elemP )

       call ComputeLocalMassMatrix(elemP )

       call InitMatrixElementShape(elemP )

       !write(*,*) 'i=',i,'Qdof = ',elemP%Qdof, ' icur =', elemP%deg_cur

    enddo

    ! assempling matrix block together
    call InitGlobalMatrixShape(gridL)


  end subroutine Prepare_Local_grid

  !> initiate the solution on the local grid and setting of the BC
  subroutine Prepare_Local_grid_solution( degP )
    integer, intent(in) :: degP
    class(element), pointer :: elemP, elem1, elem2
    real, dimension (:,:), allocatable:: xi, wi, Qphi,w
    integer :: dof, dofP, Fdof, Qdof
    integer :: i, j, k, l, ib


    dofP = DOFtriang(degP)

    ! setting of the initial approximation
    elemP => gridL%elem(1)
    Fdof = max(elemP%Qdof, maxval( elemP%face(fGdof, :)  ) )

    allocate( xi(1:Fdof, 1:2) )
    allocate( wi(1:ndim, 1:Fdof) )
    allocate( w(1:ndim, 1:dofP) )
    allocate( Qphi(1:dofP, 1:Fdof) )

    ! we go over elements of the local grid
    do i=1, gridL%nelem
       elemP => gridL%elem(i)
       Qdof = elemP%Qdof

       elem1 => grid%elem( elemP%psplit )  ! pointer to the mother element
       dof = elem1%dof

      ! setting of the patch basis functions in integ nodes of the element from the patch
       call BarycCoord(elem1, Qdof, elemP%xi(0, 1:Qdof, 1:2),  xi(1:Qdof, 1:2) )

       call Eval_phi_Qnode(elem1, dofP, Qdof, xi(1:Qdof, 1:nbDim), Qphi(1:dofP, 1:Qdof) )

       ! solution in integ nodes
       do k=1, ndim
          wi(k, 1:Qdof) = matmul( elem1%w(0, (k-1)*dof +1: k*dof),  Qphi(1:dof, 1:Qdof) )

          call IntegrateVectorB( elemP, dofP, wi(k, 1:Qdof), w(k, 1:dofP) )

          elemP%w(0, (k-1)*dofP +1: k*dofP) = matmul(elemP%MassInv%Mb(1:dofP, 1:dofP),  w(k, 1:dofP) )
       enddo

       !call PlotElemSolution3D( 40, elemP)

       !do l=1,Qdof
       !   write(26, *) elemP%xi(0,l, 1:2), wi(:, l)
       !enddo

       ! for the setting of the BC
       allocate(elemP%wSS(1:elemP%flen, 1: elemP%face(fGdof, 1), 1:ndim ) ) ! the same degree

       ! do j=1, elemP%flen
       !    if(elemP%face(neigh, j) < 0) then
       !       Qdof = elemP%face(fGdof, j)


       !       ! setting of the patch basis functions in integ nodes of the element from the patch
       !       call BarycCoord(elem1, Qdof, elemP%xi(j, 1:Qdof, 1:2),  xi(1:Qdof, 1:2) )

       !       call Eval_phi_Qnode(elem1, dofP, Qdof, xi(1:Qdof, 1:nbDim), Qphi(1:dofP, 1:Qdof) )

       !       ! solution in integ nodes
       !       do k=1, ndim
       !          !wi(k, 1:Qdof) = matmul( elem%w(0, (k-1)*dof +1: k*dof),  Qphi(1:dof, 1:Qdof) )
       !          elemP%wSS(j, 1:Qdof, k) = matmul( elem1%w(0, (k-1)*dof +1: k*dof),  Qphi(1:dof, 1:Qdof) )
       !       enddo

       !       do l=1,Qdof
       !          write(50, *) elemP%xi(j,l, 1:2), elemP%wSS(j, l, 1:ndim)
       !       enddo
       !       write(50,*) '############'

       !    endif

       ! enddo

    enddo


    ! setting of the "BC" for all boundary segments  of the local grid
    do ib =1, gridL%nbelm
       elemP => gridL%elem(gridL%b_edge(ib )%itc )
       j = gridL%b_edge(ib)%jtc
       Qdof = elemP%face(fGdof, j)

       !Qdof = elemP%Qdof

       !write(*,'(a10,30i7)') ':::::',ib,elemP%i, gridL%b_edge(ib)%ibc, j, gridL%b_edge(ib)%lbn(:)

       if( gridL%b_edge(ib)%lbn(4) < 0) then ! edge on the boundary, we used Dirichlet BC

          if(state%modelName == 'scalar' .or.state%modelName == '2eqs' ) then
             do k=1,Qdof
                call Exact_Scalar(elemP%xi(j, k, 1:nbDim), elemP%wSS(j,k,1:ndim), state%time%ctime )
             enddo

          else
             stop 'NON-implemeted edr54ede45sfe54dre in local_problem.f90'
          endif

          ! FOR NEUMANN BC it is necessary to use extrapolation !!!!
       else
          !setting of BC from the average

          elem1 =>  grid%elem(gridL%b_edge(ib )%lbn(3) )   ! global elements sharing the edge
          elem2 =>  grid%elem(gridL%b_edge(ib )%lbn(4) )


          !computing the average

          ! setting of the patch basis functions in integ nodes of the element from the patch
          call BarycCoord(elem1, Qdof, elemP%xi(j, 1:Qdof, 1:2),  xi(1:Qdof, 1:2) )

          call Eval_phi_Qnode(elem1, dofP, Qdof, xi(1:Qdof, 1:nbDim), Qphi(1:dofP, 1:Qdof) )

          ! solution in integ nodes
          do k=1, ndim
             !wi(k, 1:Qdof) = matmul( elem%w(0, (k-1)*dof +1: k*dof),  Qphi(1:dof, 1:Qdof) )
             elemP%wSS(j, 1:Qdof, k) = matmul( elem1%w(0, (k-1)*dof +1: k*dof),  Qphi(1:dof, 1:Qdof) )
          enddo

          ! setting of the patch basis functions in integ nodes of the element from the patch
          call BarycCoord(elem2, Qdof, elemP%xi(j, 1:Qdof, 1:2),  xi(1:Qdof, 1:2) )

          call Eval_phi_Qnode(elem2, dofP, Qdof, xi(1:Qdof, 1:nbDim), Qphi(1:dofP, 1:Qdof) )

          ! solution in integ nodes
          do k=1, ndim
             !wi(k, 1:Qdof) = matmul( elem%w(0, (k-1)*dof +1: k*dof),  Qphi(1:dof, 1:Qdof) )
             elemP%wSS(j, 1:Qdof, k) = (elemP%wSS(j, 1:Qdof, k) &
                  + matmul( elem2%w(0, (k-1)*dof +1: k*dof),  Qphi(1:dof, 1:Qdof) )) / 2
          enddo

       endif

       do l=1,Qdof
          write(50, *) elemP%xi(j,l, 1:2), elemP%wSS(j, l, 1:ndim)
       enddo
       write(50,*) '############'

    enddo  ! ib =1,gridL%nbelm



    deallocate(xi, wi, w,Qphi)




  end subroutine Prepare_Local_grid_solution



  !> initiate the solution on the local grid and setting of the BC
  subroutine Prepare_Local_grid_solution_Rhw( degP, dofP, Rhw )
    integer, intent(in) :: degP, dofP
    real, dimension(1:ndim, 1:dofP), intent(inout) :: Rhw
    class(element), pointer :: elemP, elem1, elem2
    real, dimension (:,:), allocatable:: xi, wi, Qphi,w
    integer :: dof, Fdof, Qdof
    integer :: i, j, k, l, ib


    ! setting of the initial approximation
    elemP => gridL%elem(1)
    Fdof = max(elemP%Qdof, maxval( elemP%face(fGdof, :)  ) )

    ! we go over elements of the local grid
    do i=1, gridL%nelem
       elemP => gridL%elem(i)

       do k=1, ndim
          elemP%w(0, (k-1)*dofP +1: k*dofP) = Rhw(k, 1:dofP)
       enddo

       !call PlotElemSolution3D( 40, elemP)

       ! for the setting of the BC
       allocate(elemP%wSS(1:elemP%flen, 1: elemP%face(fGdof, 1), 1:ndim ) ) ! the same degree
       do j=1, elemP%flen
          Qdof =  elemP%face(fGdof, 1)

          call Eval_w_Edge(elemP, j, elemP%wSS(j, 1:Qdof, 1:ndim), .false.)
          !do l=1,Qdof
          !   write(50, *) elemP%xi(j,l, 1:2), elemP%wSS(j, l, 1:ndim)
          !enddo
          !write(50,*) '############'

          ib =  -elemP%face(neigh, j) ! index of the corresponding boundary segment

          if(gridL%b_edge(ib)%ibc <= 0) then
             ! Dirichlet BC, we set the BC
             if(state%modelName == 'scalar' .or.state%modelName == '2eqs' ) then
                do k=1,Qdof
                   call Exact_Scalar(elemP%xi(j, k, 1:nbDim), elemP%wSS(j,k,1:ndim), state%time%ctime )
                enddo

             else
                stop 'NON-implemeted edr54ede45sfe54dre in local_problem.f90'
             endif
          endif

          !do l=1,Qdof
          !   write(55, *) elemP%xi(j,l, 1:2), elemP%wSS(j, l, 1:ndim)
          !enddo
          !write(55,*) '############'

       enddo


    enddo



    ! setting of the "BC" for all boundary segments  of the local grid
    ! do ib =1, gridL%nbelm
    !    elemP => gridL%elem(gridL%b_edge(ib )%itc )
    !    j = gridL%b_edge(ib)%jtc
    !    Qdof = elemP%face(fGdof, j)

    !    !Qdof = elemP%Qdof

    !    !write(*,'(a10,30i7)') ':::::',ib,elemP%i, gridL%b_edge(ib)%ibc, j, gridL%b_edge(ib)%lbn(:)

    !    if( gridL%b_edge(ib)%lbn(4) < 0) then ! edge on the boundary, we used Dirichlet BC

    !       if(state%modelName == 'scalar' .or.state%modelName == '2eqs' ) then
    !          do k=1,Qdof
    !             call Exact_Scalar(elemP%xi(j, k, 1:nbDim), elemP%wSS(j,k,1:ndim), state%time%ctime )
    !          enddo

    !       else
    !          stop 'NON-implemeted edr54ede45sfe54dre in local_problem.f90'
    !       endif

    !       ! FOR NEUMANN BC it is necessary to use extrapolation !!!!
    !    else
    !       !setting of BC from the average

    !       call Eval_w_Edge(elemP, j, elemP%wSS(j, 1:Qdof, 1:ndim), .false.)


    !    endif

    !    do l=1,Qdof
    !       write(50, *) elemP%xi(j,l, 1:2), elemP%wSS(j, l, 1:ndim)
    !    enddo
    !    write(50,*) '############'

    ! enddo  ! ib =1,gridL%nbelm



    !deallocate(xi, wi, w,Qphi)

  end subroutine Prepare_Local_grid_solution_Rhw

end module local_problem2

