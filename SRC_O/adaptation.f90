!> mesh adaptation
module mesh_adaptation
  use main_data  ! contains "type(mesh) ::  grid"   for computation
  use ama_L2interpol
  use geometry
  use mesh_oper
  use problem_oper
  use set_solution
  use stdgm_mod

  public:: SaveOrigMesh
  public:: RemoveOrigMesh
  !public:: AdaptMesh_AMA
  public:: ReprepareProblem
  public:: RecomputeElemSol
  public:: RecomputeElemSolDerefine
  public:: PassElem2Elem
  public:: DeallocateElem
!  public:: Mesh2Amesh
!  public:: Amesh2Mesh
!  public:: PlotAmesh1

contains

  !> used for the adaptation of time dependent problems
  !> save the actual mesh 'grid' into 'gridS' for a possible interpolation
  !> of an unsucessful steps
  subroutine SaveOrigMesh( )
    integer :: i, Td, k
    class(element), pointer :: elem, elemS

    state%time%ttime_save = state%time%ttime
    state%timeprn_save = state%timeprn
    state%isol_save = state%isol

    state%space%gridS_allocated = .true.

    allocate (gridS)

    gridS%nelem = grid%nelem
    gridS%npoin = grid%npoin

    allocate(gridS%x( 1:gridS%npoin,1:nbDim) )
    gridS%x( 1:gridS%npoin,1:nbDim) = grid%x( 1:gridS%npoin,1:nbDim)

    allocate(gridS%elem( 1:gridS%nelem) )

    do i=1, gridS%nelem
       elem  => grid%elem(i)
       elemS => gridS%elem(i)

       elemS%i    = elem%i
       elemS%flen = elem%flen
       elemS%deg  = elem%deg
       elemS%dof  = elem%dof
       elemS%Tdeg  = elem%Tdeg
       elemS%Tdof  = elem%Tdof
       elemS%TQnum  = elem%TQnum
       elemS%Qdeg = elem%Qdeg
       elemS%Qnum = elem%Qnum

       allocate(elemS%face(1:nei_i, 1:elemS%flen))
       elemS%face(idx:nei_i, 1:elemS%flen) = elem%face(idx:nei_i, 1:elemS%flen)

       allocate(elemS%xc(1:2) )
       elemS%xc(1:2) = elem%xc(1:2)
       elemS%area = elem%area

       elemS%deg_cur = elem%deg_cur

       allocate(elemS%F)
       elemS%F%deg = elem%F%deg
       elemS%F%dof = elem%F%dof

       allocate(elemS%F%F(1:elem%F%dof,1:nbDim))
       elemS%F%F(1:elem%F%dof, 1:nbDim) = elem%F%F(1:elem%F%dof, 1:nbDim)

       elemS%F%iFlin = elem%F%iFlin

       if(elemS%F%iFlin ) then
          allocate(elemS%F%D1F0(1:nbDim,1:nbDim) )

          elemS%F%JF0 = elem%F%JF0
          elemS%F%D1F0(1:2,1:2) = elem%F%D1F0(1:2,1:2)
       endif

       elemS%HGnode = elem%HGnode
       if(elemS%HGnode) then
          allocate(elemS%HGvertex(1:3))
          elemS%HGvertex(1:3) = elem%HGvertex(1:3)
       endif

    enddo

    if(state%time%disc_time == 'STDG') then
       !print*,' STDGM for SaveOrigMesh in adaptation.f90 not implemented'

       Td = 1
       do i=1, gridS%nelem
          elem  => grid%elem(i)
          elemS => gridS%elem(i)

          allocate(elemS%w(0:Td, 1:ndim*elem%dof) )
          ! store w at last time level, not the whole wST !
          !
          call Transfer_wST_to_w_Elem(elem , 0, elem%TQnum)
          elemS%w(0:Td, 1:ndim*elem%dof) = elem%w(0:Td, 1:ndim*elem%dof)

          !do k=1,ndim
          !   !print*,'###',k,elem%dof, size(elem%wSTfinAD, 1), size(elem%wSTfinAD,2)
          !   elemS%w(0, (k-1)*elem%dof + 1: k*elem%dof) = elem%wSTfinAD(k,1:elem%dof)
          !enddo

       enddo

    else

       Td = state%time%deg+1
       do i=1, gridS%nelem
          elem  => grid%elem(i)
          elemS => gridS%elem(i)

          allocate(elemS%w(0:Td, 1:ndim*elem%dof) )
          elemS%w(0:Td, 1:ndim*elem%dof) = elem%w(0:Td, 1:ndim*elem%dof)

       enddo
    endif

  end subroutine SaveOrigMesh


  !> remove mesh gridS used for the adaptation of time dependent problems
  subroutine RemoveOrigMesh( )
    integer :: i
    class(element), pointer :: elemS

    state%space%gridS_allocated = .false.

    deallocate(gridS%x)

    do i=1, gridS%nelem
       elemS => gridS%elem(i)
       deallocate(elemS%face, elemS%xc )

       if(elemS%F%iFlin ) deallocate(elemS%F%D1F0 )
       deallocate(elemS%F%F)
       deallocate(elemS%F)

       if(elemS%HGnode)  deallocate(elemS%HGvertex )
    enddo

    do i=1, gridS%nelem
       deallocate(gridS%elem(i)%w  )
    enddo

    deallocate(gridS%elem  )
    deallocate (gridS)

  end subroutine RemoveOrigMesh


  !> repreparation of the problem, outgoing grid% structure has to be same as
  !> in subroutine PrepareProblem in problem.f90, arrays and values from
  !> no refined elements are reused
  !>
  !> passing of datas from the old "type(mesh) :: gridO" to the new "type(mesh) :: gridN"
  !> elem%dof, elem%deg, .....
  !> deleting the elements of the old grid
  subroutine ReprepareProblem(gridN, gridO)
    type(mesh), intent(inout) :: gridN
    type(mesh), intent(inout) :: gridO
    class(element), pointer :: elem
    logical :: SimplyPass
    integer :: i,j, iold, ib

    !! general datas
    gridN%nbc = gridO%nbc

    gridO%elem(:)%dealloc = .true.

    ! for AMA done in Amesh2Mesh
    if(state%space%adapt%adapt_type /= 'AMA' .and. state%space%adapt%adapt_type /= 'ANI' &
         .and. state%space%adapt%adapt_type /= 'Ahp'.and. state%space%adapt%adapt_type /= 'Ihp') then
       !! setting of degrees of polynomial approximations
       do i=1,gridN%nelem
          elem => gridN%elem(i)
          elem%deg = elem%deg + elem%psplit  ! p-adaptation
          elem%Tdeg = state%time%deg

!          elem%TQnum = state%time%max_Tdof    !temporarily
          elem%TQnum = state%time%Qnum

          call elem%initElementDof()
       enddo
    endif

    state%space%max_dof = maxval(gridN%elem(:)%dof)

    !print*,'Setting of DOF for edges (used in the numerical quadratures)'
    call SetEdgeDegrees(gridN)

    !print *,'@@ 1'

    if ( .not. allocated( gridN%x ) ) stop 'x not alloc in ReprepareProblem'


    do i=1,gridN%nelem
       elem => gridN%elem(i)
       SimplyPass = .true.

       if(state%space%adapt%adapt_type == 'AMA' .or. state%space%adapt%adapt_type == 'ANI' &
            .or. state%space%adapt%adapt_type == 'Ahp' .or. state%space%adapt%adapt_type == 'Ihp'&
            .or. elem%RGtype /= 'N' .or. elem%psplit /= 0 ) then
          SimplyPass = .false.
       else
          do j=1,elem%flen
             if(elem%face(neigh,j) > 0) then
                if(gridN%elem(elem%face(neigh,j))%RGtype /='N' .or. &
                     gridN%elem(elem%face(neigh,j))%psplit /= 0  )  then
                   SimplyPass = .false.
                   goto 10
                endif
             endif
          enddo
       endif

10     continue
       !print*,'-----------------------------'
       !print *,'@@ 2', elem%i, SimplyPass

       !print*,'elem =',i,gridN%nelem,'  SimplyPass =',SimplyPass, &
       !     "elem%RGtype = ",elem%RGtype
       if(SimplyPass) then
          ! geometry of the element and its neighbours are the same,
          ! we pass simply all the arguments
          ! polynomial degrees of neighbours may differ
          call PassElem2Elem(elem, gridO%elem(elem%i))

          !print *,'@@ 3', elem%deg, elem%w(0,:)

          gridO%elem(elem%i)%dealloc = .false.
          elem%to_recompute = .false.
       else
          !print*,' new element', elem%i, elem%deg
          !call ComputeElementGeometry(gridN, elem)

          ! curved element
          if ( elem%deg_cur > 1 ) then
             print*, ' FR CURVED ELEMENTS WERE NOT controlled! '
             ib = elem%ibcur
             call elem%computeElementGeometry( gridN%x, gridN%b_edge(ib)%x_inn  )
          else
             call elem%computeElementGeometry( gridN%x )
          endif

!          call elem%computeElementGeometry( gridN%x, gridN%b_edge(ib)%x_inn  )

          !print *,'@@ 4'
          call SetElementQuadraturesDegrees(elem )
          !print *,'@@ 5'
          call PrepareOneElement(elem )
          !print *,'@@ 6'
          call ComputeLocalMassMatrix(elem )
          !print *,'@@ 7'
          call InitMatrixElementShape(elem )
          !print *,'@@ 8'

          ! initialization of the state vector from the previous grid
          ! for AMA and ANI and Ahp alreeady done in a different way
          if(state%space%adapt%adapt_type /= 'AMA' .and. state%space%adapt%adapt_type /= 'ANI' &
               .and. state%space%adapt%adapt_type /= 'Ahp' .and. state%space%adapt%adapt_type /= 'Ihp') then

             if(state%time%disc_time == 'STDG' ) then

                allocate(elem%w(0:1, 1:ndim*elem%dof) )  ! for  ComputeTerms( )
                allocate(elem%wST(1:ndim, 1:elem%dof, 1:elem%Tdof ) )
                !NEW for adaptation
                allocate(elem%rhsST(1:ndim, 1:elem%dof_plus, 1:elem%Tdof_plus ))
                allocate( elem%wSTfinAD(1:ndim, 1:elem%dof) ) !
             else

                allocate(elem%w(0:state%time%deg+1, 1:ndim* elem%dof) )

             endif

             ! solution recomputed in AdvancedInterpolDGsolution
             elem%to_recompute = .true.

             ! OLD variant see "WRTQ"
             ! if(elem%RGtype  == 'N' .and. elem%psplit == 0 ) then
             !    elem%w(0:state%time%deg+1, 1:ndim* elem%dof) = &
             !         gridO%elem(elem%i)%w(0:state%time%deg+1, 1:ndim* elem%dof)
             !    !print*,'ReprepareProblem :',i, elem%i, elem%RGtype

             ! elseif(elem%RGtype  == 'R' .or. elem%RGtype  == 'G' .or. &
             !      (elem%RGtype  == 'N' .and. elem%psplit /= 0) )then

             !    call RecomputeElemSol(elem, gridO%elem(elem%i))
             !    !print *,'@@ 9'
             ! elseif(elem%RGtype  == 'D' ) then
             !    !print *,'@@ 10'
             !    call RecomputeElemSolDerefine(elem, gridO)
             ! else
             !    !write(*,'(a6,i5,12es12.4)') '@@ 11',elem%i, elem%w(0,1:3)
             ! endif

          endif ! (state%space%adapt%adapt_type /= 'AMA' .and. state%space%adapt%adapt_type /= 'ANI' and 'Ahp')

       endif  ! (SimplyPass)

       elem%i = i
    enddo

    ! NEW variant of the recomputation "WRTQ"
    call AdvancedInterpolDGsolution(gridN, gridO )

    state%space%max_Qdof = maxval(gridN%elem(:)%Qdof) ! usually fGdof <= Qdof, hence this is OK
    !!!state%dof = ndim * sum(gridN%elem(:)%dof)

    state%space%h = maxval(gridN%elem(:)%diam)
    grid%h = state%space%h

    write(*, '(a38, i7, a2, i8, a2, i10, a3, f6.3,a2)') &
         '#Previous size (nelem, dof, nonzero):', gridO%nelem, ', ', &
         state%nsize,', ', state%nonzero,' (=',&
         1.*state%nonzero/state%nsize /state%nsize*100.,'%)'

    ! assempling matrix block together
    call InitGlobalMatrixShape(gridN)
    !print *,'@@ 10'
    ! setting of boundary conditions
    call SetConstBC(gridN)
    !print *,'@@ 11'

    if(state%space%adapt%adapt_type == 'AMA' .or. state%space%adapt%adapt_type == 'ANI' &
         .or. state%space%adapt%adapt_type == 'Ahp' .or. state%space%adapt%adapt_type == 'Ihp') return

    ! removing old elements which are not further used
    do i=1,gridO%nelem
       !print *,'@@ 11a',i
       if(gridO%elem(i)%dealloc ) call DeallocateElem(gridO%elem(i))
    enddo

    ! removing old bound_edge
    do i=1,gridO%nbelm
       !print *,'@@ 11b', i

       if(gridO%b_edge(i)%icurv >0) deallocate(gridO%b_edge(i)%x_inn)
       if(gridO%b_edge(i)%BC == 0) deallocate(gridO%b_edge(i)%x_div)
    enddo
    ! removing arrays of old grid
    deallocate(gridO%x)
    deallocate(gridO%xcur)
    deallocate(gridO%elem)
    deallocate(gridO%b_edge)
    !print *,'@@ 12'

    if (state%space%adapt%adapt_method == 'ALG2') then  !state%time%iter_SC has to be initialized in compute. f90 in order
                                            !not to loose the value from previous adapt cycle
       state%estim(:, :)  = 0.

       state%stop_rem = .true.
       state%stop_alg = .true.
       state%stop_lin = .true.

       state%stop_rem_glob = .true.
       state%stop_alg_glob = .true.
       state%stop_lin_glob = .true.
    endif

    if (state%space%adapt%adapt_method == 'ALG') then
       state%estim(:, :)  = 0.

       state%stop_rem = .true.
       state%stop_alg = .true.
       state%stop_lin = .true.

       state%stop_rem_glob = .true.
       state%stop_alg_glob = .true.
       state%stop_lin_glob = .true.

       ! reinitialization of the time step when Stopping Criteria based on AEE are satisfied
       state%time%iter_SC = -1
    endif

    ! multigrid
    if( state%MGsolver ) call InitMG( )



  end subroutine ReprepareProblem


  subroutine DeallocateGrid(gridA)
    type(mesh), intent(inout) :: gridA
    integer :: i

    print*,'## subroutine DeallocateGrid(gridA)'


    do i=1,grid%nelem
       call DeallocateElem(gridA%elem(i))
    enddo

    do i=1,gridA%nbelm
       if(gridA%b_edge(i)%icurv >0) deallocate(grid%b_edge(i)%x_inn)
       if(gridA%b_edge(i)%BC == 0) deallocate(grid%b_edge(i)%x_div)
    enddo


    deallocate(gridA%x)
    deallocate(gridA%xcur)
    deallocate(gridA%elem)
    deallocate(gridA%b_edge)
  end subroutine DeallocateGrid

  !> recomputation of the solution from several elements of "type(mesh) :: grid"
  !> (old grid) to new larger
  !> "type(element) :: elem", "elem" is a derefinement of elems from grid
  subroutine RecomputeElemSolDerefine(elem, grid)
    type(element), intent(inout) :: elem
    type(mesh), intent(in) :: grid
    class(element), pointer :: Delem
    real, dimension (:,:), allocatable :: x, xi! barycentric coordinates
    real, dimension (:,:), allocatable :: Ephi  ! test functions of large elem
    real, dimension (:,:), pointer ::    phi    ! test functions of small elements
    real, dimension (:,:,:), allocatable :: wi   ! integral of solution over Ks
    real, dimension (:,:,:), allocatable :: void
    real, dimension (:), allocatable :: wloc   ! solution in integ. nodes
    real, dimension (:), allocatable :: weights ! weights
    real, dimension (:),   allocatable :: vec  ! local array
    real, dimension (1:3, 1:3) :: loc_lam      ! local barycentric coordinates
    integer :: dof, Qdof,Qnum, Edof, i,j,k, l, kst
    integer :: RGhis, index, subel


    RGhis = elem%RGhistory
    subel = elem%RGreclen

    !print*,'@@@',elem%i, elem%RGrecomput(1:subel)

    Qdof = maxval(grid%elem( elem%RGrecomput(1:subel))%Qdof )
    allocate(x(1:Qdof, 1:3) )
    allocate(xi(1:Qdof, 1:3) )

    allocate(weights(1:Qdof) )

    Edof = elem%dof
    allocate(Ephi(1:Edof, 1:Qdof) )

    allocate(wi(0:state%time%deg+1, 1:Edof, 1:ndim) )
    wi(:,:,:) = 0.

    allocate(wloc(1:Qdof) )


    do i=1,subel
       Delem => grid%elem( elem%RGrecomput(i))

       dof = Delem%dof
       Qdof = Delem%Qdof
       Qnum = Delem%Qnum

       phi => state%space%V_rule(Qnum)%phi(1:dof,1:Qdof)

       x(1:Qdof, 1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)
       x(1:Qdof, 3) = 1. - x(1:Qdof, 1) - x(1:Qdof, 2)

       loc_lam(1:3, 1) = state%space%adapt%RGred(i, 2, 1:3)
       loc_lam(1:3, 2) = state%space%adapt%RGred(i, 3, 1:3)
       loc_lam(1:3, 3) = state%space%adapt%RGred(i, 1, 1:3)


       ! integration nodes of element in the barycentric coordinates with respect
       ! to the old element
       do l=1,Qdof
          xi(l, 1:3) = matmul( loc_lam(1:3, 1:3), x(l, 1:3) )
          !if(elem%i == 1) write(*,'(4es12.4, a4)' ) x(l, 1:nbDim),xi(l, 1:nbDim),'  ,,'

          !write(*,'(3es12.4,6i5)') xi(l, 1:3), l, i

       enddo

       allocate(void(1:Edof, 1:nbDim, 1:Qdof) )

       ! evaluation of basis function in "transformated" integ. node
       call PHI_orthonormal(Qdof, nbDim, xi(1:Qdof,1:nbDim), elem%type, Edof, Ephi(1:Edof, 1:Qdof), &
           void(1:Edof, 1:nbDim, 1:Qdof) )

       deallocate(void)


       ! integration  of the solution in "transformated" integ. node
       do j= 0, state%time%deg+1
          do k=1, ndim
             kst = dof*(k-1) + 1

             ! evaluation of the solution on old(small) element in integ nodes
             do l=1,Qdof
                wloc(l) = dot_product(Delem%w(j,kst:kst+dof-1), phi(1:dof,l) )
             enddo


             ! integration of w * Ephi over old (small element)
             call Eval_V_Weights(Delem, weights(1:Qdof) )

             do l=1,Edof

                wi(j,l,k) = wi(j,l,k) + &
                     dot_product(weights(1:Qdof)*wloc(1:Qdof), Ephi(l, 1:Qdof) )

             enddo
          enddo
       enddo
    enddo


    do j= 0, state%time%deg+1
       do k=1, ndim
          kst = Edof*(k-1) + 1

          elem%w(j,kst:kst+Edof-1) = matmul(elem%MassInv%Mb, wi(j,1:Edof,k) )
       enddo

    enddo

!    print*,'Stopped (23) in adaptation.f90'
!    stop


    deallocate(x, xi, Ephi, wi, weights, wloc)

  end subroutine RecomputeElemSolDerefine

  !> recomputation of the solution from "type(element) :: old_elem" to
  !> "type(element) :: new_elem", "new_elem" is a refinement of old_elem
  !> having the same polynomial degree
  subroutine RecomputeElemSol(new_elem, old_elem)
    type(element), intent(inout) :: old_elem
    type(element), intent(inout):: new_elem
    real, dimension (:,:), allocatable :: x, xi! barycentric coordinates
    real, dimension (:,:), allocatable :: phi  ! values of test functions in integ. nodes
    real, dimension (:,:), allocatable :: wi   ! solution in integ. nodes
    real, dimension (:,:,:), allocatable :: void
    real, dimension (:),   allocatable :: vec  ! local array
    real, dimension (1:3, 1:3) :: loc_lam      ! local barycentric coordinates
    integer :: dof, Qdof,Qnum, old_dof, i,j,k, l, kst
    integer :: index

    !if(new_elem%deg /= old_elem%deg) then
    !   print*,'Attention, different degrees in "RecomputeSolution: in adaptation.f90'
    !endif

    old_dof = old_elem%dof

    dof = new_elem%dof
    Qnum = new_elem%Qnum
    Qdof = new_elem%Qdof

    allocate(x(1:Qdof, 1:3) )
    allocate(xi(1:Qdof, 1:3) )
    allocate(phi(1:old_dof, 1:Qdof) )
    allocate(wi(1:Qdof, 1:ndim) )
    allocate(vec(1:dof*ndim) )

    x(1:Qdof, 1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)
    x(1:Qdof, 3) = 1. - x(1:Qdof, 1) - x(1:Qdof, 2)


    !write(*,'(a6, 3i5,2es14.6,a3,2i5)') '!!!!!',new_elem%i,&
    !     new_elem%Qnum, new_elem%Qdof, &
    !     new_elem%xc(:),new_elem%RGtype, new_elem%RGhistory, new_elem%RGindex

    !print*,'Old:', old_elem%i, old_elem%xc(:)
    !print*,'New:', new_elem%i, new_elem%xc(:)


    if(new_elem%RGtype == 'N' ) then  ! no h-refinement, only p-refinement
       loc_lam(1:3, 1) = state%space%adapt%RGred(0, 2, 1:3)
       loc_lam(1:3, 2) = state%space%adapt%RGred(0, 3, 1:3)
       loc_lam(1:3, 3) = state%space%adapt%RGred(0, 1, 1:3)

    elseif(new_elem%RGtype == 'R' ) then  ! red refinement
       !!loc_lam(1:3, 1:3) = state%RG%red(new_elem%RGlevel)%lambda(new_elem%RGindex, 1:3, 1:3)
       if(new_elem%RGhistory ==0) then
          index = new_elem%RGindex
       else
          index = 4
       endif
       !print*,'RG index',new_elem%i, index
       loc_lam(1:3, 1) = state%space%adapt%RGred(index, 2, 1:3)
       loc_lam(1:3, 2) = state%space%adapt%RGred(index, 3, 1:3)
       loc_lam(1:3, 3) = state%space%adapt%RGred(index, 1, 1:3)

    elseif(new_elem%RGtype == 'G') then  ! green refinement

       print*,'Not yet reimplemented, RGhistory has different meaning'
       stop
       !loc_lam(1:3, 1) = state%space%adapt%RGgreen(RGlevel, new_elem%RGindex, 2, 1:3)
       !loc_lam(1:3, 2) = state%space%adapt%RGgreen(RGlevel, new_elem%RGindex, 3, 1:3)
       !loc_lam(1:3, 3) = state%space%adapt%RGgreen(RGlevel, new_elem%RGindex, 1, 1:3)
    else
       print*,'Mismatch in RecomputeElemSol'
       stop
    endif

    ! integration nodes of new_element in the barycentric coordinates with respect
    ! to the old element
    do l=1,Qdof
       xi(l, 1:3) = matmul( loc_lam(1:3, 1:3), x(l, 1:3) )
    enddo


    !if(new_elem%RGtype == 'N' ) then
    !   print*,'------------------------'
    !   do l=1,Qdof
    !      write(*,'(6es12.4)') xi(l, 1:3) , x(l, 1:3)
    !   enddo
    !   !stop
    !endif

    allocate(void(1:old_dof, 1:nbDim, 1:Qdof) )
    ! evaluation of basis function in "transformated" integ. node
    call PHI_orthonormal(Qdof, nbDim, xi(1:Qdof, 1:nbDim), new_elem%type, old_dof, &
         phi(1:old_dof, 1:Qdof), void(1:old_dof, 1:nbDim, 1:Qdof) )
    deallocate(void)

    ! evaluation of the solution in "transformated" integ. node
    do j= 0, state%time%deg+1
       do k=1, ndim
          kst = old_dof*(k-1) + 1
          do i=1,Qdof
             wi(i,k) = dot_product(old_elem%w(j,kst:kst+old_dof-1), phi(1:old_dof,i) )
          enddo
       enddo

       vec(:) = 0.
       call EvalVectorB(new_elem, wi, dof, vec)

       do k=1, ndim
          kst = dof*(k-1) + 1
          new_elem%w(j,kst:kst+dof-1) = matmul(new_elem%MassInv%Mb, vec(kst:kst+dof-1) )
       enddo

    enddo
    deallocate(xi, phi, wi, vec)

  end subroutine RecomputeElemSol

  !> passing pointers and values from "type(element) :: old_elem" to
  !> "type(element) :: new_elem"
  subroutine PassElem2Elem(new_elem, old_elem )
    class(element), intent(inout) :: old_elem
    class(element), intent(inout):: new_elem
    logical :: rset
    integer :: j, nsl, nsl1

    allocate(new_elem%xc(1:nbDim) )
    new_elem%xc = old_elem%xc
    new_elem%flen = old_elem%flen
    new_elem%type = old_elem%type
    new_elem%area  = old_elem%area
    new_elem%diam  = old_elem%diam
    new_elem%r_ins  = old_elem%r_ins
    new_elem%limit_par  = old_elem%limit_par
    ! curved variable has to be reordered

    new_elem%F  => old_elem%F
    new_elem%n  => old_elem%n
    new_elem%dn  => old_elem%dn
    new_elem%nc  => old_elem%nc
    new_elem%dnc  => old_elem%dnc

    new_elem%Qnum = old_elem%Qnum
    new_elem%Qdeg = old_elem%Qdeg
    new_elem%Qdof = old_elem%Qdof
    new_elem%Tdeg = old_elem%Tdeg
    new_elem%Tdof = old_elem%Tdof
    new_elem%TQnum = old_elem%TQnum
    new_elem%Tdof_plus = old_elem%Tdof_plus

    new_elem%CP = old_elem%CP   ! Poincare constant

    new_elem%RTNflux => old_elem%RTNflux
    new_elem%res_func => old_elem%res_func   !residual function made up from residual vector


    allocate(new_elem%eta(1:max_eta, 1:ndim) )

    allocate(new_elem%estimFNCD(1:ndim) )
    allocate(new_elem%estimFNCA(1:ndim) )
    allocate(new_elem%estimFNCT(1:ndim) )

    ! the corresponding element was changed??
    rset = .false.
    if(new_elem%type /= old_elem%type .or. new_elem%flen /= old_elem%flen) then
       rset = .true.
    else
       do j=1,new_elem%flen
          if( new_elem%face(fdeg,j) /= old_elem%face(fdeg,j) ) rset = .true.
       enddo
    endif

    if(rset) then
       call SetElementQuadraturesDegrees(new_elem )
       print*,'Adaptation.f90, rset variable is .true. !!!'
    else
       new_elem%face(fGnum, 1:new_elem%flen) = old_elem%face(fGnum, 1:new_elem%flen)
       new_elem%face(fGdof, 1:new_elem%flen) = old_elem%face(fGdof, 1:new_elem%flen)
       new_elem%face(fTdeg, 1:new_elem%flen) = old_elem%face(fTdeg, 1:new_elem%flen)
       new_elem%face(fTdof, 1:new_elem%flen) = old_elem%face(fTdof, 1:new_elem%flen)
    endif

    call ComputeIntegNode(new_elem)

    ! for the same $p$ only !!!!!!!!!!!!!
    new_elem%Mass  => old_elem%Mass
    new_elem%MassInv  => old_elem%MassInv
    new_elem%Stiff  => old_elem%Stiff

    new_elem%w  => old_elem%w
    new_elem%wS  => old_elem%wS
    new_elem%MGw  => old_elem%MGw

    new_elem%wc  => old_elem%wc  ! used for storing computational solution based on AEE ALG

    ! STDG
    if(state%time%disc_time == 'STDG' ) then
       new_elem%wST  => old_elem%wST
       new_elem%rhsST  => old_elem%rhsST
    endif

    ! number of off-diagonal block is different
    if(new_elem%type /= old_elem%type .or. new_elem%flen /= old_elem%flen) then
       call DeleteMatrixElementShape(old_elem)
       call InitMatrixElementShape(new_elem)

       print*,'Adaptation.f90, rset variable is .true. !!! (2)'

    else
       ! resize of matrix block if dof of neighbouring element was changed
       new_elem%block  => old_elem%block
       new_elem%ILU  => old_elem%ILU
       new_elem%vec  => old_elem%vec

       new_elem%blockST  => old_elem%blockST
       new_elem%rhsST  => old_elem%rhsST
       new_elem%wSTfin   => old_elem%wSTfin
       new_elem%wSTfinAD => old_elem%wSTfinAD

       do j=1,new_elem%flen
          if( new_elem%face(neigh,j) > 0 .and. &
               new_elem%face(fdof,j) /= old_elem%face(fdof,j) ) then

             print*,'Adaptation.f90, rset variable is .true. !!! (3)'

             nsl = new_elem%dof * ndim
             nsl1= new_elem%face(fdof,j) * ndim

             call DeleteMblock(new_elem%block(j) )
             call InitMblock(new_elem%block(j), nsl, nsl1 )

             call DeleteMblock(new_elem%ILU(j) )
             call InitMblock(new_elem%ILU(j), nsl, nsl1 )
          endif
       enddo

    endif

  end subroutine PassElem2Elem

  !> passing pointers and values from "type(element) :: old_elem" to
  !> "type(element) :: new_elem"
  subroutine DeallocateElem(elem )
    type(element), intent(inout) :: elem
    integer :: j

    !print*,'&&&&', elem%i, elem%xc(:)
    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z0'
    deallocate(elem%F%F)

    if( elem%F%iFlin ) then
    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z1'
       deallocate(elem%F%D1F0)
    else
    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z2'
       deallocate(elem%F%V%JF, elem%F%V%D1F)
    endif

    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z3'
    if(elem%ibcur > 0) then
       do j=1,elem%flen
    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z4'
          deallocate(elem%F%E(j)%JF, elem%F%E(j)%D1F)
       enddo
    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z5'
       deallocate(elem%nc, elem%dnc )
    endif

    deallocate(elem%F)

    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z6'
    deallocate(elem%n)
    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z7'
    deallocate(elem%dn)

    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z8'
    deallocate(elem%Mass%Mb)
    deallocate(elem%Mass)
    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z9'
    deallocate(elem%MassInv%Mb)
    deallocate(elem%MassInv)
    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z10'
    deallocate(elem%Stiff%Mb)
    deallocate(elem%Stiff)

    deallocate(elem%w)

    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z11'
    if (state%time%disc_time /= 'STDG') then
       deallocate(elem%block(0)%Mb, elem%ILU(0)%Mb)
       do j=1,elem%flen
          !print*,'DeallocateElem, elem%i =', elem%i, 'type Z12'
          if(elem%face(neigh,j) >0) deallocate(elem%block(j)%Mb, elem%ILU(j)%Mb)
          !print*,'DeallocateElem, elem%i =', elem%i, 'type Z12a'
       enddo
       deallocate(elem%block, elem%ILU)

    else

       deallocate(elem%wST, elem%rhsST, elem%wSTfinAD ,  elem%wSTfin )

       deallocate(elem%block(0)%Mb, elem%blockST(0)%Mb, elem%ILU(0)%Mb)
       do j=1,elem%flen
          !print*,'DeallocateElem, elem%i =', elem%i, 'type Z12x'
          if(elem%face(neigh,j) >0)  &
               deallocate(elem%block(j)%Mb, elem%ILU(j)%Mb, elem%blockST(j)%Mb )
       enddo
       deallocate(elem%block, elem%ILU, elem%blockST)
     endif

    deallocate(elem%face)

    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z15'
    if(max_MG > 0) deallocate(elem%MGw)
    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z16'
    deallocate(elem%vec)
    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z17'
    deallocate(elem%xc, elem%xi )
    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z18'
    deallocate(elem%eta)
    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z19'


    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z20'
    if( size(elem%ibc) >0 ) deallocate(elem%ibc, elem%tbc)

    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z21'
    deallocate(elem%RTNflux)
    !print*,'DeallocateElem, elem%i =', elem%i, 'type Z22'

    if ( state%space%adapt%adapt_method == 'ALG' .or. state%space%adapt%adapt_method == 'ALG2') then

       deallocate(elem%res_func)
       deallocate(elem%wc)
    endif

    if(allocated(elem%estimFNCD) ) deallocate(elem%estimFNCD)
    if(allocated(elem%estimFNCA) ) deallocate(elem%estimFNCA)
    if(allocated(elem%estimFNCT) ) deallocate(elem%estimFNCT)
    !deallocate(elem%estimFNCA)
    !deallocate(elem%estimFNCT)




    ! AMA technique
    if(state%space%adapt%adapt_type == 'ANI' .or. state%space%adapt%adapt_type == 'Ahp'.or. state%space%adapt%adapt_type == 'Ihp') &
         deallocate ( elem%rgabc )

    ! pNeu
    if( state%space%estim_space == 'pNeu') deallocate(elem%lifting)

  end subroutine DeallocateElem


end module mesh_adaptation
