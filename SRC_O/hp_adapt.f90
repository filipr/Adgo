!> hp mesh adaptation, with hanging nodes
module hp_adaptation
  use data_mod, old_grid => grid  !mesh oper IS NOT USED for global grid (adaptation)
  use data_mod, new_grid => gridN !mesh oper IS NOT USED for global grid (adaptation)

  use geometry
  use mesh_oper
  use problem_oper
  use mesh_adaptation

  !> global variables for arrays of elem%face and elem%vec
  integer :: hg_face = 1     ! 1st index of elem%HGface(*, : ) = faces
  integer :: hg_level = 2    ! 2nd index of elem%HGface(*, : ) = HG level
  integer :: hg_index = 3   ! 3rd index of elem%HGface(*, : ) = local index in HGlevel
  integer :: max_hg_face =3  ! maximum of the previous indexes (only for allocation)


  type, public ::  hp_node
     logical :: remove                         ! remove this node
     logical :: HG                             ! HG node?
     integer :: i                              ! index of node
     real, dimension (:), allocatable   :: x              ! coordinates
     integer, dimension(1:2) :: xcur           ! curved info
     integer :: per                            ! per = 0 => no periodicity, > 0 periodic par
  end type hp_node

  type, public ::  hp_element
     logical :: remove                         ! remove thois element
     integer :: i                              ! index of element
     integer :: ie                              ! index of element
     integer :: type                           ! =3 for triang, 4 for quad
     integer :: flen                           ! number of nodes including hanging
     integer :: deg                   ! degree of polynomial approximation
     integer :: Tdeg                  ! degree of polynomial approximation
     integer, allocatable, dimension(:,:) :: face ! indexes of vertices, neighs,....
     integer, allocatable, dimension(:) :: pnew   ! indexes of new vertices
     logical :: HGnode                            !  HG structure is allocated and defined
     integer, allocatable, dimension(:) :: HGvertex ! for elements with HG nodes: verteces
     integer, allocatable, dimension(:,:) :: HGface ! for elements with HG nodes: faces
     integer :: hsplit                             ! type of h-ref of h-deref
     integer :: psplit                             ! type of p-ref
     integer :: RGhistory                           ! number of history levels of RG ref
     integer :: RGlevel                             ! level of refinement
     type(RGref), allocatable, dimension(:) :: RG     ! data for red green (de)refinement
     character :: RGtype                       ! =N (none) =R (red), =G (green) refinement
     integer :: RGindex                        ! index within level ref (see type(RGref)
     integer, dimension(1:4) :: RGrecomput
     integer :: RGreclen
     integer :: per_boun                       ! index to hp_edge if periodic edge
  end type hp_element

  type, public :: hp_edge
     logical :: remove        ! remove this element
     integer :: lbn(4)        ! indexes of end nodes
     integer :: ibc           ! index of boundary element
     integer :: itc           ! index of adjacent element
     integer :: jtc           ! inner index of adjacent element
     integer :: mid_node      ! new middle node if any
  end type hp_edge

  !> mesh   \f$ {\cal T}_h = \{ K \}_{K\in {\cal T}_h} \f$
  type, public :: hp_mesh
     integer :: nelem, npoin, nbelm, nbc, max_nelem, max_npoin, max_nbelm
     real :: xper(2,2)                       ! for periodicity
     integer :: iper(2,2)                    ! for periodicity
     integer :: periodic
     type(hp_node), pointer, dimension(:) :: node      ! nodes coordinates
     type(hp_element), pointer, dimension(:) :: elem     ! elements
     type(hp_edge), pointer, dimension(:) :: b_edge ! boundary edges
  end type hp_mesh



  public:: AdaptMesh_hp_HG
  public:: EstimateNewGrid
  public:: RedHGrefinement
  public:: RedHGderefinement
  public:: FindMiddleHGNode
  public:: PassGrid2hpGrid
  public:: PasshpGrid2Grid
  public:: DRefineGrid
  public:: RemoveHangingNodes
  public:: PlotMesh_hp
  public:: DeletehpGrid
contains

  !> hp mesh adaptation
  subroutine AdaptMesh_hp_HG( )
    type(hp_mesh) :: hp_grid
    character(len=15) :: gridfile

    call EstimateNewGrid(old_grid, hp_grid)

    call PassGrid2hpGrid(old_grid, hp_grid)

    !if(state%space%adapt%adapt_level == 1) call PlotMesh_hp(hp_grid, 0 )

    call DRefineGrid(hp_grid)

    allocate(new_grid)
    call PasshpGrid2Grid(hp_grid, new_grid)


    !!write(*,'(a23,a8,i7,a4,i7)') ' # HG adapt:, ', &
    !!     'npoin:',old_grid%npoin,'  =>', new_grid%npoin, &
    !write(*,'(a11,a8,i7,a4,i7,a13,i7,a1)')' # HP adapt ', &
    !     'nelem:',old_grid%nelem,'  =>', new_grid%nelem,'   (removed ', &
    !     count(hp_grid%elem(1:hp_grid%nelem)%remove),')'
    !!write(*,'(a23,a8,i7,a4,i7)')' #                     ', &
    !!     'nbelm:',old_grid%nbelm,'  =>', new_grid%nbelm
    !write(*,'(a11,a8,i7,a4,i7)')' #                     ', &
    !     ' DOF :',ndim*sum(old_grid%elem(:)%dof),'  =>', &
    !     ndim*sum(new_grid%elem(:)%dof)


    !print*,' DeletehpGrid'
    !call PlotMesh_hp(hp_grid, state%space%adapt%adapt_level)
    call DeletehpGrid(hp_grid)

    !print*,'WriteMesh'
    !gridfile = 'grid0'
    !call WriteMesh(gridfile,new_grid)

    !call PlotMesh(grid, 'meshA')

    !print*,'Seeking of the neighbours'
    call SeekNeighbours(new_grid)

    !print*,'Reseeking of the curved boundary'
    new_grid%curved_deg = old_grid%curved_deg
    call SeekCurvedBoundary(new_grid)

    !print*,'Passing of data from grid to gridN, deallocation of grid'
    call ReprepareProblem(new_grid, old_grid)

    old_grid => new_grid
    !print*,'AdaptMesh_hp_HG -- done'

  end subroutine AdaptMesh_hp_HG

  !> for an element without HG nodes (elem%HGnode = .false.) initiate the
  !> corresponding arrays
  subroutine InitiateHGElem(elem)
    type(hp_element), intent(inout) :: elem  ! element
    integer ::  max_flen, i

    max_flen = 4 * 2**(state%space%adapt%max_HGlevel+1)

    if(.not. elem%HGnode) then ! elem has not any hanging node yet

       elem%HGnode = .true.

       allocate(elem%HGvertex(1:elem%type+1) )

       do i=1,elem%type+1
          elem%HGvertex(i) = i
       enddo

       allocate(elem%HGface(1:max_hg_face, 1:max_flen) ) ! max number of HG nodes

       elem%HGface(:,:) = 0
       do i=1,elem%type
          elem%HGface(hg_face,  i) = i
          elem%HGface(hg_level, i) = 0
          elem%HGface(hg_index, i) = 1
       enddo
    else
       print*,'Element has already HG structure'
       stop

    endif
  end subroutine InitiateHGElem

  !> include the new HG "node" in the sequence of the element "type(hp_element):: elem"
  !> between ip1, ip2
  subroutine IncludeHGnode(elem, iface, ip1, ip2, node)
    type(hp_element), intent(inout) :: elem  ! element
    integer, intent(in) :: ip1, ip2 ! insert new node "node" between ip1, ip2
    integer, intent(in) :: node  ! index of the HG node
    integer, intent(out) :: iface  ! index of face where HG node has to be added
    integer :: i

    do i=1, elem%flen
       if(elem%face(idx, i) == ip1) then
          iface = i
          goto 100
       endif
    enddo
100 continue
    if(ip2 /= elem%face(idx, mod(iface, elem%flen)+1)) then
       print*,'Mischmatch in hp_adapt.f90'
       print*,ip1,ip2, iface
       print*,'elem =',elem%i, elem%flen
       write(*,'(a6,6i5)') 'idx:',elem%face(idx,1 :elem%flen)
       stop
    endif

    ! shifting the arrays
    !if(elem%i == 14) then
    !   write(*,*) 'iface = ',iface
    !   write(*,'(a6,8i4)') '***1*:',elem%face(idx,1:elem%flen)
    !endif

    if(iface < elem%flen) then
       ! 1=idx, 2=neigh !!!, 3=nei_i not used anymore
       elem%face(1:2,iface+2:elem%flen+1) = elem%face(1:2,iface+1:elem%flen)
       ! 1 = hg_face, 2= hg_level, 3=hg_index
       elem%HGface(1:3, iface+2:elem%flen+1) = elem%HGface(1:3, iface+1:elem%flen)
    endif

    elem%flen = elem%flen+1
    !if(elem%i == 14) write(*,'(a6,8i4)') '***2*:',elem%face(idx,1:elem%flen)

    ! data of the HG node/ segment
    elem%face(idx, iface+1) = node

    !if(elem%i == 14) write(*,'(a6,8i4)') '***3*:',elem%face(idx,1:elem%flen)

    do i=1,elem%type+1
       if(elem%HGvertex(i) > iface) elem%HGvertex(i)= elem%HGvertex(i) + 1
    enddo


    elem%HGface(hg_face, iface+1) = elem%HGface(hg_face, iface)

    elem%HGface(hg_level, iface) = elem%HGface(hg_level, iface)+1
    elem%HGface(hg_level, iface+1) = elem%HGface(hg_level, iface)


    elem%HGface(hg_index, iface) = 2*elem%HGface(hg_index, iface) - 1
    elem%HGface(hg_index, iface+1) = elem%HGface(hg_index, iface) + 1

  end subroutine IncludeHGnode

  !> adaptation of the "type(hp_mesh) :: hp_grid"
  subroutine DRefineGrid(hp_grid)
    type(hp_mesh), intent(inout) :: hp_grid
    type(hp_element), pointer :: elem, elem1
    type(hp_node), pointer :: node1, node2, Newnode, Newnode2
    integer :: i,j, k1, k2, kk, ip1, ip2, ipc, l
    integer :: cbp1, cbp2, ibp1, ibp2
    logical :: periodic

    if(grid%periodic > 0) then
       print*,'Periodic grid, HG refinement not implemented !!!!!'
       !stop
    endif


    ! seeking of node for refinement,
    ! if a node is missing then a new one is included as a hanging one
    do i=1,hp_grid%nelem
       elem => hp_grid%elem(i)

       if(elem%hsplit == 4) then
          allocate( elem%pnew(1:elem%type))
          elem%pnew(1:elem%type) = 0

          if(.not. elem%HGnode) call InitiateHGElem(elem)

          ! seeking the middle edge nodes for each edge
          do j=1,elem%type
             k1 = elem%HGvertex(j)
             k2 = elem%HGvertex(j+1)

             if( k2 == k1 + 1) then              ! new node has to be inserted
                hp_grid%npoin = hp_grid%npoin+1

                if(hp_grid%npoin > hp_grid%max_npoin ) then
                   print*,'Too much #npoin in hp_adapt.f90'
                endif

                Newnode => hp_grid%node(hp_grid%npoin)
                allocate(Newnode%x(1:nbDim) )

                Newnode%i = hp_grid%npoin
                Newnode%remove = .false.
                Newnode%per = 0

                ip1 = elem%face(idx,k1)
                ip2 = elem%face(idx,mod(k1,elem%flen)+1)

                node1 => hp_grid%node(ip1)
                node2 => hp_grid%node(ip2)


                ! coorddinates of new node
                if(node1%xcur(1)> 0 .and.  node2%xcur(1)> 0 ) then ! curved edge

                   Newnode%x(:) = (node1%x(:) +node2%x(:))/2.

                   call SeekCurvedNodeLocal(Newnode%x, node1%xcur, node2%xcur, &
                        2*Distance(node1%x(:), node2%x(:)), ipc)

                   Newnode%xcur(1) = node1%xcur(1)
                   Newnode%xcur(2) = ipc

                   Newnode%x(:) = curvbound%CBP(node1%xcur(1))%pnt(ipc,1:nbDim)


                else ! non-curved edge

                   Newnode%x(:) = (node1%x(:) +node2%x(:))/2.
                   Newnode%xcur(1:nbDim) = 0
                endif

                elem%pnew(j) = hp_grid%npoin

                call IncludeHGnode(elem, kk, ip1, ip2, hp_grid%npoin)
                Newnode%HG = .true.

                !if(node1%per > 0 .and. node2%per > 0)  then
                !   print*,'#######################################'
                !   write(*,'(a6,14i5)') 'pHG.', hp_grid%npoin -1,  hp_grid%npoin,ip1, ip2
                !   write(*,'(a6,14i5)') 'pHG..',elem%i, elem%flen,  elem%face(idx, 1:elem%flen)
                !   write(*,'(a6,12es12.4)') 'pHG...', &
                !        hp_grid%node(elem%face(idx,k1))%x(:), hp_grid%node(elem%face(idx,k2))%x(:)
                !endif

                !periodic ??
                periodic = .false.
                if(node1%per > 0 .and. node2%per > 0) then
                   periodic = .true.
                   Newnode2 => hp_grid%node(hp_grid%npoin+1)
                   allocate(Newnode2%x(1:nbDim) )

                   Newnode2%i = hp_grid%npoin+1
                   Newnode2%remove = .false.

                   Newnode2%per = Newnode%i
                   Newnode%per = Newnode2%i
                   Newnode2%x(:)  = Newnode%x(:) &
                        + hp_grid%node(node1%per)%x(:)  - hp_grid%node(ip1)%x(:)

                   Newnode2%HG = .true.

                   !write(102,*) Newnode%x(:)
                   !write(102,*) Newnode2%x(:)
                   !write(102,*)
                   !!stop
                   !write(*,'(a6,9i5)') '???',elem%face(neigh,1:elem%flen)
                endif


                !if(elem%i == 4869) then
                !   write(*,'(a5,20i5)') 'elem:',elem%i,elem%flen,elem%RGlevel,&
                !        elem%RGhistory
                !   write(*,'(a5,20i5)') 'idx:',elem%face(idx,1:elem%flen)
                !   if(elem%HGnode) then
                !      write(*,'(a5,20i5)') 'HGver:',elem%HGvertex(1:elem%type+1)
                !      write(*,'(a5,20i5)') 'HGfac:',elem%HGface(hg_face,1:elem%flen)
                !      write(*,'(a5,20i5)') 'HGlev:',elem%HGface(hg_level,1:elem%flen)
                !      write(*,'(a5,20i5)') 'HGidx:',elem%HGface(hg_index,1:elem%flen)
                !   endif
                !   do l=1,elem%flen
                !      print*,hp_grid%node(elem%face(idx,l))%x(1:nbDim),elem%face(idx,l)
                !   enddo
                !   print*,'----------------------include HG node----------------'
                !endif




                if(elem%face(neigh,k1) > 0) then
                   ! include the new node in the neigh elem sequence
                   elem1 => hp_grid%elem(elem%face(neigh, k1))

                   if(.not. elem1%HGnode) call InitiateHGElem(elem1)

                   if(periodic ) then
                      !write(*,'(a6,3i5,20es12.4)') &
                      !     '!!!',kk,ip2, ip1, hp_grid%node(ip2)%x(:), hp_grid%node(ip1)%x(:)
                      !write(*,'(a6,3i5,20es12.4)') &
                      !     '!!!',kk,node2%per, node1%per, &
                      ! hp_grid%node(node2%per)%x(:), hp_grid%node(node1%per)%x(:)
                      !stop

                      ! also the boundary nodes
                      kk = elem%per_boun
                      hp_grid%b_edge(kk)%mid_node = hp_grid%npoin

                      !print*,'###  b_edge1',elem%i, kk, hp_grid%npoin

                      ! the second node newnode2
                      hp_grid%npoin = hp_grid%npoin + 1
                      kk = elem1%per_boun
                      hp_grid%b_edge(kk)%mid_node = hp_grid%npoin

                      !print*,'###  b_edge2',elem1%i, kk, hp_grid%npoin

                      call IncludeHGnode(elem1, kk, node2%per, node1%per, hp_grid%npoin)

                      !write(*,'(a6,14i5)') 'pHG1', node2%per, node1%per, &
                      !     hp_grid%npoin -1,  hp_grid%npoin
                      !write(*,'(a6,14i5)') 'pHG2',elem%i, elem1%i, elem1%flen,  &
                      !     elem1%face(idx,1:elem1%flen)
                      !write(*,'(a6,12es12.4)') 'pHG3', &
                      !     hp_grid%node(elem%face(idx,k1))%x(:), hp_grid%node(elem%face(idx,k2))%x(:)

                   else
                      call IncludeHGnode(elem1, kk, ip2, ip1, hp_grid%npoin)

                   endif



                   !if(elem1%i == 4869) then
                   !   write(*,'(a5,4i5,a2,i5)') 'elem:',elem1%i,elem1%flen,elem1%RGlevel,&
                   !        elem1%RGhistory,'||',elem%i
                   !   write(*,'(a5,20i5)') 'idx:',elem1%face(idx,1:elem1%flen)
                   !   if(elem1%HGnode) then
                   !      write(*,'(a5,20i5)') 'HGver:',elem1%HGvertex(1:elem1%type+1)
                   !      write(*,'(a5,20i5)') 'HGfac:',elem1%HGface(hg_face,1:elem1%flen)
                   !      write(*,'(a5,20i5)') 'HGlev:',elem1%HGface(hg_level,1:elem1%flen)
                   !      write(*,'(a5,20i5)') 'HGidx:',elem1%HGface(hg_index,1:elem1%flen)
                   !   endif
                   !   do l=1,elem1%flen
                   !      print*,hp_grid%node(elem1%face(idx,l))%x(1:nbDim),elem1%face(idx,l)
                   !   enddo
                   !   print*,'----------------------include HG node----------------'
                   !endif



                else ! boundary face, adding a new boundary face
                   kk = -elem%face(neigh,k1)
                   hp_grid%b_edge(kk)%mid_node = hp_grid%npoin
                   !write(*,*) 'Bound face:',kk, elem%i, hp_grid%b_edge(kk)%mid_node
                endif

                elem%face(neigh, k1+1) = elem%face(neigh, k1)

                if(elem%face(neigh, k1) >0) then
                   elem1%face(neigh,kk+1) = elem1%face(neigh,kk)
                endif


             elseif( k2 > k1 + 1) then ! hanging node is present, we use it

                elem%pnew(j) = FindMiddleHGNode(elem, j)

                if( .not. hp_grid%node(elem%pnew(j))%HG ) then
                   print*,'Troubles in DRefineGrid with HG nodes'
                   write(*,'(a6,4i5,16es12.4)') '@@@',elem%i, elem%pnew(j), k1, k2, &
                        hp_grid%node(elem%face(idx,k1))%x(:), hp_grid%node(elem%face(idx,k2))%x(:)
                   stop

                endif

                hp_grid%node(elem%pnew(j))%HG = .false.

             else
                print*,'Troubles in DRefineGrid in "hp_adapt.f90"'
                stop
             endif

          enddo
       endif
    enddo

    !print*,'@@@@ C3'

    ! mesh refinement
    do i=1,hp_grid%nelem
       elem => hp_grid%elem(i)
       if(elem%hsplit == 4) then
          !write(*,'(a6,6i4)') 'elem:', elem%i,  elem%pnew(1:elem%type)
          !print*,'???', elem%pnew(1:elem%type)
          call RedHGrefinement(hp_grid, elem)
       endif
    enddo

    ! mesh derefinement
    do i=1,hp_grid%nelem
       elem => hp_grid%elem(i)
       if(elem%hsplit == -4 .and. elem%RGhistory > 0 ) then
          !write(*,'(a6,6i4)') 'elem:', elem%i,  elem%pnew(1:elem%type)
          !print*,'???', elem%pnew(1:elem%type)
          call RedHGderefinement(hp_grid, elem)
       endif
    enddo

    call RemoveHangingNodes(hp_grid)

    !!if(state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) stop

  end subroutine DRefineGrid

  !> we remove the haning nodes from the elements lists
  subroutine RemoveHangingNodes(hp_grid)
    type(hp_mesh), intent(inout) :: hp_grid
    type(hp_element), pointer :: elem
    logical, allocatable, dimension(:) :: rem
    logical :: HGnoderemoved
    integer :: i,j, k, l

    allocate(rem(1:maxval(hp_grid%elem(:)%flen) ) )
    HGnoderemoved = .false.

    do i=1, hp_grid%nelem
       elem => hp_grid%elem(i)
       if(elem%HGnode .and. .not. elem%remove ) then
          rem(:) = .false.

          ! we seek the no-corner nodes which are not HG, they will be removed
          k = 1
          do j=1,elem%flen
             !write(*,'(a6,3i5,3i3, 2es12.4,10i3)') &
             !     'HG???',i,elem%i,elem%face(idx,j),elem%flen, j,k, &
             !     hp_grid%node(elem%face(idx,j))%x(1:nbDim), &
             !     elem%HGvertex(1:elem%type+1),elem%HGface(1:3,j)


             if(j == elem%HGvertex(k) ) then
                k = k + 1
             else
                if(.not. hp_grid%node(elem%face(idx,j))%HG) then
                   rem(j) = .true.
                   hp_grid%node(elem%face(idx,j))%remove = .true.
                   HGnoderemoved = .true.

                   !print*,hp_grid%node(elem%face(idx,j))%x(:),elem%i,&
                   !     hp_grid%node(elem%face(idx,j))%i, '>>>'

                endif
             endif
          enddo

          if(count(rem(1:elem%flen)) > 0) then
             ! shifting the arrays
             !write(*,'(a4,40i4)') '!!!!',elem%i,elem%flen,elem%face(idx,1:elem%flen),&
             !     elem%HGface(1,1:elem%flen),elem%HGface(2,1:elem%flen),&
             !     elem%HGface(3,1:elem%flen),elem%HGvertex(1:4)

             !print*,'....',elem%flen, rem(:)
             k = 1
             do
                !print*,'@@@@',k,elem%flen, rem(k)
                if(k > elem%flen) exit
                if(rem(k)) then
                   !if(elem%i == 1235) then
                   !   write(*,'(a8,20i5)') '   ...',k, elem%HGface(hg_level,k-1:k), &
                   !        elem%HGface(hg_index,k-1:k)
                   !endif

                   elem%HGface(hg_level,k-1) =  elem%HGface(hg_level,k-1) - 1
                   elem%HGface(hg_index,k-1) =  elem%HGface(hg_index,k) /2
                   !elem%HGface(hg_index,k-1) =  max(1, elem%HGface(hg_index,k-1) /2)

                   !if(elem%HGface(hg_index,k) /2 &
                   !     /= max(1, elem%HGface(hg_index,k-1) /2)) print*,'??????',elem%i

                   !if(elem%i == 1235) then
                   !   write(*,'(a8,20i5)') '   ...',k, elem%HGface(hg_level,k-1:k), &
                   !        elem%HGface(hg_index,k-1:k)
                   !endif

                   do j = k, elem%flen - 1
                      elem%face(idx,j) =  elem%face(idx,j+1)
                      elem%HGface(hg_face:hg_index, j) = elem%HGface(hg_face:hg_index, j+1)
                      rem(j) = rem(j+1)
                   enddo

                   do l=2,elem%type+1
                      if(elem%HGvertex(l) > k) elem%HGvertex(l) = elem%HGvertex(l)-1
                   enddo
                   elem%flen = elem%flen - 1

                   !write(*,'(a4,40i4)')'!!!!',elem%i,elem%flen,elem%face(idx,1:elem%flen),&
                   !     elem%HGface(1,1:elem%flen),elem%HGface(2,1:elem%flen), &
                   !     elem%HGface(3,1:elem%flen),   elem%HGvertex(1:4)

                endif
!!!print*,'....',k,elem%flen, rem(:)
                k = k+1
             enddo

             !write(*,'(a4,40i5)') 'elem',elem%i,elem%flen, size(elem%face(idx,:)), &
             !     size(elem%HGface(idx,:))
             !write(*,'(a4,40i5)') 'idx',elem%face(idx,1:elem%flen)
             !write(*,'(a4,40i5)') 'hged',elem%HGface(1,1:elem%flen)
             !write(*,'(a4,40i5)') 'hgle',elem%HGface(2,1:elem%flen)
             !write(*,'(a4,40i5)') 'hgi',elem%HGface(3,1:elem%flen)
             !write(*,'(a4,40i5)') 'hgver',elem%HGvertex(1:4)
             !print*,'-----------'

             ! element is still HG ?
             if(elem%flen == elem%type) then
                elem%HGnode = .false.
                deallocate(elem%HGface, elem%HGvertex)
             endif
          endif
       endif
    enddo

    ! formal shifting of nodes
    i = 0
    if(HGnoderemoved) then
       do j=1,hp_grid%npoin
          if(.not. hp_grid%node(j)%remove) then
             i = i+ 1
             hp_grid%node(j)%i = i
          endif
          !print*,'@@@',hp_grid%node(j)%remove,j,hp_grid%node(j)%i
       enddo
    endif
  end subroutine RemoveHangingNodes

  !> 'red' derifinement of the element with HG nodes
  subroutine RedHGderefinement(hp_grid, elem)
    type(hp_mesh), intent(inout) :: hp_grid
    type(hp_element), intent(inout) :: elem
    type(hp_element), pointer :: elem1, elem2
    integer :: j, j1, k1, k2, l, l1
    integer :: lev, subel, flen, nsize, len, ifirst
    integer :: iel, ied, iei, iej
    integer, dimension(1:6) :: el, ed, ei, ej

    el(1:6) = (/ 1, 2, 2, 3, 3, 1/)
    ed(1:6) = (/ 1, 3, 1, 3, 1, 3/)
    ei(1:6) = (/ 1, 2, 1, 2, 1, 2/)
    ej(1:6) = (/ 1, 1, 2, 2, 3, 3/)

    if(elem%type /= 3) then
       print*,' Only triangles in derefinement are impelmented'
       stop
    endif

    ! removing of "center" element will change the ""hanginity" of its nodes
    hp_grid%node(elem%face(idx,1:elem%flen))%HG &
         = .not. hp_grid%node(elem%face(idx,1:elem%flen))%HG

    lev = elem%RGhistory
    subel = elem%RG(lev)%subel

    elem%RGreclen = subel
    elem%RGrecomput(1:subel) = elem%RG(lev)%daughter(1:subel)


    elem%RGhistory = elem%RGhistory - 1
    elem%RGlevel = elem%RGlevel - 1
    elem%RGtype = 'D'


    if(subel == 4) nsize = 3  ! first three triangles will be removed
    flen = 0
    do j=1,nsize
       elem1 => hp_grid%elem(elem%RG(lev)%daughter(j))

       flen = flen + elem1%flen
       !write(*,'(a3,12i5)') '..,',elem1%i, elem1%face(idx,1:elem1%flen), &
       !      elem1%face(neigh,1:elem1%flen)

       elem1%remove = .true.
       elem1%RGlevel = 0
    enddo
    flen = flen - elem%type

    elem%flen = flen
    deallocate(elem%face)
    allocate(elem%face(idx:neigh,1:flen) )
    elem%face(:,:) = 0

    elem%HGnode = .true.
    allocate(elem%HGface(hg_face:hg_index,1:flen) )
    elem%HGface(:,:) = 0

    allocate(elem%HGvertex(1:elem%type+1) )
    elem%HGvertex(1:elem%type) = 0
    elem%HGvertex(elem%type+1) = flen + 1

    j = 1  ! node
    do l=1,6   ! new element has "6 edges"
       iel = el(l)   ! elem
       ied = ed(l)   ! edge
       iei = ei(l)   ! sub HG segment
       iej = ej(l)   ! sub HG vertex

       if(mod(l,2) == 1) elem%HGvertex(iej) = j

       elem1 => hp_grid%elem(elem%RG(lev)%daughter(iel))

       if(.not. elem1%HGnode) then
          elem%face(idx,j) = elem1%face(idx, ied)
          elem%HGface(hg_face,j) = iej
          elem%HGface(hg_level,j) = 1
          elem%HGface(hg_index,j) = iei

          !print*,'~~~~~~~~~~~~~~',iei, j, l, elem%HGface(hg_index,j)
          j = j + 1
       else
          ifirst = elem1%HGvertex(ied)
          len = elem1%HGvertex(ied+1) - ifirst -1

          elem%face(idx, j: j+ len) = elem1%face(idx, ifirst:ifirst+ len)

          elem%HGface(hg_face, j: j+ len) = iej
          elem%HGface(hg_level, j: j+ len) = elem1%HGface(hg_level, ifirst:ifirst+len) + 1
          elem%HGface(hg_index, j: j+ len) = elem1%HGface(hg_index, ifirst:ifirst+len) &
               + (iei-1)* 2**(elem%HGface(hg_level, j: j+ len)-1)

          !print*,'???',j,j+len,ifirst, elem1%HGvertex(ied+1)-1, ifirst+len
          !print*,'~~~',elem1%HGface(hg_level, ifirst:ifirst+len), &
          !     elem1%HGface(hg_index, ifirst:ifirst+len), &
          !     elem%HGface(hg_level, j: j+ len), elem%HGface(hg_index, j: j+ len)
          j = j + len + 1
       endif

    enddo

    !if(elem%i == 3863 .or. elem%i == 4869) then
    !   write(*,'(a6,15i4)') 'DEREf',elem%i,j,j1,hp_grid%nelem,elem1%face(1,1:elem1%flen)
    !   if(elem1%HGnode) then
    !      write(*,'(a6,15i4)') 'HG ver',elem1%HGvertex(1:elem1%type)
    !      write(*,'(a6,15i4)') 'hg_fac',elem1%HGface(hg_face, 1:elem1%flen)
    !      write(*,'(a6,15i4)') 'hg_lev',elem1%HGface(hg_level, 1:elem1%flen)
    !      write(*,'(a6,15i4)') 'hg_ind',elem1%HGface(hg_index, 1:elem1%flen)
    !   endif
    !   print*,'--------------------------------'
    !endif


    !write(*,'(a6,i4,a1,30i5)') '..face',elem%i,'|',elem%face(idx,1:elem%flen), &
    !     elem%HGface(hg_face,1:elem%flen)

    !write(*,'(a6,i4,a1,30i5)')'HGface',elem%i, '|',elem%HGface(hg_level,1:elem%flen),&
    !     elem%HGface(hg_index,1:elem%flen), elem%HGvertex(1:elem%type+1)
    !print*

  end subroutine RedHGderefinement

  !> 'red' splitting of the element with HG nodes
  subroutine RedHGrefinement(hp_grid, elem)
    type(hp_mesh), intent(inout) :: hp_grid
    type(hp_element), intent(inout) :: elem
    type(hp_element), pointer :: elem1, elem2
    integer :: j, j1, k1, k2, l, l1, max_flen

    max_flen = 4 * 2**(state%space%adapt%max_HGlevel+1)

    ! we change elem%pnew(1:elem%type) from the indexes of node to their orders
    !if(elem%i <= 6) print*,'OLD',elem%pnew(1:elem%type)

    j1 = 1
    !write(*,'(a6,20i4)') 'IDX',elem%i,elem%face(idx,1:elem%flen)
    !write(*,'(a6,20i4)') 'pnew',elem%pnew(:)
    do j=1,elem%flen
       !print*,'???',elem%i,j, j1, elem%face(idx,j),elem%pnew(j1)
       if(elem%face(idx,j) == elem%pnew(j1)) then
          elem%pnew(j1) = j
          j1 = j1+1
          if(j1 > elem%type) exit ! last pnew node found
       endif
    enddo
    if(j1 /= elem%type+1) then
       print*,'pnew index does not found',j1,elem%type
       stop
    endif

    !if(elem%i <= 6) print*,'NEW',elem%pnew(1:elem%type)

    ! let's go to refinement
    if(elem%type == 3) then !triangle
       elem%RGhistory = elem%RGhistory + 1
       elem%RGlevel = elem%RGlevel + 1

       do j=1,elem%type ! corner elements
          j1 = mod(j+1,elem%type) + 1
          ! j is the first index, j1 is the last(third) one

          hp_grid%nelem = hp_grid%nelem + 1

          if(hp_grid%nelem > hp_grid%max_nelem ) then
             print*,'Too much #nelem in hp_adapt.f90'
          endif

          elem1 =>  hp_grid%elem(hp_grid%nelem)

          elem1%remove = .false.
          elem1%i = elem%i
          elem1%type = elem%type
          elem1%deg = elem%deg
          elem1%Tdeg = elem%Tdeg
          elem1%hsplit = 0
          elem1%psplit = elem%psplit

          elem1%RGhistory = 0
          elem1%RGlevel = elem%RGlevel
          elem1%RGtype = 'R'
          elem1%RGindex = j

          elem%RG(elem%RGhistory)%daughter(j) = hp_grid%nelem

          elem1%flen = elem%pnew(j) - elem%HGvertex(j) &  !first edge
               + elem%HGvertex(j1+1) -elem%pnew(j1) + 1   !third and the second edge

          allocate(elem1%face(1:3,1:max_flen) )
          elem1%face(2:3,1:max_flen) = 0

          elem1%HGnode = .false.
          if(elem1%flen > elem1%type) then   ! HG elem
             elem1%HGnode = .true.
             allocate(elem1%HGvertex(1:elem%type+1) )
             allocate(elem1%HGface(1:3,1:max_flen) )
          endif

          ! first edge
          k1 = 1
          k2 = 1+elem%pnew(j) - elem%HGvertex(j)

          elem1%face(1,k1:k2) = elem%face(1, elem%HGvertex(j):elem%pnew(j) )

          if(elem1%HGnode) then
             elem1%HGvertex(1) = k1
             elem1%HGvertex(2) = k2
             elem1%HGvertex(3) = k2+1
             elem1%HGvertex(4) = elem1%flen + 1

             do l=k1,k2-1
                l1 = elem%HGvertex(j) + l - k1
                elem1%HGface(hg_face, l) = 1
                elem1%HGface(hg_level, l) = elem%HGface(hg_level, l1) - 1
                elem1%HGface(hg_index, l) = elem%HGface(hg_index, l1)

                !!print*,'~~~~',l,l1
             enddo
             elem1%HGface(hg_face, k2) = 2
             elem1%HGface(hg_level, k2) = 0
             elem1%HGface(hg_index, k2) = 1
          endif

          ! second edge, never HG node
          ! third edge
          k1 = k2+1
          k2 = elem1%flen  !k1 + elem%HGvertex(j1+1) -elem%pnew(j1)

          if(elem%i <= -6) print*,'?e2',k1,k2, elem%pnew(j1), elem%HGvertex(j1+1) - 1

          elem1%face(1,k1:k2) = elem%face(1, elem%pnew(j1):elem%HGvertex(j1+1)-1 )

          if(elem1%HGnode) then
             do l=k1,k2
                l1 = elem%pnew(j1) + l - k1
                elem1%HGface(hg_face, l) = 3
                elem1%HGface(hg_level, l) = elem%HGface(hg_level, l1) - 1
                elem1%HGface(hg_index, l) = elem%HGface(hg_index, l1) &
                     - 2**elem1%HGface(hg_level, l)
             enddo
          endif


          !if(elem%i == 3863 .or. elem%i == 4869) then
          !   write(*,'(a6,15i4)') 'new e',elem%i,j,j1,hp_grid%nelem,elem1%face(1,1:elem1%flen)
          !   if(elem1%HGnode) then
          !      write(*,'(a6,15i4)') 'HG ver',elem1%HGvertex(1:elem1%type)
          !      write(*,'(a6,15i4)') 'hg_fac',elem1%HGface(hg_face, 1:elem1%flen)
          !      write(*,'(a6,15i4)') 'hg_lev',elem1%HGface(hg_level, 1:elem1%flen)
          !      write(*,'(a6,15i4)') 'hg_ind',elem1%HGface(hg_index, 1:elem1%flen)
          !   endif
          !   do l=1,elem%flen
          !      print*,hp_grid%node(elem%face(idx,l))%x(1:nbDim),elem%face(idx,l)
          !   enddo
          !
          !   print*,'--------------------- HG refinement-----------'
          !endif
       enddo

       ! central element replaces to old one, no hanging node
       if(elem%HGnode) deallocate(elem%HGvertex, elem%HGface)
       elem%HGnode = .false.
       elem%remove = .false.
       elem%hsplit = 0

       !if(elem%i == 422) print*,'!!!!',elem%i,elem%HGnode

       ! RG structures
       elem%RG(elem%RGhistory)%daughter(4) =  elem%i
       elem%RGtype = 'R'

       do j=1,elem%type
          elem%face(idx,j) = elem%face(idx,elem%pnew(j))
       enddo
       elem%flen = elem%type


       !write(*,'(a6,15i4)') 'new e',elem%i,j,j1,hp_grid%nelem,elem%face(1,1:elem%flen)

       !if(elem%i == 6) stop

    else if(elem%type == 4) then  !quadrilaterall
       print*,'Not YET implemented'
       stop

    else
       print*,'Problem in RedHGrefinement in "hp_adapt.f90"'
       stop
    endif

  end subroutine RedHGrefinement


  !> finding the middle node of elem's iface-th face,
  !> where is a HG node
  function FindMiddleHGNode(elem, iface)
    integer :: FindMiddleHGNode
    type(hp_element), intent(in) :: elem
    integer , intent(in) :: iface
    integer :: i

    do i=elem%HGvertex(iface), elem%flen
       if(elem%HGface(hg_face,i) == iface) then
          if(elem%HGface(hg_index,i) == 2**(elem%HGface(hg_level,i)-1)) then
             FindMiddleHGNode = elem%face(idx,i+1)
!!!FindMiddleHGNode = i+1
             return
          endif
       endif
    enddo
    print*,"Middle HG node doesn't found"
    stop

  end function FindMiddleHGNode

  !> passing the "type(mesh):: grid" to "type(hp_mesh) :: hp_grid"
  subroutine PassGrid2hpGrid(grid, hp_grid)
    type(mesh), intent(in) :: grid
    type(hp_mesh), intent(inout) :: hp_grid
    class(element), pointer :: elem
    type(hp_element), pointer :: hp_elem
    integer :: i, ii, j, ip, k,k1, k2, max_flen
    integer :: l, is1, is2, ip1, ip2, ie, iie

    max_flen = 4 * 2**(state%space%adapt%max_HGlevel+1)
    !print*,'Max_flen =',max_flen
    hp_grid%xper(:,:) = grid%xper(:,:)
    hp_grid%iper(:,:) = grid%iper(:,:)
    hp_grid%periodic = grid%periodic

    ! nodes
    hp_grid%npoin = grid%npoin
    do i=1,hp_grid%npoin
       allocate(hp_grid%node(i)%x(1:nbDim) )
       hp_grid%node(i)%x(:) = grid%x(i,:)
       hp_grid%node(i)%xcur(:) = grid%xcur(i,:)
       hp_grid%node(i)%remove = .false.
       hp_grid%node(i)%HG = .false.
       hp_grid%node(i)%i = i
       hp_grid%node(i)%per = 0
       !!if( hp_grid%node(i)%xcur(1) > 0) print*,'((((((((',i,  hp_grid%node(i)%xcur(1:nbDim)
    enddo


    ! elements
    hp_grid%nelem = grid%nelem
    do i=1,grid%nelem
       elem => grid%elem(i)
       hp_elem => hp_grid%elem(i)

       hp_elem%remove = .false.
       hp_elem%i = i
       hp_elem%type = elem%type
       hp_elem%flen = elem%flen
       hp_elem%deg = elem%deg
       hp_elem%Tdeg = elem%Tdeg
       hp_elem%hsplit = elem%hsplit
       hp_elem%psplit = elem%psplit

       hp_elem%RGlevel = elem%RGlevel
       hp_elem%RGhistory = elem%RGhistory
       hp_elem%RGindex = elem%RGindex
       hp_elem%RGtype = 'N'
       hp_elem%RGreclen = 0
       hp_elem%per_boun = 0

       ip = hp_elem%RGhistory
       if(hp_elem%hsplit >0) ip = ip + 1

       if(ip >0) then
          allocate(hp_elem%RG(1: ip))

          do j=1,hp_elem%RGhistory
             hp_elem%RG(j)%subel = elem%RG(j)%subel

             allocate(hp_elem%RG(j)%daughter(1:hp_elem%RG(j)%subel))
             hp_elem%RG(j)%daughter(:) = elem%RG(j)%daughter(:)
          enddo

          if(ip > hp_elem%RGhistory) then
             j = hp_elem%RGhistory + 1
             hp_elem%RG(j)%subel = hp_elem%hsplit !! four daughter elements, red refinement

             allocate(hp_elem%RG(j)%daughter(1:hp_elem%RG(j)%subel))
             hp_elem%RG(j)%daughter(:) = 0
          endif

       endif

       allocate(hp_elem%face(1:2,1:max_flen) ) !!! max value of additional HG nodes
       hp_elem%face(idx,1:hp_elem%flen) = elem%face(idx,1:hp_elem%flen)
       hp_elem%face(neigh,1:hp_elem%flen) = elem%face(neigh,1:hp_elem%flen)
       !hp_elem%face(nei_i,1:hp_elem%flen) = elem%face(nei_i,1:hp_elem%flen)

       do j=1,hp_elem%flen
          if(hp_elem%face(neigh,j) < 0) hp_grid%node(hp_elem%face(idx,j))%HG = .true.
       enddo

       ! element with HanGing node(s) ?
       hp_elem%HGnode = elem%HGnode
       if(hp_elem%HGnode) then
          allocate(hp_elem%HGvertex(1:hp_elem%type+1) )
          hp_elem%HGvertex(1:hp_elem%type) = elem%HGvertex(1:hp_elem%type)
          hp_elem%HGvertex(hp_elem%type+1) = hp_elem%flen+1

          allocate(hp_elem%HGface(1:3, 1:max_flen) )  ! max number of HG nodes
          ! (1, :) = index of the edge 1..type
          ! (2, :) = level of HG node refinement
          ! (3, :) = local index of HG node refinement
          hp_elem%HGface(hg_face, 1:hp_elem%flen) = elem%HGface(hg_face, 1:hp_elem%flen)
          k = 1
          do j=1,hp_elem%flen

             ip = log(1.0001*elem%HGface(2, j) ) / log(2.)

             hp_elem%HGface(hg_level, j) = ip
             hp_elem%HGface(hg_index, j) = elem%HGface(hg_level, j) - 2**ip + 1

             ! we seek which node is hanging one
             if(j == hp_elem%HGvertex(k)) then
                k = k+ 1
             else
                if(hp_grid%node(hp_elem%face(idx,j))%HG ) then
                   print*,'Troubles in initiation of HG nodes'
                   stop
                endif
                hp_grid%node(hp_elem%face(idx,j))%HG = .true.
             endif

          enddo
       endif

       !if(i == 3863 .or. i == 4869) then
       !   write(*,'(a5,20i5)') 'elem:',hp_elem%i,hp_elem%flen,hp_elem%RGlevel,&
       !        hp_elem%RGhistory
       !   write(*,'(a5,20i5)') 'idx:',hp_elem%face(idx,1:hp_elem%flen)
       !   if(elem%HGnode) then
       !      write(*,'(a5,20i5)') 'HGver:',hp_elem%HGvertex(1:hp_elem%type+1)
       !      write(*,'(a5,20i5)') 'HGfac:',hp_elem%HGface(hg_face,1:hp_elem%flen)
       !      write(*,'(a5,20i5)') 'HGlev:',hp_elem%HGface(hg_level,1:hp_elem%flen)
       !      write(*,'(a5,20i5)') 'HGidx:',hp_elem%HGface(hg_index,1:hp_elem%flen)
       !   endif
       !   print*,'--------------------------------------'
       !endif
    enddo

    ! boundary faces
    hp_grid%nbelm = grid%nbelm
    do i=1,grid%nbelm
       hp_grid%b_edge(i)%remove = .false.
       hp_grid%b_edge(i)%lbn(:) = grid%b_edge(i)%lbn(:)
       hp_grid%b_edge(i)%ibc = grid%b_edge(i)%ibc
       hp_grid%b_edge(i)%itc = grid%b_edge(i)%itc
       hp_grid%b_edge(i)%jtc = grid%b_edge(i)%jtc
       hp_grid%b_edge(i)%mid_node = 0
    enddo



    ! seeking of periodicity pairs
    do l=1, hp_grid%periodic
       is1 = 0
       is2 = 0
       ip1 = 0
       ip2 = 0
       do i=1,grid%nbelm
          if(is1 == 0 .and. grid%b_edge(i)%ibc == grid%iper(l,1) ) is1 = i
          if(grid%b_edge(i)%ibc == grid%iper(l,1) ) ip1 = ip1 + 1

          if(is2 == 0 .and. grid%b_edge(i)%ibc == grid%iper(l,2) ) is2 = i
          if(grid%b_edge(i)%ibc == grid%iper(l,2) ) ip2 = ip2 + 1

       enddo
       write(*,'(a4,7i5)') 'per', l,grid%iper(l,1:2), is1, is2, ip1, ip2

       do k = 0, ip1 - 1
          i = is1 + k
          ii = is2 + ip2 - k - 1

          !print*,'@@@', k,i, ii


          ie  = grid%b_edge( i)%itc
          iie = grid%b_edge(ii)%itc

          hp_grid%elem(ie)%per_boun = i
          hp_grid%elem(iie)%per_boun = ii

          if(k ==0 ) then
             !first node only for the first segment
             k1  = grid%b_edge( i)%lbn(1)
             k2  = grid%b_edge(ii)%lbn(2)
             hp_grid%node(k1)%per = k2
             hp_grid%node(k2)%per = k1
          endif
          ! second node
          k1  = grid%b_edge( i)%lbn(2)
          k2  = grid%b_edge(ii)%lbn(1)
          hp_grid%node(k1)%per = k2
          hp_grid%node(k2)%per = k1
       enddo

    end do

    !do i=1,hp_grid%nelem
    !   if(hp_grid%elem(i)%per_boun > 0) print*,'#####',i,hp_grid%elem(i)%per_boun
    !enddo

    do i=1,hp_grid%npoin
       if( hp_grid%node(i)%per > 0) then
          if( i > hp_grid%node(i)%per ) then
             write(101,*) hp_grid%node(i)%x(:)
             write(101,*) hp_grid%node(hp_grid%node(i)%per)%x(:)
             write(101,*)
          endif
       endif
    enddo

    !print*,' End of periodicity'
    !stop

  end subroutine PassGrid2hpGrid


  !> passing the "type(hp_mesh) :: hp_grid" to "type(mesh):: grid"
  subroutine PasshpGrid2Grid(hp_grid, grid)
    type(hp_mesh), intent(in) :: hp_grid
    type(mesh), intent(inout) :: grid
    class(element), pointer :: elem
    type(hp_element), pointer :: hp_elem
    integer :: i,j,ie, ip, l

    grid%xper(:,:) = hp_grid%xper(:,:)
    grid%iper(:,:) = hp_grid%iper(:,:)
    grid%periodic = hp_grid%periodic


    ! nodes
    !grid%npoin = hp_grid%npoin
    grid%npoin = count(.not. hp_grid%node(1:hp_grid%npoin)%remove)

    !print*,'~~~~~',hp_grid%npoin, grid%npoin

    allocate(grid%x(1:grid%npoin, 1:nbDim))
    allocate(grid%xcur(1:grid%npoin, 1:nbDim))

    ip = 0
    do i=1,hp_grid%npoin
       if(.not. hp_grid%node(i)%remove) then
          ip = ip+1
          grid%x(ip, :) = hp_grid%node(i)%x(:)
          grid%xcur(ip, :) = hp_grid%node(i)%xcur(:)
       endif
    enddo


    ! elements
    ! reindexing
    ie = 0
    do i=1,hp_grid%nelem
       hp_elem => hp_grid%elem(i)
       if(.not. hp_elem%remove) then
          ie = ie + 1
          hp_elem%ie = ie
       endif
    enddo

    ! creating new list
    grid%nelem = count(.not. hp_grid%elem(1:hp_grid%nelem)%remove)
    allocate(grid%elem(1:grid%nelem))

    !print*,'???',grid%nelem, hp_grid%nelem, ie
    ie = 0
    do i=1,hp_grid%nelem
       hp_elem => hp_grid%elem(i)

       if(.not. hp_elem%remove) then
          ie = ie + 1
          elem => grid%elem(ie)

          elem%i = hp_elem%i
          elem%type = hp_elem%type
          elem%flen = hp_elem%flen
          elem%deg = hp_elem%deg
          elem%Tdeg = hp_elem%Tdeg
          elem%hsplit = hp_elem%hsplit
          elem%psplit = hp_elem%psplit

          ! only for counting the total DOF
          elem%dof = (elem%deg + elem%psplit + 1)*(elem%deg+ elem%psplit+2)/2

          elem%RGlevel = hp_elem%RGlevel
          elem%RGhistory = hp_elem%RGhistory
          elem%RGtype = hp_elem%RGtype
          elem%RGindex = hp_elem%RGindex


          !if((ie >= 2687 .and. ie <=2689) .or. ie == 2009 .or. ie == 1053 .or. &
          !     ie == 446) then
          !   write(*,'(a6,4i5, 2e12.4)') &
          !        '@@@@',ie,elem%RGlevel, elem%RGhistory, elem%RGindex, &
          !        hp_grid%node(hp_elem%face(1,1))%x(1:nbDim)
          !endif

          !print*,'@@@',elem%RGhistory,elem%RGtype
          if(elem%RGhistory >0) then
             allocate(elem%RG(1: elem%RGhistory))

             do j=1,hp_elem%RGhistory
                elem%RG(j)%subel = hp_elem%RG(j)%subel

                allocate(elem%RG(j)%daughter(1:elem%RG(j)%subel))

                elem%RG(j)%daughter(:) = hp_grid%elem(hp_elem%RG(j)%daughter(:))%ie
             enddo
          endif

          if(elem%RGtype == 'D') then
             elem%RGreclen = hp_elem%RGreclen
             elem%RGrecomput(1:elem%RGreclen) = hp_elem%RGrecomput(1:elem%RGreclen)
             elem%RGrecomput(1:4) = hp_elem%RGrecomput(1:4)
          endif

          allocate(elem%face(1:max_face,elem%flen) )

          !elem%face(idx,1:elem%flen) = hp_elem%face(idx,1:elem%flen)
          do j=1,elem%flen
             elem%face(idx,j) = hp_grid%node(hp_elem%face(idx,j))%i
          enddo

          !write(*,'(a6,40i5)') &
          !     '@@@@',ie,elem%flen, elem%face(idx, 1:elem%flen)

          !elem%face(neigh,1:elem%flen) = hp_elem%face(neigh,1:elem%flen)
          !elem%face(nei_i,1:elem%flen) = hp_elem%face(nei_i,1:elem%flen)

          elem%HGnode = hp_elem%HGnode
          if(elem%HGnode) then
             allocate(elem%HGvertex(1:elem%type) )
             elem%HGvertex(1:elem%type) = hp_elem%HGvertex(1:elem%type)

             allocate(elem%HGface(1:2, 1:elem%flen) )

             elem%HGface(1, 1:elem%flen) = hp_elem%HGface(hg_face, 1:elem%flen)

             do j=1,elem%flen
                ip = 2**hp_elem%HGface(hg_level, j) + hp_elem%HGface(hg_index, j) - 1
                elem%HGface(2, j) = ip

                if(hp_elem%HGface(hg_level, j) > state%space%adapt%max_HGlevel) then
                   print*,'Too high HGlevel',ie,i,  hp_elem%i
                   do l=1,hp_elem%flen
                      print*,hp_grid%node(hp_elem%face(1,l))%x(1:nbDim),hp_elem%face(1,l)
                   enddo
                   write(*,'(a6,20i5)') 'HGver:',hp_elem%HGvertex(1:elem%type)
                   write(*,'(a6,20i5)') 'hgfac:',hp_elem%HGface(hg_face, 1:elem%flen)
                   write(*,'(a6,20i5)') 'hglev:',hp_elem%HGface(hg_level, 1:elem%flen)
                   write(*,'(a6,20i5)') 'hgind:',hp_elem%HGface(hg_index, 1:elem%flen)
                   stop
                endif

             enddo

          endif
       endif

       !if(ie==518 .or. ie==551  .or.ie==2373 .or.ie==2374 .or.ie==2375&
       !  .or.ie==2307 .or. ie==2308 .or.ie==2194.or.ie==2403.or.ie==2404.or.ie==2405) then
       !   print*,'--------------------------------'
       !   if( elem%HGnode) then
       !
       !      write(*,'(a6,20i5)') 'HGver:',hp_elem%HGvertex(1:elem%type)
       !      write(*,'(a6,20i5)') 'hgfac:',hp_elem%HGface(hg_face, 1:elem%flen)
       !      write(*,'(a6,20i5)') 'hglev:',hp_elem%HGface(hg_level, 1:elem%flen)
       !      write(*,'(a6,20i5)') 'hgind:',hp_elem%HGface(hg_index, 1:elem%flen)
       !   endif
       !      print*

       !   write(*,'(a6,20i5)') 'elem:',elem%i,elem%flen, ie
       !   write(*,'(a6,20i5)') 'face:',elem%face(1,:)
       !   if( elem%HGnode) then
       !      write(*,'(a6,20i5)') 'HGfa:',elem%HGface(1,:)
       !      write(*,'(a6,20i5)') 'HGfa:',elem%HGface(2,:)
       !      write(*,'(a6,20i5)') 'HGver:',elem%HGvertex(:)
       !   endif
       !   do l=1,elem%flen
       !      write(*,'(a6,2es14.6,i5)') 'x:',grid%x(elem%face(1,l),1:nbDim), l
       !   enddo
       !   print*,'----------- o u t ---------------------'
       !endif

    enddo

    !print*,'ie =',ie,sum(grid%elem(:)%dof)

    ! boundary faces
    grid%nbelm = hp_grid%nbelm + count(hp_grid%b_edge(1:hp_grid%nbelm)%mid_node >0)
    ie = count(hp_grid%elem(hp_grid%b_edge(1:hp_grid%nbelm)%itc)%remove)
    if(mod(ie,2) >0) then
       print*,'Trouble (4) in PasshpGrid'
       stop
    endif

    grid%nbelm = grid%nbelm - ie/2
    allocate(grid%b_edge(1:grid%nbelm))

    ip = 1
    do i=1,hp_grid%nbelm
       if(.not. hp_grid%b_edge(i)%remove) then
          if(hp_grid%b_edge(i)%mid_node >0) then ! two new segments replace a new one
             grid%b_edge(ip)%lbn(1) = hp_grid%node(hp_grid%b_edge(i)%lbn(1))%i
             grid%b_edge(ip)%lbn(2) = hp_grid%node(hp_grid%b_edge(i)%mid_node)%i
             grid%b_edge(ip)%ibc = hp_grid%b_edge(i)%ibc
             grid%b_edge(ip)%itc = hp_grid%b_edge(i)%itc
             grid%b_edge(ip)%jtc = 0 !hp_grid%b_edge(i)%jtc
             ip = ip+1

             !write(*,'(a4,5i5, 2es14.6)') &
             !     'bbb', i,ip-1, grid%b_edge(ip-1)%lbn(1:2), grid%b_edge(ip-1)%itc, &
             !     grid%x(grid%b_edge(ip-1)%lbn(1), :)

             grid%b_edge(ip)%lbn(1) = hp_grid%node(hp_grid%b_edge(i)%mid_node)%i
             grid%b_edge(ip)%lbn(2) = hp_grid%node(hp_grid%b_edge(i)%lbn(2))%i
             grid%b_edge(ip)%ibc = hp_grid%b_edge(i)%ibc
             grid%b_edge(ip)%itc = hp_grid%b_edge(i)%itc
             grid%b_edge(ip)%jtc = 0 !hp_grid%b_edge(i)%jtc
             ip = ip+1

             !write(*,'(a4,5i5, 2es14.6)') &
             !     'bbb', i,ip-1, grid%b_edge(ip-1)%lbn(1:2),  grid%b_edge(ip-1)%itc, &
             !     grid%x(grid%b_edge(ip-1)%lbn(1), :)

          elseif(hp_grid%elem(hp_grid%b_edge(i)%itc)%remove) then

             !write(*,'(a4,4i5, 2es14.6)') &
             !     'aaa', i,ip, hp_grid%node(hp_grid%b_edge(i)%lbn(1:2) )%i, &
             !     hp_grid%node(hp_grid%b_edge(i)%lbn(1) )%x(:)
             !write(*,'(a4,4i5, 2es14.6)') &
             !     'aaa', i+1,ip, hp_grid%node(hp_grid%b_edge(i+1)%lbn(1:2) )%i, &
             !     hp_grid%node(hp_grid%b_edge(i+1)%lbn(1) )%x(:)

             hp_grid%b_edge(i+1)%remove = .true.
             grid%b_edge(ip)%lbn(1) = hp_grid%node(hp_grid%b_edge(i)%lbn(1) )%i
             grid%b_edge(ip)%lbn(2) = hp_grid%node(hp_grid%b_edge(i+1)%lbn(2) )%i
             grid%b_edge(ip)%ibc = hp_grid%b_edge(i)%ibc
             !grid%b_edge(ip)%itc = hp_grid%b_edge(i)%itc
             !grid%b_edge(ip)%jtc = hp_grid%b_edge(i)%jtc

             !write(*,'(a4,4i5, 2es14.6)') &
             !     'bbb', i,ip, grid%b_edge(ip)%lbn(1:2), &
             !     grid%x(grid%b_edge(ip)%lbn(1), :)
             !write(*,'(a4,4i5, 2es14.6)') &
             !     'bbb', i,ip, grid%b_edge(ip)%lbn(1:2), &
             !     grid%x(grid%b_edge(ip)%lbn(2), :)


             ip = ip+1


          else ! original segment
             grid%b_edge(ip)%lbn(1:2) = hp_grid%node(hp_grid%b_edge(i)%lbn(1:2) )%i
             grid%b_edge(ip)%ibc = hp_grid%b_edge(i)%ibc
             grid%b_edge(ip)%itc = hp_grid%b_edge(i)%itc
             grid%b_edge(ip)%jtc = hp_grid%b_edge(i)%jtc
             ip = ip+1
          endif

       endif
    enddo


    !do i=1,grid%nelem
    !   elem => grid%elem(i)
    !   if(elem%RGhistory > 0) then
    !      do j=1,elem%RGhistory
    !         write(*,'(a3,2i5,a1,30i5)') '###',i,j,'|',elem%RG(j)%daughter(:)
    !      enddo
    !   endif
    !   if(elem%RGtype == 'D') &
    !        write(*,'(a3,2i5,a1,30i5)') '///',i,0,'|',elem%RGrecomput(1:4)
    !
    !enddo



  end subroutine PasshpGrid2Grid


  !> deallocate the arrays of hp_grid
  subroutine DeletehpGrid(hp_grid)
    type(hp_mesh), intent(inout) :: hp_grid
    type(hp_element), pointer :: elem
    integer :: i,j

    !do i=1,hp_grid%npoin
    !   deallocate ( hp_grid%node(i)%x, hp_grid%node(i)%xcur )
    !enddo
    deallocate ( hp_grid%node )

    do i=1,hp_grid%nelem
       elem => hp_grid%elem(i)

       deallocate(elem%face)
       if(allocated(elem%pnew)) deallocate(elem%pnew)

       !if(elem%RGlevel > 0) then
       if(allocated(elem%RG) ) then
          do j=1,size(elem%RG(:) )
             deallocate(elem%RG(j)%daughter)
          enddo
          deallocate(elem%RG)
       endif

       if(elem%HGnode) deallocate(elem%HGface, elem%HGvertex)
    enddo

    deallocate(hp_grid%elem )

    deallocate(hp_grid%b_edge )

  end subroutine DeletehpGrid


  !> estimates the new grid, maximal number of nodes, elements, bound_edges
  subroutine EstimateNewGrid(grid, hp_grid)
    type(mesh), intent(in) :: grid
    type(hp_mesh), intent(out) :: hp_grid
    integer :: i


    hp_grid%max_npoin = grid%npoin + count(grid%elem(:)%hsplit == 4)*3

    hp_grid%max_nelem = grid%nelem + count(grid%elem(:)%hsplit > 0)*3

    hp_grid%max_nbelm = grid%nbelm + count(grid%elem(grid%b_edge(:)%itc)%hsplit ==4)

    allocate(hp_grid%node(1:hp_grid%max_npoin) )
    allocate(hp_grid%elem(1:hp_grid%max_nelem) )
    allocate(hp_grid%b_edge(1:grid%nbelm) ) !!!

    !!print*,'New mesh:',hp_grid%max_npoin, hp_grid%max_nelem, hp_grid%max_nbelm

  end subroutine EstimateNewGrid

  !> plotting of the hp-mesh, mesh with local polynomial degrees of approximations
  subroutine PlotMesh_hp(hp_grid, isol)
    type(hp_mesh), intent(in) :: hp_grid
    integer, intent (inout)  :: isol
    character(len=15) :: solfile
    character(len=1) :: ch1
    character(len=2) :: ch2
    character(len=3) :: ch3
    character(len=4) :: ch4
    character(len=5) :: ch5
    real :: xc(1:nbDim)
    integer :: i,j, flen, lnd, chlen , ip

    integer :: detail = 0

    solfile = 'hpmesh'

    open(unit = 11, file = 'mesh.gp', status = 'replace', action = 'write')
    write(11,*) 'set terminal postscript eps'

    if(isol < 0) isol = 100000 + isol

    if(isol >= 0 .and. isol <= 9) then
       write( ch1, '(i1)' ) isol
       solfile(7:7) = ch1
       chlen = 7
    elseif(isol >= 10 .and. isol <= 99) then
       write( ch2, '(i2)' ) isol
       solfile(7:8) = ch2
       chlen = 8

    elseif(isol >= 100 .and. isol <= 999) then
       write( ch3, '(i3)' ) isol
       solfile(7:9) = ch3
       chlen = 9

    elseif(isol >= 1000 .and. isol <= 9999) then
       write( ch4, '(i4)' ) isol
       solfile(7:10) = ch4
       chlen = 10

    elseif(isol >= 10000 .and. isol <= 99999) then
       write( ch5, '(i5)' ) isol
       solfile(7:11) = ch5
       chlen = 11

    endif

    write(11,*) 'set output "',solfile(1:chlen),'.eps" '


    open(unit = 10, file = solfile, status = 'replace', action = 'write')

    ip = 0
    do i=1, hp_grid%nelem
       if(.not. hp_grid%elem(i)%remove) then
          ip = ip + 1
          xc(:) = 0.
          flen = hp_grid%elem(i)%flen
          do j = 0, flen
             lnd = hp_grid%elem(i)%face(idx,mod(j,flen)+1)

             write(10,'(2es14.6, i10)')  hp_grid%node(lnd)%x(:), lnd
             if(j < flen) xc(:) = xc(:) + hp_grid%node(lnd)%x(:)
          end do
          write(10,*) "#### elem = ",hp_grid%elem(i)%i, 'new_i = ',ip
          write(10,*)

          xc(:) = xc(:) / flen
          write(11,'(a11,i2,a5,es14.6,a3,es14.6, a7)') &
               "set label '",hp_grid%elem(i)%deg,"' at ",xc(1)," , ",xc(2), " center"
       endif
    end do
    close(10)


   if(detail ==0) then
      !write(11,*) 'plot [-0.2:1.2][-0.5:0.5] "mesh" w l'
      write(11,*) 'plot  "',solfile(1:chlen),'" w l'
   else
      write(11,*) 'set size 1.0, 0.5'
      write(11,*) 'set size ratio -1'

      write(11,*) 'set multiplot'
      write(11,*) 'set size 0.8, 0.8'

      write(11,*) 'set origin -0.25, 0.'
      write(11,*) 'unset key'
      write(11,*) 'set title "mesh"'
      write(11,*) 'plot  "',solfile,'" w l'

      write(11,*) 'set origin 0.33, 0.'
      write(11,*) 'unset key'
      write(11,*) 'set title "detail" '
      write(11,*) 'p [-0.5:1.5][-1:1]"',solfile,'" w l'

      write(11,*) 'unset multiplot'
   endif

   close(unit = 11, status = 'keep')
   !call system('gnuplot -persist "mesh.gp"')
 end subroutine PlotMesh_hp

end module hp_adaptation



