!*******************************************
!> red green refinement module
!*******************************************
module red_green
use main_data
!use structure
use  hp_adaptation

implicit none

public :: AdaptMesh_hp_RG
public :: PassGrid2hpGrid_RG
public :: RefineGrid_RG
contains

  !> hp adaptation with RG refinement (without derefinement)
  subroutine AdaptMesh_hp_RG()
    type(hp_mesh) :: hp_grid
    character(len=15) :: gridfile

    call EstimateNewGrid(old_grid, hp_grid)

    call PassGrid2hpGrid_RG(old_grid, hp_grid)

    !print*,'RefineGrid_RG'
    call RefineGrid_RG(hp_grid)

    !print*,'PasshpGrid2Grid'
    allocate(new_grid)
    call PasshpGrid2Grid(hp_grid, new_grid)

    !print*,'DeletehpGrid'
    call DeletehpGrid(hp_grid)

    !print*,'WriteMesh'
    !gridfile = 'grid0'
    !call WriteMesh(gridfile,new_grid)

    !call PlotMesh(grid, 'meshA')

    !print*,'Seeking of the neighbours'
    call new_grid%seekNeighbours()

    !print*,'Reseeking of the curved boundary'
    new_grid%curved_deg = old_grid%curved_deg
    call SeekCurvedBoundary(new_grid)

    !print*,'Passing of data from grid to gridN, deallocation of grid'
    call ReprepareProblem(new_grid, old_grid)

    old_grid => new_grid
    print*,'AdaptMesh_hp_RG -- done'

  end subroutine AdaptMesh_hp_RG


  !> passing the "type(mesh):: grid" to "type(hp_mesh) :: hp_grid"
  !> RG variant without HG nodes
  subroutine PassGrid2hpGrid_RG(grid, hp_grid)
    type(mesh), intent(in) :: grid
    type(hp_mesh), intent(inout) :: hp_grid
    class(element), pointer :: elem
    type(hp_element), pointer :: hp_elem
    integer :: i,ii, j, ip, k,k1, k2, max_flen
    integer :: l, is1, is2, ip1, ip2, ie, iie

    max_flen = 3  !4 * 2**(state%space%adapt%max_HGlevel+1)
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
    hp_grid%elem(:)%HGnode = .false.
    hp_grid%elem(:)%remove = .false.
    hp_grid%elem(:)%RGhistory = 0
    hp_grid%elem(:)%RGtype = 'N'
    hp_grid%elem(:)%flen = 3
    hp_grid%elem(:)%type = 3

    hp_grid%nelem = grid%nelem
    do i=1,grid%nelem
       elem => grid%elem(i)
       hp_elem => hp_grid%elem(i)

       hp_elem%i = i
       hp_elem%type = elem%type
       hp_elem%flen = elem%flen
       hp_elem%deg = elem%deg
       hp_elem%hsplit = elem%hsplit
       hp_elem%psplit = elem%psplit

       hp_elem%RGlevel = elem%RGlevel
       !hp_elem%RGhistory = elem%RGhistory
       !hp_elem%RGindex = elem%RGindex
       !hp_elem%RGtype = 'N'
       !hp_elem%RGreclen = 0
       hp_elem%per_boun = 0

       ! ip = hp_elem%RGhistory
       ! if(hp_elem%hsplit >0) ip = ip + 1

       ! if(ip >0) then
       !    allocate(hp_elem%RG(1: ip))

       !    do j=1,hp_elem%RGhistory
       !       hp_elem%RG(j)%subel = elem%RG(j)%subel

       !       allocate(hp_elem%RG(j)%daughter(1:hp_elem%RG(j)%subel))
       !       hp_elem%RG(j)%daughter(:) = elem%RG(j)%daughter(:)
       !    enddo

       !    if(ip > hp_elem%RGhistory) then
       !       j = hp_elem%RGhistory + 1
       !       hp_elem%RG(j)%subel = hp_elem%hsplit !! four daughter elements, red refinement

       !       allocate(hp_elem%RG(j)%daughter(1:hp_elem%RG(j)%subel))
       !       hp_elem%RG(j)%daughter(:) = 0
       !    endif

       ! endif

       allocate(hp_elem%face(idx:nei_i, 1:max_flen) )
       hp_elem%face(idx,1:hp_elem%flen) = elem%face(idx,1:hp_elem%flen)
       hp_elem%face(neigh,1:hp_elem%flen) = elem%face(neigh,1:hp_elem%flen)
       hp_elem%face(nei_i,1:hp_elem%flen) = elem%face(nei_i,1:hp_elem%flen)

       !do j=1,hp_elem%flen
       !   if(hp_elem%face(neigh,j) < 0) hp_grid%node(hp_elem%face(idx,j))%HG = .true.
       !enddo

       ! element with HanGing node(s) ?
       hp_elem%HGnode = elem%HGnode
       if(hp_elem%HGnode) then
          print*,'Troubles in redgreen.f90'
          stop
       endif

       ! if(hp_elem%HGnode) then
       !    allocate(hp_elem%HGvertex(1:hp_elem%type+1) )
       !    hp_elem%HGvertex(1:hp_elem%type) = elem%HGvertex(1:hp_elem%type)
       !    hp_elem%HGvertex(hp_elem%type+1) = hp_elem%flen+1

       !    allocate(hp_elem%HGface(1:3, 1:max_flen) )  ! max number of HG nodes
       !    ! (1, :) = index of the edge 1..type
       !    ! (2, :) = level of HG node refinement
       !    ! (3, :) = local index of HG node refinement
       !    hp_elem%HGface(hg_face, 1:hp_elem%flen) = elem%HGface(hg_face, 1:hp_elem%flen)
       !    k = 1
       !    do j=1,hp_elem%flen

       !       ip = log(1.0001*elem%HGface(2, j) ) / log(2.)

       !       hp_elem%HGface(hg_level, j) = ip
       !       hp_elem%HGface(hg_index, j) = elem%HGface(hg_level, j) - 2**ip + 1

       !       ! we seek which node is hanging one
       !       if(j == hp_elem%HGvertex(k)) then
       !          k = k+ 1
       !       else
       !          if(hp_grid%node(hp_elem%face(idx,j))%HG ) then
       !             print*,'Troubles in initiation of HG nodes'
       !             stop
       !          endif
       !          hp_grid%node(hp_elem%face(idx,j))%HG = .true.
       !       endif

       !    enddo
       ! endif

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

  end subroutine PassGrid2hpGrid_RG


  !> adaptation of the "type(hp_mesh) :: hp_grid" by Red-Green refinement
  subroutine RefineGrid_RG(hp_grid)
    type(hp_mesh), intent(inout) :: hp_grid
    type(hp_element), pointer :: elem, elem1
    type(hp_node), pointer :: node1, node2, Newnode, Newnode2
    integer :: i,j, jj, k1, k2, kk, ip1, ip2, ipc, l
    integer :: cbp1, cbp2, ibp1, ibp2, nelem_old
    logical :: periodic

    if(grid%periodic > 0) then
       print*,'Periodic grid, HG refinement not implemented !!!!!'
       !stop
    endif

    ! array for possible new nodes
    do i=1,hp_grid%nelem
       elem => hp_grid%elem(i)
       allocate( elem%pnew(1:elem%type))
       elem%pnew(1:elem%type) = 0
    enddo

    ! seeking of node for refinement,
    ! if a node is missing then a new one is included as a hanging one
    do i=1,hp_grid%nelem
       elem => hp_grid%elem(i)

       if(elem%hsplit == 4) then
          ! seeking the middle edge nodes for each edge
          do j=1,elem%type
             if( elem%pnew(j) == 0) then          ! new node has to be inserted
                jj = mod(j,elem%type) + 1

                hp_grid%npoin = hp_grid%npoin+1

                if(hp_grid%npoin > hp_grid%max_npoin ) then
                   print*,'Too much #npoin in redgreen.f90'
                endif

                Newnode => hp_grid%node(hp_grid%npoin)
                allocate(Newnode%x(1:nbDim) )

                Newnode%i = hp_grid%npoin
                Newnode%remove = .false.
                Newnode%per = 0

                ip1 = elem%face(idx,  j)
                ip2 = elem%face(idx, jj)

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


                ! adding the node to the opposite element
                if(elem%face(neigh, j) > 0) then
                   k1 = elem%face(neigh, j)
                   k2 = elem%face(nei_i, j)
                   elem1 => hp_grid%elem(k1)

                   if(periodic ) then
                      ! also the boundary nodes
                      kk = elem%per_boun
                      hp_grid%b_edge(kk)%mid_node = hp_grid%npoin

                      !print*,'###  b_edge1',elem%i, kk, hp_grid%npoin

                      ! the second node newnode2
                      hp_grid%npoin = hp_grid%npoin + 1
                      kk = elem1%per_boun
                      hp_grid%b_edge(kk)%mid_node = hp_grid%npoin

                      !print*,'###  b_edge2',elem1%i, kk, hp_grid%npoin

                   else  ! NON periodic element
                      elem1%pnew(k2) =  elem%pnew(j)
                   endif

                else ! boundary face, adding a new boundary face
                   kk = -elem%face(neigh,j)
                   hp_grid%b_edge(kk)%mid_node = hp_grid%npoin
                   !write(*,*) 'Bound face:',kk, elem%i, hp_grid%b_edge(kk)%mid_node
                endif  ! end of if( elem%pnew(j) == 0) then

             endif



             !!elem%face(neigh, k1+1) = elem%face(neigh, k1)

             !!if(elem%face(neigh, k1) >0) then
             !!   elem1%face(neigh,kk+1) = elem1%face(neigh,kk)
             !!endif

             !write(*,'(a8,4i5,2es12.4)') '!!!?',elem%i, elem%pnew(:)
             !if(elem%face(neigh,j) > 0) &
             !     write(*,'(a8,4i5,2es12.4)') '!!!*',elem1%i,elem1%pnew(:)


          enddo  !end of  do j=1,elem%type
       endif  !end of if(elem%hsplit == 4)
    enddo  !end of do i=1,hp_grid%nelem

    !print*,'@@@@ C3'

    ! mesh refinement
    nelem_old = hp_grid%nelem
    do i=1, nelem_old
       elem => hp_grid%elem(i)
       if(elem%hsplit == 4) then

          call Elem_ref_red(hp_grid, elem)

       else if(elem%hsplit > 0 .and. elem%hsplit < 4 ) then

          call Elem_ref_green(hp_grid, elem)

       endif
    enddo


    ! do i=1, hp_grid%nelem
    !    elem => hp_grid%elem(i)
    !    write(98,*) hp_grid%node(elem%face(idx,1))%x(:)
    !    write(98,*) hp_grid%node(elem%face(idx,2))%x(:)
    !    write(98,*) hp_grid%node(elem%face(idx,3))%x(:)
    !    write(98,*) hp_grid%node(elem%face(idx,1))%x(:)
    !    write(98,'(x)')
    ! enddo


    !!if(state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) stop


  end subroutine RefineGrid_RG


  !> 'red' splitting of the element without HG nodes
  subroutine Elem_ref_red(hp_grid, elem)
    type(hp_mesh), intent(inout) :: hp_grid
    type(hp_element), intent(inout) :: elem
    type(hp_element), pointer :: elem1
    integer :: j1, j2, j3, k1, k2, k3

    if(hp_grid%nelem +3 > hp_grid%max_nelem ) then
       print*,'Too much #nelem in hp_adapt.f90'
    endif

    k1 = elem%face(idx,1)
    k2 = elem%face(idx,2)
    k3 = elem%face(idx,3)

    j1 = elem%pnew(1)
    j2 = elem%pnew(2)
    j3 = elem%pnew(3)


    hp_grid%nelem = hp_grid%nelem + 1
    elem1 =>  hp_grid%elem(hp_grid%nelem)
    allocate(elem1%face(idx:nei_i, 1:3) )
    elem1%face(idx, 1) = k1
    elem1%face(idx, 2) = j1
    elem1%face(idx, 3) = j3
    elem1%deg  = elem%deg
    elem1%psplit  = elem%psplit
    elem1%i  = hp_grid%nelem
    elem1%RGlevel = elem%RGlevel + 1


    hp_grid%nelem = hp_grid%nelem + 1
    elem1 =>  hp_grid%elem(hp_grid%nelem)
    allocate(elem1%face(idx:nei_i, 1:3) )
    elem1%face(idx, 1) = k2
    elem1%face(idx, 2) = j2
    elem1%face(idx, 3) = j1
    elem1%deg  = elem%deg
    elem1%psplit  = elem%psplit
    elem1%i  = hp_grid%nelem
    elem1%RGtype = 'R'
    elem1%RGlevel = elem%RGlevel + 1

    hp_grid%nelem = hp_grid%nelem + 1
    elem1 =>  hp_grid%elem(hp_grid%nelem)
    allocate(elem1%face(idx:nei_i, 1:3) )
    elem1%face(idx, 1) = k3
    elem1%face(idx, 2) = j3
    elem1%face(idx, 3) = j2
    elem1%deg  = elem%deg
    elem1%psplit  = elem%psplit
    elem1%i  = hp_grid%nelem
    elem1%RGtype = 'R'
    elem1%RGlevel = elem%RGlevel + 1

    elem%face(idx, 1) = j1
    elem%face(idx, 2) = j2
    elem%face(idx, 3) = j3
    elem%RGtype = 'R'

  end subroutine Elem_ref_red


  !> 'green' splitting of the element without HG nodes
  subroutine Elem_ref_green(hp_grid, elem)
    type(hp_mesh), intent(inout) :: hp_grid
    type(hp_element), intent(inout) :: elem
    type(hp_element), pointer :: elem1, elem2
    integer :: jj, j1, j2, j3, k1, k2, k3

    if(hp_grid%nelem +2 > hp_grid%max_nelem ) then
       print*,'Too much #nelem in hp_adapt.f90'
    endif

    j1 = elem%hsplit
    j2 = mod(j1, 3) + 1
    j3 = mod(j2, 3) + 1

    jj = elem%pnew(elem%hsplit)

    k1 = elem%face(idx,j1)
    k2 = elem%face(idx,j2)
    k3 = elem%face(idx,j3)


    hp_grid%nelem = hp_grid%nelem + 1
    elem1 =>  hp_grid%elem(hp_grid%nelem)
    allocate(elem1%face(idx:nei_i, 1:3) )
    elem1%face(idx, 1) = k1
    elem1%face(idx, 2) = jj
    elem1%face(idx, 3) = k3
    elem1%deg  = elem%deg
    elem1%psplit  = elem%psplit
    elem1%i  = hp_grid%nelem
    elem1%RGtype = 'G'
    elem1%RGlevel = elem%RGlevel + 1

    elem%face(idx, 1) = jj
    elem%face(idx, 2) = k2
    elem%face(idx, 3) = k3
    elem%RGtype = 'G'

  end subroutine Elem_ref_green

end module red_green
