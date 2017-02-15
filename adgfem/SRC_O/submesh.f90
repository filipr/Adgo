!> create a conforming sub-triangulation by splitting of elements having hanging nodes
module submesh
  use paramets
  use ama_L2interpol
  use data_mod ! contains type(mesh) :: gridN for computation
  use mesh_oper


  implicit none

  integer, parameter,private :: max_subelems = 4  ! maximal number of subtriangles of the given elem
  integer, dimension(:,:), allocatable :: submes_idx ! index of sub-elems forming original elem


  public :: Create_Conforming_Subgrid
  public :: Set_Sim_elems
  public :: Set_Sim_nelem
  public :: Set_estim_from_Subgrid
contains

  !> create a conforming sub-triangulation grid_Sim by splitting of elements of grid_HG
  !> having hanging nodes
  subroutine Create_Conforming_Subgrid(grid_HG, grid_Sim)
    class(mesh), intent(in) :: grid_HG
    class(mesh), pointer, intent(out) :: grid_Sim
    class(element), pointer :: elem, elemHG
    integer :: i, j, k, nelem, Tdeg

    allocate(grid_Sim)

    ! estimate of the number of elements of the submesh
    call Set_Sim_nelem(grid_HG, nelem)

    grid_Sim%nelem = nelem
    grid_Sim%npoin = grid_HG%npoin  ! number of points rests the same
    grid_Sim%nbelm = grid_HG%nbelm  ! number of boundary elements rests the same
    grid_Sim%nbc   = grid_HG%nbc    !

    !print*,'343d', size(grid_Sim%xper, 1)

    grid_Sim%xper(:,:) = grid_HG%xper(:,:)
    grid_Sim%iper(:,:) = grid_HG%iper(:,:)

    ! mesh nodes
    allocate(grid_Sim%x(1:grid_Sim%npoin, 1:nbDim) )
    allocate(grid_Sim%xcur(1:grid_Sim%npoin, 1:nbDim) )
    grid_Sim%x(1:grid_Sim%npoin, 1:nbDim) = grid_HG%x(1:grid_Sim%npoin, 1:nbDim)


    ! elements
    call Set_Sim_elems(grid_HG, grid_Sim)

    ! boundary segments
    allocate( grid_Sim%b_edge( 1:grid_Sim%nbelm ) )
    do i=1,grid_Sim%nbelm
       grid_Sim%b_edge(i)%lbn(1:2) = grid_HG%b_edge(i)%lbn(1:2)
       grid_Sim%b_edge(i)%ibc      = grid_HG%b_edge(i)%ibc
    enddo

    grid_Sim%ElemSupports  = .false.


    ! seeking of neighbours
    call grid_Sim%seekNeighbours( )

    ! setting of elem%ibc
    call SetConstBC(grid_Sim, .true.)


    ! courved boundaries
    if (nbDim==2) then
       !!!call ReadCurvedBoundary(prof_file)
       call SeekCurvedBoundary(grid_Sim)
    endif

    call grid_sim%computeGeometry( )


    ! call SetEdgeDegrees(grid_Sim)
    call grid_Sim%setEdgeQuadratureDegrees()

    do i=1, grid_Sim%nelem
       elem => grid_Sim%elem(i)
       call SetElementQuadraturesDegrees(elem )
       call PrepareOneElement(elem )
       call ComputeLocalMassMatrix(elem)

       allocate( elem%wSS(1:1, 0:2, 1:elem%dof*ndim) )
       !print*,'??????????',size(elem%wSS, 1),size(elem%wSS, 2),size(elem%wSS, 3)
       !print*,'??????????',elem%i, elem%iBC(:)
    enddo

    !call grid_Sim%plot('mesh00')

    !print*,' recomputation of the solution  only in elem%ws(0:1, :)'
    !Tdeg = state%time%deg
    !state%time%deg = -1
    call AdvancedInterpolDGsolution(grid_Sim, grid_HG, 2 )
    !state%time%deg = Tdeg
    !print*,' recomputation of the solution  only in elem%ws(0:1, :)'

    do i=1, grid_Sim%nelem
       elem => grid_Sim%elem(i)
       allocate(elem%w(0:0, 1:elem%dof *ndim) )
       elem%w(0, 1:elem%dof *ndim) = elem%wSS(1, 0, 1:elem%dof *ndim)  ! elem%wSS(1,0,:) is identical
    enddo


    !do i=1, grid_Sim%nelem
    !   elem => grid_Sim%elem(i)
    !   call PlotElemFunction3D(100+state%space%adapt%adapt_level, elem,  elem%dof, elem%wh(0, 1:elem%dof))
    !enddo

    !if(state%space%adapt%adapt_level == 6) stop 'dedt37d3hw2ouis'

  end subroutine Create_Conforming_Subgrid


  !> setting of elements of the subgrid
  subroutine Set_Sim_elems(grid_HG, grid_Sim)
    class(mesh), intent(in) :: grid_HG
    class(mesh), pointer, intent(out) :: grid_Sim
    class(element), pointer :: elem, elem_HG
    integer :: i, j, k, nelem, ie
    integer :: k1, k2, k3, khg, khg1,khg2, khg3, i1, i2, i3


    allocate(grid_Sim%elem(1:grid_Sim%nelem) )

    !allocations of the face array
    do i=1, grid_Sim%nelem
       elem => grid_Sim%elem(i)
       elem%i = i
       elem%type = 3   ! the resulting subrid is simplicial
       elem%flen = 3

       allocate( elem%face(1:max_face, 1:elem%flen) )
    enddo


    allocate(submes_idx(1:grid_HG%nelem, 0:max_subelems) )
    submes_idx(:,:) = 0

    ! creation of the grid_Sim from grid_HG
    ie = 0
    do i=1, grid_HG%nelem
       elem_HG => grid_HG%elem(i)

       !print*,'###',i, elem_HG%i, elem_HG%HGnode

       if(elem_HG%HGnode ) then
          !print*, 'HANGING nodes to be DONE'


          if(elem_HG%flen == 4) then
             ! one HG node on a one edge, k1,k2,k3 - vertexes, khg - hanging node

             if( elem_HG%HGvertex(2) - elem_HG%HGvertex(1) == 2) then
                i1 = 1
             elseif( elem_HG%HGvertex(3) - elem_HG%HGvertex(2) == 2) then
                i1 = 2
             else
                i1 = 3
             endif
             i2 = mod(i1, 3 ) + 1
             i3 = mod(i2, 3 ) + 1

             k1 = elem_HG%face(idx, elem_HG%HGvertex(i1) )
             k2 = elem_HG%face(idx, elem_HG%HGvertex(i2) )
             k3 = elem_HG%face(idx, elem_HG%HGvertex(i3) )

             do k=1, elem_HG%flen
                khg = elem_HG%face(idx, k)
                if(khg /= k1 .and. khg /= k2 .and. khg /= k3 ) goto 10
             enddo
10           continue

             !write(*,*) grid_Sim%x(k1, 1:2),k1
             !write(*,*) grid_Sim%x(k2, 1:2),k2
             !write(*,*) grid_Sim%x(k3, 1:2),k3
             !write(*,*) grid_Sim%x(khg, 1:2),khg, elem_HG%face(idx, :)
             !write(*,*) '##', elem_HG%HGvertex(:)
             !stop

             ! two elements
             do j=1, 2
                ie = ie + 1
                submes_idx(i,0) = submes_idx(i,0) + 1
                k = submes_idx(i,0)
                submes_idx(i,k) = ie

                elem =>  grid_Sim%elem(ie)
                if(j == 1) then
                   elem%face(idx, 1) = k1;   elem%face(idx, 2) = khg;  elem%face(idx, 3) = k3;

                elseif( j== 2) then
                   elem%face(idx, 1) = k3;   elem%face(idx, 2) = khg;   elem%face(idx, 3) = k2
                endif

                elem%deg = elem_HG%deg
                call elem%initElementDof()

                ! solution will be set later
                !allocate(elem%w(0:0, 1:ndim*elem%dof) )
                allocate(elem%wS(0:1, 1:ndim*elem%dof) )
             enddo

          elseif(elem_HG%flen == 5) then
             ! two HG node on a one edge, k1,k2,k3 - vertexes, khg, khg1 - hanging nodes

             if( elem_HG%HGvertex(2) - elem_HG%HGvertex(1) == 1) then
                i1 = 1
             elseif( elem_HG%HGvertex(3) - elem_HG%HGvertex(2) == 1) then
                i1 = 2
             else
                i1 = 3
             endif
             i2 = mod(i1, 3 ) + 1
             i3 = mod(i2, 3 ) + 1

             k1 = elem_HG%face(idx, elem_HG%HGvertex(i1) )
             k2 = elem_HG%face(idx, elem_HG%HGvertex(i2) )
             k3 = elem_HG%face(idx, elem_HG%HGvertex(i3) )

             do k=1, elem_HG%flen
                if(k == i2) khg  = elem_HG%face(idx,  elem_HG%HGvertex(i2)+1)
                if(k == i3) khg1 = elem_HG%face(idx,  elem_HG%HGvertex(i3)+1)
             enddo

             !write(*,*) grid_Sim%x(k1, 1:2),k1,i1
             !write(*,*) grid_Sim%x(k2, 1:2),k2,i2
             !write(*,*) grid_Sim%x(k3, 1:2),k3,i3
             !write(*,*) grid_Sim%x(khg, 1:2),khg, elem_HG%face(idx, :)
             !write(*,*) grid_Sim%x(khg1, 1:2),khg1
             !write(*,*) '##', elem_HG%HGvertex(:)
             !stop

             ! 3 elements
             do j=1, 3
                ie = ie + 1
                submes_idx(i,0) = submes_idx(i,0) + 1
                k = submes_idx(i,0)
                submes_idx(i,k) = ie

                elem =>  grid_Sim%elem(ie)
                if(j== 1) then
                   elem%face(idx, 1) = k1 ;  elem%face(idx, 2) = k2;   elem%face(idx, 3) = khg

                elseif(j== 2) then
                   elem%face(idx, 1) = k1 ;  elem%face(idx, 2) = khg;  elem%face(idx, 3) = khg1

                elseif(j== 3) then
                   elem%face(idx, 1) = k3 ;  elem%face(idx, 2) = khg1; elem%face(idx, 3) = khg
                endif

                elem%deg = elem_HG%deg
                call elem%initElementDof()

                ! solution will be set later
                !allocate(elem%w(0:0, 1:ndim*elem%dof) )
                allocate(elem%wS(0:1, 1:ndim*elem%dof) )
             enddo

          elseif(elem_HG%flen == 6) then
             ! 3 HG nodes on a esch edge, k1,k2,k3 - vertexes, khg, khg1 - hanging nodes

             k1   = elem_HG%face(idx,1)
             khg1 = elem_HG%face(idx,2)
             k2   = elem_HG%face(idx,3)
             khg2 = elem_HG%face(idx,4)
             k3   = elem_HG%face(idx,5)
             khg3 = elem_HG%face(idx,6)


             !write(*,*) grid_Sim%x(k1, 1:2),k1,i1
             !write(*,*) grid_Sim%x(k2, 1:2),k2,i2
             !write(*,*) grid_Sim%x(k3, 1:2),k3,i3
             !write(*,*) grid_Sim%x(khg1, 1:2),khg1
             !write(*,*) grid_Sim%x(khg2, 1:2),khg2
             !write(*,*) grid_Sim%x(khg3, 1:2),khg3
             !write(*,*) elem_HG%face(idx, :)
             !write(*,*) '##', elem_HG%HGvertex(:)
             !stop 'de38ewis'

             ! four elements
             do j=1, 4
                ie = ie + 1
                submes_idx(i,0) = submes_idx(i,0) + 1
                k = submes_idx(i,0)
                submes_idx(i,k) = ie

                elem =>  grid_Sim%elem(ie)
                if(j == 1) then
                   elem%face(idx, 1) = k1;   elem%face(idx, 2) = khg1;  elem%face(idx, 3) = khg3;

                elseif(j == 2) then
                   elem%face(idx, 1) = k2;   elem%face(idx, 2) = khg2;  elem%face(idx, 3) = khg1;

                elseif(j == 3) then
                   elem%face(idx, 1) = k3;   elem%face(idx, 2) = khg3;  elem%face(idx, 3) = khg2;

                elseif( j== 4) then
                   elem%face(idx, 1) = khg1; elem%face(idx, 2) = khg2;  elem%face(idx, 3) = khg3;
                endif

                elem%deg = elem_HG%deg
                call elem%initElementDof()

                ! solution will be set later
                !allocate(elem%w(0:0, 1:ndim*elem%dof) )
                allocate(elem%wS(0:1, 1:ndim*elem%dof) )
             enddo


          else
             print*,'elem_HG%flen=', elem_HG%flen, elem_HG%xc
             stop 'OTHER case not implemented #####'
          endif

       else ! no HG node
          ie = ie + 1
          elem => grid_Sim%elem(ie)
          if(elem_HG%flen > 3) stop 'Troubles r4f47hj'

          ! sub-elem froming elem_HG
          submes_idx(i,0) = submes_idx(i,0) + 1
          k = submes_idx(i,0)
          submes_idx(i,k) = ie

          ! copying of the element
          elem%face(idx, 1:3) = elem_HG%face(idx, 1:3)
          elem%deg = elem_HG%deg
          call elem%initElementDof()

          ! solution
          !allocate(elem%w(0:0, 1:ndim*elem%dof) )
          allocate(elem%wS(0:1, 1:ndim*elem%dof) )
          !elem%w(0, 1:ndim*elem%dof) = elem_HG%w(0, 1:ndim*elem%dof)
       endif
    enddo


  end subroutine Set_Sim_elems

  !> compute the number of elements of the simplicial submesh
  subroutine Set_Sim_nelem(grid_HG, nelem)
    class(mesh), intent(in) :: grid_HG
    class(element), pointer :: elem
    integer :: i, j, k, nelem

    nelem = 0
    do i = 1, grid_HG%nelem
       elem => grid_HG%elem(i)
       if(elem%HGnode) then
          k = 0
          do j=1, elem%flen
             if( elem%HGface(2, j) == 2 .or. elem%HGface(2, j) == 3) then
                ! edge with one HG node
                 k = k+1
              elseif( elem%HGface(2, j) >  3) then
                 stop "more than one HG node per edge, not implemented"
              endif
           enddo
           if(k == 0) then
              nelem = nelem + 1
              print*,' streng, HG node?'

           elseif(k == 2) then
              nelem = nelem + 2

           elseif(k == 4) then
              nelem = nelem + 3

           elseif(k == 6) then
              nelem = nelem + 4
           else
              stop " some troubles????"
           endif

           !print*,'verify', k

       else
          nelem = nelem + 1
       endif
    end do

    print*,' ## The simplicial grid has nelem = ', nelem

  end subroutine Set_Sim_nelem


  !> setting of the error estims and regul parameters from grid_Sim to grid_HG
  subroutine Set_estim_from_Subgrid(grid_HG, grid_Sim)
    class(mesh), intent(in) :: grid_HG
    class(mesh), pointer, intent(out) :: grid_Sim
    class(element), pointer :: elem, elem_HG
    integer :: i, j, ie


    do i=1, grid_HG%nelem
       elem_HG => grid_HG%elem(i)

       elem_HG%eta(:,:) = 0.

       do j=1, submes_idx(i,0)
          ie = submes_idx(i,j)
          elem => grid_Sim%elem(ie)

          elem_HG%eta(:,:) = elem_HG%eta(:,:) + elem%eta(:,:)**2
       enddo

       elem_HG%eta(:,:) = sqrt(elem_HG%eta(:,:))
    enddo

  end subroutine Set_estim_from_Subgrid


  !> deallocation of the simplicial submesh
  subroutine Delete_Submesh(grid_Sim)
    class(mesh), pointer, intent(out) :: grid_Sim
    class(element), pointer :: elem
    integer :: i, j

    do i=1, grid_Sim%nelem
       elem => grid_Sim%elem(i)

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
          deallocate(elem%nc, elem%dnc)

       endif

       if(state%MGsolver .and. max_MG > 0) deallocate(elem%MGw)

       !do j=1,elem%flen
       !if(elem%face(neigh, j) <= 0) then
       deallocate(elem%iBC)
       deallocate(elem%tBC)
       !endif
       !enddo

       deallocate( elem%xc)
       deallocate(elem%F)

       !print*,'DeallocateElem, elem%i =', elem%i, 'type Z6'
       deallocate(elem%n)
       !print*,'DeallocateElem, elem%i =', elem%i, 'type Z7'
       deallocate(elem%dn)

       deallocate(elem%face)
       deallocate(elem%w)
       !deallocate(elem%wS)
       deallocate(elem%wSS)

       deallocate( elem%eta)
       deallocate(elem%RTNflux )

       deallocate(elem%lifting )
       deallocate(elem%xi)

       deallocate(elem%Mass)
       deallocate(elem%MassInv)
       deallocate(elem%Stiff)

    end do

    deallocate(grid_Sim%elem)

    deallocate(grid_Sim%x)
    deallocate(grid_Sim%xcur)

    deallocate(submes_idx)

  end subroutine Delete_Submesh

end module submesh
