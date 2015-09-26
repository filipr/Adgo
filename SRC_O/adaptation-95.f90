!> mesh adaptation
module mesh_adapt95
  use main_data  ! contains "type(mesh) ::  grid"   for computation
  use structure  ! contains "type(Amesh):: Agrid"   for adaptation
  use init
  use angen

  public :: Amesh2Mesh
  public :: AdaptMesh_AMA
  public :: Mesh2Amesh

contains

  !> read the mesh form gridfile in the "grid" file
  !subroutine AdaptMesh(Agrid)
  ! NOT USED ALREADY *.f95
  subroutine AdaptMesh_AMA( )
    type(Amesh), pointer :: Agrid
    type(Omesh), pointer :: Ogrid

    allocate(Agrid)
    allocate(Ogrid)

    if(grid%max_el_len /= 3) then
       print*,'AMA for triangles without HG nodes, grid%max_el_len = ',grid%max_el_len
       stop
    endif

    print*,'Mesh adaptation will be done'

    print*,'Transformatin grid ==> Agrid'
    call Mesh2Amesh(grid, Agrid)

    print*,'Start of ANGENER'
    call ANGENER(Agrid, Ogrid)

 !   print*,'Plotting of the mesh'
  !  call PlotAmesh(Agrid)

    allocate(gridN)
    print*,'Transformatin Agrid ==> grid'
    call Amesh2Mesh(Agrid, gridN)

    !print*,'Deleting Agrid'
    call RmvObj(Agrid)

    !! 1st neighbours, 2nd curved
    call SeekNeighbours(gridN)

    gridN%curved_deg = grid%curved_deg
    call SeekCurvedBoundary(gridN)

    print*,'Passing of data from grid to gridN, deallocation of grid'
    call ReprepareProblem(gridN, grid)

    grid => gridN

  end subroutine AdaptMesh_AMA


  !> transformation of type(mesh):: grid to type(Amesh):: Agrid
  subroutine Mesh2Amesh(grid,Agrid)
    type(mesh), intent(in) :: grid
    type(Amesh), pointer   :: Agrid

    integer :: prmt  = 50
    integer :: prfl  = 51

    type(node), pointer      :: CrntN   => null()
    type(edge), pointer      :: CrntEd  => null()
    type(Aelement), pointer  :: CrntE   => null()
    type(neighbour), pointer :: LastE   => null()

    real(kind=dp), dimension(2) :: v1
    real(kind=dp), dimension(2) :: v2

    integer, dimension(3) :: node_idx
    integer :: i
    integer :: j
    integer :: k
    integer :: pibc

    Agrid%npoin = grid%npoin
    Agrid%nelem = grid%nelem
    Agrid%nbelm = grid%nbelm

    Agrid%xper(:,:) = grid%xper(:,:)
    Agrid%iper(:,:) = grid%iper(:,:)
    Agrid%periodic = grid%periodic

    Agrid%ndim = ndim
    Agrid%numel = 200
    Agrid%ityp = 1
    Agrid%ifv  = 2
    Agrid%pos  = 0.02
    Agrid%p    = 10000.
    !Agrid%p    = 10.
    Agrid%eps  = 1E+15

    !creating the DLL of nodes:
    !  - setting the index of node
    !  - setting the coordinates
    !  - presuming each node inner
    !  - at last setting the actually highest index of vertices to 'npoin'
    do i = 1,Agrid%npoin
       call AddObj(CrntN)
       CrntN%ind = i
       CrntN%crd(1:nbDim) = grid%x(i,1:nbDim)
       if (i == 1) Agrid%nodes => CrntN
    end do
    Agrid%iNode = Agrid%npoin
    Agrid%LastNd => CrntN

    !creating the DLL of elements:
    !  - setting the index of Aelement
    !  - setting the pointers to the 3 vertices (by searching for the index
    !    of the node in DLL)
    !  - right away inserting the given Aelement into the vertex's DLL of Aelements
    !    containing the particular vertex (pointer 'els' always pointing to the
    !    last object!)
    !  - setting the INE-pointer to the place in the node's 'els' list
    !  - at last setting the actually highest index of Aelements to 'nelem' and
    !    setting the pointer to the last Aelement
    do i = 1,Agrid%nelem
       call AddObj(CrntE)
       CrntE%ind = i
       if (i == 1) Agrid%elmnts => CrntE
       if(grid%elem(i)%flen /= 3) then
          print *,'Only triangles without hanging nodes implemented in "Mesh2Amesh"'
          stop
       else
          node_idx(1:3) = grid%elem(i)%face(idx,1:3)
       end if

       do j = 1,3
          CrntE%vtx(j)%p => SearchNode(Agrid,node_idx(j))
          call AddObj(CrntE%vtx(j)%p%els)
          CrntE%vtx(j)%p%els%elem => CrntE
          CrntE%INE(j)%p => CrntE%vtx(j)%p%els
       end do
    end do
    Agrid%iElem = Agrid%nelem
    Agrid%LastEl => CrntE

    !creating the DLL of boundary segments:
    !  - setting the two pointers to the nodes forming the segment
    !  - setting the first node's 'edg' pointer on this edge
    !  - setting the boundary component variable 'ibc' for edge
    !  - setting inner = .false. for the particular nodes
    pibc = 0 !> previous ibc
    do i = 1,Agrid%nbelm
       call AddObj(CrntEd)
       if (i == 1) Agrid%edges => CrntEd
       node_idx(1:nbDim) = grid%b_edge(i)%lbn(1:nbDim)
       node_idx(3) = grid%b_edge(i)%ibc
       do j =1,2
          CrntN => SearchNode(Agrid,node_idx(j))
          CrntEd%nds(j)%p => CrntN
          if (j == 1) CrntN%edg => CrntEd
       end do
       CrntEd%ibc = node_idx(3)
       !> new boundary component, first node should be fixed:
       if (node_idx(3) /= pibc) CrntEd%nds(1)%p%fixed = .true.
       do j = 1,2
          if (CrntEd%nds(j)%p%inner) CrntEd%nds(j)%p%inner = .false.
       end do
       pibc = node_idx(3)
    end do

    !  - rewinding each node's 'els' pointer
    !  - reordering each node's 'els' DLL so that elements are counterclockwise and
    !  - setting each inner node's 'els' into a circle
    CrntN => Agrid%nodes
    do
       if (.not.associated(CrntN)) exit
       do
          if (.not.associated(CrntN%els%prv)) exit
          CrntN%els => CrntN%els%prv
       end do

       call ReorderEls(CrntN)

       if (CrntN%inner) then
          LastE => CrntN%els
          do
             if (.not.associated(LastE%nxt)) exit
             LastE => LastE%nxt
          end do
          CrntN%els%prv => LastE
          LastE%nxt     => CrntN%els
       end if
       CrntN => CrntN%nxt
    end do

    !setting each element's DLL of neighbours
    CrntE => Agrid%elmnts
    do
       if (.not.associated(CrntE)) exit
       call SetNei(CrntE)
       CrntE => CrntE%nxt
    end do

    call InitNCB(Agrid)
    call InitPer(Agrid)
    !> - setting nodes fixed at right angles
    CrntEd => Agrid%edges
    do
       if (.not.associated(CrntEd)) exit
       if (.not.CrntEd%curved .and. associated(CrntEd%nxt)) then
          v1 = CrntEd%nds(1)%p%crd - CrntEd%nds(2)%p%crd
          v2 = CrntEd%nxt%nds(2)%p%crd - CrntEd%nxt%nds(1)%p%crd
          if (abs(dot_product(v1,v2)) .le. 1.d+15) CrntEd%nds(2)%p%fixed = .true.
       end if
       CrntEd => CrntEd%nxt
    end do

    i = Agrid%nelem
    j = Agrid%ndim
    allocate(Agrid%w(i,j))
    do k = 1,i
       call Eval_aver_w_Elem(grid%elem(k), Agrid%w(k,1:j))
       write(700,'(4i5,4es12.4)' ) k, grid%elem(k)%face(idx,1:3) ,Agrid%w(k,:)
    end do
  end subroutine Mesh2Amesh

  !> transformation of 'type(Amesh):: Agrid' to 'type(mesh):: grid'
  subroutine Amesh2Mesh(Agrid, grid)
    type(Amesh), intent(in) :: Agrid
    type(mesh), intent(out) :: grid

    type(node), pointer      :: CrntN   => null()
    type(edge), pointer      :: CrntEd  => null()
    type(Aelement), pointer  :: CrntE   => null()
    type(neighbour), pointer :: LastE   => null()
    class(element), pointer   :: elem

    integer, dimension(3) :: node_idx
    integer :: i
    integer :: j
    integer :: k


    !passing agruments from Agrid to grid

    ! counting of mesh nodes
   CrntN => Agrid%nodes
   i = 1
   do
      if (.not.associated(CrntN)) then
         exit
      end if
      CrntN%ind = i
      i = i+ 1
      CrntN => CrntN%nxt
   enddo
   grid%npoin = i-1


   ! counting of mesh elements
   CrntE => Agrid%elmnts
   i = 1
   do
      if (.not.associated(CrntE)) then
         exit
      end if
      i = i+ 1
      CrntE => CrntE%nxt
   enddo
   grid%nelem = i - 1


   ! counting of boundary segments
   CrntEd => Agrid%edges
   i = 1
   do
      if (.not.associated(CrntEd)) then
         exit
      end if
      i = i+ 1
      CrntEd => CrntEd%nxt
   enddo
   grid%nbelm = i - 1

   grid%xper(:,:) = Agrid%xper(:,:)
   grid%iper(:,:) = Agrid%iper(:,:)
   grid%periodic = Agrid%periodic

   !! nodes
   allocate(grid%x(grid%npoin,2) )
   allocate(grid%xcur(1:grid%npoin, 1:nbDim))

   CrntN => Agrid%nodes
   i = 1
   do
      if (.not.associated(CrntN)) then
         exit
      end if
      grid%x(i,1:nbDim) = CrntN%crd(1:nbDim)
      !write(111,*) grid%x(i,1:nbDim)
      i = i+ 1
      CrntN => CrntN%nxt
   enddo

   ! elements
   allocate(grid%elem(1:grid%nelem) )
   CrntE => Agrid%elmnts
   i = 1
   do
      if (.not.associated(CrntE)) then
         exit
      end if
      elem => grid%elem(i)
      elem%flen = 3
      elem%type = 3
      elem%HGnode = .false.

      allocate(elem%face(1:max_face, 1:elem%flen))

      elem%face(idx,1) = CrntE%vtx(1)%p%ind
      elem%face(idx,2) = CrntE%vtx(2)%p%ind
      elem%face(idx,3) = CrntE%vtx(3)%p%ind


      if(elem%face(idx,1) > grid%npoin .or.  &
           elem%face(idx,2) > grid%npoin .or.  &
           elem%face(idx,3) > grid%npoin ) then
         write(*,'(a7,6i6)') '####',i,grid%npoin, elem%face(idx,1:3)
      endif

      elem%i = i
      elem%deg = state%space%deg
      call InitElemDof( elem )

      elem%hsplit = 0
      elem%psplit = 0
      !elem%RGtype = CrntE%RGtype
      !elem%RGlevel = CrntE%RGlevel
      !elem%RGindex = CrntE%RGindex
      !elem%RGmother = CrntE%RGmother
      !elem%RGdaughter(:)  = CrntE%RGdaughter(:)

      allocate(elem%w(0:state%time%deg+1, 1:ndim* elem%dof) )

      call SetOneElementConstantIC(elem, Agrid%wp(i, 1:ndim) )

      !write(*,'(a6,i5,12es12.4)') '<< 10',elem%i, elem%w(0,1:3)
      !do j=0,3
      !   write(112,*) grid%x(elem%face(idx,mod(j,3)+1), 1:nbDim)
      !enddo
      !write(112,*)

      ! taking values of the solution
      allocate(grid%elem(i)%wS(1:1,1:ndim))

      ! passiong of the solution from the Amesh to mesh
      do j=1,ndim
      ! HERE to insert AMA
         !grid%elem(i)%wS(1, j) =  CrntE%w?????
         grid%elem(i)%wS(1, j) = state%BC(1)%ww(j) ! doscasne, jen aby to probehlo
      enddo

      i = i+ 1
      CrntE => CrntE%nxt
   enddo



   ! boundary segments
   allocate(grid%b_edge(1:grid%nbelm) )
   CrntEd => Agrid%edges
   i = 1

   do
      if (.not.associated(CrntEd)) then
         exit
      end if

      if(i > grid%nbelm) then
         write(*,*) 'Troubles in Amesh2Mesh in adaptation.f90'
         write(*,*) CrntEd%nds(1)%p%crd(1:nbDim), i
         write(*,*) CrntEd%nds(2)%p%crd(1:nbDim), grid%nbelm
         write(*,*)
      else
         grid%b_edge(i)%lbn(1) = CrntEd%nds(1)%p%ind
         grid%b_edge(i)%lbn(2) = CrntEd%nds(2)%p%ind
         grid%b_edge(i)%ibc = CrntEd%ibc

         !write(100+state%space%adapt%adapt_level,*) CrntEd%nds(1)%p%crd(1:nbDim), i
         !write(100+state%space%adapt%adapt_level,*) CrntEd%nds(2)%p%crd(1:nbDim), grid%nelem
         !write(100+state%space%adapt%adapt_level,*)
      endif
      i = i+ 1
      CrntEd => CrntEd%nxt
   enddo

 end subroutine Amesh2Mesh

 !> the  Agrid in Amesh type
 subroutine PlotAmesh1(Agrid)
   type(Amesh), pointer :: Agrid
   type(Aelement), pointer :: elmnt
   integer :: i
   integer :: detail = 0

   open(unit = 10, file = 'mesh', status = 'replace', action = 'write')
   elmnt => Agrid%elmnts
   do
      if (.not.associated(elmnt)) exit
      do i = 0,3
         write(10,'(2es14.6, i10)') elmnt%vtx(mod(i,3)+1)%p%crd,elmnt%vtx(mod(i,3)+1)%p%ind
      end do
      write(10,*) "#### elem = ", elmnt%ind
      write(10,*)
      elmnt => elmnt%nxt
   end do
   close(unit = 10, status = 'keep')

   open(unit = 11, file = 'mesh.gp', status = 'replace', action = 'write')
   write(11,*) 'set terminal postscript eps'
   write(11,*) 'set output "mesh.eps" '

   if(detail ==0) then
      !write(11,*) 'plot [-0.2:1.2][-0.5:0.5] "mesh" w l'
      write(11,*) 'plot  "mesh" w l'
   else
      write(11,*) 'set size 1.0, 0.5'
      write(11,*) 'set size ratio -1'

      write(11,*) 'set multiplot'
      write(11,*) 'set size 0.8, 0.8'

      write(11,*) 'set origin -0.25, 0.'
      write(11,*) 'unset key'
      write(11,*) 'set title "mesh"'
      write(11,*) 'p "mesh" w l'

      write(11,*) 'set origin 0.33, 0.'
      write(11,*) 'unset key'
      write(11,*) 'set title "detail" '
      write(11,*) 'p [-0.5:1.5][-1:1] "mesh"  w l'

      write(11,*) 'unset multiplot'
   endif

   close(unit = 11, status = 'keep')
   call system('gnuplot -persist "mesh.gp"')

 end subroutine PlotAmesh1

end module mesh_adapt95
