!> mesh operations: reading, plotting, computing geometry
module mesh_oper
!  use matrix_oper
  use mesh_mod
  use main_data  !, g_grid => grid !mesh operation IS NOT USED for global grid (adaptation)
  use f_mapping
  use geometry, only: GetFaceIndexes
!  use mesh_oper3D

  implicit none

 ! public:: ReadMesh
  public:: WriteMesh
!  public:: PlotMesh
  public:: SeekCurvedNodeGlobal
  public:: SeekCurvedNodeLocal
  public:: SeekCurvedBoundary
  public:: SplineCurvedBoundary
!  public:: SeekNeighbours
  !public:: SetCurved         ! seeking of element adjacent to curved boundary

 ! public:: ComputeGeometry
 ! public:: ComputeElementGeometry

  public:: ReadCurvedBoundary

  public:: WriteMeshMedit
  public:: WriteMeshMatlab

  public:: SeekElemSupports

  public :: FindTriangleCoords
  public :: FindTriangleIndexs
  public :: BarycCoordOne


contains
!  !> read the mesh from gridfile in the "grid" file
!  subroutine ReadMesh(gridfile, grid)
!    character(len=*), intent(in) :: gridfile
!    type(mesh), intent(inout), target :: grid
!    class(element), pointer :: elem
!    integer :: ifile=12
!    integer, parameter :: maxlen=20;
!    integer:: i, lenloc, iloc(1:maxlen), counter
!    real, dimension(1:3) :: vec1, vec2, vec3
!    real, dimension(1:2, 1:nbDim) :: xminmax
!    real :: tmp
!
!    open(ifile, file=gridfile, status='OLD')
!    read(ifile,*) grid%npoin, grid%nelem, grid%nbelm, grid%nbc
!    read(ifile,*) grid%xper(1,1:nbDim),grid%iper(1,1:2),&
!          grid%xper(2, 1:nbDim),grid%iper(2,1:2)
!
!    ! periodicity
!    grid%periodic = 0
!    if(grid%xper(1,1) /= 0. .or. grid%xper(1,2) /= 0.) then
!       grid%periodic = 1
!       if(grid%xper(2,1) /= 0. .or. grid%xper(2,2) /= 0.) grid%periodic = 2
!    endif
!    if(nbDim == 3) print*,'Attention for periodicity in 3D in mesh.f90 !!!!'
!
!    allocate(grid%x(1:grid%npoin,1:nbDim) )
!    allocate(grid%xcur(1:grid%npoin,1:2) )  ! idx = 1 CBP, idx =2 pnt
!
!    do i=1,grid%npoin
!       read(ifile,*) grid%x(i,1:nbDim)
!    enddo
!
!    do i=1,nbDim
!       xminmax(1, i) = minval(grid%x(:, i) )
!       xminmax(2, i) = maxval(grid%x(:, i) )
!    end do
!
!    grid%diam = VectorNorm( xminmax(2,:) - xminmax(1,:) )
!
!    !a0 = loc(grid)
!    allocate( grid%elem(1:grid%nelem) )
!
!    !a1 = loc(grid%elem(1))
!    !a2 = loc(grid%elem(2))
!    !a4 = loc(grid%elem(grid%nelem))
!
!    !print*,'@@@@',a0,a1,a2, a4
!    !print*, a2-a1, (a4-a1)/(grid%nelem-1), a4-a1
!    !stop
!
!    do i=1, grid%nelem
!       elem => grid%elem(i)
!       elem%i = i
!
!       read(ifile, *) lenloc, iloc(1:lenloc)
!       if(abs(lenloc) > maxlen) then
!          print*,'Dimensional error in ReadMesh'
!          stop
!          if(lenloc >= 3 .and. lenloc <= 4) then ! triangle or quadrilateral
!             elem%type = lenloc
!             elem%flen = lenloc
!             allocate(elem%face(1:max_face, 1:elem%flen))
!             elem%face(idx,1:lenloc) = iloc(1:elem%flen)
!             elem%HGnode = .false.
!
!          elseif(lenloc >= 9) then  ! element with hanging nodes
!             elem%type = iloc(1)
!             elem%flen = iloc(2)
!             allocate(elem%face(1:max_face, 1:elem%flen))
!
!             elem%face(idx,1:iloc(2)) = iloc(3:2+elem%flen)
!
!             elem%HGnode = .true.
!             allocate(elem%HGvertex(1:iloc(1)) )
!             elem%HGvertex(1:iloc(1)) = iloc(3+elem%flen: 2+elem%flen+elem%type)
!
!             call SetHGelem(elem, grid%npoin, grid%x)
!          else
!             print*,'Bad data in init file ',gridfile,',  lenloc = ',lenloc
!             stop
!          endif
!       elseif(nbDim == 3) then   !3D case,  tetrahedra, pyramid or hexahedra
!          elem%type = lenloc
!          if (lenloc == 4) then ! tetrahedra
!             elem%flen = lenloc
!             allocate(elem%face(1:max_face, 1:elem%flen))
!             elem%face(nbnode,:) = 3
!
!             ! check orientation of tetrahedra
!             vec1(1:nbDim) = grid%x(iloc(2),1:nbDim) - grid%x(iloc(1),1:nbDim)
!             vec2(1:nbDim) = grid%x(iloc(3),1:nbDim) - grid%x(iloc(1),1:nbDim)
!             ! print*,'mmm ',grid%x(iloc(2),1:nbDim),grid%x(iloc(1),1:nbDim),vec1(1:nbDim)
!             vec3(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
!             vec3(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
!             vec3(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
!             vec1(1:nbDim) = grid%x(iloc(4),1:nbDim) - grid%x(iloc(1),1:nbDim)
!             if ((vec1(1)*vec3(1)+vec1(2)*vec3(2)+vec1(3)*vec3(3))<0.0) then
!                tmp=iloc(2)
!                iloc(2)=iloc(3)
!                iloc(3)=tmp
!                counter = counter +1
!             endif
!             !end check orientation of tetrahedra
!
!          elseif (lenloc == 5) then  ! pyramids
!             elem%flen = lenloc
!             allocate(elem%face(1:max_face, 1:elem%flen))
!             elem%face(nbnode,1:4) = 3
!             elem%face(nbnode,5) = 4
!          else                        ! hexahedra ??
!             elem%flen = lenloc-2
!             allocate(elem%face(1:max_face, 1:elem%flen))
!             elem%face(nbnode,:) = 4
!          end if
!          elem%face(idx,1:lenloc) = iloc(1:elem%flen)
!          elem%HGnode = .false.
!
!       endif
!
!       ! for red-green refinement
!       elem%RGtype = 'N'
!       elem%RGhistory = 0
!       elem%RGlevel = 0
!       elem%RGindex = 0
!    enddo
!
!    allocate(grid%b_edge(1:grid%nbelm) )
!    do i=1,grid%nbelm
!       read(ifile,*)  grid%b_edge(i)%lbn(1:2),grid%b_edge(i)%ibc
!    enddo
!    close(ifile)
!
!  end subroutine ReadMesh    endif
!
!       if(nbDim == 2) then
!          if(lenloc >= 3 .and. lenloc <= 4) then ! triangle or quadrilateral
!             elem%type = lenloc
!             elem%flen = lenloc
!             allocate(elem%face(1:max_face, 1:elem%flen))
!             elem%face(idx,1:lenloc) = iloc(1:elem%flen)
!             elem%HGnode = .false.
!
!          elseif(lenloc >= 9) then  ! element with hanging nodes
!             elem%type = iloc(1)
!             elem%flen = iloc(2)
!             allocate(elem%face(1:max_face, 1:elem%flen))
!
!             elem%face(idx,1:iloc(2)) = iloc(3:2+elem%flen)
!
!             elem%HGnode = .true.
!             allocate(elem%HGvertex(1:iloc(1)) )
!             elem%HGvertex(1:iloc(1)) = iloc(3+elem%flen: 2+elem%flen+elem%type)
!
!             call SetHGelem(elem, grid%npoin, grid%x)
!          else
!             print*,'Bad data in init file ',gridfile,',  lenloc = ',lenloc
!             stop
!          endif
!       elseif(nbDim == 3) then   !3D case,  tetrahedra, pyramid or hexahedra
!          elem%type = lenloc
!          if (lenloc == 4) then ! tetrahedra
!             elem%flen = lenloc
!             allocate(elem%face(1:max_face, 1:elem%flen))
!             elem%face(nbnode,:) = 3
!
!             ! check orientation of tetrahedra
!             vec1(1:nbDim) = grid%x(iloc(2),1:nbDim) - grid%x(iloc(1),1:nbDim)
!             vec2(1:nbDim) = grid%x(iloc(3),1:nbDim) - grid%x(iloc(1),1:nbDim)
!             ! print*,'mmm ',grid%x(iloc(2),1:nbDim),grid%x(iloc(1),1:nbDim),vec1(1:nbDim)
!             vec3(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
!             vec3(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
!             vec3(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
!             vec1(1:nbDim) = grid%x(iloc(4),1:nbDim) - grid%x(iloc(1),1:nbDim)
!             if ((vec1(1)*vec3(1)+vec1(2)*vec3(2)+vec1(3)*vec3(3))<0.0) then
!                tmp=iloc(2)
!                iloc(2)=iloc(3)
!                iloc(3)=tmp
!                counter = counter +1
!             endif
!             !end check orientation of tetrahedra
!
!          elseif (lenloc == 5) then  ! pyramids
!             elem%flen = lenloc
!             allocate(elem%face(1:max_face, 1:elem%flen))
!             elem%face(nbnode,1:4) = 3
!             elem%face(nbnode,5) = 4
!          else                        ! hexahedra ??
!             elem%flen = lenloc-2
!             allocate(elem%face(1:max_face, 1:elem%flen))
!             elem%face(nbnode,:) = 4
!          end if
!          elem%face(idx,1:lenloc) = iloc(1:elem%flen)
!          elem%HGnode = .false.
!
!       endif
!
!       ! for red-green refinement
!       elem%RGtype = 'N'
!       elem%RGhistory = 0
!       elem%RGlevel = 0
!       elem%RGindex = 0
!    enddo
!
!    allocate(grid%b_edge(1:grid%nbelm) )
!    do i=1,grid%nbelm
!       read(ifile,*)  grid%b_edge(i)%lbn(1:2),grid%b_edge(i)%ibc
!    enddo
!    close(ifile)
!
!  end subroutine ReadMesh


  !> read the mesh form gridfile in the "grid" file
  subroutine WriteMesh(gridfile, grid)
    character(len=*), intent(in) :: gridfile
    type(mesh), intent(in), target :: grid
    class(element), pointer :: elem
    integer :: ifile=12
    integer:: i, lenloc

    open(ifile, file=gridfile, status='UNKNOWN')
    write(ifile,*) grid%npoin, grid%nelem, grid%nbelm, grid%nbc
    write(ifile,*) grid%xper(1,1:nbDim),grid%iper(1,1:2),&
          grid%xper(2, 1:nbDim),grid%iper(2,1:2)

    do i=1,grid%npoin
       write(ifile,*) grid%x(i,1:nbDim)
    enddo

    do i=1, grid%nelem
       elem => grid%elem(i)

       if(elem%HGnode) then
          lenloc = 2 + elem%flen + elem%type
          write(ifile, *) lenloc, elem%type, elem%flen, elem%face(idx,1:elem%flen), &
               elem%HGvertex(1:elem%type)
       else
          write(ifile, *) elem%type, elem%face(idx,1:elem%flen)
       endif
    enddo

    do i=1,grid%nbelm
       write(ifile,*)  grid%b_edge(i)%lbn(1:2),grid%b_edge(i)%ibc
    enddo
    close(ifile)

  end subroutine WriteMesh

  !> write the file gridfile.mesh in form readable by Medit software
  subroutine WriteMeshMedit(grid, output_type )
    class(mesh), intent(in) :: grid
    character(len=*), intent(in) :: output_type
    character(len=50) :: mesh_name
    character(len=3) :: ch3
    integer :: text_size, is, num_size, i


    if(state%space%adapt%adapt_level > 0) then
        is = int(log(1.*state%space%adapt%adapt_level)/log(10.))
    else
        is = 0
    endif

    num_size = 3
    write( ch3, '(i3)' ) state%space%adapt%adapt_level

    if(output_type == 'T') then
         mesh_name = 'Tot-000    '
         text_size = 4
         mesh_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         mesh_name(num_size+text_size+2:num_size+text_size+6) = '.mesh'
    elseif (output_type == 'D') then
         mesh_name = 'Disc-000    '
         text_size = 5
         mesh_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         mesh_name(num_size+text_size+2:num_size+text_size+6) = '.mesh'
    elseif (output_type == 'A') then
         mesh_name = 'Alg-000    '
         text_size = 4
         mesh_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         mesh_name(num_size+text_size+2:num_size+text_size+6) = '.mesh'
    elseif (output_type == 'EstT') then
         mesh_name = 'EstTot-000    '
         text_size = 7
         mesh_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         mesh_name(num_size+text_size+2:num_size+text_size+6) = '.mesh'
    elseif (output_type == 'EstD') then
         mesh_name = 'EstDisc-000    '
         text_size = 8
         mesh_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         mesh_name(num_size+text_size+2:num_size+text_size+6) = '.mesh'
    elseif (output_type == 'EstA') then
         mesh_name = 'EstAlg-000    '
         text_size = 7
         mesh_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            mesh_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         mesh_name(num_size+text_size+2:num_size+text_size+6) = '.mesh'
    else
      Print*, 'Unknown output_type. Possible choices: T, D, A, EstT, EstD, EstA'
    endif

    !mesh_name = gridfile

    open(20+state%time%iter, file=mesh_name, status='UNKNOWN', position = 'append')
    write(20+state%time%iter,'(a22)') 'MeshVersionFormatted 1'
    write(20+state%time%iter,'(a1)') ' '
    write(20+state%time%iter,'(a11)') 'Dimension 2'
    write(20+state%time%iter,'(a1)') ' '
    write(20+state%time%iter,'(a8)') 'Vertices'
    write(20+state%time%iter,'(i6)') grid%npoin
    write(20+state%time%iter,'(a1)') ' '

    do i=1,grid%npoin
       write(20+state%time%iter,*) grid%x(i,1:nbDim), 0
    enddo

    write(20+state%time%iter,'(a1)') ' '
    write(20+state%time%iter,'(a9)') 'Triangles'
    write(20+state%time%iter,'(i6)') grid%nelem
    write(20+state%time%iter,'(a1)') ' '

    do i=1,grid%nelem

       if( grid%elem(i)%HGnode) then
         write(20+state%time%iter,*) grid%elem(i)%face(idx,grid%elem(i)%HGvertex(1)), &
                                grid%elem(i)%face(idx,grid%elem(i)%HGvertex(2)), &
                                grid%elem(i)%face(idx,grid%elem(i)%HGvertex(3)), 0   ! for triangles!!!
                                else
         write(20+state%time%iter,*) grid%elem(i)%face(idx,1), grid%elem(i)%face(idx,2), grid%elem(i)%face(idx,3), 0   ! for triangles!!!
       endif

    enddo

    write(20+state%time%iter,'(a1)') ' '
    write(20+state%time%iter,'(a3)') 'End'

    close(20+state%time%iter)

  end subroutine WriteMeshMedit


  !> write the files elements-***.txt and vertices-***.txt for visualization by Matlab
  subroutine WriteMeshMatlab( grid )
    class(mesh), intent(in) :: grid
    character(len=50) :: elem_name, vert_name
    character(len=3) :: ch3
    integer :: text_size, is, num_size, i


    if(state%space%adapt%adapt_level > 0) then
        is = int(log(1.*state%space%adapt%adapt_level)/log(10.))
    else
        is = 0
    endif

    num_size = 3
    write( ch3, '(i3)' ) state%space%adapt%adapt_level

    elem_name = 'elements-000   '
    text_size = 9
    elem_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
    if (stop_crit == 'L') then
       elem_name(num_size+text_size+1:num_size+text_size+2) = 'L'
    elseif (stop_crit == 'G') then
       elem_name(num_size+text_size+1:num_size+text_size+2) = 'G'
    elseif (stop_crit == 'N') then
       elem_name(num_size+text_size+1:num_size+text_size+2) = 'N'
    endif
    elem_name(num_size+text_size+2:num_size+text_size+5) = '.txt'



    vert_name = 'coordinates-000   '
    text_size = 12
    vert_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
    if (stop_crit == 'L') then
       vert_name(num_size+text_size+1:num_size+text_size+2) = 'L'
    elseif (stop_crit == 'G') then
       vert_name(num_size+text_size+1:num_size+text_size+2) = 'G'
    elseif (stop_crit == 'N') then
       vert_name(num_size+text_size+1:num_size+text_size+2) = 'N'
    endif
    vert_name(num_size+text_size+2:num_size+text_size+5) = '.txt'

    open(20+state%time%iter, file=elem_name, status='UNKNOWN', position = 'append')
    do i=1,grid%nelem

       if( grid%elem(i)%HGnode) then
         write(20+state%time%iter,*) grid%elem(i)%face(idx,grid%elem(i)%HGvertex(1)), &
                                grid%elem(i)%face(idx,grid%elem(i)%HGvertex(2)), &
                                grid%elem(i)%face(idx,grid%elem(i)%HGvertex(3))   ! for triangles!!!
                                else
         write(20+state%time%iter,*) grid%elem(i)%face(idx,1), grid%elem(i)%face(idx,2), grid%elem(i)%face(idx,3)   ! for triangles!!!
       endif

    enddo

    close(20+state%time%iter)


    open(20+state%time%iter, file=vert_name, status='UNKNOWN', position = 'append')
    do i=1,grid%npoin
       write(20+state%time%iter,*) grid%x(i,1:nbDim)
    enddo

    close(20+state%time%iter)


  end subroutine WriteMeshMatlab


  !> setting of hanging nodes
  subroutine SetHGelem(elem, npoin, x)
    type(element), intent(inout) :: elem
    integer, intent(in) :: npoin
    real, dimension(1:npoin,1:nbDim), intent(in) :: x
    integer :: i, j, k, i1, i2, k1, k2, l1, l2
    integer :: iface, idiff, ilevel
    real :: r0, r1, r2, rpos1, rpos2

    if(.not. elem%HGnode) then
       print*,'No hanging node'
       stop
    endif

    allocate(elem%HGface(1:2, 1:elem%flen) )

    iface = 0
    do i =1,elem%type
       i1 = elem%HGvertex(i)
       i2 = elem%HGvertex( mod(i,elem%type )+ 1)

       idiff = i2 - i1
       if(i == elem%type) idiff = idiff + elem%flen

       k1 = elem%face(idx, i1)
       k2 = elem%face(idx, i2)

       !write(*,'(a5,5i5,a1,2i5,a1,8i5)') &
       !     '&&&&&', elem%i,i,i1,i2,idiff,'|',k1,k2,'|', elem%HGvertex(:)

       if(idiff == 1) then ! face without hanging node
          iface = iface + 1
          elem%HGface(1, iface) = i
          elem%HGface(2, iface) = 1
       elseif(idiff >= 2 .and. idiff <= state%space%adapt%HG) then ! face with hanging nodes
          do j=1, idiff
             iface = iface + 1
             elem%HGface(1, iface) = i

             l1 = elem%face(idx, iface)
             l2 = elem%face(idx, mod(iface,elem%flen) + 1)

             if(dot_product(x(k1,1:nbDim)-x(k2,1:nbDim), x(l1,1:nbDim)-x(l2,1:nbDim)) <= 0.95* &
                  (dot_product(x(k1,1:nbDim)-x(k2,1:nbDim), x(k1,1:nbDim)-x(k2,1:nbDim) )* &
                  dot_product(x(l1,1:nbDim)-x(l2,1:nbDim), x(l1,1:nbDim)-x(l2,1:nbDim)))**0.5) then
                print*,'hanging node is not on a straight segment: u.v <= |u| |v|'
                write(*,'(a6,20i5)') 'elem =',elem%i,elem%flen,elem%face(idx,1:elem%flen)
                write(*,'(2es14.6,5i5)')x(k1,1:nbDim), k1,j,iface,idiff
                write(*,'(2es14.6,5i5)')x(l1,1:nbDim), l1
                write(*,'(2es14.6,5i5)')x(l2,1:nbDim), l2
                write(*,'(2es14.6,5i5)')x(k2,1:nbDim), k2
                print*,dot_product(x(k1,1:nbDim)-x(k2,1:nbDim), x(l1,1:nbDim)-x(l2,1:nbDim)), &
                     (dot_product(x(k1,1:nbDim)-x(k2,1:nbDim), x(k1,1:nbDim)-x(k2,1:nbDim) )* &
                  dot_product(x(l1,1:nbDim)-x(l2,1:nbDim), x(l1,1:nbDim)-x(l2,1:nbDim)))**0.5
                stop
             endif

             r0 = dot_product(x(k2,1:nbDim)-x(k1,1:nbDim),x(k2,1:nbDim)-x(k1,1:nbDim))**0.5
             r1 = dot_product(x(l1,1:nbDim)-x(k1,1:nbDim),x(l1,1:nbDim)-x(k1,1:nbDim))**0.5
             r2 = dot_product(x(l2,1:nbDim)-x(k1,1:nbDim),x(l2,1:nbDim)-x(k1,1:nbDim))**0.5
             rpos1 = r1/r0
             rpos2 = r2/r0

             !write(*,'(a6,4i5,a1,2es12.4)') 'diff',j,iface,l1,l2,'|', rpos1,rpos2

             ! seeking the level
             do k=1,state%space%adapt%max_HGlevel
                !write(*,'(a20,i5,4es12.4)') '...',k, rpos2-rpos1,  - 2.**(-k), &
                !     ((rpos2-rpos1) - 2.**(-k))/(2.**(-k))
                if( abs(((rpos2-rpos1) - 2.**(-k))/(2.**(-k))) <= 1E-3) then
                   ilevel = k
                   goto 10
                endif
             enddo
             print*,'Face in sharing face does not found (1)',elem%i,elem%flen
             write(*,*) x(elem%face(idx,elem%HGvertex(1)), 1:nbDim)
             write(*,*) x(elem%face(idx,elem%HGvertex(2)), 1:nbDim)
             write(*,*) x(elem%face(idx,elem%HGvertex(3)), 1:nbDim)


             write(*,'(a6,20i5)') '****', i, i1,i2, idiff, k1, k2
             write(*,'(a6,20i5)') 'elem:', elem%face(idx,:)
             write(*,'(a6,20i5)') 'HGel:', elem%HGface(1,:), elem%HGface(2,:)
             write(*,'(a6,20i5)') 'Hver:', elem%HGvertex(:)

             stop
10           continue

             !write(*,'(a6,4i5,a1,2es12.4,i5)') 'diff',j,iface,l1,l2,'|', rpos1,rpos2,ilevel

             !!elem%HGlevel = max(elem%HGlevel, ilevel)
             ! seeking the segment within level
             do k=0,2**ilevel-1
                if(abs(rpos1 - k* 2.**(-ilevel))/(2.**(-ilevel)) <= 1E-3) then
                   elem%HGface(2, iface) = 2.**ilevel + k
                   goto 20
                endif
             enddo
             print*,'Face in sharing face does not found (2)',elem%i,elem%flen
             write(*,*) x(elem%face(idx,elem%HGvertex(1)), 1:nbDim)
             write(*,*) x(elem%face(idx,elem%HGvertex(2)), 1:nbDim)
             write(*,*) x(elem%face(idx,elem%HGvertex(3)), 1:nbDim)


             write(*,'(a6,20i5)') '****', i, i1,i2, idiff, k1, k2
             write(*,'(a6,20i5)') 'elem:', elem%face(idx,:)
             write(*,'(a6,20i5)') 'HGel:', elem%HGface(1,:), elem%HGface(2,:)
             write(*,'(a6,20i5)') 'Hver:', elem%HGvertex(:)

             stop
20           continue


             !write(*,'(a3,4i3,6es10.2)') '###',i1, i2, k1,k2, x(k1,1:nbDim),x(k2,1:nbDim),rpos1,rpos2
             !write(*,'(a3,4i3,4es10.2, 2i3)') '###',j, iface, l1,l2, x(l1,1:nbDim),x(l2,1:nbDim),elem%HGface(1:2,iface)
             !write(*,*)


          enddo
       else
          print*,'Troubles in SetHGelem'
          stop
       endif

    enddo

  end subroutine SetHGelem

!  !> plot the mesh in the file 'meshfile' visualizable by gnuplot
!  subroutine PlotMesh(grid, meshfile)
!    type(mesh), intent(inout) :: grid
!    character(len=*), intent(in) :: meshfile
!    integer:: i,j, k, ifile = 11
!
!    open(ifile, file=meshfile)
!
!    !!print*,'@@@',grid%nelem
!
!    do i=1, grid%nelem
!       do j=0,grid%elem(i)%flen
!          k = grid%elem(i)%face(idx,mod(j,grid%elem(i)%flen) +1)
!          write(ifile,'(2es14.6,2i6)') grid%x(k,1:nbDim),0,i
!       enddo
!       write(ifile,'(x)')
!       write(ifile,'(x)')
!    enddo
!    close(ifile)
!
!  end subroutine PlotMesh

  ! seeking the corresponding "curved" nodes in  "type(curved_boundary) :: curvbound"
  ! between nodes p1 and p2
  subroutine SeekCurvedNodeLocal(xi, xcur1, xcur2,  rmin, node)
    real, dimension(1:nbDim), intent(in) :: xi      ! coordinates
    integer, dimension(1:2), intent(in) :: xcur1, xcur2   ! idx = 1 CBP, idx =2 pnt
    real, intent(in) :: rmin      ! minimal distance
    integer, intent(out) :: node
    integer :: k
    integer :: kstart, kend, klen
    real :: dist, dist1

    if(xcur1(1) /= xcur2(1)) then
       print*,'Different components in SeekCurvedNodeLocal'
       stop
    endif

    klen = size(curvbound%CBP(xcur1(1))%pnt, 1) -1

    kstart = xcur1(2)
    if(xcur1(2) < xcur2(2)) then
       kend = xcur2(2)
    elseif(xcur1(2) > xcur2(2)) then
       if(curvbound%closed(xcur1(1))) then
          kend = klen + xcur2(2)
       else
          print*,'TRouble in SeekCurvedNodeLocal', curvbound%closed(xcur1(1)), xcur1(1)
          stop
       endif
    else
       print*,'Nodes in SeekCurvedNodeLocal are identical'
       print*,' xcur1:', xcur1(1:2)
       print*,' xcur2:', xcur2(1:2)
       stop
    endif

    dist = 1E+50

    do k=kstart, kend
       dist1 = Distance(xi(1:nbDim), curvbound%CBP(xcur1(1))%pnt(mod(k-1,klen)+1,1:nbDim) )
       if(dist1 < dist ) then
          node = mod(k-1,klen)+1
          dist = dist1
       endif
    enddo

    if( dist > rmin) then ! node is not on curved boundary
       print*,'Warning, the founded node in SeekCurvedNodeLocal is far !!!'
    endif

    if(node == xcur1(2) .or. node == xcur2(2) )then
       print*,'Inner node is identical with the boundary one',&
            ' in SeekCurvedNodeLocal'
       write(*,'(a5,2i5,2e12.4)') 'BN1: ',xcur1(1:2), &
            curvbound%CBP(xcur1(1))%pnt(xcur1(2),1:nbDim)
       write(*,'(a5,2i5,2e12.4)') 'BN2: ',xcur2(1:2), &
            curvbound%CBP(xcur2(1))%pnt(xcur2(2),1:nbDim)
       write(*,'(a5,2i5,2e12.4)') 'IN : ',0,node, xi(1:nbDim)
       stop
    endif

  end subroutine SeekCurvedNodeLocal


  ! seeking the corresponding "curved" nodes in  "type(curved_boundary) :: curvbound"
  subroutine SeekCurvedNodeGlobal(xi, xcur, rmin)
    real, dimension(1:nbDim), intent(in) :: xi      ! coordinates
    integer, dimension(1:2), intent(out) :: xcur   ! idx = 1 CBP, idx =2 pnt
    real, intent(in) :: rmin      ! minimal distance
    integer :: i,k
    real :: dist, dist1

    xcur(1:2) = 0

    dist = 1E+50
    do i=1,curvbound%nbp    !  we seek its coresponding vertices from profiles
       do k=1, size(curvbound%CBP(i)%pnt, 1)
          dist1 = Distance(xi(1:nbDim), curvbound%CBP(i)%pnt(k,1:nbDim) )
          if(dist1 < dist ) then
             xcur(1) = i
             xcur(2) = k
             dist = dist1
          endif
       enddo
    enddo

    if( dist > rmin) then ! node is not on curved boundary
       xcur(1:2) = 0
    endif

    !print*,'#####',xcur(:), xi


  end subroutine SeekCurvedNodeGlobal

  !> seeking of the elements \f$ K\in {\cal T}_h\f$ adjacent to NON-polygonal
  !> boundary, using file 'prof_file'
  subroutine SeekCurvedBoundary(grid)
    type(mesh), intent(inout) :: grid
    !integer :: ifcur=13
    integer :: ip1, ip2, ib, i, j, k, k1, l
    real :: dist0, t1, t2
    real :: xc(10,2), delta
    integer :: jnode(10)
    integer :: counter_curved
    integer:: curved_deg

    call cpu_time(t1)

    curved_deg = grid%curved_deg-1

    counter_curved=0
    grid%xcur(1:grid%npoin,1:2) = -1

    !open(ifcur, file='curved-nodes', status='UNKNOWN')

    grid%b_edge(1:grid%nbelm)%icurv = 0

    grid%elem(1:grid%nelem)%ibcur = -1
    grid%elem(1:grid%nelem)%jcur = 0
    grid%elem(1:grid%nelem)%deg_cur = 1

    if(curved_deg >0 .and. curvbound%nbp > 0 ) then
       do ib=1,grid%nbelm   !  for each boundary element

          ip1 = grid%b_edge(ib)%lbn(1)
          ip2 = grid%b_edge(ib)%lbn(2)
          dist0 = Distance(grid%x(ip1, 1:nbDim),  grid%x(ip2, 1:nbDim))

          ! first node, not yet sought for previous b_edge?
          if(grid%xcur(ip1,1) == -1) &
               call SeekCurvedNodeGlobal(grid%x(ip1,:), grid%xcur(ip1,:), dist0 * 0.1)
          !if(ip1 == 8) print*,'!!!?',ip1,ip2,grid%xcur(ip1,:), dist0
          !if(ip2 == 8) print*,'..!?',ip1,ip2,grid%xcur(ip2,:), dist0

          ! second node
          call SeekCurvedNodeGlobal(grid%x(ip2,:), grid%xcur(ip2,:), dist0 * 0.1)

          ! segment on curved boundary
          if(grid%xcur(ip1,1) > 0 .and.  grid%xcur(ip2,1) >0 ) then
             if(grid%xcur(ip1,1) /= grid%xcur(ip2,1) )then
                print*,'Different boundary parts in SeekCurvedBoundary',ip1,ip2
                print*,grid%x(ip1,1:nbDim), grid%xcur(ip1,1:2)
                print*,grid%x(ip2,1:nbDim), grid%xcur(ip2,1:2)
                stop
             else
                !write(*,'(a3,32i5)')'!!!',ip1, ip2, ip(:),jp(:), &
               !grid%xcur(ip1,1:2), grid%xcur(ip2,1:2),&
               !ip(1) - grid%xcur(ip1,1), ip(2) - grid%xcur(ip2,1),&
               !jp(1) - grid%xcur(ip1,2), jp(2) - grid%xcur(ip2,2)


                if(curved_deg <= 2)then
                   ! maximal, Q_3 approximation, 2 inner nodes
                   do l=1,curved_deg
                      delta = 1.0 * l/(curved_deg + 1)
                      xc(l, 1:nbDim) = delta * grid%x(grid%b_edge(ib)%lbn(2), 1:nbDim) &
                           + (1.- delta) * grid%x(grid%b_edge(ib)%lbn(1), 1:nbDim)

                      call SeekCurvedNodeLocal(xc(l,:), grid%xcur(ip1,:), &
                           grid%xcur(ip2,:), 2*dist0, jnode(l))

                   enddo

                   !write(*,'(a6,10i6)') 'NEW',jnode(1:curved_deg)
                   grid%b_edge(ib)%icurv=curved_deg

                   ! number of inserted nodes = curved_deg
                   allocate(grid%b_edge(ib)%x_inn(1:curved_deg,1:nbDim) )

                   do l=1,curved_deg

                      !write(*,'(a6,6i6,20es12.4)') &
                      !     'cur:',l, ip1, ip2, grid%xcur(ip1,1), grid%xcur(ip1,2), jnode(l), &
                      !     grid%x(ip1,:), &
                      !     curvbound%CBP(grid%xcur(ip1,1))%pnt(grid%xcur(ip1,2),1:nbDim)


                      grid%b_edge(ib)%x_inn(l,1:nbDim) = &
                           curvbound%CBP(grid%xcur(ip1,1))%pnt(jnode(l),1:nbDim)

                      !write(ifcur,*) grid%b_edge(ib)%x_inn(l,1:nbDim), ib,l,jnode(l), &
                      !     grid%b_edge(ib)%lbn(1:nbDim)
                   enddo

                   counter_curved= counter_curved+1;

                   ! seeking of the curved edge of the corresponding element
                   i = grid%b_edge(ib)%itc
                   do j=1,grid%elem(i)%flen
                      k = grid%elem(i)%face(idx,j)
                      k1 = grid%elem(i)%face(idx,mod(j,grid%elem(i)%flen)+1)
                      if(k == grid%b_edge(ib)%lbn(1) .and. &
                           k1 == grid%b_edge(ib)%lbn(2) ) then
                         if(grid%elem(i)%ibcur > 0) then
                            print*,'Element ',i,' has probably more than 1 curved edge!!'
                            stop
                         endif

                         grid%elem(i)%ibcur = ib
                         grid%elem(i)%jcur = j
                         grid%elem(i)%deg_cur = grid%b_edge(ib)%icurv+1
                      endif
                   enddo

                   ! setting of the curved data for all nodes
                   !grid%xcur(grid%b_edge(ib)%lbn(1:2),1) = ip(1:nbDim)
                   !grid%xcur(grid%b_edge(ib)%lbn(1:2),2) = jp(1:nbDim)

                else
                   write(*, '(a23,i1,a26)') &
                        "Not yet implemented P_",curved_deg+1, &
                        " boundary approximation !!"
                   stop
                endif
             endif

          endif
       enddo
    endif
    !close(ifcur)

    !do i=1,grid%npoin
    !   if(grid%xcur(i,1) > 0) print*,')))))))))',i,  grid%xcur(i, 1:nbDim)
    !enddo


    grid%num_curv = counter_curved
    call cpu_time(t2)
    !if(state%space%adapt%max_adapt_level == 0 .or. state%space%adapt%adapt_level < 0 ) &
          write(*,'(a5,i1,a32,i5,a16)') &
          ' # P_',curved_deg+1,' approximation of the boundary, ', &
          grid%num_curv,' curved elements'

  end subroutine SeekCurvedBoundary

  !> construction of the spline
  subroutine SplineCurvedBoundary(grid )
    type(mesh), intent(inout) :: grid
    real, allocatable, dimension(:,:) :: a,b,c,d,e
    integer :: ib, l, l0, i, j, k, k1
    integer :: ifcur = 33
    integer :: counter_curved=0
    real:: t1, t2

    call cpu_time(t1)

    open(ifcur, file='curved-nodes', status='UNKNOWN')

    grid%b_edge(1:grid%nbelm)%icurv = 0

    ! pocet elementu
    do ib=1,grid%nbelm
       if(grid%b_edge(ib)%ibc == 3 .or. grid%b_edge(ib)%ibc == 8) then
          counter_curved = counter_curved+1
       endif
    enddo

    allocate(a(1:counter_curved,1:nbDim))
    allocate(b(1:counter_curved,1:nbDim))
    allocate(c(1:counter_curved,1:nbDim))
    allocate(d(1:counter_curved,1:nbDim))
    allocate(e(1:counter_curved,1:nbDim))

    ! ulozeni vrcholu do d
    l = 0
    do ib=1,grid%nbelm
       if(grid%b_edge(ib)%ibc == 3 .or. grid%b_edge(ib)%ibc == 8) then
          l = l+1
          d(l,1:nbDim) = grid%x(grid%b_edge(ib)%lbn(1), 1:nbDim)

       endif
    enddo

    ! usporadani vrcholu
    l0 = 1
    if (d(1,1) /=  1.) then
       do l = 1,counter_curved
          if (d(l,1) == 1.) exit
          e(l,1:nbDim) = d(l,1:nbDim)
       enddo
       l0 = l
       do l = l0,counter_curved
          d(l-l0+1,1:nbDim) = d(l,1:nbDim)
       enddo
       do l = 1,l0-1
          d(l-l0+1+counter_curved,1:nbDim) = e(l,1:nbDim)
       enddo
    endif

    ! konstrukce splinu
    ! sestaveni matice soustavy rovnic
    ! diagonala matice soustavy rovnic a
    a = 4.

    ! prava strana soustavy rovnic c
    do l = 2,counter_curved-1
       c(l,1:nbDim) = 6.*d(l+1,1:nbDim) - 12.*d(l,1:nbDim) + 6.*d(l-1,1:nbDim)
    enddo
    c(counter_curved,1:nbDim) = 6.*d(1,1:nbDim) - 12.*d(counter_curved,1:nbDim) + 6.*d(counter_curved-1,1:nbDim)
    ! okrajove podminky, 2. derivace rovny nule
    b(1,1:nbDim) = 0.

    ! reseni soustavy rovnic
    ! eliminace poddiagonaly matice soustavy rovnic
    do l = 3,counter_curved
       a(l,1:nbDim) = a(l,1:nbDim) - 1./a(l-1,1:nbDim)
       c(l,1:nbDim) = c(l,1:nbDim) - c(l-1,1:nbDim)/a(l-1,1:nbDim)
    enddo

    ! zpetna substituce do b
    b(counter_curved,1:nbDim) = c(counter_curved,1:nbDim)/a(counter_curved,1:nbDim)
    do l = counter_curved-1,2,-1
       b(l,1:nbDim) = (c(l,1:nbDim) - b(l+1,1:nbDim))/a(l,1:nbDim)
    enddo

    ! vypocet koeficientu a, c
    do l = 1,counter_curved-1
       a(l,1:nbDim) = (b(l+1,1:nbDim) - b(l,1:nbDim))/6.
       c(l,1:nbDim) = - b(l,1:nbDim)/3. - b(l+1,1:nbDim)/6. - d(l,1:nbDim) + d(l+1,1:nbDim)
    enddo
    a(counter_curved,1:nbDim) = - b(counter_curved,1:nbDim)/6.
    c(counter_curved,1:nbDim) = - b(counter_curved,1:nbDim)/3. - d(counter_curved,1:nbDim) +d(1,1:nbDim)
    ! vypocet koeficientu b
    b = b/2.

    ! puvodni poradi koeficientu
    if (l0 /= 1) then
       do l = 1,counter_curved-l0+1
          e(l,1:nbDim) = a(l,1:nbDim)
       enddo
       do l = 1,l0-1
          a(l,1:nbDim) = a(l-l0+1+counter_curved,1:nbDim)
       enddo
       do l = l0,counter_curved
          a(l,1:nbDim) = e(l-l0+1,1:nbDim)
       enddo
       do l = 1,counter_curved-l0+1
          e(l,1:nbDim) = b(l,1:nbDim)
       enddo
       do l = 1,l0-1
          b(l,1:nbDim) = b(l-l0+1+counter_curved,1:nbDim)
       enddo
       do l = l0,counter_curved
          b(l,1:nbDim) = e(l-l0+1,1:nbDim)
       enddo
       do l = 1,counter_curved-l0+1
          e(l,1:nbDim) = c(l,1:nbDim)
       enddo
       do l = 1,l0-1
          c(l,1:nbDim) = c(l-l0+1+counter_curved,1:nbDim)
       enddo
       do l = l0,counter_curved
          c(l,1:nbDim) = e(l-l0+1,1:nbDim)
       enddo
       do l = 1,counter_curved-l0+1
          e(l,1:nbDim) = d(l,1:nbDim)
       enddo
       do l = 1,l0-1
          d(l,1:nbDim) = d(l-l0+1+counter_curved,1:nbDim)
       enddo
       do l = l0,counter_curved
          d(l,1:nbDim) = e(l-l0+1,1:nbDim)
       enddo
    endif

    l = 0
    do ib=1,grid%nbelm
       if(grid%b_edge(ib)%ibc == 3 .or. grid%b_edge(ib)%ibc == 8) then
          l = l+1

          !! list of the nodes is here
          !print*,ib,grid%x(grid%b_edge(ib)%lbn(1), 1:nbDim ),'   !@'
          !print*,ib,grid%x(grid%b_edge(ib)%lbn(2), 1:nbDim ),'   !@'
          !print*,

          ! number of inserted nodes = 2
          grid%b_edge(ib)%icurv = 2

          allocate(grid%b_edge(ib)%x_inn(1:nbDim,1:nbDim) )

          ! setting of the inserted nodes is here
          grid%b_edge(ib)%x_inn(1,1:nbDim) = a(l,1:nbDim)/27. + b(l,1:nbDim)/9. + c(l,1:nbDim)/3. + d(l,1:nbDim)

          grid%b_edge(ib)%x_inn(2,1:nbDim) = a(l,1:nbDim)*(8./27.) + b(l,1:nbDim)*(4./9.) + c(l,1:nbDim)*(2./3.) + d(l,1:nbDim)

          do l0=1,2
             write(ifcur,*) grid%b_edge(ib)%x_inn(l0,1:nbDim), grid%b_edge(ib)%lbn(1:nbDim)
          enddo

          ! seeking of the curved edge of the corresponding element
          print*,'NOT YET CHECKED in SplineCurvedBoundary !!!!'
          i = grid%b_edge(ib)%itc
          do j=1,grid%elem(i)%flen
             k = grid%elem(i)%face(idx,j)
             k1 = grid%elem(i)%face(idx,mod(j,grid%elem(i)%flen)+1)
             if(k == grid%b_edge(ib)%lbn(1) .and. &
                  k1 == grid%b_edge(ib)%lbn(2) ) then
                if(grid%elem(i)%ibcur > 0) then
                   print*,'Element ',i,' has probably more than 1 curved edge!!'
                   stop
                endif

                grid%elem(i)%ibcur = ib
                grid%elem(i)%jcur = j
                grid%elem(i)%deg_cur = grid%b_edge(ib)%icurv+1
             endif
          enddo


       endif
    enddo

    deallocate(a,b,c,d,e)

    close(ifcur)

    grid%num_curv = counter_curved

    call cpu_time(t2)
    write(*,'(a47,i5,a16)') &
         ' # Cubic spline approximation of the boundary, ', &
         grid%num_curv,' curved elements'

  end subroutine SplineCurvedBoundary


!  !> seeking of neighbouring of elements \f$ K, K',\ K\cap K'\not=\emptyset\f$
!  subroutine SeekNeighbours(grid)
!    type(mesh), intent(inout) :: grid
!    integer, dimension(:,:), allocatable :: cyc
!    integer, dimension(:), allocatable :: len_cyc
!    integer, parameter :: max_cyc = 20
!    integer :: i, j, k, k1, l, ni, nj, is1, is2, ip1, ip2, ii, ie, iie, jj, ib, ibound
!
!    if (nbDim == 3) then
!       ! only tetrahedra
!       grid%n3elem = sum(grid%elem(:)%type-3)
!    else
!       grid%n3elem = sum(grid%elem(:)%type-2)
!    endif
!    grid%max_el_len = maxval( grid%elem(:)%type)
!
!    ! to each vertex we create a list of elements sharing this node
!    allocate(cyc(grid%npoin, max_cyc))
!    allocate(len_cyc(grid%npoin))
!
!    len_cyc(1:grid%npoin) = 0
!    do  i=1, grid%nelem
!       do j=1, grid%elem(i)%flen
!          k = grid%elem(i)%face(idx,j)
!
!          len_cyc(k) = len_cyc(k) + 1
!          cyc(k,len_cyc(k)) = i
!       enddo
!       grid%elem(i)%per = 0
!    enddo
!
!    ! inicialization  of arrays for neighbours
!
!    do i=1, grid%nelem
!       grid%elem(i)%face(neigh,:) = -10
!       grid%elem(i)%face(nei_i,:) = -10
!    enddo
!
!    do i=1, grid%nelem
!       ibound = 0
!       do j=1,grid%elem(i)%flen
!          if(grid%elem(i)%face(neigh,j) == -10) then
!             k = grid%elem(i)%face(idx,j)
!             k1 = grid%elem(i)%face(idx,mod(j,grid%elem(i)%flen)+1)
!             do l=1,len_cyc(k)
!                ni = cyc(k,l)
!                if(i /= ni) then
!                   do nj =1,grid%elem(ni)%flen
!                      if(k1 == grid%elem(ni)%face(idx,nj) .and.   &
!                           k== grid%elem(ni)%face(idx,mod(nj,grid%elem(ni)%flen)+1))&
!                           then
!                         grid%elem(i)%face(neigh,j) = ni
!                         grid%elem(ni)%face(neigh,nj) = i
!
!                         grid%elem(i)%face(nei_i,j) = nj
!                         grid%elem(ni)%face(nei_i,nj) = j
!
!                         goto 100
!                      endif
!                   enddo
!                endif
!             enddo
!100          continue
!          endif
!          ibound = ibound + 1
!       enddo
!       if(ibound > 0) then ! boundary element
!          allocate(grid%elem(i)%iBC(1:grid%elem(i)%flen) )
!          grid%elem(i)%iBC(:) = 0
!
!          allocate(grid%elem(i)%tBC(1:grid%elem(i)%flen) )
!          grid%elem(i)%tBC(:) = 0
!       endif
!    enddo
!
!    deallocate(cyc,len_cyc)
!
!    !do i= 2,6
!    !   write(*,'(15i5)') i, grid%elem(i)%face(idx,:), grid%elem(i)%face(neigh,:), &
!    !        grid%elem(i)%face(nei_i,:)
!    !enddo
!
!    ! seeking of element adjacent to boundary elements
!    grid%b_edge(1:grid%nbelm)%itc = -1
!    do i=1,grid%nelem
!       do j=1,grid%elem(i)%flen
!          if(grid%elem(i)%face(neigh,j) < 0 ) then
!             k = grid%elem(i)%face(idx,j)
!             k1 = grid%elem(i)%face(idx,mod(j,grid%elem(i)%flen)+1)
!             do ib=1,grid%nbelm
!                if(grid%b_edge(ib)%itc < 0) then
!                   !   print*,'?',k,k1,grid%b_edge(ib)%lbn(1:2)
!                   if(grid%b_edge(ib)%lbn(1) == k .and. &
!                        grid%b_edge(ib)%lbn(2) == k1 ) then
!                      grid%b_edge(ib)%itc = i
!                      grid%b_edge(ib)%jtc = j
!
!                      grid%elem(i)%face(neigh,j) = -ib
!                      goto 200
!                   endif
!                endif
!             enddo
!             print*,'Adjacent element to a boundary edge doesnt found'
!             print*,'elem:',i,j,grid%elem(i)%flen, k, k1
!             print*,'elem:',i,0,grid%elem(i)%face(idx,:)
!             !print*,'elem:',18,0,grid%elem(18)%face(idx,:)
!             print*, grid%x(k,1:nbDim)
!             print*, grid%x(k1,1:nbDim)
!             stop
!          endif
!200       continue
!          ! print*,'####','elem :',i,j,'nodes:',k,k1, '  iBe=',ib
!       enddo
!    enddo
!
!    do l=1, grid%periodic
!       is1 = 0
!       is2 = 0
!       ip1 = 0
!       ip2 = 0
!       do i=1,grid%nbelm
!          if(is1 == 0 .and. grid%b_edge(i)%ibc == grid%iper(l,1) ) is1 = i
!          if(grid%b_edge(i)%ibc == grid%iper(l,1) ) ip1 = ip1 + 1
!
!          if(is2 == 0 .and. grid%b_edge(i)%ibc == grid%iper(l,2) ) is2 = i
!          if(grid%b_edge(i)%ibc == grid%iper(l,2) ) ip2 = ip2 + 1
!
!       enddo
!       write(*,'(a4,7i5)') 'per', l,grid%iper(l,1:2), is1, is2, ip1, ip2
!
!       do k = 0, ip1 - 1
!          i = is1 + k
!          ii = is2 + ip2 - k - 1
!
!          !print*,'@@@', k,i, ii
!
!          ie  = grid%b_edge( i)%itc
!          iie = grid%b_edge(ii)%itc
!
!          j =  grid%b_edge( i)%jtc
!          jj = grid%b_edge(ii)%jtc
!
!          !write(*,'(a5,8i5)') '@@@', k,i, ii, ie, iie, j, jj
!
!          !write(101,*) (grid%x(grid%elem( ie)%face(idx, j), 1:nbDim) + &
!          !     grid%x(grid%elem( ie)%face(idx, mod(j,3)+1), 1:nbDim) )/ 2
!          !write(101,*) (grid%x(grid%elem(iie)%face(idx, mod(jj,3)+1), 1:nbDim) + &
!          !     grid%x(grid%elem(iie)%face(idx, jj), 1:nbDim) ) / 2.
!          !write(101,*)
!
!          grid%elem( ie)%per = iie
!          grid%elem(iie)%per = ie
!
!          grid%elem( ie)%face(neigh,  j ) = iie
!          grid%elem(iie)%face(neigh, jj ) = ie
!
!
!          grid%elem( ie)%face(nei_i,  j ) = jj
!          grid%elem(iie)%face(nei_i, jj ) = j
!       enddo
!
!    end do
!    !!print*,' End of periodicity'
!
!
!  end subroutine SeekNeighbours


  !> seeking of neighbouring of elements \f$ K, K',\ K\cap K'\not=\emptyset\f$
  subroutine SeekNeighbours3D(grid)
    type(mesh), intent(inout) :: grid
    integer, dimension(:,:), allocatable :: cyc
    integer, dimension(:), allocatable :: len_cyc
    integer, parameter :: max_cyc = 20
    integer :: i, j, k, l, ni, nj, is1, is2, ip1, ip2, ii, ie, iie, jj, ib, ibound
    integer, dimension(1:4) :: nodes, nodes1   ! nodes of face

    if (nbDim == 3) then
       ! only tetrahedra
       grid%n3elem = sum(grid%elem(:)%type-3)
    else
       grid%n3elem = sum(grid%elem(:)%type-2)
    endif
    grid%max_el_len = maxval( grid%elem(:)%type)

    ! to each vertex we create a list of elements sharing this node
    allocate(cyc(grid%npoin, max_cyc))
    allocate(len_cyc(grid%npoin))

    len_cyc(1:grid%npoin) = 0
    do  i=1, grid%nelem
       do j=1, grid%elem(i)%flen
          k = grid%elem(i)%face(idx,j)

          len_cyc(k) = len_cyc(k) + 1
          cyc(k,len_cyc(k)) = i
       enddo
    enddo

    ! inicialization  of arrays for neighbours

    do i=1, grid%nelem
       grid%elem(i)%face(neigh,:) = -10
    enddo

    do i=1, grid%nelem
       ibound = 0
       do j=1,grid%elem(i)%flen
          if(grid%elem(i)%face(neigh,j) == -10) then
             call GetFaceIndexes(grid%elem(i)%face,nodes,j,grid%elem(i)%flen)
             !k = grid%elem(i)%face(idx,j)
             !k1 = grid%elem(i)%face(idx,mod(j,grid%elem(i)%flen)+1)
             do l=1,len_cyc(nodes(1))
                ni = cyc(nodes(1),l)
                if(i /= ni) then
                   do nj =1,grid%elem(ni)%flen
                      call GetFaceIndexes(grid%elem(ni)%face,nodes1,nj,grid%elem(ni)%flen)
                      do k=1, grid%elem(i)%face(nbnode,j)
                         if(((grid%elem(i)%face(nbnode,j) == 2) .and. &
                              ((nodes(k) == nodes1(2)) .and. &
                              (nodes(mod(k,grid%elem(i)%face(nbnode,j))+1) == nodes1(1)))) &
                              .or. &
                              ((grid%elem(i)%face(nbnode,j) == 3) .and. &
                              ((nodes(k) == nodes1(3)) .and. &
                              (nodes(mod(k,grid%elem(i)%face(nbnode,j))+1) == nodes1(2))&
                              .and. &
                              (nodes(mod(k+1,grid%elem(i)%face(nbnode,j))+1) == nodes1(1))))&
                              .or. &
                              ((grid%elem(i)%face(nbnode,j) == 4) .and. &
                              ((nodes(k) == nodes1(4)) .and. &
                              (nodes(mod(k,grid%elem(i)%face(nbnode,j))+1) == nodes1(3)) .and. &
                              (nodes(mod(k+1,grid%elem(i)%face(nbnode,j))+1) == nodes1(2)) .and. &
                              (nodes(mod(k+2,grid%elem(i)%face(nbnode,j))+1) == nodes1(1))))) then

                            grid%elem(i)%face(neigh,j) = ni
                            grid%elem(ni)%face(neigh,nj) = i

                            grid%elem(i)%face(nei_i,j) = nj
                            grid%elem(ni)%face(nei_i,nj) = j

                            grid%elem(i)%face(rot,j) = k-1
                            grid%elem(ni)%face(rot,nj) = k-1

                            goto 100
                         endif
                      enddo
                   enddo
                endif
             enddo
100          continue
          endif
          ibound = ibound + 1
       enddo
       if(ibound > 0) then ! boundary element
          allocate(grid%elem(i)%iBC(1:grid%elem(i)%flen) )
          grid%elem(i)%iBC(:) = 0

          allocate(grid%elem(i)%tBC(1:grid%elem(i)%flen) )
          grid%elem(i)%tBC(:) = 0
       endif
    enddo

    deallocate(cyc,len_cyc)

    !do i= 2,6
    !   write(*,'(15i5)') i, grid%elem(i)%face(idx,:), grid%elem(i)%face(neigh,:), &
    !        grid%elem(i)%face(nei_i,:)
    !enddo

    ! seeking of element adjacent to boundary elements
    grid%b_edge(1:grid%nbelm)%itc = -1
    do i=1,grid%nelem
       do j=1,grid%elem(i)%flen
          if(grid%elem(i)%face(neigh,j) < 0 ) then
             call GetFaceIndexes(grid%elem(i)%face,nodes,j,grid%elem(i)%flen)
             !k = grid%elem(i)%face(idx,j)
             !k1 = grid%elem(i)%face(idx,mod(j,grid%elem(i)%flen)+1)
             do ib=1,grid%nbelm
                if(grid%b_edge(ib)%itc < 0) then
                   !   print*,'?',k,k1,grid%b_edge(ib)%lbn(1:2)
                   !                    print*,grid%elem(i)%face(nbnode,j)
                   do k=1, grid%elem(i)%face(nbnode,j)
                      if ((((grid%elem(i)%face(nbnode,j) == 2) .and. &
                           ((nodes(k) == grid%b_edge(ib)%lbn(2)) .and. &
                           (nodes(mod(k,grid%elem(i)%face(nbnode,j))+1) == grid%b_edge(ib)%lbn(1)))) .or. &
                           ((grid%elem(i)%face(nbnode,j) == 3) .and. &
                           ((nodes(k) == grid%b_edge(ib)%lbn(3)) .and. &
                           (nodes(mod(k,grid%elem(i)%face(nbnode,j))+1) == grid%b_edge(ib)%lbn(2)) .and. &
                           (nodes(mod(k+1,grid%elem(i)%face(nbnode,j))+1) == grid%b_edge(ib)%lbn(1)))) .or. &
                           ((grid%elem(i)%face(nbnode,j) == 4) .and. &
                           ((nodes(k) == grid%b_edge(ib)%lbn(4)) .and. &
                           (nodes(mod(k,grid%elem(i)%face(nbnode,j))+1) == grid%b_edge(ib)%lbn(3)) .and. &
                           (nodes(mod(k+1,grid%elem(i)%face(nbnode,j))+1) == grid%b_edge(ib)%lbn(2)) .and. &
                           (nodes(mod(k+2,grid%elem(i)%face(nbnode,j))+1) == grid%b_edge(ib)%lbn(1))))) .or. &
                           (((grid%elem(i)%face(nbnode,j) == 2) .and. &
                           ((nodes(k) == grid%b_edge(ib)%lbn(1)) .and. &
                           (nodes(mod(k,grid%elem(i)%face(nbnode,j))+1) == grid%b_edge(ib)%lbn(2)))) .or. &
                           ((grid%elem(i)%face(nbnode,j) == 3) .and. &
                           ((nodes(k) == grid%b_edge(ib)%lbn(1)) .and. &
                           (nodes(mod(k,grid%elem(i)%face(nbnode,j))+1) == grid%b_edge(ib)%lbn(2)) .and. &
                           (nodes(mod(k+1,grid%elem(i)%face(nbnode,j))+1) == grid%b_edge(ib)%lbn(3)))) .or. &
                           ((grid%elem(i)%face(nbnode,j) == 4) .and. &
                           ((nodes(k) == grid%b_edge(ib)%lbn(1)) .and. &
                           (nodes(mod(k,grid%elem(i)%face(nbnode,j))+1) == grid%b_edge(ib)%lbn(2)) .and. &
                           (nodes(mod(k+1,grid%elem(i)%face(nbnode,j))+1) == grid%b_edge(ib)%lbn(3)) .and. &
                           (nodes(mod(k+2,grid%elem(i)%face(nbnode,j))+1) == grid%b_edge(ib)%lbn(4)))))) then
                         !                    if(grid%b_edge(ib)%lbn(1) == k .and. &
                         !                         grid%b_edge(ib)%lbn(2) == k1 ) then
                         grid%b_edge(ib)%itc = i
                         grid%b_edge(ib)%jtc = j

                         grid%elem(i)%face(neigh,j) = -ib
                         goto 200
                      endif
                   enddo
                endif
             enddo

             print*,nodes(:),'Adjacent element to a boundary edge doesnt found'
             stop
          endif
200       continue
          ! print*,'####','elem :',i,j,'nodes:',k,k1, '  iBe=',ib
       enddo
    enddo

    do l=1, grid%periodic
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

          j =  grid%b_edge( i)%jtc
          jj = grid%b_edge(ii)%jtc

          !write(*,'(a5,8i5)') '@@@', k,i, ii, ie, iie, j, jj

          !write(101,*) (grid%x(grid%elem( ie)%face(idx, j), 1:nbDim) + &
          !     grid%x(grid%elem( ie)%face(idx, mod(j,3)+1), 1:nbDim) )/ 2
          !write(101,*) (grid%x(grid%elem(iie)%face(idx, mod(jj,3)+1), 1:nbDim) + &
          !     grid%x(grid%elem(iie)%face(idx, jj), 1:nbDim) ) / 2.
          !write(101,*)


          grid%elem( ie)%face(neigh,  j ) = iie
          grid%elem(iie)%face(neigh, jj ) = ie


          grid%elem( ie)%face(nei_i,  j ) = jj
          grid%elem(iie)%face(nei_i, jj ) = j
       enddo

    end do
    !!print*,' End of periodicity'


  end subroutine SeekNeighbours3D


!  !> computing of geometrical properties of each element \f$ K\in{\cal T}_h\f$,
!  !> the mapping \f$ F:\hat{K}\to K,\ K\in{\cal T}_h\f$
!  subroutine ComputeElementGeometry(grid, elem)
!    type(mesh), intent(inout) :: grid
!    type(element), intent(inout):: elem  ! elem = element
!    real :: diam, sum_inner_edges, rlen
!    real, allocatable, dimension (:,:) :: x, x0
!    !real, allocatable, dimension (:,:) :: xi, Fxi
!    integer:: flen, Fdeg, Fdof
!    integer:: j, k, k1, k0, ib, l, m, tmpint
!    integer:: l1, l2, l3, l4
!    real:: t1, t2
!    !    3D structures
!    integer, dimension(1:4) :: nodes   ! nodes of face
!    real, dimension(1:3) :: vec1, vec2
!
!
!    flen = elem%flen               ! number of faces of elem
!    Fdeg = elem%deg_cur           ! degree of curved element
!    if(Fdeg <= 1 ) Fdeg =  1      ! polygonal element Fdeg = 1
!
!    allocate(elem%F)
!    elem%F%deg = Fdeg
!
!    allocate(elem%xc(1:nbDim) )
!
!    diam = 0.
!    sum_inner_edges = 0.
!
!    ! computing of outer normals
!    allocate( elem%n(1:flen,1:nbDim))
!    allocate( elem%dn(1:flen))
!    do k=1,flen
!       if (nbDim==2) then
!          k1 = mod(k,flen) + 1
!          elem%n(k,1) = grid%x(elem%face(idx,k1),2) - grid%x(elem%face(idx,k),2)
!
!          elem%n(k,2) = grid%x(elem%face(idx,k),1)  - grid%x(elem%face(idx,k1),1)
!
!          elem%dn(k) = dot_product(elem%n(k,:),elem%n(k,:))**0.5
!       else
!          call GetFaceIndexes(elem%face,nodes,k,elem%flen)
!          if (elem%flen == 4) then  !tetrahedra
!             vec1(1:nbDim) = grid%x(nodes(2),1:nbDim) - grid%x(nodes(1),1:nbDim)
!             vec2(1:nbDim) = grid%x(nodes(3),1:nbDim) - grid%x(nodes(1),1:nbDim)
!             elem%n(k,1) = .5*(vec1(2)*vec2(3) - vec1(3)*vec2(2))
!             elem%n(k,2) = .5*(vec1(3)*vec2(1) - vec1(1)*vec2(3))
!             elem%n(k,3) = .5*(vec1(1)*vec2(2) - vec1(2)*vec2(1))
!          else
!             print*,'Only tetrahedra implemented in ComputeElementGeometry in 3d'
!             stop
!          endif
!
!          elem%dn(k) = dot_product(elem%n(k,:),elem%n(k,:))**0.5
!       endif
!    enddo
!
!    ! computing of mappings F: K^ -> K
!    if(elem%type == 3 .and. nbDim == 2) then  ! triangles
!       Fdof = (Fdeg + 1)*(Fdeg+2)/2
!       elem%F%dof = Fdof
!       ib = elem%ibcur     ! numbed of curved edge
!       j = elem%jcur       ! index of element edge adjacent to curved edge
!
!       !print*,'#####',elem%i,ib,j
!
!       if(j>0 .and. elem%HGnode) j = elem%HGface(1, j)
!
!       allocate(elem%F%F(1:Fdof,1:nbDim))
!       allocate(x(1:Fdof,1:nbDim))
!
!       if(elem%HGnode) then  ! element with hanging node
!          do k=1,3
!             x(k,1:nbDim) = grid%x(elem%face(idx,elem%HGvertex(k) ),1:nbDim)
!          enddo
!       else !element without hanging node
!          do k=1,3
!             x(k,1:nbDim) = grid%x(elem%face(idx,k),1:nbDim)
!          enddo
!       endif
!
!       elem%area = Area(x(1,:), x(2,:), x(3,:) )
!       ! evaluation of a "barycentre" of elem
!       do l=1,2
!          elem%xc(l) = sum(x(1:3,l) )/3
!       enddo
!
!       if(Fdeg == 2) then
!          ! piecewise quadratic approximation of the boundary, one node is inserted
!          ! computing of edge centres
!          do l=1,3
!             x(l+3, 1:nbDim) = (x(l,1:nbDim) + x(mod(l,3)+1,1:nbDim) )/2
!          enddo
!
!          ! a node on a curved edge is replaced by an inserted node
!
!          x(j+3,1:nbDim) = grid%b_edge(ib)%x_inn(1,1:nbDim)
!
!       elseif(Fdeg == 3) then
!
!          ! piecewise cubic approximation of the boundary
!          ! two nodes are inserted
!
!          ! computing of nodes on  edges
!          do l=1,3
!             x(2*l+3-1, 1:nbDim) = (2.*x(l,1:nbDim) + 1.*x(mod(l,3)+1,1:nbDim) )/3
!             x(2*l+3  , 1:nbDim) = (1.*x(l,1:nbDim) + 2.*x(mod(l,3)+1,1:nbDim) )/3
!          enddo
!          x(10, 1:nbDim) = (x(1,1:nbDim) + x(2,1:nbDim) + x(3,1:nbDim) )/3
!
!          ! two nodes on a curved edge is replaced by an inserted node
!          x(2*j+3-1,1:nbDim) = grid%b_edge(ib)%x_inn(1,1:nbDim)
!          x(2*j+3  ,1:nbDim) = grid%b_edge(ib)%x_inn(2,1:nbDim)
!
!          ! correction for central node, 1st possibility
!          !x(10,1) = sum(x(1:9,1))/9.
!          !x(10,2) = sum(x(1:9,2))/9.
!
!          ! correction for central node, 2nd possibility
!          l1 = 2*(j-1)+1+3
!          l2 = 2*(mod(j,3))+2+3
!          l3 = 2*(j-1)+2+3
!          l = mod(j,3) + 1
!          l4 = 2*(mod(l,3))+1+3
!          x(10,1:nbDim) = ( x(l1,1:nbDim) + x(l2,1:nbDim) + x(l3,1:nbDim) + x(l4,1:nbDim))/4.
!
!       elseif(Fdeg > 3) then
!          print*,'P_k, k> 3 approximation of boundary is not implemeneted'
!          stop
!       endif
!
!    elseif(elem%type == 4  .and. nbDim == 2) then
!       ! quadrilaterall elements
!       Fdof = (Fdeg+1)**2
!       elem%F%dof = Fdof
!       ib = elem%ibcur     ! numbed of curved edge
!       j = elem%jcur       ! index of element edge adjacent to curved edge
!       if(j>0 .and. elem%HGnode) j = elem%HGface(1, j)
!
!       allocate(elem%F%F(1:Fdof,1:nbDim))
!       allocate(x0(1:4,1:nbDim))
!       allocate(x(1:Fdof,1:nbDim))
!
!       if(elem%HGnode) then  ! element with hanging node
!          do k=1,4
!             x0(k,1:nbDim) = grid%x(elem%face(idx,elem%HGvertex(k) ),1:nbDim)
!          enddo
!       else !element without hanging node
!          do k=1,4
!             x0(k,1:nbDim) = grid%x(elem%face(idx,k),1:nbDim)
!          enddo
!       endif
!       elem%area = Area(x0(1,:), x0(2,:), x0(3,:) ) + Area(x0(3,:), x0(4,:), x0(1,:) )
!       ! evaluation of a "barycentre" of elem
!       do l=1,2
!          elem%xc(l) = sum(x0(1:4,l) )/4
!       enddo
!
!       do k0=1,Fdeg+1
!          do k1=1,Fdeg+1
!             k = (k1-1)*(Fdeg+1) + k0
!             t1 = 1.*(k0-1)/Fdeg
!             t2 = 1.*(k1-1)/Fdeg
!             !x(k,1:nbDim) = grid%x(elem%face(idx,1), 1:nbDim)  &
!             !     + (grid%x(elem%face(idx,2), 1:nbDim) - grid%x(elem%face(idx,1), 1:nbDim)) * t1  &
!             !     + (grid%x(elem%face(idx,4), 1:nbDim) - grid%x(elem%face(idx,1), 1:nbDim)) * t2  &
!             !     + (grid%x(elem%face(idx,3), 1:nbDim) - grid%x(elem%face(idx,2), 1:nbDim) &
!             !     -  grid%x(elem%face(idx,4), 1:nbDim) + grid%x(elem%face(idx,1), 1:nbDim)) * t1 * t2
!
!             x(k,1:nbDim) = x0(1, 1:nbDim)  &
!                  + (x0(2, 1:nbDim) - x0(1, 1:nbDim)) * t1  &
!                  + (x0(4, 1:nbDim) - x0(1, 1:nbDim)) * t2  &
!                  + (x0(3, 1:nbDim) - x0(2, 1:nbDim) -  x0(4, 1:nbDim) + x0(1, 1:nbDim)) * t1 * t2
!          enddo
!       enddo
!
!
!       if(Fdeg == 2) then
!          ! piecewise quadratic approximation of the boundary
!          ! one node is inserted
!
!          ! a node on a curved edge is replaced by an inserted node
!          if(j == 1) then
!             l = 2
!          elseif(j == 2) then
!             l = 6
!          elseif(j == 3) then
!             l = 8
!          elseif(j == 4) then
!             l = 4
!          else
!             print*,'#####Troubles in computing of mappings F !!'
!             stop
!          endif
!          x(l,1:nbDim) = grid%b_edge(ib)%x_inn(1,1:nbDim)
!
!          !a small modification ??
!          !x(5,1:nbDim) = (x(2,1:nbDim) + x(4,1:nbDim)+x(6,1:nbDim) + x(8,1:nbDim))/4.
!
!       elseif(Fdeg == 3) then
!          ! piecewise cubic approximation of the boundary
!          ! two nodes are inserted
!
!          ! a node on a curved edge is replaced by an inserted node
!          if(j == 1) then
!             l = 2
!             l1 = 3
!          elseif(j == 2) then
!             l = 8
!             l1 = 12
!          elseif(j == 3) then
!             l = 15
!             l1 = 14
!          elseif(j == 4) then
!             l = 9
!             l1 = 5
!          else
!             print*,'#####Troubles in computing of mappings F !!'
!             stop
!          endif
!          x(l,1:nbDim) = grid%b_edge(ib)%x_inn(1,1:nbDim)
!          x(l1,1:nbDim) = grid%b_edge(ib)%x_inn(2,1:nbDim)
!
!       elseif(Fdeg > 3) then
!          print*,'P_k, k> 3 approx. of boundary is not implemeneted'
!          stop
!       endif
!
!       deallocate( x0 )
!
!    elseif(elem%type == 4 .and. nbDim == 3) then
!       ! tetrahedra elements
!       Fdof = (Fdeg + 1)*(Fdeg+2)*(Fdeg+3)/6
!       elem%F%dof = Fdof
!       ib = elem%ibcur     ! numbed of curved edge
!       j = elem%jcur       ! index of element edge adjacent to curved edge
!
!       !print*,'#####',elem%i,ib,j
!
!       if(j>0 .and. elem%HGnode) j = elem%HGface(1, j)
!
!       allocate(elem%F%F(1:Fdof,1:nbDim))
!       allocate(x(1:Fdof,1:nbDim))
!
!       if(elem%HGnode) then  ! element with hanging node
!          do k=1,elem%flen
!             x(k,1:nbDim) = grid%x(elem%face(idx,elem%HGvertex(k) ),1:nbDim)
!          enddo
!       else !element without hanging node
!          do k=1,elem%flen
!             x(k,1:nbDim) = grid%x(elem%face(idx,k),1:nbDim)
!          enddo
!       endif
!       elem%area = Area3D(x(1,:), x(2,:), x(3,:), x(4,:) )
!       ! evaluation of a "barycentre" of elem
!       do l=1,nbDim
!          elem%xc(l) = sum(x(1:elem%flen,l))/elem%flen
!       enddo
!
!       if(Fdeg == 2) then
!          ! piecewise quadratic approximation of the boundary, one node is inserted on edge
!          ! computing of edge centres
!          tmpint=5
!          do k=1,3
!            do l=k+1,4
!              x(tmpint, 1:nbDim) = (x(l,1:nbDim) + x(k,1:nbDim) )/2
!              tmpint = tmpint + 1
!            enddo
!          enddo
!
!          ! a node on a curved edge is replaced by an inserted node
!!           x(j+3,1:2) = grid%b_edge(ib)%x_inn(1,1:2)
!
!       elseif(Fdeg == 3) then
!
!          ! piecewise cubic approximation of the boundary
!          ! two nodes are inserted to each edge
!
!          ! computing of nodes on  edges
!          tmpint=5
!          do k=1,3
!            do l=k+1,4
!              x(tmpint, 1:nbDim) = (2.*x(l,1:nbDim) + 1.*x(k,1:nbDim) )/3
!              x(tmpint+1, 1:nbDim) = (1.*x(l,1:nbDim) + 2.*x(k,1:nbDim) )/3
!              tmpint = tmpint + 2
!            enddo
!          enddo
!
!          do k=1,2
!            do l=k+1,3
!              do m=l+1,4
!                x(tmpint, 1:nbDim) = (x(l,1:nbDim) + x(k,1:nbDim)  + x(m,1:nbDim))/3
!                tmpint = tmpint + 1
!              enddo
!            enddo
!          enddo
!
!          ! two nodes on a curved edge is replaced by an inserted node
!!           x(2*j+3-1,1:2) = grid%b_edge(ib)%x_inn(1,1:2)
!!           x(2*j+3  ,1:2) = grid%b_edge(ib)%x_inn(2,1:2)
!
!          ! correction for central node, 1st possibility
!          !x(10,1) = sum(x(1:9,1))/9.
!          !x(10,2) = sum(x(1:9,2))/9.
!
!          ! correction for central node, 2nd possibility
!!           l1 = 2*(j-1)+1+3
!!           l2 = 2*(mod(j,3))+2+3
!!           l3 = 2*(j-1)+2+3
!!           l = mod(j,3) + 1
!!           l4 = 2*(mod(l,3))+1+3
!!           x(10,1:2) = ( x(l1,1:2) + x(l2,1:2) + x(l3,1:2) + x(l4,1:2))/4.
!
!       elseif(Fdeg > 3) then
!          print*,'P_k, k> 3 approximation of boundary is not implemeneted'
!          stop
!       endif
!    else
!       print*,'Sorry, only triangles and quadrilateralls are implemented'
!          stop
!    endif
!
!    call SetF(elem, Fdof, x(:,:) )
!
!    !write(*,*) 'F1:',elem%F%F(:,1)
!    !write(*,*) 'F2:',elem%F%F(:,2)
!!     write(*,*) 'F3:',elem%F%F(:,3)
!
!    ! computing the sum of length of all inner edges
!    diam = 0.
!    do k =1,elem%type-1
!       do l=k+1,elem%type
!          rlen= dot_product(x(l,1:nbDim) - x(k,1:nbDim),x(l,1:nbDim) - x(k,1:nbDim))**0.5
!          diam = max(diam, rlen)
!       enddo
!    enddo
!
!    do k=1,flen
!       rlen = VectorNorm(elem%n(k,:) )
!       diam = max(diam, rlen)
!       !if(elem%face(neigh,k) > 0) ! sum of all edges (innner and boundary)
!       sum_inner_edges = sum_inner_edges + rlen
!    enddo
!    elem%diam = diam
!    state%space%h = max(state%space%h, diam)
!
!    if(sum_inner_edges <= 0.) sum_inner_edges = elem%diam  ! only for singular case
!    elem%r_ins = 2*elem%area/sum_inner_edges
!
!    ! parameter for limiting
!    !elem%limit_par = 1./(elem%area**0.75)/elem%diam
!    !elem%limit_par = 1./(elem%area**0.75)/sum_inner_edges
!    if (nbDim == 3) then
!       elem%limit_par = 1./(elem%area)**(1./6.)/sum_inner_edges
!    else
!       elem%limit_par = 1./elem%area/sum_inner_edges
!    end if
!
!    deallocate( x )
!
!  end subroutine ComputeElementGeometry

!  !> computing of geometry \f$ \forall K\in{\cal T}_h\f$
!  !> after refinement, only refined elements have to be recomputed
!  subroutine ComputeGeometry(grid)
!    type(mesh), intent(inout) :: grid
!    integer:: i
!
!    !open(45,file='nodes')
!
!    state%space%h = 0.
!    state%space%domain_volume = 0.
!    do i=1,grid%nelem
!       call ComputeElementGeometry(grid, grid%elem(i))
!
!       state%space%domain_volume = state%space%domain_volume + grid%elem(i)%area
!
!       !if(.not. grid%elem(i)%F%iFlin) &
!       !if(grid%elem(i)%HGnode) &
!       !     call CheckElement(grid%elem(i), 45)
!    enddo
!    grid%h = state%space%h
!!    close(45)
!
!!    print*,'# Geometry computed'
!
!
!    !do i=1,grid%nelem
!    !      do j=1,grid%elem(i)%flen
!    !         if(grid%elem(i)%face(neigh,j) > 0) then
!    !            ii = grid%elem(i)%face(neigh,j)
!    !            write(98,*) grid%elem(i)%xc(1:nbDim)
!    !            write(98,*) grid%elem(ii)%xc(1:nbDim)
!    !            write(98,*)
!    !            write(98,*)
!    !         else
!    !            write(97,*) grid%elem(i)%xc(1:nbDim)
!!
!!             endif
!!
!!          enddo
!!       enddo
!!
!!       do k=1,grid%nbelm
!!          i = grid%b_edge(k)%itc
!!          j = grid%b_edge(k)%jtc
!!          if(grid%elem(i)%face(neigh,j) > 0) then
!!             ii = grid%elem(i)%face(neigh,j)
!!             write(95,*) grid%elem(i)%xc(1:nbDim)
!!             write(95,*) grid%elem(ii)%xc(1:nbDim)
!!             write(95,*)
!!             write(95,*)
!!          else
!!             write(94,*) grid%elem(i)%xc(1:nbDim)
!!
!!          endif
!!
!!       enddo
!!
!!       print*,'# Geometry computed, test done'
!!    stop
!
!
!
!  end subroutine ComputeGeometry


  subroutine ReadCurvedBoundary(profiles)
    character(len=*), intent(in) :: profiles
    integer :: prfl  = 51
    integer:: i, j, k

    open(unit = prfl, file = profiles,   status = 'old', action='read')

    read(prfl,*) curvbound%nbp
    if (curvbound%nbp > 0) then
       allocate(curvbound%CBP(1:curvbound%nbp))
       allocate(curvbound%closed(1:curvbound%nbp) )
       do i = 1, curvbound%nbp
          read(prfl,*) k
          allocate(curvbound%CBP(i)%pnt(1:k,1:2))
          do j = 1,k
             read(prfl,*) curvbound%CBP(i)%pnt(j,:)
          end do

          curvbound%closed(i) = .false.

          !print*,'!!!', curvbound%CBP(i)%pnt(1,:)
          !print*,'!!!', curvbound%CBP(i)%pnt(k,:)

          !print*,'@@@',Distance(curvbound%CBP(i)%pnt(1,:), curvbound%CBP(i)%pnt(k,:) )/  &
          !     (VectorNorm(curvbound%CBP(i)%pnt(1,:)) +  &
          !     VectorNorm(curvbound%CBP(i)%pnt(k,:) ))

          if(Distance(curvbound%CBP(i)%pnt(1,:), curvbound%CBP(i)%pnt(k,:) )/  &
               (VectorNorm(curvbound%CBP(i)%pnt(1,:)) +  &
               VectorNorm(curvbound%CBP(i)%pnt(k,:) )) <= 1E-05   ) then
             curvbound%closed(i) = .true.
          endif

       end do
    end if

    close(prfl)
    !print*,'# Reading of curved nodes from file "', profiles
  end subroutine ReadCurvedBoundary

  !> seek he list of elements having a common vertex with elem
  subroutine SeekElemSupports(gridN)
    type(mesh), intent(inout), target  :: gridN
    class(element), pointer :: elem
    integer, dimension(:), allocatable :: V  ! number  of elements corresponding to a vertex
    integer, dimension(:), allocatable :: E  ! number  of elements corresponding to an element
    integer, dimension(:,:), allocatable :: suppV  ! list of elements corresponding to a vertex
    integer, dimension(:,:), allocatable :: suppE  ! list of elements corresponding to an element
    integer :: i, j, k, l, m, ie
    integer :: maxdeg = 20

    !!print*,'###,   SeekElemSupports(gridN) called',gridN%ElemSupports

    allocate( V(1:gridN%npoin), suppV(1:gridN%npoin, 1:maxdeg) )
    V(:) = 0

    ! seeking of elements having a common vertex
    do i=1,gridN%nelem
       elem => gridN%elem(i)
       do j=1,elem%flen
          k = elem%face(idx, j)
          V(k) = V(k) + 1
          if(V(k) > maxdeg ) then
             write(*,*) 'Dimension full in SeekElemSupports'
             stop
          else
             suppV(k, V(k)) = i
          endif

       enddo
    enddo

    ! graphical verification
    !do i=1,gridN%npoin
    !!   write(*,'(a4,30i5)' ) '###',i,V(i), suppV(i,1:V(i))
    !   write(2000+elem%i, *) gridN%x(i,:)!
    !
    !   do j=1, V(i)
    !      elem => gridN%elem(suppV(i,j))
    !
    !      write(2000+i, *) gridN%x(i,:)
    !      write(2000+i, *) elem%xc(1:2)
    !      write(2000+i, '(x)' )
    !   enddo
    !enddo

    !print*,'#################################'

    allocate( E(1:gridN%nelem), suppE(1:gridN%nelem, 1: 2*maxdeg) )
    E(:) = 0
    ! creating of list of elemnts
    do i=1,gridN%nelem
       elem => gridN%elem(i)
       do j=1,elem%flen
          k = elem%face(idx, j)

          !write(*,'(a4,30i5)' ) '???',j,k,V(k), suppV(k,1:V(k))

          do l=1, V(k)
             ie = suppV(k, l) ! index of element in the list

             if(ie == i) goto 100 ! ie is identical with i
             do m = 1, E(i)   ! checking if ie is already inserted
                if(suppE(i, m) == ie) goto 100
             enddo
             ! ie is not inserted, we include it
             E(i) = E(i) + 1
             if(E(i) > 2*maxdeg ) then
                write(*,*) 'Dimension full in SeekElemSupports at 2'
                stop
             else
                suppE(i, E(i) ) = ie
                !write(*,'(a6,20i5)') '!!!',i,k,E(i), suppE(i, E(i) )
             endif

100          continue
          enddo
       enddo
    enddo

    ! filling of arrays elem%supp
    do i=1,gridN%nelem
       elem => gridN%elem(i)
       elem%isupp = E(i)
       allocate(elem%supp(E(i), 1:2) )

       !elem%supp(0, 1:1) = i
       elem%supp(1:E(i), 1) = suppE(i, 1:E(i) )

       ! elements sharing  a face       => elem%supp(j, 2) = 1
       ! elements sharing only a vertex => elem%supp(j, 2) = 2
       elem%supp(1:E(i), 2) = 2
       do j=1, E(i)
          if(elem%supp(j, 1) == elem%face(neigh, 1)  &
               .or. elem%supp(j, 1) == elem%face(neigh, 2)  &
               .or. elem%supp(j, 1) == elem%face(neigh, 3) ) elem%supp(j,2) = 1
       enddo

    enddo


    deallocate(V, E, suppV, suppE)

    gridN%ElemSupports  = .true.

    ! graphical verification
    !do i=1,gridN%nelem
    !   elem => gridN%elem(i)
    !   write(*,'(a4,30i5)' ) '###',i,elem%isupp,elem%supp(:,1)
    !   write(1000+elem%i, *) elem%xc(1:2)
    !   do j=1, elem%isupp
    !      write(1000+elem%i, *) elem%xc(1:2)
    !      write(1000+elem%i, *) gridN%elem(elem%supp(j,1) )%xc(1:2)
    !      write(1000+elem%i, '(x)' )
    !   enddo
    !enddo
    !
    !stop

  end subroutine SeekElemSupports




  !> seek he list of elements having a common vertex
  !> abs(V(i) ) is the number of elements in the list, for
  !> boundary vertexes V(i) < 0
  !> the list is oriented counter clock-wise
  subroutine SeekVertexSupports(gridN, maxdeg, inner, V, suppV)
    class(mesh), intent(inout), target  :: gridN
    class(element), pointer :: elem, elem1
    integer, intent(in) :: maxdeg
    logical, dimension(1:gridN%npoin) :: inner  ! inner nodes
    integer, dimension(1:gridN%npoin) :: V  ! number of elements corresponding to a vertex
    integer, dimension(1:gridN%npoin,1:maxdeg, 1:2) :: suppV  ! list of elements corresponding to a vertex
    integer :: i, j, k, l, ip, ie, ii, jp, jp1, itc
    integer, dimension(:), allocatable :: Lsupp

    V(:) = 0

    ! seeking of elements having a common vertex
    do i=1,gridN%nelem
       elem => gridN%elem(i)

       do j=1,elem%flen
          k = elem%face(idx, j)
          V(k) = V(k) + 1
          if(V(k) > maxdeg ) then
             write(*,*) 'Dimension full in SeekVertexSupports'
             stop
          else
             suppV(k, V(k),1) = i
             suppV(k, V(k),2) = j
          endif

       enddo
    enddo

    allocate(Lsupp(1:2) )

    inner(:) = .true.
    ! we distinguish the boundary and interior nodes
    do i=1,gridN%nbelm
       j =  gridN%b_edge(i)%lbn(1)
       k =  gridN%b_edge(i)%lbn(2)
       inner(j) = .false. ! some elements are carried out two times, doesn't mind
       inner(k) = .false.

       ! reorder the orientation counter-clock-wise
       ! first we set the first triangle
       ip = j
       itc = gridN%b_edge(i)%itc ! must be the first

       do ii=1,V(ip)
          ie =  suppV(ip, ii, 1)
          if(ie == itc) then
             Lsupp(1:2) =  suppV(ip, 1, 1:2)
             suppV(ip,  1, 1:2) = suppV(ip, ii, 1:2)
             suppV(ip, ii, 1:2) =  Lsupp(1:2)
             goto 50
          endif
       enddo
50     continue

       ! we redorder the other elements
       do ii=1,V(ip) - 1
          ie =  suppV(ip, ii, 1) ! actual element
          elem => gridN%elem(ie)

          j = suppV(ip, ii, 2)    ! inner index of vertex ip
          jp1 = mod(j,3) + 1
          jp = mod(jp1,3) + 1

          do l=ii+1, V(ip)
             k = suppV(ip, l,1)

             if(elem%face(neigh, jp) == k) then
                !  print*,'triangle found, switch', l, ii + 1
                Lsupp(1:2) =  suppV(ip, l, 1:2)
                suppV(ip, l, 1:2) = suppV(ip, ii+1, 1:2)
                suppV(ip, ii+1, 1:2) =  Lsupp(1:2)
                goto 10
             endif
          enddo
10        continue
       enddo

    enddo


    ! correct orientation
     do ip=1,gridN%npoin

        if(inner(ip) ) then ! inner node
           do ii=1,V(ip) - 1
              ie =  suppV(ip, ii, 1) ! actual element

              elem => gridN%elem(ie)
              j = suppV(ip, ii, 2)    ! inner index of vertex ip
              jp1 = mod(j,3) + 1
              jp = mod(jp1,3) + 1

              do l=ii+1, V(ip)
                 k = suppV(ip, l,1)

                 if(elem%face(neigh, jp) == k) then
                    !print*,'triangle found, switch', l, ii + 1
                    Lsupp(1:2) =  suppV(ip, l, 1:2)
                    suppV(ip, l, 1:2) = suppV(ip, ii+1, 1:2)
                    suppV(ip, ii+1, 1:2) =  Lsupp(1:2)
                    goto 20
                 endif
              enddo
20            continue

              !write(*,'(a10,20i5)') '########',suppV(ip, 1:V(ip), 1)
              !write(*,'(a10,20i5)') '########',suppV(ip, 1:V(ip), 2)
              !print*,'------------------------'
           enddo

        else  ! boundary nodes

        endif

        ! graphical verification
        ! if(mod(ip,5) == 3) then
        !    write(88,'(10es12.4)') gridN%x(ip,:)
        !    do ii=1,V(ip)
        !       ie =  suppV(ip, ii, 1) ! actual element
        !       elem => gridN%elem(ie)
        !       write(88,'(10es12.4)') elem%xc(:)
        !    enddo
        !    write(88,'(10es12.4)') gridN%x(ip,:)
        !    write(88,'(x)')
        ! endif

     enddo

     deallocate(Lsupp)

  end subroutine SeekVertexSupports


  !> split the computational domain to a cartesian grid and associate to each
  !> square cell the list of triangles whose barycenter is contained
  subroutine FindTriangleIndexs(gridO, nsq, rmx, max_tri, ntri, itri)
    type(mesh), intent(in), target	::  gridO
    class(element), pointer :: elem
    integer, intent(out) :: nsq ! domain will be split into nsq x nsq square cells
    integer, intent(out) :: max_tri  ! usually nelem
    real, dimension(1:2, 1:2), intent(out) :: rmx  ! frame of the computational domain
    !integer, dimension(1:nsq, 1:nsq, 1: max_tri), allocatable :: itri
    !integer, dimension(1:nsq, 1:nsq), allocatable :: ntri
    integer, dimension(:, :, :), allocatable, intent(out) :: itri
    integer, dimension(:, :), allocatable, intent(out) :: ntri
    integer :: ie, i, j, k

    ! spliting of the domain onto nsq x nsq square cells
!    print*, 'FR: real-> integer added in FindTriangleIndexs OK?', nsq
    nsq = int((gridO%nelem/2)**0.5)
    !!!!nsq = 15

    !max_tri = 250
    max_tri = gridO%nelem

    allocate( ntri(1:nsq, 1:nsq), itri(1:nsq, 1:nsq, 1: max_tri) )

    ! frame of the computational domain
    rmx(1, 1) = minval(gridO%x(:, 1) )
    rmx(1, 2) = maxval(gridO%x(:, 1) )
    rmx(2, 1) = minval(gridO%x(:, 2) )
    rmx(2, 2) = maxval(gridO%x(:, 2) )

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

  end subroutine FindTriangleIndexs



  !> find the traingle from gridO containing the node xp(1:2) (return the link elemO)
  !> and its barycentric coordinates x0 with respect to elemO
  subroutine  FindTriangleCoords(gridO, elemO, xp, x0, nsq, max_tri, ntri, itri, rmx  )
    class(mesh), intent(in), target :: gridO
    class(element), pointer, intent(inout) :: elemO
    class(element), pointer :: elemO_dist
    real, dimension(1:2), intent(in) :: xp
    real, dimension(1:2), intent(out) :: x0
    integer, intent(in) :: nsq, max_tri
    integer, dimension(1:nsq, 1:nsq), intent(in):: ntri
    integer, dimension(1:nsq, 1:nsq, 1:max_tri), intent(in):: itri
    real, dimension(1:2, 1:2), intent(in):: rmx
    real, dimension(:,:), allocatable :: xt !
    real, dimension(:), allocatable :: x0_bak
    integer :: si, sj, ssi, ssj, ss_plus, ii, jj, l
    integer :: itest
    real:: xtest(1:2), dist, dist_min

    allocate(x0_bak(1:2))
    dist_min = 1E+20

    !itest = 0
    !xtest(1) = -0.94444445766666663
    !xtest(2) = 0.22222219900000004
    !if(abs(xp(1) - xtest(1) ) < 1E-4 .and. abs(xp(2) - xtest(2)) <1E-4) itest = 1
    !if(itest == 1) print*,'#####',xp(1:2)

    allocate(xt(1:3, 1:2) ) ! coordinates of vertices of elem

    !write(98,*) xi(ip, 1:2), elem%i, ip
    ! starting relative coordinates of
    ssi = int(nsq * (xp(1) - rmx(1, 1) ) /  (rmx(1,2) - rmx(1, 1) )) + 1
    ssj = int(nsq * (xp(2) - rmx(2, 1) ) /  (rmx(2,2) - rmx(2, 1) )) + 1


    do ss_plus = 0, nsq
       do ii = -ss_plus, ss_plus
          do jj = -ss_plus, ss_plus
             ! checking of the correct square cell
             if(max(abs(ii), abs(jj)) == ss_plus) then
                si = ssi + ii
                sj = ssj + jj

                if(si >= 1 .and. si <= nsq .and. sj >= 1 .and. sj <= nsq) then
                   ! we go through the list of elements

                   do l=1, ntri(si, sj)
                      elemO => gridO%elem( itri(si, sj, l) )

                      if(elemO%HGnode) then
                         xt(1, 1:2) = gridO%x(elemO%face(idx, elemO%HGvertex(1)), 1:2)
                         xt(2, 1:2) = gridO%x(elemO%face(idx, elemO%HGvertex(2)), 1:2)
                         xt(3, 1:2) = gridO%x(elemO%face(idx, elemO%HGvertex(3)), 1:2)
                      else
                         xt(1, 1:2) = gridO%x(elemO%face(idx, 1), 1:2)
                         xt(2, 1:2) = gridO%x(elemO%face(idx, 2), 1:2)
                         xt(3, 1:2) = gridO%x(elemO%face(idx, 3), 1:2)
                      endif

                      call BarycCoordOne(xt(1:3, 1:2), xp(1:2), x0(1:2))

                      !if(itest > 0) then
                      !   write(*,'(a4,4i4,5es12.4)')'...',ss_plus, si, sj,  itri(si, sj, l), &
                      !        xp(1:2),  x0(1:2), 1. - x0(1) - x0(2)
                      !endif

                      ! the element found ??
                      if(x0(1) >= 0. .and. x0(2) >= 0. .and. x0(1) + x0(2) <= 1.) then
                         !if(itest > 0) then
                         !   print*,xt(1, 1:2), elemO%i
                         !   print*,xt(2, 1:2)
                         !   print*,xt(3, 1:2)
                         !endif
                         goto 10
                      else
                         dist = 0.
                         if(x0(1) < 0.) dist = dist + x0(1)**2
                         if(x0(2) < 0.) dist = dist + x0(2)**2
                         if(x0(1) + x0(2) > 1.) dist = dist + (x0(1) + x0(2) - 1.)**2
                         if(dist < dist_min) then
                            dist_min = dist
                            x0_bak(1:2) = x0(1:2)
                            elemO_dist => elemO
                         endif

                      endif
                   enddo ! l=1, ntri(si, sj)
                endif
             endif
          enddo ! jj
       enddo ! ii
    enddo ! ss_plus = 0, nsq
    !write(*,'(a15, 2es12.4,a55,es12.4)') &
    !     '### NOde xp=',xp(1:2),'  does not found, we insert the closets with dist = ',  dist_min
    !write(181,*) xp(1:2)
    !print*,'###',xp(:)
    !print*,'###',x0_bak(:), dist_min
    !print*,'###',elemO_dist%xc(:)
    !stop
    elemO => elemO_dist
    x0(1:2) = x0_bak(1:2)
10  continue

    deallocate(xt, x0_bak)

  end subroutine FindTriangleCoords


  !> xp  real physical cordinates of triangle elem
  !> Fx  real physical cordinates of ONE node (inside or outside) of  elem
  !> xi barycentric coordinates of Fx  with respect elem
  subroutine BarycCoordOne(xp, Fx,  xi )
    real, dimension(1:3, 1:2), intent(in) :: xp
    real, dimension(1:2), intent(in) :: Fx
    real, dimension(1:2), intent(inout) :: xi
    real :: D, Dx, Dy
    integer :: i

    D = (xp(2, 1) - xp(1, 1)) * (xp(3, 2) - xp(1, 2)) - (xp(2, 2) - xp(1, 2)) * (xp(3, 1) - xp(1, 1))
    ! cyclus of nodes
    !do i=1, Qdof1
    Dx = (Fx(   1) - xp(1, 1))*(xp(3, 2) - xp(1, 2)) - (Fx(   2) - xp(1, 2))*(xp(3, 1) - xp(1, 1))
    Dy = (xp(2, 1) - xp(1, 1))*(Fx(   2) - xp(1, 2)) - (xp(2, 2) - xp(1, 2))*(Fx(   1) - xp(1, 1))

    xi( 1) = Dx / D
    xi( 2) = Dy / D
    !enddo

  end subroutine BarycCoordOne

end module mesh_oper



