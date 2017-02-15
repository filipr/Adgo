!> mesh operations: reading, plotting, computing geometry
module mesh_oper3D
  use matrix_oper
  use main_data, g_grid => grid !mesh operation IS NOT USED for global grid (adaptation)

  public:: PlotMesh3D
  public:: GetFaceIndexes

contains

  !> plot the mesh in the file 'mesh' visualizable by Techplot
  subroutine PlotMesh3D(grid)
    type(mesh), intent(inout) :: grid
    integer:: i, ifile = 11

    open(ifile, file='mesh3D.plt')

    write(ifile,*) 'TITLE      = Adgfem Mesh Data'
    write(ifile,*) 'VARIABLES  =  "x0"  "x1"  "x2"'
    write(ifile,*) 'ZONE N=',grid%npoin,', E=',grid%nelem,', F=FEPOINT, ET=TETRAHEDRON'

    do i=1, grid%npoin
      write(ifile,'(3es14.6)') grid%x(i,1:nbDim)
    enddo

    do i=1, grid%nelem
      write(ifile,*) grid%elem(i)%face(idx,1:grid%elem(i)%flen)
    enddo
    close(ifile)

  end subroutine PlotMesh3D
!copied to geom!!!
!  subroutine GetFaceIndexes(faces, nodes, j, flen)
!    integer, dimension(:,:), intent(in) :: faces
!    integer, dimension(:), intent(out) :: nodes
!    integer, intent(in) :: flen
!    integer, intent(in) :: j
!
!    if (nbDim == 2) then
!     nodes(1)= faces(idx,j)
!     nodes(2)= faces(idx,mod(j,flen)+1)
!    else
!      if(flen .ne. 4) then !tetrahedra
!        print*,'Only tetrahedra implmented in GetFaceIndexes'
!        stop
!      endif
!        if (j == 1) then
!          nodes(1)= faces(idx,2)
!          nodes(2)= faces(idx,3)
!          nodes(3)= faces(idx,4)
!        elseif (j == 2) then
!          nodes(1)= faces(idx,3)
!          nodes(2)= faces(idx,1)
!          nodes(3)= faces(idx,4)
!        elseif (j == 3) then
!          nodes(1)= faces(idx,1)
!          nodes(2)= faces(idx,2)
!          nodes(3)= faces(idx,4)
!        else
!          nodes(1)= faces(idx,3)
!          nodes(2)= faces(idx,2)
!          nodes(3)= faces(idx,1)
!        endif
!    endif
!
!  end subroutine GetFaceIndexes

end module mesh_oper3D



