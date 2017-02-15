module mesh3D_mod
   use element_mod
   use lapack_oper
   use mesh_mod


implicit none

   type, public, extends(AbstrMesh_t) :: Mesh3D_t
      class( Element_t ), allocatable, dimension(:) :: elem     ! elements
      integer :: nelem, npoin, nbelm, nbc, n3elem, max_el_len
      integer :: curved_deg, num_curv
      real :: xper(2,2)                       ! for periodicity
      integer :: iper(2,2)                    ! for periodicity
      integer :: periodic
      real, allocatable, dimension(:,:) :: x      ! nodes coordinates
      integer, allocatable, dimension(:,:) :: xcur ! info for curved nodes
      type(Bound_edge_t), allocatable, dimension(:) :: b_edge ! boundary edges
      real :: diam                            ! diameter of the computational domain
      real :: h                               !diam of the largest element

      real :: domain_volume                    ! meassure of the computational domain

      integer :: stop_adaptation
      integer :: adapt_level               ! level of mesh adaptation
      integer :: max_adapt_level           ! maximal level of mesh adaptation ??
      logical :: adapt, adaptR             ! adapt the mesh ?
      real :: tol_min, tol_max             ! tolerances for mesh adaptation


   contains
      procedure :: init => initMesh3D
      procedure :: read => readMesh3D
      procedure :: allocElem => allocElemMesh3D

   end type Mesh3D_t

contains

!!!!!!!!!!!!!!!!!! 3D MESH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> allocates datas of the mesh
   subroutine initMesh3D( this, nelem, npoin, nbelm, nbc, nbdim )
      class(Mesh3D_t) :: this
      integer, intent(in) :: nelem, npoin, nbelm, nbc, nbdim
      integer :: i

      this%nelem = nelem
      this%npoin = npoin
      this%nbelm = nbelm
      this%nbc = nbc


      allocate(this%x(1:npoin,1:nbDim) )
      allocate(this%xcur(1:npoin,1:2) ) ! idx = 1 CBP, idx =2 pnt

      ! allocate elements of the mesh of the desired type
      print*, 'calling allocElem'
      call this%allocElem(nelem)
     ! allocate( ElementHG_t :: this%elem(1:nelem) )



      !init elem
      print*, 'Init elements after the gridfile is read!'
!      do i = 1, nelem
!!        call  this%elem%initElement(i)
!         call this%elem(i)%init(i)
!      enddo

   !   allocate(this%b_edge(1:this%nbelm) )


   end subroutine initMesh3D



     !> read the mesh from gridfile in the "grid" file
  subroutine readMesh3D(grid, gridfile)
    class(Mesh3D_t), intent(inout), target :: grid
    character(len=*), intent(in) :: gridfile
    type(Element_t), pointer :: elem
    integer :: ifile=12
    integer, parameter :: maxlen=20;
    integer:: i, lenloc, iloc(1:maxlen), counter
    real, dimension(1:3) :: vec1, vec2, vec3
    real, dimension(1:2, 1:nbDim) :: xminmax
    !real, dimension(1:2, 1:2) :: xminmax
    real :: tmp
    integer :: npoin, nelem, nbelm, nbc
    !real, pointer :: p_h

    print*, 'readMesh3D called - not tested!!!!!!!!!!'
    stop

    open(ifile, file=gridfile, status='OLD')
    read(ifile,*) npoin, nelem, nbelm, nbc

    !allocates datas of the mesh
    call grid%initMesh(nelem, npoin, nbelm, nbc, nbdim )

    read(ifile,*) grid%xper(1,1:nbDim),grid%iper(1,1:2),&
          grid%xper(2, 1:nbDim),grid%iper(2,1:2)

    ! periodicity
    grid%periodic = 0
    if(grid%xper(1,1) /= 0. .or. grid%xper(1,2) /= 0.) then
       grid%periodic = 1
       if(grid%xper(2,1) /= 0. .or. grid%xper(2,2) /= 0.) grid%periodic = 2
    endif
    if(nbDim == 3) print*,'Attention for periodicity in 3D in mesh.f90 !!!!'

    do i=1,grid%npoin
       read(ifile,*) grid%x(i,1:nbDim)
    enddo


    do i=1,nbDim
       xminmax(1, i) = minval( grid%x(:, i) )
       xminmax(2, i) = maxval( grid%x(:, i) )
    end do

    grid%diam = VectorNorm( xminmax(2,:) - xminmax(1,:) )

    do i = 1, nelem

!         !3D case,  tetrahedra, pyramid or hexahedra
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

    enddo

!
!    !FR allocation of elem should be in  a subroutine
!!    do i=1, grid%nelem
!!       elem => grid%elem(i)
!!       elem%i = i
!!
!!       read(ifile, *) lenloc, iloc(1:lenloc)
!!       if(abs(lenloc) > maxlen) then
!!          print*,'Dimensional error in ReadMesh'
!!          stop
!!       endif
!!!TODO not yet implemented
!! !      call elem%setFaces( lenloc, iloc )
!!
!!
!!    enddo
!
!    allocate(grid%b_edge(1:grid%nbelm) )
!    do i=1,grid%nbelm
!       read(ifile,*)  grid%b_edge(i)%lbn(1:2),grid%b_edge(i)%ibc
!    enddo
    close(ifile)

  end subroutine readMesh3D


   subroutine allocElemMesh3D(this, nelem)
      class ( Mesh3D_t ), intent( inout ) :: this
      integer, intent( in ) :: nelem
      print*, '3DMesh not tested'
      stop 'allocElemMesh3D - This one should never be called!'

   end subroutine allocElemMesh3D

end module mesh3D_mod
