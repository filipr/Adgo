!> definition of the computational grids
module mesh_mod
   use AMAdata
 ! use matrix_oper
  use main_data
 ! use f_mapping

  use element_mod
  use f_mapping
  !use mesh_oper
  use lapack_oper

  implicit none

  public :: EikonalTopology
!  public :: mesh


   !> boundary edge, from input and the curved representation, if any
   type, public :: Bound_edge_t
     integer :: lbn(4)        ! indexes of end nodes
     integer :: ibc           ! index of boundary element
     integer :: itc           ! index of adjacent element
     integer :: jtc           ! inner index of adjacent element
     integer :: BC            ! type of boundary condition
     integer :: inout         ! 0=inlet , 1=outlet, -1=solid walls
     integer:: icurv          ! icurv = 0, polygonal, icurv =k, P_{k+1} approx
     real, allocatable, dimension(:,:) :: x_inn ! coordinates of inner nodes
     real, allocatable, dimension(:,:) :: x_div ! coordinates of integ nodes for c_L,..

   contains
      procedure :: copy => copyBound_edge

   end type Bound_edge_t


   type, public,abstract :: Abstrmesh

   contains
      procedure(allocElem), deferred :: allocElem
      procedure(initAbstrMesh), deferred :: init

   end type Abstrmesh

   abstract interface
   subroutine allocElem( this , nelem )
      import :: Abstrmesh
      class ( Abstrmesh ), intent( inout ) :: this
      integer, intent( in ) :: nelem
   end subroutine allocElem

   subroutine initAbstrMesh( this, nelem, npoin, nbelm, nbc, nbdim )
      import :: Abstrmesh
      class(Abstrmesh), intent(inout) :: this
      integer, intent(in) :: nelem, npoin, nbelm, nbc, nbdim
   end subroutine initAbstrMesh
   end interface

   ! 2D mesh
   type, public, extends( Abstrmesh ) :: mesh
      class( element ), pointer, dimension(:) :: elem     ! elements
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
      logical :: submesh                   ! does the mesh have a partition? (specified in .ini file)
      real :: tol_min, tol_max             ! tolerances for mesh adaptation

      class(element), dimension(:), pointer :: elemP  ! for the solution of the local problems
      integer, allocatable,  dimension(:,:) :: elemP_idx  ! indexes of neighbouring  elements
      integer     :: elemP_idx_actual       ! index to actually considered element
      integer     :: flenP                  ! the maximal number of neighbouring triangles
      logical :: ElemSupports               ! the arrays elem%supp are already allocated

      ! variable and arrays for eikonal equations
      integer :: max_nbp
      integer, dimension(:,:), allocatable :: loc_ibp
      contains
!      public:: SaveOrigMesh
!      public:: RemoveOrigMesh
!      public:: AdaptMesh_AMA
!      public:: ReprepareProblem
!      public:: RecomputeElemSol
!      public:: RecomputeElemSolDerefine
!      public:: PassElem2Elem
!      public:: DeallocateElem
!      public:: Mesh2Amesh
!      public:: Amesh2Mesh
!      public:: PlotAmesh1
  !   procedure :: setMesh ! dependent on the method used
  !   procedure :: delMesh
  !   procedure :: copyMesh
  !    procedure :: writeMesh
      !passGrid to another grid (hp)
     ! AdaptMesh_hp_HG EstimateNewGrid RedHGrefinementRedHGderefinement FindMiddleHGNode PassGrid2hpGrid PasshpGrid2Grid DRefineGrid RemoveHangingNodes PlotMesh_hp DeletehpGrid

   !   generic, public :: assignment(=) => copyMesh


      procedure :: allocElem => allocElemMesh

      procedure :: init => initMesh
      procedure :: read => readmesh
      procedure :: plot => plotmesh
      procedure :: plotSubmesh
      procedure :: plotVertexFunction3D
      procedure :: seekNeighbours => seekNeighboursNew
      procedure :: setMesh
      procedure :: shiftMesh ! moves the mesh and magnifies it
!      procedure :: SeekCurvedBoundary => seekCurvedBoundaryNew
      procedure :: SetConvexSubmesh
      procedure :: SetTwoConvexSubmeshes
      procedure :: SetSubmesh ! gives iteger Values to elements according to submesh file
      procedure :: computeGeometry => computeGeometryNew
      procedure :: ComputeEikonalTopology => ComputeEikonalTopologyPPP
!      procedure :: copy => copyMesh !not used

      procedure :: setNewElementDofs ! p-adaptation DWR
      procedure :: setEdgeQuadratureDegrees ! new for SetEdgeDegrees

   end type mesh

   type, public, extends( mesh ) :: MeshHG_t
      integer :: max_HGlevel               ! maximal implemented level of hanging nodes          !!moved to %space in ADGo
      integer :: HG                        ! maximal number of hanging nodes per edge
      contains

      procedure :: allocElem => allocElementsHG

   end type MeshHG_t

   type, public, extends( mesh ) :: MeshRG_t
      integer :: max_RGlevel               ! maximal implemented level of RD
      real, dimension(:,:,:), allocatable :: RGred ! red green refinement
      real, dimension(:,:,:,:), allocatable :: RGgreen ! red green refinement

      contains

      procedure :: allocElem => allocElementsRG

   end type MeshRG_t

   type, public, extends( mesh ) :: MeshAMA_t
     ! integer, dimension(:,:), allocatable :: AMAhistory   !history of AMA computation
     !type(AMAmesh), allocatable :: AMA


      contains
      procedure :: allocElem => allocElementsAMA

   end type MeshAMA_t

!   type, public, extends( MeshHG_t ) :: MeshHGhp_t
!   end type MeshHGhp_t
!
!   type, public, extends( MeshHG_t ) :: MeshHGh_t
!   end type MeshHGh_t
!
!   type, public, extends( MeshHG_t ) :: MeshHGp_t
!   end type MeshHGp_t
!
!   type, public, extends( MeshAMA_t) :: MeshAMAhp_t
!   end type MeshAMAhp_t
!
!   type, public, extends( MeshAMA_t) :: MeshAMAh_t
!   end type MeshAMAh_t
!
!   type, public, extends( MeshAMA_t) :: MeshAMAp_t
!   end type MeshAMAp_t




   contains

   !> copy Bound_edge_t from old Bound_edge_t
   subroutine copyBound_edge( this, old)
      class( Bound_edge_t ), intent(inout) :: this
      class( Bound_edge_t ), intent(in) :: old
      integer :: i, bound1, bound2

      do i =1,4
         this%lbn(i) = old%lbn(i)
      enddo
      this%ibc = old%ibc
      this%itc = old%itc
      this%jtc = old%jtc
      this%BC = old%BC
      this%inout = old%inout
      this%icurv = this%icurv

      bound1 = size(old%x_inn(:,1))
      bound2 = size(old%x_inn(1,:))
      if( allocated(this%x_inn) ) &
         deallocate(this%x_inn)
      allocate( this%x_inn(1:bound1,1:bound2), source = old%x_inn(1:bound1,1:bound2) )

      if( allocated(this%x_div) ) &
         deallocate(this%x_div)
      bound1 = size(old%x_div(:,1))
      bound2 = size(old%x_div(1,:))
      allocate( this%x_div(1:bound1,1:bound2), source = old%x_div(1:bound1,1:bound2) )

   end subroutine copyBound_edge


   subroutine allocElemMesh(this, nelem)
      class ( mesh ), intent( inout ) :: this
      integer, intent( in ) :: nelem

      print*, 'Allocating elems - noAdapt'
      allocate( element :: this%elem(1:nelem) )

   end subroutine allocElemMesh

   !> not used
   !> copies the whole structure of oldMesh to new mesh this
!   subroutine copyMesh(this, oldMesh)
!      class (mesh), intent(inout) :: this
!      class (mesh), intent(in) :: oldMesh
!      integer :: i
!
!      print*, 'Copying meshes'
!
!      ! allocate new elements and set basic settings, alloc b_edge, x, x_cur
!      call this%init(oldMesh%nelem, oldMesh%npoin, oldMesh%nbelm, oldMesh%nbc, nbdim)
!      this%curved_deg = oldMesh%curved_deg
!      this%num_curv = oldMesh%num_curv
!      this%xper(1:2,1:2) = oldMesh%xper(1:2,1:2)
!      this%iper(1:2,1:2) = oldMesh%iper(1:2,1:2)
!      this%periodic = oldMesh%periodic
!
!      allocate( this%x( this%npoin, 1:nbDim) )
!      allocate( this%xcur(1:this%npoin, 1:nbDim) )
!
!      this%x( this%npoin, 1:nbDim) = oldMesh%x( this%npoin, 1:nbDim)
!      this%xcur(1:this%npoin, 1:nbDim) = oldMesh%xcur(1:this%npoin, 1:nbDim)
!
!      allocate( this%b_edge(1:this%nbelm) )
!      do i = 1, this%nbelm
!         call this%b_edge(i)%copy( oldMesh%b_edge(i) )
!      end do !i
!
!      this%diam = oldMesh%diam
!      this%h = oldMesh%h
!      this%domain_volume = oldMesh%domain_volume
!      this%stop_adaptation = oldMesh%stop_adaptation
!      this%adapt_level = oldMesh%adapt_level
!      this%max_adapt_level = oldMesh%max_adapt_level
!      this%adapt = oldMesh%adapt
!      this%adaptR = oldMesh%adaptR
!      this%tol_min = oldMesh%tol_min
!      this%tol_max = oldMesh%tol_max
!
!      ! elemP and other P variables are not copied now
!
!      do i =1, this%nelem
!         call this%elem(i)%copy( oldMesh%elem(i), state%time%disc_time )
!      end do
!
!   end subroutine copyMesh

   !> allocates datas of the mesh
   subroutine initMesh( this, nelem, npoin, nbelm, nbc, nbdim)
      class(mesh), intent(inout) :: this  ! target ???
      integer, intent(in) :: nelem, npoin, nbelm, nbc, nbdim
   !   integer, intent(in), optional :: HG
      integer :: i

      this%nelem = nelem
      this%npoin = npoin
      this%nbelm = nbelm
      this%nbc = nbc

      allocate(this%x(1:npoin,1:nbDim) )
      allocate(this%xcur(1:npoin,1:2) ) ! idx = 1 CBP, idx =2 pnt

      allocate( this%b_edge( 1:this%nbelm ) )

      !print*, 'INit mesh xcur allocated!'

      ! allocate elements of the mesh of the desired type
      call this%allocElem(nelem)

      select type (this)
         type is (MeshHG_t )
            this%max_HGlevel = state%space%adapt%max_HGlevel
            this%HG = state%space%adapt%HG
      end select


   end subroutine initMesh



  !> read the mesh from gridfile in the "grid" file
  subroutine readmesh( this, gridfile)
    class(mesh), intent(inout), target :: this
    character(len=*), intent(in) :: gridfile
  !  integer, intent(in), optional :: HG
    class(element), pointer :: elem
    integer :: ifile=12
    integer, parameter :: maxlen=20;
    integer:: i, lenloc, iloc(1:maxlen), counter, last_face
    real, dimension(1:3) :: vec1, vec2, vec3
    real, dimension(1:2, 1:nbDim) :: xminmax
    !real, dimension(1:2, 1:2) :: xminmax
    real :: tmp, xi(1:2)
    integer :: npoin, nelem, nbelm, nbc
    !real, pointer :: p_h

!
!    if ( .not. present(HG) ) then
!      select type(this)
!         type is (MeshHG_t)
!            stop 'Problem in readMesh'
!      end select
!    endif


    open(ifile, file=gridfile, status='OLD')
    read(ifile,*) npoin, nelem, nbelm, nbc

    !allocate data of the mesh
    call this%init(nelem, npoin, nbelm, nbc, nbdim )

    read(ifile,*) this%xper(1,1:nbDim),this%iper(1,1:2),&
          this%xper(2, 1:nbDim),this%iper(2,1:2)



    ! periodicity
    this%periodic = 0
    if(this%xper(1,1) /= 0. .or. this%xper(1,2) /= 0.) then
       this%periodic = 1
       if(this%xper(2,1) /= 0. .or. this%xper(2,2) /= 0.) this%periodic = 2
    endif
    if(nbDim == 3) print*,'Attention for periodicity in 3D in mesh.f90 !!!!'

    do i=1,this%npoin
       read(ifile,*) this%x(i,1:nbDim)
    enddo


    do i=1,nbDim
       xminmax(1, i) = minval( this%x(:, i) )
       xminmax(2, i) = maxval( this%x(:, i) )
    end do

    this%diam = VectorNorm( xminmax(2,:) - xminmax(1,:) )

    select type (this)
    ! mesh with HG nodes
    type is (MeshHG_t)
       do i = 1, nelem
         elem => this%elem(i)

         read(ifile, *) lenloc, iloc(1:lenloc)


         ! nbdim = 2 otherwise mesh3D should be used
         if(lenloc >= 3 .and. lenloc <= 4) then ! triangle or quadrilateral
            call elem%init(i, lenloc, lenloc, iloc(1:lenloc) )

         elseif(lenloc >= 9) then  ! element with hanging nodes
            !last_face = iloc(2) + 2
            call elem%init(i, iloc(1), iloc(2), iloc(3:lenloc) )

            select type (elem)
               type is (ElementHG_t)
                  call elem%setHG( this%npoin, this%x, this%max_HGlevel, this%HG )
               class default
                  stop 'stopping in readMesh - hanging nodes are possible only for ElementHG_t'
            end select

         else
                print*,'Bad data in init file ',gridfile,',  lenloc = ',lenloc
                stop
         endif
       enddo
    !mesh with NO HG nodes
    class default
      print*, '# No HG nodes'
       do i = 1, nelem
         elem => this%elem(i)
         read(ifile, *) lenloc, iloc(1:lenloc)
         !if (lenloc >= 9) stop 'this mesh shouldnt have HG nodes, in readMesh'

         ! nbdim = 2 otherwise mesh3D should be used
            call elem%init(i, lenloc, lenloc, iloc(1:lenloc) )
       enddo
    end select

    if ( .not. allocated(this%b_edge) ) &
      stop 'Stopping in readMesh - b_edge is not allocated'

     do i=1,this%nbelm
        read(ifile,*)  this%b_edge(i)%lbn(1:2), this%b_edge(i)%ibc

    !    if(this%b_edge(i)%ibc == 4) then
    !       xi(1:2) = this%x(this%b_edge(i)%lbn(1) , 1:2)
    !       if(xi(1) >= 0. .and. xi(1) <= 1.) then
    !          if(xi(2) > 0) write(88, *) xi(1:2)!, this%b_edge(i)%ibc
    !          if(xi(2) < 0) write(89, *) xi(1:2)!, this%b_edge(i)%ibc
    !       endif

    !    endif

     enddo

    ! stop 'prof preparation'
!!!TODO not yet implemented
!! !      call elem%setFaces( lenloc, iloc )
!!

    close(ifile)

    this%ElemSupports  = .false.

  end subroutine readmesh

  !> set the mesh from given paramets
  !> use after init() subroutine !!!
  ! an exmple of use follows FR
  subroutine setMesh( this, xper, iper, x, lenloc, iloc )
    class(mesh), intent(inout) :: this ! targer ???
    real, dimension(1:2,1:nbDim), intent(in) :: xper
    integer, dimension(1:2,1:nbDim), intent(in) :: iper
    real, dimension(1:this%npoin,1:nbDim), intent(in) :: x
    integer, dimension(1:this%nelem), intent(in) :: lenloc
    integer, dimension(:,:), intent(in) :: iloc

    class(element), pointer :: elem
!    integer, parameter :: maxlen=20;
!    integer:: i, lenloc, iloc(1:maxlen), counter, last_face
!    real, dimension(1:3) :: vec1, vec2, vec3
    real, dimension(1:2, 1:nbDim) :: xminmax
!    real :: tmp, xi(1:2)
    integer :: npoin, nelem, i
    !real, pointer :: p_h

    print*, '####### B_edge????'
!          this%b_edge(i)%lbn(1:2), this%b_edge(i)%ibc


    !allocate data of the mesh
!    call this%init(nelem, npoin, nbelm, nbc, nbdim )

    npoin = this%npoin
    nelem = this%nelem


    ! periodicity
    this%xper(1:2,1:nbDim) = xper(1:2,1:nbDim)
    this%iper(1:2,1:nbDim) = iper(1:2,1:nbDim)
    this%periodic = 0
    if(this%xper(1,1) /= 0. .or. this%xper(1,2) /= 0.) then
       this%periodic = 1
       if(this%xper(2,1) /= 0. .or. this%xper(2,2) /= 0.) this%periodic = 2
    endif
    if(nbDim == 3) stop 'Attention for periodicity in 3D in mesh.f90 !!!!'

    ! coordinates of the vertices of the mesh
    this%x(1:this%npoin,1:nbDim) = x(1:this%npoin,1:nbDim)
    do i=1,nbDim
       xminmax(1, i) = minval( this%x(:, i) )
       xminmax(2, i) = maxval( this%x(:, i) )
    end do
    this%diam = VectorNorm( xminmax(2,:) - xminmax(1,:) )

      ! init elements
    select type (this)
    ! mesh with HG nodes
    type is (MeshHG_t)
       do i = 1, nelem
         elem => this%elem(i)
         ! nbdim = 2 otherwise mesh3D should be used
         if(lenloc(i) >= 3 .and. lenloc(i) <= 4) then ! triangle or quadrilateral
            call elem%init(i, lenloc(i), lenloc(i), iloc(i,1:lenloc(i) ) )

         elseif(lenloc(i) >= 9) then  ! element with hanging nodes
            !last_face = iloc(2) + 2
            call elem%init(i, iloc(i,1), iloc(i,2), iloc(i,3:lenloc(i)) )

            select type (elem)
               type is (ElementHG_t)
                  call elem%setHG( this%npoin, this%x, this%max_HGlevel, this%HG )
               class default
                  stop 'stopping in readMesh - hanging nodes are possible only for ElementHG_t'
            end select
         else
            print*,'Bad data in init file ',gridfile,',  lenloc(i) = ',lenloc(i)
            stop
         endif
       enddo
    !mesh with NO HG nodes
    class default
!      print*, '# No HG nodes', nelem
       do i = 1, nelem
         elem => this%elem(i)
         ! nbdim = 2 otherwise mesh3D should be used
         call elem%init(i, lenloc(i), lenloc(i), iloc(i,1:lenloc(i)) )
       enddo
    end select

    if ( .not. allocated(this%b_edge) ) &
      stop 'Stopping in readMesh - b_edge is not allocated'

!    do i=1,this%nbelm
!      read(ifile,*)  this%b_edge(i)%lbn(1:2), this%b_edge(i)%ibc
!
!    !    if(this%b_edge(i)%ibc == 4) then
!    !       xi(1:2) = this%x(this%b_edge(i)%lbn(1) , 1:2)
!    !       if(xi(1) >= 0. .and. xi(1) <= 1.) then
!    !          if(xi(2) > 0) write(88, *) xi(1:2)!, this%b_edge(i)%ibc
!    !          if(xi(2) < 0) write(89, *) xi(1:2)!, this%b_edge(i)%ibc
!    !       endif
!
!    !    endif
!
!     enddo



    this%ElemSupports  = .false.

  end subroutine setMesh
  !      nelem = 4
!      npoin = 5
!      nbelm = 4
!      call this%grid%init(nelem, npoin, nbelm, 0, nbDim)
!
!      eps = this%eps
!      x(1, 1:nbDim) = this%xy_coord(1:nbDim)
!      x(2, 1:nbDim) = (/ this%xy_coord(1) - eps , this%xy_coord(2) - eps/)
!      x(3, 1:nbDim) = (/ this%xy_coord(1) + eps , this%xy_coord(2) - eps/)
!      x(4, 1:nbDim) = (/ this%xy_coord(1) + eps , this%xy_coord(2) + eps/)
!      x(5, 1:nbDim) = (/ this%xy_coord(1) - eps , this%xy_coord(2) + eps/)
!
!
!
!      xper(1:2,1:nbDim) = 0.0
!      iper(1:2,1:nbDim) = 0
!
!      allocate( lenloc(1:nelem), source = 3 )
!      allocate( iloc(1:nelem, 1:3), source = 0 )
!
!      if (nelem /= 4 ) &
!         stop 'problem in initGrid_point_value'
!      iloc( 1 , 1:lenloc(1) ) = (/ 1,2,3 /)
!      iloc( 2 , 1:lenloc(2) ) = (/ 1,3,4 /)
!      iloc( 3 , 1:lenloc(3) ) = (/ 1,4,5 /)
!      iloc( 4 , 1:lenloc(4) ) = (/ 1,5,2 /)
!
!
!      call this%grid%setMesh( xper(1:2,1:nbDim), iper(1:2,1:nbDim), &
!            x(1:npoin,1:nbDim), lenloc(1:nelem), iloc(1:nelem,1:3) )


      !set boundary edges
!      do i = 1,nbelm
!         this%grid%b_edge(i)%lbn(1:2) = (/ i+1, mod(i, nbelm) + 2 /)
!         this%grid%b_edge(i)%ibc = this%grid%elem(i)%i
!      end do


  !> plot the mesh in the file 'meshfile' visualizable by gnuplot
  subroutine plotmesh( this, meshfile )
    class(mesh), intent(inout) :: this
    character(len=*), intent(in) :: meshfile
    integer:: i,j, k, ifile = 11

    open(ifile, file=meshfile)

    !!print*,'@@@',this%nelem

    do i=1, this%nelem
       do j=0,this%elem(i)%flen
          k = this%elem(i)%face( idx, mod( j,this%elem(i)%flen ) +1 )
          write(ifile,'(2es14.6,2i6)') this%x(k,1:nbDim),0,i
       enddo
       write(ifile,'(x)')
       write(ifile,'(x)')
    enddo
    close(ifile)

  end subroutine plotmesh

  !> plot the mesh in the file 'meshfile' visualizable by gnuplot
  !> elements with %iSubmesh == iPart are plotted
  subroutine plotSubmesh( this, iPart, meshfile )
    class(mesh), intent(inout) :: this
    integer, intent(in) :: iPart ! elements with %iSubmesh == iPart are plotted
    character(len=*), intent(in) :: meshfile
    integer:: i,j, k, ifile = 11

    open(ifile, file=meshfile)

    !!print*,'@@@',this%nelem

    do i=1, this%nelem
      if (this%elem(i)%iSubMesh == iPart ) then
       do j=0,this%elem(i)%flen
          k = this%elem(i)%face( idx, mod( j,this%elem(i)%flen ) +1 )
          write(ifile,'(2es14.6,2i6)') this%x(k,1:nbDim),0,i
       enddo
       write(ifile,'(x)')
       write(ifile,'(x)')
      end if
    enddo
    close(ifile)

  end subroutine plotSubmesh

   !> shift a point to the given coordinates
   !> then magnifies the mesh around it according to the parameter eps
   !> only the coordinates of the mesh points are changed
   subroutine shiftMesh( this, iPoint, coord, eps)
      class(mesh), intent(inout) :: this
      integer, intent(in) :: iPoint
      real, dimension(1:nbDim), intent(in) :: coord
      real, intent(in) :: eps
      real, dimension(1:nbDim) :: xPoint
      integer :: i

      xPoint(1:nbDim) = this%x( iPoint, 1:nbDim )

!      shift(1:nbDim) = coord(1:nbDim) - this%x( iPoint, 1:nbDim )
      ! shift the chosen point
      this%x( iPoint, 1:nbDim ) = coord(1:nbDim)

      ! shift the other point and change their distance from iPoint by multiplying by eps
      do i=1, iPoint - 1
         this%x(i,1:nbDim) = this%x(iPoint,1:nbDim) + &
            eps * ( this%x(i,1:nbDim) - xPoint(1:nbDim) )
      end do
      do i = iPoint + 1, this%npoin
         this%x(i,1:nbDim) = this%x(iPoint,1:nbDim) + &
            eps * ( this%x(i,1:nbDim) - xPoint(1:nbDim) )
      end do

   end subroutine shiftMesh

  !> seeking of neighbouring of elements \f$ K, K',\ K\cap K'\not=\emptyset\f$
  subroutine seekNeighboursNew(this)
    class( mesh ), intent(inout) :: this
    class(element), pointer :: elem
    integer, dimension(:,:), allocatable :: cyc
    integer, dimension(:), allocatable :: len_cyc
    integer, parameter :: max_cyc = 20
    integer :: i, j, k, k1, l, ni, nj, is1, is2, ip1, ip2, ii, ie, iie, jj, ib, ibound, j1, j2, k2

    !print* , 'seekNeighboursNew, ??? - suma:', sum( this%elem(1:this%nelem)%type - 2 )

    if (nbDim == 3) then
       ! only tetrahedra
       this%n3elem = sum(this%elem(:)%type-3)
    else
       this%n3elem = sum(this%elem(:)%type-2)
    endif
    this%max_el_len = maxval( this%elem(:)%type)

    ! to each vertex we create a list of elements sharing this node
    allocate(cyc(this%npoin, max_cyc))
    allocate(len_cyc(this%npoin))

    len_cyc(1:this%npoin) = 0
    do  i=1, this%nelem
       do j=1, this%elem(i)%flen
          k = this%elem(i)%face(idx,j)

          len_cyc(k) = len_cyc(k) + 1
          cyc(k,len_cyc(k)) = i
       enddo
       this%elem(i)%per = 0
    enddo

    ! inicialization  of arrays for neighbours

    do i=1, this%nelem
       this%elem(i)%face(neigh,:) = -10
       this%elem(i)%face(nei_i,:) = -10
    enddo

    do i=1, this%nelem
       ibound = 0
       do j=1,this%elem(i)%flen
          if(this%elem(i)%face(neigh,j) == -10) then
             k = this%elem(i)%face(idx,j)
             k1 = this%elem(i)%face(idx,mod(j,this%elem(i)%flen)+1)
             do l=1,len_cyc(k)
                ni = cyc(k,l)
                if(i /= ni) then
                   do nj =1,this%elem(ni)%flen
                      if(k1 == this%elem(ni)%face(idx,nj) .and.   &
                           k== this%elem(ni)%face(idx,mod(nj,this%elem(ni)%flen)+1))&
                           then
                         this%elem(i)%face(neigh,j) = ni
                         this%elem(ni)%face(neigh,nj) = i

                         this%elem(i)%face(nei_i,j) = nj
                         this%elem(ni)%face(nei_i,nj) = j

                         goto 100
                      endif
                   enddo
                endif
             enddo
100          continue
          endif
          ibound = ibound + 1
       enddo
       if(ibound > 0) then ! boundary element
          !print*,'Alloc:', this%elem(i)%i
          allocate(this%elem(i)%iBC(1:this%elem(i)%flen) )
          this%elem(i)%iBC(:) = 0

          allocate(this%elem(i)%tBC(1:this%elem(i)%flen) )
          this%elem(i)%tBC(:) = 0
       endif
    enddo


    !do i= 2,6
    !   write(*,'(15i5)') i, this%elem(i)%face(idx,:), this%elem(i)%face(neigh,:), &
    !        this%elem(i)%face(nei_i,:)
    !enddo

    ! seeking of element adjacent to boundary elements
    this%b_edge(1:this%nbelm)%itc = -1
    do i=1,this%nelem
       do j=1,this%elem(i)%flen
          if(this%elem(i)%face(neigh,j) < 0 ) then
             k = this%elem(i)%face(idx,j)
             k1 = this%elem(i)%face(idx,mod(j,this%elem(i)%flen)+1)
             do ib=1,this%nbelm
                if(this%b_edge(ib)%itc < 0) then
                   !   print*,'?',k,k1,this%b_edge(ib)%lbn(1:2)
                   if(this%b_edge(ib)%lbn(1) == k .and. &
                        this%b_edge(ib)%lbn(2) == k1 ) then
                      this%b_edge(ib)%itc = i
                      this%b_edge(ib)%jtc = j

                      this%elem(i)%face(neigh,j) = -ib
                      goto 200
                   endif
                endif
             enddo
             print*,'Adjacent element to a boundary edge doesnt found'
             print*,'elem:',i,j,this%elem(i)%flen, k, k1
             print*,'elem:',i,0,this%elem(i)%face(idx,:)
             !print*,'elem:',18,0,this%elem(18)%face(idx,:)
             print*, this%x(k,1:nbDim)
             print*, this%x(k1,1:nbDim)
             stop
          endif
200       continue
          ! print*,'####','elem :',i,j,'nodes:',k,k1, '  iBe=',ib
       enddo
    enddo

    do l=1, this%periodic
       is1 = 0
       is2 = 0
       ip1 = 0
       ip2 = 0
       do i=1,this%nbelm
          if(is1 == 0 .and. this%b_edge(i)%ibc == this%iper(l,1) ) is1 = i
          if(this%b_edge(i)%ibc == this%iper(l,1) ) ip1 = ip1 + 1

          if(is2 == 0 .and. this%b_edge(i)%ibc == this%iper(l,2) ) is2 = i
          if(this%b_edge(i)%ibc == this%iper(l,2) ) ip2 = ip2 + 1

       enddo
       write(*,'(a4,7i5)') 'per', l,this%iper(l,1:2), is1, is2, ip1, ip2

       do k = 0, ip1 - 1
          i = is1 + k
          ii = is2 + ip2 - k - 1

          !print*,'@@@', k,i, ii

          ie  = this%b_edge( i)%itc
          iie = this%b_edge(ii)%itc

          j =  this%b_edge( i)%jtc
          jj = this%b_edge(ii)%jtc

          !write(*,'(a5,8i5)') '@@@', k,i, ii, ie, iie, j, jj

          !write(101,*) (this%x(this%elem( ie)%face(idx, j), 1:nbDim) + &
          !     this%x(this%elem( ie)%face(idx, mod(j,3)+1), 1:nbDim) )/ 2
          !write(101,*) (this%x(this%elem(iie)%face(idx, mod(jj,3)+1), 1:nbDim) + &
          !     this%x(this%elem(iie)%face(idx, jj), 1:nbDim) ) / 2.
          !write(101,*)

          this%elem( ie)%per = iie
          this%elem(iie)%per = ie

          this%elem( ie)%face(neigh,  j ) = iie
          this%elem(iie)%face(neigh, jj ) = ie


          this%elem( ie)%face(nei_i,  j ) = jj
          this%elem(iie)%face(nei_i, jj ) = j
       enddo

    end do


    !print*,' End of periodicity'

    ! arrays for the eikonal equation
    if(state%modelName == 'pedes' ) then


       ! to each vertex we create a list of nodes sharing this node
       len_cyc(1:this%npoin) = 0
       cyc(:,:) = 0

       do  ie=1, this%nelem
          elem => this%elem(ie)


          do j1=1, elem%flen
             k1 = elem%face(idx,j1)

             j2= mod( j1, elem%flen) + 1
             k2 = elem%face(idx,j2)

             ! k1 and k2 is the pair of nodes with the common edge
             do l=1, len_cyc(k1)
                if(k2 == cyc(k1, l) ) goto 110 ! this node is already in the list
             enddo
             ! we add the node k2 to the k1's list
             len_cyc(k1) = len_cyc(k1) + 1
             cyc(k1,len_cyc(k1)) = k2

110          continue
          enddo
       enddo


       ! storing of the array len_cyc, cyc
       this%max_nbp = maxval(len_cyc(1:this%npoin) )
       allocate(this%loc_ibp(1:this%npoin, 0: this%max_nbp) )
       this%loc_ibp(:,:) = 0

       do  i=1, this%npoin
          j = len_cyc(i)
          this%loc_ibp(i, 0 ) = j
          this%loc_ibp(i, 1: j) = cyc(i, 1:j)

          !write(*,'(a6,i5,a2,40i5)') 'ibp:',i,'::',  this%loc_ibp(i, 0:j )
       enddo

       call this%ComputeEikonalTopology( )
       !print*,'### EIKONAL TOLOPOGY computed'
    endif


    deallocate(cyc,len_cyc)

  end subroutine seekNeighboursNew

  !> computing of geometry \f$ \forall K\in{\cal T}_h\f$
  !> after refinement, only refined elements have to be recomputed
  subroutine ComputeGeometryNew( this )
    class( mesh ), intent(inout) :: this
    integer:: i, flen,ib, fdeg

  !    print*, 'ComputeGeometryNew'

    !open(45,file='nodes')
    this%h = 0.
    this%domain_volume = 0.

    do i=1,this%nelem
!      flen = this%elem(i)%flen
!      ibcur = this%elem(i)%ibcur
      fdeg = this%elem(i)%deg_cur

!      if( fdeg <= 1 ) Fdeg =  1      ! polygonal element Fdeg = 1
!      print*, 'Fdeg:', fdeg
!       call this%elem(i)%computeElementGeometry( this%x( this%elem(i)%face( idx, 1:this%elem(i)%flen ) , 1:nbDim ), &
!            this%b_edge( this%elem(i)%ibcur )%x_inn(:,1:nbDim) )
      !curved element
      if (fdeg > 1) then
         ib =this%elem(i)%ibcur
         call this%elem(i)%computeElementGeometry( this%x , this%b_edge(ib)%x_inn )
      else
         call this%elem(i)%computeElementGeometry( this%x )
      endif
       !this%elem(i)%area = 25.
       this%domain_volume = this%domain_volume + this%elem(i)%area
       this%h = max( this%h, this%elem(i)%diam )

       !if(.not. grid%elem(i)%F%iFlin) &
       !if(grid%elem(i)%HGnode) &
       !     call CheckElement(grid%elem(i), 45)
    enddo
    ! the following (h, domain_volume) is done in state%space%copyGeometry
 !   grid%h = state%space%h
!    close(45)

!    print*,'# Geometry computed'


    !do i=1,grid%nelem
    !      do j=1,grid%elem(i)%flen
    !         if(grid%elem(i)%face(neigh,j) > 0) then
    !            ii = grid%elem(i)%face(neigh,j)
    !            write(98,*) grid%elem(i)%xc(1:nbDim)
    !            write(98,*) grid%elem(ii)%xc(1:nbDim)
    !            write(98,*)
    !            write(98,*)
    !         else
    !            write(97,*) grid%elem(i)%xc(1:nbDim)
!
!             endif
!
!          enddo
!       enddo
!
!       do k=1,grid%nbelm
!          i = grid%b_edge(k)%itc
!          j = grid%b_edge(k)%jtc
!          if(grid%elem(i)%face(neigh,j) > 0) then
!             ii = grid%elem(i)%face(neigh,j)
!             write(95,*) grid%elem(i)%xc(1:nbDim)
!             write(95,*) grid%elem(ii)%xc(1:nbDim)
!             write(95,*)
!             write(95,*)
!          else
!             write(94,*) grid%elem(i)%xc(1:nbDim)
!
!          endif
!
!       enddo
!
!       print*,'# Geometry computed, test done'
!    stop
  end subroutine ComputeGeometryNew



  !> DO NOT USE UNTIL THE PROBLEM WITH ALLOCATION IS SOLVED IN THE GCC COMPILER (bugzilla Bug 65359)
  !> allocation of HG meshes
   subroutine allocElementsHG( this, nelem)
      class( MeshHG_t), intent(inout) :: this
      integer, intent(in) :: nelem

         allocate( ElementHG_t :: this%elem(1:nelem) )

   end subroutine allocElementsHG

   !> allocation of RG meshes
   subroutine allocElementsRG( this, nelem)
      class( MeshRG_t), intent(inout) :: this
      integer, intent(in) :: nelem

      allocate( ElementRG_t :: this%elem(1:nelem) )

   end subroutine allocElementsRG


   !> allocation of AMA meshes
   subroutine allocElementsAMA( this, nelem)
      class( MeshAMA_t), intent(inout) :: this
      integer, intent(in) :: nelem

      allocate( ElementAMA_t :: this%elem(1:nelem) )

   end subroutine allocElementsAMA


   !> compute the topology for the solution of the eikonal equations
   subroutine ComputeEikonalTopologyPPP(this )
    class( mesh ), intent(inout) :: this
    real, dimension(:,:), allocatable :: x
    integer, dimension(:,:), allocatable :: loc_ibp, lnd, ibp
    integer :: nelem, npoin, max_nbp, i, j, k

    nelem = this%nelem
    npoin = this%npoin
    max_nbp = this%max_nbp

    allocate( x(1:npoin, 1:2) )
    x(1:npoin, 1:2) = this%x(1:npoin, 1:2)

    allocate( lnd(1:nelem, 1:3) )
    do i=1,nelem
       lnd(i, 1:3) = this%elem(i)%face(idx, 1:3)
    enddo

    allocate(loc_ibp(1:npoin, 0: max_nbp) )
    loc_ibp(1:npoin, 0: max_nbp) = this%loc_ibp(1:npoin, 0: max_nbp)


    allocate( ibp(1:npoin, 1:1) )
    ibp(:,:) = -1
    do i=1,this%nbelm
       j = this%b_edge(i)%ibc

       ibp( this%b_edge(i)%lbn(1) , 1) = 0
       ibp( this%b_edge(i)%lbn(2) , 1) = 0
    enddo

    do i=1,this%nbelm
       j = this%b_edge(i)%ibc
       do k=1, state%numBC
          if(state%BC(k)%ibc == j) then
             if(state%BC(k)%inout == 1) then
                ibp( this%b_edge(i)%lbn(1) , 1) = 1
                ibp( this%b_edge(i)%lbn(2) , 1) = 1
             endif
          endif
       enddo
    enddo

    !do i=1,npoin
    !   if(ibp(i, 1) == 0) write(20+state%space%adapt%adapt_level, *) x(i, 1:2)
    !   if(ibp(i, 1) == 1) write(40+state%space%adapt%adapt_level, *) x(i, 1:2)
    !enddo

    ! OUTPUT
    ! open(22, file="eikonal.data", status='replace')
    ! write(22, *) npoin, nelem,  max_nbp
    ! do i=1,npoin
    !    write(22, *) x(i, 1:2)
    ! enddo
    ! do i=1,nelem
    !    write(22, *) lnd(i, 1:3)
    ! enddo
    ! do i=1, npoin
    !    write(22, *) loc_ibp(i, 0: max_nbp  )
    ! enddo
    ! do i=1, npoin
    !    write(22, *) ibp(i, 1)
    ! enddo

    ! close(22)

    call EikonalTopology(npoin, nelem,  max_nbp, &
         x(1:npoin, 1:2), lnd(1:nelem, 1:3), loc_ibp(1:npoin, 0: max_nbp  ), ibp(1:npoin, 1) )


    !stop'ed37dehd37d3de3'
    deallocate( x,  loc_ibp, lnd, ibp )


   end subroutine ComputeEikonalTopologyPPP

   ! !> Prepares topology for the solution of the eikonal equation
   ! subroutine EikonalTopology(npoin, nelem, max_nbp, x, lnd, loc_ibp, ibp)
   !   integer, intent(in) :: npoin ! number of the grid nodes
   !   integer, intent(in) :: nelem ! number of the grid triangles
   !   integer, intent(in) :: max_nbp ! maximal number of neighbouring nodes
   !   real, dimension(1:npoin, 1:2), intent(in) :: x  ! coordinates of the nodes
   !   integer, dimension(1:nelem, 1:3), intent(in) :: lnd ! indexes of nodes forming each triangle
   !   integer, dimension(1:npoin, 0:max_nbp), intent(in) ::loc_ibp !loc_ibp(i,0) = number of neigh nodes
   !                                                                !loc_ibp(i,1:) = indexes of ^^^^^^^

   !   integer, dimension(1:npoin, 1), intent(in) ::ibp !loc_ibp(i) = -1 ==> internal node
   !                                                    !loc_ibp(i) =  0 ==> node on boundary NOT OUTFLOW
   !                                                    !loc_ibp(i) =  1 ==> node on boundary  OUTFLOW


   ! end subroutine EikonalTopology

   !> Prepares topology for the solution of the eikonal equation
   subroutine EikonalTopology(npoin, nelem, max_nbp, x, lnd, loc_ibp, ibp)
     use vertqueue ! for queue of vertices
     use BRAlgorithm
     integer, intent(in) :: npoin ! number of the grid nodes
     integer, intent(in) :: nelem ! number of the grid triangles
     integer, intent(in) :: max_nbp ! maximal number of neighbouring nodes
     real, dimension(1:npoin, 1:2), intent(in) :: x  ! coordinates of the nodes
     integer, dimension(1:nelem, 1:3), intent(in) :: lnd ! indexes of nodes forming each triangle
     integer, dimension(1:npoin, 0:max_nbp), intent(in) ::loc_ibp !loc_ibp(i,0) = number of neigh nodes
                                                                  !loc_ibp(i,1:) = indexes of ^^^^^^^

     integer, dimension(1:npoin, 1), intent(in) ::ibp !loc_ibp(i) = -1 ==> internal node
                                                      !loc_ibp(i) =  0 ==> node on boundary NOT OUTFLOW
                                                      !loc_ibp(i) =  1 ==> node on boundary  OUTFLOW

     !*******************START-MODIFIED-KUBERA**********************************
     ! prepare sigmao set
     call PrepareSigmao(npoin, max_nbp, loc_ibp(1:npoin, 0: max_nbp), ibp(1:npoin, 1))


     !print*,' prepare set PKPi  -set of edges around Pi and other infos  '
     call PrepareBKPi(npoin, nelem,  max_nbp, x(1:npoin, 1:2), lnd(1:nelem, 1:3),ibp(1:npoin, 1))


     ! prepare potential and other vertex oriented arrays
     call PrepareVertexInfo(npoin)

     !print*,' creates vertex queue for BR iteration'
     call QueueCreate(npoin)

     !print*,'*******************END-MODIFIED-KUBERA*********************************'

   end subroutine EikonalTopology

   subroutine SetSubmesh( this, linesFile )
      class( mesh ), intent(inout) :: this
      character(len=50), intent(in) :: linesFile
      character(len=50) :: submeshType
      integer :: iFile = 39
      integer :: nLines

      this%submesh = .false.

      open( iFile, file=linesFile, status='old')
      read( iFile, *) nLines
      if (nLines > 0) then
         read( iFile, *) submeshType
         ! there is a partition of the mesh
         this%submesh = .true.
      else
         submeshType = 'NONE'
      endif
      close( iFile )

!      print*, 'subMeshType: ' , submeshType


      select case (submeshType)
      case ('convex')
         write(*,*) '# Convex submesh.'
         call this%SetConvexSubmesh( linesFile )
      ! two submeshes - for primal and also for dual problem
      case ('convex2')
         write(*,*) '# Two convex submeshes.'
         call this%SetTwoConvexSubmeshes( linesFile )
      case default
         write(*,*) '# Unknown type of submesh in SetSubmesh, elem%iSubmesh = 0'
      end select


   end subroutine SetSubmesh

   !> sets the grid%elem(:)%iSubmesh to 1 - inside the Convex submesh, 2 - outside
   subroutine SetConvexSubmesh( this, linesFile )
      class( mesh ), intent(inout) :: this
      character(len=50), intent(in) :: linesFile
      integer :: iFile = 39
      integer :: nLines, i, j, ie
      real, dimension(:,:,:), allocatable :: points
      real, dimension(:,:), allocatable :: normal
      real, dimension(:), allocatable :: coef
      real, dimension(1:nbDim) :: x_inside
      character(len=50) :: submeshType
      logical :: found
      real :: eps

      eps = 1.0E-10

      open( iFile, file=linesFile, status='old')
      read( iFile, *) nLines
      read( iFile, *) submeshType

      if ( nLines > 0 ) then

         allocate( points(1:nLines, 1:2, 1:nbDim), source = 0.0 )
         allocate( normal(1:nLines, 1:nbDim), source = 0.0 )
         allocate( coef(1:nLines), source = 0.0 )

         do i = 1, nLines
            read( iFile, *) j, points(i,1,1:nbDim) , points(i,2,1:nbDim)
            !normal vector to the given line
            normal(i,1:nbDim) = (/ points(i,2,2) - points(i,1,2), points(i,1,1) - points(i,2,1) /)
            normal(i, 1:nbDim) = normal(i, 1:nbDim) / VectorNorm( normal(i, 1:nbDim)  )
            !print*, i , 'th line:', points(i,1:2,1:nbDim)
         end do



         ! find a point inside the submesh - a barycenter of a triangle in the submesh
         found = .false.
         i = 2
         do while ( (.not. found) .and. i <= nLines )
            do j = 1,2
               if ( any( abs( points(1, 1, 1:nbDim) - points(i, j, 1:nbDim) ) > eps ) ) then
                  if ( any( abs( points(1, 2, 1:nbDim) - points(i, j, 1:nbDim) ) > eps ) ) then
                     x_inside(1:nbDim) = &
                      ( points(1, 1, 1:nbDim) + points(1, 2, 1:nbDim) + points(i, j, 1:nbDim) ) / 3.0

                     !print*, 'x_inside:', x_inside(1:nbDim)
                     found = .true.
                  end if
               end if
            end do !j
            i = i+1
         end do


         !control
         if (.not. found) &
            stop 'Problem in setSubMesh - no point inside the submesh was found!'

         ! compute the general equation for the i-th line
         ! n1*x + n2*y + coef == 0
         do i = 1,nLines
            coef(i) = (-1)*dot_product( normal(i,1:nbDim), points(i,1,1:nbDim) )
            !we need the normal to be outer
            if ( dot_product( normal(i,1:nbDim), x_inside(1:nbDim) ) + coef(i) >= 0.0 ) then
               normal(i,1:nbDim) = (-1)*normal(i,1:nbDim)
               coef(i) = (-1)*coef(i)
            end if
         end do

         do ie=1,this%nelem
            i = 1
            found = .false.
            do while (.not. found .and. i <= nLines)
               ! is the element outside?
               !print*, 'elem:', ie, 'nline=', i , dot_product( normal(i,1:nbDim), this%elem(ie)%xc(1:nbDim) ) + coef(i)
               if ( dot_product( normal(i,1:nbDim), this%elem(ie)%xc(1:nbDim) ) + coef(i) >= 0.0 ) then
                  this%elem(ie)%iSubMesh = 2
                  found = .true.
               end if
               i = i+1
            end do

            !inside
            if (.not. found) &
               this%elem(ie)%iSubMesh = 1

         end do


         deallocate( points, normal, coef )

      else
         stop 'the submesh file does not specify the support of the target functional'
      endif
      close(iFile)

   end subroutine SetConvexSubmesh

   !> sets the grid%elem(:)%iSubmesh to -1 - inside the PRIMAL Convex submesh, 1 inside DUAL Convex Submesh, 2 - outside
   subroutine SetTwoConvexSubmeshes( this, linesFile )
      class( mesh ), intent(inout) :: this
      character(len=50), intent(in) :: linesFile
      integer :: iFile = 39
      integer :: nLines, i, j, ie, pLines, dLines
      real, dimension(:,:,:), allocatable :: points
      real, dimension(:,:), allocatable :: normal
      real, dimension(:), allocatable :: coef
      real, dimension(1:nbDim) :: x_inside
      character(len=50) :: submeshType
      logical :: found
      real :: eps

      eps = 1.0E-10

      open( iFile, file=linesFile, status='old')
      read( iFile, *) nLines
      read( iFile, *) submeshType
      read( iFile, *) pLines, dLines

      if (pLines + dLines /= nLines) &
         stop 'wrong number of lines in TwoConvex Submeshes'

!      print* , 'pLines:' , pLines

      if ( pLines > 0 ) then

         allocate( points(1:pLines, 1:2, 1:nbDim), source = 0.0 )
         allocate( normal(1:pLines, 1:nbDim), source = 0.0 )
         allocate( coef(1:pLines), source = 0.0 )

         do i = 1, pLines
            read( iFile, *) j, points(i,1,1:nbDim) , points(i,2,1:nbDim)
            !normal vector to the given line
            normal(i,1:nbDim) = (/ points(i,2,2) - points(i,1,2), points(i,1,1) - points(i,2,1) /)
            normal(i, 1:nbDim) = normal(i, 1:nbDim) / VectorNorm( normal(i, 1:nbDim)  )
            !print*, i , 'th line:', points(i,1:2,1:nbDim)
         end do

         ! find a point inside the submesh - a barycenter of a triangle in the submesh
         found = .false.
         i = 2
         do while ( (.not. found) .and. i <= pLines )
            do j = 1,2
               if ( any( abs( points(1, 1, 1:nbDim) - points(i, j, 1:nbDim) ) > eps ) ) then
                  if ( any( abs( points(1, 2, 1:nbDim) - points(i, j, 1:nbDim) ) > eps ) ) then
                     x_inside(1:nbDim) = &
                      ( points(1, 1, 1:nbDim) + points(1, 2, 1:nbDim) + points(i, j, 1:nbDim) ) / 3.0

!                     print*, 'x_inside:', x_inside(1:nbDim)
                     found = .true.
                  end if
               end if
            end do !j
            i = i+1
         end do


         !control
         if (.not. found) &
            stop 'Problem in setSubMesh - no point inside the submesh was found!'

         ! compute the general equation for the i-th line
         ! n1*x + n2*y + coef == 0
         do i = 1,pLines
            coef(i) = (-1)*dot_product( normal(i,1:nbDim), points(i,1,1:nbDim) )
            !we need the normal to be outer
            if ( dot_product( normal(i,1:nbDim), x_inside(1:nbDim) ) + coef(i) >= 0.0 ) then
               normal(i,1:nbDim) = (-1)*normal(i,1:nbDim)
               coef(i) = (-1)*coef(i)
            end if
         end do

         do ie=1,this%nelem
            i = 1
            found = .false.
            do while (.not. found .and. i <= pLines)
               ! is the element outside?
!               print*, 'elem:', ie, 'nline=', i , dot_product( normal(i,1:nbDim), this%elem(ie)%xc(1:nbDim) ) + coef(i)
               if ( dot_product( normal(i,1:nbDim), this%elem(ie)%xc(1:nbDim) ) + coef(i) >= 0.0 ) then
                  this%elem(ie)%iSubMesh = 2
                  found = .true.
               end if
               i = i+1
            end do

            !inside
            if (.not. found) &
               this%elem(ie)%iSubMesh = -1


         end do


         deallocate( points, normal, coef )
      else
         stop 'the submesh file does not specify the support of the primal submesh'
      endif

!      do ie=1,this%nelem
!         print*, ie, this%elem(ie)%iSubMesh
!      end do


      ! DUAL SUBMESH
      if ( dLines > 0 ) then
!         print*, 'dLines:' , dLines

         allocate( points(1:dLines, 1:2, 1:nbDim), source = 0.0 )
         allocate( normal(1:dLines, 1:nbDim), source = 0.0 )
         allocate( coef(1:dLines), source = 0.0 )

         do i = 1, dLines
            read( iFile, *) j, points(i,1,1:nbDim) , points(i,2,1:nbDim)
            !normal vector to the given line
            normal(i,1:nbDim) = (/ points(i,2,2) - points(i,1,2), points(i,1,1) - points(i,2,1) /)
            normal(i, 1:nbDim) = normal(i, 1:nbDim) / VectorNorm( normal(i, 1:nbDim)  )
            !print*, i , 'th line:', points(i,1:2,1:nbDim)
         end do

         ! find a point inside the submesh - a barycenter of a triangle in the submesh
         found = .false.
         i = 2
         do while ( (.not. found) .and. i <= dLines )
            do j = 1,2
               if ( any( abs( points(1, 1, 1:nbDim) - points(i, j, 1:nbDim) ) > eps ) ) then
               if ( any( abs( points(1, 2, 1:nbDim) - points(i, j, 1:nbDim) ) > eps ) ) then
                     x_inside(1:nbDim) = &
                      ( points(1, 1, 1:nbDim) + points(1, 2, 1:nbDim) + points(i, j, 1:nbDim) ) / 3.0

                     !print*, 'x_inside:', x_inside(1:nbDim)
                     found = .true.
                  end if
               end if
            end do !j
            i = i+1
         end do
         !control
         if (.not. found) &
            stop 'Problem in setSubMesh - no point inside the DUAL submesh was found!'

         ! compute the general equation for the i-th line
         ! n1*x + n2*y + coef == 0
         do i = 1,dLines
            coef(i) = (-1)*dot_product( normal(i,1:nbDim), points(i,1,1:nbDim) )
            !we need the normal to be outer
            if ( dot_product( normal(i,1:nbDim), x_inside(1:nbDim) ) + coef(i) >= 0.0 ) then
               normal(i,1:nbDim) = (-1)*normal(i,1:nbDim)
               coef(i) = (-1)*coef(i)
            end if
         end do

         do ie=1,this%nelem
            i = 1
            found = .false.
            !otherwise the element it already in the primal subdomain ( NO INTERSECTION IS ALLOWED )
            if ( this%elem(ie)%iSubMesh /= -1 ) then
               do while (.not. found .and. i <= dLines)
                  ! is the element outside?
                  !print*, 'elem:', ie, 'nline=', i , dot_product( normal(i,1:nbDim), this%elem(ie)%xc(1:nbDim) ) + coef(i)
                  if ( dot_product( normal(i,1:nbDim), this%elem(ie)%xc(1:nbDim) ) + coef(i) >= 0.0 ) then
                     this%elem(ie)%iSubMesh = 2
                     found = .true.
                  end if
                  i = i+1
               end do

               !inside
               if (.not. found) &
                  this%elem(ie)%iSubMesh = 1

!               if (.not. found) &
!                  print*, 'found elem and set -1'

            endif !iSubMesh /=1

         end do


         deallocate( points, normal, coef )

      else
         stop 'the submesh file does not specify the support of the primal submesh'
      endif

!     do ie=1,this%nelem
!         print*, ie, this%elem(ie)%iSubMesh
!      end do

      close(iFile)

   end subroutine SetTwoConvexSubmeshes

   !> increase(p_plus > 0)/decrease(p_plus<0) element polynomial degrees and dofs
   !> + SetElementQuadraturesDegrees, setEdgeQuadratureDegrees
   subroutine setNewElementDofs( this, p_plus )
      class(mesh), intent(inout) :: this
      integer, intent(in) :: p_plus
      class( element ), pointer :: elem
      integer :: i

      do i = 1, this%nelem
         elem => this%elem(i)
         !p+p_plus
         if ( (elem%deg + p_plus <= MaxDegreeImplemented) .and. &
               ( elem%deg + p_plus >= 0 ) ) then
            elem%deg = elem%deg + p_plus
            call elem%initElementDof()
         else
            print*, 'elem deg: new/old', elem%deg + p_plus, elem%deg
            stop 'Problem in setNewElementDofs'
         endif
       enddo

       state%space%max_dof = maxval(this%elem(:)%dof)
       !'Setting of DOF for edges (used in the numerical quadratures)'
      call this%setEdgeQuadratureDegrees()

   end subroutine setNewElementDofs


!> setting of the degrees of edge quadratures
  subroutine setEdgeQuadratureDegrees(this)
    class(mesh), intent(inout) :: this
    integer :: i, j, ii

    do i=1,this%nelem
       do j=1,this%elem(i)%flen
          if(this%elem(i)%face(neigh,j) > 0) then
             ii = this%elem(i)%face(neigh,j)
             ! maximal degree for face quadratures
             this%elem(i)%face(fdeg,j) = max(this%elem(i)%deg, this%elem(ii)%deg )
             ! dof from opposite element for matrix shape
             this%elem(i)%face(fdof,j) = this%elem(ii)%dof

             !if(state%modelName == 'scalar' .or.state%modelName == '2eqs' .and. state%space%adapt%adapt_method == 'RTN') &   ! for dual problem
             !     this%elem(i)%face(fdof,j) = this%elem(ii)%dof_plus
             this%elem(i)%face(fTdeg,j) = this%elem(ii)%Tdeg
             this%elem(i)%face(fTdof,j) = this%elem(ii)%Tdof

             !NEW FR
             this%elem(i)%face(fdof_plus,j) = this%elem(ii)%dof_plus

!             if(state%SP) then
!	        this%elem(i)%face(fdegP,j) = max(this%elem(i)%degP, this%elem(ii)%degP )
!                this%elem(i)%face(fdofP,j) = this%elem(ii)%dofP
!             endif

          else
             this%elem(i)%face(fdeg,j) = this%elem(i)%deg
             this%elem(i)%face(fdof,j) = this%elem(i)%dof

!             if(state%SP) then
!	        this%elem(i)%face(fdegP,j) = this%elem(i)%degP
!                this%elem(i)%face(fdofP,j) = this%elem(i)%dofP
!             endif

             !if(state%modelName == 'scalar' .or.state%modelName == '2eqs' .and. state%space%adapt%adapt_method == 'RTN') &   ! for dual problem
             !     this%elem(i)%face(fdof,j) = this%elem(i)%dof_plus
          endif
       enddo
    enddo

  end subroutine setEdgeQuadratureDegrees

 ! 3D plot a function goven by values in vertices
  subroutine plotVertexFunction3D(this, ifile, func )
    class(mesh), intent(inout) :: this
    integer, intent(in) :: ifile        ! number of the file
    real, dimension(1:this%npoin) :: func
    integer :: i, npoin

    npoin = this%npoin

    do i = 1, npoin
       write(ifile, *) this%x(i, 1:nbDim), func(i)
       write(ifile,'(x)')
    enddo
    write(ifile,'(x)')

  end subroutine plotVertexFunction3D

end module mesh_mod
