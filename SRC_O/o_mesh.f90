module mesh_mod
   use AMAdata
 ! use matrix_oper
  use main_data
 ! use f_mapping

  use element_mod
  !use mesh_oper
  use lapack_oper

  implicit none


  private

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
      procedure(init), deferred :: init

   end type Abstrmesh

   abstract interface
   subroutine allocElem( this , nelem )
      import :: Abstrmesh
      class ( Abstrmesh ), intent( inout ) :: this
      integer, intent( in ) :: nelem
   end subroutine allocElem

   subroutine init( this, nelem, npoin, nbelm, nbc, nbdim )
      import :: Abstrmesh
      class(Abstrmesh), intent(inout) :: this
      integer, intent(in) :: nelem, npoin, nbelm, nbc, nbdim
   end subroutine init
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
      real :: tol_min, tol_max             ! tolerances for mesh adaptation

      class(element), dimension(:), pointer :: elemP  ! for the solution of the local problems
      integer, allocatable,  dimension(:,:) :: elemP_idx  ! indexes of neighbouring  elements
      integer     :: elemP_idx_actual       ! index to actually considered element
      integer     :: flenP                  ! the maximal number of neighbouring triangles
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
      procedure :: seekNeighbours => seekNeighboursNew
!      procedure :: SeekCurvedBoundary => seekCurvedBoundaryNew
      procedure :: computeGeometry => computeGeometryNew

      procedure :: copy => copyMesh

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

   !> copies the whole structure of oldMesh to new mesh this
   subroutine copyMesh(this, oldMesh)
      class (mesh), intent(inout) :: this
      class (mesh), intent(in) :: oldMesh
      integer :: i

      print*, 'Copying meshes'

      ! allocate new elements and set basic settings, alloc b_edge, x, x_cur
      call this%init(oldMesh%nelem, oldMesh%npoin, oldMesh%nbelm, oldMesh%nbc, nbdim)
      this%curved_deg = oldMesh%curved_deg
      this%num_curv = oldMesh%num_curv
      this%xper(1:2,1:2) = oldMesh%xper(1:2,1:2)
      this%iper(1:2,1:2) = oldMesh%iper(1:2,1:2)
      this%periodic = oldMesh%periodic

      allocate( this%x( this%npoin, 1:nbDim) )
      allocate( this%xcur(1:this%npoin, 1:nbDim) )

      this%x( this%npoin, 1:nbDim) = oldMesh%x( this%npoin, 1:nbDim)
      this%xcur(1:this%npoin, 1:nbDim) = oldMesh%xcur(1:this%npoin, 1:nbDim)

      allocate( this%b_edge(1:this%nbelm) )
      do i = 1, this%nbelm
         call this%b_edge(i)%copy( oldMesh%b_edge(i) )
      end do !i

      this%diam = oldMesh%diam
      this%h = oldMesh%h
      this%domain_volume = oldMesh%domain_volume
      this%stop_adaptation = oldMesh%stop_adaptation
      this%adapt_level = oldMesh%adapt_level
      this%max_adapt_level = oldMesh%max_adapt_level
      this%adapt = oldMesh%adapt
      this%adaptR = oldMesh%adaptR
      this%tol_min = oldMesh%tol_min
      this%tol_max = oldMesh%tol_max

      ! elemP and other P variables are not copied now

      do i =1, this%nelem
         call this%elem(i)%copy( oldMesh%elem(i), state%time%disc_time )
      end do

   end subroutine copyMesh

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
    real :: tmp
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

    !allocate datas of the mesh
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
    enddo

!!!TODO not yet implemented
!! !      call elem%setFaces( lenloc, iloc )
!!

    close(ifile)


  end subroutine readmesh

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

  !> seeking of neighbouring of elements \f$ K, K',\ K\cap K'\not=\emptyset\f$
  subroutine seekNeighboursNew(this)
    class( mesh ), intent(inout) :: this
    integer, dimension(:,:), allocatable :: cyc
    integer, dimension(:), allocatable :: len_cyc
    integer, parameter :: max_cyc = 20
    integer :: i, j, k, k1, l, ni, nj, is1, is2, ip1, ip2, ii, ie, iie, jj, ib, ibound

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
          allocate(this%elem(i)%iBC(1:this%elem(i)%flen) )
          this%elem(i)%iBC(:) = 0

          allocate(this%elem(i)%tBC(1:this%elem(i)%flen) )
          this%elem(i)%tBC(:) = 0
       endif
    enddo

    deallocate(cyc,len_cyc)

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
    !!print*,' End of periodicity'


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

















!!!!!!!!!!!!!!!!!!  HG   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!§§
   subroutine allocElementsHG( this, nelem)
      class( MeshHG_t), intent(inout) :: this
      integer, intent(in) :: nelem

         allocate( ElementHG_t :: this%elem(1:nelem) )

   end subroutine allocElementsHG

!!!!!!!!!!!!!!!!!!!  RG   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!§§
   subroutine allocElementsRG( this, nelem)
      class( MeshRG_t), intent(inout) :: this
      integer, intent(in) :: nelem

      allocate( ElementRG_t :: this%elem(1:nelem) )

   end subroutine allocElementsRG


!!!!!!!!!!!!!!!!!!!  AMA   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!§§
   subroutine allocElementsAMA( this, nelem)
      class( MeshAMA_t), intent(inout) :: this
      integer, intent(in) :: nelem

      allocate( ElementAMA_t :: this%elem(1:nelem) )

   end subroutine allocElementsAMA





end module mesh_mod
