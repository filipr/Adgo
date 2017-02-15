module target_functional_mod
  use blocks_integ
  use lapack_oper
  use mesh_mod
  use mesh_oper
  use model_oper
  use paramets
  use st_interpol
  use stdgm_mod
!  use sort_mod

   implicit none

   type, public,abstract :: Abstr_Target_functional_t

   contains
      procedure(initTarFunc), deferred :: init
      procedure(findSupp), deferred :: findSupp
!      procedure(dualRHS_scalar), deferred :: dualRHS_scalar
!      final :: delete_Target_functional
   end type Abstr_Target_functional_t

   abstract interface
   subroutine initTarFunc( this, xy_coord, grid )
      import :: Abstr_Target_functional_t
      import :: nbDim
      import :: mesh
      class(  Abstr_Target_functional_t ), intent(inout) :: this
      real, dimension(1:nbDim), intent(in) :: xy_coord
      class( mesh ), intent(in) :: grid
   end subroutine initTarFunc

   subroutine findSupp( this, grid )
      import ::  Abstr_Target_functional_t
      import ::  mesh
      class(  Abstr_Target_functional_t ), intent(inout) :: this
      class( mesh ), intent(in) :: grid
   end subroutine findSupp

   end interface

   type, extends( Abstr_Target_functional_t ) :: Target_functional_t
      integer :: id ! id of the functional
      character(len=20) :: name ! name of the target quantity
      real :: Ju ! target quantity of the computed solution
      real :: Ju_exact
!      real, dimension(:,:), allocatable :: dJ_phi ! array - for linear J equals J(phi), for nonlinear J'[u](phi)
      logical :: linear ! target functional is linear -> no need of J'(u)(.)
      logical :: boundary ! target functional involves boundary integration
      logical :: time_dependent ! target functional involves time integration
      integer, allocatable, dimension(:) :: supp ! support of the target functional J
      integer :: isupp ! # of elements in support, i.e. size of the supp array
      real :: vol_supp ! volume of the support of the functional
      real, allocatable, dimension(:,:) :: xy_coord ! coordinates of the point value in tarFunc for id=3, or other coords
      integer :: coord_i ! the coordinate from 1:ndim we are interested in
      class( mesh ), allocatable :: grid ! only for point val
      real :: eps                              ! diameter for pointvalue of tarFunc for id=3



   contains
      procedure :: clean => cleanTargetFunctional ! nulify all parts which change after adaptation
      procedure :: init => init_tarFunc
      procedure :: initGrid => initGrid_tarFunc ! initializes the dual mesh - support of the functional
      procedure :: findSupp => findSupp_tarFunc
      procedure :: dualRHS_scalar => dualRHS_scalar_tarFunc
      procedure :: computeJu_exact
      procedure :: computeJu
      procedure :: delete => delete_Target_functional

   end type Target_functional_t

   type, extends( Target_functional_t ) :: Point_value_t
!     real :: eps                              ! diameter for pointvalue of tarFunc for id=3
!      real :: eps2 ! second epsilon - eps -- supp of J, eps2 - diameter of the circle around pointval

   contains

    procedure :: init => init_point_value
    procedure :: initGrid => initGrid_point_value  ! initializes the dual mesh - support of the functional
    procedure :: findSupp => findSupp_point_value
    procedure :: findSupp_point_value
    procedure :: dualRHS_scalar => dualRHS_scalar_pointVal
    procedure :: computeJu_exact => computeJu_exact_pointVal
    procedure :: computeJu => computeJu_pointVal


   end type Point_value_t

   !> integration of the solution over a subdomain \f$ \Omega_s \f$ of \f$ \Omega\f$ 
   type, extends( Target_functional_t ) :: U_over_subdomain_t

   contains
      procedure :: computeJu => computeJu_u_over_subdomain
      procedure :: computeJu_exact => computeJu_exact_u_over_subdomain
      procedure :: dualRHS_scalar => dualRHS_scalar_u_over_subdomain
      procedure :: findSupp => findSupp_u_over_subdomain
      procedure :: init => init_u_over_subdomain

   end type U_over_subdomain_t

      !> integration of the solution over a subdomain \f$ \Omega_s\f$  of \f$ \Omega\f$ 
   type, extends( Target_functional_t ) :: dudx_t
      integer :: dx ! 1:nbDim - derivative with respect to x (1) or y (2)

   contains
      procedure :: computeJu => computeJu_dudx
      procedure :: computeJu_exact => computeJu_exact_dudx
      procedure :: dualRHS_scalar => dualRHS_scalar_dudx
      procedure :: findSupp => findSupp_dudx
      procedure :: init => init_dudx

!      procedure :: setRHS => setRHS_dudx

   end type dudx_t


contains

   subroutine init_tarFunc( this, xy_coord, grid )
      class( Target_functional_t ), intent(inout) :: this
      real, dimension(1:nbDim), intent(in) :: xy_coord
      class( mesh ), intent(in) :: grid

      stop 'init_tarFunc should not be called - abstract'

   end subroutine init_tarFunc

   subroutine initGrid_tarFunc( this, gridfile )
      class( Target_functional_t ), intent(inout) :: this
      character(len=50), intent(in) :: gridfile     ! file with dual mesh

      stop 'initGrid_tarFunc should not be called - abstract'

   end subroutine initGrid_tarFunc

   subroutine findSupp_tarFunc( this , grid)
      class( Target_functional_t ), intent(inout) :: this
      class( mesh ), intent(in) :: grid

      stop 'findSupp_tarFunc should not be called - abstract'

   end subroutine findSupp_tarFunc

   subroutine cleanTargetFunctional( this )
      class( Target_functional_t ) :: this

      this%Ju = 0.0 ! target quantity of the computed solution
      this%Ju_exact = 0.0 ! ??
      if (allocated(this%supp)) &
         deallocate(this%supp)

      this%isupp = 0
      !this%vol_supp = ??
      !grid ! only for point val

   end subroutine cleanTargetFunctional

   subroutine delete_Target_functional( this )
    class( Target_functional_t ) :: this

    if (allocated( this%supp ) ) &
      deallocate( this%supp )
    if (allocated( this%xy_coord ) ) &
      deallocate( this%xy_coord )
    if (allocated( this%grid ))then
      deallocate( this%grid )
      print*, 'Control deallocation of DWR%grid.'
    endif

   end subroutine delete_Target_functional

   subroutine computeJu_exact( this, grid )
      class( Target_functional_t ), intent(inout) :: this
      class(mesh), intent(in) :: grid

      stop 'computeJu_exact should not be called - abstract'

   end subroutine computeJu_exact

   subroutine computeJu( this, grid )
      class( Target_functional_t ), intent(inout) :: this
      class(mesh), intent(in) :: grid

      stop 'computeJu should not be called - abstract'

   end subroutine computeJu

   !> compute dualRHS in a point, depends on the choice of J
   !> this functions does not watch for support - it is called only for \f$ x \in supp(J)\f$ 
   function dualRHS_scalar_tarFunc( this, x, t) result(f)
    class( Target_functional_t ),intent(in) :: this
    real, dimension(1:nbDim), intent(in) :: x
    real, intent(in) :: t
    real, dimension(1:ndim) :: f

    stop 'dualRHS_scalar_tarFunc - only abstract'

   end function dualRHS_scalar_tarFunc





!!! POINT VALUE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine init_point_value( this, xy_coord, grid)
      class( Point_value_t ), intent(inout) :: this
      real, dimension(1:nbDim), intent(in) :: xy_coord
      class( mesh ), intent(in) :: grid
      character(len=50) :: gridfile     ! file with dual mesh

      this%linear = .true.
      this%boundary = .false.
      this%time_dependent = .false.
      this%name = 'PointValue'
      this%eps = 0.1!0.01 * grid%diam

!      this%eps = 0.02 !max( 0.1 , 0.5 * h)
      print*, 'Epsilon =', this%eps
!      print*, 'h,eps:', h, this%eps
      allocate( this%xy_coord(1,1:nbDim) , source = 0.0 )
      this%xy_coord(1,1:2) = xy_coord(1:2)

      print*, 'The coord_i set to 1 for ndim>1'
      this%coord_i = 1

      gridfile = '../Grids/gridPointVal.grid'

      call this%initGrid(gridfile)
!      stop 'AFTER initGrid'

   end subroutine init_point_value


   subroutine initGrid_point_value( this, gridfile )
      class( Point_value_t ), intent(inout) :: this
      character(len=50), intent(in) :: gridfile     ! file with dual mesh
      integer :: nelem, npoin,nbelm
      real :: eps
      real, dimension(1:5,1:nbDim) :: x
      real, dimension(1:2,1:nbDim) :: xper
      integer, dimension(1:2,1:nbDim) :: iper
      integer, dimension(:), allocatable :: lenloc
      integer, dimension(:,:), allocatable :: iloc
      integer :: i

      allocate( mesh :: this%grid )

      print*, 'TODO : deinit target Grid in the end of computation!'
      ! read the template square
      call this%grid%read( gridfile )

!      print*, 'xy_coord and eps:', this%xy_coord(1,1:nbDim), this%eps
      call this%grid%shiftMesh( 1, this%xy_coord(1,1:nbDim), this%eps )
      call this%grid%seekNeighbours()

      ! CURVED ???
      this%grid%b_edge(1:this%grid%nbelm)%icurv = 0
      this%grid%elem(1:this%grid%nelem)%ibcur = -1
      this%grid%elem(1:this%grid%nelem)%jcur = 0
      this%grid%elem(1:this%grid%nelem)%deg_cur = 1

      call this%grid%computeGeometry()

      call this%grid%plot('../DWR/meshTargetJ.gnu')

      ! controls !

!      deallocate( lenloc, iloc )

   end subroutine initGrid_point_value

   !> marks the elements of the support of J, counts them and computes also volume of this support
   !> we also compute the value of this%eps from the parameters of the mesh
   subroutine findSupp_point_value( this, grid )
      class( Point_value_t ), intent(inout) :: this
      class( mesh ), intent(in) :: grid
      class( element ), pointer :: elem
      integer :: nelem, i, j, dof, Tdof, k
      integer, dimension(:,:), allocatable :: temp_supp
      real :: dist
      logical :: inSupp
      integer :: nsq, max_tri
      real, dimension(1:2,1:2) :: rmx
      integer, dimension(:, :, :), allocatable:: itri
      integer, dimension(:, :), allocatable :: ntri
      real, dimension(1:2) :: x_bary ! not used - only for the subroutine call

!      integer, dimension(1:3) :: A = (/ 1, 5, 2 /)
!      integer, dimension(4) :: B = (/ 1, -50, 2, 0 /)
!      integer, dimension(:), allocatable :: C
!
!
!      call heapsort(A)
!      write(*,*)'Sorted array :',A
!      C = MergeArrays(A,B)
!      print*, C



      type(intersect)			:: inter
      integer, dimension(:), allocatable	:: triangles! numbers of triangles in gridN which have
      integer :: NumTri

   end subroutine findSupp_point_value


      !> marks the elements of the support of J, counts them and computes also volume of this support
   subroutine findSupp_point_value_old( this, grid )
      class( Point_value_t ), intent(inout) :: this
      class( mesh ), intent(in) :: grid
      class( element ), pointer :: elem
      integer :: nsq, max_tri
      real, dimension(1:2,1:2) :: rmx
      integer, dimension(:, :, :), allocatable:: itri
      integer, dimension(:, :), allocatable :: ntri
      real, dimension(1:2) :: x_bary ! not used - only for the subroutine call

      print*, ' HERE findSupp_point_value_old!'

      ! grid, nsq = size of the square, frames of the domain, nelem, # of elems in the i,j square, indices of the elements in squares

      call FindTriangleIndexs(grid, nsq, rmx, max_tri, ntri, itri )

      call FindTriangleCoords( grid, elem, this%xy_coord(1,1:nbDim), x_bary(1:2), &
            nsq, max_tri, ntri, itri, rmx(1:2,1:2) )

      this%vol_supp = elem%area
      this%isupp = 1

      if (allocated(this%supp) ) &
         deallocate(this%supp)

      allocate( this%supp(1:this%isupp), source = 0 )
      this%supp(1:this%isupp) = elem%i

      ! => the point xy_coord is in the element elem !
      this%eps = elem%diam

   end subroutine findSupp_point_value_old

   !> computes the dual rhs in a point
   !> used when \f$  1/|B_{\epsilon}| \int_{B_{\epsilon}} u dx  \f$ 
   !> together with findSupp_point_value_old
   function dualRHS_scalar_pointVal( this, x, t) result(f)
    class( Point_value_t ),intent(in) :: this
    real, dimension(1:nbDim), intent(in) :: x
    real, intent(in) :: t
    real, dimension(1:ndim) :: f

    f(1:ndim) = 0.0

    if (t == state%time%FinTime) then
      f(this%coord_i) = 1.0 / (this%vol_supp)
    else
      stop 'dualRHS_scalar_pointVal: CONTROL nonzero only in final time'
    endif

   end function dualRHS_scalar_pointVal

   subroutine computeJu_exact_pointVal_old( this )
      class( Point_value_t ), intent(inout) :: this
!      class(mesh), intent(in) :: grid
      real, dimension(1:ndim) :: w

      print*, 'computeJu_exact_pointVal trule the value u(x,y) not integral!'

      if ( state%time_dependent ) then
         stop 'computeJu_exact_pointVal not implemented for time-dependent problems!'
      endif

      call Set_Model_Data( this%xy_coord(1,1:nbDim), state%time%FinTime , w, 1 )

      this%Ju_exact = w( this%coord_i )

   end subroutine computeJu_exact_pointVal_old

   !> computes the value J(u)
   subroutine computeJu_pointVal_old( this, grid )
      class( Point_value_t ), intent(inout) :: this
      class(mesh), intent(in) :: grid
      class( element ), pointer :: elem
      real, allocatable, dimension(:) :: wi
      real, allocatable, dimension(:,:) :: f
      integer :: i, j, Qdof, dof
      real :: local_Ju
      real :: time

      print*, 'computeJu_pointVal: Ju ~ on the choice of the pointVal approximation!!!'
      print*, 'F@R: computeJu_pointVal - set time ( temporarily finTime )'
      time = state%time%finTime

      this%Ju = 0.0

      do i = 1, this%isupp
         elem => grid%elem( this%supp(i) )
         Qdof = elem%Qdof

         allocate(wi(1:Qdof), source = 0.0)
         ! eval w in integ nodes, 0 - endpoint
         wi(1:Qdof) = Eval_whST_iCoord_Elem( elem, elem%TQnum, 0, this%coord_i)
         ! we use only wi(:,this%coord_i)
         !          call IntegrateVectorFunction2( elem, wi(1:Qdof,1:ndim), local_Ju)

         allocate( f(1:Qdof, 1:ndim) , source = 0.0 )
         do j = 1,Qdof
            f(j,1:ndim) = this%dualRHS_scalar(elem%xi(0,j,1:nbDim), time)
         end do !j

         local_Ju = EvalL2ScalarProduct( elem, f(1:Qdof, this%coord_i), wi(1:Qdof) )

         this%Ju = this%Ju + local_Ju

         deallocate(f, wi)

      end do

   end subroutine computeJu_pointVal_old

   ! computes the integral of the exact solution (one of its dimension from 1:ndim) over the elements of the support of J
   subroutine computeJu_exact_pointVal( this, grid )
      class( Point_value_t ), intent(inout) :: this
      class(mesh), intent(in) :: grid
      class( element ), pointer :: elem
      real :: w

      stop 'not done'
      print*, 'computeJu_exact_pointVal trule the value u(x,y) not integral!'
!
!      if ( state%time_dependent ) then
!         stop 'computeJu_exact_pointVal not implemented for time-dependent problems!'
!      endif
!
!      call Set_Model_Data( this%xy_coord, state%time%FinTime , w, 1 )
!
!      do i = 1, this%isupp
!         elem => grid%elem( this%supp(i) )
!
!         do j = 1, elem%Qdof
!            call Set_Model_Data( this%xy_coord, state%time%FinTime , w, 1 )
!             elem%xi(0, j, 1:nbDim)
!
!         end do !j
!
!      enddo
!      this%Ju_exact = w( this%coord_i )

   end subroutine computeJu_exact_pointVal

   !> computes the value J(u)
   subroutine computeJu_pointVal( this, grid )
      class( Point_value_t ), intent(inout) :: this
      class(mesh), intent(in) :: grid
      class( element ), pointer :: elem
      real, allocatable, dimension(:) :: wi
      real, allocatable, dimension(:,:) :: f
      integer :: i, j, Qdof, dof
      real :: local_Ju
      real :: time

      print*, 'computeJu_pointVal: Ju ~ on the choice of the pointVal approximation!!!'

      print*, 'F@R: computeJu_pointVal - set time ( temporarily finTime )'
      time = state%time%finTime

      this%Ju = 0.0

      do i = 1, this%isupp
         elem => grid%elem( this%supp(i) )
         Qdof = elem%Qdof

         allocate(wi(1:Qdof), source = 0.0)
         ! eval w in integ nodes, 0 - endpoint
         wi(1:Qdof) = Eval_whST_iCoord_Elem( elem, elem%TQnum, 0, this%coord_i)
         ! we use only wi(:,this%coord_i)
         !          call IntegrateVectorFunction2( elem, wi(1:Qdof,1:ndim), local_Ju)

         allocate( f(1:Qdof, 1:ndim) , source = 0.0 )
         do j = 1,Qdof
            f(j,1:ndim) = this%dualRHS_scalar(elem%xi(0,j,1:nbDim), time)
         end do !j

         local_Ju = EvalL2ScalarProduct( elem, f(1:Qdof, this%coord_i), wi(1:Qdof) )

         this%Ju = this%Ju + local_Ju

         deallocate(f, wi)

      end do

   end subroutine computeJu_pointVal












!!! U OVER SUBDOMAIN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   subroutine init_u_over_subdomain( this, xy_coord, grid)
      class( U_over_subdomain_t ), intent(inout) :: this
      real, dimension(1:nbDim), intent(in) :: xy_coord
      class( mesh ), intent(in) :: grid

      this%linear = .true.
      this%boundary = .false.
      this%time_dependent = .false.
      this%name = 'UoverSubdomain'

      if (ndim > 1) &
         print*, 'The coord_i set to 1 for ndim>1'
      this%coord_i = 1

!      print*, 'U over subdomain is now used for the average of the solution over the subdomain!'
!      print*, 'change in dualRHS_scalar_u_over_subdomain !!!'

!      allocate( this%xy_coord(1:2) , source = 0.0 )
!      this%xy_coord(1:2) = xy_coord(1:2)

   end subroutine init_u_over_subdomain

   !> marks the elements of the support of J, counts them and computes also volume of this support
   !> support = element which have iSubmesh == 1
   subroutine findSupp_u_over_subdomain( this, grid )
      class( U_over_subdomain_t ), intent(inout) :: this
      class( mesh ), intent(in) :: grid
      class( element ), pointer :: elem
      integer :: nelem, i, j
      integer, dimension(:), allocatable :: temp_supp
      real :: area


      if ( grid%elem(1)%iSubMesh == 0 ) then
         stop 'The subdomain partition was not done yet in findSupp_u_over_subdomain!'
      endif

      nelem = grid%nelem
      j = 0
      area = 0.0

      allocate( temp_supp(1:nelem), source = -1 )

      do i = 1, nelem
         elem => grid%elem(i)
         if (elem%iSubMesh == 1) then
            j = j+1
            temp_supp(j) = i
            area = area + elem%area
        end if
      end do !i

      this%isupp = j

      if (allocated(this%supp) ) &
         deallocate(this%supp)

      allocate( this%supp(1:this%isupp), source = 0 )
      this%supp(1:this%isupp) = temp_supp(1:this%isupp)

      deallocate( temp_supp )

      this%eps = elem%r_ins
      this%vol_supp = area

   end subroutine findSupp_u_over_subdomain

   function dualRHS_scalar_u_over_subdomain( this, x, t) result(f)
    class( U_over_subdomain_t ),intent(in) :: this
    real, dimension(1:nbDim), intent(in) :: x
    real, intent(in) :: t
    real, dimension(1:ndim) :: f

    f(1:ndim) = 0.0

!    print*, 'HERE dualRHS_scalar_u_over_subdomain'

    if (t == state%time%FinTime) then
!      f(1:ndim) = 1.0
      ! average value
      f( this%coord_i ) = 1.0 / this%vol_supp

    else
      stop 'dualRHS_scalar_u_over_subdomain: CONTROL nonzero only in final time'
    endif

   end function dualRHS_scalar_u_over_subdomain

   subroutine computeJu_exact_u_over_subdomain( this, grid )
      class( U_over_subdomain_t ), intent(inout) :: this
      class(mesh), intent(in) :: grid
      class( element ), pointer :: elem
      real, allocatable, dimension(:,:) :: wi
      real, allocatable, dimension(:,:) :: f
      integer :: i, j, Qdof, dof
      real :: local_Ju
      real :: time

      ! do we know the exact solution
      if ( state%model%known_sol ) then

         !      print*, 'F@R: computeJu_exact - set time ( temporarily finTime )'
         time = state%time%finTime

         this%Ju_exact = 0.0

         do i = 1, this%isupp
            elem => grid%elem( this%supp(i) )
            Qdof = elem%Qdof

            allocate(wi(1:Qdof,1:ndim), source = 0.0)
            allocate( f(1:Qdof, 1:ndim) , source = 0.0 )
            ! eval w in integ nodes, 0 - endpoint
            do j = 1, Qdof
               ! integ nodes, time, w, ityp = 1 (exact solution)
               ! we use only wi(:,this%coord_i)
               call Set_Model_Data( elem%xi(0,j,1:nbDim) , time , wi(j,1:ndim), 1 )

               f(j,1:ndim) = this%dualRHS_scalar(elem%xi(0,j,1:nbDim), time)
            enddo

            local_Ju = EvalL2ScalarProduct( elem, f(1:Qdof, this%coord_i), wi(1:Qdof, this%coord_i) )
            this%Ju_exact = this%Ju_exact + local_Ju

            deallocate(f, wi)

         end do !i
      ! convection problem isca 29
      else if ( state%model%icase == 29 ) then
         this%Ju_exact = 0.114
      ! peaks
      else if ( state%model%icase == 64 ) then
         ! connected with twoSquaresIn250
         this%Ju_exact = 4.6769935E-3
         ! connected with twoSquaresIn164
!         this%Ju_exact = 1.379437E-2
      ! unknown solution
      ! LL shaped peak
      else if ( state%model%icase == 65 ) then
         this%Ju_exact = 2.138 * 0.001
      ! CROSS domain problem, Ainsworth, Rankin 2012
      ! exact = 0.01630471454734821 *25 (25: area of the subdomain = 0.04)
      else if ( state%model%icase == 70 ) then ! domain cross
         this%Ju_exact = 0.40761786368370525 ! sent from prof Rankin
      else
         this%Ju_exact = 0.0
         print*, '# DWR: The exact solution is not a priori known for icase' , state%model%icase,'. Ju_exact cannot be computed!'
      endif

   end subroutine computeJu_exact_u_over_subdomain

   subroutine computeJu_u_over_subdomain( this, grid )
      class( U_over_subdomain_t ), intent(inout) :: this
      class(mesh), intent(in) :: grid
      class( element ), pointer :: elem
      real, allocatable, dimension(:) :: wi
      real, allocatable, dimension(:,:) :: f
      integer :: i, j, Qdof, dof
      real :: local_Ju
      real :: time

!      write(debug,*) 'F@R: computeJu_pointVal - set time ( temporarily finTime )'
      time = state%time%finTime

      this%Ju = 0.0

      do i = 1, this%isupp
         elem => grid%elem( this%supp(i) )
         Qdof = elem%Qdof

         allocate(wi(1:Qdof), source = 0.0)
         allocate( f(1:Qdof, 1:ndim) , source = 0.0 )
         ! eval w in integ nodes, 0 - endpoint
         wi(1:Qdof) = Eval_whST_iCoord_Elem( elem, elem%TQnum, 0, this%coord_i)
         ! we use only wi(:,this%coord_i)
         !   call IntegrateVectorFunction2( elem, wi(1:Qdof,1:ndim), local_Ju)
         do j = 1,Qdof
            f(j,1:ndim) = this%dualRHS_scalar(elem%xi(0,j,1:nbDim), time)
         end do !j

         local_Ju = EvalL2ScalarProduct( elem, f(1:Qdof, this%coord_i), wi(1:Qdof) )
         this%Ju = this%Ju + local_Ju

         deallocate(f, wi)

      end do

   end subroutine computeJu_u_over_subdomain








!!!! DUDX_T !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   subroutine init_dudx( this, xy_coord, grid)
      class( dudx_t ), intent(inout) :: this
      real, dimension(1:nbDim), intent(in) :: xy_coord
      class( mesh ), intent(in) :: grid

      this%linear = .true.
      this%boundary = .false.
      this%time_dependent = .false.
      this%name = 'dudx'
      ! du/dx -- 1
      this%dx = 1

      if (ndim > 1) &
         print*, 'The coord_i set to 1 for ndim>1'
      this%coord_i = 1

   end subroutine init_dudx

   !> marks the elements of the support of J, counts them and computes also volume of this support
   !> support = element which have iSubmesh == 1
   !> copied fro u_over_subdomain_t
   subroutine findSupp_dudx( this, grid )
      class( dudx_t ), intent(inout) :: this
      class( mesh ), intent(in) :: grid
      class( element ), pointer :: elem
      integer :: nelem, i, j
      integer, dimension(:), allocatable :: temp_supp
      real :: area

      if ( grid%elem(1)%iSubMesh == 0 ) then
         stop 'The subdomain partition was not done yet in findSupp_u_over_subdomain!'
      endif

      nelem = grid%nelem
      j = 0
      area = 0.0

      allocate( temp_supp(1:nelem), source = -1 )

      do i = 1, nelem
         elem => grid%elem(i)
         if (elem%iSubMesh == 1) then
            j = j+1
            temp_supp(j) = i
            area = area + elem%area
        end if
      end do !i

      this%isupp = j

      if (allocated(this%supp) ) &
         deallocate(this%supp)

      allocate( this%supp(1:this%isupp), source = 0 )
      this%supp(1:this%isupp) = temp_supp(1:this%isupp)

      deallocate( temp_supp )

      this%eps = elem%r_ins
      this%vol_supp = area

   end subroutine findSupp_dudx

   function dualRHS_scalar_dudx( this, x, t) result(f)
    class( dudx_t ),intent(in) :: this
    real, dimension(1:nbDim), intent(in) :: x
    real, intent(in) :: t
    real, dimension(1:ndim) :: f

    f(1:ndim) = 0.0

    if (t == state%time%FinTime) then
!      f(1:ndim) = 1.0
      ! average value
      ! ok but must be multiplied by dw/dx
      f( this%coord_i ) =  - 1.0 !/ this%vol_supp

    else
      stop 'dualRHS_scalar_u_over_subdomain: CONTROL nonzero only in final time'
    endif


   end function dualRHS_scalar_dudx

   subroutine computeJu_exact_dudx( this, grid )
      class( dudx_t ), intent(inout) :: this
      class(mesh), intent(in) :: grid
      class( element ), pointer :: elem
      real, allocatable, dimension(:,:) :: wi
      real, allocatable, dimension(:,:) :: f
      integer :: i, j, Qdof, dof
      real :: local_Ju
      real :: time

      ! do we know the exact solution
      if ( state%model%known_sol ) then

         time = state%time%finTime

         this%Ju_exact = 0.0

         do i = 1, this%isupp
            elem => grid%elem( this%supp(i) )
            Qdof = elem%Qdof

            allocate(wi(1:Qdof,1:ndim), source = 0.0)
            allocate( f(1:Qdof, 1:ndim) , source = 0.0 )
            ! eval w in integ nodes, 0 - endpoint
            do j = 1, Qdof
               ! integ nodes, time, w, ityp = 1 (exact solution)
               ! we use only wi(:,this%coord_i) , 8 - derivative du/dx
               call Set_Model_Data( elem%xi(0,j,1:nbDim) , time , wi(j,1:ndim), 8 )

               f(j,1:ndim) = this%dualRHS_scalar(elem%xi(0,j,1:nbDim), time)
            enddo

            local_Ju = EvalL2ScalarProduct( elem, f(1:Qdof, this%coord_i), wi(1:Qdof, this%coord_i) )
            this%Ju_exact = this%Ju_exact + local_Ju

            deallocate(f, wi)

         end do !i
      ! convection problem isca 68
      else if ( state%model%icase == 68 ) then
         ! connected with submesh square8TwoTriangles
         this%Ju_exact = 1.585090814E-3
      else
         this%Ju_exact = 0.0
         print*, '# DWR: The exact solution is not a priori known for icase' , state%model%icase,'. Ju_exact cannot be computed!'
      endif

   end subroutine computeJu_exact_dudx

   subroutine computeJu_dudx( this, grid )
      class( dudx_t ), intent(inout) :: this
      class(mesh), intent(in) :: grid
      class( element ), pointer :: elem
      real, allocatable, dimension(:,:) :: wi
      real, allocatable, dimension(:,:) :: f
      integer :: i, j, Qdof, dof
      real :: local_Ju
      real :: time

      if (ndim > 1) &
         stop 'Problem in computeJu_dudx wi( 1 ,1:Qdof) with ndim>1!'

!      if ( state%dual ) &
!         stop 'the primal solution is not saved in elem%wST'

      write(debug,*) 'F@R: computeJu_pointVal - set time ( temporarily finTime )'
      time = state%time%finTime

      this%Ju = 0.0

      do i = 1, this%isupp
         elem => grid%elem( this%supp(i) )
         Qdof = elem%Qdof

         allocate(wi(1:ndim, 1:Qdof), source = 0.0)
         allocate( f(1:Qdof, 1:ndim) , source = 0.0 )

         ! eval dw/dx in integ nodes, 0 - endpoint
         wi( 1:ndim ,1:Qdof) = evalSTfunInIntTime_spaceDer( elem, ndim, elem%dof, elem%Tdof, &
            elem%wST(1:ndim,1:elem%dof,1:elem%Tdof), 0, elem%tQnum, this%dx )

         ! we use only wi(:,this%coord_i)
         !          call IntegrateVectorFunction2( elem, wi(1:Qdof,1:ndim), local_Ju)
         do j = 1,Qdof
            f(j,1:ndim) = this%dualRHS_scalar(elem%xi(0,j,1:nbDim), time)
         end do !j

         local_Ju = EvalL2ScalarProduct( elem, f(1:Qdof, this%coord_i), &
            wi(this%coord_i, 1:Qdof) )

         this%Ju = this%Ju + local_Ju

         deallocate(f, wi)

      end do

   end subroutine computeJu_dudx

end module target_functional_mod
