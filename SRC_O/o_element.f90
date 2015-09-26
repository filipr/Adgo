module element_mod
   use lapack_oper
   use paramets
   use geometry

   use elemental_estims_mod


   implicit none




  !> values of Jacobian and inverse in integ. nodes (volume, edges)
  type, public :: Der_Fmapp
     integer :: Qdof
     real, allocatable, dimension(:) :: JF        ! Jacobian of F in integ nodes
     real, allocatable, dimension(:,:,:) :: D1F   ! inverse of Jac. matrix in integ nodes

   contains

      procedure, public :: copy => copyDer_Fmapp

  end type Der_Fmapp

  !> mapping of reference triangle or square on element elem \f$ F:\hat{K}\to K\f$
  type, public :: Fmapping
     integer :: deg           ! degree of polynomial mapping
     integer :: dof           ! degree of fredoms
     logical :: iFlin         !  true for linear mapping F (triangles or paralelogram)
     real :: JF0                              ! const Jacobian of F
     real, allocatable, dimension(:,:) :: D1F0    ! const inverse of Jac. matrix
     real, allocatable, dimension(:,:) :: F       ! Langrangian nodes of element K
     type(Der_Fmapp)  :: V                         ! JF and D1F only for curved: volume
     type(Der_Fmapp), allocatable, dimension(:) :: E   ! JF and D1F only for curved: edges (1:elem%flen)

  contains

    ! procedure, public :: init => initFmapping
     procedure, public :: copy => copyFmapping

  end type Fmapping

    !> structures for red green (de)refinement
  type, public :: RGref
     integer :: subel   ! number of daughters elements
     integer, allocatable, dimension(:) :: daughter ! links to refined sub-elements
  end type RGref


!   type, public, abstract :: AbstrElement_t
!
!
!   contains
!      procedure(init), deferred :: init
!
!   end type AbstrElement_t
!
!   abstract interface
!   subroutine init( this, i, type, flen, iloc)
!      import :: AbstrElement_t
!      class( AbstrElement_t) :: this
!      integer, intent(in) :: i
!      integer, intent(in) :: type
!      integer, intent(in) :: flen
!      integer, dimension(:), intent(in) :: iloc
!   end subroutine init
!   end interface

   !> basic element type, when adaptation needed, use some of its descendants
   type, public :: element
      integer :: i
      integer :: type                           ! =3 for triang, 4 for quad
      integer :: flen                           ! number of nodes including hanging
      integer, allocatable, dimension(:,:) :: face ! indexes of vertices, neighs,....
      real, allocatable, dimension(:) :: xc      ! barycentre of the element
      real, allocatable, dimension(:,:, :) :: xi    ! integ nodes 0:flen,
      integer :: isupp                        ! number  of elements having a common vertex with elem
      integer, allocatable, dimension(:,:) :: supp  ! list of elements having a common vertex with elem
      real :: area                              ! area of an element
      real :: diam                              ! diameter of an element
      real :: r_ins                             ! radius of inscribed circle
      real :: limit_par                         ! parameter for limiting

      integer :: ibcur ! global index of the curved edge
      integer :: jcur  ! = 0 if polygonal, > 0 curved ib, index of the curved edge in the element
      !
      integer :: deg_cur
      integer :: dof_cur
      type(Fmapping), pointer :: F              ! mapping: K^ -> K
      real, pointer, dimension(:,:) :: n        ! vector of outer normal - 1:flen, 1:nbDim
      real, pointer, dimension(:) :: dn         ! size of vector of outer normal
      real, pointer, dimension(:,:) :: nc ! outer normals on curved edge in integ nodes - 1:Qdof, 1:nbDim
      real, pointer, dimension(:) :: dnc  ! dS "increase" on curved edge in integ nodes
      integer :: deg                   ! degree of polynomial approximation
      integer :: degP			! degree of polynomial approx. for pressure
      integer :: dofP			! DOF per component for l (degree for pressure)
      integer :: dof, dof_plus, dof_fix ! DOF per one components for p and p+1
      logical :: deg_plus              ! if true, compute fluxes for p+1, otherwise
      integer :: Tdeg, Tdof, Tdof_plus     ! degree of polynomial approximation in time
      integer::  TQnum                 ! index of time quadrature, number of quadrature nodes
      integer :: Qnum                  ! index of implicit volume quadrature
      integer :: Qdeg                  ! degree of implicit volume quadrature
      integer :: Qdof                  ! # of nodes of ^^^^^^^^^^^^^^

      type(Mblock), pointer  :: Mass                ! matrix block, mass matrix
      type(Mblock), pointer  :: MassInv             ! matrix block, inverse of mass matrix
      type(Mblock), pointer  :: Stiff               ! matrix block, stiff matrix
      type(Mblock), pointer, dimension (:) :: block ! matrix blocks (0:rlen)
      type(Mblock), pointer, dimension (:) :: blockST ! matrix blocks (0:rlen)
      type(Mblock), pointer, dimension (:) :: ILU   ! preconditioner matrix blocks (0:rlen)
      type(Mblock), pointer, dimension (:,:) :: SPblock ! matrix blocks (0:rlen)
      type(Mblock), pointer, dimension (:) :: XXX   ! preconditioner matrix blocks (0:rlen)
      real, pointer, dimension(:, :) :: wc               ! computational approximation (solution) vector satisfying stop. criteria based on AEE
      real, pointer, dimension(:,:) :: w               ! solution vector
      real, pointer, dimension(:,:,:) :: wSP               ! solution vector for SP - 1:3 (iV1,iV2,iP), 1:2 (time levels we need to remember), 1:max(elem%dof,elem%dofP)
      real, pointer, dimension(:,:,:) :: wST             ! solution vector for ST DGM - 1:ndim, 1:elem%dof, 1:elem%Tdof
      real, pointer, dimension(:,:) :: wSTfin       ! solution vector for ST DGM at final time in space integ nodes
      real, pointer, dimension(:,:) :: wSTfinAD       ! solution vector for ST DGM at final time - coefficients
      real, pointer, dimension(:,:,:) :: rhsST             ! RHS  ST DGM - 1:ndim, 1:elem%dof, 1:elem%Tdof
      real, pointer, dimension(:,:) :: wS              ! temporary used array
      real, pointer, dimension(:,:,:) :: wSS            ! temporary used array
      real, pointer, dimension(:,:) :: vec             ! vectors: right hand sides, ...
      real, pointer, dimension(:,:) :: res_func    !residual function made up from residual vector
      real, pointer, dimension(:,:,:) :: MGw           ! solution vector for MG


      integer, allocatable, dimension(:) :: iBC          ! index of boundary component, if any
      integer, allocatable, dimension(:) :: tBC          ! type of boundary condition, if any
      integer       ::  ncv                          ! sequence of element's dof
      integer       ::  ncvP                          ! sequence of element's dofP

      integer :: hsplit                             ! type of h-ref of h-deref
      integer :: psplit                             ! type of p-ref
      logical :: to_recompute                  ! recompute the solution


      real :: CK, Cb, CTn, CKo, Cbo
      real :: CP                                !Poincare constant


      integer :: HO_deg                         ! degree of the HO reconstruction
      real :: max_eigenvals
      real, dimension(:,:), allocatable :: lifting  ! lifting operator

      ! COMMENTED BEFORE
        !  type ( Elemental_estims_t ), allocatable :: estims

      real, pointer, dimension(:,:,:) :: RTNflux     ! RTN fluxes with time history
      integer :: RTNflux_deg
      real, pointer, dimension(:,:) :: RTNphi     ! RTN fluxes with time history
      real :: errTot, errDisc, errAlg                ! total error u - u_h^i, discretization error u - u_h, algebraic error u_h - u_h^i on the element
      real :: estTot, estDisc, estAlg, estRem        ! estimate of total error u - u_h^i, discretization error u - u_h, algebraic error u_h - u_h^i on the element
      real, allocatable, dimension(:) :: estimFNCD   ! Discretization component of Flux NonConformity estimator
      real, allocatable, dimension(:) :: estimFNCA   ! Algebraic component of Flux NonConformity estimator
      real, allocatable, dimension(:) :: estimFNCT   ! Total Flux NonConformity estimator
      integer :: MGdeg                               ! MG level for element
      integer :: MGdof                               ! DOFs for current MG level
      integer :: MGncv                               ! sequence of element's dof

      real :: estimA, estimS, estimT, estimST, estim_loc ! estimates of "residal errors"
      real :: jumpsJh                           ! \sum_{G \in \partial K} |G|^{-1} \int_G [u_h]^2 dS
      real :: rezid, rezid1                     ! reziduum for artificial viscosity
      real :: reg, reg1, reg2, regT0, regT1, regT2          ! regularity of the solution
      !real, pointer, dimension(:)  :: TEST
      real :: errL2, errH1, errL8,interL8, interLq, interH1
      real, dimension(1:2, 1:3) :: errSnorm	  ! STDGM L2/H1 norms in three different nodes 1-endpoint, 2-node of the time quad, 3-1/2


      real, allocatable, dimension(:,:) :: eta    ! array of different error indicators

      logical :: RTNfluxeval                    !whether or not RTNflux for elem at state%time%iter is already evaluated
      integer :: per                           ! > 0 periodic neighbours element

      !HG ###############xxxxxxxxx###
      ! NEEDED FOR OTHER TYPES OF ADAPT TOO
      logical :: HGnode                         ! HG structure is allocated and defined
      integer, allocatable, dimension(:) :: HGvertex ! for elements with HG nodes: vertices
      integer, allocatable, dimension(:,:) :: HGface ! for elements with HG nodes: faces

      !RG ##########################
      integer :: RGhistory                           ! levels of history  of RG refinement
      integer :: RGlevel                             ! level of refinement
      type(RGref), pointer, dimension(:) :: RG     ! data for red green (de)refinement - used in HG TOO !!!

      integer :: RGreclen                        ! number of elements in RGrecomputes
      integer, dimension(1:4) :: RGrecomput      ! (old) index of elems for recomputation
      character :: RGtype                       ! =N (none) =R (red), =G (green) refinement
      integer :: RGindex                        ! index within level ref (see type(RGref)

      logical :: dealloc                        ! elem (of old mesh) will be deallocated


      !AMA ##########################
      real    :: ama_p                               ! change of deg
      real, dimension(:),  allocatable :: rgabc ! metric for AMA adaptation
      real :: amaL                             !  ama ... size of the dominant eigenvalue
      real :: amaR                             !ama < 1 ... ratio between the small and the big eigenvals
      real, dimension(1:2) :: amaD             ! = ( cos(alpha), sin(alpha) ) orientation of AMA ellipse




!      type(FormulElem) :: formul

!      type(Aposteriori_Est) :: apostEst
!      type(Apriori_Est) :: aprEst
!      type(MultigridElem) :: multigrid


     contains

     procedure :: allocBlocks
     procedure :: allocVectors
     procedure :: init => initElem
     procedure :: computeElementGeometry => computeElementGeometry2D
     procedure :: computeF => ComputeF2D
     procedure :: computeDF => ComputeDF2D
     procedure :: copy => copyElement
     procedure :: copyBlocks
     procedure :: initElementDof
     procedure :: setF => SetF2D

!     procedure :: copyVectors
   !   procedure :: writeElem
   !   procedure :: setFaces

   !   generic, public :: assignment(=) => copyElem

   end type element

      ! DO NOT USE UNTIL THE PROBLEM WITH ALLOCATION IS SOLVED IN THE GCC COMPILER (bugzilla Bug 65359)
    type, public, extends( element ) :: ElementHG_t
!      logical :: HGnode                         ! HG structure is allocated and defined
!      integer, allocatable, dimension(:) :: HGvertex ! for elements with HG nodes: vertices
!      integer, allocatable, dimension(:,:) :: HGface ! for elements with HG nodes: faces

   contains
      procedure:: init => initElemHG
      procedure :: setHG ! set HG element
   end type ElementHG_t

   type, public, extends( element ) :: ElementAMA_t
!      real    :: ama_p                               ! change of deg
!      real, dimension(:),  allocatable :: rgabc ! metric for AMA adaptation
!      real :: amaL                             !  ama ... size of the dominant eigenvalue
!      real :: amaR                             !ama < 1 ... ratio between the small and the big eigenvals
!      real, dimension(1:2) :: amaD             ! = ( cos(alpha), sin(alpha) ) orientation of AMA ellipse
      contains

      procedure:: init => initElemAMA

   end type ElementAMA_t


   type, public, extends( element ) :: ElementRG_t
!     integer :: RGhistory                           ! levels of history  of RG refinement
!     integer :: RGlevel                             ! level of refinement
!     !type(RGref), pointer, dimension(:) :: RG     ! data for red green (de)refinement -> structure on following 2 lines
!     integer :: subel   ! number of daughters elements
!     integer, allocatable, dimension(:) :: daughter ! links to refined sub-elements
!
!     integer :: RGreclen                        ! number of elements in RGrecomputes
!     integer, dimension(1:4) :: RGrecomput      ! (old) index of elems for recomputation
!     character :: RGtype                       ! =N (none) =R (red), =G (green) refinement
!     integer :: RGindex                        ! index within level ref (see type(RGref)

     contains
     procedure:: init => initElemRG

   end type ElementRG_t

   contains

   subroutine copyFmapping( this, old_map, flen)
      class( Fmapping ), intent(inout) :: this
      class( Fmapping ), intent(in) :: old_map
      integer, intent(in), optional :: flen
      integer :: i

      this%deg = old_map%deg
      this%dof = old_map%dof

      allocate( this%F(1:this%dof,1:nbDim))
      this%F(1:this%dof, 1:nbDim) = old_map%F(1:this%dof, 1:nbDim)

      this%iFlin = old_map%iFlin

      !linear mapping - constant Jacobian
      if (this%iFlin) then
         this%JF0 = old_map%JF0
         allocate(this%D1F0(1:nbDim,1:nbDim) )
         this%D1F0(1:nbDim,1:nbDim) = old_map%D1F0(1:nbDim, 1:nbDim)
      !nonlinear mapping - curved edges
      ! element with NON constant Jacobian
      ! Jacobian and D1F in Gauss integ. nodes of  curved element
      else
         if (.not. present(flen)) &
            stop 'copyFmapping for elements with nonlinear mapping should be called with elem%flen.'

         print*, 'copyFmapping: The curved edges are NOT tested!'
         allocate( this%E(1:flen) )

         call this%V%copy( old_map%V )
         do i = 1, flen
            call this%E(i)%copy( old_map%E(i) )
         enddo
      endif

   end subroutine copyFmapping

   !> copy Der_Fmapp
   subroutine copyDer_Fmapp( this, old_Der )
      class(Der_Fmapp), intent(inout) :: this
      class(Der_Fmapp), intent(in) :: old_Der

      this%Qdof = old_Der%Qdof
      allocate( this%D1F(1:this%Qdof, 1:nbDim, 1:nbDim) )
      this%D1F(1:this%Qdof, 1:nbDim, 1:nbDim) = &
         old_Der%D1F(1:this%Qdof, 1:nbDim, 1:nbDim)
      allocate( this%JF(1:this%Qdof) )
      this%JF(1:this%Qdof) = old_Der%JF(1:this%Qdof)

   end subroutine copyDer_Fmapp

   !> initialization of basic elem data - noAdapt
   subroutine initElem( this, i, type, flen, iloc )
      class( element) :: this
      integer, intent(in) :: i
      integer, intent(in) :: type
      integer, intent(in) :: flen
      integer, dimension(:), intent(in) :: iloc

      !stop 'abstract initElem this should not be called'

      this%i = i
      this%type = type
      this%flen = flen

      this%HGnode = .false.

      allocate( this%face(1:max_face, 1:flen) )
      this%face(idx,1:flen) = iloc(1:flen)

   end subroutine initElem

   !> passing pointers and values from "type(element) :: old_elem" to
   subroutine copyElement( new_elem , old_elem, disc_time )
      class( element ), intent(in) :: old_elem
      class( element ), intent(out) :: new_elem
      character(len=20), intent(in) :: disc_time

      logical :: rset
      integer :: flen, Qdof, Fdof, j, nsl, nsl1


   !> passing pointers and values from "type(element) :: old_elem" to
  !> "type(element) :: new_elem"
!  subroutine PassElem2Elem(new_elem, old_elem )


    new_elem%i = old_elem%i
    new_elem%flen = old_elem%flen
    flen = new_elem%flen
    new_elem%type = old_elem%type

    new_elem%deg = old_elem%deg
    new_elem%degP = old_elem%degP
    new_elem%Tdeg = old_elem%Tdeg
    call new_elem%initElementDof()

    new_elem%Qnum = old_elem%Qnum
    new_elem%Qdeg = old_elem%Qdeg
    new_elem%Qdof = old_elem%Qdof
    new_elem%TQnum = old_elem%TQnum

    new_elem%CP = old_elem%CP   ! Poincare constant

    ! copy FACES
    allocate( new_elem%face(1:max_face, 1:flen) )
    new_elem%face(1:max_face, 1:flen) = old_elem%face(1:max_face, 1:flen)

    allocate( new_elem%xc(1:nbDim) )
    new_elem%xc = old_elem%xc


    Fdof = max(new_elem%Qdof, maxval( new_elem%face(fGdof, :)  ) )
    allocate(new_elem%xi(0:flen, 1:Fdof, 1:nbDim ) )
    new_elem%xi(0:flen, 1:Fdof, 1:nbDim ) = &
      old_elem%xi(0:flen, 1:Fdof, 1:nbDim )
    new_elem%isupp = old_elem%isupp
    allocate( new_elem%supp(1:new_elem%isupp, 1:2) ) ! second argument = 1 or 2 - share a face?
    new_elem%supp(1:new_elem%isupp, 1:2) = old_elem%supp(1:new_elem%isupp, 1:2)
    new_elem%area  = old_elem%area
    new_elem%diam  = old_elem%diam
    new_elem%r_ins  = old_elem%r_ins
    new_elem%limit_par  = old_elem%limit_par

    new_elem%ibcur = old_elem%ibcur
    new_elem%jcur = old_elem%jcur
    new_elem%deg_cur = old_elem%deg_cur
    new_elem%dof_cur = old_elem%dof_cur
    ! copy the Fmapping
    call new_elem%F%copy( old_elem%F, flen)

    allocate (new_elem%n(1:flen, 1:nbDim), new_elem%dn(1:flen) )
    new_elem%n(1:flen, 1:nbDim)  = old_elem%n(1:flen, 1:nbDim)
    new_elem%dn(1:flen)  = old_elem%dn(1:flen)

    !curved edge
    if (new_elem%ibcur > 0) then
      Qdof = new_elem%face(fGdof, new_elem%jcur )
      allocate( new_elem%nc(1:Qdof, 1:nbDim), new_elem%dnc(1:Qdof) )
      new_elem%nc(1:Qdof, 1:nbDim)  = old_elem%nc(1:Qdof, 1:nbDim)
      new_elem%dnc(1:Qdof)  = old_elem%dnc(1:Qdof)
    endif


    if ( associated(old_elem%RTNflux) ) &
      print*, 'RTNflux and res_func is not copied in elem%copy!'

    ! TODO allocate estimators
!    allocate(new_elem%eta(1:max_eta, 1:ndim) )
!    allocate(new_elem%estimFNCD(1:ndim) )
!    allocate(new_elem%estimFNCA(1:ndim) )
!    allocate(new_elem%estimFNCT(1:ndim) )

    !TODO:
    print*, 'etas, rset..., Are not copied in elem%copyElement!'
!    ! the corresponding element was changed??
!    rset = .false.
!    if(new_elem%type /= old_elem%type .or. new_elem%flen /= old_elem%flen) then
!       rset = .true.
!    else
!       do j=1,new_elem%flen
!          if( new_elem%face(fdeg,j) /= old_elem%face(fdeg,j) ) rset = .true.
!       enddo
!    endif
!
!    if(rset) then
!       call SetElementQuadraturesDegrees(new_elem )
!       print*,'Adaptation.f90, rset variable is .true. !!!'
!    else

!    endif



    !allocate element solution arrays and matrix blocks:
    ! TODO: p_plus ??? FERROR
    call new_elem%allocVectors( disc_time )
    call new_elem%copyBlocks( old_elem, disc_time )


   end subroutine copyElement

   !> allocate elem% blocks and solution vectors
   subroutine allocBlocks( this, disc_time)
      class( element ), intent(inout) :: this
      character(len=20), intent(in) :: disc_time
      integer :: Tdof, dof, flen, nsl, nslST,i , j, nsl1, nslST1

      dof = this%dof
      flen = this%flen
      nsl = dof * ndim

      allocate( this%Mass, this%MassInv, this%Stiff )
      call InitMblock( this%Mass, dof, dof )
      call InitMblock( this%MassInv, dof, dof )
      call InitMblock( this%Stiff, dof, dof )

      allocate( this%ILU(0:flen), this%block(0:flen) )
      call InitMblock(this%block(0), nsl, nsl )

!      allocate( this%XXX(), this%SPblock() )

      if( disc_time /= 'STDG' ) then
         call InitMblock(this%ILU(0), nsl, nsl )
         do i=1,nsl
             this%ILU(0)%Mb(i,i) = 1.
          enddo

         do j=1,flen
             if(this%face(neigh,j) > 0) then
                nsl1 = this%face(fdof,j) * ndim
                call InitMblock(this%block(j), nsl, nsl1)
                call InitMblock(this%ILU(j), nsl, nsl1)
             endif
         enddo

      else !  STDGM
         nslST = dof * this%Tdof * ndim

         call InitMblock(this%ILU(0), nslST, nslST )
         do i=1,nslST
             this%ILU(0)%Mb(i,i) = 1.
         enddo
         allocate( this%blockST( 0:flen) )
         call InitMblock(this%blockST(0), nslST, nslST)

         ! off diagonal blocks
         do j=1,flen
            if(this%face(neigh,j) > 0) then
                nsl1 = this%face(fdof,j) *  ndim
                nslST1 = nsl1 * this%face(fTdof,j)
                call InitMblock(this%block(j), nsl, nsl1)
                call InitMblock(this%blockST(j), nslST, nslST1)
                call InitMblock(this%ILU(j), nslST, nslST1)
             endif
         enddo
      endif

   end subroutine allocBlocks

   !> allocate and copy elem% blocks from old_elem
   subroutine copyBlocks( this , old_elem, disc_time)
      class( element ), intent(inout) :: this
      class( element ), intent(in) :: old_elem
      character(len=20), intent(in) :: disc_time
      integer :: Tdof, dof, flen, nsl, nslST, j

      dof = this%dof
      flen = this%flen
      nsl = dof * ndim

      allocate( this%Mass, this%MassInv, this%Stiff )
      call this%Mass%CopyMblock(old_elem%Mass)
      call this%MassInv%CopyMblock( old_elem%MassInv )
      call this%Stiff%CopyMblock( old_elem%Stiff )

      allocate( this%ILU(0:flen), this%block(0:flen) )
      call this%block(0)%CopyMblock( old_elem%block(0) )

      call this%ILU(0)%CopyMblock( old_elem%ILU(0) )

      if( disc_time /= 'STDG' ) then
         do j=1,flen
             if(this%face(neigh,j) > 0) then
                call this%block(j)%CopyMblock( old_elem%block(j) )
                call this%ILU(j)%CopyMblock( old_elem%ILU(j) )
             endif
         enddo
      else !  STDGM

         allocate( this%blockST( 0:flen) )
         call this%blockST(0)%CopyMblock( old_elem%blockST(0) )
         ! off diagonal blocks
         do j=1,flen
            if(this%face(neigh,j) > 0) then
                call this%block(j)%CopyMblock( old_elem%block(j) )
                call this%ILU(j)%CopyMblock( old_elem%ILU(j) )
                call this%blockST(j)%CopyMblock( old_elem%blockST(j) )
             endif
         enddo
      endif

   end subroutine copyBlocks

   !> allocate elem%  solution vectors
   subroutine allocVectors(this, disc_time)
      class( element ), intent(inout) :: this
      character(len=20), intent(in) :: disc_time
      integer :: Tdof, dof, flen, nsl, nslST

      dof = this%dof
      flen = this%flen
      nsl = dof * ndim
      allocate( this%w(0:this%Tdeg+1, 1:ndim*this%dof), source = 0.0 )

      if ( disc_time /= 'STDG' ) then
         allocate( this%vec(1:max_vecs, 1:this%dof_plus * ndim), source = 0.0 )
      else !  STDGM
         nslST = dof * this%Tdof * ndim
         ! TODO, FERROR: shouldnt be used Tdof_plus ???
         allocate( this%vec(1:max_vecs, 1 : this%dof_plus * ndim *this%Tdof), source = 0.0)
         allocate( this%wST(1:ndim, 1:dof, 1:this%Tdof), source = 0.0 )
         allocate( this%rhsST(1:ndim, 1:dof, 1:this%Tdof), source = 0.0 )
         allocate( this%wSTfin(1:ndim, 1:this%Qdof), source = 0.0 )
         allocate( this%wSTfinAD(1:ndim, 1:dof), source = 0.0 )
      endif

   end subroutine allocVectors



    !> computing of geometrical properties of each element \f$ K\in{\cal T}_h\f$,
  !> the mapping \f$ F:\hat{K}\to K,\ K\in{\cal T}_h\f$
  subroutine computeElementGeometry2D( this, grid_x, x_inn)
    class (element), intent(inout) :: this  ! elem = element
    real, dimension(:,:), intent(in) :: grid_x
    real, dimension(:,:), intent(in), optional :: x_inn
    real :: diam, sum_inner_edges, rlen
    real, allocatable, dimension (:,:) :: x, x0
    !real, allocatable, dimension (:,:) :: xi, Fxi
    integer:: flen, Fdeg, Fdof
    integer:: j, k, k1, k0, ib, l, m, tmpint
    integer:: l1, l2, l3, l4
    real:: t1, t2
    !    3D structures
    integer, dimension(1:4) :: nodes   ! nodes of face
    real, dimension(1:3) :: vec1, vec2


   !  print*, 'ComputeElementGeometry2D'


    flen = this%flen               ! number of faces of elem
    Fdeg = this%deg_cur           ! degree of curved elemenet
    if(Fdeg <= 1 ) Fdeg =  1      ! polygonal element Fdeg = 1

    allocate(this%F)
    this%F%deg = Fdeg


    if (allocated(this%xc) ) then
      print*, 'elem(',this%i,')%xc: already allocated'
    else
      allocate( this%xc(1:nbDim) )
    endif


    !  print*, 'elem(',this%i,')%xc: ', size(this%xc)
    diam = 0.
    sum_inner_edges = 0.

    ! computing of outer normals
    allocate( this%n(1:flen,1:nbDim))
    allocate( this%dn(1:flen))

    do k=1,flen
          k1 = mod(k,flen) + 1
          this%n(k,1) = grid_x(this%face(idx,k1),2) - grid_x(this%face(idx,k),2)

          this%n(k,2) = grid_x(this%face(idx,k),1)  - grid_x(this%face(idx,k1),1)

          this%dn(k) = dot_product(this%n(k,:),this%n(k,:))**0.5
    enddo

    ! computing of mappings F: K^ -> K
    if(this%type == 3) then  ! triangles
       Fdof = (Fdeg + 1)*(Fdeg+2)/2
       this%F%dof = Fdof

       !print*,'#####',this%i,ib,j
       ! print*, 'HG'
       ! if(j>0 .and. this%HGnode) j = this%HGface(1, j)

       allocate(this%F%F(1:Fdof,1:nbDim))
       allocate(x(1:Fdof,1:nbDim))

!       if(this%HGnode) then  ! element with hanging node
!          do k=1,3
!             x(k,1:nbDim) = grid_x(this%face(idx,elem%HGvertex(k) ),1:nbDim)
!          enddo
!       else !element without hanging node
          do k=1,3
             x(k,1:nbDim) = grid_x(this%face(idx,k),1:nbDim)
          enddo
!       endif

       this%area = Area(x(1,:), x(2,:), x(3,:) )
       ! evaluation of a "barycentre" of elem
       do l=1,2
          this%xc(l) = sum(x(1:3,l) )/3
       enddo

       if ( present(x_inn) ) then
         ib = this%ibcur     ! numbed of curved edge
         j = this%jcur       ! index of element edge adjacent to curved edge

         if(Fdeg == 2) then
          ! piecewise quadratic approximation of the boundary, one node is inserted
          ! computing of edge centres
          do l=1,3
             x(l+3, 1:nbDim) = (x(l,1:nbDim) + x(mod(l,3)+1,1:nbDim) )/2
          enddo

          ! a node on a curved edge is replaced by an inserted node
!          x(j+3,1:nbDim) = grid%b_edge(ib)%x_inn(1,1:nbDim)
          x(j+3,1:nbDim) = x_inn(1, 1:nbDim)
         elseif(Fdeg == 3) then

          ! piecewise cubic approximation of the boundary
          ! two nodes are inserted

          ! computing of nodes on  edges
          do l=1,3
             x(2*l+3-1, 1:nbDim) = (2.*x(l,1:nbDim) + 1.*x(mod(l,3)+1,1:nbDim) )/3
             x(2*l+3  , 1:nbDim) = (1.*x(l,1:nbDim) + 2.*x(mod(l,3)+1,1:nbDim) )/3
          enddo
          x(10, 1:nbDim) = (x(1,1:nbDim) + x(2,1:nbDim) + x(3,1:nbDim) )/3

          ! two nodes on a curved edge is replaced by an inserted node
          x(2*j+3-1,1:nbDim) = x_inn(1,1:nbDim)
          x(2*j+3  ,1:nbDim) = x_inn(2,1:nbDim)

          ! correction for central node, 1st possibility
          !x(10,1) = sum(x(1:9,1))/9.
          !x(10,2) = sum(x(1:9,2))/9.

          ! correction for central node, 2nd possibility
          l1 = 2*(j-1)+1+3
          l2 = 2*(mod(j,3))+2+3
          l3 = 2*(j-1)+2+3
          l = mod(j,3) + 1
          l4 = 2*(mod(l,3))+1+3
          x(10,1:nbDim) = ( x(l1,1:nbDim) + x(l2,1:nbDim) + x(l3,1:nbDim) + x(l4,1:nbDim))/4.

         elseif(Fdeg > 3) then
          print*,'P_k, k> 3 approximation of boundary is not implemeneted'
          stop
         endif
       else
         if (Fdeg > 1) stop 'Problem in ComputeElementGeometry - Curved edge!'
       endif !present




!       ! quadrilateral elements
!    elseif(this%type == 4) then
!       Fdof = (Fdeg+1)**2
!       this%F%dof = Fdof
!       ib = this%ibcur     ! numbed of curved edge
!       j = this%jcur       ! index of element edge adjacent to curved edge
!       if(j>0 .and. this%HGnode) j = this%HGface(1, j)
!
!       allocate(this%F%F(1:Fdof,1:nbDim))
!       allocate(x0(1:4,1:nbDim))
!       allocate(x(1:Fdof,1:nbDim))
!
!       if(this%HGnode) then  ! element with hanging node
!          do k=1,4
!             x0(k,1:nbDim) = grid_x(elem%face(idx,elem%HGvertex(k) ),1:nbDim)
!          enddo
!       else !element without hanging node
!          do k=1,4
!             x0(k,1:nbDim) = grid_x(elem%face(idx,k),1:nbDim)
!          enddo
!       endif
!       this%area = Area(x0(1,:), x0(2,:), x0(3,:) ) + Area(x0(3,:), x0(4,:), x0(1,:) )
!       ! evaluation of a "barycentre" of elem
!       do l=1,2
!          this%xc(l) = sum(x0(1:4,l) )/4
!       enddo
!
!       do k0=1,Fdeg+1
!          do k1=1,Fdeg+1
!             k = (k1-1)*(Fdeg+1) + k0
!             t1 = 1.*(k0-1)/Fdeg
!             t2 = 1.*(k1-1)/Fdeg
!             !x(k,1:nbDim) = grid_x(elem%face(idx,1), 1:nbDim)  &
!             !     + (grid_x(elem%face(idx,2), 1:nbDim) - grid_x(elem%face(idx,1), 1:nbDim)) * t1  &
!             !     + (grid_x(elem%face(idx,4), 1:nbDim) - grid_x(elem%face(idx,1), 1:nbDim)) * t2  &
!             !     + (grid_x(elem%face(idx,3), 1:nbDim) - grid_x(elem%face(idx,2), 1:nbDim) &
!             !     -  grid_x(elem%face(idx,4), 1:nbDim) + grid_x(elem%face(idx,1), 1:nbDim)) * t1 * t2
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


    else
       print*, 'Qudarilatterals are not implemented in o_element ComputeGeometry2D'
       print*,'Sorry, only triangles and quadrilateralls are implemented'
          stop
    endif

    call this%setF( Fdof, x(:,:) )

    !write(*,*) 'F1:',elem%F%F(:,1)
    !write(*,*) 'F2:',elem%F%F(:,2)
!     write(*,*) 'F3:',elem%F%F(:,3)

    ! computing the sum of length of all inner edges
    diam = 0.
    do k =1,this%type-1
       do l=k+1,this%type
          rlen= dot_product(x(l,1:nbDim) - x(k,1:nbDim),x(l,1:nbDim) - x(k,1:nbDim))**0.5
          diam = max(diam, rlen)
       enddo
    enddo

    do k=1,flen
       rlen = VectorNorm(this%n(k,:) )
       diam = max(diam, rlen)
       !if(this%face(neigh,k) > 0) ! sum of all edges (innner and boundary)
       sum_inner_edges = sum_inner_edges + rlen
    enddo
    this%diam = diam
      ! FR this has to be set after the subroutine call
!    this%h = max(this%h, diam)

    if(sum_inner_edges <= 0.) sum_inner_edges = this%diam  ! only for singular case
    this%r_ins = 2*this%area/sum_inner_edges

    ! parameter for limiting
    !this%limit_par = 1./(this%area**0.75)/this%diam
    !this%limit_par = 1./(this%area**0.75)/sum_inner_edges
    !2D
    this%limit_par = 1./this%area/sum_inner_edges


    deallocate( x )

  end subroutine ComputeElementGeometry2D



 !> setting of mappings \f$ F: \hat{K}\to K\quad \forall K\in{\cal T}_h\f$, i.e.
  !> mapping of the reference element onto the physical one
  !>
  !> \f$ K= \f$ elem,  given by Langrangian nodes within  TRIANGLE or QUADRILATERALL,
  !> \f$ F \f$ can be linear, quadratic, cubic, ....,
  !> F(l,1:nbDim) means  the \f$(x,y)\f$ coordinates of the l\f$^{\rm th}\f$
  !> Langrangian nodes
  subroutine SetF2D( this , nod, x)
    class( element ), intent(inout) :: this   ! elem = element
    integer, intent(in) :: nod
    real, dimension(1:nod, 1:nbDim), intent(in) :: x
    integer :: deg
    real, dimension(:, :), pointer :: F

    F => this%F%F(1:this%F%dof, 1:nbDim)

    if(nod /= this%F%dof) then
       print*,'Inconsistency in SetF:',nod, this%F%dof
       stop
    endif

    this%F%iFlin = .false.

    deg = this%F%deg

    if (nbDim .ne. 2) stop 'SetF2D in 3D!'

    if(this%type == 3) then   !  TRIANGLES
       if(deg == 1) then
          F(1, 1:2) = x(1,1:2)
          F(2, 1:2) = x(2,1:2) - x(1,1:2)
          F(3, 1:2) = x(3,1:2) - x(1,1:2)

          this%F%iFlin = .true.   ! linear mapping

       elseif(deg == 2) then
          F(1, 1:2) = x(1,1:2)
          F(2, 1:2) = -(x(2,1:2) -4*x(4,1:2) +3*F(1,1:2))
          F(3, 1:2) = -(x(3,1:2) -4*x(6,1:2) +3*F(1,1:2))
          F(4, 1:2) = 2*(x(2,1:2) + F(1,1:2) - 2*x(4,1:2))
          F(6, 1:2) = 2*(x(3,1:2) + F(1,1:2) - 2*x(6,1:2))
          F(5, 1:2) = 4*(x(5,1:2) - F(1,1:2) -(F(2,1:2)+F(3,1:2))/2  &
               - (F(4,1:2)+F(6,1:2))/4 )
       elseif(deg == 3) then
          F(1,1:2) = x(1,1:2)
          F(2,1:2) = -11./2.*x(1,1:2)+x(2,1:2)+9*x(4,1:2)-9./2.*x(5,1:2)
          F(3,1:2) = -11./2.*x(1,1:2)+x(3,1:2)-9./2.*x(8,1:2)+9*x(9,1:2)
          F(4,1:2) = 9*x(1,1:2)-9./2.*x(2,1:2)-45./2.*x(4,1:2)+18*x(5,1:2)
          F(5,1:2) = 18*x(1,1:2)-45./2.*x(4,1:2)+9./2.*x(5,1:2)-9./2.*x(6,1:2) &
               -9./2.*x(7,1:2)+9./2.*x(8,1:2)-45./2.*x(9,1:2)+27*x(10,1:2)
          F(6,1:2) = 9*x(1,1:2)-9./2.*x(3,1:2)+18*x(8,1:2)-45./2.*x(9,1:2)
          F(7,1:2) = -9./2.*x(1,1:2)+9./2.*x(2,1:2)+27./2.*x(4,1:2)-27./2.*x(5,1:2)
          F(8,1:2) = -27./2.*x(1,1:2)+27*x(4,1:2)-27./2.*x(5,1:2)+27./2.*x(6,1:2) &
               +27./2.*x(9,1:2)-27*x(10,1:2)
          F(9,1:2) = -27./2.*x(1,1:2)+27./2.*x(4,1:2)+27./2.*x(7,1:2)-27./2.*x(8,1:2)&
               +27*x(9,1:2)-27*x(10,1:2)
          F(10,1:2) = -9./2.*x(1,1:2)+9./2.*x(3,1:2)-27./2.*x(8,1:2)+27./2.*x(9,1:2)
       else
          print*,'Maximal mapping P_3 on triangles implemented'
          stop
       endif

    elseif(this%type == 4 ) then   !  QUADRILATERALS

       if(deg == 1) then
          F(1, 1:2) = x(1,1:2)
          F(2, 1:2) = x(2,1:2) - x(1,1:2)
          F(3, 1:2) = x(3,1:2) - x(1,1:2)
          F(4, 1:2) = x(4,1:2) - x(2,1:2)  - x(3,1:2) + x(1,1:2)

          if(dot_product(F(4,1:2), F(4,1:2)) <= 1D-14) this%F%iFlin = .true.   ! linear mapping
       elseif(deg == 2) then
          F(1, 1:2) = x(1,1:2)
          F(2, 1:2) = -3*x(1,1:2)+4*x(2,1:2)-x(3,1:2)
          F(3, 1:2) = -3*x(1,1:2)+4*x(4,1:2)-x(7,1:2)
          F(4, 1:2) = 9*x(1,1:2)-12*x(2,1:2)+3*x(3,1:2)-12*x(4,1:2) &
               +16*x(5,1:2)-4*x(6,1:2)+3*x(7,1:2)-4*x(8,1:2)+x(9,1:2)
          F(5, 1:2) = 2*x(1,1:2)-4*x(2,1:2)+2*x(3,1:2)
          F(6, 1:2) = -6*x(1,1:2)+12*x(2,1:2)-6*x(3,1:2)+8*x(4,1:2) &
               -16*x(5,1:2)+8*x(6,1:2)-2*x(7,1:2)+4*x(8,1:2)-2*x(9,1:2)
          F(7, 1:2) = 4*x(1,1:2)-8*x(2,1:2)+4*x(3,1:2)-8*x(4,1:2)+ &
               16*x(5,1:2)-8*x(6,1:2)+4*x(7,1:2)-8*x(8,1:2)+4*x(9,1:2)
          F(8, 1:2) = 2*x(1,1:2)-4*x(4,1:2)+2*x(7,1:2)
          F(9, 1:2) = -6*x(1,1:2)+8*x(2,1:2)-2*x(3,1:2)+12*x(4,1:2) &
               -16*x(5,1:2)+4*x(6,1:2)-6*x(7,1:2)+8*x(8,1:2)-2*x(9,1:2)

       elseif(deg == 3) then
          F(1,1:2) = x(1,1:2)
          F(2,1:2) = -5.5*x(1,1:2) +9*x(2,1:2) -4.5*x(3,1:2) +x(4,1:2)
          F(3,1:2) = -5.5*x(1,1:2) +9*x(5,1:2) -4.5*x(9,1:2) +x(13,1:2)
          F(4,1:2) = 121./4.*x(1,1:2) -49.5*x(2,1:2) +99./4.*x(3,1:2)-5.5*x(4,1:2) &
               -49.5*x(5,1:2) +81*x(6,1:2) -40.5*x(7,1:2) +9*x(8,1:2) &
               +99./4.*x(9,1:2) -40.5*x(10,1:2) +81./4.*x(11,1:2) -4.5*x(12,1:2) &
               -5.5*x(13,1:2) +9*x(14,1:2) -4.5*x(15,1:2) +x(16,1:2)
          F(5,1:2) = 9*x(1,1:2) -22.5*x(2,1:2) +18*x(3,1:2) -4.5*x(4,1:2)
          F(6,1:2) = -49.5*x(1,1:2) +495./4.*x(2,1:2) -99*x(3,1:2) +99./4.*x(4,1:2) &
               +81*x(5,1:2) -202.5*x(6,1:2) +162*x(7,1:2) -40.5*x(8,1:2)  &
               -40.5*x(9,1:2) +405./4.*x(10,1:2) -81*x(11,1:2) +81./4.*x(12,1:2) &
               +9*x(13,1:2) -22.5*x(14,1:2) +18*x(15,1:2) -4.5*x(16,1:2)
          F(7,1:2) = 81*x(1,1:2) -202.5*x(2,1:2) +162*x(3,1:2) -81./2.*x(4,1:2) &
               -405./2.*x(5,1:2) +2025./4.*x(6,1:2) -405*x(7,1:2) +405./4.*x(8,1:2) &
               +162*x(9,1:2) -405*x(10,1:2) +324*x(11,1:2) -81*x(12,1:2) &
               -81./2.*x(13,1:2) +405./4.*x(14,1:2) -81*x(15,1:2) +81./4.*x(16,1:2)
          F(8,1:2) = 9*x(1,1:2) -45./2.*x(5,1:2) +18*x(9,1:2) -4.5*x(13,1:2)
          F(9,1:2) = -99./2.*x(1,1:2) +81*x(2,1:2) -81./2.*x(3,1:2) +9*x(4,1:2) &
               +495./4.*x(5,1:2) -405./2.*x(6,1:2) +405./4.*x(7,1:2) -45./2.*x(8,1:2) &
               -99*x(9,1:2) +162*x(10,1:2) -81*x(11,1:2) +18*x(12,1:2) &
               +99./4.*x(13,1:2) -81./2.*x(14,1:2) +81./4.*x(15,1:2) -4.5*x(16,1:2)
          F(10,1:2) = -4.5*x(1,1:2) +27./2.*x(2,1:2) -27./2.*x(3,1:2) &
               +9./2.*x(4,1:2)
          F(11,1:2) = 99./4.*x(1,1:2) -297./4.*x(2,1:2) +297./4.*x(3,1:2) &
               -99./4.*x(4,1:2) -81./2.*x(5,1:2) +243./2.*x(6,1:2) -243./2.*x(7,1:2) &
               +81./2.*x(8,1:2) +81./4.*x(9,1:2) -243./4.*x(10,1:2) +243./4.*x(11,1:2) &
               -81./4.*x(12,1:2) -4.5*x(13,1:2) +27./2.*x(14,1:2) -27./2.*x(15,1:2) &
               +9./2.*x(16,1:2)
          F(12,1:2) = -81./2.*x(1,1:2) +243./2.*x(2,1:2) -243./2.*x(3,1:2) &
               +81./2*x(4,1:2) +405./4.*x(5,1:2) -1215./4.*x(6,1:2) +1215./4.*x(7,1:2) &
               -405./4.*x(8,1:2) -81*x(9,1:2) +243*x(10,1:2) -243*x(11,1:2) &
               +81*x(12,1:2) +81./4.*x(13,1:2) -243./4.*x(14,1:2) +243./4.*x(15,1:2) &
               -81./4.*x(16,1:2)
          F(13,1:2) = 81./4.*x(1,1:2) -243./4.*x(2,1:2) +243./4.*x(3,1:2) &
               -81./4.*x(4,1:2) -243./4.*x(5,1:2) +729./4.*x(6,1:2) -729./4.*x(7,1:2)&
               +243./4.*x(8,1:2)+ 243./4.*x(9,1:2)-729./4. *x(10,1:2) &
               +729./4.*x(11,1:2) -243./4.*x(12,1:2) -81./4.*x(13,1:2) &
               +243./4.*x(14,1:2) -243./4.*x(15,1:2) +81./4.*x(16,1:2)
          F(14,1:2) = -4.5*x(1,1:2) +27./2.*x(5,1:2) -27./2.*x(9,1:2)+ 9./2.*x(13,1:2)
          F(15,1:2) = 99./4.*x(1,1:2) -81./2.*x(2,1:2) +81./4.*x(3,1:2) &
               -4.5*x(4,1:2) -297./4.*x(5,1:2) +243./2.*x(6,1:2) -243./4.*x(7,1:2) &
               +27./2.*x(8,1:2)+297./4.*x(9,1:2) -243./2.*x(10,1:2) +243./4.*x(11,1:2)&
               -27./2.*x(12,1:2) -99./4.*x(13,1:2) +81./2.*x(14,1:2) -81./4.*x(15,1:2)&
               +9./2.*x(16,1:2)
          F(16,1:2) = -81./2.*x(1,1:2) +405./4.*x(2,1:2)- 81*x(3,1:2) +81./4.*x(4,1:2)&
               +243./2.*x(5,1:2) -1215./4.*x(6,1:2) +243*x(7,1:2) -243./4.*x(8,1:2)  &
               -243./2.*x(9,1:2)+ 1215./4.*x(10,1:2) -243*x(11,1:2) +243./4.*x(12,1:2) &
               +81./2.*x(13,1:2) -405./4.*x(14,1:2) +81*x(15,1:2) -81./4.*x(16,1:2)
       else
          print*,'Maximal mapping P_3 on quadrilateralls implemented'
          stop
       endif
    else
       print*,'Only triangles or quadrilaterals implemented in SetF'
       stop
    endif

  end subroutine SetF2D

 !> evaluation of \f$ F(x_i)\in R^2,\ x_i\in \hat{K},\ i=1..nod,  \f$
  !>
  !> xi= \f$x_i\in \hat{K}\f$ are arbitrary nodes within reference elements,
  !> \f$ F \f$ mapping of reference element \f$\hat{K}\f$ to the actual one \f$K\f$=elem,
  !> F can be linear, quadratic, cubic, ....
  subroutine ComputeF2D( this, nod, xi, Fx )
    class( element), intent(in):: this   ! elem = element
    integer, intent(in) :: nod            ! number of nodes
    real, dimension(1:nod, 1:nbDim), intent(in) :: xi
    real, dimension(1:nod, 1:nbDim), intent(out) :: Fx
    integer :: i, deg

    real, dimension(:, :), pointer :: F

    F => this%F%F(1:this%F%dof, 1:nbDim)
    deg = this%F%deg

    if(this%type == 3) then  ! TRIANGLES
          do i=1, 2
             ! absolute term
             Fx(1:nod, i) = F(1, i)

             if(deg >= 1) then
                ! linear terms
                Fx(1:nod, i) = Fx(1:nod, i)  &
                     + F(2, i)*xi(1:nod, 1) + F(3, i)*xi(1:nod, 2)
                if(deg >= 2) then
                   ! quadratic terms
                   Fx(1:nod, i) = Fx(1:nod, i)  &
                        + F(4, i) * xi(1:nod, 1) * xi(1:nod, 1) &
                        + F(5, i) * xi(1:nod, 1) * xi(1:nod, 2) &
                        + F(6, i) * xi(1:nod, 2) * xi(1:nod, 2)
                   if(deg >= 3) then
                      ! cubic terms
                      Fx(1:nod, i) = Fx(1:nod, i)  &
                           + F(7, i) * xi(1:nod, 1) * xi(1:nod, 1) * xi(1:nod, 1) &
                           + F(8, i) * xi(1:nod, 1) * xi(1:nod, 1) * xi(1:nod, 2) &
                           + F(9, i) * xi(1:nod, 1) * xi(1:nod, 2) * xi(1:nod, 2) &
                           + F(10,i) * xi(1:nod, 2) * xi(1:nod, 2) * xi(1:nod, 2)
                      if(deg >= 4) then
                         ! 4-th order terms
                         Fx(1:nod, i) = Fx(1:nod, i)  &
                              + F(11,i)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1) &
                              + F(12,i)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2) &
                              + F(13,i)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2) &
                              + F(14,i)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,2) &
                              + F(15,i)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,2)
                         if(deg >=5) then
                            print*,'Only polynoms up to order 4 are implemented',&
                                 'in computeF for triangles'
                            stop
                         endif
                      endif
                   endif
                endif
             endif
          enddo

    elseif(this%type == 4) then   ! QUADRILATERALS
          do i=1, 2
             ! absolute term
             Fx(1:nod, i) = F(1, i)

             if(deg >= 1) then
                ! bilinear terms
                Fx(1:nod, i) = Fx(1:nod, i) + F(2, i)*xi(1:nod, 1) &
                     + F(3, i)*xi(1:nod, 2) + F(4, i)*xi(1:nod, 1)*xi(1:nod, 2)
                if(deg >= 2) then
                   ! biquadratic terms
                   Fx(1:nod, i) = Fx(1:nod, i)  &
                        + F(5, i)*xi(1:nod, 1)*xi(1:nod, 1) &
                        + F(6, i)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 2) &
                        + F(7, i)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 2)*xi(1:nod, 2) &
                        + F(8, i)*xi(1:nod, 2)*xi(1:nod, 2) &
                        + F(9, i)*xi(1:nod, 1)*xi(1:nod, 2)*xi(1:nod, 2)
                   if(deg >= 3) then
                      ! bicubic terms
                      Fx(1:nod, i) = Fx(1:nod, i)  &
                           + F(10, i)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 1) &
                           + F(11, i)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 2) &
                           + F(12, i)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 2) &
                           *xi(1:nod, 2) &
                           + F(13, i)*xi(1:nod, 1)**3 *xi(1:nod, 2)**3 &
                           + F(14, i)*xi(1:nod, 2)*xi(1:nod, 2)*xi(1:nod, 2) &
                           + F(15, i)*xi(1:nod, 1)*xi(1:nod, 2)*xi(1:nod, 2)*xi(1:nod, 2) &
                           + F(16, i)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 2)*xi(1:nod, 2) &
                           *xi(1:nod, 2)
                      !print*,'# Nontested  Eval_QPolynoms for deg = 3 !!'
                      if(deg >=4) then
                         print*,'Only Qpolynoms up to order 3 are implemented in ComputeF_4'
                         stop
                      endif
                   endif
                endif
             endif
          enddo

    else
          print*,'Only triagles and quarilaterals are implemented in ComputeF'
          stop
    endif


  end subroutine ComputeF2D


  !> evaluation of \f$ \frac{D F(x_i)}{D \hat{x}} \in R^{2\times 2},
  !> \ x_i\in \hat{K},\ i=1..nod,  \f$
  !>
  !> xi= \f$x_i\in \hat{K}\f$ are arbitrary nodes within reference elements,
  !> \f$ F \f$ mapping of reference element \f$\hat{K}\f$ to the actual one \f$K\f$=elem,
  !> F can be linear, quadratic, cubic, ....
  subroutine ComputeDF2D( this, nod, xi, DF  )
    class( element ), intent(in):: this   ! elem = element
    integer, intent(in) ::  nod
    real, dimension(1:nod, 1:nbDim), intent(in) :: xi
    real, dimension(1:nod, 1:nbDim, 1:nbDim ), intent(out) :: DF
    integer :: deg, i

    real, dimension(:, :), pointer :: F

    F => this%F%F(1:this%F%dof, 1:nbDim)

    deg = this%F%deg

    !if(this%i == 15 ) then
    !   do i=1,nod
    !      write(99,*) xi(i, 1:2)
    !   enddo
    !endif

    if(this%type == 3) then  ! TRIANGLES
          if(deg >= 1) then
             ! linnear mapping
             DF(1:nod,1,1) = F(2,1)
             DF(1:nod,1,2) = F(3,1)
             DF(1:nod,2,1) = F(2,2)
             DF(1:nod,2,2) = F(3,2)

             if(deg >= 2) then
                ! quadratic mapping
                DF(1:nod,1,1) = DF(1:nod,1,1) + 2*F(4,1)*xi(1:nod,1)+ F(5,1)*xi(1:nod,2)
                DF(1:nod,1,2) = DF(1:nod,1,2) + 2*F(6,1)*xi(1:nod,2)+ F(5,1)*xi(1:nod,1)
                DF(1:nod,2,1) = DF(1:nod,2,1) + 2*F(4,2)*xi(1:nod,1)+ F(5,2)*xi(1:nod,2)
                DF(1:nod,2,2) = DF(1:nod,2,2) + 2*F(6,2)*xi(1:nod,2)+ F(5,2)*xi(1:nod,1)

                if(deg >= 3) then
                   ! cubic mapping
                   DF(1:nod,1,1) = DF(1:nod,1,1) + 3*F(7,1)*xi(1:nod,1)*xi(1:nod,1) &
                        + 2*F(8,1)*xi(1:nod,1)*xi(1:nod,2) &
                        + F(9,1)*xi(1:nod,2)*xi(1:nod,2)
                   DF(1:nod,1,2) = DF(1:nod,1,2) + 3*F(10,1)*xi(1:nod,2)*xi(1:nod,2) &
                        + 2*F(9,1)*xi(1:nod,1)*xi(1:nod,2) &
                        + F(8,1)*xi(1:nod,1)*xi(1:nod,1)
                   DF(1:nod,2,1) = DF(1:nod,2,1) + 3*F(7,2)*xi(1:nod,1)*xi(1:nod,1) &
                        + 2*F(8,2)*xi(1:nod,1)*xi(1:nod,2) &
                        + F(9,2)*xi(1:nod,2)*xi(1:nod,2)
                   DF(1:nod,2,2) = DF(1:nod,2,2) + 3*F(10,2)*xi(1:nod,2)*xi(1:nod,2) &
                        + 2*F(9,2)*xi(1:nod,1)*xi(1:nod,2) &
                        + F(8,2)*xi(1:nod,1)*xi(1:nod,1)
                   !      print*,'P_3 approximation of boundary is implemented for', &
                   !                  ' Jacobian but not tested !!!'

                   if(deg >=4) then
                      print*,'ComputeDF_3 for deg > 3 is not implemented'
                      stop
                   endif
                endif
             endif
          endif

    elseif(this%type == 4) then   ! QUADRILATERALS

          if(deg >= 1) then
             ! bi-linnear mapping
             DF(1:nod,1,1) = F(4,1)*xi(1:nod,2) + F(2,1)
             DF(1:nod,1,2) = F(4,1)*xi(1:nod,1) + F(3,1)

             DF(1:nod,2,1) = F(4,2)*xi(1:nod,2) + F(2,2)
             DF(1:nod,2,2) = F(4,2)*xi(1:nod,1) + F(3,2)

             if(deg >= 2) then
                ! bi-quadratic mapping
                DF(1:nod,1,1) = DF(1:nod,1,1) + &
                     2*F(5,1)*xi(1:nod,1)+ &
                     2*F(6,1)*xi(1:nod,1)*xi(1:nod,2) + &
                     2*F(7,1)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2) + &
                     F(9,1)*xi(1:nod,2)*xi(1:nod,2)

                DF(1:nod,1,2) = DF(1:nod,1,2) + &
                     F(6,1)*xi(1:nod,1)*xi(1:nod,1)+ &
                     2*F(7,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2) + &
                     2*F(8,1)*xi(1:nod,2) + &
                     2*F(9,1)*xi(1:nod,1)*xi(1:nod,2)

                DF(1:nod,2,1) = DF(1:nod,2,1) + &
                     2*F(5,2)*xi(1:nod,1)+ &
                     2*F(6,2)*xi(1:nod,1)*xi(1:nod,2) + &
                     2*F(7,2)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2) + &
                     F(9,2)*xi(1:nod,2)*xi(1:nod,2)

                DF(1:nod,2,2) = DF(1:nod,2,2) + &
                     F(6,2)*xi(1:nod,1)*xi(1:nod,1)+ &
                     2*F(7,2)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2) + &
                     2*F(8,2)*xi(1:nod,2) + &
                     2*F(9,2)*xi(1:nod,1)*xi(1:nod,2)

                if(deg >= 3) then
                   ! bi-cubic mapping
                   DF(1:nod,1,1) = DF(1:nod,1,1) + &
                        3*F(10,1)*xi(1:nod,1)*xi(1:nod,1)+ &
                        3*F(11,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2) + &
                        3*F(12,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2) + &
                        3*F(13,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2) &
                        *xi(1:nod,2)  + &
                        F(15,1)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,2) + &
                        2*F(16,1)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,2)


                   DF(1:nod,1,2) = DF(1:nod,1,2) + &
                        F(11,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)+ &
                        2*F(12,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2) + &
                        3*F(13,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2)&
                        *xi(1:nod,2)+ &
                        3*F(14,1)*xi(1:nod,2)*xi(1:nod,2) + &
                        3*F(15,1)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,1) + &
                        3*F(16,1)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,1)*xi(1:nod,1)

                   DF(1:nod,2,1) = DF(1:nod,2,1) + &
                        3*F(10,2)*xi(1:nod,1)*xi(1:nod,1)+ &
                        3*F(11,2)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2) + &
                        3*F(12,2)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2) + &
                        3*F(13,2)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2) &
                        *xi(1:nod,2)  + &
                        F(15,2)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,2) + &
                        2*F(16,2)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,2)


                   DF(1:nod,2,2) = DF(1:nod,2,2) + &
                        F(11,2)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)+ &
                        2*F(12,2)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2) + &
                        3*F(13,2)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2)&
                        *xi(1:nod,2)+ &
                        3*F(14,2)*xi(1:nod,2)*xi(1:nod,2) + &
                        3*F(15,2)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,1) + &
                        3*F(16,2)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,1)*xi(1:nod,1)

                   !print*,'#P_3 approximation of Q boundary is implemented for', &
                   !     ' Jacobian but not tested !!!'
                   if(deg >=4) then
                      print*,'Compute_DF4 for deg > 3 is not implemented'
                      stop
                   endif
                endif
             endif
          endif

    else
          print*,'Only triagles and quarilaterals are implemented in ComputeDF'
          stop

    endif

  end subroutine ComputeDF2D

  !> draw a results of mapping \f$ F: \hat{K}\to K\f$ to file 'fort.(ifile)'
  !> visualizable by gnuplot, \f$ K= \f$ elem
  subroutine CheckElement2D( this, ifile)
    class( element ), intent(inout):: this   ! elem = element
    integer, intent(in) :: ifile
    real, dimension(:,:), allocatable :: xi, Fxi
    integer :: i,j, itest, it

    itest = 20

    if(this%type == 3) then
       allocate(xi(1:(itest+1)*(itest+2)/2,1:nbDim) )
       allocate(Fxi(1:(itest+1)*(itest+2)/2,1:nbDim) )
       it = 0
       do i=0,itest
          do j=0,itest-i
             it = it + 1
             xi(it,1) = 1.*i/itest
             xi(it,2) = 1.*j/itest
          enddo
       enddo

    elseif(this%type == 4) then
       allocate(xi(1:(itest+1)**2,1:nbDim) )
       allocate(Fxi(1:(itest+1)**2,1:nbDim) )
       it = 0
       do i=0,itest
          do j=0,itest
             it = it + 1
             xi(it,1) = 1.*i/itest
             xi(it,2) = 1.*j/itest
          enddo
       enddo
    else
       ! no work
       allocate(xi(1:1,1:nbDim) )
       allocate(Fxi(1:1,1:nbDim) )
    endif

    call this%computeF( it, xi(1:it,1:nbDim), Fxi(1:it,1:nbDim) )

    do i=1,it
       if(xi(i,1)*xi(1,2)*(1-xi(i,1) - xi(i,2)) < 1E-12) then
          write(ifile,*) Fxi(i,1:nbDim),xi(i,1:nbDim),i, xi(i,1)*xi(1,2)*(1-xi(i,1) - xi(i,2))
       endif
    enddo
    write(ifile,'(x)')

    deallocate(xi,Fxi)
  end subroutine CheckElement2D












!!!!!!!!!!!!!!!!!!! HG ELEMENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine initElemHG( this, i, type, flen, iloc)
      class( ElementHG_t) :: this
      integer, intent(in) :: i
      integer, intent(in) :: type
      integer, intent(in) :: flen
      integer, dimension(:), intent(in) :: iloc

      this%i = i
      this%type = type
      this%flen = flen

      allocate( this%face(1:max_face, 1:flen) )
      this%face(idx,1:flen) = iloc(1:flen)

         ! element with hanging nodes
         if ( flen > type ) then
            this%HGnode = .true.
            allocate( this%HGvertex(1:type ) )
            this%HGvertex(1:type) = iloc(1+flen: flen+type)
            ! call this%setHG( this%npoin, this%x, HG)
            ! call SetHGelem(elem, grid%npoin, grid%x)
         else
            this%HGnode = .false.
         endif


   end subroutine initElemHG


   !> setting of hanging nodes
   subroutine setHG( this, npoin, x, max_HGlevel, HG)
    class(ElementHG_t), intent(inout) :: this
    integer, intent(in) :: npoin
    real, dimension(1:npoin,1:nbDim), intent(in) :: x
    integer, intent(in) :: max_HGlevel  ! max implemented HG level (probably = 2)
    integer, intent(in) :: HG
    integer :: i, j, k, i1, i2, k1, k2, l1, l2
    integer :: iface, idiff, ilevel
    real :: r0, r1, r2, rpos1, rpos2



    if(.not. this%HGnode) then
       print*,'No hanging node'
       stop
    endif


    allocate(this%HGface(1:2, 1:this%flen) )

    iface = 0
    do i =1,this%type
       i1 = this%HGvertex(i)
       i2 = this%HGvertex( mod(i,this%type )+ 1)

       idiff = i2 - i1
       if(i == this%type) idiff = idiff + this%flen

       k1 = this%face(idx, i1)
       k2 = this%face(idx, i2)

       !write(*,'(a5,5i5,a1,2i5,a1,8i5)') &
       !     '&&&&&', this%i,i,i1,i2,idiff,'|',k1,k2,'|', this%HGvertex(:)

       if(idiff == 1) then ! face without hanging node
          iface = iface + 1
          this%HGface(1, iface) = i
          this%HGface(2, iface) = 1
       elseif(idiff >= 2 .and. idiff <= HG) then ! face with hanging nodes
          do j=1, idiff
             iface = iface + 1
             this%HGface(1, iface) = i

             l1 = this%face(idx, iface)
             l2 = this%face(idx, mod(iface,this%flen) + 1)

             if(dot_product(x(k1,1:nbDim)-x(k2,1:nbDim), x(l1,1:nbDim)-x(l2,1:nbDim)) <= 0.95* &
                  (dot_product(x(k1,1:nbDim)-x(k2,1:nbDim), x(k1,1:nbDim)-x(k2,1:nbDim) )* &
                  dot_product(x(l1,1:nbDim)-x(l2,1:nbDim), x(l1,1:nbDim)-x(l2,1:nbDim)))**0.5) then
                print*,'hanging node is not on a straight segment: u.v <= |u| |v|'
                write(*,'(a6,20i5)') 'elem =',this%i,this%flen,this%face(idx,1:this%flen)
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
             print*, 'PROBLEMMMM'
             do k=1,max_HGlevel
                !write(*,'(a20,i5,4es12.4)') '...',k, rpos2-rpos1,  - 2.**(-k), &
                !     ((rpos2-rpos1) - 2.**(-k))/(2.**(-k))
                if( abs(((rpos2-rpos1) - 2.**(-k))/(2.**(-k))) <= 1E-3) then
                   ilevel = k
                   goto 10
                endif
             enddo
             print*,'Face in sharing face does not found (1)',this%i,this%flen
             write(*,*) x(this%face(idx,this%HGvertex(1)), 1:nbDim)
             write(*,*) x(this%face(idx,this%HGvertex(2)), 1:nbDim)
             write(*,*) x(this%face(idx,this%HGvertex(3)), 1:nbDim)


             write(*,'(a6,20i5)') '****', i, i1,i2, idiff, k1, k2
             write(*,'(a6,20i5)') 'elem:', this%face(idx,:)
             write(*,'(a6,20i5)') 'HGel:', this%HGface(1,:), this%HGface(2,:)
             write(*,'(a6,20i5)') 'Hver:', this%HGvertex(:)

             stop
10           continue

             !write(*,'(a6,4i5,a1,2es12.4,i5)') 'diff',j,iface,l1,l2,'|', rpos1,rpos2,ilevel

             !!this%HGlevel = max(this%HGlevel, ilevel)
             ! seeking the segment within level
             do k=0,2**ilevel-1
                if(abs(rpos1 - k* 2.**(-ilevel))/(2.**(-ilevel)) <= 1E-3) then
                   this%HGface(2, iface) = 2.**ilevel + k
                   goto 20
                endif
             enddo
             print*,'Face in sharing face does not found (2)',this%i,this%flen
             write(*,*) x(this%face(idx,this%HGvertex(1)), 1:nbDim)
             write(*,*) x(this%face(idx,this%HGvertex(2)), 1:nbDim)
             write(*,*) x(this%face(idx,this%HGvertex(3)), 1:nbDim)


             write(*,'(a6,20i5)') '****', i, i1,i2, idiff, k1, k2
             write(*,'(a6,20i5)') 'elem:', this%face(idx,:)
             write(*,'(a6,20i5)') 'HGel:', this%HGface(1,:), this%HGface(2,:)
             write(*,'(a6,20i5)') 'Hver:', this%HGvertex(:)

             stop
20           continue


             !write(*,'(a3,4i3,6es10.2)') '###',i1, i2, k1,k2, x(k1,1:nbDim),x(k2,1:nbDim),rpos1,rpos2
             !write(*,'(a3,4i3,4es10.2, 2i3)') '###',j, iface, l1,l2, x(l1,1:nbDim),x(l2,1:nbDim),this%HGface(1:2,iface)
             !write(*,*)


          enddo
       else
          print*,'Troubles in elem%setHG'
          stop
       endif

    enddo
   end subroutine setHG

!!!!!!!!!!!!!!! RG ELEMENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine initElemRG( this, i, type, flen, iloc)
      class( ElementRG_t) :: this
      integer, intent(in) :: i
      integer, intent(in) :: type
      integer, intent(in) :: flen
      integer, dimension(:), intent(in) :: iloc

      this%i = i
      this%type = type
      this%flen = flen

      allocate( this%face(1:max_face, 1:flen) )
      this%face(idx,1:flen) = iloc(1:flen)

      ! for red-green refinement
      this%RGtype = 'N'
      this%RGhistory = 0
      this%RGlevel = 0
      this%RGindex = 0

      this%HGnode = .false.

   end subroutine initElemRG

!!!!!!!!!!!!!!!! AMA ELEMENT !!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine initElemAMA( this, i, type, flen, iloc)
      class( ElementAMA_t) :: this
      integer, intent(in) :: i
      integer, intent(in) :: type
      integer, intent(in) :: flen
      integer, dimension(:), intent(in) :: iloc

      this%i = i
      this%type = type
      this%flen = flen

      allocate( this%face(1:max_face, 1:flen) )
      this%face(idx,1:flen) = iloc(1:flen)

      this%HGnode = .false.


   end subroutine initElemAMA


    !> initialization  of elem%dof, elem%dof_plus elem%dof_fix for a given element
  subroutine initElementDof( this )
    class( element ), intent(inout):: this
    integer :: deg

    deg = this%deg

    this%dof  =  DOFtriang(deg)
    this%deg_plus = .false.

    this%dof_plus = DOFtriang( deg+1 )   !p+ := p + 1
    !this%dof_plus = DOFtriang( deg+2 )   !p+ := p + 2
    !this%dof_plus = DOFtriang( deg+3 )   !p+ := p + 3
    !this%dof_plus = DOFtriang( deg+4 )   !p+ := p + 4

    this%dof_fix = this%dof

    ! STDG
    this%Tdof      = this%Tdeg + 1
    this%Tdof_plus = this%Tdeg + 2

    ! saddle-point problems (SP)
    deg = this%degP
    if(deg >= 0)   this%dofP  =  DOFtriang(deg)


  end subroutine initElementDof



end module element_mod










!!
!!   subroutine writeElem(this)
!!      class(Element) :: this
!!
!!      print*, 'elem i = ', this%i
!!      print*, 'elem tp =', this%type
!!
!!   end subroutine writeElem
!
!!   subroutine setFaces( elem , lenloc, iloc )
!!      class(Element) :: elem
!!      integer, intent(in) :: lenloc
!!      integer, dimension(:), intent(in) :: iloc
!!
!!!ONLY for 2D now
!!!    if(nbDim == 2) then
!!          if(lenloc >= 3 .and. lenloc <= 4) then ! triangle or quadrilateral
!!             elem%type = lenloc
!!             elem%flen = lenloc
!!
!!             allocate( elem%face(1:elem%flen) )
!!
!!             !FR SET FACE
!!
!!             elem%face(1:lenloc)%idx = iloc(1:elem%flen)
!!             elem%HGnode = .false.
!!
!!          elseif(lenloc >= 9) then  ! element with hanging nodes
!!
!!               print*, 'HG nodes - Not implemented yet'
!!!             elem%type = iloc(1)
!!!             elem%flen = iloc(2)
!!!
!!!             allocate( elem%face( 1:elem%flen) )
!!!
!!!             elem%face( 1:iloc(2) )%idx = iloc(3:2+elem%flen)
!!!
!!!             elem%HGnode = .true.
!!!
!!!             allocate(elem%HGvertex(1:iloc(1)) )
!!!             elem%HGvertex(1:iloc(1)) = iloc(3+elem%flen: 2+elem%flen+elem%type)
!!
!!
!!           ! not yet implemented FRFR
!!           !  call elem%SetHG(grid%npoin, grid%x)
!!          else
!!            ! print*,'Bad data in init file ',gridfile,',  lenloc = ',lenloc
!!               print*,'Bad data in init file -  lenloc = ',lenloc
!!             stop
!!          endif
!!!       elseif(nbDim == 3) then   !3D case,  tetrahedra, pyramid or hexahedra
!!!          elem%type = lenloc
!!!          if (lenloc == 4) then ! tetrahedra
!!!             elem%flen = lenloc
!!!             allocate(elem%face(1:max_face, 1:elem%flen))
!!!             elem%face(nbnode,:) = 3
!!!
!!!             ! check orientation of tetrahedra
!!!             vec1(1:nbDim) = grid%x(iloc(2),1:nbDim) - grid%x(iloc(1),1:nbDim)
!!!             vec2(1:nbDim) = grid%x(iloc(3),1:nbDim) - grid%x(iloc(1),1:nbDim)
!!!             ! print*,'mmm ',grid%x(iloc(2),1:nbDim),grid%x(iloc(1),1:nbDim),vec1(1:nbDim)
!!!             vec3(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
!!!             vec3(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
!!!             vec3(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
!!!             vec1(1:nbDim) = grid%x(iloc(4),1:nbDim) - grid%x(iloc(1),1:nbDim)
!!!             if ((vec1(1)*vec3(1)+vec1(2)*vec3(2)+vec1(3)*vec3(3))<0.0) then
!!!                tmp=iloc(2)
!!!                iloc(2)=iloc(3)
!!!                iloc(3)=tmp
!!!                counter = counter +1
!!!             endif
!!!             !end check orientation of tetrahedra
!!!
!!!          elseif (lenloc == 5) then  ! pyramids
!!!             elem%flen = lenloc
!!!             allocate(elem%face(1:max_face, 1:elem%flen))
!!!             elem%face(nbnode,1:4) = 3
!!!             elem%face(nbnode,5) = 4
!!!          else                        ! hexahedra ??
!!!             elem%flen = lenloc-2
!!!             allocate(elem%face(1:max_face, 1:elem%flen))
!!!             elem%face(nbnode,:) = 4
!!!          end if
!!!          elem%face(idx,1:lenloc) = iloc(1:elem%flen)
!!!          elem%HGnode = .false.
!!!
!!!       endif
!!
!!       ! for red-green refinement - FRFR not yet implemented STRUCTURE?
!!!       elem%RGtype = 'N'
!!!       elem%RGhistory = 0
!!!       elem%RGlevel = 0
!!!       elem%RGindex = 0
!!
!!   end subroutine setFaces
!

