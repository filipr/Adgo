module element3D_mod
   use element_mod

   implicit none


   type, public, extends (AbstrElement_t) :: Element3D_t
     integer :: i
     integer :: type                           ! =3 for triang, 4 for quad
     integer :: flen                           ! number of nodes including hanging
     !type(FaceElem_t), allocatable, dimension(:) :: face ! indexes of vertices, neighs,....
     integer, allocatable, dimension(:,:) :: face ! indexes of vertices, neighs,....
     real, allocatable, dimension(:,:, :) :: xi    ! integ nodes
     integer :: isupp                        ! number  of elements having a common vertex with elem
     integer, allocatable, dimension(:,:) :: supp  ! list of elements having a common vertex with elem
     real :: area                              ! area of an element
     real :: diam                              ! diameter of an element
     real :: r_ins                             ! radius of inscribed circle
     real :: limit_par                         ! parameter for limiting

   !  type ( Elemental_estims_t ), allocatable :: estims

!!     real :: estimA, estimS, estimT, estimST, estim_loc ! estimates of "residal errors"
!!     real :: jumpsJh                           ! \sum_{G \in \partial K} |G|^{-1} \int_G [u_h]^2 dS
!!     real :: rezid, rezid1                     ! reziduum for artificial viscosity
!     real :: reg, regT0, regT1, regT2          ! regularity of the solution
!     !real, pointer, dimension(:)  :: TEST
!     real :: errL2, errH1, errL8,interL8, interLq
! !    real, dimension(1:2, 1:3) :: errSnorm	  ! STDGM L2/H1 norms in three different nodes 1-endpoint, 2-node of the time quad, 3-1/2
     integer :: ibcur, jcur, deg_cur, dof_cur  ! = 0 if polygonal, > 0 curved ib
!
     type(Fmapping_t), pointer :: F              ! mapping: K^ -> K
     real, allocatable, dimension(:,:) :: n        ! vector of outer normal
     real, allocatable, dimension(:) :: dn         ! size of vector of outer normal
!     real, allocatable, dimension(:,:) :: nc ! outer normals on curved edge in integ nodes
!     real, allocatable, dimension(:) :: dnc  ! dS "increase" on curved edge in integ nodes
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
!     type(Mblock), pointer, dimension (:,:) :: SPblock ! matrix blocks (0:rlen)
     type(Mblock), pointer, dimension (:) :: XXX   ! preconditioner matrix blocks (0:rlen)
     real, allocatable, dimension(:, :) :: wc               ! computational approximation (solution) vector satisfying stop. criteria based on AEE
     real, allocatable, dimension(:,:) :: w               ! solution vector
!     real, allocatable, dimension(:,:,:) :: wSP               ! solution vector for SP - 1:3 (iV1,iV2,iP), 1:2 (time levels we need to remember), 1:max(elem%dof,elem%dofP)
     real, allocatable, dimension(:,:,:) :: wST             ! solution vector for ST DGM - 1:ndim, 1:elem%dof, 1:elem%Tdof
!     real, allocatable, dimension(:,:) :: wSTfin       ! solution vector for ST DGM at final time in space integ nodes
     real, allocatable, dimension(:,:) :: wSTfinAD       ! solution vector for ST DGM at final time - coefficients
     real, allocatable, dimension(:,:,:) :: rhsST             ! RHS  ST DGM - 1:ndim, 1:elem%dof, 1:elem%Tdof
!  !?   real, allocatable, dimension(:,:) :: wS              ! temporary used array
!   !?  real, pointer, dimension(:,:,:) :: wSS            ! temporary used array
     real, allocatable, dimension(:,:) :: vec             ! vectors: right hand sides, ...
     real, pointer, dimension(:,:) :: res_func    !residual function made up from residual vector
     real, pointer, dimension(:,:,:) :: MGw           ! solution vector for MG

!     real, pointer, dimension(:,:,:) :: RTNflux     ! RTN fluxes with time history
!     integer :: RTNflux_deg
!     real, pointer, dimension(:,:) :: RTNphi     ! RTN fluxes with time history
!     real :: errTot, errDisc, errAlg                ! total error u - u_h^i, discretization error u - u_h, algebraic error u_h - u_h^i on the element
!     real :: estTot, estDisc, estAlg, estRem        ! estimate of total error u - u_h^i, discretization error u - u_h, algebraic error u_h - u_h^i on the element
!     real, allocatable, dimension(:) :: estimFNCD   ! Discretization component of Flux NonConformity estimator
!     real, allocatable, dimension(:) :: estimFNCA   ! Algebraic component of Flux NonConformity estimator
!     real, allocatable, dimension(:) :: estimFNCT   ! Total Flux NonConformity estimator
!     integer :: MGdeg                               ! MG level for element
!     integer :: MGdof                               ! DOFs for current MG level
!     integer :: MGncv                               ! sequence of element's dof
     integer, allocatable, dimension(:) :: iBC          ! index of boundary component, if any
     integer, allocatable, dimension(:) :: tBC          ! type of boundary condition, if any
     integer       ::  ncv                          ! sequence of element's dof
     integer       ::  ncvP                          ! sequence of element's dofP

     integer :: hsplit                             ! type of h-ref of h-deref
     integer :: psplit                             ! type of p-ref
     logical :: to_recompute                  ! recompute the solution

!     real    :: ama_p                               ! change of deg
!     integer :: RGhistory                           ! levels of history  of RG refinement
!     integer :: RGlevel                             ! level of refinement
!     type(RGref), pointer, dimension(:) :: RG     ! data for red green (de)refinement
!     integer :: RGreclen                        ! number of elements in RGrecomputes
!     integer, dimension(1:4) :: RGrecomput      ! (old) index of elems for recomputation
!     character :: RGtype                       ! =N (none) =R (red), =G (green) refinement
!     integer :: RGindex                        ! index within level ref (see type(RGref)
!     logical :: dealloc                        ! elem (of old mesh) will be deallocated
!     real, allocatable, dimension(:,:) :: eta    ! array of different error indicators
     real :: CK, Cb, CTn, CKo, Cbo
     real :: CP                                !Poincare constant
!     logical :: RTNfluxeval                    !whether or not RTNflux for elem at state%time%iter is already evaluated
!     real, dimension(:),  allocatable :: rgabc ! metric for AMA adaptation
!     real :: amaL                             !  ama ... size of the dominant eigenvalue
!     real :: amaR                             !ama < 1 ... ratio between the small and the big eigenvals
!     real, dimension(1:2) :: amaD             ! = ( cos(alpha), sin(alpha) ) orientation of AMA ellipse
     integer :: per                           ! > 0 periodic neighbours element

     real :: max_eigenvals

   contains
      procedure :: init => initElem3D
      procedure :: computeF => ComputeF3D
      procedure :: computeDF => ComputeDF3D
      procedure :: setF => SetF3D



   end type Element3D_t

   contains

   subroutine initElem3D( this, i, type, flen, iloc )
      class( Element3D_t) :: this
      integer, intent(in) :: i
      integer, intent(in) :: type
      integer, intent(in) :: flen
      integer, dimension(:), intent(in) :: iloc

      stop 'abstract initElem3D this should not be called'
   end subroutine initElem3D


  !> setting of mappings \f$ F: \hat{K}\to K\quad \forall K\in{\cal T}_h\f$, i.e.
  !> mapping of the reference element onto the physical one
  !>
  !> \f$ K= \f$ elem,  given by Langrangian nodes within  TETRAHEDRA,
  !> \f$ F \f$ can be linear, quadratic, cubic, ....,
  !> F(l,1:nbDim) means  the \f$(x,y)\f$ coordinates of the l\f$^{\rm th}\f$
  !> Langrangian nodes
  subroutine SetF3D( this , nod, x)
    class( Element3D_t ), intent(inout) :: this   ! elem = element
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

    if(this%type == 4 .and. nbDim == 3) then   !  TETRAHEDRA

       if(deg == 1) then
          F(1, 1:nbDim) = x(1,1:nbDim)
          F(2, 1:nbDim) = x(2,1:nbDim) - x(1,1:nbDim)
          F(3, 1:nbDim) = x(3,1:nbDim) - x(1,1:nbDim)
          F(4, 1:nbDim) = x(4,1:nbDim) - x(1,1:nbDim)

          this%F%iFlin = .true.   ! linear mapping
       else
          print*,'Maximal mapping P_1 on tetrahedra implemented'
          stop
       endif
     else
         stop 'problem in SetF3D'
     endif

  end subroutine SetF3D


   !> evaluation of \f$ F(x_i)\in R^2,\ x_i\in \hat{K},\ i=1..nod,  \f$
  !>
  !> xi= \f$x_i\in \hat{K}\f$ are arbitrary nodes within reference elements,
  !> \f$ F \f$ mapping of reference element \f$\hat{K}\f$ to the actual one \f$K\f$=elem,
  !> F can be linear, quadratic, cubic, ....
  subroutine ComputeF3D( this, nod, xi, Fx )
    class( Element3D_t), intent(in):: this   ! elem = element
    integer, intent(in) :: nod            ! number of nodes
    real, dimension(1:nod, 1:nbDim), intent(in) :: xi
    real, dimension(1:nod, 1:nbDim), intent(out) :: Fx
    integer :: i, deg

    real, dimension(:, :), pointer :: F

    F => this%F%F(1:this%F%dof, 1:nbDim)
    deg = this%F%deg

    if(nbDim == 3) then ! 3D elements

       if(this%type == 4) then  ! TRETRAHEDRA
          do i=1, nbDim
             ! absolute term
             Fx(1:nod, i) = F(1, i)

             if(deg >= 1) then
                ! linear terms
                Fx(1:nod, i) = Fx(1:nod, i)  &
                     + F(2, i)*xi(1:nod, 1) + F(3, i)*xi(1:nod, 2) + F(4, i)*xi(1:nod, 3)

                !    if(deg >= 2) then
                !     ! quadratic terms
                !     Fx(1:nod, i) = Fx(1:nod, i)  &
                !             + F(4, i) * xi(1:nod, 1) * xi(1:nod, 1) &
                !             + F(5, i) * xi(1:nod, 1) * xi(1:nod, 2) &
                !             + F(6, i) * xi(1:nod, 2) * xi(1:nod, 2)
                if(deg >=2) then
                   print*,'Only polynoms up to order 1 are implemented',&
                        'in computeF for tetrahedra'
                   stop
                endif
                !                    endif
             endif
          enddo
       else
          print*,'Only tetrahedra are implemented in ComputeF'
          stop
       endif
    else
      stop 'problem in ComputeF3D'
    endif

  end subroutine ComputeF3D

  !> evaluation of \f$ \frac{D F(x_i)}{D \hat{x}} \in R^{2\times 2},
  !> \ x_i\in \hat{K},\ i=1..nod,  \f$
  !>
  !> xi= \f$x_i\in \hat{K}\f$ are arbitrary nodes within reference elements,
  !> \f$ F \f$ mapping of reference element \f$\hat{K}\f$ to the actual one \f$K\f$=elem,
  !> F can be linear, quadratic, cubic, ....
  subroutine ComputeDF3D( this, nod, xi, DF  )
    class( Element3D_t ), intent(in):: this   ! elem = element
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



    if (nbDim==3) then  ! 3D case

       if(this%type == 4) then  ! TETRAHEDRA
          if(deg >= 1) then
             ! linnear mapping
             do i=1,3
                DF(1:nod,i,1) = F(2,i)
                DF(1:nod,i,2) = F(3,i)
                DF(1:nod,i,3) = F(4,i)
             enddo

             !if(deg >= 2) then
             !   ! quadratic mapping
             !   DF(1:nod,1,1) = DF(1:nod,1,1) + 2*F(4,1)*xi(1:nod,1)+ F(5,1)*xi(1:nod,2)
             !   DF(1:nod,1,2) = DF(1:nod,1,2) + 2*F(6,1)*xi(1:nod,2)+ F(5,1)*xi(1:nod,1)
             !   DF(1:nod,2,1) = DF(1:nod,2,1) + 2*F(4,2)*xi(1:nod,1)+ F(5,2)*xi(1:nod,2)
             !   DF(1:nod,2,2) = DF(1:nod,2,2) + 2*F(6,2)*xi(1:nod,2)+ F(5,2)*xi(1:nod,1)
             if(deg >=2) then
                print*,'ComputeDF_3 for deg > 1 is not implemented'
                stop
             endif
             !endif
          endif
       else
          print*,'Only tetrahedra for 3D are implemented in ComputeDF'
          stop
       endif
    else
      stop 'Problem in ComputeDF3D'
    endif
   end subroutine ComputeDF3D

end module element3D_mod
