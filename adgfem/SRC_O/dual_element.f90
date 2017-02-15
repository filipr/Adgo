!module dual_element_mod
!
!implicit none
!
! !> basic element type, when adaptation needed, use some of its descendants
!   type :: Dual_element_t
!      integer :: i
!!      integer :: deg                   ! degree of polynomial approximation
!!      integer :: Tdeg    ! degree of polynomial approximation in time
!!      real, pointer, dimension(:,:) :: w               ! solution vector
!
!      real, pointer, dimension(:,:,:) :: primal_wST             ! solution vector for ST DGM - 1:ndim, 1:elem%dof, 1:elem%Tdof
!
!      real, pointer, dimension(:,:,:) :: zST        ! try
!
!!      real, pointer, dimension(:,:) :: wSTfin       ! solution vector for ST DGM at final time in space integ nodes
!!      real, pointer, dimension(:,:) :: wSTfinAD       ! solution vector for ST DGM at final time - coefficients
!!      real, pointer, dimension(:,:,:) :: rhsST             ! RHS  ST DGM - 1:ndim, 1:elem%dof, 1:elem%Tdof
!!
!!      real, pointer, dimension(:) :: vec             ! vectors: right hand sides, only one dimension
!
!     contains
!
!!     procedure :: allocVectors
!     procedure :: init => initDualElem
!
!
!   end type Dual_element_t
!
!contains
!
!   subroutine initDualElem( this )
!      class( Dual_element_t ) :: this
!
!      print*, 'init Dual elem'
!      print*, 'wST - pointers, therefore we dont allocate it'
!
!   end subroutine initDualElem
!
!end module
