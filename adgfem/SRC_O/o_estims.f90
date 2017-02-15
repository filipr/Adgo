!> module for the error estimations
module estims_mod
   use paramets

   implicit none

!   type, public :: SpaceEstims_t
!      real, allocatable, dimension(:,:) :: estim
!
!   contains
!      procedure :: init => initSpaceEstims
!
!   end type SpaceEstims_t
!
!   type, public, extends( SpaceEstims_t ) :: DWR_estim_t
!!      integer :: id !type of the DWR target functional (will be specified in paramets)
!!      character(len=20) :: name
!!      real :: Ju ! target quantity of the computed solution
!!      logical :: lin_functional ! target functional is linear -> no need of J'(u)(.)
!
!      !     integer, dimension(:), allocatable :: supp !support of the tarFunc
!!      integer :: plus ! the difference if the primal and dual approximation degree
!
!!      real :: eps                              ! diameter for pointvalue of tarFunc for id=3
!!      real, dimension(1:2) :: xy_coord ! coordinates of the point value in tarFunc for id=3
!
!   contains
!      procedure :: init => initDWR_estims
!
!   end type DWR_estim_t


contains

!   subroutine initSpaceEstims( this, Lq)
!      class(SpaceEstims_t), intent(inout) :: this
!      real, optional, intent(in) :: Lq
!
!      stop 'General initSpaceEstims is not implemented, should be abstract in future.'
!
!   end subroutine initSpaceEstims
!
!   subroutine initDWR_estims( this, Lq)
!      class( DWR_estim_t ), intent(inout) :: this
!      real, optional, intent(in) :: Lq
!
!      print*, 'initDWR_estims not used- DWR method is in separate structure'
!   end subroutine initDWR_estims





end module estims_mod
