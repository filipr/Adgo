!> module for the error estimations
module estims_mod
   use paramets

   implicit none

   type, public :: SpaceEstims_t
      real, allocatable, dimension(:,:) :: estim

   contains
      procedure :: init => initSpaceEstims

   end type SpaceEstims_t

   type, public, extends( SpaceEstims_t ) :: DWR_t
      integer :: id !type of the DWR target functional (will be specified in paramets)
      character(len=20) :: name
      real :: Ju ! target quantity of the computed solution
      logical :: lin_functional ! target functional is linear -> no need of J'(u)(.)

      !     integer, dimension(:), allocatable :: supp !support of the tarFunc
      integer :: plus ! the difference if the primal and dual approximation degree

      real :: eps                              ! diameter for pointvalue of tarFunc for id=3
      real, dimension(1:2) :: xy_coord ! coordinates of the point value in tarFunc for id=3

   contains
      procedure :: init => initDWR

   end type DWR_t


contains

   subroutine initSpaceEstims( this, Lq)
      class(SpaceEstims_t), intent(inout) :: this
      real, optional, intent(in) :: Lq

      stop 'General initSpaceEstims is not implemented, should be abstract in future.'

   end subroutine initSpaceEstims

   subroutine initDWR( this, Lq )
      class(DWR_t), intent(inout) :: this
      real, optional, intent(in) :: Lq

      if (.not. present(Lq) ) &
         stop 'Lq must be present when calling initDWR'


      !dual_deg - primal deg
      this%Ju = 0
      this%plus = 0
      !setting the tarFunc

      this%id = floor(Lq)
      select case(this%id)
         case(1)
            write(*,'(a50)') '  # Target functional: L2-norm of the solution - NOT IMPLEMENTED'
            this%lin_functional = .false.
            stop
         case(2)
            write(*,'(a50)') '  # Target functional: H1-seminorm of the solution - NOT IMPLEMENTED'
            this%lin_functional = .false.
            stop
         case(3)
            write(*,'(a50)') '  # Target functional: Point value of the solution in (0.5, 0.5)'
            this%xy_coord(1:2) = (0.5, 0.5)
            this%lin_functional = .true.
            this%name = 'PointValue'
            this%eps = 0.1
         case(4:)
            write(*,'(a50)') '  # Unknown type of target functional! STOPPING!'
            stop
      end select

   end subroutine initDWR




end module estims_mod
