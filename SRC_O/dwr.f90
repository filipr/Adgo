module dwr_oper
  use matrix_oper
  use lapack_oper
  use main_data

  implicit none

!   type :: DWR_type
!     integer :: id !type of the DWR target functional (will be specified in paramets)
!     !  type(Newton_type), pointer :: Newton !Newton arrays for dual problem
!     real :: Ju ! target quantity of the computed solution
!     logical :: lin_functional ! target functional is linear -> no need of J'(u)(.)
!     integer, dimension(:), allocatable :: supp !support of the tarFunc
!     integer :: plus ! the difference if the primal and dual approximation degree
!     integer :: nsize ! size of the Dual problem
!     !integer :: tarFunc_id                  ! type of the DWR estimates
!     ! integer :: Tdeg, deg
!     !integer :: ttime, ctime ! not necessary for stationary problems
!     integer ::  nonzero                      ! # of nonzero  elements
!     ! integer :: max_Qdof                      ! maximal of elem(i)%Qdof
!     ! integer :: max_dof                       ! maximal of elem(i)%dof
!     ! integer :: max_Tdof                      ! maximal of elem(i)%Tdof
!     ! real :: h                               ! maximal diameter
!     real :: eps                              ! diameter for pointvalue of tarFunc for id=3
!     real, dimension(1:2) :: xy_coord ! coordinates of the point value in tarFunc for id=3
!
!   end type DWR_type

!   public :: PrepareDualProblem


   public :: ComputeTarFuncPhi
   public :: ComputeJu
!   public :: FindSupp

  ! type :: TargetFunc
!   end type TargetFunc

   contains

!   subroutine PrepareDualProblem( )
!
!   ! if (state%DWR%lin_functional ) then
!
!         select case(state%DWR%id)
!            case(1)
!               print*, 'Not implemented'
!               stop
!            case(2)
!               print*, 'Not implemented'
!               stop
!
!            case(3)
!               print*, 'Point value in ', state%DWR%xy_coord(1), ' , ', state%DWR%xy_coord(2)
!               state%DWR%eps = 0.05
!               print*, 'epsilon set to ', state%DWR%eps
!               call FindSupp(state%DWR%eps)
!
!
!            case(4:)
!               print*, 'Not implemented in PrepareDualProblem'
!         end select
!
!
!!      else
!!
!!      endif
!
!
!
!
!
!   end subroutine PrepareDualProblem

!   subroutine InitDualMesh
!
!   !periodicity
!
!   grid%x, grid%xcur
!
!   end subroutine InitDualMesh



   subroutine ComputeTarFuncPhi()




   end subroutine ComputeTarFuncPhi

   subroutine ComputeJu()



   end subroutine ComputeJu


!   subroutine FindSupp(grid,eps)
!      class(mesh), intent(in) :: grid
!      real, intent(in) :: eps
!      class(element), pointer :: elem
!      integer :: i, j
!      integer, dimension(:), allocatable :: support
!
!
!
!      allocate (support(1:grid%nelem))
!      ! support(:) = -1
!      j = 1
!      do i = 1, grid%nelem
!
!         if ( Distance( grid%elem(i)%xc, state%DWR%xy_coord ) < eps ) then
!            support(j) = i
!            j = j+1
!         endif
!
!      enddo
!
!      if ( allocated(state%DWR%supp) ) then
!         deallocate( state%DWR%supp)
!      endif
!
!      allocate( state%DWR%supp(1: j-1 ) )
!      state%DWR%supp(1:j-1) = support(1:j-1)
!
!      deallocate( support )
!
!   end subroutine FindSupp

end module dwr_oper
