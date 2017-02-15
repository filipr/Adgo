module solution_mod
   use mesh_mod
   use nonlinear_mod
   use define_state
   use stdgm_mod


implicit none

   public :: allocateZST
   public :: copyZST_fromLongVector
   public :: CopyZST_toLongVector
   public :: Transfer_w_to_WS_elem
   public :: Solution_L8norm

contains

   subroutine allocateZST( grid , plus)
      class( mesh ), intent(inout) :: grid
      logical, intent(in) :: plus
      class(element), pointer :: elem
      integer :: i, deg, dof, Tdof

      do i = 1, grid%nelem
         elem => grid%elem(i)
         if (plus) then
            deg = elem%deg + state%space%plusDeg
            dof = DOFtriang( deg )
         else
            dof = elem%dof
         endif
         Tdof = elem%Tdof

         if ( associated(elem%zST) ) then
            print*, 'allocateZST:', size(elem%zST)
            stop
            deallocate(elem%zST)
         endif

         allocate( elem%zST(1:ndim,1:dof,1:Tdof), source = 0.0 )
       end do ! i

   end subroutine allocateZST

   !> distributes  elem%zST to the global vector x to
   subroutine copyZST_toLongVector( grid,  nsize, x )
      class( mesh ) , intent(in) :: grid
      integer, intent(in) :: nsize
      real, dimension(1:nsize), intent(out) :: x
      class( element ), pointer :: elem
      integer :: i,j,k,kk, dof, elemDof

      kk = 0

      do i = 1, grid%nelem
        elem => grid%elem(i)
        dof = elem%dof

        if (elem%ncv /= kk + 1) then
          print*, 'Problem while copying update from zST (solution.f90)'
          stop
        endif
        elemDof = elem%dof * ndim * elem%Tdof

        do k = 1, elem%Tdof
          do j = 1, ndim
             x(kk + 1 : kk + dof) = elem%zST(j,1:dof,k)

             kk = kk + dof
          enddo !j
        enddo !k

      enddo !i
      nullify( elem )

   end subroutine copyZST_toLongVector

   !> distributes the global solution vector x to elem%zST and elem%w(0,:)
   subroutine CopyZST_fromLongVector( grid,  nsize, x )
      class( mesh ) , intent(inout) :: grid
      integer, intent(in) :: nsize
      real, dimension(1:nsize), intent(in) :: x
!      class( DWR_t ) , intent(in) :: DWR
      class( element ), pointer :: elem
      integer :: i,j,k,kk, dof, elemDof

      kk = 0

      do i = 1, grid%nelem
        elem => grid%elem(i)
        dof = elem%dof

        if (elem%ncv /= kk + 1) then
          print*, 'Problem while copying update to wST (solution.f90)'
          stop
        endif
        elemDof = elem%dof * ndim * elem%Tdof

        do k = 1, elem%Tdof
          do j = 1, ndim
             elem%zST(j,1:dof,k) = x(kk + 1 : kk + dof)
             kk = kk + dof
          enddo !j
        enddo !k
        !  compute the solution in the endpoints

!        call Transfer_wST_to_w_Elem( elem , 0, elem%TQnum)
        elem%w(0, 1:dof*ndim) = Transfer_funST_to_fun( &
            elem%zST(1:ndim, 1:dof, 1:elem%Tdof), &
            dof, elem%Tdof, 0, elem%TQnum)

      enddo !i

      nullify( elem )

   end subroutine CopyZST_fromLongVector

   !copy elem%w(0,:) to wS
   subroutine Transfer_w_to_WS_elem( elem )
     class( element ), intent(inout) :: elem
     integer :: j, dof

     dof = elem%dof

     if ( allocated( elem%wS ) ) then
         deallocate( elem%wS )
     endif

     allocate( elem%wS(1:ndim, 1:elem%dof) )

     do j = 1, ndim
         elem%wS(j, 1:dof) = elem%w( 0, (j-1)*dof + 1 : j*dof )
     enddo !j

   end subroutine Transfer_w_to_WS_elem


   !> compute (pseudo maximal norm) - in space integration nodes
   !> in the time interval endpoint
  function Solution_L8norm(grid) result(norm)
   class( mesh ), intent(in) :: grid
   real :: norm

   class( element ), pointer :: elem
   real, allocatable, dimension(:,:)  :: wi
   real :: temp
   integer :: i

      norm = 0.0

      do i = 1,grid%nelem
         elem => grid%elem(i)
         allocate( wi(1:elem%Qdof, 1:ndim), source = 0.0 )

         call Eval_whST_Elem( elem, 1, 0, wi)
         temp = maxval( abs( wi(:,:) ) )
         norm = max( norm, temp)
         deallocate( wi )
      end do
  end function Solution_L8norm







end module solution_mod
