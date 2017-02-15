! NOT USED
!> basis geometrical entities which can be stored in parallel mode

module elemental_mod
  use paramets

  implicit none

   private
  !> elemental arrays - 1 dimensional
   type, public :: Elemental1_t
      real, allocatable, dimension (:) :: x  ! block of matrix

      contains
!      procedure :: addElemental1
!      procedure :: subElemental1
!      procedure :: assignElemental1
      procedure :: length
      procedure :: printme
      procedure :: delete => deleteElemental1
      ! final :: deleteElemental1

!      generic, public :: operator(+) => addElemental1
!      generic, public :: operator(-) => subElemental1
!      generic, public :: assignment(=) => assignElemental1
!     procedure :: CopyMblocK
    end type Elemental1_t

    interface Elemental1_t
       procedure :: constructFromNumber
       procedure :: constructFromArray
    end interface

    interface assignment (=)
      procedure :: copyElemental1
      procedure :: copyArrayToElemental1
!!      procedure :: copyElemental1ToArray
    end interface

    public :: assignment(=)



    ! 3 dimensional arrays
    !> elemental arrays - 3 dimensional
   type, public :: Elemental3_t
      real, allocatable, dimension (:,:,:) :: x  ! block of matrix

      contains

!      procedure :: addElemental1
!      procedure :: subElemental1
!      procedure :: assignElemental1
      procedure :: getSize => getSize3
      procedure :: printme => printme3
      procedure :: copyTo1Darray => copyElemental3to1Darray ! watch for the ordering
      procedure :: delete => deleteElemental3

!      generic, public :: operator(+) => addElemental1
!      generic, public :: operator(-) => subElemental1
!      generic, public :: assignment(=) => assignElemental1
!     procedure :: CopyMblocK
    end type Elemental3_t

    interface Elemental3_t
       procedure :: constructFromNumber3
       procedure :: constructFromArray3
    end interface

    interface assignment (=)
      procedure :: copyElemental3
      procedure :: copyArrayToElemental3
!!      procedure :: copyElemental1ToArray
    end interface

!    public :: assignment(=)

!   !> elemental arrays - 2 dimensional
!   type, public :: Elemental2_t
!      real, allocatable, dimension (:,:)  :: x  ! block of matrix
!
!      contains
!      procedure :: init => initElemental2
!      final :: del => deleteElemental2
!   end type Elemental2_t
!
!   !> elemental arrays - 3 dimensional
!   type, public :: Elemental3_t
!      real, allocatable, dimension (:,:,:)  :: x  ! block of matrix
!
!      contains
!      procedure :: init => initElemental3
!      final :: del => deleteElemental2
!  end type Elemental3_t

   public :: copy3Darrayto1Darray
   public :: copy1DarrayTo3Darray


contains

   subroutine printme(this)
      class( Elemental1_t ), intent(inout) :: this

      write(*,*) 'Elemental_t%x:', this%x

   end subroutine printme

   elemental function length( this)
      class( Elemental1_t ), intent(in) :: this
      real :: length
      length = size( this%x(:) )
   end function length

!   subroutine initElemental3( this, deg1, deg2,deg3)
!      class( Elemental3_t ), intent(inout) :: this
!      integer, intent(in) :: deg1
!      integer, intent(in) :: deg2
!      integer, intent(in) :: deg3
!
!      allocate(this%x(1:deg1,1:deg2,1:deg3), source = 0.0)
!
!   end subroutine initElemental3

    !> n = length of the vector
    function constructFromNumber(n) result(this)
        type( Elemental1_t ) :: this
        integer, intent(in) :: n
        if (allocated(this%x)) &
            deallocate(this%x)

        allocate(this%x(1:n), source=0.0)
    end function

    function constructFromArray(array) result(this)
        type( Elemental1_t ) :: this
        real, dimension(:), intent(in) :: array
        allocate(this%x(1:size(array)))
        this%x = array

    end function

!    function constructCopy( other ) result(this)
!        type( Elemental1_t ), intent(in) :: other
!        type( Elemental1_t ) :: this
!
!        allocate( this%x(1: other%length()) )
!        this%x = other%x
!    end function
!

   subroutine copyElemental1(this,other)
      class(Elemental1_t),intent(inout) :: this
      class(Elemental1_t),intent(in) :: other

      if (allocated(this%x)) then
         deallocate(this%x)
      end if

      allocate( this%x(1: other%length() ) )
      this%x = other%x
    end subroutine

    subroutine copyArrayToElemental1(this,array)
      class(Elemental1_t),intent(inout) :: this
      real, dimension(:),intent(in) :: array

      call copyElemental1(this, Elemental1_t( array ))
    end subroutine


   subroutine deleteElemental1(this)
      class( Elemental1_t ) :: this
      if (allocated(this%x)) &
            deallocate(this%x)

    end subroutine deleteElemental1





!!!!!!!!!!!! 3 !!!!!!!!!!!!!!!!!!!
   subroutine printme3(this)
      class( Elemental3_t ), intent(inout) :: this

      write(*,*) 'Elemental3_t%x:', this%x

   end subroutine printme3

   function getSize3( this) result( length )
      class( Elemental3_t ), intent(in) :: this
      real, dimension(1:3) :: length
      length(1) = size( this%x(:,1,1) )
      length(2) = size( this%x(1,:,1) )
      length(3) = size( this%x(1,1,:) )
   end function getSize3


   subroutine deleteElemental3(this)
        class( Elemental3_t ) :: this
        if (allocated(this%x)) then
            deallocate(this%x)
        end if
    end subroutine deleteElemental3

          !> n = length of the vector
    function constructFromNumber3(n,m,k) result(this)
        type( Elemental3_t ) :: this
        integer, intent(in) :: n
        integer, intent(in) :: m
        integer, intent(in) :: k

        if (allocated(this%x)) then
         deallocate(this%x)
        end if

        allocate(this%x(1:n,1:m,1:k), source=0.0)
    end function constructFromNumber3

    function constructFromArray3(array) result(this)
        type( Elemental3_t ) :: this
        real, dimension(:,:,:), intent(in) :: array
        integer :: k,l,m

        k = size(array(:,1,1))
        l = size(array(1,:,1))
        m = size(array(1,1,:))
        allocate( this%x(1:k,1:l,1:m), &
                   source = array(1:k,1:l,1:m) )

    end function constructFromArray3


   subroutine copyElemental3(this,other)
      class(Elemental3_t),intent(inout) :: this
      class(Elemental3_t),intent(in) :: other
      integer :: k,l,m

      if (allocated(this%x)) then
         deallocate(this%x)
      end if
      k = size(other%x(:,1,1))
      l = size(other%x(1,:,1))
      m = size(other%x(1,1,:))

      allocate( this%x(1:k,1:l,1:m), &
                source = other%x(1:k,1:l,1:m) )

    end subroutine copyElemental3

    subroutine copyArrayToElemental3(this,array)
      class(Elemental3_t),intent(inout) :: this
      real, dimension(:,:,:),intent(in) :: array
      integer :: k,l,m

      if (allocated(this%x)) then
         deallocate(this%x)
      end if
      k = size(array(:,1,1))
      l = size(array(1,:,1))
      m = size(array(1,1,:))

      allocate( this%x(1:k,1:l,1:m), &
                source = array(1:k,1:l,1:m) )

    end subroutine copyArrayToElemental3


    ! copy Elemental3_t to 1D array , permutation sets the ordering
    ! standard permutation for STDGM (/ 3,1,2 /)
    ! watch for the ordering 1:Tdof, 1:ndim, 1:dof - 3,1,2
    function copyElemental3to1Darray( this, nsize) result(array)
      class(Elemental3_t),intent(in) :: this
      integer, intent(in) :: nsize
!      integer, dimension(1:3), intent(in) :: permutation ! order of the cycles
      real, dimension(1:nsize) :: array
      integer, dimension(1:3) :: length
      integer :: i, k
      integer :: m, n , p, nn

      length = getSize3( this )

      m = length(1)
      n = length(2)
      p = length(3)

      if ( size( array(:) ) /= m*n*p ) &
         stop 'wrong size of the array in copyElemental3to1Darray'

      nn = 0

      do i = 1, p ! start with the 3rd index (TDOF)
         do k = 1, m ! 1st index - ndim
            array(nn+1:nn+n) = this%x(k,1:n,i )
         end do !k
         nn = nn + n
      end do ! n

    end function copyElemental3to1Darray


    ! copy3D array to 1D array
    ! standard permutation for STDGM (/ 3,1,2 /)
    ! watch for the ordering 1:Tdof, 1:ndim, 1:dof - 3,1,2
    function copy3DarrayTo1Darray( array3, nsize) result(array)
      real, dimension(:,:,:), intent(in) :: array3
      integer, intent(in) :: nsize
      real, dimension(1:nsize) :: array
      integer, dimension(1:3) :: length
      integer :: i, k
      integer :: m, n , p, nn

      m = size( array3(:,1,1) )
      n = size( array3(1,:,1) )
      p = size( array3(1,1,:) )

      if ( nsize /= m*n*p ) then
         print*, nsize, m,n,p
         stop 'wrong size of the array in copyElemental3to1Darray'
      endif
      nn = 0

      do i = 1, p ! start with the 3rd index (TDOF)
         do k = 1, m ! 1st index ndim ,
            array(nn+1:nn+n) = array3(k,1:n,i ) !2nd index - dof
         end do !k
         nn = nn + n
      end do ! n

    end function copy3DarrayTo1Darray

    ! distribute 1D array to 3D array
    ! the ordering it typical for the STDG setting (/ 3,1,2 /)
    ! 1:ndim, 1:dof, 1:Tdof is ordered Tdof->ndim->dof
    function copy1DarrayTo3Darray( array1, n1, n2, n3 ) result (array3)
      real, dimension(:), intent(in) :: array1
      integer, intent(in) :: n1, n2, n3
      real, dimension(1:n1,1:n2,1:n3) :: array3

      integer :: i, k
      integer :: nn
      nn = 0

      do i = 1, n3 ! start with the 3rd index (TDOF)
         do k = 1, n1 ! 1st index ndim ,
            array3(k, 1:n2,i) = array1(nn+1:nn+n2) !2nd index - dof
         end do !k
         nn = nn + n2
      end do ! n

      if ( size(array1(:)) /= nn) then
         print*,  n1,n2,n3, nn, size(array1(:))
         stop 'wrong sizes in copy1DarrayTo3Darray'
      endif

    end function copy1DarrayTo3Darray

end module elemental_mod
