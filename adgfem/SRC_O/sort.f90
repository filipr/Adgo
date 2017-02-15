module sort_mod

implicit none




public :: MergeArrays
!public :: MergeSort
private :: heapsort_int
private :: heapsort_real
private :: shiftdown_int
private :: shiftdown_real

interface heapsort
   module procedure heapsort_int, heapsort_real
end interface

contains

   subroutine heapsort_real(a)
      real, intent(in out) :: a(0:)
      integer :: start, n, bottom
      real :: temp

      n = size(a)
      do start = (n - 2) / 2, 0, -1
        call shiftdown_real(a, start, n);
      end do

      do bottom = n - 1, 1, -1
        temp = a(0)
        a(0) = a(bottom)
        a(bottom) = temp;
        call shiftdown_real(a, 0, bottom)
      end do

   end subroutine heapsort_real

   subroutine heapsort_int(a)
      integer, intent(in out) :: a(0:)
      integer :: start, n, bottom
      integer :: temp

      n = size(a)
      do start = (n - 2) / 2, 0, -1
        call shiftdown_int(a, start, n);
      end do

      do bottom = n - 1, 1, -1
        temp = a(0)
        a(0) = a(bottom)
        a(bottom) = temp;
        call shiftdown_int(a, 0, bottom)
      end do

   end subroutine heapsort_int

   subroutine shiftdown_real(a, start, bottom)
     real, intent(in out) :: a(0:)
     integer, intent(in) :: start, bottom
     integer :: child, root
     real :: temp

     root = start
     do while(root*2 + 1 < bottom)
       child = root * 2 + 1

       if (child + 1 < bottom) then
         if (a(child) < a(child+1)) child = child + 1
       end if

       if (a(root) < a(child)) then
         temp = a(child)
         a(child) = a (root)
         a(root) = temp
         root = child
       else
         return
       end if
     end do

   end subroutine shiftdown_real

   subroutine shiftdown_int(a, start, bottom)
     integer, intent(in out) :: a(0:)
     integer, intent(in) :: start, bottom
     integer :: child, root
     integer :: temp

     root = start
     do while(root*2 + 1 < bottom)
       child = root * 2 + 1

       if (child + 1 < bottom) then
         if (a(child) < a(child+1)) child = child + 1
       end if

       if (a(root) < a(child)) then
         temp = a(child)
         a(child) = a (root)
         a(root) = temp
         root = child
       else
         return
       end if
     end do

   end subroutine shiftdown_int


   ! Sorts and merges two INTEGER arrays into one, with no repetition
   function MergeArrays(A,B) result (C)
      integer, dimension(:), intent(in) :: A
      integer, dimension(:), intent(in) :: B
      integer, allocatable, dimension(:) :: C
      integer :: NA,NB,NC,N
      integer :: I,J,K
      integer, dimension(:), allocatable :: temp, tempA, tempB

      NA = size(A)
      NB = size(B)
      N = NA + NB

      allocate( tempA(1:NA) , source = A)
      allocate( tempB(1:NB) , source = B)
      allocate( temp(1:N) , source = 0 )

      call heapsort(tempA)
      call heapsort(tempB)

      I = 1; J = 1; K = 2;
      if ( tempA(1) <= tempB(1)) then
         temp(1) = tempA(1)
         I = I + 1
      else
         temp(1) = tempB(1)
         J=J+1
      endif

      do while(I <= NA .and. J <= NB)
         if ( tempA(I) <= tempB(J)) then
            if ( temp(K-1) /= tempA(I) ) then
               temp(K) = tempA(I)
               K = K+1
            endif
            I = I+1
         else
            if ( temp(K-1) /= tempB(J) ) then
               temp(K) = tempB(J)
               K = K+1
            endif
            J = J+1
         endif
!         K = K + 1
      enddo
      do while (I <= NA)
         if ( temp(K-1) /= tempA(I) ) then
               temp(K) = tempA(I)
               K = K+1
         endif
         I = I + 1
      enddo
      do while (J <= NB)
         if ( temp(K-1) /= tempB(J) ) then
               temp(K) = tempB(J)
               K = K+1
         endif
         J = J + 1
      enddo
      NC = K-1

      if (allocated(C)) deallocate(C)
      allocate(C(1:NC))
      C(1:NC) = temp(1:NC)
      deallocate(temp,tempA,tempB)

   end function MergeArrays


end module sort_mod
