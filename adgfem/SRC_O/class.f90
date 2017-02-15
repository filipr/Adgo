module oop

   implicit none

   private

   public :: Mblock_N


   type, public :: Mblock_N
      integer :: n
      real, allocatable, dimension( : , : ) :: Mb

      procedure :: addMblock_N
      procedure :: subMblock_N
      procedure :: assignMblock_N
      procedure :: delMblock_N


      generic, public :: operator(+) => addMblock_N
      generic, public :: operator(-) => subMblock_N
      generic, public :: assignment(=) => assignMblock_N
   end type Mblock_N

   interface Mblock_N
      module procedure initMblock_N
   end interface

   contains

   type(Mblock_N) function initMblock_N(n)
      integer, intent(in) :: n

      allocate ( initMblock_N%Mb(n,n) )
      initMblock_N%Mb(:,:) = 0.
      initMblock_N%n = n
   end function initMblock_N

end module oop
