module pokus_mod

implicit none

   integer :: try_int = 8

   contains

   subroutine realSub(ndimL, Qdof, w, f_s, x )
    integer, intent(in) :: Qdof, ndimL
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:try_int,1:ndimL), intent(inout) :: f_s
    real, dimension(1:Qdof,1 :try_int), intent(in) :: x

    integer, parameter :: ityp = 0   ! function f_s
    integer :: i, k, l

    do l=1,ndimL
       do i=1,Qdof
          do k=1,try_int
             f_s(i, k, l) =5
          enddo
       enddo
    enddo

   end subroutine realSub

   subroutine subWInter( AbstrSubroutine, i )
      interface
         subroutine AbstrSubroutine(ndimL, Qdof, w, f_s, x )
            integer, intent(in) :: Qdof, ndimL
            real, dimension(:,:), intent(in):: w !state  w in #Qdof nodes
            real, dimension(:,:,:), intent(inout) :: f_s
            real, dimension(:,:), intent(in) :: x
         end subroutine
      end interface
      integer, intent(in) :: i

      integer :: k, j,m
      real, dimension(:,:), allocatable :: w
      real, dimension(:,:,:), allocatable :: f_s
      real, dimension(:,:), allocatable :: x


      do k = 1,i
         print*, 'HERE I AM'
      enddo
      j = 4

      allocate( w(1:j,1:j) , f_s(1:j,1:j, 1:j), x(1:j, 1:j) )

      call AbstrSubroutine( j, m ,w(1:2,1:2), f_s(:,:,:), x(:,:) )



   end subroutine subWInter

end module pokus_mod
