module AGMGwrap
! This module is dedicated to be used with AGMG by Y.Notay. It contains 
! >example_seq< test subroutine from original source files.

public example_seq

implicit none

! Subroutine from AGMG distr.
 subroutine example_seq
!
!  Solves the discrete Laplacian on the unit square by simple call to agmg.
!  The right-hand-side is such that the exact solution is the vector of all 1.
!
        implicit none
        real ,allocatable :: a(:),f(:),x(:)
        integer,allocatable :: ja(:),ia(:)
        integer :: n,iter,iprint,nhinv
        real  :: tol
!
!       set inverse of the mesh size (feel free to change)
        nhinv=500
!
!       maximal number of iterations
        iter=50
!
!       tolerance on relative residual norm
        tol=1.e-6
!
!       unit number for output messages: 6 => standard output
        iprint=6
!
!       generate the matrix in required format (CSR)
!
!         first allocate the vectors with correct size
            N=(nhinv-1)**2
            allocate (a(5*N),ja(5*N),ia(N+1),f(N),x(N))
!         next call subroutine to set entries
            call uni2d(nhinv-1,f,a,ja,ia)
!
!       call agmg
!         argument 5 (ijob)  is 0 because we want a complete solve
!         argument 7 (nrest) is 1 because we want to use flexible CG
!                            (the matrix is symmetric positive definite)
!
         call dagmg(N,a,ja,ia,f,x,0,iprint,1,iter,tol)
!
!      uncomment the following lines to write solution on disk for checking
!
!       open(10,file='sol.out',form='formatted')
!       write(10,'(e22.15)') f(1:n)
!       close(10)
      END subroutine example_seq
!----------------------------------------------------------------------
      subroutine uni2d(m,f,a,ja,ia)
!
! Fill a matrix in CSR format corresponding to a constant coefficient
! five-point stencil on a square grid
!
      implicit none
      real  :: f(*),a(*)
      integer :: m,ia(*),ja(*)
      integer :: k,l,i,j
      real , parameter :: zero=0,cx=-1.,cy=-1., cd=4.
!
      k=0
      l=0
      ia(1)=1
      do i=1,m
        do j=1,m
          k=k+1
          l=l+1
          a(l)=cd
          ja(l)=k
          f(k)=zero
          if(j < m) then
             l=l+1
             a(l)=cx
             ja(l)=k+1
            else
             f(k)=f(k)-cx
          end if
          if(i < m) then
             l=l+1
             a(l)=cy
             ja(l)=k+m
            else
             f(k)=f(k)-cy
          end if
          if(j > 1) then
             l=l+1
             a(l)=cx
             ja(l)=k-1
            else
             f(k)=f(k)-cx
          end if
          if(i >  1) then
             l=l+1
             a(l)=cy
             ja(l)=k-m
            else
             f(k)=f(k)-cy
          end if
          ia(k+1)=l+1
        end do
      end do
!
      return
      end subroutine uni2D

end module
