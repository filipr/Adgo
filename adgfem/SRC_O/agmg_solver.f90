! agmg library neccessary !
module agmg_solver

use matrix_oper
use mtxform

implicit none

public:: example_seq
public:: SolveViaAGMG

CONTAINS

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

  subroutine SolveViaAGMG(nsize,x,b)
  integer,intent(in) :: nsize
  real,dimension(:),allocatable :: Ax
  integer,dimension(:),allocatable :: Ai, Ap
  real,dimension(:),intent(inout) :: x,b
  integer :: iter,iprint,ijob,nrest
  real :: tol
  !
  !call DG2CROWS(Ap,Ai,Ax)
  call DG2HB(Ap,Ai,Ax)
  !
  !maximal number of iterations
  iter=50
  !
  !tolerance on relative residual norm
  tol=1.e-6
  !
  !unit number for output messages: 6 => standard output
  !nonpositive input suppresses all messegess, but error messagess printed on standardd output
  iprint=6
  !
  !ijob usage or remark
  ! 0 performs setup + solve + memory release, no initial guess
  ! 10 performs setup + solve + memory release, initial guess in x
  ! 1 performs setup only(preprocessing: prepares all parameters for subsequent solves)
  ! 2 solves only (based on previous setup), no initial guess
  ! 12 solves only (based on previous setup), initial guess in x
  ! 3 the vector returned in x is not the solution of the linear system,
  !   but the result of the action of the multigrid
  !   preconditioner on the right hand side in f
  ! -1 erases the setup and releases internal memory
  ! 100,110,101,102,112 same as, respectively, 0,10,1,2,12 but
  !  use the TRANSPOSE of the input matrix
  ! 2,3,12,102,112 require that one has previously called AGMG
  ! with ijob=1 or ijob=101
  ! 1,101 the preconditioner defined when calling AGMG 
  ! is entirely kept in internal memory; hence:
  ! 3 the input matrix is not accessed
  ! 2,12,102,112 the input matrix is used only to perform matrix vector product
  ! within the main iterative solution process
  ijob = 100 !110

  ! NREST restat parameter for GCR (an alternative implementtion of GMRES)
  ! nonpositive value is convered to 10 (suggested dafault)
  ! nrest=1 flexible CG is used instead of GCR, simplification are performed.
  ! For symmetric and positive definite matrix only
  nrest = 10

  call dagmg(nsize,Ax,Ai,Ap,b,x,ijob,iprint,nrest,iter,tol)

  end subroutine

end module agmg_solver
