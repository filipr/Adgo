!> subroutines for MulritGrid (MG) solvers
module pMultiGrid
  use main_data
  use paramets
  use geometry
  use matrix_oper
  use matrix_oper_int
!  use gmres_solver
!  use lin_solvers
  use kind_parameters
  use sparse_utils  ! module for sparse matrix operations

  implicit none

  type,public :: MGSetting
    integer :: Vcycle = 2
    integer :: presmooth = 1
    integer :: postsmooth = 1
  end type MGSetting

  public:: example_seq


  public:: InitMG
  public:: pMGSolver
  public:: bJacobi
  public:: MGbMVprod
  public:: SolveViaUMFPACK		! direct solver via UMFPACK
  public:: pMGLinSolverRecu
  public:: pMGLinSolverRecu4LINSOL
  private:: LINMGLVL_mod


  ! NOT USED FOR WHILE
  public:: MGcopyElem
  public:: RestrictElem
  private:: NOT_USED_LINMGLVL			! replaced by LINMGLVL_mod


  ! IN DEVELOPMENT
  public:: MGbMViLUprod
  public:: pMGLinSolverRecuILU
  private:: LINMGLVL_mod_ILU

  public:: test_umfpack			! will be *DELETED*
  private:: MGrestVec, MGprolVec_new, MGprolVec

  ! fcn related to matrix
  public:: CompareMatrix
  public:: DajHBij
  public:: GetHBij
  public:: DG2HB
  public:: DG2CROWS 			! --- NOT TESTED YET --- !

  private:: SORT1, PIKSRT, TWRIT
  private:: write_umfpack_info
  private:: XXX 			! for testing *ONLY*

  public:: umfpack_wrapper
  public:: L2Res

  ! TO CHECK:
  public:: pMGSolverRecu



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
!pMGLinSolverRecu4LINSOL
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

  subroutine SolveViaAGMG(b)
  real,dimension(:),allocatable :: Ax
  integer,dimension(:),allocatable :: Ai, Ap
  real,dimension(state%nsize) :: x,y,bb
  real,dimension(state%nsize),intent(in) :: b
  integer :: iter,tol,iprint
  !
  !call DG2CROWS(Ap,Ai,Ax)
  call DG2HB(Ap,Ai,Ax)
  !
  !maximal number of iterations
  iter=100
  !
  !tolerance on relative residual norm
  tol=1.e-6
  !
  !unit number for output messages: 6 => standard output
  iprint=6
  !
bb = b
  call dagmg(state%nsize,Ax,Ai,Ap,b,x,100,iprint,1,iter,tol)
  print*,iter
  print*,norm2(b),norm2(bb)
  call bMVprod(y,x,state%nsize)
  print*, norm2(b-y)/norm2(b),tol
print*, norm2(bb-y)/norm2(bb),tol
pause
  end subroutine


  subroutine InitMG()

    grid%elem(:)%MGdof = grid%elem(:)%dof
    grid%elem(:)%MGdeg = grid%elem(:)%deg
    grid%elem(:)%MGncv = grid%elem(:)%ncv

  end


  subroutine XXX

  integer :: N = 3
  real, dimension(9) :: Ax = (/1, 0, 0, 0, 1, 0, 0, 0, 1/)
  integer, dimension(9) :: Ai = (/1, 2, 3, 1, 2, 3, 1, 2, 3/)
  integer,dimension(4) :: Ap = (/1, 4, 7, 10/)

  !real(kind=8), dimension(3) :: Ax = (/1,1,1/)
  !integer, dimension(3) :: Ai = (/1,2,3/)
  !integer,dimension(3+1) :: Ap = (/1, 2, 3, 4/)

  real(kind=8), dimension(3) :: b = (/1, 2, 3/)
  real(kind=8), dimension(3) :: x

  print*,'Ax',Ax
  print*,'Ai',Ai
  print*,'Ap',Ap
  print*,'b',b
  print*,'x',x

  write (unit=6, fmt=*) "input matrix:"
  call csc1_spy (6, Ap, Ai)          ! 1-based

  ! 0-based
  Ai=Ai-1
  Ap=Ap-1
  call umfpack_wrapper(Ax,Ai,Ap,x,b,1)

  print*, x

  end subroutine XXX

  subroutine SORT1(A,N,iout)
  integer,dimension(:,:),intent(inout) ::    A     !Table to be sorted
  integer,intent(in) :: N,iout

  if (iout.eq.1) then
    print *,' '
    print *,'Table to be sorted:'
    call TWRIT(N,A)
  endif


  call PIKSRT(N,A)

  if (iout.eq.1) then
    print *,' '
    print *,'Sorted table (straight insertion):'
    call TWRIT(N,A)
  endif

  end subroutine SORT1

  ! Sorting an array with straight insertion method *
  subroutine PIKSRT(N,ARR)
  integer,intent(in) :: N
  integer,intent(inout) :: ARR(3,N)
  integer :: i,j,a(3)

  do j=2,N
    a = ARR(1:3,j)
    do i=j-1,1,-1
      if ( ARR(1,i)<=a(1)) goto 10
      ARR(1:3,i+1)=ARR(1:3,i)
    end do
    i=0
10  ARR(1:3,i+1 )= a
  end do
  return

  end subroutine PIKSRT


  !write table of size N to standard output
  subroutine TWRIT(N,ARR)
  integer::N
  integer :: ARR(3,N)
  integer :: i

  write(*,"(I6,I6,I6)") (ARR(1,i),ARR(2,i),ARR(3,i),i=1,N)
  return
  end subroutine TWRIT


  !> Store DG-matrix in Harwell-Boeing (compressed column storage) format
  !> with 1-based indexing
  subroutine DG2HB(Ap,Ai,Ax)
  integer,dimension(:),allocatable,intent(inout) :: Ap, Ai
  real,dimension(:),allocatable,intent(inout) :: Ax
  integer :: i,j,k,shiftidx,last
  integer :: allsize, facesize, ptrsize, ptridx
  type(element),pointer :: elem, elem1
  integer :: dof,ndof,ndof1
  integer,dimension(:,:),allocatable :: array
  real,dimension(:,:),allocatable :: accum
  integer :: iii


  allocate( accum(maxval(grid%elem%MGdof)*ndim, maxval(grid%elem%MGdof)*ndim) )
  allocate( array(3, 1+maxval(grid%elem%flen) ) )

  if( 0.eq.1) then
    print*,'#elem',grid%nelem
    do i=1,grid%nelem,1
      print*,'elem',i
      do j=1,grid%elem(i)%flen,1
        k = grid%elem(i)%face(neigh,j)
        print*,'sused',k
      enddo
    enddo
  end if

  ! evaluate sizes of arrays Ap,Ai,Ax
  allsize = 0
  ptrsize = 1	! first after, ie. we have 0 columns

  do i=1,grid%nelem,1
    elem => grid%elem(i)
    ptrsize = ptrsize + ndim*elem%MGdof
    facesize = elem%MGdof
    do j=1,elem%flen,1
      k = elem%face(neigh,j)
      if (k>0) then
        facesize = facesize + grid%elem(k)%MGdof
      endif
    enddo
    allsize = allsize + ndim**2*elem%MGdof*facesize
  enddo

  ! allocate arrays
  allocate(Ap(ptrsize), Ai(allsize), Ax(allsize))

  Ap=0
  Ai=0
  Ax=0.

  !print*,ptrsize-1,'?=',state%nsize

  ! print*,'fill arrays with data'
  shiftidx = 1
  ptridx = 0

  do i=1,grid%nelem,1
    elem => grid%elem(i)
    dof = elem%MGdof
    ndof = ndim*dof

    ! Array "array" is used to get information about column structure, ie.
    ! which block in column is above and/or bellow mtx diagonal
    ! array(1,:) stores information for faces
    array(1,:) = (/ i, elem%face(neigh, 1:elem%flen) /)
    ! array(2,:) stores information where corresponding flux matrix block
    ! is stored, ie. elem%block(*)
    array(2,:) = (/ (k, k=0,elem%flen,1) /)
    array(3,:) = (/ 0, elem%face(nei_i, 1:elem%flen) /)

    call SORT1(array(:,1:elem%flen+1),elem%flen+1,0)

    !print*,'set data for indices array Ai'
    last = 0;
    do j=1,elem%flen+1,1
      if ( array(1,j).gt.0 ) then
        iii = array(1,j)
        elem1 => grid%elem(iii)
        ndof1 = ndim*elem1%MGdof
        Ai(shiftidx+last : shiftidx+last+ndof1-1) = (/ (k, k=elem1%MGncv,elem1%MGncv+ndof1-1,1) /)
        last = last+ndof1
      endif
    enddo
    !print*,'copy Ai for all 2:ndof columns'
    do k=2,ndof,1
      Ai(shiftidx+(k-1)*last : shiftidx+k*last-1) = Ai(shiftidx : shiftidx+last-1)
    enddo

    ! evaluate diagonal matrix block
    accum(1:ndof,1:ndof) = elem%block(0)%Mb(1:ndof, 1:ndof)

    if ( eta.ne.0. ) then
      do k = 0,ndof-1,dof
        accum(k+1:k+dof,k+1:k+dof) = accum(k+1:k+dof,k+1:k+dof) + &
          eta * elem%Mass%Mb(1:dof, 1:dof)
      enddo
    endif

    !print*,'set data for values array Ax'
    do k=1,ndof,1
      ptridx = ptridx+1
      Ap(ptridx)=shiftidx

      do j=1,elem%flen+1,1

        if ( array(1,j).eq.i ) then ! equivalent (array(2,j)==0)
          !print*,'diagonal block'
          Ax(shiftidx:shiftidx+ndof-1) = accum(1:ndof,k)
          shiftidx = shiftidx+ndof
        elseif ( array(1,j).gt.0 ) then
          !print*,'offdiagonal block in column'
          elem1 => grid%elem( array(1,j) )
          ndof1 = ndim*elem1%MGdof
          Ax(shiftidx:shiftidx+ndof1-1) = grid%elem( array(1,j) )%block( array(3,j) )%Mb(1:ndof1,k)
          shiftidx = shiftidx+ndof1
        endif

      enddo !j=1,4,1
    enddo !k=1,ndof,1pMGLinSolverRecu4LINSOL
  enddo !i=1,grid%nelem,1

  Ap(ptrsize) = allsize+1

  end subroutine DG2HB


  !> Store DG-matrix in COMPRESSED ROW STORAGE format
  !> with 1-based indexing
  subroutine DG2CROWS(Ap,Ai,Ax)
  integer,dimension(:),allocatable,intent(inout) :: Ap, Ai
  real,dimension(:),allocatable,intent(inout) :: Ax
  integer :: i,j,k,shiftidx,last
  integer :: allsize, facesize, ptrsize, ptridx
  type(element),pointer :: elem, elem1
  integer :: dof,ndof,ndof1
  integer,dimension(:,:),allocatable :: array
  real,dimension(:,:),allocatable :: accum
  integer :: iii


  allocate( accum(maxval(grid%elem%MGdof)*ndim, maxval(grid%elem%MGdof)*ndim) )
  allocate( array(3, 1+maxval(grid%elem%flen) ) )

  if( 0.eq.1) then
    print*,'#elem',grid%nelem
    do i=1,grid%nelem,1
      print*,'elem',i
      do j=1,grid%elem(i)%flen,1
        k = grid%elem(i)%face(neigh,j)
        print*,'sused',k
      enddo
    enddo
  end if

  ! evaluate sizes of arrays Ap,Ai,Ax
  allsize = 0
  ptrsize = 1	! first after, ie. we have 0 columns

  do i=1,grid%nelem,1
    elem => grid%elem(i)
    ptrsize = ptrsize + ndim*elem%MGdof
    facesize = elem%MGdof
    do j=1,elem%flen,1
      k = elem%face(neigh,j)
      if (k>0) then
        facesize = facesize + grid%elem(k)%MGdof
      endif
    enddo
    allsize = allsize + ndim**2*elem%MGdof*facesize
  enddo

  ! allocate arrays
  allocate(Ap(ptrsize), Ai(allsize), Ax(allsize))

  Ap=0
  Ai=0
  Ax=0.

  !print*,ptrsize-1,'?=',state%nsize

  ! print*,'fill arrays with data'
  shiftidx = 1
  ptridx = 0

  do i=1,grid%nelem,1
    elem => grid%elem(i)
    dof = elem%MGdof
    ndof = ndim*dof

    ! Array "array" is used to get information about column structure, ie.
    ! which block in column is above and/or bellow mtx diagonal
    ! array(1,:) stores information for faces
    array(1,:) = (/ i, elem%face(neigh, 1:elem%flen) /)
    ! array(2,:) stores information where corresponding flux matrix block
    ! is stored, ie. elem%block(*)
    array(2,:) = (/ (k, k=0,elem%flen,1) /)
    array(3,:) = (/ 0, elem%face(nei_i, 1:elem%flen) /)

    call SORT1(array(:,1:elem%flen+1),elem%flen+1,0)

    !print*,'set data for indices array Ai'
    last = 0;
    do j=1,elem%flen+1,1
      if ( array(1,j).gt.0 ) then
        iii = array(1,j)
        elem1 => grid%elem(iii)
        ndof1 = ndim*elem1%MGdof
        Ai(shiftidx+last : shiftidx+last+ndof1-1) = (/ (k, k=elem1%MGncv,elem1%MGncv+ndof1-1,1) /)
        last = last+ndof1
      endif
    enddo
    !print*,'copy Ai for all 2:ndof columns'
    do k=2,ndof,1
      Ai(shiftidx+(k-1)*last : shiftidx+k*last-1) = Ai(shiftidx : shiftidx+last-1)
    enddo

    ! evaluate diagonal matrix block
    accum(1:ndof,1:ndof) = elem%block(0)%Mb(1:ndof, 1:ndof)

    if ( eta.ne.0. ) then
      do k = 0,ndof-1,dof
        accum(k+1:k+dof,k+1:k+dof) = accum(k+1:k+dof,k+1:k+dof) + &
          eta * elem%Mass%Mb(1:dof, 1:dof)
      enddo
    endif

    !print*,'set data for values array Ax'
    do k=1,ndof,1
      ptridx = ptridx+1
      Ap(ptridx)=shiftidx

      do j=1,elem%flen+1,1


        if ( array(1,j).eq.i ) then ! equivalent (array(2,j)==0)
          !print*,'diagonal block'
          Ax(shiftidx:shiftidx+ndof-1) = accum(k,1:ndof)
          shiftidx = shiftidx+ndof
        elseif ( array(1,j).gt.0 ) then
          !print*,'offdiagonal block in column'
          elem1 => grid%elem( array(1,j) )
          ndof1 = ndim*elem1%MGdof
          Ax(shiftidx:shiftidx+ndof1-1) = grid%elem( array(1,j) )%block( array(2,j) )%Mb(k,1:ndof1)
          shiftidx = shiftidx+ndof1
        endif

      enddo !j=1,4,1 pMGLinSolverRecu4LINSOL
    enddo !k=1,ndof,1
  enddo !i=1,grid%nelem,1

  Ap(ptrsize) = allsize+1

  end subroutine DG2CROWS


  subroutine write_umfpack_info(text,info)
    real(kind=double),dimension(90),intent(in) :: info
    character(len=*),intent(in) :: text

    write(*,fmt="(a,            /,  &
                 &a12, f5.0,      /,  &
                 &a12, es10.2,a6, /,  &
                 &a42,            /,  &
                 &a18, f10.2,a5,  /,  &
                 &a18, f10.2,a5,  /,  &
                 &a18, es10.2,    /,  &
                 &a18, f10.0,     /,  &
                 &a18, f10.0)")       & ! END OF FORMAT
   text,                  & ! Data to print from here
   "   status:  ",info(1) ,               &
   "   time:    ",info(16)," (sec)",      &
   "   estimates (upper bound) for numeric LU:",           &
   "   size of LU:    ",(info(21)*info(4))/2**20," (MB)",  &
   "   memory needed: ",(info(22)*info(4))/2**20," (MB)",  &
   "   flop count:    ", info(23),        &
   "   nnz (L):       ", info(24),        &
   "   nnz (U):       ", info(25)
  end subroutine write_umfpack_info





  !> Direct solver via UMFPACK subroutine
  !> input:
  !> N-by-N matrix in Harwell-Boeing format 0-based indexing, i.e. [Ax],[Ai],[Ap]
  !> right hand side [b]
  !> output:
  !> vector [x] such, that Ax=b
  subroutine umfpack_wrapper(Ax,Ai,Ap,x,b,iout)
  ! input: NxN matrix in Harwell/Boeing format 0-base indexing, solution, RHS
  real(kind=8),dimension(:),intent(in)    :: Ax
  integer,dimension(:),intent(in)         :: Ai, Ap
  real(kind=8),dimension(:),intent(in)    :: b
  real(kind=8),dimension(:),intent(inout) :: x
  integer                                 :: N
  ! umfpack control variables
  real(kind=8),dimension(20)              :: control
  real(kind=8),dimension(90)              :: info
  integer                                 :: numeric,symbolic,status,sys,filenum

  integer :: iout

  ! file numbering for symbolic analysis and numerical factorization, i.e.
  ! s<filenum>.umf and n<filenum>.umf, respectively
  filenum = 0

  N = size(Ap) -1

  call umf4def(control)

  ! print control parameters. set control (1) to 1 to print
  ! error messages only
  !control(1) = 2
  control(1) = 1
  call umf4pcon(control)

  ! pre-order and symbolic analysis
  call umf4sym(N, N, Ap, Ai, Ax, symbolic, control, info)

  ! check umf4sym error condition
  if (info(1) < 0) then
    print *, "Error occurred in umf4sym: ", info(1)
    stop
  else
    if(iout.eq.1) call write_umfpack_info("symbolic analysis:",info)
  endif

  ! numeric factorization
  call umf4num(Ap, Ai, Ax, symbolic, numeric, control, info)

  ! check umf4num error condition
  if (info(1) < 0) then
    print *, "Error occurred in umf4num: ", info(1)
    stop
  else
    if(iout.eq.1) call write_umfpack_info("numeric factorization:",info)
  endif

  ! save the symbolic analysis to the file s<filenum>.umf
  ! note that this is not needed until another matrix is
  ! factorized, below.
  call umf4ssym(symbolic, filenum, status)

  if (status < 0) then
    print *, "Error occurred in umf4ssym: ", status
    stop
  endif

  ! save the LU factors to the file n<filenum>.umf
  call umf4snum(numeric, filenum, status)
  if (status < 0) then
    print *, "Error occurred in umf4snum: ", status
    stop
  endif

  ! free the symbolic analysis
  call umf4fsym(symbolic)

  ! free the numeric factorization
  call umf4fnum(numeric)

  ! load the numeric factorization back in (filename: n<filenum>.umf)
  call umf4lnum(numeric, filenum, status)
  if (status < 0) then
    print *, "Error occurred in umf4lnum: ", status
    stop
  endif

  ! solve Ax=b, without iterative refinement
  sys = 0
  call umf4sol(sys, x, b, numeric, control, info)
  if (info(1) < 0) then
    print *, "Error occurred in umf4sol: ", info(1)
    stop
  else
    if(iout.eq.1) call write_umfpack_info("solving Ax=b:",info)
  endif

  ! free the numeric factorizationpMGLinSolverRecu4LINSOL
  call umf4fnum(numeric)

  ! No LU factors (symbolic or numeric) are in memory at this point.

  ! print final statistics
  call umf4pinf(control, info)
  if(iout.eq.1) call write_umfpack_info("final statistic:",info)

  end subroutine umfpack_wrapper


  function L2Res(x,b,nsize,prod)
    real :: L2Res
    ! input
    integer,intent(in) :: nsize
    real,dimension(1:nsize),intent(in) :: x	! solution
    real,dimension(1:nsize),intent(in) :: b 	! RHS
    ! local variable
    real,dimension(1:nsize) :: Ax		! Ax = A*x, approximation of RHS

    interface
      subroutine prod(b,x,nsize)
        integer, intent (in):: nsize
        real, dimension(1:nsize), intent(in) :: x
        real, dimension(1:nsize), intent(inout) :: b
      end subroutine prod
     end interface

     call prod(Ax,x,nsize)

     ! res = b-A*x = b-Ax
     L2Res = norm2(b-Ax)

  end function L2Res



  subroutine SolveViaUMFPACK(N,x,b)
  integer,intent(in):: N
  real,dimension(state%nsize),intent(out) :: x
  real,dimension(state%nsize),intent(in) :: b
  integer,dimension(:),allocatable :: Ap,Ai
  real(kind=8),dimension(:),allocatable :: Ax

  integer,parameter :: iout = 0
  integer:: iwrite=6

!print*,'@@@####### u1',grid%elem(5)%face(idx,:)
  call DG2HB(Ap,Ai,Ax)
!print*,'@@@####### u2',grid%elem(5)%face(idx,:)
  ! -- Make sense only for small matrices
!open(iwrite,file='matica.txt')
!  write (unit=iwrite, fmt=*) "input matrix:"
!  call csc1_spy (iwrite, Ap, Ai)          ! 1-based
!close(iwrite)
!pause
  ! 0-based
  Ap = Ap-1
  Ai = Ai-1

!print*,N,size(Ap),size(Ai),size(Ax),state%nonzero,state%nsize
!print*,'Ax',Ax(1:3)
!print*,'DG',grid%elem(1)%block(0)%Mb(1:3,1)+eta*grid%elem(1)%Mass%Mb(1:3,1)
!pause

  call umfpack_wrapper(Ax,Ai,Ap,x,b, iout)


  end subroutine SolveViaUMFPACK


  ! -- ONLY FOR TESTING PURPOUSES
  subroutine test_umfpack(yy,bb)
  ! see UMFPACK manual for correct data setting
  integer :: iread=5                  ! File-unit numbers
  integer :: iwrite=6                 ! File-unit numbers
  integer :: fid
  integer :: ios
  integer :: numeric, symbolic, status, sys, filenum ! UMFPACK-variables
  integer,dimension(:),allocatable :: Ap,Ai,Bp,Bi
  real(kind=8),dimension(:),allocatable :: Ax
  real(kind=8),dimension(:),allocatable :: Bx
  real(kind=8),dimension(20)            :: control
  real(kind=8),dimension(90)            :: info
  character(len=72)                          :: title
  character(len=8)			     :: key
  real,dimension(state%nsize),intent(out) :: yy
  real,dimension(state%nsize),intent(in) :: bb
  real,dimension(state%nsize) :: rhs
  integer :: N
  CHARACTER(LEN=20), PARAMETER :: FMT = "(I10)"
  integer :: NZ
  real,dimension(3):: y,b
  integer,parameter :: iout = 1;


  grid%elem(:)%MGdof = grid%elem(:)%dof
  grid%elem(:)%MGncv = grid%elem(:)%ncv

  call DG2HB(Ap,Ai,Ax)

  if( 1.eq.0 ) then
    N = 3
    NZ = 4

    deallocate(Ap,Ai,Ax)
    allocate(Ax(NZ),Ai(NZ),Ap(N+1))
    Ax = (/ 1.,1.,sqrt(2.),1. /)
    Ai=(/ 1,2,1,3 /)
    Ap=(/ 1,2,3,5 /)

    y(:) = 0.
    b = (/4, 2, 1/)
  endif

  write(6,FMT) size(Ap)
  write(6,FMT) size(Ai)
  write(6,FMT) size(Ax)
  write(6,FMT) N
  write(6,FMT) Ap(N)
  print*,'N',N

  write (unit=iwrite, fmt=*) "input matrix:"
  call csc1_spy (iwrite, Ap, Ai)          ! 1-based

pause


  fid = 100

  open (unit=fid, file='a.rua',iostat=ios)
  if (ios.eq.0) then
    title = 'DG-matrix'
    key = '--'
    call rua_write(fid,Ap,Ai,Ax,title,key)
    close (unit=fid)
  else
    print*,'!!! Problem in writting Harwell/Boeing format into file'
    stop
  endif

!  call WriteMatrixA(state%eta(1))

  allocate(Bi(size(Ai)), Bp(size(Ap)), Bx(size(Ax)))

  if(1.eq.1) then
    open (unit=fid, file='a.rua',iostat=ios)
    if (ios.eq.0) then
      title = 'DG-matrix'
      key = '--'
      call rua_read(fid,Bp,Bi,Bx,title)
      close(unit=fid)
    else
      print*,'!!! Problem in reading Harwell/Boeing format from file'
      stop
    endif
  endif


  if(all(Ai.eq.Bi)) then
    print*,'Inidices OK'
  else
    print*,'Inidices WRONG'
  endif

  if(all(Ax.eq.Bx)) then
    print*,'Values OK'
  else
    print*,'Values WRONG'
  endif

  if(all(Ap.eq.Bp)) then
    print*,'Pointers OK'
  else
    print*,'Pointers WRONG'
  endif

  Ap = Ap-1
  Ai = Ai-1

if( 1.eq.0 ) then
  write(*,'(/A/)'),'@@@ Varianta 1 @@@'

! convert from 1-based to 0-based indexing

  call umf4def(control)

  ! print control parameters. set control (1) to 1 to print
  ! error messages only
  !control(1) = 2
  control(1) = 1
  call umf4pcon(control)

  ! pre-order and symbolic analysis
  call umf4sym(N, N, Ap, Ai, Ax, symbolic, control, info)

  ! check umf4sym error condition
  if (info(1) < 0) then
    print *, "Error occurred in umf4sym: ", info(1)
    stop
  else
    call write_umfpack_info("symbolic analysis:",info)
  endif

  ! numeric factorization
  call umf4num(Ap, Ai, Ax, symbolic, numeric, control, info)

  ! check umf4num error condition
  if (info(1) < 0) then
    print *, "Error occurred in umf4num: ", info(1)
    stop
  else
    call write_umfpack_info("numeric factorization:",info)
  endif

  ! save the symbolic analysis to the file s0.umf
  ! note that this is not needed until another matrix is
  ! factorized, below.
  filenum = 0
  call umf4ssym(symbolic, filenum, status)
  if (status < 0) then
    print *, "Error occurred in umf4ssym: ", status
    stop
  endif

  ! save the LU factors to the file n0.umf
  call umf4snum(numeric, filenum, status)
  if (status < 0) then
    print *, "Error occurred in umf4snum: ", status
    stop
  endif

  ! free the symbolic analysis
  call umf4fsym(symbolic)

  ! free the numeric factorization
  call umf4fnum(numeric)

  ! load the numeric factorization back in (filename: n0.umf)
  call umf4lnum(numeric, filenum, status)
  if (status < 0) then
    print *, "Error occurred in umf4lnum: ", status
    stop
  endif

  ! solve Ax=b, without iterative refinement
  sys = 0
  call umf4sol(sys, y, b, numeric, control, info)
  if (info(1) < 0) then
    print *, "Error occurred in umf4sol: ", info(1)
    stop
  else
    call write_umfpack_info("solving Ax=b:",info)
  endif
print*,'@@@',numeric
  ! free the numeric factorization
  call umf4fnum(numeric)

  ! No LU factors (symbolic or numeric) are in memory at this point.

  ! print final statistics
  call umf4pinf(control, info)
  call write_umfpack_info("final statistic:",info)

endif



pause('MOJ')

print*,'Ai',Ai
print*,'Ap',Ap
!print*,'Ax',Ax

  call umfpack_wrapper(Ax,Ai,Ap,yy,bb,iout)
pause('MOJ KONIEC')



if (1.eq.1) then
write(*,'(/A/)'),'@@@ Varianta 2 @@@'
Bp = Bp-1
Bi = Bi-1
! varianta 2

  call umf4def(control)

  control(1) = 1
  call umf4pcon(control)

  ! pre-order and symbolic analysis
  call umf4sym(N, N, Bp, Bi, Bx, symbolic, control, info)

  ! check umf4sym error condition
  if (info(1) < 0) then
    print *, "Error occurred in umf4sym: ", info(1)
    stop
  else
    call write_umfpack_info("symbolic analysis:",info)
  endif

  ! numeric factorization
  call umf4num(Bp, Bi, Bx, symbolic, numeric, control, info)
print*,'@@@',numeric
  ! check umf4num error condition
  if (info(1) < 0) then
    print *, "Error occurred in umf4num: ", info(1)
    stop
  else
    call write_umfpack_info("numeric factorization:",info)
  endif

  ! save the symbolic analysis to the file s0.umf
  ! note that this is not needed until another matrix is
  ! factorized, below.
  filenum = 100

  call umf4ssym(symbolic, filenum, status)

  if (status < 0) then
    print *, "Error occurred in umf4ssym: ", status
    stop
  endif

  ! save the LU factors to the file n0.umf
  call umf4snum(numeric, filenum, status)
  if (status < 0) then
    print *, "Error occurred in umf4snum: ", status
    stop
  endif
print*,'@@@',numeric
  ! free the symbolic analysis
  call umf4fsym(symbolic)

  ! free the numeric factorization
  call umf4fnum(numeric)
print*,'@@@',numeric
  ! load the numeric factorization back in (filename: n0.umf)
  call umf4lnum(numeric, filenum, status)
  if (status < 0) then
    print *, "Error occurred in umf4lnum: ", status
    stop
  endif
print*,'@@@_a',numeric
  ! solve Ax=b, without iterative refinement
  sys = 0
  call umf4sol(sys, yy, bb, numeric, control, info)
  if (info(1) < 0) then
    print *, "Error occurred in umf4sol: ", info(1)
    stop
  else
    call write_umfpack_info("solving Ax=b:",info)
  endif

print*,'@@@_o',numeric
  ! free the numeric factorization
  call umf4fnum(numeric)

  ! No LU factors (symbolic or numeric) are in memory at this point.

  ! print final statistics
  call umf4pinf(control, info)
  call write_umfpack_info("final statistics:",info)

endif

call MGbMVprod(rhs,yy,state%nsize)
print*,' '
print*,'***********************'
print*,'norm = ', norm2(bb-rhs)

  end subroutine test_umfpack

  subroutine pMGLinSolverRecuILU(nsize, eta, b, x, rezid, tot_iter, &
         precond_update, not_conv,exact_solu)
    use matrix_oper_int, mx_eta => eta
    !use lin_solvers, only: ComputeBlockILUPrecond
    integer :: nsize                               ! size of the algebraic problem
    real, intent(in):: eta
    real, dimension(1:nsize), intent(inout):: x    ! solution
    real, dimension(1:nsize), intent(inout):: b    ! RHS
    real, intent(inout) :: rezid                   ! reziduum
    integer, intent(inout) :: tot_iter             ! number of iterations
    logical, intent(in) :: precond_update          ! = .false. preconditioner is not update
    integer, intent(inout) :: not_conv             ! convergency

    character(len=*),intent(in):: exact_solu       ! GMRES, bJacobi

    real,dimension(1:nsize):: rr
    !character(len=1) :: precond  ! type of preconditioner: ' ', 'D', 'L', 'J'
    !integer:: iout, iter
    !real :: size(-1:3)
    !integer:: restart = 30 !30    ! 50  ! GMRES restarted after 'restart' iterations  !45
    !integer:: nloops = 5  !5     ! 10   ! maximal number of restarted cycles !100	  !40
    !real :: t0, t1, t2
    integer  :: max_deg
    integer,parameter :: min_deg = 1
    integer,parameter :: presmooth = 1
    integer,parameter :: postsmooth = 1
    integer,parameter :: Vcycle = 2

    max_deg = maxval( grid%elem(:)%deg)

    grid%elem(:)%MGdof = grid%elem(:)%dof
    grid%elem(:)%MGncv = grid%elem(:)%ncv


    !call ComputeBlockILUPrecond( )


!print*,'### ZACIATOK LINMGLVL', max_deg,min_deg,state%space%adaptation
!pause('ZACIATOK LINMGLVL')

    !call LINMGLVL(x, b, nsize, max_deg, min_deg, Vcycle, presmooth, postsmooth,exact_solu)
    call LINMGLVL_mod_ILU(x, b, nsize, max_deg, min_deg, Vcycle, presmooth, postsmooth,exact_solu,(/grid%elem(:)%deg/))

!print*,'### KONIEC LINMGLVL', max_deg,min_deg,state%space%adaptation
!pause('KONIEC LINMGLVL')

    call bMVprod(rr,x,nsize)    ! rr = Ax
    rr(1:nsize) = rr(1:nsize) - b(1:nsize)                ! rr = Ax - b

    rezid = (dot_product( rr(1:nsize), rr(1:nsize) ) )**0.5

  end subroutine pMGLinSolverRecuILU


    subroutine pMGLinSolverRecu4LINSOL(nsize, eta, b, x, max_deg,exact_solu,rezid)

    use matrix_oper_int, mx_eta => eta
    !use lin_solvers, only: ComputeBlockDiagPrecond
    integer :: nsize                               ! size of the algebraic problem
    real, intent(in):: eta
    real, dimension(1:nsize), intent(inout):: x    ! solution
    real, dimension(1:nsize), intent(inout):: b    ! RHS

    real, intent(inout) :: rezid                   ! reziduum

    character(len=*),intent(in):: exact_solu       ! GMRES, bJacobi

    real,dimension(1:nsize):: rr
    !character(len=1) :: precond  ! type of preconditioner: ' ', 'D', 'L', 'J'
    !integer:: iout, iterhttp://ona.idnes.cz/diskuse.aspx?iddiskuse=A140324_110857_dieta_pet
    !real :: size(-1:3)
    !integer:: restart = 30 !30    ! 50  ! GMRES restarted after 'restart' iterations  !45
    !integer:: nloops = 5  !5     ! 10   ! maximal number of restarted cycles !100	  !40
    !real :: t0, t1, t2
    integer  :: max_deg
    integer,parameter :: min_deg = 1
    integer,parameter :: presmooth = 1
    integer,parameter :: postsmooth = 1
    integer,parameter :: Vcycle = 1


!print*,'### ZACIATOK LINMGLVL', max_deg,min_deg,state%space%adaptation
!pause('ZACIATOK LINMGLVL')

    !call LINMGLVL(x, b, nsize, max_deg, min_deg, Vcycle, presmooth, postsmooth,exact_solu)
    call LINMGLVL_mod(x, b, nsize, max_deg, min_deg, Vcycle, presmooth, postsmooth,exact_solu,(/grid%elem(:)%deg/))

!print*,'### KONIEC LINMGLVL', max_deg,min_deg,state%space%adaptation
!pause('KONIEC LINMGLVL')

    call bMVprod(rr,x,nsize)    ! rr = Ax
    rr(1:nsize) = rr(1:nsize) - b(1:nsize)                ! rr = Ax - b

    rezid = (dot_product( rr(1:nsize), rr(1:nsize) ) )**0.5

  end subroutine pMGLinSolverRecu4LINSOL


  subroutine pMGLinSolverRecu(nsize, eta, b, x, rezid, tot_iter, &
         precond_update, not_conv,exact_solu)
    use matrix_oper_int, mx_eta => eta
    !use lin_solvers, only: ComputeBlockDiagPrecond
    integer :: nsize                               ! size of the algebraic problem
    real, intent(in):: eta
    real, dimension(1:nsize), intent(inout):: x    ! solution
    real, dimension(1:nsize), intent(inout):: b    ! RHS
    real, intent(inout) :: rezid                   ! reziduum
    integer, intent(inout) :: tot_iter             ! number of iterations
    logical, intent(in) :: precond_update          ! = .false. preconditioner is not update
    integer, intent(inout) :: not_conv             ! convergency

    character(len=*),intent(in):: exact_solu       ! GMRES, bJacobi

    real,dimension(1:nsize):: rr
    !character(len=1) :: precond  ! type of preconditioner: ' ', 'D', 'L', 'J'
    !integer:: iout, iterhttp://ona.idnes.cz/diskuse.aspx?iddiskuse=A140324_110857_dieta_pet
    !real :: size(-1:3)
    !integer:: restart = 30 !30    ! 50  ! GMRES restarted after 'restart' iterations  !45
    !integer:: nloops = 5  !5     ! 10   ! maximal number of restarted cycles !100	  !40
    !real :: t0, t1, t2
    integer  :: max_deg
    integer,parameter :: min_deg = 1
    integer,parameter :: presmooth = 1
    integer,parameter :: postsmooth = 1
    integer,parameter :: Vcycle = 1

    max_deg = maxval( grid%elem(:)%deg)

    grid%elem(:)%MGdof = grid%elem(:)%dof
    grid%elem(:)%MGncv = grid%elem(:)%ncv


!print*,'### ZACIATOK LINMGLVL', max_deg,min_deg,state%space%adaptation
!pause('ZACIATOK LINMGLVL')

    !call LINMGLVL(x, b, nsize, max_deg, min_deg, Vcycle, presmooth, postsmooth,exact_solu)
    call LINMGLVL_mod(x, b, nsize, max_deg, min_deg, Vcycle, presmooth, postsmooth,exact_solu,(/grid%elem(:)%deg/))

!print*,'### KONIEC LINMGLVL', max_deg,min_deg,state%space%adaptation
!pause('KONIEC LINMGLVL')

    call bMVprod(rr,x,nsize)    ! rr = Ax
    rr(1:nsize) = rr(1:nsize) - b(1:nsize)                ! rr = Ax - b

    rezid = (dot_product( rr(1:nsize), rr(1:nsize) ) )**0.5

  end subroutine pMGLinSolverRecu



  subroutine pMGSolverRecu()
    integer :: i
    type(Newton_type), pointer :: Newton
    Newton => state%nlSolver

    !Newton%b = (/(2*i,i=1,state%nsize,1)/)
    Newton%x = 0.
    print*,' '
    print*,' Newton-Multigrid procedure'
    print*,' '

    ! simple init of MG structure part
    grid%elem(:)%MGdof = grid%elem(:)%dof
    grid%elem(:)%MGncv = grid%elem(:)%ncv
    call NOT_USED_LINMGLVL(Newton%x, Newton%b,state%nsize,5,1,1,1,1,'GMRES')
  end


  ! Get element of DG matrix on position [i,j]
  !!!!!! wihtout Mass matrix !!!!!!
  function GetDGij(i,j)
    real :: GetDGij
    integer, intent(in) :: i,j
    type(element),pointer :: elemi,elemj
    integer :: cumsum,k

    if( (min(i,j)<1).or.(max(i,j)>state%nsize) ) then
      print*,'Wrong index - GetDJij'
      STOP
    end if

    elemi => null()
    elemj => null()

    ! find element for row "i"
    ! find element for column "j"
    cumsum = 0
    do k=1,grid%nelem,1

      cumsum = cumsum + ndim*grid%elem(k)%dof

      if( .not.associated(elemi) ) then
        if( cumsum.ge.i) then
          elemi => grid%elem(k)
          !print*,'I: elem,i,j',k,i,j
          if( associated(elemj) ) EXIT
        end if
      end if

      if( .not.associated(elemj) ) then
        if( cumsum.ge.j) then
          elemj => grid%elem(k)
          !print*,'J: elem,i,j',k,i,j
          if( associated(elemi) ) EXIT
        end if
      end if

    end do

    if( elemi%i.eq.elemj%i) then ! diagonal matrix block
        ! mass matrix

        ! flux matrix
        GetDGij = elemi%block(0)%Mb( 1+i-elemi%ncv , 1+j-elemj%ncv )
    elseif( any( elemi%face(neigh,1:elemi%flen).eq.elemj%i ) ) then ! nonzero off-diagonal matrix block
        do k=1,elemi%flen,1
          if( elemi%face(neigh,k).eq.elemj%i) then
            EXIT
          end if
        end do
        GetDGij = elemi%block(k)%Mb( 1+i-elemi%ncv , 1+j-elemj%ncv )
    else ! zero matrix block
      GetDGij = 0
    end if

  end function GetDGij


  subroutine CompareMatrix
  integer :: i,j
  integer,allocatable,dimension(:) :: Ap
  integer,allocatable,dimension(:) :: Ai
  real,allocatable,dimension(:) :: Ax
  integer :: N


  grid%elem(:)%MGdeg = grid%elem(:)%deg
  grid%elem(:)%MGdof = grid%elem(:)%dof
  grid%elem(:)%MGncv = grid%elem(:)%ncv

  call DG2hb(Ap,Ai,Ax)
print*,'@@ eta',eta,state%time%tau(1)

print*,'Ax',Ax(1:3)
print*,'C=',grid%elem(1)%block(0)%Mb(1:3,1)
print*,'M=',grid%elem(1)%Mass%Mb(1:3,1), maxval(  grid%elem(1)%Mass%Mb )
pause
    do i=1,state%nsize,1
      do j=1,state%nsize,1
        print*,i,j,state%nsize
        if( GetDGij(i,j).eq.GetHBij(Ap,Ai,Ax,i,j) ) then

        else
          print*,'PROBLEM',i,j
          PAUSE
        end if
      end do
    end do
  end subroutine CompareMatrix


  function DajHBij(i,j)
  real :: DajHBij
  integer,intent(in) :: i,j
  integer,allocatable,dimension(:) :: Ap
  integer,allocatable,dimension(:) :: Ai
  real,allocatable,dimension(:) :: Ax
  integer :: N

  grid%elem(:)%MGdeg = grid%elem(:)%deg
  grid%elem(:)%MGdof = grid%elem(:)%dof
  grid%elem(:)%MGncv = grid%elem(:)%ncv

  call DG2HB(Ap,Ai,Ax)

  DajHBij = GetHBij(Ap,Ai,Ax,i,j)

  end function

  ! Get element of HB matrix on position [i,j]
  ! HB = Harwell-Boeing matrix is described by arrays: Ap,Ai,Ax
  function GetHBij(Ap,Ai,Ax,i,j)
    real :: GetHBij
    integer,dimension(:),intent(in) :: Ap,Ai
    real,dimension(:),intent(in) :: Ax
    integer,intent(in) :: i,j
    integer :: k, nzidx

    nzidx = -1

!    print*,'i=',i
!    print*,Ai(Ap(j):Ap(j+1)-1)


    do k=Ap(j),Ap(j+1)-1,1
      if( Ai(k).eq.i) then
        nzidx = k
        EXIT
      end if
    end do

    if( nzidx.gt.-1) then
      GetHBij = Ax( nzidx )
    else
      GetHBij = 0
    end if

  end function GetHbij


! * NOT USED * !
  recursive subroutine NOT_USED_LINMGLVL(x,b,nsize,lvl,lowest,ncycles,presmooth, postsmooth,exact_solu)
    use gmres_solver

    integer,intent(in) :: nsize,presmooth,postsmooth,lowest,lvl, ncycles
    real, dimension(:) :: b
    real, dimension(:), intent(inout) :: x

    class(element), pointer :: elem
    real :: Tstart, Tend
    integer :: i,restsize,prolsize
    real, dimension(:),allocatable :: y,cor,res,prol,u,v

    character(len=7),intent(in)::exact_solu

    ! GMRES
    real :: rezid
    integer :: not_conv
    integer :: restart = 30
    integer :: nloops = 50
    integer :: iter

    integer :: iout = 1

    if (nsize /= size(x)) then
      print*,'ERROR in LINMGLVL. Different vector size as expected'
      STOP
    end if


    allocate( y(1:nsize) )
    y(1:nsize) = 0.

    if (lvl>lowest) then      ! RESTRICTION

      !print*,'####',lvl,lowest

      ! PRESMOOTHING
      if (iout.eq.1) then
        print*,'--- PRESMOOTHING RESIDUUM ---  on level',lvl,'\ '
        write(*,'(A6,I2)') 'LEVEL:',lvl
      endif

      do i=1,presmooth,1
        call bJacobi(x,b,nsize,y)
        x = y
        call MGbMVprod(y,x,nsize)
        if( iout.gt.1) write(*,"(I3,F41.15)") i,norm2(b-y)
      end do

      ! Niekolko dalsich iteracii Jacobiho metody
      allocate(u(1:nsize))
      allocate(v(1:nsize))
      u = x
      do i=1,9,1
        call bJacobi(u,b,nsize,v)
        u = v
        call MGbMVprod(v,u,nsize)
        if( iout.eq.1)  write(*,"(A1,I2,A1,F40.15)") '(',i+presmooth,')',norm2(b-v)
      end do
      deallocate(u,v)

      ! NEW LEVEL

      allocate(cor(1:nsize))

      call MGrestVec(b-y, nsize, cor(1:nsize), restsize)

      allocate(res(1:restsize))
      res(1:restsize) = cor(1:restsize)
      cor(1:restsize) = 0.

!      write(*,"(10F10.5)") cor(1:10)

      if (lvl>lowest+1) then
        do i=1,ncycles,1
          call NOT_USED_LINMGLVL(cor(1:restsize),res(1:restsize),restsize, &
                        lvl-1,lowest,ncycles,presmooth,postsmooth, exact_solu)
        enddo
      else
        call NOT_USED_LINMGLVL(cor(1:restsize),res(1:restsize),restsize, &
                      lvl-1,lowest,ncycles,presmooth,postsmooth, exact_solu)
      endif

      allocate(prol(1:nsize))

      prolsize = nsize
      call MGProlVec(cor(1:restsize),restsize,prol(1:prolsize),prolsize)


      if( iout.gt.1) print*,'Malo by byt rovnake?',prolsize,nsize
      if( iout.eq.1) print*,'--- POSTSMOOTHING RESIDUUM --- on level',lvl,'/'

      if (1.eq.1) then
        ! POSTSMOOTHING
        if( iout.gt.1) print*,x(1:10)
        x(1:nsize) = x(1:nsize) + prol(1:prolsize)
        if( iout.gt.1) print*,x(1:10)

        do i=1,postsmooth,1
          call bJacobi(x,b,nsize,y)
          x = y
          if( iout.gt.1) print*,x(1:10)
          call MGbMVprod(y,x,nsize)
          if( iout.eq.1) write(*,"(I3,F41.15)") i,norm2(b-y)
        end do
      else
        ! POSTSMOOTHING
      allocate(u(1:nsize))
        do i=1,postsmooth,1
          call bJacobi(prol(1:nsize),b-y,nsize,u)
          prol(1:nsize) = u
          call MGbMVprod(u,prol(1:nsize),nsize)
          write(*,"(I3,F41.15)") i,norm2(b-y-u)
        end do
        x(1:nsize) = x(1:nsize) + prol(1:prolsize)
        call MGbMVprod(y,x,nsize)
        write(*,"(F41.15)") norm2(b-y)
      end if

    else ! lvl == lowest
        if( iout.eq.1) print*,'--- "EXACT" RESIDUUM ---       on level',lvl,'V'
        if( iout.eq.1) write(*,"(A6,I2)") 'LEVEL:',lvl

!	call gmres(nsize, x, b, restart*nloops, state%linSolver%tol,  &
!            bMVprod, bMViLUprod, restart,  state%linSolver%tol, iout, iter, rezid, &
!            not_conv)

        if( exact_solu.eq.'GMRES' ) then
          call gmres(nsize, x(1:nsize), b(1:nsize), restart*nloops, state%linSolver%tol, &
                   MGbMVprod, MGbMViLUprod, 30,  state%linSolver%tol, iout, iter, rezid, &
                   not_conv)
        elseif(exact_solu.eq.'bJacobi') then
          do i=1,100,1
            call bJacobi(x,b,nsize,y)
             x = y
            call MGbMVprod(y,x,nsize)
            if ( (i==1).or.(i==100) ) then
              if( iout.eq.1) write(*,"(I5,F38.15)") i,norm2(b-y)
            end if
          end do
       elseif(exact_solu.eq.'bJacobi10') then
          do i=1,10,1
            call bJacobi(x,b,nsize,y)
             x = y
            call MGbMVprod(y,x,nsize)
            if ( (i==1).or.(i==10) ) then
              if( iout.eq.1) write(*,"(I5,F38.15)") i,norm2(b-y)
            end if
          end do
      elseif(exact_solu.eq.'bJacobi1') then
          do i=1,1,1
            call bJacobi(x,b,nsize,y)
             x = y
            call MGbMVprod(y,x,nsize)
            if ( (i==1).or.(i==1) ) then
              if( iout.eq.1) write(*,"(I5,F38.15)") i,norm2(b-y)
            end if
          end do
        else
          print*,'UNKNOWN solver type in LINMGLVL - STOP'
          stop
        endif

    endif


    deallocate(y)
    if ( allocated(cor) ) then
      deallocate(cor,res)
    end if

  end subroutine NOT_USED_LINMGLVL


recursive subroutine LINMGLVL_mod(x,b,nsize,lvl,lowest,ncycles,presmooth, postsmooth,exact_solu,degrees)
    use gmres_solver

    integer,intent(in) :: nsize,presmooth,postsmooth,lowest,lvl, ncycles
    real, dimension(:) :: b
    real, dimension(:), intent(inout) :: x

    class(element), pointer :: elem
    real :: Tstart, Tend
    integer :: i,restsize,prolsize
    real, dimension(:),allocatable :: y,cor,res,prol,u,v

    character(len=*),intent(in)::exact_solu
    integer,dimension(grid%nelem),intent(in) :: degrees

    ! GMRES
    real :: rezid
    integer :: not_conv
    integer :: restart = 30
    integer :: nloops = 50
    integer :: iter

    integer :: iout = 0

    if (nsize /= size(x)) then
      print*,'ERROR in LINMGLVL. Different vector size as expected'
      STOP
    end if


    allocate( y(1:nsize) )
    y(1:nsize) = 0.

    if (lvl>lowest) then      ! RESTRICTION

      !print*,'####',lvl,lowest

      ! PRESMOOTHING
      if (iout.eq.1) then
        print*,'--- PRESMOOTHING RESIDUUM ---  on level',lvl,'\ '
        write(*,"(A6,I2)") 'LEVEL:',lvl
      endif

      do i=1,presmooth,1
        call bJacobi(x,b,nsize,y)
        x = y

        if( iout.eq.1) then
          call MGbMVprod(y,x,nsize)
          write(*,"(I3,F41.15)") i,norm2(b-y)
        end if

      end do

      ! NEW LEVEL

      allocate(cor(1:nsize))

      call MGbMVprod(y,x,nsize)
      call MGrestVec(b-y, nsize, cor(1:nsize), restsize)

      allocate(res(1:restsize))
      res(1:restsize) = cor(1:restsize)
      cor(1:restsize) = 0.

      if (lvl>lowest+1) then
        do i=1,ncycles,1
	  call LINMGLVL_mod(cor(1:restsize),res(1:restsize),restsize, &
                        lvl-1,lowest,ncycles,presmooth,postsmooth, exact_solu, &
                        (/grid%elem(:)%MGdeg/) )
        enddo
      else
        call LINMGLVL_mod(cor(1:restsize),res(1:restsize),restsize, &
                      lvl-1,lowest,ncycles,presmooth,postsmooth, exact_solu, &
                      (/grid%elem(:)%MGdeg/) )
      endif

      allocate(prol(1:nsize))

      prolsize = nsize
!      call MGProlVec(cor(1:restsize),restsize,prol(1:prolsize),prolsize)
      call MGProlVec_new(cor(1:restsize),restsize,prol(1:prolsize),prolsize, degrees)


      if( iout.eq.1) print*,'Malo by byt rovnake?',prolsize,nsize
      if( iout.eq.1) print*,'--- POSTSMOOTHING RESIDUUM --- on level',lvl,'/'

      ! POSTSMOOTHING
      if( iout.eq.1) print*,x(1:10)
      x(1:nsize) = x(1:nsize) + prol(1:prolsize)
      if( iout.eq.1) print*,x(1:10)

      do i=1,postsmooth,1
        call bJacobi(x,b,nsize,y)
        x = y

        if( iout.eq.1) then
          call MGbMVprod(y,x,nsize)
          write(*,"(I3,F41.15)") i,norm2(b-y)
        end if

      end do

    else ! lvl == lowest
        if( iout.eq.1) print*,'--- "EXACT" RESIDUUM ---       on level',lvl,'V'
        if( iout.eq.1) write(*,"(A6,I2)") 'LEVEL:',lvl

!	call gmres(nsize, x, b, restart*nloops, state%linSolver%tol,  &
!            bMVprod, bMViLUprod, restart,  state%linSolver%tol, iout, iter, rezid, &
!            not_conv)

        if( exact_solu.eq.'GMRES' ) then
          call gmres(nsize, x(1:nsize), b(1:nsize), restart*nloops, state%linSolver%tol, &
                   MGbMVprod, MGbMViLUprod, 30,  state%linSolver%tol, iout, iter, rezid, &
                   not_conv)

        elseif( exact_solu.eq.'bJacobi') then
          do i=1,100,1
            call bJacobi(x,b,nsize,y)
             x = y

            if( iout.eq.1) then
              if ( (i==1).or.(i==100) ) then
                call MGbMVprod(y,x,nsize)
                write(*,"(I5,F38.15)") i,norm2(b-y)
              end if
            end if

          end do

        elseif( exact_solu.eq.'UMFPACK') then
          call SolveViaUMFPACK(x,b)

        elseif(exact_solu.eq.'bJacobi100') then
          do i=1,100,1
            call bJacobi(x,b,nsize,y)
             x = y
            call MGbMVprod(y,x,nsize)
            if ( (i==1).or.(i==100) ) then
              if( iout.eq.1) write(*,"(I5,F38.15)") i,norm2(b-y)
            end if
          end do

       elseif(exact_solu.eq.'bJacobi10') then
          do i=1,10,1
            call bJacobi(x,b,nsize,y)
             x = y
            call MGbMVprod(y,x,nsize)
            if ( (i==1).or.(i==10) ) then
              if( iout.eq.1) write(*,"(I5,F38.15)") i,norm2(b-y)
            end if
          end do

       elseif(exact_solu.eq.'bJacobi1') then
          do i=1,1,1
            call bJacobi(x,b,nsize,y)
             x = y
            call MGbMVprod(y,x,nsize)
            if ( (i==1).or.(i==1) ) then
              if( iout.eq.1) write(*,"(I5,F38.15)") i,norm2(b-y)
            end if
          end do

        elseif( exact_solu.eq.'ILU') then
            call MGbMViLUprod(y,b,nsize)
             x = y
            call MGbMVprod(y,x,nsize)
            if ( (i==1).or.(i==100) ) then
              if( iout.eq.1) write(*,"(I5,F38.15)") i,norm2(b-y)
            end if
        else
          print*,'UNKNOWNsnehurka karlin mff cuni cz solver type in LINMGLVL - STOP'
          stop
        endif

    endif


    deallocate(y)
    if ( allocated(cor) ) then
      deallocate(cor,res)
    end if

  end subroutine LINMGLVL_mod


recursive subroutine LINMGLVL_mod_ILU(x,b,nsize,lvl,lowest,ncycles,presmooth, postsmooth,exact_solu,degrees)
    use gmres_solver

    integer,intent(in) :: nsize,presmooth,postsmooth,lowest,lvl, ncycles
    real, dimension(:) :: b
    real, dimension(:), intent(inout) :: x

    class(element), pointer :: elem
    real :: Tstart, Tend
    integer :: i,restsize,prolsize
    real, dimension(:),allocatable :: y,cor,res,prol,u,v

    character(len=*),intent(in)::exact_solu
    integer,dimension(grid%nelem),intent(in) :: degrees

    ! GMRES
    real :: rezid
    integer :: not_conv
    integer :: restart = 30
    integer :: nloops = 50
    integer :: iter

    integer :: iout = 0

    if (nsize /= size(x)) then
      print*,'ERROR in LINMGLVL. Different vector size as expected'
      STOP
    end if


    allocate( y(1:nsize) )
    y(1:nsize) = 0.

    if (lvl>lowest) then      ! RESTRICTION

      !print*,'####',lvl,lowest

      ! PRESMOOTHING
      if (iout.eq.1) then
        print*,'--- PRESMOOTHING RESIDUUM ---  on level',lvl,'\ '
        write(*,"(A6,I2)") 'LEVEL:',lvl
      endif

      ! * only "1" presmoothing iteration *
      do i=1,1,1
        if( 1.eq.1) then
          call MGbMVprod(y,x,nsize)
          write(*,"(I3,F41.15)") i,norm2(b-y)
        end if

        call MGbMViLUprod(y,b,nsize)
        x = y

        if( 1.eq.1) then
          call MGbMVprod(y,x,nsize)
          write(*,"(I3,F41.15)") i,norm2(b-y)
          print*,norm2(b-y)
        end if


      end do

      ! NEW LEVEL

      allocate(cor(1:nsize))

      call MGbMVprod(y,x,nsize)
      call MGrestVec(b-y, nsize, cor(1:nsize), restsize)

      allocate(res(1:restsize))
      res(1:restsize) = cor(1:restsize)
      cor(1:restsize) = 0.

      if (lvl>lowest+1) then
        do i=1,ncycles,1
	  call LINMGLVL_mod_ILU(cor(1:restsize),res(1:restsize),restsize, &
                        lvl-1,lowest,ncycles,presmooth,postsmooth, exact_solu, &
                        (/grid%elem(:)%MGdeg/) )
        enddo
      else
        call LINMGLVL_mod_ILU(cor(1:restsize),res(1:restsize),restsize, &
                      lvl-1,lowest,ncycles,presmooth,postsmooth, exact_solu, &
                      (/grid%elem(:)%MGdeg/) )
      endif

      allocate(prol(1:nsize))

      prolsize = nsize
!      call MGProlVec(cor(1:restsize),restsize,prol(1:prolsize),prolsize)
      call MGProlVec_new(cor(1:restsize),restsize,prol(1:prolsize),prolsize, degrees)


      if( iout.eq.1) print*,'Malo by byt rovnake?',prolsize,nsize
      if( iout.eq.1) print*,'--- POSTSMOOTHING RESIDUUM --- on level',lvl,'/'

      ! POSTSMOOTHING
      if( iout.eq.1) print*,x(1:10)
      x(1:nsize) = x(1:nsize) + prol(1:prolsize)
      if( iout.eq.1) print*,x(1:10)

      ! * only "1" postsmoothing iteration *
      do i=1,1,1

        call MGbMViLUprod(y,b,nsize)
        x = y

        if( iout.eq.1) then
          call MGbMVprod(y,x,nsize)
          write(*,"(I3,F41.15)") i,norm2(b-y)
        end if

      end do

    else ! lvl == lowest
        if( iout.eq.1) print*,'--- "EXACT" RESIDUUM ---       on level',lvl,'V'
        if( iout.eq.1) write(*,"(A6,I2)") 'LEVEL:',lvl

!	call gmres(nsize, x, b, restart*nloops, state%linSolver%tol,  &
!            bMVprod, bMViLUprod, restart,  state%linSolver%tol, iout, iter, rezid, &
!            not_conv)

        if( exact_solu.eq.'GMRES' ) then
          call gmres(nsize, x(1:nsize), b(1:nsize), restart*nloops, state%linSolver%tol, &
                   MGbMVprod, MGbMViLUprod, 30,  state%linSolver%tol, iout, iter, rezid, &
                   not_conv)
        elseif( exact_solu.eq.'ILU') then
          do i=1,1,1
            call MGbMViLUprod(y,b,nsize)
             x = y

            if( iout.eq.1) then
              if ( (i==1).or.(i==100) ) then
                call MGbMVprod(y,x,nsize)
                write(*,"(I5,F38.15)") i,norm2(b-y)
              end if
            end if

          end do

        elseif( exact_solu.eq.'UMFPACK') then
          call SolveViaUMFPACK(x,b)

       elseif(exact_solu.eq.'bJacobi100') then
          do i=1,100,1
            call bJacobi(x,b,nsize,y)
             x = y
            call MGbMVprod(y,x,nsize)
            if ( (i==1).or.(i==100) ) then
              if( iout.eq.1) write(*,"(I5,F38.15)") i,norm2(b-y)
            end if
          end do

       elseif(exact_solu.eq.'bJacobi10') then
          do i=1,10,1
            call bJacobi(x,b,nsize,y)
             x = y
            call MGbMVprod(y,x,nsize)
            if ( (i==1).or.(i==10) ) then
              if( iout.eq.1) write(*,"(I5,F38.15)") i,norm2(b-y)
            end if
          end do

      elseif(exact_solu.eq.'bJacobi1') then
          do i=1,1,1
            call bJacobi(x,b,nsize,y)
             x = y
            call MGbMVprod(y,x,nsize)
            if ( (i==1).or.(i==1) ) then
              if( iout.eq.1) write(*,"(I5,F38.15)") i,norm2(b-y)
            end if
          end do
        else
          print*,'UNKNOWN solver type in LINMGLVL - STOP'
          stop
        endif

    endif


    deallocate(y)
    if ( allocated(cor) ) then
      deallocate(cor,res)
    end if

  end subroutine LINMGLVL_mod_ILU


  !> copy of elem%w -> elem%MGw
  subroutine MGcopyElem(elem)
    type(element), intent(inout) :: elem

    elem%MGw(1, MGv, 1:elem%dof * ndim) = elem%w(0, 1:elem%dof * ndim)
  end subroutine MGcopyElem

  !> restriction of elem%MGw
  subroutine RestrictElem(elem)
    type(element), intent(inout) :: elem
    integer :: i, k, k1, k2, deg, deg1, dof, dof1

    deg = elem%MGdeg
    deg1 = deg - 1  ! KONTROLA
    dof  = (deg  + 1)*(deg  + 2)/2
    dof1 = (deg1 + 1)*(deg1 + 2)/2

    k = 1
    do i=1,ndim
       k1 = (i-1)*dof*ndim
       k2 = (k1 -1 + dof1)*ndim
       elem%MGw(elem%deg - deg1 + 1, MGv, k:k+dof1 -1 * ndim) &
            = elem%MGw(elem%deg - deg1, MGv, k1:k2)

       k = k + dof1*ndim
    enddo

  end subroutine RestrictElem

  !> Block Jacobi method for linear problem (\f$ (A + \eta\M)x= b\f$
  !>
  !> Consider eq \f$ Ax=b \f$.
  !> \f$ A \f$ is decomposed on diagonal component \f$ D \f$ and remaining part \f$ R \f$, i.e. \f$ A=D+R \f$.
  !> \f$ Ax=b -> (D+R)x=b \ldots Dx=b-Rx \ldots Dx^{k+1}=b-Rx^k \ldots x^{k+1}=inv(D)(b-Rx^k) \f$
subroutine bJacobi(xin, b, nsize, xout)

  use matrix_oper_int
  use lapack_oper
  integer :: nsize
  real, dimension(1:nsize), intent(in):: xin, b
  real, dimension(1:nsize), intent(out):: xout
  class(element), pointer :: elem, elem1
  integer :: i,j,k,is, is1,ndof,ndof1,dof,dof1
  real, dimension(:,:), allocatable :: accuM
  integer :: maxsize

  maxsize = maxval(grid%elem%MGdof)*ndim
  allocate(accuM(1:maxsize,1:maxsize))

  accuM(:,:)=0.
  xout(:)=0.

  do i=1,grid%nelem
    elem => grid%elem(i)
		dof = elem%MGdof
    ndof = ndim * dof
    is = elem%MGncv

    ! off-diagonal blocks
    do j=1,elem%flen
        k = elem%face(neigh,j)

        if(k > 0) then
          elem1 => grid%elem(k)
					dof1 = elem1%MGdof
          ndof1 = ndim * dof1
          is1 = elem1%MGncv

          xout(is : is+ndof-1) = xout(is : is+ndof-1) &
                  + matmul( elem%block(j)%Mb(1:ndof, 1:ndof1), xin(is1 : is1+ndof1-1) )
        endif

    enddo ! j=1,elem%flen

    xout(is : is+ndof-1) = b(is : is+ndof-1) - xout(is : is+ndof-1)

    ! inverse of diagonal blocks
    if (eta /= 0.) then
      accuM(1:ndof,1:ndof) = elem%block(0)%Mb(1:ndof, 1:ndof)
       do k = 0,ndof-1,dof
           accuM(k+1:k+dof,k+1:k+dof) &
		= accuM(k+1:k+dof,k+1:k+dof) + eta*elem%Mass%Mb(1:dof,1:dof)
       enddo

    else
      ! diagonal block
       accuM(1:ndof,1:ndof) = elem%block(0)%Mb(1:ndof, 1:ndof)
    endif

    call MblockInverse(ndof,accuM(1:ndof,1:ndof))

    xout(is:is+ndof-1) = matmul(accuM(1:ndof,1:ndof),xout(is:is+ndof-1))

  enddo ! i=1,grid%nelem

  deallocate(accum)

end subroutine bJacobi


subroutine MGbMVprod(b,x,nsize)
! Block matrix-vector product in pMG cycle.
    use matrix_oper_int
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    class(element), pointer:: elem,elem1 ! one element
!    integer :: i,j,k, ndof, ndof1, is, is1 !, l_accum
    integer :: i,j,k, dof, ndof, ndof2, is,is1, ndof1, is2
    real, dimension(:),allocatable :: accum
    !real, dimension(:),allocatable, save:: accum

    ! allocate accum once to accomodate for the largest dof.
    allocate(accum(maxval(grid%elem%dof) * ndim ) )
    !l_accum= maxval(grid%elem%dof) * ndim
    !if(size(accum) <= l_accum) then
    !   deallocate(accum)
    !   allocate(accum(1:l_accum))
    !endif

    do i=1,grid%nelem
       elem => grid%elem(i)

       dof = elem%MGdof  ! povodne ndof1
       ndof = ndim*dof   != elem%MGdof * ndim
       is = elem%MGncv

       if (eta /= 0.) then
         do k = 0,ndof-1,dof
           accum(k+1:k+dof) = eta * matmul(elem%Mass%Mb(1:dof, 1:dof), &
                x(is+k: is+k+dof-1))
         enddo
         ! diagonal block
         accum(1:ndof) = accum(1:ndof) &
              + matmul(elem%block(0)%Mb(1:ndof, 1:ndof), x(is: is+ndof-1) )
       else
         ! diagonal block
         accum(1:ndof) = matmul(elem%block(0)%Mb(1:ndof, 1:ndof), x(is: is+ndof-1) )
       endif


       !! off-diagonal blocks
       do j=1,elem%flen
          k = grid%elem(i)%face(neigh,j)

          if(k > 0) then
             elem1 => grid%elem(k)
             ndof1 = elem1%MGdof * ndim
             is1 = elem1%MGncv

             accum(1:ndof) = accum(1:ndof) &
                  + matmul(elem%block(j)%Mb(1:ndof, 1:ndof1), x(is1: is1+ndof1-1) )
          endif
       enddo
       b(is: is+ndof-1) = accum(1:ndof)
    enddo

    deallocate(accum)

  end subroutine MGbMVprod



!> ILU preconditioning: \f$ b = (LU)^{-1}x \f$,  \f$ LU  \f$ is the incomplete LU block
  !> preconditioner having the same structure as matrix
  subroutine MGbMViLUprod(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    real, dimension(:), allocatable :: y
    class(element), pointer:: elem, elem1 ! one element
    type(Mblock) :: Loc

    integer :: i, ndof, is, j, i1, ndof1, is1, ii

    call InitMblock(Loc, grid%elem(1)%MGdof * ndim, grid%elem(1)%MGdof * ndim)

    allocate(y(1:nsize) )

    !! L solution
    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof = elem%MGdof * ndim
       is = elem%MGncv

       y(is: is+ndof-1) = x(is: is+ndof-1)

       do j=1,elem%flen
          i1 = elem%face(neigh,j)
          if(i1 > 0 .and. i1 < i) then
             elem1 => grid%elem(i1)
             ndof1 = elem1%MGdof * ndim
             is1 = elem1%MGncv

             y(is: is+ndof-1) = y(is: is+ndof-1) &
                  - matmul(elem%ILU(j)%Mb(1:ndof, 1:ndof1), y(is1: is1+ndof1-1) )

             if(is1 > is) print*,'Problem MGbMViLUprod!: L solution',is,is1, i,i1

          endif
       enddo
    enddo

    !! U solution
    do ii=1,grid%nelem
       i = grid%nelem - ii + 1

       elem => grid%elem(i)
       ndof = elem%MGdof * ndim
       is = elem%MGncv

       do j=1,elem%flen
          i1 = elem%face(neigh,j)

          if( i1 > i) then
             elem1 => grid%elem(i1)
             ndof1 = elem1%MGdof * ndim
             is1 = elem1%MGncv

             y(is: is+ndof-1) = y(is: is+ndof-1) &
                  - matmul(elem%ILU(j)%Mb(1:ndof, 1:ndof1), b(is1: is1+ndof1-1) )

             if(is1 < is) print*,'Problem MGbMViLUprod! U solution',is,is1, i,i1
          endif
       enddo

       if(ndof .ne. size(Loc%Mb,1)) then
          deallocate (Loc%Mb)
          call InitMblock(Loc, ndof, ndof)
       endif

       Loc%Mb(1:ndof,1:ndof) = grid%elem(i)%ILU(0)%Mb(1:ndof,1:ndof)
       call MblockInverse(ndof, Loc%Mb)

       b(is: is+ndof-1) = matmul(Loc%Mb(1:ndof,1:ndof), y(is: is+ndof-1) )

    enddo


    deallocate (Loc%Mb)
    deallocate(y)

  end subroutine MGbMViLUprod


subroutine MGrestVec(x, nsize, xout, outsize)
! pMG-restriction operator.
! Vector x(1:nsize) is restricted on vector xout(1:outsize), xout(outsize+1:nsize) is garbage.

! ..POTREBUJEM globalne deklarovane: d, ndim

  ! INPUT
  real, dimension(1:nsize), intent(in) :: x	! ..MOZE BYT intent(inout), ROZMYSLIET..
  integer, intent(in) :: nsize			! old problem size
!  integer,dimension(grid%nelem),intent(inout) :: degrees
  ! OUTPUT
  real, dimension(1:nsize), intent(out) :: xout
  integer, intent(out) :: outsize		! new problem size, i.e. restricted problem size
  integer :: MGnsize
  ! INNER
  class(element), pointer :: elem
  integer :: i, k, deg, dof, ndof, MGdof, is, MGdeg, MGndof
  integer :: d 				! global dimension 2 or 3
  !integer :: ndim = 4
  integer :: MGminP = 1				! minimal polynomial degree of aproximation used in pMG algorithm

  d = nbDim

  MGnsize = 0

  do i=1,grid%nelem
    elem => grid%elem(i)

    deg = elem%MGdeg		! actual polynomial degree on current element (will be reduced)
    dof = elem%MGdof		! #DOF on element
    ndof = ndim*dof		! size of problem on the current element, = \frac{d+2}{d!} \Pi_{j=1}^{d}(deg+j)
    is = elem%MGncv		! position in state vector for current element and MG-level

    if (deg == MGminP) then	! we cannot decrease pol.degree on current element, because we reached minimum pol.degree
      ! MGdof = dof		! unchanged
      ! MGdeg = deg		! unchanged
      ! MGndof = ndof		! unchanged, = ndim*dof

      ! update stored      ! update stored information for current element information for current element
      ! elem%MGdof = dof
      ! elem%MGdeg = deg
      elem%MGncv = MGnsize+1	!first record (in xout) for current element

      ! copying WHOLE relevant part of vector:
print*,MGnsize,ndof,size(xout),size(x)
print*,grid%nelem,nsize,MGnsize,elem%deg,elem%MGdeg,is
!pause

      !xout(1 : MGnsize+ndof ) = x( is : is+ndof-1 ) !FIXME??
      xout(MGnsize+1 : MGnsize+ndof ) = x( is : is+ndof-1 )

      MGnsize = MGnsize + ndof	! yet last record in xout

    else ! we decrease polynomial degree
      MGndof = ndof*deg/(deg+d)	! #DOF is reduced
      MGdeg = deg-1		! aprox. degree is reduced= \frac{d+2}{d!} \Pi_{j=1}^{d}(deg+j)
      MGdof = MGndof/ndim

      ! update stored information for current element
      elem%MGdof = MGdof
      elem%MGdeg = MGdeg
      elem%MGncv = MGnsize+1	!first record (in xout) for current element

      ! copying RESTRICTED relevant parts of vector:
      do k=1,ndim
        xout( MGnsize+1 : MGnsize+MGdof ) = x( is + (k-1)*dof : is + (k-1)*dof + MGdof -1 )
        MGnsize = MGnsize + MGdof
      enddo

    endif
  enddo ! i=1:grid%nelem

  outsize = MGnsize

end subroutine MGrestVec


subroutine MGprolVec_new(x, nsize, xout, outsize, degrees)
! pMG prolongation operator
! Vector x(1:nsize) is prolongated on vector xout(1:state%nsize), xout(outsize+1:state%nsize) is garbage.

!..POTREBUJEM globalne deklarovane: d, ndim, state%nsize

  ! INPUT
  real, dimension(1:nsize), intent(in) :: x
  integer, intent(in) :: nsize			! old problem size
  ! OUTPUT
! tu je problem, malo by tam byt outsize namiesto 1:state%nsize
  real, dimension(1:outsize), intent(out) :: xout
  integer, intent(inout) :: outsize		! new problem size, i.e. restricted problem size
  integer :: MGnsize
  integer,dimension(outsize),intent(in) :: degrees
  ! INNER
  class(element), pointer :: elem
  integer :: i, k, deg, dof, ndof, is, MGdeg, MGdof, MGndof, maxdeg
  integer :: d				! global dimension 2 or 3
  !integer :: ndim = 4
  !integer :: MGmaxP = 1			! maximal polynomial degree of aproximation used in pMG algorithm

  xout(:) = 0.
  MGnsize = 0
  d = nbDim

  do i=1,grid%nelem,1
    elem => grid%elem(i)

    maxdeg = elem%deg		! ?default polynomial degree
    deg = elem%MGdeg		! actual polynomial degree on current element (will be raised)
    dof = elem%MGdof		! #DOF on element
    ndof = ndim*dof		! size of problem on the current element, = \frac{d+2}{d!} \Pi_{j=1}^{d}(deg+j)
    is = elem%MGncv		! position in state vector for current element and MG-level

!    if (deg == maxdeg) then	! maximal pol. degree is reached
    if (deg == degrees(i)) then	! maximal pol. degree is reached
      ! MGdof = dof		! unchanged
      ! MGdeg = deg		! unchanged
      ! MGndof = ndof 		! unchanged, = ndim*dof

      ! update stored information for current element
      ! elem%MGdof = dof
      ! elem%MGdeg = deg
      elem%MGncv = MGnsize+1	! first record (in xout) for current element

      ! copying WHOLE relevant part of vector:
      xout( MGnsize+1 : MGnsize+ndof ) = x( is : is+ndof-1 )

      MGnsize = MGnsize + ndof	! yet last record in xout

    else
      MGdeg = deg+1
      MGndof = ndof*(MGdeg+d)/MGdeg
      MGdof = MGndof/ndim

      elem%MGdeg = MGdeg
      elem%MGdof = MGdof
      elem%MGncv = MGnsize+1	!first record (in xout) for current element

!print*,'#elem',grid%nelem,'deg=',deg,nsize,outsize

      do k=1,ndim,1
        xout( MGnsize+1 : MGnsize+dof ) = x( is+(k-1)*dof : is+k*dof - 1 )
        MGnsize = MGnsize + MGdof
      enddo

    endif

  enddo !i=1:grid%elem

  outsize = MGnsize

end subroutine MGprolVec_new


subroutine MGprolVec(x, nsize, xout, outsize)
! pMG prolongation operator
! Vector x(1:nsize) is prolongated on vector xout(1:state%nsize), xout(outsize+1:state%nsize) is garbage.

!..POTREBUJEM globalne deklarovane: d, ndim, state%nsize

  ! INPUT
  real, dimension(1:nsize), intent(in) :: x
  integer, intent(in) :: nsize			! old problem size
  ! OUTPUT
! tu je problem, malo by tam byt outsize namiesto 1:state%nsize
  real, dimension(1:state%nsize), intent(out) :: xout
  integer, intent(out) :: outsize		! new problem size, i.e. restricted problem size
  integer :: MGnsize
  ! INNER
  class(element), pointer :: elem
  integer :: i, k, deg, dof, ndof, is, MGdeg, MGdof, MGndof, maxdeg
  integer :: d				! global dimension 2 or 3
  !integer :: ndim = 4
  !integer :: MGmaxP = 1			! maximal polynomial degree of aproximation used in pMG algorithm

  xout(:) = 0.
  MGnsize = 0
  d = nbDim

  do i=1,grid%nelem,1
    elem => grid%elem(i)

    maxdeg = elem%deg		! ?default polynomial degree
    deg = elem%MGdeg		! actual polynomial degree on current element (will be raised)
    dof = elem%MGdof		! #DOF on element
    ndof = ndim*dof		! size of problem on the current element, = \frac{d+2}{d!} \Pi_{j=1}^{d}(deg+j)
    is = elem%MGncv		! position in state vector for current element and MG-level

    if (deg == maxdeg) then	! maximal pol. degree is reached
      ! MGdof = dof		! unchanged
      ! MGdeg = deg		! unchanged
      ! MGndof = ndof 		! unchanged, = ndim*dof

      ! update stored information for current element
      ! elem%MGdof = dof
      ! elem%MGdeg = deg
      elem%MGncv = MGnsize+1	! first record (in xout) for current element

      ! copying WHOLE relevant part of vector:
      xout( MGnsize+1 : MGnsize+ndof ) = x( is : is+ndof-1 )

      MGnsize = MGnsize + ndof	! yet last record in xout

    else
      MGdeg = deg+1
      MGndof = ndof*(MGdeg+d)/MGdeg
      MGdof = MGndof/ndim

      elem%MGdeg = MGdeg
      elem%MGdof = MGdof
      elem%MGncv = MGnsize+1	!first record (in xout) for current element

      do k=1,ndim,1
        xout( MGnsize+1 : MGnsize+dof ) = x( is+(k-1)*dof : is+k*dof - 1 )
        MGnsize = MGnsize + MGdof
      enddo

    endif

  enddo !i=1:grid%elem

  outsize = MGnsize

end subroutine MGprolVec


!> linear MG solution
  subroutine pMGSolver( )
	integer :: i, nsize, outsize, outsize2,outsize3, outsize4, dummysize
	real, dimension(:), allocatable :: MGX,MGB,MGY, &
					   MGB2,MGX2,MGY2, &
				 	   MGB3,MGX3,MGY3, &
				           MGB4,MGX4,MGY4, &
				           MGB5,MGX5,MGY5, &
				           MGB6,MGX6,MGY6, &
					   DUMMY
	type(Newton_type), pointer :: Newton
        class(element), pointer :: elem
	real :: Tstart,Tend;

	! initialize: init solution, RHS
	Newton => state%nlSolver

	print*,' '
	print*,' Newton-Multigrid procedure'

	print*,' '

	! simple initof MG structure part
	grid%elem(:)%MGdof = grid%elem(:)%dof
	grid%elem(:)%MGncv = grid%elem(:)%ncv

	allocate(MGX(1:state%nsize))
	allocate(MGB(1:state%nsize))
	allocate(MGY(1:state%nsize))

	MGX = Newton%x
	MGB = (/(2*i,i=1,state%nsize,1)/)

	call cpu_time(Tstart)
!
! 1. iterace
!
	do i=1,3,1
		call bJacobi(MGX,MGB,state%nsize,MGY)
		MGX=MGY
		call MGbMVprod(MGY,MGX,state%nsize)
		print*,i,'resi fine = ',norm2(MGB-MGY)
	end do

	print*,MGX(1:10)
	call MGbMVprod(MGY,MGX,state%nsize)
	print*,i,'resi fine = ',norm2(MGB-MGY)

!
! 1. restrikce
!
	allocate(MGB2(1:state%nsize))
	MGB2 = MGB-MGY;
	print*,i,'resi fine = ',norm2(MGB2)
	call MGrestVec(MGB2, state%nsize, MGY, outsize)

	print*,'original'
	write(*,'(7F10.5)') MGB2(1:20)

	print*,'restrikce'
	write(*,'(7F10.5)') MGY(1:20)

	deallocate(MGB2)
	allocate(MGX2(1:outsize))
	allocate(MGB2(1:outsize))
	allocate(MGY2(1:outsize))

	MGB2 = MGY(1:outsize)
	MGX2(:) = 0.

	call MGbMVprod(MGY2,MGX2,outsize)
!	print*,'MGY2',MGY2(1:20)
!	print*,'MGX2',MGX2(1:20)
!	print*,size(MGX2),size(MGB2),size(MGY2),outsize
!	print*,norm2(MGX2)
!	print*,norm2(MGY2)
!	print*,'xxxxxx resi cors = ',norm2(MGB2(1:outsize)-MGY2(1:outsize)), norm2(MGB2(1:outsize)), norm2(MGY2(1:outsize))

!
! 2. iterace
!

	print*,'MGB',MGB(1:10)
	print*,'MGB2',MGB2(1:10)

	print*,' *** LEVEL',grid%elem(1)%MGdeg,'***'
	do i=1,100,1
	  call bJacobi(MGX2,MGB2,outsize,MGY2)
	  MGX2=MGY2
	  call MGbMVprod(MGY2,MGX2,outsize)
	  print*,i,'resi cors = ',norm2(MGB2-MGY2)
	end do


	print*,'MGX2'
	write(*,'(3F10.5)') MGX2(1:20)
	call MGbMVprod(MGY2,MGX2,outsize)
	print*,i,'resi fine = ',norm2(MGB2-MGY2)
!
! 1. prolongace
!

allocate(DUMMY(1:state%nsize))

	call MGprolVec(MGX2,outsize,DUMMY,dummysize)
	print*,state%nsize,dummysize

!	print*,'original'
!	write(*,'(7F10.5)') MGX2(1:20)

!	print*,'prolongace'
!	write(*,'(7F10.5)') DUMMY(1:20)

!	print*,size(MGX),outsize,dummysize

	MGX(1:dummysize) = MGX(1:dummysize) + DUMMY(1:dummysize)

	do i=1,100,1
	  !print*,i,'solu norm = ',norm2(MGX)
	  call bJacobi(MGX,MGB,dummysize,MGY)
	  MGX = MGY
	  call MGbMVprod(MGY,MGX,dummysize)
	  !print*,i,'resi cors = ',norm2(MGB-MGY)
	end do
print*,i,'resi cors = ',norm2(MGB-MGY)
	print*,'MGX = '
	write(*,'(6F10.5)') MGX(1:20)

RETURN



		allocate(MGB3(1:outsize))
		MGB3 = MGB2-MGY2;
		print*,i,'resi fine = ',norm2(MGB3)
		call MGrestVec(MGB3, outsize, MGY2, outsize2)

!print*,'original = ',MGB2(1:20)
!print*,'restrikc = ',MGY(1:20)
print*,'original'
write(*,'(7F10.5)') MGB3(1:20)

print*,'restrikce'
write(*,'(7F10.5)') MGY2(1:20)

		deallocate(MGB3)
		allocate(MGX3(1:outsize2))
		allocate(MGB3(1:outsize2))
		allocate(MGY3(1:outsize2))

		MGB3 = MGY2(1:outsize2)
		MGX3(:) = 0.

		call MGbMVprod(MGY3,MGX3,outsize2)
print*,'MGY3',MGY3(1:20)
print*,'MGX3',MGX3(1:20)
print*,size(MGX3),size(MGB3),size(MGY3),outsize2
print*,norm2(MGX3)
print*,norm2(MGY3)
		print*,'xxxxxx resi cors = ',norm2(MGB3(1:outsize2)-MGY3(1:outsize2)), norm2(MGB3(1:outsize2)), norm2(MGY3(1:outsize2))


print*,' *** LEVEL',grid%elem(1)%MGdeg,'***'

do i=1,3,1
!	  	print*,i,'solu norm = ',norm2(MGX)
			call bJacobi(MGX3,MGB3,outsize2,MGY3)
			MGX3=MGY3
		call MGbMVprod(MGY3,MGX3,outsize2)
		print*,i,'resi cors = ',norm2(MGB3-MGY3)
		end do


pause

		allocate(MGB4(1:outsize2))
		MGB4 = MGB3-MGY3;
		print*,i,'resi fine = ',norm2(MGB4)
		call MGrestVec(MGB4, outsize2, MGY3, outsize3)

!print*,'original = ',MGB2(1:20)
!print*,'restrikc = ',MGY(1:20)
print*,'original'
write(*,'(7F10.5)') MGB4(1:20)

print*,'restrikce'
write(*,'(7F10.5)') MGY3(1:20)

		deallocate(MGB4)
		allocate(MGX4(1:outsize3))
		allocate(MGB4(1:outsize3))
		allocate(MGY4(1:outsize3))

		MGB4 = MGY3(1:outsize3)
		MGX4(:) = 0.

		call MGbMVprod(MGY4,MGX4,outsize3)
print*,'MGY4',MGY4(1:20)
print*,'MGX4',MGX4(1:20)
print*,size(MGX4),size(MGB4),size(MGY4),outsize3
print*,norm2(MGX4)
print*,norm2(MGY4)
		print*,'xxxxxx resi cors = ',norm2(MGB4(1:outsize3)-MGY4(1:outsize3)), norm2(MGB4(1:outsize3)), norm2(MGY4(1:outsize3))

print*,' *** LEVEL',grid%elem(1)%MGdeg,'***'

! SOLVE EXACT
do i=1,50,1
!	  	print*,i,'solu norm = ',norm2(MGX)
			call bJacobi(MGX4,MGB4,outsize3,MGY4)
			MGX4=MGY4
		call MGbMVprod(MGY4,MGX4,outsize3)
		print*,i,'resi cors = ',norm2(MGB4-MGY4)
		end do

! PROLONGATION
allocate(DUMMY(1:state%nsize))

call MGprolVec(MGX4,outsize3,DUMMY,dummysize)
print*,outsize2,dummysize

print*,'original'
write(*,'(7F10.5)') MGX4(1:20)

print*,'prolongace'
write(*,'(7F10.5)') DUMMY(1:20)

	print*,size(MGX3),outsize2,dummysize



	print*,'Norma korekce',norm2(DUMMY)
	print*,DUMMY(1:5)
	pause


	print*,'PREresiduum'
	call MGbMVprod(MGY3,MGX3,outsize2)
	print*,i,'resi cors = ',norm2(MGB3-MGY3)


	print*,'POSTresiduum1'
	call MGbMVprod(MGY3,MGX3(1:outsize2) - DUMMY(1:dummysize),outsize2)
	print*,i,'resi cors = ',norm2(MGB3-MGY3)

	print*,'POSTresiduum2'
	call MGbMVprod(MGY3,MGX3(1:outsize2) + DUMMY(1:dummysize),outsize2)
	print*,i,'resi cors = ',norm2(MGB3-MGY3)



	MGX3(1:outsize2) = MGX3(1:outsize2) - DUMMY(1:dummysize)


	print*,'POSTresiduum0'
	call MGbMVprod(MGY3,MGX3,outsize2)
	print*,i,'resi cors = ',norm2(MGB3-MGY3)


	pause

	do i=1,3,1
	  !print*,i,'solu norm = ',norm2(MGX)
	  call bJacobi(MGX3,MGB3,outsize2,MGY3)
	  MGX3=MGY3
	  call MGbMVprod(MGY3,MGX3,outsize2)
	  print*,i,'resi cors = ',norm2(MGB3-MGY3)
	end do


call MGprolVec(MGX3,outsize2,DUMMY,dummysize)
print*,outsize,dummysize

print*,'original'
write(*,'(7F10.5)') MGX3(1:20)

print*,'prolongace'
write(*,'(7F10.5)') DUMMY(1:20)

	print*,size(MGX2),outsize,dummysize

	MGX2(1:outsize) = MGX2(1:outsize) - DUMMY(1:dummysize)

	do i=1,3,1
	  !print*,i,'solu norm = ',norm2(MGX)
	  call bJacobi(MGX2,MGB2,outsize,MGY2)
	  MGX2 = MGY2
	  call MGbMVprod(MGY2,MGX2,outsize)
	  print*,i,'resi cors = ',norm2(MGB2-MGY2)
	end do





call MGprolVec(MGX2,outsize,DUMMY,dummysize)
print*,state%nsize,dummysize

print*,'original'
write(*,'(7F10.5)') MGX2(1:20)

print*,'prolongace'
write(*,'(7F10.5)') DUMMY(1:20)

	print*,size(MGX),outsize,dummysize

	MGX(1:dummysize) = MGX(1:dummysize) + DUMMY(1:dummysize)

	do i=1,3,1
	  !print*,i,'solu norm = ',norm2(MGX)
	  call bJacobi(MGX,MGB,dummysize,MGY)
	  MGX = MGY
	  call MGbMVprod(MGY,MGX,dummysize)
	  print*,i,'resi cors = ',norm2(MGB-MGY)
	end do




		call cpu_time(Tend)

		print*,Tend-Tstart,' s'

  end subroutine pMGSolver


!  subroutine ComputeBlockqLUPrecond()
!    integer:: i,k,j
!
!    do i=2,grid%nelem,1
!      do k=1,i-1,1
!    end do
!  end subroutine ComputeBlockqLUPrecond



end module pMultiGrid
