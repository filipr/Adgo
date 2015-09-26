!> subroutines for MulritGrid (MG) solvers
module dir_solver
  use main_data
  use paramets
  use geometry
  use mtxform
  use matrix_oper
!  use matrix_oper_int
!  use gmres_solver
!  use lin_solvers
  use kind_parameters
  use sparse_utils  ! module for sparse matrix operations

  implicit none

  public:: SolveViaUMFPACK		! grid%elem%MG*** need to be initialized

  public:: umfpack_wrapper		! DEMO
  private:: write_umfpack_info		! used by umfpack_wrapper

  public:: test_umfpack_simple		! simple test of umfpack_wrapper
  public:: test_umfpack_DG		! umfpack wrapper for "real" DG matrix

CONTAINS


  subroutine SolveViaUMFPACK(nsize,x,b)
  integer,intent(in) :: nsize
  real,dimension(nsize),intent(out) :: x	!FIXME: velkost?
  real,dimension(nsize),intent(in) :: b		!FIXME: velkost?
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

!print*,size(Ap),size(Ai),size(Ax),state%nonzero,state%nsize
!print*,Ap(1)
!print*,Ap(state%nsize+1)
!print*,''
!print*,'Ap = ',Ap
!print*,''
!print*,'Ai = ',Ai
!print*,'Ax',Ax(1:3)
!print*,'DG',grid%elem(1)%block(0)%Mb(1:3,1)+eta*grid%elem(1)%Mass%Mb(1:3,1)
!pause

  call umfpack_wrapper(Ax,Ai,Ap,x,b, iout)  

  end subroutine SolveViaUMFPACK


  !> Direct solver via UMFPACK subroutine - DEMO
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

  ! >> FILE NUMBERING << for symbolic analysis and numerical factorization, i.e.
  ! s<filenum>.umf and n<filenum>.umf, respectively
  filenum = 0

  N = size(Ap) -1

  call umf4def(control)	

  ! >> PRINT CONTROL PARAMETERS <<
  ! set control (1) to 1 to print error messages only
  !control(1) = 2
  control(1) = 1
  call umf4pcon(control)

  ! >> PRE-ORDER AND SYMBOLIC ANALYSIS <<
  call umf4sym(N, N, Ap, Ai, Ax, symbolic, control, info)
  ! check umf4sym error condition
  if (info(1) < 0) then
    print *, "Error occurred in umf4sym: ", info(1)
    stop
  else
    if(iout.eq.1) call write_umfpack_info("symbolic analysis:",info)
  endif

  ! NUMERIC FACTORIZATION
  call umf4num(Ap, Ai, Ax, symbolic, numeric, control, info)
  ! check umf4num error condition
  if (info(1) < 0) then
    print *, "Error occurred in umf4num: ", info(1)
    stop
  else
    if(iout.eq.1) call write_umfpack_info("numeric factorization:",info)
  endif

  ! >> SAVE THE SYMBOLIC ANALYSIS << to the file s<filenum>.umf
  ! note that this is not needed until another matrix is
  ! factorized, below.
  call umf4ssym(symbolic, filenum, status)
  if (status < 0) then
    print *, "Error occurred in umf4ssym: ", status
    stop
  endif

  ! >> SAVE THE LU FACTORS << to the file n<filenum>.umf
  call umf4snum(numeric, filenum, status)
  if (status < 0) then
    print *, "Error occurred in umf4snum: ", status
    stop
  endif    

  ! FREE THE SYMBOLIC ANALYSIS
  call umf4fsym(symbolic)
  ! free the numeric factorization
  call umf4fnum(numeric)

  ! >> LOAD THE NUMERIC FACTORIZATION << back in (filename: n<filenum>.umf)
  call umf4lnum(numeric, filenum, status)
  if (status < 0) then
    print *, "Error occurred in umf4lnum: ", status
    stop
  endif

  ! >> SOLVE << Ax=b, without iterative refinement
  sys = 0
  call umf4sol(sys, x, b, numeric, control, info)
  if (info(1) < 0) then
    print *, "Error occurred in umf4sol: ", info(1)
    stop
  else
    if(iout.eq.1) call write_umfpack_info("solving Ax=b:",info)
  endif
   
  ! FREE THE NUMERIC FACTORIZATION
  call umf4fnum(numeric)
  
  ! >>>  
  ! No LU factors (symbolic or numeric) are in memory at this point.
  ! <<<
    
  ! print final statistics
  call umf4pinf(control, info)
  if(iout.eq.1) call write_umfpack_info("final statistic:",info)

  end subroutine umfpack_wrapper


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







  ! UMFPACK test for simple matix
  subroutine test_umfpack_simple
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

  end subroutine test_umfpack_simple



  ! -- ONLY FOR TESTING PURPOUSES
  subroutine test_umfpack_DG(yy,bb)
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

  end subroutine test_umfpack_DG


end module dir_solver
