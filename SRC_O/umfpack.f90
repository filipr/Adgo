!> interface with UMFPACK 

module umfpack_iterface
  use main_data

  use kind_parameters
  use sparse_utils  ! module for sparse matrix operations

  implicit none

  public:: solve_umfpack

contains 
  !> solution of My = r by imfpack
  subroutine solve_umfpack( M, y, z )
    real, dimension (1:Mshape%nonzero), intent(in) :: M
    real, dimension (1:Mshape%nsize), intent(in) :: z
    real, dimension (1:Mshape%nsize), intent(inout) :: y

    integer :: n, nz, i, j, p, istat, ncol, nrow       ! Fortran-variables
    integer :: numeric, symbolic, status, sys, filenum ! UMFPACK-variables
    integer :: iread=5,iwrite=6                        ! File-unit numbers
    integer,dimension(:),allocatable           :: Ap,Ai
    real(kind=double),dimension(:),allocatable :: Ax,x,b,r
    real(kind=double),dimension(20)            :: control
    real(kind=double),dimension(90)            :: info
    real(kind=double)                          :: aij, xj
    character(len=72)                          :: title
    logical                                    :: rhs

    integer :: k,l, irepeat
    real :: val

    print*,'@@',M(1)
    !----------------------------------------------------------------
    ! read the Harwell/Boeing matrix
    !----------------------------------------------------------------
!    print*,'Matrix input' 
    !call HBsize_read(iread,nrow,ncol,nz,rhs)
    !if (nrow /= ncol) then
    !   print *,"Can handle square matrices only! STOP."
    !   stop
    !endif
    nrow = Mshape%nsize
    ncol = Mshape%nsize
    nz = Mshape%nonzero

    n = ncol

    allocate(Ap(n+1),Ai(nz),Ax(nz),x(n),b(n),r(n),stat=istat)
    if (istat /= 0) then
       print *,"Error allocating matrices/vectors. STOP"
       stop
    endif

!    call rua_read(iread,Ap,Ai,Ax,title)

    title = "test matrix"
    Ap(1:ncol+1) = Mshape%irwst(1:ncol+1)
    Ai(1:nz) = Mshape%idx(1:nz)
    
    ! transform from the row ordered matrix to the colum ordered matrix
    do i=1,n
       do j=Ap(i), Ap(i+1) - 1  ! index 
          k = Ai(j)
          
          do l=Ap(k), Ap(k+1) - 1
             if(Ai(l) == i) then
                Ax(j) = M(l)
             endif
          enddo
       enddo
    enddo
    ! reorder in to the "ascending order" 
    do i=1,n
       do irepeat = 1, Ap(i+1) - Ap(i)
          l = 0
          do j=Ap(i), Ap(i+1) - 2  ! index 
             if(Ai(j) > Ai(j+1) ) then
                k = Ai(j)
                Ai(j) = Ai(j+1)
                Ai(j+1) = k
                
                val = Ax(j)
                Ax(j) = Ax(j+1)
                Ax(j+1) = val
                
                l = l + 1
             endif
          enddo
          if(l == 0) goto 15
       enddo
15     continue
    enddo

    !Ax(1:nz) = M(1:nz)

    !print*,'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
    !open(10, file='output.file', status='UNKNOWN')
    !call rua_write(10,Ap,Ai,Ax,title,'AA')
    !close(10)
    !print*,'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'


    !write(32,'(a15)') 'Mshape:'
    !write(32,'(200i5)') Ap(:)
    
    !j = 6
    !do i=1, Mshape%nonzero/j
    !   write(32,'(2i5,a1,6i5)') j*(i-1)+1,i*j,'|', Ai(j*(i-1)+1:i*j)
    !enddo
    !write(32,'(2i5,a1,6i5)') j*(i-1)+1,i*j,'|', Ai(j*(i-1)+1:)

    !write(32,'(a5)') 'Ax:'
    !j = 3
    !do i=1, nz/j
    !   write(32,'(2i5,a1,3es22.14)') j*(i-1)+1,i*j,'|', Ax(j*(i-1)+1:i*j)
    !enddo
    !   write(32,'(2i5,a1,3es22.14)') j*(i-1)+1,nz,'|', Ax(j*(i-1)+1:nz)

    !write(32,'(a5)') 'rhs:'
    !j = 3
    !do i=1, ncol/j
    !   write(32,'(2i5,a1,3es22.14)') j*(i-1)+1,i*j,'|', z(j*(i-1)+1:i*j)
    !enddo
    !   write(32,'(2i5,a1,3es22.14)') j*(i-1)+1,n,'|', z(j*(i-1)+1:n)

    !write(32,'(a5)') 'x:'
    !j = 3
    !do i=1, ncol/j
    !   write(32,'(2i5,a1,3es22.14)') j*(i-1)+1,i*j,'|', y(j*(i-1)+1:i*j)
    !enddo
    !write(32,'(2i5,a1,3es22.14)') j*(i-1)+1,n,'|', y(j*(i-1)+1:n)

    ! write the matrix shape
    !write(unit=iwrite,fmt=*) "input matrix:"
    !call csc1_spy(iwrite,Ap,Ai)          ! 1-based

    !----------------------------------------------------------------
    ! create the right-hand-side, assume x(i) = 1 + i/n
    !----------------------------------------------------------------

    !do i = 1,n
    !   x(i) = 1.0 + real(i)/real(n)
    !enddo

    ! b = A*x
    !call csc1_vector_mult(Ap,Ai,Ax,y,b)

    !write(32,'(a15)') 'new Ax:'
    !j = 3
    !do i=1, Mshape%nsize/j
    !   write(32,'(2i5,a1,3es22.14)') j*(i-1)+1,i*j,'|', b(j*(i-1)+1:i*j)
    !enddo
    !write(32,'(2i5,a1,3es22.14)') j*(i-1)+1,n,'|', b(j*(i-1)+1:n)


    !write(32,'(a15)') 'Ay - z:'
    !print*,'------------------------'
    !j = 3
    !do i=1, ncol/j
    !   write(32,'(2i5,a1,3es22.14)') j*(i-1)+1,i*j,'|', b(j*(i-1)+1:i*j) - z(j*(i-1)+1:i*j)
    !enddo
    !write(32,'(2i5,a1,3es22.14)') j*(i-1)+1,ncol,'|', b(j*(i-1)+1:ncol) - z(j*(i-1)+1:ncol)


    ! RHS
    b(1:ncol) = z(1:ncol)
    
    x(1:ncol) = y(1:ncol)

    !----------------------------------------------------------------
    ! convert from 1-based to 0-based
    !----------------------------------------------------------------
    Ap = Ap - 1
    Ai = Ai - 1
    

    !call resid(n, Ap, Ai, Ax, y, b, r)
    

    !----------------------------------------------------------------
    ! factor the matrix and save to a file
    !----------------------------------------------------------------
    ! set default parameters


    call umf4def(control)

    ! print control parameters.  set control (1) to 1 to print
    ! error messages only
    !control(1) = 2
    control(1) = 1
    call umf4pcon(control)

    ! pre-order and symbolic analysis
    call umf4sym(n, n, Ap, Ai, Ax, symbolic, control, info)
    
    ! print statistics computed so far
    ! call umf4pinf(control, info) could also be done.
    !write(unit=iwrite,fmt="(a18,            /,  &
    !     &a12, f5.0,      /,  &
    !     &a12, es10.2,a6, /,  &
    !     &a42,            /,  &
    !     &a18, f10.2,a5,  /,  &
    !     &a18, f10.2,a5,  /,  &
    !     &a18, es10.2,    /,  &
    !     &a18, f10.0,     /,  &
    !     &a18, f10.0)")       & ! END OF FORMAT
    !     "symbolic analysis:",                  & ! Data to print from here
    !     "   status:  ",info(1) ,               &
    !     "   time:    ",info(16)," (sec)",      &
    !     "   estimates (upper bound) for numeric LU:",           &
    !     "   size of LU:    ",(info(21)*info(4))/2**20," (MB)",  &
    !     "   memory needed: ",(info(22)*info(4))/2**20," (MB)",  &
    !     "   flop count:    ", info(23),        &
    !     "   nnz (L):       ", info(24),        &
    !     "   nnz (U):       ", info(25)
    ! 
    ! check umf4sym error condition
    if (info(1) < 0) then
       print *, "Error occurred in umf4sym: ", info(1)
       stop
    endif
    
    ! numeric factorization
    call umf4num(Ap, Ai, Ax, symbolic, numeric, control, info)
    
    ! print statistics for the numeric factorization
    ! call umf4pinf(control, info) could also be done.
    !write(unit=iwrite,fmt="(a22,            /,  &
    !     &a12, f5.0,      /,  &
    !     &a12, es10.2,    /,  &
    !     &a32,            /,  &
    !     &a18, f10.2, a5, /,  &
    !     &a18, f10.2, a5, /,  &
    !     &a18, es10.2,    /,  &
    !     &a18, f10.0,     /,  &
    !     &a18, f10.0)")       & ! END OF FORMAT
    !     "numeric factorization:",              & ! Data to print from here
    !     "   status:  ",info(1),                &
    !     "   time:    ",info(66),               &
    !     "   actual numeric LU statistics:",    &
    !     "   size of LU:    ",(info(41)*info(4))/2**20," (MB)",  &
    !     "   memory needed: ",(info(42)*info(4))/2**20," (MB)",  &
    !     "   flop count:    ",info(43),         &
    !     "   nnz (L):       ",info(44),         &
    !     "   nnz (U):       ",info(45)
    
    ! check umf4num error condition
    if (info(1) < 0) then
       print *, "Error occurred in umf4num: ", info(1)
       stop
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
    
    ! No LU factors (symbolic or numeric) are in memory at this point.
    
    !----------------------------------------------------------------
    ! load the LU factors back in, and solve the system
    !----------------------------------------------------------------
    
    ! At this point the program could terminate and load the LU
    ! factors (numeric) from the n0.umf file, and solve the
    ! system (see below).  Note that the symbolic object is not
    ! required.
    
    ! load the numeric factorization back in (filename: n0.umf)
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
    endif
    
    ! free the numeric factorization
    call umf4fnum(numeric)
    
    ! No LU factors (symbolic or numeric) are in memory at this point.
    
    ! print final statistics
!    print *,"Final Statistics"
    call umf4pinf(control, info)
    
    ! print the residual.  x(i) should be 1 + i/n
!    call resid(n, Ap, Ai, Ax, x, b, r)
    
    !----------------------------------------------------------------
    ! load the symbolic analysis back in, and factorize a new matrix
    !----------------------------------------------------------------
    ! Again, the program could terminate here, recreate the matrix,
    ! and refactorize.  Note that umf4sym is not called.
    
    ! load the symbolic factorization back in (filename: s0.umf)
!    call umf4lsym(symbolic, filenum, status)
!    if (status < 0) then
!       print *, "Error occurred in umf4lsym: ", status
!       stop
!    endif
    
!    ! arbitrarily change the values of the matrix but not the pattern
!    do p = 1, nz
!       Ax(p) = Ax(p) + 3.14159 / 100.0
!    enddo
    
!    ! numeric factorization of the modified matrix
!    call umf4num(Ap, Ai, Ax, symbolic, numeric, control, info)
!    if (info(1) < 0) then
!       print *, "Error occurred in umf4num: ", info(1)
!       stop
!    endif
    
!    ! free the symbolic analysis
!    call umf4fsym(symbolic)
    
!    ! create a new right-hand-side, assume x(i) = 7 - i/n
!    do i = 1,n
!       b(i) = 0
!    enddo
!    ! b = A*x, with the modified matrix A (note that A is now 0-based)
!    do j = 1,n
!       xj = j
!       xj = 7 - xj / n
!       do p = Ap(j) + 1, Ap(j+1)
!          i    = Ai(p) + 1
!          aij  = Ax(p)
!          b(i) = b(i)  + aij * xj
!       enddo
!    enddo
    
!    !----------------------------------------------------------------
!    ! solve Ax=b, with iterative refinement
!    !----------------------------------------------------------------
!    print *," "
!    print *,"START   W I T H  Iterative refinement:"
    
!    sys = 0
!    call umf4solr(sys, Ap, Ai, Ax, x, b, numeric, control, info)
!    if (info(1) < 0) then
!       print *, "Error occurred in umf4solr: ", info(1)
!       stop
!    endif
    
!    ! print the residual.  x(i) should be 7 - i/n
!    call resid(n, Ap, Ai, Ax, x, b, r)
!    print *," "
!    print *,"END     With Iterative refinement:"
    
!    !----------------------------------------------------------------
!    ! solve Ax=b, without iterative refinement, broken into steps
!    !----------------------------------------------------------------
    
!    print *," "
!    print *,"START   WITHOUT  Iterative refinement:"
!    ! the factorization is PAQ=LU, PRAQ=LU, or P(R\A)Q=LU.
    
!    ! x = R*b (or x=R\b, or x=b, as appropriate)
!    call umf4scal(x, b, numeric, status)
!    if (status < 0) then
!       print *, "Error occurred in umf4scal: ", status
!       stop
!    endif
    
!    ! solve P'Lr=x for r (using r as workspace)
!    sys = 3
!    call umf4sol(sys, r, x, numeric, control, info)
!    if (info(1) < 0) then
!       print *, "Error occurred in umf4sol: ", info(1)
!       stop
!    endif
    
!    ! solve UQ'x=r for x
!    sys = 9
!    call umf4sol(sys, x, r, numeric, control, info)
!    if (info(1) < 0) then
!       print *, "Error occurred in umf4sol: ", info(1)
!       stop
!    endif
    
!    ! free the numeric factorization
!    call umf4fnum(numeric)
    
!    ! print the residual.  x(i) should be 7 - i/n
!    call resid(n, Ap, Ai, Ax, x, b, r)
!    print *," "
!    print *,"END     WITHOUT  Iterative refinement:"
    

!    write(*,'(a5)') 'x:'
!    j = 3
!    do i=1, n/j
!       write(*,'(2i5,a1,3es22.14)') j*(i-1)+1,i*j,'|', x(j*(i-1)+1:i*j) 
!    enddo
!       write(*,'(2i5,a1,3es22.14)') j*(i-1)+1,n,'|', x(j*(i-1)+1:n) 
!
!    write(*,'(a5)') 'y:'
!    j = 3
!    do i=1, n/j
!       write(*,'(2i5,a1,3es22.14)') j*(i-1)+1,i*j,'|', y(j*(i-1)+1:i*j)
!    enddo
!       write(*,'(2i5,a1,3es22.14)') j*(i-1)+1,n,'|', y(j*(i-1)+1:n)


!    write(*,'(a5)') 'x-y:'
!    j = 3
!    do i=1, n/j
!       write(*,'(2i5,a1,3es22.14)') j*(i-1)+1,i*j,'|', x(j*(i-1)+1:i*j) - y(j*(i-1)+1:i*j)
!    enddo
!       write(*,'(2i5,a1,3es22.14)') j*(i-1)+1,n,'|', x(j*(i-1)+1:n) - y(j*(i-1)+1:n)


    y(1:n) = x(1:n)

    deallocate(Ap,Ai,Ax,x,b,r,stat=istat)

!    print*,'umpfack _ stopped'
!    stop

  end subroutine solve_umfpack


  !=======================================================================
  !== resid ==============================================================
  !=======================================================================
  !> Compute the residual, r = Ax-b, its max-norm, and print the max-norm
  !> Note that A is zero-based.
  
  subroutine resid1(n, Ap, Ai, Ax, x, b, r)
    integer,intent(in)                          :: n
    integer                                     :: j, i, p
    integer,dimension(:),intent(in)             :: Ap,Ai
    real(kind=double)                           :: rmax, rmin, aij
    real(kind=double),dimension(:),intent(in)   :: Ax,x,b
    real(kind=double),dimension(:),intent(out)  :: r
    
    r = 0.0    ! initialize r(i)=0.0 for i=1,n (array expression).
    do j = 1,n
       do p = Ap(j) + 1, Ap(j+1)

          i    = Ai(p) + 1
          aij  = Ax(p)
          r(i) = r(i)  + aij * x(j)
       enddo
    enddo
    
    rmax = 0.0
    rmin = 0.0
    do i = 1, n
       rmax = max(rmax, (  r(i)-b(i)  ))  ! This is not exactly a norm
       rmin = min(rmin, (  r(i)-b(i)  ))  ! but max positive/negative error
    enddo
    
    print *,"Max. pos. error (A*x-b): ", rmax
    print *,"Max. neg. error (A*x-b): ", rmin
    
    return
  end subroutine resid1

end module umfpack_iterface
