!> basic LAPACK matrix operations for vectors and matrices
module lapack_oper
  implicit none

  public:: VectorNorm         ! returns a norm of vector
  public:: Distance           ! returns a distance of two nodes = vector
  public:: SDot               ! scalar product
  public:: Area               ! returns of area forming by three nodes
  public:: Area3D             ! returns of volumes forming by four nodes
  public:: MblockInverse      ! computation of the inverse of a block
  public:: SchurComplements   ! solving of saddle point problem by Schur complements
  public:: SchurComplementsNEW ! solving of saddle point problem by Schur complements without inversion
  public:: factorial

  public :: SolveLocalMatrixProblem
  public :: SolveLocalTrasposedMatrixProblem
contains


  function VectorNorm(x)
    real :: VectorNorm
    real, dimension(:), intent(in) :: x

    VectorNorm = sqrt(dot_product(x,x))

  end function VectorNorm

  function Distance(p,q)
    real :: Distance
    real, dimension(:), intent(in) :: p,q
    if(size(p) /= size(q) ) then
       print*,'Different dimension of p and q in Distance(p,q)'
       Distance = 0.
    else
       Distance = sqrt(dot_product(q-p,q-p))
    endif
  end function Distance

  function SDot(p,q)
    real :: SDot
    real, dimension(:), intent(in) :: p,q
    if(size(p) /= size(q) ) then
       print*,'Different dimension of p and q in SDot(p,q)'
       SDot = 0.
    else
       SDot = dot_product(p, q)
    endif
  end function SDot

  function Area(p,q,r)
    ! compute an area of triangle
    real :: Area
    real, dimension(:), intent(in) :: p,q,r
    if(size(p) /= 2 .or. size(q) /=2 .or. size(r) /=2 ) then
       print*,'Bad dimension of p, q or r in Area(p,q,r)'
       Area = 0.
    else
       Area = p(1)*(q(2) - r(2)) + q(1)*(r(2) - p(2)) + r(1)*(p(2)-q(2) )
       Area = Area /2.
    endif

    if(Area <= 0) then
       print*,'Negative area, see file "area"  ', Area
       open(10, file='area', status='UNKNOWN')
       write(10,*) p(:)
       write(10,*) q(:)
       write(10,*) r(:)
       stop
    endif

  end function Area

  function Area3D(p,q,r,s)
    ! compute an area of triangle
    real :: Area3D
    real, dimension(:), intent(in) :: p,q,r,s
    real, dimension(1:3) :: v1, v2, v3
    if(size(p) /= 3 .or. size(q) /=3 .or. size(r) /=3 .or. size(s) /=3 ) then
       print*,'Bad dimension of p, q, r or s in Area(p,q,r,s)'
       Area3D = 0.
    else
       v1(:)=q(:)-p(:)
       v2(:)=r(:)-p(:)
       v3(:)=s(:)-p(:)
       Area3D = v1(1)*v2(2)*v3(3)+v1(2)*v2(3)*v3(1)+v1(3)*v2(1)*v3(2) &
             -v1(1)*v2(3)*v3(2)-v1(2)*v2(1)*v3(3)-v1(3)*v2(2)*v3(1)
       Area3D = Area3D /6.

! if (area3D < 1.0e-13) then
! print*,'#### ',Area3D
! v1(:) = q(:)-p(:)
! print*, '## ',sum(v1(:)*v1(:))**0.5
! v1(:) = r(:)-p(:)
! print*, '## ',sum(v1(:)*v1(:))**0.5
! v1(:) = s(:)-p(:)
! print*, '## ',sum(v1(:)*v1(:))**0.5
! v1(:) = r(:)-q(:)
! print*, '## ',sum(v1(:)*v1(:))**0.5
! v1(:) = s(:)-q(:)
! print*, '## ',sum(v1(:)*v1(:))**0.5
! v1(:) = s(:)-r(:)
! print*, '## ',sum(v1(:)*v1(:))**0.5
! end if

    endif

    if(Area3D <= 0) then
       print*,'Negative area'
       stop
    endif

  end function Area3D


  !> evaluatiotn of the positive and negative part of the matrix
  subroutine  Matrix_plus_minus(n, P, Ppm)
    integer, intent(in) :: n
    real, dimension(1:n,1:n), intent(in) :: P
    real, dimension(1:2,1:n,1:n), intent(out) :: Ppm
    !external:: dgetri, dgetrf          ! subroutines from LAPACK
    real, dimension(:), allocatable :: ident, work
    integer :: info, iwork

  end subroutine Matrix_plus_minus

  !> evaluation of the inverse (\f$ A^{-1} \f$ ) to \f$ A\in R^{n\times n} \f$
  subroutine  MblockInverseOrig(n, A)
    integer, intent(in) :: n
    real, dimension(1:n,1:n), intent(inout) :: A
    external:: dgetri, dgetrf          ! subroutines from LAPACK
    real, dimension(:), allocatable :: ident, work
    integer :: info, iwork

    !do i=1,n
    !   write(*,'(i5,100es12.4)') i, A(i,:)
    !enddo
    !write(*,*)

    iwork = 100 * 30
    allocate(ident(1:n), work(1:iwork) )
    ident(:) = 1.

    !print*,'@@@', state%space%max_dof

    call DGETRF(n, n, A, n, ident, info )
    if(info /= 0 ) print*,'Problem 1 in MblockInverse in matrix.f90 ', info
    if(info /= 0 ) stop

    call DGETRI(n, A, n, ident, work,  iwork, info )
    if(info /= 0 ) print*,'Problem 2 in MblockInverse in matrix.f90 ', info
    if(info /= 0 ) stop

    deallocate(ident, work)

  end subroutine MblockInverseOrig


  !> evaluation of the inverse (\f$ A^{-1} \f$ ) to \f$ A\in R^{n\times n} \f$
  subroutine  MblockInverse(n, A)
    integer, intent(in) :: n
    real, dimension(1:n,1:n), intent(inout) :: A
    real, dimension(:, :), allocatable :: AA
    integer :: i

    allocate(AA(1:n, 1:n) )
    AA(:,:) = 0.
    do i=1,n
       AA(i,i) = 1.
    enddo

    call SolveLocalMatrixProblem(n, A(1:n, 1:n), n, AA(1:n, 1:n) )

    A(1:n, 1:n) = AA(1:n, 1:n)

    deallocate(AA)

  end subroutine MblockInverse


  !> evaluation of the inverse (\f$ A^{-1} \f$ ) to \f$ A\in R^{n\times n} \f$
  !> DOES NOT VEROFIED !!!!!!!!!!!!!!!!!!!!!
  subroutine  MblockInverseSym(n, A)
    integer, intent(in) :: n
    real, dimension(1:n,1:n), intent(inout) :: A
    external:: dgetri, dgetrf          ! subroutines from LAPACK
    real, dimension(:), allocatable :: ident, work
    integer :: info, iwork, i

    !do i=1,n
    !   write(*,'(i5,100es12.4)') i, A(i,:)
    !enddo
    !write(*,*)

    iwork = 100 * 30
    allocate(ident(1:n), work(1:iwork) )
    ident(:) = 1.

    !print*,'@@@', state%space%max_dof

    !!call DGETRF(n, n, A, n, ident, info )

    call DPOTRF( 'L', n, A, n, INFO )

    if(info /= 0 ) print*,'Problem 1 in MblockInverse in matrix.f90 ', info
    if(info /= 0 ) stop

    do i=1, n
       write(14,'(a6, i5,30es14.6)') 'EDW',i,A(i, 1:n)
    enddo
    stop



    call DGETRI(n, A, n, ident, work,  iwork, info )
    if(info /= 0 ) print*,'Problem 2 in MblockInverse in matrix.f90 ', info
    if(info /= 0 ) stop

    deallocate(ident, work)

  end subroutine MblockInverseSym



  !> solve the problem  (\f$ A x = b \f$ ) to \f$ A\in R^{n\times n} \f$,
  !> \f$ b\f$ represents m RHS, solution overwrites b
  subroutine SolveLocalMatrixProblem(n, A, m, b) !!, prec)
    integer, intent(in) :: n, m
    real, dimension(1:n,1:n), intent(in) :: A
    real, dimension(1:n, 1:m), intent(inout) :: b
    !integer, intent(in), optional :: prec

    external:: DGESV          ! subroutines from LAPACK
    real, dimension(:), allocatable :: ipiv
    real, dimension(:,:,:), allocatable :: bb
    real, dimension(:,:), allocatable :: AA
    real :: val
    integer :: k, i, info

    !if( present(prec)) print*,'WEDE$#################@@@@SSEE', prec

    !print*,'SolveLocalMatrixProblem started ', n,  m

    allocate(ipiv(1:n))
    do k=1, n
       ipiv(k) = k
    enddo


    !if( present(prec) ) then
    allocate( AA(1:n, 1:n)  )
    AA(1:n, 1:n) = A(1:n, 1:n)  ! the original matrix is not changed !!
    !  bb(1:n, 1:m, 1 ) = b(1:n, 1:m)
    !endif

    !print*,'###', prec1

    call DGESV(n, m, AA(1:n, 1:n), n, ipiv, b(1:n, 1:m), n, info)
    if(INFO /= 0) print*,'Warning in lap_sub.f90 LAPACK: DGESV,  INFO = ',info

!     ! iterative refinement
!     ! DOES NOT WORK TOO MUCH !!!!
!         if( present(prec) ) then
!        ! residum
       !val= 0.
       !do k=1, m
       !   bb(1:n, k, 2) = bb(1:n, k, 1) - matmul(AA(1:n, 1:n), b(1:n, k))
       !   val = max(val, VectorNorm(bb(1:n, k, 2)) )
       !   write(*,'(a5,2i5,40es12.4)') '#@REZ',n,k, VectorNorm(bb(1:n, k, 2)) , bb(1:min(n,6), k, 2)
       !enddo

!        if( val > 1E-10) then
!           ! an improvement of the accuracy  necessary

!           write(*,'(a5,2i5,40es12.4)') '#@REZ',n,0, val !, VectorNorm(bb(1:n, 1, 2)) , bb(1:min(n,6), 1, 2)
!           do i=1,1  !3

!              ! actual solution
!              bb(1:n, 1:m, 3 ) = b(1:n, 1:m)

!              ! refreshing of the matrix
!              A(1:n, 1:n) = AA(1:n, 1:n)

!              ! new RHS is the original residuum
!              b(1:n, 1:m) = bb(1:n, 1:m, 2)

!              call DGESV(n, m, A(1:n, 1:n), n, ipiv, b(1:n, 1:m), n, info)
!              if(INFO /= 0) print*,'Warning in lap_sub.f90 LAPACK: DGESV, (22)  INFO = ',info

!              b(1:n, 1:m) =  b(1:n, 1:m) + bb(1:n, 1:m, 3)

!              val= 0.
!              do k=1, m
!                 bb(1:n, k, 2) = matmul(AA(1:n, 1:n), b(1:n, k)) -  bb(1:n, k, 1)
!                 val = max(val, VectorNorm(bb(1:n, k, 2)) )
!                ! write(*,'(a5,2i5,40es12.4)') '#@REZ',n,k, VectorNorm(bb(1:n, k, 2)) , bb(1:min(n,6), k, 2)
!              enddo
!              write(*,'(a5,2i5,40es12.4)') '#@REZ',n,i, val
!           enddo
!           stop
!        endif

!        deallocate( AA, bb)
!     endif

! 10  continue
    deallocate(ipiv)
    deallocate(AA)

  end subroutine SolveLocalMatrixProblem

!> solve the problem  (\f$ A x = b \f$ ) to \f$ A\in R^{n\times n} \f$,
  !> \f$ b\f$ represents m RHS, solution overwrites b
  subroutine SolveLocalTrasposedMatrixProblem(n, A, m, b)
    integer, intent(in) :: n, m
    real, dimension(1:n,1:n), intent(in) :: A
    real, dimension(1:n, 1:m), intent(inout) :: b

    external:: DGESV          ! subroutines from LAPACK
    real, dimension(:), allocatable :: ipiv
    real, dimension(:,:,:), allocatable :: bb
    real, dimension(:,:), allocatable :: AA
    real :: val
    integer :: k, i, info

    allocate(ipiv(1:n))
    do k=1, n
       ipiv(k) = k
    enddo

    call DGESV(n, m, transpose(A(1:n, 1:n)), n, ipiv, b(1:n, 1:m), n, info)
    if(INFO /= 0) print*,'Warning in lap_sub.f90 LAPACK: DGESV,  INFO = ',info

    deallocate(ipiv)

  end subroutine SolveLocalTrasposedMatrixProblem



  !> solve the overdetermined linear algebraic problem \f$ Ax \approx b \f$
  !> by the least squares technique, multiplication by the transposed matrix
  !> \f$ n \ge m \f$, number of RHS is k
  subroutine SolveMatrixProblemLeastSquares(n, m, A, k, b, x)
    integer, intent(in) :: n, m, k
    real, dimension(1:n,1:m), intent(inout) :: A
    real, dimension(1:n,1:k), intent(inout) :: b
    real, dimension(1:m,1:k), intent(inout) :: x
    real, dimension(:, :), allocatable :: ATA

    if(n < m) stop 'too little equations in SolveMatrixProblemLeastSquares'

    allocate(ATA(1:m, 1:m) )

    ATA(1:m, 1:m) = matmul( transpose (A(1:n, 1:m) ), A(1:n, 1:m) )
    x(1:m, 1:k)   = matmul( transpose (A(1:n, 1:m) ), b(1:n, 1:k) ) ! RHS will be overwritten

    call SolveLocalMatrixProblem(m, ATA(1:m, 1:m), k, x(1:m, 1:k))


    deallocate(ATA)
  end subroutine SolveMatrixProblemLeastSquares

  !> evaluation of the inverse (\f$ A^{-1} \f$ ) to \f$ A\in R^{n\times n} \f$
  subroutine  MblockInverseOLD(n, A)
    integer, intent(in) :: n
    real, dimension(1:n,1:n), intent(inout) :: A

    real, dimension(:, :), allocatable :: L, U
    real, dimension(:), allocatable :: f, b, av
    integer :: dof, j, j1, k, k1

    dof = n

    ! TODO optimize

    allocate(L(1:dof, 1:dof), U(1:dof, 1:dof))
    allocate(av(1:dof), f(1:dof), b(1:dof))

    ! unit matrix
    L(1:dof, 1:dof) = 0.
    U(1:dof, 1:dof) = 0.
    do j=1, dof
       L(j,j) = 1.
    enddo

    ! let us compute the inverse elem%Mass by LU decomposition
    do k=1,dof
       do j=k,dof
          U(k,j) = A(k,j) - sum(L(k,1:k-1)*U(1:k-1,j) )
       enddo
       if(k /= dof) then
          do k1=k+1,dof
             L(k1,k) = (A(k1,k) - sum(L(k1,1:k-1)*U(1:k-1,k) )) /U(k,k)
          enddo
       endif
    enddo


    ! backward solution: LU a = e_k,  e_k is a canonical basis
    do k=1,dof
       f(1:dof) = 0.
       f(k) = 1.

       b(1) = f(1)/L(1,1)
       do j=2, dof
          b(j) = (f(j) - sum(b(1:j-1)*L(j,1:j-1) ) ) / L(j,j)
       enddo

       av(dof) = b(dof)/U(dof,dof)

       do j1 = 1, dof-1
          j = dof -j1
          av(j) = (b(j) - sum(av(j+1:dof) * U(j,j+1:dof) ) )/U(j,j)
       enddo

       A(1:dof, k) = av(1:dof)
    enddo

    deallocate(L, U, f, b, av)

  end subroutine MblockInverseOLD

  !> solving of a saddle point problem by Schur complements
  !  (  A    B ) ( x ) = ( f )
  !  ( -B^T  0 ) ( y ) = ( g )
  subroutine SchurComplements(n, m, A, B, ip, f, g, x, y  )
    integer, intent(in) :: n, m
    real, dimension(1:n, 1:n), intent(inout) :: A
    real, dimension(1:n, 1:m), intent(in) :: B
    integer, intent(in) :: ip   ! number of RHS
    real, dimension(1:n, 1:ip), intent(in) :: f
    real, dimension(1:m, 1:ip), intent(in) :: g
    real, dimension(1:n, 1:ip), intent(inout) :: x
    real, dimension(1:m, 1:ip), intent(inout) :: y
    external:: DGESV          ! subroutines from LAPACK
    real, dimension(:,:), allocatable :: A1, BTA1B, AA
    real, dimension(:,:), allocatable :: ff, gg
    !integer, dimension(:), allocatable :: ipiv
    integer :: i, j, k,info
    logical :: iprint

    iprint = .false.
    !iprint = .true.

    ! checking the symmetry of A
    do i=1,n
       do j= i+1, n
          if( abs(a (i,j) - a(j, i)) > 0. ) &
               write(*,'(a6,2i5,4es18.10)') 'WW#@', i,j, abs(a (i,j) - a(j, i)), a(i,j), a(j,i)
       enddo
    enddo

    if(iprint) then
       print*,'SCHUR COMPLEMENTS STARTED'
       call WriteArray(A, n, n, 1, n, 1, n)
       print*,'------------------------------------   A'

       call WriteArray(B, m, m, 1, m, 1, m)
       print*,'------------------------------------   B'
    endif

    allocate(A1(1:n, 1:n), BTA1B(1:m, 1:m), AA(1:n, 1:n) )

    A1(1:n, 1:n) = A(1:n, 1:n )
    AA(1:n, 1:n) = A(1:n, 1:n )   ! ONLY for verification
    call MblockInverse(n, A1(1:n, 1:n ) )
    !call MblockInverseSym(n, A1(1:n, 1:n ) )

    if(iprint) then
       call WriteArray(A1, n, n, 1, n, 1, n)
       print*,'------------------------------------   A^{-1}'
    endif

    BTA1B(1:m, 1:m) = matmul( transpose(B(1:n, 1:m)), &
         matmul(A1(1:n, 1:n ), B(1:n, 1:m) ) )


    if(iprint) then
       call WriteArray(BTA1B, m, m, 1, m, 1, m)
       print*,'------------------------------------ BTA1B'
    endif

    allocate(gg(1:m, 1:ip), ff(1:n, 1:ip) )

    do i=1,ip
       gg(1:m, i) = matmul(transpose(B(1:n, 1:m)), &
            matmul(A1(1:n, 1:n ),  f(1:n,i ))) +  g(1:m, i )
    enddo

    if(iprint) then
       write(*,'(a6,50es11.3)') 'RHSA1:', gg(:, 1)
       write(*,'(a6,50es11.3)') 'RHSA2:', gg(:, 2)
       print*
    endif

    !allocate(ipiv(1:max(m, n) ))
    !do i=1, max(m, n)
    !   ipiv(i) = i
    !enddo

    !call DGESV(m, ip, BTA1B(1:m, 1:m), m, ipiv(1:m), gg(1:m, 1:ip),  m, info)
    !if(INFO /= 0) print*,'Warning in lap_sup.f90 (1)LAPACK: DGESV,  INFO = ',info

    call SolveLocalMatrixProblem(m, BTA1B(1:m, 1:m), ip, gg(1:m, 1:ip) )


    if(iprint) then
       write(*,'(a6,50es11.3)') 'solA1:', gg(:, 1)
       write(*,'(a6,50es11.3)') 'solA2:', gg(:, 2)
       print*,'-----------------'
    endif

    ! second system
    do i=1,ip
       ff(1:n, i) = f(1:n, i ) -  matmul( B(1:n, 1:m), gg(1:m, i) )
    enddo

    if(iprint) then
       write(*,'(a6,50es11.3)') 'RHSB1:', ff(:, 1)
       write(*,'(a6,50es11.3)') 'RHSB2:', ff(:, 2)
       print*,'-----------------'
    endif

    !call DGESV(n, ip, A(1:n, 1:n), n, ipiv(1:n), ff(1:n, 1:ip),  n, info)
    !if(INFO /= 0) print*,'Warning in lap_sup.f90 (2)LAPACK: DGESV,  INFO = ',info

    call SolveLocalMatrixProblem(n, A(1:n, 1:n), ip, ff(1:n, 1:ip) )


    if(iprint) then
       write(*,'(a6,50es11.3)') 'solB1:', ff(:, 1)
       write(*,'(a6,50es11.3)') 'solB2:', ff(:, 2)
       print*,'-----------------'
    endif

    ! output of the solution
    x(1:n, 1:ip) = ff(1:n, 1:ip)
    y(1:m, 1:ip) = gg(1:m, 1:ip)

    do i=1,ip
       ff(1:n, i) = matmul(AA(1:n, 1:n ), x(1:n, i) ) + matmul( B(1:n, 1:m), y(1:m, i) ) - f(1:n, i )
       gg(1:m, i) = - matmul(  transpose (B(1:n, 1:m) ), x(1:n, i) ) - g(1:m, i )

       if(VectorNorm ( ff(1:n, i)) > 1E-6 .or. VectorNorm ( gg(1:m, i)) > 1E-6) then
          write(999,*) n+m, n, m
          do k=1, n
             do j=1,n
                write(999,*) AA(k,j), k,j
             enddo
             do j=1,m
                write(999,*) B(k,j), k,j+n
             enddo
          enddo
          do k=1,m
             do j=1, n
                write(999,*) - B(j, k) , k+n, j
             enddo
             do j=1, m
                write(999,*) 0., k+n, j+ n
             enddo
          enddo
          print*,'Matrix written in fort.999'
          stop

       endif


       ! if(VectorNorm ( ff(1:n, i)) > 1E-10) &
       !      write(201,'(i4,a4,es12.4,a2,600es8.1)') i,'ff', VectorNorm ( ff(1:n, i)), '|',ff(1:1, i)
       ! if(VectorNorm ( gg(1:m, i)) > 1E-10) &
       !      write(201,'(i4,a4,es12.4,a2,600es8.1)') i, 'gg', VectorNorm(gg(1:m, i) ), '|', gg(1:1, i)

       ! if(VectorNorm ( ff(1:n, i)) > 1E-08) &
       !      write(301,'(i4,a4,es12.4,a2,600es8.1)') i,'ff', VectorNorm ( ff(1:n, i)), '|',ff(1:1, i)
       ! if(VectorNorm ( gg(1:m, i)) > 1E-08) &
       !      write(301,'(i4,a4,es12.4,a2,600es8.1)') i, 'gg', VectorNorm(gg(1:m, i) ), '|', gg(1:1, i)

       ! if(VectorNorm ( ff(1:n, i)) > 1E-06) &
       !      write(401,'(i4,a4,es12.4,a2,600es8.1)') i,'ff', VectorNorm ( ff(1:n, i)), '|',ff(1:1, i)
       ! if(VectorNorm ( gg(1:m, i)) > 1E-06) &
       !      write(401,'(i4,a4,es12.4,a2,600es8.1)') i, 'gg', VectorNorm(gg(1:m, i) ), '|', gg(1:1, i)
    enddo
    !write(201,'(x)')
    !write(201,*) '############################################################'

    !do i=1,ip
    !   write(99,'(i4,a4,es12.4,a2,600es8.1)') i,'ff', VectorNorm ( ff(1:n, i)), '|',x(1:, i)
    !   write(99,'(i4,a4,es12.4,a2,600es8.1)') i, 'gg', VectorNorm(gg(1:m, i) ), '|', y(1:, i)
    !enddo


    deallocate(A1, BTA1B, gg, ff, AA )

  end subroutine SchurComplements


  !> solving of a saddle point problem by Schur complements
  !  (  A    B ) ( x ) = ( f )
  !  ( -B^T  0 ) ( y ) = ( g )
  subroutine SchurComplementsNEW(ilev, irp, nelem, Fdeg, xc, n, m, A, B, ip, f, g, x, y  )
    integer, intent(in) :: ilev, irp, n, m, nelem, Fdeg
    real, dimension(1:2), intent(in) :: xc
    real, dimension(1:n, 1:n), intent(inout) :: A
    real, dimension(1:n, 1:m), intent(in) :: B
    integer, intent(in) :: ip   ! number of RHS
    external:: DGESV          ! subroutines from LAPACK
    real, dimension(1:n, 1:ip), intent(in) :: f
    real, dimension(1:m, 1:ip), intent(in) :: g
    real, dimension(1:n, 1:ip), intent(inout) :: x
    real, dimension(1:m, 1:ip), intent(inout) :: y
    real, dimension(:,:), allocatable :: BTA1B, AA
    real, dimension(:,:), allocatable :: ff, gg , zz
    integer, dimension(:), allocatable :: ipiv
    real, dimension(:,:), allocatable :: D_AF, xx
    real, dimension(:), allocatable :: D_R, D_C, WORK , FERR, BERR
    integer, dimension(:), allocatable ::  IWORK
    integer :: LWORK
    real ::  RCOND
    character*1 :: EQUED
    integer :: i, j, k,info, irp_test, F1, F2, FFFF
    real:: val
    logical :: iprint

    irp_test = 0 !15

    allocate(ipiv(1:max(m, n) ))
    do i=1, max(m, n)
       ipiv(i) = i
    enddo

    !print*,'SCHUR COMPLEMENTS NEW STARTED'

    iprint = .false.
    !if(irp_test == irp )     iprint = .true.

    ! checking the symmetry of A
    !if(iprint) then
    if(irp == irp_test) then
       print*,'SCHUR COMPLEMENTS STARTED', n, m
       if(ip < 0) then
          call WriteArray(A, n, n, 1, n, 1, n)
          print*,'------------------------------------   A'

          !call WriteArray(B, m, m, 1, m, 1, m)
          !print*,'------------------------------------   B'
       else
          call WriteArray(A, n, n, 1, n, 1, n)
          F1 = Fdeg+1
          F2 = (Fdeg+1)*(Fdeg+3) - 3 *F1
          FFFF = nelem*10 + Fdeg
          print*,'DED',FFFF, nelem, Fdeg,F1, F2
          do i=1,n
             if(FFFF == 61) write(*,1061) i, ':',A(i,:)
             if(FFFF == 62) write(*,1062) i, ':',A(i,:)
             if(FFFF == 63) write(*,1063) i, ':',A(i,:)
             if(FFFF == 64) write(*,1064) i, ':',A(i,:)
             if(FFFF == 65) write(*,1065) i, ':',A(i,:)
             if(FFFF == 66) write(*,1066) i, ':',A(i,:)
             if(FFFF == 67) write(*,1067) i, ':',A(i,:)

             if(mod(i,F1+F2) == F1) write(*,'(x)')
             if(mod(i,F1+F2) == 0) write(*,'(x)')
          enddo
       endif
    endif
    !stop

1061 format(i2,a2,6(2es12.4,' |', 2es12.4,' |'))
1062 format(i2,a2,6(3es12.4,' |', 6es12.4,' |'))
1063 format(i2,a2,6(4es12.4,' |',12es12.4,' |'))
1064 format(i2,a2,6(5es12.4,' |',20es12.4,' |'))
1065 format(i2,a2,6(6es12.4,' |',20es12.4,' |'))
1066 format(i2,a2,6(7es12.4,' |',42es12.4,' |'))
1067 format(i2,a2,6(8es12.4,' |',56es12.4,' |'))



    allocate( BTA1B(1:m, 1:m), AA(1:n, 1:n) )

    AA(1:n, 1:n) = A(1:n, 1:n )   ! storing for computations

    allocate( gg(1:m, 1:ip), ff(1:n, 1:ip) , zz(1:n, 1:m+ip) )


    ! FIRST system setting :
    ! matrix   BTA1B = B^T A^{-1} B
    ! RHS       gg = g + B^T A^{-1} f
    zz(1:n,   1: m   ) =  B(1:n, 1:m )
    zz(1:n, m+1: m+ip) =  f(1:n, 1:ip )

    ! scalling
    !do i=1, n
    !   val = AA(i,i)
    !   print*,'iii',i,val
    !   !AA(i, 1:n) = AA(i, 1:n) / val
    !   !zz(i, 1:m+ip) = zz(i, 1:m+ip) / val
    !enddo


    ! ! simple subroutine
    call DGESV(n, m+ip, AA(1:n, 1:n), n, ipiv(1:n), zz(1:n, 1:m+ip),  n, info)
    if(INFO /= 0) print*,'Warning in lap_sup.f90 (0)LAPACK: DGESV,  INFO = ',info, n,m+ip, irp!, grid%x(irp, :)

    ! ! ! more sophisticated subroutine
    ! LWORK = 4 * n
    ! allocate(D_AF(1:n, 1:n), D_R(1:n), D_C(1:n), xx(1:n, 1:m+ip) )
    ! allocate(FERR(1:m+ip), BERR(1:m+ip), WORK(LWORK), IWORK(n) )

    ! !call DGESVX( 'E', 'N', n, m+ip, AA(1:n, 1:n), n, D_AF(1:n, 1:n), n, ipiv(1:n),  &
    ! !     EQUED, D_R(1:n), D_C(1:n), zz(1:n, 1:m+ip), n, xx(1:n, 1:m+ip), n, RCOND, FERR, BERR, &
    ! !     WORK, IWORK, INFO )

    ! call DSYSVX( 'N', 'U', n, m+ip, AA(1:n, 1:n), n, D_AF(1:n, 1:n), n, ipiv(1:n),  &
    !      !EQUED, D_R(1:n), D_C(1:n),  &
    !      zz(1:n, 1:m+ip), n, xx(1:n, 1:m+ip), n, RCOND, FERR, BERR, &
    !      WORK, LWORK, IWORK, INFO )

    ! zz(1:n, 1:m+ip) = xx(1:n, 1:m+ip)


    ! if(irp == 10) &
    !       write(132,'(es10.2,a3,es10.2)')  1./RCOND, '&', FERR(1)!,EQUED

    ! if(irp == irp_test) then
    !    do i=1,m+ip
    !       write(*,'(a6,2i5,3es14.6,a5)') 'DED',info, i, FERR(i), BERR(i), RCOND !,EQUED
    !    enddo
    ! endif


    ! if(INFO /= 0) print*,'Warning in lap_sup.f90 (0)LAPACK: DGESVx,  INFO = ',info, n,m+ip, irp!, grid%x(irp, :)


    ! !if(irp == irp_test) then
    ! !if(irp == irp_test) then
    ! if(irp == 0) then !15) then
    !    do k=1,1 !m
    !       val = VectorNorm(matmul( A(1:n, 1:n ), zz(1:n,k)) -B(1:n,k))
    !       if(val > 1E-10) write(91,'(a6,3i5,300es12.4)') 'Schu1:',ilev, irp, k,val, xc(1:2)
    !       if(irp == irp_test) write(*,'(a6,3i5,300es12.4)') 'Schu1:',ilev, irp, k,val, xc(1:2)
    !       !write(92,'(100es14.6)')  zz(1:n,k)

    !       write(90,'(a6,2i5,10es14.6)') 'DED',info, Fdeg, FERR(1), BERR(1), RCOND, val, xc(1:2)
    !       write(89,'(2i5,10es12.4)') info, Fdeg, A(1,1), A(3*Fdeg+1, 3*Fdeg+1), &
    !            RCOND, val, FERR(1), BERR(1), xc(1:2)
    !       !write(89,*)'------------------'
    !    enddo
    !    !do k=1,ip
    !    !   val = VectorNorm(matmul( A(1:n, 1:n ), zz(1:n,m+k))-f(1:n,k))
    !    !   if(val > 1E-10) write(91,'(a6,3i5,300es12.4)') 'Schu2:',ilev, irp, k, val, xc(1:2)
    !    !   if(irp == irp_test) write(*,'(a6,3i5,300es12.4)') 'Schu2:',ilev, irp, k, val, xc(1:2)
    !    !enddo

    !    !print*, 'x = ',xc(1:2)
    !    !stop
    ! endif


    ! deallocate(D_AF, D_R, D_C, xx, FERR, BERR, WORK, IWORK)
    ! ! more sophisticated subroutine

    !call SolveLocalMatrixProblem(n, AA(1:n, 1:n), m+ip, zz(1:n, 1:m+ip)) !, .true. )

    ! RHS
    do i=1,ip
       gg(1:m, i) = matmul(transpose(B(1:n, 1:m)), zz(1:n, m+ i)) +  g(1:m, i )
    enddo

    ! MATRIX
    BTA1B(1:m, 1:m) = matmul( transpose(B(1:n, 1:m)), zz(1:n, 1:m) )

    !!call WriteArray(BTA1B, m, m, 1, m, 1, m)


    !endif

    ! simple subroutine
    call DGESV(m, ip, BTA1B(1:m, 1:m), m, ipiv(1:m), gg(1:m, 1:ip),  m, info)
    if(INFO /= 0) print*,'Warning in lap_sup.f90 (1)LAPACK: DGESV,  INFO = ',info, m,ip,irp!, grid%x(irp, :)

    ! ! more sophisticated subroutine
    ! allocate(D_AF(1:m, 1:m), D_R(1:m), D_C(1:m), xx(1:m, 1:ip) )
    ! allocate(FERR(1:ip), BERR(1:ip), WORK(4*m), IWORK(m) )

    ! call DGESVX( 'E', 'N', m, ip, BTA1B(1:m, 1:m), m, D_AF(1:m, 1:m), m, ipiv(1:m),  &
    !      EQUED, D_R(1:m), D_C(1:m), gg(1:m, 1:ip), m, xx(1:m, 1:ip), m, RCOND, FERR, BERR, &
    !      WORK, IWORK, INFO )

    ! !do i=1,m+ip
    ! !   write(*,'(a6,2i5,4es14.6)') 'DED',info, i, FERR(i), BERR(i), RCOND
    ! !enddo

    !  if(INFO /= 0) print*,'Warning in lap_sup.f90 (0)LAPACK: DGESVx,  INFO = ',info, n,m+ip, irp!, grid%x(irp, :)

    !  gg(1:m, 1:ip) = xx(1:m, 1:ip)
    !  deallocate(D_AF, D_R, D_C, xx, FERR, BERR, WORK, IWORK)
    ! ! more sophisticated subroutine


    !call SolveLocalMatrixProblem(m, BTA1B(1:m, 1:m), ip, gg(1:m, 1:ip)) !, .true. )


    if(iprint) then
       write(*,'(a6,500es11.3)') 'solA1:', gg(:, 1)
       write(*,'(a6,500es11.3)') 'solA2:', gg(:, 2)
       print*,'-----------------'
    endif

    ! second system
    do i=1,ip
       ff(1:n, i) = f(1:n, i ) -  matmul( B(1:n, 1:m), gg(1:m, i) )
    enddo

    if(iprint) then
       write(*,'(a6,500es11.3)') 'RHSB1:', ff(:, 1)
       write(*,'(a6,500es11.3)') 'RHSB2:', ff(:, 2)
       print*,'-----------------'
    endif

    AA(1:n, 1:n) = A(1:n, 1:n )

    call DGESV(n, ip, AA(1:n, 1:n), n, ipiv(1:n), ff(1:n, 1:ip),  n, info)
    if(INFO /= 0) print*,'Warning in lap_sup.f90 (2)LAPACK: DGESV,  INFO = ',info, n, ip,irp!, grid%x(irp, :)

    !call SolveLocalMatrixProblem(n, AA(1:n, 1:n), ip, ff(1:n, 1:ip) )


    if(iprint) then
       write(*,'(a6,500es11.3)') 'solB1:', ff(:, 1)
       write(*,'(a6,500es11.3)') 'solB2:', ff(:, 2)
       print*,'-----------------'
    endif

    ! output of the solution
    x(1:n, 1:ip) = ff(1:n, 1:ip)
    y(1:m, 1:ip) = gg(1:m, 1:ip)



    !do i=1,ip

    !    ff(1:n, i) = matmul(A(1:n, 1:n ), x(1:n, i) ) + matmul( B(1:n, 1:m), y(1:m, i) ) - f(1:n, i )
    !    gg(1:m, i) = - matmul(  transpose (B(1:n, 1:m) ), x(1:n, i) ) - g(1:m, i )

       ! FOR MATLAB

       ! if(VectorNorm ( ff(1:n, i)) > 1E-6 .or. VectorNorm ( gg(1:m, i)) > 1E-6) then
       !    write(999,*) n+m, n, m
       !    do k=1, n
       !       do j=1,n
       !          write(999,*) AA(k,j), k,j
       !       enddo
       !       do j=1,m
       !          write(999,*) B(k,j), k,j+n
       !       enddo
       !    enddo
       !    do k=1,m
       !       do j=1, n
       !          write(999,*) - B(j, k) , k+n, j
       !       enddo
       !       do j=1, m
       !          write(999,*) 0., k+n, j+ n
       !       enddo
       !    enddo
       !    print*,'Matrix written in fort.999'
       !    stop

       ! endif


        !write(*,'(i4,a4,es12.4,a2,600es8.1)') i,'ff', VectorNorm ( ff(1:n, i)), '|',ff(1:, i)
        !write(*,'(i4,a4,es12.4,a2,600es8.1)') i, 'gg', VectorNorm(gg(1:m, i) ), '|', gg(1:, i)

       !  if(VectorNorm ( ff(1:n, i)) > 1E-10) &
       !       write(201,'(i4,a4,es12.4,a2,600es8.1)') i,'ff', VectorNorm ( ff(1:n, i)), '|',ff(1:1, i)
       !  if(VectorNorm ( gg(1:m, i)) > 1E-10) &
       !       write(201,'(i4,a4,es12.4,a2,600es8.1)') i, 'gg', VectorNorm(gg(1:m, i) ), '|', gg(1:1, i)

       ! if(VectorNorm ( ff(1:n, i)) > 1E-08) &
       !      write(301,'(i4,a4,es12.4,a2,600es8.1)') i,'ff', VectorNorm ( ff(1:n, i)), '|',ff(1:1, i)
       ! if(VectorNorm ( gg(1:m, i)) > 1E-08) &
       !      write(301,'(i4,a4,es12.4,a2,600es8.1)') i, 'gg', VectorNorm(gg(1:m, i) ), '|', gg(1:1, i)

       ! if(VectorNorm ( ff(1:n, i)) > 1E-06) &
       !      write(401,'(i4,a4,es12.4,a2,600es8.1)') i,'ff', VectorNorm ( ff(1:n, i)), '|',ff(1:1, i)
       ! if(VectorNorm ( gg(1:m, i)) > 1E-06) &
       !      write(401,'(i4,a4,es12.4,a2,600es8.1)') i, 'gg', VectorNorm(gg(1:m, i) ), '|', gg(1:1, i)
     !enddo
    !write(201,'(x)')
    !write(201,*) '############################################################'

    !do i=1,ip
    !   write(99,'(i4,a4,es12.4,a2,600es8.1)') i,'ff', VectorNorm ( ff(1:n, i)), '|',x(1:, i)
    !   write(99,'(i4,a4,es12.4,a2,600es8.1)') i, 'gg', VectorNorm(gg(1:m, i) ), '|', y(1:, i)
    !enddo

    deallocate(ipiv)
    deallocate(AA, BTA1B, gg, ff, zz )

  end subroutine SchurComplementsNEW

  !> write matrix or its frame on  the screene
  subroutine WriteArray( A, n1, n2, nx1, kx1, ny1, ky1)
    integer, intent(in) :: n1, n2,nx1, kx1, ny1, ky1
    real, dimension(1:n1, 1:n2), intent(in) :: A
    integer :: nx, kx, ny, ky
    integer:: i

    nx = nx1
    kx = kx1
    ny = ny1
    ky = ky1

    if(kx > n1) then
       kx = n1
       nx = max(1, n1 - 20)
    endif
    if(ky > n2) then
       ky = n2
       ny = max(1, n2 - 20)
    endif

    write(*,'(a28,i5,a3,i5)')'------- Matrix block -------',n1,'x',n2
    write(*,'(a2,20i11)') 'c:',ny, ny+1, ny+2, ny+3, ny+4, ny+5, ny+6, ny+7, ny+8, ny+9, ny+10, ny+11, ny+12, ny+13, ny+14,  ny+15
!    write(*,'(a2,20i11)') 'c:',13,14,15,16,17,18,19,20
    write(*,'(a10,2i5,a3,2i5)') 'frame: ',nx,kx,'x',ny,ky

    do i=nx, kx
    !do i=13, size(A%Mb,1)
    !do i=1, 3
       !write(*,'(40es9.1)') A(i,:)
       !write(*,'(i2,a2,40es11.3)') i,': ',A(i,13:)
       write(*,'(i2,a2,500es11.3)') i,': ',A(i,ny:ky)
    enddo
    write(*,*)' End of matrix block'

  end subroutine WriteArray

  function factorial(n)
    real :: factorial
    integer, intent(in) :: n
    integer :: i

    factorial = 1.
    do i = 2, n
       factorial = factorial * i
    enddo

  end function factorial




  !> is point x left to the lines x1, x2  !!! orientation has to be kept !!!
  function LeftHandSidePlane(x, x1, x2, eps)
    real:: LeftHandSidePlane  ! =1 inside, =0 outside, closeness to edges
    real, dimension(1:2), intent(in) :: x, x1, x2
    real, intent(in) :: eps
    real :: r, rl, pi


    pi = 2* acos(0.)

    rl = sqrt(dot_product(x2 - x1, x2 - x1))
    r = determinant(x(1:2), x1(1:2), x2(1:2) ) / rl

    if(r > eps) then
       LeftHandSidePlane = 1.

    elseif(r < -eps) then
       LeftHandSidePlane = 0.

    else
       !LeftHandSidePlane = abs(r)
       LeftHandSidePlane = (sin( r/ (2*eps)*pi  ) + 1 ) / 2

       !write(*,'(8es12.4)') x(1:2)
       !write(*,'(8es12.4)') x1(1:2)
       !write(*,'(8es12.4)') x2(1:2)
       !write(*,'(8es12.4)') r, LeftHandSidePlane, eps
       !stop "i9r443jd43o"

    endif
  end function LeftHandSidePlane


  !> x ... investigated node
  !> x1,x2,x3 ... trianle vertices
  function InsideTriangle(x, x1, x2, x3, eps, closed)
    real:: InsideTriangle  ! =1 inside, =0 outside, closeness to edges
    real, dimension(1:2), intent(in) :: x, x1, x2, x3
    real, intent(in) :: eps
    logical, intent(in) :: closed   ! closed polygon?
    real :: r1,r2,r3

    r1 = LeftHandSidePlane(x, x1, x2, eps)
    r2 = LeftHandSidePlane(x, x2, x3, eps)
    r3 = 1E+20
    if(closed) r3 = LeftHandSidePlane(x, x3, x1, eps)

    !InsideTriangle = r1*r2*r3
    InsideTriangle = min (r1, r2, r3 )
  end function InsideTriangle

  !> x ... investigated node
  !> x1,x2,x3 ... trianle vertices
  function InsideQuadrilaterall(x, x1, x2, x3, x4, eps, closed)
    real:: InsideQuadrilaterall  ! =1 inside, =0 outside, closeness to edges
    real, dimension(1:2), intent(in) :: x, x1, x2, x3, x4
    real, intent(in) :: eps
    logical, intent(in) :: closed   ! closed polygon?
    real :: r1,r2,r3, r4

    r1 = LeftHandSidePlane(x, x1, x2, eps)
    r2 = LeftHandSidePlane(x, x2, x3, eps)
    r3 = LeftHandSidePlane(x, x3, x4, eps)
    r4 = 1E+20
    if(closed) r4 = LeftHandSidePlane(x, x4, x1, eps)

    ! InsideQuadrilaterall = r1*r2*r3*r4
    InsideQuadrilaterall = min(r1, r2, r3, r4)
  end function InsideQuadrilaterall


  function determinant(x1, x2, x3 )
    real :: determinant
    real, dimension(1:2), intent(in) :: x1, x2, x3

    determinant = x1(1) * (x2(2) - x3(2) ) + x2(1) * (x3(2) - x1(2) )  + x3(1) * (x1(2) - x2(2) )

  end function determinant

end module lapack_oper

