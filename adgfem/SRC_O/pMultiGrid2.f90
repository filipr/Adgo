!> subroutines for MulritGrid (MG) solvers
MODULE pMultiGrid2

  use main_data
  use data_mod
  use paramets
  use geometry
  use matrix_oper
  use matrix_oper_int
  use lapack_oper
  use dir_solver
  use f_mapping
  use ama_L2interpol

  implicit none

  type,public :: MGSetting
     integer :: Vcycle = 1
     integer :: presmooth = 1
     integer :: postsmooth = 1
  end type MGSetting

  type(MGSetting) :: MGRun

  public:: L2Res

  public:: pMGLinSolverRecu4LINSOL
  public:: pMGLinSolverRecu
  public:: LINMGLVL_mod
  public:: RestrictElem
  public:: MGcopyElem
  public:: bJacobi
  public:: bGS
  public:: MGbMViLUprod
  public:: MGrestVec
  public:: MGprolVec_new
  public:: MGprolVec
  public:: pMGSolver

  public::  Compute_MG_BlockDiagPrecond

CONTAINS





  FUNCTION L2Res(x,b,nsize,prod)

    real :: L2Res
    ! input
    integer,intent(in) :: nsize
    real,dimension(1:nsize),intent(in) :: x	! solution
    real,dimension(1:nsize),intent(in) :: b ! RHS
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

  END FUNCTION L2Res



  SUBROUTINE pMGLinSolverRecu4LINSOL(nsize, eta, b, x, max_deg,exact_solu,rezid,SMOOTHER)

    integer :: nsize                               ! size of the algebraic problem
    real, intent(in):: eta
    real, dimension(1:nsize), intent(inout):: x    ! solution
    real, dimension(1:nsize), intent(inout):: b    ! RHS
    real, intent(inout) :: rezid                   ! reziduum
    character(len=*),intent(in):: exact_solu       ! GMRES, bJacobi

    !real,dimension(1:nsize):: rr
    real,dimension(:), allocatable :: rr
    !character(len=1) :: precond  ! type of preconditioner: ' ', 'D', 'L', 'J'
    !integer:: iout, iter
    !real :: size(-1:3)
    !integer:: restart = 30 !30    ! 50  ! GMRES restarted after 'restart' iterations  !45
    !integer:: nloops = 5  !5     ! 10   ! maximal number of restarted cycles !100	  !40
    !real :: t0, t1, t2
    integer  :: max_deg
    integer,parameter :: min_deg = 1
    integer:: presmooth
    integer:: postsmooth
    integer:: Vcycle

    interface
       subroutine SMOOTHER(x,b,n)
         real,dimension(n),intent(inout):: x
         real,dimension(n),intent(in):: b
         integer,intent(in):: n
       end subroutine SMOOTHER
    end interface


    allocate( rr(1:nsize) )

    presmooth = MGRun%presmooth
    postsmooth = MGRun%postsmooth
    Vcycle = MGRun%Vcycle

    !write(*,"(//a80)") '/----------------------------------------------------------------------------\ '
    write(*,"(1X,A23,6X,2I3,L3)") '----- START of MG-cycle', max_deg, min_deg, state%space%adaptation
    !pause

    !
    call LINMGLVL_mod(x, b, nsize, max_deg, min_deg, Vcycle, presmooth, postsmooth, &
         exact_solu,(/grid%elem(:)%deg/),SMOOTHER)
    !
    !write(*,"(1X,A21,8X,2I3,L3)") '----- END of MG-cycle',max_deg, min_deg, state%space%adaptation
    !
    call bMVprod(rr,x,nsize)    ! rr = Ax
    rr(1:nsize) = rr(1:nsize) - b(1:nsize)                ! rr = Ax - b
    !
    rezid = (dot_product( rr(1:nsize), rr(1:nsize) ) )**0.5
    !
    !write(*,"(/1X,A29,19X,F30.15/49XF30.15)") '----- FINAL MG-cycle reziduum',rezid,L2Res(x,b,nsize,bMVprod)
    !write(*,"(a80//)") '\------------------------------------------------------------------------------/'

    !pause

    !STOP

    deallocate(rr)

  END SUBROUTINE pMGLinSolverRecu4LINSOL



  SUBROUTINE pMGLinSolverRecu(nsize, eta, b, x, rezid, tot_iter, &
       not_conv,exact_solu)

    use matrix_oper_int, mx_eta => eta
    !use lin_solvers, only: ComputeBlockDiagPrecond

    integer :: nsize                               ! size of the algebraic problem
    real, intent(in):: eta
    real, dimension(1:nsize), intent(inout):: x    ! solution
    real, dimension(1:nsize), intent(inout):: b    ! RHS
    real, intent(inout) :: rezid                   ! reziduum
    integer, intent(inout) :: tot_iter             ! number of iterations
!    logical, intent(in) :: precond_update          ! = .false. preconditioner is not update
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
    integer,parameter :: Vcycle = 1

    max_deg = maxval( grid%elem(:)%deg)

    grid%elem(:)%MGdof = grid%elem(:)%dof
    grid%elem(:)%MGncv = grid%elem(:)%ncv

    !print*,'### ZACIATOK LINMGLVL', max_deg,min_deg,state%space%adaptation
    !pause('ZACIATOK LINMGLVL')

    !call LINMGLVL(x, b, nsize, max_deg, min_deg, Vcycle, presmooth, postsmooth,exact_solu)
    call LINMGLVL_mod(x, b, nsize, max_deg, min_deg, Vcycle, presmooth, postsmooth,exact_solu,(/grid%elem(:)%deg/),bJacobi)

    !print*,'### KONIEC LINMGLVL', max_deg,min_deg,state%space%adaptation
    !pause('KONIEC LINMGLVL')

    call bMVprod(rr,x,nsize)    ! rr = Ax
    rr(1:nsize) = rr(1:nsize) - b(1:nsize)                ! rr = Ax - b

    rezid = (dot_product( rr(1:nsize), rr(1:nsize) ) )**0.5

  END SUBROUTINE pMGLinSolverRecu



  RECURSIVE SUBROUTINE LINMGLVL_mod(x,b,nsize,lvl,lowest,ncycles,presmooth, postsmooth,exact_solu, &
       degrees,SMOOTHER)
    use gmres_solver

    integer,intent(in) :: nsize,presmooth,postsmooth,lowest,lvl, ncycles
    real, dimension(:) :: b
    real, dimension(:), intent(inout) :: x

    class(element), pointer :: elem
    real :: Tstart, Tend
    integer :: i,restsize,prolsize
    real, dimension(:),allocatable :: y,cor,res,prol,u,v

    character(len=*),intent(in)::exact_solu
    integer,dimension(1:grid%nelem),intent(in) :: degrees

    ! GMRES
    real :: rezid
    integer :: not_conv
    integer :: restart = 30
    integer :: nloops = 50
    integer :: iter, k, dof

    !integer :: iout = 0 !1
    integer :: iout = 1 !1


    interface
       subroutine SMOOTHER(x,b,nsize)
         real,dimension(nsize),intent(inout):: x
         real,dimension(nsize),intent(in):: b
         integer,intent(in):: nsize
       end subroutine SMOOTHER
    end interface


    if( nsize /= size(x)) then
       print*,'ERROR in LINMGLVL. Different vector size as expected'
       STOP
    end if

    allocate( y(1:nsize) )
    y(1:nsize) = 0.

    if( lvl>lowest) then      ! RESTRICTION


       do  i=1,presmooth,1
          !call bJacobi(x,b,nsize
          !call bGS(x,b,nsize)

          call SMOOTHER(x,b,nsize)

          if( iout.eq.1) then
             call MGbMVprod(grid,y,x,nsize)
             write(*,"(1X,A36,7X,I1,1X,A3,I3,es16.6)") &
                  '----- PRESMOOTHING RESIDUUM on level',lvl,'\  ',i,norm2(b-y)
          end if
       end do

       ! NEW LEVEL
       allocate(cor(1:nsize))

       call MGbMVprod(grid,y,x,nsize)
       call MGrestVec(b-y, nsize, cor(1:nsize), restsize)

       !write(66,'(a18, i6, 500es12.4)') 'b-x',nsize, b-y
       !write(66,'(a18, i6, 500es12.4)') 'cor',restsize, cor(:)
       !write( 66,*) '.,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,..'

       !write( 66,*) '#######################################'
       !dof = (degrees(1) + 1) *(degrees(1) + 2) / 2
       !write(66,'(a8, 8i10)') 'VRE',lvl, nsize, restsize, dof
       !write(66,'(a18, 2a18)') '  ', 'REZID   ', 'RESTRICT   '
       !do k=1,nsize
       !   write(66,'(a8, 2i5,30es18.10)') 'VMG',lvl, k,b(k)-y(k), cor(k)
       !   if(mod(k, dof*ndim) == 0) write(66,*)'-------- elem%i', k / dof /ndim
       !   if(mod(k, dof) == 0) write(66,*)
       !enddo


       allocate(res(1:restsize))
       res(1:restsize) = cor(1:restsize)
       cor(1:restsize) = 0.

       !write(*,'(a8, 30es12.4)') 'VMG',res(1:15)
       !write(*,'(a8, 30es12.4)') 'VMG',cor(1:15)
       !print*,'####################################################'


       if( lvl>lowest+1) then
          do  i=1,ncycles,1
             call LINMGLVL_mod(cor(1:restsize),res(1:restsize),restsize, &
                  lvl-1,lowest,ncycles,presmooth,postsmooth, exact_solu, &
                  (/grid%elem(:)%MGdeg/),SMOOTHER )
          end do
       else
          call LINMGLVL_mod(cor(1:restsize),res(1:restsize),restsize, &
               lvl-1,lowest,ncycles,presmooth,postsmooth, exact_solu, &
               (/grid%elem(:)%MGdeg/),SMOOTHER )
       end if

       allocate(prol(1:nsize))

       prolsize = nsize
       call MGprolVec_new(cor(1:restsize),restsize,prol(1:prolsize),prolsize, degrees)

       !write(66,'(a18, i6, 500es12.4)') 'cor',restsize, cor(:)
       !write(66,'(a18, i6, 500es12.4)') 'cor',prolsize, prol(:)
       !write( 66,*) '............................................................'

       !write( 66,*) '#######################################'
       !write(66,'(a8, 6i10)') 'VPR',lvl,  restsize, prolsize, nsize
       !write(66,'(a18, 2a18)') '  ', 'REZID   ', 'PROLONG   '
       !do k=1,nsize
       !   write(66,'(a8, 2i5,30es18.10)') 'VMG',lvl, k, cor(k), prol(k)
       !   if(mod(k, dof*ndim) == 0) write(66,*)'-------- elem%i', k / dof /ndim
       !   if(mod(k, dof) == 0) write(66,*)
       !enddo

       !if( iout.eq.1) print*,'Malo by byt rovnake?',prolsize,nsize
       !if( iout.eq.1) print*,'--- POSTSMOOTHING RESIDUUM --- on level',lvl,'/'
       ! POSTSMOOTHING
       !if( iout.eq.1) print*,x(1:10)

       if( iout.eq.1) then
          call MGbMVprod(grid,y,x,nsize)
          write(*,"(1X,A37,6X,I1,1X,A3,I3,es16.6)") &
               '----- POSTSMOOTHING RESIDUUM on level',lvl,'   ',-1,norm2(b-y)
       end if

       x(1:nsize) = x(1:nsize) + prol(1:prolsize)

       if( iout.eq.1) then
          call MGbMVprod(grid,y,x,nsize)
          write(*,"(1X,A37,6X,I1,1X,A3,I3,es16.6)") &
               '----- POSTSMOOTHING RESIDUUM on level',lvl,'  .',0,norm2(b-y)
       end if

       do  i=1,postsmooth,1
          !call bJacobi(x,b,nsize)
          !call bGs(x,b,nsize)
          call SMOOTHER(x,b,nsize)
          if( iout.eq.1) then
             call MGbMVprod(grid,y,x,nsize)
             write(*,"(1X,A37,6X,I1,1X,A3,I3,es16.6)") &
                  '----- POSTSMOOTHING RESIDUUM on level',lvl,'  /',i,norm2(b-y)
          end if
       end do

    else ! lvl == lowest

       ! "EXACT" SOLVER
       if( exact_solu.eq.'GMRES' ) then

          !call WriteMblock(grid%elem(1)%ILU(1) )

          call  Compute_MG_BlockDiagPrecond()

          !call WriteMblock(grid%elem(1)%ILU(0) )

          !stop

          !call gmres(nsize, x(1:nsize), b(1:nsize), restart*nloops, 1.E-10, &
          !     MGbMVprod, MGbMViLUprod, 30,  1.E-10, iout, iter, rezid, &
          !     not_conv)

          !iout = 0
          call gmres(nsize, x(1:nsize), b(1:nsize), restart*nloops, 5.E-1, &
                                !MGbMVprod, bMVnull, 30,  1.E-3, iout, iter, rezid, &
               grid, MGbMVprod, MGbMVprodDiag, 30,  1.E-3, iout, iter, rezid, &
               not_conv)
          !iout = 1

          if( iout.eq.1) then
             call MGbMVprod(y,x,nsize)
             write(*,"(1X,A31,12X,I1,1X,A3,I3,es16.6)") &
                  '---- "EXACT" SOL gmres on level',lvl,' V ',1,norm2(b-y)
          end if

          !stop

          !
          !PAUSE
          !
       elseif( exact_solu.eq.'bGS') then
          do  i=1,10,1
             call bGS(x,b,nsize)
             if( iout.eq.1) then
                !if( (i==1).or.(i==100) ) then
                call MGbMVprod(y,x,nsize)
                write(*,"(1X,A31,12X,I1,1X,A3,I3,es16.6)") &
                     '----- "EXACT" SOLUTION on level',lvl,' V ',1,norm2(b-y)
                !end if
             end if
          end do


       elseif( exact_solu.eq.'bJacobi') then
          do  i=1,10,1
             call bJacobi(x,b,nsize)
             if( iout.eq.1) then
                if( (i==1).or.(i==100) ) then
                   call MGbMVprod(y,x,nsize)
                   write(*,"(1X,A31,12X,I1,1X,A3,I3,es16.6)") &
                        '----- "EXACT" SOLUTION on level',lvl,' V ',1,norm2(b-y)
                end if
             end if
          end do

       elseif( exact_solu.eq.'UMFPACK' ) then
          call SolveViaUMFPACK(nsize,x,b)

       elseif( exact_solu.eq.'bJacobi100' ) then
          do  i=1,100,1
             call bJacobi(x,b,nsize)
             call MGbMVprod(y,x,nsize)
             if( (i==1).or.(i==100) ) then
                if( iout.eq.1) write(*,"(I5,F38.15)") i,norm2(b-y)
             end if
          end do

       elseif( exact_solu.eq.'bJacobi10') then
          do  i=1,10,1
             call bJacobi(x,b,nsize)
             call MGbMVprod(y,x,nsize)
             if( (i==1).or.(i==10) ) then
                if( iout.eq.1) write(*,"(I5,F38.15)") i,norm2(b-y)
             end if
          end do

       elseif( exact_solu.eq.'bJacobi1') then
          do  i=1,1,1
             call bJacobi(x,b,nsize)
             call MGbMVprod(y,x,nsize)
             if( (i==1).or.(i==1) ) then
                if( iout.eq.1) write(*,"(I5,F38.15)") i,norm2(b-y)
             end if
          end do

       elseif( exact_solu.eq.'ILU' ) then
          call MGbMViLUprod(y,b,nsize)
          x = y
          call MGbMVprod(y,x,nsize)
          if( (i==1).or.(i==100) ) then
             if( iout.eq.1) write(*,"(I5,F38.15)") i,norm2(b-y)
          end if

       else
          print*,'UNKNOWN solver type in LINMGLVL -- STOP'
          STOP
       endif

    endif


    deallocate(y)
    !print*,'####$$$',lvl, allocated(cor), allocated(prol), allocated(u), allocated(v)
    if( allocated(cor) ) deallocate(cor,res)
    if( allocated(prol)) deallocate(prol)
    if( allocated(u)) deallocate(u)
    if( allocated(v)) deallocate(v)


  END SUBROUTINE LINMGLVL_mod



  !> copy of elem%w -> elem%MGw
  SUBROUTINE MGcopyElem(elem)

    type(element), intent(inout) :: elem

    elem%MGw(1, MGv, 1:elem%dof * ndim) = elem%w(0, 1:elem%dof * ndim)

  END SUBROUTINE MGcopyElem



  !> restriction of elem%MGw
  SUBROUTINE RestrictElem(elem)
    type(element), intent(inout) :: elem
    integer :: i, k, k1, k2, deg, deg1, dof, dof1

    deg = elem%MGdeg
    deg1 = deg - 1  ! KONTROLA
    dof  = (deg  + 1)*(deg  + 2)/2
    dof1 = (deg1 + 1)*(deg1 + 2)/2

    k = 1
    do  i=1,ndim
       k1 = (i-1)*dof*ndim
       k2 = (k1 -1 + dof1)*ndim
       elem%MGw(elem%deg - deg1 + 1, MGv, k:k+dof1 -1 * ndim) &
            = elem%MGw(elem%deg - deg1, MGv, k1:k2)

       k = k + dof1*ndim
    end do

  END SUBROUTINE RestrictElem



  SUBROUTINE SolveVia10bJacobi(nsize,x,b)

    integer,intent(in):: nsize
    real,dimension(nsize),intent(inout) :: x
    real,dimension(nsize),intent(in) :: b

    call bJacobi(x,b,nsize)

  END SUBROUTINE SolveVia10bJacobi



  !> Block Jacobi method for linear problem (\f$ (A + \eta\M)x= b\f$
  !>
  !> Consider eq \f$ Ax=b \f$.
  !> \f$ A \f$ is decomposed on diagonal component \f$ D \f$ and remaining part \f$ R \f$, i.e. \f$ A=D+R \f$.
  !> \f$ Ax=b -> (D+R)x=b \ldots Dx=b-Rx \ldots Dx^{k+1}=b-Rx^k \ldots x^{k+1}=inv(D)(b-Rx^k) \f$
  SUBROUTINE bJacobi(xin, b, nsize)

    integer,intent(in) :: nsize
    real,dimension(1:nsize),intent(in):: b
    real,dimension(1:nsize),intent(inout):: xin
    real,dimension(1:nsize):: xout
    class(element), pointer :: elem, elem1
    integer :: i,j,k,MGis,MGis1,MGndof,MGndof1,MGdof,MGdof1,ki,kj,dof,dof1
    real, dimension(:,:),allocatable :: accuM
    integer :: maxsize

    maxsize = maxval(grid%elem(:)%MGdof)*ndim
    allocate(accuM(1:maxsize,1:maxsize))

    xout(:)=0.

    do  i=1,grid%nelem,1

       elem => grid%elem(i)
       MGdof = elem%MGdof
       MGndof = ndim * MGdof
       MGis = elem%MGncv
       dof = elem%dof

       ! off-diagonal blocks for i-th element
       do  j=1,elem%flen,1
          k = elem%face(neigh,j)

          if( k > 0) then
             elem1 => grid%elem(k)
             MGdof1 = elem1%MGdof
             MGndof1 = ndim * MGdof1
             MGis1 = elem1%MGncv
             dof1 = elem1%dof

             do  ki=0,ndim-1,1
                do  kj=0,ndim-1,1
                   xout(MGis+ki*MGdof : MGis+(ki+1)*MGdof-1) = xout(MGis+ki*MGdof : MGis+(ki+1)*MGdof-1) &
                        + matmul( &
                        elem%block(j)%Mb(ki*dof+1:ki*dof+MGdof, kj*dof1+1:kj*dof1+MGdof1), &
                        xin(MGis1+kj*MGdof: MGis1+(kj+1)*MGdof1-1) &
                        )
                end do
             end do
          end if
       end do ! j=1,elem%flen

       xout(MGis : MGis+MGndof-1) = b(MGis : MGis+MGndof-1) - xout(MGis : MGis+MGndof-1)

       ! INVERSE of diagonal block
       ! -- initialization of accuM
       accuM(:,:)=0.

       do  ki=0,ndim-1,1
          do  kj=0,ndim-1,1
             accuM(ki*MGdof+1:ki*MGdof+MGdof, kj*MGdof+1:kj*MGdof+MGdof) &
                  = elem%block(0)%Mb(ki*dof+1:ki*dof+MGdof, kj*dof+1:kj*dof+MGdof)
          end do
       end do

       if( eta /= 0.) then
          do  k = 0,MGndof-1,MGdof
             accuM(k+1:k+MGdof,k+1:k+MGdof) &
                  = accuM(k+1:k+MGdof,k+1:k+MGdof) + eta*elem%Mass%Mb(1:MGdof,1:MGdof)
          end do
       end if

       call MblockInverse(MGndof,accuM(1:MGndof,1:MGndof))

       xout(MGis:MGis+MGndof-1) = matmul(accuM(1:MGndof,1:MGndof),xout(MGis:MGis+MGndof-1))

    end do ! i=1,grid%nelem

    xin = xout

    deallocate(accuM)

  END SUBROUTINE bJacobi



  !> block Gauss-Seidel method ??
  SUBROUTINE bGS(x, b, nsize)
    integer,intent(in) :: nsize
    real, dimension(1:nsize), intent(in):: b
    real, dimension(1:nsize), intent(inout):: x
    real, dimension(1:nsize):: x1
    real, dimension(:),allocatable:: suma
    class(element), pointer :: elem, elem1
    integer :: i,j,k,MGis,MGis1,MGndof,MGndof1,MGdof,MGdof1,ki,kj,dof,dof1
    real, dimension(:,:),allocatable :: accuM
    integer :: maxsize

    maxsize = maxval(grid%elem(:)%MGdof)*ndim
    allocate(accuM(1:maxsize,1:maxsize))

    accuM(:,:)=0.
    x1(:) = 0.

    allocate(suma( ndim*grid%elem(1)%dof ))

    do  i=1,grid%nelem
       elem => grid%elem(i)
       MGdof = elem%MGdof
       MGndof = ndim * MGdof
       MGis = elem%MGncv
       dof = elem%dof

       if( size(suma,1)/=MGndof) then
          deallocate(suma)
          allocate(suma(MGndof))
       end if
       suma=0.

       ! off-diagonal blocks j<i
       do  j=1,elem%flen
          k = elem%face(neigh,j)

          if( (0<k).AND.(k<i)) then
             elem1 => grid%elem(k)
             MGdof1 = elem1%MGdof
             MGndof1 = ndim * Mgdof1
             MGis1 = elem1%MGncv
             dof1 = elem1%dof

             do  ki=0,ndim-1,1
                do  kj=0,ndim-1,1
                   suma(ki*MGdof+1 : ki*MGdof+MGdof) = suma(ki*MGdof+1 : ki*MGdof+MGdof) + &
                        matmul( &
                        elem%block(j)%Mb(ki*dof+1 : ki*dof+MGdof, kj*dof1+1 : kj*dof1+MGdof1), &
                        x1(MGis1+kj*MGdof: MGis1+(kj+1)*MGdof1-1) &
                        )
                end do
             end do
          end if
       end do ! j=1,elem%flen

       ! off-diagonal blocks j>i
       do  j=1,elem%flen
          k = elem%face(neigh,j)

          if( k>i) then
             elem1 => grid%elem(k)
             MGdof1 = elem1%MGdof
             MGndof1 = ndim * dof1
             MGis1 = elem1%MGncv
             dof1 = elem1%dof

             do  ki=0,ndim-1,1
                do  kj=0,ndim-1,1
                   suma(ki*MGdof+1 : ki*MGdof+MGdof) = suma(ki*MGdof+1 : ki*MGdof+MGdof) + &
                        matmul( &
                        elem%block(j)%Mb(ki*dof+1 : ki*dof+MGdof, kj*dof1+1 : kj*dof1+MGdof1), &
                        x(MGis1+kj*MGdof: MGis1+(kj+1)*MGdof1-1) &
                        )
                end do
             end do
          end if
       end do ! j=1,elem%flen

       suma = b(MGis : MGis+MGndof-1) - suma

       ! INVERSE of diagonal block
       ! -- initialization of accuM
       accuM(:,:)=0.

       do  ki=0,ndim-1,1
          do  kj=0,ndim-1,1
             accuM(ki*MGdof+1:ki*MGdof+MGdof, kj*MGdof+1:kj*MGdof+MGdof) &
                  = elem%block(0)%Mb(ki*dof+1:ki*dof+MGdof, kj*dof+1:kj*dof+MGdof)
          end do
       end do

       if( eta /= 0.) then
          do  k = 0,MGndof-1,MGdof
             accuM(k+1:k+MGdof,k+1:k+MGdof) &
                  = accuM(k+1:k+MGdof,k+1:k+MGdof) + eta*elem%Mass%Mb(1:MGdof,1:MGdof)
          end do
       end if

       call MblockInverse(MGndof,accuM(1:MGndof,1:MGndof))

       x1(MGis:MGis+MGndof-1) = matmul(accuM(1:MGndof,1:MGndof),suma)

    end do ! i=1,grid%nelem

    x=x1

    deallocate(accuM,suma)

  END SUBROUTINE bGS



  !> ILU preconditioning: \f$ b = (LU)^{-1}x \f$,  \f$ LU  \f$ is the incomplete LU block
  !> preconditioner having the same structure as matrix
  SUBROUTINE MGbMViLUprod(b,x,nsize)

    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    real, dimension(:), allocatable :: y
    class(element), pointer:: elem, elem1 ! one element
    type(Mblock) :: Loc

    integer :: i, ndof, is, j, i1, ndof1, is1, ii

    call InitMblock(Loc, maxval( grid%elem(:)%MGdof ) * ndim, maxval( grid%elem(:)%MGdof ) * ndim)

    allocate(y(1:nsize) )

    !! L solution
    do  i=1,grid%nelem,1
       elem => grid%elem(i)
       ndof = elem%MGdof * ndim
       is = elem%MGncv

       y(is: is+ndof-1) = x(is: is+ndof-1)

       do  j=1,elem%flen
          i1 = elem%face(neigh,j)
          if( i1 > 0 .and. i1 < i) then
             elem1 => grid%elem(i1)
             ndof1 = elem1%MGdof * ndim
             is1 = elem1%MGncv

             y(is: is+ndof-1) = y(is: is+ndof-1) &
                  - matmul(elem%XXX(j)%Mb(1:ndof, 1:ndof1), y(is1: is1+ndof1-1) )

             if( is1 > is) print*,'Problem MGbMViLUprod!: L solution',is,is1, i,i1

          end if
       end do
    end do

    !! U solution
    do  ii=1,grid%nelem
       i = grid%nelem - ii + 1

       elem => grid%elem(i)
       ndof = elem%MGdof * ndim
       is = elem%MGncv

       do  j=1,elem%flen
          i1 = elem%face(neigh,j)
          if( i1 > i) then
             elem1 => grid%elem(i1)
             ndof1 = elem1%MGdof * ndim
             is1 = elem1%MGncv

             y(is: is+ndof-1) = y(is: is+ndof-1) &
                  - matmul(elem%XXX(j)%Mb(1:ndof, 1:ndof1), b(is1: is1+ndof1-1) )

             if(is1 < is) print*,'Problem MGbMViLUprod! U solution',is,is1, i,i1

          end if
       end do

       if( ndof .ne. size(Loc%Mb,1)) then
          deallocate (Loc%Mb)
          call InitMblock(Loc, ndof, ndof)
       end if

       Loc%Mb(1:ndof,1:ndof) = grid%elem(i)%XXX(0)%Mb(1:ndof,1:ndof)
       call MblockInverse(ndof, Loc%Mb)

       b(is: is+ndof-1) = matmul(Loc%Mb(1:ndof,1:ndof), y(is: is+ndof-1) )

    end do


    deallocate (Loc%Mb)
    deallocate(y)

  END SUBROUTINE MGbMViLUprod


  !> pMG-restriction operator.
  !> Vector x(1:nsize) is restricted on vector xout(1:outsize), xout(outsize+1:nsize) is garbage.
  SUBROUTINE MGrestVec(x, nsize, xout, outsize)

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

    do  i=1,grid%nelem
       elem => grid%elem(i)

       deg = elem%MGdeg		! actual polynomial degree on current element (will be reduced)
       dof = elem%MGdof		! #DOF on element
       ndof = ndim*dof		! size of problem on the current element, = \frac{d+2}{d!} \Pi_{j=1}^{d}(deg+j)
       is = elem%MGncv		! position in state vector for current element and MG-level

       if( deg == MGminP) then	! we cannot decrease pol.degree on current element, because we reached minimum pol.degree
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
          pause

          !xout(1 : MGnsize+ndof ) = x( is : is+ndof-1 ) !FIXME??
          xout(MGnsize+1 : MGnsize+ndof ) = x( is : is+ndof-1 )

          MGnsize = MGnsize + ndof	! yet last record in xout

       else ! we decrease polynomial degree

          !MGndof = ndof*deg/(deg+d)	! #DOF is reduced
          !MGdeg = deg-1		! aprox. degree is reduced= \frac{d+2}{d!} \Pi_{j=1}^{d}(deg+j)
          !MGdof = MGndof/ndim
          if( nbDim==2) then
             MGdeg = deg-1
             MGdof = (MGdeg+1)*(MGdeg+2)/2
             MGndof = ndim*MGdof
          else
             print*,'Not implemented for nbDim=',nbDim,'--STOP'
             STOP
          end if

          ! update stored information for current element
          elem%MGdof = MGdof
          elem%MGdeg = MGdeg
          elem%MGncv = MGnsize+1	!first record (in xout) for current element

          ! copying RESTRICTED relevant parts of vector:
          do  k=1,ndim,1
             xout( MGnsize+1 : MGnsize+MGdof ) = x( is + (k-1)*dof : is + (k-1)*dof + MGdof -1 )
             MGnsize = MGnsize + MGdof
          end do
       end if
    end do ! i=1:grid%nelem

    outsize = MGnsize

  END SUBROUTINE MGrestVec



  SUBROUTINE MGprolVec_renew(xin, insize, xout, outsize, degrees)
    real,dimension(1:insize),intent(in):: xin
    integer,intent(in):: insize			! old problem size
    integer, intent(inout):: outsize		! on input exceptet size of prolongation, on output real size
    real,dimension(1:outsize),intent(out):: xout
    integer,dimension(grid%nelem),intent(in):: degrees

    ! INNER
    class(element), pointer :: elem
    integer :: i, k, deg, dof, ndof, is, MGdeg, MGdof, MGndof, Mgnsize
    !integer :: MGmaxP = 1			! maximal polynomial degree of aproximation used in pMG algorithm

    MGnsize = 0
    xout(:) = 0.

    do  i=1,grid%nelem,1
       elem => grid%elem(i)

       deg = elem%MGdeg		! actual polynomial degree on current element (will be raised)
       dof = elem%MGdof		! #DOF on element
       ndof = ndim*dof		  ! size of problem on the current element, = \frac{d+2}{d!} \Pi_{j=1}^{d}(deg+j)
       is = elem%MGncv		  ! position in state vector for current element and MG-level

       if( deg == degrees(i)) then	! maximal pol. degree is reached
          ! MGdof = dof		! unchanged
          ! MGdeg = deg		! unchanged
          ! MGndof = ndof 		! unchanged, = ndim*dof

          ! update stored information for current element
          ! elem%MGdof = dof
          ! elem%MGdeg = deg
          elem%MGncv = MGnsize+1	! first record (in xout) for current element

          ! copying WHOLE relevant part of vector:
          xout( MGnsize+1 : MGnsize+ndof ) = xin( is : is+ndof-1 )

          MGnsize = MGnsize + ndof	! yet last record in xout
       else
          if( nbDim==2) then
             MGdeg = deg+1
             MGdof = (MGdeg+1)*(MGdeg+2)/2
             MGndof = ndof*MGdof
          else
             print*,'Not implemented for nbDim=',nbDim,'--STOP'
             STOP
          end if
          ! update MG values
          elem%MGdeg = MGdeg
          elem%MGdof = MGdof
          elem%MGncv = MGnsize+1	!first record (in xout) for current element

          do  k=1,ndim,1
             xout( MGnsize+1 : MGnsize+dof ) = xin( is+(k-1)*dof : is+k*dof - 1 )
             MGnsize = MGnsize + MGdof
          end do
       end if
    end do !i=1:grid%elem

    if( MGnsize/=outsize) then
       print*,'Something is WRONG: sizes after prolongation are not equal -- STOP'
       STOP
    end if

    outsize = MGnsize

  END SUBROUTINE MGprolVec_renew


  !> pMG prolongation operator
  !> Vector x(1:nsize) is prolongated on vector xout(1:state%nsize), xout(outsize+1:state%nsize) is garbage
  subroutine MGprolVec_new(x, nsize, xout, outsize, degrees)
    ! INPUT
    real, dimension(1:nsize), intent(in) :: x
    integer, intent(in) :: nsize			! old problem size
    ! OUTPUT
    ! tu je problem, malo by tam byt outsize namiesto 1:state%nsize
    integer, intent(inout) :: outsize		! new problem size, i.e. restricted problem size
    real, dimension(1:outsize), intent(out) :: xout
    integer :: MGnsize
    integer,dimension(outsize),intent(in) :: degrees
    ! INNER
    class(element), pointer :: elem, elem1
    integer :: i, k, deg, dof, ndof, is, MGdeg, MGdof, MGndof, maxdeg
    integer :: l, Qnum, Qdof
    real, dimension(:,:), pointer :: phi0 ! the test functions
    real, dimension(:,:), allocatable :: wi, Fx
    real, dimension(:,:,:,:), allocatable :: Dw
    !integer :: d				! global dimension 2 or 3
    !integer :: ndim = 4
    !integer :: MGmaxP = 1			! maximal polynomial degree of aproximation used in pMG algorithm

    xout(:) = 0.
    MGnsize = 0
    !d = nbDim



    ! TEST OF A NEW RECONSTRUCTION
    goto 10

    print*,'! TEST OF A NEW RECONSTRUCTION'
    do i=1, grid%nelem
       elem => grid%elem(i)
       deg = elem%MGdeg
       dof = elem%MGdof

       is = elem%MGncv		! position in state vector for current element and MG-level
       allocate(elem%wS(1:ndim, 1:dof) )

       do  k=1,ndim,1
          elem%wS(k, 1:dof) = x( is+(k-1)*dof : is+k*dof - 1 )

          !!write(*,'(a3, 4i5, 30es12.4)') 'WQS',i,k, is, dof, elem%wS(k, 1:dof)

       enddo

       ! vypis
       ! Qnum = state%space%Qdeg(deg+1, 1)
       ! Qdof = state%space%V_rule(Qnum)%Qdof

       ! phi0 => state%space%V_rule(Qnum)%phi(1:dof, 1:Qdof)
       ! allocate( wi(k, 1:Qdof) )
       ! do k=1, ndim
       !    wi(k, 1:Qdof) = matmul(elem%wS(k, 1 : dof), phi0(1:dof, 1:Qdof) )
       ! enddo

       ! allocate(Fx(1:Qdof, 1:2) )
       ! call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )
       ! do l=1,Qdof
       !    write(70+i, *) Fx(l,1:2), wi(:, l)
       ! enddo


       ! deallocate(wi, Fx)
       ! konec vypis
    enddo

    do i=1, grid%nelem
       elem => grid%elem(i)
       !deg = elem%MGdeg
       !dof = elem%MGdof

       MGdeg = elem%MGdeg+1
       MGdof = (MGdeg+1)*(MGdeg+2)/2
       Qnum = state%space%Qdeg(deg+1, 1)


       allocate(Dw(1:ndim,0:0, 0:0, 1: MGdof) )
       call  LeastSquareInterpolationL2(elem, ndim, .true. , Qnum, 0, MGdof, &
            Dw(1:ndim,0:0, 0:0, 1: MGdof), 0  )

       ! vypis
       ! do  k=1,ndim,1
       !    write(*,'(a3, 4i5, 30es12.4)') '##S',i,k, is, dof, Dw(k, 0,0, 1:MGdof)

       ! enddo
       ! Qdof = state%space%V_rule(Qnum)%Qdof

       ! phi0 => state%space%V_rule(Qnum)%phi(1:MGdof, 1:Qdof)
       ! allocate( wi(k, 1:Qdof) )
       ! do k=1, ndim
       !    wi(k, 1:Qdof) = matmul(Dw(k, 0, 0, 1 : MGdof), phi0(1:MGdof, 1:Qdof) )
       ! enddo

       ! allocate(Fx(1:Qdof, 1:2) )
       ! call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )
       ! do l=1,Qdof
       !    write(80+i, *) Fx(l,1:2), wi(:, l)
       ! enddo


       ! deallocate(wi, Fx)
       ! end vypis

       deallocate(Dw)
    enddo


    ! DEALLOCATE elem%wS
    do i=1, grid%nelem
       elem => grid%elem(i)
       deallocate(elem%wS)
    enddo

    ! END OF THE TEST OF A NEW RECONSTRUCTION

    10 continue

    do  i=1,grid%nelem,1
       elem => grid%elem(i)

       maxdeg = elem%deg		! ?default polynomial degree
       deg = elem%MGdeg		! actual polynomial degree on current element (will be raised)
       dof = elem%MGdof		! #DOF on element
       ndof = ndim*dof		! size of problem on the current element, = \frac{d+2}{d!} \Pi_{j=1}^{d}(deg+j)
       is = elem%MGncv		! position in state vector for current element and MG-level

       !    if (deg == maxdeg) then	! maximal pol. degree is reached
       if( deg == degrees(i)) then	! maximal pol. degree is reached
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

          !MGdeg = deg+1
          !MGndof = ndof*(MGdeg+d)/MGdeg
          !MGdof = MGndof/ndim
          if( nbDim==2) then
             MGdeg = deg+1
             MGdof = (MGdeg+1)*(MGdeg+2)/2
             MGndof = ndof*MGdof
          else
             print*,'Not implemented for nbDim=',nbDim,'--STOP'
             STOP
          end if

          elem%MGdeg = MGdeg
          elem%MGdof = MGdof
          elem%MGncv = MGnsize+1	!first record (in xout) for current element

          do  k=1,ndim,1
             xout( MGnsize+1 : MGnsize+dof ) = x( is+(k-1)*dof : is+k*dof - 1 )
             MGnsize = MGnsize + MGdof
          end do
       end if
    end do !i=1:grid%elem

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
       call bJacobi(MGX,MGB,state%nsize)
       call MGbMVprod(grid,MGY,MGX,state%nsize)
       print*,i,'resi fine = ',norm2(MGB-MGY)
    end do

    print*,MGX(1:10)
    call MGbMVprod(grid,MGY,MGX,state%nsize)
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

    call MGbMVprod(grid, MGY2,MGX2,outsize)
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
       call bJacobi(MGX2,MGB2,outsize)
       call MGbMVprod(grid,MGY2,MGX2,outsize)
       print*,i,'resi cors = ',norm2(MGB2-MGY2)
    end do


    print*,'MGX2'
    write(*,'(3F10.5)') MGX2(1:20)
    call MGbMVprod(grid,MGY2,MGX2,outsize)
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
       call bJacobi(MGX,MGB,dummysize)
       call MGbMVprod(grid,MGY,MGX,dummysize)
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

    call MGbMVprod(grid,MGY3,MGX3,outsize2)
    print*,'MGY3',MGY3(1:20)
    print*,'MGX3',MGX3(1:20)
    print*,size(MGX3),size(MGB3),size(MGY3),outsize2
    print*,norm2(MGX3)
    print*,norm2(MGY3)
    print*,'xxxxxx resi cors = ',norm2(MGB3(1:outsize2)-MGY3(1:outsize2)), norm2(MGB3(1:outsize2)), norm2(MGY3(1:outsize2))


    print*,' *** LEVEL',grid%elem(1)%MGdeg,'***'

    do i=1,3,1
       !	  	print*,i,'solu norm = ',norm2(MGX)
       call bJacobi(MGX3,MGB3,outsize2)
       call MGbMVprod(grid, MGY3,MGX3,outsize2)
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

    call MGbMVprod(grid,MGY4,MGX4,outsize3)
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
       call bJacobi(MGX4,MGB4,outsize3)
       call MGbMVprod(grid,MGY4,MGX4,outsize3)
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
    call MGbMVprod( grid,MGY3,MGX3,outsize2)
    print*,i,'resi cors = ',norm2(MGB3-MGY3)


    print*,'POSTresiduum1'
    call MGbMVprod(grid,MGY3,MGX3(1:outsize2) - DUMMY(1:dummysize),outsize2)
    print*,i,'resi cors = ',norm2(MGB3-MGY3)

    print*,'POSTresiduum2'
    call MGbMVprod(grid,MGY3,MGX3(1:outsize2) + DUMMY(1:dummysize),outsize2)
    print*,i,'resi cors = ',norm2(MGB3-MGY3)



    MGX3(1:outsize2) = MGX3(1:outsize2) - DUMMY(1:dummysize)


    print*,'POSTresiduum0'
    call MGbMVprod(grid,MGY3,MGX3,outsize2)
    print*,i,'resi cors = ',norm2(MGB3-MGY3)


    pause

    do i=1,3,1
       !print*,i,'solu norm = ',norm2(MGX)
       call bJacobi(MGX3,MGB3,outsize2)
       call MGbMVprod(grid,MGY3,MGX3,outsize2)
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
       call bJacobi(MGX2,MGB2,outsize)
       call MGbMVprod(grid,MGY2,MGX2,outsize)
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
       call bJacobi(MGX,MGB,dummysize)
       call MGbMVprod(grid,MGY,MGX,dummysize)
       print*,i,'resi cors = ',norm2(MGB-MGY)
    end do




    call cpu_time(Tend)

    print*,Tend-Tstart,' s'

  end subroutine pMGSolver




  !> evaluation of block diagonal preconditioner,
  !> LU decomposition of elem%block(0) using LAPACK
  subroutine Compute_MG_BlockDiagPrecond()
    use matrix_oper_int
    class(element), pointer:: elem
    integer :: i, ki, kj, kki, kkj, Mki, Mkj,dof,   MGdof, MGndof

    print*,'###', eta, 1./state%time%tau(1)

    do i=1,grid%nelem
       elem => grid%elem(i)
       dof = elem%dof
       MGdof = elem%MGdof
       MGndof = MGdof * ndim

       do  kki=0,ndim-1,1
          do  kkj=0,ndim-1,1

             ki = kki * dof
             kj = kkj * dof

             Mki = kki * MGdof
             Mkj = kkj * MGdof

             ! flux matrix
             elem%ILU(0)%Mb(Mki+1: Mki + MGdof,Mkj+1: Mkj + MGdof ) &
                  = elem%block(0)%Mb(ki+1: ki + MGdof,kj+1: kj + MGdof )
             !print*,'###',Mki+1, Mki + MGdof,':',Mkj+1, Mkj + MGdof,'||', &
             !     ki+1,  ki + MGdof,':', kj+1,  kj + MGdof


             ! adding of mass matrix
             if (eta /= 0. .and. ki == kj ) then
                elem%ILU(0)%Mb(Mki+1: Mki+MGdof, Mki+1: Mki+MGdof) = &
                     elem%ILU(0)%Mb(Mki+1 : Mki+MGdof, Mki+1: Mki+MGdof) &
                     + eta * elem%Mass%Mb(1:MGdof,1:MGdof)
             end if
          enddo
       enddo
       call MblockInverse(MGndof, elem%ILU(0)%Mb(1:MGndof, 1:MGndof))

    enddo
  end subroutine Compute_MG_BlockDiagPrecond

end module pMultiGrid2
