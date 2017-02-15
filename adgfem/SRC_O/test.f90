MODULE test
use main_data
use geometry
use paramets
use problem_oper
!use pMultiGrid2
use lin_solvers

implicit none

private:: FictGrid
public:: RunMGTest



CONTAINS



SUBROUTINE FictGrid(NELEM,DEG,X)
  integer,intent(in):: NELEM,DEG
  real,dimension(:),allocatable,intent(inout):: X
  integer:: DOF,NCV
  integer:: i,j,k
  !
  DOF = (DEG+1)*(DEG+2)/2
  NCV = 1
  grid%nelem = NELEM
  !
  allocate( X(ndim*NELEM*DOF) )
  !
  do  i=1,NELEM,1
      grid%elem(i)%ncv = NCV
      grid%elem(i)%deg = DEG
      grid%elem(i)%dof = DOF
      do  j=1,ndim,1
          X(NCV:NCV+DOF-1) = (/(k+(i-1)*DOF,k=1,DOF,1)/)
          NCV = NCV+DOF
      end do
  end do
  !
  !print*,X
  !
  !print*,'--STOP'
  !STOP
END SUBROUTINE FictGrid



SUBROUTINE RunMGTest
  integer:: nelem,deg,dof,nsize,outsize,nsize2
  real,allocatable,dimension(:):: xin,xou,xii
  integer,allocatable,dimension(:):: degrees

  !ndim=1
  nelem = 3
  deg = 3
  dof = (deg+1)*(deg+2)/2
  nsize = ndim*nelem*dof

  allocate( degrees(nelem) )
  degrees(:) = deg

  call FictGrid(nelem,deg,xin)
  call initMG()

  if( size(xin,1) /=nsize) then
      print*,'Nespravna velkost --STOP'
  end if

  allocate( xou(nsize) )

  call MGrestVec(xin,nsize,xou,outsize)
  print*,nsize
  write(*,"(10F8.1)") xin
  print*,''
  print*,outsize
  write(*,"(6F8.1)") xou
  !
  allocate( xii(nsize) )
  nsize2=nsize
  call MGprolVec_renew(xou, outsize, xii, nsize2, degrees)
  print*,nsize, nsize2
  write(*,"(10F8.1)") xii
  !
  !print*,'--STOP'
  !STOP
  !
END SUBROUTINE RunMGTest


SUBROUTINE RunILUSolu(x,b,nsize)
  real,dimension(nsize),intent(in):: b
  real,dimension(nsize),intent(inout):: x
  integer:: nsize
  real,dimension(nsize)::y

  call ComputeBlockILUPrecond()
  !print*,Newton%b(1:5)
  !print*,Newton%x(1:5)
  call bMViLUprod(x,b,nsize)
  !print*,Newton%b(1:5)
  !print*,Newton%x(1:5)
  call bMVprod(y,x,nsize)
  print*,'Residuum in iLU solver:',norm2(b-y)
  !
  !print*,'--STOP'
  !STOP
  !
END SUBROUTINE RunILUSolu

END MODULE test
