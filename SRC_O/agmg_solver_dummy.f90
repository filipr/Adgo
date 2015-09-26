module agmg_solver

private:: MSSg
public:: SolveViaAGMG
public::  dagmg

CONTAINS

subroutine MSSG(text)
character(*),intent(in) :: text

print*,'--- YOU TRY TO CALL AGMG SUBROUTINE, BUT CODE WAS BUILT WITH "DUMMY" SUBROUTINE >> ',text,' <<'
print*,'--- STOP'
STOP

end subroutine MSSG


subroutine dagmg(nsize,Ax,Ai,Ap,b,x,ijob,iprint,nrest,iter,tol)
integer,intent(in) :: nsize
real,dimension(:),allocatable :: Ax
integer,dimension(:),allocatable :: Ai, Ap
real,dimension(:),intent(inout) :: x,b
integer :: iter,iprint,ijob,nrest
real:: tol

call MSSG('dagmg')

end subroutine dagmg


subroutine SolveViaAGMG(nsize,x,b)
integer,intent(in) :: nsize
real,dimension(:),intent(inout) :: x,b
integer :: iter,tol,iprint,ijob

call MSSG('SolveViaAGMG')

end subroutine SolveViaAGMG

end module agmg_solver
