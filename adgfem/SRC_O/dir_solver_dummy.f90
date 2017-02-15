module dir_solver
! this is dummy module, which shoud replace module stored in >dir_solver.f90<
! whenever no UMFPACK

public:: SolveViaUMFPACk
public:: umfpack_wrapper
public:: test_umfpack_simple
public:: test_umfpack_DG

private:: MSSG

CONTAINS

subroutine MSSG(text)
character(*),intent(in) :: text

print*,'--- YOU TRY TO CALL UMFPACK SUBROUTINE, BUT CODE WAS BUILT WITH "DUMMY" SUBROUTINE >> ',text,' <<'
print*,'--- STOP'
STOP

end subroutine MSSG


subroutine SolveViaUMFPACK(nsize,x,b)
integer,intent(in) :: nsize
real,dimension(:),intent(in) :: x
real,dimension(:),intent(in) :: b

call MSSG( 'SolveViaUMFPACK' )
end subroutine SolveViaUMFPACK


subroutine umfpack_wrapper(Ax,Ai,Ap,x,b,iout)
! input: NxN matrix in Harwell/Boeing format 0-base indexing, solution, RHS
real(kind=8),dimension(:),intent(in)    :: Ax
integer,dimension(:),intent(in)         :: Ai, Ap
real(kind=8),dimension(:),intent(in)    :: b
real(kind=8),dimension(:),intent(inout) :: x
integer :: iout

call MSSG( 'umfpack_wrapper' )
end subroutine umfpack_wrapper


subroutine test_umfpack_simple
call MSSG( 'test_umfpack_simple' )
end subroutine test_umfpack_simple


subroutine test_umfpack_DG
call MSSG( 'test_umfpack_DG' )
end subroutine test_umfpack_DG


end module dir_solver
