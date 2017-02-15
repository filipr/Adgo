module main_data
  use geometry
  use state_mod
  use mesh_mod

  implicit none

  real(8), parameter :: pi=3.1415926535d0
  !> computational grids : paralelizable
  !type(mesh) :: grid
  !type(mesh), pointer :: grid, gridN, gridS, gridD

  !> problem data : non-paralelizable


  type( solver_state ) :: state

  !> coordinates of nodes on curved boundary part
  type(curved_boundary) :: curvbound

  !> matrix shape for conforming FEM only
  type(MatrixShape) :: Mshape

   !> global file variables
  class (mesh), allocatable :: grid, gridN, gridS, gridD



  character(len=50) :: gridfile     ! file with mesh
  character(len=50) :: solfile      ! solution with basis coefficients
  character(len=50) :: rsolfile     ! solution in Lagrange nodes
  character(len=50) :: convfile     ! hisotry of computation
  character(len=50) :: prof_file    ! nodes lying on curved parts of boundary


  contains

end module main_data

module helpers
  use main_data
  implicit none


  public:: saveoct
  public:: Write_state_settings


contains
  subroutine saveoct(name,mat,uplo)
  character(*),intent(in):: name
  character(*),intent(in):: uplo
  real:: mat(:,:)
  integer:: i,j
  open(27,file=name,action='write')
  if (uplo == 'N') then
    do i = 1,size(mat,1)
      write(27,'(1000E26.16)') mat(i,:)
    end do
  else if (uplo == 'L') then
    do i = 1,size(mat,1)
      write(27,'(1000E26.16)') mat(i,:i),(0d0,j=i+1,size(mat,2))
    end do
  else if (uplo == 'U') then
    do i = 1,size(mat,1)
      write(27,'(1000E26.16)') (0d0,j=1,i-1),mat(i,i:)
    end do
  else if (uplo == 'SL') then
    do i = 1,size(mat,1)
      write(27,'(1000E26.16)') mat(i,:i),mat(i+1:,i)
    end do
  else if (uplo == 'SU') then
    do i = 1,size(mat,1)
      write(27,'(1000E26.16)') mat(:i-1,i),mat(i,i:)
    end do
  end if
  close(27)
end subroutine saveoct

  !> write the setting of state
  subroutine Write_state_settings( )
    integer :: i

    write(*,*) '######################### Write_state_settings( ) '
    write(*,*) state%model%kappa , state%model%kappa1, state%model%Pr

  ! scalar case
    if(state%modelName == 'scalar' .or.state%modelName == '2eqs') &
         write(*,*)  state%model%icase,  state%model%param1

  ! default degree of approximation of the solution
   write(*,*)  state%space%deg
   write(*,*) state%time%time_method

   write(*,*) state%time%time_method, state%time%deg

   write(*,*) state%time%cn , state%time%stdgm , state%time%tdg,  state%time%stdgm
   write(*,*) state%time%deg

   write(*,*) 'Newton ', state%nlSolver%max_iter, state%nlSolver%tol, state%nlSolver%tol2, &
       state%nlSolver%min_update

  ! penalty parameter
   write(*,*) state%space%m_IPG, state%space%C_W   !!!!!, state%space%pen_deg


   write(*,*) state%type_IC

   write(*,*) 'max_iter' , state%time%maxiter


   write(*,*) state%space%adapt%max_adapt_level, state%space%adapt%adapt_method

  ! stopping tolerance for steady-state reziduum
   write(*,*) 'SSres',state%conv_rez, state%EcD_tol, state%EcL_tol, state%EcM_tol !

  ! stopping final time
   write(*,*) 'TTime',state%time%FinTime
   write(*,*) 'taus',state%time%tau_new,  state%time%tau_fixed, state%time%tau_fixed_size

   write(*,*)'Mtol', state%linSolver%tol, state%linSolver%tol_fixed, state%MGsolver

   write(*,'(a31,i7,a7,e11.4,a11,e11.4)') &
        ' # Stopping criteria: max_iter=',state%time%maxiter, &
        ', time=', state%time%FinTime,', reziduum=',state%conv_rez

  ! output time
   write(*,*) '# Output files with time equidistance  = ', state%time%OutTime

  ! reading if init number of sol* files
  write(*,*) '# Initial values of sol* file = ',state%isol


  ! maximal CFL number
  write(*,*) '# Maximum CFL number = ',abs(state%time%CFL)

  ! tolerance for BDF
  write(*,*) '# Tolerance for adaptive BDF = ',state%time%BDFtol, state%time%tau_choice

  print*,'# Reynolds number/viscosity = ',state%model%Re,'/',state%model%Re1

  ! reading of boundary conditions
  do i=1,state%numBC

     write(*,*) state%BC(i)%ibc, state%BC(i)%inout, state%BC(i)%w(1:ndim)
     write(*,*) state%BC(i)%ibc, state%BC(i)%inout, state%BC(i)%ww(1:ndim)
  enddo

  if(ndim >= 4 .and. ndim <= nbDim + 4) then
     write(*,'(a38,6es9.2)')' # Far field: rho, v, alpha, p, M, E: ', &
          state%rho_infty, state%v_infty, state%alpha_infty/3.1415926*180, &
          state%p_infty, state%v_infty/(state%model%kappa*state%p_infty/state%rho_infty)**0.5, &
          state%p_infty/state%model%kappa1-0.5*state%rho_infty*state%v_infty**2.
     write(*,*)
  endif

  write(*,'(a30,4es9.2)') '# Stabilization parameters  :',  &
       state%ST_Vp, state%ST_Vc, state%ST_Ep, state%ST_Ec

  write(*,*) state%space%adapt%tol_max, state%space%adapt%tol_min, state%space%adapt%Lq

  write(*,*)'# ---------------------------------------------------------'


  end subroutine Write_state_settings



end module helpers
