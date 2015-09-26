!> conforming P_1 FEM for scalar

module conforming_FEM
  use gmres_solver
  use main_data
  use model_oper
  use problem_oper
  use inviscid_fluxes
  use matrix_oper_int
  use euler_problem
  use umfpack_iterface

  implicit none

  public:: ConformingFEM
  public:: ComputeStiffFEM
  public:: ComputeMassFEM
  public:: ComputeConvectionFEM
  public:: ComputeReactionFEM
  public:: ComputeSourceTerms
  public:: ComputeRHSmass
  public:: ComputeDphi
  public:: ComputePhi
  public:: SetMshape
  public:: WriteMshape
  public:: WriteMatrix
  public:: InitDirichletBC
  public:: SetDirichletBC

contains

  subroutine ConformingFEM(convfile)
    type(MatrixShape) :: Fshape  ! full Mshape with trivial lines !!!
    character(len=*), intent(in) :: convfile
    real, dimension(:,:), allocatable :: phi
    real, dimension(:,:,:), allocatable :: Dphi
    real, dimension(:), pointer :: A
    real, dimension(:), allocatable :: Mass, Flux
    real, dimension(:), allocatable :: RHS, RHSmass, rez, x, b
    real, dimension(:,:), allocatable :: w
    real, dimension(:,:), allocatable :: FullM
    integer, dimension(:), allocatable :: ibnd
    character(len=1) :: ch
    character(len=3) :: yes_conv ! convergence?
    character(len=50) :: command_name, tri_name, sol_name, exa_name, err_name, est_name
    integer :: it_hours, it_min, it_sec
    integer :: npoin, nelem, i,j
    real ::  t_sta, t_end, norm_res
    integer :: iter, ifile
    integer :: nit, nvec, iout, not_converge
    real :: norm8, err8, new_err, lognew_err, logdiff, errL2, errH1, logerr, norm, normL2, normH1

    ifile = 10

    ! GMRES settings
    nvec = 50
    nit = 10*nvec
    iout = 0

    print*
    print*,' # Computation starts '
    print*

    npoin = grid%npoin
    nelem = grid%nelem

    call SetMshape(Fshape, grid)
    call WriteMshape(Fshape)

    allocate(phi(1:3, 1:3) )
    call ComputePhi(phi)

    allocate(Dphi(1:nelem, 1:3, 1:nbDim) )
    call ComputeDphi(grid, Dphi)


    allocate(Mass(1:Fshape%nonzero) )
    allocate(Flux(1:Fshape%nonzero) )
    Mass(:) = 0.
    Flux(:) = 0.

    !!!allocate(FullM(1:Fshape%nsize, 1:Fshape%nsize) )


    allocate(state%A(1:Fshape%nonzero) )
    A => state%A(1:Fshape%nonzero)

    ! diffusive term
    if(state%model%Re > 0.) then
       call ComputeStiffFEM(grid, Fshape, Dphi, Flux)
       Flux(:) = Flux(:) / state%model%Re
    endif

    !write(*,'(a6,10es12.4)') 'diff:',Flux(1:7)

    ! convective term
    call ComputeConvectionFEM(grid, Fshape, Dphi, phi, Flux)

    !write(*,'(a6,10es12.4)') 'conv:',Flux(1:7)

    ! reactive term
    call ComputeReactionFEM(grid, Fshape, phi, Flux)

    !write(*,'(a6,10es12.4)') 'reac:',Flux(1:7)

    ! mass matrix
    call ComputeMassFEM(grid, Fshape, phi, Mass)


    allocate(RHS(1:npoin) )
    allocate(RHSmass(1:npoin) )
    allocate(w(0:state%time%deg+1, 1:npoin) )
    allocate(rez(1:npoin) )
    allocate(ibnd(1:npoin) )

    allocate(x(1:npoin), b(1:npoin) )

    call InitSolutionVector(grid, w(0, 1:npoin) )

    ! setting of the Dirichlet boundary conditions
    call InitDirichletBC(grid, ibnd(1:npoin) )

    do i=1,state%time%deg+1
       w(i,1:npoin) = w(0,1:npoin)
    enddo

    ! variangle t_sta is overwritten when solution is written on output
    call cpu_time(t_sta)

    ! iteartive loop
    do iter = 1,state%time%maxiter

       ch = ' '
       ! solution at this iteration whould be saved?
       if(state%time%OutTime > 0.  .and. state%space%adapt%adapt_level >= 0 ) then
          if(state%timeprn + state%time%tau(1) >= state%time%OutTime &
            .or. iter == state%time%maxiter )  ch = '*'
       endif


       call PrepareNewTimeStep(iter)

       do i=1,state%time%deg_actual  !!! EBDF
          j = state%time%deg_actual +1  - i  !! EBDF
          w(j+1, 1:npoin ) = w(j,1 :npoin )
       enddo
       w(1, 1:npoin ) = w(0, 1:npoin )

       !write(38,'(a5,2i5,12es12.4)') '###',iter, state%time%deg_actual, &
       !     w(0:4,18 )
       state%time%ctime = state%time%ttime + state%time%tau(1)

       ! seting of coefficients \alpha, \beta for ABDF
       call SetABDF(state%BDF, state%time%deg_actual, state%time%tau(1:state%time%deg_actual+1) )


       ! total matrix
       A(1:Fshape%nonzero) = state%time%alpha(0)/state%time%tau(1) * Mass(1:Fshape%nonzero) &
            + Flux(1:Fshape%nonzero)


       ! computation of the BDF mass vector
       !write(*,'(a7,20es12.4)') 'w_0 :',w(0,1:5)
       w(0,1:npoin) = 0.
       do i=1,state%time%deg
          !write(*,'(a7,20es12.4)') 'w_0 :',w(0,1:5), state%time%alpha(i)

          w(0,1:npoin) = w(0,1:npoin) - state%time%alpha(i) * w(i,1:npoin)
       enddo

       RHSmass(:) = 0
       call ComputeRHSmass(grid,  phi, w(0, 1:npoin), RHSmass(1:npoin) )
       !write(*,'(a7,20es12.4)') 'w_0 :',w(0,1:5)
       !write(*,'(a7,20es12.4)') 'RhsM:',RHSmass(1:5)


       ! RHS - source terms
       RHS(:) = 0
       call ComputeSourceTerms(grid,  phi, RHS(1:npoin) )


       ! total RHS
       RHS(1:npoin) = RHS(1:npoin) + RHSmass(1:npoin)/ state%time%tau(1)

       !write(*,'(a6,10es12.4)') 'A:',A(1:7)
       !write(*,'(a6,10es12.4)') 'RHS:',RHS(1:7)


       ! initial guess of the solution
       w(0,1:npoin) = w(1,1:npoin)

       ! setting of the Dirichlet boundary conditions
       call SetDirichletBC(grid,  w(0,1:npoin) )


       call RemoveTrivialElements(Fshape, A, w(0,:), RHS(:), ibnd(:) )
       !write(*,'(a7,20i5)') 'ibnd :',ibnd(:)

       !call FillFullMatrix(Fshape, A, FullM)
       !call MblockInverse(Fshape%nsize, FullM )

       !do i=1,npoin
       !   write(*,'(a2,i5,100es10.2)') '..',i,FullM(i, 1:npoin)
       !enddo


       !x(1:npoin) = matmul(FullM(1:npoin, 1:npoin), RHS(1:npoin) )

       !do i=1,npoin
       !   write(*,'(a2,100es10.2)') 'sol',x( 1:npoin)
       !enddo


       !write(*,'(a15)') 'Fshape:'
       !write(*,'(200i5)') Fshape%irwst(:)
       !j = 6
       !do i=1, Fshape%nonzero/j
       !   write(*,'(2i5,a1,14i5)') j*(i-1)+1,i*j,'|', Fshape%idx(j*(i-1)+1:i*j)
       !enddo
       !write(*,'(2i5,a1,14i5)') j*(i-1)+1,i*j,'|', Fshape%idx(j*(i-1)+1:)

       !write(*,'(a5)') 'Ax:'
       !j = 6
       !do i=1, Fshape%nonzero/j
       !   write(*,'(2i5,a1,6es11.3)') j*(i-1)+1,i*j,'|', A(j*(i-1)+1:i*j)
       !enddo
       !write(*,'(2i5,a1,6es11.3)') j*(i-1)+1,i*j,'|', A(j*(i-1)+1:)


       call ReduceProblem(Fshape, A, w(0,:), RHS(:), x, b, ibnd(:))

       !write(31,'(a15)') 'Mshape:'
       !write(31,'(200i5)') Mshape%irwst(:)

       !j = 6
       !do i=1, Mshape%nonzero/j
       !   write(31,'(2i5,a1,6i5)') j*(i-1)+1,i*j,'|', Mshape%idx(j*(i-1)+1:i*j)
       !enddo
       !write(31,'(2i5,a1,6i5)') j*(i-1)+1,i*j,'|', Mshape%idx(j*(i-1)+1:)

       !write(31,'(a5)') 'p Ax:'
       !j = 3
       !do i=1, Mshape%nonzero/j
       !   write(31,'(2i5,a1,3es22.14)') j*(i-1)+1,i*j,'|', A(j*(i-1)+1:i*j)
       !enddo
       !write(31,'(2i5,a1,3es22.14)') j*(i-1)+1,i*j,'|', A(j*(i-1)+1:Mshape%nonzero)

       !do i=1,npoin
       !   write(31,'(a5,i5,4es12.4)') '###',i,w(0,i), x(i), RHS(i), b(i)
       !enddo


       !!call Normalize(Mshape, A, RHS(:) )

       !call Writematrix(Mshape, A)


       !call prodFEM(rez(1:npoin), x(1:npoin), npoin)
       !rez(:) = rez(:) - RHS(:)
       !write(*,'(a7,20es12.4)') 'w_0 :',w(0,1:5)
       !write(*,'(a7,20es12.4)') 'rhs :',RHS(1:5)
       !write(*,'(a7,60es10.2)') 'rez :',rez(1:5)

       ! solution of the corresponding linear algebraic solver
       !iout = 1


       !call Normalize(Mshape, A, b(:) )

!       call gmres(npoin, x(1:Mshape%nsize), b(1:Mshape%nsize), &
!            nit, state%linSolver%tol, prodFEM, bMVnull, nvec,  &
!            0.01*state%linSolver%tol, iout, state%linSolver%iter, state%linSolver%residuum, not_converge)


       !write(31,*)' RHS:'
       !j = 3
       !do i=1, Mshape%nsize/j
       !   write(31,'(2i5,a1,3es22.14)') j*(i-1)+1,i*j,'|', b(j*(i-1)+1:i*j)
       !enddo
       !write(31,'(2i5,a1,3es22.14)') j*(i-1)+1,npoin,'|', b(j*(i-1)+1:npoin)

       !write(31,*)' sol:'
       !j = 3
       !do i=1, Mshape%nsize/j
       !   write(31,'(2i5,a1,3es22.14)') j*(i-1)+1,i*j,'|', x(j*(i-1)+1:i*j)
       !enddo
       !write(31,'(2i5,a1,3es22.14)') j*(i-1)+1,npoin,'|', x(j*(i-1)+1:npoin)

       !call prodFEM(rez(1:npoin), x(1:npoin), npoin)

       !write(31,'(a15)') 'new Ax:'
       !j = 3
       !do i=1, Mshape%nsize/j
       !   write(31,'(2i5,a1,3es22.14)') j*(i-1)+1,i*j,'|', rez(j*(i-1)+1:i*j)
       !enddo
       !write(31,'(2i5,a1,3es22.14)') j*(i-1)+1,npoin,'|', rez(j*(i-1)+1:npoin)

       !write(*,'(a7,120es12.4)') 'sol1 :',x(100:120)

       call solve_umfpack(A(:), x(:), b(:) )

       !call gmres(npoin, w(0,:), RHS, nit, state%linSolver%tol, prodFEM, diagFEM, nvec,  &
       !     0.01*state%linSolver%tol, iout, state%linSolver%iter, state%linSolver%residuum, not_converge)

       !write(*,'(a7,120es12.4)') 'sol :',x(100:120)

       ! possing of he solution from w(*,*) to elem structure
       w(0,:) = x(:)

       call PassSol2Elem(grid, w(0,:) )

       deallocate(Mshape%idx, Mshape%irwst)

       state%time%iter = state%time%iter + 1
       state%time%ttime = state%time%ttime + state%time%tau(1)
       state%timeprn = state%timeprn + state%time%tau(1)

       ! convergence to the steady state solution (if exists)
       call ComputeConvErrors(norm, errL2, norm8, err8)


       ! new version, based on SS residuum
       norm_res = 0.
       if(state%time%iter == 1 .or. state%err(err_0) == 0) state%err(err_0) = max(1E-15, norm_res)
       state%err(SSnew) = norm_res/state%err(err_0)

       state%err(SSL8) = (errL2/norm)**0.5

       ! former stopping criterium
       errL2 = (errL2/norm)**(0.5)/state%time%tau(1)
       state%err(L2) = errL2


       logerr = -999.999
       if(errL2 > 0) logerr = log10(errL2)


       if( state%modelName == 'scalar' .or.state%modelName == '2eqs' ) then
          call ComputeL2H1Error(errL2, errH1, normL2, normH1)
       endif


       !call ComputeEoFC()


       open(ifile, file=convfile, status='OLD', position = 'append')
       !write(ifile,'(i6,16es16.8,i5,i8, 8es16.8, i6)') &
       write(ifile,'(i6,8es16.8,i5,i8, 30es16.8)') &
            state%time%iter, logerr, state%time%ttime, state%time%tau(1), &                 ! 1..4
            state%max_eigenvals*state%time%tau(1), norm_res, state%err(SSL8),&  ! 5..7
            !coeffs(1:5),  state%EcDLM(1:3), &                                ! 8..15
            state%linSolver%tol, state%linSolver%residuum,  state%linSolver%iter, state%linSolver%iter_tot, & ! 16..19
            state%err(SSnew), & !  & !diff , &                                     !20
            1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
            !te-state%start_time, te-ts, tA1e-tA1s, tA2e-tA2s,  &     ! 21..27
            !(tA1e-tA1s)/(te-ts), (tA2e-tA2s)/(te-ts), (tA1e-tA1s+ tA2e-tA2s)/(te-ts), &
            !state%no_refused_steps
       close(ifile)


       ! screen output
       if(mod(iter,20) == 1) then
          if( state%modelName == 'scalar' .or.state%modelName == '2eqs' ) then
             print*,'iter (niter)    tau       time     Ax=b&
                  & res(iter)  ||res||  ||u-u_h||  diff(u_h)'
             print*,'---------------------------------------------------&
                  &----------------------------'
          else
             print*,'iter (niter)    tau       time     Ax=b&
                  & res(iter)  || res ||  l_CFL  %CPU cDLM'
             print*,'---------------------------------------------------&
                  &----------------------------'
          endif
          !call ComputeCF_Pk()
       endif

       yes_conv = '...'
       !if(state%EcDLM(1) <= state%EcD_tol)  yes_conv(1:1) = 'D'
       !if(state%EcDLM(2) <= state%EcL_tol)  yes_conv(2:2) = 'L'
       !if(state%EcDLM(3) <= state%EcM_tol)  yes_conv(3:3) = 'M'

       if(iter == 1 .or. mod(iter,1) == 0 .or. ch == '*')  then
          if( state%modelName == 'scalar' .or.state%modelName == '2eqs' ) then
             write(*,'(i5,a1,i6,a1,es10.2,es11.3,a1,es9.2,a1,i4,a1,es9.2, es12.5,&
                  & es9.2)') iter,'(',state%time%iter,')',state%time%tau(1),&
                  state%time%ttime,ch, state%linSolver%residuum,'(',state%linSolver%iter,')', &
                  !state%err(SSnew), & !norm_res <-> sum1
                  norm_res, &
                  !state%max_eigenvals*state%time%tau(1),&
                  errL2,state%err(SSL8)!  &
                  !' ', (tA2e-tA2s)/(te-ts)!, &
                  !yes_conv
          else
             print*,'NOT IMPLEMENTED'
             !write(*,'(i5,a1,i6,a1,es10.2,es11.3,a1,es9.2,a1,i4,a1, es10.3, es9.2,&
             !     & a1, f4.2, a4)') iter,'(',state%time%iter,')',state%time%tau(1),&
             !     state%time%ttime,ch, state%linSolver%residuum,'(',state%linSolver%iter,')', &
             !     state%err(SSnew), & !norm_res <-> sum1
             !     !norm_res, &
             !     state%max_eigenvals*state%time%tau(1), ' ', (tA2e-tA2s)/(te-ts), &
             !     yes_conv
          endif
       endif


       if(ch == '*') then
          state%isol = state%isol + 1
          call SetFileNames(.false., command_name, tri_name, sol_name, exa_name, err_name, est_name)
          if(state%space%adapt%max_adapt_level >0) call OutputDGFEMtri(tri_name)
          call OutputDGFEMsol(sol_name)
          state%timeprn = state%timeprn - state%time%OutTime

       endif


       call cpu_time(t_end)
       if(t_end - t_sta > 300) then
          !  computation takes more than 5 minuts, we save the achieved results
          call WriteResults('Gsol.bak')

          if(state%time%OutTime <= 0.) then
             call SetFileNames(.false., command_name, tri_name, sol_name, exa_name, err_name, est_name)
             call OutputDGFEMsol(sol_name)
          endif

          it_hours = t_end/3600
          it_min = (t_end - it_hours*3600)/60
          it_sec = t_end - it_hours*3600 - it_min*60
          write(*,'(a40,i3,a1,i2,a1,i2,a11)') &
               'Results saved in file "Gsol.bak" after ', &
               it_hours,':',it_min,':',it_sec,' (hh:mm:ss)'
          t_sta = t_end
       endif


       ! steady state stopping criterion
       if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
          if(state%err(SSL8)/errL2 < state%conv_rez .or. state%time%ttime > state%time%FinTime) goto 100
       else
          print*,'NOT implemented'
          stop
       !   if((state%err(SSnew) <=  state%conv_rez &
       !   !if((norm_res <=  state%conv_rez &
       !        .and.  yes_conv(1:3) =='DLM')  &
       !        .or. state%time%ttime > state%time%FinTime)      goto 100
       endif

    enddo
100 continue

    deallocate(RHS, RHSmass, w) !, residue )
    deallocate(Dphi, phi, A, Mass, Flux)


    if(state%type_IC .eq. 4 .or. state%modelName == 'scalar' .or.state%modelName == '2eqs' ) then
       call ComputeL2H1Error(errL2, errH1, normL2, normH1)
       print*,'Error in L_2 norm = ',errL2
       print*,'Error in H_1 semi-norm = ',errH1

       open(ifile, file='order.dat', status='UNKNOWN', position = 'append')
       write(ifile,'(i6,2es14.6, 3i4, 25es14.6)') &
            grid%nelem, state%space%h, state%time%tau(1),state%space%deg, state%time%deg, &
            grid%curved_deg, errL2, errH1, &
            (errL2**2 + errH1**2/state%model%Re)**0.5, &
            1. , t_end-state%start_time, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
       close(ifile)
    endif

    print*
    print*,'# Iteration process finished'
    print*

    if( state%err(SSnew) > 1E+15) then
       print*,'Too high error, computation stopped'
       stop
    endif



  end subroutine ConformingFEM


  subroutine FillFullMatrix(Fshape, A, FullM)
    type(MatrixShape) :: Fshape  ! full Mshape with trivial lines !!!
    real, dimension(1:Fshape%nonzero), intent(in) :: A
    real, dimension(1:Fshape%nsize, 1:Fshape%nsize ), intent(inout) ::  FullM

    integer :: i,j, k

    FullM(:,:) = 0.
    do i=1,Fshape%nsize
       do j=Fshape%irwst(i), Fshape%irwst(i+1) - 1
          k = Fshape%idx(j)
          FullM(i,k) = A(j)
       enddo
    enddo
  end subroutine FillFullMatrix

  !> removing of zero elements
  subroutine ReduceProblem(Fshape, A, w, RHS, x, b, ibnd)
    !!!type(MatrixShape) :: Mshape  ! real Mshape, global array
    type(MatrixShape) :: Fshape  ! full Mshape with trivial lines !!!
    real, dimension(1:Fshape%nonzero), intent(inout) :: A
    real, dimension(1:Fshape%nsize), intent(in) ::  w, RHS
    real, dimension(1:Fshape%nsize), intent(inout) ::  x, b
    integer, dimension(1:Fshape%nsize), intent(inout) ::  ibnd
    integer :: i,j, ipoc
    real :: val

    !print*,'###',count(A(1:Fshape%nonzero) /= 0.)

    Mshape%nsize = Fshape%nsize
    !Mshape%nonzero = Fshape%nonzero
    Mshape%nonzero = count(A(1:Fshape%nonzero) /= 0.)

    allocate(Mshape%irwst(1:Mshape%nsize + 1) )
    allocate(Mshape%idx(1:Mshape%nonzero) )

    ipoc = 0
    Mshape%irwst(1) =  1
    do i=1,Fshape%nsize

       do j=Fshape%irwst(i), Fshape%irwst(i+1) - 1
          if( A(j) /= 0) then
             ipoc = ipoc + 1
             A(ipoc) = A(j)

             Mshape%idx(ipoc) = Fshape%idx(j)
          endif
       enddo
       Mshape%irwst(i+1) = ipoc + 1
    enddo

    x(:) = w(:)
    b(:) = RHS(:)


  end subroutine ReduceProblem

  !> recomputation of the solution from vector to elem
  subroutine PassSol2Elem(grid0, w )
    type(mesh), intent(in) :: grid0
    real, dimension(1:grid0%npoin), intent(inout) :: w
    class(element), pointer :: elem

    real, dimension(1:ndim, 1:3) :: wi, wj
    integer :: i,j,k, dof

    do i=1,grid0%nelem
       elem => grid0%elem(i)
       dof = elem%dof

       do j=1,elem%flen
          k = elem%face(idx,j)
          wi(1:ndim, j) = w(k)
          !!write(99,*) grid0%x(k, 1:nbDim), w(k), dof
       enddo

       call Lagr2BasisDiff(elem, 1, wi(1:ndim, 1:3), elem%Qnum, dof, wj(1:ndim, 1:dof) )

       do k=1,ndim
          elem%w(0, (k-1)*dof + 1: k*dof) = wj(k, 1:dof)
       enddo

    enddo

  end subroutine PassSol2Elem


  subroutine Normalize(Mshape, A, b )
    type(MatrixShape), intent(in) :: Mshape
    real, dimension(1:Mshape%nonzero), intent(inout) :: A
    real, dimension(1:Mshape%nsize), intent(inout) ::  b
    integer :: i,j
    real :: val

    ! moving the given term from the left to the right
    do i=1,Mshape%nsize
       val = A(Mshape%irwst(i))
       if(val /= 0.) then
          b(i) = b(i) /val
          do j= Mshape%irwst(i), Mshape%irwst(i+1)-1
             A(j) = A(j) / val
          enddo
       endif
    enddo

  end subroutine Normalize


  !> change of matrix and vectors due to the Dirichlet boundary condition
  subroutine RemoveTrivialElements(Mshape, A, x, b, ibnd )
    type(MatrixShape), intent(in) :: Mshape
    real, dimension(1:Mshape%nonzero), intent(inout) :: A
    real, dimension(1:Mshape%nsize), intent(inout) :: x, b
    integer, dimension(1:Mshape%nsize), intent(in) :: ibnd
    integer :: i,j,k

    ! moving the given term from the left to the right
    do i=1,Mshape%nsize
       do j= Mshape%irwst(i), Mshape%irwst(i+1)-1
          k = Mshape%idx(j)

          if(ibnd(k) == 0) then
             b(i) = b(i) - A(j)*x(k)
             A(j) = 0.
          endif
       enddo
    enddo

    ! trivial lines for the known value
    do i=1,Mshape%nsize
       if(ibnd(i) == 0) then
          A(Mshape%irwst(i)) = 1.    ! first element in row is the diagonal one

          A(Mshape%irwst(i)+1: Mshape%irwst(i+1)-1) = 0.

          b(i) = x(i)
       endif
    enddo

  end subroutine RemoveTrivialElements


  !> setting of the Dirichlet BC
  subroutine InitDirichletBC(grid0, ibnd)
    type(mesh), intent(in) :: grid0
    integer, dimension(1:grid0%npoin), intent(inout) :: ibnd
    integer :: ib

    ibnd(:) = 1
    do ib=1,grid0%nbelm
       ibnd(grid0%b_edge(ib)%lbn(1:nbDim)) = 0
    enddo

  end subroutine InitDirichletBC

  !> setting of the Dirichlet BC
  subroutine SetDirichletBC(grid0,  w )
    type(mesh), intent(in) :: grid0
    real, dimension(1:grid0%npoin), intent(inout) :: w
    real, dimension(1:nbDim) :: xi
    real, dimension(1:ndim) :: fx

    integer :: nbelm, npoin, ib,j,k, ibc

    nbelm = grid0%nbelm

    do ib=1,nbelm
       !if(grid0%b_edge(ib)%BC == 1) then ! Dirichlet boundary condition

       do j = 1,2
          k = grid0%b_edge(ib)%lbn(j)
          xi(1:nbDim) = grid0%x(k, 1:nbDim)
          call Set_Model_Data(xi(1:nbDim), state%time%ctime, fx(1:ndim), 1)
          w(k) = fx(1)
          !write(*,*) xi(:), w(k)
       enddo
    enddo

  end subroutine SetDirichletBC

  !> setting of the vector from elem%w
  subroutine InitSolutionVector(grid0, w )
    type(mesh), intent(in) :: grid0
    real, dimension(1:grid0%npoin), intent(inout) :: w
    class(element), pointer :: elem
    real, dimension(:,:), pointer :: phi
    real, dimension(:), allocatable :: q, qi
    real, dimension(:), allocatable :: weights
    integer :: npoin, i,j,k,l, Qnum, Qdof,  dof

    npoin = grid0%npoin

    Qdof = 3
    Qnum = 1  ! piecewise linear
    Qdof = state%space%L_rule(Qnum)%Qdof
    phi => state%space%L_rule(Qnum)%phi(1:state%space%max_dof, 1:Qdof)

    allocate(q(1:state%space%max_dof), qi(1:Qdof) )

    allocate(weights(1:npoin))
    weights(:) = 0.
    w(:) = 0.

    do i=1,grid0%nelem
       elem => grid0%elem(i)

       dof = elem%dof
       q(1:dof) = elem%w(0,1:dof)   ! only first component

       ! evaluation of w in the Langrangian nodes
       do l=1, Qdof
          qi(l) = dot_product( q(1:dof), phi(1:dof,l) )
       enddo

       do j=1,elem%flen
          k = elem%face(idx,j)
          w(k) = w(k) + qi(j)
          weights(k)  = weights(k) +1
       enddo
    enddo

    w(1:npoin) = w(1:npoin)/weights(1:npoin)

    do i=1,npoin
       write(66,*) grid0%x(i,1:nbDim), w(i)
    enddo

    deallocate(weights, q, qi)

  end subroutine InitSolutionVector


  !> evaluation of linear convection term
  subroutine ComputeConvectionFEM(grid0, Mshape, Dphi, phi, A)
    type(MatrixShape), intent(inout) :: Mshape
    type(mesh), intent(in) :: grid0
    real, dimension(1:grid0%nelem, 1:3, 1:nbDim), intent(inout) :: Dphi
    real, dimension(1:3, 1:3), intent(inout) :: phi
    real, dimension(1:Mshape%nonzero), intent(inout) :: A
    real, dimension(1:nbDim) :: xi, b
    class(element), pointer :: elem
    real, dimension(:,:), pointer :: x
    integer :: i,j,j1, k, k1, nij, l, l1
    real :: val

    x => grid0%x(1:grid0%npoin, 1:nbDim)

    do i=1,grid0%nelem
       elem => grid0%elem(i)
       do j=1,elem%flen                 ! index of row
          k = elem%face(idx, j)
          do j1=1,elem%flen             ! index of column
             k1 = elem%face(idx, j1)

             nij = Mshape%nij(i,j,j1)

             val = 0.
             do l=1, 3  ! over integ nodes, midle of edges
                l1 = mod(l,3) + 1
                xi(1:nbDim) = (x(elem%face(idx,l), 1:nbDim) + x(elem%face(idx,l1), 1:nbDim))/2
                b(1) = Scal_Lin(1, xi(1), xi(2), state%time%ctime)
                b(2) = Scal_Lin(2, xi(1), xi(2), state%time%ctime)

                val = val + (Dphi(i, j1, 1)*b(1) + Dphi(i, j1, 2)*b(2))*phi(j,l)

             enddo

             !write(21, *) xi(:), val * elem%area/3

             !if(k==5 .and. k1 == 5) write(*,'(a5,6i5,2es12.4)')'conv', i,j,j1,k,k1,nij, &
             !     val, A(nij) + val * elem%area/3

             A(nij) = A(nij) + val * elem%area/3
          enddo
       enddo
    enddo

  end subroutine ComputeConvectionFEM

  !> evaluation of the source term
  subroutine ComputeSourceTerms(grid0, phi, f )
    type(mesh), intent(in) :: grid0
    real, dimension(1:3, 1:3), intent(inout) :: phi
    real, dimension(1:grid0%npoin), intent(inout) :: f
    real, dimension(1:nbDim) :: xi
    class(element), pointer :: elem
    real, dimension(:,:), pointer :: x
    integer :: i,j,j1, k, k1, nij, l, l1
    real, dimension(1:ndim) :: fx
    real :: val

    x => grid0%x(1:grid0%npoin, 1:nbDim)

    do i=1,grid0%nelem
       elem => grid0%elem(i)
       do j=1,elem%flen                 ! index of row
          k = elem%face(idx, j)

          val = 0.
          do l=1, 3  ! over integ nodes, midle of edges
             l1 = mod(l,3) + 1
             xi(1:nbDim) = (x(elem%face(idx,l), 1:nbDim) + x(elem%face(idx,l1), 1:nbDim))/2

             call Set_Model_Data(xi(1:nbDim), state%time%ctime, fx(1:ndim), 2)

             val = val + fx(1)*phi(j,l)
             !if(k == 1) write(*,'(a5,12es12.4)') 'fx :',xi(1:nbDim), fx(1),phi(j,l),val
          enddo

          f(k) = f(k) + val * elem%area/3
       enddo

    enddo

    !do i=1,grid0%npoin
    !   print*,'rhs:',i, f(i)
    !enddo
  end subroutine ComputeSourceTerms

  !> evaluation of the RHS mass term
  subroutine ComputeRHSmass(grid0, phi, w, f )
    type(mesh), intent(in) :: grid0
    real, dimension(1:3, 1:3), intent(inout) :: phi
    real, dimension(1:grid0%npoin), intent(in) :: w
    real, dimension(1:grid0%npoin), intent(inout) :: f
    real, dimension(1:3) :: wi
    class(element), pointer :: elem
    integer :: i,j,j1, k, k1, l, l1
    real :: val

    do i=1,grid0%nelem
       elem => grid0%elem(i)
       do j=1,elem%flen                 ! index of row
          j1 = mod(j, 3) + 1
          k = elem%face(idx, j)
          k1 = elem%face(idx, j1)

          wi(j) = (w(k) + w(k1))/2.  ! values in middle of edges

          !if(k == 7)
       enddo

       do j=1,elem%flen                 ! index of row
          k = elem%face(idx, j)

          val = 0.
          do l=1, 3  ! over integ nodes, midle of edges
             val = val + wi(l)*phi(j,l)

          enddo

          !if(k == 7 )  write(*,'(a5,12es12.4)') 'fx :',wi(1:3),  val

          f(k) = f(k) + val * elem%area/3
       enddo

    enddo

    !do i=1,grid0%npoin
    !   print*,'rhsMass:',i, f(i)
    !enddo
  end subroutine ComputeRHSmass

  !> evaluation of linear reactive term
  subroutine ComputeReactionFEM(grid0, Mshape, phi, A)
    type(MatrixShape), intent(inout) :: Mshape
    type(mesh), intent(in) :: grid0
    real, dimension(1:3, 1:3), intent(inout) :: phi
    real, dimension(1:Mshape%nonzero), intent(inout) :: A
    real, dimension(1:nbDim) :: xi
    class(element), pointer :: elem
    real, dimension(:,:), pointer :: x
    integer :: i,j,j1, k, k1, nij, l, l1
    real :: val, c

    x => grid0%x(1:grid0%npoin, 1:nbDim)

    do i=1,grid0%nelem
       elem => grid0%elem(i)
       do j=1,elem%flen                 ! index of row
          k = elem%face(idx, j)
          do j1=1,elem%flen             ! index of column
             k1 = elem%face(idx, j1)

             nij = Mshape%nij(i,j,j1)

             val = 0.
             do l=1, 3  ! over integ nodes, midle of edges
                l1 = mod(l,3) + 1
                xi(1:nbDim) = (x(elem%face(idx,l), 1:nbDim) + x(elem%face(idx,l1), 1:nbDim))/2
                c = Scal_Lin(0, xi(1), xi(2), state%time%ctime)

                val = val + c*phi(j1, l)*phi(j,l)

             enddo

             !if(k==5 .and. k1 == 5) write(*,'(a5,6i5,2es12.4)')'reac', i,j,j1,k,k1,nij, &
             !     val, A(nij) + val * elem%area/3

             A(nij) = A(nij) + val * elem%area/3
          enddo
       enddo

    enddo
  end subroutine ComputeReactionFEM


  !> evaluation of the stiff matrix - linear diffusion term
  subroutine ComputeStiffFEM(grid0, Mshape, Dphi, A)
    type(MatrixShape), intent(inout) :: Mshape
    type(mesh), intent(in) :: grid0
    real, dimension(1:grid0%nelem, 1:3, 1:nbDim), intent(inout) :: Dphi
    real, dimension(1:Mshape%nonzero), intent(inout) :: A
    class(element), pointer :: elem
    integer :: i,j,j1, k, k1, nij

    do i=1,grid0%nelem
       elem => grid0%elem(i)
       do j=1,elem%flen
          k = elem%face(idx, j)
          do j1=1,elem%flen
             k1 = elem%face(idx, j1)

             nij = Mshape%nij(i,j,j1)
             A(nij) = A(nij) + &
                  (Dphi(i, j, 1)*Dphi(i, j1, 1) + Dphi(i, j, 2)*Dphi(i, j1, 2)) * elem%area

          enddo
       enddo

    enddo
  end subroutine ComputeStiffFEM

  !> evaluation of the mass matrix - linear diffusion term
  subroutine ComputeMassFEM(grid0, Mshape, phi, A)
    type(MatrixShape), intent(inout) :: Mshape
    type(mesh), intent(in) :: grid0
    real, dimension(1:3, 1:3), intent(inout) :: phi
    real, dimension(1:Mshape%nonzero), intent(inout) :: A
    class(element), pointer :: elem
    integer :: i,j,j1, k, k1, l, nij
    real :: val

    do i=1,grid0%nelem
       elem => grid0%elem(i)
       do j=1,elem%flen
          k = elem%face(idx, j)
          do j1=1,elem%flen
             k1 = elem%face(idx, j1)

             val = 0.
             do l=1,3
                val = val + phi(j, l) * phi(j1,l)
             enddo

             nij = Mshape%nij(i,j,j1)

             A(nij) = A(nij) + val * elem%area /3

          enddo
       enddo

    enddo
  end subroutine ComputeMassFEM


  !> evaluation of the shape functions in middle of edges
  subroutine ComputePhi( phi)
    real, dimension(1:3, 1:3), intent(inout) :: phi

    phi(1,1:3) = (/ 0.5, 0. , 0.5  /)
    phi(2,1:3) = (/ 0.5, 0.5, 0.   /)
    phi(3,1:3) = (/ 0. , 0.5, 0.5  /)

  end subroutine ComputePhi

  !> evaluation of the gradient of the shape functions
  subroutine ComputeDphi(grid0, Dphi)
    type(mesh), intent(in) :: grid0
    real, dimension(1:grid0%nelem, 1:3, 1:nbDim), intent(inout) :: Dphi
    class(element), pointer :: elem
    real, dimension(:,:), pointer :: x

    integer :: npoin, nelem
    integer :: i

    nelem = grid0%nelem
    npoin = grid0%npoin

    x=> grid0%x(1:npoin, 1:nbDim)

    do i=1,nelem
       elem=> grid0%elem(i)

       call Gradient(x(elem%face(idx,1:3), 1:nbDim), Dphi(i, 1:3, 1:nbDim) )

       !write(*,'(a5,6es12.4)')'Dphix', Dphi(i,1:3,1)
       !write(*,'(a5,6es12.4)')'Dphiy',Dphi(i,1:3,2)
       !print*
    enddo

  end subroutine ComputeDphi

  !> Computing of the gradient
  subroutine Gradient(x, Dphi)
    real, dimension(1:3,1:nbDim), intent(in) :: x  ! vertex coordinates
    real, dimension(1:3,1:nbDim), intent(out) :: Dphi   ! derivative of test functions
    real :: D, Dx, Dy
    integer :: j, j1, j2

    D = x(1,1)*(x(2,2) - x(3,2) ) + x(2,1)*(x(3,2) - x(1,2) ) + x(3,1)*(x(1,2) - x(2,2) )

    do j=1,3
       j1 = mod(j, 3) + 1
       j2 = mod(j1, 3) + 1

       Dphi(j,1) = (x(j1,2) - x(j2, 2) ) / D
       Dphi(j,2) = (x(j2,1) - x(j1, 1) ) / D
    enddo

  end subroutine Gradient



  !> setting of the matrix shape for P_1 conforming FEM
  subroutine SetMshape(Mshape, grid0)
    type(MatrixShape), intent(inout) :: Mshape
    type(mesh), intent(in) :: grid0
    class(element), pointer :: elem
    integer, parameter :: max_cyc = 20
    integer, dimension(:), allocatable :: len_cyc
    integer, dimension(:,:), allocatable :: cyc

    integer :: npoin, nelem, nonzero
    integer :: i, j, k, j1, k1, l

    nelem = grid0%nelem
    npoin = grid0%npoin


    ! to each vertex we create a list of nodes sharing this node
    allocate(cyc(1:npoin, 1:max_cyc))
    allocate(len_cyc(npoin))

    len_cyc(1:npoin) = 0

    do  i=1, nelem
       elem => grid0%elem(i)
       if(elem%HGnode .or. elem%flen /= 3) then
          print*,'mesh contains HG node, it is not supported'
          print*,'element is not triangle'
          stop
       endif

       do j=1, elem%flen
          k = elem%face(idx,j)

          do j1=1, elem%flen
             k1 = elem%face(idx,j1)

             if(k == k1) goto 100
             do l = 1, len_cyc(k1)
                if( cyc(k1, l) ==  k) goto 100
             enddo

             len_cyc(k1) = len_cyc(k1) + 1
             cyc(k1,len_cyc(k1)) = k

             !write(*,'(a5,5i5,a1,20i5)') &
             !     '...:',i,j,k, j1, k1, '|',k1, len_cyc(k1), cyc(k1, 1: len_cyc(k1) )

100          continue
          enddo
       enddo
    enddo

    nonzero = sum(len_cyc(1:npoin)) + npoin
    Mshape%nsize = npoin
    Mshape%nonzero = nonzero
    allocate(Mshape%idx(1:nonzero) )
    allocate(Mshape%irwst(1:npoin+1) )
    allocate(Mshape%nij(1:nelem, 1:3, 1:3) )


    nonzero = 0
    Mshape%irwst(1) = 1
    do i=1,npoin
       !write(*,'(a5,20i5)') 'cyc:',i, len_cyc(i), cyc(i, 1: len_cyc(i) )
       Mshape%idx(nonzero + 1) = i
       Mshape%idx(nonzero + 2 : nonzero + len_cyc(i) + 1) = cyc(i, 1: len_cyc(i) )

       !write(*,'(a2,50i3)') 'Q',Mshape%idx(:)
       nonzero = nonzero + len_cyc(i) + 1
       Mshape%irwst(i+1) = nonzero + 1

    enddo

    if(nonzero /= Mshape%nonzero) then
       print*,' Mishmatch size in Set Mshape'
       stop
    endif

    do i=1,nelem
       elem => grid0%elem(i)
       do j=1, elem%flen
          k = elem%face(idx,j)

          do j1=1, elem%flen
             k1 = elem%face(idx,j1)

             if( k == k1) then
                Mshape%nij(i, j, j1) = Mshape%irwst(k)   ! diagonal element
             else
                do l=1, len_cyc(k)
                  ! write(*,'(a5,5i5,a2,5i5)') '####',i,j,j1,k,k1,'|',l, cyc(k,l),Mshape%irwst(k) + l
                   if(cyc(k,l) == k1)  then
                      Mshape%nij(i, j, j1) = Mshape%irwst(k) + l
                      goto 200
                   endif
                enddo
200             continue

             endif
          enddo
       enddo
       !write(*,'(a5,3(3i5,a1))') 'nij',Mshape%nij(i,1, 1:3),'|', Mshape%nij(i,2, 1:3), &
       !     '|',Mshape%nij(i,3, 1:3),'|'
       !if(i==3) stop
    enddo

    deallocate(cyc, len_cyc)

  end subroutine SetMshape


  !> write the matrix shape
  subroutine WriteMshape(Mshape)
    type(MatrixShape), intent(inout) :: Mshape
    integer :: ifile = 15
    integer :: i, j

    open(ifile, file='Mshape', status='unknown')
    do i=1, Mshape%nsize
       do j=Mshape%irwst(i), Mshape%irwst(i+1)-1
          write(15,*) i,  Mshape%idx(j) , j
       enddo
    enddo
    close(ifile)
  end subroutine WriteMshape


  !> write the matrix
  subroutine WriteMatrix(Mshape, A)
    type(MatrixShape), intent(inout) :: Mshape
    real, dimension(1:Mshape%nonzero), intent(in) ::A
    real, dimension(:), allocatable :: row
    integer :: i,j

    allocate(row(1:Mshape%nsize) )

    print*,'-------------matrix------------------'
    do i=1, Mshape%nsize
       row(:) = 0.

       do j=Mshape%irwst(i), Mshape%irwst(i+1)-1
          row(Mshape%idx(j) ) = A(j)
       enddo

       write(*,'(i5,a1,30es12.4)') i,':',row(:)
    enddo
    print*,'-------------matrix------------------'

    deallocate(row)

  end subroutine WriteMatrix


end module conforming_FEM
