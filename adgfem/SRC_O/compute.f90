!> Navier-Stokes equations
module compute_oper
  use main_data
  use io_sub
  use solve_problem
  use marking
  use estimates
  use errorDual
  use errorFlux
  use hp_adaptation
  !!use pMultiGrid !MAKE_MG
  use helmholtz_estim
  use dwr_mod
  use dual_estim
!  use alg_estim
  use st_interpol
  use anisotropic
  !use AMA_estims
  use ama_hp_interpol
  use error_subs
  use estimates
  use mesh_oper
  use newton_mod
  use lin_solvers
  use loc_rav_tho_ned
  use solution_mod
  use stdgm_mod

  implicit none

  public:: finishComputationProcess
  public:: ResetComputation
  public:: PerformOneTimeStep
  public:: PerformOneSTDGMstep
  !public:: PerformOneSTDGMstep_OLD
  public:: PerformOneSTDGMstep_Andrerson
  public:: PassToNewTimeStep
  public:: WriteOutputFiles
  public:: WriteOutputError
  public:: WriteOutputScreen
  public:: WriteOutputsAfterSolveProblem

contains

  !> resetting of the computation, i.e., start from the beginning (\f$ t=0 \f$),
  !> recomputation of the initial condition to the actual grid
  subroutine ResetComputation( )
    print*,'#---------------------------------------------------'
    print*,'# Reseting of the computation, original IC'
    print*,'#---------------------------------------------------'


    if(state%type_IC > 0) then
       state%time%ttime = 0.
       state%time%ctime = 0.
       state%time%iter = 0
       state%timeprn = 0.

       call SetElementsIC()
    endif

  end subroutine ResetComputation


  !> perform the one time step on the given grid
  !>
  !> evaluate the flux matrix/vector, solve by explicit, semi-implicit or implicit
  !> time discretization, for the implicit the Newton-like method is used
  subroutine PerformOneTimeStep( time_prepare, time_solve, time_estim )
    real, intent(inout) ::  time_prepare, time_solve, time_estim
    !type(NonlinearSol_t), pointer :: Newton
    class(element), pointer :: elem
    real, dimension(:,:), allocatable :: cDLM
    real:: t1, t2, rel_res, tt, res_max_val
    !real :: mu, mu1
    !real :: normL2, normH1
    integer :: iter, newton_count, elemDof
    integer :: i, j, k, l, i1, i2, kk
    logical :: vector_update, loc_implicitly, update
    integer :: imp
    character(len=7) :: Newtonx, MatrixA, RHSideb, MatrixB, MatrixC
    real,dimension(:),allocatable :: y,z !MAKE_MG
    !integer:: deg,dof,nelem,nsize,outsize
    integer :: vypis
    !type(MGSetting) :: stateMG
    logical :: deg_plus

    deg_plus = .false.

    associate ( Newton => state%nlSolver )

    ! seting of coefficients \alpha, \beta for ABDF
    i = state%time%deg_actual
    call state%time%SetBDF( i, state%time%tau(1:i+1) )

    !print*,'Compute "Mass RHS", vector M w^{k} and extrapolate from the old levels'
    call SetMassVector( )

    time_prepare = 0.
    time_solve  = 0.
    time_estim = 0.

    !!write(*,'(a6, 3es12.4)') 'times:',time_prepare, time_solve , time_estim

    !Newton => state%nlSolver
    allocate(cDLM(0:Newton%max_iter, 1:5))
    Newton%norm_res = 1.

    !if(ndim == 4) call DirectComputeDrag_Lift(cDLM(0, 1:5 ) )
    !open(44, file='cDLM', status='UNKNOWN', position = 'append')
    !write(44,'(x)')
    !write(44,'(2i5,15es11.3)') state%time%iter, Newton%iter, &
    !     state%time%iter + 1.*0., cDLM(Newton%iter, 1:3)
    !close(44)

    ! minimal update of the flux matrix after imp time steps
    imp = Newton%min_update
    if(state%modelName == 'scalar' .or. state%modelName == '2eqs') imp = max(1, imp)

    if( state%modelName == 'porous' )  imp = max(1, imp)
    !if( state%modelName == 'porous' ) &
    !     print*,"if( state%modelName == 'porous' )  imp = max(1, imp)"

    ! pseudo-time stepping
    if(state%time%time_method == 'P')  imp = 1000

    ! for DUA error estimates
    grid%elem(:)%CK = 0.
    grid%elem(:)%Cb = 0.
    grid%elem(:)%CKo = 0.
    grid%elem(:)%Cbo = 0.

    ! Has the big vector the right size? If not (e.g., after mesh adapation) then resize
    if(size(Newton%b(:) ) /= state%nsize  )then
       deallocate( Newton%b, Newton%b1, Newton%x, Newton%rr )

       allocate(Newton%b(1:state%nsize), Newton%b1(1:state%nsize), &
            Newton%x(1:state%nsize), Newton%rr(1:state%nsize) )
    endif

    if(state%time%time_method == 'I' .and. imp > 0) then
       if( mod(state%time%iter_loc,imp) == 1 .or. imp == 1)  state%nlSolver%implicitly = .true.

    elseif(state%time%time_method == 'I' .and. imp == 0) then
       state%nlSolver%implicitly = .true.

    endif


    if(state%time%time_method == 'I') then
       open(54, file='newton-update', status='unknown', position='append')
       open(55, file='newton', status='unknown', position='append')
       !open(56, file='implicit', status='unknown', position='append')
       write(55,*) '  '
       write(55,*) '## it  i  l    ritr     ritr2   ', &
            '   |F(x)|/N  |F(x*)|/N   |F(x)|   |F(x*)|      rel_res ', &
            '     gmr_tol    lambda     theta       etaS    algeb   LAiter IM Pr'
    endif

    if( state%time%cn ) call PrepareCrankNicolson( )

    state%time%keep_tau = .true.

    Newton%newton_count = 0
    Newton%iter  = 0
    do iter = 1, Newton%max_iter
       call cpu_time(t1)

       if(imp == 0)  state%nlSolver%implicitly= .true.

       loc_implicitly = state%nlSolver%implicitly
!       state%linSolver%precond_update = .false.
       update = .false.

       if( state%nlSolver%implicitly ) then
          !print*,'Compute matrix C(w), vector q(w)'
          !write(600,*)'# Compute matrix C(w), vector q(w)'
          update = .true.

          !write(*,*) "Newtonova iterace 1, iter, w", iter, grid%elem(1)%w(0,1:ndim)  ! wet steam testing
          !write(*,*) "Newtonova iterace 1", iter

          !call cpu_time(tt)
          !write(*,'(a30, 2es12.4)') 'CT_starts 1 :', tt - state%start_time

          call ComputeTerms(deg_plus)

          !call cpu_time(tt)
          !write(*,'(a30, 2es12.4)') 'CT_ends   1 :', tt - state%start_time
          !print*

          !write(*,*) "@@@@@@@@"   ! wet steam testing
          !call WriteMatrixA(eta)
          !write(*,*) "$$$$$$$$"

          Newton%newton_count = 0
          Newton%updates = Newton%updates + 1
          state%linSolver%precond_update = .true.

          if(state%time%time_method == 'I' .or. state%time%time_method == 'S') then
             if(state%modelName == 'NSe') then
                !print*, 'AA2'
                !write(*,*) "GGGGGGGG, iter", iter
                !call WriteMatrixA(eta)   !  wet steam testing
                !write(*,*) "HHHHHHHH"

                call SetF_q_Cw_fast( )

                !write(*,*) 'TTTTTTTTT, iter', iter, eta
                !call WriteMatrixA(eta)   ! neni NaN wet steam testing
                !write(*,*) "XXXXXXXX"
             else

                !print*, 'AA1'
                !print*
                !print*,'Compute matrix C(w), vector q(w) - (2)'
                !write(600,*) '# Compute matrix C(w), vector q(w) - (2)'
                state%nlSolver%implicitly = .false.

                !write(*,*) "Newtonova iterace 2", iter

                !call cpu_time(tt)
                !write(*,'(a30, 2es12.4)') 'CT_starts 2 :', tt - state%start_time

                call ComputeTerms(deg_plus )

                !call cpu_time(tt)
                !write(*,'(a30, 2es12.4)') 'CT_ends   2 :', tt - state%start_time
                !print*


                state%nlSolver%implicitly = .true.
             endif
          endif

          vector_update = .true.
          if(state%time%time_method == 'I') state%nlSolver%implicitly = .false.

       else   !!!if( .not. state%nlSolver%implicitly .and. Newton%iter == 1) then

          !write(*,*) "if ( .not. state%nlSolver%implicitly .and. Newton%iter == 1 ))", .not. state%nlSolver%implicitly, Newton%iter   ! wet steam testing

          ! for iter > 1, array b(:) = F(x^k) was already computed
          if(iter == 1) then
             !print*,'Compute  ONLY vector q(w)'
             !write(600,*) '# Compute  ONLY vector q(w)'

             !write(*,*) "Newtonova iterace 3", iter

             !call cpu_time(tt)
             !write(*,'(a30, 2es12.4)') 'CT_starts 3 :', tt - state%start_time

             call ComputeTerms(deg_plus )

             !call cpu_time(tt)
             !write(*,'(a30, 2es12.4)') 'CT_ends   3 :', tt - state%start_time
             !print*

             vector_update = .true.
          else
             !write(*,*) 'YYYYYYYYYY, iter', iter  ! wet steam testig
             !call WriteMatrixA(eta)   !  wet steam testig
             !write(*,*) "QQQQQQQQQQ, iter", iter   ! wet steam testig
             vector_update = .false.
             Newton%b(1:state%nsize) = Newton%b1(1:state%nsize)
             Newton%res = Newton%res0
          endif


       endif ! end of if(state%nlSolver%implicitly)


       !eta = state%time%alpha(0)/state%time%tau(1)
       eta = 1./state%time%tau(1)



       if(state%time%time_method == 'E') then

          ! explicit discretization, direct solution of the system via
          ! multiplication of the Mass Matrix inversion
          call SolveBlockDiagonalSystem(eta, Newton%norm_res)

          print*,'VDVDVDVD',Newton%norm_res

       else if(state%time%time_method == 'P') then

          ! pseudo-time stepping  NOT WORKING

           ! if(iter == 1) then
           !    do i=1,grid%nelem
           !       elem => grid%elem(i)
           !       do j=1,elem%flen
           !          l = elem%face(neigh, j)
           !          if(l > 0) elem%max_eigenvals = max(elem%max_eigenvals, grid%elem(l)%max_eigenvals )
           !       enddo

           !    enddo
           ! endif


          if(vector_update) then
             call FillVector(Newton%b(1:state%nsize), eta ) ! value of eta is not used, only its sign
             !write(*,'(a6,i5, 40es12.4)') 'flux:',0,Newton%b(1:state%nsize)
             !Newton%res  = VectorNorm(Newton%b)
             Newton%res  = VectorPrecondNorm(Newton%b)
          endif


          ! size of the pseudo-time step
          Newton%lambda = 0.02 / state%max_eigenvals

          !print*,'pseud tau =',  Newton%lambda
          !pseudo-time update
          Newton%res0  = 0
          k = 0
          do i=1,grid%nelem
             elem => grid%elem(i)

             !Newton%lambda = 0.005 / elem%max_eigenvals
             Newton%lambda = 0.01 / (elem%max_eigenvals *  state%max_eigenvals)**0.5

             j = elem%dof

             do l=1,ndim
                i1 = (l-1)*j + 1
                i2 = l*j


                elem%w(0,i1:i2) = elem%w(0,i1:i2)  &
                     !+ Newton%lambda *  Newton%b(k+1:k+j)
                     + Newton%lambda * matmul(elem%MassInv%Mb(1:j,1:j),  Newton%b(k+1:k+j) )

                k = k+j

             enddo
          enddo

          !stop

          Newton%lambda_old = Newton%lambda

          !print*,'Compute ONLY vector q(w)', state%nlSolver%implicitly
          !write(600,*) '# Compute ONLY vector q(w)', state%nlSolver%implicitly

          !write(*,*) "Newtonova iterace 4", iter
          !call cpu_time(tt)
          !write(*,'(a30, 2es12.4)') 'CT_starts 4 :', tt - state%start_time

          call ComputeTerms(deg_plus )

          !call cpu_time(tt)
          !write(*,'(a30, 2es12.4)') 'CT_ends   4 :', tt - state%start_time
          !print*


          call FillVector(Newton%b1(1:state%nsize),eta ) !include time derivative term

          !!!  VD  VD
          Newton%res0 = VectorPrecondNorm( Newton%b1 )
          !!!Newton%res0  = (Newton%res0 / state%nsize)**0.5

          Newton%theta  =  Newton%res0  / Newton%res

          !write(*,'(a6,i5,20es18.10)') 'Pseud:',iter,  Newton%res,  Newton%res0,  Newton%theta, Newton%lambda
          !print*

       else  ! fully implicit or semi-implicit methods

          if(vector_update) then ! vector_update == .TRUE. for wet steam case
             call FillVector(Newton%b(1:state%nsize), eta ) ! value of eta is not used, only its sign
             !Newton%res  = VectorNorm(Newton%b)
             Newton%res  = VectorPrecondNorm(Newton%b)
          endif

          !do i=1,grid%nelem
          !   elem => grid%elem(i)
          !   write(*,'(a8,i5,200es12.4)') 'RHS:', elem%i, elem%vec(rhs,1:9)
          !enddo
          !stop 'rde3se4de4wd3d43e in compute.f90'

          !write(*,*) 'JJJJJJJJJJ, iter', iter, eta, state%time%tau(:)
          !call WriteMatrixA(eta)   ! objevi se NaN wet steam testing
          !write(*,*) "SSSSSSSSS"
          !write(*,'(a6,500es12.4)') 'RHS',Newton%b(1:state%nsize)
          !stop


          if(iter == 1) Newton%res_ref = Newton%res

          ! first exit based on the absolute value of the residuum vector || F(w) ||
          if(state%time%time_method == 'I' .and. iter > 1  &
               .and. Newton%res/state%nsize < Newton%tol) goto 200

          Newton%iter  = Newton%iter + 1
          Newton%Aiter = Newton%Aiter  + 1

          Newton%lambda_old = 0.
          Newton%x(1:state%nsize) = 0.

          if (.not. state%linSolver%tol_fixed) then

             if(iter == 1) then
                !state%linSolver%tol = 0.5
                state%linSolver%tol = 0.25
             else
                !state%linSolver%tol = min(0.5, 0.5*(Newton%res / Newton%res1)**1.5)
                !state%linSolver%tol = min(0.5, 0.25*(Newton%res / Newton%res1)**2.0)
                !state%linSolver%tol = 0.25**Newton%iter
                state%linSolver%tol = min(0.25, abs(Newton%res - state%linSolver%residuum) / Newton%res1)

                !write(*,'(a4,2es12.4,a1,4es12.4)') '@@@@', &
                !     state%linSolver%tol,  min(0.5, 0.25*(Newton%res / Newton%res1)**2.0), &
                !     '|', abs(Newton%res - state%linSolver%residuum) / Newton%res1, &
                !     Newton%res, state%linSolver%residuum ,  Newton%res1
                state%linSolver%tol = max(state%linSolver%tol, 1E-6)
             endif
          endif

          !print*,'Solution of linear algebraic system'
          eta = state%time%alpha(0)/state%time%tau(1)
          if(state%time%cn ) eta = eta*2  ! Crank-Nicolson

          call cpu_time(t2)
          time_prepare = time_prepare + t2 - t1

          !write(*,'(a6, 3es12.4)') 'timesA',time_prepare, time_solve , time_estim

!          if ( state%space%adapt%adapt_method == 'ALG' .or. &
!                   (state%space%adapt%adapt_method == 'ALG2' .and. ( .not. state%first_GMRES_conv) ) ) then
!
!             if (state%space%adapt%adapt_level > 0 .and. state%time%iter_loc == 1) then  !computation of momentums in the first iteration step in an adaptation
!                 !Print*, '*************V cyklu*********'
!                 if (.not. state%loc_RTN(state%space%deg)%defined ) then
!                    Print*, 'Error, loc_RTN not defined yet!!'
!                    stop
!                 endif
!
!                 do i=1, grid%nelem
!                    call InitElemRTNMomentumCD(grid%elem(i), state%space%deg, SetRTNdof(state%space%deg) )
!                 enddo
!             endif
!          endif

          if( state%MGsolver ) then
             !!!call InitMG() moved

              ! *** TEST iLU as a solver ***
              !call RunILUSolu(Newton%x,Newton%b,state%nsize)
              !print*,L2Res(Newton%x,Newton%b,state%nsize,bMVprod)

              ! *** COMPARE Jacobi and G-S ***
              !allocate(y(state%nsize))
              !do  i=1,50,1
              !    call bJacobi(Newton%x,Newton%b,state%nsize,y)
              !    Newton%x = y
              !    !call bGS(Newton%x,Newton%b,state%nsize)
              !    print*,i,L2Res(Newton%x,Newton%b,state%nsize,bMVprod)
              !end do
              !deallocate(y)

              !print*,'--STOP'
              !STOP

              ! *** ILU for MG *** part of InitrMG
              !do i=1,grid%nelem
              !   call InitMGMatrixElementShape(grid%elem(i))
              !end do
              !call ComputeBlockMGILUPrecond()


             !call cpu_time(t1)
             call SolveBlockLinearProblem(state%nsize, eta, Newton%b, Newton%x, &
                  state%linSolver%residuum, state%linSolver%iter, state%linSolver%lin_solver_not_conv)
             !call cpu_time(t2)

            !\linear multigrid
          else
             !write(*,'(a6,4es14.6)')'AAA0', eta, Newton%res, state%linSolver%tol

             !MatrixA = 'MatrixA'
             !open(21, file=MatrixA, status='UNKNOWN', position = 'append')
             !do i=1,grid%nelem
             !   do j=1, grid%elem(i)%dof
             !      write(21,'(15es14.6)') grid%elem(i)%Mass%Mb(j, 1:grid%elem(i)%dof)
             !   enddo
             !enddo
             !write(21,'(a10)') '**********'
             !close(21)



             !MatrixB = 'MatrixB'
             !open(23, file=MatrixB, status='UNKNOWN', position = 'append')
             !do i=1,grid%nelem
             !   do j=1, grid%elem(i)%dof
             !      write(23,'(15es14.6)') grid%elem(i)%block(0)%Mb(j, 1:grid%elem(i)%dof)
             !   enddo
             !enddo
             !write(23,'(a10)') '**********'
             !close(23)

             !MatrixC = 'MatrixC'
             !open(24, file=MatrixC, status='UNKNOWN', position = 'append')
             !do i=1,grid%nelem
             !   do j=1, grid%elem(i)%flen
             !      k = grid%elem(i)%face(neigh,j)
             !      if (k > 0) then
             !         do l=1, grid%elem(i)%dof
             !            write(24,'(15es14.6)') grid%elem(i)%block(j)%Mb(l, 1:grid%elem(k)%dof)
             !         enddo
             !      endif
             !   enddo
             !enddo
             !write(24,'(a10)') '**********'
             !close(24)


             ! wet steam
             !write(*,*) "---------------------------------------------"
             !write(*,*) "PerformOneTimeStep, Newton iteration: ", iter
             !write(*,*) 'NNNNNN, iter', iter
             !call WriteMatrixA(eta)   ! NaN wet steam testing
             !write(*,*) 'LLLLLLL'

             !do i=1,6
             !   write(*,'(a5,i5,100es12.4)') 'MASS', i, grid%elem(1)%Mass%Mb(i,:)
             !enddo

             !call WriteMatrixA(eta)

             !!call cpu_time(tt)
             !!print*,'CT,',tt,' impl LSb'

             call SolveBlockLinearProblem(state%nsize, eta, Newton%b, Newton%x, &
                  state%linSolver%residuum, state%linSolver%iter, &
                  state%linSolver%lin_solver_not_conv)

             !!call cpu_time(tt)
             !!print*,'CT,',tt,' impl LSe'



!             open(73, file = 'linSolver2', action="write" )
!
!             print*, 'HERE AI', eta
!             write(73, * ) state%nsize
!             write(73, * ) eta
!             write(73, * ) Newton%b
!             write(73, * ) Newton%x
!
!
!             close(73)

             !stop 'afterSolveBlockLin problem'

             !!write(*,'(a6,5es14.6)')'AAA1', eta, time_prepare, time_solve, t1, t2

             !RHSideb = 'RHSideb'
             !open(22, file=RHSideb, status='UNKNOWN', position = 'append')
             !do i=1,state%nsize
             !   write(22,'(es14.6)') Newton%b(i)
             !enddo
             !write(22,'(a10)') '**********'
             !close(22)

          endif

          call cpu_time(t1)
          time_solve = time_solve + t1 - t2

          if(state%time%time_method == 'S') then ! semi-implicit scheme (NO damping !!)

             call bMVprod(Newton%rr, Newton%x, state%nsize)    ! rr = Ax

             Newton%rr(:) = Newton%rr(:) - Newton%b(:)                ! rr = Ax - b

             ! update of the solution
             k = 0
             do i=1,grid%nelem
                elem => grid%elem(i)
                j = elem%dof*ndim
                elem%w(0,1:j) = elem%w(0,1:j)  + Newton%x(k+1:k+j)

                ! components of th residual vector
                elem%vec( res_vec, 1:j) = - Newton%rr(k+1:k+j)

                k = k+j
             enddo


          else  ! fully implicit method

             ! residuum of the linear algebraic system
             eta = state%time%alpha(0)/state%time%tau(1)
             if(state%time%cn ) eta = eta*2  ! Crank-Nicolson

             !write(*,*) 'Newton%x(:) 1', Newton%x(:)  ! wet steam
             call bMVprod(Newton%rr, Newton%x, state%nsize)    ! rr = Ax
             !write(*,*) 'Newton%x(:) 2', Newton%x(:)  ! wet steam

             Newton%rr(:) = Newton%rr(:) - Newton%b(:)                ! rr = Ax - b


             ! some manipulation for the fail of teh computation
             if(state%linSolver%lin_solver_not_conv > 0) then
                open(53, file='GMRES_faild', status='UNKNOWN', position='append')
                res_max_val = 0.
                do i = 1, grid%nelem
                   elem => grid%elem(i)

                   kk = elem%ncv - 1
                   elemDof = elem%dof * ndim
                   elem%vec( res_vec, 1:elemDof) =  Newton%rr(kk+1:kk+ elemDof)

                   res_max_val = max( res_max_val, &
                        sqrt( dot_product(elem%vec( res_vec, 1:elemDof), elem%vec( res_vec, 1:elemDof) )))
                enddo

             endif

             if(state%linSolver%lin_solver_not_conv > 0 .and. state%linSolver%residuum > 1.) then
                print*,'Bad BDF iteration', state%linSolver%residuum
                state%nlSolver%implicitly = .true.
                return
             endif


             !rel_res = VectorPrecondNorm(Newton%rr(:) ) / Newton%res
             !rel_res = state%linSolver%residuum / VectorPrecondNorm(Newton%b(:) )

             ! seeking of optimal damping parameter
             Newton%lambda = 1.0   ! initialization of lambda (= damping factor)

             do l=1,10    ! iterations, seeking the optimal damping factor
             !do l=1,1    ! iterations, seeking the optimal damping factor
                ! update of the solution

                !?? write(*,*) 'w', elem%w(0,1:ndim)  ! wet steam

                k = 0
                do i=1,grid%nelem
                   elem => grid%elem(i)
                   j = elem%dof*ndim
                   !write(*,*) 'Newton%x(k+1:k+j)', Newton%x(k+1:k+j)  ! wet steam
                   elem%w(0,1:j) = elem%w(0,1:j)  &
                        + (Newton%lambda-Newton%lambda_old) * Newton%x(k+1:k+j)

                   ! components of th residual vector
                   elem%vec( res_vec, 1:j) = - Newton%rr(k+1:k+j)


                   k = k+j
                enddo
                Newton%lambda_old = Newton%lambda

                !print*,'Compute ONLY vector q(w)', state%nlSolver%implicitly
                !write(600,*) '# Compute ONLY vector q(w)', state%nlSolver%implicitly

                !write(*,*) "Newtonova iterace 5", iter
                !write(*,*) 'Newtonova iterace 5, w ', elem%w(0,1:ndim)  ! wet steam

                !!call cpu_time(tt)
                !!print*,'CT,',tt,' impl 4b'

                !call cpu_time(tt)
                !write(*,'(a30, 2es12.4)') 'CT_starts 5 :', tt - state%start_time

                call ComputeTerms(deg_plus )

                !call cpu_time(tt)
                !write(*,'(a30, 2es12.4)') 'CT_ends   5 :', tt - state%start_time
                !print*

                !!print*,'CT,',tt,' impl 4e'
                !write(*,*) 'sfs;ldfjsf'

                call FillVector(Newton%b1(1:state%nsize),eta ) !include time derivative term

                Newton%res0 = VectorPrecondNorm( Newton%b1 )
                Newton%theta  =  Newton%res0  / Newton%res


                !write(*,'(a8,i5,20es14.6)') 'newton:',iter,  Newton%res,  Newton%res0,  Newton%theta

                !call cpu_time(t2)

                write(55,'(3i4,13es11.3, i5,2l3,F10.1)') state%time%iter,Newton%iter,l, &   ! 1..3
                     state%time%iter + 1.*(Newton%iter-1)/(Newton%max_iter+2), &      ! 4
                     state%space%adapt%adapt_level + 3.*state%nlSolver%Aiter /(Newton%max_iter+2)/state%time%maxiter, &
                     Newton%res/state%nsize,  Newton%res0/state%nsize, &         ! 6..7
                     Newton%res,  Newton%res0, &                                 ! 8..9
                     rel_res, state%linSolver%tol,  Newton%lambda, Newton%theta, &        !10 .. 13
                     state%L_estim(resA:resT),  &                                !14 .. 17
                     state%linSolver%iter, loc_implicitly, &
                     state%linSolver%precond_update, t2- time_solve

                if(loc_implicitly) then
                   write(54,'(3i4,12es11.4, i5,2l3)') state%time%iter,Newton%iter,l, &   ! 1..3
                        state%time%iter + 1.*(Newton%iter-1)/(Newton%max_iter+2), &      ! 4
                        state%space%adapt%adapt_level+3.*state%nlSolver%Aiter/(Newton%max_iter+2)/state%time%maxiter, &
                        Newton%res/state%nsize,  Newton%res0/state%nsize, &         ! 6..7
                        Newton%res,  Newton%res0, &                                 ! 8..9
                        rel_res, state%linSolver%tol, Newton%lambda, Newton%theta,  &
                        state%L_estim(resS), state%L_estim(resA), &
                        state%linSolver%iter,loc_implicitly,state%linSolver%precond_update

                endif

                Newton%newton_count = Newton%newton_count + 1

                !if(loc_implicitly) print*,'update was applied'

                ! residuum is not decreasing, reduce lambda
                if(Newton%theta >= 1 .and. Newton%res0 > 1E-8 ) then
                !if(Newton%theta > 1 - Newton%lambda/4  .and. Newton%res0 > 1E-8 )then
                   !Newton%lambda = min(mu1, 0.5* Newton%lambda)
                   Newton%lambda = 0.5* Newton%lambda
                else
                   Newton%lambda1 = min(1., 1.5*Newton%lambda )
                   !
                   !if(Newton%lambda1 >= 4 * Newton%lambda) then  ! lambda was to short
                   !   Newton%lambda = Newton%lambda1
                   !else
                      !if(Newton%theta > 1./4) then
                   !   state%nlSolver%implicitly = .true.
                   !endif
                   !Newton%lambda = Newton%lambda1

                   goto 15
                   !endif
                endif

                !if(Newton%lambda < 1E-2) goto 15  ! too small lambda, NOT CONVERGE
                if(Newton%lambda < 2E-1) goto 15  ! too small lambda, NOT CONVERGE

             enddo   !l=1,10    ! iterations, seeking the optimal damping factor

             !!state%nlSolver%implicitly = .true. ! last Newton iteration was not sucessfull
15           continue

             if(Newton%newton_count >= 8) state%nlSolver%implicitly = .true.
             !if(newton_count >= 4) state%nlSolver%implicitly = .true.
             !!!if(Newton%theta > 0.75) state%nlSolver%implicitly = .true.
             if(Newton%theta > 0.5) state%nlSolver%implicitly = .true.

             !if(Newton%lambda < 0.9) state%nlSolver%implicitly = .true. !lambda is small, C(w) update
             if(state%linSolver%lin_solver_not_conv > 0) state%nlSolver%implicitly = .true.

             !if(state%nlSolver%implicitly ) write(83,'(a8,i5, a14,i5,a8,es12.4,a8,l2 )')&
             !     'iter = ',state%time%iter, '### Ncount =',newton_count, &
             !     ',  theta = ', Newton%theta,',  non conv', state%linSolver%lin_solver_not_conv


             ! nonlinear algebraic residuum criterion
             if(Newton%non_alg_stop == 'aRES') then
                !!if(Newton%tol2 > 0. .and. (state%space%adapt%adapt_method == 'RES' .or. state%space%adapt%adapt_method == 'ANI' &
                !!  .or. state%space%adapt%adapt_method == 'Ahp') ) then

                call cpu_time(t2)
                time_prepare = time_prepare + t2 - t1

                !!call cpu_time(tt)
                !!print*,'CT,',tt,' rez 5b'

                !call cpu_time(tt)
                !write(*,'(a30, 2es12.4)') 'CT_starts 6 :', tt - state%start_time

                call RezidErrorEstimates( .true., .false. )

                !call cpu_time(tt)
                !write(*,'(a30, 2es12.4)') 'CT_ends   6 :', tt - state%start_time

                !!call cpu_time(tt)
                !!print*,'CT,',tt,' rez 5e'

                call cpu_time(t1)
                time_estim = time_estim + t1 - t2

                open(91, file='criter', status='unknown', position='append')
                if(Newton%iter == 1) write(91,'(x)')
                if(Newton%Aiter == 1) write(91,'(x)')
                write(91,'(2i6, 50es12.4)' ) state%time%iter, Newton%iter, &
                     state%time%iter + 1.*(Newton%iter-1)/(Newton%max_iter+2), &
                     state%space%adapt%adapt_level + 3.* Newton%Aiter /(Newton%max_iter+2)/state%time%maxiter, &
                     state%linSolver%tol, &   !5
                     state%L_estim(resA:resST),  &
                     state%estim(resA_ST,1), state%estim(resA_ST_loc,1),   &      ! 10..11
                     Newton%norm_res,  Newton%res0, Newton%res0/state%nsize, &      ! 12..14
                     Newton%tol, Newton%tol2,  t2 - state%start_time              ! 15..17


                close(91)

                !if(state%L_estim(resA)/state%G_estim(15) < Newton%tol2 ) goto 200
                ! VD
                !print*
                !write(*,'(a26,3i5,20es14.6)') '!AW! iter tol resA/resST ',state%time%iter, Newton%iter, 0, &
                !     Newton%tol2,  state%estim(resA_ST,1)  !, & !  Newton%res/state%nsize, Newton%tol, &
                !!state%L_estim(resA), state%L_estim(resS) !, state%L_estim(resST)
                !print*,'###################################################################'

                if(state%estim(resA_ST,1) < Newton%tol2 ) then
                !if(state%estim(resA_ST_loc,1) < Newton%tol2 ) then
                   call cpu_time(t2)
                   time_prepare = time_prepare + t2 - t1
                   !!write(*,'(a6, 3es12.4)') 'timesE',time_prepare, time_solve , time_estim
                   goto 200
                endif

             else  ! Newton%non_alg_stop == 'precL2'

                ! other techniques
                !if(Newton%res0/Newton%res_ref < Newton%tol2 ) goto 200

                !if(ndim == 4) &
                !     call DirectComputeDrag_Lift(cDLM(Newton%iter, 1:5 )  )
                !diffs = ( (cDLM(Newton%iter, 1)  - cDLM(Newton%iter-1, 1) )**2 &
                !     + (cDLM(Newton%iter, 2)  - cDLM(Newton%iter-1, 2) )**2  &
                !     + (cDLM(Newton%iter, 3)  - cDLM(Newton%iter-1, 3) )**2 )**0.5
                !
                !open(44, file='cDLM', status='UNKNOWN', position = 'append')
                !write(44,'(2i5,15es11.3)') state%time%iter, Newton%iter, &
                !     state%time%iter + 1.*(Newton%iter)/(Newton%max_iter+1), &
                !     cDLM(Newton%iter, 1:3), &
                !     abs(cDLM(Newton%iter, 1)  - cDLM(Newton%iter-1, 1) ), &
                !     abs(cDLM(Newton%iter, 2)  - cDLM(Newton%iter-1, 2) ), &
                !     abs(cDLM(Newton%iter, 3)  - cDLM(Newton%iter-1, 3) ), diffs
                !close(44)

                if(Newton%res0 < Newton%tol * state%nsize )  then
                   call cpu_time(t2)
                   time_prepare = time_prepare + t2 - t1
                   goto 200
                endif

             endif  ! end of other technique
          endif  ! end of semi-implicit / implicit

       endif    ! end of explicit

       ! go to the next time step
       Newton%res1 = Newton%res

       call cpu_time(t2)
       time_prepare = time_prepare + t2 - t1
       !write(*,'(a6,4es14.6)')'AAA3', time_prepare, time_solve, t1, t2
       !!write(*,'(a6, 3es12.4)') 'timesF',time_prepare, time_solve , time_estim

    enddo  ! end of Newton iterations: Newton%iter = 1, Newton%max_iter

200 continue

    write(55,'(3i4,12es11.3, i5,2l3)') state%time%iter,Newton%iter,l, &   ! 1..3
         state%time%iter + 1.*(Newton%iter+0)/(Newton%max_iter+2), &      ! 4
         state%space%adapt%adapt_level + 3.*(state%nlSolver%Aiter+1) /(Newton%max_iter+2)/state%time%maxiter, &
         Newton%res0/state%nsize,  Newton%res0/state%nsize, &         ! 6..7
         Newton%res0,  Newton%res0

    if(loc_implicitly) &
         write(54,'(3i4,12es11.4, i5,2l3)') state%time%iter,Newton%iter,l, &   ! 1..3
         state%time%iter + 1.*(Newton%iter-1)/(Newton%max_iter+2), &      ! 4
         state%space%adapt%adapt_level+3.*state%nlSolver%Aiter/(Newton%max_iter+2)/state%time%maxiter, &
         Newton%res/state%nsize,  Newton%res0/state%nsize, &         ! 6..7
         Newton%res,  Newton%res0, &                                 ! 8..9
         rel_res, state%linSolver%tol, Newton%lambda, Newton%theta,  &
         state%L_estim(resS), state%L_estim(resA),state%linSolver%iter,loc_implicitly,state%linSolver%precond_update


    !write(*,'(a6,3i5,20es12.4)') '????',state%time%iter, Newton%iter, 0, &
    !     state%estim(min_resA_ST,1), Newton%tol2 ,  Newton%res/state%nsize, Newton%tol


    ! space, time, space-time error estimates
    if( state%space%estim_space == 'RES' .or. state%time%estim_time == 'tRES') then

       !call cpu_time(tt)
       !write(*,'(a30, 2es12.4)') 'CT_starts 7 :', tt - state%start_time

       call RezidErrorEstimates( .false., .false. )

       !call cpu_time(tt)
       !write(*,'(a30, 2es12.4)') 'CT_ends   7 :', tt - state%start_time
       !print*


       call cpu_time(t1)
       time_estim = time_estim + t1 - t2
    endif
    !!write(*,'(a6, 3es12.4)') 'timesX',time_prepare, time_solve , time_estim

    if(state%time%time_method == 'S') then
       state%nlSolver%implicitly = .false.
       !print*,'Computation of the residuum'

       !write(*,*) "Newtonova iterace 6", iter
       call ComputeTerms(deg_plus )

       state%nlSolver%implicitly = .true.
       Newton%norm_res = EvalSSresidExplicit( )
       !Newton%norm_res = EvalSSresid( )

       !!call CheckResiduum( )

    elseif(state%time%time_method == 'I') then
       loc_implicitly = state%nlSolver%implicitly
       state%nlSolver%implicitly  = .false.

       !write(*,*) "Newtonova iterace 7", iter
       call ComputeTerms(deg_plus )

       Newton%norm_res = EvalSSresidExplicit( )
       state%nlSolver%implicitly = loc_implicitly
    endif

    !write(*,*) ' #### Total number of newton updates = ', Newton%updates, Newton%norm_res
    !write(111,*) state%time%iter+1, Newton%norm_res
    if(state%time%time_method == 'I') then
       write(54,*) ' #### Total number of newton updates = ', Newton%updates
       write(55,*) ' '
       close(54)
       close(55)
       !close(56)
    endif

    deallocate(cDLM)

    !if(mod(state%time%iter, 3) == 0) then
    !   state%isol = state%isol + 1
    !   call WriteProgressOutput( 'A' )
    !endif

    state%CPU_prepare = state%CPU_prepare  + time_prepare
    state%CPU_solve   = state%CPU_solve    + time_solve
    state%CPU_estim   = state%CPU_estim    + time_estim

    !VD
    !write(*,'(a6, 4es12.4)') 'timXXX',time_prepare, time_solve , time_estim, &
    !     time_prepare+time_solve+time_estim
    !write(*,'(a6, 3es12.4)') 'timXXX',state%CPU_prepare, state%CPU_solve , state%CPU_estim

    end associate ! Newton

  end subroutine PerformOneTimeStep


  !> perform the one time step of STDGM on the given grid
  subroutine PerformOneSTDGMstep()
    class(element), pointer :: elem
!    real, dimension(:,:), allocatable :: cDLM
    real:: t1, t2, rel_res, lambda_loc
    real :: tt, tt1
    real :: res_max_val, val, lam_old
    integer :: iter
    integer :: i, j, k, l, m, kk , ll, max_l
    integer :: elemDof
    logical :: vector_update, loc_implicitly, update
    integer :: imp
    character(len=7) :: Newtonx, MatrixA, RHSideb, MatrixB, MatrixC
    real :: some_criter, etaD
    integer :: restart, nloops, nDamping, iterDWR
    logical :: deg_plus, newtonDone
    logical :: performNewtonIter
    real :: time, lost

    ! WARNING Newton is asociated to state%nlSolver i.e. equals in this subroutine !!!
    associate ( Newton => state%nlSolver )
    select type ( Newton )
      type is ( Newton_t )

    !print*," !if(Newton%non_alg_stop == 'aRES') deg_plus = .true.  !!!! "

    call state%cpuTime%startPrepareTime()

!    allocate(cDLM(0:Newton%max_iter, 1:5))
    Newton%norm_res = 1.
    Newton%lambda = 1.0   ! initialization of lambda (= damping factor)
    Newton%lambda1 = 1.0   ! initialization of lambda (= damping factor)

    ! minimal update of the flux matrix after imp time steps
    imp = Newton%min_update
    if(state%modelName == 'scalar' .or.state%modelName == '2eqs') imp = max(1, imp)

    if( state%modelName == 'porous' )  imp = max(1, imp)
    !if( state%modelName == 'porous' ) &
    !     print*,"if( state%modelName == 'porous' )  imp = max(1, imp)"


    ! Has the big vector the right size? If not (e.g., after mesh adapation) then resize
    if(size(Newton%b(:) ) /= state%nsize  )then
       deallocate( Newton%b, Newton%b1, Newton%x, Newton%rr )
       allocate(Newton%b(1:state%nsize), Newton%b1(1:state%nsize), &
            Newton%x(1:state%nsize), Newton%rr(1:state%nsize) )
    endif

    open(54, file='newton-update', status='unknown', position='append')
    open(55, file='newton', status='unknown', position='append')

    write(55,*) '  '
    write(55,*) '## it  i  l    ritr     ritr2   ', &
         '   |F(x)|/N  |F(x*)|/N   |F(x)|   |F(x*)|      rel_res ', &
         '     gmr_tol    lambda     theta       etaS    algeb   LAiter IM Pr'

    Newton%newton_count = 0
    Newton%iter  = 0

    ! clean and init DWR and also DWR%aDWR when aDWR algebraic criterion is used
    if (state%space%estim_space == 'DWR') then
      call DWR%update(grid)
      call allocateZST( grid, DWR%deg_plus )
    end if

    deg_plus = .false.
    if(Newton%non_alg_stop == 'aRES') deg_plus = .true.

    state%time%keep_tau = .true.
    iter = 1
    call state%cpuTime%addPrepareTime( )

    ! NEWTON CYCLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(state%nlSolver%non_alg_stop == 'aDWR') then
      ! outer DWR cycle
      iterDWR = 0
      do while ( DWR%estimNL > DWR%aDWR%nlTol )
         ! update the tolerance for Newton
         state%nlSolver%tol = DWR%aDWR%nlTol * DWR%aDWR%C_L
         ! Newton loop with tolerance state%nlSolver%tol
         call NewtonSolve( grid, imp, iter, deg_plus, DWR)
         call computeDWRestimates( DWR, grid )
         call DWR%writeNlErrorFile (grid%nelem)
         ! disable infinite cycle
         iterDWR = iterDWR + 1
         if ( iterDWR > 10 ) &
            stop 'Too much of outer DWR iterations'
      end do ! while DWR
    else
      call NewtonSolve( grid, imp, iter, deg_plus)
    endif

    call state%cpuTime%startPrepareTime()

    do i = 1, grid%nelem
       elem => grid%elem(i)
       call Transfer_wST_to_w_Elem( elem , 0, elem%TQnum)
    enddo

    !call Newton_output(55, iter+1, loc_implicitly )
    !if(loc_implicitly) call Newton_output(54, iter+1, loc_implicitly )

    loc_implicitly = state%nlSolver%implicitly
    state%nlSolver%implicitly  = .false.

    Newton%norm_res = EvalSSresidExplicit( )
    state%nlSolver%implicitly = loc_implicitly

    write(54,*) ' #### Total number of newton updates = ', Newton%updates
    write(55,*) ' '
    close(54)
    close(55)

!    deallocate(cDLM)

   class default
         stop 'nonNewton type of nlSolver in Newton method.'
   end select
   end associate

    call state%cpuTime%addPrepareTime( )

    !print*, 'End PerformOneSTDGMstep'
  end subroutine PerformOneSTDGMstep

  !> perform the one time step of STDGM on the given grid
  !> using the Anderson acceleration
  subroutine PerformOneSTDGMstep_Andrerson( time_prepare, time_solve, time_estim )
    real, intent(inout) ::  time_prepare, time_solve, time_estim
    class(element), pointer :: elem
    real, dimension(:,:), allocatable :: cDLM
    real:: t1, t2, rel_res, lambda_loc
    real :: tt, tt1
    real :: res_max_val, val, lam_old
    !real :: mu, mu1
    !real :: normL2, normH1
    integer :: iter, newton_count
    integer :: i, j, k, l, m, kk , ll, max_l
    integer :: elemDof
    logical :: vector_update, loc_implicitly, update
    integer :: imp
    character(len=7) :: Newtonx, MatrixA, RHSideb, MatrixB, MatrixC
    real :: some_criter, etaD
    integer :: restart, nloops
    logical :: deg_plus
    real, dimension(:), allocatable :: alpha   ! for Anderson acceleration
    real, dimension(:,:), allocatable :: dd, dd1
    integer, parameter :: AA_F = 1  ! update d^{k+1} = w^{k+1} - w^{k}
    integer, parameter :: AA_w = 2  ! result of update w^{k+1}
    integer, parameter :: AA_m = 2   ! the maximal value AA_F, ...
    integer :: mk
    associate ( Newton => state%nlSolver )
      !! WARNING Newton is asociated to state%nlSolver i.e. equals in this subroutine !!!

    deg_plus = .false.
    if(Newton%non_alg_stop == 'aRES') deg_plus = .true.
    !print*," !if(Newton%non_alg_stop == 'aRES') deg_plus = .true.  !!!! "

    time_prepare = 0.
    time_solve  = 0.
    time_estim = 0.

    ! Anderson acceleration
    Newton%AA_max = 1  ! NO acceleration
    Newton%AA_max = 4

    allocate( Newton%AA_x(1:AA_m, 0:Newton%AA_max, 1:state%nsize),  Newton%AA_update( 1:state%nsize) )

    Newton%AA_x = 0.
    allocate( dd(0:Newton%AA_max-1, 0:Newton%AA_max-1 ), dd1(1:Newton%AA_max, 1:Newton%AA_max ) )
    allocate(alpha (1:Newton%AA_max+1 ))
    dd = 0.

    allocate(cDLM(0:Newton%max_iter, 1:5))
    Newton%norm_res = 1.
    Newton%lambda = 1.0   ! initialization of lambda (= damping factor)
    Newton%lambda1 = 1.0   ! initialization of lambda (= damping factor)

    !open(63, file = 'stdgm-comp', status = 'UNKNOWN', position = 'append')

    !if(ndim == 4) call DirectComputeDrag_Lift(cDLM(0, 1:5 ) )
    !open(44, file='cDLM', status='UNKNOWN', position = 'append')
    !write(44,'(x)')
    !write(44,'(2i5,15es11.3)') state%time%iter, Newton%iter, &
    !     state%time%iter + 1.*0., cDLM(Newton%iter, 1:3)
    !close(44)

    ! minimal update of the flux matrix after imp time steps
    imp = Newton%min_update
    if(state%modelName == 'scalar' .or.state%modelName == '2eqs') imp = max(1, imp)

    if( state%modelName == 'porous' )  imp = max(1, imp)
    !if( state%modelName == 'porous' ) &
    !     print*,"if( state%modelName == 'porous' )  imp = max(1, imp)"


    ! Has the big vector the right size? If not (e.g., after mesh adapation) then resize
    if(size(Newton%b(:) ) /= state%nsize  )then
       deallocate( Newton%b, Newton%b1, Newton%x, Newton%rr )
       !write(63,*) 'Changing size of Newton%b'

       allocate(Newton%b(1:state%nsize), Newton%b1(1:state%nsize), &
            Newton%x(1:state%nsize), Newton%rr(1:state%nsize) )
    endif

      !WHY? Do we need this in STDGM too?
  !  if(imp > 0) then
  !     if( mod(state%time%iter_loc,imp) == 1 .or. imp == 1)  state%nlSolver%implicitly = .true.
  !  elseif(imp == 0) then
  !     state%nlSolver%implicitly = .true.
  !  endif

    open(54, file='newton-update', status='unknown', position='append')
    open(55, file='newton', status='unknown', position='append')

    write(55,*) '  '
    write(55,*) '## it  i  l    ritr     ritr2   ', &
         '   |F(x)|/N  |F(x*)|/N   |F(x)|   |F(x*)|      rel_res ', &
         '     gmr_tol    lambda     theta       etaS    algeb   LAiter IM Pr'


    newton_count = 0
    Newton%iter  = 0
    if ( Newton%non_alg_stop == 'aDWR' ) then
       !print*, 'cleanDWR called! '
      call DWR%clean()
      call DWR%aDWR%init()
    endif

    if (state%space%estim_space == 'DWR') &
         call InitDualProblem( DWR, DWR%deg_plus )

    if ( Newton%non_alg_stop == 'aDWR' ) then
       !print*, 'InitDualLinSolver in PerformOneSTDGMstep'
       call InitDualLinSolver( DWR, grid )

    endif


    state%time%keep_tau = .true.
    do iter = 1, Newton%max_iter

       !!!if(state%modelName == 'pedes' ) call Pedestrian_velocity_limiting_STDG( )
       !if(state%time_dependent  .and.  iter >=  0.75*Newton%max_iter ) state%time%keep_tau = .false.

       call cpu_time(t1)

       if(imp == 0)  state%nlSolver%implicitly = .true.

       loc_implicitly = state%nlSolver%implicitly
!       state%linSolver%precond_update = .false.
       update = .false.

       if( state%nlSolver%implicitly ) then

          !write(600,*)'# Compute matrix C(w), vector q(w)'

          update = .true.
          !call cpu_time(tt1)

          call ComputeSTDGM_Terms(deg_plus )

          !call WriteMatrixST_Blocks( 0. )
          !stop "0i034iwokjd3olwk"


          !call cpu_time(tt)
          !write(*,'(a30, 2es12.4)') 'CT_ends   1 :', tt - state%start_time,  tt - tt1

         ! print*, 'STDGMTerms done with implicitly = .true., iter =', iter
          newton_count = 0
          Newton%updates = Newton%updates + 1
          !state%linSolver%precond_update = .true.
          !print*, 'precond update set:', state%linSolver%precond_update, iter

          state%nlSolver%implicitly = .false.

          !call cpu_time(tt1)
          call ComputeSTDGM_Terms(deg_plus )
          !call cpu_time(tt)
          !write(*,'(a30, 2es12.4)') 'CT_ends   2 :', tt - state%start_time ,  tt - tt1
          state%nlSolver%implicitly = .true.

          vector_update = .true.
          state%nlSolver%implicitly = .false.

          ! computing with the matrices from previous time step ?
       else

          ! for iter > 1, array b(:) = F(x^k)
          if(iter == 1) then

             !call cpu_time(tt1)
             call ComputeSTDGM_Terms(deg_plus )
             !call cpu_time(tt)
             !write(*,'(a30, 2es12.4)') 'CT_ends   3 :', tt - state%start_time,  tt - tt1
            vector_update = .true.

          else

             vector_update = .false.
             Newton%b(1:state%nsize) = Newton%b1(1:state%nsize)
             Newton%res = Newton%res0
          endif

       endif ! end of if(state%nlSolver%implicitly)

       eta = 1./state%time%tau(1)


       if(vector_update) then  !filling elem%rhsST into global rhs

          call FillVectorST(Newton%b(1:state%nsize), eta ) ! value of eta is not used, only its sign

          !Newton%res  = VectorNorm(Newton%b)
          Newton%res = VectorPrecondNorm(Newton%b)
       endif


       if(iter == 1) Newton%res_ref = Newton%res

       ! first exit based on the absolute value of the residuum vector || F(w) ||
       if(iter > 1 .and. Newton%res/state%nsize < Newton%tol) then
       !   print*, 'goto200 - Newton%res/state%nsize < Newton%tol',Newton%res/state%nsize, Newton%tol
       !   print*, 'FR Commented newton stop - it should stop iterating because |F|<eps'
       !   print*,
          goto 200
       endif

       Newton%iter  = Newton%iter + 1
       Newton%Aiter = Newton%Aiter  + 1

       Newton%lambda_old = 0.
       Newton%x(1:state%nsize) = 0.

       if (.not. state%linSolver%tol_fixed) then
          if(iter == 1) then
             state%linSolver%tol = 0.25
          else
             state%linSolver%tol = min(0.25, abs(Newton%res - state%linSolver%residuum) / Newton%res1)
             state%linSolver%tol = max(state%linSolver%tol, 1E-6)
          endif
       endif


       call cpu_time(t2)
       time_prepare = time_prepare + t2 - t1

       !print*, 'Before Solve Linear Problem'

       !call WriteMatrixST_Blocks(eta)
       !call test_bMVprodST()
       !call Write_rhsST()

       !call cpu_time(tt1)
       !write(*,'(a30, 2es12.4)') 'CT_starts LA :', tt1 - state%start_time

       ! set number of iterations and number of nloops of lin alg. solver
       if (Newton%non_alg_stop == 'aDWR') then
         ! set number of iter after which restart is done
         restart = DWR%aDWR%restart_primal

!         print*, 'restart: ' , restart
         nloops = 1
         DWR%aDWR%iter = iter

       else
         restart = 45 !30
         nloops = 10 !5 ! 10
       endif


       if(state%model%varying_time_term) eta = 0. ! time deriv terms already included

       call cpu_time(tt1)
       !print*,'before LA'
       call SolveBlockLinearSTDGMProblem(state%nsize, eta, Newton%b, Newton%x, &
            state%linSolver%residuum, state%linSolver%iter, &
            state%linSolver%lin_solver_not_conv, restart, nloops)

       !!!!!!Newton%x(:) =  Newton%b(:) * 0.01

       call cpu_time(t1)
       time_solve = time_solve + t1 - t2
       !print*,'after LA', t1- t2

       !write(*,'(a40, a4, f8.2, f12.4)') '#CPU# CT_ends   LA :','    ', &
       !     t1 - state%start_time , t1 - tt1


       call bMVprodST(Newton%rr, Newton%x, state%nsize)    ! rr = Ax

       Newton%rr(:) = Newton%rr(:) - Newton%b(:)

       ! some manipulation for the fail of the computation
       if(state%linSolver%lin_solver_not_conv > 0) then
          open(53, file='GMRES_failed', status='UNKNOWN', position='append')
          res_max_val = 0.
          do i = 1, grid%nelem
             elem => grid%elem(i)

             kk = elem%ncv - 1
             elemDof = elem%dof * ndim * elem%Tdof
             elem%vec( res_vec, 1:elemDof) =  Newton%rr(kk+1:kk+ elemDof)

             res_max_val = max( res_max_val, &
                  sqrt( dot_product(elem%vec( res_vec, 1:elemDof), elem%vec( res_vec, 1:elemDof) )))
          enddo

       endif

       ! seeking of optimal damping parameter
       Newton%lambda_old = 0.

       Newton%lambda = 1.0   ! initialization of lambda (= damping factor)
       !Newton%lambda = Newton%lambda1
       !max_l = 10
       !max_l = 4
       max_l = 1

       !if(state%linSolver%lin_solver_not_conv > 0) then
       !   max_l = 1
       !   Newton%lambda = 0.1
       !endif

       do l=1,max_l   ! iterations, seeking the optimal damping factor
          ! update of the solution

          kk = 0
          do i = 1, grid%nelem
             !!lambda_loc = Newton%lambda
             elem => grid%elem(i)
             if (elem%ncv /= kk + 1) then
                print*, 'Problem while copying update to wST (compute.f90)'
                stop
             endif
             elemDof = elem%dof * ndim * elem%Tdof
             elem%vec( res_vec, 1:elemDof) =  Newton%rr(kk+1:kk+ elemDof)


             do k = 1, elem%Tdof
                do j = 1, ndim

                   elem%wST(j,1:elem%dof,k) = elem%wST(j,1:elem%dof,k) &
                        + (Newton%lambda - Newton%lambda_old) &
                        * Newton%x(kk + 1 : kk + elem%dof)

                   kk = kk + elem%dof
                   !write(*,'(a10, 2i5, 300es12.4)') 'b_sol:',elem%i, elemDof,  elem%wST(j,1:elem%dof,k)
                enddo !j
             enddo !k

             if(state%linSolver%lin_solver_not_conv > 0) then
                if(sqrt( dot_product(elem%vec( res_vec, 1:elemDof), elem%vec( res_vec, 1:elemDof))) > &
                     0.5 * res_max_val ) then
                   write(53,*), elem%xc(:), elem%i, &
                        sqrt(dot_product(elem%vec(res_vec, 1:elemDof), elem%vec( res_vec, 1:elemDof))),&
                        res_max_val, elem%wST(1:ndim,1:elem%dof,:)
                endif
             endif

          enddo !i


          if(state%linSolver%lin_solver_not_conv > 0 ) close(53)

          lam_old = Newton%lambda_old
          Newton%lambda_old = Newton%lambda



          !   print*,'Compute ONLY vector q(w)', state%nlSolver%implicitly
          !write(600,*) '# Compute ONLY vector q(w)', state%nlSolver%implicitly
          !   print*, 'after update wST implicitly:', state%nlSolver%implicitly

          !call cpu_time(tt1)
          !write(*,'(a30, 2es12.4)') 'CT_starts 4 :', tt1 - state%start_time


          call ComputeSTDGM_Terms(deg_plus )
          !	print*, 'after STDGMTerms'

          !call cpu_time(tt)
          !write(*,'(a30, 2es12.4)') 'CT_ends   4 :', tt - state%start_time ,  tt - tt1
          !write(*,'(a6,2i5)') '###', iter
          !write(14, *) state%time%iter_loc + 1.*iter/30, tt - tt1,  tt - state%start_time

          call FillVectorST(Newton%b1(1:state%nsize), eta )

          Newton%res0 = VectorPrecondNorm( Newton%b1 )

          !Newton%res0 = dot_product(Newton%b1(:), Newton%b1(:) )**0.5
          Newton%theta  =  Newton%res0  / Newton%res
          !write(63,*) 'Newton%res0', Newton%res0
          !write(63,*) 'Newton%res', Newton%res
          !write(63,*) 'Newton%theta:', Newton%theta

          !if(l==1) print*
          !write(*,'(a8, 3i5, 1es10.2, a1, 30es12.4)') 'newton:',state%time%iter, Newton%iter, l, &
          !    Newton%lambda,'|', &
          !    Newton%res0, Newton%res, Newton%theta, state%linSolver%residuum


          if(Newton%iter > 1) call Newton_output(grid, 55, l, loc_implicitly )

          if(loc_implicitly) call Newton_output(grid, 54, l, loc_implicitly )

          newton_count = newton_count + 1
          !      print*, '(related to seeking lambda???) newton count = ', newton_count

          ! residuum is not decreasing, reduce lambda
          if(Newton%theta >= 1 .and. Newton%res0 > 1E-8 ) then
             !if(Newton%theta > 1 - Newton%lambda/4  .and. Newton%res0 > 1E-8 )then
             !Newton%lambda = min(mu1, 0.5* Newton%lambda)
             !Newton%lambda = Newton%lambda / 2.0
             Newton%lambda = Newton%lambda / 1.5
          else
             !Newton%lambda1 = min(1., 2.0*Newton%lambda )
             Newton%lambda1 = min(1., 5.*Newton%lambda )

             goto 15
          endif

          ! too small lambda, NOT CONVERGE
          if(Newton%lambda < 1E-2 .and. state%modelName /= 'porous') goto 15
          !if(Newton%lambda < 1E-3) goto 15  ! too small lambda, NOT CONVERGE

       enddo   !l=1,10    ! iterations, seeking the optimal damping factor

       !!state%nlSolver%implicitly = .true. ! last Newton iteration was not sucessfull
15     continue


       !! ANDERSON ACCELERATION
       ! maximal number of used vectors
       mk = min( Newton%iter, Newton%AA_max) - 1
       ! mk = 0

       ! storing for the vectors for the Anderson acceleration from the previous iterations
       do l =  mk, 1, -1
          Newton%AA_x(1:AA_m, l, 1:state%nsize) = Newton%AA_x(1:AA_m, l-1, 1:state%nsize)
       enddo

       ! shifting of the dd matrix terms
       if(mk > 0) then
          dd(1:mk, 1:mk) = dd( 0:mk-1, 0:mk-1)
       endif

       !! new values for the  Anderson acceleration:
       !  update d^{k+1} = w^{k+1} - w^{k} = G(w^k) - w^k
       Newton%AA_x(AA_F, 0, 1:state%nsize) =  Newton%x(1:state%nsize) * Newton%lambda

       ! result of updated w^{k+1}= G(w^k)
       call FillVectorST_from_wST( Newton%AA_x(AA_w, 0, 1:state%nsize) )

       ! Anderson update
       if(mk == 0) then
          ! only computation of the scalar product
          dd(0,0) = dot_product(Newton%AA_x(AA_F, 0, :), Newton%AA_x(AA_F, 0, : )  )
          !write(*,'(a8,i5, 20es12.4)') 'dd:',0, dd(0, :)

       else

          ! do i=0, mk
          !    write(*,'(a8,i5, 200es12.4)') 'F(x^k):',i,Newton%AA_x(AA_F, i, :)
          ! enddo
          ! print*
          ! do i=0, mk
          !    write(*,'(a8,i5, 200es12.4)') 'w(x^k):',i,Newton%AA_x(AA_w, i, :)
          ! enddo
          ! print*

          !do i=0, 0 ! mk
          i = 0
          do j=i, mk
             dd(i,j) = dot_product(Newton%AA_x(AA_F, i, :), Newton%AA_x(AA_F, j, : )  )
             dd(j,i) = dd(i,j)
          enddo
          !enddo
          !stop"7e39ud3o"

          !do j=i, mk
          !   write(*,'(a8,i5, 20es12.4)') 'dd:',j, dd(j, :)
          !enddo

          dd1(1:mk+1, 1:mk+1) = dd(0:mk, 0:mk)
          alpha(1:mk+1) = 1.
          call SolveLocalMatrixProblem(mk+1, dd1(1:mk+1, 1:mk+1), 1, alpha(1:mk+1) )
          val = sum (alpha(1:mk+1))
          alpha(1:mk+1) = alpha(1:mk+1) / val
          !write(*,'(a8,i5, 20es12.4)') 'alpha:',99,  alpha(1:mk+1)

          Newton%AA_update = 0.
          do i=0,mk
             Newton%AA_update(:) = Newton%AA_update(:) + alpha(i+1) *  Newton%AA_x(AA_w, i, :)
             !write(*,'(a10, 300es12.4)') 'update:', Newton%AA_update(1:10)
          enddo

          ! setting the new solution
          call Fill_wST_from_VectorST( Newton%AA_update(:) )

          call ComputeSTDGM_Terms(deg_plus )

          call FillVectorST(Newton%b1(1:state%nsize), eta )
          ! the actual residuum
          Newton%res0 = VectorPrecondNorm( Newton%b1 )

          !write(*,'(a8, 3i5, 1es10.2, a1, 30es12.4)') 'anderson:',state%time%iter, Newton%iter, -1, &
          !     Newton%lambda,'|', &
          !     Newton%res0, Newton%res, Newton%theta, state%linSolver%residuum

          !print*,'___________________________________________________________________________'
          !stop "i30id03ik3ew"

       endif


       !if(state%nlSolver%implicitly ) write(83,'(a8,i5, a14,i5,a8,es12.4,a8,l2 )')&
       !     'iter = ',state%time%iter, '### Ncount =',newton_count, &
       !     ',  theta = ', Newton%theta,',  non conv', state%linSolver%lin_solver_not_conv


       !					state%isol = state%isol + 1
       !
       !					call WriteProgressOutput( 'STE' )

       !     print*, 'Before algebraic residuum criterion'

       !call WriteMatrixLU()


       if(newton_count >= 8) state%nlSolver%implicitly = .true.
       !if(newton_count >= 4) state%nlSolver%implicitly = .true.
       ! !!if(Newton%theta > 0.75) state%nlSolver%implicitly = .true.
       if(Newton%theta > 0.5) state%nlSolver%implicitly = .true.

       !if(Newton%lambda < 0.9) state%nlSolver%implicitly = .true. !lambda is small, C(w) update
       if(state%linSolver%lin_solver_not_conv > 0) state%nlSolver%implicitly = .true.




       ! nonlinear algebraic residuum criterion
       if(Newton%non_alg_stop == 'aRES') then

          call cpu_time(t2)
          time_prepare = time_prepare + t2 - t1

          ! print*, 'PerformOneSTDGMstep calling RezidErrorEstimate(true)'
          !call RezidErrorEstimates( .true., .false. )  !onlyAS =.true.  -- almost the same
          call RezidErrorEstimates( .false., .false. )

           do i=1, grid%nelem
              elem => grid%elem(i)
              if(elem%i == -1) then ! 556
                 !if(i==1) print*
                 write(*,'(a2, i5,i3, i2, 2es9.1, a1,  30es13.5))') &
                      'N:', elem%i,Newton%iter, l, &
                      elem%eta(resA:resS,1), '|', &
                      elem%wST(1,1:elem%dof,1) ,elem%wST(1,1:elem%dof,2)
              endif
           enddo


          call cpu_time(t1)
          time_estim = time_estim + t1 - t2
          !write(*,'(a6, 3es12.4)') 'timesD',time_prepare, time_solve , time_estim

          if(state%time%iter == -10) then
             do i=1,grid%nelem
                elem => grid%elem(i)
                write(100+Newton%iter , *) i, elem%xc, max(1E-15, elem%eta(resA,1)), &
                     max(1E-15, elem%eta(resS,1)), &
                     max(1E-15, elem%eta(resA,1)/ max(1E-15, elem%eta(resS,1)) )

                write(200+Newton%iter , *) i, elem%xc, max(1E-15, elem%eta(resA,1))
                write(200+Newton%iter , *) i, elem%xc, max(1E-15, elem%eta(resS,1))
                write(200+Newton%iter , '(x)' )
                write(200+Newton%iter , '(x)' )
                write(200+Newton%iter , '(x)' )

             enddo
          endif

          ! if(Newton%iter == 1) &
          !      write(*,'(a4,a4,5a12,a7,2a12)') '  ', 'iter', '|F(x^k)| ','|F(x^k+1)|', &
          !      'eta_A', 'eta_S', &
          !      'rel. err. ', '  LA it',&
          !      'LA resid', '   CPU(s) '

          ! write(*,'(a4,2i4,5es12.4,a2,i5,2es12.4,l2)') 'NEW:', Newton%iter, l-1, &
          !       Newton%res,Newton%res0, &
          !       state%L_estim(resA), state%L_estim(resST), &
          !       state%L_estim(resA)/ state%L_estim(resST),' |', &
          !       state%linSolver%iter, state%linSolver%residuum, t1 -  state%start_time,&
          !       loc_implicitly

          open(91, file='criter', status='unknown', position='append')
          if(Newton%iter == 1) write(91,'(x)')
          if(Newton%Aiter == 1) write(91,'(x)')
          call Newton_output(grid, 91, l, loc_implicitly )
          close(91)


          if(state%estim(resA_ST,1) < Newton%tol2 ) then
             !if(state%estim(resA_ST_loc,1) < Newton%tol2 ) then
             call cpu_time(t2)
             time_prepare = time_prepare + t2 - t1
             !!write(*,'(a6, 3es12.4)') 'timesE',time_prepare, time_solve , time_estim
             goto 200  ! Newton method has achieved the prescribed tolerance
          endif

       else if (Newton%non_alg_stop == 'aDWR') then
         print*, 'Compute aDWR estimate (elseif in Perform)'
         stop

       else ! Newton%non_alg_stop == 'rezL2'
          ! other techniques

          if(Newton%res0 < Newton%tol * state%nsize )  then
             !write(*,*) 'Reziduum:' , Newton%res0,  '<' , Newton%tol * state%nsize
             call cpu_time(t2)
             time_prepare = time_prepare + t2 - t1
             goto 200  ! Newton method has achieved the prescribed tolerance
          endif

          !write(*,*) 'Reziduum:' , Newton%res0,  '>' , Newton%tol * state%nsize

       endif  ! end of other technique
       !  endif  ! end of semi-implicit / implicit


       ! go to the next Newton iteration
       !write(63,*) '####################################################################'

       !write(*,*) 'Going to the' , iter+1 , '-th Newton iteration'
       !stop 'FR stopped in compute.f90'

       Newton%res1 = Newton%res

       call cpu_time(t2)
       time_prepare = time_prepare + t2 - t1
       !write(*,'(a6,4es14.6)')'AAA3', time_prepare, time_solve, t1, t2
       !!write(*,'(a6, 3es12.4)') 'timesF',time_prepare, time_solve , time_estim

    enddo  ! end of Newton iterations: Newton%iter = 1, Newton%max_iter


200 continue

    do i = 1, grid%nelem
       elem => grid%elem(i)

       call Transfer_wST_to_w_Elem( elem , 0, elem%TQnum)
    enddo

    !call Newton_output(55, iter+1, loc_implicitly )
    !if(loc_implicitly) call Newton_output(54, iter+1, loc_implicitly )


    if( state%space%estim_space == 'RES' .or. state%time%estim_time == 'tRES') then
       ! P++       print*,' estimates of the time and space-time error already done '
       !call RezidErrorEstimates( .false., .false. )
       call cpu_time(t1)
       time_estim = time_estim + t1 - t2
    endif

    loc_implicitly = state%nlSolver%implicitly
    state%nlSolver%implicitly  = .false.

    !call cpu_time(tt1)
    !write(*,'(a30, 2es12.4)') 'CT_starts 5 :', tt1 - state%start_time

    ! P++  print*,' not necessary'
    !call ComputeSTDGM_Terms(deg_plus )

    !call cpu_time(tt)
    !write(*,'(a30, 2es12.4)') 'CT_ends   5 :', tt - state%start_time ,  tt - tt1
    !write(*,'(a6,2i5)') '###', iter
    !write(15, *) state%time%iter_loc + 1.*iter/30, tt - tt1,  tt - state%start_time


    Newton%norm_res = EvalSSresidExplicit( )
    state%nlSolver%implicitly = loc_implicitly

    write(54,*) ' #### Total number of newton updates = ', Newton%updates
    write(55,*) ' '
    close(54)
    close(55)


    deallocate(cDLM)

    !if(mod(state%time%iter, 3) == 0) then
    !   state%isol = state%isol + 1
    !   call WriteProgressOutput( 'A' )
    !endif

    state%CPU_prepare = state%CPU_prepare  + time_prepare
    state%CPU_solve   = state%CPU_solve    + time_solve
    state%CPU_estim   = state%CPU_estim    + time_estim
    !!write(*,'(a6, 3es12.4)') 'timXXX',state%CPU_prepare, state%CPU_solve , state%CPU_estim


    deallocate(  Newton%AA_x,  Newton%AA_update, dd, dd1, alpha)

    end associate

    !print*, 'End PerformOneSTDGMstep'

  end subroutine PerformOneSTDGMstep_Andrerson


!  !> perform the one time step of STDGM on the given grid
!  subroutine PerformOneSTDGMstep_OLD( time_prepare, time_solve, time_estim )
!    real, intent(inout) ::  time_prepare, time_solve, time_estim
!    class(element), pointer :: elem
!    real, dimension(:,:), allocatable :: cDLM
!    real:: t1, t2, rel_res, lambda_loc
!    real :: tt, tt1
!    real :: res_max_val, val, lam_old
!    !real :: mu, mu1
!    !real :: normL2, normH1
!    integer :: iter, newton_count
!    integer :: i, j, k, l, m, kk , ll, max_l
!    integer :: elemDof
!    logical :: vector_update, loc_implicitly, update
!    integer :: imp
!    character(len=7) :: Newtonx, MatrixA, RHSideb, MatrixB, MatrixC
!    real :: some_criter, etaD
!    integer :: restart, nloops
!    logical :: deg_plus
!
!
!    associate ( Newton => state%nlSolver )
!      !! WARNING Newton is asociated to state%nlSolver i.e. equals in this subroutine !!!
!
!    deg_plus = .false.
!    if(Newton%non_alg_stop == 'aRES') deg_plus = .true.
!    !print*," !if(Newton%non_alg_stop == 'aRES') deg_plus = .true.  !!!! "
!
!    time_prepare = 0.
!    time_solve  = 0.
!    time_estim = 0.
!
!    allocate(cDLM(0:Newton%max_iter, 1:5))
!    Newton%norm_res = 1.
!    Newton%lambda = 1.0   ! initialization of lambda (= damping factor)
!    Newton%lambda1 = 1.0   ! initialization of lambda (= damping factor)
!
!    !open(63, file = 'stdgm-comp', status = 'UNKNOWN', position = 'append')
!
!    !if(ndim == 4) call DirectComputeDrag_Lift(cDLM(0, 1:5 ) )
!    !open(44, file='cDLM', status='UNKNOWN', position = 'append')
!    !write(44,'(x)')
!    !write(44,'(2i5,15es11.3)') state%time%iter, Newton%iter, &
!    !     state%time%iter + 1.*0., cDLM(Newton%iter, 1:3)
!    !close(44)
!
!    ! minimal update of the flux matrix after imp time steps
!    imp = Newton%min_update
!    if(state%modelName == 'scalar' .or.state%modelName == '2eqs') imp = max(1, imp)
!
!    if( state%modelName == 'porous' )  imp = max(1, imp)
!    !if( state%modelName == 'porous' ) &
!    !     print*,"if( state%modelName == 'porous' )  imp = max(1, imp)"
!
!
!    ! Has the big vector the right size? If not (e.g., after mesh adapation) then resize
!    if(size(Newton%b(:) ) /= state%nsize  )then
!       deallocate( Newton%b, Newton%b1, Newton%x, Newton%rr )
!       !write(63,*) 'Changing size of Newton%b'
!
!       allocate(Newton%b(1:state%nsize), Newton%b1(1:state%nsize), &
!            Newton%x(1:state%nsize), Newton%rr(1:state%nsize) )
!    endif
!
!      !WHY? Do we need this in STDGM too?
!  !  if(imp > 0) then
!  !     if( mod(state%time%iter_loc,imp) == 1 .or. imp == 1)  state%nlSolver%implicitly = .true.
!  !  elseif(imp == 0) then
!  !     state%nlSolver%implicitly = .true.
!  !  endif
!
!    open(54, file='newton-update', status='unknown', position='append')
!    open(55, file='newton', status='unknown', position='append')
!
!    write(55,*) '  '
!    write(55,*) '## it  i  l    ritr     ritr2   ', &
!         '   |F(x)|/N  |F(x*)|/N   |F(x)|   |F(x*)|      rel_res ', &
!         '     gmr_tol    lambda     theta       etaS    algeb   LAiter IM Pr'
!
!
!    newton_count = 0
!    Newton%iter  = 0
!    if ( Newton%non_alg_stop == 'aDWR' ) then
!       !print*, 'cleanDWR called! '
!      call DWR%clean()
!      call DWR%aDWR%init()
!    endif
!
!    if (state%space%estim_space == 'DWR') &
!         call InitDualProblem( DWR, DWR%deg_plus )
!
!    if ( Newton%non_alg_stop == 'aDWR' ) then
!       !print*, 'InitDualLinSolver in PerformOneSTDGMstep'
!       call InitDualLinSolver( DWR, grid )
!
!    endif
!
!
!    state%time%keep_tau = .true.
!    do iter = 1, Newton%max_iter
!
!       !!!if(state%modelName == 'pedes' ) call Pedestrian_velocity_limiting_STDG( )
!       !if(state%time_dependent  .and.  iter >=  0.75*Newton%max_iter ) state%time%keep_tau = .false.
!
!       call cpu_time(t1)
!
!       if(imp == 0)  state%nlSolver%implicitly = .true.
!
!       loc_implicitly = state%nlSolver%implicitly
!!       state%linSolver%precond_update = .false.
!       update = .false.
!
!       if( state%nlSolver%implicitly ) then
!
!          !write(600,*)'# Compute matrix C(w), vector q(w)'
!
!          update = .true.
!          !call cpu_time(tt1)
!
!          call ComputeSTDGM_Terms(deg_plus )
!
!          !call cpu_time(tt)
!          !write(*,'(a30, 2es12.4)') 'CT_ends   1 :', tt - state%start_time,  tt - tt1
!
!         ! print*, 'STDGMTerms done with implicitly = .true., iter =', iter
!          newton_count = 0
!          Newton%updates = Newton%updates + 1
!          !state%linSolver%precond_update = .true.
!          !print*, 'precond update set:', state%linSolver%precond_update, iter
!
!          state%nlSolver%implicitly = .false.
!
!          !call cpu_time(tt1)
!          call ComputeSTDGM_Terms(deg_plus )
!          !call cpu_time(tt)
!          !write(*,'(a30, 2es12.4)') 'CT_ends   2 :', tt - state%start_time ,  tt - tt1
!          state%nlSolver%implicitly = .true.
!
!          vector_update = .true.
!          state%nlSolver%implicitly = .false.
!
!          ! computing with the matrices from previous time step ?
!       else
!
!          ! for iter > 1, array b(:) = F(x^k)
!          if(iter == 1) then
!
!             !call cpu_time(tt1)
!             call ComputeSTDGM_Terms(deg_plus )
!             !call cpu_time(tt)
!             !write(*,'(a30, 2es12.4)') 'CT_ends   3 :', tt - state%start_time,  tt - tt1
!            vector_update = .true.
!
!          else
!
!             vector_update = .false.
!             Newton%b(1:state%nsize) = Newton%b1(1:state%nsize)
!             Newton%res = Newton%res0
!          endif
!
!       endif ! end of if(state%nlSolver%implicitly)
!
!       eta = 1./state%time%tau(1)
!
!
!       if(vector_update) then  !filling elem%rhsST into global rhs
!
!          call FillVectorST(Newton%b(1:state%nsize), eta ) ! value of eta is not used, only its sign
!
!          !Newton%res  = VectorNorm(Newton%b)
!          Newton%res = VectorPrecondNorm(Newton%b)
!       endif
!
!       !write(121, *) newton%b(1:state%nsize)
!
!       if(iter == 1) Newton%res_ref = Newton%res
!
!       ! first exit based on the absolute value of the residuum vector || F(w) ||
!       if(iter > 1 .and. Newton%res/state%nsize < Newton%tol) then
!       !   print*, 'goto200 - Newton%res/state%nsize < Newton%tol',Newton%res/state%nsize, Newton%tol
!       !   print*, 'FR Commented newton stop - it should stop iterating because |F|<eps'
!       !   print*,
!          goto 200
!       endif
!
!       Newton%iter  = Newton%iter + 1
!       Newton%Aiter = Newton%Aiter  + 1
!
!       Newton%lambda_old = 0.
!       Newton%x(1:state%nsize) = 0.
!
!       if (.not. state%linSolver%tol_fixed) then
!          if(iter == 1) then
!             state%linSolver%tol = 0.25
!          else
!             state%linSolver%tol = min(0.25, abs(Newton%res - state%linSolver%residuum) / Newton%res1)
!             state%linSolver%tol = max(state%linSolver%tol, 1E-6)
!          endif
!       endif
!
!
!       call cpu_time(t2)
!       time_prepare = time_prepare + t2 - t1
!
!       !print*, 'Before Solve Linear Problem'
!
!       !call WriteMatrixST_Blocks(eta)
!       !call test_bMVprodST()
!       !call Write_rhsST()
!
!       !call cpu_time(tt1)
!       !write(*,'(a30, 2es12.4)') 'CT_starts LA :', tt1 - state%start_time
!
!       ! set number of iterations and number of nloops of lin alg. solver
!       if (Newton%non_alg_stop == 'aDWR') then
!         ! set number of iter after which restart is done
!         restart = DWR%aDWR%restart_primal
!
!!         print*, 'restart: ' , restart
!         nloops = 1
!         DWR%aDWR%iter = iter
!
!       else
!         restart = 45 !30
!         nloops = 10 !5 ! 10
!       endif
!
!
!       if(state%model%varying_time_term) eta = 0. ! time deriv terms already included
!
!       call cpu_time(tt1)
!       !print*,'before LA'
!       call SolveBlockLinearSTDGMProblem(state%nsize, eta, Newton%b, Newton%x, &
!            state%linSolver%residuum, state%linSolver%iter, &
!            state%linSolver%lin_solver_not_conv, restart, nloops)
!
!       call cpu_time(t1)
!       time_solve = time_solve + t1 - t2
!       !print*,'after LA', t1- t2
!
!       !write(*,'(a40, a4, f8.2, f12.4)') '#CPU# CT_ends   LA :','    ', &
!       !     t1 - state%start_time , t1 - tt1
!
!
!       call bMVprodST(Newton%rr, Newton%x, state%nsize)    ! rr = Ax
!
!       Newton%rr(:) = Newton%rr(:) - Newton%b(:)
!
!       ! some manipulation for the fail of the computation
!       if(state%linSolver%lin_solver_not_conv > 0) then
!          open(53, file='GMRES_failed', status='UNKNOWN', position='append')
!          res_max_val = 0.
!          do i = 1, grid%nelem
!             elem => grid%elem(i)
!
!             kk = elem%ncv - 1
!             elemDof = elem%dof * ndim * elem%Tdof
!             elem%vec( res_vec, 1:elemDof) =  Newton%rr(kk+1:kk+ elemDof)
!
!             res_max_val = max( res_max_val, &
!                  sqrt( dot_product(elem%vec( res_vec, 1:elemDof), elem%vec( res_vec, 1:elemDof) )))
!          enddo
!
!       endif
!
!       !write(131, *) newton%x(1:state%nsize)
!
!       ! seeking of optimal damping parameter
!       Newton%lambda_old = 0.
!
!       !Newton%lambda = 1.0   ! initialization of lambda (= damping factor)
!       Newton%lambda = Newton%lambda1
!
!       max_l = 10
!
!       !if(state%linSolver%lin_solver_not_conv > 0) then
!       !   max_l = 1
!       !   Newton%lambda = 0.1
!       !endif
!
!       do l=1,max_l   ! iterations, seeking the optimal damping factor
!          ! update of the solution
!          kk = 0
!          do i = 1, grid%nelem
!             !!lambda_loc = Newton%lambda
!             elem => grid%elem(i)
!             if (elem%ncv /= kk + 1) then
!                print*, 'Problem while copying update to wST (compute.f90)'
!                stop
!             endif
!             elemDof = elem%dof * ndim * elem%Tdof
!             elem%vec( res_vec, 1:elemDof) =  Newton%rr(kk+1:kk+ elemDof)
!
!
!             do k = 1, elem%Tdof
!                do j = 1, ndim
!
!                   elem%wST(j,1:elem%dof,k) = elem%wST(j,1:elem%dof,k) &
!                        + (Newton%lambda - Newton%lambda_old) &
!                        * Newton%x(kk + 1 : kk + elem%dof)
!
!                   kk = kk + elem%dof
!                   !write(*,'(a10, 2i5, 300es12.4)') 'b_sol:',elem%i, elemDof,  elem%wST(j,1:elem%dof,k)
!                enddo !j
!             enddo !k
!
!             if(state%linSolver%lin_solver_not_conv > 0) then
!                if(sqrt( dot_product(elem%vec( res_vec, 1:elemDof), elem%vec( res_vec, 1:elemDof))) > &
!                     0.5 * res_max_val ) then
!                   write(53,*), elem%xc(:), elem%i, &
!                        sqrt(dot_product(elem%vec(res_vec, 1:elemDof), elem%vec( res_vec, 1:elemDof))),&
!                        res_max_val, elem%wST(1:ndim,1:elem%dof,:)
!                endif
!             endif
!
!          enddo !i
!
!
!          if(state%linSolver%lin_solver_not_conv > 0 ) close(53)
!
!          lam_old = Newton%lambda_old
!          Newton%lambda_old = Newton%lambda
!
!
!
!          !   print*,'Compute ONLY vector q(w)', state%nlSolver%implicitly
!          !write(600,*) '# Compute ONLY vector q(w)', state%nlSolver%implicitly
!          !   print*, 'after update wST implicitly:', state%nlSolver%implicitly
!
!          !call cpu_time(tt1)
!          !write(*,'(a30, 2es12.4)') 'CT_starts 4 :', tt1 - state%start_time
!
!
!          call ComputeSTDGM_Terms(deg_plus )
!          !	print*, 'after STDGMTerms'
!
!          !call cpu_time(tt)
!          !write(*,'(a30, 2es12.4)') 'CT_ends   4 :', tt - state%start_time ,  tt - tt1
!          !write(*,'(a6,2i5)') '###', iter
!          !write(14, *) state%time%iter_loc + 1.*iter/30, tt - tt1,  tt - state%start_time
!
!          call FillVectorST(Newton%b1(1:state%nsize), eta )
!
!          Newton%res0 = VectorPrecondNorm( Newton%b1 )
!
!          !Newton%res0 = dot_product(Newton%b1(:), Newton%b1(:) )**0.5
!          Newton%theta  =  Newton%res0  / Newton%res
!          !write(63,*) 'Newton%res0', Newton%res0
!          !write(63,*) 'Newton%res', Newton%res
!          !write(63,*) 'Newton%theta:', Newton%theta
!
!          !if(l==1) print*
!          !write(*,'(a8, 3i5, 3es10.2, a1, 30es12.4)') 'newton:',state%time%iter, Newton%iter, l, &
!          !    Newton%lambda, lam_old,  Newton%lambda - lam_old,'|', &
!          !    Newton%res0, Newton%res, Newton%theta, state%linSolver%residuum
!
!
!          if(Newton%iter > 1) call Newton_output(grid, 55, l, loc_implicitly )
!
!          if(loc_implicitly) call Newton_output(grid, 54, l, loc_implicitly )
!
!          newton_count = newton_count + 1
!          !      print*, '(related to seeking lambda???) newton count = ', newton_count
!
!          ! residuum is not decreasing, reduce lambda
!          if(Newton%theta >= 1 .and. Newton%res0 > 1E-8 ) then
!             !if(Newton%theta > 1 - Newton%lambda/4  .and. Newton%res0 > 1E-8 )then
!             !Newton%lambda = min(mu1, 0.5* Newton%lambda)
!             !Newton%lambda = Newton%lambda / 2.0
!             Newton%lambda = Newton%lambda / 1.5
!          else
!             !Newton%lambda1 = min(1., 2.0*Newton%lambda )
!             Newton%lambda1 = min(1., 5.*Newton%lambda )
!
!             goto 15
!          endif
!
!          ! too small lambda, NOT CONVERGE
!          if(Newton%lambda < 1E-2 .and. state%modelName /= 'porous') goto 15
!          !if(Newton%lambda < 1E-3) goto 15  ! too small lambda, NOT CONVERGE
!
!       enddo   !l=1,10    ! iterations, seeking the optimal damping factor
!
!       !!state%nlSolver%implicitly = .true. ! last Newton iteration was not sucessfull
!15     continue
!
!
!       if(newton_count >= 8) state%nlSolver%implicitly = .true.
!       !if(newton_count >= 4) state%nlSolver%implicitly = .true.
!       ! !!if(Newton%theta > 0.75) state%nlSolver%implicitly = .true.
!       if(Newton%theta > 0.5) state%nlSolver%implicitly = .true.
!
!       !if(Newton%lambda < 0.9) state%nlSolver%implicitly = .true. !lambda is small, C(w) update
!       if(state%linSolver%lin_solver_not_conv > 0) state%nlSolver%implicitly = .true.
!
!
!       !if(state%nlSolver%implicitly ) write(83,'(a8,i5, a14,i5,a8,es12.4,a8,l2 )')&
!       !     'iter = ',state%time%iter, '### Ncount =',newton_count, &
!       !     ',  theta = ', Newton%theta,',  non conv', state%linSolver%lin_solver_not_conv
!
!
!       !					state%isol = state%isol + 1
!       !
!       !					call WriteProgressOutput( 'STE' )
!
!       !     print*, 'Before algebraic residuum criterion'
!
!       !call WriteMatrixLU()
!
!
!
!       ! nonlinear algebraic residuum criterion
!       if(Newton%non_alg_stop == 'aRES') then
!
!          call cpu_time(t2)
!          time_prepare = time_prepare + t2 - t1
!
!          ! print*, 'PerformOneSTDGMstep calling RezidErrorEstimate(true)'
!          !call RezidErrorEstimates( .true., .false. )  !onlyAS =.true.  -- almost the same
!          call RezidErrorEstimates( .false., .false. )
!
!           do i=1, grid%nelem
!              elem => grid%elem(i)
!              if(elem%i == -1) then ! 556
!                 !if(i==1) print*
!                 write(*,'(a2, i5,i3, i2, 2es9.1, a1,  30es13.5))') &
!                      'N:', elem%i,Newton%iter, l, &
!                      elem%eta(resA:resS,1), '|', &
!                      elem%wST(1,1:elem%dof,1) ,elem%wST(1,1:elem%dof,2)
!              endif
!           enddo
!
!
!          call cpu_time(t1)
!          time_estim = time_estim + t1 - t2
!          !write(*,'(a6, 3es12.4)') 'timesD',time_prepare, time_solve , time_estim
!
!          if(state%time%iter == -10) then
!             do i=1,grid%nelem
!                elem => grid%elem(i)
!                write(100+Newton%iter , *) i, elem%xc, max(1E-15, elem%eta(resA,1)), &
!                     max(1E-15, elem%eta(resS,1)), &
!                     max(1E-15, elem%eta(resA,1)/ max(1E-15, elem%eta(resS,1)) )
!
!                write(200+Newton%iter , *) i, elem%xc, max(1E-15, elem%eta(resA,1))
!                write(200+Newton%iter , *) i, elem%xc, max(1E-15, elem%eta(resS,1))
!                write(200+Newton%iter , '(x)' )
!                write(200+Newton%iter , '(x)' )
!                write(200+Newton%iter , '(x)' )
!
!             enddo
!          endif
!
!          ! if(Newton%iter == 1) &
!          !      write(*,'(a4,a4,5a12,a7,2a12)') '  ', 'iter', '|F(x^k)| ','|F(x^k+1)|', &
!          !      'eta_A', 'eta_S', &
!          !      'rel. err. ', '  LA it',&
!          !      'LA resid', '   CPU(s) '
!
!          ! write(*,'(a4,2i4,5es12.4,a2,i5,2es12.4,l2)') 'NEW:', Newton%iter, l-1, &
!          !       Newton%res,Newton%res0, &
!          !       state%L_estim(resA), state%L_estim(resST), &
!          !       state%L_estim(resA)/ state%L_estim(resST),' |', &
!          !       state%linSolver%iter, state%linSolver%residuum, t1 -  state%start_time,&
!          !       loc_implicitly
!
!          open(91, file='criter', status='unknown', position='append')
!          if(Newton%iter == 1) write(91,'(x)')
!          if(Newton%Aiter == 1) write(91,'(x)')
!          call Newton_output(grid, 91, l, loc_implicitly )
!          close(91)
!
!
!          if(state%estim(resA_ST,1) < Newton%tol2 ) then
!             !if(state%estim(resA_ST_loc,1) < Newton%tol2 ) then
!             call cpu_time(t2)
!             time_prepare = time_prepare + t2 - t1
!             !!write(*,'(a6, 3es12.4)') 'timesE',time_prepare, time_solve , time_estim
!             goto 200  ! Newton method has achieved the prescribed tolerance
!          endif
!
!       else if (Newton%non_alg_stop == 'aDWR') then
!          print*, 'Compute aDWR estimate (elseif in Perform)'
!          stop
!
!       else ! Newton%non_alg_stop == 'rezL2'
!          ! other techniques
!
!          if(Newton%res0 < Newton%tol * state%nsize )  then
!             !write(*,*) 'Reziduum:' , Newton%res0,  '<' , Newton%tol * state%nsize
!             call cpu_time(t2)
!             time_prepare = time_prepare + t2 - t1
!             goto 200  ! Newton method has achieved the prescribed tolerance
!          endif
!
!          !write(*,*) 'Reziduum:' , Newton%res0,  '>' , Newton%tol * state%nsize
!
!       endif  ! end of other technique
!       !  endif  ! end of semi-implicit / implicit
!
!
!       ! go to the next Newton iteration
!       !write(63,*) '####################################################################'
!
!       !write(*,*) 'Going to the' , iter+1 , '-th Newton iteration'
!       !stop 'FR stopped in compute.f90'
!
!       Newton%res1 = Newton%res
!
!       call cpu_time(t2)
!       time_prepare = time_prepare + t2 - t1
!       !write(*,'(a6,4es14.6)')'AAA3', time_prepare, time_solve, t1, t2
!       !!write(*,'(a6, 3es12.4)') 'timesF',time_prepare, time_solve , time_estim
!
!    enddo  ! end of Newton iterations: Newton%iter = 1, Newton%max_iter
!
!
!200 continue
!
!    do i = 1, grid%nelem
!       elem => grid%elem(i)
!
!       call Transfer_wST_to_w_Elem( elem , 0, elem%TQnum)
!    enddo
!
!    !call Newton_output(55, iter+1, loc_implicitly )
!    !if(loc_implicitly) call Newton_output(54, iter+1, loc_implicitly )
!
!
!    if( state%space%estim_space == 'RES' .or. state%time%estim_time == 'tRES') then
!       ! P++       print*,' estimates of the time and space-time error already done '
!       !call RezidErrorEstimates( .false., .false. )
!       call cpu_time(t1)
!       time_estim = time_estim + t1 - t2
!    endif
!
!    loc_implicitly = state%nlSolver%implicitly
!    state%nlSolver%implicitly  = .false.
!
!    !call cpu_time(tt1)
!    !write(*,'(a30, 2es12.4)') 'CT_starts 5 :', tt1 - state%start_time
!
!    ! P++  print*,' not necessary'
!    !call ComputeSTDGM_Terms(deg_plus )
!
!    !call cpu_time(tt)
!    !write(*,'(a30, 2es12.4)') 'CT_ends   5 :', tt - state%start_time ,  tt - tt1
!    !write(*,'(a6,2i5)') '###', iter
!    !write(15, *) state%time%iter_loc + 1.*iter/30, tt - tt1,  tt - state%start_time
!
!
!    Newton%norm_res = EvalSSresidExplicit( )
!    state%nlSolver%implicitly = loc_implicitly
!
!    write(54,*) ' #### Total number of newton updates = ', Newton%updates
!    write(55,*) ' '
!    close(54)
!    close(55)
!
!
!    deallocate(cDLM)
!
!    !if(mod(state%time%iter, 3) == 0) then
!    !   state%isol = state%isol + 1
!    !   call WriteProgressOutput( 'A' )
!    !endif
!
!    state%CPU_prepare = state%CPU_prepare  + time_prepare
!    state%CPU_solve   = state%CPU_solve    + time_solve
!    state%CPU_estim   = state%CPU_estim    + time_estim
!    !!write(*,'(a6, 3es12.4)') 'timXXX',state%CPU_prepare, state%CPU_solve , state%CPU_estim
!
!    end associate
!
!    !print*, 'End PerformOneSTDGMstep'
!
!  end subroutine PerformOneSTDGMstep_OLD
!

  !> Passing to the new time step, after a sucessfull time step,
  !> we compute various indications
  subroutine PassToNewTimeStep( )
    real :: norm, errL2, norm8, err8,t
    real :: errH1, normL2, normH1 , temp
    integer :: i

    state%time%iter = state%time%iter + 1
    state%time%ttime = state%time%ttime + state%time%tau(1)
    state%timeprn = state%timeprn + state%time%tau(1)

    state%time%ctime =  state%time%ttime  ! for the evaluation of the errors

    ! estimate of the time error, state%err(Terr_loc) evaluated in ProposeNewTimeStep in time.f90
    state%err(Terr) = state%err(Terr) +  state%err(Terr_loc)  ! NOT  * state%time%tau(1), Terr_loc = O(\tau^{n+1})


    if (state%time%disc_time == 'STDG') then
       ! VD the following was moved to subroutine UpdateElementW in time.f90
       !do i = 1, grid%nelem
       !   call Eval_wSTfin_Elem (grid%elem(i))

       !   ! call Transfer_wST_to_w_Elem(grid%elem(i), 0, 1)
       !   ! w(rhs,1:grid%elem(i)%dof) * V_rule%phi
       !enddo
       !  print*, 'Eval_wSTfin_Elem in PassToNewTimeStep'

       ! convergence to the steady state solution (if exists)
       call ComputeConvErrors(norm, errL2, norm8, err8)

       ! estimate of the time error, state%err(Terr_loc) evaluated in ProposeNewTimeStep in time.f90
       !state%err(Terr) = state%err(Terr) +  state%err(Terr_loc)  ! NOT  * state%time%tau(1), Terr_loc = O(\tau^{n+1})

       ! new version, based on SS residuum
       if(state%time%iter == 1 .or. state%err(err_0) == 0) &
            state%err(err_0) = max(1E-15, state%nlSolver%norm_res)

       state%err(SSnew) = state%nlSolver%norm_res/state%err(err_0)

       if ( abs(norm) > 0 ) then
          state%err(SSL8) = (errL2/norm)**0.5
       else
          print*, 'NORM is equal to zero in PassToNewTimeStep (cannot divide L2)' ,  errL2,norm
          state%err(SSL8) = 1.0
       endif

       ! former stopping criterium, now in state%err(SSL8)
       !!!state%err(SSL2) = (errL2/norm)**(0.5)/state%time%tau(1)

       ! convergence to the steady state solution in conservative variables
       if(state%modelName == 'NSe' ) then
          !if(state%time%iter <= 2) print*, 'NS-equations not implemented for STDGM'
          !stop
          call DirectComputeDrag_Lift(state%cDLM(state%time%iter, 1:5 )  )
          call ComputeEoFC()
       endif


       if( state%modelName == 'NSe' .and.  (state%type_IC .eq. 6 .or. state%type_IC .eq. 7) ) then
          !print*, 'NS-equations not implemented for STDGM'
          !stop
          state%err(L2_old) = state%err(L2)
          state%err(H1_old) = state%err(H1)

          !call ComputeL2H1ErrorOLD(state%err(L2), state%err(H1), normL2, normH1)
          call ComputeL2H1Error(state%err(L2), state%err(H1), state%err(H1_discrete), normL2, normH1)
          !FILIP this may be use instead of ComputeL2H1Error
!          state%errSnorm(:,:) = 0.
!          !i = 1
!          t = 1.
!          call ComputeL2H1ErrorST(t,state%errSnorm(1,1), state%errSnorm(2,1),&
!                  normL2, normH1)

          call SpaceTimeErrors( )

       elseif( state%modelName == 'scalar' .or. state%modelName == '2eqs' &
            .or. state%modelName == 'porous') then

          state%err(L2_old) = state%err(L2)
          state%err(H1_old) = state%err(H1)


          !FILIP New
          state%errSnorm(:,:) = 0.
          !i = 1
          t = 1.
          call ComputeL2H1ErrorST(t,state%errSnorm(1,1), state%errSnorm(2,1),&
                  normL2, normH1)
          !call ComputeL2H1Error(state%err(L2), state%err(H1), state%err(H1_discrete), normL2, normH1)

          state%err(L2) = state%errSnorm(1,1)
          state%err(H1) = state%errSnorm(2,1)

          ! Warning - errSnorm(1,i) and errSnorm(2,i) are probably not in the same time moment
          !i = 2 - time quadrature nodes
          !wont work if R_rule replaced by T-rule
          ! used to compute the error in Radau nodes -> superconvergence
          do i = 1, state%time%max_Tdof - 1
            t = state%time%T_rule(state%time%max_Tdof)%lambda(i)
            call ComputeL2H1ErrorST(t, errL2, errH1, normL2, normH1)
            state%errSnorm(1,2) =max(state%errSnorm(1,2), errL2)
            state%errSnorm(2,2) =max(state%errSnorm(2,2), errH1)

          enddo

          !i = 3 - arbitrary nodes of the time interval
          do i = 1, 11
            t = (i - 1.) / 10.
            call ComputeL2H1ErrorST(t, errL2, errH1, normL2, normH1)
            state%errSnorm(1,3) =max(state%errSnorm(1,3), errL2)
            state%errSnorm(2,3) =max(state%errSnorm(2,3), errH1)

          enddo

          call SpaceTimeErrors( )

       end if

    ! not STDGM
    else

       ! convergence to the steady state solution (if exists)
       call ComputeConvErrors(norm, errL2, norm8, err8)

       ! stopping criterium for ADIGMA
       !if(state%time%iter == 1) state%err(err_0) = max(1E-15, errL2**(0.5)/state%time%tau(1) )
       !new_err = errL2**(0.5)/state%time%tau(1)
       !state%err(SSnew) = new_err/state%err(err_0)

       ! new version, based on SS residuum
       if(state%time%iter == 1 .or. state%err(err_0) == 0) &
            state%err(err_0) = max(1E-15, state%nlSolver%norm_res)
       state%err(SSnew) = state%nlSolver%norm_res/state%err(err_0)

       !print*
       !write(*,'(a8,i5,20es12.4)') '#@@@!',state%time%iter, state%err(SSnew), state%nlSolver%norm_res, &
       !     state%err(err_0)
       !print*

       state%err(SSL8) = (errL2/norm)**0.5

       ! former stopping criterium , now in state%err(SSL8)
       !!!!state%err(SSL2) = (errL2/norm)**(0.5)/state%time%tau(1)



       ! convergence to the steady state solution in conservative variables
       if(ndim >= 4 .and. ndim <= 6 .and. nbDim == 2) then
          call DirectComputeDrag_Lift(state%cDLM(state%time%iter, 1:5 )  )

          call ComputeEoFC()
       endif


       if( ndim == 4  .and.  (state%type_IC .eq. 6 .or. state%type_IC .eq. 7) ) then
          state%err(L2_old) = state%err(L2)
          state%err(H1_old) = state%err(H1)

          !call ComputeL2H1ErrorOLD(state%err(L2), state%err(H1), normL2, normH1)
          call ComputeL2H1Error(state%err(L2), state%err(H1), state%err(H1_discrete), normL2, normH1)
          !write(*, '(a10,10es16.8)') 'ERRORS:', errL2, errH1, normL2, normH1


          call SpaceTimeErrors( )

          !write(81,'(i5,20es12.4)') state%time%iter,state%err(L2), state%err(H1), &
          !     state%errSTnorm(L2H1)**0.5, state%errSTnorm(L2L2eH1)**0.5

       endif

       if( state%modelName == 'scalar' .or. state%modelName == '2eqs' &
            .or. state%modelName == 'porous') then
          state%err(L2_old) = state%err(L2)
          state%err(H1_old) = state%err(H1)

          call ComputeL2H1Error(state%err(L2), state%err(H1), state%err(H1_discrete), normL2, normH1)

          !write(*,'(a8,10es12.4)') ,'EDERSW',state%err(L2), state%err(H1)

          ! if piecewise constant approximation in time then
          !state%time%ctime = state%time%ttime - state%time%tau(1)
          !call ComputeL2H1Error(state%err(L2_old), state%err(H1_old), normL2, normH1)
          !state%time%ctime = state%time%ttime

          if(state%space%adapt%adapt_method == 'RTN') then
             state%errSTnorm(L2H1) = state%errSTnorm(L2H1) + state%time%tau(1) &
                  *(state%err(H1_old)**2 + state%err(H1)**2 + state%err(H1_old)*state%err(H1))/3

             !call ComputeH1ErrorTimeDer( ) ! error over (0,T) is accumulated in state%errSTnorm(L2L2eH1)
             call ComputeH1ErrorTimeDerDual(1, errH1 ) ! error over (t_{k-1}, t_k)
             !call ComputeH1ErrorTimeDerDualMore(1, errH1 ) ! error over (t_{k-1}, t_k)
             state%errSTnorm(L2L2eH1) = state%errSTnorm(L2L2eH1) + errH1

          elseif(state%space%adapt%adapt_method == 'DUA') then

             state%errSTnorm(L2H1) = state%errSTnorm(L2H1) + state%time%tau(1) &
                  *(state%err(H1_old)**2 + state%err(H1)**2 + state%err(H1_old)*state%err(H1))/3


             call ComputeFluxError( ) ! error over (t_{k-1}, t_k) is stored in state%estim

             call SetRHSDualSTerror( )

             state%errSTnorm(L2L2eH1) = sum( state%estim(eN1p : eN3p, 1) )
             state%errSTnorm(L2F) = sum( state%estim(eN1  : eN2 , 1) )
             state%errSTnorm(NC) = state%estim( NC1n, 1)

             !write(100,'(6es14.6)')  &
             !     state%time%ttime,(2*state%errSTnorm(L2L2eH1))**0.5, state%errSTnorm(L2F)**0.5, &
             !     state%errSTnorm(NC)**0.5

!          elseif (state%space%adapt%adapt_method == 'ALG') then
!
!             if (stop_crit == 'N') then
!                state%errSTnorm(L2H1) = state%err(H1_old)**2   ! CHANGED!!! elem%w(1, :) is considered as a comput sol
!             else
!                state%errSTnorm(L2H1) = state%err(H1_old)**2   ! 'L' or 'G' stop crit, computational solution was obtained at the previous step
!             endif
!
!          elseif (state%space%adapt%adapt_method == 'ALG2') then
!
!             !if (stop_crit == 'N') then
!             !   state%errSTnorm(L2H1) = state%err(H1_old)**2
!             !else
!             if (.not. state%first_GMRES_conv) then
!                state%errSTnorm(L2H1) = state%err(H1_old)**2   ! 'L' or 'G' stop crit, computational solution was obtained at the previous step
!             endif

             !endif

          else !if (.not. state%time%disc_time == 'STDG') then

             !state%errSTnorm(L2H1) = state%errSTnorm(L2H1) + state%time%tau(1) &
             !     *(state%err(H1_old)**2 + state%err(H1)**2 + state%err(H1_old)*state%err(H1))/3
             !
             !state%errSTnorm(L2L2eH1) = state%errSTnorm(L2L2eH1) + state%time%tau(1) &
             !     *( ( state%err(L2_old)**2 + state%model%Re1 * state%err(H1_old)**2 ) &
             !     +  ( state%err(L2)**2 + state%model%Re1 * state%err(H1)**2 ) &
             !     +  ( state%err(L2)**2 + state%model%Re1 * state%err(H1)**2 )**0.5  &
             !     * ( state%err(L2_old)**2 + state%model%Re1 * state%err(H1_old)**2 )**0.5 )/3

             call SpaceTimeErrors( )
          endif
       end if

    endif !not STDGM

  end subroutine PassToNewTimeStep

  ! used to be a part of WriteOutputScreen
  !> control various criteria to end the computation process
  function finishComputationProcess() result (finish)
    logical :: finish
    character(len=3) :: yes_conv ! convergence?

     ! aerodynamical coefficients were already converged?
    if(ndim > 1) then
       yes_conv = '...'
       if(state%EcDLM(1) <= state%EcD_tol)  yes_conv(1:1) = 'D'
       if(state%EcDLM(2) <= state%EcL_tol)  yes_conv(2:2) = 'L'
       if(state%EcDLM(3) <= state%EcM_tol)  yes_conv(3:3) = 'M'
    endif

     ! end of the computation ??
    finish = .false.
    if(state%time%iter_loc == state%time%maxiter) finish = .true.


    ! steady state stopping criterion
    if(ndim > 1) then
       if((state%err(SSnew) <=  state%conv_rez &
            .and.  yes_conv(1:3) =='DLM')  &
            .or. state%time%ttime >= 0.9999999 * state%time%FinTime) finish = .true.

    else
       if(state%err(SSL8) <= state%conv_rez &
            .or. state%time%ttime >= 0.9999999 * state%time%FinTime) &
            finish = .true.

!       if ( state%space%estim_space == 'DWR' ) &
!        print*, 'state%err', state%err(SSL8), 'conv_rez:' , state%conv_rez

       if( state%nlSolver%norm_res > 1E+15) then
          print*,'Too high error, computation stopped'
          !call SetFileNames(.false., command_name, tri_name, sol_name)
          !call OutputDGFEMsol(sol_name)
          stop
       endif

    endif  ! if ndim



  end function finishComputationProcess


  !> write outputs to the screen
  ! finishComputationProcess should be called after this subroutine - they used to be together
  subroutine WriteOutputScreen(iter, time_prepare, time_solve, time_estim)
    integer, intent(in) :: iter
    real, intent(in) ::  time_prepare, time_solve, time_estim
!    logical, intent(inout) :: finish
    integer :: ifile = 11
    integer :: it_hours, it_min, it_sec, i, is
    character(len=1) :: ch
    character(len=3) :: yes_conv ! convergence?
    logical :: Aname
    real:: t_end, ratio

    ! aerodynamical coefficients were already converged?
    if(ndim > 1) then
       yes_conv = '...'
       if(state%EcDLM(1) <= state%EcD_tol)  yes_conv(1:1) = 'D'
       if(state%EcDLM(2) <= state%EcL_tol)  yes_conv(2:2) = 'L'
       if(state%EcDLM(3) <= state%EcM_tol)  yes_conv(3:3) = 'M'
    endif

    ! output of files tri*xxxxx and/or sol*xxxxx
    ch = ' '

    ! no adaptation
    if( state%space%adapt%max_adapt_level == 0 ) then
       if(state%time%ttime >= 0.999999*state%time%FinTime  ) ch = '*'  ! Final time ?
       if( iter >= state%time%maxiter ) ch = '*'      ! Last iteration?
    endif


    ! screen output
    if(  (state%time%time_method == 'E' .and. mod(iter,200) == 1) .or. &
         (state%time%time_method /= 'E' .and. mod(iter, 20) == 1) )then
       if( state%modelName == 'scalar' .or. state%modelName == '2eqs' .or. state%modelName == 'porous' ) then
          print*,'iter (niter)    tau       time     Ax=b&
               & res(iter)   ||res||  ||u-u_h||   diff(u_h)'
          print*,'---------------------------------------------------&
               &------------------------------'
       else
          print*,'iter (niter)    tau       time     Ax=b&
               & res(iter)   || res ||  l_CFL  %CPU cDLM'
          print*,'---------------------------------------------------&
               &-----------------------------'
       endif
    endif

    if(iter == 1 .or. ch == '*' .or. &
         (state%time%time_method == 'E' .and. mod(iter,10) == 0) .or. &
         (state%time%time_method /= 'E' .and. mod(iter, 1) == 0)  )  then
       if( state%modelName == 'scalar' .or.state%modelName == '2eqs'.or.state%modelName == 'porous') then
          write(*,'(i5,a1,i6,a1,es10.2,es11.3,a1,es9.2,a1,i4,a1,i2,es9.2, es12.5,&
               & es9.2)') iter,'(',state%time%iter,')',state%time%tau(1),&
               state%time%ttime,ch, state%linSolver%residuum,'(',state%linSolver%iter,')', &
               state%nlSolver%iter, state%err(SSnew), state%err(L2),state%err(SSL8)
       else
          write(*,'(i5,a1,i6,a1,es10.2,es11.3,a1,es9.2,a1,i4,a1, i2, es10.3, es9.2,&
               & a1, f4.2, a4)') iter,'(',state%time%iter,')',state%time%tau(1),&
               state%time%ttime,ch, state%linSolver%residuum,'(',state%linSolver%iter,')', &
               state%nlSolver%iter,   state%err(SSnew), &
               state%max_eigenvals*state%time%tau(1), ' ', &
               time_solve/(time_solve+time_prepare+1E-8), &
               yes_conv
       endif
    endif



  end subroutine WriteOutputScreen


  !> write outputs to the convfile
  subroutine WriteOutputFiles(convfile, iter, t_sta, time_prepare, time_solve, time_estim)
    character(len=*), intent(in) :: convfile
    integer, intent(in) :: iter
    real, intent(inout) :: t_sta
    real, intent(in) ::  time_prepare, time_solve, time_estim
    integer :: ifile = 11
    integer :: it_hours, it_min, it_sec, i, is
    character(len=1) :: ch
    character(len=3) :: yes_conv ! convergence?
    logical :: Aname
    real:: t_end, ratio


    call cpu_time(t_end)


    ! history of convergence to the file convfile
    open(ifile, file=convfile, status='OLD', position = 'append')
    write(ifile,'(i6,16es14.6, i5,i8, es14.6, i5,  5es12.4, i6, 25es14.6, i8, 6es16.6)') &
         state%time%iter, state%err(SSL8)/state%time%tau(1), state%time%ttime, state%time%tau(1), &   ! 1..4
         state%max_eigenvals*state%time%tau(1), state%nlSolver%norm_res, state%err(SSL8),  & ! 5..7
         state%cDLM(state%time%iter, 1:5 ),  state%EcDLM(1:3), &                  ! 8..15
         state%linSolver%tol, state%linSolver%residuum,  state%linSolver%iter, state%linSolver%iter_tot, &    !16..19
         state%err(SSnew), state%nlSolver%iter,    & !                            !20..21
         t_end-state%start_time, time_prepare, time_solve,  &                !22..24
         time_prepare/(time_prepare + time_solve + 1E-08), &                 !25
         time_solve/(time_prepare + time_solve+ 1E-08), &                    !26
         state%no_refused_steps, &                                           !27
         state%err(L2), state%err(H1), &                                     !28..29
         state%errSTnorm(L2H1)**0.5,  state%errSTnorm(L2L2eH1)**0.5, &       !30..31
         state%errSTnorm(L2F)**0.5, state%errSTnorm(NC)**0.5, &              !32..33
         state%errSTnorm(L2F)**0.5+ state%errSTnorm(NC)**0.5, &              !34
         state%estim(total,1), abs(state%estim( DFnS: DFnT, 1))**0.5, &           !35..37
         state%L_estim(resA:resST),                   &                      !38..41
         abs(state%estim(resA:resST, 1))**0.5, &                                  !42..45
         state%estim(resA_ST,1),state%estim(resA_ST_loc,1), &                !46..47
         state%estim(min_resT_S,1),state%estim(max_resT_S, 1),&              !48..49
         state%err(Terr_loc), state%err(Terr), state%err(Terr_loc)/(1E-12+state%L_estim(resS)), &!50..52
         state%nlSolver%Aiter, state%CPU_prepare, state%CPU_solve, state%CPU_estim, & ! 53..56
         state%CPU_constaint, state%CPU_estim2, time_estim                           ! 57..59

    close(ifile)

    !print*,' solution at this iteration will be saved?',state%time%OutTime
     if(state%time%OutTime > 0. ) then
10      continue
        is = state%isol+1

        if(state%time%ttime -state%time%tau(1) < state%time%OutTime * is &
             .and. state%time%OutTime * is <= state%time%ttime &
             .and. state%time%OutTime * is <= state%time%FinTime*1.0001) then
           state%isol = is
           ratio = ( state%time%OutTime * is - (state%time%ttime - state%time%tau(1)) ) / state%time%tau(1)

           !print*,'ede4erede',ratio
           call WriteProgressOutput_time( 'TSA', ratio)
           !write(*,'(a8, i8,5es12.4)') '>>..<<',is, ratio,  state%time%OutTime * is , state%time%FinTime*1.0001
           goto 10
        endif
        !stop
     endif




    !  computation takes more than 5 minuts, we save the achieved results
    if(t_end - t_sta > 300 .and. .not. state%time_dependent) then
       call WriteResults('Gsol.bak')
       if(ndim >= 4) call ComputeSkinFriction()

       if(state%time%OutTime <= 0. ) call WriteProgressOutput( 'ST', .false., 1. )

       it_hours = t_end/3600
       it_min = (t_end - it_hours*3600)/60
       it_sec = t_end - it_hours*3600 - it_min*60
       write(*,'(a40,i3,a1,i2,a1,i2,a11)') &
            'Results saved in file "Gsol.bak" after ', &
            it_hours,':',it_min,':',it_sec,' (hh:mm:ss)'
       t_sta = t_end
    endif

  end subroutine WriteOutputFiles

  !> write the errors and thier estimates in the file "order.dat"
  subroutine WriteOutputError( )
    class(element), pointer :: elem
    real :: G_estim, G_err, G_err2, t_end, tol_scale, tau
    integer :: i, ifile


    ifile = 11
    open(ifile, file='order.dat', status='UNKNOWN', position = 'append')
    call cpu_time(t_end)


    state%err(L8) =  maxval(grid%elem(:)%errL8)
    state%err(XX) = (state%err(L2)**2 + state%model%Re1*state%err(H1)**2)**0.5

    if (state%time%disc_time == 'STDG') state%errSTnorm(Snorm1:Snorm3) =  state%errSnorm(1,1:3)**2

    !write(*,'(a6,60es16.6)') '!!!!',state%estim(1:max_eta, 1)
    !write(*,'(a6,60es16.6)') '!!!!',state%estim(1:max_eta, 1)**0.5

    !write(*,'(a6,60es16.6)') '!!!!',state%err(1:max_errS)

!    print*, 'State estim:',state%estim(1:max_eta, 1)

    tau = state%time%tau(1)
    if(state%time%tau_choice == 'fixed') tau = state%time%tau_fixed_size

    if(state%modelName == 'NSe' ) then
       state%err(CFD_cD) = state%cDLM(state%time%iter, 1)
       state%err(CFD_cL) = state%cDLM(state%time%iter, 2)
       state%err(CFD_cM) = state%cDLM(state%time%iter, 3)
    endif


    write(ifile,'(i6,2es14.6,i9, 3i4, 2i7, es11.3, 60es16.6, 2i7, 3es12.4)') &
         grid%nelem,  state%space%h, tau,  state%nsize, &                 ! 1..4
         maxval(grid%elem(:)%deg ), state%time%deg, grid%curved_deg, &    ! 5..7
         state%time%iter, state%space%adapt%adapt_level, state%time%ttime,&! 8..10
         state%err(1:max_errS),                           &               ! 11..30
         max(state%errSTnorm(L8L2), state%errSTloc(L8L2) ),     &         ! 31
         sqrt(state%errSTnorm(2:max_errSTnorm) ),                 &       ! 32..40
         sqrt( abs(state%estim(1:max_eta, 1)) ), &                        ! 41..60
         t_end-state%start_time, state%space%domain_volume, &             ! 61..62
         state%CPU_prepare, state%CPU_solve, state%CPU_estim, &           ! 63..65
         state%CPU_estim2, state%CPU_adapt,  1.*state%num_limits /6.,&    ! 66 .. 69
         state%space%adapt%tol_max, &                                     ! 69
         state%space%C_W, state%linSolver%iter_tot, state%space%m_IPG, &  ! 70..72
         state%space%adapt%tol_max / sqrt(state%time%FinTime), &          ! BDF   73
         state%T_estim(resST)/ sqrt( state%time%ttime - state%time%ttime_save ), &  !BDF  74
         state%space%adapt%tol_max *sqrt(( state%time%tau(1))/state%time%FinTime)   ! STDG 75
    close(ifile)

    if(state%modelName == 'scalar') then
       write(*,'(a8,3(a8,es12.4, a2))') 'Errors:', &
            'L^2-norm',  state%err(1), ',', &
            'H^1-semi',  state%err(2), ',', &
            'J_h jumps',  state%err(18), ','
    endif


  end subroutine WriteOutputError

  ! Called after SolveProblemAD to write multiple outputs
  subroutine WriteOutputsAfterSolveProblem( outputType, Aname )
    character(len=*), intent(in) :: outputType
    logical, intent(in) :: Aname
    integer :: ifile
    real :: tt

   ! ALL - from SolveProblemAD
   ! output the final estimates of the error
!    if (state%space%estim_space == 'DWR') then
!         print*, 'WriteOutputError not done for DWR method!'
!         print*, '-----------------------------------------'
!    else
      call WriteOutputError( )
!    endif

    open (11, file=convfile, status='OLD', position = 'append')
    write(11, '(x)' )
    close(11)

    if( (state%modelName == 'scalar' .or. state%modelName == '2eqs') .and. state%type_IC .eq. 5) then
      !print*,'Output for LEVEL SET METHOD in file "data.vofx"'
      call WriteOutput_LevelSet()
    endif


    ! NO ADAPT - from Compute_NoAdapt after SolveProblemAD
    if ( outputType == 'NoAdapt' ) then

       !if(state%space%estim_space == 'pNeu') call WriteProgressOutput('A', Aname)
       if(ndim <= 2)call WriteProgressOutput( 'SEA', Aname, 1. )
       if(ndim >  2)call WriteProgressOutput( 'SA' , Aname, 1.)
       !if(state%type_IC == 7) call WriteProgressOutput( 'E' )

    ! ADAPT STAT - from Compute_stat
    elseif ( outputType == 'AdaptStat' ) then

!       if (state%space%estim_space == 'DWR') then
!         print*, 'WriteOutputError not done for DWR method!'
!         print*, '-----------------------------------------'
!       else
!         call WriteOutputError( )
!       endif



       if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
         if ( state%space%estim_space == 'DWR') then
            call WriteProgressOutput( 'STEAD' , Aname, 1.)
         else
            call WriteProgressOutput( 'STEA' , Aname, 1.)
         endif
       else
         if ( state%space%estim_space == 'DWR') then
            call WriteProgressOutput( 'STAD', Aname, 1. )
         else
            call WriteProgressOutput( 'STA', Aname, 1. )
         endif
       endif

       ! output of the physical coefficients
       if(state%modelName == 'NSe') then
          call cpu_time(tt)
          ifile = 11
          open(ifile, file='AD_coeffs_loc.dat', status='UNKNOWN', position = 'append')
          write(ifile,'(2i6,i9,  30es12.4)') &
               state%space%adapt%adapt_level, grid%nelem,  state%nsize,  tt-state%start_time, & !1..4
               state%space%adapt%tol_max, state%space%adapt%tol_min, &                          !5..6
               state%cDLM(state%time%iter, 1:5 ), state%err(interLq: interH1), &                !7..14
               sqrt(state%estim(1:4,1))                                                        !15..18
          close(ifile)
       endif

    !ADAPT NONSTAT - not used
    else if( outputType == 'AdaptNonStat' ) then
      !nothing special

        if(state%tri_solA) call WriteProgressOutput( 'STA', Aname, 1. )
    else
       stop 'Do nothing in WriteOutputsAfterSolveProblem. Should be called with NoAdapt or AdaptStat or AdaptNonStat!'
    endif


  end subroutine WriteOutputsAfterSolveProblem




end module compute_oper
