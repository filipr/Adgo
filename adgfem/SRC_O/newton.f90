module newton_mod

   use dual_problem_mod
   use dwr_mod
   use element_mod
   use estimates
   use euler_problem
   use matrix_oper
   use matrix_oper_int
   use main_data
   use mesh_mod
   use paramets

implicit none

   public :: EstimateNewtonADwr
   public :: EstimateNewtonARes
   public :: Newton_output
   public :: NewtonSolve
   public :: performOneNewtonStdgmStep
   public :: setDampingInNewton

contains

   ! seeking of optimal damping parameter in the Newton method,
   ! this should help the convergence of Newton when the approximation is still rough
   subroutine setDampingInNewton( newton, grid, deg_plus, eta, loc_implicitly, res_max_val, l)
      class( Newton_t ), intent ( inout ) :: newton
      class( mesh ), intent(inout) :: grid
      logical, intent(in) :: deg_plus
      real, intent(in) :: eta
      logical, intent(in) :: loc_implicitly
      real, intent(in) :: res_max_val
      integer, intent(out) :: l ! number of damping iterations
      class(element), pointer :: elem
      integer :: max_l, kk, elemDof, i, k, j


       newton%lambda = 1.0   ! initialization of lambda (= damping factor)
       max_l = 10


       do l=1,max_l   ! iterations, seeking the optimal damping factor
          ! update of the solution
          kk = 0
          do i = 1, grid%nelem
             !!lambda_loc = newton%lambda
             elem => grid%elem(i)
             if (elem%ncv /= kk + 1) then
                print*, 'Problem while copying update to wST (compute.f90)'
                stop
             endif
             elemDof = elem%dof * ndim * elem%Tdof
             elem%vec( res_vec, 1:elemDof) =  newton%rr(kk+1:kk+ elemDof)


             do k = 1, elem%Tdof
                do j = 1, ndim

                   elem%wST(j,1:elem%dof,k) = elem%wST(j,1:elem%dof,k) &
                        + (newton%lambda - newton%lambda_old) &
                        * newton%x(kk + 1 : kk + elem%dof)

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

          newton%lambda_old = newton%lambda

          call ComputeSTDGM_Terms(deg_plus )

          call FillVectorST( newton%b1(1:state%nsize), eta )

          newton%res0 = VectorPrecondNorm( newton%b1 )


          newton%theta  =  newton%res0  / newton%res

          if(newton%iter > 1) call Newton_output(grid, 55, l, loc_implicitly )

          if(loc_implicitly) call Newton_output( grid, 54, l, loc_implicitly )

          newton%newton_count = newton%newton_count + 1

          ! residuum is not decreasing, reduce lambda
          if(newton%theta >= 1 .and. newton%res0 > 1E-8 ) then
             !if(newton%theta > 1 - newton%lambda/4  .and. newton%res0 > 1E-8 )then
             !newton%lambda = min(mu1, 0.5* newton%lambda)
             newton%lambda = 0.5* newton%lambda
          else
             newton%lambda1 = min(1., 1.5*newton%lambda )
             ! residuum decreases, quit iterations
             goto 15
          endif

          if(newton%lambda < 1E-2) goto 15  ! too small lambda, NOT CONVERGE
          !if(newton%lambda < 2E-1) goto 15  ! too small lambda, NOT CONVERGE

       enddo   !l=1,10    ! iterations, seeking the optimal damping factor

15     continue

   end subroutine setDampingInNewton


   !> perform one step of the Newton method, i.e.
   !> solve the linear alg. system given by C(u)d^k = -F(u)
   !> solution of the problem is done by GMRES and its parameters: nloops and restart are set here !
   !> but other methods may be implemented in SolveBlockLinearSTDGMProblem later
   !> the tolerance on the linear algebraic error is saved in state%linSolver%tol
   !> When Newton%non_alg_stop == 'aDWR' then the DWR_t structure is needed to compute the estimates of the alg. error
   subroutine performOneNewtonStdgmStep( newton, grid, imp, iter, &
      deg_plus, newtonDone, res_max_val, loc_implicitly, DWR )
      class( NonlinearSol_t ), intent ( inout ) :: newton ! this should the Newton type
      class( mesh ), intent(inout) :: grid
      integer, intent(in) :: imp
      integer, intent(in) :: iter
      logical, intent(in) :: deg_plus
      logical, intent(out) :: newtonDone
      real, intent(out) :: res_max_val
      logical, intent(out) :: loc_implicitly
      type( DWR_t), intent(inout), optional :: DWR

      logical :: update
      class(element), pointer :: elem
      integer :: max_l, l, kk, elemDof, i, k, j
      integer :: restart, nloops
      logical :: vector_update
      !real :: t1, t2, time_prepare, tt, tt1

      newtonDone = .false.

      call state%cpuTime%startPrepareTime()
      if(imp == 0)  newton%implicitly = .true.

      if ( state%nlSolver%non_alg_stop == 'aDWR') then
         print*, 'There may be problems with implicitly!'
         !print*, 'Probably in the matrix there are values with deg_plus!'
         print*, 'Set implicitly = .true.'
         newton%implicitly = .true.
      endif

      loc_implicitly = newton%implicitly
      update = .false.

      if( newton%implicitly ) then
          update = .true.
          call ComputeSTDGM_Terms(deg_plus )

          newton%newton_count = 0
          newton%updates = Newton%updates + 1

          newton%implicitly = .false.
          call ComputeSTDGM_Terms(deg_plus )
          newton%implicitly = .true.

          vector_update = .true.
          newton%implicitly = .false.

      ! computing with the matrices from previous time step ?
      else

          ! for iter > 1, array b(:) = F(x^k)
          if(iter == 1) then
             call ComputeSTDGM_Terms(deg_plus )
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
          ! TODO - what does eta do ???
          call FillVectorST(newton%b(1:state%nsize), eta ) ! value of eta is not used, only its sign
          newton%res = VectorPrecondNorm(Newton%b)
       endif

       if(iter == 1) newton%res_ref = Newton%res

       ! first exit based on the absolute value of the residuum vector || F(w) ||
       if(iter > 1 .and. newton%res/state%nsize < newton%tol &
                   .and. state%nlSolver%non_alg_stop /= 'aDWR' ) then
       !   print*, 'goto200 - Newton%res/state%nsize < Newton%tol',Newton%res/state%nsize, Newton%tol
          !goto 200
          newtonDone = .true.
       endif

       newton%iter  = newton%iter + 1
       newton%Aiter = newton%Aiter  + 1

       newton%lambda_old = 0.
       newton%x(1:state%nsize) = 0.

       if (.not. state%linSolver%tol_fixed) then
          if(iter == 1) then
             state%linSolver%tol = 0.25
          else
             state%linSolver%tol = min(0.25, abs(newton%res - state%linSolver%residuum) / newton%res1)
             state%linSolver%tol = max(state%linSolver%tol, 1E-6)
          endif
       endif

       call state%cpuTime%addPrepareTime()

       !call WriteMatrixST_Blocks(eta)
       !call test_bMVprodST()
       !call Write_rhsST()

       ! set number of iterations and number of nloops of lin alg. solver
       if (newton%non_alg_stop == 'aDWR') then
         ! set number of iter after which restart is done
         if ( present(DWR) ) then
            restart = DWR%aDWR%restart_primal
            nloops = 1
            DWR%aDWR%iter = DWR%aDWR%iter + 1

            ! set the right tolerance
            state%linSolver%tol = DWR%aDWR%linTol * DWR%aDWR%C_L
            !print*, 'state%linSolver%tol set to : ' , state%linSolver%tol
            ! SOLVE THE SYSTEM
            call PrimalLinSolver( DWR, grid, newton, state%linSolver%tol )

         else
            stop 'performOneNewtonStdgmStep must be called with DWR, if Newton%non_alg_stop == aDWR '
         endif
       else
         ! parameters for the linear (GMRES) solver
         restart = 45 !30
         nloops = 10 !5 ! 10

         if(state%model%varying_time_term) eta = 0. ! time deriv terms already included

         ! SOLVE THE SYSTEM
         call state%cpuTime%startSolveTime()
         call SolveBlockLinearSTDGMProblem(state%nsize, eta, newton%b, newton%x, &
            state%linSolver%residuum, state%linSolver%iter, &
            state%linSolver%lin_solver_not_conv, restart, nloops)

         call bMVprodST(newton%rr, newton%x, state%nsize)    ! rr = Ax
         newton%rr(:) = Newton%rr(:) - Newton%b(:)

         call state%cpuTime%addSolveTime()
       endif

       call state%cpuTime%startPrepareTime()

       ! some manipulation for the fail of the computation
       if(state%linSolver%lin_solver_not_conv > 0) then
          open(53, file='GMRES_failed', status='UNKNOWN', position='append')
          res_max_val = 0.
          do i = 1, grid%nelem
             elem => grid%elem(i)

             kk = elem%ncv - 1
             elemDof = elem%dof * ndim * elem%Tdof
             elem%vec( res_vec, 1:elemDof) =  newton%rr(kk+1:kk+ elemDof)
             res_max_val = max( res_max_val, &
                  sqrt( dot_product(elem%vec( res_vec, 1:elemDof), elem%vec( res_vec, 1:elemDof) )))
          enddo
          ! FR added
          close(53)
       endif

       call state%cpuTime%addPrepareTime()

   end subroutine performOneNewtonStdgmStep

   !> write the outputs of the Newton method into the file 'ifile'
   subroutine Newton_output(grid, ifile, l, loc_implicitly )
    class( mesh ), intent(in) :: grid
    integer, intent(in) :: ifile
    integer, intent(in) :: l
    logical, intent(in) :: loc_implicitly
    real :: t2

    associate ( Newton => state%nlSolver )

      call cpu_time(t2)


      write(ifile,'(3i4,14es13.5, i5,2l3,F10.1,2i7,es12.4)') &
           state%time%iter,Newton%iter,l, &                             ! 1..3
           state%time%iter + 1.*(Newton%iter-1)/(Newton%max_iter+2), &      ! 4
           state%estim(resA_ST,1), &                                       ! 5
           Newton%res/state%nsize,  Newton%res0/state%nsize, &             ! 6..7
           Newton%res,  Newton%res0, &                                     ! 8..9
           state%linSolver%tol,  Newton%lambda, Newton%theta, &            !10 .. 12
           state%L_estim(resA:resST), state%time%tau(1), &                 !13 .. 17
           state%linSolver%iter, loc_implicitly, &                         !18 .. 19
           state%linSolver%precond_update, t2- state%start_time, &         !20 .. 21
           grid%nelem, state%nsize,state%time%ttime+ state%time%tau(1)     !22 .. 24

    end associate

   end subroutine Newton_output

   !> Estimate the nonlinear error in Newton for aRES
   !> it should be called only when Newton%non_alg_stop == 'aRES'
   subroutine EstimateNewtonARes( Newton, time_prepare, time_estim, done)
      class( Newton_t ), intent ( in ) :: Newton
      real, intent(inout) :: time_prepare
      real, intent(inout) :: time_estim
      logical, intent(inout) :: done
      real :: t1, t2

      call cpu_time(t2)
      time_prepare = time_prepare + t2 - t1

      ! print*, 'PerformOneSTDGMstep calling RezidErrorEstimate(true)'
      !call RezidErrorEstimates( .true., .false. )  !onlyAS =.true.  -- almost the same
      call RezidErrorEstimates( .false., .false. )

      call cpu_time(t1)
      time_estim = time_estim + t1 - t2
      !write(*,'(a6, 3es12.4)') 'timesD',time_prepare, time_solve , time_estim

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


      if(state%estim(resA_ST,1) < Newton%tol2 ) then
       !if(state%estim(resA_ST_loc,1) < Newton%tol2 ) then
       call cpu_time(t2)
       time_prepare = time_prepare + t2 - t1
       !!write(*,'(a6, 3es12.4)') 'timesE',time_prepare, time_solve , time_estim

       done = .true.

       !goto 200  ! Newton method has achieved the prescribed tolerance
      endif

   end subroutine EstimateNewtonARes


   !> Estimate the nonlinear error in Newton for aDWR
   !> AND compute the dual problem, then compute nl estimate again
   !> it should be called only when Newton%non_alg_stop == 'aDWR'
   subroutine EstimateNewtonADwr( grid, DWR, iter, newtonDone)
      class( mesh ), intent(inout) :: grid
      type( DWR_t), intent(inout) :: DWR
      integer, intent(in) :: iter
      logical, intent(out) :: newtonDone
      real :: etaD, some_criter, t2, t1, dualTol

      ! a) estimate nlError for the 1st time
      ! if we are the 1st mesh then no ZST is computed
      if ( grid%adapt_level > 0 .or. iter > 1 ) then
         call computeNonlinDWRestimates( DWR, grid )

         if ( DWR%estimNL < state%nlSolver%tol ) then
            newtonDone = .true.
         else
             newtonDone = .false.
         endif

         print*, 'After 1st nonlin est! estNL <? TOL' , DWR%estimNL , state%nlSolver%tol
      else
         newtonDone = .true.
      endif

      if (newtonDone) then
            ! b) Dual problem
            call UpdateDualProblem( DWR )
            ! solve the dual problem with tol = C_A * estimNL
            dualTol = DWR%aDWR%linTol * DWR%aDWR%C_L
            call DualLinSolver( DWR, grid, iter, dualTol)

            !DWR%aDWR%iter_lin_dual =  -10  !DWR%linSolver_iter !DWR%aDWR%iter_lin_dual  +
            !call PlotSolDual( DWR%aDWR%iter_lin_dual )
            call PlotSolDual( iter )

            ! c) estimate nlError for the 2nd time
            call computeNonlinDWRestimates( DWR, grid )

            if ( DWR%estimNL < state%nlSolver%tol ) then
               newtonDone = .true.
            else
               newtonDone = .false.
            endif
            print*, 'After second nonlin est! estNL <? TOL' ,  DWR%estimNL , state%nlSolver%tol
      end if

   end subroutine EstimateNewtonADwr

   !> solve the Newton problem F(w) = 0
   !> with tolerance given in state%nlSolver%tol with respect to one of the
   !> estimating methods: aRES, aDWR, rezL2
   subroutine NewtonSolve( grid, imp, iter, deg_plus, DWR)
      class( mesh ), intent(inout) :: grid
      integer, intent(inout) :: imp ! integer which is used to define implicitly - Should be changed
      integer, intent(inout) :: iter ! in - previous iterations , out - total iterations
      logical, intent(in) :: deg_plus
      type( DWR_t), intent(inout), optional :: DWR
      logical :: newtonDone
      real :: time_prepare, time_estim, t1, t2
      integer :: nDamping
      real :: res_max_val, time, lost
      logical :: loc_implicitly

      newtonDone = .false.

      associate( nlSolver=> state%nlSolver)

      ! MAIN NEWTON CYCLE
      do while ( (.not. newtonDone) .and. (iter <= nlSolver%max_iter) )
          ! Three parts:
          !    1) solve the linear system
          !    2) find the damping
          !    3) estimate the error and decide whether to quit - newtonDone = .true.

          ! 1) SOLVE THE LINEAR SYSTEM
          if(nlSolver%non_alg_stop == 'aDWR') then
            if (present(DWR)) then
               ! inner linear primal problem step
               call performOneNewtonStdgmStep( nlSolver, grid, imp,&
                  iter, deg_plus, newtonDone, res_max_val, loc_implicitly, DWR )
            else
               stop 'performOneNewtonStdgmStep should be called with DWR in its &
                     argument for DWR method'
            endif
          else
            call performOneNewtonStdgmStep( nlSolver, grid, imp, iter, &
                  deg_plus, newtonDone, res_max_val, loc_implicitly )
          endif

          ! 2) seeking of optimal damping parameter
          call state%cpuTime%startSolveTime()
          select type ( nlSolver )
            type is ( Newton_t )
               call setDampingInNewton( nlSolver, grid, deg_plus, &
               eta, loc_implicitly, res_max_val, nDamping)
            class default
               stop 'other type of nlSolver'
          end select

          call state%cpuTime%addSolveTime()
          call state%cpuTime%startEstimTime()

          if(nlSolver%newton_count >= 8) nlSolver%implicitly = .true.
          ! !!if(Newton%theta > 0.75) nlSolver%implicitly = .true.
          if(nlSolver%theta > 0.5) nlSolver%implicitly = .true.
          !if(Newton%lambda < 0.9) nlSolver%implicitly = .true. !lambda is small, C(w) update
          if(state%linSolver%lin_solver_not_conv > 0) nlSolver%implicitly = .true.
          !if(nlSolver%implicitly ) write(83,'(a8,i5, a14,i5,a8,es12.4,a8,l2 )')&
          !     'iter = ',state%time%iter, '### Ncount =',newton_count, &
          !     ',  theta = ', Newton%theta,',  non conv', state%linSolver%lin_solver_not_conv
          !					state%isol = state%isol + 1
          !					call WriteProgressOutput( 'STE' )
          !call WriteMatrixLU()

          ! 3 ) ESTIMATE THE ALGEBRACIC (NONLINEAR) ERROR
          select type ( nlSolver )
            type is ( Newton_t )

             ! nonlinear algebraic residuum criterion
             if(nlSolver%non_alg_stop == 'aRES') then
               call EstimateNewtonARes( nlSolver, time_prepare, time_estim, newtonDone )

             else if (nlSolver%non_alg_stop == 'aDWR') then
               call DWR%J%computeJu(grid)
               call DWR%J%computeJu_exact(grid)
               call DWR%writeNlErrorFile(grid%nelem)

               call EstimateNewtonADwr( grid, DWR, iter, newtonDone)
               !call DWR%writeNlErrorFile(grid%nelem)

             else if (nlSolver%non_alg_stop == 'rezL2') then
                if(nlSolver%res0 < nlSolver%tol * state%nsize )  then
                   !write(*,*) 'Reziduum:' , nlSolver%res0,  '<' , nlSolver%tol * state%nsize
                   !call cpu_time(t2)
                   !time_prepare = time_prepare + t2 - t1
                   ! Newton method has achieved the prescribed tolerance
                   newtonDone = .true.
                endif
                !write(*,*) 'Reziduum:' , nlSolver%res0,  '>' , nlSolver%tol * state%nsize
             else
               stop 'Unknown method for nonlinear algebraic error estimation in NewtonSolve!'
             endif  ! end of other technique
            class default
               stop 'other type of nlSolver'
          end select

          call state%cpuTime%addEstimTime( )

          ! Write Newton output
          open(91, file='criter', status='unknown', position='append')
          if(nlSolver%iter == 1) write(91,'(x)')
          if(nlSolver%Aiter == 1) write(91,'(x)')
          if(state%time%iter == 1) print*, 'TODO: Repair NewtonOutput for other methods than aRES!'
          call Newton_output(grid, 91, nDamping, loc_implicitly )
          close(91)

          nlSolver%res1 = nlSolver%res
          iter = iter + 1
      end do ! while
      ! end of newton step, save the solution, write down and  go on

      if (( iter > nlSolver%max_iter) .and. (.not. newtonDone) ) &
         print*, 'Newton solver ended without reaching the tolerance', &
                 'achieving max. number of iterations:' , nlSolver%max_iter
      end associate ! nlSolver

   end subroutine NewtonSolve

end module newton_mod
