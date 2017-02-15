module dual_problem_mod
!  use compute_oper
  use ama_L2interpol
  use data_mod
  use dwr_mod
  use element_mod
  use euler_problem
  use lin_solvers
  use main_data
  use matrix_oper_int
  use mesh_mod
  use solution_mod

implicit none

   public :: InitDualLinSolver
   public :: DualLinSolver
   public :: PrimalLinSolver
   public :: prepareDualProblem
   public :: SolveDualProblem_stat
   public :: updateProblemPlus
   public :: WriteAdwrErrorsScreen
   public :: WriteOutputScreenDP
   public :: WriteDWRErrorsScreen


contains

  !> solve \f$ C(u(T))^T z = J(\phi) \f$ for nonlinear stationary problem
  subroutine SolveDualProblem_stat(DWR, grid, convfile)
    type( DWR_t ), intent(inout) :: DWR
    class( mesh ), intent(inout) :: grid
    character(len=*), intent(in) :: convfile
    class( element ), pointer :: elem
    integer :: i, j, kk, k, elemDof, dof, degP, dofP, R_type, Qnum, nsize
    real :: res_max_val
!    logical :: precond_update
    character(len=50) :: plot_file =  'plot_sol.gnu' !'../DWR/plot_sol.gnu'
    integer :: iPlot = 39
    real, allocatable, dimension(:,:,:,:) :: Dw ! for the least square reconstruction
    real :: l_norm
    real :: residuum
    integer :: lin_solver_not_conv, restart, nloops

    call state%cpuTime%startEstimTime()

    if ( state%time_dependent ) &
      stop 'DWR method is not implemented for time-dependent problems YET!'

    if (allocated(DWR%b) ) &
      stop 'stop DWR%b already allocated here'

      ! change in future to dual_size add control for nonlinear problems size can be ok
    nsize = state%nsize
    allocate(DWR%b(1:nsize), DWR%b1(1:nsize), &
            DWR%x(1:nsize), DWR%rr(1:nsize) )


    eta = 0.0

    !filling elem%rhsST into global rhs
!    print*, 'RHS has to be prepared for DWR%fillRHSvector'
!    call FillVectorST( DWR%b(1:state%nsize), eta ) ! value of eta is not used, only its sign
    ! DWR%rhs -> DWR%b
    ! call with plus = false , because the dof is increased now for DWR_P
    call DWR%fillRHSvector( grid, .false. )

    !call  FillVector(Newton%b(1:state%nsize), eta )
    DWR%x(1:state%nsize) = 0.

    !print*, 'Nullify the number of lin. iterations before dual problem.'
    DWR%linSolver_iter = 0
    !maybe globaly ?
    DWR%residuum = 0 !state%linSolver%residuum
    lin_solver_not_conv = state%linSolver%lin_solver_not_conv

    ! FR  restart and nloops should be set globally
    !set DWR_PARAMETERS
    restart = 50
    nloops = 25

    call SolveBlockLinearSTDGMProblem_Dual( nsize, eta, DWR%b, DWR%x, &
            DWR%residuum, DWR%linSolver_iter, &
            lin_solver_not_conv, DWR%linSolver_tol, restart, nloops)

    !should be changed
    call bMVprodST(DWR%rr, DWR%x, nsize)    ! rr = Ax
    DWR%rr(:) = DWR%rr(:) - DWR%b(:)

    if( lin_solver_not_conv > 0 ) then
      print*, 'lin solver for DUAL problem has not converge yet '
          ! does not do anything !
          open(53, file='GMRES_DP_failed', status='UNKNOWN', position='append')
          res_max_val = 0.
          do i = 1, grid%nelem
             elem => grid%elem(i)

             kk = elem%ncv - 1
             elemDof = elem%dof * ndim * elem%Tdof
             elem%vec( res_vec, 1:elemDof) =  DWR%rr(kk+1:kk+ elemDof)

             res_max_val = max( res_max_val, &
                  sqrt( dot_product(elem%vec( res_vec, 1:elemDof), elem%vec( res_vec, 1:elemDof) )))
          enddo
          close(53)
    endif

    ! copy the solution vector back to elements WST, W, and wS for the least squares
    call CopyZST_fromLongVector( grid, nsize, DWR%x(1:nsize) )


    ! PLOT solution
    ! splot 'plot_sol.gnu' with lines notitle
    open( iPlot, file = plot_file, action="write", status="replace" )
    do i = 1, grid%nelem
      elem => grid%elem(i)
      ! from elem%w( 0, ...)
      call PlotElemSolution3D(iPlot, elem)
    end do
    print*, 'Dual solution saved to plot_sol.gnu'
    close(iPlot)
      ! number of iter = 1, used only with aDWR
    call WriteOutputScreenDP(DWR, 1)

    call state%cpuTime%addEstimTime()

!    l_norm = Solution_L8norm( grid )
!    print*, ' L8 norm of Z: ', l_norm
!    print*,

  end subroutine SolveDualProblem_stat


  !> solve and estimate the linear primal problem until
  !> DWR%estimLP < tol
  !> only for aDWR method - uses DWR%aDWR class, which may not be initialized in other cases !
  subroutine PrimalLinSolver(DWR, grid, newton, tol)
    type( DWR_t ), intent(inout) :: DWR
    class( mesh ), intent(inout) :: grid
    class( NonlinearSol_t ), intent(inout) :: newton
!    integer, intent(in) :: iter ! this subroutine is called for the iter-th time
    real, intent(in) :: tol ! given tolerance
    class( element ), pointer :: elem
    real :: tt, res_max_val, eta1
    integer :: elemDof, kk, i, nsize
    integer :: iPlot = 68
    integer :: lin_solver_not_conv
    integer :: restart, nloops
    logical :: done
    integer :: it, maxIt

    eta1 = 0
    nsize = state%nsize

    ! maximal number of iterations
    maxIt = 50
    it = 1
    done = .false.

    ! only one loop for the alg iterative solver is done
    ! TODO - better working with linear solver ???
    nloops = 1

    do while (.not. done .and. ( it <= maxIt) )
       !Solve
       call state%cpuTime%startSolveTime()
       call SolveBlockLinearSTDGMProblem( nsize, eta1, newton%b, newton%x, &
            state%linSolver%residuum, DWR%aDWR%iter_lin_primal, &
            state%linSolver%lin_solver_not_conv, DWR%aDWR%restart_primal, nloops)

       call state%cpuTime%addSolveTime()
       call state%cpuTime%startPrepareTime()

       !print*, 'iter in PrimalLinSolver:' , DWR%aDWR%iter_lin_primal, state%linSolver%iter
   !    time_solve = time_solve + t1 - t2
       !should be changed

       call bMVprodST(newton%rr, newton%x, state%nsize)    ! rr = Ax
       newton%rr(:) = Newton%rr(:) - Newton%b(:)
       call state%cpuTime%addPrepareTime()

       ! estimate -> DWR%estimLP
       call computePrimalLinearDWRestimates( DWR, grid, newton%rr)


       ! is the error measured with respect to the target quantity
       ! under the given tolerance ?
       if ( DWR%estimLP < tol ) &
         done = .true.

       it = it + 1
    end do !while

    ! PLOT solution
    ! splot 'plot_sol.gnu' with lines notitle
!    open( iPlot, file = plot_file, action="write", status="replace" )
!    do i = 1, grid%nelem
!      elem => grid%elem(i)
!      ! from elem%w( 0, ...)
!      call PlotElemSolution3D(iPlot, elem)
!    end do
!    print*, 'Primal solution saved to plot_sol.gnu'
!    close(iPlot)

      ! the only use of iter
     !call WriteOutputScreenDP(DWR, iter )

  end subroutine PrimalLinSolver




  !> solve and estimate the linear dual problem until
  !> DWR%estimLD < DWR%aDWR%linTol
  !> for the problem \f$  C(u(T))^T z = J(\phi) \f$
  !> only for aDWR method - uses DWR%aDWR class, which may not be initialized in other cases !
  subroutine DualLinSolver(DWR, grid, iter, tol)
    type( DWR_t ), intent(inout) :: DWR
    class( mesh ), intent(inout) :: grid
    integer, intent(in) :: iter ! this subroutine is called for the iter-th time
    real, intent(in) :: tol
    class( element ), pointer :: elem
    real :: t1, res_max_val, eta1
    integer :: elemDof, kk, i, nsize
    integer :: iPlot = 68
    integer :: lin_solver_not_conv
    integer :: restart, nloops
    logical :: done
    integer :: it, maxIt

    eta1 = 0
    nsize = state%nsize

    ! maximal number of iterations
    maxIt = 50
    it = 1
    done = .false.

    ! only one loop for the alg iterative solver is done
    ! TODO - better working with linear solver ???
    nloops = 1

    do while (.not. done .and. ( it <= maxIt) )
      !Solve
      call SolveBlockLinearSTDGMProblem_Dual( nsize, eta1, DWR%b, DWR%x, &
              DWR%residuum, DWR%aDWR%iter_lin_dual, &
              lin_solver_not_conv, DWR%linSolver_tol, DWR%aDWR%restart_dual, nloops)

       !should be changed
       call bMVprodST_Dual(DWR%rr, DWR%x, nsize)    ! rr = A^Tx

       DWR%rr(:) = DWR%rr(:) - DWR%b(:)

       ! copy the solution vector back to elements WST, W, and wS for the least squares
       call CopyZST_fromLongVector( grid, nsize, DWR%x(1:nsize) )

      !Estimate
      call computeDualLinearDWRestimates(DWR, grid)

      if ( DWR%estimLD < tol ) &
         done = .true.

      it = it + 1
    end do !while

!    if( lin_solver_not_conv > 0 ) then
!      print*, 'lin solver for DUAL problem does not converge'
!      ! does not do anything !
!      open(53, file='GMRES_DP_failed', status='UNKNOWN', position='append')
!      res_max_val = 0.
!      do i = 1, grid%nelem
!         elem => grid%elem(i)
!
!         kk = elem%ncv - 1
!         elemDof = elem%dof * ndim * elem%Tdof
!         elem%vec( res_vec, 1:elemDof) =  DWR%rr(kk+1:kk+ elemDof)
!
!         res_max_val = max( res_max_val, sqrt( dot_product( &
!            elem%vec( res_vec, 1:elemDof), elem%vec( res_vec, 1:elemDof) )) )
!      enddo
!      close(53)
!   endif

    ! PLOT solution
    ! splot 'plot_sol.gnu' with lines notitle
!    open( iPlot, file = plot_file, action="write", status="replace" )
!    do i = 1, grid%nelem
!      elem => grid%elem(i)
!      ! from elem%w( 0, ...)
!      call PlotElemSolution3D(iPlot, elem)
!    end do
!    print*, 'Dual solution saved to plot_sol.gnu'
!    close(iPlot)

      ! the only use of iter
     call WriteOutputScreenDP(DWR, iter )

  end subroutine DualLinSolver


  !> updates the size of arrays and matrix blocks when the polynomial degree globally changes
  !> we have to call computeTerms afterwards to fill the blockST...
  subroutine updateProblemPlus( grid, plus )
    class( mesh ), intent( inout ) :: grid
!    type( DWR_t ), intent( in ) :: DWR
    logical, intent(in) :: plus ! true - plus state%space%plusDeg, false - minus %plusDeg
    class( element ), pointer :: elem
    integer :: i, nelem
    integer :: p_plus

!    print*,
!    print*, 'updateProblemPlus called',  plus
!    print*,

    if (state%time_dependent) then
      stop 'updateProblemPlus time dependent - not implemented'
      ! reinit the time stepping - reverse from (0,T) to (T,0)
      !call state%time%reverseTime()
    endif

    nelem = grid%nelem

    ! increase/decrease the polynomial degree
    if ( plus ) then
         p_plus = state%space%plusDeg
    else
         p_plus = - state%space%plusDeg
    endif

    ! the size of the system changes, also dofs and quadratures
    if ( p_plus /= 0 ) then

      call grid%setNewElementDofs( p_plus )

      ! update the sizes of mesh arrays and matrices
      call UpdateMesh_p( grid )
    endif

  end subroutine updateProblemPlus

  !> prepare the structures for computation of the dual problem
  !> it should be called before every solution process, i.e. after every mesh adaptation
  subroutine PrepareDualProblem( DWR )
    class( DWR_t ), intent(inout) :: DWR
    class(element), pointer :: elem
    integer :: i,j
    integer, pointer :: exp1, exp2
    character(len=20) :: loc_adapt_space
    character(len=50) :: plot_file_primal = 'plot_sol_primal.gnu' !'../DWR/plot_sol_primal.gnu'
    integer :: iPlotPrimal = 41
    logical :: impl_loc

    call state%cpuTime%startEstimTime()

    if (state%time%quadType /= 'Radau') &
      stop 'Radau quadrature must be used - to compute C(u(T))'

    ! When preparing the matrix for
    if ( DWR%J%boundary ) &
      stop 'We have to change the matrix C(u) for boundary target functionals!'

    ! print primal solution
    open( iPlotPrimal, file = plot_file_primal, action="write", status="replace" )
    do i = 1, grid%nelem
      elem => grid%elem(i)
      call Transfer_wST_to_w_Elem( elem , 0, elem%TQnum)
      call PlotElemSolution3D(iPlotPrimal, elem)
    end do
    print*, 'Primal solution saved to plot_sol_primal.gnu'
    close(iPlotPrimal)

    ! do nothing for deg_plus == .false.
    if ( DWR%deg_plus ) then
       ! compute the larger matrix and allocate longer vectors for p_plus
       call updateProblemPlus( grid, .true. )
    endif

    ! we need to compute the matrix at u(T)
    impl_loc = state%nlSolver%implicitly
    state%nlSolver%implicitly = .true.
    ! ComputeSTDGM_Terms( ) call
    call ComputeSTDGM_Terms( .false. )
    ! put back
    state%nlSolver%implicitly = impl_loc

    call DWR%J%findSupp( grid )
    ! RHS must be computed after ComputeSTDGM_Terms( )
    !print*, 'setRHS_new called in PrepareDualProblem'
    !call DWR%J%setRHS( grid )
    call DWR%setRHS_new( grid )

    if ( DWR%J%isupp == 0 ) then
        print*, 'epsilon: ', DWR%J%eps, 'grid%h/2', grid%h / 2.0
       stop 'zero support of target func'
    endif

    if (state%time_dependent) &
      stop 'we need to compute IC for dual problem in PrepareDualProblem!'

    call state%cpuTime%addEstimTime()

   end subroutine PrepareDualProblem

  !> write outputs to the screen
  !> finishComputationProcess should be called after this subroutine - they used to be together
  subroutine WriteOutputScreenDP( DWR, iter )
    class( DWR_t ), intent(in) :: DWR
    integer :: ifile = 11
    integer :: it_hours, it_min, it_sec, i, is
    character(len=1) :: ch
    character(len=3) :: yes_conv ! convergence?
    logical :: Aname
    real:: t_end, ratio
    real :: err1, err2, err3
    integer linIter, nlIter, iter

!    iter = 1
    linIter = DWR%linSolver_iter ! DWR%aDWR%lin_iter_dual ???
    ! dual problem is linear
    nlIter = 0 !state%nlSolver%iter
    err1 = -1.0 !state%err(SSnew)
    err2 = -1.0 !state%err(L2)
    err3 = -1.0 !state%err(SSL8)


    ! output of files tri*xxxxx and/or sol*xxxxx
    ch = ' '

    ! no adaptation
    if( state%space%adapt%max_adapt_level == 0 ) then
       if(state%time%ttime >= 0.999999*state%time%FinTime  ) ch = '*'  ! Final time ?
       if( iter >= state%time%maxiter ) ch = '*'      ! Last iteration?
    endif

    ! screen output
    if( state%modelName == 'scalar' .or.state%modelName == '2eqs' ) then
       print*,'iter    tau       time     Ax=b&
            & res(iter)  nlIter  ||res||  ||u-u_h||   diff(u_h)'
       print*,'---------------------------------------------------&
            &------------------------------'
    else
       stop 'WriteOutputScreenDP not done for nonscalar problems!'
    endif

    if(iter == 1 .or. ch == '*') then
       if( state%modelName == 'scalar' .or.state%modelName == '2eqs' ) then
         write(*,'(i5,es10.2,es11.3,a1,es9.2,a1,i4,a1,i6,a2,es9.2, es9.2,&
               & es9.2)') &
               iter, state%time%tau(1),&
               state%time%ttime,ch, DWR%residuum,'(',linIter,')', &
               nlIter,'  ', err1, err2, err3
       else
         stop 'WriteOutputScreenDP not done for nonscalar problems!'
       endif
    endif

  end subroutine WriteOutputScreenDP

  !> write on screen the table of dual error estimates
  subroutine WriteDWRErrorsScreen( DWR )
   class( DWR_t ), intent(in) :: DWR

    print*,'  J(u_h)    J(u)       Error     eta_S    eta_Sabs     eta_A'
    print*,'----------------------------------------------------------------'
    write(*,'( es9.2, a1, es9.2, a2 , es9.2, a2 ,  es9.2, a2, es9.2, a2 , es9.2 )') &
            DWR%J%Ju , ' ', DWR%J%Ju_exact , '  ', DWR%J%Ju - DWR%J%Ju_exact , '  ', &
            sqrt(state%estim(dwrS, 1:ndim)) , '  ', sqrt(state%estim(dwrS_abs, 1:ndim)), '  ', &
            sqrt( state%estim(dwrA, 1:ndim) )
    print*,'----------------------------------------------------------------'
    print*, 'eta_dS    eta_dSabs     eta_dA'
    print*,'----------------------------------------------------------------'
    write(*,'( es9.2, a2, es9.2, a2 , es9.2)') &
            sqrt( state%estim(dwr_dualS, 1:ndim)) , '  ', &
            sqrt( state%estim(dwr_dualS_abs, 1:ndim) ) , '  ', &
            sqrt( state%estim(dwr_dualA, 1:ndim) )
    print*,'----------------------------------------------------------------'
    print*, 'Ju:' , DWR%J%Ju , ' Ju_exact:', DWR%J%Ju_exact
    print*,'----------------------------------------------------------------'

  end subroutine WriteDWRErrorsScreen

    !> write on screen the table of dual error estimates
  subroutine WriteAdwrErrorsScreen( DWR )
   class( DWR_t ), intent(in) :: DWR

    print*,'  Error     eta_Discr    etaNL     etaA    etaAD'
    print*,'----------------------------------------------------------------'
    write(*,'( es9.2, a1, es9.2, a2 , es9.2, a2 ,  es9.2, a2, es9.2)') &
            DWR%J%Ju - DWR%J%Ju_exact , '  ', &
            DWR%estimDiscr, '  ', DWR%estimNL, '  ', &
            DWR%estimLP, '  ', DWR%estimLD
    print*,'----------------------------------------------------------------'

  end subroutine WriteAdwrErrorsScreen


   !> update the structures for computation of the dual problem
   !> it should be called before lin_problem solve
   !> used only for aDWR method for estimation of the alg errors
  subroutine UpdateDualProblem( DWR )
    class( DWR_t ), intent(inout) :: DWR
    class(element), pointer :: elem
    integer :: i,j
    integer, pointer :: exp1, exp2
    logical :: time_dependent
    character(len=20) :: loc_adapt_space
    character(len=50) :: plot_file_primal = 'plot_sol_primal.gnu' !'../DWR/plot_sol_primal.gnu'
    integer :: iPlotPrimal = 41
    logical :: impl_loc
    logical :: deg_plus

    deg_plus = .false.  !??
!    print*, 'DEG_PLUS IN UPDATEDUALPROBLEM????'

!    print*, 'UpdateDualProblem( DWR )'

    ! do nothing for deg_plus == .false.
    if ( DWR%deg_plus ) then
       ! compute the larger matrix and allocate longer vectors for p_plus
!       print*, 'updateProblemPlus'
       call updateProblemPlus( grid, .true. )
    endif

    ! we need to compute the matrix at u(T)
    impl_loc = state%nlSolver%implicitly
    state%nlSolver%implicitly = .true.
    ! ComputeSTDGM_Terms( ) call:
    ! print*, 'calling compute terms in UpdateDualProblem. Although expensive, necessary for nonl problems.'
    ! FR: calling compute terms in UpdateDualProblem. Although expensive, necessary for nonl problems.
    call ComputeSTDGM_Terms( deg_plus )
    ! put back
    state%nlSolver%implicitly = impl_loc

!    call InitDualLinSolver( DWR, grid )

   end subroutine UpdateDualProblem


end module dual_problem_mod
