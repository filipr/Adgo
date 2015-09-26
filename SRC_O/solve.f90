!> main iterative loop, solution of the unstationary compressible
!> Navier-Stokes equations

module solve_problem
  !  use gmres_solver
  use matrix_oper
  use main_data
  use error_subs
  use problem_oper
  use inviscid_fluxes
  use project_estimation
  use apost_estimation
  use euler_problem
  use dual_estim
  use alg_estim

  implicit none

  public:: InitSolveProblem
  public :: InitSolveProblemALG2
!  public:: SolveProblem
  public:: SolveCDTDGProblem



contains

  !> initiation of the SolveProblem subroutine
  !> called only once in the beginning of the computation
  subroutine InitSolveProblem(convfile)
    character(len=*), intent(in) :: convfile
    real ::  normL2, normH1
    integer :: i, ifile = 11

    state%space%diam = grid%diam

    state%timeprn = 0.

    state%err(SSnew) = 1.
    state%err(SSL8) = 1.

    state%err(Terr) = 0. !  estimate of the time error in compute.f90

    state%space%adapt%adapt_level = 0

    state%linSolver%consistency = 0.
    state%linSolver%iter_tot = 0

    state%no_refused_steps = 0

    state%time%unchanged_time_steps = 0

    if((state%type_IC == 0  .or. state%type_IC == -2) .and. state%time%iter > 0) then
       !if(ndim > 1) &
       allocate(state%cDLM(1:(state%time%iter+(state%space%adapt%max_adapt_level+1)*state%time%maxiter), 1:5))
       state%cDLM(:,:) = 0.

       call ReadConvFile(convfile)
    else
       open(ifile, file=convfile, status='replace') ! deleting of the old file
       !if(ndim > 1) & ! troubles in output for ndim == 1
       allocate(state%cDLM(1:max(15, state%space%adapt%max_adapt_level+1)*state%time%maxiter, 1:5 ) )
       state%cDLM(:,:) = 0.
       close (ifile)   !  file will be opened/closed  at each time step
    endif
    allocate( state%EcDLM (1:3) )
    state%EcDLM(:) = 0.

    state%itime =  0
    state%CPU_prepare = 0.
    state%CPU_solve = 0.
    state%CPU_estim = 0.
    state%CPU_estim2 = 0.
    state%CPU_constaint = 0.

    ! adaptive choice of time step adaptive BDF
    !if(state%time%BDFtol > 0.) then
       ! EABDF
       !if(state%time%deg > 3) then
       !   print*,'Maximum third order of ABDF is implemented'
       !   stop
       !endif
    !endif

    if(state%time%tdg ) then
       ! allocate fake bdf
       call state%time%AllocateBDF( state%time%deg+1)
       state%time%deg_actual = state%time%deg
    elseif(state%time%disc_time /= 'STDG') then
       call state%time%AllocateBDF(state%time%deg)
    endif

    !if (state%time%disc_time == 'STDG') then
    !   open(63, file = 'stdgm-comp', status = 'replace')
    !   close(63)
    !endif

    state%max_eigenvals = 0.    ! maximal eigenvals*|Gi|/|Ki|

    !currently same for stdgm and other methods
    !if (state%time%disc_time == 'STDG') then
   !    do i=1,grid%nelem    ! setting of elem%w(:,:)
   !      call InitElementW(grid%elem(i) )
   !      call InitElementEigenvals(grid%elem(i) )
   !    enddo
   ! else

    do i=1,grid%nelem    ! setting of elem%w(:,:)
       call InitElementW(grid%elem(i) )
       if(.not. state%SP) then
          call InitElementEigenvals(grid%elem(i) )
       endif
    enddo
   ! endif

    !print*,'Maximal eigenvalue:', state%max_eigenvals

    ! the following array stores the values of time steps tau_{k-l}, l=1,..,Tdeg,
    allocate(state%time%tau(1:state%time%deg +1) )

    ! the following array stores a short history for the smooth time adaptation
    allocate(state%time%rat_tau(1:4) )
    state%time%rat_tau(:) = 1.

    if(state%time%tau_fixed .or. state%type_IC == -2) then
       state%time%tau(1:state%time%deg+1) = state%time%tau_new
    else
       print*, 'InitSolveProblem - is it always right to put tau := 0.1 directly?'
       state%time%tau(1:state%time%deg+1) = 0.1
    endif

    ! arrays for dual error estimates
    allocate(state%estim(1:max_eta, 1:ndim), source = 0.0 )
!    state%estim(:, :)  = 0.
    write(debug,*) 'loc_estim temporarily used for RTNst'
    allocate(state%loc_estim(1:max_eta, 1:ndim), source = 0.0 )
!    state%loc_estim(:, :)  = 0.

    allocate(state%L_estim(1:max_eta), source = 0.0 )  !  previously  1:16  !!!!!!

    allocate(state%T_estim(1:max_eta), source = 0.0 )
!    state%L_estim(:)  = 0.
!    state%T_estim(:)  = 0.

    ! arrays for RTN aposteriori error estimates
    allocate(state%eta(1:ndim))
    allocate(state%etaRn_sum(1:ndim), state%etaDFn_sum(1:ndim))
    allocate(state%etaNC1_sum(1:ndim), state%etaNC2_sum(1:ndim))
    allocate(state%etaOsc_sum(1:ndim))


    state%eta(:) = 0  !initialization
    state%etaRn_sum(:) = 0
    state%etaDFn_sum(:) = 0
    state%etaNC1_sum(:) = 0
    state%etaNC2_sum(:) = 0
    state%etaOsc_sum(:) = 0

    state%err(H1_old) = 0.
    state%err(H1) = 0.
    ! state%errSTnorm(L2H1) = 0.
    ! state%errSTnorm(L2L2eH1) = 0.
    ! state%errSTnorm(L2F) = 0.
    ! state%errSTnorm(NC) = 0.

    ! error in the initial condition
    if( ndim == 4  .and.  state%type_IC .eq. 6) then
       !call ComputeL2H1ErrorOLD(state%err(IC_L2),state%err(IC_H1), normL2, normH1)
    endif

    if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then ! both following subroutines are necessary
       !call ComputeL2H1ErrorOLD(state%err(IC_L2),state%err(IC_H1), normL2, normH1)
       call ComputeL2H1Error(state%err(IC_L2), state%err(IC_H1), state%err(H1_discrete), normL2, normH1)
    endif

    state%err(L2) = state%err(IC_L2)
    state%err(H1) = state%err(IC_H1)

    ! inicialization for DUA, allocation of elem%RTNflux, computation for t_k = 0.
    if( state%space%adapt%adapt_method == 'DUA') call ComputeDualEstim( )

    ! inicialization for ALG, computation with initial conditions (as first step)
    if ( state%space%adapt%adapt_method == 'ALG') then
       state%stop_alg = .true.
       state%stop_lin = .true.
       state%stop_alg_glob = .true.
       state%stop_lin_glob = .true.
       call ComputeAlgEstim( )
    endif

    state%space%adapt%stop_adaptation = 0

  end subroutine InitSolveProblem


  !> initiation of the SolveProblem subroutine
  !> after solving discretization solution u_h
  !> the same as InitSolveProblem except for allocations of arrays
  subroutine InitSolveProblemALG2(convfile)
    character(len=*), intent(in) :: convfile
    real ::  normL2, normH1
    integer :: i, ifile = 11

    state%space%diam = grid%diam

    state%timeprn = 0.

    state%err(SSnew) = 1.
    state%err(SSL8) = 1.

    state%err(Terr) = 0. !  estimate of the time error in compute.f90

    state%space%adapt%adapt_level = 0

    state%linSolver%consistency = 0.
    state%linSolver%iter_tot = 0

    state%no_refused_steps = 0

    state%time%unchanged_time_steps = 0

    if((state%type_IC == 0  .or. state%type_IC == -2) .and. state%time%iter > 0) then
       !if(ndim > 1) &
       !!!allocate(state%cDLM(1:(state%time%iter+(state%space%adapt%max_adapt_level+1)*state%time%maxiter), 1:5))   ! done already before u_h was computed
       call ReadConvFile(convfile)
    else
       open(ifile, file=convfile, status='replace') ! deleting of the old file
       !if(ndim > 1) & ! troubles in output for ndim == 1
       !!!allocate(state%cDLM(1:max(15, state%space%adapt%max_adapt_level+1)*state%time%maxiter, 1:5 ) )   ! done already before u_h was computed
       close (ifile)   !  file will be opened/closed  at each time step
    endif

    state%itime =  0
    state%CPU_prepare = 0.
    state%CPU_solve = 0.
    state%CPU_estim = 0.
    state%CPU_estim2 = 0.
    state%CPU_constaint = 0.

    if(state%time%tdg ) then
       ! allocate fake bdf
       !!!call AllocateABDF(state%BDF, state%time%deg+1)  ! done already before u_h was computed
       state%time%deg_actual = state%time%deg
    else
       !!!call AllocateABDF(state%BDF, state%time%deg)   ! done already before u_h was computed
    endif

    state%max_eigenvals = 0.    ! maximal eigenvals*|Gi|/|Ki|

    do i=1,grid%nelem    ! setting of elem%w(:.:)
       call InitElementW(grid%elem(i) )
       call InitElementEigenvals(grid%elem(i) )
    enddo

    ! the following array stores the values of time steps tau_{k-l}, l=1,..,Tdeg,
    allocate(state%time%tau(1:state%time%deg +1) )      ! done already before u_h was computed

    ! the following array stores a short history for the smooth time adaptation
    allocate(state%time%rat_tau(1:4) )   ! done already before u_h was computed
    state%time%rat_tau(:) = 1.

    if(state%time%tau_fixed .or. state%type_IC == -2) then
       state%time%tau(1:state%time%deg+1) = state%time%tau_new
    else
       state%time%tau(1:state%time%deg+1) = 0.1
    endif

    ! arrays for dual error estimates
    !!!allocate(state%L_estim(1:16) )   ! 16 is given manually !!!!!!   ! done already before u_h was computed

    ! arrays for RTN aposteriori error estimates
    !!!allocate(state%eta(1:ndim))
    !!!allocate(state%etaRn_sum(1:ndim), state%etaDFn_sum(1:ndim))           ! done already before u_h was computed
    !!!allocate(state%etaNC1_sum(1:ndim), state%etaNC2_sum(1:ndim))
    !!!allocate(state%etaOsc_sum(1:ndim))

    state%eta(:) = 0  !initialization
    state%etaRn_sum(:) = 0
    state%etaDFn_sum(:) = 0
    state%etaNC1_sum(:) = 0
    state%etaNC2_sum(:) = 0
    state%etaOsc_sum(:) = 0

    !!!allocate(state%estim(1:max_eta, 1:ndim) )   ! done already before u_h was computed
    state%estim(:, :)  = 0.

    state%err(H1_old) = 0.
    state%err(H1) = 0.
    state%errSTnorm( : ) = 0.

    ! error in the initial condition
    if( ndim == 4  .and.  state%type_IC .eq. 6) then
       !call ComputeL2H1ErrorOLD(state%err(IC_L2),state%err(IC_H1), normL2, normH1)
    endif

    if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then ! both following subroutines are necessary
       !call ComputeL2H1ErrorOLD(state%err(IC_L2),state%err(IC_H1), normL2, normH1)
       call ComputeL2H1Error(state%err(IC_L2), state%err(IC_H1), state%err(H1_discrete), normL2, normH1)
    endif

    state%err(L2) = state%err(IC_L2)
    state%err(H1) = state%err(IC_H1)

    ! inicialization for DUA, allocation of elem%RTNflux, computation for t_k = 0.
    if( state%space%adapt%adapt_method == 'DUA') call ComputeDualEstim( )

    ! inicialization for ALG, computation with initial conditions (as first step)
    if ( state%space%adapt%adapt_method == 'ALG2') then
       state%stop_alg = .true.
       state%stop_lin = .true.
       state%stop_alg_glob = .true.
       state%stop_lin_glob = .true.
       call ComputeAlgEstim( )
    endif

  end subroutine InitSolveProblemALG2


!  ! Solve constant-diffusion time-discontinuous galerkin discrete problem
  subroutine SolveCDTDGProblem(convfile)
    character(len=*), intent(in) :: convfile
    character(len=1) :: ch
    character(len=3) :: yes_conv ! convergence?
    character(len=50) :: command_name, tri_name, sol_name, exa_name, err_name, est_name
    integer :: ifile=84!, icDcL = 81
    real, dimension(:,:), allocatable ::  b, X, residue
    real, allocatable:: X0(:)
    integer :: it, iter, i,j, iflag
    real :: sum1in, sum1, errL2, errH1, logerr, norm
    real :: errL2sup, errH1int, normL2, normH1
    real :: norm_res, norm_res2!    external dgmres;
    real :: ts, tA1s, tA1e, tA2s, tA2e, te,te1, t_sta, t_end
    integer :: it_hours, it_min, it_sec
    real :: diff = 1D-10 ! diference ||w - w'|| computed by 2 different schemes
    real :: tau_new     ! computed new time step for adaptive BDF
    integer :: istep, NOT_CONV1, NOT_CONV2
    logical :: refuseStep
    real, dimension(1:5) :: coeffs   ! cD, cL, cM, cD_p, cL_p
    real, dimension(:,:),allocatable:: KK
    real :: err8_0 ! initial residual
    real :: norm8, err8
    real :: new_err, lognew_err, logdiff
    real :: y, resid
    integer :: k, l, dof, kst, inode
    type(Gauss_rule):: G_rule

    call Create_Gauss_rule(G_rule, state%time%deg+1)
    allocate(KK(state%time%deg+1,state%time%deg+1))

    allocate(b(1:state%nsize,state%time%deg+1) )
    allocate(X(1:state%nsize,state%time%deg+1), X0(1:state%nsize) )
    allocate(residue(1:state%nsize,state%time%deg+1))

    if(state%space%adapt%max_adapt_level == 0 .or. &
         (state%space%adapt%max_adapt_level > 0 .and. state%space%adapt%adapt_level <= 0) )then

       if (.not. state%time%tau_fixed) then
          call SetTimeStep()
          state%time%tau_new = state%time%tau(1)
       else
          if(state%space%adapt%max_adapt_level > 0 .and. state%space%adapt%adapt_level < 0) then
             state%time%tau(1) = state%time%tau_fixed_size/100
             state%time%tau_new = state%time%tau(1)
          else
             state%time%tau(1) = state%time%tau_fixed_size
             state%time%tau_new = state%time%tau(1)
          endif
          print *,'using explicit time step in CDTDG = ',state%time%tau_new
       endif

       state%time%tau_old = state%time%tau_new
       !if( state%space%adapt%max_adapt_level > 0 .and. state%space%adapt%adapt_level < 0)  then
       !   state%time%tau_old = state%time%tau_new/100
       !endif


       state%time%tau = state%time%ttime + state%time%tau_old*(G_rule%lambda - 1)
       do i=1,grid%nelem
          call ExactElementW(grid%elem(i), state%time%tau)
       enddo
       !print*,'!!!!',state%time%tau(:)
    endif


    print*
    print*,'# Iteration process starts'
    print*

    errL2sup = 0.
    errH1int = 0.

    call cpu_time(t_sta)
    do iter = 1,state%time%maxiter
       call cpu_time(ts)

       ch = ' '
       ! no adaptation
       if( state%space%adapt%max_adapt_level == 0 ) then
          if(state%time%ttime + state%time%tau(1) >= state%time%FinTime  ) ch = '*'  ! Final time ?
          if( iter >= state%time%maxiter ) ch = '*'      ! Last iteration?
       endif

       ! solution at this iteration whould be saved?
       if(state%time%OutTime > 0. .and. state%space%adapt%adapt_level >= 0) then
          if(state%timeprn + state%time%tau(1) >= state%time%OutTime)  ch = '*'
       endif

       ! assemble X0
       k = 0
       do i=1,grid%nelem
         j = grid%elem(i)%dof*ndim
         X0(k+1: k+j) = grid%elem(i)%w(0,1:j)
         k = k+j
       enddo

       !write(*,'(a5,i5,21es12.4)') 'X0',inode, state%time%ctime,X0(:)

       call cpu_time(tA1s)

       do inode = 1,state%time%deg + 1
         state%time%tau = state%time%tau_old * (G_rule%lambda - 1) - state%time%tau_new * G_rule%lambda(inode)

         ! seting of fake coefficients \alpha, \beta for ABDF, to actually perform a
         ! lagrange extrapolation
         call state%time%FakeBDFLagrangeExtr( state%time%deg_actual, state%time%tau)
         state%time%ctime = state%time%ttime + state%time%tau_new * G_rule%lambda(inode)

         !print*,'Compute Elements Terms'
         if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then   ! scalar equation
           call ComputeElementsTerms(Set_f_s_scalar, Set_A_s_scalar, Set_Ppm_scalar, &
                Set_R_s_scalar, Set_K_sk_scalar, Set_S_scalar, Set_DS_scalar)
         else
           print*,'SolveCDTDGProblem not implemented for non-scalar equations'
         endif
         !print*,'ComputeElementsTerms -- done'

         call FillVectorsOLD(b(:,inode),  X(:,inode), 0.) ! FIXME je tohle OK?

         !write(*,'(a5,i5,21es12.4)') 'RHS',inode, state%time%ctime,b(:,inode)
         !write(*,'(a5,i5,21es12.4)') 'X',inode, state%time%ctime,X(:,inode)

       end do

       !print *,'extrap = ',sum(X,1)

       call cpu_time(tA1e)

       call cpu_time(tA2s)
       call PrepareTimeElementLinearProblem(state%time%tau_new, G_rule%weights, &
            G_rule%lambda, b, KK, X0)

       state%linSolver%iter = 0
       call SolveMixedSylvesterProblem(X, b, KK, state%linSolver%residuum, state%linSolver%iter, &
            NOT_CONV1)
       call cpu_time(tA2e)

       !write(*,'(a5,3es12.4)') '  B0',grid%elem(2000)%w(0,:)

       !print *,'solved = ',sum(X,1)

       do inode = 1,state%time%deg+1
         k = 0
         do i=1,grid%nelem
           j = grid%elem(i)%dof*ndim
           grid%elem(i)%w(inode,1:j) = X(k+1:k+j,inode)
           grid%elem(i)%w(0,1:j) = grid%elem(i)%w(inode,1:j)
           k = k+j
         end do
         state%time%ctime = state%time%ttime + state%time%tau_new * G_rule%lambda(inode)

         !call ComputeL2H1ErrorOLD(errL2, errH1, normL2, normH1)
         call ComputeL2H1Error(errL2, errH1,  state%err(H1_discrete),normL2, normH1)
         errL2sup = max (errL2sup, errL2)
         errH1int = sqrt(errH1int**2 + state%time%tau_new * G_rule%weights(inode) * errH1**2)

         !write(*,'(a5,3es12.4)') '  B1',grid%elem(2000)%w(0,:)

       end do

       state%time%tau = state%time%tau_new * (G_rule%lambda - 1)

       ! seting of fake coefficients \alpha, \beta for ABDF, to actually perform a
       ! lagrange extrapolation
       state%time%ctime = state%time%ttime + state%time%tau_new
       call state%time%FakeBDFLagrangeExtr( state%time%deg, state%time%tau)

       !print*,'Compute Elements Terms'
       !write(*,'(a5,3es12.4)') '  b2',grid%elem(2000)%w(0,:)

       if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then   ! scalar equation
         call ComputeElementsTerms(Set_f_s_scalar, Set_A_s_scalar, Set_Ppm_scalar,&
              Set_R_s_scalar, Set_K_sk_scalar, Set_S_scalar, Set_DS_scalar)
       else
         print*,'SolveCDTDGProblem not implemented for ndim=',ndim
       endif

       !write(*,'(a5,3es12.4)') '  B2',grid%elem(2000)%w(0,:)


       state%time%iter = state%time%iter + 1
       state%time%ttime = state%time%ttime + state%time%tau_new
       state%timeprn = state%timeprn + state%time%tau_new

       state%time%tau(1) = state%time%tau_new

       ! convergence to the steady state solution (if exists)
       call ComputeConvErrorsTDG(norm, errL2, norm8, err8)

       ! stopping criterium for ADIGMA
       if(state%time%iter == 1) state%err(err_0) = max(1E-15, errL2**(0.5)/state%time%tau(1) )
       new_err = errL2**(0.5)/state%time%tau_new
       state%err(SSnew) = new_err/state%err(err_0)

       !if(state%time%iter == 1) err8_0 = max(1E-15, err8 )
       !state%err(SSL8) = err8/err8_0
       ! NEW VALUE for SSerrL8
       state%err(SSL8) = (errL2/norm)**0.5

       lognew_err = -999.999
       if(new_err > 0) lognew_err = log10(new_err)

       ! former stopping criterium
       errL2 = (errL2/norm)**(0.5)/state%time%tau(1)
       state%err(L2) = errL2

       ! logarithmus of the convergence error
       logerr = -999.999
       if(errL2 > 0) logerr = log10(errL2)

       ! convergence to the steady state solution in conservative variables
       ! call ComputeResConserVar(state%w, state%wold, state%rez_conser)

       if(ndim == 4) then
          call DirectComputeDrag_Lift(coeffs(1:5) )
          state%cDLM(state%time%iter, 1:3 ) = coeffs(1:3)
       endif

       call ComputeEoFC()

       call cpu_time(te)

       !write(*,'(a5,3es12.4)') '  C0',grid%elem(2000)%w(0,:)

       !write(16,'(a6,6es12.4)') 'times:', &
       !     (tA1s -ts )/(te - ts), (tA1e -tA1s)/(te - ts), &
       !     (tA2e -tA1e)/(te - ts), &
       !     (te1 -tA2e )/(te - ts), (te -te1)/(te - ts), te - ts

       !write(17,'(a6,6es12.4)') 'times:', &
       !     (tA1s -ts ), (tA1e -tA1s), &
       !     (tA2e -tA1e), &
       !     (te1 -tA2e ), (te -te1), te - ts

       open(ifile, file=convfile, status='OLD', position = 'append')
       write(ifile,'(i6,16e14.6,i5, 7e14.6)') &
            state%time%iter, logerr, state%time%ttime, state%time%tau_new, &                 ! 1..4
            state%max_eigenvals*state%time%tau_new,  state%err(SSnew), lognew_err, &  ! 5..7
            coeffs(1:5),  state%EcDLM(1:3), &                                ! 8..15
            state%linSolver%tol, state%linSolver%residuum, state%linSolver%iter, &                        ! 16..18
            logdiff, te-ts, tA1e-tA1s, tA2e-tA2s,  &                         ! 19..25
            (tA1e-tA1s)/(te-ts), (tA2e-tA2s)/(te-ts),  (tA1e-tA1s+ tA2e-tA2s)/(te-ts)
       close(ifile)


       !call ComputeErrorRingleb(errL2)
       !write(71,*)state%time%iter, errL2

       ! screen output
       if(mod(iter,20) == 1) then
          print*,'iter (niter)    tau       time     Ax=b&
               & res(iter)    errL2sup    errH1int  conv_error  l_CFL  %CPU cDLM'
!               & res(iter)    errL2sup    errH1int  conv_error  l_CFL  %CPU cDLM'
          print*,'---------------------------------------------------&
               &----------------------------'
          !call ComputeCF_Pk()
       endif

       yes_conv = '...'
       if(state%EcDLM(1) <= state%EcD_tol)  yes_conv(1:1) = 'D'
       if(state%EcDLM(2) <= state%EcL_tol)  yes_conv(2:2) = 'L'
       if(state%EcDLM(3) <= state%EcM_tol)  yes_conv(3:3) = 'M'


       if(iter == 1 .or. mod(iter,1) == 0 .or. ch =='*')  then
          write(*,'(i5,a1,i6,a1,es10.2,es11.3,a1,es9.2,a1,i4,a1, es12.4, es9.2,&
                  & a1, f4.2, a4)') iter,'(',state%time%iter,')',state%time%tau(1),&
                  state%time%ttime,ch, state%linSolver%residuum,'(',state%linSolver%iter,')', state%err(SSnew), &
                  state%max_eigenvals*state%time%tau_new, ' ',&
                  (tA1e-tA1s+ tA2e-tA2s)/(te-ts), yes_conv
       endif

       if(ch == '*') then
          state%isol = state%isol + 1
          call SetFileNames(.false., command_name, tri_name, sol_name, exa_name, err_name, &
               est_name, .false.)
          if(state%space%adapt%max_adapt_level >0) call OutputDGFEMtri(tri_name)
          call OutputDGFEMsol(sol_name)
          state%timeprn = state%timeprn - state%time%OutTime
       endif


       call cpu_time(t_end)
       if(t_end - t_sta > 300) then
          !  computation takes more than 5 minuts, we save the achieved results
          call WriteResults('Gsol.bak')

          it_hours = t_end/3600
          it_min = (t_end - it_hours*3600)/60
          it_sec = t_end - it_hours*3600 - it_min*60
          write(*,'(a40,i3,a1,i2,a1,i2,a11)') &
               'Results saved in file "Gsol.bak" after ', &
               it_hours,':',it_min,':',it_sec,' (hh:mm:ss)'
          t_sta = t_end
       endif

       if((state%err(SSnew) <=  state%conv_rez .and.  yes_conv(1:3) =='DLM')  &
            .or. state%time%ttime >= state%time%FinTime .or. &
            state%err(SSnew) > 1E+15)    goto 100


       !write(*,'(a5,3es12.4)') '  D0',grid%elem(2000)%w(0,:)

    enddo

100 continue

    deallocate(b, X, residue )


    if(state%type_IC .eq. 4 .or. state%modelName == 'scalar' .or.state%modelName == '2eqs' ) then
       state%time%ctime = state%time%ttime
       !call ComputeL2H1ErrorOLD(errL2, errH1, normL2, normH1)
       call ComputeL2H1Error(errL2, errH1,  state%err(H1_discrete),normL2, normH1)
       print*,'Error in L_2 norm = ',errL2
       print*,'Error in H_1 semi-norm = ',errH1

       print*,'Error in sup L_2 norm = ',errL2sup
       print*,'Error in int H_1 semi-norm = ',errH1int

       open(ifile, file='order.dat', status='UNKNOWN', position = 'append')
       write(ifile,'(i6,2es14.6, 3i4, 5es14.6)') &
            grid%nelem, state%space%h, state%time%tau_new,state%space%deg, state%time%deg, &
            grid%curved_deg, errL2, errH1, errL2sup, errH1int, t_end-state%start_time

       close(ifile)
    endif

    print*
    print*,'# Iteration process finished'
    print*

    if( state%err(SSnew) > 1E+15) then
       print*,'Too high error, computation stopped'
       !stop
    endif

  end subroutine SolveCDTDGProblem

end module solve_problem
