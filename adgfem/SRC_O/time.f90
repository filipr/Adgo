!> setting of the size of the time steps
module time_sets
  use main_data
  use error_subs
  use inviscid_fluxes
  use stdgm_mod

  implicit none

  public:: SetTimeStep
  public:: SetAdaptTimeStep

  public:: UpdateElementW
!  public:: PrepareNewTimeStep
  public:: PrepareNewTimeStepAD
  public:: ExtrapolatedEstimate
  public:: LocalErrorEstimate
  public:: ProposeNewTimeStep

contains

  !> setting the size of the time step \f$\tau_k \f$ by ad hoc incerese of
  !> CFL number up to value 'state.CFL'
  subroutine SetTimeStep()
    real ::  local_CFL
    real :: alpha = 5.0
    !real :: alpha = 2.5
    !real :: alpha = 1.5
    !real :: alpha = 1.0
    !real :: alpha = 0.75
    !real :: alpha = 0.25
    !real :: alpha = 0.05
    real :: x

    !print*, 'SETTIMESTEP called!'

    if( state%time%tau_choice == 'cfl' ) then
       ! fixed initiation for state%time%BDFtol == 0.
       local_CFL = state%time%CFL

    elseif( state%time%tau_choice == 'exp' ) then
       ! a priori given increase
       local_CFL = state%time%CFL

       local_CFL = local_CFL - (state%time%CFL-1.0)*exp(-alpha*state%time%ttime)
       !local_CFL = local_CFL - (state%time%CFL-0.5)*exp(-alpha*state%time%ttime)
       !local_CFL = local_CFL - (state%time%CFL-0.2)*exp(-alpha*state%time%ttime)
       !local_CFL = local_CFL - (state%time%CFL-0.1)*exp(-alpha*state%time%ttime)
       !local_CFL = local_CFL - (state%time%CFL-0.05)*exp(-alpha*state%time%ttime)

    elseif( state%time%tau_choice == 'adapt' ) then  !ONLY FEW (1 or 2) iterations
       local_CFL = 1.0
       !local_CFL = min(0.5, state%time%BDFtol)

       if(state%space%adapt%adapt_level < 0) local_CFL = 0.01

       !if(state%time_dependent) local_CFL = 0.5
       if(state%time_dependent) local_CFL = 0.01

       !if(state%modelName == 'scalar' .or.state%modelName == '2eqs')  state%time%tau(1) =  state%time%BDFtol
       !state%time%tau(1) = local_CFL * state%time%tau(1)

    endif


    state%time%tau(1) = local_CFL/state%max_eigenvals

    !print*,'Time step?',state%time%iter, state%time%tau(1)

    ! if(state%time%BDFtol > 0. ) then
    !    local_CFL = 1.
    !    !state%time%tau(1) = 0.5 / state%max_eigenvals
    !    !state%time%tau(1) = 1.0/state%max_eigenvals

    !    !state%time%tau(1) = min(0.1, state%time%BDFtol) / state%max_eigenvals
    !    !state%time%tau(1) = min(1.0, state%time%BDFtol) / state%max_eigenvals

    ! endif

    ! !print*,'Time step?',state%time%iter, state%time%tau(1), state%max_eigenvals

  end subroutine SetTimeStep

  !> adaptive setting the size of the time step \f$\tau_k \f$ by
  !> EABDF method
  subroutine SetAdaptTimeStep(time_deg, Tdeg_size, taus, estL, tau_new, refuseStep, non_conv)
    integer, intent(in) :: time_deg,  Tdeg_size
    integer, intent(in) :: non_conv
    real, dimension(1: Tdeg_size), intent(in) :: taus
    real, intent(inout) :: estL  ! estimate of the local error
    real, intent(out) :: tau_new
    logical, intent(out) :: refuseStep
    real, dimension(:), pointer  :: rat
    !real, dimension(7) :: factorial
    !real :: theta, est, coef, e1, e2, coef1
    real :: est, coef, coef1, facmin, facmax,rmin, rmax, fac, err_aver, min_cfl
    real :: kP, kI, kD, tau_min
    integer :: i

    !print*, ' SetAdaptTimeStep'

    if( state%max_eigenvals <= 0.) state%max_eigenvals = 1. ! NO convection case

    if(state%modelName == 'porous' ) state%max_eigenvals = 10. ! NO convection case

    !factorial(1:7) = (/ 1., 2., 6., 24., 120., 720., 5040. /)
    refuseStep= .false.

    facmin = 0.9
    !if (state%time%disc_time == 'STDG') facmin = 0.98
    if (state%time_dependent .and. state%time%disc_time == 'STDG') facmin = 0.95

    if (state%modelName == 'pedes' ) facmin = 0.9


    facmax = 25.
    if( state%ST_Vc >0 .or.  state%ST_Ec> 0.) facmax = 5.
    if( state%type_IC .eq. 10) facmax = 5.
    !if(state%time_dependent) facmax = 1.25
    if(state%time_dependent) facmax = 5.0  ! Porous media
    if(state%time_dependent) facmax = 2.0  ! Porous media

    rmin = 0.98
    rmax = 1.02

    ! local error estimates
    est =  estL


    ! heurictic approach JCP 2011
    if(state%time%estim_time == 'loc') then
       coef = (state%time%BDFtol/est )**(1./(time_deg+1))
       coef1 = coef
       coef = coef *0.85
    endif

    ! approach based on the RES space and time error estims
    !!if(state%space%adapt%adapt_method == 'RES') then  OLD version
    if(state%time%estim_time == 'tRES') then
       !coef from estL is overwritten !!
       if (state%time%disc_time == 'STDG') then
          !coef = (state%time%BDFtol / state%estim(max_resT_S,1) )**(1./(time_deg+1))
          coef = (state%time%BDFtol * state%L_estim(resS)/state%L_estim(resT))**(1./(time_deg+1))
          coef1 = coef
          coef = coef * 0.9
          !!!if(coef >= rmin .and. coef <= rmax ) coef = 1.
          !if( coef > 1 .and. coef <= rmax ) coef = 1.

          coef = max(coef, 0.1) ! minimal decreasing of the time step
          !!coef = max(coef, 0.25) ! minimal decreasing of the time step
       else
!          print*, 'BDFtol, L_estims, L_estims, timedeg: ' , &
!             state%time%BDFtol , state%L_estim(resS) , state%L_estim(resT), time_deg
         !print*, 'resS, resT ' , resS , resT
          coef = (state%time%BDFtol * state%L_estim(resS)/state%L_estim(resT))**(1./(time_deg+1))
          coef1 = coef
          coef = coef * 0.95
          if(coef >= rmin .and. coef <= rmax ) coef = 1.
       endif
       !write(*,'(a6,8es12.4)') '#@!WSE', state%L_estim(resT), state%L_estim(resS), &
       !     state%estim(min_resT_S, 1) , &
       !     state%L_estim(resT) / state%L_estim(resS), coef1, coef



       !open(22, file='ABDF', status='UNKNOWN', position = 'append')
       !write(22,'(2i5,20es14.6)') state%time%iter, time_deg, &
       !     state%L_estim(resT)/state%L_estim(resS), state%L_estim(resT), &
       !     state%L_estim(resS), coef1, coef, taus(1), taus(1)*coef
       !close(22)

       !coef =  state%time%BDFtol / state%estim(min_resT_S, 1)
       !coef =  state%time%BDFtol / (state%err(Terr_loc)/state%L_estim(resS) )

       !err_aver = state%estim(resS, 1) + state%L_estim(resS)**2 *state%time%tau(1)
       !err_aver = (err_aver / (state%time%ttime  + state%time%tau(1) ))**0.5
       !coef =  state%time%BDFtol / (state%err(Terr_loc)/ err_aver )

       if(state%linSolver%lin_solver_not_conv > 0) coef = min(coef, facmin)
    endif ! RES


    if(coef > facmax ) coef = facmax
    if(coef1 < facmin) then   ! coef1 is correct !!!!
       refuseStep = .true.
       !coef = coef * facmin
       !coef = 0.75
    endif

    !if(coef >= rmin .and. coef <= rmax ) coef = 1.

    !tau_new = min(taus(1) * coef , 1E+12/state%max_eigenvals)
    tau_new = min(taus(1) * coef , 1E+15)

    !print*,'e37ye38hds', tau_new, taus(1) * coef , 1E+12/state%max_eigenvals, state%max_eigenvals

    min_cfl = 1E-02
    if(state%modelName == 'scalar' .or. state%modelName == 'porous') min_cfl = 0.   ! no restriction

    !print*, 'Maximal eigenvalue2:', state%max_eigenvals, tau_new
    if(tau_new  <  min_cfl/state%max_eigenvals) then
       tau_new = min_cfl/state%max_eigenvals
       refuseStep = .false.
    endif

    !print*, 'Maximal eigenvalue2:', state%max_eigenvals, tau_new

    if (state%modelName == 'pedes' ) then
       tau_min = 2E-1
       state%time%tau_new = max (state%time%tau_new, tau_min)
       if(state%time%tau(1) <= tau_min * 1.00000001) refuseStep = .false.
    endif



  end subroutine SetAdaptTimeStep

  !> smooths 'coeff' in the SetAdaptTimeStep
  function SmoothTimeCoeff(x)
    real:: SmoothTimeCoeff
    real:: x, SM
    real:: q1L,  q1R, q2L, q2R,  q3,pi

    pi = 3.14159267

    q1L = 0.4
    q2L = 0.8

    q1R = 0.5
    q2R = 2.

    if(x >= 1. - q1L .and. x <= 1. + q1R) then
       SM = 1.

    elseif(x >= 1. + q1R .and. x <= 1. + q2R) then
       SM = q2R * (sin(pi*0.5*(x - q2R - 1.) /(q2R - q1R) ) + 1) + 1

    elseif(x >= 1. - q2L .and. x <= 1. + q1L) then
       SM = q2L * (sin(pi*0.5*(x + q2L - 1.) /(q2L - q1L) ) - 1) + 1

    else
       SM = x
    endif
    SmoothTimeCoeff = SM

  end function SmoothTimeCoeff

  subroutine PrepareNewTimeStepAD(iter)
    integer, intent(in) :: iter
    integer :: i, j
    real :: tt

    call state%cpuTime%startPrepareTime()

    ! setting of the time degree approximation and stroring hisotroy solution
    if( state%time%disc_time /= 'STDG') then

       if(state%time_dependent) then
          state%time%deg_actual = min(state%time%iter - 0 ,  state%time%deg)
       else
          state%time%deg_actual = min(state%time%iter_loc-0, state%time%deg)
       endif

       state%time%deg_actual = max(1, state%time%deg_actual)

       !print*,'#### Tdeg',state%time%deg_actual,state%time%iter_loc, state%time%iter

       ! storing of old time steps
       do i=1,state%time%deg_actual  !!! EBDF
          j = state%time%deg_actual +1  - i  !! EBDF
          state%time%tau(j+1) = state%time%tau(j)
       enddo

    else
       state%time%deg_actual = state%time%deg
    endif ! not STDGM

    ! storing of w on old time levels

    do i=1,grid%nelem
       !print*, 'UpdateElementW ERR'
	   call UpdateElementW(grid%elem(i) )
    enddo




    !print*,' computing of maximal eigenvalue in the first iteration ',state%time%tau_choice, state%time%keep_tau
    if(iter == 1)  then
       state%max_eigenvals = 0.
       do i=1,grid%nelem
          call InitElementEigenvals(grid%elem(i) )
       enddo
       !print*,'Maximal eigenvalue:', state%max_eigenvals

       ! NO convection case
       if(.not. state%model%convective) state%max_eigenvals = 1.  ! NO convection case
       if(state%modelName == 'porous' ) state%max_eigenvals = 10.
       !print*,'Maximal eigenvalue:', state%max_eigenvals
    endif


    ! setting of the new time step obtained from  the previous iteration
    if(state%time%tau_choice == 'fixed') then
       if(state%space%adapt%max_adapt_level > 0) then
          if(state%space%adapt%adapt_level <  0) then
             state%time%tau(1:state%time%deg+1) = state%time%tau_fixed_size/100
          else
             state%time%tau(1:state%time%deg+1) = state%time%tau_fixed_size
          endif
       endif

    elseif(state%time%tau_choice == 'cfl' .or. state%time%tau_choice == 'exp') then
       call SetTimeStep()

    elseif(state%time%tau_choice == 'adapt') then

       if( (state%type_IC == 0 .or. state%type_IC == -2 .or. state%type_IC == -3) .and.  &
            state%space%adapt%adapt_level == 0 .and.  state%time%iter_loc <= 1 ) then
          ! time step is known from the previous computation from 'G*.rsol'
          state%time%tau(1) =  state%time%tau_old

       elseif((state%time%iter <= 2 .and. state%time%estim_time == 'loc' ) .or. &
            (state%time%iter <= 1 .and. state%time%estim_time == 'tRES' )  ) then
          call SetTimeStep()
          !print*,'BAABBB####@@@', state%time%tau(1), state%time%tau_new

       elseif( state%space%adapt%adapt_level == 0 .and. state%time%iter_loc <= 1) then
          call SetTimeStep()
          !print*,'BBBBBB####@@@', state%time%tau(1), state%time%tau_new
          !  state%time%tau(1) =  5E-03


       elseif( state%time%keep_tau .and. state%time%iter_loc <= 1) then
          ! Multi-time step in time dependent computation was Sucesfull
          !print*,' keep the current time step',state%time%tau_new
          state%time%tau(1) = state%time%tau_new

       elseif( .not. state%time%keep_tau .and. state%time%iter_loc <= 1  ) then
          ! Multi-time step in time dependent computation was NOT Sucesfull
          !state%time%tau(1) = state%time%tau_new / 4
          !state%time%tau(1) = state%time%tau(1) / 2
          state%time%tau(1) = state%time%tau(1) * 0.85
          !print*,'####@@@', state%time%tau(1), state%time%tau_new

       else
          !print*,'state%time%tau_new was given in SetAdaptTimeStep'
          state%time%tau(1) = state%time%tau_new
          if(state%space%adapt%adapt_level < 0) state%time%tau(1) = state%time%tau(1) /10.

       endif
    endif


    !print*, state%time%tau(1), state%time%FinTime

    ! to reach the state%time%FinTime exactly
    if(state%time%ttime + state%time%tau(1) > state%time%FinTime ) &
          state%time%tau(1) = state%time%FinTime - state%time%ttime


    !print*,'  SetTimeStep - computation with the stabilization '
    !if( state%time%iter_loc <= 1 .and. ( state%ST_Vc >0 .or.  state%ST_Ec> 0.) .and. &
    !     state%time%tau_choice /= 'fixed' )  &
    !     state%time%tau(1) = 0.1/state%max_eigenvals

    ! clearing of this value
    state%max_eigenvals = 0.

    call state%cpuTime%addPrepareTime()

  !print*, 'End of PrepareNewTimeStep'
  end subroutine PrepareNewTimeStepAD

  !> restore STDGM solution before the repetition of the refused time step
  subroutine  RestoreSTDGMsolution( )
    class(element), pointer:: elem
    integer :: i, j

    do i=1,grid%nelem
       elem => grid%elem(i)
       elem%wST(1:ndim, 1:elem%dof, 1:elem%Tdof) = elem%wSS(1:ndim, 1:elem%dof, 1:elem%Tdof)
    enddo
  end subroutine RestoreSTDGMsolution

  !> deallocation of the arrays with the stored solution
  subroutine  DeallocateStoredSolution( )
    integer :: i

    do i=1,grid%nelem
       if (allocated( grid%elem(i)%wSS) ) &
         deallocate ( grid%elem(i)%wSS)
    enddo

  end subroutine DeallocateStoredSolution

  !> updating of the element solutions
  subroutine UpdateElementW(elem )
    type(element):: elem
    integer :: i, j

!    if(state%SP) then ! saddle point
!       	do i=1,state%time%deg
!	    j = state%time%deg  - i  !! EBDF
!	    elem%wSP(wV1,j+1, : ) = elem%wSP(wV1,j, : )
!	    elem%wSP(wV2,j+1, : ) = elem%wSP(wV2,j, : )
!	    elem%wSP(wP, j+1, : ) = elem%wSP(wP, j, : )
!         enddo
!    else ! NOT  saddle point
       if( state%time%disc_time /= 'STDG') then
          ! storing of old time steps
          do i=1,state%time%deg_actual  !!! EBDF
             j = state%time%deg_actual +1  - i  !! EBDF
             elem%w(j+1, : ) = elem%w(j, : )
          enddo
       else
          call Eval_wSTfin_Elem (elem )
          ! storing of the actual solution
          allocate(elem%wSS(1:ndim, 1:elem%dof, 1:elem%Tdof) )
          elem%wSS(1:ndim, 1:elem%dof, 1:elem%Tdof) = elem%wST(1:ndim, 1:elem%dof, 1:elem%Tdof)
       endif
       elem%w(1, : ) = elem%w(0, : )
!    endif !SP


  end subroutine UpdateElementW

  subroutine ExtrapolatedEstimate(errL2)
    real, intent(out) :: errL2
    real, dimension(:), allocatable :: wi
    real :: norm, errE, norm8, err8
    integer :: i,j, dimdof

    dimdof = (MaxDegreeImplemented +1)*(MaxDegreeImplemented+2)/2
    dimdof = dimdof * ndim

    allocate(wi(1:dimdof) )

    errL2 = 0.
    norm = 0.
    err8 = 0.
    norm8 = 0.

    do i = 1, grid%nelem
       dimdof = grid%elem(i)%dof * ndim

       !print*,MaxDegreeImplemented**2* ndim, ndim


       ! EABDF extrapolation of wi from old time levels
       wi(1:dimdof) = 0.
       do j=1,state%time%deg+1
          wi(1:dimdof) = wi(1:dimdof)  &
               + state%time%Bextrap(j) * grid%elem(i)%w(j,1:dimdof)
       enddo
       ! FIXME: k cemu je BDF%Bextrap

       call ComputeElementConvErrors(grid%elem(i), &
            grid%elem(i)%w(0,:), wi(:), norm, errE, norm8, err8)

       !if(i == 1) then
       !   print*,'dimdof =',dimdof
       !   write(*,'(a3,10es12.4)') 'wE:',wi(1:8)
       !   write(*,'(a3,8es12.4)') 'w1:',grid%elem(i)%w(1,1:8)
       !   write(*,'(a3,8es12.4)') 'w2:',grid%elem(i)%w(2,1:8)
          !write(*,'(a3,10es12.4)') 'wi:',grid%elem(i)%w(0,1:8)
          !write(*,'(a3,10es12.4)') 'wi:',grid%elem(i)%w(0,11:nbDim0)
          !write(*,'(a3,10es12.4)') 'wE:',wi(11:nbDim0)
          !write(*,'(a3,10es12.4)') 'wi:',grid%elem(i)%w(0,21:30)
          !write(*,'(a3,10es12.4)') 'wE:',wi(21:30)
          !write(*,'(a3,10es12.4)') 'wi:',grid%elem(i)%w(0,31:40)
          !write(*,'(a3,10es12.4)') 'wE:',wi(31:40)
       !endif

       errL2 = errL2 + errE
    enddo

    !print*,'DIFF =',errL2

    deallocate(wi)

  end subroutine ExtrapolatedEstimate

  !> estimate of the local discretization error with respect to the time
  !> using the approximation of the high order derivatives
  subroutine LocalErrorEstimate(errL2, errL8)
    real, intent(out) :: errL2, errL8
    real, dimension(:), allocatable :: wi
    real :: norm2, norm8, e1, e2
    integer :: i,j, dimdof, k

    dimdof = (MaxDegreeImplemented +1)*(MaxDegreeImplemented+2)/2
    dimdof = dimdof * ndim

    allocate(wi(1:dimdof) )

    errL2 = 0.
    errL8 = 0.

    do i = 1, grid%nelem
       dimdof = grid%elem(i)%dof * ndim

       ! evaluation of the local error estimates
       wi(1:dimdof) = 0.
       do j=0,state%time%deg+1
          wi(1:dimdof) = wi(1:dimdof)  &
               + state%time%delta(j) * grid%elem(i)%w(j,1:dimdof)

       enddo

       call ComputeElementConvErrors(grid%elem(i), &
            grid%elem(i)%w(0,:), wi(:), norm2, e1, norm8, e2) ! e1, e2 not used

       errL2 = errL2 + norm2
    enddo

    errL2 = (errL2 ** 0.5)*abs(state%time%gamm)
    errL8 = norm8*abs(state%time%gamm)

    deallocate(wi)

  end subroutine LocalErrorEstimate

  !>  estimate the time discretization error and propose new size of the time step
  subroutine ProposeNewTimeStep(refuseStep )
    logical, intent (inout) :: refuseStep
    real :: estL2, estL8
    integer :: NOT_CONV, Tdeg_actual, size_Tdeg

    Tdeg_actual = state%time%deg_actual
    !!if (state%time%disc_time == 'STDG')   Tdeg_actual =  Tdeg_actual + 1
    !print*,'!! EXTRAPOLATED  BDF, \sum_{l=1}^{state%time%deg_actual} \beta_l wold_l',Tdeg_actual
    !call ExtrapolatedEstimate(estL2)

    !  estL2 = L^2 norm,  estL8 = L^{\infty} norm
    if (state%time%disc_time /= 'STDG') &
       call LocalErrorEstimate(estL2, estL8)

    ! estimate of the time error
    state%err(Terr_loc) = estL2 !* state%time%tau(1)**2


    !print*,'###', state%time%deg_actual, state%time%deg, size(state%time%tau, 1)
    !estL2 = estL8

    size_Tdeg = size(state%time%tau, 1)
    call SetAdaptTimeStep(Tdeg_actual, size_Tdeg, state%time%tau(1:size_Tdeg), &
         estL2, state%time%tau_new, refuseStep, NOT_CONV)


    if(refuseStep) then
       write(*,'(a15,a5,es12.4,a5,es12.4,a2,3es12.4)') &
            'Refused tau:','old=', state%time%tau(1),'new=',state%time%tau_new, &
            ':', state%L_estim(resS), state%L_estim(resT), state%L_estim(resT)/state%L_estim(resS)
       !write(*,'(1x)')

       if(.not. state%time_dependent .and. state%modelName == 'NSe' .and.  &
            state%time%tau(1) < 1E-10 ) then
          ! an exception
          refuseStep = .false.

       else

          state%time%tau(1) = state%time%tau_new
          state%no_refused_steps = state%no_refused_steps + 1
       endif

    endif

  end subroutine ProposeNewTimeStep

end module time_sets
