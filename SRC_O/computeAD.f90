!> new (more compact)
!> main iterative loop, solution of the unstationary compressible
!> Navier-Stokes equations

module computeAD_oper
  use main_data
  use io_sub
  use solve_problem
  use marking
  use estimates
  use errorDual
  use errorFlux
  use hp_adaptation
  use red_green
  !use pMultiGrid !MAKE_MG

  use helmholtz_estim
  use neumann_estim
  use dual_estim
  use alg_estim
  use st_interpol
  use anisotropic
  use AMA_estims
  use ama_hp_interpol
  use error_subs
  use mesh_oper
  use loc_rav_tho_ned
  use compute_oper
!  use compute_operSP
  use higher_order_estims
  use rtn_st_mod
  use higher_order_local

  implicit none

  public:: ErrorEstims
  public:: ErrorEstims_time_dependent
  public:: ErrorEstims_time_independent
  public:: SolveProblemAD
  public:: RecomputeBackMesh
  public:: AdaptationAD
!  public:: ComputeAD
  public:: SolveDualProblem

contains

  !> several type of error estimates
  subroutine ErrorEstims( finish )
    logical, intent(in) :: finish
    integer :: i
    real :: t1, t2

    call cpu_time(t1)

    !print*, 'Calling ErrorEstims'
    if(state%time_dependent) then
       call ErrorEstims_time_dependent(finish )

    else
       call ErrorEstims_time_independent( )
    end if  !if(state%time_dependent

    call cpu_time(t2)
    state%CPU_estim2 = state%CPU_estim2 + t2 - t1


  end subroutine ErrorEstims

  !> several type of error estimates
  !> TIME-DEPENDENT COMPUTATIONS
  subroutine ErrorEstims_time_dependent(finish )
    logical, intent(in) :: finish
    real ::  r_tol, r_est, r_estOLD, Lq
    integer :: i, tQnum

    r_tol = 0.
    r_est = 0.

    if( state%space%estim_space == 'RES' ) then
       ! STDGM
       if( state%time%disc_time == 'STDG')  then
          do i=resA, resST  ! i.e., resA, resS, resT, resST
             state%estim(i, 1) = state%estim(i, 1) + state%L_estim(i)**2 !*state%time%tau(1)
             state%T_estim(i)  = state%T_estim(i)  + state%L_estim(i)**2 !*state%time%tau(1)
          enddo

          ! estimated error and its tolerance
          r_tol = state%space%adapt%tol_max  &
               * (( state%time%ttime - state%time%ttime_save ) /state%time%FinTime)**0.5
          r_est = state%T_estim(resST)**0.5

       else ! BDF time discretization
          call RezidErrorEstimates( .false.  , .true.)
          do i=resA, resST  ! i.e., resA, resS, resT, resST
             state%estim(i, 1) = state%estim(i, 1) + state%L_estim(i)**2 *state%time%tau(1)
             state%T_estim(i)  = state%T_estim(i)  + state%L_estim(i)**2 *state%time%tau(1)
          enddo

          ! estimated error and its tolerance
          r_tol = state%space%adapt%tol_max /state%time%FinTime**0.5
          r_est =( state%T_estim(resST)/( state%time%ttime - state%time%ttime_save ))**0.5
          print*,'Not tested FWE'

       endif

       ! screen output
       if( finish ) then
          print*
          write(*,'(a40, l2, es12.4, a3, es12.4)') &
               '# hp-iso adapt ?,estim_tot, tol:', &
               state%space%adaptation, r_est, ' < ', r_tol
          write(*,'(a20,4es12.4)') 'Time local estims', state%T_estim(resA:resST)**0.5
          write(*,'(a20,4es12.4)') 'Time global estims',state%estim(resA:resST, 1)**0.5
          print*
       endif

       ! error estimate is bigger than the tolerance, we need a remeshing
       state%space%adapt%stop_adaptation = 1
       if(r_est > r_tol) then
          if( state%space%adapt%adapt_type == 'Ihp') then
             call IsotropicMetric( )
          elseif( state%space%adapt%adapt_type == 'Ahp') then
             call  AnisotInterpolEstimates( )
          else
             print*,'VERIFY edet53dei in computeAD.f90'
          endif

          state%space%adapt%stop_adaptation = 0 ! remeshing necessary
          print*
          write(*, *) ' ###  Stopping criterion was not achieved'

       endif


    elseif( state%space%estim_space == 'inter') then ! interpolation error estimates
       Lq = state%space%adapt%Lq

       r_estOLD = state%err(interL8)                      ! L^\infinity- norm
       if(Lq >=  1.) r_estOLD = state%err(interLq)        ! L^q - norm
       if(Lq <= -1.) r_estOLD = state%err(interH1)        ! H^1 - semi-norm

       call AnisotInterpolEstimates( )
       ! estimated error and its tolerance
       r_tol = state%space%adapt%tol_max

       r_est = state%err(interL8)                ! L^\infinity- norm
       if(Lq >=  1.) r_est = state%err(interLq)  ! L^q - norm
       if(Lq <= -1.) r_est = state%err(interH1)  ! H^1 - semi-norm

       ! screen output
       if( finish ) then
          print*
          write(*,'(a40, l2, es12.4, a3, es12.4)') &
               '# hp-iso adapt ?,estim_tot, tol:', &
               state%space%adaptation, r_est, ' < ', r_tol
       endif

       ! state%space%adapt%stop_adaptation is set in  AnisotInterpolEstimates( )

    elseif ( state%space%estim_space == 'RTNst' ) then
       if( state%time%disc_time == 'STDG')  then
         ! number of time nodes for integration

         ! minimally tdof +1 is absolutely necessary because of the Radau polynomial which is in P^{q+1}
         !not needed now the variable is global now
         !tQnum = state%time%max_Tdof + 1

         call ComputeRTNstEstim( state%loc_estim(1:max_eta, 1:ndim) , state%time%Qnum )

         state%estim(:, :) = state%estim(:, :) + state%loc_estim(1:max_eta,:)
         write(debug,*) 'Probably dont use RTNall since it has to be squarerooted in WriteOutputErrors'
         state%estim(RTNall, :) = sqrt( state%estim(RTNeta,:) ) + state%estim(RTNjump,:)
!         print*, 'probably problem in RTNjump - no squareroot should be used ???'

!         write(*,'(a40, l2, es12.4, a3, es12.4)') &
!            '# hp-iso adapt ?,estim_tot, tol:', &
!            state%space%adaptation, r_est, ' < ', r_tol
!         write(*,'(a20,5es12.4)') 'Time local estims', state%loc_estim(2:6,1)**0.5
!         write(*,'(a20,5es12.4)') 'Time global estims',state%estim(1:6, 1)

         write(debug,*) 'What is r_tol in ErrorEstims_time_dependent for???'
         r_tol = 1.0

       else
         stop ' Space-time RTN estimates are possible only for STDGM'
       endif

    ! Screen output
    print*,
    write(*,'(a42, a6, es12.4, a3, es12.4, a2)') &
         '## error estimates: method, estim, tol:  ', &
         state%space%estim_space, r_est, ' < ', r_tol,'?'

    elseif( state%space%estim_space == 'DUA') then
       call ComputeDualEstim( )
       print*, 'DUA'
       !write(*,'(a28,6es12.4)') 'state%estim(1:Osc,1:ndim):', state%estim(1:Osc,1)

    else
       stop ' UNKNOWN TYPE in ErrorEstims in compuateAD.f90 for TIME-DEPENDENT'
    endif



    if( r_est <= r_tol )  state%space%adapt%stop_adaptation = 1

  end subroutine ErrorEstims_time_dependent


  !> several type of error estimates
  !> TIME-INDEPENDET COMPUTATIONS
  subroutine ErrorEstims_time_independent( )
    real ::  r_tol, r_est, Lq, r_estOLD
    real :: cF, xF
    integer :: i, tQnum

    r_tol = 0.
    r_est = 0.

    !print*,'#####.. check ', state%space%adapt%stop_adaptation
    state%space%adapt%stop_adaptation = 0


    if( state%space%estim_space == 'RES' ) then ! residual error estimates

       if( state%time%disc_time == 'STDG') then
          do i=resA, resST  ! i.e., resA, resS, resT, resST
             !state%estim(i, 1) = state%estim(i, 1) + state%L_estim(i)**2 !*state%time%tau(1)
             !state%T_estim(i)  = state%T_estim(i)  + state%L_estim(i)**2 !*state%time%tau(1)

             state%estim(i, 1) = state%L_estim(i)**2 /state%time%tau(1)**2
             state%T_estim(i)  = state%L_estim(i)**2 /state%time%tau(1)**2
          enddo
       else! BDF time discretization
          call RezidErrorEstimates( .false.  , .true.)
          do i=resA, resST  ! i.e. resA, resS, resT, resST
             state%estim(i, 1) = state%L_estim(i)**2
          enddo
       endif

       ! estimated error and its tolerance
       r_est = state%estim(resS,1)**0.5
       r_tol = state%space%adapt%tol_max


    elseif(state%space%estim_space == 'RTN') then
       call ComputeApostEstim( )
       ! write(*,'(a28,6es12.4)')'state%estim(1:Osc,1:ndim):', state%estim(1:Osc,1)

    elseif( state%space%estim_space == 'DUA') then
       call ComputeDualEstim( )
       !write(*,'(a28,6es12.4)') 'state%estim(1:Osc,1:ndim):', state%estim(1:Osc,1)

    elseif( state%space%estim_space == 'RTNst' ) then
       ! number of time nodes for integration
       stop 'The RTNst estimates work only for time dependent problems!'
!       call ComputeRTNstEstim( state%estim(1:max_eta, 1:ndim) , tQnum )

    elseif( state%space%estim_space == 'HEL') then
       call ComputeHelmholtzEstim( )
       ! write(*,'(a28,6es12.4)')'state%estim(1:Osc,1:ndim):', state%estim(1:Osc,1)

    elseif( state%space%estim_space == 'HO_rec') then
       !!call ComputeHigherOrderEstims( )
       call ComputeHigherOrderEstims_OPT( )
       !call ComputeHO_LocalProblems( )

       ! estimated error and its tolerance
       r_est = state%estim(HO_estim_H1_p2,1)**0.5
       r_tol = state%space%adapt%tol_max


    elseif( state%space%estim_space == 'inter') then ! interpolation error estimates
       Lq = state%space%adapt%Lq

       r_estOLD = state%err(interL8)                      ! L^\infinity- norm
       if(Lq >=  1.) r_estOLD = state%err(interLq)        ! L^q - norm
       if(Lq <= -1.) r_estOLD = state%err(interH1)        ! H^1 - semi-norm

       ! only for settings arrays state%estim(:, HO_rec_p2_*)
       call ComputeHigherOrderEstims_OPT( )

       call AnisotInterpolEstimates( )
       ! estimated error and its tolerance
       r_tol = state%space%adapt%tol_max

       r_est = state%err(interL8)                ! L^\infinity- norm
       if(Lq >=  1.) r_est = state%err(interLq)  ! L^q - norm
       if(Lq <= -1.) r_est = state%err(interH1)  ! H^1 - semi-norm

       ! special action, need to be verified
       !if(r_est > r_estOLD .and. r_est > 1.1 * r_estOLD)  state%space%adapt%stop_adaptation = -1


    elseif( state%space%estim_space == 'pNeu') then  ! p-robust EE based on local Neumann problems

       r_estOLD = state%estim(P_tot,1)**0.5
       call ComputeLocalNeumannEstim( )

       ! estimated error and its tolerance
       r_est = state%estim(P_tot,1)**0.5
       r_tol = state%space%adapt%tol_max

       if(r_est > r_tol ) then
          if( state%space%adapt%adapt_type == 'Ihp' ) then
             cF = 0.95
             xF = 0.9
             if(r_est > cF*r_estOLD) then
                state%space%adapt%tol_min = state%space%adapt%tol_min*xF
                print*,'Security factor switch ON', state%space%adapt%tol_min
             endif

             print*,'REMOVe D'
             call ComputeHigherOrderEstims_OPT( )
             !call IsotropicMetric( )

          else
             stop 'No verified other variants in  ErrorEstims_time_independent'
          endif
       endif

    elseif( state%space%estim_space == 'DWR' ) then
       print*, 'subroutine ErrorEstims in computeAD. Still not implemented for DWR method. Using RES.'


    else
       stop ' UNKNOWN TYPE in ErrorEstims in compuateAD.f90 for TIME-INDEPENDENT'

    endif


    ! Screen output
    print*,
    write(*,'(a42, a6, es12.4, a3, es12.4, a2)') &
         '## error estimates: method, estim, tol:  ', &
         state%space%estim_space, r_est, ' < ', r_tol,'?'

    if( r_est <= r_tol )  state%space%adapt%stop_adaptation = 1
    !print*,'#####.. check2', state%space%adapt%stop_adaptation

  end subroutine ErrorEstims_time_independent




  subroutine ErrorEstimsOLD( )
    integer :: i, tQnum
    real :: t1, t2

    call cpu_time(t1)

    print*, 'Calling ErrorEstims'

    !if(state%modelName == 'scalar' .or.state%modelName == '2eqs')  then
    if(state%space%estim_space == 'RTN') then
       call ComputeApostEstim( )
       ! write(*,'(a28,6es12.4)')'state%estim(1:Osc,1:ndim):', state%estim(1:Osc,1)

    elseif( state%space%estim_space == 'DUA') then
       call ComputeDualEstim( )
       !write(*,'(a28,6es12.4)') 'state%estim(1:Osc,1:ndim):', state%estim(1:Osc,1)
!
!    elseif( state%space%estim_space == 'RTNst' ) then
!      ! number of time nodes for integration
!       tQnum = state%time%max_Tdof ! number of integration nodes
!       call ComputeRTNstEstim( state%estim(1:max_eta, 1:ndim) , tQnum )

    elseif( state%space%estim_space == 'HEL') then
       call ComputeHelmholtzEstim( )
       ! write(*,'(a28,6es12.4)')'state%estim(1:Osc,1:ndim):', state%estim(1:Osc,1)

    elseif( state%space%estim_space == 'HO_rec') then
       !!call ComputeHigherOrderEstims( )
       call ComputeHigherOrderEstims_OPT( )

    elseif( state%space%estim_space == 'inter') then
       call AnisotInterpolEstimates( )

    elseif( state%space%estim_space == 'pNeu') then
      call ComputeLocalNeumannEstim( )

    elseif( state%space%estim_space == 'DWR' ) then
         print*, 'subroutine ErrorEstims in computeAD. Still not implemented for DWR method. Using RES.'


    elseif( state%space%estim_space == 'RES' .or. state%space%estim_space == 'DWR' ) then ! .or. state%space%estim_space == 'ANI' &
       !.or. state%space%estim_space == 'Ahp' ) then

       !print*, ' else Error estimates for RES already computed in PerformOneTimeStep'
       !if(state%time%time_method /= 'I'.or. state%nlSolver%non_alg_stop /= 'aRES' .or. &
       !     (state%time%tau_choice == 'adapt' .and.  state%time%estim_time /= 'tRES') ) then
       !   if( state%time%disc_time /= 'STDG') then
       !      call RezidErrorEstimates( .false. , .true.)
       !   else
       ! print*,'computeAD - errorEstims calling RezidErrorEst (false)!! NOT IMPLEMENTED'
       ! for  state%time%disc_time == 'STDG'already computed
       if( state%time%disc_time /= 'STDG')  &
            call RezidErrorEstimates( .false.  , .true.)
       !   endif
       !endif

       if(state%time_dependent) then

          if( state%time%disc_time /= 'STDG') then
             do i=resA, resST  ! i.e., resA, resS, resT, resST
                state%estim(i, 1) = state%estim(i, 1) + state%L_estim(i)**2 *state%time%tau(1)
                state%T_estim(i)  = state%T_estim(i)  + state%L_estim(i)**2 *state%time%tau(1)
             enddo
          else
             do i=resA, resST  ! i.e., resA, resS, resT, resST
                state%estim(i, 1) = state%estim(i, 1) + state%L_estim(i)**2 !*state%time%tau(1)
                state%T_estim(i)  = state%T_estim(i)  + state%L_estim(i)**2 !*state%time%tau(1)
             enddo
          endif
       else
          if( state%time%disc_time /= 'STDG') then
             do i=resA, resST  ! i.e. resA, resS, resT, resST
                state%estim(i, 1) = state%L_estim(i)**2
             enddo
          else
             do i=resA, resST  ! i.e., resA, resS, resT, resST
                !state%estim(i, 1) = state%estim(i, 1) + state%L_estim(i)**2 !*state%time%tau(1)
                !state%T_estim(i)  = state%T_estim(i)  + state%L_estim(i)**2 !*state%time%tau(1)

                state%estim(i, 1) = state%L_estim(i)**2 /state%time%tau(1)**2
                state%T_estim(i)  = state%L_estim(i)**2 /state%time%tau(1)**2
             enddo
          endif
       endif

       !open(99, file='estims', status='UNKNOWN', position = 'append')
       !write(99,'(a6,i5,30es12.4)') 'ETAs:', state%time%iter, &
       !     state%time%tau(1), state%L_estim(1:4), state%estim(1:4, 1)**0.5
       !close(99)
    else
       stop ' UNKNOWN TYPE in ErrorEstims in compuateAD.f90'

    endif

    call cpu_time(t2)
    state%CPU_estim2 = state%CPU_estim2 + t2 - t1


  end subroutine ErrorEstimsOLD

  !> Perform some time steps on the given grid
  subroutine SolveProblemAD(convfile)
    use matrix_oper_int
    character(len=*), intent(in) :: convfile
    class(element), pointer :: elem
    integer :: iter, istep, i, j, k
    real :: t_sta, time_prepare, time_solve, time_estim, diff, diff1
    logical :: refuseStep, finish, loc_implicitly!!, CrankNicolson

    refuseStep = .false.
    state%linSolver%residuum = 0.

    !!CrankNicolson = state%time%cn

    call state%nlSolver%InitNLSolverSettings( state%time%time_method )

    print*
    write(*,'(a31, 2i2,a2, a12, a8, a8, es12.4)')' # SolveProblemAD starts (Tdeg=', &
         state%time%deg_actual, state%time%deg,') ', state%linSolver%name,state%time%disc_time, &
         ', tau = ',state%time%tau(1)
    print*


    call cpu_time(t_sta) ! this variable is overwritten when solution is written on output

    state%T_estim(:) = 0.
    state%errSTloc(:) = 0.
    state%time_AD_start = state%time%ttime

    do iter = 1,state%time%maxiter
       !call cpu_time(ts)
       state%time%iter_loc = iter

       !print*,'PrepareNewTimeStepAD(iter)'
       call PrepareNewTimeStepAD(iter)

       !! switch to Crank/Nicolson
       !if(state%time%deg_actual == 1 .and. state%time%deg > 1 .and. state%time_dependent) then
       !   print*,'switch to Crank-Nicolson'
       !   state%time%cn = .true.
       !endif

       state%linSolver%iter = 0

       !! cycles repeated if the time step is refused
       do istep = 1, 5 ! 10
          state%time%ctime = state%time%ttime + state%time%tau(1)


          if (state%SP) then  ! saddle point
             !print*,' perform one time step', state%time%tau(:)
             stop 'call PerformOneTimeStepSP( time_prepare, time_solve, time_estim )'

          else ! not
             ! for BDF, not for STDGM

             if( state%time%disc_time /= 'STDG') then
                !print*,' perform one time step', state%time%tau(:)

                !write(*,*)
                !write(*,*) "---------------------------------------------"  ! wet steam test
                !write(*,*) 'w', grid%elem(1)%w(0,1:ndim)
                !write(*,*) "SolveProblemAD calls PerformOneTimeStep; time iter", iter        ! wet steam test
                !write(*,*) "---------------------------------------------"  ! wet steam test
                !write(*,*)
                call PerformOneTimeStep( time_prepare, time_solve, time_estim )

             else  ! STDGM
                !print*,' perform one STDGM step', state%time%tau(:)
                call PerformOneSTDGMstep( time_prepare, time_solve, time_estim )
                !print*, 'End of PerformOneSTDGMstep in SolveproblemAD', 'state%space%adapt%adapt_method', state%space%adapt%adapt_method
             endif
             !print*,'DESWD',grid%elem(1)%w(0,:)
          endif

          !if(state%time%tau_choice /= 'fixed') then
          if(state%time%tau_choice == 'adapt') then
             state%time%tau_new = state%time%tau(1)
             !print*,' adaptive choice of the time step'
             call ProposeNewTimeStep(refuseStep )  ! for fixed time only error estim

             !  if(.not. refuseStep) goto 90
          endif

          if(.not. refuseStep) goto 90

       enddo

90     continue

       !print*
       !print*
       !print*,'Tests of the projection ONLY in computeAD.f90'
       !call SetElementsIC()
       !print*
       !print*
       !print*

       !loc_implicitly = state%nlSolver%implicitly
       !state%nlSolver%implicitly = .false.
       !call ComputeSTDGM_Terms( )
       !state%nlSolver%implicitly = loc_implicitly

       ! write(*,'(a4,20es12.4)') '***',0., grid%elem(121)%rhsST(1, :, 1)
       ! write(*,'(a4,20es12.4)') '***',0., grid%elem(121)%rhsST(1, :, 2)
       ! write(*,'(a4,20es12.4)') '***',0., grid%elem(121)%rhsST(1, :, 3)
       ! write(*,'(a4,20es12.4)') '***',0., grid%elem(121)%wST(:, :, :)
       ! print*,'------------------', state%time%ctime

       !print*,' call JumpsEvaluation( ) '
       call JumpsEvaluation( )

       !print* ,' one time step was perfomed, we go its evaluation'
       call PassToNewTimeStep( )

       !!state%time%cn = CrankNicolson

       !loc_implicitly = state%nlSolver%implicitly
       !state%nlSolver%implicitly = .false.
       !state%time%ttime = state%time%ttime - state%time%tau(1)
       !call ComputeSTDGM_Terms( )
       !state%time%ttime = state%time%ttime + state%time%tau(1)
       !state%nlSolver%implicitly = loc_implicitly

       ! write(*,'(a4,20es12.4)') '*A*',0., grid%elem(121)%rhsST(1, :, 1)
       ! write(*,'(a4,20es12.4)') '*A*',0., grid%elem(121)%rhsST(1, :, 2)
       ! write(*,'(a4,20es12.4)') '*A*',0., grid%elem(121)%rhsST(1, :, 3)
       ! write(*,'(a4,20es12.4)') '***',0., grid%elem(121)%wST(:, :, :)
       ! print*,'------------------', state%time%ctime

       state%time%ctime = state%time%ttime  ! ctime was changed, we keep the actual value for Estims
       !  print*, 'calling ErrorEstims in SolveProblemAD after passToNewTimeStep'
       !write(*,'(a10,3es12.4)') 'Errors:',state%err(L2), state%err(H1), state%err(H1_discrete)

       call WriteOutputScreen (iter, time_prepare, time_solve, time_estim,  finish)

       if(state%time_dependent) call ErrorEstims( finish  )

       call WriteOutputFiles (convfile, iter, t_sta, time_prepare, time_solve, time_estim )

70     continue

       !!print*,'Supress output here - ehdtedeiu'
       !!call WriteProgressOutput( 'S', .false. )
       !!state%isol = state%isol + 1

       if(finish ) goto 100

    enddo  ! go to next time step

100 continue
    state%time_AD_end = state%time%ttime

    if(.not. state%time_dependent) call ErrorEstims( finish )

    !print*, 'Before AnisotInterpolEstimates'
    !
    ! computation of the metric and iterpol error estims
    !if(state%space%adapt%max_adapt_level > 0 .and. &
    !     (state%space%adapt%adapt_type == 'Ahp' .or. state%space%adapt%adapt_type == 'ANI' ) ) &
    !     call AnisotInterpolEstimates( )

    ! output the final estimates of the error
    call WriteOutputError( )

    open (11, file=convfile, status='OLD', position = 'append')
    write(11, '(x)' )
    close(11)

    if( state%MGsolver ) call DeInitMG( )


    !open(91, file='../estimST', status='unknown', position='append')
    !write(91,*) state%time%tau(1), state%space%h, &
    !     state%L_estim(resS), state%G_estim(resT), 1./state%time%tau(1), &
    !     state%L_estim(resS)/ state%time%tau(1)**0.5, &
    !     state%L_estim(resT)/ state%time%tau(1)**0.5
    !close(91)

    !print*,'# Iteration process SolveProblemAD finished  '

    if( (state%modelName == 'scalar' .or. state%modelName == '2eqs') .and. state%type_IC .eq. 5) then
       !print*,'Output for LEVEL SET METHOD in file "data.vofx"'
       call WriteOutput_LevelSet()
    endif

    !print*,'_____________________________________________________________________________'

    !print*,'##### END of SolveProblemAD'

  end subroutine SolveProblemAD

   !FR not done
    !> Perform some time steps on the given grid
  subroutine SolveDualProblem(convfile)
    use matrix_oper_int
    character(len=*), intent(in) :: convfile
    class(element), pointer :: elem
    integer :: iter, istep, i, j, k
    real :: t_sta, time_prepare, time_solve, time_estim, diff, diff1
    logical :: refuseStep, finish, loc_implicitly!!, CrankNicolson

    refuseStep = .false.
    state%linSolver%residuum = 0.

    !!CrankNicolson = state%time%cn

    !??? when to use this and how - DWR%newton vs %newton
!    if(size(Newton%b(:) ) /= state%nsize  )then
!       deallocate( Newton%b, Newton%b1, Newton%x, Newton%rr )
!
!       allocate(Newton%b(1:state%nsize), Newton%b1(1:state%nsize), &
!            Newton%x(1:state%nsize), Newton%rr(1:state%nsize) )
!    endif


!    call InitNewtonSettingsAD( state%DWR%Newton )

    print*
    write(*,'(a50, 2i3,a1, 2a12, a8, es12.4)')' # Iteration process SolveProblemAD starts (Tdeg=', &
         state%time%deg_actual, state%time%deg,')', state%linSolver%name,state%time%disc_time, &
         ', tau = ',state%time%tau(1)
    print*


    call cpu_time(t_sta) ! this variable is overwritten when solution is written on output

    state%T_estim(:) = 0.
    state%errSTloc(:) = 0.
    state%time_AD_start = state%time%ttime

    do iter = 1,state%time%maxiter
       !call cpu_time(ts)
       state%time%iter_loc = iter

       !print*,'PrepareNewTimeStepAD(iter)'
       call PrepareNewTimeStepAD(iter)

       !! switch to Crank/Nicolson
       !if(state%time%deg_actual == 1 .and. state%time%deg > 1 .and. state%time_dependent) then
       !   print*,'switch to Crank-Nicolson'
       !   state%time%cn = .true.
       !endif

       state%linSolver%iter = 0

       !! cycles repeated if the time step is refused
       do istep = 1, 5 ! 10
          state%time%ctime = state%time%ttime + state%time%tau(1)

          ! for BDF, not for STDGM
          if( state%time%disc_time /= 'STDG') then
             !print*,' perform one time step', state%time%tau(:)
             call PerformOneTimeStep( time_prepare, time_solve, time_estim )

          else  ! STDGM
             !print*,' perform one STDGM step', state%time%tau(:)
             call PerformOneSTDGMstep( time_prepare, time_solve, time_estim )
             !print*, 'End of PerformOneSTDGMstep in SolveproblemAD', 'state%space%adapt%adapt_method', state%space%adapt%adapt_method
          endif

          !if(state%time%tau_choice /= 'fixed') then
          if(state%time%tau_choice == 'adapt') then
             state%time%tau_new = state%time%tau(1)
             !print*,' adaptive choice of the time step'
             call ProposeNewTimeStep(refuseStep )  ! for fixed time only error estim

             !  if(.not. refuseStep) goto 90
          endif

          if(.not. refuseStep) goto 90

       enddo

90     continue

       !loc_implicitly = state%nlSolver%implicitly
       !state%nlSolver%implicitly = .false.
       !call ComputeSTDGM_Terms( )
       !state%nlSolver%implicitly = loc_implicitly

       ! write(*,'(a4,20es12.4)') '***',0., grid%elem(121)%rhsST(1, :, 1)
       ! write(*,'(a4,20es12.4)') '***',0., grid%elem(121)%rhsST(1, :, 2)
       ! write(*,'(a4,20es12.4)') '***',0., grid%elem(121)%rhsST(1, :, 3)
       ! write(*,'(a4,20es12.4)') '***',0., grid%elem(121)%wST(:, :, :)
       ! print*,'------------------', state%time%ctime

       ! one time step was perfomed, we go its evaluation
       call PassToNewTimeStep( )

       !!state%time%cn = CrankNicolson

       !loc_implicitly = state%nlSolver%implicitly
       !state%nlSolver%implicitly = .false.
       !state%time%ttime = state%time%ttime - state%time%tau(1)
       !call ComputeSTDGM_Terms( )
       !state%time%ttime = state%time%ttime + state%time%tau(1)
       !state%nlSolver%implicitly = loc_implicitly

       ! write(*,'(a4,20es12.4)') '*A*',0., grid%elem(121)%rhsST(1, :, 1)
       ! write(*,'(a4,20es12.4)') '*A*',0., grid%elem(121)%rhsST(1, :, 2)
       ! write(*,'(a4,20es12.4)') '*A*',0., grid%elem(121)%rhsST(1, :, 3)
       ! write(*,'(a4,20es12.4)') '***',0., grid%elem(121)%wST(:, :, :)
       ! print*,'------------------', state%time%ctime

       state%time%ctime = state%time%ttime  ! ctime was changed, we keep the actual value for Estims
       !  print*, 'calling ErrorEstims in SolveProblemAD after passToNewTimeStep'

       call WriteOutputScreen (iter, time_prepare, time_solve, time_estim, finish)

       call ErrorEstims( finish )

       call WriteOutputFiles (convfile, iter, t_sta, time_prepare, time_solve, time_estim )

70     continue

       if(finish ) goto 100

    enddo  ! go to next time step

100 continue
    state%time_AD_end = state%time%ttime

    call JumpsEvaluation( )

    ! computation of the metric and iterpol error estims
    if(state%space%adapt%max_adapt_level > 0 .and. &
         (state%space%adapt%adapt_type == 'Ahp' .or. state%space%adapt%adapt_type == 'ANI' ) ) &
         call AnisotInterpolEstimates( )

    ! output the final estimates of the error
    call WriteOutputError( )

    open (11, file=convfile, status='OLD', position = 'append')
    write(11, '(x)' )
    close(11)

    if( state%MGsolver ) call DeInitMG( )


    !open(91, file='../estimST', status='unknown', position='append')
    !write(91,*) state%time%tau(1), state%space%h, &
    !     state%L_estim(resS), state%G_estim(resT), 1./state%time%tau(1), &
    !     state%L_estim(resS)/ state%time%tau(1)**0.5, &
    !     state%L_estim(resT)/ state%time%tau(1)**0.5
    !close(91)

    !print*,'# Iteration process SolveProblemAD finished  '

    if( (state%modelName == 'scalar' .or. state%modelName == '2eqs') .and. state%type_IC .eq. 5) then
       print*,'Output for LEVEL SET METHOD in file "data.vofx"'
       call WriteOutput_LevelSet()
    endif

    print*,'_____________________________________________________________________________'

    !print*,'##### END of SolveDualProblem'

  end subroutine SolveDualProblem

  !> call the approxiate mesh adaptation technique, create new mesh: grid
  subroutine AdaptationAD(  )
    logical :: identical_AMA_grids
    logical :: metric
    integer :: indexi
    real ::  r_tol, r_est

    print*,'_____________________________________________________________________________'

    !call WriteResults('Gsol.bak')
    if(state%space%adapt%adapt_type == 'Ahp' .or. state%space%adapt%adapt_type == 'Ihp') then ! ANGENER based on HO interpol


       ! if(state%space%estim_space == 'pNeu') then
       !    indexi = P_tot
       ! elseif(state%space%estim_space == 'HO_rec') then
       !    indexi = HO_estim_H1_p2

       ! elseif(state%space%estim_space == 'RES') then
       !    if(.not. state%time_dependent) indexi = resS
       !    if(      state%time_dependent) indexi = resST

       ! else
       !    print*,'undefined value in  AdaptationAD(  ), stopping!'
       !    stop
       ! endif


       metric = .false.   !metric = .true. ==> metric computed by ANGENER F77
       identical_AMA_grids = .false.

       ! already  done in compute1
       !!if(.not. metric) call AnisotInterpolEstimates( )


       ! if(state%space%adapt%adapt_type == 'Ihp' .or. state%space%adapt%adapt_method == 'ANI' ) then
       !    state%space%adapt%stop_adaptation = 0

       !    write(*,*)
       !    if(state%time_dependent) then
       !       if( state%time%disc_time /= 'STDG') then
       !          r_tol = state%space%adapt%tol_max /state%time%FinTime**0.5
       !          r_est =( state%T_estim(indexi)/( state%time%ttime - state%time%ttime_save ))**0.5
       !          print*,'Not tested FWE'
       !       else
       !          r_tol = state%space%adapt%tol_max * (( state%time%ttime - state%time%ttime_save ) / state%time%FinTime)**0.5
       !          r_est = state%T_estim(indexi)**0.5
       !       endif

       !       write(*,'(a40, l2, es12.4, a3, es12.4)') &
       !            '# hp-iso adapt ?,estim_tot, tol:', &
       !            state%space%adaptation, r_est, ' < ', r_tol
       !       write(*,'(a20,4es12.4)') 'Time local estims', state%T_estim(resA:resST)**0.5
       !       write(*,'(a20,4es12.4)') 'Time global estims',state%estim(resA:resST, 1)**0.5
       !    else
       !       r_est = state%estim(indexi,1)**0.5
       !       r_tol = state%space%adapt%tol_max
       !       write(*,'(a40, l2, es12.4, a3, es12.4)') &
       !            '# hp-iso adapt ?,estim_tot, tol:', &
       !            state%space%adaptation, r_est, ' < ', r_tol
       !    endif

       !    write(*,*)


       !    if( r_est <= r_tol )  state%space%adapt%stop_adaptation = 1

       !    if( state%space%adapt%adapt_type == 'Ihp' .and. state%space%adapt%stop_adaptation == 0 )  then
       !       if( state%space%estim_space == 'HO_rec') then
       !          ! already done, SHOULD be rewritten


       !          ! OLD VARIANT using elem%reg, elem%regT2
       !          !call SetHigherOrder_HP_adapt( )

       !          ! NEW variant - not working, direct setting in ComputeHigherOrderElem_HP
       !          ! NEWEST variant - direct setting in ComputeHigherOrderElem_OPT

       !       else
       !          !call IsotropicMetricOLD( )
       !          call IsotropicMetric( )
       !       endif

       !    endif


       ! endif

       !if(state%space%adapt%stop_adaptation > 0 ) then
       !   print*
       !   print*,' # Sucesfull multi-time step'
       !   print*
       !endif

       !print*,'SW3x'
       !pause

       !if( state%space%adapt%stop_adaptation <= 0 ) then  !!!.and. state%time%ttime < state%time%FinTime) &
       !print*
       !   write(*, *) ' ###  Stopping criterion was not achieved'
       call AdaptMesh_Angener(metric, identical_AMA_grids)
       !endif


       !print*,'SW3y'
       !pause

       if(identical_AMA_grids) state%space%adapt%stop_adaptation = 11

    elseif ( state%space%adapt%adapt_type == 'HG' .or. state%space%adapt%adapt_type == 'RG') then

       ! !!call PlotMesh(grid, 'mesh')

       call MarkElements( )

       ! if(.not. state%time_dependent) then

       !    if(state%space%estim_space == 'pNeu' ) then
       !       write(*,*)
       !       write(*,'(a40, l2, es12.4, a3, es12.4)') &
       !            '# hp-iso adapt ?,estim_tot, tol:', &
       !            state%space%adaptation, state%estim(P_tot, 1)**0.5, &
       !            ' < ', state%space%adapt%tol_max
       !       write(*,*)
       !       if( state%estim(P_tot, 1) <= state%space%adapt%tol_max**2 ) &
       !             state%space%adapt%stop_adaptation = 1

       !    else

       !       !if(state%space%adapt%adapt_method == 'AMA') then
       !       !!write(*,'(a25, l2, 3es12.4, a4, es12.4)') '# hp-isotropic adapt  ?', &

       !       write(*,'(a40, l2, 3es12.4, a3, es12.4)') &
       !            '# hp-iso adapt ?,reS, resJ, resSJ, tol:', &
       !            state%space%adaptation,  state%L_estim(resS),  state%err(EjumpsG), &
       !            (state%L_estim(resS)**2 + state%err(EjumpsG)**2)**0.5, ' <?', state%space%adapt%tol_max

       !       if(.not. state%space%adaptation) state%space%adapt%stop_adaptation = 2

       !       if(state%L_estim(resS)**2 + state%err(EjumpsG)**2 < state%space%adapt%tol_max**2 ) &
       !            state%space%adapt%stop_adaptation = 1

       !    endif
       ! endif

       !if(state%space%adapt%stop_adaptation <= 0) then

          ! hp variant
          if(state%space%adapt%adapt_type == 'HG') then
             call AdaptMesh_hp_HG( )

          elseif(state%space%adapt%adapt_type == 'RG') then
             call AdaptMesh_hp_RG( )

          endif

          !! plotting of recomputed solution
          !state%space%adapt%adapt_level = state%space%adapt%adapt_level + 1
          !call WriteProgressOutput( 'S' )
          !state%space%adapt%adapt_level = state%space%adapt%adapt_level - 1

          ! if(state%space%adapt%adapt_level > 0) then
          !    is = int(log(1.*state%space%adapt%adapt_level)/log(10.))
          ! else
          !    is = 0
          ! endif

          ! num_size = 1
          ! text_size = 9
          ! write( ch1, '(i1)' ) state%space%adapt%adapt_level

          ! newgridfile = 'new.grid-0'
          ! newgridfile(num_size+text_size-is:num_size+text_size) = ch1(num_size-is:num_size)
          ! call WriteMesh(newgridfile, grid)

          if(state%time%deg > 2) print*,'### !! change state%time%deg_actual = min(2, Tdeg)'
          !state%time%deg_actual = 1
          state%time%deg_actual = min(2, state%time%deg)

       !end if

    else
       print*,'UNKNOWN state%space%adapt%adapt_method, stopped in AdaptationAD'
       stop

    end if
    print*,'_____________________________________________________________________________'

  end subroutine AdaptationAD


  !> used for the adaptation of time dependent problems
  !> recomputation back in time, interpolation results from gridS on new grid
  subroutine RecomputeBackMesh( ) !!lev)
    !integer, intent(in) :: lev
    integer :: i !, lev

    !lev = state%time%recompute_back
    !print*
    write(*,'(a55,i2)') &
         ' Unsuccessfull multi-time step, BACK recomputation No.',state%time%recompute_back
    !print*

    state%time%ttime = state%time%ttime_save
    state%time%ctime = state%time%ttime_save
    state%timeprn = state%timeprn_save
    state%isol = state%isol_save

    if(state%time%ttime > 0.) then
       call AdvancedInterpolDGsolution(grid, gridS )
       !call SimpleInterpolDGsolution(grid, gridS )
       !print*,'##### AdvancedInterpolDGsolution(grid, gridS )'
    else
       call SetElementsIC()

       if(state%time%disc_time == 'STDG') then  ! ST DGM
          do i=1,grid%nelem    ! setting of elem%w(:,:)
             call InitElementW(grid%elem(i) )
             !call InitElementEigenvals(grid%elem(i) )
          enddo
       end if

       !print*,'#####   call SetElementsIC()'
    endif

    !recomputation back of ST errors and estims: resA, resS, resT, resST
    state%estim(resA:resST, 1) = state%estim(resA:resST, 1) &
         - state%T_estim(resA:resST)

    state%errSTnorm(L2L2:)  = state%errSTnorm(L2L2:) - state%errSTloc(L2L2:)
    !print*,' RecomputeBackMesh DONE'

  end subroutine RecomputeBackMesh



  !> main procedure, compute the solution including adaptation - Moved to o_scompute
  !>
  !> three types of computation
  !> 1) NO ADAPTATION, steady as well as unsteady
  !> 2)    ADAPTATION, steady
  !> 3)    ADAPTATION, unsteady, reinterpolation back is necessary
!  subroutine ComputeAD(convfile, solfile)
!    character(len=*), intent(in) :: convfile, solfile
!    !character(len=50) :: command_name, tri_name, sol_name
!    character(len=50) :: newgridfile
!    integer :: i, i0, j, max_recomput_back
!    character(len=1) :: ch1
!    integer :: text_size, is, num_size
!    logical :: Aname
!
!    ! COMPUTATION   1) NO ADAPTATION, steady as well as unsteady
!    if(state%space%adapt%max_adapt_level == 0) then
!       Aname = .false.
!
!       !do i=0, state%space%adapt%max_adapt_level
!       !state%space%adapt%adapt_level = i
!
!
!       if(state%time%OutTime > 0.) then
!          !if (state%time%disc_time .ne. 'STDG') then
!
!          call WriteProgressOutput( 'TS', Aname )
!          !elseif (state%isol == 0) then
!          !   call WriteProgressOutput( 'TS' )
!          !endif
!       else
!          call WriteProgressOutput( 'T', Aname )
!       endif
!
!       !call WriteProgressOutput( 'E' , Aname)
!
!       call SolveProblemAD(convfile)
!
!       call WriteResults(solfile)
!
!       if(state%space%estim_space == 'pNeu' .or. state%space%estim_space == 'HO_rec') &
!            call WriteProgressOutput('A', Aname)
!
!       if(ndim <= 2)call WriteProgressOutput( 'SE', Aname )
!       if(ndim >  2)call WriteProgressOutput( 'S' , Aname)
!       !if(state%type_IC == 7) call WriteProgressOutput( 'E' )
!
!       ! end of COMPUTATION   1)
!    else
!
!       if(.not. state%time_dependent) then
!          ! COMPUTATION   2) ADAPTATION, steady
!          Aname = .true.
!
!          ! only for initialization
!          if(state%space%estim_space == 'RES' .or. state%space%estim_space == 'Ahp') then
!             state%time%deg_actual = 1
!             call RezidErrorEstimates( .true., .false. )
!             call JumpsEvaluation( )
!          endif
!
!          do i= 0,  state%space%adapt%max_adapt_level
!             state%space%adapt%adapt_level = i
!!!!call WriteProgressOutput( 'TS' )
!
!             call SolveProblemAD(convfile)
!
!             if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
!                call WriteProgressOutput( 'STEA' , Aname)
!             else
!                call WriteProgressOutput( 'STA', Aname )
!             endif
!!!!call WriteProgressOutput('A')
!            if ( state%space%estim_space == 'DWR' ) then
!               call PrepareDualProblem( )
!               print*, 'Solve dualProblem not implemented. Stop.'
!               stop
!
!               call SolveDualProblem( convfile )
!            endif
!
!
!             if(i < state%space%adapt%max_adapt_level .and.  state%space%adapt%stop_adaptation==0 ) &
!                  call AdaptationAD( )
!
!             !! SMAZ
!
!             !!state%space%adapt%adapt_level = state%space%adapt%adapt_level + 1
!
!             !!call WriteProgressOutput( 'ST', Aname )
!             !!state%space%adapt%adapt_level = state%space%adapt%adapt_level - 1
!             !! SMAZ END
!
!             if(state%space%adapt%stop_adaptation > 0) goto 10
!
!          enddo! state%space%adapt%max_adapt_level
!
!10        continue
!
!          if(state%space%adapt%stop_adaptation == 1) then
!             write(*,*) ' # Error criterion achieved, no further adaptation'
!          elseif(state%space%adapt%stop_adaptation == 2) then
!             write(*,*) ' # No element marked for adaptation'
!          elseif(state%space%adapt%stop_adaptation == 11) then
!             ! written in anisot.f90
!             !write(*,*) 'The new hp-grid is (almost) identical with the previous ones'
!          endif
!
!       else
!          ! COMPUTATION   3) ADAPTATION, unsteady
!          Aname = .true.
!
!          ! only for initialization
!          if(state%space%estim_space == 'RES' .or. state%space%estim_space == 'Ahp') then
!             state%time%deg_actual = 1
!             if( state%time%disc_time /= 'STDG') then
!                call RezidErrorEstimates( .true., .false. )
!             endif
!             call JumpsEvaluation( )
!          endif
!
!          state%space%adapt%adapt_level = 0
!          if(state%tri_solA)      call WriteProgressOutput( 'ST', Aname)
!          if(state%time%OutTime > 0.)  call WriteProgressOutput( 'ST', .false. )
!          !print*,'SW1'
!          !pause
!
!          max_recomput_back =  8  !!5
!          do i= 1,  state%space%adapt%max_adapt_level + 1
!             !state%space%adapt%adapt_level = i
!             !call WriteProgressOutput( 'TS' )
!
!             call SaveOrigMesh( )
!
!             do j=1, max_recomput_back  ! inner loop within one time interval
!                state%time%recompute_back = j
!                ! REM TEST
!                !state%isol = state%isol + 1
!                !call WriteProgressOutput( 'STE' )
!
!!!!!if(state%isol == 6) stop
!
!                !print*,'SW2'
!                !pause
!                call SolveProblemAD(convfile)
!
!                ! REM TEST
!                !state%isol = state%isol + 1
!                !call WriteProgressOutput( 'STE' )
!
!                !if(state%modelName == 'scalar' .or.state%modelName == '2eqs') call WriteProgressOutput( 'STE' )
!                !if(ndim >  1) call WriteProgressOutput( 'ST' )
!                ! !!call WriteProgressOutput('A')
!
!                !print*,'before adaptation '
!                !pause
!
!!!! mesh from the initial grid is always recomputed
!                !!if(i <= 2 .and. j == 1) state%space%adapt%stop_adaptation = 0
!
!
!                call AdaptationAD( )
!
!                !print*,'after adaptation '
!                !pause
!
!                !!stop
!
!                !state%isol = state%isol + 1
!                !call WriteProgressOutput( 'ST', Aname )
!
!                if(state%space%adapt%stop_adaptation <= 0 .and. j < max_recomput_back) then
!                   call RecomputeBackMesh(  )
!
!                   ! saving of the initial condition on the new mesh
!                   if(state%space%adapt%adapt_level == 0 .and. state%tri_solA) &
!                        call WriteProgressOutput( 'ST', Aname )
!                   if(state%space%adapt%adapt_level == 0 .and. state%time%OutTime > 0.)  &
!                        call WriteProgressOutput( 'ST', .false. )
!                   state%time%keep_tau = .false.
!                else
!                   state%errSTnorm(L8L2) = max(state%errSTnorm(L8L2),state%errSTloc(L8L2))
!                   print*
!                   print*,'Sucesfull multi-time step'
!                   print*
!                   state%err(interLq) = 1E+30
!                   state%err(interL8) = 1E+30
!
!                   state%space%adapt%adapt_level = state%space%adapt%adapt_level + 1
!
!                   if(state%space%gridS_allocated ) call RemoveOrigMesh( )
!                   state%time%keep_tau = .true.
!
!                   state%space%adapt%stop_adaptation = 1  ! only for a possible exit
!
!                   if(state%tri_solA) call WriteProgressOutput( 'ST', Aname )
!
!                   goto 20
!                endif
!
!                !print*,'SW4'
!                !pause
!
!             enddo
!
!20           continue
!
!             !state%space%adapt%adapt_level = state%space%adapt%adapt_level + 1
!             !if(state%time%OutTime == 0.) then
!             !   state%isol = state%isol + 1
!             !   call WriteProgressOutput( 'ST', Aname )
!             !endif
!
!             ! Final Time achieved
!             if(state%time%ttime >= state%time%FinTime .and. state%space%adapt%stop_adaptation > 0) goto 30
!
!          enddo! state%space%adapt%max_adapt_level
!30        continue
!
!       endif
!
!
!    endif
!
!  end subroutine ComputeAD



  subroutine CleanMemory ( )
    integer :: i,j
    !type(Newton_type), pointer:: Newton

    !Newton => state%nlSolver

    do i=1,state%numBC
       deallocate(state%BC(i)%w)
       deallocate(state%BC(i)%ww)
    enddo

    deallocate (state%BC )

    do i = 1, maxGrule
       deallocate (state%space%G_rule(i)%weights, state%space%G_rule(i)%lambda)
       deallocate (state%space%G_rule(i)%phi, state%space%G_rule(i)%Dphi)
       !j = i + maxVrule
       !deallocate (state%space%V_rule(j)%weights, state%space%V_rule(j)%lambda)
       !deallocate (state%space%V_rule(j)%phi, state%space%V_rule(j)%Dphi)
    enddo
    deallocate (state%space%G_rule)


    do i = 1, maxVrule + maxGrule
       if( state%space%V_rule(i)%def) then
          deallocate (state%space%V_rule(i)%weights, state%space%V_rule(i)%lambda)
          deallocate (state%space%V_rule(i)%phi, state%space%V_rule(i)%Dphi)
       endif
    enddo
    deallocate (state%space%V_rule)

    do i = 1, maxTrule
       deallocate (state%time%T_rule(i)%weights, state%time%T_rule(i)%lambda)
       deallocate (state%time%T_rule(i)%phi, state%time%T_rule(i)%Dphi)
    enddo
    deallocate (state%time%T_rule)

!    if (state%time%disc_time == 'STDG') then
!       do i = 1, maxTrule
!         ! deallocate (state%time%T_rule(i)%weights, state%time%T_rule(i)%lambda)
!          !!deallocate (state%time%T_rule(i)%phi, state%time%T_rule(i)%Dphi)
!       enddo
!       deallocate (state%time%T_rule)
!    endif

    do i = 0, maxLrule
       deallocate ( state%space%L_rule(i)%lambda)
       deallocate (state%space%L_rule(i)%phi, state%space%L_rule(i)%psi, state%space%L_rule(i)%psi_pow)

       j = i+maxVrule
       deallocate ( state%space%L_rule(j)%lambda)
       deallocate (state%space%L_rule(j)%phi)

    enddo
    deallocate (state%space%L_rule)

    deallocate(state%space%Qdeg, state%space%ldeg)

    deallocate(state%err)
    if (state%time%disc_time == 'STDG') deallocate( state%errSTnorm, state%errSTloc, state% errSnorm)

    deallocate( state%nlSolver%x )
    deallocate( state%nlSolver%b )
    deallocate( state%nlSolver%b1 )
    deallocate( state%nlSolver%rr )

    deallocate(state%space%GScoeff)
    if ( allocated( state%space%adapt%RGred )) &
      deallocate ( state%space%adapt%RGred )
    if ( allocated( state%space%adapt%RGgreen )) &
      deallocate ( state%space%adapt%RGgreen )

    deallocate(state%estim, state%T_estim, state%L_estim, state%eta, state%time%tau, state%time%rat_tau)

    deallocate(state%cDLM)

    if(allocated( state%space%adapt%AMAhistory) ) deallocate( state%space%adapt%AMAhistory)



    call DeallocateGrid(grid)

    do i = 0, MaxRTNImplemented
       deallocate(state%RTN(i)%Vnode, state%RTN(i)%phi, state%RTN(i)%Dphi, state%RTN(i)%MM, &
            state%RTN(i)%MMinv)
       ! MISSING
    enddo
    deallocate( state%loc_RTN)


  end subroutine CleanMemory

end module computeAD_oper
