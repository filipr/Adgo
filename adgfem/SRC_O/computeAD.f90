!> new (more compact)
!> main iterative loop, solution of the unstationary compressible
!> Navier-Stokes equations

module computeAD_oper
  use main_data
  use data_mod
  use dual_problem_mod
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
  use st_interpol
  use anisotropic
  use ama_hp_interpol
  use ama_hp_resid
  use error_subs
  use errp_estimates_mod
  use mesh_oper
  use loc_rav_tho_ned
  use compute_oper
  use higher_order_estims
  use rtn_st_mod
  use higher_order_local
  use regularity
  use solution_mod
  use higher_order_recovery

!  use mesh_adapt95
!  use hp_ama_main

  implicit none

  public:: ErrorEstims
  public:: ErrorEstims_time_dependent
  public:: ErrorEstims_time_independent
  public:: SolveProblemAD
  public:: RecomputeBackMesh
  public:: AdaptationAD

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
    real ::  r_tol, r_est, r_estOLD, Lq, val
    real :: fac_ref_min,  fac_ref_max
    integer :: i, tQnum
    integer :: spaceLevels, timeLevels

    call state%cpuTime%startEstimTime()

    r_tol = 0.
    r_est = 0.

    fac_ref_min = 0.1
    fac_ref_max = 1.1

    if( state%space%estim_space == 'RES' ) then
       ! STDGM
       if( state%time%disc_time == 'STDG')  then
          do i=resA, resST  ! i.e., resA, resS, resT, resST
             state%estim(i, 1) = state%estim(i, 1) + state%L_estim(i)**2 !*state%time%tau(1)
             state%T_estim(i)  = state%T_estim(i)  + state%L_estim(i)**2 !*state%time%tau(1)
          enddo

          ! estimated error and its tolerance
          r_tol = state%space%adapt%tol_max  &
               * sqrt(( state%time%ttime - state%time%ttime_save ) /state%time%FinTime)
          r_est = sqrt(state%T_estim(resST) )

          ! estimated error and its tolerance for unit length interval
          !r_tol = state%space%adapt%tol_max  &
          !     * sqrt(( state%time%ttime - state%time%ttime_save ) * state%time%FinTime )
          !r_est = sqrt(state%T_estim(resST) )

          !r_tol = state%space%adapt%tol_max  ! L^\infty case ??????????

       else ! BDF time discretization
          call RezidErrorEstimates( .false.  , .true.)
          do i=resA, resST  ! i.e., resA, resS, resT, resST
             state%estim(i, 1) = state%estim(i, 1) + state%L_estim(i)**2 *state%time%tau(1)
             state%T_estim(i)  = state%T_estim(i)  + state%L_estim(i)**2 *state%time%tau(1)
          enddo

          ! estimated error and its tolerance
          r_tol = state%space%adapt%tol_max /state%time%FinTime**0.5
          r_est =( state%T_estim(resST)/( state%time%ttime - state%time%ttime_save ))**0.5
          if( state%time%iter_loc <= 1) print*,'Not tested FWE'

       endif

       ! screen output
       if( finish ) then
          print*
          write(*,'(a40, l2, es12.4, a3, es12.4, a11, es12.4)') &
               '# hp-iso adapt ?,estim_tot, tol:', &
               state%space%adaptation, r_est, ' < ', r_tol, ', ratio = ', r_est/r_tol
          !write(*,'(a20,4es12.4)') 'Time one-step estims', state%L_estim(resA:resST)
          write(*,'(a20,4es12.4)') 'Time local estims', sqrt(state%T_estim(resA:resST))
          write(*,'(a20,4es12.4)') 'Time global estims', sqrt(state%estim(resA:resST, 1))


          if(state%time%iter== 1) r_est =2 * r_tol   ! enforce refinement on the first level

          ! error estimate is bigger than the tolerance, we need a remeshing
          state%space%adapt%stop_adaptation = 1

          if( (r_est >  fac_ref_max * r_tol .or. r_est <   fac_ref_min * r_tol ) &
               .and. state%space%adapt%max_adapt_level > 0) then

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

          call ComputeRTNstEstim( state%loc_estim(1:max_eta, 1:ndim) , state%time%Qnum )
          state%estim(:, :) = state%estim(:, :) + state%loc_estim(1:max_eta,:)

          ! EXPORT USED TO SOLVE THE DUAL 3D PROBLEM WITH FENICS
          ! how many inner values in each elem is computed (degree of Lagrange nodes )
          spaceLevels = 1 !4* (state%space%deg + 1 )           ! number of equidistantly distributed time levels + 1 , i.e. # of nodes
          timeLevels  = 2! 4* (state%time%deg  + 1 ) + 1       ! +1 error in the setting, i.e. 5 ~ 4 intervals
          call exportRTNst( grid, spaceLevels, timeLevels )

          print*, 'Total RTN estimate:', sqrt( state%estim( RTNeta, 1) + state%estim( RTNjump, 1 ) ) , &
          sqrt( state%estim( RTNeta, 1) ) / sqrt( state%estim( RTNeta, 1) + state%estim( RTNjump, 1 ) ) , '%'
          print*, 'Upper bound(norm):', sqrt( state%estim( RTNfluxnorm, 1) + state%estim( RTNjump, 1 ) ) , &
          sqrt( state%estim( RTNfluxnorm, 1) ) / sqrt( state%estim( RTNfluxnorm, 1) + state%estim( RTNjump, 1 ) ) , '%'
          !write(debug,*) 'Probably dont use RTNall since it has to be squarerooted in WriteOutputErrors'
          !state%estim(RTNall, :) = sqrt( state%estim(RTNeta,:) ) + state%estim(RTNjump,:)

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

    elseif ( state%space%estim_space == 'DWR' ) then
       stop 'ErrorEstims_time_dependent called. NOT IMPLEMENTED FOR DWR METHOD'

    else
       stop ' UNKNOWN TYPE in ErrorEstims in compuateAD.f90 for TIME-DEPENDENT'
    endif


    ! too large error
    if( r_est <=  fac_ref_max * r_tol )  then
       state%space%adapt%stop_adaptation = 1

       ! too small error
       if( r_est < fac_ref_min * r_tol )  then
          state%space%adapt%stop_adaptation = 3
       endif

    endif

    call state%cpuTime%addEstimTime()

  end subroutine ErrorEstims_time_dependent


  !> several type of error estimates
  !> TIME-INDEPENDET COMPUTATIONS
  subroutine ErrorEstims_time_independent( )
    real ::  r_tol, r_est, Lq, r_estOLD
    real :: l_norm
    real :: cF, xF
    integer :: i, tQnum
    real :: t1, t2

    r_tol = 0.
    r_est = 0.


    !print*, 'ErrorEstims_time_independent called !!!'
    state%space%adapt%stop_adaptation = 0

    if( state%space%estim_space == 'RES' ) then ! residual error estimates

       if( state%time%disc_time == 'STDG') then
          do i=resA, resSr  ! i.e., resA, resS, resT, resST, resSr
             !state%estim(i, 1) = state%estim(i, 1) + state%L_estim(i)**2 !*state%time%tau(1)
             !state%T_estim(i)  = state%T_estim(i)  + state%L_estim(i)**2 !*state%time%tau(1)

             state%estim(i, 1) = state%L_estim(i)**2 /state%time%tau(1) !**2
             state%T_estim(i)  = state%L_estim(i)**2 /state%time%tau(1) !**2
          enddo


       else! BDF time discretization

          call RezidErrorEstimates( .false.  , .true.)

          do i=resA, resSr  ! i.e. resA, resS, resT, resST, resSr
             state%estim(i, 1) = state%L_estim(i)**2
          enddo
       endif

       ! estimated error and its tolerance
       r_est = sqrt(state%estim(resS,1))
       r_tol = state%space%adapt%tol_max

       if(state%space%adapt%adapt_type == 'Ahp' .or. state%space%adapt%adapt_type == 'Ihp')  then
          ! marking elements for the refinement
          !call Set_Mesh_regularity( )

          ! isotropic hp refinement using the maximal top refinement, defines new metrices
          !call IsotropicMetric_SimpleOrdering( )


          ! print*,'ATTENTION, ComputeLocalNeumannEstim overwrites state%estim(:,:) '
          ! do i =1, grid%nelem
          !    grid%elem(i)%errL8 = grid%elem(i)%eta(resS, 1)
          ! enddo
          ! call ComputeLocalNeumannEstim( )
          ! !write(*,'(a8, 50es12.4)') 'estimsP:', sqrt( abs(state%estim(1:max_eta, 1)) )
          ! do i =1, grid%nelem
          !    grid%elem(i)%eta(resS, 1) = grid%elem(i)%errL8
          ! enddo


          !print*,'Vertex based HO reconstruction'
          !call ComputeHO_Recovery( )
          !call ComputeHO_LocalProblems( )

          ! setting of the metric
          call AnisotInterpolEstimates( )
          !call AnisotResidEstimates( )
          state%estim(11, 1) = r_est**2
          state%estim(2, 1) = r_est**2

       endif

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
       !print*,'Higher-order _reconstruction subroutines'

       !print*,'Higher-order _reconstruction subroutines AnisotInterpolEstimates( )'
       !call AnisotInterpolEstimates( )

       !print*,'Higher-order _reconstruction subroutines'
       !!call ComputeHigherOrderEstims( )
       call ComputeHigherOrderEstims_OPT( )


       ! print*,'Vertex based HO reconstruction'
       !call ComputeHO_LocalProblems( )

       !print*,"'! estimated error and its tolerance"
       r_est = state%estim(HO_estim_H1_p2,1)**0.5
       r_tol = state%space%adapt%tol_max



    elseif( state%space%estim_space == 'inter') then ! interpolation error estimates
       Lq = state%space%adapt%Lq

       r_estOLD = state%err(interL8)                      ! L^\infinity- norm
       if(Lq >=  1.) r_estOLD = state%err(interLq)        ! L^q - norm
       if(Lq <= -1.) r_estOLD = state%err(interH1)        ! H^1 - semi-norm

       call AnisotInterpolEstimates( )

       ! only for settings arrays state%estim(:, HO_rec_p2_*)
       !call ComputeHigherOrderEstims_OPT( )

       ! print*,'ATTENTION, RezidErrorEstimates overwrites state%estim(:,:) 87d3jd3'
       ! call RezidErrorEstimates( .false.  , .true.)
       ! do i=resA, resSr  ! i.e. resA, resS, resT, resST, resSr
       !    state%estim(i, 1) = state%L_estim(i)**2
       ! enddo
       ! r_est = state%estim(2, 1)


       ! print*,'ATTENTION, ComputeLocalNeumannEstim overwrites state%estim(:,:) 87d3jd3'
       ! do i =1, grid%nelem
       !    grid%elem(i)%errL8 = grid%elem(i)%eta(resS, 1)
       ! enddo
       ! call ComputeLocalNeumannEstim( )
       ! !write(*,'(a8, 50es12.4)') 'estimsP:', sqrt( abs(state%estim(1:max_eta, 1)) )
       ! do i =1, grid%nelem
       !    grid%elem(i)%eta(resS, 1) = grid%elem(i)%errL8
       ! enddo
       ! state%estim(2, 1) = r_est

       ! print*, 'RES / pNEU error estimate:', sqrt(state%estim(resS,1)), sqrt(state%estim(P_tot,1))

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
          if( state%space%adapt%adapt_type == 'IhpOLD' ) then ! not used for Ihp
             cF = 0.95
             xF = 0.9
             !if(r_est > cF*r_estOLD .and. state%space%adapt%adapt_level > 1) then
             !   state%space%adapt%tol_min = state%space%adapt%tol_min*xF
             !   print*,'Security factor switch ON', state%space%adapt%tol_min, r_estOLD, r_est
             !endif

             !print*,'REMOVe D'
             !call ComputeHigherOrderEstims_OPT( )

             call IsotropicMetric( )

             !call IsotropicMetric_SimpleOrdering( )   ! variant adapting only top 10% of elements

          elseif( state%space%adapt%adapt_type == 'HG' ) then
             !print*,' HG-refinement'

          elseif(state%space%adapt%adapt_type == 'Ahp' .or. state%space%adapt%adapt_type == 'Ihp')  then

             ! setting of the metric
             call AnisotInterpolEstimates( )


          else
             stop 'No verified other variants in  ErrorEstims_time_independent'
          endif
       endif

    elseif( state%space%estim_space == 'DWR' ) then
       ! for aDWR the estimates are already computed
       if (state%nlSolver%non_alg_stop /= 'aDWR') then
          call DWR%clean()
          call PrepareDualProblem( DWR )
          call SolveDualProblem_stat( DWR, grid, convfile )
          DWR%dualProblem_computed = .true.
          call computeDWRestimates( DWR, grid )
       endif

       !      print*,'------------  AASD3 ' , DWR%J%coord_i
       r_est = sqrt( abs( state%estim( DWR%eta_index, DWR%J%coord_i ) ))
       r_tol = state%space%adapt%tol_max

       if(r_est > r_tol ) then
          if( state%space%adapt%adapt_type == 'Ihp' .or.  state%space%adapt%adapt_type == 'Ahp' ) then
             call AnisotInterpolEstimates( )
             !call IsotropicMetric( )
             !call IsotropicMetric_SimpleOrdering( )   ! variant adapting only top 10% of elements

          elseif( state%space%adapt%adapt_type == 'HG' ) then
             !             print*,' HG-refinement'
          else
             stop 'No verified other variants in  ErrorEstims_time_independent'
          endif

      endif
      !print*, 'now DWR%del_DWR was called but it should be called only when the computation is done.'
      !call DWR%delete()

    ! estimate based on the Ritz reconstruction
    elseif ( state%space%estim_space == 'ERRp' ) then

       call ComputeERRpEstimate( grid, convfile )
       r_est = sqrt( sum(state%estim(resS, 1:ndim)) )
       r_tol = state%space%adapt%tol_max

    else
       stop ' UNKNOWN TYPE in ErrorEstims in computeAD.f90 for TIME-INDEPENDENT'
    endif

    ! Screen output
    print*,
    write(*,'(a42, a6, es12.4, a3, es12.4, a2)') &
         '## Error estimates: method, estim, tol:  ', &
         state%space%estim_space, r_est, ' < ', r_tol,'?'

    if( r_est <= r_tol )  state%space%adapt%stop_adaptation = 1
    !print*,'#####.. check2', state%space%adapt%stop_adaptation

  end subroutine ErrorEstims_time_independent


!  subroutine ErrorEstimsOLD( )
!    integer :: i, tQnum
!    real :: t1, t2
!
!    call cpu_time(t1)
!
!    print*, 'Calling ErrorEstims'
!
!    !if(state%modelName == 'scalar' .or.state%modelName == '2eqs')  then
!    if(state%space%estim_space == 'RTN') then
!       call ComputeApostEstim( )
!       ! write(*,'(a28,6es12.4)')'state%estim(1:Osc,1:ndim):', state%estim(1:Osc,1)
!
!    elseif( state%space%estim_space == 'DUA') then
!       call ComputeDualEstim( )
!       !write(*,'(a28,6es12.4)') 'state%estim(1:Osc,1:ndim):', state%estim(1:Osc,1)
!!
!!    elseif( state%space%estim_space == 'RTNst' ) then
!!      ! number of time nodes for integration
!!       tQnum = state%time%max_Tdof ! number of integration nodes
!!       call ComputeRTNstEstim( state%estim(1:max_eta, 1:ndim) , tQnum )
!
!    elseif( state%space%estim_space == 'HEL') then
!       call ComputeHelmholtzEstim( )
!       ! write(*,'(a28,6es12.4)')'state%estim(1:Osc,1:ndim):', state%estim(1:Osc,1)
!
!    elseif( state%space%estim_space == 'HO_rec') then
!       !!call ComputeHigherOrderEstims( )
!       call ComputeHigherOrderEstims_OPT( )
!
!    elseif( state%space%estim_space == 'inter') then
!       call AnisotInterpolEstimates( )
!
!    elseif( state%space%estim_space == 'pNeu') then
!      call ComputeLocalNeumannEstim( )
!
!    elseif( state%space%estim_space == 'DWR' ) then
!         print*, 'subroutine ErrorEstims in computeAD. Still not implemented for DWR method. Using RES.'
!
!
!    elseif( state%space%estim_space == 'RES' .or. state%space%estim_space == 'DWR' ) then ! .or. state%space%estim_space == 'ANI' &
!       !.or. state%space%estim_space == 'Ahp' ) then
!
!       !print*, ' else Error estimates for RES already computed in PerformOneTimeStep'
!       !if(state%time%time_method /= 'I'.or. state%nlSolver%non_alg_stop /= 'aRES' .or. &
!       !     (state%time%tau_choice == 'adapt' .and.  state%time%estim_time /= 'tRES') ) then
!       !   if( state%time%disc_time /= 'STDG') then
!       !      call RezidErrorEstimates( .false. , .true.)
!       !   else
!       ! print*,'computeAD - errorEstims calling RezidErrorEst (false)!! NOT IMPLEMENTED'
!       ! for  state%time%disc_time == 'STDG'already computed
!       if( state%time%disc_time /= 'STDG')  &
!            call RezidErrorEstimates( .false.  , .true.)
!       !   endif
!       !endif
!
!       if(state%time_dependent) then
!
!          if( state%time%disc_time /= 'STDG') then
!             do i=resA, resSr  ! i.e., resA, resS, resT, resST
!                state%estim(i, 1) = state%estim(i, 1) + state%L_estim(i)**2 *state%time%tau(1)
!                state%T_estim(i)  = state%T_estim(i)  + state%L_estim(i)**2 *state%time%tau(1)
!             enddo
!          else
!             do i=resA, resSr  ! i.e., resA, resS, resT, resST
!                state%estim(i, 1) = state%estim(i, 1) + state%L_estim(i)**2 !*state%time%tau(1)
!                state%T_estim(i)  = state%T_estim(i)  + state%L_estim(i)**2 !*state%time%tau(1)
!             enddo
!          endif
!       else
!          if( state%time%disc_time /= 'STDG') then
!             do i=resA, resSr  ! i.e. resA, resS, resT, resST
!                state%estim(i, 1) = state%L_estim(i)**2
!             enddo
!          else
!             do i=resA, resSr  ! i.e., resA, resS, resT, resST
!                !state%estim(i, 1) = state%estim(i, 1) + state%L_estim(i)**2 !*state%time%tau(1)
!                !state%T_estim(i)  = state%T_estim(i)  + state%L_estim(i)**2 !*state%time%tau(1)
!
!                state%estim(i, 1) = state%L_estim(i)**2 /state%time%tau(1)**2
!                state%T_estim(i)  = state%L_estim(i)**2 /state%time%tau(1)**2
!             enddo
!          endif
!       endif
!
!       !open(99, file='estims', status='UNKNOWN', position = 'append')
!       !write(99,'(a6,i5,30es12.4)') 'ETAs:', state%time%iter, &
!       !     state%time%tau(1), state%L_estim(1:4), state%estim(1:4, 1)**0.5
!       !close(99)
!    else
!       stop ' UNKNOWN TYPE in ErrorEstims in compuateAD.f90'
!
!    endif
!
!    call cpu_time(t2)
!    state%CPU_estim2 = state%CPU_estim2 + t2 - t1
!
!
!  end subroutine ErrorEstimsOLD

  !> Perform some time steps on the given grid
  subroutine SolveProblemAD(convfile, finish)
    use matrix_oper_int
    use pedes_averaging
    character(len=*), intent(in) :: convfile
    logical, intent(out) :: finish
    class(element), pointer :: elem
    integer :: iter, istep, i, j, k
    real :: tt, t_sta, time_prepare, time_solve, time_estim, diff, diff1
    logical :: refuseStep, loc_implicitly!!, CrankNicolson
    logical :: pedes_use_default_initc, do_not_repeat_time_step
    real :: time, lost

    call state%cpuTime%startPrepareTime()


    pedes_use_default_initc = .true.   ! only for the pedestriang flow with eikonal equations
    if(state%modelName == 'pedes' ) call  pedes_w_bar_alloc(grid%nelem)

    !print *,' call Clear_Elem_Estim_locL( )'
    call Clear_Elem_Estim_locL( )

    refuseStep = .false.
    state%linSolver%residuum = 0.

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

    call state%cpuTime%addPrepareTime()

    do iter = 1,state%time%maxiter

       state%time%iter_loc = iter

       !print*,'PrepareNewTimeStepAD(iter)'
       call PrepareNewTimeStepAD(iter)
       !print*,'PrepareNewTimeStepAD(iter)', state%time%tau(1)

       if(state%modelName == 'pedes' ) then
          call PedestrianEikonalEquation(pedes_use_default_initc )
          !call PlotEikonalVelocityVectors(1000+state%time%iter, 2000+state%time%iter)
       endif

       state%linSolver%iter = 0

       do_not_repeat_time_step = .false.

       !! cycles repeated if the time step is refused
       do istep = 1, 5 ! 10
          state%time%ctime = state%time%ttime + state%time%tau(1)

          !elem => grid%elem(150)
          !do k=1,ndim
          !   do i=1, elem%Tdof
          !      write(*,'(a8,3i5,200es12.4)') 'ELEM:',elem%i,k,i, elem%wST(k, 1:elem%dof, i)
          !   enddo
          !enddo

          if( state%time%disc_time /= 'STDG') then
             !print*,' perform one time step', state%time%tau(:)

             call PerformOneTimeStep( time_prepare, time_solve, time_estim )

          else  ! STDGM

             if( state%modelName == 'porous' ) then
                call PerformOneSTDGMstep_Andrerson( time_prepare, time_solve, time_estim )
                !call PerformOneSTDGMstep( time_prepare, time_solve, time_estim )
                !call PerformOneSTDGMstep_OLD( time_prepare, time_solve, time_estim )

             else
                call PerformOneSTDGMstep()
                ! OLD - should be removed in future
                time_prepare = state%cpuTime%prepare
                time_solve   = state%cpuTime%solve
                time_estim   = state%cpuTime%estim
             endif

          endif

          !!!if(ndim > 1) call AvoidNonphysicalValues( )
          call state%cpuTime%startPrepareTime()

          ! adding  of the elem%estims to elem%estimsL
          if(state%time_dependent) then
             do i=1,grid%nelem
                grid%elem(i)%estim_locL = sqrt( grid%elem(i)%estim_locL**2  &
                     + grid%elem(i)%eta(resST, 1)**2)
             enddo
          endif

          if(state%time%tau_choice == 'fixed' ) then

          !elseif( (state%linSolver%lin_solver_not_conv > 0 .or. state%linSolver%residuum > 1.) &
          !     .or. (  state%nlSolver%iter >=  state%nlSolver%max_iter) )then

             ! if( state%linSolver%lin_solver_not_conv > 0 .and. state%linSolver%residuum > 1.) then
             !    state%time%tau_new = state%time%tau(1) * 0.2
             ! else
             !    state%time%tau_new = state%time%tau(1) * 0.75
             ! endif

             ! state%time%tau(1) = state%time%tau_new
             ! write(*,'(a20,a12,es12.4,a20, es12.4, a10, 2i6 )') &
             !      'Refused time step:', 'residuum =', state%linSolver%residuum, &
             !      ', NEW time step = ', state%time%tau_new,  &
             !      ', iters:',state%nlSolver%iter,  state%nlSolver%max_iter
             ! refuseStep = .true.

          else
             !if(state%time%tau_choice /= 'fixed') then
             if(state%time%tau_choice == 'adapt') then

                if(state%time%keep_tau) then
                   state%time%tau_new = state%time%tau(1)
                   !print*,' adaptive choice of the time step',state%L_estim(resS), state%L_estim(resT)
                   !print*,' Adaptive choice of the time step' ,  state%time%tau_new
                   call ProposeNewTimeStep(refuseStep )  ! for fixed time only error estim
                   !print*,' adaptive choice of the time step' ,  state%time%tau_new


                else
                   state%time%tau_new = state%time%tau(1) * 0.5
                   state%time%tau(1) = state%time%tau_new
                   !  if(.not. refuseStep) goto 90
                endif
             endif
          endif

          if(.not. refuseStep .or. do_not_repeat_time_step) goto 90

          if (state%modelName == 'pedes') then
             if( abs( (state%time%tau_new  - state%time%tau(1)) /  state%time%tau(1)) < 1E-5) &
                  do_not_repeat_time_step = .true.
          endif


          ! recomputation of the elem%estims
          if(state%time_dependent) then
             do i=1,grid%nelem
                grid%elem(i)%estim_locL = sqrt( grid%elem(i)%estim_locL**2 &
                     -grid%elem(i)%eta(resST, 1)**2)
             enddo
          endif

          !do i=1,grid%nelem
          !   call PlotElem_D_Function3D(10+iter*10 + istep, grid%elem(i),  grid%elem(i)%dof, &
          !        grid%elem(i)%w(0, 1:grid%elem(i)%dof) )
          !enddo
          ! restoring of the the solution from the values before time refusing
          if( state%time%disc_time == 'STDG') call RestoreSTDGMsolution( )

          call state%cpuTime%addPrepareTime()

       enddo ! do istep = 1,5


90     continue
       call state%cpuTime%addPrepareTime()
       call state%cpuTime%startPrepareTime()

       if( state%time%disc_time == 'STDG') call DeallocateStoredSolution( )

       !print*,' call JumpsEvaluation( ) '
       call JumpsEvaluation( )

       !print* ,' one time step was perfomed, we go its evaluation'
       call PassToNewTimeStep( )

       !!state%time%cn = CrankNicolson

       state%time%ctime = state%time%ttime  ! ctime was changed, we keep the actual value for Estims
       !  print*, 'calling ErrorEstims in SolveProblemAD after passToNewTimeStep'
       !write(*,'(a10,3es12.4)') 'Errors:',state%err(L2), state%err(H1), state%err(H1_discrete)

       finish = finishComputationProcess()
       call WriteOutputScreen(iter, time_prepare, time_solve, time_estim)

       ! pedestrian flow, the domain is empty
       if (state%modelName == 'pedes' )  call Pedes_Empty_domain (finish)


       if(finish) then
         if(ndim >= 4) call  ComputeSkinFriction()
       endif

       call state%cpuTime%addPrepareTime()

       if(state%time_dependent) call ErrorEstims( finish  )

       ! convfile, progress output
       call WriteOutputFiles (convfile, iter, t_sta, time_prepare, time_solve, time_estim )

70     continue

       if(finish ) goto 100

    enddo  ! go to next time step

100 continue
    state%time_AD_end = state%time%ttime

    if( state%MGsolver ) call DeInitMG( )

    if(state%modelName == 'pedes' ) call  pedes_w_bar_DEalloc( )

!    if( (state%modelName == 'scalar' .or. state%modelName == '2eqs') .and. state%type_IC .eq. 5) then
!       !print*,'Output for LEVEL SET METHOD in file "data.vofx"'
!       call WriteOutput_LevelSet()
!    endif

    !print*,'##### END of SolveProblemAD'

  end subroutine SolveProblemAD



  !> call the approxiate mesh adaptation technique, create new mesh: grid
  subroutine AdaptationAD(  )
    logical :: identical_AMA_grids
    logical :: metric
    integer :: indexi, i
    real ::  r_tol, r_est, regularity
    class(element), pointer :: elem
    real :: tt , t2

    call cpu_time(tt)
    ! new structure for cpu time measuring
    call state%cpuTime%startAdaptTime( )

    ! ONLY TEST
 !    do i = 1, grid%nelem
 !       elem => grid%elem(i)
 !       !call Set_Elem_Coeffs_Decay(elem)
 !       !call Set_Elem_Regul_Estim_Decay(elem)

 !       call ElementEdgeJumpsAllProj(elem, regularity)

 !    enddo
 !    ! ONLY TEST END

   ! test of the regularity


!    print*,'______________________________AdaptationAD____________________________________'
    print*
    write(*,'(a40, i4,a18)') '___________________________AdaptationAD:', state%space%adapt%adapt_level+1, '__________________'

    !call WriteResults('Gsol.bak')
    if(state%space%adapt%adapt_type == 'Ahp' .or. state%space%adapt%adapt_type == 'Ihp') then ! ANGENER based on HO interpol

       metric = .false.   !metric = .true. ==> metric computed by ANGENER F77
       identical_AMA_grids = .false.

       ! already  done in compute1
       !!if(.not. metric) call AnisotInterpolEstimates( )

       !   write(*, *) ' ###  Stopping criterion was not achieved'

        !write(*, *) ' DEVELOPMENT OF NEW VARIANT OF ANGENER'
        !call AdaptMesh_AMA( )     ! f95 given by Gabriel Patho
        call AdaptMesh_Angener(metric, identical_AMA_grids)    ! ANGENER 77 - variant

       if(identical_AMA_grids) state%space%adapt%stop_adaptation = 11

    elseif ( state%space%adapt%adapt_type == 'HG' .or. state%space%adapt%adapt_type == 'RG') then
       ! !!call PlotMesh(grid, 'mesh')

       if (state%space%estim_space == 'DWR') then
          print*,'WE USE MarkTOPElements (DWR) !!!!!!!!!!!!!!!!'
          print*

          call MarkTopElements( )
!          print*,'WE USE MarkElements not MarkTOPElements (e.g., for SISC) !!!!!!!!!!!!!!!!'
!          print*
!
!          call MarkElements( )

       else
          print*
          print*,'WE USE MarkElements not MarkTOPElements (e.g., for SISC) !!!!!!!!!!!!!!!!'
          print*

          call MarkElements( )
       endif


       !call MarkTopElements( )

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

    call cpu_time(t2)
    ! new version
    call state%cpuTime%addAdaptTime( )
    state%CPU_adapt = state%CPU_adapt + t2 - tt

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

end module computeAD_oper
