!> module distinguishing different mode of the computations: steady/unsteady, adaptation/non-adaptation
module compute_mod
  use computeAD_oper
  use dual_problem_mod

  implicit none

  public :: ComputeADGo
  public :: Compute_NoAdapt
  public :: Compute_Stat
  public :: Compute_NonStat

contains

  !> main procedure, compute the solution including adaptation
  !>
  !> three types of computation
  !> 1) NO ADAPTATION, steady as well as unsteady
  !> 2)    ADAPTATION, steady
  !> 3)    ADAPTATION, unsteady, reinterpolation back is necessary
  subroutine ComputeADGo( convfile, solfile, adaptation, time_dependent)
    character(len=*), intent(in) :: convfile, solfile
    logical, intent(in) :: adaptation
    logical, intent(in) :: time_dependent
    integer :: ifile
    real :: tt

    if(state%model%varying_time_term .and. state%time%disc_time /= 'STDG') then
       print*,'state%model%varying_time_term == .true. & state%time%disc_time /= "STDG" not tested'
       stop
    endif

    if ( ( .not. adaptation ) .or. state%space%adapt%max_adapt_level == 0 ) then
       print*, '# Compute_NoAdapt called'
       call Compute_NoAdapt( convfile, solfile )

    elseif ( .not. time_dependent ) then
       print*, '# Compute_Stat called'
       call Compute_Stat( convfile, solfile )

    else
       print*, '# Compute_NonStat called'
       call Compute_NonStat( convfile, solfile )
    endif

    ! saving of the DG solution in the coefficients of the orthogonal basis
    call WriteResults(solfile)
    call WriteMesh('FinalGrid.grid', grid)


    ! output of the physical coefficients
    if(.not. time_dependent .and. state%modelName == 'NSe') then
       call cpu_time(tt)

       ifile = 11
       open(ifile, file='AD_coeffs.dat', status='UNKNOWN', position = 'append')
       write(ifile,'(2i6,i9,  30es12.4)') &
                state%space%adapt%adapt_level, grid%nelem,  state%nsize,  tt-state%start_time, & !1..4
               state%space%adapt%tol_max, state%space%adapt%tol_min, &                          !5..6
               state%cDLM(state%time%iter, 1:5 ), state%err(interLq: interH1), &                !7..13
               sqrt(state%estim(1:4,1))                                                        !14..17
      close(ifile)

       open(ifile, file='AD_coeffs_loc.dat', status='UNKNOWN', position = 'append')
       write(ifile, '(x)')
       close(ifile)

    endif

  end subroutine ComputeADGo

  !> steady as well as unsteady computation WITHOUT space adaptation
  subroutine Compute_NoAdapt( convfile, solfile )
    character(len=*), intent(in) :: convfile, solfile
    logical :: Aname, finish

    Aname = .false.

    !do i=0, state%space%adapt%max_adapt_level
    !state%space%adapt%adapt_level = i

    if(state%time%OutTime > 0.) then
       !if (state%time%disc_time .ne. 'STDG') then
       call WriteProgressOutput( 'TS', Aname, 1. )
       !elseif (state%isol == 0) then
       !   call WriteProgressOutput( 'TS' )
       !endif
    else
       call WriteProgressOutput( 'T', Aname, 1. )
    endif

    !call WriteProgressOutput( 'E' , Aname)

    call SolveProblemAD(convfile, finish)

    if(.not. state%time_dependent) call ErrorEstims( finish )

    call WriteOutputsAfterSolveProblem('NoAdapt', Aname)

  end subroutine Compute_NoAdapt

  !> STEADY computation WITH space adaptation
  subroutine Compute_Stat( convfile, solfile )
    character(len=*), intent(in) :: convfile, solfile
    integer :: i, ifile
    logical :: Aname, finish
    real :: time, lost

    Aname = .true.

    ! only for initialization
    if(state%space%estim_space == 'RES') then ! .or. state%space%estim_space == 'Ahp') then
       state%time%deg_actual = 1
       call Clear_Elem_Estim_locL( )
       call RezidErrorEstimates( .false., .false. )
       call JumpsEvaluation( )
    endif

    do i= 0,  state%space%adapt%max_adapt_level
       state%space%adapt%adapt_level = i
       ! !!call WriteProgressOutput( 'TS' )

       call SolveProblemAD( convfile, finish )

       ! !!call WriteProgressOutput('A')
       call ErrorEstims( finish )

       ! write various outputs
       call WriteOutputsAfterSolveProblem('AdaptStat', Aname)

       if(i < state%space%adapt%max_adapt_level .and.  state%space%adapt%stop_adaptation==0 ) &
            call AdaptationAD( )

       if(state%space%adapt%stop_adaptation > 0) goto 10

    enddo! state%space%adapt%max_adapt_level

10  continue

    if(state%space%adapt%stop_adaptation == 1) then
       write(*,*) ' # Error criterion achieved, no further adaptation'

    elseif(state%space%adapt%stop_adaptation == 2) then
       write(*,*) ' # No element marked for adaptation'

    elseif(state%space%adapt%stop_adaptation == 11) then
       ! written in anisot.f90
       !write(*,*) 'The new hp-grid is (almost) identical with the previous ones'
    endif

  end subroutine Compute_Stat


  !> UNSTEADY computation WITH space adaptation
  subroutine Compute_NonStat( convfile, solfile )
    character(len=*), intent(in) :: convfile, solfile
    character(len=50) :: newgridfile
    integer :: i, i0, j, k, max_recomput_back
    character(len=1) :: ch1
    integer :: text_size, is, num_size
    logical :: Aname, finish
    real :: tt, t2
    logical :: deg_plus

    deg_plus = .true.
    Aname = .true.

    ! only for initialization
    if(state%space%estim_space == 'RES' .or. state%space%estim_space == 'Ahp') then
       state%time%deg_actual = 1

       if( state%time%disc_time == 'STDG') then
          call ComputeSTDGM_Terms(deg_plus )
       else
          call ComputeTerms(deg_plus )
       endif

       call Clear_Elem_Estim_locL( )
       call RezidErrorEstimates( .false., .false. )
       !endif
       call JumpsEvaluation( )
    endif


    state%space%adapt%adapt_level = 0
    if(state%tri_solA)      call WriteProgressOutput( 'ST', Aname, 1.)
    if(state%time%OutTime > 0.)  call WriteProgressOutput( 'ST', .false., 1. )
    !print*,'SW1'
    !pause

    max_recomput_back =  1 !3 !5  !!8  !!5
    do i= 1,  state%space%adapt%max_adapt_level + 1
       !state%space%adapt%adapt_level = i
       !call WriteProgressOutput( 'TS' )

       call SaveOrigMesh( )

       do j=1, max_recomput_back  ! inner loop within one time interval
          state%time%recompute_back = j
          ! REM TEST
          !state%isol = state%isol + 1
          !call WriteProgressOutput( 'STE' )

!!!!if(state%isol == 6) stop

          !print*,'SW2'
          !pause
          call SolveProblemAD(convfile, finish)
          ! Aname not used
          call WriteOutputsAfterSolveProblem('AdaptNonStat', Aname )

          ! do k=1, 5 !grid%nelem
          !    write(*,'(a8, 2i5,30es12.4)') 'estims:', grid%elem(k)%i,grid%elem(k)%dof, &
          !         grid%elem(k)%estim_loc, grid%elem(k)%eta(1:4, 1)
          ! enddo

          !!print* , "call WriteMesh('FinalGrid.grid', grid)"
          !!call WriteMesh('FinalGrid.grid', grid)


          ! REM TEST
          !state%isol = state%isol + 1
          !call WriteProgressOutput( 'STE' )

          !if(state%modelName == 'scalar' .or.state%modelName == '2eqs') call WriteProgressOutput( 'STE' )
          !if(ndim >  1) call WriteProgressOutput( 'ST' )
          ! !!call WriteProgressOutput('A')

          !print*,'before adaptation '
          !pause

!!! mesh from the initial grid is always recomputed
          !!if(i <= 2 .and. j == 1) state%space%adapt%stop_adaptation = 0


          call cpu_time(tt)
          if(state%space%adapt%stop_adaptation == 0) call AdaptationAD( )

          if(state%space%adapt%stop_adaptation == 3) then
             print*
             print*
             print*
             write(*,*) ' #  Error criterion achieved, too small error, further adaptation'
             print*
             print*
             print*
             call AdaptationAD( )
          endif

          ! already in AdaptationAD ? ! new structure in state%cpuTime
          call cpu_time(t2)
          state%CPU_adapt = state%CPU_adapt + t2 - tt

          !state%isol = state%isol + 1
          !call WriteProgressOutput( 'ST', Aname )

          !print*,'dtge5sred',state%space%adapt%stop_adaptation, j, max_recomput_back, &
          !     state%space%adapt%adapt_level,  state%tri_solA

          if(state%space%adapt%stop_adaptation <= 0 .and. j < max_recomput_back) then
             call RecomputeBackMesh(  )

             ! saving of the initial condition on the new mesh
             !if(state%space%adapt%adapt_level == 0 .and. state%tri_solA) &
             !     call WriteProgressOutput( 'STA', Aname, 1. )

             if(state%space%adapt%adapt_level == 0 .and. state%time%OutTime > 0.)  &
                  call WriteProgressOutput( 'ST', .false., 1. )

             state%time%keep_tau = .false.


             ! NOT YET TESTED
             !if(state%time%recompute_back == 4) then
             !   state%space%adapt%tol_max = state%space%adapt%tol_max * 0.75
             !   print*,' The tolerance state%space%adapt%tol_max was changed to:', &
             !        state%space%adapt%tol_max
             !endif

          else
             state%errSTnorm(L8L2) = max(state%errSTnorm(L8L2),state%errSTloc(L8L2))
             !print*
             !print*,'Sucesfull multi-time step'
             !print*
             state%err(interLq) = 1E+30
             state%err(interL8) = 1E+30

             state%space%adapt%adapt_level = state%space%adapt%adapt_level + 1

             if(state%space%gridS_allocated ) call RemoveOrigMesh( )
             state%time%keep_tau = .true.

             state%space%adapt%stop_adaptation = 1  ! only for a possible exit

             !if(state%tri_solA) call WriteProgressOutput( 'STA', Aname, 1. )

             goto 20
          endif

          !print*,'SW4'
          !pause

       enddo

20     continue

       !state%space%adapt%adapt_level = state%space%adapt%adapt_level + 1
       !if(state%time%OutTime == 0.) then
       !   state%isol = state%isol + 1
       !   call WriteProgressOutput( 'ST', Aname )
       !endif

       ! Final Time achieved
       if(state%time%ttime >= state%time%FinTime .and. state%space%adapt%stop_adaptation > 0) goto 30

    enddo! state%space%adapt%max_adapt_level
30  continue

  end subroutine Compute_NonStat

end module compute_mod
