!> module distinguishing different mode of the computations: steady/unsteady, adaptation/non-adaptation
module compute_mod
  use computeAD_oper

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


  end subroutine ComputeADGo

  !> steady as well as unsteady computation WITHOUT space adaptation
  subroutine Compute_NoAdapt( convfile, solfile )
    character(len=*), intent(in) :: convfile, solfile
    logical :: Aname

    Aname = .false.

    !do i=0, state%space%adapt%max_adapt_level
    !state%space%adapt%adapt_level = i

    if(state%time%OutTime > 0.) then
       !if (state%time%disc_time .ne. 'STDG') then
       call WriteProgressOutput( 'TS', Aname )
       !elseif (state%isol == 0) then
       !   call WriteProgressOutput( 'TS' )
       !endif
    else
       call WriteProgressOutput( 'T', Aname )
    endif

    !call WriteProgressOutput( 'E' , Aname)

    call SolveProblemAD(convfile)

    call WriteResults(solfile)

    !if(state%space%estim_space == 'pNeu') call WriteProgressOutput('A', Aname)

    if(ndim <= 2)call WriteProgressOutput( 'SEA', Aname )
    if(ndim >  2)call WriteProgressOutput( 'S' , Aname)
    !if(state%type_IC == 7) call WriteProgressOutput( 'E' )



  end subroutine Compute_NoAdapt

  !> STEADY computation WITH space adaptation
  subroutine Compute_Stat( convfile, solfile )
    character(len=*), intent(in) :: convfile, solfile
    integer :: i
    logical :: Aname

    Aname = .true.

    if ( state%space%estim_space /= 'DWR' ) then

       ! only for initialization
       if(state%space%estim_space == 'RES') then ! .or. state%space%estim_space == 'Ahp') then
          state%time%deg_actual = 1
          call RezidErrorEstimates( .true., .false. )
          call JumpsEvaluation( )
       endif

       do i= 0,  state%space%adapt%max_adapt_level
          state%space%adapt%adapt_level = i
          ! !!call WriteProgressOutput( 'TS' )

          call SolveProblemAD(convfile)

          if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
             call WriteProgressOutput( 'STEA' , Aname)
          else
             call WriteProgressOutput( 'STA', Aname )
          endif

          ! !!call WriteProgressOutput('A')
          if ( state%space%estim_space == 'DWR' ) then
             call PrepareDualProblem( )
             print*, 'Solve dualProblem not implemented. Stop.'
             stop

             call SolveDualProblem( convfile )
          endif


          if(i < state%space%adapt%max_adapt_level .and.  state%space%adapt%stop_adaptation==0 ) &
               call AdaptationAD( )

          !! print*,' SMAZ edtkcehcek'
          !! state%space%adapt%adapt_level = state%space%adapt%adapt_level + 1
          !! call WriteProgressOutput( 'ST', Aname )
          !! state%space%adapt%adapt_level = state%space%adapt%adapt_level - 1
          !! SMAZ END

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
    else ! DWR error estimation

       print*, 'call initDualProblem commented'
       ! call initDualProblem()

       do i= 0,  state%space%adapt%max_adapt_level
          state%space%adapt%adapt_level = i
          ! !!call WriteProgressOutput( 'TS' )

          call SolveProblemAD(convfile)

          if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
             call WriteProgressOutput( 'STEA' , Aname)
          else
             call WriteProgressOutput( 'STA', Aname )
          endif


          call PrepareDualProblem( )
          print*, 'Solve dualProblem not implemented. Stop.'
          stop

          call SolveDualProblem( convfile )



          if(i < state%space%adapt%max_adapt_level .and.  state%space%adapt%stop_adaptation==0 ) &
               call AdaptationAD( )

          if(state%space%adapt%stop_adaptation > 0) goto 20

       enddo! state%space%adapt%max_adapt_level

   20  continue

       if(state%space%adapt%stop_adaptation == 1) then
          write(*,*) ' # Error criterion achieved, no further adaptation'
       elseif(state%space%adapt%stop_adaptation == 2) then
          write(*,*) ' # No element marked for adaptation'
       elseif(state%space%adapt%stop_adaptation == 11) then
          ! written in anisot.f90
          !write(*,*) 'The new hp-grid is (almost) identical with the previous ones'
       endif

    endif

  end subroutine Compute_Stat


  !> UNSTEADY computation WITH space adaptation
  subroutine Compute_NonStat( convfile, solfile )
    character(len=*), intent(in) :: convfile, solfile
    character(len=50) :: newgridfile
    integer :: i, i0, j, max_recomput_back
    character(len=1) :: ch1
    integer :: text_size, is, num_size
    logical :: Aname


    Aname = .true.

    ! only for initialization
    if(state%space%estim_space == 'RES' .or. state%space%estim_space == 'Ahp') then
       state%time%deg_actual = 1

       if( state%time%disc_time /= 'STDG') then
          call RezidErrorEstimates( .true., .false. )
       endif
       call JumpsEvaluation( )
    endif


    state%space%adapt%adapt_level = 0
    if(state%tri_solA)      call WriteProgressOutput( 'ST', Aname)
    if(state%time%OutTime > 0.)  call WriteProgressOutput( 'ST', .false. )
    !print*,'SW1'
    !pause

    max_recomput_back =  8  !!5
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
          call SolveProblemAD(convfile)

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


          if(state%space%adapt%stop_adaptation == 0) call AdaptationAD( )

          !print*,'after adaptation '
          !pause

          !!stop

          !state%isol = state%isol + 1
          !call WriteProgressOutput( 'ST', Aname )

          !print*,'dtge5sred',state%space%adapt%stop_adaptation, j, max_recomput_back, &
          !     state%space%adapt%adapt_level,  state%tri_solA

          if(state%space%adapt%stop_adaptation <= 0 .and. j < max_recomput_back) then
             call RecomputeBackMesh(  )

             ! saving of the initial condition on the new mesh
             if(state%space%adapt%adapt_level == 0 .and. state%tri_solA) &
                  call WriteProgressOutput( 'ST', Aname )
             if(state%space%adapt%adapt_level == 0 .and. state%time%OutTime > 0.)  &
                  call WriteProgressOutput( 'ST', .false. )
             state%time%keep_tau = .false.
          else
             state%errSTnorm(L8L2) = max(state%errSTnorm(L8L2),state%errSTloc(L8L2))
             print*
             print*,'Sucesfull multi-time step'
             print*
             state%err(interLq) = 1E+30
             state%err(interL8) = 1E+30

             state%space%adapt%adapt_level = state%space%adapt%adapt_level + 1

             if(state%space%gridS_allocated ) call RemoveOrigMesh( )
             state%time%keep_tau = .true.

             state%space%adapt%stop_adaptation = 1  ! only for a possible exit

             if(state%tri_solA) call WriteProgressOutput( 'ST', Aname )

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
