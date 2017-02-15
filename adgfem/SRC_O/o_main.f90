program adgo
  use compute_mod
  use cpu_time_mod
  use data_mod
  use define_state

!  use elemental_mod

  use main_data
  use mesh_oper, only : ReadCurvedBoundary, SeekCurvedBoundary
  use mesh_mod
  use model_mod
  use problem_oper
  use solver_mod

  implicit none

  real :: t2

  open(debug, file = debug_file, action="write", status="UNKNOWN" )

  call cpu_time(state%start_time)
  ! new structure for CPU times
  allocate( state%cpuTime )
  call state%cpuTime%initCpuTime()
  call state%cpuTime%startPrepareTime()

  call readMainData()

  !print*,' PREPROCESSING'
  ! general association, fixed during all computations and adaptations
  call state%init_Problem()

  ! subroutine initMesh is called in readMesh
  call grid%read( gridfile )

  call grid%plot('mesh0')
  if (nbDim == 3) then
   stop 'Mesh is not implemented for 3D'
    !call PlotMesh3D(grid)
  endif
  print*,'# Mesh plotted'

  call grid%seekNeighbours()

  if (nbDim==2) then
     call ReadCurvedBoundary(prof_file)
     call SeekCurvedBoundary(grid)
     !call grid%SeekCurvedBoundary()
     !  call SplineCurvedBoundary( )
     !print*,'Curved boundary found'
  endif

  ! Init the DWR technique
  if (state%space%estim_space == 'DWR') then
!    DWR = DWR_t( DWR%id, grid)
    call DWR%init( grid )
  endif

!to prepare problem
!  call grid%computeGeometry()
!
!  write(debug,*) 'H and domain volume:', state%space%h , state%space%domain_volume
!  call state%space%copyGeometry( grid%h , grid%domain_volume )

 ! call mesh%prepare()

  !write(debug,*) 'PrepareProblem + ReprepareProblem - PUT TOGETHER'
  !Prepare problem is called only once
  call PrepareProblem(grid, rsolfile)
  call state%cpuTime%addPrepareTime()

  print*,'# Problem prepared for computation within ', state%cpuTime%prepare, ' s'
  print*,'# ------------------------------------------------------'

  ! this is called only here and never again
  call InitSolveProblem(convfile)

  !call Reconstruction_test(state%space )
  !stop 'stopped after rconstruction test'

  call ComputeADGo( convfile, solfile, state%space%adaptation, state%time_dependent )

  !call Reconstruction_test2D( )
  !stop 'stopped after reconstruction test'

  call state%cpuTime%printCpuTimes()


  !write(*,'(a20, 16es12.4)') 'CPU_times:', t2-state%start_time, &
  !     state%CPU_prepare, state%CPU_solve, state%CPU_prepare+ state%CPU_solve, &
  !     state%CPU_estim, state%CPU_estim2, &
  !     state%CPU_estim/ (t2-state%start_time), state%CPU_estim2 / (t2-state%start_time)
  call CleanMemory ( )

  close(debug)

end program adgo
