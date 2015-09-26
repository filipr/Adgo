program adgo
  ! use oop
  use compute_mod
  use data_mod
  use define_state
  use main_data
  use mesh_oper, only : ReadCurvedBoundary, SeekCurvedBoundary
  use mesh_mod
  use model_mod
  use problem_oper
  use solver_mod
  !use solve_problem



  implicit none

  real :: t2
!  real, pointer :: aa, bb
!
!  allocate(aa)
!  aa = 3.5
!  print*,  aa
!
!  bb => aa
!
!  print* , bb
!  deallocate(aa)
!  bb = 0.0
!  print*, bb


!  print*, 'HELLO WORLD'
!
!  stop 'POZOR na copyElement'


  open(debug, file = debug_file, action="write", status="replace" )

  call cpu_time(state%start_time)

  call readMainData()

  ! PREPROCESSING
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

!to prepare problem
!  call grid%computeGeometry()
!
!  write(debug,*) 'H and domain volume:', state%space%h , state%space%domain_volume
!  call state%space%copyGeometry( grid%h , grid%domain_volume )

 ! call mesh%prepare()

  write(debug,*) 'state%space%QDEG(:,1) ' , state%space%Qdeg(:, 2)
  write(debug,*) 'state%space%QDEG(:,2) ' , state%space%Qdeg(:, 1)

  write(debug,*) 'PrepareProblem + ReprepareProblem - PUT TOGETHER'
  write(debug,*) 'Prepare problem is called only once'
  call PrepareProblem(grid, rsolfile)

  call cpu_time(t2)
  print*,'# Problem prepared for computation within ',t2-state%start_time, ' s'
  print*,'# ------------------------------------------------------'

  ! this is called only here and never again
  call InitSolveProblem(convfile)
 !print*,'InitSolveProblem -- done'

  write(debug,*) 'main - where to put Set_R_s ??? NOT USED?'
  !state%Set_R_s = 'Set_R_s_scalar'


  !call Reconstruction_test(state%space )
  !stop 'stopped after rconstruction test'


  call ComputeADGo( convfile, solfile, state%space%adaptation, state%time_dependent )

  !call Reconstruction_test2D( )
  !stop 'stopped after reconstruction test'

  call cpu_time(t2)

  print*,'# ADGFEM finished after ',t2-state%start_time, ' s'

  call CleanMemory ( )

  close(debug)

end program
