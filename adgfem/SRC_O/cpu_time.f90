module cpu_time_mod

implicit none

type,public :: Cpu_t
   integer :: itime                         ! = 0, no CPU time outputs
   real :: start_time, end_time             ! initial CPU time and end
   real :: timeprn                          ! time for printing output
   real :: prepare, solve , estim           ! CPU time necessary for particular opeartions
   real :: estim2, adapt                    ! CPU time necessary for particular opeartions
   real :: constaint                        ! ?
   real :: t1                               ! local starting time
   real :: adaptStart, estimStart, prepareStart, solveStart ! starting times , clean in addTime
   logical :: prepareIsRunning
   logical :: adaptIsRunning
   logical :: estimIsRunning
   logical :: solveIsRunning

   ! MAYBE add also loc_times - for each adapt level ?

   contains

   procedure :: initCpuTime
   procedure :: addAdaptTime
   procedure :: startAdaptTime
   procedure :: addEstimTime
   procedure :: startEstimTime
   procedure :: addPrepareTime
   procedure :: startPrepareTime
   procedure :: addSolveTime
   procedure :: startSolveTime

   procedure :: cleanTimes
   procedure :: printCpuTimes
   procedure :: totalTime


end type cpu_t


contains

   subroutine initCpuTime(this)
      class( Cpu_t ) :: this
      real :: t

      call cpu_time(t)
      this%start_time = t

      this%prepare = 0.0
      this%adapt   = 0.0
      this%estim   = 0.0
      this%estim2  = 0.0
      this%solve   = 0.0
      this%end_time= 0.0

      this%adaptStart = 0.0
      this%estimStart = 0.0
      this%prepareStart = 0.0
      this%solveStart = 0.0

      this%adaptIsRunning = .false.
      this%estimIsRunning = .false.
      this%prepareIsRunning = .false.
      this%solveIsRunning = .false.

   end subroutine initCpuTime

   !> start clock for time adapt time
   subroutine startAdaptTime(this)
      class( Cpu_t ) :: this

      if (this%adaptIsRunning) then
         print*, 'Adapt CPU time is already running!'
         stop
      else if (this%estimIsRunning) then
         print*, 'ANOTHER CPU RUNNING: Estim in Adapt!'
      else if (this%prepareIsRunning) then
         print*, 'ANOTHER CPU RUNNING: Prepare in Adapt!'
         stop
      else if (this%solveIsRunning) then
         print*, 'ANOTHER CPU RUNNING: Solve in Adapt!'
         stop
      else
         call cpu_time(this%adaptStart)
         this%adaptIsRunning = .true.
      endif

   end subroutine startAdaptTime

   !> add CPU time spent on adaptation
   subroutine addAdaptTime(this)
      class( Cpu_t ) :: this
      real :: t

      if (this%adaptIsRunning) then
         call cpu_time(t)
         this%adapt = this%adapt + t - this%adaptStart
         this%adaptStart = 0.0
         this%adaptIsRunning = .false.
      else
         print*, 'Cannot add CPU AdaptTime! Adapt is not running!'
         stop
      endif

   end subroutine addAdaptTime

   !> start clock for time estimation time
   subroutine startEstimTime(this)
      class( Cpu_t ) :: this

      if (this%estimIsRunning) then
         print*, 'Estim CPU time is already running!'
      else if (this%adaptIsRunning) then
         print*, 'ANOTHER CPU RUNNING: Adapt in Estim!'
         stop
      else if (this%prepareIsRunning) then
         print*, 'ANOTHER CPU RUNNING: Prepare in Estim!'
         stop
      else if (this%solveIsRunning) then
         print*, 'ANOTHER CPU RUNNING: Solve in Estim!'
         stop
      else
         call cpu_time(this%estimStart)
         this%estimIsRunning = .true.
      endif

   end subroutine startEstimTime

   !> add CPU time spent on estimation
   subroutine addEstimTime(this)
      class( Cpu_t ) :: this
      real :: t

      if (this%estimIsRunning) then
         call cpu_time(t)
         this%estim = this%estim + t - this%estimStart
         this%estimStart = 0.0
         this%estimIsRunning = .false.
      else
         print*, 'Cannot add CPU EstimTime! Estim is not running!'
         stop
      endif

   end subroutine addEstimTime

   !> start clock for time preparation time
   subroutine startPrepareTime(this)
      class( Cpu_t ) :: this
!      real, intent(out) :: start_time
      real :: t

      if (this%prepareIsRunning) then
         print*, 'Prepare CPU time is already running!'
         stop
      else if (this%estimIsRunning) then
         print*, 'ANOTHER CPU RUNNING: Estim in Prepare!'
      else if (this%adaptIsRunning) then
         print*, 'ANOTHER CPU RUNNING: Adapt in Prepare!'
         stop
      else if (this%solveIsRunning) then
         print*, 'ANOTHER CPU RUNNING: Solve in Prepare!'
         stop
      else
         call cpu_time(this%prepareStart)
         this%prepareIsRunning = .true.
      endif

   end subroutine startPrepareTime

   !> add CPU time spent on preparation
   subroutine addPrepareTime(this)
      class( Cpu_t ) :: this
      real :: t

      if (this%prepareIsRunning) then
         call cpu_time(t)
         this%prepare = this%prepare + t - this%prepareStart
         this%prepareStart = 0.0
         this%prepareIsRunning = .false.
      else
         print*, 'Cannot add CPU PrepareTime! Prepare is not running!'
         stop
      endif

   end subroutine addPrepareTime

   !> start clock for time solve time
   subroutine startSolveTime(this)
      class( Cpu_t ) :: this

      if (this%solveIsRunning) then
         print*, 'Solve CPU time is already running!'
         stop
      else if (this%estimIsRunning) then
         print*, 'ANOTHER CPU RUNNING: Estim in Solve!'
      else if (this%adaptIsRunning) then
         print*, 'ANOTHER CPU RUNNING: Adapt in Solve!'
         stop
      else if (this%prepareIsRunning) then
         print*, 'ANOTHER CPU RUNNING: Prepare in Solve!'
         stop
      else
         call cpu_time(this%solveStart)
         this%solveIsRunning = .true.
      endif

   end subroutine startSolveTime

   !> add CPU time spent on Solution of the algebraic systems
   subroutine addSolveTime(this)
      class( Cpu_t ) :: this
      real :: t

      if (this%solveIsRunning) then
         call cpu_time(t)
         this%solve = this%solve + t - this%solveStart
         this%solveStart = 0.0
         this%solveIsRunning = .false.
      else
         print*, 'Cannot add CPU SolveTime! Solve is not running!'
         stop
      endif

   end subroutine addSolveTime


!   subroutine addSolveTime(this, start_time)
!      class( Cpu_t ) :: this
!      real, intent(in) :: start_time
!      real :: t
!
!      call cpu_time(t)
!      this%solve = this%solve + t - start_time
!
!   end subroutine addSolveTime
!
!   subroutine addEstimTime(this, start_time)
!      class( Cpu_t ) :: this
!      real, intent(in) :: start_time
!      real :: t
!
!      call cpu_time(t)
!      this%estim = this%estim + t - start_time
!
!   end subroutine addEstimTime
!
!   subroutine addAdaptTime(this, start_time)
!      class( Cpu_t ) :: this
!      real, intent(in) :: start_time
!      real :: t
!
!      call cpu_time(t)
!      this%adapt = this%adapt + t - start_time
!
!   end subroutine addAdaptTime


   subroutine cleanTimes(this)
      class( Cpu_t ) :: this

      this%prepare = 0.0
      this%adapt   = 0.0
      this%estim   = 0.0
      this%estim2  = 0.0
      this%solve   = 0.0
      this%end_time= 0.0

      this%adaptStart = 0.0
      this%estimStart = 0.0
      this%prepareStart = 0.0
      this%solveStart = 0.0

   end subroutine cleanTimes

   subroutine printCpuTimes( this)
      class( Cpu_t ) :: this
      real :: t, t_total, t_lost

      call cpu_time(t)
      t_total = t-this%start_time
      t_lost = t_total - this%totalTime()

      print*,'# ADGFEM finished! CPU TIMES: '
      print*, 'Prepare | Solve |Estimate | Adapt | TOTAL | Lost'
      write(*,'(a1, f7.2,a1, f7.2, a1, f7.2, a1, f7.2, a1, f7.2, a1, f7.2)') &
         ' ', this%prepare, ' ', this%solve,' ', this%estim, ' ', this%adapt, &
         ' ', t_total, ' ', t_lost

   end subroutine printCpuTimes

   function totalTime( this) result (t)
      class( Cpu_t ) :: this
      real :: t

      t = this%prepare + this%estim + this%adapt + this%solve + this%estim2

   end function totalTime

end module cpu_time_mod
