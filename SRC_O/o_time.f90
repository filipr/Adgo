module time_mod
!   use tdg_mod
 !  use mblock_mod
 use geometry
 use integration

   implicit none

   private

   type, public :: Time_t
      character(len=20) :: disc_time ! duplicated from state, should replace it in future
      character :: time_method             ! 'E'= explicit, 'I' = implicit,'S'=semi-implicit! duplicated from state, should replace it in future
      !FERROR should be loaded in readMainData()
      character(len=20) :: estim_time
      !character(len=20) :: adapt_time

      integer :: deg

      real :: tau_old                          ! time step from previous comput
      real, dimension(:), allocatable :: tau       ! time steps
      real :: tau_new                          ! new time step
      real, dimension(:), allocatable :: rat_tau   ! ratio tau, array for time step controll
      integer :: unchanged_time_steps       ! number of unchanged number of time steps
      integer :: max_Tdof                      ! maximal of elem(i)%Tdof

      logical :: stdgm
      logical:: tdg, cn                        ! DG method for time ?; Crank-Nicolson?

      real :: ttime                        ! totaL time
      real :: ctime                        ! current time
      real :: FinTime                      ! Final  time
      real :: OutTime                      ! output time

      real :: ttime_save, timeprn_save     ! total time save for possible adaptation !stays in state

      integer :: iter                      ! index of actual time step
      integer :: iter_loc                  ! local index of time step (within one adapt. loop)
      integer :: maxiter                   ! maximum number of time steps
      integer :: iter_SC                   ! index of a time step when Stopping Criteria based on AEE are satisfied

      logical :: keep_tau
      integer :: recompute_back

      logical :: tau_fixed                 ! time step is fixed in *.ini
      real :: tau_fixed_size               ! size of ^^^^^^^^^^^^^^^^
      character(len=20) :: tau_choice      ! time step choice type
      real :: CFL                          ! maximal allowed CFL number

   ! FERROR : BDFtol is used also for stdgm? ??? YES
     real :: BDFtol                       ! tolerance for ABDF

     !for STDG only
      character(len=10) :: quadType
      class(Time_rule), dimension(:), allocatable :: T_rule     ! integ rules for time intervals
!!!!! BDF !!!!!!!!!!!!
      integer :: deg_actual               ! actual order of BDF time discretization = min(iter, Tdeg)
      integer :: Qnum ! the degree (number of nodes) of time quadrature used
      !( it has to be synchornized for RTNst estims and also the RES methods needs it to be set right because of the deg_plus !!!


      !FOR ABDF type in BDF type only
      integer :: max_deg      ! maximal implemented degree
      !integer :: deg          ! degree of this scheme
      !integer :: deg2         ! degree of the second scheme for the error estim.
      real, allocatable, dimension(:) :: alpha
      real, allocatable, dimension(:) :: Bextrap
      real, allocatable, dimension(:) :: extrap
      real, allocatable, dimension(:) :: delta
      real  :: gamm

   contains
!     procedure, deferred ::  performOneTimeStep

      procedure :: init => initTime
      procedure :: initTimeStepAdapt
      procedure :: AllocateBDF
      procedure :: FakeBDFLagrangeExtr
      procedure :: SetBDF

      procedure :: createTimeQuadRulesAndBasis
      procedure :: write_TrulePhi      ! should be moved to T-quadrature

   end type Time_t

   type, public, extends(Time_t) :: TimeTDG_t
!      character(len=10) :: quadType
!      type(Time_rule), dimension(:), allocatable :: T_rule     ! integ rules for time intervals
      type(Mblock), allocatable  :: refTimeMatrix   !the maximal (Tdof = max_Tdof+1) time part of the Mass Matrix in STDGM method, +1 - Tdof_plus
      type(Mblock), allocatable  :: StiffTimeMatrix   !the maximal (Tdof = max_Tdof+1) time part of the Mass Matrix in STDGM method, +1 - Tdof_plus

      contains

!      procedure :: createTimeQuadRulesAndBasis
!      procedure :: write_TrulePhi


   end type TimeTDG_t

   type, public, extends(Time_t) :: TimeBDF_t

     !integer :: deg_actual


     !moved from ABDF
!     integer :: max_deg      ! maximal implemented degree
!!     integer :: deg          ! degree of this scheme
!     !integer :: deg2         ! degree of the second scheme for the error estim.
!!     real, pointer, dimension(:) :: alpha
!     real, pointer, dimension(:) :: Bextrap
!     real, pointer, dimension(:) :: extrap
!     real, pointer, dimension(:) :: delta
!     real  :: gamm


   contains

   end type TimeBDF_t

!
!   abstract interface
!      subroutine performOneTimeStep(this)
!
!         implicit none
!         import :: DiscTime_t
!         class(DiscTime_t) :: this
!      end subroutine
!   end interface

   contains
   !> initialization of state%time, should be used after the ini file is read
   !> quadratures, basis, method logicals
   subroutine initTime( this )
      class( Time_t ), intent( inout ) :: this
      integer :: i


      this%cn = .false.  ! OLD techniques
      this%tdg = .false. ! OLD techniques
      this%stdgm  = .false.


      if ( this%deg > MaxTimeDegree) then
       print*, 'Desired time degree is greater than implemented. Stopping'
       stop
      endif

      select type(this)
      type is (TimeTDG_t)
         this%stdgm = .true.
         this%quadType = 'Radau' !'Gauss'
         !write(*,*) '# Gauss quadrature is used for the time integration. (Change to Radau in initTime)'
         write(*,*) '# Radau quadrature is used for the time integration. (Change to Gauss in initTime)'

         call this%createTimeQuadRulesAndBasis()

      type is (TimeBDF_t)
        ! stop 'initTime - BDF not yet done in ADGO'
        this%quadType = 'Gauss'
         write(*,*) '# Gauss quadrature is used for the time integration. (Change to Gauss-Radau in initTime)'
         call this%createTimeQuadRulesAndBasis()

      class default
         stop 'Not defined type of time disc.'
      end select

      this%keep_tau = .false.




      !set time basis functions

!    print*,'begin in SDERFTRE', maxTrule
!     integ rules for time intervals, connected with the GAUSS at this moment
!    do iq=1, maxTrule
!       !print*,'*********************',iq,maxTrule
!       call SetTrulePhi(state%time%T_rule(iq))
!       !do i=1,state%time%T_rule(iq)%Qdof
!       !   write(100+iq, *) state%time%T_rule(iq)%lambda(i),state%time%T_rule(iq)%phi(:,i)
!       !   write(200+iq, *) state%time%T_rule(iq)%lambda(i),state%time%T_rule(iq)%Dphi(:,i)
!       !enddo
!
!       !if(iq == 5) stop
!    enddo
!    !print*,'stopped in SDERFTRE'
!    !stop
!
!    !writing all state%time%T_rule(:)%Phi to file TRule_Phi
!    !call Write_TrulePhi()


   end subroutine initTime

 !> adaptation of the time step, should be called in the beginning when the ini file is read
 subroutine initTimeStepAdapt( this, val )
   class( Time_t ), intent( inout ) :: this
   real, intent( in ) :: val


   if(this%tau_choice == 'fixed') then
      this%tau_fixed = .true.
      this%tau_new = val
      this%tau_fixed_size = val
      !if (this%tau_new <= 0.D+00) this%tau_fixed = .false.
      !this%tau_fixed_size = this%tau_new

    write(*,'(a25,es12.4)') &
          '  # Fixed time step tau = ',this%tau_new

   elseif(this%tau_choice == 'cfl') then
      this%CFL = val
    write(*,'(a20,es12.4)') &
          '  # Fixed CFL number = ',this%CFL

   elseif(this%tau_choice == 'exp') then
      this%CFL = val

    write(*,'(a40,es12.4)') &
          '  # Exponentially increasing  CFL number = ',this%CFL

   elseif(this%tau_choice == 'adapt') then
      this%BDFtol = val

      if(this%estim_time == 'loc') then
         write(*,'(a57,es12.4)') &
              ' # Adaptively chosen time step based on LOC estim: tol =',this%BDFtol
      elseif(this%estim_time == 'tRES') then
         write(*,'(a57,es12.4)') &
              ' # Adaptively chosen time step based on RES estim: tol =',this%BDFtol
      else
         print*,' Unknown type of the adaptive choice of the time step:  ',this%estim_time
         print*, ' Only techniques "loc" and "tRES" implemented !'
         stop
      endif
  else
     print*,' Unknown type of the choice of the time step:',this%tau_choice
     print*, 'Possibilities are: fixed, cfl, exp, adapt'
  endif


 end subroutine initTimeStepAdapt


  !> allocate BDF method
  subroutine AllocateBDF( this, max_deg)
    class( Time_t ), intent(inout) :: this
    integer, intent(in) :: max_deg

    this%max_deg = max_deg
    allocate(this%alpha(0:max_deg+1) )
    allocate(this%delta(0:max_deg+1) )   ! for estimation of the local error

    allocate(this%extrap(1:max_deg) )
    !allocate(BDF%Bextrap(1:max_deg+1) )
  end subroutine AllocateBDF

    !> set up a fake BDF formula to simply perform Langrange extrapolation
  subroutine FakeBDFLagrangeExtr( this, deg, tau)
    class (Time_t), intent(inout) :: this
    integer, intent(in) :: deg
    real, dimension (1:deg+1), intent(in) :: tau
    integer:: n,i

    this%deg = deg + 1
    n = size(tau)
    ! FIXME: This is the classical formula. Go via DD eventually!
    do i = 1,n
       ! the leading minus is due to minus used in BDF extrapolation
       this%alpha(i) = -product(-tau/(tau(i) - tau),mask=tau/=tau(i))
    end do
  end subroutine FakeBDFLagrangeExtr


  !> set coefficients \f$ \alpha_k,\ \beta_k,\ \dots \f$ of ABDF method
  subroutine SetBDF( this, deg, tau )
    class( Time_t ), intent(inout) :: this
    integer, intent(in) :: deg
    real, dimension (1:deg+1), intent(in) :: tau

    real :: theta0, theta1, theta2, A0, A1, A2, B0, B1, C0, C_BDF, C_EXT
    real, dimension(:,:), allocatable :: lag
    real :: factorial, prod
    integer :: i,j,k, version

    if(deg > this%max_deg) then
       print*,'Error in SetBDF, try to set degree ',deg,&
            ', but there is allocated the degree ',this%max_deg
       stop
    endif

   ! NEW version
   this%deg = deg

   ! Lagrangian interpolation \sum_{i=0}^{deg(+1)} w(i,:) L_i(t)

   allocate(lag(0:deg+1, 0:deg+1) )

   factorial = 1.
   do i=0, deg+1
       lag(i,i) = 0.
       do j= i+1, deg + 1
          lag(i, j) = sum(tau(i+1:j) ) ! lag(i,j) = t_{k-i} - t_{k-j}
          lag(j, i) = - lag(i, j)
       enddo
       if(i> 0) factorial = factorial * i   ! factorial = (deg+1)!
   enddo


!   if ( .not. allocated( this%tau ) ) then
!      stop 'SetBDF only for bdf method, time%tau not allocated'
!   endif
!
!   print*, 'tau:', size( this%tau )
!   print*, this%tau(1)
!
!
!   print*, 'size of alpha', size( this%alpha(:) ), 'and deg: ', deg


   do i= 0, deg
    ! alpha(i) = (L_i)' (t_k)
    this%alpha(i) = 0.

    do j= 0, deg
       if(j /= i) then
          prod = 1./lag(i,j)

          do k=0, deg
             if(k /= j .and. k/=i )  prod = prod * lag(0,k)/(lag(i,k))
          enddo

          this%alpha(i) = this%alpha(i) + prod
       endif
    enddo

    this%alpha(i) = this%alpha(i) * tau(1)
   enddo


   ! coefficients delta for local error estimate
   ! delta(i) = (deg+1)-th derivative of L_i, = const

   lag(0:deg+1, 0:deg+1) = lag(0:deg+1, 0:deg+1) / tau(1)

   do i=0, deg+1
    lag(i,i) = 1.
   enddo

   do i=0, deg + 1
    this%delta(i) =  factorial / product(lag(i, 0:deg+1) )
   enddo

   ! coefficient gamma
   this%gamm = -1.
   do i=1,deg
    this%gamm = this%gamm * sum(tau(1:i)) /(tau(1) * (i+1))
   enddo


   !write(*,'(a7,i5,es12.4,a1,8es12.4)') 'alpha',deg,sum(this%alpha(0:deg)), '|', &
   !     this%alpha(0:deg)
   !write(*,'(a7,i5,es12.4,a1,8es12.4)') 'delta',deg,sum(this%delta(0:deg+1)),'|', &
   !     this%delta(0:deg+1)
   !write(*,'(a7,i5,2es12.4)') 'gamma',deg, this%gamm, 1./(deg+1)
   !print*,'----------------------'


   this%deg = deg
   select case (deg)
   case(:1)   ! EBDF extrapolated BDF
    this%extrap(1) = 1.

   case(2)
    theta0 = tau(1)/tau(2)
    theta1 = tau(2)/tau(3)

    this%extrap(1) = theta0 + 1.
    this%extrap(2) = -theta0

   case(3)
    theta0 = tau(1)/tau(2)
    theta1 = tau(2)/tau(3)
    theta2 = tau(3)/tau(4)

    A0 = theta0 + 1
    A1 = theta1*A0 + 1
    A2 = theta2*A1 + 1

    B0 = theta1 + 1
    B1 = theta2*B0 + 1

    this%extrap(1) = A0*A1/B0
    this%extrap(2) = -theta0*A1
    this%extrap(3) = theta0*theta1**2 *A0/B0

   case(4:)

    if( this%iter <= 5) &
         print*,'m-step extrapolated BDF scheme for m>3 is not implemented'
    stop

   end select

  end subroutine SetBDF

!!!!!!!!!!!    TDG                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !> set Time integ Quad rules and time basis functions Phi and Dphi
   subroutine createTimeQuadRulesAndBasis( this )
      class( Time_t ), intent( inout ) :: this
      integer :: i

      select case ( this%quadType )
      case('Radau')
         !global variable!!!!
         maxTrule = maxRrule
         !print*, 'Radau in createTimeQuadRule'
         allocate( RadauTime_rule :: this%T_rule(1:maxTrule) )
         do i=1,maxTrule
            call this%T_rule(i)%createTimeRule(i)
            call this%T_rule(i)%SetTruleBasis()
         end do

      case('Gauss')
         maxTrule = maxGrule
         allocate( GaussTime_rule :: this%T_rule(1:maxTrule) )
         do i=1,maxTrule
            call this%T_rule(i)%createTimeRule(i)
            call this%T_rule(i)%SetTruleBasis()
         end do
      end select

      !    writing T_rule(i)%Phi to file TRule_Phi
      !    call this%T_rule(3)%WriteTruleBasis()
      !    stop


   end subroutine createTimeQuadRulesAndBasis




  subroutine write_TrulePhi(this)
   class( Time_t ), intent(in) :: this
   integer :: i,j , Qdeg


   open (58, file="T_rulePhi", action ="write", status="replace")

   do j = 1, maxTrule
    !  T_rule => state%time%T_rule(j)
    !  Qdeg = T_rule%Qdeg
    !  Qdof = T_rule%Qdof

   write(58,*) , 'Trule:' , this%T_rule(j)%Qdeg
   do i = 1, MaxTimeDegree + 2
      write(58, * ) , this%T_rule(j)%phi(i,:)
   enddo
   write(58,*) '---------------'

   enddo !j
   close(58)

  end subroutine write_TrulePhi




end module time_mod
