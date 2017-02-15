module dwr_mod
   use dwr_alg_mod
   use elemental_mod
   use target_functional_mod

   implicit none

   type :: DWR_t
      integer :: id !type of the DWR target functional (will be specified in paramets)
      !  type(Newton_type), pointer :: Newton !Newton arrays for dual problem
      logical :: RESindicators ! if .true. : the adaptation of the mesh uses the  RES error indicators instead of DWR
      logical :: deg_plus ! dual problem solution for p + p_plus
      integer :: R_type ! type of the reconstruction , 0 - L2, 1 - H1, -1 - H-1, 2 - reconstructiton based on a(u,vp+) = res(vp+)
      integer :: eta_index ! which eta should be used to adaptation of the mesh (temporarily set in the program), should be in ini
      class( Target_functional_t ), allocatable :: J ! target functional
      logical :: dualProblem_computed
      real :: linSolver_tol ! tolerance of the linsolver for dual problem
      integer :: linSolver_iter ! number of lin Solver iterations
      real :: residuum
      real, dimension(:), allocatable ::  x, b, b1, rr  ! arrays for computation
      type( Elemental3_t ), dimension(:), allocatable :: rhs
      type( Elemental3_t ), dimension(:), allocatable :: dualRes ! the residual vector of the dual problem
      type( DWR_alg_t ) :: aDWR ! class of all needed variables for the aDWR method for error estimation of the alg. errors

      logical :: PU ! use the partition of unity for localization of the estimates
      type( Elemental1_t ), dimension(:), allocatable :: pu_estim ! used for local estims using PU
      real, dimension(:), allocatable :: vertex_estim ! used for local vertex estims using PU, 1:npoin
      integer, dimension(:), allocatable :: num_of_vertex_elem ! used for local vertex estims using PU, 1:npoin

      real :: estimLP, estimLD, estimNL ! linear primal and dual estimate, nonlinear estimate
      real :: estimDiscr ! discretization estimate
      real :: J_uhn ! actual value of the target functional

   contains
      procedure :: delete => del_DWR
      procedure :: init => initDWR
      procedure :: update => updateDWR
      procedure :: clean => cleanDWR
      procedure :: setRHS_new

      procedure :: fillRHSvector ! %rhs -> %b
      procedure :: fillDualResVector ! dualResidual
      procedure :: distributePuEstimToVertices
      procedure :: initDualLinSolver
      procedure :: etaFromVertices
      procedure :: writeNlErrorFile
!      procedure :: newtonDone

   end type DWR_t

!   interface DWR_t
!       procedure :: initDWR
!      !procedure :: otherWayOfInitDWR - e.g. copy constructor
!   end interface DWR_t

!   public :: PrepareDualProblem
    public :: InitDualProblem
    public :: computeDualLinearDWRestimates
    public :: computePrimalLinearDWRestimates

   !> main structure for DWR error estimates
   type( DWR_t ), allocatable :: DWR

   contains
    !> init basic DWR structure
    !> should be called only once - in the beginning of computation
    !  function initDWR( id, grid) result(this) - it does new allocation???
    subroutine initDWR( this, grid )
      class(DWR_t), intent(inout) :: this
      class( mesh ), intent(in) :: grid
      integer :: i
      integer :: iFile = 41
      real, dimension(1:nbDim) :: xy_coord

      if ( state%time%deg > 1 ) &
            stop 'in initDWR - stationary problems work only for q=0'

      this%dualProblem_computed = .false.

      !set DWR_PARAMETERS
      this%PU = .false. !.true. ! use partition of unity to localization
      if (this%PU) &
         print*, ' # DWR: Vertex oriented partition of unity is used for mesh adaptation.'

      !set DWR_PARAMETERS
      ! type of reconstruction
      this%R_type = 1 ! 1 - H1, 0 - L2, 2 - Ritz
      if (this%R_type <= 1) then
         print*, '# DWR: We use the LS reconstruction.'
      else if (this%R_type==2) then
         print*, '# DWR: We use the reconstruction based on a(up,vp+) = res(vp+) forall vp+ \in Shpp!'
      else
         stop 'Unknown type of reconstrution.'
      endif

      !set DWR_PARAMETERS
      ! tolerance of the lin solver for dual problem
      this%linSolver_tol = state%linSolver%tol
      print*, '# DWR: Linear solver tolerance for DUAL problem set to ' , this%linSolver_tol

      !set eta index
      !set DWR_PARAMETERS
!      this%eta_index = dwrS
!      write(*,*) '# DWR: PRIMAL residual estimator used for mesh adaptation!'
!      this%eta_index = dwr_dualS
!      write(*,*) '# DWR: DUAL residual estimator used for mesh adaptation!'
       this%eta_index = dwr_aver
      write(*,*) '# DWR: AVERAGE residual estimator used for mesh adaptation!'


      ! init the target functional this%J
      select case(this%id)
         case(1)
            write(*,'(a50)') '  # Target functional: L2-norm of the solution - NOT IMPLEMENTED'
!            this%linear = .false.
            stop
         case(2)
            write(*,'(a50)') '  # Target functional: H1-seminorm of the solution - NOT IMPLEMENTED'
!            this%linear = .false.
            stop
         case(3)
            write(*,'(a50)') '  # Target functional: Point value of the solution in (0.5, 0.5)'
            ! FERROR coordinates should be set in the .ini file in future
            xy_coord(1:2) = (0.75, 0.25)
            allocate( Point_value_t :: this%J )
            call this%J%init( xy_coord, grid )
         case(4)
            write(*,'(a100)') '  # Target functional: Integral of the solution u(T) over subdomain specified in the submesh file.'
            write(*,'(a100)') '  # For now it approximates the value in the middle (divides by the area of the support)!'
            ! FERROR coordinates should be set in the .ini file in future
            xy_coord(1:2) = (0.0, 0.0) ! not used here
            allocate( U_over_subdomain_t :: this%J )
            call this%J%init( xy_coord, grid )
         case(5)
            write(*,'(a100)') '  # Target functional: Integral of the solution du/dx over subdomain specified in the submesh file.'
            ! FERROR coordinates should be set in the .ini file in future
            xy_coord(1:2) = (0.0, 0.0) ! not used here
            allocate( dudx_t :: this%J )
            call this%J%init( xy_coord, grid )

         case(6:)
            write(*,'(a50)') '  # Unknown type of target functional! STOPPING!'
            stop
      end select

      this%estimDiscr = 1.0
      ! temporal for export
      this%aDWR%file_error_new = "aDWR_nlErrors"

      if ( state%nlSolver%non_alg_stop == 'aDWR' ) &
         call this%aDWR%init()
      if ( state%nlSolver%non_alg_stop == 'aDWR' .and. DWR%deg_plus) &
         stop 'the DWR_P method is nnot compatible with aDWR algebraic estimates!'


   end subroutine initDWR

   ! final subroutine, cleaning the DWR structure
   ! should be called in the end of the computation
   subroutine del_DWR(this)
     class( DWR_t ) :: this
     integer :: i,l

     call this%clean()

     if (allocated(this%J)) then
         call this%J%delete()
         deallocate(this%J)
     endif

   end subroutine del_DWR

   !> update DWR structures in the beginning of each solution procedure
   !> should be called after each mesh adaptation
   subroutine updateDWR( this, grid )
      class(DWR_t), intent(inout) :: this
      class( mesh ), intent(in) :: grid
      integer :: i

      if ( state%nlSolver%non_alg_stop == 'aDWR' ) then
       !print*, 'cleanDWR called! '
         call this%clean()
         call this%aDWR%update()
      endif

      !this%estimDiscr = 1.0
      this%estimNL = 1.0
      this%estimLD = 1.0
      this%estimLP = 1.0

      !moved from InitDualProblem( this, this%deg_plus )
      if (state%time%quadType /= 'Radau') &
         stop 'Radau quadrature must be used - to compute C(u(T))'
      if (.not. allocated(this%J)) &
         stop  'DWR%J is not allocated !!!'
      if ( this%J%boundary ) &
         stop 'We have to change the matrix C(u) for boundary target functionals!'

      call this%J%findSupp( grid )

      if ( this%J%isupp == 0 ) then
         print*, 'epsilon: ', DWR%J%eps, 'grid%h/2', grid%h / 2.0
         stop 'zero support of target func'
      endif

      call this%setRHS_new( grid )

      if (this%deg_plus) &
         print*, 'Control updateDWR for DWR_P - length of ZST will be decreased in DWRestims'

      if (state%time_dependent) &
         stop 'we need to compute IC for dual problem in updateDWR!'

      if ( state%nlSolver%non_alg_stop == 'aDWR' ) then
       !print*, 'InitDualLinSolver in PerformOneSTDGMstep'
         call this%initDualLinSolver( grid )
      endif

   end subroutine updateDWR

   ! subroutine cleaning the structures in DWR,
   ! which should be called after each adaptation
   subroutine cleanDWR(this)
      class( DWR_t ) :: this
      integer :: i , l

      if (allocated(this%b)) &
         deallocate( this%b )
      if (allocated(this%x)) &
         deallocate( this%x )
      if (allocated(this%b1)) &
         deallocate( this%b1 )
      if (allocated(this%rr )) &
         deallocate( this%rr )

      if (allocated( this%rhs ) ) then
         l = size( this%rhs(:) )
         do i = 1,l
            call this%rhs(i)%delete()
         enddo
         deallocate( this%rhs )
      endif

      if (allocated( this%dualRes ) ) then
         l = size( this%dualRes(:) )
         do i = 1,l
            call this%dualRes(i)%delete()
         enddo
         deallocate( this%dualRes )
      endif

      this%linSolver_iter = 0
      this%residuum = 0
      this%dualProblem_computed = .false.

      if (allocated( this%pu_estim ) ) then
         l = size( this%pu_estim(:) ) ! = nelem
         do i = 1,l
            call this%pu_estim(i)%delete()
         enddo
         deallocate( this%pu_estim )
      endif
      if (allocated(this%vertex_estim)) &
         deallocate(this%vertex_estim)
      if (allocated(this%num_of_vertex_elem)) &
         deallocate(this%num_of_vertex_elem)

      call this%J%clean()

   end subroutine cleanDWR


   !> new version: save the dual_rhs to the vector b(1:nsize) DOES NOT use elem%rhsST
   !> set RHS for the dual problem used for DWR error estimation method
   subroutine setRHS_new(this, grid)
    class( DWR_t ), intent(inout) :: this
    class( mesh ), intent(in) :: grid
    class( element ), pointer :: elem
    integer :: i,j, dof, Qdof, k,l, nelem, i_elem, Tdof
    real :: time
    real, dimension(:,:), allocatable :: f


    if ( allocated( this%rhs ) ) then
!      print*,
!      print*, 'Rhs is already allocated -> do nothing in setRHS_new'
!      print*,
    else
       if ( this%J%isupp == 0 ) then
          print*, 'J epsilon:', this%J%eps
          stop 'zero support of target func'
       endif
       associate( JJ => this%J )
   !    select type ( JJ )
   !         type is ( dudx_t )
   !            if ( JJ%dx /= 1 .and. JJ%dx /= 2 ) &
   !               stop 'DWR%J%dx is not set correctly in setRHS_new'
   !    end select

       nelem = grid%nelem
       allocate( this%rhs(1:grid%nelem) )

       do i = 1, nelem
         elem => grid%elem(i)
         dof = elem%dof_plus
         Tdof = elem%Tdof
         this%rhs(i) = Elemental3_t( ndim, dof,Tdof )
       enddo

       ! TIME dependent version
       if (state%time_dependent) then
         stop 'set RHS not implemented for nonstationary problems'

         print*, 'F@R: setRHS - set time (temporarily finTime)'
         time = state%time%finTime

       ! TIME independent
       else
         time = state%time%finTime

         do i = 1, this%J%isupp
           i_elem = this%J%supp(i)
           elem => grid%elem( i_elem )
           Qdof = elem%Qdof

           if ( elem%deg_plus ) then
             dof = elem%dof_plus
           else
             dof = elem%dof
           end if


           allocate( f(1:Qdof, 1:ndim) , source = 0.0 )
           do j = 1,Qdof
             f(j,1:ndim) = this%J%dualRHS_scalar(elem%xi(0,j,1:nbDim), time)
           end do !j
           ! TODO CONTROL SIZE OF VEC
           ! computes \int_{K_{elem}} f \phi_{R} \ dx \f$, for phi_1:dof

            ! set rhs
   !        call EvalVectorB( elem, f(1:Qdof,1:ndim), dof, elem%vec(rhs,1:ndim*dof) )
           select type ( JJ )
            class is ( dudx_t )
               call EvalVectorDphiDx_2dim( elem, f(1:Qdof, 1:ndim), &
                  dof, this%rhs(i_elem)%x(1:ndim, 1:dof, 1), JJ%dx )
            class default
               call EvalVectorB_2dim( elem, f(1:Qdof,1:ndim), &
                  dof, this%rhs(i_elem)%x(1:ndim, 1:dof, 1) )
           end select

           deallocate(f)

           if (elem%Tdof > 1) &
               stop 'problem in setRHS_new'

         enddo ! i
       endif ! time dependent

    end associate ! JJ

    endif ! allocated rhs

   end  subroutine setRHS_new

   ! fill the global dual res vector to DWR%dualRes(:,:,:)
   subroutine fillDualResVector( this , grid, nsize, b )
      class( DWR_t ), intent(inout) :: this
      class( mesh ), intent(in) :: grid
      integer, intent(in) :: nsize ! size of the global vector
      real, dimension(1:nsize), intent(in) :: b ! the vector containing A^T*zST
      class( element ), pointer :: elem
      integer :: i, k, l
      integer :: ivec, kvec, dof_plus, ndof, nelem, Tdof

      nelem = grid%nelem
      if ( allocated(this%dualRes) ) then
          deallocate(this%dualRes)
      endif
      allocate( this%dualRes(1:nelem) )

      ivec = 0
      do i=1,nelem
         elem => grid%elem(i)
         dof_plus = elem%dof_plus
         Tdof = elem%Tdof
         ndof = dof_plus * Tdof * ndim
         !allocate the array
         this%dualRes(i) = Elemental3_t( ndim, dof_plus,Tdof )

         kvec = ivec

         ! save to dualRes the dual residual vector
         ! Ritz - need this vector to rhs for the local problem used for computation of zSTplus
         do l = 1, Tdof
            do k = 1, ndim
               this%dualRes(i)%x(k,1:dof_plus,l) = &
                  this%rhs(i)%x(k,1:dof_plus,l) - b(kvec+1:kvec + dof_plus)
               kvec = kvec + dof_plus
            enddo !k
         enddo !l
         ivec = ivec + ndof
      end do !i

   end subroutine fillDualResVector

   !> fills the 3-dim RHS array this%rhs(:,:,:) to vector of the RHS this%b
   subroutine fillRHSvector( this, grid, plus )
      class( DWR_t ) :: this
      class( mesh ), intent(in) :: grid
      logical, intent(in) :: plus
      class( element ), pointer :: elem
      integer :: nsize, i, nelem, k, l, dof, ndof, kvec, ivec

      nelem = size( this%rhs )
      ivec = 0

      if ( .not. allocated(this%rhs) ) &
         print*, 'DWR rhs is not allocated!'

      if (allocated(this%b)) &
         deallocate( this%b )

      if (plus) then
        nsize = sum( ndim*grid%elem(:)%dof_plus * grid%elem(:)%Tdof )
        !print*, 'nsize', nsize
        allocate( this%b(1:nsize), source = 0.0 )

        do i = 1,nelem
            elem => grid%elem(i)
            dof = elem%dof_plus
            ndof = ndim * dof * elem%Tdof
            kvec = ivec

            do l = 1, elem%Tdof
               do k = 1, ndim
                  this%b(kvec+1:kvec + dof) = this%rhs(i)%x(k,1:dof,l )
                  kvec = kvec + dof
               end do ! k
            enddo !l
            ivec = ivec + ndof

         end do

      else
         nsize = sum( ndim*grid%elem(:)%dof * grid%elem(:)%Tdof )
         allocate( this%b(1:nsize), source = 0.0 )
         do i = 1,nelem
            elem => grid%elem(i)
            dof = elem%dof
            ndof = ndim * dof * elem%Tdof
            kvec = ivec
            do l = 1, elem%Tdof
               do k = 1, ndim
                  this%b(kvec+1:kvec + dof) = this%rhs(i)%x(k,1:dof,l )
                  kvec = kvec + dof
               end do ! k
            enddo !l
            ivec = ivec + ndof
         end do

      endif

   end subroutine fillRHSvector


     !> init C(u(T))^T z = J(\phi) for nonlinear stationary problem
  !> called only with aDWR
  subroutine initDualLinSolver( this, grid)
    class( DWR_t ), intent(inout) :: this
    class( mesh ), intent(in) :: grid
    class( element ), pointer :: elem
    integer :: i, j, kk, k, elemDof, dof, degP, dofP, R_type, Qnum, nsize
    real :: t1, t2, time_prepare, time_solve, time_estim, res_max_val
    character(len=50) :: plot_file =  'plot_sol.gnu' !'../DWR/plot_sol.gnu'
    integer :: iPlot = 39
    real, allocatable, dimension(:,:,:,:) :: Dw ! for the least square reconstruction
    real :: l_norm
    real :: residuum
    integer :: lin_solver_not_conv

    ! TODO:
    ! CPU_TIME

    if ( state%time_dependent ) &
      stop 'DWR method is not implemented for time-dependent problems YET!'

    time_prepare = 0.
    time_solve  = 0.
    time_estim = -1.

    if ( this%deg_plus ) then
      nsize = sum( grid%elem(:)%dof_plus )
    else
      ! change in future to dual_size add control for nonlinear problems size can be ok
      nsize = state%nsize
    endif

    write( debug,*) 'Allocation of DWR%x should be done better!'

    if ( .not. allocated(this%b) ) &
      allocate(this%b(1:nsize), source = 0.0)
    if ( .not. allocated(this%b1) ) &
      allocate(this%b1(1:nsize) , source = 0.0)
    if ( .not. allocated(this%x) ) &
      allocate(this%x(1:nsize), source = 0.0 )
    if ( .not. allocated(this%rr) ) &
      allocate(this%rr(1:nsize), source = 0.0)

    call cpu_time(t1)

    eta = 0.0

    ! not needed - fills primal RHS
    !    call FillVectorST( this%b(1:state%nsize), eta ) ! value of eta is not used, only its sign

    !filling dwr%rhs into global rhs vector
    ! false - already p+1
    call this%fillRHSvector( grid, this%deg_plus )
!    DWR%x(1:state%nsize) = 0.

    call cpu_time(t2)
    time_prepare = time_prepare + t2 - t1

    !print*, 'Nullify the number of lin. iterations before dual problem.'
    this%linSolver_iter = 0
    !maybe globaly ?
    this%residuum = -1 !state%linSolver%residuum
    lin_solver_not_conv = 0

  end subroutine initDualLinSolver


   ! set the global vertex-wise array of the local estimates
   ! sum of pu_estims from the corresponding elements
   subroutine distributePuEstimToVertices( this, grid )
      class(DWR_t), intent(inout) :: this
      class( mesh ), intent(in) :: grid
      class( element ), pointer :: elem
      integer :: i, nelem, npoin,j , glIndex

      nelem = grid%nelem
      npoin = grid%npoin

      if (allocated( this%vertex_estim ) ) &
         deallocate( this%vertex_estim )
      allocate( this%vertex_estim(1:npoin), source = 0.0 )

      ! TRY to use in PU disrtibution to elements
      if (allocated( this%num_of_vertex_elem ) ) &
         deallocate( this%num_of_vertex_elem )
      allocate( this%num_of_vertex_elem(1:npoin), source = 0 )

      do i = 1, nelem
         elem => grid%elem(i)
         !go through vertices
         do j = 1,3
            glIndex = elem%face(idx, j)
            this%vertex_estim(glIndex) = this%vertex_estim(glIndex) + &
               this%pu_estim(i)%x(j)
         enddo

         ! number of elements per vertex
         do j = 1, elem%flen
            glIndex = elem%face(idx, j)
            this%num_of_vertex_elem(glIndex) = this%num_of_vertex_elem(glIndex) + 1
         enddo
      enddo

   end subroutine distributePuEstimToVertices

   function etaFromVertices(this, elem) result( estim )
      class(DWR_t), intent(inout) :: this
      type(element), intent(in) :: elem
      class( element ), pointer :: neighbor
      real, dimension(1:ndim)  :: estim
      integer :: i,j, npoints, glIndex

      if (ndim > 1) &
         stop 'EtaFromVertices works only for scalar case'

      npoints = elem%flen
      estim(:) = 0.0


      do i = 1, npoints
         ! get global index of the elem vertex
         glIndex = elem%face(idx, i)
         estim(1) = estim(1) &
            + ( this%vertex_estim(glIndex) / this%num_of_vertex_elem(glIndex) )
      enddo

   end function etaFromVertices

        ! write one line intp the file for errors from state%estim()
   ! now it is cleaned after every adaptation
   subroutine writeNlErrorFile( this, nelem)
      class( DWR_t ) :: this
      integer, intent(in) :: nelem
!      integer, intent(in) :: iter ! number of the line
!      integer, intent(in) :: iter_lin_primal !
!      integer, intent(in) :: iter_lin_dual ! number of the lin. alg. iterations
      integer :: iFile = 42

      ! SOME NUMBERS NEED to be squarerooted !!!!!!!!!

!      call this%J%computeJu( grid )
      open( iFile, file = this%aDWR%file_error_new, action="write", position="append" )
        write(iFile,'(i6, i6, i6, i6, i6, 8es18.10)') &
         state%space%adapt%adapt_level + 1, nelem, this%aDWR%iter, &
         this%aDWR%iter_lin_primal, this%aDWR%iter_lin_dual, &
         this%J%Ju_exact, this%J%Ju, &
         this%estimDiscr, this%estimNL, &
         this%estimLP, this%estimLD
      close(iFile)

  end subroutine writeNlErrorFile


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! END OF DWR PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !> init the structures for computation of the dual problem
  !> it should be called before every solution process, i.e. after every mesh adaptation
  !> called only once per mesh level
  subroutine InitDualProblem( DWR, plus )
    class( DWR_t ), intent(inout) :: DWR
    logical, intent(in) :: plus ! dual problem computed for deg_plus
    class(element), pointer :: elem
    integer :: i,j
    integer, pointer :: exp1, exp2
    logical :: time_dependent
    character(len=20) :: loc_adapt_space
    character(len=50) :: plot_file_primal = 'plot_sol_primal.gnu' !'../DWR/plot_sol_primal.gnu'
    integer :: iPlotPrimal = 41
    logical :: impl_loc
    integer :: deg, dof, Tdof

    if (state%time%quadType /= 'Radau') &
      stop 'Radau quadrature must be used - to compute C(u(T))'
    ! When preparing the matrix for

    if (.not. allocated(DWR%J)) &
      stop  'DWR%J is not allocated !!!'

    if ( DWR%J%boundary ) &
      stop 'We have to change the matrix C(u) for boundary target functionals!'

    call DWR%J%findSupp( grid )

    if ( DWR%J%isupp == 0 ) then
        print*, 'epsilon: ', DWR%J%eps, 'grid%h/2', grid%h / 2.0
       stop 'zero support of target func'
    endif

    call DWR%setRHS_new( grid )

    if (plus) &
      print*, 'Control InitDualProblem for DWR_P - length of ZST will be decreased in DWRestims'

    ! allocate zST
    do i = 1, grid%nelem
        elem => grid%elem(i)
        if (plus) then
          deg = elem%deg + state%space%plusDeg
          dof = DOFtriang( deg )
        else
          dof = elem%dof
        endif
        Tdof = elem%Tdof

        if ( associated(elem%zST) ) then
          print*, 'InitDualProblem:', size(elem%zST)
          stop
          deallocate(elem%zST)
        endif

        allocate( elem%zST(1:ndim,1:dof,1:Tdof), source = 0.0 )

    end do ! i

    if (state%time_dependent) &
      stop 'we need to compute IC for dual problem in InitDualProblem!'

   end subroutine InitDualProblem


  !> ! primal linear algebraic estimate for DWR method
  !> if the dual solution is not ready yet use norm of the residual
  subroutine computePrimalLinearDWRestimates( DWR, grid, primalRes)
    type( DWR_t ), intent(inout) :: DWR
    class( mesh ), intent(inout) :: grid
    real, dimension(:), intent(in) :: primalRes
    class( element ), pointer :: elem
    integer :: i, k, dof, Tdof, ndof, kvec, ivec, l
    integer :: nelem
    real :: estLinP
    real :: normZ
    real temp

    call state%cpuTime%startEstimTime()

    temp = 0.0 ! try DWR%x ~ zST
    estLinP = 0.0

    normZ = norm2( DWR%x)
    if (normZ < 1.E-15 ) then
      ! dual solution zST was not computed yet
      print*, 'PrimalLinearDWRestimates - zST is not ready, we use ||res|| instead'
      print*, 'Should it be multiplied somehow? '
      estLinP = norm2( primalRes)

    else
       nelem = grid%nelem
       ivec = 0
       do i=1,nelem
            elem => grid%elem(i)
            Tdof = elem%Tdof
            dof = elem%dof
            ndof = dof * Tdof * ndim

            kvec = ivec
            do l = 1, Tdof
               do k = 1, ndim
                  elem%eta(dwrLinP,1) = dot_product( &
                     primalRes(kvec+1:kvec+dof), grid%elem(i)%zST(k,1:dof,l) )
                  estLinP = estLinP + elem%eta(dwrLinP,1)
                  kvec = kvec + dof
               enddo !k
            enddo !l
            ivec = ivec + ndof
         end do !i
         temp = dot_product( primalRes(1:state%nsize), DWR%x(1:state%nsize) )
         if ( abs(temp - estLinP) > 1.E-9)  &
            print*, 'Difference between DWR%x and zST (should be the same):' , temp - estLinP
    endif

    DWR%estimLP = abs( estLinP )

    call state%cpuTime%addEstimTime()
    !print*, 'estLinP = ' , estLinP

  end subroutine computePrimalLinearDWRestimates

  !> dual linear algebraic estimate for the DWR method
  !> compute \eta_I^* = ( C^T*z_h - J(u_h) , u_h )
  subroutine computeDualLinearDWRestimates( DWR, grid)
    type( DWR_t ), intent(inout) :: DWR
    class( mesh ), intent(inout) :: grid
    class( element ), pointer :: elem
    integer :: i, k, dof, Tdof, ndof, kvec, ivec, l
    integer :: nelem
    real :: estLinD

    estLinD = 0.0

    nelem = grid%nelem
    ivec = 0
    do i=1,nelem
         elem => grid%elem(i)
         Tdof = elem%Tdof
         dof = elem%dof
         ndof = dof * Tdof * ndim

         kvec = ivec
         do l = 1, Tdof
            do k = 1, ndim
               elem%eta(dwrLinD,1) = dot_product( &
                  DWR%rr(kvec+1:kvec+dof), grid%elem(i)%wST(k,1:dof,l) )
!                  b(kvec+1:kvec+dof) - DWR%rhs(i)%x(k,1:dof,l), grid%elem(i)%wST(k,1:dof,l) )
               estLinD = estLinD + elem%eta(dwrLinD,1)
               kvec = kvec + dof
            enddo !k
         enddo !l
         ivec = ivec + ndof
      end do !i

      DWR%estimLD = abs( estLinD )
      !print*, 'estLinD = ' , estLinD

  end subroutine computeDualLinearDWRestimates




end module dwr_mod
