!> general error estimation subroutines
module estimates
  use ama_L2interpol
  use dwr_mod
  use dual_problem_mod
  use main_data  ! contains type(mesh) :: grid for computation
  use euler_problem
  use apost_estimation
  use project_estimation
  use solution_mod
  use eval_jumps

  implicit  none

  public:: computeNonlinDWRestimates
  public:: ComputeDWRestimates
  public:: DualDWRrezidErrorEstimates
  public:: DWRrezidErrorEstimates
  public:: FluxVectorDifference
  public:: PWpolynomialReconstComputeTerms
  public:: reconstructSolution
  public:: RezidErrorEstimates
  public:: SolutionDifference
  public:: SetVectorsFields
  public:: Set_Elem_Regul_Estim_Decay
  public:: RitzReconstruction


contains

  !> clear the arrays elem%estim_locL
  subroutine Clear_Elem_Estim_locL( )
    class(element), pointer :: elem
    integer :: i

    do i = 1, grid%nelem
       elem => grid%elem(i)
       elem%estim_locL = 0.
    enddo
  end subroutine Clear_Elem_Estim_locL

  !> perform the error estimates using the dual norm
  !> including (non-)linear algebraic error
  subroutine RezidErrorEstimates( onlyAS, Ttime_updated )
    logical, intent(in) :: onlyAS   ! only space and algebraic estimates
    logical, intent(in) :: Ttime_updated ! Ttime was already updated by %tau(1)
    class(element), pointer :: elem, elem1
    real, dimension(:), allocatable :: L_estim
    real, dimension(:,:), allocatable :: wi
    real, dimension(:,:,:), allocatable :: Dwi
    real :: machine_tol, rmax, t0, t1, t2, ttime, val, val1, weight
    integer :: i, j, k, ndof, ndofP, itest, imax, ipoc
    logical :: loc_implicitly

    !print*,'####  RezidErrorEstimates  start',onlyAS, Ttime_updated
    !!call cpu_time(t0)

    !itest = 360
    itest = -480

    ttime = state%time%ttime
    if(Ttime_updated) state%time%ttime = state%time%ttime - state%time%tau(1)

    allocate(L_estim(1:max_eta) )
    loc_implicitly = state%nlSolver%implicitly

    state%nlSolver%implicitly = .false.
    grid%elem(:)%deg_plus = .true.


    ! setting of fields elem%vec(rhs,*), elem%vec(rhsT,*) for error estimate
    call SetVectorsFields( onlyAS )

    call cpu_time(t0)


    state%estim(max_resT_S,:) = 0.
    state%estim(min_resT_S,:) =  1E+50
    state%estim(min_resT_S_loc,:) =  1E+50
    state%estim(resA_ST,:) = 1E+50

    rmax = 0.

    call cpu_time(t1)

    L_estim(:) = 0.   ! total value of the residuum

    do i=1, grid%nelem
       elem => grid%elem(i)
       ! NOT MULTIPLIED,1/tau in project.f90 removed
       !elem%vec(rhs,:) = elem%vec(rhs,:) * state%time%tau(1)

       ! the following should be performed for the STDGM approach
       !elem%vec(rhsT,:) = elem%vec(rhsT,:) / state%time%tau(1)

       !!!call EnergyReziduumElemEstimate(elem)  ! element residuum
       !do j=1,4

       if( state%time%disc_time /= 'STDG') then
          call DualElemEstimate(elem, 3, onlyAS)  ! 3 => element residuum in X- norm
       else
          !call ST_DualElemEstimate_Var2(elem,3)
          !call ST_DualElemEstimate(elem, 2 )  ! 2 => element residuum in the H^1-norm (for inviscid??)
          call ST_DualElemEstimate(elem, 3 )  ! 3 => element residuum in X- norm
       endif

       !if(elem%i < 10 .or. elem%i == 1234) &
       !     write(*,'(a10, 2i5, 300es12.4)') 'etas:',elem%i, elem%dof,  elem%eta(1:5,:)

       ! limitation along the shock waves
       val = 1.
       if(state%type_IC == 8) val = 2*elem%area/elem%diam  ! .or. abs(elem%xc(1) - 1.) > 2E-2) then

       !if(val > 1E-2) then
       if(state%type_IC /= 8  .or. abs(elem%xc(1) - 1.) > 5E-2) then

          L_estim( resA) = L_estim( resA) + elem%eta(resA, 1)**2
          L_estim( resS) = L_estim( resS) + elem%eta(resS, 1)**2
          L_estim( resT) = L_estim( resT) + elem%eta(resT, 1)**2
          L_estim(resST) = L_estim(resST) + elem%eta(resST, 1)**2
          L_estim(resSr) = L_estim(resSr) + elem%eta(resSr, 1)**2

       endif

       ! Verfurth approach
       !print*,'$$$', elem%eta(resT, 1), elem%eta(resS, 1)
       if(elem%eta(resS, 1) > 0.) &
            state%estim(min_resT_S_loc,1) = min(state%estim(min_resT_S_loc,1),  &
            elem%eta(resT, 1) /  elem%eta(resS, 1) )

       !if(elem%eta(resT, 1) > rmax) then
       !   rmax = elem%eta(resT, 1)
       !   imax = i
       !endif

       state%estim(max_resT_S,1) = max(state%estim(max_resT_S,1),  &
            elem%eta(resT, 1) / max(1E-15, elem%eta(resS, 1)) )
       !!state%estim(resA_ST,:) = 1E+50

    enddo

    !print*, 'stopped in estimates.f90 after ST_DualElemEstimate'
    ! stop


    state%L_estim(1:max_eta) = sqrt(L_estim(1:max_eta))


    machine_tol = 1.E-01
    state%estim(resA_ST_loc,1) = 0.

    ! local algebraic criterion
    do i=1,grid%nelem
       elem => grid%elem(i)

       if( elem%eta(resST, 1) >  machine_tol * state%L_estim(resST) / grid%nelem**0.5 ) then

          state%estim(resA_ST_loc,1) = max(state%estim(resA_ST_loc,1), &
               elem%eta(resA, 1)/ elem%eta(resST, 1))

          !write(198,*) elem%xc(:), elem%eta(resA, 1), elem%eta(resST, 1) ,elem%eta(resA, 1)/ elem%eta(resST, 1), elem%i
       endif
    enddo


    ! steady-state approach
    !print*, 'here' , state%L_estim(resA) , state%L_estim(resS)
    state%estim(resA_ST,1) = state%L_estim(resA) /  state%L_estim(resS)

    ! STDG approach
    state%estim(min_resT_S, 1) = state%L_estim(resT) / state%L_estim(resS)
    !print*,'#### WERTY',state%estim(min_resT_S, 1) , state%estim(min_resT_S_loc,1), state%estim(max_resT_S,1)

    grid%elem(:)%deg_plus = .false.

    state%nlSolver%implicitly  = loc_implicitly

    !write(*,'(a10,i5, 20es12.4)') '##RDE342ed', state%nlSolver%iter, state%estim(1:9, 1)
    !write(*,'(a10,i5, 20es12.4)') '##RDE342ed', state%nlSolver%iter, state%L_estim(1:4)

    deallocate(L_estim)

    call cpu_time(t2)
    state%CPU_constaint = state%CPU_constaint + t2 - t1

    !write(*,'(a10,i5, 20es12.4)') '##RDE342ed', state%nlSolver%iter, state%estim(1:9, 1)
    !write(*,'(a10,i5, 20es12.4)') '##RDE342ed', state%nlSolver%iter, state%L_estim(1:4)

    ! computing of quantities for the mesh refinement
    if(.not. onlyAS) then

       !call JumpsEvaluation( )

       !if(.not. grid%ElemSupports) & ! create the list of elements sharing at least a vertex with elem
       !     call SeekElemSupports(grid)

       ! averaging of the estimate
       ! setting of 'elem%estim_loc' including elem%eta(resST, 1) of neighbours elements
       ! elem%estim_loc**2 = sum_{K\in N(K)} elem%eta(resST, 1)**2
       do i=1,grid%nelem
          elem => grid%elem(i)
          if( state%time%disc_time /= 'STDG') then
             elem%estim_loc = elem%eta(resS, 1)**2    ! Verfurth approach
          else
             if(state%time_dependent) then  !!!  VERIFY ????
                elem%estim_loc = elem%eta(resST, 1)**2    ! STDGM approach
             else
                elem%estim_loc = elem%eta(resST, 1)**2 /state%time%tau(1)   ! STDGM approach
             endif

          endif

          ipoc = 0

          !if(i == 1) write(*,'(a8,20es12.4)') '#@#@@#@',elem%estim_loc

          weight = 0.0
          if(weight > 0.) then
             ! only neighbouring elements
             ! FOR HP_STAEDY
             do j=1,elem%flen
                k = elem%face(neigh,j)

                ! all elements sharing at least a vertex
                !  FOR ST_ESTIMS
                !do j=1,elem%isupp
                !   k = elem%supp(j,1)

                if(k > 0) then
                   elem1 => grid%elem(k)
                   if( state%time%disc_time /= 'STDG') then
                      elem%estim_loc = elem%estim_loc + weight*elem1%eta(resS, 1)**2  ! Verfurth
                      !elem%estim_loc = max(elem%estim_loc , elem1%eta(resS, 1)**2)   ! Verfurth
                   else
                      elem%estim_loc = elem%estim_loc + weight*elem1%eta(resST, 1)**2   ! STDGM approach
                   endif

                   ipoc = ipoc + 1
                endif
             enddo
             ! !!  elem%estim_loc = elem%estim_loc**0.5  ! we store the square
             elem%estim_loc = elem%estim_loc / (1. + weight*ipoc)
          endif !if(weight > 0.) then

          !if(i == 1) write(*,'(a8,20es12.4)') '#@#@@#@2',elem%estim_loc,elem%estim_locL

          !!write(*,'(a6,i5,6es12.4)')'EST:xd',elem%i, elem%eta(resS, 1), elem%estim_loc, 1.*ipoc, (1. + weight*ipoc)
          elem%estim_loc = sqrt( elem%estim_loc)  ! + elem%jumpsJh)  !!elem%rezid)

          ! storing of several time levs - IF USED MUST BE DONE IN DIFFERENT WAY
          ! OUTSIDE OF  RezidErrorEstimates
          ! ALREADY DONE IN COMPUTEad.F90, BUT NOT TESTED !!!!
          !elem%estim_locL = sqrt( elem%estim_locL**2 + elem%estim_loc**2)
          !elem%estim_loc = elem%estim_locL  ! used for adaptivity

          !if(i == 1) write(*,'(a8,20es12.4)') '#@#@@#@3',elem%estim_loc

       enddo ! do i=1,grid%nelem


       !val = 0.; val1 = 0.
       !do i = 1, grid%nelem
       !  elem => grid%elem(i)
       !  val  = val  + elem%estim_loc**2
       !  val1 = val1 + elem%eta(resST,1)**2
       !enddo


       ! will be deallocated at the end of grid  of arrays allocated in SeekElemSupports
       !do i=1,grid%nelem
       !   elem => grid%elem(i)
       !   deallocate(elem%supp)
       !enddo

    endif

    state%time%Ttime = ttime

    !print*,'####  RezidErrorEstimates  END'

    !if( state%time%cn ) then
    !   do i=1,grid%nelem
    !      elem => grid%elem(i)
    !   enddo
    !endif

    call cpu_time(t2)
    !write(*,'(a40, 2l6, f12.4)') &
    !      '#CPU#  RezidErrorEstimates  ends',onlyAS, Ttime_updated, t2 - t0

    !do k=1, 5 ! grid%nelem
    !   write(*,'(a8, 2i5,30es12.4)') 'est_Loc:', grid%elem(k)%i,grid%elem(k)%dof, &
    !        grid%elem(k)%estim_loc, grid%elem(k)%eta(1:4, 1)
    !enddo

  end subroutine RezidErrorEstimates

  !> setting of field elem%vec(rhs,*), elem%vec(rhsT,*) for error estimate
  subroutine SetVectorsFields(onlyAS )
    logical, intent(in) :: onlyAS   ! only space and algebraic estimates
    integer :: ityp, ityp1, itest, i
    class(element), pointer :: elem, elem1
    logical :: deg_plus

    deg_plus = .true.
    !itest = 25

    !ityp = 1  ! backward Euler, STDGM approach
    !ityp = 2  ! Crank-Nicolson, STDGM approach
    !ityp = 3  ! backward Euler, Verfurth approach
    !ityp = 4 ! Crank-Nicolson, Verfurth approach

    ityp1 = 0    ! pw polynomial  approximation in time
    !ityp1 = 1   !STDGM approach
    !ityp1 = 2   !Verfurth approach

    if( state%time%cn ) then
       ityp = 2* ityp1
    else
       ityp = 2* ityp1 - 1
    endif

    !state%num_flux = .false.      ! use physical fluxes instead of numerical ones


    !print*,'SetVectorFields, ityp == ',ityp
    if( state%time%disc_time == 'STDG') then

       !if(state%time%iter <= 1) print*,'####, ATTENTION in estimates.f90, ttime was updated?', state%time%ttime
!         print*, 'SetVectorsFields calling ComputeSTDGM_Terms with implicitly = ', state%nlSolver%implicitly

       !print *,'Already DONE'

       !call ComputeSTDGM_Terms( deg_plus )

       ! do i=1, 5 !grid%nelem
       !    elem => grid%elem(i)
       !    write(*,'(a10, 2i5, 300es12.4)') 'b_sol:',elem%i, elem%dof,  elem%wST(:,1:elem%dof,:)
       ! enddo


    elseif(ityp == -1) then  ! backward Euler, piecewise linear time reconstruction
       ! PREDELAT !!!! VD

       !itest = 25
       !write(*,'(a6,i5,16es12.4)') 'w0', itest, 0., grid%elem(itest)%w(0, 1:2), &
       !     grid%elem(itest)%w(1, 1:2),grid%elem(itest)%w(2, 1:2), grid%elem(itest)%w(3, 1:2)
       !print*,'-----------, vector'

       call PWpolynomialReconstComputeTerms(onlyAS ) !flux vector:  %vec(rhs, * ) = \int_{I_m} c_h(w_h, * )

    elseif(ityp == 1) then  ! backward Euler, STDGM approach
       call ComputeTerms(deg_plus )     ! flux vector:  %vec(rhs, * ) = c_h(w_h, * )

       call TimeDerivativeVector( ) ! BDF term: %vec(rhs,*) += \sum_n \alpha_n (w^{k-n},*)/ tau

       call SolutionDifference( ) ! for estimT:  %vec(rhsT, * ) = (w^{k} - w^{k-1}, *)

    elseif( ityp==0 .or. ityp == 2 ) then ! Crank-Nicolson

       call CrankNicolsonVectorFlux( ) ! %vec(rhs, * ) = ( c_h(w_h^k, * ) + c_h(w_h^{k-1}, * ))/2

       call TimeDerivativeVector( ) ! BDF term: %vec(rhs,*) += \sum_n \alpha_n (w^{k-n},*)/ tau

       call SolutionDifference( ) ! for estimT:  %vec(rhsT, * ) = (w^{k} - w^{k-1}, *)


    elseif(ityp == 3) then ! backward Euler, Verfurth approach
       !elem => grid%elem(itest)
       !!call PWpolynomialReconstComputeTerms( ) !flux vector:  %vec(rhs, * ) = \int_{I_m} c_h(w_h, * )

       !write(197,150) 'PW rhsT:',itest, elem%vec(rhsT, :)
       !write(197,150) 'PW rhs :',itest, elem%vec(rhs , :)
       !!call DualElemEstimate(elem, 3)
       !write(197,'(a6,20es12.4)') 'etas:',elem%eta(resA, 1), elem%eta(resS, 1), elem%eta(resT, 1), elem%eta(resST, 1)
       !write(197,*) '--------------------'


       ! vec(rhsT, * ) = 1/\tau_k \int_k < F(w_h) - F(\tilde{w}_h), \phi_i > dt
       ! vec(rhs, * ) = 1//\tau_m (w_h^k - w_h^{k-1}, \phi) + < F(w_h) , \phi_i >

       !approach based on Verfurth, MUST BE CALLED BEFORE SPACE ESTIMATE !!!
       call FluxVectorDifference( ) ! for estimT:  %vec(rhsT, * ) = c_h(w_h, * ) - c_h(W_h, * )

       ! vec(rhs, * ) already evaluated in FluxVectorDifference( )
       !!!call ComputeTerms( )     ! flux vector:  %vec(rhs, * ) = c_h(w_h, * )

       call TimeDerivativeVector( ) ! BDF term: %vec(rhs,*) += \sum_n \alpha_n (w^{k-n},*)/ tau

       !write(197,150) 'Ver rhsT:',itest, elem%vec(rhsT, :)
       !write(197,150) 'Ver rhs :',itest, elem%vec(rhs , :)
       !call DualElemEstimate(elem, 3)
       !write(197,'(a6,20es12.4)') 'etas:',elem%eta(resA, 1), elem%eta(resS, 1), elem%eta(resT, 1), elem%eta(resST, 1)
       !write(197,*) '--------------------'

150 format(a10,i5, 4(3es12.4 ' |',3es12.4,' ## ') )

    elseif(ityp == 4) then ! Crank-Nicolson, Verfurth approach
       ! vec(rhsT, * ) = 1/\tau_k \int_k < F(w_h) - F(\tilde{w}_h), \phi_i > dt
       ! vec(rhs, * ) = 1//\tau_m (w_h^k - w_h^{k-1}, \phi) + < F(w_h) , \phi_i >

       !approach based on Verfurth, MUST BE CALLED BEFORE SPACE ESTIMATE !!!
       call FluxVectorDifference( ) ! for estimT:  %vec(rhsT, * ) = c_h(w_h, * ) - c_h(W_h, * )

       call CrankNicolsonVectorFlux( ) ! %vec(rhs, * ) = ( c_h(w_h^k, * ) + c_h(w_h^{k-1}, * ))/2

       call TimeDerivativeVector( ) ! BDF term: %vec(rhs,*) += \sum_n \alpha_n (w^{k-n},*)/ tau

    else
       print*,'Unknown ityp in subroutine SetVectorsFields, estimated.f90'
       stop

    endif

    !do i=1,1 !grid%nelem
    !!do i=18,18
    !   elem => grid%elem(i)
    !   !write(*,'(a6,i5, 80es12.4)') 'rhsT:',i, elem%vec(rhsT, :)
    !   write(*,'(a6,i5, 80es12.4)') 'rhs :',i, elem%vec(rhs , :)
    !enddo



  end subroutine SetVectorsFields

  !> reconstruct piecewise linear in time solution
  !> \f$\tilde{w}(t) = w^{k-1} + \frac{t - t_{k-1}}{\tau_k} (w^k - w^{k-1})\f$
  !> \f$ \%vec(rhs, \varphi ) = \int_{I_k} c_h(\tilde{w}(t), \varphi ) d t\f$,
  !> \f$ \varphi \in S_{hp}^{+}\f$
  subroutine PWpolynomialReconstComputeTerms(onlyAS )
    logical, intent(in) :: onlyAS   ! only space and algebraic estimates
    class(element), pointer :: elem
    type(Gauss_rule), pointer :: G_rule
    class(Time_rule), pointer :: T_rule
    real, dimension(:,:), allocatable :: Lag_coef
    integer :: time_deg, Tdeg
    integer :: Gnum, Gdof , j, i, l, l1, ndof, ndofP, itest
    real :: t, ctime_store
    logical ::  deg_plus

    deg_plus = .true.
!!! TWO CHANGES, 1) time_deg = Tdeg , 2)  Gnum = time_deg  !!!!!!!!!

    if(.not. onlyAS) then
       ctime_store = state%time%ctime

       Tdeg = state%time%deg_actual  ! degree of actual BDF rule
       allocate(Lag_coef(0:1, 0:Tdeg) )  ! index1 = 0 -> Lagr functiom, =1 -> its derivative


       !time_deg = 3   ! maximal degree of test functions in time
       time_deg = Tdeg   ! maximal degree of test functions in time  ????


       !itest = 360
       itest = 25

       ! rule of integration in time
       !!Gnum = time_deg  + 1  ! ???
       Gnum = time_deg   ! ???

       !TODO why are we using G_rule - may not work if T_rule is Radau type
       G_rule => state%space%G_rule(Gnum)
       T_rule => state%time%T_rule(Gnum)



       Gdof = G_rule%Qdof

       ! preparing of arrays
       do i=1,grid%nelem
          elem => grid%elem(i)
          ndof = elem%dof * ndim
          ndofP = elem%dof_plus * ndim
          ! indexes of elem%wS: i= -Tdeg:0 ... storing of w^k, w^{k-1}, ..., w^{k-Tdeg},
          ! i=1..time_deg  computing of c(w_h, \vp_*) \phi_i
          allocate(elem%wS( -Tdeg:time_deg, 1:ndofP ) )

          ! storring of the actual solution
          do j=0, Tdeg
             elem%wS(-j, 1:ndof) = elem%w(j, 1:ndof)
          end do

          elem%wS( 1:time_deg, 1:ndofP) = 0. ! arrays for computing of c(w_h, \vp_*) \phi_i
          elem%vec(rhs, 1:ndofP) = 0.
       enddo

       ! integration over the time interval
       do j=1, Gdof
          ! actual time
          t = G_rule%lambda(j)
          state%time%ctime = state%time%ttime + state%time%tau(1) * t   ! %ttime is not yet updated

          call SetLagrCoeffs(Tdeg,  t, Lag_coef(0:1, 0:Tdeg) )

          ! pw polynomial reconstruction of w and dw/dt at t
          do i=1,grid%nelem
             elem => grid%elem(i)
             ndof = elem%dof * ndim

             ! evaluation of w and dw/dt at t using the Lagrangian interpolation
             elem%w(0,1:ndof) = elem%wS(-Tdeg , 1:ndof) * Lag_coef(0, 0)  ! w
             elem%w(1,1:ndof) = elem%wS(-Tdeg , 1:ndof) * Lag_coef(1, 0)  ! dw/dt
             do l=1, Tdeg
                elem%w(0,1:ndof) = elem%w(0,1:ndof) + elem%wS(-Tdeg + l, 1:ndof) * Lag_coef(0, l)
                elem%w(1,1:ndof) = elem%w(1,1:ndof) + elem%wS(-Tdeg + l, 1:ndof) * Lag_coef(1, l)
             enddo

          enddo

          !elem => grid%elem(itest)
          !write(*,*) state%time%ttime+state%time%tau(1), elem%vec(rhs,1:6)

          !print*,'########## bef ComuteTerms', state%time%iter
          !print*,'CHECK ComputeTerms(deg_plus )'
          call ComputeTerms(deg_plus )      ! flux vector:  %vec(rhs, * ) = c_h(w_h, * )
          !print*,'########## aft ComuteTerms', state%time%iter

          !elem%vec(rhsT,:) = elem%vec(rhs,1:6)

          !! REMOVE
          !do i=1,grid%nelem
          !   elem => grid%elem(i)
          !   elem%vec(rhs,:) = 0.
          !enddo

          call AddTimeDerivativeVector( ) ! %vec(rhs, * ) = %vec(rhs, * ) + (D_t w, *)

          !elem => grid%elem(itest)
          !write(283,*) state%time%ttime+state%time%tau(1), elem%vec(rhs,1:6)
          !write(284,*) state%time%ttime+state%time%tau(1), -(elem%vec(rhs,1:6)-elem%vec(rhsT,1:6))


          ! summning over time integ. nodes for each element
          do i=1,grid%nelem
             elem => grid%elem(i)
             ndofP = elem%dof_plus * ndim


             !if(i == itest) then
             !   write(791,'(a6,30es14.6)') 'new', elem%vec(rhs, 1:ndofP)
             !   write(791,'(a6,30es14.6)') 'DwO', (elem%wS(0, 1:) - elem%wS(-1, 1:) )/state%time%tau(1)
             !   write(791,'(a6,30es14.6)') 'DwN', elem%w(1, 1:)
             !   write(791,*) '--------------------------------'
             !endif

             !! adding of the BDF term: %vec(rhs,*) += \sum_n \alpha_n (w^{k-n},*)/ tau
             !elem%vec(rhs, 1:ndofP) = elem%vec(rhs, 1:ndofP) + elem%vec(rhsT, 1:ndofP)

             do l = 1, time_deg   ! c_h(w_h, \varphi_{1:ndofP}) \phi_l  in integ node "j"
                elem%ws(l, 1:ndofP) =  elem%ws(l, 1:ndofP) &
                     + G_rule%weights(j) * elem%vec(rhs, 1:ndofP) * T_rule%phi(l, j)
             enddo

             !if(i == itest) then
             !   l = 6
             !   write(186,'(a8,3i5, 20es14.6)') 'wk..',state%time%iter, 0,0, elem%ws(0,1:6)
             !   write(186,'(a8,3i5, 20es14.6)') 'wk-1',state%time%iter, 0,0, elem%ws(-1,1:6)
             !   write(186,'(a8,3i5, 20es14.6)') 'wi  :', state%time%iter, itest, 0, elem%w(0,1:l)
             !   write(186,'(a8,3i5, 20es14.6)') 'wi_t:', state%time%iter, 0,  0, elem%w(1,1:l)
             !   do l1=1,time_deg
             !      write(186,'(a8,3i5, 20es14.6)') 'ws i:', state%time%iter, j, l1, elem%ws(l1, 1:l)
             !   end do
             !   write(186,*) '--------------------'
             !   write(187,*) state%time%iter + G_rule%lambda(j), elem%ws(1, 1:6)
             !   write(188,*) state%time%ttime+state%time%tau(1), elem%ws(0,1:6)
             !   write(189,*) state%time%ttime, elem%ws(-1,1:6)
             !   write(190,*) state%time%ttime+state%time%tau(1)* G_rule%lambda(j), elem%w(0,1:6)

             !endif

          enddo ! i=1,grid%nelem
       enddo ! j=1,Gdof
       !do not mupltiply with the size of the time step, done in the totat summing in compute.f90


       !write(187,*) '#####    '

       !elem => grid%elem(itest)
       !write(381,*) state%time%ttime+state%time%tau(1), elem%ws(1,1:6)
       !write(382,*) state%time%ttime+state%time%tau(1), elem%ws(2,1:6)
       !write(383,*) state%time%ttime+state%time%tau(1), elem%ws(3,1:6)


       ! refreshing of the actual solution
       ! time derivative in elem%vec(rhsT,:) need not be stored any more
       do i=1,grid%nelem
          elem => grid%elem(i)
          ndof   = elem%dof * ndim
          ndofP = elem%dof_plus * ndim

          elem%w(0, 1:ndof) = elem%wS( 0, 1:ndof)  ! restoring of the solution from last 2 levels
          elem%w(1, 1:ndof) = elem%wS(-1, 1:ndof)

          ! summing of squares of T_rule%phi(l,:), phi_l(t) are L^2 orthogonal
          elem%vec(rhsT, 1:ndofP) =  0.
          do l=1, time_deg
             elem%vec(rhsT, 1:ndofP) = elem%vec(rhsT, 1:ndofP) + elem%wS(l, 1:ndofP)**2

             !if(i == itest) write(19,'(a3,i5,30es14.6)') '!!!', l, &
             !     VectorNorm( elem%vec(rhsT, 1:ndofP) ), &
             !     VectorNorm( elem%wS(l, 1:ndofP) ), elem%wS(l, 1:ndofP)

          enddo

          elem%vec(rhsT, 1:ndofP) = elem%vec(rhsT, 1:ndofP)**0.5

          !!if(i == itest) write(*,'(a3,30es14.6)') '!!!',elem%vec(rhsT, 1:6) - elem%vec(rhs, 1:6)
          !if(i == itest) write(19,*)'-------------------------------'

          deallocate(elem%wS)
       enddo


       deallocate(Lag_coef)
       state%time%ctime = ctime_store
    endif  ! if (.not. onlyAS)

    ! space  estimate
    call ComputeTerms(deg_plus )     ! flux vector:  %vec(rhs, * ) = c_h(w_h, * )

    call TimeDerivativeVector( ) ! BDF term: %vec(rhs,*) += \sum_n \alpha_n (w^{k-n},*)/ tau

    itest = -981
    if(itest > 0 .and. itest <= grid%nelem) then
       elem => grid%elem(itest)

       !write(*,'(a4,i5,8es12.4)') &
       !     'PWS',elem%i, elem%w(0:3, 1), elem%w(0:3, 2)

       write(*,'(a4,i5, 60es12.4)') &
            'PWS',elem%i, elem%vec(rhs,:)

       write(*,'(a4,i5, 60es12.4)') &
            'PWT',elem%i, elem%vec(rhsT,:)
    endif




    !elem => grid%elem(itest)
    !write(63,'(a5,200es12.4)') 'rhs',elem%vec(rhs, :)
    !write(63,'(a5,200es12.4)') 'rhsT',elem%vec(rhsT, :)
    !write(63,*) '/////////////////'
    !stop
  end subroutine PWpolynomialReconstComputeTerms

  !> evaluate   elem%vec(rhs, i ) =
  !> \f$ \frac12 \left(c_h(w_h^k, \varphi_i ) + c_h(w_h^k, \varphi_i ) \right)\f$
  !> for \f$\varphi_i \in S_{hp}\f$
  subroutine CrankNicolsonVectorFlux( )
    class(element), pointer :: elem
    integer ::  i, j, k, ndof, ndofP
    logical :: deg_plus

    deg_plus = .true.

    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof  = elem%dof * ndim
       ndofP = elem%dof_plus * ndim

       allocate(elem%wS(1:1,1:elem%dof * ndim) )

       elem%wS(1,1:elem%dof * ndim)  = elem%w(0,1:ndof)

       elem%w(0,1:elem%dof * ndim) = elem%w(1, 1:ndof)
    enddo

    call ComputeTerms(deg_plus )     ! flux vector:  %vec(rhs, * ) = c_h(w_h^{k-1}, * )
    do i=1,grid%nelem
       elem => grid%elem(i)
       elem%vec(rhsOLD,1:ndofP)  = elem%vec(rhs,1:ndofP)

       elem%w(0,1:ndof) = elem%wS(1,1:ndof)
       deallocate(elem%wS )
    enddo

    call ComputeTerms(deg_plus )     ! flux vector:  %vec(rhs, * ) = c_h(w_h^{k }, * )

    do i=1,grid%nelem        ! Crank-Nicolson (F(w^k) + F(w^{k-1} ) /2
       elem => grid%elem(i)
       elem%vec(rhs,1:ndofP)  = (elem%vec(rhs,1:ndofP) + elem%vec(rhsOLD,1:ndofP) )/2.
    enddo
  end subroutine CrankNicolsonVectorFlux



  !> compute the vector: elem\%vec(i) =
  !> \f$ \frac{1}{\tau_m} \int_{I_m}
  !> ( F({\bf w}_h) - F(\tilde{{\bf w}}_h), \varphi_i)_{L^2(K)} \, {\rm d}t\f$,
  !> \f$i=1,\dots, DOF^+_K,\  K \f$ = elem \f$ \in T_h\f$,
  !> \f$ {\bf w} \f$ is piecewise constant in time
  !> \f$ \tilde{{\bf w}} \f$ is piecewise affine in time
  subroutine FluxVectorDifference( )
    class(element), pointer :: elem      ! elem = element
    type(Gauss_rule), pointer :: G_rule
    real, dimension(:,:), allocatable :: wk_wk1
    integer :: i, j, dof, ndof, ndofP, Qdof, Gnum, Gdof
    real :: t
    integer :: itest
    logical :: deg_plus

    deg_plus = .true.

    itest = -360


    ! saving of elem%w(0,*)
    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof  = elem%dof * ndim

       allocate(elem%wS(1:1,1:elem%dof * ndim) )
       elem%wS(1, 1:ndof )  = elem%w(0,1:ndof)

       ndofP = elem%dof_plus * ndim
       elem%vec(rhsT, 1:ndofP) = 0.
    enddo

    ! integration in time
    Gnum = 2   ! 15 is the maximal one !!!!!!!!!
    G_rule => state%space%G_rule(Gnum)
    Gdof = G_rule%Qdof

    do j=1, Gdof
       ! actual time
       t = G_rule%lambda(j)
       !!!!t = 1.

       state%time%ctime = state%time%ttime + state%time%tau(1) * t   ! %ttime is not yet updated

       do i=1,grid%nelem
          elem => grid%elem(i)
          ndof  = elem%dof * ndim
          ! piecewise affine reconstructin in time
          elem%w(0,1:ndof) = (1. - t) * elem%w(1,1:ndof) + t * elem%wS(1, 1:ndof)

       enddo

       call ComputeTerms( deg_plus)     ! flux vector:  %vec(rhs, * ) = c_h(w_h^{k-1}, * )

       do i=1,grid%nelem
          elem => grid%elem(i)
          ndofP = elem%dof_plus * ndim

           if(i == itest)  write(*,'(i5,200es14.6)') j-1,elem%vec(rhsT, 13:15),elem%vec(rhs, 13:15)!,elem%w(0,1:ndof)

          elem%vec(rhsT, 1:ndofP) = elem%vec(rhsT, 1:ndofP) &   ! weighted average
               + G_rule%weights(j) * elem%vec(rhs, 1:ndofP)
           if(i == itest)  write(*,'(i5,200es14.6)') j,elem%vec(rhsT, 13:15),elem%vec(rhs, 13:15)!,elem%w(0,1:ndof)
       enddo

    enddo

    ! backward setting of the actual solution to elem%w
    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof  = elem%dof * ndim
       elem%w(0,1:ndof) = elem%wS(1,1:ndof)
       deallocate(elem%wS )

    enddo


    ! for t= t^k
    call ComputeTerms(deg_plus )     ! flux vector:  %vec(rhs, * ) = c_h(w_h^{k-1}, * )

    do i=1,grid%nelem
       elem => grid%elem(i)
       ndofP = elem%dof_plus * ndim

       elem%vec(rhsT, 1:ndofP) = elem%vec(rhsT, 1:ndofP)  - elem%vec(rhs, 1:ndofP)
       if(i == itest)  write(*,'(i5,200es14.6)') 88,elem%vec(rhs, 13:15)
       if(i == itest)  write(*,'(i5,200es14.6)') 99,elem%vec(rhsT, 13:15)

       !write(*,'(a6,i5, 20es12.4)') 'rhsT:',i, elem%vec(rhsT, 1:ndofP)
       !write(*,'(a6,i5, 20es12.4)') 'rhsT:',i, elem%vec(rhs , 1:ndofP)
       !print*
    enddo
    !!print*,'^^^^^^^^^^^^^^^^^^^^ Ver'

  end subroutine FluxVectorDifference

 !> compute the vector: elem\%vec(i) =
 !> \f$ ({\bf w}_h^k - {\bf w}_h^{k-1}, \varphi_i)_{L^2(K)} \f$
 !> \f$i=1,\dots, DOF^+_K,\  K \f$ = elem \f$ \in T_h\f$
  subroutine SolutionDifference( )
    class(element), pointer :: elem      ! elem = element
    real, dimension(:,:), allocatable :: wk_wk1
    integer :: i, dof, dofA, Qdof

    do i=1,grid%nelem
       elem => grid%elem(i)

       dofA = dof
       if(elem%deg_plus) dofA = elem%dof_plus
       Qdof = elem%Qdof

       elem%vec(rhsT,:) = 0.
       allocate( wk_wk1(1:Qdof, 1:ndim) )

       call Eval_w_w_Elem(elem, wk_wk1)
       call EvalVectorB(elem, wk_wk1, dofA, elem%vec(rhsT,1:dofA*ndim) )
       ! NOT : "normalization:, in order to be in agreement with elem%vec(rhs, 1:dofA*ndim)
       !!elem%vec(rhsT,1:dofA*ndim) = elem%vec(rhsT,1:dofA*ndim) /state%time%tau(1)

       !if( i <= 2) write(*,'(a6,i5, 25es12.4)') 'w-w:',i, &
       !     wk_wk1(:,1), 999999., elem%vec(rhsT,1:dofA*ndim)

       deallocate(wk_wk1)
    enddo
    !print*,'_________________________________________________________'

  end subroutine SolutionDifference

  !> reconstruct the solution wST/zST_plus from the space P^{p+1}
  !> using the LeastSquareInterpolationWeighted, it uses the array elem%wS
  !> primal == .true. -> wST , .false. -> zST
  !> R_type = -1 -> H1 interpolation, 0 -> L2 interpolation
  subroutine reconstructSolution( grid, primal, R_type )
    class( mesh ), intent(inout) :: grid
    logical, intent(in) :: primal
    integer, intent(in) :: R_type
    integer :: nelem, degP, dofP, Qnum, i
    class( element ), pointer :: elem
    real, allocatable, dimension(:,:,:,:) :: Dw ! for the least square reconstruction
    logical :: flag ! problem flag


    flag = .true.

    if(.not. grid%ElemSupports) & ! create the list of elements sharing at least a vertex with elem
            call SeekElemSupports(grid)

    nelem = grid%nelem

    if (R_type == 2) then
      !need Aplus, allocated zST(1:dof_plus) , res(vp+) forall vp+ \in Shpp
      !update matrix and rhsvector
      ! primal - rhs = res(u_h)(vp+) vp+ \in Shpp
      ! dual   - rhs = res*(z_h)(vp+) vp+ \in Shpp
      ! saved to wSTplus or zSTplus the reconstruction for all elements
      call RitzReconstruction( grid, primal, flag, DWR )

    else

       if (state%time%deg > 0) &
         stop 'reconstructSolution works only for stationary problems - q = 0'

       if (primal) then
         do i = 1, nelem
            elem => grid%elem(i)
            ! rewrites - need no allocation check
            elem%wS = Transfer_funST_to_funS( elem%wST( 1:ndim, 1:elem%dof, 1:elem%Tdof ) , &
                      elem%dof, elem%Tdof, 0, elem%TQnum )
          enddo !i
       else
         do i = 1, nelem
            elem => grid%elem(i)

            elem%wS = Transfer_funST_to_funS( elem%zST( 1:ndim, 1:elem%dof, 1:elem%Tdof ) , &
                      elem%dof, elem%Tdof, 0, elem%TQnum )
          enddo !i
       endif

       !primal solution
       if (primal) then
       ! reconstruct the solution I_h(w)
          do i = 1, nelem
            elem => grid%elem(i)
            degP = elem%deg + state%space%plusDeg  ! ?
            dofP = DOFtriang(degP) ! ?
            Qnum = elem%Qnum

            ! only 1:ndim,0,0,1:dofP is used - the others are for the derivatives
            allocate( Dw(1:ndim, 0:degP, 0:degP, 1:dofP ) , source = 0.0 )
            call LeastSquareInterpolationWeighted(elem, ndim, .false., Qnum, degP, dofP, Dw, R_type  )

            if ( associated( elem%wSTplus ) ) then
               deallocate( elem%wSTplus )
            endif
            allocate( elem%wSTplus( 1:ndim, 1:dofP, 1:elem%Tdof ) )
            elem%wSTplus( 1:ndim, 1:dofP, 1 ) = Dw(1:ndim, 0, 0, 1:dofP)

            !CONTROL - not zero somewhere
            if ( norm2(Dw(1, 0, 0, elem%dof + 1 : dofP) ) > 1E-9 ) then
               flag = .false.
            endif

            deallocate( Dw )
          enddo !i

       !dual
       else
          ! reconstruct the solutions I_h(z)
          do i = 1, nelem
            elem => grid%elem(i)
            degP = elem%deg + state%space%plusDeg  ! ?
            dofP = DOFtriang(degP) ! ?
            Qnum = elem%Qnum

            ! only 1:ndim,0,0,1:dofP is used - the others are for the derivatives
            allocate( Dw(1:ndim, 0:degP, 0:degP, 1:dofP ) , source = 0.0 )
            call LeastSquareInterpolationWeighted(elem, ndim, .false., Qnum, degP, dofP, Dw, R_type  )

            if ( associated( elem%zSTplus ) ) then
               deallocate( elem%zSTplus )
            endif
            allocate( elem%zSTplus( 1:ndim, 1:dofP, 1:elem%Tdof ) )
            elem%zSTplus( 1:ndim, 1:dofP, 1 ) = Dw(1:ndim, 0, 0, 1:dofP)

            !CONTROL - not zero somewhere
            if ( norm2(Dw(1, 0, 0, elem%dof + 1 : dofP) ) > 1E-9 ) then
               flag = .false.
            endif

   !         if (i==1) then
   !
   !            print*, 'zST: ', elem%zST, 'size: ', size(elem%zST)
   !            print*, 'zST_plus:' , elem%zSTplus , 'size: ',size(elem%zSTplus)
   !
   !         end if

            deallocate( Dw )
          enddo !i

       endif

       !deallocate WS
       do i = 1, nelem
         elem => grid%elem(i)
         deallocate( elem%wS )
       enddo !i

    endif ! R_type

    if ( flag ) then
      print*, 'primal:', primal, 'R_type(L2/H1/H-1:', R_type
      stop 'Problem in reconstructSolution - the reconstructed solution is almost (1E-9) zero at ALL triangles.'
    endif

  end subroutine reconstructSolution

  !> computes the nonlinear algebraic estimate for DWR method
  !> output: DWR%estimNL = res(u_h)(z_h) = F(u_h) * z_h
  subroutine computeNonlinDWRestimates( DWR, grid)
      type( DWR_t ), intent(inout) :: DWR
      class( mesh ), intent(inout) :: grid
      class( element ), pointer :: elem
      integer :: i
      real, dimension(1:ndim) :: loc_eta, alg_estim
      logical :: loc_implicitly

      loc_implicitly = state%nlSolver%implicitly
      state%nlSolver%implicitly = .false.
      ! fill elem%vec(rhs,:)  - residual
      call ComputeSTDGM_Terms( .false. )
      state%nlSolver%implicitly = loc_implicitly

      alg_estim(1:ndim) = 0.0
      do i=1, grid%nelem
         elem => grid%elem(i)
         !compute the algebraic error
         loc_eta(1:ndim) = DWRElemEstim_alg(elem)
         alg_estim(1:ndim) = alg_estim(1:ndim) + abs( loc_eta(1:ndim) )
      !         elem%eta(dwrA,1:ndim)  = abs( loc_eta(1:ndim) )
       enddo

       DWR%estimNL = abs( alg_estim(1) )

       if ( state%nlSolver%non_alg_stop == 'aDWR' ) then
         ! C_L should not be here (safety parameter is used directly when computing)
         DWR%aDWR%linTol = DWR%estimNL * DWR%aDWR%C_A
       endif

  end subroutine computeNonlinDWRestimates



   !> compute the LS reconstruction from primal and dual solution
  !> compute the weighted residual estimates:
  !> 1. global estimate ~ J(u) - J(u_h)
  !> 2. local indicators for the mesh refinement
  !> discretization error estimates
  ! STATIONARY PROBLEMS ONLY
  subroutine computeDWRestimates( DWR, grid )
    type( DWR_t ), intent(inout) :: DWR
    class( mesh ), intent(inout) :: grid
    class( element ), pointer :: elem
    integer :: i, j, kk, k, elemDof, dof, Qnum, nelem, mdof
    integer :: iPlot
    character(len=50) :: plot_file

    call state%cpuTime%startEstimTime()

    nelem = grid%nelem
    !write(debug,*) 'alokace zST : nezkracuji ho'

    ! Put the whole block to DWREstimate subroutine

    ! compute the Dual estimator
    ! depending on the dual residual r*(u_h,z_h)( I_h^{p+1 }u_h  - u_h )
    ! multiply b = ( J(phi) - zST^T * C^+ ) )
    ! multiply b * (u_plus - u)
    call DualDWRrezidErrorEstimates( grid, DWR, .true. )
    ! update back to the smaller system - done in DualDWRRezidErrorEstimates

    ! compute projection OR reconstruction of the solution zSTplus
    if ( DWR%deg_plus ) then

      !  change the size of zST to elem%dof!
      do i = 1,grid%nelem
         elem => grid%elem(i)
         elem%zSTplus => elem%zST
         nullify( elem%zST )
         ! projection of z_{p+1} to z_p
         allocate( elem%zST(1:ndim, 1:elem%dof, 1:elem%Tdof) , source = 0.0 )
         ! due to the OG of basis function this is the projection prom P^{p+1} to P^p
         !FERROR : WOULD NOT WORK FOR CURVED ELEMENTS
         elem%zST(1:ndim, 1:elem%dof, 1:elem%Tdof) = &
               elem%zSTplus(1:ndim, 1:elem%dof, 1:elem%Tdof)
      end do !i
    else
       ! reconstruct the solutions I_h(z)
       call reconstructSolution( grid, .false. , DWR%R_type )
    endif


    ! we compute the standard RES indicators and use them to mesh adaptation
    ! TRUE -> use standard RES error indicators for mesh adaptation
    if (DWR%RESindicators) then
       !print*, 'ttime:' , state%time%ttime
       !write(*,'(a10, 16es12.4)') 'ESTIMSS1:',  sqrt(state%estim(1:5,1)), grid%elem(1)%estim_loc, &
       !      grid%elem(1)%eta(2, 1)

       ! saves the indicators to elem(:)%estim_loc
!       call RezidErrorEstimates( .true. , .true. )
       call RezidErrorEstimates( .false. , .true. )

       !write(*,'(a10, 16es12.4)') 'ESTIMSS2:',  sqrt(state%estim(1:5,1)), grid%elem(1)%estim_loc, &
       !      grid%elem(1)%eta(2, 1)

!         print*, 'ESTIMATES:' , state%estim( :, : )
       state%estim( dwrS, 1:ndim ) = state%estim( max_resT_S , 1:ndim )
       !write(*,'(a10, 16es12.4)') 'ESTIMSS3:',  sqrt(state%estim(1:5,1)), grid%elem(1)%estim_loc, &
       !      grid%elem(1)%eta(2, 1)

       ! the total RES estimate
       state%estim( dwr_etaS, 1:ndim ) = state%L_estim(resS)**2 /  state%time%tau(1)

       ! fill the elem%eta(resST, 1) array
       do i = 1, nelem
          ! copy the estimate to eta_S
          grid%elem(i)%eta(dwrS, DWR%J%coord_i ) = grid%elem(i)%estim_loc
          !            print*, 'estim_loc:', grid%elem(i)%estim_loc
       enddo

    ! standard way of DWR
    else
       ! primal rezidual DWR estimate
       ! reconstruct or project the solution u ~ deg_plus
       ! computed from the PRIMAL residual
       call DWRrezidErrorEstimates( grid, DWR, .true. )

       do i = 1, nelem
            elem => grid%elem(i)

            ! FR PROBLEM: here we set eta_aver = |eta| + |eta*|
            ! due to the theory , it may be better to set eta_aver = |eta + eta*|
            elem%eta( dwr_aver, DWR%J%coord_i) = &
               0.5 * ( elem%eta(dwr_sign, DWR%J%coord_i) * elem%eta(dwrS, DWR%J%coord_i) &
                     + elem%eta(dwr_dual_sign, DWR%J%coord_i) * elem%eta(dwr_dualS, DWR%J%coord_i) )
!            elem%eta( dwr_aver, DWR%J%coord_i) = &
!               0.5 * ( elem%eta(dwrS, DWR%J%coord_i) + elem%eta(dwr_dualS, DWR%J%coord_i) )
!            print*, 'rel difference between primal and dual eta:', &
!               elem%eta(dwrS, DWR%J%coord_i) , elem%eta(dwr_dualS, DWR%J%coord_i), &
!               (elem%eta(dwrS, DWR%J%coord_i) - elem%eta(dwr_dualS, DWR%J%coord_i)) / elem%eta(dwrS, DWR%J%coord_i)

            elem%estim_loc = elem%eta( DWR%eta_index , DWR%J%coord_i )
       enddo

    endif

    call DWR%J%computeJu( grid )
    state%estim( dwr_Juh, 1:ndim ) =  DWR%J%Ju**2.0 ! has to be sqrt in future - its done in WriteOutputError
    call DWR%J%computeJu_exact( grid )

    state%estim( dwrE ,1:ndim ) = 0.0
    state%estim( dwr_aver,1:ndim ) = ( 0.5*( sqrt( state%estim(dwrS,1:ndim) ) &
               + sqrt( state%estim(dwr_dualS,1:ndim) ) ) )**2.0
    state%estim( dwr_aver_abs,1:ndim ) = ( 0.5*( sqrt( state%estim(dwrS_abs,1:ndim) ) &
               + sqrt( state%estim(dwr_dualS_abs,1:ndim) ) ) )**2.0

    ! Ju_exact is not right if the exact solution is unknown
    state%estim( dwrE, DWR%J%coord_i ) = abs( DWR%J%Ju - DWR%J%Ju_exact )**2.0 ! its squarerooted in WriteOutputEstims

    ! stopping criterion in aDWR
    DWR%estimDiscr = 0.5*( sqrt( state%estim(dwrS,1) ) &
               + sqrt( state%estim(dwr_dualS,1) ) )
    !print*, 'estim DISCR: ' , DWR%estimDiscr

    call WriteDWRErrorsScreen( DWR )

    ! write the errors into the file "aDWR_errors" (DWR%aDWR%file_error_name)
    ! this file is used to write tables and graphs of the DWR errors
    if (state%nlSolver%non_alg_stop == 'aDWR') then
      call WriteAdwrErrorsScreen(DWR)
      call DWR%aDWR%writeErrorFile( )
      ! C_L should NOT be in this tolerance
      DWR%aDWR%nlTol = DWR%estimDiscr * DWR%aDWR%C_A ! * DWR%aDWR%C_L
    endif

    iPlot = 93
    plot_file = 'plot_eta.gnu' !'../DWR/plot_eta.gnu'
    ! PLOT estimate
    ! splot 'plot_sol.gnu' with lines notitle
    open( iPlot, file = plot_file, action="write", status="replace" )
    do i = 1, grid%nelem
      elem => grid%elem(i)
!         print*, 'dwrS eta: ' , elem%eta(dwrS, 1)
      call PlotElemFunction3D( iPlot, elem,  1, elem%eta(dwrS, 1))
    end do
    print*, 'Space error indicator saved to plot_eta.gnu'
    close(iPlot)

    ! PLOT dual estimator
    plot_file = 'plot_dual_eta.gnu' !'../DWR/plot_eta.gnu'
    open( iPlot, file = plot_file, action="write", status="replace" )
    do i = 1, grid%nelem
      elem => grid%elem(i)
      call PlotElemFunction3D( iPlot, elem,  1, elem%eta(dwr_dualS, 1))
    end do
    print*, 'Space error indicator saved to plot_dual_eta.gnu'
    close(iPlot)


    if (state%nlSolver%non_alg_stop /= 'aDWR') then
       DWR%estimDiscr = sqrt(state%estim( dwr_aver,1))
       DWR%estimNL = sqrt(state%estim( dwrA, 1))
       DWR%estimLD = sqrt( state%estim( dwr_dualA, 1) )
       !DWR%estimLD = 10.0
       DWR%estimLP = 0.0
       DWR%aDWR%iter_lin_primal = state%linSolver%iter !, sqrt(state%estim( dwrS,1))
       DWR%aDWR%iter = state%nlSolver%iter !, sqrt(state%estim( dwrA, 1))
       DWR%aDWR%iter_lin_dual = DWR%linSolver_iter !, sqrt(state%estim( dwr_aver, 1))

       call DWR%writeNlErrorFile( grid%nelem)
    end if

    call state%cpuTime%addEstimTime()


!   stop 'End of ComputeDWREstimates'
  end subroutine ComputeDWRestimates

  !> DWR error computation similar to rezid error estimates
  !> the maximum is substituted by the dual solution and reconstruction
  !> perform the error estimates using the dual norm
  !> including (non-)linear algebraic error
  !> ONLY STATIONARY CASE - q=0 in STDG
  subroutine DWRrezidErrorEstimates(grid, DWR, onlyAS )
    class( mesh ), intent(inout) :: grid
    type( DWR_t ), intent(inout) :: DWR
    logical, intent(in) :: onlyAS   ! only space and algebraic estimates
    class(element), pointer :: elem, elem1
    real, dimension(1:ndim) :: alg_estim, space_estim, space_estimABS, loc_eta
    integer :: i, iPlot
    character(len=50) :: plot_file
    logical :: loc_implicitly
    real :: t1, t2, ttime
    integer :: kk = 0
    real, allocatable, dimension(:) :: one
    logical :: deg_plus

    deg_plus = .true.

    allocate( one(1:ndim), source = 1.0 )

    if ( .not. onlyAS ) &
      stop 'DWRrezidErrorEstimates not implemented for onlyAS = .false.'

    state%space%adapt%stop_adaptation=0

    ttime = state%time%ttime

    loc_implicitly = state%nlSolver%implicitly

    state%nlSolver%implicitly = .false.
    grid%elem(:)%deg_plus = .true.

    ! setting of fields elem%vec(rhs,*), elem%vec(rhsT,*) for error estimate
    ! deg_plus = .true.
    call ComputeSTDGM_Terms( .true. )

    call cpu_time(t1)

    if( state%time%disc_time /= 'STDG') &
          stop 'DWRrezidErrorEstimates only for STDGM'

    if ( state%time_dependent ) then
      stop 'DWR method is not implemented for time-dependent problems YET!'

    !STATIONARY
    else
      alg_estim(1:ndim) = 0.0
      space_estim(1:ndim) = 0.0 ! sum
      space_estimABS(1:ndim) = 0.0 ! sqrt of **2.

      if (state%time%deg > 0 ) &
         stop 'DWRElemEstim are implemented only for q == 0!'

      if (DWR%PU) then
         !print*, 'DWR%PU we need p+2 long residual, dont forget to put it back'
         !print*, 'Omezit volani ComputeSTDGM_Terms'
         call updateProblemPlus( grid, .true. )
         ! we need the larger matrix
         state%nlSolver%implicitly = .false.
         !grid%elem(:)%deg_plus = .true.
         call ComputeSTDGM_Terms(deg_plus )
         call updateProblemPlus( grid, .false. )
         !grid%elem(:)%deg_plus = .true.
         call ComputeSTDGM_Terms(deg_plus )
      endif

      ! compute the vertex oriented estimates using PU
      if (DWR%PU) then

         allocate( DWR%pu_estim(1: grid%nelem) )

         do i = 1,grid%nelem
            DWR%pu_estim(i) = Elemental1_t(3)
            elem => grid%elem(i)
            DWR%pu_estim(i)%x(1:3) = DWRElemEstim_PU(elem)
         enddo
         call DWR%distributePuEstimToVertices( grid )

         iPlot = 56
         plot_file = "vertex_estims.gnu"
         open( iPlot, file = plot_file, action="write", status="replace" )
         call grid%plotVertexFunction3D( iPlot, DWR%vertex_estim(1:grid%npoin) )
         print*, 'Vertex estims saved to vertex_estims.gnu'
         close(iPlot)
      endif


      do i=1, grid%nelem
         elem => grid%elem(i)

         ! we use elem%vec( rhs, :) where the vector a(u,\phi_i) from 1 to ndof_plus is saved

         ! compute the space error
         if (DWR%PU) then
            loc_eta(1:ndim) = DWR%EtaFromVertices( elem ) ! distribute estims from vertices to elements
         else
            loc_eta(1:ndim) = DWRElemEstim_space(elem)
         endif

         space_estim(1:ndim) = space_estim(1:ndim) + loc_eta(1:ndim)
         space_estimABS(1:ndim) = space_estimABS(1:ndim) + abs( loc_eta(1:ndim) )
         elem%eta(dwrS,1:ndim) = abs( loc_eta(1:ndim) )
         ! signum of eta( dwrS ), can be used when we try to improve the localization
         elem%eta(dwr_sign,1:ndim) = sign( one(1:ndim), loc_eta(1:ndim) )

         !compute the algebraic error
         loc_eta(1:ndim) = DWRElemEstim_alg(elem)
         !          alg_estim(1:ndim) = alg_estim(1:ndim) + loc_eta(1:ndim)
         alg_estim(1:ndim) = alg_estim(1:ndim) + abs( loc_eta(1:ndim) )
         elem%eta(dwrA,1:ndim)  = abs( loc_eta(1:ndim) )
         !          if ( minval(loc_eta(1:ndim)) < 0.0 ) then
         !            print*, '(PROBLEM?) The local algebraic error indicator is negative on element(', &
         !               elem%i , ') and equals:', loc_eta(1:ndim)
         !          endif
       enddo

       state%estim( dwrS, 1:ndim ) = space_estim(1:ndim)**2.0 ! has to be sqrt in future - its done in WriteOutputError
       state%estim( dwrS_abs, 1:ndim) = space_estimABS(1:ndim)**2.0 ! has to be sqrt in future - its done in WriteOutputError
       state%estim( dwrA ,1:ndim ) = alg_estim(1:ndim)**2.0 ! has to be sqrt in future - its done in WriteOutputError

!       print*, 'alg_estim old = ' , alg_estim(1:ndim)

    endif !stationary

    state%time%ttime = ttime
    state%nlSolver%implicitly  = loc_implicitly
    grid%elem(:)%deg_plus = .false.

    call cpu_time(t2)
    state%CPU_constaint = state%CPU_constaint + t2 - t1

    deallocate(one)

  end subroutine DWRrezidErrorEstimates


  !> including (non-)linear algebraic error
  !> ONLY STATIONARY CASE - q=0 in STDG
  !> compute the DUAL estimator depending on the dual residual r*(u_h,z_h)( I_h^{p+1 }u_h  - u_h )
  subroutine DualDWRrezidErrorEstimates( grid, DWR, onlyAS )
    class( mesh ), intent(inout) :: grid
    type( DWR_t ), intent(inout) :: DWR
    logical, intent(in) :: onlyAS   ! only space and algebraic estimates
    class(element), pointer :: elem, elem1
    integer :: i, iPlot
    character(len=50) :: plot_file
    logical :: loc_implicitly
    real :: t1, t2, ttime
    real, allocatable, dimension(:) :: x,b
    integer :: ivec, kvec, ndof, dof, dof_plus, mdof, plusDeg, Tdof, k, l, nsize, dd
    real, dimension(1:ndim) :: loc_eta, loc_etaA
    real, dimension(1:ndim) :: space_estim, space_estimABS, alg_estim
    integer :: j, dofA
    real, allocatable, dimension(:) :: wi
    real, pointer, dimension(:,:,:) :: temp
    logical :: impl
    integer :: kk = 0
    real, allocatable, dimension(:) :: one

    allocate( one(1:ndim) , source = 1.0 )

    if ( .not. onlyAS ) &
      stop 'DualDWRrezidErrorEstimates not implemented for onlyAS = .false.'

    if (state%time%deg > 0 ) &
      stop 'DualDWRElemEstim are implemented only for q == 0!'

!    state%space%adapt%stop_adaptation = 0
    ! for ComputeTerms
    impl = state%nlSolver%implicitly
    state%nlSolver%implicitly = .true.

!    call cpu_time(t1)

    if( state%time%disc_time /= 'STDG') &
          stop 'DWRrezidErrorEstimates only for STDGM'

    if ( state%time_dependent ) then
      stop 'DWR method is not implemented for time-dependent problems YET!'

    !STATIONARY
    else
      plusDeg = state%space%plusDeg

      ! we have to increase the pol degrees, for deg_plus it is already done!
      if ( .not. DWR%deg_plus ) then
         call updateProblemPlus( grid, .true. )
         call ComputeSTDGM_Terms( .false. )
      endif

      ! we need to fill elem%blockPlusNow - Ritz reconstruction
      ! fill blockPlus

      if (state%nlSolver%non_alg_stop == 'aDWR' .and. DWR%PU) &
            stop 'aDWR is not implemented with PU!'
      if ( state%nlSolver%non_alg_stop /= 'aDWR') &
            call CopyBlocksSTtoBlockPlus()

      ! add zeros to the end of zST -> size ~ dof_plus

!      nsize = sum( grid%elem(:)%dof_plus * ndim )
      nsize = sum( grid%elem(:)%dof * ndim * grid%elem(:)%Tdof )

      allocate( x(1:nsize), source = 0.0 )
      allocate( b(1:nsize), source = 0.0 )


      ! fill zST to x
      ! ALGEB - use DWR%x -> x, just enlarge the size
      ivec = 0
      do i=1,grid%nelem
         elem => grid%elem(i)
         !dof_plus = elem%dof_plus
         ! lower dof ~ p

         mdof = DOFtriang( elem%deg - plusDeg )
         ! ~ p+1
         dof = elem%dof
         Tdof = elem%Tdof
         ndof = dof * Tdof * ndim
         kvec = ivec

         do l = 1, Tdof
            do k = 1, ndim
                ! SHOULD BE ALWAYS ZST ACCORDING TO THEORY
!                if (DWR%deg_plus) then
                  x(kvec+1:kvec + mdof) =  elem%zST(k,1:mdof, l)
                  x(kvec+1+mdof:kvec+dof) = 0.0
!                else
!                  x(kvec+1:kvec + mdof) =  elem%zST(k,1:mdof, l)
!                  x(kvec+1+mdof:kvec+dof) = 0.0
!!                  x(kvec+1:kvec + dof) =  elem%zSTplus(k,1:dof, l)
!                endif
               kvec = kvec + dof
            enddo !k
         enddo !l

         ivec = ivec + ndof
      end do !i

      ! multiply zST^T * C -> b
      ! p has to be increased , otherwise we multiply only with a part from the matrix
      call bMVprodST_Dual( b, x, nsize )
!      print*, 'x :' norm2(x)
!      print*, 'b :', norm2(b)
      ! control b - the first columns should be almost zero


      ! DECREASE THE POLYNOMIAL DEGREE NOW
!      print*, 'DECREASE THE POLYNOMIAL DEGREE NOW '
!      print*, 'WATCH: We may not allocate everything, since mesh adaptation is coming now.'
!      UPDATE IN FUTURE: We may not allocate everything, since mesh adaptation is coming now.
      call updateProblemPlus( grid, .false. )

      ! instead of wSS , use DWR%dualRes(:)%x(:,:,:) - used to Ritz reconstruction
      call DWR%fillDualResVector( grid, nsize, b(1:nsize) )
      deallocate( x, b )

      ! smaller matrix again
      ! we need to compute ComputeTerms only for deg_plus == false

      ! Not needed - matrix is already computed
      !call ComputeSTDGM_Terms( .false. )

      ! reconstruct the primal solution wST -> wSTplus
      call reconstructSolution( grid, .true. , DWR%R_type )
      ! it should be already computed for deg_plus
!      print*, 'control rhsST - must be always called after ComputeTerms! (for the right dof)'

      grid%elem(:)%deg_plus = .true.

      alg_estim(1:ndim) = 0.0
      space_estim(1:ndim) = 0.0 ! sum
      space_estimABS(1:ndim) = 0.0 ! sqrt of **2.

      do i = 1, grid%nelem
         elem => grid%elem(i)

         dof = elem%dof
         dofA = elem%dof_plus
         Tdof = elem%Tdof
         ! put to new subroutine on every elem
         allocate( wi(1:dofA) )

         do j = 1,ndim
            ! DWR_PARAMETERS
            ! CONTROL whether (wSTplus - wST) or wSTplus is used, differs by the alg error
            wi(1:dof) = elem%wSTplus( j, 1:dof,1 ) - elem%wST(j,1:dof,1)
            wi(dof+1:dofA) = elem%wSTplus(j,dof+1:dofA,1)

            loc_eta(j) = dot_product( wi(1:dofA), DWR%dualRes(i)%x(j,1:dofA,1)  )
            loc_etaA(j) = dot_product( elem%wST(j,1:dof,1), DWR%dualRes(i)%x(j,1:dof,1) )
            !print*, 'dual eta:', loc_eta, loc_etaA

         enddo ! j

!         deallocate( elem%wSS )
         deallocate( wi )

         space_estim(1:ndim) = space_estim(1:ndim) + loc_eta(1:ndim)
         space_estimABS(1:ndim) = space_estimABS(1:ndim) + abs( loc_eta(1:ndim) )
         elem%eta( dwr_dualS, 1:ndim ) = abs( loc_eta(1:ndim) )
         ! the signum of eta_dualS - may be used when we try to improve the localization
         elem%eta(dwr_dual_sign,1:ndim) =  sign( one(1:ndim) , loc_eta(1:ndim) )

         alg_estim(1:ndim) = alg_estim(1:ndim) + abs( loc_etaA(1:ndim) )
         elem%eta(dwr_dualA, 1:ndim)  = abs( loc_etaA(1:ndim) )

      end do !i

      state%estim( dwr_dualS, 1:ndim ) = space_estim(1:ndim)**2.0 ! has to be sqrt in future - its done in WriteOutputError
      state%estim( dwr_dualS_abs, 1:ndim) = space_estimABS(1:ndim)**2.0 ! has to be sqrt in future - its done in WriteOutputError
      state%estim( dwr_dualA ,1:ndim ) = alg_estim(1:ndim)**2.0 ! has to be sqrt in future - its done in WriteOutputError

      grid%elem(:)%deg_plus = .false.

    endif !stationary

    !put back
    state%nlSolver%implicitly = impl

!    call cpu_time(t2)
!    state%CPU_constaint = state%CPU_constaint + t2 - t1
     deallocate( one )

  end subroutine DualDWRrezidErrorEstimates


  subroutine Set_Elem_Regul_Estim_Decay(elem)
   type(element), target, intent(inout):: elem ! elem = element
    real, dimension(:), pointer :: weights
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:,:), pointer :: Dphi
    real, dimension(:,:), allocatable :: wi, wii
    real, dimension(:,:,:), allocatable :: wwP   ! projections
    integer :: i, Qnum, Qdof, dof, dofL, deg, k, kst, l, ifile, ifile1, ifile2
    integer :: min_deg, ni
    real, dimension (:, :), allocatable :: ss
    real, dimension (:, :), allocatable :: eta_loc
    real :: val, val1, order, sumy, sumyi, sumx, sumx2, sumy2, sumxy, eta_resS
    logical :: iprint, singularity, loc_implicitly

    iprint = .false.
    !iprint = .true.

    call Detect_apriori_known_singularity(elem, singularity)
    if(singularity) iprint = .true.

    if( mod(elem%i, 1) == 0) iprint = .true.

    !if(elem%i == 1 .or. elem%i == 101 .or. elem%i == 111) iprint = .true.   ! mesh LL-shape.uns.365.grid
    !if(elem%i == 2 .or. elem%i == 38 ) iprint = .true.  ! mesh LL-shape.uns.100.grid


    ! storing of array eta
    allocate(eta_loc(1:max_eta, 1:ndim) )
    eta_loc(1:max_eta, 1:ndim) = elem%eta(1:max_eta, 1:ndim)

    dof = elem%dof


    ! projection of the solution of the element
    allocate(wwP(0:elem%deg, 1:dof, 1:ndim) )
    !if(iprint) &
    !   write(*,'(a6,3i5, 50es12.4)') 'proj W:',elem%i, elem%deg, deg, elem%w(0, 1:dof)
    do k=1, ndim
       wwP(elem%deg, 1:dof, k) = elem%w(0, (k-1)*dof + 1 : k*dof)
    enddo


    !write(*,'(a6,3i5, 50es12.4)') 'proj W:',elem%i, elem%deg, deg, wwP(elem%deg, 1:dof, 1)
    do deg = elem%deg, 0, - 1
       call Energy_Elem_deg_projection(elem, deg, wwP(deg, 1:dof, 1:ndim) )
    !!   if(iprint) then
    !!   write(*,'(a6,3i5, 50es12.4)') 'proj W:',elem%i, elem%deg, deg, wwP(deg, 1:dof, 1)
    !!   endif
    enddo

    !write(*,'(a12, 3i5, l3,40es12.4)') 'rhs:',elem%i, elem%deg, &
    !     elem%dof, state%nlSolver%implicitly, elem%vec(rhs,:)

    loc_implicitly = state%nlSolver%implicitly
    state%nlSolver%implicitly = .false.

    elem%deg_plus = .true.


    ! files outputs
    if(singularity) then
       ifile = state%space%adapt%adapt_level + 550
    else
       ifile = state%space%adapt%adapt_level + 500
    endif
    ifile1 = ifile + 100
    ifile2 = ifile + 200

    ! least squares
    sumx = 0.; sumy = 0.; sumx2 = 0.; sumy2 = 0.; sumxy = 0.; ni = 0

    min_deg = max(0, elem%deg -2)
    !min_deg = 0.
    do deg = elem%deg, min_deg, - 1
       do k=1, ndim
          elem%w(0, (k-1)*dof + 1 : k*dof)  = wwP(deg, 1:dof, k)
       enddo
       call Compute_ONLY_ONE_ELEMENT_Terms(elem )

       eta_resS = elem%eta(resS, 1)

       if( state%time%disc_time /= 'STDG') then
          call DualElemEstimate(elem, 3, .true.)  ! 3 => element residuum in X- norm
       else
          !call ST_DualElemEstimate_Var2(elem,3)
          !call ST_DualElemEstimate(elem, 2 )  ! 2 => element residuum in the H^1-norm (for inviscid??)
          call ST_DualElemEstimate(elem, 3 )  ! 3 => element residuum in X- norm
       endif

       !if(elem%deg == 2) &
            write(*,'(a12, 3i5, l3,2es12.4,a2, 40es12.4)') 'rhs:',elem%i, deg, (deg+1)*(deg+2)/2, &
            state%nlSolver%implicitly, elem%eta(resA:resS, 1), '|',elem%vec(rhs,:)

       ! least squares:
       if(deg <= elem%deg .and. deg >= 2 ) then
          ni = ni + 1
          sumx = sumx + log(1. * deg)
          sumx2 = sumx2 + log(1. * deg)**2
          sumy = sumy + log(elem%eta(resS, 1))
          sumy2 = sumy2 + log(elem%eta(resS, 1))**2
          sumxy = sumxy + log(1. * deg) * log(elem%eta(resS, 1))
         ! write(*,'(a8, i5, 30es12.4)') 'L R 23:',ni, 1.*deg, elem%eta(resS, 1), sumx, sumx2, sumy, sumy2, sumxy
       endif

       ! direct evaluation
       val = 1.
       val1 = 1.
       if(deg < elem%deg .and. deg >= 2 ) then
          val = 1.5 + log(1./sqrt(elem%area) *  eta_resS / elem%eta(resS, 1)) / log(1.*(deg -1 ) / deg)
          val1 = 1.5 + log( eta_resS / elem%eta(resS, 1) ) / log(1.*(deg -1 ) / deg)
          !print*,'#E@S@', eta_resS / elem%eta(resS, 1) , 1.*(deg -1 ) / deg,  &
          !     log( eta_resS / elem%eta(resS, 1)), log(1.*(deg -1 ) / deg)


          write(ifile1, *) elem%i, elem%deg, elem%deg - 0.1 * (elem%deg - deg), &
               elem%eta(resS, 1),val, val1

          !print*,'vals:', val, val1

       endif

       write(ifile, *) elem%i, elem%deg, elem%deg - 0.1 * (elem%deg - deg), &
               elem%eta(resS, 1),val, val1
    enddo

    ! least squares
    if( ni * sumx2-sumx*sumx /= 0.) then
       order =  - (ni*sumxy -sumx*sumy)/(ni*sumx2-sumx*sumx) + 1.5
    else
       order = 1.
    endif

    write(ifile2, *) elem%i, elem%deg, elem%deg, order, elem%xc

    !if(elem%deg == 2) print*,'---', order
    !if(ni == 3) stop "93u30ojd3ow,xwsaxzwsq]'xs21"

    write(ifile, '(x)')
    write(ifile1, '(x)')

    ! RESTORING varaibles and array back
    state%nlSolver%implicitly = loc_implicitly
    do k=1, ndim
       elem%w(0, (k-1)*dof + 1 : k*dof)  = wwP(elem%deg, 1:dof, k)
    enddo

    elem%eta(1:max_eta, 1:ndim) = eta_loc(1:max_eta, 1:ndim)
    deallocate(eta_loc)

  end subroutine Set_Elem_Regul_Estim_Decay


  !> compute the Ritz reconstruction based on the local problem
  !> \f$ a(uPlus, \varphi) = res(\varphi) \forall \varphi \in S_{hp}\f$
  !> primal / dual solution is used
  !> if return flag == .true. then the reconstruction is zero on whole domain
  !> it needs to the matrix elem%blockPlus,and RHS saved in elem%rhsST for primal and DWR%rhs for dual
  subroutine RitzReconstruction( grid , primal, flag, DWR)
    class( mesh ), intent(inout) :: grid
    logical, intent(in) :: primal
    logical, intent(out) :: flag
    type( DWR_t ), intent(in), optional :: DWR
    class(element), pointer :: elem
    integer :: nelem, i, j , k, nsize, dof, Tdof
    real, dimension(:), allocatable :: b ! rhs
    logical :: loc_implicitly
    logical :: deg_plus

    deg_plus = .true.

    ! we suppose that DWR%deg_plus = false ! problem for DWR_P, where it is TRUE

    nelem = grid%nelem

    flag = .true.

    print*, 'Ritz reconstruction primal:', primal

    if (.not. allocated( grid%elem(1)%blockPlus ) ) &
      stop 'blockPlus is not allocated in RitzReconstruction'


    if (primal) then

      ! fill the RHS
      loc_implicitly = state%nlSolver%implicitly
      state%nlSolver%implicitly = .false.
      !grid%elem(:)%deg_plus = .true.
      call ComputeSTDGM_Terms( deg_plus )
      !grid%elem(:)%deg_plus = .false.
      state%nlSolver%implicitly = loc_implicitly

      ! not done yet
      do i = 1, nelem
         elem => grid%elem(i)
         dof = elem%dof_plus
         Tdof = elem%Tdof
         nsize = size( grid%elem(i)%blockPlus%Mb(:,1) )
         allocate( b(1:nsize), source = 0.0 )

         b(1:nsize) = copy3Darrayto1Darray( elem%rhsST( 1:ndim,1:dof,1:Tdof), nsize )
!         if (i==1) then
!!            print*, 'MATRIX::' , elem%blockPlus%Mb(:,:)
!         if ( abs( sum(elem%rhsST( 1:ndim,1:dof,1:Tdof) - DWR%dualRes(i)%x(:,:,:))) > 0.00001 ) &
!            print*, 'difference in primal and dual RHS::', elem%rhsST( 1:ndim,4:dof,1:Tdof) - DWR%dualRes(i)%x(:,4:dof,:)
!            print*, 'elem vec:' , rhs, elem%vec(rhs,:)
!            stop
!         endif
!
         call SolveLocalMatrixProblem(nsize, elem%blockPlus%Mb(1:nsize,1:nsize), 1, b(1:nsize) )

         !rhs -> wSTplus
         if ( associated( elem%wSTplus ) ) &
               deallocate( elem%wSTplus )
         allocate( elem%wSTplus( 1:ndim, 1:dof, 1:Tdof ) )
         elem%wSTplus(1:ndim, 1:dof, 1:Tdof) = copy1DarrayTo3Darray( &
                                       b(1:nsize), ndim, dof, Tdof)
         !CONTROL - the reconstruction is not zero somewhere
         if ( norm2(b(:)) > 1E-9 ) &
               flag = .false.
         deallocate( b )
      enddo

    else ! dual !will be never called for DWR_P
      if ( .not. present(DWR) ) &
         stop 'RitzReconstraction fith primal=.false. cannot be called without DWR, needed for the rhs!'

      do i = 1, nelem
         elem => grid%elem(i)
         dof = elem%dof_plus
         Tdof = elem%Tdof
         nsize = size( grid%elem(i)%blockPlus%Mb(:,1) )
         allocate( b(1:nsize), source = 0.0 )

         b(1:nsize) = DWR%dualRes(i)%copyTo1DArray( nsize )
         call SolveLocalTrasposedMatrixProblem(nsize, elem%blockPlus%Mb(1:nsize,1:nsize), 1, b(1:nsize) )

         !rhs -> zSTplus
         if ( associated( elem%zSTplus ) ) &
               deallocate( elem%zSTplus )
         allocate( elem%zSTplus( 1:ndim, 1:dof, 1:Tdof ) )
         elem%zSTplus(1:ndim, 1:dof, 1:Tdof) = copy1DarrayTo3Darray( &
                                       b(1:nsize), ndim, dof, Tdof)
         !CONTROL - the reconstruction is not zero somewhere
         if ( norm2(b(:)) > 1E-9 ) &
               flag = .false.
         deallocate( b )
      enddo
    endif

  end subroutine RitzReconstruction


end module estimates



