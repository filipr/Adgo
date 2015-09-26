!> general error estimation subroutines
module estimates
  use main_data  ! contains type(mesh) :: grid for computation
  use euler_problem
  use apost_estimation
  use project_estimation

  implicit  none

  public:: RezidErrorEstimates
  public:: JumpsEvaluation
  public:: SolutionDifference
  public:: PWpolynomialReconstComputeTerms

contains



  !> perform the error estimates using the dual norm
  !> including (non-)linear algebraic error
  subroutine RezidErrorEstimates( onlyAS, Ttime_updated )
    logical, intent(in) :: onlyAS   ! only space and algebraic estimates
    logical, intent(in) :: Ttime_updated ! Ttime was already updated by %tau(1)
    class(element), pointer :: elem, elem1
    real, dimension(:), allocatable :: L_estim
    real, dimension(:,:), allocatable :: wi
    real, dimension(:,:,:), allocatable :: Dwi
    real :: machine_tol, rmax, t1, t2, ttime, val
    integer :: i, j, k, ndof, ndofP, itest, imax, ipoc
    logical :: loc_implicitly

    !print*,'####  RezidErrorEstimates  start',onlyAS, Ttime_updated

    !itest = 360
    itest = -480

    ttime = state%time%ttime
    if(Ttime_updated) state%time%ttime = state%time%ttime - state%time%tau(1)

    allocate(L_estim(1:max_eta) )
    loc_implicitly = state%nlSolver%implicitly

    state%nlSolver%implicitly = .false.
    grid%elem(:)%deg_plus = .true.

    !itest = 25
    !write(*,'(a6,i5,16es12.4)') 'w0', itest, 0., grid%elem(itest)%w(0, 1:2), &
    !     grid%elem(itest)%w(1, 1:2),grid%elem(itest)%w(2, 1:2), grid%elem(itest)%w(3, 1:2)
    !print*,'-----------, rezid',onlyAS

    ! setting of fields elem%vec(rhs,*), elem%vec(rhsT,*) for error estimate
    call SetVectorsFields( onlyAS )

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
          call ST_DualElemEstimate(elem, 2 )  ! 2 => element residuum in the H^1-norm (for inviscid??)
          !call ST_DualElemEstimate(elem, 3 )  ! 3 => element residuum in X- norm
       endif

       ! limitation along the shock waves
       val = 1.
       if(state%type_IC == 8) val = 2*elem%area/elem%diam  ! .or. abs(elem%xc(1) - 1.) > 2E-2) then

       !if(val > 1E-2) then
       if(state%type_IC /= 8  .or. abs(elem%xc(1) - 1.) > 5E-2) then

          L_estim( resA) = L_estim( resA) + elem%estimA**2
          L_estim( resS) = L_estim( resS) + elem%estimS**2
          L_estim( resT) = L_estim( resT) + elem%estimT**2
          L_estim(resST) = L_estim(resST) + elem%estimST**2

       endif


       !if(elem%i < 5) &
       !write(*,'(a12,i5,10es12.4)') 'RES estims', &
       !     elem%i, elem%estimA, elem%estimS, elem%estimT, elem%estimST

       !if(elem%estimS /  elem%estimT  < state%estim(min_resT_S,1)) k = i

       !if(elem%estimS /  elem%estimT  < state%estim(min_resT_S,1)) k = i
       !
       ! Verfurth approach
       !print*,'$$$', elem%estimT, elem%estimS
       if(elem%estimS > 0.) &
            state%estim(min_resT_S_loc,1) = min(state%estim(min_resT_S_loc,1),  &
            elem%estimT /  elem%estimS )

       !if(elem%estimT > rmax) then
       !   rmax = elem%estimT
       !   imax = i
       !endif

       state%estim(max_resT_S,1) = max(state%estim(max_resT_S,1), elem%estimT / max(1E-15, elem%estimS) )
       !!state%estim(resA_ST,:) = 1E+50

    enddo

  ! print*, 'stopped in estimates.f90 after ST_DualElemEstimate'
  ! stop

    state%L_estim(1:max_eta) = L_estim(1:max_eta)**0.5

    !write(*,'(a6,4es12.4,a2,3es12.4)') &
    !     'EWQDRT', state%L_estim(1:4),'|',state%L_estim(1)/state%L_estim(2)

    machine_tol = 1.E-01
    state%estim(resA_ST_loc,1) = 0.
    ! local algebraic criterion
    do i=1,grid%nelem
       elem => grid%elem(i)

       if( elem%estimST >  machine_tol * state%L_estim(resST) / grid%nelem**0.5 ) then

          state%estim(resA_ST_loc,1) = max(state%estim(resA_ST_loc,1), elem%estimA/ elem%estimST)

          !write(198,*) elem%xc(:), elem%estimA, elem%estimST ,elem%estimA/ elem%estimST, elem%i
       endif
    enddo


    ! steady-state approach
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
       call SeekElemSupports(grid)  ! create the list of elements sharing at least a vertex with elem

       ! averaging of the estimate
       ! setting of 'elem%estim_loc' including elem%estimST of neighbours elements
       ! elem%estim_loc**2 = sum_{K\in N(K)} elem%estimST**2
       do i=1,grid%nelem
          elem => grid%elem(i)
          if( state%time%disc_time /= 'STDG') then
             elem%estim_loc = elem%estimS**2    ! Verfurth approach
          else
             elem%estim_loc = elem%estimST**2    ! STDGM approach
          endif

          ipoc = 1

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
                   elem%estim_loc = elem%estim_loc + elem1%estimS**2   ! Verfurth approach
                   !elem%estim_loc = max(elem%estim_loc , elem1%estimS**2)   ! Verfurth approach
                else
                   elem%estim_loc = elem%estim_loc + elem1%estimST**2   ! STDGM approach
                endif

                ipoc = ipoc + 1
             endif
          enddo
          ! !!  elem%estim_loc = elem%estim_loc**0.5  ! we store the square
          elem%estim_loc = elem%estim_loc / (ipoc*4)  !  *4 for neighbouring elements works without
       enddo

       ! deallocation of arrays allocated in SeekElemSupports
       do i=1,grid%nelem
          elem => grid%elem(i)
          deallocate(elem%supp)
       enddo

    endif

    state%time%Ttime = ttime

    !print*,'####  RezidErrorEstimates  END'

    !if( state%time%cn ) then
    !   do i=1,grid%nelem
    !      elem => grid%elem(i)
    !   enddo
    !endif
  end subroutine RezidErrorEstimates

  !> setting of field elem%vec(rhs,*), elem%vec(rhsT,*) for error estimate
  subroutine SetVectorsFields(onlyAS )
    logical, intent(in) :: onlyAS   ! only space and algebraic estimates
    integer :: ityp, ityp1, itest, i
    class(element), pointer :: elem, elem1

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
         !print*, 'SetVectorsFields calling ComputeSTDGM_Terms with implicitly = ', state%nlSolver%implicitly
       call ComputeSTDGM_Terms( )

    elseif(ityp == -1) then  ! backward Euler, piecewise linear time reconstruction
       ! PREDELAT !!!! VD

       !itest = 25
       !write(*,'(a6,i5,16es12.4)') 'w0', itest, 0., grid%elem(itest)%w(0, 1:2), &
       !     grid%elem(itest)%w(1, 1:2),grid%elem(itest)%w(2, 1:2), grid%elem(itest)%w(3, 1:2)
       !print*,'-----------, vector'

       call PWpolynomialReconstComputeTerms(onlyAS ) !flux vector:  %vec(rhs, * ) = \int_{I_m} c_h(w_h, * )

    elseif(ityp == 1) then  ! backward Euler, STDGM approach
       call ComputeTerms( )     ! flux vector:  %vec(rhs, * ) = c_h(w_h, * )

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
       !write(197,'(a6,20es12.4)') 'etas:',elem%estimA, elem%estimS, elem%estimT, elem%estimST
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
       !write(197,'(a6,20es12.4)') 'etas:',elem%estimA, elem%estimS, elem%estimT, elem%estimST
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
             call ComputeTerms( )      ! flux vector:  %vec(rhs, * ) = c_h(w_h, * )
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
    call ComputeTerms( )     ! flux vector:  %vec(rhs, * ) = c_h(w_h, * )

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

    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof  = elem%dof * ndim
       ndofP = elem%dof_plus * ndim

       allocate(elem%wS(1:1,1:elem%dof * ndim) )

       elem%wS(1,1:elem%dof * ndim)  = elem%w(0,1:ndof)

       elem%w(0,1:elem%dof * ndim) = elem%w(1, 1:ndof)
    enddo

    call ComputeTerms( )     ! flux vector:  %vec(rhs, * ) = c_h(w_h^{k-1}, * )
    do i=1,grid%nelem
       elem => grid%elem(i)
       elem%vec(rhsOLD,1:ndofP)  = elem%vec(rhs,1:ndofP)

       elem%w(0,1:ndof) = elem%wS(1,1:ndof)
       deallocate(elem%wS )
    enddo

    call ComputeTerms( )     ! flux vector:  %vec(rhs, * ) = c_h(w_h^{k }, * )

    do i=1,grid%nelem        ! Crank-Nicolson (F(w^k) + F(w^{k-1} ) /2
       elem => grid%elem(i)
       elem%vec(rhs,1:ndofP)  = (elem%vec(rhs,1:ndofP) + elem%vec(rhsOLD,1:ndofP) )/2.
    enddo
  end subroutine CrankNicolsonVectorFlux


  !>  calculation of inter-element jumps
  !> elem%rez = \f$\int_{\partial K} [u_h]^2 dS \f$,
  !> jumpL = \f$ \sum_{\Gamma \in \partial K} |\Gamma|^{-1} \int_{\Gamma} [u_h]^2 dS\f$
  subroutine JumpsEvaluation( )
    class(element), pointer :: elem
    real :: jumps, jumps2, jumpL, rez_old, jumps1, jumpL1, jumpL2
    real, dimension(:), allocatable :: jumpV, jumpV1, jumpV2
    real :: wwi
    integer :: i, j, ifile


    allocate(jumpV(1:ndim), jumpV1(1:ndim), jumpV2(1:ndim))

    ! calculating of jumps
    jumps = 0.    
    jumps1 = 0.   
    jumps2 = 0.   

    state%time%ctime = state%time%ttime

    do i=1,grid%nelem
       elem => grid%elem(i)

       elem%rezid = 0.
       elem%rezid1 = 0.

       jumpL = 0.
       jumpL1 = 0.
       jumpL2 = 0.

       ! lifting operator for the dicrete gradient in pNeu for SIPG and NIPG
       if( state%space%estim_space == 'pNeu') elem%lifting(:,:) = 0.

       do j = 1, elem%flen
          wwi = 1.
          if(elem%face(neigh,j) > 0) wwi = 0.5
          !k = elem%face(neigh,j)

          !call IntegElementEdgeJump(elem,  j, jumpV)
          call IntegElementEdgeJumpProj(elem,  j, jumpV, jumpV1, jumpV2)

          ! the whole vector
          elem%rezid = elem%rezid + sum(jumpV(:) )
          elem%rezid1 = elem%rezid1 + sum(jumpV1(:) )


          !if(elem%i < 3) print*,'###',elem%i, j, elem%rezid


          ! only the density
          !elem%rezid  = elem%rezid  + jumpV(1)
          !elem%rezid1 = elem%rezid1 + jumpV1(1)

          jumpL = jumpL + sum(jumpV(:) ) / elem%dn(j) ! penalty par 1/|\Gamma|
          jumpL1 = jumpL1 + sum(jumpV1(:) ) / elem%dn(j) ! penalty par 1/|\Gamma|

          jumpL2 = jumpL2 + wwi * sum(jumpV2(:) ) / elem%dn(j) ! penalty par 1/|\Gamma|

          !call ElementEdgeJumpIndicator(elem, j)
          !jumpL = jumpL + (elem%rezid - rez_old) / elem%dn(j) ! penalty par 1/|\Gamma|

       enddo

       jumps = jumps + jumpL  ! elem%rezid is already the square of || [u_h] ||
       !jumps1 = jumps1 + jumpL1  ! elem%rezid is already the square of || [u_h] ||
       jumps1 = jumps1 + elem%rezid  ! elem%rezid is already the square of || [u_h] ||


       jumps2 = jumps2 + jumpL2


       elem%jumpsJh = jumpL    ! for estimates || . ||_X^2 + || . ||_J^2


       !if(elem%i == 100) &
       !write(*,'(a6,2i5,12es12.4)') ' jumps:',elem%i, elem%deg, elem%rezid, elem%rezid1
       !write(94,'(7es12.4,a20)') elem%xc(1:2), elem%rezid**0.5,elem%rezid1**0.5, &
       !      elem%rezid1**0.5/ elem%rezid**0.5, &
       !      elem%rezid/ elem%area/ elem%diam**(2*max(1,elem%deg)-1)*elem%diam**4. , &
       !      (elem%rezid1**0.5/ elem%rezid**0.5) / &
       !      (elem%rezid/ elem%area/ elem%diam**(2*max(1,elem%deg)-1)*elem%diam**4. ), &
       !      'estimates.f90'

    enddo
    !state%G_jumps = jumps**0.5
    !state%G_jumps1 = jumps1**0.5
    state%err(Ejumps)  = jumps1**0.5
    state%err(EjumpsG) = jumps**0.5
    state%err(EjumpsPi) = jumps2**0.5

    deallocate(jumpV, jumpV1, jumpV2)

  end subroutine JumpsEvaluation

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

       call ComputeTerms( )     ! flux vector:  %vec(rhs, * ) = c_h(w_h^{k-1}, * )

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
    call ComputeTerms( )     ! flux vector:  %vec(rhs, * ) = c_h(w_h^{k-1}, * )

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



end module estimates



