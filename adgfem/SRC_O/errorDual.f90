!> compute the computational error in different norms if the exact solution is known
module errorDual
  !use mesh_oper
  use main_data
  !use eval_sol
  use io_sub
  use model_oper
  !use set_solution
  use problem_oper
  use euler_problem
  implicit none

  public:: ComputeH1ErrorTimeDerDual
  public:: ComputeH1ErrorTimeDerDualMore
  public:: ElementTimeDerError
  public:: ElementRHSoscill
contains

  !> compute the dual problem \f$ -\Delta \psi = g,\ \psi=0 \mbox{ on }\partial\Omega\f$,
  !> if ityp == 1 then \f$ g = \partial_t (u(t) - u_h(t) )\f$
  !> if ityp == 2 then \f$ g = f(t^n) - f(t) \f$
  !> output  is \f$ \| \nabla \psi\|_{L^2(\Omega)}^2 \f$
  !> problem is solved in \f$ s_{hp}\f$
  subroutine ComputeH1ErrorTimeDerDual( ityp, errorDualInteg )
    integer, intent(in) :: ityp
    real, intent(inout) :: errorDualInteg
    class(element), pointer :: elem
    real, dimension(:), allocatable :: x, b
    type(Gauss_rule), pointer :: G_rule
    logical :: loc_implicitly
    real :: LAresiduum
    integer :: LAiter
    real :: normL2, normH1, errorDual
    integer :: i, j,k, l, dof, dof1, Gnum, Gdof
    integer :: NOT_CONV

    !errorDualInteg = 0.
    !return

    ! NO enrichment of the space
    do i=1,grid%nelem
       elem => grid%elem(i)
       allocate(elem%wS(0:1, 1:elem%dof*ndim) )
    enddo

    ! solution of the dual problem
    call ClearMatrixBlocks()
    call ClearVectorBlocks()

    loc_implicitly = state%nlSolver%implicitly
    state%nlSolver%implicitly = .true.

    state%RHS_presented = .false.
    state%homogenDirichlet = .true.

    ! assembling of the flux matrix with Laplace problem with u=0 on \gom, fixed in time
    do i = 1, grid%nelem
       elem => grid%elem(i)
       call ComputeOneElementTerms(Set_f_s_Laplace, Set_A_s_Laplace, Set_Ppm_Laplace,&
            Set_R_s_Laplace, Set_K_sk_Laplace, Set_S_empty, Set_DS_empty,  elem )
    enddo

    ! allocation of vectors for the right-hans side and the solution
    allocate(x(1:state%nsize), b(1:state%nsize) )

    Gnum = 2   ! 15 is the maximal one !!!!!!!!!
    G_rule => state%space%G_rule(Gnum)
    Gdof = G_rule%Qdof


    errorDualInteg = 0.
    do l=1, Gdof     ! cyclus over Gauss integ nodes in time
       call ClearVectorBlocks()

       ! actual time
       state%time%ctime = state%time%ttime - state%time%tau(1) * G_rule%lambda(l)

       ! assembling of the dual right-hand side:  g = \partial_t(u-u_h) or g = f(t) - f^n
       do i = 1, grid%nelem
          elem => grid%elem(i)

          if(ityp == 1) then
             call ElementTimeDerError(elem)

          elseif(ityp == 2) then
             call ElementRHSoscill(elem)  ! OSCILATION terms

          else
             print*,' UNKNOWN ityp in  ComputeH1ErrorTimeDerDual at errorDual.f90'
             stop
          endif
       enddo

       !eta = 0.
       call FillVector(b(1:state%nsize), 0. )

       x(1:state%nsize) = 0.


       LAiter = 0
       ! precond_update moved to linSolver
       state%linSolver%precond_update = .true.
       call SolveBlockLinearProblem(state%nsize, eta, b, x, &
            LAresiduum, LAiter, NOT_CONV)

       if(ityp == 1) print*,'dual time der,  ','LA solved', LAresiduum, LAiter
       if(ityp == 2) print*,'dual RHS oscil, ','LA solved', LAresiduum, LAiter

       !write(*,'(a10,es14.6, i5)')'LA solved', LAresiduum, LAiter

       k = 0
       do i=1,grid%nelem
          elem => grid%elem(i)
          j = elem%dof*ndim
          elem%wS(0,1:j) = x(k+1:k+j)
          k = k+j
       enddo

       errorDual = 0.
       ! evaluation of the error
       do i=1,grid%nelem
          elem => grid%elem(i)

          call ComputeH1NormElement(elem, elem%dof, normL2, normH1)
          errorDual = errorDual + normH1    ! errH1 is already squared !!!

       enddo

       errorDualInteg = errorDualInteg + errorDual * G_rule%weights(l)

    enddo ! l= 1, Gdof

    errorDualInteg = errorDualInteg*state%time%tau(1)

100 continue


    ! backward settings
    state%time%ctime =  state%time%ttime
    state%nlSolver%implicitly = loc_implicitly
    state%RHS_presented = .true.
    state%homogenDirichlet = .false.

    do i=1,grid%nelem
       elem => grid%elem(i)
       deallocate(elem%wS)
    enddo


    deallocate(x, b)


  end subroutine ComputeH1ErrorTimeDerDual



  !> compute the dual problem \f$ -\Delta \psi = g,\ \psi=0 \mbox{ on }\partial\Omega\f$,
  !> if ityp == 1 then \f$ g = \partial_t (u(t) - u_h(t) ) \f$,
  !> if ityp == 2 then \f$ g = f(t^n) - f(t) \f$,
  !> output  is \f$ \| \nabla \psi\|_{L^2(\Omega)}^2 \f$
  !> problem is solved in \f$ s_{hp}^+\f$
  subroutine ComputeH1ErrorTimeDerDualMore( ityp, errorDualInteg )
    integer, intent(in) :: ityp
    real, intent(inout) :: errorDualInteg
    class(element), pointer :: elem
    real, dimension(:), allocatable :: x, b
    type(Gauss_rule), pointer :: G_rule
    logical :: loc_implicitly
    real :: LAresiduum
    integer :: LAiter
    real :: normL2, normH1, errorDual
    integer :: i, j,k, l, dof, dof1, Gnum, Gdof
    integer :: NOT_CONV

    !errorDualInteg = 0.
    !return


    ! enrichment of the space
    do i=1,grid%nelem
       elem => grid%elem(i)

       if(elem%dof /= elem%dof_fix) print*,'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
       allocate(elem%wS(0:1, 1:elem%dof_plus*ndim) )

       elem%dof = elem%dof_plus
       elem%deg = elem%deg + 1
    enddo

    call grid%setEdgeQuadratureDegrees()

    ! new matrix shape
    state%init_only_matrix_blocks = .true.
    do i=1,grid%nelem
       elem => grid%elem(i)
       call DeleteMatrixElementShape(elem)
       ! dof_plus = .false.
       call InitMatrixElementShape(elem, .false.)
    end do

    ! assempling matrix block together
    call InitGlobalMatrixShape(grid)

    ! solution of the dual problem
    call ClearMatrixBlocks()
    call ClearVectorBlocks()


    loc_implicitly = state%nlSolver%implicitly
    state%nlSolver%implicitly = .true.

    state%RHS_presented = .false.  ! here problem ???
    state%homogenDirichlet = .true.

    ! assembling of the flux matrix with Laplace problem with u=0 on \gom, fixed in time
    do i = 1, grid%nelem
       elem => grid%elem(i)
       call ComputeOneElementTerms(Set_f_s_Laplace, Set_A_s_Laplace, Set_Ppm_Laplace,&
            Set_R_s_Laplace, Set_K_sk_Laplace, Set_S_empty, Set_DS_empty, elem )
    enddo

    ! allocation of vectors for the right-hans side and the solution
    allocate(x(1:state%nsize), b(1:state%nsize) )

    Gnum = 3   ! 15 is the maximal one !!!!!!!!!
    G_rule => state%space%G_rule(Gnum)
    Gdof = G_rule%Qdof


    errorDualInteg = 0.
    do l=1, Gdof     ! cyclus over Gauss integ nodes in time
       call ClearVectorBlocks()

       ! actual time
       state%time%ctime = state%time%ttime - state%time%tau(1) * G_rule%lambda(l)

       ! assembling of the dual right-hand side:  g = \partial_t(u-u_h) or g = f(t) - f^n
       do i = 1, grid%nelem
          elem => grid%elem(i)

          if(ityp == 1) then
             call ElementTimeDerError(elem)

          elseif(ityp == 2) then
             call ElementRHSoscill(elem)  ! OSCILATION terms

          else
             print*,' UNKNOWN ityp in  ComputeH1ErrorTimeDerDual at errorDual.f90'
             stop
          endif
       enddo

       eta = 0.
       call FillVector(b(1:state%nsize), eta )

       x(1:state%nsize) = 0.

       if(ityp == 1) print*,'dual time der'
       if(ityp == 2) print*,'dual RHS oscil'


       LAiter = 0
       ! precond_update moved to linSolver
       state%linSolver%precond_update = .true.
       call SolveBlockLinearProblem(state%nsize, eta, b, x, &
            LAresiduum, LAiter, NOT_CONV)

       !write(*,'(a10,es14.6, i5)')'LA solved', LAresiduum, LAiter

       k = 0
       do i=1,grid%nelem
          elem => grid%elem(i)
          j = elem%dof*ndim
          elem%wS(0,1:j) = x(k+1:k+j)
          k = k+j
       enddo

       errorDual = 0.
       ! evaluation of the error
       do i=1,grid%nelem
          elem => grid%elem(i)

          call ComputeH1NormElement(elem, elem%dof_plus, normL2, normH1)
          errorDual = errorDual + normH1    ! errH1 is already squared !!!

       enddo

       errorDualInteg = errorDualInteg + errorDual * G_rule%weights(l)

    enddo ! l= 1, Gdof

    errorDualInteg = errorDualInteg*state%time%tau(1)

100 continue


    ! backward settings
    state%time%ctime =  state%time%ttime
    state%nlSolver%implicitly = loc_implicitly
    state%RHS_presented = .true.
    state%homogenDirichlet = .false.

    !grid%elem(1:grid%nelem)%dof = grid%elem(1:grid%nelem)%dof_fix

    do i=1,grid%nelem
       elem => grid%elem(i)
       elem%dof = elem%dof_fix
       elem%deg = elem%deg - 1
       deallocate(elem%wS)
    enddo

    call grid%setEdgeQuadratureDegrees()

    ! new matrix shape
    do i=1,grid%nelem
       elem => grid%elem(i)
       call DeleteMatrixElementShape(elem)
       ! dof_plus = .false.
       call InitMatrixElementShape(elem, .false.)
    end do

    ! assempling matrix block together
    call InitGlobalMatrixShape(grid)
    ! setting of boundary conditions
    state%init_only_matrix_blocks = .false.

    !call ClearVectorBlocks()

    deallocate(x, b)


  end subroutine ComputeH1ErrorTimeDerDualMore





  !> compute the dual problem \f$ -\Delta \psi = g,\ \psi=0 \mbox{ on }\partial\Omega\f$,
  !> if ityp == 1 then \f$ g = \partial_t (u(t) - u_h(t) )\f$,
  !> if ityp == 2 then \f$ g = f(t^n) - f(t)\f$,
  !> output  is \f$ \| \nabla \psi\|_{L^2(\Omega)}^2 \f$
  subroutine ComputeH1ErrorTimeDerDual0( ityp, errorDualInteg )
    integer, intent(in) :: ityp
    real, intent(inout) :: errorDualInteg
    class(element), pointer :: elem
    real, dimension(:), allocatable :: x, b
    type(Gauss_rule), pointer :: G_rule
    logical :: loc_implicitly
    real :: LAresiduum
    integer :: LAiter
    real :: normL2, normH1, errorDual
    integer :: i, j,k, l, dof, dof1, Gnum, Gdof
    integer :: NOT_CONV

    errorDualInteg = 0.
    !return

    !write(*,'(a6,35es14.6)') 'iii',grid%elem(1)%w(0,:)
    !print*,'ComputeH1ErrorTimeDerDual( ) - begin'
    ! enrichment of the space
    do i=1,grid%nelem
       elem => grid%elem(i)
       dof = elem%dof
       dof1 = elem%dof_plus

       !write(300,'(a6,i5,24es12.4)' ) 'w^k  :', i, elem%w(0,:)

       elem%w(0,dof1+1:dof1+dof) = elem%w(0,1:dof)
       elem%w(1,dof1+1:dof1+dof) = elem%w(1,1:dof)
       elem%w(0:1, dof+1 :dof1) = 0.
       !WWW elem%deg = elem%deg + 1

       !elem%dof = (elem%deg + 1)*(elem%deg + 2) / 2
       elem%dof = elem%dof_plus


       ! not remove, values elem%w are necessary for the evaluation of the time derivative
       ! of the approximate solution !!!!!
       !elem%w(0:1,1:dof1) = 0.
    enddo

    call grid%setEdgeQuadratureDegrees()

    ! new matrix shape
    state%init_only_matrix_blocks = .true.
    do i=1,grid%nelem
       elem => grid%elem(i)
       call DeleteMatrixElementShape(elem)
       ! dof_plus = .false.
       call InitMatrixElementShape(elem, .false.)
    end do

    ! assempling matrix block together
    !print*,'???',state%nsize
    call InitGlobalMatrixShape(grid)
    !print*,'???',state%nsize

    ! solution of the dual problem
    call ClearMatrixBlocks()
    call ClearVectorBlocks()

    loc_implicitly = state%nlSolver%implicitly
    state%nlSolver%implicitly = .true.

    state%RHS_presented = .false.
    state%homogenDirichlet = .true.


    ! assembling of the flux matrix with Laplace problem with u=0 on \gom, fixed in time
    do i = 1, grid%nelem
       elem => grid%elem(i)
       call ComputeOneElementTerms(Set_f_s_Laplace, Set_A_s_Laplace, Set_Ppm_Laplace,&
            Set_R_s_Laplace, Set_K_sk_Laplace, Set_S_empty, Set_DS_empty, elem )
    enddo


    allocate(x(1:state%nsize), b(1:state%nsize) )

    Gnum = 2   ! 15 is the maximal one !!!!!!!!!
    G_rule => state%space%G_rule(Gnum)
    Gdof = G_rule%Qdof

    !write(*,*)' @@@@',0, state%time%ctime, state%time%ttime

    errorDualInteg = 0.
    do l=1, Gdof     ! cyclus over Gauss integ nodes in time
       call ClearVectorBlocks()

       ! assembling of the dual right-hand side:  f = \partial_t(u-u_h)

       ! actual time
       state%time%ctime = state%time%ttime - state%time%tau(1) * G_rule%lambda(l)
       !state%time%ctime = state%time%ttime  ! old version qith Gnum = 1

       do i = 1, grid%nelem
          elem => grid%elem(i)

          if(ityp == 1) then
             call ElementTimeDerError(elem)
          elseif(ityp == 2) then
             call ElementRHSoscill(elem)  ! OSCILATION terms
          else
             print*,' UNKNOWN ityp in  ComputeH1ErrorTimeDerDual at errorDual.f90'
             stop
          endif
       enddo


       eta = 0.
       call FillVector(b(1:state%nsize), eta )

       !call WriteMblock(grid%elem(1)%block(0) )
       !write(*,'(a4,40es11.3)') 'b: ', grid%elem(1)%vec(rhs,:)

       x(1:state%nsize) = 0.

       state%linSolver%tol = 1E-4

       !print*,'Solution of linear algebraic system'

       LAiter = 0

       !write(*,*) '###',size(b, 1), size(x,1)
         ! precond_update moved to linSolver
       state%linSolver%precond_update = .true.
       call SolveBlockLinearProblem(state%nsize, eta, b, x, &
            LAresiduum, LAiter, NOT_CONV)
       !write(*,'(a10,es14.6, i5)')'LA solved', LAresiduum, LAiter
       !write(*,*) 'ityp = ',ityp, ' in errorDual.f90'

       k = 0
       do i=1,grid%nelem
          elem => grid%elem(i)
          j = elem%dof*ndim
          !WWW
          elem%w(0,1:j) = x(k+1:k+j)
          k = k+j

          !write(*,'(a6,24es12.4)' ) 'w^k  :', elem%w(0,:)
       enddo

       !state%isol = state%isol + 100
       !call WriteProgressOutput( 'TS' )
       !state%isol = state%isol - 100


       errorDual = 0.
       ! evaluation of the error
       do i=1,grid%nelem
          elem => grid%elem(i)
          call ComputeH1NormElement(elem,elem%dof, normL2, normH1)
          errorDual = errorDual + normH1    ! errH1 is already squared !!!

          ! approximate solution is given back in order to compute ElementTimeDerError
          ! correctly
          dof = elem%dof_fix
          dof1 = elem%dof_plus
          !WWW
          elem%w(0,1:dof) = elem%w(0,dof1+1:dof1+dof)
       enddo

       errorDualInteg = errorDualInteg + errorDual * G_rule%weights(l)

       !write(*,'(a6,i5,8es14.6)')' errD',l, &
       !     state%time%ctime, errorDual, G_rule%weights(l),errorDualInteg

    enddo ! l= 1, Gdof

    !state%errSTnorm(L2L2eH1) = state%errSTnorm(L2L2eH1) + errorDualInteg*state%time%tau(1)
    errorDualInteg = errorDualInteg*state%time%tau(1)

    ! backward settings
    state%time%ctime =  state%time%ttime
    state%nlSolver%implicitly = loc_implicitly
    state%RHS_presented = .true.
    state%homogenDirichlet = .false.


    ! OLD
    !errorDual = errorDual * state%time%tau(1)
    !state%errSTnorm(L2L2eH1) = state%errSTnorm(L2L2eH1) + errorDual

    !print*,'@@@',state%errSTnorm(L2L2eH1) - errorDual,state%errSTnorm(L2L2eH1), errorDual
    ! backward settings
    do i=1,grid%nelem
       elem => grid%elem(i)
       !WWW elem%deg = elem%deg - 1
       elem%dof = elem%dof_fix
       dof = elem%dof
       dof1 = elem%dof_plus
       elem%w(0,1:dof) = elem%w(0,dof1+1:dof1+dof)
       elem%w(1,1:dof) = elem%w(1,dof1+1:dof1+dof)

       !write(400,'(a6,i5,24es12.4)' ) 'w^k  :', i, elem%w(0,:)

    enddo

    call grid%setEdgeQuadratureDegrees()
    ! new matrix shape
    do i=1,grid%nelem
       elem => grid%elem(i)
       call DeleteMatrixElementShape(elem)
       ! dof_plus = .false.
       call InitMatrixElementShape(elem, .false.)
    end do

    ! assempling matrix block together
    call InitGlobalMatrixShape(grid)
    ! setting of boundary conditions
    state%init_only_matrix_blocks = .false.

    deallocate(x, b)

    !write(*,'(a6,35es14.6)') '|||',grid%elem(1)%w(0,:)

    !print*,'ComputeH1ErrorTimeDerDual( ) - end,  state%errSTnorm(L2L2eH1)=',state%errSTnorm(L2L2eH1)**0.5

  end subroutine ComputeH1ErrorTimeDerDual0







  !> compute \f$ \int_K \partial_t (u_{h} - u(t))   \phi_i dx \f$, \f$ t = state%time%ctime \f$
  subroutine ElementTimeDerError(elem)
    type(element):: elem
    real, dimension(:,:), allocatable :: xi, x, ut, uht
    integer ::  dof,  Qdof, Qnum,  l

    dof = elem%dof
    Qdof = elem%Qdof
    Qnum = elem%Qnum

    allocate( ut(1:Qdof, 1:ndim) )
    allocate( xi(1:Qdof, 1:nbDim))
    allocate( x(1:Qdof, 1:nbDim))
    allocate( uht(1:Qdof, 1:nbDim))


    ! time derivative of the approximate solution in integ. nodes
    elem%dof = elem%dof_fix
    call Eval_w_t_Elem(elem, uht)
    elem%dof = dof
    !elem%dof = elem%dof_plus

    ! setting of the function \f$ f \f$ in integ nodes x(:, 1:nbDim)
    xi(1:Qdof,1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)
    call ComputeF(elem, Qdof, xi(1:Qdof,1:nbDim), x(1:Qdof, 1:nbDim) )

    ! time derivative of the exact solution in integ. nodes
    do l=1,Qdof
       call TimeDer_Exact_Scalar( x(l,1:nbDim), ut(l,1:ndim), state%time%ctime )

       !ut(l, :) = 2*(x(l,1)*(1-x(l,1)) + x(l,2)*(1-x(l,2)) )
       !ut(l, :)  = 0.
       !write(200+state%time%iter,*) x(l,1:2), ut(l,1), uht(l,1),ut(l,1)- uht(l,1)
       !write(200+int(300*state%time%ctime),*) x(l,1:2), ut(l,1), uht(l,1),ut(l,1)- uht(l,1)
    enddo

    ut(1:Qdof,1:ndim) = ut(1:Qdof,1:ndim) - uht(1:Qdof,1:ndim)

    ! evaluation of the vector
    !if(dot_product(elem%vec(rhs, :), elem%vec(rhs, :)) > 0 ) &
    !     write(*,'(a4,i5, 40es12.4)') 'eD1:',elem%i, elem%vec(rhs, :)

    call EvalVectorB(elem, ut(1:Qdof,1:ndim), dof, elem%vec(rhs,1:ndim*dof) )

    deallocate(ut, uht, x, xi)

  end subroutine ElementTimeDerError


  !> compute \f$ \int_K  (f^{n} - f(t)  \phi_i dx \f$, \f$ t = state%time%ctime \f$
  subroutine ElementRHSoscill(elem)
    type(element):: elem
    real, dimension(:,:), allocatable :: xi, x, fn, ft
    integer ::  dof,  Qdof, Qnum,  l

    dof = elem%dof
    Qdof = elem%Qdof
    Qnum = elem%Qnum

    allocate( xi(1:Qdof, 1:nbDim))
    allocate( x(1:Qdof, 1:nbDim))
    allocate( fn(1:Qdof, 1:ndim) )
    allocate( ft(1:Qdof, 1:nbDim))

    ! setting of the function \f$ f \f$ in integ nodes x(:, 1:nbDim)
    xi(1:Qdof,1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)
    call ComputeF(elem, Qdof, xi(1:Qdof,1:nbDim), x(1:Qdof, 1:nbDim) )

    ! time derivative of the exact solution in integ. nodes
    do l=1,Qdof
       call RHS_Scalar( x(l,1:nbDim), fn(l,1:ndim), state%time%ttime )
       call RHS_Scalar( x(l,1:nbDim), ft(l,1:ndim), state%time%ctime )

       !ut(l, :) = 2*(x(l,1)*(1-x(l,1)) + x(l,2)*(1-x(l,2)) )
       !ut(l, :)  = 0.
       !write(200+state%time%iter,*) x(l,1:2), ut(l,1), uht(l,1),ut(l,1)- uht(l,1)
       !write(*,'(20es14.6)') x(l,1:2), fn(l,1), ft(l,1),fn(l,1)- ft(l,1)
    enddo

    ft(1:Qdof,1:ndim) = fn(1:Qdof,1:ndim) - ft(1:Qdof,1:ndim)

    ! evaluation of the vector
    call EvalVectorB(elem, ft(1:Qdof,1:ndim), dof, elem%vec(rhs,1:ndim*dof) )

    deallocate(fn, ft, x, xi)

  end subroutine ElementRHSoscill



end module errorDual
