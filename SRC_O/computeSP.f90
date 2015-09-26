!> main iterative loop, solution of the saddle point problems
!> as incompressible Navier-Stokes equations

module compute_operSP

  use main_data
  use io_sub
  use solve_problem
  use marking
  use estimates
  !use errorDual
  !use errorFlux
  use hp_adaptation
  !!use pMultiGrid !MAKE_MG
  !use helmholtz_estim
  use test

  !use dual_estim
  !use alg_estim
  use st_interpol
  use anisotropic
  use AMA_estims
  use ama_hp_interpol
  use error_subs
  use mesh_oper
  !use loc_rav_tho_ned

  implicit none


contains

!> evaluation of the state vector \f$w\f$ on element \f$ elem \f$  in
!> integ nodes, i.e. recomputation of elem%w into wi in volume integ. nodes
subroutine Eval_w_ElemSP(elem, wi)
    type(element), intent(in):: elem      ! elem = element
    real, dimension(1:elem%Qdof,1:nbdim), intent(inout):: wi ! w in integ. nodes
    real, dimension(:,:), pointer:: phi ! local store arrays
    integer :: dof, Qdof, k, kst, i

    dof = elem%dof
    Qdof = elem%Qdof

    phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)

    do i=1,Qdof
       wi(i,1) = dot_product(elem%wSP(wV1,0,1:dof), phi(1:dof,i) )
    enddo

    do i=1,Qdof
       wi(i,2) = dot_product(elem%wSP(wV2,0,1:dof), phi(1:dof,i) )
    enddo
    ! Faster ?
    !wi(1:Qdof, k) = matmul(elem%w(0,kst:kst+dof-1), phi(1:dof,1:Qdof) )

end subroutine Eval_w_ElemSP

!> evaluation of the mass RHS,
!> extrapolation from the old time levels,
subroutine SetMassVectorSP()
    class(element), pointer :: elem
    integer :: ie, i

    do ie = 1, grid%nelem
       elem => grid%elem(ie)

       call SetElementMassVectorSP(elem)

       if (state%modelName == 'incNS') then
		if (ie == 1)  print*, 'CHECK 01'
          elem%wSP(wV1,0,:) = elem%wSP(wV1,1,:)
       endif
    enddo
end subroutine SetMassVectorSP

subroutine SetElementMassVectorSP(elem)
    type(element):: elem
    real, dimension(:,:), allocatable :: wi
    integer ::  dof, dofP,  Qdof, i

    dof = elem%dof
    dofP = elem%dofP
    Qdof = elem%Qdof

    elem%wSP(wV1,0,:) = 0.
    elem%wSP(wV2,0,:) = 0.
    elem%wSP(wP,0,:) = 0.
    if (elem%i == 1) print*, 'CHECK 02'
    do i=1,state%time%deg
       elem%wSP(wV1,0,:) = elem%wSP(wV1,0,:) - state%time%alpha(i) * elem%wSP(wV1,i,:)
       elem%wSP(wV2,0,:) = elem%wSP(wV2,0,:) - state%time%alpha(i) * elem%wSP(wV2,i,:)
       elem%wSP(wP,0,:) = elem%wSP(wP,0,:) - state%time%alpha(i) * elem%wSP(wP,i,:)
       !if(elem%i == 10)
       !write(*,'(a3,2i5,12es12.4)') &
       !     ',@@',i,state%time%deg, state%time%alpha(i), elem%w(i,1) ,elem%w(0,1)
    enddo
    if (elem%i == 1) print*, 'elem%vec was not allocated'
    allocate(elem%vec(rhsM,ndim*dof))
    elem%vec(rhsM,:) = 0.

    ! setting of w in integration nodes
    allocate(wi(1:Qdof,1:ndim))
    call Eval_w_ElemSP(elem, wi(1:Qdof,1:nbdim) )
    wi(1:dof,ndim) = 0. !pressure part of the vector

    ! evaluation of the vector
    call EvalVectorB(elem, wi(1:Qdof,1:ndim), dof, elem%vec(rhsM,1:ndim*dof) )

    deallocate(wi)

  end subroutine SetElementMassVectorSP

  !> perform the one time step on the given grid
  !>
  !> evaluate the flux matrix/vector, solve by explicit, semi-implicit or implicit
  !> time discretization, for the implicit the Newton-like method is used
  subroutine PerformOneTimeStepSP( time_prepare, time_solve, time_estim )
    real, intent(inout) ::  time_prepare, time_solve, time_estim
    type(Newton_type), pointer :: Newton
    class(element), pointer :: elem
    real, dimension(:,:), allocatable :: cDLM
    real:: t1, t2, rel_res
    !real :: mu, mu1
    !real :: normL2, normH1
    integer :: iter, newton_count
    integer :: i, j, k, l, i1, i2
    logical :: precond_update, vector_update, loc_implicitly, update
    integer :: imp
    character(len=7) :: Newtonx, MatrixA, RHSideb, MatrixB, MatrixC
    real,dimension(:),allocatable :: y,z !MAKE_MG
    !integer:: deg,dof,nelem,nsize,outsize
    integer :: vypis
    !type(MGSetting) :: stateMG

    ! seting of coefficients \alpha, \beta for ABDF
    call SetABDF(state%BDF, state%time%deg_actual, state%time%tau(1:state%time%deg_actual+1) )

    print*,'Compute "Mass RHS", vector M w^{k} and extrapolate from the old levels'
    print*,'TO DONE'
    call SetMassVectorSP( )

    time_prepare = 0.
    time_solve  = 0.
    time_estim = 0.

    !!write(*,'(a6, 3es12.4)') 'times:',time_prepare, time_solve , time_estim

    Newton => state%nlSolver
    allocate(cDLM(0:Newton%max_iter, 1:5))
    Newton%norm_res = 1.

    !if(ndim == 4) call DirectComputeDrag_Lift(cDLM(0, 1:5 ) )
    !open(44, file='cDLM', status='UNKNOWN', position = 'append')
    !write(44,'(x)')
    !write(44,'(2i5,15es11.3)') state%time%iter, Newton%iter, &
    !     state%time%iter + 1.*0., cDLM(Newton%iter, 1:3)
    !close(44)

    ! minimal update of the flux matrix after imp time steps
    imp = Newton%min_update
    imp = 1
    print*,'Fixed update'


    ! for DUA error estimates
    !grid%elem(:)%CK = 0.
    !grid%elem(:)%Cb = 0.
    !grid%elem(:)%CKo = 0.
    !grid%elem(:)%Cbo = 0.

    ! Has the big vector the right size? If not (e.g., after mesh adapation) then resize
    if(size(Newton%b(:) ) /= state%nsize  )then
       deallocate( Newton%b, Newton%b1, Newton%x, Newton%rr )

       allocate(Newton%b(1:state%nsize), Newton%b1(1:state%nsize), &
            Newton%x(1:state%nsize), Newton%rr(1:state%nsize) )
    endif

    if(state%time%time_method == 'I' .and. imp > 0) then
       if( mod(state%time%iter_loc,imp) == 1 .or. imp == 1)  state%nlSolver%implicitly = .true.
    elseif(state%time%time_method == 'I' .and. imp == 0) then
       state%nlSolver%implicitly = .true.
    endif


    if(state%time%time_method == 'I') then
       open(54, file='newton-update', status='unknown', position='append')
       open(55, file='newton', status='unknown', position='append')
       !open(56, file='implicit', status='unknown', position='append')
       write(55,*) '  '
       write(55,*) '## it  i  l    ritr     ritr2   ', &
            '   |F(x)|/N  |F(x*)|/N   |F(x)|   |F(x*)|      rel_res ', &
            '     gmr_tol    lambda     theta       etaS    algeb   LAiter IM Pr'
    endif

    newton_count = 0
    Newton%iter  = 0
    do iter = 1, Newton%max_iter
       call cpu_time(t1)

       if(imp == 0)  state%nlSolver%implicitly= .true.

       loc_implicitly = state%nlSolver%implicitly
       precond_update = .false.
       update = .false.

       if( state%nlSolver%implicitly ) then
          !print*,'Compute matrix C(w), vector q(w)'
          !write(600,*)'# Compute matrix C(w), vector q(w)'
          update = .true.
          call ComputeTermsSP( )

          newton_count = 0
          Newton%updates = Newton%updates + 1
          precond_update = .true.

          if(state%time%time_method == 'I' .or. state%time%time_method == 'S') then
             if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
                !print*
                !print*,'Compute matrix C(w), vector q(w) - (2)'
                !write(600,*) '# Compute matrix C(w), vector q(w) - (2)'
                state%nlSolver%implicitly = .false.
                call ComputeTermsSP( )
                state%nlSolver%implicitly = .true.
             else
                call SetF_q_Cw_fast( )
             endif
          endif

          vector_update = .true.
          if(state%time%time_method == 'I') state%nlSolver%implicitly = .false.

       else   !!!if( .not. state%nlSolver%implicitly .and. Newton%iter == 1) then

          ! for iter > 1, array b(:) = F(x^k) was already computed
          if(iter == 1) then
             !print*,'Compute  ONLY vector q(w)'
             !write(600,*) '# Compute  ONLY vector q(w)'
             call ComputeTermsSP( )
             vector_update = .true.
          else
             vector_update = .false.
             Newton%b(1:state%nsize) = Newton%b1(1:state%nsize)
             Newton%res = Newton%res0
          endif


       endif ! end of if(state%nlSolver%implicitly)


       !eta = state%time%alpha(0)/state%time%tau(1)
       eta = 1./state%time%tau(1)


       if(state%time%time_method == 'E') then
          ! explicit discretization, direct solution of the system via
          ! multiplication of the Mass Matrix inversion
          call SolveBlockDiagonalSystem(eta, Newton%norm_res)

          print*,'VDVDVDVD',Newton%norm_res

       else if(state%time%time_method == 'P') then
          ! pseudo-time stepping  NOT WORKING

           ! if(iter == 1) then
           !    do i=1,grid%nelem
           !       elem => grid%elem(i)
           !       do j=1,elem%flen
           !          l = elem%face(neigh, j)
           !          if(l > 0) elem%max_eigenvals = max(elem%max_eigenvals, grid%elem(l)%max_eigenvals )
           !       enddo

           !    enddo
           ! endif


          if(vector_update) then
             call FillVector(Newton%b(1:state%nsize), eta ) ! value of eta is not used, only its sign
             !write(*,'(a6,i5, 40es12.4)') 'flux:',0,Newton%b(1:state%nsize)
             !Newton%res  = VectorNorm(Newton%b)
             Newton%res  = VectorPrecondNorm(Newton%b)
          endif


          ! size of the pseudo-time step
          Newton%lambda = 0.02 / state%max_eigenvals

          !print*,'pseud tau =',  Newton%lambda
          !pseudo-time update
          Newton%res0  = 0
          k = 0
          do i=1,grid%nelem
             elem => grid%elem(i)

             !Newton%lambda = 0.005 / elem%max_eigenvals
             Newton%lambda = 0.01 / (elem%max_eigenvals *  state%max_eigenvals)**0.5

             j = elem%dof

             do l=1,ndim
                i1 = (l-1)*j + 1
                i2 = l*j


                elem%w(0,i1:i2) = elem%w(0,i1:i2)  &
                     !+ Newton%lambda *  Newton%b(k+1:k+j)
                     + Newton%lambda * matmul(elem%MassInv%Mb(1:j,1:j),  Newton%b(k+1:k+j) )

                k = k+j

             enddo
          enddo

          !stop

          Newton%lambda_old = Newton%lambda

          !print*,'Compute ONLY vector q(w)', state%nlSolver%implicitly
          !write(600,*) '# Compute ONLY vector q(w)', state%nlSolver%implicitly
          call ComputeTermsSP( )
          call FillVector(Newton%b1(1:state%nsize),eta ) !include time derivative term

          !!!  VD  VD
          Newton%res0 = VectorPrecondNorm( Newton%b1 )
          !!!Newton%res0  = (Newton%res0 / state%nsize)**0.5

          Newton%theta  =  Newton%res0  / Newton%res

          write(*,'(a6,i5,20es18.10)') 'Pseud:',iter,  Newton%res,  Newton%res0,  Newton%theta, Newton%lambda
          !print*

       else  ! fully implicit or semi-implicit methods

          if(vector_update) then
             call FillVector(Newton%b(1:state%nsize), eta ) ! value of eta is not used, only its sign
             !Newton%res  = VectorNorm(Newton%b)
             Newton%res  = VectorPrecondNorm(Newton%b)
          endif


          if(iter == 1) Newton%res_ref = Newton%res

          ! first exit based on the absolute value of the residuum vector || F(w) ||
          if(state%time%time_method == 'I' .and. iter > 1  &
               .and. Newton%res/state%nsize < Newton%tol) goto 200

          Newton%iter  = Newton%iter + 1
          Newton%Aiter = Newton%Aiter  + 1

          Newton%lambda_old = 0.
          Newton%x(1:state%nsize) = 0.

          if (.not. state%linSolver%tol_fixed) then
             if(iter == 1) then
                !state%linSolver%tol = 0.5
                state%linSolver%tol = 0.25
             else
                !state%linSolver%tol = min(0.5, 0.5*(Newton%res / Newton%res1)**1.5)
                !state%linSolver%tol = min(0.5, 0.25*(Newton%res / Newton%res1)**2.0)
                !state%linSolver%tol = 0.25**Newton%iter
                state%linSolver%tol = min(0.25, abs(Newton%res - state%linSolver%residuum) / Newton%res1)

                !write(*,'(a4,2es12.4,a1,4es12.4)') '@@@@', &
                !     state%linSolver%tol,  min(0.5, 0.25*(Newton%res / Newton%res1)**2.0), &
                !     '|', abs(Newton%res - state%linSolver%residuum) / Newton%res1, &
                !     Newton%res, state%linSolver%residuum ,  Newton%res1
                state%linSolver%tol = max(state%linSolver%tol, 1E-6)
             endif
          endif

          !print*,'Solution of linear algebraic system'
          eta = state%time%alpha(0)/state%time%tau(1)

          call cpu_time(t2)
          time_prepare = time_prepare + t2 - t1

          !!write(*,'(a6, 3es12.4)') 'timesA',time_prepare, time_solve , time_estim

          !call WriteMblock(grid%elem(1)%block(0) )       !call WriteMatrixA(0.)

          if ( state%space%adapt%adapt_method == 'ALG' .or. &
                   (state%space%adapt%adapt_method == 'ALG2' .and. ( .not. state%first_GMRES_conv) ) ) then

             if (state%space%adapt%adapt_level > 0 .and. state%time%iter_loc == 1) then  !computation of momentums in the first iteration step in an adaptation
                 !Print*, '*************V cyklu*********'
                 if (.not. state%loc_RTN(state%space%deg)%defined ) then
                    Print*, 'Error, loc_RTN not defined yet!!'
                    stop
                 endif

                 do i=1, grid%nelem
                    call InitElemRTNMomentumCD(grid%elem(i), state%space%deg, SetRTNdof(state%space%deg) )
                 enddo
             endif
          endif

          if( state%MGsolver ) then
             !!!call InitMG() moved

              ! *** TEST iLU as a solver ***
              !call RunILUSolu(Newton%x,Newton%b,state%nsize)
              !print*,L2Res(Newton%x,Newton%b,state%nsize,bMVprod)

              ! *** COMPARE Jacobi and G-S ***
              !allocate(y(state%nsize))
              !do  i=1,50,1
              !    call bJacobi(Newton%x,Newton%b,state%nsize,y)
              !    Newton%x = y
              !    !call bGS(Newton%x,Newton%b,state%nsize)
              !    print*,i,L2Res(Newton%x,Newton%b,state%nsize,bMVprod)
              !end do
              !deallocate(y)

              !print*,'--STOP'
              !STOP

              ! *** ILU for MG *** part of InitrMG
              !do i=1,grid%nelem
              !   call InitMGMatrixElementShape(grid%elem(i))
              !end do
              !call ComputeBlockMGILUPrecond()

             !call cpu_time(t1)
             call SolveBlockLinearProblem(state%nsize, eta, Newton%b, Newton%x, &
                  state%linSolver%residuum, state%linSolver%iter, precond_update, state%linSolver%lin_solver_not_conv)
             !call cpu_time(t2)

            !\linear multigrid
          else
             !write(*,'(a6,4es14.6)')'AAA0', eta, Newton%res, state%linSolver%tol

             !MatrixA = 'MatrixA'
             !open(21, file=MatrixA, status='UNKNOWN', position = 'append')
             !do i=1,grid%nelem
             !   do j=1, grid%elem(i)%dof
             !      write(21,'(15es14.6)') grid%elem(i)%Mass%Mb(j, 1:grid%elem(i)%dof)
             !   enddo
             !enddo
             !write(21,'(a10)') '**********'
             !close(21)



             !MatrixB = 'MatrixB'
             !open(23, file=MatrixB, status='UNKNOWN', position = 'append')
             !do i=1,grid%nelem
             !   do j=1, grid%elem(i)%dof
             !      write(23,'(15es14.6)') grid%elem(i)%block(0)%Mb(j, 1:grid%elem(i)%dof)
             !   enddo
             !enddo
             !write(23,'(a10)') '**********'
             !close(23)

             !MatrixC = 'MatrixC'
             !open(24, file=MatrixC, status='UNKNOWN', position = 'append')
             !do i=1,grid%nelem
             !   do j=1, grid%elem(i)%flen
             !      k = grid%elem(i)%face(neigh,j)
             !      if (k > 0) then
             !         do l=1, grid%elem(i)%dof
             !            write(24,'(15es14.6)') grid%elem(i)%block(j)%Mb(l, 1:grid%elem(k)%dof)
             !         enddo
             !      endif
             !   enddo
             !enddo
             !write(24,'(a10)') '**********'
             !close(24)

             call SolveBlockLinearProblem(state%nsize, eta, Newton%b, Newton%x, &
                  state%linSolver%residuum, state%linSolver%iter, precond_update, state%linSolver%lin_solver_not_conv)
             !!write(*,'(a6,5es14.6)')'AAA1', eta, time_prepare, time_solve, t1, t2

             !RHSideb = 'RHSideb'
             !open(22, file=RHSideb, status='UNKNOWN', position = 'append')
             !do i=1,state%nsize
             !   write(22,'(es14.6)') Newton%b(i)
             !enddo
             !write(22,'(a10)') '**********'
             !close(22)

          endif

          call cpu_time(t1)
          time_solve = time_solve + t1 - t2

          if(state%time%time_method == 'S') then ! semi-implicit scheme (NO damping !!)

             call bMVprod(Newton%rr, Newton%x, state%nsize)    ! rr = Ax

             Newton%rr(:) = Newton%rr(:) - Newton%b(:)                ! rr = Ax - b

             ! update of the solution
             k = 0
             do i=1,grid%nelem
                elem => grid%elem(i)
                j = elem%dof*ndim
                elem%w(0,1:j) = elem%w(0,1:j)  + Newton%x(k+1:k+j)

                ! components of th residual vector
                elem%vec( res_vec, 1:j) = - Newton%rr(k+1:k+j)

                k = k+j
             enddo


          else  ! fully implicit method

             ! residuum of the linear algebraic system
             eta = state%time%alpha(0)/state%time%tau(1)

             call bMVprod(Newton%rr, Newton%x, state%nsize)    ! rr = Ax

             Newton%rr(:) = Newton%rr(:) - Newton%b(:)                ! rr = Ax - b

             !rel_res = VectorPrecondNorm(Newton%rr(:) ) / Newton%res
             !rel_res = state%linSolver%residuum / VectorPrecondNorm(Newton%b(:) )

             ! seeking of optimal damping parameter
             Newton%lambda = 1.0   ! initialization of lambda (= damping factor)

             do l=1,10    ! iterations, seeking the optimal damping factor
             !do l=1,1    ! iterations, seeking the optimal damping factor
                ! update of the solution
                k = 0
                do i=1,grid%nelem
                   elem => grid%elem(i)
                   j = elem%dof*ndim
                   elem%w(0,1:j) = elem%w(0,1:j)  &
                        + (Newton%lambda-Newton%lambda_old) * Newton%x(k+1:k+j)

                   ! components of th residual vector
                   elem%vec( res_vec, 1:j) = - Newton%rr(k+1:k+j)


                   k = k+j
                enddo
                Newton%lambda_old = Newton%lambda

                !print*,'Compute ONLY vector q(w)', state%nlSolver%implicitly
                !write(600,*) '# Compute ONLY vector q(w)', state%nlSolver%implicitly
                call ComputeTermsSP( )
                call FillVector(Newton%b1(1:state%nsize),eta ) !include time derivative term

                Newton%res0 = VectorPrecondNorm( Newton%b1 )
                Newton%theta  =  Newton%res0  / Newton%res


                !write(*,'(a8,i5,20es14.6)') 'newton:',iter,  Newton%res,  Newton%res0,  Newton%theta

                !call cpu_time(t2)

                write(55,'(3i4,13es11.3, i5,2l3,F10.1)') state%time%iter,Newton%iter,l, &   ! 1..3
                     state%time%iter + 1.*(Newton%iter-1)/(Newton%max_iter+2), &      ! 4
                     state%space%adapt%adapt_level + 3.*state%nlSolver%Aiter /(Newton%max_iter+2)/state%time%maxiter, &
                     Newton%res/state%nsize,  Newton%res0/state%nsize, &         ! 6..7
                     Newton%res,  Newton%res0, &                                 ! 8..9
                     rel_res, state%linSolver%tol,  Newton%lambda, Newton%theta, &        !10 .. 13
                     state%L_estim(resA:resT),  &                                !14 .. 17
                     state%linSolver%iter, loc_implicitly, &
                     precond_update, t2- time_solve

                if(loc_implicitly) then
                   write(54,'(3i4,12es11.4, i5,2l3)') state%time%iter,Newton%iter,l, &   ! 1..3
                        state%time%iter + 1.*(Newton%iter-1)/(Newton%max_iter+2), &      ! 4
                        state%space%adapt%adapt_level+3.*state%nlSolver%Aiter/(Newton%max_iter+2)/state%time%maxiter, &
                        Newton%res/state%nsize,  Newton%res0/state%nsize, &         ! 6..7
                        Newton%res,  Newton%res0, &                                 ! 8..9
                        rel_res, state%linSolver%tol, Newton%lambda, Newton%theta,  &
                        state%L_estim(resS), state%L_estim(resA), &
                        state%linSolver%iter,loc_implicitly,precond_update

                endif

                newton_count = newton_count + 1

                !if(loc_implicitly) print*,'update was applied'
                !print*,'newton count = ', newton_count

                ! residuum is not decreasing, reduce lambda
                if(Newton%theta >= 1 .and. Newton%res0 > 1E-8 ) then
                !if(Newton%theta > 1 - Newton%lambda/4  .and. Newton%res0 > 1E-8 )then
                   !Newton%lambda = min(mu1, 0.5* Newton%lambda)
                   Newton%lambda = 0.5* Newton%lambda
                else
                   Newton%lambda1 = min(1., 1.5*Newton%lambda )
                   !
                   !if(Newton%lambda1 >= 4 * Newton%lambda) then  ! lambda was to short
                   !   Newton%lambda = Newton%lambda1
                   !else
                      !if(Newton%theta > 1./4) then
                   !   state%nlSolver%implicitly = .true.
                   !endif
                   !Newton%lambda = Newton%lambda1

                   goto 15
                   !endif
                endif

                !if(Newton%lambda < 1E-2) goto 15  ! too small lambda, NOT CONVERGE
                if(Newton%lambda < 2E-1) goto 15  ! too small lambda, NOT CONVERGE

             enddo   !l=1,10    ! iterations, seeking the optimal damping factor

             !!state%nlSolver%implicitly = .true. ! last Newton iteration was not sucessfull
15           continue

             if(newton_count >= 8) state%nlSolver%implicitly = .true.
             !if(newton_count >= 4) state%nlSolver%implicitly = .true.
             !!!if(Newton%theta > 0.75) state%nlSolver%implicitly = .true.
             if(Newton%theta > 0.5) state%nlSolver%implicitly = .true.

             !if(Newton%lambda < 0.9) state%nlSolver%implicitly = .true. !lambda is small, C(w) update
             if(state%linSolver%lin_solver_not_conv > 0) state%nlSolver%implicitly = .true.

             !if(state%nlSolver%implicitly ) write(83,'(a8,i5, a14,i5,a8,es12.4,a8,l2 )')&
             !     'iter = ',state%time%iter, '### Ncount =',newton_count, &
             !     ',  theta = ', Newton%theta,',  non conv', state%linSolver%lin_solver_not_conv


             ! nonlinear algebraic residuum criterion
             if(Newton%non_alg_stop == 'aRES') then
                !!if(Newton%tol2 > 0. .and. (state%space%adapt%adapt_method == 'RES' .or. state%space%adapt%adapt_method == 'ANI' &
                !!  .or. state%space%adapt%adapt_method == 'Ahp') ) then

                call cpu_time(t2)
                time_prepare = time_prepare + t2 - t1

                call RezidErrorEstimates( .true., .false. )

                call cpu_time(t1)
                time_estim = time_estim + t1 - t2

                open(91, file='criter', status='unknown', position='append')
                if(Newton%iter == 1) write(91,'(x)')
                if(Newton%Aiter == 1) write(91,'(x)')
                write(91,'(2i6, 50es12.4)' ) state%time%iter, Newton%iter, &
                     state%time%iter + 1.*(Newton%iter-1)/(Newton%max_iter+2), &
                     state%space%adapt%adapt_level + 3.* Newton%Aiter /(Newton%max_iter+2)/state%time%maxiter, &
                     state%linSolver%tol, &   !5
                     state%L_estim(resA:resST),  &
                     state%estim(resA_ST,1), state%estim(resA_ST_loc,1),   &      ! 10..11
                     Newton%norm_res,  Newton%res0, Newton%res0/state%nsize, &      ! 12..14
                     Newton%tol, Newton%tol2,  t2 - state%start_time              ! 15..17


                close(91)

                !if(state%L_estim(resA)/state%G_estim(15) < Newton%tol2 ) goto 200
                ! VD
                !print*
                !write(*,'(a26,3i5,20es14.6)') '!AW! iter tol resA/resST ',state%time%iter, Newton%iter, 0, &
                !     Newton%tol2,  state%estim(resA_ST,1)  !, & !  Newton%res/state%nsize, Newton%tol, &
                !!state%L_estim(resA), state%L_estim(resS) !, state%L_estim(resST)
                !print*,'###################################################################'

                if(state%estim(resA_ST,1) < Newton%tol2 ) then
                !if(state%estim(resA_ST_loc,1) < Newton%tol2 ) then
                   call cpu_time(t2)
                   time_prepare = time_prepare + t2 - t1
                   !!write(*,'(a6, 3es12.4)') 'timesE',time_prepare, time_solve , time_estim
                   goto 200
                endif

             else  ! Newton%non_alg_stop == 'precL2'

                ! other techniques
                !if(Newton%res0/Newton%res_ref < Newton%tol2 ) goto 200

                !if(ndim == 4) &
                !     call DirectComputeDrag_Lift(cDLM(Newton%iter, 1:5 )  )
                !diffs = ( (cDLM(Newton%iter, 1)  - cDLM(Newton%iter-1, 1) )**2 &
                !     + (cDLM(Newton%iter, 2)  - cDLM(Newton%iter-1, 2) )**2  &
                !     + (cDLM(Newton%iter, 3)  - cDLM(Newton%iter-1, 3) )**2 )**0.5
                !
                !open(44, file='cDLM', status='UNKNOWN', position = 'append')
                !write(44,'(2i5,15es11.3)') state%time%iter, Newton%iter, &
                !     state%time%iter + 1.*(Newton%iter)/(Newton%max_iter+1), &
                !     cDLM(Newton%iter, 1:3), &
                !     abs(cDLM(Newton%iter, 1)  - cDLM(Newton%iter-1, 1) ), &
                !     abs(cDLM(Newton%iter, 2)  - cDLM(Newton%iter-1, 2) ), &
                !     abs(cDLM(Newton%iter, 3)  - cDLM(Newton%iter-1, 3) ), diffs
                !close(44)

                if(Newton%res0 < Newton%tol * state%nsize )  then
                   call cpu_time(t2)
                   time_prepare = time_prepare + t2 - t1
                   goto 200
                endif

             endif  ! end of other technique
          endif  ! end of semi-implicit / implicit

       endif    ! end of explicit

       ! go to the next time step
       Newton%res1 = Newton%res

       call cpu_time(t2)
       time_prepare = time_prepare + t2 - t1
       !write(*,'(a6,4es14.6)')'AAA3', time_prepare, time_solve, t1, t2
       !!write(*,'(a6, 3es12.4)') 'timesF',time_prepare, time_solve , time_estim

    enddo  ! end of Newton iterations: Newton%iter = 1, Newton%max_iter

200 continue

    write(55,'(3i4,12es11.3, i5,2l3)') state%time%iter,Newton%iter,l, &   ! 1..3
         state%time%iter + 1.*(Newton%iter+0)/(Newton%max_iter+2), &      ! 4
         state%space%adapt%adapt_level + 3.*(state%nlSolver%Aiter+1) /(Newton%max_iter+2)/state%time%maxiter, &
         Newton%res0/state%nsize,  Newton%res0/state%nsize, &         ! 6..7
         Newton%res0,  Newton%res0

    if(loc_implicitly) &
         write(54,'(3i4,12es11.4, i5,2l3)') state%time%iter,Newton%iter,l, &   ! 1..3
         state%time%iter + 1.*(Newton%iter-1)/(Newton%max_iter+2), &      ! 4
         state%space%adapt%adapt_level+3.*state%nlSolver%Aiter/(Newton%max_iter+2)/state%time%maxiter, &
         Newton%res/state%nsize,  Newton%res0/state%nsize, &         ! 6..7
         Newton%res,  Newton%res0, &                                 ! 8..9
         rel_res, state%linSolver%tol, Newton%lambda, Newton%theta,  &
         state%L_estim(resS), state%L_estim(resA),state%linSolver%iter,loc_implicitly,precond_update


    !write(*,'(a6,3i5,20es12.4)') '????',state%time%iter, Newton%iter, 0, &
    !     state%estim(min_resA_ST,1), Newton%tol2 ,  Newton%res/state%nsize, Newton%tol


    ! space, time, space-time error estimates
    if( state%space%estim_space == 'RES' &  !!!.or. state%space%estim_space == 'ANI' &
         .or. state%space%estim_space == 'Ahp' ) then
       call RezidErrorEstimates( .false., .false. )
       call cpu_time(t1)
       time_estim = time_estim + t1 - t2
    endif
    !!write(*,'(a6, 3es12.4)') 'timesX',time_prepare, time_solve , time_estim

    if(state%time%time_method == 'S') then
       state%nlSolver%implicitly = .false.
       !print*,'Computation of the residuum'
       call ComputeTermsSP( )
       state%nlSolver%implicitly = .true.
       Newton%norm_res = EvalSSresidExplicit( )
       !Newton%norm_res = EvalSSresid( )

       !!call CheckResiduum( )

    elseif(state%time%time_method == 'I') then
       loc_implicitly = state%nlSolver%implicitly
       state%nlSolver%implicitly  = .false.

       call ComputeTermsSP( )
       Newton%norm_res = EvalSSresidExplicit( )
       state%nlSolver%implicitly = loc_implicitly
    endif

    !write(*,*) ' #### Total number of newton updates = ', Newton%updates, Newton%norm_res
    !write(111,*) state%time%iter+1, Newton%norm_res
    if(state%time%time_method == 'I') then
       write(54,*) ' #### Total number of newton updates = ', Newton%updates
       write(55,*) ' '
       close(54)
       close(55)
       !close(56)
    endif

    deallocate(cDLM)

    !if(mod(state%time%iter, 3) == 0) then
    !   state%isol = state%isol + 1
    !   call WriteProgressOutput( 'A' )
    !endif

    state%CPU_prepare = state%CPU_prepare  + time_prepare
    state%CPU_solve   = state%CPU_solve    + time_solve
    state%CPU_estim   = state%CPU_estim    + time_estim

    !VD
    !write(*,'(a6, 4es12.4)') 'timXXX',time_prepare, time_solve , time_estim, &
    !     time_prepare+time_solve+time_estim
    !write(*,'(a6, 3es12.4)') 'timXXX',state%CPU_prepare, state%CPU_solve , state%CPU_estim


  end subroutine PerformOneTimeStepSP

!> evaluation of the matrix and vector blocks
 subroutine ComputeTermsSP( )
   integer :: i
   ! clearing of the appropriate arrays
   if(state%nlSolver%implicitly) call ClearMatrixBlocksSP()
   call ClearVectorBlocksSP()

   ! setting of the corresponding fluxes
   if(state%modelName == 'incNS') then
      call ComputeElementsTerms(Set_f_s_scalar, Set_A_s_scalar, Set_Ppm_scalar, &
           Set_R_s_scalar, Set_K_sk_scalar, Set_S_scalar, Set_DS_scalar)


   else
      print*,'Compute Elements Terms not implemented for ', state%modelName
   endif

   !do i=1,grid%nelem
   !   write(*,'(a4,i5,20es12.4)') 'w: ',i, grid%elem(i)%w(0,:)
   !   call WriteMblock(grid%elem(i)%block(0) )
   !enddo


 end subroutine ComputeTermsSP


  !> clearing the matrix blocks
  subroutine ClearMatrixBlocksSP( )
    class(element), pointer :: elem
    integer :: i, j

    do i = 1, grid%nelem
       elem => grid%elem(i)

       elem%SPblock(bVV,0)%Mb(:,:) = 0. ! diagonal blocks
	elem%SPblock(bVP,0)%Mb(:,:) = 0.
	elem%SPblock(bPV,0)%Mb(:,:) = 0.

       do j = 1, elem%flen         ! off-diagonal blocks
          if(elem%face(neigh,j) > 0)  then
		elem%SPblock(bVV,j)%Mb(:,:) = 0.
		elem%SPblock(bVP,j)%Mb(:,:) = 0.
		elem%SPblock(bPV,j)%Mb(:,:) = 0.
	  endif
       enddo
    enddo
  end subroutine ClearMatrixBlocksSP

  !> clearing the vector blocks
  subroutine ClearVectorBlocksSP( )
    integer :: i

    do i = 1, grid%nelem
       grid%elem(i)%vec(rhs,:) = 0. ! clearing of rhs terms, e.g., BC
    enddo

  end subroutine ClearVectorBlocksSP

end module compute_operSP
