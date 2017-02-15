!> various problem operations
module problem_oper
  use main_data
  use f_mapping
  use mesh_oper
  use blocks_integ
  use model_oper
  use integration
  use rav_tho_nedelec
  use io_sub
  use eval_sol
  use set_solution
  use basis
  use loc_rav_tho_ned
 ! use dwr_oper
  implicit none



  !public:: InitProblem - new version Init_Problem
  public:: PrepareProblem

  public:: TimeDerivativeVector
  public:: AddTimeDerivativeVector

  public:: PrepareOneElement
!  public:: SetElementQuadraturesDegrees
  public:: InitMatrixElementShape
  public:: DeleteMatrixElementShape
  public:: InitGlobalMatrixShape
  public:: GetBlocksSizeElem
  public:: GetBlocksSizeElemST
  public:: ComputeLocalMassMatrix

  public:: InitMG
  public:: DeInitMG
  public:: InitMGMatrixElementShape

  public:: ComputeRefTimeMatrix
  public:: UpdateMesh_p
  public :: PlotSolPrimal
  public :: PlotSolDual
!  public:: SetEdgeDegrees
contains


 !> compute the vector: elem\%vec(i) + =
 !> \f$\frac{1}{\tau_k} ({\bf w}_h^k - {\bf w}_h^{k-1}, \varphi_i)_{L^2(K)} \f$
 !> \f$i=1,\dots, DOF^+_K,\  K \f$ = elem \f$ \in T_h\f$, similarly for higher order BDF
 subroutine TimeDerivativeVector( )
   real, dimension(:,:), pointer :: phi      ! pointer to the test function
   real, dimension(:,:), allocatable :: wi   ! w^k - w^{k-1} in integ nodes
   class(element), pointer :: elem      ! elem = element
   integer :: dof, dofA, Qdof
   integer :: i, k, kst, n, ie

   associate ( time => state%time)
   select type ( time )
      class is ( TimeBDF_t )

         allocate(wi(1:state%space%max_Qdof,1:ndim) )

         do ie=1,grid%nelem
            elem => grid%elem(ie)
            dof = elem%dof
            dofA = dof
            if(elem%deg_plus) dofA = elem%dof_plus

            Qdof = elem%Qdof

            phi => state%space%V_rule(elem%Qnum)%phi(1:dofA, 1:Qdof)

            ! setting of w^k - w^{k-1} in integration nodes
            do k=1,ndim              ! k = index of component of w
               kst = dof*(k-1) + 1
               do i=1,Qdof
                  wi(i,k) = 0.
                  do n=0,time%deg
                     !if(elem%i == 2 .and. k==1 .and. i == 1) &
                     !     write(*,'(a4,3i5,12es12.4)') '@@@', n, time%deg,dofA, &
                     !     time%alpha(n)!, elem%w(n,kst:kst+dof-1), elem%vec(rhs,1)

                     wi(i,k) = wi(i,k) - time%alpha(n) &
                          * dot_product(elem%w(n,kst:kst+dof-1), phi(1:dof,i) )
                  enddo
               enddo
               !if(elem%i == 18) then
               !  write(197,*) 'k =   ', k,'     elem%Qnum = ',elem%Qnum
               !   write(197,'(a6, 80es12.4)') 'wi-ii:', wi(:,k)
               !   write(197,'(a6, 80es12.4)') 'wi-ii:', wi(:,k)/state%time%tau(1)
               !   write(197,*) sum(wi(:,k) )/time%tau(1)
               !endif

            enddo

            wi(1:Qdof, 1:ndim) = wi(1:Qdof, 1:ndim) / time%tau(1)

            call EvalVectorB(elem, wi(1:Qdof, 1:ndim), dofA, elem%vec(rhs, 1:dofA*ndim) )

         enddo

         deallocate(wi)

      class default
         stop 'Subroutine TimeDerivativeVector is implemented for BDF only.'
   end select
   end associate

 end subroutine TimeDerivativeVector


 !> compute the vector: elem\%vec(i) + =
 !> \f$\frac{1}{\tau_k} ({\bf w}_h^k - {\bf w}_h^{k-1}, \varphi_i)_{L^2(K)} \f$
 !> \f$i=1,\dots, DOF^+_K,\  K \f$ = elem \f$ \in T_h\f$, similarly for higher order BDF
 !> ONLY ONE ELEMENT
 subroutine TimeDerivativeVector_ONLY_ONE_ELEM(elem )
   class(element), target, intent(inout) :: elem      ! elem = element
   real, dimension(:,:), pointer :: phi      ! pointer to the test function
   real, dimension(:,:), allocatable :: wi   ! w^k - w^{k-1} in integ nodes
   integer :: dof, dofA, Qdof
   integer :: i, k, kst, n, ie

   associate ( time => state%time)
   select type ( time )
      class is ( TimeBDF_t )

         allocate(wi(1:state%space%max_Qdof,1:ndim) )

         dof = elem%dof
         dofA = dof
         if(elem%deg_plus) dofA = elem%dof_plus

         Qdof = elem%Qdof

         phi => state%space%V_rule(elem%Qnum)%phi(1:dofA, 1:Qdof)

         ! setting of w^k - w^{k-1} in integration nodes
         do k=1,ndim              ! k = index of component of w
            kst = dof*(k-1) + 1
            do i=1,Qdof
               wi(i,k) = 0.
               do n=0,time%deg
                  !if(elem%i == 2 .and. k==1 .and. i == 1) &
                  !     write(*,'(a4,3i5,12es12.4)') '@@@', n, time%deg,dofA, &
                  !     time%alpha(n)!, elem%w(n,kst:kst+dof-1), elem%vec(rhs,1)

                  wi(i,k) = wi(i,k) - time%alpha(n) &
                       * dot_product(elem%w(n,kst:kst+dof-1), phi(1:dof,i) )
               enddo
            enddo
            !if(elem%i == 18) then
            !  write(197,*) 'k =   ', k,'     elem%Qnum = ',elem%Qnum
            !   write(197,'(a6, 80es12.4)') 'wi-ii:', wi(:,k)
            !   write(197,'(a6, 80es12.4)') 'wi-ii:', wi(:,k)/state%time%tau(1)
            !   write(197,*) sum(wi(:,k) )/time%tau(1)
            !endif

         enddo

         wi(1:Qdof, 1:ndim) = wi(1:Qdof, 1:ndim) / time%tau(1)

         call EvalVectorB(elem, wi(1:Qdof, 1:ndim), dofA, elem%vec(rhs, 1:dofA*ndim) )

         deallocate(wi)

      class default
         stop 'Subroutine TimeDerivativeVector_ONLY_ONE_ELEM is implemented for BDF only.'
      end select
    end associate

  end subroutine TimeDerivativeVector_ONLY_ONE_ELEM


 !> compute the vector: elem\%vec(i) + =
 !> \f$ ( \partial{\bf w}_{h\tau}(\bar{t}), \varphi_i)_{L^2(K)} \f$
 !> \f$i=1,\dots, DOF^+_K,\  K \f$ = elem \f$ \in T_h\f$,
 !>  \f${\bf w}_{h\tau}(\bar{t}) \f$ ( pw polynomial) stored in elem\%w(1,*) !!!!!
 subroutine AddTimeDerivativeVector( )
   real, dimension(:,:), pointer :: phi      ! pointer to the test function
   real, dimension(:,:), allocatable :: wi   ! w^k - w^{k-1} in integ nodes
   class(element), pointer :: elem      ! elem = element
   integer :: dof, dofA, Qdof
   integer :: i, k, kst, ie

   allocate(wi(1:state%space%max_Qdof,1:ndim) )

   do ie=1,grid%nelem
      elem => grid%elem(ie)
      dof = elem%dof
      dofA = dof
      if(elem%deg_plus) dofA = elem%dof_plus

      Qdof = elem%Qdof

      phi => state%space%V_rule(elem%Qnum)%phi(1:dofA, 1:Qdof)

      ! setting of w^k - w^{k-1} in integration nodes
      do k=1,ndim              ! k = index of component of w
         kst = dof*(k-1) + 1
         do i=1,Qdof
            wi(i,k) = -dot_product(elem%w(1,kst:kst+dof-1), phi(1:dof,i) ) ! already time der
            ! --------^  minus is correct, since the flux vector is on RHS (?)
         enddo
      enddo

      call EvalVectorB(elem, wi(1:Qdof, 1:ndim), dofA, elem%vec(rhs, 1:dofA*ndim) )
   enddo

   deallocate(wi)

 end subroutine AddTimeDerivativeVector


  subroutine WriteQuadratures()
    integer :: i, j, ifile

    ifile = 12
    open(ifile, file='Gnodes.tex', status='UNKNOWN')

    write(ifile,*) '{\footnotesize '
    !write(ifile,*) '{\scriptsize '
    write(ifile,*) '\begin{tabular}{l|ccc|ccc}'
    write(ifile,*) '\hline'
    write(ifile,*) ' $G_k$ & $j$ & $w_j$  & $x_j$ &$j$ & $w_j$  & $x_j$\\'
    write(ifile,*) '\hline'
    write(ifile,*) '\hline'

    !do i=1,maxGrule
    do i=1,12
       !write(ifile,'( a4,i2,a20)') '$G_{',i,'} $ & & & & & \\'

       do j = 1, state%space%G_rule(i)%Qdof/2
          write(ifile,'(a4,i2,a3, 2(a3,i3,a3,f16.14,a3,f16.14),a3)') &
               '$G_{',i,'} $', &
               ' &', 2*(j-1)+1,'&', state%space%G_rule(i)%weights(2*(j-1)+1),' & ', &
               state%space%G_rule(i)%lambda(2*(j-1)+1), ' & ', &
               2*j,' & ', state%space%G_rule(i)%weights(2*j),' & ', &
               state%space%G_rule(i)%lambda(2*j),' \\'
       enddo
       do j = state%space%G_rule(i)%Qdof/2 * 2 + 1, state%space%G_rule(i)%Qdof
          write(ifile,'(a4,i2,a3, 1(a3,i3,a3,f16.14,a3,f16.14),a3)') &
               '$G_{',i,'} $', &
               ' &', j,'&', state%space%G_rule(i)%weights(j),' & ', &
               state%space%G_rule(i)%lambda(j), ' \\'

       enddo
       write(ifile,*) '\hline'
    enddo

    write(ifile,*) '\end{tabular}'
    write(ifile,*) '}'

    close(ifile)

    ! volume quadratures
    open(ifile, file='Vnodes.tex', status='UNKNOWN')

    write(ifile,*) '{\footnotesize '
    !write(ifile,*) '{\scriptsize '
    write(ifile,*) '\begin{tabular}{l|ccccc}'
    write(ifile,*) '\hline'
    write(ifile,*)'$D_k$ & $j$ & $w_j$ & $\lambda_{j,1}$ & $\lambda_{j,2}$ & $\lambda_{j,3}$ \\'
    write(ifile,*) '\hline'
    write(ifile,*) '\hline'

    !do i=1,maxVrule
    do i=1,7
       !write(ifile,'( a4,i2,a20)') '$G_{',i,'} $ & & & & & \\'

       do j = 1, state%space%V_rule(i)%Qdof
          write(ifile,'(a4,i2,a5, i3,a3, 4(f16.14,a3))') &
               '$D_{',i,'} $ &', &
               j, '&', &
               state%space%V_rule(i)%weights(j),' & ', &
               state%space%V_rule(i)%lambda(j,1), ' & ', &
               state%space%V_rule(i)%lambda(j,2), ' & ', &
               state%space%V_rule(i)%lambda(j,3), ' \\'
       enddo
       write(ifile,*) '\hline'
    enddo

    write(ifile,*) '\end{tabular}'
    write(ifile,*) '}'

    close(ifile)




!    call Test_V_rules( )
    !stop " call Test_V_rules( ) stOPPED"

    print*,'# WriteQuadratures --- done'
  end subroutine WriteQuadratures

! replaced by init_Problem which is a part of state
!  !> initialization of the problem,
!  !> should be used one time at the beginning of computation
!  !>
!  !> creation of volume and edge quadrature rules
!  !> evaluation of reference test functions at quadrature rules
!  subroutine InitProblem()
!    type(Gauss_rule), pointer :: G_rule
!    integer :: i, j, j1
!
!    state%init_only_matrix_blocks = .false.
!    state%RHS_presented = .true.
!    state%homogenDirichlet = .false.
!    state%num_flux = .true.
!
!    !MinDegreeImplemented = 1
!    !MaxDegreeImplemented = 3
!    !MaxDegreeImplemented = 5
!    !MaxDegreeImplemented = 6
!    !MaxDegreeImplemented = 7
!    !MaxDegreeImplemented = 8
!    !MaxDegreeImplemented = 10
!    !MaxDegreeImplemented = 11  doesn't work, more accurate V_rule required !!!
!
!    if(state%space%deg > MaxDegreeImplemented) then
!       write(*,*) 'The default state%space%deg=',state%space%deg,&
!            ' is higher than MaxDegreeImplemented =', MaxDegreeImplemented
!       stop
!    endif
!
!     ! array for the computation of the hierarchical projection error
!    allocate(state%space%ldeg(-1: MaxDegreeImplemented ) )
!    do i=-1, MaxDegreeImplemented
!       state%space%ldeg(i) =  DOFtriang(i)
!    enddo
!
!
!    ! maximal allocated V and G rules
!    !maxVrule = 23
!    !state%maxVrule = 20
!    !state%maxGrule = 15
!    !maxTrule = state%maxGrule  ! associated with Gauss quadratures
!    !maxRrule = 12 !Radau quadrature rule
!
!    ! for given degree of polynomial approximation gives the degree of Dunavant/Wandzura quadrature
!    allocate(state%space%Qdeg(0:10, 1:2) )
!
!    if(MaxDegreeImplemented > 10) then
!       print*,'Adding of Dunavant/Wandzura quadrature necessary in problem.f90'
!       stop
!    endif
!
!    ! degree of approx: elem%deg: 0   1   2   3   4   5   6   7   8   9  10
!    state%space%Qdeg(0:10, 1)   =    (/ 2,  5,  8, 12, 14, 17, 21, 22, 23, 23, 23 /)
!
!    !maximal level for MG
!    !FIXME: bude este upravene
!    if(state%MGsolver) state%MaxMGlevel = 9
!    !state%MGsolver = .false.
!    !state%MGsolver = .true
!    QnumOffset = maxVrule   !  for the old assessment
!    ! allocation of the V_rules:
!    !                            1 ... state%maxVrule  for triangles
!    !           state%maxVrule + 1 ... state%maxVrule + state%maxGrule  for quadrilaterals
!    allocate(state%space%V_rule(1:maxVrule + maxGrule ) )
!    state%space%V_rule(:)%def = .false.
!
!    ! preparing of volume integration quadratures
!    !do i=1,state%maxVrule
!    do j = 0, MaxDegreeImplemented
!       i = state%space%Qdeg(j, 1)
!       if(.not. state%space%V_rule(i)%def)  call Create_volume_rule(state%space%V_rule(i), i)
!
!       !if(minval(state%space%V_rule(i)%weights(:) ) < 0. ) print *,'???? w',i
!       !if(minval(state%space%V_rule(i)%lambda(:,:)) < 0. ) print *,'???? l',i
!    enddo
!
!    ! testing of the accuracy of quadrature rules on monomials
!    !call Test_V_rules( )
!    !stop
!
!
!    ! preparing of edge integration quadratures
!    allocate(state%space%G_rule(1:maxGrule) )
!    do i=1, maxGrule
!       call Create_Gauss_rule(state%space%G_rule(i), i)
!       call Init_Legenre_pols(state%space%G_rule(i), MaxDegreeImplemented )
!    end do
!
!    ! LEGENDRE polynomials - TEST
!    ! i = 12
!    ! do j=1,state%space%G_rule(i)%Qdof
!    !    write(21, *) state%space%G_rule(i)%lambda(j), state%space%G_rule(i)%Leg_phi(0:MaxDegreeImplemented, j)
!   !  enddo
!    ! G_rule => state%space%G_rule(i)
!    ! do j=0,MaxDegreeImplemented
!    !    do j1=0,MaxDegreeImplemented
!    !       print*,'EDS',j,j1, dot_product(G_rule%weights(:), &
!    !            G_rule%Leg_phi(j, 1:) * G_rule%Leg_phi(j1, 1:) )
!    !    enddo
!    !    print*, G_rule%Leg_phi(j,0)
!    ! enddo
!    ! stop
!
!    ! Bi-Gauss volume integration for quadrilaterals, only weights
!    do i=1, maxGrule
!       call Add_4volume_rule(i,state%space%G_rule(i), state%space%V_rule(i+maxVrule))
!    enddo
!
!    ! quadrature rules for integration over time intervals
!    !MaxTimeDegree = 4
!
!    if (state%time%deg > MaxTimeDegree) then
!       print*, 'Desired time degree is greater than implemented. Stopping'
!       stop
!    endif
!
!
!
!
!    ! preparing of integ rules for time intervals (Gauss)
!    allocate(state%time%T_rule(1:maxTrule) )
!    do i=1,maxTrule
!       call Create_time_rule(state%time%T_rule(i), i)
!    enddo
!
!    if (state%time%disc_time == 'STDG') then
!       !Radau quadrature
!       allocate(state%time%T_rule(1:maxRrule) )
!       do i=1,maxRrule
!          call Create_Radau_rule(state%time%T_rule(i), i)
!       enddo
!    endif
!
!
!    ! Space errors
!    allocate(state%err(1:max_errS) )
!    state%err(:) = 0.
!    state%err(interLq) = 1E+30
!    state%err(interL8) = 1E+30
!
!    ! Space-time errors
!    allocate(state%errSTnorm(1:max_errSTnorm), state%errSTloc(1:max_errSTnorm))
!
!    state%errSTnorm(:) = 0.
!    state%errSTloc(:)  = 0.
!
!    if (state%time%disc_time == 'STDG') then
!       allocate( state%errSnorm( 1:2,1:3) )     !L2/H1 norms in 3 time nodes
!       state%errSnorm(:,:) = 0.
!    endif
!    !call  WriteQuadratures()
!
!
!    ! preparing of volume integration quadratures
!    ! has to be < = QnumOffset, >QnumOffset are for quadrilateralls !!!
!    !maxLrule = MaxDegreeImplemented
!    allocate(state%space%L_rule(0:maxLrule + maxVrule ) )
!
!    do i=0,  maxLrule
!       call Create_L_rule(state%space%L_rule(i), i)
!       call Create_L_rule(state%space%L_rule(i+maxVrule), i+maxVrule)
!    enddo
!
!    ! maximal level of HANGING NODES = 2**max_HGlevel
!    state%space%adapt%max_HGlevel = 2
!    state%space%adapt%HG = 2**(state%space%adapt%max_HGlevel+1)-1
!
!    ! TWO POSSIBILITIES:
!    !call SetBasisFunctions(PHI_Taylor) ! old version INCOMPATIBLE !!!!!!!!!!!
!
!    ! symbolic Gram-Schmidt, for deg <= 5
!    !call SetBasisFunctions(PHI_orthonormal_OLD)
!
!    ! numerical Gram-Schmidt, general deg
!    call ComputeGramSchmidtCoefficients( )
!    call SetBasisFunctions(PHI_orthonormal)
!
!
!    !do i=1,21
!    !   do j =1, state%space%V_rule(15)%Qdof
!    !      write(400+i,*) state%space%V_rule(15)%lambda(j, 1:nbDim), state%space%V_rule(15)%phi(i, j)
!    !   enddo
!    !enddo
!
!    !stop
!
!
!    ! check ComputeElementConvErrors in euler.f90 and
!    ! AddElementMassMatrix in problem.f90
!
!    ! initialization of RTN finite elements
!    !MaxRTNImplemented = 3
!    !MaxRTNImplemented = MaxDegreeImplemented
!
!    ! for local RTN functions
!    allocate( state%loc_RTN(0: MaxDegreeImplemented ) )
!    state%loc_RTN(:)%defined = .false.  ! will be defined during the computation
!
!
!    ! standard RTN functions using reference element
!    allocate( state%RTN(0: MaxRTNImplemented) )
!    do i = 0, MaxRTNImplemented
!       state%RTN(i)%deg = i
!       state%RTN(i)%Qnum = state%space%Qdeg(i, 1) !!! OLD max(1, 2*i)
!       state%RTN(i)%Gnum = i+1
!       call Init_RTN_FE(state%RTN(i), &
!            state%space%V_rule(state%RTN(i)%Qnum), state%space%G_rule(state%RTN(i)%Gnum) )
!
!       allocate(state%RTN(i)%Vnode(1:maxVrule) )
!       state%RTN(i)%Vnode(1:maxVrule)%def = .false.
!       state%RTN(i)%Vnode(1:maxVrule)%defder = .false.
!    enddo
!
!    call InitRGrefinement ( )
!
!    state%print = .false.
!
!    ! initialization for Newton method
!    !allocate(state%nlSolver)  first allocation in main.f90
!    allocate(state%nlSolver%b(1:1), state%nlSolver%b1(1:1), &
!               state%nlSolver%x(1:1), state%nlSolver%rr(1:1) )
!
!    ! initialization of the time step when Stopping Criteria based on AEE are satisfied
!    state%time%iter_SC = -1
!
!    !history of AMA adaptations
!    if(state%space%adapt%adapt_type == 'ANI' .or. state%space%adapt%adapt_type == 'Ahp'.or. state%space%adapt%adapt_type == 'Ihp') &
!         allocate ( state%space%adapt%AMAhistory(0:max(1, state%space%adapt%max_adapt_level), 1:5) )
!
!    !print*,'# Problem initialized'
!
!    state%space%gridS_allocated = .false.
!    state%time%keep_tau = .false.
!
!  end subroutine InitProblem



  ! replaced in o_adaption.f90
  !> initialization of the red green refinement
  subroutine InitRGrefinement()
    integer :: split = 4
    !real, dimension(1:3, 1:3) :: Rlambda

    !given in paramets.f90
    !max_RGlevel = 15  !! defined in state.f90

    allocate(state%space%adapt%RGred(0:split, 1:3, 1:3) )

    state%space%adapt%RGred(0, 1, 1:3) = (/0., 0., 1. /)
    state%space%adapt%RGred(0, 2, 1:3) = (/1., 0., 0. /)
    state%space%adapt%RGred(0, 3, 1:3) = (/0., 1., 0. /)

    call SplitTriangToFour(state%space%adapt%RGred(0, 1:3, 1:3), state%space%adapt%RGred( 1:split, 1:3, 1:3) )

    allocate(state%space%adapt%RGgreen(1:3, 1:nbDim, 1:3, 1:3) ) !!! ONLY for triangles

    !state%RG%green(0)%lambda(1, 1, 1:3) = (/0., 0., 1. /)
    !state%RG%green(0)%lambda(1, 2, 1:3) = (/1., 0., 0. /)
    !state%RG%green(0)%lambda(1, 3, 1:3) = (/0., 1., 0. /)

    ! split edge 1
    state%space%adapt%RGgreen(1, 1, 1, 1:3) = (/0.5, 0., 0.5 /)
    state%space%adapt%RGgreen(1, 1, 2, 1:3) = (/1., 0., 0. /)
    state%space%adapt%RGgreen(1, 1, 3, 1:3) = (/0., 1., 0. /)

    state%space%adapt%RGgreen(1, 2, 1, 1:3) = (/0., 0., 1. /)
    state%space%adapt%RGgreen(1, 2, 2, 1:3) = (/0.5, 0., 0.5 /)
    state%space%adapt%RGgreen(1, 2, 3, 1:3) = (/0., 1., 0. /)


    ! split edge 2
    state%space%adapt%RGgreen(2, 1, 1, 1:3) = (/0.5, 0.5, 0. /)
    state%space%adapt%RGgreen(2, 1, 2, 1:3) = (/0., 1., 0. /)
    state%space%adapt%RGgreen(2, 1, 3, 1:3) = (/0., 0., 1. /)

    state%space%adapt%RGgreen(2, 2, 1, 1:3) = (/1., 0., 0. /)
    state%space%adapt%RGgreen(2, 2, 2, 1:3) = (/0.5, 0.5, 0. /)
    state%space%adapt%RGgreen(2, 2, 3, 1:3) = (/0., 0., 1. /)

    ! split edge 3
    state%space%adapt%RGgreen(3, 1, 1, 1:3) = (/0., 0.5, 0.5 /)
    state%space%adapt%RGgreen(3, 1, 2, 1:3) = (/0., 0., 1. /)
    state%space%adapt%RGgreen(3, 1, 3, 1:3) = (/1., 0., 0. /)

    state%space%adapt%RGgreen(3, 2, 1, 1:3) = (/0., 1., 0. /)
    state%space%adapt%RGgreen(3, 2, 2, 1:3) = (/0., 0.5, 0.5 /)
    state%space%adapt%RGgreen(3, 2, 3, 1:3) = (/1., 0., 0. /)



  end subroutine InitRGrefinement

  ! replaced in o_adaption.f90
  !> split triangle with barycentric coordinates x to four triangles with zi)
  subroutine SplitTriangToFour(x, zi)
    real, dimension(1:3, 1:3), intent(in) :: x
    real, dimension(1:4, 1:3, 1:3), intent(out) :: zi


    zi(1, 1, 1:3) = x(1,1:3)
    zi(1, 2, 1:3) = (x(1,1:3) + x(2,1:3))/2
    zi(1, 3, 1:3) = (x(3,1:3) + x(1,1:3))/2

    zi(2, 1, 1:3) = x(2,1:3)
    zi(2, 2, 1:3) = (x(2,1:3) + x(3,1:3))/2
    zi(2, 3, 1:3) = (x(1,1:3) + x(2,1:3))/2

    zi(3, 1, 1:3) = x(3,1:3)
    zi(3, 2, 1:3) = (x(3,1:3) + x(1,1:3))/2
    zi(3, 3, 1:3) = (x(2,1:3) + x(3,1:3))/2

    zi(4, 1, 1:3) = (x(1,1:3) + x(2,1:3))/2
    zi(4, 2, 1:3) = (x(2,1:3) + x(3,1:3))/2
    zi(4, 3, 1:3) = (x(3,1:3) + x(1,1:3))/2

  end subroutine SplitTriangToFour





  !> preparation of one element \f$ K\in{\cal T}_h\f$; Jacobi matrixes,
  !> preparation of volume quadratures, setting of JF=\f$J_F\f$ (Jacobian) and
  !> D1F = \f$ \left(\frac{D\,F}{D\, x}\right)\f$ in integ. nodes
  !> setting of Poincare constant for element
  subroutine PrepareOneElement( elem )
    type(element), intent(inout):: elem   ! elem = element

!    real, allocatable, dimension (:,:) :: xi, F
!    real, allocatable, dimension (:,:,:) :: DF
    integer :: Gnum,  Qdeg, Qdof, Qnum, ie, je, ic, jc, l, Fdof
    real :: t, n(2)

    call elem%Set_elem_d_gamma( state%space%m_IPG )


    ! setting of Poincare constant
    elem%CP = 1./Pi                ! for CONVEX element

    !Fdof = state%RTN(Fdeg)%dof
    Fdof = SetRTNdof(state%space%deg)
    if( state%space%estim_space == 'pNeu' .or. state%space%estim_space == 'DUA' &
         .or. state%space%estim_space == 'RTN' .or. state%space%estim_space == 'RTNst')  &
         allocate(elem%RTNflux(0:1, 1:ndim, 1:Fdof) )


    elem%errTot = 0.
    elem%errDisc = 0.
    elem%errAlg = 0.
    elem%estTot = 0.
    elem%estDisc = 0.
    elem%estAlg = 0.

    elem%errL2 = 0.
    elem%errH1 = 0.
    elem%errL8 = 0.
    elem%interLq = 0.
    elem%interL8 = 0.
    elem%interH1 = 0.
    elem%estim_loc = 0.

    !! setting of Jacobian and inverse of Jacobi matrix
    ! different for elems with curved edges - non constant Jacobian
    if(elem%F%iFlin ) then
       call ComputeJacobian_constant( elem )
    else
       call ComputeJacobian_curved( elem )
    endif

    call ComputeIntegNode(elem )

    if ( allocated( elem%eta ) ) deallocate(elem%eta)
    allocate(elem%eta(1:max_eta, 1:ndim), source = 0.0 )
    !elem%eta(1:max_eta, 1:ndim) = 0.

    ! MULTIGRID method
    if(state%MGsolver .and. max_MG > 0) &
         allocate(elem%MGw(1:state%MaxMGlevel, 1:max_MG, 1:ndim*elem%dof) )
    elem%MGdeg = elem%deg


    ! AMA technique
    if( state%space%adapt%adapt_type == 'ANI'  &
         .or. state%space%adapt%adapt_type == 'Ahp' &
         .or. state%space%adapt%adapt_type == 'Ihp') then
         if ( allocated( elem%rgabc ) ) deallocate(elem%rgabc)
         allocate ( elem%rgabc(1:3) )
    endif
    ! initial setting for correct function of:
    !  valgrind --leak-check=full ../SRC/AAdgfem shock-ST-AMAhp5a.ini
    elem%errL8 = 0.

    ! lifting operator for the dicrete gradient in pNeu for SIPG and NIPG [Ern, Vohralik, SINUM 15]
    if( state%space%estim_space == 'pNeu') allocate(elem%lifting(1:ndim, 1:nbDim) )
    !print*,'elem%lifting(1:ndim, 1:nbDim) )', elem%i


    elem%ama_p = 0.
    elem%reg2 = 0.

  end subroutine PrepareOneElement

  !> setting of the space dependent parameter, at this moment only scalar
  !> is set to elem%xi(:,:, 0)
  !> used for the Forchheimer flow, setting of nonconstant permeability
  subroutine Setting_of_space_variable_coeffs(gridA)
    type(mesh), intent(inout) :: gridA
    class(element), pointer :: elem
    type(Lagrang_rule), pointer :: L_rule
    real, dimension(:), pointer :: xi
    real :: val1,val2,val3
    real, dimension(:), allocatable :: vals
    integer :: i, je, jc, l, k, Gnum, Ldof

    print*,' subroutine Setting_of_space_variable_coeffs called'

    if(state%modelName == 'scalar' .and.  &
         (state%model%idiff == 13  ) ) then
       ! OLD COULD BE REMOVED
       do i=1, gridA%nelem
          elem => gridA%elem(i)

          ! volume nodes
          do l=1,elem%Qdof
             xi => elem%xi(0, l, 1:2)
             if( state%model%idiff == 13) then
                call  Set_porous_media_A2(xi(2), -xi(1), val1, val2, val3)
                elem%xi(0, l, 0) = val1
             else
                stop "No other implemented (1)"
             endif
          enddo

          ! edge integ nodes
          do je=1, elem%flen
             jc = je
             if(elem%HGnode) jc = elem%HGface(1, je)
             Gnum = elem%face(fGnum,je)

             do l=1,Gnum
                xi => elem%xi(je, l, 1:2)

                if( state%model%idiff == 13) then
                   call  Set_porous_media_A2(xi(2), -xi(1), val1, val2, val3)
                   elem%xi(je, l, 0) = val1
                else
                   stop "No other implemented (2)"
                endif
             enddo
          enddo
       enddo

    elseif(state%modelName == 'porous') then

       allocate( vals(1:iRe) )

       do i=1, gridA%nelem
          elem => gridA%elem(i)

          do je=0, elem%flen
             jc = je

             if(je ==0) then      ! volume integ nodes
                Gnum  = elem%Qdof
                xi => elem%xi(0, l, 1:2)

             else                 ! edge integ nodes
                if(elem%HGnode) jc = elem%HGface(1, je)
                Gnum = elem%face(fGnum,jc)
             endif

             do l=1,Gnum
                xi => elem%xi(jc, l, 1:2)

                if( state%model%idiff == 1) then  ! TEST CASE
                   !call  Set_porous_media_Hraz(xi(1), xi(2), val1, val2, val3)
                   elem%xi(jc, l, 3) =  1.

                elseif( state%model%idiff == 2) then  ! Hraz (damp) by M. Kuraz - only test case
                   call  Set_porous_media_Hraz_OLD(xi(1), xi(2), val1, val2, val3)
                   elem%xi(jc, l, 3) = val1

                elseif( state%model%idiff == 3) then  ! Hraz (damp) by M. Kuraz
                   call  Set_porous_media_Hraz(xi(1), xi(2), iRe, vals(1:iRe) )
                   elem%xi(jc, l, 2+1: 2+iRe-1) = vals(1: iRe-1)

                   !write(102, *)elem%xi(jc, l, 1: 2+iRe-1), sum(elem%xi(jc, l, 2+1: 2+iRe-1) )

                   !call  Set_porous_media_Hraz(elem%xc(1), elem%xc(2), iRe, vals(1:iRe) )
                   !elem%xi(jc, l, 2+1: 2+iRe-1) = vals(1: iRe-1)

                   !if(jc == 0) &
                   !     write(101, *) elem%xc(:), vals(:)

                   if(je == 0 .and. l==1) then
                      elem%Ldeg = elem%deg
                      L_rule => state%space%L_rule(elem%Ldeg)
                      Ldof = L_rule%Qdof

                      allocate(elem%Lxi(1:Ldof, 1:nbDim+iRe-1 ) ) !3rd component: 1:2=coords, 3,4,.. .=coeffs
                      
                      call ComputeF(elem, Ldof, L_rule%lambda(1:Ldof, 1:2), elem%Lxi(1:Ldof, 1:nbDim) )

                      do k = 1,Ldof
                         call  Set_porous_media_Hraz(elem%Lxi(k,1), elem%Lxi(k,2), iRe, vals(1:iRe) )
                         elem%Lxi(k, 2+1: 2+iRe-1) = vals(1: iRe-1) 
                         
                         !write(103,*) Fxi(i, 1:2), vals(i, 1:2+iRe)
                      enddo
                   endif


                   
                elseif( state%model%idiff == 4) then
                   ! no action

                elseif( state%model%idiff == 13) then
                   call  Set_porous_media_A2(xi(2), -xi(1), val1, val2, val3)
                   elem%xi(jc, l, 3) = val1

                else
                   stop "No other implemented (3)"
                endif
             enddo
          enddo ! je = 0, flen
       enddo  ! i=1,gridA%nelem
       deallocate( vals)
    endif

    !stop ' subroutine Setting_of_space_variable_coeffs'

  end subroutine Setting_of_space_variable_coeffs


  !> compute the volume and edge integ nodes.
  subroutine ComputeIntegNode(elem)
    type(element), intent(inout) :: elem
    real, dimension(:, :), allocatable :: xi
    integer :: Fdof, Qdof, Qnum, Gnum
    integer je, jc, l
    real :: t

    ! setting of integ nodes

    Fdof = max(elem%Qdof, maxval( elem%face(fGdof, :)  ) )

    if ( allocated( elem%xi ) ) deallocate(elem%xi)

    allocate(elem%xi(0:elem%flen, 1:Fdof, 1:nbDim+iRe-1 ) ) !3rd component: 1:2=coords, 3,4,.. .=coeffs


    ! volume integ nodes
    Qdof = elem%Qdof
    Qnum = elem%Qnum
    call ComputeF(elem, Qdof,  state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim), &
         elem%xi(0, 1:Qdof,1:nbDim) )

    !do l=1,Qdof
    !   write(*, *) elem%i, state%space%V_rule(Qnum)%lambda(l,1:nbDim), &
    !     elem%xi(0, l,1:nbDim)
    !enddo


    ! edge integ nodes
    do je=1, elem%flen
       jc = je
       if(elem%HGnode) jc = elem%HGface(1, je)
       Gnum = elem%face(fGnum,je)

       allocate(xi(1:Gnum, 1:nbDim) )

       if(elem%type == 3) then ! triangle

          do l=1,Gnum
             ! Gauss integration node on (0,1)
             t = state%space%G_rule(Gnum)%lambda(l)
             if(elem%HGnode) t = ResizeHG(t, elem%HGface(2, je))

             ! xi(1:nbDim) ... integration node on reference element
             if(jc == 1) xi(l,1:nbDim) = (/ t, 0./)
             if(jc == 2) xi(l,1:nbDim) = (/ 1-t, t/)
             if(jc == 3) xi(l,1:nbDim) = (/ 0., 1- t/)
          enddo

       elseif(elem%type == 4) then ! quadrilateral

          do l=1,Gnum
             t = state%space%G_rule(Gnum)%lambda(l)
             if(elem%HGnode) t = ResizeHG(t, elem%HGface(2, je))

             ! xi(1:nbDim) ... barycentric coordinates of integration node
             if(jc == 1) xi(l,1:nbDim) = (/ t    , 0.   /)
             if(jc == 2) xi(l,1:nbDim) = (/ 1.   , t    /)
             if(jc == 3) xi(l,1:nbDim) = (/ 1.-t , 1.   /)
             if(jc == 4) xi(l,1:nbDim) = (/ 0.   , 1.-t  /)
          enddo

       else
          print*,'Only triang and quad in "PrepareOneElement" implemented (2)'
          stop
       endif

       call ComputeF(elem, Gnum, xi(1:Gnum, 1:nbDim), elem%xi(je, 1:Gnum, 1:nbDim) )
       deallocate( xi)
    enddo

    ! VOLUME INTEG NODES
    !do l=1, elem%Qdof
    !   write(80,*) elem%xi(0,l,1:2)
    !enddo

    ! EDGE INTEG NODES
    !do je=1,elem%flen
    !   do l=1, elem%face(fGdof, je)
    !      write(90,*) elem%xi(je,l,1:2)
    !   enddo
    !enddo
  end subroutine ComputeIntegNode

! moved to o_mesh.f90 setEdgeQuadratureDegrees
  !> setting of the degrees of edge quadratures
!  subroutine SetEdgeDegrees(grid)
!    use data_mod, g_grid => grid !mesh operation IS NOT USED for global grid (adaptation)
!    type(mesh), intent(inout) :: grid
!    integer :: i, j, ii
!
!    do i=1,grid%nelem
!       do j=1,grid%elem(i)%flen
!          if(grid%elem(i)%face(neigh,j) > 0) then
!             ii = grid%elem(i)%face(neigh,j)
!             ! maximal degree for face quadratures
!             grid%elem(i)%face(fdeg,j) = max(grid%elem(i)%deg, grid%elem(ii)%deg )
!             ! dof from opposite element for matrix shape
!             grid%elem(i)%face(fdof,j) = grid%elem(ii)%dof
!
!             !if(state%modelName == 'scalar' .or.state%modelName == '2eqs' .and. state%space%adapt%adapt_method == 'RTN') &   ! for dual problem
!             !     grid%elem(i)%face(fdof,j) = grid%elem(ii)%dof_plus
!             grid%elem(i)%face(fTdeg,j) = grid%elem(ii)%Tdeg
!             grid%elem(i)%face(fTdof,j) = grid%elem(ii)%Tdof
!
!             !NEW FR
!             grid%elem(i)%face(fdof_plus,j) = grid%elem(ii)%dof_plus
!
!!             if(state%SP) then
!!	        grid%elem(i)%face(fdegP,j) = max(grid%elem(i)%degP, grid%elem(ii)%degP )
!!                grid%elem(i)%face(fdofP,j) = grid%elem(ii)%dofP
!!             endif
!
!          else
!             grid%elem(i)%face(fdeg,j) = grid%elem(i)%deg
!             grid%elem(i)%face(fdof,j) = grid%elem(i)%dof
!!             if(state%SP) then
!!	        grid%elem(i)%face(fdegP,j) = grid%elem(i)%degP
!!                grid%elem(i)%face(fdofP,j) = grid%elem(i)%dofP
!!             endif
!
!             !if(state%modelName == 'scalar' .or.state%modelName == '2eqs' .and. state%space%adapt%adapt_method == 'RTN') &   ! for dual problem
!             !     grid%elem(i)%face(fdof,j) = grid%elem(i)%dof_plus
!          endif
!       enddo
!    enddo
!
!  end subroutine SetEdgeDegrees



  !> preparation of the problem
  !> should be used after InitProblem or after each change of degree of approximation
  !> but NOT the degree of boundary approximation (if so use InitProblem
  !> at first)
  !>
  !> setting of degree of polynomial approximation for each element,
  !> either by default or by the solution file "rsolfile"
  !>
  !> initialization of matrix elements shape, setting of volume quadrature
  !> evaluation of f_mapping at integ nodes,
  !> computation of local mass matrices, ....
  subroutine PrepareProblem(grid, rrsolfile)
    use data_mod, g_grid => grid !mesh operation IS NOT USED for global grid (adaptation)
    class(mesh), intent(inout), target :: grid
    character(len=*), intent(in) :: rrsolfile
    !FR rsolfile changed to rrsolfile due to ambiguous reference
    class(element), pointer :: elem
    integer :: i

    call grid%computeGeometry()

    call state%space%copyGeometry( grid%h , grid%domain_volume )
    write(debug,*) 'H and domain volume:', state%space%h , state%space%domain_volume


    ! initiation of $p$ vector (degree of polynomial approximation)
    if(state%type_IC == 0 .or. state%type_IC == -3) then
       ! initial condition from an input file
       if(state%type_IC ==  0 ) call ElementReadResults(rrsolfile)
       if(state%type_IC == -3 ) call ElementReadResults_dgm(  )

       if(state%time%disc_time == 'STDG' ) then
          do i=1,grid%nelem
             elem => grid%elem(i)

             elem%Tdeg = state%time%deg
             call elem%initElementDof()

             allocate( elem%wST(1:ndim, 1:elem%dof, 1:elem%Tdof ) )
             allocate( elem%wSTfinAD(1:ndim, 1:elem%dof) )
             allocate( elem%rhsST(1:ndim, 1:elem%dof_plus, 1:elem%Tdof_plus ) )

          enddo
          print*,'  PrepareProblem just in implementation for STDGM'
          !stop
       endif

    else
       ! default initial conditions
       ! same degree (=state%space%deg)  of polynomial approx on each element
       state%nsize =  ndim * (state%space%deg+1)*(state%space%deg+2)/2 * grid%nelem
       state%time%ttime = 0.
       state%time%iter = 0.
       state%time%tau_old = 0.

       do i=1,grid%nelem
          elem => grid%elem(i)
          elem%deg  = state%space%deg
          elem%Tdeg = state%time%deg

!          elem%degP  = state%space%degP

          ! testing of the "oscillating" p CHANGE
          !if(elem%i == 1) then
          !   print*,'##################'
          !   print*,'ATTENTION IN PrepareProblem'
          !   print*,'##################'
          !endif

          !if( mod(elem%i, 2) == 0 ) elem%deg = elem%deg + 1

          !if( mod(elem%i, 3) == 0 ) elem%deg = elem%deg + 1
          !if( mod(elem%i, 3) == 1 ) elem%deg = elem%deg + 2

          call elem%initElementDof()

          ! the following command was moved to subroutine InitElemDof( elem )
          !elem%Tdof = elem%Tdeg + 1

!          if( state%SP ) then
!             if (state%modelName == 'incNS') then
!
!                allocate(elem%wSP(wV1:wP,0:state%time%deg,1:max(elem%dof,elem%dofP)))
!
!             endif
!          else
          if(state%time%disc_time == 'STDG' ) then

             allocate(elem%w(0:1, 1:ndim*elem%dof) )  ! used for saving wST in specific time
             allocate(elem%wST(1:ndim, 1:elem%dof, 1:elem%Tdof ) )
             allocate( elem%wSTfinAD(1:ndim, 1:elem%dof) )
             !NEW
             allocate(elem%rhsST(1:ndim, 1:elem%dof_plus, 1:elem%Tdof_plus ))
          else

             if( (state%modelName == 'scalar' .or.state%modelName == '2eqs') &
                  .and. state%space%adapt%adapt_method == 'RTN') then  ! for dual problem
                ! size w(*, dof+1: dof_plus) is for the evaluation of the error in dual norm
                !allocate(elem%w(0:state%time%deg+1, &
                !     elem%dof + elem%dof_plus) )
                allocate(elem%w(0:state%time%deg+1, 1:ndim*elem%dof_plus) )
             else
                allocate(elem%w(0:state%time%deg+1, 1:ndim*elem%dof) )
             endif

          endif

       enddo
    endif


    state%space%max_dof = maxval(grid%elem(:)%dof)
    state%time%max_Tdof = maxval(grid%elem(:)%Tdof)

    call grid%setEdgeQuadratureDegrees()

    !prepare the Time part of MassMatrix and set the tQnum of time quadrature used
    associate ( time => state%time )
    select type (time)
      class is ( TimeTDG_t )

       if ( state%space%estim_space == 'RTNst' ) then
         ! F@R Should it differ for nonlinear problems ???
         ! + 1 is needed for the Radau reconstruction
         time%Qnum = time%max_Tdof + 1
         write(*,*) 'TEMPORARILY!!!!  time%Qnum set in problem to max_Tdof + 1!!!'
       else
         ! + 1 is needed only for the RES deg_plus
         time%Qnum = time%max_Tdof  + 1
         if(state%model%varying_time_term)  time%Qnum = time%max_Tdof  + 2
         write(*, '(a46,i3)') '# time%Qnum set in problem to max_Tdof + ? = ', time%Qnum
       endif

       !NEW for adaptation - Tdof_plus -> max_Tdof + 1
       allocate( time%refTimeMatrix, time%StiffTimeMatrix)
       call InitMblock( time%refTimeMatrix, time%max_Tdof + 1, time%max_Tdof + 1)
       call InitMblock( time%StiffTimeMatrix, time%max_Tdof + 1, time%max_Tdof + 1)

       call ComputeRefTimeMatrix(time%refTimeMatrix, time%StiffTimeMatrix, time%max_Tdof + 1)

    end select
    end associate


    ! initiation of element setting
    do i=1, grid%nelem
       elem => grid%elem(i)
       call SetElementQuadraturesDegrees(elem )

       call PrepareOneElement(elem )

       call ComputeLocalMassMatrix(elem )
       ! dof_plus = .false.
       call InitMatrixElementShape(elem, .false. )

    enddo

    state%space%max_Qdof = maxval(grid%elem(:)%Qdof)  ! usually fGdof <= Qdof, hence this is OK


    if(state%space%adapt%adapt_type == 'ANI' .or. state%space%adapt%adapt_type == 'Ahp' &
         .or. state%space%adapt%adapt_type == 'Ihp') &
         call EvalAMAparams(grid)

    ! assempling matrix block together
    call InitGlobalMatrixShape(grid)


    ! setting of boundary conditions
    call SetConstBC(grid)

    ! setting of the space dependent coefficients, if necessary
    call Setting_of_space_variable_coeffs(grid )

    ! setting of initial conditions, filling arays elem()%w(0,:)
    if(state%type_IC == 0) then
       ! initial condition from an input file done in InitStateVector

    else
       ! setting of initial state vector
       if(state%type_IC == -1) then
          ! reading of initial values from file 'resultsx' from ANGENER
          call SetElementsIC_AMA()

       elseif(state%type_IC == -2 .or. state%type_IC == -3) then
          ! reading of initial values from file 'dgm.sol' from ExplDGM
          call SetElementsIC_DGM()

       else
          ! setting of initial values from the BC
          call SetElementsIC()

       endif
    endif

    ! setting of the correct numbers of output siles
    if( (state%type_IC== 0 .or. state%type_IC== -1 .or. state%type_IC== -2) .and. state%time%OutTime > 0. ) &
       state%isol = int( state%time%ttime / state%time%OutTime)

    ! multigrid
    !    print*, 'Multigrid:' , state%MGsolver
    if( state%MGsolver ) call InitMG( )

    !print*,'# Problem prepared for computation'
  end subroutine PrepareProblem


SUBROUTINE InitMG()
  integer :: i

  grid%elem(:)%MGdof = grid%elem(:)%dof
  grid%elem(:)%MGdeg = grid%elem(:)%deg
  grid%elem(:)%MGncv = grid%elem(:)%ncv

  !do i=1,grid%nelem
  !   call InitMGMatrixElementShape(grid%elem(i))
  !end do

END SUBROUTINE InitMG

SUBROUTINE DeInitMG()
  integer :: i, j
  class(element), pointer :: elem
  grid%elem(:)%MGdof = grid%elem(:)%dof
  grid%elem(:)%MGdeg = grid%elem(:)%deg
  grid%elem(:)%MGncv = grid%elem(:)%ncv


  print*
  print*,' ANDREJ, repare DeInitMG, please :-) '
  print*
  return


  do i=1,grid%nelem
     elem => grid%elem(i)
     do j=0,elem%flen
        deallocate(elem%XXX(j)%Mb)
     enddo

     deallocate(elem%XXX)
  end do

END SUBROUTINE DeInitMG


  ! FIXME: dealocation is needed
  subroutine InitMGMatrixElementShape(elem)
    type(element),  intent(inout):: elem   ! elem = element
    integer :: dof, nsl, nsl1,  i,j, nslST, nslST1

    ! allocation of blocks
    dof = elem%MGdof
    nsl = dof * ndim

    ! blocks of the (main) sparse preconditioner matrix
    !allocate( elem%ILU( 0:0) )
    allocate( elem%XXX( 0:elem%flen) )

    call InitMblock(elem%XXX(0), nsl, nsl)
    do i=1,nsl
       elem%XXX(0)%Mb(i,i) = 1.
    enddo

    ! off diagonal blocks
    do j=1,elem%flen
       if ( elem%face(neigh,j) > 0 ) then
          nsl1 = elem%face(fdof,j) * ndim
          call InitMblock(elem%XXX(j), nsl, nsl1)
        endif
    enddo

  end subroutine InitMGMatrixElementShape


  !> allocates %block, %blockST, %ILU %vec
  !> would not work if elem%Tdof differs through the mesh (elem%Tdof_plus on edges)
  !> initialization of state vectors, matrix blocks for each element \f$ K\in{\cal T}_h\f$
  subroutine InitMatrixElementShape(elem, deg_plus)
    type(element),  intent(inout):: elem   ! elem = element
    logical, intent(in) :: deg_plus
    integer :: dof, dofP, dofP1, nsl, nsl1,  i,j, nslST, nslST1, nslPlus

!    if(state%SP) then ! saddle point problems
!
!     ! blocks of the (main) sparse matrix
!     allocate( elem%SPblock(1:Bdim , 0:elem%flen) )
!
!     if(state%modelName == 'incNS') then ! incompressible  Navier-Stokes equations
!        ! diagonal blocks
!        nsl = nbDim*dof   ! velocity
!        dofP = elem%dofP
!
!        call InitMblock(elem%SPblock(bVV,0), nsl, nsl)
!        call InitMblock(elem%SPblock(bVP,0), nsl, dofP)
!        call InitMblock(elem%SPblock(bPV,0), dofP, nsl)
!
!         ! off diagonal blocks
!          do j=1,elem%flen
!             if(elem%face(neigh,j) > 0) then
!                nsl1 = elem%face(fdof,j) * nbDim
!                dofP1 = elem%face(fdofP,j)
!
!                call InitMblock(elem%SPblock(bVV,j), nsl, nsl1)
!                call InitMblock(elem%SPblock(bVP,j), nsl, dofP1)
!                call InitMblock(elem%SPblock(bPV,j), dofP, nsl1)
!
!             endif
!          enddo
!
!     else
!        print*, 'UNKNOWN type', state%modelName, ' in InitMatrixElementShape(elem)'
!        stop
!     endif
!
!    else    ! not state%SP (standard convection-diffusion)


      if (state%time%disc_time /= 'STDG') then
               ! allocation of blocks
          if (deg_plus) then
            dof = elem%dof_plus
          else
            dof = elem%dof
          endif
          nsl = dof * ndim

          if(.not. state%init_only_matrix_blocks) then
             ! right hand sides
             !allocate( elem%vec(1:max_vecs, 1:nsl) )
             allocate( elem%vec(1:max_vecs, 1:elem%dof_plus * ndim) , source = 0.0 )  ! for aposteriori
          endif

          ! blocks of the (main) sparse matrix
          allocate( elem%block( 0:elem%flen) )

          ! blocks of the (main) sparse preconditioner matrix
          !allocate( elem%ILU( 0:0) )
          allocate( elem%ILU( 0:elem%flen) )

          ! diagonal blocks
          call InitMblock(elem%block(0), nsl, nsl)

          call InitMblock(elem%ILU(0), nsl, nsl)
          do i=1,nsl
             elem%ILU(0)%Mb(i,i) = 1.
          enddo

          ! off diagonal blocks
          do j=1,elem%flen
             if(elem%face(neigh,j) > 0) then
                if ( deg_plus ) then
                  nsl1 = elem%face(fdof_plus,j) * ndim
                else
                  nsl1 = elem%face(fdof,j) * ndim
                endif

                call InitMblock(elem%block(j), nsl, nsl1)
                call InitMblock(elem%ILU(j), nsl, nsl1)
             endif

          enddo

          ! ST DGM method inicialization
       else
          ! allocation of blocks
          if (deg_plus) then
            dof = elem%dof_plus
            nslST = dof * ndim * elem%Tdof_plus
          else
            dof = elem%dof
            nslST = dof * ndim * elem%Tdof
          endif
          nsl = dof * ndim

          ! vector
          allocate( elem%vec(1:res_vec, 1 : elem%dof_plus * ndim *elem%Tdof), source = 0.0 )

          ! blocks of the (main) sparse matrix
          allocate( elem%block(   0:elem%flen) )
          allocate( elem%blockST( 0:elem%flen) )
          allocate( elem%ILU(     0:elem%flen) )

          !diagonal blocks
          call InitMblock(elem%block(0), nsl, nsl)
          call InitMblock(elem%blockST(0), nslST, nslST)
          call InitMblock(elem%ILU(0), nslST, nslST)

          do i=1,nslST
             elem%ILU(0)%Mb(i,i) = 1.
          enddo

          ! off diagonal blocks
          do j=1,elem%flen
             if(elem%face(neigh,j) > 0) then
               if ( deg_plus ) then
                  nsl1 = elem%face(fdof_plus,j) * ndim
                  nslST1 = nsl1 * ( elem%Tdof_plus )
                else
                  nsl1 = elem%face(fdof,j) * ndim
                  nslST1 = nsl1 * elem%face(fTdof,j)
                endif

                call InitMblock(elem%block(j), nsl, nsl1)
                call InitMblock(elem%blockST(j), nslST, nslST1)
                call InitMblock(elem%ILU(j), nslST, nslST1)
             endif
          enddo


       endif !stdgm

!    endif ! state%SP   saddle-point problem

  end subroutine InitMatrixElementShape

  !> deallocations of state vectors, matrix blocks for each element \f$ K\in{\cal T}_h\f$
  subroutine DeleteMatrixElementShape(elem)
    type(element),  intent(inout):: elem   ! elem = element
    integer :: j

    if(.not. state%init_only_matrix_blocks) &
         deallocate( elem%vec )

    deallocate(elem%block(0)%Mb, elem%ILU(0)%Mb)
    do j=1,elem%flen
       if(elem%face(neigh,j) > 0) deallocate(elem%block(j)%Mb, elem%ILU(j)%Mb)
    enddo

    !!deallocate( elem%block( 0:elem%flen), elem%ILU( 0:elem%flen) )
  end subroutine DeleteMatrixElementShape

  !> store the degrees of frredom of the problem
  subroutine EvalAMAparams(gridL)
    type(mesh), intent(inout) :: gridL
    integer :: i

    i = max(0, state%space%adapt%adapt_level)
    state%space%adapt%AMAhistory(i, 1) = gridL%npoin
    state%space%adapt%AMAhistory(i, 2) = gridL%nelem
    state%space%adapt%AMAhistory(i, 3) = sum(gridL%elem(:)%dof*ndim)

    !!write(*,'(a8,i5,3i8)') '@@@@@@@@',i, state%space%adapt%AMAhistory(i, 1:3)

  end subroutine EvalAMAparams

  subroutine PlotSolPrimal( k )
      integer, intent(in) :: k ! add this number to the name of the file
      class(element), pointer :: elem
      character(len=50) :: plot_file_primal = 'plot_sol_primal_' !'../DWR/plot_sol_primal.gnu'
      integer :: iPlotPrimal = 41
      integer :: i

      if (k<10) then
         write(plot_file_primal, "(A18,I1,A4)" ) 'plot_sol_primal_00', k, '.gnu'
      elseif (k<100) then
         write(plot_file_primal, "(A17,I2,A4)" ) 'plot_sol_primal_0', k, '.gnu'
      else
         write(plot_file_primal, "(A16,I3,A4)" ) 'plot_sol_primal_', k, '.gnu'
      endif

      ! print primal solution
      open( iPlotPrimal, file = plot_file_primal, action="write", status="replace" )
      do i = 1, grid%nelem
         elem => grid%elem(i)
         call Transfer_wST_to_w_Elem( elem , 0, elem%TQnum)
         call PlotElemSolution3D(iPlotPrimal, elem)
       end do
       print*, 'Primal solution saved to ',  plot_file_primal !, '.gnu'
       close(iPlotPrimal)

  end subroutine PlotSolPrimal

  subroutine PlotSolDual( k )
      integer, intent(in) :: k ! add this number to the name of the file
      class(element), pointer :: elem
      character(len=50) :: plot_file_dual = 'plot_sol_dual_' !'../DWR/plot_sol_primal.gnu'
      integer :: iPlot = 41
      integer :: i

      if (k<10) then
         write(plot_file_dual, "(A16,I1,A4)" ) 'plot_sol_dual_00', k, '.gnu'
      elseif (k<100) then
         write(plot_file_dual, "(A15,I2,A4)" ) 'plot_sol_dual_0', k, '.gnu'
      else
         write(plot_file_dual, "(A14,I3,A4)" ) 'plot_sol_dual_', k, '.gnu'
      endif

      ! print primal solution
      open( iPlot, file = plot_file_dual, action="write", status="replace" )
      do i = 1, grid%nelem
         elem => grid%elem(i)

         elem%w(0,:) = Transfer_funST_to_fun( elem%zST(1:ndim, 1:elem%dof, 1:elem%Tdof) , &
               elem%dof, elem%Tdof, 0, elem%TQnum )
         call PlotElemSolution3D(iPlot, elem)
       end do
       print*, 'Dual solution saved to ', plot_file_dual !,'.gnu'
       close(iPlot)

  end subroutine PlotSolDual



  !> initiate Global Matrix Shape, neccesary to carried out after each h of p adaptation
  subroutine InitGlobalMatrixShape(grid)
    use data_mod, grid_g => grid
    type(mesh), intent(inout) :: grid
    class(element), pointer :: elem   ! elem = element
    integer :: i, j, nsize, nsizeP
!    logical :: deg_plus
    integer :: loc_nsize, loc_nonzero

!    if ( size( grid%elem(1)%blockST(0)%Mb(1,:)) == &
!         grid%elem(1)%dof * grid%elem(1)%Tdof * ndim ) then
!       deg_plus = .false.
!       print*, 'InitGlobalMatrixShape with deg_plus = FALSE'
!    else if ( size( grid%elem(1)%blockST(0)%Mb(1,:)) == &
!         grid%elem(1)%dof_plus * grid%elem(1)%Tdof_plus * ndim ) then
!       deg_plus = .true.
!       print*, 'InitGlobalMatrixShape with deg_plus = TRUE'
!    else
!      stop 'wrong size of elem%blockST in InitGlobalMatrixShape'
!    endif


    !    InitGlobalMatrixShape always with deg_plus = false nsize - smaller

    if (state%time%disc_time == 'STDG') then

      state%nonzero = 0
      nsize = 0

      state%nsize = sum(grid%elem(:)%dof*ndim*grid%elem(:)%Tdof)

      !print*, 'Problem size is: (InitGlobalMatrixShape)', state%nsize

      do i=1,grid%nelem
         elem => grid%elem(i)
         ! necessary for serial block matrix-element multiplication
         elem%ncv = nsize + 1     ! number within the sequence of element's dof

         call GetBlocksSizeElemST( elem, loc_nsize, loc_nonzero, .false. )
         nsize = nsize + loc_nsize
         state%nonzero  = state%nonzero + loc_nonzero
      enddo

      if ( state%nsize /= nsize ) then
         print*, 'problem in InitGlobalMatrixShape: state%nsize /= nsize :' , &
            state%nsize, nsize
         stop
      endif

    else
       state%nsize = sum(grid%elem(:)%dof*ndim)

       state%nonzero = 0
       !open(14,file='matrix_shape',status='unknown')
       nsize = 0

       do i=1,grid%nelem
          elem => grid%elem(i)

          ! necessary for serial block matrix-element multiplication
          elem%ncv = nsize + 1     ! number within the sequence of element's dof
          call GetBlocksSizeElem( elem, loc_nsize, loc_nonzero, .false. )
          nsize = nsize + loc_nsize
          state%nonzero  = state%nonzero + loc_nonzero

       enddo

    endif ! stdgm

    !close(14)
    !if(state%time%iter_loc <= 0) &
    if(.not. state%local_problem ) then
!       if (state%SP) then
!          write(*, '(a68, i7, a2, i8, a2, i8, a2, i10, a2, i10, a2, i10, a3, f6.3, a20, f6.3, a2)') &
!               '# Problem size (nelem, dof, dofP, nonzero, nonzeroVP, nonzeroPV):', grid%nelem, ', ', &
!               state%nsize,', ', state%nsizeP, ', ', state%nonzero, ', ', state%nonzeroVP, ', ', state%nonzeroPV,' (=',&
!               1.*state%nonzero/state%nsize /state%nsize*100.,'% and for VP blocks ', &
!               1.*state%nonzeroVP/state%nsize /state%nsizeP*100.,'%)'
!          ! (1.*state%nonzero/state%nsize) *100.,'%)'
!       else
          write(*, '(a38, i7, a2, i8, a2, i10, a3, f6.3,a2)') &
               '# Problem size (nelem, dof, nonzero):', grid%nelem, ', ', &
               state%nsize,', ', state%nonzero,' (=',&
               1.*state%nonzero/state%nsize /state%nsize*100.,'%)'
          ! (1.*state%nonzero/state%nsize) *100.,'%)'
!       endif SP
       !!print*,'#######',state%time%iter_loc
    endif

  end subroutine InitGlobalMatrixShape


  !> computing of the mass matrix, stiff matrix and the inverse of mass matrix
  !> for each element
  !>
  !> mass: \f$ M_{R,C} = \int_{elem} \phi_R\, \phi_C\ dx,\quad \f$
  !> stiff: \f$ S_{R,C} = \int_{elem} \nabla\phi_R\cdot \nabla\phi_C\ dx,\quad \f$
  !> inverse mass: \f$ M^{-1}\f$
  subroutine ComputeLocalMassMatrix(elem)
    type(element), intent (inout) :: elem
    real, dimension(:), allocatable :: ident
    integer :: dof, dof_plus, Qdof, Qnum


    dof = elem%dof
    !dof_plus for aposteriori estimates
    dof_plus = elem%dof_plus

    Qdof = elem%Qdof

    Qnum = state%space%Qdeg(elem%deg+1,1)

    ! identity vector
    allocate(ident(1:Qdof) )
    ident(1:Qdof) = 1.

    ! local mas matrix
    allocate(elem%Mass)
    call InitMblock(elem%Mass,dof_plus,dof_plus)
    call IntegrateBlockBBmass(elem, Qnum, dof_plus, dof_plus , elem%Mass%Mb(1:dof_plus,1:dof_plus) )
    ! normalization
    !elem%Mass%Mb(:,:) = elem%Mass%Mb(:,:) /elem%area

    ! inverse of local mass matrix
    allocate(elem%MassInv)
    call InitMblock(elem%MassInv,dof,dof)

    elem%MassInv%Mb(1:dof, 1:dof) = elem%Mass%Mb(1:dof, 1:dof)
    call MblockInverse(dof, elem%MassInv%Mb)

    ! stiff matrix
    allocate(elem%Stiff)
    call InitMblock(elem%Stiff,dof,dof)
    call IntegrateBlockD2(elem, elem%dof, ident(1:Qdof), elem%Stiff%Mb(1:dof,1:dof) )

    deallocate(ident)

 end subroutine ComputeLocalMassMatrix


 subroutine ComputeRefTimeMatrix( A , B, Tdof)
   type(Mblock), intent(inout) :: A, B
   integer, intent(in) :: Tdof
   integer :: i, j, l, TQnum
   class(Time_rule), pointer :: T_rule
   real :: val, val1, val2

   !print*,'   TQnum = Tdof ! number of quadrature nodes '
   TQnum = Tdof ! number of quadrature nodes

   A%Mb(:,:) = 0.

   T_rule => state%time%T_rule(TQnum)

      !print*, 'TQNUm: ', TQnum , Tdof

      do i = 1,Tdof

         do j = 1, Tdof
            !	print*, 'dphi = T_rule%Dphi(j,1:TQnum)'
            A%Mb(i,j) = dot_product(T_rule%weights(1:TQnum), T_rule%Dphi(j,1:TQnum) *T_rule%phi(i,1:TQnum))&
                 + T_rule%phi(i,-1)*T_rule%phi(j,-1)

            !write(*,'(a8,2i5, 30es12.4)') 'refTim',i,j,  A%Mb(i,j), &
            !     dot_product(T_rule%weights(1:TQnum), T_rule%Dphi(j,1:TQnum) *T_rule%phi(i,1:TQnum)), &
            !      T_rule%phi(i,-1)*T_rule%phi(j,-1)

            !print*, 't-phi', T_rule%phi(i,-1)*T_rule%phi(j,-1)

            B%Mb(i,j) = dot_product( T_rule%weights(1:TQnum), &
                 !T_rule%phi(j,1:TQnum) * T_rule%phi(i,1:TQnum) &
                  T_rule%Dphi(j,1:TQnum) * T_rule%Dphi(i,1:TQnum) )

            ! TEST of orthogonality
            !val = dot_product(T_rule%weights(1:TQnum) , T_rule%phi(j,1:TQnum) * T_rule%phi(i,1:TQnum) )
            !val1 = dot_product(T_rule%weights(1:TQnum) , T_rule%Dphi(j,1:TQnum) * T_rule%phi(i,1:TQnum) )
            !val2 = dot_product(T_rule%weights(1:TQnum) , T_rule%Dphi(j,1:TQnum) * T_rule%Dphi(i,1:TQnum) )
            !write(*,'(a6,2i5,20es16.8)') 'ORT', i,j, val, val1, val2, B%Mb(i,j), A%Mb(i,j), T_rule%phi(i,-1)
         enddo
         !write(*,'(a6,i5,20es16.8)') 'ORT', i,  B%Mb(i,:)
      enddo

      !print*, 'Still running'
      !stop

 end subroutine ComputeRefTimeMatrix



 !> test the accuracy of quadrature rules on monomials
 subroutine Test_V_rules( )
   integer :: i, i1, l
   integer :: a,b, deg
   real :: area, coef, quad, x, y, val, exact, err
   logical :: line


   ! negative weights or nodes outsides of element
   do i1 = 0, MaxDegreeImplemented
      i = state%space%Qdeg(i1, 1)

      do l = 1,state%space%V_rule(i)%Qdof

         if(minval(state%space%V_rule(i)%lambda(l,:)) < 0. ) &
              write(*,'(a6,2i5,3es12.4)') 'negQUA',i,l,state%space%V_rule(i)%lambda(l,:)
         if(maxval(state%space%V_rule(i)%lambda(l,:)) > 1. ) &
              write(*,'(a6,2i5,3es12.4)') 'negQUA',i,l,state%space%V_rule(i)%lambda(l,:)
         if(state%space%V_rule(i)%weights(l) < 0. ) &
              write(*,'(a6,2i5,3es12.4)') 'negQUA',i,l,state%space%V_rule(i)%weights(l)
      enddo
   enddo
   print*,'######## positivity of quadratures checked'


   area = 0.5D+00

   !do deg = 34, 36
   do deg = 0, 36
      do b = 0, deg

         a = deg - b
         write ( *, '(a6,i3,a1,i2,a1,i2,a28)') &
              'deg=', a+b,'(',a,'+',b,') --------------------------------'
         !
         !  Multiplying X**A * Y**B by COEF will give us an integrand
         !  whose integral is exactly 1.  This makes the error calculations easy.
         !
         coef = real ( a + b + 2, kind = 8 ) * real ( a + b + 1, kind = 8 )
         do i = 1, b
            coef = coef * real ( a + i, kind = 8 ) / real ( i, kind = 8 )
         end do

         line = .false.
         !do i = 1,2
         do i1 = 1, MaxDegreeImplemented
            i = state%space%Qdeg(i1, 1)

            quad = 0.0D+00

            do l = 1,state%space%V_rule(i)%Qdof

               x = state%space%V_rule(i)%lambda(l,1)
               y = state%space%V_rule(i)%lambda(l,2)

               if ( a == 0 .and. b == 0 ) then
                  val = coef
               else if ( a == 0 .and. b /= 0 ) then
                  val = coef        * y**b
               else if ( a /= 0 .and. b == 0 ) then
                  val = coef * x**a
               else if ( a /= 0 .and. b /= 0 ) then
                  val = coef * x**a * y**b
               end if

               quad = quad + state%space%V_rule(i)%weights(l) * val

            end do

            quad = area * quad

            exact = 1.0D+00
            err = abs ( exact - quad )

            !if(i >= 20) then
               if(abs(err) < 1E-10) then
                  write ( *, '(a6,i3,a1,i2,a1,i2,a8,i8,2x,es18.10,2x,es14.6)' ) &
                       'deg=', a+b,'(',a,'+',b,') Qnum=', i, quad, err
                  line = .true.
               else
                  write ( *, '(a6,i3,a1,i2,a1,i2,a8,i8,2x,es18.10,2x,es14.6, a6)' ) &
                       'deg=', a+b,'(',a,'+',b,') Qnum=', i, quad, err,'!!!!!!'
               endif
            !endif



            !deallocate ( wtab )
            !deallocate ( xytab )

      end do
      !if(line) write ( *,*)
    end do

  end do


 end subroutine Test_V_rules

 !> computes nsize and nzero for one element
 subroutine GetBlocksSizeElemST( elem, nsize, nonzero, deg_plus )
   class( element ), intent(in) :: elem
   integer, intent(out) :: nsize
   integer, intent(out) :: nonzero
   logical, intent(in) :: deg_plus
   integer :: j

   if (deg_plus) then
      nsize = elem%dof_plus*elem%Tdof_plus*ndim
   else
      nsize = elem%dof*elem%Tdof*ndim
   endif
   nonzero = size( elem%blockST(0)%Mb(1,:)) * size(elem%blockST(0)%Mb(:,1))

   do j=1,elem%flen
         if(elem%face(neigh,j) > 0) &
            nonzero  = nonzero &
            + size(elem%blockST(j)%Mb(1,:)) * size(elem%blockST(j)%Mb(:,1))
   enddo

 end subroutine GetBlocksSizeElemST

  !> computes nsize and nzero for one element
 subroutine GetBlocksSizeElem( elem, nsize, nonzero, deg_plus )
   class( element ), intent(in) :: elem
   integer, intent(out) :: nsize
   integer, intent(out) :: nonzero
   logical, intent(in) :: deg_plus
   integer :: j

   if (deg_plus) then
      nsize = elem%dof_plus * ndim
   else
      nsize = elem%dof * ndim
   endif

   nonzero = size(elem%block(0)%Mb(1,:)) * size(elem%block(0)%Mb(:,1))

   do j=1,elem%flen
      if(elem%face(neigh,j) > 0) &
         nonzero  = nonzero &
            + size(elem%block(j)%Mb(1,:)) * size(elem%block(j)%Mb(:,1))
   enddo

 end subroutine GetBlocksSizeElem


   !> update all needed array sizes and matrix sizes when the %elem deg was changed
   !> called after call grid%setNewElementDofs( p_plus )
   !> the solution vectors are not filled , it is need to call ComputeSTDGM_Terms with implicitly = .true. and .false. after this subroutine
   subroutine UpdateMesh_p( grid)
      class( mesh ), intent(inout) :: grid
      class( element ), pointer :: elem
      integer :: i, nelem, dof, Tdof
      real, pointer, dimension(:,:,:) :: temp
      integer :: dd, min_dim

      nelem = grid%nelem
      !print*, 'We do not deallocate elem%rhsST in UpdateMesh_p when reducing the system!'
      !print*, 'CONTROL whether it doesnt break the computation!'

   do i = 1, nelem
      elem => grid%elem(i)
      dof = elem%dof
      Tdof = elem%Tdof

      call SetElementQuadraturesDegrees(elem )

      ! parts of prepareOneElement:
      ! jacobian in integration nodes for nonlinear F
      if ( .not. elem%F%iFlin ) then
         ! probably leads to problems with allocation
         call ComputeJacobian_curved( elem )
      end if

      call ComputeIntegNode(elem )

      ! solution w, wST, rhsST, wSTfinAD
      ! deallocate
      temp => elem%wST

      ! prolong/shrink the vector wST , copy data
      nullify( elem%wST )
      allocate( elem%wST(1:ndim,1:dof,1:Tdof), source = 0.0 )
      dd = size( temp(1,:,1) )
      min_dim = min( dd, dof )
      if (dd > dof) then
         if (temp(1,min_dim + 1,1) /= 0 ) &
              stop 'some data is lost in UpdateMesh_p'
      endif

      elem%wST(1:ndim,1:min_dim ,1:Tdof) =  temp(1:ndim,1:min_dim ,1:Tdof)
      deallocate( temp )

      deallocate(elem%w)
      deallocate(elem%vec)
      ! DO NOT deallocate elem%wST, it is used!
      deallocate( elem%wSTfinAD )


      ! we can still use the residual
      ! WE NEVER REDUCE THE SIZE OF rhsST !!! NEEDED FOR PU-DWR ESTIMATES
      if ( size(elem%rhsST(:,:,:)) < ndim*elem%dof_plus*elem%Tdof_plus) then
       deallocate( elem%rhsST )
       ! rhsST is used in ComputeSTDGM_Terms for deg_plus
       allocate( elem%rhsST(1:ndim,1:elem%dof_plus,1:elem%Tdof_plus), source = 0.0 )
      else
         ! DO NOTHING
      endif
   !          call elem%deallocateSolutionVectors()

      ! allocate
      allocate( elem%w(0:1, 1:ndim*dof), source = 0.0 )
      allocate( elem%wSTfinAD(1:ndim, 1:dof), source = 0.0 )

      ! deallocate all matrix blocks:
      ! Mass, MassInv,Stiff,block, ILU, blockST
      call elem%deallocateMatrixBlocks( )
      ! recompute - needed only if p_plus > 1
      call ComputeLocalMassMatrix(elem )
      ! allocates also elem%vec - to dof_plus ?? PROBLEM ???
      call InitMatrixElementShape(elem, .false.)
   end do

   state%space%max_Qdof = maxval( grid%elem(:)%Qdof )
   ! assempling matrix block together
   call InitGlobalMatrixShape(grid)
   ! setting of boundary conditions - needed because integ point have changed
   call SetConstBC( grid )

 end subroutine UpdateMesh_p

! ! Plot solution wST in a given time which is number in [0,1] relative to the actual time step
! subroutine PlotSolST( relTime )
!      integer, intent(in) :: relTime ! add this number to the name of the file
!      class(element), pointer :: elem
!      character(len=50) :: plot_file = 'plot_sol_intime' !'../DWR/plot_sol_primal.gnu'
!      integer :: iPlotPrimal = 41
!      integer :: i
!
!!      write(plot_file, "(A16,I1,A4)" ) 'plot_sol_', , '.gnu'
!
!      ! print primal solution
!      open( iPlotPrimal, file = plot_file, action="write", status="replace" )
!      do i = 1, grid%nelem
!         elem => grid%elem(i)
!         call Transfer_wST_to_w_Elem( elem , 0, elem%TQnum)
!         call PlotElemSolution3D(iPlotPrimal, elem)
!       end do
!       print*, 'Primal solution saved to ',  plot_file_primal !, '.gnu'
!       close(iPlotPrimal)
!
!  end subroutine PlotSolST


end module problem_oper
