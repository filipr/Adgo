!> compute the computational error in different norms if the exact solution is known
module error_subs
  use mesh_oper
  use main_data
  use eval_sol
  use model_oper
  use set_solution
  use stdgm_mod
  use tquadrature_mod

  implicit none

  public:: ComputeH1NormElement
  public:: ComputeL2H1Error
  public:: ComputeL2H1ErrorElement
  public:: ComputeL2H1ErrorST
  public:: ComputeL2H1ErrorElementST
  public:: ComputeL2H1ErrorOLD
  public:: ComputeL2H1ErrorElementOLD
  public:: ComputeH1ErrorTimeDer
  public:: ComputeTotDiscAlgErrorElement
  public:: ComputeDiscAlgError
  public:: H1seminormOfErrorInTimeDerElem
  public:: SpaceTimeErrors
  public:: SpaceTimeErrorsElem
  public:: SpaceTimeErrorsElemST


  public:: ComputeElementConvErrors
  public:: ComputeConvErrors
  public:: ComputeConvErrorsTDG
  public:: ComputeEoFC

contains

 !> compute the \f$ L^2(\Omega) \f$-norm of
 !>  the difference of \f$ {\bf w}_h\f$
 !> and the exact solution
 subroutine ComputeL2H1Error(errL2, errH1, errH1_discrete, normL2, normH1)
   real, intent (inout) :: errL2, errH1, errH1_discrete, normL2, normH1
   class(element), pointer :: elem
   real :: err_locL2, norm_locL2, err_locH1, norm_locH1, errH1_dis_loc
   integer :: i

   errL2 = 0.
   normL2 = 0.
   errH1 = 0.
   normH1 = 0.
   errH1_discrete = 0.

   do i=1,grid%nelem
      elem => grid%elem(i)
      call  ComputeL2H1ErrorElement(elem, err_locL2, err_locH1, errH1_dis_loc, norm_locL2, norm_locH1)
      !call ComputeL2H1ErrorElementOLD(elem, err_locL2, err_locH1, norm_locL2, norm_locH1)

      errL2 = errL2 + err_locL2
      normL2 = normL2 + norm_locL2
      errH1 = errH1 + err_locH1
      normH1 = normH1 + norm_locH1
      errH1_discrete = errH1_discrete +  errH1_dis_loc

      elem%errL2 = err_locL2**0.5
      elem%errH1 = err_locH1**0.5

      if(state%time%iter == 0) elem%eta(IC,1) = elem%errL2

     ! if(elem%i > 125)
     ! write(*,'(a6,i5,12es12.4)') ':::',elem%i,   err_locH1, errH1 , err_locH1**0.5, errH1**0.5

      !write(550+state%space%adapt%adapt_level ,'(i5,12es12.4)'),elem%i, elem%errL2, errL2**0.5
      !if(elem%errL2 > 1E-3 .and. state%time%iter_loc == 2) &
      !     write(525 ,'(a6,2i5,12es12.4)'),'######', state%space%adapt%adapt_level, elem%i, &
      !     elem%errL2, errL2**0.5, elem%xc(:)

      !if(elem%i == 68) then
      !   write(515,'(a4,i5,80es12.4)') 'er!1',elem%i,err_locL2, errL2
      !   write(515,'(a4,i5,80es12.4)') 'er!1',elem%i,err_locL2**0.5, errL2**0.5
      !   write(515,*)'-------------------------------------------------'
      !endif
      !write(*,'(a4,i5,80es12.4)') 'er!1',elem%i,err_locL2, err_locH1, errH1_dis_loc
   enddo

  ! write(*,'(a4,i8,80es12.4)') 'er!1',grid%nelem, errL2**0.5, errH1**0.5,  err_locL2**0.5, &
  !      err_locH1**0.5,  (err_locL2 * grid%nelem)**0.5,  (err_locH1 * grid%nelem)**0.5


   errL2 = errL2**0.5
   errH1 = errH1**0.5
   normL2 = normL2**0.5
   normH1 = normH1**0.5
   errH1_discrete= errH1_discrete**0.5

 end subroutine ComputeL2H1Error

 !> compute the square of the \f$ L^2(K) \f$-norm of the difference of \f$ {\bf w}_h\f$
 !> and the exact solution, \f$ K \f$ = elem
 subroutine ComputeL2H1ErrorElement(elem, errL2, errH1, errH1_discrete, normL2, normH1)
   type(element), intent(inout):: elem   ! elem = element
   real, intent (inout) :: errL2, normL2, errH1, normH1, errH1_discrete
   real, dimension(:), allocatable :: wrc1, wrc2   ! w recomputed  in integ nodes
   real, dimension(:,:), allocatable :: wExact     ! exact solution  in integ nodes
   real, dimension(:,:,:), allocatable :: DwExact  ! Der of exact solution in integ nds
   real, dimension(:,:), allocatable :: wi, Fx      ! Der of exact solution in integ nds
   real, dimension(:,:,:), allocatable :: Dwi      ! Der of exact solution in integ nds
   integer :: j, k, l, dof, Qnum, Qdof
   real :: val

   errL2 = 0.
   errH1 = 0.
   normL2 = 0.
   normH1 = 0.
   errH1_discrete = 0.

   dof = elem%dof
   Qnum = elem%Qnum
   Qdof = elem%Qdof

   allocate(wrc1(1:Qdof), wrc2(1:Qdof) )
   allocate(wExact(1:Qdof, 1:ndim) )
   allocate(DwExact(1:Qdof, 1:ndim, 1:nbDim))
   allocate(wi(1:Qdof, 1:ndim))
   allocate(Dwi(1:Qdof, 1:ndim, 1:nbDim))


   !allocate(Fx(1:Qdof, 1:2) )

   !if(elem%i == 1) write(*,'(a8, i2, 2es14.6)') '??errs',j,state%time%ctime, state%time%ttime


   ! setting of the exact solution in integ nodes

   call SetExactSolutionQnodes(elem, state%space%V_rule(Qnum), &

        wExact(1:Qdof, 1:ndim), DwExact(1:Qdof, 1:ndim, 1:nbDim))

   !call ComputeF(elem, Qdof, state%space%V_rule(elem%Qnum)%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )


   ! setting of the approximate solution in integ nodes
   call Eval_w_Elem(elem, wi(1:Qdof, 1:ndim) )
   call Eval_Dw_Elem(elem, Dwi(1:Qdof, 1:ndim, 1:nbDim) )

   elem%errL8 = 0.

   ! setting of w in integration nodes
   !do k=1,1              ! k = index of component of w
   do k=1,ndim              ! k = index of component of w

      ! L2-norm
      ! setting of the exact solution in integ nodes: array wrc2
      wrc1(1:Qdof) = (wExact(1:Qdof, k) - wi(1:Qdof, k) )**2
      wrc2(1:Qdof) = (wExact(1:Qdof, k) )**2

      ! L^{\infty}-norm
      elem%errL8 = max( elem%errL8, maxval( wrc1(1:Qdof) )**0.5)

      call IntegrateFunction(elem, wrc1(1:Qdof), val  )
      errL2 = errL2 + val

!      if(elem%i == 1) &
!           write(*,'(a4,8es14.6)')'-ER--',state%time%ctime, val, wrc1(1:6)

      !do l=1,Qdof
      !   write(80+state%time%iter,*) Fx(l, 1:2), wExact(l, :), wi(l, :), abs(wExact(l, :)- wi(l, :))
      !enddo

      !   write(*,'(a4,i5,80es12.4)') 'er!i',elem%i, wi(:,1)
      !   write(*,'(a4,i5,80es12.4)') 'er!E',elem%i, wExact(:,1)
      !   write(*,'(a4,i5,80es12.4)') 'er!e',elem%i, wrc1(:)
      !!   write(515,*) '-------------------------------------------------'
      !endif

      call IntegrateFunction(elem, wrc2(1:Qdof), val  )
      normL2 = normL2 + val

      ! H1-seminorm
      ! setting of the exact solution in integ nodes

      do j=1,Qdof
         wrc1(j) = dot_product( Dwi(j, k, 1:nbDim)- DwExact(j, k, 1:nbDim), &
              Dwi(j, k, 1:nbDim) - DwExact(j, k, 1:nbDim) )
         wrc2(j) = dot_product( DwExact(j, k, 1:nbDim), DwExact(j, k, 1:nbDim))
      enddo

      call IntegrateFunction(elem, wrc1(1:Qdof), val  )
      errH1 = errH1 + val
      !print*,'###',errH1, val

      call IntegrateFunction(elem, wrc2(1:Qdof), val  )
      normH1 = normH1 + val

      ! discrete H1-seminorm with the  lifting operator for the dicrete gradient
      ! in pNeu for SIPG and NIPG [Ern, Vohralik, SINUM 15]
      if( state%space%estim_space == 'pNeu')  then
         do j=1,Qdof
            Dwi(j, k, 1:nbDim) = Dwi(j, k, 1:nbDim) + state%space%m_IPG*elem%lifting(k, 1:nbDim) !opposite sign!

            wrc1(j) = dot_product( Dwi(j, k, 1:nbDim)- DwExact(j, k, 1:nbDim), &
                 Dwi(j, k, 1:nbDim) - DwExact(j, k, 1:nbDim) )
         enddo
         call IntegrateFunction(elem, wrc1(1:Qdof), val  )
         errH1_discrete = errH1_discrete + val
      endif
   end do

   !if(elem%i <= 4 .or. elem%i >= grid%nelem - 1) &
   !     write(*,'(a8,6es12.4)') 'Erro:',elem%errL8, errL2**0.5, errH1**0.5, errH1_discrete  !, normL2, normH1

   deallocate(wrc1, wrc2, wExact, DwExact, wi, Dwi )

 end subroutine ComputeL2H1ErrorElement


  !> compute the \f$ H^1(\Omega) \f$-seminorm of the discretization error and algebraic error
  subroutine ComputeDiscAlgError(DiscErr, AlgErr)
    real, intent (inout) :: DiscErr, AlgErr
    class(element), pointer :: elem
    real :: err_locDisc, err_locAlg
    integer :: i

    DiscErr = 0.
    AlgErr = 0.

    do i = 1, grid%nelem
       elem => grid%elem(i)
       call ComputeTotDiscAlgErrorElement(elem)

       DiscErr = DiscErr + elem%errDisc**2
       AlgErr = AlgErr + elem%errAlg**2

    enddo

    DiscErr = DiscErr**0.5
    AlgErr = AlgErr**0.5

  end subroutine ComputeDiscAlgError


  !> compute the \f$ H^1(K) \f$-seminorm of the total error, discretization error, and algebraic error, \f$ K \f$ = elem
  subroutine ComputeTotDiscAlgErrorElement(elem)
    type(element), intent(inout):: elem   ! elem = element
    real, dimension(:), allocatable :: wrc1, wrc2, wrc3   ! Tot, Disc, and Alg errors in integ nodes
    real, dimension(:,:), allocatable :: wExact     ! exact solution  in integ nodes
    real, dimension(:,:,:), allocatable :: DwExact  ! Der of exact solution in integ nds
    real, dimension(:,:,:), allocatable :: Dwi      ! Der of solution when GMRES converged in integ nds
                                                        ! - ALG in w(0,:), ALG2 in wc(0,:)
    real, dimension(:,:,:), allocatable :: Dwh      ! Der of computational solution based on AEE ALG in integ nds
                                                        ! - ALG in wc(0,:), ALG2 in w(1,:)
    integer :: k, dof, Qnum, Qdof, j

    dof = elem%dof
    Qnum = elem%Qnum
    Qdof = elem%Qdof

    allocate(wrc1(1:Qdof), wrc2(1:Qdof), wrc3(1:Qdof) )
    allocate(wExact(1:Qdof, 1:ndim) )
    allocate(DwExact(1:Qdof, 1:ndim, 1:nbDim) )
    allocate(Dwi(1:Qdof, 1:ndim, 1:nbDim) )
    allocate(Dwh(1:Qdof, 1:ndim, 1:nbDim) )


    ! setting of the eact solution in integ nodes
    call SetExactSolutionQnodes(elem, state%space%V_rule(elem%Qnum), &
         wExact(1:Qdof, 1:ndim), DwExact(1:Qdof, 1:ndim, 1:nbDim))

    if (state%space%adapt%adapt_method == 'ALG') then
      call Eval_Dw_Elem(elem, Dwi(1:Qdof, 1:ndim, 1:nbDim) )  ! w(0,:) = u_h

      call Eval_DwAEE_Elem(elem, Dwh(1:Qdof, 1:ndim, 1:nbDim) )  ! wc(0,:) = u_h^i

      if (elem%i < 5) then
         !print*, '******D u_h po spocteni adapt cyklus 0 pri max_adapt 0*********'
         !print*, 'Dwi(1:Qdof, 1, 1)', Dwi(1:Qdof, 1, 1)
         !print*, 'Dwi(1:Qdof, 1, 2)', Dwi(1:Qdof, 1, 2)
         !print*, ' '
         !print*, '******D u_h^i po spocteni adapt cyklus 0 pri max_adapt 0*********'
         !print*, 'Dwh(1:Qdof, 1, 1)', Dwh(1:Qdof, 1, 1)
         !print*, 'Dwh(1:Qdof, 1, 2)', Dwh(1:Qdof, 1, 2)
         !pause
      endif

    elseif (state%space%adapt%adapt_method == 'ALG2') then
       call Eval_DwAEE_Elem(elem, Dwi(1:Qdof, 1:ndim, 1:nbDim) )  ! wc(0,:) = u_h

       call Eval_Dw_Elem_time(elem, Dwh(1:Qdof, 1:ndim, 1:nbDim) )  ! w(1,:) = u_h^i

       if (elem%i < 5) then
         !print*, '******D u_h po spocteni adapt cyklus 0 pri max_adapt 0*********'
         !print*, 'Dwi(1:Qdof, 1, 1)', Dwi(1:Qdof, 1, 1)
         !print*, 'Dwi(1:Qdof, 1, 2)', Dwi(1:Qdof, 1, 2)
         !print*, ' '
         !print*, '******D u_h^i po spocteni adapt cyklus 0 pri max_adapt 0*********'
         !print*, 'Dwh(1:Qdof, 1, 1)', Dwh(1:Qdof, 1, 1)
         !print*, 'Dwh(1:Qdof, 1, 2)', Dwh(1:Qdof, 1, 2)
         !pause
      endif

    else
       print*, 'Sub. ComputeTotDiscAlgErrorElement is used only with adapt_method ALG and ALG2.'
       stop
    endif

    do k=1,ndim

       do j=1,Qdof
          wrc1(j) = dot_product( Dwh(j, k, 1:nbDim)- DwExact(j, k, 1:nbDim), &
               Dwh(j, k, 1:nbDim) - DwExact(j, k, 1:nbDim) )
          wrc2(j) = dot_product( Dwi(j, k, 1:nbDim)- DwExact(j, k, 1:nbDim), &
               Dwi(j, k, 1:nbDim) - DwExact(j, k, 1:nbDim) )
          wrc3(j) = dot_product( Dwi(j, k, 1:nbDim)- Dwh(j, k, 1:nbDim), &
               Dwi(j, k, 1:nbDim) - Dwh(j, k, 1:nbDim) )
       enddo

       call IntegrateFunction(elem, wrc1(1:Qdof), elem%errTot  )
       call IntegrateFunction(elem, wrc2(1:Qdof), elem%errDisc )
       call IntegrateFunction(elem, wrc3(1:Qdof), elem%errAlg )

       elem%errTot = elem%errTot**0.5
       elem%errDisc = elem%errDisc**0.5
       elem%errAlg = elem%errAlg**0.5


    enddo !k



    deallocate(wrc1, wrc2, wrc3, wExact, DwExact, Dwi, Dwh )

  end subroutine ComputeTotDiscAlgErrorElement

 !> Setting of the exact solution in volume integration nodes
 subroutine SetExactSolutionQnodes(elem, V_rule,  wi, Dwi)
   type(element), intent(inout):: elem   ! elem = element
   type(volume_rule), intent(in) :: V_rule
   real, dimension(1:V_rule%Qdof, 1:ndim), intent(out) :: wi  !exact solution in integ nds
   real, dimension(1:V_rule%Qdof, 1:ndim, 1:nbDim), intent(out) :: Dwi  ! derivatives in ^^^
   real, dimension(:,:), allocatable :: xi, Fx , wS    !
   logical :: grad_analytical
   integer:: Qdof, dof, ndof

   Qdof = V_rule%Qdof

   allocate(xi(1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim) )

   xi(1:Qdof, 1:nbDim) = V_rule%lambda(1:Qdof,1:nbDim)   !

   call ComputeF(elem, Qdof, xi(1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim) )

   call Exact_Sol(Qdof, Fx(1:Qdof, 1:nbDim), wi(1:Qdof, 1:ndim), state%time%ctime  )

   call Exact_Sol_Der(Qdof, Fx(1:Qdof, 1:nbDim), Dwi(1:Qdof, 1:ndim, 1:nbDim), state%time%ctime,&
        grad_analytical )

!   write(*,'(a6,i5,20es16.8)') 'Dwi1',elem%i,state%time%ctime,Dwi(1:3,1,1)
!   write(*,'(a6,i5,20es16.8)') 'Dwi2',elem%i,state%time%ctime,Dwi(1:3,1,2)


   ! if analytical relation for Dw is not given, we copute it as a gradient of a projection
   if(.not. grad_analytical) then  ! MAKE HIGHER ORDER !!!!!
      dof = elem%dof
      ndof = dof * ndim

      ! storring of the numerical solution
      allocate(wS(0:0,1:ndof))
      wS(0, 1:ndof) = elem%w(0,1:ndof)

      ! setting of the exact solution into elem%w at time state%time%ctime
      call SetOneElementIC(elem, .false.)

      call Eval_Dw_Elem(elem, Dwi(1:Qdof, 1:ndim, 1:nbDim) )

!   write(*,'(a6,i5,20es16.8)') 'Dwi1',elem%i,state%time%ctime,Dwi(1:3,1,1)
!   write(*,'(a6,i5,20es16.8)') 'Dwi2',elem%i,state%time%ctime,Dwi(1:3,1,2)

      ! re-freshing of the numerical solution
      elem%w(0, 1:ndof) = wS(0,1:ndof)
      deallocate(wS )

   end if


   deallocate(xi, Fx)

 end subroutine SetExactSolutionQnodes

 !> Setting of the exact solution in arbitrary  integration nodes
 subroutine SetExactSolutionArbitraryNodes(elem, Qdof, xi, wi, Dwi)
   type(element), intent(inout):: elem   ! elem = element
   real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi  !exact solution in integ nds
   real, dimension(1:Qdof, 1:ndim), intent(out) :: wi  !exact solution in integ nds
   real, dimension(1:Qdof, 1:ndim, 1:nbDim), intent(out) :: Dwi  ! derivatives in ^^^
   real, dimension(:,:), allocatable :: Fx, wS     !
   logical :: grad_analytical
   integer:: Qdof, dof, ndof


   allocate(Fx(1:Qdof, 1:nbDim) )


   call ComputeF(elem, Qdof, xi(1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim) )

   call Exact_Sol(Qdof, Fx(1:Qdof, 1:nbDim), wi(1:Qdof, 1:ndim), state%time%ctime  )

   call Exact_Sol_Der(Qdof, Fx(1:Qdof, 1:nbDim), Dwi(1:Qdof, 1:ndim, 1:nbDim), state%time%ctime,&
        grad_analytical )

!   write(*,'(a6,i5,20es16.8)') 'Dwi1',elem%i,state%time%ctime,Dwi(1:3,1,1)
!   write(*,'(a6,i5,20es16.8)') 'Dwi2',elem%i,state%time%ctime,Dwi(1:3,1,2)


   ! if analytical relation for Dw is not given, we copute it as a gradient of a projection
   if(.not. grad_analytical) then  ! MAKE HIGHER ORDER !!!!!
      dof = elem%dof
      ndof = dof * ndim

      ! storring of the numerical solution
      allocate(wS(0:0,1:ndof))
      wS(0, 1:ndof) = elem%w(0,1:ndof)

      ! setting of the exact solution into elem%w at time state%time%ctime
      call SetOneElementIC(elem, .false.)

      call Eval_Dw_Elem(elem, Dwi(1:Qdof, 1:ndim, 1:nbDim) )

!   write(*,'(a6,i5,20es16.8)') 'Dwi1',elem%i,state%time%ctime,Dwi(1:3,1,1)
!   write(*,'(a6,i5,20es16.8)') 'Dwi2',elem%i,state%time%ctime,Dwi(1:3,1,2)

      ! re-freshing of the numerical solution
      elem%w(0, 1:ndof) = wS(0,1:ndof)
      deallocate(wS )

   end if


   deallocate(Fx)

 end subroutine SetExactSolutionArbitraryNodes

 !> compute the \f$ L^2(\Omega) \f$-norm of the difference of \f$ {\bf w}_h\f$
 !> and the exact solution
 subroutine ComputeL2H1ErrorOLD(errL2, errH1, normL2, normH1)
   real, intent (inout) :: errL2, errH1, normL2, normH1
   class(element), pointer :: elem
   real :: err_locL2, err_locH1, norm_locL2, norm_locH1
   integer :: i

   errL2 = 0.
   errH1 = 0.
   normL2 = 0.
   normH1 = 0.

   do i=1,grid%nelem
      elem => grid%elem(i)

      call ComputeL2H1ErrorElementOLD(elem, err_locL2, err_locH1, norm_locL2, norm_locH1)

      errL2 = errL2 + err_locL2
      errH1 = errH1 + err_locH1

      normL2 = normL2 + norm_locL2
      normH1 = normH1 + norm_locH1

      elem%errL2 = err_locL2**0.5
      elem%errH1 = err_locH1**0.5
   enddo

   errL2 = errL2**0.5
   errH1 = errH1**0.5

   normL2 = normL2**0.5
   normH1 = normH1**0.5
 end subroutine ComputeL2H1ErrorOLD




 !> compute the square of the \f$ L^2(K) \f$-norm of the difference of \f$ {\bf w}_h\f$
 !> and the exact solution, \f$ K \f$ = elem
 subroutine ComputeL2H1ErrorElementOLD(elem, errL2, errH1, normL2, normH1)
   type(element), intent(inout):: elem   ! elem = element
   real, intent (inout) :: errL2, errH1, normL2, normH1
   real, dimension(:), allocatable :: w_h   ! w recomputed  in integ nodes
   integer :: i, k, dof

   errL2 = 0.
   errH1 = 0.
   normL2 = 0.
   normH1 = 0.

   dof = elem%dof

   allocate(w_h(1:ndim*dof))

   ! storring of the numerical solution
   w_h(1:ndim*dof) = elem%w(0,1:ndim*dof)

   ! setting of the exact solution into elem%w
   call SetOneElementIC(elem, .false.)
   !!!call SetOneElementAnalyticalIC(elem)  !!! SJEDNOTIT !!!!!

   ! difference of numerical and exact solution
   elem%w(0,1:ndim*dof) =  w_h(1:ndim*dof)  - elem%w(0,1:ndim*dof)

   do k=1,1              ! k = index of component of w
   !do k=1,ndim              ! k = index of component of w

      do i=(k-1)*dof +1, k*dof
         errL2 = errL2 + elem%w(0,i) &
              * dot_product(elem%Mass%Mb(i,1:dof), elem%w(0,(k-1)*dof+1:k*dof) )
         errH1 = errH1 + elem%w(0,i) &
              * dot_product(elem%Stiff%Mb(i,1:dof), elem%w(0,(k-1)*dof+1:k*dof) )

         normL2 = normL2  + w_h(i) * &
              dot_product(elem%Mass%Mb(i,1:dof), w_h((k-1)*dof+1:k*dof) )
         normH1 = normH1  + w_h(i) * &
              dot_product(elem%Stiff%Mb(i,1:dof), w_h((k-1)*dof+1:k*dof) )
      enddo

   enddo

   ! numerical solution back to elem%w(0,:)
   elem%w(0,1:ndim*dof) = w_h(1:ndim*dof)

   !!if(elem%i == 1) write(*,'(a8,6es12.4)') 'Erro:',errL2, errH1, normL2, normH1

   deallocate(w_h)

 end subroutine ComputeL2H1ErrorElementOLD


 !> compute the square of the \f$ H^1(K) \f$-seminorm
 !> from the array elem%wS(0,:) !!!!!!!!!!, \f$ K \f$ = elem
 subroutine ComputeH1NormElement(elem, dof, normL2, normH1)
   type(element), intent(inout):: elem   ! elem = element
   integer, intent(in) :: dof ! dof of the array elem%dofS
   real, intent (inout) :: normL2, normH1
   real, dimension(:,:,:), allocatable :: Der  ! derivative of test functions
   real, dimension(:,:,:), allocatable :: Dwi   ! w recomputed  in integ nodes
   real, dimension(:), allocatable :: func
   integer :: i, k, kst, Qdof
   real :: val

   normL2 = 0.
   normH1 = 0.

   Qdof = elem%Qdof

   allocate(Der(1:dof, 1:nbDim, 1:Qdof) )
   allocate(Dwi(1:Qdof, 1:ndim, 1:nbDim), func(1:Qdof))

   call Eval_Dphi(elem, dof, Der)

   do k=1,ndim              ! k = index of component of Dw
      kst = dof*(k-1) + 1
      do i=1,Qdof
         Dwi(i, k, 1) = dot_product(elem%wS(0,kst:kst+dof-1), Der(1:dof, 1, i) )
         Dwi(i, k, 2) = dot_product(elem%wS(0,kst:kst+dof-1), Der(1:dof, 2, i) )
      enddo

      func(1:Qdof) = Dwi(1:Qdof, k, 1)**2 + Dwi(1:Qdof, k, 2)**2
      call IntegrateFunction(elem, func, val)
      normH1 = normH1 + val

   enddo
   deallocate(Der)

   deallocate(Dwi, func)

 end subroutine ComputeH1NormElement

 !> compute the error \f$ \| {\bf w}(t,x) - {\bf w}_{h\tau}(t,x) \|\f$,
 !> \f$ {\bf w}_{h\tau}\f$ is piecewise polynomial in space as well as time
 !> from BDF reconstruction
 subroutine SpaceTimeErrors( )
   class(element), pointer :: elem
   class(Time_rule), pointer :: T_rule
   real, dimension(:,:,:), allocatable :: Lag_coef
   integer :: Tdeg
   integer :: Gnum, Gdof , j, i
   real :: errL2_loc, errH1_loc, errL2, errH1
	real, dimension(:), allocatable :: errL8, errL8_loc
   real :: t

   associate ( time => state%time )
   select type (time)
      class is ( TimeBDF_t )
         Tdeg = time%deg_actual  ! degree of actual BDF rule

         ! for the integration in time
         Gnum = Tdeg  + 2  ! ???
!         print*, 'T_rule is not allocated for BDF methods. Therefore we use space%G_rule. FR'
!         print*, 'it is now - should be changed'

         T_rule => state%time%T_rule(Gnum)
         Gdof = T_rule%Qdeg

         ! index1 Gaus integ nodes, index2 = 0 -> Lagr functiom, =1 -> its derivative
         allocate(Lag_coef(1:Gdof,0:1, 0:Tdeg) )

         ! integration over the time interval
         do j=1, Gdof
            ! actual time
            t = T_rule%lambda(j)
            !t = 1. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            time%ctime = time%ttime + time%tau(1) * (t -1)
            call SetLagrCoeffs(Tdeg,  t, Lag_coef(j, 0:1, 0:Tdeg) )
         enddo

         errL2 = 0.
         errH1 = 0.
         ! evaluation of errors of pw polynomial approximations
         do i=1,grid%nelem
            elem => grid%elem(i)
            call SpaceTimeErrorsElem(elem, Tdeg, Gnum, Lag_coef(1:Gdof, 0:1, 0:Tdeg), &
                 errL2_loc, errH1_loc )
            errL2 = errL2 + errL2_loc
            errH1 = errH1 + errH1_loc
         enddo

         state%errSTnorm(L2L2) = state%errSTnorm(L2L2) + time%tau(1) * errL2
         state%errSTnorm(L2H1) = state%errSTnorm(L2H1) + time%tau(1)* (errH1 + errL2)
         state%errSTnorm(L2H10)= state%errSTnorm(L2H10) + time%tau(1) * errH1
         state%errSTnorm(L2L2eH1) = state%errSTnorm(L2L2eH1) + time%tau(1) * (errL2 + state%model%Re1* errH1)
         !state%errSTnorm(L8L2) = max( state%errSTnorm(L8L2) , maxval( errL8(1:Gdof) ) )

         state%errSTloc(L2L2) = state%errSTloc(L2L2) + time%tau(1) * errL2
         state%errSTloc(L2H1) = state%errSTloc(L2H1) + time%tau(1) * (errH1 + errL2)
         state%errSTloc(L2H10)= state%errSTloc(L2H10) + time%tau(1) * errH1
         state%errSTloc(L2L2eH1) = state%errSTloc(L2L2eH1) + time%tau(1) * (errL2 + state%model%Re1* errH1)
         state%errSTloc(L8L2) = max(state%errSTloc(L8L2), errL2**0.5)

         deallocate(Lag_coef)

         !write(198,*) time%iter,time%ttime, errL2**0.5, errH1**0.5, errL2, errH1, state%errSTnorm(L2H1)**0.5, &
         !     state%errSTnorm(L2L2eH1)**0.5


      class is ( TimeTDG_t )
          !open(63, file='stdgm-comp' , status='unknown', position='append')
          !write(63,*) 'SpaceTimeErrors in stdgm'
          !write(63,*) ''
          !close(63)

          Gnum = time%Qnum ! We need to integrate function of deg = 2*tdeg
          Gdof = time%T_rule(Gnum)%Qdeg

          errL2 = 0.
          errH1 = 0.
          errL2_loc = 0.
          errH1_loc = 0.

          allocate( errL8_loc(1:Gdof))
          allocate( errL8(1:Gdof))
          errL8( 1:Gdof) = 0.
          errL8_loc( 1:Gdof) = 0.
          do i = 1, grid%nelem
             elem => grid%elem(i)
             call SpaceTimeErrorsElemST( elem, Gdof, errL2_loc, errH1_loc , errL8_loc(1:Gdof) )
             errL2 = errL2 + errL2_loc
             errH1 = errH1 + errH1_loc
             errL8( 1:Gdof) = errL8 + errL8_loc( 1:Gdof)
             elem%errL8 = maxval(errL8( 1:Gdof))
          enddo !i

          state%errSTnorm(L2L2) = state%errSTnorm(L2L2) + time%tau(1) * errL2
          state%errSTnorm(L2H1) = state%errSTnorm(L2H1) + time%tau(1)* (errH1 + errL2)
          state%errSTnorm(L2H10)= state%errSTnorm(L2H10) + time%tau(1) * errH1
          state%errSTnorm(L2L2eH1) = state%errSTnorm(L2L2eH1) + time%tau(1) * (errL2 + state%model%Re1* errH1)
          !state%errSTnorm(L8L2) = max( state%errSTnorm(L8L2) , maxval( errL8(1:Gdof) ) )

          state%errSTloc(L2L2) = state%errSTloc(L2L2) + time%tau(1) * errL2
          state%errSTloc(L2H1) = state%errSTloc(L2H1) + time%tau(1) * (errH1 + errL2)
          state%errSTloc(L2H10)= state%errSTloc(L2H10) + time%tau(1) * errH1
          state%errSTloc(L2L2eH1) = state%errSTloc(L2L2eH1) + time%tau(1) * (errL2 + state%model%Re1* errH1)
          state%errSTloc(L8L2) = max(state%errSTloc(L8L2),  maxval( errL8(1:Gdof) )**0.5 )



          !	print*, 'state%errSTnorm(L2L2)' , state%errSTnorm(L2L2)
          !	print*, 'state%errSTnorm(L2H1)' , state%errSTnorm(L2H1)
          !	print*, 'state%errSTnorm(L2H10)', state%errSTnorm(L2H10)
          !	print*, 'state%errSTnorm(L8L2)' , state%errSTnorm(L8L2)

          deallocate( errL8_loc )
          deallocate( errL8 )

      class default
         stop 'unknown time_disc method in SpaceTimeErrors'
   end select
   end associate

   !write(99,'(2(a5,6es12.4))') 'GLOB',state%errSTnorm(1:6),'LOC ',state%errSTloc(1:6)




end subroutine SpaceTimeErrors


 !> compute the error \f$ \| {\bf w}(t,x) - {\bf w}_{h\tau}(t,x) \|\f$,
 !> \f$ {\bf w}_{h\tau}\f$ is piecewise polynomial in space as well as time
 !> from BDF reconstruction on one element
 subroutine SpaceTimeErrorsElem(elem, Tdeg, Gnum, Lag_coef, errL2, errH1 )
   type(element), intent(inout) :: elem
   integer, intent(in) :: Tdeg    ! degree of Lagrangian reconstruction
   integer, intent(in) :: Gnum    ! degree of Gauss quandrature
   real, dimension(1:state%space%G_rule(Gnum)%Qdof, 0:1, 0:Tdeg), intent(in) :: Lag_coef
   real, intent(out) :: errL2, errH1
   type(Gauss_rule), pointer :: G_rule
   real, dimension(:,:), allocatable :: wi, wE      ! approx solution in integ nds
   real, dimension(:,:,:), allocatable :: Dwi, DwE   ! Der of approxt solution in integ nds
   real, dimension(:), allocatable :: wrc1, wrc2 ! temporary arrays
   integer:: Qdof, Gdof, i, j, k, ndof
   real :: t, val, errL2t, errH1t

   G_rule => state%space%G_rule(Gnum)
   Gdof = G_rule%Qdof

   Qdof = elem%Qdof
   ndof = elem%dof*ndim

   allocate(wi(1:Qdof, 1:ndim),  wE(1:Qdof, 1:ndim))
   allocate(Dwi(1:Qdof, 1:ndim, 1:nbDim),  DwE(1:Qdof, 1:ndim, 1:nbDim))
   allocate(wrc1(1:Qdof),  wrc2(1:Qdof))

   errL2 = 0.
   errH1 = 0.
   ! integration over the time interval
   do j=1, Gdof
      ! actual time
      t = G_rule%lambda(j)

      !t = 1. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      state%time%ctime = state%time%ttime + state%time%tau(1) * (t - 1)

      ! setting of the exact solution in integ nodes
      call SetExactSolutionQnodes(elem, state%space%V_rule(elem%Qnum), &
           wE(1:Qdof, 1:ndim), DwE(1:Qdof, 1:ndim, 1:nbDim))

      ! setting of the polynomial reconstruction of w_h and D w_h in time integ nodes
      call Eval_wht_Elem(Tdeg, Lag_coef(j,0,0:Tdeg), elem, &
           wi(1:Qdof,1:ndim), Dwi(1:Qdof,1:ndim, 1:nbDim ) )

      !if(elem%i == 1 .and. j==1 ) then
         !write(197,'(a6,2i5,20es16.8)') 'wi(t)',elem%i,j,state%time%ctime,wi(1:3,1)
         !write(197,'(a6,2i5,20es16.8)') 'wE(t)',elem%i,j,state%time%ctime,wE(1:3,1)
         !write(198,'(a6,2i5,20es16.8)') 'Dwi(t)',elem%i,j,state%time%ctime,Dwi(1:3,1, 1),Dwi(1:3,1,2)
         !write(198,'(a6,2i5,20es16.8)') 'DwE(t)',elem%i,j,state%time%ctime,DwE(1:3,1, 1),DwE(1:3,1,2)
      !   write(*,*)'-------------------------------------------------',state%time%ctime
      !   write(*,'(a4,i5,80es12.4)') 'er*i',elem%i, wi(:,1)
      !   write(*,'(a4,i5,80es12.4)') 'er*E',elem%i, wE(:,1)
      !   write(*,'(a4,i5,80es12.4)') 'er*e',elem%i,  wi(1:Qdof,1:ndim) - wE(1:Qdof,1:ndim)
      !endif


      ! error = approximate - exact
      wi(1:Qdof,1:ndim) = wi(1:Qdof,1:ndim) - wE(1:Qdof,1:ndim)
      Dwi(1:Qdof,1:ndim, 1:nbDim) = Dwi(1:Qdof,1:ndim, 1:nbDim) - DwE(1:Qdof,1:ndim, 1:nbDim)

      errL2t = 0.
      errH1t = 0.
      ! evaluation of the error
      !do k=1,1              ! k = index of component of w
      do k=1,ndim              ! k = index of component of w

         ! L2-norm
         wrc1(1:Qdof) =  wi(1:Qdof, k)**2
         !wrc2(1:Qdof) = (wEt(1:Qdof, k) )**2

         call IntegrateFunction(elem, wrc1(1:Qdof), val  )
         errL2t = errL2t + val

         !call IntegrateFunction(elem, wrc2(1:Qdof), val  )
         !normL2 = normL2 + val

         ! H1-seminorm
         do i=1,Qdof
            wrc1(i) = dot_product( Dwi(i, k, 1:nbDim), Dwi(i, k, 1:nbDim) )
            !wrc2(j) = dot_product( DwE(j, k, 1:nbDim), DwE(j, k, 1:nbDim))
         enddo

         call IntegrateFunction(elem, wrc1(1:Qdof), val  )
         errH1t = errH1t + val

         !call IntegrateFunction(elem, wrc2(1:Qdof), val  )
         !normH1 = normH1 + val
      end do  !k=1,ndim

      !if(elem%i == 1) write(*,'(a8,6es12.4)') '....:',errL2t, errH1t

      errL2 = errL2 + G_rule%weights(j) * errL2t
      errH1 = errH1 + G_rule%weights(j) * errH1t


   enddo

      !if(elem%i == 1) write(*,'(a8,6es12.4)') 'Er_ht:',errL2, errH1


   deallocate(wi, Dwi, wE, DwE)


 end subroutine SpaceTimeErrorsElem

 !> compute the error \f$ \| {\bf w}(t,x) - {\bf w}_{h\tau}(t,x) \|\f$,
 !> \f$ {\bf w}_{h\tau}\f$ is piecewise polynomial in space as well as time
 ! wouldn't work if Gnum /= T_rule%Qdof
 subroutine SpaceTimeErrorsElemST(elem, Gnum, errL2, errH1, errL8)
   type(element), intent(inout) :: elem
   integer, intent(in) :: Gnum    ! degree of time-Gauss quandrature
   real, intent(out) :: errL2, errH1
   real, dimension(:), intent(out) :: errL8
   class(Time_rule), pointer :: T_rule
   real, dimension(:,:), allocatable :: wi, wE      ! approx solution in integ nodes
   real, dimension(:,:,:), allocatable :: Dwi, DwE   ! Der of approxt solution in integ nodes
   real, dimension(:), allocatable :: wrc1, wrc2 ! temporary arrays
   integer:: Qdof, Gdof, i, j, k, ndof
   real :: t, val, errL2t, errH1t,cTime


   cTime = state%time%ctime

   T_rule => state%time%T_rule(Gnum)
      Gdof = T_rule%Qdeg

      if (Gnum /=T_rule%Qdeg) then
         print*, 'Problem in SpaceTimeErrorsElemST, Gnum <> T_rule%Qdof'
         stop
      endif

      Qdof = elem%Qdof
      ndof = elem%dof*ndim


      allocate(wi(1:Qdof, 1:ndim),  wE(1:Qdof, 1:ndim))
      allocate(Dwi(1:Qdof, 1:ndim, 1:nbDim),  DwE(1:Qdof, 1:ndim, 1:nbDim))
      allocate(wrc1(1:Qdof),  wrc2(1:Qdof))

      errL2 = 0.
      errH1 = 0.
      errL8(1:Gdof) = 0.
      ! integration over the time interval
      do j=1, Gdof
      !do j= Gdof, 1, -1
         ! actual time
         t = T_rule%lambda(j)

         !t = 1. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         state%time%ctime = state%time%ttime + state%time%tau(1) * (t - 1)

         !if(elem%i == 1) write(*,'(a8, i2, 2es14.6)') 'STerrs',j,state%time%ctime, state%time%ttime

         ! setting of the exact solution in integ nodes
         call SetExactSolutionQnodes(elem, state%space%V_rule(elem%Qnum), &
              wE(1:Qdof, 1:ndim), DwE(1:Qdof, 1:ndim, 1:nbDim))

         ! setting of the space-time polynomial of whST and D whST in time integ nodes
         call Eval_whST_Elem(elem, Gdof, j, wi(1:Qdof,1:ndim), Dwi(1:Qdof,1:ndim,1:nbDim)  )

   !      if(elem%i == 1  ) then
   !         write(*,'(a6,2i5,20es16.8)') 'wi(t)',elem%i,j,state%time%ctime,wi(1:3,1)
   !         write(*,'(a6,2i5,20es16.8)') 'wE(t)',elem%i,j,state%time%ctime,wE(1:3,1)
   !         write(*,'(a6,2i5,20es16.8)') 'Dwi(t)',elem%i,j,state%time%ctime,Dwi(1:3,1, 1),Dwi(1:3,1,2)
   !         write(*,'(a6,2i5,20es16.8)') 'DwE(t)',elem%i,j,state%time%ctime,DwE(1:3,1, 1),DwE(1:3,1,2)
   !         write(*,*)'-------------------------------------------------',state%time%ctime
   !      !   write(*,'(a4,i5,80es12.4)') 'er*i',elem%i, wi(:,1)
   !      !   write(*,'(a4,i5,80es12.4)') 'er*E',elem%i, wE(:,1)
   !      !   write(*,'(a4,i5,80es12.4)') 'er*e',elem%i,  wi(1:Qdof,1:ndim) - wE(1:Qdof,1:ndim)
   !      endif


         ! error = approximate - exact
         wi(1:Qdof,1:ndim) = wi(1:Qdof,1:ndim) - wE(1:Qdof,1:ndim)
         Dwi(1:Qdof,1:ndim, 1:nbDim) = Dwi(1:Qdof,1:ndim, 1:nbDim) - DwE(1:Qdof,1:ndim, 1:nbDim)

         errL2t = 0.
         errH1t = 0.
         ! evaluation of the error
         !do k=1,1              ! k = index of component of w
         do k=1,ndim              ! k = index of component of w


            ! L2-norm
            wrc1(1:Qdof) =  wi(1:Qdof, k)**2
            !wrc2(1:Qdof) = (wEt(1:Qdof, k) )**2


            call IntegrateFunction(elem, wrc1(1:Qdof), val  )
            errL2t = errL2t + val

            !if(elem%i == 1) &
            !  write(*,'(a4,8es14.6)')'-ST--',state%time%ctime, val, wrc1(1:6)

            !call IntegrateFunction(elem, wrc2(1:Qdof), val  )
            !normL2 = normL2 + val

            ! H1-seminorm
            do i=1,Qdof
               wrc1(i) = dot_product( Dwi(i, k, 1:nbDim), Dwi(i, k, 1:nbDim) )
               !wrc2(j) = dot_product( DwE(j, k, 1:nbDim), DwE(j, k, 1:nbDim))
            enddo

            call IntegrateFunction(elem, wrc1(1:Qdof), val  )
            errH1t = errH1t + val

            !call IntegrateFunction(elem, wrc2(1:Qdof), val  )
            !normH1 = normH1 + val
         end do  !k=1,ndim

         !if(elem%i == 1) write(*,'(a8,6es12.4)') '....:',errL2t, errH1t

         errL2 = errL2 + T_rule%weights(j) * errL2t
         errH1 = errH1 + T_rule%weights(j) * errH1t
         errL8(j) = errL2t



      enddo  !j=1,Gdof

         !if(elem%i == 1) write(*,'(a8,6es12.4)') 'Er_ht:',errL2, errH1


      deallocate(wi, Dwi, wE, DwE)

      state%time%ctime = cTime


 end subroutine SpaceTimeErrorsElemST



 !> computation of
 subroutine H1seminormOfErrorInTimeDerElem(elem, t, rvalue)
   type(element), intent(in) :: elem
   real, dimension(:), pointer :: weights
   real, dimension(:,:,:), allocatable :: Der
   real, dimension(:,:,:), pointer:: Dphi ! pointers to test functions
   real, dimension(:,:,:), allocatable :: wi
   real, dimension(:,:,:), allocatable :: utgrad
   real, dimension(:,:), allocatable :: intercalc
   real, dimension(:,:), allocatable :: x, Fx
   real, dimension(1:ndim), intent(out) :: rvalue
   real :: t
   integer :: k, Qdof, Qnum, ist, l, j, i


    rvalue(1:ndim) = 0
    Qdof = elem%Qdof
    Qnum = elem%Qnum
    weights => state%space%V_rule(Qnum)%weights(1:Qdof)

    allocate(x(1:Qdof, 1:2))
    x(1:Qdof, 1:2) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:2)
    allocate(Fx(1:Qdof, 1:2))
    call ComputeF(elem, Qdof, x, Fx)

    allocate(utgrad(1:ndim, 1:2, 1:Qdof) )
    do l=1, Qdof
       call Der_TimeDer_Exact_Scalar(Fx(l, 1:2), utgrad(1:ndim, 1:2, l), t)
    enddo !l

    allocate (Der(1:elem%dof,1:2,1:Qdof))

    ! derivatives of the reference test functions
    Dphi => state%space%V_rule(Qnum)%Dphi(1:elem%dof, 1:2, 1:Qdof)
    do j=1, elem%dof
       Der(j, 1, 1:Qdof) = &
            elem%F%D1F0(1,1) * Dphi(j, 1, 1:Qdof) &
            +elem%F%D1F0(1,2)* Dphi(j, 2, 1:Qdof)

       Der(j, 2, 1:Qdof) = &
            elem%F%D1F0(2,1) * Dphi(j, 1, 1:Qdof) &
            +elem%F%D1F0(2,2)* Dphi(j, 2, 1:Qdof)
    enddo !j

    allocate( wi(0:1, 1:ndim, 1:state%space%max_dof) )
    allocate( intercalc(1:ndim, 1:2) )
    !intercalc(1:ndim, 1:2) = 0. ! CHECK

    do k=1, ndim
       ist = (k-1)*elem%dof + 1
          wi(0, k, 1:elem%dof) = elem%w(0,ist:ist+elem%dof-1)
          wi(1, k, 1:elem%dof) = elem%w(1,ist:ist+elem%dof-1)
          do l=1, Qdof
             do i=1, 2
             intercalc(k, i) = utgrad(k, i, l) &
                  - dot_product(wi(0, k, 1:elem%dof), Der(1:elem%dof, i, l))&
                  * (1/(state%time%tau(1)) )&
                  + dot_product(wi(1, k, 1:elem%dof), Der(1:elem%dof, i, l))&
                  * (1/(state%time%tau(1)) )
          enddo !i
          rvalue(k) = rvalue(k) + weights(l) &
               * dot_product(intercalc(k, 1:2), intercalc(k, 1:2)) * elem%F%JF0 * 0.5
       enddo !l

    enddo !k

    deallocate(Der, wi, utgrad, intercalc, x, Fx)

  end subroutine H1seminormOfErrorInTimeDerElem

  !> compute the \f$ H^1 \f$-seminorm of the time derivative
  subroutine ComputeH1ErrorTimeDer( )
    class(element), pointer :: elem
    real, dimension(1:ndim, 0:1) :: errgradtime
    integer :: i

    do i=1,grid%nelem
       elem => grid%elem(i)
       call H1seminormOfErrorInTimeDerElem(elem, state%time%ttime, errgradtime(1:ndim, 0))
       call H1seminormOfErrorInTimeDerElem(elem, state%time%ttime - state%time%tau(1), &
            errgradtime(1:ndim, 1))

       state%errSTnorm(L2L2eH1) = state%errSTnorm(L2L2eH1) + &
            state%time%tau(1)*(errgradtime(1, 0) + errgradtime(1, 1)) / 2  ! ONLY for ndim = 1
    enddo !i

  end subroutine ComputeH1ErrorTimeDer


  !> compute errE\f$= \int_{K} (w_1 - w_2)^2\ dx \f$
  !> and normE \f$=\int_{K} (w_2)^2\ dx \f$
  !>
  !> elem = \f$ K\in {\cal T}_h\f$,
  !> w1, w2 = \f$ w_1, w_2 \in {\bf S}_{hp}\f$
  subroutine ComputeElementConvErrors(elem, w1, w2, normE, errE, norm8, err8)
    type(element):: elem
    real, dimension(1:elem%dof*ndim), intent (in) :: w1, w2
    real, intent(inout) :: normE, errE    ! norm and error in L^2 norm
    real, intent(inout) :: norm8, err8    ! norm and error in L^{\infty} norm

    real, dimension(:), allocatable :: wrc1, wrc2   ! w recomputed  in integration nodes
    real, dimension(:,:),   pointer:: phi ! local store arrays
    integer ::  j, k, kst, dof, Qnum, Qdof
    real :: val

    dof = elem%dof
    Qnum = elem%Qnum
    Qdof = elem%Qdof

    allocate(wrc1(1: Qdof), wrc2(1:Qdof) )

    phi => state%space%V_rule(Qnum)%phi(1:dof, 1:Qdof)

    normE = 0.
    errE = 0.

    !if(elem%F%iFlin ) then
    !   ! only for ORTHONORMAL basis
    !   do k=1,ndim              ! k = index of component of w
    !      kst = dof*(k-1) + 1
    !
    !      errE = errE +  elem%area*(5 -elem%type) &
    !           * dot_product(w1(kst: k*dof)-w2(kst: k*dof), w1(kst: k*dof)-w2(kst: k*dof) )
    !
    !      normE = normE +  elem%area*(5 -elem%type) &
    !           * dot_product(w2(kst: k*dof), w2(kst: k*dof) )
    !   enddo
    !else
       ! setting of w in integration nodes

       do k=1,ndim              ! k = index of component of w
          kst = dof*(k-1) + 1
          do j=1,Qdof
             wrc1(j) = dot_product(w1((k-1)*dof+1: k*dof), phi(1:dof,j) )
             wrc2(j) = dot_product(w2((k-1)*dof+1: k*dof), phi(1:dof,j) )

             !if(j == 1) then
             !if(elem%i == 1) then
                !write(*,'(a4,12es11.3)') 'w2',w2((k-1)*dof+1: k*dof)
                !write(*,'(a4,12es11.3)') 'phi',phi(1:dof, j)
             !   write(*,'(a4,12es11.3)') 'wrc1',wrc1(j)
             !   write(*,'(a4,12es11.3)') 'wrc2',wrc2(j)
             !endif

             err8 = max(err8, abs(wrc1(j) - wrc2(j)) )
             norm8 = max(norm8, abs(wrc2(j)) )

             wrc1(j) = (wrc1(j) - wrc2(j))**2
             wrc2(j) = wrc2(j)**2
          enddo

          call IntegrateFunction(elem, wrc1(1:Qdof), val  )
          errE = errE + val

          call IntegrateFunction(elem, wrc2(1:Qdof), val  )
          !if(val > 1E-2) write(*,'(a3,i5, 8es11.3)' ) 'wrc',elem%i, val, wrc1(1:Qdof)!, normE

          normE = normE + val
       end do
    !endif

    deallocate(wrc1, wrc2 )

  end subroutine ComputeElementConvErrors

  subroutine ComputeConvErrors(norm, errL2, norm8, err8)
    real, intent(out) :: norm, errL2
    real :: normE, errE
    real :: norm8, err8
    integer :: i

    errL2 = 0.
    norm = 0.
    err8 = 0.
    norm8 = 0.

    do i = 1, grid%nelem
       call ComputeElementConvErrors(grid%elem(i), &
            grid%elem(i)%w(0,:), grid%elem(i)%w(1,:),normE, errE, norm8, err8)
       errL2 = errL2 + errE
       norm = norm + normE
       !write(500+state%time%iter, *) i, grid%elem(i)%xc(:), errE, normE, err8, norm8
    enddo
  end subroutine ComputeConvErrors

  subroutine ComputeConvErrorsTDG(norm, errL2, norm8, err8)
    real, intent(out) :: norm, errL2
    real :: normE, errE
    real :: norm8, err8
    integer :: i

    errL2 = 0.
    norm = 0.
    err8 = 0.
    norm8 = 0.

    if(state%time%iter <=2) print*,'@@@', ' CHECK ComputeConvErrorsTDG in euler.f90'
    do i = 1, grid%nelem
       call ComputeElementConvErrors(grid%elem(i), &
            grid%elem(i)%w(0,:), &
            !sum(grid%elem(i)%wold, 1) / size(grid%elem(i)%wold, 1), &
            sum(grid%elem(i)%w, 1) / size(grid%elem(i)%w, 1), &
            normE, errE, norm8, err8)
       errL2 = errL2 + errE
       norm = norm + normE
    enddo
  end subroutine ComputeConvErrorsTDG


  !> evaluation of the difference of cDLM_max - cDLM_min within last 10 % of time steps
  !> cDLM(:,1) = cD, cDLM(:,2) = cL, cDLM(:,3) = cM
  !> EcDLM = cDLM_max - cDLM_min
  subroutine ComputeEoFC()
    real :: cDLM_max, cDLM_min
    integer :: start_it, fin_it
    integer :: i,k

    !start_it = aint(0.9*state%time%iter)
    !start_it = aint(0.95*state%time%iter)
    start_it = state%time%iter - state%time%iter_loc + aint(0.9*state%time%iter_loc)

    start_it = max(start_it,1)
    fin_it = state%time%iter

    do k=1,3
       cDLM_max = state%cDLM(start_it,k)
       cDLM_min = state%cDLM(start_it,k)

       if (start_it > 0) then
          do i=start_it,fin_it
             if (state%cDLM(i,k) > cDLM_max) then
                cDLM_max=state%cDLM(i,k)
             else if (state%cDLM(i,k) < cDLM_min) then
                cDLM_min=state%cDLM(i,k)
             endif
          enddo
          state%EcDLM(k)  = cDLM_max-cDLM_min
       else !1st iteration
          state%EcDLM(k) = 1.
       endif
    enddo
  end subroutine ComputeEoFC

   !> compute the square of the \f$ L^2(\Omega) \f$-norm of
 !>  the difference of \f$ {\bf w}_h\f$
 !> and the exact solution in timeNode node of the Time quadrature formula
  subroutine ComputeL2H1ErrorST(timeNode, errL2, errH1, normL2, normH1)
    real, intent (in) :: timeNode
    real, intent (inout) :: errL2, errH1, normL2, normH1
    class(element), pointer :: elem
    class(Time_rule), pointer :: T_rule
    real :: err_locL2, norm_locL2, err_locH1, norm_locH1, cTime
    integer :: i, Tnum, dof

    Tnum = state%time%Qnum ! We need to integrate function of deg = 2*tdeg - works with Gauss rule

    T_rule => state%time%T_rule(Tnum)

       dof = T_rule%Qdeg


       cTime = state%time%ctime
       state%time%ctime = state%time%ttime + state%time%tau(1) * (timeNode -1)

       errL2 = 0.
       normL2 = 0.
       errH1 = 0.
       normH1 = 0.
       !eval the basis functions phi in time in timeNode

       ! new Trule basis is evaluated in this subroutine
       !dof, Qdof, xi - nodes in (0,1), phi, dphi

!       call Eval_LegendrePolynomialsOneNode( state%time%max_Tdof + 1, timeNode, &
!         T_rule%Phi(1:(state%time%max_Tdof + 1) , -2), T_rule%DPhi(1: state%time%max_Tdof + 1 , -2) )

       call Eval_LegendrePolynomialsOneNode( state%time%Qnum , timeNode, &
         T_rule%Phi(1:(state%time%Qnum) , -2), T_rule%DPhi(1:state%time%Qnum, -2) )

       do i=1,grid%nelem
          elem => grid%elem(i)
          call  ComputeL2H1ErrorElementST(elem, err_locL2, err_locH1, &
               norm_locL2, norm_locH1)

          errL2 = errL2 + err_locL2
          !      print*, 'errL2 i', errL2
          normL2 = normL2 + norm_locL2
          errH1 = errH1 + err_locH1
          normH1 = normH1 + norm_locH1

          !FERROR elem error should be probably computed here, not in bdf version
          !  elem%errL2 = err_locL2**0.5
          !  elem%errH1 = err_locH1**0.5

          !  if(state%time%iter == 0) elem%eta(IC,1) = elem%errL2

       enddo


       errL2 = errL2**0.5
       errH1 = errH1**0.5
       normL2 = normL2**0.5
       normH1 = normH1**0.5

       ! print*, 'errL2', errL2
       state%time%ctime = cTime


  end subroutine ComputeL2H1ErrorST

 !> compute the square of the \f$ L^2(K) \f$-norm of the difference of \f$ {\bf w}_h\f$
 !> and the exact solution, \f$ K \f$ = elem in time = state%time%ctime
 subroutine ComputeL2H1ErrorElementST(elem, errL2, errH1, normL2, normH1)
   type(element), intent(inout):: elem   ! elem = element
   real, intent (inout) :: errL2, normL2, errH1, normH1
   real, dimension(:), allocatable :: wrc1, wrc2   ! w recomputed  in integ nodes
   real, dimension(:,:), allocatable :: wExact     ! exact solution  in integ nodes
   real, dimension(:,:,:), allocatable :: DwExact  ! Der of exact solution in integ nds
   real, dimension(:,:), allocatable :: wi
   real, dimension(:,:,:), allocatable :: Dwi
   integer :: j, k, dof, Qnum, Qdof
   real :: val

   errL2 = 0.
   errH1 = 0.
   normL2 = 0.
   normH1 = 0.
   !elem%errL8 = 0.

   dof = elem%dof
   Qnum = elem%Qnum
   Qdof = elem%Qdof

   allocate(wrc1(1:Qdof), wrc2(1:Qdof) )
   allocate(wExact(1:Qdof, 1:ndim) )
   allocate(DwExact(1:Qdof, 1:ndim, 1:nbDim))
   allocate(wi(1:Qdof, 1:ndim))
   allocate(Dwi(1:Qdof, 1:ndim, 1:nbDim))

   ! setting of the exact solution in integ nodes
   call SetExactSolutionQnodes(elem, state%space%V_rule(elem%Qnum), &
        wExact(1:Qdof, 1:ndim), DwExact(1:Qdof, 1:ndim, 1:nbDim))


   ! setting of the approximate solution in integ nodes
   !call Eval_w_Elem(elem, wi(1:Qdof, 1:ndim) )
   !call Eval_Dw_Elem(elem, Dwi(1:Qdof, 1:ndim, 1:nbDim) )

   ! setting of the space-time polynomial of whST and D whST in time defined before, saved in %Phi(:,-2)
   call Eval_whST_Elem(elem, elem%TQnum, -2, wi(1:Qdof,1:ndim), Dwi(1:Qdof,1:ndim,1:nbDim)  )

   !elem%errL8 = 0.

   ! setting of w in integration nodes
   do k=1,ndim              ! k = index of component of w

      ! L2-norm
      ! setting of the exact solution in integ nodes: array wrc2
      wrc1(1:Qdof) = (wExact(1:Qdof, k) - wi(1:Qdof, k) )**2
      wrc2(1:Qdof) = (wExact(1:Qdof, k) )**2

      ! L^{\infty}-norm
      !elem%errL8 = max( elem%errL8, maxval( wrc1(1:Qdof) )**0.5)

      call IntegrateFunction(elem, wrc1(1:Qdof), val  )
      errL2 = errL2 + val

     ! if(elem%i == 8) then
     !    write(*,*)'-------------------------------------------------',state%time%ctime
    !     write(*,'(a4,i5,80es12.4)') 'er!i',elem%i, wi(:,1)
    !     write(*,'(a4,i5,80es12.4)') 'er!E',elem%i, wExact(:,1)
     !    write(*,'(a4,i5,80es12.4)') 'er!e',elem%i, wrc1(:)
     ! endif

      call IntegrateFunction(elem, wrc2(1:Qdof), val  )
      normL2 = normL2 + val

      ! H1-seminorm
      ! setting of the exact solution in integ nodes

      do j=1,Qdof
         wrc1(j) = dot_product( Dwi(j, k, 1:nbDim)- DwExact(j, k, 1:nbDim), &
              Dwi(j, k, 1:nbDim) - DwExact(j, k, 1:nbDim) )
         wrc2(j) = dot_product( DwExact(j, k, 1:nbDim), DwExact(j, k, 1:nbDim))
      enddo

      call IntegrateFunction(elem, wrc1(1:Qdof), val  )
      errH1 = errH1 + val
      !print*,'###',errH1, val

      call IntegrateFunction(elem, wrc2(1:Qdof), val  )
      normH1 = normH1 + val
   end do
   !if(elem%i == 8)  write(*,'(a8,6es12.4)') 'Erro:',errL2, errH1  !, normL2, normH1

   deallocate(wrc1, wrc2, wExact, DwExact, wi, Dwi )

 end subroutine ComputeL2H1ErrorElementST

end module error_subs
