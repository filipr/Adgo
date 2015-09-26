!> estimates for AMA adaptation
module higher_order_estims

  use main_data  ! contains type(mesh) :: grid for computation
  use problem_oper
  use euler_problem
  use estimates
  use plot_geom
  use eval_sol
  use ama_hp_interpol
  use AMA_estims
  use local_problem
  use local_problem2

  implicit none

  public:: ComputeHigherOrderEstims
  public:: ComputeHigherOrderEstims_OPT
  public:: SetHigherOrder_HP_adapt
  public:: ComputeHigherOrderElem
  public:: ComputeHigherOrderElem_HP
  public:: ComputeHigherOrderElem_OPT
  public:: Decide_hp_type_adapt

contains

  !> compute the error estimates by a higher order reconstruction
  subroutine ComputeHigherOrderEstims( )
    class(element), pointer :: elem
    real, dimension(:,:), allocatable :: wp ! array  with solution in vertexes
    integer :: i, j, k,  imt, imt1, is, dof
    logical :: loc_implicitly
    character(len=15) :: file1, file2
    character(len=5) :: ch5

    if(nbDim /=2 ) then
       print*,' ## ComputeHigherOrderEstims only for nbDim = 2!'
       stop
    endif

    ! create the list of elements sharing at least a vertex with elem
    call SeekElemSupports(grid)  

    ! setting of the array elem%wS
    do i=1,grid%nelem
       elem => grid%elem(i)
       dof = elem%dof
       allocate(elem%wS( 1:ndim, 1:elem%dof ) )
       do k=1,ndim
          elem%wS(k, 1:dof) = elem%w(0, (k-1)*dof+1 : k*dof) ! first component (or density)
       enddo
    enddo


    state%estim( : , :) = 0.
    do i=1,grid%nelem
       elem => grid%elem(i)
       call ComputeHigherOrderElem(elem)

       !call ComputeHigherOrderElem_HP(elem)

       state%estim(1:max_eta, : ) = state%estim(1:max_eta, : ) + elem%eta(1: max_eta, : ) 

    enddo

    write(*,'(a25, 4es12.4)') 'HO_rec error estims', & 
         state%estim(HO_estim_L2_p1: HO_estim_H1_p1, :)**0.5, &
         state%estim(HO_estim_L2_p2: HO_estim_H1_p2, :)**0.5

    do i=1,grid%nelem
       elem => grid%elem(i)
       deallocate(elem%wS)
       deallocate(elem%supp)
    enddo

  end subroutine ComputeHigherOrderEstims

  !> compute the error estimates by a higher order reconstruction
  subroutine ComputeHigherOrderEstims_OPT( )
    class(element), pointer :: elem
    real, dimension(:,:), allocatable :: wp ! array  with solution in vertexes
    integer :: i, j, k,  imt, imt1, is, dof, dofP
    integer :: ideg
    character(len=15) :: file1, file2
    character(len=5) :: ch5

    if(nbDim /=2 ) then
       print*,' ## ComputeHigherOrderEstims only for nbDim = 2!'
       stop
    endif

    ! create the list of elements sharing at least a vertex with elem
    call SeekElemSupports(grid)  

    ! setting of the array elem%wS
    do i=1,grid%nelem
       elem => grid%elem(i)
       dofP = DOFtriang(elem%deg + 1)
       allocate(elem%wS( 1:2*ndim, 1: dofP ) )
       
       elem%wS(:,:) = 0.

       dof = elem%dof
       do k=0,ndim-1
          elem%wS(k+1, 1:dof) = elem%w(0, k*dof+1 : k*dof+dof)
       enddo

       ! for storing of the results to ho-local.f90
       !allocate(elem%wSS(1:ndim, 1:dofP, 0:2 ) )

       !elem%eta( : , 1 :ndim  ) = 0.
    enddo

    !print*,'REMOVe A'
    state%estim( : , :) = 0.

    !do ideg = -1, 0   ! reconstruction of function which is polynom of degree  elem%deg + ideg 
    !do ideg = -1, 1
   ! do ideg = -2, 0

      ! do i=1,grid%nelem
      !    elem => grid%elem(i)
      !    dof = elem%dof
      !    dofP = DOFtriang(elem%deg + ideg)

          ! if(ideg <= 0) then
          !    elem%wS(:,:) = 0.
          !    if(dofP > 0) then
          !       do k=0,ndim-1
          !          elem%wS(k+1, 1:dofP) = elem%w(0, k*dof+1 : k*dof+dofP) ! first component (or density)
          !       enddo
          !    endif

          ! else 
          !    do k=1,ndim
          !       elem%wS(k, 1:dofP) = elem%wS(ndim + k, 1:dofP)
          !    enddo

          ! endif

      ! enddo


       do i=1,grid%nelem
          elem => grid%elem(i)

          !if(elem%deg + ideg >= 0) then
          call ComputeHigherOrderElem_OPT(elem, 0)  ! ideg
          !endif

          !call ComputeHigherOrderElem_HP(elem)
          
          !if(elem%i == 1)
          !print*,'REMOVe B',i, grid%nelem
          state%estim(1:max_eta, : ) = state%estim(1:max_eta, : ) + elem%eta(1: max_eta, : ) 
          !write(*,'(a10,2i5,6es12.5)') 'recovery:',elem%i,HO_recovery,elem%eta(HO_recovery, k),  state%estim(HO_recovery,:)
       enddo

     !  if(ideg == 0) &
        write(*,'(a25, 4es12.4)') 'HO_rec error estimS', & 
             state%estim(HO_estim_H1_p2, :)**0.5, &
             state%estim(HO_recovery, :)**0.5, &
             state%estim(HO_rec_estim, :)**0.5

       ! write(*,'(a25, 6es12.4)') 'HO_rec TRNC arrays', & 
       !      state%estim(HO_trunc_L2_p0: HO_trunc_H1_p0, :)**0.5, &
       !      state%estim(HO_trunc_L2_p1: HO_trunc_H1_p1, :)**0.5, &
       !      state%estim(HO_trunc_L2_p2: HO_trunc_H1_p2, :)**0.5
    !enddo

    do i=1,grid%nelem
       elem => grid%elem(i)
       deallocate(elem%wS)
       deallocate(elem%supp)
    enddo


    !!

    ! !variant of high order Riemann metric
    ! file1 = 'metrixA00000'
    ! file2 = 'metrixS00000'
    ! !file2 = 'rgabcxA00000'

    ! is = 0
    ! if(state%space%adapt%adapt_level > 0) is = int(log(1. * state%space%adapt%adapt_level)/log(10.)) 

    ! write( ch5, '(i5)' ) state%space%adapt%adapt_level  ! change the format if num_size /= 5 !!!
    ! file1(12-is: 12)  = ch5(5-is:5)
    ! file2(12-is: 12)  = ch5(5-is:5)

    
    ! imt = 24
    ! open(imt, file=file1, status='UNKNOWN')
    ! do i=1,grid%nelem
    !    elem => grid%elem(i)
    !    !    if( mod(i, 3) == 1) &
    !    !if( dot_product(elem%xc(:), elem%xc(:))**0.5 < 4E-2) &
    !    !    if(abs(elem%xc(1) - 1.5) < 0.25 .and. elem%xc(2) > 1.75 ) &
    !    call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )

    ! enddo
    ! close(imt)

    if(state%space%adapt%max_adapt_level > 0) call SmoothMetric( ) 


    ! ! variant of high order Riemann metric
    ! imt = 24
    ! open(imt, file=file2, status='UNKNOWN')
    ! do i=1,grid%nelem !,2
    !   elem => grid%elem(i)
    !   !if( dot_product(elem%xc(:), elem%xc(:))**0.5 < 4E-2) &
    !        call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )
    ! enddo
    ! close(imt)




  end subroutine ComputeHigherOrderEstims_OPT


  !> set HP type adaptation for the HO_recontruction
  subroutine SetHigherOrder_HP_adapt( )
    class(element), pointer :: elem
    real, dimension(:,:), allocatable :: wp ! array  with solution in vertexes
    integer :: i, j, k,  imt, imt1, is, dof
    logical :: loc_implicitly
    character(len=15) :: file1, file2
    character(len=5) :: ch5

    ! do i=1,grid%nelem
    !    elem => grid%elem(i)

    !    call Decide_hp_type_adapt(elem)
    !    !call Decide_hp_type_adapt_simple(elem)

    ! enddo


    call IsotropicMetric( )

    !variant of high order Riemann metric
    file1 = 'metrixA00000'
    file2 = 'metrixS00000'
    !file2 = 'rgabcxA00000'

    is = 0
    if(state%space%adapt%adapt_level > 0) is = int(log(1. * state%space%adapt%adapt_level)/log(10.)) 

    write( ch5, '(i5)' ) state%space%adapt%adapt_level  ! change the format if num_size /= 5 !!!
    file1(12-is: 12)  = ch5(5-is:5)
    file2(12-is: 12)  = ch5(5-is:5)


   ! imt = 24
   ! open(imt, file=file1, status='UNKNOWN')
   ! do i=1,grid%nelem
   !    elem => grid%elem(i)
   !    !    if( mod(i, 3) == 1) &
   !    !if( dot_product(elem%xc(:), elem%xc(:))**0.5 < 4E-2) &
   !         !    if(abs(elem%xc(1) - 1.5) < 0.25 .and. elem%xc(2) > 1.75 ) &
   !         call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )
      
   ! enddo
   ! close(imt)

   !  call SmoothMetric( ) 


   !  ! variant of high order Riemann metric
   !  imt = 24
   !  open(imt, file=file2, status='UNKNOWN')
   !  do i=1,grid%nelem !,2
   !    elem => grid%elem(i)
   !    !if( dot_product(elem%xc(:), elem%xc(:))**0.5 < 4E-2) &
   !         call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )
   !  enddo
   !  close(imt)


  end subroutine SetHigherOrder_HP_adapt


  !> compute the error estimates by a higher order reconstruction for one element
  subroutine ComputeHigherOrderElem_HP(elem)
    class(element), intent(inout) :: elem
    real, dimension(:,:), allocatable :: estim
    real :: norm, h_K, scale, ratio, loc_tol, glob_tol, min_dof
    real :: h_opt, lam_max, weight, fac, a, f
    integer :: k, ndimL, ideg, degP, degP1, ityp, min_idx, i, n
    logical :: iprint

    iprint = .false.
    !if(dot_product(elem%xc, elem%xc)**0.5 < 0.25) iprint = .true.
    

    ndimL = ndim

    h_K = elem%diam


   ! (ityp==1) scale = 1, (ityp==2) scale = 1./N, (ityp==3) scale = |K|/|Omega|
    !ityp = 1
    !ityp = 2  
    ityp = 3   

    scale = 1.
    if(ityp == 2) scale = 1./(1.05 * grid%nelem)**0.5
    if(ityp == 3) scale = (elem%area/state%space%domain_volume)**0.5

    loc_tol = scale *  state%space%adapt%tol_min

    !print*,'ATTENTION !!!!!!!!!!!! EDSGETSYEH'
    !loc_tol = 1

    ! temporal storing of the high order derivatives in array elem%wSS(:, :, : )
    ! first index : component of the state vector
    ! second index: order of derivatives 0 = p, 1 = p+1, 2 = p+2
    ! third index :  which partial derivative
    ! elem%wSS(k,i,j)  = \frac{ \partial^{deg+i} w^k }{ partial x^{j} partial y^{deg+i - j} }
    allocate(elem%wSS(1:ndimL, 0:2, 0:elem%deg+ 2 ) )   
    

    ! first component: 1 = error estim, 2 = N_estim, 3 = DOF_estim
    ! first component: 1 = h_opt, 2 = N_estim, 3 = DOF_estim
    ! second component = ideg
    allocate(estim(1:3, 0:2) ) 

    call  Eval_All_Derivatives(elem, ndimL) 

    k = 1 ! component of the vector

    n = 180

    min_dof = 1E+30

    do ideg = 0,2
       degP1 = elem%deg + ideg
       degP = degP1 - 1
       
       if(degP >= 0) then
       
          ! H^\mu seminorm of the HO reconstruction \mu = p_K + ideg
          norm = dot_product(elem%wSS(k, ideg, 0:degP1 ), elem%wSS(k, ideg, 0:degP1 ) )**0.5
          !norm = max_val(abs(

          !norm = 0
          do i=1,n
             a =  pi * i / n

             f = DirectionalDerivative(degP1, elem%wSS(k, ideg, 0:degP1), cos(a), sin(a))
       
            ! norm = max(f,norm)
             
          enddo


          !estim(1, ideg) = h_K**degP * norm
          !estim(2, ideg) = (estim(1, ideg) / loc_tol)**(2./degP)
          !estim(3, ideg) = (degP + 1)*(degP + 2) * estim(2, ideg)

          glob_tol = state%space%adapt%tol_min * 0.1

          estim(1, ideg) = (glob_tol**2 *(2*degP +4)/(state%space%domain_volume * norm**2))**(0.5/degP1)
          estim(2, ideg) = elem%area / estim(1, ideg)**2
          estim(3, ideg) = (degP + 1)*(degP + 2) * estim(2, ideg)
       else
          estim(3, ideg) = 1E+20
       endif

       if(estim(3, ideg) < min_dof) then
          min_dof = estim(3, ideg)
          min_idx = ideg
       endif

       if(iprint) then
          write(*,'(a6,3i5,16es12.4)') &
               'DE@WS',elem%i, ideg, degP, elem%xc, norm,loc_tol, estim(1:3, ideg)
       endif


    enddo ! ideg =0,2

    if(iprint) then
       write(*,'(a6,3i5,16es12.4)') &
            '####',elem%i, min_idx, min_idx+elem%deg -1, elem%xc, 0.,loc_tol, estim(1:3, min_idx)
       print*
    endif

    !stop

    !ratio = elem%diam / (sqrt(3.) * estim(1, ideg))

    elem%psplit = min_idx - 1

    if(elem%psplit == 1) then
       ratio = 1.
    else
       ratio = (elem%estimST/ loc_tol )**(1./(elem%deg- elem%psplit))
       ratio = min(ratio, 2.)
    endif

    h_opt = elem%diam / ratio 
    !h_opt = 2. * elem%area**0.5 / ratio 

    !write(*,'(a6,8es12.4)') 'DE#??', elem%diam,  elem%area**0.5,   elem%area**0.5/ elem%diam, &
    !      2*elem%area/ elem%diam, 2*elem%area/ elem%diam / elem%diam
    
    lam_max = 3./h_opt**2
    
    weight = 1.0
    elem%ama_p = weight * elem%psplit
    
    
    
    elem%rgabc(1) = lam_max 
    elem%rgabc(2) = 0.
    elem%rgabc(3) = lam_max

    deallocate(elem%wSS)


  end subroutine ComputeHigherOrderElem_HP



  !> compute the error estimates by a higher order reconstruction for one element
  subroutine ComputeHigherOrderElem_OPT(elem, ideg1)
    class(element), intent(inout) :: elem
    integer, intent(in) :: ideg1
    type(volume_rule), pointer ::V_rule
    real, dimension(:,:,:), allocatable :: Rhw
    real, dimension(:,:), allocatable :: RhwL
    real, dimension(:,:,:,:), allocatable :: Dw
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:,:), allocatable :: Dphi
    real, dimension(:), allocatable :: weights
    real, dimension(:,:,:), allocatable :: wi
    real, dimension(:,:), allocatable :: qq
    real, dimension(:,:), allocatable :: MM
    real, dimension(:,:), allocatable :: wExact     ! exact solution  in integ nodes
    real, dimension(:,:,:), allocatable :: DwExact  ! Der of exact solution in integ nds

    real, dimension(:,:), allocatable :: xi, Fx ! real physical coordinates
    integer :: i, j, j1,j2, k, l, degR, dof, degP, degPP, dofP, dofPP, dofPpj, Qnum, Qdof, degA,iff, iff1
    integer:: dofL, dofU, ielem, R_type, i_type, i_val, ityp, ideg
    real :: m_val, val, ratio, weight, h_opt, lam_max, lam_min, loc_tol, scale
    real :: lam_maxA, lam_minA, max_a  ! anisotropy
    real :: deref_level, max_ratio, min_ratio
    real, dimension(:,:), allocatable :: estim_neumann
    integer :: pps1, pps2, ik

    allocate(estim_neumann(1:max_eta, 1:ndim) ) 
    estim_neumann(1:max_eta, 1:ndim) = elem%eta(1:max_eta, 1:ndim)
    !print*,'   subroutine ComputeHigherOrderElem_OPT(elem) starts :', elem%i, ideg
    elem%eta(1:max_eta, 1:ndim) = 0.

    dof = elem%dof

    ! maximal possible dof
    !degP = min(elem%deg + 2, MaxDegreeImplemented)
    degP = min(elem%deg + 1, MaxDegreeImplemented)

    !degPP = elem%deg
    degPP = min(elem%deg + 1, MaxDegreeImplemented)

    dofP  = DOFtriang( degP)
    dofPP = DOFtriang( degPP)

    dofPpj = DOFtriang( elem%deg + 1)

    ! array for the reconstruction - solution of the local problem
    allocate(RhwL( 1:ndim,  1:dofPP) )

    !print*,'Before SolveLocalProblem(elem,',degP, dofPP
    !print*,'SKIPPED !!!!! '
    !call SolveLocalProblem(elem, degPP, dofPP, RhwL( 1:ndim, 1:dofPP) ) 
    !call SolveLocalProblem2(elem, degPP, dofPP, RhwL( 1:ndim, 1:dofPP) ) 
    !write(*,'(a6,300es12.4)') 'RhwL::', RhwL(:,:) 


    ! array for the the gradient recovery
    allocate(Rhw(1:ndim,  1:dofP, 0:2) )
    Rhw(:,:,:) = 0.
        
    ! increase of quadrature
    if(elem%deg == MaxDegreeImplemented ) then
       print*,'$$## alert, possible troubles in ama-hp_metric.f90'
       Qnum = elem%Qnum
    else
       Qnum = state%space%Qdeg(min(10, elem%deg + 2), 1) 
    endif

    ! quadrature rule for deg, deg+1 as well as deg+2
    V_rule => state%space%V_rule(Qnum)
    Qdof = V_rule%Qdof

    ! quadrature weights
    allocate(weights(1:Qdof) )
    call Eval_V_Weights_plus(elem, V_rule, weights(1:Qdof))

    ! basis functions in integ nodes
    phi => V_rule%phi(1:dofP, 1:Qdof)

    ! derivatives of basis functions in integ nodes
    allocate(Dphi(1:dofP, 1:nbDim, 1:Qdof) )
    call  Eval_Dphi_plus(elem, V_rule, dofP, Dphi(1:dofP, 1:nbDim, 1:Qdof) )

    ! exact solution
    allocate(wExact(1:Qdof, 1:ndim) )
    allocate(DwExact(1:Qdof, 1:ndim, 1:nbDim))

    !call SetExactSolutionQnodes(elem, V_rule, &
    !     wExact(1:Qdof, 1:ndim), DwExact(1:Qdof, 1:ndim, 1:nbDim))

    ! array for HO reconstruction
    allocate(Dw (1:ndim, 0:degP, 0:degP, 1:dofP) )   

    ! h-refinement
    !allocate(xi(1:Qdof, 1:3), Fx(1:Qdof, 1:2) )
    allocate(Fx(1:Qdof, 1:2) )

    ! projection  and its derivatives  in integ nodes
    allocate(wi(1:15, 1:ndim, 1:Qdof) )

    ! physical coordinates (for the exact solution)
    !!!!call ComputeF(elem, Qdof, xi(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )
    call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )
    
    ! exact solution in integ nodes on subelements
    call SetExactSolutionArbitraryNodes(elem, Qdof, V_rule%lambda(1:Qdof, 1:2),  &
         wExact(1:Qdof, 1:ndim), DwExact(1:Qdof, 1:ndim, 1:nbDim))

    ! plot of the exact solution  - begin
    ! allocate(MM(1:dofP, 1:dofP+3) )
    ! do k=1,dofP
    !    do j=1,dofP
    !       MM(k,j) =  dot_product( weights(1:Qdof), V_rule%phi(k, 1:Qdof)* V_rule%phi(j, 1:Qdof))
    !    enddo
    !    MM(k, dofP+1) = dot_product( weights(1:Qdof), V_rule%phi(k, 1:Qdof)* wExact(1:Qdof, 1))
    !    MM(k, dofP+2) = dot_product( weights(1:Qdof), V_rule%phi(k, 1:Qdof)* DwExact(1:Qdof, 1, 1))
    !    MM(k, dofP+3) = dot_product( weights(1:Qdof), V_rule%phi(k, 1:Qdof)* DwExact(1:Qdof, 1, 2))
    ! enddo

    ! call  SolveLocalMatrixProblem(dofP, MM(1:dofP, 1:dofP), 3, MM(1:dofP, dofP+1: dofP+3) )

    ! call PlotElem_D_Function3D(10, elem, dofP, MM( 1:dofP, dofP+1) )
    ! call PlotElem_D_Function3D(11, elem, dofP, MM( 1:dofP, dofP+2) )
    ! call PlotElem_D_Function3D(12, elem, dofP, MM( 1:dofP, dofP+3) )
    
    ! Dw(1, 1, 0, 1:dofP) =  MM( 1:dofP, dofP+1)
    ! !Dw(1, 2, 0, 1:dofP) =  MM( 1:dofP, dofP+2)
    ! !Dw(1, 3, 0, 1:dofP) =  MM( 1:dofP, dofP+3)
    ! deallocate(MM)
    ! plot of the exact solution  - end


    ! ! projection approximate solution in integ nodes on elements
    ! do k=1, ndim
    !    ik = (k-1)*dof 
    !    ! w_h
    !    wi(4, k, 1:Qdof) = matmul(elem%w(0, ik + 1: ik + dof), phi(1:dof, 1:Qdof) )
    !    ! D_x w_h,  Dy w_h
    !    wi(5, k, 1:Qdof) = matmul(elem%w(0, ik + 1: ik + dof), Dphi(1:dof, 1, 1:Qdof) )
    !    wi(6, k, 1:Qdof) = matmul(elem%w(0, ik + 1: ik + dof), Dphi(1:dof, 2, 1:Qdof) )
    ! enddo


    !if(elem%i == 1) print*,'ECHO HERED VFRTSTEGS'
    !elem%eta(:,:) = 2. * elem%eta(:,:)
    !wi(:,:,:) = 0.

    do ideg = 0, -2, -1
       if(elem%deg + ideg >= 0) then
          ! order of the HO reconstruction
          degP = elem%deg + ideg + 1    
          !degP = elem%deg + 1 
          
          if(degP > MaxDegreeImplemented) stop 'too high degP in ho-estims.f90'
          
          ! dof for HO reconstruction
          dofP  = DOFtriang( degP)

          ! projection of the solution 
          degPP = elem%deg + ideg      ! used for the projection of the approximate solution
          dofPP = DOFtriang( degPP)
          
          ! type of the stencil for the reconstruction
          !R_type = 1   ! L^2-norm
          R_type = -1  ! H^1 - norm

          if(ideg == -1) R_type = R_type * 10 ! reconstruction from \Pi^{p-1} u_h
          if(ideg == -2) R_type = R_type * 20 ! reconstruction from \Pi^{p-2} u_h
          if(ideg == +1) R_type = R_type * 50 ! reconstruction from \Pi^{p+1} u_h

          !if(elem%i <= 1) print*,'##EDE', elem%i, i_type, j2,'|',elem%face(neigh, :),'|', &
          !     R_type, degP,dofPP

          if(ideg == 0) then
             ! higher order reconstruction via weighted least square method
             call LeastSquareInterpolationWeighted(elem, ndim, .false., Qnum, degP, dofP, &
                  Dw (1:ndim, 0:degP, 0:degP, 1:dofP), R_type )

             ! graphical output - begin
             !  call PlotElem_D_Function3D(30, elem, dofP, Dw(1, 0, 0, 1:dofP) )

             !  call PlotElem_D_Function3D(20, elem, elem%dof, elem%w( 0, 1:elem%dof) )

             ! ! ! error = exact - approximate
             !  Dw(1, 1, 1, 1:dofP) = Dw(1, 1, 0, 1:dofP) 
             !  Dw(1, 1, 1, 1:elem%dof) = Dw(1, 1, 1, 1:elem%dof) -  elem%w( 0, 1:elem%dof)
             !  call PlotElem_D_Function3D(120, elem, elem%dof, Dw(1, 1, 1, 1:elem%dof)) 
             
             ! ! ! truncation = exact - HO reconstruction
             !  Dw(1, 2, 1, 1:dofP) =  Dw(1, 1, 0, 1:dofP) -  Dw(1, 0, 0, 1:dofP)
             !  call PlotElem_D_Function3D(130, elem, elem%dof, Dw(1, 2, 1, 1:elem%dof)) 

             ! ! ! truncation = approximate - HO reconstruction
             !  Dw(1, 2, 2, 1:dofP) = Dw(1, 0, 0, 1:dofP) 
             !  Dw(1, 2, 2, 1:elem%dof) = Dw(1, 2, 2, 1:elem%dof) -  elem%w( 0, 1:elem%dof)
             !  call PlotElem_D_Function3D(230, elem, elem%dof, Dw(1, 2, 2, 1:elem%dof)) 

             ! graphical output - end


           ! gradient recovery 
             !call LeastSquareInterpolationNodes(elem, Rhw(1:ndim,  1:dof, 1:2) )

             ! recovery in integ nodes on elements
             !do k=1, ndim
             !   wi(12, k, 1:Qdof) = matmul(Rhw(k, 1:dof, 1), phi(1:dof, 1:Qdof) )
             !   wi(13, k, 1:Qdof) = matmul(Rhw(k, 1:dof, 2), phi(1:dof, 1:Qdof) )
             !enddo

             ! polynomial preserving reconstruction 
             call LeastSquareInterpolationPP(elem, dofP, Rhw(1:ndim,  1:dofP, 0:2) )

             !call PlotElemFunction3D(4000 , elem, dofP, Rhw(1, 1:dofP, 0) ) 
             !call PlotElemFunction3D(4001 , elem, dofP, Rhw(1, 1:dofP, 1) ) 
             !call PlotElemFunction3D(4002 , elem, dofP, Rhw(1, 1:dofP, 2) ) 

             !call PlotElemFunction3D(20, elem, dofP, Rhw(1, 1:dofP, 0) )

             ! eastSquareInterpolation
             !!RhwL(1:ndim, 1:dofP) = Rhw(1:ndim,  1:dofP, 0)  ! BC from LeastSquareInterpolationPP
             !RhwL(1:ndim, 1:dofP) = Dw(1:ndim, 0, 0, 1:dofP) ! BC from LeastSquareInterpolationWeighted
             !call SolveLocalProblem2(elem, degP, dofP, RhwL( 1:ndim, 1:dofP) ) 

             !call PlotElemFunction3D(10, elem, dofP, RhwL(1, 1:dofP) )

             ! recovery in integ nodes on elements
             do k=1, ndim
                ! recontruction of Dw, use dof in ^^^^^^^^^^^^^^^^^
                wi(11, k, 1:Qdof) = matmul(Rhw(k, 1:dof, 0),  phi(1:dof,  1:Qdof) )
                wi(12, k, 1:Qdof) = matmul(Rhw(k, 1:dof, 1),  phi(1:dof,  1:Qdof) )
                wi(13, k, 1:Qdof) = matmul(Rhw(k, 1:dof, 2),  phi(1:dof,  1:Qdof) )
                
                !! recontruction of w, use dofP in ^^^^^^^^^^^^^^^^^
                !wi(11, k, 1:Qdof) = matmul(Rhw(k, 1:dofP, 0),  phi(1:dofP,  1:Qdof) )
                !wi(12, k, 1:Qdof) = matmul(Rhw(k, 1:dofP, 0), Dphi(1:dofP, 1, 1:Qdof) )
                !wi(13, k, 1:Qdof) = matmul(Rhw(k, 1:dofP, 0), Dphi(1:dofP, 2, 1:Qdof) )


                ! solution of the local problem 
                !wi(11, k, 1:Qdof) = matmul(RhwL( k, 1:dofP),  phi(1:dofP,    1:Qdof) )
                !wi(12, k, 1:Qdof) = matmul(RhwL( k, 1:dofP), Dphi(1:dofP, 1, 1:Qdof) )
                !wi(13, k, 1:Qdof) = matmul(RhwL( k, 1:dofP), Dphi(1:dofP, 2, 1:Qdof) )

                !write(*,'(a6,2i5,200es12.4)') 'WE..',-99, dofP, RhwL(1, k, 1:dofP)
                !print*,',,..'

             enddo

             !do l=1,Qdof
             !   write(50, *) Fx(l, 1:2), wi(11:13,1,l)
             !enddo

             !elem%wSS(1:ndim, 1:dofP, 0) = Dw(1:ndim, 0, 0, 1:dofP)    ! weighted least square
             !elem%wSS(1:ndim, 1:dofP, 0:2) = Rhw(1:ndim,  1:dofP, 0:2) ! polynomial preserving
             !elem%wSS(1:ndim, 1:dofP, 0) = RhwL(1:ndim,  1:dofP)  ! solution of the local problem
          endif

          ! the recontruction of degree elem%deg+1 is used for the reconstruction of degree elem%deg + 2

          !if(ideg == 0) then
          !   if(dofPpj /= dofP) stop 'It is correct in #ESEDE#@#E# ?'
          !   elem%wS(ndim +1 : 2 * ndim,  1:dofPpj) =  Dw (1:ndim, 0, 0, 1:dofPpj)
          !endif

          ! if(elem%i == 1) then
          !    write(*,'(a6,2i5,100es12.4)') 'ws 1',ideg,dofP,  elem%wS(1, : )
          !    !write(*,'(a6,2i5,100es12.4)') 'wS 2',ideg,dofP,  elem%wS(2, : )
          !    write(*,'(a6,2i5,100es12.4)') 'save',ideg,dofPpj,  Dw (1, 0, 0, 1:dofP) !dofPpj)
          ! endif


          ! HO reconstruction in integ nodes on elements
          do k=1, ndim
             ! w_rec
             wi(1, k, 1:Qdof) = matmul(Dw(k, 0, 0, 1:dofP), phi(1:dofP, 1:Qdof) )
             ! D_x w_rec, D_y w_rec
             wi(2, k, 1:Qdof) = matmul(Dw(k, 0, 0, 1:dofP), Dphi(1:dofP, 1, 1:Qdof) )
             wi(3, k, 1:Qdof) = matmul(Dw(k, 0, 0, 1:dofP), Dphi(1:dofP, 2, 1:Qdof) )
          enddo


          if(dofPP > 0) then
             ! projection approximate solution in integ nodes on elements
             do k=1, ndim
                ik = (k-1)*dof 
                ! w_h
                wi(4, k, 1:Qdof) = matmul(elem%w(0, ik + 1: ik + dofPP), phi(1:dofPP, 1:Qdof) )
                ! D_x w_h,  Dy w_h
                wi(5, k, 1:Qdof) = matmul(elem%w(0, ik + 1: ik + dofPP), Dphi(1:dofPP, 1, 1:Qdof) )
                wi(6, k, 1:Qdof) = matmul(elem%w(0, ik + 1: ik + dofPP), Dphi(1:dofPP, 2, 1:Qdof) )
             enddo
          else
             wi(4:6, 1:ndim, 1:Qdof) = 0.
          endif

          ! evaluation of the error estimates
          !!iff = (ideg - 1) * 4 + HO_estim_L2_p1 - 1     
          if(ideg == -2 ) iff = HO_estim_L2_p0 - 1
          if(ideg == -1 ) iff = HO_estim_L2_p1 - 1
          if(ideg ==  0 ) iff = HO_estim_L2_p2 - 1


          do k=1,ndim
             ! error estimate in the L2 norm
             elem%eta(iff+1, k) =  & !elem%eta(iff+1, k) 
                  + dot_product(weights(1:Qdof),  &
                  ( wi(1, k, 1:Qdof)- wi(4, k, 1:Qdof))**2 )

             ! error estimate in the H1 semi-norm
             elem%eta(iff+2, k) =  & ! elem%eta(iff+2, k) 
                  + dot_product(weights(1:Qdof),  &
                  ( wi(2, k, 1:Qdof)- wi(5, k, 1:Qdof))**2 + ( wi(3, k, 1:Qdof)- wi(6, k, 1:Qdof))**2 )

             !!if(elem%i == 1) write(*,'(a10,6i5)') 'HO  iff =',iff+1, iff+2

             if(ideg ==  0 ) then
                ! truncation error  in the L2 norm
                elem%eta(iff+3, k) = & ! elem%eta(iff+3, k) 
                     + dot_product(weights(1:Qdof),  &
                     ( wExact(1:Qdof, k)- wi(1, k, 1:Qdof))**2 )

                ! truncation error in the H1 semi-norm
                elem%eta(iff+4, k) = & ! elem%eta(iff+4, k) 
                     + dot_product(weights(1:Qdof),  &
                     (wi(2, k, 1:Qdof)- DwExact(1:Qdof, k,1) )**2 + (wi(3, k, 1:Qdof)- DwExact(1:Qdof,k,2))**2)

                ! recovery error in the H1 semi-norm
                elem%eta(iff+5, k) = & ! elem%eta(iff+5, k) 
                     + dot_product(weights(1:Qdof),  &
                     (DwExact(1:Qdof, k,1) - wi(12, k, 1:Qdof) )**2 + (DwExact(1:Qdof,k,2) - wi(13, k, 1:Qdof))**2)
                if(elem%eta(iff+5, k) > 1E-12) &
                     write(94,'(3es12.5, 2i5)') elem%xc(:), elem%eta(iff+5, k), elem%i,iff+5

                
                ! estimate by the recovery error in the H1 semi-norm
                elem%eta(iff+6, k) = & ! elem%eta(iff+5, k) 
                     + dot_product(weights(1:Qdof),  &
                     (wi(5,k, 1:Qdof) - wi(12, k, 1:Qdof) )**2 + (wi(6,k, 1:Qdof) - wi(13, k, 1:Qdof))**2)
                !  do l=1,Qdof
                !     write(60, *) Fx(l,1:2),  wExact(l, k   ), wi(4, k, l), wi(1, k, l), wi(11, k, l)
                ! !    write(61, *) Fx(l,1:2), DwExact(l, k, 1), wi(5, k, l), wi(2, k, l), wi(12, k, l)
                ! !    write(62, *) Fx(l,1:2), DwExact(l, k, 2), wi(6, k, l), wi(3, k, l), wi(13, k, l)

                !     write(63, *) Fx(l,1:2),  wExact(l, k   ), wi(4, k, l)-wExact(l, k),  &
                !          wi(1, k, l)-wExact(l, k), wi(11, k, l)-wExact(l, k), &
                !          abs(wi(4, k, l)-wExact(l, k)), abs(wi(11, k, l)-wExact(l, k))
                !  enddo
                !stop 'EDESEUJEU$%'
             endif

          enddo
       else

          if(ideg == -2 ) iff = HO_estim_L2_p0 - 1
          if(ideg == -1 ) iff = HO_estim_L2_p1 - 1
          if(ideg ==  0 ) iff = HO_estim_L2_p2 - 1
          elem%eta(iff + 1 : iff+2, 1:ndim) = 1E+20

       endif ! if (elem%deg + ideg >= 0) 

    enddo

    ! array for hp-parameters
    allocate(qq(-2 : 0, 1:5) )
    
    m_val = 1E+11
    
    iff1 = HO_estim_H1_p2  ! method HO reconstruction
    !iff1 = HO_trunc_H1_p2  ! method: solution local problems

    !if(elem%i == 1) print*,'REMOVe C'
    elem%estimST = elem%eta(iff1, 1)**0.5  ! HO reconstruction
    !elem%estimST =  estim_neumann(P_tot,1) ! pNeu - p_robust

    !!elem%estimST = elem%eta(HO_estim_H1_p2, 1)**0.5
    !!elem%estimST = elem%eta(HO_trunc_H1_p2, 1)**0.5
     
    !do j=-1, 1
    do j=-2, 0
       if(j == -2 ) iff = HO_estim_H1_p0 
       if(j == -1 ) iff = HO_estim_H1_p1 
       if(j ==  0 ) iff = HO_estim_H1_p2 
       !iff = iff1 + j*4

       qq(j, 1) = (elem%eta(iff, 1) / elem%eta(iff1, 1))**0.5

       !!if(elem%i == 1) write(*,'(a10,6i5)') 'adapt iff =',iff, iff1

       if(elem%deg +j  > 0 ) then
          qq(j, 2) = qq(j,1) **(2.0/(elem%deg + j))  ! originally 0.5 bad!!!
       else
          qq(j, 2) = 1E+10
       endif
       
       ! the necessary number of degrees of freedom 
       qq(j,3) = qq(j,2)**2 * DOFtriang(elem%deg + j)
       
       if(qq(j,3)  < m_val) then
          m_val = qq(j, 3)
          i_val = j
       endif
       
       
       ! if(elem%i <= -2 .or. dot_product(elem%xc(:), elem%xc(:))**0.5 < -1E-2 )  then
       ! !   if(j== - 2) write(*,'(a20, i5, 4es14.6)') '****************', elem%i, elem%xc(:), &
       ! !        dot_product(elem%xc(:), elem%xc(:))**0.5 
       !    write(*,'(a6,4i5, 5es12.4, i5)') 'regul',j, elem%deg+j, iff,iff1, elem%eta(iff, 1)**0.5, &
       !         qq(j,1), qq(j,2), qq(j,3), m_val, i_val
       !    if(j==0) print*,'----------------------------------'
       ! endif
    enddo


    ! anisotropic adaptation
    if(state%space%adapt%adapt_type == 'Ahp') then
       ! evaluation of derivatives of degree elem%deg+i_val + 1
       ! aray Dw (1:ndim, 0:degP, 0:degP, 1:dofP) contains polynomial reconstruction of degree elem%deg+1
       ! DO NOT CHANGE the previous, OTHERWISE call again  LeastSquareInterpolationWeighted
       ! restoring the values
       degP = min(elem%deg + 1, MaxDegreeImplemented)  
       dofP  = DOFtriang( degP)

       degA =  elem%deg + i_val 
       call SetElementAnisotropy(elem, degA, degP, dofP, Dw(1:ndim, 0:degP, 0:degP, 1:dofP ), &
            lam_maxA, lam_minA, max_a )

       !stop ' in EDWSEDR%$F& if(state%space%adapt%adapt_type == Ahp) then'
    endif

    !if(elem%i <= 2 .or. dot_product(elem%xc(:), elem%xc(:))**0.5 < 0.1 )  then

    ! the reduction of DOF for p-refinement is less then 2%, we keep the original p
    !if(i_val == 0 .and. qq(-1, 3) < 1.0 *  m_val) then
    !   i_val = -1 
    !endif
    !print*,'HERE', qq(-1, 3), m_val, i_val 
    !      print*
    !      print*
    !   endif
    !endif

    !if(elem%i == 1) print*,' REMOVE XX'
    !elem%eta(1:max_eta, 1:ndim) = estim_neumann(1:max_eta, 1:ndim)


    deallocate(estim_neumann)

    deref_level = 1.0 

    max_ratio = 2.0    ! maximal ratio of refinement
    min_ratio = 0.  !0.9    ! maximal ratio of DErefinement
    
    ! (ityp==1) scale = 1, (ityp==2) scale = 1./N, (ityp==3) scale = |K|/|Omega|
    !ityp = 1
    ityp = 2  
    !ityp = 3   
    
    scale = 1.
    !if(ityp == 2) scale = 1./(1.05 * grid%nelem)**0.5
    if(ityp == 2) scale = 1./(1.0 * grid%nelem)**0.5
    if(ityp == 3) scale = (elem%area/state%space%domain_volume)**0.5
    
    loc_tol = scale *  state%space%adapt%tol_min
    
    
    ! we decide p-refinement, p-kept, p-derefinement,  one-parameter variant
    elem%psplit = 0
    ratio = 1.

    if( elem%estimST > loc_tol  ) then
       ! NEW variant 
       if(i_val == 0 .and. elem%deg < MaxDegreeImplemented - 1) then !p-refinement
          elem%psplit = 1

       elseif(i_val == -1) then
          ratio = (elem%estimST/ loc_tol )**(1./elem%deg)
          !ratio = 2.
          !ratio = 1.5

       elseif(i_val == -2) then
          elem%psplit = - 1
          ratio = (elem%estimST/ loc_tol )**(1./(elem%deg-1))
          !ratio = 2.

       else !if(i_val == 0 .and. elem%deg == MaxDegreeImplemented - 1
          ratio = (elem%estimST/ loc_tol )**(1./elem%deg)
          !ratio = 1.5
       endif
       

    elseif(elem%estimST <= deref_level * loc_tol  ) then
       ! NEW VARIANT
       elem%psplit = i_val + 1
       if(elem%psplit + elem%deg > MaxDegreeImplemented - 1)  elem%psplit = 0.

       ratio = (elem%estimST/ (deref_level *loc_tol) )**(1./(elem%deg+ i_val + 1) )

    endif
    elem%estim_loc = ratio
    !ratio = ratio * 1.25
    
    
    if(elem%psplit == -1) then
       ratio = min( ratio, 5.)
    else
       ratio = min(ratio, max_ratio)
       ratio = max(ratio, min_ratio)
    endif

    elem%estimA = ratio

    h_opt = elem%diam / ratio 
    !h_opt = 2. * elem%area**0.5 / ratio 

    !write(*,'(a6,8es12.4)') 'DE#??', elem%diam,  elem%area**0.5,   elem%area**0.5/ elem%diam, &
    !      2*elem%area/ elem%diam, 2*elem%area/ elem%diam / elem%diam
    
    lam_max = 3./h_opt**2
    
    weight = 1.25
    elem%ama_p = weight * elem%psplit
    
    
    if(state%space%adapt%max_adapt_level > 0) then    
       if(state%space%adapt%adapt_type == 'Ahp') then
          lam_min = lam_max
          ratio = (lam_maxA/lam_minA)**0.5 ! different meaning prof previous
          !write(*, '(a6,6es12.4)')'EDE', lam_max, lam_min, lam_maxA/lam_minA, lam_maxA, lam_minA, max_a
          lam_max = lam_max * ratio
          lam_min = lam_min / ratio
          !write(*, '(a6,6es12.4)')'EDE', lam_max, lam_min, lam_max/ lam_min, lam_maxA, lam_minA, max_a
          !print*,'..'
          ! new metric
          elem%rgabc(1) =  lam_max * cos(max_a)**2 + lam_min * sin(max_a)**2
          elem%rgabc(2) = (lam_max - lam_min) * cos(max_a) * sin(max_a)
          elem%rgabc(3) =  lam_max * sin(max_a)**2 + lam_min * cos(max_a)**2

       else
          elem%rgabc(1) = lam_max 
          elem%rgabc(2) = 0.
          elem%rgabc(3) = lam_max
       endif
    endif


    ! draw of the reconstructed solution
    ! if(elem%i == 1 .or. elem%i == 30) then
    !    write(*,'(a6,i3,30es12.4)') 'estims', elem%i, &
    !         elem%eta(HO_estim_H1_p0, 1)**0.5, &
    !         elem%eta(HO_estim_H1_p1, 1)**0.5, &
    !         elem%eta(HO_estim_H1_p2, 1)**0.5, ratio, elem%ama_p
    
    !    !do l=1,Qdof
    !    !   write(101 + ideg, *) Fx(l, 1:2), wi(1,:,l), wi(4, 1, 1:Qdof)
    !    !enddo
    ! endif
    
  
    !endif  ! ideg == 0



    ! SMAZ - begin
    !if(elem%i == 1) then
    ! allocate(Fx(1:Qdof, 1:2) )
    ! call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )
    ! do l=1,Qdof
    !    write(100,*) Fx(l, 1:2), dot_product(elem%w(0, 1:dof), phi(1:dof, l) ) , &
    !         dot_product(Dw(ideg, 1, 0, 0, 1:dofP), phi(1:dofP, l) ) , &
    !         dot_product(elem%w(0, 1:dof), phi(1:dof, l) ) &
    !         - dot_product(Dw(ideg, 1, 0, 0, 1:dofP), phi(1:dofP, l) ) 
    ! enddo
    ! deallocate(Fx)
    !endif
    ! SMAZ -end

    deallocate(Fx,  Dphi, qq, Rhw, RhwL) 
    deallocate(weights, wExact, DwExact, Dw, wi)

  end subroutine ComputeHigherOrderElem_OPT


  !> setting of the optimal anisotropy for the given degree degA
  subroutine SetElementAnisotropy(elem, degA, degP, dofP, Dw, lam_max, lam_min, max_a )
    class(element), intent(inout) :: elem
    integer,intent(in) :: degA           ! degree of the sought anisotropy
    integer,intent(in) :: degP, dofP     ! degree of the recnstruction
    real, dimension(1:ndim, 0:degP, 0:degP, 1:dofP), intent(inout) :: Dw !derivatives: basis coeffs
    real,intent(out) :: lam_max, lam_min, max_a ! element anisotropy
    real, dimension(:,:,:,:), allocatable :: Dwi   ! all derivatives of w in integ node
    real, dimension(:,:,:), allocatable :: Dphi ! derivative of the test functions
    real, dimension(:,:), allocatable :: MassInv ! Inversion of local mass matrix of order deg+1!!
    real, dimension(:), allocatable :: ident, vec ! temporary array
    real, dimension(:,:), allocatable :: Kparams
    type(volume_rule), pointer :: V_rule
    real :: val
    integer :: degAA, dofAA, Qnum, Qdof
    integer :: i, j, k, l, ix, jx

    integer :: variant = 1   ! 1 = using reconstruction from HO, =2 using reconstruction from interpol

    if(variant == 1) then
       !print*,'ATTENTION here, degAA = degA + 1 !!!!!!!!!!!!!!!!!'
       !degAA = degA
       degAA = degA + 1
       dofAA = DOFtriang( degAA )

       !print*,'EDED',elem%i, elem%deg, degA, degAA, degP, dofP
       if(degAA > degP) stop 'troubles in SetElementAnisotropy EDRETSH'

       ! computations od derivatives of degree degAA = degA + 1
       ! increase of quadrature
       if(degP == MaxDegreeImplemented ) then
          print*,'$$ alert, possible troubles in ho-estims.f90 EDHTA'
          Qnum = elem%Qnum
       else
          Qnum = state%space%Qdeg(degP, 1)
       endif

       V_rule => state%space%V_rule(Qnum)
       Qdof = V_rule%Qdof

       !local mass matrix of order dofP x dofP
       allocate(MassInv(1:dofP, 1:dofP) )
       MassInv(:,:) = 0.

       ! identity vector
       allocate(ident(1:Qdof) )
       ident(:) = 1.

       ! evaluation of the inverse of the mass matrix
       call IntegrateBlockBBplus(elem, Qnum, dofP, ident(1:Qdof), MassInv(1:dofP,1:dofP) )
       call MblockInverse(dofP, MassInv(1:dofP, 1:dofP) )

       allocate(vec(1:dofP) )

       allocate(Dphi(1:dofP, 1:nbDim, 1:Qdof ) )
       call  Eval_Dphi_plus(elem, V_rule, dofP,  Dphi(1:dofP, 1:nbDim, 1:Qdof) )

       allocate(Dwi(1:ndim, 0:degP, 0:degP, 1:Qdof) )
       Dwi = 0.

       do k=1,ndim
          do i=1,degAA   ! i = i-th derivative
             do j=0,i  !   d^i w /dx^j dy^(i-j)

                ! evaluation of the partial derivative in integ nodes
                do l=1, Qdof  ! integ node

                   ! jx which (i-1)th derivative has to be integrate
                   ! ix according which variable: ix = 1 according x_1, ix = 2 according x_2
                   if(j==0) then
                      jx = 0
                      ix = 1
                   elseif(j==1 .and. degAA >2) then
                      jx = 0
                      ix = 2
                   else
                      jx = j-1
                      ix = 2
                   endif

                   Dwi(k, i, j, l) = dot_product(Dw(k,i-1, jx, 1:dofP), Dphi(1:dofP, ix, l) )
                enddo ! l

                ! evaluation of the basis coefficients
                call IntegrateVectorBplus(elem, Qnum, dofP, Dwi(k, i, j, 1:Qdof), vec(1:dofP) )

                Dw(k, i,j, 1:dofP) = matmul(MassInv(1:dofP,1:dofP), vec(1:dofP) )
                !MassS(:,:) = Mass(:,:)
                !call SolveLocalMatrixProblem(dofP, MassS(1:dofP, 1:dofP), 1, vec(1:dofP) )
                !Dw(k, i,j, 1:dofP) =  vec(1:dofP)
             enddo ! j
          enddo ! i
       enddo ! k


       ! write(*,*) '-----', elem%i, '  ------------------------------'
       ! k = 1
       ! do i=1,degAA   ! i = i-th derivative
       !    do j=0,i  !   d^i w /dx^j dy^(i-j)
       !       write(*,'(a6,3i5,120es12.4)') 'Dwi:',i,j,Qdof,Dwi(k, i, j, 1:Qdof)
       !    enddo
       ! enddo

       ! write(*,*) '-----------------------------------'

       ! do i=0,degAA  ! i = i-th derivative
       !    do j=0,i  !   d^i w /dx^j dy^(i-j)
       !       write(*,'(a6,3i5,120es12.4)') 'Dw:',i,j,dofP,Dw(k, i, j, 1:dofP)
       !    enddo
       ! enddo

       ! write(*,*) '########################################'


       ! degree of the function is degA, we compute derivaties of degree degA+1 = degAA
       allocate(elem%wSS(1:ndim, 1:1, 0:degA+1) )

       ! averaging of the higher-order derivatives
       do k = 1, ndim
          i = degAA   ! i = i-th derivative
          do j=0, i  !   d^i w /dx^j dy^(i-j)
             val = dot_product(V_rule%weights(1:Qdof), Dwi(k, i, j, 1:Qdof)) ! value in integ nodes
             Dwi(k, i, j, 1:Qdof) = val
             ! evaluation of the basis coefficients
             call IntegrateVectorBplus(elem, Qnum, dofP, Dwi(k, i, j, 1:Qdof), vec(1:dofP) )

             elem%wSS(k, 1, j) = dot_product(MassInv(1,1:dofP), vec(1:dofP) ) ! 1st basis coefficient

             !write(*,'(a6,3i5,120es12.4)') 'Dwi:',i,j,Qdof,elem%wSS(k, 0, j), val, Dwi(k, i, j, 1:Qdof)
          enddo

       enddo


       deallocate(Dwi, MassInv, vec, ident, Dphi)
    

    elseif(variant == 2) then

       degAA = elem%deg + 1
       allocate(elem%wSS(1:ndim, 0:2, 0:elem%deg+ 2 ) )
       call Eval_All_Derivatives(elem, ndim)

    else
       stop 'UNKNOWN variant $%DS%#G in ho-estims.f90'

    endif



    allocate( Kparams(1:1,1:4) )
    ! computing of the anisotropy
    if(state%space%adapt%Lq > -0.01) then ! L^q-nrom
       call FindAnisotropyEIp(elem, ndim, degAA, elem%wSS(1:ndim, 1, 0:degAA), &
            Kparams(1,1:4) ) !!, Kparams(ideg,2), Kparams(ideg,3), Kparams(ideg,4) )
       
    else   ! H^1-seminorm
       call FindAnisotropyEIp_H1(elem, ndim, degAA, elem%wSS(1:ndim, 1, 0:degAA), &
            Kparams(1,1:4) ) !!, Kparams(ideg,2), Kparams(ideg,3), Kparams(ideg,4) )
    endif
    
    lam_max = Kparams(1, 1)
    lam_min = Kparams(1, 2)
    max_a   = Kparams(1, 3)
    
    
    
    deallocate(elem%wSS, Kparams) 


  end subroutine SetElementAnisotropy


  !> compute the error estimates by a higher order reconstruction for one element
  subroutine ComputeHigherOrderElem(elem)
    class(element), intent(inout) :: elem
    type(volume_rule), pointer ::V_rule
    real, dimension(:,:,:,:,:), allocatable :: Dw
    !real, dimension(:,:), allocatable :: phi
    real, dimension(:), allocatable :: weights
    real, dimension(:,:,:), allocatable :: wwi
    !real, dimension(:,:,:,:), allocatable :: wi
    !real, dimension(:,:,:), allocatable :: Dphi
    real, dimension(:,:), allocatable :: xi,  loc_lam     ! exact solution  in integ nodes
    real, dimension(:,:,:), allocatable :: Fxi
    real, dimension(:,:), allocatable :: wExact     ! exact solution  in integ nodes
    real, dimension(:,:,:), allocatable :: DwExact  ! Der of exact solution in integ nds
    real, dimension(:,:), allocatable :: MM !, MMM
    real, dimension(:,:), allocatable :: proj
    real, dimension(:,:,:), allocatable :: phiR, derR
    real, dimension(:,:,:,:), allocatable :: DphiR

    real, dimension(:,:), allocatable :: Fx ! real physical coordinates
    integer :: i, j, j1,j2, k, l, degP, dofP, dofPP, Qnum, dof, Qdof, val, ideg, iff
    integer:: dofL, dofU, ielem, R_type, i_type

    !print*,'   subroutine ComputeHigherOrderElem(elem) starts :', elem%i

    dof = elem%dof

    ! increase of quadrature
    if(elem%deg == MaxDegreeImplemented ) then
       print*,'$$## alert, possible troubles in ama-hp_metric.f90'
       Qnum = elem%Qnum
    else
       Qnum = state%space%Qdeg(elem%deg + 2, 1) 
    endif

    ! quadrature rule for deg, deg+1 as well as deg+2
    V_rule => state%space%V_rule(Qnum)
    Qdof = V_rule%Qdof

    ! quadrature weights
    allocate(weights(1:Qdof) )
    call Eval_V_Weights_plus(elem, V_rule, weights(1:Qdof))

    ! maximal possible dof
    degP = min(elem%deg + 2, MaxDegreeImplemented)
    dofPP = DOFtriang( degP)

    ! basis functions in integ nodes
    !phi => V_rule%phi(1:dofPP, 1:Qdof)

    ! derivatives of basis functions in integ nodes
    !allocate(phi(1:dofPP, 1:Qdof) )
    !allocate(Dphi(1:dofPP, 1:nbDim, 1:Qdof) )
    !call  Eval_Dphi_plus(elem, V_rule, dofPP, Dphi(1:dofPP, 1:nbDim, 1:Qdof) )

    ! exact solution
    allocate(wExact(1:Qdof, 1:ndim) )
    allocate(DwExact(1:Qdof, 1:ndim, 1:nbDim))

    !call SetExactSolutionQnodes(elem, V_rule, &
    !     wExact(1:Qdof, 1:ndim), DwExact(1:Qdof, 1:ndim, 1:nbDim))

    ! array for HO reconstruction
    allocate(Dw (0:3, 1:ndim, 0:degP, 0:degP, 1:dofPP) )   

    ! solution and its derivatives  in integ nodes
    !print*,'SE#',Qdof, 3*8*ndim*Qdof
    !allocate(wi(0:2, 1:8, 1:ndim, 1:Qdof) )

    ! h-refinement
    allocate(xi(1:Qdof, 1:3), loc_lam(1:3, 1:3), Fxi(0:3, 1:Qdof, 1:2), )
    allocate(phiR(0:3,1:dofPP, 1:Qdof), DphiR(0:3, 1:dofPP, 1:nbDim, 1:Qdof), derR(1:dofPP, 1:nbDim, 1:Qdof) )
    allocate(MM(1:dofPP, 1:dofPP+ndim) )  !, MMM(1:dofPP, 1:dofPP+ndim)  )

    ! projection  and its derivatives  in integ nodes
    allocate(wwi(1:8, 1:ndim, 1:Qdof) )


    !elem%eta( HO_hp_Hp_Pm:  HO_hp_Hp_Pp , 1 :ndim  ) = 0.
    elem%eta( : , 1 :ndim  ) = 0.

    ! computation of the HO reconstructions over subelements
    weights(1:Qdof) = weights(1:Qdof) !/4  ! weights for subelements

    ! array for the projection of piecewise polynomial functions on subelements to the original element
    allocate( proj(1:dofPP, 1:ndim))
    proj(:,:) = 0.
    !do i_type = 0, 3  ! i_type = 0 center sub-elements, i_type =1,2,3  vertex subelements
    do i_type = 0, 0  ! i_type = 0 center sub-elements, i_type =1,2,3  vertex subelements

       ! integ nodes for sub-elements, not used 
       ielem = i_type
       if(ielem == 0) ielem = 4
       loc_lam(1:3, 1) = state%space%adapt%RGred(ielem, 2, 1:3)
       loc_lam(1:3, 2) = state%space%adapt%RGred(ielem, 3, 1:3)
       loc_lam(1:3, 3) = state%space%adapt%RGred(ielem, 1, 1:3)

       do l=1,Qdof
          !xi(l, 1:3) = matmul( loc_lam(1:3, 1:3), V_rule%lambda(l, 1:3) )
          xi(l, 1:3) = V_rule%lambda(l, 1:3) 
       enddo

       ! physical coordinates (only for CHECK)
       call ComputeF(elem, Qdof, xi(1:Qdof, 1:2), Fxi(i_type, 1:Qdof, 1:2) )

       !do l=1,Qdof
       !   write(120+ielem, *) xi(l, 1:2), Fxi(i_type, l, 1:2)
       !enddo

       ! exact solution in integ nodes on subelements
       call SetExactSolutionArbitraryNodes(elem, Qdof, xi(1:Qdof, 1:2),  &
            wExact(1:Qdof, 1:ndim), DwExact(1:Qdof, 1:ndim, 1:nbDim))

       ! reference basis functions in integ nodes on subelements
       call PHI_orthonormal(Qdof, nbDim, xi(1:Qdof, 1:nbDim), 3, dofPP, &
            phiR(i_type, 1:dofPP, 1:Qdof), derR(1:dofPP, 1:nbDim, 1:Qdof) )

       ! derivatives of REAL basis functions in integ nodes,  only NON_CURVED elements
       do j=1, dofPP
          DphiR(i_type,j,1,1:Qdof) = elem%F%D1F0(1,1)*derR(j, 1, 1:Qdof) + elem%F%D1F0(1,2)*derR(j, 2, 1:Qdof)
          DphiR(i_type,j,2,1:Qdof) = elem%F%D1F0(2,1)*derR(j, 1, 1:Qdof) + elem%F%D1F0(2,2)*derR(j, 2, 1:Qdof)
       enddo
       
       
       ! approximate solution in integ nodes on subelements
       do k=1, ndim
          ! w_h
          wwi(4, k, 1:Qdof) = matmul(elem%w(0, (k-1)*dof + 1: k* dof), phiR(i_type,1:dof, 1:Qdof) )
          ! D_x w_h,  Dy w_h
          wwi(5, k, 1:Qdof) = matmul(elem%w(0, (k-1)*dof + 1: k* dof), DphiR(i_type,1:dof, 1, 1:Qdof) )
          wwi(6, k, 1:Qdof) = matmul(elem%w(0, (k-1)*dof + 1: k* dof), DphiR(i_type,1:dof, 2, 1:Qdof) )
       enddo



       ! order of the reconstruction
       do ideg = 1, 2
          !!!degP = elem%deg + ideg
          degP = elem%deg + 1

          if(degP > MaxDegreeImplemented) stop 'too high degP in ho-estims.f90'

          ! dof for HO reconstruction
          dofP = DOFtriang( degP)
          
          ! type of the stencil for the reconstruction
          R_type = i_type
          j1 = mod(i_type, 3) + 1
          j2 = mod( j1 , 3) + 1
          if(i_type > 0 ) then 
             if( elem%face(neigh,i_type) < 0 .and. elem%face(neigh,j2) < 0) R_type = 0
          endif

          R_type = 0   ! neighbours reconstruction
          R_type = -1  ! elements sharing a node reconstruction


          if(ideg == 2) R_type = 10  ! reconstruction from \Pi^{p-1} u_h, neighbours 
          if(ideg == 2) R_type = -10  ! reconstruction from \Pi^{p-1} u_h,  elements sharing a node

          !if(elem%i <= 1) print*,'##EDE', elem%i, i_type, j2,'|',elem%face(neigh, :),'|', &
          !     R_type, degP,dofPP
          
          ! higher order reconstruction via least square method
          call LeastSquareInterpolationH1(elem, ndim, .false., Qnum, degP, dofP, &
               Dw (i_type, 1:ndim, 0:degP, 0:degP, 1:dofP), R_type )
          
          !call LeastSquareInterpolationL2(elem, ndim, .false., Qnum, degP, dofP, &
          !     Dw(i_type, 1:ndim, 0:degP, 0:degP, 1:dofP),  R_type  )


          ! HO reconstruction in integ nodes on (sub)-elements
          do k=1, ndim
             ! w_rec
             wwi(1, k, 1:Qdof) = matmul(Dw(i_type, k, 0, 0, 1:dofP), phiR(i_type,1:dofP, 1:Qdof) )
             ! D_x w_rec, D_y w_rec
             wwi(2, k, 1:Qdof) = matmul(Dw(i_type, k, 0, 0, 1:dofP), DphiR(i_type,1:dofP, 1, 1:Qdof) )
             wwi(3, k, 1:Qdof) = matmul(Dw(i_type, k, 0, 0, 1:dofP), DphiR(i_type,1:dofP, 2, 1:Qdof) )
          enddo
          

          ! evaluation of the error estimates
          iff = (ideg - 1) * 4 + HO_estim_L2_p1 - 1 
          do k=1,ndim
             ! error estimate in the L2 norm
             elem%eta(iff+1, k) = elem%eta(iff+1, k) + dot_product(weights(1:Qdof),  &
                  ( wwi(1, k, 1:Qdof)- wwi(4, k, 1:Qdof))**2 )
             
             ! error estimate in the H1 semi-norm
             elem%eta(iff+2, k) = elem%eta(iff+2, k) + dot_product(weights(1:Qdof),  &
                  ( wwi(2, k, 1:Qdof)- wwi(5, k, 1:Qdof))**2 + ( wwi(3, k, 1:Qdof)- wwi(6, k, 1:Qdof))**2 )
             
             ! truncation error  in the L2 norm
             elem%eta(iff+3, k) = elem%eta(iff+3, k) + dot_product(weights(1:Qdof),  &
                  ( wExact(1:Qdof, k)- wwi(1, k, 1:Qdof))**2 )
             
             ! truncation error in the H1 semi-norm
             elem%eta(iff+4, k) = elem%eta(iff+4, k) + dot_product(weights(1:Qdof),  &
                  (wwi(2, k, 1:Qdof)- DwExact(1:Qdof, k,1) )**2 + (wwi(3, k, 1:Qdof)- DwExact(1:Qdof,k,2))**2)
          enddo

       enddo ! ideg =1,2


       ! parameters for adaptations
       elem%reg = (elem%eta(HO_estim_H1_p1, 1) / elem%eta(HO_estim_H1_p2, 1))**0.5
       elem%regT0 = 2.5 * elem%diam**0.5    ! case K 1.5 ! case L 2.0
       elem%regT2 = elem%diam**(-0.5) 
       
       elem%estimST = elem%eta(HO_estim_H1_p1, 1)**0.5

       !if(elem%i == 2)  write(*,'(a6,3i5, 10es16.8)') 'REGS:', &
       !     elem%i,ideg,HO_estim_H1_p1,elem%estimST,elem%reg, elem%regT0, elem%regT2

       
       ! draw of the reconstructed solution
       if(elem%i <= -1) then
           !!write(*,'(a6,30es12.4)') 'DW', Dw(i_type, 1, 0, 0, 1:dofPP)
           do l=1,Qdof
              write(100, *) Fxi(i_type, l, 1:2), wwi(1,:,l), phiR(i_type,1:dofPP, l)
           enddo
        endif


       ! ! type of possibles adaptations

       ! ! Dw(i_type, k, 0, 0, 1:dofP) contains the p+2 order reconstruction
       ! ! H-refinement
       ! iff = HO_hp_Hp_Pm ! iff+1 =  HO_hp_H0_P0, iff+2 =  HO_hp_H0_Pp
       ! do ideg = 0, 2    ! 0 => p_K - 1 ,    1 => p_K,    2 => p_k + 1
       !     dofP = DOFtriang( elem%deg - 1 + ideg)

       !   ! HO reconstruction in integ nodes on subelements
       !    do k=1, ndim
       !       ! D_x w_rec, D_y w_rec
       !       wwi(7, k, 1:Qdof) = matmul(Dw(i_type, k, 0, 0, 1:dofP), DphiR(i_type,1:dofP, 1, 1:Qdof) )
       !       wwi(8, k, 1:Qdof) = matmul(Dw(i_type, k, 0, 0, 1:dofP), DphiR(i_type,1:dofP, 2, 1:Qdof) )
           
       !       elem%eta( iff+ideg, k) =  elem%eta( iff+ideg, k) + dot_product(weights(1:Qdof), &
       !            (wwi(2, k, 1:Qdof)-wwi(7, k, 1:Qdof))**2 +  (wwi(3, k, 1:Qdof)-wwi(8, k, 1:Qdof))**2 )

       !    enddo
       ! enddo ! ideg = 0,2


       ! ! H-fixed, projection of the reconstruction from the sub-elements to orig element
       ! ! L^2 scalar product of the reconstruction with the basis functions
       ! do k=1,ndim
       !    do i=1,dofPP
       !       proj(i,k) = proj(i,k) + dot_product(weights(1:Qdof), wwi(1, k,1:Qdof)*phiR(i_type,i,1:Qdof))
       !    enddo
       ! enddo

    enddo ! i_type=0,0(4)

    ! H -fixed, projection of the reconstruction from the sub-elements to orig element 
    
    ! iff = HO_hp_H0_Pm ! iff+1 =  HO_hp_H0_P0, iff+2 =  HO_hp_H0_Pp
    ! do ideg = 0, 2    ! 0 => p_K - 1 ,    1 => p_K,    2 => p_k + 1
    !    dofP = DOFtriang( elem%deg - 1 + ideg)

    !    ! projection of to the original element
    !    ! mass matrix
    !    MM(1:dofP,1:dofP) = elem%Mass%Mb(1:dofP, 1:dofP)
    !    ! RHS
    !    MM(1:dofP,dofP+1: dofP + ndim) = proj(1:dofP,1:ndim) 
    !    call SolveLocalMatrixProblem(dofP, MM(1:dofP, 1:dofP), ndim, MM(1:dofP, dofP+1:dofP+ndim) )

    !    ! integration over subelements
    !    do i_type = 0, 3  ! i_type = 0 center sub-elements, i_type =1,2,3  vertex subelements

    !       !print*,' HO reconstruction in integ nodes on subelements',ideg,i_type
    !       do k=1, ndim
    !          ! w_rec
    !          wwi(1, k, 1:Qdof) = matmul(Dw(i_type, k, 0, 0, 1:dofPP), phiR(i_type,1:dofPP, 1:Qdof) )
    !          ! D_x w_rec, D_y w_rec
    !          wwi(2, k, 1:Qdof) = matmul(Dw(i_type, k, 0, 0, 1:dofPP), DphiR(i_type,1:dofPP, 1, 1:Qdof) )
    !          wwi(3, k, 1:Qdof) = matmul(Dw(i_type, k, 0, 0, 1:dofPP), DphiR(i_type,1:dofPP, 2, 1:Qdof) )
    !       enddo

    !       ! L^2 -projection in integ nodes on subelements
    !       do k=1, ndim
    !          ! w_rec
    !          wwi(4, k, 1:Qdof) = matmul(MM(1:dofP, dofP+ k), phiR(i_type,1:dofP, 1:Qdof) )
    !          ! D_x w_rec, D_y w_rec
    !          wwi(5, k, 1:Qdof) = matmul(MM(1:dofP, dofP+ k), DphiR(i_type,1:dofP, 1, 1:Qdof) )
    !          wwi(6, k, 1:Qdof) = matmul(MM(1:dofP, dofP+ k), DphiR(i_type,1:dofP, 2, 1:Qdof) )

    !          elem%eta(iff+ideg, k) = elem%eta(iff+ideg, k) + dot_product(weights(1:Qdof),  &
    !               ( wwi(2, k, 1:Qdof)- wwi(5, k, 1:Qdof))**2  +&
    !               ( wwi(3, k, 1:Qdof)- wwi(6, k, 1:Qdof))**2 )
    !       enddo

    !       ! draw of the reconstructed solution
    !       if(elem%i <= 1) then
    !          do l=1,Qdof
    !             write(100*(ideg+2), *) Fxi(i_type, l, 1:2), wwi(4,:,l)
    !          enddo
    !       endif

    !    enddo ! i_type=0,3

    ! enddo ! ideg = 0,2




       ! ! subelement mass matrix
       ! do i=1,dofPP
       !    do j=i,dofPP
       !       MM(i,j) = dot_product(weights(1:Qdof), phiR(i, 1:Qdof) * phiR(j, 1:Qdof) )
       !       MM(j,i) = MM(i,j)
       !    enddo

       !    ! RHS : scalar product wwi and phiR
       !    do k=1,ndim
       !       MM(i,dofPP + k) = dot_product(weights(1:Qdof), phiR(i, 1:Qdof) * wwi(1, k, 1:Qdof) )
       !    enddo
       !    !write(*,'(a6,i5,40es12.4)') 'MM',i, MM(i,:)
       ! enddo
       
       ! if(elem%i <= 2) then
       !    do l=1,Qdof
       !       write(100+ielem, *) Fxi(l, 1:2), wwi(1,:,l)
       !       !    write(120+ielem, *) xi(l, 1:2), phiR(1:dofPP, l)
       !       !    write(130+ielem, *) xi(l, 1:2), DphiR(1:dofPP, 1,l)
       !       !    write(140+ielem, *) xi(l, 1:2), DphiR(1:dofPP, 2,l)
       !    enddo
       ! endif

       ! iff = HO_hp_Hp_Pm ! iff+1 =  HO_hp_H0_P0, iff+2 =  HO_hp_H0_Pp
       ! do ideg = 0, 2
       !    dofP = DOFtriang( elem%deg - 1 + ideg)

       !    !print*,'##EDE',ideg, dofL, dofU
       !    MMM(1:dofPP, 1:dofPP+ndim) = MM(1:dofPP, 1:dofPP+ndim)

       !    ! evaluation of the basis coefficients of the projection
       !    call SolveLocalMatrixProblem(dofP, MMM(1:dofP, 1:dofP), ndim, MMM(1:dofP, dofPP+1:dofPP+ndim) )

       !    !write(*,'(a6,100es12.4)') 'COEFS:',MMM(1:dofP, dofPP+1:dofPP+ndim)

       !    do k=1,ndim
       !       ! D \Pi w_rec_x, D \Pi_rec_y
       !       wwi(4, k, 1:Qdof) = matmul(MMM(1:dofP, dofPP + k), DphiR(1:dofP, 1, 1:Qdof) )
       !       wwi(5, k, 1:Qdof) = matmul(MMM(1:dofP, dofPP + k), DphiR(1:dofP, 2, 1:Qdof) )

       !       wwi(6, k, 1:Qdof) = matmul(MMM(1:dofP, dofPP + k), phiR(1:dofP, 1:Qdof) )

       !       if(elem%i <= 2) then
       !          !print*,'EDE#@#$#',ideg,dofP, dofPP
       !          do l=1,Qdof
       !             write(200+10*ideg, *) Fxi(l, 1:2), wwi(6,1,l)
       !          enddo
       !       endif

       !       !write(*,'(a6,100es12.4)') 'Dwi_x:', wwi(2, k, 1:Qdof)
       !       !write(*,'(a6,100es12.4)') 'Dwi_y:', wwi(3, k, 1:Qdof)

       !       ! h-kept, p-decrease, kept, increase
       !       elem%eta( iff+ideg, k) =  elem%eta( iff+ideg, k) + dot_product(weights(1:Qdof), &
       !            (wwi(2, k, 1:Qdof)-wwi(4, k, 1:Qdof))**2 +  (wwi(3, k, 1:Qdof)-wwi(5, k, 1:Qdof))**2 )

       !       !if(elem%i <= 1)print*,'EDE#@#$#',ideg,dofP, dofPP, &
       !       !     elem%eta( iff+ideg, k) , dot_product(weights(1:Qdof), &
       !       !     (wwi(2, k, 1:Qdof)-wwi(4, k, 1:Qdof))**2 +  (wwi(3, k, 1:Qdof)-wwi(5, k, 1:Qdof))**2 )

       !    enddo ! k=1,ndim

       !    ! area is one fourth, weights are for the whole triangle
       !    elem%eta( iff+ideg, 1:ndim) = elem%eta( iff+ideg, 1:ndim) / 4  
       ! enddo  ! ideg =0,2

       !if(elem%i == 1) write(*,'(a6,i5,30es12.4)') 'Hp_P*',ielem, elem%eta(HO_hp_Hp_Pm: HO_hp_Hp_Pp, 1)


    !stop 'stoped HERE WSERTG'

    deallocate(xi, Fxi, loc_lam, phiR, DphiR, derR, MM)  !, MMM)

    ! SMAZ - begin
    !if(elem%i == 1) then
    ! allocate(Fx(1:Qdof, 1:2) )
    ! call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )
    ! do l=1,Qdof
    !    write(100,*) Fx(l, 1:2), dot_product(elem%w(0, 1:dof), phi(1:dof, l) ) , &
    !         dot_product(Dw(ideg, 1, 0, 0, 1:dofP), phi(1:dofP, l) ) , &
    !         dot_product(elem%w(0, 1:dof), phi(1:dof, l) ) &
    !         - dot_product(Dw(ideg, 1, 0, 0, 1:dofP), phi(1:dofP, l) ) 
    ! enddo
    ! deallocate(Fx)
    !endif
    ! SMAZ -end

    deallocate(weights, wExact, DwExact, Dw, wwi)  ! , Dphi )

  end subroutine ComputeHigherOrderElem



  !> decide which hp adaptation will be performed
  subroutine  Decide_hp_type_adapt(elem)
     class(element), intent(inout) :: elem
    ! type(volume_rule), pointer ::V_rule
    ! real, dimension(:,:), allocatable :: errs    ! estimated error of the each case of the refinement
    ! integer, dimension(:), allocatable :: idof ! DOF  of the each case of the refinement
    ! real :: max_dec, h_opt, ratio, lam_max, scale, loc_tol, alpha
    ! integer :: deg, dof, l, k, iff, i_dec, itest, ityp

    ! itest = 0 !12

    ! if(elem%i == 1 .or.elem%i == 2 .or.elem%i == 12) & ! .or.elem%i == 20 .or. elem%i == 30 ) &
    !      itest = elem%i
    ! if(dot_product(elem%xc(1:2), elem%xc(1:2))**0.5 < 0.1) itest = elem%i

    ! allocate(errs(0:10, 1:2), idof(0:10) )
    ! k = 1  ! solution component

    ! ! actual state
    ! deg = elem%deg
    ! errs(0,1) = elem%eta(HO_estim_H1_p2, k)**0.5 
    ! idof(0) = DOFtriang( deg)

    ! ! possible cases of hp adaptation
    ! iff = HO_hp_H0_Pm - 1
    ! do l = 1, 6   ! at this moment 6 cases
    !    errs(l,1) = elem%eta(iff + l , k)**0.5

    !    select case(iff + l)

    !    case(HO_hp_H0_Pm)
    !       idof(l) = DOFtriang( deg - 1)
    !    case(HO_hp_H0_P0)
    !       idof(l) = DOFtriang( deg   )
    !    case(HO_hp_H0_Pp)
    !       idof(l) = DOFtriang( deg + 1)

    !    case(HO_hp_Hp_Pm)
    !       idof(l) = DOFtriang( deg - 1) * 4
    !    case(HO_hp_Hp_P0)
    !       idof(l) = DOFtriang( deg   )  * 4
    !    case(HO_hp_Hp_Pp)
    !       idof(l) = DOFtriang( deg + 1) * 4
    !    end select
    ! enddo

    ! !errs(0,1) = elem%eta(HO_hp_H0_P0, k)**0.5 
    ! errs(0,1) = elem%eta(HO_estim_H1_p2, k)**0.5 
    ! idof(0) = DOFtriang( deg)

    ! if(elem%i == itest) then
    !    write(*,'(a8,i5,i10,es12.4,i5,2es12.4)') 'bas', 0,idof(0), errs(0,1), elem%i, elem%xc(:)
    !    print*,'-------------------------------'
    ! endif

    ! max_dec = 0.
    ! i_dec = 0
    ! do l=1, 6
    !    if(iff+l /= HO_hp_H0_P0 .and. errs(l,1) < errs(0,1) )  then
    !       errs(l,2) = - log(errs(l,1)/ errs(0,1)) / log(1. * idof(l) / idof(0) )
    !    else
    !       errs(l,2) = -888.
    !    endif
    !    if(max_dec < errs(l,2)) then
    !       max_dec = errs(l,2)
    !       i_dec = l
    !    endif

    !    if(elem%i == itest) then
    !       if(l <= 3) then
    !          alpha = log(errs(l,1)/errs(l+3,1))/ log(2.)
    !          loc_tol = state%space%adapt%tol_min /(1.05 * grid%nelem)**0.5
             

    !          write(*,'(a8,i5, i10,2es12.4,a2,5es12.4)') 'HP',l, idof(l), errs(l,1:2), &
    !               '|',errs(l,1)*idof(l),alpha, idof(l) *(errs(l,1) / loc_tol)**(2./alpha)
    !       else
    !          write(*,'(a8,i5, i10,2es12.4,a2,5es12.4)') 'HP',l, idof(l), errs(l,1:2), &
    !               '|',errs(l,1)*idof(l)
             
    !       endif
    !       if(l==3) print*
    !    endif
    ! enddo

    ! if(elem%i == itest) then
    !    print*
    !    write(*,'(a8,i5, i10,5es12.4)') 'best',i_dec, idof(i_dec), errs(i_dec,1:2)
    !    print*,'######################################################'
    ! endif


    ! ! setting of the adaptation

    ! ! (ityp==1) scale = 1, (ityp==2) scale = 1./N, (ityp==3) scale = |K|/|Omega|
    ! !ityp = 1
    ! ityp = 2  
    ! !ityp = 3   

    ! scale = 1.
    ! if(ityp == 2) scale = 1./(1.05 * grid%nelem)**0.5
    ! if(ityp == 3) scale = (elem%area/state%space%domain_volume)**0.5

    ! ! setting of the local tolerance
    ! loc_tol = scale *  state%space%adapt%tol_min


    ! elem%psplit = 0
    ! ratio = 1.

    ! if(elem%eta(HO_estim_H1_p1, 1)**0.5  > loc_tol)  then
    !    select case(iff + i_dec)

    !    case(HO_hp_H0_Pm)
    !       elem%psplit = -1

    !    case(HO_hp_H0_P0)

    !    case(HO_hp_H0_Pp)
    !       elem%psplit = +1

    !    case(HO_hp_Hp_Pm)
    !       elem%psplit = -1
    !       ratio = 2.

    !    case(HO_hp_Hp_P0)
    !       ratio = 2.

    !    case(HO_hp_Hp_Pp)
    !       elem%psplit = +1
    !       ratio = 2.

    !    end select

    ! endif

    ! h_opt = elem%diam / ratio 

    ! lam_max = 3./h_opt**2

    ! elem%ama_p = elem%psplit

    ! elem%rgabc(1) = lam_max 
    ! elem%rgabc(2) = 0.
    ! elem%rgabc(3) = lam_max


    ! if(elem%i == itest) then
    !    do l=1, 6
    !       write(500+l, *) idof(0), errs(0,1)
    !       write(500+l, *) idof(l), errs(l,1)
    !       write(500+l,'(x)')
    !    enddo
    ! endif

    ! !stop 'stoped in  Decide_hp_type_adapt'

    ! deallocate(errs, idof)


  end subroutine Decide_hp_type_adapt


  !> decide which hp adaptation will be performed
  subroutine  Decide_hp_type_adapt_simple(elem)
    class(element), intent(inout) :: elem
    type(volume_rule), pointer ::V_rule
    real, dimension(:,:), allocatable :: errs    ! estimated error of the each case of the refinement
    integer, dimension(:), allocatable :: idof ! DOF  of the each case of the refinement
    real :: max_dec, h_opt, ratio, lam_max, scale, loc_tol, alpha, ratio_min, weight, max_ratio
    integer :: deg, dof, l, k,  i_dec, itest, ityp
    integer :: iA, iB, iC, iD, iE, iFF, iG, iH, iI

    itest = 0 !12

    if(elem%i == 1 .or.elem%i == 2 .or.elem%i == 12) & ! .or.elem%i == 20 .or. elem%i == 30 ) &
         itest = elem%i
    if(dot_product(elem%xc(1:2), elem%xc(1:2))**0.5 < 0.1) itest = elem%i



    ! setting of the adaptation

    ! (ityp==1) scale = 1, (ityp==2) scale = 1./N, (ityp==3) scale = |K|/|Omega|
    !ityp = 1
    ityp = 2  
    !ityp = 3   

    scale = 1.
    if(ityp == 2) scale = 1./(1.05 * grid%nelem)**0.5
    if(ityp == 3) scale = (elem%area/state%space%domain_volume)**0.5

    ! setting of the local tolerance
    loc_tol = scale *  state%space%adapt%tol_min

    elem%psplit = 0
    ratio = 1.

    ratio_min = 0.1

    max_ratio = 2.

    if(elem%estimST > loc_tol  ) then
             
       !if(elem%reg < par_reg .and. elem%deg < MaxDegreeImplemented - 2 ) then ! p-refinement
       if(elem%reg < elem%regT0 .and. elem%deg < MaxDegreeImplemented - 2) then !p-refinement
          elem%psplit = 1
          iA = iA + 1
          !ratio = 0.75
          !ratio = (elem%estimST/ loc_tol )**(1./(elem%deg+1) )
          
          !write(100 + 10*state%space%adapt%adapt_level + 1, *) elem%xc(:),ratio,elem%psplit
          
       else 
          ! h-refinement
          ratio = (elem%estimST/ loc_tol )**(1./elem%deg)
          !ratio = ratio ! case J
          ratio = ratio * 1.25  !1.25  !case E, J  1.5 case K
          
          !!ratio = max(ratio, 1.5)  ! case K
          
          !ratio = 2.
          
          !WSW!
          !if(Algol1) ratio = 2.
          
          iB = iB + 1
          !write(100 + 10*state%space%adapt%adapt_level + 2, *) elem%xc(:),ratio,elem%psplit, &
          !     elem%estimST/ loc_tol, (elem%estimST/ loc_tol )**(1./elem%deg)
          
       endif
       
    elseif(elem%estimST <= ratio_min * loc_tol  ) then
       
       if( elem%reg < elem%regT0 .and. elem%deg < MaxDegreeImplemented - 2) then
          elem%psplit = 1    !p-refinement
          ratio = (elem%estimST/ (ratio_min *loc_tol) )**(1./(elem%deg+1) )
          
          iC = iC + 1
          !write(100 + 10*state%space%adapt%adapt_level + 3, *) elem%xc(:),ratio,elem%psplit
          !   
          !elseif(elem%eta(P_potP, 1) <  ratio_max * loc_tol  .and.  elem%deg > 1) then 
       elseif( elem%reg > elem%regT2   .and.  elem%deg > 1) then 
          ! p-derefinement
          elem%psplit = -1 
          ratio = (elem%estimST/ (ratio_min * loc_tol) )**(1./(elem%deg-1))
          
          iD = iD + 1
          !write(100 + 10*state%space%adapt%adapt_level + 4, *) elem%xc(:),ratio,elem%psplit
          
       else
          ! h-(de)refinement
          ratio = (elem%estimST/ (ratio_min * loc_tol) )**(1./elem%deg)
          iE = iE + 1
          !write(100 + 10*state%space%adapt%adapt_level + 5, *) elem%xc(:),ratio,elem%psplit
          
       endif
       
    else
       !write(100 + 10*state%space%adapt%adapt_level + 6, *) elem%xc(:),ratio,elem%psplit
    endif

    !if( dot_product(elem%xc(1:2), elem%xc(1:2) )**0.5 < 1E-3) then
    !   write(77,'(a4,3i5,10es12.4)') 'IAMA',state%space%adapt%adapt_level, elem%i, elem%deg, &
    !        elem%estimST, loc_tol, elem%reg, ratio, min(ratio, 2.5)**1.2
    !endif
    
    ratio = min( ratio, max_ratio)
    
    h_opt = elem%diam / ratio 
    !h_opt = 2. * elem%area**0.5 / ratio 
    
    !write(*,'(a6,8es12.4)') 'DE#??', elem%diam,  elem%area**0.5,   elem%area**0.5/ elem%diam, &
    !      2*elem%area/ elem%diam, 2*elem%area/ elem%diam / elem%diam
    
    lam_max = 3./h_opt**2
    
    weight = 1.0
    elem%ama_p = weight * elem%psplit


    !write(*,'(a6,8es12.4)') '#ISO#',lam_max, lam_actual
    !lam_max = min(lam_max, lam_actual * 10)    ! at most 10^{1/2} times smaller element
    !lam_max = max(lam_max, lam_actual / 10)    ! at most 10^{1/2} times bigger element


    elem%rgabc(1) = lam_max 
    elem%rgabc(2) = 0.
    elem%rgabc(3) = lam_max
    !write(*,'(a6,i5,8es12.4)') '#ISO#',elem%i,elem%rgabc(:), lam_max, lam_actual,elem%estimST


    write(*,'(a10,10i7)') 'AdatG:',iA, iB, iC, iD, iE, iFF, iG, iH, iI
    if(state%space%adapt%adapt_level == 0) then
        write(84,'(a10,10a12)') '***:',&
             'E>T p+', ' E>T h+',' E<T p+',' E<T p-', 'E<T p0', ' ', ' ', ' ', ''
        write(84,'(a10,10a12)') '***:','iA', 'iB', 'iC','iD','iE', 'iFF', 'iG ', 'iH ', ' iI'
     endif
    
    write(84,'(a10,10i12)') 'AdatG:',iA, iB, iC, iD, iE, iFF, iG, iH, iI


    if(elem%i == itest) then
       do l=1, 6
          write(500+l, *) idof(0), errs(0,1)
          write(500+l, *) idof(l), errs(l,1)
          write(500+l,'(x)')
       enddo
    endif

    !stop 'stoped in  Decide_hp_type_adapt'

    deallocate(errs, idof)


  end subroutine Decide_hp_type_adapt_simple



 !> test for recontructions
 subroutine Reconstruction_test2D( )
   class(element), pointer :: elem
   real, dimension(:,:), pointer :: phi
   type(volume_rule), pointer :: V_rule
   real :: val, max_val1
   integer :: i,j,k, l, Qnum, Qdof, dof, dofM


   do i=1,1  ! grid%nelem
      elem => grid%elem(i)
      write(*,'(a9, 2i5,30es12.4)') 'ww:', i, elem%dof, elem%w(0,:)

      dof = elem%dof
      dofM = elem%dof - (elem%deg +1)


      do k=1, elem%Qnum

         Qnum = k
         V_rule => state%space%V_rule(Qnum)
      

         Qdof = V_rule%Qdof


         phi => V_rule%phi(1:dof, 1:Qdof)


         max_val1 = 0.
         do l=1,Qdof
            val = abs( dot_product(phi(1:dof, l), elem%w(0,1:dof) ) &
                 - dot_product(phi(1:dofM, l), elem%w(0,1:dofM) ) )
            write(10 + Qnum,*) V_rule%lambda(l, 1:2), &
                 dot_product(phi(1:dof, l), elem%w(0,1:dof) ), &
                 dot_product(phi(1:dofM, l), elem%w(0,1:dofM) ), val

            max_val1 = max(max_val1, val)
         enddo

         print*,'###',dof, k, Qnum, Qdof, max_val1
         
      enddo

   enddo


 end subroutine Reconstruction_test2D

end module Higher_order_estims
