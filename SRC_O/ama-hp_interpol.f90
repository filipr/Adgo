!> subroutines for hp method based on interpolation error
module ama_hp_interpol

  use main_data  ! contains type(mesh) :: grid for computation
  use problem_oper
  use euler_problem
  use estimates
  use lapack_oper
  use plot_geom
  use eval_sol
  use geometry
  use marking
  use ama_L2interpol

  implicit none


  public:: AnisotInterpolEstimates
  public:: SetQuantities4Metric
  public:: Eval_hp_Metric
  public:: Eval_hp_MetricElem
  public:: Set_hp_metric_hpREZ
  public:: Set_hp_metric_Inter
  public:: Eval_All_Derivatives
  !public:: LeastSquareInterpolationH1
  !public:: LeastSquareInterpolationL2
  public:: EvalSolutionVertexes
  public:: FindAnisotropyEIp
  public:: FindAnisotropyEIp_H1
  public:: FindAnisotropyEIpOLD
  public:: SmoothMetric

  public:: TriangleInsideEllipse

  public:: IsotropicMetric
contains

  !> perform the anisotropic error estimates using interpolation error estimates
  subroutine AnisotInterpolEstimates( )
    class(element), pointer :: elem
    !real, dimension(:,:), allocatable :: wp ! array  with solution in vertexes
    integer:: ndimL   ! number of quantities for metric evaluation
    integer :: i, j, k,  imt, imt1, is
    logical :: loc_implicitly
    character(len=15) :: file1, file2
    character(len=5) :: ch5

    if(nbDim /=2 ) then
       print*,' ## AnisotErrorEstimates only for nbDim = 2!'
       stop
    endif

    !allocate( wp(1:grid%npoin, 1:ndim) )
    !
    !call EvalSolutionVertexes(wp(1:grid%npoin, 1:ndim) )

    call SeekElemSupports(grid)  ! create the list of elements sharing at least a vertex with elem

    call SetQuantities4Metric(ndimL )

    !stop

    state%err(algeb) = 0.

    state%num_limits = 0
    call Eval_hp_Metric( ndimL )  !!!wp(1:grid%npoin, 1:ndim) )

    !deallocate(wp)

    file1 = 'metrixA00000'
    file2 = 'metrixS00000'

    is = 0
    if(state%space%adapt%adapt_level > 0) is = int(log(1. * state%space%adapt%adapt_level)/log(10.))

    write( ch5, '(i5)' ) state%space%adapt%adapt_level  ! change the format if num_size /= 5 !!!
    !!!write( ch5, '(i5)' ) state%time%iter  ! change the format if num_size /= 5 !!!
    file1(12-is: 12)  = ch5(5-is:5)
    file2(12-is: 12)  = ch5(5-is:5)


    ! ! variant of high order Riemann metric
    !   imt = 24
    !   open(imt, file=file1, status='UNKNOWN')
    !   do i=1,grid%nelem
    !      elem => grid%elem(i)
    !  !    if( mod(i, 3) == 1) &
    ! !     if( dot_product(elem%xc(:), elem%xc(:))**0.5 < 4E-3) &
    ! ! !    if(abs(elem%xc(1) - 1.5) < 0.25 .and. elem%xc(2) > 1.75 ) &
    !           call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )

    !   enddo
    !   close(imt)

     !print*,'STOPPED IN ESRTWGQ'
     !stop

    call SmoothMetric( )


    ! variant of high order Riemann metric
    !imt = 24
    !open(imt, file=file2, status='UNKNOWN')
    !do i=1,grid%nelem,2
    !   elem => grid%elem(i)
    !   call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )
    !enddo
    !close(imt)

    ! deallocation of arrays allocated in SeekElemSupports
    do i=1,grid%nelem
       elem => grid%elem(i)
       deallocate(elem%supp, elem%wS)
    enddo

    !print*,'stoped in NEW ama-hp_interpol.f90',  state%err(algeb)
    !stop

  end subroutine AnisotInterpolEstimates

  !> setting of quantities for evaluation of the metric in temporary array elem%wS
  subroutine SetQuantities4Metric(ndimL )
    integer, intent(inout):: ndimL   ! number of quantities for metric evaluation
    class(element), pointer :: elem
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:), allocatable :: wi, q
    integer :: i,l,k
    integer :: dof, Qnum, Qdof

    ndimL = 1
    !if(state%time_dependent) ndimL = 2

    !ndimL = 3
    !!!ndimL = grid%curved_deg
    !if(state%space%adapt%max_adapt_level == 0) ndimL = 1

    if(ndimL > 1) then
       do i=1,grid%nelem
          elem => grid%elem(i)
          allocate(elem%wSS(1:1, 1:1,1:elem%dof*ndim ) )
       enddo


       ! recomputation of the stored solution in gridS to the actual grid
       call AdvancedInterpolDGsolution(grid, gridS, .true. )
    endif

    !write(241, *) grid%nelem, ndim, state%time%ttime, state%time%tau(1), state%time%iter ! new
    !write(242, *) grid%nelem, ndim, state%time%ttime, state%time%tau(1), state%time%iter ! new

    do i=1,grid%nelem
       elem => grid%elem(i)
       allocate(elem%wS( 1:ndimL, 1:elem%dof ) )

       elem%wS(1, 1 : elem%dof) = elem%w(0, 1 : elem%dof) ! first component (or density)
       !if(ndimL > 1) elem%wS(2, 1 : elem%dof) = elem%w(0, 1 : elem%dof) ! first component (or density)
       if(ndimL > 1) elem%wS(2, 1 : elem%dof) = elem%wSS(1,1,  1 : elem%dof) !solution from the last multistep
       if(ndimL > 2) elem%wS(3, 1 : elem%dof) = 2*elem%wS(1, 1 : elem%dof) - elem%wS(2, 1 : elem%dof) ! extrapolation ahead


       if(ndimL > 1) deallocate(elem%wSS)

       !write(*,'(a8,i5,200es12.4)') '@@@@!',elem%i, elem%wS(1, 1 : elem%dof)
       !write(*,'(a8,i5,200es12.4)') '@@@@!',elem%i, elem%wS(2, 1 : elem%dof)

       ! dof = elem%dof
       ! Qnum = elem%deg
       ! Qdof = state%space%L_rule(Qnum)%Qdof
       ! allocate( q(1:ndimL, 1:Qdof) )

       ! phi => state%space%L_rule(Qnum)%phi(1:dof, 1:Qdof)
       ! do k=1,ndimL
       !    do l=1,Qdof
       !       q(k, l) = dot_product( elem%wS(k, 1 : dof), phi(1:dof,l) )
       !    enddo
       !    write(240+k,*) elem%deg,  0, q(k, 1:Qdof)
       ! enddo

       ! deallocate(q)

    enddo


    !write(241,*) 1., 1., 1., 0., 1.
    !write(242,*) 1., 1., 1., 0., 1.

  end subroutine SetQuantities4Metric

  !> evaluate the Riemann metric for hp-mesh
  subroutine Eval_hp_Metric( ndimL )    !!!wp)
    !!real, dimension(1:grid%npoin, 1:ndim), intent(inout) :: wp
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    class(element), pointer :: elem
    integer :: i
    real :: Lq, interOLD, Lqq

    !Lq =  0.   ! error in the L^{\infty}-norm
    !Lq = -1.   ! error in the H^1-seminorm
    !Lq =  2.   ! error in the L^2-norm
    Lq = state%space%adapt%Lq
    Lqq = Lq

    if(Lq <= 0.001 ) Lqq = 2.

    if(Lq >= 1.) then
       interOLD = state%err(interLq)
    elseif(Lq <= -1.) then
       interOLD = state%err(interH1)
    else
       interOLD = state%err(interL8)
    endif

    state%err(interLq) =  0.
    state%err(interL8) =  0.
    state%err(interH1) =  0.

    do i=1,grid%nelem
       elem => grid%elem(i)
       call Eval_hp_MetricElem(elem, ndimL)    !!!, wp)
       state%err(interLq) = state%err(interLq) + elem%interLq
       state%err(interH1) = state%err(interH1) + elem%interH1
       state%err(interL8) = max(state%err(interL8), elem%interL8 )

       !state%space%adapt%stop_adaptation = state%space%adapt%stop_adaptation + abs(elem%psplit)

       !if(VectorNorm(elem%xc(1:2) )  < 0.1) &
       !     write(*,*) 'elemA', elem%i, elem%psplit,  elem%ama_p

    enddo

    state%space%adapt%stop_adaptation = 0

    state%err(interH1) = state%err(interH1) ** (0.5)
    state%err(interLq) = state%err(interLq) ** (1./Lqq)

    ! the following NOT NECESSARY ?
    if(  state%space%estim_space == 'inter') then !!state%space%adapt%adapt_method == 'Ahp') then
       ! H1-seminorm
       if(Lq <= -0.99) then
          if(state%err(interH1) <= state%space%adapt%tol_max ) state%space%adapt%stop_adaptation = 1
          if(state%err(interH1) > interOLD .and. state%err(interH1) < 1.1 * interOLD) &
               state%space%adapt%stop_adaptation = -1

       ! Lq-norm
       elseif(Lq >= 1.) then
          if(state%err(interLq) <= state%space%adapt%tol_max ) state%space%adapt%stop_adaptation = 1
          if(state%err(interLq) > interOLD .and. state%err(interLq) < 1.1 * interOLD) &
               state%space%adapt%stop_adaptation = -1

       ! L^infty norm
       else
          if(state%err(interL8) <= state%space%adapt%tol_max ) state%space%adapt%stop_adaptation = 1
          if(state%err(interL8) > interOLD .and. state%err(interL8) < 1.1 * interOLD) &
               state%space%adapt%stop_adaptation = -1
       endif
    endif
    ! END of the following NOT NECESSARY ?


    !write(*,*) ' ###  end subroutine Eval_hp_Metric, state%space%adapt%stop_adaptation =', &
    !     state%space%adapt%stop_adaptation
    ! print*,
    ! print*,'Estimates of the interpolation error (# metric limiting = ', &
    !      state%num_limits ,  int(1.*state%num_limits/(ndimL *3)),':'
    ! !write(*,'(a40, 4es12.4, f4.0,i5)') ' ###  estim L2, Loo, tol, Lq OLD (stop):', &
    ! !     state%err(interLq), state%err(interL8), state%space%adapt%tol_max,&
    ! !     interOLD, state%space%adapt%Lq,  state%space%adapt%stop_adaptation

    ! write(*,'(a40, 4es12.4, f4.0,i5)') ' ###  estim L2, Loo, H1: ', &
    !      state%err(interLq), state%err(interL8),state%err(interH1)

    ! write(*,'(a40, 2es12.4, f4.0,i5)') ' ###  tol, Lq OLD (stop):', &
    !      state%space%adapt%tol_max, interOLD, state%space%adapt%Lq,  state%space%adapt%stop_adaptation

    !stop

    !print*
    !print*,'Number of metric limiting = ',  state%num_limits ,  int(1.*state%num_limits/(ndimL *3))
    !print*


  end subroutine Eval_hp_Metric

  !> evaluate the Riemann metric for hp-mesh for one element
  subroutine Eval_hp_MetricElem(elem, ndimL)    !!!, wp)
     !!real, dimension(1:grid%npoin, 1:ndim), intent(inout) :: wp
    class(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    real :: elem_esti
    integer :: i, j, k
    integer :: deg, deg1, dof, Qnum, Qdof



    if(elem%flen /= 3) then
       write(*,*)' subroutine Eval_hp_MetricElem only for triangles without HG nodes'
       stop
    endif

    ! temporal storing of the high order derivatives in array elem%wSS(:, :, : )
    ! first index : component of the state vector
    ! second index: order of derivatives 0 = p, 1 = p+1, 2 = p+2
    ! third index :  which partial derivative
    ! elem%wSS(k,i,j)  = \frac{ \partial^{deg+i} w^k }{ partial x^{j} partial y^{deg+i - j} }
    allocate(elem%wSS(1:ndimL, 0:2, 0:elem%deg+ 2 ) )

    call Eval_All_Derivatives(elem, ndimL)   !!!, wp)

    !call Set_hp_metric_hpREZ( elem )

    if(state%space%estim_space == 'inter') then ! state%space%adapt%adapt_method == 'Ahp' ) then
       call Set_hp_metric_Inter( elem, ndimL )

    else if(state%space%estim_space == 'RES') then
       call Set_hp_metric_hpREZ( elem, ndimL )

    else
       stop 'unknown estim_space in ama-hp_interpol.f90 deyde6a'
    endif


    deallocate( elem%wSS )

    !print*,'end subroutine Eval_hp_MetricElem'
  end subroutine Eval_hp_MetricElem

  !>  evaluation of metric from directional derivatives of degree elem%deg, elem%deg+1, elem%deg+2
  subroutine Set_hp_metric_Inter( elem, ndimL )
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    real :: lam_max, lam_min, max_a, lam_actual
    !real :: diam, diamT, elem_estim, est0, est0T,
    real :: h_min, h_max, lambda_min, lambda_max
    real, dimension(:,:), allocatable :: Kparams
    real :: ratio, epsilon1, ropt !, tol_loc
    integer :: ideg, degP, i, j,  pK, ip, iopt, ideg_min, ideg_max
    real :: weight
    logical :: ins

    ! approach based on the Riemann metric generated by high order interpol error estimate
    !tol_loc = state%space%adapt%tol_min

    !tol_loc = state%space%adapt%tol_min / grid%nelem**0.5

    epsilon1 = 1E+12
    h_max = state%space%diam / 4.
    !h_max = state%space%diam / 2.

    h_min = h_max / epsilon1**0.5
    lambda_min = 1./h_max**2
    lambda_max = 1./h_min**2

    elem%psplit = 0
    pK = max(1,elem%deg)

    ! Kparams: parameters of the proposed metric, first index is the degree of interpolation,
    ! second index: 1= lambda_max, 2=lambda_min, 3= max_a, 4=max_f (maximal interp error)
    ! 5=dof/area
    allocate(Kparams(0:3, 1:5) )

    ropt = 1E+50

    ! testing of the directional derivatives of degree elem%deg, elem%deg+1, elem%deg+2
    if(state%space%adapt%adapt_space == 'AMAh' .or. state%space%adapt%adapt_space == 'IMAh')  then
       ideg_min = 1
       ideg_max = 1
    else  !!! (state%space%adapt%adapt_space == 'AMAhp')
       ideg_min = 0
       ideg_max = 2
    endif

    do ideg = ideg_min, ideg_max   ! hp-AMA
    !do ideg = 1, 1   ! only h-AMA
       degP = elem%deg + ideg   ! given degree of the directional derivative

       if(degP > MaxDegreeImplemented -1) then
          Kparams(2, 1:5) = Kparams(1, 1:5)

       else
          if(state%space%adapt%Lq > -0.01) then ! L^q-nrom
             !call FindAnisotropyEIpOLD(elem, degP, elem%wSS(1:ndim, ideg, 0:degP), tol_loc, &
             !     Kparams(ideg, 1), Kparams(ideg, 2), Kparams(ideg, 3), Kparams(ideg, 4) )
             
             
             if(ideg == 1) &  ! ONLY evaluation of the interpolation error estimate in the H^1-seminorm
                  call FindAnisotropyEIp_H1(elem, ndimL, degP, elem%wSS(1:ndimL, ideg, 0:degP), &
                  Kparams(ideg,1:4) ) !!, Kparams(ideg,2), Kparams(ideg,3), Kparams(ideg,4) )


             call FindAnisotropyEIp(elem, ndimL, degP, elem%wSS(1:ndimL, ideg, 0:degP), &
                  Kparams(ideg,1:4) ) !!, Kparams(ideg,2), Kparams(ideg,3), Kparams(ideg,4) )

          else   ! H^1-seminorm
             if(degP > 1)  then ! at least P_1 approximation (degP = p+1) for H^1-seminorm
                call FindAnisotropyEIp_H1(elem, ndimL, degP, elem%wSS(1:ndimL, ideg, 0:degP), &
                     Kparams(ideg,1:4) ) !!, Kparams(ideg,2), Kparams(ideg,3), Kparams(ideg,4) )
             else
                Kparams(ideg, 1:2) = 1.5*ropt  ! => we skip this case
             endif
          endif

          !if(mod(elem%i, 50) == 1) &
          !if(elem%i == -1) &
               ! write(*,'(a4,3i5,3es16.8,a2,20es12.4)') &
               ! 'DD:',ideg, degP, DOFtriang(degP), Kparams(ideg, 1:3), '|', &
               ! 1./(Kparams(ideg, 1) * Kparams(ideg, 2) )**0.5, &
               ! DOFtriang(degP)*(Kparams(ideg, 1) * Kparams(ideg, 2))**0.5

          ! parameter of the "quality" of hp-anisotropic approaximation
          Kparams(ideg, 5) = DOFtriang(degP)  *( Kparams(ideg, 1) * Kparams(ideg, 2))**0.5
          !Kparams(ideg, 5) = DOFtriang(degP-1)*( Kparams(ideg, 1) * Kparams(ideg, 2))**0.5

          if(Kparams(ideg, 5) < ropt .and. degP > 1) then ! for degP=1 only first order derivatives
             ropt = Kparams(ideg, 5)
             iopt = ideg
          endif
       endif
    enddo  ! ideg

    !write(*,'(a4,3i5,8es16.8)') &
    !    'DD:',elem%i, degP, iopt, Kparams(iopt, 1:3), Kparams(iopt, 5), ropt


    !!! NO ADITIONAL P-ADAPTATION
    !!!! if(state%space%adapt%adapt_level >= 10) iopt = 1

    !!call TriangleInsideEllipse(elem, Kparams(1, 1:5), ins )

    ! drawing possible choices
    !call SetOptimal_hp_Metric(elem, Kparams(0:3, 1:5))

    !write(200+state%space%adapt%adapt_level, *) elem%xc,  Kparams(0:2, 4)

    !if(iopt ==0) write(300+state%space%adapt%adapt_level, *) elem%xc,  Kparams(0:2, 4)
    !if(iopt ==1) write(400+state%space%adapt%adapt_level, *) elem%xc,  Kparams(0:2, 4)
    !if(iopt ==2) write(500+state%space%adapt%adapt_level, *) elem%xc,  Kparams(0:2, 4)

    !if(mod(elem%i, 50) == 1) print*,'##', elem%i, iopt, elem%xc(:)

    ! already the MaxDegreeImplemented
    if(elem%deg + elem%psplit >= MaxDegreeImplemented -1 .and. iopt == 2) then
       iopt = 1
    endif

    ! minimal degree = 1
    if(elem%deg == 1 .and. iopt == 0) iopt = 1

    !if(VectorNorm(elem%xc(1:2) )  < 0.1) &
    !     write(*,*) 'elem', elem%i, iopt


    ! if the proposed "quality" of approximation is very close to the "quality" of actual degree
    if( iopt /= 1 .and. Kparams(iopt, 5) >= 0.95 * Kparams(1, 5)) then
    !if( iopt /= 1 .and. Kparams(iopt, 5) >= 0.70 * Kparams(1, 5)) then
       !write(*,'(a8,2i5,3es12.4)') &
       !  'opt',elem%deg, elem%deg + iopt -1, Kparams(iopt, 5),  Kparams(1, 5) , Kparams(iopt, 5)/  Kparams(1, 5)
       iopt = 1
    endif

    ! setting of the p-adaptation
    elem%psplit = iopt - 1


    ! setting the metric from the best candidate: pK-1, pK, pK+1
    !lam_max = Kparams(iopt, 1)
    !lam_min = Kparams(iopt, 2)
    !max_a   = Kparams(iopt, 3)

    ! setting the metric from the best candidate: pK-1, pK, pK+1, averaging with pK
    !weight = 0.9
    !weight = 0.75
    weight = 0.85
    !if(state%modelName == 'scalar' .or.state%modelName == '2eqs') weight = 0.5

    lam_max = weight * Kparams(iopt, 1) + (1. - weight) * Kparams(1, 1)
    lam_min = weight * Kparams(iopt, 2) + (1. - weight) * Kparams(1, 2)
    max_a   = Kparams(iopt, 3)

    if( state%space%adapt%adapt_space == 'IMAh' .or. state%space%adapt%adapt_space == 'IMAhp') then
       lam_min = lam_max
       max_a = 0.
    endif


    elem%ama_p = weight * elem%psplit

    !if(elem%deg + elem%ama_p <  1) then
    !   write(*,'(a20,2es12.4,2i6,4es12.4)') &
    !        '##D#DE#D#DSW', elem%deg + elem%ama_p, elem%ama_p, iopt, elem%deg, elem%xc
    !endif


    ! THE OPTIMAL METRIC
    !lam_max = Kparams(3, 1)
    !lam_min = Kparams(3, 2)
    !max_a =   Kparams(3, 3)
    !elem%ama_p = Kparams(3, 4)
    !!!par = Kparams(3, 5)

    !lam_actual = 1. / elem%diam**2.0
    lam_actual = (2*elem%area/elem%diam)**(-2)

    ! limitation of the aspect ratio
    ratio = lam_min / lam_max
    !ratio = max(ratio, 1E-4)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    ratio = max(ratio, 2E-5)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    !ratio = max(ratio, 2.5E-4)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    !!ratio = max(ratio, 1E-3)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    !!ratio = max(ratio, 1E-2)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    !ratio = max(ratio, diamT/ diam *0.2)   ! ratio may increase maximally 5-times
    if(state%modelName == 'scalar' .or.state%modelName == '2eqs') ratio = max(ratio, 1E-3)

    ! limitation of the refinement at one level

    !write(*,'(a6,8es12.4)') '####',lam_max, lam_actual
    lam_max = min(lam_max, lam_actual * 10)    ! at most 10^{1/2} times smaller element
    lam_max = max(lam_max, lam_actual / 10)    ! at most 10^{1/2} times bigger element


    ! new lam_min after limitation of the aspect ratio and lam_max
    lam_min = lam_max * ratio

    ! limitation from the point of view of reasonability
    lam_max = max(lambda_min, min(lambda_max, lam_max) )
    lam_min = max(lambda_min, min(lambda_max, lam_min) )


    ! new metric
    elem%rgabc(1) =  lam_max * cos(max_a)**2 + lam_min * sin(max_a)**2
    elem%rgabc(2) = (lam_max - lam_min) * cos(max_a) * sin(max_a)
    elem%rgabc(3) =  lam_max * sin(max_a)**2 + lam_min * cos(max_a)**2

    ! if(abs(elem%xc(1) - 1.) < 0.05) then
    !    write(*,'(a6,i5,3es14.6, a2, 4es14.6)') &
    !         '####',elem%i,lam_max**(-0.5), lam_min**(-0.5),  lam_min**(-0.5)/  lam_max**(-0.5), &
    !         '|',h_max, h_min
    !    !'####',elem%i,lam_max, lam_min, max_a, '|',elem%rgabc(:),  elem%ama_p
    ! endif

    !ip = elem%deg + elem%psplit
    !write(3000+state%space%adapt%adapt_level*10 + ip, *) elem%xc(:), ip

    deallocate(Kparams)
  end subroutine Set_hp_metric_Inter



  !> setting of the optimalmetric from three candiates
  subroutine SetOptimal_hp_Metric(elem, Kparams )
    type(element), intent(inout) :: elem
    real, dimension(0:3, 1:5), intent(inout) :: Kparams
    integer :: ideg, ip, degP, nstep, i, n, l
    real, dimension(1:2) :: xi
    real :: xi_tol, lam_max, lam_min, max_a, rho, t, pi, pdeg, par
    real, dimension(:,:), allocatable :: coeffs, vals
    integer, dimension(:), allocatable :: i_ang
    pi = 2 * asin(1.)

    ! graphical checking
    !xi(1:2) = (/ 0., 0.55/)
    xi(1:2) = (/ 0., 0./)
    xi_tol = 0.01
    !if(state%space%adapt%adapt_level == 2) xi_tol = 0.1

    do ideg = 0, 2
       lam_max = Kparams(ideg, 1)
       lam_min = Kparams(ideg, 2)
       max_a   = Kparams(ideg, 3)
       degP = elem%deg + ideg

       elem%rgabc(1) =  lam_max * cos(max_a)**2 + lam_min * sin(max_a)**2
       elem%rgabc(2) = (lam_max - lam_min) * cos(max_a) * sin(max_a)
       elem%rgabc(3) =  lam_max * sin(max_a)**2 + lam_min * cos(max_a)**2

       ip = 1000 + 10*(state%space%adapt%adapt_level) + ideg + 1
       !if( elem%i == -18) then
       ! if(VectorNorm(elem%xc(1:2) - xi(1:2) ) <= xi_tol )then
       !      call DrawEllips(ip, elem%rgabc(1:3), elem%xc(1:2) )

       !      write(99,'(2i5,16es12.4)') &
       !           elem%i,elem%deg+ideg-1, &
       !           -999.,lam_max, lam_min, &
       !           lam_max/lam_min, &
       !           max_a, max_a / pi , &
       !           !cos(max_a), sin(max_a), &
       !           Kparams(ideg, 5)
       !      if(ideg ==2) write(99,'(x)')
       !   endif
    enddo

    return


    ! interpolation between ideg=0 .. 2
    allocate(vals(0:2, 1:3) )    ! second index is quantity:h_max, rho, alpha
    allocate(coeffs(0:2, 1:3) )  ! second index is quantity:h_max, rho, alpha
    allocate(i_ang(1:3) )

    !vals(:,:)=0.
    !coeffs(:,:)= 0.
    do ideg = 0, 2
       vals(ideg, 1) = Kparams(ideg, 1)
       vals(ideg, 2) = Kparams(ideg, 1) / Kparams(ideg, 2)
       vals(ideg, 3) = Kparams(ideg, 3)
    enddo

    ! modification of angles
    !!if(elem%i == 10 .and. state%space%adapt%adapt_level == 1) then
    !do ideg = 0,2
    !   write(*,'(i5, 30es12.4)') ideg, vals(ideg, :)
    !enddo
    !write(*,'(a6,2i5, 20es12.4)') '##@@', 0, 0, &
    !        vals(:, 3),  maxval(vals(:, 3)) - minval(vals(:, 3)) - pi/2
    n = 1
23  continue
    if( maxval(vals(:, 3)) - minval(vals(:, 3))  > pi/2. .and. n <= 3) then
       i_ang(1:3) = minloc(vals, 1)   ! index 1, minimum in columns
       ! !!write(*,'(4i5)') -9999, i_ang(:)
       l = i_ang(3) - 1
       vals(l,3) = vals(l,3) + pi

       !write(*,'(a6, 2i5, 20es12.4)') '##@@', n, l, &
       !     vals(:, 3),  maxval(vals(:, 3)) - minval(vals(:, 3)) - pi/2

       n = n + 1
       goto 23
    endif
       !stop
    !!endif

    ! coefficients of the quadratic interpolation
    do i = 1, 3
       coeffs(0, i) = vals(0, i)
       coeffs(1, i) = -3./2*vals(0, i) + 2.*vals(1, i) - vals(2, i) /2.
       coeffs(2, i) =  1./2*vals(0, i) -    vals(1, i) + vals(2, i) /2.
    enddo

    ! seeking of the minimum, Kparams(4, :) contains the optimal metric
    Kparams(3, 1:5) = Kparams(1, 1:5)
    Kparams(3, 4) = 0.

    nstep = 10
    do i =0, nstep
       t = 2. *i / nstep
       pdeg = elem%deg + t

       if(pdeg <= MaxDegreeImplemented - 1) then
          lam_max = coeffs(0,1) + coeffs(1, 1) * t + coeffs(2, 1) * t * t
          rho     = coeffs(0,2) + coeffs(1, 2) * t + coeffs(2, 2) * t * t
          max_a   = coeffs(0,3) + coeffs(1, 3) * t + coeffs(2, 3) * t * t

          rho = max(rho, 1.)

          lam_min = lam_max / rho

          !write(*,'(a6,2i5,30es12.4)') '####',elem%i, i,t, &
          !     rho, vals(0:2, 2)

          par = DOFtriangR(pdeg - 1.) * ( lam_max * lam_min )**0.5
          if(par < Kparams(3, 5)) then
             Kparams(3, 1) = lam_max
             Kparams(3, 2) = lam_min
             Kparams(3, 3) = max_a
             Kparams(3, 4) = t - 1
             Kparams(3, 5) = par
          endif

          elem%rgabc(1) =  lam_max * cos(max_a)**2 + lam_min * sin(max_a)**2
          elem%rgabc(2) = (lam_max - lam_min) * cos(max_a) * sin(max_a)
          elem%rgabc(3) =  lam_max * sin(max_a)**2 + lam_min * cos(max_a)**2

          ip = 2000 + 100*(state%space%adapt%adapt_level) + i + 1
          !if( elem%i == -18 .and. state%space%adapt%adapt_level == 0) then
          ! print*,'###',elem%i,VectorNorm(elem%xc(1:2) - xi(1:2) ), xi_tol
          ! if(VectorNorm(elem%xc(1:2) - xi(1:2) ) <= xi_tol )then
          !    call DrawEllips(ip, elem%rgabc(1:3), elem%xc(1:2) )

          !    write(99,'(2i5,16es12.4)') &
          !         elem%i,elem%deg, pdeg, &
          !         lam_max, lam_min, &
          !         lam_max/lam_min, &
          !         max_a, max_a / pi, &
          !         !cos(max_a), sin(max_a), &
          !         par
          !    if(i == nstep) write(99,'(x)')
          ! endif
       endif

    enddo

    deallocate(coeffs, vals, i_ang)

  end subroutine SetOptimal_hp_Metric

  !>  evaluation of metric from high order derivatives of degree elem%deg+1
  !> using info from the residuall error estimates
  subroutine Set_hp_metric_hpREZ( elem, ndimL)
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    real ::  lam_max, lam_min, max_a, max_f
    real ::  diam, diamT, elem_estim, epsilon1, est0, est0T, h_min, h_max, lambda_min, lambda_max
    real :: ratio, loc_tol,  weight, scale , ratioE, h_opt, lam_actual
    real, dimension(:), allocatable :: Kparams
    integer :: ideg, degP, i, j,  pK, ityp

    ! quantity for error control, depends on the method and the problem
    if(state%time_dependent) then
       if( state%time%disc_time /= 'STDG') then
          elem_estim = (elem%estim_loc + elem%jumpsJh)**0.5     !  [ VD, MATCOM ??]
       else
          elem_estim = elem%estimST
       endif
    else
       if( state%time%disc_time /= 'STDG') then
          elem_estim = (elem%estim_loc + elem%jumpsJh)**0.5     !  [ VD, MATCOM ??]
       else
          elem_estim = elem%estimS / state%time%tau(1)
       endif
    endif


    !ityp = 1
    !ityp = 2
    ityp = 3

    scale = 1.

    allocate(Kparams(1:4) )

    ! approach based on the Riemann metric generated by high order interpol error estimate
    if(ityp == 2) scale = 1./(1.05 * grid%nelem)**0.5
    if(ityp == 3) scale = (elem%area/state%space%domain_volume)**0.5

    loc_tol = state%space%adapt%tol_min * scale
    !print*,'Verify tolerances in Set_hp_metric_hpREZ(, ama-hp_interpol.f90'

    if(state%time_dependent)  then
       if( state%time%disc_time /= 'STDG') then
          loc_tol = scale *  state%space%adapt%tol_min / state%time%FinTime**0.5
          print*,'Not tested FDG'
       else
          loc_tol = scale *  state%space%adapt%tol_min *(( state%time%ttime - state%time%ttime_save )/ state%time%FinTime)**0.5
       endif
    endif

    epsilon1 = 1E+8
    h_max = state%space%diam / 4.
    if( state%type_IC .eq. 8) h_max = 0.25  ! TRY to relax, fixed for paper stdg_est
    h_min = h_max / epsilon1**0.5
    lambda_min = 1./h_max**2
    lambda_max = 1./h_min**2

    elem%psplit = 0
    pK = max(1,elem%deg)
    elem%reg = elem%rezid / elem%area / elem%diam**(2*pK)* elem%diam**3

    elem%regT0 = elem%diam**(-2)
    elem%regT1 = elem%diam**(-3)
    elem%regT2 = elem%diam**(-4)

    !write(100+state%space%adapt%adapt_level, *) elem%xc(:), elem_estim, elem%estim_loc**0.5, elem%jumpsJh**0.5, &
    !     elem_estim  / loc_tol , &
    !     elem%reg, elem%regT0, elem%regT1, elem%regT2, elem%diam, 2*elem%area / elem%diam

    if(state%space%adapt%adapt_space == 'ANIh') then
       elem%psplit = 0
       !!lam_actual = (2.*elem%area/elem%diam)**2

    else  !state%space%adapt%adapt_space == 'ANIhp'

       !if(elem_estim > loc_tol ) then
          if(elem%reg <= elem%regT0 .and. elem%deg < MaxDegreeImplemented-1 ) then
             ! solution is regular => p-refinement
             elem%psplit = 1
             !          write(41,*) elem%xc(:), elem%reg, elem_estim, loc_tol
          !else
          !   ! solution isn't regular => h-refinement
          !
          !   if(elem%reg >= elem%regT2 .and. elem%deg > MinDegreeImplemented) then
          !      elem%psplit = -1
          !      !             write(51,*) elem%xc(:), elem%reg, elem_estim, loc_tol
          !   else
          !      !             write(61,*) elem%xc(:), elem%reg, elem_estim, loc_tol
          !   endif
          !endif

       !elseif(elem_estim < tol_loc_min ) then

          !elseif(elem%reg >= elem%regT2 .and. elem%deg > MinDegreeImplemented) then
          elseif(elem%reg >= 1.  .and. elem%deg > MinDegreeImplemented) then
             elem%psplit = -1
             !          write(71,*) elem%xc(:), elem%reg, elem_estim, loc_tol
          else
             !          write(81,*) elem%xc(:), elem%reg, elem_estim, loc_tol
          endif
       !else

          !       write(91,*) elem%xc(:), elem%reg, elem_estim, loc_tol

       !endif
    endif

    if(elem%deg + elem%psplit >= MaxDegreeImplemented -1 ) elem%psplit = 0

    ideg = elem%psplit + 1
    degP = elem%deg + ideg   ! given degree of the directional derivative

    !call FindAnisotropyEIpOLD(elem,ndimL, degP, elem%wSS(1:ndimL, ideg, 0:degP), loc_tol, &
    !     lam_max, lam_min,  max_a, max_f)

    call FindAnisotropyEIp(elem,ndimL, degP, elem%wSS(1:ndimL, ideg, 0:degP), Kparams(1:4))

    lam_max = Kparams(1)
    lam_min = Kparams(2)
    max_a   = Kparams(3)

    ratioE = (elem_estim/ loc_tol )**(1./elem%deg)
    ratioE = min(ratioE, 2.5)
    ratioE = ratioE**1.2


    !diam = elem%diam
    !est0 = 2. / diamT**2
    !est0T = 2. / elem%diam**2

    ! limitation of the aspect ratio
    ratio = lam_min / lam_max
    ratio = max(ratio, 1E-4)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    !ratio = max(ratio, diamT/ diam *0.2)   ! ratio may increase maximally 5-times

    !write(*,'(a6,i5,6es12.4)') '#@!',elem%i,ratio, ratioE, elem%diam, 2*elem%area / elem%diam, &
    !     elem%diam /( 2*elem%area / elem%diam)

    ! logical, but works?
    !!diamT = elem%diam
    diamT = 2 * elem%area / elem%diam
    h_opt = diamT / ratioE   ! opposite definition of ratioE
    h_opt = min(h_opt, h_max)
    lam_max = 1./h_opt**2
    lam_min = lam_max * ratio

    !write(99,*) elem%xc(:), elem%estimST, ratioE, diamT, h_opt

    ! works, but nonlogical
    !diamT = elem%diam
    !h_opt = diamT / ratioE
    !!lam_min = (1.5/h_opt)**2
    !!lam_max =  lam_min / ratio
    !lam_max = 3./h_opt**2
    !lam_min = lam_max * ratio

    ! isotropic metric
    !elem%psplit = 0
    !lam_max = est0 * 2**0.5
    !lam_min = est0 * 2**0.5

    !lam_max = 50 * 4.**state%space%adapt%adapt_level
    !lam_min = 50 * 4.**state%space%adapt%adapt_level

    ! a priori isotropic
    !loc_tol = loc_tol * 10
    !lam_max = ( elem_estim / loc_tol )
    !lam_max = (1./diamT *( elem_estim / loc_tol)**(1./(elem%deg+elem%psplit)))**2
    !lam_min = lam_max

    !if(elem%i ==1) write(78,*) state%space%adapt%adapt_level, loc_tol, grid%nelem
    !if(elem%face(idx,1) == 1 .or.elem%face(idx,2) == 1 .or.elem%face(idx,3) == 1 ) &
    !     write(77,'(2i5, 4es12.4, a2, 10es12.4)') &
    !     state%space%adapt%adapt_level, elem%i, elem_estim, loc_tol, &
    !     elem_estim / loc_tol, state%G_estim(7), '|', &
    !     diam, &
    !     diamT, (1./diamT *( elem_estim / loc_tol)**(1./(elem%deg+elem%psplit)))**(-1), &
    !     1./diamT**2, (1./diamT *( elem_estim / loc_tol)**(1./(elem%deg+elem%psplit)))**2

    !  keeping the maximal and the minimal eigenvalues == minimal and maximal edge
    !write(76,'(2i5, 10es12.4)') &
    !     state%space%adapt%adapt_level, elem%i, lam_max, lam_min, ratio, lambda_max, lambda_min

    lam_max = max(lambda_min, min(lambda_max, lam_max) )
    lam_min = max(lambda_min, min(lambda_max, lam_min) )

    elem%rgabc(1) =  lam_max * cos(max_a)**2 + lam_min * sin(max_a)**2
    elem%rgabc(2) = (lam_max - lam_min) * cos(max_a) * sin(max_a)
    elem%rgabc(3) =  lam_max * sin(max_a)**2 + lam_min * cos(max_a)**2

    weight = 1  !!0.75
    elem%ama_p = weight * elem%psplit

    !call DrawEllips(ip3, elem%rgabc(1:3), elem%xc(1:2) )
    !write(ip4,*) elem%xc, elem%wS(k, :)

    deallocate(Kparams)
  end subroutine Set_hp_metric_hpREZ

  !> seek the optimal anisotropy of the interpolation errorfunction
  !> \f$ E^I_p \f$, i.e.
  !> the direction having the maximal directional derivative of order degP
  !> direction is max_a, size of the derivative lam_max, size of the derivative in
  !> the perpendicular direction is lam_min
  !> paper anisot_hp
  subroutine FindAnisotropyEIp( elem, ndimL, degP, der, Kparams) !  lam_max, lam_min,  max_a, max_f)
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    integer, intent(in) :: degP  ! degree of the directional derivative
    real, dimension(1:ndimL,  0:degP), intent(in) :: der  ! partial derivatives of degree degP
    real, dimension(1:4), intent(inout) :: Kparams !  lam_max, lam_min,  max_a, max_f
    ! real, intent(out) :: lam_max, lam_min,  max_a ! characterization of the metric
    ! real, intent(out) :: max_f
    real, dimension(1:2) :: xi
    real, dimension(:,:), allocatable :: Kpar_loc
    real, dimension(:,:), allocatable :: D, R, RT, M
    real, dimension(:,:), allocatable :: Fx
    real, dimension(:), allocatable :: wi
    real, dimension(:,:), allocatable :: du_val
    real ::  lam_max, lam_min,  max_a, max_f
    real :: beta, beta_max, beta_min, min_f, tol_loc
    real :: pi, fac, f, a, rho, f_bound, f_der, scale, area
    real :: f1, f2, a1, area_min, rho_opt, max_f_opt
    integer :: n, k, j, i, ip, ifig, Qdof, du_max, itest
    real :: Lq, Lqq, rN, area1
    real :: max_u, KK, u, r_lim, h_min

    itest = 25

    ! limitation of the maximal derivative
    !max_u = 1E+12
    max_u = 1E+10
    !max_u = 1E+08
    KK = 10.


    tol_loc = state%space%adapt%tol_min
    !if(state%space%adapt%stop_adaptation == -1) tol_loc = tol_loc * 0.75
    !if(state%space%adapt%stop_adaptation == -1) tol_loc = tol_loc * 0.5
    !if(state%space%adapt%stop_adaptation == -1) tol_loc = tol_loc * 0.25

    if(( state%modelName == 'scalar' .or.state%modelName == '2eqs' ) &
         .and. state%space%adapt%stop_adaptation == -1) tol_loc = tol_loc * (0.5)**state%time%recompute_back

    !if( state%modelName == 'NSe'.and. state%space%adapt%stop_adaptation == -1) tol_loc = tol_loc * (0.75)**state%time%recompute_back

    ! if(state%space%adapt%stop_adaptation == -1.and. elem%i == 1) then
    !    print*
    !    print*,state%space%adapt%stop_adaptation
    !    print*
    !    print*,'###',  state%err(interLq), state%space%adapt%tol_max, state%space%adapt%Lq
    !    print*
    !    print*
    ! endif


    ! seeking the element size, orientation and the aspect ratio
    n = 180  ! 360 ! number of test directions
    pi = 2 * asin(1.)

    ifig = 10

    ! contains the local maxima of the directional derivatives
    !allocate(du_val(1:4*degP, 1:2) )
    !du_max = 0

    allocate(Kpar_loc(1:ndimL, 1:3) )

    fac = factorial(degP)

    !k = 1 ! from the density only
    do k=1, ndimL

       max_f = 0.
       ! we seek the direct with the maximal directional derivative
       do i=1,n
          a =  pi * i / n

          f = DirectionalDerivative(degP, der(k, 0:degP), cos(a), sin(a))
          !beta = ( tol_loc / f )**(1./ degP)

          if( f  > max_f )then
             max_f = f
             max_a = a
             !beta_max = beta
          endif

       enddo

       ! size in the perpendicular direction
       a = max_a + pi /2
       min_f = DirectionalDerivative(degP, der(k, 0:degP), cos(a), sin(a) )

       if(min_f > 0.) then
          rho = max_f / min_f  ! initial aspect ratio
       else
          rho = 1E+5
       endif


       ! limitation of the maximal derivative
!       r_lim = 1.
!       if( state%modelName == 'NSe'  .and. max_f > max_u) then
          !u = max_f
          !max_f = (KK * max_u * u - max_u *max_u) / (u + max_u * (KK-2) )

          ! Var 1
          !!r_lim = max_f / u * (0.8)
          !!if( state%space%adapt%stop_adaptation == -1) r_lim = max_f / u * (0.8)**state%time%recompute_back

          ! Var 2
          !r_lim = max_f / u * (0.8)**state%time%recompute_back

!          r_lim = 0.


!          state%num_limits = state%num_limits + 1
          !write(*,'(a5,3i5,10es12.4)') 'WEW',elem%i, degP, k, u, max_f, rho, r_lim
!       endif

       !if(max_f > 1E+10) then
       !   write(*,'(a5,3i5,10es12.4)') 'WEW',elem%i, degP, k, max_f, max_a
       !   write(54,*) elem%xc(:), degP, k, max_f, 'WEW'
       !endif



       !write(*,'(a8,3i5,20es12.4)') '####@@!',elem%i, degP, k, max_f, min_f, rho

       goto 10 ! a simplification of the algorithm from APNUM

       ! setting of the upper estimate of type
       ! DirectDeriv M= max_f ( x^T Q D Q^T x )**(degP / 2)
       ! D = diag (1, rho**(-2/degP) )

       allocate(D(1:2, 1:2), R(1:2, 1:2), RT(1:2, 1:2), M(1:2, 1:2) )

       R(1,1) = cos(max_a)
       R(1,2) = -sin(max_a)
       R(2,1) = sin(max_a)
       R(2,2) = cos(max_a)

       RT(1:2, 1:2) = transpose(R(1:2, 1:2) )

       D(2,1) = 0.
       D(1,2) = 0.

       f1 = max_f
       f2 = rho

       D(1,1) = 1.
       D(2,2) = rho ** (-2./degP)
       M(1:2, 1:2) = matmul(R(1:2, 1:2), matmul(D(1:2, 1:2), RT(1:2, 1:2)))

       !if(elem%i == itest) then
       !   ifig = ifig + 1
       !   call DrawEstimate(ifig, elem%xc(1:2), max_f, degP, M(1:2, 1:2) )
       !endif

       !area = (rho / max_f**2)**(-1./degP)
       !write(*,'(a6, 3i5,8es12.4)') &
       !     '???',elem%i, 0, 0, max_f, min_f, rho, area

       area_min = 1E+30

       ! seeking of optimal values f_max and rho satisfying the estimate
       ! for several values of f_max, we seek the minimal rho
       do j = 0, 10
          max_f = f1 * (1.01**j)
          !rho = f2
          rho = 1E+6

          D(2,2) = rho ** (-2./degP)
          M(1:2, 1:2) = matmul(R(1:2, 1:2), matmul(D(1:2, 1:2), RT(1:2, 1:2)))

          do i=1,n
             !a = 2 * pi * i / n
             a =  pi * i / n
             xi(1) = cos(a)
             xi(2) = sin(a)

             ! directional derivative in direction a
             f_der = DirectionalDerivative(degP, der(k, 0:degP), cos(a), sin(a) )

             !f  = M(1,1) * cos(a)**2 + ( M(1,2) + M(2,1) )*cos(a) *sin(a) &
             !     + M(2,2)*sin(a)**2
             f = dot_product(xi(1:2), matmul(M(1:2, 1:2), xi(1:2) ) )

             f_bound =  max_f * abs(f)**(degP/2.)

             if(f_der > 1.001 * f_bound .and. a /=max_a) then
                ! estimate is violated,
                ! we have to modify the anisotropy of the interpolation error function

                !! decrease of the ratio of the interpolation error function
                scale = ((f_der / max_f)**(2./degP) - (cos(a - max_a))**2) &
                     /(sin(a - max_a))**2
                scale = scale**(degP/2.)
                min_f = max_f * scale

                !! increase of the size of the interpolation error function
                !scale = f_der / f_bound
                !max_f = max_f * scale
                !min_f = min_f * scale

                ! refreshing of the matrix on the right-hand-side
                rho = max_f / min_f
                D(2,2) = rho ** (-2./degP)
                M(1:2, 1:2) = matmul(R(1:2, 1:2), matmul(D(1:2, 1:2), RT(1:2, 1:2)))

             endif
          enddo  ! i=1,n

          area = (rho / max_f**2)**(-1./degP)
          !write(*,'(a6, 3i5,8es12.4)') &
          !     '???',elem%i, j, i, max_f, min_f, rho, area

          if(area < area_min) then
             area_min = area
             max_f_opt = max_f
             rho_opt = rho
          else
             ! we expect that the dependence of area on max_f is strictly convex
             ! so when are starts to increase, the optimal value was found
             goto 20
          endif

          !if(elem%i == itest)  then
          !   ifig = ifig + 1
          !   call DrawEstimate(ifig, elem%xc(1:2), max_f, degP, M(1:2, 1:2) )
          !endif

       enddo   ! j=1,5

       deallocate(D, R, RT, M)

20     continue
       ! the optimal anisotropy of interpolation error function found
       max_f = max_f_opt
       rho = rho_opt
       min_f = max_f / rho
       area = (rho / max_f**2)**(-1./degP)

       !write(*,'(a6, 3i5,8es12.4)') &
       !     '???',elem%i, -99, -99, max_f, min_f, rho, area
       !write(*,*) '(__________________________________)'

10     continue

       !if(elem%i == itest) stop
       Lq = state%space%adapt%Lq
       if(Lq >= 1.) then  ! q < \infty
          ! balance || e||_K <= tol (|K|/|Omega|)^(1/q)
          area = (tol_loc * rho**0.5 / ( c_p_q (degP, Lq) * max_f))**(2./degP)
          area = area * state%space%domain_volume**(-2./(Lq * degP) )
          ! !! area1 = area

          ! balance || e||_K <= tol (1 / N))^(1/q)
          !rN = grid%nelem**(1./Lq)
          !area = tol_loc * rho**0.5 / ( rN * c_p_q (degP, Lq) * max_f)
          !area = area**(2.*Lq /(Lq * degP + 2) )

          !write(*,'(a6,6es14.6)') '@@@@@',area1, area, elem%area
       else
          area = (tol_loc * rho**0.5 /  max_f)**(2./degP)
       endif

       !write(*,'(a5,2i5,8es12.4)')  '!!!', elem%i, degP, &
       !     tol_loc, rho**0.5 , c_p_q (degP, Lq),  max_f, (2./degP)

       ! ????
       !beta_min = (rho**(1./degP) * area/pi)**0.5
       beta_min = (rho**(1./degP) * area/1.)**0.5
       beta_max = beta_min / rho**(1./degP)

       ! new tests of limitation
       r_lim = 1.
       if( state%modelName == 'NSe' ) then
          !h_min = 5E-4
          !h_min = 1E-3
          h_min = state%model%Re1 * degP
          !h_min = 5E-4 * degP  ! should be the same as ^^^ (state%model%Re1 * degP )

          if(state%type_IC == 8) h_min = 0.005
          !!!if(beta_max < 5E-3) print*,'####', beta_max , h_min

          if( beta_max < h_min) then
             !write(*,'(a4,4i5,8es12.4)') '###',elem%i,degP,k,  state%num_limits, beta_min, beta_max,&
             !     beta_min/ beta_max, max_f
             beta_max = h_min
             beta_min =   beta_max * rho**(1./degP)
             r_lim = 0.
             state%num_limits = state%num_limits + 1

          endif
       endif

       !beta_min = max(beta_min, 1E-3)
       !beta_max = max(beta_max, 1E-3)

       !if(elem%i == itest) then
       !    print*,'************************************'
       !    print*,'***', beta_max, beta_min, max_f
       !    print*,'***', (tol_loc / max_f )**(1./ degP), &
       !         (tol_loc / min_f )**(1./ degP)
       !write(192,'(a5,2i5,12es10.2)') '???', elem%i, degP, &
       !     beta_max, beta_min, max_f, min_f, max_f / min_f,  beta_max/  beta_min
       !        lam_max, lam_min, max_a, a
       !endif

       !beta_min = ( tol_loc / min_f )**(1./ degP)

       lam_max = 1./beta_max**2
       lam_min = 1./beta_min**2

       Kpar_loc(k,1) = lam_max
       Kpar_loc(k,2) = lam_min
       Kpar_loc(k,3) = max_a

       if(elem%i == -1) &
          write(*,'(a6,3i5,3es14.6, a2,20es12.4)') 'WW',elem%i,degP,k,Kpar_loc(k,1:3),'|' !,  der(k, 0:degP)
    enddo  ! end of k


    ! metric compositions
    if(ndimL > 1) then

       ! SMAZ the following
       !Kparams(1) = lam_max
       !Kparams(2) = lam_min
       !Kparams(3) = max_a
       Kparams(4) = max_f


       call MetricComposition(ndimL, Kpar_loc(1:ndimL, 1:3),  Kparams(1:3), elem%xc(1:2), elem%i, degP )

    else
       Kparams(1) = lam_max
       Kparams(2) = lam_min
       Kparams(3) = max_a
       Kparams(4) = max_f

    endif

     if(elem%i == -1) then
        write(*,'(a6,3i5,3es14.6, a2,20es12.4)') 'QQ',elem%i,degP,0,Kparams(1:3),'|' !,  der(k, 0:degP)
        print*,'#########################################################'
     endif
     !!if(degP == elem%deg+ 2) stop


    ! estimates of the interpolation error in the L^{\infty} norm
    if(degP == elem%deg + 1) then

       !do k=1, ndim
       do k=1, 1 ! from the first component only (the actual solution)
          elem%interL8 = 0.

          if(Lq >= 1.) then
             Lqq = Lq
          else
             Lqq = 2.
          endif

          Qdof = elem%Qdof

          allocate( Fx(1:Qdof, 1:nbDim), wi(1:Qdof) )
          !integration nodes on K
          call ComputeF(elem, Qdof, state%space%V_rule(elem%Qnum)%lambda(1:Qdof,1:nbDim), &
               Fx(1:Qdof, 1:nbDim) )

          ! interpolation error function in integ. nodes
          do i=1,Qdof
             xi(1:2) =Fx(i, 1:2) - elem%xc(1:2)
             wi(i) = DirectionalDerivative(degP, der(k, 0:degP), xi(1),xi(2) )

             ! limitation of the maximal derivative
             wi(i) = wi(i) * r_lim

             ! q-th power
             wi(i) = (abs(wi(i)) )**Lqq
          enddo
          ! error  in the Lq norm
          call IntegrateFunction(elem, wi, elem%interLq)

          ! NOT POWERED !!!
          !elem%interLq = elem%interLq**(1./Lqq)

          deallocate(Fx, wi)

          ! interpolation error function estimate in L^\infty norm
          do i=1,3  ! vertexes of triangle
             xi(1:2) = grid%x(elem%face(idx, i), 1:2) - elem%xc(1:2)

             f = DirectionalDerivative(degP, der(k, 0:degP), xi(1),xi(2) )

             ! limitation of the maximal derivative
             f = f * r_lim

             elem%interL8 = max ( elem%interL8, f )
             !write(*,'(a6,3i5,6es12.4)') &
             !     '??',elem%i, i, elem%face(idx, i), xi(1:2), abs(f),elem%interL8
          enddo ! i

          ! size of the non-counted area
          if(r_lim == 0.)  state%err(algeb) = state%err(algeb) + elem%area

       enddo  ! k
    endif

    !deallocate(du_val)
    deallocate(Kpar_loc)

  end subroutine FindAnisotropyEIp

  !> seek the optimal anisotropy of the gradient of the interpolation error function
  !> \f$ |\nabla E^I_p| \f$, i.e.
  !> the direction having the maximal directional derivative of order degP
  !> direction is max_a, size of the derivative lam_max, size of the derivative in 
  !> the perpendicular direction is lam_min
  !> paper anisot_hp
  subroutine FindAnisotropyEIp_H1( elem, ndimL, degP, der, Kparams) !  lam_max, lam_min,  max_a, max_f)
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    integer, intent(in) :: degP  ! degree of the directional derivative
    real, dimension(1:ndimL,  0:degP), intent(in) :: der  ! partial derivatives of degree degP
    real, dimension(1:4), intent(inout) :: Kparams !  lam_max, lam_min,  max_a, max_f
    ! real, intent(out) :: lam_max, lam_min,  max_a ! characterization of the metric
    ! real, intent(out) :: max_f
    real, dimension(1:2) :: xi
    integer :: degPP
    real, dimension(:,:), allocatable :: derPP, alpha
    real, dimension(:,:), allocatable :: Kpar_loc
    real, dimension(:,:), allocatable :: D, R, RT, M
    real, dimension(:,:), allocatable :: Fx
    real, dimension(:), allocatable :: wi
    real, dimension(:,:), allocatable :: du_val
    real ::  lam_max, lam_min,  max_a, max_f
    real :: beta, beta_max, beta_min, min_f, tol_loc
    real :: pi, fac, f, a, rho, f_bound, f_der, scale, area
    real :: f1, f2, a1, area_min, rho_opt, max_f_opt
    integer :: n, k, j, i, ip, ifig, Qdof, du_max, itest, l, ll, deg0
    real :: Lq, Lqq, rN, area1
    real :: max_u, KK, u, r_lim, h_min

    r_lim = 1.
    itest = 25

    ! minmal degree = 1
    if(degP <= 1) then
       print*,'   ... TROUBLE in FindAnisotropyEIp_H1', degP-1
       stop
    endif


    ! computations of the coefficients of the magnitude of the gradient of interp. function
    deg0  = degP - 1
    degPP = 2*deg0
    allocate(derPP(1:ndimL, 0:degPP) )
    allocate(alpha(1:8, 0:degPP) )
    alpha(:,:) = 0.

    do k=1, ndimL
       ! paper AM15 - coefficients alpha_l, l=0,...,degP, OPPOSITE ORDERING
       do l=0,degP
          alpha(1, l) = der(k, degP-l)  
       enddo

       do l=0, degP
          alpha(1, l) = alpha(1,l) / factorial(l) / factorial(degP - l) 
          !alpha(1,l) = 1. + l   ! TEST ONLY
       enddo

       ! paper AM15 - coefficients beta^{(1,2)}_l, l=0,...,degP - 1
       do l=0, degP-1
          alpha(2, l) = (l + 1) * alpha(1, l+1)
          alpha(3, l) = (degP - l) * alpha(1, l)
       enddo

       ! paper AM15 - coefficients gamma^{(1,2)}_l, l=0,...,degPP
       do l=0, degP-1
          do ll = 0, l
             alpha(4, l)         = alpha(4, l)         + alpha(2,ll)           * alpha(2,l - ll)
             alpha(5, l)         = alpha(5, l)         + alpha(3,ll)           * alpha(3,l - ll)

             if(degPP-l > l) then
                alpha(4, degPP-l) = alpha(4, degPP-l) + alpha(2,degP -1 - ll) * alpha(2,degP-1 -l + ll)
                alpha(5, degPP-l) = alpha(5, degPP-l) + alpha(3,degP -1 - ll) * alpha(3,degP-1 -l + ll)
             endif

             !write(*,'(a6,2i5,a3,200es12.4)') '???:',l,ll,'|',alpha(4, :)
          enddo

       enddo
       ! paper AM15 - final coefficients, back to ORIGINAL ORDERING
       do l=0, degPP
          derPP(k, degPP-l) = alpha(4, l) + alpha(5, l)
       enddo

    enddo

    !print*,'------------------------'
    !k = 1
    !write(*,'(a6,4i6)') 'params',degP-1, degP, degPP
    !write(*,'(a6,200es12.4)') 'DER:', der(k, :)
    !write(*,'(a6,200es12.4)') 'a:', alpha(1, :)
    !write(*,'(a6,200es12.4)') 'b1:',alpha(2, :)
    !write(*,'(a6,200es12.4)') 'b2:',alpha(3, :)
    !write(*,'(a6,200es12.4)') 'c1:',alpha(4, :)
    !write(*,'(a6,200es12.4)') 'c2:',alpha(5, :)
    !write(*,'(a6,200es12.4)') '|v|^2:',derPP(k, 0:degPP)
    !stop

    deallocate(alpha)

    if(state%space%adapt%Lq > -0.01) goto 300

    ! EVALUATION of the appropriate metrix
    
    ! limitation of the maximal derivative
    !max_u = 1E+12
    max_u = 1E+10
    !max_u = 1E+08
    KK = 10.


    tol_loc = state%space%adapt%tol_min
    !if(state%stop_adaptation == -1) tol_loc = tol_loc * 0.75
    !if(state%stop_adaptation == -1) tol_loc = tol_loc * 0.5
    !if(state%stop_adaptation == -1) tol_loc = tol_loc * 0.25

    if(( state%modelName == 'scalar' .or.state%modelName == '2eqs' ) &
         .and. state%space%adapt%stop_adaptation == -1)  &
         tol_loc = tol_loc * (0.5)**state%time%recompute_back

    !if( state%modelName == 'NSe'.and. state%stop_adaptation == -1) tol_loc = tol_loc * (0.75)**state%recompute_back

    ! seeking the element size, orientation and the aspect ratio
    n = 180  ! 360 ! number of test directions
    pi = 2 * asin(1.) 

    ifig = 10

    ! contains the local maxima of the directional derivatives
    !allocate(du_val(1:4*degP, 1:2) )
    !du_max = 0

    allocate(Kpar_loc(1:ndimL, 1:3) )

    fac = factorial(degP)

    !k = 1 ! from the density only
    do k=1, ndimL

       max_f = 0.
       ! we seek the direct with the maximal directional derivative
       do i=1,n
          a =  pi * i / n

          !f = DirectionalDerivative(degP, der(k, 0:degP), cos(a), sin(a))
          f = DirectionalDerivative_H1(degPP, derPP(k, 0:degPP), cos(a), sin(a))
          !beta = ( tol_loc / f )**(1./ degP) 
       
          if( f  > max_f )then
             max_f = f
             max_a = a
             !beta_max = beta
          endif

       enddo

       ! size in the perpendicular direction
       a = max_a + pi /2
       !min_f = DirectionalDerivative(degP, der(k, 0:degP), cos(a), sin(a) )
       min_f = DirectionalDerivative_H1(degPP, derPP(k, 0:degPP), cos(a), sin(a) )
       !!!!!!!!!!!       min_f = max_f ! ONLY TEST

       if(min_f > 0.) then
          rho = max_f / min_f  ! initial aspect ratio
       else
          rho = 1E+5
       endif


       ! limitation of the maximal derivative
!       r_lim = 1.
!       if( state%modelName == 'NSe'  .and. max_f > max_u) then
          !u = max_f
          !max_f = (KK * max_u * u - max_u *max_u) / (u + max_u * (KK-2) )

          ! Var 1
          !!r_lim = max_f / u * (0.8)
          !!if( state%stop_adaptation == -1) r_lim = max_f / u * (0.8)**state%recompute_back

          ! Var 2
          !r_lim = max_f / u * (0.8)**state%recompute_back

!          r_lim = 0.


!          state%num_limits = state%num_limits + 1
          !write(*,'(a5,3i5,10es12.4)') 'WEW',elem%i, degP, k, u, max_f, rho, r_lim
!       endif

       !if(max_f > 1E+10) then
       !   write(*,'(a5,3i5,10es12.4)') 'WEW',elem%i, degP, k, max_f, max_a
       !   write(54,*) elem%xc(:), degP, k, max_f, 'WEW'
       !endif



       !write(*,'(a8,3i5,20es12.4)') '####@@!',elem%i, degP, k, max_f, min_f, rho

       goto 10 ! a simplification of the algorithm from APNUM 

!        ! setting of the upper estimate of type
!        ! DirectDeriv M= max_f ( x^T Q D Q^T x )**(degP / 2)
!        ! D = diag (1, rho**(-2/degP) )
       
!        allocate(D(1:2, 1:2), R(1:2, 1:2), RT(1:2, 1:2), M(1:2, 1:2) )
       
!        R(1,1) = cos(max_a)
!        R(1,2) = -sin(max_a)
!        R(2,1) = sin(max_a)
!        R(2,2) = cos(max_a)
       
!        RT(1:2, 1:2) = transpose(R(1:2, 1:2) )
       
!        D(2,1) = 0.
!        D(1,2) = 0.
       
!        f1 = max_f
!        f2 = rho
       
!        D(1,1) = 1.
!        D(2,2) = rho ** (-2./degP)
!        M(1:2, 1:2) = matmul(R(1:2, 1:2), matmul(D(1:2, 1:2), RT(1:2, 1:2)))
       
!        !if(elem%i == itest) then
!        !   ifig = ifig + 1
!        !   call DrawEstimate(ifig, elem%xc(1:2), max_f, degP, M(1:2, 1:2) )
!        !endif
       
!        !area = (rho / max_f**2)**(-1./degP)
!        !write(*,'(a6, 3i5,8es12.4)') &
!        !     '???',elem%i, 0, 0, max_f, min_f, rho, area
       
!        area_min = 1E+30
       
!        ! seeking of optimal values f_max and rho satisfying the estimate
!        ! for several values of f_max, we seek the minimal rho 
!        do j = 0, 10
!           max_f = f1 * (1.01**j)
!           !rho = f2
!           rho = 1E+6
          
!           D(2,2) = rho ** (-2./degP)
!           M(1:2, 1:2) = matmul(R(1:2, 1:2), matmul(D(1:2, 1:2), RT(1:2, 1:2)))
          
!           do i=1,n
!              !a = 2 * pi * i / n
!              a =  pi * i / n
!              xi(1) = cos(a)
!              xi(2) = sin(a)
             
!              ! directional derivative in direction a
!              f_der = DirectionalDerivative(degP, der(k, 0:degP), cos(a), sin(a) )
             
!              !f  = M(1,1) * cos(a)**2 + ( M(1,2) + M(2,1) )*cos(a) *sin(a) &
!              !     + M(2,2)*sin(a)**2
!              f = dot_product(xi(1:2), matmul(M(1:2, 1:2), xi(1:2) ) )
             
!              f_bound =  max_f * abs(f)**(degP/2.)
             
!              if(f_der > 1.001 * f_bound .and. a /=max_a) then
!                 ! estimate is violated, 
!                 ! we have to modify the anisotropy of the interpolation error function
                
!                 !! decrease of the ratio of the interpolation error function
!                 scale = ((f_der / max_f)**(2./degP) - (cos(a - max_a))**2) &
!                      /(sin(a - max_a))**2
!                 scale = scale**(degP/2.)
!                 min_f = max_f * scale 
                
!                 !! increase of the size of the interpolation error function
!                 !scale = f_der / f_bound
!                 !max_f = max_f * scale
!                 !min_f = min_f * scale
                
!                 ! refreshing of the matrix on the right-hand-side
!                 rho = max_f / min_f 
!                 D(2,2) = rho ** (-2./degP)
!                 M(1:2, 1:2) = matmul(R(1:2, 1:2), matmul(D(1:2, 1:2), RT(1:2, 1:2)))
                
!              endif
!           enddo  ! i=1,n
          
!           area = (rho / max_f**2)**(-1./degP)
!           !write(*,'(a6, 3i5,8es12.4)') &
!           !     '???',elem%i, j, i, max_f, min_f, rho, area
          
!           if(area < area_min) then
!              area_min = area
!              max_f_opt = max_f
!              rho_opt = rho
!           else
!              ! we expect that the dependence of area on max_f is strictly convex
!              ! so when are starts to increase, the optimal value was found
!              goto 20
!           endif

!           !if(elem%i == itest)  then
!           !   ifig = ifig + 1
!           !   call DrawEstimate(ifig, elem%xc(1:2), max_f, degP, M(1:2, 1:2) )
!           !endif
          
!        enddo   ! j=1,5
       
!        deallocate(D, R, RT, M)
       
! 20     continue
!        ! the optimal anisotropy of interpolation error function found
!        max_f = max_f_opt
!        rho = rho_opt
!        min_f = max_f / rho
!        area = (rho / max_f**2)**(-1./degP)
       
!        !write(*,'(a6, 3i5,8es12.4)') &
!        !     '???',elem%i, -99, -99, max_f, min_f, rho, area
!        !write(*,*) '(__________________________________)'
       
10     continue

       !write(*,'(a6, 3i5,8es12.4)') &
       !     '???',elem%i, -99, -99, max_f, max_a, min_f, rho, area

       !stop

       !if(elem%i == itest) stop
       Lq = state%space%adapt%Lq
       if(Lq <= -0.99) then  ! H^1-semi-norm
          ! balance || e||_K <= tol (|K|/|Omega|)^(1/q)
          area = ((tol_loc**2 * rho**0.5 *degP *pi**deg0 ) /  max_f)**(1./deg0)
          area = area /  state%space%domain_volume**(1./deg0 )
          ! !! area1 = area
          
          ! balance || e||_K <= tol (1 / N))^(1/q)
          !rN = grid%nelem**(1./Lq)
          !area = tol_loc * rho**0.5 / ( rN * c_p_q (degP, Lq) * max_f)
          !area = area**(2.*Lq /(Lq * degP + 2) )
          
          !write(*,'(a6,6es14.6)') '@@@@@', area, elem%area
       else
          stop 'TRoubles in datas in ama-hp_interpol.f90 EDREYTS'

       endif
       
       !write(*,'(a5,2i5,8es12.4)')  '!!!', elem%i, degP, &
       !     tol_loc, rho**0.5 , c_p_q (degP, Lq),  max_f, (2./degP)

       ! ????
       !beta_min = (rho**(1./degP) * area/pi)**0.5
       beta_min = (rho**(1./(2*deg0)) * area/1.)**0.5
       beta_max = beta_min / rho**(1./(2*deg0))
       
       ! new tests of limitation
       r_lim = 1.
       if( state%modelName == 'NSe' ) then
          !h_min = 5E-4
          !h_min = 1E-3
          h_min = state%model%Re1 * degP
          !h_min = 5E-4 * degP  ! should be the same as ^^^ (state%model%Re1 * degP )

          if(state%type_IC == 8) h_min = 0.005
          !!!if(beta_max < 5E-3) print*,'####', beta_max , h_min

          if( beta_max < h_min) then
             !write(*,'(a4,4i5,8es12.4)') '###',elem%i,degP,k,  state%num_limits, beta_min, beta_max,&
             !     beta_min/ beta_max, max_f
             beta_max = h_min
             beta_min =   beta_max * rho**(1./degP)
             r_lim = 0.
             state%num_limits = state%num_limits + 1
          
          endif
       endif

       !beta_min = max(beta_min, 1E-3)
       !beta_max = max(beta_max, 1E-3)

       !if(elem%i == itest) then
       !    print*,'************************************'
       !    print*,'***', beta_max, beta_min, max_f
       !    print*,'***', (tol_loc / max_f )**(1./ degP), &
       !         (tol_loc / min_f )**(1./ degP)  
       !write(192,'(a5,2i5,12es10.2)') '???', elem%i, degP, &
       !     beta_max, beta_min, max_f, min_f, max_f / min_f,  beta_max/  beta_min
       !        lam_max, lam_min, max_a, a
       !endif
       
       !beta_min = ( tol_loc / min_f )**(1./ degP) 

       !write(*,'(a6,6es14.6,a3,4es12.4)') &
       !     '@@@@@', area, elem%area, beta_max, beta_min, elem%diam, elem%area/elem%diam,'|', &
       !     elem%errH1
       
       lam_max = 1./beta_max**2
       lam_min = 1./beta_min**2

       Kpar_loc(k,1) = lam_max
       Kpar_loc(k,2) = lam_min
       Kpar_loc(k,3) = max_a

       if(elem%i == -1) &
          write(*,'(a6,3i5,3es14.6, a2,20es12.4)') 'WW',elem%i,degP,k,Kpar_loc(k,1:3),'|' !,  der(k, 0:degP)
    enddo  ! end of k


    ! metric compositions
    if(ndimL > 1) then

       ! SMAZ the following
       !Kparams(1) = lam_max
       !Kparams(2) = lam_min
       !Kparams(3) = max_a
       Kparams(4) = max_f
       
       
       call MetricComposition(ndimL, Kpar_loc(1:ndimL, 1:3),  Kparams(1:3), elem%xc(1:2), elem%i, degP )

    else
       Kparams(1) = lam_max
       Kparams(2) = lam_min
       Kparams(3) = max_a
       Kparams(4) = max_f

    endif
 
     if(elem%i == -1) then
        write(*,'(a6,3i5,3es14.6, a2,20es12.4)') 'QQ',elem%i,degP,0,Kparams(1:3),'|' !,  der(k, 0:degP)
        print*,'#########################################################'
     endif
     !!if(degP == elem%deg+ 2) stop

     deallocate(Kpar_loc)

300  continue

    ! estimates of the interpolation error in the L^{\infty} norm
    if(degP == elem%deg + 1) then
       
       !do k=1, ndim
       do k=1, 1 ! from the first component only (the actual solution)
          elem%interL8 = 0.

          if(Lq >= 1.) then 
             Lqq = Lq
          else
             Lqq = 2.    
          endif

          Qdof = elem%Qdof
          
          allocate( Fx(1:Qdof, 1:nbDim), wi(1:Qdof) )
          !integration nodes on K
          call ComputeF(elem, Qdof, state%space%V_rule(elem%Qnum)%lambda(1:Qdof,1:nbDim), &
               Fx(1:Qdof, 1:nbDim) )
          
          !L^q -norm
          ! interpolation error function in integ. nodes
          do i=1,Qdof
             xi(1:2) =Fx(i, 1:2) - elem%xc(1:2)
             wi(i) = DirectionalDerivative(degP, der(k, 0:degP), xi(1),xi(2) )

             ! limitation of the maximal derivative
             wi(i) = wi(i) * r_lim

             ! q-th power
             wi(i) = (abs(wi(i)) )**Lqq
          enddo
          ! error  in the Lq norm
          call IntegrateFunction(elem, wi, elem%interLq)
          
          ! NOT POWERED !!!
          !elem%interLq = elem%interLq**(1./Lqq)
          
          ! interpolation error function estimate in L^\infty norm
          do i=1,3  ! vertexes of triangle
             xi(1:2) = grid%x(elem%face(idx, i), 1:2) - elem%xc(1:2)

             f = DirectionalDerivative(degP, der(k, 0:degP), xi(1),xi(2) )

             ! limitation of the maximal derivative
             f = f * r_lim

             elem%interL8 = max ( elem%interL8, f )
             !write(*,'(a6,3i5,6es12.4)') &
             !     '??',elem%i, i, elem%face(idx, i), xi(1:2), abs(f),elem%interL8
          enddo ! i 

          ! size of the non-counted area
          !if(r_lim == 0.)  state%err(algeb) = state%err(algeb) + elem%area

          ! H^1-seminorm
           do i=1,Qdof
             xi(1:2) =Fx(i, 1:2) - elem%xc(1:2)
             wi(i) = DirectionalDerivative_H1(degPP, derPP(k, 0:degPP), xi(1),xi(2) )

             ! limitation of the maximal derivative
             wi(i) = wi(i) * r_lim

             ! q-th power
             wi(i) = (abs(wi(i)) )
          enddo
          ! error  in the Lq norm
          call IntegrateFunction(elem, wi, elem%interH1)
          ! NOT POWERED !!!
          !elem%interH1 = elem%interH1**(0.5)

          !!write(*,'(a6,2i5, 4es14.6)') 'EDE',elem%i, degPP,elem%interLq, elem%interL8, elem%interH1
          deallocate(Fx, wi)
       

       enddo  ! k
    endif

    !deallocate(du_val)
    deallocate(derPP)

  end subroutine FindAnisotropyEIp_H1


  !> composition of several metrices by an intersection of ellipses
  subroutine MetricComposition(ndimL, parL,  parG, xc, ie, degP)
    integer, intent(in):: ndimL
    real, dimension(1:ndimL, 1:3), intent(in) :: parL
    real, dimension(1:3), intent(inout) :: parG
    real, dimension(1:2), intent(in) :: xc  ! elemen barycentre
    integer, intent(in):: ie, degP  ! element index, pol. degree
    real, dimension(:,:), allocatable :: rgabc
    real, dimension(:,:), allocatable :: A
    real, dimension(:,:), allocatable :: M,Q, ZZ,ZZ1
    real, dimension(:), allocatable :: D,E,w, work
    external :: DSYEV,  DPTEQR
    integer :: lwork, info
    integer :: k, ikf, itest

    itest = -1  !-146

    ikf = 100 + degP*10

    if(ie == itest) write(ikf-1,*) xc(:)

    allocate( rgabc(1:ndimL, 1:3) )

    k = 1
    rgabc(k,1) =  parL(k,1) * cos(parL(k,3))**2 + parL(k,2) * sin(parL(k,3))**2
    rgabc(k,2) = (parL(k,1) - parL(k,2)) * cos(parL(k,3)) * sin(parL(k,3))
    rgabc(k,3) =  parL(k,1) * sin(parL(k,3))**2 + parL(k,2) * cos(parL(k,3))**2

    if(ie == itest)  call DrawEllips(ikf+1, rgabc(k, 1:3), xc(1:2))

    if(ie == itest) write(*,'(a6, 6es12.4)' ) 'init', rgabc(1,1:3), xc(1:2)

    do k=2, ndimL
       rgabc(k,1) =  parL(k,1) * cos(parL(k,3))**2 + parL(k,2) * sin(parL(k,3))**2
       rgabc(k,2) = (parL(k,1) - parL(k,2)) * cos(parL(k,3)) * sin(parL(k,3))
       rgabc(k,3) =  parL(k,1) * sin(parL(k,3))**2 + parL(k,2) * cos(parL(k,3))**2

       if(ie == itest)  call DrawEllips(ikf+k, rgabc(k, 1:3), xc(1:2))

       call INTERSECTION_METRIC(rgabc(1,1), rgabc(1,2), rgabc(1,3), rgabc(k,1), rgabc(k,2), rgabc(k,3))
       if(ie == itest) write(*,'(a6, 6es12.4)' ) 'init', rgabc(1,1:3), rgabc(k,1:3)

    enddo

    if(ie == itest)  call DrawEllips(ikf, rgabc(1, 1:3), xc(1:2))


    ! setting of eigenvalues and eigenvectors
    ! lwork = 68
    ! allocate(A(1:2, 1:2), w(1:2), work(1:lwork) )
    ! A(1,1) = rgabc(1, 1)
    ! A(1,2) = rgabc(1, 2)
    ! A(2,1) = 0. !rgabc(1, 2)
    ! A(2,2) = rgabc(1, 3)
    ! call DSYEV( 'V', 'U', 2, A, 2, W, WORK, LWORK, INFO )
    ! if(info /= 0) then
    !    print*,'Troubles in DSYEV in ama-hp_interpol.f90', INFO
    !    stop
    ! endif

    ! parG(1) = w(2)
    ! parG(2) = w(1)
    ! parG(3) = acos(A(1,1))
    ! if(A(1,1) < 0.) parG(3) =pi - acos(A(1,1))

    ! if(ie == itest) &
    !      write(*,'(a6,3i5,3es12.4, a2,20es12.4)') 'WW',ie,degP,k,parG(1:3),'%'

    ! if(ie == itest) then
    !    print*
    !    write(*,'(a8,6es12.4)') '@@@', w(1:2), work(1:2)
    !    write(*,'(a8,6es12.4)') '@@@', A(1,:)
    !    write(*,'(a8,6es12.4)') '@@@', A(2,:)
    !    print*,'#########################################################'
    ! endif
    ! deallocate(A, w,  work)

    ! second possibility

    ! setting of eigenvalues and eigenvectors
    ! allocate(M(1:2,1:2), Q(1:2, 1:2) , ZZ(1:2, 1:2))
    ! M = 0.
    ! M(1,1) = parL(1,1)
    ! M(2,2) = parL(1,2)
    ! Q(1,1) = cos(parL(1,3))
    ! Q(1,2) =-sin(parL(1,3))
    ! Q(2,1) = sin(parL(1,3))
    ! Q(2,2) = cos(parL(1,3))
    ! ZZ = matmul(Q, matmul(M, transpose(Q) ))
    ! write(*,'(a8,6es12.4)') '@ZZ', ZZ(1,:)
    ! write(*,'(a8,6es12.4)') '@ZZ', ZZ(2,:)


    lwork = 8
    allocate(D(1:2),E(1:1), w(1:2), A(1:2, 1:2), work(1:lwork) )
    D(1) = rgabc(1, 1)
    D(2) = rgabc(1, 3)
    E(1) = rgabc(1, 2)
    call DPTEQR( 'I', 2, D, E, A, 2, WORK, INFO )
    !call DSYEV( 'V', 'U', 2, A, 2, W, WORK, LWORK, INFO )
    if(info /= 0) then
       print*,'Troubles in DSYEV in ama-hp_interpol.f90', INFO
       stop
    endif

    parG(1) = D(1)
    parG(2) = D(2)
    if(A(2,1) >= 0) then
       parG(3) = acos(A(1,1))
    else
       parG(3) = acos(-A(1,1))
    endif


    ! M = 0.
    ! M(1,1) = D(1)
    ! M(2,2) = D(2)

    if(ie == itest) &
         write(*,'(a6,3i5,3es14.6, a2,20es12.4)') '*W',ie,degP,k,parG(1:3),'%'

    if(ie == itest) then
       print*
       write(*,'(a8,6es12.4)') '*DD', D(1:2)
       write(*,'(a8,6es12.4)') '@EV', A(1,:)
       write(*,'(a8,6es12.4)') '@EV', A(2,:)
       print*

    ! ZZ1 = matmul(A , matmul(M, transpose(A)  ))
    ! write(*,'(a8,6es12.4)') '@ZZ1', ZZ1(1,:)
    ! write(*,'(a8,6es12.4)') '@ZZ1', ZZ1(2,:)

    ! print*, parG(1) * cos(parG(3))**2 + parG(2) * sin(parG(3))**2
    ! print*, (parG(1) - parG(2)) * cos(parG(3)) * sin(parG(3))
    ! print*, parG(1) * sin(parG(3))**2 + parG(2) * cos(parG(3))**2

    endif
    deallocate(A, w,  work, D, E)

    deallocate(rgabc )
  end subroutine MetricComposition


  function c_p_q(p, q)
    real :: c_p_q
    real, intent(in) :: q
    integer, intent(in) :: p  ! p is already p+1
    real :: pi
    pi = 2 * asin(1.)

    c_p_q = (1./pi)**(q*p/2 + 1)
    c_p_q = c_p_q *2 * pi /(q*pi +2 )
    c_p_q = c_p_q ** (1./q)

  end function c_p_q

  !> compute the directional derivative in the direction ( cos(a), sin(a) )
  function DirectionalDerivative(degP, der, a1, a2)
    real :: DirectionalDerivative
    integer, intent(in) :: degP
    real, dimension(0:degP), intent(in) :: der  ! partial derivatives of degree degP
    real, intent(in) :: a1, a2
    integer :: j
    real :: f, fac

    f = 0.
    do j=0, degP
       fac = factorial(j) * factorial(degP - j)
       f = f + a2**j * a1**(degP - j) * der(j)  / fac
    enddo

    DirectionalDerivative = abs(f)

  end function DirectionalDerivative

  !> compute the directional derivative in the direction ( cos(a), sin(a) ), without facorial factors
  function DirectionalDerivative_H1(degP, der, a1, a2)
    real :: DirectionalDerivative_H1
    integer, intent(in) :: degP
    real, dimension(0:degP), intent(in) :: der  ! partial derivatives of degree degP
    real, intent(in) :: a1, a2
    integer :: j
    real :: f, fac

    f = 0.
    do j=0, degP
       f = f + a2**j * a1**(degP - j) * der(j)
    enddo
    
    DirectionalDerivative_H1 = abs(f)

  end function DirectionalDerivative_H1



  !> draw the right and side of the estimate
  subroutine DrawEstimate(ifig, xc, max_f, degP, M )
    integer, intent(in) :: ifig, degP
    real, intent(in) :: max_f
    real, dimension(1:2), intent(in) :: xc
    real, dimension(1:2, 1:2), intent(in) :: M
    real :: a, pi, f, xi(1:2 ), tol
    integer :: i, n
    pi = 3.14159267

    n = 100

    tol = 0.1

    do i=1,n

       a = 2 * pi * i / n
       xi(1) = cos(a)
       xi(2) = sin(a)

       f = dot_product(xi(1:2), matmul(M(1:2, 1:2), xi(1:2) ) )
       f = max_f * f **( degP/2. )
       write(ifig, *) xc(1) + cos(a) * f, xc(2) + sin(a) * f, &
            xc(1) + cos(a) * (tol/f)**(1./degP), xc(2) + sin(a) * (tol/f)**(1./degP)
    end do
  end subroutine DrawEstimate


    subroutine INTERSECTION_METRIC(ra1,rb1,rc1,ra2,rb2,rc2)
      real, intent(inout) :: ra1,rb1,rc1,ra2,rb2,rc2
      real:: a1,a2,b1,b2,c1,c2,p1,p2,p3,p4,d,rmat1(2),rmat2(2)
      real:: rlam1,rlam2,pi1,pi2,pi3,pi4
      real:: rmat112, rmat121, rmat212, rmat221

      a1 = ra1
      b1 = rb1
      c1 = rc1
      a2 = ra2
      b2 = rb2
      c2 = rc2

      p3 = 1.D+00
      p4 = 1.D+00

      d = (a1*c2 - a2*c1)**2 -4.D+00*(a1*b2-b1*a2)*(b1*c2-c1*b2)
      if((a1*b2-b1*a2) .ne. 0 .and. d .ge. 0.) then
         p2 = (a2*c1-a1*c2 - d**0.5)/2.D+00/(a1*b2-b1*a2)
      elseif(a2*c1-a1*c2 .ne. 0. .and. d .ge. 0. ) then
         p2 = (b1*c2-c1*b2)/(a2*c1-a1*c2)
      else
!         print *,'zero polynom'
         if(a1 .le. a2) then
            ra1 = ra2
            rb1 = rb2
            rc1 = rc2
         endif
         return
      endif

      if(p2*a2 + b2 .ne. 0.) then
         p1 = -(b2*p2 + c2)/(p2*a2 + b2)
      else
         if(a1 .le. a2) then
            ra1 = ra2
            rb1 = rb2
            rc1 = rc2
         endif
         return
      endif

      rmat1(1) = p1**2*a1+2*p1*p3*b1+p3**2*c1
      rmat112 = p2*p1*a1+p2*p3*b1+p4*p1*b1+p4*p3*c1
      rmat121 = p2*p1*a1+p2*p3*b1+p4*p1*b1+p4*p3*c1
      rmat1(2) = p2**2*a1+2*p2*p4*b1+p4**2*c1

      rmat2(1) = p1**2*a2+2*p1*p3*b2+p3**2*c2
      rmat212 = p2*p1*a2+p2*p3*b2+p4*p1*b2+p4*p3*c2
      rmat221 = p2*p1*a2+p2*p3*b2+p4*p1*b2+p4*p3*c2
      rmat2(2) = p2**2*a2+2*p2*p4*b2+p4**2*c2

      rlam1 = max(rmat1(1),rmat2(1) )
      rlam2 = max(rmat1(2),rmat2(2) )

      if(rlam1 .le. 0. .or. rlam2 .le. 0) then
         print *,'nonpositive eigenvalues'
         stop
      endif

      d = p1*p4 - p2*p3
      if(d == 0.) then
         print *,'zero deter.'
         stop
      endif
      pi1 = p4/d
      pi2 = -p2/d
      pi3 = -p3/d
      pi4 = p1/d

      ra1 = pi1**2*rlam1+pi3**2*rlam2
      rb1 = pi1*rlam1*pi2+pi3*rlam2*pi4
      rc1 = pi2**2*rlam1+pi4**2*rlam2

      return
    end subroutine INTERSECTION_METRIC


  !> seek the direction having the maximal directional derivative of order degP
  !> direction is max_a, size of the derivative lam_max, size of the derivative in
  !> the perpendicular direction is lam_min
  !> paper AM2013
  subroutine FindAnisotropyEIpOLD(elem, ndimL, degP, der, tol_loc_max, lam_max, lam_min,  max_a, max_f)
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    integer, intent(in) :: degP  ! degree of the directional derivative
    real, dimension(1:ndimL,  0:degP), intent(in) :: der  ! partial derivatives of degree degP
    real, intent(in) :: tol_loc_max   ! local tolerance for the element
    real, intent(out) :: lam_max, lam_min,  max_a ! characterization of the metric
    real, intent(out) :: max_f
    real, dimension(1:2) :: xi
    real :: beta, beta_max, beta_min, min_f
    real :: pi, fac, f, a
    integer :: n, k, j, i, ip, itest

    itest = 25

    ! seeking the element size, orientation and the aspect ratio
    !n = 100 ! number of test directions
    n = 180 ! number of test directions
    pi = 2 * asin(1.)

    fac = factorial(degP)

    !do k=1, ndim
    do k=1, 1 ! from the density only
       max_f = 0.
       do i=1,n
          a = 2 * pi * i / n

          !f = 0.
          !do j=0, degP
          !   fac = factorial(j) * factorial(degP - j)
          !   f = f + sin(a)**j * cos(a)**(degP - j) * der(k, j)  / fac
          !enddo
          !f = abs(f)

          f = DirectionalDerivative(degP, der(k, 0:degP), cos(a), sin(a))

          !beta = ( tol_loc_max / f )**(1./ degP)

          if( f  > max_f )then
             max_f = f
             max_a = a
             !beta_max = beta
          endif

          !xi(1:2) = (/ 0.5, 0.95/)
          !if(VectorNorm(elem%xc(1:2) - xi(1:2) ) <= 0.030)  then
          !   ip = 10*(state%space%adapt%adapt_level+1) + 6
          !  if(i == 1) then
          !      write(ip+0,'(x)' )
          !      write(ip+1,'(x)' )
          !      write(ip+2,'(x)' )
          !   endif
          !   write(ip+0,*) a, f
          !   write(ip+1,*) elem%xc(1) + cos(a) * f   , elem%xc(2) + sin(a) * f
          !   write(ip+2,*) elem%xc(1) + cos(a) * beta, elem%xc(2) + sin(a) * beta
          !endif

       enddo


       ! size in the parallel direction
       a = max_a + pi /2
       !f = 0.
       !do j=0, degP
       !   fac = factorial(j) * factorial(degP - j)
       !   !f = f + sin(a)**j * cos(a)**(degP - j) * der(k, j) !(/fac) !!! CHYBA ???
       !f = f + sin(a)**j * cos(a)**(degP - j) * der(k, j) / fac !!! OK
       !enddo
       !!min_f = abs(f) / fac   !!! CHYBA ???
       !min_f = abs(f)   !!! OK

       min_f = DirectionalDerivative(degP, der(k, 0:degP), cos(a), sin(a) )


       beta_max = ( tol_loc_max / max_f )**(1./ degP)
       beta_min = ( tol_loc_max / min_f )**(1./ degP)


       lam_max = 1./beta_max**2
       lam_min = 1./beta_min**2

       !if(elem%i == itest) then
       !    print*,'************************************'
       !    print*,'***', beta_max, beta_min, max_f
       !    print*,'***', (tol_loc_max / max_f )**(1./ degP), &
       !         (tol_loc_max / min_f )**(1./ degP)
       !    write(191,'(a5,2i5,12es10.2)') '???', elem%i, degP, &
       !         beta_max, beta_min, max_f, min_f, max_f / min_f,  beta_max/  beta_min
       !        lam_max, lam_min, max_a, a
       !endif

       !write(*,'(a3,i2,a22,6es10.2)') 'p = ', p, ', max_f, min_f, angle:',max_f, min_f, max_a, a, &
       !     lam_max, lam_min


       !write(88,*) elem%xc(:), elem%wS(k, 0:degP)

    end do

    ! estimates of the interpolation error in the L^{\infty} norm
    if(degP == elem%deg + 1) then

       !do k=1, ndim
       do k=1, 1 ! from the density only
          elem%interL8 = 0.

          do i=1,3  ! vertexes of triangle
             xi(1:2) = grid%x(elem%face(idx, i), 1:2) - elem%xc(1:2)
             f = 0.
             do j=0, degP
                fac = factorial(j) * factorial(degP - j)
                f = f +xi(2)**j * xi(1)**(degP - j) * der(k, j)  / fac
             enddo ! j

             elem%interL8 = max ( elem%interL8, abs(f) )
             !write(*,'(a6,3i5,6es12.4)') &
             !     '??',elem%i, i, elem%face(idx, i), xi(1:2), abs(f),elem%interL8
          enddo ! i
       enddo  ! k
    endif

  end subroutine FindAnisotropyEIpOLD


  !> evaluate all partial derivatives of \f$ w \f$ on elem
  subroutine Eval_All_Derivatives(elem, ndimL)   !!!, wp)
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    !real, dimension(1:grid%npoin, 1:ndim), intent(inout) :: wp
    real, dimension(:,:,:,:), allocatable :: Dw    ! all derivatives of w: bases coefficients
    real, dimension(:,:,:,:), allocatable :: Dwi   ! all derivatives of w in integ node
    real, dimension(:,:,:), allocatable :: Dphi ! derivative of the test functions
    real, dimension(:,:), allocatable :: MassInv ! Inversion of local mass matrix of order deg+1!!
    !real, dimension(:,:), allocatable :: Mass, MassS ! local mass matrix of order deg+1!!
    real, dimension(:), allocatable :: ident, vec ! temporary array
    integer :: Qnum, Qdof, dof, dofP, deg, degP
    integer :: ideg, i, j, k, l,  ix, jx

    deg = elem%deg
    dof = elem%dof
    !Qnum = elem%Qnum + Qdiff

    ! increase of quadrature
    if(deg == MaxDegreeImplemented ) then
        print*,'$$ alert, possible troubles in ama-hp_metric.f90'
        Qnum = elem%Qnum
     else
        Qnum = state%space%Qdeg(deg+1, 1)
     endif

    !!if(Qnum > 20 ) print*,'#####',elem%deg, elem%Qnum, Qnum, maxVrule

    !if(Qnum > maxVrule ) then
    !   Qnum = maxVrule
    !   Qdiff = Qnum - elem%Qnum
    !   print*,'$$ alert, possible troubles in ama-hp_metric.f90'
    !endif

    Qdof = state%space%V_rule(Qnum)%Qdof

    !!!!!Qdof = elem%Qdof OLD!!

    do ideg = 0,2  ! we test the interpolation errors of order deg, deg+1, deg+2

       degP = deg + ideg

       if(deg <= MaxDegreeImplemented) then

          dofP = DOFtriang( degP )

          !local mass matrix of order dofP x dofP
          allocate(MassInv(1:dofP, 1:dofP) )
          MassInv(:,:) = 0.
          !allocate(Mass(1:dofP, 1:dofP), MassS(1:dofP, 1:dofP) )
          !Mass(:,:) = 0.

          ! identity vector
          allocate(ident(1:Qdof) )
          ident(:) = 1.

          ! evaluation of the inverse of the mass matrix
          !call IntegrateBlockBBplus(elem, Qnum, dofP, ident(1:Qdof), Mass(1:dofP,1:dofP) )

          call IntegrateBlockBBplus(elem, Qnum, dofP, ident(1:Qdof), MassInv(1:dofP,1:dofP) )
          call MblockInverse(dofP, MassInv(1:dofP, 1:dofP) )


          ! Dw(k,i,j, l) => k-th component of w, (i,j) d^i/(dx^j dy^(i-j)), l coeff or integ node
          allocate(Dw (1:ndimL, 0:degP, 0:degP, 1:dofP) )
          allocate(Dwi(1:ndimL, 0:degP, 0:degP, 1:Qdof) )
          Dw = 0.
          Dwi = 0.

          ! higher order reconstruction via Oswald interpolation
          !call OswaldInterpolation(elem, ndimL,degP, dofP,  wp(1:grid%npoin, 1:ndim), &
          !     MassInv(1:dofP, 1:dofP), Dw (1:ndimL, 0:degP, 0:degP, 1:dofP) )

          ! higher order reconstruction via least square method
          call LeastSquareInterpolationH1(elem, ndimL, .false., Qnum, degP, dofP, &
               Dw (1:ndimL, 0:degP, 0:degP, 1:dofP), 0 )

          ! do i=0,degP
          !    do j=0,degP
          !       write(*,'(a3,2i5,200es14.6)') '@@@',i,j,Dw(1,i,j,:)
          !    enddo
          ! enddo
          !stop

          !call LeastSquareInterpolationL2(elem, ndimL, .false., Qnum, degP, dofP, &
          !     Dw (1:ndimL, 0:degP, 0:degP, 1:dofP), 0 )

          allocate(vec(1:dofP) )

          allocate(Dphi(1:dofP, 1:nbDim, 1:Qdof ) )
          call  Eval_Dphi_plus(elem, state%space%V_rule(Qnum), dofP,  Dphi(1:dofP, 1:nbDim, 1:Qdof) )

          do k=1,ndimL

             do i=1,degP   ! i = i-th derivative
                do j=0,i  !   d^i w /dx^j dy^(i-j)

                   ! evaluation of the partial derivative in integ nodes
                   do l=1, Qdof  ! integ node

                      ! jx which (i-1)th derivative has to be integrate
                      ! ix according which variable: ix = 1 according x_1, ix = 2 according x_2
                      if(j==0) then
                         jx = 0
                         ix = 1
                      elseif(j==1 .and. degP >2) then
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

          !write(99,*) '-----', elem%i, '  ------------------------------'
          !k = 1
          !do i=1,degP   ! i = i-th derivative
          !   do j=0,i  !   d^i w /dx^j dy^(i-j)
          !      write(99,'(a6,3i5,120es12.4)') 'Dwi:',i,j,Qdof,Dwi(k, i, j, 1:Qdof)
          !   enddo
          !enddo

          !write(99,*) '-----------------------------------'

          !do i=0,degP   ! i = i-th derivative
          !   do j=0,i  !   d^i w /dx^j dy^(i-j)
          !      write(99,'(a6,3i5,120es12.4)') 'Dw:',i,j,dofP,Dw(k, i, j, 1:dofP)
          !   enddo
          !enddo

          !write(99,*) '########################################'

          do k=1, ndimL
             do i=0, degP
                ! elem%wSS contains only the derivatives of maximal order
                !Dwi(k,degP, i, 1:Qdof) is almost constant
                elem%wSS(k, ideg, i) = sum(Dwi(k,degP, i, 1:Qdof) )/Qdof
             enddo
             !write(100 + 10*state%space%adapt%adapt_level + ideg,'(120es12.4)') &
             !write(*,'(120es12.4)') &
             !     elem%xc(:), elem%wSS(k, ideg, 0: degP), &
             !     1.*ideg, 1.*degP, 1.* elem%i
          enddo
       else
          ! this degree is not implemented, hence too high derivative
          elem%wSS(1:ndimL, ideg, 0:degP) = 1.E+60
       endif

       deallocate(Dw, Dwi,  vec, Dphi, ident)
       deallocate( MassInv)
       !deallocate( Mass, MassS)

    end do  ! ideg
    !!stop


  end subroutine Eval_All_Derivatives

  !> evaluate the Oswald interpolation on elem of order /f$ P^{degP} /f$
  ! UNCORRECTED ndim --> ndimL
  subroutine OswaldInterpolation(elem, degP, dof1, wp, MassInv, Dw  )
    type(element), intent(inout) :: elem
    !!integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    integer, intent(in) :: degP, dof1
    real, dimension(1:grid%npoin, 1:ndim), intent(inout) :: wp
    real, dimension(1:ndim, 0:degP, 0:degP, 1:dof1), intent(inout) :: Dw
    real, dimension(1:dof1, 1:dof1), intent(in) :: MassInv
    class(element), pointer :: elem1
    type(Lagrang_rule), pointer :: L_rule
    real, dimension(:,:), pointer :: phi ! the test functions
    real, dimension(:,:), allocatable :: Fx ! real physical coordinates
    real, dimension(:,:), allocatable :: wi ! solution in Lagrangian integ nodes
    real, dimension(:), allocatable :: wii   ! solution in Lagrangian integ nodes of neighbour
    real, dimension(:), allocatable :: f, qi  ! temporary array
    real, dimension(:,:), allocatable :: Lpsi    ! Lagr. test functions in V_rule integ nodes
    integer :: Ldof, Lnum, Qdof, dof, dofii
    integer :: i, i1, j, k, l, l1, l2

    ! evaluation of a "high-order projection" of w
    dof = elem%dof
    Qdof = elem%Qdof

    Lnum = degP
    L_rule => state%space%L_rule(Lnum)
    Ldof = L_rule%Qdof
    phi => L_rule%phi(1:dof1, 1:Ldof)

    allocate(Fx(1:Ldof, 1:2) )
    call ComputeF(elem, Ldof, L_rule%lambda(1:Ldof, 1:2), Fx(1:Ldof, 1:2) )

    !do i=1, dof1
    !   write(*,'(a6,3i5,30es12.4)') 'L_rule:',Lnum, i,Ldof, phi(i, 1:Ldof)
    !enddo

    allocate( wi(1:ndim,1:Ldof), wii(1:ndim) )
    do k=1, ndim
       !!WS!! wi(k, 1:Ldof) = matmul(elem%w(0, (k-1)*dof +1 : k*dof), phi(1:dof, 1:Ldof) )
       wi(k, 1:Ldof) = matmul(elem%wS(k, 1 : dof), phi(1:dof, 1:Ldof) )
    enddo

    !write(*,*) '---------------------------'
    !do k=1,ndim
    !!   write(*,'(a6,i5,30es12.4)') 'wi:',k, wi(k, 1:Ldof)
    !   do l=1, Ldof
    !      write(98,*) L_rule%lambda(l, 1:2), wi(k, l), Fx(l, 1:2)
    !   enddo
    !enddo

    ! averaging at vertices => we replace wi(k, vertex) by wp(vertex , k)
    do i=1, elem%flen
       k = elem%face(idx, i)         ! global index of vertex
       j = VertexIndex(Lnum, i)      ! local index of vertex as Lagrangian node
       wi(1:ndim, j) = wp(k, 1:ndim)

       !k = 1
       !write(*,'(a6,i5,30es12.4)') 'wiV:',k, wi(k, 1:Ldof)
    enddo

    !do l=1, Ldof
    !   write(97,*) L_rule%lambda(l, 1:2), wi(1, l), Fx(l, 1:2)
    !enddo
    !write(*,*) '---------------------------'

    ! averaging at edges => average from the neighbouring element
    do i=1, elem%flen
       k = elem%face(neigh, i)         ! global index of vertex
       if(k > 0) then ! inner face
          elem1 => grid%elem(k)
          dofii = elem1%dof
          i1 = elem%face(nei_i, i)

          do l = 1, Lnum-1             ! only interior Lagr. nodes on edges
             l1 = EdgeIndex(Lnum,  i, l)       ! local index of vertex as Lagrangian node
             l2 = EdgeIndex(Lnum, i1, Lnum -l) ! from opposite side

             do k=1, ndim
                !!WS!! wii(k) = dot_product(elem1%w(0, (k-1)*dofii +1 : k*dofii), phi(1:dofii, l2) )
                wii(k) = dot_product(elem1%wS(k, 1 : dofii), phi(1:dofii, l2) )
             enddo

             !write(*,'(a6,3i5,9es12.4)') &
             !     '####',i,l,l1, L_rule%lambda(l1, 1:2), 1 - sum(L_rule%lambda(l1, 1:2)), wi(1,l1)
             !write(*,'(a6,3i5,9es12.4)') &
             !     '####',i,l,l2, L_rule%lambda(l2, 1:2), 1 - sum(L_rule%lambda(l2, 1:2)), wii(1)
             !print*,'---'

             wi(1:ndim, l1) = ( wi(1:ndim, l1) + wii(1:ndim) ) / 2

             !write(*,'(a6,i5,30es12.4)') 'wiE:',1, wi(1, 1:Ldof)
          enddo
       end if
    enddo


    !do l=1, Ldof
    !   write(96,*) L_rule%lambda(l, 1:2), wi(1, l), Fx(l, 1:2)
    !   !! wi(1,l) = Fx(l, 1)**2 + 0.05 * Fx(l, 2)**2   ! test case
    !enddo
    !write(*,*) '---------------------------'

    ! wi(1:ndim, 1:Ldof) contains the values in (elem%deg +1)-Lagr. nodes,
    ! we recompute them as a polynomial fucntion
    allocate( Lpsi(1:Ldof, 1:Qdof)  ) ! Lagrang. test functions in V_rule integ nodes
    allocate( f(1:Qdof), qi(1:dof1) )

    call Eval_L_Direct(L_rule, Qdof, state%space%V_rule(elem%Qnum)%lambda(1:Qdof,1:3), Lpsi(1:Ldof, 1:Qdof) )
    !call Eval_L_Direct(L_rule, Ldof, state%space%L_rule(Lnum)%lambda(1:Ldof,1:2), Lpsi(1:Ldof, 1:Ldof) )

    !do k=1,Qdof
    !   !write(*,'(a6,i5, 30es12.4)') 'Lpsi:',k, Lpsi(k, 1:Qdof)
    !   write(95, *) state%space%V_rule(elem%Qnum)%lambda(k,1:2),Lpsi(1:Ldof, k)
    !enddo
    !stop

    deallocate(Fx)
    allocate(Fx(1:Qdof, 1:2) )
    call ComputeF(elem, Qdof, state%space%V_rule(elem%Qnum)%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )

    do k=1, ndim
       f(1:Qdof) = matmul(wi(k,1:Ldof), Lpsi(1:Ldof, 1:Qdof) )

       qi(:) = 0.
       call IntegrateVectorB(elem, dof1, f(1:Qdof), qi(1:dof1) )

       ! zero order derivative = a reconstruction of  function w
       Dw(k, 0, 0, 1:dof1)  = matmul(MassInv(1:dof1, 1:dof1), qi(1:dof1) )
    enddo

    !do l=1, Qdof
    !   write(95,*) state%space%V_rule(elem%Qnum)%lambda(l, 1:2), f(l), Fx(l, 1:2)
    !enddo

    deallocate( Fx, wi, wii, f, qi )

  end subroutine OswaldInterpolation

  !> evaluate the solution in vertexes by Oswald interpolation
  ! UNCORRECTED ndim --> ndimL
  subroutine EvalSolutionVertexes(wp)
    real, dimension(1:grid%npoin, 1:ndim), intent(inout) :: wp
    class(element), pointer :: elem
    type(Lagrang_rule), pointer :: L_rule
    real, dimension(:,:), pointer :: phi   ! Lagrangian interpolation
    real, dimension(:,:), allocatable :: wi
    real, dimension(:), allocatable :: supp
    integer :: dof, Qdof, QLdof
    integer :: i, j, k, l, ix, jx

    allocate(supp(1:grid%npoin))
    supp(:) = 0.
    wp(:,:) = 0.

    L_rule => state%space%L_rule(1)  ! 1 => P_1 Lagrangian functions
    QLdof = L_rule%Qdof
    phi => L_rule%phi(1:state%space%max_dof, 1:QLdof)

    !do i=1, state%space%max_dof
    !   write(*,'(a6,2i5,30es12.4)') 'L_rule:', i,QLdof, phi(i, 1:QLdof)
    !enddo

    allocate( wi(1:ndim,1:QLdof) )

    do i=1,grid%nelem
       elem => grid%elem(i)
       dof = elem%dof

       do k=1, ndim
          !!WS!!wi(k, 1:QLdof) = matmul(elem%w(0, (k-1)*dof +1 : k*dof), phi(1:dof, 1:QLdof) )
          wi(k, 1:QLdof) = matmul(elem%wS(k, 1 : dof), phi(1:dof, 1:QLdof) )
          !write(*,'(a6,2i5,30es12.4)') 'wi', i, dof, wi(k,1:QLdof)
       enddo

       wp(elem%face(idx, 1), 1:ndim) = wp(elem%face(idx, 1), 1:ndim) + wi(1:ndim, 1)
       wp(elem%face(idx, 2), 1:ndim) = wp(elem%face(idx, 2), 1:ndim) + wi(1:ndim, 2)
       wp(elem%face(idx, 3), 1:ndim) = wp(elem%face(idx, 3), 1:ndim) + wi(1:ndim, 3)
       supp(elem%face(idx, 1:3)) = supp(elem%face(idx, 1:3)) + 1.
    enddo


    do i=1,grid%npoin
       wp(i, 1:ndim)  = wp(i, 1:ndim) / supp(i)

       !!write(90,*) grid%x(i, 1:2), wp(i, 1:ndim)

    enddo

  end subroutine EvalSolutionVertexes

  !> smoothing of the metric stored ????????
  subroutine SmoothMetric( )
    class(element), pointer :: elem, elem1
    real :: est, estT
    real, dimension(1:2) :: grad
    real, dimension(:,:), allocatable :: rgabc
    real :: weigth
    integer :: i, j, k, ipoc
    logical :: smooth_max, smooth_max1

    smooth_max = .false.
    smooth_max1 = .false.

    ipoc = 3;   weigth = 1.0     ! ESCO scalar case
    !ipoc = 2;   weigth = 0.75     ! ESCO NS case ??
    if(state%modelName == 'NSe')  &
         !!ipoc = 3;   weigth = 0.75  ! ESCO NS case ??
         ipoc = 2;   weigth = 0.25  ! ESCO NS case ??
    !ipoc = 1;   weigth = 0.5     ! ESCO NS case ??

    if( state%space%estim_space == 'RES') & !!state%space%adapt%adapt_method == 'ANI') &
         ipoc = 2;   weigth = 0.4

    if(state%space%adapt%adapt_type == 'Ihp' ) then
       !ipoc = 2;   weigth = 1.0 ! case K
       !ipoc = 2;   weigth = 0.8 ! case K
       ipoc = 1;   weigth = 0.4  !case E, J   ! make no sense for smooth_max = .true.
       smooth_max = .true.
       !smooth_max = .true.
    endif

    if(state%space%adapt%adapt_type == 'Ahp' .and. state%space%estim_space /= 'inter' ) then
       !ipoc = 2;   weigth = 1.0 ! case K
       !ipoc = 2;   weigth = 0.8 ! case K
       ipoc = 1;   weigth = 0.4  !case E, J   ! make no sense for smooth_max = .true.
       smooth_max = .false.
       !smooth_max = .true.
    endif

    !weigth = 0.5
    !weigth = 0.25

    allocate( rgabc(1:grid%nelem, 1:4) )

    do k=1, ipoc
       rgabc(:,:) = 0.

       do i=1,grid%nelem
          elem => grid%elem(i)


          rgabc(i,1:3) = elem%rgabc(1:3)
          rgabc(i,4) = 1.


          if(smooth_max1) then

             do j=1,elem%flen
                if(elem%face(neigh,j) > 0) then
                   elem1 => grid%elem(elem%face(neigh,j))

                   rgabc(i,1:3) =  max(elem1%rgabc(1:3),  rgabc(i,1:3) )
                endif
             enddo

          else
             do j=1,elem%flen
                if(elem%face(neigh,j) > 0) then
                   elem1 => grid%elem(elem%face(neigh,j))

                   rgabc(i,1:3) = rgabc(i,1:3) + elem1%rgabc(1:3) * weigth
                   rgabc(i,4) = rgabc(i,4) + weigth
                endif
             enddo

          endif
       enddo

       do i=1,grid%nelem
          elem => grid%elem(i)
          if(smooth_max) then
             !elem%rgabc(1:3) = max(elem%rgabc(1:3),  0.6 * rgabc(i,1:3) )
             elem%rgabc(1:3) = max(elem%rgabc(1:3),  rgabc(i,1:3)/ rgabc(i, 4 ) )
          else
             elem%rgabc(1:3) = rgabc(i,1:3) / rgabc(i, 4)
          endif
       enddo
    enddo

    !do i=1,grid%nelem
    !   elem => grid%elem(i)
    !
    !   est  = elem%amaL
    !   estT = elem%amaR * est
    !   grad(1:2) = elem%amaD(1:2)
    !
    !   elem%rgabc(1) = est * grad(1)**2 + estT * grad(2)**2
    !   elem%rgabc(2) = (est - estT) * grad(1) * grad(2)
    !   elem%rgabc(3) = est * grad(2)**2 + estT * grad(1)**2
    !enddo

    deallocate( rgabc )
  end subroutine SmoothMetric



 subroutine DrawEllips(ieli, rgabc, xc)
   integer, intent(in) :: ieli
   real, dimension(1:3), intent(in) :: rgabc
   real, dimension(1:2), intent(in) :: xc
   integer :: i, num
   real :: x,y,t, xs, ys, rnorm

   !print*,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
   num = 100
   do i=0,num
      t = 1.*i/num*6.283185307
      x = cos(t)
      y = sin(t)
      rnorm = ( x*x *rgabc(1) + 2*x*y *rgabc(2) + y*y *rgabc(3) )**0.5
      xs = x/rnorm + xc(1)
      ys = y/rnorm + xc(2)
      write(ieli,*) xs,ys !, i, xi,yi, ieli
   enddo
   write(ieli, *)'## '
   write(ieli, *)

 end subroutine DrawEllips

 !> test if triangle elem is inside ellipse given by  Kparams
 subroutine TriangleInsideEllipse(elem, Kparams, ins )
   type(element), intent(inout) :: elem
   real, dimension(1:5),intent(in) :: Kparams
   logical, intent(out) :: ins
   real, dimension(:,:), allocatable :: M
   real, dimension(: ), allocatable :: xi
   real :: lam_max, lam_min, max_a, rr
   integer :: l

   allocate( M (1:2, 1:2), xi(1:2) )

   ! setting of matrix M defining the metric
   lam_max = Kparams( 1)
   lam_min = Kparams( 2)
   max_a   = Kparams( 3)

   M(1,1) =  lam_max * cos(max_a)**2 + lam_min * sin(max_a)**2
   M(1,2) = (lam_max - lam_min) * cos(max_a) * sin(max_a)
   M(2,1) = M(1,2)
   M(2,2) = lam_max * sin(max_a)**2 + lam_min * cos(max_a)**2

   ins = .true.
   ! we test if all vertices are inside ellipse
   do l=1,elem%flen
      xi(1:2) = grid%x(elem%face(idx, l), 1:2) - elem%xc(1:2)
      rr = dot_product( xi(1:2), matmul(M(1:2, 1:2), xi(1:2) ) )
      if(rr > 1) ins = .false.
   enddo

   deallocate(M, xi)


 end subroutine TriangleInsideEllipse

 !> set the Metric for isotropic refinement
  subroutine IsotropicMetricOLD( )
    class(element), pointer :: elem
    integer :: i, j, k,  imt, imt1, is, ityp, pK
    !!!real :: !tol_min, tol_max,
    real :: loc_tol, h_opt, lam_actual, lam_max
    real :: scale, weight, ratio, ratio_min
    logical :: hp_adapt
    integer :: iA, iB, iC, iD, iE, iFF, iG, iH, iI

    hp_adapt = .true.
    if(state%space%adapt%adapt_space == 'IMAh') hp_adapt = .false.  ! only h-adaptation

    !if(                    estim < tol * ratio_min ) then derefinement
    !if( tol * ratio_min <  estim < tol             ) then nothing
    ratio_min = 0.1
    if(state%time_dependent) ratio_min = 0.4

    ! (ityp==1) scale = 1, (ityp==2) scale = 1./N, (ityp==3) scale = |K|/|Omega|
    !ityp = 1

    ityp = 2
    !ityp = 3

    scale = 1.
    if(ityp == 2) scale = 1./(1.05 * grid%nelem)**0.5


     if(state%space%estim_space /= 'pNeu')  call  MarkElements( )

    iA = 0; iB =0;  iC= 0;iD = 0; iE =0;  iFF=0;  iG= 0; iH = 0; iI = 0
    do i=1,grid%nelem
       elem => grid%elem(i)

       if(ityp == 3) scale = (elem%area/state%space%domain_volume)**0.5
       loc_tol = scale *  state%space%adapt%tol_min
       if(state%time_dependent)  then
          if( state%time%disc_time /= 'STDG') then
             loc_tol = scale *  state%space%adapt%tol_min / state%time%FinTime**0.5
             print*,'Not tested FDG'
          else
             loc_tol = scale *  state%space%adapt%tol_min *(( state%time%ttime - state%time%ttime_save )/ state%time%FinTime)**0.5
          endif
       endif



       !print*,'########$$$',loc_tol

       !lam_actual = 3./elem%diam**2

       elem%hsplit = 0
       elem%psplit = 0

       ratio = 1.

       if(.not. hp_adapt .and. state%time_dependent) then
          ratio = (elem%estimST/ loc_tol )**(1./elem%deg)

       else
          !write(*,'(a6,i5,6es12.4)') 'EDWS',elem%i, &
          !     elem%estimST, loc_tol, elem%reg, elem%regT0 , elem%regT2

          if(elem%estimST > loc_tol  ) then

             if(elem%reg <  elem%regT0 .and. elem%deg < MaxDegreeImplemented - 2 &

                  .and. hp_adapt) then
                ! p-refinement
                elem%psplit = 1
                ratio = 0.75

                iC = iC + 1
             else
                ! h-refinement
                !elem%hsplit = 1
                ratio = (elem%estimST/ loc_tol )**(1./elem%deg)

                if(elem%reg >  elem%regT2 .and. elem%deg > 2 .and. hp_adapt) then
                   ! p-derefinement
                   !elem%psplit = -1
                   !ratio = (elem%estimST/ loc_tol )**(1./elem%deg-1)
                   iI = iI + 1
                else
                   iFF = iFF + 1
                endif
             endif

          else if(elem%estimST < ratio_min * loc_tol  ) then
             !elem%hsplit = -1
             ! h-derefinement
             ratio = (elem%estimST/ (0.1 * loc_tol) )**(1./elem%deg)

             if( elem%reg <  elem%regT0 .and. hp_adapt) then

                if(elem%deg < MaxDegreeImplemented - 2) then

                   ! p-refinement
                   !ratio = ratio * 0.5   ! 0.75
                   !elem%psplit = 1

                   iA = iA + 1
                endif

                !write(*,'(a5,i5,8es12.4)') &
                !     'Eppp',elem%i, elem%estimST, 0.1 * loc_tol, ratio, elem%reg ,  elem%regT0
             elseif( elem%reg >  elem%regT2 .and. hp_adapt) then
                iG = iG + 1
                elem%psplit = -1
                !ratio = (elem%estimST/ 0.1 * loc_tol )**(1./elem%deg-1)

             else
                iD = iD + 1

                !write(*,'(a5,i5,8es12.4)') &
                !     'ESDF',elem%i, elem%estimST, 0.1 * loc_tol, ratio, elem%reg ,  elem%regT0
             endif
          elseif (hp_adapt) then
             !ratio = (elem%estimST/ (0.85 * loc_tol) )**(1./elem%deg)
             if( elem%reg <  elem%regT0 .and. elem%deg < MaxDegreeImplemented - 3) then
                ! h->p substitution
                !ratio = ratio * 0.5 !0.75
                !elem%psplit = 1

                iB = iB + 1


             elseif(  elem%reg >  elem%regT2 .and. elem%deg > 2) then
                !print*,'p-DEFRE'
                !elem%psplit = -1
                !ratio = (elem%estimST/ loc_tol )**(1./elem%deg-1)

                iH = iH + 1
             else
                iE = iE + 1


             endif

          endif

       endif


       ratio = min(ratio, 5.0)
       !ratio = ratio**1.2

       !if(elem%hsplit /= 0 .or. elem%psplit /= 0) &
       !     write(*,'(3es12.4, 2i5, 3es12.4)') elem%xc(:), ratio, elem%hsplit, elem%psplit, &
       !     elem%estimST, tol_max *scale

       h_opt = elem%diam / ratio

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
    end do


    write(*,'(a10,10i7)') 'AdatG:',iA, iB, iC, iD, iE, iFF, iG, iH, iI
    if(state%space%adapt%adapt_level == 0) &
         write(84,'(a10,10a7)') '***:','iA', ' iB', ' iC', ' iD', ' iE', ' iFF', ' iG', ' iH', ' iI'
    write(84,'(a10,10i7)') 'AdatG:',iA, iB, iC, iD, iE, iFF, iG, iH, iI

    !print*
    !print*,'Sum', sum(grid%elem(:)%estimST**2)**0.5, tol_max
    !print*


  end subroutine IsotropicMetricOLD

 !> set the Metric for isotropic refinement
  subroutine IsotropicMetric( )
    class(element), pointer :: elem
    integer :: i, j, k,  imt, imt1, is, ityp, pK
    !!!real :: !tol_min, tol_max, 
    real :: loc_tol, h_opt, lam_actual, lam_max
    real :: scale, weight, ratio, deref_level
    logical :: hp_adapt
    integer :: iA, iB, iC, iD, iE, iFF, iG, iH, iI
    character(len=15) :: file1, file2
    character(len=5) :: ch5
    real :: max_ratio, min_ratio, sum
    logical :: Algol1, Algol3
 
    !Algol1 = .true.
    Algol1 = .false.

    !Algol3 = .true.
    Algol3 = .false.

    !!! if (Algol1 == .false.  .and. Algol3 == .false. ) THEN  Algol2

    !max_ratio = 4.0 
    !max_ratio = 2.0 ! case E, J   ! originally 2.5 for cas J
    max_ratio = 2.0    ! maximal ratio of refinement
    min_ratio = 0.  !0.9    ! maximal ratio of DErefinement

    hp_adapt = .true.
    if(state%space%adapt%adapt_space == 'IMAh') hp_adapt = .false.  ! only h-adaptation

    !if(                     estim < tol * deref_level ) then derefinement
    !if( tol * deref_level < estim < tol             ) then nothing

    !deref_level = 1.0  ! 0.85  !-0.1  !0.2
    !deref_level = 0.1  !case E, J
    deref_level = 0.1  

    if(Algol3) deref_level = 1.0  ! test for HO_rec  

    if(Algol1) deref_level = -0.1

    if(state%time_dependent) deref_level = 0.4

    ! (ityp==1) scale = 1, (ityp==2) scale = 1./N, (ityp==3) scale = |K|/|Omega|
    !ityp = 1
    ityp = 2  
    !ityp = 3   

    scale = 1.
    if(ityp == 2) scale = 1./(1.05 * grid%nelem)**0.5


    ! It is necessary ????
    !if(state%space%estim_space /= 'pNeu' .and. state%space%estim_space /= 'HO_rec') call  MarkElements( )

    iA = 0; iB =0;  iC= 0;iD = 0; iE =0;  iFF=0;  iG= 0; iH = 0; iI = 0
    do i=1,grid%nelem
       elem => grid%elem(i)

       if(ityp == 3) scale = (elem%area/state%space%domain_volume)**0.5
       loc_tol = scale *  state%space%adapt%tol_min
       if(state%time_dependent)  then
          if( state%time%disc_time /= 'STDG') then
             loc_tol = scale *  state%space%adapt%tol_min / state%time%FinTime**0.5
             print*,'Not tested FDG'
          else
             loc_tol = scale *  state%space%adapt%tol_min *(( state%time%ttime - state%time%ttime_save )/ state%time%FinTime)**0.5
          endif
       endif


       ! we decide p-refinement, p-kept, p-derefinement,  one-parameter variant
       elem%psplit = 0
       ratio = 1.

       if(.not. hp_adapt .and. state%time_dependent) then
          ratio = (elem%estimST/ loc_tol )**(1./elem%deg)

       else
          !write(*,'(a6,i5,6es12.4)') 'EDWS',elem%i, &
          !     elem%estimST, loc_tol, elem%reg, elem%regT0 , elem%regT2

          !sum = sum + 
          !write(50+state%space%adapt%adapt_level) sum

          if(elem%estimST > loc_tol  ) then
             
             if(elem%reg < elem%regT0 .and. elem%deg < MaxDegreeImplemented - 2) then !p-refinement
                elem%psplit = 1
                iA = iA + 1
                !ratio = 0.75
                !ratio = (elem%estimST/ loc_tol )**(1./(elem%deg+1) )
                
                !write(100 + 10*state%space%adapt%adapt_level + 1, '(3es14.6,i5,10es14.6)') & 
                !     elem%xc(:),ratio,elem%psplit, &
                !     elem%estimST/ loc_tol, (elem%estimST/ loc_tol )**(1./elem%deg), &
                !     elem%reg / elem%regT0 , elem%reg / elem%regT2

             else 
                ! h-refinement
                ratio = (elem%estimST/ loc_tol )**(1./elem%deg)
                !ratio = ratio 
                ratio = ratio * 1.25  !1.25  !case E, J, K

                !WSW!
                if(Algol1) ratio = 2.
                if(Algol3) ratio = (elem%estimST/ loc_tol )**(1./elem%deg)

                iB = iB + 1
                !write(100 + 10*state%space%adapt%adapt_level + 2,  '(3es14.6,i5,10es14.6)') & 
                !     elem%xc(:),ratio,elem%psplit, &
                !     elem%estimST/ loc_tol, (elem%estimST/ loc_tol )**(1./elem%deg), &
                !     elem%reg / elem%regT0 , elem%reg / elem%regT2

             endif
             
          elseif(elem%estimST <= deref_level * loc_tol  ) then

             !if(Algol3) deref_level = deref_level / 2.

             if( elem%reg < elem%regT0 .and. elem%deg < MaxDegreeImplemented - 2) then
                elem%psplit = 1    !p-refinement
                ratio = (elem%estimST/ (deref_level *loc_tol) )**(1./(elem%deg+1) )

                iC = iC + 1
                !write(100 + 10*state%space%adapt%adapt_level + 3,  '(3es14.6,i5,10es14.6)') & 
                !     elem%xc(:),ratio,elem%psplit, &
                !     elem%estimST/ loc_tol, (elem%estimST/ loc_tol )**(1./elem%deg), &
                !     elem%reg / elem%regT0 , elem%reg / elem%regT2

             elseif( elem%reg > elem%regT2   .and.  elem%deg > 1) then 
                ! p-derefinement
                elem%psplit = -1 
                ratio = (elem%estimST/ (deref_level * loc_tol) )**(1./(elem%deg-1))

                iD = iD + 1
                !write(100 + 10*state%space%adapt%adapt_level + 4,  '(3es14.6,i5,10es14.6)') & 
                !     elem%xc(:),ratio,elem%psplit, &
                !     elem%estimST/ loc_tol, (elem%estimST/ loc_tol )**(1./elem%deg), &
                !     elem%reg / elem%regT0 , elem%reg / elem%regT2

             else
                 ! h-(de)refinement
                ratio = (elem%estimST/ (deref_level * loc_tol) )**(1./elem%deg)
                iE = iE + 1
                !write(100 + 10*state%space%adapt%adapt_level + 5,  '(3es14.6,i5,10es14.6)') & 
                !     elem%xc(:),ratio,elem%psplit, &
                !     elem%estimST/ loc_tol, (elem%estimST/ loc_tol )**(1./elem%deg), &
                !     elem%reg / elem%regT0 , elem%reg / elem%regT2

             endif
             
          else
             iFF = iFF + 1
             !write(100 + 10*state%space%adapt%adapt_level + 6,  '(3es14.6,i5,10es14.6)') &
             !     elem%xc(:),ratio,elem%psplit, &
             !        elem%estimST/ loc_tol, (elem%estimST/ loc_tol )**(1./elem%deg), &
             !        elem%reg / elem%regT0 , elem%reg / elem%regT2
          endif
       endif   ! if(.not. hp_adapt 

       !if( dot_product(elem%xc(1:2), elem%xc(1:2) )**0.5 < 1E-3) then
       !   write(77,'(a4,3i5,10es12.4)') 'IAMA',state%space%adapt%adapt_level, elem%i, elem%deg, &
       !        elem%estimST, loc_tol, elem%reg, ratio, min(ratio, 2.5)**1.2
       !endif

       ratio = min( ratio, max_ratio)
       if(Algol3) ratio = max(ratio, min_ratio)

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
       !write(*,'(a6,i5,8es12.4)') '#ISO#',elem%i,elem%rgabc(:), lam_max, lam_actual,elem%estimST
    end do


    ! write(*,'(a10,10i7)') 'AdatG:',iA, iB, iC, iD, iE, iFF, iG, iH, iI
    ! if(state%space%adapt%adapt_level == 0) then
    !     write(84,'(a10,10a12)') '***:',&
    !          'E>T p+', ' E>T h+',' E<T p+',' E<T p-', 'E<T p0', ' ', ' ', ' ', ''
    !     write(84,'(a10,10a12)') '***:','iA', 'iB', 'iC','iD','iE', 'iFF', 'iG ', 'iH ', ' iI'
    !  endif
    
    ! write(84,'(a10,10i12)') 'AdatG:',iA, iB, iC, iD, iE, iFF, iG, iH, iI

    !print*
    !print*,'Sum', sum(grid%elem(:)%estimST**2)**0.5, tol_max
    !print*

    file1 = 'metrixI00000'
    file2 = 'metrixS00000'

    is = 0
    if(state%space%adapt%adapt_level > 0) is = int(log(1. * state%space%adapt%adapt_level)/log(10.)) 

    write( ch5, '(i5)' ) state%space%adapt%adapt_level  ! change the format if num_size /= 5 !!!
    !!!write( ch5, '(i5)' ) state%time%iter  ! change the format if num_size /= 5 !!!
    file1(12-is: 12)  = ch5(5-is:5)
    file2(12-is: 12)  = ch5(5-is:5)

    ! !variant of high order Riemann metric
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

    call SmoothMetric( ) 


    ! ! variant of high order Riemann metric
    ! imt = 24
    ! open(imt, file=file2, status='UNKNOWN')
    ! do i=1,grid%nelem !,2
    !    elem => grid%elem(i)
    !    !if( dot_product(elem%xc(:), elem%xc(:))**0.5 < 4E-2) &
    !    call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )
    ! enddo
    ! close(imt)


  end subroutine IsotropicMetric



end module ama_hp_interpol
