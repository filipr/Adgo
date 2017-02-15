!> parameters for the ama hp adaptation
module ama_hp_interpol_params
  real :: ama_err_max
  real :: ama_err_min
  real :: ama_err_aver
  real :: ama_err_total
  real :: ama_err_total_reduced
  real :: ama_err_total_old
  real :: ama_target_tol

  real, parameter, private :: ama_r_max = 4.0   ! the maximal allowed reduction of the element area
  real, parameter, private :: ama_c_max = 4.0   ! the maximal allowed increase  of the element area
  integer, dimension(:), allocatable :: ama_iest
  real, dimension(:), allocatable :: ama_est

end module ama_hp_interpol_params

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
  use regularity

  implicit none


  public:: AnisotInterpolEstimates
  public:: SetQuantities4Metric
  public:: Eval_hp_Metric
  public:: Eval_hp_MetricElem
  public:: Set_hp_metric_hpREZ
  public:: Set_hp_metric_Inter
  public:: Eval_All_Derivatives
  public:: EvalSolutionVertexes
  public:: FindAnisotropyEIp
  public:: FindAnisotropyEIp_H1
  public:: SmoothMetric

  public:: TriangleInsideEllipse

  public:: IsotropicMetric
  public:: IsotropicMetric_SimpleOrdering
contains

  !> perform the anisotropic error estimates using interpolation error estimates
  subroutine AnisotInterpolEstimates( )
    class(element), pointer :: elem
    !real, dimension(:,:), allocatable :: wp ! array  with solution in vertexes
    integer:: ndimL   ! number of quantities for metric evaluation
    integer :: i, j, k,  imt, imt1, is
    logical :: loc_implicitly
    character(len=15) :: file1, file2, file3
    character(len=5) :: ch5

    if(nbDim /=2 ) then
       print*,' ## AnisotErrorEstimates only for nbDim = 2!'
       stop
    endif


    !allocate( wp(1:grid%npoin, 1:ndim) )
    !
    !call EvalSolutionVertexes(wp(1:grid%npoin, 1:ndim) )

    if(.not. grid%ElemSupports) &
         call SeekElemSupports(grid)  ! create the list of elements sharing at least a vertex with elem

    call SetQuantities4Metric(ndimL )


    ! New variant of the HO recontruction
    ! it allocates elem%wSD, it is necessary to DEALLOCATE
    !call Compute_means_values_DW(grid, ndimL, 1)

    !stop

    state%err(algeb) = 0.

    state%num_limits = 0
    call Eval_hp_Metric( ndimL )  !!!wp(1:grid%npoin, 1:ndim) )

    !deallocate(wp)

    file1 = 'metrixA00000'
    file2 = 'metrixS00000'
    file3 = 'meshixS00000'

    is = 0
    if(state%space%adapt%adapt_level > 0) is = int(log(1. * state%space%adapt%adapt_level)/log(10.))

    write( ch5, '(i5)' ) state%space%adapt%adapt_level  ! change the format if num_size /= 5 !!!
    !!!write( ch5, '(i5)' ) state%time%iter  ! change the format if num_size /= 5 !!!
    file1(12-is: 12)  = ch5(5-is:5)
    file2(12-is: 12)  = ch5(5-is:5)
    file3(12-is: 12)  = ch5(5-is:5)


    ! variant of high order Riemann metric
   !  imt = 24
   !  open(imt, file=file1, status='UNKNOWN')
   !  do i=1,grid%nelem
   !    elem => grid%elem(i)
   !    !if( sqrt(dot_product(elem%xc(:), elem%xc(:)) ) < 0.1) &
   !    ! ! !    if(abs(elem%xc(1) - 1.5) < 0.25 .and. elem%xc(2) > 1.75 ) &

   !    !if(elem%xc(1) > 0.30 .and. elem%xc(1) < 0.38 .and.  &
   !    !     elem%xc(2) > -30 .and. elem%xc(2) < 0.1) &

   !    if(elem%xc(1) > -1.30 .and. elem%xc(1) < - 0.5 .and.  &
   !         elem%xc(2) > 0.5 .and. elem%xc(2) < 10.1) &
   !    !if( mod(i, 3) == 1) &
   !         call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )

   ! enddo
   ! close(imt)

    !    print*,'STOPPED IN ESRTWGQ'
    !    stop

    if(state%modelName == 'pedes' ) call Refine_a_priori( )

    call SmoothMetric( )


    ! ! variant of high order Riemann metric
     imt = 24
     open(imt, file=file2, status='UNKNOWN')
     do i=1,grid%nelem !,2
        elem => grid%elem(i)
        !if( sqrt(dot_product(elem%xc(:), elem%xc(:)) ) < 0.1) &
             call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )
     enddo
     close(imt)

     ! mesh
     open(imt, file=file3, status='UNKNOWN')
     do i=1,grid%nelem !,2
        elem => grid%elem(i)
        write(imt, '(2es14.6)') grid%x(elem%face(idx, 1), 1:2)
        write(imt, '(2es14.6)') grid%x(elem%face(idx, 2), 1:2)
        write(imt, '(2es14.6)') grid%x(elem%face(idx, 3), 1:2)
        write(imt, '(2es14.6)') grid%x(elem%face(idx, 1), 1:2)
        write(imt,'(x)')
     enddo
     close(imt)

    ! deallocation of arrays allocated in SeekElemSupports
    do i=1,grid%nelem
       elem => grid%elem(i)
       !write(*,'(a6, i5, 6es12.4)') 'ede3',i, elem%xc, elem%rgabc
       !deallocate(elem%supp)
       deallocate(elem%wS)
    enddo

    ! !print*,'stoped in NEW ama-hp_interpol.f90',  state%err(algeb)
    ! !stop

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
    if(state%time_dependent) ndimL = 2
    if(state%space%adapt%max_adapt_level == 0) ndimL = 1  ! necessary for ndimL = 2

    !ndimL = 3
    !!!ndimL = grid%curved_deg

    if(ndimL > 1) then
       do i=1,grid%nelem
          elem => grid%elem(i)
          allocate(elem%wSS(1:1, 1:1,1:elem%dof*ndim ) )
       enddo


       ! recomputation of the stored solution in gridS to the actual grid
       call AdvancedInterpolDGsolution(grid, gridS, 1 )
    endif

    !write(241, *) grid%nelem, ndim, state%time%ttime, state%time%tau(1), state%time%iter ! new
    !write(242, *) grid%nelem, ndim, state%time%ttime, state%time%tau(1), state%time%iter ! new

    do i=1,grid%nelem
       elem => grid%elem(i)
       !stop 'allocate elem%wS'
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

       ! if(elem%i == 1) print*,'fort.11 written HEDt53d'
       ! if(elem%xc(1) < 5.) &
       ! call PlotElemFunction3D(11, elem,  elem%dof, elem%wS(1, 1:elem%dof) )

    enddo


    !write(241,*) 1., 1., 1., 0., 1.
    !write(242,*) 1., 1., 1., 0., 1.

  end subroutine SetQuantities4Metric

  !> evaluate the Riemann metric for hp-mesh
  subroutine Eval_hp_Metric( ndimL )    !!!wp)
    use ama_hp_interpol_params
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

    !evaluate of the parameters for the ama_hp_adaptation
    if(state%space%estim_space == 'RES' .or. state%space%estim_space == 'DWR' .or. &
         state%space%estim_space == 'pNeu' ) call Eval_AMA_paramets( )

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

    ! setting, if next adaptation should be carried out
    if(  state%space%estim_space == 'inter') then !!state%space%adapt%adapt_method == 'Ahp') then
       ! H1-seminorm
       if(Lq <= -0.99) then
          if(state%err(interH1) <= state%space%adapt%tol_max ) state%space%adapt%stop_adaptation = 1
          !if(state%err(interH1) > interOLD .and. state%err(interH1) < 1.1 * interOLD) &
          !     state%space%adapt%stop_adaptation = -1

       ! Lq-norm
       elseif(Lq >= 1.) then
          if(state%err(interLq) <= state%space%adapt%tol_max ) state%space%adapt%stop_adaptation = 1
          !if(state%err(interLq) > interOLD .and. state%err(interLq) < 1.1 * interOLD) &
          !     state%space%adapt%stop_adaptation = -1

       ! L^infty norm
       else
          if(state%err(interL8) <= state%space%adapt%tol_max ) state%space%adapt%stop_adaptation = 1
          !if(state%err(interL8) > interOLD .and. state%err(interL8) < 1.1 * interOLD) &
          !     state%space%adapt%stop_adaptation = -1
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
    elem%wSS = 0.

    call Eval_All_Derivatives(elem, ndimL)

    !print*,'....', elem%i
    ! does not work accuratly, Eval_All_Derivatives seems to be better
    ! 21 / 11 / 2016 - not nice grid observed for porous media, 
    !strong local refinement, bad angles
    !call Eval_All_Derivatives_by_Projection(elem, ndimL)
    !print*,'...B', elem%i


    elem%interLq = 0.
    elem%interL8 = 0.
    elem%interH1 = 0.

    if(state%space%estim_space == 'inter' .or.  state%space%estim_space == 'HO_rec') then
       call Set_hp_metric_Inter( elem, ndimL )

    else if(state%space%estim_space == 'RES' .or. state%space%estim_space == 'DWR' .or. &
         state%space%estim_space == 'pNeu' ) then

       call Set_hp_metric_hpREZ( elem, ndimL )
       !call Set_hp_metric_Inter( elem, ndimL )

       !if(elem%i == 1) print*,'Interpol subroutine is CALLED !!!!'
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
    real :: lam_max, lam_min, max_a, lam_actual, lam_actualT
    !real :: diam, diamT, elem_estim, est0, est0T,
    real :: h_min, h_max, lambda_min, lambda_max, lam_max_actual_ratio
    real, dimension(:,:), allocatable :: Kparams
    real :: ratio, epsilon1, ropt !, tol_loc
    integer :: ideg, degP, i, j,  pK, ip, iopt, ideg_min, ideg_max
    real :: weight, regularity, pKdiff
    logical :: ins, singularity
    integer :: deg_opt

    ! approach based on the Riemann metric generated by high order interpol error estimate
    !tol_loc = state%space%adapt%tol_min

    !tol_loc = state%space%adapt%tol_min / grid%nelem**0.5

    epsilon1 = 1E+12
    h_max = min( state%space%diam / 4., 4.)
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

    call SetOptimalAnisotropyEIp( elem, ndimL, iopt, Kparams(0:3, 1:5)  )


   ! modification of optimal  polynomial degree using information of the singularity
    if(state%space%adapt%adapt_space /= 'AMAh' .and. state%space%adapt%adapt_space /= 'IMAh')  then

       ! detection of the singularity
       singularity = .false.
       !call Detect_apriori_known_singularity(elem, singularity)

       ! estimate of the local reglarity of the solution
       regularity = -1.
       call ElementEdgeJumpsAllProj(elem, regularity)

       ! the resulting polynomial degree
       pK = elem%deg + iopt - 1

       ! change of the optimal meric based on the estimate of the regularity
       if(regularity > 0) then
          pKdiff = pK + 1  - regularity - 2

          if(pKdiff  > 1.) then
             iopt = 0
             write( (state%space%adapt%adapt_level+1)*10 + 1, *) elem%xc,  pKdiff, pK, regularity
          elseif(pKdiff  > 0.) then

             if(iopt == 2) then
                iopt = 1
                write( (state%space%adapt%adapt_level+1)*10 + 2, *) elem%xc,  pKdiff, pK, regularity
             endif

          endif

       endif
    endif  !if(state%space%adapt%adapt_space /= 'AMAh' .and. state%space%adapt%adapt_space /= 'IMAh')

    ! setting of the p-adaptation
    elem%psplit = iopt - 1


    !! detection of the singularity
    !singularity = .false.
    !!call Detect_apriori_known_singularity(elem, singularity)
    !if(singularity) then
    !   deg_opt = 2   ! has to be changed !!
    !
    !   iopt = max( 0, deg_opt -  elem%deg + 1)   ! choice of the tested metric
    !   elem%psplit =  deg_opt - elem%deg         ! setting of psplit (chnage of elem%deg)
    !   if(elem%psplit < 0) Kparams(iopt, 1:2) =  Kparams(iopt, 1:2) * 2.
    !endif



    ! setting the metric from the best candidate: pK-1, pK, pK+1, averaging with pK
    !!weight = 0.9
    !!weight = 0.75
    !weight = 0.85
    weight = 1.0
    !if(state%modelName == 'scalar' .or.state%modelName == '2eqs') weight = 0.5

    lam_max = weight * Kparams(iopt, 1) + (1. - weight) * Kparams(1, 1)
    lam_min = weight * Kparams(iopt, 2) + (1. - weight) * Kparams(1, 2)
    max_a   = Kparams(iopt, 3)

    elem%ama_p = weight * elem%psplit

    if( state%space%adapt%adapt_space == 'IMAh' .or. state%space%adapt%adapt_space == 'IMAhp') then
       lam_min = lam_max
       max_a = 0.
    endif



    !if(elem%deg + elem%ama_p <  1) then
    !   write(*,'(a20,2es12.4,2i6,4es12.4)') &
    !        '##D#DE#D#DSW', elem%deg + elem%ama_p, elem%ama_p, iopt, elem%deg, elem%xc
    !endif

    ! limitation from the point of view of reasonability, too large and too small elements
    lam_max = max(lambda_min, min(lambda_max, lam_max) )
    lam_min = max(lambda_min, min(lambda_max, lam_min) )

    !if(elem%i <=10) write(*,'(a6, 6es12.4)') 'lim 1:', lam_max,  lam_min, 1./sqrt(lam_max),1./sqrt(lam_min)

    ! THE OPTIMAL METRIC
    !lam_max = Kparams(3, 1)
    !lam_min = Kparams(3, 2)
    !max_a =   Kparams(3, 3)
    !elem%ama_p = Kparams(3, 4)
    !!!par = Kparams(3, 5)

    !lam_actual = 1. / elem%diam**2.0
    lam_actual = (2*elem%area/elem%diam)**(-2)
    lam_actualT = (elem%diam)**(-2)
    ! limitation of the aspect ratio
    ratio = lam_min / lam_max

    !if(elem%xc(2) > 0 .and. elem%xc(2) < 1. .and. elem%xc(1) > 0.5 .and. elem%xc(1) < 0.8) &
    !if( elem%xc(1) > 19 .and. elem%xc(1) < 21) &
    !if( elem%xc(1) < 1 ) &
    !if(ratio < 2E-5) &
    !     write(*,'(a5,2i5,8es12.4)')  '!!ded!', elem%i, degP, ratio, lam_min, lam_max, max_a

    !ratio = max(ratio, 1E-4)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    ratio = max(ratio, 2E-5)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    !ratio = max(ratio, 2.5E-4)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    !!ratio = max(ratio, 1E-3)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    !!ratio = max(ratio, 1E-2)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    !ratio = max(ratio, diamT/ diam *0.2)   ! ratio may increase maximally 5-times
    if(state%modelName == 'scalar' .or.state%modelName == '2eqs' .or.state%modelName == 'porous' ) &
         ratio = max(ratio, 1E-3)


     ! if( abs(elem%xc(1)) < 0.05 .and. abs(elem%xc(2)) <  0.05) then
     !    write(*,'(a6,2(3es12.4,a2), 3es12.4)') '####',  &
     !         lam_max, lam_min, ratio,'|', &
     !         lam_actual, lam_actualT, lam_actualT / lam_actual, '|',lam_actualT / lam_actual / ratio
     ! endif

    ! limitation of the refinement at one level
    !lam_max_actual_ratio = 10  !4. !10
    lam_max_actual_ratio = 4. !10
    !if(.not. state%time_dependent) lam_max_actual_ratio = 50.

    lam_max = min(lam_max, lam_actual * lam_max_actual_ratio)  ! at most r^{1/2} times smaller element
    lam_max = max(lam_max, lam_actual / lam_max_actual_ratio)  ! at most r^{1/2} times bigger element

    !if(elem%i <=10) write(*,'(a6, 6es12.4)') 'lim 2:', lam_max,  lam_min, 1./sqrt(lam_max),1./sqrt(lam_min)

    ! limitation of the ratio
    ratio = max(ratio,  lam_actualT / lam_actual / 2)


    ! new lam_min after limitation of the aspect ratio and lam_max
    lam_min = lam_max * ratio

    !if(elem%i <=10) write(*,'(a6, 6es12.4)') 'lim 3:', lam_max,  lam_min, 1./sqrt(lam_max),1./sqrt(lam_min)

    !if(elem%i < 5) print*,'ATTENTION 3d5edeud3hd38de3'
    !lam_max = lam_actual * 2
    !lam_min = lam_actualT * 2
    !elem%ama_p = 0.

    !if( abs(elem%xc(1)) < 0.05 .and. abs(elem%xc(2)) <  0.05) then
    !   write(*,'(a6,2(3es12.4,a2), 3es12.4)') '####',  &
    !        lam_max, lam_min, ratio,'|', &
    !        lam_actual, lam_actualT, lam_actualT / lam_actual, '|',lam_actualT / lam_actual / ratio
    !   print*,'--------'
    !endif


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




  !> evaluation of the parameters for the ama_hp_adaptation
  subroutine Eval_AMA_paramets( )
    use ama_hp_interpol_params
    real :: ratio 

    ama_err_max = maxval(grid%elem(:)%estim_loc)
    ama_err_min = minval(grid%elem(:)%estim_loc)
    ama_err_aver = sum(grid%elem(:)%estim_loc) / grid%nelem

    ama_err_total = sqrt( dot_product(grid%elem(:)%estim_loc, grid%elem(:)%estim_loc) )
    ama_err_total_reduced = 0.

    ratio = 10
    if(state%space%adapt%adapt_level == 0) then
       !ama_target_tol = state%space%adapt%tol_max
       ama_target_tol =  ama_err_total / ratio
    else
       !target_tol = ama_err_total / 10.
       !target_tol = min( state%space%adapt%tol_max, ama_err_total / 5.)
       if(ama_err_total / ratio < ama_target_tol ) then
          ama_target_tol = ama_err_total / ratio
       endif
       if(ama_err_total / ama_err_total_old > 0.9) &
            ama_target_tol = ama_target_tol * 0.5
    endif


    ! setting of the target tolerance
    print*,'#### taget tol = state%space%adapt%tol_min !!!'
    if(state%time_dependent) then
       ama_target_tol = state%space%adapt%tol_min * &
            sqrt(( state%time%ttime - state%time%ttime_save ) /state%time%FinTime)
    else
       ama_target_tol = state%space%adapt%tol_min

    endif


    open(91, file='Eval_AMA_paramets_history', status='UNKNOWN', position='append')
    write(91,'( i5, 7es12.4,2i8)') state%space%adapt%adapt_level, &
         ama_err_min, ama_err_aver, ama_err_max, ama_err_total, &
         ama_target_tol, state%space%adapt%tol_min, state%space%adapt%tol_max,&
         grid%nelem, state%nsize
    close(91)

    ama_err_total_old = ama_err_total

  end subroutine Eval_AMA_paramets


  !> evaluation of the parameters for the ama_hp_adaptation
  subroutine Sort_maximal_elements( )
    use ama_hp_interpol_params
    integer :: i

    allocate(ama_iest(1:grid%nelem), ama_est(1:grid%nelem) )
    do i=1,grid%nelem
       ama_iest(i) = i
       ama_est(i) = grid%elem(i)%estim_loc

    enddo

    call order_estims(grid%nelem, ama_iest(1:grid%nelem), ama_est(1:grid%nelem))


  end subroutine Sort_maximal_elements


  subroutine Sort_maximal_elements_dealloc( )
    use ama_hp_interpol_params
    deallocate(ama_iest, ama_est)
  end subroutine Sort_maximal_elements_dealloc

  !>  evaluation of metric from high order derivatives of degree elem%deg+1
  !> using info from the residuall error estimates
  subroutine Set_hp_metric_hpREZ( elem, ndimL)
    use ama_hp_interpol_params
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    real :: lam_max, lam_min, max_a, lam_actual, lam_actualT
    !real :: diam, diamT, elem_estim, est0, est0T,
    real :: h_min, h_max, lambda_min, lambda_max, lam_max_actual_ratio
    real, dimension(:,:), allocatable :: Kparams
    real :: ratio, epsilon1, ropt !, tol_loc
    integer :: ideg, degP, i, j,  ji, pK, ip, iopt, ideg_min, ideg_max, imt
    real :: weight, Lq, pKdiff, regularity, factorX
    logical :: ins, singularity, iprint
    real :: area, factor, tols, area_new
    integer :: deg_opt

    iprint = .false.
    !if(elem%xc(1) > 0.5 .and. elem%xc(1) < 0.7  .and. &
    !     elem%xc(2) > 0. .and. elem%xc(2) < 0.3  )iprint = .true.


    ! ama_target_tol given in Eval_AMA_paramets( )

!    if( state%space%estim_space == 'DWR') elem%estim_loc = elem%eta( dwrS, 1)

    ! error equidistribution
    !print*,'##DE##',/state%space%adapt%Lq
    if(abs(-1. - state%space%adapt%Lq) < 1E-5) then
       Lq = 0.5  ! H^1 seminorm
    elseif(abs(state%space%adapt%Lq) < 1E-5) then
       Lq = 0.  ! L^\infty norm
    else
       Lq = 1./state%space%adapt%Lq   ! L^q norm
    endif

    !tols = ama_target_tol * (elem%area/  state%space%domain_volume)**Lq
    tols = ama_target_tol * sqrt(1./grid%nelem)   ! equidistribution
    if(elem%i == 1) print*,'tols: =  ', tols, Lq
    
    factor = elem%estim_loc / tols
    factor = factor**(2./(max(1E+00, 1.*elem%deg+1)))
    
    if((elem%i == 1 .or. elem%i == grid%nelem) ) & !.and.  target_tol < state%space%adapt%tol_max) &
         write(*,'(a12,4es12.4,a2,8es12.4)') &
         'AMA_pars ii:', ama_err_min, ama_err_aver, ama_err_max,ama_err_total,'|', &
         ama_target_tol, state%space%adapt%tol_max
    
    !write(100+state%space%adapt%adapt_level,'(15es12.4)') &
    !     elem%xc(:), tols, elem%estim_loc, elem%estim_loc / tols, factor

    factorX = factor
    if(state%time_dependent) then
       factor = max(0.01, min(factor, 4.) )
       
    else
       factor = max(0.1, min(factor, 4.) )
    endif

    area = elem%area / factor
    !area = elem%diam**2 / factor

    !epsilon1 = 1E+15
    epsilon1 = 1E+10

    h_max = min( state%space%diam / 4., 4.)
    !h_max = state%space%diam / 2.


    h_min = h_max / epsilon1**0.5
    lambda_min = 1./h_max**2
    lambda_max = 1./h_min**2


    !if(factor > 1.) &
    !     write(*,'(a8, 3i5, 20es12.4)') 'YUTelem',elem%i, elem%deg, elem%psplit, elem%rgabc(:), factor,factorX
    !if(elem%i <= 3) write(*,*) 'a testing of the variant mixing HO_rec, REZ and inter, comments !EDF'
    !EDF
    elem%psplit = 0

    ! Kparams: parameters of the proposed metric, first index is the degree of interpolation,
    ! second index: 1= lambda_max, 2=lambda_min, 3= max_a, 4=max_f (maximal interp error)
    ! 5=dof/area
    allocate(Kparams(0:3, 1:5) )
    Kparams = 1.E+20

    !write(10+10*state%space%adapt%adapt_level + 5,'(300es12.4)') &
    !     elem%xc, factor, area, elem%area 

    call SetOptimalAnisotropy_AreaFixed( elem, ndimL, iopt, area, Kparams(0:3, 1:5)  )


    !EDF  iopt = elem%psplit + 1
    if(iprint) then
       do j=0, 2
          write(*,'(a4,4i5,8es12.4)') &
               'DD:',elem%i, elem%deg, iopt, j,elem%estim_loc / tols,  factorX, factor, Kparams(j, 1:5)
          write(150+state%space%adapt%adapt_level, '(2es12.4, 4i5,8es12.4)') &
               elem%xc,elem%i, elem%deg, iopt, j,elem%estim_loc / tols,  factorX, factor, Kparams(j, 1:5)
       enddo
       print*
       write(150+state%space%adapt%adapt_level, '(x)' )
    endif


    ! modification of optimal  polynomial degree using information of the singularity
    if(state%space%adapt%adapt_space /= 'AMAh' .and. state%space%adapt%adapt_space /= 'IMAh')  then

       ! detection of the singularity
       singularity = .false.
       !if(elem%i == 1) print*,'Detect_apriori_known_singularity activated !!!'
       !call Detect_apriori_known_singularity(elem, singularity)

       ! estimate of the local reglarity of the solution
       regularity = -1.
       !call ElementEdgeJumpsAllProj(elem, regularity)
       !elem%reg2 = regularity

       ! the resulting polynomial degree
       pK = elem%deg + iopt - 1

       ! change of the optimal meric based on the estimate of the regularity
       if(regularity > 0) then
          !pKdiff = pK + 1  - regularity - 4
          !pKdiff = pK + 1  - regularity - 2
          !pKdiff = pK + 1  - regularity - 3

          pKdiff = -10.  ! no regularity detection

          ! print *,'###S#D', elem%i, pK, regularity, pKdiff

          if(pKdiff  > 1.) then
             iopt = 0
             !!write( (state%space%adapt%adapt_level+1)*10 + 1, *) elem%xc,  pKdiff, pK, regularity
          elseif(pKdiff  > 0.) then

             if(iopt == 2) then
                iopt = 1
              !!  write( (state%space%adapt%adapt_level+1)*10 + 2, *) elem%xc,  pKdiff, pK, regularity
             endif

          endif

       endif
    endif  !if(state%space%adapt%adapt_space /= 'AMAh' .and. state%space%adapt%adapt_space /= 'IMAh')


    !EDF  iopt = elem%psplit + 1
    if(iprint ) then
       do j=0, 2
          write(*,'(a4,4i5,8es12.4)') &
               'EE:',elem%i, elem%deg, iopt, j,elem%estim_loc / tols,  factorX, factor, Kparams(j, 1:5)
          write(150+state%space%adapt%adapt_level, '(2es12.4, 4i5,8es12.4)') &
               elem%xc,elem%i, elem%deg, iopt, j,elem%estim_loc / tols,  factorX, factor, Kparams(j, 1:5)
       enddo
       print*
       write(150+state%space%adapt%adapt_level, '(x)' )
    endif

    ! setting of the p-adaptation
    elem%psplit = iopt - 1

    if(singularity) then
       deg_opt = 1   ! has to be changed !!
       !deg_opt = 2   ! has to be changed !!

       iopt = max( 0, deg_opt -  elem%deg + 1)   ! choice of the tested metric
       elem%psplit =  deg_opt - elem%deg         ! setting of psplit (chnage of elem%deg)
    !   if(elem%psplit < 0) Kparams(iopt, 1:2) =  Kparams(iopt, 1:2) * 2.

    !   Kparams(iopt, 1:2)  = Kparams(iopt, 1:2) * 4
    !   !write(*,'(a8, 4i5, 30es12.4)' ) 'Sing5s2',elem%i, elem%deg, elem%psplit, iopt, factor
    !else
    !   !write(*,'(a8, 4i5, 30es12.4)' ) 'Reg 5s2',elem%i, elem%deg, elem%psplit, iopt, factor

    !endif

    !if(elem%xc(1) > 0.95 .and. elem%xc(2) > 0.95) then
    !   do j=0, 2
    !      write(*,'(a4,4i5,8es12.4)') &
    !           'DD:',elem%i, elem%deg, iopt, j,elem%estim_loc / tols,  factor, Kparams(j, 1:5)
    !   enddo
    !   print*
    endif


    ! plotting of the metric
    !do ji=0, 3
    do ji =1, -1  ! NO PLOTTING

       j = ji
       if(ji == 3) j = iopt

       lam_max = Kparams(j, 1)
       lam_min = Kparams(j, 2)
       max_a   = Kparams(j, 3)

       elem%rgabc(1) =  lam_max * cos(max_a)**2 + lam_min * sin(max_a)**2
       elem%rgabc(2) = (lam_max - lam_min) * cos(max_a) * sin(max_a)
       elem%rgabc(3) =  lam_max * sin(max_a)**2 + lam_min * cos(max_a)**2

       imt = 10 * (state%space%adapt%adapt_level+1) + ji + 1
       !if( sqrt( (elem%xc(1))**2 + (elem%xc(2))**2 ) < 0.1 )&
       !if( singularity) &
            !if(elem%i >= 4 .and. elem%i <= 5) &
       !write(*,'(a8,3i5, 8es12.4)') 'RGABC =', ji, j, elem%i, elem%rgabc(1:3),Kparams(j,:)
       !if( dot_product(elem%xc, elem%xc) > 1) &
       if( elem%xc(2)  > 0.5) &
            call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )
    enddo

    !if(singularity) &
    !     write(*,'(a10, 2i5,40es12.4)') 'rgabc 5is',elem%i, elem%deg, elem%ama_p, elem%rgabc(1:3)


    ! setting the metric from the best candidate: pK-1, pK, pK+1, averaging with pK

    weight = 1.0
    !if(state%modelName == 'scalar' .or.state%modelName == '2eqs') weight = 0.5

    lam_max = weight * Kparams(iopt, 1) + (1. - weight) * Kparams(1, 1)
    lam_min = weight * Kparams(iopt, 2) + (1. - weight) * Kparams(1, 2)
    max_a   = Kparams(iopt, 3)

    elem%ama_p = weight * elem%psplit

    !if(singularity)  elem%ama_p = 1.5 * elem%psplit

    !if(singularity) &
    !     write(*,'(a10, 2i5,40es12.4)') 'rgabc 5gs',elem%i, elem%deg, elem%ama_p, elem%rgabc(1:3), &
    !     lam_max, lam_min, max_a

    !if(elem%deg + elem%ama_p <  1) then
    !   write(*,'(a20,2es12.4,2i6,4es12.4)') &
    !        '##D#DE#D#DSW', elem%deg + elem%ama_p, elem%ama_p, iopt, elem%deg, elem%xc
    !endif

    ! limitation from the point of view of reasonability, too large and too small elements
    lam_max = max(lambda_min, min(lambda_max, lam_max) )
    lam_min = max(lambda_min, min(lambda_max, lam_min) )

    if(lam_min < lambda_min .or. lam_min > lambda_max .or. &
         lam_min < lambda_min .or. lam_min > lambda_max ) then
       write(*,'(a10, 40es12.4)') 'lams limit',lam_min,lam_max,  lambda_min, lambda_max 
    endif


    !if(singularity) &
    !write(101,'(a10, 2i5,40es12.4)') 'rgabc 5gA',elem%i, elem%deg, elem%ama_p,  &
    !     lam_max, lam_min, max_a

    !lam_actual = 1. / elem%diam**2.0
    lam_actual = (2*elem%area/elem%diam)**(-2)
    lam_actualT = (elem%diam)**(-2)
    ! limitation of the aspect ratio
    ratio = lam_min / lam_max

    if(ratio < 1E-3) then
        write(*,'(a10, 40es12.4)') 'ratio limit', ratio, lam_min, lam_max
     endif

    !if(singularity) &
    !     write(*,'(a10, 2i5,40es12.4)') 'rgabc 5gB',elem%i, elem%deg, elem%ama_p, elem%rgabc(1:3), &
    !     lam_max, lam_min, max_a, lam_actual, lam_actualT


    !!ratio = max(ratio, 1E-4)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    !ratio = max(ratio, 2E-5)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    !ratio = max(ratio, 2.5E-4)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    ratio = max(ratio, 1E-3)  ! maximal aspect ration (=1./ratio) !! (should be positivity)

    !ratio = max(ratio, 1E-2)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    !ratio = max(ratio, 0.25)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    !ratio = max(ratio, 1.0)  ! maximal aspect ration (=1./ratio) !! (should be positivity)

    !ratio = max(ratio, diamT/ diam *0.2)   ! ratio may increase maximally 5-times

    !ratio = max(ratio, 1.)
    if(state%modelName == 'scalar' .or.state%modelName == '2eqs' .or.state%modelName == 'porous') &
         ratio = max(ratio, 1E-3)

    if(iprint ) then
       do j=iopt, iopt
          write(150+state%space%adapt%adapt_level, '(2es12.4, 3i5, a5,8es12.4)') &
               elem%xc,elem%i, elem%deg, iopt, 'K',elem%estim_loc / tols,  factorX, factor, lam_max, lam_min, max_a, -1., lam_actual
       enddo
       print*
       write(150+state%space%adapt%adapt_level, * ) 'ratio = ', ratio
    endif

     ! if( abs(elem%xc(1)) < 0.05 .and. abs(elem%xc(2)) <  0.05) then
     !    write(*,'(a6,2(3es12.4,a2), 3es12.4)') '####',  &
     !         lam_max, lam_min, ratio,'|', &
     !         lam_actual, lam_actualT, lam_actualT / lam_actual, '|',lam_actualT / lam_actual / ratio
     ! endif

    ! limitation of the refinement at one level
    !lam_max_actual_ratio = 10  !4. !10
    lam_max_actual_ratio = 4. ! =10  => not so nice monotone convergence, final grids are ROUGH !!!
    !lam_max_actual_ratio = 6. ! =10  => not so nice monotone convergence, final grids are ROUGH !!!
    !if(.not. state%time_dependent) lam_max_actual_ratio = 50.

    lam_max = min(lam_max, lam_actual * lam_max_actual_ratio)  ! at most r^{1/2} times smaller element
    lam_max = max(lam_max, lam_actual / lam_max_actual_ratio)  ! at most r^{1/2} times bigger element

    if(iprint ) then
       do j=iopt, iopt
          write(150+state%space%adapt%adapt_level, '(2es12.4, 3i5, a5,8es12.4)') &
               elem%xc,elem%i, elem%deg, iopt, 'X',elem%estim_loc / tols,  factorX, factor, lam_max, lam_min, max_a, -1., lam_actual
       enddo
       print*
       write(150+state%space%adapt%adapt_level, * ) 'ratio = ', ratio
    endif


    !if(singularity) &
    !     write(*,'(a10, 2i5,40es12.4)') 'rgabc 5gX',elem%i, elem%deg, elem%ama_p, elem%rgabc(1:3), &
    !     lam_max, lam_min, max_a

    !if(elem%i <=10) write(*,'(a6, 6es12.4)') 'lim 2:', lam_max,  lam_min, 1./sqrt(lam_max),1./sqrt(lam_min)

    ! limitation of the ratio
    !ratio = max(ratio,  lam_actualT / lam_actual / 2)


    ! new lam_min after limitation of the aspect ratio and lam_max
    lam_min = lam_max * ratio

    if(iprint ) then
       do j=iopt, iopt
          write(*,'(a4,4i5,8es12.4)') &
               'GG:',elem%i, elem%deg, iopt, j,elem%estim_loc / tols,  factorX, factor, Kparams(j, 1:5)
          write(150+state%space%adapt%adapt_level, '(2es12.4, 4i5,8es12.4)') &
               elem%xc,elem%i, elem%deg, iopt, j,elem%estim_loc / tols,  factorX, factor, lam_max, lam_min, max_a, -1., ratio
       enddo
       print*
       write(150+state%space%adapt%adapt_level, '(x)')
    endif

    !if(singularity) &
    !     write(*,'(a10, 2i5,40es12.4)') 'rgabc 5gY',elem%i, elem%deg, elem%ama_p, elem%rgabc(1:3), &
    !     lam_max, lam_min, max_a

    if( state%space%adapt%adapt_space == 'IMAh' .or. state%space%adapt%adapt_space == 'IMAhp') then
       lam_min = lam_max
       max_a = 0.
    endif


    !if(elem%i <= 3) print*,'A priori refinement j03id032dwoqw'
    !lam_max = 25 * 2**(state%space%adapt%adapt_level+1)
    !lam_min = lam_max ! / 200
    !max_a = 0.   !pi/2  !0.
    !elem%ama_p = 0

    if(iprint ) then
       do j=iopt, iopt
          write(*,'(a4,4i5,8es12.4)') &
               'II:',elem%i, elem%deg, iopt, j,elem%estim_loc / tols,  factorX, factor, Kparams(j, 1:5)
          write(150+state%space%adapt%adapt_level, '(2es12.4, 4i5,8es12.4)') &
               elem%xc,elem%i, elem%deg, iopt, j,elem%estim_loc / tols,  factorX, factor, lam_max, lam_min, max_a
       enddo
       print*
       write(150+state%space%adapt%adapt_level, *)'------------------------------'
    endif

    ! new metric
    elem%rgabc(1) =  lam_max * cos(max_a)**2 + lam_min * sin(max_a)**2
    elem%rgabc(2) = (lam_max - lam_min) * cos(max_a) * sin(max_a)
    elem%rgabc(3) =  lam_max * sin(max_a)**2 + lam_min * cos(max_a)**2

    elem%rgabc(1:3) = elem%rgabc(1:3)

    !if(singularity) &
    !     write(*,'(a10, 2i5,40es12.4)') 'rgabc 5gZ',elem%i, elem%deg, elem%ama_p, elem%rgabc(1:3), &
    !     lam_max, lam_min, max_a

    ! if(abs(elem%xc(1) - 1.) < 0.05) then
    !    write(*,'(a6,i5,3es14.6, a2, 4es14.6)') &
    !         '####',elem%i,lam_max**(-0.5), lam_min**(-0.5),  lam_min**(-0.5)/  lam_max**(-0.5), &
    !         '|',h_max, h_min
    !    !'####',elem%i,lam_max, lam_min, max_a, '|',elem%rgabc(:),  elem%ama_p
    ! endif

    !ip = elem%deg + elem%psplit
    !write(3000+state%space%adapt%adapt_level*10 + ip, *) elem%xc(:), ip

    area_new = 1./sqrt(elem%rgabc(1) * elem%rgabc(3))
    !if(factorX > 2.) then
    !   write(*,'(a8,3i5, 30es12.4)') 'New_Area:',elem%i, elem%deg, int(elem%ama_p), &
    !        elem%area,  elem%area/factorx, area_new,  area_new/ (elem%area/factorx), factorX, factor
    !endif
    if(factorX > 2. .and. area_new/ (elem%area/factorx) > 1.5) then
    else
       ama_err_total_reduced = ama_err_total_reduced + elem%estim_loc**2
    endif

    !if((elem%i == grid%nelem) ) & !.and.  target_tol < state%space%adapt%tol_max) &
    !     write(*,'(a12,4es12.4,a2,8es12.4)') &
    !     'AMA_pars ii:', ama_err_min, ama_err_aver, ama_err_max,ama_err_total,'|', &
    !     sqrt(ama_err_total_reduced)
    
    deallocate(Kparams)

  end subroutine Set_hp_metric_hpREZ



  !> seeks the best anisotropy frmo the candiates elem%deg + 0,1,2
  subroutine SetOptimalAnisotropyEIp( elem, ndimL, iopt, Kparams )
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    integer, intent(inout) :: iopt
    real, dimension(0:3,1:5), intent(inout) :: Kparams
    integer :: ideg_min, ideg_max, ideg, degP
    real :: ropt
    logical :: singularity

    ropt = 1E+50

    ! testing of the directional derivatives of degree elem%deg, elem%deg+1, elem%deg+2
    if(state%space%adapt%adapt_space == 'AMAh' .or. state%space%adapt%adapt_space == 'IMAh')  then
       ideg_min = 1
       ideg_max = 1
    else  !!! (state%space%adapt%adapt_space == 'AMAhp')
       ideg_min = 0
       ideg_max = 2
    endif

    if(elem%deg == 0 .and. ideg_min == 0) ideg_min = 1  ! we avoid the pathological case

    do ideg = ideg_min, ideg_max   ! hp-AMA
    !do ideg = 1, 1   ! only h-AMA
       degP = elem%deg + ideg   ! given degree of the directional derivative

       !if(degP == 1) degP = 2  ! for P_0 approximation ###P_0

       if(degP > MaxDegreeImplemented -1) then
          Kparams(2, 1:5) = Kparams(1, 1:5)

       else
          if(state%space%adapt%Lq > -0.01) then ! L^q-nrom

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
          Kparams(ideg, 5) = DOFtriang(degP)  * sqrt( Kparams(ideg, 1) * Kparams(ideg, 2))
          !Kparams(ideg, 5) = DOFtriang(degP-1)*( Kparams(ideg, 1) * Kparams(ideg, 2))**0.5

          if(Kparams(ideg, 5) < ropt .and. degP >= 1) then ! for degP=1 only first order derivatives
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

    if(elem%i == -995) &
         write(*, '(a6,2i5,8es12.4)') 'KKpae:',elem%i,  iopt, Kparams(1, 1:5)
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


  end subroutine SetOptimalAnisotropyEIp


  !> seeks the best anisotropy frmo the candiates elem%deg + 0,1,2
  subroutine SetOptimalAnisotropy_AreaFixed( elem, ndimL, iopt, area, Kparams)
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    integer, intent(inout) :: iopt
    real, intent(in) :: area  !   the prescibed element area
    real, dimension(0:3,1:5), intent(inout) :: Kparams
    integer :: ideg_min, ideg_max, ideg, degP
    real :: ropt, areaR
    logical :: singularity

    ropt = 1E+50

    ! testing of the directional derivatives of degree elem%deg, elem%deg+1, elem%deg+2
    if(state%space%adapt%adapt_space == 'AMAh' .or. state%space%adapt%adapt_space == 'IMAh')  then
       ideg_min = 1
       ideg_max = 1
    else  !!! (state%space%adapt%adapt_space == 'AMAhp')
       ideg_min = 0
       ideg_max = 2
    endif

    if(elem%deg == 0 .and. ideg_min == 0) ideg_min = 1  ! we avoid the pathological case

    ! pathological case for errors in H^1 seminorm
    if( abs( state%space%adapt%Lq + 1.0) < 1E-6 .and. elem%deg == 1 .and. ideg_min == 0) ideg_min = 1
    ! ^^^^   LEPSI udelat ve VZTAHU k MinDegreeImplemented !!!

    do ideg = ideg_min, ideg_max   ! hp-AMA
       degP = elem%deg + ideg   ! given degree of the directional derivative

       ! modification of the area in order to keep the same DOF density
       areaR = area * (degP *(degP+1))/ ( (elem%deg+1) *(elem%deg +2))

       if(degP > MaxDegreeImplemented -1) then
          Kparams(2, 1:5) = Kparams(1, 1:5)

       else

          call FindAnisotropy_AreaFixed(elem, ndimL, degP, elem%wSS(1:ndimL, ideg, 0:degP), &
               areaR,  Kparams(ideg,1:5) )

          !write(10+10*state%space%adapt%adapt_level + ideg,'(300es12.4)') &
          !     elem%xc, area, areaR, elem%area 

          ! Kparams(ideg, 5) contains the estimate of the interpolation error for the given area
          if(Kparams(ideg, 5) < ropt .and. degP >= 1) then ! for degP=1 only first order derivatives
             ropt = Kparams(ideg, 5)
             iopt = ideg
          endif

          !if(elem%xc(1) > 0.75 .and. elem%xc(2) > 0.75) &
          !   write(*,'(a10, 4i5, 40es12.4)') 'elem A', elem%i, ideg, degP,iopt, ropt, Kparams(ideg,1:5),&
          !   1./sqrt(Kparams(ideg,1) * Kparams(ideg,2) )

       endif
    enddo  ! ideg


    ! already the MaxDegreeImplemented
    if(elem%deg + elem%psplit >= MaxDegreeImplemented -1 .and. iopt == 2) iopt = 1

    ! minimal degree = 1
    if(elem%deg == 1 .and. iopt == 0) iopt = 1


    ! a priori known singularity
    !call Detect_apriori_known_singularity(elem, singularity)
    !if(singularity) iopt = 0

    !!if(singularity) print*,'Singul:', elem%xc(:), iopt

    !if(VectorNorm(elem%xc(1:2) )  < 0.1) &
    !write(*,*) 'elem', elem%i, iopt, ropt

    ! if the proposed "quality" of approximation is very close to the "quality" of actual degree
    !if( iopt /= 1 .and. Kparams(iopt, 5) >= 0.95 * Kparams(1, 5)) then
    !if( iopt /= 1 .and. Kparams(iopt, 5) >= 0.70 * Kparams(1, 5)) then
       !write(*,'(a8,2i5,3es12.4)') &
       !  'opt',elem%deg, elem%deg + iopt -1, Kparams(iopt, 5),  Kparams(1, 5) , Kparams(iopt, 5)/  Kparams(1, 5)
       !iopt = 1
    !endif

  end subroutine SetOptimalAnisotropy_AreaFixed

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


  !> limitation of the anisotropy due to smoothing
  subroutine  Anisotropy_limitation ( lam_max,  lam_min)
    real, intent(inout) ::  lam_max,  lam_min
    real :: maximal_lambda
    real :: maximal_sigma
    real :: sigma

    maximal_sigma = 5E+03
    maximal_lambda = (4 / state%space%diam  * 5.E+05)**2

    sigma =  lam_max / lam_min

    if(sigma <=  maximal_sigma * pi / 2) then
       sigma =  maximal_sigma * sin( sigma /  maximal_sigma)
    else
       sigma =  maximal_sigma 
    endif

    if(lam_max <=  maximal_lambda * pi / 2) then
       lam_max = maximal_lambda * sin(lam_max /  maximal_lambda)
    else
       lam_max = maximal_lambda
    endif

    lam_min = lam_max / sigma

  end subroutine Anisotropy_limitation


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
    integer :: n, k, j, i, ip, ifig, Qdof, du_max, itest, astop
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

    if(( state%modelName == 'scalar' .or.state%modelName == '2eqs' .or.state%modelName == 'porous' ) &
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

       max_f = -10.
       ! we seek the direct with the maximal directional derivative
       do i=1, n, 2  ! equidistance 2 degrees is enough
          a =  pi * i / n

          f = DirectionalDerivative(degP, der(k, 0:degP), cos(a), sin(a))
          !beta = ( tol_loc / f )**(1./ degP)

          if( f  > max_f )then
             max_f = f
             max_a = a
             !beta_max = beta
          endif

          !astop = 0
          !if( abs( elem%xc(1) - 0.63) < 0.1 .and.  abs( elem%xc(2) - 0.1) < 0.1 &
          !     .and. state%space%adapt%adapt_level >= 3 ) then
          !if(elem%i == 995) then
          !   write(99,'(20es16.8)') elem%xc(:),a,f,max_f, max_a
          !   write(99,'(20es16.8)') elem%xc(1) +  cos(a) * f, elem%xc(2) +  sin(a) * f
          !   !astop = 1
          !endif

       enddo

       !if(elem%i == 995) stop'dedyu83dy832djiws'

       ! size in the perpendicular direction
       a = max_a + pi /2
       min_f = DirectionalDerivative(degP, der(k, 0:degP), cos(a), sin(a) )

       if(min_f > 0.) then
          rho = max_f / min_f  ! initial aspect ratio
       else
          if(max_f <=0.) then
             rho = 1.
          else
             rho = 1E+5
          endif

       endif


       !if(astop  == 1) then
       !   print*,'edte54de3e38d43h48yd3ije382dy3xns4xn4p5x',max_f/min_f
       !   stop
       !endif


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

       ! setting of the optimal area of the element
       if(max_f <= 0.) then
          area = state%space%domain_volume
       else
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

             if(elem%i == -599) write(*,'(a6,6es14.6)') '@@@@@',area1, area, elem%area
          else
             area = (tol_loc * rho**0.5 /  max_f)**(2./degP)
          endif
       endif

       !if(rho > 1E+3) &
       ! write(*,'(a5,2i5,8es12.4)')  '!!!', elem%i, degP, max_f, min_f, rho, max_a
       !     tol_loc, rho**0.5 , c_p_q (degP, Lq),  max_f, (2./degP)

       !write(*,'(a5,2i5,8es12.4)')  '!!!', elem%i, degP, rho, area
       ! ????
       !beta_min = (rho**(1./degP) * area/pi)**0.5
       beta_min = (rho**(1./degP) * area/1.)**0.5
       beta_max = beta_min / rho**(1./degP)

       ! new tests of limitation
       r_lim = 1.
       if( state%modelName == 'NSe' .or. state%modelName == 'pedes' ) then

          if( state%modelName == 'NSe' ) then
             !h_min = 5E-4
             !h_min = 1E-3
             h_min = state%model%Re1 * degP
             !h_min = 5E-4 * degP  ! should be the same as ^^^ (state%model%Re1 * degP )

             !if(state%type_IC == 8 .or. state%model%Re1 == 0) h_min = 0.005
             if(state%type_IC == 8 .or. state%model%Re1 == 0) h_min = 2E-03

             if(state%type_IC /= 8 .and. state%model%Re1 == 0) h_min = 1E-04
             !if(beta_max < 5E-3) print*,'####', beta_max , h_min

          elseif( state%modelName == 'pedes' ) then
             !h_min = 0.02
             !h_min = state%space%diam / 400
             h_min = state%space%diam / 500

          endif

          if( beta_max < h_min) then
             !write(*,'(a10,4i5,8es12.4)')'#!!E#E##',elem%i,degP,k, state%num_limits, beta_min,beta_max,&
             !     beta_min/ beta_max, max_f, h_min, state%model%Re1
             !write(97, *) elem%xc(:)
             beta_max = h_min
             beta_min =   beta_max * rho**(1./degP)

             r_lim = 0.
             if( state%modelName == 'pedes' )  r_lim = 0.1
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
    enddo  ! end of k=1, ndimL


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
    if(degP == elem%deg + 1 ) then

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

          !if(abs(elem%xc(1) - 20) < 1. ) &
               !if(elem%interLq > 0) &
          !     write(*,'(a8, i5, 30es12.4)') '4d43de3',degP, elem%xc(:), wi(1:3), elem%interLq, r_lim

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


  !> seek the optimal anisotropy of the interpolation error function
  !> \f$ E^I_p \f$, for the FIX AREA i.e.
  !> the direction having the maximal directional derivative of order degP
  !> direction is max_a, size of the derivative lam_max, size of the derivative in
  !> the perpendicular direction is lam_min
  !> paper anisot_hp
  subroutine FindAnisotropy_AreaFixed( elem, ndimL, degP, der, area, Kparams)
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    integer, intent(in) :: degP  ! degree of the directional derivative
    real, intent(inout) :: area  !   the prescibed element area
    real, dimension(1:ndimL,  0:degP), intent(inout) :: der  ! partial derivatives of degree degP
    real, dimension(1:5), intent(inout) :: Kparams !  lam_max, lam_min,  max_a, max_f, err
    real, dimension(1:2) :: xi
    real, dimension(:,:), allocatable :: Kpar_loc
    real, dimension(:,:), allocatable :: Fx
    real, dimension(:), allocatable :: wi
    real, dimension(:,:), allocatable :: derPP

    real :: pi, rho, lam_max, lam_min,  max_a, max_f, min_f, r_lim
    real :: a, beta_max, beta_min, err_estim, power
    real :: f, fac, h_min
    real :: Lq, Lqq, c_p, area_const, rho_max
    integer :: deg0, degPP
    integer :: i, k, n, Qdof, ifig, itest, itype
    logical :: singularity

    itest = 25

    ! seeking the element size, orientation and the aspect ratio
    n = 180  ! 360 ! number of test directions
    pi = 2 * asin(1.)

    ifig = 10

    ! L^q, q\in[1,\infty]- norms
    if(state%space%adapt%Lq > -0.01) then
       degPP = degP

       allocate(derPP (1:ndimL, 0:degPP) )
       derPP (1:ndimL, 0:degPP) = der (1:ndimL, 0:degPP)

       itype = 0

    elseif( abs( state%space%adapt%Lq + 1.0) < 1E-6) then
       ! H^1-semi-norm
       ! print*,'ede',degP
       ! area = 0.01
       ! der(1, :) = 0.
       ! der(1, 0) = 6.
       ! der(1, 1) = 0.
       ! der(1, 2) = 0.
       ! der(1, 3) = 0.

       deg0  = degP - 1
       degPP = 2*deg0
       allocate(derPP(1:ndimL, 0:degPP) )

       call Set_Der_coeffs_H1(ndimL, degP, degPP, der (1:ndimL, 0:degP), derPP(1:ndimL, 0:degPP))

       !print*,'--------------------------------------------------'
       !write(*,'(a8, i5, 30es12.4)') 'der P:',degP, der (1:ndimL, 0:degP)
       !write(*,'(a8, i5, 30es12.4)') 'der PP:',degPP, derPP (1:ndimL, 0:degPP)

       itype = 1

    else
       print*, "UNKNOWN (semi-)norm in  FindAnisotropy_AreaFixed", state%space%adapt%Lq
    endif


    allocate(Kpar_loc(1:ndimL, 1:5) )

    fac = factorial(degP)

    ! we go over all quantities taken into account
    do k=1, ndimL

       max_f = -10.
       ! we seek the direct with the maximal directional derivative
       do i=0, n, 2  ! equidistance 2 degrees is enough
          a =  pi * i / n   ! angle


          if(itype == 0) then
             f = DirectionalDerivative(degPP, derPP(k, 0:degPP), cos(a), sin(a))

             !if(elem%xc(1) > 0.95 .and. elem%xc(2) > 0.95) then
             !   write(*,'(a8, 30es12.4)') 'DF max_F f:',f, derPP(k, 0:degPP)
             !endif

          elseif(itype == 1) then
             f = DirectionalDerivative_H1(degPP, derPP(k, 0:degPP), cos(a), sin(a))

             ! if(elem%xc(1) > 0.95 .and. elem%xc(2) > 0.95) then
             !    write(*,'(a8, 30es12.4)') 'DF max_F f:',f, derPP(k, 0:degPP)
             ! endif

          else
             stop "yd83yhd3iu7ju;7[]"
          endif

          if( f  > max_f )then
             max_f = f
             max_a = a
          endif

       enddo

       ! size in the perpendicular direction
       a = max_a + pi /2
       if(itype == 0) then
          min_f = DirectionalDerivative(degPP, derPP(k, 0:degPP), cos(a), sin(a) )

       elseif(itype == 1) then
          min_f = DirectionalDerivative_H1(degPP, derPP(k, 0:degPP), cos(a), sin(a) )

       else
          stop "yd83yhd33w32"
       endif

       if(state%space%adapt%adapt_type == 'Ihp' ) then
          min_f = max_f
          max_a = 0.
       endif

       if(min_f > 0.) then
          rho = max_f / min_f  ! aspect ratio
       else
          if(max_f <=0.) then  ! derivatives are very very small, hence isotropic metric
             rho = 1.
          else
             rho = 1E+15
          endif

       endif

       !rho = min(rho, 1E+15)  ! in order to avoid vanishing derivatives in one direction
       
       ! limitation of the aspect ratio
       !rho_max = 1E+3

       rho_max = 20.
       if(itype == 0) rho_max = rho_max**(degP)        ! maximal ratio h_max/h_min
       if(itype == 1) rho_max = rho_max**(2*deg0) 

       rho = min(rho, rho_max)  

       !area_const = 1.  ! the simplicity
       !area_const = pi  ! ellipse
       area_const = 3.* sqrt(3.) / 4  ! isoscele triangle

       ! setting of the element sizes in both directions
       if(itype == 0) then
          beta_min = sqrt(rho**(1./degP) * area/ area_const )
          beta_max = beta_min / rho**(1./degP)

       elseif(itype == 1) then
          beta_min = sqrt(rho**(1./(2*deg0)) * area /area_const  )
          beta_max = beta_min / rho**(1./(2*deg0))
       else
          stop "yd83yhd3"
       endif

       !if(elem%xc(1) > 0.75 .and. elem%xc(2) > 0.75) then
       !    write(*,'(a15, 3i5, 30es12.4)') &
       !         'ify  !!', elem%i, elem%deg, degP,  max_f, rho, &
       !         beta_min, beta_max, beta_min / beta_max, area, elem%area
       ! endif


       ! evaluation of the estimate of the interpolation error
       Lq = state%space%adapt%Lq
       if(itype == 0) then
          if(Lq >=1.) then    ! q < \infty
             power = degP / 2. + 1./Lq
             err_estim =  c_p_q (degP, Lq) * max_f / sqrt(rho) * area**power

          else !q = infty
             power = degP / 2.
             err_estim =   max_f / sqrt(rho) * area**power
          endif

       elseif(itype == 1) then ! H^1-seminorm
          c_p = pi**(-degP+1)/ degP
          power = degP
          err_estim =  sqrt(c_p * max_f / sqrt(rho) * area**power )



       else
          stop "yd83yhd38jo2;11"

       endif


       ! detection of the singularity
       !singularity = .false.
       !call Detect_apriori_known_singularity(elem, singularity)

       ! if(elem%i < 5 .or. singularity) then
       !write(*,'(a25, 2i5, 30es12.4)') &
       !     'Verify 9eu39jud93oi !!!!!', elem%i,degP, err_estim, area, max_f, min_f, &
       !     max_f / sqrt(rho)
       !      if(degP == elem%deg + 2) print*
       !   endif

       !call PlotElemFunction3D(71, elem,  elem%dof, elem%w(0, 1:elem%dof) )

       !if(state%space%adapt%adapt_level == 1 .and. elem%i == 2) &
       !     stop '###S#@SD#LKDS_#@LKPLDS#@PLDSP@SIS#NH+}D~{§§§'

       goto 99
       ! ! new tests of limitation
       ! r_lim = 1.
       ! if( state%modelName == 'NSe' .or. state%modelName == 'pedes' ) then

       !    if( state%modelName == 'NSe' ) then
       !       !h_min = 5E-4
       !       !h_min = 1E-3
       !       h_min = state%model%Re1 * degP
       !       !h_min = 5E-4 * degP  ! should be the same as ^^^ (state%model%Re1 * degP )

       !       !if(state%type_IC == 8 .or. state%model%Re1 == 0) h_min = 0.005
       !       if(state%type_IC == 8 .or. state%model%Re1 == 0) h_min = 2E-03

       !       if(state%type_IC /= 8 .and. state%model%Re1 == 0) h_min = 1E-04
       !       !if(beta_max < 5E-3) print*,'####', beta_max , h_min

       !    elseif( state%modelName == 'pedes' ) then
       !       !h_min = 0.02
       !       !h_min = state%space%diam / 400
       !       h_min = state%space%diam / 500

       !    endif

       !    if( beta_max < h_min) then
       !       !write(*,'(a10,4i5,8es12.4)')'#!!E#E##',elem%i,degP,k, state%num_limits, beta_min,beta_max,&
       !       !     beta_min/ beta_max, max_f, h_min, state%model%Re1
       !       !write(97, *) elem%xc(:)
       !       beta_max = h_min
       !       beta_min =   beta_max * rho**(1./degP)

       !       r_lim = 0.
       !       if( state%modelName == 'pedes' )  r_lim = 0.1
       !       state%num_limits = state%num_limits + 1

       !    endif
       ! endif
99     continue

       !beta_min = max(beta_min, 1E-3)
       !beta_max = max(beta_max, 1E-3)

       lam_max = 1./beta_max**2
       lam_min = 1./beta_min**2


       !beta_max = lam_max
       !beta_min = lam_min
       !call Anisotropy_limitation ( lam_max,  lam_min)

       !if( abs(beta_max - lam_max) > 1E-1 .or. &
       !     abs(beta_min - lam_min) > 1E-1 .or. &
       !     abs(beta_max/beta_min - lam_max/lam_min) > 1E-1 ) then
       !   write(*,'(a6,5es14.6)') 'lim 0',  beta_max,  beta_min,  beta_max/ beta_min
       !   write(*,'(a6,5es14.6)') 'lim 1',  lam_max,  lam_min,  lam_max/ lam_min
       !   print*
       !endif

       Kpar_loc(k,1) = lam_max
       Kpar_loc(k,2) = lam_min
       Kpar_loc(k,3) = max_a
       Kpar_loc(k,5) = err_estim

       if(elem%i == -1) &
            write(*,'(a6,3i5,5es14.6, a2,20es12.4)') 'WW',elem%i,degP,k,Kpar_loc(k,1:5),'|' !,  der(k, 0:degP)
    enddo  ! end of k=1, ndimL


    ! metric compositions
    if(ndimL > 1) then

       ! SMAZ the following
       !Kparams(1) = lam_max
       !Kparams(2) = lam_min
       !Kparams(3) = max_a
       Kparams(4) = max_f
       Kparams(5) = maxval(Kpar_loc(1:ndimL, 5))

       call MetricComposition(ndimL, Kpar_loc(1:ndimL, 1:3),  Kparams(1:3), elem%xc(1:2), elem%i, degP )

    else
       Kparams(1) = lam_max
       Kparams(2) = lam_min
       Kparams(3) = max_a
       Kparams(4) = max_f
       Kparams(5) = err_estim
    endif

    if(elem%i == -1) then
       write(*,'(a6,3i5,3es14.6, a2,20es12.4)') 'QQ',elem%i,degP,0,Kparams(1:3),'|' !,  der(k, 0:degP)
       print*,'#########################################################'
    endif
    !!if(degP == elem%deg+ 2) stop


    ! SKIPPED
    ! estimates of the interpolation error in the L^{\infty} norm
    if(degP == elem%deg + 1 .and. degP == elem%deg + 1000000 ) then

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
             wi(i) = DirectionalDerivative(degPP, derPP(k, 0:degPP), xi(1),xi(2) )

             ! limitation of the maximal derivative
             wi(i) = wi(i) * r_lim

             ! q-th power
             wi(i) = (abs(wi(i)) )**Lqq
          enddo
          ! error  in the Lq norm
          call IntegrateFunction(elem, wi, elem%interLq)

          !if(abs(elem%xc(1) - 20) < 1. ) &
               !if(elem%interLq > 0) &
          !     write(*,'(a8, i5, 30es12.4)') '4d43de3',degP, elem%xc(:), wi(1:3), elem%interLq, r_lim

          ! NOT POWERED !!!
          !elem%interLq = elem%interLq**(1./Lqq)

          deallocate(Fx, wi)

          ! interpolation error function estimate in L^\infty norm
          do i=1,3  ! vertexes of triangle
             xi(1:2) = grid%x(elem%face(idx, i), 1:2) - elem%xc(1:2)

             f = DirectionalDerivative(degPP, derPP(k, 0:degPP), xi(1),xi(2) )

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
    deallocate(Kpar_loc, derPP)

  end subroutine FindAnisotropy_AreaFixed

  !> setting of coefficients of the derivatives for H^1-seminorm in order to use the
  !> same procedure above
  subroutine Set_Der_coeffs_H1(ndimL, degP, degPP, der, derPP)
    integer, intent(in) :: ndimL, degP, degPP
    real, dimension (1:ndimL, 0:degP) :: der
    real, dimension (1:ndimL, 0:degPP) :: derPP
    real, dimension(:,:), allocatable :: alpha
    integer :: k, l, ll

    allocate(alpha(1:5, 0:degPP) )

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

    deallocate(alpha)

  end subroutine Set_Der_coeffs_H1


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
       !print*,'   ... TROUBLE in FindAnisotropyEIp_H1', degP-1
       return
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

    if(( state%modelName == 'scalar' .or.state%modelName == '2eqs' .or.state%modelName == 'porous' ) &
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

    if(ie == itest)  write(*,'(a6, 6es12.4)' ) 'init', rgabc(1,1:3), xc(1:2)

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


    lwork = 10
    allocate(D(1:2),E(1:1), w(1:2), A(1:2, 1:2), work(1:lwork) )
    D(1) = rgabc(1, 1)
    D(2) = rgabc(1, 3)
    E(1) = rgabc(1, 2)
    A = 0.
    work = 0.

    !write(*,'(a8, 30es12.4)') "M 1:", D(1), E(1)
    !write(*,'(a8, 30es12.4)') "M 2:", E(1), D(2)

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

    !if(degP == 2) print*,'#@#$#!#', degP, der(:)

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
       !if(a1 > 0 .and. a2 >0) &
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



  !> evaluate all partial derivatives of \f$ w \f$ on elem
  subroutine Eval_All_Derivatives(elem, ndimL)   !!!, wp)
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    !real, dimension(1:grid%npoin, 1:ndim), intent(inout) :: wp
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:,:,:), allocatable :: Dw    ! all derivatives of w: bases coefficients
    real, dimension(:,:,:,:), allocatable :: Dwi   ! all derivatives of w in integ node
    real, dimension(:,:,:), allocatable :: Dphi ! derivative of the test functions
    real, dimension(:,:), allocatable :: MassInv ! Inversion of local mass matrix of order deg+1!!
    !real, dimension(:,:), allocatable :: Mass, MassS ! local mass matrix of order deg+1!!
    real, dimension(:), allocatable :: ident, vec ! temporary array
    integer :: Qnum, Qdof, dof, dofP, deg, degP
    integer :: ideg, ideg_min, ideg_max, i, j, k, l,  ix, jx

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


    !if(Qnum > maxVrule ) then
    !   Qnum = maxVrule
    !   Qdiff = Qnum - elem%Qnum
    !   print*,'$$ alert, possible troubles in ama-hp_metric.f90'
    !endif

    Qdof = state%space%V_rule(Qnum)%Qdof
    phi => state%space%V_rule(Qnum)%phi(1:dof, 1:Qdof)

    !!!!!Qdof = elem%Qdof OLD!!

    if(state%space%adapt%adapt_space == 'AMAh' .or. state%space%adapt%adapt_space == 'IMAh')  then
       ideg_min = 1
       ideg_max = 1
    else  !!! (state%space%adapt%adapt_space == 'AMAhp')
       ideg_min = 0
       ideg_max = 2
    endif

    if(deg == 0 .and. ideg_min == 0) ideg_min = 1 ! otherwise a patological case

    do ideg = ideg_min, ideg_max  ! we evaluate the derivatives of order deg+ideg

       degP = deg + ideg
       !###P_0
       !if(degP == 1) degP = 2

       !if(degP >= 8 ) print*,'#####',elem%deg, ideg, degP,  MaxDegreeImplemented

       if(degP == -1) then
          ! direct computation of the second order derivatives
          !?? allocate(Dw (1:ndimL, 0:degP, 0:degP, 1:dofP) )

          call P_0_to_P2_Interpolation(elem, ndimL, .false., Qnum, degP, dofP, &
               Dw (1:ndimL, 0:degP, 0:degP, 1:dofP), -1 )

       else
          if(degP <= MaxDegreeImplemented) then

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
             !call LeastSquareInterpolationH1(elem, ndimL, .false., Qnum, degP, dofP, &
             !     Dw (1:ndimL, 0:degP, 0:degP, 1:dofP), 0 )

             !call LeastSquareInterpolationL2(elem, ndimL, .false., Qnum, degP, dofP, &
             !     Dw (1:ndimL, 0:degP, 0:degP, 1:dofP), 0 )

             ! higher order reconstruction via least square method
             call LeastSquareInterpolationWeighted(elem, ndimL, .false., Qnum, degP, dofP, &
                  Dw (1:ndimL, 0:degP, 0:degP, 1:dofP), -1 )  ! last argument=0 -> L2, = -1 -> H1

             ! components of 1..dofP - degP - 1 has no influence on the highest order derivatives
             Dw (1:ndimL, 0, 0, 1:dofP - degP - 1) = 0.

             if(ideg == 0 ) then
                if(elem%i == 1) print*,'ENFORSING Element SOlUtiON'
                Dw(1, 0, 0, 1:dofP) = elem%w(0, 1:dofP)
             endif


             !if(elem%i == 1) print*,'fort.11 written HEDt53d'
             !if(elem%xc(1) > 19. .and. elem%xc(1) < 21) &
             !     !print*, elem%i , elem%xc
             !     call PlotElemFunction3D(11, elem,  elem%dof, elem%wS(1, 1:elem%dof) )

             !if(elem%i == 1) print*,'fort.21 written HEDt53d'
             !if(elem%xc(1) < 5.) &
             !if(elem%i == -599) &
             !if(elem%xc(1) > 0.95 .and. elem%xc(2) > 0.95) then
             !    call PlotElemFunction3D(500+10*state%space%adapt%adapt_level+ideg, &
             !         elem,  dofP, Dw(1, 0, 0, 1:dofP) )
             !    !!!print*,'DF:',  Dw (1:ndimL, 0:degP, 0:degP, 1:dofP)
             ! endif



             ! do i=0,degP
             !    do j=0,degP
             !       write(*,'(a3,2i5,200es14.6)') '@@@',i,j,Dw(1,i,j,:)
             !    enddo
             ! enddo
             !stop

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

             if(elem%i == -995) then
             !if(elem%xc(1) > 0.95 .and. elem%xc(2) > 0.95) then
                !if( abs(elem%xc(1) -20 ) < 0.5 ) then
                write(*,*) '-----', elem%i, '  ----', elem%xc
                k = 1
                do i=1,degP   ! i = i-th derivative
                   do j=0,i  !   d^i w /dx^j dy^(i-j)
                      write(*,'(a6,3i5,500es12.4)') 'Dwi:',i,j,Qdof,Dwi(k, i, j, 1:Qdof)
                   enddo
                enddo

                write(*,*) '-----------------------------------'

                do i=0,degP   ! i = i-th derivative
                   do j=0,i  !   d^i w /dx^j dy^(i-j)
                      write(*,'(a6,3i5,500es12.4)') 'Dw:',i,j,dofP,Dw(k, i, j, 1:dofP)
                   enddo
                enddo

                write(*,*) 'DF ########################################'
                !stop 'u39j3o'
             endif

             do k=1, ndimL
                do i=0, degP
                   ! elem%wSS contains only the derivatives of maximal order
                   !Dwi(k,degP, i, 1:Qdof) is almost constant
                   elem%wSS(k, ideg, i) = sum(Dwi(k,degP, i, 1:Qdof) )/Qdof / phi(1,1)

                   if( sqrt( (elem%xc(1) - 0.5)**2 + (elem%xc(2) -0.5)**2 ) < -0.05 ) &
                        write(*,'(2es12.4,4i5, 2es12.4,a2, 120es12.4)') &
                        elem%xc(:), elem%i, ideg, degP, i, &
                        Dwi(k,degP, i, 1:2),'|',elem%wSS(k, ideg, i)

                enddo
                ! if( sqrt( (elem%xc(1) - 0.5)**2 + (elem%xc(2) -0.5)**2 ) < 0.05 ) print*
                !write(100 + 10*state%space%adapt%adapt_level + ideg,'(120es12.4)') &
                !if( elem%xc(1) > 19 .and. elem%xc(1) < 21) &
                !if(elem%i == -599) &
                !if(elem%xc(1) > 0.95 .and. elem%xc(2) > 0.95) &
                !     write(*,'(2es12.4,3i5, a8, 120es12.4)') &
                !     elem%xc(:), elem%i, ideg, degP, &
                !     '| nd |',elem%wSS(k, ideg, 0: degP)*phi(1, 1)

             enddo

             deallocate(Dw, Dwi,  vec, Dphi, ident)
             deallocate( MassInv)

          else
             ! this degree is not implemented, hence too high derivative
             elem%wSS(1:ndimL, ideg, 0:degP) = 1.E+60
          endif

          !deallocate( Mass, MassS)

       endif ! if(degP == 1)
    end do  ! ideg

    !if(elem%xc(1) > 0.95 .and. elem%xc(2) > 0.95)  &
    !     print*,' . .     . .'
    !!stop


  end subroutine Eval_All_Derivatives




  !> evaluate all partial derivatives of \f$ w \f$ on elem
  subroutine Eval_All_Derivatives_by_Projection(elem, ndimL)   !!!, wp)
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    type(volume_rule), pointer :: V_rule
     !real, dimension(1:grid%npoin, 1:ndim), intent(inout) :: wp
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:,:,:), allocatable :: Dw    ! all derivatives of w: bases coefficients
    real, dimension(:,:,:,:), allocatable :: Dwi   ! all derivatives of w in integ node
    real, dimension(:,:,:), allocatable :: Dww   ! all derivatives of w: bases coefficients
    real, dimension(:,:), allocatable :: psi ! canonical basis functions
    real, dimension(:,:), allocatable :: Fx  ! physical coordinates
    real, dimension(:), allocatable :: weights ! Inversion of local mass matrix of order deg+1!!
    real, dimension(:,:), allocatable :: Mass ! Local mass matrix of order deg+1!!
    real, dimension(:,:), allocatable :: MassInv ! Inversion of local mass matrix of order deg+1!!
    !real, dimension(:,:), allocatable :: Mass, MassS ! local mass matrix of order deg+1!!
    real, dimension(:,:,:), allocatable ::  vec  ! temporary array
    real, dimension(:), allocatable ::    wi ! temporary array
    real :: val, norm
    integer :: Qnum, Qdof, dof, dofP, deg, degP, degS
    integer :: ideg, ideg_min, ideg_max, i, j, k, l,  ix, jx

    deg = elem%deg
    dof = elem%dof
    !Qnum = elem%Qnum + Qdiff

    ! increase of quadrature
    if(deg == MaxDegreeImplemented ) then
       print*,'$$ alert, possible troubles in ama-hp_metric.f90 (2)'
       Qnum = elem%Qnum
    else
       Qnum = state%space%Qdeg(deg+1, 1)
    endif


    !if(Qnum > maxVrule ) then
    !   Qnum = maxVrule
    !   Qdiff = Qnum - elem%Qnum
    !   print*,'$$ alert, possible troubles in ama-hp_metric.f90'
    !endif

    V_rule => state%space%V_rule(Qnum)
    Qdof = V_rule%Qdof

    !!!!!Qdof = elem%Qdof OLD!!

    if(state%space%adapt%adapt_space == 'AMAh' .or. state%space%adapt%adapt_space == 'IMAh')  then
       ideg_min = 1
       ideg_max = 1
    else  !!! (state%space%adapt%adapt_space == 'AMAhp')
       ideg_min = 0
       ideg_max = 2
    endif

    if(deg == 0 .and. ideg_min == 0) ideg_min = 1 ! otherwise a patological case

    ! Dw(k,i,j, l) => k-th component of w, (i,j) d^i/(dx^j dy^(i-j)), l coeff or integ node
    dofP = DOFtriang( deg + ideg_max)

    allocate(Dww (ideg_min: ideg_max, 1:ndimL, 1:dofP) )
    Dww = 0.

    do ideg = ideg_min, ideg_max  ! we evaluate the derivatives of order deg+ideg

       degP = deg + ideg
       !###P_0
       !if(degP == 1) degP = 2

       !if(degP >= 8 ) print*,'#####',elem%deg, ideg, degP,  MaxDegreeImplemented

       if(degP == -1) then
          ! direct computation of the second order derivatives
          !?? allocate(Dw (1:ndimL, 0:degP, 0:degP, 1:dofP) )

          call P_0_to_P2_Interpolation(elem, ndimL, .false., Qnum, degP, dofP, &
               Dw (1:ndimL, 0:degP, 0:degP, 1:dofP), -1 )

       else
          if(degP <= MaxDegreeImplemented) then

             dofP = DOFtriang( degP )

             ! Dw(k,i,j, l) => k-th component of w, (i,j) d^i/(dx^j dy^(i-j)), l coeff or integ node
             allocate(Dw (1:ndimL, 0:degP, 0:degP, 1:dofP) )
             Dw = 0.

             ! higher order reconstruction via least square method
             call LeastSquareInterpolationWeighted(elem, ndimL, .false., Qnum, degP, dofP, &
                  Dw (1:ndimL, 0:degP, 0:degP, 1:dofP), -1 )  ! last argument=0 -> L2, = -1 -> H1

             !if(ideg == 1) &
             !     call LeastSquareRecovery_Dw_mean_values(elem, ndimL, .false., Qnum, degP, dofP, &
             !     Dw (1:ndimL, 0:degP, 0:degP, 1:dofP), -1 )


             Dww (ideg, 1:ndimL, 1:dofP) = Dw (1:ndimL, 0, 0, 1:dofP)

             deallocate(Dw)
          else
             ! this degree is not implemented, hence too high derivative
             elem%wSS(1:ndimL, ideg, 0:degP) = 1.E+60
          endif

          !deallocate( Mass, MassS)

       endif ! if(degP == 1)
    end do  ! ideg


    ! computation of the HO derivatives by the projection
    dofP = DOFtriang( deg + ideg_max)   ! the maximal DOF
    allocate( psi(1:dofP, 1:Qdof) )

    allocate(Fx(1:Qdof, 1:2) )
    call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )

    Fx(1:Qdof, 1) = (Fx(1:Qdof, 1) - elem%xc(1)) / elem%diam
    Fx(1:Qdof, 2) = (Fx(1:Qdof, 2) - elem%xc(2)) / elem%diam

    ! canonical basis functions (x- xc)**j (y - yc)**(i-j)
    k = 0
    do i=0, degP
       do j=0, i
          k = k + 1
          if( i + j == 0) then
             psi(k, 1:Qdof) = 1.

          elseif( j == 0) then
             psi(k, 1:Qdof) = Fx(1:Qdof, 1)**(i - j)

          elseif( j - i == 0) then
             psi(k, 1:Qdof) = Fx(1:Qdof, 2)**j

          else
             psi(k, 1:Qdof) = Fx(1:Qdof, 1)**(i - j) * Fx(1:Qdof, 2)**j

          endif

          psi(k, 1:Qdof) =  psi(k, 1:Qdof) !/ elem%diam**i
          !do l=1, Qdof
          !   write(100+k, *) Fx(l, 1:2) + elem%xc(1:2), psi(k,l)
          !enddo

       enddo
    enddo

    ! mass matrix
    allocate(weights(1:Qdof) )
    call Eval_V_Weights_plus(elem, V_rule, weights(1:Qdof) )

    allocate( Mass(1:dofP, 1:dofP) )
    allocate( vec (ideg_min: ideg_max,  1:dofP, 1:ndimL ) )
    allocate( wi (1:Qdof) )

    phi => V_rule%phi(1:dofP, 1:Qdof)

    !print*,'#####', elem%F%JF0

    do i=1,dofP
       do j=1,i
          !if(state%space%adapt%adapt_level >= 7) then
          !   write(*,'(a8,i5, 300es12.4)') '### we', i, weights(1:Qdof)
          !   write(*,'(a8,i5, 300es12.4)') '### ps', j, psi(j, 1:Qdof)
          !endif


          Mass(i,j) = dot_product(weights(1:Qdof) *  psi(i, 1:Qdof), psi(j, 1:Qdof) ) !!*elem%F%JF0
          Mass(j,i) = Mass(i,j)
       enddo

       do k = 1, ndimL

          do ideg=ideg_min, ideg_max

             wi(1:Qdof) = matmul( Dww (ideg, k, 1:dofP), phi(1:dofP, 1:Qdof) )

             vec(ideg, i, k) = dot_product(weights(1:Qdof) *  psi(i, 1:Qdof), wi(1:Qdof) ) !!*elem%F%JF0
          enddo
       enddo
    enddo


    !do l=1, dofP
    !    write(*,'(i5, 30es12.4)') l, Mass(l, 1:dofP)
    ! enddo
    ! print*
    ! do ideg=ideg_min, ideg_max
    !    write(*,'(i5, 30es12.4)') ideg, Dww (ideg, 1, 1:dofP)
    ! enddo
    ! print*
    ! do ideg=ideg_min, ideg_max
    !    write(*,'(i5, 30es12.4)') ideg, vec(ideg, :, 1)
    ! enddo
    ! print*,'------------------------------------------------'

    ! evaluation of the function in canonical basis

    do ideg = ideg_min, ideg_max
       degP = deg + ideg
       if(degP <= MaxDegreeImplemented) then
          dofP = DOFtriang( degP )
          allocate( MassInv(1:dofP, 1:dofP) )

          MassInv(1:dofP, 1:dofP) = Mass(1:dofP, 1:dofP)

          call SolveLocalMatrixProblem(dofP, MassInv(1:dofP, 1:dofP), ndimL, vec(ideg, 1:dofP,1:ndimL))

          !do k=1,ndimL
          !   write(*,'(a4, i5, 30es12.4)') 'sols:', ideg, vec(ideg, 1:dofP, 1)
          !enddo


          degS = dofP - degP
          norm =  1./ (elem%diam)**degP
          ! settingf of the HO derivatives
          do k=1,ndimL
             do i=0, degP   ! D^{degP}{ Dx^{degP-i}}
                val = vec(ideg, degS + i, k) * norm * factorial(degP-i) * factorial(i)

                elem%wSS(k, ideg, i) = val / phi(1, 1)  ! constant value of the derivative

                !write(*,'(a6,3i5,500es12.4)') 'Dw:',degP,dofP, i, val, elem%wSS(k, ideg, i),norm,  &
                !     factorial(degP-i), factorial(i)
             enddo

             !if(elem%xc(1) > 0.95 .and. elem%xc(2) > 0.95) &
             !     write(*,'(2es12.4,3i5, a8, 120es12.4)') &
             !     elem%xc(:), elem%i, ideg, degP, &
             !     '| Md |',elem%wSS(k, ideg, 0: degP)*phi(1, 1)
          enddo

          deallocate(MassInv)

       endif
    enddo  ! do ideg =

    !if(elem%xc(1) > 0.95 .and. elem%xc(2) > 0.95) &
    !     print*,'-----------------------------------'

    !stop '73d93djjsdyu923w'

    deallocate(Mass, vec, wi)

  end subroutine Eval_All_Derivatives_by_Projection

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
    real :: weight
    integer :: i, j, k, ipoc
    logical :: smooth_max, smooth_max1

    smooth_max = .false.
    smooth_max1 = .false.

    !ipoc = 1;   weight = 0.5     ! ESCO scalar case
    ipoc = 2;   weight = 0.75     ! ESCO NS case ??
    !ipoc = 3;   weight = 0.4

    if(state%modelName == 'NSe'.and. state%space%estim_space == 'inter' ) then
       !!ipoc = 3;   weight = 0.75  ! ESCO NS case ??
       !ipoc = 2;   weight = 0.25  ! ESCO NS case ??
       ipoc = 1;   weight = 0.5     ! ANGENER with almost non_obtuse elements
    endif

    if(state%modelName == 'pedes'.and. state%space%estim_space == 'inter' )  then
       ipoc = 1;   weight = 0.5     ! ANGENER with almost non_obtuse elements
    endif

    !!state%space%adapt%adapt_method == 'ANI') &
    if( state%space%estim_space == 'RES' .or. state%space%estim_space == 'DWR' .or.  &
         state%space%estim_space == 'pNeu' ) then
       !!ipoc = 3;   weight = 0.25
       ipoc = 3;   weight = 0.4
       !ipoc = 2;   weight = 1.0
       !ipoc = 1;   weight = 1.0
       !ipoc = 1;   weight = 0.5
       !smooth_max = .true.

       if(state%modelName == 'NSe') then
          ipoc = 2;   weight = 1.0
       endif

    endif

   if(state%time_dependent) then
      ipoc = 1;   weight = 0.85
   endif


    if(state%space%adapt%adapt_type == 'Ihp' ) then
       !ipoc = 1;   weight = 1.0 ! case K
       !ipoc = 2;   weight = 0.8 ! case K
       !ipoc = 1;   weight = 0.4  !case E, J   ! make no sense for smooth_max = .true.
       ipoc = 2;   weight = 1.0  ! battery, RES
       smooth_max = .true.
       !smooth_max = .true.
    endif


    !weight = 0.5
    !weight = 0.25

    print*, "Smoothing parameters: ipoc,  weight ", ipoc,  weight

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

                   rgabc(i,1:3) = rgabc(i,1:3) + elem1%rgabc(1:3) * weight
                   rgabc(i,4) = rgabc(i,4) + weight
                endif
             enddo

          endif
       enddo

       do i=1,grid%nelem
          elem => grid%elem(i)

          rgabc(i,1:3) = rgabc(i,1:3)/ rgabc(i, 4 )

          if(smooth_max) then  ! we do not change the strongest metric characterized by lam_1 * lam_2

             if(elem%rgabc(1) * elem%rgabc(3) > rgabc(i,1)*rgabc(i, 3) ) then
                !elem%rgabc(1:3) = elem%rgabc(1:3)
             else
                elem%rgabc(1:3) = rgabc(i,1:3)
             endif

          else
             elem%rgabc(1:3) = rgabc(i,1:3)
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
  !> needed precomputed input arguments - elem%estim_loc, elem%reg, elem%regT0, regT2, elem%deg
  !> parameters Algol1,Algol2,Algol3, hp_adapt
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

     !Algol2 is default
    !Algol3 = .true.
    Algol3 = .false.

    print *,'!! Subroutine IsotropicMetric, Algol1, Algol3 :', Algol1, Algol3

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
             loc_tol = scale *  state%space%adapt%tol_min &
                  *(( state%time%ttime - state%time%ttime_save )/ state%time%FinTime)**0.5
          endif
       endif


       ! we decide p-refinement, p-kept, p-derefinement,  one-parameter variant
       elem%psplit = 0
       ratio = 1.

       if(.not. hp_adapt .and. state%time_dependent) then
          ratio = (elem%estim_loc/ loc_tol )**(1./elem%deg)

       else
          !write(*,'(a6,i5,6es12.4)') 'EDWS',elem%i, &
          !     elem%estim_loc, loc_tol, elem%reg, elem%regT0 , elem%regT2

          !sum = sum +
          !write(50+state%space%adapt%adapt_level) sum

          if(elem%estim_loc > loc_tol  ) then

             if(elem%reg < elem%regT0 .and. elem%deg < MaxDegreeImplemented - 2) then !p-refinement
                elem%psplit = 1
                iA = iA + 1
                !ratio = 0.75
                !ratio = (elem%estim_loc/ loc_tol )**(1./(elem%deg+1) )

                !write(100 + 10*state%space%adapt%adapt_level + 1, '(3es14.6,i5,10es14.6)') &
                !     elem%xc(:),ratio,elem%psplit, &
                !     elem%estim_loc/ loc_tol, (elem%estim_loc/ loc_tol )**(1./elem%deg), &
                !     elem%reg / elem%regT0 , elem%reg / elem%regT2

             else
                ! h-refinement
                ratio = (elem%estim_loc/ loc_tol )**(1./elem%deg)
                !ratio = ratio
                !ratio = ratio * 1.1  ! 1.25  !1.25  !case E, J, K

                !WSW!
                if(Algol1) ratio = 2.
                if(Algol3) ratio = (elem%estim_loc/ loc_tol )**(1./elem%deg)

                iB = iB + 1
                !write(100 + 10*state%space%adapt%adapt_level + 2,  '(3es14.6,i5,10es14.6)') &
                !write(*,'(a6,6es12.4)') 'h_ref:',&
                !     elem%xc(:),sqrt(dot_product(elem%xc, elem%xc) ), elem%diam, ratio
                !!     elem%estim_loc/ loc_tol, (elem%estim_loc/ loc_tol )**(1./elem%deg), &
                !!     elem%reg / elem%regT0 , elem%reg / elem%regT2

             endif

          elseif(elem%estim_loc <= deref_level * loc_tol  ) then

             !if(Algol3) deref_level = deref_level / 2.

             if( elem%reg < elem%regT0 .and. elem%deg < MaxDegreeImplemented - 2) then
                elem%psplit = 1    !p-refinement
                ratio = (elem%estim_loc/ (deref_level *loc_tol) )**(1./(elem%deg+1) )

                iC = iC + 1
                !write(100 + 10*state%space%adapt%adapt_level + 3,  '(3es14.6,i5,10es14.6)') &
                !     elem%xc(:),ratio,elem%psplit, &
                !     elem%estim_loc/ loc_tol, (elem%estim_loc/ loc_tol )**(1./elem%deg), &
                !     elem%reg / elem%regT0 , elem%reg / elem%regT2


             !elseif( elem%reg > elem%regT2   .and.  elem%deg > 1) then
             !   ! p-derefinement
             !   elem%psplit = -1
             !   ratio = (elem%estim_loc/ (deref_level * loc_tol) )**(1./(elem%deg-1))
             !
             !   iD = iD + 1
             !   !write(100 + 10*state%space%adapt%adapt_level + 4,  '(3es14.6,i5,10es14.6)') &
             !   !     elem%xc(:),ratio,elem%psplit, &
             !   !     elem%estim_loc/ loc_tol, (elem%estim_loc/ loc_tol )**(1./elem%deg), &
             !   !     elem%reg / elem%regT0 , elem%reg / elem%regT2
             !

             else
                 ! h-(de)refinement
                ratio = (elem%estim_loc/ (deref_level * loc_tol) )**(1./elem%deg)
                iE = iE + 1
                !write(100 + 10*state%space%adapt%adapt_level + 5,  '(3es14.6,i5,10es14.6)') &
                !     elem%xc(:),ratio,elem%psplit, &
                !     elem%estim_loc/ loc_tol, (elem%estim_loc/ loc_tol )**(1./elem%deg), &
                !     elem%reg / elem%regT0 , elem%reg / elem%regT2

             endif

          else
             iFF = iFF + 1
             !write(100 + 10*state%space%adapt%adapt_level + 6,  '(3es14.6,i5,10es14.6)') &
             !     elem%xc(:),ratio,elem%psplit, &
             !        elem%estim_loc/ loc_tol, (elem%estim_loc/ loc_tol )**(1./elem%deg), &
             !        elem%reg / elem%regT0 , elem%reg / elem%regT2
          endif
       endif   ! if(.not. hp_adapt

       !if( dot_product(elem%xc(1:2), elem%xc(1:2) )**0.5 < 1E-3) then
       !   write(77,'(a4,3i5,10es12.4)') 'IAMA',state%space%adapt%adapt_level, elem%i, elem%deg, &
       !        elem%estim_loc, loc_tol, elem%reg, ratio, min(ratio, 2.5)**1.2
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
       !write(*,'(a6,i5,8es12.4)') '#ISO#',elem%i,elem%rgabc(:), lam_max, lam_actual,elem%estim_loc
    end do


    ! write(*,'(a10,10i7)') 'AdatG:',iA, iB, iC, iD, iE, iFF, iG, iH, iI
    ! if(state%space%adapt%adapt_level == 0) then
    !     write(84,'(a10,10a12)') '***:',&
    !          'E>T p+', ' E>T h+',' E<T p+',' E<T p-', 'E<T p0', ' ', ' ', ' ', ''
    !     write(84,'(a10,10a12)') '***:','iA', 'iB', 'iC','iD','iE', 'iFF', 'iG ', 'iH ', ' iI'
    !  endif

    ! write(84,'(a10,10i12)') 'AdatG:',iA, iB, iC, iD, iE, iFF, iG, iH, iI

    !print*
    !print*,'Sum', sum(grid%elem(:)%estim_loc**2)**0.5, tol_max
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
    !  open(imt, file=file1, status='UNKNOWN')
    !  do i=1,grid%nelem
    !     elem => grid%elem(i)
    !     !    if( mod(i, 3) == 1) &
    !     if( sqrt(dot_product(elem%xc(:), elem%xc(:))) < 0.1) &
    !     !    if(abs(elem%xc(1) - 1.5) < 0.25 .and. elem%xc(2) > 1.75 ) &
    !     call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )

    !  enddo
    !  close(imt)

    call SmoothMetric( )


    ! ! ! variant of high order Riemann metric
    ! imt = 24
    ! open(imt, file=file2, status='UNKNOWN')
    ! do i=1,grid%nelem !,2
    !    elem => grid%elem(i)
    !    !    !if( dot_product(elem%xc(:), elem%xc(:))**0.5 < 4E-2) &
    !    if( sqrt(dot_product(elem%xc(:), elem%xc(:))) < 0.1) &
    !         call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )
    ! enddo
    ! close(imt)


  end subroutine IsotropicMetric


  !> set the Metric for isotropic refinement, refinement only for the TOP elements
  !> refinement uses elem%estim_loc
  subroutine IsotropicMetric_SimpleOrdering( )
    class(element), pointer :: elem
    integer :: i, j, k,  imt, imt1, is, ityp, pK
    real :: tol_max
    real :: loc_tol, h_opt, lam_actual, lam_max
    real :: scale, weight, ratio, deref_level
    real :: limit_min, limit_max
    logical :: hp_adapt
    integer :: iA, iB, iC, iD, iE, iFF, iG, iH, iI
    character(len=15) :: file1, file2
    character(len=5) :: ch5
    real :: max_ratio, min_ratio, sum
    real , dimension(:), allocatable :: est
    integer , dimension(:), allocatable :: iest
    logical :: only_h_adapt, only_p_adapt
    real :: ratio_max, ratio_min, estim_max,  estim_min,  estim_tot

    max_ratio = 5.

    only_h_adapt = .false.
    only_p_adapt = .false.
    if(state%space%adapt%adapt_space == 'HGh') only_h_adapt = .true. ! only h-refinement
    if(state%space%adapt%adapt_space == 'HGp') only_p_adapt = .true. ! only p-refinement


    ! ordering of the elements according the values of estimator
    allocate(iest(1:grid%nelem), est(1:grid%nelem) )
    do i=1,grid%nelem
       iest(i) = i
       est(i) = grid%elem(i)%estim_loc

       elem => grid%elem(i)
       elem%psplit = 0
       elem%hsplit = 0
    enddo

    estim_max = maxval(grid%elem(:)%estim_loc)
    estim_min = minval(grid%elem(:)%estim_loc)
    estim_tot = sqrt(sum(grid%elem(:)%estim_loc**2))
    tol_max = state%space%adapt%tol_max

    write(*,'(a40, 6es12.4)')'IMSO estim_max, estim_min, estim_tot, tol_max:', &
         estim_max, estim_min, estim_tot, tol_max

    call order_estims(grid%nelem, iest(1:grid%nelem), est(1:grid%nelem))


    ! variant with elemets with error at lest 0.1 % of the errors
    !limit_max = 0.5 * estim_max
    !limit_min =   5 * estim_min

    !ratio_max = 1.0
    !ratio_min = 0.0


    ! variant with 10% of elemets ordered with inceasing error
    limit_max = 0.
    !limit_min = 5 * estim_min
    limit_min = 5 * estim_max   ! NO COARSENING

    ratio_max = 0.1    ! first 10% of elements refined
    ratio_min = 1 + 0.95  !0.9    ! last  10% of elements de-refined



    do j=1, grid%nelem !!  int(0.1 * grid%nelem)
       i = iest(j)
       elem => grid%elem(i)


       ! refinement
       if(elem%estim_loc >= limit_max .and. j <= int(ratio_max * grid%nelem) ) then
          if(only_h_adapt) then
             elem%hsplit = 4
          elseif(only_p_adapt) then
             elem%psplit = 1

          else
             if(elem%reg <= elem%regT0 .and. elem%deg < MaxDegreeImplemented-1 ) then

                elem%psplit = 1
                !write(100 + state%space%adapt%adapt_level + 1,*) elem%xc(:), grid%x(elem%face(idx, 1:3), 1:2)
                !write(*,'(a5,6es12.4)') 'p+', elem%xc(:), elem%estim_loc, elem%reg, elem%regT0

             else
                elem%hsplit = 4
                !write(200 + state%space%adapt%adapt_level + 1,*) elem%xc(:), grid%x(elem%face(idx, 1:3), 1:2)
                !write(*,'(a5,6es12.4)') 'h+', elem%xc(:), elem%estim_loc, elem%reg, elem%regT0

             endif
          endif

       ! de-refinement
       elseif(elem%estim_loc < limit_min .and. j >= int(ratio_min * grid%nelem) ) then

          if(only_h_adapt) then
             elem%hsplit = -4
          elseif(only_p_adapt) then
             elem%psplit = -1

          else
             if(elem%reg <= elem%regT0) then
                elem%hsplit = -1
                !write(300 + state%space%adapt%adapt_level + 1,*) elem%xc(:), grid%x(elem%face(idx, 1:3), 1:2)

             else
                if(elem%deg > MinDegreeImplemented)  elem%hsplit = -1  ! elem%psplit = -1
                !write(400 + state%space%adapt%adapt_level + 1,*) elem%xc(:), grid%x(elem%face(idx, 1:3), 1:2)
                !write(*,'(a5,6es12.4)') 'p-', elem%xc(:), elem%estim_loc, elem%reg, elem%regT0

             endif
          endif

       endif

       !write(100 + 10*state%space%adapt%adapt_level + 5 + elem%hsplit/4,*) elem%xc(:)
    enddo


    deallocate(iest, est)


    do i=1, grid%nelem
       elem => grid%elem(i)

       ratio = 1.
       if(elem%hsplit == 4) ratio = 2.  !2.5 ! 2.5

       ratio = min( ratio, max_ratio)
       !if(Algol3) ratio = max(ratio, min_ratio)

       ! diam based setting
       !h_opt = elem%diam / ratio
       !lam_max = 3./h_opt**2

       ! area based setting
       lam_max = 3*sqrt(3.)/4 / (elem%area / ratio /ratio)

       !if(elem%hsplit == 4) &
       !     write(*,'(a6,8es12.4)') 'DE#??', h_opt, elem%diam,  elem%area**0.5, elem%area**0.5/ elem%diam, &
       !     2*elem%area/ elem%diam, 2*elem%area/ elem%diam / elem%diam


       weight = 1.   !5
       elem%ama_p = weight * elem%psplit

       elem%rgabc(1) = lam_max
       elem%rgabc(2) = 0.
       elem%rgabc(3) = lam_max

       !write(*,'(a6,2i9,8es12.4)') '#ISO#',elem%i,elem%psplit, elem%ama_p , elem%rgabc(:), elem%estim_loc
       !stop

       !write(*,'(a6,i5,8es12.4)') '#ISO#',elem%i,elem%rgabc(:), lam_max, lam_actual,elem%estim_loc
    end do


    ! write(*,'(a10,10i7)') 'AdatG:',iA, iB, iC, iD, iE, iFF, iG, iH, iI
    ! if(state%space%adapt%adapt_level == 0) then
    !     write(84,'(a10,10a12)') '***:',&
    !          'E>T p+', ' E>T h+',' E<T p+',' E<T p-', 'E<T p0', ' ', ' ', ' ', ''
    !     write(84,'(a10,10a12)') '***:','iA', 'iB', 'iC','iD','iE', 'iFF', 'iG ', 'iH ', ' iI'
    !  endif

    ! write(84,'(a10,10i12)') 'AdatG:',iA, iB, iC, iD, iE, iFF, iG, iH, iI

    !print*
    !print*,'Sum', sum(grid%elem(:)%estim_loc**2)**0.5, tol_max
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


    ! ! ! variant of high order Riemann metric
    ! imt = 24
    ! open(imt, file=file2, status='UNKNOWN')
    ! do i=1,grid%nelem !,2
    !    elem => grid%elem(i)
    !    !if( dot_product(elem%xc(:), elem%xc(:))**0.5 < 4E-2) &
    !    call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )
    ! enddo
    ! close(imt)


  end subroutine IsotropicMetric_SimpleOrdering


  !> prescribe minimal refinement around the obtacles
  subroutine Refine_a_priori( )
    class(element), pointer :: elem
    integer :: i, j
    integer :: NN
    real, dimension(:,:), allocatable :: xc
    real, dimension(:,:), allocatable :: rc
    real ::  h_min, lam_min, rlen, hh, lam, h_max, lam_max, alpha

    NN = 1

    allocate( xc(1:NN, 1:2), rc(1:NN, 1:2) )  ! rc(*, 1) r_min, rc(*,2) r_max

    xc(1, 1:2) = (/ 35, 5 /)  ! center  (/ 35, 5 /)
    rc(1,1) = 2.0
    rc(1,2) = 4.0

    h_max = 0.5
    lam_min = 3./h_max**2
    lam = lam_min


    !h_min = 0.5
    !lam_max = 3./h_min**2

    do i=1,grid%nelem
       elem => grid%elem(i)

       do j=1, NN
          rlen = sqrt(dot_product(elem%xc(1:2) - xc(j, 1:2), elem%xc(1:2) - xc(j, 1:2)))

          if(rlen <= rc(j, 2)) then
             alpha = (rlen - rc(j, 1) ) / (rc(j, 2) - rc(j,1) )
             hh = h_max *( 1 + alpha)
             lam = 3./hh**2
             !print*,'ede3d', rlen, h_max, hh, alpha

             ! a priori refinement
             !if( elem%rgabc(1) < lam .and.  elem%rgabc(3) < lam ) then
                !elem%rgabc(1) = lam *( 1- alpha)  +  alpha * elem%rgabc(1)
                !elem%rgabc(2) =                      alpha * elem%rgabc(2)
                !elem%rgabc(3) = lam *( 1- alpha)  +  alpha * elem%rgabc(3)

                elem%rgabc(1) = lam
                elem%rgabc(2) = 0.
                elem%rgabc(3) = lam
             !endif

             !elem%rgabc(1) = max(lam, elem%rgabc(1))
             !elem%rgabc(3) = max(lam, elem%rgabc(3))

             !elem%rgabc(1) = min(lam_max, elem%rgabc(1))
             !elem%rgabc(3) = min(lam_max, elem%rgabc(3))
             !endif
          else
          !   elem%rgabc(1) = lam / 50
          !   elem%rgabc(2) = 0.
          !   elem%rgabc(3) = lam / 50
          endif

       enddo
    enddo

    deallocate(xc, rc)

  end subroutine Refine_a_priori

end module ama_hp_interpol
