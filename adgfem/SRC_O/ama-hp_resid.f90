!> subroutines for hp method based on residuals estimates
module ama_hp_resid

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
  use ama_hp_interpol

  implicit none


  public:: AnisotResidEstimates
  public:: Eval_hp_resid_Metric
  public:: Eval_hp_resid_MetricElem
  public:: Set_Optimal_HO_degree
  public:: Eval_HO_Derivatives
contains

  !> perform the anisotropic error estimates using residual (RES)
  !> or dual weighted residual (DWR) error estimates 
  subroutine AnisotResidEstimates( )
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

    if(.not. grid%ElemSupports) &
         call SeekElemSupports(grid)  ! create the list of elements sharing at least a vertex with elem

    call SetQuantities4Metric(ndimL )

    !stop

    state%err(algeb) = 0.

    state%num_limits = 0
    call Eval_hp_resid_Metric( ndimL )  !!!wp(1:grid%npoin, 1:ndim) )

    !deallocate(wp)

    file1 = 'metrixA00000'
    file2 = 'metrixS00000'

    is = 0
    if(state%space%adapt%adapt_level > 0) is = int(log(1. * state%space%adapt%adapt_level)/log(10.))

    write( ch5, '(i5)' ) state%space%adapt%adapt_level  ! change the format if num_size /= 5 !!!
    !!!write( ch5, '(i5)' ) state%time%iter  ! change the format if num_size /= 5 !!!
    file1(12-is: 12)  = ch5(5-is:5)
    file2(12-is: 12)  = ch5(5-is:5)


    ! variant of high order Riemann metric
    imt = 24
    open(imt, file=file1, status='UNKNOWN')
    do i=1,grid%nelem !, 2
       elem => grid%elem(i)
       ! !!write(*,'(a8,i5,8es12.4)') 'ede4r',elem%i, elem%rgabc(1:3), elem%xc(1:2)
       !  !    if( mod(i, 3) == 1) &
       ! !     if( dot_product(elem%xc(:), elem%xc(:))**0.5 < 4E-3) &
       ! ! !    if(abs(elem%xc(1) - 1.5) < 0.25 .and. elem%xc(2) > 1.75 ) &
       !if(elem%xc(1) < -10 .and. elem%xc(1) < -1 .and.  &
       !        elem%xc(2) > -30 .and. elem%xc(2) < 30) &
       call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )

    enddo
    close(imt)

    !print*,'STOPPED IN ESRTWGQ'
    !stop

    if(state%modelName == 'pedes' ) call Refine_a_priori( )

    call SmoothMetric( )


    ! ! variant of high order Riemann metric
    ! imt = 24
    ! open(imt, file=file2, status='UNKNOWN')
    ! do i=1,grid%nelem !,2
    !    elem => grid%elem(i)
    !    call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )
    ! enddo
    ! close(imt)

    ! deallocation of arrays allocated in SeekElemSupports
    do i=1,grid%nelem
       elem => grid%elem(i)
       !write(*,'(a6, i5, 6es12.4)') 'ede3',i, elem%xc, elem%rgabc
       !deallocate(elem%supp)
       deallocate(elem%wS)
    enddo

    !print*,'stoped in ama-hp_resid.f90'
    !stop

  end subroutine AnisotResidEstimates



  !> evaluate the Riemann metric for hp-mesh
  subroutine Eval_hp_resid_Metric( ndimL )    !!!wp)
    use ama_hp_interpol_params
    !!real, dimension(1:grid%npoin, 1:ndim), intent(inout) :: wp
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    class(element), pointer :: elem
    integer :: i, dof, k, dofL

    !evaluate of the parameters for the ama_hp_adaptation
    if(state%space%estim_space == 'RES' .or. state%space%estim_space == 'DWR') call Eval_AMA_paramets( ) 

    ! sorting of elements according the maximal estim_loc
    call  Sort_maximal_elements( ) 

    ! filling of the arrays by the L^2-projection
    call Assemble_wSS_arrays( ) ! allocations elem%wSS(1, 0:1, 1:dof*ndim)

    do i=1,grid%nelem
       elem => grid%elem(i)
       call Eval_hp_resid_MetricElem(elem, ndimL)   
    enddo

    !deallocation of arrays
    call Sort_maximal_elements_dealloc( )
    !stop 'end subroutine Eval_hp_resid_MetricElem'


    call Deallocate_wSS_arrays( ) 

  end subroutine Eval_hp_resid_Metric


  !> evaluate the Riemann metric for hp-mesh for one element
  subroutine Eval_hp_resid_MetricElem(elem, ndimL)  
    use ama_hp_interpol_params
    class(element), intent(inout) :: elem
    integer, intent(in):: ndimL   ! number of quantities for metric evaluation
    real, dimension(:,:,:,:), allocatable :: Dw    ! all derivatives of w: bases coefficients
    real, dimension(:,:,:,:), allocatable :: Dwi   ! all derivatives of w: integ nodes
    real, dimension(:), allocatable :: Kparams 
    real, dimension(:,:), allocatable :: xi, yi
    real :: slope_h, slope_p, err_h_deref, err_p_deref
    real :: area, target_tol , Lq, factor, tols, lam_max
    real :: elem_reg, reg_param
    integer :: i, j, k, rn
    integer :: deg, dof, degP, dofP, Qdof
    integer :: itype, is
    character(len=13) :: p_split, h_split
    character(len=5) :: ch5
    logical :: singularity
 
    !itype = 0  ! a priori knowledge of the singularity
    !itype = 1  ! equidistribution
    !itype = 2  ! top elements
    !itype = 3  ! top elements with RES strategy
    !itype = 4  ! top elements with RES strategy + hp-decay (SISC16) indication
    itype = 5  ! top elements with RES strategy + full hp-decay indication
    
    elem%interLq = 0.
    elem%interL8 = 0.

    if(elem%flen /= 3) then
       write(*,*)' subroutine Eval_hp_MetricElem only for triangles without HG nodes'
       stop
    endif

    Qdof = elem%Qdof
    degP = elem%deg + 1
    dofP = DOFtriang( degP )

    ! Dw(k,i,j, l) => k-th component of w, (i,j) d^i/(dx^j dy^(i-j)), l coeff or integ node
    allocate(Dw (1:ndimL, 0:degP, 0:degP, 1:dofP) )
    allocate(Dwi(1:ndimL, 0:degP, 0:degP, 1:Qdof) )

    ! evaluation of the all derivatives up to the degree degP from the HO reconstruction
    call Eval_HO_Derivatives(elem, ndimL, degP, dofP, Dw (1:ndimL, 0:degP, 0:degP, 1:dofP), &
        Dwi(1:ndimL, 0:degP, 0:degP, 1:Qdof) )

    ! setting of the optimal degree from the HO reconstruction
    call Set_Optimal_HO_degree(elem, ndimL, degP, dofP, Dw(1:ndimL, 0:degP, 0:degP, 1:dofP)) 


    ! regularity parameter :: elem_reg <= reg_param  ==> p-refinement
    reg_param = 0.2

    elem_reg = ( elem%eta(resS, 1) - elem%eta(resSr, 1)) / max(1E-15, elem%eta(resSr, 1))

    target_tol = state%space%adapt%tol_max

    ! setting of the target area
    if(itype == 0) then  !  a priori knowledge of the singularity

       if(elem%i == 1) print*,'! full hp-derefinement indication in ama-hp_resid'

       call Detect_apriori_known_singularity(elem, singularity) 

       factor = 1. 
       elem%psplit = 0

       ! refinement
       do j=1, int(0.1 * grid%nelem)

          if(elem%i == ama_iest(j)  ) then
             if( .not. singularity) then 
                elem%psplit = 1
                !!factor = 1.1
             else
                factor = 4.
             endif
             write(*,'(a12,3i5,4es12.4)') 'marked APR', elem%i, elem%deg, elem%psplit, factor
          endif
       enddo

    else if(itype == 1) then ! error equidistribution
       target_tol = state%space%adapt%tol_min
       !target_tol = ama_err_total / 10.
       !target_tol = min( state%space%adapt%tol_max, ama_err_total / 5.)
       
       Lq = 1./state%space%adapt%Lq
       tols = target_tol * (elem%area/  state%space%domain_volume)**Lq
       factor = elem%estim_loc / tols
       factor = factor*(2./(max(1E-2, 1.*elem%deg)))

    elseif(itype == 2) then  ! marking of the top elements
       factor = 1. 
       do j=1, int(0.1 * grid%nelem)
          if(elem%i == ama_iest(j) .and. elem%psplit <= 0. ) factor = 4.
       enddo

    elseif(itype == 3) then  ! marking of the top elements
       factor = 1. 
       elem%psplit = 0

       ! refinement
       do j=1, int(0.1 * grid%nelem)

          if(elem%i == ama_iest(j)  ) then
             if( elem_reg < reg_param &
                  .and. elem%deg < MaxDegreeImplemented-2 ) then
                elem%psplit = 1
                !!factor = 1.1
             else
                factor = 4.
             endif
             !write(*,'(a12,3i5,4es12.4)') 'marked', elem%i, elem%deg, elem%psplit, factor
          endif
       enddo

          
       ! derefinement
       do j= int(0.1 * grid%nelem)+1,  int(0.9 * grid%nelem)
          if(elem%i == ama_iest(j)  ) then
             if( elem_reg < reg_param &
                  .and. elem%deg < MaxDegreeImplemented-2 ) then
                !elem%psplit = 1
                !factor = 0.8
             !else
             !   factor = 0.75
             endif
             !write(*,'(a12,3i5,4es12.4)') 'marked', elem%i, elem%deg, elem%psplit, factor
          endif
       enddo

       ! derefinement
       do j= int(0.9 * grid%nelem)+1,  grid%nelem
          if(elem%i == ama_iest(j)  ) then
             if( elem_reg < reg_param &
                  .and. elem%deg < MaxDegreeImplemented-2 ) then
                !elem%psplit = 1
                !factor = 0.5
             else
                !factor = 0.75
             endif
             !write(*,'(a12,3i5,4es12.4)') 'marked', elem%i, elem%deg, elem%psplit, factor
          endif
       enddo

    elseif(itype == 4) then  ! marking of the top elements + hp-decay regularity indication

       if(elem%i == 1) print*,'! hp-derefinement indication in ama-hp_resid'

       call Regularity_hp_derefinement(elem, err_h_deref, err_p_deref)

       !slope_h = (err_h_deref - elem%eta(P_tot, 1) ) / (1.5*(elem%deg +1)*(elem%deg + 2) )
       !slope_p = (err_p_deref - elem%eta(P_tot, 1) ) / ( elem%deg +1 )


       rn = 3
       allocate(xi(1:rn, 1:2) )
       
       ! p-1, h
       xi(1, 1) = err_p_deref                            ! error
       xi(1, 2) = 1.*(elem%deg + 0) * (elem%deg + 1) / 2    ! DOF
       
       ! p, 2 * h 
       xi(2, 1) = err_h_deref                            ! error
       xi(2, 2) = 1.*(elem%deg + 1) * (elem%deg + 2) / 8    ! DOF

       ! p, h
       xi(3, 1) = elem%eta(resS, 1)    ! NONSENCE !!!!!!                     ! error
       xi(3, 2) = 1.*(elem%deg + 1) * (elem%deg + 2) / 2    ! DOF


       ! EOC for h-refinement
       slope_h = log(xi(2, 1) / xi(3, 1)) /( xi(3, 2)**(1./3) - xi(2, 2)**(1./3) )

       ! EOC for p-refinement
       slope_p = log(xi(1, 1) / xi(3, 1)) /( xi(3, 2)**(1./3) - xi(1, 2)**(1./3) )


       ! write(100 + state%space%adapt%adapt_level, *)  xi(3, 2), xi(3, 1)
       ! write(100 + state%space%adapt%adapt_level, *)  xi(1, 2), xi(1, 1)
       ! write(100 + state%space%adapt%adapt_level, *) '   '

       ! write(200 + state%space%adapt%adapt_level, *)  xi(3,2), xi(3, 1)
       ! write(200 + state%space%adapt%adapt_level, *)  xi(2,2), xi(2, 1)
       ! write(200 + state%space%adapt%adapt_level, *) '   '          
  
       !write(*,'(a7, i5, 6es12.4)') 'regul:',elem%i, elem%eta(resS, 1), err_h_deref, err_p_deref,slope_h, slope_p

       !if(slope_h >  slope_p) & 
       !     !if(err_h_deref <  err_p_deref) & 
       !     write(*,'(a8, 2i5,3es12.4,3(a2, 3es12.4))') &
       !     'decayL',elem%i, elem%dof, elem%eta(resS, 1), err_h_deref, err_p_deref, &
       !     '|',xi(1:3, 1), '|', xi(1:3, 2)
            

       deallocate(xi)


       factor = 1. 
       elem%psplit = 0

       ! refinement
       do j=1, int(0.1 * grid%nelem)

          if(elem%i == ama_iest(j)  ) then
             if( slope_p >= slope_h .and. elem%deg < MaxDegreeImplemented-2 ) then
                elem%psplit = 1
                !!factor = 1.1
             else
                factor = 4.
             endif
             !write(*,'(a12,3i5,4es12.4)') 'marked', elem%i, elem%deg, elem%psplit, factor
          endif
       enddo

    elseif(itype == 5) then  ! marking of the top elements + full hp-decay regularity indication

       if(elem%i == 1) print*,'! full hp-derefinement indication in ama-hp_resid'


       call Regularity_full_hp_derefinement(elem, slope_h, slope_p)

        !if(slope_h <  slope_p) then
       !    write(*,'(a7, i5,es12.4,a2, 16es12.4)') 'regul:',elem%i, elem%eta(resS, 1), '|', &
       !         err_p_deref, err_h_deref, err_hp_deref, slope_p, slope_h, xi(1:3, 2)
       ! !endif

       !if(elem%i == 5) stop '8yr39jdo3w'
       p_split = 'p_split-00000'
       h_split = 'h_split-00000'

       is = int(log(1.*state%space%adapt%adapt_level+0.1)/log(10.)) 
       write( ch5, '(i5)' ) state%space%adapt%adapt_level
       p_split(13-is:13) = ch5(5-is : 5)
       h_split(13-is:13) = ch5(5-is : 5)
 
       open(11, file=p_split, status='UNKNOWN', position='append')
       open(12, file=h_split, status='UNKNOWN', position='append')
       if( slope_p >= slope_h ) write(11, *) elem%xc
       if( slope_p <  slope_h ) write(12, *) elem%xc
       close(11)
       close(12)


       factor = 1. 
       elem%psplit = 0

       ! refinement
       do j=1, int(0.1 * grid%nelem)

          if(elem%i == ama_iest(j)  ) then
             if( slope_p >= slope_h .and. elem%deg < MaxDegreeImplemented-2 ) then
                elem%psplit = 1
                !!factor = 1.1
             else
                factor = 4.
             endif
             !write(*,'(a12,3i5,4es12.4)') 'marked', elem%i, elem%deg, elem%psplit, factor
          endif
       enddo

    else
       stop "UNKNOWN itype in ama-hp_resid.f90 73d3ie"
    endif

    factor = min(factor, 4.)
    elem%reg2 = factor
    area = elem%area / factor

    !if(factor > 2.) write(state%space%adapt%adapt_level * 10 + 10 ,*) elem%xc
    !if(elem%psplit == 1 ) write(state%space%adapt%adapt_level * 10 + 11 , *) elem%xc
    !if(elem%psplit == -1 ) write(state%space%adapt%adapt_level * 10 + 12 , *) elem%xc


    if((elem%i == 1 .or. elem%i == grid%nelem) ) & !.and.  target_tol < state%space%adapt%tol_max) &
         write(*,'(a12,4es12.4,a2,8es12.4)') &
         'AMA_pars ii:', ama_err_min, ama_err_aver, ama_err_max,ama_err_total,'|', &
         target_tol, state%space%adapt%tol_max
    
    !write(*,'(15es12.4)') ,  elem%xc(:), tols, elem%eta(resS, 1), elem%eta(resS, 1) / tols, factor


    if(itype == 3 .or. itype == 4  .or. itype == 5) then
       ! area based setting
       lam_max = 3*sqrt(3.)/4 / (elem%area / factor)

       elem%ama_p = elem%psplit ! * 1.5

       elem%rgabc(1) = lam_max
       elem%rgabc(2) = 0.
       elem%rgabc(3) = lam_max

    else

       !print*,' setting of the anisotropy'
       allocate(Kparams(1:5) )
       call FindAnisotropy_AreaFixed( elem, ndimL, degP, Dwi(1:ndimL, degP, 0:degP, 1), &
            area, Kparams(1:5)) 

       !print*,' setting of the metric for the anisotropy'
       call Setting_Metric_Anisotropy(elem, Kparams(1:5) )
       
       deallocate(  Kparams)

    endif

    deallocate(Dw, Dwi)

  end subroutine Eval_hp_resid_MetricElem

  !> setting of the metric for the anisotropy
  subroutine Setting_Metric_Anisotropy(elem, Kparams )
    type(element), intent(inout) :: elem
    real, dimension(1:5), intent(in) :: Kparams
    real :: h_min, h_max, lambda_min, lambda_max
    real :: lam_max, lam_min, max_a,  epsilon1, ratio
    real :: lam_actual, lam_actualT, lam_max_actual_ratio

    ! metric limitation
    h_max = min( state%space%diam / 4., 4.)
    !h_max = state%space%diam / 2.
    epsilon1 = 1E+12

    h_min = h_max / epsilon1**0.5

    lambda_min = 1./h_max**2
    lambda_max = 1./h_min**2

    ! actual metric
    lam_max = Kparams(1)
    lam_min = Kparams(2)
    max_a   = Kparams(3)

    elem%ama_p = elem%psplit


    ! metric limitation
    lam_max = max(lambda_min, min(lambda_max, lam_max) )
    lam_min = max(lambda_min, min(lambda_max, lam_min) )
    
    ratio = lam_min / lam_max
    ratio = max(ratio, 2E-5)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
 
    ! limitation of not too fast refinement
    lam_actual = (2*elem%area/elem%diam)**(-2)
    lam_actualT = (elem%diam)**(-2)
    lam_max_actual_ratio = 4. !10
    !if(.not. state%time_dependent) lam_max_actual_ratio = 50.

    lam_max = min(lam_max, lam_actual * lam_max_actual_ratio)  ! at most r^{1/2} times smaller element
    lam_max = max(lam_max, lam_actual / lam_max_actual_ratio)  ! at most r^{1/2} times bigger element

    ! the consequent limitation of the ratio
    ratio = max(ratio,  lam_actualT / lam_actual / 2)

    !if(lam_max >=  lam_actual * lam_max_actual_ratio) &
    !     write(*,'(a8,3i5, 3es12.4,a2,30es12.4)') 'limits:',elem%i, elem%deg, elem%deg+ elem%psplit, &
    !     lam_min /ratio, lam_min, ratio, '|',  &
    !     lam_actual, lam_max_actual_ratio, lam_actual * lam_max_actual_ratio, lam_max

    ! new lam_min after limitation of the aspect ratio and lam_max
    lam_min = lam_max * ratio


    elem%rgabc(1) =  lam_max * cos(max_a)**2 + lam_min * sin(max_a)**2
    elem%rgabc(2) = (lam_max - lam_min) * cos(max_a) * sin(max_a)
    elem%rgabc(3) =  lam_max * sin(max_a)**2 + lam_min * cos(max_a)**2

  end subroutine Setting_Metric_Anisotropy


  !> setting of the optimal polynomial dregree using the approach [Dol, Solin: AMC 2016]
  subroutine Set_Optimal_HO_degree(elem, ndimL, degP, dofP, Dw)
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL, degP, dofP  ! ndimL, number of quantities for metric evaluation
    ! Dw(k,i,j, l) => k-th component of w, (i,j) d^i/(dx^j dy^(i-j)), l coeff or integ node
    real, dimension(1:ndimL, 0:degP, 0:degP, 1:dofP) , intent(inout) :: Dw

    type(volume_rule), pointer ::V_rule
    real, dimension(:), allocatable :: weights
    real, dimension(:,:,:), allocatable :: wi
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:,:), allocatable :: Dphi
    real, dimension(:,:), allocatable :: qq
    real :: m_val
    integer :: Qdof, i, j, k, ik, ideg, degPP, dofPP, dof
    integer :: iff, iff1, i_val


    !k = 1
    !do i=0,degP   ! i = i-th derivative
    !   do j=0,i  !   d^i w /dx^j dy^(i-j)
    !      write(*,'(a6,3i5,120es12.4)') 'Dw:',i,j,dofP,Dw(k, i, j, 1:dofP)
          
          !call PlotElemFunction3D(100+i*10+j, elem,  dofP, Dw(1, i, j, 1:dofP) )
          
    !   enddo
    !enddo
    !write(*,*) "################################################"


    !call PlotElemFunction3D(100, elem,  dofP, Dw(1, 0, 0, 1:dofP) )

    dof = elem%dof

    V_rule => state%space%V_rule(elem%Qnum)
    Qdof = V_rule%Qdof

    ! quadrature weights
    allocate(weights(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof))
    
    ! basis functions in integ nodes
    phi => V_rule%phi(1:dofP, 1:Qdof)

    ! derivatives of basis functions in integ nodes
    allocate(Dphi(1:dofP, 1:nbDim, 1:Qdof) )
    call  Eval_Dphi_plus(elem, V_rule, dofP, Dphi(1:dofP, 1:nbDim, 1:Qdof) )

    ! projection  and its derivatives  in integ nodes
    allocate(wi(1:15, 1:ndim, 1:Qdof) )


    ! HO reconstruction in integ nodes on elements
    do k=1, ndim
       ! w_rec
       wi(1, k, 1:Qdof) = matmul(Dw(k, 0, 0, 1:dofP), phi(1:dofP, 1:Qdof) )
       
       ! D_x w_rec, D_y w_rec
       wi(2, k, 1:Qdof) = matmul(Dw(k, 0, 0, 1:dofP), Dphi(1:dofP, 1, 1:Qdof) )
       wi(3, k, 1:Qdof) = matmul(Dw(k, 0, 0, 1:dofP), Dphi(1:dofP, 2, 1:Qdof) )
    enddo
    
    do ideg = 0, -2, -1
       ! projection of the solution 
       degPP = elem%deg + ideg      ! used for the projection of the approximate solution
       dofPP = DOFtriang( degPP)

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

       if(ideg == -2 )  iff = res_HO_p2 
       if(ideg == -1 )  iff = res_HO_p1
       if(ideg ==  0 )  iff = res_HO_p0

       ! ATTENTION, MISCHMATCH ndim AND ndimL
       do k=1,ndim
          ! error estimate in the H1 semi-norm
          elem%eta(iff, k) =  & ! elem%eta(iff+2, k) 
               + dot_product(weights(1:Qdof),  &
               ( wi(2, k, 1:Qdof)- wi(5, k, 1:Qdof))**2 + ( wi(3, k, 1:Qdof)- wi(6, k, 1:Qdof))**2 )
       enddo

       !write(*,'(a6,4i5, 500es12.4)') 'etas',ideg, degPP, dofPP, iff, elem%eta(iff, 1)**0.5, &
       !     wi(2, 1, 1:5),  wi(5, 1, 1:5)

    enddo

    !stop "73yrh39ide"

    iff1 = res_HO_p0

    m_val = 1E+30
    ! array for hp-parameters
    allocate(qq(-2 : 0, 1:5) )
    do j=-2, 0
       if(j == -2 )  iff = res_HO_p2 
       if(j == -1 )  iff = res_HO_p1
       if(j ==  0 )  iff = res_HO_p0

       qq(j, 1) = sqrt(elem%eta(iff, 1) / elem%eta(iff1, 1) )

       !!if(elem%i == 1) write(*,'(a10,6i5)') 'adapt iff =',iff, iff1

       if(elem%deg +j  > 0 ) then
          qq(j, 2) = qq(j,1) **(2.0/(elem%deg + j))  
       else
          qq(j, 2) = 1E+10
       endif
       
       ! the necessary number of degrees of freedom 
       qq(j,3) = qq(j,2) * DOFtriang(elem%deg + j)
       
       if(qq(j,3)  < m_val) then
          m_val = qq(j, 3)
          i_val = j
       endif

       !write(*,'(a6,4i5, 5es12.4, i5)') 'regul',j, elem%deg+j, iff,iff1, elem%eta(iff, 1)**0.5, &
       !     qq(j,1), qq(j,2), qq(j,3), m_val, i_val
       !if(j==0) print*,'----------------------------------'
       
    enddo

    ! setting of elem%psplit
    elem%psplit = i_val + 1    ! -1, 0, 1

    deallocate(Dphi, wi, weights, qq)

  end subroutine Set_Optimal_HO_degree

  !> evaluate all HO partial derivatives of \f$ w \f$ of degree degP on elem
  subroutine Eval_HO_Derivatives(elem, ndimL, degP, dofP, Dw, Dwi)
    type(element), intent(inout) :: elem
    integer, intent(in):: ndimL, degP, dofP  ! ndimL, number of quantities for metric evaluation
    ! Dw(k,i,j, l) => k-th component of w, (i,j) d^i/(dx^j dy^(i-j)), l coeff or integ node
    real, dimension(1:ndimL, 0:degP, 0:degP, 1:dofP) , intent(inout) :: Dw 
    real, dimension(1:ndimL, 0:degP, 0:degP, 1:elem%Qdof) , intent(inout) :: Dwi 

    real, dimension(:,:,:), allocatable :: Dphi ! derivative of the test functions
    real, dimension(:,:), allocatable :: MassElem ! Inversion of local mass matrix of order deg+1!!
    !real, dimension(:,:), allocatable :: Mass, MassS ! local mass matrix of order deg+1!!
    integer :: Qnum, Qdof, dof
    integer :: ideg, ideg_min, ideg_max, i, j, k, l,  ix, jx

    
    Qdof = elem%Qdof

    if(degP == -1) then
       ! direct computation of the second order derivatives
       print*,'VERIFY deduy83'
       call P_0_to_P2_Interpolation(elem, ndimL, .false., Qnum, degP, dofP, &
            Dw (1:ndimL, 0:degP, 0:degP, 1:dofP), -1 )

    else
       if(degP <= MaxDegreeImplemented) then

          ! copy of the local mass matrix used in the solution of Ax = b
          allocate(MassElem(1:dofP, 1:dofP) )

          ! Dw(k,i,j, l) => k-th component of w, (i,j) d^i/(dx^j dy^(i-j)), l coeff in integ node
          Dw = 0.
          Dwi = 0.


          ! higher order reconstruction via least square method
          call LeastSquareInterpolationWeighted(elem, ndimL, .false., elem%Qnum, degP, dofP, &
               Dw (1:ndimL, 0:degP, 0:degP, 1:dofP), -1 )  ! last argument=0 -> L2, = -1 -> H1

          !if(elem%i == 1) print*,'fort.11 written HEDt53d'
          !if(elem%xc(1) > 19. .and. elem%xc(1) < 21) &
          !     !print*, elem%i , elem%xc
          !call PlotElemFunction3D(11, elem,  elem%dof, elem%wS(1, 1:elem%dof) )
          
          !if(elem%i == 1) print*,'fort.21 written HEDt53d'
          !if(elem%xc(1) < 5.) &
          !if(elem%i == -599) &
          !call PlotElemFunction3D(21, elem,  dofP, Dw(1, 0, 0, 1:dofP) )
          
          
          allocate(Dphi(1:dofP, 1:nbDim, 1:Qdof ) )
          call  Eval_Dphi(elem, dofP,  Dphi(1:dofP, 1:nbDim, 1:Qdof) )

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
                   call IntegrateVectorB(elem, dofP, Dwi(k, i, j, 1:Qdof), Dw(k, i,j, 1:dofP) )
                   
                   !local mass matrix of order dofP x dofP
                   MassElem(1:dofP, 1:dofP) = elem%Mass%Mb(1:dofP, 1:dofP)
                   call SolveLocalMatrixProblem(dofP, MassElem(1:dofP, 1:dofP), 1, Dw(k, i,j, 1:dofP) )

                enddo ! j
             enddo ! i
          enddo ! k
          
          if(elem%i <= -10) then
             !if( abs(elem%xc(1) -20 ) < 0.5 ) then
             write(*,*) '-----', elem%i, '  ----', elem%xc
             k = 1
             do i=1,degP   ! i = i-th derivative
                do j=0,i  !   d^i w /dx^j dy^(i-j)
                   write(*,'(a6,3i5,120es12.4)') 'Dwi:',i,j,Qdof,Dwi(k, i, j, 1:Qdof)
                enddo
             enddo
             
             write(*,*) '-----------------------------------'
             
             do i=0,degP   ! i = i-th derivative
                do j=0,i  !   d^i w /dx^j dy^(i-j)
                   write(*,'(a6,3i5,120es12.4)') 'Dw:',i,j,dofP,Dw(k, i, j, 1:dofP)

                   call PlotElemFunction3D(100+i*10+j, elem,  dofP, Dw(1, i, j, 1:dofP) )

                enddo
             enddo
             
             write(*,*) '########################################'
          endif
          
          
       else
          stop '! this degree is not implemented, hence too high derivative f437yfh3'
          
       endif
       
       deallocate( MassElem, Dphi)

    endif ! if(degP == 1)

  end subroutine Eval_HO_Derivatives



  !>  assembling arrays elem%wSS with the projections
  !> the temporary arrays contains the following values
  !> elem%wSS(1, 0, :) = actual solution
  !> elem%wSS(1, 1, :) = p-1 projection
  subroutine Assemble_wSS_arrays( ) 
    class(element), pointer :: elem
    integer :: i, dof, k, dofL
    do i=1,grid%nelem
       elem => grid%elem(i)
       dof = elem%dof
       allocate( elem%wSS(1, 0:1, 1:dof*ndim) )
       elem%wSS = 0.

       ! original solution
       elem%wSS(1, 0, 1:dof*ndim) = elem%w(0, 1:dof*ndim) 
       
       ! p-1 projection
       dofL = dof - elem%deg - 1
       do k=1, ndim
          elem%wSS(1, 1, (k-1)*dof + 1: (k-1)*dof + dofL) = elem%w(0, (k-1)*dof + 1: (k-1)*dof + dofL) 
       enddo
       !write(*,'(a8,i5, 40es12.4)') '%w :', elem%i, elem%w(0,:)
       !write(*,'(a8,i5, 40es12.4)') '%wS0 :', dof,  elem%wSS(1, 0,:)
       !write(*,'(a8,i5, 40es12.4)') '%wS1 :', dofL, elem%wSS(1,1,:)
       !print*,'---'

    enddo
  end subroutine Assemble_wSS_arrays

  !> deallocation of the temporary arrays
  subroutine Deallocate_wSS_arrays( ) 
    class(element), pointer :: elem
    integer :: i
    do i=1,grid%nelem
       elem => grid%elem(i)
       deallocate( elem%wSS )
    end do
  end subroutine Deallocate_wSS_arrays

end module ama_hp_resid
