!> estimates for AMA adaptation
module AMA_estims

  use main_data  ! contains type(mesh) :: grid for computation
  use problem_oper
  use euler_problem
  use estimates
  use plot_geom
  use eval_sol
  use ama_hp_interpol

  implicit none

  public:: AnisotErrorEstimates
  public:: AnisotElemEstimate
  public:: EvalHessians
  public:: EvalMaximumLangAMA
  public:: EvalMaximumLangAMABeta
contains

  !> perform the anisotropic error estimates using the Hessian
  !> the commented alternative is the dual norm including (non-)linear algebraic error
  subroutine AnisotErrorEstimates( )
    class(element), pointer :: elem
    real, dimension(:,:), allocatable :: wp ! array  with solution in vertexes
    integer :: i, j, k,  imt, imt1, is
    logical :: loc_implicitly
    character(len=15) :: file1, file2
    character(len=5) :: ch5

    if(nbDim /=2 ) then
       print*,' ## AnisotErrorEstimates only for nbDim = 2!'
       stop
    endif

    ! variant of high order Riemann metric
    !allocate( wp(1:grid%npoin, 1:ndim) )
    !call EvalSolutionVertexes(wp(1:grid%npoin, 1:ndim) )
    !call Eval_hp_Metric(wp(1:grid%npoin, 1:ndim) )
    !deallocate(wp)
    !call SmoothMetric( )
    ! end of variant of high order Riemann metric


    file1 = 'metrixA00000'
    file2 = 'metrixS00000'
    !file2 = 'rgabcxA00000'

    is = 0
    if(state%space%adapt%adapt_level > 0) is = int(log(1. * state%space%adapt%adapt_level)/log(10.))

    write( ch5, '(i5)' ) state%space%adapt%adapt_level  ! change the format if num_size /= 5 !!!
    file1(12-is: 12)  = ch5(5-is:5)
    file2(12-is: 12)  = ch5(5-is:5)


    ! variant of high order Riemann metric
    !imt = 24
    !open(imt, file=file1, status='UNKNOWN')
    !do i=1,grid%nelem
    !   elem => grid%elem(i)
    !   call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )
    !enddo
    !close(imt)
    !return
    ! end of variant of high order Riemann metric


    !imt = 24
    !open(imt, file=file1, status='UNKNOWN')
    !imt1 = 25
    !open(imt1, file=file2, status='UNKNOWN')



    !!call RezidErrorEstimates( )  ! evaluate elem%estimS, SetvectorFields, etc.

    call EvalHessians ( ) ! elem%rgabc contains the second order derivative terms

    do i=1,grid%nelem
       elem => grid%elem(i)
       do j=3,3
          call AnisotElemEstimate(elem, j)  ! element residuum in X- norm
       enddo

       !call DrawEllips(imt, elem%rgabc(1:3), elem%xc(1:2) )

    enddo

    call SmoothMetric1( )

    !do i=1,grid%nelem
    !   elem => grid%elem(i)

    !   !write(*,* ) elem%i,  elem%rga, elem%rgabc(2), elem%rgabc(3)

    !   call DrawEllips(imt1, elem%rgabc(1:3), elem%xc(1:2) )
    !enddo

    !close(imt)
    !close(imt1)

  end subroutine AnisotErrorEstimates

  !> smoothing of the metric stored ????????
  subroutine SmoothMetric1( )
    class(element), pointer :: elem, elem1
    real :: est, estT
    real, dimension(1:2) :: grad
    real, dimension(:,:), allocatable :: rgabc
    real :: weigth
    integer :: i, j, k, ipoc

    ipoc = 3
    weigth = 1.0
    !weigth = 0.5
    !weigth = 0.25

    allocate( rgabc(1:grid%nelem, 1:4) )

    do k=1, ipoc
       rgabc(:,:) = 0.

       do i=1,grid%nelem
          elem => grid%elem(i)

          rgabc(i,1:3) = elem%rgabc(1:3)
          rgabc(i,4) = 1.

          do j=1,elem%flen
             if(elem%face(neigh,j) > 0) then
                elem1 => grid%elem(elem%face(neigh,j))

                rgabc(i,1:3) = rgabc(i,1:3) + elem1%rgabc(1:3) * weigth
                rgabc(i,4) = rgabc(i,4) + weigth
             endif
          enddo
       enddo

       do i=1,grid%nelem
          elem => grid%elem(i)
          elem%rgabc(1:3) = rgabc(i,1:3) / rgabc(i, 4)
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
  end subroutine SmoothMetric1


  !> estimation of the element reziduum,
  !> \f$ \eta_K := \max_{\| \varphi\|_X=1}
  !> \frac{|c_h({\bf w}, \varphi_i) |}{\|\varphi\|_X}\f$,
  !> \f$\varphi\in P^{p_K+1}(K)\f$; \
  !> ityp = 1  \f$ \Rightarrow \| v \|_X = \| v \|_{L^2} \f$,
  !> ityp = 2  \f$ \Rightarrow \| v \|_X = \| v \|_{H^1} \f$,
  !> ityp = 3  \f$ \Rightarrow \| v \|_X^2 = \|u\|_0^2 + \varepsilon |u|_1^2 \f$,
  !> ityp = 4  \f$ \Rightarrow \| v \|_X^2 = \varepsilon ( |u|_1^2 + J_h^{\sigma}(u,u) )\f$,
  subroutine AnisotElemEstimate(elem, ityp)
    type(element) :: elem
    integer, intent(in) :: ityp
    type(volume_rule), pointer :: V_rule
    real, dimension(:,:), allocatable :: S, D, alpha, beta, Dpsi, DpsiT
    real, dimension(:), allocatable :: ident, grad, gradT
    real, dimension(:,:), allocatable :: Drw, Fx
    real, dimension(:,:,:), allocatable :: Dphi  ! derivatives of the test functions
    real, dimension(:,:,:), allocatable :: Dwi   ! gradient of the solution
    real :: elem_estim
    real :: est, estT, ratio, ratioM, diam, diamT, est0, est0T, lambda1
    real :: lambda_max, lambda_min, h_max, h_min, epsilon1, tol_loc, lam_actual
    integer :: i, dofA, dof1, dof, itest, j, k
    integer :: pK, Qnum, Qdof
    character(len=3):: metrictype ! SHX (X-norm), SHP(H1-norm), HES(Hessian metric)

    metrictype='HES'

    ! CH tol_loc snizovat postupne
    tol_loc = state%space%adapt%tol_min / grid%nelem**0.5
    epsilon1 = 1E+8
    h_max = state%space%diam / 2.
    h_min = h_max / epsilon1**0.5
    lambda_min = 1./h_max**2
    lambda_max = 1./h_min**2


    !write(*,'(a10,4es14.6)') 'h_min,max:', h_min, h_max, state%space%diam
    !write(*,'(a10,4es14.6)') 'l_min,max:', lambda_min, lambda_max

    elem%psplit = 0
    pK = max(1,elem%deg)
    elem%reg = elem%rezid / elem%area / elem%diam**(2*pK)* elem%diam**3

    elem%regT0 = elem%diam**(-2)
    elem%regT1 = elem%diam**(-3)
    elem%regT2 = elem%diam**(-4)

    ! quantity meassuring the (isotropic) size of the error [ VD, MATCOM ??]
    elem_estim = (elem%estim_loc + elem%jumpsJh)**0.5

    if(elem_estim > tol_loc ) then
       if(elem%reg <= elem%regT0 .and. elem%deg < MaxDegreeImplemented-1 ) then
          ! solution is regular => p-refinement
          elem%psplit = 1
!          write(41,*) elem%xc(:), elem%reg, elem_estim, tol_loc
       else
          ! solution isn't regular => h-refinement

          if(elem%reg >= elem%regT2 .and. elem%deg > MinDegreeImplemented) then
             elem%psplit = -1
!             write(51,*) elem%xc(:), elem%reg, elem_estim, tol_loc
          else
!             write(61,*) elem%xc(:), elem%reg, elem_estim, tol_loc
          endif
       endif

    elseif(elem_estim < tol_loc ) then

       if(elem%reg >= elem%regT2 .and. elem%deg > MinDegreeImplemented) then
          elem%psplit = -1
!          write(71,*) elem%xc(:), elem%reg, elem_estim, tol_loc
       else
!          write(81,*) elem%xc(:), elem%reg, elem_estim, tol_loc
       endif
    else

!       write(91,*) elem%xc(:), elem%reg, elem_estim, tol_loc

    endif

    !!elem%psplit = 0


!    write(*,'(a8,2i5,12es12.4)') 'regul:',elem%i, elem%psplit, elem%reg,&
!         elem_estim- tol_loc, elem%regT0, elem%regT2

    !stop

    ! the following is not used up to A3
    dof = elem%dof
    dofA = elem%dof_plus
    dof1 = 1

    Qnum = elem%Qnum
    Qdof = elem%Qdof

    allocate (Dwi(1:Qdof, 1:ndim, 1:nbDim),  Drw(1:Qdof, 1:2) )
    ! Drw(1:Qdof, 1)  ... derivative in the direction of the orientation of the element
    ! Drw(1:Qdof, 2)  ... derivative in the direction perpendicular to the orientation of the element

    ! evaluation od the "Stiff" matrix in the appropriate  norm
    ! S(i,j) = (( \phi_i, \phi_j)),  ((. , . )) generate the norm
    allocate(S(1:dofA, 1:dofA), D(1:dofA, 1:dofA) )
    S(:,:) = 0.   ! scalar product in X-norm
    D(:,:) = 0.   ! H^1 scalar product (stiff matrix)

    allocate(ident(1:elem%Qdof) )
    ident(:) = 1.
    call IntegrateBlockD2(elem, dofA, ident(:), D(1:dofA, 1:dofA) )
    deallocate(ident)

    if(ityp == 4) then
       !print*,' penalty'
       call IntegrateBlockJumps(elem, dofA, S(1:dofA, 1:dofA) )
    endif

    if(ityp >= 2) then  !print*,'semi-H_1 norm'
       allocate(ident(1:elem%Qdof) )
       ident(:) = 1.
       call IntegrateBlockD2(elem, dofA, ident(:), S(1:dofA, 1:dofA) )
       ! print*,'
       deallocate(ident)
    endif

    if(ityp == 3) then  !! "normalization: |u|_1 -> \sqrt{\varepsilon}| u |_1
       S(1:dofA, 1:dofA) =  S(1:dofA, 1:dofA) * state%model%Re1

       ! adding of the L2 norm
       do i=dof1,dofA
          S(i,i) = S(i,i) + 2*elem%area !*(2**0.5) !! sqrt(2) is the "size of convection"
       enddo
    endif

    if(ityp == 1 .or. ityp == 2 ) then  !! L_2 norm
       do i=dof1,dofA
          S(i,i) = S(i,i) + 2*elem%area
       enddo
    endif

    if(ityp == 4) then  !! "normalization: |u|_1 -> \sqrt{\varepsilon}| u |_1
       if(state%model%Re > 0.) S(1:dofA, 1:dofA) =  S(1:dofA, 1:dofA) / state%model%Re
    endif

    if(ityp == 5) then  !!???????? norm
       do i=dof1,dofA
          S(i,i) = 1.
       enddo
    endif


!    ! evaluation of the corresponding maxima over a discrete sets
    allocate(alpha(1:ndim, dof1:dofA) )
    allocate(beta(1:ndim, dof1:dofA) )


    ! "space discretization error", alpha contains the coefficients maximazing error estims A1 begin
    call EvalMaximumLangAMA(dofA, dof1, dofA, S(1:dofA, 1:dofA), &
         elem%vec(rhs, 1:dofA*ndim), est, alpha(1:ndim, dof1:dofA) )

    !A1 end
    if (metrictype=='SHX') then
       ! "space discretization orthogonal error", beta contains the coefficients maximazing error estims
       call EvalMaximumLangAMABeta(dofA, dof1, dofA, S(1:dofA, 1:dofA), &
            elem%vec(rhs, 1:dofA*ndim), estT, alpha, beta)
    endif

    !A2 begin
    if (metrictype=='SHP') then
       ! "space discretization orthogonal error", beta contains the coefficients maximazing error estims
       call EvalMaximumLangAMABeta2(dofA, dof1, dofA, S(1:dofA, 1:dofA), D(1:dofA, 1:dofA), &
            elem%vec(rhs, 1:dofA*ndim), estT, alpha(1:ndim, dof1:dofA), beta(1:ndim, dof1:dofA) )
    endif
    !A2 end

    ! gradients of the shape functions
    allocate(Dphi(1:dofA, 1:nbDim, 1:Qdof) )
    call Eval_Dphi(elem, dofA,  Dphi(1:dofA, 1:nbDim, 1:Qdof) )

    allocate(Dpsi(1:Qdof, 1:nbDim) )  ! Dpsi = gradient of the maximazing function from Shp
    allocate(grad(1:nbDim) )          ! averaging of Dpsi
    allocate(DpsiT(1:Qdof, 1:nbDim) )  ! DpsiT = gradient of the maximazing function from Shp \ psi
    allocate(gradT(1:nbDim) )          ! averaging of DpsiT

    V_rule => state%space%V_rule(Qnum)

    call  Eval_Dw_Elem(elem, Dwi(1:Qdof, 1:ndim, 1:nbDim) )

    !do k=1,ndim   ! FOR ALL COMPONENT CHANGES OF  grad(*) NECESSARY !!!!!
    do k=1,1
       ! A1 approach based on the seeking of constrained extrama in orthogonal subspaces of S_hp^+
       Dpsi(1:Qdof, 1) = matmul(alpha(k, 1:dofA) , Dphi(1:dofA, 1, 1:Qdof) )
       Dpsi(1:Qdof, 2) = matmul(alpha(k, 1:dofA) , Dphi(1:dofA, 2, 1:Qdof) )

       DpsiT(1:Qdof, 1) = matmul(beta(k, 1:dofA) , Dphi(1:dofA, 1, 1:Qdof) )
       DpsiT(1:Qdof, 2) = matmul(beta(k, 1:dofA) , Dphi(1:dofA, 2, 1:Qdof) )


       ! A2 alternative approach base on the directional derivatives of u_h
       !Dpsi(1:Qdof, 1:2) =  Dwi(1:Qdof, k, 1:2 )

       do i=1,Qdof
          if(Dpsi(i,1) < 0)  Dpsi(i,1:2) = -  Dpsi(i,1:2)
          if(DpsiT(i,1) < 0) DpsiT(i,1:2) = - DpsiT(i,1:2)
       enddo


       grad(1) = dot_product(V_rule%weights(1:Qdof),  Dpsi(1:Qdof, 1) )
       grad(2) = dot_product(V_rule%weights(1:Qdof),  Dpsi(1:Qdof, 2) )

       gradT(1) = dot_product(V_rule%weights(1:Qdof),  DpsiT(1:Qdof, 1) )
       gradT(2) = dot_product(V_rule%weights(1:Qdof),  DpsiT(1:Qdof, 2) )

       ! grad(1:2) is the direction perpendicular to the orientation of the element
       grad(1:2) = grad(1:2)/( dot_product(grad(:), grad(:)))**0.5 !normalization, grad(1)=cos, grad(2)=sin
       gradT(1:2) = gradT(1:2)/( dot_product(gradT(:), gradT(:)))**0.5 !normalization

       ! A2 alternative approach base on the directional derivatives of u_h
       !Drw(1:Qdof,1) = Dwi(1:Qdof, k, 1)* grad(2) - Dwi(1:Qdof, k, 2)* grad(1)
       !Drw(1:Qdof,2) = Dwi(1:Qdof, k, 1)* grad(1) + Dwi(1:Qdof, k, 2)* grad(2)

    enddo


    ! A1 approach based on the seeking of constrained extrama in orthogonal subspaces of S_hp^+
    ! inverse of the aspect ratio of the new element
    if (metrictype=='SHP') then
       ratio = (estT / est)**1.0
    endif

    if (metrictype=='SHX') then
       ratio = (estT / est)**1.0
    endif

    ! A2  alternative approach base on the directional derivatives of u_h
    ! inverse of the ratio = ratio of the direction derivative of w_h in the orientation of the element
    ! and the perpendicular direction
    !ratio = abs ( dot_product(V_rule%weights(1:Qdof),  Drw(1:Qdof, 1) ) / &
    !     dot_product(V_rule%weights(1:Qdof),  Drw(1:Qdof, 2) )  )

    !write(*,'(a6,10es14.6)') 'AAA', dot_product(V_rule%weights(1:Qdof),  Drw(1:Qdof, 1) ) , &
    !     dot_product(V_rule%weights(1:Qdof),  Drw(1:Qdof, 2) ), ratio


    ! the previous is not used, old versions
    ! A3 approach partly based on the Hessian matrix
    ! the following defines metric of the actual (old) grid, element size est, element ratio est0 /est0T


    diam = elem%diam
    diamT = 2*elem%area / diam
    est0 = 1.5 / diamT**2
    est0T = 2*3**0.5 / elem%diam**2

    lam_actual = (2*elem%area/elem%diam)**(-2)

    ! approach based on the Hessian matrix only for A3
    !A1
    if (metrictype=='HES') then
       call EvalRatioAngle(elem, lambda1, ratio, grad(1), grad(2) )
       !all
    endif

    !ch jine 0.001, 0.0001
    ratio = max(ratio, 0.001)  ! maximal aspect ration (=1./ratio) !! (should be positivity)
    !ratio = max(ratio, diamT/ diam *0.2)   ! ratio may increase maximally 5-times

    ! ryzi hodnoty est, estT, grad(1,2) k dispozici ????

    ! est = lambda_1 - optimal value
    est = (1./diamT *( elem_estim / tol_loc)**(1./(elem%deg+elem%psplit)))**2

    ! maximal enlargement / de-enlargement
    est  = max(2.E-02 * lam_actual , min(50.E+00 * lam_actual,  est) )
    !est  = max(2.E-01 * est0 , min(5.E+00 * est0,  est) )
    !estT = max(2.E-01 * est0T, min(5.E+00 * est0T, est) )

    estT = est  * ratio

    !  keeping the maximal and the minimal eigenvalues == minimal and maximal edge
    est  = max(lambda_min, min(lambda_max, est) )
    estT = max(lambda_min, min(lambda_max, estT) )

    elem%amaL = est
    elem%amaR = estT / est
    elem%amaD(1:2) = grad(1:2)

    elem%rgabc(1) = est * grad(1)**2 + estT * grad(2)**2
    elem%rgabc(2) = (est - estT) * grad(1) * grad(2)
    elem%rgabc(3) = est * grad(2)**2 + estT * grad(1)**2

    !write(*,'(i5, 3es14.6, a2, 4(2es14.6, a2), 10es14.6)') &
    !     elem%i,  elem%rgabc(1), elem%rgabc(2), elem%rgabc(3), '|',est, estT, '|', &
    !     est0, est0T, '|', ratio, ratio, '|', grad(1:2), '|', elem%xc(:)


    deallocate (S, D, alpha, beta, Dwi, Drw, Dphi, Dpsi)
    deallocate (grad)


  end subroutine AnisotElemEstimate

  !> evaluate the "optimal" aspect ratio and the triangle orientation form the Hessian matrix
  !> stored in elem%rgabc(1), elem%rgabc(2)(off-diagonal), elem%rgabc(3)
  subroutine EvalRatioAngle(elem, lambda, ratio, cosalpha, sinalpha)
    type(element), intent(inout) :: elem
    real, intent(out):: lambda, ratio, cosalpha, sinalpha
    external :: DSTEQR
    CHARACTER ::  COMPZ
    INTEGER   ::  INFO, LDZ, N
    REAL, dimension(:),  allocatable:: D, E, WORK
    REAL, dimension(:,:),  allocatable:: Z
    REAL, dimension(:,:),  allocatable:: M,DD,R,R1,S
    integer :: i

    COMPZ = 'I'    ! only for tridiagonal matrix
    N = 2
    LDZ = N
    allocate( D(1:N), E(1:N-1), Z(1:LDZ, 1:N), WORK(1: max(1,2*N-2)) )

    allocate(M(1:N, 1:N), DD(1:N, 1:N), R(1:N, 1:N), R1(1:N, 1:N), S(1:N, 1:N) )

    D(1) = elem%rgabc(1)
    D(2) = elem%rgabc(3)
    E(1) = elem%rgabc(2)


    call DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
    if(INFO /= 0) then
       print*,'Trouble in EvalRatioAngle in calling the LAPACK subroutine DSTEQR, file ama-estims.f90'
       stop
    endif

    DD(:,:) = 0.
    DD(1,1) = D(1)
    DD(2,2) = D(2)
    M(1,1) = elem%rgabc(1)
    M(1,2) = elem%rgabc(2)
    M(2,1) = elem%rgabc(2)
    M(2,2) = elem%rgabc(3)

    R(:,:) = Z(1:2, 1:2)
    R1(:,:) = transpose(R(:,:) )
    S(:,:)  = matmul(R(:,:), matmul(DD(:,:), R1(:,:) ) )

    !do i=1, N
    !   write(*,'(4(a3, 2es12.4))') 'M',M(i,:), 'S',S(i, :), 'R',R(i, :)
    !enddo
    !print*,'------------------------'
    !write(*,'(4(a3, 2es12.4))') 'Mv', matmul(M(:,:), Z(1:2,1)), 'lv',D(1)*Z(1:2, 1)
    !write(*,'(4(a3, 2es12.4))') 'Mv', matmul(M(:,:), Z(1:2,2)), 'lv',D(2)*Z(1:2, 2)

    !print*,'INFO=', INFO
    !write(*,'(4(a3, 2es12.4))') 'M',elem%rgabc(1), elem%rgabc(2),'D',D(1), 0., 'V',Z(1, 1), Z(1,2)
    !write(*,'(4(a3, 2es12.4))') 'M',elem%rgabc(2), elem%rgabc(3),'D',0., D(2), 'V',Z(2, 1), Z(2,2)

    DD(:,:) = abs(DD(:,:) )
    S(:,:)  = matmul(R(:,:), matmul(DD(:,:), R1(:,:) ) )
    elem%rgabc(1) = S(1,1)
    elem%rgabc(2) = S(1,2)
    elem%rgabc(3) = S(2,2)

    deallocate(M,DD,R,R1,S)

    lambda = max( abs(D(1)), abs(D(2)) )
    ratio =  min( abs(D(1)), abs(D(2)) ) / lambda
    cosalpha = Z(1, 1)
    sinalpha = Z(2, 1)

    if(abs(D(1)) <  abs(D(2)) ) then
       cosalpha = Z(1, 2)
       sinalpha = Z(2, 2)
    endif

    !write(*,'(a6,i5,20es12.4)') 'Ani:',elem%i, lambda, ratio, cosalpha, sinalpha, abs(D(1))/ abs(D(2))

    deallocate(D, E, WORK, Z)
  end subroutine EvalRatioAngle


  ! evaluate the appropriate maximum over a discrete set
  subroutine EvalMaximumLangAMA(dof_tot, dof_L, dof_U, S, vec, estim, alpha)
    integer, intent(in) :: dof_tot, dof_L, dof_U ! dof: total, lowe, upper
    real, dimension(1:dof_tot, 1:dof_tot), intent(in) :: S
    real, dimension(1:dof_tot*ndim), intent(in) :: vec
    real, intent(inout) :: estim
    real, dimension(1:ndim, dof_L:dof_U), intent(out) :: alpha
    real, dimension(:, :), allocatable :: Sinv
    real :: lambda
    integer :: i, k, k1, k2

    if(dof_L > dof_tot .or. dof_U < dof_L ) then
       print*,'bad arrays bounds in project.f90:EvalMaximumLang',dof_tot, dof_L, dof_U
       stop
    endif

    allocate(Sinv(dof_L:dof_U, dof_L:dof_U) )

    ! "enrichment" of the diagonal following from the Lagrange multiplier
    Sinv(dof_L:dof_U, dof_L:dof_U ) = S(dof_L:dof_U, dof_L:dof_U )
    do i = dof_L, dof_U
       Sinv(i,i) = 2*Sinv(i,i)
    enddo

    call MblockInverse(dof_U - dof_L +1 , Sinv(dof_L:dof_U, dof_L:dof_U ))

    alpha(:,:) = 0.

    estim = 0.
    do k=1,ndim
       k1 = (k-1) * dof_tot + dof_L
       k2 = (k-1) * dof_tot + dof_U

       alpha(k, dof_L:dof_U) = matmul(Sinv(dof_L:dof_U, dof_L:dof_U), vec(k1:k2) )
       !write(*,'(a6,12es12.4)') 'alpha:',alpha(k, dof_L:dof_U)
       lambda = ( dot_product(alpha(k, dof_L:dof_U), &
            matmul(S(dof_L:dof_U,dof_L:dof_U), alpha(k, dof_L:dof_U) )) )**0.5

       !!!write(197,*) '///',k,k1,k2,dof_L,dof_U, lambda

       if(lambda > 0) then
          alpha(k, dof_L:dof_U) = alpha(k, dof_L:dof_U) / lambda
          estim = estim + abs(dot_product(vec(k1:k2), alpha(k, dof_L:dof_U)) )

       !else
          ! alpha = 0., estim = 0.
       endif
    enddo

    deallocate(Sinv)


  end subroutine EvalMaximumLangAMA



subroutine EvalMaximumLangAMABeta(dof_tot, dof_L, dof_U, S, vec, estim, alpha, beta)
    integer, intent(in) :: dof_tot, dof_L, dof_U ! dof: total, lowe, upper
    real, dimension(dof_L:dof_U, dof_L:dof_U), intent(in) :: S
    real, dimension(1:dof_tot*ndim), intent(in) :: vec
    real, intent(inout) :: estim !, estimSBeta
    real, dimension(1:ndim, dof_L:dof_U), intent(in) :: alpha
    real, dimension(dof_L:dof_U, 1:ndim):: alphaT ! transposed matrix
    real, dimension(1:ndim, dof_L:dof_U), intent(inout):: beta
    real, dimension(:, :), allocatable :: Sinv
    real, dimension(1:ndim,dof_L:dof_U) ::A, B, D
    real, dimension(dof_L:dof_U, 1:ndim) ::C
    real :: lambda, mu !,estimSBeta
    integer :: i, k, k1, k2, j

    if(dof_L > dof_tot .or. dof_U < dof_L ) then
       print*,'bad arrays bounds in project.f90:EvalMaximumLang',dof_tot, dof_L, dof_U
       stop
    endif

    allocate(Sinv(dof_L:dof_U, dof_L:dof_U) )

    ! "enrichment" of the diagonal following from the Lagrange multiplier
    Sinv(dof_L:dof_U, dof_L:dof_U ) = S(dof_L:dof_U, dof_L:dof_U )
    do i = dof_L, dof_U
       Sinv(i,i) = 2*Sinv(i,i)
    enddo

    call MblockInverse(dof_U - dof_L +1 , Sinv(dof_L:dof_U, dof_L:dof_U ))

    do i=dof_L, dof_U
       do j=1,ndim
          alphaT(i,j)=alpha(j,i)
       enddo
    enddo


    beta(:,:)=0
    C(:,:)=0
    B(:,:)=0

    estim = 0.

    do k=1,ndim
       k1 = (k-1) * dof_tot + dof_L
       k2 = (k-1) * dof_tot + dof_U
       mu=(- dot_product(alphaT( dof_L:dof_U, k), &
            (matmul(S(dof_L:dof_U, dof_L:dof_U),matmul(Sinv(dof_L:dof_U, dof_L:dof_U), vec(k1:k2))) )))/ &
            (dot_product(alphaT(dof_L:dof_U, k), matmul(matmul(S(dof_L:dof_U, dof_L:dof_U), Sinv(dof_L:dof_U, dof_L:dof_U)),&
            matmul(S(dof_L:dof_U, dof_L:dof_U), alpha(k, dof_L:dof_U)))))

       A(k, dof_L:dof_U)=matmul(alpha(k,dof_L:dof_U), S(dof_L:dof_U, dof_L:dof_U))
       B(k, dof_L:dof_U)=vec(k1:k2) +mu *A(k, dof_L:dof_U)
       D(k, dof_L:dof_U)=matmul(B(k, dof_L:dof_U), Sinv(dof_L:dof_U, dof_L:dof_U))
       do i=dof_L, dof_U
          do j=1, ndim
             C(i,j)=D(j,i)
          enddo
       enddo
       lambda=(dot_product(C(dof_L:dof_U, k), matmul(S(dof_L:dof_U, dof_L:dof_U), D(k, dof_L:dof_U))))**0.5
       beta(k, dof_L:dof_U)=-1/lambda*matmul(Sinv(dof_L:dof_U, dof_L:dof_U),B(k, dof_L:dof_U))
       estim = estim + abs(dot_product(vec(k1:k2), beta(k, dof_L:dof_U)) )
    enddo

  end subroutine EvalMaximumLangAMABeta


  !> evaluate \f$ \max_{\psi'\in S_{hp}^+, \nabla\psi'\cdot\nabla\psi =0} c_h(w_h,\psi')\f$,
  !> where \f$\psi\f$ is the function maximazing
  !> \f$ \max_{\psi\in S_{hp}^+} c_h(w_h,\psi)\f$ given by vector alpha
subroutine EvalMaximumLangAMABeta2(dof_tot, dof_L, dof_U, S, D, vec, estim, alpha, beta)
    integer, intent(in) :: dof_tot, dof_L, dof_U ! dof: total, lowe, upper
    real, dimension(1:dof_tot, 1:dof_tot), intent(in) :: S  ! X-scalar product
    real, dimension(1:dof_tot, 1:dof_tot), intent(in) :: D  ! H1-scalar product
    real, dimension(1:dof_tot*ndim), intent(in) :: vec
    real, intent(inout) :: estim !, estimSBeta
    real, dimension(1:ndim, dof_L:dof_U), intent(in) :: alpha
    !real, dimension(dof_L:dof_U, 1:ndim):: alphaT ! transposed matrix
    real, dimension(1:ndim, dof_L:dof_U), intent(inout):: beta
    real, dimension(:, :), allocatable :: Sinv
    !real, dimension(1:ndim,dof_L:dof_U) ::A, B, D
    real, dimension(:), allocatable :: Da, Sd, SDa
    real :: lam1, lam2
    integer :: i, k, k1, k2, j

    if(dof_L > dof_tot .or. dof_U < dof_L ) then
       print*,'bad arrays bounds in project.f90:EvalMaximumLang',dof_tot, dof_L, dof_U
       stop
    endif

    allocate(Da(dof_L:dof_U),  SDa(dof_L:dof_U), Sd(dof_L:dof_U) )
    allocate(Sinv(dof_L:dof_U, dof_L:dof_U) )

    ! "enrichment" of the diagonal following from the Lagrange multiplier
    Sinv(dof_L:dof_U, dof_L:dof_U ) = S(dof_L:dof_U, dof_L:dof_U )
    do i = dof_L, dof_U
       Sinv(i,i) = 2*Sinv(i,i)
    enddo

    call MblockInverse(dof_U - dof_L +1 , Sinv(dof_L:dof_U, dof_L:dof_U ))

    !do i=dof_L, dof_U
    !   do j=1,ndim
    !      alphaT(i,j)=alpha(j,i)
    !   enddo
    !enddo


    beta(:,:)=0

    estim = 0.

    do k=1,ndim
       k1 = (k-1) * dof_tot + dof_L
       k2 = (k-1) * dof_tot + dof_U
       Da(dof_L:dof_U)  = matmul(D(dof_L:dof_U, dof_L:dof_U), alpha(k, dof_L:dof_U) )
       Sd(dof_L:dof_U)  = matmul(Sinv(dof_L:dof_U, dof_L:dof_U), vec(k1:k2) )
       SDa(dof_L:dof_U) = matmul(Sinv(dof_L:dof_U, dof_L:dof_U), Da(dof_L:dof_U) )

       lam2 = - dot_product(  Sd(dof_L:dof_U), Da(dof_L:dof_U) ) &
            / dot_product( SDa(dof_L:dof_U), Da(dof_L:dof_U) )

       Sd(dof_L:dof_U) = -vec(k1:k2) - lam2* Da(dof_L:dof_U)
       Sda(dof_L:dof_U) =  matmul(Sinv(dof_L:dof_U, dof_L:dof_U), Sd(dof_L:dof_U) )

       lam1 = dot_product(Sda(dof_L:dof_U), matmul( S(dof_L:dof_U, dof_L:dof_U ), Sda(dof_L:dof_U) ))

       beta(k, dof_L:dof_U) = Sda(dof_L:dof_U) / lam1**0.5

       estim = estim + abs( dot_product(vec(k1:k2), beta(k, dof_L:dof_U)) )
    enddo

    deallocate(Da, Sd, SDa, Sinv)
  end subroutine EvalMaximumLangAMABeta2



  subroutine AnisotErrorEstimatesBeta( )
    class(element), pointer :: elem, elem1
    real, dimension(:), allocatable :: L_estim
    real, dimension(:,:), allocatable :: wi
    real, dimension(:,:,:), allocatable :: Dwi
    real :: machine_tol, rmax
    integer :: i, j, k, ndof, ndofP, itest, imax
    logical :: loc_implicitly

    !itest = 360
    itest = -480

    allocate(L_estim(1:max_eta) )
    loc_implicitly = state%nlSolver%implicitly

    state%nlSolver%implicitly = .false.
    grid%elem(:)%deg_plus = .true.

    ! setting of fields elem%vec(rhs,*), elem%vec(rhsT,*) for error estimate
    call SetVectorsFields(.false.)

    state%estim(max_resT_S,:) = -1E+50
    state%estim(min_resT_S,:) =  1E+50
    state%estim(min_resT_S_loc,:) =  1E+50
    state%estim(resA_ST,:) = 1E+50

    rmax = 0.

    L_estim(:) = 0.   ! total value of the residuum
    do i=1,grid%nelem
       elem => grid%elem(i)
       ! NOT MULTIPLIED,1/tau in project.f90 removed
       !elem%vec(rhs,:) = elem%vec(rhs,:) * state%time%tau(1)

       ! the following should be performed for the STDGM approach
       !elem%vec(rhsT,:) = elem%vec(rhsT,:) / state%time%tau(1)

       !!!call EnergyReziduumElemEstimate(elem)  ! element residuum
       !do j=1,4
       call AnisotElemEstimateBeta(elem, 3)  ! 3 => element residuum in X- norm
       L_estim( resA) = L_estim( resA) + elem%estimA**2
       L_estim( resS) = L_estim( resS) + elem%estimS**2
       L_estim( resT) = L_estim( resT) + elem%estimT**2
       L_estim(resST) = L_estim(resST) + elem%estimST**2


       ! Verfurth approach
       state%estim(min_resT_S_loc,1) = min(state%estim(min_resT_S_loc,1),  &
            elem%estimT /  elem%estimS )


    enddo
    print*,'CODE STOPPED HERE in ama-estims.f90'

    stop

    state%L_estim(1:max_eta) = L_estim(max_eta)**0.5

    !write(*,*) '####', imax, rmax

    machine_tol = 1.E-01
    state%estim(resA_ST_loc,1) = 0.
    ! local algebraic criterion
    do i=1,grid%nelem
       elem => grid%elem(i)

       if( elem%estimST >  machine_tol * state%L_estim(resST) / grid%nelem**0.5 .and. &
            elem%estimA/ elem%estimST >  state%estim(resA_ST_loc,1)) then

          state%estim(resA_ST_loc,1) = max(state%estim(resA_ST_loc,1), elem%estimA/ elem%estimST)

       endif
    enddo

    ! global algebraic criterion
    ! STDG approach
    !state%estim(resA_ST,1) = state%L_estim(resA)/state%L_estim(resST)

    ! Verfurth approach
    !state%estim(resA_ST,1) = state%L_estim(resA) / ( state%L_estim(resS) + state%L_estim(resT) )

    ! steady-state approach
    state%estim(resA_ST,1) = state%L_estim(resA) /  state%L_estim(resS)

    ! Verfurth approach
    !state%estim(min_resT_S, 1) = state%L_estim(resST)/state%L_estim(resS)

    ! STDG approach
    state%estim(min_resT_S, 1) = state%L_estim(resT)/state%L_estim(resS)

    grid%elem(:)%deg_plus = .false.

    ! setting of 'elem%estim_loc' including elem%estimST of neighbours elements
    ! elem%estim_loc**2 = sum_{K\in N(K)} elem%estimST**2
    do i=1,grid%nelem
       elem => grid%elem(i)
       !elem%estim_loc = elem%estimST**2    ! STDGM approach
       elem%estim_loc = elem%estimS**2    ! Verfurth approach
       do j=1,elem%flen
          k = elem%face(neigh,j)
          if(k > 0) then
             elem1 => grid%elem(k)
             !elem%estim_loc = elem%estim_loc + elem1%estimST**2   ! STDGM approach
             elem%estim_loc = elem%estim_loc + elem1%estimS**2   ! Verfurth approach
             !elem%estim_loc = max(elem%estim_loc , elem1%estimS**2)   ! Verfurth approach
          endif
       enddo
       !!!  elem%estim_loc = elem%estim_loc**0.5  ! we store the square
    enddo

    state%nlSolver%implicitly  = loc_implicitly
    deallocate(L_estim)

    !if( state%time%cn ) then
    !   do i=1,grid%nelem
    !      elem => grid%elem(i)
    !   enddo
    !endif

  end subroutine AnisotErrorEstimatesBeta


!>subroutine to compute special metrics
 subroutine AnisotElemEstimateBeta(elem, ityp)
    type(element) :: elem
    integer, intent(in) :: ityp
    type(volume_rule), pointer :: V_rule
    real, dimension(:,:), allocatable :: S, alpha, Dpsi, beta, DpsiT
    real, dimension(:), allocatable :: ident, grad, gradT
    real, dimension(:,:,:), allocatable :: Dphi
    integer :: i, dofA, dof1, dof, itest, j, k
    integer :: Qnum, Qdof

    dof = elem%dof
    dofA = elem%dof_plus
    dof1 = 1

    Qnum = elem%Qnum
    Qdof = elem%Qdof

    ! evaluation od the "Stiff" matrix in the appropriate  norm
    ! S(i,j) = (( \phi_i, \phi_j)),  ((. , . )) generate the norm
    allocate(S(1:dofA, 1:dofA) )
    S(:,:) = 0.

    if(ityp == 4) then
       !print*,' penalty'
       call IntegrateBlockJumps(elem, dofA, S(1:dofA, 1:dofA) )
    endif

    if(ityp >= 2) then  !print*,'semi-H_1 norm'
       allocate(ident(1:elem%Qdof) )
       ident(:) = 1.
       call IntegrateBlockD2(elem, dofA, ident(:), S(1:dofA, 1:dofA) )
       ! print*,'
       deallocate(ident)
    endif

    if(ityp == 3) then  !! "normalization: |u|_1 -> \sqrt{\varepsilon}| u |_1
       S(1:dofA, 1:dofA) =  S(1:dofA, 1:dofA) * state%model%Re1

       ! adding of the L2 norm
       do i=dof1,dofA
          S(i,i) = S(i,i) + 2*elem%area !*(2**0.5) !! sqrt(2) is the "size of convection"
       enddo
    endif

    if(ityp == 1 .or. ityp == 2 ) then  !! L_2 norm
       do i=dof1,dofA
          S(i,i) = S(i,i) + 2*elem%area
       enddo
    endif

    if(ityp == 4) then  !! "normalization: |u|_1 -> \sqrt{\varepsilon}| u |_1
       if(state%model%Re > 0.) S(1:dofA, 1:dofA) =  S(1:dofA, 1:dofA) / state%model%Re
    endif

    if(ityp == 5) then  !!???????? norm
       do i=dof1,dofA
          S(i,i) = 1.
       enddo
    endif


    Qdof = elem%Qdof

    allocate(alpha(1:ndim, dof1:dofA) )
    allocate(beta(1:ndim, dof1:dofA) )

    call EvalMaximumLangAMA(dofA, dof1, dofA, S(1:dofA, 1:dofA), &
         elem%vec(rhs, 1:dofA*ndim), elem%estimS, alpha)

    ! "space discretization orthogonal error", beta contains the coefficients maximazing error estims
    call EvalMaximumLangAMABeta(dofA, dof1, dofA, S(1:dofA, 1:dofA), &
         elem%vec(rhs, 1:dofA*ndim), elem%estimS, alpha, beta)

    ! gradients of the shape functions
    allocate(Dphi(1:dofA, 1:nbDim, 1:Qdof) )
    call Eval_Dphi(elem, dofA,  Dphi(1:dofA, 1:nbDim, 1:Qdof) )

    allocate(Dpsi(1:Qdof, 1:nbDim) )  ! Dpsi = gradient of the maximazin function
    allocate(grad(1:nbDim) )          ! averaging of Dpsi
    V_rule => state%space%V_rule(Qnum)

    do k=1,ndim
       Dpsi(1:Qdof, 1) = matmul(alpha(k, 1:dofA) , Dphi(1:dofA, 1, 1:Qdof) )
       Dpsi(1:Qdof, 2) = matmul(alpha(k, 1:dofA) , Dphi(1:dofA, 2, 1:Qdof) )

       grad(1) = dot_product(V_rule%weights(1:Qdof),  Dpsi(1:Qdof, 1) )
       grad(2) = dot_product(V_rule%weights(1:Qdof),  Dpsi(1:Qdof, 2) )

       write(110,*) elem%xc(1:2),elem%i

       write(120,*) elem%xc(1:2)
       write(120,*) elem%xc(1:2) + grad(1:2)/( dot_product(grad(:), grad(:) )**0.5)*elem%diam/4
       write(120,'(x)')
       write(130,*) elem%xc(1:2), elem%estimS
    enddo


    ! gradients of the 'orthogonal function'
    allocate(DpsiT(1:Qdof, 1:nbDim) )  ! DpsiT = gradient of the maximazin function
    allocate(gradT(1:nbDim) )          ! averaging of DpsiT
    V_rule => state%space%V_rule(Qnum)

    do k=1,ndim
       DpsiT(1:Qdof, 1) = matmul(beta(k, 1:dofA) , Dphi(1:dofA, 1, 1:Qdof) )
       DpsiT(1:Qdof, 2) = matmul(beta(k, 1:dofA) , Dphi(1:dofA, 2, 1:Qdof) )

       gradT(1) = dot_product(V_rule%weights(1:Qdof),  Dpsi(1:Qdof, 1) )
       gradT(2) = dot_product(V_rule%weights(1:Qdof),  Dpsi(1:Qdof, 2) )

       write(110,*) elem%xc(1:2),elem%i

       write(120,*) elem%xc(1:2)
       write(120,*) elem%xc(1:2) + gradT(1:2)/( dot_product(gradT(:), gradT(:) )**0.5)*elem%diam/4
       write(120,'(x)')
       write(130,*) elem%xc(1:2), elem%estimS
    enddo


    elem%rgabc(1) = elem%estimS / (state%space%adapt%tol_min )**0.5
    elem%rgabc(2) = 0.
    elem%rgabc(3) = elem%estimS / (state%space%adapt%tol_min )**0.5

	!elem%rgabc(1) = elem%estimS / (state%space%adapt%tol_min )**0.5
    !elem%rgabc(2) = 0.
    !elem%rgabc(3) = elem%estimS / (state%space%adapt%tol_min )**0.5



    deallocate(S, alpha, Dphi, Dpsi, grad, DpsiT, gradT)

  end subroutine AnisotElemEstimateBeta


  !> evaluate Hassian matrix for each element, stored in elem%rgabc(1), elem%rgabc(2)(off-diagonal), elem%rgabc(3)
  subroutine EvalHessians ( )
    class(element), pointer :: elem
    real, dimension(:), allocatable :: wp, we, supp
    real, dimension(:,:), allocatable :: wi, dw, d2w, xp, Dwi
    character*5 quantity
    integer :: i,j,k, j1, k1
    integer:: itest

    if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
       quantity = 'RO'
    else
       quantity = 'RO'  !'M'
    endif

    allocate(wi(1:1, 1:ndim), we(1: grid%nelem) , wp(1: grid%npoin), supp(1: grid%npoin) )

    wp(:) = 0.
    supp(:) = 0.

    do i=1,grid%nelem
       elem => grid%elem(i)
       call Eval_aver_w_Elem(elem, wi(1, 1:ndim) )
       call ComputeQuantity(1, ndim, wi(1, 1:ndim), we(i), quantity)

       do j=1,elem%flen
          k = elem%face(idx, j)
          wp(k) = wp(k) + we(i) * elem%area
          supp(k) = supp(k) + elem%area
       enddo

       !write(88,*) elem%xc, we(i)
    enddo

    do i=1,grid%npoin
       wp(i) = wp(i) / supp(i)
       !write(89,*) grid%x(i,:), wp(i)
    enddo

    allocate(dw(1:grid%nelem, 1:2), xp(1:3, 1:3) )
    do i=1,grid%nelem
       elem => grid%elem(i)
       xp(1:3, 1:2) = grid%x(elem%face(idx, 1:3), 1:2)
       xp(1:3, 3) = wp(elem%face(idx, 1:3) )

       call EvalGradP1( xp(1:3, 1:3), dw(i, 1:2) )

       !write(87,*) xp(1, 1:3)
       !write(87,*) xp(2, 1:3)
       !write(87,*) xp(3, 1:3)
       !write(87,*) xp(1, 1:3)
       !write(87,* ) '  '
       !write(87,* ) '  '
       !write(87,* ) '  '
       !write(86,*) elem%xc(1:2), dw(i, 1:2)
    enddo

    allocate( Dwi(1:ndim,1:nbDim))
    do i=1,grid%nelem
       elem => grid%elem(i)
       if(elem%deg >= 1) then
          call Eval_aver_Dw_Elem(elem, Dwi(1:ndim, 1:nbDim) )
          dw(i,1:2)  = Dwi(1, 1:2)  !! ONLY DENSITY AT THIS MOMENT
       endif
       !write(85,*) elem%xc(1:2), dw(i, 1:2)

    enddo
    deallocate(Dwi)

    !itest = 2 ! 2

    allocate(d2w(1:grid%npoin, 1:5) )
    d2w(:,:) = 0.
    do i=1,grid%nelem
       elem => grid%elem(i)

       !if(i==16) write(*,'(a6,10i5)') '???',elem%i, i,  elem%face(idx, :),  elem%face(neigh, :)

       do j=1,elem%flen
          k = elem%face(idx, j)
          j1 = mod(j,  elem%flen) + 1

          d2w(k, 1) = d2w(k, 1) + dw(i, 1) * elem%n(j1, 1) /2 !/ elem%dn(j1)
          d2w(k, 2) = d2w(k, 2) + dw(i, 1) * elem%n(j1, 2) /2  !/ elem%dn(j1)
          d2w(k, 3) = d2w(k, 3) + dw(i, 2) * elem%n(j1, 1) /2  !/ elem%dn(j1)
          d2w(k, 4) = d2w(k, 4) + dw(i, 2) * elem%n(j1, 2) /2  !/ elem%dn(j1)
          d2w(k, 5) = d2w(k, 5) + elem%area /3.!!!!!!!!!!!!!* 4

          !if(k==itest) write(*,'(a1,4i5,20es12.4)') 'i',i, j, k, j1,  &
          !     dw(i, 1), elem%n(j1, 1),  dw(i, 1) * elem%n(j1, 1) ,  d2w(k, 1)

          if(elem%face(neigh, j) <= 0  )  then  ! adding of the boundary edge
             d2w(k, 1) = d2w(k, 1) + dw(i, 1) * elem%n(j, 1) /2  !/ elem%dn(j)
             d2w(k, 2) = d2w(k, 2) + dw(i, 1) * elem%n(j, 2) /2  !/ elem%dn(j)
             d2w(k, 3) = d2w(k, 3) + dw(i, 2) * elem%n(j, 1) /2  !/ elem%dn(j)
             d2w(k, 4) = d2w(k, 4) + dw(i, 2) * elem%n(j, 2) /2  !/ elem%dn(j)

             !if(k==itest) write(*,'(a1, 4i5,20es12.4)') 'b',i, j, k, elem%face(neigh, j),  &
             !     dw(i, 1), elem%n(j, 1), dw(i, 1) * elem%n(j, 1) ,  d2w(k, 1)

             k1 = elem%face(idx, j1)
             d2w(k1, 1) = d2w(k1, 1) + dw(i, 1) * elem%n(j, 1) /2  !/ elem%dn(j)
             d2w(k1, 2) = d2w(k1, 2) + dw(i, 1) * elem%n(j, 2) /2  !/ elem%dn(j)
             d2w(k1, 3) = d2w(k1, 3) + dw(i, 2) * elem%n(j, 1) /2  !/ elem%dn(j)
             d2w(k1, 4) = d2w(k1, 4) + dw(i, 2) * elem%n(j, 2) /2  !/ elem%dn(j)

             !if(k1==itest) write(*,'(a1, 4i5,20es12.4)') 'b',i, j, k1, elem%face(neigh, j),  &
             !     dw(i, 1), elem%n(j, 1), dw(i, 1) * elem%n(j, 1) ,  d2w(k1, 1)
          endif

       enddo
    enddo

    do i=1,grid%npoin

       !if(i==itest) write(*,'(a1, 4i5,20es12.4)') 'T',i, 0, 0,0,d2w(i, 1:4) ,d2w(i, 1:4) / d2w(i, 5)


       d2w(i, 1:4) = d2w(i, 1:4) / d2w(i, 5)
       !write(84,*) grid%x(i, 1:2), d2w(i, 1:5)
    enddo

    do i=1,grid%nelem
       elem => grid%elem(i)
       elem%rgabc(1) = 0.
       elem%rgabc(2) = 0.
       elem%rgabc(3) = 0.

       k1 = 0.
       do j=1,elem%flen
          k = elem%face(idx, j)
          if(d2w(k, 5) > 0.) then
             elem%rgabc(1) = elem%rgabc(1) + d2w(k, 1)
             elem%rgabc(2) = elem%rgabc(2) + d2w(k, 2)
             elem%rgabc(3) = elem%rgabc(3) + d2w(k, 4)
             k1 = k1 + 1
          endif
       enddo
       elem%rgabc(1) = elem%rgabc(1) / k1
       elem%rgabc(2) = elem%rgabc(2) / k1

       elem%rgabc(3) = elem%rgabc(3) / k1
       !write(83, *) elem%xc(:), elem%rgabc(1), elem%rgabc(2), elem%rgabc(3)
    enddo

    deallocate(wi, wp, we, supp, dw, d2w, xp)
    !stop

  end subroutine EvalHessians

  ! evaluate the gradient of piecewise linear function on element

  subroutine EvalGradP1(xp, dw )
    real, dimension(1:3, 1:3), intent(in) :: xp
    real, dimension(1:2), intent(inout) :: dw
    real, dimension(1:3):: det

    det(3) = (xp(1,1) - xp(2,1))*xp(3,2) + (xp(2,1) - xp(3,1))*xp(1,2) + (xp(3,1) - xp(1,1))*xp(2,2)
    det(1) = (xp(1,3) - xp(2,3))*xp(3,2) + (xp(2,3) - xp(3,3))*xp(1,2) + (xp(3,3) - xp(1,3))*xp(2,2)
    det(2) = (xp(1,1) - xp(2,1))*xp(3,3) + (xp(2,1) - xp(3,1))*xp(1,3) + (xp(3,1) - xp(1,1))*xp(2,3)

    dw(1:2) = det(1:2)/ det(3)

    !write(*, '(20es14.6)' ) xp(1, 1:3)
    !write(*, '(20es14.6)' ) xp(2, 1:3)
    !write(*, '(20es14.6)' ) xp(3, 1:3)
    !write(*, '(20es14.6)' ) det(1:3)
    !write(*, '(20es14.6)' ) dw(1:2)

    !stop

  end subroutine EvalGradP1

end module AMA_estims
