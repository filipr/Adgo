!> definition of porous media flow model
module modelPorous
  use main_data
  use f_mapping
  use mesh_oper
  use define_state
  use blocks_integ
  use lapack_oper
  use porous_fnc

  implicit none

  public:: Set_f_s_empty
  public:: Set_A_s_empty
  public:: Set_Ppm_empty
  public:: Set_R_s_porous
  public:: Set_K_sk_porous

  contains

  !> empty convective  terms
  subroutine Set_f_s_empty(ndimL, nbDim, Qdof, w, f_s, x, ie )
    integer, intent(in) :: Qdof, ndimL, nbDim
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
    real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
    integer, intent(in) :: ie


    f_s = 0.

  end subroutine Set_f_s_empty

  !> empty convective  terms
  subroutine Set_A_s_empty(ndimL, nbDim, Qdof, w, A_s, xi, ie)
    integer, intent(in) :: Qdof, nbdim, ndimL
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: A_s
    ! matrices A_s in  -- " --
    real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
    integer, intent(in) :: ie

    A_s = 0.

  end subroutine Set_A_s_empty

  !> empty convective  terms
  subroutine Set_Ppm_empty(ndimL, nbDim, Qdof, w, n, xi, Ppm, one_over_area, elem)
    integer, intent(in) :: Qdof, ndimL, nbDim
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: Ppm
    ! matrices Ppm in  -- " --
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: n   ! outer normal
    real, dimension(1:Qdof, 1:nbDim),intent(in) ::  xi          ! node on the edge?
    real, intent(in), optional :: one_over_area !
    type( element ), intent( inout ), optional :: elem !not used


    Ppm = 0.

  end subroutine Set_Ppm_empty


  !> compute viscous fluxes R_s, s=1,2 for the porous media model
  !> in integ nodes
  subroutine Set_R_s_porous(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
    integer, intent(in) :: ndimL, nbDim, iRe, Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi !physical cooedinates
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    !real, intent(in) :: Re_1                     ! inverse of Reynolds number
    real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
    real, dimension(:, :), allocatable :: K_sk
    integer :: k

    R_s(:, :, :) = 0.

    allocate(K_sk( 1:nbDim,  1:nbDim) )

    do k=1,Qdof
       
       call Eval_Diff_Porous_Coeffs(w(k,1), Dw(k, 1, 1:nbDim), K_sk(1:nbDim, 1:nbDim), &
            Re_1(1:iRe, k), 0, xi(k, 1:nbDim) )
       
       R_s(k, 1:nbDim, 1) =  matmul ( K_sk( 1:nbDim, 1:nbDim), Dw(k, 1, 1:nbDim) )
       
    enddo

    deallocate (K_sk)

  end subroutine Set_R_s_porous

  !> compute "matrices" 1x1 K_sk, s,k=1,2 for porous media model in integ nodes
  subroutine Set_K_sk_porous(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, K_sk, xi)
    integer, intent(in) :: ndimL, nbDim, iRe, Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi !physical cooedinates
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,ndimL), intent(inout) :: K_sk
    integer :: k

    K_sk(:, :, :, : ,:) = 0.

    do k=1,Qdof
       call Eval_Diff_Porous_Coeffs(w(k,1), Dw(k, 1, 1:nbDim), K_sk(k, 1:nbDim, 1:nbDim, 1, 1), &
            Re_1(1:iRe, k), 0, xi(k, 1:nbDim) )

    enddo
    
    
    !do k=1, Qdof
    !   if( abs( xi(k, 2) -0.) < 1) &
    !        write(63,'(30es12.4)') &
    !        xi(k, 1:2), K_sk(k, 1, 1, 1, 1), K_sk(k, 2, 2, 1, 1), K_sk(k, 1, 2, 1, 1), Re_1(2:iRe, k)
    !enddo

  end subroutine Set_K_sk_porous

  !> compute the matrix in front of the time derivative term in integ nodes
  subroutine Set_Time_Matrix_porous(elem, ndimL,  Qdof, wi, xi, TA, wR)
    type(element), intent (inout) :: elem
    integer, intent(in) :: ndimL,  Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: wi  !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:2+iRe), intent(in):: xi  !coodinates + participation to components
    real, dimension(1:Qdof, 1:ndimL, 1:ndimL), intent(inout) :: TA ! output matrix 
    real, dimension(1:ndimL, 1:elem%dof_plus), optional :: wR
    real, dimension(:, :), allocatable :: ww
    real :: rK, xc(2)
    integer :: i, l, dof
    real, dimension(:,:), pointer :: phi

    if(  state%model%idiff == 3 ) then
       do l=1, Qdof
          rK = 0.
          do i = 1, 3  ! SET THE NUMBER of materials
             ! elem%xi(0, l, 2+i) is the coefficent of affiliation to  i-th material component
             if(elem%xi(0, l, 2+i) > 0.) & 
                  rK = rK +  capacity(wi(l,1), i, xi(l, 2) ) * xi(l, 2+i)
             !rK = rK +  capacity(wi(l,1), i, elem%xi(0,l, 2) ) * elem%xi(0, l, 2+i)

             !write(*, '(a8,i5, 6es12.4)') 'conduct:', i, u, uu, xi(2), Re_1(i+1), rK
          enddo
          TA(l,1, 1) = rK
          !if(rk < 0.) &
          !write(21, '(23es12.4)' ) elem%xi(0,l,1:2), rK,  elem%xi(0,l,3:)

          !if(state%time%iter >= 1) then
          !   write(66,'(30es12.4)')  &
          !        elem%xi(0, l, 1:2), rK, wi(l, 1),elem%xi(0, l, 2+i)
          !   if(elem%i == grid%nelem) stop '9ue93jdo3dmzd39u393i'
          !endif

       enddo

    elseif(  state%model%idiff == 4 ) then
       ! C1
       !TA(1:Qdof,1, 1) = 1. !* wi(1:Qdof, 1)

       ! C2
       !TA(1:Qdof,1, 1) = 2. * wi(1:Qdof, 1)

       ! C3
       !TA(1:Qdof,1, 1) = 2. * state%time%ctime / wi(1:Qdof, 1)

       ! C4
       TA(1:Qdof,1, 1) = 2  * wi(1:Qdof, 1)

       !!print*,'eleme:', elem%i, state%time%ctime,  TA(1:2,1, 1)
    endif


    ! xc(1) = 10.5; xc(2) =0.25
    ! do l=1, Qdof
    !    !if( abs( elem%xi(0, l, 1) -12.) < 3 .and. abs( elem%xi(0, l, 2) -0.) < 3) &
    !    !if(TA(l,1,1) < 0.) &
    !    if( dot_product(elem%xc- xc , elem%xc- xc) < 2.) &
    !        write(62,'(30es12.4)') elem%xi(0,l, 1:2), TA(l,1,1)
    ! enddo


    ! ! projection of the integ nodes given functions onto polynomial function
    dof = elem%dof_plus
    !dof = elem%dof
    !dof = DoFtriang( max(0, elem%deg - 1) )
    allocate(ww(1:dof, 1:2) )

    !!if( abs( TA(1, 1, 1) - TA(2, 1,1 ) ) > 0.15 ) then
    !if(elem%i == 843) then
    !   write(*,'(a12, 2i5, 300es12.4)') ' TA:', elem%i, Qdof,  TA(1:Qdof, 1, 1)
    !endif

    ww(:,:) = 0.
    call IntegrateVectorB(elem, dof, TA(1:Qdof, 1, 1), ww(1:dof, 1) )

    !write(*,'(a12, 2i5, 300es12.4)') ' w,phi:', elem%i, dof,  ww(1:dof, 1)

    !do l=1,dof
    !   ww(l, 2) = dot_product(elem%MassInv%Mb(l,1:dof), ww(1:dof, 1) )
    !enddo
    ww(1:dof, 2) = ww(1:dof, 1)

    call SolveLocalMatrixProblem(dof, elem%mass%Mb(1:dof, 1:dof), 1, ww(1:dof, 2))

    if( present(wR) ) then
       wR = 0.
       wR(1, 1:dof) = ww(1:dof, 2)
    endif


    !write(*,'(a12, 2i5, 300es12.4)') 'M^-1 w,phi:', elem%i, dof,  ww(2, 1:dof)

    !phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)
    !TA(1:Qdof, 1, 1) = matmul(ww(1:dof, 2) , phi(1:dof,1:Qdof) )


    deallocate(ww)

    ! do l=1, Qdof
    !    !if( abs( elem%xi(0, l, 1) -12.) < 3 .and. abs( elem%xi(0, l, 2) -0.) < 3) &
    !   !if(TA(l,1,1) < 0.) &

    !    if( dot_product(elem%xc- xc , elem%xc- xc) < 2.) &
    !        write(64,'(30es12.4)') elem%xi(0,l, 1:2), TA(l,1,1)
    ! enddo


    !if( abs( TA(1, 1, 1) - TA(2, 1,1 ) ) > 0.15 ) then
    !if(elem%i == 843) then
    !   write(*,'(a12, 2i5, 300es12.4)') ' TA:', elem%i, dof,  TA(1:Qdof, 1, 1)
    !   print*,"_____________________",  abs( TA(1, 1, 1) - TA(2, 1,1 ) )
    !   !stop '9e39ud93o'
    !endif

  end subroutine Set_Time_Matrix_porous




  !> evaluate "exact" solution, used for the IC and BC
  subroutine Exact_Porous(Qdof, x, wi, t)
    integer, intent(in) :: Qdof
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: x
    real, dimension(1:Qdof, 1:ndim), intent(out) :: wi
    real, intent(in) :: t
    integer :: imod
    real, dimension(:,:), allocatable :: rr
    real :: u_max, u_min, r_max, tt, t_max
    integer :: l

    allocate( rr(1:2, 1:Qdof) )
    !call Set_Model_Data(x, t, wi, 1)
    imod = state%model%iexact

    !print*,'imod = ', imod, 'idiff = ', state%model%idiff

    ! imod = 1 damp (hraz) 

    select case (imod)
    case(1)
       wi(:, 1) = x(:, 1) 

    case(2)
       u_min = 0.

       wi(:, :) = u_min
       u_max = 30. - u_min
       r_max = 0.5
       
       do l=1,Qdof
          ! distance from the left boundary
          if(x(l, 1) <= 0.) then
             rr(1, l) = abs( x(l, 2) )
          elseif( x(l, 2) < 0. .and. x(l, 1) < -15./13 * x(l, 2) ) then
             rr(1,l) = sqrt( dot_product(x(l, 1:2), x(l, 1:2) ) )
          else
             rr(1,l) =  (15./13*x(l, 1) - x(l, 2) ) / sqrt(394./ 169) 
          endif
          !write(21, *) x(l, 1:2), rr(1,l)

          if( rr(1, l) < r_max) then
             wi( l, 1) = u_max * ( cos(rr(1, l )/r_max *pi ) + 1 )/ 2. + u_min
          endif

          t_max = -0.5
          if(t < t_max) then
             tt = sin(t * pi /2 /t_max)
             wi(:, 1) = wi(:, 1) * tt
          endif

          !write(33, *) x(l, 1:2), wi(l, 1), rr(1,l)
       enddo

    case(3) ! test nonlinear case: 2u u_t - (u^2 u_x)_x = 0
       ! C1
       !wi(1:Qdof, 1) = x(1:Qdof, 1)**2 + t

       ! C2
       !wi(1:Qdof, 1) = x(1:Qdof, 1) + t

       ! C3
       !wi(1:Qdof, 1) = exp( t * x(1:Qdof, 1)) 

       ! C4
       wi(1:Qdof, 1) = exp( t +  x(1:Qdof, 1) )


       !wi(1:Qdof, 1) = x(1:Qdof, 1) 
       !if(t == 0) wi(1:Qdof, 1) = wi(1:Qdof, 1) *0.5
    case(4:)
 
        stop 'UNKNOWN type in Exact_Porous'

    end select


  end subroutine Exact_Porous



  !> evaluation of diffusion coefficients and their derivatives
  !> \f$ K_{s,k}^{i,j},\ s,k=1,2 (\mbox{space dimension}),\ i,j=1,\dots, ndim\f$,
  !> \f$ ider =0 \Rightarrow K(u),\f$ or \f$ ider =0 \Rightarrow K(|\nabla u|),\f$
  !> \f$ ider =1 \Rightarrow \frac{\rm d}{{\rm d} u} K(u) \f$ or
  !> \f$ ider =1 \Rightarrow \frac{\rm d}{{\rm d} |\nabla u|} K(|\nabla u|) \f$
  subroutine Eval_Diff_Porous_Coeffs(u, Du, K_sk, Re_1, ider, xi)
    real, intent(in) :: u            ! solution
    real, dimension(1:nbDim), intent(in) :: Du ! derivative of the solution
    real, dimension(1:nbDim), intent(in) :: xi ! physical coordinate
    real, dimension(1:nbDim, 1:nbDim), intent(inout) :: K_sk ! output diffusion matrix
    real, dimension(1:iRe), intent(in) :: Re_1         ! viscosity
    integer, intent(in) :: ider      ! =0 => K(u), =1 => d K(u)/ d u
    integer :: i, j, imod, nn     ! IMOD
    real :: m, uu, rK, rKp, val1, val2, val3
    real :: viscos, compress, permeab, a0, a1

    imod = state%model%idiff
    !imod = 1    ! Laplace, linear diffusion
    !imod = 2    ! linear diffusion with different coeficients
    !imod = 3    ! nonlinear diffusion, atan
    !imod = 4    ! nonlinear diffusion, atan, anisotrop
    !imod = 5    ! Kacur: degenerate parabolic problem (Eymard, Hilhorst, Vohralik 2006)
    !imod = 6    ! Barenblatt, porus media flow, Radu et all 2008
    !imod = 7    ! NONLINEAR elliptic [Houston, Sulli, Robson 2007]
    !imod = 8    ! NONLINEAR elliptic [Houston, Sulli, Robson 2007] second

    K_sk(:, :) = 0.
    select case (imod)
    case(0)   ! no diffusion
       K_sk(:, :) = 0.

    case(1)     ! linear diffusion
       if(ider == 0) then
          ! functions
          K_sk(1, 1) = Re_1(1)
          K_sk(2, 2) = Re_1(1)
       else
          ! derivatives
          ! K_sk = 0.
       endif

    case(2)     ! non-linear diffusion, original test problem

       permeab = Re_1(2)  ! uses the precomputed value
       
       viscos = 1.3E-03
       compress = 5E-10

       a0 = viscos / permeab
       a1 = 550 /sqrt(permeab )
       !a1 = 0.
       
       if(( a0 + sqrt(a0*a0 + 4 * a1 * uu)) <=  0.) &
            write(*,'(a8,4es12.4)') 'permeaB:', permeab, a0, a1

       uu = sqrt( Du(1)*Du(1) + Du(2)*Du(2) )

       rK = Re_1(1) * 2./ ( a0 + sqrt(a0*a0 + 4 * a1 * uu)) 

       rKp = - Re_1(1) * 1./ ( a0 + sqrt(a0*a0 +4 * a1 * uu))**2  / sqrt(a0*a0+ 4. * a1*  uu) * 4 *a1

       !write(20, '(6es12.4)') xi(1:2), rK / compress

       if(ider == 0) then
          ! functions
          K_sk(1, 1) = rK / compress
          K_sk(2, 2) = K_sk (1,1)

          !print*,'#DE#DE#',  rK / compress
       else
          ! derivatives
          K_sk(1, 1) = rKp / compress
          K_sk(2, 2) = K_sk (1,1)

          print* ,'NOT NECESSARY'
          print*,'#E#E#E:', rKp / compress
       endif


    case(3)     ! non-linear diffusion, original test problem

       ! size of the gradient of the head preasure
       uu = sqrt(dot_product( Du(1:2) , Du(1:2) )  )

       rK = 0.
       do i = 1, 3  ! SET THE NUMBER of materials
          if(Re_1(i+1) >0.) &
               rK = rK +  forch_conduct(u, uu, i, xi(2) ) * Re_1(i+1)

          ! if(xi(2) > -0.2 .and. xi(2) <= 0.0 .and. abs(xi(1) -28) <= 0.25) then
          !    write(*, '(a8,i5, 16es12.4)') 'conduct:', &
          !         i, u, uu, xi(1:2), Re_1(i+1), forch_conduct(u, uu, i, xi(2) ), rK
          !    if(i == 3) print*
          !    if(i == 3) write(65, *) xi(1:2), rK
          ! endif

       enddo
       
       !write(22, '(30es12.4)' ) xi(1:2), rK

       if(ider == 0) then
          ! functions
          K_sk(1, 1) = rK 
          K_sk(2, 2) = K_sk (1,1)

          !print*,'#DE#DE#',  rK / compress
       else
          ! derivatives
          !K_sk(1, 1) = rKp / compress
          !K_sk(2, 2) = K_sk (1,1)

          print* ,'NOT NECESSARY'
       endif

       ! print*,'(ID#(UJD(#JU#OEJK#O'
       ! rK = 30.
       ! nn = 10000
       ! do i = -nn, nn
       !    uu = 1.*i /nn * rK
       !    do j=1, 3
       !       write(90+j,*) -uu, capacity(0., j, uu), forch_conduct(0., 0., j, uu) !vangen(0., 1, uu)/soilpar(1)%Ths*soilpar(1)%Ss
       !    enddo

       ! enddo
       ! stop '3u9439iewjs'

   case(4)     ! test nonlinear case: 2u u_t - (u^2 u_x)_x = 0
       if(ider == 0) then
          ! functions
          ! C1
          !K_sk(1, 1) = 1./2 
          
          ! C2
          !K_sk(1, 1) = u * u 

          ! C3
          !K_sk(1, 1) = xi(1)**2  / u 

          ! C4
          K_sk(1, 1) =  u 


          K_sk(2, 2) = 0.
       else
          ! derivatives
          ! K_sk = 0.
       endif



    case(13)     !porous media flow,  Forchheimer 2-term law
       ! ider
       !call  Set_porous_media_A2(xi(2), -xi(1), val1, val2, val3)

       !if( val1/Re_1(2) > 1E+3 .or. val1/Re_1(2) < 1E-3) then
       !   write(*,'(a8,8es12.4)') 'diffgfee:', val1, Re_1(2), val1/Re_1(2), xi(1:2)
       !endif

       !permeab = val1
       permeab = Re_1(2)  ! uses the precomputed value

       viscos = 1.3E-03
       compress = 5E-10

       a0 = viscos / permeab
       a1 = 550 /sqrt(permeab )
       !a1 = 0.
       
       if(( a0 + sqrt(a0*a0 + 4 * a1 * uu)) <=  0.) &
            write(*,'(a8,4es12.4)') 'permeaB:', permeab, a0, a1

       uu = sqrt( Du(1)*Du(1) + Du(2)*Du(2) )

       rK = Re_1(1) * 2./ ( a0 + sqrt(a0*a0 + 4 * a1 * uu)) 

       rKp = - Re_1(1) * 1./ ( a0 + sqrt(a0*a0 +4 * a1 * uu))**2  / sqrt(a0*a0+ 4. * a1*  uu) * 4 *a1

       write(20, '(6es12.4)') xi(1:2), rK / compress

       if(ider == 0) then
          ! functions
          K_sk(1, 1) = rK / compress
          K_sk(2, 2) = K_sk (1,1)

          !print*,'#DE#DE#',  rK / compress
       else
          ! derivatives
          K_sk(1, 1) = rKp / compress
          K_sk(2, 2) = K_sk (1,1)

          print*,'#E#E#E:', rKp / compress
       endif


    case(14:)
       stop 'UNKNOWN TYPE in Eval_Diff_Porous_Coeffs'

    end select

  end subroutine Eval_Diff_Porous_Coeffs


 !> setting of Hraz by Michal Kuraz, OLD SETTING
  subroutine Set_porous_media_Hraz_OLD(xii, yii, val1, val2, f)
    real, intent(in) :: xii, yii
    real, intent(inout) :: val1, val2, f
    real, dimension(:,:), allocatable :: x
    real, dimension(:,:), allocatable :: pq
    real, dimension (:), allocatable :: xi
    real ::  qq1, qq2, qq
    integer :: ityp
    integer :: i
    !real:: eps = 1E-5
    real:: eps = 1E-2


    allocate(pq(1:3, 1:3) )

    ! coefficient pq(:, 1 ) = a0 =   K / mu * (1/kappa)  Darcy
    ! coefficient pq(:, 2 ) = a1 =   K / mu * (1/kappa)  Darcy
    ! K = permeability of the medium = 1E-12 ..  1E-15
    ! mu = viscosity of the fluid   = 1.3E-3  (for water)
    ! 1/kappa = compressibility 5E-10 for water
    
    pq(1, 1:3) = (/ 1E-14, 1E-14, 0. /)
    pq(2, 1:3) = (/ 1E-12, 1E-12, 0. /)
    pq(3, 1:3) = (/ 1E-10, 1E-10, 0. /)


    allocate(x(1:16, 1:2) )
    
    x( 1, 1:2) = (/-7. ,  -12. /)
    x( 2, 1:2) = (/ 14.8 ,-12. /)
    x( 3, 1:2) = (/ 15.2 ,-12. /)
    x( 4, 1:2) = (/ 37. , -12. /)
    x( 5, 1:2) = (/ 37. ,  0. /)
    x( 6, 1:2) = (/ 30.  ,  0./)
    x( 7, 1:2) = (/ 17. ,  15./)
    x( 8, 1:2) = (/ 13. ,  15. /)
    x( 9, 1:2) = (/  0. ,  -0.0/)
    x(10, 1:2) = (/ -7. ,  0./)
    x(11, 1:2) = (/ 15.2 , -5. /)
    x(12, 1:2) = (/ 18. ,  0.0/)
    x(13, 1:2) = (/ 16. , 14.0 /)
    x(14, 1:2) = (/ 14. , 14. /)
    x(15, 1:2) = (/ 12. , 0. /)
    x(16, 1:2) = (/ 14.8, -5. /)


    allocate(xi(1:2) )
    xi(1) = xii
    xi(2) = yii


    ! NEW variant with boundary smoothing
    val1 = 0.;    val2 = 0;   f = 0.;

    ! left bottom 1
    qq = InsideQuadrilaterall(xi, x(1,:), x(2,:), x(16,:), x(10,:), eps, .false.)
    ityp = 1; val1 = val1 + qq * pq(ityp, 1);   val2 = val2 + qq * pq(ityp, 2); f = f + qq * pq(ityp, 3)

    ! left bottom 2
    qq = InsideTriangle(xi, x(10,:), x(16,:), x(15,:), eps, .true.)
    ityp = 1; val1 = val1 + qq * pq(ityp, 1);   val2 = val2 + qq * pq(ityp, 2); f = f + qq * pq(ityp, 3)

    ! right bottom 1
    qq = InsideQuadrilaterall(xi, x(5,:), x(11,:), x(3,:), x(4,:), eps, .false.)
    ityp = 1; val1 = val1 + qq * pq(ityp, 1);   val2 = val2 + qq * pq(ityp, 2); f = f + qq * pq(ityp, 3)

    ! right bottom 2
    qq = InsideTriangle(xi, x(11,:), x(5,:), x(12,:), eps, .true.)
    ityp = 1; val1 = val1 + qq * pq(ityp, 1);   val2 = val2 + qq * pq(ityp, 2); f = f + qq * pq(ityp, 3)

    ! top left
    qq = InsideQuadrilaterall(xi, x(9,:), x(15,:), x(14,:), x(8,:), eps, .false.)
    ityp = 2; val1 = val1 + qq * pq(ityp, 1);   val2 = val2 + qq * pq(ityp, 2); f = f + qq * pq(ityp, 3)

    ! top top
    qq = InsideQuadrilaterall(xi, x(8,:), x(14,:), x(13,:), x(7,:), eps, .false.)
    ityp = 2; val1 = val1 + qq * pq(ityp, 1);   val2 = val2 + qq * pq(ityp, 2); f = f + qq * pq(ityp, 3)

    ! top right
    qq = InsideQuadrilaterall(xi, x(7,:), x(13,:), x(12,:), x(6,:), eps, .false.)
    ityp = 2; val1 = val1 + qq * pq(ityp, 1);   val2 = val2 + qq * pq(ityp, 2); f = f + qq * pq(ityp, 3)

    ! core top 
    qq = InsideQuadrilaterall(xi, x(12,:), x(13,:), x(14,:), x(15,:), eps, .true.)
    ityp = 3; val1 = val1 + qq * pq(ityp, 1);   val2 = val2 + qq * pq(ityp, 2); f = f + qq * pq(ityp, 3)

    ! core center
    qq = InsideQuadrilaterall(xi, x(11,:), x(12,:), x(15,:), x(16,:), eps, .true.)
    ityp = 3; val1 = val1 + qq * pq(ityp, 1);   val2 = val2 + qq * pq(ityp, 2); f = f + qq * pq(ityp, 3)

    ! core tube
    qq = InsideQuadrilaterall(xi, x(3,:), x(11,:), x(16,:), x(2,:), eps, .false.)
    ityp = 3; val1 = val1 + qq * pq(ityp, 1);   val2 = val2 + qq * pq(ityp, 2); f = f + qq * pq(ityp, 3)




    ! if(qq1 <= 0. .and. qq2 <= 0. ) then
    !write(99,*) xi(1), xi(2),  val1

    !  elseif(qq1 >= 1. .or. qq2 >= 1. ) then

    ! !    write(22,*) -xi(2), xi(1),  qq1, qq2

    !  else
    !     write(*,*) -xi(2), xi(1), qq1,qq2, qq
    !  endif

    deallocate(x,  xi, pq)

  end subroutine Set_porous_media_Hraz_OLD


 !> setting of Hraz by Michal Kuraz, returns the relative affiliation to each component
  subroutine Set_porous_media_Hraz(xii, yii, iRe, vals ) 
    real, intent(in) :: xii, yii
    integer, intent(in) :: iRe
    real, dimension(1: iRe),intent(inout) :: vals
    real, dimension(:,:), allocatable :: x
    real, dimension(:,:), allocatable :: pq
    real, dimension (:), allocatable :: xi
    real ::  qq1, qq2, qq, norm
    integer :: ityp
    integer :: i
    !real:: eps = 1E-5
    !real:: eps = 1E-3
    !real:: eps = 1E-2
    !real:: eps = 1E-1
    real:: eps = 5E-1


    allocate(pq(1:3, 1:3) )

    ! coefficient pq(:, 1 ) = a0 =   K / mu * (1/kappa)  Darcy
    ! coefficient pq(:, 2 ) = a1 =   K / mu * (1/kappa)  Darcy
    ! K = permeability of the medium = 1E-12 ..  1E-15
    ! mu = viscosity of the fluid   = 1.3E-3  (for water)
    ! 1/kappa = compressibility 5E-10 for water
    
    pq(1, 1:3) = (/ 1E-14, 1E-14, 0. /)
    pq(2, 1:3) = (/ 1E-12, 1E-12, 0. /)
    pq(3, 1:3) = (/ 1E-10, 1E-10, 0. /)


    allocate(x(1:18, 1:2) )
    
    x( 1, 1:2) = (/  -7.5, -11.  /)
    x( 2, 1:2) = (/ 14.8 , -11. /)
    x( 3, 1:2) = (/ 15.2 , -11. /)
    x( 4, 1:2) = (/ 37.5 , -11. /)
    x( 5, 1:2) = (/ 37.5 ,  0. /)
    x( 6, 1:2) = (/ 30.  ,  0./)
    x( 7, 1:2) = (/ 17. ,  15./)
    x( 8, 1:2) = (/ 13. ,  15. /)
    x( 9, 1:2) = (/  0. ,  -0.0/)
    x(10, 1:2) = (/ -7.5,  0./)
    x(11, 1:2) = (/ 15.2 , -5. /)
    x(12, 1:2) = (/ 18. ,  0.0/)
    x(13, 1:2) = (/ 16. , 14.0 /)
    x(14, 1:2) = (/ 14. , 14. /)
    x(15, 1:2) = (/ 12. , 0. /)
    x(16, 1:2) = (/ 14.8, -5. /)

    x(17, 1:2) = (/ 0. , -11. /)
    x(18, 1:2) = (/ 30., -11. /)


    allocate(xi(1:2) )
    xi(1) = xii
    xi(2) = yii


    vals(:) = 0.


!     x(10, 1:2) = (/ -0.1 , 0. /)
!     x(1, 1:2) = (/ 0.1 , 20.1 /)
!     xi(1:2) =  (/ 0.0025 , -10.002 /)


    ! left bottom 1
    qq = InsideQuadrilaterall(xi, x(10,:), x(1,:), x(17,:), x(9,:), eps, .false.)
    ityp = 3;   vals(ityp) = vals(ityp) + qq

    ! left bottom 2
    qq = InsideTriangle(xi, x(2,:), x(16,:), x(17,:), eps, .false.)
    ityp = 3;   vals(ityp) = vals(ityp) + qq

    ! left bottom 3
    qq = InsideQuadrilaterall(xi, x(17,:), x(16,:), x(15,:), x(9,:), eps, .true.)
    ityp = 3;   vals(ityp) = vals(ityp) + qq


    ! right bottom 1
    qq = InsideQuadrilaterall(xi, x(6,:), x(18,:), x(4,:), x(5,:), eps, .false.)
    ityp = 3;   vals(ityp) = vals(ityp) + qq

    ! right bottom 2
    qq = InsideTriangle(xi, x(18,:), x(11,:), x(3,:), eps, .false.)
    ityp = 3;   vals(ityp) = vals(ityp) + qq

    ! right bottom 3
    qq = InsideQuadrilaterall(xi, x(6,:), x(12,:), x(11,:), x(18,:), eps, .true.)
    ityp = 3;   vals(ityp) = vals(ityp) + qq

    ! top left
    qq = InsideQuadrilaterall(xi, x(9,:), x(15,:), x(14,:), x(8,:), eps, .false.)
    if(xi(1) < x(9,1)) qq = 0.
    ityp = 1;   vals(ityp) = vals(ityp) + qq

    ! top top
    qq = InsideQuadrilaterall(xi, x(8,:), x(14,:), x(13,:), x(7,:), eps, .false.)
    ityp = 1;   vals(ityp) = vals(ityp) + qq

    ! top right
    qq = InsideQuadrilaterall(xi, x(7,:), x(13,:), x(12,:), x(6,:), eps, .false.)
    if(xi(1) > x(6,1)) qq = 0.
    ityp = 1;   vals(ityp) = vals(ityp) + qq

    ! core top 
    qq = InsideQuadrilaterall(xi, x(12,:), x(13,:), x(14,:), x(15,:), eps, .true.)
    ityp = 2;   vals(ityp) = vals(ityp) + qq

    ! core center
    qq = InsideQuadrilaterall(xi, x(11,:), x(12,:), x(15,:), x(16,:), eps, .true.)
    ityp = 2;   vals(ityp) = vals(ityp) + qq

    ! core tube
    qq = InsideQuadrilaterall(xi, x(3,:), x(11,:), x(16,:), x(2,:), eps, .false.)
    ityp = 2;   vals(ityp) = vals(ityp) + qq


    norm = sum(vals(:) )
    vals(:) = vals(:) / norm


    ! if(qq1 <= 0. .and. qq2 <= 0. ) then
    !write(99,*) xi(1), xi(2),  val1

    !  elseif(qq1 >= 1. .or. qq2 >= 1. ) then

    ! !    write(22,*) -xi(2), xi(1),  qq1, qq2

    !  else
    !     write(*,*) -xi(2), xi(1), qq1,qq2, qq
    !  endif

    deallocate(x,  xi, pq)

  end subroutine Set_porous_media_Hraz


end module modelPorous

