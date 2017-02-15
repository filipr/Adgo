
!> definition of shallow water equations
module modelSWE
  
  use data_mod
  use main_data
  use f_mapping
  use mesh_oper
  use blocks_integ
  use matrix_oper_int
  use pedes_averaging

  implicit none
  
  
  public:: Set_f_s_swe
  public:: Set_A_s_swe
  public:: Set_Ppm_swe

  public:: Set_R_s_swe
  public:: Set_K_sk_swe

  public:: Set_S_swe
  public:: Set_DS_swe

contains



  !> compute the swe flow fluxes \f$ f_s({\bf w}),\quad s=1,2\f$
  subroutine Set_f_s_swe(ndimL, nbDim, Qdof, w, f_s, xi, ie)
    !Qdof ... pocet konecnych objemu
    !nbDim ... dimenze rychlosti = 2
    !ndimL ... dimenze vystupu = 3
    integer, intent(in) :: Qdof, nbDim, ndimL
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in integ nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s ! fluxes f_s in  -- " --
    real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
    integer, intent(in) :: ie  ! index of the element
    real, dimension(:), allocatable :: u,v,p,gterm
    real, dimension(:,:,:,:), allocatable :: A_s
    integer:: i
    real :: h_minimal_allowed, fac

    h_minimal_allowed = state%model%Pr
 
    allocate( u(1:Qdof), v(1:Qdof), gterm(1:Qdof) )


    do i=1, Qdof  
       !i konecny objem (trojuhelnik)
  
       !print*,'####', h_minimal_allowed, i, Qdof
       if(w(i, 1) <= 0.) then
          !if(w(i, 1) <= h_minimal_allowed) then
          f_s(i, :, :) = 0.

       else
          !call Set_limit_factor(w(i,1), fac)
          fac = 1.

          u(i) = w(i, 2) / w(i, 1)
          v(i) = w(i, 3) / w(i, 1)

          gterm(i) = 0.5*9.81*w(i, 1)**2

          !f_s(konecny objem,[1,2]..s, [1,2,3]..dim of output)
          f_s(i, 1, 1) = w(i, 2)
          f_s(i, 1, 2) = w(i, 2) * u(i) + gterm(i)
          f_s(i, 1, 3) = w(i, 2) * v(i)

          f_s(i, 2, 1) = w(i, 3)
          f_s(i, 2, 2) = w(i, 3) * u(i)
          f_s(i, 2, 3) = w(i, 3) * v(i) + gterm(i)

          !if(w(i,1) < 1E-5) write(*,'(a8,40es12.4)') 'de f:s',w(i,1:3), u(i), v(i),p(i),f_s(i, 1, 1:3)

          f_s(i, :, :) = f_s(i, :, :) * fac
       endif


    enddo

    !!endif

    deallocate(u,v,gterm )

  end subroutine Set_f_s_swe


  !> compute matrices \f$ A_s = \frac{D{\bf f_s}({\bf w})}{D{\bf w}},\quad s=1,2\f$
  !> for swe flow fluxes
  subroutine Set_A_s_swe(ndimL, nbDim, Qdof, w, A_s, xi, ie)
    !Qdof ... pocet konecnych objemu
    !nbDim ... dimenze rychlosti = 2
    !ndimL ... dimenze vystupu = 3
    integer, intent(in) :: ndimL, nbDim, Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: A_s
                                               ! matrices A_s in  -- " --
    real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
    integer, intent(in) :: ie  ! index of the element
    real, dimension(:), allocatable :: u,v,h, pp, uv, u2,v2,ggterm
    integer :: i
    real :: h_minimal_allowed, fac

    !TODO -- podivat se na state%model%Pr
    h_minimal_allowed = state%model%Pr


    allocate( u(1:Qdof), v(1:Qdof), h(1:Qdof), ggterm(1:Qdof) )
    allocate( uv(1:Qdof), u2(1:Qdof), v2(1:Qdof) )
    
       
    do i=1,Qdof
       h(i) = w(i,1)
       !if( h(i) > 0. ) h(i) = max(h(i), h_minimal_allowed)

       if(h(i) > 0.   ) then

          !call Set_limit_factor(h(i), fac)
          fac = 1.

          !if(h(i) > h_minimal_allowed) then
          u(i) = w(i,2)/h(i)
          v(i) = w(i,3)/h(i)
       else
          u(i) = 0.
          v(i) = 0.
       endif


       if(h(i) <= 0.) then
       !if(h(i) <= h_minimal_allowed) then
          A_s(i, :,:,:) = 0.
       else

!          pp(i) =  state%model%p0 *  h(i)**state%model%kappa1 ! IS NOT PRESSURE !!!!
          ggterm(i)=9.81*w(i,1)
          uv(i) = u(i)*v(i)
          u2(i) = u(i)*u(i)
          v2(i) = v(i)*v(i)
          
          !A_s(element, kolikataMatice, index1, index2)
          A_s(i, 1, 1, 1) = 0.
          A_s(i, 1, 1, 2) = 1.
          A_s(i, 1, 1, 3) = 0.

          A_s(i, 1, 2, 1) = -u2(i) + ggterm(i)
          A_s(i, 1, 2, 2) =  2 * u(i)
          A_s(i, 1, 2, 3) = 0.

          A_s(i, 1, 3, 1) = -uv(i)
          A_s(i, 1, 3, 2) = v(i)
          A_s(i, 1, 3, 3) = u(i)


          A_s(i, 2, 1, 1:2) = 0.
          A_s(i, 2, 1, 3) = 1.

          A_s(i, 2, 2, 1) = -uv(i)
          A_s(i, 2, 2, 2) = v(i)
          A_s(i, 2, 2, 3) = u(i)

          A_s(i, 2, 3, 1) = -v2(i) + ggterm(i)
          A_s(i, 2, 3, 2) = 0.
          A_s(i, 2, 3, 3) = 2 * v(i)

          !if(w(i,1) < 1E-5) write(*,'(a8,40es12.4)') 'de A:s',w(i,1:3), u(i),v(i),pp(i),u2(i),uv(i),v2(i)
          A_s(i, :,:,:) = A_s(i, :,:,:) * fac

       endif
    enddo

    deallocate(u,v,h, ggterm, uv, u2, v2 )

  end subroutine Set_A_s_swe



  !> compute matrices
  !> \f$ P^{\pm} = \left(\frac{D({\bf f_1}({\bf w})n_1+{\bf f_2}({\bf w})n_2}{D{\bf w}}
  !>  \right)^{\pm}\f$
  !> for the swe flow equation
  subroutine Set_Ppm_swe(ndimL, nbDim, Qdof, w, n, xi, Ppm, one_over_area, elem)
    !Qdof ... pocet konecnych objemu
    !nbDim ... dimenze rychlosti = 2
    !ndimL ... dimenze vystupu = 3
    integer, intent(in) :: Qdof, ndimL, nbDim
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: Ppm
                                               ! matrices Ppm in  -- " --
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: n   ! outer normal
    real, dimension(1:Qdof, 1:nbDim),intent(in) ::  xi                    ! node on the edge?
    real, intent(in), optional :: one_over_area
    type(element), intent(inout), optional :: elem
    real, dimension(:,:,:), allocatable :: PPP
    real, dimension(:,:), allocatable :: t, t1, q, q1  !,  DD, DDp, DDm
    real, dimension(:), allocatable   :: dp, dm, nv
    real :: rlen, eta, u, v,pp, ggterm, p0, aa, kappa, kappa1, val, fac, fac_r
    real :: eta_minimal_allowed
    integer :: ie, i, j, itest1, iprint
    real :: tt, tt1, velocity, max_velocity, h

    itest1 = -grid%nelem

    call cpu_time(tt)

    if(ndimL /= 3) stop 'Bad dimension (ndim) in subroutine Set_Ppm_swe !!'

    allocate( t(1: 3 ,1: 3 ), t1(1: 3 ,1: 3 ), q(1: 3 ,1: 3 ), q1(1: 3 ,1: 3 ) )
    allocate( PPP( 1:2, 1: 3 ,1: 3 ) )

    allocate( dp(1: 3 ), dm(1: 3 ), nv(1:nbDim))



    if (.not. present(elem) .or. .not. present(one_over_area) ) &
      stop 'elem and one_over_area must be present in Set_Ppm_swe!'

    kappa = state%model%kappa  
    kappa1 = state%model%kappa1  
    p0 =  state%model%p0 
    eta_minimal_allowed = state%model%Pr
    max_velocity = 10.

    iprint = 0

    do ie=1,Qdof

       !normal vector size + normalized normal vector
       rlen = (n(ie,1)*n(ie, 1) + n(ie, 2)*n(ie, 2) )**0.5
       nv(1:nbDim) = n(ie, 1:nbDim)/rlen

       ! height of water
       h = w(ie, 1)

       if(h <= 0. )  then
       !if(h <= h_minimal_allowed) then
          Ppm(ie,1:2, 1: 3 , 1: 3 ) = 0.

       else

          !call Set_limit_factor(h, fac)
          fac = 1.

          !transformation 
          u = ( w(ie, 2) * nv(1) + w(ie, 3) * nv(2) ) / h
          v = (-w(ie, 2) * nv(2) + w(ie, 3) * nv(1) ) / h

          
          !aa = c in the 'On a numerical flux for the shallow water equations' 
          aa = sqrt(9.81 * h)


          ! eigenvalues
          ! TODO, rlen does not belong into eigenvalues, but would pop up later because of the method, right? 
          dp(1) = (u - aa) *rlen
          dp(2) =    u     *rlen
          dp(3) = (u + aa) *rlen


          !calculating the negative (dm) and the positive (dp)  parts of eigenvalues
          do i=1, 3 
             dm(i)=0.
             if(dp(i) < 0.)then
                dm(i)=dp(i)
                dp(i)=0.
             endif
          enddo

          elem%max_eigenvals = max( abs(dp(1)), abs(dp(3) ) ) *one_over_area
          if(h > 0.01) &
               state%max_eigenvals = max(state%max_eigenvals, elem%max_eigenvals  )


          ! matrix with eigenvectors
          t(1, 1) = 1. 
          t(1, 2) = 0.     
          t(1, 3) = 1.   

          t(2, 1) = u-aa
          t(2, 2) = 0.  !TODO 0 instead of 1 
          t(2, 3) = u+aa  

          t(3, 1) = v      
          t(3, 2) = 1.  
          t(3, 3) = v    

          ! its inverse          
          t1(1, 1) = (aa+u)/(2*aa)  
          t1(1, 2) = (-1.)/(2*aa) 
          t1(1, 3) = 0.   

          t1(2, 1) = -v
          t1(2, 2) = 0.  
          t1(2, 3) = 1.  

          t1(3, 1) = (aa-u)/(2*aa)   
          t1(3, 2) = (1.)/(2*aa)    
          t1(3, 3) =  0.    


          !transformation matrixes
          q(1:3,1:3) = 0. 
          q(1,1) = 1.
          
          q(2,2) =  nv(1)
          q(2,3) = nv(2)
          
          q(3,2) = -nv(2)
          q(3,3) = nv(1)

          !its inverse
          q1(1:3,1:3 ) = q(1:3,1:3) 

          q1(2,3) = -nv(2)
          
          q1(3,2) =  nv(2)



          !PPP are matrices A1+, A1-
          do i=1, 3 
             do j=1, 3 
                PPP(1, i, j) = sum(t(i, 1: 3 ) * dp(1: 3 ) * t1(1: 3 , j) )
                PPP(2, i, j) = sum(t(i, 1: 3 ) * dm(1: 3 ) * t1(1: 3 , j) )
             enddo
          enddo


          ! adding of the adiabatic terms
          ! TODO, rlen does not belong into this term, but would pop up later because of the method, right? 
          val = - 0.5 * 9.81 * h * rlen

          
          !TODO HOW WILL THIS TERM LOOK LIKE?
          PPP(1:2, 2, 1) = PPP(1:2, 2, 1) + val


          ! limiting for the vanishing density
          PPP(1:2, 1: 3 , 1: 3 ) = PPP(1:2, 1: 3 , 1: 3 ) * fac


          ! final transformation Q^{-1}*A1^{+-}*Q
          Ppm(ie,1, 1: 3 , 1: 3 ) = &
               matmul(q1(1: 3 , 1: 3 ), matmul(PPP(1, 1: 3 , 1: 3 ), q(1: 3 , 1: 3 ) ) )

          Ppm(ie,2, 1: 3 , 1: 3 ) = &
               matmul(q1(1: 3 , 1: 3 ), matmul(PPP(2, 1: 3 , 1: 3 ), q(1: 3 , 1: 3 ) ) )


       endif ! if(h <= 0.) 
    enddo

    deallocate( t, t1, q, q1, dp, dm, nv, PPP) 


  end subroutine Set_Ppm_swe

  !> compute matrices 4x4 K_sk, s,k=1,2 for N.-S. equations
  !> in integ nodes
  subroutine Set_K_sk_swe(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, K_sk, xi)
    integer, intent(in) :: ndimL, nbDim, iRe, Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes, not USED
    real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    !real, intent(in) :: Re_1                     ! inverse of Reynolds number
    real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: K_sk
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates

    K_sk(1:Qdof,1:nbDim,1:nbDim,1:ndimL,1:ndimL) = 0.

  end subroutine Set_K_sk_swe


  !> compute viscous fluxes R_s, s=1,2 for N.-S. equations
  !> in integ nodes
  subroutine Set_R_s_swe(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
    integer, intent(in) :: Qdof, nbDim, iRe, ndimL
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    !real, intent(in) :: Re_1                     ! inverse of Reynolds number
    real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates

    R_s(:, :, :) = 0.

  end subroutine Set_R_s_swe



  !> compute reactive terms S in integ nodes
  subroutine Set_S_swe(ndimL, nbDim, Qdof, xi, w, Dw, S)
    integer, intent(in) :: ndimL, nbDim, Qdof
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi!optimal velocity,SOLUTION OF THE EIKONAL EQUATION
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL), intent(inout) :: S
!    real, dimension(1:ndimL), intent(in) :: 
    integer :: i,j,k
    real :: V_eta, h, rlen
    real :: mu(1:2),gradZ(1:2)
    real :: h_minimal_allowed , h_actual, penalty

    !TODO look at state%model%Pr
    h_minimal_allowed = state%model%Pr

    !Find out bottom topology
    gradZ =BottomTopology()   

    !S(element, 1-3 slozka output)
    S(1:Qdof, :) = 0.
    !return

    do k=1,Qdof
       h = w(k, 1)

       h_actual = h
       
       !if(eta_actual < eta_minimal_allowed)  then
       !   if(eta_actual < 0.75*eta_minimal_allowed) then
       !      eta_actual = 0.
       !   else
       !      eta_actual = 4 * eta_actual - 3*eta_minimal_allowed
       !   endif
       !
       !endif
       !write(10,*) eta, eta_actual 
       
       !TODO Shouldn't there be a 'then' instead of '&'?
       if(h_actual > 0.) & 
            !S(k, 2:3) = (eta_actual* xi(k, 1:2) - w(k, 2:3) ) / state%model%tau  
            !TODO plus or minus sign? -- probably + and then changed into - by the line after enddo.
            !TODO divide by state%model%tau? Is it something which would pop up later in the method?
            S(k, 2:3) = (-9.81*h_actual*gradZ(1:2) )

       
    enddo
    
    ! the reactive term is on the right-hand side
    S(1:Qdof, 1:ndimL) = -S(1:Qdof, 1:ndimL)

    

  end subroutine Set_S_swe

  !> compute derivative of the reactive terms S in integ nodes
  subroutine Set_DS_swe(ndimL, nbDim, Qdof, xi, w, Dw, DS)
    integer, intent(in) :: ndimL, nbDim, Qdof
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:ndimL), intent(inout) :: DS
!    real, dimension(1:2), intent(in) :: 
    integer :: i,j,k
    real :: V_eta, eta, rlen, penalty
    real :: mu(1:2),gradZ(1:2)
    real :: h_minimal_allowed, eta_actual

    h_minimal_allowed = state%model%Pr


    DS(1:Qdof, :,:) = 0.
    !return
    
    !Find out bottom topology
    gradZ =BottomTopology()  

    do k=1,Qdof
       eta = w(k, 1)

       eta_actual = eta

       !if(eta_actual < eta_minimal_allowed)  then
       !   if(eta_actual < 0.75*eta_minimal_allowed) then
       !      eta_actual = 0.
       !   else
       !      eta_actual = 4 * eta_actual - 3*eta_minimal_allowed
       !   endif
       !endif
       !write(11,*) eta, eta_actual 

       if(eta_actual > 0.) then
          !TODO plus or minus sign?
          !TODO divide by state%model%tau?
          DS(k, 2:3, 1) = -9.81*gradZ(1:2)
          
       endif


    enddo

    ! the reactive term is on the right-hand side
    DS(1:Qdof, 1:ndimL, 1:ndimL)  = -DS(1:Qdof, 1:ndimL, 1:ndimL) 
  end subroutine Set_DS_swe



  function BottomTopology()  !Doplnit a dopsat ???
!    real, dimension(1:2), intent(inout) :: 
     real :: BottomTopology(1:2)
     integer :: i
   
     do i=1,2
       BottomTopology(i) = 0
     enddo
  end function BottomTopology

  ! function Barotropic_pressure(ndimL, U)
  !   real :: Barotropic_pressure
  !   integer, intent(in) :: ndimL
  !   real, dimension(1:ndimL), intent(in) :: U

  !   Barotropic_pressure = state%model%p0 * U(1)**state%model%kappa 

  ! end function Barotropic_pressure

  !> setting of BC using the solution of linearized Riemann problem
  subroutine SWE_BC_LRP(Qdof, ndimL, wi, wD, n, press_extrap, xc)
    integer, intent(in) :: Qdof, ndimL
    real, dimension(1:Qdof,1:ndimL), intent(in) :: wi
    real, dimension(1:Qdof,1:ndimL), intent(inout) ::  wD
    real, dimension(1:nbDim), intent(in) :: n, xc
    real,  intent(in) ::  press_extrap
    real, dimension(:,:), allocatable :: TT, TT1, QQ, QQ1  !,  DD, DDp, DDm
    real, dimension(:), allocatable   :: alpha, beta, omega, nv, dp, dm
    real, dimension(:), allocatable   :: qi, qD, qBC

    integer :: ie,j, iprint
    real :: kappa, kappa1, p0, h, aa
    real :: size, pp, u, v, h_minimal_allowed

    if(ndimL /= 3) then
       print*,'Bad implementation of SWE_BC_LRP'
       stop
    endif


    allocate( TT(1: 3 ,1: 3 ), TT1(1: 3 ,1: 3 ), QQ(1: 3 ,1: 3 ), QQ1(1: 3 ,1: 3 ) )
    allocate(alpha(1:3), beta(1:3), omega(1:3) )
    allocate( dp(1: 3 ), dm(1: 3 ), nv(1:nbDim))
    allocate( qi(1:3), qD(1:3), qBC(1:3)  )

    kappa = state%model%kappa  
    kappa1 = state%model%kappa1  
    p0 =  state%model%p0 

    size = sqrt(dot_product(n(1:nbDim), n(1:nbDim) ))
    nv(1:nbDim) = n(1:nbDim)/size


    !transformation matrixes
    QQ(1:3,1:3) = 0.
    QQ(1,1) = 1.
    QQ(2,2) =  nv(1)
    QQ(2,3) = nv(2)
    QQ(3,2) = -nv(2)
    QQ(3,3) = nv(1)
    
    QQ1(1:3,1:3 ) = QQ(1:3,1:3)
    QQ1(2,2) =  nv(1)
    QQ1(2,3) = -nv(2)
    QQ1(3,2) =  nv(2)
    QQ1(3,3) =  nv(1)

    iprint = 0
    h_minimal_allowed = state%model%Pr
    do ie=1,Qdof

       ! physical quantities
       h = wi(ie, 1)
       if(h > 0.) then

          ! h = max(h, h_minimal_allowed)
          ! if(h <=  h_minimal_allowed .and. iprint == 0 ) then
          !    !write(*,'(a40, 8es12.4)') 'Density <=0,  Set_BC_LRP, xc, wi:', xc(:), wi(ie, :)
          !    iprint = iprint + 1
             
          !    open( 94, file='bad_elems', status='UNKNOWN',position='append')
          !    write(94, *) xc(:), state%time%iter
          !    close(94)
          !    !stop
          ! endif
       
          ! rotation forward
          qi(1:3) = matmul(QQ(1:3, 1:3), wi(ie, 1:3) )
          qD(1:3) = matmul(QQ(1:3, 1:3), wD(ie, 1:3) )

          !transformation 
          u = qi(2) / qi(1) 
          v = qi(3) / qi(1)
          pp =  p0 *  h**kappa  !Bude potreba upravit.
          !aa = (kappa * pp / h)**0.5   ! too slow for kappa * pp / h = 1 ????
          aa = sqrt(kappa * pp / h)


          ! eigenvalues
          dp(1) = (u - aa) 
          dp(2) =    u     
          dp(3) = (u + aa) 

          !if(abs( xc(2) - 5) < 0.1 .and. dp(1) > 0.) &
          !     write(*,'(a20, 15es12.4)') 'h,u,v,pp,aa:',h,u,v,pp,aa,kappa, dp(1:3)
          !write(*,'(a20, 15es12.4)') 'wi::',wi(ie, 1:3)
          !write(*,'(a20, 15es12.4)') 'wD::',wD(ie, 1:3)

          ! matrix with eigenvectors
          TT(1, 1) = 1.     ;   TT(1, 2) = 0.  ;     TT(1, 3) =   1.   ;
          TT(2, 1) = u- aa  ;   TT(2, 2) = 0.  ;     TT(2, 3) =   u+aa ;
          TT(3, 1) = v      ;   TT(3, 2) = 1.  ;     TT(3, 3) =   v    


          ! its inverse
          TT1(1, 1) = aa + u    ;   TT1(1, 2) = -1. ;     TT1(1, 3) =   0.   ;   
          TT1(2, 1) = -2.*aa*v  ;   TT1(2, 2) =  0. ;     TT1(2, 3) = 2.*aa  ;
          TT1(3, 1) = aa - u    ;   TT1(3, 2) =  1. ;     TT1(3, 3) =  0.;   

          TT1(1:3, 1:3) = TT1(1:3, 1:3) / (2 * aa)

          ! transformation
          alpha(1:3) = matmul(TT1(1:3, 1:3) , qi(1:3) )
          beta(1:3) = matmul(TT1(1:3, 1:3) , qD(1:3) )

          do j=1,ndimL
             if(dp(j) >= 0) then
                omega(j) = alpha(j)
             else
                omega(j) = beta(j)
             endif
          enddo

          ! reverse transformation
          qBC(1:3) = matmul( TT(1:3, 1:3), omega(1:3) )

          ! backward rotation
          wD(ie, 1:3) = matmul(QQ1(1:3, 1:3), qBC(1:3) )

          !write(*,'(a20, 15es12.4)') 'wD ... ::',wD(ie, 1:3)

          !write(*,'(a20, 15es12.4)') 'qi::',qi( 1:3), alpha(1:3)
          !write(*,'(a20, 15es12.4)') 'qD::',qD( 1:3), beta(1:3) 
          !write(*,'(a20, 15es12.4)') 'qBC::',qBC( 1:3), omega(1:3)
          !print*

       else
          wD(ie, 1:3) = 0.
          !h = max(h, h_minimal_allowed)
          !if(h <=  h_minimal_allowed .and. iprint == 0 ) then
          !   !write(*,'(a40, 8es12.4)') 'Density <=0,  Set_BC_LRP, xc, wi:', xc(:), wi(ie, :)
          !   iprint = iprint + 1

          !   open( 94, file='bad_elems', status='UNKNOWN',position='append')
          !   write(94, *) xc(:), state%time%iter
          !   close(94)
          !   !stop
          !endif

       endif
    end do ! ie =1,Qdof


    deallocate(TT, TT1, QQ, QQ1, alpha, beta, omega, qi, qD, qBC, dp, dm)

  end subroutine SWE_BC_LRP

  !> setting of BC using the extrapolation based on the number of characteristics
  subroutine SWE_BC_Extrapolated(Qdof, ndimL, wi, wD, n, press_extrap, xc)
    integer, intent(in) :: Qdof, ndimL
    real, dimension(1:Qdof,1:ndimL), intent(in) :: wi
    real, dimension(1:Qdof,1:ndimL), intent(inout) ::  wD
    real, dimension(1:nbDim), intent(in) :: n, xc
    real,  intent(in) ::  press_extrap
    integer :: i,j
    real :: kappa, kappa1, p0, h
    real :: size, nn(2), vn, p, c, pD

    if(ndimL /= 3) then
       print*,'Bad implementation of SWE_BC_Characteristic'
       stop
    endif

    kappa = state%model%kappa  
    kappa1 = state%model%kappa1  
    p0 =  state%model%p0 

    size = (dot_product(n, n))**0.5
    nn(1:nbDim) = n(1:nbDim)/size

    do i=1, Qdof
       h = wi(i,1)
       vn = dot_product(nn(1:nbDim), wi(i,2:3) )/ h

       p = p0 *  h**kappa
       c = (kappa * p / h)**0.5

       !write(*,'(a2,7es10.3)' ) '>>',wi(i,1:ndim), p
       if(vn .lt. 0) then
          if(vn .lt. -c) then
             ! supersonic inlet
             ! no action
          else
             ! subsonic inlet
             wD(i,1) = h 
          endif
       else
          if(vn .gt. c) then
             ! supersonic outlet
             wD(i,1:ndimL) = wi(i,1:ndimL)
          else
             ! subsonic outlet
             !!!pD = state%model%kappa1*(wD(i,4)-dot_product(wD(i,2:3),wD(i,2:3))/wD(i,1)/2)
             pD = p0 *  wD(i,1)**kappa

             ! NEW correction in te sence of the mean value of pressure
             if(press_extrap > 0.)  then
                !write(22,'(8es12.4)' ) &
                !     xc(1:nbDim),pD, p, press_extrap, pD + p - press_extrap,vn

                pD = pD + p - press_extrap
             endif

             wD(i,2:3) = wi(i,2:3)
             !!!wD(i,4) = pD/state%model%kappa1 + dot_product(wi(i,2:3),wi(i,2:3))/wi(i,1)/2
             wD(i,1) = (pD / p0)**(1./kappa)

             ! do nothing
             !wD(i,1:ndimL) = wi(i,1:ndimL)

          endif
       endif
    enddo
  end subroutine SWE_BC_Extrapolated

end module modelSWE
