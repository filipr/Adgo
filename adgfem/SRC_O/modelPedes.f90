module pedes_averaging
  real, dimension(:,:), allocatable :: w_bar
  integer, dimension(:), allocatable :: limit_w_bar
  public:: pedes_w_bar_alloc
  public:: pedes_w_bar_DEalloc

contains

  subroutine pedes_w_bar_alloc(nelem )
    integer, intent(in) :: nelem
    allocate(w_bar(1:nelem, 1:3) )
    allocate(limit_w_bar(1:nelem) )
    limit_w_bar(:) = 0
  end subroutine pedes_w_bar_alloc


  subroutine pedes_w_bar_DEalloc( )
    deallocate(w_bar )
    deallocate(limit_w_bar )
  end subroutine pedes_w_bar_DEalloc

end module pedes_averaging


!> definition of 2D models which are simulated: Euler, Navier-Stokes
module modelPedes
  
  use data_mod
  use main_data
  use f_mapping
  use mesh_oper
  use blocks_integ
  use matrix_oper_int
  use pedes_averaging

  implicit none
  
  
  public:: Set_f_s_pedes
  public:: Set_A_s_pedes
  public:: Set_Ppm_pedes

  public:: Set_R_s_empty
  public:: Set_K_sk_empty

  public:: Set_S_pedes
  public:: Set_DS_pedes

contains

  !> limiting factor for pedestrian flow fluxes
  !> \f$ fac \in [0,1] \f$, 
  !> \f$ \rho <= fac_r \rho_{\min} \ \Rightarrow \ fac = 0 \ \Rightarrow\ \f$ no flux
  subroutine Set_limit_factor(rho, fac)
    real, intent(in) :: rho
    real, intent(out) :: fac
    real :: fac_r

    !fac_r = 0.125    
    fac_r = 0.25    !work reasonably
    !fac_r = 0.50

    if(rho <=  fac_r * state%model%Pr) then
       fac = 0.

    else if(rho >=  state%model%Pr) then
       fac = 1.

    else
       fac = -cos( pi *(rho - fac_r* state%model%Pr)/ ((1.-fac_r) * state%model%Pr)) /2 + 0.5
       !print*,'####', rho,  state%model%Pr, fac
    endif

    !!fac = 1.

  end subroutine Set_limit_factor

  !> evaluate the amount of penalty for gradient penalization
  !> penalty (rho_min) = 0 & derivative of penalty at rho_min = 0.
  subroutine Set_grad_penalty(rho, penalty)
    real, intent(in) :: rho
    real, intent(out) :: penalty
    real :: alpha, rho_min
    real :: max_arg , val_arg

    alpha = 100. 

    max_arg = -8.

    rho_min = state%model%Pr

    if(rho >=  rho_min) then
       penalty = 0.

    else
       
       !val_arg = alpha *(rho - rho_min)
       !if(val_arg <  max_arg ) val_arg = max_arg
       !penalty = exp(val_arg ) + val_arg - 1.


       penalty = (max(0., rho) - rho_min)**2 * alpha
        !penalty = exp(-alpha *(rho - rho_min) ) + alpha*(rho - rho_min) - 1.

    endif

    !if(penalty > 1.) print*,' Grad penalty yt74u:', rho, penalty

  end subroutine Set_grad_penalty


  !> compute the pedestrian flow fluxes \f$ f_s({\bf w}),\quad s=1,2\f$
  subroutine Set_f_s_pedes(ndimL, nbDim, Qdof, w, f_s, xi, ie)
    integer, intent(in) :: Qdof, nbDim, ndimL
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in integ nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s ! fluxes f_s in  -- " --
    real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
    integer, intent(in) :: ie  ! index of the element
    real, dimension(:), allocatable :: u,v,p
    real, dimension(:,:,:,:), allocatable :: A_s
    integer:: i
    real :: rho_minimal_allowed, fac

    rho_minimal_allowed = state%model%Pr
 
    allocate( u(1:Qdof), v(1:Qdof),  p(1:Qdof) )

    ! if( limit_w_bar(ie) == 1) then
    !    ! linearization around the mean value)
    !    allocate(A_s(1:Qdof,1:nbDim,1:ndimL,1:ndimL)) 

    !    ! setting of the matrix from the mean values
    !    call Set_A_s_pedes(ndimL, nbDim, Qdof, w, A_s(1:Qdof,1:nbDim,1:ndimL,1:ndimL), &
    !         xi(1:Qdof,1 :nbDim), ie)

    !    ! linearized matrix from mean value times the standard vector
    !    do i=1, Qdof
    !       f_s(i, 1, 1:ndimL) = matmul(A_s(i, 1, 1:ndimL, 1:ndimL), w(i, 1:ndimL) )
    !       f_s(i, 2, 1:ndimL) = matmul(A_s(i, 2, 1:ndimL, 1:ndimL), w(i, 1:ndimL) )
    !    enddo
    !    deallocate(A_s)

    ! else ! standard setting

    do i=1, Qdof

       !print*,'####', rho_minimal_allowed, i, Qdof
       if(w(i, 1) <= 0.) then
          !if(w(i, 1) <= rho_minimal_allowed) then
          f_s(i, :, :) = 0.

       else
          call Set_limit_factor(w(i,1), fac)

          u(i) = w(i, 2) / w(i, 1)
          v(i) = w(i, 3) / w(i, 1)
          p(i) = state%model%p0 *  w(i, 1)**state%model%kappa


          f_s(i, 1, 1) = w(i, 2)
          f_s(i, 1, 2) = w(i, 2) * u(i) + p(i)
          f_s(i, 1, 3) = w(i, 2) * v(i)

          f_s(i, 2, 1) = w(i, 3)
          f_s(i, 2, 2) = w(i, 3) * u(i)
          f_s(i, 2, 3) = w(i, 3) * v(i) + p(i)

          !if(w(i,1) < 1E-5) write(*,'(a8,40es12.4)') 'de f:s',w(i,1:3), u(i), v(i),p(i),f_s(i, 1, 1:3)

          f_s(i, :, :) = f_s(i, :, :) * fac
       endif


    enddo

    !!endif

    deallocate(u,v,p )

  end subroutine Set_f_s_pedes


  !> compute matrices \f$ A_s = \frac{D{\bf f_s}({\bf w})}{D{\bf w}},\quad s=1,2\f$
  !> for pedestrian flow fluxes
  subroutine Set_A_s_pedes(ndimL, nbDim, Qdof, w, A_s, xi, ie)
    integer, intent(in) :: ndimL, nbDim, Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: A_s
                                               ! matrices A_s in  -- " --
    real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
    integer, intent(in) :: ie  ! index of the element
    real, dimension(:), allocatable :: u,v,rho, pp, uv, u2,v2
    integer :: i
    real :: rho_minimal_allowed, fac

    rho_minimal_allowed = state%model%Pr


    allocate( u(1:Qdof), v(1:Qdof), rho(1:Qdof), pp(1:Qdof) )
    allocate( uv(1:Qdof), u2(1:Qdof), v2(1:Qdof) )
    
    ! ! setting of the density and velocity in integ nodes
    ! if( limit_w_bar(ie) == 1) then
    !    ! linearization around the mean value)
    !    rho(1:Qdof) = w_bar(ie, 1)
    !    u(1:Qdof) = w_bar(ie, 2) / w_bar(ie, 1)
    !    v(1:Qdof) = w_bar(ie, 3) / w_bar(ie, 1)
       
    ! else ! standard setting
       
    do i=1,Qdof
       rho(i) = w(i,1)
       !if( rho(i) > 0. ) rho(i) = max(rho(i), rho_minimal_allowed)

       if(rho(i) > 0.   ) then

          call Set_limit_factor(rho(i), fac)
          
          !if(rho(i) > rho_minimal_allowed) then
          u(i) = w(i,2)/rho(i)
          v(i) = w(i,3)/rho(i)
       else
          u(i) = 0.
          v(i) = 0.
       endif

       !enddo
       !endif

       !do i=1,Qdof

       if(rho(i) <= 0.) then
       !if(rho(i) <= rho_minimal_allowed) then
          A_s(i, :,:,:) = 0.
       else

          pp(i) =  state%model%p0 *  rho(i)**state%model%kappa1 ! IS NOT PRESSURE !!!!
          uv(i) = u(i)*v(i)
          u2(i) = u(i)*u(i)
          v2(i) = v(i)*v(i)
          
          A_s(i, 1, 1, 1) = 0.
          A_s(i, 1, 1, 2) = 1.
          A_s(i, 1, 1, 3) = 0.

          A_s(i, 1, 2, 1) = -u2(i) + pp(i)
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

          A_s(i, 2, 3, 1) = -v2(i) + pp(i)
          A_s(i, 2, 3, 2) = 0.
          A_s(i, 2, 3, 3) = 2 * v(i)

          !if(w(i,1) < 1E-5) write(*,'(a8,40es12.4)') 'de A:s',w(i,1:3), u(i),v(i),pp(i),u2(i),uv(i),v2(i)
          A_s(i, :,:,:) = A_s(i, :,:,:) * fac

       endif
    enddo

    deallocate(u,v,rho, pp, uv, u2, v2 )

  end subroutine Set_A_s_pedes



  !> compute matrices
  !> \f$ P^{\pm} = \left(\frac{D({\bf f_1}({\bf w})n_1+{\bf f_2}({\bf w})n_2}{D{\bf w}}
  !>  \right)^{\pm}\f$
  !> for the pedestrian flow equation
  subroutine Set_Ppm_pedes(ndimL, nbDim, Qdof, w, n, xi, Ppm, one_over_area, elem)
    !class(mesh), intent(in) :: grid
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
    real :: rlen, rho, u, v,pp, p0, aa, kappa, kappa1, val, fac, fac_r
    real :: rho_minimal_allowed
    integer :: ie, i, j, itest1, iprint
    real :: tt, tt1, velocity, max_velocity

    itest1 = -grid%nelem

    call cpu_time(tt)

    if(ndimL /= 3) stop 'Bad dimension (ndim) in subroutine Set_Ppm_pedes !!'

    allocate( t(1: 3 ,1: 3 ), t1(1: 3 ,1: 3 ), q(1: 3 ,1: 3 ), q1(1: 3 ,1: 3 ) )
    allocate( PPP( 1:2, 1: 3 ,1: 3 ) )

    allocate( dp(1: 3 ), dm(1: 3 ), nv(1:nbDim))

    !allocate( DDp(1: 3 , 1: 3 ), DDm(1: 3 , 1: 3 ) )

    !tt1 = tt
    !call cpu_time(tt)
    !if(elem%i == itest1) print*,'MM,',tt,' P+-1', tt - tt1


    !print*, 'Set_Ppm_pedes', Qdof

    if (.not. present(elem) .or. .not. present(one_over_area) ) &
      stop 'elem and one_over_area must be present in Set_Ppm_pedes!'

    kappa = state%model%kappa  
    kappa1 = state%model%kappa1  
    p0 =  state%model%p0 
    rho_minimal_allowed = state%model%Pr
    max_velocity = 10.

    iprint = 0

    do ie=1,Qdof

       rlen = (n(ie,1)*n(ie, 1) + n(ie, 2)*n(ie, 2) )**0.5
       nv(1:nbDim) = n(ie, 1:nbDim)/rlen

       ! physical quantities
       rho = w(ie, 1)

       if(rho <= 0. )  then
       !if(rho <= rho_minimal_allowed) then
          Ppm(ie,1:2, 1: 3 , 1: 3 ) = 0.

       else

          call Set_limit_factor(rho, fac)

          !rho = max(rho, rho_minimal_allowed)
          !if(rho <= rho_minimal_allowed .and. iprint == 0 ) then
          !   !write(*,'(a40,i5, 8es12.4)') 'Density <=0, Set_Ppm_pedes, i, xc, xi ', &
          !   !     elem%i,elem%xc(:), w(ie, :)
          !   iprint = iprint + 1
          !
          !   open( 94, file='bad_elems', status='UNKNOWN',position='append')
          !   write(94, *) elem%xc(:), state%time%iter
          !   close(94)
          !   !stop
          !endif


          !transformation 
          u = ( w(ie, 2) * nv(1) + w(ie, 3) * nv(2) ) / rho
          v = (-w(ie, 2) * nv(2) + w(ie, 3) * nv(1) ) / rho

          !u = ( w(ie, 2) * nv(1) + w(ie, 3) * nv(2) ) / max( rho, 0.5 * state%model%Pr) 
          !v = (-w(ie, 2) * nv(2) + w(ie, 3) * nv(1) ) / max( rho, 0.5 * state%model%Pr) 

          !if(rho <= rho_minimal_allowed) print*,'????', rho, u, v
          !if(rho <= 0.5 * state%model%Pr )  then
          !if(rho <= state%model%Pr .and. elem%i == 2214 .and. (u > 5. .or. v > 5.)) then
          !      write(*,'(a8, i5, 3(4es12.4, a3))') &
          !        '????', elem%i, state%model%Pr, rho, u, v, '||', w(ie, 1:3)
          !   stop "de3u39jiw"
          !endif


          !velocity = sqrt(u*u + v*v) 
          !if( velocity > max_velocity) then
          !   u = u / velocity * max_velocity
          !   v = v / velocity * max_velocity
          !endif

          pp =  p0 *  rho**kappa
          !aa = (kappa * pp / rho)**0.5   ! too slow for kappa * pp / rho = 1 ????
          aa = sqrt(kappa * pp / rho)

          !write(*,'(a20, 15es12.4)') 'rho,u,v,pp,aa:',rho,u,v,pp,aa,kappa

          ! eigenvalues
          dp(1) = (u - aa) * rlen
          dp(2) =    u     * rlen
          dp(3) = (u + aa) * rlen

          !if(w(ie,1) < 1E-5) write(*,'(a8,40es12.4)') 'de P+-',w(ie,1:3), u,v, aa, dp(:)
         

          !calculating the negative (dm) and the positive(dp)  parts of eigenvalues
          do i=1, 3 
             dm(i)=0.
             if(dp(i) < 0.)then
                dm(i)=dp(i)
                dp(i)=0.
             endif
          enddo

          elem%max_eigenvals = max( abs(dp(1)), abs(dp(3) ) ) *one_over_area
          if(rho > 0.01) &
               state%max_eigenvals = max(state%max_eigenvals, elem%max_eigenvals  )

          !tt1 = tt
          !call cpu_time(tt)
          !if(elem%i == itest1) print*,'MM,',tt,' P+-2', tt - tt1, ie



          ! matrix with eigenvectors
          t(1, 1) = 1.     
          t(1, 2) = 0.  
          t(1, 3) =   1.   

          t(2, 1) = u- aa  
          t(2, 2) = 0.  
          t(2, 3) =   u+aa 

          t(3, 1) = v      
          t(3, 2) = 1.  
          t(3, 3) =   v    


          ! its inverse
          t1(1, 1) = aa + u    
          t1(1, 2) = -1.  
          t1(1, 3) =   0.   

          t1(2, 1) = -2.*aa*v  
          t1(2, 2) =  0.  
          t1(2, 3) = 2.*aa  

          t1(3, 1) = aa - u    
          t1(3, 2) =  1.  
          t1(3, 3) =  0.    


          t1(1:3, 1:3) = t1(1:3, 1:3) / (2 * aa)


          !tt1 = tt
          !call cpu_time(tt)
          !if(elem%i == itest1) print*,'MM,',tt,' P+-b', tt - tt1, ie

          !transformation matrixes
          q(1:3,1:3) = 0.
          q(1,1) = 1.
          q(2,2) =  nv(1)
          q(2,3) = nv(2)
          q(3,2) = -nv(2)
          q(3,3) = nv(1)

          q1(1:3,1:3 ) = q(1:3,1:3)
          q1(2,2) =  nv(1)
          q1(2,3) = -nv(2)
          q1(3,2) =  nv(2)
          q1(3,3) =  nv(1)


          !write(*,'(a10,40es12.4)') 'ede3ws',matmul(t, t1)
          !write(*,'(a10,40es12.4)') 'ede3ws',matmul(t1, t)


          !  multiplication of matrices: pp=t*dp*t1, where dp is diagonal

          !write(*,'(a20, 6es12.4)') 'Eigenvals+=',dp(1: 3 )
          !write(*,'(a20, 6es12.4)') 'Eigenvals-=',dm(1: 3 )
          !print*,'...'

          !DDm(:,:) = 0.
          !DDp(:,:) = 0.
          !do i=1, 3 
          !   DDm(i,i) = dm(i)
          !   DDp(i,i) = dp(i)
          !enddo

          !tt1 = tt
          !call cpu_time(tt)
          !if(elem%i == itest1) print*,'MM,',tt,' P+-c', tt - tt1, ie


          do i=1, 3 
             do j=1, 3 
                PPP(1, i, j) = sum(t(i, 1: 3 ) * dp(1: 3 ) * t1(1: 3 , j) )
                PPP(2, i, j) = sum(t(i, 1: 3 ) * dm(1: 3 ) * t1(1: 3 , j) )
             enddo
          enddo

          !PPP(1, 1: 3 , 1: 3 ) = matmul(t(1: 3 , 1: 3 ),  matmul(DDp(1: 3 , 1: 3 ), t1(1: 3 , 1: 3 ) ) )

          !PPP(2, 1: 3 , 1: 3 ) = matmul(t(1: 3 , 1: 3 ),  matmul(DDm(1: 3 , 1: 3 ), t1(1: 3 , 1: 3 ) ) )

          !  do i=1,ndim
          !     write(*,'(a3,i5,3es12.4,a2,3es12.4)') 'DD',i, DDp( i, 1:ndim),'|', DDm( i, 1:ndim)
          !  enddo

          !  print*
          !  do i=1,ndim
          !     write(*,'(a3,i5,3es12.4,a2,3es12.4)') 'tt',i, t( i, 1:ndim),'|', t1( i, 1:ndim)
          !  enddo

          !  print*
          ! do i=1,ndim
          !     write(*,'(a3,i5,3es12.4,a2,3es12.4)') 'QQ',i, q( i, 1:ndim),'|', q1( i, 1:ndim)
          !  enddo

          !  print*

          !  do i=1,ndim
          !     write(*,'(a3,i5,3es12.4,a2,3es12.4)') 'PP',i, PPP(1,i, 1:ndim),'|', PPP(2, i, 1:ndim)
          ! enddo

          !tt1 = tt
          !call cpu_time(tt)
          !if(elem%i == itest1) print*,'MM,',tt,' P+-3', tt - tt1, ie


          ! adding of the adiabatic terms
          val = - 0.5 *kappa1 * p0 * rho**kappa1 *rlen

          PPP(1:2, 2, 1) = PPP(1:2, 2, 1) + val


          ! limiting for the vanishing density
          PPP(1:2, 1: 3 , 1: 3 ) = PPP(1:2, 1: 3 , 1: 3 ) * fac

          !print*
          ! do i=1,ndim
          !   write(*,'(a3,i5,3es12.4,a2,3es12.4)') 'PPd',i, PPP(1, i, 1:ndim),'|', PPP(2, i, 1:ndim)
          !enddo

          ! final transformation
          Ppm(ie,1, 1: 3 , 1: 3 ) = &
               matmul(q1(1: 3 , 1: 3 ), matmul(PPP(1, 1: 3 , 1: 3 ), q(1: 3 , 1: 3 ) ) )

          Ppm(ie,2, 1: 3 , 1: 3 ) = &
               matmul(q1(1: 3 , 1: 3 ), matmul(PPP(2, 1: 3 , 1: 3 ), q(1: 3 , 1: 3 ) ) )


          !tt1 = tt
          !call cpu_time(tt)
          !if(elem%i == itest1) print*,'MM,',tt,' P+-4', tt - tt1, ie


          !print*,'...,'
       endif ! if(rho <= 0.) 
    enddo

    deallocate( t, t1, q, q1, dp, dm, nv, PPP) 
    !deallocate(DDm, DDp)
    
    !tt1 = tt
    !call cpu_time(tt)
    !if(elem%i == itest1) print*,'MM,',tt,' P+-X', tt - tt1


  end subroutine Set_Ppm_pedes

  !> compute matrices 4x4 K_sk, s,k=1,2 for N.-S. equations
  !> in integ nodes
  subroutine Set_K_sk_empty(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, K_sk, xi)
    integer, intent(in) :: ndimL, nbDim, iRe, Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes, not USED
    real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    !real, intent(in) :: Re_1                     ! inverse of Reynolds number
    real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: K_sk
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates

    K_sk(1:Qdof,1:nbDim,1:nbDim,1:ndimL,1:ndimL) = 0.

  end subroutine Set_K_sk_empty


  !> compute viscous fluxes R_s, s=1,2 for N.-S. equations
  !> in integ nodes
  subroutine Set_R_s_empty(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
    integer, intent(in) :: Qdof, nbDim, iRe, ndimL
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    !real, intent(in) :: Re_1                     ! inverse of Reynolds number
    real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates

    R_s(:, :, :) = 0.

  end subroutine Set_R_s_empty



  !> compute reactive terms S in integ nodes
  subroutine Set_S_pedes(ndimL, nbDim, Qdof, xi, w, Dw, S)
    integer, intent(in) :: ndimL, nbDim, Qdof
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi!optimal velocity,SOLUTION OF THE EIKONAL EQUATION
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL), intent(inout) :: S
    integer :: i,j,k
    real :: V_rho, rho, rlen
    real :: mu(1:2)
    real :: rho_minimal_allowed , rho_actual, penalty

    rho_minimal_allowed = state%model%Pr

    S(1:Qdof, :) = 0.
    !return

    do k=1,Qdof
       rho = w(k, 1)

       rho_actual = rho
       
       !if(rho_actual < rho_minimal_allowed)  then
       !   if(rho_actual < 0.75*rho_minimal_allowed) then
       !      rho_actual = 0.
       !   else
       !      rho_actual = 4 * rho_actual - 3*rho_minimal_allowed
       !   endif
       !
       !endif
       !write(10,*) rho, rho_actual 
       
       if(rho_actual > 0.) & 
            S(k, 2:3) = (rho_actual* xi(k, 1:2) - w(k, 2:3) ) / state%model%tau 

       
    enddo
    
    ! the reactive term is on the right-hand side
    S(1:Qdof, 1:ndimL) = -S(1:Qdof, 1:ndimL)

    

  end subroutine Set_S_pedes

  !> compute derivative of the reactive terms S in integ nodes
  subroutine Set_DS_pedes(ndimL, nbDim, Qdof, xi, w, Dw, DS)
    integer, intent(in) :: ndimL, nbDim, Qdof
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:ndimL), intent(inout) :: DS
    integer :: i,j,k
    real :: V_rho, rho, rlen, penalty
    real :: mu(1:2)
    real :: rho_minimal_allowed, rho_actual

    rho_minimal_allowed = state%model%Pr


    DS(1:Qdof, :,:) = 0.
    !return

    do k=1,Qdof
       rho = w(k, 1)

       rho_actual = rho

       !if(rho_actual < rho_minimal_allowed)  then
       !   if(rho_actual < 0.75*rho_minimal_allowed) then
       !      rho_actual = 0.
       !   else
       !      rho_actual = 4 * rho_actual - 3*rho_minimal_allowed
       !   endif
       !endif
       !write(11,*) rho, rho_actual 

       if(rho_actual > 0.) then
          !DS(k, 2:3, 1) = rho_actual * xi(k,1:2)  / state%model%tau  ! ERROR !!!
          DS(k, 2:3, 1) = xi(k,1:2)  / state%model%tau 
          DS(k, 2, 2) = - 1. / state%model%tau
          DS(k, 3, 3) = - 1. / state%model%tau
       endif


    enddo

    ! the reactive term is on the right-hand side
    DS(1:Qdof, 1:ndimL, 1:ndimL)  = -DS(1:Qdof, 1:ndimL, 1:ndimL) 
  end subroutine Set_DS_pedes





  function Barotropic_pressure(ndimL, U)
    real :: Barotropic_pressure
    integer, intent(in) :: ndimL
    real, dimension(1:ndimL), intent(in) :: U

    Barotropic_pressure = state%model%p0 * U(1)**state%model%kappa 

  end function Barotropic_pressure

  !> setting of BC using the solution of linearized Riemann problem
  subroutine Pedestrian_BC_LRP(Qdof, ndimL, wi, wD, n, press_extrap, xc)
    integer, intent(in) :: Qdof, ndimL
    real, dimension(1:Qdof,1:ndimL), intent(in) :: wi
    real, dimension(1:Qdof,1:ndimL), intent(inout) ::  wD
    real, dimension(1:nbDim), intent(in) :: n, xc
    real,  intent(in) ::  press_extrap
    real, dimension(:,:), allocatable :: TT, TT1, QQ, QQ1  !,  DD, DDp, DDm
    real, dimension(:), allocatable   :: alpha, beta, omega, nv, dp, dm
    real, dimension(:), allocatable   :: qi, qD, qBC

    integer :: ie,j, iprint
    real :: kappa, kappa1, p0, rho, aa
    real :: size, pp, u, v, rho_minimal_allowed

    if(ndimL /= 3) then
       print*,'Bad implementation of Pedestrian_BC_LRP'
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
    rho_minimal_allowed = state%model%Pr
    do ie=1,Qdof

       ! physical quantities
       rho = wi(ie, 1)
       if(rho > 0.) then

          ! rho = max(rho, rho_minimal_allowed)
          ! if(rho <=  rho_minimal_allowed .and. iprint == 0 ) then
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
          pp =  p0 *  rho**kappa
          !aa = (kappa * pp / rho)**0.5   ! too slow for kappa * pp / rho = 1 ????
          aa = sqrt(kappa * pp / rho)


          ! eigenvalues
          dp(1) = (u - aa) 
          dp(2) =    u     
          dp(3) = (u + aa) 

          !if(abs( xc(2) - 5) < 0.1 .and. dp(1) > 0.) &
          !     write(*,'(a20, 15es12.4)') 'rho,u,v,pp,aa:',rho,u,v,pp,aa,kappa, dp(1:3)
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
          !rho = max(rho, rho_minimal_allowed)
          !if(rho <=  rho_minimal_allowed .and. iprint == 0 ) then
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

  end subroutine Pedestrian_BC_LRP

  !> setting of BC using the extrapolation based on the number of characteristics
  subroutine Pedestrian_BC_Extrapolated(Qdof, ndimL, wi, wD, n, press_extrap, xc)
    integer, intent(in) :: Qdof, ndimL
    real, dimension(1:Qdof,1:ndimL), intent(in) :: wi
    real, dimension(1:Qdof,1:ndimL), intent(inout) ::  wD
    real, dimension(1:nbDim), intent(in) :: n, xc
    real,  intent(in) ::  press_extrap
    integer :: i,j
    real :: kappa, kappa1, p0, rho
    real :: size, nn(2), vn, p, c, pD

    if(ndimL /= 3) then
       print*,'Bad implementation of Pedestrian_BC_Characteristic'
       stop
    endif

    kappa = state%model%kappa  
    kappa1 = state%model%kappa1  
    p0 =  state%model%p0 

    size = (dot_product(n, n))**0.5
    nn(1:nbDim) = n(1:nbDim)/size

    do i=1, Qdof
       rho = wi(i,1)
       vn = dot_product(nn(1:nbDim), wi(i,2:3) )/ rho

       p = p0 *  rho**kappa
       c = (kappa * p / rho)**0.5

       !write(*,'(a2,7es10.3)' ) '>>',wi(i,1:ndim), p
       if(vn .lt. 0) then
          if(vn .lt. -c) then
             ! supersonic inlet
             ! no action
          else
             ! subsonic inlet
             wD(i,1) = rho 
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
  end subroutine Pedestrian_BC_Extrapolated

end module modelPedes
