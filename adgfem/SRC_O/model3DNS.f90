!> definition of 3D models which are simulated: Euler, Navier-Stokes
module model3DNS

  use main_data
  use f_mapping
  use mesh_oper
  use blocks_integ

  implicit none

  public:: Set_f_s_Euler3D
  public:: Set_A_s_Euler3D
  public:: Set_Ppm_Euler3D
  public:: ComputeEigenValsVec3D
  public:: Set_R_s_NS3D
  public:: Set_K_sk_NS3D

contains
  !> compute the Euler fluxes \f$ f_s({\bf w}),\quad s=1,3\f$ 
  subroutine Set_f_s_Euler3D(ndim, Qdof, w, f_s)
    integer, intent(in) :: Qdof, ndim
    real, dimension(1:Qdof, 1:ndim), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndim), intent(inout) :: f_s
                                               ! matrices A_s in  -- " --
    real, dimension(:), allocatable :: u,v,z,p ! velocity components and pressure
    
    allocate( u(1:Qdof), v(1:Qdof), z(1:Qdof),  p(1:Qdof) )


    u(1:Qdof) = w(1:Qdof, 2) / w(1:Qdof, 1)
    v(1:Qdof) = w(1:Qdof, 3) / w(1:Qdof, 1)
    z(1:Qdof) = w(1:Qdof, 4) / w(1:Qdof, 1)
    p(1:Qdof) = state%model%kappa1 * (w(1:Qdof,5)  &
         - (u(1:Qdof)*w(1:Qdof, 2) + v(1:Qdof)*w(1:Qdof, 3) + z(1:Qdof)*w(1:Qdof, 4) )/2. )


    f_s(1:Qdof, 1, 1) = w(1:Qdof, 2)
    f_s(1:Qdof, 1, 2) = w(1:Qdof, 2)*u(1:Qdof) + p(1:Qdof)
    f_s(1:Qdof, 1, 3) = w(1:Qdof, 2)*v(1:Qdof)
    f_s(1:Qdof, 1, 4) = w(1:Qdof, 2)*z(1:Qdof)
    f_s(1:Qdof, 1, 5) = (w(1:Qdof,5) + p(1:Qdof) ) * u(1:Qdof)
    
    f_s(1:Qdof, 2, 1) = w(1:Qdof,3)
    f_s(1:Qdof, 2, 2) = w(1:Qdof, 3)*u(1:Qdof) 
    f_s(1:Qdof, 2, 3) = w(1:Qdof, 3)*v(1:Qdof) + p(1:Qdof)
    f_s(1:Qdof, 2, 4) = w(1:Qdof, 3)*z(1:Qdof) 
    f_s(1:Qdof, 2, 5) = (w(1:Qdof,5) + p(1:Qdof) ) * v(1:Qdof)

    f_s(1:Qdof, 3, 1) = w(1:Qdof,3)
    f_s(1:Qdof, 3, 2) = w(1:Qdof, 4)*u(1:Qdof) 
    f_s(1:Qdof, 3, 3) = w(1:Qdof, 4)*v(1:Qdof)
    f_s(1:Qdof, 3, 4) = w(1:Qdof, 4)*z(1:Qdof) + p(1:Qdof) 
    f_s(1:Qdof, 3, 5) = (w(1:Qdof,5) + p(1:Qdof) ) * z(1:Qdof)


    deallocate(u,v,z,p )

  end subroutine Set_f_s_Euler3D

  !> compute matrices \f$ A_s = \frac{D{\bf f_s}({\bf w})}{D{\bf w}},\quad s=1,2,3\f$ 
  !> for Euler fluxes 
  subroutine Set_A_s_Euler3D(ndim, Qdof, stateW, A_s)
    integer, intent(in) :: Qdof, ndim
    real, dimension(1:Qdof, 1:ndim), intent(in):: stateW !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndim,1:ndim), intent(inout) :: A_s
    ! matrices A_s in  -- " --
    real :: kappa,kappa1,kappa3
    real, dimension(:), allocatable :: u,v,w,rho,e,uv,uw,vw,kappaE,u2,v2,w2,u2v2w2
    integer :: i
    

    allocate( u(1:Qdof), v(1:Qdof), w(1:Qdof), rho(1:Qdof), e(1:Qdof) )
    allocate( uv(1:Qdof),uw(1:Qdof),vw(1:Qdof), kappaE(1:Qdof), &
         u2(1:Qdof), v2(1:Qdof), w2(1:Qdof), u2v2w2(1:Qdof) )
    
    kappa  = state%model%kappa
    kappa1 = state%model%kappa1
    kappa3 = kappa - 3
    
    rho(1:Qdof)= stateW(1:Qdof,1)
    u(1:Qdof)  = stateW(1:Qdof,2)/rho(1:Qdof)
    v(1:Qdof)  = stateW(1:Qdof,3)/rho(1:Qdof)
    w(1:Qdof)  = stateW(1:Qdof,4)/rho(1:Qdof)
    e(1:Qdof)  = stateW(1:Qdof,5)
    uv(1:Qdof) = u(1:Qdof)*v(1:Qdof)
    uw(1:Qdof) = u(1:Qdof)*w(1:Qdof)
    vw(1:Qdof) = v(1:Qdof)*w(1:Qdof)
    u2(1:Qdof) = u(1:Qdof)*u(1:Qdof)
    v2(1:Qdof) = v(1:Qdof)*v(1:Qdof)
    w2(1:Qdof) = w(1:Qdof)*w(1:Qdof)
    kappaE(1:Qdof) = kappa*e(1:Qdof)
    
    u2v2w2(1:Qdof) = u2(1:Qdof) + v2(1:Qdof) + w2(1:Qdof)
    
    A_s(1:Qdof,1,1,1) = 0
    A_s(1:Qdof,1,1,2) = 1.0
    A_s(1:Qdof,1,1,3) = 0
    A_s(1:Qdof,1,1,4) = 0
    A_s(1:Qdof,1,1,5) = 0
    
    A_s(1:Qdof,1,2,1) =  kappa1*u2v2w2(1:Qdof)/2.0 - u2(1:Qdof)
    A_s(1:Qdof,1,2,2) = -kappa3*u(1:Qdof)
    A_s(1:Qdof,1,2,3) = -kappa1*v(1:Qdof)
    A_s(1:Qdof,1,2,4) = -kappa1*w(1:Qdof)
    A_s(1:Qdof,1,2,5) = kappa1
    
    A_s(1:Qdof,1,3,1) = -uv(1:Qdof)
    A_s(1:Qdof,1,3,2) = v(1:Qdof)
    A_s(1:Qdof,1,3,3) = u(1:Qdof)
    A_s(1:Qdof,1,3,4) = 0
    A_s(1:Qdof,1,3,5) = 0
    
    A_s(1:Qdof,1,4,1) = -uw(1:Qdof)
    A_s(1:Qdof,1,4,2) = w(1:Qdof)
    A_s(1:Qdof,1,4,3) = 0
    A_s(1:Qdof,1,4,4) = u(1:Qdof)
    A_s(1:Qdof,1,4,5) = 0
    
    A_s(1:Qdof,1,5,1) = u(1:Qdof)*(kappa1*u2v2w2(1:Qdof) - kappaE(1:Qdof)/rho(1:Qdof))
    A_s(1:Qdof,1,5,2) = kappaE(1:Qdof)/rho(1:Qdof) - kappa1*(u2(1:Qdof) + u2v2w2(1:Qdof)/2.0)
    A_s(1:Qdof,1,5,3) = -kappa1*uv(1:Qdof)
    A_s(1:Qdof,1,5,4) = -kappa1*uw(1:Qdof)
    A_s(1:Qdof,1,5,5) = kappa*u(1:Qdof)
    
    A_s(1:Qdof,2,1,1) = 0
    A_s(1:Qdof,2,1,2) = 0
    A_s(1:Qdof,2,1,3) = 1.0
    A_s(1:Qdof,2,1,4) = 0
    A_s(1:Qdof,2,1,5) = 0
    
    A_s(1:Qdof,2,2,1) = -uv(1:Qdof)
    A_s(1:Qdof,2,2,2) = v(1:Qdof)
    A_s(1:Qdof,2,2,3) = u(1:Qdof)
    A_s(1:Qdof,2,2,4) = 0
    A_s(1:Qdof,2,2,5) = 0
    
    A_s(1:Qdof,2,3,1) =  kappa1*u2v2w2(1:Qdof)/2.0 - v2(1:Qdof)
    A_s(1:Qdof,2,3,2) = -kappa1*u(1:Qdof)
    A_s(1:Qdof,2,3,3) = -kappa3*v(1:Qdof)
    A_s(1:Qdof,2,3,4) = -kappa1*w(1:Qdof)
    A_s(1:Qdof,2,3,5) = kappa1
    
    A_s(1:Qdof,2,4,1) = -vw(1:Qdof)
    A_s(1:Qdof,2,4,2) = 0
    A_s(1:Qdof,2,4,3) = w(1:Qdof)
    A_s(1:Qdof,2,4,4) = v(1:Qdof)
    A_s(1:Qdof,2,4,5) = 0
    
    A_s(1:Qdof,2,5,1) = v(1:Qdof)*(kappa1*u2v2w2(1:Qdof) - kappaE(1:Qdof)/rho(1:Qdof))
    A_s(1:Qdof,2,5,2) = -kappa1*uv(1:Qdof)
    A_s(1:Qdof,2,5,3) = kappaE(1:Qdof)/rho(1:Qdof) - kappa1*(v2(1:Qdof) + u2v2w2(1:Qdof)/2.0)
    A_s(1:Qdof,2,5,4) = -kappa1*vw(1:Qdof)
    A_s(1:Qdof,2,5,5) = kappa*v(1:Qdof)
    
    A_s(1:Qdof,3,1,1) = 0
    A_s(1:Qdof,3,1,2) = 0
    A_s(1:Qdof,3,1,3) = 0
    A_s(1:Qdof,3,1,4) = 1.0
    A_s(1:Qdof,3,1,5) = 0
    
    A_s(1:Qdof,3,2,1) = -uw(1:Qdof)
    A_s(1:Qdof,3,2,2) = w(1:Qdof)
    A_s(1:Qdof,3,2,3) = 0
    A_s(1:Qdof,3,2,4) = u(1:Qdof)
    A_s(1:Qdof,3,2,5) = 0
    
    A_s(1:Qdof,3,3,1) = -vw(1:Qdof)
    A_s(1:Qdof,3,3,2) = 0
    A_s(1:Qdof,3,3,3) = w(1:Qdof)
    A_s(1:Qdof,3,3,4) = v(1:Qdof)
    A_s(1:Qdof,3,3,5) = 0
    
    A_s(1:Qdof,3,4,1) =  kappa1*u2v2w2(1:Qdof)/2.0 - w2(1:Qdof)
    A_s(1:Qdof,3,4,2) = -kappa1*u(1:Qdof)
    A_s(1:Qdof,3,4,3) = -kappa1*v(1:Qdof)
    A_s(1:Qdof,3,4,4) = -kappa3*w(1:Qdof)
    A_s(1:Qdof,3,4,5) = kappa1
    
    A_s(1:Qdof,3,5,1) = w*(kappa1*u2v2w2(1:Qdof) - kappaE(1:Qdof)/rho(1:Qdof))
    A_s(1:Qdof,3,5,2) = -kappa1*uw(1:Qdof)
    A_s(1:Qdof,3,5,3) = -kappa1*vw(1:Qdof)
    A_s(1:Qdof,3,5,4) = kappaE(1:Qdof)/rho - kappa1*(w2(1:Qdof) + u2v2w2(1:Qdof)/2.0)
    A_s(1:Qdof,3,5,5) = kappa*w(1:Qdof)


    deallocate(u,v,w,rho,e,uv,uw,vw,kappaE,u2,v2,w2,u2v2w2)

  end subroutine Set_A_s_Euler3D




  !> compute matrices 
  !> \f$ P^{\pm}  =  \left(\frac{D({\bf f_1}({\bf w})n_1+{\bf f_2}({\bf w})n_2}{D{\bf w}}
  !>  \right)^{\pm}\f$ 
  !> for the Euler equation
  subroutine Set_Ppm_Euler3D(ndim, Qdof, w, n, xi, Ppm, one_over_area, elem_i)
    integer, intent(in) :: Qdof, ndim, elem_i
    real, dimension(1:Qdof, 1:ndim), intent(in):: w !state  w in #Qdof nodes
    real, intent(in) :: one_over_area
    real, dimension(1:Qdof,1:2,1:ndim,1:ndim), intent(inout) :: Ppm
                                               ! matrices Ppm in  -- " --
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: n   ! outer normal 
    real, dimension(1:nbDim),intent(in) ::  xi                    ! node on the edge?
    real, dimension(1:ndim,1:ndim) :: t, t1
    real, dimension(1:ndim) :: dp, dm
    real, dimension (1:nbDim) :: nv
    real :: rlen
    integer :: ie, i, j

    do ie=1,Qdof
       call ComputeEigenValsVec3D(w(ie, 1:ndim), n(ie, 1:nbDim),xi(1:nbDim), &
            t1, t, dp, elem_i)


       do i=1,ndim
          state%max_eigenvals = max(state%max_eigenvals, &
               abs(dp(i)) *one_over_area )
       enddo
       
!       write(31,*)'Eigenvalues ',ie,dp(1),  one_over_area, state%max_eigenvals
       
       !calculating the negative (dm) and the positive(dp)  parts of eigenvalues
       do i=1,ndim                                                       
          dm(i) = 0.                                                         
          if(dp(i).lt.0.)then                                           
             dm(i) = dp(i)                                                
             dp(i) = 0.                                                   
          endif
       enddo
       !  multiplication of matrices: pp=t*dp*t1, where dp is diagonal
       !        print*,'Eigenvals+=',dp(1:ndim)
       !        print*,'Eigenvals-=',dm(1:ndim)
       do i=1,ndim
          do j=1,ndim
             Ppm(ie,1,i,j) = sum(t(i, 1:ndim) * dp(1:ndim) * t1(1:ndim, j) )
             Ppm(ie,2,i,j) = sum(t(i, 1:ndim) * dm(1:ndim) * t1(1:ndim, j) )
          enddo
       enddo
    enddo
    
  end subroutine Set_Ppm_Euler3D

  !> compute a descomposition of 
  !> \f$ \frac{D{\bf f_1}({\bf w})}{D{\bf w}} n_1 +
  !> \frac{D{\bf f_2}({\bf w})}{D{\bf w}} n_2 \f$ on positive and negative parts
  subroutine ComputeEigenValsVec3D(w, n, xi, t1, t, dp, elem_i)
    real :: w(1:ndim), n(1:nbDim), xi(1:nbDim), &
         t(1:ndim,1:ndim), t1(1:ndim,1:ndim), dp(1:ndim)
    real :: Q(1:ndim,1:ndim), Q1(1:ndim,1:ndim),&
         TT(1:ndim,1:ndim), TT1(1:ndim,1:ndim)
    real :: r, u,v,z,p,c, h, c2, c22, nv(nbDim), rlen, kap, kap1, n12, vv, tmp
    integer :: i, j, k, elem_i, ii
    integer :: ifile = 11 
    kap = state%model%kappa
    kap1 = state%model%kappa1
    
    !  r-density, u,v-velocity, p-pressure, kappa-poisson constant
    !  c-the local speed of sound, h-enthalpy          
    
    rlen = sqrt(sum(n(1:nbDim)*n(1:nbDim)))
    nv(1:nbDim) = n(1:nbDim)/rlen
    
    n12 = sqrt(nv(1)*nv(1) + nv(2)*nv(2))
    r=w(1)
    if (n12 .ne. 0.) then
       u = (w(2)*nv(1) + w(3)*nv(2) + w(4)*nv(3))/r
       v = (-w(2)*nv(2) + w(3)*nv(1))/(r*n12)
       z = (-w(2)*nv(1)*nv(3) - w(3)*nv(2)*nv(3) + w(4)*n12*n12)/(r*n12)
    else
       u = sign(1.,nv(3))*w(4)/r
       v = w(3)/r
       z = -sign(1.,nv(3))*w(2)/r
    endif
    
    vv=(u*u + v*v + z*z)/2.0
    
    
    p=kap1*(w(5)-r*vv)
    !     p=kap1*(w(5)-0.5*(sum(w(2:4)*w(2:4))/r))
    ! print*,'Ahoj',abs(p-kap1*(w(5)-r*vv)),n(:),rlen
    ! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
    if (abs(p-kap1*(w(5)-0.5*(sum(w(2:4)*w(2:4))/r))) > 0.0001) then
       print*,'problem with p in ComputeEigenValsVec3D'
    end if
    !print*,'EigenvalsVec=',r,u,v,p, kap*p/r 
    if( kap*p/r .le. 0. ) then
       open(ifile, file ='bad_ele', status = 'UNKNOWN')	   
       
       do ii = 1, grid%elem(elem_i)%flen	  
          write(ifile,*)  grid%x(grid%elem(elem_i)%face(idx,ii),1:nbDim)
       enddo
       write(ifile,*) grid%x(grid%elem(elem_i)%face(idx,1),1:nbDim) 	
       close(ifile)
       
       print*,'nonpositive square of speed of sound'
       print*,'pressure =', p
       print*,'w () = ', w(1:ndim)
       print*,'x =',xi(1:nbDim),elem_i
       
       !if(r .le. 0) then
       !   r= 0.001
       !endif
       !if(p .le. 0) then
       !   p= 0.001
       !endif
       stop
    endif
    
    c=sqrt(kap*p/r)
    h=c*c/kap1+(u*u+v*v)/2
    c2=c*c
    c22=2.*c2
    !  t and t1 - transformation matrices
    Q(:,:)=0.0
    Q1(:,:) = 0.0
    if (n12 .ne. 0.0) then
       Q(1,1) = 1.0
       Q(2,2) = nv(1)
       Q(2,3) = nv(2)
       Q(2,4) = nv(3)
       Q(3,2) = -nv(2)/n12
       Q(3,3) = nv(1)/n12
       Q(4,2) = -nv(1)*nv(3)/n12
       Q(4,3) = -nv(2)*nv(3)/n12
       Q(4,4) = n12
       Q(5,5) = 1.0
    else
       Q(1,1) = 1.0
       Q(2,2) = 0.0
       Q(2,3) = 0.0
       Q(2,4) = sign(1.,nv(3))
       Q(3,2) = 0.0
       Q(3,3) = 1.0
       Q(4,2) = -sign(1.,nv(3))
       Q(4,3) = 0.0
       Q(4,4) = 0.0
       Q(5,5) = 1.0
    endif
    
    Q1(1,1) = 1.0
    Q1(2,2) = 1.0
    Q1(3,3) = 1.0
    Q1(4,4) = 1.0
    Q1(5,5) = 1.0
    
    do i = 2,4
       if (Q(i,i) == 0.0) then
          j=i+1
          do while ((Q(j,i) == 0.0) .and. (j<5))
             j=j+1
          enddo
          if (j == 5) then
             ! mmmmmmmmmmmmm
             print*,'Problem in computation of eigenval Vec 3D '
          else
             do k=2,4
                !           swap(Q(i,k),Q(j,k))
                tmp=Q(i,k)
                Q(i,k)=Q(j,k)
                Q(j,k)=tmp
                !           swap(Q1(i,k),Q1(j,k))
                tmp=Q1(i,k)
                Q1(i,k)=Q1(j,k)
                Q1(j,k)=tmp
             enddo
          endif
       endif

       do j=4,2,-1
          Q1(i,j) =Q1(i,j)/Q(i,i)
       enddo

       do j=4,i,-1
          Q(i,j) =Q(i,j)/Q(i,i)
       enddo

       do j=2,4
          if (j .ne. i) then
             do k=4,2,-1
                Q1(j,k) = Q1(j,k)- Q1(i,k)*Q(j,i)
             enddo
             do k=4,i,-1
                Q(j,k) = Q(j,k) - Q(i,k)*Q(j,i)
             enddo
          endif
       enddo
    enddo

    if (n12 .ne. 0.) then
       Q(1,1) = 1.0
       Q(2,2) = nv(1)
       Q(2,3) = nv(2)
       Q(2,4) = nv(3)
       Q(3,2) = -nv(2)/n12
       Q(3,3) = nv(1)/n12
       Q(4,2) = -nv(1)*nv(3)/n12
       Q(4,3) = -nv(2)*nv(3)/n12
       Q(4,4) = n12
       Q(5,5) = 1.0
    else
       Q(1,1) = 1.0
       Q(2,2) = 0.0
       Q(2,3) = 0.0
       Q(2,4) = sign(1.,nv(3))
       Q(3,2) = 0.0
       Q(3,3) = 1.0
       Q(4,2) = -sign(1.,nv(3))
       Q(4,3) = 0.0
       Q(4,4) = 0.0
       Q(5,5) = 1.0
    endif
    
    TT(1,1) = 1.
    TT(1,2) = 1.
    TT(1,3) = 1.
    TT(1,4) = 1.
    TT(1,5) = 1.
    TT(2,1) = u - c
    TT(2,2) = u
    TT(2,3) = u
    TT(2,4) = u
    TT(2,5) = u + c
    TT(3,1) = v
    TT(3,2) = v
    TT(3,3) = v - c
    TT(3,4) = v
    TT(3,5) = v
    TT(4,1) = z
    TT(4,2) = z
    TT(4,3) = z
    TT(4,4) = z - c
    TT(4,5) = z
    TT(5,1) = vv + c2/kap1 - u*c
    TT(5,2) = vv
    TT(5,3) = vv - v*c
    TT(5,4) = vv - z*c
    TT(5,5) = vv + c2/kap1 + u*c
    
    TT1(1,1) = (kap1*vv + u*c)/2.0
    TT1(1,2) = -(c + kap1*u)/2.0
    TT1(1,3) = -v*kap1/2.0
    TT1(1,4) = -z*kap1/2.0
    TT1(1,5) = kap1/2.0
    TT1(2,1) = c2 - c*(v + z) - kap1*vv
    TT1(2,2) = u*kap1
    TT1(2,3) = c + v*kap1
    TT1(2,4) = c + z*kap1
    TT1(2,5) = -kap1
    TT1(3,1) = v*c
    TT1(3,2) = 0.0
    TT1(3,3) = -c
    TT1(3,4) = 0.0
    TT1(3,5) = 0.0
    TT1(4,1) = z*c
    TT1(4,2) = 0.0
    TT1(4,3) = 0.0
    TT1(4,4) = -c
    TT1(4,5) = 0.0
    TT1(5,1) = (kap1*vv - u*c)/2.0
    TT1(5,2) = (c - kap1*u)/2.0
    TT1(5,3) = -v*kap1/2.0
    TT1(5,4) = -z*kap1/2.0
    TT1(5,5) = kap1/2.0
    
    TT1(:,:) = TT1(:,:)/c2
    

    t(:,:)=0.0
    t1(:,:)=0.0

    do i=1,5
       do j=1,5
          t(i,j) = sum(Q1(i,1:ndim)*TT(1:ndim,j))
          t1(i,j) = sum(TT1(i,1:ndim)*Q(1:ndim,j))
       enddo
    enddo
    ! dp ... eigenvalues of the matrix n1*a(w)+n2*b(w)

    dp(2) = u
    dp(3) = dp(2)
    dp(4) = dp(2)
    dp(1) = dp(2)-c
    dp(5) = dp(2)+c

    dp(:) = dp(:)*rlen

    ! print*,'**************************************************************************'
    ! print*,TT(1,:)
    ! print*,TT1(:,1),sum(TT1(:,1)*dp(:))
    ! print*,w(:)
    ! print*,nv(:)
    ! print*,u,v,z,p,r,vv,c
    ! print*,'T'
    !   do i=1,5
    !     print*,t(i,:)
    !   enddo
    ! print*,'T1'
    !   do i=1,5
    !     print*,t1(i,:)
    !   enddo
    ! print*
    ! 
    ! print*,'T'
    !   do i=1,5
    !     print*,t(i,:)
    !   enddo
    ! print*,'T1'
    !   do i=1,5
    !     print*,t1(i,:)
    !   enddo
    ! print*
    
  end subroutine ComputeEigenValsVec3D

  !> compute 3D viscous fluxes R_s, s=1,3,3 for N.-S. equations
  !> in integ nodes
  subroutine Set_R_s_NS3D(Qdof, w, Dw, Re_1, R_s, xi)
    integer, intent(in) :: Qdof
    real, dimension(1:Qdof, 1:ndim), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndim, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    !real, intent(in) :: Re_1                     ! inverse of Reynolds number
    real, dimension(1:Qdof, 1:nbDim, 1:ndim), intent(inout) :: R_s
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
    real, dimension(:), allocatable :: u, v, z, oRe, e

    allocate( u(1:Qdof), v(1:Qdof), z(1:Qdof), e(1:Qdof), oRe(1:Qdof) )

    u(1:Qdof) = w(1:Qdof, 2)   / w(1:Qdof, 1)
    v(1:Qdof) = w(1:Qdof, 3)   / w(1:Qdof, 1)
    z(1:Qdof) = w(1:Qdof, 4)   / w(1:Qdof, 1)
    e(1:Qdof) = w(1:Qdof, 5)   / w(1:Qdof, 1)
    oRe(1:Qdof) = Re_1(1:Qdof) / w(1:Qdof, 1)

    !p(1:Qdof) = state%model%kappa1 * (w(1:Qdof,4)  &
    !     - (u(1:Qdof)*w(1:Qdof, 2) + v(1:Qdof)*w(1:Qdof, 3) )/2. )

    print*,'!!!!!', 'NOT UPDATED for 3D !!!! in model3D.f90 !!!'
    stop


    R_s(1:Qdof, 1:nbDim, 1) = 0.

    R_s(1:Qdof, 1, 2) = 2./3* oRe(1:Qdof) &
         *(2*(Dw(1:Qdof, 2, 1) - u(1:Qdof)*Dw(1:Qdof, 1, 1) )&
         - (Dw(1:Qdof, 3, 2) - v(1:Qdof)*Dw(1:Qdof, 1, 2) ) )

    R_s(1:Qdof, 1, 3) = oRe(1:Qdof) &
         *((Dw(1:Qdof, 3, 1) - v(1:Qdof)*Dw(1:Qdof, 1, 1) )&
         + (Dw(1:Qdof, 2, 2) - u(1:Qdof)*Dw(1:Qdof, 1, 2) ) )


    R_s(1:Qdof, 1, 4) = u(1:Qdof) * R_s(1:Qdof, 1, 2) + v(1:Qdof) * R_s(1:Qdof, 1, 3) &
         + state%model%kappa/state%model%Pr *  oRe(1:Qdof) &
         * ( Dw(1:Qdof, 4, 1) - e(1:Qdof)*Dw(1:Qdof, 1, 1) &
         - (u(1:Qdof) * Dw(1:Qdof, 2, 1) + v(1:Qdof) * Dw(1:Qdof, 3, 1)) &
         + (u(1:Qdof) * u(1:Qdof) + v(1:Qdof) *v(1:Qdof) ) *  Dw(1:Qdof, 1, 1) )



    R_s(1:Qdof, 2, 2) = oRe(1:Qdof) &
         *((Dw(1:Qdof, 3, 1) - v(1:Qdof)*Dw(1:Qdof, 1, 1) )&
         + (Dw(1:Qdof, 2, 2) - u(1:Qdof)*Dw(1:Qdof, 1, 2) ) )

    R_s(1:Qdof, 2, 3) = 2./3* oRe(1:Qdof) &
         *(2*(Dw(1:Qdof, 3, 2) - v(1:Qdof)*Dw(1:Qdof, 1, 2) )&
         - (Dw(1:Qdof, 2, 1) - u(1:Qdof)*Dw(1:Qdof, 1, 1) ) )


    R_s(1:Qdof, 2, 4) = u(1:Qdof) * R_s(1:Qdof, 2, 2) + v(1:Qdof) * R_s(1:Qdof, 2, 3) &
         + state%model%kappa/state%model%Pr *  oRe(1:Qdof) &
         * ( Dw(1:Qdof, 4, 2) - e(1:Qdof)*Dw(1:Qdof, 1, 2) &
         - (u(1:Qdof) * Dw(1:Qdof, 2, 2) + v(1:Qdof) * Dw(1:Qdof, 3, 2)) &
         + (u(1:Qdof) * u(1:Qdof) + v(1:Qdof) *v(1:Qdof) ) *  Dw(1:Qdof, 1, 2) )


    deallocate(u, v, e, oRe)

  end subroutine Set_R_s_NS3D
  


  !> compute matrices 5x5 K_sk, s,k=1,2,3 for N.-S. equations
  !> in integ nodes
  subroutine Set_K_sk_NS3D(Qdof, w, Re_1, K_sk, xi)
    integer, intent(in) :: Qdof
    real, dimension(1:Qdof, 1:ndim), intent(in):: w !state  w in #Qdof nodes
    real, intent(in) :: Re_1                     ! inverse of Reynolds number
    real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndim,1:ndim), intent(inout) :: K_sk
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates

    real :: ro, v1,v2,v3, hlpK51,hlpK52,hlpK53,hlpK54,hlpK55, one_over_Rew1, gamPr
    integer :: i


    !initialization
    K_sk(1:Qdof,1:nbDim,1:nbDim,1:ndim,1:ndim) = 0.

    gamPr = state%model%kappa/state%model%Pr

    do i=1,Qdof
       ro = w(i,1)
       v1 = w(i,2)/w(i,1)
       v2 = w(i,3)/w(i,1)
       v3 = w(i,4)/w(i,1)

       one_over_Rew1 = Re_1/w(i,1)

       hlpK51 = -(v1*v1+v2*v2+v3*v3)*one_over_Rew1 &
                  + one_over_Rew1*gamPr*(-w(i,5)/ro + v1*v1 + v2*v2 + v3*v3)
       hlpK52 = (1.-gamPr)*v1*one_over_Rew1
       hlpK53 = (1.-gamPr)*v2*one_over_Rew1
       hlpK54 = (1.-gamPr)*v3*one_over_Rew1
       hlpK55 =  gamPr*one_over_Rew1

       K_sk(i,1,1,2,1) = -4.0/3.0*v1*one_over_Rew1
       K_sk(i,1,1,2,2) = 4.0/3.0*one_over_Rew1
       K_sk(i,1,1,3,1) = -v2*one_over_Rew1
       K_sk(i,1,1,3,3) = one_over_Rew1
       K_sk(i,1,1,4,1) = -v3*one_over_Rew1
       K_sk(i,1,1,4,4) =  one_over_Rew1
       K_sk(i,1,1,5,1) = hlpK51-1.0/3.0*v1*v1*one_over_Rew1
       K_sk(i,1,1,5,2) = hlpK52 + 1.0/3.0*v1*one_over_Rew1
       K_sk(i,1,1,5,3) = hlpK53
       K_sk(i,1,1,5,4) = hlpK54
       K_sk(i,1,1,5,5) = hlpK55

       K_sk(i,2,2,2,1) = -v1*one_over_Rew1
       K_sk(i,2,2,2,2) =  one_over_Rew1
       K_sk(i,2,2,3,1) = -4.0/3.0*v2*one_over_Rew1
       K_sk(i,2,2,3,3) = 4.0/3.0*one_over_Rew1
       K_sk(i,2,2,4,1) = -v3*one_over_Rew1
       K_sk(i,2,2,4,4) = one_over_Rew1
       K_sk(i,2,2,5,1) = hlpK51-1.0/3.0*v2*v2*one_over_Rew1
       K_sk(i,2,2,5,2) = hlpK52
       K_sk(i,2,2,5,3) = hlpK53 + 1.0/3.0*v2*one_over_Rew1
       K_sk(i,2,2,5,4) = hlpK54
       K_sk(i,2,2,5,5) = hlpK55

       K_sk(i,3,3,2,1) = -v1*one_over_Rew1
       K_sk(i,3,3,2,2) =  one_over_Rew1
       K_sk(i,3,3,3,1) = -v2*one_over_Rew1
       K_sk(i,3,3,3,3) =  one_over_Rew1
       K_sk(i,3,3,4,1) = -4.0/3.0*v3*one_over_Rew1
       K_sk(i,3,3,4,4) =  4.0/3.0*one_over_Rew1
       K_sk(i,3,3,5,1) = hlpK51-1.0/3.0*v3*v3*one_over_Rew1
       K_sk(i,3,3,5,2) = hlpK52
       K_sk(i,3,3,5,3) = hlpK53
       K_sk(i,3,3,5,4) = hlpK54 + 1.0/3.0*v3*one_over_Rew1
       K_sk(i,3,3,5,5) = hlpK55

       K_sk(i,1,2,2,1) =  2.0/3.0*v2*one_over_Rew1
       K_sk(i,1,2,2,3) = -2.0/3.0*one_over_Rew1
       K_sk(i,1,2,3,1) = -v1*one_over_Rew1
       K_sk(i,1,2,3,2) =  one_over_Rew1
       K_sk(i,1,2,5,1) = -1.0/3.0*v1*v2*one_over_Rew1
       K_sk(i,1,2,5,2) =  v2*one_over_Rew1
       K_sk(i,1,2,5,3) = -2.0/3.0*v1*one_over_Rew1

       K_sk(i,1,3,2,1) =  2.0/3.0*v3*one_over_Rew1
       K_sk(i,1,3,2,3) = -2.0/3.0*one_over_Rew1
       K_sk(i,1,3,4,1) = -v1*one_over_Rew1
       K_sk(i,1,3,4,2) =  one_over_Rew1
       K_sk(i,1,3,5,1) = -1.0/3.0*v1*v3*one_over_Rew1
       K_sk(i,1,3,5,2) =  v3*one_over_Rew1
       K_sk(i,1,3,5,4) = -2.0/3.0*v1*one_over_Rew1

       K_sk(i,2,1,2,1) = -v2*one_over_Rew1
       K_sk(i,2,1,2,3) =  one_over_Rew1
       K_sk(i,2,1,3,1) =  2.0/3.0*v1*one_over_Rew1
       K_sk(i,2,1,3,2) = -2.0/3.0*one_over_Rew1
       K_sk(i,2,1,5,1) = -1.0/3.0*v1*v2*one_over_Rew1
       K_sk(i,2,1,5,2) = -2.0/3.0*v2*one_over_Rew1
       K_sk(i,2,1,5,3) =  v1*one_over_Rew1

       K_sk(i,2,3,3,1) =  2.0/3.0*v3*one_over_Rew1
       K_sk(i,2,3,3,3) = -2.0/3.0*one_over_Rew1
       K_sk(i,2,3,4,1) = -v2*one_over_Rew1
       K_sk(i,2,3,4,3) =  one_over_Rew1
       K_sk(i,2,3,5,1) = -1.0/3.0*v2*v3*one_over_Rew1
       K_sk(i,2,3,5,3) =  v3*one_over_Rew1
       K_sk(i,2,3,5,4) = -2.0/3.0*v2*one_over_Rew1

       K_sk(i,3,1,2,1) = -v3*one_over_Rew1
       K_sk(i,3,1,2,4) =  one_over_Rew1
       K_sk(i,3,1,4,1) =  2.0/3.0*v1*one_over_Rew1
       K_sk(i,3,1,4,2) = -2.0/3.0*one_over_Rew1
       K_sk(i,3,1,5,1) = -1.0/3.0*v1*v3*one_over_Rew1
       K_sk(i,3,1,5,2) = -2.0/3.0*v3*one_over_Rew1
       K_sk(i,3,1,5,4) =  v1*one_over_Rew1

       K_sk(i,3,2,3,1) = -v3*one_over_Rew1
       K_sk(i,3,2,3,4) =  one_over_Rew1
       K_sk(i,3,2,4,1) =  2.0/3.0*v2*one_over_Rew1
       K_sk(i,3,2,4,3) = -2.0/3.0*one_over_Rew1
       K_sk(i,3,2,5,1) = -1.0/3.0*v2*v3*one_over_Rew1
       K_sk(i,3,2,5,3) = -2.0/3.0*v3*one_over_Rew1
       K_sk(i,3,2,5,4) =  v2*one_over_Rew1

    enddo
  end subroutine Set_K_sk_NS3D

end module model3DNS
