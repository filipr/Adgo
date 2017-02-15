!> definition of 2D models which are simulated: Euler, Navier-Stokes
module model2DNS

  use data_mod
  use main_data
  use f_mapping
  use mesh_oper
  use blocks_integ
  use matrix_oper_int

  implicit none

  public:: Exact_Static_Shock_Vortex
  public:: Exact_Ringleb
  public:: Exact_Sod
  public:: Exact_Isent_Vortex
  public:: Exact_Gaussian_pulse
  public:: Exact_SteadyChannel
  public:: DMR_IC

  public:: SetEulerFluxes1

  public:: Set_f_s_Euler
  public:: Set_f_s_WS
  public:: Set_A_s_Euler
  public:: Set_A_s_WS
  public:: Set_Ppm_Euler
  public:: Set_Ppm_WS
  public:: ComputeEigenValsVec
  public:: ComputeEigenValsVec_WS
  public:: Set_R_s_NS
  public:: Set_R_s_WS
  public:: Set_K_sk_NS
  public:: Set_K_sk_WS

  public:: RotationForward
  public:: RotationBackward
  public:: SetBCCharacteristic
  public:: SetBCexactRiemann
  public:: ExactRiemannSolver
  public:: ReprepareBCCharacteristic
  public:: TimeDependentIOBC

contains

  !> mirror BC, \f$ {\bf w} \to \tilde{\bf w} \f$,
  !> \f$  {\bf w}\cdot{\bf n} = - \tilde{\bf w} \cdot{\bf n} \f$,
  !> \f$  {\bf w}\cdot{\bf t} =  \tilde{\bf w} \cdot{\bf t} \f$,
  !> nc is the UNIT normal !!!
  subroutine Mirror_W(ndim, Qdof, wi, nc)
    integer, intent(in) :: ndim, Qdof
    real, dimension(1:Qdof,1:ndim), intent(inout) :: wi          !  w in integ nodes
    real, dimension(1:Qdof,1:nbDim), intent(in) :: nc    ! normals in integ nodes
    real :: vn
    integer :: j

    if(ndim /= 4 .and. ndim /= 3) then
       print*,'Mirror_W only for dim = 3 or 4'
       stop
    endif

    ! v = v - 2*(v . un) un = v - 2*(v . n) n /|n|^2
    do j=1,Qdof
       vn = dot_product(wi(j, 2:3), nc(j, 1:nbDim)) !!/dot_product(nc(j,1:nbDim), nc(j,1:nbDim))
       wi(j,2:3) = wi(j, 2:3) - 2* vn * nc(j, 1:nbDim)
    enddo

  end subroutine Mirror_W

  !> mirror BC, \f$ {\bf w} \to \tilde{\bf w} \f$,
  !> \f$  {\bf w}\cdot{\bf n} = 0 \f$,
  !> \f$  {\bf w}\cdot{\bf t} =  \tilde{\bf w} \cdot{\bf t} \f$
  !> nc is the UNIT normal !!!
  subroutine UpdateMirror(ndim, Qdof, wi, nc)
    integer, intent(in) :: ndim, Qdof
    real, dimension(1:Qdof,1:ndim), intent(inout) :: wi          !  w in integ nodes
    real, dimension(1:Qdof,1:nbDim), intent(in) :: nc    ! normals in integ nodes
    real :: vn
    integer :: j

    if(ndim /= 4 .and. ndim /= 3) then
       print*,'UpdateMirror only for dim = 3 or 4'
       stop
    endif

    ! v = v - (v . un) un = v - (v . n) n /|n|^2
    do j=1,Qdof
       vn = dot_product(wi(j, 2:3), nc(j, 1:nbDim)) /dot_product(nc(j,1:nbDim), nc(j,1:nbDim))
       wi(j,2:3) = wi(j, 2:3) -  vn * nc(j, 1:nbDim)
    enddo

  end subroutine UpdateMirror

 !> Evaluate exact solution \f$wi\f$ at node \f$x\f$ for
 !> the steady state solution in a channel with a source term
 subroutine Exact_SteadyChannel(x, wi)
   real, dimension(1:nbDim), intent(in) :: x
   real, dimension(1:ndim), intent(out) :: wi
   real:: rho, v1, v2, pressure, energy, temperature
   real :: p0, p1

   rho = 1.
   p0 = 2.
   p1 = 0.01
   v1 = state%model%Re*p1/2. * x(2)*( 1. - x(2) )
   v2 = 0.

   pressure = p0  - p1 * x(1)
   temperature =  pressure / state%model%kappa1 / rho
   energy = rho * temperature + rho * (v1 * v1 + v2 * v2)/2.

   wi(1) = rho
   wi(2) = rho * v1
   wi(3) = rho * v2
   wi(4) = energy

 end subroutine Exact_SteadyChannel

 !> evaluate RHS  \f$ f \f$ at node \f$x\f$
 !> for  the steady state solution in a channel with a source term
 subroutine RHS_SteadyChannel(x, f, t)
   real, dimension(1:nbDim), intent(in) :: x
   real, dimension(1:ndim), intent(out) :: f
   real, intent(in) :: t

   real:: rho, v1, v2, pressure, energy, temperature, dedx
   real :: p0, p1

   rho = 1.
   p0 = 2.
   p1 = 0.01
   v1 = state%model%Re*p1/2. * x(2)*( 1. - x(2) )
   v2 = 0.

   pressure = p0  - p1 * x(1)
   temperature =  pressure / state%model%kappa1 / rho
   energy = rho * temperature + rho * (v1 * v1 + v2 * v2)/2.

   dedx = - p1/state%model%kappa1 + rho *(state%model%Re * p1/2)**2 *x(2)*(1-x(2))*(1-2*x(2))

   f(1:3) = 0.
   !f(4) = v1*(dedx - p1) - state%model%Re * (p1/2)**2*((1-2*x(2))**2 - 2*x(2)*(1-x(2) ))
   f(4) = -0.25/state%model%kappa1 * p1*p1*state%model%Re &
        *(state%model%kappa*(1-2*x(2))**2 - 1 + 6*x(2)*(1- x(2) ))

 end subroutine RHS_SteadyChannel


 !> Evaluate exact solution \f$wi\f$ at node \f$x\f$ for Ringleb flow problem
 subroutine Exact_Gaussian_pulse(x, wi)
   real, dimension(1:nbDim), intent(in) :: x
   real, dimension(1:ndim), intent(out) :: wi

   real:: epsilon = 0.1
   real :: r, val

   r = dot_product(x(1:nbDim), x(1:nbDim) )

   val = epsilon * exp(-log(2.)/(log(exp(1.) ))* r)

   wi(1) = state%rho_infty + val
   wi(2) = wi(1) * state%v_infty * cos(state%alpha_infty)
   wi(3) = wi(1) * state%v_infty * sin(state%alpha_infty)
   wi(4) = (state%p_infty +val)/state%model%kappa1 &
        + 0.5*(wi(2)*wi(2) + wi(3)*wi(3))/wi(1)

   !print*, state%rho_infty, state%v_infty, state%p_infty
   !write(88,*) x(1:nbDim), wi(1:ndim)
 end subroutine Exact_Gaussian_pulse


 !> static shock placed at X0 = 0
 subroutine Exact_Static_Shock_Vortex(x, wi, ctime)
   real, dimension(1:nbDim), intent(in) :: x
   real, dimension(1:ndim), intent(out) :: wi
   real, intent(in) :: ctime
   real :: x0, kappa, r1, r2, u1, u2, p1,p2, coef, coef2, M_in
   real, dimension(1:2) :: xr, xi
   real :: rr, uc, rc, uu, beta, epsilon, bar_v1, bar_v2, f, theta, eta
   real:: p_l, rho_l, v1_l, v2_l
   real :: kP, kM, kK, C1, C2
   real :: rr1, rr2, Mn, vm
   real :: ve, p_s, r_s, u_s, xx_s, x_s
   integer :: i

   ! smoothing parameter
   !ve = 5E-3  TOO small
   !ve = 0.99E-4
   !ve = max (1E-4, 10.**(-state%time%recompute_back -1) )
   !ve = 1E-3
   !ve = 2E-4   ! ATAN smoothing
   !ve = 5E-4   ! CUBIC smoothing

   ve = state%model%Re1 !/ 2  ! CUBIC smoothing
   !if(ve <= 1E-4) ve = 2E-4

   !ve = 1E-4
   !ve = 1E-5 TOO much

   xr(1) = 0.5
   xr(2) = 1.0

   !uc = 0.25
   uc = 0.5
   rc = 0.075
   C1 = uc / rc
   C2 = 1. / (2 * rc * rc)

   x0 = 1.
   M_in = 1.1588


   kappa = state%model%kappa

   r1 = 1.
   u1 = kappa**0.5 * M_in
   p1 = 1.

   !furst
   kP  = (kappa + 1)/ 2.
   kM  = (kappa - 1)/ 2.
   kK  = (kappa - 1) /kappa

   coef = kP * M_in**2 /(1 + kM * M_in**2)
   coef2 = (kappa * M_in**2 - kM)/ kP

   r2 = r1 * coef
   u2 = u1 / coef
   p2 = p1 * coef2

   if(x(1) < x0) then
      r_s = r1
      u_s = u1
      p_s = p1
   else
      r_s = r2
      u_s = u2
      p_s = p2

   endif

   if(state%space%adapt%adapt_level == 0) then
      ! ATAN smoothing of the shock
      !xx_s = (x(1) - x0)/ve
      !r_s = (r1 + r2) / 2 + (r2 - r1)/pi*atan(xx_s)
      !u_s = (u1 + u2) / 2 + (u2 - u1)/pi*atan(xx_s)
      !p_s = (p1 + p2) / 2 + (p2 - r1)/pi*atan(xx_s)

      ! cubic smoothing of he shock
      if( abs(x(1) - x0) < ve) then
         x_s = x(1) - x0
         xx_s = x_s *(x_s**2 - 3 * ve**2) / (4 * ve**3)
         r_s = (r1 + r2) /2 + (r1- r2 ) *  xx_s
         u_s = (u1 + u2) /2 + (u1- u2 ) *  xx_s
         p_s = (p1 + p2) /2 + (p1- p2 ) *  xx_s
      endif


   endif

   !if(x(1) < x0) then
   !   write(54,*) x(1), r1, u1, p1, r_s, u_s, p_s
   !else
   !   write(54,*) x(1), r2, u2, p2, r_s, u_s, p_s
   !endif


   !C1 = 0.   ! no vertex
   if(x(1) < x0) then
      ! adding of vortex
      xi(1:2) = x(1:2) - xr(1:2)
      rr = dot_product(xi(1:2), xi(1:2) )**0.5

      uu = C1 * rr * exp(- C2 * rr*rr)

      p_l =  p_s**kK - kK*(C1**2 /(4*C2)) *exp(-2*C2*rr**2)
      p_l = p_l**(1./kK)

      rho_l = r_s * (p_l /p_s)**(1./kappa)

      v1_l = u_s - uu * xi(2)/rr
      v2_l = 0. + uu * xi(1)/rr


      ! furst begin
      !rr1 = 0.075
      !rr2 = 0.175
      !vm = 0.5
      !Mn = vm/(kappa * p1/ rho1)

      ! furst end

      wi(1) = rho_l
      wi(2) = rho_l * v1_l
      wi(3) = rho_l * v2_l
      wi(4) = p_l /(kappa -1) + 0.5 * rho_l * (v1_l**2 + v2_l**2)

      !wi(1) = r1
      !wi(2) = r1 * u1
      !wi(3) = r1 * 0.
      !wi(4) = p1 /(kappa -1) + 0.5 * r1 * u1**2
   else
      ! NO stationary SHOCK
      !r2 = r1
      !u2 = u1
      !p2 = p1

      wi(1) = r_s
      wi(2) = r_s * u_s
      wi(3) = 0.
      wi(4) = p_s /(kappa -1) + 0.5 * r_s * u_s**2
   endif


 end subroutine Exact_Static_Shock_Vortex

 ! evaluate the angle bewteen vector x(1:2) and x_1 axis
 function Eval_Angle(x)
   real:: Eval_Angle
   real, dimension(1:2), intent(in) :: x
   real :: Pi2

   Pi2 = asin(1.)

   if( dot_product(x,x) <= 0.) then
      print*,'Zero angle'
      Eval_Angle = 0.
   else
      if(x(1) == 0) then
         if(x(2) > 0) then
            Eval_Angle = Pi2
         else
            Eval_Angle = 3* Pi2
         end if
      elseif(x(2) == 0) then
         if(x(1) > 0) then
            Eval_Angle = 0.
         else
            Eval_Angle = 2* Pi2
         end if
      else
         Eval_Angle = atan( x(2) / x(1) )

         if(x(1) >  0. .and. x(2) > 0) then
            Eval_Angle = Eval_Angle + 0.

         elseif(x(1) <  0. .and. x(2) > 0) then
            Eval_Angle = Eval_Angle + 2 * Pi2

         elseif(x(1) <  0. .and. x(2) < 0) then
            Eval_Angle = Eval_Angle + 2 *Pi2

         elseif(x(1) >  0. .and. x(2) < 0) then
            Eval_Angle = Eval_Angle + 4* Pi2
         endif
      endif

   endif

 end function Eval_Angle

 !> Evaluate exact solution \f$wi\f$ at node \f$x\f$ for isentropic vortex
 subroutine Exact_Isent_Vortex(x, wi, ctime)
   real, dimension(1:nbDim), intent(in) :: x
   real, dimension(1:ndim), intent(out) :: wi
   real, intent(in) :: ctime
   real, dimension(1:nbDim) :: xi, xii
   real:: xx, radius2, r2, epsilon, bar_v1, bar_v2, bar_rho, bar_p
   real:: f, kappa, rho_l, v1_l, v2_l, p_l
   real:: eta, theta
   integer :: l


   epsilon = 5.
   kappa = state%model%kappa

   xx =  mod (5. + ctime , 10.) ! [xx,xx] position of the center of the vertex

   xi(1:2) = x(1:2)  ! coordinates

   ! shift of the period
   do l=1,nbDim ! = 2
      if(abs(xi(l) - xx) > 5.) then
         if(xi(l) > xx) then
            xi(l) = xi(l) - 10
         elseif(xi(l) < xx) then
            xi(l) = xi(l) + 10
         endif
      endif
   enddo

   !if(x(1) > 9.9 .and. x(2) > 9.9)  write(*,'(a6,12es12.4)') '???', &
   !     x(1:2), xx, xi(1:2)

   xi(1:2) = xi(1:2) - xx  ! relative coordinates with respect to the vortex center
   radius2 = dot_product(xi(1:2), xi(1:2) )

   !if(x(1) > 9.9 .and. x(2) > 9.9)  write(*,'(a6,12es12.4)') '?..', &
   !     x(1:2), xx, xi(1:2), radius2**0.5

   bar_v1 = epsilon/2/pi*exp(0.5*(1-radius2))*(-xi(2))
   bar_v2 = epsilon/2/pi*exp(0.5*(1-radius2))*xi(1)

   f = -(kappa-1)*epsilon*epsilon/8/kappa/pi/pi*exp(1-radius2)

   v1_l  = bar_v1  + state%BC(1)%w(2)
   v2_l  = bar_v2  + state%BC(1)%w(3)

   theta = state%BC(1)%w(4) / state%BC(1)%w(1) + f
   eta   = state%BC(1)%w(4) / (state%BC(1)%w(1)**kappa)

   rho_l = (theta / eta)**(1./(kappa - 1.) )
   p_l   = rho_l * theta

   wi(1) = rho_l
   wi(2) = v1_l*rho_l
   wi(3) = v2_l*rho_l
   wi(4) = p_l/(kappa-1)+0.5*rho_l*(v1_l*v1_l+v2_l*v2_l)

   !write(82,'(20e12.4)') xx, r,v,p, state%BC(1)%ww(1:ndim)
   !write(*,'(a6,12e12.4)') 'VORTEX:',state%BC(1)%w(1:ndim), theta, eta, wi(1:4)

   !write(*,'(a6,12e12.4)') 'VORTEX:', theta, eta, wi(1:4)


   !print*,';-------------------'

 end subroutine Exact_Isent_Vortex


 !> Evaluate exact solution \f$wi\f$ at node \f$x\f$ for Ringleb flow problem
 subroutine Exact_Ringleb(x, wi)
   real, dimension(1:nbDim), intent(in) :: x
   real, dimension(1:ndim), intent(out) :: wi

   real:: cspeed

   cspeed = cspeed_ringleb(state%model%kappa, x(1), x(2) )
   call state_ringleb(state%model%kappa, x(1), cspeed, wi(1:ndim))

   !write(*,'(a3,7e12.4)') 'Ri:',x(1:nbDim), cspeed, wi(1:ndim)
 end subroutine Exact_Ringleb


  function cspeed_ringleb(gamma,  x, y)
    real :: cspeed_ringleb
    real, intent(in) :: gamma, x, y
    real :: tol, clow, chigh, cnew, err
    integer ::  nit

    tol = 1.e-15
    clow  = .5000000
    chigh = .9999999

    nit = 0
    err = 1.0
    do nit = 1,100
       cnew = (clow+chigh)/2.
       if( fun_ringleb(gamma,x,y,cnew)*fun_ringleb(gamma,x,y,chigh) > 0.) then
          chigh = cnew
       else
          clow  = cnew
       endif

       err  = abs(chigh-clow)
       if (err < tol) goto 100
    enddo
100 continue

    if(err > tol) then
       print*,"Error tolerance was not achieved in bisection iteration: ",err,'>',tol
       stop
    endif
    cspeed_ringleb = cnew
  end function cspeed_ringleb


  function fun_ringleb(gamma, x,  y, c)
    real :: fun_ringleb
    real, intent(in) :: gamma, x,  y, c
    real:: jval, gb, rho, q2

    gb = gamma - 1
    rho = c**(2./gb)
    jval = 1./c + 1./(3.*c**3) + 1./(5.*c**5) -  0.5*log((1+c)/(1.-c) )
    q2 = 2.*(1.-c*c)/gb

    fun_ringleb = (x-jval/2.)**2 + y*y - 1./(4.*rho*rho*q2*q2)
  end function fun_ringleb


  subroutine state_ringleb(gamma, x, c, soln)
    real, intent(in) :: gamma, x, c
    real, dimension(1:ndim), intent(out) :: soln
    real gb, jval, rho, qval, kval, u,v,p,en,q2

    gb = gamma - 1.

    rho = c**(2./gb)
    jval = 1./c + 1./(3.*c**3) + 1./(5.*c**5) -  0.5*log((1+c)/(1.-c) )
    q2   = 2.*(1.-c*c)/gb
    qval = sqrt(q2)
    kval = sqrt( 2./(1./q2 -2.*rho*(x - jval/2.)))

    v  = q2/kval
    u  = sqrt( abs(q2 - v*v) )
    p  = c*c*rho/gamma
    en = p/gb + 0.5*rho*(u*u + v*v)

    soln(1) = rho
    soln(2) = u*rho
    soln(3) = v*rho
    soln(4) = en
  end subroutine state_ringleb


 !> Evaluate exact solution \f$wi\f$ at node \f$x\f$ for Sod tube
 subroutine Exact_Sod(x, wi, ctime)
   real, dimension(1:nbDim), intent(in) :: x
   real, dimension(1:ndim), intent(out) :: wi
   real, intent(in) :: ctime
   real :: ti, rk, pl, pr, rl, rr, cl, cr, p2, rleft, vc, vs, x1, xx, xx0, p, v, r
   xx0 = 0.
   if(ctime == 0.) then
      if(x(1) < xx0) then
         wi(1:ndim) = state%BC(1)%ww(1:ndim)
      else
         wi(1:ndim) = state%BC(2)%ww(1:ndim)
      endif
   else
      ti = ctime

      rk = state%model%kappa
      pl = 1.
      pr = 0.1

      rl = 1.0
      rr = 0.125

      cl = (rk*pl/rl)**0.5
      cr = (rk*pr/rr)**0.5

      p2 = 0.303130209

      rleft = (2./rk/(rk-1))**0.5*(p2/pr -1)/(1+(rk+1)/(rk-1)*p2/pr)**0.5

      vc = cr*rleft
      vs = cr**2*(p2/pr-1)/rk/vc

      x1 = xx0 +  ((rk +1)/2*vc-cl)*ti
      xx = x(1)

      if(xx < (xx0-cl*ti)) then
         p = pl
         v = 0.
         r = rl

      elseif((xx >=  xx0-cl*ti) .and. (  xx < x1)) then
         p = pl*(1-(rk-1)/(rk+1)*(1+(xx-xx0)/cl/ti))**(2*rk/(rk-1))
         v = 2/(rk+1)*((xx-xx0)/ti + cl)
         r =rl*(1-(rk-1)/(rk+1)*(1+(xx-xx0)/cl/ti))**(2/(rk-1))

      elseif((xx >= x1).and.( xx < xx0 + vc*ti)) then
         p = p2
         v = vc
         r = rl*(p2/pl)**(1./rk)

      elseif((xx >= xx0+vc*ti).and.(xx < xx0 + vs*ti)) then
         p = p2
         v = vc
         r = rr*(1+(rk+1)/(rk-1)*p2/pr)/((rk+1)/(rk-1) + p2/pr)

      elseif(xx >= xx0+vs*ti)then
         p = pr
         v = 0.
         r = rr
      else
         print *,'problems in  Exact_Sod'

      endif

      wi(1) = r
      wi(2) = r * v
      wi(3) = 0.
      wi(4) = (rk - 1) * p + 0.5*r * v *v

   endif
 end subroutine Exact_Sod


 !> Evaluate IC \f$wi\f$ at node \f$x\f$ for double Mach reflection
 subroutine DMR_IC(x, wi, t)
   real, dimension(1:nbDim), intent(in) :: x
   real, dimension(1:ndim), intent(out) :: wi
   real, intent(in) :: t

   real, dimension(1:ndim) :: qL, qR
   real :: steep
   integer :: i

   !print*,'###',state%space%h

!   print*,(3.+3.*0.012)/state%space%h, 6/state%space%h, state%space%h
   steep = min(4., 2. + 2.*t/0.001 )/state%space%h
!   steep = min(20., 2. + 30.*t/0.001 )/state%space%h
!   steep = 2./state%space%h
   !steep = 1./state%space%h

   do i=1,state%numBC
      if(state%BC(i)%inout .eq. 1) then
         qL(1:ndim) = state%BC(i)%ww(1:ndim)
      elseif(state%BC(i)%inout .eq. 0) then
         qR(1:ndim) = state%BC(i)%ww(1:ndim)
      else
         stop 'Trouble in model2DNS.f90, DMR_IC'
      endif
   enddo

   wi(1:ndim) = qL(1:ndim) + (qR(1:ndim) - qL(1:ndim)) &
              * (tanh(steep*(1./6. + (x(2)+20.*t) / 3**0.5 - x(1)) ) + 1)/2

 end subroutine DMR_IC




  subroutine Prepare_A_s_Euler(max_Qdof, A_s)
    ! matrices A_s for Euler equations in integration nodes
    ! setting of constant elements
    integer, intent(in) :: max_Qdof
    real, dimension(1:max_Qdof,1:nbDim,1:ndim,1:ndim), intent(inout) :: A_s
    integer :: i

    A_s(1:max_Qdof,1:nbDim,1:ndim,1:ndim) = 0.

    do i=1,max_Qdof
       A_s(i, 1, 1, 2) = 1.
       A_s(i, 1, 2, 4) = state%model%kappa - 1.
       A_s(i, 2, 1, 3) = 1.
       A_s(i, 2, 3, 4) = state%model%kappa - 1.
    enddo

  end subroutine Prepare_A_s_Euler



  !> compute the mach number in integ nodes from wi
  subroutine EvalMachNumber(Qdeg, w, mach)
    integer, intent(in) :: Qdeg
    real, dimension(1:Qdeg, 1:ndim), intent(in):: w !state  w
    real, dimension(1:Qdeg), intent(out):: mach
    integer :: l
    real :: r,u,v,p

    do l=1,Qdeg
       r = w(l,1)
       u = w(l,2)/w(l,1)
       v = w(l,3)/w(l,1)
       p = state%model%kappa1*(w(l,4) - r*(u*u + v*v)/2)

       mach = ((u*u + v*v)/(state%model%kappa * p /r))**0.5
    enddo
  end subroutine EvalMachNumber

  !> compute the Euler fluxes \f$ {\bf f}_1({\bf w})\f$ and
  !> \f$ {\bf f}_2({\bf w})\f$  w in integ nodes
  !> IT IS NEEDED??
  subroutine SetEulerFluxes1(Qdeg, w, f1, f2)
    integer, intent(in) :: Qdeg
    real, dimension(1:Qdeg, 1:ndim), intent(in):: w !state  w
    real, dimension(1:Qdeg, 1:ndim), intent(out):: f1, f2 !state  w
    integer :: l
    real :: r,u,v,p

    do l=1,Qdeg
       r = w(l,1)
       u = w(l,2)/w(l,1)
       v = w(l,3)/w(l,1)
       p = state%model%kappa1*(w(l,4) - r*(u*u + v*v)/2)

       f1(l,1) = w(l,2)
       f1(l,2) = r*u*u + p
       f1(l,3) = r*u*v
       f1(l,4) = (w(l,4) + p)*u

       f2(l,1) = w(l,3)
       f2(l,2) = r*u*v
       f2(l,3) = r*v*v + p
       f2(l,4) = (w(l,4) + p)*v
    enddo

  end subroutine SetEulerFluxes1

  !> compute the Euler fluxes \f$ f_s({\bf w}),\quad s=1,2\f$
  subroutine Set_f_s_Euler(ndimL, nbDim, Qdof, w, f_s, xi, ie)
    integer, intent(in) :: Qdof, nbDim, ndimL
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
                                               ! fluxes f_s in  -- " --
    real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
    integer, intent(in) :: ie
    real, dimension(:), allocatable :: u,v,p

    allocate( u(1:Qdof), v(1:Qdof),  p(1:Qdof) )

    u(1:Qdof) = w(1:Qdof, 2) / w(1:Qdof, 1)
    v(1:Qdof) = w(1:Qdof, 3) / w(1:Qdof, 1)
    p(1:Qdof) = state%model%kappa1 * (w(1:Qdof,4)  &
         - (u(1:Qdof)*w(1:Qdof, 2) + v(1:Qdof)*w(1:Qdof, 3) )/2. )


    f_s(1:Qdof, 1, 1) = w(1:Qdof, 2)
    f_s(1:Qdof, 1, 2) = w(1:Qdof, 2)*u(1:Qdof) + p(1:Qdof)
    f_s(1:Qdof, 1, 3) = w(1:Qdof, 2)*v(1:Qdof)
    f_s(1:Qdof, 1, 4) = (w(1:Qdof,4) + p(1:Qdof) ) * u(1:Qdof)

    f_s(1:Qdof, 2, 1) = w(1:Qdof,3)
    f_s(1:Qdof, 2, 2) = w(1:Qdof, 3)*u(1:Qdof)
    f_s(1:Qdof, 2, 3) = w(1:Qdof, 3)*v(1:Qdof) + p(1:Qdof)
    f_s(1:Qdof, 2, 4) = (w(1:Qdof,4) + p(1:Qdof) ) * v(1:Qdof)


    deallocate(u,v,p )

  end subroutine Set_f_s_Euler

  !> compute the Euler fluxes for wet steam case \f$ f_s({\bf w}),\quad s=1,2\f$
   subroutine Set_f_s_WS(ndimL, nbDim,Qdof, w, f_s, xi, ie)
     integer, intent(in) :: Qdof, nbDim, ndimL    ! Qdof: number of integration nodes, ndimL: number of quantities for metric evaluation
     real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
     real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
                                                ! fluxes f_s in  -- " --
     real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi   ! integ nodes
     integer, intent(in) :: ie
     real, dimension(:), allocatable :: u,v,p

     allocate( u(1:Qdof), v(1:Qdof) )

     u(1:Qdof) = w(1:Qdof, 2) / w(1:Qdof, 1)
     v(1:Qdof) = w(1:Qdof, 3) / w(1:Qdof, 1)

     f_s(1:Qdof, :, :) = 0.
     call Set_f_s_Euler(4, nbDim, Qdof, w, f_s(1:Qdof, 1:2, 1:4), xi, ie)

     f_s(1:Qdof, 1, 5) = w(1:Qdof,5)*u(1:Qdof)
     f_s(1:Qdof, 1, 6) = w(1:Qdof,6)*u(1:Qdof)
     f_s(1:Qdof, 1, 7) = w(1:Qdof,7)*u(1:Qdof)
     f_s(1:Qdof, 1, 8) = w(1:Qdof,8)*u(1:Qdof)

     f_s(1:Qdof, 2, 5) = w(1:Qdof,5)*v(1:Qdof)
     f_s(1:Qdof, 2, 6) = w(1:Qdof,6)*v(1:Qdof)
     f_s(1:Qdof, 2, 7) = w(1:Qdof,7)*v(1:Qdof)
     f_s(1:Qdof, 2, 8) = w(1:Qdof,8)*v(1:Qdof)

     !write(*,*) '***************'
     !write(*,*) 'w(1,:) ',  w(1,:) ! wet steam test
     !write(*,*) '---------------------'
     !write(*,*) 'f_s(1,2,:) ', f_s(1,2,:)  !
     !write(*,*) '***************'
     !call WriteMatrixA(eta)

     deallocate(u,v)

   end subroutine Set_f_s_WS

  !> compute matrices \f$ A_s = \frac{D{\bf f_s}({\bf w})}{D{\bf w}},\quad s=1,2\f$
  !> for Euler fluxes
  subroutine Set_A_s_Euler(ndimL, nbDim, Qdof, w, A_s, xi, ie)
    integer, intent(in) :: ndimL, nbDim, Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: A_s
                                               ! matrices A_s in  -- " --
    real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
    integer, intent(in) :: ie
    real :: kappa, k1, k3
    real, dimension(:), allocatable :: u,v,rho,e,uv,ke,u2,v2
    integer :: i


    allocate( u(1:Qdof), v(1:Qdof), rho(1:Qdof), e(1:Qdof) )
    allocate( uv(1:Qdof), ke(1:Qdof), u2(1:Qdof), v2(1:Qdof) )

    kappa = state%model%kappa
    k3=state%model%kappa-3
    k1=state%model%kappa1

    rho(1:Qdof) = w(1:Qdof,1)
    u(1:Qdof) = w(1:Qdof,2)/rho(1:Qdof)
    v(1:Qdof) = w(1:Qdof,3)/rho(1:Qdof)
    e(1:Qdof) = w(1:Qdof,4)
    uv(1:Qdof)= u(1:Qdof)*v(1:Qdof)
    u2(1:Qdof)= u(1:Qdof)*u(1:Qdof)
    v2(1:Qdof)=v(1:Qdof)*v(1:Qdof)
    ke(1:Qdof) = kappa*e(1:Qdof)

    A_s(1:Qdof, 1, 1, 1) = 0.

    A_s(1:Qdof, 1, 1, 2) = 1.
    A_s(1:Qdof, 1, 1, 3:4) = 0.

    A_s(1:Qdof,1,2,1) = k3*u2(1:Qdof)/2+k1*v2(1:Qdof)/2;
    A_s(1:Qdof,1,2,2) = -k3*u(1:Qdof);
    A_s(1:Qdof,1,2,3) = -k1*v(1:Qdof);
    A_s(1:Qdof, 1, 2, 4) = kappa - 1.

    A_s(1:Qdof,1,3,1) = -uv(1:Qdof);
    A_s(1:Qdof,1,3,2) = v(1:Qdof);
    A_s(1:Qdof,1,3,3) = u(1:Qdof);
    A_s(1:Qdof, 1, 3, 4) = 0.

    A_s(1:Qdof,1,4,1) = -ke(1:Qdof)*u(1:Qdof)/rho(1:Qdof) &
         +k1*u(1:Qdof)*(u2(1:Qdof)+v2(1:Qdof));
    A_s(1:Qdof,1,4,2) = ke(1:Qdof)/rho(1:Qdof)-k1*(3*u2(1:Qdof)+v2(1:Qdof))/2;
    A_s(1:Qdof,1,4,3) = -k1*uv(1:Qdof);
    A_s(1:Qdof,1,4,4) = kappa*u(1:Qdof);

    A_s(1:Qdof, 2, 1, 1:2) = 0.
    A_s(1:Qdof, 2, 1, 3) = 1.
    A_s(1:Qdof, 2, 1, 4) = 0.

    A_s(1:Qdof,2,2,1) = -uv(1:Qdof);
    A_s(1:Qdof,2,2,2) = v(1:Qdof);
    A_s(1:Qdof,2,2,3) = u(1:Qdof);
    A_s(1:Qdof, 2, 2, 4) = 0.

    A_s(1:Qdof,2,3,1) = k3*v2(1:Qdof)/2+k1*u2(1:Qdof)/2;
    A_s(1:Qdof,2,3,2) = -k1*u(1:Qdof);
    A_s(1:Qdof,2,3,3) = -k3*v(1:Qdof);
    A_s(1:Qdof, 2, 3, 4) = kappa - 1.

    A_s(1:Qdof,2,4,1) = -ke(1:Qdof)*v(1:Qdof)/rho(1:Qdof) &
         +k1*v(1:Qdof)*(u2(1:Qdof)+v2(1:Qdof));
    A_s(1:Qdof,2,4,2) = -k1*uv(1:Qdof);
    A_s(1:Qdof,2,4,3) = ke(1:Qdof)/rho(1:Qdof)-k1*(u2(1:Qdof)+3*v2(1:Qdof))/2;
    A_s(1:Qdof,2,4,4) = kappa*v(1:Qdof);


    deallocate(u,v,rho,e,uv,ke,u2,v2 )
  end subroutine Set_A_s_Euler

  !> compute matrices \f$ A_s = \frac{D{\bf f_s}({\bf w})}{D{\bf w}},\quad s=1,2\f$
  !> for wet steam
  subroutine Set_A_s_WS(ndimL, nbDim, Qdof, w, A_s, xi, ie)
    integer, intent(in) :: Qdof, nbDim, ndimL
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: A_s
                                               ! matrices A_s in  -- " --
    real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
    integer, intent(in) :: ie

    real :: kappa, k1, k3
    real, dimension(:), allocatable :: u,v,rho,omega,Q2,Q1,Q0
    integer :: i

    allocate(u(1:Qdof),v(1:Qdof),rho(1:Qdof),omega(1:Qdof),Q2(1:Qdof),Q1(1:Qdof),Q0(1:Qdof))

    rho(1:Qdof) = w(1:Qdof,1)
    u(1:Qdof) = w(1:Qdof,2)/rho(1:Qdof)
    v(1:Qdof) = w(1:Qdof,3)/rho(1:Qdof)
    omega(1:Qdof) = w(1:Qdof,5)/rho(1:Qdof)
    Q2(1:Qdof)= w(1:Qdof,6)/rho(1:Qdof)
    Q1(1:Qdof)= w(1:Qdof,7)/rho(1:Qdof)
    Q0(1:Qdof)=w(1:Qdof,8)/rho(1:Qdof)

    A_s(:,:,:,:) = 0.


    call Set_A_s_Euler(4, nbDim, Qdof, w, A_s(1:Qdof,1:nbDim,1:4,1:4), xi, ie)

    A_s(1:Qdof, 1, 1, 5) = - omega(1:Qdof)*u(1:Qdof)
    A_s(1:Qdof, 1, 1, 6) = - Q2(1:Qdof)*u(1:Qdof)
    A_s(1:Qdof, 1, 1, 7) = - Q1(1:Qdof)*u(1:Qdof)
    A_s(1:Qdof, 1, 1, 8) = - Q0(1:Qdof)*u(1:Qdof)

    A_s(1:Qdof, 1, 2, 5) =  omega(1:Qdof)
    A_s(1:Qdof, 1, 2, 6) =  Q2(1:Qdof)
    A_s(1:Qdof, 1, 2, 7) =  Q1(1:Qdof)
    A_s(1:Qdof, 1, 2, 8) =  Q0(1:Qdof)

    A_s(1:Qdof, 1, 3, 5:8) = 0.

    A_s(1:Qdof, 1, 4, 5:8) = 0.

    A_s(1:Qdof, 1, 5, 1:4) = 0.
    A_s(1:Qdof, 1, 5, 5) = u(1:Qdof)
    A_s(1:Qdof, 1, 5, 6:8) = 0.

    A_s(1:Qdof, 1, 6, 1:5) = 0.
    A_s(1:Qdof, 1, 6, 6) = u(1:Qdof)
    A_s(1:Qdof, 1, 6, 7:8) = 0.

    A_s(1:Qdof, 1, 7, 1:6) = 0.
    A_s(1:Qdof, 1, 7, 7) = u(1:Qdof)
    A_s(1:Qdof, 1, 7, 8) = 0.

    A_s(1:Qdof, 1, 8, 1:7) = 0.
    A_s(1:Qdof, 1, 7, 8) = u(1:Qdof)


    A_s(1:Qdof, 2, 1, 5) = - omega(1:Qdof)*v(1:Qdof)
    A_s(1:Qdof, 2, 1, 6) = - Q2(1:Qdof)*v(1:Qdof)
    A_s(1:Qdof, 2, 1, 7) = - Q1(1:Qdof)*v(1:Qdof)
    A_s(1:Qdof, 2, 1, 8) = - Q0(1:Qdof)*v(1:Qdof)

    A_s(1:Qdof, 2, 2, 5:8) = 0.

    A_s(1:Qdof, 2, 3, 5) =  omega(1:Qdof)
    A_s(1:Qdof, 2, 3, 6) =  Q2(1:Qdof)
    A_s(1:Qdof, 2, 3, 7) =  Q1(1:Qdof)
    A_s(1:Qdof, 2, 3, 8) =  Q0(1:Qdof)

    A_s(1:Qdof, 2, 4, 5:8) = 0.

    A_s(1:Qdof, 2, 5, 1:4) = 0.
    A_s(1:Qdof, 2, 5, 5) = v(1:Qdof)
    A_s(1:Qdof, 2, 5, 6:8) = 0.

    A_s(1:Qdof, 2, 6, 1:5) = 0.
    A_s(1:Qdof, 2, 6, 6) = v(1:Qdof)
    A_s(1:Qdof, 2, 6, 7:8) = 0.

    A_s(1:Qdof, 2, 7, 1:6) = 0.
    A_s(1:Qdof, 2, 7, 7) = v(1:Qdof)
    A_s(1:Qdof, 2, 7, 8) = 0.


    !write(*,*) 'A_s_Euler:'  ! wet steam writing
    !do i=1,ndim
    !  write(*,*) A_s(1,1,i,:)
    !  write(*,*) '-----------------'
    !end do

    deallocate(u,v,rho,omega,Q2,Q1,Q0)

  end subroutine Set_A_s_WS


  !> compute matrices 4x4 K_sk, s,k=1,2 for N.-S. equations
  !> in integ nodes
  subroutine Set_K_sk_NS(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, K_sk, xi)
    integer, intent(in) :: ndimL, nbDim, iRe, Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes, not USED
    real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    !real, intent(in) :: Re_1                     ! inverse of Reynolds number
    real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: K_sk
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
    real :: v1, v2, E, one_over_Rew1, gamPr
    integer :: i

    !initialization
    K_sk(1:Qdof,1:nbDim,1:nbDim,1:ndimL,1:ndimL) = 0.

    do i=1,Qdof
       one_over_Rew1 = Re_1(1,1)/w(i,1)
       v1 = w(i,2)/w(i,1)
       v2 = w(i,3)/w(i,1)
       E = w(i,4)/w(i,1)
       gamPr = state%model%kappa/state%model%Pr

       !s=1	: k=1
       K_sk(i,1,1,2,1) = -4./3*v1
       K_sk(i,1,1,2,2) = 4./3

       K_sk(i,1,1,3,1) = -v2
       K_sk(i,1,1,3,3) = 1.

       K_sk(i,1,1,4,1) = -(4./3*v1*v1 + v2*v2 + gamPr*(E - v1*v1 - v2*v2) )

       K_sk(i,1,1,4,2) = (4./3 - gamPr) * v1
       K_sk(i,1,1,4,3) = (1. - gamPr) *v2
       K_sk(i,1,1,4,4) = gamPr

       !s=1 : k=2
       K_sk(i,1,2,2,1) = 2./3*v2
       K_sk(i,1,2,2,3) = -2./3

       K_sk(i,1,2,3,1) = -v1
       K_sk(i,1,2,3,2) = 1.

       K_sk(i,1,2,4,1) = -1./3*v1*v2
       K_sk(i,1,2,4,2) = v2
       K_sk(i,1,2,4,3) = -2./3*v1

       !s=2 : k=1
       K_sk(i,2,1,2,1) = -v2
       K_sk(i,2,1,2,3) = 1.

       K_sk(i,2,1,3,1) = 2./3*v1
       K_sk(i,2,1,3,2) = -2./3

       K_sk(i,2,1,4,1) = -1./3*v1*v2
       K_sk(i,2,1,4,2) = -2./3*v2
       K_sk(i,2,1,4,3) = v1

      !s=2	: k=2
       K_sk(i,2,2,2,1) = -v1
       K_sk(i,2,2,2,2) = 1.

       K_sk(i,2,2,3,1) = -4./3*v2
       K_sk(i,2,2,3,3) = 4./3

       K_sk(i,2,2,4,1) = -(v1*v1 + 4./3*v2*v2 + gamPr*(E - v1*v1 - v2*v2))
       K_sk(i,2,2,4,2) = (1. - gamPr)*v1
       K_sk(i,2,2,4,3) = (4./3 - gamPr)*v2
       K_sk(i,2,2,4,4) = gamPr

       K_sk(i,1:nbDim,1:nbDim,1:ndimL,1:ndimL) &
            = K_sk(i,1:nbDim,1:nbDim,1:ndimL,1:ndimL) * one_over_Rew1
    enddo
  end subroutine Set_K_sk_NS


  !> compute matrices 8x8 K_sk, s,k=1,2 for wet steam equations
  !> in integ nodes
  subroutine Set_K_sk_WS(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, K_sk, xi)
    integer, intent(in) :: ndimL, nbDim, iRe, Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes, not USED
    real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    !real, intent(in) :: Re_1                     ! inverse of Reynolds number
    real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: K_sk
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
    integer :: i,j,k,l

    !> wet steam values of K_sk are equal to 0
    K_sk(1:Qdof,1:nbDim,1:nbDim,1:ndimL,1:ndimL) = 0.

    !> set K_sk for NS equations
    call Set_K_sk_NS(4, nbDim, iRe, Qdof, w(1:Qdof, 1:4), Dw(1:Qdof,1:4,1:nbDim), Re_1, &
         K_sk(1:Qdof,1:nbDim,1:nbDim,1:4,1:4), xi)

    !! wet steam writing for testing
    !write(*,*) "------K_sk-------"
    !do i=1,ndimL
    !  !do j=1,ndimL
    !     write(*,*) K_sk(1,1,1,i,:)
    !     write(*,*)
    !  !end do
    !end do
    !write(*,*) "---------------"

  end subroutine Set_K_sk_WS


  !> compute viscous fluxes R_s, s=1,2 for N.-S. equations
  !> in integ nodes
  subroutine Set_R_s_NS(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
    integer, intent(in) :: Qdof, nbDim, iRe, ndimL
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    !real, intent(in) :: Re_1                     ! inverse of Reynolds number
    real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
    real, dimension(:), allocatable :: u, v, oRe, e

    allocate( u(1:Qdof), v(1:Qdof), e(1:Qdof), oRe(1:Qdof) )

    u(1:Qdof) = w(1:Qdof, 2)   / w(1:Qdof, 1)
    v(1:Qdof) = w(1:Qdof, 3)   / w(1:Qdof, 1)
    e(1:Qdof) = w(1:Qdof, 4)   / w(1:Qdof, 1)
    oRe(1:Qdof) = Re_1(1,1:Qdof) / w(1:Qdof, 1)

    !p(1:Qdof) = state%model%kappa1 * (w(1:Qdof,4)  &
    !     - (u(1:Qdof)*w(1:Qdof, 2) + v(1:Qdof)*w(1:Qdof, 3) )/2. )


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

    !write(*,*) '===== Dw(1,1,:) ======='
    !write(*,*) Dw(1,1,:)
    !write(*,*) '--------------'
    !write(*,*) R_s(1,2,:)
    !write(*,*) '================='

    deallocate(u, v, e, oRe)

  end subroutine Set_R_s_NS


  !> compute viscous fluxes R_s, s=1,2 for wet steam equations
  !> in integ nodes
  subroutine Set_R_s_WS(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
    integer, intent(in) :: Qdof, iRe, nbDim, ndimL
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates


    R_s(1:Qdof, 1:nbDim, 1:ndimL) = 0.
    ! set R_s for wet steam equations
    !R_s(1:Qdof,1,5:8) = 0.
    !R_s(1:Qdof,2,5:8) = 0.

    ! set R_s fof NS equations
    call Set_R_s_NS(4, nbDim, iRe, Qdof, w(1:Qdof,1:4), Dw(1:Qdof,1:4,1:nbDim), Re_1, &
         R_s(1:Qdof,1:nbDim,1:4), xi )


    ! wet steam writing for testing
    !write(*,*) '===== R_s ======='
    !write(*,*) R_s(1,1,:)
    !write(*,*) '--------------'
    !write(*,*) R_s(1,2,:)
    !write(*,*) '================='

  end subroutine Set_R_s_WS


  !> compute reactive terms S,  for empty equation in integ nodes
  subroutine Set_S_empty(ndimL, nbDim, Qdof, xi, w, Dw, S)
    integer, intent(in) :: ndimL, nbDim, Qdof
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL), intent(inout) :: S
    integer :: i,j,k

    S(1:Qdof, :) = 0.
    !do k=1,Qdof
    !   S(k, 1) = Eval_Reaction_Coeffs(w(k,1), Dw(k,1,1:nbDim), xi(k, 1:nbDim), 0 )
    !enddo

  end subroutine Set_S_empty

  !> compute derivative of the reactive terms S,  for empty equation in integ nodes
  subroutine Set_DS_empty(ndimL, nbDim, Qdof, xi, w, Dw, DS)
    integer, intent(in) :: ndimL, nbDim, Qdof
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:ndimL), intent(inout) :: DS
    integer :: i,j,k

    DS(1:Qdof, :,:) = 0.
    !do k=1,Qdof
    !   DS(k, 1) = Eval_Reaction_Coeffs(w(k,1), Dw(k,1,1:nbDim), xi(k, 1:nbDim), 1 )
    !enddo

  end subroutine Set_DS_empty


  !> compute a descomposition of
  !> \f$ \frac{D{\bf f_1}({\bf w})}{D{\bf w}} n_1 +
  !> \frac{D{\bf f_2}({\bf w})}{D{\bf w}} n_2 \f$ on positive and negative parts
  subroutine ComputeEigenValsVec(grid, ndimL, w, n, xi, t1, t, dp, elem)
    class(mesh), intent(in) :: grid
    integer, intent(in) :: ndimL
    real :: w(1:ndimL), n(1:nbDim), xi(1:nbDim), t(1:4,1:4), t1(1:4,1:4), dp(1:4)
    real :: r, u,v,p,c, h, c2, c22, nv(2), rlen, kap, kap1
    type(element), intent(in) :: elem
    integer :: ii
    integer :: ifile = 11
    kap = state%model%kappa
    kap1 = state%model%kappa1

    !  r-density, u,v-velocity, p-pressure, kappa-poisson constant
    !  c-the local speed of sound, h-enthalpy

    rlen = sqrt(n(1)*n(1) + n(2)*n(2) )
    nv(1:nbDim) = n(1:nbDim)/rlen


    r=w(1)
    u=w(2)/r
    v=w(3)/r
    p=kap1*(w(4)-0.5*r*(u*u+v*v))

    !print*,'EigenvalsVec=',r,u,v,p, kap*p/r

    if( kap*p/r .le. 0.) then
       open(ifile, file ='bad_ele', status = 'UNKNOWN', position='append')

       do ii = 1, elem%flen
          write(ifile,*)  grid%x(elem%face(idx,ii),1:nbDim)
       enddo
       write(ifile,*) grid%x(elem%face(idx,1),1:nbDim)
       write(ifile,'(x)')
       close(ifile)

       ! print*,'nonpositive square of speed of sound (1)'
       ! print*,'pressure =', p
       ! print*,'w () = ', w(1:ndimL)
       ! print*,'x =',xi(1:nbDim),elem%i
       ! print*,'xc=',elem%xc(1:nbDim)

       write(*,'(a20,2(a4,es10.2),a9,i5,a9,i2, a16)') &
            'nonpositive p/ro:', &
            ' p =', p, &
            ' r =', r, &
            ', elem%i=',elem%i,  &
            ', Newton= ', state%nlSolver%iter , &
            ', file "bad_ele"'
       if(r .le. 0) then
          r= 0.001
       endif
       if(p .le. 0) then
          p= 0.001
       endif
       !stop
    endif

    c=sqrt(kap*p/r)
    h=c*c/kap1+(u*u+v*v)/2.
    c2=c*c
    c22=2.*c2
    !  t and t1 - transformation matrices

    t(1,1) = 1.
    t(1,2) = 0.
    t(1,3) = 1./c22
    t(1,4) = 1./c22
    t(2,1) = u
    t(2,2) = nv(2)
    t(2,3) = (u+c*nv(1))/c22
    t(2,4) = (u-c*nv(1))/c22
    t(3,1) = v
    t(3,2) = -nv(1)
    t(3,3) = (v+c*nv(2))/c22
    t(3,4) = (v-c*nv(2))/c22
    t(4,1) = (u*u+v*v)/2.
    t(4,2) = nv(2)*u-nv(1)*v
    t(4,3) = (h+c*(nv(1)*u+nv(2)*v))/c22
    t(4,4) = (h-c*(nv(1)*u+nv(2)*v))/c22

    t1(1,1) = 1-kap1/c2*(u*u/2+v*v/2)
    t1(1,2) = kap1/c2*u
    t1(1,3) = kap1/c2*v
    t1(1,4) = -kap1/c2
    t1(2,1) = nv(1)*v-nv(2)*u
    t1(2,2) = nv(2)
    t1(2,3) = -nv(1)
    t1(2,4) = 0.
    t1(3,1) = -c*(nv(1)*u+nv(2)*v)+kap1*(u*u/2+v*v/2)
    t1(3,2) = c*nv(1)-kap1*u
    t1(3,3) = c*nv(2)-kap1*v
    t1(3,4) = kap1
    t1(4,1) = c*(nv(1)*u+nv(2)*v)+kap1*(u*u/2+v*v/2)
    t1(4,2) = -c*nv(1)-kap1*u
    t1(4,3) = -c*nv(2)-kap1*v
    t1(4,4) = kap1

    ! dp ... eigenvalues of the matrix n1*a(w)+n2*b(w)

    dp(1)=n(1)*u+n(2)*v
    dp(2)=dp(1)
    dp(3)=dp(1)+c*rlen
    dp(4)=dp(1)-c*rlen


  end subroutine ComputeEigenValsVec

  !> compute a descomposition of
  !> \f$ \frac{D{\bf f_1}({\bf w})}{D{\bf w}} n_1 +
  !> \frac{D{\bf f_2}({\bf w})}{D{\bf w}} n_2 \f$ on positive and negative parts
  !> for the wet steam equations
  subroutine ComputeEigenValsVec_WS(grid, ndimL, w, n, xi, t1, t, dp, elem)
    class(mesh), intent(in) :: grid
    integer, intent(in) :: ndimL
    real :: w(1:ndimL), n(1:nbDim), xi(1:nbDim), t(1:ndimL,1:ndimL), t1(1:ndimL,1:ndimL), dp(1:ndimL)
    real :: r, u,v,p,c, h, c2, c22, nv(2), rlen, kap, kap1, omega, L, p_ws  ! wet steam - omega = mass fraction of liquid water, L - latent heat of condensation, p_ws - wet steam pressure
    type(element), intent(in) :: elem
    integer :: ii, i,j
    integer :: ifile = 11
    kap = state%model%kappa
    kap1 = state%model%kappa1

    !  r-density, u,v-velocity, p-pressure, kappa-poisson constant
    !  c-the local speed of sound, h-enthalpy

    rlen = sqrt(n(1)*n(1) + n(2)*n(2) )
    nv(1:nbDim) = n(1:nbDim)/rlen


    r=w(1)
    u=w(2)/r
    v=w(3)/r

    p=kap1*(w(4)-0.5*r*(u*u+v*v))
    ! wet steam
    omega=w(5)/r
    L = 2260.E3  ! L = 2260 kJ/kg (Wikipedia) ! should be define in module wet_steam_paramets.f90
    p_ws = kap1*(1-omega)/(1+omega*kap1)*(w(4)-0.5*r*(u*u+v*v) + r*omega*L)
    !write(*,*) "w(i)", (w(i), i=1,ndimL)
    !write(*,*) "omega, p, p_ws", omega, p, p_ws

    !print*,'EigenvalsVec=',r,u,v,p, kap*p/r

    if( kap*p/r .le. 0. ) then
       open(ifile, file ='bad_ele', status = 'UNKNOWN')

       do ii = 1, elem%flen
          write(ifile,*)  grid%x(elem%face(idx,ii),1:nbDim)
       enddo
       write(ifile,*) grid%x(elem%face(idx,1),1:nbDim)
       write(ifile,'(x)')
       close(ifile)

       print*,'nonpositive square of speed of sound (2)'
       print*,'pressure =', p
       print*,'w () = ', w(1:ndimL)
       print*,'x =',xi(1:nbDim),elem%i
       print*,'xc=',elem%xc(1:nbDim)

       if(r .le. 0) then
          r= 0.001
       endif
       if(p .le. 0) then
          p= 0.001
       endif
       stop
    endif

    c=sqrt(kap*p/r)
    h=c*c/kap1+(u*u+v*v)/2.
    c2=c*c
    c22=2.*c2
    !  t and t1 - transformation matrices

    t(1,1) = 1.
    t(1,2) = 0.
    t(1,3) = 1./c22
    t(1,4) = 1./c22
    t(2,1) = u
    t(2,2) = nv(2)
    t(2,3) = (u+c*nv(1))/c22
    t(2,4) = (u-c*nv(1))/c22
    t(3,1) = v
    t(3,2) = -nv(1)
    t(3,3) = (v+c*nv(2))/c22
    t(3,4) = (v-c*nv(2))/c22
    t(4,1) = (u*u+v*v)/2.
    t(4,2) = nv(2)*u-nv(1)*v
    t(4,3) = (h+c*(nv(1)*u+nv(2)*v))/c22
    t(4,4) = (h-c*(nv(1)*u+nv(2)*v))/c22
    ! wet steam part
    t(1:4,5:ndimL) = 0.
    t(5,1:4) = 0.
    t(5,5) = 1.
    t(5,6:ndimL) = 0.
    t(6,1:5) = 0.
    t(6,6) = 1.
    t(6,7:ndimL) = 0.
    t(7,1:6) = 0.
    t(7,7) = 1.
    t(7,8) = 0.
    t(8,1:7) = 0.
    t(8,8) = 1.

    t1(1,1) = 1.-kap1/c2*(u*u/2.+v*v/2.)
    t1(1,2) = kap1/c2*u
    t1(1,3) = kap1/c2*v
    t1(1,4) = -kap1/c2
    t1(2,1) = nv(1)*v-nv(2)*u
    t1(2,2) = nv(2)
    t1(2,3) = -nv(1)
    t1(2,4) = 0.
    t1(3,1) = -c*(nv(1)*u+nv(2)*v)+kap1*(u*u/2+v*v/2)
    t1(3,2) = c*nv(1)-kap1*u
    t1(3,3) = c*nv(2)-kap1*v
    t1(3,4) = kap1
    t1(4,1) = c*(nv(1)*u+nv(2)*v)+kap1*(u*u/2+v*v/2)
    t1(4,2) = -c*nv(1)-kap1*u
    t1(4,3) = -c*nv(2)-kap1*v
    t1(4,4) = kap1

    ! wet steam part
    t1(1:5,4:ndimL) = 0.
    t1(5,1:4) = 0.
    t1(5,5) = 1.
    t1(5,6:ndimL) = 0.
    t1(6,1:5) = 0.
    t1(6,6) = 1.
    t1(6,7:ndimL) = 0.
    t1(7,1:6) = 0.
    t1(7,7) = 1.
    t1(7,8) = 0.
    t1(8,1:7) = 0.
    t1(8,8) = 1.

    ! dp ... eigenvalues of the matrix n1*a(w)+n2*b(w)

    dp(1)=n(1)*u+n(2)*v
    dp(2)=dp(1)
    dp(3)=dp(1)+c*rlen
    dp(4)=dp(1)-c*rlen
    ! wet steam part
    dp(5:ndimL)=u*nv(1)+v*nv(2)

    !! wet steam writing for testing
    !do i=1,ndimL
    !  write(*,*) '==== t, t1 ======'
    !  do j=1, ndimL
    !    write(*,*) t(i,j), t1(i,j)
    !  end do
    !  write(*,*)
    !  write(*,*) dp(i)
    !  write(*,*) '=========='
    !end do

  end subroutine ComputeEigenValsVec_WS


  !> compute matrices
  !> \f$ P^{\pm} = \left(\frac{D({\bf f_1}({\bf w})n_1+{\bf f_2}({\bf w})n_2}{D{\bf w}}
  !>  \right)^{\pm}\f$
  !> for the Euler equation
  subroutine Set_Ppm_Euler(ndimL, nbDim, Qdof, w, n, xi, Ppm, one_over_area, elem)
    !class(mesh), intent(in) :: grid
    integer, intent(in) :: Qdof, ndimL, nbDim
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: Ppm
                                               ! matrices Ppm in  -- " --
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: n   ! outer normal
    real, dimension(1:Qdof, 1:nbDim),intent(in) ::  xi                    ! node on the edge?
    real, intent(in), optional :: one_over_area
    type(element), intent(inout), optional :: elem
    real, dimension(1:ndimL,1:ndimL) :: t, t1
    real, dimension(1:ndimL) :: dp, dm
    real, dimension (1:nbDim) :: nv
    real :: rlen
    integer :: ie, i, j

    !print*, 'Set_Ppm_Euler', Qdof

    if (.not. present(elem) .or. .not. present(one_over_area) ) &
      stop 'elem and one_over_area must be present in Set_Ppm_Euler!'

    do ie=1,Qdof
       call ComputeEigenValsVec(grid,ndimL, w(ie, 1:ndimL), n(ie, 1:nbDim), &
            xi(ie,1:nbDim),&
            t1, t, dp, elem)

      ! print*, 't = ', t
      ! print*, 'dp = ', dp

       do i=1,ndimL
          elem%max_eigenvals = abs(dp(i)) *one_over_area
          state%max_eigenvals = max(state%max_eigenvals, elem%max_eigenvals  )
       enddo

       !calculating the negative (dm) and the positive(dp)  parts of eigenvalues
       do i=1,ndimL
          dm(i)=0.
          if(dp(i).lt.0.)then
             dm(i)=dp(i)
             dp(i)=0.
          endif
       enddo

       !  multiplication of matrices: pp=t*dp*t1, where dp is diagonal

       !print*,'Eigenvals+=',dp(1:ndimL)
       !print*,'Eigenvals-=',dm(1:ndimL)
       do i=1,ndimL
          do j=1,ndimL
             Ppm(ie,1,i,j) = sum(t(i, 1:ndimL) * dp(1:ndimL) * t1(1:ndimL, j) )
             Ppm(ie,2,i,j) = sum(t(i, 1:ndimL) * dm(1:ndimL) * t1(1:ndimL, j) )

          enddo
       enddo

    enddo

  end subroutine Set_Ppm_Euler

  !> compute matrices
  !> for the wet steam equations
  subroutine Set_Ppm_WS(ndimL, nbDim, Qdof, w, n, xi, Ppm, one_over_area, elem)
    !class(mesh), intent(in) :: grid
    integer, intent(in) :: Qdof, ndimL, nbDim
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: Ppm
                                               ! matrices Ppm in  -- " --
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: n   ! outer normal
    real, dimension(1:Qdof, 1:nbDim),intent(in) ::  xi                    ! node on the edge?
    real, intent(in), optional :: one_over_area
    type(element), intent(inout), optional :: elem
    real, dimension(1:ndimL,1:ndimL) :: t, t1
    real, dimension(1:ndimL) :: dp, dm
    real, dimension (1:nbDim) :: nv
    real :: rlen
    integer :: ie, i, j

    if (.not. present(elem) .or. .not. present(one_over_area) ) &
      stop 'elem and one_over_area must be present in Set_Ppm_Euler!'


    do ie=1,Qdof
       call ComputeEigenValsVec_WS(grid,ndimL, w(ie, 1:ndimL), n(ie, 1:nbDim), &
            xi(ie,1:nbDim),&
            t1, t, dp, elem)



       do i=1,ndimL
          elem%max_eigenvals = abs(dp(i)) *one_over_area
          state%max_eigenvals = max(state%max_eigenvals, elem%max_eigenvals  )
!          write(*,*) t(i,:), dp(i)
       enddo
!       ! Only for testing
!       do i=1,ndimL
!         do j=1,ndimL
!           write(*,*) i, j, t(i,j), t1(i,j), dp(i)
!         end do
!       end do
!       stop
!       ! Only for testing
!       write(31,*)'Eigenvalues ',ie,dp(1),  one_over_area, state%max_eigenvals

       !calculating the negative (dm) and the positive(dp)  parts of eigenvalues
       do i=1,ndimL
          dm(i)=0.
          if(dp(i).lt.0.)then
             dm(i)=dp(i)
             dp(i)=0.
          endif
       enddo

       !  multiplication of matrices: pp=t*dp*t1, where dp is diagonal

       !print*,'Eigenvals+=',dp(1:ndimL)
       !print*,'Eigenvals-=',dm(1:ndimL)
       do i=1,ndimL
          do j=1,ndimL
             Ppm(ie,1,i,j) = sum(t(i, 1:ndimL) * dp(1:ndimL) * t1(1:ndimL, j) )
             Ppm(ie,2,i,j) = sum(t(i, 1:ndimL) * dm(1:ndimL) * t1(1:ndimL, j) )

          enddo
       enddo

       !create P
       !call Matrix_plus_minus(ndimL, P, Ppm(ie, 1:2, 1:ndimL, 1:ndimL) )
    enddo

  end subroutine Set_Ppm_WS


  subroutine Set_Ppm_Euler_Slip( ndimL, nbDim, e_Qdof, w_ein, n_e, Ppm )
    ! compute matrix Pp on a slip boundary
    integer, intent(in) :: e_Qdof, ndimL, nbDim
    real, dimension(1:e_Qdof, 1:ndimL), intent(in):: w_ein !state  w in #Qdof nodes
    real, dimension(1:e_Qdof,2:3,1:ndimL), intent(inout) :: Ppm
                                               ! matrices Ppm in  -- " --
    real, dimension(1:e_Qdof,1:nbDim), intent(in) :: n_e      ! outer normal
    real :: kappa, kappa1
    real :: v(2), vv, nn(2)
    integer :: ie, i, j

    kappa = state%model%kappa
    kappa1 = state%model%kappa1

    do ie=1,e_Qdof
       v(1) = w_ein(ie,2)/w_ein(ie,1)
       v(2) = w_ein(ie,3)/w_ein(ie,1)
       vv = dot_product(v, v)
       nn(1:nbDim) = n_e(ie,1:nbDim)*kappa1
       !Ppm(ie,1,1:ndimL) = 0.    !  first line of Pp
       !Ppm(ie,4,1:ndimL) = 0.    !  last  line of Pp

       !print*,'SLIP ne:',n_e(ie,1:nbDim),nn(1:nbDim)

       do j=1,2
          Ppm(ie,j+1,1) = vv*nn(j)/2
          Ppm(ie,j+1,2) = -v(1)*nn(j)
          Ppm(ie,j+1,3) = -v(2)*nn(j)
          Ppm(ie,j+1,4) = nn(j)
       enddo
       !do j=1,4
       !   print*,'FW =',Ppm(ie,j,1:ndimL)
       !enddo
       !print*,'*****************',ie
    enddo
  end subroutine Set_Ppm_Euler_Slip


  subroutine RotationForward(ndimL, wi, n, qi)
    integer, intent(in) :: ndimL
    real, dimension(1:2) :: n
    real, dimension(1:ndimL) :: wi, qi

    qi(1) = wi(1);
    qi(2) = wi(2)*n(1) + wi(3)*n(2);
    qi(3) = -wi(2)*n(2) + wi(3)*n(1);
    qi(4) = wi(4);
  end subroutine RotationForward

  subroutine RotationBackward(ndimL, wi, n, qi)
    integer, intent(in) :: ndimL
    real, dimension(1:2) :: n
    real, dimension(1:ndimL) :: wi, qi

    qi(1) = wi(1);
    qi(2) = wi(2)*n(1) - wi(3)*n(2);
    qi(3) = wi(2)*n(2) + wi(3)*n(1);
    qi(4) = wi(4);
  end subroutine RotationBackward


  !> repreparation of BC, a use of the old (characteristic) approach
  subroutine ReprepareBCCharacteristic(Qdof, ndimL, wi, wD, n, press_extrap, xc)
    integer, intent(in) :: Qdof, ndimL
    real, dimension(1:Qdof,1:ndimL), intent(in) :: wi
    real, dimension(1:Qdof,1:ndimL), intent(inout) ::  wD
    real, dimension(1:nbDim), intent(in) :: n, xc
    real,  intent(in) ::  press_extrap
    integer :: i,j
    real :: size, nn(2), vn, p, c, pD

    if(ndimL /= 4) then
       print*,'Bad implementation of ReprepareCharacteristic'
       stop
    endif
    size = (dot_product(n, n))**0.5
    nn(1:nbDim) = n(1:nbDim)/size

    do i=1, Qdof
       vn = dot_product(nn(1:nbDim), wi(i,2:3))/wi(i,1)

       p = state%model%kappa1 *(wi(i,4)-dot_product(wi(i,2:3),wi(i,2:3))/wi(i,1)/2 )
       c = (state%model%kappa * p /wi(i,1))**0.5D+00

       !write(*,'(a2,7es10.3)' ) '>>',wi(i,1:ndim), p
       if(vn .lt. 0) then
          if(vn .lt. -c) then
             ! supersonic inlet
             ! no action
          else
             ! subsonic inlet
             wD(i,4) = p/state%model%kappa1 + dot_product(wD(i,2:3),wD(i,2:3))/wD(i,1)/2
          endif
       else
          if(vn .gt. c) then
             ! supersonic outlet
             wD(i,1:ndimL) = wi(i,1:ndimL)
          else
             ! subsonic outlet
             pD = state%model%kappa1*(wD(i,4)-dot_product(wD(i,2:3),wD(i,2:3))/wD(i,1)/2)

             ! NEW correction in te sence of the mean value of pressure
             if(press_extrap > 0.)  then
                !write(22,'(8es12.4)' ) &
                !     xc(1:nbDim),pD, p, press_extrap, pD + p - press_extrap,vn

                pD = pD + p - press_extrap
             endif

             wD(i,1:3) = wi(i,1:3)
             wD(i,4) = pD/state%model%kappa1 + dot_product(wi(i,2:3),wi(i,2:3))/wi(i,1)/2
          endif
       endif
    enddo

  end subroutine ReprepareBCCharacteristic



  !> rsetting of BC, a use of the exact Riemann problem
  subroutine SetBCexactRiemann(Qdof, ndimL, wi, wD, w_BC, n, xi )
    integer, intent(in) :: Qdof, ndimL
    real, dimension(1:Qdof,1:ndimL), intent(in) :: wi, wD
    real, dimension(1:Qdof,1:ndimL), intent(out) ::  w_BC
    real, dimension(1:nbDim), intent(in) :: n, xi
    integer :: i

    if(ndimL /= 4) then
       print*,'Bad implementation of SetBCExactRiemann'
       stop
    endif

    do i=1, Qdof
       call ExactRiemannSolver(ndimL,wi(i,1:ndimL), wD(i,1:ndimL), n(1:nbDim),&
            w_BC(i,1:ndimL))
    enddo

  end subroutine SetBCexactRiemann



  !>     Calculates the 2D Euler flux using the exact Riemann Solver
  !>
  !>     NB:  Input parameters are the left (UL) and right (UR) state
  !>     vectors of conserved variables,  with
  !>     U = (rho, rho u, rho v, rho w, rho E), and the side-area
  !>     vector, S, where S = d(s1,s2), where d is
  !>     the area of the cell face and (s1,s2) is the outward
  !>     pointing cell-face normal vector, pointing from cell L to cell R.
  !>     The solution is rotated to the S vector frame of reference.
  !>     The EXACTRS subroutine also needs the xphi and PRESSURE subroutines
  !>     listed below.
  !>     The flux is returned in vector US.
  subroutine ExactRiemannSolver(ndimL, UL, UR, S, US)
    integer, intent(in) :: ndimL
    real, dimension(1:ndimL), intent(in) :: UL, UR
    real, dimension(1:nbDim), intent(in) :: S
    real, dimension(1:ndimL), intent(out) :: US

    integer I,conv;
    real :: s1,s2,RGAMMA, SM,cl,cr,ql,qr,pl,pr,d,Rrl,Rrr, &
         ps,ut,vt,uu,vv, ML,MR,rconst,lconst,sv,newmr,newml,t1,t2,pi,qi,ri, &
         RSTAR,CSTAR,EEXP,SEXP,alpha,eps1,eps2,pn,ei,rdd;

    RGAMMA = state%model%kappa;
    conv = 0;
    eps1 = 1.0E-6;
    eps2 = 1.0E-6;
    alpha = 1.0E0;
    d = (S(1)*S(1)+S(2)*S(2))**0.5;
    rdd = 1.0E0/(d+1.0E-20);
    s1 = S(1)*rdd;
    s2 = S(2)*rdd;
    Rrl = 1.0E0/UL(1);
    Rrr = 1.0E0/UR(1);
    ql = (s1*UL(2)+s2*UL(3))*Rrl;
    qr = (s1*UR(2)+s2*UR(3))*Rrr;


    !print*,'##',s1,s2,Rrl,Rrr,ql,qr;

    pl = pressure(ndimL, UL);
    pr = pressure(ndimL, UR);
    cl = (RGAMMA*pl*Rrl)**0.5;
    cr = (RGAMMA*pr*Rrr)**0.5;
    ps = (pl+pr)*0.5;

    !print*
    !print*,'pl=',pl
    !print*,'pr=',pr
    !print*,'cl=',cl
    !print*,'cr=',cr

    !!print*,'##',pl,pr,cl,cr,ps,ql;

    !    START OF EXACT RIEMANN SOLVER
    rconst = UR(1)*pr**0.5;
    lconst = UL(1)*pl**0.5;
    MR = rconst*xphi(ps/pr);
    ML = lconst*xphi(ps/pl);


    !print*,'MR=',MR, rconst
    !print*,'ML=',ML, lconst

    !     Godunov's iteration
10  continue
    I=0;
20  continue
    pn = (pr/MR+pl/ML+ql-qr)/(1.0E0/ML+1.0E0/MR);
    pn = alpha*max(eps1,pn)+(1.0E0-alpha)*ps;
    newmr = rconst*xphi(pn/pr);
    newml = lconst*xphi(pn/pl);
    t1 = abs(newmr-MR);
    t2 = abs(newml-ML);
    if( max(t1,t2) < eps1  .or.  alpha < eps2 ) conv=1;
    ps = pn;
    MR = newmr;
    ML = newml;
    I = I+1;

    !write(*,'(a6,i5,6es14.6)')'####',I,ps,MR, ML, s1, s2

    if( I <20 .and. conv == 0 )goto 20

    alpha = alpha*0.8;
    if( conv == 0) goto 10
    SM = (pl-pr+ML*ql+MR*qr)/(ML+MR);

    !!print*,'##@!!',SM,pl,pr,ML,ql,MR,qr

    !write(*,'(a6,6es20.12)')'?????', SM, ps/pr
    !print*,'ps=',ps
    !print*,'ql=',ql
    !print*,'qr=',qr

!    Now find what is happening at the interface
    if (SM < 0.0E0) then
       !      RIGHT OF THE CONTACT
       ut = UR(2)*Rrr-qr*s1;
       vt = UR(3)*Rrr-qr*s2;
       if ((ps/pr) > 1.0E0) then
          !      RIGHT TRAVELLING SHOCK
          sv = qr+MR*Rrr;
          if (sv < 0.0E0)then
             !         SUPERSONIC FROM RIGHT TO LEFT
             ri = UR(1);
             qi = qr;
             pi = pr;
          else
             !         BETWEEN CONTACT AND RIGHT MOVING SHOCK
             ri = MR/(sv-SM);
             qi = SM;
             pi = ps;
          endif
       else
          !      RIGHT TRAVELLING EXPANSION
          EEXP=(qr+cr);
          if (EEXP < 0.0E0) then
             !         SUPERSONIC FROM RIGHT TO LEFT
             ri = UR(1);
             qi = qr;
             pi = pr;
          else
             RSTAR=(ps*(UR(1)**RGAMMA)/pr)**(1.0E0/RGAMMA);
             !!CSTAR=RGAMMA*ps/RSTAR**0.5;
             CSTAR=(RGAMMA*ps/RSTAR)**0.5;
             SEXP=(CSTAR+SM);
             if (SEXP > 0.0E0)then
                !            BETWEEN CONTACT AND START OF EXPANSION
                ri = RSTAR;
                qi = SM;
                pi = ps;
             else
                !            IN THE EXPANSION WAVE
                qi=(((RGAMMA-1.0E0)*qr*0.5)-cr)* 2.0E0/(RGAMMA+1.0E0);
                ri=1.0E0+(RGAMMA-1.0E0)*(qi-qr)/(2.0E0*cr);
                ri=UR(1)*ri*(2.0E0/(RGAMMA-1.0E0));
                pi=(pr/(UR(1)**RGAMMA))*(ri**RGAMMA);
             endif
          endif
       endif
    else
       !      LEFT OF THE CONTACT
       ut = UL(2)*Rrl-ql*s1;
       vt = UL(3)*Rrl-ql*s2;

!!       write(*,'(a6,6es18.10)')'!ZZ!',ps/pl, ps/pl - 1.0E0, ps/pl - 1.0D0,ps/pl - 1.

       if (ps/pl >= 1.0E0) then
          !        LEFT TRAVELLING SHOCK
          sv = ql-ML*Rrl;

          if (sv > 0.0E0)then
             !          SUPERSONIC FROM LEFT TO RIGHT
             qi = ql;
             ri = UL(1);
             pi = pl;
          else
             !          BETWEEN SHOCK AND CONTACT
             pi = ps;
             qi = SM;
             ri =  ML/(SM-sv);
          endif
       else
          !       LEFT TRAVELLING EXPANSION
          SEXP=ql-cl;

!!          write(*,'(a6,6es20.12)')'!AA!',ps/pl, SEXP

          if (SEXP > 0.0E0)then
             !          Supersonic from left to right
             ri = UL(1);
             pi = pl;
             qi = ql;
          else
             RSTAR=(ps*(UL(1)**RGAMMA)/pl)**(1.0E0/RGAMMA);
             !!CSTAR=RGAMMA*ps/RSTAR**0.5;
             CSTAR=(RGAMMA*ps/RSTAR)**0.5;
             EEXP=SM-CSTAR;
             !       if(ielem .eq. 79) write (*,*) 'SEXP=',SEXP,EEXP

!          write(*,'(a6,6es14.6)')'!..!',ps/pl, SEXP,EEXP
!          print*,'RSTAR=',RSTAR
!          print*,'CSTAR=',CSTAR
!          print*,'ps=',ps
!          print*,'cl=',cl
!          print*,'ql=',ql
!          print*,'ga=',RGAMMA

             if (EEXP > 0.0E0)then
                !             In the expansion wave
                qi=((RGAMMA-1.0E0)*ql*0.5+cl)*2.0E0/(RGAMMA+1.0E0);
                ri=1.0E0-((RGAMMA-1.0E0)/(2.0E0*cl))*(qi-ql);

!                print*,'ri:',ri,qi-ql

                ri=UL(1)*ri**(2.0E0/(RGAMMA-1.0E0));
                pi=(pl/(UL(1)**RGAMMA))*(ri**RGAMMA);

!                print*,'qi=',qi,ql
!                print*,'ri=',ri
!                print*,'pi=',pi

             else
                !             BETWEEN EXPANSION AND CONTACT
                ri = RSTAR;
                qi = SM;
                pi = ps;
                !    if(ielem .eq. 79) write (*,*) 'ri,qi,pi=',ri,qi,pi
             endif
          endif
       endif
    endif


    uu = s1*qi+ut;
    vv = s2*qi+vt;

    !     if(ielem .eq. 79) write (*,*) ri,uu,vv,pi,qi,'&&&&&'
    !     Find the flux at the interface and scale by the face area, d
    US(1) = ri;
    US(2) = ri*uu;
    US(3) = ri*vv;
    ei = pi/(RGAMMA-1.0E0)+0.5*ri*(uu*uu+vv*vv);
    US(4) = ei;

  end subroutine ExactRiemannSolver

  function xphi(x)
    real :: xphi
    real, intent(in) :: x

    if (x < 0.999) then
        xphi =  state%model%kappa1*(1.0D0-x) /(2.3664319*(1.0E0- x**0.14285714) );
    else
        xphi = 1.2*x+0.2**0.5 ;
    endif
    return
  end function xphi


  function pressure(ndimL, U)
    integer, intent(in) :: ndimL
    real :: pressure
    real, dimension(1:ndimL), intent(in) :: U
    pressure = state%model%kappa1*(U(4)-(U(2)*U(2)+U(3)*U(3) )/2.0E0/U(1) );
    return
  end function pressure
  subroutine TimeDependentIOBC(Qdof, ndimL, wD)
    ! change of the BC state vectors for time-depended inflow/outflow BC
    integer, intent(in) :: Qdof, ndimL
    real, dimension(1:Qdof,1:ndimL), intent(inout) :: wD
    real, dimension(1:Qdof) :: v1, v2, p

    real :: pi = 3.14159265359
    real :: omega = 0.5      ! omega = frequency
    real :: degree_max = 10   ! alpha_max = amplitude of angle
    real :: alpha


    ! OSCILLATIONS OF ANGLE OF ATTACK AT INFLOW/OUTFLOW
    ! actual angle amplitude
    !alpha = degree_max/180*pi * sin(2*pi*omega * (state%Ttime+state%time%tau(1)))
    !!print*,'$$$',state%time%iter,state%Ttime+state%time%tau(1), alpha

    !v1(1:Qdof) = wD(1:Qdof,2) * cos(alpha) - wD(1:Qdof,3) * sin(alpha)
    !v2(1:Qdof) = wD(1:Qdof,2) * sin(alpha) + wD(1:Qdof,3) * cos(alpha)
    !wD(1:Qdof,2) = v1(1:Qdof)
    !wD(1:Qdof,3) = v2(1:Qdof)

    ! OSCILLATIONS OF OUTLET PRESSURE
    p(1:Qdof) = (wD(1:Qdof,4) - &
         (wD(1:Qdof,2)*wD(1:Qdof,2) + wD(1:Qdof,3)*wD(1:Qdof,3))/2/wD(1:Qdof,1)) &
         *state%model%kappa1

    p(1:Qdof) = p(1:Qdof)*(1+0.1*sin(2*pi*omega*(state%time%ttime +state%time%tau(1)) ))

    wD(1:Qdof,4) = p(1:Qdof)/state%model%kappa1 + &
         (wD(1:Qdof,2)*wD(1:Qdof,2) + wD(1:Qdof,3)*wD(1:Qdof,3))/2/wD(1:Qdof,1)

    !print*,'$$$',state%time%iter,state%Ttime+state%time%tau(1), p(1:nbDim)

  end subroutine TimeDependentIOBC

  !> rsetting of BC, a use of the linearized Riemann problem
  subroutine SetBCCharacteristic(grid, Qdof, ndimL, wi, wD, w_BC, n, xi, elem )
    class(mesh), intent(in) :: grid
    integer, intent(in) :: Qdof, ndimL
    real, dimension(1:Qdof,1:ndimL), intent(in) :: wi, wD
    real, dimension(1:Qdof,1:ndimL), intent(out) ::  w_BC
    real, dimension(1:nbDim), intent(in) :: n
    type(element), intent(in) :: elem
    real, dimension(1:nbDim) :: nn, xi
    real, dimension(1:ndimL,1:ndimL) :: t, t1
    real, dimension(1:ndimL) :: qi, qD, qj, dp, alpha, beta, omega
    real :: size, wRP(ndimL)
    integer :: i,j

    if(ndimL /= 4) then
       print*,'Bad implementation of SetBCCharacteristic'
       stop
    endif
    size = (dot_product(n, n))**0.5
    nn(1:nbDim) = n(1:nbDim)/size

    do i=1, Qdof
       wRP(1:ndimL) = wi(i, 1:ndimL)
       call ComputeEigenValsVec(grid, ndimL, wRP(1:ndimL),nn(1:nbDim),xi(1:nbDim), t1, t, dp, elem)

       call RotationForward(ndimL, wi(i,1:ndimL), nn, qi);
       call RotationForward(ndimL, wD(i,1:ndimL), nn, qD);

       do j=1,ndimL
          alpha(j) = dot_product(t1(j,1:ndimL), qi(1:ndimL) )
          beta(j) = dot_product(t1(j,1:ndimL), qD(1:ndimL) )
       enddo

       do j=1,ndimL
          if(dp(j) >= 0) then
             omega(j) = alpha(j)
          else
             omega(j) = beta(j)
          endif
       enddo

       do j=1,ndimL
          qj(j) = dot_product(t(j,1:ndimL), omega(1:ndim) )
       enddo
       call RotationBackward(ndimL, qj, nn, w_BC(i,1:ndimL) )

    enddo

  end subroutine SetBCCharacteristic

  !> Transform of the conservative variables to the physical ones for NSe
  subroutine Transform_W2Q_NSe(Qdof, wi, qi)
    integer, intent(in) :: Qdof
    real, dimension(1:4, 1:Qdof), intent(in) :: wi
    real, dimension(1:4, 1:Qdof), intent(inout) :: qi

    qi(1, 1:Qdof) = wi(1, 1:Qdof)
    qi(2, 1:Qdof) = wi(2, 1:Qdof) / wi(1, 1:Qdof)
    qi(3, 1:Qdof) = wi(3, 1:Qdof) / wi(1, 1:Qdof)

    qi(4, 1:Qdof) =  state%model%kappa1 * ( wi(4, 1:Qdof) &
         - 0.5* (wi(2, 1:Qdof)*wi(2, 1:Qdof) + wi(3, 1:Qdof)*wi(3, 1:Qdof) ) / wi(1, 1:Qdof) )

  end subroutine Transform_W2Q_NSe

  !> Transform of the physical variables to the conservative ones for NSe
  subroutine Transform_Q2W_NSe(Qdof, qi, wi)
    integer, intent(in) :: Qdof
    real, dimension(1:4, 1:Qdof), intent(in) :: qi
    real, dimension(1:4, 1:Qdof), intent(inout) :: wi

    wi(1, 1:Qdof) = qi(1, 1:Qdof)
    wi(2, 1:Qdof) = qi(2, 1:Qdof) * qi(1, 1:Qdof)
    wi(3, 1:Qdof) = qi(3, 1:Qdof) * qi(1, 1:Qdof)

    wi(4, 1:Qdof) =    qi(4, 1:Qdof)/ state%model%kappa1  &
         + 0.5 * (qi(2, 1:Qdof)*qi(2, 1:Qdof) + qi(3, 1:Qdof)*qi(3, 1:Qdof) ) * wi(1, 1:Qdof)

  end subroutine Transform_Q2W_NSe

end module model2DNS
