!> parameters and functions for wet steam source terms
module wet_steam_paramets
  use eval_sol

  implicit none
  public
 
  public :: RHS_WS

 
  ! constants
  real, parameter :: Pii = 3.14159265359
  real, parameter :: kb = 1.3806488E-23        ! Boltzmann's constant [J K^-1]
  real, parameter :: Rv = 461.52               ! gas constant of water vapor  [J K^-1 kg^-1]
  real, parameter :: mv = 2.9904610E-26        ! mass of vapor molecule [kg]
  
  ! nedouzavrene veliciny
  !----------------------
  real ::  L = 2260.E3  ! latent heat of condensation [kJ kg^-1] (Wikipedia)
  real :: r_l = 960.    ! water density (material property) at approx. 100°C [kg/m^3]

  ! http://www.peacesoftware.de/einigewerte/wasser_dampf_e.html - saturated steam T = 400
  real ::  cp = 2.2160192812928E3     ! specific heat of vapor at constant pressure [J K^-1 kg^-1]
  !real ::  cv = 1.6415527886029E3     ! specific heat of vapor at constant volume [J K^-1 kg^-1]
  !real ::  cv = cp-Rv
  real :: eta_v = 1.3192492103853E-5  ! dynamic viscosity of water vapor [Pa s]
  real :: lambda_v = 0.0283471820364  ! heat conduction coefficient of vater wapor [W m^-1 K^-1]
  real :: p_s = 245.753186304E3       ! pressure of saturated water vapor [Pa], should be a function of...
  real :: theta_s = 390.              ! temperature of saturated water vapor [K], should be a function of... - estimated from state equatin
  
  real :: sigma = 58.85E-3            ! water surface tension at 370 K [N m^-1], http://en.wikipedia.org/wiki/Surface_tension 
  ! correction values or functions (??needed??)
  real :: beta  ! 
  real :: p_cor !

  contains      ! source terms for wet steam equations

    !> Computes density of vapor in wet steam flow (not material property)
    subroutine Rho_vap(omega,rho,rho_v)
      real, intent(in) :: omega,rho
      real, intent(inout) :: rho_v

      rho_v = (1-omega)*rho
    end subroutine Rho_vap

    !> vapor temperature
    subroutine VaporTheta(p,rho_v,theta_v)
      real, intent(in) :: p, rho_v
      real, intent(inout) :: theta_v
 
      theta_v = p/(Rv*rho_v)
    end subroutine VaporTheta

    !> vapor pressure
    subroutine VaporPressure(kappa,omega,E,u,v,rho,L,p_v)
      real, intent(in) :: kappa,omega,e,u,v,rho,L
      real, intent(inout) :: p_v

      p_v = (kappa-1)*(1-omega)/(1+omega*(kappa-1))*(E-0.5*(u**2 + v**2)+rho*L*omega)
    end subroutine VaporPressure

    

    !> specific heat of vapor at constant volume [J K^-1 kg^-1]
    subroutine Cv_vapor(Rv,cp,cv)
      real, intent(in) :: Rv, cp
      real, intent(inout) :: cv
      
      cv = cp-Rv
    end subroutine Cv_vapor

    !!> Computes density of water in wet steam flow (not material property)
    !subroutine Rho_liq(omega,rho,rho_l)
    !  real, intent(in) :: omega,rho
    !  real, intent(inout) :: rho_l

    !  rho_l = omega*rho
    !end subroutine Rho_liq

    !> Computes the average water drop radius
    subroutine AverDropRadius(omega,Q2,Q0,r_aver)
      real, intent(in) :: omega, Q2, Q0
      real, intent(inout) :: r_aver
      
      if(omega <= 1.E-6) then
        r_aver = 0.
      else
        if(Q0 <=0.) then
          write(*,*) 'Problem: Q0 <= 0'
          stop
        endif
        r_aver = sqrt(Q2/Q0)
      endif
    end subroutine AverDropRadius

    !> Computes the critical drop radius
    subroutine CritRadius(theta_v,p_v,p_s,rho_l,sigma,r_crit)
      real, intent(in) :: theta_v,p_v,p_s,rho_l,sigma
      real, intent(inout) :: r_crit
      
      real :: S
      

      S = log(p_v/p_s)
      if (S <= 0.) then
        write(*,*) "pv < ps, can't compute S"
        stop
      end if
      
      r_crit = 2*sigma/(r_l*Rv*theta_v*S)

    end subroutine CritRadius

    !> Computes Knudsen's number
    subroutine KnudsenNumber(ni_v,theta_v,p_v,r_aver,Kn)
      real, intent(in) :: ni_v,theta_v,p_v,r_aver
      real, intent(inout) :: Kn

      if(r_aver <=0.) then
        write(*,*) 'Problem in Knudsen: r_aver <= 0'
        stop
      endif
      Kn = ni_v*sqrt(2*Pii*Rv*theta_v)/(4*r_aver*p_v)
    end subroutine KnudsenNumber

    !> Computes the parametrisation of time derivation of the water drop radius
    subroutine RadGrowth(lambda_v,theta_s,theta_v,r_aver,r_crit,r_l,L,Kn,r_dot)
      real, intent(in) :: lambda_v,theta_s,theta_v,r_aver,r_crit,r_l,L,Kn
      real, intent(inout) :: r_dot

      if(r_aver <=0.) then
        write(*,*) 'Problem in RadGrowth: r_aver <= 0'
        stop
      endif
      r_dot = lambda_v*(theta_s-theta_v)*(r_aver-r_crit)/(L*r_l*(1.+3.18*Kn)*r_aver*r_aver)         
    end subroutine RadGrowth

    !> Computes nucleation rate J
    subroutine NuclRate(sigma,rho_v,r_l,r_crit,theta_v,J)
      real, intent(in) :: sigma,rho_v,r_l,r_crit,theta_v
      real, intent(inout) :: J

      J = sqrt(2*sigma/(Pii*mv**3))*rho_v*rho_v/r_l*exp(-4*Pii*r_crit*r_crit*sigma/(3*kb*theta_v))
    end subroutine NuclRate


    !> Computes the right hand side vector of source terms for wet steam
    subroutine RHS_WS(elem,Qdof,x,f)
      type(element), intent(in) :: elem
      real, dimension(1:nbDim), intent(in) :: x
      integer, intent(in) :: Qdof
      real, dimension(1:Qdof,1:ndim), intent(inout) :: f

      real, dimension(:,:), allocatable :: wi   
      real, dimension(1:Qdof) :: omega, Q2, Q1, Q0, rho, u,v, theta_vapor
      real :: rho_v,r_crit,theta_v,theta_s,J,p_v,p_s,Kn,r_dot,r_aver
      integer :: i

      allocate(wi(1:Qdof,1:ndim))
      call Eval_w_Elem(elem,wi(1:Qdof,1:ndim))    !> evaluation of the state vector \f$w\f$ on element \f$ elem \f$  in 
                                                  !> integ nodes, i.e. recomputation of elem%w into wi in volume integ. nodes
      
      rho(:) = wi(:,1)
      u(:) = wi(:,2)
      v(:) = wi(:,3)
      !theta_vapor(:) = (wi(:)/rho(:) - ((u(:)**2+v(:)**2)/2.))/cv Predelat
      omega(:) = wi(:,5)/rho(:)
      Q2(:) = wi(:,6)/rho(:)
      Q1(:) = wi(:,7)/rho(:)
      Q0(:) = wi(:,8)/rho(:)
      
      do l=1,Qdof
        if(theta_vapor(l) <= 0.) then
          write(*,*) 'Absolute temperature lower than 0 in element'
          write(*,*) elem%i
        end if
      end do


      f(1:Qdof,1:4) = 0.  ! no source terms for NSe
      !! sourcei terms for wet steam part
      !do i=1,Qdof
      !  f(i,5) = 4./3.*Pii*CritRad(theta_v,p_v,p_s,rho_l,sigma,r_crit)**3*r_l*NuclRate() + 4*rho(:)*Pii*Q2(:)*RadGrowth()*Rho_l()
      !  f(i,6) = CritRad(theta_v,p_v,p_s,rho_l,sigma,r_crit)**2*NuclRate() + 2*w(:,1)*Q1(:)*RadGrowth()*Rho_l()
      !  f(i,7) = CritRad(theta_v,p_v,p_s,rho_l,sigma,r_crit)*NuclRate() + w(:,1)*Q0(:)*RadGrowth()
      !  f(i,8) = NuclRate()
      !end do
      write(*,*) 'SMRK'
      stop
      deallocate(wi)

    end subroutine RHS_WS


end module wet_steam_paramets
