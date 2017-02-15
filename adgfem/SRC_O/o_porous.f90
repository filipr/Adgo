!> algebraic functions for the porous media flow model by M. Kuraz
!> data parameters
module porous_fnc
  public :: init_porous
  public :: vangen
  
  type, private :: soilpar_str
    real :: alpha  ! Van Genuchten 
    real :: n      ! Van Genuchten
    real :: m      ! Van Genuchten
    real :: Thr    ! theta_r
    real :: Ths    ! theta_s
    real :: Ss     ! storativity
    real :: Ks     ! saturated conductivity
  end type soilpar_str
  
  type(soilpar_str), private, dimension(:), allocatable :: soilpar
  
  !> dynamic viscosity
  real, parameter, private :: mu=1.62e-7
  !> density of water
  real, parameter, private :: rho=1000.0
  !> "configuration of the void space" exact value definition should be improved
  !real, parameter, private :: beta=1.0
  ! VD
  real, parameter, private :: beta=0.0
  
  
  contains
  
    !> assigns soil parameters and allocates structures, at this moment this function is just simple initialization of values for the case study of dam seepage, can be always improved:)
    subroutine init_porous_coeffs()
      implicit none
      integer :: layers
      
      layers = 3 ! in current case study we have 3 different materials
      
      allocate(soilpar(layers))
      
      !!!!!!!!!!!!!!!!!!!!!!
      !units: m, days
      !!!!!!!!!!!!!!!!!!!!!
      
      !layer 1 = gravel
      ! original values
      !soilpar(1)%alpha = 100.0
      !soilpar(1)%n = 2.0
      !soilpar(1)%m = 0.626

      ! New values by Michal
      soilpar(1)%alpha = 2.0
      soilpar(1)%n = 1.41
      soilpar(1)%m = 0.291

      soilpar(1)%Ks = 7.128
      soilpar(1)%ths = 0.43
      soilpar(1)%thr = 0.01   ! 0.
      soilpar(1)%Ss = 1.0E-02
      
      !layer 2 = clay
      soilpar(2)%alpha = 0.8
      soilpar(2)%n = 1.2
      soilpar(2)%m = 0.167
      soilpar(2)%Ks = 0.048
      soilpar(2)%ths = 0.38
      soilpar(2)%thr = 0.06
      soilpar(2)%Ss = 1.0E-02
      
      !layer 3 = silt clay
      soilpar(3)%alpha = 2.0
      soilpar(3)%n = 1.41
      soilpar(3)%m = 0.291
      soilpar(3)%Ks = 0.108
      soilpar(3)%ths = 0.45
      soilpar(3)%thr = 0.067
      soilpar(3)%Ss = 1.0E-02
      
      !soilpar(1:3)%Ss = 0. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

      ! soilpar(1:3)%alpha = 100.0
      ! soilpar(1:3)%n = 2.0
      ! soilpar(1:3)%m = 0.626
      ! soilpar(1:3)%Ks = 7.128
      ! soilpar(1:3)%ths = 0.43
      ! soilpar(1:3)%thr = 0.0
      ! soilpar(1:3)%Ss = 1.0E-02
      
      ! soilpar(1:3)%alpha = 2.0
      ! soilpar(1:3)%n = 1.41
      ! soilpar(1:3)%m = 0.291
      ! soilpar(1:3)%Ks = 0.108
      ! soilpar(1:3)%ths = 0.45
      ! soilpar(1:3)%thr = 0.067
      ! soilpar(1:3)%Ss = 1.0E-02
    
       ! soilpar(1:3)%alpha = 2.0
       ! soilpar(1:3)%n = 4.
       ! soilpar(1:3)%m = 2.
       ! soilpar(1:3)%Ks = 1.0
       ! soilpar(1:3)%ths = 0.5
       ! soilpar(1:3)%thr = 0.1
       ! soilpar(1:3)%Ss = 1.0E-00
    
    end subroutine init_porous_coeffs


    !> \brief Van Genuchten relation \f$ \theta = f(pressure) \f$
    !>  \f$ \theta_e = \frac{1}{(1+(\alpha*h)^n)^m} \f$
    !> water content is considered as absolute value not the relative one \n
    !> see \f$ \theta_e = \frac{\theta - \theta_r}{\theta_s-\theta_r} \f$
    !>
    function vangen(u, layer, z) result(theta)
      implicit none
      !> solution
      real, intent(in) :: u
      !> material id
      integer, intent(in) :: layer
      !> geodetic head
      real, intent(in) :: z
      !> resulting water content
      real :: theta

      real :: a,n,m, theta_e, h
      
      h = u - z

      
      a = soilpar(layer)%alpha
      n = soilpar(layer)%n
      m = soilpar(layer)%m
      

      if (h >=0.0) then
	  theta = soilpar(layer)%Ths
	  RETURN
      else
	  theta_e = 1/(1+(a*(abs(h)))**n)**m
	  theta = theta_e*( soilpar(layer)%Ths - soilpar(layer)%Thr ) +soilpar(layer)%Thr
      end if

    end function vangen
    
    !> \brief so-called retention water capacity, it is a derivative of the retention curve function
    !> \f$ E(h) = C(h) + \frac{\theta(h)}{\theta_s}S_s \f$
    !> where
    !> \f$ C(h) = \left\{ \begin{array}{l l}\frac{m n \alpha  (-h \alpha )^{-1+n}}{\left(1+(-h \alpha )^n\right)^{1+m}}(\theta_s - \theta_r), 
    !> & \quad \mbox{$\forall$ $h \in (-\infty, 0 )$}\\ 0, & \quad \mbox{$\forall$ $h \in \langle 0, + \infty )$}\\ \end{array} \right. \f$
    !> and 
    !> \f$ \theta(h) = \left\{ \begin{array}{l l} \frac{\theta_s -\theta_r}{(1+(-\alpha h)^n_{vg})^m_{vg}} + \theta_r,  
    !> & \quad \mbox{$\forall$ $h \in (-\infty, 0 )$}\\ \theta_S, & \quad \mbox{$\forall$ $h \in \langle 0, + \infty )$}\\ \end{array} \right. \f$
    !>
    function capacity(u, layer, z) result(E)

      implicit none
      !> solution
      real, intent(in) :: u
      !> material id
      integer, intent(in) :: layer
      !> geodetic head
      real, intent(in) :: z
      !> resulting system capacity (elasticity)
      real :: E

      real :: C, a, m, n, tr, ts, h
          
      h = u - z

      if (h < 0) then
	a = soilpar(layer)%alpha
	n = soilpar(layer)%n
	m = soilpar(layer)%m
	tr = soilpar(layer)%Thr
	ts = soilpar(layer)%Ths
	C = a*m*n*(-tr + ts)*(-(a*h))**(-1 + n)*(1 + (-(a*h))**n)**(-1 - m)
      else
	E = soilpar(layer)%Ss
	RETURN
      end if

      E = C + vangen(u, layer, z)/soilpar(layer)%Ths*soilpar(layer)%Ss

      ! VD
      !E = soilpar(layer)%Ss
      !write(*,'(a8, 20es12.4)') 'params:', soilpar(layer)%Ss, E , E/ soilpar(layer)%Ss
      !write(*,'(a8, 20es12.4)') 'params:', soilpar(layer)%Ss, a, n, m, tr, ts, soilpar(layer)%Ks

    end function capacity
    
    
    !> \brief Mualem's function for unsaturated hydraulic conductivity with van Genuchten's water content substitution
    !> \f$   K(h) = \left\{ \begin{array}{l l} K_s\frac{\left( 1- (-\alpha h)^{n_{vg}m_{vg}} \left( 1+ (-\alpha h)^{n_{vg}} 
    !> \right)^{-m_{vg}} \right)^2}{\left(1+(-\alpha h)^{n_{vg}} \right)^{\frac{m_{vg}}{2}}},  &  \mbox{$\forall$  
    !> $h \in$ $(-\infty,0)$}\\ K_s, & \mbox{$\forall$   $h \in$ $\langle 0, +\infty)$}\\ \end{array} \right. \f$
    !>
    function conduct(u, layer, z) result(K)
    
      implicit none
      !> solution
      real, intent(in) :: u
      !> material id
      integer, intent(in) :: layer
      !> geodetic head
      real, intent(in) :: z
      !> resulting hydraulic conductivity
      real :: K
      real :: a,n,m,h
      real :: K1, ah

      h = u - z
      
      if (h >  0) then
         K = soilpar(layer)%Ks
      else
         a = soilpar(layer)%alpha
         n = soilpar(layer)%n
         m = soilpar(layer)%m
         
         K =  (1 - (-(a*h))**(m*n)/(1 + (-(a*h))**n)**m)**2/(1 + (-(a*h))**n)**(m/2.0) * soilpar(layer)%Ks
         !ah = -a*h

         !print*,'#DE#:',( 1 + (ah)**n)**(m)
         !print*,'#DE#:',( 1 + (ah)**n )**(m/2), ah, n*m
         !print*,'????:',  ( 1  + ah**(n*m) / ( 1 + (ah)**n)**(m) )**2 

         !K1 = ( 1 - ah**(n*m) / ( 1 + ah**n)**m )**2 / ( 1 + ah**n )**(m/2) * soilpar(layer)%Ks
         
         !if( K /soilpar(layer)%Ks  < 1E-10) &
         !     write(*,'(a8, 26es12.4)') ' K = ', K, K1, soilpar(layer)%Ks,  K/  soilpar(layer)%Ks, &
         !     h, u, z
      end if

      ! VD
      !K = soilpar(layer)%Ks
      
    end function conduct
    
    function forch_conduct(u, gradu, layer, z) result(kappa)
      implicit none
      !> solution
      real, intent(in) :: u
      !> L2 norm of the solution gradient
      real, intent(in) :: gradu
      !> material id
      integer, intent(in) :: layer
      !> geodetic head
      real, intent(in) :: z
      
      !> resulting kappa conductivity for Forchheimer equation
      real :: kappa
      
      real :: K
      
      K = conduct(u, layer, z)

      
      ! ERROR ?, exchanged a_0 and a_1 !!!
      !kappa = 2.0/(1.0/K + sqrt( (rho*beta/(mu*K))*(rho*beta/(mu*K)) + 4.0/K*gradu))

      ! better computer arithmetic, works for K= 0
      !kappa = 2.0 * K /(1.0 + sqrt( (rho*beta/(mu))*(rho*beta/(mu)) + 4.0*K*gradu))

      ! VD CORRECTED
      kappa = 2.0 * K /(1.0 + sqrt(1 +  4. * (rho*beta/ mu)*K *gradu))
    
    end function forch_conduct

end module porous_fnc


!> definition of model of porous media flow
module porous_mod

 ! use main_data
  use model_mod
!  use f_mapping
!  use mesh_oper
!  use define_state
!  use blocks_integ
!  use model3DNS
!  use model2DNS
!  use modelTurb2e
!  use modelFE
!  use modelLaplace
  use porous_fnc

  implicit none

 type, EXTENDS( Model_t ), public :: Porous_t
!   real :: kappa, kappa1
!   real :: Pr


 contains
   procedure :: init => initPorous

 end type Porous_t


 contains
 !> initialization of the NS case
 !> isca, t2 not used
  subroutine initPorous( this, Re, isca, t2)
    class (Porous_t), intent(inout) :: this
    real, intent(in) ::Re
    integer, intent(in), optional :: isca
    real, intent(in), optional :: t2

    !stop 'initPorous not implemented yet'
    !??? ndim?
    this%ndim = 1 ! is not sufficient for the global setting, set in readModelData !!
    this%convective = .false.

    this%subdomainRHS = .false.
    this%rhsTestFunDerivative = 0 ! rhs: \int f*phi dx
    ! is the exact solution a priori known?
    this%known_sol = .false.

    this%Re = Re
    this%icase = isca
    this%param1 = t2

    this%ireac = 0

    if( Re > 0.) then
       this%Re1  = 1./ this%Re
    else
       this%Re1 = 0.
    endif

    ! initialization of the data coefficients
    call init_porous_coeffs()
    
    print*,' # Porous media flow problem, Forchheimer'
    
    !has to be modified !!!
    select case (this%icase)
    case(1)   ! damp (HRAZ),  linear test case
       this%idiff = 2
       this%iconv = 0
       this%iexact = 2
       this%conv_rat = 1.0

       !this%varying_time_term = .true. !  here only for the test
 
    case(2)   ! Forcheimer  damp (HRAZ),  Forchheimer 2-term law, only test case
       this%idiff = 2
       this%iconv = 0
       this%iexact = 2
       this%conv_rat = 1.0

       this%varying_time_term = .true. ! term in fron of the time derivative
 
    case(3)   ! Forcheimer  damp (HRAZ),  Forchheimer 2-term law, REAL data by M. Kuraz
       this%idiff = 3
       this%iconv = 0
       this%iexact = 2  !change later
       this%conv_rat = 1.0

       this%varying_time_term = .true. ! term in fron of the time derivative
 
    case(4)   ! test nonlinear case:  2u u_t - (u^2 u_x)_x = 0
       this%idiff = 4
       this%iconv = 0
       this%iexact = 3  
       this%conv_rat = 1.0

       this%varying_time_term = .true. ! term in fron of the time derivative
 
    case(5:71)
       print*,' UNKNOWN TYPE of scalar%icase !!!'
       print*,'Adddefinitions in:'
       print*,'                    o_porous.f90:    subroutine initPorous '
       print*,'                    modelPorous.f90: function Eval_Diff_Porous_Coeffs'
       print*,'                    problem.f90:      subroutine Setting_of_space_variable_coeffs'
       stop

    case(72)   ! porous media flow,  Forchheimer 2-term law
       this%idiff = 13
       this%iconv = 0
       this%iexact = 65
       this%conv_rat = 1.0
    case(73:)
       print*,' UNKNOWN TYPE of this%isca !!!'
       stop
       
    end select
    


!    this%kappa = 1.4
!    this%kappa1 = this%kappa - 1.
!    this%Pr = 0.72   ! Prandtl number (almost all gases)
!
!
!    if ( Re == 0.) then
!        print*,' # Compressible Euler equations'
!        this%Re1 = 0.
!     elseif  ( Re > 0.) then
!        this%Re1 = 1./this%Re
!        print*,' # Compressible Navier-Stokes equations, Re=',this%Re
!
!     else
!        print*,'# Reynolds number is negative',this%Re,' STABILIZATION !!!'
!        !stop
!     endif

  end subroutine initPorous

end module porous_mod



