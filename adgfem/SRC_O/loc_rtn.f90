!> Raviart-Thomas-Nedelec (RTN) finite elements constructed from the
!> local basis on each element
module loc_rav_tho_ned
  use geometry
  use integration
  use lapack_oper

  implicit none

  !> local RTN basis 
  type :: basis_rtn_fe
     logical :: defined
     integer :: Fdeg
     integer :: Fdof   ! degree of freedom
     !type(loc_rtn_fe), dimension(:), pointer :: phi
     integer, dimension(:,:), allocatable :: ipsi
     !integer, pointer, dimension(:,:) :: ipsi
  end type basis_rtn_fe

  public :: SetRTNdof
  public :: Init_Loc_RTN
contains
  
  !> return dof of RTN finite element of degref Fdeg
  function SetRTNdof(Fdeg)
    integer :: SetRTNdof
    integer, intent(in) :: Fdeg

    SetRTNdof = (Fdeg+1)*(Fdeg+3)
    
  end function SetRTNdof

  !> init basis of the loc_RTN which consists of the composition of standard
  !> DG FE basis functions on elem in the form:
  !>
  !> \f$\psi = (\varphi_i, \varphi_j)^{\rm T} + (x_1,x_2)^{\rm T}\varphi_k\f$,
  !> \f$\varphi_i, \varphi_j \f$ are polynomials from  \f$S_p(K)\f$ of degree 
  !> \f$ \leq Fdeg\f$ and 
  !> \f$\varphi_k, \varphi_j \f$ is polynomial from  \f$S_p(K)\f$ of degree 
  !> \f$ = Fdeg\f$. Hence we need only indexes
  subroutine  Init_Loc_RTN(loc_RTN, Fdeg)
    type(basis_rtn_fe), intent(inout) :: loc_RTN
    integer, intent(in) :: Fdeg
    integer :: Fdof, i, dof

    Fdof = SetRTNdof(Fdeg)
    if(nbDim == 2) dof = (Fdeg+1)*(Fdeg+2) /2 
    if(nbDim == 3) dof = (Fdeg+1)*(Fdeg+2)*(Fdeg+2) /6 

    loc_RTN%Fdeg = Fdeg
    loc_RTN%Fdof = Fdof

    allocate(loc_RTN%ipsi(1:Fdof,1:2) )


    ! test functions with \phi_i /= 0
    do i=1,dof
       loc_RTN%ipsi(i, 1) = 1
       loc_RTN%ipsi(i, 2) = i
    enddo
    
    ! test functions with \phi_j /= 0
    do i=1,dof
       loc_RTN%ipsi(dof+i, 1) = 2
       loc_RTN%ipsi(dof+i, 2) = i
    enddo
    
    ! test functions with \phi_k /= 0
    do i=1,Fdeg+1
       loc_RTN%ipsi(2*dof+i, 1) = 3
       loc_RTN%ipsi(2*dof+i, 2) = dof-Fdeg-1 + i
    enddo

    !write(*,'(a8,80i3)') 'ipsi(1):',loc_RTN%ipsi(:,1)
    !write(*,'(a8,80i3)') 'ipsi(2):',loc_RTN%ipsi(:,2)

    loc_RTN%defined = .true.
    
  end subroutine Init_Loc_RTN


end module loc_rav_tho_ned
