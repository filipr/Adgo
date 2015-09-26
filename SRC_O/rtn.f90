!> Raviart-Thomas-Nedelec (RTN) finite element
module rav_tho_nedelec
  use geometry
  use integration
  use lapack_oper

  implicit none

  !> RTN basis functions in V_rule integration nodes of \f$ S_{hp} \f$
  !> derivatives of RTN basis functions in V_rule integration nodes of \f$ S_{hp} \f$
  type :: rtn_Vnode
     logical :: def
     real, dimension(:,:,:), allocatable :: phi    ! basis functions
     logical :: defder
     real, dimension(:,:,:,:), allocatable :: Dphi !derivatives of basis functions
  end type rtn_Vnode

  !> RTN finite element
  type :: rtn_fe
     integer :: deg
     integer :: dof   ! degree of freedom
     integer :: Qnum  ! index of volume quadrature
     integer :: Qdof  ! numbger of volume integ nodes
     integer :: Gnum  ! index of edge quadrature
     integer :: Gdof  ! numbger of edge integ nodes
     real, dimension(:,:,:,:), pointer :: phi        ! basis functions
     real, dimension(:,:,:,:,:), pointer :: Dphi     ! their derivatives
     real, dimension(:,:), allocatable :: MM     !  Momentum matrix
     real, dimension(:,:), allocatable :: MMinv  !  inverse of Momentum matrix
     type(rtn_Vnode), dimension(:), pointer :: Vnode  ! basis function in different V_rules
  end type rtn_fe


  public :: Compute_RTN_Momentums
  public :: EvalMomentumFace
  public :: EvalMomentumVolume
  public :: EvalMomentumVolumeD
  public :: Init_RTN_FE
  public :: Init_RTN_FE_Vnode
  public :: Init_RTN_FE_Vnode_Der
  public :: Set_RTN_integ_node

contains
  !> inicialization of RTN elements
  subroutine Init_RTN_FE(RTN, V_rule, G_rule)
    type(rtn_fe), intent(inout) :: RTN
    type(rtn_fe) :: RTNtest      !TEST
    type(volume_rule), intent(in) :: V_rule
    type(Gauss_rule), intent(in) :: G_rule
    real, dimension(:,:,:), allocatable :: x
    real, dimension(:), allocatable :: thetai     !TEST
    real, dimension(:,:,:), allocatable :: psihat, psi    !TEST 2
    real, dimension(:,:), allocatable :: B    !TEST 2
    integer :: dof, Qdof, Gdof, mdof, deg
    integer :: i,j,l, it, k

    Qdof = V_rule%Qdof
    RTN%Qdof = Qdof
    RTNtest%Qdof = Qdof            !TEST

    Gdof = G_rule%Qdof
    RTN%Gdof = Gdof
    RTNtest%Gdof = Gdof            !TEST


    allocate(x(0:3, 1:max(Qdof, Gdof), 1:2))
    !--------- ^^^  0= volume, 1,2,3= 1st, 2nd, 3rd face

    call Set_RTN_integ_node(V_rule, G_rule, Qdof, Gdof, x)

    deg = RTN%deg
    dof = (deg+1)*(deg+3)   !pro N=2
    RTN%dof = dof
    allocate (RTN%phi(1:dof, 1:2, 0:3, 1:max(Qdof, Gdof) ) )
    !---------------------------- ^^^  0= volume, 1,2,3= 1st, 2nd, 3rd face

    RTN%phi(1:dof, 1:2, 0:3, 1:max(Qdof, Gdof) ) = 0    !initialization

    allocate (RTNtest%phi(1:dof, 1:2, 0:3, 1:max(Qdof, Gdof) ) )   !TEST

    RTNtest%phi(1:dof, 1:2, 0:3, 1:max(Qdof, Gdof) ) = 0.    !TEST

    allocate (psihat(1:2, 0:3, 1:max(Qdof, Gdof) ), psi(1:2, 0:3, 1:max(Qdof, Gdof) ), B(1:2, 1:2) )   !TEST 2

    psihat(1:2, 0:3, 1:max(Qdof, Gdof) ) = 0.    !TEST 2
    psi(1:2, 0:3, 1:max(Qdof, Gdof) ) = 0.    !TEST 2
    B(1,1) = 0.
    B(1,2) = -1.
    B(2,:) = 1.
    do l=1,2
       !write(*,'(es12.4)') B(l, 1:2)
    enddo

    do j=0, 3                     !TEST 2
       mdof = Gdof
       if(j == 0) mdof = Qdof

       psihat(1, j, 1:mdof) = x(j, 1:mdof, 1)   !x^
       psihat(2, j, 1:mdof) = x(j, 1:mdof, 2)   !y^

       do l=1, mdof
          psi(1:2, j, l) = MATMUL(B, psihat(1:2, j, l))
          !psi(2, j, 1:mdof) =
       enddo
    enddo                          !TEST 2


    do j=0, 3
       mdof = Gdof
       if(j == 0) mdof = Qdof

       it = 0  !poradi bazove funkce       !UPRAVA
       do k=0, deg                         !  |
          do i=0, k                        !  V
             it = it + 1
             RTN%phi(it, 1, j, 1:mdof) = x(j, 1:mdof, 1)**(k-i) * x(j, 1:mdof, 2)**i
             RTN%phi(it, 2, j, 1:mdof) = 0.
             !print*,'RTN1:',deg,j,it

             it = it + 1
             RTN%phi(it, 1, j, 1:mdof) = 0.
             RTN%phi(it, 2, j, 1:mdof) = x(j, 1:mdof, 1)**(k-i) * x(j, 1:mdof, 2)**i
             !print*,'RTN2:',deg,j,it
          enddo
       enddo

       do k=0, deg
          it = it + 1
          RTN%phi(it, 1, j, 1:mdof) = x(j, 1:mdof, 1)**(deg-k+1) * x(j, 1:mdof, 2)**k   ! x_1*(x_1**(deg-k))*(x_2**k)
          RTN%phi(it, 2, j, 1:mdof) = x(j, 1:mdof, 1)**(deg-k) * x(j, 1:mdof, 2)**(k+1)     ! x_2*(x_1**(deg-k))*(x_2**k)
          !print*,'RTN3:',deg,j,it

       enddo      ! KONEC UPRAV

    enddo


    do j=0, 3
       mdof = Gdof
       if(j == 0) mdof = Qdof

       it = 0
       do k=0, deg
          do i=0, k
             it = it + 1
             RTNtest%phi(1, 1, j, 1:mdof) = RTNtest%phi(1, 1, j, 1:mdof) + it*RTN%phi(it, 1, j, 1:mdof)

             it = it + 1
             RTNtest%phi(1, 2, j, 1:mdof) = RTNtest%phi(1, 2, j, 1:mdof) + it*RTN%phi(it, 2, j, 1:mdof)
          enddo
       enddo

       do k=0, deg
          it = it + 1
          RTNtest%phi(1, 1, j, 1:mdof) = RTNtest%phi(1, 1, j, 1:mdof) + it*RTN%phi(it, 1, j, 1:mdof)
          RTNtest%phi(1, 2, j, 1:mdof) = RTNtest%phi(1, 2, j, 1:mdof) + it*RTN%phi(it, 2, j, 1:mdof)
       enddo

    enddo

    !   print*,'###jake je j:',j

    !   do it =2, dof
    !      RTNtest%phi(it, 1, j, 1:mdof) = 0.
    !      RTNtest%phi(it, 2, j, 1:mdof) = 0.
    !   enddo

       !do j=0, 3
          !mdof = Gdof
          !if(j == 0) mdof = Qdof

          !write(*,'(a3,12es12.4)') '??x', RTNtest%phi(1,1, j, 1:mdof)   !TEST
          !write(*,'(a3,12es12.4)') '??y', RTNtest%phi(1,2, j, 1:mdof)   !TEST

          !do l=1,dof
             !if(j == 0) write(*,'(a3,12es12.4)') '??x', RTN%phi(l,1, j, 1:mdof)
             !if(j == 0) write(*,'(a3,12es12.4)') '??y', RTN%phi(l,2, j, 1:mdof)
             !if(j == 0) print*
          !enddo
       !enddo
    !else
    !   print*,'RTN not yet implemeted '
    !   stop
    !endif

    ! computation of derivatives of RTN basis functions in V_rule integ nodes
    allocate (RTN%Dphi(1:dof, 1:2, 0:3, 1:max(Qdof, Gdof), 1:2 ) )
                                                          !^^^-derivative with respect to the first/second variable
    it = 0  !order of basis functions
    do k=0, deg
       do i=0, k
          it = it + 1
          RTN%Dphi(it, 1, 0, 1:Qdof, 1) = (k-i) * x(0, 1:Qdof, 1)**(k-i-1) * x(0, 1:Qdof, 2)**i
          RTN%Dphi(it, 2, 0, 1:Qdof, 1) = 0.

          RTN%Dphi(it, 1, 0, 1:Qdof, 2) = x(0, 1:Qdof, 1)**(k-i) * i * x(0, 1:Qdof, 2)**(i-1)
          RTN%Dphi(it, 2, 0, 1:Qdof, 2) = 0.

          it = it + 1
          RTN%Dphi(it, 1, 0, 1:Qdof, 1) = 0.
          RTN%Dphi(it, 2, 0, 1:Qdof, 1) = (k-i) * x(0, 1:Qdof, 1)**(k-i-1) * x(0, 1:Qdof, 2)**i

          RTN%Dphi(it, 1, 0, 1:Qdof, 2) = 0.
          RTN%Dphi(it, 2, 0, 1:Qdof, 2) = x(0, 1:Qdof, 1)**(k-i) * i * x(0, 1:Qdof, 2)**(i-1)
       enddo
    enddo

    do k=0, deg
       it = it + 1
       RTN%Dphi(it, 1, 0, 1:Qdof, 1) = (deg-k+1) * x(0, 1:Qdof, 1)**(deg-k) * x(0, 1:Qdof, 2)**k
       RTN%Dphi(it, 2, 0, 1:Qdof, 1) = (deg-k) * x(0, 1:Qdof, 1)**(deg-k-1) * x(0, 1:Qdof, 2)**(k+1)

       RTN%Dphi(it, 1, 0, 1:Qdof, 2) = x(0, 1:Qdof, 1)**(deg-k+1) * k * x(0, 1:Qdof, 2)**(k-1)
       RTN%Dphi(it, 2, 0, 1:Qdof, 2) = x(0, 1:Qdof, 1)**(deg-k) * (k+1) * x(0, 1:Qdof, 2)**k
    enddo

          !do l=1,dof
             !write(*,'(a7,12es12.4)') '??x, dx', RTN%Dphi(l,1, 0, 1:Qdof, 1)
             !write(*,'(a7,12es12.4)') '??y, dx', RTN%Dphi(l,2, 0, 1:Qdof, 1)

             !write(*,'(a7,12es12.4)') '??x, dy', RTN%Dphi(l,1, 0, 1:Qdof, 2)
             !write(*,'(a7,12es12.4)') '??y, dy', RTN%Dphi(l,2, 0, 1:Qdof, 2)
             !print*
          !enddo
          !stop


    ! computation of RTN momentums
    allocate( RTN%MM(1:dof, 1:dof) )

    allocate( RTNtest%MM(1:dof, 1:dof) )    !TEST
    RTNtest%MM(:, :) = 0.

    call Compute_RTN_Momentums(V_rule, G_rule, RTN%deg, dof, Qdof, Gdof,  &
         RTN%phi(1:dof, 1:2, 0:3, 1:max(Qdof, Gdof) ), RTN%MM(1:dof, 1:dof) )

    call Compute_RTN_Momentums(V_rule, G_rule, RTN%deg, dof, Qdof, Gdof,  &
         RTNtest%phi(1:dof, 1:2, 0:3, 1:max(Qdof, Gdof) ), RTNtest%MM(1:dof, 1:dof) )  !TEST

    !do l=1,dof
    !   !write(*,'(10es12.4)') RTN%MM(l, 1:dof)
    !enddo

    allocate( RTN%MMinv(1:RTN%dof, 1:RTN%dof) )

    RTN%MMinv(1:RTN%dof, 1:RTN%dof) = RTN%MM(1:RTN%dof, 1:RTN%dof)
    call MblockInverse(RTN%dof, RTN%MMinv(1:RTN%dof, 1:RTN%dof) )

    allocate( thetai(1:dof) )

    thetai(1:dof) = MATMUL(RTN%MMinv, RTNtest%MM(1:dof, 1))

    !write(*,'(10es12.4)') thetai(1:dof)

    deallocate(x, thetai, RTNtest%phi, RTNtest%MM,  psihat)
    !stop

  end subroutine Init_RTN_FE



  !> setting of RTN elements in integ nodes for a given V_rule
  subroutine Init_RTN_FE_Vnode(RTN, V_rule)
    type(rtn_fe), intent(inout) :: RTN
    type(volume_rule), intent(in) :: V_rule
    real, dimension(:,:), allocatable :: x
    !real, dimension(:), allocatable :: thetai     !TEST
    integer :: dof, Qdof, deg, Qdeg
    integer :: i, it, k

    Qdeg = V_rule%Qdeg

    Qdof = V_rule%Qdof

    RTN%Vnode(Qdeg)%def = .true.


    allocate(x(1:Qdof, 1:2))
    x(1:Qdof, 1:2) = V_rule%lambda(1:Qdof,1:2)

    deg = RTN%deg
    dof = (deg+1)*(deg+3)   !pro N=2

    allocate (RTN%Vnode(Qdeg)%phi(1:dof, 1:2, 1:Qdof ) )

    RTN%Vnode(Qdeg)%phi(1:dof, 1:2, 1:Qdof ) = 0.    !initialization

    it = 0  !poradi bazove funkce
    do k=0, deg
       do i=0, k
          it = it + 1
          RTN%Vnode(Qdeg)%phi(it, 1, 1:Qdof) = x(1:Qdof, 1)**(k-i) * x(1:Qdof, 2)**i
          RTN%Vnode(Qdeg)%phi(it, 2, 1:Qdof) = 0.

          it = it + 1
          RTN%Vnode(Qdeg)%phi(it, 1, 1:Qdof) = 0.
          RTN%Vnode(Qdeg)%phi(it, 2, 1:Qdof) = x(1:Qdof, 1)**(k-i) * x(1:Qdof, 2)**i
       enddo
    enddo

    do k=0, deg
       it = it + 1
       RTN%Vnode(Qdeg)%phi(it, 1, 1:Qdof) = x(1:Qdof, 1)**(deg-k+1) * x(1:Qdof, 2)**k
       RTN%Vnode(Qdeg)%phi(it, 2, 1:Qdof) = x(1:Qdof, 1)**(deg-k)   * x(1:Qdof, 2)**(k+1)
    enddo

    deallocate(x)

  end subroutine Init_RTN_FE_Vnode


  !> setting of derivatives of RTN basis functions in integ nodes for a given V_rule
  subroutine Init_RTN_FE_Vnode_Der(RTN, V_rule)
    type(rtn_fe), intent(inout) :: RTN
    type(volume_rule), intent(in) :: V_rule
    real, dimension(:,:), allocatable :: x
    integer :: dof, Qdof, deg, Qdeg, it, k, i

    Qdeg = V_rule%Qdeg
    Qdof = V_rule%Qdof

    RTN%Vnode(Qdeg)%defder = .true.

    allocate(x(1:Qdof, 1:2))
    x(1:Qdof, 1:2) = V_rule%lambda(1:Qdof,1:2)

    deg = RTN%deg
    dof = (deg+1)*(deg+3)   !pro N=2

    allocate (RTN%Vnode(Qdeg)%Dphi(1:dof, 1:2, 1:Qdof, 1:2 ) )

    RTN%Vnode(Qdeg)%Dphi(1:dof, 1:2, 1:Qdof, 1:2 ) = 0.    !initialization

    it = 0  !order of basis functions
    do k=0, deg
       do i=0, k
          it = it + 1
          if(k-i /= 0) &
               RTN%Vnode(Qdeg)%Dphi(it, 1, 1:Qdof, 1) = (k-i) * x(1:Qdof, 1)**(k-i-1) * x(1:Qdof, 2)**i
          RTN%Vnode(Qdeg)%Dphi(it, 2, 1:Qdof, 1) = 0.

          if(i /= 0) &
               RTN%Vnode(Qdeg)%Dphi(it, 1, 1:Qdof, 2) = x(1:Qdof, 1)**(k-i) * i * x(1:Qdof, 2)**(i-1)
          RTN%Vnode(Qdeg)%Dphi(it, 2, 1:Qdof, 2) = 0.

          it = it + 1
          RTN%Vnode(Qdeg)%Dphi(it, 1, 1:Qdof, 1) = 0.
          RTN%Vnode(Qdeg)%Dphi(it, 2, 1:Qdof, 1) = (k-i) * x(1:Qdof, 1)**(k-i-1) * x(1:Qdof, 2)**i

          RTN%Vnode(Qdeg)%Dphi(it, 1, 1:Qdof, 2) = 0.
          RTN%Vnode(Qdeg)%Dphi(it, 2, 1:Qdof, 2) = x(1:Qdof, 1)**(k-i) * i * x(1:Qdof, 2)**(i-1)
       enddo
    enddo

    do k=0, deg
       it = it + 1
       RTN%Vnode(Qdeg)%Dphi(it, 1, 1:Qdof, 1) = (deg-k+1) * x(1:Qdof, 1)**(deg-k) * x(1:Qdof, 2)**k
       RTN%Vnode(Qdeg)%Dphi(it, 2, 1:Qdof, 1) = (deg-k) * x(1:Qdof, 1)**(deg-k-1) * x(1:Qdof, 2)**(k+1)

       RTN%Vnode(Qdeg)%Dphi(it, 1, 1:Qdof, 2) = x(1:Qdof, 1)**(deg-k+1) * k * x(1:Qdof, 2)**(k-1)
       RTN%Vnode(Qdeg)%Dphi(it, 2, 1:Qdof, 2) = x(1:Qdof, 1)**(deg-k) * (k+1) * x(1:Qdof, 2)**k
    enddo

    !do i=1,dof
    !write(*,'(a7,12es12.4)') '??x, dx', RTN%Vnode(Qdeg)%Dphi(i, 1, 1:Qdof, 1)
    !write(*,'(a7,12es12.4)') '??y, dx', RTN%Vnode(Qdeg)%Dphi(i, 2, 1:Qdof, 1)

    !write(*,'(a7,12es12.4)') '??x, dy', RTN%Vnode(Qdeg)%Dphi(i, 1, 1:Qdof, 2)
    !write(*,'(a7,12es12.4)') '??y, dy', RTN%Vnode(Qdeg)%Dphi(i, 2, 1:Qdof, 2)
    !print*
    !enddo


    deallocate(x)
    !stop
  end subroutine Init_RTN_FE_Vnode_Der


  !> setting the volume and edge integ. nodes
  subroutine Set_RTN_integ_node(V_rule, G_rule, Qdof, Gdof, xi)
    type(volume_rule), intent(in) :: V_rule
    type(Gauss_rule), intent(in) :: G_rule
    integer, intent(in) :: Qdof, Gdof
    real, dimension(0:3,1:max(Qdof,Gdof),1:nbDim), intent(inout) :: xi
    real, dimension(1:nbDim) :: x0, a0
    real :: t
    integer :: ie, l

    ! volume integ. nodes

    !TODO: problem - all of G_rules on edge have to be the same

    if ( Qdof /= V_rule%Qdof .or. Gdof /= G_rule%Qdof ) &
      stop 'Problem in Set_RTN_integ_node - # of nodes is different from Qdof'

    xi(0,1:Qdof, 1:2) = V_rule%lambda(1:Qdof,1:2)

    ! face integ. nodes
    do ie = 1, 3   ! loop through edges
       if(ie == 1 ) then
          x0(1:2) = (/0., 0./)
          a0(1:2) = (/1., 0./)
       elseif(ie == 2 ) then
          x0(1:2) = (/1., 0./)
          a0(1:2) = (/-1., 1./)
       elseif(ie == 3) then
          x0(1:2) = (/0., 1./)
          a0(1:2) = (/0., -1./)
       endif

       do l=1, Gdof  !loop though 1D Gauss quadrature rule
          t = G_rule%lambda(l)
          xi(ie,l,1:2) = x0(1:2) +  a0(1:2) * t
       enddo ! l
    enddo ! ie
  end subroutine Set_RTN_integ_node


  !> computation of moments of RTN elements
  !> deg = degree of RTN element
  !> dof = degree of freedoms = number of momentums
  !> Qdof = degree of integ. nodes
  subroutine Compute_RTN_Momentums(V_rule, G_rule, deg, dof, Qdof, Gdof, phi, MM )
    type(volume_rule), intent(in) :: V_rule
    type(Gauss_rule), intent(in) :: G_rule
    integer, intent(in) :: deg, dof, Qdof, Gdof
    real, dimension(1:dof, 1:2, 0:3, 1:max(Qdof, Gdof)), intent(in):: phi
    real, dimension(1:dof, 1:dof) :: MM
    real, dimension(:,:), allocatable :: qi
    real, dimension(:,:,:), allocatable :: x
    real, dimension(1:2) :: nn
    integer :: i,j, ib, it, deg_interior

    allocate( qi(1:max(Qdof, Gdof), 1:2  ) )

  ! if(deg == 1) then                   !UPRAVA
       ! P_1 RTN functions
    ib = 0   ! index of momentum

    ! face momentums
    do i=1,3 ! loop over triagle edges
       ! outer normal
       if(i == 1 ) nn(1:2) = (/ 0., -1./)
       if(i == 2 ) nn(1:2) = (/ 1.,  1./)
       if(i == 3 ) nn(1:2) = (/-1.,  0./)

       do j=1, deg+1  ! deg+1 degrees of freedom over face   !UPRAVA
          ib = ib + 1

          call EvalMomentumFace(i, deg, ib, Gdof, G_rule, qi(1:Gdof, 1))

          !if (ib == 6) then
          !   write(*,'(a5,i1,a8,2F8.2)') 'hrana', i, 'nn(1:2)', nn(1:2)  !TEST
          !   write(*,'(a16,2es12.4)') 'G_rule%weights', G_rule%weights(1:Gdof)  !TEST
          !   write(*,'(a12,2es12.4)') 'qi(1:Gdof)', qi(1:Gdof,1)  !TEST
          !endif


          do it = 1, dof
             MM(ib, it) = dot_product(G_rule%weights(1:Gdof) * qi(1:Gdof,1),  &
                  phi(it, 1, i,  1:Gdof) * nn(1) &
                  + phi(it, 2, i,  1:Gdof) * nn(2) )
              !  if (ib == 6) then
              !    write(*,'(a2,i2,a3,12es12.4)') 'it', it, '??x', phi(it,1, i, 1:Gdof)  !TEST
               !   write(*,'(a2,i2,a3,12es12.4)') 'it', it, '??y', phi(it,2, i, 1:Gdof)  !TEST
                !endif

          enddo

          !if (ib == 6) then
          !   write(*,'(a12,es12.4)') 'MM(ib,5)=', MM(ib,5)  !TEST
          !endif

       enddo
    enddo

    allocate(x(0:3, 1:max(Qdof, Gdof), 1:2))                !change
    !--------- ^^^  0= volume, 1,2,3= 1st, 2nd, 3rd face

    call Set_RTN_integ_node(V_rule, G_rule, Qdof, Gdof, x)

    deg_interior = deg*(deg+1)  ! deg*(deg+1) degrees of freedom over element   !UPRAVA

    ! volumes momentums
    do j=1, deg_interior
       ib = ib + 1

       call EvalMomentumVolume(deg, ib, Qdof, 0, x(0:3, 1:Qdof, 1:2), qi(1:Qdof, 1:2))
       do it = 1, dof

          MM(ib, it) = dot_product( 0.5* V_rule%weights(1:Qdof),   &
               phi(it, 1, 0, 1:Qdof) * qi(1:Qdof, 1) &
               + phi(it, 2, 0, 1:Qdof) * qi(1:Qdof, 2) )

       enddo

    enddo

   ! else                               !UPRAVA
   !    print*,'deg > 1 in Compute_RTN_Momentums not yet implemented'
   !    stop
   ! endif

    deallocate(qi, x)
  end subroutine Compute_RTN_Momentums


  !> evaluation of the test function for the "edge momentums"
  subroutine EvalMomentumFace(ie, deg, ib, Qdof, G_rule,  qi)
    type(Gauss_rule), intent(in) :: G_rule
    integer, intent(in) :: deg, ib, Qdof, ie
    real, dimension(1:Qdof), intent(inout) :: qi
    integer :: k

    ! better approach ???
    k = mod(ib, deg+1)  ! i.e., k == ib .or. ib == (deg+1)+k .or. ib == 2*(deg+1)+k
    if( k==0) k = deg+1
    qi(1:Qdof) = (G_rule%lambda(1:Qdof))**(k-1)

    return

    ! original version:
    do k=1, deg+1
       if(ib == k .or. ib == (deg+1)+k .or. ib == 2*(deg+1)+k) then
          !print*,'$$$$$',k,ib,ie
          if(ie == 1 ) then
             qi(1:Qdof) = (G_rule%lambda(1:Qdof))**(k-1)
          elseif(ie == 2 ) then
             qi(1:Qdof) = (G_rule%lambda(1:Qdof))**(k-1)
             !qi(1:Qdof) = (1. - G_rule%lambda(1:Qdof))**(k-1)  ! new test

          elseif(ie == 3) then
             if (k == 1) then
                qi(1:Qdof) = 1.
             else
                !qi(1:Qdof) = (G_rule%lambda(1:Qdof))**(k-1) ! new test
                qi(1:Qdof) = (1. - G_rule%lambda(1:Qdof))**(k-1)
             endif
          endif

          goto 101
       endif
    enddo

   ! else                   !UPRAVA
   !    print*,'deg > 1 in EvalMomentumFace not yet implemented'
   !    stop
   ! endif
101 continue
end subroutine EvalMomentumFace

  !> evaluation of the test function for the "volume momentums"
  subroutine EvalMomentumVolume(deg, ib, Qdof, j, x,  qi)
    !type(volume_rule), intent(in) :: V_rule
    integer, intent(in) :: deg, ib, Qdof, j
    real, dimension(0:3,1:Qdof,1:2), intent(inout) :: x
    real, dimension(1:Qdof, 1:2 ) :: qi
    integer :: k, i, poc, poradi

    poradi = ib - 3*(deg+1)
    poc = 0
    do k = 0, deg-1
       do i = 0, k
          poc = poc + 1
          if (poradi == poc) then
             qi(1:Qdof, 1) = x(j, 1:Qdof, 1)**(k-i) * x(j, 1:Qdof, 2)**i !V_rule%lambda(1:Qdof,1)**(k-i) * V_rule%lambda(1:Qdof,2)**i
             qi(1:Qdof, 2) = 0.
             goto 102
          endif
          poc = poc + 1
          if (poradi == poc) then
             qi(1:Qdof, 1) = 0.
             qi(1:Qdof, 2) = x(j, 1:Qdof, 1)**(k-i) * x(j, 1:Qdof, 2)**i  !V_rule%lambda(1:Qdof,1)**(k-i) * V_rule%lambda(1:Qdof,2)**i
             goto 102
          endif
       enddo
    enddo

102 continue
  end subroutine EvalMomentumVolume


  !> evaluation of the test function for the "volume momentums"
  subroutine EvalMomentumVolumeD(deg, ib, Qdof, x,  qi)
    !type(volume_rule), intent(in) :: V_rule
    integer, intent(in) :: deg, ib, Qdof
    real, dimension(1:Qdof,1:2), intent(inout) :: x
    real, dimension(1:Qdof, 1:2 ) :: qi
    integer :: k, i, poc, poradi
    real :: factor

    !x(1:Qdof, 1:2) = x(1:Qdof, 1:2) -1./3

    poradi = ib - 3*(deg+1)

    if (poradi <= 0) then
      print*, 'problem in EvalMomentumVolumeD: poradi:', poradi
      stop
    endif

    poc = 0
    do k = 0, deg-1
       do i = 0, k
          poc = poc + 1
          if (poradi == poc) then
             qi(1:Qdof, 1) = x(1:Qdof, 1)**(k-i) * x(1:Qdof, 2)**i !V_rule%lambda(1:Qdof,1)**(k-i) * V_rule%lambda(1:Qdof,2)**i
             qi(1:Qdof, 2) = 0.
             goto 102
          endif
          poc = poc + 1
          if (poradi == poc) then
             qi(1:Qdof, 1) = 0.
             qi(1:Qdof, 2) = x(1:Qdof, 1)**(k-i) * x(1:Qdof, 2)**i  !V_rule%lambda(1:Qdof,1)**(k-i) * V_rule%lambda(1:Qdof,2)**i
             goto 102
          endif
       enddo
    enddo

102 continue

    ! normalization w.r.t. polynomial degree p
    factor = factorial(k) / factorial(k-i)/factorial(i)
    factor = factor * (k+1) *(k+2)

    !! VD
    qi(1:Qdof, 1:2) = qi(1:Qdof, 1:2) * factor

    !print*,'###EDE',ib,poradi, k, i, k-i, factor

  end subroutine EvalMomentumVolumeD


end module rav_tho_nedelec
