!> local Raviart-Thomas-Nedelec (RTN) finite elements
!> subroutines using global structures
module eval_rav_tho_ned
  use geometry
  use basis
  use integration
  use lapack_oper
  use loc_rav_tho_ned
  use define_state
  use main_data
  use f_mapping

  implicit none

  public :: Eval_Loc_RTN_Elem
  public :: Eval_Loc_RTN_ElemGE
  public :: Eval_Loc_RTN_Div_Elem
  public :: Eval_Loc_RTN_Edge
  public :: EvalMomentumFaceHG
  public :: Eval_Loc_RTN_Edge_HG

  public :: EvalRTNMomentumEdgeTestFunction
contains


  !> evaluate local basis RTN functions at volumes integration nodes
  subroutine Eval_Loc_RTN_Elem(elem, loc_RTN, psi )
    type(element), intent(in) :: elem
    type(basis_rtn_fe), intent(in) :: loc_RTN
    real, dimension(1:loc_RTN%Fdof, 1:nbDim, 1:elem%Qdof), intent(out) :: psi
    real, dimension(:,:), pointer:: phi        ! DGFE test functions
    real, dimension(:,:), allocatable :: Fx
    integer :: Fdof, Qdof, dof, Qnum
    integer :: i, j, indx, iphi

    Fdof = loc_RTN%Fdof
    !dof = elem%dof
    dof = maxval(loc_RTN%ipsi(:, 2))

    Qnum = elem%Qnum
    Qdof = elem%Qdof

    ! setting of physical coordinates of integration nodes
    allocate(Fx(1:Qdof, 1:nbDim))
    call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof, 1:nbDim), &
         Fx(1:Qdof, 1:nbDim) )

    ! "renormalization" of the coordinates
    do i=1,Qdof
       Fx(i, 1:nbDim) = (Fx(i, 1:nbDim) - elem%xc(1:nbDim))/elem%diam
    enddo

    ! pointer to DGFE functions
    phi => state%space%V_rule(Qnum)%phi(1:dof,1:Qdof)

    psi(:,:,:) = 0.

    do i=1,Fdof
       indx  = loc_RTN%ipsi(i, 1)
       iphi = loc_RTN%ipsi(i, 2)

       if(indx <= nbDim) then
          psi(i, indx, 1:Qdof) = phi(iphi, 1:Qdof)

       else
          psi(i, 1, 1:Qdof) = phi(iphi, 1:Qdof) * Fx(1:Qdof,1)
          psi(i, 2, 1:Qdof) = phi(iphi, 1:Qdof) * Fx(1:Qdof,2)

       endif

    enddo

    ! only for output
    !call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof, 1:nbDim),&
    !     Fx(1:Qdof, 1:nbDim) )
    !do j=1,Fdof
    !   do i=1,Qdof
    !      write(100*elem%i+50+j,*) Fx(i,1:nbDim), psi(j,1:2,i)
    !   enddo
    !enddo


    deallocate(Fx)
  end subroutine Eval_Loc_RTN_Elem


  !> evaluate local basis RTN functions
  !> (2nd idx =1,2) and its divergence (2nd idx = 3)
  !> at GEneral volumes integration nodes
  !> of quadrature V_rule%Qnum != elem%Qnum
  subroutine Eval_Loc_RTN_ElemGE(elem, V_rule, loc_RTN, psi )
    type(element), intent(in) :: elem
    type(volume_rule), target, intent(in) :: V_rule
    type(basis_rtn_fe), intent(in) :: loc_RTN
    real, dimension(1:loc_RTN%Fdof, 1:3, 1:V_rule%Qdof) :: psi
    real, dimension(:,:), pointer:: phi        ! DGFE test functions
    real, dimension(:,:,:), allocatable :: Dphi
    real, dimension(:,:), allocatable :: Fx
    integer :: Fdof, Qdof, dof
    integer :: i, j, indx, iphi

    Fdof = loc_RTN%Fdof
    Qdof = V_rule%Qdof

    ! setting of physical coordinates of integration nodes
    allocate(Fx(1:Qdof, 1:nbDim))
    call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof, 1:nbDim), &
         Fx(1:Qdof, 1:nbDim) )

    ! "renormalization" of the coordinates
    do i=1,Qdof
       Fx(i, 1:nbDim) = (Fx(i, 1:nbDim) - elem%xc(1:nbDim))/elem%diam
    enddo

    ! pointer to DGFE functions
    !dof = DOFtriang( MaxDegreeImplemented )
    dof = maxval(loc_RTN%ipsi(:, 2))
    phi => V_rule%phi(1:dof,1:Qdof)

    ! their derivatives of  DGFE functions
    allocate( Dphi(1:dof, 1:nbDim, 1:Qdof) )
    call Eval_Dphi_plus(elem, V_rule, dof, Dphi(1:dof, 1:nbDim, 1:Qdof) )

    psi(:,:,:) = 0.

    do i=1,Fdof
       indx  = loc_RTN%ipsi(i, 1)
       iphi = loc_RTN%ipsi(i, 2)

       if(indx <= nbDim) then
          psi(i, indx, 1:Qdof) =  phi(iphi, 1:Qdof)     ! 1st or 2nd comp of \Psi
          psi(i,   3,  1:Qdof) = Dphi(iphi, indx, 1:Qdof)  ! div of \Psi

       else
          psi(i, 1, 1:Qdof) = phi(iphi, 1:Qdof) * Fx(1:Qdof,1)
          psi(i, 2, 1:Qdof) = phi(iphi, 1:Qdof) * Fx(1:Qdof,2)

          ! chain rule  (x/h \phi(x))' = D\phi(x) + x/h phi(x)
          psi(i, 3, 1:Qdof) = &
               Dphi(iphi, 1, 1:Qdof) * Fx(1:Qdof,1) + phi(iphi,1:Qdof)/elem%diam &
               +Dphi(iphi, 2, 1:Qdof) * Fx(1:Qdof,2) + phi(iphi,1:Qdof)/elem%diam
       endif

    enddo

    ! ! graphical verification
     ! call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof, 1:nbDim),&
     !      Fx(1:Qdof, 1:nbDim) )
     ! do j=1,Fdof
     !    do i=1,Qdof
     !       write(100*elem%i+50+j,*) Fx(i,1:nbDim), psi(j,1:3,i)
     !    enddo
     ! enddo
     ! stop


    deallocate(Fx, Dphi)
  end subroutine Eval_Loc_RTN_ElemGE


  !> evaluate divergence local basis RTN functions at volumes integration nodes
  subroutine Eval_Loc_RTN_Div_Elem(elem, loc_RTN, Divpsi )
    type(element), intent(in) :: elem
    type(basis_rtn_fe), intent(in) :: loc_RTN
    real, dimension(1:loc_RTN%Fdof, 1:elem%Qdof) :: Divpsi
    real, dimension(:,:),  pointer :: phi        ! DGFE test functions
    real, dimension(:,:,:), allocatable:: Dphi        ! DGFE test functions
    real, dimension(:,:), allocatable :: Fx
    integer :: Fdof, Qdof, dof, Qnum
    integer :: i, j, indx, iphi

    Fdof = loc_RTN%Fdof
    dof = elem%dof
    Qnum = elem%Qnum
    Qdof = elem%Qdof

    ! setting of physical coordinates of integration nodes
    allocate(Fx(1:Qdof, 1:nbDim))
    call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof, 1:nbDim), &
         Fx(1:Qdof, 1:nbDim) )

    ! "renormalization" of the coordinates
    do i=1,Qdof
       Fx(i, 1:nbDim) = (Fx(i, 1:nbDim) - elem%xc(1:nbDim))/elem%diam
    enddo


    ! pointer to DGFE functions
    phi => state%space%V_rule(Qnum)%phi(1:dof,1:Qdof)

    ! their derivatives of  DGFE functions
    allocate( Dphi(1:dof, 1:nbDim, 1:Qdof) )
    call Eval_Dphi(elem, dof, Dphi(1:dof, 1:nbDim, 1:Qdof) )


    do i=1,Fdof
       indx  = loc_RTN%ipsi(i, 1)
       iphi = loc_RTN%ipsi(i, 2)

       if(indx <= nbDim) then
          Divpsi(i, 1:Qdof) = Dphi(iphi, indx, 1:Qdof)

       else ! chain rule  (x/h \phi(x))' = D\phi(x) + x/h phi(x)
          Divpsi(i, 1:Qdof) = &
               Dphi(iphi, 1, 1:Qdof) * Fx(1:Qdof,1) + phi(iphi,1:Qdof)/elem%diam  + &
               Dphi(iphi, 2, 1:Qdof) * Fx(1:Qdof,2) + phi(iphi,1:Qdof)/elem%diam
       endif

    enddo

    deallocate(Dphi)

    ! only for output
    !call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof, 1:nbDim),&
    !     Fx(1:Qdof, 1:nbDim) )
    !do j=1,Fdof
    !   do i=1,Qdof
    !      write(300*elem%i+50+j,*) Fx(i,1:nbDim), Divpsi(j,i)
    !   enddo
    !enddo


    deallocate(Fx)
  end subroutine Eval_Loc_RTN_Div_Elem


  !> evaluate local basis RTN functions at edge  integration nodes
  subroutine Eval_Loc_RTN_Edge(elem, loc_RTN, ie, psi )
    type(element), intent(in) :: elem
    type(basis_rtn_fe), intent(in) :: loc_RTN
    integer, intent(in) :: ie  ! index of the edge
    real, dimension(1:loc_RTN%Fdof, 1:nbDim, 1:elem%face(fGdof, ie)), intent(out) :: psi
    real, dimension(:,:), allocatable :: phi        ! DGFE test functions
    real, dimension(:,:), allocatable :: x, Fx
    integer :: Fdof, Qdof, dof, Qnum
    integer :: i, j, indx, iphi

    Fdof = loc_RTN%Fdof
    dof = elem%dof
    Qnum = elem%face(fGnum, ie)
    Qdof = elem%face(fGdof, ie)

    !if (elem%i == 6) then
    !    print*, 'v------------In Eval_Loc_RTN_Edge sub.----------------v'
    !    write(*,'(a29,5i3)') 'ie, Fdof, dof, Gdof, Qnum:', ie, Fdof, dof, Qdof, Qnum
    !    write(*,'(a39,2i3)') 'elem%HGface(1, ie), elem%HGface(2, ie):', elem%HGface(1, ie), elem%HGface(2, ie)
    !    print*, ' '
    !endif

    ! setting of physical coordinates of integration nodes
    allocate(x(1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim))
    call SetEdgeNodes(Qnum, ie, 3, x(1:Qdof, 1:nbDim) )
    call ComputeF(elem, Qdof, x(1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim) )

    !if (elem%i == 6) then
    !    print*, 'v------------In Eval_Loc_RTN_Edge sub.----------------v'
    !    write(*,'(a15,8es14.6)') 'x(1:Qdof, 1):', x(1:Qdof, 1)
    !    write(*,'(a15,8es14.6)') 'x(1:Qdof, 2):', x(1:Qdof, 2)
    !    print*, ' '
    !    write(*,'(a15,8es14.6)') 'Fx(1:Qdof, 1):', Fx(1:Qdof, 1)
    !    write(*,'(a15,8es14.6)') 'Fx(1:Qdof, 2):', Fx(1:Qdof, 2)
    !   print*, ' '
    !endif

    ! "renormalization" of the coordinates
    do i=1,Qdof
       !!write(11,*)x(i,1:2), Fx(i,1:2),ie
       Fx(i, 1:nbDim) = (Fx(i, 1:nbDim) - elem%xc(1:nbDim))/elem%diam
    enddo

    !if (elem%i == 6) then
    !    print*, 'v------------In Eval_Loc_RTN_Edge sub.----------------v'
    !    write(*,'(a15,8es14.6)') 'Fx(1:Qdof, 1) renorm:', Fx(1:Qdof, 1)
    !    write(*,'(a15,8es14.6)') 'Fx(1:Qdof, 2) renorm:', Fx(1:Qdof, 2)
    !   print*, ' '
    !endif

    ! DGFE functions on the edge
    allocate(phi(1:dof, 1:Qdof) )
    call Eval_Phi_Edge(elem, dof, ie, phi, .false.)     ! it is wrong for elements containing HG nodes
    !phi = state%space%G_rule(Qnum)%phi(elem%type, ie, 1, 1:dof, 1:Qdof)

    !if (elem%i == 6) then
    !    print*, 'v------------In Eval_Loc_RTN_Edge sub.----------------v'
    !    do j=1, dof
    !    write(*,'(a5,i3,a10, 8es14.6)') 'phi(', j, ', 1:Qdof):', phi(j, 1:Qdof)
    !    enddo
    !    print*, ' '
    !   pause
    !endif


    psi(:,:,:) = 0.

    do i=1,Fdof
       indx  = loc_RTN%ipsi(i, 1)
       iphi = loc_RTN%ipsi(i, 2)

       if(indx <= nbDim) then
          psi(i, indx, 1:Qdof) = phi(iphi, 1:Qdof)

       else
          psi(i, 1, 1:Qdof) = phi(iphi, 1:Qdof) * Fx(1:Qdof,1)
          psi(i, 2, 1:Qdof) = phi(iphi, 1:Qdof) * Fx(1:Qdof,2)

       endif

    enddo

    ! only for output
    !call ComputeF(elem, Qdof, x(1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim) )
    !do j=1,Fdof
    !   do i=1,Qdof
    !      write(100*elem%i+50+j,*) Fx(i,1:nbDim), psi(j,1:2,i)
    !   enddo
    !enddo


    deallocate(Fx)
  end subroutine Eval_Loc_RTN_Edge


  !> evaluate local basis RTN functions at edge integration nodes including HG nodes
  subroutine Eval_Loc_RTN_Edge_HG(elem, loc_RTN, ie, HGidx, psi, opposite )
    type(element), intent(in) :: elem
    type(basis_rtn_fe), intent(in) :: loc_RTN
    integer, intent(in) :: ie  ! index of the (sub)edge
    integer, intent(in) :: HGidx   ! type of the (sub)edge
    logical, intent(in) :: opposite   ! opposite orientation of integ. nodes on (sub)edge
    real, dimension(1:loc_RTN%Fdof, 1:nbDim, 1:elem%face(fGdof, ie)) :: psi
    real, dimension(:,:), allocatable :: phi        ! DGFE test functions
    real, dimension(:,:), allocatable :: x, Fx
    integer :: Fdof, Qdof, dof, Qnum
    integer :: i, j, indx, iphi

    Fdof = loc_RTN%Fdof
    !dof = elem%dof
    dof = maxval(loc_RTN%ipsi(:, 2))

    Qnum = elem%face(fGnum, ie)
    Qdof = elem%face(fGdof, ie)

    ! setting of physical coordinates of integration nodes
    allocate(x(1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim))
    !write(*,'(a15, i3, a3, i3)') 'elem%HGface(1,', ie,')', elem%HGface(1, ie)
    !write(*,'(a6, i3)') 'HGidx', HGidx
    !write(*,'(a5, i3)') 'Qnum', Qnum
    !write(*,'(a5, i3)') 'Qdof', Qdof
    !write(*,'(a18, 6es14.6)') 'x(1:Qdof, 1)', x(1:Qdof, 1)
    !write(*,'(a18, 6es14.6)') 'x(1:Qdof, 2)', x(1:Qdof, 2)
    if(.not. elem%HGnode) then
       call SetEdgeNodes(Qnum, ie, 3, x(1:Qdof, 1:nbDim) )
    else
       call SetEdgeNodesHG(HGidx, Qnum, elem%HGface(1, ie), 3, x(1:Qdof, 1:nbDim) )
    endif

    !if (elem%i == 6) then
        !write(*,*) 'ie, x(1:Qdof, 1)', ie, x(1:Qdof, 1)
        !write(*,*) 'ie, x(1:Qdof, 2)', ie, x(1:Qdof, 2)
        !print*, ' '
    !endif

    call ComputeF(elem, Qdof, x(1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim) )

    ! "renormalization" of the coordinates
    do i=1,Qdof
       !!write(11,*)x(i,1:2), Fx(i,1:2),ie
       Fx(i, 1:nbDim) = (Fx(i, 1:nbDim) - elem%xc(1:nbDim))/elem%diam
    enddo

    ! DGFE functions on the edge
    allocate(phi(1:dof, 1:Qdof) )

    ! if opposite then reorder
    if(opposite) then
       call Eval_Phi_Edge(elem, dof, ie, phi, .true.)
    else
       call Eval_Phi_Edge(elem, dof, ie, phi, .false.)
    endif

    !if (elem%i == 6) then
    !  do i=1,dof
    !   write(*,*) 'phi(1:dof, 1:Qdof)', i, phi(i, 1:Qdof)
    !    write(*,*) 'phi(1:dof, 1:Qdof)', i, phi(i, 1:Qdof)
    !  enddo
    !endif

    psi(:,:,:) = 0.

    do i=1,Fdof
       indx  = loc_RTN%ipsi(i, 1)
       iphi = loc_RTN%ipsi(i, 2)

       if(indx <= nbDim) then
          psi(i, indx, 1:Qdof) = phi(iphi, 1:Qdof)

       else
          psi(i, 1, 1:Qdof) = phi(iphi, 1:Qdof) * Fx(1:Qdof,1)
          psi(i, 2, 1:Qdof) = phi(iphi, 1:Qdof) * Fx(1:Qdof,2)

       endif

    enddo

    ! only for output
    !call ComputeF(elem, Qdof, x(1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim) )
    !do j=1,Fdof
    !   do i=1,Qdof
    !      write(100*elem%i+50+j,*) Fx(i,1:nbDim), psi(j,1:2,i)
    !   enddo
    !enddo


    deallocate(Fx)
  end subroutine Eval_Loc_RTN_Edge_HG


  !> evaluation of the local RTN functions for the "edge momentums" including HG nodes
  subroutine EvalMomentumFaceHG(ie, HGidx, deg, ib, Qdof, G_rule,  qi)
    type(Gauss_rule), intent(in) :: G_rule
    integer, intent(in) :: ie      ! index of the edge 1, 2, 3
    integer, intent(in) :: HGidx   ! type of the subedge
    integer, intent(in) :: deg     ! degree of the flux reconstruction
    integer, intent(in) :: ib      ! index of the momentum
    integer, intent(in) :: Qdof    ! number of inte nodes
    real, dimension(1:Qdof), intent(inout) :: qi
    integer :: k,j
    real :: t

    k = mod(ib, deg+1)  ! i.e., k == ib .or. ib == (deg+1)+k .or. ib == 2*(deg+1)+k
    if( k==0) k = deg+1

    if(HGidx == 1) then
       qi(1:Qdof) = (G_rule%lambda(1:Qdof))**(k-1)
    else
       do j=1,Qdof
          qi(j) = ResizeHG( G_rule%lambda(j), HGidx)

       enddo

       !write(*,'(a4,i5,20es12.4)') 'hg G',HGidx, G_rule%lambda(:)
       !write(*,'(a4,i5,20es12.4)') 'hg q',HGidx, qi(:)
       !write(*,'(a4,i5,20es12.4)') 'hg q',HGidx, qi(:)**(k-1)

       qi(1:Qdof) = qi(1:Qdof)**(k-1)

    endif

  end subroutine EvalMomentumFaceHG

   ! eval EDGE test function for RTN iMoment-th Momentum problem
   subroutine EvalRTNMomentumEdgeTestFunction(ie, iMoment, G_rule, qi)
      integer, intent(in) :: ie ! index of the edge in a triangle
      integer, intent(in) :: iMoment
      type(Gauss_rule), intent(in) :: G_rule
      !integer, intent(in) :: Qdof ! number of
      real, dimension(1:G_rule%Qdof), intent(inout) :: qi ! test functions in Grule integ nodes

      qi(1:G_rule%Qdof) = ( G_rule%lambda( 1:G_rule%Qdof ) )**(iMoment-1)

   end subroutine EvalRTNMomentumEdgeTestFunction


end module eval_rav_tho_ned
