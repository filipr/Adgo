!> different evaluations of the solution jumps

module eval_jumps
  use mesh_oper
  use main_data
  use data_mod
  use integration
  use blocks_integ
  use model_oper
  use basis
  use eval_sol


  implicit none

  public:: ElementEdgeJump_time
  public:: ElementEdgeJump
  public:: ElementEdgeJumpProj

  public:: IntegElementEdgeJump
  public:: IntegElementEdgeJump_time
  public:: IntegElementEdgeJumpDer

  public:: ElementEdgeJumpIndicator
  public:: ElementEdgeDerJumpIndicator
  
  public::  ElementEdgeJumpsAllProj
contains
  !> evaluate \f$ [{\bf w_h}] \f$ in integ nodes
  !>
  !> \f$ \Gamma \f$ is an face shared by elements 'elem' and 'elem1'
  subroutine ElementEdgeJump(elem,  ie,  jump)
    type(element), intent(in):: elem ! elem = element
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    real, dimension(1:elem%face(fGdof,ie), 1: ndim), intent(out) :: jump
    class(element), pointer ::   elem1  ! elem1 = neigh element
    real, dimension(:,:), allocatable :: wi     ! w recomputed  in integ nodes
    real, dimension(:,:), allocatable :: wii    ! opposite or boundary state
    real, dimension(:), pointer :: xi
    integer ::  dof1, ie1, Qnum, Qdof
    integer :: l, ii

    !! seting of degree of the Gauss quadrature
    Qnum = elem%face(fGnum,ie)
    Qdof = state%space%G_rule(Qnum)%Qdof
    if(Qdof /= elem%face(fGdof,ie)) print*,'## Trouble in ElementEdgeJump'

    allocate(wi(1:Qdof, 1:ndim ), wii(1:Qdof, 1:ndim ) )
    call Eval_w_Edge(elem, ie, wi(1:Qdof, 1:ndim), .false.)

    !print*, 'ElementEdgeJump wi1:', wi(1:Qdof, 1:ndim)

    ii = elem%face(neigh, ie)
    if( ii > 0) then  !! inner face
       elem1 => grid%elem(ii)
       ie1 = elem%face(nei_i,ie)

       if(Qdof /= elem1%face(fGdof,ie1)) print*,'## Trouble in ElementEdgeJump (2)'

       call Eval_w_Edge(elem1, ie1, wii(1:Qdof, 1:ndim), .true.)

       !print*, 'ElementEdgeJump elem(',ii,') wi1:', wii(1:Qdof, 1:ndim)


    else  !! boundary face

       if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
          do l=1,Qdof
             !xi => grid%b_edge(-elem%face(neigh,ie))%x_div(l, 1:nbDim)
             !call Exact_Scalar(xi(1:nbDim), wii(l,1:ndim), state%time%ctime )


             call Exact_Scalar(grid%b_edge(-elem%face(neigh,ie))%x_div(l, 1:nbDim), &
                  wii(l,1:ndim), state%time%ctime )
          enddo
       else
          ! NOT YET IMPLEMENTED !!!!!!!!!
          wii(1:Qdof, 1:ndim) = wi(1:Qdof, 1:ndim)
       endif
    endif

    jump(1:Qdof,1:ndim) = wi(1:Qdof,1:ndim) - wii(1:Qdof,1:ndim)

    deallocate(wi, wii)

  end subroutine ElementEdgeJump


  !> evaluate \f$ [{\bf w_h}] \f$ in integ nodes
  !>
  !> \f$ \Gamma \f$ is an face shared by elements 'elem' and 'elem1'  at t_{k-1}
  subroutine ElementEdgeJump_time(elem,  ie,  jump)
    type(element), intent(in):: elem ! elem = element
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    real, dimension(1:elem%face(fGdof,ie), 1: ndim), intent(out) :: jump
    class(element), pointer ::   elem1  ! elem1 = neigh element
    real, dimension(:,:), allocatable :: wi     ! w recomputed  in integ nodes
    real, dimension(:,:), allocatable :: wii    ! opposite or boundary state
    real, dimension(:), pointer :: xi
    integer ::  dof1, ie1, Qnum, Qdof
    integer :: l, ii

    !! seting of degree of the Gauss quadrature
    Qnum = elem%face(fGnum,ie)
    Qdof = state%space%G_rule(Qnum)%Qdof
    if(Qdof /= elem%face(fGdof,ie)) print*,'## Trouble in ElementEdgeJump'

    allocate(wi(1:Qdof, 1:ndim ), wii(1:Qdof, 1:ndim ) )
    call Eval_w_Edge_time(elem, ie, wi(1:Qdof, 1:ndim), .false.)

    ii = elem%face(neigh, ie)
    if( ii > 0) then  !! inner face
       elem1 => grid%elem(ii)
       ie1 = elem%face(nei_i,ie)

       if(Qdof /= elem1%face(fGdof,ie1)) print*,'## Trouble in ElementEdgeJump (2)'

       call Eval_w_Edge_time(elem1, ie1, wii(1:Qdof, 1:ndim), .true.)

    else  !! boundary face

       if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
          do l=1,Qdof
             xi => grid%b_edge(-elem%face(neigh,ie))%x_div(l, 1:nbDim)
             call Exact_Scalar(xi(1:nbDim), wii(l,1:ndim), state%time%ttime -state%time%tau(1) )
          enddo
       else
          ! NOT YET IMPLEMENTED !!!!!!!!!
          wii(1:Qdof, 1:ndim) = wi(1:Qdof, 1:ndim)
       endif
    endif

    jump(1:Qdof,1:ndim) = wi(1:Qdof,1:ndim) - wii(1:Qdof,1:ndim)


    deallocate(wi, wii)

  end subroutine ElementEdgeJump_time


  !> evaluate \f$ [{\bf w_h}] \f$ and \f$ [{\bf \Pi^{p-1} w_h}] \f$ in integ nodes
  !>
  !> \f$ \Gamma \f$ is an face shared by elements 'elem' and 'elem1'
  subroutine ElementEdgeJumpProj(elem,  ie,  jump, jump1, grid_L)
    type(element), target, intent(in):: elem ! elem = element
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    real, dimension(1:elem%face(fGdof,ie), 1: ndim), intent(out) :: jump, jump1
    class(mesh), pointer, intent(in), optional :: grid_L
    class(element), pointer ::   elem1  ! elem1 = neigh element
    real, dimension(:,:,:), allocatable :: wi     ! w recomputed  in integ nodes
    real, dimension(:,:,:), allocatable :: wii    ! opposite or boundary state
    real, dimension(:,:), allocatable :: nc       ! outer normals
    real, dimension(:), pointer :: xi
    integer ::  dof1, ie1, Qnum, Qdof
    integer :: l, ii

    !! seting of degree of the Gauss quadrature
    Qnum = elem%face(fGnum,ie)
    Qdof = state%space%G_rule(Qnum)%Qdof
    if(Qdof /= elem%face(fGdof,ie)) print*,'## Trouble in ElementEdgeJump'

    allocate(wi(1:2,1:Qdof, 1:ndim ), wii(1:2,1:Qdof, 1:ndim ) )
    call Eval_w_EdgeProj(elem, ie, wi(1, 1:Qdof, 1:ndim), wi(2, 1:Qdof, 1:ndim), .false.)

    ii = elem%face(neigh, ie)
    if( ii > 0) then  !! inner face
       if(present(grid_L)) then
          elem1 => grid_L%elem(ii)
       else
          elem1 => grid%elem(ii)
       endif

       ie1 = elem%face(nei_i,ie)

       if(Qdof /= elem1%face(fGdof,ie1)) print*,'## Trouble in ElementEdgeJump (2)'

       call Eval_w_EdgeProj(elem1, ie1, wii(1,1:Qdof, 1:ndim), wii(2,1:Qdof, 1:ndim), .true.)

    else  !! boundary face

       if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
          do l=1,Qdof
             if(present(grid_L)) then
                xi => elem%xi(ie, l, 1:nbDim)
             else
                xi => grid%b_edge(-elem%face(neigh,ie))%x_div(l, 1:nbDim)
             endif

             call Exact_Scalar(xi(1:nbDim), wii(1, l,1:ndim), state%time%ctime )
          enddo
          wii(2, 1:Qdof ,1:ndim) = wii(1, 1:Qdof ,1:ndim)

       elseif(state%modelName == 'porous') then
          do l=1,Qdof
             if(present(grid_L)) then
                xi => elem%xi(ie, l, 1:nbDim)
             else
                xi => grid%b_edge(-elem%face(neigh,ie))%x_div(l, 1:nbDim)
             endif

             call Exact_Scalar(xi(1:nbDim), wii(1, l,1:ndim), state%time%ctime )
          enddo
          wii(2, 1:Qdof ,1:ndim) = wii(1, 1:Qdof ,1:ndim)

       elseif( state%type_IC == 8 ) then
          ! IGNORING THE VIOLATION IN BOUNDARY CONDITIONS FOR SHOCK-VORTEX INTERACTION
          wii(1:2,1:Qdof, 1:ndim) = wi(1:2,1:Qdof, 1:ndim)

       else
          ! extrapolation
          wii(1:2,1:Qdof, 1:ndim) = wi(1:2,1:Qdof, 1:ndim)

          if(elem%iBC(ie) == 0 ) then ! SOLID WALLS
             if(state%model%Re == 0.) then        ! Euler equations
                ! setting of outer normals in integration nodes
                allocate(nc(1:Qdof, 1:nbDim) )
                if(elem%ibcur > 0) then  ! UNIT normal
                   nc(1:Qdof,1) = elem%nc(1:Qdof,1)/ elem%dnc(1:Qdof)
                   nc(1:Qdof,2) = elem%nc(1:Qdof,2)/ elem%dnc(1:Qdof)
                else
                   nc(1:Qdof,1) = elem%n(ie,1) / elem%dn(ie)
                   nc(1:Qdof,2) = elem%n(ie,2) / elem%dn(ie)
                endif

                ! mirror state vector
                call Mirror_W(ndim, Qdof, wii(1,1:Qdof, 1:ndim), nc)
                call Mirror_W(ndim, Qdof, wii(2,1:Qdof, 1:ndim), nc)

                ! projection to the tangential component
                !call UpdateMirror(ndim, Qdof, wii(1,1:Qdof, 1:ndim), nc)
                !call UpdateMirror(ndim, Qdof, wii(2,1:Qdof, 1:ndim), nc)

                !if(elem%i > 340 .and. elem%i < 360) then
                !if(elem%i == 507 ) then
                !   print*,elem%xc(:)
                !   do l= 1,Qdof
                !      xi => grid%b_edge(-elem%face(neigh,ie))%x_div(l, 1:nbDim)
                !      write(*,'(20es12.4)') xi(1:2), xi(1:2) + elem%n(ie,1:2),&
                !           xi(1:2) + wi(1, l,2:3)/100., xi(1:2) + wii(1, l,2:3)/100.
                !   enddo
                !   do l= 1,Qdof
                !      write(*,'(2i5,20es12.4)') elem%i,l,elem%n(ie,1:2), &
                !           wi(1, l,2:3) ,wii(1, l,2:3), &
                !           dot_product(wi(1, l,2:3) +wii(1, l,2:3), nc(l,1:2))
                !   enddo
                !endif

                deallocate(nc)

             elseif(state%model%Re > 0.) then     ! Navier-Stokes equations
                wii(1:2,1:Qdof, 2:3) = 0.                      ! no-slip BC
                !wii(1:2,1:Qdof, 2:3) = - wi(1:2,1:Qdof, 2:3)   ! mirror BC
             endif

          endif ! end of SOLID WALLS
       endif
    endif

    jump(1:Qdof,1:ndim) = wi(1, 1:Qdof,1:ndim) - wii(1, 1:Qdof,1:ndim)
    jump1(1:Qdof,1:ndim) = wi(2, 1:Qdof,1:ndim) - wii(2, 1:Qdof,1:ndim)

    deallocate(wi, wii)

  end subroutine ElementEdgeJumpProj

  !> evaluate \f$ \int_{\Gamma} [{\bf w_h}]^2\, dS \f$,
  !>
  !> \f$ \Gamma \f$ is an face shared by elements 'elem' and 'elem1'
  subroutine IntegElementEdgeJump(elem,  ie, jump)
    type(element), intent(in):: elem ! elem = element
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    real, dimension(1:ndim), intent(inout) :: jump
    real, dimension(:,:), allocatable :: wi     ! w recomputed  in integ nodes
    integer ::   Qnum, Qdof

    !! seting of degree of the Gauss quadrature
    Qnum = elem%face(fGnum,ie)
    Qdof = state%space%G_rule(Qnum)%Qdof
    if(Qdof /= elem%face(fGdof,ie)) print*,'## Trouble in ElementEdgeJump'

    allocate(wi(1:Qdof, 1:ndim ) )
    call ElementEdgeJump(elem,  ie, wi)

    wi(1:Qdof,1:ndim) = wi(1:Qdof,1:ndim)**2 * elem%dn(ie)

    jump(1:ndim) =  matmul(state%space%G_rule(Qnum)%weights(1:Qdof), wi(1:Qdof,1:ndim) )

    deallocate(wi)

  end subroutine IntegElementEdgeJump


  !> evaluate only the lifting
  subroutine Eval_Lifting(grid_L)
    class(mesh), pointer, intent(out) :: grid_L
    class(element), pointer :: elem
    real, dimension(:,:), allocatable :: wi,wi1
    real, dimension(:), allocatable :: jump2
    real ::weight
    integer :: i, ie, Qnum, Qdof

    if( state%space%estim_space /= 'pNeu')then
       stop 'Eval_Lifting only for pNeu'
    endif

    do i=1, grid_L%nelem
       elem => grid_L%elem(i)
       elem%lifting(:,:) = 0.

       do ie=1, elem%flen
          Qnum = elem%face(fGnum,ie)
          Qdof = state%space%G_rule(Qnum)%Qdof
          if(Qdof /= elem%face(fGdof,ie)) print*,'## Trouble in ElementEdgeJump'

          allocate(wi(1:Qdof, 1:ndim ), wi1(1:Qdof, 1:ndim ) )
          call ElementEdgeJumpProj(elem,  ie, wi, wi1, grid_L)
          
          allocate(jump2(1:ndim) )
          ! P^0 projection
          jump2(1:ndim) =  matmul(state%space%G_rule(Qnum)%weights(1:Qdof), wi(1:Qdof,1:ndim) )

          ! lifting operator for the dicrete gradient in pNeu for SIPG and NIPG [Ern, Vohralik, SINUM 15]
          weight = 1.0
          if(elem%face(neigh,ie) > 0) weight = 0.5
          elem%lifting(1:ndim,1)=elem%lifting(1:ndim, 1) + weight*elem%n(ie,1) * jump2(1:ndim) / elem%area
          elem%lifting(1:ndim,2)=elem%lifting(1:ndim, 2) + weight*elem%n(ie,2) * jump2(1:ndim) / elem%area
          deallocate(wi, wi1, jump2)
       end do
    enddo

  end subroutine Eval_Lifting


  !> evaluate \f$ \int_{\Gamma} [{\bf w_h}]^2\, dS \f$,
  !> evaluate \f$ \int_{\Gamma} [{\bf \Pi^{p-1} w_h}]^2\, dS \f$ and
  !> evaluate \f$ \int_{\Gamma} (\Pi^0 [w_h])^2\, dS \f$,
  !>
  !> \f$ \Gamma \f$ is an face shared by elements 'elem' and 'elem1'
  subroutine IntegElementEdgeJumpProj(elem,  ie, jump, jump1, jump2)
    type(element), intent(inout):: elem ! elem = element
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    real, dimension(1:ndim), intent(inout) :: jump,jump1, jump2
    real, dimension(:,:), allocatable :: wi, wi1     ! w recomputed  in integ nodes
    real :: weight
    integer ::   Qnum, Qdof

    !print*,'IntegElementEdgeJumpProj(elem,  ie, jump, jump1, jump2)'

    !! seting of degree of the Gauss quadrature
    Qnum = elem%face(fGnum,ie)
    Qdof = state%space%G_rule(Qnum)%Qdof
    if(Qdof /= elem%face(fGdof,ie)) print*,'## Trouble in ElementEdgeJump'

    allocate(wi(1:Qdof, 1:ndim ), wi1(1:Qdof, 1:ndim ) )
    call ElementEdgeJumpProj(elem,  ie, wi, wi1)

    ! P^0 projection
    jump2(1:ndim) =  matmul(state%space%G_rule(Qnum)%weights(1:Qdof), wi(1:Qdof,1:ndim) )

    ! lifting operator for the dicrete gradient in pNeu for SIPG and NIPG [Ern, Vohralik, SINUM 15]
     if( state%space%estim_space == 'pNeu')then
        weight = 1.0
        if(elem%face(neigh,ie) > 0) weight = 0.5
        elem%lifting(1:ndim,1)=elem%lifting(1:ndim, 1) + weight*elem%n(ie,1) * jump2(1:ndim) / elem%area
        elem%lifting(1:ndim,2)=elem%lifting(1:ndim, 2) + weight*elem%n(ie,2) * jump2(1:ndim) / elem%area
     endif
    !print*,'#DE#ESWW##',elem%lifting(1:ndim, 1:2)

    ! L^2-norm of the P^0 projection
    jump2(1:ndim) = jump2(1:ndim)**2 * elem%dn(ie)


    wi(1:Qdof,1:ndim) = wi(1:Qdof,1:ndim)**2 * elem%dn(ie)
    wi1(1:Qdof,1:ndim) = wi1(1:Qdof,1:ndim)**2 * elem%dn(ie)

    jump(1:ndim) =  matmul(state%space%G_rule(Qnum)%weights(1:Qdof), wi(1:Qdof,1:ndim) )
    jump1(1:ndim) =  matmul(state%space%G_rule(Qnum)%weights(1:Qdof), wi1(1:Qdof,1:ndim) )

    deallocate(wi, wi1)

  end subroutine IntegElementEdgeJumpProj



  !> evaluate \f$ \int_{t_{k-1}}^{t_k}\int_{\Gamma} [{\bf w_h}]^2\, dS \, dt\f$,
  !> weight =  \f$ tau_k\f$ is included
  !>
  !> \f$ \Gamma \f$ is an face shared by elements 'elem' and 'elem1'
  subroutine IntegElementEdgeJump_time(elem,  ie, jump)
    type(element), intent(in):: elem ! elem = element
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    real, dimension(1:ndim), intent(inout) :: jump
    real, dimension(:, :,:), allocatable :: wi     ! w recomputed  in integ nodes
    integer ::   Qnum, Qdof

    state%time%ctime = state%time%ttime

    !! seting of degree of the Gauss quadrature
    Qnum = elem%face(fGnum,ie)
    Qdof = state%space%G_rule(Qnum)%Qdof
    if(Qdof /= elem%face(fGdof,ie)) print*,'## Trouble in ElementEdgeJump'

    allocate(wi(1:3, 1:Qdof, 1:ndim ) )
    call ElementEdgeJump(elem,  ie, wi(1, 1:Qdof, 1:ndim) )
    call ElementEdgeJump_time(elem,  ie, wi(2, 1:Qdof, 1:ndim) )

    wi(3, 1:Qdof,1:ndim) = (wi(1, 1:Qdof,1:ndim)**2 + wi(2, 1:Qdof,1:ndim)**2 &
         + wi(1, 1:Qdof,1:ndim)* wi(2, 1:Qdof,1:ndim) ) /3

    wi(3, 1:Qdof,1:ndim) = wi(3, 1:Qdof,1:ndim) * elem%dn(ie)

    jump(1:ndim) =  matmul(state%space%G_rule(Qnum)%weights(1:Qdof), wi(3, 1:Qdof,1:ndim) )&
         *state%time%tau(1)

    deallocate(wi)

  end subroutine IntegElementEdgeJump_time



  !> evaluate \f$ \int_{\Gamma} [\nabla {\bf w}\cdot{\bf n}]^2\, dS \f$,
  !>
  !> \f$ \Gamma \f$ is an face shared by elements 'elem' and 'elem1'
  subroutine IntegElementEdgeJumpDer(elem,  ie, jump)
    type(element), intent(in):: elem ! elem = element
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    real, dimension(1:ndim), intent(inout) :: jump
    class(element), pointer ::   elem1  ! elem1 = neigh element
    real, dimension(:,:), allocatable :: Dw     ! Dw recomputed  in integ nodes
    real, dimension(:,:,:), allocatable :: Dwi     ! Dw recomputed  in integ nodes
    real, dimension(:,:,:), allocatable :: Dwii    ! Dw opposite or boundary state
    real, dimension(:), pointer :: xi
    integer ::  dof, dof1, ie1, Qnum, Qdof
    integer :: l, ii

    dof = elem%dof

    !! seting of degree of the Gauss quadrature
    Qnum = elem%face(fGnum,ie)
    Qdof = state%space%G_rule(Qnum)%Qdof
    if(Qdof /= elem%face(fGdof,ie)) print*,'## Trouble in ElementEdgeJumpDer'

    allocate(Dwi(1:Qdof, 1:ndim,1:nbDim ), Dwii(1:Qdof, 1:ndim,1:nbDim ) )
    call Eval_Dw_Edge(elem, ie, Dwi(1:Qdof, 1:ndim, 1:nbDim), .false.)

    ii = elem%face(neigh, ie)
    if( ii > 0) then  !! inner face
       elem1 => grid%elem(ii)
       ie1 = elem%face(nei_i,ie)

       if(Qdof /= elem1%face(fGdof,ie1)) print*,'## Trouble in ElementEdgeJumpDer (2)'

       call Eval_Dw_Edge(elem1, ie1, Dwii(1:Qdof, 1:ndim, 1:nbDim), .true.)

    else  !! boundary face

       !if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
       !   do l=1,Qdof
       !      xi => grid%b_edge(-elem%face(neigh,ie))%x_div(l, 1:nbDim)
       !      call Exact_Scalar(xi(1:nbDim), wii(l,1:ndim), state%time%ctime )
       !   enddo
       !else
          ! NOT YET IMPLEMENTED !!!!!!!!!
          Dwii(1:Qdof, 1:ndim, 1:nbDim) = Dwi(1:Qdof, 1:ndim, 1:nbDim)
       !endif
    endif


    !if(elem%i == 1) then
    !   write(*,'(a6,i5,12es14.6)') 'Dwi S', elem%i,Dwi(:,1,1)
    !   write(*,'(a6,i5,12es14.6)') 'Dwi S', elem%i,Dwi(:,1,2)
    !   write(*,'(a6,i5, 12es14.6)') 'Dwii S',ii,Dwii(:,1,1)
    !   write(*,'(a6,i5, 12es14.6)') 'Dwii S',ii,Dwii(:,1,2)
    !endif


    allocate(Dw(1:Qdof,1:ndim) )

    Dw(1:Qdof,1:ndim) = &
         (Dwi(1:Qdof,1:ndim, 1) - Dwii(1:Qdof,1:ndim, 1))*elem%n(ie, 1)  +  &
         (Dwi(1:Qdof,1:ndim, 2) - Dwii(1:Qdof,1:ndim, 2))*elem%n(ie, 2)

    !if(elem%i == 1) &
    !   write(*,'(a6,i5,12es14.6)') 'Dw S', elem%i, Dw(:,1), elem%n(ie,1:2)

    Dw(1:Qdof,1:ndim) = Dw(1:Qdof,1:ndim)**2

    !if(elem%i == 1) &
    !   write(*,'(a6,i5,12es14.6)') 'Dw A', elem%i, Dw(:,1)


    jump(1:ndim) =  matmul(state%space%G_rule(Qnum)%weights(1:Qdof), Dw(1:Qdof,1:ndim) ) &
         / elem%dn(ie)    ! since elem%n has length = h, it is squared

    deallocate( Dw, Dwi, Dwii)

  end subroutine IntegElementEdgeJumpDer

  !> evaluate \f$ \int_{\Gamma} [{\bf w}]^2\, dS \f$,
  !>
  !> \f$ \Gamma \f$ is an face shared by elements 'elem' and 'elem1'
  subroutine ElementEdgeJumpIndicator(elem,  ie)
    type(element), intent(inout):: elem ! elem = element
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    class(element), pointer ::   elem1  ! elem1 = neigh element
    real, dimension(:,:), allocatable :: phi, phi1 ! test functions
    real, dimension(:), allocatable :: wi       ! w recomputed  in integ nodes
    real, dimension(:,:,:), allocatable :: wiw
    real, dimension(:,:), allocatable :: wB     ! boundary state
    real, dimension(1:nbDim) :: xi
    real :: val

    integer ::  dof, dof1, ie1, kst, kst1, Qnum, Qdof
    integer :: l, l1, ii

    dof = elem%dof

    !! seting of degree of the Gauss quadrature
    Qnum = elem%face(fGnum,ie)
    Qdof = state%space%G_rule(Qnum)%Qdof

    allocate(wi(1:Qdof ) )
!    allocate(qL(1:Qdof ), qR(1:Qdof ) )

    allocate(phi(1:dof, 1:Qdof))
    call Eval_Phi_Edge(elem, dof, ie, phi, .false.)

    if(state%modelName == 'pedes' ) allocate(wiw(1:3, 1:3, 1:Qdof) )

    ii = elem%face(neigh, ie)

    if( ii > 0) then  !! inner face

       elem1 => grid%elem(ii)

       dof1 = elem1%dof
       ie1 = elem%face(nei_i,ie)

       allocate(phi1(1:dof1, 1:Qdof))
       call Eval_Phi_Edge(elem1, dof1, ie1, phi1, .true.)

       !print*,'##',ie,ie1,elem%flen, Qnum, Qdof, state%space%G_rule(Qnum)%Qdeg

       ! evaluation of (w_ie^+ - w_ie^-) in integ. nodes
       do l=1, Qdof       ! ndim ONLY the DENSITY
          wi(l) = (dot_product(phi(1:dof ,l),  elem%w(0,1: dof) ) &
               - dot_product(phi1(1:dof1 ,l), elem1%w(0,1: dof1) ) )**2
       enddo !! l

       ! if(state%modelName == 'pedes' ) then
       !    do l=1, Qdof       ! indicators = jumps of the velocity
       !       ! left values
       !       wiw(1, 1, l) = dot_product(phi(1:dof ,l),  elem%w(0, 0*dof + 1: 1*dof) ) ! rho
       !       wiw(1, 2, l) = dot_product(phi(1:dof ,l),  elem%w(0, 1*dof + 1: 2*dof) ) ! rho*u
       !       wiw(1, 3, l) = dot_product(phi(1:dof ,l),  elem%w(0, 2*dof + 1: 3*dof) ) ! rho*v
       !       wiw(1, 2:3, l)  = wiw(1, 2:3, l) / wiw(1, 1, l)                           ! u, v  
             
       !       ! right values
       !       wiw(2, 1, l) = dot_product(phi1(1:dof1 ,l),  elem1%w(0, 0*dof1 + 1: 1*dof1) ) ! rho
       !       wiw(2, 2, l) = dot_product(phi1(1:dof1 ,l),  elem1%w(0, 1*dof1 + 1: 2*dof1) ) ! rho*u
       !       wiw(2, 3, l) = dot_product(phi1(1:dof1 ,l),  elem1%w(0, 2*dof1 + 1: 3*dof1) ) ! rho*v
       !       wiw(2, 2:3, l)  = wiw(2, 2:3, l) / wiw(2, 1, l)                           ! u, v  
             
       !       ! square of the difference
       !       wiw(3, 1:3, l) = ( wiw(1, 1:3, l) - wiw(2, 1:3, l) )**2
             
       !       ! inserting in the array the square of thye magnitude of the velocity
       !       !wi(l) = wiw(3,2,l) + wiw(3,3,l)

       !    enddo !! l
       ! endif  ! pedes


       deallocate(phi1)
    else  !! boundary face

       ! evaluation of w_ie^L in integ. nodes
       do l=1, Qdof       ! ndim ONLY the DENSITY
          wi(l) = dot_product(phi(1:dof ,l),  elem%w(0,1: dof) )
       enddo

       !boundary state    ! ONLY first component (density) !!!!!!!
       allocate(wB(1:Qdof, 1:ndim ) )

       if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
          do l=1,Qdof
             xi(1:nbDim) = grid%b_edge(-elem%face(neigh,ie))%x_div(l, 1:nbDim)
             call Exact_Scalar(xi(1:nbDim), wB(l,1:ndim), state%time%ctime )
          enddo
       else
          ! NOT YET IMPLEMENTED !!!!!!!!!
          wB(1:Qdof, 1) = wi(1:Qdof)
       endif

       ! ONLY first component (density) !!!!!!!
       wi(1:Qdof) = (wi(1:Qdof) - wB(1:Qdof, 1) )**2

       deallocate(wB)
    endif

    elem%rezid = elem%rezid + dot_product(wi(:), state%space%G_rule(Qnum)%weights(:) ) &
         * elem%dn(ie)

    !if(elem%i == 1) then
    !   write(*,'(a8, 3i5, 20es12.4)') 'jump:', elem%i, ie, l, elem%rezid, &
    !        dot_product(wi(:), state%space%G_rule(Qnum)%weights(:) ) * elem%dn(ie), &
    !        wi(1),  wiw(3, 1:3, 1)
    !endif

    deallocate(wi, phi)
    if(state%modelName == 'pedes' ) deallocate(wiw )

  end subroutine ElementEdgeJumpIndicator

  !> evaluate \f$ \int_{\partial K} [{\bf v}]^2\, dS \f$,
  subroutine ElementJumpVelocityIndicator(elem) 
     type(element), intent(inout):: elem 
     integer :: j,k

     elem%rezid = 0.

     do j = 1, elem%flen
        k = elem%face(neigh,j)
        if(k > 0) then
           !elem1 => grid%elem(k)
           !call ElementEdgeJumpIndicator(elem,  j)
           call ElementEdgeJumpVelocityIndicator(elem,  j)
           !call ElementEdgeDerJumpIndicator(elem, elem1, j)
        endif
     enddo
   end subroutine ElementJumpVelocityIndicator


  !> evaluate \f$ \int_{\Gamma} [{\bf v}]^2\, dS \f$,
  !>
  !> \f$ \Gamma \f$ is an face shared by elements 'elem' and 'elem1'
  subroutine ElementEdgeJumpVelocityIndicator(elem,  ie)
    type(element), intent(inout):: elem ! elem = element
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    class(element), pointer ::   elem1  ! elem1 = neigh element
    real, dimension(:,:), allocatable :: phi, phi1 ! test functions
    real, dimension(:), allocatable :: wi       ! w recomputed  in integ nodes
    real, dimension(:,:,:), allocatable :: wiw
    real, dimension(:,:), allocatable :: wB     ! boundary state
    real, dimension(1:nbDim) :: xi
    real :: val

    integer ::  dof, dof1, ie1, kst, kst1, Qnum, Qdof
    integer :: l, l1, ii

    dof = elem%dof

    !! seting of degree of the Gauss quadrature
    Qnum = elem%face(fGnum,ie)
    Qdof = state%space%G_rule(Qnum)%Qdof

    allocate(wi(1:Qdof ) )
!    allocate(qL(1:Qdof ), qR(1:Qdof ) )

    allocate(phi(1:dof, 1:Qdof))
    call Eval_Phi_Edge(elem, dof, ie, phi, .false.)

    allocate(wiw(1:3, 1:3, 1:Qdof) ) ! velocity array

    ii = elem%face(neigh, ie)

    if( ii > 0) then  !! inner face

       elem1 => grid%elem(ii)

       dof1 = elem1%dof
       ie1 = elem%face(nei_i,ie)

       allocate(phi1(1:dof1, 1:Qdof))
       call Eval_Phi_Edge(elem1, dof1, ie1, phi1, .true.)

       !print*,'##',ie,ie1,elem%flen, Qnum, Qdof, state%space%G_rule(Qnum)%Qdeg

       ! evaluation of (w_ie^+ - w_ie^-) in integ. nodes
       do l=1, Qdof       ! ndim ONLY the DENSITY
          wi(l) = (dot_product(phi(1:dof ,l),  elem%w(0,1: dof) ) &
               - dot_product(phi1(1:dof1 ,l), elem1%w(0,1: dof1) ) )**2
       enddo !! l

       ! if(state%modelName == 'pedes' ) then
       do l=1, Qdof    !!!   indicators = jumps of the velocity
          ! left values
          wiw(1, 1, l) = dot_product(phi(1:dof ,l),  elem%w(0, 0*dof + 1: 1*dof) ) ! rho
          wiw(1, 2, l) = dot_product(phi(1:dof ,l),  elem%w(0, 1*dof + 1: 2*dof) ) ! rho*u
          wiw(1, 3, l) = dot_product(phi(1:dof ,l),  elem%w(0, 2*dof + 1: 3*dof) ) ! rho*v
          if(wiw(1, 1, l) > state%model%Pr) then
             wiw(1, 2:3, l)  = wiw(1, 2:3, l) / wiw(1, 1, l)                           ! u, v  
          else
             wiw(1, 2:3, l)  = 0.
          endif

          ! right values
          wiw(2, 1, l) = dot_product(phi1(1:dof1 ,l),  elem1%w(0, 0*dof1 + 1: 1*dof1) ) ! rho
          wiw(2, 2, l) = dot_product(phi1(1:dof1 ,l),  elem1%w(0, 1*dof1 + 1: 2*dof1) ) ! rho*u
          wiw(2, 3, l) = dot_product(phi1(1:dof1 ,l),  elem1%w(0, 2*dof1 + 1: 3*dof1) ) ! rho*v

          if(wiw(2, 1, l) > state%model%Pr) then
             wiw(2, 2:3, l)  = wiw(2, 2:3, l) / wiw(2, 1, l)                           ! u, v  
          else
             wiw(2, 2:3, l)  = 0.
          endif


          ! square of the difference
          wiw(3, 1:3, l) = ( wiw(1, 1:3, l) - wiw(2, 1:3, l) )**2

          ! inserting in the array the square of thye magnitude of the velocity
          wi(l) = wiw(3,2,l) + wiw(3,3,l)

       enddo !! l
       ! endif  ! pedes


       deallocate(phi1)
    else  !! boundary face

       ! evaluation of w_ie^L in integ. nodes
       do l=1, Qdof       ! ndim ONLY the DENSITY
          wi(l) = dot_product(phi(1:dof ,l),  elem%w(0,1: dof) )
       enddo

       !boundary state    ! ONLY first component (density) !!!!!!!
       allocate(wB(1:Qdof, 1:ndim ) )

       if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
          do l=1,Qdof
             xi(1:nbDim) = grid%b_edge(-elem%face(neigh,ie))%x_div(l, 1:nbDim)
             call Exact_Scalar(xi(1:nbDim), wB(l,1:ndim), state%time%ctime )
          enddo
       else
          ! NOT YET IMPLEMENTED !!!!!!!!!
          wB(1:Qdof, 1) = wi(1:Qdof)
       endif

       ! ONLY first component (density) !!!!!!!
       wi(1:Qdof) = (wi(1:Qdof) - wB(1:Qdof, 1) )**2

       deallocate(wB)
    endif

    elem%rezid = elem%rezid + dot_product(wi(:), state%space%G_rule(Qnum)%weights(:) ) &
         * elem%dn(ie)

    !if(elem%i == 1) then
    !   write(*,'(a8, 3i5, 20es12.4)') 'jump:', elem%i, ie, l, elem%rezid, &
    !        dot_product(wi(:), state%space%G_rule(Qnum)%weights(:) ) * elem%dn(ie), &
    !        wi(1),  wiw(3, 1:3, 1)
    !endif

    deallocate(wi, phi)
    !if(state%modelName == 'pedes' ) deallocate(wiw )

  end subroutine ElementEdgeJumpVelocityIndicator



  !> evaluate \f$ \int_{\Gamma} [{\bf Der w}]^2\, dS \f$,
  !>
  !> \f$ \Gamma \f$ is an face shared by elements 'elem' and 'elem1'
  subroutine ElementEdgeDerJumpIndicator(elem, elem1, ie)
    type(element), intent(inout):: elem, elem1  ! elem = element, elem1 = neigh element
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    real, dimension(:,:), allocatable :: w, w1, wa   ! w recomputed  in integ nodes
    real, dimension(:,:,:), allocatable :: Dw, Dw1   ! Dw recomputed  in integ nodes
    real, dimension(:), allocatable :: val           ! w recomputed  in integ nodes

    integer :: Qdof, ie1, l, l1


    ie1 = elem%face(nei_i,ie)
    Qdof = elem%face(fGdof,ie)

    if(Qdof .ne. elem1%face(fGdof,ie1)) then
       print*,'Incompatible Gdof in ElementEdgeDerJumpIndicator'
       print*, elem%face(fGdof,ie), elem1%face(fGdof,ie1)
       stop
    endif

    allocate(w(1:Qdof, 1:ndim ) )
    allocate(w1(1:Qdof, 1:ndim ) )

    call Eval_w_Edge(elem, ie, w, .false.)
    call Eval_w_Edge(elem1, ie1, w1, .true.)

    allocate(wa(1:Qdof, 1:ndim ) )

    call Eval_aver_w_Edge(elem, elem1, ie,Qdof, wa)

    allocate(Dw(1:Qdof, 1:ndim, 1:nbDim ) )
    allocate(Dw1(1:Qdof, 1:ndim, 1:nbDim ) )

    call Eval_Dw_Edge(elem, ie, Dw, .false.)


    call Eval_Dw_Edge(elem1, ie1, Dw1, .true.)

    ! evaluation of \nabla (w_ie^+ - w_ie^-) in integ. nodes, ONLY the DENSITY

    !print*,'$$$$$$$$$$    ', Qdof, elem%face(fGdof,ie), elem1%face(fGdof,ie1)
    !
    !print*,'Dwx :', Dw(1:Qdof, 1, 1)
    !print*,'Dwx :', Dw1(1:Qdof, 1, 1)
    !print*
    !
    !print*,'Dwy :', Dw(1:Qdof, 1, 2)
    !print*,'Dwy :', Dw1(1:Qdof, 1, 2)
    !print*

    allocate(val(1:Qdof) )

    ! jump of the size gradient
    !val(1:Qdof) = (Dw(1:Qdof, 1, 1)-Dw1(1:Qdof, 1, 1) ) &
    !     * (Dw(1:Qdof, 1, 1)-Dw1(1:Qdof, 1, 1) ) &
    !     + (Dw(1:Qdof, 1, 2)-Dw1(1:Qdof, 1, 2) ) &
    !     * (Dw(1:Qdof, 1, 2)-Dw1(1:Qdof, 1, 2) )

    ! jump of the grandient multiplied by momentum
    !val(1:Qdof) = (Dw(1:Qdof, 1, 1)-Dw1(1:Qdof, 1, 1) ) * wa(1:Qdof, 2 ) &
    !     + (Dw(1:Qdof, 1, 2)-Dw1(1:Qdof, 1, 2) ) * wa(1:Qdof, 3)


    ! density indicator:
    !val(1:Qdof)= (w1(1:Qdof,1) - w(1:Qdof,1)) &
    !     * (elem%n(ie,1) *  wa(1:Qdof, 2 )/ wa(1:Qdof, 1 ) &
    !     +  elem%n(ie,2) *  wa(1:Qdof, 3 )/ wa(1:Qdof, 1 ) )/  elem%dn(ie)


    !if(dot_product(val(:), val(:) )> 1E-4) then
    !   print*,'rho_i:', w(1:Qdof,1)
    !   print*,'rho_j:', w1(1:Qdof,1)
    !   print*,'nn=  ', elem%n(ie,1:nbDim)
    !   print*,'v_1 :', wa(1:Qdof, 2 )/ wa(1:Qdof, 1)
    !   print*,'v_2 :', wa(1:Qdof, 3 )/ wa(1:Qdof, 1)
    !   print*,'v*n',  (elem%n(ie,1) *  wa(1:Qdof, 2 )/ wa(1:Qdof, 1 ) &
    !        +  elem%n(ie,2) *  wa(1:Qdof, 3 )/ wa(1:Qdof, 1 ) )/  elem%dn(ie)
    !   print*,'Val :', val(1:Qdof)
    !print*,'-------------------',dot_product(val(:), val(:) ),'*****', &
    !     (elem%xc(:)+elem1%xc(:) )/2
    !endif


    ! nabla density indicator:
    if(ndim == 4) then
       val(1:Qdof)= (Dw1(1:Qdof,1,1) + Dw(1:Qdof,1,1)) * wa(1:Qdof, 2)/ wa(1:Qdof, 1) &
            + (Dw1(1:Qdof,1,2) + Dw(1:Qdof,1,2)) * wa(1:Qdof, 3)/ wa(1:Qdof, 1)
    else
       Dw(1:Qdof,1,1:nbDim) = Dw1(1:Qdof,1,1:nbDim) + Dw(1:Qdof,1,1:nbDim)
       val(1:Qdof)= Dw(1:Qdof,1,1) * Dw(1:Qdof,1,1)
    endif

    do l=1,Qdof
       val(l) = max(0., val(l) )
    enddo


    !!if(dot_product(val(:), val(:) )> 1E+2) then
    !   print*,'Dwx:', (Dw1(1:Qdof,1,1) + Dw(1:Qdof,1,1))
    !   print*,'Dwy:', (Dw1(1:Qdof,1,2) + Dw(1:Qdof,1,2))
    !
    !   print*,'v_1 :', wa(1:Qdof, 2 )/ wa(1:Qdof, 1)
    !   print*,'v_2 :', wa(1:Qdof, 3 )/ wa(1:Qdof, 1)
    !   print*,'###: ',( Dw1(1:Qdof,1,1) + Dw(1:Qdof,1,1)) * wa(1:Qdof, 2)/ wa(1:Qdof, 1) &
    !     + (Dw1(1:Qdof,1,2) + Dw(1:Qdof,1,2)) * wa(1:Qdof, 3)/ wa(1:Qdof, 1)
    !   print*,'Val :', val(1:Qdof)
    !print*,'-------------------',dot_product(val(:), val(:) ),'*****', &
    !     (elem%xc(:)+elem1%xc(:) )/2
    !endif

    elem%rezid = elem%rezid + elem%dn(ie) &
         * dot_product(val(1:Qdof), state%space%G_rule(elem%face(fGnum,ie))%weights(1:Qdof) )

    !print*,'Rez =', elem%i, ie, elem1%i, elem%rezid,  &
    !      dot_product(Dw(1:Qdof,1,1), state%space%G_rule(elem%face(fGnum,ie))%weights(1:Qdof) )

    !print*,'-------------------',dot_product(val(:), val(:) ),'*****', &
    !     (elem%xc(:)+elem1%xc(:) )/2


    deallocate(Dw, Dw1, w, w1, wa, val)

  end subroutine ElementEdgeDerJumpIndicator

  !> evaluate the jumps of the solution and all its projections
  subroutine  ElementEdgeJumpsAllProj(elem, regularity)
    type(element), target, intent(in):: elem ! elem = element
    real,  intent(inout) :: regularity
    class(element), pointer ::   elem1 
    real, dimension(:), pointer :: weights
    real, dimension(:,:), allocatable :: jumps
    real, dimension(:,:), allocatable :: wi, wii
    real, dimension(:,:,:), allocatable :: wwP   ! projections
    real, dimension(:,:), allocatable:: phi ! test functions
    integer :: ii, ie, ie1, Qnum, Qdof, dof, dofL, deg, k, kst, l, ifile, ifile2, ifile3
    real :: sumy, sumyi, sumx, sumx2, sumy2, sumxy
    integer :: icase, ni, min_deg, max_deg
    real, dimension (:, :), allocatable :: ss
    real :: val, val1
    logical :: iprint, singularity

    iprint = .false.
    !iprint = .true.
    
    !call Detect_apriori_known_singularity(elem, singularity)
    !if(singularity) iprint = .true.

    allocate(jumps(0:elem%deg, 0:ndim ) ) 
    jumps(:,:) = 0.

    dof = elem%dof

    ! projection of the solution of the element
    allocate(wwP(0:elem%deg, 1:dof, 1:ndim) )  
    wwP = 0.
    
    ! least squares for decays of energy projections in H^1-seminorm
    sumx = 0.; sumy = 0.; sumx2 = 0.; sumy2 = 0.; sumxy = 0.; ni = 0

    ! projections of all polynomial  degrees
    !!!do icase = 1, 2  
    icase = 1  ! =1 projection in the L^2-norm, =2 energy projection

    do deg = elem%deg, 0, - 1

       if(icase == 2)  call Energy_Elem_deg_projection(elem, deg, wwP(deg, 1:dof, 1:ndim) )

       if(icase == 1) then
          dofL =  (deg + 1) * (deg + 2) / 2
          
          do k=1, ndim
             wwP(deg, 1:dofL, k) = elem%w(0, (k-1)*dof + 1 :  (k-1)*dof + dofL) 
          enddo
       endif
       
       !if(iprint) then
       !   write(*,'(a6,5i5, 50es12.4)') 'proj W:',elem%i, elem%deg, deg, dof, dofL,wwP(deg, 1:dof, 1)
       !endif
    enddo
    !if(iprint)  print*

    ! maximal projection degree used for least squares
    max_deg = elem%deg

    do ie = 1, elem%flen
       ii = elem%face(neigh, ie)

       if(ii > 0) then ! interior edge
          !! seting of degree of the Gauss quadrature
          Qnum = elem%face(fGnum,ie)
          Qdof = state%space%G_rule(Qnum)%Qdof
          if(Qdof /= elem%face(fGdof,ie)) print*,'## Trouble in ElementEdgeJumpsAllProj'
          
          weights => state%space%G_rule(Qnum)%weights(1:Qdof)
          
          allocate(wi(1:Qdof, 1:ndim), wii(1:Qdof, 1:ndim) )
          
          ! neighbouring element
          elem1 => grid%elem(ii)
          
          ! maximal projection degree used for least squares
          max_deg = min(max_deg, elem1%deg)
             
          !neighbouring elements in integ nodes
          ie1 = elem%face(nei_i,ie)
             
          ! value of the solution from outside 
          call Eval_w_Edge(elem1, ie1, wii(1:Qdof, 1:ndim), .true.)
          
          ! test function of the edge from the given element
          allocate(phi(1:dof, 1:Qdof))
          call Eval_Phi_Edge(elem, dof, ie, phi(1:dof, 1:Qdof), .false.)
          
          ! solution from inside
          do k=1,ndim              ! k = index of component of w
             kst = dof*(k-1) + 1
             
             do l=1,Qdof
                wi(l,k) = dot_product(elem%w(0,kst:kst+dof-1), phi(1:dof,l) )
             enddo
          enddo
          
          ! solution on the edge, mean value of the approximate solution
          wii(1:Qdof,1:ndim) = ( wii(1:Qdof,1:ndim) + wi(1:Qdof,1:ndim) ) / 2

          do deg = elem%deg, 0, - 1
             
             dofL = (deg + 1) * (deg + 2) / 2
             
             do k=1,ndim              ! k = index of component of w
                kst = dof*(k-1) + 1
                
                do l=1,Qdof
                   !wi(l,k) = dot_product(elem%w(0,kst:kst+dofL-1), phi(1:dofL,l) )
                   wi(l,k) = dot_product(wwP(deg, 1:dofL, k), phi(1:dofL,l) )
                enddo
             enddo

             ! square of jumps in integ nodes
             wi(1:Qdof, 1:ndim) = ( wi(1:Qdof,1:ndim) - wii(1:Qdof,1:ndim))**2
             
             ! integrals over edge
             do k=1,ndim
                jumps(deg, k) = jumps(deg, k) + &
                     dot_product(wi(1:Qdof, k), weights(1:Qdof)) * elem%dn(ie)
             enddo
             jumps(deg, 0) = jumps(deg, 0) + sum( jumps(deg, 1:ndim)) 
             
          enddo

          deallocate(phi, wi, wii)

       endif
    enddo ! ie = 1, elem%flen

    jumps(0:elem%deg, 0:ndim) = sqrt( jumps(0:elem%deg, 0:ndim) )


    ! estimate of the regularity
    !allocate(ss(1:elem%deg, 1:2) )
    !ss(:,:)= 0.
    !do deg = elem%deg, 2, - 1
    !   ss(deg, 1) = 0.5+(log(jumps(deg, 0) / jumps(deg-1, 0)) - log(elem%diam) ) / log(1.*(deg -1) / deg)
    !   ss(deg, 2) = 0.5+(log(jumps(deg, 0) / jumps(deg-1, 0)) ) / log(1.*(deg -1) / deg)
    !enddo
    
    min_deg = 1
    max_deg = max_deg 

    do deg = max_deg, min_deg, - 1
       ! least squares:
       if(deg <= elem%deg .and. deg >= 1 ) then
          ni = ni + 1
          sumx = sumx + log(1. * deg)
          sumx2 = sumx2 + log(1. * deg)**2
          sumy = sumy + log(jumps(deg,0))
          sumy2 = sumy2 + log(jumps(deg,0) )**2
          sumxy = sumxy + log(1. * deg) * log(jumps(deg,0))
       endif
       !if(singularity .and. state%space%adapt%adapt_level == 4) &
       !     write(*,'(a8, i5, 30es12.4)') 'L R 23:',ni, 1.*deg, jumps(deg,0), sumx, sumx2, sumy, sumy2, sumxy
    enddo
       
    if( ni * sumx2-sumx*sumx /= 0.) then
       regularity =  - (ni*sumxy -sumx*sumy)/(ni*sumx2-sumx*sumx) + 0.5
    else
       regularity = -1.
    endif
    
    if(singularity) then
       ifile = state%space%adapt%adapt_level + 150
    else
       ifile = state%space%adapt%adapt_level + 100
    endif
    ifile2 = ifile + 100
    
    if(iprint) then
       !write(10, *) elem%i, elem%deg, jumps(0:elem%deg, 0) 
       do deg = elem%deg, 0, -1
          write(ifile, *) elem%i, elem%deg, elem%i - 0.1 * (elem%deg - deg), elem%xc(1:2), &
               jumps(deg, 0)
          !if(deg >=2) &
          !     write(ifile2, *) elem%i, elem%deg, elem%i - 0.1 * (elem%deg - deg), elem%xc(1:2), &
          !     jumps(deg, 0), ss(deg, 1:2),  regularity
       enddo
       write(ifile, '(x)') 
       write(ifile, '(x)') 
       write(ifile, '(x)') 
       
       !write(ifile2, '(x)') 
       !write(ifile2, '(x)') 
       !write(ifile2, '(x)') 
       
       deg = elem%deg
       write(ifile2, *) elem%i, elem%deg, elem%i - 0.1 * (elem%deg - deg), elem%xc(1:2), &
            regularity, max_deg, min_deg, ni
       
       ! do deg = elem%deg, 0, -1
       !    k= deg*(deg+1)/2 + 1
       !    l = (deg+1)*(deg+2)/ 2
       !    val =  sqrt( dot_product( wwP(elem%deg, k:l, 1),  wwP(elem%deg, k:l, 1)) )
       
       !    write(500+elem%i, *) elem%i, elem%deg, elem%deg - 0.1 * (elem%deg - deg), &
       !         val
       ! enddo
       !  write(500+elem%i, '(x)') 
    endif
    
    !deallocate(ss)
    deallocate( wwP, jumps)
    
    !enddo  ! icase=1,2
     
  end subroutine ElementEdgeJumpsAllProj

end module eval_jumps
