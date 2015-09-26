!> main subroutines for the solution of the unstationary compressible
!> Navier-Stokes equations

module euler_problem
  use main_data
  use problem_oper
  use lin_solvers
  use inviscid_fluxes
  use time_sets
  use stdgm_mod
  use matrix_oper_int


  implicit none

  public:: ComputeTerms
  public:: ComputeElementsTerms
  public:: ComputeOneElementFluxes
  public:: ComputeOneElementTerms

  public:: ComputeSTDGM_Terms


  public:: FillVector
  public:: FillVectorST
  public:: PrepareCrankNicolson

  public:: ElementEdgeJump_time
  public:: ElementEdgeJump
  public:: ElementEdgeJumpProj

  public:: IntegElementEdgeJump
  public:: IntegElementEdgeJump_time
  public:: IntegElementEdgeJumpDer

  public:: ElementEdgeJumpIndicator
  public:: ElementEdgeDerJumpIndicator
  public:: SmoothResid
  public:: EquivResid
  public:: EquivResid1

  public:: ReadConvFile

  public:: CheckResiduum

  public:: ClearMatrixBlocks
  public:: ClearMatrixBlocks_STDGM
  public:: ClearVectorBlocks
  public:: ClearVectorBlocks_STDGM

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
  subroutine ElementEdgeJumpProj(elem,  ie,  jump, jump1)
    type(element), intent(in):: elem ! elem = element
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    real, dimension(1:elem%face(fGdof,ie), 1: ndim), intent(out) :: jump, jump1
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
       elem1 => grid%elem(ii)
       ie1 = elem%face(nei_i,ie)

       if(Qdof /= elem1%face(fGdof,ie1)) print*,'## Trouble in ElementEdgeJump (2)'

       call Eval_w_EdgeProj(elem1, ie1, wii(1,1:Qdof, 1:ndim), wii(2,1:Qdof, 1:ndim), .true.)

    else  !! boundary face

       if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
          do l=1,Qdof
             xi => grid%b_edge(-elem%face(neigh,ie))%x_div(l, 1:nbDim)
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

    deallocate(wi, phi)

  end subroutine ElementEdgeJumpIndicator



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


  !> evaluate \f$ \int_{K} [{\bf D w}]^2\, dx \f$,
  subroutine ElementGradIndicator(elem)
    type(element), intent(inout):: elem  ! elem = element
    real, dimension(:), allocatable :: Mach   ! w recomputed  in integ nodes
    real, dimension(:,:), allocatable :: wi   ! w recomputed  in integ nodes
    real, dimension(:,:,:), allocatable :: Dwi   ! w recomputed  in integ nodes
    real, dimension(:), allocatable     :: weights     !> weights in integ nodes
    !real, dimension(:,:,:), allocatable :: Der       ! Der recomputed  in integ nodes
    real, dimension(1:nbDim) :: D_rho, velocity
    integer :: Qdof, dof, i, k
    real :: ha

    Qdof = elem%Qdof
    dof = elem%dof

    ! setting of integ. weights
    allocate(weights(1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof) )

    allocate(wi(1:Qdof, 1:ndim), Mach(1:Qdof))
    call Eval_w_Elem(elem, wi)
    call EvalmachNumber(Qdof, wi, Mach)

    allocate(Dwi(1:Qdof, 1:ndim, 1:nbDim))
    call Eval_Dw_Elem(elem, Dwi)

    do i=1,2
       D_rho(i) = dot_product(weights(1:Qdof), Dwi(1:Qdof, 1, i) ) /elem%area

       velocity(i) = dot_product(weights(1:Qdof), wi(1:Qdof, i+1)/wi(1:Qdof, 1) ) &
            /elem%area
    enddo

    elem%rezid = max(0., dot_product(D_rho(1:nbDim), velocity(1:nbDim) ) )

    !elem%rezid = min(10., dot_product(D_rho(1:nbDim), velocity(1:nbDim) ) )

    ! Mach number over 1, necessary condition
    !if( dot_product(weights(1:Qdof), Mach(1:Qdof) ) /elem%area < 1.)  elem%rezid =0.
    !if( maxval( Mach(1:Qdof) ) < 1.)  elem%rezid =0.

    ha = 0.
    !if(elem%HGnode) then
    !   do i=1, elem%ftype
    !   enddo
    !else
    do i=1,elem%flen
       ha = max(ha, abs((D_rho(1)*elem%n(i,2) - D_rho(2)*elem%n(i,1)) &
            /dot_product(D_rho(1:nbDim),D_rho(1:nbDim) )**0.5 ) )
    enddo

    !if(elem%rezid > 0) &
    !     write(200+state%time%iter,*)&
    !     elem%xc(1:nbDim), elem%rezid, dot_product(D_rho(1:nbDim), D_rho(1:nbDim) ), &
    !     dot_product(velocity(1:nbDim), velocity(1:nbDim) ), &
    !     elem%rezid*(ha/(elem%deg +1)**2), ha, ha/elem%diam, &
    !     dot_product(weights(1:Qdof), Mach(1:Qdof) ) /elem%area,elem%rezid

    ! amount of "numerical viscosity"
    !if(elem%rezid > 5E-01) then
    if(elem%rezid > 2E-00) then

       !elem%rezid = abs(state%model%Re) * 2. * elem%rezid * ( ha/(elem%deg + 1))**2.
       !elem%rezid = abs(state%model%Re) * ( ha/(elem%deg + 1))**2.

       !elem%rezid = abs(state%model%Re)  *    2.      * ( ha/2)**2. /(elem%deg + 1)
       elem%rezid = abs(state%model%Re) * elem%rezid * ( ha/2)**2. /(elem%deg + 1)
       !elem%rezid = abs(state%model%Re) * elem%rezid * ( ha/2)**2.
    else
       elem%rezid = 0.
    endif

    deallocate(Dwi, weights)
    !deallocate(Der)

  end subroutine ElementGradIndicator

  subroutine SmoothResid()
    real, dimension(:), allocatable :: rez
    integer, dimension(:), allocatable :: coun
    integer :: i,j, ii

    allocate(rez(1:grid%nelem), coun(1:grid%nelem) )
    rez(:) = 0.
    coun(:) = 0

    do i=1,grid%nelem
       rez(i) = rez(i) + 4.*grid%elem(i)%rezid
       coun(i) = coun(i) + 4
       do j=1,grid%elem(i)%flen
          if(grid%elem(i)%face(neigh,j) > 0) then
             ii = grid%elem(i)%face(neigh,j)
             rez(ii) = rez(ii) + grid%elem(i)%rezid
             coun(ii) = coun(ii) + 1
          endif
       enddo
    enddo

    do i=1,grid%nelem
       grid%elem(i)%rezid = rez(i) / coun(i)

       !if(rez(i) > 0.) then
       !   do j=0,3
       !      write(200+state%time%iter,*) &
       !           grid%x(grid%elem(i)%face(idx,mod(j,3)+1), 1:nbDim ), grid%elem(i)%rezid
       !   enddo
       !   write(200+state%time%iter,*) '   '
       !   write(200+state%time%iter,*) '   '
       !   write(200+state%time%iter,*) '   '
       !endif

    enddo

    deallocate(rez, coun)

  end subroutine SmoothResid

  !> "equilibriate the residuum" for smoother resolution,
  !> resulting distribution is constant on subdomains
  subroutine EquivResid()
    class(element), pointer :: elem, elem1
    integer :: i,j, k
    logical :: change

    do k=1, 100
       change = .false.

       do i=1,grid%nelem
          elem => grid%elem(i)

          if(elem%rezid > 0.) then
             do j=1,elem%flen

                if(elem%face(neigh,j) > 0) then
                   elem1 => grid%elem(elem%face(neigh,j))

                   if(elem1%rezid > 0.) then
                      ! both neighbouring elements were detected, we take maximum over them
                      elem%rezid = max(elem%rezid, elem1%rezid)
                      elem1%rezid = elem%rezid
                      change = .true.
                   endif

                endif
             enddo
          endif
       enddo
       if(.not. change) goto 10
    enddo

    10 continue

  end subroutine EquivResid


  !> "equilibriate the residuum" for smoother resolution
  !> resulting distribution should me "monotone"
  subroutine EquivResid1()
    class(element), pointer :: elem, elem1
    integer :: i,j, k, ichange
    real :: rmin, rmax
    logical :: change

    do k=1, 90
       change = .false.
       ichange = 0

       do i=1,grid%nelem
          elem => grid%elem(i)

          if(elem%rezid > 0.) then
             rmin = 1E+20
             rmax = 0.
             do j=1,elem%flen

                if(elem%face(neigh,j) > 0) then
                   elem1 => grid%elem(elem%face(neigh,j))

                   rmin = min(rmin, elem1%rezid)
                   rmax = max(rmax, elem1%rezid)

                endif
             enddo

             if(elem%rezid > rmax) then

                elem%rezid = rmax
                change = .true.

                ichange = ichange + 1

             elseif(elem%rezid < rmin) then

                !write(700+k,*) elem%xc(1:nbDim), elem%rezid,elem%i
                !write(700+k,*) elem%xc(1:nbDim), rmin
                !write(700+k,*) '  '
                !write(700+k,*) '  '

                elem%rezid = rmin
                change = .true.

                ichange = ichange + 1

             endif

          endif
       enddo
       if(.not. change) goto 10

       !do i=1,grid%nelem
       !   if(grid%elem(i)%rezid > 0.) &
       !        write(500+k,*) grid%elem(i)%xc(1:nbDim), grid%elem(i)%rezid
       !enddo

    enddo

    10 continue

  end subroutine EquivResid1


  !> setting of the artificial viscosity in elem%vec(aRe, *) from elem%rezid
  !> \f$ P_1 \f$ continuous reconstruction
  subroutine SetArtificialViscosity( )
    class(element), pointer :: elem
    real, dimension(:), allocatable :: rez
    integer, dimension(:), allocatable :: coun
    real, dimension(:), allocatable :: w_lag  ! viscosity in Lagr. node
    real, dimension(:), allocatable :: w_bas  ! viscosity in basis coeff
    integer :: i,j, k, Ndeg

    allocate(rez(1:grid%npoin), coun(1:grid%npoin) )

    Ndeg = 1 ! only P_1 lagr. approximation
    allocate(w_lag(1:(Ndeg+1)*(Ndeg+2)/2) )
    allocate(w_bas(1: state%space%max_dof ) )

    rez(:) = 0.
    coun(:) = 0

    do i=1,grid%nelem
       elem => grid%elem(i)
       do j=1,elem%flen
          k = elem%face(idx, j)
          rez(k) = rez(k) + elem%rezid
          coun(k) = coun(k) + 1
       enddo
    enddo

    do k=1,grid%npoin
       rez(k) = rez(k) / coun(k)
       !if(rez(k) > 0.) write(300,*) grid%x(k, 1:nbDim), rez(k)
    enddo

    do i=1,grid%nelem
       elem => grid%elem(i)
       if(elem%HGnode) then
          w_lag(1:3) = rez(elem%face(idx, elem%HGvertex(1:3) ) )
       else
          w_lag(1:3) = rez(elem%face(idx, 1:3) )
       endif

       !if(sum(w_lag(1:3) )/= 0.) then
       !   if(elem%HGnode) then
       !      do j=0,3
       !         write(100+state%time%iter,*) &
       !              grid%x(elem%face(idx,elem%HGvertex(mod(j,3)+1)), 1:nbDim ),&
       !              w_lag(mod(j,3)+1)
       !      enddo
       !   else
       !      do j=0,3
       !         write(100+state%time%iter,*) &
       !              grid%x(elem%face(idx,mod(j,3)+1), 1:nbDim ), w_lag(mod(j,3)+1)
       !      enddo
       !   endif
       !   write(100+state%time%iter,*)
       !   write(100+state%time%iter,*)
       !endif

       call Lagr2BasisDiffOne(elem, Ndeg,  w_lag, elem%Qnum, elem%dof, w_bas(1:elem%dof) )

       elem%vec(aRe, 1:elem%dof) = w_bas(1:elem%dof)

    enddo

    deallocate(rez, coun)

  end subroutine SetArtificialViscosity


  subroutine ReadConvFile(convfile)
    character(len=*), intent(in) :: convfile
    integer :: dummy_i
    real, dimension(1:6) :: dummy_r     !5 - new_err, 6 - log10(new_err)
    integer :: ifile = 11
    integer :: i,j,k

    open(ifile, file=convfile, status='OLD')

    do i= 1, state%time%iter
       read(ifile, *, end = 100) dummy_i, dummy_r(1:6),  state%cDLM(i,1:5)
       if(i==1) state%err(err_0) = dummy_r(5)
    enddo
100 state%time%iter = i - 1
    close(ifile)

  end subroutine ReadConvFile

  !> clearing the matrix blocks
  subroutine ClearMatrixBlocks( )
    class(element), pointer :: elem
    integer :: i, j

    do i = 1, grid%nelem
       elem => grid%elem(i)

       elem%block(0)%Mb(:,:) = 0. ! diagonal blocks

       do j = 1, elem%flen         ! off-diagonal blocks
          if(elem%face(neigh,j) > 0)  elem%block(j)%Mb(:,:) = 0.
       enddo
    enddo
  end subroutine ClearMatrixBlocks

  !> clearing the vector blocks
  subroutine ClearVectorBlocks( )
    integer :: i

    do i = 1, grid%nelem
       grid%elem(i)%vec(rhs,:) = 0. ! clearing of rhs terms, e.g., BC
    enddo

  end subroutine ClearVectorBlocks

    !> clearing the matrix blocks elem%blockST
  subroutine ClearMatrixBlocks_STDGM( )
   class(element), pointer :: elem
   integer i,j

   do i = 1, grid%nelem
      elem => grid%elem(i)
      elem%blockST(0)%Mb(:,:) = 0.

      do j = 1, elem%flen
         if (elem%face(neigh,j) > 0) elem%blockST(j)%Mb(:,:) = 0.
      enddo
   enddo
  end subroutine ClearMatrixBlocks_STDGM

    !> clearing the vector blocks elem%rhsST
  subroutine ClearVectorBlocks_STDGM( )
   integer :: i

   do i = 1, grid%nelem
      !F@R change back to 0
      grid%elem(i)%rhsST(1:ndim,:,:) = 0. ! clearing of rhs terms
   enddo

  end subroutine ClearVectorBlocks_STDGM


  !> evaluation of the mass RHS,
  !> extrapolation from the old time levels,
  subroutine SetMassVector( )
    class(element), pointer :: elem
    integer :: ie, i

    do ie = 1, grid%nelem
       elem => grid%elem(ie)

       call SetElementMassVector(elem)

       if (.not. state%time%tdg) then
          if(ndim == 4) then
             elem%w(0,:) = elem%w(1,:)
          else
             elem%w(0,:)  = 0.

             do i=1,state%time%deg
                elem%w(0,:) = elem%w(0,:) + state%time%extrap(i) * elem%w(i,:)
             enddo
          endif
       endif
    enddo
  end subroutine SetMassVector

  !>  evaluation of inviscid and viscous fluxes through element
  !> \f$ \int_{\partial K} \vec{F}(w_h, \nabla w_h) \cdot \vec{n}\, dS \f$
  subroutine ComputeOneElementFluxes(elem, Set_f_s, Set_R_s)
    type(element), intent(inout) :: elem
    interface
      subroutine Set_R_s(ndimL, nbDim, Qdof, w, Dw, Re_1, R_s)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
      end subroutine Set_R_s
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
      end subroutine Set_f_s
    end interface

    integer :: i,ie, j,k, dofA

    elem%vec(rhs,:) = 0.

    !elem%deg_plus = .true.
    elem%deg_plus = .false.

    dofA = elem%dof
    if(elem%deg_plus) dofA = elem%dof_plus
    !dofA = 1

    ! volume terms
    call ElementInvVisVolumes(elem, Set_f_s, Set_R_s, dofA )

    !print*,'computing of edges edges', elem%i,dofA
    do j = 1, elem%flen
       k = elem%face(neigh,j)

       call ElementInvVisFlux(elem, j, Set_f_s, Set_R_s, dofA )
    enddo

    ! adding of the source term
    if(state%RHS_presented .and. &
         ( state%modelName == 'scalar' .or.state%modelName == '2eqs' .or. &
         (state%modelName == 'NSe' .and. state%type_IC .eq. 9) ) ) call ElementRHS(elem )

    !write(200+state%space%adapt%adapt_level,'(30es12.4)') elem%xc(:), abs(elem%vec(rhs,:) )
    !write(*,'(30es12.4)') elem%xc(:), elem%vec(rhs,1:dofA)

    ! ||phi||_{L^2(K)} = 1
    elem%estimS = (dot_product(elem%vec(rhs,1:dofA), elem%vec(rhs,1:dofA) ) &
         / elem%area)**0.5

    !elem%estimS = (dot_product(elem%vec(rhs,1:dofA), elem%vec(rhs,1:dofA) ) )**0.5

    !print*,'Stooped in ComputeOneElementFluxes'
    !stop

    elem%deg_plus = .false.

  end subroutine ComputeOneElementFluxes

  !>  evaluation of all integrals and filling the matrix for ONE element,
  !>  the appropriate parts of
  !> matrix \f$ {\bf C}( {\bf w}) \f$, vectors \f$ {\bf q}( {\bf w}) \f$,
  !> \f$ {\bf m}( {\bf w}) \f$
  subroutine ComputeOneElementTerms(Set_f_s, Set_A_s, Set_Ppm, Set_R_s, Set_K_sk, Set_S, Set_DS, elem)
   interface
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
      end subroutine Set_f_s
      subroutine Set_A_s(ndimL, nbDim, Qdof, w, A_s, xi)
         integer, intent(in) :: Qdof, nbdim, ndimL
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: A_s
         ! matrices A_s in  -- " --
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
      end subroutine
      subroutine Set_Ppm( ndimL, nbDim, Qdof, w, n, xi, Ppm, one_over_area, elem)
         import :: element
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: Ppm
                                               ! matrices Ppm in  -- " --
         real, dimension(1:Qdof, 1:nbDim), intent(in) :: n   ! outer normal
         real, dimension(1:Qdof, 1:nbDim),intent(in) ::  xi                    ! node on the edge?
         real, intent(in), optional :: one_over_area
         type(element), intent(inout), optional :: elem
      end subroutine
      subroutine Set_R_s(ndimL, nbDim, Qdof, w, Dw, Re_1, R_s)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
      end subroutine Set_R_s
      subroutine Set_K_sk(ndimL, nbDim,Qdof, w, Dw, Re_1, K_sk)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,ndimL), intent(inout) :: K_sk
       end subroutine Set_K_sk
       subroutine Set_S(ndimL, nbDim, Qdof, xi, w, Dw, S)
          integer, intent(in) :: ndimL, nbDim, Qdof
          real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
          real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
          real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
          real, dimension(1:Qdof, 1:ndimL), intent(inout) :: S
       end subroutine
       subroutine Set_DS(ndimL, nbDim, Qdof, xi, w, Dw, DS)
          integer, intent(in) :: ndimL, nbDim, Qdof
          real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
          real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
          real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
          real, dimension(1:Qdof, 1:ndimL, 1:ndimL), intent(inout) :: DS
       end subroutine
    end interface
    type(element), intent(inout) :: elem
    class(element), pointer :: elem1
    integer :: ie, j,k, i
    real :: Re_1, Re, val, Mnorm
    integer :: i_min, i_max
    real, dimension(:), allocatable :: x1,x2,x3
    logical :: explicit_convection

    i = elem%i

    !explicit_convection = .true.   ! only for scalar equation (ndim = 1) makes sense
    explicit_convection = .false.

    if( explicit_convection ) then !scalar equation, explicit inviscid terms

       ! print*, 'Scalar linear convection-reaction'
       call LinearInviscidVolumes(elem )

       !print*,'explicit inviscid volumes',i
       call ExplicitInviscidVolumes(elem, Set_f_s_scalar)

       !print*,'explicit of inviscid edges', elem%i
       do j = 1, elem%flen
          k = elem%face(neigh,j)

          if(k > 0) then
             elem1 => grid%elem(k)
             !print*,'Inner edge',j,k,state%time%tau(:)
             call ExplicitInviscidEdge(elem, elem1, j, Set_NumFlux_scalar)

          else
             !print*,'Boundary edge',j,k,state%time%tau(:)
             call ExplicitInviscidEdge(elem, elem, j,  Set_NumFlux_scalar)
          endif

       enddo

    else ! fully (semi-)implicit

       !write(*,'(a6,30es12.4)') 'vec:',elem%vec(rhs,:)
       !write(*,'(a8,i3,50es12.4)') 'vec:',i, elem%vec(rhs,:)

       ! print*,'Matrix ani A'
       !call WriteMblock_Screene(elem%block(0) )       !call WriteMatrixA(0.)
       !write(*,'(a6,200es12.4)') 'wsol:',grid%elem(1)%w(0, :)



       !write(*,'(a35,i5,2es12.4)')'computing of inviscid volumes',elem%i, elem%xc(:)
       call ElementInviscidVolumes(elem, Set_f_s, Set_A_s)
       !write(*,'(a8,i3,50es12.4)') 'vec:',i, elem%vec(rhs,:)

       !call WriteMblock_Screene(elem%block(0) )       !call WriteMatrixA(0.)

       if(state%ST_Vc >0 .and. state%model%Re == 0.) then
          !print*,'computing volume stabilization for inviscid flows'
          call ElementVolumesStabil(elem )
       endif



       !write(*,'(a8,2i3,50es12.4)') 'vec:',i, 0, elem%vec(rhs,:)

       !print*,'computing of inviscid edges', elem%i
       do j = 1, elem%flen
          k = elem%face(neigh,j)
          
          if(k > 0) then
             ! setting of the correct neighbour for the local problem
             if(state%local_problem)  then 
                elem1 => gridL%elem(k)
             else ! standard case
                elem1 => grid%elem(k)
             endif

             !print*,'Inner edge',j,k,elem%i, elem1%i
             call ElementInviscidInnerEdge(elem, elem1, j, Set_Ppm, Set_f_s)
             !write(*,'(a8,i3,50es12.4)') 'vec:',i, elem%vec(rhs,:)

          elseif(elem%iBC(j) == 0 ) then
             !print*,'FIXED WALLS ',elem%i,j,k,elem%iBC(j)
             if(state%modelName == 'scalar' .or. state%modelName == '2eqs') then
                !print*,' scalar equation, pure extrapolation wii = wi'
                call ElementInviscidNeumannlEdge(elem, j, Set_Ppm, Set_f_s)
             else
                ! Euler and/or Navier-Stokes equations
                !if(state%model%Re > 0.) then  !viscous case
                !   call ElementInviscidIOEdge(elem, j, Set_Ppm)
                !else  ! inviscid BC

                ! variant from the pressure linearization
                !call ElementInviscidWallEdge(elem, j)

                ! mirror BC
                call ElementInviscidWall_2_Edge(elem, j, Set_Ppm)
                !endif

             endif
             !write(*,'(a8,i3,50es12.4)') 'vec:',i, elem%vec(rhs,:)

          else
             !print*,'Input/output  edges',elem%i,j,k, 'iBC', elem%iBC(j)
             call ElementInviscidIOEdge(elem, j, Set_Ppm)
             !write(*,'(a8,i3,50es12.4)') 'vec:',i, elem%vec(rhs,:)

          endif

          !if(state%local_problem)  return

          !if(elem%i == 1) then
          !   print*,'EDGE : = ',j
          !   call WriteMblock_Screene(elem%block(0) )
          !endif

          !if(elem%i == 1) write(*,'(a5,50es12.4)') 'q:',  elem%vec(rhs, 1:4)
          !call WriteMblock(elem%block(0) )

          !if(i >= 99 .and. i <= 100) &
          !write(*,'(a6,i5,30es12.4)') 'resIe:',elem%i,elem%vec(rhs, :)

          !write(*,'(a8,2i3,50es12.4)') 'vec:',i, j, elem%vec(rhs,:)

       enddo
    endif



    if(state%model%Re < 0. .or. state%ST_Vc <  -1.E-05 ) then

       !print*,' shock capturing via viscous terms'
       call ElementViscousVolumes(elem,  Set_R_s_Laplace, Set_K_sk_Laplace)
       !call ElementViscousVolumes(elem,  Set_R_s, Set_K_sk)

       do j = 1, elem%flen
          k = elem%face(neigh,j)

          if(k > 0) then
            ! setting of the correct neighbour for the local problem
             if(state%local_problem)  then 
                elem1 => gridL%elem(k)
             else ! standard case
                elem1 => grid%elem(k)
             endif

             ! viscous edges + penalty
             call ElementViscousInnerEdge(elem, elem1, j, Set_R_s_Laplace, Set_K_sk_Laplace)

             ! only penalty
             !call ElementViscousInnerPenalty(elem, elem1, j,Set_R_s_Laplace, Set_K_sk_Laplace)

             !call ElementViscousInnerEdge(elem, elem1, j,  Set_R_s, Set_K_sk)
          else ! boundary edge
             !!if(elem%tBC(j) == 0 )  then  ! inlet or walls
             !print*,' Inlet or solid wall',elem%i,j,elem%iBC(j)
             ! HERE
!SC           call ElementViscousBoundEdge(elem, j, Set_R_s_Laplace, Set_K_sk_Laplace)

          endif
       enddo
    endif


    !write(*,'(a10,i5,100es12.4)') 'ED RHS a', 1, grid%elem(1)%vec(rhs, :)

    !print*,' computing of viscous and penalty terms for viscous flow',elem%i
    if(state%model%Re > 0.) then
       !print*,'viscous volume',elem%i
       call ElementViscousVolumes(elem,  Set_R_s, Set_K_sk)

       !if(i >= 99 .and. i <= 100) &
       !   write(*,'(a6,i5,30es12.4)') 'resVV:',elem%i,elem%vec(rhs, :)

       !print*,'computing of viscous edges', elem%i
       do j = 1, elem%flen
          k = elem%face(neigh,j)

          if(k > 0) then
             ! setting of the correct neighbour for the local problem
             if(state%local_problem)  then 
                elem1 => gridL%elem(k)
             else ! standard case
                elem1 => grid%elem(k)
             endif


             !print*,'Inner edges',elem%i,j,k
             call ElementViscousInnerEdge(elem, elem1, j, Set_R_s, Set_K_sk)

          else ! boundary edge
             ! pocita cleny i pres vystup, testovano na MTC3, nema zadny vliv
             ! konvergence je naprosto stejna
             !!if(elem%tBC(j) == 0 )  then  ! inlet or walls
             !print*,' Inlet or solid wall',elem%i,j,elem%iBC(j)
             call ElementViscousBoundEdge(elem, j, Set_R_s, Set_K_sk)
             !!endif
          endif
          !call WriteMblock(elem%block(0) )
       enddo

    endif


    !print*, 'reaction terms (at this moment only for scalar)'
    if( (state%modelName == 'scalar' .or.state%modelName == '2eqs') .and. state%model%ireac > 0) &
       call ElementReactionVolumes(elem, Set_S, Set_DS)


    !print*,' adding of the source term', elem%i
    if(state%RHS_presented .and. &
         ( state%modelName == 'scalar' .or.state%modelName == '2eqs' .or. &
         (state%modelName == 'NSe' .and. state%type_IC .eq. 9) ) ) call ElementRHS(elem )
         !.or. state%modelName == 'wet_steam'  ) ) call ElementRHS(elem )
    !write(*,'(a8,i3,50es12.4)') 'vec:',i, elem%vec(rhs,:)


  end subroutine ComputeOneElementTerms

  !>  evaluation of all integrals and filling the matrix per elements
  !>  matrix \f$ {\bf C}( {\bf w}) \f$, vectors \f$ {\bf q}( {\bf w}) \f$,
  !> \f$ {\bf m}( {\bf w}) \f$
  subroutine ComputeElementsTerms(Set_f_s, Set_A_s, Set_Ppm, Set_R_s, Set_K_sk, Set_S, Set_DS)
    interface
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
      end subroutine Set_f_s
      subroutine Set_A_s(ndimL, nbDim, Qdof, w, A_s, xi)
         integer, intent(in) :: Qdof, nbdim, ndimL
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: A_s
         ! matrices A_s in  -- " --
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
      end subroutine
      subroutine Set_Ppm( ndimL, nbDim, Qdof, w, n, xi, Ppm, one_over_area, elem)
         import :: element
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: Ppm
                                               ! matrices Ppm in  -- " --
         real, dimension(1:Qdof, 1:nbDim), intent(in) :: n   ! outer normal
         real, dimension(1:Qdof, 1:nbDim),intent(in) ::  xi                    ! node on the edge?
         real, intent(in), optional :: one_over_area
         type(element), intent(inout), optional :: elem
      end subroutine
      subroutine Set_R_s(ndimL, nbDim, Qdof, w, Dw, Re_1, R_s)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
      end subroutine Set_R_s
      subroutine Set_K_sk(ndimL, nbDim,Qdof, w, Dw, Re_1, K_sk)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,ndimL), intent(inout) :: K_sk
       end subroutine Set_K_sk
       subroutine Set_S(ndimL, nbDim, Qdof, xi, w, Dw, S)
          integer, intent(in) :: ndimL, nbDim, Qdof
          real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
          real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
          real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
          real, dimension(1:Qdof, 1:ndimL), intent(inout) :: S
       end subroutine
       subroutine Set_DS(ndimL, nbDim, Qdof, xi, w, Dw, DS)
          integer, intent(in) :: ndimL, nbDim, Qdof
          real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
          real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
          real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
          real, dimension(1:Qdof, 1:ndimL, 1:ndimL), intent(inout) :: DS
       end subroutine
    end interface
    class(element), pointer :: elem, elem1
    integer :: i,ie, j,k
    real :: Re_1, Re, val, Mnorm
    integer :: i_min, i_max
    real, dimension(:), allocatable :: x1,x2,x3

    do i = 1, grid%nelem
       elem => grid%elem(i)

       !if(state%modelName == 'NSe') call CheckPhysicalProperties(elem)

       !print*,'computing of jump indicator', i
       elem%rezid = 0.

       if(state%model%Re < 0.) call ElementGradIndicator(elem)

       !if(elem%xc(1) > 0.2) then
       if(state%ST_Vc >0 .or.  state%ST_Ec> 0. .or. state%ST_Vc <= -1.E-05 ) then
          do j = 1, elem%flen
             k = elem%face(neigh,j)
             if(k > 0) then
                !elem1 => grid%elem(k)
                call ElementEdgeJumpIndicator(elem,  j)
                !call ElementEdgeDerJumpIndicator(elem, elem1, j)
             endif
          enddo
          Re = elem%rezid
          elem%rezid = elem%rezid * elem%limit_par
          !!elem%rezid = (elem%rezid ) **0.5

          val = 0.
          if(elem%rezid > 1.) then
             val = elem%rezid
             val = val**0.5

             if(state%modelName == 'scalar' .or.state%modelName == '2eqs') val = 1.
             !val = 1. ! used for MTC2 not M=2.00 ????

             if(state%ST_Vc < -1.E-05) val = 1.

          else
             !elseif(elem%rezid > 0.4) then
             !val = (sin((elem%rezid -0.85)*3.14159267/2/0.15) + 1)/2 ! BAD
             !val = (sin((elem%rezid -0.8)*3.14159267/2/0.2) + 1)/2
             !val = (sin((elem%rezid -0.7)*3.14159267/2/0.3) + 1)/2
             !val = (sin((elem%rezid -0.6)*3.14159267/2/0.4) + 1)/2
             val = (sin((elem%rezid -0.5)*3.14159267/2/0.5) + 1)/2
          endif

          !write(1000+state%time%iter,*) elem%xc(:), val, state%nlSolver%iter, Re, elem%rezid

          !if(abs(elem%xc(1) - 0.2) <= 0.1 .and. abs(elem%xc(2) - 0.35) <= 0.1 )&
          !     write(600,'(2es12.4,2i5,4es12.4)') &
          !     elem%xc(:),elem%i, state%nlSolver%iter, state%time%iter + 1.*state%nlSolver%iter/30, &
          !     Re, elem%rezid, val

          !if(elem%rezid> 1e-9) write(100,*) elem%xc(1:nbDim), elem%rezid

          elem%rezid = val

          !elem%rezid = 0.

          ! new version
          if(state%ST_Vc < -1.E-05) elem%rezid = val * elem%diam**2 * abs(state%ST_Vc )


       endif

       ! this cycle can not be connected with the following ones thanks the jump indic
    enddo

    do i=1,2
       call SmoothResid()
    enddo

    !print *,'Residuum for stabilization computed'

    ! new version
    if(state%ST_Vc < -1.E-05) then
       do i=1,2
          call SmoothResid()
       enddo

       call SetArtificialViscosity( )
    endif


    if(state%model%Re < 0.) then
       !print*,'Smoothing procedure'
       !call EquivResid()
       !call EquivResid1()

       do i=1,3
          call SmoothResid()
       enddo

       call SetArtificialViscosity( )
    endif


    if(state%modelName == 'NSe' .or. state%modelName == 'wet_steam') call EvalOutputPressure()

!    call omp_set_num_threads( 4 )
!! !$OMP PARALLEL DO !I is private by default
    do i = 1, grid%nelem
       elem => grid%elem(i)

       !if( elem%rezid > 1E-01) &
       !if( elem%i == 99) &
       !     write( *,'(a4, i5, 5es14.6)') '??',elem%i, elem%xc(:), elem%rezid

       !write(*,*) '------------------------------'  ! wet steam
       !write(*,*) 'elem(i)', i
       !write(*,*) 'w', elem%w(0,1:ndim)
       !write(*,*) "-----------------------------------------------------"
       !write(*,*) "ComputeElementsTerms calls ComputeOneElementTerms for element number", i  ! wet steam
       call ComputeOneElementTerms(Set_f_s, Set_A_s, Set_Ppm, Set_R_s, Set_K_sk, &
            Set_S, Set_DS, elem )

       !print* , 'after block:', elem%block(0)%Mb(1,1)

    enddo
!!  !$OMP END PARALLEL DO

!    close(state%time%iter+100)

    !print*,'ComputeElementsTerms - done'
  end subroutine ComputeElementsTerms

  subroutine SaveMatrixBlocks(elem)
    type(element), intent(in) :: elem
    integer :: ifile, i,j,k

    ifile = 12

    open(ifile, file='blocks', status='unknown', position = 'append')
    write(ifile,'(3i6,2es12.4)') elem%i,elem%flen,count(elem%face(neigh,:) >0),elem%xc(:)
    write(ifile, '(a5,40es12.4)') 'rhs:',elem%vec(rhs,:)
    !write(ifile, *)
    write(ifile,'(3i6)') elem%i,0,0
    do k=1,elem%dof*ndim
       write(ifile, '(i5,40es12.4)') k,elem%block(0)%Mb(k,:)
    enddo

    do j=1,elem%flen
       if(elem%face(neigh,j) >0) then
          write(ifile, *)
          write(ifile,'(3i6,2es12.4)') elem%i,j,elem%face(neigh,j),&
               grid%x(elem%face(idx,j),1:nbDim)

          do k=1,elem%dof*ndim
             write(ifile, '(i5,40es12.4)') k,elem%block(j)%Mb(k,:)
          enddo
       endif
    enddo
    !write(ifile,*)
    !write(ifile,*)

    close(ifile)
  end subroutine SaveMatrixBlocks

  !> Setting of the right-hand side and init solution into a global vector
  subroutine FillVectorsOLD(b, x, etaL)
    real, dimension(:), intent(inout) :: b, x
    real,intent(in):: etaL
    class(element), pointer:: elem ! one element
    integer  :: ie, j, k, Row, ivec, ndof

    ivec = 1

    do ie=1,grid%nelem
      elem => grid%elem(ie)
      ndof = elem%dof * ndim
      x(ivec:ivec+ndof-1) = elem%w(0,1:ndof)

      if (etaL == 0.) then
         b(ivec:ivec+ndof-1) = elem%vec(rhs,1:ndof)
      else
         b(ivec:ivec+ndof-1) = elem%vec(rhs,1:ndof) + etaL * elem%vec(rhsM,1:ndof)
      endif
      ivec = ivec + ndof
   end do

   !return
   !! normalization:  b := b/|K_i|,   A := A/|K_i|
   !ivec = 1

   !do ie=1,grid%nelem
   !   elem => grid%elem(ie)
   !   ndof = elem%dof * ndim
   !
   !   b(ivec:ivec+ndof-1) = b(ivec:ivec+ndof-1)  / elem%area
   !   ivec = ivec + ndof
   !
   !   elem%block(0)%Mb(:,:) = elem%block(0)%Mb(:,:)  / elem%area
   !   do j=1,elem%flen
   !      if(elem%face(neigh,j) > 0) &
   !           elem%block(j)%Mb(:,:) = elem%block(j)%Mb(:,:)  / elem%area
   !   enddo
   !end do

 end subroutine FillVectorsOLD

 !> Setting of the right-hand side and init solution into a global vector
 subroutine FillVector(b, etaL)
   real, dimension(:), intent(inout) :: b
   real, intent(in) :: etaL
   class(element), pointer:: elem ! one element
   integer  :: ie, ivec, ndof

   ivec = 1

   do ie=1,grid%nelem
      elem => grid%elem(ie)
      ndof = elem%dof * ndim

      !write(*,'(a6,20es12.4)') 'CN FV', elem%vec(rhs,1:ndof)
      !write(*,'(a6,20es12.4)') 'CN FV', elem%vec(rhsOLD,1:ndof)

      if( etaL > 0) then
         ! with mass matrix terms
         call SetElementMassVectorS(elem, b(ivec:ivec+ndof-1))

         b(ivec:ivec+ndof-1) = - state%time%alpha(0)/state%time%tau(1)* b(ivec:ivec+ndof-1)

         if(state%time%cn) then
            ! Crank Nicolson
            b(ivec:ivec+ndof-1) = 2. *b(ivec:ivec+ndof-1) &
                 + elem%vec(rhs,1:ndof) +   elem%vec(rhsOLD,1:ndof)    &
                 + 2. * elem%vec(rhsM,1:ndof) / state%time%tau(1)
         else

            b(ivec:ivec+ndof-1) = b(ivec:ivec+ndof-1) &
                 + elem%vec(rhs,1:ndof) +   elem%vec(rhsM,1:ndof) / state%time%tau(1)
         endif

      else
         ! without mass matrix terms
         b(ivec:ivec+ndof-1) = elem%vec(rhs,1:ndof)
      endif

      ivec = ivec + ndof
   end do
 end subroutine FillVector

 !> Setting of the right-hand side and init solution into a global vector for STDG method
 subroutine FillVectorST(b, etaL)
   real, dimension(:), intent(inout) :: b
   real, intent(in) :: etaL
   class(element), pointer:: elem ! one element
   integer  :: ie, ivec, dof, ndof, l, k, kvec

   ivec = 0

   !write(59,'(a6,100es12.4)' ) '#bb#',b(:)
   do ie=1,grid%nelem
      elem => grid%elem(ie)
      dof = elem%dof
      ndof = dof * elem%Tdof * ndim

      b(ivec+1:ivec+ndof) = 0 !elem%rhsST()
      kvec = ivec

      do l = 1, elem%Tdof
         do k = 1, ndim

            b(kvec+1:kvec + dof) =  elem%rhsST(k,1:elem%dof, l)

            !write(59, '(5i5, 200es12.4)')elem%i, l , k , kvec+1, kvec+dof, b(kvec+1:kvec + dof)

            kvec = kvec + dof
         enddo !k
      enddo !l

      ivec = ivec + ndof
   end do
   !write(59,'(a6,100es12.4)' ) '#bb#',b(:)
   !write(59,*)'----------------',state%nsize,dof, ndof
 end subroutine FillVectorST


 !> evaluation of \f$ F(w^{k-1}) \f$ for the Crank-Nicolson,
 !> \f$ w^{k-1} \f$ is the initial approximation for \f$ w^{k} \f$ hence we use it
 subroutine PrepareCrankNicolson( )
   class(element), pointer:: elem ! one element
   logical ::  loc_implicitly
   real :: loc_ctime
   integer  :: ie, ivec, ndof

   loc_implicitly = state%nlSolver%implicitly
   loc_ctime = state%time%ctime

   state%nlSolver%implicitly = .false.
   state%time%ctime = state%time%ttime   ! time level $t_{k-1}$ !!

   !print*,'Storing of vec(rhsOLD,*)', state%nlSolver%implicitly, state%time%ctime

   call ComputeTerms( )

   do ie=1,grid%nelem
      elem => grid%elem(ie)
      ndof = elem%dof * ndim
      ! Crank-Nicolson, we store the vector F(w^k) = F(w^{k-1})
      elem%vec(rhsOLD,1:ndof)  = elem%vec(rhs,1:ndof)

   end do

   state%nlSolver%implicitly = loc_implicitly
   state%time%ctime = loc_ctime

 end subroutine PrepareCrankNicolson

 !> evaluation of the matrix and vector blocks
 subroutine ComputeTerms( )
   integer :: i

   ! clearing of the appropriate arrays
   if(state%nlSolver%implicitly) call ClearMatrixBlocks()
   call ClearVectorBlocks()

   ! setting of the corresponding fluxes
   if(state%modelName == 'scalar') then        ! 2D scalar equation
      call ComputeElementsTerms(Set_f_s_scalar, Set_A_s_scalar, Set_Ppm_scalar, &
           Set_R_s_scalar, Set_K_sk_scalar, Set_S_scalar, Set_DS_scalar)

   elseif(state%modelName == 'NSe' ) then    ! 2D Euler and Navier-Stokes equations
      call ComputeElementsTerms(Set_f_s_Euler, Set_A_s_Euler, Set_Ppm_Euler, &
           Set_R_s_NS, Set_K_sk_NS, Set_S_empty, Set_DS_empty)

   elseif(state%modelName == '2eqs') then    ! 2 scalar equations
      call ComputeElementsTerms(Set_f_s_scalar, Set_A_s_scalar, Set_Ppm_scalar, &
           Set_R_s_scalar, Set_K_sk_scalar, Set_S_scalar, Set_DS_scalar)

   elseif(state%modelName == 'wet_steam' ) then    ! 2D wet_stem
      call ComputeElementsTerms(Set_f_s_WS, Set_A_s_WS, Set_Ppm_WS, &
           Set_R_s_WS, Set_K_sk_WS, Set_S_empty, Set_DS_empty)


  ! elseif(nbDim == 2 .and. ndim == 6) then    ! RANS - k-omega  model
  !    call ComputeElementsTerms(Set_f_s_Turb2e, Set_A_s_Turb2e, Set_Ppm_Turb2e, &
  !         Set_R_s_Turb2e, Set_K_sk_Turb2e, Set_S_empty, Set_DS_empty)

  ! elseif(nbDim == 3 .and. ndim == 5) then    ! 3D Euler and Navier-Stokes equations
  !    call ComputeElementsTerms(Set_f_s_Euler3D, Set_A_s_Euler3D, Set_Ppm_Euler3D, &
  !         Set_R_s_NS3D, Set_K_sk_NS3D, Set_S_empty, Set_DS_empty )


   else
      print*,'Compute Elements Terms not implemented for ndim=',ndim
   endif

   !do i=1,grid%nelem
   !   !write(*,'(a4,i5,20es12.4)') 'iw: ',i, grid%elem(i)%w(0,:)
   !   !print*, 'EEDE',size (grid%elem(i)%block(0)%Mb(:,:), 1)
   !   call WriteMblock(grid%elem(i)%block(0) )
   !   write(*,'(a10,i5,100es12.4)') 'ED RHS', i, grid%elem(i)%vec(rhs, :)
   !enddo

!  stop 'stopped in the end of ComputeTerms'

 end subroutine ComputeTerms

 !> if block flux matrix \f$ {\bf C}({\bf w}) \f$ and vector \f$ {\bf q}({\bf w}) \f$
 !> are given from implicit variant of ComputeTerms, we put
 !> \f$ {\bf q}({\bf w}) :=  {\bf q}({\bf w})  - {\bf C}({\bf w}){{\bf w}} \f$
 !> instead of explicit variant of ComputeTerms
 subroutine SetF_q_Cw_fast( )
    class(element), pointer:: elem,elem1 ! one element
    integer :: i, j, k, ndof, ndof1

    real, dimension(:),allocatable:: accum

    ! allocate accum once to accomodate for the largest dof.
    allocate(accum(maxval(grid%elem%dof) * ndim))

    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof = elem%dof  * ndim

       ! diagonal block
       accum(1:ndof) = matmul(elem%block(0)%Mb(1:ndof, 1:ndof), elem%w(0,1:ndof) )

       !! off-diagonal blocks
       do j=1,elem%flen
          k = grid%elem(i)%face(neigh,j)

          if(k > 0) then
             elem1 => grid%elem(k)
             ndof1 = elem1%dof * ndim

             accum(1:ndof) = accum(1:ndof)  &
                  + matmul(elem%block(j)%Mb(1:ndof, 1:ndof1), elem1%w(0,1:ndof1) )
          endif
       enddo

       elem%vec(rhs,1:ndof) = elem%vec(rhs,1:ndof) - accum(1:ndof)

    enddo

    deallocate(accum)

 end subroutine SetF_q_Cw_fast

 !> evaluate \f$ \| {\bf w}\| \f$
 function L2_norm( )
   real :: L2_norm
   class(element), pointer :: elem
   integer :: i,j,k

   L2_norm = 0.
   k = 0
   do i=1,grid%nelem
      elem => grid%elem(i)
      j = elem%dof*ndim
      L2_norm =  L2_norm + dot_product(elem%w(0,1:j), elem%w(0,1:j))
      k = k+j
   enddo
   L2_norm = max(L2_norm, 1.)
 end function L2_norm

 subroutine CheckResiduum( )
   class(element), pointer :: elem
   real, dimension(:), allocatable :: accum
   integer:: i, l, ndof, ifile1, ifile2, k, k1, k2

   allocate(accum(1:state%space%max_dof) )

   state%nlSolver%implicitly = .false.
   call ComputeTerms( )
   state%nlSolver%implicitly = .true.

   ifile1 = 20 + state%time%iter
   ifile2 = 30 + state%time%iter

   write(ifile1,*) '# File generated by CheckResiduum in euler.f90'
   write(ifile2,*) '# File generated by CheckResiduum in euler.f90'

   do i=1,grid%nelem
      elem => grid%elem(i)
      ndof = elem%dof*ndim
      accum(1:ndof) = elem%vec(rhs,1:ndof)

      do k=1,ndim
         k1 = (k-1)*elem%dof + 1
         k2 = k*elem%dof

         !print*,'...',k1,k2,state%time%tau(1)
         !print*,'???',accum(k1:k2)
         !print*,'mmm',(elem%w(0,k1:k2) - elem%w(1,k1:k2))/state%time%tau(1)
         !print*,'###',matmul(elem%Mass%Mb(1:elem%dof,1:elem%dof), &
         !     (elem%w(0,k1:k2) - elem%w(1,k1:k2))/state%time%tau(1) )

         accum(k1:k2) = accum(k1:k2)  &
              - matmul(elem%Mass%Mb(1:elem%dof,1:elem%dof), &
              (elem%w(0,k1:k2) - elem%w(1,k1:k2))/state%time%tau(1) )

         !print*,'!!!',accum(k1:k2)
         !print*

      enddo

      do l=1,ndof
         write(ifile1,'(20es14.6)') elem%xc(1:2), accum(l), elem%vec(rhs,l), &
              matmul(elem%Mass%Mb(1:elem%dof,1:elem%dof), &
              (elem%w(0,k1:k2) - elem%w(1,k1:k2))/state%time%tau(1) )
      enddo

      write(ifile2,*) elem%xc(1:2), dot_product(accum(1:ndof), accum(1:ndof) ), &
           dot_product(elem%vec(rhs,1:ndof), elem%vec(rhs,1:ndof))

   enddo

   print*,'# subroutine CheckResiduum terminated'
 end subroutine CheckResiduum

 !> evaluation of the matrix and vector blocks for STDGM method
 subroutine ComputeSTDGM_Terms( )
   class(element), pointer :: elem
   class(Time_rule), pointer :: T_rule
   !real, dimension(:,:), pointer :: phi
   integer :: i, j, k, l, r, kk
   integer :: alpha, Qdeg, dof, Tdof, s_dim, f_dof, f_Tdof, f_dim
   integer :: wTdof, wdof  !NEW for adaptation - smaller than dof/Tdof if elem%deg_plus=TRUE
   real :: cTime
   real :: val, val_vec
   real :: local_eta
   logical ::iprint
   integer :: dofA

!   print*, 'ComputeSTDGM_Terms'



   cTime = state%time%ctime

   local_eta = 1. / state%time%tau(1)


!   !F@R control phi(Tdof + 1) ma koreny v int uzlech

!     OLD
!   if (grid%elem(1)%deg_plus) then
!      Qdeg = state%time%max_Tdof + 1
!   else
!      Qdeg = state%time%max_Tdof + 1
!   endif

   !NEW
   Qdeg = state%time%Qnum

   T_rule => state%time%T_rule(Qdeg)

   associate ( time => state%time)
   select type ( time )
      class is ( TimeTDG_t )

      ! rhsST = q(w)
      if(state%nlSolver%implicitly) then
         !print*, 'ComputeSTDGM_Terms with implicitly = TRUE called'
         call ClearMatrixBlocks_STDGM()
         call ClearVectorBlocks_STDGM()

         do alpha = 1, Qdeg ! temporarily max_Tdof =  max time quadrature nodes

            do i = 1, grid%nelem
               elem => grid%elem(i)

               !	G_rule => state%space%G_rule(elem%TQnum)
               if (Qdeg /= elem%TQnum) then
               !F@R Verify if it is OK, some nodes could be in wrong position
                !  print*, 'Problem in ComputeSTDGM_Terms 1'
                !  print*,'Verify !!!'
               endif
               !call Transfer_wST_to_w_Elem(elem , alpha, elem%TQnum)
               ! save the wST space-time solution in quadrature index alpha to w
               call Transfer_wST_to_w_Elem(elem , alpha, Qdeg)
               !print*, 'still running', elem%w

            enddo

            !we have to run ComputeTerms() in the time-quadrature nodes
            !cTime = state%time%ctime
            state%time%ctime = state%time%ttime + state%time%tau(1) * T_rule%lambda(alpha)
            !   print*, 'ctime:' , state%time%ctime, cTime
            !   print*, '--------------','ctime:', state%time%ctime, state%nlSolver%implicitly
            call ComputeTerms( )
            state%time%ctime = cTime

            !print*, 'vec = ' ,  	grid%elem(2)%vec(rhs, 1 : grid%elem(2)%dof)

            do i =1, grid%nelem
               elem => grid%elem(i)
               Tdof = elem%Tdof
               dof = elem%dof
               s_dim = ndim *dof


               do l = 1, Tdof
                  !diag blocks of blockST
                  do r = 1, Tdof

                     val = T_rule%phi(l,alpha)*T_rule%phi(r,alpha) * T_rule%weights(alpha)  !phi_l(t_alpha) * phi_r(t_alpha)* w_alpha

                     elem%blockST(0)%Mb( (l-1)*s_dim + 1 : l*s_dim, (r-1)*s_dim + 1 : r*s_dim) = &
                          elem%blockST(0)%Mb((l-1)*s_dim + 1 : l*s_dim, (r-1)*s_dim + 1 : r*s_dim) &
                          + elem%block(0)%Mb(1:s_dim, 1:s_dim) * val
                  enddo !r

                  !offdiagonal blocks of blockST
                  do j = 1, elem%flen
                     if(elem%face(neigh,j) > 0) then
                        f_dof = elem%face(fdof,j)
                        f_Tdof = elem%face(fTdof,j)
                        f_dim = ndim * f_dof
                        do r = 1, f_Tdof
                           val = T_rule%phi(l,alpha)*T_rule%phi(r,alpha) * T_rule%weights(alpha)  !phi_l(t_alpha) * phi_r(t_alpha)* w_alpha
                           elem%blockST(j)%Mb((l-1)*s_dim + 1 : l*s_dim , (r-1)*f_dim + 1 : r*f_dim) =  &
                              elem%blockST(j)%Mb((l-1)*s_dim + 1 : l*s_dim , (r-1)*f_dim + 1 : r*f_dim)  &
                              + elem%block(j)%Mb(1:s_dim, 1:f_dim) * val
                        enddo !r=1,f_Tdof
                     endif
                  enddo !j


                  !vec -> rhsST
                  val_vec = T_rule%phi(l,alpha) * T_rule%weights(alpha) ! for now G_rule%weights used for time quadrature

                 ! write(*,'(a8,12es12.4)') 'elem%vec', elem%vec(rhs,1:3)

                  do k = 1, ndim
                     elem%rhsST(k, 1:dof, l) = elem%rhsST(k, 1:dof, l)  &
                        + ( elem%vec(rhs,(k-1)*dof + 1 : k*dof) * val_vec )
                  enddo !k

               enddo !l = 1,Tdof
            enddo !i

         end do !alpha =  1,Qdeg

         !adding the timejump part
         !TO ensure that in wSTfin is really the value from previous time interval
         do i = 1, grid%nelem
            call Elem_wSTfinToRhsST( grid%elem(i) , state%nlSolver%implicitly)
         enddo !i




      !implicitly = .false. => rshST = F(w)
      else

         !print*, 'ComputeSTDGM_Terms with implicitly = FALSE called'
         iprint = .false.
        ! if (grid%elem(1)%deg_plus) iprint = .true.

         if (iprint) print*, 'ComputeSTDGM_Terms with implicitly = FALSE called', grid%elem(1)%deg_plus
         if (iprint) print*, 'size of vec ' , size(grid%elem(1)%vec(rhs,:))

         call ClearVectorBlocks_STDGM()  ! rhsST = 0

         if  (iprint) then

            do l = 1,grid%elem(1)%Tdof_plus
               print*, 'rhsST', 'Tdof:',l, grid%elem(1)%rhsST(1, 1:grid%elem(1)%dof_plus, l)
            enddo
            print*, '____________'
         endif

         cTime = state%time%ctime  ! storing of the actual value

         do alpha = 1, Qdeg ! temporarily max_Tdof =  max time quadrature nodes

            do i = 1, grid%nelem
               elem => grid%elem(i)
               !F@R control
               call Transfer_wST_to_w_Elem(elem , alpha, Qdeg)
               !call Transfer_wST_to_w_Elem(elem , alpha, elem%TQnum)
            enddo

            !we have to run ComputeTerms() in the time-quadrature nodes
            state%time%ctime = state%time%ttime + state%time%tau(1) * state%time%T_rule(Qdeg)%lambda(alpha)
            !print*, 'cTime=', state%cTime
            !     print*
            !      print*, '-----------------------','ctime:', state%time%ctime, state%nlSolver%implicitly
            !      write(*,'(a6,12es12.4)') 'w(t):',grid%elem(1)%w(0,:)*2**0.5
            call ComputeTerms( )

            !if (iprint) print*, 'vec  rhs' , grid%elem(1)%vec(rhs,:)

            do i =1, grid%nelem
               elem => grid%elem(i)


               dofA = elem%dof
               if(elem%deg_plus) dofA = elem%dof_plus

               !NEW for adaptation
               if(elem%deg_plus) then
                  dof = elem%dof_plus
                  Tdof = elem%Tdof_plus

               else
                  dof = elem%dof
                  Tdof = elem%Tdof
               endif
   !kontrola
          !     do k  = elem%dof + 1, dof
          !        write*,
          !     end

               s_dim = ndim * dof

               !         write(*,'(a8,12es12.4)') 'Elem%vec', elem%vec(rhs,1:3)

               !vec -> rhsST

               do k = 1,ndim
                  do l = 1, Tdof
                     val_vec = T_rule%phi(l,alpha) * T_rule%weights(alpha)
                     elem%rhsST(k, 1:dof, l) = elem%rhsST(k, 1:dof, l)  &
                          +  val_vec *  elem%vec( rhs, (k-1)*dof + 1 : k*dof )   	! dofA
                     ! !!! EXP +  val_vec *  elem%vec( rhs, (k-1)*dof + 1 : k*dof )

                     if (i == 1 .and. iprint) then
                        print*, 'alpha' , alpha
                        print*,'val_vec in STDGM_Terms',  val_vec
                        print*, elem%vec( rhs, (k-1)*dof + 1 : k*dof )
                        print*, 'rhsST', 'Tdof:',l, elem%rhsST(k, 1:dof, l)
                        print*
                     endif
                  enddo !l

               enddo !k

            enddo !i
         end do !alpha =  1,Qdeg

         ! putting back the original value
         !state%time%ctime = cTime

         !    write(*,'(a15, 12es12.4)') 'wST:' , grid%elem(1)%wST(1,:,:)
         !   write(*,'(a15, 12es12.4)') 'F^m: aft vec ' , grid%elem(1)%rhsST(1,:,:)

         do i =1, grid%nelem
            elem => grid%elem(i)

            !phi => state%space%V_rule(elem%Qnum)%phi

            !NEW for adaptation
            wdof = elem%dof
            wTdof = elem%Tdof

               !NEW for adaptation
            if(elem%deg_plus) then
              dof = elem%dof_plus
              Tdof = elem%Tdof_plus

            else
              dof = elem%dof
              Tdof = elem%Tdof
            endif


            !Tdof = elem%Tdof
            !dof = elem%dof

          !  iprint=.false.

            !if( state%time%iter >=118 .and. state%space%adapt%adapt_level >=33 .and. &
            !     elem%i >=650 .and. elem%i <= 650 ) &

         !   if(elem%deg_plus)   iprint =.true.

            do k = 1, ndim

               ! if(iprint) then
               !    write(*,'(a6,5i5)'),'!@@@!',elem%i, dof, Tdof, wdof, wTdof
               !    !write(*,'(a4,i5,30es12.4)')'wST',j,  elem%wST(k,1:wdof,1:wTdof)
               !    do l=1,Tdof
               !       write(*,'(a4,i5,30es12.4)')'refTime',l, state%time%refTimeMatrix%Mb(l,:)
               !    enddo
               !    do j=1,dof
               !       write(*,'(a4,i5,30es12.4)')'MASS',j, elem%Mass%Mb(j, 1:wdof)
               !    enddo
               !    print*,'### Trouble HERE'

               !    print*, 'Size of Mass:' , size(elem%Mass%Mb(:,1)) ,'x', size(elem%Mass%Mb(1,:))
               !    stop
               ! endif

               do l = 1,Tdof
                  do j = 1,dof   ! dofA
                     !NEW for adaptation

                     val_vec = local_eta * dot_product( time%refTimeMatrix%Mb(l,1:wTdof) , &
                          matmul( elem%Mass%Mb(j, 1:wdof), elem%wST(k,1:wdof,1:wTdof) ))

                   !  if(iprint) print*, l,'from mass', val_vec


                     elem%rhsST(k, j, l) = elem%rhsST(k, j, l) - val_vec
                  enddo !j
               enddo !l

            enddo !k
         enddo !i


         !  write(*,'(a15, 12es12.4)') 'F^m: aft mass' , grid%elem(1)%rhsST(1,:,:)


         ! the w_{m-1}^- part:
         do i = 1, grid%nelem
           call Elem_wSTfinToRhsST( grid%elem(i) , state%nlSolver%implicitly)
         enddo !i

      endif   !state%nlSolver%implicitly
      class default
      stop 'For TimeTDG_t only!!!'
   end select
   end associate ! time

   ! putting back the original value
   state%time%ctime = cTime

 !  print*, '________________________________________'
 !  print*, 'End of ComputeSTDGM_Terms with implicitly = ', state%nlSolver%implicitly


 end subroutine ComputeSTDGM_Terms


end module euler_problem




