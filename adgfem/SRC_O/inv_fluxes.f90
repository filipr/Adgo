!> subroutine of evaluation of Euler and Navier-Stokes fluxes,
!> the main iterative loop
module inviscid_fluxes
  use main_data
  use problem_oper
!  use matrix_oper_int
  use set_solution
  use stdgm_mod
  use mesh_mod
!  use define_state
!  use wet_steam_paramets
!  use porous_fnc

  implicit none

  public:: InitElementW
  public:: ElementInviscidVolumes
  public:: ElementInviscidInnerEdge
  public:: ElementInviscidIOEdge
  public:: ElementInviscidWallEdge
  public:: ElementInviscidWall_2_Edge
  public:: ElementVolumesStabil

  public:: ElementViscousVolumes
  public:: ElementViscousInnerEdge
  public:: ElementViscousBoundEdge

  public:: ElementInvVisFlux

  public :: ElementRHS
  public :: ElementSubdomainRHS

  !public:: AddElementMassMatrix
  public:: SetElementMassVector
  public:: SetElementMassVectorS

  public:: DirectComputeDrag_Lift

  public :: ExplicitInviscidEdge
  public :: ExplicitInviscidVolumes
  public :: LinearInviscidVolumes

  public:: ComputeSkinFriction
  public:: ElementReactionVolumes

  public :: Elem_wSTfinToRhsST

contains




  !> evaluation of inviscid edge integrals in solid walls
  subroutine ElementInviscidWallEdge(elem, ie)
    type(element), intent(inout):: elem          ! elem = element
    integer, intent (in) :: ie                   ! inner index of edge, 1,2,3, (4)
    real, dimension(:,:), allocatable :: wi        ! w in integ nodes
    real, dimension(:,:), allocatable :: nc    ! normals in integ nodes
    real, dimension(:,:), allocatable :: f_s     ! explicit discretization
    real, dimension(:), allocatable :: ident     ! explicit discretization
    real, dimension(:,:,:), allocatable :: Ppm   ! matrices Ppm in integ nodes
    integer ::  dof, dofA, Qdof, kst
    integer :: i, k, k1, row, col

    ! stop for wet steam case
    if(state%modelName == "wet_steam") then
      write(*,*) "ElementInviscidWallEdge is not suitable for wet steam"
      stop
    end if

    dof = elem%dof
    dofA = dof
    if(elem%deg_plus) dofA = elem%dof_plus

    Qdof = elem%face(fGdof,ie)    ! "num" integration nodes

    allocate(wi(1:Qdof,1:ndim))

    call Eval_w_Edge(elem, ie, wi, .false.)

    ! setting of outer normals in integration nodes
    allocate(nc(1:Qdof, 1:nbDim) )

    if(elem%ibcur > 0) then
       nc(1:Qdof,1:nbDim) = elem%nc(1:Qdof,1:nbDim)
    else
       nc(1:Qdof,1) = elem%n(ie,1)
       nc(1:Qdof,2) = elem%n(ie,2)
    endif


    ! evaluation of matrix
    allocate(Ppm(1:Qdof, 1:ndim, 1:ndim) )
    if(ndim == 4 .and. nbDim == 2) then
       call Set_Ppm_Euler_Slip(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), nc(1:Qdof,1:nbDim), &
            Ppm(1:Qdof, 2:3, 1:ndim))
    elseif(ndim > 4 .and. ndim <= 6 .and. nbDim == 2) then
       call Set_Ppm_Turb2e_Slip(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), nc(1:Qdof,1:nbDim), &
            Ppm(1:Qdof, 2:3, 1:ndim))
    endif

    if(state%nlSolver%implicitly) then
       do k=2,3                   ! k = index of component of w, only second and third !!
          row = (k-1)*dof

          do k1=1,ndim          ! k1 = index of component of w
             col = (k1-1)*dof

             call IntegrateEdgeBlockBB(elem, ie, elem, ie, Ppm(1:Qdof, k, k1), &
                  elem%block(0)%Mb(row+1:row+dof, col+1:col+dof) )
          enddo !k1
       enddo ! k

    else ! explicit discretization

       allocate(f_s(1:Qdof,1:ndim), ident(1:Qdof))
       ident(:) = 1.

       Ppm(1:Qdof, 1, 1:ndim) = 0.
       Ppm(1:Qdof, 4:ndim, 1:ndim) = 0.

       do i=1,Qdof
          f_s(i,1:ndim) = matmul(Ppm(i, 1:ndim, 1:ndim), wi(i, 1:ndim) )

       enddo

       call EvalEdgeVectorB(elem, ie, ident(1:Qdof), -f_s(1:Qdof,1:ndim), &
            dofA, elem%vec(rhs, 1:ndim*dofA ) )

    endif


    deallocate(Ppm )
    deallocate(nc)
    deallocate(wi)

  end subroutine ElementInviscidWallEdge



  !> evaluation of inviscid edge integrals for Neumann BC for scalar equation
  subroutine ElementInviscidNeumannlEdge(elem, ie, Set_Ppm, Set_f_s)
    type(element), intent(inout):: elem          ! elem = element
    integer, intent (in) :: ie                   ! inner index of edge, 1,2,3, (4)
    interface
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
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x, ie )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
         integer, intent(in) :: ie
      end subroutine Set_f_s
    end interface

    real, dimension(:,:), allocatable :: wi        ! w in integ nodes
    real, dimension(:,:), allocatable :: nc    ! normals in integ nodes
    real, dimension(:,:), allocatable :: f_s     ! explicit discretization
    real, dimension(:), allocatable :: ident     ! explicit discretization
    real, dimension(:,:,:,:), allocatable :: Ppm   ! matrices Ppm in integ nodes
    integer ::  dof, dofA, Qdof, kst
    integer :: i


    dof = elem%dof
    dofA = dof
    if(elem%deg_plus) dofA = elem%dof_plus

    Qdof = elem%face(fGdof,ie)    ! "num" integration nodes

    allocate(wi(1:Qdof,1:ndim))

    call Eval_w_Edge(elem, ie, wi, .false.)

    !write(33,*) elem%xc(:)

    ! setting of outer normals in integration nodes
    allocate(nc(1:Qdof, 1:nbDim) )

    if(elem%ibcur > 0) then
       nc(1:Qdof,1:nbDim) = elem%nc(1:Qdof,1:nbDim)
    else
       nc(1:Qdof,1) = elem%n(ie,1)
       nc(1:Qdof,2) = elem%n(ie,2)
    endif


    ! evaluation of matrix
    allocate(Ppm(1:Qdof, 1:nbDim, 1:ndim, 1:ndim) )
    call Set_Ppm(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), nc(1:Qdof,1:nbDim),  &
         elem%xi(ie, 1:Qdof, 1:nbDim), &
         Ppm(1:Qdof, 1:nbDim, 1:ndim, 1:ndim), 1./elem%area, elem  )

    ! exrapolation wii = wi
    Ppm(1:Qdof,1,1:ndim, 1:ndim) = Ppm(1:Qdof,1,1:ndim, 1:ndim) + Ppm(1:Qdof,2,1:ndim, 1:ndim)


    if(state%nlSolver%implicitly) then
       ! evaluation of matrix terms - diagonal block
       call EvalEdgeBlockBB(elem, ie, elem, ie, Ppm(1:Qdof, 1, 1:ndim, 1:ndim), &
            elem%block(0)%Mb(1:ndim *dof, 1:ndim *dof) )

    else ! explicit discretization

       allocate(f_s(1:Qdof,1:ndim), ident(1:Qdof))
       ident(:) = 1.

       do i=1,Qdof
          f_s(i,1:ndim) = matmul(Ppm(i, 1, 1:ndim, 1:ndim), wi(i, 1:ndim) )
       enddo

       call EvalEdgeVectorB(elem, ie, ident(1:Qdof), -f_s(1:Qdof,1:ndim), &
            dofA, elem%vec(rhs, 1:ndim*dofA ) )

    endif


    deallocate(Ppm )
    deallocate(nc)
    deallocate(wi)

  end subroutine ElementInviscidNeumannlEdge

  !> boundary condition on fixed wall taking into account the eikonal velocity
  subroutine Set_W_pedes(elem, ie, Qdof, wi, nc, pedesW)
    class(element), intent(in) :: elem
    integer, intent(in) :: ie, Qdof
    real, dimension(1:Qdof,1:ndim), intent(inout) :: wi          !  w in integ nodes
    real, dimension(1:Qdof,1:nbDim), intent(in) :: nc    ! normals in integ nodes
    logical, intent(in)  :: pedesW
    real :: vn, val, vec, rel
    integer :: j

    !print*, '! v = v - 2*(v . un) un = v - 2*(v . n) n /|n|^2'
    do j=1,Qdof
       vec = sqrt( dot_product( elem%xi(ie, j, 1:2),  elem%xi(ie, j, 1:2) ) ) ! eikonal velocity
       val = sqrt( dot_product(wi(j, 2:3), wi(j, 2:3) ))                      ! actual velocity
       vn = dot_product(wi(j, 2:3), nc(j, 1:nbDim))    ! normal velocity
       rel = sqrt( (elem%xc(1) - 35)**2 + (elem%xc(2) - 5)**2)

       ! velocity in the direction of eikonal and size of the actual
       ! NEW
       if(vec > 0) then
          wi(j, 2:3)  =  elem%xi(ie, j, 1:2) / vec * val
       else
          wi(j,2:3) = 0.
       endif


       ! write(94, *) elem%xc(1:2), elem%xi(ie, j, 1:2)
       ! write(94, *) elem%xc(1:2) + elem%xi(ie, j, 1:2)
       ! write(94, *) '  '

       ! write(93, *) elem%xc(1:2)
       ! write(93, *) elem%xc(1:2) + wi(j, 2:3)/ wi(j,1)*100
       ! write(93, *) '  '

       if(pedesW) then
          ! we prescribe the direction velocity as the eikonal one
          ! the magnitude rest the same

          !!!wi(j, 2:3) =  elem%xi(ie, j, 1:2) / vec * val

          !write(92, *) elem%xc(1:2)
          !write(92, *) elem%xc(1:2) + wi(j, 2:3)/ wi(j,1)*100
          !write(92, *) '  '

       else
          ! NEW
!          wi(j,2:3) = wi(j, 2:3) - 2* vn * nc(j, 1:nbDim)

          !write(91, *) elem%xc(1:2)
          !write(91, *) elem%xc(1:2) + wi(j, 2:3)/ wi(j,1)*100
          !write(91, *) '  '
       endif
    enddo

  end subroutine Set_W_pedes


  !> evaluation of inviscid edge integrals on solid walls using the
  !> mirror boundary conditions
  subroutine ElementInviscidWall_2_Edge(elem, ie, Set_Ppm)
    type(element), intent(inout):: elem      ! elem = element
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    interface
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
    end interface
    real, dimension(:,:), allocatable :: wi, wii     !  w in integ nodes
    real, dimension(:,:), allocatable :: f_s         !  inviscid fluxes
    real, dimension(:), allocatable :: ident      !  indentity vector
    real, dimension(:,:), allocatable :: nc    ! normals in integ nodes
    real, dimension(:,:,:,:), allocatable :: Ppm   ! matrices Ppm in integ nodes
    real, dimension(1:nbDim) :: xi
    logical :: pedesW

    integer ::  dof, dofA, Qdof
    integer :: j, k, k1, row, col, i

    dof = elem%dof
    dofA = dof
    if(elem%deg_plus) dofA = elem%dof_plus

    Qdof = elem%face(fGdof,ie)  ! "num" integration nodes

    pedesW = .false.
    ! the outlet (right) part of the circle obstacle
    !if( elem%xc(1) > 35 .and. elem%xc(1) < 38 .and. abs(elem%xc(2) - 5) < 3) pedesW = .true.
    ! NEW
    !if(state%modelName == 'pedes' ) pedesW = .true.


    ! setting of outer normals in integration nodes
    allocate(nc(1:Qdof, 1:nbDim) )
    if(elem%ibcur > 0) then
       nc(1:Qdof,1:nbDim) = elem%nc(1:Qdof,1:nbDim)
    else
       nc(1:Qdof,1) = elem%n(ie,1)
       nc(1:Qdof,2) = elem%n(ie,2)
    endif


    allocate(wi(1:Qdof,1:ndim), wii(1:Qdof,1:ndim))
    call Eval_w_Edge(elem, ie, wi, .false.)
    wii(1:Qdof,1:ndim) = wi(1:Qdof,1:ndim)

    ! using BC for the modification of wi
    if(pedesW) then
       ! NEW
       call Set_W_pedes(elem, ie, Qdof, wii(1:Qdof, 1:ndim), nc(1:Qdof, 1:2), pedesW)

    elseif(state%model%Re == 0. .and. ndim >= 4) then
       call UpdateMirror(ndim, Qdof, wii, nc)

    elseif(state%modelName == 'pedes') then
       call UpdateMirror(ndim, Qdof, wii, nc)

    else
!       wii(1:Qdof, 2:3) = 0.
    endif

    !write(77,'(4es12.4)')  elem%xc(1:nbDim)
    !write(77,'(4es12.4)')  elem%xc(1:nbDim)+ 0.05* wi(1, 2:3)
    !write(77,*)


    allocate(Ppm(1:Qdof, 1:nbDim, 1:ndim, 1:ndim) )
    call Set_Ppm(ndim, nbDim, Qdof, wii(1:Qdof,1:ndim), nc(1:Qdof,1:nbDim), &
         elem%xi(ie, 1:Qdof, 1:nbDim),  &
         Ppm(1:Qdof,1:nbDim, 1:ndim, 1:ndim), 1./elem%area, elem )


    !do j=1,Qdof
    !   do k=1,ndim
    !      write(*,'(4es11.4, a3, 4es11.4)') Ppm(j,1,k, 1:ndim),' | ', &
    !           Ppm(j,2,k, 1:ndim)
    !   enddo
    !   write(*,*) 'Elem = ',elem%i,',  Qdof = ', j
    !   print*,'-------------------------------------------------'
    !enddo


    if(elem%ibcur > 0) then  ! UNIT normal
       nc(1:Qdof,1) = elem%nc(1:Qdof,1)/ elem%dnc(1:Qdof)
       nc(1:Qdof,2) = elem%nc(1:Qdof,2)/ elem%dnc(1:Qdof)
    else
       nc(1:Qdof,1) = elem%n(ie,1) / elem%dn(ie)
       nc(1:Qdof,2) = elem%n(ie,2) / elem%dn(ie)
    endif

    !!elem%block(0)%Mb(:,:) = 0.
    !!elem%vec(rhs, :) = 0.


    if(state%nlSolver%implicitly) then

       if(pedesW) then
          ! the linearization is P+ + P-
          ! NEW
          !Ppm(1:Qdof, 1, 1:ndim, 1:ndim) = Ppm(1:Qdof, 1, 1:ndim, 1:ndim) &
          !     + Ppm(1:Qdof, 2, 1:ndim, 1:ndim)

       else ! standard type of linearization reflecting the mirror

          Ppm(1:Qdof, 1, 1:ndim, 1:ndim) = Ppm(1:Qdof, 1, 1:ndim, 1:ndim) &
               + Ppm(1:Qdof, 2, 1:ndim, 1:ndim)


          ! adding of "mirror" part of impermeable BC
          if(state%model%Re <= 0.) then     ! inviscid

             do k=1,ndim
                ! column 2
                Ppm(1:Qdof, 1, k, 2) = Ppm(1:Qdof, 1, k, 2) &
                     - 2 * Ppm(1:Qdof, 2, k, 2) * nc(1:Qdof,1) * nc(1:Qdof,1) &
                     - 2 * Ppm(1:Qdof, 2, k, 3) * nc(1:Qdof,1) * nc(1:Qdof,2)

                ! column 3
                Ppm(1:Qdof, 1, k, 3) = Ppm(1:Qdof, 1, k, 3) &
                     - 2 * Ppm(1:Qdof, 2, k, 2) * nc(1:Qdof,2) * nc(1:Qdof,1) &
                     - 2 * Ppm(1:Qdof, 2, k, 3) * nc(1:Qdof,2) * nc(1:Qdof,2)
             enddo

          else   ! viscous

             do k=1,ndim
                ! column 2
                Ppm(1:Qdof, 1, k, 2) = Ppm(1:Qdof, 1, k, 2) - 2. * Ppm(1:Qdof, 2, k, 2)

                ! column 3
                Ppm(1:Qdof, 1, k, 3) = Ppm(1:Qdof, 1, k, 3) - 2. * Ppm(1:Qdof, 2, k, 3)
             enddo

          endif

       endif

       ! evaluation of matrix terms - diagonal block
       call EvalEdgeBlockBB(elem, ie, elem, ie, Ppm(1:Qdof, 1, 1:ndim, 1:ndim), &
            elem%block(0)%Mb(1:ndim *dof, 1:ndim *dof) )

    else     ! explicitly
       allocate(f_s(1:Qdof,1:ndim), ident(1:Qdof))
       ident(:) = 1.

       ! f_s = P^+ w_i
       do i=1,Qdof
          f_s(i,1:ndim) = matmul(Ppm(i, 1, 1:ndim, 1:ndim), wi(i, 1:ndim) )
       enddo

       ! f_s = P^i ^w_i,  ^w_i is the mirror of wi

       !write(*,'(a4,6es12.4)') 'wi :',wi(1, 1:4), nc(1, 1:nbDim)

       if(state%modelName == 'pedes' ) then
          ! BC from the eikonal velocity
          !call Set_W_pedes(elem, ie, Qdof, wi(1:Qdof, 1:ndim), nc(1:Qdof, 1:nbDim), pedesW )
          call Mirror_W(ndim, Qdof, wi(1:Qdof, 1:ndim), nc(1:Qdof, 1:nbDim))

       elseif(state%model%Re <= 0.) then     ! inviscid
          call Mirror_W(ndim, Qdof, wi(1:Qdof, 1:ndim), nc(1:Qdof, 1:nbDim))

       else ! viscous
          wi(1:Qdof, 2:3) = - wi(1:Qdof, 2:3)
       endif

       !write(*,'(a4,4es12.4)') 'wii:',wi(1, 1:4)

       do i=1,Qdof
          f_s(i,1:ndim) = f_s(i,1:ndim) +matmul(Ppm(i, 2, 1:ndim, 1:ndim), wi(i, 1:ndim) )
       enddo

       call EvalEdgeVectorB(elem, ie, ident(1:Qdof), -f_s(1:Qdof,1:ndim), &
            dofA, elem%vec(rhs, 1:ndim*dofA ) )

       !do i=1,elem%dof * ndim
       !   write(*,'(a6,16es14.6)') 'Iin', elem%vec(rhs, i), &
       !        dot_product(elem%block(0)%Mb(i, 1:ndim*dof), elem%w(0,1:ndim*dof) ), &
       !        elem%vec(rhs, i) + &
       !        dot_product(elem%block(0)%Mb(i, 1:ndim*dof), elem%w(0,1:ndim*dof) ), &
       !        f_s(1:Qdof, i),  dot_product(Ppm(1, 1, i, 1:ndim), wi(1, 1:ndim) )
       !enddo
       !write(*,*) '-------------------'

       deallocate(ident, f_s)
    endif


    deallocate(wi)
    deallocate(Ppm )

    deallocate(nc)

  end subroutine ElementInviscidWall_2_Edge

  !> extrapolate pressure from interior at outlet
  subroutine EvalOutputPressure()
    class(element), pointer :: elem
    real, dimension(:,:), allocatable :: weights
    real, dimension(:), allocatable :: pressure, identity, density
    real, dimension(:,:), allocatable :: wi
    real :: val, val1, valD
    integer :: ib, ie, je, j, Qdof, iBC

    allocate(weights(1:state%numBC, 1:3) )
    allocate(wi(1:state%space%max_Qdof,1:ndim))
    allocate(pressure(1:state%space%max_Qdof) )
    allocate(identity(1:state%space%max_Qdof) )
    allocate(density(1:state%space%max_Qdof) )
    identity(:) = 1.
    weights(:,:) = 0.

    do ib=1,grid%nbelm
       !if(grid%elem(grid%b_edge(ib)%itc)%tBC(grid%b_edge(ib)%jtc) >= 0 .and. &
       !     grid%elem(grid%b_edge(ib)%itc)%tBC(grid%b_edge(ib)%jtc) <= 2 ) then
       !if(grid%b_edge(ib)%BC == 1) then

       elem => grid%elem(grid%b_edge(ib)%itc)
       je = grid%b_edge(ib)%jtc
       Qdof = elem%face(fGdof,je)  ! "num" integration nodes
       iBC = elem%iBC(je)

       if(1 <= iBC .and. iBC <= state%numBC) then
          !write(*,'(a6,8i5)') '?????',ib,elem%i, je, Qdof, iBC
          !write(*,'(6es10.3)' ) elem%w(0,1:6)
          !write(*,'(6es10.3)' ) elem%w(0,7:12)
          !write(*,'(6es10.3)' ) elem%w(0,13:18)
          !write(*,'(6es10.3)' ) elem%w(0,19:24)

          call Eval_w_Edge(elem, je, wi(1:Qdof, 1:ndim), .false.)

          do j=1,Qdof
             density(j) = wi(j,1)
             pressure(j) = state%model%kappa1 *(wi(j,4) &
                  - dot_product(wi(j,2:3),wi(j,2:3))/wi(j,1)/2 )
             !write(71,*) elem%xc(1:nbDim), pressure(j)
             !write(*,'(a2,7es10.3)' ) '<<',wi(j,1:4), pressure(j)
          enddo

          call  IntegrateFunctionEdge(elem, je, pressure(1:Qdof), val)
          call  IntegrateFunctionEdge(elem, je, density(1:Qdof), valD)
          call  IntegrateFunctionEdge(elem, je, identity(1:Qdof), val1)
          weights(iBC,1) = weights(iBC,1) + val
          weights(iBC,2) = weights(iBC,2) + valD
          weights(iBC,3) = weights(iBC,3) + val1

          !write(*,'(a3,26es10.3)') '...', val/val1, pressure(1:Qdof)

       endif
    enddo

    do ie=1,state%numBC
       if(weights(ie, 3) > 0.) state%BC(ie)%press_extrap = weights(ie, 1) / weights(ie, 3)
       if(weights(ie, 3) > 0.) state%BC(ie)%rho_extrap = weights(ie, 2) / weights(ie, 3) ! desnity
       !write(*,'(a8,i5,3es14.6)' ) &
       !     '@@@@@',ie,state%BC(ie)%press_extrap, weights(ie, 1), weights(ie, 2)

    enddo

    deallocate(weights, pressure, wi, identity, density)

  end subroutine EvalOutputPressure
  ! wet steam equations
  !> Evaluation of inviscid fluxes on I/O edges
  subroutine ElementInviscidIOEdge(elem, ie, Set_Ppm)
    type(element), intent(inout):: elem      ! elem = element
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    interface
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
    end interface
    real, dimension(:,:), allocatable :: wi          !  w in integ nodes
    real, dimension(:,:), pointer:: nc           ! outer normsl in integ nodes
    real, dimension(:,:), pointer:: phi          ! local store arrays
    real, dimension(:,:), allocatable :: wD, w_BC, f_s ! w recomputed  in integ nodes
    real, dimension(:,:,:,:), allocatable :: Ppm   ! matrices Ppm in integ nodes
    real, dimension(:), allocatable :: ident, temp ! local store arrays
    real, dimension(1:nbDim) :: xi
    integer ::  dof, dofA, Qdof, ndimL
    integer :: i, j, k, k1, row, col, ibb

    dof = elem%dof
    dofA = dof
    if(elem%deg_plus) dofA = elem%dof_plus

    Qdof = elem%face(fGdof,ie)  ! "num" integration nodes

    allocate(wi(1:Qdof,1:ndim))
    call Eval_w_Edge(elem, ie, wi, .false.)

    ! setting of outer normals in integration nodes
    allocate(nc(1:Qdof, 1:nbDim) )
    if(elem%ibcur > 0) then
       nc(1:Qdof,1:nbDim) = elem%nc(1:Qdof,1:nbDim)
    else
     ! no courved edge on inlet Outlet
       nc(1:Qdof,1) = elem%n(ie,1)
       nc(1:Qdof,2) = elem%n(ie,2)
    endif

    allocate(Ppm(1:Qdof, 1:nbDim, 1:ndim, 1:ndim) )

    call Set_Ppm(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), nc(1:Qdof,1:nbDim), &
         elem%xi(ie, 1:Qdof, 1:nbDim),  &
         Ppm(1:Qdof,1:nbDim, 1:ndim, 1:ndim), 1./elem%area, elem )

    !if(elem%i == -8) then
    !   print*,'----------------sedt54aeXXXXX---'
    !   do i=1,Qdof
    !      write(*,'(10es16.8)') elem%xi(ie, i, 1:2), wi(i,1), Ppm(i, 1:2, 1, 1)
    !   enddo
    !endif

    if(state%nlSolver%implicitly) then
       ! evaluation of matrix terms - diagonal block
       call EvalEdgeBlockBB(elem, ie, elem, ie, Ppm(1:Qdof, 1, 1:ndim, 1:ndim), &
            elem%block(0)%Mb(1:ndim *dof, 1:ndim *dof) )
    endif


    ! outside of Omega ==> RHS
    ! setting of boundary conditions
    allocate( wD(1:Qdof, 1:ndim), w_BC(1:Qdof, 1:ndim) )
    allocate(temp(1:Qdof) )


    if(state%modelName == 'NSe' ) then
       if(state%local_problem)  then ! local problem for EE, BC given apriori in local_problem2.f90
          w_BC(1:Qdof, 1:ndim) =  elem%wSS(ie, 1:Qdof, 1:ndim)

       elseif(state%type_IC == 2) then     ! double Mach reflection

          do k=1,Qdof
             xi(1:nbDim) = grid%b_edge(-elem%face(neigh,ie))%x_div(k, 1:nbDim)
             call DMR_IC(xi(1:nbDim), w_BC(k,1:ndim), state%time%ctime )
          enddo

       elseif(state%type_IC == 4) then     ! RINGLEB FLOW problem

          do k=1,Qdof
             xi(1:nbDim) = grid%b_edge(-elem%face(neigh,ie))%x_div(k, 1:nbDim)
             call Exact_Ringleb(xi(1:nbDim), w_BC(k,1:ndim) )
          enddo

       elseif(state%type_IC == 9) then     ! steady state solution in the channel

          do k=1,Qdof
             xi(1:nbDim) = grid%b_edge(-elem%face(neigh,ie))%x_div(k, 1:nbDim)
             call Exact_SteadyChannel(xi(1:nbDim), w_BC(k,1:ndim) )
          enddo

       else ! external flow or channel

          if(elem%iBC(ie) == 0 ) then  ! viscous impermeable edge
             w_BC(1:Qdof, 1) = wi(1:Qdof, 1)
             w_BC(1:Qdof, 2) = 0.
             w_BC(1:Qdof, 3) = 0.
             w_BC(1:Qdof, 4) = wi(1:Qdof, 4) &
                  - 0.5*(wi(1:Qdof,2)**2+ wi(1:Qdof,3)**2)/wi(1:Qdof,1)

          else ! inlet/outlet, viscous or inviscid

             do k=1,ndim
                wD(1:Qdof,k) = state%BC(elem%iBC(ie))%ww(k)   ! includes wet steam case, BC are given from *.ini
             enddo


             ! repreparation for viscous flows
             !if(state%model%Re > 0. ) then

             if(  elem%tBC(ie) == 2) then
                ! repreparation for channel outlet
                call ReprepareBCCharacteristic(Qdof, ndim, wi(1:Qdof,1:ndim), &
                     wD(1:Qdof,1:ndim), elem%n(ie,:), &
                     state%BC(elem%iBC(ie))%press_extrap, elem%xc(:) )
                     !0., elem%xc(:) )

                !write(171,*)grid%x(grid%b_edge(-elem%face(neigh,ie))%lbn(1), 1:2), &
                !     pressure(4,wD(1, 1:ndim))

                !!   w_BC(1:Qdof,1:ndim) = wD(1:Qdof,1:ndim)
                ! do nothing
                !wD(1:Qdof,1:ndim) = wi(1:Qdof,1:ndim)
             endif
             !else

             !if(elem%tBC(ie) == 0 .or. elem%tBC(ie) == 1) then

             ! far-field BC
             ! approach based on the solution of the exact Riemann problem
             call SetBCexactRiemann(Qdof, ndim, wi(1:Qdof,1:ndim), &
                  wD(1:Qdof,1:ndim), w_BC(1:Qdof,1:ndim), elem%n(ie,:),  &
                  elem%xc(1:nbDim) )

                !write(*,'(a8,8es12.4)') '#i W_D',wD(1, 1:ndim),  pressure(4,wD(1, 1:ndim))
                !write(*,'(a8,8es12.4)') '#i W_i',wi(1, 1:ndim),  pressure(4,wi(1, 1:ndim))
                !write(*,'(a8,8es12.4)') '#i W_BC',w_BC(1, 1:ndim),  pressure(4,w_BC(1, 1:ndim))
                !print*

                ! approach based on the linearized Riemann problem
                !   call SetBCCharacteristic(Qdof, ndim, wi(1:Qdof,1:ndim), &
                !        wD(1:Qdof,1:ndim), w_BC(1:Qdof,1:ndim), elem%n(ie,:), &
                !         elem%xc(:), elem )

            ! elseif(elem%tBC(ie) == 2) then
                ! channel outlet

                ! !write(*,'(a8,8es12.4)') '#o W_D',wD(1, 1:ndim),  pressure(4,wD(1, 1:ndim))

                ! ! approach based on the physical extrapolation: sub-/super-sonic in/out-let
                ! call ReprepareBCCharacteristic(Qdof, ndim, wi(1:Qdof,1:ndim), &
                !      wD(1:Qdof,1:ndim), elem%n(ie,:), &
                !      state%BC(elem%iBC(ie))%press_extrap, elem%xc(:) )
                ! !w_BC(1:Qdof,1:ndim) = wD(1:Qdof,1:ndim)

                ! call SetBCexactRiemann(Qdof, ndim, wi(1:Qdof,1:ndim), &
                !      wD(1:Qdof,1:ndim), w_BC(1:Qdof,1:ndim), elem%n(ie,:),  &
                !      elem%xc(1:nbDim) )

                ! !if(elem%i == 60) then
                ! !   write(*,'(a8,8es12.4)') '#o W_i',wi(1, 1:ndim),  pressure(4,wi(1, 1:ndim)) , &
                ! !        state%BC(elem%iBC(ie))%press_extrap, 3.01044889
                ! !   write(*,'(a8,8es12.4)') '#o W_BC',w_BC(1, 1:ndim) , pressure(4,w_BC(1, 1:ndim))
                ! !   print*, elem%i
                ! !endif

             !endif

          ! elseif(ndim > 4 .and. ndim <= 6) then ! turbulence model
          !    ! (exact) characteristic BC
          !    !print*,'#############',ndim
          !    ndimL = 4
          !    call SetBCexactRiemann(Qdof, ndimL, wi(1:Qdof,1:ndimL), &
          !         wD(1:Qdof,1:ndimL), w_BC(1:Qdof,1:ndimL), elem%n(ie,:),  &
          !         elem%xc(1:nbDim) )

          !    ! extrapolation for turbulence
          !    w_BC(1:Qdof, ndimL+1: ndim) = wD(1:Qdof, ndimL+1:ndim)

          endif
       endif

    elseif(state%modelName == 'pedes' ) then
       ! setting of BC conditions based using the extrapolation
       !call Pedestrian_BC_Extrapolated(Qdof, ndim, wi(1:Qdof,1:ndim), &
       !        wD(1:Qdof,1:ndim), elem%n(ie,:), &
       !        state%BC(elem%iBC(ie))%press_extrap, elem%xc(:) )


       wD(:,1) = state%BC(elem%iBC(ie))%ww(1)
       wD(:,2) = state%BC(elem%iBC(ie))%ww(2)
       wD(:,3) = state%BC(elem%iBC(ie))%ww(3)

       !write(*,'(a4, 200es12.4)') 'OBCP:',wD(:, 1)
       !write(*,'(a4, 200es12.4)') 'OBCP:',wD(:, 2)
       !write(*,'(a4, 200es12.4)') 'OBCP:',wD(:, 3)
       !write(*,*) '  '

       !!wD(j,1) = wi(j,1)  ! wD density is kept from BC, velocity is free
       !wD(1:Qdof, 2) = wi(1:Qdof, 2) !/ wi(1:Qdof, 1) *  wD(1:Qdof, 1)
       !wD(1:Qdof, 3) = wi(1:Qdof, 3) ! / wi(1:Qdof, 1) *  wD(1:Qdof, 1)


       ! completly free output
       !wD(1:Qdof,1:ndim) = wi(1:Qdof,1:ndim)

       !write(*,'(a4, 200es12.4)') 'OBCP:',wi(:, 1)
       !write(*,'(a4, 200es12.4)') 'OBCP:',wi(:, 2)
       !write(*,'(a4, 200es12.4)') 'OBCP:',wi(:, 3)
       !write(*,*) '  '
       !write(*,'(a4, 200es12.4)') 'OBCP:',wD(:, 1)
       !write(*,'(a4, 200es12.4)') 'OBCP:',wD(:, 2)
       !write(*,'(a4, 200es12.4)') 'OBCP:',wD(:, 3)
       !write(*,*) ' -------------------------------------------- '

       ! setting of BC conditions based using the solution of linearized Riemann problem
       call Pedestrian_BC_LRP(Qdof, ndim, wi(1:Qdof,1:ndim), &
            wD(1:Qdof,1:ndim), elem%n(ie,:), &
            state%BC(elem%iBC(ie))%press_extrap, elem%xc(:) )
       !0., elem%xc(:) )

       ! write(171,*)grid%x(grid%b_edge(-elem%face(neigh,ie))%lbn(1), 1:2), &
       !      pressure(4,wD(1, 1:ndim))
       ! write(*,'(2es12.4, 3(a2,3es12.4))') &
       !      grid%x(grid%b_edge(-elem%face(neigh,ie))%lbn(1), 1:2), &
       !      '|', wi(1,1:ndim), &
       !      '|', state%BC(elem%iBC(ie))%ww(1:ndim),  &
       !      '|', wD(1, 1:ndim)

       w_BC(1:Qdof,1:ndim) = wD(1:Qdof,1:ndim)

       ! do nothing
       !w_BC(1:Qdof,1:ndim) = wi(1:Qdof,1:ndim)

    elseif(state%modelName == 'swe' ) then


       wD(:,1) = state%BC(elem%iBC(ie))%ww(1)
       wD(:,2) = state%BC(elem%iBC(ie))%ww(2)
       wD(:,3) = state%BC(elem%iBC(ie))%ww(3)


       ! setting of BC conditions based using the solution of linearized Riemann problem
       call Pedestrian_BC_LRP(Qdof, ndim, wi(1:Qdof,1:ndim), &
            wD(1:Qdof,1:ndim), elem%n(ie,:), &
            state%BC(elem%iBC(ie))%press_extrap, elem%xc(:) )


       w_BC(1:Qdof,1:ndim) = wD(1:Qdof,1:ndim)


    elseif(state%modelName == 'scalar' .or.state%modelName == '2eqs' ) then

       ! BC for scalar equation
       !print*, 'BC for scalar equation'
       do k=1,Qdof

          !xi(1:nbDim) = grid%b_edge(-elem%face(neigh,ie))%x_div(k, 1:nbDim)
          xi(1:nbDim) = elem%xi(ie, k, 1:nbDim)
          !print*, size(elem%xi(:,1,1)),size(elem%xi(1,:,1))
          !print*, xi(1:nbDim), w_BC(k,1:ndim), state%time%ctime
          call Exact_Scalar(xi(1:nbDim), w_BC(k,1:ndim), state%time%ctime )

          !write(*,'(a6,6e12.4)') 'wi:',xi(1:nbDim), w_BC(k,1:1),wi(k,1:1)
          !write(98,*) xi(1:nbDim), w_BC(k,1:1),wi(k,1:1)
       enddo


    elseif(state%modelName == "wet_steam") then  ! wet steam equations

       call SetBCexactRiemann(Qdof, 4, wi(1:Qdof,1:4), &
            wD(1:Qdof,1:4), w_BC(1:Qdof,1:4), elem%n(ie,:),  &   ! intent(out) :: w_BC
            elem%xc(1:nbDim) )

       ibb = state%BC(elem%iBC(ie))%inout    ! ie - inner edge 1,2,3

       if( ibb == 0 ) then  ! inlet
          w_BC(1:Qdof,5:8) = wD(1:Qdof, 5:8)  ! ???
          !!write(21, *)  elem%xc(:)
       elseif (ibb == 1 .or. ibb == 2) then  ! outlet
          w_BC(1:Qdof,5:8) = wi(1:Qdof, 5:8)  ! ???
          !!write(22,*) elem%xc(:)
       else
          print*,'@@ee@@',elem%iBC(ie), ibb
          stop
       endif

    else
       print*,'inv_fluxes: IO BC not implemented for ndim=',ndim
       stop
    endif

    allocate(f_s(1:Qdof,1:ndim), ident(1:Qdof))
    ident(:) = 1.

    do k=1,Qdof
       f_s(k,1:ndim) = matmul(Ppm(k,2,1:ndim, 1:ndim), w_BC(k, 1:ndim) )
    enddo

    if(.not. state%nlSolver%implicitly) then
       do k=1,Qdof
          f_s(k,1:ndim) = f_s(k,1:ndim) + matmul(Ppm(k,1,1:ndim, 1:ndim), wi(k, 1:ndim) )
       enddo
    endif


    !if(elem%i == 2)
    !write(*,'(a6,2i5,20es14.6)') 'CfluX=', elem%i,ie, -f_s(1:Qdof, 1:ndim)



    call EvalEdgeVectorB(elem, ie, ident(1:Qdof), -f_s(1:Qdof,1:ndim), &
         dofA, elem%vec(rhs, 1:ndim*dofA ) )

    !write(*,'(a4,2i5,42es10.2)') 'invP',elem%i, ie, Ppm(:,1,1:ndim, 1:ndim), &
    !     Ppm(:,2,1:ndim, 1:ndim), wi(:, 1:ndim), -f_s(1:Qdof,1:ndim)
    !write(*,'(a4,i5,l5,120es10.2)') 'invM',elem%i, state%nlSolver%implicitly, &
    !     Ppm(:,2,1:ndim, 1:ndim),w_BC(:,1:ndim), &
    !     elem%vec(rhs, 1:ndim*dofA)

   deallocate(wi, wD, w_BC)
   deallocate(nc)
   deallocate(Ppm )
   deallocate(temp)


  end subroutine ElementInviscidIOEdge

  !> evaluate of inviscid and viscous volumes via the form
  !> \f$ \int_{K} \vec{F(w, \nabla w)}\cdot \nabla \phi_i\, dS,\quad
  !> i=1,\dots, dofA\f$
  subroutine ElementInvVisVolumes(elem, Set_f_s, Set_R_s, dofA)
    type(element), intent(inout):: elem  ! elem = element,
    interface
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x, ie )
        integer, intent(in) :: Qdof, ndimL, nbDim
        real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
        real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
        real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
        integer, intent(in) :: ie
      end subroutine Set_f_s
      subroutine Set_R_s(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
        integer, intent(in) :: ndimL, nbDim, iRe, Qdof
        real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
        real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
        real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
        !real, intent(in) :: Re_1                     ! inverse of Reynolds number
        real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
      end subroutine Set_R_s
    end interface

    integer, intent(in) :: dofA          ! how many test functions to evaluate
    real, dimension(:,:), allocatable :: wi           ! solution  in integ nodes
    real, dimension(:,:,:), allocatable :: Dphi       ! test functions  in integ nodes
    real, dimension(:,:,:), allocatable :: Dwi        ! Der solution  in integ nodes
    real, dimension(:,:,:), allocatable :: f_s, R_s   ! fluxes in integ nodes
    real, dimension(:,:), allocatable :: Re_1
    real, dimension(:), allocatable ::  weights
    integer ::  Qdof, Qnum, k, i, row

    Qdof = elem%Qdof
    Qnum = elem%Qnum

    allocate(wi(1:Qdof,1:ndim) )
    call Eval_w_Elem(elem, wi(1:Qdof,1:ndim) )


    ! setting of fluxes f_s in integration nodes
    allocate(f_s(1:Qdof, 1:nbDim, 1:ndim))
    call Set_f_s(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), f_s(1:Qdof, 1:nbDim, 1:ndim), &
         elem%xi(0,1:Qdof, 1:nbDim), elem%i)

    ! setting of fluxes R_s in integration nodes
    if(state%model%Re > 0.) then
       allocate(  Dwi(1:Qdof,1:ndim,1:nbDim) )
       call Eval_Dw_Elem(elem, Dwi(1:Qdof, 1:ndim, 1:nbDim) )

       allocate(Re_1(1:iRe, 1:Qdof) )
       Re_1(1,1:Qdof) = 1./state%model%Re

       allocate(R_s(1:Qdof,1:nbDim, 1:ndim))
       call Set_R_s(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof, 1:ndim, 1:nbDim), Re_1, &
            R_s, elem%xi(0, 1:Qdof, 1:nbDim) )


       ! adding of terms
       f_s(1:Qdof, 1:nbDim, 1:ndim) = f_s(1:Qdof, 1:nbDim, 1:ndim) - R_s(1:Qdof, 1:nbDim, 1:ndim)

       deallocate(Dwi, R_s, Re_1)
    endif

    ! evaluation test functions on the edge
    allocate(Dphi( 1:dofA, 1:nbDim, 1:Qdof) )
    call Eval_Dphi(elem, dofA, Dphi)

    allocate(weights( 1:Qdof) )
    call Eval_V_Weights(elem, weights(1:Qdof)  )

    do k=1,ndim
       row = (k-1)*dofA

       do i=1,dofA
          elem%vec(rhs, row+i) = elem%vec(rhs, row+i) &
               + dot_product(weights(1:Qdof),  &
               f_s(1:Qdof, 1, k) * Dphi(i, 1, 1:Qdof) &
               + f_s(1:Qdof, 2, k) * Dphi(i, 2, 1:Qdof) )

          !if(i == 1) &
          !     write(*,'(a3,2i5,60es10.2)') 'VVV',elem%i,i, &
          !     weights(1:Qdof), Dphi(i, 1, 1:Qdof), Dphi(i, 2, 1:Qdof), &
          !     f_s(1:Qdof, 1, k),  f_s(1:Qdof, 2, k), &
          !     elem%vec(rhs,row+i)

       enddo
    enddo

    deallocate(wi, f_s, weights, Dphi)

  end subroutine ElementInvVisVolumes

  !> evaluate of inviscid and viscous fluxes through edge via the form
  !> \f$ \int_{\Gamma} \vec{F(w, \nabla w)}\cdot \vec{n} \phi_i\, dS,\quad
  !> i=1,\dots, dofA,\ \Gamma\subset \partial K\f$
  subroutine ElementInvVisFlux(elem, ie, Set_f_s, Set_R_s, dofA)
    type(element), intent(inout):: elem  ! elem = element,
    integer, intent (in) :: ie           ! inner index of edge, 1,2,3, (4)
    interface
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x, ie )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
         integer, intent(in) :: ie
      end subroutine Set_f_s
      subroutine Set_R_s(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
      end subroutine Set_R_s
    end interface
    integer, intent(in) :: dofA          ! how many test functions to evaluate
    real, dimension(:,:), allocatable :: nc           ! outer normal in integ nodes
    real, dimension(:,:), allocatable :: wi           ! solution  in integ nodes
    real, dimension(:,:), allocatable :: phi          ! test functions  in integ nodes
    real, dimension(:,:,:), allocatable :: Dwi        ! Der solution  in integ nodes
    real, dimension(:,:,:), allocatable :: f_s, R_s   ! fluxes in integ nodes
    real, dimension(:,:), allocatable :: Re_1
    real, dimension(:), allocatable :: temp
    integer ::  Qdof, Qnum, k, i, row

    Qdof = elem%face(fGdof,ie)  ! "num" integration nodes
    Qnum = elem%face(fGnum,ie)

    allocate(wi(1:Qdof,1:ndim) )
    call Eval_w_Edge( elem, ie,  wi(1:Qdof,1:ndim),  .false.)

    ! setting of outer normals in integration nodes
    allocate(nc(1:Qdof, 1:nbDim) )
    if(elem%ibcur > 0) then
       nc(1:Qdof,1:nbDim) = elem%nc(1:Qdof,1:nbDim)
    else
     ! no courved edge on inlet Outlet
       nc(1:Qdof,1) = elem%n(ie,1)
       nc(1:Qdof,2) = elem%n(ie,2)
    endif

    ! setting of fluxes f_s in integration nodes
    allocate(f_s(1:Qdof, 1:nbDim, 1:ndim))
    call Set_f_s(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), f_s(1:Qdof, 1:nbDim, 1:ndim), &
         elem%xi(ie ,1:Qdof, 1:nbDim), elem%i)

    ! setting of fluxes R_s in integration nodes
    if(state%model%Re > 0.) then
       allocate(  Dwi(1:Qdof,1:ndim,1:nbDim) )
       call Eval_Dw_Edge(elem, ie,  Dwi, .false.)

       allocate(Re_1(1:iRe, 1:Qdof) )
       Re_1(1,1:Qdof) = 1./state%model%Re

       allocate(R_s(1:Qdof,1:nbDim, 1:ndim))
       call Set_R_s(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof, 1:ndim, 1:nbDim), Re_1, &
            R_s(1:Qdof,1:nbDim, 1:ndim),  elem%xi(ie ,1:Qdof, 1:nbDim) )


    ! adding of terms
       f_s(1:Qdof, 1:nbDim, 1:ndim) = f_s(1:Qdof, 1:nbDim, 1:ndim) - R_s(1:Qdof, 1:nbDim, 1:ndim)

       deallocate(Dwi, R_s, Re_1)
    endif

    ! evaluation test functions on the edge
    allocate(phi( 1:dofA, 1:Qdof) )
    call Eval_Phi_Edge(elem, dofA, ie, phi, .false.)

    allocate(temp( 1:Qdof) )


    do k=1,ndim
       row = (k-1)*dofA
       temp(1:Qdof) = state%space%G_rule(Qnum)%weights(1:Qdof) &
            *(f_s(1:Qdof, 1, k) *  nc(1:Qdof,1) + f_s(1:Qdof, 2, k) *  nc(1:Qdof,2) )

       !write(*,'(a3,i5,20es10.2)') 'wi ',elem%i,wi(1:Qdof, 1:ndim)
       !write(*,'(a3,i5,20es10.2)') 'n1 ',elem%i,nc(1,1:nbDim)
       !write(*,'(a3,i5,20es10.2)') 'n2 ',elem%i,nc(2,1:nbDim)
       !write(*,'(a3,i5,20es10.2)') 'f1 ',elem%i,f_s(1:Qdof, 1, k)
       !write(*,'(a3,i5,20es10.2)') 'f2 ',elem%i,f_s(1:Qdof, 2, k)
       !write(*,'(a3,i5,20es10.2)') 'f.n',k,f_s(1:Qdof, 1, k) *  nc(1:Qdof,1) + f_s(1:Qdof, 2, k) *  nc(1:Qdof,2)
       !write(*,'(a3,i5,20es10.2)') '!!!',elem%i,state%space%G_rule(Qnum)%weights(1:Qdof)

       do i=1,dofA
          elem%vec(rhs, row+i) = elem%vec(rhs, row+i) &
               - dot_product(temp(1:Qdof), phi(i, 1:Qdof) )

          !if(i == 1) &
          !     write(*,'(a3,2i5,20es10.2)') '???',elem%i,i,phi(i, 1:Qdof),&
          !     temp(1:Qdof),dot_product(temp(1:Qdof), phi(i, 1:Qdof) ),  elem%vec(rhs,row+i)

       enddo
    enddo

    deallocate(wi, f_s, nc, temp, phi)

  end subroutine ElementInvVisFlux

  !> evaluate of inviscid edge integrals over inner edges
  subroutine ElementInviscidInnerEdge(elem, elem1, ie, Set_Ppm, Set_f_s)
    type(element), intent(inout):: elem, elem1  ! elem = element, elem1 = neigh elem
    class(element), pointer    :: elemP, elemP1  ! local pointers
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    interface
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
       end subroutine Set_Ppm
       subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x, iee )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
         integer, intent(in) :: iee
      end subroutine Set_f_s

    end interface

    real, dimension(:,:), allocatable :: wi, wii ! w in integ nodes
    real, dimension(:,:), allocatable :: nc      ! outer normal in integ nodes
    real, dimension(:,:,:,:), allocatable :: Ppm ! matrices Ppm in integ nodes
    real, dimension(:,:), allocatable :: f_s     !   vectors f_s
    real, dimension(:,:,:), allocatable :: ff_s  !   vectors f_s
    real, dimension(:), allocatable :: ident     !   identity
    real :: param, lambda_LF, tt, tt1
!    real, dimension(1:4) :: qL, qR
    integer ::  dof, dofA, dof1, ie1, Qdof, k, i, itest1

    itest1 = -grid%nelem

    dof = elem%dof
    dofA = dof
    if(elem%deg_plus) dofA = elem%dof_plus

    dof1 = elem1%dof

    ie1 = elem%face(nei_i,ie)

    Qdof = elem%face(fGdof,ie)

    allocate(wi(1:Qdof,1:ndim) )

    call Eval_aver_w_Edge(elem, elem1, ie, Qdof, wi(1:Qdof,1:ndim))

    ! for local problem, linearization with respect to the original solution
    !if(state%local_problem) then
    !   elemP  => grid%elem( abs(elem%i) )
    !   elemP1 => grid%elem( abs(elem1%i) )
    !   call Eval_aver_w_Edge(elemP, elemP1, ie, Qdof, wi(1:Qdof,1:ndim))
    !endif


    !  if (state%nlSolver%implicitly ) print*, 'WI: ', wi(1:Qdof,1:ndim)

    !if(elem%i == 2214 .and. state%time%iter_loc >= 19) then
    !   print*
    !   write(*,'(a4,1i5,20es12.4)') 'w_i ',elem%i, elem%w(0, :)
    !   write(*,'(a4,1i5,20es12.4)') 'w_ii',elem1%i, elem1%w(0, :)
    !   do i=1,ndim
    !      write(*,'(a4,2i5,20es12.4)') 'w_av',elem%i,elem1%i,wi(:, i)
    !   enddo
    !   print*,'-----------------------'
    !endif

    state%print = .false.

    ! setting of outer normals in integration nodes
    allocate(nc(1:Qdof, 1:nbDim) )
    !if(elem%ibcur > 0) then
    !   nc(1:Qdof,1:nbDim) = elem%nc(1:Qdof,1:nbDim)
    !else
     ! no courved edge on inner edges
       nc(1:Qdof,1) = elem%n(ie,1)
       nc(1:Qdof,2) = elem%n(ie,2)
    !endif

       !print*,'ATTENTION ED6etd6eja76e   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       !nc(1:Qdof,1) = 0.
       !nc(1:Qdof,2) = 1.


    allocate(Ppm(1:Qdof, 1:nbDim, 1:ndim, 1:ndim) )

    !!call cpu_time(tt1)
      ! Vijajasyndaram
    call Set_Ppm( ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), nc(1:Qdof,1:nbDim), &
         elem%xi(ie, 1:Qdof, 1:nbDim),  &
         Ppm(1:Qdof,1:2, 1:ndim, 1:ndim), 1./elem%area, elem )

    !!call cpu_time(tt)
    !!if(elem%i == itest1) print*,'IN,',tt,' Ppm', tt - tt1

    !if(elem%i == -8) then
    !   print*,'----------------sedt54ae---'
    !   do i=1,Qdof
    !      write(*,'(10es16.8)') elem%xi(ie, i, 1:2), wi(i,1), Ppm(i, 1:2, 1, 1)
    !   enddo
    !endif

!    if (state%nlSolver%implicitly) then
!      print*, 'PPM:' , Ppm(1:Qdof, 1, 1:ndim, 1:ndim)
!      stop
!    endif

    ! if( elem%i == 11) then
    !    write(*,'(a10,10es16.8)') 'ni:, wi:',nc(1,1:nbDim), wi(1,1:ndim)
    !    print*
    !    do i=1,ndim
    !       write(*,'(a3,i5,10es16.8)') 'Pp+',i, Ppm(1,1, i, 1:ndim), &
    !            dot_product(Ppm(1,1, i, 1:ndim), wi(1,1:ndim) )
    !    enddo
    !    print*
    !    do i=1,ndim
    !       write(*,'(a3,i5,10es16.8)') 'Pp-',i, Ppm(1,2, i, 1:ndim), &
    !            dot_product(Ppm(1,2, i, 1:ndim), wi(1,1:ndim) )
    !    enddo
    !    print*
    !    do i=1,ndim
    !       write(*,'(a3,i5,10es16.8)') 'PP',i, Ppm(1,1, i, 1:ndim) + Ppm(1,2, i, 1:ndim), &
    !            dot_product(  Ppm(1,1, i, 1:ndim) + Ppm(1,2, i, 1:ndim),  wi(1,1:ndim) )
    !    enddo
    !    print*
    !    !  Ppm(1,1, 1:ndim, 1:ndim)= Ppm(1,1, 1:ndim, 1:ndim) + Ppm(1,2, 1:ndim, 1:ndim)

    !    !  do i=1,ndim
    !    !     write(*,'(a3,i5,10es16.8)') 'PP',i, Ppm(1,1, i, 1:ndim)
    !    !  enddo
    !    !  print*
    !    !   do i=1,ndim
    !    !     write(*,'(a3,i5,10es16.8)') 'fs',i, dot_product(Ppm(1,1, i, 1:ndim), wi(1, 1:ndim) )
    !    !  enddo


    !    !  print*,'...........'
    !    !  stop
    ! endif



    ! Lax-Friedrichs, setting of matrixes A_s in integration nodes
    !call Set_A_s_Euler(ndim, Qdof, wR(1:Qdof,1:ndim),  &
    !     state%A_s(1:Qdof,1:nbDim,1:ndim,1:ndim) , elem%xi(0,1:Qdof, 1:nbDim))

    !
    !Ppm(1:Qdof,1,1:ndim,1:ndim)= (state%A_s(1:Qdof,1,1:ndim,1:ndim)*elem%n(ie,1) &
    !     +  state%A_s(1:Qdof,2,1:ndim,1:ndim)*elem%n(ie,2))/2.
    !
    !call Set_A_s_Euler(ndim, Qdof, wC(1:Qdof,1:ndim),  &
    !     state%A_s(1:Qdof,1:nbDim,1:ndim,1:ndim), elem%xi(0,1:Qdof, 1:nbDim) )
    !
    !Ppm(1:Qdof,2,1:ndim,1:ndim)= (state%A_s(1:Qdof,1,1:ndim,1:ndim)*elem%n(ie,1) &
    !     +  state%A_s(1:Qdof,2,1:ndim,1:ndim)*elem%n(ie,2))/2.


    !lambda_LF = state%time%tau(1)/ elem%dn(ie)
    !lambda_LF = 100./ elem%dn(ie)

    ! stabilization
     if(.not. state%local_problem) then
        param = (elem%rezid + elem1%rezid)/2.* elem%dn(ie)**state%ST_Ep * state%ST_Ec &
             * elem%dn(ie)  ! Ppm is already multiplied by |n| = dn
     else
        param = 0.
     endif


    !if(elem%i < elem1%i) &
    !     print*,'I', elem%i, elem1%i, elem%dn(ie),& ! elem%rezid, elem1%rezid,
    !     param/elem%dn(ie)

    do k=1,ndim                 ! k = index of component of w
       !!write(state%time%iter+51,*) elem%xc(1:nbDim), Ppm(1:Qdof, 1, k, k), param ,k,elem%i
       Ppm(1:Qdof, 1, k, k) = Ppm(1:Qdof, 1, k, k) + param !!+ 1./lambda_LF
       Ppm(1:Qdof, 2, k, k) = Ppm(1:Qdof, 2, k, k) - param !!- 1./lambda_LF
    enddo

    if(state%nlSolver%implicitly) then
!!    elem%block(0)%Mb(:,:) = 0.
!!    elem%block(ie)%Mb(:,:) = 0.
!!    elem%vec(rhs,:) = 0.

       ! evaluation of matrix terms - diagonal block
        !            print* , 'before1:', ie, elem%block(0)%Mb(1,1)
        !            print*, 'PPM:' , Ppm(1:Qdof, 1, 1:ndim, 1:ndim)

       !print*,'wes34ewsd',ie, dof, ndim
       !print* , 'before1:', ie , elem%block(0)%Mb(1,1)
       call EvalEdgeBlockBB(elem, ie, elem, ie, Ppm(1:Qdof, 1, 1:ndim, 1:ndim), &
            elem%block(0)%Mb(1:ndim *dof, 1:ndim *dof) )
       !print*, elem%xc(:)
       !print*, elem1%xc(:)
       !print* , 'after1:', ie , elem%block(0)%Mb(1,1)

       ! evaluation of matrix terms - off-diagonal block
       !!!!if(.not.  state%only_diag) &
       call EvalEdgeBlockBB(elem, ie, elem1, ie1, Ppm(1:Qdof, 2, 1:ndim, 1:ndim), &
            elem%block(ie)%Mb(1:ndim*dof, 1:ndim*dof1) )

    else  ! explicitly, i.e,  state%nlSolver%implicitly = .false.
       allocate(f_s(1:Qdof,1:ndim), wii(1:Qdof,1:ndim), ident(1:Qdof))
       ident(:) = 1.

       call Eval_w_Edge(elem,  ie,  wi,  .false.)
       call Eval_w_Edge(elem1, ie1, wii, .true.)
       !!if(state%local_problem) wii(1:Qdof,1:ndim) = wi(1:Qdof,1:ndim)

       if(state%num_flux) then  ! numerical flux used, standard way
          do i=1,Qdof
             f_s(i,1:ndim) = matmul(Ppm(i,1,1:ndim, 1:ndim), wi(i, 1:ndim) ) &
                  + matmul(Ppm(i,2,1:ndim, 1:ndim), wii(i, 1:ndim) )

          enddo
       else         ! we used pysical flux, only for error estimates
          do i=1,Qdof
             f_s(i,1:ndim) = matmul(Ppm(i,1,1:ndim, 1:ndim), wi(i, 1:ndim) ) &
                  + matmul(Ppm(i,2,1:ndim, 1:ndim), wii(i, 1:ndim) )
          enddo

          !if(elem%i == 8 .and. ie == 2) &
          !     write(*,'(a6,2i5,20es12.4)') 'Nflux:',elem%i,ie,f_s(:,1)

          allocate(ff_s(1:Qdof,1:nbDim, 1:ndim) )

          call Set_f_s(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), ff_s(1:Qdof, 1:nbDim, 1:ndim),&
               elem%xi(ie,1:Qdof, 1:nbDim), elem%i)

          f_s(1:Qdof,1:ndim) = ff_s(1:Qdof, 1, 1:ndim) * elem%n(ie, 1)  &
               + ff_s(1:Qdof, 2, 1:ndim) * elem%n(ie, 2)  !  ONLY STRAIGHT EDGES

          !if(elem%i == 8 .and. ie == 2) then
          !write(*,'(a6,2i5,20es12.4)') 'Fflux:',elem%i,ie,f_s(:,1)
          !   write(*,'(a6,2i5,20es12.4)') 'w i:',elem%i,ie,wi(:,:),wii(:,:)
          !   write(*,'(a6,2i5,20es12.4)') 'F n :',elem%i,1,ff_s(1:Qdof, 1, 1:ndim), elem%n(ie, 1)
          !   write(*,'(a6,2i5,20es12.4)') 'F n :',elem%i,2,ff_s(1:Qdof, 2, 1:ndim), elem%n(ie, 2)
          !   print*
          !endif

          deallocate(ff_s)
       endif

       !if(elem%i >= 99 .and. elem%i <= 100) then
       !   !write(*,'(a4,i5,12es10.2)') 'VOL ',elem%i, wi(:, 1:ndim)
       !   write(22,'(a4,2i5,120es14.6)') 'fn ',elem%i, ie, f_s(:,1:ndim)
       !   !write(*,'(a4,i5,12es10.2)') '    ',dofA, elem%vec(rhs, 1:dofA)
       !endif

       !do i=1,ndim
       !   write(*,'(a4,2i5,20es12.4)') 'f_s ',elem%i,i,f_s(:, i)
       !enddo

       !if( elem%i == 11 )write(*,'(a8,2i5,20es12.4)') 'CfluA=', elem%i,ie, -f_s(1, 1:ndim)

       call EvalEdgeVectorB(elem, ie, ident(1:Qdof), -f_s(1:Qdof,1:ndim), &
            dofA, elem%vec(rhs, 1:ndim*dofA ) )

       deallocate(f_s, wii, ident)
    endif

    deallocate(Ppm )
    deallocate (wi)
    deallocate (nc)
  end subroutine ElementInviscidInnerEdge

  !> computation of inviscid volume integrals on \f$elem=K\f$, i.e.,
  !> \f$-\int_K \sum_{s=1}^2 (A_s(w) \phi_C) \partial_s \phi_R dx \f$
  subroutine ElementInviscidVolumes(elem, Set_f_s,  Set_A_s)
    type(element):: elem
    ! Compute inviscid fluxes f_s and matrices A_s(w)
    interface
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x, ie )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
         integer, intent(in) :: ie
      end subroutine Set_f_s
      subroutine Set_A_s(ndimL, nbDim, Qdof, w, A_s, xi, ie)
         integer, intent(in) :: Qdof, nbdim, ndimL
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: A_s
         ! matrices A_s in  -- " --
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
         integer, intent(in) :: ie
      end subroutine
    end interface
    real, dimension(:,:), allocatable :: wi ! w recomputed  in integ nodes
    real, dimension(:,:,:,:), allocatable :: A_s ! matrices A_s
    real, dimension(:,:,:), allocatable :: f_s !   vectors f_s
    real :: param
    integer ::  Qdof, dof, dofA
    integer :: i

    dof = elem%dof
    dofA = dof
    if(elem%deg_plus) dofA = elem%dof_plus

    Qdof = elem%Qdof

    ! setting of the state vector in integration nodes
    allocate(wi(1:Qdof,1:ndim) )
    call Eval_w_Elem(elem, wi(1:Qdof,1:ndim) )


    if(state%nlSolver%implicitly) then
       ! setting of matrixes A_s in integration nodes
       allocate(A_s(1:Qdof,1:nbDim,1:ndim,1:ndim))
       call Set_A_s(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), &
            A_s(1:Qdof, 1:nbDim, 1:ndim, 1:ndim), elem%xi(0,1:Qdof, 1:nbDim), elem%i )

       ! for DUA error estimates
       if(.not. state%local_problem .and. (state%modelName == 'scalar'.or.state%modelName == '2eqs'))then
          param =  maxval(abs(A_s(1:Qdof,1:nbDim, 1:ndim, 1:ndim )) )
          elem%Cbo = max(elem%Cbo, param)
          elem%Cb  = max(elem%Cb,  grid%diam * param)
       end if

       ! evaluation of matrix terms
       call EvalBlockBD(elem,  -A_s(1:Qdof,1:nbDim,1:ndim,1:ndim), &
            elem%block(0)%Mb(1:ndim*dof, 1:ndim*dof ) )

       deallocate(A_s)
    else ! explicitly, i.e,  state%nlSolver%implicitly = .false.

       ! setting of fluxes f_s in integration nodes
       allocate(f_s(1:Qdof, 1:nbDim, 1:ndim))
       call Set_f_s(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), f_s(1:Qdof, 1:nbDim, 1:ndim), &
            elem%xi(0,1:Qdof, 1:nbDim), elem%i)

       call EvalVectorD(elem, f_s(1:Qdof,1:nbDim,1:ndim), dofA, elem%vec(rhs, 1:ndim*dofA ) )


       ! allocate(A_s(1:Qdof,1:nbDim,1:ndim,1:ndim))
       ! call Set_A_s(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), &
       !      A_s(1:Qdof, 1:nbDim, 1:ndim, 1:ndim), elem%xi(0,1:Qdof, 1:nbDim) )


       ! !if(elem%i >= 99 .and. elem%i <= 100) then
       ! !   !write(*,'(a4,i5,12es10.2)') 'VOL ',elem%i, wi(:, 1:ndim)
       ! if( abs(f_s(1,2, 1) ) > 1E-8) then
       !    do i=1,Qdof
       !       write(*,'(a4,i5,120es14.6)') 'f_s ',elem%i,  &
       !            f_s(i,1,1:ndim) - matmul(A_s(i, 1, 1:ndim, 1:ndim), wi(i, 1:ndim) ),  &
       !            f_s(i,2,1:ndim) - matmul(A_s(i, 2, 1:ndim, 1:ndim), wi(i, 1:ndim) )
       !    enddo
       ! endif

       ! deallocate(A_s)

       deallocate(f_s)
    endif

    deallocate(wi)
    !stop

  end subroutine ElementInviscidVolumes


  !> compute \f$ \int_K  param) \nabla \phi_j  \nabla d\phi_i dx \f$
  subroutine ElementVolumesStabil(elem)
    type(element):: elem
    real, dimension(:), allocatable :: Stab  ! w recomputed  in integration nodes
    real, dimension(:,:,:), allocatable :: Dwi  ! Dw recomputed  in integration nodes
    real, dimension(:,:,:), allocatable :: R_s  ! R_s(w) recomputed  in integration nodes
    integer :: i, k, k1, dof, dofA, Qdof
    integer :: row, col

    dof = elem%dof
    dofA = dof
    if(elem%deg_plus) dofA = elem%dof_plus

    Qdof = elem%Qdof

    allocate(Stab (1:Qdof) )

    Stab(1:Qdof) = elem%rezid * elem%diam**state%ST_Vp * state%ST_Vc

    if(state%nlSolver%implicitly) then
       do k=1,ndim                 ! k = index of component of w
          row = (k-1)*dof
          call IntegrateBlockD2(elem, elem%dof, Stab(1:Qdof), &
               elem%block(0)%Mb(row+1:row+dof, row+1:row+dof ) )
       enddo
    else
       ! evaluation of the gradient of w in integ nodes
       allocate(Dwi(1:Qdof, 1:ndim, 1:nbDim))
       call Eval_Dw_Elem(elem, Dwi(1:Qdof, 1:ndim, 1:nbDim) )

       allocate(R_s(1:Qdof, 1:nbDim, 1:ndim))
       do k = 1,ndim
          R_s(1:Qdof, 1, k) = Stab(1:Qdof) * Dwi(1:Qdof, k, 1)
          R_s(1:Qdof, 2, k) = Stab(1:Qdof) * Dwi(1:Qdof, k, 2)
       enddo

       call EvalVectorD(elem,-R_s(1:Qdof,1:nbDim,1:ndim), dofA, elem%vec(rhs, 1:ndim*dofA ) )

       deallocate(Dwi, R_s)
    endif

    deallocate(Stab)

  end subroutine ElementVolumesStabil


  !> compute \f$ \int_K \sum_{s=1}^2 K_sk(w) d\phi_j/d x_s   d\phi_i/d x_k dx \f$
  !> \f$ Re = \f$ Reynolds number, real or artificial
  subroutine ElementViscousVolumes(elem,  Set_R_s, Set_K_sk)
    type(element):: elem
    ! Compute viscous fluxes and matrices K_sk(w), s,k=1,2
    interface
      subroutine Set_R_s(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
      end subroutine Set_R_s
      subroutine Set_K_sk(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, K_sk, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,ndimL), intent(inout) :: K_sk
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
       end subroutine Set_K_sk
    end interface
    real, dimension(:,:), allocatable :: wi  ! w recomputed  in integration nodes
    real, dimension(:,:), allocatable :: Re_1  ! inverse of Reynolds number in integ nodes
    real, dimension(:,:,:,:,:), allocatable :: Ksk  ! viscous terms
    real, dimension(:,:,:), allocatable :: Dwi  ! Dw recomputed  in integration nodes
    real, dimension(:,:,:), allocatable :: R_s  ! R_s(w) recomputed  in integration nodes
    real :: param
    integer ::  dof,dofA, Qdof, k, i

    dof = elem%dof
    dofA = elem%dof
    if(elem%deg_plus) dofA = elem%dof_plus

    Qdof = elem%Qdof

    ! Reynolds number in integ nodes.
    allocate(Re_1(1:iRe, 1:Qdof) )

    if(state%model%Re > 0.) then
       Re_1(1, 1:Qdof) = 1./state%model%Re
    else
       call Eval_func_Elem(elem, elem%vec(aRe, 1:dof), Re_1(1, 1:Qdof) )
    endif

    ! setting the precomputed values, otherwise == 1
    !if(state%modelName == 'scalar' .and.  &     ! POROUS MEDIA FLOW
    !     (state%model%idiff == 13  ) ) then
    Re_1(2:iRe, 1:Qdof) = transpose( elem%xi(0, 1:Qdof, 2+1:2+iRe-1) )
    !write(*,'(a8, 30es12.4)') 'Re_1:', Re_1(2:iRe, 1)
    !endif


    ! solution in integ nodes
    allocate(wi(1:Qdof,1:ndim) )
    call Eval_w_Elem(elem, wi(1:Qdof,1:ndim) )

    ! derivative of the solution in integ nodes
    allocate(Dwi(1:Qdof, 1:ndim, 1:nbDim))

    if(state%nlSolver%implicitly) then

       ! diffusion of the scalar problem can have diffusion depending on the gradient of the solution
       if(state%modelName == 'scalar' .or. state%modelName == '2eqs' .or. state%modelName == 'porous') &
            call Eval_Dw_Elem(elem, Dwi(1:Qdof, 1:ndim, 1:nbDim) )


       allocate(Ksk (1:Qdof,1:nbDim, 1:nbDim, 1:ndim, 1:ndim ) )
       call Set_K_sk(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof, 1:ndim, 1:nbDim), &
            Re_1(1:iRe, 1:Qdof), &
            Ksk(1:Qdof,1:nbDim, 1:nbDim, 1:ndim, 1:ndim ), elem%xi(0, 1:Qdof, 1:nbDim) )

       ! if(state%time%iter >= 4) then
       !     do k=1, Qdof
       !        write(64,'(30es12.4)')  &
       !             elem%xi(0, k, 1:2), Ksk(k, 1, 1, 1, 1), Ksk(k, 2, 2, 1, 1), wi(k, 1),&
       !             Re_1(2:iRe, k)
       !     enddo
       !     if(elem%i == grid%nelem) stop '9ue93jdo3dmzd39u393i'
       ! endif

       !write(*,'(a8,6es12.4)') 'Ksk:',  Ksk(1,1:nbDim, 1:nbDim, 1:ndim, 1:ndim )

       ! for DUA error estimates
       if(.not. state%local_problem .and. (state%modelName == 'scalar'.or.state%modelName == '2eqs'))then
          param = maxval(abs(Ksk(1:Qdof,1:nbDim, 1:nbDim, :, :)))
          elem%CKo = max(elem%CKo, param )
          elem%CK  = max(elem%CK,  grid%diam / elem%diam * param )
!          print*, 'param in KSK', param
       endif

       if(state%ST_Vc >0 ) then ! volume  stabilization
          !if(elem%i == 1 .and. state%time%iter < 5) &
          !  print*,'Volume stabil for NS not implemented in inv_fluxes.f90  ERD'
!
          param = elem%rezid * elem%diam**state%ST_Vp * state%ST_Vc
          do k=1,ndim
             Ksk(1:Qdof,1,1,k,k) = Ksk(1:Qdof,1,1,k,k) + param
             Ksk(1:Qdof,2,2,k,k) = Ksk(1:Qdof,2,2,k,k) + param
          enddo
       endif

       ! evaluation of matrix terms
       call EvalBlockDD(elem, Ksk(1:Qdof,1:nbDim,1:nbDim,1:ndim,1:ndim), &
            elem%block(0)%Mb(1:ndim*dof, 1:ndim*dof ) )

       deallocate(Ksk)
    else
       ! evaluation of the gradient of w in integ nodes
       call Eval_Dw_Elem(elem, Dwi(1:Qdof, 1:ndim, 1:nbDim) )

       !if(elem%i >= 99 .and. elem%i <= 100) then
       !   write(*,'(a6, 2i5, 30es12.4)') 'Dwi:',elem%i,1,  Dwi(1:Qdof, 1:ndim, 1)
       !   write(*,'(a6, 2i5, 30es12.4)') 'Dwi:',elem%i,2,  Dwi(1:Qdof, 1:ndim, 2)
       !   write(*,'(a6, 2i5, 30es12.4)') 'Dwi:',elem%i,0, &
       !        (Dwi(1:Qdof, 1, 1)**2 + Dwi(1:Qdof, 1, 2)**2 )**0.5
       !endif

       ! setting of fluxes R_s in integration nodes
       allocate(R_s(1:Qdof, 1:nbDim, 1:ndim))
       call Set_R_s(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim),Dwi(1:Qdof, 1:ndim, 1:nbDim), Re_1, &
            R_s(1:Qdof, 1:nbDim, 1:ndim), elem%xi(0, 1:Qdof, 1:nbDim))

       if(state%ST_Vc >0 ) then
          param = elem%rezid * elem%diam**state%ST_Vp * state%ST_Vc
          do k=1,ndim
             R_s(1:Qdof,1,k) = R_s(1:Qdof,1,k) + param*Dwi(1:Qdof, k, 1)
             R_s(1:Qdof,2,k) = R_s(1:Qdof,2,k) + param*Dwi(1:Qdof, k, 2)
          enddo
       endif

       !do i=1,Qdof
       !   R_s(i,1,1:ndim) = matmul(Ksk(i,1,1,1:ndim, 1:ndim), Dwi(i, 1:ndim, 1) ) &
       !        + matmul(Ksk(i,1,2,1:ndim, 1:ndim), Dwi(i, 1:ndim, 2) )
       !
       !   R_s(i,2,1:ndim) = matmul(Ksk(i,2,1,1:ndim, 1:ndim), Dwi(i, 1:ndim, 1) ) &
       !        + matmul(Ksk(i,2,2,1:ndim, 1:ndim), Dwi(i, 1:ndim, 2) )
       !enddo

       call EvalVectorD(elem, -R_s(1:Qdof,1:nbDim,1:ndim), dofA, elem%vec(rhs, 1:ndim*dofA ) )

       !!deallocate(KsK)
       deallocate(R_s)
    endif

    deallocate(wi, Dwi, Re_1 )

  end subroutine ElementViscousVolumes


  !> Computing of "viscous" edge integrals for inner edges
  subroutine ElementViscousInnerEdge(elem, elem1, ie, Set_R_s, Set_K_sk)
    type(element), intent(inout):: elem, elem1  ! elem = element, elem1 = neigh element
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    ! Compute viscous fluxes and matrices K_sk(w), s,k=1,2
    interface
      subroutine Set_R_s(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
      end subroutine Set_R_s
      subroutine Set_K_sk(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, K_sk, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,ndimL), intent(inout) :: K_sk
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
       end subroutine Set_K_sk
    end interface
    !!!real, intent(in) :: Re_1       ! inverse of Reynolds number
    real, dimension(:,:), allocatable :: wi          ! w recomputed  in integ nodes
    real, dimension(:,:), allocatable :: Re_1    ! inverse of Reynolds number in integ nodes
    real, dimension(:,:,:,:,:), allocatable :: Ksk   ! matrices Ksk in integ nodes
    real, dimension(:), allocatable :: penal         ! local store arrays
    real, dimension(:,:), allocatable :: wii  ! w recomputed  in integration nodes
    real, dimension(:,:,:), allocatable :: Dwi  ! Dw recomputed  in integration nodes
    real, dimension(:,:,:), allocatable :: R_s  ! R_s(w) recomputed  in integration nodes
    real, dimension(:,:), allocatable :: Rflux  ! viscous flux
    real, dimension(:), allocatable :: ident         ! local store arrays

    real :: val, param
    real, dimension(1:ndim) :: qL, qR

    integer ::  dof, dofA, dof1, ie1, Qnum, Qdof, ndof, ndof1
    integer :: k, row, row1, i

    dof = elem%dof
    dofA = dof
    if(elem%deg_plus) dofA = elem%dof_plus

    dof1 = elem1%dof

    ndof = ndim*dof
    ndof1 = ndim*dof1

    ie1 = elem%face(nei_i,ie)

    !! seting of degree of the Gauss quadrature
    Qdof = elem%face(fGdof,ie)
    if(Qdof .ne. elem1%face(fGdof,ie1) ) then
       print*, 'Troubles in degrees of Gauss formulaes', elem%i, elem1%i
       print*,'##', Qdof, elem1%face(fGdof,ie1)
       print*,'##',elem%xc(:)
       print*,'##',elem1%xc(:)
       stop
    endif

    ! Reynolds number in integ nodes.
    allocate(Re_1(1:iRe, 1:Qdof) )
    if(state%model%Re > 0.) then
       Re_1(1,1:Qdof) = 1./state%model%Re
    else
       call Eval_func_Edge(elem, ie, elem%vec(aRe, 1:dof), Re_1(1,1:Qdof), .false. )
    endif

    ! interior penalty
    allocate(penal(1:Qdof) )
    !penal(1:Qdof) = state%space%sigma*Re_1(1:Qdof)  !!  elem%dn(ie)/elem%dn(ie)= 1
    ! VD sigma
    penal(1:Qdof) = Re_1(1, 1:Qdof) * max(elem%d_gamma , elem1%d_gamma)  !!  elem%dn(ie)/elem%dn(ie)= 1

    !penal(1:Qdof) = state%space%sigma*Re_1/elem%dn(ie)**(state%space%pen_deg-1)

    !if(state%modelName == 'scalar' .and.  &     ! POROUS MEDIA FLOW
    !     (state%model%idiff == 13  ) ) then
    Re_1(2:iRe, 1:Qdof) = transpose( elem%xi(ie, 1:Qdof, 2+1:2+iRe-1) ) ! setting the precomputed values
    !endif


    if(state%nlSolver%implicitly) then
       allocate( wi(1:Qdof,1:ndim), Dwi(1:Qdof,1:ndim,1:nbDim) )
       call Eval_w_Edge(elem, ie,  wi, .false.)

       if(state%modelName == 'scalar' .or. state%modelName == '2eqs' .or. state%modelName == 'porous') &
            call Eval_Dw_Edge(elem,  ie,  Dwi(1:Qdof,1:ndim,1:nbDim),  .false.)

       allocate(Ksk(1:Qdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim) )
       call Set_K_sk(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof,1:ndim,1:nbDim), Re_1, &
            Ksk(1:Qdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim), elem%xi(ie, 1:Qdof, 1:nbDim) )

       ! if(state%time%iter >= 4) then
       !     do k=1, Qdof
       !        write(65,'(30es12.4)')  &
       !             elem%xi(ie, k, 1:2), Ksk(k, 1, 1, 1, 1), Ksk(k, 2, 2, 1, 1), wi(k, 1),&
       !             Re_1(2:iRe, k)
       !     enddo
       !     if(elem%i == grid%nelem) stop '9ue93j......39u393i'
       ! endif


      !!if(.not.  state%only_diag)  then ! ORIGINAL SCHEME
        Ksk(:,:,:,:,:) = - Ksk(:,:,:,:,:)/2.    ! 1/2 in < . > operator
       !!else
       !!   Ksk(:,:,:,:,:) = - Ksk(:,:,:,:,:)
       !!endif


       ! ORIGINAL TERMS from Green's formula
       ! evaluation of matrix terms - diagonal block

       call EvalEdgeBlockDB(elem, ie, elem, ie, &
            Ksk(1:Qdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim), &
            elem%block(0)%Mb(1:ndof, 1:ndof) )

       ! penalty
       call EvalEdgeBlockDiagBB(elem, ie, elem, ie, penal(1:Qdof), elem%block(0)%Mb(1:ndof, 1:ndof)  )

       ! penalty
       if(.not.  state%only_diag) &
       call EvalEdgeBlockDiagBB(elem, ie, elem1,ie1,-penal(1:Qdof), elem%block(ie)%Mb(1:ndof, 1:ndof1) )

       if(state%space%m_IPG .ne. 0) then  ! not for IIPG
          Ksk(:,:,:,:,:) = -state%space%m_IPG* Ksk(:,:,:,:,:)    ! SIPG, NIPG, IIPG

          ! diagonal block,  transposition already included
          call EvalEdgeBlockBD(elem, ie, elem, ie, &
               Ksk(1:Qdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim), &
               elem%block(0)%Mb(1:ndof, 1:ndof))

          !!off diagonal block, transposition already included
          if(.not.  state%only_diag) &
               call EvalEdgeBlockBD(elem, ie, elem1, ie1, &
               Ksk(1:Qdof, 1:nbDim, 1:nbDim,1:ndim, 1:ndim),&
               elem%block(ie)%Mb(1:ndof, 1:ndof1) )

       endif

       ! ORIGINAL TERMS from Green's formula
       ! evaluation of matrix terms - off-diagonal block

       if(.not.  state%only_diag)  then
          ! evaluation of wi in integ. nodes from opposite side
          call Eval_w_Edge(elem1, ie1, wi, .true.)

          if(state%modelName == 'scalar' .or.state%modelName == '2eqs') &
               call Eval_Dw_Edge(elem1,  ie1,  Dwi(1:Qdof,1:ndim,1:nbDim),  .true.)


          ! setting the precomputed values
          Re_1(2:iRe, 1:Qdof) = transpose( elem1%xi(ie1, 1:Qdof, 2+1:2+iRe-1) )

          call Set_K_sk(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof,1:ndim,1:nbDim), Re_1, &
               Ksk(1:Qdof, 1:nbDim, 1:nbDim,1:ndim, 1:ndim),  elem1%xi(ie1, 1:Qdof, 1:nbDim))

          Ksk(:,:,:,:,:) = - Ksk(:,:,:,:,:)/2.    ! 1/2 in < . > operator

          call EvalEdgeBlockDB(elem, ie, elem1, ie1, Ksk(1:Qdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim), &
               elem%block(ie)%Mb(1:ndof, 1:ndof1) )
       endif


      deallocate(wi, Dwi, Ksk)

    else ! explicit terms
       allocate( wi(1:Qdof,1:ndim), wii(1:Qdof,1:ndim),  Dwi(1:Qdof,1:ndim,1:nbDim) )
       allocate(R_s(1:Qdof,1:nbDim, 1:ndim), Rflux(1:Qdof, 1:ndim), ident(1:Qdof))
       ident(:) = 1.


       ! compute < R_s(w) Dw > n_s

       ! from opposite side
       call Eval_w_Edge(elem1, ie1,  wii, .true.)
       call Eval_Dw_Edge(elem1,  ie1,  Dwi(1:Qdof,1:ndim,1:nbDim),  .true.)

       !if(state%local_problem)  then
       !   call Eval_w_Edge(elem1, ie1,  wii, .false.)
       !   call Eval_Dw_Edge(elem1,  ie1,  Dwi(1:Qdof,1:ndim,1:nbDim),  .false.)
       !endif


       call Set_R_s(ndim, nbDim, iRe, Qdof, wii(1:Qdof,1:ndim), Dwi(1:Qdof, 1:ndim, 1:nbDim), Re_1, &
            R_s(1:Qdof,1:nbDim, 1:ndim),  elem1%xi(ie1, 1:Qdof, 1:nbDim) )

       Rflux(1:Qdof, 1:ndim) = &
            R_s(1:Qdof, 1, 1:ndim) * elem%n(ie,1) + R_s(1:Qdof, 2, 1:ndim) * elem%n(ie,2)

       ! from actual side
       call Eval_w_Edge(elem, ie,  wi, .false.)
       call Eval_Dw_Edge(elem,  ie,  Dwi(1:Qdof,1:ndim,1:nbDim),  .false.)

       call Set_R_s(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof, 1:ndim, 1:nbDim), Re_1, &
            R_s(1:Qdof,1:nbDim, 1:ndim), elem%xi(ie, 1:Qdof, 1:nbDim))


       !if(.not. state%local_problem) then  ! ORIGINAL SCHEME
       Rflux(1:Qdof, 1:ndim) = Rflux(1:Qdof, 1:ndim) &
            +R_s(1:Qdof, 1, 1:ndim) * elem%n(ie,1) + R_s(1:Qdof, 2, 1:ndim) * elem%n(ie,2)
       Rflux(1:Qdof, 1:ndim) = Rflux(1:Qdof, 1:ndim) / 2.  ! 1/2 in < . > operator

       !else  !!solution of the local problem
       !   Rflux(1:Qdof, 1:ndim) = &
       !        R_s(1:Qdof, 1, 1:ndim) * elem%n(ie,1) + R_s(1:Qdof, 2, 1:ndim) * elem%n(ie,2)
       !endif


       ! alternative version
       !call Eval_aver_R_s_Edge(Set_R_s, elem, ie, Rflux)

       ! term < K(w) Dw> n \phi = < R(w, Dw) > n \phi
       call EvalEdgeVectorB(elem, ie, ident(1:Qdof),  Rflux(1:Qdof, 1:ndim), &
            dofA, elem%vec(rhs, 1:ndim*dofA ) )

       ! interior penalty, jump of the solution
       wi(1:Qdof, 1:ndim) = wii(1:Qdof, 1:ndim) - wi(1:Qdof, 1:ndim)

       !!!write(*,'(a8,2i5, 20es12.4)') 'IntPen:', elem%i, elem1%i, wi(1:Qdof, 1:ndim)

       call EvalEdgeVectorB(elem, ie, penal(1:Qdof), wi(1:Qdof, 1:ndim), &
            dofA, elem%vec(rhs,1:ndim*dofA) )


       ! STABILIZATION TERMS  from Green's formula
       if(state%space%m_IPG .ne. 0) then  ! not for IIPG
          if(ndim > 1 .and. elem%i == 1) then
             print*,'Does the array wi contains the correect values??'
             print*,'May be the following :::', '38ur93jd3woidsw'
          endif
          !call Eval_w_Edge(elem, ie,  wi, .false.)  ! Is this command correct?

          ! setting the precomputed values
          Re_1(2:iRe, 1:Qdof) = transpose( elem%xi(ie, 1:Qdof, 2+1:2+iRe-1) )

          allocate(Ksk(1:Qdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim) )
          call Set_K_sk(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof,1:ndim,1:nbDim), Re_1, &
               Ksk(1:Qdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim),  elem%xi(ie, 1:Qdof, 1:nbDim) )

          Ksk(:,:,:,:,:) = state%space%m_IPG* Ksk(:,:,:,:,:)/2.    ! SIPG, NIPG, IIPG

          call EvalEdgeVectorD(elem, ie, Ksk(1:Qdof,1:nbDim,1:nbDim,1:ndim, 1:ndim), &
               wi(1:Qdof, 1:ndim), dofA,  elem%vec(rhs, 1:ndim*dofA ) )
          deallocate(Ksk)
       endif

       deallocate(R_s, Dwi, wi, wii, ident)

    endif   ! explicitly, implicitly

    deallocate(penal, Re_1)

  end subroutine ElementViscousInnerEdge

  !> Computing of ONLY PENALTY integrals for inner edges (for stabilization)
  subroutine ElementViscousInnerPenalty(elem, elem1, ie, Set_R_s, Set_K_sk)
    type(element), intent(inout):: elem, elem1  ! elem = element, elem1 = neigh element
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    ! Compute viscous fluxes and matrices K_sk(w), s,k=1,2
    interface
      subroutine Set_R_s(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
      end subroutine Set_R_s
      subroutine Set_K_sk(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, K_sk, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,ndimL), intent(inout) :: K_sk
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
       end subroutine Set_K_sk
    end interface
    !!!real, intent(in) :: Re_1       ! inverse of Reynolds number
    real, dimension(:,:), allocatable :: wi          ! w recomputed  in integ nodes
    real, dimension(:,:), allocatable :: Re_1    ! inverse of Reynolds number in integ nodes
    real, dimension(:), allocatable :: penal         ! local store arrays
    real, dimension(:,:), allocatable :: wii  ! w recomputed  in integration nodes

    real :: val, param
    real, dimension(1:ndim) :: qL, qR

    integer ::  dof, dofA, dof1, ie1, Qnum, Qdof, ndof
    integer :: k, row, row1, i, ndof1

    dof = elem%dof
    dofA = dof
    if(elem%deg_plus) dofA = elem%dof_plus

    ndof  = ndim * dof
    dof1  = elem1%dof
    ndof1 = ndim * dof1

    ie1 = elem%face(nei_i,ie)

    !! seting of degree of the Gauss quadrature
    Qdof = elem%face(fGdof,ie)
    if(Qdof .ne. elem1%face(fGdof,ie1) ) then
       print*, 'Troubles in degrees of Gauss formulaes'
       stop
    endif

    ! Reynolds number in integ nodes.
    allocate(Re_1(1:iRe, 1:Qdof) )
    if(state%model%Re > 0.) then
       Re_1(1, 1:Qdof) = 1./state%model%Re
    else
       call Eval_func_Edge(elem, ie, elem%vec(aRe, 1:dof), Re_1(1,1:Qdof), .false. )
    endif

    ! interior penalty
    allocate(penal(1:Qdof) )
    !penal(1:Qdof) = state%space%sigma*Re_1(1:Qdof)  !!  elem%dn(ie)/elem%dn(ie)= 1
    ! VD sigma
    penal(1:Qdof) = Re_1(1,1:Qdof) * max(elem%d_gamma , elem1%d_gamma)  !!  elem%dn(ie)/elem%dn(ie)= 1

    !penal(1:Qdof) = state%space%sigma*Re_1/elem%dn(ie)**(state%space%pen_deg-1)



    if(state%nlSolver%implicitly) then
       ! penalty
       call EvalEdgeBlockDiagBB(elem, ie, elem, ie, penal(1:Qdof), elem%block(0)%Mb(1:ndof, 1:ndof)  )

       ! penalty
       if(.not.  state%only_diag) &
            call EvalEdgeBlockDiagBB(elem,ie,elem1,ie1,-penal(1:Qdof),elem%block(ie)%Mb(1:ndof,1:ndof1))


    else ! explicit terms
       allocate( wi(1:Qdof,1:ndim), wii(1:Qdof,1:ndim) )

       ! compute < R_s(w) Dw > n_s
       ! from opposite side
       call Eval_w_Edge(elem1, ie1,  wii, .true.)

       ! from actual side
       call Eval_w_Edge(elem, ie,  wi, .false.)

       ! interior penalty, jump of the solution
       wi(1:Qdof, 1:ndim) = wii(1:Qdof, 1:ndim) - wi(1:Qdof, 1:ndim)

       call EvalEdgeVectorB(elem, ie, penal(1:Qdof), wi(1:Qdof, 1:ndim), &
            dofA, elem%vec(rhs,1:ndim*dofA) )

       deallocate(wi, wii)

    endif   ! explicitly, implicitly

    deallocate(penal, Re_1)

  end subroutine ElementViscousInnerPenalty

  !> Computing of "viscous" edge integrals for boundary edges
  subroutine ElementViscousBoundEdge(elem,  ie, Set_R_s, Set_K_sk)
    type(element), intent(inout):: elem        ! elem = element
    integer, intent (in) :: ie                 ! inner index of edge, 1,2,3, (4)
    ! Compute viscous fluxes and matrices K_sk(w), s,k=1,2
    interface
      subroutine Set_R_s(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
      end subroutine Set_R_s
      subroutine Set_K_sk(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, K_sk, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,ndimL), intent(inout) :: K_sk
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
       end subroutine Set_K_sk
    end interface
    !real, intent(in) :: Re_1                   ! inverse of Reynolds number
    real, dimension(:,:), allocatable :: wi   ! w recomputed  in integ nodes
    real, dimension(:,:), allocatable :: Re_1    ! inverse of Reynolds number in integ nodes
    real, dimension(:,:), allocatable :: wB, w_BC, wD      ! wB recomputed  in integ nodes
    real, dimension(:,:,:,:,:), allocatable :: Ksk  ! matrices Ksk in integ nodes
    real, dimension(:), allocatable :: penal        ! local store arrays
    real, dimension(:), allocatable :: dnc    ! increment of face in integ nodes
    real, dimension(:,:,:), allocatable :: Dwi  ! Dw recomputed  in integration nodes
    real, dimension(:,:,:), allocatable :: R_s  ! R_s(w) recomputed  in integration nodes
    real, dimension(:,:), allocatable :: Rflux  ! viscous flux
    real, dimension(:), allocatable :: ident         ! local store arrays
    real, dimension(1:nbDim) :: xi
    real ::  rho_theta

    integer :: dof, dofA, kst, kst1, Qnum, Qdof, ndimL, ndof
    integer :: l, k, i, k1, row, col, ibb

    dof = elem%dof
    dofA = dof
    if(elem%deg_plus) dofA = elem%dof_plus

    ndof = ndim * dof

    !! seting of degree of the Gauss quadrature
    Qdof = elem%face(fGdof,ie)

    ! Reynolds number in integ nodes
    allocate(Re_1(1:iRe,1:Qdof) )
    if(state%model%Re > 0.) then
       Re_1(1,1:Qdof) = 1./state%model%Re
    else
       call Eval_func_Edge(elem, ie, elem%vec(aRe, 1:dof), Re_1(1,1:Qdof), .false. )
    endif

    ! boundary penalty
    allocate(penal(1:Qdof) )
    !penal(1:Qdof) = state%space%sigma*Re_1(1:Qdof)      !!  elem%dn(ie)/elem%dn(ie)= 1
    ! VD sigma
    penal(1:Qdof) = Re_1(1,1:Qdof) * elem%d_gamma  !!  elem%dn(ie)/elem%dn(ie)= 1


    !if(state%modelName == 'scalar' .and.  &     ! POROUS MEDIA FLOW
    !     (state%model%idiff == 13  ) ) then
    Re_1(2:iRe,1:Qdof) = transpose( elem%xi(ie, 1:Qdof, 2+1:2+iRe-1) ) ! setting the precomputed values
    !endif

    ! stronger penalty for boundary layers !!!
    !if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
    !   penal(1:Qdof) = state%space%sigma*max(1E-5, state%model%Re1 ) !!  elem%dn(ie)/elem%dn(ie)= 1
    !endif

    ! evaluation of w_ie in integ nodes
    allocate(wi(1:Qdof,1:ndim) )
    call Eval_w_Edge(elem, ie, wi, .false.)


    !if(elem%iBC(ie) == 0 ) then
    !   do j=1,Qdof
    !      wi(j,4) = wi(j,4) - 0.5*(wi(j,2)*wi(j,2)+ wi(j,3)*wi(j,3))/wi(j,1)
    !   enddo
    !   wi(1:Qdof,2:3) = 0.          ! velocities are zero
    !endif

    ! setting of outer normals in integration nodes
    allocate(dnc(1:Qdof) )
    if(elem%ibcur > 0) then
       dnc(1:Qdof) = elem%dnc(1:Qdof)
    else
     ! no courved edge on inlet Outlet
       dnc(1:Qdof) = elem%dn(ie)
    endif


    ! vector wB for BOUNDARY conditions
    allocate(wB(1:Qdof, 1:ndim), w_BC(1:Qdof, 1:ndim), wD(1:Qdof,1:ndim) )

    if(state%local_problem)  then ! local problem for EE, BC given apriori in local_problem2.f90
       wB(1:Qdof, 1:ndim) =  elem%wSS(ie, 1:Qdof, 1:ndim)

    elseif(state%modelName == 'NSe') then    ! Navier-Stokes equations
       if( state%type_IC == 9 ) then ! exact Steady state solution
          do k=1,Qdof
             xi(1:nbDim) = grid%b_edge(-elem%face(neigh,ie))%x_div(k, 1:nbDim)
             call Exact_SteadyChannel(xi(1:nbDim), wB(k,1:ndim) )
             !write(*,'(3e12.4,a5)') xi(1:nbDim), wB(k,1:1), '  Ri:'
          enddo
       else

          if(elem%iBC(ie) == 0 ) then
             !print*,'FIXED WALLS ', ie
             wB(1:Qdof,1) = wi(1:Qdof,1)  ! density is extrapolated
             wB(1:Qdof,2:3) = 0.          ! velocities are zero
             !do j=1,Qdof
             ! extrapolate from values on  Gamma_W
             !wB(1:Qdof,4) = wi(1:Qdof,4) - 0.5*(wi(1:Qdof,2)**2 + wi(1:Qdof,3)**2)/wi(1:Qdof,1)
             wB(1:Qdof,4) = wi(1:Qdof,4)  ! energy is extrapolated (Neumann BC for 4th component)

             !other case (60)
             !w_BC(j,4) = wi(j,1)*theta_D(j)
             !enddo

             ! viscous shock-vortex interaction
             if(state%type_IC == 8) wB(1:Qdof,1:4) = wi(1:Qdof,1:4)

             !       else if (state%BC(elem%iBC(ie))%inout == 0) then ! only inlet
          else   ! inlet or outlet
             do k=1,ndim
                wB(1:Qdof,k) = state%BC(elem%iBC(ie))%ww(k)
             enddo
             ! new
             !call ReprepareBCCharacteristic(Qdof, ndim, wi(1:Qdof,1:ndim), &
             !     wB(1:Qdof,1:ndim), elem%n(ie,:), &
             !     state%BC(elem%iBC(ie))%press_extrap, elem%xc(:) )

             ! (exact) characteristic BC
             call SetBCexactRiemann(Qdof, ndim, wi(1:Qdof,1:ndim), &
                  wB(1:Qdof,1:ndim), w_BC(1:Qdof,1:ndim), elem%n(ie,:),  &
                  elem%xc(1:nbDim) )
             wB(1:Qdof,1:ndim) = w_BC(1:Qdof,1:ndim)

             ! old
             !do j=1,Qdof
             !   ! extrapolation from values on Gamma_I
             !   rho_theta = wi(j,4) - 0.5*(wi(j,2)*wi(j,2)+ wi(j,3)*wi(j,3) )/wi(j,1)
             !
             !   wB(j,4) = rho_theta + 0.5*(wB(j,2)*wB(j,2) + wB(j,3)*wB(j,3))/wB(j,1)
             !enddo
             !       else
             !          print*,'problem in inv_fluxes.f90'
             !          stop
             ! NOTHING on OUTLET
          endif
       endif

    ! elseif(ndim >= 5 .and. ndim <= 6) then    ! RANS
    !    if(elem%iBC(ie) == 0 ) then
    !       !print*,'FIXED WALLS  ', ie
    !       wB(1:Qdof,1) = wi(1:Qdof,1)  ! density is extrapolated
    !       wB(1:Qdof,2:3) = 0.          ! velocities are zero
    !       !do j=1,Qdof
    !       ! extrapolate from values on  Gamma_W
    !       !wB(1:Qdof,4) = wi(1:Qdof,4) - 0.5*(wi(1:Qdof,2)**2 + wi(1:Qdof,3)**2)/wi(1:Qdof,1)
    !       wB(1:Qdof,4) = wi(1:Qdof,4)  ! energy is extrapolated (Neumann BC for 4th component)

    !       ! turbulence extrapolation
    !       wB(1:Qdof,5:ndim) = wi(1:Qdof,5:ndim )
    !    else   ! inlet or outlet
    !       !print*,'INLET OUTLET ', ie
    !       do k=1,ndim
    !          wB(1:Qdof,k) = state%BC(elem%iBC(ie))%ww(k)
    !       enddo
    !       ! new
    !       !call ReprepareBCCharacteristic(Qdof, ndim, wi(1:Qdof,1:ndim), &
    !       !     wB(1:Qdof,1:ndim), elem%n(ie,:), &
    !       !     state%BC(elem%iBC(ie))%press_extrap, elem%xc(:) )

    !       ! (exact) characteristic BC
    !       ndimL = 4
    !       call SetBCexactRiemann(Qdof, ndimL, wi(1:Qdof,1:ndimL), &
    !            wB(1:Qdof,1:ndimL), w_BC(1:Qdof,1:ndimL), elem%n(ie,:),  &
    !            elem%xc(1:nbDim) )
    !       wB(1:Qdof,1:ndimL) = w_BC(1:Qdof,1:ndimL)

    !       ! turbulence - extrapolation
    !       wB(1:Qdof, 5:ndim) = wi(1:Qdof,5:ndim)
    !    endif

    elseif(state%modelName == 'scalar' .or.state%modelName == '2eqs' ) then   ! BC for scalar equation
       if(elem%iBC(ie) == 0 ) then !   !print*,'FIXED WALLS  ', ie, NEUMANN BOUNDARY CONDITION
          wB(1:Qdof,1:ndim) = wi(1:Qdof,1:ndim)  ! density is extrapolated

       else ! Dirichlet boundary condition
          if(state%homogenDirichlet) then
             wB(:,:) = 0.

          else
             do k=1,Qdof
                xi(1:nbDim) = grid%b_edge(-elem%face(neigh,ie))%x_div(k, 1:nbDim)
                call Exact_Scalar(xi(1:nbDim), wB(k,1:ndim), state%time%ctime )
                !write(42,'(30e12.4)') xi(1:2), wB(k,1)
             enddo

          endif

       endif

    elseif(state%modelName == 'wet_steam') then    ! wet_steam equations

          if(elem%iBC(ie) == 0 ) then
             !print*,'FIXED WALLS ', ie
             wB(1:Qdof,1) = wi(1:Qdof,1)  ! density is extrapolated
             wB(1:Qdof,2:3) = 0.          ! velocities are zero
             !do j=1,Qdof
             ! extrapolate from values on  Gamma_W
             !wB(1:Qdof,4) = wi(1:Qdof,4) - 0.5*(wi(1:Qdof,2)**2 + wi(1:Qdof,3)**2)/wi(1:Qdof,1)
             wB(1:Qdof,4) = wi(1:Qdof,4)  ! energy is extrapolated (Neumann BC for 4th component)

             ! viscous shock-vortex interaction
             if(state%type_IC == 8) wB(1:Qdof,1:4) = wi(1:Qdof,1:4)
             !       else if (state%BC(elem%iBC(ie))%inout == 0) then ! only inlet


             ! wet steam part
             wB(1:Qdof,5:8) = wi(1:Qdof,5:8)  ! do nothing

          ! Same as in ElementInviscidIOEdge !??
          else   ! inlet or outlet
             do k=1,ndim
                wB(1:Qdof,k) = state%BC(elem%iBC(ie))%ww(k)
             enddo

             ! (exact) characteristic BC
             call SetBCexactRiemann(Qdof, 4, wi(1:Qdof,1:4), &
                  wB(1:Qdof,1:4), w_BC(1:Qdof,1:4), elem%n(ie,:),  &
                  elem%xc(1:nbDim) )
             wB(1:Qdof,1:4) = w_BC(1:Qdof,1:4)

             ibb = state%BC(elem%iBC(ie))%inout
             if(ibb == 0) then  ! inlet
               wB(1:Qdof,5:8) = wD(1:Qdof,5:8)  ! data from wet.ini
             elseif(ibb == 1) then  ! outlet
               wB(1:Qdof,5:8) = wi(1:Qdof,5:8)  ! do nothing
             endif

          endif

    elseif(state%modelName == 'porous' ) then   ! BC for porous media flow
       if(elem%iBC(ie) == 0 ) then !   !print*,'FIXED WALLS  ', ie, NEUMANN BOUNDARY CONDITION
          wB(1:Qdof,1:ndim) = wi(1:Qdof,1:ndim)  ! density is extrapolated

       else ! Dirichlet boundary condition
          if(state%homogenDirichlet) then
             wB(:,:) = 0.

          else

             call Exact_Porous(Qdof, grid%b_edge(-elem%face(neigh,ie))%x_div(1:Qdof, 1:nbDim), &
                  wB(1:Qdof,1:ndim), state%time%ctime )

             !do k=1,Qdof
             !   xi(1:nbDim) = grid%b_edge(-elem%face(neigh,ie))%x_div(k, 1:nbDim)
             !   !call Exact_Scalar(xi(1:nbDim), wB(k,1:ndim), state%time%ctime )
             !   write(22,'(30e12.4)') xi(1:2), wB(k,1)
             !enddo

          endif

       endif
    else
       print*,'Viscous BC not implemented for ndim=',ndim
    endif

    ! fixed walls
    !if(elem%iBC(ie) == 0 ) penal(:) = 1. * penal(:)

    !elem%block(0)%Mb(:,:) = 0.
    !elem%vec(rhs, :) = 0.

    !if(elem%i == 1) write(*,'(a6,3e12.4,i5)') '!!!???', elem%xc(1:nbDim), wB(1,1:1), elem%i


    if(state%nlSolver%implicitly) then
       allocate( Dwi(1:Qdof,1:ndim,1:nbDim) )

       if(state%modelName == 'scalar' .or. state%modelName == '2eqs' .or. state%modelName == 'porous') &
            call Eval_Dw_Edge(elem,  ie,  Dwi(1:Qdof,1:ndim,1:nbDim),  .false.)

       allocate(Ksk(1:Qdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim) )
       call Set_K_sk(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof,1:ndim,1:nbDim), Re_1,&
            Ksk(1:Qdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim), elem%xi(ie, 1:Qdof, 1:nbDim) )
       Ksk(:,:,:,:,:) = - Ksk(:,:,:,:,:)

       !SC   if(state%model%Re < 0. .or. state%ST_Vc <  -1.E-05) then ! only for stabilization
       !SC      ! ORIGINAL TERMS from Green's formula, diagonal block
       !SC   call EvalEdgeBlockDB(elem, ie, elem, ie, Ksk(1:Qdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim),&
       !SC       elem%block(0)%Mb(1:ndim*dof, 1:ndim*dof) )
       !SC   else

       ! fixed walls
       if(ndim >= 4 .and. elem%iBC(ie) == 0 ) &
            Ksk(1:Qdof, 1:nbDim, 1:nbDim, 4:ndim , 1:ndim) = 0.

       ! ORIGINAL TERMS from Green's formula, diagonal block
       call EvalEdgeBlockDB(elem, ie, elem, ie, &
            Ksk(1:Qdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim), &
            elem%block(0)%Mb(1:ndim*dof, 1:ndim*dof) )

       if(ndim >= 4 .and. elem%iBC(ie) == 0) then  !!! fixed walls
          ! penalty, EvalEdgeVectorB not called since v = 0
          call EvalEdgeBlockDiagBB23(elem, ie, elem, ie, penal(1:Qdof), elem%block(0)%Mb(1:ndof, 1:ndof))

       else
          ! penalty
          call EvalEdgeBlockDiagBB(elem, ie, elem, ie, penal(1:Qdof), elem%block(0)%Mb(1:ndof, 1:ndof) )

          !write(*,'(a4,40es11.3)') 'e0: ', grid%elem(1)%vec(rhs,:)

          ! boundary penalty
          call EvalEdgeVectorB(elem, ie, penal(1:Qdof), wB(1:Qdof, 1:ndim), &
               dof, elem%vec(rhs,1:ndim*dof) )

          !write(*,'(a4,40es11.3)') 'e1: ', grid%elem(1)%vec(rhs,:)

       endif
       !SC    endif

       ! STABILIZATION TERMS  from Green's formula
       if(state%space%m_IPG .ne. 0) then ! not IIPG
          Ksk(:,:,:,:,:) = -state%space%m_IPG* Ksk(:,:,:,:,:)    ! SIPG, NIPG, IIPG

          ! evaluation of matrix terms - diagonal block
          call EvalEdgeBlockBD(elem, ie, elem, ie, &
               Ksk(1:Qdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim), &
               elem%block(0)%Mb(1:ndim*dof, 1:ndim*dof) )

          ! evaluating of the Vector
          call EvalEdgeVectorD(elem, ie, Ksk(1:Qdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim), &
               wB(1:Qdof, 1:ndim),  dof, elem%vec(rhs,1:ndim*dof) )
       endif

       deallocate(Dwi, Ksk)

    else ! explicit discretization

       allocate( R_s(1:Qdof,1:nbDim, 1:ndim), Rflux(1:Qdof, 1:ndim), ident(1:Qdof))
       ident(:) = 1.

       !if(elem%i == 1) write(*,'(a5, 8e12.4)') 'Mvec', elem%vec(rhs,:)

       ! evaluation of <R(w, DW) n>
       allocate( Dwi(1:Qdof,1:ndim,1:nbDim) )
       call Eval_Dw_Edge(elem,  ie,  Dwi(1:Qdof,1:ndim,1:nbDim),  .false.)

       call Set_R_s(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof, 1:ndim, 1:nbDim), Re_1, &
            R_s(1:Qdof,1:nbDim, 1:ndim) , elem%xi(ie, 1:Qdof, 1:nbDim) )

       ! fixed walls for the Navier-Stokes equations
       if(ndim >=4 .and. ndim <= 6 .and. elem%iBC(ie) == 0 )&
            R_s(1:Qdof, 1:nbDim, 4:ndim) = 0.

       if(elem%ibcur > 0) then  ! for curved element with non-constant outer normals
          do i=1,Qdof
             Rflux(i, 1:ndim) = &
                  R_s(i, 1, 1:ndim) * elem%nc(i,1) + R_s(i, 2, 1:ndim) * elem%nc(i,2)
          enddo
       else
          Rflux(1:Qdof, 1:ndim) = &
               R_s(1:Qdof, 1, 1:ndim) *elem%n(ie,1) + R_s(1:Qdof, 2, 1:ndim) *elem%n(ie,2)
       endif

       !if(elem%i == 2) write(*,'(a6,2i5,20es14.6)') 'RfluX=',elem%i,ie ,Rflux(1:Qdof, 1:ndim)

       !if(elem%iBC(ie) == 0 ) then
       !   write(*,'(a6,i5,30es12.4)') 'resVV:',elem%i,elem%vec(rhs, :)
       !   write(*,'(a6,2i5,20es14.6)')' fluX=',elem%i,ie ,wi(1:Qdof, 1:ndim)
       !endif

       call EvalEdgeVectorB(elem, ie, ident(1:Qdof),  Rflux(1:Qdof, 1:ndim),&
            dofA, elem%vec(rhs, 1:ndim*dofA ) )

       ! STABILIZATION TERMS  from Green's formula
       if(state%space%m_IPG .ne. 0) then  ! not for IIPG
          allocate(Ksk(1:Qdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim) )
          call Set_K_sk(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwi(1:Qdof,1:ndim,1:nbDim), Re_1, &
               Ksk(1:Qdof, 1:nbDim, 1:nbDim, 1:ndim, 1:ndim), elem%xi(ie, 1:Qdof, 1:nbDim) )
          Ksk(:,:,:,:,:) = state%space%m_IPG* Ksk(:,:,:,:,:)    ! SIPG, NIPG, IIPG
       endif

       ! boundary penalty, jump of the solution
       wi(1:Qdof, 1:ndim) = wB(1:Qdof, 1:ndim) - wi(1:Qdof, 1:ndim)

       ! no penalty for density and energy
       if(ndim >= 4 .and. ndim <= 6 .and. elem%iBC(ie) == 0 ) then
          wi(1:Qdof, 1) = 0.
          wi(1:Qdof, nbDim+2:ndim) = 0.
       endif

       !if(elem%iBC(ie) == 0 ) then
       !   write(*,'(a6,i5,30es12.4)') 'resVV:',elem%i,elem%vec(rhs, :)
       !   write(*,'(a6,2i5,20es14.6)')'JfluX=',elem%i,ie ,wi(1:Qdof, 1:ndim)
       !endif


       call EvalEdgeVectorB(elem, ie, penal(1:Qdof), wi(1:Qdof, 1:ndim), &
            dofA, elem%vec(rhs,1:ndim*dofA) )

       ! STABILIZATION TERMS  from Green's formula
       if(state%space%m_IPG .ne. 0) then  ! not for IIPG

          call EvalEdgeVectorD(elem, ie, Ksk(1:Qdof,1:nbDim,1:nbDim,1:ndim, 1:ndim), &
               wi(1:Qdof, 1:ndim), dofA, elem%vec(rhs, 1:ndim*dofA ) )

          deallocate(Ksk)
       endif

       deallocate(R_s, Dwi, ident)

    endif   ! explicitly, implicitly

    deallocate(wB, w_BC, wi, dnc)
    deallocate(penal, Re_1)

  end subroutine ElementViscousBoundEdge





  !> initiate w(*,i):= w(0,i)
  subroutine InitElementW (elem )
    type(element), intent(inout) :: elem
    integer :: k, nlev,dof,Tdof

!
!    if (state%SP) then
!	if(state%modelName == 'incNS') then
!             !do nothing ???
!	endif
!    else
    if(state%time%disc_time == 'STDG') then  ! ST DGM
       !used only for initial condition and when the mesh is recomputed
       dof = elem%dof
       Tdof = elem%Tdof

       do k=1,ndim
          elem%wST(k, 1:dof, 1) = elem%w(0, (k-1)*dof+1:k*dof )
       enddo
       elem%wST(1:ndim, 1:dof, 2:Tdof) = 0.

       ! setting of the initial time level
       call Eval_wSTfin_Elem( elem )

    else
       nlev = state%time%deg+1

       do k=1, nlev
          elem%w(k,:) = elem%w(0,:)
       enddo

    endif

  end subroutine InitElementW

  !> compute the term (matrix in front of the timederivative term
  subroutine Compute_Time_deriv(elem, deg_plus,  sdim, T_mat )
    type(element), intent (inout) :: elem
    logical, intent(in) :: deg_plus
    integer, intent(in) :: sdim    ! size of the space block
    real, dimension(1:sdim, 1:sdim), intent(inout) :: T_mat  ! resulting space block

    type(volume_rule), pointer :: V_rule
    type(Lagrang_rule), pointer :: L_rule
    real, dimension(:,:), pointer :: phi

    real, dimension(:,:), allocatable :: wi ! w recomputed  in integ nodes
    real, dimension(:,:), allocatable :: wR ! reconstructed wR
    real, dimension(:), allocatable :: wbar ! average of w
    real, dimension(:,:,:), allocatable :: TA ! matrix
    real, dimension(:,:,:), allocatable :: TAp ! matrix
    real, dimension(:,:), allocatable :: Tf ! matrix
    !real, dimension(:,:), allocatable :: vals
    real :: param, rK, alpha
    integer ::  Qdof, dof, dofA
    integer :: i, l, deg1, Ldof, k, ist


    dof = elem%dof
    dofA = dof
    !if(elem%deg_plus) dofA = elem%dof_plus
    if(deg_plus) dofA = elem%dof_plus

    Qdof = elem%Qdof
    V_rule => state%space%V_rule(elem%Qnum)
    !!Qdof = V_rule%Qdof

    ! setting of the state vector in integration nodes
    allocate(wi(1:Qdof,1:ndim) )
    call Eval_w_Elem(elem, wi(1:Qdof,1:ndim) )
    !!call Eval_w_Elem_plus(elem, V_rule, wi(1:Qdof,1:ndim) )

    ! setting of matrixes TA in integration nodes
    allocate(TA(1:Qdof, 1:ndim,1:ndim))

    ! constant matrix:
    TA=0.
    do i=1, ndim
       TA(:,i, i) = 1.
    enddo

    if(state%modelName == 'porous' ) then
       !allocate(wbar(1:ndim) )

       !if( state%nlSolver%implicitly ) then  ! avereging for the Jacobian
       !   alpha = 1.0
       !   call Eval_aver_w_Elem(elem, wbar)
       !   do i=1, ndim
       !      wi(1:Qdof, i) =  alpha *  wi(1:Qdof, i) + (1. - alpha) * wbar(i)
       !   enddo
       !endif
       !deallocate( wbar)

       !allocate( wR(1:ndim, 1:elem%dof_plus) )

       call Set_Time_Matrix_porous(elem, ndim,  Qdof, wi(1:Qdof, 1:ndim), elem%xi(0,1:Qdof, 1:2+iRe),&
            TA(1:Qdof, 1:ndim, 1:ndim) ) !, wR)

       !do l=1,Qdof
       !   write(22, '(23es12.4)' ) elem%xi(0,l,1:2), TA(l, 1,1)
       !   write(*,'(a6,i5, 200es12.4)') 'TA:',l,wi(l, 1:ndim), elem%xi(0,l, 1:2+iRe-1), TA(l,1,1)
       !enddo
       !print*,'--------33e3'


       !if( abs( elem%xc(1) -12.) < 3 .and. abs( elem%xc(2) -0.) < 3) &
       ! if( dot_product(elem%xc- grid%elem(1428)%xc , elem%xc- grid%elem(1428)%xc) < 2.) &
       !call PlotElemFunction3D(63, elem, elem%dof_plus, wR(1, 1:elem%dof_plus) )

       !deallocate(wR)

       !if(elem%i == grid%nelem)  stop "9eu393id3"

       ! NEW VARIANT USING CONTINUOUS PROJECTION,
       !ITS DEGREE DEFINED IN  Setting_of_space_variable_coeffs in problem.f90
       deg1 = elem%Ldeg !+ 1
       L_rule => state%space%L_rule(deg1)
       Ldof = L_rule%Qdof

       ! w in Lagrang. integ nodes
       phi => L_rule%phi(1:dof, 1:Ldof)

       do k=1, ndim
          ist = (k-1)*dof + 1
          wi(1:Ldof, k) =  matmul( elem%w(0,ist:ist+dof-1), phi(1:dof, 1:Ldof))
       enddo

       ! setting of matrixes TA in Lagrang integration nodes
       allocate(TAp(1:Ldof, 1:ndim,1:ndim))

       call Set_Time_Matrix_porous(elem, ndim,  Ldof, wi(1:Ldof, 1:ndim), elem%Lxi(1:Ldof, 1:2+iRe), &
            TAp(1:Ldof, 1:ndim, 1:ndim))

       !do i=1,Ldof
       !   write(*,'(a6,i5, 200es12.4)') 'TMp:',i,wi(i, 1:ndim), vals(i, 1:2+iRe-1), TAp(i,1,1)
       !   write(23,*) Fxi(i, 1:2), TAp(i,1,1)
       !enddo
       !print*,'################################### u39j3o'


       ! recomputation to V_rule integ nodes
       call Lagr2QnodesDiff(elem, deg1,  TAp(1:Ldof, 1:ndim, 1), V_rule,  TA(1:Qdof, 1:ndim, 1) )

       !do i=1,Qdof
       !   write(*,'(a6,i5, 200es12.4)') 'TMp:',i,wi(i, 1:ndim), vals(i, 1:2+iRe-1), TAp(i,1,1)
       !   write(24,*) elem%xi(0,i, 1:2),vals(i, 1)
       !enddo

       !TA(1:Qdof, 1:ndim, 1) = vals(1:Qdof, 1:ndim)

       deallocate(TAp)
       !if(elem%i == grid%nelem) stop "u39udj3jd3odj3pod3o"

    endif


    ! evaluation of matrix terms
    T_mat = 0.

    call EvalBlockBB_dofA(elem,  TA(1:Qdof, 1:ndim, 1:ndim), dofA, T_mat(1:ndim*dofA, 1:ndim*dofA ) )

    deallocate(wi, TA)
    !deallocate(vals )

  end subroutine Compute_Time_deriv

  ! doesn't work if elem%dof (elem%Qdof) changes in time
  ! imp = state%nlSolver%implicitly - influences the sign of added part (.true. -> '-' , .false. -> '+')
  !> STDG: Evaluation of vector Q - \f$q^{i,j} \phi_m^{j}(t_{m-1}+) (w_{m-1}^{-},\varphi_m^i)\f$
  subroutine Elem_wSTfinToRhsST(elem , imp)
   type(element), intent (inout) :: elem
   logical, intent (in) :: imp
   class(Time_rule), pointer :: T_rule
   type(volume_rule), pointer  :: V_rule
   real, dimension(1:elem%Qdof) :: weights
   real, dimension(:,:), allocatable :: T_mat

!   real, dimension(:,:), pointer :: phi
   integer :: Qdof, Tdof, TQnum, dof, i, j, k, l
   integer :: wdof                   !dof/Tdof for the solution wST can be smaller than dof/Tdof for test func if elem%deg_plus
   real :: val, val_new

   Qdof = elem%Qdof
!   if (Qdof /= size(elem%wSTfin(1,:))) then
!      print*, 'Problem in Elem_wSTfinToRhsST', Qdof, size(elem%wSTfin(1,:))
!      stop
!   endif

   wdof = elem%dof
   !wTdof = elem%Tdof

   !NEW for estimates
   ! FR: Upraveno podle Vitkova pristupu v ComputeTerms OK??
   if(elem%deg_plus .and. .not. imp) then
      dof = elem%dof_plus
      Tdof = elem%Tdof_plus
   else
      dof = elem%dof
      Tdof = elem%Tdof
   endif


   if (Qdof == state%space%V_rule(elem%Qnum)%Qdof) then
      Qdof = elem%Qdof
   else
      stop
      print* , 'Problem in wSTfinToRhsST'
   endif

   V_rule => state%space%V_rule(elem%Qnum)
   T_rule => state%time%T_rule(elem%TQnum)

   associate ( phi => state%space%V_rule(elem%Qnum)%phi( 1:dof,1:Qdof)  )

!   phi => V_rule%phi  ! , &
!             )

   !, T_rule => state%time%T_rule(elem%TQnum), &
   !            phi => V_rule%phi( 1:dof,1:Qdof) )

   !   V_rule => state%space%V_rule(elem%Qnum)
   !   !FR problem when we change T_rule in ComputeSTDGM_Terms, maybe NOT because we are using only the value of phi(t_{m-1}^+) which is the same
   !   T_rule => state%time%T_rule(elem%TQnum)
   !
   !   phi => V_rule%phi( 1:dof,1:Qdof)

      !  print*, 'Elem_wSTfinToRhsST called', dof

      ! print*, 'eta = (1/tau)', 1./ state%time%tau(1)

      call Eval_V_Weights(elem, weights(1:Qdof) )

      !the sign +-( (1./ state%time%tau(1)) * T_rule%phi(1:Tdof,-1) * val_new ) depends on implictly=F/T
      if (imp) then

         stop "NEVER CALLED !!!!"
!         if(elem%deg_plus) then
!            !FERROR
!            !if(elem%i == 1) &
!            ! we never call implicitly=TRUE && deg_plus
!            !print*, 'inv_fluxes/wSTfinToRhsST called with implicitly=TRUE && deg_plus, Not used!'
!            !stop
!         endif

         do k = 1,ndim
            do i = 1, dof
   !!! Takhle nelze počítat - V_rule%weights se mění v závislosti na F: K_ref -> K
               ! val = elem%area * sum( V_rule%weights(1:Qdof) * elem%wSTfin(k,1:Qdof) * V_rule%phi(i,1:Qdof) )

               ! print*, 'wSTfin_ToRhsST val driv:',  val
               !val = sum( weights(1:Qdof) * elem%wSTfin(k,1:Qdof) * V_rule%phi(i,1:Qdof) )

               !NEW for adapt
               val_new = dot_product( weights(1:Qdof),  &
                    matmul(elem%wSTfinAD(k,1:dof), phi(1:dof, 1:Qdof)) * V_rule%phi(i,1:Qdof) )
               elem%rhsST(k,i,1:Tdof) = elem%rhsST(k,i,1:Tdof)  &
                    - ( (1./ state%time%tau(1)) * T_rule%phi(1:Tdof,-1) * val_new )
            enddo !i
         enddo !k
      !imp = FALSE
      else

         if( .not. state%model%varying_time_term) then ! original subroutine
            do k = 1,ndim
               do i = 1, dof
                  !val = sum( weights(1:Qdof) * elem%wSTfin(k,1:Qdof) * V_rule%phi(i,1:Qdof) )
                  val_new = dot_product( weights(1:Qdof),  &
                       matmul(elem%wSTfinAD(k,1:wdof), phi(1:wdof, 1:Qdof)) * V_rule%phi(i,1:Qdof) )
                  !if (abs(val - val_new) > 1E-15) then
                  !   print*, 'problem in whST to rhst' , val-val_new
                  !   stop
                  !endif
                  !               if ( .not. associated( elem%rhsST ) ) then
                  !                  stop 'not alloc'
                  !               else
                  !                  print*, 'size RHS:', size( elem%rhsST )
                  !                  print*, state%time%tau, val_new
                  !                  print*, T_rule%phi(1:Tdof,-1)
                  !               endif
                  elem%rhsST(k,i,1:Tdof) = elem%rhsST(k,i,1:Tdof)  &
                       + ( (1./ state%time%tau(1)) * T_rule%phi(1:Tdof,-1) * val_new )
               enddo !i
            enddo !k

         else !state%model%varying_time_term
            allocate( T_mat(1: ndim*dof, 1:ndim*dof) )

            ! setting of the solution from the t_{m-1}^-
            do k = 1,ndim
               elem%w(0, (k-1)*wdof +1 : k*wdof) = elem%wSTfinAD(k,1:wdof)
            enddo

            ! computing of the matrix in front of the time derivative term at t_{m-1}^+
            !call Compute_Time_deriv(elem, elem%deg_plus, ndim*dof, T_mat(1:ndim*dof, 1:ndim*dof) )

            ! Time Penalty Var1
            T_mat(1:ndim*dof,1:ndim*dof)  = 0.
            do l=1,ndim*dof
               T_mat(l,l) = 1.
            enddo

            ! setting of the terms
            do k = 1,ndim
               do l = 1, Tdof

                  val_new = T_rule%phi(l,-1) / state%time%tau(1)

                  do j=1, dof
                     elem%rhsST(k, j, l) = elem%rhsST(k, j, l)  +  val_new  &
                          * dot_product( T_mat( (k-1)* dof + j, (k-1)*dof+1:(k-1)*dof + wdof ), &
                          elem%wSTfinAD(k, 1:wdof) )
                  enddo
               enddo  ! l=1,Tdof

            enddo !k
            deallocate( T_mat)

         endif !  state%model%varying_time_term

      endif
   end associate


  end subroutine Elem_wSTfinToRhsST

  !> compute the maximal corresponsing eigenvalues of inviscid system
  subroutine InitElementEigenvals(elem )
    type(element):: elem
    real :: v2, p, c, maxLam
    integer :: dof, j

    !! seeking of maximal eigenvalue
!    if(state%SP) then ! saddle point
!       if(state%modelName == 'incNS') then
!          maxLam = 1.
!          if(elem%i ==1)  print*,'Check in InitElementEigenvals(elem )'
!       else
!          print*,'Init Element W not implemented for ',state%modelName
!          stop
!       endif
!
!    else

       if(state%modelName == 'NSe' .or. state%modelName == 'wet_steam') then ! Euler and Navier-Stokes equations and wet steam equations

          dof = elem%dof
          v2 = 0.
          do j=2,nbDim + 1
             v2 = v2 + (elem%w(0,(j-1)*dof+1)/elem%w(0,1))**2
          enddo
          p = state%model%kappa1 *(elem%w(0,3*dof+1) - 0.5*elem%w(0,1)*v2)
          c = (p*state%model%kappa/elem%w(0,1) )**0.5
          maxLam = (c + v2**0.5)*elem%diam/elem%area

       elseif( state%modelName == 'swe') then ! SWE flows

          if(elem%i == 1) print*,' Verify  InitElementEigenvals 49uj4jfd'

          dof = elem%dof
          v2 = 0.
          if(elem%w(0,1) > 0) then
             do j=2,nbDim + 1
                v2 = v2 + (elem%w(0,(j-1)*dof+1)/elem%w(0,1))**2
             enddo

             ! pressure
             p = state%model%p0 * elem%w(0,1)**state%model%kappa

             ! speed of sound
             c = (p*state%model%kappa/elem%w(0,1) )**0.5

             ! maximal eigenvalue
             maxLam = (c + v2**0.5)*elem%diam/elem%area
          else
             maxLam = 0.
          endif


       elseif(state%modelName == 'pedes' ) then ! pedestrian flows

          dof = elem%dof
          v2 = 0.
          if(elem%w(0,1) > 0) then
             do j=2,nbDim + 1
                v2 = v2 + (elem%w(0,(j-1)*dof+1)/elem%w(0,1))**2
             enddo

             ! pressure
             p = state%model%p0 * elem%w(0,1)**state%model%kappa

             ! speed of sound
             c = (p*state%model%kappa/elem%w(0,1) )**0.5

             ! maximal eigenvalue
             maxLam = (c + v2**0.5)*elem%diam/elem%area
          else
             maxLam = 0.
          endif


       elseif(state%modelName == 'scalar' .or. state%modelName == '2eqs') then ! scalar equations
          maxLam = elem%diam/elem%area   !! ???

       elseif(state%modelName == 'porous') then ! porous media flow
          maxLam = elem%diam/elem%area   !! ???

       else
          print*,'Init Element W (2) not implemented for model = ',state%modelName
          stop
       endif

!    endif

    !print*,'### inv_fl',elem%i, state%max_eigenvals, maxLam

    state%max_eigenvals = max(state%max_eigenvals, maxLam )

    !print*,'InitElementW',state%max_eigenvals
  end subroutine InitElementEigenvals

  subroutine ExactElementW(elem, time)
    type(element):: elem
    real,intent(in):: time(:)
    !real:: w_save(size(elem%w))
    real:: w_save(ndim * state%space%max_dof)
    integer:: k


    w_save(1:ndim * elem%dof) = elem%w(0,1:ndim * elem%dof)

    !print*,'#### CHECK in ExactElemW in inv_fluxes.f90:','size(elem%w,1) =',size(elem%w,1),'??'
    do k = 1,size(elem%w,1)-1
       state%time%ctime = time(k)
       call SetOneElementIC(elem, .false.)
       elem%w(k,1:ndim * elem%dof) = elem%w(0,1:ndim * elem%dof)
    end do

    elem%w(0,1:ndim * elem%dof) = w_save(1:ndim * elem%dof)

  end subroutine ExactElementW

! FR NOT USED ???
!  !> evaluate elem%block := elem%block + eta * MassMatrix
!  subroutine AddElementMassMatrix (elem, eta )
!    type(element):: elem
!    real, intent(in) :: eta
!    integer :: k, dof, il, it, dof1
!
!    dof = elem%dof
!
!
!    ! multiply by tau
!    !if(elem%i == 10) print*,'@@@',state%time%alpha(0)
!
!    ! Add mass matrix
!    !if(elem%F%iFlin) then
!    !   ! only for ORTHOGONAL
!    !   do k=1,ndim*dof
!    !      elem%block(0)%Mb(k, k) =  elem%block(0)%Mb(k, k) &
!    !           + state%time%alpha(0) * elem%area*(5-elem%type)
!    !   enddo
!    !else
!       do k=1,ndim
!          il = (k-1)*dof + 1   ! lower bound of block
!          it = k*dof           ! upper bound of block
!
!          elem%block(0)%Mb(il:it, il:it ) =  elem%block(0)%Mb(il:it, il:it) &
!               + state%time%alpha(0)*eta *elem%Mass%Mb(1:dof,1:dof)
!
!       enddo
!    !endif
!
!
!  end subroutine AddElementMassMatrix

  !> compute \f$ \int_K \hat{w}  \phi_i dx \f$
  subroutine SetElementMassVector(elem)
    type(element):: elem
    real, dimension(:,:), allocatable :: wi
    integer ::  dof,  Qdof, i

    dof = elem%dof
    Qdof = elem%Qdof

    elem%w(0,1:ndim*dof) = 0.
    do i=1,state%time%deg
       elem%w(0,:) = elem%w(0,:) - state%time%alpha(i) * elem%w(i,:)
       !if(elem%i == 10)
       !write(*,'(a3,2i5,12es12.4)') &
       !     ',@@',i,state%time%deg, state%time%alpha(i), elem%w(i,1) ,elem%w(0,1)
    enddo

    elem%vec(rhsM,1:ndim*dof) = 0.

    ! setting of w in integration nodes
    allocate(wi(1:Qdof,1:ndim))
    call Eval_w_Elem(elem, wi(1:Qdof,1:ndim) )

    ! evaluation of the vector
    call EvalVectorB(elem, wi(1:Qdof,1:ndim), dof, elem%vec(rhsM,1:ndim*dof) )

    deallocate(wi)

  end subroutine SetElementMassVector

  !> compute \f$ \int_K {\bf w}  \phi_i dx \f$
  subroutine SetElementMassVectorS(elem, vecM)
    type(element):: elem
    real, dimension(1:ndim*elem%dof), intent(inout) :: vecM
    real, dimension(:,:), allocatable :: wi
    integer ::  dof,  Qdof, i

    dof = elem%dof
    Qdof = elem%Qdof

    vecM(:) = 0.

    ! setting of w in integration nodes
    allocate(wi(1:Qdof,1:ndim))
    call Eval_w_Elem(elem, wi(1:Qdof,1:ndim) )

    ! evaluation of the vector
    call EvalVectorB(elem, wi(1:Qdof,1:ndim), dof, vecM(1:ndim*dof) )

    deallocate(wi)

  end subroutine SetElementMassVectorS

  !> compute \f$ \int_K f  \phi_i dx \f$
  subroutine ElementRHS(elem)
    type(element):: elem
    real, dimension(:,:), allocatable :: f !!!xi, x
    integer ::  dof,  dofA, Qdof, Qnum,  l

    dof = elem%dof
    dofA = dof
    if(elem%deg_plus) dofA = elem%dof_plus

    ! solution in integ nodes
    !allocate(wi(1:Qdof,1:ndim) )
    !call Eval_w_Elem(elem, wi(1:Qdof,1:ndim) )

    Qdof = elem%Qdof
    !Qnum = elem%Qnum

    allocate( f(1:Qdof, 1:ndim) )
    !allocate( xi(1:Qdof, 1:nbDim))
    !allocate( x(1:Qdof, 1:nbDim))

    ! setting of the function \f$ f \f$ in integ nodes x(:, 1:nbDim)
    !xi(1:Qdof,1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)

    !call ComputeF(elem, Qdof, xi(1:Qdof,1:nbDim), x(1:Qdof, 1:nbDim) )

    do l=1,Qdof
       if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
          !call RHS_Scalar( x(l,1:nbDim), f(l,1:ndim), state%time%ctime )
          call RHS_Scalar( elem%xi(0, l, 1:nbDim), f(l,1:ndim), state%time%ctime )
       elseif(ndim == 4 .and. state%type_IC .eq. 9 ) then
          !call RHS_SteadyChannel( x(l,1:nbDim), f(l,1:ndim), state%time%ctime )
          call RHS_SteadyChannel(elem%xi(0, l, 1:nbDim), f(l,1:ndim), state%time%ctime )
       endif
    enddo

    ! wet steam source terms
    !if(state%modelName == 'wet_steam') then
    !   call RHS_WS(elem,elem%xi(0,1:Qdof,1:nbdim),f(1:Qdof,1:ndim))
    !end if

    ! evaluation of the vector
    call EvalVectorB(elem, f(1:Qdof,1:ndim), dofA, elem%vec(rhs,1:ndim*dofA) )
    !write(*,'(a6,i5,26es12.4)') 'f(x_i)',Qdof, f(:,1)
    !write(*,'(a6,i5,6es12.4)') 'in_RHS',elem%i, elem%vec(rhs,1:ndim*dofA)

    deallocate(f) !!, x, xi)

  end subroutine ElementRHS


    !> compute \f$ \int_K f  \phi_i dx \f$
    !> output = elem%vec(rhs,1:ndim*dofA)
  subroutine ElementSubdomainRHS(elem)
    type(element):: elem
    real, dimension(:,:), allocatable :: f !!!xi, x
    integer ::  dof,  dofA, Qdof, Qnum,  l , k
    real, allocatable, dimension(:,:,:) :: g

    if (ndim > 1) &
      stop ' ElementSubdomainRHS there may be problem with computing f due to ndim>1!'

    dof = elem%dof
    dofA = dof
    if(elem%deg_plus) dofA = elem%dof_plus


    ! is the element in the support of the PRIMAL RHS ?
    if ( elem%iSubMesh == -1 ) then

!       print*, elem%i , 'element iSubmesh == -1'


       Qdof = elem%Qdof

       allocate( f(1:Qdof, 1:ndim) )

       ! setting of the function \f$ f \f$ in integ nodes x(:, 1:nbDim)
       !xi(1:Qdof,1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)
       !call ComputeF(elem, Qdof, xi(1:Qdof,1:nbDim), x(1:Qdof, 1:nbDim) )

       do l=1,Qdof
          if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
             call RHS_Scalar( elem%xi(0, l, 1:nbDim), f(l,1:ndim), state%time%ctime )
          elseif(ndim == 4 .and. state%type_IC .eq. 9 ) then
             call RHS_SteadyChannel(elem%xi(0, l, 1:nbDim), f(l,1:ndim), state%time%ctime )
          endif
       enddo

       ! evaluation of the vector
       k = state%model%rhsTestFunDerivative
       ! standart computation f*phi dx
       if ( k == 0 ) then
            call EvalVectorB(elem, f(1:Qdof,1:ndim), dofA, elem%vec(rhs,1:ndim*dofA) )
       ! derivatives of the test functions in the rhs
       else
            allocate( g(1:Qdof, 1:nbDim, 1:ndim), source = 0.0 )
            !dPhi/dx or dPhi/dy
            if (k==1 .or. k == 2) then
               g( 1:Qdof , k , 1:ndim ) = f(1:Qdof, 1:ndim)
               !  !> integrate \f$ Vector(R) =
               ! \int_{K_{elem}} (f_1 \partial_1 \phi_{R} + f_2 \partial_2 \phi_{R} )\ dx \f$
               ! \f$ \forall R \f$
               call EvalVectorD( elem, g(1:Qdof, 1:nbDim, 1:ndim), dofA, elem%vec(rhs,1:ndim*dofA) )
            !dPhi/dy
!            elseif (k==2) then
!              ! call EvalVectorD(elem, func, dofA, Vector)
!            !unknown
            else
               stop 'unknown choice of k in ElementSubdomainRHS'
            endif
            deallocate(g)
       endif

       !write(*,'(a6,i5,26es12.4)') 'f(x_i)',Qdof, f(:,1)
       !write(*,'(a6,i5,6es12.4)') 'in_RHS',elem%i, elem%vec(rhs,1:ndim*dofA)

       deallocate(f) !!, x, xi)

    else
      ! DO NOTHING - do not set 0.0 (for implicitly=True there are more inputs to rhs)
      ! DO NOT !!!! elem%vec(rhs,1:ndim*dofA) = 0.0
    endif

  end subroutine ElementSubdomainRHS

  !> evaluation of flow coefficients: \f$ c_D,\ c_L, c_M, c_{D,p}, c_{L,p} \f$
  subroutine DirectComputeDrag_Lift(coeffs)
    real, dimension (1:5), intent(inout) :: coeffs
    real :: CD, CD_p, CL, CL_p, CM ! values of drag, lift and moment c.
    class(element), pointer:: elem
    real, dimension(:,:), allocatable :: nc          ! normals in integ nodes
    real, dimension(:,:), allocatable :: x_div  ! local arrays of modif. int. nds

    real, dimension(:,:,:), allocatable :: ST ! stress tensor) in integration nodes
    real, dimension(:), allocatable :: press        ! - pressure  in integration nodes
    real, dimension(:, :), allocatable :: wi        ! w in integ nodes
    real, dimension(:, :, :), allocatable :: Dwi    ! Dw in integ nodes
    real, dimension(:), allocatable :: temp_CD, temp_CL, temp_CM,temp_CD_p,temp_CL_p

    integer :: ib, ie, j, k, l, kst, dof, Qnum, Qdof,  OK, i, jbc
    real :: t1, t2, new_CD, new_CL, new_CD_p, new_CL_p
    real, dimension(2) :: x_ref

    call cpu_time(t1)


    allocate(ST(1:state%space%max_Qdof, 1:nbDim, 1:nbDim) )

    allocate(wi(1:state%space%max_Qdof, 1:ndim) )
    allocate(Dwi(1:state%space%max_Qdof, 1:ndim, 1:nbDim ) )
    allocate(press(1:state%space%max_Qdof) )

    allocate(temp_CD(1:state%space%max_Qdof) )
    allocate(temp_CL(1:state%space%max_Qdof) )
    allocate(temp_CM(1:state%space%max_Qdof) )
    allocate(temp_CD_p(1:state%space%max_Qdof) )
    allocate(temp_CL_p(1:state%space%max_Qdof) )

    OK = 0  !test OK

    !inicialization
    x_ref(1)=0.25
    x_ref(2)=0.0
    CD = 0.0
    CL = 0.0
    CM = 0.0
    CD_p = 0.0
    CL_p = 0.0

    do ib =1,grid%nbelm
       if(grid%b_edge(ib)%BC == 0 ) then  ! impermeable walls = profiles

          OK = 1
          elem => grid%elem(grid%b_edge(ib)%itc)
          ie = grid%b_edge(ib)%jtc

          Qnum = elem%face(fGnum,ie)
          Qdof = elem%face(fGdof,ie)

          allocate(x_div(1:Qdof, 1:nbDim))
          x_div(1:Qdof, 1) = grid%b_edge(ib)%x_div(1:Qdof,1) - x_ref(1)
          x_div(1:Qdof, 2) = grid%b_edge(ib)%x_div(1:Qdof,2) - x_ref(2)

          ! setting of vector solution w in integ nodes
          call Eval_w_Edge(elem, ie, wi(1:Qdof, 1:ndim),.false. )

          ! setting of derivatives of vector solution w in integ nodes
          call Eval_Dw_Edge(elem, ie, Dwi(1:Qdof, 1:ndim,1:nbDim),.false. )

          !do l=1,Qdof
          !   Dwi(l, 1:ndim,1:nbDim) = Dwi(l, 1:ndim,1:nbDim) / &
          !        state%space%G_rule(Qnum)%weights(l)
          !enddo

          ! setting of stress tensor
          call Set_Stress_Tensor(ndim, Qdof, wi(1:Qdof,1:ndim), &
               Dwi(1:Qdof,1:ndim,1), Dwi(1:Qdof,1:ndim,2), &
               ST(1:Qdof,1:nbDim,1:nbDim), press(1:Qdof) )

          allocate(nc(1:Qdof, 1:nbDim) )

          ! setting of outer normals in integration nodes
          if(elem%ibcur > 0) then
             nc(1:Qdof,1:nbDim) = elem%nc(1:Qdof,1:nbDim)
          else
             nc(1:Qdof,1) = elem%n(ie,1)
             nc(1:Qdof,2) = elem%n(ie,2)
          endif


          ! temps
          temp_CD(1:Qdof) = ST(1:Qdof,1,1)*nc(1:Qdof,1) + ST(1:Qdof,1,2)*nc(1:Qdof,2)
          temp_CL(1:Qdof) = ST(1:Qdof,2,1)*nc(1:Qdof,1) + ST(1:Qdof,2,2)*nc(1:Qdof,2)
          temp_CM(1:Qdof) = x_div(1:Qdof,1)*temp_CL(1:Qdof) &
               - x_div(1:Qdof,2)*temp_CD(1:Qdof)

          temp_CD_p(1:Qdof) = -press(1:Qdof)*nc(1:Qdof,1)
          temp_CL_p(1:Qdof) = -press(1:Qdof)*nc(1:Qdof,2)

          !write(*,'(a8,2i5,4es12.4,a2,20es12.4)') 'cDLM_ed3a',elem%i, Qdof, &
          !     maxval(press(1:Qdof)), minval(press(1:Qdof)), &
          !     maxval(press(1:Qdof))/ minval(press(1:Qdof)), &
          !     (maxval(press(1:Qdof))- minval(press(1:Qdof)))/elem%diam, '|', elem%xc(:)


          ! if(state%time%iter_loc > 25) then
          !    write(100*(state%space%adapt%adapt_level+1) + state%time%iter_loc,'(20es12.4)') &
          !         elem%xc(:), &
          !         maxval(press(1:Qdof)), minval(press(1:Qdof)), &
          !         maxval(press(1:Qdof))/ minval(press(1:Qdof)), &
          !         (maxval(press(1:Qdof))- minval(press(1:Qdof)))/elem%diam

          !    do j=1,Qdof
          !       write(1000+100*(state%space%adapt%adapt_level+1) + state%time%iter_loc,'(20es12.4)') &
          !            x_div(j, 1:2)+x_ref(1:2), &
          !            temp_CD_p(j), temp_CL_p(j), temp_CD(j), temp_CL(j), temp_CM(j), &
          !            press(j),nc(j,1:2), nc(j,1:2)/sqrt(dot_product(nc(j,1:2), nc(j,1:2))), &
          !            dot_product(temp_CD(1:Qdof), state%space%G_rule(Qnum)%weights(1:Qdof)), &
          !            dot_product(temp_CL(1:Qdof), state%space%G_rule(Qnum)%weights(1:Qdof))

          !       ! write(2000+100*(state%space%adapt%adapt_level+1) + state%time%iter_loc,'(20es12.4)') &
          !       !      x_div(j, 1:2)+x_ref(1:2)
          !       ! write(2000+100*(state%space%adapt%adapt_level+1) + state%time%iter_loc,'(20es12.4)') &
          !       !      x_div(j, 1:2)+x_ref(1:2) + nc(j,1:2)
          !       ! write(2000+100*(state%space%adapt%adapt_level+1) + state%time%iter_loc,'(x)')


          !       ! write(3000+100*(state%space%adapt%adapt_level+1) + state%time%iter_loc,'(20es12.4)') &
          !       !      x_div(j, 1:2)+x_ref(1:2)
          !       ! write(3000+100*(state%space%adapt%adapt_level+1) + state%time%iter_loc,'(20es12.4)') &
          !       !      x_div(j, 1:2)+x_ref(1:2) +elem%n(ie,1:2)
          !       ! write(3000+100*(state%space%adapt%adapt_level+1) + state%time%iter_loc,'(x)')
          !    enddo
          ! endif



          ! multiply by weights
          CD = CD - dot_product(temp_CD(1:Qdof), state%space%G_rule(Qnum)%weights(1:Qdof))
          CL = CL - dot_product(temp_CL(1:Qdof), state%space%G_rule(Qnum)%weights(1:Qdof))
          CM = CM - dot_product(temp_CM(1:Qdof), state%space%G_rule(Qnum)%weights(1:Qdof))

          CD_p = CD_p &
               - dot_product(temp_CD_p(1:Qdof), state%space%G_rule(Qnum)%weights(1:Qdof))
          CL_p = CL_p &
               - dot_product(temp_CL_p(1:Qdof), state%space%G_rule(Qnum)%weights(1:Qdof))


          deallocate(nc, x_div)

          !write(*,'(a2,i5,4es14.6)') '!!',ib, &
          !     dot_product(temp_CD(1:Qdof), state%space%G_rule(Qnum)%weights(1:Qdof)),CD,&
          !     dot_product(temp_CL(1:Qdof), state%space%G_rule(Qnum)%weights(1:Qdof)),CL

       endif  ! end of impermeable walls
    enddo ! end for all boundary edges


    !close(ifile)

    !modification transformation with respect to the angle of attack
    new_CD =  cos(state%alpha_infty)*CD + sin(state%alpha_infty)*CL
    new_CL = -sin(state%alpha_infty)*CD + cos(state%alpha_infty)*CL

    new_CD_p =  cos(state%alpha_infty)*CD_p + sin(state%alpha_infty)*CL_p
    new_CL_p = -sin(state%alpha_infty)*CD_p + cos(state%alpha_infty)*CL_p

    !write(*,'(a3,4es14.6)')'???',CD,CL,new_CD,new_CL

    CD= new_CD *2.0 / (state%rho_infty * state%v_infty * state%v_infty)
    CL= new_CL *2.0 / (state%rho_infty * state%v_infty * state%v_infty)
    CM =    CM *2.0 / (state%rho_infty * state%v_infty * state%v_infty)

    !write(*,'(a3,4es14.6)')'???',CD,CL,new_CD,new_CL

    CD_P= new_CD_P * 2.0/ (state%rho_infty * state%v_infty * state%v_infty)
    CL_P= new_CL_p * 2.0/ (state%rho_infty * state%v_infty * state%v_infty)

    !print*,'!!!!!!!!!',state%rho_infty, state%v_infty, state%alpha_infty
    !print*,'final: ',CD, CL, CM

    if (OK == 0) print*, 'No boundary segments on fixed walls  ',  &
         'in DirectCompute CD_CL'

    coeffs(1) = CD
    coeffs(2) = CL
    coeffs(3) = CM
    coeffs(4) = CD_p
    coeffs(5) = CL_p

    deallocate(ST, wi, Dwi, press)
    deallocate(temp_CD, temp_CL, temp_CM, temp_CD_p, temp_CL_p )

    call cpu_time(t2)

    !!if(state%itime > 0) &
    !print*,'# DirectComputeDrag_Lift finished after ',t2-t1, ' s'
    !print*,'CCCC:',coeffs(1:5)
    !stop

  end subroutine DirectComputeDrag_Lift


  !> compute the stress tensor and the pressure in integ nodes:
  !>
  !> \f$ \tau_{ij} = -p \delta_{ij} +
  !> \frac{1}{Re}\left(\frac{\partial v_i}{\partial x_j}
  !> + \frac{\partial v_j}{\partial x_i}
  !> -\frac23\mbox{div } \vec{v}\delta_{ij}\right)\f$
  subroutine Set_Stress_Tensor(ndimL, Qdof, w, w_x1, w_x2, ST, press)
    integer, intent(in) :: Qdof, ndimL
    real, dimension(1:Qdof, 1:ndim), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndim), intent(in):: w_x1 !state  dw/dx in #Qdof nodes
    real, dimension(1:Qdof, 1:ndim), intent(in):: w_x2 !state  dw/dy in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:nbDim), intent(inout) :: ST  ! stress tensor in integ nodes

    real, dimension(1:Qdof), intent(inout) :: press
    real, dimension(:), allocatable :: dv1dx, dv1dy, dv2dx, dv2dy

    real :: kappa1, Re

    if(nbDim > 2)  print*,'inv_fluxes.f90 Set_Stress_Tensor: attention'

    allocate( dv1dx(1:Qdof), dv1dy(1:Qdof), dv2dx(1:Qdof), dv2dy(1:Qdof) )

    kappa1 = state%model%kappa1
    Re = state%model%Re


    ! for each node
    ! setting partial derivatives of velocity field
    dv1dx(1:Qdof) = w_x1(1:Qdof,2)/w(1:Qdof,1) &
         - w_x1(1:Qdof,1)*w(1:Qdof,2)/(w(1:Qdof,1)*w(1:Qdof,1))

    dv1dy(1:Qdof) = w_x2(1:Qdof,2)/w(1:Qdof,1) &
         - w_x2(1:Qdof,1)*w(1:Qdof,2)/(w(1:Qdof,1)*w(1:Qdof,1))

    dv2dx(1:Qdof) = w_x1(1:Qdof,3)/w(1:Qdof,1) &
         - w_x1(1:Qdof,1)*w(1:Qdof,3)/(w(1:Qdof,1)*w(1:Qdof,1))

    dv2dy(1:Qdof) = w_x2(1:Qdof,3)/w(1:Qdof,1) &
         - w_x2(1:Qdof,1)*w(1:Qdof,3)/(w(1:Qdof,1)*w(1:Qdof,1))

    ! setting pressure from energy
    press(1:Qdof) = kappa1*(w(1:Qdof,4) &
         - 0.5*(w(1:Qdof,2)*w(1:Qdof,2) + w(1:Qdof,3)*w(1:Qdof,3))/w(1:Qdof,1))

    ! setting stress tensor
    if(Re > 0.) then
       ST(1:Qdof,1,1) = -press(1:Qdof) - 2./(3.*Re)*(dv1dx(1:Qdof) + dv2dy(1:Qdof))&
            + 2.*dv1dx(1:Qdof)/Re
       ST(1:Qdof,1,2) = (dv1dy(1:Qdof) + dv2dx(1:Qdof))/Re

       ST(1:Qdof,2,1) = ST(1:Qdof,1,2)

       ST(1:Qdof,2,2) = -press(1:Qdof) - 2./(3.*Re)*(dv1dx(1:Qdof) + dv2dy(1:Qdof))&
            + 2.*dv2dy(1:Qdof)/Re
    else

       ST(1:Qdof,1,1) = -press(1:Qdof)
       ST(1:Qdof,1,2) = 0.0
       ST(1:Qdof,2,1) = 0.0
       ST(1:Qdof,2,2) = -press(1:Qdof)
    endif

    deallocate( dv1dx, dv1dy, dv2dx, dv2dy )

  end subroutine Set_Stress_Tensor


  !> computation of inviscid volume integrals on \f$elem=K\f$, i.e.,
  !> \f$-\int_K \sum_{s=1}^2 (f_s(w) \partial_s \phi_R dx \f$
  subroutine ExplicitInviscidVolumes(elem,  Set_f_s)
    type(element):: elem
    ! Compute matrices A_s(w)
    interface
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x, ie )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
         integer, intent(in) :: ie
      end subroutine Set_f_s
    end interface
    real, dimension(:,:), allocatable :: wi ! w recomputed  in integ nodes
    real, dimension(:,:,:), allocatable ::f_s ! vectors f_s
    integer ::  Qdof, dof

    dof = elem%dof
    Qdof = elem%Qdof

    ! setting of the state vector in integration nodes
    allocate(wi(1:Qdof,1:ndim) )
    call Eval_w_Elem(elem, wi(1:Qdof,1:ndim) )

    ! setting of matrixes A_s in integration nodes
    allocate(f_s(1:Qdof,1:nbDim,1:ndim))
    call Set_f_s(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), f_s(1:Qdof, 1:nbDim, 1:ndim), &
         elem%xi(0,1:Qdof, 1:nbDim), elem%i)

    ! evaluation of matrix terms
    call EvalVectorD(elem,  f_s(1:Qdof,1:nbDim,1:ndim), dof, elem%vec(rhs,1:ndim*dof) )

    deallocate(wi, f_s)
  end subroutine ExplicitInviscidVolumes



  !> explicit evaluation of inviscid edge integrals
  subroutine ExplicitInviscidEdge(elem, elem1, ie, Set_NumFlux)
    type(element), intent(inout):: elem, elem1  ! elem = element, elem1 = neigh elem
    integer, intent (in) :: ie                  ! inner index of edge, 1,2,3, (4)
    interface
      subroutine Set_NumFlux(ndimL, nbDim, Qdof, wi, wj, nc, xi, H, area_1, ie )
          integer, intent(in) :: Qdof, ndimL, nbDim, ie
          real, dimension(1:Qdof, 1:ndimL), intent(in):: wi, wj ! state  w in integ nodes
          real, dimension(1:Qdof, 1:nbDim), intent(in):: nc        ! outer normal in integ nodes
          real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
          real, dimension(1:Qdof,1:ndimL), intent(inout) :: H   ! numer flux H in  -- " --
          real, intent(in) :: area_1
      end subroutine
    end interface
    real, dimension(:,:), allocatable :: wi,wi1  ! w in integ nodes
    real, dimension(:,:), allocatable :: nc      ! outer normal in integ nodes
    real, dimension(:,:), allocatable :: H       ! numerical flux in integ nodes
!    real, dimension(1:nbDim) :: xi
    real, dimension(:,:), allocatable :: xi !FR changed
    real :: param, lambda_LF
    integer ::   dof, dof1, ie1, Qdof, k

    write(debug, *) 'Problem in arguments when calling Set_NumFlux!!! CONTROL: xi(1:nbDim) changed to xi(1:Qdof, 1:nbDim)'

    dof = elem%dof
    ie1 = elem%face(nei_i,ie)
    Qdof = elem%face(fGdof,ie)

    allocate(wi(1:Qdof,1:ndim) )
    call Eval_w_Edge(elem, ie, wi(1:Qdof,1:ndim), .false.)

    allocate( xi(1:Qdof, 1:nbDim) )


    allocate(nc(1:Qdof, 1:nbDim) )     ! setting of outer normals in integ nodes
    allocate(wi1(1:Qdof,1:ndim) )  ! outside solution in integ nodes
    if( elem%face(neigh,ie) > 0) then
       ! inner edge
       call Eval_w_Edge(elem1, ie1, wi1(1:Qdof,1:ndim), .true.)
       nc(1:Qdof,1) = elem%n(ie,1)
       nc(1:Qdof,2) = elem%n(ie,2)

    else
       ! extrapolation
       wi1(1:Qdof,1:ndim)= wi(1:Qdof,1:ndim)
       ! Dirichlet BC
       do k=1,Qdof
          xi(k, 1:nbDim) = grid%b_edge(-elem%face(neigh,ie))%x_div(k, 1:nbDim)
          if(state%time%tdg) then
             call Exact_Scalar(xi(k,1:nbDim), wi1(k,1:ndim), state%time%ctime)
          else
             call Exact_Scalar(xi(k,1:nbDim), wi1(k,1:ndim), state%time%ttime) ! EXPLICIT !!!
             !call Exact_Scalar(xi(1:nbDim), wi1(k,1:ndim), state%time%ctime) ! IMPLICIT !!!
          endif
       enddo

       !write(*,'(a4,8es12.4)') 'wi: ',wi(1:Qdof, 1), elem%w(0,:)
       !write(*,'(a4,4es12.4)') 'wD: ',wi1(1:Qdof, 1)
       !write(*,'(a4,4es12.4)') 'dif:',wi(1:Qdof, 1) - wi1(1:Qdof, 1)
       !print*,'@@@',state%time%ttime, state%time%ctime

       ! outer normal
       if(elem%ibcur > 0) then
          nc(1:Qdof,1:nbDim) = elem%nc(1:Qdof,1:nbDim)
       else
          ! polygonal edge
          nc(1:Qdof,1) = elem%n(ie,1)
          nc(1:Qdof,2) = elem%n(ie,2)
       endif

    endif


    allocate(H(1:Qdof,1:ndim) )
    ! numerical flux  - ndimL, nbDim, Qdof, wi, wj, nc, xi, H, area_1, ie
    call Set_NumFlux( ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), wi1(1:Qdof,1:ndim), nc(1:Qdof,1:nbDim), &
         xi(1:Qdof, 1:nbDim), H(1:Qdof, 1:ndim), 1./elem%area,elem%i )

    !if( elem%face(neigh,ie) <= 0) then
    !   do k=1,Qdof
    !      write(*,'(a4,4es12.4)') ' H:', &
    !           grid%b_edge(-elem%face(neigh,ie))%x_div(k, 1:nbDim), &
    !           H(k, 1), grid%b_edge(-elem%face(neigh,ie))%x_div(k, 1)**2/2 &
    !           *sum(elem%n(ie,1:nbDim))
    !   enddo
    !endif

    !! stabilization
    !param = (elem%rezid + elem1%rezid)/2.* elem%dn(ie)**state%ST_Ep * state%ST_Ec &
    !     * elem%dn(ie)  ! Ppm is already multiplied by |n| = dn

    !if(elem%i < elem1%i) &
    !     print*,'I', elem%i, elem1%i, elem%dn(ie),& ! elem%rezid, elem1%rezid,
    !     param/elem%dn(ie)

    !do k=1,ndim                 ! k = index of component of w
    !   !!write(state%time%iter+51,*) elem%xc(1:nbDim), Ppm(1:Qdof, 1, k, k), param ,k,elem%i
    !   H(1:Qdof, k) = H(1:Qdof, k) + param !!+ 1./lambda_LF
    !enddo

    ! evaluation of vector terms
    call ExplEdgeB(elem, ie, -H(1:Qdof, 1:ndim),  dof, elem%vec(rhs,1:ndim *dof ) )

    deallocate(H )
    deallocate (wi, wi1)
    deallocate (nc, xi)
  end subroutine ExplicitInviscidEdge


  !> computation of linear inviscid volume integrals on \f$elem=K\f$, i.e.,
  !> \f$ \int_K \sum_{s=1}^2 (A_s(w) \phi_C) \partial_s \phi_R dx \f$
  subroutine LinearInviscidVolumes(elem)
    type(element):: elem
    real, dimension(:,:,:,:), allocatable :: conv ! linear convection in integ nodes
    real, dimension(:,:,:), allocatable :: reac   ! linear reaction in integ nodes
    real, dimension(:,:), allocatable :: xi, x
    integer ::  Qdof, dof,  Qnum
    integer :: i,j

    dof = elem%dof
    Qnum = elem%Qnum
    Qdof = elem%Qdof

    ! evaluation of integration nodes
    allocate( xi(1:Qdof, 1:nbDim))
    allocate( x(1:Qdof, 1:nbDim))

    xi(1:Qdof,1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)
    call ComputeF(elem, Qdof, xi(1:Qdof,1:nbDim), x(1:Qdof, 1:nbDim) )


    ! CONVECTION
    ! setting of the convection in integration nodes
    allocate(conv(1:Qdof, 1:nbDim, 1:ndim, 1:ndim) )

    call LinearConvection(Qdof, ndim, x(1:Qdof,1:nbDim), state%time%ctime, &
         conv(1:Qdof,1:nbDim, 1:ndim, 1:ndim))

    call EvalBlockDB(elem,  conv(1:Qdof,1:nbDim,1:ndim,1:ndim), &
               elem%block(0)%Mb(1:ndim*dof, 1:ndim*dof ) )


    ! REACTION
    ! setting of the reaction in integration nodes
    allocate(reac(1:Qdof, 1:ndim, 1:ndim) )

    call LinearReaction(Qdof, ndim, x(1:Qdof,1:nbDim), state%time%ctime, &
         reac(1:Qdof,1:ndim,1:ndim))


    call EvalBlockBB(elem,  reac(1:Qdof, 1:ndim,1:ndim), &
               elem%block(0)%Mb(1:ndim*dof, 1:ndim*dof ) )

    !call WriteMblock(elem%block(0) )

    deallocate(conv, reac)
  end subroutine LinearInviscidVolumes

  !> Setting of names for tri* sol*
  !> if command = .true. then setting of the commmand for plotting
  subroutine SetCFFileNames(CF_name, CF_name1)
    character(len=50),intent(inout) :: CF_name, CF_name1
    character(len=5) :: ch5
    integer :: num_size, text_size, long_text_size, file_size, inum
    integer :: is

    long_text_size = 10
    text_size = 4
    num_size = 5
    file_size = text_size + num_size

    !print*,'###',state%space%adapt%max_adapt_level,  state%time%OutTime

    if(state%space%adapt%max_adapt_level == 0 .or. state%time%OutTime > 0.) then
       inum = state%isol
       CF_name  = 'CFp-00000    '
       CF_name1 = 'CFc-00000    '
    else
       inum = state%space%adapt%adapt_level
       if(state%space%adapt%adapt_level < 0) inum = 10**num_size + state%space%adapt%adapt_level
       CF_name  = 'CFpA00000    '
       CF_name1 = 'CFcA00000    '
    endif

    if(inum > 0) then
       is = int(log(1.*inum)/log(10.))
    else
       is = 0
    endif

    !print*,'!!!',inum,is, num_size+text_size-is, num_size+text_size, num_size-is, num_size

    write( ch5, '(i5)' ) inum  ! change the format if num_size /= 5 !!!
    CF_name(num_size+text_size-is:num_size+text_size) = ch5(num_size-is: num_size)
    CF_name1(num_size+text_size-is:num_size+text_size) = ch5(num_size-is: num_size)
  end subroutine SetCFFileNames

!> compute and visualise the skin friction coefficient  along bound impermeable walls
  subroutine ComputeSkinFriction()
    class(element), pointer:: elem
    real, dimension(:,:), allocatable :: nc       ! normals in integ nodes
    real, dimension(:,:), allocatable :: tc       ! tangents in integ nodes
    real, dimension(:,:), allocatable :: wi       ! w in integ nodes
    real, dimension(:,:,:), allocatable :: Dwi    ! Dw in integ nodes
    real, dimension(:,:,:), allocatable :: ST     ! stress tensor in integ nodes
    real, dimension(:), allocatable :: press      ! - pressure  in integ nodes
    character(len=50) :: CF_name, CF_name1

    integer :: ie, ib, l, Qdof,Qnum
    logical :: wall_present
    integer :: ifile = 15, ifile1 = 16
    real ::  skin, skin_aver


    allocate(wi(1:state%space%max_Qdof, 1:ndim) )
    allocate(Dwi(1:state%space%max_Qdof, 1:ndim, 1:nbDim ) )
    allocate(ST(1:state%space%max_Qdof, 1:nbDim, 1:nbDim) )
    allocate(press(1:state%space%max_Qdof) )

    call SetCFFileNames(CF_name, CF_name1)
    open(ifile, file=CF_name, status ='UNKNOWN')
    open(ifile1, file=CF_name1, status ='UNKNOWN')

    ! verify if impermeable walls are present
    wall_present = .false.

!    print*,'# rho_infty, v_infty',ro_infty,v_infty
    do ib =1,grid%nbelm
       if(grid%b_edge(ib)%BC == 0 ) then
          wall_present = .true.

          elem => grid%elem(grid%b_edge(ib)%itc)
          ie = grid%b_edge(ib)%jtc

          Qdof = elem%face(fGdof,ie)
          Qnum = elem%face(fGnum,ie)

          ! setting of vector solution w in integ nodes
          call Eval_w_Edge(elem, ie, wi(1:Qdof, 1:ndim),.false. )

          ! setting of derivatives of vector solution w in integ nodes
          call Eval_Dw_Edge(elem, ie, Dwi(1:Qdof, 1:ndim,1:nbDim),.false. )

          ! setting of stress tensor
          call Set_Stress_Tensor(ndim, Qdof, wi(1:Qdof,1:ndim), &
               Dwi(1:Qdof,1:ndim,1), Dwi(1:Qdof,1:ndim,2), &
               ST(1:Qdof,1:nbDim,1:nbDim), press(1:Qdof) )

          allocate(nc(1:Qdof, 1:nbDim), tc(1:Qdof, 1:nbDim) )

          ! setting of outer normals in integration nodes
          if(elem%ibcur > 0) then
             nc(1:Qdof,1) = elem%nc(1:Qdof,1) / elem%dnc(1:Qdof)
             nc(1:Qdof,2) = elem%nc(1:Qdof,2) / elem%dnc(1:Qdof)
          else
             nc(1:Qdof,1) = elem%n(ie,1) / elem%dn(ie)
             nc(1:Qdof,2) = elem%n(ie,2) / elem%dn(ie)
          endif

          tc(1:Qdof,1) =  nc(1:Qdof,2)
          tc(1:Qdof,2) =  -nc(1:Qdof,1)

          !nn(1:nbDim) = grid%b_edge(ie)%n(1:nbDim)
          !len = (dot_product(nn(1:nbDim), nn(1:nbDim)))**0.5
          !nn(1:nbDim) = nn(1:nbDim)/len
          !tt(1) = nn(2)
          !tt(2) = -nn(1)

          ! only viscous part of stress tensor
          ST(1:Qdof,1,1) = ST(1:Qdof,1,1) - press(1:Qdof)
          ST(1:Qdof,2,2) = ST(1:Qdof,2,2) - press(1:Qdof)

          skin_aver = 0.
          do l=1, Qdof

             skin = tc(l,1)*(ST(l,1,1) * nc(l,1) + ST(l,1,2) * nc(l,2)) &
                  + tc(l,2)*(ST(l,2,1) * nc(l,1) + ST(l,2,2) * nc(l,2))
             skin = 2*skin/state%rho_infty/state%v_infty/state%v_infty

             skin_aver = skin_aver + skin * state%space%G_rule(Qnum)%weights(l)

             write(ifile,*) grid%b_edge(ib)%x_div(l,1:nbDim), skin
          enddo

          write(ifile,*) '   '

          write(ifile1,*) grid%b_edge(ib)%x_div((Qdof+1)/2,1:nbDim), skin_aver

          deallocate(nc, tc)
       end if
    end do

    if(wall_present) then
    !   print*, '# CF coefficient computed in files: CF, CFaver'
    else
       print*, 'No boundary segment on fixed wall in grid for ComputeCF'
    endif

    deallocate(ST, wi, Dwi, press)

    close(ifile)
    close(ifile1)

  end subroutine ComputeSkinFriction


  !> computation of reaction volume integrals on \f$elem=K\f$, i.e.,
  !> \f$-\int_K (S(w) \phi_C \phi_R dx \f$
  subroutine ElementReactionVolumes(elem, Set_S,  Set_DS)
    type(element):: elem
    ! Compute inviscid fluxes f_s and matrices A_s(w)
    interface
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

    real, dimension(:,:), allocatable :: wi ! w recomputed  in integ nodes
    real, dimension(:,:,:), allocatable :: Dwi ! Dw recomputed  in integ nodes
    real, dimension(:,:,:), allocatable :: DS ! matrices A
    real, dimension(:,:), allocatable :: S !   vectors S
    real :: param
    integer ::  Qdof, dof, dofA
    integer :: i

    dof = elem%dof
    dofA = dof
    if(elem%deg_plus) dofA = elem%dof_plus

    Qdof = elem%Qdof

    ! setting of the state vector in integration nodes
    allocate(wi(1:Qdof,1:ndim) )
    call Eval_w_Elem(elem, wi(1:Qdof,1:ndim) )

    allocate(Dwi(1:Qdof, 1:ndim, 1:nbDim) )
    call Eval_Dw_Elem(elem, Dwi(1:Qdof, 1:ndim, 1:nbDim) )


    if(state%nlSolver%implicitly) then
       ! setting of matrixes A_s in integration nodes
       allocate(DS(1:Qdof, 1:ndim, 1:ndim))
       call Set_DS(ndim, nbDim, Qdof, elem%xi(0,1:Qdof, 1:nbDim), wi(1:Qdof,1:ndim), &
            Dwi(1:Qdof, 1:ndim, 1:nbDim), DS(1:Qdof, 1:ndim, 1:ndim) )

       ! evaluation of matrix terms
       call EvalBlockBB(elem,  +DS(1:Qdof, 1:ndim,1:ndim), &
            elem%block(0)%Mb(1:ndim*dof, 1:ndim*dof ) )

       deallocate(DS)
    else ! explicitly, i.e,  state%nlSolver%implicitly = .false.

       ! setting of fluxes f_s in integration nodes
       allocate( S(1:Qdof, 1:ndim) )
       call Set_S(ndim, nbDim, Qdof, elem%xi(0,1:Qdof, 1:nbDim), wi(1:Qdof,1:ndim), &
            Dwi(1:Qdof, 1:ndim, 1:nbDim), S(1:Qdof, 1:ndim) )

       call EvalVectorB(elem, -S(1:Qdof,1:ndim), dofA, elem%vec(rhs, 1:ndim*dofA ) )


       !if(elem%i >= 5 .and. elem%i <= 5) then
       !!   !   !write(*,'(a4,i5,12es10.2)') 'VOL ',elem%i, wi(:, 1:ndim)
       !   !write(*,'(a4,2i5,120es14.6)') 'f_s ',elem%i, 1, S(:,1)
       !   write(*,'(a4,i5,i2,120es10.2)') 'f_s ',elem%i, 2, S(:,2),elem%xc(:), &
       !        wi(1,1),elem%xi(0,1, 1), wi(1,2), wi(1,2)/wi(1,1),wi(1,1)* elem%xi(0,1, 1) - wi(1,2)
       !   !write(*,'(a4,2i5,120es14.6)') 'f_s ',elem%i, 3, S(:,3)
       !   !print*
       !!   write(*,'(a4,i5,12es10.2)') '    ',dofA, elem%vec(rhs, 1:dofA)
       !endif

       deallocate(S)
    endif

    deallocate(wi, Dwi)
    !stop

  end subroutine ElementReactionVolumes



  !> computation of reaction volume integrals on \f$elem=K\f$, i.e.,
  !> \f$-\int_K (S(w) \phi_C \phi_R dx \f$
  subroutine ElementGradientStabilization(elem)  !!, Set_S,  Set_DS)
    type(element):: elem
    ! ! Compute inviscid fluxes f_s and matrices A_s(w)
    ! interface
    !    subroutine Set_S(ndimL, nbDim, Qdof, xi, w, Dw, S)
    !       integer, intent(in) :: ndimL, nbDim, Qdof
    !       real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
    !       real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    !       real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    !       real, dimension(1:Qdof, 1:ndimL), intent(inout) :: S
    !    end subroutine
    !    subroutine Set_DS(ndimL, nbDim, Qdof, xi, w, Dw, DS)
    !       integer, intent(in) :: ndimL, nbDim, Qdof
    !       real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
    !       real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    !       real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    !       real, dimension(1:Qdof, 1:ndimL, 1:ndimL), intent(inout) :: DS
    !    end subroutine
    ! end interface

    real, dimension(:,:), allocatable :: wi ! w recomputed  in integ nodes
    real, dimension(:,:,:), allocatable :: Dwi ! Dw recomputed  in integ nodes
    real, dimension(:,:,:, :), allocatable :: DS ! matrices A
    real, dimension(:,:), allocatable :: S !   vectors S
    real, dimension(:), allocatable :: penalty !   vectors penal
    real :: param, val
    integer ::  Qdof, dof, dofA
    integer :: i, k, l

    dof = elem%dof
    dofA = dof
    if(elem%deg_plus) dofA = elem%dof_plus

    Qdof = elem%Qdof

    ! setting of the state vector in integration nodes
    allocate(wi(1:Qdof,1:ndim) )
    call Eval_w_Elem(elem, wi(1:Qdof,1:ndim) )

    allocate(Dwi(1:Qdof, 1:ndim, 1:nbDim) )
    call Eval_Dw_Elem(elem, Dwi(1:Qdof, 1:ndim, 1:nbDim) )


    ! adding of the penalty leading to piecewise constant solution for vacuum state
    allocate(penalty(1:Qdof) )
    do i=1,Qdof
       call Set_grad_penalty(wi(i, 1), penalty(i))
    enddo
    val = maxval(penalty(:) )
    penalty(:) = val

    !if(maxval(penalty(:)) > 0. ) &
    !     write(100+state%time%iter , *) elem%xc(:), minval(wi(:, 1) ), maxval(penalty(:))

    if(state%nlSolver%implicitly) then
       ! setting of the linearization of penalty terms in integ nodes
       allocate(DS(1:Qdof, 1:nbDim, 1:ndim, 1:ndim))
       DS = 0.

       do i=1,Qdof
          do k= 1,ndim
             do l=1,nbDim
                DS(i, l, k, k) = penalty(i) * Dwi(i, k, l)
             enddo
          enddo
       enddo

       ! evaluation of matrix terms
       call EvalBlockDB(elem,  -DS(1:Qdof, 1:nbDim, 1:ndim,1:ndim), &
            elem%block(0)%Mb(1:ndim*dof, 1:ndim*dof ) )

       ! if(maxval(penalty(:)) > 0. )  then
       !    write(*,'(a4,2i5,120es14.6)') 'e(r)',elem%i, 0, penalty(:)
       !    write(*,'(a4,2i5,120es14.6)') 'DS ',elem%i, 1, DS(:,1, 1,1)
       !    write(*,'(a4,2i5,120es14.6)') 'DS ',elem%i, 2, DS(:, 1, 2, 2)
       !    write(*,'(a4,2i5,120es14.6)') 'DS ',elem%i, 3, DS(:,1, 3, 3)
       !    write(*,'(a4,2i5,120es14.6)') 'DS ',elem%i, 1, DS(:, 2, 1,1)
       !    write(*,'(a4,2i5,120es14.6)') 'DS ',elem%i, 2, DS(:, 2, 2, 2)
       !    write(*,'(a4,2i5,120es14.6)') 'DS ',elem%i, 3, DS(:,2, 3, 3)

       ! endif


       deallocate(DS)
    else ! explicitly, i.e,  state%nlSolver%implicitly = .false.

       ! setting of penalty in integration nodes
       allocate( S(1:Qdof, 1:ndim) )

       S = 0.

       do i=1,Qdof
          do k = 1,ndim
             S(i, k) = penalty(i) * dot_product(Dwi(i, k, 1:nbDim), Dwi(i, k, 1:nbDim) )
          enddo
       enddo

       call EvalVectorB(elem, S(1:Qdof,1:ndim), dofA, elem%vec(rhs, 1:ndim*dofA ) )


       ! if(maxval(penalty(:)) > 0. )  then
       !    write(*,'(a4,2i5,120es14.6)') 'e(r)',elem%i, 0, penalty(:)
       !    write(*,'(a4,2i5,120es14.6)') 'f_s ',elem%i, 1, S(:,1)
       !    write(*,'(a4,2i5,120es14.6)') 'f_s ',elem%i, 2, S(:,2)
       !    write(*,'(a4,2i5,120es14.6)') 'f_s ',elem%i, 3, S(:,3)
       ! endif

       deallocate(S)
    endif

    deallocate(wi, Dwi)
    !stop

  end subroutine ElementGradientStabilization



end module inviscid_fluxes
