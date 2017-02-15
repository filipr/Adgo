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
  use mesh_mod
  use eval_jumps
  use modelPorous
  use higher_order_local

  implicit none

  public:: ComputeTerms
  public:: ComputeElementsTerms
  public:: ComputeOneElementFluxes
  public:: ComputeOneElementTerms

  public:: ComputeSTDGM_Terms


  public:: FillVector
  public:: FillVectorST
  public:: PrepareCrankNicolson

  public:: SmoothResid
  public:: SetF_q_Cw_fast
  public:: EquivResid
  public:: EquivResid1

  public:: ReadConvFile

  public:: CheckResiduum

  public:: ClearMatrixBlocks
  public:: ClearMatrixBlocks_STDGM
  public:: ClearVectorBlocks
  public:: ClearVectorBlocks_STDGM

  public:: AvoidNonphysicalValues
  public:: PedestrianFlow_AvoidVacuum
  public:: PedestrianFlow_MinimalDensity
  public:: PedestrianEikonalEquation

  ! porous media flow
  public:: ComputeCapacityConductivity
contains



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
      grid%elem(i)%rhsST(1:ndim,1:grid%elem(i)%dof_plus,1:grid%elem(i)%Tdof_plus) = 0. ! clearing of rhs terms
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
      subroutine Set_R_s(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
       end subroutine Set_R_s
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x, ie )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
         integer, intent(in) :: ie
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
         (state%modelName == 'NSe' .and. state%type_IC .eq. 9) ) ) then

!         print*, 'subdomainRHS:' , state%model%subdomainRHS

         if ( state%model%subdomainRHS ) then
!            print*, 'ComputeOneElementFluxes' , state%nlSolver%implicitly
            call ElementSubdomainRHS( elem )
         else
            call ElementRHS(elem )
         endif

    endif
    !write(200+state%space%adapt%adapt_level,'(30es12.4)') elem%xc(:), abs(elem%vec(rhs,:) )
    !write(*,'(30es12.4)') elem%xc(:), elem%vec(rhs,1:dofA)

    ! ||phi||_{L^2(K)} = 1
    elem%eta(resS, 1) = (dot_product(elem%vec(rhs,1:dofA), elem%vec(rhs,1:dofA) ) &
         / elem%area)**0.5

    !elem%eta(resS, 1) = (dot_product(elem%vec(rhs,1:dofA), elem%vec(rhs,1:dofA) ) )**0.5

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
       end subroutine Set_A_s
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
       subroutine Set_S(ndimL, nbDim, Qdof, xi, w, Dw, S)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL), intent(inout) :: S
       end subroutine Set_S
       subroutine Set_DS(ndimL, nbDim, Qdof, xi, w, Dw, DS)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:ndimL), intent(inout) :: DS
       end subroutine Set_DS
    end interface
    type(element), intent(inout) :: elem
    class(element), pointer :: elem1
    integer :: ie, j,k, i
    real :: Re_1, Re, val, Mnorm, tt, tt1
    integer :: i_min, i_max, itest, itest1
    real, dimension(:), allocatable :: x1,x2,x3
    logical :: explicit_convection

    itest = -11
    !if( abs(elem%xc(1) - 32) < 2) itest = elem%i
    itest1 = -grid%nelem

    i = elem%i

    !explicit_convection = .true.   ! only for scalar equation (ndim = 1) makes sense
    explicit_convection = .false.

    !print*,' computing of convective (inviscid) terms',state%model%convective, elem%i
    if(state%model%convective) then

       if( explicit_convection ) then !scalar equation, explicit inviscid terms
          stop "explicit_convection commented"

          ! ! print*, 'Scalar linear convection-reaction'
          ! call LinearInviscidVolumes(elem )

          ! !print*,'explicit inviscid volumes',i
          ! call ExplicitInviscidVolumes(elem, Set_f_s_scalar)

          ! !print*,'explicit of inviscid edges', elem%i
          ! do j = 1, elem%flen
          !    k = elem%face(neigh,j)

          !    if(k > 0) then
          !       elem1 => grid%elem(k)
          !       !print*,'Inner edge',j,k,state%time%tau(:)
          !       call ExplicitInviscidEdge(elem, elem1, j, Set_NumFlux_scalar)

          !    else
          !       !print*,'Boundary edge',j,k,state%time%tau(:)
          !       call ExplicitInviscidEdge(elem, elem, j,  Set_NumFlux_scalar)
          !    endif

          ! enddo

       else ! fully (semi-)implicit

          !write(*,'(a6,30es12.4)') 'vec:',elem%vec(rhs,:)
          if(elem%i == itest) then
             print*, state%nlSolver%implicitly,elem%xc
             write(*,'(a8,i3,500es12.4)') 'vecA :',i, elem%vec(rhs,:)
          endif

          ! print*,'Matrix ani A'
          !call WriteMblock_Screene(elem%block(0) )       !call WriteMatrixA(0.)
          !write(*,'(a6,200es12.4)') 'wsol:',grid%elem(1)%w(0, :)

          !write(*,'(a35,i5,2es12.4)')'computing of inviscid volumes',elem%i, elem%xc(:)
          call ElementInviscidVolumes(elem, Set_f_s, Set_A_s)

          if(state%ST_Vc >0 .and. state%model%Re == 0.) then
             !print*,'computing volume stabilization for inviscid flows'
             call ElementVolumesStabil(elem )
          endif

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

                !print*,'Inner edge',j,k,elem%i, elem1%i, elem1%xc
                call ElementInviscidInnerEdge(elem, elem1, j, Set_Ppm, Set_f_s)
                !write(*,'(a8,i3,50es12.4)') 'vec:',i, elem%vec(rhs,:)
                if(elem%i == itest) write(*,'(a8,i3,500es12.4)') 'vecC :',i, elem%vec(rhs,:)

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

                   if(elem%i == itest) write(*,'(a8,i3,500es12.4)') 'vecD :',i, elem%vec(rhs,:)

                endif
                !write(*,'(a8,i3,50es12.4)') 'vec:',i, elem%vec(rhs,:)

             else
                !print*,'Input/output  edges',elem%i,j,k, 'iBC', elem%iBC(j)
                call ElementInviscidIOEdge(elem, j, Set_Ppm)
                !write(*,'(a8,i3,50es12.4)') 'vec:',i, elem%vec(rhs,:)

                if(elem%i == itest) write(*,'(a8,i3,500es12.4)') 'vecE :',i, elem%vec(rhs,:)

             endif

             !if(state%local_problem)  return

             !if(elem%i == itest1) then
             !   print*,'EDGE : = ',j
             !   call WriteMblock_Screene(elem%block(0) )
             !endif

          enddo

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

       end if  ! eximplicit convection

    endif   ! state%model%convection

    !tt1 = tt
    !call cpu_time(tt)
    !if(elem%i == itest1) print*,'BP,',tt,' Inv7', tt - tt1


    !write(*,'(a10,i5,100es12.4)') 'ED RHS a', 1, grid%elem(1)%vec(rhs, :)

    !print*,' computing of viscous and penalty terms for viscous flow',elem%i
    if(state%model%Re > 0.) then
       !if (elem%i == 1) print*,'viscous volume',elem%i
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
             if(elem%iBC(j) == 0 .and. (state%modelName == 'scalar'.or.state%modelName /= '2eqs')) then
                !HOMOGENEOUS Neumann
             else
                ! Dirchlet BC or fixed wall for NSe
                call ElementViscousBoundEdge(elem, j, Set_R_s, Set_K_sk)
                !write(43, *) elem%xc, elem%i, j
             endif


             !!endif
          endif
          !call WriteMblock(elem%block(0) )
       enddo

    endif

    !tt1 = tt
    !call cpu_time(tt)
    !if(elem%i == itest1) print*,'BP,',tt,' Inv8', tt - tt1


    !print*, 'reaction terms (at this moment only for scalar)'
    if( ((state%modelName == 'scalar' .or.state%modelName == '2eqs') .and. state%model%ireac > 0) &
         .or. state%modelName == 'pedes' ) &
         call ElementReactionVolumes(elem, Set_S, Set_DS)

    !print*, ' adding of the penalty leading to piecewise constant solution for vacuum state'
    !if( state%modelName == 'pedes' ) &
    !     call ElementGradientStabilization(elem)

    !print*,' adding of the source term', elem%i
    if(state%RHS_presented .and. &
         ( state%modelName == 'scalar' .or.state%modelName == '2eqs' .or. &
         (state%modelName == 'NSe' .and. state%type_IC .eq. 9) ) ) then
       !.or. state%modelName == 'wet_steam'  ) ) call ElementRHS(elem )

       !         print*, 'subdomainRHS:' , state%model%subdomainRHS

       if ( state%model%subdomainRHS ) then
          !         print*, 'ComputeOneElementTerms: impl' , state%nlSolver%implicitly
          call ElementSubdomainRHS( elem )
       else
          call ElementRHS(elem )
       endif
    endif

    !write(*,'(a8,i3,50es12.4)') 'vec:',i, elem%vec(rhs,:)

    if(elem%i == itest) write(*,'(a8,2i3,500es12.4)') 'vecZ :',i, size(elem%vec, 2) , elem%vec(rhs,:)

    !call WriteMblock_Screene(elem%block(0) )
    !call WriteMblock_Screene(elem%block(1) )
    !call WriteMblock_Screene(elem%block(2) )
    !call WriteMblock_Screene(elem%block(3) )
    !if(elem%i == 50) stop "d73d3d3d3"

    !tt1 = tt
    !call cpu_time(tt)
    !if(elem%i == itest1) then
    !   print*,'BP,',tt,' Inv9', tt - tt1
    !   print*,'...............................................'
    !endif

  end subroutine ComputeOneElementTerms

  !>  evaluation of all integrals and filling the matrix per elements
  !>  matrix \f$ {\bf C}( {\bf w}) \f$, vectors \f$ {\bf q}( {\bf w}) \f$,
  !> \f$ {\bf m}( {\bf w}) \f$
  subroutine ComputeElementsTerms(Set_f_s, Set_A_s, Set_Ppm, Set_R_s, Set_K_sk, Set_S, Set_DS)
    interface
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x, i)
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
         integer, intent(in) :: i
      end subroutine Set_f_s
      subroutine Set_A_s(ndimL, nbDim, Qdof, w, A_s, xi, i)
         integer, intent(in) :: Qdof, nbdim, ndimL
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: A_s
         ! matrices A_s in  -- " --
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
         integer, intent(in) :: i
      end subroutine
      subroutine Set_Ppm( ndimL, nbDim, Qdof, w, n, xi, Ppm, one_over_area, elem)
         import :: element
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: Ppm
                                               ! matrices Ppm in  -- " --
         real, dimension(1:Qdof, 1:nbDim), intent(in) :: n   ! outer normal
         real, dimension(1:Qdof, 1:nbDim),intent(in) ::  xi         ! node on the edge?
         real, intent(in), optional :: one_over_area
         type(element), intent(inout), optional :: elem
      end subroutine
      subroutine Set_R_s(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1     ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
      end subroutine Set_R_s
      subroutine Set_K_sk(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, K_sk, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1      ! inverse of Reynolds number
         real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,ndimL), intent(inout) :: K_sk
         real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
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
    real :: Re_1, Re, val, Mnorm, tt, tt1, val1, val2, val3
    integer :: i_min, i_max
    real, dimension(:), allocatable :: x1,x2,x3

    !!call cpu_time(tt)
    !!print*,'CT,',tt,' sub 1'

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

          if(state%modelName == 'pedes' ) then
             call ElementJumpVelocityIndicator(elem)
             !elem%rezid = (elem%rezid / elem%diam)**2
             elem%rezid = elem%rezid * elem%limit_par
          endif


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

    if(state%ST_Vc >0 .or.  state%ST_Ec> 0. ) then
       do i=1,2
          call SmoothResid()
       enddo
    endif

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


    !!call cpu_time(tt)
    !!print*,'CT,',tt,' sub 7'

    if(state%modelName == 'NSe' .or. state%modelName == 'wet_steam') call EvalOutputPressure()

    !call cpu_time(tt)
    !print*,'CT,',tt,' sub 8', state%nlSolver%implicitly,  grid%elem(1)%deg_plus
    !tt1 = tt

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

       !call cpu_time(tt)
       !print*,'CT,',tt,' sub .', elem%i

       ! do j=1,elem%Qdof
       !    call Set_Battery(elem%xi(0, j,1), elem%xi(0,j,2), val1, val2, val3)
       !    write(84, *) elem%xi(0, j,1), elem%xi(0,j,2), val1, val2, val3
       ! enddo
       ! do ie =1, elem%flen
       !    do j=1,elem%face(fGdof, ie)
       !       call Set_Battery(elem%xi(ie, j,1), elem%xi(ie,j,2), val1, val2, val3)
       !       write(85, *) elem%xi(ie, j,1), elem%xi(ie,j,2), val1, val2, val3
       !    enddo
       ! enddo

    enddo

    !stop "Battery writting 73d39iud"


!!  !$OMP END PARALLEL DO

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
            kvec = kvec + dof

            ! if(state%time%recompute_back >= 2 .and. &
            !      dot_product(elem%rhsST(k,1:elem%dof, l), elem%rhsST(k,1:elem%dof, l)) > 1E-15  ) &
            !      write(*,'(a4, 5i5,200es12.4)')'b:',elem%i,l,k, kvec+1, kvec+dof, b(kvec+1:kvec + dof)

         enddo !k
      enddo !l

      ivec = ivec + ndof
   end do
   !write(59,'(a6,100es12.4)' ) '#bb#',b(:)
   !write(59,*)'----------------',state%nsize,dof, ndof
 end subroutine FillVectorST

 !> Setting of the large vector from the actual STDG solution
 subroutine FillVectorST_from_wST(b)
   real, dimension(:), intent(inout) :: b
   class(element), pointer:: elem ! one element
   integer  :: kk, i, k, j

   kk = 0
   do i = 1, grid%nelem
      elem => grid%elem(i)
      do k = 1, elem%Tdof
         do j = 1, ndim
            b( kk + 1 : kk + elem%dof) = elem%wST(j,1:elem%dof,k)
            kk = kk + elem%dof
         enddo !j
      enddo !k
   enddo !i
   if(kk /= size(b, 1)) stop "Trouble in  FillVectorST_from_wST"

 end subroutine FillVectorST_from_wST

 !> Setting of the large vector from the actual STDG solution
 subroutine Fill_wST_from_VectorST(b)
   real, dimension(:), intent(inout) :: b
   class(element), pointer:: elem ! one element
   integer  :: kk, i, k, j

   kk = 0
   do i = 1, grid%nelem
      elem => grid%elem(i)
      do k = 1, elem%Tdof
         do j = 1, ndim
            elem%wST(j,1:elem%dof,k) = b( kk + 1 : kk + elem%dof)
            kk = kk + elem%dof
         enddo !j
      enddo !k
   enddo !i
   if(kk /= size(b, 1)) stop "Trouble in  FillVectorST_from_wST"

 end subroutine Fill_wST_from_VectorST

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

   call ComputeTerms(.false. )

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
 subroutine ComputeTerms( deg_plus )
   logical, intent(in) :: deg_plus
   integer :: i
   real :: tt, tt1

   call cpu_time(tt1)
   !write(*,'(a30, 2es12.4)') 'CT_starts :', tt - state%start_time

   grid%elem(:)%deg_plus = deg_plus

   ! clearing of the appropriate arrays
   if(state%nlSolver%implicitly) call ClearMatrixBlocks()
   call ClearVectorBlocks()

   !print*,'! limiting of velocity'
   !if(state%modelName == 'pedes' ) call Pedestrian_velocity_limiting( )
   !if(state%modelName == 'pedes' ) call Pedestrian_velocity_limiting2( )


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

   elseif(state%modelName == 'pedes' ) then    ! Pedestrian flow
      call ComputeElementsTerms(Set_f_s_pedes, Set_A_s_pedes, Set_Ppm_pedes, &
           Set_R_s_empty, Set_K_sk_empty, Set_S_pedes, Set_DS_pedes)

   elseif(state%modelName == 'swe' ) then    ! Shallow water equations
      call ComputeElementsTerms(Set_f_s_swe, Set_A_s_swe, Set_Ppm_swe, &
           Set_R_s_swe, Set_K_sk_swe, Set_S_swe, Set_DS_swe)

   elseif(state%modelName == 'porous' ) then    ! porous media flow
      !!call ComputeCapacityConductivity( )

      call ComputeElementsTerms(Set_f_s_empty, Set_A_s_empty, Set_Ppm_empty, &
           Set_R_s_porous, Set_K_sk_porous, Set_S_empty, Set_DS_empty)


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

   !call cpu_time(tt)
   !write(*,'(a30, 2es12.4)') 'CT_ends :', tt - state%start_time , tt - tt1

 end subroutine ComputeTerms


 !> evaluation of the matrix and vector blocks
 subroutine Compute_ONLY_ONE_ELEMENT_Terms(elem )
   class(element), target, intent(inout) :: elem
   integer :: j

   ! clearing of the appropriate arrays
   if(state%nlSolver%implicitly) then
      elem%block(0)%Mb(:,:) = 0. ! diagonal blocks
      do j = 1, elem%flen         ! off-diagonal blocks
         if(elem%face(neigh,j) > 0)  elem%block(j)%Mb(:,:) = 0.
      enddo
   endif

   elem%vec(rhs,:) = 0.

   !print*,'! limiting of velocity'
   !if(state%modelName == 'pedes' ) call Pedestrian_velocity_limiting( )
   !if(state%modelName == 'pedes' ) call Pedestrian_velocity_limiting2( )


   ! setting of the corresponding fluxes
   if(state%modelName == 'scalar') then        ! 2D scalar equation
      call ComputeOneElementTerms(Set_f_s_scalar, Set_A_s_scalar, Set_Ppm_scalar, &
           Set_R_s_scalar, Set_K_sk_scalar, Set_S_scalar, Set_DS_scalar, elem)

   elseif(state%modelName == 'NSe' ) then    ! 2D Euler and Navier-Stokes equations
      call ComputeOneElementTerms(Set_f_s_Euler, Set_A_s_Euler, Set_Ppm_Euler, &
           Set_R_s_NS, Set_K_sk_NS, Set_S_empty, Set_DS_empty, elem)

   elseif(state%modelName == '2eqs') then    ! 2 scalar equations
      call ComputeOneElementTerms(Set_f_s_scalar, Set_A_s_scalar, Set_Ppm_scalar, &
           Set_R_s_scalar, Set_K_sk_scalar, Set_S_scalar, Set_DS_scalar, elem)

   elseif(state%modelName == 'wet_steam' ) then    ! 2D wet_stem
      call ComputeOneElementTerms(Set_f_s_WS, Set_A_s_WS, Set_Ppm_WS, &
           Set_R_s_WS, Set_K_sk_WS, Set_S_empty, Set_DS_empty, elem)

   elseif(state%modelName == 'pedes' ) then    ! Pedestrian flow
      call ComputeOneElementTerms(Set_f_s_pedes, Set_A_s_pedes, Set_Ppm_pedes, &
           Set_R_s_empty, Set_K_sk_empty, Set_S_pedes, Set_DS_pedes, elem)

   elseif(state%modelName == 'swe' ) then    ! Shallow water equations
      call ComputeOneElementTerms(Set_f_s_swe, Set_A_s_swe, Set_Ppm_swe, &
           Set_R_s_swe, Set_K_sk_swe, Set_S_swe, Set_DS_swe, elem)

   elseif(state%modelName == 'swe' ) then    ! Shallow water equations
      call ComputeOneElementTerms(Set_f_s_empty, Set_A_s_empty, Set_Ppm_empty, &
           Set_R_s_porous, Set_K_sk_porous, Set_S_empty, Set_DS_empty, elem)


  ! elseif(nbDim == 2 .and. ndim == 6) then    ! RANS - k-omega  model
  !    call ComputeElementTerms(Set_f_s_Turb2e, Set_A_s_Turb2e, Set_Ppm_Turb2e, &
  !         Set_R_s_Turb2e, Set_K_sk_Turb2e, Set_S_empty, Set_DS_empty)

  ! elseif(nbDim == 3 .and. ndim == 5) then    ! 3D Euler and Navier-Stokes equations
  !    call ComputeElementTerms(Set_f_s_Euler3D, Set_A_s_Euler3D, Set_Ppm_Euler3D, &
  !         Set_R_s_NS3D, Set_K_sk_NS3D, Set_S_empty, Set_DS_empty )


   else
      print*,'Compute_ONLY_ONE_ELEMENT_Terms not implemented for state%modelName',state%modelName
   endif

   ! adding of time derivative terms
   call TimeDerivativeVector_ONLY_ONE_ELEM(elem )

   !do i=1,grid%nelem
   !   !write(*,'(a4,i5,20es12.4)') 'iw: ',i, grid%elem(i)%w(0,:)
   !   !print*, 'EEDE',size (grid%elem(i)%block(0)%Mb(:,:), 1)
   !   call WriteMblock(grid%elem(i)%block(0) )
   !   write(*,'(a10,i5,100es12.4)') 'ED RHS', i, grid%elem(i)%vec(rhs, :)
   !enddo

!  stop 'stopped in the end of ComputeTerms'

   !call cpu_time(tt)
   !write(*,'(a30, 2es12.4)') 'CT_ends :', tt - state%start_time , tt - tt1

 end subroutine Compute_ONLY_ONE_ELEMENT_Terms

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
   logical :: deg_plus

   deg_plus = .false.

   allocate(accum(1:state%space%max_dof) )

   state%nlSolver%implicitly = .false.
   call ComputeTerms(deg_plus )
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
 !> if implicitly == False -> vector F(u) saved in elem%vec(rhs,:)
 !> if implicitly == True  -> sets the matrix C(*,*) elem%vec(rhs = q(u),
 subroutine ComputeSTDGM_Terms( deg_plus )
   logical, intent (in) :: deg_plus
   class(element), pointer :: elem
   class(Time_rule), pointer :: T_rule
   real, dimension(:,:), allocatable :: T_mat
   integer :: i, j, k, l, r, kk, rr, ll
   integer :: alpha, Qdeg, dof, Tdof, s_dim, f_dof, f_Tdof, f_dim
   integer :: wTdof, wdof  !NEW for adaptation - smaller than dof/Tdof if elem%deg_plus=TRUE
   integer :: m, mm, n, nn
   real :: cTime, tt, tt1
   real :: val, val_vec
   real :: local_eta
   logical ::iprint
   integer :: dofA, itest

   ! elem => grid%elem(100)
   ! allocate(smaz(1:elem%Tdof * elem%dof * ndim, 1:elem%Tdof * elem%dof * ndim ) )
   ! allocate(smaz1(1:elem%Tdof * elem%dof * ndim, 1:elem%Tdof * elem%dof * ndim ) )
   ! smaz = 0.
   ! smaz1 = 0.

   !print*, 'ComputeSTDGM_Terms with implicitly = ' ,state%nlSolver%implicitly, 'deg_plus:', deg_plus
   call cpu_time(tt1)

   cTime = state%time%ctime

   local_eta = 1. / state%time%tau(1)

    ! adding of the term in front of the time derivative
   if(state%model%varying_time_term) then
      s_dim =  maxval( grid%elem(:)%dof_plus)*ndim
      allocate( T_mat(1: s_dim, 1:s_dim) )
   endif

   !   !F@R control phi(Tdof + 1) ma koreny v int uzlech

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

         state%linSolver%precond_update = .true.

         do alpha = 1, Qdeg ! temporarily max_Tdof =  max time quadrature nodes
            !print*,'#E#E#:',alpha,Qdeg
            do i = 1, grid%nelem
               elem => grid%elem(i)
               if (Qdeg /= elem%TQnum) then
                  !F@R Verify if it is OK, some nodes could be in wrong position
                  stop 'Verify if it is OK, some nodes could be in wrong position'
               endif
               ! save the wST space-time solution in quadrature index alpha to w
               call Transfer_wST_to_w_Elem(elem , alpha, Qdeg)
            enddo

            !we have to run ComputeTerms() in the time-quadrature nodes
            !cTime = state%time%ctime
            state%time%ctime = state%time%ttime + state%time%tau(1) * T_rule%lambda(alpha)
            !   print*, 'ctime:' , state%time%ctime, cTime
            !   print*, '--------------','ctime:', state%time%ctime, state%nlSolver%implicitly
            call ComputeTerms(deg_plus )
            state%time%ctime = cTime


            do i =1, grid%nelem
               elem => grid%elem(i)
               Tdof = elem%Tdof
               dof = elem%dof
               s_dim = ndim *dof

               do l = 1, Tdof
                  !diag blocks of blockST
                  do r = 1, Tdof

                     val = T_rule%phi(l,alpha)*T_rule%phi(r,alpha) * T_rule%weights(alpha)  !phi_l(t_alpha) * phi_r(t_alpha)* w_alpha

                     elem%blockST(0)%Mb((l-1)*s_dim + 1 : l*s_dim, (r-1)*s_dim +1 : r*s_dim) = &
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


                  ! ! VD: the following is not necessary, elem%rhsST is set in implicitly = F
                  ! !vec -> rhsST
                  ! !write(*,'(a8,12es12.4)') 'elem%vec', elem%vec(rhs,1:3)
                  ! ! for now G_rule%weights used for time quadrature
                  ! val_vec = T_rule%phi(l,alpha) * T_rule%weights(alpha)

                  ! do k = 1, ndim
                  !    elem%rhsST(k, 1:dof, l) = elem%rhsST(k, 1:dof, l)  &
                  !         + ( elem%vec(rhs,(k-1)*dof + 1 : k*dof) * val_vec )
                  ! enddo !k

               enddo !l = 1,Tdof


               ! adding of the term in front of the time derivative - time deriv term itself
               if(state%model%varying_time_term) then

                  ! computing of the matrix in front of the time derivative term at t_alpha
                  !print*,'Compute_Time_deriv A', s_dim, elem%dof
                  call Compute_Time_deriv(elem, .false., s_dim,  T_mat(1:s_dim,1:s_dim) )

                  !val = T_mat(1,1)  ! in implicitly mode,  T_mat(1,1) =  T_mat(i,i), i=1,2,3,...
                  !if( sqrt(T_mat(1,1)**2 + T_mat(2,2)**2 + T_mat(3,3)**2 ) > 1E-1) then
                  !   do l=1,s_dim
                  !      write(*,'(a8,3i5, 120es12.4)') 'T_mat', elem%i, s_dim, l, T_mat(l,1:s_dim)
                  !   enddo
                  !   print*
                  !endif

                  ! Time Penalty Var2
                  !T_mat(1:ndim*dof,1:ndim*dof)  = 0.
                  !do l=1,ndim*dof
                  !   !T_mat(l,l) = 1.
                  !   T_mat(l,l) = val
                  !enddo

                  do l=1,Tdof
                     do r = 1, Tdof
                        ll = (l -1) * s_dim
                        rr = (r -1) * s_dim

                        val = T_rule%phi(l,alpha)*T_rule%Dphi(r,alpha) * T_rule%weights(alpha)* local_eta

                        elem%blockST(0)%Mb( ll + 1 : l*s_dim, rr +1 : r*s_dim)  &
                             = elem%blockST(0)%Mb( ll + 1 : l*s_dim, rr +1 : r*s_dim)  &
                             + val * T_mat(1:s_dim,1:s_dim)

                        !!if(elem%i == 100) smaz( ll + 1 : l*s_dim, rr +1 : r*s_dim) = &
                        !smaz( ll + 1 : l*s_dim, rr +1 : r*s_dim) + val * T_mat(1:s_dim,1:s_dim)

                     enddo  ! r
                  enddo ! l


               endif ! state%model%varying_time_term

            enddo !i

         end do !alpha =  1,Qdeg

         ! itest = 1721
         ! elem => grid%elem(itest)
         ! do r=3, 0, -1
         !    if(r == 0 .or. elem%face(neigh, r) > 0) then
         !       do l=1, elem%dof*elem%Tdof
         !          write(*,'(a8,2i5,300es11.3)') 'bST:',r, l,elem%blockST(r)%Mb(l, :)
         !       enddo
         !       print*
         !    endif

         ! enddo
         ! print*,'_________________________________'

         ! adding of the term in front of the time derivative - jump term,  w_{m-1}^+ part:
         if(state%model%varying_time_term) then
            do i = 1, grid%nelem
               elem => grid%elem(i)

               Tdof = elem%Tdof
               dof = elem%dof
               s_dim = ndim *dof

               ! save the wST space-time solution in quadrature index alpha to w
               call Transfer_wST_to_w_Elem(elem , -1, Qdeg)

               ! computing of the matrix in front of the time derivative term at t_{m-1}^+
               !print*,'Compute_Time_deriv B', s_dim, elem%dof
               call Compute_Time_deriv(elem, .false., s_dim, T_mat(1:s_dim,1:s_dim) )

               ! Time Penalty Var0
               T_mat(1:s_dim,1:s_dim)  = 0.
               do l=1,s_dim
                  T_mat(l,l) = 1.
               enddo


               !if( dot_product( T_mat(1, 1:s_dim),  T_mat(1, 1:s_dim)) > 1E-2) then
               !   do l=1, s_dim
               !      write(*,'(a6, i5, i3, 20es10.2)') 'T_mat', l, s_dim, T_mat(l, :)
               !   enddo
               !   print*,'I:', i, elem%xc(:), dot_product( T_mat(1, 1:s_dim),  T_mat(1, 1:s_dim))
               !   print*
               !   write(91, *)  elem%xc(:), dot_product( T_mat(1, 1:s_dim),  T_mat(1, 1:s_dim))
               !endif

               do l=1,Tdof
                  do r = 1, Tdof
                     ll = (l -1) * s_dim
                     rr = (r -1) * s_dim

                     val = T_rule%phi(l,-1)*T_rule%phi(r,-1) * local_eta

                     elem%blockST(0)%Mb( ll + 1 : l*s_dim, rr +1 : r*s_dim)  &
                          = elem%blockST(0)%Mb( ll + 1 : l*s_dim, rr +1 : r*s_dim)  &
                          + val * T_mat(1:s_dim,1:s_dim)

                     !if(elem%i == 100) then
                     !   write(*,'(a12, 5i5, 20es12.4)') &
                     !        'd30i430i40', l,r,ll, rr, T_rule%phi(l,-1)*T_rule%phi(r,-1)
                     !endif

                     ! if(elem%i == 100) smaz(ll + 1 : l*s_dim, rr +1 : r*s_dim) = &
                     !      smaz(ll + 1 : l*s_dim, rr +1 : r*s_dim)+ val * T_mat(1:s_dim,1:s_dim)

                     ! if(elem%i == 100) smaz1(ll + 1 : l*s_dim, rr +1 : r*s_dim) = &
                     !      smaz1(ll + 1 : l*s_dim, rr +1 : r*s_dim)+ val * T_mat(1:s_dim,1:s_dim)


                  enddo  ! r
               enddo ! l
            enddo ! i=1,grid%nelem

         endif !(state%model%varying_time_term

         ! itest = 1721
         ! elem => grid%elem(itest)
         ! do r=0, 0
         !    do l=1, elem%dof*elem%Tdof
         !       write(*,'(a8,2i5,300es11.3)') 'bST:',r, l,elem%blockST(r)%Mb(l, :)
         !    enddo
         !    print*
         ! enddo
         ! print*,'##############################################'


         !stop "9u49u49of4"

         ! ! VD: the following is not necessary, elem%rhsST is set in implicitly = F
         !adding the timejump part
         !do i = 1, grid%nelem
         !   elem => grid%elem(i)
         !   call Elem_wSTfinToRhsST( grid%elem(i) , state%nlSolver%implicitly)
         !enddo !i

         !implicitly = .false. => rshST = F(w)


!          ! comparison of matrices
!          if(state%time%recompute_back <= 2 ) then
!             do i = 100, 100 !grid%nelem
!                elem => grid%elem(i)
!                if(mod(i, 100) /= 0) goto 20
!                Tdof = elem%Tdof
!                dof = elem%dof
!                s_dim = ndim *dof


!                if(state%model%varying_time_term) then
!                   do l=1, Tdof * dof
!                      write(*,'(a8, 4i5,200es12.4)')'B_ST1:',elem%i,dof, Tdof,l, smaz1(l, :) !elem%blockST(0)%Mb(l, :)
!                   enddo
!                   print*
!                   !do l=1, Tdof * dof
!                   !   write(*,'(a8, 4i5,200es12.4)')'B_ST:',elem%i,dof, Tdof,l, smaz(l, :) !elem%blockST(0)%Mb(l, :)
!                   !enddo
!                   deallocate(smaz, smaz1)
!                else
!                   allocate(T_mat(1:dof*Tdof*ndim, 1:dof*Tdof*ndim) )
!                   do m =1,Tdof
!                      mm = (m-1) * s_dim
!                      do n = 1,Tdof
!                         nn = (n-1) * s_dim
!                         do k = 0,s_dim-1,dof
!                            do l = 0, s_dim-1, dof

!                               !val =  time%refTimeMatrix%Mb(m,n)

!                               val = 0. * dot_product(T_rule%weights(1:Qdeg), &
!                                    T_rule%Dphi(n,1:Qdeg) *T_rule%phi(m,1:Qdeg)) &
!                                    +   T_rule%phi(m,-1)*T_rule%phi(n,-1)

!                               print*,'d30i430i40', val, time%refTimeMatrix%Mb(m,n)
!                               !     val-time%refTimeMatrix%Mb(m,n) , &
!                               !     dot_product(T_rule%weights(1:Qdeg), &
!                               !     T_rule%Dphi(n,1:Qdeg) *T_rule%phi(m,1:Qdeg))

!                               T_mat(mm + k + 1 : mm+ k + dof, nn + l + 1 : nn+ l + dof) &
!                                    = local_eta *  elem%Mass%Mb(1:dof, 1:dof)  * val


!                            enddo
!                         enddo
!                      enddo
!                   enddo


!                   do l=1,Tdof*dof
!                      write(*,'(a8, 4i5,200es12.4)')'B_ST:',elem%i,dof, Tdof, l,  T_mat(l, :) !&
!                      !elem%blockST(0)%Mb(l, :) + T_mat(l, :)
!                   enddo

!                   deallocate(T_mat)
!                endif
! 20             continue
!             enddo
!          if(state%time%recompute_back >= 2 ) stop "(*&^%$##%^&*()"
!          endif


      else

!         print*, 'ComputeSTDGM_Terms with implicitly = FALSE called'
!         print*,
         iprint = .false.

!         if (grid%elem(1)%deg_plus) iprint = .true.
!
!         if (iprint) print*, 'ComputeSTDGM_Terms with implicitly = FALSE called', grid%elem(1)%deg_plus
!         if (iprint) print*, 'size of vec ' , size(grid%elem(1)%vec(rhs,:))

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

               call Transfer_wST_to_w_Elem(elem , alpha, Qdeg)
               !call Transfer_wST_to_w_Elem(elem , alpha, elem%TQnum)
            enddo

            !we have to run ComputeTerms() in the time-quadrature nodes
            state%time%ctime = state%time%ttime + state%time%tau(1) * state%time%T_rule(Qdeg)%lambda(alpha)

            !      print*, '-----------------------','ctime:', state%time%ctime, state%nlSolver%implicitly
            !      write(*,'(a6,12es12.4)') 'w(t):',grid%elem(1)%w(0,:)*2**0.5
            call ComputeTerms( deg_plus )

            do i =1, grid%nelem
               elem => grid%elem(i)

               dofA = elem%dof
               if(elem%deg_plus) dofA = elem%dof_plus

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

               s_dim = ndim * dof

               !vec -> rhsST
               do k = 1,ndim
                  do l = 1, Tdof
                     val_vec = T_rule%phi(l,alpha) * T_rule%weights(alpha)

                     elem%rhsST(k, 1:dof, l) = elem%rhsST(k, 1:dof, l)  &
                          +  val_vec *  elem%vec( rhs, (k-1)*dof + 1 : k*dof )   	! dofA

                  enddo !l
               enddo !k


               ! adding of the term in front of the time derivative - time deriv term itself
               if(state%model%varying_time_term) then

                  ! computing of the matrix in front of the time derivative term at t_alpha
                  !print*,'Compute_Time_deriv C', s_dim, elem%dof_plus
                  call Compute_Time_deriv(elem, elem%deg_plus, s_dim,  T_mat(1:s_dim,1:s_dim) )

                  do k = 1,ndim
                     do l = 1, Tdof

                        val_vec = T_rule%phi(l,alpha) * T_rule%weights(alpha) * local_eta

                        do j=1, dof
                           elem%rhsST(k, j, l) = elem%rhsST(k, j, l)  &
                                -  val_vec  &
                                * dot_product(T_rule%Dphi(1:wTdof,alpha) , &
                                matmul( T_mat( (k-1)* dof + j, (k-1)*dof+1:(k-1)*dof + wdof ),&
                                elem%wST(k, 1:wdof, 1:wTdof) ) )
                        enddo
                     enddo  ! l=1,Tdof
                  enddo ! k=1,ndim
               endif ! state%model%varying_time_term

            enddo !i
         end do !alpha =  1,Qdeg

         !    write(*,'(a15, 12es12.4)') 'wST:' , grid%elem(1)%wST(1,:,:)
!            write(*,'(a15, 12es12.4)') 'F^m: aft vec ' , grid%elem(1)%rhsST(1,:,:)

         ! adding the time derivative and the the w_{m-1}^+ part:
         do i =1, grid%nelem
            elem => grid%elem(i)

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

            ! original subroutine for the adding the time derivative and the the w_{m-1}^+ part
            if(.not.  state%model%varying_time_term ) then
               do k = 1, ndim
                  do l = 1,Tdof
                     do j = 1,dof   ! dofA
                        !NEW for adaptation
                        val_vec = local_eta * dot_product( time%refTimeMatrix%Mb(l,1:wTdof) , &
                             matmul( elem%Mass%Mb(j, 1:wdof), elem%wST(k,1:wdof,1:wTdof) ))

                        elem%rhsST(k, j, l) = elem%rhsST(k, j, l) - val_vec

                     enddo !j
                  enddo !l

               enddo !k
            endif

            ! adding of the term in front of the time derivative -  jump term,  w_{m-1}^+ part:
            if(state%model%varying_time_term) then
               call Transfer_wST_to_w_Elem(elem , -1, Qdeg)

               ! computing of the matrix in front of the time derivative term at t_{m-1}^+
               s_dim = ndim * dof
               !print*,'Compute_Time_deriv D', s_dim, elem%dof_plus
               !call Compute_Time_deriv(elem, elem%deg_plus, s_dim, T_mat(1:s_dim, 1:s_dim) )

               ! Time Penalty Var1
               T_mat(1:s_dim,1:s_dim)  = 0.
               do l=1,s_dim
                  T_mat(l,l) = 1.
               enddo

               do k = 1,ndim
                  do l = 1, Tdof

                     val_vec = T_rule%phi(l,-1) * local_eta

                     do j=1, dof
                        elem%rhsST(k, j, l) = elem%rhsST(k, j, l)  &
                             -  val_vec  &
                             * dot_product(T_rule%phi(1:wTdof,-1) , &
                             matmul( T_mat( (k-1)* dof + j, (k-1)*dof+1:(k-1)*dof + wdof ),&
                             elem%wST(k, 1:wdof, 1:wTdof) ) )
                     enddo
                  enddo  ! l=1,Tdof
               enddo ! k=1,ndim

            endif   !state%model%varying_time_term

            !if(elem%i >= 22 .and. state%time%recompute_back == 2  ) then
            !   do l=1, size(elem%rhsST, 3)
            !      write(*,'(a10,2i5, 300es12.4)') 'rhsST E', alpha, l,&
            !           elem%rhsST(1, :, l)
            !   enddo
            !endif

         enddo !i

!         write(*,'(a15, 12es12.4)') 'F^m: aft mass' , grid%elem(1)%rhsST(1,:,:)

         ! the w_{m-1}^- part:   "state%model%varying_time_term" given inside
         do i = 1, grid%nelem
            elem => grid%elem(i)
            call Elem_wSTfinToRhsST( grid%elem(i) , state%nlSolver%implicitly)


            !if(elem%i >= 22 .and. state%time%recompute_back == 2 ) then
            !   do l=1, size(elem%rhsST, 3)
            !      write(*,'(a10,2i5, 300es12.4)') 'rhsST F', alpha, l,&
            !           elem%rhsST(1, :, l)
            !   enddo
            !endif
         enddo !i

         !if(state%time%recompute_back == 2 )stop "3d3ed399"

      endif   !state%nlSolver%implicitly
      class default
      stop 'ComputeSTDGM_Terms, For TimeTDG_t only!!!'
   end select
 end associate ! time

 ! putting back the original value
 state%time%ctime = cTime

 if(state%model%varying_time_term) deallocate( T_mat)

 call cpu_time(tt)
 !  print*, '________________________________________'
 !write(*,'(a36, l4, l4, f8.2, f12.4)') &
 !     '#CPU# ComputeSTDGM: implic, dof++: ', state%nlSolver%implicitly, elem%deg_plus, &
 !     tt - state%start_time, tt - tt1
 end subroutine ComputeSTDGM_Terms


 !> try to avoin nonphysical solution, negative density, pressure, etc.
 subroutine Pedes_Empty_domain (finish)
   logical, intent(inout) :: finish
   class(element), pointer :: elem ! elem = element
   real, dimension(:,:), pointer:: phi ! local store arrays
   real, dimension(:), pointer:: weights ! local store arrays
   real, dimension(:), allocatable :: wi ! w in integ. nodes
   integer :: i, dof, Qdof
   real :: density, density_tot, density_tot1

   density_tot = 0.
   density_tot1 = 0.
   !if (state%time%disc_time == 'STDG') then
   do i=1,grid%nelem
      elem => grid%elem(i)

      dof = elem%dof
      Qdof = elem%Qdof

      phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)
      weights => state%space%V_rule(elem%Qnum)%weights(1:Qdof)

      allocate(wi(1:Qdof) ) ! w in integ. nodes

      wi(1:Qdof) = matmul(elem%w(0,1:dof), phi(1:dof, 1:Qdof) )
      density = dot_product(wi(1:Qdof), weights(1:Qdof) )

      if(density > 0.)             density_tot = density_tot + density * elem%area
      if(density > state%model%Pr) density_tot1 = density_tot1 + density * elem%area
      deallocate(wi)
   enddo

   !if( density_tot < state%model%Pr * state%space%domain_volume) finish = .true.
   if( density_tot < 2. ) finish = .true.

   ! output file
   open(11, file = 'pedes_history', status = 'unknown', position = 'append')
   if(state%time%iter_loc == 1) &
        write(11, '(a6, 8a12)') 'iter', 'ttime', 'dens_tot', 'dens_tot1', 'dens_min'

   write(11,'(i6, 8es12.4)')  state%time%iter_loc, state%time%ttime, &
        density_tot, density_tot1, state%model%Pr * state%space%domain_volume
   close(11)


 end subroutine Pedes_Empty_domain

 !> try to avoin nonphysical solution, negative density, pressure, etc.
 subroutine AvoidNonphysicalValues( )
   class(element), pointer :: elem ! elem = element
   integer :: i

    if(state%modelName == 'pedes' ) then
       if (state%time%disc_time == 'STDG') then
          do i=1,grid%nelem
             elem => grid%elem(i)
             call PedestrianFlow_AvoidVacuumSTDGM(elem)
          enddo
       else
          do i=1,grid%nelem
             elem => grid%elem(i)
             call PedestrianFlow_AvoidVacuum(elem)
          enddo

       endif
    endif

 end subroutine AvoidNonphysicalValues

  !> test the possible vocuum and its elimination by the adding of an average
  subroutine PedestrianFlow_AvoidVacuumSTDGM(elem)
    type(element), intent(inout) :: elem
    class(Time_rule), pointer :: T_rule
    real :: rho_min, rho_mean, rho_minG, rho_meanG, theta
    integer :: it, alpha, k, dof, Qdeg, Qdof,  Tdof

    dof = elem%dof
    theta = 1.

    Qdeg = elem%TQnum   !!state%time%Qnum
    Tdof = elem%Tdof
    T_rule => state%time%T_rule(Qdeg)

    do it = 1, 3
       rho_minG = 10000.
       rho_meanG = 0.

       do alpha = 1, Qdeg
          call Transfer_wST_to_w_Elem(elem , alpha, Qdeg)

          call PedestrianFlow_MinimalDensity(elem, rho_min, rho_mean)

          !if(elem%i == 42 ) &
          !  call PlotElem_D_Function3D(elem%i*10+alpha, elem,  elem%dof, elem%w(0, 1:dof) )

          !!print*,'###ede635dej', elem%w(0,1), rho_min, rho_mean, rho_mean / sqrt(2.)
          rho_minG = min(rho_minG, rho_min)
          rho_meanG =  rho_meanG +  rho_mean *  T_rule%weights(alpha)
       enddo


       if(rho_minG <  state%model%Pr) then
          write(*,'(a10,2i5, 2es12.4, a2, 6es12.4)')  &
            'rho_max:',elem%i, it, elem%xc(:), '|',state%model%Pr, rho_minG, rho_meanG, theta

          !print*,'ED#ED', size(elem%wST, 1), size(elem%wST, 2), size(elem%wST, 3)
          !print*,':::::',Qdeg, Tdof

          elem%wST(1:ndim,2:dof,1:1   ) = 0.75 * elem%wST(1:ndim,2:dof, 1:1)
          elem%wST(1:ndim,2:dof,2:Tdof) = 0.75 * elem%wST(1:ndim,2:dof,2:Tdof)

          !if(elem%i == 42 ) then
          !   do alpha = 1, Qdeg
          !      call Transfer_wST_to_w_Elem(elem , alpha, Qdeg)
          !      call PlotElem_D_Function3D(1000+elem%i*10+alpha, elem,  elem%dof, elem%w(0, 1:dof) )
          !   enddo
          !endif

       !stop '  PedestrianFlow_AvoidVacuumSTDGM dey34td3eij'
       else
          goto 100
       endif

    enddo
100 continue

  end subroutine PedestrianFlow_AvoidVacuumSTDGM

  !> test the possible vocuum and its elimination by the adding of an average
  subroutine PedestrianFlow_AvoidVacuum(elem)
    type(element), intent(inout) :: elem
    real :: rho_min, rho_mean, theta
    integer :: it, k, dof

    dof = elem%dof
    theta = 1.

    do it = 1, 1
       call PedestrianFlow_MinimalDensity(elem, rho_min, rho_mean)

       !!print*,'###ede635dej', elem%w(0,1), rho_min, rho_mean, rho_mean / sqrt(2.)

       if(rho_min <  state%model%Pr .or. it > 1) &
            write(*,'(a10,2i5, 2es12.4, a2, 6es12.4)')  &
            'rho_max:',elem%i, it, elem%xc(:), '|',state%model%Pr, rho_min, rho_mean, theta


       ! P_0 approximation
       elem%w(0,         2 :   dof ) = 0.
       elem%w(0,   dof + 2 : 2*dof ) = 0.
       elem%w(0, 2*dof + 2 : 3*dof ) = 0.

       if(rho_min >=  state%model%Pr) goto 100

       !call PlotElem_D_Function3D(981, elem,  elem%dof, elem%w(0, 1:dof) )
       !call PlotElem_D_Function3D(982, elem,  elem%dof, elem%w(0, dof+1: 2*dof) )
       !call PlotElem_D_Function3D(983, elem,  elem%dof, elem%w(0, 2*dof+ 3:dof) )
       !call PlotElem_D_Function3D(984, elem,  elem%dof, elem%w(0, dof+1: 2*dof) / elem%w(0, 1:dof))
       !call PlotElem_D_Function3D(985, elem,  elem%dof, elem%w(0, 2*dof+ 3:dof) / elem%w(0, 1:dof))


       if(rho_mean == rho_min ) then
          theta = 1.
       else
          theta = min(1., rho_mean /(rho_mean - rho_min ) )
       endif

       !theta = (rho_mean - state%model%Pr)  / (rho_mean - rho_min)
       theta = 0.75

       !!! modification of the density and the momentum
       ! ONLY DENSITY !!
       do k = 1, 1 ! ndim

       !   elem%w(0, (k-1)*dof + 1 ) = theta * elem%w(0, (k-1)*dof + 1 ) + (1.-theta)*rho_mean /sqrt(2.)

       !   elem%w(0, (k-1)*dof + 2 :k*dof ) = theta * elem%w(0, (k-1)*dof + 2: k*dof )
       enddo

       !call PlotElem_D_Function3D(991, elem,  elem%dof, elem%w(0, 1:dof) )
       !call PlotElem_D_Function3D(992, elem,  elem%dof, elem%w(0, dof+1: 2*dof) )
       !call PlotElem_D_Function3D(993, elem,  elem%dof, elem%w(0, 2*dof+ 3:dof) )
       !call PlotElem_D_Function3D(994, elem,  elem%dof, elem%w(0, dof+1: 2*dof) / elem%w(0, 1:dof))
       !call PlotElem_D_Function3D(995, elem,  elem%dof, elem%w(0, 2*dof+ 3:dof) / elem%w(0, 1:dof))

       ! stop "PedestrianFlow_AvoidVacuum"

    enddo

100 continue


  end subroutine PedestrianFlow_AvoidVacuum

  !> compute the minimal density for the pedestrian flow in volume and edge integ nodes
  subroutine PedestrianFlow_MinimalDensity(elem, rho_min, rho_mean)
    type(element), intent(inout) :: elem
    real, intent(inout) :: rho_min    ! minimal value
    real, intent(inout) :: rho_mean   ! mean value
    real, dimension(:,:), allocatable :: wi ! w recomputed  in integ nodes
    integer :: ie, Qdof

    rho_min = 1E+30

    Qdof = max(elem%Qdof, maxval( elem%face(fGdof,:) ) )

    allocate(wi(1:Qdof,1:ndim) )

    call Eval_w_Elem(elem, wi(1:elem%Qdof,1:ndim) )
    rho_min = minval(wi(1:Qdof, 1) )

    rho_mean = dot_product(wi(1:Qdof, 1), state%space%V_rule(elem%Qnum)%weights(1:elem%Qdof) )

    do ie = 1, elem%flen
       call Eval_w_Edge(elem, ie,  wi(1:elem%face(fGdof, ie), 1:ndim), .false. )
       rho_min = min(rho_min, minval (wi(1:elem%face(fGdof, ie), 1) ) )
    enddo

    deallocate(wi)

  end subroutine PedestrianFlow_MinimalDensity

  !> deallocation of arrays for the pedestrian eikonal equation
  subroutine DeallocatePedestrianEikonal( )
    use BRAlgorithm
    use vertqueue

    call BRDeallocate()

    call QueueDeallocate()

  end subroutine DeallocatePedestrianEikonal

  !> solution of the pedestrian eikonal equation
  subroutine PedestrianEikonalEquation(use_default_initc )
    logical, intent(inout) :: use_default_initc
    real, dimension(:,:), allocatable :: density   ! second component is the weght
    real, dimension(:,:), allocatable :: velocity, velocity_nodes
    integer, dimension(:,:), allocatable :: lnd
    class(element), pointer :: elem
    type(volume_rule), pointer :: V_rule
    type(Gauss_rule), pointer :: G_rule
    type(Bound_edge_t), pointer :: b_edge
    real :: rho, val, val1, val2, xi(2), ti(2)
    integer :: npoin, nelem, nbelm, max_nbp, i, ib, j, j1, k, k1, k2, k3, Qdof


    npoin = grid%npoin
    nelem = grid%nelem
    nbelm = grid%nbelm
    max_nbp = grid%max_nbp

    allocate(density(1:npoin, 1:2), velocity(1:nelem, 1:2), velocity_nodes(1:npoin, 0:2) )

    density(:, :) = 0.
    do i=1,nelem
       elem => grid%elem(i)
       rho = elem%w(0, 1) / sqrt(0.5)  ! mean value of the desnity on the element
       do j=1, elem%flen
          k = elem%face(idx, j)

          density(k, 1) = density(k, 1) + rho * elem%area
          density(k, 2) = density(k, 2) +  elem%area
       enddo
    enddo

    ! the "normalization"
    density(1:npoin, 1) =  density(1:npoin, 1)  / density(1:npoin, 2)

   ! ! OUTPUT
    ! open(22, file="eikonal-density.data", status='replace' )
    ! write(22, *) npoin
    ! do i=1,npoin
    !    write(22, *) density(i, 1)
    ! enddo
    ! close(22)

    ! local store
    allocate( lnd(1:nelem, 1:3) )
    do i=1,nelem
       lnd(i, 1:3) = grid%elem(i)%face(idx, 1:3)
    enddo


    call EikonalVelocity(npoin, nelem, max_nbp, grid%x(1:npoin, 1:2), lnd(1:nelem, 1:3), &
         grid%loc_ibp(1:npoin, 0: max_nbp), &
         density(1:npoin, 1), velocity(1:nelem, 1:2),  use_default_initc)


    ! recomputation of the eikonal velocity to integ nodes
    call EikonalVelocityIntegNodes(nelem, npoin, velocity(1:nelem, 1:2), velocity_nodes(1:npoin, 0:2) )


    use_default_initc = .false.   ! for the next time the precomputed quantities will be used


    deallocate(density, velocity, lnd, velocity_nodes)

    !stop'e4de6de3iudh33'

  end subroutine PedestrianEikonalEquation


  !> eikonal velocity: recomputation piecewise linear approximation onto integ nodes
  subroutine EikonalVelocityIntegNodes(nelem, npoin, velocity, velocity_nodes)
    integer, intent(in) :: nelem, npoin
    real, dimension(1:nelem, 1:2), intent(in) :: velocity
    real, dimension(1:npoin, 0:2), intent(out) :: velocity_nodes
    class(element), pointer :: elem
    type(Bound_edge_t), pointer :: b_edge
    type(volume_rule), pointer :: V_rule
    type(Gauss_rule), pointer :: G_rule
    integer :: i,j,k, k1,k2,k3, ib, j1, j2, l,l1, n1,n2,n3, Qdof
    real ::  val, val1, val2, rho
    real, dimension(:,:), pointer :: xl
    real, dimension(:), allocatable :: xi, ti, xll
    real, dimension(:,:), allocatable :: velP2

    allocate(xi(1:2), ti(1:2), velP2(1:6, 1:2) )


    ! recomputation of the optimal velocity into nodes by least squares
    velocity_nodes(:,:) = 0.
    do i = 1,nelem
       elem => grid%elem(i)
       do j=1, elem%flen
          k = elem%face(idx, j)
          velocity_nodes(k, 0)  = velocity_nodes(k, 0)  + elem%area


          velocity_nodes(k, 1)  = velocity_nodes(k, 1)  + elem%area * velocity(i, 1)
          velocity_nodes(k, 2)  = velocity_nodes(k, 2)  + elem%area * velocity(i, 2)

          ! if(elem%i == 1 .and. j == 1 .and.  state%time%iter_loc <= 3) &
          !      print*, ' given direction of the eikonal velocity !!!!  in euler.f90'
          ! rho = elem%w(0, 1) * sqrt(2.)
          ! val = 2 * exp(-7.5 * (rho / 9.)**2)

          ! velocity_nodes(k, 1)  = velocity_nodes(k, 1)  + elem%area * val
          ! velocity_nodes(k, 2)  = velocity_nodes(k, 2)  + 0.
          !print*,'rho =', rho, val
       enddo
    enddo

    ! division by the weight
    velocity_nodes(1:npoin, 1) = velocity_nodes(1:npoin, 1) / velocity_nodes(1:npoin, 0)
    velocity_nodes(1:npoin, 2) = velocity_nodes(1:npoin, 2) / velocity_nodes(1:npoin, 0)

    !do i=1,npoin
    !!   velocity_nodes(i, 1:2) = grid%x(i, 1:2)/100.
    !   if( grid%x(i,1) < 0.1) &
    !  write(*,*) grid%x(i, 1:2), velocity_nodes(i, 1:2)
    !enddo


    !return


    ! modification for boundary node, velocity has to point inside of the domain
    do ib=1,grid%nbelm
       b_edge => grid%b_edge(ib)

       elem => grid%elem(b_edge%itc)  ! adjacent triangle
       j = b_edge%jtc

       if(elem%iBC(j) == 0) then ! fixed walls
          k1 = b_edge%lbn(1)
          k2 = b_edge%lbn(2)
          k3 = 0

          ! tangential direction from k1
          xi(1:2) = grid%x(k2, 1:2) - grid%x(k1, 1:2)
          val1 = dot_product(velocity_nodes(k1, 1:2),  xi(1:2)  )
          if(val1 > 0) then  ! the correct element
             k3 = k1
             ti(1:2) = xi(1:2)
          endif

          ! tangential direction from k2
          xi(1:2) = grid%x(k1, 1:2) - grid%x(k2, 1:2)
          val2 = dot_product(velocity_nodes(k2, 1:2),  xi(1:2)  )
          if(val2 > 0) then  ! the correct element
             k3 = k2
             ti(1:2) = xi(1:2)
          endif


          if(k3 > 0) then
             ! velocity in the outer normal direction
             val = dot_product(velocity_nodes(k3, 1:2),   elem%n(j, 1:2)  )

             if(val > 0) then ! eikonal velocity points out to the domain

                ! unit tangential direction
                val = dot_product(velocity_nodes(k3, 1:2),  ti(1:2)  )
                if(val < 0.) ti(1:2) = -ti(1:2) ! opposite orientation

                ti(1:2) = ti(1:2)/ sqrt(dot_product(xi(1:2), xi(1:2))) !unit direction

                ! velocity magnitude
                val = sqrt(dot_product(velocity_nodes(k3, 1:2), velocity_nodes(k3, 1:2)))

                ! projection to the tangential direction
                val1 = 1.
                !!!if(abs(elem%xc(2) - 5) < 3.) &
                !!!     val1 = max( 0.5, min( 0.5 + 0.5*(elem%xc(1) - 33), 1.) )

                velocity_nodes(k3, 1:2) = val * ti(1:2) * val1 ! reduction of the velocity


             endif ! if(val < 0) ! eikonal velocity points out to the domain

          !else
             ! stagnation point

          endif ! if k3 > 0
       endif  ! if fixed walls
    enddo

    !stop 'de3de'

    !return

    ! reduction of the velocity in nodes closed to the stagantion point
    ! xi(1:2) = (/ -33., 5. /)
    ! val1 = 2.0
    ! do i=1,grid%npoin
    !    val = sqrt(dot_product(xi(1:2) - grid%x(i, 1:2), xi(1:2) - grid%x(i, 1:2)) )
    !    if(val <= val1) then
    !       velocity_nodes(i, 1:2) = velocity_nodes(i, 1:2)  * sqrt(val /val1)
    !    endif

    ! enddo


    ! ! reduction for the stagnation elements
    ! do ib=1,grid%nbelm
    !    b_edge => grid%b_edge(ib)

    !    elem => grid%elem(b_edge%itc)  ! adjacent triangle
    !    j = b_edge%jtc

    !    if(elem%iBC(j) == 0) then ! fixed walls
    !       k1 = b_edge%lbn(1)
    !       k2 = b_edge%lbn(2)

    !       ! tangential direction from k1
    !       xi(1:2) = grid%x(k2, 1:2) - grid%x(k1, 1:2)
    !       val1 = dot_product(velocity_nodes(k1, 1:2),  xi(1:2)  )

    !       ! tangential direction from k2
    !       xi(1:2) = grid%x(k1, 1:2) - grid%x(k2, 1:2)
    !       val2 = dot_product(velocity_nodes(k2, 1:2),  xi(1:2)  )


    !       if(val1 < 0. .and. val2 < 0) then  ! the stagnation point
    !          velocity_nodes(k1, 1:2) = velocity_nodes(k1, 1:2) / 4.
    !          velocity_nodes(k2, 1:2) = velocity_nodes(k2, 1:2) / 4.
    !          !print*,'EDEREDE',k1,k2, velocity_nodes(k1, 1:2)
    !          !print*,'EDEREDE',k1,k2, velocity_nodes(k1, 1:2)/ 4.
    !       endif
    !    endif
    ! enddo


    ! piecewise linear interpolation into volume integ nodes
    do i = 1,nelem
       elem => grid%elem(i)
       V_rule => state%space%V_rule(elem%Qnum)
       Qdof = elem%Qdof
       elem%xi(0,1:Qdof, 1:2) = 0.

       do j=1,elem%flen
          j1 = mod(j,  elem%flen) + 1
          !j2 = mod(j1, elem%flen) + 1
          k = elem%face(idx, j1)
          elem%xi(0,1:Qdof, 1) = elem%xi(0,1:Qdof, 1) + velocity_nodes(k, 1) * V_rule%lambda(1:Qdof,j)
          elem%xi(0,1:Qdof, 2) = elem%xi(0,1:Qdof, 2) + velocity_nodes(k, 2) * V_rule%lambda(1:Qdof,j)
       enddo

    enddo

    ! piecewise linear interpolation into edge integ nodes on the output
    do ib = 1, grid%nbelm
       b_edge => grid%b_edge(ib)
       !!!if(b_edge%inout == 1) then ! outlet

       elem => grid%elem(b_edge%itc)  ! adjacent triangle
       j = b_edge%jtc
       k1 = b_edge%lbn(1)
       k2 = b_edge%lbn(2)

       G_rule =>  state%space%G_rule(elem%face(fGnum, j))
       ! face integ nodes
       Qdof = G_rule%Qdof
       do k=1, Qdof
          ti = elem%xi(j,k, 1:2)
          val = G_rule%lambda(k)

          xi(1:2) = val * grid%x(k2, 1:2)  + (1.-val) * grid%x(k1, 1:2)
          !!write(*,'(a8,2i5,8es12.4)')'ed43ds',i,k,xi(1:2),elem%xi(j,k, 1:2), elem%xi(j,k, 1:2)-xi(1:2)
          elem%xi(j,k, 1:2) = val * velocity_nodes(k2, 1:2)  + (1.-val) * velocity_nodes(k1, 1:2)

       enddo

    enddo

    !!!!
    return

    ! the following P@ improvement causes some troubles with convergence

    ! P2 improvement for element with the stagnation point
    do ib=1,grid%nbelm
       b_edge => grid%b_edge(ib)

       elem => grid%elem(b_edge%itc)  ! adjacent triangle
       j = b_edge%jtc

       if(elem%iBC(j) == 0 .and. elem%ibcur > 0 ) then ! fixed walls
          k1 = b_edge%lbn(1)
          k2 = b_edge%lbn(2)

          ! tangential direction from k1
          xi(1:2) = grid%x(k2, 1:2) - grid%x(k1, 1:2)
          val1 = dot_product(velocity_nodes(k1, 1:2),  xi(1:2)  )

          ! tangential direction from k2
          xi(1:2) = grid%x(k1, 1:2) - grid%x(k2, 1:2)
          val2 = dot_product(velocity_nodes(k2, 1:2),  xi(1:2)  )


          if(val1 < 0. .and. val2 < 0) then  ! the stagnation point
             j1 = mod(j, 3) + 1
             j2 = mod(j1, 3) + 1
             k3 = elem%face(idx, j2)

             !print*,'EDEREDE',k1,k2, val1, val2, (grid%x(k1, :) + grid%x(k2, :)) / 2
             !write(59,*) (grid%x(k1, :) + grid%x(k2, :)) / 2
             n1 = j2
             n2 = j
             n3 = j1

             velP2(2*n1-1 , 1:2) = velocity_nodes(k1, 1:2)
             velP2(2*n1   , 1:2) = 0.
             !velP2(2*n1   , 1:2) = (velocity_nodes(k1, 1:2) + velocity_nodes(k2, 1:2) ) /2
             velP2(2*n2-1 , 1:2) = velocity_nodes(k2, 1:2)
             velP2(2*n2   , 1:2) = (velocity_nodes(k2, 1:2) + velocity_nodes(k3, 1:2) ) /2
             velP2(2*n3-1 , 1:2) = velocity_nodes(k3, 1:2)
             velP2(2*n3   , 1:2) = (velocity_nodes(k3, 1:2) + velocity_nodes(k1, 1:2) ) /2

             ! volume integ nodes
             V_rule => state%space%V_rule(elem%Qnum)
             Qdof = elem%Qdof
             elem%xi(0,1:Qdof, 1:2) = 0.

             allocate( xll(1:Qdof) )

             xl => V_rule%lambda(1:Qdof,1:3)

             do l =1, 3
                l1 = mod(l, 3) + 1
                !vertex
                xll(1:Qdof) = (2*xl(1:Qdof, l) - 1) * xl(1:Qdof, l)

                elem%xi(0,1:Qdof,1) = elem%xi(0,1:Qdof,1) + xll(1:Qdof) * velP2(2*l - 1, 1)
                elem%xi(0,1:Qdof,2) = elem%xi(0,1:Qdof,2) + xll(1:Qdof) * velP2(2*l - 1, 2)

                ! edge middle
                xll(1:Qdof) = 4 * xl(1:Qdof,l) * xl(1:Qdof,l1)

                elem%xi(0,1:Qdof,1) = elem%xi(0,1:Qdof,1) + xll(1:Qdof) * velP2(2*l , 1)
                elem%xi(0,1:Qdof,2) = elem%xi(0,1:Qdof,2) + xll(1:Qdof) * velP2(2*l , 2)

             enddo


             ! Boundary integ nodes
             G_rule =>  state%space%G_rule(elem%face(fGnum, j))
             ! face integ nodes
             Qdof = G_rule%Qdof
             do k=1, Qdof
                val = G_rule%lambda(k)
                xi(1:2) = val * grid%x(k2, 1:2)  + (1.-val) * grid%x(k1, 1:2)

                elem%xi(j,k, 1:2) = 2 * val * (val -0.5) * velocity_nodes(k2, 1:2)  &
                     + 2. * (1.-val) * (0.5 - val) * velocity_nodes(k1, 1:2)

             enddo
             deallocate(xll)

          endif  ! if stagnation point
       endif  ! fixed walls
    enddo


    deallocate(xi, ti, velP2)

  end subroutine EikonalVelocityIntegNodes

  !> plot the eikonal velocity arrays
  subroutine  PlotEikonalVelocityVectors(ifile, ifile1)
    integer, intent(in) :: ifile, ifile1
    class(element), pointer:: elem
    real, dimension(:,:), allocatable :: xi, Fx
    integer :: i, j, k, k1, k2, Qdof, Qnum
    type(Gauss_rule), pointer :: G_rule


    write(ifile, *) '##   PlotEikonalVelocityVectors(ifile) in euler.f90'

    do i=1,grid%nelem
       elem => grid%elem(i)

       ! plot of the vectors
       !write(ifile,*) elem%xc(:)! , elem%xi(0,1, 1:nbDim)
       !write(ifile,*) elem%xc(:)+elem%xi(0,1, 1:nbDim)! , elem%xi(0,1, 1:nbDim)
       !write(ifile,*)

       Qdof = elem%Qdof
       Qnum = elem%Qnum

       allocate( Fx(1:Qdof,1:nbDim) )

       call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim), &
            Fx(1:Qdof,1:nbDim) )

       do j=1,Qdof
          !if(abs(elem%xi(0,j, 2) - 5) < 5E-2) &
          !     write(ifile+1000,*) Fx(j,1:2), elem%xi(0,j, 1:nbDim), &
          !     sqrt(dot_product(elem%xi(0,j,:), elem%xi(0,j,:)))
          write(ifile,*) Fx(j, 1:2)
          write(ifile,*) Fx(j, 1:2) +elem%xi(0,j, 1:nbDim)
          write(ifile,*) Fx(j, 1:2) +elem%xi(0,j, 1:nbDim)*0.95
          write(ifile,*) '  '
       enddo
       deallocate(Fx)

       do k=1,elem%flen
          if(elem%face(neigh,k) < 0. .and.  elem%iBC(k) == 0) then ! fixed walls
             k1 = mod(k, 3) + 1
             k2 = mod(k1, 3) + 1

             G_rule =>  state%space%G_rule(elem%face(fGnum, k))
             Qdof = G_rule%Qdof
             allocate(xi (1:Qdof, 1:3), Fx (1:Qdof, 1:2) )
             xi = 0.
             xi(1:Qdof, k) = G_rule%lambda(1:Qdof)
             xi(1:Qdof, k2 ) = 1. - G_rule%lambda(1:Qdof)

             call ComputeF(elem, Qdof, xi(1:Qdof, 1:2),  Fx(1:Qdof,1:nbDim) )

             do j=1,Qdof
                write(ifile1,*) Fx(j, 1:2)
                write(ifile1,*) Fx(j, 1:2) +elem%xi(k,j, 1:nbDim)
                write(ifile1,*) Fx(j, 1:2) +elem%xi(k,j, 1:nbDim)*0.95
                write(ifile1,*) '  '
             enddo

             deallocate(xi, Fx)
          endif
       enddo


      !if(elem%xc(1) > 32 .and. elem%xc(1) < 34 .and. elem%xc(2) > 4 .and. elem%xc(2) < 6) then
!       write(*,'(a8,40i6)') 'ede43',deg,dof,QQdof, V_rule%Qdof, Qdof, size(V_rule%phi, 1), size(V_rule%phi, 2)
!         write(*,'(a8,40es12.4)') 'EDED43', elem%xi(0, 1:Qdof, 2)
!       write(*,'(a8,40es12.4)') 'ede43',sqrt( elem%xi(0, 1:QQdof, 1)**2 + elem%xi(0, 1:QQdof, 2)**2 )
!       write(*,'(a8,40es12.4)') 'ede43', wwi(1, 1:dof)
!       write(*,'(a8,40es12.4)') 'ede43', wi(1, 1:dof)
!       write(*,'(a8,40es12.4)') 'ede43',qExact(1:Qdof, 1)
!       stop
!    endif

       !if(i ==5) stop '487y7y790ik'
    enddo
  end subroutine PlotEikonalVelocityVectors


  !> solves the Eikonal equation
  subroutine EikonalVelocity(npoin, nelem, max_nbp,x,lnd, loc_ibp, density, velocity, use_default_initc)
     use BRAlgorithm
     implicit none
     integer, intent(in) :: npoin ! number of the grid nodes
     integer, intent(in) :: nelem ! number of the grid triangles
     integer, intent(in) :: max_nbp ! maximal number of neighbouring nodes
     real, dimension(1:npoin, 1:2), intent(in) :: x  ! coordinates of the nodes
     integer, dimension(1:nelem, 1:3), intent(in) :: lnd ! indexes of nodes forming each triangle

     integer, dimension(1:npoin, 0:max_nbp), intent(in) ::loc_ibp !loc_ibp(i,0) = number of neigh nodes
                                                                  !loc_ibp(i,1:) = indexes of ^^^^^^^

     real, dimension(1:npoin), intent(in) :: density  ! density at vertexes - INPUT
     real, dimension(1:nelem, 1:2), intent(out) :: velocity  ! velocity on elements - OUTPUT
     logical ,intent(in):: use_default_initc ! set which type of initial condition for BR algorithm is used
     integer :: i


     ! write(200, *) npoin, nelem,max_nbp
     ! do i=1,npoin
     !    write(200, *) x(i, 1:2)
     ! enddo
     ! do i=1,nelem
     !    write(200, *) lnd(i, 1:3)
     ! enddo
     ! do i=1, npoin
     !    write(200, *) loc_ibp(i, 0: max_nbp  )
     ! enddo
     ! do i=1, npoin
     !    write(200, *) density(i)
     ! enddo


     !print*,'initialize BR-algorithm, use_default_initc = ',use_default_initc
     call InitializeBRAlgorithm(npoin,density,use_default_initc)

     !print*,'solve eikonal equation'
     call BRSolve(npoin,max_nbp,x,loc_ibp)

     !compute velocity from potential phi
     call ComputeVelocity(npoin, nelem,  max_nbp, x(1:npoin, 1:2), lnd(1:nelem, 1:3), &
          density(1:npoin),velocity(1:nelem,1:2))


     ! EXPORT DATA -comment !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !for testing purposes, export phi (phi.dta), x-velocity (vx.dta), y-velocity (vy.dta) for gnuplot    pm3d
     !call Export(npoin, nelem,  max_nbp, x(1:npoin, 1:2), lnd(1:nelem, 1:3),velocity(1:nelem,1:2))


     !stop 'ede7d38e3dj3o'



   end subroutine EikonalVelocity


 !> try to avoin nonphysical behaiviour of the velocity for small density
 subroutine Pedestrian_velocity_limiting2( )
   use pedes_averaging
   class(element), pointer :: elem, elem1 ! elem = element
   real, dimension(:,:), pointer:: phi ! local store arrays
   real, dimension(:), pointer:: weights ! local store arrays
   real, dimension(:,:), allocatable :: wi ! w in integ. nodes
   real, dimension(:,:), allocatable :: wS, w ! new wS
   real, dimension(:), allocatable :: w_min, velocity
   integer :: i, j, k, dof, Qdof
   real :: density, density_tot, density_tot1, vel(1:2)

   do i=1,grid%nelem
      elem => grid%elem(i)

      dof = elem%dof
      Qdof = elem%Qdof

      phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)
      weights => state%space%V_rule(elem%Qnum)%weights(1:Qdof)

      allocate(wi(1:3, 1:Qdof) ) ! w in integ. nodes
      allocate(w_min(1:3) ) ! minimal value in integ nodes

      call Eval_min_w_Elem(elem, w_min(1:ndim) )

      do k=1,ndim
         wi(k, 1:Qdof) = matmul(elem%w(0,(k-1)*dof+ 1: k*dof), phi(1:dof, 1:Qdof) )
      enddo

      density = dot_product(wi(1, 1:Qdof), weights(1:Qdof) )

      !if(elem%xc(1) < 0.3 .and. state%time%iter_loc >= 15) &

      !!if(density < 0.) then
         !write(*,'(a8,i5,50es12.4)') 'limiting2', elem%i, density, w_min(1), elem%xc(:)
      !!elem%w(0,  1:ndim*dof)  = 0.

      ! small density, velocity is avereged and taken as the eikonal velocity
      !!else

      if(density < state%model%Pr) then
      !if(w_min(1) < state%model%Pr) then
         !allocate(velocity (1:Qdof) )
         !velocity(1:Qdof) = sqrt(  wi(2, 1:Qdof)**2 +  wi(3, 1:Qdof)**2)
         !write(*,'(a4,50es12.4)') 'velo',velocity(1:Qdof)

         ! eikonal velocity, DOES NOT work for P_0 (FVM), although it is similar to the following one
         !wi(2, 1:Qdof) =  abs(wi(1, 1:Qdof)) * elem%xi(0,1:Qdof, 1)
         !wi(3, 1:Qdof) =  abs(wi(1, 1:Qdof)) * elem%xi(0,1:Qdof, 2)

         ! eikonal velocity, works for P_0 (FVM)
         !vel(1) = dot_product(elem%xi(0,1:Qdof, 1), weights(1:Qdof) )
         !vel(2) = dot_product(elem%xi(0,1:Qdof, 2), weights(1:Qdof) )

         !wi(2, 1:Qdof) =  density * vel(1)
         !wi(3, 1:Qdof) =  density * vel(2)

         !wi(2, 1:Qdof) =  density * elem%xi(0,1:Qdof, 1) * (density / state%model%Pr)**2
         !wi(3, 1:Qdof) =  density * elem%xi(0,1:Qdof, 2) * (density / state%model%Pr)**2

         !do j=1,Qdof
         !   if(velocity(j) > 0.) then
         !      wi(2, j) =  wi(1, j) * wi(2, j) /velocity(j) * 2
         !      wi(3, j) =  wi(1, j) * wi(3, j) /velocity(j) * 2
         !   else
         !      wi(2:3, j) = 0.
         !   endif
         !
         !enddo

         !write(*,'(a8,i5, 30es16.8)')  'lim 2:',elem%i, elem%w(0, 1:dof )
         !write(*,'(a8,i5, 30es16.8)')  'lim 2:',dof, elem%w(0, dof+1:2*dof )
         !write(*,'(a8,i5, 30es16.8)')  'lim 2:',-1, elem%w(0, 2*dof+1:3*dof )
         !print*

         ! for P0 DGM
         elem%w(0,    dof + 1 : 3*dof) = 0.

         elem%w(0,          1 ) =  density / phi(1,1)
         elem%w(0,    dof + 1 ) = dot_product(wi(2, 1:Qdof), weights(1:Qdof) ) / phi(1,1)
         elem%w(0, 2* dof + 1 ) = dot_product(wi(3, 1:Qdof), weights(1:Qdof) ) / phi(1,1)

         !write(*,'(a8,i5, 30es16.8)')  'lim 2:',elem%i, elem%w(0, 1:dof )
         !write(*,'(a8,i5, 30es16.8)')  'lim 2:',dof, elem%w(0, dof+1:2*dof )
         !write(*,'(a8,i5, 30es16.8)')  'lim 2:',-1, elem%w(0, 2*dof+1:3*dof )
         !print*,'------------------------------------'

         ! allocate(wS(2:3, 1:dof), w(2:3, 1:dof) ) ! w in integ. nodes

         ! call IntegrateVectorB(elem, dof, wi(2, 1:Qdof), wS(2, 1:dof) )
         ! call IntegrateVectorB(elem, dof, wi(3, 1:Qdof), wS(3, 1:dof) )

         ! w(2,1:dof) = matmul(elem%MassInv%Mb(1:dof, 1:dof), wS(2, 1:dof) )
         ! w(3,1:dof) = matmul(elem%MassInv%Mb(1:dof, 1:dof), wS(3, 1:dof) )

         ! !write(*,'(a4,50es12.4)') 'limiting2', density
         ! !write(*,'(a4,50es12.4)') 'ws 2', elem%w(0,dof + 1: 2*dof)
         ! !write(*,'(a4,50es12.4)') 'ws 2', w(2,1:dof)
         ! !print*
         ! !write(*,'(a4,50es12.4)') 'ws 3', elem%w(0,2*dof + 1: 3*dof)
         ! !write(*,'(a4,50es12.4)') 'ws 3', w(3,1:dof)
         ! !print*,'--------------------------'

         ! !if(elem%i == 2214 .or. elem%i == 2201) then
         ! !   write(*,'(a4,i5,4es12.4,a4, 40es12.4)') 'w23', elem%i, &
         ! !        w(2:3,1:dof), elem%w(0,2)/elem%w(0,1),  elem%w(0,3)/elem%w(0,1), &
         ! !        '|',elem%w(0,1:3)
         ! !   !'|',wi(1, 1:Qdof)
         ! !endif


         ! !elem%w(0,    dof + 1 : 2*dof) = w(2,1:dof)
         ! !elem%w(0, 2* dof + 1 : 3*dof) = w(3,1:dof)

         ! !elem%w(0,    dof + 2 : 2*dof) = 0.
         ! !elem%w(0, 2* dof + 2 : 3*dof) = 0.

         ! !!!!elem%w(0,  1:ndim*dof)  = 0.
         !  deallocate(wS , w, velocity)

      endif
      deallocate(wi, w_min)
   enddo


 end subroutine Pedestrian_velocity_limiting2


 !> try to avoin nonphysical behaiviour of the velocity for small density
 subroutine Pedestrian_velocity_limiting( )
   use pedes_averaging
   class(element), pointer :: elem, elem1 ! elem = element
   type(volume_rule), pointer :: V_rule
   real, dimension(:,:,:), allocatable :: velocity, Dphi
   integer, dimension(:), allocatable :: ic
   real, dimension(:,:), allocatable :: w, xc, DwR, ww
   real, dimension(:,:,:), allocatable :: wi
   real, dimension(:,:), pointer:: phi ! test functions
   real :: diff, norm, normR, det, valx, valy
   real, dimension(:,:), allocatable :: v_proj, v_proj_nodes
   integer :: i, j, k, l, nei, dof, Qdof


   allocate(w(1:ndim, 0:2 ) )
   !allocate(velocity(1:grid%nelem, 1:ndim, 0:3), xc(1:3, 1:2), ic(1:3), DwR(1:ndim, 1:2) )

   ! w(1:ndim, 0:2)        second index =0, w, =1, D_x w, =2 D_y
   ! velocity(:,1:ndim, 0:2)  3rd index =0, w, =1, D_x w, =2 D_y

   !allocate(v_proj(0:2, 1:2) )  ! reconstructed projection of the velocity, P_0 and P_1
   limit_w_bar(:) = 0

   if(state%time%iter_loc == 1) print*,' Pedestrian_velocity_limiting( ) is NOT USED'
   return

   ! evaluation of the mean values of w and Dw on triangles
   do i=1,grid%nelem
      elem => grid%elem(i)

      call ElementJumpVelocityIndicator(elem)
      elem%rezid = elem%rezid / elem%diam

      if(elem%rezid > 1E-00) then
         !limit_w_bar(i) = 1
         call Eval_aver_w_Elem(elem, w(1:ndim, 0) )
         w_bar(i, 1:ndim) = w(1:ndim, 0)
         !write(*, *) elem%xc, elem%rezid
         !write(100+state%time%iter*10, *) elem%xc, elem%rezid, w_bar(i, 2)/ w_bar(i,1), w_bar(i,:), elem%i
      endif



      ! call Eval_aver_w_Elem(elem, w(1:ndim, 0) )
      ! call Eval_aver_Dw_Elem(elem, w(1:ndim, 1:2) )

      ! if(w(1,0) > 0) then
      !    !mean value of the density and velocity
      !    velocity(i, 1, 0) = w(1,0)
      !    velocity(i, 2, 0) = w(2,0) / w(1,0)
      !    velocity(i, 3, 0) = w(3,0) / w(1,0)

      !    !mean value of the derivative of the density and velocity with respect to x_1
      !    velocity(i, 1, 1) =  w(1,1)
      !    velocity(i, 2, 1) = (w(2,1) * w(1,0) - w(2,0) * w(1,1) ) / w(1,0) /w(1,0)
      !    velocity(i, 3, 1) = (w(3,1) * w(1,0) - w(3,0) * w(1,1) ) / w(1,0) /w(1,0)

      !    !mean value of the derivative of the density and velocity with respect to x_2
      !    velocity(i, 1, 2) =  w(1,2)
      !    velocity(i, 2, 2) = (w(2,2) * w(1,0) - w(2,0) * w(1,2) ) / w(1,0) /w(1,0)
      !    velocity(i, 3, 2) = (w(3,2) * w(1,0) - w(3,0) * w(1,2) ) / w(1,0) /w(1,0)


      ! else
      !    velocity(i,:,:)= 0.
      ! endif

      ! Qdof = elem%Qdof
      ! allocate(wi(1:Qdof, 1:ndim, 0:2))
      ! call Eval_w_Elem(elem, wi(1:Qdof, 1:ndim, 0))
      ! call Eval_Dw_Elem(elem, wi(1:Qdof, 1:ndim, 1:2) )

      ! ! velocity in integ nodes
      ! do l=1,Qdof
      !    if(wi(l,1,0) > 0) then
      !       wi(l, 2:3, 0) = wi(l, 2:3, 0) / wi(l, 1, 0)
      !    else
      !       wi(l, 2:3, 0) = 0.
      !    endif
      !    !write(70,*) elem%xc(:), wi(l, 2:3, 0)
      ! enddo

      ! ! square of the difference of the velocity and the mean value
      ! wi(1:Qdof, 2, 0) = (wi(1:Qdof, 2, 0) -  velocity(i, 2, 0))**2
      ! wi(1:Qdof, 3, 0) = (wi(1:Qdof, 3, 0) -  velocity(i, 3, 0))**2


      ! ! integration of v - \bar{v} over K
      ! call IntegrateFunction(elem, wi(1:Qdof, 2, 0) , valx)
      ! call IntegrateFunction(elem, wi(1:Qdof, 3, 0) , valy)

      ! ! detection for the limiting
      ! velocity(i, 1, 3) = sqrt(valx + valy) / (elem%area / elem%diam)

      ! !if(elem%xc(2) > 4 .and. elem%xc(2) < 6) then
      !    ! write(71,*) elem%xc(:), 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
      !    ! write(71,*) elem%xc(:), velocity(i, 1:3, 0), velocity(i, 1:3, 1), velocity(i, 1:3, 2)
      !    ! write(71,'(x)')
      !    ! write(71,'(x)')

      !    ! write(72,*) elem%xc(:), 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
      !    ! write(72,*) elem%xc(:), valx, valy, sqrt(valx + valy), sqrt(valx + valy) / (elem%area / elem%diam)
      !    ! write(72,'(x)')
      !    ! write(72,'(x)')
      ! !endif

      ! deallocate( wi)
   enddo

   !stop 'deu73ed29p3w'

   ! ! evaluation of the gradients from neighbouring mean values
   ! do i=1,grid%nelem
   !    elem => grid%elem(i)
   !    ! nei = 0
   !    ! do j=1, elem%flen
   !    !    k = elem%face(neigh, j)
   !    !    if(k > 0) then
   !    !       elem1 => grid%elem(k)
   !    !       nei = nei + 1
   !    !       ic(j) = k
   !    !       xc(j, 1:2) = elem1%xc(1:2)
   !    !    endif
   !    ! enddo

   !    !if( nei == 3) then
   !    !   call GradientReconstruct(xc(1:3, 1:2), velocity(ic(1),1:ndim, 0:2),  &
   !    !        velocity(ic(2),1:ndim, 0:2), velocity(ic(3),1:ndim, 0:2), DWR(1:ndim, 1:2) )

   !    !   ! difference between mean values of the velocity gradient and its reconstruction
   !    !   diff = sqrt( dot_product(velocity(i, 2, 1:2) - DwR(2, 1:2), velocity(i, 2, 1:2) - DwR(2, 1:2))&
   !    !        +  dot_product(velocity(i, 3, 1:2) - DwR(3, 1:2), velocity(i, 3, 1:2) - DwR(3, 1:2)) )

   !    !  ! the sizes of the velocity gradient
   !    !   norm = sqrt( dot_product(velocity(i, 2, 1:2), velocity(i, 2, 1:2) )&
   !    !        +  dot_product(velocity(i, 3, 1:2), velocity(i, 3, 1:2) ) )

   !    !   normR = sqrt( dot_product( DwR(2, 1:2),  DwR(2, 1:2))&
   !    !        +  dot_product( DwR(3, 1:2),  DwR(3, 1:2)) )


   !    !   if(diff > 0.05 .and. norm > 0.0 * normR .and. norm > 0) then

   !    if(velocity(i, 1, 3) > 0.1) then
   !          ! limiting
   !          !if(diff / norm > 0.5) then

   !       !write(75,*) elem%xc(:), diff/norm, norm, normR, diff, norm - normR, norm/normR, &
   !       !        diff / abs(norm - normR),           DwR(2:3, 1), DwR(2:3, 2)


   !       V_rule => state%space%V_rule(elem%Qnum)
   !       Qdof = V_rule%Qdof

   !       dof = elem%dof

   !       allocate(v_proj_nodes(0:2, 1:Qdof)) ! reconstructed projection of the velocity in integ nodes

   !       phi => V_rule%phi(1:dof,1:Qdof)
   !       allocate(Dphi(1:3, 1:nbDim, 1:Qdof))  ! first three functions are enough
   !       call  Eval_Dphi(elem, 3,  Dphi(1:3, 1:nbDim, 1:Qdof) )

   !       ! P_0 reconstruction
   !       v_proj(0, 1:2) = velocity(i, 2:3, 0) /  phi(1,1)

   !       ! P_1 reconstruction
   !       ! Dphi(2:3, 1:2, 1:Qdof) = constants
   !       det = Dphi(2, 1, 1) * Dphi(3, 2, 1) - Dphi(2, 2, 1) * Dphi(3, 1, 1)

   !       v_proj(1, 1:2) = velocity(i, 2:3, 1) *  Dphi(3, 2, 1) - velocity(i, 2:3, 2) *  Dphi(3, 1, 1)
   !       v_proj(2, 1:2) = velocity(i, 2:3, 2) *  Dphi(2, 1, 1) - velocity(i, 2:3, 1) *  Dphi(2, 2, 1)

   !       v_proj(1, 1:2) =  v_proj(1, 1:2) / det
   !       v_proj(2, 1:2) =  v_proj(2, 1:2) / det

   !       ! P_0 velocity in integ node
   !       v_proj_nodes(1,1:Qdof) = v_proj(0, 1) * phi(1,1:Qdof)
   !       v_proj_nodes(2,1:Qdof) = v_proj(0, 2) * phi(1,1:Qdof)

   !       ! P_1 velocity in integ node
   !       !v_proj_nodes(1,1:Qdof) = v_proj_nodes(1,1:Qdof) &
   !       !     + v_proj(1, 1) * phi(2,1:Qdof) + v_proj(2, 1) * phi(3,1:Qdof)

   !       !v_proj_nodes(2,1:Qdof) = v_proj_nodes(2,1:Qdof) &
   !       !     + v_proj(1, 2) * phi(2,1:Qdof) + v_proj(2, 2) * phi(3,1:Qdof)

   !       ! density in integ nodes
   !       v_proj_nodes(0,1:Qdof) = matmul(elem%w(0, 1:dof), phi(1:dof,1:Qdof) )

   !       ! P_0 density
   !       v_proj_nodes(0,1:Qdof) = velocity(i, 1, 0)

   !       !write(*,'(a8, 30es12.4)') 'rho',  v_proj_nodes(0,1:Qdof)
   !       !write(*,'(a8, 30es12.4)') 'v_1',  v_proj_nodes(1,1:Qdof)
   !       !write(*,'(a8, 30es12.4)') 'v_2',  v_proj_nodes(2,1:Qdof)


   !       ! momentum in integ nodes = velocity * density
   !       v_proj_nodes(1,1:Qdof) = v_proj_nodes(1,1:Qdof) * v_proj_nodes(0,1:Qdof)
   !       v_proj_nodes(2,1:Qdof) = v_proj_nodes(2,1:Qdof) * v_proj_nodes(0,1:Qdof)

   !       !write(*,'(a8, 30es12.4)') 'rv_1',  v_proj_nodes(1,1:Qdof)
   !       !write(*,'(a8, 30es12.4)') 'rv_2',  v_proj_nodes(2,1:Qdof)

   !       allocate(ww(0:2, 1:dof)) ! RHS of the projection =\int_K v_proj * \phi_i \dx
   !       call IntegrateVectorB(elem, dof, v_proj_nodes(0,1:Qdof), ww(0, 1:dof) )
   !       call IntegrateVectorB(elem, dof, v_proj_nodes(1,1:Qdof), ww(1, 1:dof) )
   !       call IntegrateVectorB(elem, dof, v_proj_nodes(2,1:Qdof), ww(2, 1:dof) )

   !       ! solution in basis coeffs, projection
   !       elem%w(0,         1:   dof) = matmul(elem%MassInv%Mb(1:dof, 1:dof), ww(0, 1:dof) )
   !       elem%w(0,   dof + 1: 2*dof) = matmul(elem%MassInv%Mb(1:dof, 1:dof), ww(1, 1:dof) )
   !       elem%w(0, 2*dof + 1: 3*dof) = matmul(elem%MassInv%Mb(1:dof, 1:dof), ww(2, 1:dof) )


   !       !call PlotElem_D_Function3D(103, elem,  dof, elem%w(0, 1:  dof) )
   !       !call PlotElem_D_Function3D(104, elem,  dof, elem%w(0, dof + 1:  2*dof) )
   !       !call PlotElem_D_Function3D(105, elem,  dof, elem%w(0, 2*dof + 1:  3*dof) )

   !       !print*,'EDE#E', state%time%iter, state%nlSolver%iter, i
   !       !endif


   !       deallocate(v_proj_nodes, Dphi, ww)
   !       !stop "d4e43fd48eiu32"

   !    endif


   !       !write(72,*) elem%xc(:), velocity(i, 1:3, 1), velocity(i, 1:3, 2)
   !       !write(73,*) elem%xc(:), DwR(1:3, 1), DwR(1:3, 2)

   !    !endif

   ! enddo


   !stop ' Pedestrian_velocity_limiting'

   deallocate(w)
   !deallocate(velocity, xc, ic, DWR, v_proj)

 end subroutine Pedestrian_velocity_limiting




 !> try to avoin nonphysical behaiviour of the velocity for small density
 subroutine Pedestrian_velocity_limiting_STDG( )
   class(element), pointer :: elem, elem1 ! elem = element
   class(Time_rule), pointer :: T_rule
   type(volume_rule), pointer :: V_rule
   real, dimension(:,:,:), allocatable :: velocity, Dphi
   integer, dimension(:), allocatable :: ic
   real, dimension(:,:), allocatable :: w, xc, DwR, ww
   real, dimension(:,:,:), allocatable :: Dwi
   real, dimension(:,:), pointer:: phi ! test functions
   real :: diff, norm, normR, det, valx, valy
   real, dimension(:,:), allocatable :: v_proj, v_proj_nodes
   integer :: i, j, k, l, nei, alpha, dof, Qdof, Qdeg, Tdof


   !allocate(velocity(1:grid%nelem, 1:ndim, 0:3), w(1:ndim, 0:2 ), xc(1:3, 1:2), ic(1:3), DwR(1:ndim, 1:2) )
   ! w(1:ndim, 0:2)        second index =0, w, =1, D_x w, =2 D_y
   ! velocity(:,1:ndim, 0:2)  3rd index =0, w, =1, D_x w, =2 D_y

   !allocate(v_proj(0:2, 1:2) )  ! reconstructed projection of the velocity, P_0 and P_1

   Qdeg = maxval(grid%elem(:)%TQnum )  !!state%time%Qnum
   T_rule => state%time%T_rule(Qdeg)
   Tdof = T_rule%Qdeg
   !Tdof = elem%Tdof

  ! evaluation of the mean values of w and Dw on triangles
   do i=1,grid%nelem
      elem => grid%elem(i)

      call ElementJumpVelocityIndicator(elem)
      elem%rezid = elem%rezid / elem%diam

      if(elem%rezid > 1E-2) then
         write(*, *) elem%xc, elem%rezid
         !write(100+state%time%iter*10, *) elem%xc, elem%rezid
      endif


      ! dof = elem%dof
      ! theta = 1.

      ! Qdeg = elem%TQnum   !!state%time%Qnum
      ! Tdof = elem%Tdof
      ! T_rule => state%time%T_rule(Qdeg)

      ! do alpha = 1, Qdeg
      !    call Transfer_wST_to_w_Elem(elem , alpha, Qdeg)
      ! enddo

      ! call Eval_aver_w_Elem(elem, w(1:ndim, 0) )
      ! call Eval_aver_Dw_Elem(elem, w(1:ndim, 1:2) )

      ! if(w(1,0) > 0) then
      !    !mean value of the density and velocity
      !    velocity(i, 1, 0) = w(1,0)
      !    velocity(i, 2, 0) = w(2,0) / w(1,0)
      !    velocity(i, 3, 0) = w(3,0) / w(1,0)

      !    !mean value of the derivative of the density and velocity with respect to x_1
      !    velocity(i, 1, 1) =  w(1,1)
      !    velocity(i, 2, 1) = (w(2,1) * w(1,0) - w(2,0) * w(1,1) ) / w(1,0) /w(1,0)
      !    velocity(i, 3, 1) = (w(3,1) * w(1,0) - w(3,0) * w(1,1) ) / w(1,0) /w(1,0)

      !    !mean value of the derivative of the density and velocity with respect to x_2
      !    velocity(i, 1, 2) =  w(1,2)
      !    velocity(i, 2, 2) = (w(2,2) * w(1,0) - w(2,0) * w(1,2) ) / w(1,0) /w(1,0)
      !    velocity(i, 3, 2) = (w(3,2) * w(1,0) - w(3,0) * w(1,2) ) / w(1,0) /w(1,0)


      ! else
      !    velocity(i,:,:)= 0.
      ! endif

      ! Qdof = elem%Qdof
      ! allocate(wi(1:Qdof, 1:ndim, 0:2))
      ! call Eval_w_Elem(elem, wi(1:Qdof, 1:ndim, 0))
      ! call Eval_Dw_Elem(elem, wi(1:Qdof, 1:ndim, 1:2) )

      ! ! velocity in integ nodes
      ! do l=1,Qdof
      !    if(wi(l,1,0) > 0) then
      !       wi(l, 2:3, 0) = wi(l, 2:3, 0) / wi(l, 1, 0)
      !    else
      !       wi(l, 2:3, 0) = 0.
      !    endif
      !    !write(70,*) elem%xc(:), wi(l, 2:3, 0)
      ! enddo

      ! ! square of the difference of the velocity and the mean value
      ! wi(1:Qdof, 2, 0) = (wi(1:Qdof, 2, 0) -  velocity(i, 2, 0))**2
      ! wi(1:Qdof, 3, 0) = (wi(1:Qdof, 3, 0) -  velocity(i, 3, 0))**2


      ! ! integration of v - \bar{v} over K
      ! call IntegrateFunction(elem, wi(1:Qdof, 2, 0) , valx)
      ! call IntegrateFunction(elem, wi(1:Qdof, 3, 0) , valy)

      ! ! detection for the limiting
      ! velocity(i, 1, 3) = sqrt(valx + valy) / (elem%area / elem%diam)

      ! !if(elem%xc(2) > 4 .and. elem%xc(2) < 6) then
      !    ! write(71,*) elem%xc(:), 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
      !    ! write(71,*) elem%xc(:), velocity(i, 1:3, 0), velocity(i, 1:3, 1), velocity(i, 1:3, 2)
      !    ! write(71,'(x)')
      !    ! write(71,'(x)')

      !    ! write(72,*) elem%xc(:), 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
      !    ! write(72,*) elem%xc(:), valx, valy, sqrt(valx + valy), sqrt(valx + valy) / (elem%area / elem%diam)
      !    ! write(72,'(x)')
      !    ! write(72,'(x)')
      ! !endif

      ! deallocate( wi)
   enddo
 end subroutine Pedestrian_velocity_limiting_STDG


 !> reconstruction of the gradient from point values on triangle
 subroutine GradientReconstruct(xc, Dw1, Dw2, Dw3, DwR)
   real, dimension(1:3, 1:2), intent(in) :: xc ! nodes coordinates
   real, dimension(1:ndim, 0:2), intent(in)  :: Dw1 ! gradient of w at 1st node
   real, dimension(1:ndim, 0:2), intent(in)  :: Dw2 ! gradient of w at 2nd node
   real, dimension(1:ndim, 0:2), intent(in)  :: Dw3 ! gradient of w at 3rd node
   real, dimension(1:ndim, 1:2), intent(inout) :: DwR ! reconstructed gradient
   real :: det

   det = xc(1,1)*( xc(2,2) - xc(3,2)) + xc(2,1)*( xc(3,2) - xc(1,2)) + xc(3,1)*( xc(1,2) - xc(2,2))

   DwR(1:ndim, 1) = Dw1(1:ndim,0) * (xc(2,2) - xc(3,2) ) &
        + Dw2(1:ndim,0) * (xc(3,2) - xc(1,2)) + Dw3(1:ndim,0) * (xc(1,2) - xc(2,2) )

   DwR(1:ndim, 2) = Dw1(1:ndim,0) * (xc(3,1) - xc(2,1) ) &
        + Dw2(1:ndim,0) * (xc(1,1) - xc(3,1)) + Dw3(1:ndim,0) * (xc(2,1) - xc(1,1) )

   DWR(1:ndim, 1:2) = DWR(1:ndim, 1:2) / det

 end subroutine GradientReconstruct



 !> evaluation of the matrix and vector blocks for STDGM method
 !> if implicitly == False -> vector F(u) saved in elem%vec(rhs,:)
 !> if implicitly == True  -> sets the matrix C(*,*) elem%vec(rhs = q(u),
 subroutine ComputeCplus( deg_plus )
   logical, intent(in) :: deg_plus
   class(element), pointer :: elem
   class(Time_rule), pointer :: T_rule
   !real, dimension(:,:), pointer :: phi
   integer :: i, j, k, l, r, kk
   integer :: alpha, Qdeg, dof, Tdof, s_dim, f_dof, f_Tdof, f_dim
   integer :: wTdof, wdof  !NEW for adaptation - smaller than dof/Tdof if elem%deg_plus=TRUE
   real :: cTime, tt, tt1
   real :: val, val_vec
   real :: local_eta
   logical ::iprint
   integer :: dofA

   !   print*, 'ComputeSTDGM_Terms'
   call cpu_time(tt1)

   cTime = state%time%ctime

   local_eta = 1. / state%time%tau(1)


   !   !F@R control phi(Tdof + 1) ma koreny v int uzlech


   !NEW
   Qdeg = state%time%Qnum

   T_rule => state%time%T_rule(Qdeg)

   associate ( time => state%time)
   select type ( time )
      class is ( TimeTDG_t )

         ! zero
         call ClearMatrixBlocks_STDGM()

         call ClearVectorBlocks_STDGM()

         ! new allocation


         do alpha = 1, Qdeg ! temporarily max_Tdof =  max time quadrature nodes

            do i = 1, grid%nelem
               elem => grid%elem(i)

               !	G_rule => state%space%G_rule(elem%TQnum)
               if (Qdeg /= elem%TQnum) then
                  !F@R Verify if it is OK, some nodes could be in wrong position
               endif
               !call Transfer_wST_to_w_Elem(elem , alpha, elem%TQnum)
               ! save the wST space-time solution in quadrature index alpha to w
               call Transfer_wST_to_w_Elem(elem , alpha, Qdeg)
            enddo

            !we have to run ComputeTerms() in the time-quadrature nodes
            !cTime = state%time%ctime
            state%time%ctime = state%time%ttime + state%time%tau(1) * T_rule%lambda(alpha)
            !   print*, 'ctime:' , state%time%ctime, cTime
            !   print*, '--------------','ctime:', state%time%ctime, state%nlSolver%implicitly
            call ComputeTerms( deg_plus )
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
         do i = 1, grid%nelem
            elem => grid%elem(i)
            call Elem_wSTfinToRhsST( grid%elem(i) , state%nlSolver%implicitly)

            !write(91, *) elem%i, elem%wSTfinAD(k,1:dof)

         enddo !i


      class default
      stop 'ComputeSTDGM_Terms, For TimeTDG_t only!!!'
   end select
 end associate ! time

 ! putting back the original value
 state%time%ctime = cTime

 !call cpu_time(tt)
 !write(*,'(a30, 2es12.4)') 'CT_ends   STDGM :', tt - state%start_time ,  tt - tt1
 !write(*,'(a6,2i5)') '###'

 !  print*, '________________________________________'
 !  print*, 'End of ComputeSTDGM_Terms with implicitly = ', state%nlSolver%implicitly


  end subroutine ComputeCplus


  !> POROUS MEDIA FLOW, compute the HO continuous reconstruction of the capacity and conductivity
  subroutine ComputeCapacityConductivity( )
    class(element), pointer :: elem
    real, dimension(:,:), allocatable :: wi
    real, dimension(:,:,:), allocatable :: TA
    integer :: itype, Qdof, i

    open(23, file ='conductivity', status='unknown')

    ! element-wise (discontinuous) reconstruction
    do i=1, grid%nelem
       elem => grid%elem(i)

       Qdof = elem%Qdof

       ! solution in integ nodes
       allocate(wi(1:Qdof, 1:ndim))
       call Eval_w_Elem(elem, wi(1:Qdof,1:ndim) )

       ! setting of matrixes TA in integration nodes
       allocate(TA(1:Qdof, 1:ndim,1:ndim))

       ! array for the reconstruction
       allocate(elem%wSD(1:2, 1:ndim, 1:elem%dof_plus) )

       ! setting of capacity
       call Set_Time_Matrix_porous(elem, ndim,  Qdof, wi(1:Qdof,1:ndim), elem%xi(0,1:Qdof, 1:2+iRe), &
            TA(1:Qdof, 1:ndim,1:ndim), elem%wSD(1, 1:ndim, 1:elem%dof_plus) )


       call PlotElem_D_Function3D(120 , elem, elem%dof_plus, &
            elem%wSD(1, 1, 1:elem%dof_plus) )

       deallocate(wi, TA) ! not used acctually
    enddo

    close(23)
    ! local continuous reconstructions, from elem%wSD stored in elem%wS
    itype = 2   ! reconstruction of the conductivity and capacity
    call ComputeHO_LocalProblems( itype )


    do i=1, grid%nelem
       elem => grid%elem(i)
       deallocate(elem%wSD )
       deallocate(elem%wS )
    enddo

  end subroutine ComputeCapacityConductivity

end module euler_problem




