!> aposteriori error estimation based on the solution of local Neumann problems
module neumann_estim
  use paramets
!  use main_data
  use ama_L2interpol
  use data_mod ! contains type(mesh) :: grid for computation
  use mesh_oper
  use problem_oper
  use euler_problem
  use rav_tho_nedelec
  use loc_rav_tho_ned
  use eval_rav_tho_ned
  use dual_estim
  use inviscid_fluxes
  use f_mapping
  use blocks_integ
  use basis
  use eval_sol
  use errorDual

  implicit none

  public:: ComputeLocalNeumannEstim
  public:: LocalNeumannVertexProblem
  public:: Eval_DGphi_elem
  public:: Eval_RTNphi_elem
  public:: Eval_tau_g_elem
  public:: Integ_dot_product_RTNphi_elem
  public:: Integ_divRTNphi_DGphi_elem
  public:: ComputeLocRTNMomentumsElem2
  public:: Eval_p_rob_estims
  public:: Eval_BC_estim

contains
  !> evaluate the error estimates based on the Helmholtz decomposition
  subroutine  ComputeLocalNeumannEstim( )
    class(element), pointer :: elem, elem1
    logical, dimension(:), allocatable :: inner ! inner vertex
    integer, dimension(:), allocatable :: N  ! number of elements sharing a vertex
    integer, dimension(:,:,:), allocatable :: supp !list of corresponding elements
    real, dimension(:,:), allocatable :: regul !list of corresponding elements
    integer :: i_var, is, ismoothing
    integer :: maxdeg = 20
    integer :: i,j, dof
    real :: normF, normS, weightE, weightR

    !print*,' # starting of ComputeLocalNeumannEstim( )'
    ! create a list of elements sharing a vertex
    allocate( inner(1:grid%npoin), N(1:grid%npoin))
    allocate( supp(1:grid%npoin, 1:maxdeg, 1:2) )

    call SeekVertexSupports(grid, maxdeg, inner, N, supp)

    !write(201,*) '***********************', state%space%adapt%adapt_level
    !write(301,*) '***********************', state%space%adapt%adapt_level
    !write(401,*) '***********************', state%space%adapt%adapt_level

    ! ! graphical verification
    ! do i=1,grid%npoin
    !    write(*,'(a4,2i5,l3, 30i5)' ) '###',i,N(i),inner(i), supp(i,1: N(i) )

    !i = 15
    !do j=1, abs(N(i))
    !   elem => grid%elem(supp(i,j,1))
    !   
    !   write(2000+i, *) grid%x(i,:)
    !   write(2000+i, *) elem%xc(1:2)
    !   write(2000+i, '(x)' ) 
    !enddo
    ! enddo

    !print*,'#################################'

    do i_var = 0, 1   ! 0 => p-1, 1=> p
       do i=1,grid%nelem
          elem => grid%elem(i)
          dof = DOFtriang( MaxDegreeImplemented ) !???
          ! degree of reconstruction depends also on the neighbour
          if(i_var == 0) &
               allocate( elem%RTNphi(1:dof, 1:5) ) ! flux reconstruction on elem in DG basis functions
          elem%RTNphi(:,:) = 0.
          elem%RTNflux_deg = 0 

          if(i_var == 0) then
             elem%eta( : , :) = 0.
          else
             elem%eta(1:P_s_p4, :) = 0.
          endif

       enddo
       
       ! solution of local Neumann problem for each vertex
       do i=1,grid%npoin
          !do i=21, 21
          call LocalNeumannVertexProblem(i_var, i, inner(i), N(i), supp(i, 1:N(i),1:2)) 
       enddo

       !print*,'ATTENTION IN WERTYSDF'
       
       ! compute ESTIMATES
       state%estim(:, :)  = 0.
       do i=1,grid%nelem
          elem => grid%elem(i)
          call Eval_p_rob_estims(i_var, elem, normF, normS)
          
          state%estim(P_Flux:P_potP, :) = state%estim(P_Flux:P_potP, :) + elem%eta(P_Flux:P_potP, :)**2
          
          state%estim(P_F_p1:P_F_p4, :) = state%estim(P_F_p1:P_F_p4, :) &
               + (elem%eta(P_F_p1:P_F_p4, :)*normF)**2
          
          state%estim(P_s_p1:P_s_p4, :) = state%estim(P_s_p1:P_s_p4, :) &
               + (elem%eta(P_s_p1:P_s_p4, :)*normS)**2
          
          !state%estim(:, :) = state%estim(:, :) + elem%eta(:, :)**2
          
          !if(i==1) &
          !     print*,'##S#E', state%estim(P_tot, 1),  state%estim(P_potP, 1), &
          !     write(*,'(a6,120es9.1)') 'EDE#@', state%estim(:,1)
          !print*
          
          if(i_var == 1) deallocate( elem%RTNphi ) 
       enddo
    !print*,'##S#E', state%estim(P_pot, 1),  state%estim(P_potP, 1)

       write(*,'(a2, 6(a10, es9.2))') 'Ee',&
            'P_Rez = ', state%estim(P_rez, 1)**0.5 , &
            'P_Flux = ', state%estim(P_flux, 1)**0.5 , &
            'P_FR = ', state%estim(P_FR, 1)**0.5 , &
            'P_pot = ', state%estim(P_pot, 1)**0.5, & 
            'P_BC = ', state%estim(P_BC, 1)**0.5, & 
            'P_tot = ', state%estim(P_tot, 1)**0.5
       !print*,state%estim(P_potP, 1)**0.5  !, grid%elem(1)%eta(P_potP, 1)
    enddo ! i_var

    deallocate( inner, N, supp )
       

    ! values for adaptation
    do i = 1, grid%nelem
       elem => grid%elem(i)
       elem%estimST =  elem%eta(P_tot, 1)
       
       ! indicator of the regularity
       elem%reg =  1E-15
       elem%reg1 = 1E-15
       elem%reg2 = 1E-15
       
       ! 1st variant
       !if(elem%deg >= 2) elem%reg  = elem%eta(P_F_p1, 1) / elem%eta(P_F_p2, 1) !/  elem%diam**0.5
       !if(elem%deg >= 3) elem%reg1 = elem%eta(P_F_p2, 1) / elem%eta(P_F_p3, 1) !/  elem%diam**0.5

       ! Babuska variant
       if(elem%deg >= 1) elem%reg  = elem%eta(P_tot, 1) / elem%eta(P_potP, 1) !/  elem%diam**0.5
       !if(elem%deg >= 2) elem%reg1 = elem%eta(P_F_p1, 1) / elem%eta(P_F_p2, 1) /  elem%diam**0.5
       !if(elem%deg >= 3) elem%reg2 = elem%eta(P_F_p2, 1) / elem%eta(P_F_p3, 1) /  elem%diam**0.5

       ! the following already NECESSARY 
       !elem%regT0 = 0.75 * elem%diam**0.5       ! elem%reg <  elem%regT0  ==> p-refinement
       elem%regT0 = 0.75 * elem%diam**0.5       ! elem%reg <  elem%regT0  ==> p-refinement
       !!!elem%regT0 = 0.95 * elem%diam**0.5 ! case J      ! elem%reg <  elem%regT0  ==> p-refinement
       !elem%regT0 = 1.5*elem%diam**0.5       ! elem%reg <  elem%regT0  ==> p-refinement
       !elem%regT1 = 1.                   ! elem%reg <  elem%regT1  ==> h->p substitution
       elem%regT2 = elem%diam**(-0.5)    ! elem%reg >  elem%regT2  ==> p-DErefinement

       !elem%regT0 = elem%diam**1.0                ! elem%reg <  elem%regT0  ==> p-refinement
       !elem%regT1 = elem%diam**(0.5)    ! elem%reg <  elem%regT1  ==> h->p substitution
       !elem%regT2 = 1.  !0.9  !elem%diam**(0.0)    ! elem%reg >  elem%regT2  ==> p-DErefinement

    enddo

    
    ! smoothing of the regul parameter
    allocate(regul(1:grid%nelem, 1:12) )
    
    ismoothing = 1 !1 ! number of smoothing cycles
    weightE = 0.0  !0.05  ! weight for error estim of the neighbouring elements
    weightR = 0.25 !0.75  ! weight for regularity of the neighbouring elements

    do is = 1, ismoothing
       regul(1:grid%nelem, 1) = grid%elem(1:grid%nelem)%reg
       regul(1:grid%nelem, 2) = grid%elem(1:grid%nelem)%reg1
       regul(1:grid%nelem, 3) = grid%elem(1:grid%nelem)%reg2

       regul(1:grid%nelem, 6) = grid%elem(1:grid%nelem)%estimST

       regul(1:grid%nelem, 11) = 1.
       regul(1:grid%nelem, 12) = 1.

       do i = 1, grid%nelem
          elem => grid%elem(i)
          do j=1,elem%flen
             if(elem%face(neigh,j) > 0) then
                elem1 => grid%elem(elem%face(neigh,j))

                regul(i,1) = regul(i,1) + weightR * elem1%reg
                regul(i,2) = regul(i,2) + weightR * elem1%reg1
                regul(i,3) = regul(i,3) + weightR * elem1%reg2
                regul(i,11) = regul(i,11) + weightR

                regul(i,6) = regul(i,6) + weightE * elem1%estimST
                regul(i,12) = regul(i,12) + weightE
             endif
          enddo
       enddo

       do i = 1, grid%nelem
          elem => grid%elem(i)

          !write(*,'(a5,i5, 3es12.4)') 'DESW',i, elem%reg, elem%reg1
          elem%reg     = regul(i,1) / regul(i,11)
          !elem%reg = max ( elem%reg , regul(i,1) / regul(i,11))

          elem%reg1    = regul(i,2) / regul(i,11)
          elem%reg2    = regul(i,3) / regul(i,11)

          elem%estimST = regul(i,6) / regul(i,12)
          !write(*,'(a5,i5, 3es12.4)') 'DESW',i, elem%reg, elem%reg1
          !print*

       enddo
    enddo ! do is

    deallocate(regul)

    !print*,'Stopped in Neumann'
    !stop



100 format(a6,2es12.4,'|',14es12.4)
!    write(22,100) &
!         '$$$',state%space%h, state%time%tau(1), estimL(Hrez:Heta1, 1:ndim)



    !print*,' # finish of ComputeLocalNeumannEstim( )', state%estim(P_tot, :)**0.5

  end subroutine ComputeLocalNeumannEstim


  !> evaluate the error estimates based on the solution of local Neumann problems
  !> RTN flux and H^1 potential were already evaluated
  subroutine Eval_p_rob_estims(i_var, elem, normF, normS)
    integer, intent(in) :: i_var   ! index of variant
    class(element), intent(inout) :: elem
    real :: normF, normS
    type(volume_rule), pointer :: V_rule
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:,:), allocatable :: Dphi
    real, dimension(:,:), allocatable :: Dwi
    real, dimension(:,:), allocatable :: Dpot, DpotP
    real, dimension(:,:), allocatable :: flux
    real, dimension(:,:), allocatable :: Dflux
    real, dimension(:,:), allocatable :: Fx
    real, dimension(:,:), allocatable :: ff
    real, dimension(:), allocatable :: weights
    real, dimension(:,:), allocatable :: wExact
    real, dimension(:,:,:), allocatable :: DwExact

    integer :: Qnum, Qdof, i, Fdeg, dof, Fdof, max_dof, j, l, Fdof1
    integer :: k1, k2, ib, i1, i2
    integer :: FdofL, FdofU
    real :: pi, val, weigh  !!, normS, normF

    pi = 2* acos(0.)

    Fdeg = elem%RTNflux_deg
    Fdof = DOFtriang( Fdeg )  

    FdofU = DOFtriang( Fdeg -2 )  

    !if(elem%i == 1) print*,'EEEEEEEEEEEE', Fdeg, Fdof, FdofU

    Qnum =  state%space%Qdeg(Fdeg+1, 1) 
    V_rule => state%space%V_rule(Qnum)

    Qdof = V_rule%Qdof
    
    dof = elem%dof
    max_dof = max (dof, Fdof)

    allocate(wExact(1:Qdof, 1:ndim) )
    allocate(DwExact(1:Qdof, 1:ndim, 1:nbDim))
    call SetExactSolutionQnodes(elem, V_rule, wExact(1:Qdof, 1:ndim), DwExact(1:Qdof, 1:ndim, 1:nbDim))

    !print*,'###',dof, Fdof, max_dof
    !write(*,'(a4,30es10.2)') 'RTN1',elem%RTNphi(1:Fdof, 1)
    !write(*,'(a4,30es10.2)') 'RTN2',elem%RTNphi(1:Fdof, 2)
    !write(*,'(a4,30es10.2)') 'RTN5',elem%RTNphi(1:Fdof, 5)

    allocate(weights(1:Qdof) )
    call Eval_V_Weights_plus(elem, V_rule, weights(1:Qdof) )

    phi => V_rule%phi(1:max_dof, 1:Qdof)

    allocate(Dphi(1:max_dof, 1:2, 1:Qdof)) 
    call Eval_Dphi_plus(elem, V_rule, max_dof,  Dphi(1:max_dof, 1:2, 1:Qdof) )

    allocate( Dwi(1:Qdof, 0:2)  )
    Dwi(1:Qdof, 0) = matmul(elem%w(0,1:dof), phi(1:dof, 1:Qdof) )
    Dwi(1:Qdof, 1) = matmul(elem%w(0,1:dof), Dphi(1:dof, 1, 1:Qdof) )
    Dwi(1:Qdof, 2) = matmul(elem%w(0,1:dof), Dphi(1:dof, 2, 1:Qdof) )

    ! RTN flux
    allocate( flux(1:Qdof, 1:2)  )
    flux(1:Qdof, 1) = matmul(elem%RTNphi(1:Fdof, 1), phi(1:Fdof, 1:Qdof) )
    flux(1:Qdof, 2) = matmul(elem%RTNphi(1:Fdof, 2), phi(1:Fdof, 1:Qdof) )
 
    ! div of RTN flux
    allocate( Dflux(1:Qdof, 1:2)  )
    Dflux(1:Qdof, 1) = matmul(elem%RTNphi(1:Fdof, 1), Dphi(1:Fdof, 1, 1:Qdof) )
    Dflux(1:Qdof, 2) = matmul(elem%RTNphi(1:Fdof, 2), Dphi(1:Fdof, 2, 1:Qdof) )
 
    ! potential and its gradient
    allocate( Dpot(1:Qdof, 0:2)  )
    Dpot(1:Qdof, 0) = matmul(elem%RTNphi(1:Fdof,5), phi(1:Fdof, 1:Qdof) )

    ! OLD VARIANT but working variant
    !Dpot(1:Qdof, 1) = matmul(elem%RTNphi(1:Fdof,5), Dphi(1:Fdof, 1, 1:Qdof) )
    !Dpot(1:Qdof, 2) = matmul(elem%RTNphi(1:Fdof,5), Dphi(1:Fdof, 2, 1:Qdof) )

    ! NEW VARIANT
    Dpot(1:Qdof, 1) = matmul(elem%RTNphi(1:Fdof,3), phi(1:Fdof, 1:Qdof) )
    Dpot(1:Qdof, 2) = matmul(elem%RTNphi(1:Fdof,4), phi(1:Fdof, 1:Qdof) )

    ! projection of the potential
    allocate( DpotP(1:Qdof, 0:2)  )
    DpotP(1:Qdof, 1) = matmul(elem%RTNphi(1:FdofU,5), Dphi(1:FdofU, 1, 1:Qdof) )
    DpotP(1:Qdof, 2) = matmul(elem%RTNphi(1:FdofU,5), Dphi(1:FdofU, 2, 1:Qdof) )

    !do l=1,Fdof
    !   write(*,'(a4,i5,40es14.6)') 'DEWD',l, elem%RTNphi(l,3:4)
    !   !Dpot(l, 1:2),DpotP(l, 1:2)
    !enddo


    ! setting of source terms in integ nodes
    allocate( Fx(1:Qdof, 1:nbDim) )
    !integration nodes on K
    call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof,1:nbDim), &
         Fx(1:Qdof, 1:nbDim) )
       
    allocate( ff(1:Qdof, 1:7*ndim) )
    if(elem%i ==1 .and. ndim > 1) print*,'TROUBLE in nemuann.F0000'

    do j=1,Qdof
       call RHS_Scalar(Fx(j, 1:2), ff(j, 1:ndim), state%time%ctime)

       !if(vectornorm(Dpot(j,1:2) -Dpot(j,3:4)) > 1E-10) &
       !!     write(200,'(8es12.4)') &
       !     write(*,'(a5,2i5,8es12.4)')'pot:',elem%i,j,&
       !     Fx(j,1:2), abs(Dpot(j,1:2) -Dpot(j,3:4)),  Dpot(j,1:2),Dpot(j,3:4)

       !write(600+elem%i,*) Fx(j,1:2), Dwi(j, 1:2), flux(j,1:2), Dflux(j, 1:2), ff(j,1)
       !write(700,*) Fx(j,1:2), Dpot(j,0:2) , Dwi(j, 1:2), flux(j,1:2), Dflux(j, 1:2), ff(j,1)
    enddo

    !call DRAW_reconstruction(800, elem, 8, Fdof, 5)


    ! indicator P_flux = || \nabla u_h + \sigma_h ||
    ff(1:Qdof, 2) =  (Dwi(1:Qdof, 1) + flux(1:Qdof, 1))**2  &
         + (Dwi(1:Qdof, 2) + flux(1:Qdof, 2))**2 

    elem%eta(P_flux, 1) = dot_product(weights(1:Qdof),  ff(1:Qdof, 2) )**0.5

    
    ! indicator P_rez = h_K/pi || f - \nabla \cdot \sigma_h ||
    ff(1:Qdof, 3) =  (ff(1:Qdof, 1) - Dflux(1:Qdof, 1) - Dflux(1:Qdof, 2)  )**2

    elem%eta(P_rez, 1) = dot_product(weights(1:Qdof),  ff(1:Qdof, 3) )**0.5
    elem%eta(P_rez, 1) = elem%eta(P_Rez, 1) * elem%diam / pi


    elem%eta(P_FR, 1) = elem%eta(P_rez, 1) + elem%eta(P_flux, 1) 

    ! indicator P_pot = || \nabla u_h  - \nabla s_h ||  S_h = potential

    ! discrete gradient with the  lifting operator for the dicrete gradient 
    ! in pNeu for SIPG and NIPG [Ern, Vohralik, SINUM 15]
    Dwi(1:Qdof,1) = Dwi(1:Qdof,1) + state%space%m_IPG*elem%lifting(1, 1) ! opposite sign !
    Dwi(1:Qdof,2) = Dwi(1:Qdof,2) + state%space%m_IPG*elem%lifting(1, 2)


    ff(1:Qdof, 4) =  (Dpot(1:Qdof, 1) -  Dwi(1:Qdof, 1))**2 &
         + (Dpot(1:Qdof, 2) -  Dwi(1:Qdof, 2))**2 

    elem%eta(P_pot, 1) =  dot_product(weights(1:Qdof),  ff(1:Qdof, 4))**0.5

    !! projection "indicator"  P_pot = || \nabla u_h  - \Pi \nabla s_h ||  S_h = potential
    ! not used at this moment
    !ff(1:Qdof, 4) =  (DpotP(1:Qdof, 1) -  Dwi(1:Qdof, 1))**2 &
    !     + (DpotP(1:Qdof, 2) -  Dwi(1:Qdof, 2))**2 
    !
    !elem%eta(P_potP, 1) =  dot_product(weights(1:Qdof),  ff(1:Qdof, 4))**0.5


    !if(abs(elem%eta(P_pot,1) - elem%eta(P_potP,1)) > 1E-3) &
    !     print*,'#EDE', elem%i, elem%eta(P_pot,1), elem%eta(P_potP,1)

    !print*,'#. # ', 0., 0.,     elem%eta(P_pot,1), elem%eta(P_potP,1)

    ! BC indicator
    call  Eval_BC_estim(elem, V_rule)
    

    !print*,'### eta(P_BC) ', elem%eta(P_BC, 1)

    ! total indicator
    elem%eta(P_tot, 1) = (elem%eta(P_FR, 1)**2  + (elem%eta(P_pot,1) +elem%eta(P_BC,1) )**2)**0.5

    !!!!!elem%eta(P_sF: max_eta , : ) = 0.; normS = 1.; normF = 1.

    ! NEW variant ala Babuska
    if(i_var == 0) elem%eta(P_potP, 1) = elem%eta(P_tot, 1)
    
    if(i_var == 0) return

    return

    ! TESTING parameters
    ! indicator  || \nabla s +  \sigma_h ||
    ff(1:Qdof, 5) =   (Dpot(1:Qdof, 1) +  flux(1:Qdof, 1))**2 &
         + (Dpot(1:Qdof, 2) +  flux(1:Qdof, 2))**2 

    elem%eta(p_sF, 1) =  dot_product(weights(1:Qdof),  ff(1:Qdof, 5))**0.5

    ! indicator  || \nabla s - \nabla u ||
    ff(1:Qdof, 6) =   (Dpot(1:Qdof, 1) - DwExact(1:Qdof,1, 1))**2 &
         + (Dpot(1:Qdof, 2)  - DwExact(1:Qdof,1, 2) )**2 

    elem%eta(p_su, 1) =  dot_product(weights(1:Qdof),  ff(1:Qdof, 6))**0.5

    ! indicator  || \nabla u +sigma ||
    ff(1:Qdof, 7) =   (flux(1:Qdof, 1) + DwExact(1:Qdof,1, 1))**2 &
         + (flux(1:Qdof, 2)  + DwExact(1:Qdof,1, 2) )**2 

    elem%eta(p_FDu, 1) =  dot_product(weights(1:Qdof),  ff(1:Qdof, 7))**0.5

    ! norm  || \nabla s_h||
    ff(1:Qdof, 7) =  Dpot(1:Qdof, 1)**2 + Dpot(1:Qdof, 2)**2 
    normS =  dot_product(weights(1:Qdof),  ff(1:Qdof, 7))**0.5

    ! norm  || sigma_h||
    ff(1:Qdof, 7) = flux(1:Qdof, 1)**2 + flux(1:Qdof, 2)**2
    normF =  dot_product(weights(1:Qdof),  ff(1:Qdof, 7))**0.5


    ! indicator  || sigma - \Pi^{p-l} sigma ||
    call ElemRTNFluxProjection(elem, 4)
    elem%eta(P_F_p1 : P_F_p4, 1) = elem%eta(p_F_p1 : p_F_p4, 1) / max (1E-15, normF)

    ! indicator  ||\nabla s - \Pi^{p-l} \nabla s ||
    !do l = 1, min(Fdeg , 4)
    do l = 0, min(Fdeg+1 , 4)
       
       ! approach using gradient of the potential $ s $ evaluated in local Neumann problems
       ! FdofL = 1
       ! FdofU = DOFtriang( Fdeg - l)  

       ! ff(1:Qdof, 5) = matmul(elem%RTNphi(FdofL:FdofU,5), Dphi(FdofL:FdofU, 1, 1:Qdof) )
       ! ff(1:Qdof, 6) = matmul(elem%RTNphi(FdofL:FdofU,5), Dphi(FdofL:FdofU, 2, 1:Qdof) )

       ! ! projection of the gradient of the potential
       ! ff(1:Qdof, 7) =   (Dpot(1:Qdof, 1) - ff(1:Qdof, 5))**2 &
       !   + (Dpot(1:Qdof, 2)  - ff(1:Qdof, 6) )**2 

       !! projection of potential 
       !!FdofL = DOFtriang( Fdeg - l)  
       !!FdofU = DOFtriang( Fdeg )  

       ! approach using the potential $ s $ itself evaluated in local Neumann problems
       FdofL = 1
       FdofU = DOFtriang( Fdeg-l )  


       ff(1:Qdof, 7) = matmul(elem%RTNphi(FdofL:FdofU,5), phi(FdofL:FdofU, 1:Qdof) )
       ff(1:Qdof, 7) = (Dpot(1:Qdof, 0) - ff(1:Qdof, 7))**2

       normS = dot_product(weights(1:Qdof),   Dpot(1:Qdof, 0)**2 )**0.5

       ! value of the projection (potential or its gradient)
       if(l > 0) &
            elem%eta(P_s_p1-1+l, 1) =  &
            dot_product(weights(1:Qdof),  ff(1:Qdof, 7))**0.5 / max (1E-15, normS)

       !if(elem%i == 1) write(*,'(a4,5i5,3es18.10)')'WED',elem%i, l, Fdof, FdofL, FdofU, elem%eta(P_s_p1-1+l, 1), &
       !     dot_product(weights(1:Qdof),  ff(1:Qdof, 7))**0.5 , normS

    enddo


    ! if(elem%i <= 2 .or. (elem%i > 200 .and. elem%i < 203) .or. &
    !       (elem%i > 500 .and. elem%i < 503) ) then
    !    !write(*,'(a5,i5,4es12.4, a2, 4es12.4)') &
    !    !     'Pi:',elem%i,  elem%eta(p_F_p1 : p_F_p4, 1),'|', &
    !    !     elem%eta(p_s_p1 : p_s_p4, 1)
    !    write(*,'(a5,i5,4es12.4, a2, 6es12.4)') &
    !         'Pi:',elem%i,  elem%eta(p_s_p1 : p_s_p4, 1),'|', &
    !         elem%eta(p_s_p2,1) / elem%eta(p_s_p1,1), &
    !         elem%eta(p_s_p3,1) / elem%eta(p_s_p2,1), &
    !         elem%eta(p_s_p4,1) / elem%eta(p_s_p3,1), &
    !         elem%eta(p_s_p2,1) / elem%eta(p_s_p1,1)*elem%diam, &
    !         elem%eta(p_s_p3,1) / elem%eta(p_s_p2,1)*elem%diam, &
    !         elem%eta(p_s_p4,1) / elem%eta(p_s_p3,1)*elem%diam
            
    ! endif


    !write(99,*)  elem%i,  elem%eta(p_s_p1 : p_s_p4, 1), & ! 1,2..5
    !     elem%eta(p_s_p1,1) / elem%eta(p_s_p2,1), &       ! 6
    !     elem%eta(p_s_p2,1) / elem%eta(p_s_p3,1), &       ! 7
    !     elem%eta(p_s_p3,1) / elem%eta(p_s_p4,1), &       ! 8
    !     elem%eta(p_s_p1,1) / elem%eta(p_s_p2,1)/elem%diam, &
    !     elem%eta(p_s_p2,1) / elem%eta(p_s_p3,1)/elem%diam, &
    !     elem%eta(p_s_p3,1) / elem%eta(p_s_p4,1)/elem%diam

    !do j=1,Qdof
    !    if(elem%i >= 8 .and. elem%i <= 10) &
    !         write(400,*) Fx(j,1:2), Dpot(j,0:2) , & ! 1.. 5
    !         Dwi(j, 0:2), flux(j,1:2), Dflux(j, 1:2),    & ! 6..12
    !         ff(j,1:4),                                   & !13..16 
    !         elem%eta(1:P_tot, 1),Dwi(j, 1:2)+ flux(j,1:2),&!17..23
    !         DwExact(j, 1, 1:nbDim)      ,                & !24..25 
    !         Dflux(j, 1) +  Dflux(j, 2)                     ! 26..

    !    if(elem%i < 100) &
    !         write(500+elem%i,*) Fx(j,1:2), Dpot(j,0:2) , & ! 1.. 5
    !         Dwi(j, 0:2), flux(j,1:2), Dflux(j, 1:2),    & ! 6..12
    !         ff(j,1:4),                                   & !13..16 
    !         elem%eta(1:P_tot, 1),Dwi(j, 1:2)+ flux(j,1:2),&!17..23
    !         DwExact(j, 1, 1:nbDim)      ,                & !24..25 
    !         Dflux(j, 1) +  Dflux(j, 2)                     ! 26..

    !    write(700,*) Fx(j,1:2), Dpot(j,0:2) , & ! 1.. 5
    !         Dwi(j, 0:2), flux(j,1:2), Dflux(j, 1:2),    & ! 6..12
    !         ff(j,1:4),                                   & !13..16 
    !         elem%eta(1:P_tot, 1),Dwi(j, 1:2)+ flux(j,1:2),&!17..23
    !         DwExact(j, 1, 1:nbDim)           ,           & !24..25 
    !         Dflux(j, 1) +  Dflux(j, 2)                     ! 26..
    ! enddo

    !write(*,'( 2(a12, es12.4), a12,4es12.4)')  &
    !      'P_Rez = ', elem%eta(P_rez, 1) , &
    !      'P_Flux = ', elem%eta(P_flux, 1) , &
    !      'orthog:',&
    !      dot_product(weights(1:Qdof),  ff(1:Qdof, 2) *  phi(1, 1:Qdof)), &
    !      dot_product(weights(1:Qdof),  ff(1:Qdof, 2) *  phi(2, 1:Qdof)), &
    !      dot_product(weights(1:Qdof),  ff(1:Qdof, 2) *  phi(3, 1:Qdof))
         
    deallocate(weights, Dphi, flux, Dflux, ff, Fx, Dpot, DpotP)
       
  
   !print*,'Stopped in Neumann, subroutine Eval_p_rob_estims', elem%i
   !stop
   

  end subroutine Eval_p_rob_estims


  !> array elem%RTNphi(1:Fdof, 1) contains flux from RTN_p, p = elem%RTNflux_deg
  !> we compute error if its projection into RTN_{p-l}, l=1,2,...max_lev
  subroutine ElemRTNFluxProjection(elem, max_lev)
    class(element), intent(inout) :: elem
    integer, intent(in) :: max_lev
    type(volume_rule), pointer :: V_rule
    real, dimension(:,:,:),  allocatable:: RTNphi
    real, dimension(:,:,:),  allocatable:: flux
    real, dimension(:,:),  allocatable:: MM
    real, dimension(:),  allocatable:: weights
    real, dimension(:,:),  pointer:: phi
    integer :: FFdeg, FFdof
    integer :: Qdof, i, Fdeg, dof, Fdof, max_dof, j, l, Fdof1
    real :: pi, normS, normF
    integer :: FdofL, FdofU
   
    V_rule => state%space%V_rule(state%space%Qdeg(elem%RTNflux_deg, 1) )
    !V_rule => state%space%V_rule(elem%Qnum)

    Qdof = V_rule%Qdof

    allocate(weights(1:Qdof) )
    call Eval_V_Weights_plus(elem, V_rule, weights(1:Qdof) )

    Fdeg = elem%RTNflux_deg - 1 
    Fdof = DOFtriang( Fdeg )  
    FFdof = SetRTNdof( Fdeg )  

    ! RTN basis function in integ nodes
    allocate( RTNphi(1:FFdof, 1:3, 1:Qdof) ) 

    ! mass matrix for the projection
    allocate( MM(1:FFdof, 1:FFdof+2) ) 

    dof = elem%dof
    max_dof = max (dof, Fdof)

    phi => V_rule%phi(1:max_dof, 1:Qdof)

    ! RTN flux in integ nodes
    allocate( flux(0:2, 1:Qdof, 1:2)  )

    flux(0,1:Qdof, 1) = matmul(elem%RTNphi(1:Fdof, 1), phi(1:Fdof, 1:Qdof) )
    flux(0,1:Qdof, 2) = matmul(elem%RTNphi(1:Fdof, 2), phi(1:Fdof, 1:Qdof) )


    !print*,'#EDE',max_lev, Fdeg
    do l = 1, min(max_lev , Fdeg)

       FFdeg = Fdeg - l
       Fdof = DOFtriang( FFdeg  )  
       FFdof = SetRTNdof( FFdeg )  
    
       call Eval_RTNphi_elem(elem, V_rule, FFdeg, FFdof, RTNphi(1:FFdof, 1:3, 1:Qdof), 0  )

       !do i=1,Qdof
       !   write(97,*) V_rule%lambda(i,1:2),  RTNphi(1:FFdof, 1,i)
       !   write(98,*) V_rule%lambda(i,1:2),  RTNphi(1:FFdof, 2,i)
       !enddo
    
       do i=1,FFdof
          do j=i, FFdof
             MM(i,j) = dot_product( weights(1:Qdof), &
                  RTNphi(i, 1, 1:Qdof)*RTNphi(j, 1, 1:Qdof) &
                  + RTNphi(i, 2, 1:Qdof)*RTNphi(j, 2, 1:Qdof) )

             MM(j,i) = MM(i,j)
          enddo
          MM(i,FFdof+1) = dot_product( weights(1:Qdof), &
               RTNphi(i, 1, 1:Qdof)* flux(0, 1:Qdof, 1) &
               + RTNphi(i, 2, 1:Qdof)*flux(0, 1:Qdof, 2) )
       enddo

       ! coefficients of the projection

       !call MblockInverse(FFdof, MM(1:FFdof, 1:FFdof) )
       !MM(1:FFdof, FFdof+2) = matmul( MM(1:FFdof, 1:FFdof),  MM(1:FFdof, FFdof+1))
       
       call SolveLocalMatrixProblem(FFdof, MM(1:FFdof, 1:FFdof), 1, MM(1:FFdof, FFdof+1) ) !:FFdof+1) )
       MM(1:FFdof, FFdof+2) = MM(1:FFdof, FFdof+1)

       ! projection in integ nodes
       flux(1, 1:Qdof, 1) = matmul( MM(1:FFdof, FFdof+2),  RTNphi(1:FFdof, 1, 1:Qdof) )
       flux(1, 1:Qdof, 2) = matmul( MM(1:FFdof, FFdof+2),  RTNphi(1:FFdof, 2, 1:Qdof) )

       !write(*,'(a6,140es12.4)') 'coefOr', flux(0, 1:Qdof, 1:2)
       !write(*,'(a6,140es12.4)') 'coefNEW', flux(1, 1:Qdof, 1:2)


       ! || \sigma - \Pi_{-l} \sigma ||
       flux(2,1:Qdof, 1) =  (flux(0,1:Qdof, 1) - flux(1, 1:Qdof, 1))**2 &
            + (flux(0,1:Qdof, 2) - flux(1, 1:Qdof, 2))**2 

       ! || \Pi_{-l} \sigma ||
       !flux(2,1:Qdof, 1) =  flux(1, 1:Qdof, 1)**2  + flux(1, 1:Qdof, 2)**2 



       if(l > 0) elem%eta(P_F_p1 -1 + l, 1) =  &
            dot_product(weights(1:Qdof),  flux(2,1:Qdof, 1))**0.5 !!!/ normF

       !print*,'###',l,elem%eta(P_F_p1 -1 + l, 1), dot_product(weights(1:Qdof),  flux(2,1:Qdof, 1))**0.5, &
       !Fdeg, Fdof, FFdof,'|',elem%Qnum, V_rule%Qdeg, Qdof

       !do i=1,Qdof
       !   write(100-l,*) V_rule%lambda(i,1:2),  flux(0, i, 1:2),  flux(1, i, 1:2)
       !enddo
    enddo

    deallocate(RTNphi, flux, MM, weights)

    !stop

  end subroutine ElemRTNFluxProjection


  !> solve the local Neumann problem on the patch corresponding to a vortex of the mesh
  subroutine LocalNeumannVertexProblem(i_var, ip, inner, N, supp )
    integer, intent(in) :: i_var   ! index of variant
    integer, intent(in) :: ip      ! index of node
    logical, intent(in) :: inner   ! inner vertex?
    integer, intent(in) :: N   ! number of elememnts sharing vertex
    integer, dimension(1:N,1:2), intent(in) :: supp   !idx of elems sharing vertex
    class(element), pointer :: elem
    type(volume_rule), pointer :: V_rule
    integer :: dN, dNR, dNF         ! dimension of the Neumann problem
    integer :: dN_pt, dNR_pt, dNF_pt     ! dimension of the Neumann problem
    real, dimension (:,:), allocatable  :: A, B  ! matrix for Neumann local problem
    real, dimension (:,:), allocatable  :: A_pt, B_pt  ! matrix for Neumann local problem

    real, dimension (:,:), allocatable  :: rhs, x ! RHSs and sols for Neumann local problem
    real, dimension (:,:), allocatable  :: rhs_pt, x_pt ! RHSs and sols for Neumann local problem
    real, dimension (:,:,:), allocatable  :: DGphi ! DG basis functions
    real, dimension (:,:,:,:), allocatable  :: RTNphi ! RTN basis functions
    real, dimension (:,:), allocatable  :: tau_g , xi ! RHS
    real, dimension (:,:), allocatable  :: Fx, A_RR, A_RD, uD
    real, dimension (:,:,:), allocatable  :: sigma
    integer, dimension (:,:, :), allocatable  :: itrans
    real :: area
    integer :: i, i_n, i_p, j, j1, k, kk, iBC, ie0 !, k1, k2, l1, l2 
    integer :: F_face, F_vol, F_tot, F_size, iptest
    integer :: Qdof, Qnum, Qdeg, ie, ie1, ie2
    integer :: deg, deg1, Fdof, Rdof, info
    logical :: innerR

    iptest = 1

    !innerR = .true.  ! simplified modification for non-homogenous BC
    innerR = inner   !  a rigorous modification for non-homogenous BC

    !if(inner) then ! inner vortex

    deg = maxval(grid%elem(supp(1:N,1))%deg) ! setting of p of reconstruction

    if(i_var == 0)  deg = deg - 1   ! low order econstruction for hp-variant
    
    !if(ip <= 2) 
    !print*,'!! ATTENTION, TEST in neumann.f90 !!!! ', deg
    if(deg <= 0) return

    !deg = 2
    !deg = 1
    !deg = 0

    deg1 = deg + 1
    
    Rdof = DOFtriang( deg )  
    Fdof = SetRTNdof( deg) 

    ! setting of num_quadrature
    !Qnum = state%space%Qdeg(deg + 3, 1)
    Qnum = state%space%Qdeg(deg + 2, 1)
    Qdeg = Qnum   

    V_rule => state%space%V_rule(Qnum)
    Qdof = V_rule%Qdof

    ! setting number of degrees of freedom of Q_h^a, V_h^a

    ! RTN DOF 
    F_face = deg1         ! for one face
    F_vol  = deg * deg1   ! for one element
    F_size =  2*(deg+1) + F_vol  ! opposite is always ignored

    F_tot = F_face + F_vol

    if(inner) then
       ! inner nodes
       dNR = N * Rdof  - 1    ! Q_h^a constraint
       dNF = N * (F_face + F_vol) ! V_h^a: 1 face +  1 vol per elem

    else
       ! boundary  nodes
       dNR = N * Rdof          ! Q_h^a NO constraint
       dNF = N * (F_face + F_vol) + F_face ! open list of adjacent elemenst

       ! for potential reconstruction
       dNR_pt = N * Rdof - 1         ! Q_h^a NO constraint
       dNF_pt = N * (F_face + F_vol) - F_face ! both edgech on \partial\omega are taken off
       dN_pt  = dNR_pt + dNF_pt
    endif

    dN = dNR + dNF

    allocate(A(1:dNF, 1:dNF+2), B(1:dNF+2, 1:dNR) )
    allocate( rhs(1:dN, 1:2), x(1:dN, 1:2) )
    A(:,:) = 0.
    B(:,:) = 0.

    ! reconstruction of the potetial for boundary faces has different dimension
    if(.not. innerR ) then
       allocate(A_pt(1:dNF_pt, 1:dNF_pt+2), B_pt(1:dNF_pt+2, 1:dNR_pt) )
       allocate( rhs_pt(1:dN_pt, 1:2), x_pt(1:dN_pt, 1:2) )
       A_pt(:,:) = 0.
       B_pt(:,:) = 0.

       ! nonhomogeneous Dirichlet BC
       allocate(uD(1: F_face, 1:2)  )
    endif 
       
    !test functions on the element
    allocate(DGphi(1:N, 1:Rdof, 1:Qdof) )  

    allocate(RTNphi(1:N, 1:Fdof, 1:3, 1:Qdof) )  !2nd idx: 1,2 -components, 3 = div

    allocate(A_RR(1:Fdof, 1:Fdof+2), A_RD(1:Fdof+2, 1:Rdof ) )

    allocate(tau_g(1:5, 1:Qdof), xi(1:3,1:2) )! 1st and 2nd parts of RHS

    area = sum(grid%elem(supp(:,1))%area)

    ! itrans(:,:, 1:2) inner OR flux,  itrans(:,:,3:4)  bound AND potential
    allocate(itrans(1:N, 1: F_size, 1:4) )  

    do i=1, N
       ie = supp(i,2)   ! inner index of the vertex

       elem => grid%elem(supp(i,1))

       ! DG basis fucntion with the constraint \int_\omega_a \phi_i dx = 0
       call Eval_DGphi_elem(elem, V_rule, Rdof, DGphi(i, 1:Rdof, 1:Qdof), area, inner)

       ! RTN basis functions with \q \cdot \nn = 0 pn \partial \omega_a
       call Eval_RTNphi_elem(elem, V_rule, deg, Fdof, &
            RTNphi(i, 1:Fdof, 1:3, 1:Qdof), ie)

       ! right hand side of the Neumann problem
       ! triangle coordinates for gradient of \psi_a
       xi(1:3, 1:2) = grid%x(elem%face(idx, 1:3), 1:2)
       call Eval_tau_g_elem(elem, i_var, V_rule, tau_g(1:5, 1:Qdof), supp(i,2), xi )

       call Integ_dot_product_RTNphi_elem(elem, V_rule,  Fdof, &
            RTNphi(i, 1:Fdof, 1:2, 1:Qdof), tau_g(1:5, 1:Qdof), &
            A_RR(1:Fdof, 1:Fdof+2) )

       call Integ_divRTNphi_DGphi_elem(elem, V_rule,  Fdof, Rdof,  &
            RTNphi(i, 1:Fdof, 3, 1:Qdof), DGphi(i, 1:Rdof, 1:Qdof), &
            tau_g(3, 1:Qdof), A_RD(1:Fdof+2, 1:Rdof) )

       
        ! if(ip == 1 .and. i==1 ) then
        !  !  write(89,*) grid%nelem
        !  !  write(89,'(2i5,10es14.6)') deg, 3*(deg1+1), A_RR(1,1), A_RR(3*deg1+1, 3*deg1+1)
       ! write(*,*) 'deg = ', deg, 'RTN matrix  ', 1,'..', 3*(deg+1),' *** ||  *** ',3*(deg+1) +1,'...', Fdof
       !     write(*,'(a6,400i12)') '   ',1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20
       !      do j=1,Fdof
       !         write(*,'(2i5,400es12.4)')i,j,A_RR(j,:)
       !      enddo
       !      print*,'--------------  A_RR matrix',i, A(1,2)
       !      stop
        ! endif

       ! if(ip == 1) then
       !     do j=1,Fdof+1
       !        write(*,'(2i3,20es12.4)')i,j,A_RD(j,:)
       !     enddo
       !     print*,'--------------  A_RG matrix'
        ! endif

       call CreateTransfPairs(deg, F_size, N, i, ie, itrans(i, 1: F_size, 1:4), inner )

       !print*,'-------A1',ip, i, elem%i
       call AssembMatrixRTN_RTN(A(1:dNF, 1:dNF+2), A_RR(1:Fdof, 1:Fdof+2), &
            dNF, Fdof, F_size, itrans(i, 1: F_size, 1:2) )

       !call WriteArray(A, dNF, dNF+2, 1, dNF, 1, dNF+2)

       !print*,'-------A2',ip
         
       call AssembMatrixRTN_DG(B(1:dNF+2, 1:dNR), A_RD(1:Fdof+2, 1:Rdof), &
            dNF, dNR, Fdof, Rdof, F_size, itrans(i, 1: F_size, 1:2), i, N, inner)


       !call WriteArray(B, dNF+1, dNR, 1, dNF+1, 1, dNR)


       if(.not. innerR) then 
          ! setting of nonhomogeneous Dirichlet BC
          !if(ip == iptest) print*,'###',ip,elem%i,i,supp(i,1:2)
          if(i == 1 ) then
             !call Eval_Dir_BC(elem, deg, F_face, ie, .true., uD(1:F_face, 1) )
             call Eval_Dir_BC_L2proj(elem, deg, F_face, ie, .true., uD(1:F_face, 1) )
             !if(ip == 1) write(*,'(a10,i5,8es12.4)') 'proj 1',i,uD(1:F_face, 1)

             ! adding to the RHS
             ie0 = ie
             do k=1,F_face
                kk = (ie0 -1) * F_face + k
                !print*,'i=1',k, kk
                A_RR(1:Fdof, Fdof+2) = A_RR(1:Fdof, Fdof+2) - uD(k, 1) *  A_RR(1:Fdof, kk )
                A_RD(Fdof+2, 1:Rdof ) = A_RD(Fdof+2, 1:Rdof ) - uD(k,1) * A_RD(kk, 1:Rdof) 
             enddo
          endif

          if(i == N ) then
             !call Eval_Dir_BC(elem, deg, F_face, ie, .false., uD(1:F_face, 2) )
             call Eval_Dir_BC_L2proj(elem, deg, F_face, ie, .false., uD(1:F_face, 2) )
             !if(ip == 1) write(*,'(a10,i5,8es12.4)') 'proj N',i,uD(1:F_face, 2)

             !stop

             ! adding to the RHS
             ie1 = mod(ie , 3) + 1
             ie0 = mod(ie1, 3) + 1
             do k=1,F_face
                kk = (ie0 -1) * F_face + k
                !print*,'i=N',k, kk
                A_RR(1:Fdof, Fdof+2) = A_RR(1:Fdof, Fdof+2) - uD(k, 2) *  A_RR(1:Fdof, kk )
                A_RD(Fdof+2, 1:Rdof ) = A_RD(Fdof+2, 1:Rdof ) - uD(k,2) * A_RD(kk, 1:Rdof) 
             enddo

          endif


          !print*,'-------A3',ip
          call AssembMatrixRTN_RTN(A_pt(1:dNF_pt, 1:dNF_pt+2), A_RR(1:Fdof, 1:Fdof+2), &
               dNF_pt, Fdof, F_size, itrans(i, 1: F_size, 3:4) )
         
          !call WriteArray(A_pt, dNF_pt, dNF_pt+2, 1, dNF_pt, 1, dNF_pt+2)

          !print*,'-------A5',ip

          call AssembMatrixRTN_DG(B_pt(1:dNF_pt+2, 1:dNR_pt), A_RD(1:Fdof+2, 1:Rdof), &
               dNF_pt, dNR_pt, Fdof, Rdof, F_size, itrans(i, 1: F_size, 3:4), i, N, .true.)


          !call WriteArray(B_pt, dNF_pt+1, dNR_pt, 1, dNF_pt+1, 1, dNR_pt)

       endif

       !!if(ip == 1) 
          

       !stop 'stopped ERTGVBNHG'

       !  allocate( Fx(1:Qdof, 1:nbDim) )
       ! !integration nodes on K
       !  call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof,1:nbDim), &
       !       Fx(1:Qdof, 1:nbDim) )
       
       !  do j=1,Qdof
       !     write(70+i,*) V_rule%lambda(j, 1:nbDim), Fx(j,1:2),DGphi(i,1:Rdof,j)
       !     write(100+i,*)Fx(j,1:2), RTNphi(i,1:Fdof,1,j)
       !     write(200+i,*)Fx(j,1:2), RTNphi(i,1:Fdof,2,j)
       !     write(300+i,*)Fx(j,1:2), RTNphi(i,1:Fdof,3,j)
       !     write(400+i,*)Fx(j,1:2), tau_g(1:3,j)
       !  enddo
       
       !  deallocate(Fx)
       
       
    enddo   ! do i=1,N

       
    !if(.not. innerR ) stop 'Stopped in WQRTY'
    !print*,'#############################################################'
    !if(ip == 15 ) stop 'Stopped in WQRTY'


    !call WriteArray(A, dNF, dNF+1, 1, dNF, 1, dNF+1)
    
    ! if(ip == 15) call WriteArray(B, dNF+1, dNR, 1, dNF+1, 1, dNR)
    
    rhs(:,:) = 0.
    ! setting of RHS for flux
    rhs(    1 : dNF     , 1 ) =  A(1:dNF, dNF + 1 ) 
    rhs(dnF+1 : dNF +dNR, 1 ) =  B(dNF + 1, 1:dNR ) 

    ! setting of RHS for potential
    rhs(    1 : dNF     , 2 ) =  A(1:dNF, dNF + 2 ) 
    rhs(dnF+1 : dNF +dNR, 2 ) =  B(dNF + 2, 1:dNR ) 
    !write(*,'(a6,40es12.4)') 'rhsB', rhs(dnF+1 : dNF +dNR, 2 ) 
    
    !write(201,*) 'Schur:', ip, innerR, 2, grid%x(ip, 1:2), dNF, dNR
    !call SchurComplements(dNF, dNR, A(1:dNF, 1:dNF), -B(1:dNF, 1:dNR), 2, &
    !     rhs(1:dNF, 1:2), rhs(dNF+1:dNF+dNR, 1:2), &
    !     x(1:dNF, 1:2), x(dNF+1:dNF+dNR, 1:2)  ) 


    !call SchurComplements( &
    call SchurComplementsNEW( state%space%adapt%adapt_level, ip,  N, deg, grid%x(ip,1:2), &
    !call SchurComplements_ITER( &
         dNF, dNR, A(1:dNF, 1:dNF), -B(1:dNF, 1:dNR), 2, &
         rhs(1:dNF, 1:2), rhs(dNF+1:dNF+dNR, 1:2), &
         x(1:dNF, 1:2), x(dNF+1:dNF+dNR, 1:2)  ) 


    !write(*,'(a6,50es11.3)') 'A x1:', x(1:dN,1)
    !write(*,'(a6,50es11.3)') 'A x2:', x(1:dN,2)

    ! potential reconstruction for potential, different algebraic problem
    if(.not. innerR) then
       rhs_pt(:,:) = 0.
       ! setting of RHS for flux no necessary
       !rhs_pt(       1 : dNF_pt        , 1 ) =  A_pt(  1:dNF_pt, dNF_pt + 1 ) 
       !rhs_pt(dnF_pt+1 : dNF_pt +dNR_pt, 1 ) =  B_pt(dNF_pt + 1, 1:dNR_pt ) 
       
       ! setting of RHS for potential
       rhs_pt(   1 : dNF_pt     , 2 ) =  A_pt(1:dNF_pt, dNF_pt + 2 ) 
       rhs_pt(dnF_pt+1 : dNF_pt +dNR_pt, 2 ) =  B_pt(dNF_pt + 2, 1:dNR_pt ) 
       !write(*,'(a6,40es12.4)') 'rhsB', rhs_pt(dnF_pt+1 : dNF_pt +dNR_pt, 2 ) 

       
       !print*,'---------',ip
       !call WriteArray(A_pt, dNF_pt, dNF_pt+2, 1, dNF_pt, 1, dNF_pt+2)
       !call WriteArray(B_pt, dNF_pt+1, dNR_pt, 1, dNF_pt+1, 1, dNR_pt)

       !write(201,*) 'Schur:', ip, innerR, 1, grid%x(ip, 1:2), dNF_pt, dNR_pt
       ! call SchurComplements(dNF_pt, dNR_pt, A_pt(1:dNF_pt, 1:dNF_pt), &
       !      -B_pt(1:dNF_pt, 1:dNR_pt), 1, &
       !      rhs_pt(1:dNF_pt, 2), rhs_pt(dNF_pt+1:dNF_pt+dNR_pt, 2), &
       !      x_pt(1:dNF_pt, 2), x_pt(dNF_pt+1 : dNF_pt+dNR_pt, 2)  ) 
       
       !call SchurComplements( &
       call SchurComplementsNEW( state%space%adapt%adapt_level, -ip, N, deg, grid%x(ip,1:2), &
       !call SchurComplements_ITER( &
            dNF_pt, dNR_pt, A_pt(1:dNF_pt, 1:dNF_pt), &
            -B_pt(1:dNF_pt, 1:dNR_pt), 1, &
            rhs_pt(1:dNF_pt, 2), rhs_pt(dNF_pt+1:dNF_pt+dNR_pt, 2), &
            x_pt(1:dNF_pt, 2), x_pt(dNF_pt+1 : dNF_pt+dNR_pt, 2)  ) 

       !write(*,'(a6,50es11.3)') 'B x1:', x_pt(1:dN_pt,1)
       !write(*,'(a6,50es11.3)') 'B x2:', x_pt(1:dN_pt,2)

    endif
    
    !if(ip == 15) then
    !print*,'**********************************   AFTER SCHUR'
    !  if(ip == 15) call WriteArray(B, dNF+1, dNR, 1, dNF+1, 1, dNR)
    !    print*
    !    write(*,'(a4,50es11.3)') 'Bx2:', &
    !         matmul(transpose(B(1:dNF, 1:dNR)), x(1:dNF,2))
    ! endif

    ! setting of the solution sigma in integ. nodes
    allocate( sigma(1:N, 1:2, 1:Qdof) )
    
    do i=1, N
       elem => grid%elem(supp(i,1))
       ie = supp(i,2)   ! inner index of the vertex

       if(innerR) then
          ! both flux and potential reconstruction
          !write(*,'(a10,10i5)') 'All',i,elem%i,ie
          call AssociateElementReconstruction(elem, N, dNF, Fdof, V_rule, ip, i, ie, deg, 2, &
               F_size, itrans(i, 1: F_size, 1:2), &
               x(1:dNF, 1:2), RTNphi(i, 1:Fdof, 1:3, 1:V_rule%Qdof), F_face )
       else
          ! flux reconstruction
          !if(N==1) write(*,'(a10,10i5)') 'flux',i,elem%i,ie
          call AssociateElementReconstruction(elem, N, dNF, Fdof, V_rule, ip, i, ie, deg, 1, &
               F_size, itrans(i, 1: F_size, 1:2), &
               x(1:dNF, 1), RTNphi(i, 1:Fdof, 1:3, 1:V_rule%Qdof), F_face )

          !if(N==1) write(*,'(a10,10i5)') 'potent',i,elem%i,ie

          ! potential  reconstruction
          if(i == 1 .or. i == N ) then
             ! adding of the non-homogeneous Dirichlet BC
             call AssociateElementReconstruction(elem, N, dNF_pt, Fdof, V_rule, ip, i, ie, deg, -1, &
                  F_size, itrans(i, 1: F_size, 3:4), &
                  x_pt(1:dNF_pt, 2), RTNphi(i, 1:Fdof, 1:3, 1:V_rule%Qdof), F_face, uD(1:F_face, 1:2) )
             ! adding of the non-homogeneous Dirichlet BC  -------------------------^^^^^^^^^^^^^^^^^
          else

             call AssociateElementReconstruction(elem, N, dNF_pt, Fdof, V_rule, ip, i, ie, deg, -1, &
                  F_size, itrans(i, 1: F_size, 3:4), &
                  x_pt(1:dNF_pt, 2), RTNphi(i, 1:Fdof, 1:3, 1:V_rule%Qdof), F_face )
          endif
          
       endif

    enddo
       
    deallocate(itrans)

    deallocate(A_RR, A_RD)
    deallocate(DGphi, RTNphi )
    deallocate(rhs, x)
    deallocate(A, B)
    deallocate(tau_g, sigma)

    if(.not. innerR ) then
       deallocate(A_pt, B_pt, rhs_pt, x_pt, uD)
    endif
    
    !if(.not. inner) then
    !if(ip == 3) then 
    !if( inner) then
    !print*,'stopped in LocalNeumannVertexProblem, ip = ',ip
    !   stop 
    !endif

    !endif

  end subroutine LocalNeumannVertexProblem

  !> evaluate the DG functions in integ nodes, generally different with elem%Qdof
  !> for inner = .true. a constrained is imposed
  subroutine  Eval_DGphi_elem(elem,  V_rule, Rdof, DGphi, area, inner )
    class(element), intent(in) :: elem
    type(volume_rule), intent(in) :: V_rule    
    integer, intent(in) :: Rdof
    real, dimension(1:Rdof, 1:V_rule%Qdof), intent(inout) :: DGphi
    real, intent(in) :: area
    logical, intent(in) :: inner
    real, dimension(:), allocatable :: weights 
    real :: val
    integer :: Qdof, i
    
    Qdof = V_rule%Qdof
    
    DGphi(1:Rdof, 1:Qdof) = V_rule%phi(1:Rdof, 1:Qdof)

    if(inner) then
       ! mean value over \omega_a has to be zero
       allocate(weights(1:Qdof) )
       
       call Eval_V_Weights_plus(elem, V_rule, weights(1:Qdof) )

       ! first function is constant = |D|/ 6|K|, |D| = area of \omega_a  
       DGphi(1, 1:Qdof) =  area / elem%F%JF0 / 3  
       
       ! mean values of shape functions i=2,3, .. have to be zero, we shift them
       do i=2, Rdof  
          val = dot_product( DGphi(i, 1:Qdof), weights(1:Qdof) )
          DGphi(i, 1:Qdof) = DGphi(i, 1:Qdof) - val ! shifting
       enddo
    
       deallocate(weights )

    endif

  end subroutine Eval_DGphi_elem

  !> evaluate the function \f$ \tau := \psi_a \nabla u_h \f$ and
  !> \f$ g:=\psi_a f - \nabla u_h \cdot \nabla \psi_a\f$ in integ nodes
  !> generally different with elem%Qdof,
  !> \f$ \psi_a\f$ is the hat function
  subroutine  Eval_tau_g_elem(elem,  i_var, V_rule, tau_g, ie, xi )
    class(element), intent(in) :: elem
    integer, intent(in) :: i_var   ! index of variant
    type(volume_rule), target, intent(in) :: V_rule    
    real, dimension(1:5, 1:V_rule%Qdof), intent(inout) :: tau_g
    integer, intent(in) :: ie  ! inner index of the vertex a (for hat function)
    real, dimension(1:3, 1:2), intent(inout) :: xi
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:,:), allocatable :: Dphi, RTNphi
    real, dimension(:,:), allocatable :: Fx, MM
    real, dimension(:), allocatable :: f, weights
    real, dimension(:,:), allocatable :: uhDPi, uhR
    real :: psi, psiDx, psiDy, uh, uhDx, uhDy
    integer :: Qdof, dof, i, j, ie1, ie2, kst
    integer :: Fdeg,  FFdof, dofM
    
    Qdof = V_rule%Qdof
    dof = elem%dof

    ! test iP2
    if(i_var == 0) dof = dof - (elem%deg + 1)

    phi => V_rule%phi(1:dof, 1:Qdof)

    ! derivative of test functions
    allocate( Dphi(1:dof, 1:nbDim, 1:Qdof) ) 
    call Eval_Dphi_plus(elem, V_rule, dof, Dphi(1:dof, 1:nbDim, 1:Qdof) ) 
 
    ! physical coordinates of integ nodes for evaluating of souce term
    allocate( Fx(1:Qdof, 1:nbDim), f(1:ndim)  )
    call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof,1:nbDim), &
         Fx(1:Qdof, 1:nbDim) )
 
    ! k =1, only scalar case index of component of Dw
    !kst = dof*(k-1) + 1

    ie1 = mod(  ie, 3) + 1
    ie2 = mod( ie1, 3) + 1

    ! gradient of psi_a (constant on elem)
    psiDx = (xi(ie1, 2) - xi(ie2, 2))/elem%area/2
    psiDy = (xi(ie2, 1) - xi(ie1, 1))/elem%area/2
       
    ! test iP2
    if(i_var == -10) then  ! THIS VARIANT SHOULD be carrected with respect this one in SRC/
       allocate(uhDPi(1:6, 1:Qdof) )

       if(elem%i == 1) write(*,'(a6,40es12.4)') '%w',elem%w(0, 1:dof)
       kst = 1
       do i=1,Qdof
          ! value of \psi_a at integ node
          psi = V_rule%lambda(i, ie2)
          !psi = 1
          
          ! approximate solution
          uhDPi(1, i) = dot_product(elem%w(0,kst:kst+dof-1), phi(1:dof, i))
          
          ! gradient of approximate solution
          uhDPi(2, i) = dot_product(elem%w(0,kst:kst+dof-1), Dphi(1:dof, 1, i))
          uhDPi(3, i) = dot_product(elem%w(0,kst:kst+dof-1), Dphi(1:dof, 2, i))
       enddo

       if(elem%i == 1) write(*,'(a6,40es12.4)') 'wi',uhDPi(1, 1:Qdof)
       
       ! projection to the lower space
       dofM = dof - (elem%deg + 1)
       allocate(uhR(1:6, 1:dofM) )

       ! integration with the basis functions 
       uhR(1, 1:dofM) = matmul( phi(1:dofM, 1:Qdof), uhDPi(1,1:Qdof) )
       uhR(2, 1:dofM) = matmul( phi(1:dofM, 1:Qdof), uhDPi(2,1:Qdof) )
       uhR(3, 1:dofM) = matmul( phi(1:dofM, 1:Qdof), uhDPi(3,1:Qdof) )
       
       ! evaluation of u and its gradient in bases functions
       uhR(4, 1:dofM) = matmul( elem%MassInv%Mb(1:dofM, 1:dofM), uhR(1,1:dofM) )
       uhR(5, 1:dofM) = matmul( elem%MassInv%Mb(1:dofM, 1:dofM), uhR(2,1:dofM) )
       uhR(6, 1:dofM) = matmul( elem%MassInv%Mb(1:dofM, 1:dofM), uhR(3,1:dofM) )

       if(elem%i == 1) write(*,'(a6,40es12.4)') 'Pi w_h',uhR(1, 1:dofM)

       ! evaluation in integ nodes
       kst = 1
       do i=1,Qdof
          ! value of \psi_a at integ node
          psi = V_rule%lambda(i, ie2)
          !psi = 1
          
          ! approximate solution
          uh = dot_product(uhR(4, 1:dofM), phi(1:dofM, i))
          
          ! gradient of approximate solution
          uhDx = dot_product(uhR(5, 1:dofM), phi(1:dofM, i))
          uhDy = dot_product(uhR(6, 1:dofM), phi(1:dofM, i))

          ! potential reconstruction 
          ! tau_g(4:5,:)= R_{\pi/2) \nabla(\psi^a u_h) 
          tau_g(4, i) = -(psiDy * uh + psi * uhDy)
          tau_g(5, i) =   psiDx * uh + psi * uhDx
          
          ! flux reconstruction
          ! discrete gradient with the  lifting operator for the dicrete gradient 
          ! in pNeu for SIPG and NIPG [Ern, Vohralik, SINUM 15]
          uhDx = uhDx + state%space%m_IPG*elem%lifting(1, 1)   ! opposite sign !
          uhDy = uhDy + state%space%m_IPG*elem%lifting(1, 2)
          
          tau_g(1, i) = psi * uhDx
          tau_g(2, i) = psi * uhDy
          
          call RHS_Scalar(Fx(i, 1:2), f(1:ndim), state%time%ctime)
          
          tau_g(3, i) =  psi * f(1) - uhDx * psiDx - uhDy * psiDy
       enddo

       deallocate(uhDPi, uhR)
    else
       kst = 1
       do i=1,Qdof
          ! value of \psi_a at integ node
          psi = V_rule%lambda(i, ie2)
          !psi = 1
          
          ! approximate solution
          uh = dot_product(elem%w(0,kst:kst+dof-1), phi(1:dof, i))
          
          ! gradient of approximate solution
          uhDx = dot_product(elem%w(0,kst:kst+dof-1), Dphi(1:dof, 1, i))
          uhDy = dot_product(elem%w(0,kst:kst+dof-1), Dphi(1:dof, 2, i))
          
          ! potential reconstruction 
          ! tau_g(4:5,:)= R_{\pi/2) \nabla(\psi^a u_h) 
          tau_g(4, i) = -(psiDy * uh + psi * uhDy)
          tau_g(5, i) =   psiDx * uh + psi * uhDx
          
          ! flux reconstruction
          ! discrete gradient with the  lifting operator for the dicrete gradient 
          ! in pNeu for SIPG and NIPG [Ern, Vohralik, SINUM 15]
          uhDx = uhDx + state%space%m_IPG*elem%lifting(1, 1)   ! opposite sign !
          uhDy = uhDy + state%space%m_IPG*elem%lifting(1, 2)
          
          tau_g(1, i) = psi * uhDx
          tau_g(2, i) = psi * uhDy
          
          call RHS_Scalar(Fx(i, 1:2), f(1:ndim), state%time%ctime)
          
          tau_g(3, i) =  psi * f(1) - uhDx * psiDx - uhDy * psiDy
          
       enddo

    endif

    ! projection of tau_g(1:2, :) and tau_g(4:5, :) to RTN_{p-1}
    if(i_var == 10) then 

       !write(*,'(a6,300es12.4)') ' tau_h:', tau_g(1, 1:Qdof)
       !write(*,'(a6,300es12.4)') ' tau_h:', tau_g(2, 1:Qdof)
       !write(*,'(x)') 

       allocate(weights(1:Qdof) )
       call Eval_V_Weights_plus(elem, V_rule, weights(1:Qdof) )

       Fdeg = elem%deg !- 1
       FFdof = SetRTNdof( Fdeg )  

       ! RTN basis function in integ nodes
       allocate( RTNphi(1:FFdof, 1:3, 1:Qdof) ) 

       ! mass matrix for the projection
       allocate( MM(1:FFdof, 1:FFdof+2) ) 

       call Eval_RTNphi_elem(elem, V_rule, Fdeg, FFdof, RTNphi(1:FFdof, 1:3, 1:Qdof), 0  )

       !do i=1,Qdof
       !   write(97,*) V_rule%lambda(i,1:2),  RTNphi(1:FFdof, 1,i)
       !   write(98,*) V_rule%lambda(i,1:2),  RTNphi(1:FFdof, 2,i)
       !enddo
    
       do i=1,FFdof
          do j=i, FFdof
             MM(i,j) = dot_product( weights(1:Qdof), &
                  RTNphi(i, 1, 1:Qdof)*RTNphi(j, 1, 1:Qdof) &
                  + RTNphi(i, 2, 1:Qdof)*RTNphi(j, 2, 1:Qdof) )

             MM(j,i) = MM(i,j)
          enddo
          MM(i,FFdof+1) = dot_product( weights(1:Qdof), &
               RTNphi(i, 1, 1:Qdof)* tau_g(1, 1:Qdof) + RTNphi(i, 2, 1:Qdof)* tau_g(2, 1:Qdof) )

          MM(i,FFdof+2) = dot_product( weights(1:Qdof), &
               RTNphi(i, 1, 1:Qdof)* tau_g(4, 1:Qdof) + RTNphi(i, 2, 1:Qdof)* tau_g(5, 1:Qdof) )
       enddo

       ! coefficients of the projection
       call SolveLocalMatrixProblem(FFdof, MM(1:FFdof, 1:FFdof), 2, MM(1:FFdof, FFdof+1 : FFdof+2) ) 


       ! projection in integ nodes
       tau_g(1, 1:Qdof) = matmul( MM(1:FFdof, FFdof+1),  RTNphi(1:FFdof, 1, 1:Qdof) )
       tau_g(2, 1:Qdof) = matmul( MM(1:FFdof, FFdof+1),  RTNphi(1:FFdof, 2, 1:Qdof) )

       tau_g(4, 1:Qdof) = matmul( MM(1:FFdof, FFdof+2),  RTNphi(1:FFdof, 1, 1:Qdof) )
       tau_g(5, 1:Qdof) = matmul( MM(1:FFdof, FFdof+2),  RTNphi(1:FFdof, 2, 1:Qdof) )

       
       !write(*,'(a6,300es12.4)') 'Pi tau_h:', tau_g(1, 1:Qdof)
       !write(*,'(a6,300es12.4)') 'Pi tau_h:', tau_g(2, 1:Qdof)
       !write(*,*) '---------------------------------------------'
       !stop
       deallocate(weights, MM, RTNphi)
    endif ! if(i_var == 0) 

       !if(elem%i == 1) then
    !    write(*,'(a6,2i5,40es12.4)') 'g^a:',elem%i,ie,psiDx, psiDy,  tau_g(3, 1:Qdof)
    ! endif
    deallocate(Dphi, Fx)
 

  end subroutine Eval_tau_g_elem


  !> evaluate the RTN functions (2nd idx =1,2) and its divergence (2nd idx = 3)
  !> in integ nodes,  generally different with elem%Qdof
  subroutine  Eval_RTNphi_elem(elem, V_rule, Fdeg, Fdof, RTNphi, ie )
    class(element), intent(in) :: elem
    type(volume_rule), intent(in) :: V_rule
    integer, intent(in) :: Fdeg, Fdof
    real, dimension(1:Fdof, 1:3, 1:V_rule%Qdof), intent(inout):: RTNphi
    integer, intent(in) :: ie  ! inner index of the vertex
    type(basis_rtn_fe), pointer:: loc_RTN
    real, dimension(:,:), allocatable :: MMelem, Fx
    real, dimension(:,:,:), allocatable :: psi
    real :: val, val_sum
    integer :: Qdof, i, ii, j, k, FF
    
    Qdof = V_rule%Qdof
    
    loc_RTN => state%loc_RTN(Fdeg)
    if(.not. loc_RTN%defined ) call Init_Loc_RTN(state%loc_RTN(Fdeg), Fdeg )

    allocate(MMelem(1:Fdof, 1:Fdof) )

    ! evaluate the momentums of the local RTN basis on K
    call ComputeLocRTNMomentumsElem2(elem, V_rule, Fdeg, Fdof, MMelem) 

     ! do i=1,Fdof
     !    write(*,'(i5,20es12.4)') i,MMelem(i,:)
     !enddo
     ! print*,'---------- orig  ----------'

    MMelem(1:Fdof, 1:Fdof) = transpose( MMelem(1:Fdof, 1:Fdof) )

      ! do i=1,Fdof
      !    write(*,'(i5,200es12.4)') i,MMelem(i,:)
      ! enddo
      ! print*,'---------- transpose  ----------'

    ! momentums of the element RTN functions form identical matrix
    call MblockInverse(Fdof, MMelem )
    !call MblockInverseVAR(Fdof, MMelem )
    
    
    ! do i=1,Fdof
    !     write(*,'(i5,200es12.4)') i,MMelem(i,:)
    !  enddo
    !  print*,'---------- inverse  ----------'
    ! stop

    ! on the ie-th edge we need the value -1, opposite normal
    if(ie > 0) then
       do i=1,Fdeg+1
          ii = (ie-1)*(Fdeg+1) + i
          MMelem(ii, 1:Fdof) = - MMelem(ii, 1:Fdof)
          !print*,'### neumann.f90 #######', elem%i, i, ii
       enddo
    endif

     ! do i=1,Fdof
     !    write(*,'(i5,20es12.4)') i,MMelem(i,:)
     ! enddo
     ! print*,'---------- inverse  ----------'

    !  RTN basis functions on elem in integ nodes
    allocate(psi(1:Fdof, 1:3, 1:Qdof) )
    call Eval_Loc_RTN_ElemGE(elem, V_rule, loc_RTN, psi(1:Fdof, 1:3, 1:Qdof) )

    do i=1,Qdof
       RTNphi(1:Fdof, 1, i) = matmul(MMelem(1:Fdof, 1:Fdof), psi(1:Fdof, 1, i) )
       RTNphi(1:Fdof, 2, i) = matmul(MMelem(1:Fdof, 1:Fdof), psi(1:Fdof, 2, i) )
       RTNphi(1:Fdof, 3, i) = matmul(MMelem(1:Fdof, 1:Fdof), psi(1:Fdof, 3, i) )
    enddo

    deallocate(MMelem, psi)

    ! do i=1, Fdof
    !    write(*,'(2i5,420es12.4)') i, 1, maxval(abs(RTNphi(i, 1, 1:Qdof))), RTNphi(i, 1, 1:Qdof)
    !    write(*,'(2i5,420es12.4)') i, 2, maxval(abs(RTNphi(i, 2, 1:Qdof))), RTNphi(i, 2, 1:Qdof)
    !    write(*,'(2i5,420es12.4)') i,-3, maxval(abs(RTNphi(i, 3, 1:Qdof))), RTNphi(i, 3, 1:Qdof)
    !    print*
    ! enddo
    ! stop

    ! physical coordinates of integ nodes for evaluating of souce term
    ! allocate( Fx(1:Qdof, 1:nbDim),  )
    ! call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof,1:nbDim), &
    !      Fx(1:Qdof, 1:nbDim) )

    ! if(ie == 1 .and. Fdof == 8) then
    !    do i=1,Qdof
    !       write(41,*) Fx(i, 1:2),V_rule%lambda(i,1:nbDim),'|', RTNphi(1:Fdof, 1, i)
    !       write(42,*) Fx(i, 1:2),V_rule%lambda(i,1:nbDim),'|', RTNphi(1:Fdof, 2, i)
    !       write(43,*) Fx(i, 1:2),V_rule%lambda(i,1:nbDim),'|', RTNphi(1:Fdof, 3, i)
    !    enddo
    ! endif

    ! deallocate(Fx)
    
    !return

    ! orthognalization of the "volume" RTN functions
    allocate(psi(1:Fdof, 1:3, 1:Qdof) )
    psi(1:Fdof, 1:3, 1:Qdof) = RTNphi(1:Fdof, 1:3, 1:Qdof)
    FF = (Fdeg+1)*3 

    !do i= (Fdeg+1)*3 +1, Fdof
    !   write(*,'(a4,i5,200es12.4)') 'DWE',i, psi(i, 1, 1:8)
    !enddo

    
    !print*,'matrix bofore:            ................'
    !do i= (Fdeg+1)*3 +1, Fdof
    !   do k = (Fdeg+1)*3 +1, Fdof
    !      val = dot_product(V_rule%weights(1:Qdof),  &
    !           psi(i, 1, 1:Qdof)*psi(k, 1, 1:Qdof) + psi(i, 2, 1:Qdof)*psi(k, 2, 1:Qdof) )
    !      write(*,'(a4,2i5,200es12.4)') 'DWE',i, k, val
    !   enddo
    !enddo

    do k = FF +1, Fdof
       do i=  FF +1, k-1
          val = dot_product(V_rule%weights(1:Qdof),  &
               psi(i, 1, 1:Qdof)*psi(k, 1, 1:Qdof) + psi(i, 2, 1:Qdof)*psi(k, 2, 1:Qdof) )
          
          psi(k, 1:3, 1:Qdof) =  psi(k, 1:3, 1:Qdof) - val * psi(i, 1:3, 1:Qdof)
          !print*,'d...fres',i, k, val
       enddo
       
       val = dot_product(V_rule%weights(1:Qdof),  &
            psi(k, 1, 1:Qdof)*psi(k, 1, 1:Qdof) + psi(k, 2, 1:Qdof)*psi(k, 2, 1:Qdof) )

       psi(k, 1:3, 1:Qdof)  = psi(k, 1:3, 1:Qdof)  / val**0.5

       
       !print*,'ddedfres',i, val
    enddo

    !do i= (Fdeg+1)*3 +1, Fdof
    !   write(*,'(a4,i5,200es12.4)') 'DWE',i, psi(i, 1, 1:8)
    !enddo

     ! print*,'matrix after:            ................'
     ! do i= (Fdeg+1)*3 +1, Fdof
     !    do k = (Fdeg+1)*3 +1, Fdof
     !       val = dot_product(V_rule%weights(1:Qdof),  &
     !            psi(i, 1, 1:Qdof)*psi(k, 1, 1:Qdof) + psi(i, 2, 1:Qdof)*psi(k, 2, 1:Qdof) )
     !       write(*,'(a4,2i5,200es12.4)') 'DWE',i, k, val
     !    enddo
     ! enddo

    !!VD 
    RTNphi(FF+1:Fdof, 1:3, 1:Qdof) = psi(FF+1:Fdof, 1:3, 1:Qdof) !* elem%area**0.5

    !return

    ! a partial orogonalization of the edge momentum functions
    do k=1, FF
       do i=FF+1, Fdof
          val = dot_product(V_rule%weights(1:Qdof),  &
               psi(i, 1, 1:Qdof)*psi(k, 1, 1:Qdof) + psi(i, 2, 1:Qdof)*psi(k, 2, 1:Qdof) )
          
          psi(k, 1:3, 1:Qdof) =  psi(k, 1:3, 1:Qdof) - val * psi(i, 1:3, 1:Qdof)
       enddo
    enddo
    
    RTNphi(1:FF, 1:3, 1:Qdof) = psi(1:FF, 1:3, 1:Qdof)
    !stop

    ! "orthoNORMALIZATION"
    val_sum = 0.
    do k=1, FF
       val = dot_product(V_rule%weights(1:Qdof),  &
            psi(k, 1, 1:Qdof)*psi(k, 1, 1:Qdof) + psi(k, 2, 1:Qdof)*psi(k, 2, 1:Qdof) )
       val_sum = val_sum + val
    enddo
    val_sum = val_sum / FF

    do k= FF + 1 , Fdof
       val = dot_product(V_rule%weights(1:Qdof),  &
            psi(k, 1, 1:Qdof)*psi(k, 1, 1:Qdof) + psi(k, 2, 1:Qdof)*psi(k, 2, 1:Qdof) )

       RTNphi(k, 1:3, 1:Qdof) = psi(k, 1:3, 1:Qdof) /val *  val_sum 
    enddo

    deallocate(psi)
    
  end subroutine Eval_RTNphi_elem

  

  !> integrate the dot_product of the RTN functions over elem
  !> last colom of A_RR is the right-hand side
  subroutine  Integ_dot_product_RTNphi_elem(elem, V_rule, Fdof, RTNphi, tau_g, &
       A_RR)
    class(element), intent(in) :: elem
    type(volume_rule), intent(in) :: V_rule
    integer, intent(in) :: Fdof
    real, dimension(1:Fdof, 1:2, 1:V_rule%Qdof), intent(in):: RTNphi
    real, dimension(1:5, 1:V_rule%Qdof), intent(in):: tau_g ! RHS
    real, dimension(1:Fdof,1:Fdof+2), intent(inout) :: A_RR
    real, dimension(:), allocatable :: weights
    integer :: Qdof, i,j
    
    Qdof = V_rule%Qdof

        ! mean value over \omega_a has to be zero
    allocate(weights(1:Qdof) )

    if(elem%F%iFlin) then ! linear element, constant Jacobian
       weights(1:Qdof)  = V_rule%weights(1:Qdof) * elem%F%JF0 / 2
    else
       print*,'Curved elements in  Eval_DGphi_elem not implemented'
    endif

    do i=1,Fdof
       do j=i,Fdof
          A_RR(i,j) = dot_product( weights(1:Qdof), &
               RTNphi(i, 1, 1:Qdof) * RTNphi(j, 1, 1:Qdof) &
               + RTNphi(i, 2, 1:Qdof) * RTNphi(j, 2, 1:Qdof) )
          A_RR(j,i) = A_RR(i,j)

       enddo
       ! RHS for flux reconstruction
       A_RR(i, Fdof+1) =  - dot_product( weights(1:Qdof), &
            tau_g(1, 1:Qdof) * RTNphi(i, 1, 1:Qdof) &
            +tau_g(2, 1:Qdof) * RTNphi(i, 2, 1:Qdof) )

       ! RHS for potential reconstruction
       A_RR(i, Fdof+2) =  - dot_product( weights(1:Qdof), &
            tau_g(4, 1:Qdof) * RTNphi(i, 1, 1:Qdof) &
            +tau_g(5, 1:Qdof) * RTNphi(i, 2, 1:Qdof) )
    enddo


    deallocate(weights)

  end subroutine Integ_dot_product_RTNphi_elem

  
  !> integrate the product of div RTN functions with DG functions over elem
  !> last line is the RHS
  subroutine Integ_divRTNphi_DGphi_elem(elem, V_rule, Fdof, Rdof, divRTNphi, &
       DGphi, tau_g, A_RD )
    class(element), intent(in) :: elem
    type(volume_rule), intent(in) :: V_rule
    integer, intent(in) :: Fdof, Rdof
    real, dimension(1:Fdof, 1:V_rule%Qdof), intent(in):: divRTNphi
    real, dimension(1:Rdof, 1:V_rule%Qdof), intent(in):: DGphi
    real, dimension(1:V_rule%Qdof), intent(in):: tau_g
    real, dimension(1:Fdof+2, 1:Rdof), intent(inout) :: A_RD
    real, dimension(:), allocatable :: weights
    integer :: Qdof, i,j
    
    Qdof = V_rule%Qdof

    ! mean value over \omega_a has to be zero
    allocate(weights(1:Qdof) )

    if(elem%F%iFlin) then ! linear element, constant Jacobian
       weights(1:Qdof)  = V_rule%weights(1:Qdof) * elem%F%JF0 / 2
    else
       print*,'Curved elements in  Eval_DGphi_elem not implemented'
    endif

    do j=1,Rdof
       do i=1,Fdof
          A_RD(i,j) = dot_product( weights(1:Qdof),  divRTNphi(i, 1:Qdof) * DGphi(j, 1:Qdof)  )
       enddo
       A_RD(Fdof+1,j) = dot_product( weights(1:Qdof),  tau_g(1:Qdof) * DGphi(j, 1:Qdof)  )
       A_RD(Fdof+2,j) = 0.
    enddo

    deallocate(weights)

  end subroutine Integ_divRTNphi_DGphi_elem

  

  !> evaluate the momentums of the basis RTN function of degree Fdeg on elem
  !> and the resulting matrix is stored in MMRE
  !> for the face momentum, the nodes values are used 
  !> (compare ComputeLocRTNMomentumsElem )  better conditionality
  !> no HG nodes
  subroutine ComputeLocRTNMomentumsElem2(elem, V_rule, Fdeg, Fdof, MMRE)
    class(element), intent(in) :: elem      
    type(volume_rule), intent(in) :: V_rule
    integer, intent(in) :: Fdeg, Fdof
    real, dimension(1:Fdof, 1:Fdof), intent(inout) :: MMRE 
    type(basis_rtn_fe), pointer :: loc_RTN
    real, dimension(:,:,:), allocatable :: psi
    real, dimension(:,:), allocatable :: xi, Fxi, qi, phi
    integer :: Qdof, dof
    integer :: ie, ie1, ie2, j, ib, it, i
    integer :: indx, iphi
    real :: rlen

    loc_RTN => state%loc_RTN(Fdeg)
    dof = maxval(loc_RTN%ipsi(:, 2))

    Qdof = 3*(Fdeg + 1) ! number of nodes on faces
    allocate( xi(1: Qdof, 1:3)) ! barycentric coordinates of nodes on faces
    allocate(Fxi(1: Qdof, 1:2)) ! physical  coordinates of nodes on faces
    allocate(phi(1: dof, 1:Qdof)) ! test functions on given nodes
    allocate( psi(1:Fdof, 1:nbDim, 1: Qdof)  )! RTN test functions

    MMRE(:, :) = 0 ! first index moment, second index RTN test function
    
    ! barycentres of 
    xi(:,:) = 0.
    it = 0   

    ! face momentums, setting of the nodes on faces, where $\psi\cdot\nn$ is given
    do ie=1,3 ! loop over triagle edges (including HG nodes )
       ie1 = mod(ie , 3) + 1
       ie2 = mod(ie1, 3) + 1

       
       if(Fdeg == 0) then  ! only one node in face centers
          it = it + 1
          xi(it, ie2) = 0.5
          xi(it, ie) = 1. - xi(it,ie2)
          !write(*,'(a6,4i5, 4es12.4)') '!@#',elem%i,ie, ie1, ie2,  xi(1,:)
       else
          rlen = 1./Fdeg
          do j=0, Fdeg
             it = it + 1
             xi(it, ie2) = 1. - j * rlen
             xi(it, ie) = 1. - xi(it,ie2)
             !write(*,'(a6,4i5, 4es12.4)') '!@#',elem%i,ie, ie1, ie2,  xi(1,:)
          enddo
       endif
    enddo

    !do it=1, Qdof
    !   write(*,'(a6,i5, 7es12.4)') '!@#', it, xi(it,:),Fxi(it, 1:nbDim)
    !   write(40,*) Fxi(it, 1:nbDim),  phi(1:dof,it)
    !enddo

    call ComputeF(elem, Qdof, xi(1:Qdof, 1:nbDim), Fxi(1:Qdof, 1:nbDim) )
    ! "renormalization" of the coordinates
    do i=1,Qdof
       Fxi(i, 1:nbDim) = (Fxi(i, 1:nbDim) - elem%xc(1:nbDim))/elem%diam
    enddo

    call Eval_phi_Qnode(elem, dof, Qdof, xi(1:Qdof,1:2), phi(1:dof, 1:Qdof) )

    ! reference RTN basis functions of elem
    psi(:,:,:) = 0.
    do i=1,Fdof
       indx = loc_RTN%ipsi(i, 1)
       iphi = loc_RTN%ipsi(i, 2)

       if(indx <= nbDim) then
          psi(i, indx, 1:Qdof) = phi(iphi, 1:Qdof)
          
       else
          psi(i, 1, 1:Qdof) = phi(iphi, 1:Qdof) * Fxi(1:Qdof,1)
          psi(i, 2, 1:Qdof) = phi(iphi, 1:Qdof) * Fxi(1:Qdof,2)

       endif
    enddo

    ! MMRE contains the values of \psi \cdot \nn on the face
    ! face momentums
    do ie=1,3 ! loop over triagle edges (including HG nodes )

       ie1 = (ie-1)*(Fdeg + 1 ) + 1   ! initial and final node on edge ie
       ie2 =     ie*(Fdeg + 1 ) 

       do j=ie1, ie2  ! nodes corresponding to the momentums
          MMRE(j, 1:Fdof) = (psi(1:Fdof, 1, j) * elem%n(ie, 1) &
               + psi(1:Fdof, 2, j) * elem%n(ie, 2)) / elem%dn(ie)
       enddo !j
    enddo !ie

    deallocate(xi, Fxi, phi, psi)

    ! VOLUME MOMENTUMS
    !Qdof = elem%Qdof
    ! NV
    Qdof = V_rule%Qdof

    !print*,'#EDE', elem%Qdof, Qdof
    !stop

    allocate( psi(1:Fdof, 1:3, 1: Qdof)  )
    allocate(  qi(1:Qdof, 1:nbDim  ) )


    ! volume integ nodes
    allocate(xi(1:Qdof, 1:nbDim)) 
    !xi(1:Qdof, 1:2) = state%space%V_rule(elem%Qnum)%lambda(1:Qdof,1:2)
    ! NV
    xi(1:Qdof, 1:2) = V_rule%lambda(1:Qdof,1:2)

    !call Eval_Loc_RTN_Elem(elem, loc_RTN, psi(1:Fdof, 1:nbDim, 1:Qdof) )
    ! NV
    call Eval_Loc_RTN_ElemGE(elem, V_rule, loc_RTN, psi(1:Fdof, 1:3, 1:Qdof) )

    ib = 3*(Fdeg + 1)
    do j=1, Fdeg*(Fdeg+1)  ! Fdeg*(Fdeg+1) degrees of freedom over element 
       ib = ib + 1

       call EvalMomentumVolumeD(Fdeg, ib, Qdof, xi(1:Qdof, 1:nbDim), &
            qi(1:Qdof, 1:nbDim))

       ! write(*,'(2i5,420es12.4)') ib, 1, maxval(abs(qi(1:Qdof, 1))), qi(1:Qdof, 1)
       ! write(*,'(2i5,420es12.4)') ib, 2, maxval(abs(qi(1:Qdof, 2))), qi(1:Qdof, 2)
       ! print*

       ! normalization w.r.t. area of element
       !qi(1:Qdof, 1:2) = qi(1:Qdof, 1:2) / elem%area
       !!! VD
       qi(1:Qdof, 1:2) = qi(1:Qdof, 1:2) / elem%area  / elem%diam !!**2

       do it = 1, Fdof
          !call IntegrateFunctionsVec(elem,  psi(it, 1:nbDim, 1:Qdof), &
          !     qi(1:Qdof, 1:nbDim),  MMRE(ib, it))
          ! NV
          call IntegrateFunctionsVec_plus(elem, V_rule, psi(it, 1:nbDim, 1:Qdof), &
               qi(1:Qdof, 1:nbDim),  MMRE(ib, it))

       enddo !it

       !write(*,'(2i5,420es12.4)') ib, 1, maxval(abs(MMRE(ib, :))), MMRE(ib, :)

    enddo !j 

    !!! MMRE is not inverse !!!!
    !call MblockInverse(Fdof, MMRE )

    !do j=1, Fdof
    !   write(*,'(i5,2es12.4,a3,400es12.4)') &
    !        j, maxval(abs(MMRE(j, 1:Fdof))), sum(abs(MMRE(j, 1:Fdof)))/Fdof, '|', &
    !        MMRE(j, 1:Fdof)
    !   if(j == 3*(Fdeg+1) ) print*
    !enddo
    !stop
    deallocate(xi, qi, psi)

    !print*,'  end subroutine ComputeRTNMomentumsRealElem'
    
  end subroutine ComputeLocRTNMomentumsElem2

  !> assamble the global matrix over the whole patch \f$ \omega_a \f$
  !> with the local matrix computed on \f$ K \in \omega_a\f, 9 blocks
  !subroutine AssembMatrixRTN_RTN(A, A_RR, dNF, Fdof, i, N, ie, deg)
  subroutine AssembMatrixRTN_RTN(A, A_RR, dNF, Fdof, F_size, itrans)
    real, dimension(1:dNF, 1:dNF+2), intent(inout) :: A
    real, dimension(1:Fdof, 1:Fdof+2), intent(inout) :: A_RR
    integer, intent(in) :: dNF, Fdof !, i, N, ie, deg
    integer, intent(in) :: F_size   ! number of DTN DOF on one element
    integer, dimension(1:F_size, 1:2), intent(inout) :: itrans
    integer :: j, k !ie1, ie2, i_n, j,k

    do j=1, F_size

       if(itrans(j,2) > 0) then  ! is negative for potential reconstr.n at bound. node

          do k=1, F_size
             if(itrans(k,2) > 0) &
                  A(itrans(j,2), itrans(k,2) ) = A(itrans(j,2), itrans(k,2) ) &
                  + A_RR(itrans(j,1), itrans(k,1) ) 

             !!!print*,'WESA',j,k,itrans(j,1), itrans(k,1)
          enddo
          ! RHS for flux
          A(itrans(j,2), dNF+1 ) = A(itrans(j,2), dNF + 1 ) &
               + A_RR(itrans(j,1), Fdof+1 ) 
          
          ! RHS for potential
          A(itrans(j,2), dNF+2 ) = A(itrans(j,2), dNF + 2 ) &
               + A_RR(itrans(j,1), Fdof+2 ) 
       endif
    enddo

    !deallocate(itrans)

  end subroutine AssembMatrixRTN_RTN


  !> assamble the global matrix over the whole patch \f$ \omega_a \f$
  !> with the local matrix computed on \f$ K \in \omega_a\f, 9 blocks
  !subroutine AssembMatrixRTN_DG(B, A_RD, dNF, dNR, Fdof, Rdof, i, N, ie, deg)
  !> if (DGconstraint) we require that \f$ (q_h,1)_{\omega_a}  = 0 \forall q_h\in Q_h^a\f$
  subroutine AssembMatrixRTN_DG(B, A_RD, dNF, dNR, Fdof, Rdof, F_size, itrans, i, N, DGconstraint)
    real, dimension(1:dNF+2, 1:dNR), intent(inout) :: B
    real, dimension(1:Fdof+2, 1:Fdof), intent(inout) :: A_RD
    integer, intent(in) :: dNF, dNR, Fdof, Rdof,  i, N !, ie, deg
    integer, intent(in) :: F_size   ! number of DTN DOF on one element
    integer, dimension(1:F_size, 1:2), intent(inout) :: itrans
    logical, intent(in) :: DGconstraint
    !integer, dimension(:,:), allocatable :: itrans
    integer :: j, jj, k, kk, ii, ishift, il
    !integer :: F_face, F_vol, F_tot, deg1


    ishift = 0
    if(DGconstraint) ishift = -1  ! due to the constrain 1st DG function is not basis function

    !print*,'DGconstraint', DGconstraint, ishift

    do j=1, F_size
       if(itrans(j,2) > 0) then ! is negative for potentl reconstr.n at bound. node
          do k=1, Rdof  ! local DG index
             kk = (i-1)*Rdof + k +ishift  ! global DG index
             if(kk > 0) then
                B(itrans(j,2), kk ) = B(itrans(j,2), kk) + A_RD(itrans(j,1), k) 
                
                !  write(*,'(a10,2i5,a3,2i5)') & 'transA:',itrans(j,1), k ,'==>',itrans(j,2), kk
             else
                ! first piecewise constant function, has to be substracted
                ! from all triangles
                do ii = 2, N
                   kk = (ii-1)*Rdof 
                   B(itrans(j,2), kk ) = B(itrans(j,2), kk) - A_RD(itrans(j,1), 1) 
                enddo
                
             endif
          enddo
       end if
    enddo

    ! RHS
    do il = 1, 2
       j = Fdof + il
       jj = dNF + il
       do k=1, Rdof  ! local DG index
          kk = (i-1)*Rdof + k + ishift  ! global DG index
          if(kk > 0) then
             B(jj, kk ) = B(jj, kk) + A_RD(j, k) 
             
             !write(*,'(a10,2i5,a3,2i5)')'transB:',j, k ,'==>', jj, kk
          else
             ! first piecewise constant function, has to be substracted
             ! from all triangles
             do ii = 2, N
                kk = (ii-1)*Rdof 
                B(jj, kk ) = B(jj, kk) - A_RD(j, 1) 
             enddo
             
          endif
       enddo
    enddo

  end subroutine AssembMatrixRTN_DG


  !> solution of the local neumann problem is recomputed in integ nodes of each element
  !> and then in basis functions, irhs is number of of RHS:
  !> irhs = 2   => flux and potentials are reconstructed in the same way
  !> irhs = 1   => flux is reconstructed (for boundary nodes) 
  !> irhs = -1  => potential is reconstructed (for boundary nodes) 
  subroutine AssociateElementReconstruction(elem, N, dNF, Fdof, V_rule, ip, i, ie, deg, irhs, &
       F_size, itrans,x, RTNphi, F_face, uD)
    class(element), intent(inout) :: elem
    integer, intent(in) :: N   ! number of elements in the stencil
    integer, intent(in) :: ip  ! number of the vertex (index of the Neumann problem)
    integer, intent(in) :: i   ! index  of the element in the stencil
    integer, intent(in) :: dNF, Fdof, deg, irhs, ie
    type(volume_rule), target, intent(in) :: V_rule
    integer, intent(in) :: F_size   ! number of DTN DOF on one element
    integer, dimension(1:F_size, 1:2), intent(inout) :: itrans
    real, dimension(1:dNF, 1: abs(irhs)), intent(inout) :: x ! solution of lin. alg. problem by Schur
    real, dimension(1:Fdof, 1:3, 1:V_rule%Qdof) , intent(inout) :: RTNphi 
    integer, intent(in):: F_face
    real, dimension(1:F_face,1:2), intent(in), optional :: uD
    !integer, intent(in), optional :: ie0
    real, dimension(:,:,:), allocatable :: flux
    real, dimension(:,:), allocatable :: Fx, Fxi, xi, Fxj
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:,:), allocatable :: Dphi
    real, dimension(:,:), allocatable :: phiT,phiL
    real, dimension(:,:,:), allocatable :: DphiL, DphiLREF
    real, dimension(:,:), allocatable :: M, Mx
    real, dimension(:), allocatable :: weights, pot_coeffs, val
    integer :: ie0, ie1, ie2, i_n, j, k, kk, kk1, jg, jl, Qdof, dof
    integer :: F_vol, F_tot, deg1, deg2, dof2
    integer :: nface, it, j1, j2,l
    type(Lagrang_rule), pointer :: L_rule
    real :: si
    logical :: igraf
    integer:: ig_Fe,  ig_Fv, ig_Pe, ig_Pv = 91

    igraf = .false.
    !igraf = .true.

    !if(elem%i == 1 .and. ip == 15) igraf = .true.
    !if(ip == 1) igraf = .true.
    !if(ip == 2) igraf = .true.
    !if(ip == 12) igraf = .true.
    !if(elem%i == 1 )  igraf = .true.
    !if(elem%i >= 7 .and. elem%i <= 10)  igraf = .true.
    !if(elem%i >= 40 .and. elem%i <= 42)  igraf = .true.
    !if(elem%i >= 80 .and. elem%i <= 82)  igraf = .true.
    !if(elem%i >= 112 .and. elem%i <= 114)  igraf = .true.

    !print*,'EDEWSEEL ELEMENET ', elem%i

    ! weigths for numerical integration
    Qdof = V_rule%Qdof
    allocate(weights(1:Qdof) )

    call Eval_V_Weights_plus(elem, V_rule, weights(1:Qdof) )

    ! polynomial degress of reconstructure
    deg1 = deg + 1
    deg2 = deg + 1  ! 2   NO INFLUENCE since gradient is divergence free

    elem%RTNflux_deg =  max(elem%RTNflux_deg, deg1, deg2)

    if( elem%RTNflux_deg > MaxDegreeImplemented) then
       print*,'Troubles in AssociateElementReconstruction with max deg'
       stop
    endif

    ! setting of the flux reconstruction (solution $x$) into the integ. nodes (global => local)
    allocate(flux(1: abs(irhs)+1, 1:3, 1:Qdof) ) !2nd idx: 1,2 -components, 3 = div

    flux(:,:,:) = 0.

    do k=1, abs(irhs)

       do j=1, F_size
          jg = itrans(j,2)  ! global index
          jl = itrans(j,1)  ! local index

          if(jg > 0) then
               flux(k,1:3, 1:Qdof) = flux(k,1:3, 1:Qdof) +  RTNphi( jl, 1:3, 1:Qdof) * x(jg, k)

               !if(ip == 15 .and. k == 2) then
               !do l=3,3
               !print*,'WED'
               !if(N==1) write(*,'(a6,6i5,2(a2,2es12.4),a2,es12.4)') 'FLUX',elem%RTNflux_deg, k, j,jl,jg
               !write(*,'(a6,6i5,2(a2,2es12.4),a2,es12.4)') 'FLUX',elem%RTNflux_deg, k, j,jl,jg, l,&
               !    '|', flux(k, l, 1:2), '|',  RTNphi( jl, l, 1:2), '|', x(jg, k)
               !enddo
               !endif
            endif
       enddo

       ! adding of the non-homogeneous Dirichlet BC
       if( present(uD)) then
          ! the first face
          if(i==1) then
             ie0 = ie
             do j =1,F_face
                jl = (ie0 -1) * F_face + j
                flux(k,1:3, 1:Qdof) = flux(k,1:3, 1:Qdof) +  RTNphi( jl, 1:3, 1:Qdof) * uD(j,1)
                
                !if(N==1) write(*,'(a6,6i5,2(a2,2es12.4),a2,es12.4)') 'ADFLUX',ie0,k, j,jl
             enddo
          endif

          ! the second face
          if(i==N) then
             ie1 = mod(ie , 3) + 1
             ie0 = mod(ie1, 3) + 1
             do j =1,F_face
                jl = (ie0 -1) * F_face + j
                flux(k,1:3, 1:Qdof) = flux(k,1:3, 1:Qdof) +  RTNphi( jl, 1:3, 1:Qdof) * uD(j,2)
                
                !if(N == 1) write(*,'(a6,6i5,2(a2,2es12.4),a2,es12.4)') 'ADFLUX',ie0,k, j,jl
             enddo
          endif

       endif

       !if(ip == 15 .and. k == 2)      print*
    enddo

    !write(*,'(a6,i5,20es12.4)') 'FLUX',elem%RTNflux_deg, flux(2, 1, 1:10)
    ! write(*,'(a6,i5,20es12.4)') 'FLUX',elem%RTNflux_deg, flux(2, 2, 1:10)
    ! write(*,'(a6,i5,20es12.4)') 'FLUX',elem%RTNflux_deg, flux(2, 3, 1:10)
    ! print* !,Qdof

    ! potential  reconstruction
    ! flux(2, 1:2, 1:Qdof) contains \dzeta^a = -R_{\pi/2} \nabla s_h in integ nodes
    ! we perform -R_{\pi/2} operation
    if(irhs ==2) then
       flux(3, 1:2, 1:Qdof) = flux(2, 1:2, 1:Qdof)
       flux(2, 1, 1:Qdof) = -flux(3, 2, 1:Qdof)
       flux(2, 2, 1:Qdof) =  flux(3, 1, 1:Qdof)
    endif
    ! the same for the potential  reconstruction on boundary nodes
    if(irhs == -1) then
       flux(2, 1:2, 1:Qdof) = flux(1, 1:2, 1:Qdof)
       flux(1, 1, 1:Qdof) = -flux(2, 2, 1:Qdof)
       flux(1, 2, 1:Qdof) =  flux(2, 1, 1:Qdof)
    endif

    !write(*,'(a6,i5,120es12.4)') 'FLUX',elem%RTNflux_deg, flux(2, 1, 1:Qdof)
    ! write(*,'(a6,i5,120es12.4)') 'FLUX',elem%RTNflux_deg, flux(2, 2, 1:Qdof)
    ! write(*,'(a6,i5,120es12.4)') 'FLUX',elem%RTNflux_deg, flux(2, 3, 1:Qdof)
    ! print*,'------------------- AFTER CHANGE '

    !write(*,'(a6,i5,120es12.4)') 'FLUX',elem%RTNflux_deg, flux(1, 1, 1:Qdof)
    !write(*,'(a6,i5,120es12.4)') 'FLUX',elem%RTNflux_deg, flux(1, 2, 1:Qdof)
    !write(*,'(a6,i5,120es12.4)') 'FLUX',elem%RTNflux_deg, flux(1, 3, 1:Qdof)

    ! flux reconstruction in the DG basis functions, degree deg1 = deg + 1
    if(irhs == 2 .or. irhs == 1 .or. irhs == -1) then
       dof =  DOFtriang( deg1 )  
       phi => V_rule%phi(1:dof, 1:Qdof)
       
       allocate(M(1:dof, 1:dof+ 2* 2) ) !irhs) )  ! mass matrix and RHS of reconstructions
       M(:,:) = 0.
       do j =1, dof
          ! mass matrix = (\phi_i, \phi_j)_K
          do k=j, dof
             M(j,k) = dot_product(weights(1:Qdof), phi(j, 1:Qdof) * phi(k, 1:Qdof))
             M(k,j) = M(j,k)
          enddo
          ! RHS = (flux, \phi_i)_K
          do k=1, abs(irhs)
             M(j,dof + 2*k-1) = dot_product(weights(1:Qdof), phi(j, 1:Qdof) * flux(k,1, 1:Qdof) )
             M(j,dof + 2*k  ) = dot_product(weights(1:Qdof), phi(j, 1:Qdof) * flux(k,2, 1:Qdof) )
          enddo
       enddo
       
       ! OLD
       !call MblockInverse(dof, M(1:dof, 1:dof) )
       
       ! NEW
       !print*,'#DE$  AssociateElementReconstruction  ', elem%i,irhs, elem%xc(:)
       call SolveLocalMatrixProblem(dof, M(1:dof, 1:dof), 2*abs(irhs), M(1:dof, dof+1:dof+2*abs(irhs)) )

       ! setting of the basis coefficients of the flux reconstruction
       do k=1,abs(irhs)
          do j=1,2
             kk = (k-1)*2 + j
             kk1 = kk
             if(irhs == -1) then
                kk = 2 + j
                kk1 = j
             endif

             ! accumulation from all elements 
             ! OLD
             !elem%RTNphi(1:dof, kk) = elem%RTNphi(1:dof, kk) &
             !     + matmul(M(1:dof, 1:dof) , M(1:dof, dof+kk1) )

             ! NEW
             elem%RTNphi(1:dof, kk) = elem%RTNphi(1:dof, kk) + M(1:dof, dof+kk1) 

             !if(elem%i <=3) print*,'ATTENTION HERE !!!!!!'
          enddo
       enddo

       deallocate(M)
    endif
       
    !if(elem%i <= -27) then
    !   write(*,'(a6,3i5,30es12.4)') 'sig|_K',elem%i,i, ip, elem%RTNphi(1:10, 1)
    !   write(*,'(a6,3i5,30es12.4)') 'sig|_K',elem%i,i, ip, elem%RTNphi(1:10, 2)
    !   write(*,'(a6,3i5,30es12.4)') 'sig|_K',elem%i,i, ip, elem%RTNphi(1:10, 3)
    !   write(*,'(a6,3i5,30es12.4)') 'sig|_K',elem%i,i, ip, elem%RTNphi(1:10, 4)
    !   print*,'-------------------------'
    !endif
    

    ! RECOMPUTATION of the potential  (-R_{\pi/2} = \nabla s_h) 
    ! elem%RTNphi(1:dof, 3:4 ) contains \nabla s_h  in basis coeffs
    ! flux(2, 1:2, 1:Qdof) contains the same in integ nodes
    !!deg2 = deg + 1

    if(irhs == 2 .or. irhs == -1) then
       kk = 2
       if( irhs == -1 ) kk = 1

       dof2 =  DOFtriang( deg2 )  ! flux \in RTN_p \subset P_{p+1} thus s_h \in P_{p+2}
       allocate(pot_coeffs(1:dof2) )

       call SetPotential(elem, ie, V_rule, flux(kk, 1:2, 1:Qdof), deg2, dof2, &
            weights(1:Qdof), pot_coeffs(1:dof2) )
       
       ! accumulation of the potential from all elements  
       elem%RTNphi(1:dof2, 5) = elem%RTNphi(1:dof2, 5) + pot_coeffs(1:dof2)

       deallocate(pot_coeffs)
    endif

    

    ! the following is the graphical verification 
    if(igraf) then

       !write(*,*) '##############',ip,irhs

       ig_Fe = 80
       ig_Fv = 81

       ig_Pe = 90
       ig_Pv = 91

       !iFx = 3
       !iFy = 4

       dof2 = DOFtriang( deg1 + 1 ) 
       
       nface = 8  ! number of test points on the face
       allocate( Fx(1:Qdof, 1:nbDim), xi(3*nface, 1:3), Fxi(3*nface, 1:2)  )
       xi(:,:) = 0.
       it = 0
       do j=1,3
          j1 = mod(j,3) + 1
          j2 = mod(j1,3) + 1
          do k=0,nface-1
             it = it + 1
             xi(it, j2) = 1. - 1.*k /(nface-1)
             xi(it, j) = 1. - xi(it,j2)
          enddo
       enddo
       
       
       !integration nodes on K
       call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof,1:nbDim), &
            Fx(1:Qdof, 1:nbDim) )
       
       ! nodes on boundary of K
       call ComputeF(elem, 3*nface, xi(1:3*nface,1:nbDim), &
            Fxi(1:3*nface, 1:nbDim) )
       
       L_rule => state%space%L_rule(nface)
       
       allocate( Fxj(L_rule%Qdof, 1:2)  )
       
       ! nodes inside boundary of K
       call ComputeF(elem, L_rule%Qdof, L_rule%lambda(1:L_rule%Qdof, 1:2), &
            Fxj(1:L_rule%Qdof, 1:2) )
       
       ! DG test functions in xi
       allocate( phiT(1:dof2, 1:3*nface) ) 
       call PHI_orthonormal(3*nface, nbDim, xi(1:3*nface, 1:nbDim), 3, dof2, phiT(1:dof2, 1:3*nface) )
       
       ! DG test functions in Lagrange nodes
       allocate( phiL(1:dof2, 1:L_rule%Qdof), DphiLREF(1:dof2, 1:2, 1:L_rule%Qdof), &
            DphiL(1:dof2, 1:2, 1:L_rule%Qdof) ) 

       call PHI_orthonormal(L_rule%Qdof, nbDim, L_rule%lambda(1:L_rule%Qdof, 1:nbDim), 3, &
            dof2, phiL(1:dof2, 1:L_rule%Qdof), DphiLREF(1:dof2, 1:2, 1:L_rule%Qdof) )

       do j=1, dof2
          DphiL(j, 1, 1:L_rule%Qdof) = &
               elem%F%D1F0(1,1) * DphiLREF(j, 1, 1:L_rule%Qdof) &
               +elem%F%D1F0(1,2)* DphiLREF(j, 2, 1:L_rule%Qdof)

          DphiL(j, 2, 1:L_rule%Qdof) = &
               elem%F%D1F0(2,1) * DphiLREF(j, 1, 1:L_rule%Qdof) &
               +elem%F%D1F0(2,2)* DphiLREF(j, 2, 1:L_rule%Qdof)
       enddo
       
       ! fluxes on faces
       do j=1,Qdof
          write(50+irhs,*)  Fx(j, 1:2),  flux(1, 1:3, j)
          write(60+irhs,*)  Fx(j, 1:2),  flux(2, 1:3, j)
       enddo

       it = 0
       do j=1,3
          si = 1.
          if(elem%face(neigh, j) > elem%i) si = -1.
          do k=1,nface
             it = it + 1
             !write(83,*) xi(it, 1:3), Fxi(it, 1:2)
             if(irhs > 0) &
                  write(ig_Fe,*) Fxi(it, 1:2), & 
                  dot_product(phiT(1:dof,it), elem%RTNphi(1:dof, 1)), &
                  dot_product(phiT(1:dof,it), elem%RTNphi(1:dof, 2)), &
                  si*(dot_product(phiT(1:dof,it), elem%RTNphi(1:dof, 1))*elem%n(j,1) &
                  +dot_product(phiT(1:dof,it), elem%RTNphi(1:dof, 2))*elem%n(j,2)), elem%i

             if(irhs == 2 .or. irhs == -1) &
                  write(ig_Pe,*) Fxi(it, 1:2), & 
                  dot_product(phiT(1:dof,it), elem%RTNphi(1:dof, 3)), &
                  dot_product(phiT(1:dof,it), elem%RTNphi(1:dof, 4)), &
                  si*(-dot_product(phiT(1:dof,it), elem%RTNphi(1:dof, 3))*elem%n(j,2) &
                  +dot_product(phiT(1:dof,it), elem%RTNphi(1:dof, 4))*elem%n(j,1)), &
                  dot_product(phiT(1:dof2,it), elem%RTNphi(1:dof2, 5)), elem%i ! reconstruction s_h

          enddo
          !write(ig_Fe,'(x)')
          !write(ig_Fe,'(x)')
       enddo
       write(ig_Fe,'(x)')
       write(ig_Fe,'(x)')
       
       write(ig_Pe,'(x)')
       write(ig_Pe,'(x)')
       
       
       ! fluxes in Lagrang nodes
       ! first orientation
       it= 0
       do j=0, L_rule%Qdeg
          do k=0, L_rule%Qdeg - j
             it = it + 1
             if(irhs > 0) &
                  write(ig_Fv,*) Fxj(it, 1:2), & 
                  dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 1)), &
                  dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 2))

             if(irhs == 2 .or. irhs == -1) &
                  write(ig_Pv,*) Fxj(it, 1:2), & 
                  dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 3)), &
                  dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 4)), &
                  dot_product( DphiL(1:dof,2,it), elem%RTNphi(1:dof, 3)) &  ! divergence of flux
                  -dot_product(DphiL(1:dof,1,it), elem%RTNphi(1:dof, 4)), &
                  dot_product(phiL(1:dof2,it), elem%RTNphi(1:dof2, 5)), & ! reconstruction s_h
                  dot_product(DphiL(1:dof2,1, it), elem%RTNphi(1:dof2, 5)), & ! reconstruction s_h_x
                  dot_product(DphiL(1:dof2,2, it), elem%RTNphi(1:dof2, 5))   ! reconstruction s_h_y

             if(irhs == -1 .and. (it==1 .or. it==9 .or. it == 45) ) then
                if( abs(Fxj(it,1)*Fxj(it, 2)) < 1E-15) then ! Bondary node
                   allocate(val(1:ndim) )
                   call Exact_Scalar(Fxj(it, 1:2), val(1:ndim), state%time%ttime)
                   write(*,'(a6,3i5,4es12.4,a1,8es14.6)') &
                        'WEDx',ip,elem%i,it, L_rule%lambda(it, 1:2), Fxj(it,1:2), '|', &
                        dot_product(phiL(1:dof2,it), elem%RTNphi(1:dof2, 5)), val(1), &
                        dot_product(phiL(1:dof2,it), elem%RTNphi(1:dof2, 5)) - val(1)
                   write(138,*) &
                        Fxj(it,1:2), &
                        dot_product(phiL(1:dof2,it), elem%RTNphi(1:dof2, 5)), val(1), &
                        dot_product(phiL(1:dof2,it), elem%RTNphi(1:dof2, 5)) - val(1)
                   deallocate(val)
                endif

             endif

          enddo
          write(ig_Fv,'(x)')
          write(ig_Pv,'(x)')
       enddo
       write(ig_Fv,'(x)')
       write(ig_Pv,'(x)')
       print*,'*************       EGFDS'

       ! second orientation
       do j=0, L_rule%Qdeg
          it = nface + 1 - j 
          do k=0, L_rule%Qdeg - j
             !print*,'####',j,k,it
             if(irhs > 0) &
                  write(ig_Fv,*) Fxj(it, 1:2), & 
                  dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 1)), &
                  dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 2))

             if(irhs == 2 .or. irhs == -1) &
                  write(ig_Pv,*) Fxj(it, 1:2), & 
                  dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 3)), &
                  dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 4)), &
                  dot_product( DphiL(1:dof,2,it), elem%RTNphi(1:dof, 3)) &  ! divergence of flux
                  -dot_product(DphiL(1:dof,1,it), elem%RTNphi(1:dof, 4)), &
                  dot_product(phiL(1:dof2,it), elem%RTNphi(1:dof2, 5)), & ! reconstruction s_h
                  dot_product(DphiL(1:dof2,1, it), elem%RTNphi(1:dof2, 5)), & ! reconstruction s_h_x
                  dot_product(DphiL(1:dof2,2, it), elem%RTNphi(1:dof2, 5))   ! reconstruction s_h_y
             it = it + nface - k 
          enddo
          write(ig_Fv,'(x)')
          write(ig_Pv,'(x)')
       enddo
       
       write(ig_Fv,'(x)')
       write(ig_Pv,'(x)')
       
       ! third orientation
       do j=0, L_rule%Qdeg
          it = j + 1
          do k=0, L_rule%Qdeg - j
             !print*,'!!!!',j,k,it
             if(irhs > 0) &
                  write(ig_Fv,*) Fxj(it, 1:2), & 
                  dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 1)), &
                  dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 2))

             if(irhs == 2 .or. irhs == -1) &
                  write(ig_Pv,*) Fxj(it, 1:2), & 
                  dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 3)), &
                  dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 4)), &
                  dot_product( DphiL(1:dof,2,it), elem%RTNphi(1:dof, 3)) &  ! divergence of flux
                  -dot_product(DphiL(1:dof,1,it), elem%RTNphi(1:dof, 4)), &
                  dot_product(phiL(1:dof2,it), elem%RTNphi(1:dof2, 5)), & ! reconstruction s_h
                  dot_product(DphiL(1:dof2,1, it), elem%RTNphi(1:dof2, 5)), & ! reconstruction s_h_x
                  dot_product(DphiL(1:dof2,2, it), elem%RTNphi(1:dof2, 5))   ! reconstruction s_h_y
             it = it + nface +1 - k 
          enddo
          write(ig_Fv,'(x)')
          write(ig_Pv,'(x)')
       enddo
       
       write(ig_Pv,'(x)')
       write(ig_Fv,'(x)')
       
       
       
       deallocate(Fx, xi, Fxi, Fxj)

       deallocate(phiL, phiT, DphiL, DphiLREF)
  
    endif   ! if(igraf)   graphical verification


    deallocate(flux, weights)
    ! end of reconstructures 


    !if(ip == 15) then
    !   print*,'Stopped in  AssociateElementReconstruction'
    !   stop
    !endif

  end subroutine AssociateElementReconstruction

  !> create the pair of RTN DOF on elem (local index) and the patch \f$\omega_a \f$ (global index)
  !> RTN DOF are face momentum on two faces sharing a vertex and volume momentums
  !> itrans(:,:, 1:2)  inner OR flux
  !> itrans(:,:, 3:4)  bound AND potential
  subroutine CreateTransfPairs(deg, F_size, N, i, ie, itrans, inner)
    integer, intent(in) :: deg   ! deg of recontruction
    integer, intent(in) :: F_size   ! number of RTN DOF on one element
    integer, intent(in) :: N  ! number of elements in supp around the vertex 
    integer, intent(in) :: i  ! index  of the given elem in supp around the vertex 
    integer, intent(in) :: ie ! local index of this vertex
    integer, dimension(1:F_size, 1:4), intent(inout) :: itrans
    logical, intent(in) :: inner !  inner node
    integer :: i_n, j, ie1, ie2, deg1, F_face, F_vol, F_tot

    if(inner) then
       i_n = mod(i, N) + 1   ! closed list of elements
    else
       i_n = i + 1           ! open list of elements
    endif

    ie1 = mod( ie, 3) + 1
    ie2 = mod(ie1, 3) + 1
    
    deg1 = deg + 1
    F_face = deg1
    F_vol  = deg * deg1
    F_tot = F_face + F_vol

    if(F_size /=  2*deg1+F_vol) then  ! size of assembled submatrix
       print*,' Inconsistency in the number of DOF in neumann.f90'
       stop
    endif

    ! itrans indexes of transformation
    !> itrans( : , 1 ) local  vector / matrix 
    !> itrans( : , 2 ) global vector / matrix
    !allocate(itrans(1: F_size, 1:2) ) 

    ! local edges
    do j=1,deg1
       itrans(j       , 1) = (ie  - 1)*deg1 + j
       itrans(j + deg1, 1) = (ie2 - 1)*deg1 + j
    enddo

    ! local volumes
    do j=1,deg1*deg
       itrans(j +2*deg1, 1) = 3*deg1 + j
    enddo

    ! global edges
    do j=1,deg1
       itrans(j       , 2) = (i   - 1)*F_tot + j
       itrans(j + deg1, 2) = (i_n - 1)*F_tot + deg1 - j + 1
    enddo
    ! global volumes
    do j=1,deg1*deg
       itrans(j +2*deg1, 2) = (i   - 1)*F_tot +  F_face + j
    enddo

    !write(*,'(a4,l2,i5,a2,30i5)') 'CTF',inner, i,'|', itrans(1:F_size, 1)
    !write(*,'(a6,i5,a2,30i5)') '    ', i,'|', itrans(1:F_size, 2)
    !print*

    ! potential reconstruction for boundary nodes
    ! the first and the last edges have to be removed
    if(.not. inner) then

       itrans(1:F_size, 3:4) = itrans(1:F_size, 1:2)
       ! we remove the first edge, hence we shift the global array
       itrans(1:F_size, 4) = itrans(1:F_size, 4) - deg1
       if(i == 1) itrans(1:deg1, 3) = 0

       ! we remove the last edge
       if(i == N) itrans(deg1 + 1: 2* deg1, 3) = 0

       do j=1,F_size
          if(itrans(j, 3) == 0) itrans(j, 4) = -1
       enddo

       !write(*,'(a4,l2,i5,a2,30i5)') 'CTF2',inner, i,'|', itrans(1:F_size, 3)
       !write(*,'(a6,i5,a2,30i5)') '    ', i,'|', itrans(1:F_size, 4)
       !print*

    endif

  end subroutine CreateTransfPairs


  !> setting of potential (set in elem%RTNphi(2:dof2, 5)) from its gradient
  subroutine SetPotential(elem, ie, V_rule, flux, deg2, dof2, weights, pot_coeffs)
    class(element), intent(inout) :: elem
    integer :: ie    ! inner index of the leading vertex
    type(volume_rule), target, intent(in) :: V_rule
    real, dimension(1:2, 1:V_rule%Qdof), intent(inout) :: flux
    integer, intent(in) :: deg2  ! polynomial degree of flux, potential has degree deg + 2
    integer, intent(in) :: dof2  ! corresponding dof
    real, dimension(1:V_rule%Qdof), intent(in) :: weights
    real, dimension(1:dof2), intent(inout) :: pot_coeffs 
    real, dimension(:), allocatable :: potential
    real, dimension(:, :), allocatable :: Fx
    real, dimension(:, :), allocatable :: psi
    real, dimension(:, :), allocatable :: M
    real, dimension(:, :), allocatable :: aij, sij, zij
    real, dimension(:,:), allocatable :: xi
    real, dimension(:,:), allocatable :: phiE
    real, dimension(:,:), pointer :: phi
    real :: hs, val
    integer:: Qdof
    integer :: i,j,k, l, QEdof, j1, j2, it


    hs =  elem%diam**2  ! scaling factor
    !hs = 1.
    !print*,'ATEnTiOn !@!!!!^%%!!!'

    Qdof = V_rule%Qdof
    allocate( potential(1:Qdof)) 

    !physical coordinates of integ node
    allocate( Fx(1:Qdof, 1:nbDim)) 
    call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof,1:nbDim),  Fx(1:Qdof, 1:nbDim) )

    ! test
    !flux(:,:) = 0.
    !do l=1,Qdof
    !   flux(1, l) = 2*Fx(l,2) * Fx(l,1)
    !   flux(2, l) = Fx(l,1)**2 
    !enddo
    ! test end


    !if(elem%i <= 15) then
    !   write(*,'(a6,i5,50es12.4)') 'D_x s',elem%i, flux( 1, 1:Qdof)
    !   write(*,'(a6,i5,50es12.4)') 'D_y s',elem%i, flux( 2, 1:Qdof)
    !   !print*
    !   !write(*,'(a6,i5,30es12.4)') 'D_x s',elem%i, elem%RTNphi(1:10, 3)
    !   !write(*,'(a6,i5,30es12.4)') 'D_y s',elem%i, elem%RTNphi(1:10, 4)
    !   print*,'-----------------', deg2, dof2
    !endif



    allocate(psi(1:dof2, 1:Qdof))  ! Taylor basis functions on elements
    allocate( xi( 1:3, 1:2) )

    do l=1,Qdof
       xi(1, 1:2) = (Fx(l, 1:2) - elem%xc(1:2) ) / hs
       call phi_values(dof2, xi(1, 1:2), psi(1:dof2, l) )
    enddo


    allocate(M(1:dof2, 1:dof2 + 2) )  ! mass matrix and RHS of reconstructions
    M(:,:) = 0.
    do j =1, dof2
       ! mass matrix = (\phi_i, \phi_j)_K
       do k=j, dof2
          M(j,k) = dot_product(weights(1:Qdof), psi(j, 1:Qdof) * psi(k, 1:Qdof))
          M(k,j) = M(j,k)
       enddo
       ! RHS = (flux, \phi_i)_K
       do k=1,2
          M(j,dof2 + k ) = dot_product(weights(1:Qdof), psi(j, 1:Qdof) * flux(k, 1:Qdof) )
       enddo
    enddo
         
    !!call WriteArray(M, dof2, dof2+2, 1, dof2, 1, dof2+2)

    ! RAD
    !call MblockInverse(dof2, M(1:dof2, 1:dof2) )

    !!call WriteArray(M, dof2, dof2, 1, dof2, 1, dof2)

    call SolveLocalMatrixProblem(dof2, M(1:dof2,1:dof2), 2, M( 1:dof2, dof2 + 1: dof2 + 2) )

    allocate(aij(1:2, 1:dof2) )  ! coefficients of flux in Taylor basis

    aij(1, 1:dof2) = M( 1:dof2, dof2 + 1)
    aij(2, 1:dof2) = M( 1:dof2, dof2 + 2)

    !RAD
    !aij(1, 1:dof2) = matmul(M(1:dof2, 1:dof2), M( 1:dof2, dof2 + 1) )
    !aij(2, 1:dof2) = matmul(M(1:dof2, 1:dof2), M( 1:dof2, dof2 + 2) )

    !if(elem%i <= 15) then
    !   print*
    !   write(*,'(a6,i5,30es12.4)') 'D_x s',elem%i, aij(1, 1:dof2)
    !   write(*,'(a6,i5,30es12.4)') 'D_y s',elem%i, aij(2, 1:dof2)
    !endif

   allocate(sij(0:deg2, 0:deg2), zij(0:deg2, 0:deg2) )  ! coefficients of potential in Taylor basis

   sij(:,:) = 0.
   zij(:,:) = 0.

   l = 0
   do k=1, deg2
      do j = 0, k - 1
         i = k - j
         l = l + 1
         sij(i,j) = aij(1, l) * hs / i 
         !!print*,'#####',l,i, j, k, '/' , i
      enddo
   enddo
   !print*,'-----------------------------------'

   !!call WriteArray(zij(0:deg2, 0:deg2), deg2+1, deg2+1, 1, deg2+1, 1, deg2+1)


   l = 0
   do k=1, deg2
      do i = k - 1, 0, -1
         j = k - i
         l = l + 1
         zij(i,j) = aij(2, l) * hs / j
         !!print*,'#####',l,i, j, k, '/' , j
      enddo
   enddo

   !call WriteArray(sij(0:deg2, 0:deg2), deg2+1, deg2+1, 1, deg2+1, 1, deg2+1)

   !call WriteArray(zij(0:deg2, 0:deg2), deg2+1, deg2+1, 1, deg2+1, 1, deg2+1)


   ! averaging
   do i = 0, deg2
      do j = 1, deg2 
         if(i == 0) then
            sij(i,j) = zij(i,j)
         else
            sij(i,j) = ( sij(i,j) +  zij(i,j) ) / 2
         endif
      enddo
   enddo
   
   !!call WriteArray(sij(0:deg2, 0:deg2), deg2+1, deg2+1, 1, deg2+1, 1, deg2+1)

   ! setting of potential in integ nodes
   potential(:) = 0
   l = 0
   do k = 0, deg2
      do j= 0, k
         i = k - j
         l = l + 1
         !write(*,'(a10,6i5)') 'sum_pot:',l, i, j
         potential(1:Qdof) = potential(1:Qdof) + psi(l, 1:Qdof) * sij(i, j)
      enddo
   enddo


   ! projection of the potential into DGM basis functions

    phi => V_rule%phi(1:dof2, 1:Qdof)

    ! basis coeffs of the potential
    M(:,:) = 0.
    do j =1, dof2
       ! mass matrix = (\phi_i, \phi_j)_K
       do k=j, dof2
          M(j,k) = dot_product(weights(1:Qdof), phi(j, 1:Qdof) * phi(k, 1:Qdof))
          M(k,j) = M(j,k)
       enddo
       ! RHS = (flux, \phi_i)_K
       M(j,dof2 + 1) = dot_product(weights(1:Qdof), phi(j, 1:Qdof) * potential(1:Qdof) )
    enddo

    !!call WriteArray(M, dof2, dof2+2, 1, dof2, 1, dof2+2)

    call MblockInverse(dof2, M(1:dof2, 1:dof2) )
    M(1:dof2, dof2 + 2) =  matmul(M(1:dof2, 1:dof2) , M(1:dof2, dof2+1) )

    ! NEW  slowly
    !!M(1:dof2, dof2 + 2) = M(1:dof2, dof2 + 1)
    !!call SolveLocalMatrixProblem(dof2, M(1:dof2, 1:dof2), 1, M(1:dof2, dof2+2 : dof2+2) )
    

    !!call WriteArray(M, dof2, dof2+2, 1, dof2, 1, dof2+2)

    ! value of the potential on the boundary edge
    deallocate (xi)
    QEdof = 3

    allocate(xi(1:QEdof, 1:3), phiE(1:dof2, 1:QEdof) )
    xi(:,:)= 0.
    j1 = mod(ie,3) + 1
    j2 = mod(j1,3) + 1

    it = 0
    do k=0,QEdof-1
       it = it + 1
       xi(it, ie) = 1. - 1.*k /(QEdof-1)
       xi(it, j1) = 1. - xi(it,ie)
    enddo

    !call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof,1:nbDim),  Fx(1:Qdof, 1:nbDim) )
    !do k=1,Qdof
    !   write(62,*) Fx(k, 1:2), phi(1:dof2, k)
    !enddo

    ! value of DG basis functions on the boundary of K \subset \partial \omega_a
    call PHI_orthonormal(QEdof, nbDim, xi(1:QEdof, 1:nbDim), 3, dof2, phiE(1:dof2, 1:QEdof) )

    call ComputeF(elem, QEdof, xi(1:QEdof,1:nbDim),  Fx(1:QEdof, 1:nbDim) )

    ! setting of the trace of potential on the edge
    val = 0.
    do k=1,QEdof
       val = val +  dot_product( phiE(1:dof2, k),  M(1:dof2, dof2 + 2) )
    enddo
    val = val /QEdof

    ! shifting of the potential
    potential(1:Qdof) = potential(1:Qdof) - val

    ! reconstruction of the potential in DG basis functions
    do j = 1, dof2
       ! RHS = (flux, \phi_i)_K
       M(j,dof2 + 1) = dot_product(weights(1:Qdof), phi(j, 1:Qdof) * potential(1:Qdof) )
    enddo

    pot_coeffs(1:dof2) =  matmul(M(1:dof2, 1:dof2) , M(1:dof2, dof2+1) )

    ! NEW
    !!call SolveLocalMatrixProblem(dof2, M(1:dof2, 1:dof2), 1, M(1:dof2, dof2+1 : dof2+1) )
    !!pot_coeffs(1:dof2) = M(1:dof2, dof2 + 1)


   !  if(elem%i == 1 .and. ie == 3) then
   !     print*,'###', elem%i, ie
       
   !     do l=1,Qdof
   !        write(65, *) Fx(l, 1:2), flux(1:2, l)
   !        write(66, *) Fx(l, 1:2), potential(l),  psi(1:dof2, l)
   !     enddo
       
   !     do k=1,QEdof
   !        print*,'@@@',k, dot_product( phiE(1:dof2, k),  M(1:dof2, dof2 + 2) )
   !        write(61,*) Fx(k, 1:2),  dot_product( phiE(1:dof2, k),  M(1:dof2, dof2 + 2) ),xi(k, 1:3)
          
   !        !write(63,*) Fx(k, 1:2), phiE(1:dof2, k)
          
   !     enddo
   !  endif
   !     !stop

   ! !  elem%RTNphi(1, 5) = 0.  ! to be set

     
     deallocate(M, Fx, xi, psi)

  end subroutine SetPotential


  !> evaluate the nonhomogeneous Dirichlet BC in momentum boundary nodes
  subroutine Eval_Dir_BC(elem, deg, F_face, ie, first, uD)
    class(element), intent(in) :: elem      
    integer, intent(in) :: deg, F_face   ! F_case should be deg + 1
    integer, intent(in) ::  ie  ! inner index of boundary edge
    logical, intent(in) :: first  ! first edge, vertex is the first point
    real, dimension(1:F_face), intent(inout) ::uD
    real, dimension(:,:), allocatable :: xi, Fxi
    real, dimension(:), allocatable :: psi
    real, dimension(:,:, :), allocatable :: Dwi
    integer :: Qdof, dof
    integer :: ie0,  ie1, ie2, j, ib, it, i
    integer :: indx, iphi
    real ::  rlen, Der_u

    allocate( xi(1: F_face, 1:3)) ! barycentric coordinates of nodes on faces
    allocate(Fxi(1: F_face, 1:2)) ! physical  coordinates of nodes on faces
    allocate( Dwi(1: F_face, 1:ndim, 0:2)) ! solution and its derivatives
    allocate( psi(1: F_face ) )
    ! barycentres of 
    xi(:,:) = 0.
    it = 0   


    if( first) then  ! setting of the corresponding boundary edge
       ie0 = ie
    else
       ie1 = mod(ie , 3) + 1
       ie0 = mod(ie1, 3) + 1
    endif

    !write(*,'(a9,4i5)' ) 'evalBC',elem%i, ie, ie0
    
    ! face momentums, setting of the nodes on faces, where $\psi\cdot\nn$ is given
    !do ie=1,3 ! loop over triagle edges (including HG nodes )
    ie1 = mod(ie0 , 3) + 1
    ie2 = mod(ie1, 3) + 1

       
    if(deg == 0) then  ! only one node in face centers
       it = it + 1
       xi(it, ie2) = 0.5
       xi(it, ie0) = 1. - xi(it,ie2)
       !write(*,'(a6,4i5, 4es12.4)') '!@#',elem%i,ie, ie1, ie2,  xi(it,:)
       psi(it) = 0.5
    else
       rlen = 1./deg
       do j=0, deg
          it = it + 1
          xi(it, ie2) = 1. - j * rlen
          xi(it, ie0) = 1. - xi(it,ie2)
          !write(*,'(a6,4i5, 4es12.4)') '!@#',elem%i,ie, ie1, ie2,  xi(it,:)

          ! the hat function
          if(first) then
             psi(it) =  1. - j * rlen
          else
             psi(F_face-j) =  1. - j * rlen
          end if

       enddo
    endif

    if(it /= F_face) stop 'Troubles in QWSTYBVC'

    call ComputeF(elem, F_face, xi(1:F_face, 1:nbDim), Fxi(1:F_face, 1:nbDim) )

    do j=1, F_face
       ! ONLY FOR SCALAR !!!!
       ! exact solution (i.e., boundary condition) and its derivative
       call Exact_Scalar(Fxi(j, 1:2), Dwi(j,1:1, 0), state%time%ttime)

       call Der_Exact_Scalar(Fxi(j, 1:2), Dwi(j,1:1, 1:2), state%time%ttime)

       ! the tangential derivative
       Der_u =  (Dwi(j,1, 1) * elem%n(ie0,2) -  Dwi(j,1, 2) * elem%n(ie0,1)) / elem%dn(ie0)

       if( first) then
          uD(j) =   Dwi(j,1, 0) / elem%dn(ie0) + Der_u * psi(j) 
          !uD(F_face-j+1) =   Dwi(j,1, 0) / elem%dn(ie0) + Der_u * psi(j) 
          !uD(j) = - uD(j)
       else
          uD(j) =   Dwi(j,1, 0) / elem%dn(ie0) - Der_u * psi(j) 
          !uD(F_face-j+1) =   Dwi(j,1, 0) / elem%dn(ie0) - Der_u * psi(j) 
          !uD(j) = - uD(j)
       endif
       !print*,'##', Der_u * psi(j),  Dwi(j,1, 0) / elem%dn(ie0), uD(j)

       !write(*,'(2es12.4,a2,6es12.4,a3,6es14.6)')  Fxi(j, :), '|', Dwi(j,1, 0:2),Der_u, &
       !    psi(j), elem%dn(ie0),'!',  & ! Dwi(j,1, 0) / elem%dn(ie0) , & ! Der_u * psi(j), 
       !    uD(j) !, 4 + 2*Fxi(j,1) - 12*Fxi(j,1)**2, 4 + 2*Fxi(j,2) - 12*Fxi(j,2)**2
       !write(15,*) Fxi(j, :), uD(j)
    enddo
    !write(15,'(x)')

    deallocate(xi, Fxi, Dwi, psi)
  end subroutine Eval_Dir_BC


  !> evaluate the nonhomogeneous Dirichlet BC in momentum boundary nodes
  subroutine Eval_Dir_BC_L2proj(elem, deg, F_face, ie, first, uD)
    class(element), intent(in) :: elem      
    integer, intent(in) :: deg, F_face   ! F_case should be deg + 1
    integer, intent(in) ::  ie  ! inner index of boundary edge
    logical, intent(in) :: first  ! first edge, vertex is the first point
    type(Gauss_rule), pointer :: G_rule
    type(Time_rule), pointer :: T_rule
    real, dimension(1:F_face), intent(inout) ::uD
    real, dimension(:,:), allocatable :: xi, Fxi, Lphi
    real, dimension(:), allocatable :: psi, uDi, Li, yi
    real, dimension(:,:, :), allocatable :: Dwi
    integer :: Gdof, dof
    integer :: ie0,  ie1, ie2, j, ib, it, i, l
    integer :: indx, iphi
    real ::  rlen, Der_u, ti

    !G_rule => state%space%G_rule(max( maxGrule , deg+2) )
    G_rule => state%space%G_rule(max( maxGrule , deg+3) )  ! (deg+2) is no
    Gdof  = G_rule%Qdof

    !print*,'##EDE',deg+2, Gdof

    allocate( xi(1: Gdof, 1:3)) ! barycentric coordinates of nodes on faces
    allocate(Fxi(1: Gdof, 1:2)) ! physical  coordinates of nodes on faces
    allocate( Dwi(1: Gdof, 1:ndim, 0:2)) ! solution and its derivatives
    allocate( psi(1: Gdof ), uDi(1: Gdof ) )
    ! barycentres of 
    xi(:,:) = 0.
    it = 0   


    if( first) then  ! setting of the corresponding boundary edge
       ie0 = ie
    else
       ie1 = mod(ie , 3) + 1
       ie0 = mod(ie1, 3) + 1
    endif

    !write(*,'(a9,4i5)' ) 'G_evalBC',elem%i, ie, ie0
    
    ! face momentums, setting of the nodes on faces, where $\psi\cdot\nn$ is given
    !do ie=1,3 ! loop over triagle edges (including HG nodes )
    ie1 = mod(ie0 , 3) + 1
    ie2 = mod(ie1, 3) + 1

    do l= 1,Gdof
       ti = G_rule%lambda(l)
       xi(l, ie2) = 1. - ti
       xi(l, ie0) = 1. - xi(l,ie2)
       !!write(*,'(a6,4i5, 4es12.4)') '**@#',elem%i,ie, ie1, ie2,  xi(l,:)

       ! the hat function
       if(first) then
          psi(l) =  1. - ti
       else
          psi(l) =  ti
       end if

    enddo

    call ComputeF(elem, Gdof, xi(1:Gdof, 1:nbDim), Fxi(1:Gdof, 1:nbDim) )

    ! evaluation of the projection in the Gauss integ nodes
    do j=1, Gdof
       ! ONLY FOR SCALAR !!!!
       ! exact solution (i.e., boundary condition) and its derivative
       call Exact_Scalar(Fxi(j, 1:2), Dwi(j,1:1, 0), state%time%ttime)

       call Der_Exact_Scalar(Fxi(j, 1:2), Dwi(j,1:1, 1:2), state%time%ttime)

       ! the tangential derivative
       Der_u =  (Dwi(j,1, 1) * elem%n(ie0,2) -  Dwi(j,1, 2) * elem%n(ie0,1)) / elem%dn(ie0)

       if( first) then
          uDi(j) =   Dwi(j,1, 0) / elem%dn(ie0) + Der_u * psi(j) 
          !uDi(j) = -uDi(j) 
       else
          uDi(j) =   Dwi(j,1, 0) / elem%dn(ie0) - Der_u * psi(j) 
          !uDi(j) = -uDi(j) 
       endif
       !print*,'##', Der_u * psi(j),  Dwi(j,1, 0) / elem%dn(ie0), uDi(j)

       !write(*,'(2es12.4,a2,6es12.4,a3,6es14.6)')  Fxi(j, :), '|', Dwi(j,1, 0:2),elem%n(ie0,1:2), Der_u !, &
       !    psi(j), elem%dn(ie0),'!',  & ! Dwi(j,1, 0) / elem%dn(ie0) , & ! Der_u * psi(j), 
       !    uDi(j) !, 4 + 2*Fxi(j,1) - 12*Fxi(j,1)**2, 4 + 2*Fxi(j,2) - 12*Fxi(j,2)**2
       !write(17,*) Fxi(j, :), uDi(j)
    enddo

    allocate( Li(0: deg ) )
    allocate( yi(1: F_face ) )

    ! evaluation of the projection in the basis of Legendre function
    do i=0,deg
       Li(i) = dot_product(G_rule%weights(1:Gdof),  G_rule%Leg_phi(i, 1:Gdof)*uDi(1:Gdof) )
       Li(i) = Li(i) / G_rule%Leg_phi(i,0)
    enddo

    
    ! evaluation of the projection in the 1D Lagrange nodes (=face  momentums)
    if(deg+1 /= F_face) then
       print*,'TROUBLE in Eval_Dir_BC_L2proj'
       stop
    endif

    deallocate(xi, Fxi)
    allocate( xi(1: deg+1, 1:3)) ! barycentric coordinates of nodes on faces
    allocate(Fxi(1: deg+1, 1:2)) ! physical  coordinates of nodes on faces

    xi(:,:) = 0.

    allocate( Lphi(0:deg, 1:F_face) )
    it = 0
    if(deg == 0) then  ! only one node in face centers
       it = it + 1
       yi(it) = 0.5
    else
       rlen = 1./deg
       do j=0, deg
          it = it + 1
          if(first) then
             !yi(it) = 1. - j * rlen
             yi(it) =  j * rlen
          else
             !yi(it) = 1. - j * rlen
             yi(it) =  j * rlen
          endif
          
          xi(it, ie2) = 1. - yi(it)
          xi(it, ie0) = 1. - xi(it,ie2)

          !write(*,'(a6,4i5, 4es12.4)') '!@#',elem%i,ie, ie1, ie2,  xi(it,:)
       enddo
    endif
    
    call ComputeF(elem, deg+1, xi(1:deg+1, 1:nbDim), Fxi(1:deg+1, 1:nbDim) )
    !do j=1, deg+1
    !   write(*,'(a6,4i5, 4es12.4)') 'WE@@#',elem%i,ie, ie1, ie2,  xi(j,1:2), Fxi(j,:)
    !   !write(*,'(95, 4es12.4)') 'WE@@#',elem%i,ie, ie1, ie2,  xi(j,1:2), Fxi(j,:)
    !enddo

    !print*,'------------------------------', F_face, deg
    
    call Eval_Legenre_pols(G_rule, deg, F_face, yi(1:F_face), Lphi(0:deg, 1:F_face) )

    uD(1:F_face) = matmul(Li(0: deg ),  Lphi(0:deg, 1:F_face) ) 


    !do j=1, deg+1
    !   write(*,'(a10,4es12.4)') 'wD at xi:',Fxi(j,1:2), uD(j)
    !   write(94,'(4es12.4)') Fxi(j,1:2), uD(j)
    !enddo

    !do j=1,F_face
       !write(18,*) yi(j), uD(j)
    !   write(18,*) yi(j)*0.25, uD(j)
       !write(*,'(a6,4i5, 40es12.4)') '!@#',elem%i,ie, ie1, ie2, yi(j),Lphi(0:deg,j)
    !enddo

    !stop

    deallocate(xi, Fxi, Dwi, psi, uDi, Li, yi, Lphi)

  end subroutine Eval_Dir_BC_L2proj

  !> evaluation of the estimator measuring the violation of BC
  !> see Vohralik (aposteriori), Theorem 7.1.5, page 71
  !>  direct integration over sub-element 
  subroutine Eval_BC_estim(elem, V_rule)
    class(element), intent(inout) :: elem
    type(volume_rule), intent(in) :: V_rule
    real, dimension(:,:), allocatable :: xi, xt, phiQ, Mb, phiBC, ff
    real, dimension(:,:,:), allocatable :: DphiQ
    real, dimension(:), allocatable :: lam, bar, barQ
    integer :: i, ib, l, j, k1, k2, Qdof, Fdeg, Fdof
    real :: val, weigh

    Fdeg = elem%RTNflux_deg
    Fdof = DOFtriang( Fdeg )  

    Qdof = V_rule%Qdof
 
    allocate( ff(1:Qdof, 1:3*ndim) )

    allocate( xi(1:5, 1:2), xt(1:3, 1:2), lam(1:2), bar(1:3), barQ(1:3) )
    allocate( phiQ(1:Fdof, 1:Qdof),  DphiQ(0:Fdof, 1:2, 1:Qdof), Mb(1:Fdof, 1:Fdof+2) )
    allocate( phiBC(1:Fdof, 1:Qdof) )

    ! physical cordinates of elem
    xt(1, 1:2) = grid%x(elem%face(idx, 1), 1:2)
    xt(2, 1:2) = grid%x(elem%face(idx, 2), 1:2)
    xt(3, 1:2) = grid%x(elem%face(idx, 3), 1:2)
  

    do i=1, elem%flen
       if(elem%face(neigh, i) < 0 ) then !boundary edge
          ib =  abs(elem%face(neigh,i))
          k1 = grid%b_edge(ib)%lbn(1)
          k2 = grid%b_edge(ib)%lbn(2)

          ! coordinates of the subtriangle
          xi(1,1:2) = grid%x(k1,1:2)
          xi(2,1:2) = grid%x(k2,1:2)

          ! third point of the triangle 
          ! barycentre
          xi(3,1:2) = elem%xc(1:2)    
          weigh = 3.

          ! opposite vertex
          !i1 = mod(i, 3) + 1
          !i2 = mod(i1, 3) + 1
          !xi(3,1:2) = grid%x(elem%face(idx,i2),1:2)
          !weigh = 1.

          do l =1, Qdof
             ! coordinates of integ node
             xi(5, 1:2) = V_rule%lambda(l,1) * xi(1, 1:2) &
                  + V_rule%lambda(l,2) * xi(2, 1:2) &
                  + V_rule%lambda(l,3) * xi(3, 1:2)

             ! coordinates of the projection 
             lam(1:2) = V_rule%lambda(l,1:2) / ( V_rule%lambda(l,1) +  V_rule%lambda(l,2) )
             xi(4, 1:2) = lam(1) *  xi(1,1:2)  + lam(2) *  xi(2,1:2) 


             ! barycentric coordinates of xi(4, :) within elem
             call BarycCoordOne(xt(1:3, 1:2), xi(4, 1:2), bar(1:2) )

             ! nodes on the boundary after the projection
             call Eval_phi_Qnode(elem, Fdof, 1, bar(1:nbDim), phiBC(1:Fdof,l)  )

             ! potential at the boundary
             ff(l, 2) = dot_product(phiBC(1:Fdof, l) , elem%RTNphi(1:Fdof,5)) 
             !write(*,'(a6,i5,10es12.4)') 'weLAM',l, xi(5,1:2), bar(1:2), 1-bar(1) - bar(2)

             ! exact solution at the boundary
             call Exact_Scalar(xi(4, 1:2), ff(l, 1), state%time%ctime)


             ! exact solution (= BC) - potential
             ! extrapolation at integ node - linear interpol between boundary (=1) and barycentre(=0)
             !print *,'ATTENTION EHDTAIDYETSJ'
             !ff(l,3) = 1. & !(ff(l,1) - ff(l,2) ) &
             !ff(l,3) = lam(1)**2 & !(ff(l,1) - ff(l,2) ) &
             ff(l,3) = (ff(l,1) - ff(l,2) ) &
                  * VectorNorm(xi(5,1:2) - xi(3,1:2)) / VectorNorm(xi(4,1:2) - xi(3,1:2) )


             !write(*,'(a6,i5,10es12.4)') 'weLEA',l,xi(4, 1:2), ff(l, 1:3)


             ! nodes in the element interior 
             call BarycCoordOne(xt(1:3, 1:2), xi(5, 1:2), barQ(1:2) )
             call Eval_phi_Qnode(elem, Fdof, 1, barQ(1:nbDim), phiQ(1:Fdof,l), DphiQ(1:Fdof,1:2,l) )

             
             !write(90,'(30es16.8)')  xi(5, 1:2), ff(l,3)
             !write(91,'(30es16.8)')  xi(4, 1:2), ff(l,1) - ff(l,2) ,lam(1:2), lam(1:2)**2
             !write(91,'(30es16.8)')  xi(3, 1:2), 0.
             !write(91,'(x)') 
             !write(91,'(x)') 
             !write(91,'(x)') 

             !write(100+elem%i,'(30es16.8)')  xi(5, 1:2), ff(l,3)

          enddo
          !close(41)

          ! array ff(l, 3) contains  the values of (u_D - s_h) at integ nodes on subtriangle
          ! we compute its gradient
          
          ! projection into polynomials
          do l=1, Fdof
             do j=1,Fdof
                Mb(l,j) = dot_product(V_rule%weights(1:Qdof), phiQ(l,1:Qdof) * phiQ(j,1:Qdof))
             enddo
             Mb(l,Fdof+1) = dot_product(V_rule%weights(1:Qdof), phiQ(l,1:Qdof) * ff(1:Qdof, 3))
             !!write(*,'(a6,2i5,60es12.4)') 'Mb:DS',Fdof, l,  Mb(l,1:Fdof+1)
          enddo

          !call MblockInverse(Fdof, Mb(1:Fdof, 1:Fdof))
          !Mb(1:Fdof, Fdof+2) = matmul(Mb(1:Fdof, 1:Fdof), Mb(1:Fdof, Fdof+1) )

          call SolveLocalMatrixProblem(Fdof, Mb(1:Fdof, 1:Fdof), 1, Mb(1:Fdof, Fdof+1:Fdof+1) )
          Mb(1:Fdof, Fdof+2) = Mb(1:Fdof, Fdof+1)

          ! graident u_D - s_h at integ nodes 
          DphiQ(0, 1, 1:Qdof) = matmul( Mb(1:Fdof, Fdof+2) ,  DphiQ(1:Fdof, 1, 1:Qdof))
          DphiQ(0, 2, 1:Qdof) = matmul( Mb(1:Fdof, Fdof+2) ,  DphiQ(1:Fdof, 2, 1:Qdof))


          val = elem%area/ weigh * dot_product(V_rule%weights(1:Qdof), &
               DphiQ(0,1,1:Qdof) * DphiQ(0,1,1:Qdof) + DphiQ(0,2,1:Qdof) * DphiQ(0,2,1:Qdof) )

          elem%eta(P_BC, 1) = elem%eta(P_BC, 1) + val**0.5

          ! print*,'Qdof=',Qdof, Fdeg, Qnum
          ! write(*,'(a6,90es12.4)' ) 'xi', ff(1:Qdof, 3)
          ! write(*,'(a6,90es12.4)' ) 'Dwi',  DphiQ(0, 1, 1:Qdof)
          ! write(*,'(a6,90es12.4)' ) 'Dwi',  DphiQ(0, 2, 1:Qdof)
          ! write(*,'(a6,40es12.4)' ) 'wi', Mb(1:Fdof, Fdof+2)
          ! write(*,'(a6,40es12.4)' ) 'val', val**0.5, val, &
          !      dot_product(V_rule%weights(1:Qdof), &
          !      matmul( Mb(1:Fdof, Fdof+2) ,  phiQ(1:Fdof, 1:Qdof))), &    ! 3 !
          !      dot_product(V_rule%weights(1:Qdof), &
          !      DphiQ(0,1,1:Qdof) * DphiQ(0,1,1:Qdof) + DphiQ(0,2,1:Qdof) * DphiQ(0,2,1:Qdof) ), & !4
          !      elem%area , 1./(1./8 * 2./3)  !5
          ! print*,elem%xc,  sum(V_rule%weights(1:Qdof) )
          ! stop
       endif

    enddo
    deallocate(xi, xt, lam, bar, barQ, phiQ, phiBC, DphiQ, Mb, ff)
    !end BC according Vohralik (aposteriori), Theorem 7.1.5, page 71



  end subroutine Eval_BC_estim

  ! drawing of the reconstruction
  subroutine DRAW_reconstruction(ifile, elem, nface, Fdof, idx)
    class(element), intent(in) :: elem
    integer, intent(in) :: ifile !  file index
    integer, intent(in) :: nface ! number of test points on the face
    integer, intent(in) :: Fdof
    integer, intent(in) :: idx ! index of component
    real, dimension(:, :), allocatable :: Fxj, xi, Fxi
    integer :: i, j, j1, j2, it, k, Qdof
    type(Lagrang_rule), pointer :: L_rule
    real, dimension(:,:,:), allocatable :: u_exact
    real, dimension(:,:), allocatable :: phiT,phiL
    real, dimension(:,:,:), allocatable :: DphiL, DphiLREF


    L_rule => state%space%L_rule(nface)
    Qdof = L_rule%Qdof

    allocate( Fxj(1:Qdof, 1:2)  )

    ! nodes inside boundary of K
    call ComputeF(elem, Qdof, L_rule%lambda(1:Qdof, 1:2), Fxj(1:Qdof, 1:2) )

    allocate(u_exact(1:Qdof, 1:ndim, 0:2 ) )
    do j=1 , Qdof
       call Exact_Scalar(Fxj(j, 1:2), u_exact(j, 1:ndim, 0), state%time%ttime)
       call Der_Exact_Scalar(Fxj(j, 1:2),  u_exact(j, 1:ndim, 1:2), state%time%ttime)
    enddo


    !! DG test functions in xi

    ! DG test functions in Lagrange nodes
    allocate( phiL(1:Fdof, 1:Qdof), DphiLREF(1:Fdof, 1:2, 1:Qdof), &
         DphiL(1:Fdof, 1:2, 1:Qdof) ) 

    call PHI_orthonormal(Qdof, nbDim, L_rule%lambda(1:Qdof, 1:nbDim), 3, &
         Fdof, phiL(1:Fdof, 1:Qdof), DphiLREF(1:Fdof, 1:2, 1:Qdof) )

    do j=1, Fdof
       DphiL(j, 1, 1:Qdof) = &
            elem%F%D1F0(1,1) * DphiLREF(j, 1, 1:Qdof) &
            +elem%F%D1F0(1,2)* DphiLREF(j, 2, 1:Qdof)

       DphiL(j, 2, 1:Qdof) = &
            elem%F%D1F0(2,1) * DphiLREF(j, 1, 1:Qdof) &
            +elem%F%D1F0(2,2)* DphiLREF(j, 2, 1:Qdof)
    enddo


    ! fluxes in Lagrang nodes
    ! first orientation
    it= 0
    do j=0, L_rule%Qdeg
       do k=0, L_rule%Qdeg - j
          it = it + 1
          if(it > Qdof) stop "TROUB FRDWQKJVF 1"
          write(ifile,*) Fxj(it, 1:2), & 
               !dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 3)), &
               !dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 4)), &
               !dot_product( DphiL(1:dof,2,it), elem%RTNphi(1:dof, 3)) &  ! divergence of flux
               !-dot_product(DphiL(1:dof,1,it), elem%RTNphi(1:dof, 4)), &
               dot_product(phiL(1:Fdof,it), elem%RTNphi(1:Fdof, idx)), & ! reconstruction s_h
               dot_product(DphiL(1:Fdof,1, it), elem%RTNphi(1:Fdof, idx)), & ! reconstruction s_h_x
               dot_product(DphiL(1:Fdof,2, it), elem%RTNphi(1:Fdof, idx)), &   ! reconstruction s_h_y
               u_exact(it,1,0), u_exact(it,1,0)- dot_product(phiL(1:Fdof,it), elem%RTNphi(1:Fdof, idx)), &
               u_exact(it,1,1:2), &
               ((u_exact(it,1,1) - dot_product(DphiL(1:Fdof,1, it), elem%RTNphi(1:Fdof, idx)))**2 + &
               (u_exact(it,1,2) - dot_product(DphiL(1:Fdof,2, it), elem%RTNphi(1:Fdof, idx)))**2)**0.5
               
       enddo
       write(ifile,'(x)')
    enddo
    write(ifile,'(x)')

    ! second orientation
    do j=0, L_rule%Qdeg
       it = nface + 1 - j 
       do k=0, L_rule%Qdeg - j
          if(it > Qdof) stop "TROUB FRDWQKJVF 1"
          write(ifile,*) Fxj(it, 1:2), & 
               !!dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 3)), &
               !!dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 4)), &
               !!dot_product( DphiL(1:dof,2,it), elem%RTNphi(1:dof, 3)) &  ! divergence of flux
               !!-dot_product(DphiL(1:dof,1,it), elem%RTNphi(1:dof, 4)), &
               dot_product(phiL(1:Fdof,it), elem%RTNphi(1:Fdof, idx)), & ! reconstruction s_h
               dot_product(DphiL(1:Fdof,1, it), elem%RTNphi(1:Fdof, idx)), & ! reconstruction s_h_x
               dot_product(DphiL(1:Fdof,2, it), elem%RTNphi(1:Fdof, idx)), &   ! reconstruction s_h_y
               u_exact(it,1,0), u_exact(it,1,0)- dot_product(phiL(1:Fdof,it), elem%RTNphi(1:Fdof, idx))
          it = it + nface - k 
       enddo
       write(ifile,'(x)')
    enddo

    write(ifile,'(x)')

    ! third orientation
    do j=0, L_rule%Qdeg
       it = j + 1
       do k=0, L_rule%Qdeg - j
          if(it > Qdof) stop "TROUB FRDWQKJVF 1"
          write(ifile,*) Fxj(it, 1:2), & 
               !!dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 3)), &
               !!dot_product(phiL(1:dof,it), elem%RTNphi(1:dof, 4)), &
               !!dot_product( DphiL(1:dof,2,it), elem%RTNphi(1:dof, 3)) &  ! divergence of flux
               !!-dot_product(DphiL(1:dof,1,it), elem%RTNphi(1:dof, 4)), &
               dot_product(phiL(1:Fdof,it), elem%RTNphi(1:Fdof, idx)), & ! reconstruction s_h
               dot_product(DphiL(1:Fdof,1, it), elem%RTNphi(1:Fdof, idx)), & ! reconstruction s_h_x
               dot_product(DphiL(1:Fdof,2, it), elem%RTNphi(1:Fdof, idx)), &  ! reconstruction s_h_y
               u_exact(it,1,0), u_exact(it,1,0)-dot_product(phiL(1:Fdof,it), elem%RTNphi(1:Fdof, idx))
          it = it + nface +1 - k 
       enddo
       write(ifile,'(x)')
    enddo

    write(ifile,'(x)')



    deallocate(Fxj, u_exact, phiL,  DphiL, DphiLREF)

    !deallocate( xi, Fxi)
    !deallocate(phiT)

  end subroutine DRAW_reconstruction


  !> solving of a saddle point problem by Schur complements
  !  (  A    B ) ( x ) = ( f )
  !  ( -B^T  0 ) ( y ) = ( g )
  subroutine SchurComplements_ITER(n, m, A, B, ip, f, g, x, y  ) 
    integer, intent(in) :: n, m
    real, dimension(1:n, 1:n), intent(inout) :: A
    real, dimension(1:n, 1:m), intent(in) :: B
    integer, intent(in) :: ip   ! number of RHS
    real, dimension(1:n, 1:ip), intent(in) :: f
    real, dimension(1:m, 1:ip), intent(in) :: g
    real, dimension(1:n, 1:ip), intent(inout) :: x
    real, dimension(1:m, 1:ip), intent(inout) :: y
    real, dimension(:), allocatable :: u, rhs, dd
    integer:: restart = 30 !30    ! 50  ! GMRES restarted after 'restart' iterations  !45
    integer:: nloops = 25 !30  !5     ! 10   ! maximal number of restarted cycles !100	  !40
    real :: Mtol = 1E-12  
    real :: rezid 
    integer :: not_conv
    integer :: iter, iout = 0

    integer :: i, j, k,info, nsize
    logical :: iprint



    nsize = m + n
    allocate(u(1:nsize), rhs(1:nsize) , dd(1:nsize) )

    allocate(state%Schur_A(1:nsize, 1:nsize) )
    !allocate(state%Schur_B)

    state%Schur_A(1:n, 1:n) = A(1:n, 1:n)
    state%Schur_A(1:n, n+1: n+m) = B(1:n, 1:m)
    state%Schur_A(n+1 : n+m, 1: n) = - transpose( B(1:n, 1:m) )
    state%Schur_A(n+1 : n+m, n+1: n+m) = 0.

    do i=1, n
       dd(i) = state%Schur_A(i,i)
       state%Schur_A(i, 1:n) = state%Schur_A(i, 1:n) / dd(i)
    enddo
    
    

    do i = 1, ip
    
       u(1:n)    = x(1:n,i)
       u(n+1:n+m) = y(1:m, i)
       
       rhs(1:n)     = f(1:n, i) / dd(i)
       rhs(n+1:n+m) = g(1:m, i)
       
       
       call gmres(nsize, u, rhs, restart*nloops, 1.,  &
            bMVprod_SCHUR, bMVnull, restart,  Mtol, iout, iter, rezid, &
            not_conv)
       write(*,'(a10,i5, es14.6)' ) 'gMrES:', iter, rezid

        x(1:n ,i) = u(1:n)    
        y(1:m, i) = u(n+1:n+m) 

    enddo

    deallocate(u, rhs, state%Schur_A, dd)

    !stop 'end subroutine SchurComplements_ITER'

  end subroutine SchurComplements_ITER

end module neumann_estim

