!> estimates for AMA adaptation
module higher_order_local

  use main_data  ! contains type(mesh) :: grid for computation
  use problem_oper
  use euler_problem
  use estimates
  use plot_geom
  use eval_sol
  use ama_hp_interpol
  use AMA_estims

  implicit none

  public:: ComputeHO_LocalProblems
  public:: HO_LocalVertexProblem
  public:: Set_FE_Nodal_functions
  public:: CreateLocal_GlobalPairs
  public:: AssembGlobalMatrix

contains

  !> compute the error estimates by a higher order reconstruction
  subroutine ComputeHO_LocalProblems( )
    class(element), pointer :: elem
    logical, dimension(:), allocatable :: inner ! inner vertex
    integer, dimension(:), allocatable :: N  ! number of elements sharing a vertex
    integer, dimension(:,:,:), allocatable :: supp !list of corresponding elements
    real, dimension(:,:), allocatable :: regul !list of corresponding elements
    integer :: i, j, k, dof
    integer :: maxdeg = 20
    logical :: loc_implicitly
    character(len=15) :: file1, file2
    character(len=5) :: ch5

    if(nbDim /=2 ) then
       print*,' ## ComputeHigherOrderEstims only for nbDim = 2!'
       stop
    endif

    print*,' # starting of  ComputeHigherLocalProblems( )'
    ! create a list of elements sharing a vertex
    allocate( inner(1:grid%npoin), N(1:grid%npoin))
    allocate( supp(1:grid%npoin, 1:maxdeg, 1:2) )

    call SeekVertexSupports(grid, maxdeg, inner(1:grid%npoin), N(1:grid%npoin), &
         supp(1:grid%npoin, 1:maxdeg,1:2) )


    ! create the list of elements sharing at least a vertex with elem
    call SeekElemSupports(grid)  

    ! setting of the array elem%wS
    do i=1,grid%nelem
       elem => grid%elem(i)
       !dof = elem%dof
       dof = DOFtriang( MaxDegreeImplemented ) !???

       allocate(elem%wS( 1:ndim, 1:dof ) )
       elem%wS(:,:) = 0.

       elem%HO_deg = 0  ! degree of the HO reconstruction
       elem%eta(:,:) = 0.
    enddo

    state%estim( : , :) = 0.
    do i=1,grid%npoin
    !do i= 21, 21
    !do i= 1, 1
    !do i= 12, 12
       call HO_LocalVertexProblem(1, i, inner(i), N(i), supp(i, 1:N(i),1:2)) 
    enddo
    

    do i=1,grid%nelem
       elem => grid%elem(i)

       ! the final solution 
       !call PlotElemFunction3D(3000 , elem, elem%dof + elem%deg+ 2, &
       !     elem%wS(1, 1:elem%dof + elem%deg+ 2) )

       !call PlotElem_D_Function3D(3001 , elem, elem%dof + elem%deg+ 2, &
       !     elem%wS(1, 1:elem%dof + elem%deg+ 2) )


       call Set_EE_HO_LocalProblem(elem)
 
       state%estim(1:max_eta, : ) = state%estim(1:max_eta, : ) + elem%eta(1: max_eta, : ) 
 
       deallocate(elem%wS)

    enddo

    write(*,'(a25, 2(3es12.4, a2))') 'HO_rec error estims L2', & 
         state%err(L2), state%estim(HO_estim_L2_p2, :)**0.5, state%estim(HO_trunc_L2_p2, :)**0.5,'|',&
         state%err(H1), state%estim(HO_estim_H1_p2, :)**0.5, state%estim(HO_trunc_H1_p2, :)**0.5

    deallocate(inner, N, supp)

  end subroutine ComputeHO_LocalProblems


  !> solve the local Neumann problem on the patch corresponding to a vortex of the mesh
  subroutine HO_LocalVertexProblem(i_var, ip, inner, N, supp )
    integer, intent(in) :: i_var   ! index of variant
    integer, intent(in) :: ip      ! index of node
    logical, intent(in) :: inner   ! inner vertex?
    integer, intent(in) :: N   ! number of elememnts sharing vertex
    integer, dimension(1:N,1:2), intent(in) :: supp   !idx of elems sharing vertex
    class(element), pointer :: elem
    type(volume_rule), pointer :: V_rule
    type(Lagrang_rule), pointer :: L_rule
    real, dimension(:, :, :), allocatable :: phi
    real, dimension(:, :), allocatable :: phi_loc
    real, dimension(:, :), allocatable :: lambda
    real, dimension(:), allocatable :: weights
    real, dimension(:, :, :), allocatable :: Dphi, DL_phi
    real, dimension(:, :, :), allocatable :: Dwi
    real, dimension( :, :), allocatable :: Dhat
    real, dimension( :, :), allocatable :: ww
    real, dimension( :, :, :), allocatable :: w_BC
    real, dimension( :, :, :), allocatable :: Dwi_BC
    real, dimension(:, :), allocatable :: A_loc, A
    integer, dimension(:, :), allocatable :: itrans
    real :: val
    integer :: i, ie, ie1, ie2, deg, deg1, dof, Qnum, Qdof, dN, dofE
    integer :: j, k, l, l1, l2

    !real, dimension(:, :), allocatable :: Fx

    ! setting of local dof
    deg = maxval(grid%elem(supp(1:N,1))%deg) ! setting of p of reconstruction
    deg1 = deg + 1
    
    dof = DOFtriang( deg1 )  

    L_rule => state%space%L_rule(deg1)


    ! setting of quadrature
    Qnum = state%space%Qdeg(deg1, 1)
    V_rule => state%space%V_rule(Qnum)

    Qdof = V_rule%Qdof
    !allocate( Fx(1:Qdof, 1:nbDim) )
 

    allocate( itrans( 1:N, 1:dof) ) ! pair of the local and global indexes
    call CreateLocal_GlobalPairs(deg1, dof, N, itrans(1:N, 1:dof), inner, dN)


    allocate( phi_loc(1:dof, 1:dof) ) ! local basis functions in Lagrange nodes in correct order
    allocate( phi(1:N, 1:dof, 1:dof) )     ! local basis functions in basis coefficients

    allocate( Dwi( 0:nbDim, 1:ndim, 1:Qdof) ) ! solution and its derivatives in integ nodes

    allocate( w_BC( 1:N, 1:ndim, 1:dof) ) ! BC terms in basis coefficients
    allocate( Dwi_BC( 1:ndim, 0:nbDim, 1:Qdof) ) ! BC terms in integ nodes

    allocate( Dhat( 0:4, 1:Qdof) ) ! the hat functions and its derivatives in integ nodes

    allocate( lambda(1:dof, 1:3) )  ! barycentric coordinates of Lagrange nodes in correct order

    allocate( Dphi(1:dof, 1:nbDim, 1:Qdof) )   ! DG basis functions in integ nodes

    allocate( DL_phi(1:dof, 0:nbDim, 1:Qdof) ) ! Lagrange basis functions in integ nodes

    allocate( weights( 1:Qdof) ) ! weights for integration


    allocate( ww( 1:ndim, 1:dof) ) ! solution in DG basis functions

    !allocate( A_loc( 1:dof, 1:dof + ndim) ) ! local stiff matrix

    allocate( A( 1:dN, 1:dN + ndim) ) ! global stiff matrix
    A(:,:) = 0.

    !write(20,*) grid%x(ip, :)
    do i=1, N
       ie  = supp(i,2)   ! inner index of the vertex

       elem => grid%elem(supp(i,1))
       

       !write(*,'(a8,5i5,2es12.4)' ) 'elemW:',ip, i, ie, elem%i,dN, elem%xc
       !write(*,'(2(a8,i5))' ) 'deg = ',deg1,',   dof = ', dof
       !write(21,*) elem%xc

       ! setting of the barycentric coordinates of the Lagrange nodes in the correct order
       call  Set_FE_Nodal_functions( deg1, dof, ie,  lambda(1:dof, 1:3) )
       call  Eval_L_Direct(L_rule, dof, lambda(1:dof, 1:2), phi_loc(1:dof, 1:dof))

       !do k=1,dof
       !   write(*,'(a6, i3, 3es12.4,a1,300es12.4)') 'lam:2w',k, lambda(k, :),'|', phi_loc(1:dof, k)
       !enddo


       ! evaluation of the Lagrangian basis functions (in the correct order) in DG basis functions
       do j=1, dof
          call Lagr2BasisDiffOne(elem, deg1,  phi_loc(1:dof, j), Qnum, dof,  phi(i, j, 1:dof) )
          !write(*,'(a8,5i5,200es12.4)' ) 'phi:',ip, i, ie, elem%i,j, phi(j, :)


          !!call PlotElemFunction3D(100*i+j, elem, dof, phi(i, j,1:dof) )

       enddo

       ! derivatives of DG basis functions in integ nodes 
       call Eval_Dphi_plus(elem, V_rule, dof,  Dphi(1:dof, 1:nbDim, 1:Qdof) )


       ! Lagrange basis functions and its derivatives in integ nodes
       DL_phi(1:dof, 0, 1:Qdof) = matmul(phi(i, 1:dof, 1:dof), V_rule%phi(1:dof, 1:Qdof) )
       DL_phi(1:dof, 1, 1:Qdof) = matmul(phi(i, 1:dof, 1:dof), Dphi(1:dof, 1, 1:Qdof) )
       DL_phi(1:dof, 2, 1:Qdof) = matmul(phi(i, 1:dof, 1:dof), Dphi(1:dof, 2, 1:Qdof) )

       !  weights in integ nodes
       call Eval_V_Weights_plus(elem, V_rule, weights(1:Qdof) )
       
       ! approximate solution and its derivatives in integ nodes
       call  Eval_w_Elem_plus(elem, V_rule, Dwi(0:2, 1:ndim, 1:Qdof) )

       ! TAKES VALUES from HO RECONSTRUCTION
       !Dwi(0, 1:ndim, 1:Qdof) = matmul(elem%wSS(1:ndim, 1:dof, 0), V_rule%phi(1:dof, 1:Qdof) )
       !Dwi(1, 1:ndim, 1:Qdof) = matmul(elem%wSS(1:ndim, 1:dof, 0), Dphi(1:dof, 1, 1:Qdof) )
       !Dwi(2, 1:ndim, 1:Qdof) = matmul(elem%wSS(1:ndim, 1:dof, 0), Dphi(1:dof, 2, 1:Qdof) )

       ! variant polynomially preserving
       !Dwi(0, 1:ndim, 1:Qdof) = matmul(elem%wSS(1:ndim, 1:dof, 0), V_rule%phi(1:dof, 1:Qdof) )
       !Dwi(1, 1:ndim, 1:Qdof) = matmul(elem%wSS(1:ndim, 1:dof, 1), V_rule%phi(1:dof, 1:Qdof) )
       !Dwi(2, 1:ndim, 1:Qdof) = matmul(elem%wSS(1:ndim, 1:dof, 2), V_rule%phi(1:dof, 1:Qdof) )

       !!call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof,1:nbDim),Fx(1:Qdof, 1:nbDim) )

       ! the hat function
       ie1 = mod(  ie, 3) + 1
       ie2 = mod( ie1, 3) + 1

       ! the hat function and its (constant) gradient on elem
       Dhat(0, 1:Qdof) = V_rule%lambda(1:Qdof, ie2) 
       Dhat(1, :) = (grid%x(elem%face(idx, ie1), 2) - grid%x(elem%face(idx, ie2), 2) ) /elem%area/2
       Dhat(2, :) = (grid%x(elem%face(idx, ie2), 1) - grid%x(elem%face(idx, ie1), 1) ) /elem%area/2
       
       
       ! the functions used for the computation
       !do j=1,Qdof
       !   write(1000, *) Fx(j, 1:2), Dwi(0:2, 1, j) 
       !   write(1001, *) Fx(j, 1:2), Dhat(0:2,  j) 
       !   write(1002, *) Fx(j, 1:2), Dhat(0:2,  j) * Dwi(0:2, 1, j) 
       !enddo



       ! the setting of the stiff matrix
       do l1=1,dof
          do l2=1, l1  !dof  !only symmetric part
             k = itrans(i, l1)
             l = itrans(i, l2)
             
             if( k > 0 .and. l > 0 ) then
                val = dot_product(weights(1:Qdof), &
                     DL_phi(l1, 1, 1:Qdof) * DL_phi(l2, 1, 1:Qdof) +&
                     DL_phi(l1, 2, 1:Qdof) * DL_phi(l2, 2, 1:Qdof)  )

                ! direct setting of the global matrix
                A( k, l) = A( k, l) + val
                if( l /= k) A( l, k) = A( l, k) + val


                !A_loc(l1, l2) = dot_product(weights(1:Qdof), &
                !     DL_phi(l1, 1, 1:Qdof) * DL_phi(l2, 1, 1:Qdof) +&
                !     DL_phi(l1, 2, 1:Qdof) * DL_phi(l2, 2, 1:Qdof)  )
                !A_loc(l2, l1) = A_loc(l1, l2)
             endif

          enddo
          
          ! RHS
          if(k > 0) then
             ! RHS = D( hat * u_h) =  D hat * u_h + hat * D u_h
             do j = 1, ndim
                Dhat(3, 1:Qdof) = Dhat(1,1:Qdof) * Dwi(0,j,1:Qdof) +  Dhat(0,1:Qdof) * Dwi(1,j,1:Qdof)
                Dhat(4, 1:Qdof) = Dhat(2,1:Qdof) * Dwi(0,j,1:Qdof) +  Dhat(0,1:Qdof) * Dwi(2,j,1:Qdof)

                val =  dot_product(weights(1:Qdof), &
                     DL_phi(l1, 1, 1:Qdof) * Dhat(3, 1:Qdof) +&
                     DL_phi(l1, 2, 1:Qdof) * Dhat(4, 1:Qdof)  )

                A( k, dN+j ) = A( k, dN+j) + val
                
                !A_loc(l1, dof+j) = val
             enddo
          endif

       enddo
       
       ! adding of the BC
       if(.not. inner  ) then 

          !   BC functions in FE basis coefficients
          call Set_BC_LocalvertexProblem(elem, i, ie, N, deg1, dof, w_BC(i, 1:ndim, 1:dof) )
          Dwi_BC(:,:, :) = 0.

          ! BC FE functions in integ nodes
          do j=1, ndim
             do l=1, dof
                if(w_BC(i, j, l) /= 0.) then
                   !!write(*,'(a6, 2i5, 3es12.4)') '...BC:',i,l,w_BC(i, j, l)
                   Dwi_BC(j, 0:2, 1:Qdof) = Dwi_BC(j, 0:2,1:Qdof) + w_BC(i,j,l) * DL_phi(l, 0:2, 1:Qdof)
                endif
             enddo

             ! the BC function used for the computation
             !do l=1,Qdof
             !   write(1008, *) Fx(l, 1:2), Dwi_BC(j, 0:2, l) 
             !enddo

             ! adding of BC terms to RHS
             do l1=1,dof
                k = itrans(i, l1)
                if(k > 0) then
                   val = dot_product(weights(1:Qdof), &
                        DL_phi(l1, 1, 1:Qdof) * Dwi_BC(j, 1, 1:Qdof) +&
                        DL_phi(l1, 2, 1:Qdof) * Dwi_BC(j, 2, 1:Qdof)  )
                   
                   A( k, dN+j ) = A( k, dN+j) - val
                         
                   !!write(*,'(a8,3i5, 300es12.4)') 'BC BC BC:',l1, k, j, val

                endif
             enddo
          enddo
       endif
       !do k=1,dof
       !   write(*,'(a8,i5, 300es12.4)') 'local:',k, A_loc(k, :)
       !enddo


       !call AssembGlobalMatrix(dN, dof, A(1:dN, 1:dN+1), A_loc(1:dof, 1:dof+1),itrans(i, 1: dof) )


      ! if(i == N) stop 'ede3223s21s2'
       
    enddo  ! do i=1,N

    !do k=1,dN
    !   write(*,'(a8,i5, 300es12.4)') 'Stiff:',k, A(k, :)
    !enddo
    
    ! solution of the local problem
    call SolveLocalMatrixProblem(dN, A(1:dN, 1:dN), ndim, A(1:dN, dN+1 : dN+ndim) ) 

    !print*,'------------------------------------------'
    !do k=1,dN
    !   write(*,'(a8,i5, 300es12.4)') 'Stiff:',k, A(k, :)
    !enddo

    ! setting of the solution to particular elements
    do i=1, N
       elem => grid%elem(supp(i,1))

       ww(:,:) = 0.
       do l1=1,dof
          k = itrans(i, l1)
          
          if(k > 0) then  ! setting the DG functions from the FE ones
             do j=1, ndim
                ww(j, 1:dof) = ww(j, 1:dof) + phi(i, l1, 1:dof) * A(k, dN + j)
             enddo
          endif
       enddo

       ! adding of BC
       if(.not. inner ) then 
          do j=1, ndim
             do l=1, dof
                if(w_BC(i, j, l) /= 0.) then
                   !write(*,'(a6, 2i5, 3es12.4)') '!!!BC:',i,l,w_BC(i, j, l)
                   ww(j, 1:dof) = ww(j, 1:dof) + phi(i, l, 1:dof) * w_BC(i, j, l)
                endif
             enddo
          enddo
       endif

       ! solution of the local problem
       !call PlotElemFunction3D(1000 + 100*i, elem, dof, ww(1, 1:dof) )

       elem%HO_deg = max(elem%HO_deg, deg1) 
       do j=1,ndim
          elem%wS(j,1:dof) = elem%wS(j, 1:dof) + ww(j, 1:dof) 
       enddo


    enddo


    deallocate( phi, phi_loc, Dphi, DL_phi, itrans, Dwi, w_BC, Dwi_BC, Dhat, A)

  end subroutine HO_LocalVertexProblem

  !>  setting of barycentric coordinates of the FE nodal basis functions in the correct order
  subroutine Set_FE_Nodal_functions( deg, dof, ie, lambda )
    integer, intent (in) :: deg, dof
    integer, intent (in) :: ie  ! inner index of the central vertex node
    real, dimension(1:dof, 1:3), intent(inout) :: lambda
    integer :: k, j, l1, l2, ie0, ie1, ie2

    ie0 = ie
    ie1 = mod(ie0, 3) + 1
    ie2 = mod(ie1, 3) + 1

    ! vertex FE basis function
    lambda(1, ie0) = 0. ; lambda(1, ie1) = 0. ; lambda(1, ie2) = 1. 
    lambda(2, ie0) = 1. ; lambda(2, ie1) = 0. ; lambda(2, ie2) = 0. 
    lambda(3, ie0) = 0. ; lambda(3, ie1) = 1. ; lambda(3, ie2) = 0. 

    ! edge FE basis fucntions
    if(deg > 1) then
       k = 3

       do j=1, 3   ! we go over all three edges
          do l1= 1, deg-1  ! nodes for each edge
             k = k + 1

             lambda(k, ie0) = 1.* l1 / deg
             lambda(k, ie1) = 0.
             lambda(k, ie2) = 1. - lambda(k, ie0)  - lambda(k, ie1) 
          enddo
          ie0 = ie1
          ie1 = mod(ie0, 3) + 1
          ie2 = mod(ie1, 3) + 1
       enddo

    endif
             
    ! interior  FE basis fucntions, all Lag nodes exept the boundaries
    if(deg > 2) then

       ie0 = ie
       ie1 = mod(ie0, 3) + 1
       ie2 = mod(ie1, 3) + 1

       do l1 = 1, deg - 1
          do l2 = 1, deg -1 - l1
             k = k + 1
             lambda(k, ie0)  = 1.* l2 / deg
             lambda(k, ie1) = 1.* l1 / deg
             lambda(k, ie2) = 1. - lambda(k, ie)  - lambda(k, ie1) 
          enddo
       enddo
    endif


    if(k /= dof) stop 'inconsistency in Set_FE_Nodal_functions in the module higher_order_local'

  end subroutine Set_FE_Nodal_functions


  !> create the pair of FE DOF on elem (local index) and the patch \f$\omega_a \f$ (global index)
  !> itrans(i,k) : i element index, k, local FE index, itrans(i,k) global FE index
  !> moreover, it set dN - size of the local problem
  subroutine CreateLocal_GlobalPairs(deg, dof, N, itrans, inner, dN)
    integer, intent(in) :: deg   ! deg of recontruction
    integer, intent(in) :: dof   ! number of RTN DOF on one element
    integer, intent(in) :: N  ! number of elements in supp around the vertex 
    integer, dimension(1:N, 1:dof), intent(inout) :: itrans  ! 3rd index =1 matrix, =2 BC
    logical, intent(in) :: inner !  inner node
    integer, intent(inout) :: dN  ! global DOF
    integer :: i, j, k, l, l1, l2,  N_vertex, N_elem, l_BC
    
    itrans(:, :) = 0

    N_elem = (deg - 1)  ! DOF for one element
    if(deg > 2) N_elem = N_elem + (deg-2)*(deg-1) / 2
 

    if(inner) then
       itrans(1:N, 1) = 1   ! only the central one
       N_vertex = 1
       
    else ! boundary node, no vertex FE basis
       N_vertex = - (deg -1)  ! the first edge has no DOF
       
    endif

    l = N_vertex
    ! edge and volume FE basis functions
    do i=1,N 
       
       ! edge FE basis functions
       if(deg > 1) then
          
          ! first edge  j = 1
          if(inner .or. i > 1) then 
             k = 3                          ! local index
             l = N_vertex + (i-1)*N_elem    ! global index
             
             do l1= 1, deg-1  ! nodes for each edge
                k = k + 1
                l = l + 1
                itrans(i, k) = l
             enddo
          endif
          

          ! volume FE functions
          if(deg > 2) then  
             k = 3 + 3*(deg -1)                       ! local index
             l = N_vertex + (i-1)*N_elem + deg - 1    ! global index
             do l1 = 1, deg - 1
                do l2 = 1, deg -1 - l1
                   l = l + 1
                   k = k + 1
                   itrans(i, k) =  l
                enddo
             enddo
             
          endif
          
          if(inner .or. i < N) then
             ! third edge  j = 1
             k = 3 + 2*(deg - 1)                              ! local index
             l = N_vertex + (i-1)*N_elem +  2*(deg - 1)  + 1  ! global index
             if(deg > 2) l = l + (deg -1) * ( deg - 2) / 2
             
             if(inner .and. i == N) l = N_vertex + deg  ! the last edge is the first one
             
             do l1= 1, deg-1  ! nodes for each edge
                k = k + 1
                l = l - 1     ! opposite orientation
                itrans(i, k) = l
             enddo
          endif

          
       endif
       
       !if(i == 1) then
       !   write(*,'(a4,l2,i5,a2,30i5)') 'CTF',inner, i,'|', 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20
       !print*
       !endif

       !write(*,'(a6,l2,i5,a2,30i5)') 'CTF 1',inner, i,'|', itrans(i, 1:dof)

    enddo ! do i=1,N

    !write(*,*) 'dN = ', dN

    dN = maxval(itrans(1:N, 1:dof) )
    !write(*,*)  'dN = ', dN
    

  end subroutine CreateLocal_GlobalPairs

  !> assembling of the global matrix by the local blocks
  subroutine  AssembGlobalMatrix(dN, dof, A, A_loc, itrans )
    integer, intent(in) :: dN, dof
    real, dimension(1:dN, 1:dN+1), intent(inout) :: A
    real, dimension(1:dof, 1:dof+1), intent(inout) :: A_loc
    integer, dimension(1: dof), intent(in) ::  itrans
    integer :: i, j, k, l

    do i=1, dof
       k = itrans(i)
       if(k > 0 ) then

          do j=1,dof
             l = itrans(j)

             if(k > dN .or. l > dN) &
                  stop 'trouble detdxe53dhe3 in ho-local.f90, subroutine AssembGlobalMatrix'

             if(l > 0)  A(k,l) = A(k,l) + A_loc(i,j)

          end do
          ! RHS
          A(k,dN+1) = A(k,dN+1) + A_loc(i,dof+1)
       endif
    enddo


  end subroutine AssembGlobalMatrix

  !> settingf of the BC for the local vertex problem, array w_BC contains the FE basis 
  !> coefficients
  subroutine Set_BC_LocalvertexProblem(elem, i, ie, N, deg1, dof, w_BC )
    class(element), target, intent(in) :: elem
    integer :: i, ie, N, deg1, dof
    real, dimension(1:ndim, 1:dof), intent(inout) :: w_BC
    real, dimension(:), allocatable :: xi
    integer :: l, k, ie1, ie2
    real:: rat

    allocate(xi(1:2) )
    w_BC(:,:) = 0.

    if(state%modelName /= 'scalar' .and. state%modelName /= '2eqs' ) then
       print*, 'BC for this model non-implemented de5ded8edhe'
       stop
    endif

    if(i == 1) then ! the first triangle
       
       ! vertex
       k = 1 ! the ordering of the vertexes
       xi(1:2) = grid%x(elem%face(idx, ie), 1:2)
       call Exact_Scalar(xi(1:nbDim), w_BC(1:ndim, k), state%time%ctime )

       ! the first edge
       if(deg1 > 1) then
          ie1 = mod(ie, 3) + 1
          do l=1, deg1-1
             k = 3 + l ! the ordering of the vertexes
             rat = 1. * l / deg1

             xi(1:2) = (1- rat) * grid%x(elem%face(idx, ie), 1:2) &
                  +   rat  *  grid%x(elem%face(idx, ie1), 1:2)

             call Exact_Scalar(xi(1:nbDim), w_BC(1:ndim, k), state%time%ctime )

             w_BC(1:ndim, k) = w_BC(1:ndim, k) * (1. - rat)  ! multiplication by the hat function

          enddo
       endif

    elseif(i == N) then ! the last triangle

       ! vertex
       k = 1 ! the ordering of the vertexes
       xi(1:2) = grid%x(elem%face(idx, ie), 1:2)
       call Exact_Scalar(xi(1:nbDim), w_BC(1:ndim, k), state%time%ctime )

       ! the last edge
       if(deg1 > 1) then
          ie1 = mod(ie, 3) + 1
          ie2 = mod(ie1, 3) + 1
          do l=1, deg1-1
             k = 3 + (deg1 - 1)*3 - l +1 ! the ordering of the vertexes
             rat = 1. * l / deg1

             xi(1:2) = (1- rat) * grid%x(elem%face(idx, ie), 1:2) &
                  +   rat  *  grid%x(elem%face(idx, ie2), 1:2)

             call Exact_Scalar(xi(1:nbDim), w_BC(1:ndim, k), state%time%ctime )

             w_BC(1:ndim, k) = w_BC(1:ndim, k) * (1. - rat)  ! multiplication by the hat function

          enddo
       endif

    else ! not first and not last element, only firts vertex function
       k = 1 ! the ordering of the vertexes
       xi(1:2) = grid%x(elem%face(idx, ie), 1:2)
       call Exact_Scalar(xi(1:nbDim), w_BC(1:ndim, k), state%time%ctime )

    endif


    deallocate(xi)
  end subroutine Set_BC_LocalvertexProblem


  !> evaluate the error estimate based on the solution of the HO local problem
  subroutine Set_EE_HO_LocalProblem(elem)
    class(element), target, intent(inout) :: elem
    type(volume_rule), pointer :: V_rule
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:,:), allocatable :: Dphi
    real, dimension(:,:,:), allocatable :: wi  
    real, dimension(:,:,:), allocatable :: DwExact  ! exact solution and its deriv in integ nds
    real, dimension(:), allocatable :: weights
  
    integer :: degP, dofP, dof, Qnum, Qdof
    integer :: iff, k

    dof = elem%dof   ! dof of the approximate solution 

    degP = elem%HO_deg
    dofP = DOFtriang( degP )  ! dof of the HO reconstruction

    Qnum = state%space%Qdeg(degP, 1)
    V_rule => state%space%V_rule(Qnum)

    Qdof = V_rule%Qdof

   ! quadrature weights
    allocate(weights(1:Qdof) )
    call Eval_V_Weights_plus(elem, V_rule, weights(1:Qdof))

    phi => V_rule%phi(1:dofP, 1:Qdof)

    allocate( Dphi(1:dofP, 1:nbDim, 1:Qdof) )
    call Eval_Dphi_plus(elem, V_rule, dofP,  Dphi(1:dofP, 1:nbDim, 1:Qdof) )


   ! exact solution
    allocate(DwExact(1:Qdof, 1:ndim, 0:nbDim))  ! 3rd index: 0 solution, 1-2 its gradient

    ! exact solution in integ nodes on subelements
    call SetExactSolutionArbitraryNodes(elem, Qdof, V_rule%lambda(1:Qdof, 1:2),  &
         DwExact(1:Qdof, 1:ndim, 0), DwExact(1:Qdof, 1:ndim, 1:nbDim))

    ! projection  and its derivatives  in integ nodes
    allocate(wi(1:15, 1:ndim, 1:Qdof) )


    ! recovery in integ nodes on elements
    do k=1, ndim
       wi(1, k, 1:Qdof) = matmul(elem%wS( k, 1:dofP ),  phi(1:dofP, 1:Qdof) )
       wi(2, k, 1:Qdof) = matmul(elem%wS( k, 1:dofP ), Dphi(1:dofP, 1, 1:Qdof) )
       wi(3, k, 1:Qdof) = matmul(elem%wS( k, 1:dofP ), Dphi(1:dofP, 2, 1:Qdof) )
    enddo

    ! approximate solution 
    do k=1, ndim
       wi(4, k, 1:Qdof) = matmul(elem%w(0, (k-1)*dof+1 : k*dof ),  phi(1:dof, 1:Qdof) )
       wi(5, k, 1:Qdof) = matmul(elem%w(0, (k-1)*dof+1 : k*dof ), Dphi(1:dof, 1, 1:Qdof) )
       wi(6, k, 1:Qdof) = matmul(elem%w(0, (k-1)*dof+1 : k*dof ), Dphi(1:dof, 2, 1:Qdof) )
    enddo

    !call PlotElemFunction3D(100, elem, dof, elem%w(0, 1:dof) )
    call PlotElemFunction3D(200, elem, dofP, elem%wS(1, 1:dofP) )

    

     iff = HO_estim_L2_p2 - 1
     
     do k=1,ndim
        ! error estimate in the L2 norm
        elem%eta(iff+1, k) =  & !elem%eta(iff+1, k) 
             + dot_product(weights(1:Qdof),  &
             ( wi(1, k, 1:Qdof)- wi(4, k, 1:Qdof))**2 )

        ! error estimate in the H1 semi-norm
        elem%eta(iff+2, k) =  & ! elem%eta(iff+2, k) 
             + dot_product(weights(1:Qdof),  &
             ( wi(2, k, 1:Qdof)- wi(5, k, 1:Qdof))**2 + ( wi(3, k, 1:Qdof)- wi(6, k, 1:Qdof))**2 )

        !!if(elem%i == 1) write(*,'(a10,6i5)') 'HO  iff =',iff+1, iff+2

        ! truncation error  in the L2 norm
        elem%eta(iff+3, k) = & ! elem%eta(iff+3, k) 
             + dot_product(weights(1:Qdof),  &
             ( DwExact(1:Qdof, k, 0)- wi(1, k, 1:Qdof))**2 )

        ! truncation error in the H1 semi-norm
        elem%eta(iff+4, k) = & ! elem%eta(iff+4, k) 
             + dot_product(weights(1:Qdof),  &
             (wi(2, k, 1:Qdof)- DwExact(1:Qdof, k,1) )**2 + (wi(3, k, 1:Qdof)- DwExact(1:Qdof,k,2))**2)
     enddo

     !if(elem%i == 1) then
     !   write(*,'(a6,2i5,40es12.4)' )'w:',elem%i, dof, elem%w( 0, 1:dof )
     !   write(*,'(a6,2i5,40es12.4)' )'wS:',elem%i, dofP, elem%wS( 1, 1:dofP )
     !   print*
     !endif

    deallocate(Dphi, DwExact, wi, weights)

  end subroutine Set_EE_HO_LocalProblem

end module higher_order_local
