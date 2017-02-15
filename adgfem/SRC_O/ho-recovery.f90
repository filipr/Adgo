!> hogher order recovery on the vertex patches 
module higher_order_recovery

  use main_data  ! contains type(mesh) :: grid for computation
  use problem_oper
  use euler_problem
  use estimates
  use plot_geom
  use eval_sol
  use ama_hp_interpol
  use higher_order_local
  use eval_sol
  !use AMA_estims

  implicit none

  public:: ComputeHO_Recovery
  public:: LocalVertexProblemRecovery
  public:: Set_FE_Nodal_functionsREDO
  public:: CreateLocal_GlobalPairsREDO
  public:: AssembGlobalMatrixREDO

contains

  !> compute the error estimates by a higher order recovery of the vertex patches
  subroutine ComputeHO_Recovery( )
    class(element), pointer :: elem
    logical, dimension(:), allocatable :: inner ! inner vertex
    integer, dimension(:), allocatable :: N  ! number of elements sharing a vertex
    integer, dimension(:,:,:), allocatable :: supp !list of corresponding elements
    real, dimension(:,:), allocatable :: regul !list of corresponding elements
    integer :: i, j, k, dof, ifile, ifile1
    integer :: maxdeg = 20
    logical :: loc_implicitly
    character(len=15) :: file1, file2
    character(len=5) :: ch5

    if(nbDim /=2 ) then
       print*,' ## ComputeHO_Recovery only for nbDim = 2!'
       stop
    endif

    print*,' # starting of  ComputeHigherLocalProblems( )'
    ! create a list of elements sharing a vertex
    allocate( inner(1:grid%npoin), N(1:grid%npoin))
    allocate( supp(1:grid%npoin, 1:maxdeg, 1:2) )

    call SeekVertexSupports(grid, maxdeg, inner(1:grid%npoin), N(1:grid%npoin), &
         supp(1:grid%npoin, 1:maxdeg,1:2) )


    ! IS IT NECESSARY?
    ! create the list of elements sharing at least a vertex with elem
    if(.not. grid%ElemSupports) call SeekElemSupports(grid)  

    ! setting of the array elem%wS
    do i=1,grid%nelem
       elem => grid%elem(i)
       !dof = elem%dof
       dof = DOFtriang( MaxDegreeImplemented ) !???

       allocate(elem%wS( 1:ndim, 1:dof ) )
       elem%wS(:,:) = 0.

       elem%HO_deg = 0  ! degree of the HO reconstruction
       !A!elem%eta(:,:) = 0.
    enddo

    !!state%estim( : , :) = 0.

    do i=1,grid%npoin
    !do i= 21, 21
    !do i= 1, 1
    !do i= 12, 12
       call LocalVertexProblemRecovery(1, i, inner(i), N(i), supp(i, 1:N(i),1:2)) 
    enddo
    

    do i=1,grid%nelem
       elem => grid%elem(i)

       ifile = 400 + state%space%adapt%adapt_level
       ifile1 = 500 + state%space%adapt%adapt_level
       ! the final solution 
       !call PlotElemFunction3D(ifile , elem, elem%dof + elem%deg+ 2, &
       !     elem%wS(1, 1:elem%dof + elem%deg+ 2) )

       !if(elem%face(idx, 1) == 13 .or. &
       !     elem%face(idx, 2) == 13 .or. &
       !     elem%face(idx, 3) == 13 ) &

       !dof = DOFtriang(elem%HO_deg)
       !call PlotElem_D_Function3D(ifile , elem, dof,  elem%wS(1, 1:dof) ) 

       !!!elem%wS(1, 1:elem%dof + elem%deg+ 2) = 0.
       !elem%wS(1, 1:elem%dof ) = elem%wS(1, 1:elem%dof ) - elem%w(0, 1:elem%dof) 

       !call PlotElem_D_Function3D(ifile1, elem, dof,  elem%wS(1, 1:dof) ) 
       !elem%wS(1, 1:elem%dof ) = elem%wS(1, 1:elem%dof ) + elem%w(0, 1:elem%dof) 


       !A!call Set_EE_HO_LocalProblemREDO(elem)
 
       !A!state%estim(1:max_eta, : ) = state%estim(1:max_eta, : ) + elem%eta(1: max_eta, : ) 
 

    enddo

    !A!write(*,'(a25, 2(3es12.4, a2))') 'HO_rec error estims L2', & 
    !A!     state%err(L2), state%estim(HO_estim_L2_p2, :)**0.5, state%estim(HO_trunc_L2_p2, :)**0.5,'|',&
    !A!     state%err(H1), state%estim(HO_estim_H1_p2, :)**0.5, state%estim(HO_trunc_H1_p2, :)**0.5


    ! computation of the higher order derivatives 
    !call Compute_means_values_DW(grid, ndim, 2)



    !do i= 21, 21
    !do i= 13, 13   ! corner
    !do i= 8, 8  ! inner

    !do i= 8, 13, 5  ! inner
    do i=1,grid%npoin
       call LocalPatchProjections( i, inner(i), N(i), supp(i, 1:N(i),1:2)) 
    enddo
 
    deallocate(inner, N, supp)


    do i=1,grid%nelem
       elem => grid%elem(i)

       deallocate(elem%wS)
       !deallocate(elem%wSD)
    enddo


  
    !stop "end subroutine ComputeHO_RECOVERY  8438de"

  end subroutine ComputeHO_Recovery


  !> solve the local Neumann problem on the patch corresponding to a vortex of the mesh
  subroutine LocalVertexProblemRecovery(i_var, ip, inner, N, supp )
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
    logical :: iprint
    real, dimension(:, :), allocatable :: Fx

    iprint = .false.
    !if(ip ==  8) iprint = .true.   !inner    '../Grids/LL-shape.uns.100.grid'
    !if(ip ==  3) iprint = .true.   !boundary  '../Grids/LL-shape.uns.100.grid'
    !if(ip ==  4) iprint = .true.   !boundary  '../Grids/LL-shape.uns.100.grid'
    !if(ip == 13) iprint = .true.   !inner corner  '../Grids/LL-shape.uns.100.grid'

    !if(ip == 15) iprint = .true.   !inner corner  '../Grids/LL-shape.uns.100.grid'
    !if(ip == 27) iprint = .true.   !inner corner  '../Grids/LL-shape.uns.100.grid'
    !if(ip == 55) iprint = .true.   !inner corner  '../Grids/LL-shape.uns.100.grid'
    !if(ip ==  4) iprint = .true.   !inner corner  '../Grids/LL-shape.uns.100.grid'
    !if(ip == 50) iprint = .true.   !inner corner  '../Grids/LL-shape.uns.100.grid'
    !if(ip == 73) iprint = .true.   !inner corner  '../Grids/LL-shape.uns.100.grid'


    ! setting of local dof
    deg = maxval(grid%elem(supp(1:N,1))%deg) ! setting of p of reconstruction
    !deg1 = deg  ! we keep the degree
    deg1 = deg + 1
    ! deg1 = deg + 2  ! deg +2 is too many, we have more unknowns than input DOF data
    
    dof = DOFtriang( deg1 )  

    L_rule => state%space%L_rule(deg1)


    ! setting of quadrature
    Qnum = state%space%Qdeg(deg1, 1)
    V_rule => state%space%V_rule(Qnum)

    Qdof = V_rule%Qdof
    allocate( Fx(1:Qdof, 1:nbDim) )
 

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

    if(iprint) write(*,*) grid%x(ip, :), inner

    do i=1, N
       ie  = supp(i,2)   ! inner index of the vertex

       elem => grid%elem(supp(i,1))
       
       
       !write(*,'(a8,5i5,2es12.4)' ) 'elemW:',ip, i, ie, elem%i,dN, elem%xc
       !write(*,'(2(a8,i5))' ) 'deg = ',deg1,',   dof = ', dof
       !write(21,*) elem%xc

       ! setting of the barycentric coordinates of the Lagrange nodes in the correct order
       call  Set_FE_Nodal_functionsREDO( deg1, dof, ie,  lambda(1:dof, 1:3) )
       call  Eval_L_Direct(L_rule, dof, lambda(1:dof, 1:2), phi_loc(1:dof, 1:dof))

       ! if(iprint) then
       !    do k=1,dof
       !       write(*,'(a6, 2i3, 3es12.4,a1,300es12.4)') &
       !            'lam:2w',supp(i,1), k, lambda(k, :),'|', phi_loc(1:dof, k)
       !    enddo
       ! endif


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
       call  Eval_Dw_Elem_plus(elem, V_rule, Dwi(0:2, 1:ndim, 1:Qdof) )

       ! TAKES VALUES from HO RECONSTRUCTION
       !Dwi(0, 1:ndim, 1:Qdof) = matmul(elem%wSS(1:ndim, 1:dof, 0), V_rule%phi(1:dof, 1:Qdof) )
       !Dwi(1, 1:ndim, 1:Qdof) = matmul(elem%wSS(1:ndim, 1:dof, 0), Dphi(1:dof, 1, 1:Qdof) )
       !Dwi(2, 1:ndim, 1:Qdof) = matmul(elem%wSS(1:ndim, 1:dof, 0), Dphi(1:dof, 2, 1:Qdof) )

       ! variant polynomially preserving
       !Dwi(0, 1:ndim, 1:Qdof) = matmul(elem%wSS(1:ndim, 1:dof, 0), V_rule%phi(1:dof, 1:Qdof) )
       !Dwi(1, 1:ndim, 1:Qdof) = matmul(elem%wSS(1:ndim, 1:dof, 1), V_rule%phi(1:dof, 1:Qdof) )
       !Dwi(2, 1:ndim, 1:Qdof) = matmul(elem%wSS(1:ndim, 1:dof, 2), V_rule%phi(1:dof, 1:Qdof) )

       call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof,1:nbDim),Fx(1:Qdof, 1:nbDim) )

       ! the hat function
       ie1 = mod(  ie, 3) + 1
       ie2 = mod( ie1, 3) + 1

       ! the hat function and its (constant) gradient on elem
       Dhat(0, 1:Qdof) = V_rule%lambda(1:Qdof, ie2) 
       Dhat(1, :) = (grid%x(elem%face(idx, ie1), 2) - grid%x(elem%face(idx, ie2), 2) ) /elem%area/2
       Dhat(2, :) = (grid%x(elem%face(idx, ie2), 1) - grid%x(elem%face(idx, ie1), 1) ) /elem%area/2
       
       
       ! the functions used for the computation
       !if(iprint) then
       !   do j=1,Qdof
       !      write(1000, *) Fx(j, 1:2), Dwi(0:2, 1, j) 
       !      write(1001, *) Fx(j, 1:2), Dhat(0:2,  j) 
       !      write(1002, *) Fx(j, 1:2), Dhat(0:2,  j) * Dwi(0:2, 1, j) 
       !   enddo
       !endif


       ! the setting of the stiff matrix
       do l1=1,dof
          do l2=1, l1  !dof  !only symmetric part
             k = itrans(i, l1)
             l = itrans(i, l2)

             if( k > 0 .and. l > 0 ) then
                val = dot_product(weights(1:Qdof), &
                     DL_phi(l1, 1, 1:Qdof) * DL_phi(l2, 1, 1:Qdof) + &
                     DL_phi(l1, 2, 1:Qdof) * DL_phi(l2, 2, 1:Qdof)  )

                ! direct setting of the global matrix
                A( k, l) = A( k, l) + val
                if( l /= k) A( l, k) = A( l, k) + val


                !A_loc(l1, l2) = dot_product(weights(1:Qdof), &
                !     DL_phi(l1, 1, 1:Qdof) * DL_phi(l2, 1, 1:Qdof) +&
                !     DL_phi(l1, 2, 1:Qdof) * DL_phi(l2, 2, 1:Qdof)  )
                !A_loc(l2, l1) = A_loc(l1, l2)
                !if(iprint) write(*,'(a6,8i5)') 'idx:',ip, l1, l2, k, l 
             endif

          enddo ! l2 = 1,dof
          
          ! RHS
          if(k > 0) then
             ! RHS = D( hat * u_h) =  D hat * u_h + hat * D u_h
             do j = 1, ndim
                Dhat(3, 1:Qdof) = Dhat(1,1:Qdof) * Dwi(0,j,1:Qdof) +  Dhat(0,1:Qdof) * Dwi(1,j,1:Qdof)
                Dhat(4, 1:Qdof) = Dhat(2,1:Qdof) * Dwi(0,j,1:Qdof) +  Dhat(0,1:Qdof) * Dwi(2,j,1:Qdof)

                ! NO HAT function
                !Dhat(3, 1:Qdof) = Dwi(1,j,1:Qdof)
                !Dhat(4, 1:Qdof) = Dwi(2,j,1:Qdof)

                val =  dot_product(weights(1:Qdof), &
                     DL_phi(l1, 1, 1:Qdof) * Dhat(3, 1:Qdof) +&
                     DL_phi(l1, 2, 1:Qdof) * Dhat(4, 1:Qdof)  )

                A( k, dN+j ) = A( k, dN+j) + val
                
                !A_loc(l1, dof+j) = val
             enddo
          endif

       enddo  ! l1 = 1,dof
       
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

    !print*,'dN = ', dN, grid%x(ip, :)
    !do k=1,dN
    !   write(*,'(a8,i5, 300es12.4)') 'Stiff:',k, A(k, :)
    !enddo
    
    ! solution of the local problem on the given patch
    if(dN > 0) &
         call SolveLocalMatrixProblem(dN, A(1:dN, 1:dN), ndim, A(1:dN, dN+1 : dN+ndim) ) 

    !print*,'------------------------------------------'
    !do k=1,dN
    !   write(*,'(a8,i5, 300es12.4)') 'Stiff:',k, A(k, :)
    !enddo

    ! setting of the solution to particular elements
    do i=1, N
       elem => grid%elem(supp(i,1))

       ww(:,:) = 0.
       if(dN > 0)  then
          do l1=1,dof
             k = itrans(i, l1)
             
             if(k > 0) then  ! setting the DG functions from the FE ones
                do j=1, ndim
                   ww(j, 1:dof) = ww(j, 1:dof) + phi(i, l1, 1:dof) * A(k, dN + j)
                enddo
             endif
          enddo
       endif

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

       ! plotting of one "hat function"
       if(iprint ) then
          call PlotElemFunction3D(900, elem, dof, ww(1, 1:dof) )
       endif

       elem%HO_deg = max(elem%HO_deg, deg1) 

       do j=1,ndim
          elem%wS(j,1:dof) = elem%wS(j, 1:dof) + ww(j, 1:dof) 
       enddo


       ! plotting of a function (including the previous one computed for lower ip
       !if( sqrt(dot_product( grid%x(ip, 1:2), grid%x(ip, 1:2)) ) < 1E-3) & 
       !if( ip == 8 ) &
       !call PlotElemFunction3D(100, elem, dof, elem%wS(1, 1:dof) )


    enddo

    deallocate(Fx)
    deallocate( phi, phi_loc, Dphi, DL_phi, itrans, Dwi, w_BC, Dwi_BC, Dhat, A)
    
  end subroutine LocalVertexProblemRecovery

  !>  setting of barycentric coordinates of the FE nodal basis functions in the correct order
  subroutine Set_FE_Nodal_functionsREDO( deg, dof, ie, lambda )
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

    k = 3

    ! edge FE basis fucntions
    if(deg > 1) then

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


    if(k /= dof) stop 'Inconsistency in Set_FE_Nodal_functions in the module higher_order_local'

  end subroutine Set_FE_Nodal_functionsREDO


  !> create the pair of FE DOF on elem (local index) and the patch \f$\omega_a \f$ (global index)
  !> itrans(i,k) : i = element index, k = local FE index, itrans(i,k) = global FE index
  !> moreover, it set dN - size of the patch problem
  subroutine CreateLocal_GlobalPairsREDO(deg, dof, N, itrans, inner, dN)
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
    

  end subroutine CreateLocal_GlobalPairsREDO

  !> assembling of the global matrix by the local blocks
  subroutine  AssembGlobalMatrixREDO(dN, dof, A, A_loc, itrans )
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


  end subroutine AssembGlobalMatrixREDO

  !> settingf of the BC for the local vertex problem, array w_BC contains the FE basis 
  !> coefficients
  subroutine Set_BC_LocalvertexProblemREDO(elem, i, ie, N, deg1, dof, w_BC )
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

    else ! not first and not last element, only first vertex function
       k = 1 ! the ordering of the vertexes
       xi(1:2) = grid%x(elem%face(idx, ie), 1:2)
       call Exact_Scalar(xi(1:nbDim), w_BC(1:ndim, k), state%time%ctime )

    endif


    deallocate(xi)
  end subroutine Set_BC_LocalvertexProblemREDO


  !> evaluate the error estimate based on the solution of the HO local problem
  subroutine Set_EE_HO_LocalProblemREDO(elem)
    class(element), target, intent(inout) :: elem
    type(volume_rule), pointer :: V_rule
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:,:), allocatable :: Dphi
    real, dimension(:,:,:), allocatable :: wi  
    real, dimension(:,:,:), allocatable :: DwExact  ! exact solution and its deriv in integ nds
    real, dimension(:), allocatable :: weights
  
    integer :: degP, dofP, dof, Qnum, Qdof, dofM
    integer :: iff, k

    dof = elem%dof   ! dof of the approximate solution 

    degP = elem%HO_deg
    dofP = DOFtriang( degP )  ! dof of the HO reconstruction

    dofM = max(dof, dofP)

    Qnum = state%space%Qdeg(degP, 1)
    V_rule => state%space%V_rule(Qnum)

    Qdof = V_rule%Qdof

   ! quadrature weights
    allocate(weights(1:Qdof) )
    call Eval_V_Weights_plus(elem, V_rule, weights(1:Qdof))

   
    phi => V_rule%phi(1: dofM, 1:Qdof)

    allocate( Dphi(1:dofM, 1:nbDim, 1:Qdof) )
    call Eval_Dphi_plus(elem, V_rule, dofM,  Dphi(1:dofM, 1:nbDim, 1:Qdof) )


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
    !call PlotElemFunction3D(200, elem, dofP, elem%wS(1, 1:dofP) )

    

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

  end subroutine Set_EE_HO_LocalProblemREDO


  !> computation of the polynomial projections on the vertex patch
  subroutine  LocalPatchProjections( ip, inner, N, supp )
    integer, intent(in) :: ip      ! index of node
    logical, intent(in) :: inner   ! inner vertex?
    integer, intent(in) :: N   ! number of elememnts sharing vertex
    integer, dimension(1:N,1:2), intent(in) :: supp   !idx of elems sharing vertex
    real, dimension(:,:), allocatable :: MM ! solution matrix
    real, dimension(:,:), allocatable :: Mass ! mass matrix
    real, dimension(:,:), allocatable :: Stiff ! mass matrix
    real, dimension(:,:), allocatable :: rhs ! right-hand-side
    real, dimension(:,:), allocatable :: rhsSff ! right-hand-side
    real, dimension(:,:), allocatable :: xi ! barycentric coordinates
    real, dimension(:,:), allocatable :: Fx ! real physical coordinates
    real, dimension(:,:,:), allocatable :: phi ! the test functions
    real, dimension(:,:), allocatable :: phiW ! the test functions multiplied by weights
    real, dimension(:,:,:), allocatable :: Dphi ! Der of test functions
    real, dimension(:,:,:), allocatable :: DphiW ! Der of test functions multiplied by weights
    real, dimension(:,:), pointer :: phi0 ! the test functions
    real, dimension(:,:,:), allocatable :: Dphi0 ! the test functions
    real, dimension(:), allocatable :: weights, wi
    real, dimension(:,:), allocatable :: Dwi
    real, dimension(:,:), allocatable :: wii, qi  ! local store arrays
    real, dimension(:,:), allocatable :: estim  !difference between the HO rec and projection
    real, dimension(:), allocatable :: val
    class(element), pointer :: elem, elem1
    type(volume_rule), pointer :: V_rule

    integer :: i, ie, ie1, ie2, degP, deg1, dofP,  dof1, Qnum, Qnum1, Qdof, Qdof1, dN, dofE
    integer :: j, k, l, l1, l2, ndimL
    real :: L2weight, est_min

    ndimL = ndim

    ! setting of local dof
    !degP = minval(grid%elem(supp(1:N,1))%HO_deg) ! degree of the  available HO reconstruction
    degP = maxval(grid%elem(supp(1:N,1))%HO_deg) ! degree of the  available HO reconstruction
    !degP = maxval(grid%elem(supp(1:N,1))%deg) + 1
    dofP = DOFtriang(degP)    

    print*
    write(*,'(a8,3i5, 2es12.4, 20i5)') '###AA', &
         ip, degP, minval(grid%elem(supp(1:N,1))%HO_deg),  grid%x(ip,:),&
         grid%elem(supp(1:N,1))%HO_deg
    write(*,'(a8,3i5, 2es12.4, 20i5)') '###AA', &
         ip, degP, minval(grid%elem(supp(1:N,1))%HO_deg),  grid%x(ip,:),&
         grid%elem(supp(1:N,1))%deg

    ! increase of quadrature
    if(degP == MaxDegreeImplemented ) then
       print*,'$$ alert, possible troubles in ho-recovery.f90'
       Qnum = elem%Qnum
    else
       Qnum = state%space%Qdeg(degP+1, 1)
    endif

    V_rule => state%space%V_rule(Qnum)
    Qdof = V_rule%Qdof


    ! Mass(i,j) = \int_{K and its neighbours} \phi_i \phi_j dx
    ! rhs(i) = \int_{K and its neighbours} \phi_i w_h dx
    allocate(Mass(1:dofP, 1:dofP), Stiff(1:dofP, 1:dofP))
    allocate(rhs(1:dofP, 1:ndimL), rhsSff(1:dofP, 1:ndimL) )
 
    Mass(:,:) = 0.
    Stiff(:,:) = 0.
    rhs(:,:) = 0.
    rhsSff(:,:) = 0.
    

    allocate(phi(1:N, 1:dofP, 1:Qdof) ) ! value of the test function in integ nodes


    ! initial element
    elem => grid%elem(supp(1,1))

    ! we go over triangles is the patch
    do i=1, N
       ie  = supp(i,1)   ! inner index of the vertex
       elem1 => grid%elem(ie)
       !dof1  = elem1%dof
       dof1  = DOFtriang(elem1%HO_deg)
       
       !write(*,'(a8,4i5, 4es12.4)') '#E#E', ip, ie, dof1, DOFtriang(elem1%HO_deg), &
       !     elem1%xc, grid%x(ip,:) 

       ! MUST be the same
       Qnum1 = Qnum  !! probably more accurate than the previous one
       Qdof1 = state%space%V_rule(Qnum1)%Qdof

       ! Fx integ nodes on elem1 - real physical cordinates
       ! xi barycentric coordinates of Fx (integ nodes on elem1) with respect elem
       allocate(Fx(1:Qdof1, 1:2), xi(1:Qdof1, 1:2) )
       call ComputeF(elem1, Qdof1, state%space%V_rule(Qnum1)%lambda(1:Qdof1, 1:2), &
            Fx(1:Qdof1, 1:2) )
       call BarycCoord(elem, Qdof1, Fx(1:Qdof1, 1:2),  xi(1:Qdof1, 1:2) )

       !allocate(phi(1:dofP, 1:Qdof1) ) ! value of the test function in integ nodes
       allocate(Dphi(1:dofP, 1:nbDim, 1:Qdof1) ) ! value of Deriv test functions


       call Eval_phi_Qnode(elem, dofP, Qdof1, xi(1:Qdof1, 1:nbDim), &
            phi(i, 1:dofP, 1:Qdof1), Dphi(1:dofP, 1:nbDim, 1:Qdof1) )


       allocate(wi(1:Qdof1),  weights(1:Qdof1) )
       allocate(Dwi(1:Qdof1, 1:nbDim)  )
       call Eval_V_Weights_plus(elem1, state%space%V_rule(Qnum1), weights(1:Qdof1))

       allocate(phiW(1:dofP, 1:Qdof1) )           ! phi multiplied by weights
       allocate(DphiW(1:dofP, 1:nbDim, 1:Qdof1) ) ! Dphi multiplied by weights

       do l=1, dofP
          phiW(l, 1:Qdof1) = phi(i, l, 1:Qdof1) * weights(1:Qdof1)

          do l1=1,nbDim
             DphiW(l, l1, 1:Qdof1) = Dphi(l, l1, 1:Qdof1) * weights(1:Qdof1)
          enddo
       enddo

       ! adding of Mass and Stiff  matrices
       do l1 = 1, dofP
          do l2 = 1, dofP
             Mass(l1, l2) = Mass(l1, l2) &
                  + dot_product(phiW(l1, 1:Qdof1), phi(i, l2, 1:Qdof1) )
             
             do l=1,nbDim
                Stiff(l1, l2) = Stiff(l1, l2) &
                     + dot_product(DphiW(l1, l, 1:Qdof1), Dphi(l2, l, 1:Qdof1) )
             enddo
          enddo
       enddo
       
       ! adding of rhs
       phi0 => state%space%V_rule(Qnum1)%phi(1:dof1, 1:Qdof1)

       allocate( Dphi0(1:dof1, 1:nbDim, 1:Qdof1) )
       call Eval_Dphi_plus(elem1, state%space%V_rule(Qnum1), dof1, Dphi0(1:dof1, 1:nbDim, 1:Qdof1) )


       do k=1, ndimL
          ! wi = values of w at integ nodes of elem1
          !!WS!! wi(1:Qdof1) = matmul(elem1%w(0, (k-1)*dof1 +1 : k*dof1), &
          wi(1:Qdof1) = matmul(elem1%wS(k, 1 : dof1),  phi0(1:dof1, 1:Qdof1) )
          !write(*,'(a8, 2i5, 20es12.4)') 'de3:',elem%i,k, wi(1), (elem1%xc(1)/40)**2, elem1%xc
          ! product (wi, phi0)
          do l=1, dofP
             rhs(l, k) = rhs(l, k)  &
                  + dot_product(wi(1:Qdof1), phiW(l, 1:Qdof1) )
          enddo
          

          ! Dwi in integ nodes
          do l=1,nbDim
             !!WS!! Dwi(1:Qdof1, n) = matmul(elem1%w(0, (k-1)*dof1 +1 : k*dof1), &
             Dwi(1:Qdof1, l) = matmul(elem1%wS(k, 1 : dof1), Dphi0(1:dof1, l, 1:Qdof1))
          enddo
          
          ! product (Dwi, Dphi0)
          do l=1, dofP
             do l1=1,nbDim
                rhsSff(l, k) = rhsSff(l, k)  &
                     + dot_product(Dwi(1:Qdof1, l1), DphiW(l, l1, 1:Qdof1) )
             enddo
             
          enddo
          
       enddo

       deallocate(Fx, xi, phiW, wi, weights)
       deallocate(Dphi, DphiW, Dwi, Dphi0)
       

    enddo ! i=1,N

    ! the H1-norm (comment for the L^2-norm)
    L2weight = 0.0
    Mass(1:dofP, 1:dofP) = Mass(1:dofP, 1:dofP) +  L2weight * Stiff(1:dofP, 1:dofP)
    rhs(1:dofP, 1:ndimL) = rhs(1:dofP, 1:ndimL) +  L2weight * rhsSff(1:dofP, 1:ndimL)

    ! array for the difference between the HO recovvery and the projection to the patch
    allocate(estim(-2:1, 0:5), val(1:ndimL) )
    estim = 0.

    !array for storing of the projections
    do i=1, N
       ie  = supp(i,1)   ! index of the triangle
       elem1 => grid%elem(ie)

       allocate(elem1%wSS(-2:1, 1:dofP, 1:ndimL) )
       elem1%wSS = 0.
    end do

    ! levels of the projection
    do l = 1, -2, -1
       deg1 = degP + l - 1  !  degree of the projection
       if(deg1 >= 0) then
          dof1 = DOFtriang(deg1)
          
          allocate(MM(1:dof1, 1:dof1+ndimL ) )
          MM(1:dof1, 1:dof1) = Mass(1:dof1, 1:dof1)
          MM(1:dof1, dof1+1:dof1+ndimL) = rhs(1:dof1, 1:ndimL)
          
          call SolveLocalMatrixProblem(dof1, MM(1:dof1, 1:dof1), ndimL, &
               MM(1:dof1, dof1+1:dof1+ndimL) )
          
          ! setting of the projection to the particular elements from the patch
          allocate(qi( 1:dofP, 1:ndimL))
          allocate(wii( 1:Qdof, 1:ndimL))
          
          do i=1, N           ! we go over triangles is the patch
             ie  = supp(i,1)  ! index of the triangle
             elem1 => grid%elem(ie)

             Qnum1 = elem1%Qnum
             Qdof1 = elem1%Qdof
             phi0 =>  state%space%V_rule(Qnum1)%phi(1:max(dof1, dofP), 1:Qdof1)
             

             ! SIMPLIFY for linear elements Mass matrix is only one number !!!!
             ! local Mass matrix
             if(dof1 > size(elem1%Mass%Mb, 1) ) then
                ! necasary to compute
                allocate(weights(1:Qdof1) ) 
                call Eval_V_Weights(elem1, weights(1:Qdof1) )
                !print*,'#DE#EO#',dof1,  size(elem1%Mass%Mb, 1), elem1%dof
                do l1=1,dof1
                   MM(l1, 1:dof1) = matmul(phi0(1:dof1, 1:Qdof1), &
                        weights(1:Qdof1)*phi0(l1, 1:Qdof1) )
                enddo
                deallocate(weights)
                !do k=1,dof1
                !   write(*,'(a2,i5,20es12.4)') 'M:',k,MM(k, 1:dof1)
                !enddo

             else
                ! use of the precomputed one
                MM(1:dof1, 1:dof1) = elem1%Mass%Mb(1:dof1, 1:dof1) ! local Mass matrix
             endif

             do k=1,ndimL
                !print*,'_____________________________',Qdof,elem1%Qdof, Qnum,dof1
                
                wii(1:Qdof,k) = matmul( MM(1:dof1, dof1+k),  phi(i, 1:dof1, 1:Qdof)) 
                call IntegrateVectorBplus(elem1, Qnum, dofP, wii(1:Qdof, k), qi(1:dof1, k) )
                
                
                !write(*,'(a6,20es12.4)') 'ED rhs:',MM(1:dof1,dof1+ k)
                !write(*,'(a6,200es12.4)') 'ED wii:',wii(1:Qdof, k)
                !write(*,'(a6,20es12.4)') 'ED qi:', qi(1:dof1, k)
             enddo
             
             !do k=1,dof1
             !   write(*,'(a2,i5,20es12.4)') 'M:',k,Mass(k, 1:dof1)
             !enddo
             
             call SolveLocalMatrixProblem(dof1, MM(1:dof1, 1:dof1), ndimL, qi(1:dof1, 1:ndimL))
             
             !write(*,'(a6,20es12.4)') 'ED qi:', qi(1:dof1, 1)
             
             elem1%wSS(l, 1:dof1,1:ndimL) = qi(1:dof1, 1:ndimL)
             
             !call PlotElemFunction3D(30+l, elem1, dofP,  elem1%wSS(l, 1:dofP,1) )


             ! difference between HO_recovery (elem%wS) and the projection (elem%wSS)


             do k=1,ndimL
                wii(1:Qdof1,k) = matmul(elem1%wSS(l, 1:dof1, k),  phi0(1:dof1, 1:Qdof1) ) &
                     - matmul( elem1%wS(k, 1:dofP), phi0(1:dofP, 1:Qdof1))
             enddo

             call IntegrateSquareVectorFunction2(elem1, wii(1:Qdof1, 1:ndimL), val(1:ndimL) )
             estim(l,1) = estim(l,1) + sum(val(1:ndimL))
             estim(l,0) = estim(l,0) + elem1%area
          enddo ! i=1,N

          estim(l,1) = sqrt(estim(l,1))  ! error
          estim(l,2) = estim(l,1) / estim(l,0)**( (deg1+1.)/2.)  ! constant c_p
          estim(l,3) = state%space%adapt%tol_max * sqrt(elem%area/state%space%domain_volume)
          estim(l,4) = (estim(l,3) / estim(l,2) )**(2./(deg1+1) )  ! "new" area
          estim(l,5) = (deg1+1)*(deg1+2)/2 / estim(l,4)  ! density of DOF

          write(*,'(a8,2i5, 9es12.4)') 'estim:',l,deg1,estim(l,0:)
          deallocate(wii, qi)
          deallocate(MM)
       end if ! deg1 > 0
    end do ! l=0,1

    deallocate(Mass, Stiff, rhs, rhsSff, phi)
    do i=1, N
       ie  = supp(i,1)   ! index of the triangle
       elem1 => grid%elem(ie)
       deallocate(elem1%wSS)
    enddo

    est_min = 1E+50
    do l=-2,1
       if(estim(l,0) > 0.) then
          if(estim(l,5) <  est_min) then
             est_min = estim(l,5)
             k = l
          endif
       endif
    enddo

    write(80+state%space%adapt%adapt_level,*) grid%x(ip,:), k,  degP + k - 1 



  end subroutine LocalPatchProjections

end module higher_order_recovery
