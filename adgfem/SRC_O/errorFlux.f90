!> compute the computational error in flux quantities
module errorFlux
  !use mesh_oper
  use main_data
  !use eval_sol
  use io_sub
  use model_oper
  !use set_solution
  use problem_oper
  use euler_problem
  use eval_jumps


  implicit none

  public:: SetRHSDualSTerror
  public:: EvalDualSubElem
  public:: NodesInTwoTriangles
  public:: NodeInTriangle
  public:: ComputeFluxError
  public:: ComputeFluxErrorElem
contains

  !> prepare the RHS for the computation of the RHS for the evaluation of
  !> the space-time dual norm for STDGM in Fenics
  subroutine SetRHSDualSTerror( )
    class(element), pointer :: elem1, elem2
    integer :: i, j, k, i1, j1, l, ie, idat1, idat2
    integer :: Nx, Nt, N, Gnum, Gdof
    integer:: Mx, My, Mt
    real :: t, t1, t2, val, val2, val_N
    real, dimension (1:2) :: xi, xp, xmin, xmax

    if(nbDim /= 2 .or. ndim /= 1) then
       print*,'Case nbDim ==',nbDim,' & ndim ==',ndim, &
            ',   SetRHSDualSTerror not implemented in errorFlux.f90'
       stop
    endif

    Nx = 4  ! number of subelements for each bi-triangle - Nx * Nx squares
    Nt = 4  ! number of time sublayers for each time interval
    Nx = int(state%space%adapt%tol_min)  ! GIVEN IN *.INI FILE !!!!!!!!
    Nt = int(state%space%adapt%tol_min)


    N = (grid%nelem/2)**0.5  ! number of elements in one direction

    Mx = 0
    My = 0
    Mt = (state%time%iter_loc - 1) * Nt

    idat1 = 11
    idat2 = 12
    if(state%time%iter == 1) then
       open(idat1, file='DualRHS1', status='Unknown')
       write(idat1, *) N*Nx, N*Nx, int(state%time%FinTime/ state%time%tau(1) + 0.5)*Nt, &
            N*grid%elem(1)%diam/2**0.5, state%time%FinTime, &
            grid%elem(1)%diam/2**0.5, state%time%tau(1), grid%elem(1)%CTn, Nx, Nt

       !open(idat2, file='DualRHS2', status='Unknown')
       !write(idat2, *) N*Nx, N*Nx, int(state%time%FinTime/ state%time%tau(1) + 0.5)*Nt, &
       !     N*grid%elem(1)%diam/2**0.5, state%time%FinTime, &
       !     grid%elem(1)%diam/2**0.5, state%time%tau(1), grid%elem(1)%CTn, Nx, Nt
    else
       open(idat1, file='DualRHS1', status='OLD', position='append')
       !open(idat2, file='DualRHS2', status='OLD', position='append')
    endif



    do k=1, Nt   ! time levels
       Mt = Mt + 1                     ! next level in t-direction
       Mx = 0

       do i=1, N
          do i1 = 1, Nx                   ! subsquares of element
             Mx = Mx + 1               ! next level in x-direction

             My = 0
             do j = 1, N
                ie = 2*((i-1)*N + j) - 1
                elem1 => grid%elem(ie)          ! elem1 and elem2 form a square
                elem2 => grid%elem(ie+1)


                do l=1,2
                   xmin(l) = min(minval(grid%x(elem1%face(idx,:), l)), &
                        minval(grid%x(elem2%face(idx,:), l)) )

                   xmax(l) = max(maxval(grid%x(elem1%face(idx,:), l)), &
                        maxval(grid%x(elem2%face(idx,:), l)) )
                enddo

                do j1 = 1, Nx                ! subsquares of element
                   My = My + 1               ! next level in y-direction

                   call EvalDualSubElem(N, Nx, Nt, k, i, j, i1, j1,  &
                        xmin, xmax, elem1, elem2, Set_R_s_scalar, Set_f_s_scalar,val, val2, val_N)

                   !write(500+Mt, '(3i5, 5es12.4)') Mx, My, Mt, val, &
                   !     val *(xmax(1)-xmin(1))*(xmax(2)-xmin(2))/Nx/Nx*state%time%tau(1)/Nt, &
                   !     val_N,  val_N*(xmax(1)-xmin(1))*(xmax(2)-xmin(2))/Nx/Nx


                   !val = 1.
                   !val_N = 0.

                   write(idat1,'(3i5, 5es12.4)') Mx, My, Mt, val, &
                        val *(xmax(1)-xmin(1))*(xmax(2)-xmin(2))/Nx/Nx*state%time%tau(1)/Nt, &
                        val_N,  val_N*(xmax(1)-xmin(1))*(xmax(2)-xmin(2))/Nx/Nx

                   !write(idat2,'(3i5, 5es12.4)') Mx, My, Mt, val2, &
                   !     val2 *(xmax(1)-xmin(1))*(xmax(2)-xmin(2))/Nx/Nx*state%time%tau(1)/Nt, &
                   !     val_N,  val_N*(xmax(1)-xmin(1))*(xmax(2)-xmin(2))/Nx/Nx

                enddo !end of j1
             enddo ! end of j
          enddo  ! end of i1
       enddo   ! end of i
    enddo  ! end k=1,Nt

    close(idat1)
    !close(idat2)

  end subroutine SetRHSDualSTerror

  ! evaluation of RHS for one subelement
  subroutine EvalDualSubElem(N, Nx, Nt, k, i, j, i1, j1, xmin, xmax, elem1, elem2, &
       Set_R_s, Set_f_s, val, val2, val_N)
    interface
      subroutine Set_R_s(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1   ! inverse of Reynolds number
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

    integer, intent(in) :: N, Nx, Nt, k, i, j, i1, j1
    type(element), intent(in) :: elem1, elem2
    class(element), pointer :: elem, elemR
    integer :: elemi, elemK, facei, facei1, faceiR, faceiR1, faceiR2, iface, ifaceR
    type(Gauss_rule), pointer :: G_rule, SG_rule
    type(volume_rule), pointer :: V_rule
    real, dimension(1:2) , intent(in):: xmin, xmax
    real, dimension (0:4, 1:2) :: xp
    real, dimension (0:4, 1:3) :: xpr
    real, dimension (:,:), allocatable :: xi, xii, xiR, xiRR, wi, wE, Dwt, wiR, DivF
    real, dimension (:,:,:), allocatable :: Dwx, DwxR, f_s, R_s, f_sR, R_sR
    real, intent(out) :: val, val2, val_N
    integer :: Gdof, Gnum, Qdof, Qnum, Sdof, Snum, gi,gj, gk, ie, l, l1, nedges, n1, n2
    real, dimension (:), allocatable :: f
    real, dimension (:,:), allocatable :: Re_1
    integer, dimension(:,:), allocatable :: iedges, itri
    real:: t1, t2, t, rlen
    real :: px1,px2,px3,px4
    real :: val_f, val_S, val2_S
    logical :: inside

    ! integrals with respect to time
    Gnum = 2   ! 15 is the maximal one !!!!!!!!!
    G_rule => state%space%G_rule(Gnum)
    Gdof = G_rule%Qdof

    ! space -- edge integrals
    Snum =  max(elem1%deg, elem2%deg)*2  ! 15 is the maximal one !!!!!!!!!
    Snum = 3
    SG_rule => state%space%G_rule(Snum)
    Sdof = SG_rule%Qdof

    ! volume integrals
    Qnum = max(elem1%deg, elem2%deg)*3    ! 15 is the maximal one !!!!!!!!!
    Qnum = elem1%Qnum
    V_rule => state%space%V_rule(Qnum)
    Qdof = V_rule%Qdof

    if(Sdof > Qdof) then
       print*,'Problems in memory allocation in xi, xiR in errorFlux.f90'
       stop
    endif

    allocate( xi(1:Qdof, 1:2), xii(1:Qdof, 1:3), xiR(1:Qdof, 1:3), xiRR(1:Qdof, 1:3) )
    allocate( wi(1:Qdof, 1:ndim), wiR(1:Qdof, 1:ndim), wE(1:Qdof, 1:ndim), Dwt(1:Qdof, 1:ndim) )

    allocate(Dwx(1:Qdof,1:ndim, 1:nbDim), DwxR(1:Qdof,1:ndim, 1:nbDim) )
    allocate(DivF(1:Qdof, 1:ndim) )
    allocate(f_s(1:Qdof, 1:nbDim, 1:ndim), R_s(1:Qdof, 1:nbDim, 1:ndim) )
    allocate(f_sR(1:Qdof, 1:nbDim, 1:ndim), R_sR(1:Qdof, 1:nbDim, 1:ndim) )

    val = 0.
    val2 = 0.
    val_f = 0.
    val_S = 0.
    val2_S = 0.

    allocate(Re_1(1:iRe, 1:Qdof) )
    Re_1(1:iRe, 1:Qdof) = state%model%Re1

    allocate( f(1:ndim) )

    t1 = state%time%ttime - state%time%tau(1)
    t2 = state%time%ttime


    ! we split square in 4 triangles by diagonal cuts
    xp(1, 1) = xmin(1) + 1. * (i1-1) / Nx * (xmax(1) - xmin(1) )
    xp(1, 2) = xmin(2) + 1. * (j1-1) / Nx * (xmax(2) - xmin(2) )

    xp(2, 1) = xmin(1) + 1. * (i1-0) / Nx * (xmax(1) - xmin(1) )
    xp(2, 2) = xmin(2) + 1. * (j1-1) / Nx * (xmax(2) - xmin(2) )

    xp(3, 1) = xmin(1) + 1. * (i1-0) / Nx * (xmax(1) - xmin(1) )
    xp(3, 2) = xmin(2) + 1. * (j1-0) / Nx * (xmax(2) - xmin(2) )

    xp(4, 1) = xmin(1) + 1. * (i1-1) / Nx * (xmax(1) - xmin(1) )
    xp(4, 2) = xmin(2) + 1. * (j1-0) / Nx * (xmax(2) - xmin(2) )

    xp(0, 1) = sum(xp(1:4, 1)) / 4
    xp(0, 2) = sum(xp(1:4, 2)) / 4

    nedges = 8
    allocate(iedges(1:nedges, 1:2), itri(1:3, 1:3) )

    iedges(1:nedges, 1) = (/ 1, 2, 3, 4, 0, 0, 0, 0/)  ! index of edges of four subtriangles
    iedges(1:nedges, 2) = (/ 2, 3, 4, 1, 1, 2, 3, 4/)  ! of one small square

    itri(1, 1:3) = (/3, 1, 2 /)                   ! if lambda_i(node) = 0  then node lies on
    itri(2, 1:3) = (/3, 1, 2 /)                   !  itri(elemi,i)-th face of triangles elemi
    itri(3, 1:3) = (/2, 3, 1 /)                   !   inverse mapping to ^^^^^

    ! check of the indexes of edges
    do elemK=1,2
       if( elemK == 1) elem => grid%elem(elem1%i)
       if( elemK == 2) elem => grid%elem(elem2%i)

       do l = 1, elem%flen
          l1 = mod(l, elem%flen) + 1

          !write(*,'(a6,4i5,3es12.4)') '#!!#',elemK,l,l1,elem1%flen

          xi(1, 1:2) = (grid%x(elem%face(idx,l),1:2) + grid%x(elem%face(idx,l1),1:2) )/ 2

          call NodeInTriangle(1, grid%x(elem%face(idx,1),1:2), grid%x(elem%face(idx,2),1:2), &
               grid%x(elem%face(idx,3),1:2), xi(1:1, 1:2), xiR(1:1, 1:3), inside)

          do gi = 1,3
             if(abs( xiR(1,gi)) < 1E-5 ) then
                if(itri(elemK, gi) /= l) then
                   write(*,*) 'Troubles in seeking edges in errorFlux.f90'
                   write(*, '(2es12.4,4i5,6es12.4)') &
                        xi(1, 1:2), elemK, gi,l,itri(elemK, gi), xiR(1, 1:3)
                   stop
                endif
             endif
          enddo

       enddo
    enddo


    ! we go over subedges of subtriangles
    do ie = 1,nedges   ! all edges
       !!!do ie = 5,nedges   ! only diagonal edges within the finite volume
       rlen = 1. / Nx
       if(ie > 4) rlen = rlen / 2   ! diagonal subfaces have a half size

       l  = iedges(ie, 1)
       l1 = iedges(ie, 2)

       do gi=1,Sdof  ! integ nodes, real coordinates      ??? orientation ??
          xi(gi, 1:2) = SG_rule%lambda(gi) * xp(l,1:2) &
               + (1. - SG_rule%lambda(gi)) * xp(l1,1:2)

       end do

       call NodesInTriangleEdges(elem1, elem2, Sdof, xi(1:Sdof, 1:2), xiR(1:Sdof, 1:3), &
            inside, elemi, facei)

       if(inside)  then  !! sub-face with a jump of the flux ??
          facei1 = mod(facei, 3) + 1  ! facei-th barycentric coordinate = 0,
                                      ! face1-th is the relative position of the edge

          iface = itri(elemi,facei)
          if( elemi == 1) elem => grid%elem(elem1%i)
          if( elemi == 2) elem => grid%elem(elem2%i)

          !do gi=1,Sdof  ! integ nodes, real coordinates      ??? orientation ??
          !   write(24,'(6es12.4,l3,6i5)') &
          !        xi(gi, 1:2), xiR(gi, 1:3), xiR(gi, facei1), &
          !        inside, l,l1,ie, elemi, facei, iface
          !enddo
          !write(24, *)

          !print*,'@@@',inside, elem%i, iface, elem%face(neigh, iface), elem%face(nei_i, iface)

          if(elem%face(neigh, iface) > 0) then ! interior edge
             l = elem%face(neigh, iface)
             elemR => grid%elem(l)
             ifaceR = elem%face(nei_i, iface)

             !write(26, *) (elem%xc(:) + elemR%xc(:) ) /2


             faceiR = itri(3,ifaceR)
             faceiR1 = mod(faceiR, 3) + 1
             faceiR2 = mod(faceiR1, 3) + 1

             xiRR(1:Sdof, faceiR) = 0.
             xiRR(1:Sdof, faceiR1) = 1.- xiR(1:Sdof, facei1)
             xiRR(1:Sdof, faceiR2) = xiR(1:Sdof, facei1)



             !write(*,'(a6,15i5)' ) &
             !        '##@@',elem%i,iface, facei, facei1, ifaceR, faceiR,faceiR1, faceiR2
             !
             !do n1 = 1,Sdof
             !   write(*,'(a6,1i5,6es12.4)' ) &
             !        '##@@', n1, xiR(n1,1:3), xiRR(n1,1:3)
             !enddo

             !print*,'---------------------------'
             !do n1=1,Sdof
             !   px1 = xiR(n1, 1)* grid%x(elem%face(idx, 2), 1) &
             !        + xiR(n1, 2)* grid%x(elem%face(idx, 3), 1) &
             !        + xiR(n1, 3)* grid%x(elem%face(idx, 1), 1)
             !
             !   px2 = xiR(n1, 1)* grid%x(elem%face(idx, 2), 2) &
             !        + xiR(n1, 2)* grid%x(elem%face(idx, 3), 2) &
             !        + xiR(n1, 3)* grid%x(elem%face(idx, 1), 2)
             !
             !   px3 = xiRR(n1, 1)* grid%x(elemR%face(idx, 2), 1) &
             !        + xiRR(n1, 2)* grid%x(elemR%face(idx, 3), 1) &
             !        + xiRR(n1, 3)* grid%x(elemR%face(idx, 1), 1)
             !
             !   px4 = xiRR(n1, 1)* grid%x(elemR%face(idx, 2), 2) &
             !        + xiRR(n1, 2)* grid%x(elemR%face(idx, 3), 2) &
             !        + xiRR(n1, 3)* grid%x(elemR%face(idx, 1), 2)
             !
             !   write(*,'(10es14.6)') px1,px2,px3,px4, px1-px3,px2-px4
             !   write(58,'(10es14.6)') px1,px2,px3,px4, px1-px3,px2-px4
             !
             !enddo


             !write(52, *) xmin(1:2), xmax(1:2)
             !stop

             ! time integration
             do gk = 1, Gdof
                t = ( (k-1) + G_rule%lambda(gk))/ Nt
                state%time%ctime = t1 + t * (t2 - t1)

                !t = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                call Eval_w_Qnode(elem, Sdof, xiR(1:Sdof, 1:nbDim), t, wi(1:Sdof, 1:ndim), &
                     Dwt(1:Sdof, 1:ndim), Dwx(1:Sdof,1:ndim, 1:nbDim) )

                call Eval_w_Qnode(elemR, Sdof, xiRR(1:Sdof, 1:nbDim), t, wiR(1:Sdof, 1:ndim), &
                     Dwt(1:Sdof, 1:ndim), DwxR(1:Sdof,1:ndim, 1:nbDim) )

                !write(*,'(a6,2i5,20es12.4)') 'xiR:',elemi, Sdof, xiR(1:Sdof, 1)
                !write(*,'(a6,2i5,20es12.4)') 'xiR:',elemi, Sdof, xiR(1:Sdof, 2)
                !do l=1,Sdof
                !   write(27,'(10es12.4)') xiR(l, 1:2), wi(l, 1), Dwx(l, 1, 1:2)
                !enddo



                ! diffusive terms
                call Set_R_s(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwx(1:Qdof, 1:ndim, 1:nbDim), Re_1,&
                     R_s(1:Qdof, 1:nbDim, 1:ndim),  elem%xi(0,1:Qdof, 1:nbDim) )

                ! convective terms
                call Set_f_s(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), f_s(1:Qdof, 1:nbDim, 1:ndim), &
                     elem%xi(0,1:Qdof, 1:nbDim), elem%i )

                ! total flux
                f_s(1:Qdof, 1:nbDim, 1:ndim)  = R_s(1:Qdof, 1:nbDim, 1:ndim) &
                     - f_s(1:Qdof, 1:nbDim, 1:ndim)


                ! diffusive terms
                call Set_R_s(ndim, nbDim, iRe, Qdof, wiR(1:Qdof,1:ndim), DwxR(1:Qdof, 1:ndim, 1:nbDim),Re_1,&
                     R_sR(1:Qdof, 1:nbDim, 1:ndim),  elem%xi(0,1:Qdof, 1:nbDim) )

                ! convective terms
                call Set_f_s(ndim, nbDim, Qdof, wiR(1:Qdof,1:ndim), f_sR(1:Qdof, 1:nbDim, 1:ndim), &
                     elem%xi(0,1:Qdof, 1:nbDim), elem%i )

                ! total flux
                f_sR(1:Qdof, 1:nbDim, 1:ndim)  = R_sR(1:Qdof, 1:nbDim, 1:ndim) &
                     - f_sR(1:Qdof, 1:nbDim, 1:ndim)

                DivF(1:Qdof, 1:ndim)  = 0.
                do gi=1,nbDim
                   DivF(1:Qdof, 1:ndim)  = DivF(1:Qdof, 1:ndim)  + &
                        (f_s(1:Qdof, gi, 1:ndim) - f_sR(1:Qdof, gi, 1:ndim)) * elem%n(iface, gi)
                enddo

                ! 1/2 = average of the piecewise constant test function
                DivF(1:Qdof, 1:ndim)  = DivF(1:Qdof, 1:ndim) / 2

                do gi = 1, Sdof
                   val_S = val_S + SG_rule%weights(gi) * G_rule%weights(gk) *DivF(gi, 1) *rlen

                   if(ie >= 5) &  ! ONLY diagonal edges
                        val2_S = val2_S &
                        + SG_rule%weights(gi) * G_rule%weights(gk) *DivF(gi, 1) *rlen
                enddo

                !!write(31,'(10es12.4)') xi(1, 1:2), DivF(1:Sdof, 1)


             enddo  ! end of gk

          endif  ! interior edge
       endif ! sub-face with a jump of the flux
    enddo ! end of ie = 1,nedges

    !write(24, *) '################################'
    val_S = - val_S / ( (elem1%area + elem2%area)/ (Nx * Nx) )
    val2_S = - val2_S / ( (elem1%area + elem2%area)/ (Nx * Nx) )


!    call NodeInTriangle(5, grid%x(elem1%face(idx,1),:), grid%x(elem1%face(idx,2),:), &
!         grid%x(elem1%face(idx,3),:), xp(0:4, 1:2), xpR(0:4, 1:3), inside)

!    elemK = 1
!    if(.not. inside) then
!       ! barycentric coordinates with respect to elem2
!       call NodeInTriangle(5, grid%x(elem2%face(idx,1),:), grid%x(elem2%face(idx,2),:), &
!            grid%x(elem2%face(idx,3),:), xp(0:4, 1:2), xii(1:5, 1:3), inside)
!
!       if(.not. inside) then
!          print*,'sqaure is not in elem1 and elem2'
!          stop
!       endif
!       xpR(0:4, 1:3) = xii(1:5, 1:3)
!       elemK = 2
!    endif

!    write(21,*) '#####',elemK
!    do l=0,4
!       write(21,'(20es12.4)') xp(l, 1:2), xpR(l, 1:3)
!    enddo
!    do l=1,3
!       write(21,'(20es12.4)') grid%x(elem1%face(idx,l),:)
!    enddo
!    stop


    ! we go over subtriangles
    do l=1,4
       l1 = mod(l, 4) + 1

       do gi=1,Qdof  ! integ nodes, real coordinates      ??? orientation ??
          xi(gi, 1:2) = V_rule%lambda(gi,1) * xp(l,1:2) &
               + V_rule%lambda(gi,2) * xp(l1,1:2) &
               + V_rule%lambda(gi,3) * xp(0,1:2)

       end do

       call NodesInTwoTriangles(elem1, elem2, Qdof, xi, xiR, inside, elemi)


       ! time integration
       do gk = 1, Gdof
          t = ( (k-1) + G_rule%lambda(gk))/ Nt
          state%time%ctime = t1 + t * (t2 - t1)

          !t = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if(elemi == 1) then
             call Eval_w_Qnode(elem1, Qdof, xiR(1:Qdof, 1:nbDim), t, wi(1:Qdof, 1:ndim), &
                  Dwt(1:Qdof, 1:ndim), Dwx(1:Qdof,1:ndim, 1:nbDim) )
          else
             call Eval_w_Qnode(elem2, Qdof, xiR(1:Qdof, 1:nbDim), t, wi(1:Qdof, 1:ndim), &
                  Dwt(1:Qdof, 1:ndim), Dwx(1:Qdof,1:ndim, 1:nbDim) )
          endif
          !write(*,'(a6,2i5,20es12.4)') 'xiR:',elemi, Qdof, xiR(:, 1)
          !write(*,'(a6,2i5,20es12.4)') 'xiR:',elemi, Qdof, xiR(:, 2)


          ! diffusive terms
          call Set_R_s(ndim, nbDim, iRe, Qdof, wi(1:Qdof,1:ndim), Dwx(1:Qdof, 1:ndim, 1:nbDim), Re_1,&
               R_s(1:Qdof, 1:nbDim, 1:ndim),  elem%xi(0,1:Qdof, 1:nbDim) )

          ! convective terms
          call Set_f_s(ndim, nbDim, Qdof, wi(1:Qdof,1:ndim), f_s(1:Qdof, 1:nbDim, 1:ndim), &
               elem%xi(0,1:Qdof, 1:nbDim), elem%i )

          ! total flux
          f_s(1:Qdof, 1:nbDim, 1:ndim)  = R_s(1:Qdof, 1:nbDim, 1:ndim) &
               - f_s(1:Qdof, 1:nbDim, 1:ndim)

          if(elemi == 1) then
             call EvalDiv_F(elem1, Qnum, Qdof, xiR(1:Qdof, 1:nbDim), &
                  f_s(1:Qdof, 1:nbDim, 1:ndim), DivF(1:Qdof, 1:ndim) )
          else
             call EvalDiv_F(elem2, Qnum, Qdof, xiR(1:Qdof, 1:nbDim), &
                  f_s(1:Qdof, 1:nbDim, 1:ndim), DivF(1:Qdof, 1:ndim) )
          endif

          do gi=1,Qdof  ! integ nodes
             call RHS_Scalar(xi(gi, 1:2), f, state%time%ctime)

             val_f = val_f + V_rule%weights(gi) * G_rule%weights(gk) &
                  * (f(1) - Dwt(gi, 1) +  DivF(gi, 1) )
                  !* 1.

             !if(i==1 .and. j==1 .and. k==1 .and. state%time%iter==1) &
             !     !write(53, '(15es12.4,11i5)') xi(gi, :), t, &
             !     write(53, '(7es12.4,11i5)') xi(gi, :), t, &
             !     val_f, V_rule%weights(gi) * G_rule%weights(gk),V_rule%weights(gi), G_rule%weights(gk),&
             !     !f(1), wi(gi, 1), Dwt(gi, 1), Dwx(gi, 1, 1:2), &
             !     !f_s(gi, 1:2, 1), DivF(gi, 1), &
             !     state%time%iter, i,j,k,i1,j1, l, gk, gi, elem1%i, elem2%i
          enddo  ! end of gi
          !!stop

       enddo ! end of gk

       !stop

    enddo ! end of l


    !if(i==1 .and. j==1 .and. k==1 .and. state%time%iter==1) &
    !write(56, '(4es12.4, 8i5)') xp(0,:), t1 + (k-0.5)/Nt * (t2 - t1), &
    !     val_f/4,  state%time%iter, i,j,k, i1,j1, elem1%i, elem2%i

    !val = val_f / 4   ! we have four subelements
    val = val_S  + val_f / 4   ! we have four subelements
    val2= val2_S + val_f / 4   ! we have four subelements

    val_N = 0.

    ! adding of initial condition
    if(k == 1 .and. state%time%iter == 1) then
       t = 0.  ! initial condition
       state%time%ctime = 0.

       if(elemi == 1) then
          call Eval_w_Qnode(elem1, Qdof, xiR(1:Qdof, 1:nbDim), t, wi(1:Qdof, 1:ndim), &
               Dwt(1:Qdof, 1:ndim), Dwx(1:Qdof,1:ndim, 1:nbDim) )

       else
          call Eval_w_Qnode(elem2, Qdof, xiR(1:Qdof, 1:nbDim), t, wi(1:Qdof, 1:ndim), &
               Dwt(1:Qdof, 1:ndim), Dwx(1:Qdof,1:ndim, 1:nbDim) )
       endif

       ! exact solution = IC
       do l=1,Qdof
          call Exact_Scalar( xi(l,1:nbDim), wE(l,1:ndim), state%time%ctime )
       enddo

       val_N =  dot_product(V_rule%weights(1:Qdof), wE(1:Qdof, 1) - wi(1:Qdof, 1) )

    end if

    deallocate (f, xi, xii, xiR, xiRR, wi, wE, Dwt, Dwx, wiR, DwxR, DivF)
    deallocate(iedges, itri)

  end subroutine EvalDualSubElem


  !> nodes lie on triangles edges
  subroutine NodesInTriangleEdges(elem1, elem2, Qdof, xi, xiR, inside, elemi, facei)
    type(element), intent(in) :: elem1, elem2
    integer, intent(in) :: Qdof
    real, dimension(1:Qdof, 1:2), intent(inout) :: xi        ! real coordinates
    real, dimension(1:Qdof, 1:3), intent(inout) :: xiR       ! relative (barycentric) coords
    logical, intent(inout) :: inside
    integer, intent(inout) :: elemi, facei
    real, dimension(:, :), allocatable :: xii
    integer :: j

    allocate( xii(1:Qdof, 1:3) )

    !!write(22, '(10es12.4)') xi(1, 1:2)

    ! barycentric coordinates with respect to elem1
    call NodeInTriangle(Qdof, grid%x(elem1%face(idx,1),1:2), grid%x(elem1%face(idx,2),1:2), &
         grid%x(elem1%face(idx,3),1:2), xi(1:Qdof, 1:2), xiR(1:Qdof, 1:3), inside)


    elemi = 1

    if(.not. inside) then
       ! barycentric coordinates with respect to elem2
       call NodeInTriangle(Qdof, grid%x(elem2%face(idx,1),1:2), grid%x(elem2%face(idx,2),1:2), &
            grid%x(elem2%face(idx,3),1:2), xi(1:Qdof, 1:2), xii(1:Qdof, 1:3), inside)

       !write(*, '(10es12.4)') xi(1, 1:2), xii(1, 1:2)

       if(.not. inside) then
          print*,'edge is not in elem1 and elem2'
          do j=1,Qdof
             write(50,*) xi(j, 1:2)
          enddo
          write(51,*)  grid%x(elem1%face(idx,1),1:2)
          write(51,*)  grid%x(elem1%face(idx,2),1:2)
          write(51,*)  grid%x(elem1%face(idx,3),1:2)
          write(51,*)  grid%x(elem1%face(idx,1),1:2)
          write(51,'(x)')
          write(51,*)  grid%x(elem2%face(idx,1),1:2)
          write(51,*)  grid%x(elem2%face(idx,2),1:2)
          write(51,*)  grid%x(elem2%face(idx,3),1:2)
          write(51,*)  grid%x(elem2%face(idx,1),1:2)
          stop
       endif
       xiR(1:Qdof, 1:3) = xii(1:Qdof, 1:3)
       elemi = 2
    endif


    ! we test the barycentric coordinates if nodes lie on an edge j
    inside = .false.
    do j=1,3
       if(sum(abs(xiR(1:Qdof,j) )) < 1E-5) then
          inside = .true.
          facei = j
          goto 100
       endif
    enddo
    facei = 0
100 continue

    deallocate( xii )

    !!write(*, '(2es12.4,i5,6es12.4)') xi(1, 1:2), elemi, xiR(1, 1:3)

  end subroutine NodesInTriangleEdges

  !> nodes are in a square formed by two triangles
  subroutine NodesInTwoTriangles(elem1, elem2, Qdof, xi, xiR, inside, elemi)
    type(element), intent(in) :: elem1, elem2
    integer, intent(in) :: Qdof
    real, dimension(1:Qdof, 1:2), intent(out) :: xi        ! real coordinates
    real, dimension(1:Qdof, 1:3), intent(out) :: xiR       ! relative (barycentric) coords
    logical, intent(inout) :: inside
    integer, intent(inout) :: elemi
    real, dimension(:, :), allocatable :: xii

    allocate( xii(1:Qdof, 1:3) )

    ! barycentric coordinates with respect to elem1
    call NodeInTriangle(Qdof, grid%x(elem1%face(idx,1),:), grid%x(elem1%face(idx,2),:), &
         grid%x(elem1%face(idx,3),:), xi(1:Qdof, 1:2), xiR(1:Qdof, 1:3), inside)

    elemi = 1
    if(.not. inside) then
       ! barycentric coordinates with respect to elem2
       call NodeInTriangle(Qdof, grid%x(elem2%face(idx,1),:), grid%x(elem2%face(idx,2),:), &
            grid%x(elem2%face(idx,3),:), xi(1:Qdof, 1:2), xii(1:Qdof, 1:3), inside)

       if(.not. inside) then
          print*,'sqaure is not in elem1 and elem2'
          stop
       endif
       xiR(1:Qdof, 1:3) = xii(1:Qdof, 1:3)
       elemi = 2
    endif

    deallocate(xii)

  end subroutine NodesInTwoTriangles

  !> nodes are in one triangle
  subroutine NodeInTriangle(Qdof, x1, x2, x3, xi, lambda, inside)
    integer, intent(in) :: Qdof
    real, dimension(1:2), intent(in) :: x1, x2, x3
    real, dimension(1:Qdof, 1:2), intent(out) :: xi         ! real coordinates
    real, dimension(1:Qdof, 1:3), intent(out) :: lambda     ! relative (barycentric) coords
    logical, intent(out) :: inside
    logical :: inside1
    real :: D, Dx, Dy
    integer :: i,j


    !print*,'-----------------------', Qdof

    D  = (x1(1) - x3(1) ) * (x2(2)  - x3(2) ) - (x1(2) - x3(2) ) * (x2(1)  - x3(1) )

    do i=1,Qdof
       Dx = (xi(i,1) - x3(1) ) * (x2(2)  - x3(2) ) - (xi(i,2) - x3(2) ) * (x2(1)  - x3(1) )
       Dy = (x1(1) - x3(1) ) * (xi(i,2)  - x3(2) ) - (x1(2) - x3(2) ) * (xi(i,1)  - x3(1) )

       lambda(i,3) = Dx / D
       lambda(i,1) = Dy / D
       lambda(i,2) = 1. - lambda(i,1) - lambda(i,3)


       if(minval(lambda(i,1:3) ) >= -1E-05 .and. maxval(lambda(i,1:3) ) <= 1 + 1E-05 )then
          inside = .true.
       else
          inside = .false.
       endif

       !print*,'???',i, minval(lambda(i,1:3) ), maxval(lambda(i,1:3) ),inside, inside1


       if(i ==1 ) then
          inside1 = inside
       else
          if(inside .neqv. inside1 ) then
             print*, ' trouble in subroutine NodeInTriangle ???'
             do j=1, i
                print*,'???',j, minval(lambda(j,1:3) ), maxval(lambda(j,1:3) ),inside, inside1
             enddo
             stop
          endif
       endif
    enddo
  end subroutine NodeInTriangle


  !> compute the flux error
  !> \f$ \sum_K (c_f \|u_h -u\| + c_K \|R_s(u_h) - R_s(u) \| + c_f \|f_s(u_h) - f_s\|)\f$
  subroutine ComputeFluxError(  )
    class(element), pointer :: elem
    integer :: i

    do i=1,grid%nelem
       elem => grid%elem(i)

       if( state%modelName == 'scalar' .or.state%modelName == '2eqs') then
          call ComputeFluxErrorElem(Set_R_s_scalar, Set_f_s_scalar, elem )

          state%estim(eN1:quadra,:) = state%estim(eN1:quadra,:) + elem%eta(eN1:quadra,:)
          state%estim( NC1n   ,:) = state%estim( NC1n   ,:) + elem%eta( NC1n   ,:)

          !write(*,'(a6, i5, 12es12.4)') 'err??',i, elem%eta(eN1:eN3p, 1), elem%eta(NC1n, 1)

       else
          print*,'Case nbDim ==',nbDim,' & ndim ==',ndim,' not implemented in errorFlux.f90'
          stop
       endif
    enddo

  end subroutine ComputeFluxError

  !> compute the flux error for element
  !> \f$ c_f \|u_h -u\| + c_K \|R_s(u_h) - R_s(u) \| + c_f \|f_s(u_h) - f_s\| \f$
  subroutine ComputeFluxErrorElem(Set_R_s, Set_f_s, elem )
    interface
      subroutine Set_R_s(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
         integer, intent(in) :: ndimL, nbDim, iRe, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1    ! inverse of Reynolds number
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
    type(element), intent(inout) :: elem
    type(Gauss_rule), pointer :: G_rule
    real, dimension(:,:,:), allocatable :: wi
    real, dimension(:,:,:,:), allocatable :: Dwi
    real, dimension(:,:), allocatable :: wE, x, func, funcS
    real, dimension(:,:,:), allocatable :: DwE, R_s, R_sE, f_s, f_sE
    real, dimension(:,:), allocatable :: Re_1
    real, dimension(:), allocatable :: val
    real :: errU, errR, errF, errFR, errQ

    real :: t, CKb, hF
    integer :: Gnum, Gdof, Qnum, Qdof, j, l, i

    Qdof = elem%Qdof
    Qnum = elem%Qnum

    allocate(wi(1:3, 1:Qdof, 1:ndim) )
    call Eval_w_Elem(elem, wi(1, 1:Qdof, 1:ndim) )
    call Eval_w_Elem_time(elem, wi(2, 1:Qdof, 1:ndim) )

    allocate(Dwi(1:3, 1:Qdof, 1:ndim, 1:nbDim) )
    call Eval_Dw_Elem(elem, Dwi(1, 1:Qdof, 1:ndim, 1:nbDim) )
    call Eval_Dw_Elem_time(elem, Dwi(2, 1:Qdof, 1:ndim, 1:nbDim) )

    ! exact solutions
    allocate(wE(1:Qdof, 1:ndim) )
    allocate(DwE(1:Qdof, 1:ndim, 1:nbDim) )

    allocate(R_s(1:Qdof, 1:nbDim, 1:ndim),  f_s(1:Qdof, 1:nbDim, 1:ndim) )
    allocate(R_sE(1:Qdof, 1:nbDim, 1:ndim), f_sE(1:Qdof, 1:nbDim, 1:ndim) )

    allocate( x(1:Qdof,1:nbDim) ) ! integ nodes
    call ComputeF(elem, Qdof, state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim), x(1:Qdof, 1:nbDim) )

    allocate(func(1:Qdof, 1:ndim), funcS(1:Qdof, 1:ndim), val(1:ndim)  )

    allocate(Re_1(1:iRe, 1:Qdof) )
    Re_1(1:iRe, 1:Qdof) = state%model%Re1

    Gnum = 3   ! 15 is the maximal one !!!!!!!!!
    G_rule => state%space%G_rule(Gnum)
    Gdof = G_rule%Qdof

    errU = 0.
    errR = 0.
    errF = 0.
    errFR = 0.
    errQ = 0.

    do j=1, Gdof
       ! actual time
       t = G_rule%lambda(j)
       state%time%ctime = state%time%ttime - state%time%tau(1) * t
       state%time%ctime = state%time%ttime - state%time%tau(1) * (1 - t ) !!!!!!!!!!!!

       ! exact solution at time %ctime
       do l=1,Qdof
          call Exact_Scalar( x(l,1:nbDim), wE(l,1:ndim), state%time%ctime )
          call Der_Exact_Scalar( x(l,1:nbDim), DwE(l,1:ndim,1:nbDim), state%time%ctime )
       enddo

       ! approximate solution at time %ctime
       wi(3, 1:Qdof, 1:ndim) = t * wi(1, 1:Qdof, 1:ndim)  + (1.-t)*wi(2, 1:Qdof, 1:ndim)
       Dwi(3, 1:Qdof, :, :) = t * Dwi(1, 1:Qdof, :, :)  + (1.-t)*Dwi(2, 1:Qdof, :, :)

       ! SOLUTION TERM
       func(1:Qdof, 1:ndim) = wi(3, 1:Qdof, 1:ndim) - wE( 1:Qdof, 1:ndim)
       call IntegrateSquareVectorFunction2(elem, func, val)
       errU = errU +  G_rule%weights(j) * sum(val(:) )

       !write(71,'(3i5,30es12.4)') state%time%iter, elem%i, j, wi(3, 1:Qdof, 1:ndim)
       !write(71,'(3i5,30es12.4)') state%time%iter, elem%i, j, wE(   1:Qdof, 1:ndim)
       !write(71,'(3i5,30es12.4)') state%time%iter, elem%i, j, func(1:Qdof, 1:ndim), val
       !write(71, *) '----'

       ! DIFFUSIVE TERMS
       ! exact value
       call Set_R_s(ndim, nbDim, iRe, Qdof, wE(1:Qdof,1:ndim), DwE(1:Qdof, 1:ndim, 1:nbDim), Re_1,&
            R_sE(1:Qdof, 1:nbDim, 1:ndim), elem%xi(0,1:Qdof, 1:nbDim) )

       ! appoximate value
       call Set_R_s(ndim, nbDim, iRe, Qdof, wi(3, 1:Qdof,1:ndim), Dwi(3, 1:Qdof, 1:ndim, 1:nbDim), Re_1,&
            R_s(1:Qdof, 1:nbDim, 1:ndim), elem%xi(0,1:Qdof, 1:nbDim) )


       func(:,:) = 0
       do i=1,nbDim
          func(1:Qdof, 1:ndim) = func(1:Qdof, 1:ndim) &
               + ( R_s(1:Qdof, i, 1:ndim) - R_sE(1:Qdof, i, 1:ndim) )**2
       enddo
       call IntegrateVectorFunction2(elem, func, val)
       errR = errR +  G_rule%weights(j) * sum(val(:) )


       ! CONVECTIVE TERMS
       ! exact value
       call Set_f_s(ndim, nbDim, Qdof, wE(1:Qdof,1:ndim), f_sE(1:Qdof, 1:nbDim, 1:ndim), &
            elem%xi(0,1:Qdof, 1:nbDim), elem%i )

       ! appoximate value
       call Set_f_s(ndim, nbDim, Qdof, wi(3, 1:Qdof,1:ndim), f_s(1:Qdof, 1:nbDim, 1:ndim), &
            elem%xi(0,1:Qdof, 1:nbDim), elem%i )

       func(:,:) = 0.
       funcS(:,:) = 0.
       do i=1,nbDim
          func(1:Qdof, 1:ndim) = func(1:Qdof, 1:ndim) &
               + ( f_s(1:Qdof, i, 1:ndim) - f_sE(1:Qdof, i, 1:ndim) )**2

          funcS(1:Qdof, 1:ndim) = funcS(1:Qdof, 1:ndim) &
               + ( f_s(1:Qdof, i, 1:ndim) - f_sE(1:Qdof, i, 1:ndim) &
               - R_s(1:Qdof, i, 1:ndim) + R_sE(1:Qdof, i, 1:ndim) )**2
       enddo
       call IntegrateVectorFunction2(elem, func, val)
       errF = errF +  G_rule%weights(j) * sum(val(:) )


       funcS(:,:) = funcS(:, :) + func(:, :)
       call IntegrateVectorFunction2(elem, funcS, val)
       errFR = errFR +  G_rule%weights(j) * sum(val(:) )

       !write(*,'(a6,2i3,20es11.4)') 'wi:',elem%i, l, wi(3,:,1)
       !write(*,'(a6,2i3,20es11.4)') 'err*:',elem%i, l, errU, errR, errF

       ! quadrature errors
       ! f_s contains viscous and inviscid fluxes in integ nodes
       f_s(1:Qdof, 1:nbDim, 1:ndim) = f_s(1:Qdof, 1:nbDim, 1:ndim) - R_s(1:Qdof, 1:nbDim, 1:ndim)

       call ElementPolynProj(elem, f_s(1:Qdof, 1:nbDim, 1:ndim), f_sE(1:Qdof, 1:nbDim, 1:ndim) )
       func(:,:) = 0.


       func(:,:) = 0.
       do i=1,nbDim
          func(1:Qdof, 1:ndim) = func(1:Qdof, 1:ndim) &
               + ( f_s(1:Qdof, i, 1:ndim) - f_sE(1:Qdof, i, 1:ndim) )**2
       enddo
       call IntegrateVectorFunction2(elem, func, val)
       errQ = errQ +  G_rule%weights(j) * sum(val(:) )

       !write(61,'(a6,2i5,50es12.4)')'A',elem%i, Qdof, f_s(1:Qdof, 1:nbDim, 1:ndim)
       !write(61,'(a6,2i5,50es12.4)')'P',elem%i, Qdof, f_sE(1:Qdof, 1:nbDim, 1:ndim)
       !write(61,'(a6,2i5,50es12.4)')'D',elem%i, Qdof, f_s(1:Qdof, 1:nbDim, 1:ndim)-f_sE(1:Qdof, 1:nbDim, 1:ndim)
       !!write(*,'(a4,i5,20es12.4)' ) "eeF", elem%i, errU,errF,errR,errFR,errQ

    enddo  ! end of j=1,Gdof  integration over (t_{k-1}, t_k)

    ! length of the time interval
    errU  = errU  * state%time%tau(1)
    errF  = errF  * state%time%tau(1)
    errR  = errR  * state%time%tau(1)
    errFR = errFR * state%time%tau(1)
    errQ = errQ * state%time%tau(1)

    !write(62, *) elem%xc(:), errQ

    !elem%CTn = (elem%CK + elem%Cb)/elem%diam**2 + state%time%FinTime/state%time%tau(1)**2

    ! new variant (revision for SIAM), constant value of CTn
    if(state%time%iter_loc == 1) then
       elem%CTn = maxval(grid%elem(:)%CK + grid%elem(:)%Cb, 1) /  &
            grid%elem(1)%diam**2 + state%time%FinTime/state%time%tau(1)**2
    endif


    ! complete variant
    elem%eta(eN1:quadra , :) = 0.

    elem%eta(eN1,:) = errU  / (elem%CTn * state%time%tau(1)**2)
    elem%eta(eN2,:) = errFR / (elem%CTn * elem%diam**2)

    elem%eta(eN1p,:) =  errU / state%time%FinTime
    if(elem%CK >0) elem%eta(eN2p,:) =  errR / elem%CK
    if(elem%Cb > 0) elem%eta(eN3p,:) =  errF / elem%Cb

    elem%eta(quadra,:) =  errQ / (elem%CTn * elem%diam**2)


    ! nonconformity term
    elem%eta(NC1n, 1:ndim) = 0.

    do i = 1, elem%flen
       ! jumps of the solution
       call IntegElementEdgeJump_time(elem, i, val(1:ndim) )

       if(i <= 5) print*,'ATTENTION in errorFLUX.f90, penalty parameters depends on elem%deg'
       hF = elem%dn(i)
       CKb = elem%CKo**2 / hF + hF * elem%Cbo**2 + (state%model%Re1*state%space%sigma)**2 / hF
       CKb = CKb / elem%CTn / (elem%diam**2)

       elem%eta(NC1n, 1:ndim) = elem%eta(NC1n, 1:ndim) + CKb * val(1:ndim)

    enddo

    if(elem%i == 1 .and. state%time%iter == 1) then
       write(*,'(2(a6,es12.4))') 'h_K=',elem%diam, '\tau=',state%time%tau(1)
       write(*,'(4(a6,es12.4))') 'CK =',elem%CK,'Cb =',elem%Cb, &
            'CTF =',state%time%FinTime/state%time%tau(1)**2, 'CTn =',elem%CTn
       write(*,'(5(a6,es12.4))') &
            'CF1 =', (state%model%Re1*state%space%sigma)**2 / hF, &
            'CF2 =', elem%CKo**2 / hF, &
            'CF3 =', hF * elem%Cbo**2 , &
            'CF = ', &
            elem%CKo**2 / hF + hF * elem%Cbo**2 + (state%model%Re1*state%space%sigma)**2 / hF, &
            'CFt=',CKb
       write(*,'(a6,4es12.4)') 'rat:', 1. / elem%CTn / (elem%diam**2), &
            elem%CTn,elem%diam,  elem%diam**2
    endif


    deallocate(wi, Dwi, wE, DwE, x, func, funcS, R_s, R_sE, f_s, f_sE, val, Re_1 )

  end subroutine ComputeFluxErrorElem

end module errorFlux
