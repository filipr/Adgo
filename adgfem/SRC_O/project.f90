!> general error estimation subroutines
module project_estimation

  use main_data  ! contains type(mesh) :: grid for computation
  use paramets
  use problem_oper
  use euler_problem
  use elemental_mod

  implicit none

  public:: DualElemEstimate
  public:: DWRElemEstim_space
  public:: DWRElemEstim_alg
  public:: ReziduumElemEstimate
  public:: ST_DualElemEstimate
  public:: ST_DualElemEstimate_Var2
  public:: ProjectionElemEstimate
  public:: ElementRegularityEstimator
contains

  !> estimation of the element reziduum,
  !> \f$ \eta_K := \int_{K} \nabla \cdot \vec{F}(w_h, \nabla u_h)\, d x =
  !> \int_{\partial K} \vec{F}(w_h, \nabla w_h) \cdot \vec{n}\, dS \f$
  subroutine ReziduumElemEstimate(elem)
    type(element) :: elem

    if(state%modelName == 'NSe') then
       call ComputeOneElementFluxes(elem, Set_f_s_Euler, Set_R_s_NS)

    elseif(state%modelName == 'scalar' .or.state%modelName == '2eqs') then
       call ComputeOneElementFluxes(elem, Set_f_s_scalar, Set_R_s_scalar)

    else
       print*,'Stopped (22) in project.f90'
       stop
    endif

  end subroutine ReziduumElemEstimate




  !> estimation of the element reziduum,
  !> \f$ \eta_K := \max_{\| \varphi\|_X=1}
  !> \frac{|c_h({\bf w}, \varphi_i) |}{\|\varphi\|_X}\f$,
  !> \f$\varphi\in P^{p_K+1}(K)\f$; \
  !> ityp = 1  \f$ \Rightarrow \| v \|_X = \| v \|_{L^2} \f$,
  !> ityp = 2  \f$ \Rightarrow \| v \|_X = \| v \|_{H^1} \f$,
  !> ityp = 3  \f$ \Rightarrow \| v \|_X^2 = \|u\|_0^2 + \varepsilon |u|_1^2 \f$,
  !> ityp = 4  \f$ \Rightarrow \| v \|_X^2 = \varepsilon ( |u|_1^2 + J_h^{\sigma}(u,u) )\f$,
  !> space-time DG, time basis is ORTHOGONAL
  subroutine ST_DualElemEstimate(elem, ityp )
    type(element) :: elem
    integer, intent(in) :: ityp
    real, dimension(:,:), allocatable :: SS, SS1
    real, dimension(:), allocatable :: b, ident
    real, dimension(:, :), allocatable :: func
    real :: est1, est2
    integer :: dof, Tdof, dofP, TdofP, Qdof
    integer :: i,j,k,l
    logical :: iprint
    real ::val1, val2, val3

    !FR FIXME - control for all elements should be connected earlier
    associate( time => state%time )
    select type ( time )
    class is (TimeTDG_t)

       !iprint = .true.
       iprint = .false.
       !if(elem%i >= 122 .and. elem%i <= 122) iprint = .true.


       dof  = elem%dof
       Tdof = elem%Tdof
       dofP  = elem%dof_plus
       TdofP = elem%Tdof_plus

       Qdof = elem%Qdof

       allocate(  SS1(1:dofP, 1:dofP) )
       SS1(:,:)=0.

       !if(elem%i == 1 .and. state%time%iter == 1) print*,'subroutine ST_DualElemEstimate is not READY'

       !if(iprint) then
       !!if(elem%rhsST(1,1,1) > 1.E-5) then
       !   write(*,'(a3,6i5)') 'ELM',elem%i,dof, dofP, Tdof, TdofP
       !
       !

       !endif

       if(ityp >= 2) then  !print*,'semi-H_1 norm'
          allocate(ident(1:Qdof) )
          ident(:) = 1.
          call IntegrateBlockD2(elem, dofP, ident(1:Qdof), SS1(1:dofP, 1:dofP) )
          ! print*,'
          deallocate(ident)
       endif

       if(ityp == 3) then  !! "normalization: |u|_1 -> \sqrt{\varepsilon}| u |_1
          SS1(1:dofP, 1:dofP) =  SS1(1:dofP, 1:dofP) * state%model%Re1

          if(state%modelName == 'scalar'  &
               .and. (state%model%idiff == 9 .or. state%model%idiff == 10) ) then

             if(state%model%idiff == 9) then
                ! battery, diffusion is not constant
                call Set_Battery(elem%xc(1), elem%xc(2), val1, val2, val3)

             elseif( state%model%idiff == 10) then
                ! simplified battery, diffusion is not constant
                call Set_Battery_Simplified(elem%xc(1), elem%xc(2), val1, val2, val3)
             endif

             val3 = (val1 + val2) / 2

             SS1(1:dofP, 1:dofP) =  SS1(1:dofP, 1:dofP) * val3
             ! allocate(func(1:2, 1:elem%Qdof) )

             ! func(1, 1:elem%Qdof) = sqrt(val1)
             ! func(2, 1:elem%Qdof) = sqrt(val2)

             ! SS1(:,:) = 0.
             ! call IntegrateBlockD2F2(elem, dofP, func(1:2, 1:elem%Qdof), SS1(1:dofP, 1:dofP) )
             ! deallocate(func)
          endif

       endif

       ! adding of the L2 norm
       if(ityp >= 1 .and. ityp <= 3) then
          do i=1,dofP
             SS1(i,i) = SS1(i,i) + 2*elem%area !*(2**0.5) !! sqrt(2) is the "size of convection"
          enddo
       endif

       ! storing
       !!!!SS1(1:dofP, 1:dofP) = SS(1:dofP, 1:dofP)

       ! time-dependent test functions, L^2-orthogonal
       !SS(1:dofP, 1:dofP) = SS(1:dofP, 1:dofP) * state%time%tau(1)
       !SS(1:dofP, 1:dofP) = SS(1:dofP, 1:dofP) / state%time%tau(1)**2  ! approximation of H^1-norm (I_m)

       !if(iprint) then
       ! if(elem%i <= 1) then
       !    do l=1, dofP
       !       write(*,'(a6,i5,30es12.4)')'@@P',l, SS1(l, 1:dofP)
       !    enddo
       ! endif

       elem%eta(resA, 1) = 0.
       elem%eta(resS, 1) = 0.
       elem%eta(resT, 1) = 0.
       elem%eta(resST, 1) = 0.

       allocate( b (1:ndim*dofP),  SS(1:dofP, 1:dofP) )

       ! we USE the orthogonality of the time basis, we go over time levels
       do i=1, TdofP
          do k=1, ndim
             b( (k-1)*dofP + 1: k*dofP) = elem%rhsST(k, 1:dofP, i) * state%time%tau(1) ! scaling ?
          enddo

          ! multiplication by the norm of the time DG basis functions
          !b( 1: ndim*dofP) = b( 1: ndim*dofP) *(1.+ state%time%StiffTimeMatrix%Mb(i,i)) *state%time%tau(1)
          !b( 1: ndim*dofP) = b( 1: ndim*dofP) *(state%time%tau(1)+ state%time%StiffTimeMatrix%Mb(i,i) /state%time%tau(1) )

          ! L^2 -norm with respect to the time
          !SS(1:dofP, 1:dofP) = SS1(1:dofP, 1:dofP) * state%time%tau(1)

          ! H^1 -norm with respect to the time
          SS(1:dofP, 1:dofP) = SS1(1:dofP, 1:dofP) &
                *(time%tau(1)+ time%StiffTimeMatrix%Mb(i,i) /time%tau(1) )

          !if(elem%i <= 1) then
          !   !do k=1,dofP
          !   !   write(*,'(a10,2i5, 300es12.4)') 'b, SS:',i,k,b(k), SS(k, 1:dofP)
          !   !enddo
          !   do k=1,dofP
          !      write(*,'(a10,2i5, 300es12.4)') 'b, SS1:',i,k,b(k), SS1(k, 1:dofP)
          !   enddo
          !   print*,'3ed3434ds', time%tau(1),i, time%StiffTimeMatrix%Mb(i,i)
          !endif


          if(i <=  Tdof) then
             ! algebraic error
             call EvalMaximumLang(dofP, 1, dof,  SS(1:dofP, 1:dofP), b(1:dofP*ndim), est1)

             ! space error
             call EvalMaximumLang(dofP, 1, dofP, SS(1:dofP, 1:dofP), b(1:dofP*ndim), est2)

             elem%eta(resA, 1) = elem%eta(resA, 1) + est1**2
             elem%eta(resS, 1) = elem%eta(resS, 1) + est2**2

             elem%eta(resT, 1)  = elem%eta(resT, 1)  + est1**2
             elem%eta(resST, 1) = elem%eta(resST, 1) + est2**2

          else if(i == TdofP) then
             ! time error
             call EvalMaximumLang(dofP, 1, dof,  SS(1:dofP, 1:dofP), b(1:dofP*ndim), est1)

             ! space-time error
             call EvalMaximumLang(dofP, 1, dofP, SS(1:dofP, 1:dofP), b(1:dofP*ndim), est2)

             elem%eta(resT, 1)  = elem%eta(resT, 1)  + est1**2
             elem%eta(resST, 1) = elem%eta(resST, 1) + est2**2
          endif
          !if(iprint) &
          !write(*,'(a6,i5,30es12.4)')'EST',elem%i, elem%eta(resA, 1), elem%eta(resS, 1), elem%eta(resT, 1), elem%eta(resST, 1)
       enddo

       !if(elem%i == 1104) &
       !      write(*,'(a6,i5,30es12.4)')'EST',elem%i, elem%eta(resA, 1), elem%eta(resS, 1), elem%eta(resT, 1), elem%eta(resST, 1)

       elem%eta(resA,  1) = sqrt(elem%eta(resA, 1) )
       elem%eta(resS,  1) = sqrt(elem%eta(resS, 1) )
       elem%eta(resT,  1) = sqrt(elem%eta(resT, 1) )
       elem%eta(resST, 1) = sqrt(elem%eta(resST, 1) )

       elem%eta(resSr, 1) = elem%eta(resS, 1) ! NOT USED !!!!!


       if(iprint) then
          print*
          write(*,'(a6,i5,30es12.4)')'EST',elem%i, elem%eta(resA, 1), elem%eta(resS, 1), elem%eta(resT, 1), elem%eta(resST, 1)
       endif

       !if(elem%i > 10) stop

       deallocate( b, SS)

      class default
      stop 'STDG only'
   end select
   end associate

  end subroutine ST_DualElemEstimate


!> estimation of the element reziduum,
  !> \f$ \eta_K := \max_{\| \varphi\|_X=1}
  !> \frac{|c_h({\bf w}, \varphi_i) |}{\|\varphi\|_X}\f$,
  !> \f$\varphi\in P^{p_K+1}(K)\f$; \
  !> ityp = 1  \f$ \Rightarrow \| v \|_X = \| v \|_{L^2} \f$,
  !> ityp = 2  \f$ \Rightarrow \| v \|_X = \| v \|_{H^1} \f$,
  !> ityp = 3  \f$ \Rightarrow \| v \|_X^2 = \|u\|_0^2 + \varepsilon |u|_1^2 \f$,
  !> ityp = 4  \f$ \Rightarrow \| v \|_X^2 = \varepsilon ( |u|_1^2 + J_h^{\sigma}(u,u) )\f$,
  !> space-time DG, time basis is ORTHOGONAL

!ityp = 1 L2(0,T;L2) ; ityp = 2  L2(0,T;H1) ; ityp = 3 H1(0,T;H1-element ortogonality) ; ityp = 4 H1(0,T;H1-DG)
  subroutine ST_DualElemEstimate_Var2(elem, ityp)
    type(element) :: elem
    integer, intent(in) :: ityp
    real, dimension(:,:), allocatable :: S, SP, Sh, St,Imatrix
    real, dimension(:), allocatable :: b, bP, ident
    real :: estA, estT, estS, estST
    integer :: dof, dofP, Tdof, TdofP, Qdof,dof_tot, dof_AT
    integer :: i,j, k, pombP,pomb, k1,k2, l1,l2, m1,m2, n1,n2
    logical :: iprint

        !FR FIXME - control for all elements should be connected earlier
    associate( time => state%time )
    select type ( time )
      class is (TimeTDG_t)


       iprint = .false.

       dof  = elem%dof
       dofP  = elem%dof_plus
       Tdof = elem%Tdof
       TdofP = elem%Tdof_plus

       Qdof = elem%Qdof

       dof_AT = TdofP*dof
       dof_tot = TdofP*dofP

       allocate(Sh(1:dofP, 1:dofP)) !storing X-norm respect to space
       allocate( b(1:ndim*dof_AT), S(1:dof_AT,1:dof_AT) )
       allocate( bP(1:ndim*dof_tot),  SP(1:dof_tot, 1:dof_tot) )

       Sh(:,:)=0.
       S(:,:) = 0.
       SP(:,:) = 0.
       b(:) = 0.
       bP(:) = 0.


       if(ityp == 4) then
          !print*,' penalty'
          call IntegrateBlockJumps(elem, dofP, Sh(1:dofP, 1:dofP) )
       endif

       if(ityp >= 2) then  !'semi-H_1 norm'
          allocate(ident(1:Qdof) )
          ident(:) = 1.
          call IntegrateBlockD2(elem, dofP, ident(1:Qdof), Sh(1:dofP, 1:dofP) )
          deallocate(ident)
       endif

       if((ityp==1).or.(ityp==2)) then  ! add L_2 norm. X = L2 or X = H1
          do i = 1,dofP
             Sh(i,i) = Sh(i,i) + 2*elem%area
          enddo
       else if(ityp == 3) then  ! X = H1-element ortogonality
          ! "normalization: |u|_1 -> \sqrt{\varepsilon}| u |_1
          Sh(1:dofP, 1:dofP) =  Sh(1:dofP, 1:dofP) * state%model%Re1
          ! adding of the L2 norm
          do i=1,dofP
             Sh(i,i) = Sh(i,i) + 2*elem%area
          enddo
       else if(ityp == 4) then  ! "normalization: |u|_1 -> \sqrt{\varepsilon}| u |_1
          if(state%model%Re > 0.) Sh(1:dofP, 1:dofP) =  Sh(1:dofP, 1:dofP) / state%model%Re
       endif

       do k = 1, ndim
         pombP = (k-1)*dof_tot
         pomb = (k-1)*dof_AT
         do i = 1, TdofP
             k1 = pombP  + (i-1)*dofP + 1
             k2 = k1 + dofP - 1 != (k-1)*dof_tot + i*dofP
             bP(k1:k2) = elem%rhsST(k, 1:dofP, i) * time%tau(1)

             l1 = pomb + (i-1)*dof + 1
             l2 = l1 + dof - 1 != (k-1)*dof_AT + i*dof
             b(l1:l2) = bP(k1:k1+dof-1)
         enddo
       enddo


       if((ityp==3).or.(ityp==4)) then  ! H1(0,T;X)

         ! storing H1-norm respect to time
         allocate(St(1:TdofP, 1:TdofP), Imatrix(1:TdofP, 1:TdofP))
         St(:,:)=0.
         Imatrix(:,:)=0. !identity matrix
         do i = 1,TdofP
            Imatrix(i,i) = 1.
         enddo
         St(1:TdofP, 1:TdofP) = time%tau(1)*Imatrix(1:TdofP, 1:TdofP) &
              + time%StiffTimeMatrix%Mb(1:TdofP, 1:TdofP) /time%tau(1)

         do i = 1,TdofP
           k1 = (i-1)*dof + 1
           k2 = k1 + dof - 1  !i*dof
           m1 = (i-1)*dofP + 1
           m2 = m1 + dofP - 1 !i*dofP
     ! goto 55
           do j = i, TdofP
              l1 = (j-1)*dof + 1
              l2 = l1 + dof - 1   !j*dof
              n1 = (j-1)*dofP + 1
              n2 = n1 + dofP - 1  !j*dofP
               ! H1(0,T;X)-norm with respect to the time, p,q+1
              S(k1:k2,l1:l2) = St(i,j)*Sh(1:dof,1:dof)
               ! H1(0,T;X)-norm with respect to the time, p+1,q+1
              SP(m1:m2,n1:n2) = St(i,j)*Sh(1:dofP,1:dofP)
              if (j /= i) then
                 S(l1:l2,k1:k2) = S(k1:k2,l1:l2)
                 SP(n1:n2,m1:m2) = SP(m1:m2,n1:n2)
              endif
           enddo
     ! 55
     !     SS(k1:k2,k1:k2) = SSt(i,i)*SS1(1:dof,1:dof)
     !     SSP(m1:m2,m1:m2) = SSt(i,i)*SS1(1:dofP,1:dofP)
         enddo
         deallocate(Imatrix,St)

       else if((ityp==1).or.(ityp==2)) then  ! L2(0,T;X)
         do i = 1,TdofP
           k1 = (i-1)*dof + 1
           k2 = k1 + dof - 1
           m1 = (i-1)*dofP + 1
           m2 = m1 + dofP - 1
             ! L2(0,T;X) norm with respect to the time, S_{hp,t(q+1)}
           S(k1:k2,k1:k2) = time%tau(1)*Sh(1:dof,1:dof)
             ! L2(0,T;X)-norm with respect to the time, S_{h(p+1),t(q+1)}
           SP(m1:m2,m1:m2) = time%tau(1)*Sh(1:dofP,1:dofP)
         enddo
       endif

       deallocate(Sh)

       elem%eta(resA, 1) = 0.
       elem%eta(resS, 1) = 0.
       elem%eta(resT, 1) = 0.
       elem%eta(resST, 1) = 0.


        ! algebraic error
       call EvalMaximumLang(dof_AT, 1, dof*Tdof,  S(1:dof_AT,1:dof_AT), b(1:dof_AT*ndim), estA)
       elem%eta(resA, 1) = estA

        ! time error
       call EvalMaximumLang(dof_AT, 1, dof*TdofP,  S(1:dof_AT,1:dof_AT), b(1:dof_AT*ndim), estT)
       elem%eta(resT, 1) = estT

        ! space error
       call EvalMaximumLang(dof_tot, 1, dofP*Tdof, SP(1:dof_tot, 1:dof_tot), bP(1:dof_tot*ndim), estS)
       elem%eta(resS, 1) = estS

        ! space-time error
       call EvalMaximumLang(dof_tot, 1, dof_tot, SP(1:dof_tot, 1:dof_tot), bP(1:dof_tot*ndim), estST)
       elem%eta(resST, 1) = estST

        if(iprint) then
          print*
          write(*,'(a6,i5,30es12.4)')'EST',elem%i, elem%eta(resA, 1), elem%eta(resS, 1), elem%eta(resT, 1), elem%eta(resST, 1)
       endif

       deallocate(b, bP, S, SP)

      class default
      stop 'STDG only'
   end select
   end associate


  end subroutine ST_DualElemEstimate_Var2



  !> estimation of the element reziduum,
  !> \f$ \eta_K := \max_{\| \varphi\|_X=1}
  !> \frac{|c_h({\bf w}, \varphi_i) |}{\|\varphi\|_X}\f$,
  !> \f$\varphi\in P^{p_K+1}(K)\f$; \
  !> ityp = 1  \f$ \Rightarrow \| v \|_X = \| v \|_{L^2} \f$,
  !> ityp = 2  \f$ \Rightarrow \| v \|_X = \| v \|_{H^1} \f$,
  !> ityp = 3  \f$ \Rightarrow \| v \|_X^2 = \|u\|_0^2 + \varepsilon |u|_1^2 \f$,
  !> ityp = 4  \f$ \Rightarrow \| v \|_X^2 = \varepsilon ( |u|_1^2 + J_h^{\sigma}(u,u) )\f$,
  subroutine DualElemEstimate(elem, ityp, onlyAS )
    type(element) :: elem
    integer, intent(in) :: ityp
    logical, intent(in) :: onlyAS   ! only space and algebraic estimates
    real, dimension(:,:), allocatable :: S
    real, dimension(:), allocatable :: ident
    real, dimension(:, :), allocatable :: func
    real, dimension(:,:), pointer:: phi ! local store arrays
    integer :: i, k, dofA, dof1, dof, dofS, itest
    real ::val1, val2, val3

    !itest = 18
    !itest = 6

    dof = elem%dof
    dofA = elem%dof_plus
    dofS = elem%dof + elem%deg + 2

    !print*,'###', elem%i, dof, dofA, dofS

    dof1 = 1

    elem%eta(resA, 1)  = 0.
    elem%eta(resS, 1)  = 0.
    elem%eta(resT, 1)  = 0.
    elem%eta(resST, 1) = 0.
    elem%eta(resSr, 1)  = 0.

    !if(elem%i == 1) &
    !     write(*,'(a13,l2,4i5,200es12.4)') '## DualElemEsti',onlyAS, elem%i,dof,dofA, size(elem%vec,2) !,elem%vec(rhs,:)

    ! evaluation od the "Stiff" matrix in the appropriate  norm
    ! S(i,j) = (( \phi_i, \phi_j)),  ((. , . )) generate the norm
    allocate(S(1:dofA, 1:dofA) )
    S(:,:) = 0.

    if(ityp == 4) then
       !print*,' penalty'
       call IntegrateBlockJumps(elem, dofA, S(1:dofA, 1:dofA) )
    endif

    if(ityp >= 2) then  !print*,'semi-H_1 norm'
       allocate(ident(1:elem%Qdof) )
       ident(:) = 1.
       call IntegrateBlockD2(elem, dofA, ident(:), S(1:dofA, 1:dofA) )
       ! print*,'
       deallocate(ident)
    endif

    if(ityp == 3) then  !! "normalization: |u|_1 -> \sqrt{\varepsilon}| u |_1
       S(1:dofA, 1:dofA) =  S(1:dofA, 1:dofA) * state%model%Re1

       !do i=1, dofA
       !   write(*,'(a6, i5, 400es12.4)') 'matr:', i, S(i, :)
       !enddo

       ! THE FOLLOWING MODIFICATIONS DOES NOIT WORK, WHY ???
       if(state%modelName == 'scalar' .and. (state%model%idiff == 9 .or. state%model%idiff== 10))then


          if(state%model%idiff == 9) then
             ! battery, diffusion is not constant
             call Set_Battery(elem%xc(1), elem%xc(2), val1, val2, val3)

          elseif( state%model%idiff == 10) then
             ! simplified battery, diffusion is not constant
             call Set_Battery_Simplified(elem%xc(1), elem%xc(2), val1, val2, val3)
          endif

          val3 = (val1 + val2) / 2

          S(1:dofA, 1:dofA) =  S(1:dofA, 1:dofA) * val3

          !allocate(func(1:2, 1:elem%Qdof) )
          !func(1, 1:elem%Qdof) = sqrt( (val1 + val2 )/2 )
          !func(2, 1:elem%Qdof) = sqrt( (val1 + val2 )/2 )
          !S = 0.
          !call IntegrateBlockD2F2(elem, dofA, func(1:2, 1:elem%Qdof), S(1:dofA, 1:dofA) )

          !deallocate(func)
       endif

       !print*
       !do i=1, dofA
       !   write(*,'(a6, i5, 400es12.4)') 'matr:', i, S(i, :)
       !enddo
       !stop "====================================== jf943jf3oe"

       ! adding of the L2 norm
       do i=dof1,dofA
          S(i,i) = S(i,i) + 2*elem%area !*(2**0.5) !! sqrt(2) is the "size of convection"
       enddo
    endif

    if(ityp == 1 .or. ityp == 2 ) then  !! L_2 norm
       do i=dof1,dofA
          S(i,i) = S(i,i) + 2*elem%area
       enddo
    endif

    if(ityp == 4) then  !! "normalization: |u|_1 -> \sqrt{\varepsilon}| u |_1
       if(state%model%Re > 0.) S(1:dofA, 1:dofA) =  S(1:dofA, 1:dofA) / state%model%Re
    endif

    if(ityp == 5) then  !!???????? norm
       do i=dof1,dofA
          S(i,i) = 1.
       enddo
    endif

    !do i=1,dofA
    !    write(*,'(a4,i4,30es12.4)') ':::;',i, S(i, 1:dofA)
    ! enddo


    ! evaluation of the corresponding maxima over a discrete sets
    ! print*,"algebraic error"
    call EvalMaximumLang(dofA, dof1, dof, S(1:dofA, 1:dofA), &
         elem%vec(rhs, 1:dofA*ndim), elem%eta(resA, 1))
    ! Verfurth - not normalization,centrally done in SolveProblem1
    !!elem%eta(resA, 1) = elem%eta(resA, 1) * (state%time%tau(1)**0.5 )  ! time normalization, *tau in estimates.90

    !if(elem%i == itest) &
    !     write(197,'(a8,es12.4,a2,3i5,50es12.4)') 'algeb', elem%eta(resA, 1),'|',dofA,dof1, dof, elem%vec(rhs, 1:dof*ndim)


    !print*, "space discretization error"
    call EvalMaximumLang(dofA, dof1, dofA, S(1:dofA, 1:dofA), &
         elem%vec(rhs, 1:dofA*ndim), elem%eta(resS, 1))

    !print*, "space discretization error with dofS, for regularity"
    call EvalMaximumLang(dofA, dof1, dofS, S(1:dofA, 1:dofA), &
         elem%vec(rhs, 1:dofA*ndim), elem%eta(resSr, 1) )


    !print*,' Verfurth - not normalization,centrally done in SolveProblem1'
    !!elem%eta(resS, 1) = elem%eta(resS, 1) * (state%time%tau(1)**0.5 ) ! time normalization, *tau in estimates.90

    !if( elem%i <= 5 .or. sum( abs(elem%xc(:) ) ) < 0.4) &
    !     write(*,'(a10, 3i5, 2es12.4,a2,30es12.4)') 'estims:gtt', &
    !     elem%i, dof, dofA,elem%xc(:), '|',elem%eta(resA, 1),  elem%eta(resSr, 1) , elem%eta(resS, 1)


    ! A POSTERIORI ERROR ESTIMATES TEST
    !if(.not. onlyAS) then
    !   ! evaluation of the function maximazing the supremum
    !   allocate( xix(1:ndim, 1:dofA))

    !   xix(:,:) = 0.
    !   call EvalMaximumLang(dofA, dof1, dofA, S(1:dofA, 1:dofA), &
    !        elem%vec(rhs, 1:dofA*ndim), elem%eta(resS, 1), xix(1:ndim, dof1:dofA) )

    !   call PlotElemFunction3D(1000+state%space%adapt%adapt_level, elem, dofA, xix(1, 1:dofA) )
    !
    !   deallocate(xix)
    !endif
    ! END test for the dual function


    !if(elem%i == itest) &
    !     write(197,'(a8,es12.4,a2,3i5,50es12.4)') 'space', elem%eta(resS, 1),'|',dofA,dof1, dofA, elem%vec(rhs, 1:dofA*ndim)

    if(.not. onlyAS) then
       ! "time discretization error: q_K+1 component"
       call EvalMaximumLang(dofA, dof1, dof, S(1:dofA, 1:dofA), &
            elem%vec(rhsT, 1:dofA*ndim), elem%eta(resT, 1))

       !if(elem%i == itest) &
       !     write(197,'(a8,es12.4,a2,3i5,50es12.4)') 'time', elem%eta(resT, 1),'|',dofA,dof1, dof, elem%vec(rhsT, 1:dof*ndim)



       ! Verfurth - nor normalization,centrally done in SolveProblem1
       !!elem%eta(resT, 1) = elem%eta(resT, 1) * (state%time%tau(1)**0.5 ) ! time normalization, *tau in estimates.90

       ! STDGM approach
       !elem%eta(resT, 1) = elem%eta(resT, 1) / (state%time%tau(1)**2 + 12 )**0.5
       !elem%eta(resT, 1) = (1.* elem%eta(resT, 1)**2 + elem%eta(resA, 1)**2)**0.5

       ! "space-time discretization error: q_K+1 component"
       call EvalMaximumLang(dofA, dof1, dofA, S(1:dofA, 1:dofA), &
            elem%vec(rhsT, 1:dofA*ndim), elem%eta(resST, 1))

       !if(elem%i == itest) &
       !     write(197,'(a8,es12.4,a2,3i5,50es12.4)') 'ST', elem%eta(resST, 1),'|',dofA,dof1, dofA, elem%vec(rhsT, 1:dofA*ndim)
       !if(elem%i == itest) write(197,'(x)')

       ! Verfurth - nor normalization,centrally done in SolveProblem1
       !elem%eta(resST, 1) = elem%eta(resST, 1) * (state%time%tau(1)**0.5 ) ! time normalization, *tau in estimates.90

       ! STDGM approach
       !elem%eta(resST, 1) = elem%eta(resST, 1) / (state%time%tau(1)**2 + 12 )**0.5
       !elem%eta(resST, 1) = (1.* elem%eta(resST, 1)**2 + elem%eta(resS, 1)**2)**0.5

    endif


    deallocate(S)


  end subroutine DualElemEstimate

  !> DWR estimation of the element reziduum
  !> STATIONARY PROBLEMS ONLY
  !> computes the space estimate from elem%zST,zSTplus,vec(rhs,:)
  !> \f$ \max_{\| \varphi\|_X=1} \f$ frod Rezid method is substituted by the DUAL solution
  !> \f$ \frac{|c_h({\bf w}, \varphi_i) |}{\| zSTplus - zST \|_X}\f$, 
  !> watch if (zSTplus - zST) or zSTplus is used
  function DWRElemEstim_space(elem) result( estim )
    type(element), intent(in) :: elem
    real, dimension(1:ndim) :: estim
    integer :: j, dofA, dof
    real, allocatable, dimension(:) :: zi
!    real, allocatable, dimension(:) :: rhside
!    real, dimension(1:ndim) :: try_est

    dof = elem%dof
    dofA = elem%dof_plus

    !TODO FR CONTROL NDIM>1
     allocate( zi(1:dofA) )

!     allocate( rhside( 1:dofA*ndim ))
!     rhside(1:ndim*dofA) = Transfer_funST_to_fun( elem%rhsST(1:ndim, 1:dofA, 1:elem%Tdof_plus), &
!                           dofA, elem%Tdof_plus, 0, elem%TQnum )

!!     rhside2(1:ndim*dofA) = Transfer_funST_to_fun( elem%rhsST(1:ndim, 1:dofA, 1:elem%Tdof), &
!!                           dofA, elem%Tdof, 0, elem%TQnum )

    ! should it be here only idim not 1:ndim ????
    ! would not work for tdeg > 0 -> we need to transfer zST > z in space format

    if ( ndim > 1) &
      print*, 'NDIM>1 may lead to problems in DWRElemEstim_space'

!    print*, 'zST is not used in DWRElemEstim_Space!'
    do j = 1,ndim
       ! DWR_PARAMETERS
       ! DWR-EST - control whether zp or (zp - zh) is used
      zi(1:dof) = elem%zSTplus( j, 1:dof,1 ) - elem%zST(j,1:dof,1)
      zi(dof+1:dofA) = elem%zSTplus(j,dof+1:dofA,1)

      estim(j) = dot_product( zi(1:dofA), elem%vec( rhs, (j-1)*dofA + 1: (j-1)*dofA + dofA ) )

      ! FR: this is better , WATCH when is elem%vec( rhs, :) filled ? It seems that both is right in this case.
!      try_estim(j) = dot_product( zi(1:dofA), rhside( (j-1)*dofA + 1: (j-1)*dofA + dofA ) )
!      if ( abs(estim(j) - try_est(j)) > 1E-9 )  then
!         print*,  'estim(j) /= try_est(j): ' , estim(j) , try_est(j)
!         stop
!      endif
    enddo

    deallocate(zi)
!    deallocate(rhside)

  end function DWRElemEstim_space

  !> DWR estimation of the element reziduum
  !> STATIONARY PROBLEMS ONLY
  !> computes the ALGEBRAIC estimate from elem%zST vec(rhs,:)
  !> \f$ \max_{\| \varphi\|_X=1} \f$ from Rezid method is substituted by the DUAL solution
  !> \f$ \frac{|c_h({\bf w}, \varphi_i) |}{\| zST \|_X}\f$,
  function DWRElemEstim_alg(elem) result( estim )
    type(element), intent(in) :: elem
    real, dimension(1:ndim)  :: estim
    integer :: j, dofA,dof

    dof = elem%dof
    dofA = elem%dof_plus

    estim(1:ndim) = 0.0

!    print*, 'cdsadas = ' , elem%zST(1,1:dof,1) , elem%vec( rhs, 1:dof )

    do j = 1,ndim
    estim(j) = dot_product( elem%zST(j,1:dof,1) , elem%vec( rhs, (j-1)*dofA + 1: (j-1)*dofA+dof ) )
!      estim = estim + &
!         dot_product( elem%zST(j,1:dof,1) , elem%vec( rhs, (j-1)*dofA + 1: (j-1)*dofA+dof ) )
    enddo

!    print*, 'estim = ' , estim(1)

  end function DWRElemEstim_alg


  !> DWR estimation of the element reziduum
  !> STATIONARY PROBLEMS ONLY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> computes the space estimate from elem%zST,zSTplus,vec(rhs,:)
  !> \f$ \frac{|c_h({\bf w}, \varphi_i) |}{\| zSTplus - zST \|_X}\f$, 
  !> watch if (zSTplus - zST) or zSTplus is used
  !> DWR local error estimate using PARTITION OF UNITY
  function DWRElemEstim_PU(elem) result( pu_estim )
    type(element), intent(in) :: elem
    real, dimension(1:3) :: pu_estim ! estimate for three main vertices of the triangle
    integer :: j, dofPlus, dofPlusTwo, dof, Qnum, Qdof
    real, allocatable, dimension(:) :: zi
    real, dimension(:), allocatable :: ti ! zi in integ nodes
    real, dimension(:,:), allocatable :: wi ! zi*psi_i in basis coeffs
    real, dimension(:,:), allocatable :: hat
    type(volume_rule), pointer :: V_rule

    if ( ndim > 1) &
      print*, 'NDIM>1 may lead to problems in DWRElemEstim_spacePU'

!    print*, 'HOW TO SET V_rule(Qnum) ?', elem%Qnum
    Qnum = elem%Qnum + 1
    V_rule => state%space%V_rule(Qnum)

    if (V_rule%Qdof < DOFtriang(elem%deg+2) ) &
      stop 'problem in DWRElemEstim_PU, not enough quad. nodes!'

    ! number of the quadrature points
    Qdof = V_rule%Qdof

    ! allocate the hat functions 1:number of vertices, 1: number of intgeg points
    allocate( hat(1:3,1:Qdof), source = 0.0 )

    ! integ nodes barycentric coordinates
    ! hat func <---> with the right vertex
    ! hat function in integ nodes
    hat(1:3,1:Qdof) = transpose( V_rule%lambda(1:Qdof, 1:3) )

    dof = elem%dof
    dofPlus = elem%dof_plus

    allocate( zi(1:dofPlus) )
!    print*, 'ZI CHANGED - does not contain zSTplus'
    ! compute ZSTplus - zST
    do j = 1, 1
       ! if ndim > 1 STOP NEEDS CHANGES
       zi(1:dof) =  elem%zSTplus( j, 1:dof,1 ) - elem%zST(j,1:dof,1)
       zi(dof+1:dofPlus) = elem%zSTplus(j,dof+1:dofPlus,1)
    enddo
!    print*, 'zST:' ,  elem%zST(1,1:dof,1)
!    print*, 'diff:',  elem%zSTplus( 1, 1:dof,1 ) - elem%zST(1,1:dof,1)
!    print*, 'zi:' , zi(1:dofPlus)

    ! zSTplus-zST in integ nodes
    allocate( ti(1:Qdof), source = 0.0 )
    call EvalFuncInIntNodes( V_rule, dofPlus, zi(1:dofPlus), ti(1:Qdof) )
    deallocate(zi)

    ! hat = (zSTplus - zST) * hat
    do j = 1,3
      hat(j, 1:Qdof) = hat(j,1:Qdof) * ti(1:Qdof)
    enddo
    deallocate( ti )

    ! t -> express in basis of P^{p+2} !!!
    ! this should be p+2 or even more
    dofPlusTwo = DOFtriang( elem%deg + state%space%plusDeg + 1 )
    allocate( wi(1:3, 1:dofPlusTwo), source = 0.0 )
    call Trans_Integ_Nodes_to_Basis( elem, dofPlusTwo, 3, V_rule, &
                          hat(1:3, 1:Qdof), wi(1:3, 1:dofPlusTwo) )
    deallocate(hat)

    ! now we have in wi the basis coefficients of zSTplus - zST) * hat
    ! get together with the residual

      ! use matmul instead ;-)
      ! only for ndim = 1 !!!
    do j = 1,3
         pu_estim(j) = dot_product( wi(j,1:dofPlusTwo), elem%rhsST( 1,1:dofPlusTwo, 1 ) )
    end do

    deallocate( wi )

  end function DWRElemEstim_PU



  !> evaluate the appropriate maximum over a discrete set
  subroutine EvalMaximumLang(dof_tot, dof_L, dof_U, S, vec, estim, xix)
    integer, intent(in) :: dof_tot, dof_L, dof_U ! dof: total, lowe, upper
    real, dimension(1:dof_tot, 1:dof_tot), intent(in) :: S
    real, dimension(1:dof_tot*ndim), intent(in) :: vec
    real, intent(inout) :: estim
    real, dimension(1:ndim, dof_L:dof_U ), intent(inout), optional :: xix
    real, dimension(:, :), allocatable :: Sinv, alpha
    real :: lambda
    integer :: i, k, k1, k2



    if(dof_L > dof_tot .or. dof_U < dof_L ) then
       print*,'bad arrays bounds in project.f90:EvalMaximumLang',dof_tot, dof_L, dof_U
       stop
    endif

    allocate(Sinv(dof_L:dof_U, dof_L:dof_U) )

    ! "enrichment" of the diagonal following from the Lagrange multiplier
    Sinv(dof_L:dof_U, dof_L:dof_U ) = 2*S(dof_L:dof_U, dof_L:dof_U )
    ! correction of the bug vvv ------^^
    !do i = dof_L, dof_U
    !   Sinv(i,i) = 2*Sinv(i,i)
    !enddo

    call MblockInverse(dof_U - dof_L +1 , Sinv(dof_L:dof_U, dof_L:dof_U ))

    allocate(alpha(1:ndim, dof_L:dof_U) )
    alpha(:,:) = 0.

    estim = 0.
    do k=1,ndim
       k1 = (k-1) * dof_tot + dof_L
       k2 = (k-1) * dof_tot + dof_U

       !write(*,'(a6,12es12.4)') 'alpha:',vec(k1:k2), Sinv(dof_L:dof_U, dof_L:dof_U)

       alpha(k, dof_L:dof_U) = matmul(Sinv(dof_L:dof_U, dof_L:dof_U), vec(k1:k2) )

       ! if(dot_product(alpha(k, dof_L:dof_U), &
       !      matmul(S(dof_L:dof_U,dof_L:dof_U), alpha(k, dof_L:dof_U) )) <= 0.  ) then
       !    write(*,'(a6,120es12.4)') 'alpha:',alpha(k, dof_L:dof_U)
       !    write(*,'(a6,3i5,12es12.4)') 'alpha:',k, dof_L, dof_U, dot_product(alpha(k, dof_L:dof_U), &
       !         matmul(S(dof_L:dof_U,dof_L:dof_U), alpha(k, dof_L:dof_U) ))
       !    do i=1, dof_U
       !       write(*,'(a6,i4,120es12.4)') 'SS:',i, S(i, :)
       !    enddo

       ! endif

       lambda = sqrt(abs( dot_product(alpha(k, dof_L:dof_U), &
            matmul(S(dof_L:dof_U,dof_L:dof_U), alpha(k, dof_L:dof_U) ) ) ) )

       !!!write(197,*) '///',k,k1,k2,dof_L,dof_U, lambda

       if(lambda > 0) then
          alpha(k, dof_L:dof_U) = alpha(k, dof_L:dof_U) / lambda
          estim = estim + abs(dot_product(vec(k1:k2), alpha(k, dof_L:dof_U)) )

       !else
          ! alpha = 0., estim = 0.
       endif
    enddo

    if( present(xix) ) xix(1:ndim, dof_L:dof_U) = alpha(1:ndim, dof_L:dof_U)

    deallocate(Sinv, alpha)


  end subroutine EvalMaximumLang




  !> estimation of the regularity of the solution based on the basis coefficients
  !> expansions
  subroutine ProjectionElemEstimate(elem)
    type(element) :: elem
    real, dimension(:), allocatable :: err_p
    real, dimension(:), allocatable :: coeff_p
    real, dimension(:), allocatable :: b, a, c
    real :: s, s1, s1a, s2
    integer :: lmin, lmax
    integer :: deg, dof, ldeg, i, k

    deg = elem%deg
    dof = elem%dof

    allocate(err_p(-1:deg) )
    allocate(coeff_p(-1:deg) )

    allocate(a(0:deg), b(0:deg), c(0:deg))

    err_p(:) = 0.
    coeff_p(:) = 0.

    err_p(-1) =  dot_product(elem%w(0,1:ndim*dof), elem%w(0,1:ndim*dof) )

    do i=0, deg
       ldeg = state%space%ldeg(i) + 1

       lmin = state%space%ldeg(i-1)+1
       lmax = state%space%ldeg(i)

       !if(elem%i == 1) write(*,'(a3,8i5)') 'DeG',i,state%space%ldeg(i), lmin, lmax, ldeg, dof

       do k=1, ndim
          err_p(i) = err_p(i) + dot_product(elem%w(0,(k-1)*dof + ldeg: k*dof), &
               elem%w(0,(k-1)*dof + ldeg: k*dof) )

          coeff_p(i) = coeff_p(i) &
               + dot_product(elem%w(0,(k-1)*dof + lmin: (k-1)*dof + lmax), &
               elem%w(0,(k-1)*dof + lmin: (k-1)*dof + lmax) )
       enddo
    enddo
    err_p(-1:deg) = err_p(-1:deg)**0.5

    a(0:deg) = coeff_p(0:deg)**0.5

    if(deg > 0) then
       elem%reg = log(a(deg-1)/a(deg))
    else
       elem%reg = 1.
    endif

     if(elem%i <= -3  .or. elem%i == -55 .or. elem%i == -87 ) then

       print*,'   i  ideg    a(i)    sigma     b(i)   c(i)  a(i)/a(i-1) log(err)/log(p) '
       write(*,'(a25, i5, a4,2es12.4, a15)') &
            '------------------- elem=',elem%i,'xc=', elem%xc(:),'  ------------'
       do i=1, deg
          !a(i) = coeff_p(i)

          b(i) = (a(i) / a(i-1) )**2
          b(i) =  i*(1.-b(i))/(1+b(i))


          c(i) = log(1./a(i)/a(i) )/ (2* log(1.*i) )

          write(*,'(2i5, 2es10.2,a2,2es10.2,a2,2es10.2)') &
               i, state%space%ldeg(i), a(i), log(a(i-1)/a(i) ),' |', &
               b(i), c(i), &
               !i*(1.-b(i))/2, &
               ' |', &
               a(i)/a(i-1), &
               log(err_p(i)/err_p(i-1) ) / log(1.*(i-1)/i )


          write(99,*) &
               i, elem%xc(:), a(i), a(0), log(a(i-1)/a(i) ), &
               b(i), c(i), &
               a(i)/a(i-1), &
               log(err_p(i)/err_p(i-1) ) / log(1.*(i-1)/i )


       enddo
       print*,'------------------------------------------------'
    endif




    deallocate(err_p, coeff_p)
    deallocate(a, b, c)

  end subroutine ProjectionElemEstimate

  !> estimation of the regularity of the exact solution based on
  !>\f$ \int_{\partial K} [{\bf w}]^2\, dS \f$,
  !>
  subroutine ElementRegularityEstimator(elem)
    type(element), intent(inout):: elem ! elem = element
    class(element), pointer ::   elem1  ! elem1 = neigh element
    real, dimension(:,:), allocatable :: phi, phi1 ! test functions
    real, dimension(:, :), allocatable :: wi       ! w recomputed  in integ nodes
    real :: val, val1, valM

    integer ::  dof, dof1, ie1,  Qnum, Qdof, dofM, dofM1
    integer :: ie, l, l1, ii, k, kst, kst1

    dof = elem%dof
    dofM = elem%deg*(elem%deg+1)/2

    val  = 0.
    valM = 0.

    do ie =1, elem%flen
       ii = elem%face(neigh, ie)
       if( ii > 0) then  !! inner face
          !! seting of degree of the Gauss quadrature
          Qnum = elem%face(fGnum,ie)
          Qdof = state%space%G_rule(Qnum)%Qdof

          allocate(wi(1:2, 1:Qdof ) )

          allocate(phi(1:dof, 1:Qdof))
          call Eval_Phi_Edge(elem, dof, ie, phi, .false.)

          elem1 => grid%elem(ii)

          dof1 = elem1%dof
          dofM1 = elem1%deg*(elem1%deg+1)/2
          ie1 = elem%face(nei_i,ie)

          allocate(phi1(1:dof1, 1:Qdof))
          call Eval_Phi_Edge(elem1, dof1, ie1, phi1, .true.)


          ! evaluation of (w_ie^+ - w_ie^-) in integ. nodes
          do k=1,ndim
             kst  = (k-1)*dof + 1
             kst1 = (k-1)*dof1 + 1
             do l=1, Qdof       ! ndim ONLY the DENSITY
                wi(1, l) = (dot_product(phi(1:dof ,l),  elem%w(0,kst : kst+dof-1 ) ) &
                     - dot_product(phi1(1:dof1 ,l), elem1%w(0,kst1: kst1+dof1-1) ) )**2
                if(dofM > 0 .and. dofM1 >0) &
                     wi(2, l) = (dot_product(phi(1:dofM ,l), elem%w(0,kst : kst+dofM-1) )&
                     - dot_product(phi1(1:dofM1 ,l), elem1%w(0,kst1: kst1+dofM1-1) ) )**2
             enddo !! l
          enddo  !! k

          val  = val  + dot_product(wi(1,:), state%space%G_rule(Qnum)%weights(:) ) * elem%dn(ie)
          val1 = val1 + dot_product(wi(2,:), state%space%G_rule(Qnum)%weights(:) ) * elem%dn(ie)

          write(200+state%space%adapt%adapt_level,*) elem%xc(:),elem%eta(resST, 1), elem%errL2, &
               elem%reg,val/val1,val,val1,val/val1/(elem%diam**2)

          elem%reg = val/val1/(elem%diam**1)

          deallocate(wi, phi, phi1)


       endif
    enddo

  end subroutine ElementRegularityEstimator

end module project_estimation
