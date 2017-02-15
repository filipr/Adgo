!> setting of basis functions
module basis
  use mesh_oper
  use main_data
!  use model_oper
!  use model3D_oper
!  use eval_sol
  use blocks_integ

  implicit none

!  public:: SetBasisFunctions
  public:: Lagr2Basis
  public:: Lagr2BasisDiff
  public:: Lagr2BasisDiffOne
!  public:: SetTrulePhi
  public:: EvalTrulePhi
!  public:: Write_TrulePhi

  public:: SetEdgeNodes
  public:: SetEdgeNodesHG

  public :: PHI_orthonormal




contains

  !> recompute the polynomial function defined in Lagrangian nodes
  !> to the coefficient of basis functions of the SAME DEGREE as elem
  subroutine Lagr2Basis(elem, Ndeg,  w, Qnum, dof,  wi)
    type(element), intent(inout):: elem        ! elem = element
    integer, intent(in) :: Ndeg   ! degree of Lagrangian
    integer, intent(in) :: dof    ! degree of freedom recomputed
    integer, intent(in) :: Qnum   ! index of the used volume quadrature
    real, dimension(1:ndim, 1:(Ndeg+1)*(Ndeg+2)/2), intent(in) :: w  ! sol in Lag.
    real, dimension(1:ndim, 1:dof), intent(inout) :: wi ! recomp sol
    real, dimension(:), allocatable :: qi, f
    real, dimension(:,:), allocatable :: psi
    real, dimension(:,:), pointer :: xi
    integer :: Ndof,  Qdof, k

    Ndof = (Ndeg+1)*(Ndeg+2)/2

    Qdof = state%space%V_rule(Qnum)%Qdof

    allocate( qi(1:dof) )
    allocate( f(1:Qdof) )
    allocate( psi(1:Ndof, 1:Qdof))

    ! evaluation of the Lagrangian test function in integration nodes
    xi => state%space%V_rule(Qnum)%lambda(1:Qdof,1:3)
    !xi => state%space%L_rule(Qnum)%lambda(1:Qdof,1:3)

    ! older variant
    !call Eval_L_rule(state%space%L_rule(Ndeg), Qdof, xi, psi)

    ! better variant
    call Eval_L_Direct(state%space%L_rule(Ndeg), Qdof, xi, psi)

    do k=1, ndim
       f(1:Qdof) = matmul(w(k,1:Ndof), psi(1:Ndof, 1:Qdof) )

       qi(:) = 0.
       call IntegrateVectorB(elem, dof, f(1:Qdof), qi(1:dof) )

       wi(k,1:dof) = matmul(elem%MassInv%Mb(1:dof,1:dof), qi(1:dof) )

    enddo
    deallocate(qi,f, psi)

    !stop
  end subroutine Lagr2Basis


  !> recompute the polynomial function defined in Lagrangian nodes
  !> to the coefficient of basis functions of the DIFFERENT DEGREE as elem
  subroutine Lagr2BasisDiff(elem, Ndeg,  w, Qnum, dof,  wi)
    type(element), intent(in):: elem        ! elem = element
    integer, intent(in) :: Ndeg   ! degree of Lagrangian
    real, dimension(1:ndim, 1:(Ndeg+1)*(Ndeg+2)/2), intent(in) :: w  ! sol in Lag.
    integer, intent(in) :: Qnum   ! index of the used volume quadrature
    integer, intent(in) :: dof    ! degree of freedom to be recomputed
    real, dimension(1:ndim, 1:dof), intent(inout) :: wi ! recomp sol
    real, dimension(:), allocatable :: qi, f, weights
    real, dimension(:,:), allocatable :: psi, Mass
    real, dimension(:,:), pointer :: xi, phi
    integer :: Ndof,  Qdof, k, l

    Ndof = (Ndeg+1)*(Ndeg+2)/2

    Qdof = state%space%V_rule(Qnum)%Qdof

    allocate( qi(1:dof) )
    allocate( f(1:Qdof) )
    allocate( psi(1:Ndof, 1:Qdof))

    ! evaluation of the Lagrangian test function in integration nodes
    xi => state%space%V_rule(Qnum)%lambda(1:Qdof,1:3)
    !call Eval_L_rule(state%space%L_rule(Ndeg), Qdof, xi, psi)
    ! better variant
    call Eval_L_Direct(state%space%L_rule(Ndeg), Qdof, xi, psi)

    ! basis test function
    phi => state%space%V_rule(Qnum)%phi(1:dof, 1:Qdof)

    allocate( weights(1:Qdof) )
    if(.not. elem%F%iFlin .and. elem%Qnum == Qnum) then !linear element, constant Jacobian
       ! if Qnum /= elem%Qnum, only some approximation is used
       weights(1:Qdof)  = state%space%V_rule(elem%Qnum)%weights(1:elem%Qdof) &
            * elem%F%V%JF(1:Qdof)/ (5 -elem%type)
    else
       weights(1:Qdof)  = state%space%V_rule(Qnum)%weights(1:Qdof) * elem%F%JF0 /(5 -elem%type)
    endif

    ! computation of the mass matrix
    allocate( Mass(1:dof, 1:dof) )
    do k=1,dof
       do l=k,dof
          Mass(k,l) = dot_product(weights(1:Qdof)* phi(k, 1:Qdof), phi(l,1:Qdof) )
          Mass(l,k) = Mass(k,l)
       enddo
    enddo
    ! its inversion
    call MblockInverse(dof, Mass)

    do k=1, ndim
       f(1:Qdof) = matmul(w(k,1:Ndof), psi(1:Ndof, 1:Qdof) )
       f(1:Qdof) = f(1:Qdof) * weights(1:Qdof)
       !call PlotElemFunctionQ(400+state%time%iter, elem, 'V', Qnum, Qdof, f(1:Qdof) )

       qi(:) = 0.
       do l=1, dof
          qi(l) = qi(l) + dot_product( f(1:Qdof), phi(l, 1:Qdof) )
       enddo

       wi(k,1:dof) = matmul(Mass(1:dof,1:dof), qi(1:dof) )
    enddo
    deallocate(qi, f, psi, weights, Mass)

  end subroutine Lagr2BasisDiff


  !> recompute the polynomial function defined in Lagrangian nodes
  !> to the volume integ nodes of the DIFFERENT DEGREE as elem
  subroutine Lagr2QnodesDiff(elem, Ndeg,  w, V_rule,  wi)
    type(element), intent(in):: elem        ! elem = element
    integer, intent(in) :: Ndeg   ! degree of Lagrangian
    real, dimension(1:(Ndeg+1)*(Ndeg+2)/2, 1:ndim), intent(in) :: w  ! sol in Lag.
    type(volume_rule), target, intent(in) :: V_rule
    real, dimension( 1:elem%Qdof, 1:ndim), intent(inout) :: wi ! recomp sol
    real, dimension(:,:), allocatable :: psi 
    real, dimension(:,:), pointer :: xi 
    integer :: Ndof,  Qdof, k

    Ndof = (Ndeg+1)*(Ndeg+2)/2

    Qdof = V_rule%Qdof

    allocate( psi(1:Ndof, 1:Qdof))

    ! evaluation of the Lagrangian test function in integration nodes
    xi => V_rule%lambda(1:Qdof,1:3)

    !call Eval_L_rule(state%space%L_rule(Ndeg), Qdof, xi, psi)
    ! better variant
    call Eval_L_Direct(state%space%L_rule(Ndeg), Qdof, xi, psi)


    do k=1, ndim
       wi(1:Qdof, k) = matmul(w(1:Ndof, k), psi(1:Ndof, 1:Qdof) )
    enddo
    deallocate( psi) !, qi, weights, Mass)

  end subroutine Lagr2QnodesDiff



  !> recompute the polynomial function defined in Lagrangian nodes
  !> to the coefficient of basis functions of the DIFFERENT DEGREE as elem
  !> only one variable
  subroutine Lagr2BasisDiffOne(elem, Ndeg,  w, Qnum, dof,  wi)
    type(element), intent(in):: elem        ! elem = element
    integer, intent(in) :: Ndeg   ! degree of Lagrangian
    integer, intent(in) :: dof    ! degree of freedom to be recomputed
    integer, intent(in) :: Qnum   ! index of the used volume quadrature
    real, dimension(1:(Ndeg+1)*(Ndeg+2)/2), intent(in) :: w  ! sol in Lag.
    real, dimension(1:dof), intent(inout) :: wi ! recomp sol
    real, dimension(:), allocatable :: qi, f, weights
    real, dimension(:,:), allocatable :: psi, Mass
    real, dimension(:,:), pointer :: xi, phi
    integer :: Ndof,  Qdof, k, l

    Ndof = (Ndeg+1)*(Ndeg+2)/2

    Qdof = state%space%V_rule(Qnum)%Qdof

    allocate( qi(1:dof) )
    allocate( f(1:Qdof) )
    allocate( psi(1:Ndof, 1:Qdof))

    ! evaluation of the Lagrangian test function in integration nodes
    xi => state%space%V_rule(Qnum)%lambda(1:Qdof,1:3)

    !call Eval_L_rule(state%space%L_rule(Ndeg), Qdof, xi, psi)
    ! better variant
    call Eval_L_Direct(state%space%L_rule(Ndeg), Qdof, xi(1:Qdof,1:3), psi)

    ! basis test function
    phi => state%space%V_rule(Qnum)%phi(1:dof, 1:Qdof)

    allocate( weights(1:Qdof) )

    if(.not. elem%F%iFlin  .and. elem%Qnum == Qnum) then !linear element, constant Jacobian
       ! if Qnum /= elem%Qnum, only some approximation is used
       weights(1:Qdof)  = state%space%V_rule(elem%Qnum)%weights(1:elem%Qdof) &
            * elem%F%V%JF(1:Qdof)/ (5 -elem%type)
    else
       weights(1:Qdof)  = state%space%V_rule(Qnum)%weights(1:Qdof) * elem%F%JF0 /(5 -elem%type)
    endif

    ! computation of the mass matrix
    allocate( Mass(1:dof, 1:dof) )
    do k=1,dof
       do l=1,dof
          Mass(k,l) = dot_product(weights(1:Qdof)* phi(k, 1:Qdof), phi(l,1:Qdof) )
          Mass(l,k) = Mass(k,l)
       enddo
    enddo
    ! its inversion
    call MblockInverse(dof, Mass)

    !do k=1, ndim

    f(1:Qdof) = matmul(w(1:Ndof), psi(1:Ndof, 1:Qdof) )
    f(1:Qdof) = f(1:Qdof) * weights(1:Qdof)

    !qi(:) = 0.
    !do l=1, dof
    !   qi(l) = qi(l) + dot_product( f(1:Qdof), phi(l, 1:Qdof) )
    !enddo
    qi(1:dof) = matmul(phi(1:dof, 1:Qdof), f(1:Qdof) )

    wi(1:dof) = matmul(Mass(1:dof,1:dof), qi(1:dof) )
    !enddo
    deallocate(qi,f, psi, weights, Mass)

  end subroutine Lagr2BasisDiffOne


  !> NOT USED ANYMORE
  !> evaluation of orthonormal test functions in node \f$ (x, y) \f$
!  subroutine PHI_orthonormal_OLD(Qdof, xi, len, dof,  phi, Dphi)
!    integer, intent(in) :: Qdof          ! number of integ nodes
!    real, dimension (1:Qdof, 1:nbDim), intent(in) :: xi   ! coordinates of integ. node
!    integer, intent(in) :: len           ! triangle/square
!    integer, intent(in) :: dof           ! number of test functions
!    real, dimension(1:dof,1:Qdof) :: phi        ! value of test function
!    real, dimension(1:dof,1:nbDim,1:Qdof), optional:: Dphi  ! value of derivatives of test function
!    real :: x, y
!    integer :: l
!
!    real :: s2, s3, s5, s6, s7, s10, s11, s14, s15, s21, s30, s33, s35, s70
!
!    s2 = 2.**0.5
!    s3 = 3.**0.5
!    s5 = 5.**0.5
!    s7 = 7.**0.5
!    s11 = 11.**0.5
!    s33 = 33.**0.5
!
!    s6  = s2 * s3
!    s10 = s2 * s5
!    s14 = s2 * s7
!    s15 = s3 * s5
!    s21 = s3 * s7
!    s30 = s2 * s3 * s5
!    s35 = s5 * s7
!    s70 = s7 * s10
!
!    do l=1,Qdof
!       x = xi(l,1)
!       y = xi(l,2)
!       if(dof > 21) then
!          print*,'Implementaion in ORTH_PHI only till dof=15'
!          stop
!       else
!          if(len == 3) then ! triangle
!             if(dof>=1) phi( 1, l) = s2
!             if(dof>=2) phi( 2, l) = 6 * (x-1./3)
!             if(dof>=3) phi( 3, l) = 4 * s3 * (y +x/2 -0.5)
!             if(dof>=4) phi( 4, l) = 10. * s6 * (x *(x-0.8) + 0.1)
!             if(dof>=5) phi( 5, l) = 30. * s2 * (x *(x/2+y-0.6) -0.2*y + 0.1)
!             if(dof>=6) phi( 6, l) =      s30 * ( (x - 1)*(x-1) + 6*y*(y+x-1) )
!             if(dof>=7) phi( 7, l) = 10. * s2 * (x *(x *(7.*x-9.) +3.) - 0.2)
!             if(dof>=8) phi( 8, l) =  2. * s6 * (x *(x *(21*x-33)+13) + x*y*(42*x -24) + 2*y -1)
!             if(dof>=9) phi( 9, l) = 2. * s10 * (42*x*y*(x+y-8./7) + 6*y*(1-y) + x*(x*(7*x-15)+9)-1)
!             if(dof>=10) phi(10, l) = 2.*s14 *(x*y*(12*x+30*y-24) + y*(y*(20*y-30)+12) + (x-1)**3)
!             ! added P4 basis
!             if(dof>=11) phi(11, l) = s10*(x*(x*(x*(126*x-224)+126)-24)+1)
!             if(dof>=12) phi(12, l) = s30*(x*(x*(x*(84*x-168)+105)-22)+1+y*(x*(x*(168*x-168)+42)-2))
!             if(dof>=13) phi(13, l) = 5*s2*(x*(x*(x*(36*x-88)+69)-18)+1+y*(x*(x*(216*x-312)+102)-6) &
!                  + y*y*(x*(216*x-96)+6))
!             if(dof>=14) phi(14, l) = s70*(x*(x*(x*(9*x-28)+30)-12)+1+y*(x*(x*(108*x-228)+132)-12) &
!                  + y*y*(x*(270*x-300)+30)+y*y*y*(180*x-20))
!             if(dof>=15) phi(15, l) = 3*s10*(x*(x*(x*(x-4)+6)-4)+1+y*(x*(x*(20*x-60)+60)-20) &
!                  + y*y*(x*(90*x-180)+90)+y*y*y*(140*x-140)+70*y*y*y*y)
!!!!!!!!!!!!!!!!!!
!             ! added P5 basis
!             if(dof>=16) phi(16, l) = 2*s3*(x*(x*(x*(x*(462*x-1050)+840)-280)+35)-1)
!             if(dof>=17) phi(17, l) = 6*(x*(x*(x*(x*(330*x-810)+696)-248)+33)-1 &
!                  + y*(x*(x*(x*(660*x-960)+432)-64)+2))
!             if(dof>=18) phi(18, l) = 2*s15*(x*(x*(x*(x*(165*x-465)+462)-190)+29)-1 &
!                  + y*(x*(x*(x*(990*x-1800)+972)-168)+6) &
!                  + y*y*(x*(x*(990*x-810)+162)-6))
!             if(dof>=19) phi(19, l) = 2*s21*(x*(x*(x*(x*(55*x-185)+226)-118)+23)-1 &
!                  + y*(x*(x*(x*(660*x-1560)+1152)-264)+12) &
!                  + y*y*(x*(x*(1650*x-2250)+630)-30)+y*y*y*(x*(1100*x-400)+20))
!             if(dof>=20) phi(20, l) = 6*s3*(x*(x*(x*(x*(11*x-45)+70)-50)+15)-1 &
!                  + y*(x*(x*(x*(220*x-680)+720)-280)+20) &
!                  + y*y*(x*(x*(990*x-2070)+1170)-90) &
!                  + y*y*y*(x*(1540*x-1680)+140)+y*y*y*y*(770*x-70))
!             if(dof>=21) phi(21, l) = 2*s33*(x*(x*(x*(x*(x-5)+10)-10)+5)-1 &
!                  + y*(x*(x*(x*(30*x-120)+180)-120)+30) &
!                  + y*y*(x*(x*(210*x-630)+630)-210) &
!                  + y*y*y*(x*(560*x-1120)+560) &
!                  + y*y*y*y*(630*x-630)+252*y*y*y*y*y)
!!!!!!!!!!!!!!!!!!
!             if(present(Dphi)) then
!                if(dof>=1) Dphi( 1, 1, l) = 0.
!                if(dof>=1) Dphi( 1, 2, l) = 0.
!                if(dof>=2) Dphi( 2, 1, l) = 6.
!                if(dof>=3) Dphi( 3, 1, l) = 2*s3
!                if(dof>=2) Dphi( 2, 2, l) = 0.
!                if(dof>=3) Dphi( 3, 2, l) = 4*s3
!                if(dof>=4) Dphi( 4, 1, l) = 10*s6 * (2*x-0.8)
!                if(dof>=4) Dphi( 4, 2, l) = 0.
!                if(dof>=5) Dphi( 5, 1, l) = 30.*s2* (x + y -0.6)
!                if(dof>=5) Dphi( 5, 2, l) = 30.*s2* (x - 0.2)
!                if(dof>=6) Dphi( 6, 1, l) = 2 * s30*(x + 3*y -1)
!                if(dof>=6) Dphi( 6, 2, l) = 6 * s30*(x + 2*y -1)
!                if(dof>=7) Dphi( 7, 1, l) = 10 * s2*(x*(21*x-18) +3)
!                if(dof>=7) Dphi( 7, 2, l) = 0.
!                if(dof>=8) Dphi( 8, 1, l) = 2.* s6 *(12*y*(7*x-2) + x*(63*x - 66) + 13)
!                if(dof>=8) Dphi( 8, 2, l) = 4.* s6 *(x*(21*x - 12) + 1)
!                if(dof>=9) Dphi( 9, 1, l) = 6.*s10 *(2*y*(7*y-8) + x*(7*x-10) +28*x*y + 3)
!                if(dof>=9) Dphi( 9, 2, l) = 6.*s10 *(2*x*(7*x-8) + 4*y*(7*x -1) + 2 )
!                if(dof>=10) Dphi(10, 1, l) = 6.*s14 *((x-1)**2 + 8*y*(x-1) + 10*y*y)
!                if(dof>=10) Dphi(10, 2, l) = 24.*s14*(5*y*(y-1) + x*(x+5*y-2) + 1)
!                ! added P4 basis derivatives
!                if(dof>=11) Dphi(11, 1, l) = 6*s10*(x*(x*(84*x-112)+42)-4)
!                if(dof>=11) Dphi(11, 2, l) = 0.
!                if(dof>=12) Dphi(12, 1, l) = 2*s30*(x*(x*(168*x-252)+105)-11+y*(x*(252*x-168)+21))
!                if(dof>=12) Dphi(12, 2, l) = 2*s30*(x*(x*(84*x-84)+21)-1)
!                if(dof>=13) Dphi(13, 1, l) = 30*s2*(x*(x*(24*x-44)+23)-3+y*(x*(108*x-104)+17) &
!                     + y*y*(72*x-16))
!                if(dof>=13) Dphi(13, 2, l) = 30*s2*(x*(x*(36*x-52)+17)-1+y*(x*(72*x-32)+2))
!                if(dof>=14) Dphi(14, 1, l) = 12*s70*(x*(x*(3*x-7)+5)-1+y*(x*(27*x-38)+11) &
!                     + y*y*(45*x-25)+15*y*y*y)
!                if(dof>=14) Dphi(14, 2, l) = 12*s70*(x*(x*(9*x-19)+11)-1+y*(x*(45*x-50)+5) &
!                     + y*y*(45*x-5))
!                if(dof>=15) Dphi(15, 1, l) = 6*s10*(x*(x*(2*x-6)+6)-2+y*(x*(30*x-60)+30) &
!                     + y*y*(90*x-90)+70*y*y*y)
!                if(dof>=15) Dphi(15, 2, l) = 30*s10*(x*(x*(2*x-6)+6)-2+y*(x*(18*x-36)+18) &
!                     + y*y*(42*x-42)+28*y*y*y)
!!!!!!!!!!!!!!!!!!!!!!
!                ! added P5 basis derivatives
!                if(dof>=16) Dphi(16, 1, l) = 14*s3*(x*(x*(x*(330*x-600)+360)-80)+5)
!                if(dof>=16) Dphi(16, 2, l) = 0.
!                if(dof>=17) Dphi(17, 1, l) = 6*(x*(x*(x*(1650*x-3240)+2088)-496)+33 &
!                     + y*(x*(x*(2640*x-2880)+864)-64))
!                if(dof>=17) Dphi(17, 2, l) = 12*(x*(x*(x*(330*x-480)+216)-32)+1)
!                if(dof>=18) Dphi(18, 1, l) = 2*s15*(x*(x*(x*(825*x-1860)+1386)-380)+29 &
!                     + y*(x*(x*(3960*x-5400)+1944)-168) &
!                     + y*y*(x*(2970*x-1620)+162))
!                if(dof>=18) Dphi(18, 2, l) = 12*s15*(x*(x*(x*(165*x-300)+162)-28)+1 &
!                     + y*(x*(x*(330*x-270)+54)-2))
!                if(dof>=19) Dphi(19, 1, l) = 2*s21*(x*(x*(x*(275*x-740)+678)-236)+23 &
!                     + y*(x*(x*(2640*x-4680)+2304)-264) &
!                     + y*y*(x*(4950*x-4500)+630)+y*y*y*(2200*x-400))
!                if(dof>=19) Dphi(19, 2, l) = 8*s21*(x*(x*(x*(165*x-390)+288)-66)+3 &
!                     + y*(x*(x*(825*x-1125)+315)-15)+y*y*(x*(825*x-300)+15))
!                if(dof>=20) Dphi(20, 1, l) = 30*s3*(x*(x*(x*(11*x-36)+42)-20)+3 &
!                     + y*(x*(x*(176*x-408)+288)-56) &
!                     + y*y*(x*(594*x-828)+234)+y*y*y*(616*x-336)+154*y*y*y*y)
!                if(dof>=20) Dphi(20, 2, l) = 60*s3*(x*(x*(x*(22*x-68)+72)-28)+2 &
!                     + y*(x*(x*(198*x-414)+234)-18) &
!                     + y*y*(x*(462*x-504)+42)+y*y*y*(308*x-28))
!                if(dof>=21) Dphi(21, 1, l) = 2*s33*(x*(x*(x*(5*x-20)+30)-20)+5 &
!                     + y*(x*(x*(120*x-360)+360)-120) &
!                     + y*y*(x*(630*x-1260)+630)+y*y*y*(1120*x-1120)+630*y*y*y*y)
!                if(dof>=21) Dphi(21, 2, l) = 12*s33*(x*(x*(x*(5*x-20)+30)-20)+5 &
!                     + y*(x*(x*(70*x-210)+210)-70) &
!                     + y*y*(x*(280*x-560)+280)+y*y*y*(420*x-420)+210*y*y*y*y)
!!!!!!!!!!!!!!!!!!!!!!
!             endif
!          elseif(len == 4) then !square
!             if(dof>=1) phi( 1, l) = 1
!             if(dof>=2) phi( 2, l) = s3 * (2*x - 1)
!             if(dof>=3) phi( 3, l) = s3 * (2*y - 1)
!             if(dof>=4) phi( 4, l) = s5 * (6*x*(x-1) + 1)
!             if(dof>=5) phi( 5, l) = 3. * (2*x*(y-1) + 2*y*(x-1) + 1)
!             if(dof>=6) phi( 6, l) = s5 * (6*y*(y-1) + 1)
!             if(dof>=7) phi( 7, l) = s7 * (x*(x*(20*x - 30) + 12) - 1)
!             if(dof>=8) phi( 8, l) = s15*(2*y-1)*(6*x*(x-1) + 1)
!             if(dof>=9) phi( 9, l) = s15*(2*x-1)*(6*y*(y-1) + 1)
!             if(dof>=10) phi(10, l) = s7 * (y*(y*(20*y - 30) + 12) - 1)
!             ! added P4 basis
!             if(dof>=11) phi(11, l) = 3*(x*(x*(x*(70*x-140)+90)-20)+1)
!             if(dof>=12) phi(12, l) = s21*(x*(x*(-20*x+30)-12)+1+y*(x*(x*(40*x-60)+24)-2))
!             if(dof>=13) phi(13, l) = 5*(x*(6*x-6)+1+y*(x*(-36*x+36)-6)+y*y*(x*(36*x-36)+6))
!             if(dof>=14) phi(14, l) = s21*(y*(y*(-20*y+30)-12)+1+x*(y*(y*(40*y-60)+24)-2))
!             if(dof>=15) phi(15, l) = 3*(y*(y*(y*(70*y-140)+90)-20)+1)
!!!!!!!!!!!!!!!!!!
!             ! added P5 basis
!             if(dof>=16) phi(16, l) = s11*(x*(x*(x*(x*(252*x-630)+560)-210)+30)-1)
!             if(dof>=17) phi(17, l) = 3*s3*(x*(x*(x*(-70*x+140)-90)+20)-1 &
!                  + y*(x*(x*(x*(140*x-280)+180)-40)+2))
!             if(dof>=18) phi(18, l) = s35*(x*(x*(20*x-30)+12)-1 &
!                  + y*(x*(x*(-120*x+180)-72)+6)+y*y*(x*(x*(120*x-180)+72)-6))
!             if(dof>=19) phi(19, l) = s35*(y*(y*(20*y-30)+12)-1 &
!                  + x*(y*(y*(-120*y+180)-72)+6)+x*x*(y*(y*(120*y-180)+72)-6))
!             if(dof>=20) phi(20, l) = 3*s3*(y*(y*(y*(-70*y+140)-90)+20)-1 &
!                  + x*(y*(y*(y*(140*y-280)+180)-40)+2))
!             if(dof>=21) phi(21, l) = s11*(y*(y*(y*(y*(252*y-630)+560)-210)+30)-1)
!!!!!!!!!!!!!!!!!!
!             if(present(Dphi)) then
!                if(dof>=1) Dphi( 1, 1, l) = 0.
!                if(dof>=1) Dphi( 1, 2, l) = 0.
!                if(dof>=2) Dphi( 2, 1, l) = 2.*s3
!                if(dof>=2) Dphi( 2, 2, l) = 0.
!                if(dof>=3) Dphi( 3, 1, l) = 0.
!                if(dof>=3) Dphi( 3, 2, l) = 2.*s3
!                if(dof>=4) Dphi( 4, 1, l) = 6. * s5 * (2*x - 1)
!                if(dof>=4) Dphi( 4, 2, l) = 0.
!                if(dof>=5) Dphi( 5, 1, l) = 6. * (2*y - 1)
!                if(dof>=5) Dphi( 5, 2, l) = 6. * (2*x - 1)
!                if(dof>=6) Dphi( 6, 1, l) = 0
!                if(dof>=6) Dphi( 6, 2, l) = 6. * s5 * (2*y - 1)
!                if(dof>=7) Dphi( 7, 1, l) = 4.  * s7 * (15*x*(x-1) + 3)
!                if(dof>=7) Dphi( 7, 2, l) = 0.
!                if(dof>=8) Dphi( 8, 1, l) = 6. * s15 * (2*x*(y-1) + 2*y*(x-1) + 1)
!                if(dof>=8) Dphi( 8, 2, l) = 2. * s15 * (6*x*(x-1) + 1)
!                if(dof>=9) Dphi( 9, 1, l) = 2. * s15 * (6*y*(y-1) + 1)
!                if(dof>=9) Dphi( 9, 2, l) = 6. * s15 * (2*x*(y-1) + 2*y*(x-1) + 1)
!                if(dof>=10) Dphi(10, 1, l) = 0.
!                if(dof>=10) Dphi(10, 2, l) = 4. * s7  * (15*y*(y-1) + 3)
!                ! added P4 basis derivatives
!                if(dof>=11) Dphi(11, 1, l) = 30*(x*(x*(28*x-42)+18)-2)
!                if(dof>=11) Dphi(11, 2, l) = 0.
!                if(dof>=12) Dphi(12, 1, l) = 4*s21*(x*(-15*x+15)-3+y*(x*(30*x-30)+6))
!                if(dof>=12) Dphi(12, 2, l) = 2*s21*(x*(x*(20*x-30)+12)-1)
!                if(dof>=13) Dphi(13, 1, l) = 30*(2*x-1+y*(-12*x+6)+y*y*(12*x-6))
!                if(dof>=13) Dphi(13, 2, l) = 30*(x*(-6*x+6)-1+y*(x*(12*x-12)+2))
!                if(dof>=14) Dphi(14, 1, l) = 2*s21*(y*(y*(20*y-30)+12)-1)
!                if(dof>=14) Dphi(14, 2, l) = 4*s21*(y*(-15*y+15)-3+x*(y*(30*y-30)+6))
!                if(dof>=15) Dphi(15, 1, l) = 0.
!                if(dof>=15) Dphi(15, 2, l) = 30*(y*(y*(28*y-42)+18)-2)
!!!!!!!!!!!!!!!!!!!!!!
!                ! added P5 basis derivatives
!                if(dof>=16) Dphi(16, 1, l) = 6*s11*(x*(x*(x*(210*x-420)+280)-70)+5)
!                if(dof>=16) Dphi(16, 2, l) = 0.
!                if(dof>=17) Dphi(17, 1, l) = 60*s3*(x*(x*(-14*x+21)-9)+1+y*(x*(x*(28*x-42)+18)-2))
!                if(dof>=17) Dphi(17, 2, l) = 6*s3*(x*(x*(x*(70*x-140)+90)-20)+1)
!                if(dof>=18) Dphi(18, 1, l) = 12*s35*(x*(5*x-5)+1+y*(x*(-30*x+30)-6) &
!                     + y*y*(x*(30*x-30)+6))
!                if(dof>=18) Dphi(18, 2, l) = 6*s35*(x*(x*(-20*x+30)-12)+1+y*(x*(x*(40*x-60)+24)-2))
!                if(dof>=19) Dphi(19, 1, l) = 6*s35*(y*(y*(-20*y+30)-12)+1+x*(y*(y*(40*y-60)+24)-2))
!                if(dof>=19) Dphi(19, 2, l) = 12*s35*(y*(5*y-5)+1+x*(y*(-30*y+30)-6) &
!                     + x*x*(y*(30*y-30)+6))
!                if(dof>=20) Dphi(20, 1, l) = 6*s3*(y*(y*(y*(70*y-140)+90)-20)+1)
!                if(dof>=20) Dphi(20, 2, l) = 60*s3*(y*(y*(-14*y+21)-9)+1+x*(y*(y*(28*y-42)+18)-2))
!                if(dof>=21) Dphi(21, 1, l) = 0.
!                if(dof>=21) Dphi(21, 2, l) = s11*(y*(y*(y*(y*(252*y-630)+560)-210)+30)-1)
!!!!!!!!!!!!!!!!!!!!!!
!             endif
!
!          else
!             print*,'Only triangle or square in ORTH_PHI implemented'
!          endif
!       endif
!    end do
!
!  end subroutine PHI_orthonormal_OLD


  !> evaluation of the coefficients of the transformation to the orthonormal basis
  !> by GramSchmidt
  subroutine ComputeGramSchmidtCoefficients( )
    real, dimension(:,:), allocatable :: phi
    type(volume_rule), pointer :: V_rule
    real, dimension(1:nbDim) :: xi, xc
    integer:: Qnum, Qdof, dof
    integer:: i, k, len, l

    dof = (MaxDegreeImplemented + 1) * (MaxDegreeImplemented + 2) / 2

    allocate(state%space%GScoeff(3:4, 1:dof, 1:dof) )

    do len = 3, 4  ! triangles or quadrilaterals
       if(len == 3) then
          Qnum = maxVrule
          xc(1:nbDim) = 1./3
       else
          Qnum = maxVrule + maxGrule
          xc(1:nbDim) = 1./2
       endif

       ! setting of test functions in integ nodes
       V_rule => state%space%V_rule(Qnum)

       Qdof = V_rule%Qdof
       allocate(phi(1:dof,1:Qdof) )

       do l=1, Qdof ! values in l-th integ. node
          xi(1:nbDim) = V_rule%lambda(l,1:nbDim)
          call phi_values(dof, xi-xc, phi(1:dof, l) )

          !write(90+i,'(25es12.4)')xi(1:nbDim)-xc(1:nbDim), phi(1:nbDim1, l)
       enddo

       ! Gram-Schmidt
       do k=1,dof
          do i=1, k-1
             state%space%GScoeff(len, i, k) = dot_product(V_rule%weights(1:Qdof),  &
                  phi(k, 1:Qdof) * phi(i, 1:Qdof) ) * 0.5
          enddo

          do i=1, k-1
             phi(k, 1:Qdof) = phi(k, 1:Qdof) - phi(i, 1:Qdof) * state%space%GScoeff(len, i, k)
          enddo

          state%space%GScoeff(len, k, k) = dot_product(V_rule%weights(1:Qdof),  &
               phi(k, 1:Qdof) * phi(k, 1:Qdof) ) * 0.5

          if(state%space%GScoeff(len, k, k) <= 0) then
             print*, 'Trouble in  ComputeGramSchmidtCoefficients',dof, Qnum, Qdof
             stop
          endif

          state%space%GScoeff(len, k, k) = (state%space%GScoeff(len, k, k))**0.5

          phi(k, 1:Qdof) =  phi(k, 1:Qdof) / state%space%GScoeff(len, k, k)
       enddo

       !do l=1, Qdof ! values in l-th integ. node
       !   xi(1:nbDim) = V_rule%lambda(l,1:nbDim)
       !   write(90+len,'(25es12.4)') xi(1:nbDim), phi(1:nbDim1, l)
       !enddo

       deallocate(phi)

    enddo
  end subroutine ComputeGramSchmidtCoefficients


  !> doubled by phi_vals routine in o_basis.f90, both routines are used !!!
  !> evalaution of test functions \f$ \varphi_i,\ i=1,\dots, dof \f$ in
  !> the node \f$ xi(1:nbDim) \f$
  subroutine phi_values(dof, xi, phi, Dphi)
    integer, intent(in) :: dof    ! number of test functions
    real, dimension(1:nbDim), intent(in)  :: xi   ! coodrinates of an integ node
    real, dimension(1:dof), intent(inout) :: phi        ! value of test function
    real, dimension(1:dof,1:nbDim), optional :: Dphi   ! value of derivatives of test function
    real :: Fx, Fy, Dfx, Dfy
    integer :: i, k, l

    l = 0
    do k = 0, MaxDegreeImplemented
       do i= 0, k
          l = l + 1
          if(l <= dof) then
             Fx = xi(1)**(k - i)
             Fy = xi(2)**(i)
             phi(l) = Fx * Fy

             if(present(Dphi)) then
                if(k - i == 0) then
                   Dfx = 0
                else
                   DFx = (k -i) * xi(1)**(k - i - 1)
                endif
                if( i == 0) then
                   DFy = 0.
                else
                   DFy = i * xi(2)**(i-1)
                endif

                Dphi(l, 1) = DFx *  Fy
                Dphi(l, 2) =  Fx * DFy
             endif

          endif
       enddo
    enddo
  end subroutine phi_values

  !> doubled in PHI_orthonormalNew in o_basis.90, both routines are used !!!
  !> evaluation of orthonormal test functions in node \f$ (x, y) \f$
  subroutine PHI_orthonormal(Qdof, nbDim, xi, len, dof, phi, Dphi)
    integer, intent(in) :: Qdof           ! number of integ nodes
    integer, intent(in) :: nbDim
    real, dimension(1:Qdof,1:nbDim),intent(in) :: xi ! coordinates of integ. node
    integer, intent(in) :: len           ! triangle/square
    integer, intent(in) :: dof           ! number of test functions
    real, dimension(1:dof,1:Qdof) :: phi        ! value of test function
    real, dimension(1:dof,1:nbDim, 1:Qdof), optional :: Dphi   ! derivatives of test function
    real, dimension(1:nbDim) :: xc
    integer :: i, k, l

    if(dof > (MaxDegreeImplemented+1)*(MaxDegreeImplemented+2)/2 ) then
       print*,'Implementaion in ORTH_PHI only forl deg= MaxDegreeImplemented'
       stop
    else
       if(len == 3) then ! triangle
          xc(1:nbDim) = 1./3
       elseif(len == 4) then !square
          xc(1:nbDim) = 1./2
       else
          print*,'Only triangle or square in ORTH_PHI implemented'
          stop
       endif

       if(present(Dphi)) then
          do l=1,Qdof
             call phi_values(dof, xi(l, 1:nbDim) - xc(1:nbDim), phi(1:dof, l), &
                  Dphi(1:dof, 1:nbDim, l) )
          enddo
       else
          do l=1,Qdof
             call phi_values(dof, xi(l, 1:nbDim) - xc(1:nbDim), phi(1:dof, l) )
          enddo
       endif

       ! Gram-Schmidt
       do k=1,dof
          do i=1, k-1
             phi(k, 1:Qdof) = phi(k, 1:Qdof) - phi(i, 1:Qdof) * state%space%GScoeff(len, i, k)
          enddo
          phi(k, 1:Qdof) =  phi(k, 1:Qdof) / state%space%GScoeff(len, k, k)
       enddo

       if(present(Dphi)) then
          do k=1,dof
             do i=1, k-1
                Dphi(k, 1:nbDim, 1:Qdof) = Dphi(k, 1:nbDim, 1:Qdof) &
                     - Dphi(i, 1:nbDim, 1:Qdof) * state%space%GScoeff(len, i, k)
             enddo
             Dphi(k, 1:nbDim, 1:Qdof) =  Dphi(k, 1:nbDim, 1:Qdof) / state%space%GScoeff(len, k, k)
          enddo
       endif
    end if

  end subroutine PHI_orthonormal

  !> evaluation of "Taylor" test functions in node \f$ (x, y) \f$
  subroutine PHI_Taylor(x, y, len, dof, phi, Dphi)
    real, intent(in) :: x, y             ! coordinates of integ. node
    integer, intent(in) :: len           ! triangle/square
    integer, intent(in) :: dof           ! number of test functions
    real, dimension(1:dof) :: phi        ! value of test function
    real, dimension(1:dof,1:nbDim), optional :: Dphi   ! value of derivatives of test function
    real :: xc, yc

    if(dof > 10) then
       print*,'Implementaion in ORTH_PHI only till dof=10'
       stop
    else
       if(len == 3) then ! triangle
          xc = 1./3
          yc = 1./3
       elseif(len == 4) then !square
          xc = 1./2
          yc = 1./2
       else
          print*,'Only triangle or square in ORTH_PHI implemented'
          stop
       endif

       if(dof>=1) phi(1) = 1.
       if(dof>=2) phi(2) = (x-xc)
       if(dof>=3) phi(3) = (y-yc)
       if(dof>=4) phi(4) = (x-xc)*(x-xc)
       if(dof>=5) phi(5) = (x-xc)*(y-yc)
       if(dof>=6) phi(6) = (y-yc)*(y-yc)
       if(dof>=7) phi(7) = (x-xc)*(x-xc)*(x-xc)
       if(dof>=8) phi(8) = (x-xc)*(x-xc)*(y-yc)
       if(dof>=9) phi(9) = (x-xc)*(y-yc)*(y-yc)
       if(dof>=10) phi(10) = (y-yc)*(y-yc)*(y-yc)

       if(present(Dphi)) then
          if(dof>=1) Dphi(1, 1) = 0.
          if(dof>=2) Dphi(2, 1) = 1.
          if(dof>=3) Dphi(3, 1) = 0.
          if(dof>=4) Dphi(4, 1) = 2*(x-xc)
          if(dof>=5) Dphi(5, 1) = (y-yc)
          if(dof>=6) Dphi(6, 1) = 0.
          if(dof>=7) Dphi(7, 1) = 3.*(x-xc)*(x-xc)
          if(dof>=8) Dphi(8, 1) = 2.*(x-xc)*(y-yc)
          if(dof>=9) Dphi(9, 1) = (y-yc)*(y-yc)
          if(dof>=10) Dphi(10, 1) = 0.

          if(dof>=1) Dphi(1, 2) = 0.
          if(dof>=2) Dphi(2, 2) = 0.
          if(dof>=3) Dphi(3, 2) = 1.
          if(dof>=4) Dphi(4, 2) = 0.
          if(dof>=5) Dphi(5, 2) = (x-xc)
          if(dof>=6) Dphi(6, 2) = 2.*(y-yc)
          if(dof>=7) Dphi(7, 2) = 0.
          if(dof>=8) Dphi(8, 2) = (x-xc)*(x-xc)
          if(dof>=9) Dphi(9, 2) = 2.*(x-xc)*(y-yc)
          if(dof>=10) Dphi(10, 2) = 3.*(y-yc)*(y-yc)
       endif
    endif
  end subroutine PHI_TAYLOR


  !> NOT USED ANYMORE - replaced by ComputeBasisFunctions in o_space.f90
!  !> setting of basis functions, evaluation of test functions and their derivatives
!  !> on reference elements \f$\hat{K} \f$ in all integ node (volumes, edges)
!  !>
!  !> \f$ \hat{\varphi}_{d}(\hat{x}_{iq, l}),\
!  !> \hat{\partial}_j \hat{\varphi}_{d}(\hat{x}_{iq, l} )\f$,
!  !> \f$\hat{\varphi}_{d} \in \f$
!  !> hierachical basis, \f$ \hat{x}_{iq, l}\f$: \f$l^{\rm th}\f$  node
!  !> of \f$iq^{\rm th}\f$ quadrature rule,
!  !> volume rules for triangles and quadrilateralls, edge rules.
!  !> For hanging nodes (HG) each edge can be split in 2, 4, 8, ... equidistant pieces
!  subroutine  SetBasisFunctions(BasisType)
!    interface
!    subroutine BasisType(Qdof, nbDim, xi, len, dof, phi, Dphi)
!       integer, intent(in) :: Qdof           ! number of integ nodes
!       integer, intent(in) :: nbDim
!       real, dimension(1:Qdof,1:nbDim),intent(in) :: xi ! coordinates of integ. node
!       integer, intent(in) :: len           ! triangle/square
!       integer, intent(in) :: dof           ! number of test functions
!       real, dimension(1:dof,1:Qdof) :: phi        ! value of test function
!       real, dimension(1:dof,1:nbDim, 1:Qdof), optional :: Dphi   ! derivatives of test function
!      end subroutine BasisType
!    end interface
!    type(volume_rule), pointer :: V_rule
!    integer :: iq, iK, l, ie, len, i, j, is
!    integer :: max_deg, dof, Qdeg, Qdof
!    real ::  t
!    real, dimension (1:nbDim) :: x0, a0
!    real, dimension(:,:), allocatable :: p0
!    real, dimension(:,:,:), allocatable :: void
!    logical :: V_rule_alloc
!
!    max_deg = MaxDegreeImplemented
!    dof = DOFtriang( max_deg )
!
!    print*, '############ SetBasisFunctions'
!
!    do iq=1, maxVrule+maxGrule  !!! ERROR ???
!       if( iq <= maxVrule) then
!          ! triangular volume quadratures, only some of them are used
!          len = 3
!          V_rule_alloc = .false.
!          do j=0, max_deg
!             if(state%space%Qdeg(j, 1) == iq) V_rule_alloc = .true.
!          enddo
!
!          if(.not. V_rule_alloc) goto 100
!
!       elseif(iq >= maxVrule + 1) then
!          ! quadrilateral  volume quadratures
!          len = 4
!       else
!          ! NO quadratures
!          goto 100
!       endif
!
!       V_rule => state%space%V_rule(iq)
!       Qdeg = V_rule%Qdeg
!       Qdof = V_rule%Qdof
!
!       allocate(V_rule%phi(1:dof,1:Qdof) )
!       allocate(V_rule%Dphi(1:dof,1:nbDim,1:Qdof) )
!
!       allocate(p0(1:Qdof, 1:nbDim) )
!       p0(1:Qdof, 1:nbDim) = V_rule%lambda(1:Qdof,1:nbDim)!values in integ nodes
!
!       call BasisType(Qdof, nbDim, p0, len, dof, V_rule%phi(1:dof, 1:Qdof), &
!            V_rule%Dphi(1:dof,1:nbDim,1:Qdof) )
!
!       ! print*,'plotting of the values of test functions in integ nodes'
!       !if(iq == 15) then
!       !   do i=1,dof
!       !      do j=1,Qdof
!       !         write(20+i, *) p0(j, 1:2), V_rule%phi(i, j)
!       !      enddo
!       !   enddo
!       !endif
!
!       deallocate(p0)
!
!100    continue
!    enddo ! iq
!
!    ! EDGE QUADRATURES
!    do iq=1, maxGrule
!       Qdeg = state%space%G_rule(iq)%Qdeg
!       Qdof = state%space%G_rule(iq)%Qdof
!
!       allocate(state%space%G_rule(iq)%phi(3:4,1:4,state%space%adapt%HG, 1:dof,1:Qdof) )
!       allocate(state%space%G_rule(iq)%Dphi(3:4,1:4,state%space%adapt%HG, 1:dof, 1:nbDim,1:Qdof) )
!
!       do iK = 3,4   ! triangle or quadrilateral
!          do ie = 1, iK   ! loop through edges
!
!             if(iK == 3) then   ! edges of triangle
!                if(ie == 1 ) then
!                   x0(1:nbDim) = (/0., 0./)
!                   a0(1:nbDim) = (/1., 0./)
!                elseif(ie == 2 ) then
!                   x0(1:nbDim) = (/1., 0./)
!                   a0(1:nbDim) = (/-1., 1./)
!                elseif(ie == 3) then
!                   x0(1:nbDim) = (/0., 1./)
!                   a0(1:nbDim) = (/0., -1./)
!                endif
!             else               ! edges of quadrilateral
!                if(ie == 1 ) then
!                   x0(1:nbDim) = (/0., 0./)
!                   a0(1:nbDim) = (/1., 0./)
!                elseif(ie == 2 ) then
!                   x0(1:nbDim) = (/1., 0./)
!                   a0(1:nbDim) = (/0., 1./)
!                elseif(ie == 3) then
!                   x0(1:nbDim) = (/1., 1./)
!                   a0(1:nbDim) = (/-1., 0./)
!                elseif(ie == 4) then
!                   x0(1:nbDim) = (/0., 1./)
!                   a0(1:nbDim) = (/0., -1./)
!                endif
!             endif
!
!             is = 0
!             do i=0,state%space%adapt%max_HGlevel
!                do j=1, 2**i
!                   is = is + 1
!
!                   allocate(p0(1:Qdof, 1:nbDim) )
!                   do l=1, Qdof  !loop though 1D Gauss quadrature rule
!                      t = state%space%G_rule(iq)%lambda(l)
!
!                      p0(l, 1:nbDim) = x0(1:nbDim) +  a0(1:nbDim) * (j-1 +t)/(2**i)
!                   enddo ! l
!
!                   call BasisType(Qdof, nbDim, p0, iK, dof, &
!                        state%space%G_rule(iq)%phi( iK, ie, is, 1:dof, 1:Qdof), &
!                        state%space%G_rule(iq)%Dphi(iK, ie, is, 1:dof, 1:nbDim, 1:Qdof) )
!
!                   deallocate(p0)
!                enddo
!             enddo
!
!          enddo ! ie
!       enddo  ! iK
!    enddo ! iq
!
!    !do l=0, 3
!    !   write(*,'(10es12.4)')  state%space%G_rule(2)%phi(3,1,l+4,2:3,1)
!    !enddo
!    !write(*,'(a4,10e10.3)') '??>.', state%space%G_rule(2)%phi(3,2,1:10,2)
!
!
!    ! LANGRANGE  QUADRATURES
!    do iq=0, maxLrule  + maxVrule
!       if( iq <= maxLrule) then
!          ! triangular volume quadratures
!          len = 3
!       elseif(iq >= maxVrule) then
!          ! quadrilateral  volume quadratures
!          len = 4
!       else
!          ! NO quadratures
!          goto 200
!       endif
!
!       Qdeg = state%space%L_rule(iq)%Qdeg
!       Qdof = state%space%L_rule(iq)%Qdof
!
!       allocate(state%space%L_rule(iq)%phi(1:dof,1:Qdof) )
!       allocate(void(1:dof,1:nbDim,1:Qdof) )
!
!       allocate(p0(1:Qdof, 1:nbDim) )
!
!       p0(1:Qdof, 1:nbDim) =  state%space%L_rule(iq)%lambda(1:Qdof,1:nbDim)  ! values in integ node
!
!       call BasisType(Qdof, nbDim, p0, len, dof, state%space%L_rule(iq)%phi(1:dof, 1:Qdof), &
!            void(1:dof,1:nbDim,1:Qdof) )
!
!       deallocate(p0)
!       deallocate(void)
!
!200    continue
!
!
!    enddo ! iq
!
!!    associate ( time => state%time )
!!    !select type ( time )
!!    !  class is ( TimeTDG_t )
!!
!!         !print*,'begin in SDERFTRE', maxTrule
!!       ! integ rules for time intervals, connected with the GAUSS at this moment
!!       do iq=1, maxTrule
!!          print*,'*********************',iq, time%T_rule(iq)%Qdeg
!!          call SetTrulePhi( time%T_rule(iq))
!!          !do i=1,state%time%T_rule(iq)%Qdof
!!          !   write(100+iq, *) state%time%T_rule(iq)%lambda(i),state%time%T_rule(iq)%phi(:,i)
!!          !   write(200+iq, *) state%time%T_rule(iq)%lambda(i),state%time%T_rule(iq)%Dphi(:,i)
!!          !enddo
!!
!!          !if(iq == 5) stop
!!       enddo
!!       !print*,'stopped in SDERFTRE'
!!       !stop
!!
!!       !writing all state%time%T_rule(:)%Phi to file TRule_Phi
!!       !call Write_TrulePhi()
!!
!!
!!    !end select
!!    end associate
!
!
!  end subroutine SetBasisFunctions

 subroutine EvalTrulePhi(t, dof, Phi, DPhi)
   real, intent (in) :: t
   integer, intent (in) :: dof
   real, dimension(1:dof), intent(inout) :: Phi
   real, dimension(1:dof), intent(inout), optional :: DPhi
   real :: val, dval
   integer :: i


   !dof = size(Phi(:))

  !print*, 'dof in EvalTrulePhi' ,dof

!  if (dof < MaxTimeDegree + 2) then
!  print*, 'not counting all posible PHI(i)', dof
!  endif

!  dof = MaxTimeDegree + 2 ! 2nd +1 for Tdof_plus
! 	deallocate(T_rule%extraNode)
! 	allocate( T_rule%extraNode( 1 : dof ))

  do i=1,dof
     select case (i)
     case(1)
        val = 1.
        dval = 0.
     case(2)
        val = 12.**0.5 *(t - 0.5)
        dval = 12.**0.5
     case(3)
        val = 3.*20.**0.5 *(t*t - t  + 1./6.)
        dval = 3.*20.**0.5 * (2.*t - 1.)
     case(4)
        val = 7.**0.5 * (20. * t**(3.) - 30. * t**(2.) + 12. * t - 1. )
        dval = 7.**0.5 * (60. * t*t - 60. * t + 12.)
     case(5)
        val = 3. * ( 70. * t**(4.) - 140. * t**(3.) + 90. * t**(2.) - 20. * t + 1. )
        dval = 60.*( 14. * t**(3.) - 21.  * t**(2.) + 9.  * t - 1. )
     case(6)
        val = 11.**0.5 * ( 252. * t**(5.) - 630. * t**(4.) + 560. * t**(3.) - 210. * t**(2.) + 30. * t - 1. )
        dval = 11.**0.5 *(1260. * t**(4.) - 2520. *t**(3.) + 1680. *t**(2.) - 420. * t + 30. )
     case(7:)
        print*,' order >= 7 in SetTrulePhi in basis.f90 not implemented'
        stop
     end select

     Phi(i) = val
     if( present(DPhi) ) DPhi(i) = dval

  end do


 end subroutine EvalTrulePhi


  !> return the barycentric (cartesian coordinates on reference element)
  !> coordinates of edge integ nodes
  subroutine  SetEdgeNodes(iq, ie, iK, x )
    integer, intent(in) :: iq    ! Qnum of the Gauss quadrature
    integer, intent(in) :: ie    ! index of the edge
    integer, intent(in) :: iK  !  3 = triangle, 4 = quadrilateral
    real, dimension(1:state%space%G_rule(iq)%Qdof, 1:nbDim), intent(out) :: x
    real, dimension(1:nbDim) :: x0, a0
    integer:: Qdeg, Qdof, l
    real :: t

    Qdeg = state%space%G_rule(iq)%Qdeg
    Qdof = state%space%G_rule(iq)%Qdof

    if(iK == 3) then   ! edges of triangle
       if(ie == 1 ) then
          x0(1:nbDim) = (/0., 0./)
          a0(1:nbDim) = (/1., 0./)
       elseif(ie == 2 ) then
          x0(1:nbDim) = (/1., 0./)
          a0(1:nbDim) = (/-1., 1./)
       elseif(ie == 3) then
          x0(1:nbDim) = (/0., 1./)
          a0(1:nbDim) = (/0., -1./)
       endif
    else               ! edges of quadrilateral
       if(ie == 1 ) then
          x0(1:nbDim) = (/0., 0./)
          a0(1:nbDim) = (/1., 0./)
       elseif(ie == 2 ) then
          x0(1:nbDim) = (/1., 0./)
          a0(1:nbDim) = (/0., 1./)
       elseif(ie == 3) then
          x0(1:nbDim) = (/1., 1./)
          a0(1:nbDim) = (/-1., 0./)
       elseif(ie == 4) then
          x0(1:nbDim) = (/0., 1./)
          a0(1:nbDim) = (/0., -1./)
       endif
    endif

    do l=1, Qdof  !loop though 1D Gauss quadrature rule
       t = state%space%G_rule(iq)%lambda(l)

       !x(l, 1:nbDim) = x0(1:nbDim) +  a0(1:nbDim) * (j-1 +t)/(2**i)  ! for HG NODES
       x(l, 1:nbDim) = x0(1:nbDim) +  a0(1:nbDim) * t
    enddo ! l


  end subroutine SetEdgeNodes


  !> return the barycentric (cartesian coordinates on reference element)
  !> coordinates of edge integ nodes including edges with hanging nodes
  subroutine  SetEdgeNodesHG (HGidx, iq, ie, iK, x)
    integer, intent(in) :: iq    ! Qnum of the Gauss quadrature
    integer, intent(in) :: ie    ! index of the edge
    integer, intent(in) :: HGidx   ! type of the subedge
    integer, intent(in) :: iK  !  3 = triangle, 4 = quadrilateral
    real, dimension(1:state%space%G_rule(iq)%Qdof, 1:nbDim), intent(out) :: x
    real, dimension(1:nbDim) :: x0, a0
    integer:: Qdof, l
    real :: t

    Qdof = state%space%G_rule(iq)%Qdof

    if(iK == 3) then   ! edges of triangle
       if(ie == 1 ) then
          x0(1:nbDim) = (/0., 0./)
          a0(1:nbDim) = (/1., 0./)
       elseif(ie == 2 ) then
          x0(1:nbDim) = (/1., 0./)
          a0(1:nbDim) = (/-1., 1./)
       elseif(ie == 3) then
          x0(1:nbDim) = (/0., 1./)
          a0(1:nbDim) = (/0., -1./)
       endif
    else               ! edges of quadrilateral
       if(ie == 1 ) then
          x0(1:nbDim) = (/0., 0./)
          a0(1:nbDim) = (/1., 0./)
       elseif(ie == 2 ) then
          x0(1:nbDim) = (/1., 0./)
          a0(1:nbDim) = (/0., 1./)
       elseif(ie == 3) then
          x0(1:nbDim) = (/1., 1./)
          a0(1:nbDim) = (/-1., 0./)
       elseif(ie == 4) then
          x0(1:nbDim) = (/0., 1./)
          a0(1:nbDim) = (/0., -1./)
       endif
    endif


    !do i=0,state%space%adapt%max_HGlevel
       !do j=1, 2**i

          do l=1, Qdof  !loop though 1D Gauss quadrature rule
             t = state%space%G_rule(iq)%lambda(l)

             x(l, 1:nbDim) = x0(1:nbDim) +  a0(1:nbDim) * ResizeHG(t, HGidx) !(j-1 +t)/(2**i)
          enddo ! l

       !enddo ! j
    !enddo ! i

  end subroutine SetEdgeNodesHG

end module basis
