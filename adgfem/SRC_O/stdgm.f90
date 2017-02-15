module stdgm_mod
   use blocks_integ
   use mesh_mod

implicit none

   public :: Eval_whST_Elem
   public :: Eval_whST_iCoord_Elem
   public :: evalWSTinTimeDof

   public :: evalL2STNorm_Elem
   public :: evalL2STNormQ_Elem
   public :: evalL2STNormVecQ_Elem
   public :: evalH1L2STNorm_Elem
   public :: evalL8H1STNorm_Elem

   public :: evalSTfunInIntTime
   public :: evalSTfunInIntTime_Der
   public :: evalSTfunInIntTime_spaceDer
   public :: evalSTfunInIntTime_gradient
   public :: evalSTfunInIntTimeDof
   public :: evalSTfunInIntTimeDof_Der
   public :: evalSTfunInRealTime
   public :: evalSTfunInRealTimeDof
   public :: evalSTfunInRealTimeDofAndDt ! w and dw/dt
   public :: evalSTfunInRealTimeLagrange
   public :: evalSTfunInRealTimeLagrange_dxdt ! w, dw/dt, grad(w)
   public :: evalwSTjumpDof
   public :: evalwSTjumpInLagrangeNodes


   public :: Eval_wSTfin_Elem
   public :: Transfer_wST_to_w_Elem
   public :: Transfer_wST_to_w_Elem_Real
   public :: Transfer_wST_to_w_Patch
   public :: Transfer_funST_to_fun
   public :: Transfer_funST_to_funS

contains

! timeNode = 0 -> endpoint of the interval
  !timeNode = -1 -> starting point of the interval
  !> STDG: wST -> w in node timeNode of time quadrature, Tdeg = number of quadrature nodes
  subroutine Transfer_wST_to_w_Elem( elem , timeNode, Tdeg)
    type(element), intent (inout) :: elem
    integer, intent(in) :: timeNode
    integer, intent(in) :: Tdeg
    class(Time_rule), pointer :: T_rule
    integer :: k,dof, Tdof

    dof = elem%dof
    Tdof = elem%Tdof ! number of time basis functions on element elem

    elem%w(0,1:ndim*dof) = 0.

    T_rule => state%time%T_rule(Tdeg)
       !  print*, 'alpha:' , timeNode
       !  print*, 'elem wST in Transfer', elem%wST(1,1:dof,1:Tdof)
       !print*,'####',elem%i, dof, Tdof, timeNode
       do k = 1, ndim
          !write(*,'(120es10.2)') elem%wST(k,1:dof,1:Tdof)
          !write(*,'(120es10.2)') T_rule%phi(1:Tdof,timeNode)

          elem%w(0,dof*(k-1) + 1 :dof*k) = &
               matmul(elem%wST(k,1:dof,1:Tdof), T_rule%phi(1:Tdof,timeNode))
       enddo

      !  print*, 'elem w in Transfer', elem%w( 0,1 :dof )

  end subroutine Transfer_wST_to_w_Elem


  ! timeNode = 0 -> endpoint of the interval
  !timeNode = -1 -> starting point of the interval
  !> STDG: funST format -> fun space in node timeNode of time quadrature, Tdeg = number of quadrature nodes
  function Transfer_funST_to_fun( funST, dof, Tdof, timeNode, Tdeg) result( fun )
    real, dimension(1:ndim,1:dof,1:Tdof), intent(in) :: funST
    integer, intent(in) :: dof
    integer, intent(in) :: Tdof
    integer, intent(in) :: timeNode
    integer, intent(in) :: Tdeg
    real, dimension(1:dof*ndim):: fun

    class(Time_rule), pointer :: T_rule
    integer :: k, big_dof

    big_dof = dof * ndim

    fun(1:big_dof) = 0.0

    T_rule => state%time%T_rule(Tdeg)
       do k = 1, ndim
          fun(dof*(k-1) + 1 :dof*k) = &
               matmul( funST(k,1:dof,1:Tdof), T_rule%phi(1:Tdof,timeNode) )
       enddo

  end function Transfer_funST_to_fun

  ! timeNode = 0 -> endpoint of the interval
  !timeNode = -1 -> starting point of the interval
  !> STDG: funST format -> fun space in node timeNode of time quadrature, Tdeg = number of quadrature nodes
  !> used for elem%wS format
  function Transfer_funST_to_funS( funST, dof, Tdof, timeNode, Tdeg) result( fun )
    real, dimension(1:ndim,1:dof,1:Tdof), intent(in) :: funST
    integer, intent(in) :: dof
    integer, intent(in) :: Tdof
    integer, intent(in) :: timeNode
    integer, intent(in) :: Tdeg
    real, dimension(1:ndim, 1:dof) :: fun

    class(Time_rule), pointer :: T_rule
    integer :: k

    fun(1:ndim,1:dof) = 0.0

    T_rule => state%time%T_rule(Tdeg)
       do k = 1, ndim
          fun(k, 1:dof ) = &
               matmul( funST(k,1:dof,1:Tdof), T_rule%phi(1:Tdof,timeNode) )
       enddo

  end function Transfer_funST_to_funS

   !  transfer the solution elem%wST -> elem%w(0,:) in any time of the reference interval
   subroutine Transfer_wST_to_w_Elem_Real( elem, timeNode )
      class(element), intent(in) :: elem
      real, intent(in) :: timeNode ! time in the reference interval [0,1]

      real, allocatable, dimension(:) :: phi
      integer :: k, dof, Tdof


      dof = elem%dof
      Tdof = elem%Tdof ! number of time basis functions on element elem
      allocate( phi(1:Tdof), source=0.0 )

      ! evaluate the basis time functions in time timeNode
      call Eval_LegendrePolynomialsOneNode( Tdof, timeNode, phi(1:Tdof) )

       do k = 1, ndim
          elem%w(0,dof*(k-1) + 1 :dof*k) = &
               matmul( elem%wST(k,1:dof,1:Tdof), phi(1:Tdof) )
       enddo

      deallocate( phi )

   end subroutine Transfer_wST_to_w_Elem_Real

  ! timeNode = 0 -> endpoint of the interval
  !timeNode = -1 -> starting point of the interval
  !> STDG: transfers wST -> w on elem and its neighbors
  !> in node timeNode of time quadrature, tQnum = number of quadrature nodes
  subroutine Transfer_wST_to_w_Patch(grid, elem, timeNode, tQnum)
   class(mesh), intent(inout) :: grid
   class(element), intent (inout) :: elem
    integer, intent(in) :: timeNode
    integer, intent(in) :: tQnum

    integer :: iedge, ielem

    call Transfer_wST_to_w_Elem( elem, timeNode, tQnum)

    do iedge = 1, elem%flen
      ielem = elem%face(neigh, iedge)

      if ( ielem > 0 ) &
         call Transfer_wST_to_w_Elem( grid%elem( ielem ), timeNode, tQnum)

    end do !iedge

  end subroutine Transfer_wST_to_w_Patch






  !> STDG: evaluate \f$ w_{m}^{-}\f$ on element elem in space integration nodes
  subroutine Eval_wSTfin_Elem( elem )
    type(element), intent (inout) :: elem
    class(Time_rule), pointer :: T_rule
    type(volume_rule), pointer  :: V_rule
    integer :: Qdof, Tdof,dof, i, j, k, l
    real :: temp

    T_rule => state%time%T_rule(elem%TQnum)
    V_rule => state%space%V_rule(elem%Qnum)

    Qdof = elem%Qdof
    Tdof = elem%Tdof
    dof = elem%dof

!    elem%wSTfin(:,:) = 0.
    elem%wSTfinAD(:,:) = 0.

    do i = 1,ndim
       !new version for adaptation
       elem%wSTfinAD(i,1:dof) = matmul(elem%wST(i,1:dof,1:Tdof) , T_rule%phi(1:Tdof,0))

       do j =1,Qdof
          !new
!          elem%wSTfin(i,j) = dot_product(V_rule%phi(1:dof,j), elem%wSTfinAD(i,1:dof) )
          !old
          !elem%wSTfin(i,j) = dot_product(V_rule%phi(1:dof,j), &
          !     matmul(elem%wST(i,1:dof,1:Tdof) , T_rule%phi(1:Tdof,0)) )

       enddo

    enddo
    ! print*, 'elem%wSTfin(1,1:Qdof)', elem%wSTfin(1,1:Qdof)


    !the last parameter is not used in this case - for every Tdeg in elem%phi(:,0) is the final time of the interval
    !call Transfer_wST_to_w_Elem(elem, 0, 1)


  end subroutine Eval_wSTfin_Elem

 ! timeNode = 0 -> endpoint of the interval
  !timeNode = -1 -> starting point of the interval
 !> STDG: evaluate \f$ wh_(t)\f$ on element elem in space integration nodes in the node node of Time quadrature with tQnum nodes
  subroutine Eval_whST_Elem( elem, tQnum, node, wi, Dwi )
    type(element), intent (inout) :: elem
    integer, intent (in) :: tQnum
    integer, intent (in) :: node
    real, dimension( 1:elem%Qdof,1:ndim), intent (inout) :: wi
    real, dimension( 1:elem%Qdof,1:ndim,1:nbDim), intent (inout), optional :: Dwi
    class(Time_rule), pointer :: T_rule
    type(volume_rule), pointer  :: V_rule
    integer :: Qdof, Tdof,dof, i, j, k, l
    real, dimension( 1:elem%dof) :: temp
    real, dimension(1:elem%dof,1:nbDim,1:elem%Qdof) :: Der

    T_rule => state%time%T_rule( tQnum )
    V_rule => state%space%V_rule( elem%Qnum )

    if (present(Dwi)) &
      call Eval_Dphi(elem, elem%dof,  Der)
    Qdof = elem%Qdof
    Tdof = elem%Tdof
    dof = elem%dof

    wi(:,:) = 0.

    if (present(Dwi)) &
      Dwi(:,:,:) = 0.

    do i = 1,ndim
       do j =1,Qdof
          temp(1:dof) = matmul( elem%wST(i,1:dof,1:Tdof), T_rule%phi(1:Tdof,node) )
          wi(j,i) = dot_product(V_rule%phi(1:dof,j), temp(1:dof) )
          if (present(Dwi)) then
            do k = 1, nbdim
               Dwi(j,i,k) = dot_product(Der(1:dof,k,j), temp( 1:dof) )
            enddo
          endif
       enddo
    enddo


  end subroutine Eval_whST_Elem


  ! timeNode = 0 -> endpoint of the interval
  !timeNode = -1 -> starting point of the interval
  !> STDG: evaluate the i-th coordinate (1:ndim) \f$ wh_(t)\f$ on element elem
  !>in space integration nodes in the node node of Time quadrature with tQnum nodes
  function Eval_whST_iCoord_Elem( elem, tQnum, node, iCoord) result( wi )
    type(element), intent (in) :: elem
    integer, intent (in) :: tQnum
    integer, intent (in) :: node
    integer, intent(in) :: iCoord
    real, dimension( 1:elem%Qdof) :: wi
    class(Time_rule), pointer :: T_rule
    type(volume_rule), pointer  :: V_rule
    integer :: Qdof, Tdof,dof, j, k
    real, dimension( 1:elem%dof) :: temp
    real, dimension(1:elem%dof,1:nbDim,1:elem%Qdof) :: Der

    T_rule => state%time%T_rule( tQnum )
    V_rule => state%space%V_rule( elem%Qnum )


    Qdof = elem%Qdof
    Tdof = elem%Tdof
    dof = elem%dof

    wi(:) = 0.

    do j =1,Qdof
       temp(1:dof) = matmul( elem%wST(iCoord,1:dof,1:Tdof), T_rule%phi(1:Tdof,node) )
       wi(j) = dot_product(V_rule%phi(1:dof,j), temp(1:dof) )
    enddo

  end function Eval_whST_iCoord_Elem

   ! eval space-time function in any real time node from [0,1] in Volume integration nodes
  function evalSTfunInRealTime( elem, d, dof, Tdof, fun, node ) result(wi)
      class(element), intent(in) :: elem
      integer, intent(in) :: d, dof, Tdof
      real, dimension(1:d, 1:dof, 1:Tdof), intent(in) :: fun
      real, intent(in) :: node ! time in the reference interval [0,1]
      real, dimension(1:d, 1:elem%Qdof) :: wi

      integer :: i, Qnum, Qdof
      real, dimension(1:d, 1:dof) :: temp
      type(volume_rule), pointer :: V_rule

      temp(1:d, 1:dof) = evalSTfunInRealTimeDof( elem, d, dof, Tdof, fun(1:d,1:dof,1:Tdof), node )

      Qnum = elem%Qnum
      V_rule => state%space%V_rule( Qnum )
      Qdof = V_rule%Qdof

      do i = 1, d
         wi(i, 1:Qdof ) = matmul( temp(i,1:dof) , V_rule%phi( 1:dof,1:Qdof) )
      end do

  end function evalSTfunInRealTime

   ! eval space-time function in INode of a T_rule in Volume integration nodes
  function evalSTfunInIntTime( elem, d, dof, Tdof, fun, iNode, tQnum ) result(wi)
      class(element), intent(in) :: elem
      integer, intent(in) :: d, dof, Tdof
      real, dimension(1:d, 1:dof, 1:Tdof), intent(in) :: fun
      integer, intent(in) :: iNode ! time in the reference interval [0,1]
      integer, intent(in) :: tQnum

      real, dimension(1:d, 1:elem%Qdof) :: wi
      integer :: i, Qnum, Qdof
      real, dimension(1:d, 1:dof) :: temp
      type(volume_rule), pointer :: V_rule


      temp(1:d, 1:dof) = evalSTfunInIntTimeDof( elem, d, dof, Tdof, fun(1:d,1:dof,1:Tdof), iNode, tQnum )

      Qnum = elem%Qnum
      V_rule => state%space%V_rule( Qnum )
      Qdof = V_rule%Qdof

      do i = 1, d
         wi(i, 1:Qdof ) = matmul( temp(i,1:dof) , V_rule%phi( 1:dof,1:Qdof) )
      end do

  end function evalSTfunInIntTime

  ! eval space-time function time derivative in INode of a T_rule in Volume integration nodes
  function evalSTfunInIntTime_Der( elem, tau, d, dof, Tdof, fun, iNode, tQnum ) result(wi)
      class(element), intent(in) :: elem
      real, intent(in) :: tau
      integer, intent(in) :: d, dof, Tdof
      real, dimension(1:d, 1:dof, 1:Tdof), intent(in) :: fun
      integer, intent(in) :: iNode ! time in the reference interval [0,1]
      integer, intent(in) :: tQnum

      real, dimension(1:d, 1:elem%Qdof) :: wi
      integer :: i, Qnum, Qdof
      real, dimension(1:d, 1:dof) :: temp
      type(volume_rule), pointer :: V_rule


      temp(1:d, 1:dof) = evalSTfunInIntTimeDof_Der( elem, tau, d, dof, Tdof, &
         fun(1:d,1:dof,1:Tdof), iNode, tQnum )

      Qnum = elem%Qnum
      V_rule => state%space%V_rule( Qnum )
      Qdof = V_rule%Qdof

      do i = 1, d
         wi(i, 1:Qdof ) = matmul( temp(i,1:dof) , V_rule%phi( 1:dof,1:Qdof) )
      end do

  end function evalSTfunInIntTime_der

   ! eval space-time function space derivative in iNode of a T_rule in Volume integration nodes
  function evalSTfunInIntTime_spaceDer( elem, d, dof, Tdof, fun, iNode, tQnum, dx ) result(wi)
      class(element), intent(in) :: elem
      integer, intent(in) :: d, dof, Tdof
      real, dimension(1:d, 1:dof, 1:Tdof), intent(in) :: fun
      integer, intent(in) :: iNode ! time in the reference interval [0,1]
      integer, intent(in) :: tQnum
      integer, intent(in) :: dx ! derivative with respect to dx-th coordinate

      real, dimension(1:d, 1:elem%Qdof) :: wi
      integer :: i, Qnum, Qdof
      real, dimension(1:d, 1:dof) :: temp
      type(volume_rule), pointer :: V_rule

      real, allocatable, dimension(:,:,:) :: Der


      ! eliminate time
      temp(1:d, 1:dof) = evalSTfunInIntTimeDof( elem, d, dof, Tdof, fun(1:d,1:dof,1:Tdof), iNode, tQnum )

      Qnum = elem%Qnum
      V_rule => state%space%V_rule( Qnum )
      Qdof = V_rule%Qdof

      ! eval Dphi
      allocate( Der( 1:dof, 1:nbDim, 1:Qdof ) )
      call Eval_Dphi_plus( elem, V_rule, dof,  Der )

      do i = 1, d
         wi(i, 1:Qdof) = matmul( temp(i,1:dof) , Der(1:dof, dx, 1:Qdof) )
      end do

      deallocate(Der)

  end function evalSTfunInIntTime_spaceDer

   ! eval space-time function space derivative in iNode of a T_rule in Volume integration nodes
  function evalSTfunInIntTime_gradient( elem, d, dof, Tdof, fun, iNode, tQnum ) result(wi)
      class(element), intent(in) :: elem
      integer, intent(in) :: d, dof, Tdof
      real, dimension(1:d, 1:dof, 1:Tdof), intent(in) :: fun
      integer, intent(in) :: iNode ! time in the reference interval [0,1]
      integer, intent(in) :: tQnum

      real, dimension(1:d, 1:elem%Qdof, 1:nbDim) :: wi ! gradient in integration nodes
      integer :: i, Qnum, Qdof, j
      real, dimension(1:d, 1:dof) :: temp
      type(volume_rule), pointer :: V_rule

      real, allocatable, dimension(:,:,:) :: Der

      ! eliminate time
      temp(1:d, 1:dof) = evalSTfunInIntTimeDof( elem, d, dof, Tdof, fun(1:d,1:dof,1:Tdof), iNode, tQnum )

      Qnum = elem%Qnum
      V_rule => state%space%V_rule( Qnum )
      Qdof = V_rule%Qdof

      ! eval Dphi
      allocate( Der( 1:dof, 1:nbDim, 1:Qdof ) )
      call Eval_Dphi_plus( elem, V_rule, dof,  Der )

      do i = 1, d
         do j = 1,nbDim
         wi(i, 1:Qdof,j) = matmul( temp(i,1:dof) , Der(1:dof, j, 1:Qdof) )
         enddo
      end do

      deallocate(Der)

  end function evalSTfunInIntTime_gradient


   ! eval wST in any time of the reference interval
   function evalWSTinTimeDof( elem, node) result( wi )
      class(element), intent(in) :: elem
      real, intent(in) :: node ! time in the reference interval [0,1]

      real, dimension(1:ndim, 1:elem%dof) :: wi
      real, dimension(1:elem%Tdof) :: phi
      !class( Time_rule ), pointer :: T_rule
      integer :: k

      !T_rule => state%time%T_rule( elem%Tdeg + 1 )


      call Eval_LegendrePolynomialsOneNode( elem%Tdof, node, phi(1:elem%Tdof) )

      do k = 1, ndim
         wi(k, 1:elem%dof) = matmul( elem%wST(k, 1:elem%dof, 1:elem%Tdof) , phi(1:elem%Tdof) )
      end do

   end function evalWSTinTimeDof

   ! eval ST function in any time of the reference interval
   function evalSTfunInRealTimeDof( elem, d, dof, Tdof, fun, node ) result( wi )
      class(element), intent(in) :: elem
      integer, intent(in) :: d, dof, Tdof
      real, dimension(1:d, 1:dof, 1:Tdof), intent(in) :: fun
      real, intent(in) :: node ! time in the reference interval [0,1]

      real, dimension(1:d, 1:dof) :: wi
      real, dimension(1:Tdof) :: phi
      !class( Time_rule ), pointer :: T_rule
      integer :: k

      !T_rule => state%time%T_rule( elem%Tdeg + 1 )

      call Eval_LegendrePolynomialsOneNode( Tdof, node, phi(1:Tdof) )

      do k = 1, d
         wi(k, 1:dof) = matmul( fun(k, 1:dof, 1:Tdof) , phi(1:Tdof) )
      end do

   end function evalStfunInRealTimeDof

   ! eval ST function and its time derivative in any time of the reference interval [0,1]
   subroutine evalSTfunInRealTimeDofAndDt( elem, tau, d, dof, Tdof, fun, node, w, wDt )
      class(element), intent(in) :: elem
      real, intent(in) :: tau
      integer, intent(in) :: d, dof, Tdof
      real, dimension(1:d, 1:dof, 1:Tdof), intent(in) :: fun
      real, intent(in) :: node ! time in the reference interval [0,1]

      real, dimension(1:d, 1:dof), intent(out) :: w, wDt
      real, dimension(1:Tdof) :: phi, dPhi !
      integer :: k

      call Eval_LegendrePolynomialsOneNode( Tdof, node, phi(1:Tdof), &
            dPhi(1:Tdof) )

      do k = 1, d
         w(k, 1:dof) = matmul( fun(k, 1:dof, 1:Tdof) , phi(1:Tdof) )
         wDt(k, 1:dof) = (1.0 / tau) * &
                         matmul( fun(k, 1:dof, 1:Tdof) , dPhi(1:Tdof) )
      end do

   end subroutine evalSTfunInRealTimeDofAndDt

   ! eval ST function in any time of the reference interval
   function evalSTfunInIntTimeDof( elem, d, dof, Tdof, fun, iNode, tQnum ) result( wi )
      class(element), intent(in) :: elem
      integer, intent(in) :: d, dof, Tdof
      real, dimension(1:d, 1:dof, 1:Tdof), intent(in) :: fun
      integer, intent(in) :: iNode ! i_th node of T_rule(tQnum)
      integer, intent(in) :: tQnum
      real, dimension(1:d, 1:dof) :: wi
!      real, dimension(1:Tdof) :: phi
      class( Time_rule ), pointer :: T_rule
      integer :: k

      T_rule => state%time%T_rule( tQnum )

!      call Eval_LegendrePolynomialsOneNode( Tdof, node, phi(1:Tdof) )

      do k = 1, d
         wi(k, 1:dof) = matmul( fun( k, 1:dof, 1:Tdof ) , T_rule%phi( 1:Tdof , iNode ) )
      end do

   end function evalSTfunInIntTimeDof

   ! eval ST function derivative in any INTEGER time of the reference interval
   function evalSTfunInIntTimeDof_Der( elem, tau, d, dof, Tdof, fun, iNode, tQnum ) result( dwi )
      class(element), intent(in) :: elem
      real, intent(in) :: tau
      integer, intent(in) :: d, dof, Tdof
      real, dimension(1:d, 1:dof, 1:Tdof), intent(in) :: fun
      integer, intent(in) :: iNode ! i_th node of T_rule(tQnum)
      integer, intent(in) :: tQnum


      real, dimension(1:d, 1:dof) :: dwi
!      real, dimension(1:Tdof) :: phi
      class( Time_rule ), pointer :: T_rule
      integer :: k

      T_rule => state%time%T_rule( tQnum )

!      call Eval_LegendrePolynomialsOneNode( Tdof, node, phi(1:Tdof) )

      do k = 1, d
         dwi(k, 1:dof) = (1.0 / tau) * matmul( fun( k, 1:dof, 1:Tdof ) , T_rule%dphi( 1:Tdof , iNode ) )
      end do

   end function evalSTfunInIntTimeDof_Der


   !> computes space-time L2 norm SQUARED:
   !> \f$ norm(i) = ( \int_{I_m} \int_K {\bf f_I}^2_{h\tau} )\f$
   !> \f$ {\bf f}_{h\tau}\f$ is piecewise polynomial in space as well as time
   function evalL2STNorm_Elem(elem, tau, d, dof, Tdof, fun ) result( norm )
      class(element), intent(in) :: elem
      real, intent(in) :: tau
      integer, intent(in) :: d ! dimension
      integer, intent(in) :: dof
      integer, intent(in) :: Tdof
      real, dimension(1:d,1:dof,1:Tdof) :: fun ! function coeffs with resp to space and time basis functions
      real, dimension(1:d) :: norm

      integer :: tQnum, Qdof
      integer :: m, k
      class(Time_rule), pointer :: T_rule
      real, dimension(1:d,1:elem%Qdof) :: f_space
      real, dimension(1:d,1:Tdof) :: fTime


      ! degree of f is <= Tdof-1

      tQnum = Tdof !T_rule%Qdeg
      if ( tQnum /= state%time%T_rule( Tdof )%Qdeg) &
         stop 'Problem in T_rule in evalL2STnorm_Elem.'

      Qdof = elem%Qdof
!      ndof = elem%dof*ndim

      norm(:) = 0.

      ! integration over the time interval
      do m=1, tQnum
         f_space(1:d,1:Qdof) = ( evalSTfunInIntTime( elem, d, dof, Tdof, &
            fun(1:d, 1:dof, 1:Tdof), m, tQnum ) )**2.0

         do k = 1,d
            call IntegrateFunction(elem, f_space(k,1:Qdof), fTime(k, m) )
         end do

      end do ! m

      do k=1,d
         norm(k) =  IntegrateTimeFunction( tQnum, tau, fTime(k, 1:tQnum) )
      end do


 end function evalL2STNorm_Elem

   !> computes space-time L2 norm SQUARED:
   !> \f$ norm(i) = ( \int_{I_m} \int_K {\bf f_I}^2_{h\tau} )\f$
   !> \f$ {\bf f}_{h\tau}\f$ is piecewise polynomial in space as well as time given in integration nodes
   function evalL2STNormVecQ_Elem(elem, tau, tQnum, fun ) result( norm )
      class(element), intent(in) :: elem
      real, intent(in) :: tau
      integer, intent(in) :: tQnum
      real, dimension(1:ndim,1:nbDim,1:elem%Qdof,1:tQnum), intent(in) :: fun ! fun in integ nodes
      real, dimension(1:ndim) :: norm

      integer :: Qdof
      integer :: m, k

      real, dimension(1:ndim, 1:elem%Qdof, 1:tQnum) :: ff
      real, dimension(1:ndim, 1:tQnum) :: fTime

      Qdof = elem%Qdof
!      ndof = elem%dof*ndim

      norm(:) = 0.

      ! FOR nbDIM = 2 only
      ff(1:ndim,1:Qdof,1:tQnum) = fun(1:ndim, 1, 1:Qdof, 1:tQnum )**2.0 &
         + fun(1:ndim, 2, 1:Qdof, 1:tQnum )**2.0

      ! integration over the time interval
      do k = 1,ndim
         do m=1, tQnum
            call IntegrateFunction(elem, ff(k, 1:Qdof, m), fTime(k, m) )
         end do !m

         norm(k) =  IntegrateTimeFunction( tQnum, tau, fTime(k, 1:tQnum) )
      end do ! k

   end function evalL2STNormVecQ_Elem


   !> computes space-time L2 norm SQUARED:
   !> \f$ norm(i) = ( \int_{I_m} \int_K {\bf f_I}^2_{h\tau} )\f$
   !> \f$ {\bf f}_{h\tau}\f$ is piecewise polynomial in space and int time, given in integration nodes
   function evalL2STNormQ_Elem(elem, tau, tQnum, fun ) result( norm )
      class(element), intent(in) :: elem
      real, intent(in) :: tau
      integer, intent(in) :: tQnum
      real, dimension(1:ndim,1:elem%Qdof,1:tQnum), intent(in) :: fun ! fun in integ nodes
      real, dimension(1:ndim) :: norm

      integer :: Qdof
      integer :: m, k

      ! real, dimension(1:ndim, 1:elem%Qdof, 1:tQnum) :: ff
      real, dimension(1:ndim, 1:tQnum) :: fTime

      Qdof = elem%Qdof
      norm(:) = 0.

      !fun(1:ndim,1:Qdof,1:tQnum) = fun(1:ndim, 1:Qdof, 1:tQnum )**2.0

      ! integration over the time interval
      do k = 1,ndim
         do m=1, tQnum
            call IntegrateFunction(elem, fun(k, 1:Qdof, m)**2.0, fTime(k, m) )
         end do !m
         norm(k) =  IntegrateTimeFunction( tQnum, tau, fTime(k, 1:tQnum) )
      end do ! k

   end function evalL2STNormQ_Elem





   !> computes space-time L2 norm of the time derivative SQUARED:
   !> \f$ norm(i) = ( \int_{I_m} \int_K {\bf (f_i)'}^2_{h\tau} )\f$
   !> \f$ {\bf f}_{h\tau}\f$ is piecewise polynomial in space as well as time
   function evalH1L2STNorm_Elem(elem, tau, d, dof, Tdof, fun ) result( norm )
      class(element), intent(in) :: elem
      real, intent(in) :: tau
      integer, intent(in) :: d ! dimension
      integer, intent(in) :: dof
      integer, intent(in) :: Tdof
      real, dimension(1:d,1:dof,1:Tdof) :: fun
      real, dimension(1:d) :: norm

      integer :: tQnum, Qdof
      integer :: m, k
      class(Time_rule), pointer :: T_rule
      real, dimension(1:d,1:elem%Qdof) :: f_space
      real, dimension(1:d,1:Tdof) :: fTime


      ! degree of f is <= Tdof-1

      tQnum = Tdof !T_rule%Qdeg
      if ( tQnum /= state%time%T_rule( Tdof )%Qdeg) &
         stop 'Problem in T_rule in evalH1L2STnorm_Elem.'

      Qdof = elem%Qdof
!      ndof = elem%dof*ndim

      norm(:) = 0.

      ! integration over the time interval
      do m=1, tQnum
         f_space(1:d,1:Qdof) = ( evalSTfunInIntTime_Der( elem, tau, d, dof, Tdof, &
            fun(1:d, 1:dof, 1:Tdof), m, tQnum ) )**2.0

         do k = 1,d
            call IntegrateFunction(elem, f_space(k,1:Qdof), fTime(k, m) )
         end do

      end do ! m

      do k=1,d
         norm(k) =  IntegrateTimeFunction( tQnum, tau, fTime(k, 1:tQnum) )
      end do


 end function evalH1L2STNorm_Elem


 !> computes space-time L infinity H1 seminorm norm  SQUARED:
   !> \f$ norm(i) = ( \max_{I_m} \int_K {\bf (f_i)}^2_{h\tau} )\f$
   !> \f$ {\bf f}_{h\tau}\f$ is piecewise polynomial in space as well as time
   function evalL8H1STNorm_Elem(elem, d, dof, Tdof, fun ) result( norm )
      class(element), intent(in) :: elem
      integer, intent(in) :: d ! dimension
      integer, intent(in) :: dof
      integer, intent(in) :: Tdof
      real, dimension(1:d,1:dof,1:Tdof) :: fun
      real, dimension(1:d) :: norm

      integer :: tQnum, Qdof
      integer :: m, k
      class(Time_rule), pointer :: T_rule
      real, dimension(1:d,1:elem%Qdof, 1:nbDim) :: f_space
      real, dimension(1:d) :: norm_loc
      ! degree of f is <= Tdof-1
      norm(1:d) = -1.0

      tQnum = Tdof !T_rule%Qdeg
      if ( tQnum /= state%time%T_rule( Tdof )%Qdeg) &
         stop 'Problem in T_rule in evalH1L2STnorm_Elem.'

      Qdof = elem%Qdof
      norm(:) = 0.

      ! maximum is searched through the time interval
      do m=1, tQnum
         f_space(1:d,1:Qdof, 1:nbDim) = evalSTfunInIntTime_gradient( elem, &
                           ndim, dof, Tdof, fun(1:d,1:dof,1:Tdof), m, tQnum)

         do k = 1,d
           norm_loc(k) =  EvalVectorL2ScalarProduct( elem, nbDim, &
                  f_space(k,1:Qdof,1:nbDim), f_space(k,1:Qdof,1:nbDim) )
         end do
         norm(1:d) = max( norm(1:d), norm_loc(1:d) )

      end do ! m

 end function evalL8H1STNorm_Elem

   !> eval space-time function in any real time node from [0,1] in Lagrange integration nodes
  function evalSTfunInRealTimeLagrange( elem, d, dof, Tdof, fun, node, lagrDeg ) result(wi)
      class(element), intent(in) :: elem
      integer, intent(in) :: d, dof, Tdof
      real, dimension(1:d, 1:dof, 1:Tdof), intent(in) :: fun
      real, intent(in) :: node ! time in the reference interval [0,1]
      integer, intent(in) :: lagrDeg ! degree of the Lagrangian rule
      real, dimension(1:d, 1:((lagrDeg+1)*(lagrDeg+2)/2)) :: wi
      integer :: i, Qnum, Qdof
      real, dimension(1:d, 1:dof) :: temp
      type(Lagrang_rule), pointer :: L_rule

      temp(1:d, 1:dof) = evalSTfunInRealTimeDof( elem, d, dof, Tdof, fun(1:d,1:dof,1:Tdof), node )

      L_rule => state%space%L_rule( lagrDeg )
      Qdof = L_rule%Qdof

      do i = 1, d
         wi(i, 1:Qdof ) = matmul( temp(i,1:dof) , L_rule%phi( 1:dof,1:Qdof) )
      end do

  end function evalSTfunInRealTimeLagrange

   !> eval space-time function in any real time node from [0,1] in Lagrange integration nodes
  subroutine evalSTfunInRealTimeLagrange_dxdt( elem, d, dof, Tdof, fun, tNode, &
                                               lagrDeg , Qdof, w , wDt , wDx)
      class(element), intent(in) :: elem
      integer, intent(in) :: d, dof, Tdof
      real, dimension(1:d, 1:dof, 1:Tdof), intent(in) :: fun
      real, intent(in) :: tNode ! time in the reference interval [0,1]
      integer, intent(in) :: lagrDeg ! degree of the Lagrangian rule
      integer, intent(in) :: Qdof ! number of nodes
      real, dimension(1:d, 1:Qdof), intent(out) :: w
      real, dimension(1:d, 1:Qdof), intent(out) :: wDt
      real, dimension(1:d, 1:Qdof, 1:nbDim), intent(out) :: wDx

      integer :: i, j, Qnum
      real, dimension(1:d, 1:dof) :: wDof, wDtDof
      real, dimension(:,:,:), allocatable :: dPhiSpace
      type(Lagrang_rule), pointer :: L_rule
      real :: tau

      tau = state%time%tau(1)
      ! compute w and dw/dt in tNode
      call evalSTfunInRealTimeDofAndDt( elem, tau, d, dof, Tdof, &
                                        fun(1:d,1:dof,1:Tdof), &
                                        tNode, wDof(1:d, 1:dof), &
                                        wDtDof(1:d, 1:dof) )

      L_rule => state%space%L_rule( lagrDeg )
      if( Qdof /= L_rule%Qdof) &
         stop 'dimension problem in evalSTfunInRealTimeLagrange_dxdt'

      allocate( dPhiSpace(1:dof,1:nbDim, 1:Qdof ), source = 0.0 )
      ! compute dPhi in space - these cannot be taken directly from L_rule ( derivative of F )
      call Eval_Dphi_L_rule( lagrDeg, elem, dof, dPhiSpace )

      do i = 1, d
         w(i, 1:Qdof ) = matmul( wDof(i,1:dof) , L_rule%phi( 1:dof,1:Qdof) )
         wDt(i, 1:Qdof ) = matmul( wDtDof(i,1:dof) , L_rule%phi( 1:dof,1:Qdof) )
         do j = 1, nbDim
            wDx(i,1:Qdof, j) = matmul( wDof(i,1:dof) , dPhiSpace(1:dof, j, 1:Qdof) )
         end do
      end do

      deallocate( dPhiSpace )

  end subroutine evalSTfunInRealTimeLagrange_dxdt

   !> computes the jump in the solution
   !> we have to control whether in elem%wSTfinAD is still the solution in t_m-1 not t_m !
   function evalwSTjumpDof( elem ) result ( w_jump )
      class( element ), intent(in) :: elem
      real, dimension(1:ndim, 1:elem%dof) :: w_jump

      integer :: i, dof, Tdof
      class( Time_rule ), pointer :: T_rule

      dof = elem%dof
      Tdof = elem%Tdof
      T_rule => state%time%T_rule( Tdof )

      do i = 1, ndim
         ! eval the solution in the time moment t_{m-1}^+, phi(i,-1) = phi_i(0)
         w_jump(i, 1:dof) = matmul( elem%wST( i, 1:elem%dof, 1:elem%Tdof) , T_rule%phi(1:elem%Tdof,-1) )
      end do

      !eval the jump $\{ u_h(t_{m-1})\}$
      w_jump(1:ndim, 1:dof) = w_jump(1:ndim, 1:dof) - elem%wSTfinAD( 1:ndim, 1:dof )

   end function evalwSTjumpDof

   !> computes the jump in the solution in given set of nodes - given by lambda coords
   !> we have to control whether in elem%wSTfinAD is still the solution in t_m-1 not t_m !
   function evalwSTjumpInLagrangeNodes( elem, lagrDeg ) result ( w_jump )
      class( element ), intent(in) :: elem
      integer, intent(in) :: lagrDeg
      real, dimension(1:ndim, 1: (lagrDeg+1)*(lagrDeg+2)/2 ) :: w_jump
      integer :: i, dof, Qdof
      real, dimension(1:ndim, 1:elem%dof) :: w_jumpDof
      type(Lagrang_rule), pointer :: L_rule

      dof = elem%dof
      w_jumpDof = evalwSTjumpDof( elem )

      L_rule => state%space%L_rule( lagrDeg )
      Qdof = L_rule%Qdof

      do i = 1, ndim
         w_jump(i, 1:Qdof ) = matmul( w_jumpDof(i,1:dof) , L_rule%phi( 1:dof,1:Qdof) )
      end do
      nullify( L_rule )

   end function evalwSTjumpInLagrangeNodes

end module stdgm_mod
