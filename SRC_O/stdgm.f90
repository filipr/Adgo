module stdgm_mod
   use blocks_integ
   use data_mod

implicit none

   public :: Eval_whST_Elem
   public :: evalWSTinTimeDof

   public :: evalL2STNorm_Elem
   public :: evalL2STNormQ_Elem
   public :: evalL2STNormVecQ_Elem
   public :: evalH1L2STNorm_Elem

   public :: evalSTfunInIntTime
   public :: evalSTfunInIntTime_Der
   public :: evalSTfunInRealTime
   public :: evalSTfunInIntTimeDof
   public :: evalSTfunInIntTimeDof_Der
   public :: evalSTfunInRealTimeDof

   public :: Eval_wSTfin_Elem
   public :: Transfer_wST_to_w_Elem
   public :: Transfer_wST_to_w_Patch




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
  !> STDG: transfers wST -> w on elem and its neighbors
  !> in node timeNode of time quadrature, tQnum = number of quadrature nodes
  subroutine Transfer_wST_to_w_Patch(elem, timeNode, tQnum)
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
      !print*, 'Qdof', Qdof, 'TQnum', elem%TQnum, 'QNUM', elem%Qnum

       elem%wSTfin(:,:) = 0.
       elem%wSTfinAD(:,:) = 0.

       do i = 1,ndim
         !new version for adaptation
          elem%wSTfinAD(i,1:dof) = matmul(elem%wST(i,1:dof,1:Tdof) , T_rule%phi(1:Tdof,0))

          do j =1,Qdof
              !new
              elem%wSTfin(i,j) = dot_product(V_rule%phi(1:dof,j), elem%wSTfinAD(i,1:dof) )
             !old
             !elem%wSTfin(i,j) = dot_product(V_rule%phi(1:dof,j), &
             !     matmul(elem%wST(i,1:dof,1:Tdof) , T_rule%phi(1:Tdof,0)) )

          enddo

       enddo
      ! print*, 'elem%wSTfin(1,1:Qdof)', elem%wSTfin(1,1:Qdof)


      !the last parameter is not used in this case - for every Tdeg in elem%phi(:,0) is the final time of the interval
      !call Transfer_wST_to_w_Elem(elem, 0, 1)


  end subroutine Eval_wSTfin_Elem


 !> STDG: evaluate \f$ wh_(t)\f$ on element elem in space integration nodes in the node node of Time quadrature with tQnum nodes
  subroutine Eval_whST_Elem( elem, tQnum, node, wi, Dwi )
    type(element), intent (inout) :: elem
    integer, intent (in) :: tQnum
    integer, intent (in) :: node
    real, dimension( 1:elem%Qdof,1:ndim), intent (inout) :: wi
    real, dimension( 1:elem%Qdof,1:ndim,1:nbDim), intent (inout) :: Dwi
    class(Time_rule), pointer :: T_rule
    type(volume_rule), pointer  :: V_rule
    integer :: Qdof, Tdof,dof, i, j, k, l
    real, dimension( 1:elem%dof) :: temp
    real, dimension(1:elem%dof,1:nbDim,1:elem%Qdof) :: Der

    T_rule => state%time%T_rule( tQnum )
    V_rule => state%space%V_rule( elem%Qnum )



       call Eval_Dphi(elem, elem%dof,  Der)
       Qdof = elem%Qdof
       Tdof = elem%Tdof
       dof = elem%dof

       wi(:,:) = 0.
       Dwi(:,:,:) = 0.


       do i = 1,ndim
          do j =1,Qdof
             temp(1:dof) = matmul( elem%wST(i,1:dof,1:Tdof), T_rule%phi(1:Tdof,node) )
             wi(j,i) = dot_product(V_rule%phi(1:dof,j), temp(1:dof) )
             do k = 1, nbdim
                Dwi(j,i,k) = dot_product(Der(1:dof,k,j), temp( 1:dof) )
             enddo
          enddo
       enddo


  end subroutine Eval_whST_Elem

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



end module stdgm_mod
