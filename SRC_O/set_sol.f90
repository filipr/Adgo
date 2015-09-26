!> compute the computational error in different norms if the exact solution is known
module set_solution
  use mesh_oper
  use main_data
  use eval_sol
  use model_oper
  use basis
  use io_sub

  implicit none

  public:: SetElementsIC_AMA
  public:: SetElementsIC_DGM
  public:: SetOneElementConstantIC
  public:: SetOneElementAnalyticalICOLD
  public:: SetOneElementIC
  public:: SetElementsIC

  public:: SetConstBC
  public:: ElementPolynProj
  public:: CheckPhysicalProperties
contains

 !> setting of the constant initial conditions by wi for each element
  subroutine SetOneElementConstantIC(elem, wi)
   type(element), intent (inout) :: elem     ! element where the basis is given
   real, dimension(1:ndim), intent(in) :: wi
   integer :: dof, nsl, j, Tdeg

   !if(ndim /= 4) then
   !   print*,'Trouble in SetOneElementConstantIC',ndim
   !   stop
   !endif

   dof = elem%dof
   nsl = ndim * dof

   Tdeg = state%time%deg+1
   if(state%time%disc_time == 'STDG') Tdeg = 1

   !!!allocate( elem%w(0:state%time%deg+1, 1:nsl))
   elem%w(0: Tdeg,1:nsl) = 0.

   ! state%space%V_rule(state%space%Qdeg(0,1))%phi(1,1) is the value of the constant test functio
   do j=1,ndim
      elem%w(0:Tdeg, (j-1)*dof + 1) = wi(j) / (state%space%V_rule( state%space%Qdeg(0,1) )%phi(1,1) )
      !write(*,'(a6,2i5,40es12.4)') 'EDW',elem%i,j, elem%w(0, (j-1)*dof + 1)
   enddo

 end subroutine SetOneElementConstantIC

 !> setting of the constant initial conditions by
 !> analytical function for each element
 subroutine SetOneElementAnalyticalICOLD(elem)
   type(element), intent (inout) :: elem     ! element where the basis is given
   real, dimension(:), allocatable :: wi, qi
   real, dimension(:,:), allocatable :: xi, x
   real, dimension(:,:), allocatable :: f
   integer :: dof, nsl, k, l, Qnum, Qdof, deg

   if(ndim /= 4) then
      print*,'Trouble in SetOneElementAnalyticalIC',ndim
      stop
   endif

   deg = elem%deg
   dof = elem%dof
   nsl = ndim * dof
   Qnum = elem%Qnum
   Qdof = state%space%V_rule(Qnum)%Qdof

   !print*,'##',elem%i, deg, dof, nsl, Qnum, Qdof

   !!!allocate( elem%w(0:state%time%deg+1, 1:nsl))

   allocate( wi(1:dof) )
   allocate( qi(1:dof) )
   allocate( f(1:ndim, 1:Qdof) )
   allocate( xi(1:Qdof, 1:nbDim))
   allocate( x(1:Qdof, 1:nbDim))


   ! setting of the function \f$ f \f$ in integ nodes x(:, 1:nbDim)
   xi(1:Qdof,1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)

   call ComputeF(elem, Qdof, xi(1:Qdof,1:nbDim), x(1:Qdof, 1:nbDim) )

   if( state%type_IC .eq. 2) then ! double mach reflection
      do l=1,Qdof
         call DMR_IC(x(l,1:nbDim), f(1:4, l), state%time%ctime)
      enddo

   elseif( state%type_IC .eq. 6 ) then ! Sod shock tube
      do l=1,Qdof
         call Exact_Sod(x(l,1:nbDim), f(1:4, l), state%time%ctime )
      enddo

   elseif( state%type_IC .eq. 4) then ! Ringleb flow
      do l=1,Qdof
         call Exact_Ringleb(x(l,1:nbDim), f(1:4, l) )
      enddo

   elseif( state%type_IC .eq. 5) then ! Gaussian pulse
      do l=1,Qdof
         call Exact_Gaussian_pulse(x(l,1:nbDim), f(1:4, l) )
      enddo

   elseif( state%type_IC .eq. 9) then !  "exact" steady state in channel
      do l=1,Qdof
         call Exact_SteadyChannel(x(l,1:nbDim), f(1:4, l) )
      enddo

   else
      do l=1,Qdof
         f(1, l) = x(l,1)
         f(2, l) = x(l,2)*(1. - x(l,2))
         f(3, l) = 0.
         f(4, l) = 5. + 0.1*x(l,1) + 0.2*x(l,2)

         !!!write(*,'(4es12.4)') xi(l, 1:nbDim), x(l, 1:nbDim)
      enddo
   endif

   !!write(*,'(a2, a2, 20es14.6)') '..',elem%RGtype, f(1,:)

   ! evaluating of the coefficients of basis expansion
   do k=1, ndim
      qi(:) = 0.
      call IntegrateVectorB(elem, dof, f(k, 1:Qdof), qi(1:dof) )

      !if( elem%i == 1) write(*,'(a2, a2, 24es14.6)') ',,',elem%RGtype, qi(1:dof)

      do l=1,dof
         wi(l) = dot_product(elem%MassInv%Mb(l,1:dof), qi(1:dof) )
      enddo
      ! P_0 approximation
      !wi(2:dof) = 0.

      do l=0, state%time%deg+1
         elem%w(l, (k-1)*dof + 1: k*dof) = wi(1:dof)
      enddo

   enddo
   deallocate(wi, f, xi, x, qi)

 end subroutine SetOneElementAnalyticalICOLD


 !> setting of the initial conditions, analytical function for each element
 subroutine SetOneElementIC(elem, all_time_levels)
   type(element), intent (inout) :: elem     ! element where the basis is given
   logical, intent(in) :: all_time_levels
   real, dimension(:), allocatable :: wi, qi
   real, dimension(:,:), allocatable :: xi, x
   real, dimension(:,:), allocatable :: f
   integer :: dof, ndof, k, l, Qnum, Qdof, deg, nlev

   nlev = 0
   if(all_time_levels) nlev = state%time%deg+1
   if(state%time%disc_time == 'STDG') nlev = 0   ! ST DGM

   deg = elem%deg
   dof = elem%dof
   ndof = ndim * dof
   Qnum = elem%Qnum
   Qdof = state%space%V_rule(Qnum)%Qdof

   !if (.not. associated(elem%w)) then
   !if (size(elem%w)==0) then
   !  allocate( elem%w(0:state%time%deg+1, 1:ndof))
   !endif

   allocate( wi(1:dof) )
   allocate( qi(1:dof) )
   allocate( f(1:Qdof, 1:ndim) )
   allocate( xi(1:Qdof, 1:nbDim))
   allocate( x(1:Qdof, 1:nbDim))


   ! setting of the function \f$ f \f$ in integ nodes x(:, 1:nbDim)
   xi(1:Qdof,1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)

   call ComputeF(elem, Qdof, xi(1:Qdof,1:nbDim), x(1:Qdof, 1:nbDim) )

   call Exact_Sol(Qdof, x(1:Qdof, 1:nbDim), f(1:Qdof, 1:ndim), state%time%ctime  )

   ! evaluating of the coefficients of basis expansion
   do k=1, ndim
      qi(:) = 0.
      call IntegrateVectorB(elem, dof, f(1:Qdof, k), qi(1:dof) )

      do l=1,dof
         wi(l) = dot_product(elem%MassInv%Mb(l,1:dof), qi(1:dof) )
      enddo

      do l=0,nlev
         elem%w(l, (k-1)*dof + 1: k*dof) = wi(1:dof)
      enddo

      !!write(300+state%time%iter,*) elem%xc(1:nbDim),elem%w(0,1:3)

      !if(elem%i == 20) write(*,'(a4,i4,80es12.4)') '!!w', &
      !     elem%i,elem%xc(1:nbDim),elem%w(0,1:3)
      !if(elem%i == 20) write(*,'(a4,i4,80es12.4)') '!!f', &
      !     elem%i,f(1,:)
      !if(elem%i == 20) write(*,'(a4,i4,80es12.4)') '!!q', &
      !     elem%i,qi(:)

   enddo
   deallocate(wi, f, xi, x, qi)

 end subroutine SetOneElementIC

 !> setting of the initial conditions from the separate files
 !> for the level set methods
 subroutine SetElementsIC_LevelSet()
   class(element), pointer :: elem     ! element where the basis is given
   type(volume_rule), pointer :: V_rule
   type(Gauss_rule), pointer :: G_rule
   real, dimension(:,:), pointer:: Rphi ! pointers to test functions
   real, dimension(:,:), allocatable :: wi, w, xi
   character(len=15) :: itext
   real :: val, area
   integer :: i, j, k, max_dof, Qdof
   integer :: iflow=11, ilevel =12, ivec = 13

   allocate(wi(1:60, 1:3), w(1:60, 1:3) )
   open(iflow, file='data.flow', status = 'OLD')
   open(ilevel, file='data.vof', status = 'OLD')
   open(ivec, file='data.vec', status = 'UNKNOWN')

   read(iflow,*) itext
   read(iflow,*) itext
   read(iflow,*) itext

   read(iflow, *) state%time%FinTime, state%LevelSet_Area, state%isol
   !!!state%LevelSet_Area =  state%LevelSet_Area * 0.9

   !print*,' state%time%FinTime over-written:',  state%time%FinTime,';   state%LevelSet_Area=', state%LevelSet_Area

   read(ilevel, *)

   do i=1, grid%nelem
      !do i=1,1 !grid%nelem
      elem => grid%elem(i)

      V_rule =>  state%space%V_rule(elem%Qnum)

      ! do k=1,elem%flen
      !    do j=1,elem%face(fGnum, k)
      !       write(52,'(40es12.4)') elem%xi(k, j , 1:2)
      !    enddo
      ! enddo


      max_dof = max(maxval( elem%face(fGnum,:)), elem%Qdof,  6)
      allocate(xi(1: max_dof, 1:2) )

      ! reading the flow velocity in P2- Lagrange nodes
      !read(iflow, *) wi(1:6, 1)
      !read(iflow, *) wi(1:6, 2)

      read(iflow, *) wi(1,1), wi(3,1), wi(6,1), wi(2,1), wi(5,1), wi(4,1)  !SPECIAL ORDERING
      read(iflow, *) wi(1,2), wi(3,2), wi(6,2), wi(2,2), wi(5,2), wi(4,2)  !SPECIAL ORDERING
      read(ilevel,*) wi(1,3), wi(3,3), wi(6,3), wi(2,3), wi(5,3), wi(4,3)  !SPECIAL ORDERING


      ! ! ONLY TEST
       xi(1, 1:2) = grid%x(elem%face(idx,1), 1:2)
       xi(3, 1:2) = grid%x(elem%face(idx,2), 1:2)
       xi(6, 1:2) = grid%x(elem%face(idx,3), 1:2)
       xi(2, 1:2) = (xi(1,1:2) + xi(3,1:2))/ 2
       xi(4, 1:2) = (xi(1,1:2) + xi(6,1:2))/ 2
       xi(5, 1:2) = (xi(3,1:2) + xi(6,1:2))/ 2


       ! TEST for rotace
       !do j=1,6
       !   wi(j, 1) = -(xi(j, 2) -0.5)
       !   wi(j, 2) =  (xi(j, 1)-0.5)
       !enddo


       !do j=1,6
       do j=4,4
          val = dot_product(wi(j , 1:2), wi(j , 1:2))**0.5
          if(val > 0) val = elem%diam *0.15 / val
          !if(val > 0) val = val *0.1
          write(ivec,*) xi(j,1:2)
          write(ivec,*) xi(j,1:2) + wi(j , 1:2)*val
          write(ivec, '(x)')
       enddo
      ! ! END OF TEST

       !do j=1,6
       !   write(*,'(a6,40es12.4)') 'wi:',wi(j, 1:2), xi(j, 2)-0.5, -(xi(j,1) -0.5)
       !enddo
       !print*,'----------------------------------------------------'
      !!write(*,'(a6,40es12.4)') 'wi:',wi(1:6, 2)

      ! transforming  of the flow velocity to DG-basis
      !!call Lagr2Basis(elem, 2, wi(1:6, 1), elem%Qnum, 6, w(1:6,1) )
      !!call Lagr2Basis(elem, 2, wi(1:6, 2), elem%Qnum, 6, w(1:6,2) )
      call Lagr2BasisDiff(elem, 2, wi(1:6, 1), elem%Qnum, 6, w(1:6,1) )
      call Lagr2BasisDiff(elem, 2, wi(1:6, 2), elem%Qnum, 6, w(1:6,2) )
      call Lagr2BasisDiff(elem, 2, wi(1:6, 3), elem%Qnum, elem%dof, w(1:elem%dof,3) )

      !!write(*,'(a6,40es12.4)') 'WW:',w(1:6, 1)
      !!write(*,'(a6,40es12.4)') 'WW:',w(1:6, 2)

      ! setting of the initial condition
      elem%w(0, 1:elem%dof)  = w(1:elem%dof,3)

      ! evaluation of the flow velocity in volume integ nodes
      Qdof = elem%Qdof
      xi(1:Qdof, 1:2) = elem%xi(0,1:Qdof, 1:2)
      do j=1,Qdof
         elem%xi(0, j, 1) = dot_product(V_rule%phi(1:6,j), w(1:6,1))
         elem%xi(0, j, 2) = dot_product(V_rule%phi(1:6,j), w(1:6,2))
      enddo

      !!do j=1,Qdof
      !!   write(51,'(40es12.4)') xi(j,1:2), elem%xi(0, j , 1:2)
      !!enddo

      ! evaluation of the flow velocity in edge integ nodes
      do k=1,elem%flen
         Qdof = elem%face(fGnum,k)
         G_rule =>  state%space%G_rule(Qdof)

         if(elem%HGnode) then
            Rphi => G_rule%phi(elem%type, elem%HGface(1, k), elem%HGface(2, k), 1:6, 1:Qdof)
         else
            Rphi => G_rule%phi(elem%type, k, 1, 1:6, 1:Qdof)
         endif


         xi(1:Qdof, 1:2) = elem%xi(k,1:Qdof, 1:2)
         do j=1,Qdof
            elem%xi(k, j, 1) = dot_product(Rphi(1:6,j), w(1:6,1))
            elem%xi(k, j, 2) = dot_product(Rphi(1:6,j), w(1:6,2))
         enddo

         !!do j=1,Qdof
         !!   write(52,'(40es12.4)') xi(j,1:2), elem%xi(k, j , 1:2)
         !!enddo

      end do

      deallocate(xi)
   enddo

   close(iflow)
   close(ilevel)
   close(ivec)

   deallocate(w, wi)

   call Comp_LevelSet_Area(area, val )
   write(*,'(a50,2es14.6,es10.2,a12,es12.4)' ) &
        ' # LEVEL SETs Areas (comput prescrib rel_diff): ', &
        area, state%LevelSet_Area, abs(area-state%LevelSet_Area)/area, &
        ',  Fin Time=',state%time%FinTime

   !!call WriteOutput_LevelSet()

   !print*, 'stopped in SetElementsIC_LevelSet'
   !stop

 end subroutine SetElementsIC_LevelSet

 !> output of the level set function to the P2 Lagrang nodes
 subroutine WriteOutput_LevelSet()
   class(element), pointer :: elem     ! element where the basis is given
   type(Lagrang_rule), pointer :: L_rule
   real, dimension(:,:), pointer:: phi ! pointers to test functions
   real, dimension(:,:), allocatable :: wi, q
   character(len=15) :: itext
   real :: area, length, shift
   integer :: i,  k,l, Qnum, Qdof, ist, dof, lev
   integer :: ilevel =12

   ! modification of the level set function in order to fit the volume
   print*,'LEVEL SET: modification of the level set function in order to preserve the volume'
   call Comp_LevelSet_Area( area, length )
   write(*,*) '                rel. diff.  computed area   prescibed area   shift'
   write(*,'(a10, i2,a3, es12.4, 2es16.8)' ) &
        ' Fitting(', 0,') :', abs(state%LevelSet_Area - area)/ area, area, state%LevelSet_Area

   do lev=1, 8
      !if( abs(state%LevelSet_Area - area)/ area > 1E-2) then
      !if( abs(state%LevelSet_Area - area)/ area > 1E-3) then
      if( abs(state%LevelSet_Area - area)/ area > 1E-4) then
         shift = (area -  state%LevelSet_Area) / length
         shift = shift*0.5

         do i=1, grid%nelem
            elem => grid%elem(i)
            do k=1,ndim
               ist = (k-1)*dof + 1
               elem%w(0,ist) = elem%w(0,ist) - shift
            enddo
         enddo

         call Comp_LevelSet_Area( area, length )
         write(*,'(a10, i2,a3, es12.4, 2es16.8, 6es12.4)' ) &
              ' Fitting(',lev,') :', abs(state%LevelSet_Area - area)/ area, area, &
              state%LevelSet_Area, shift

      else
         goto 10
      end if
   enddo
10 continue

   print*,'Output for LEVEL SET METHOD in file "data.vofx"'
   allocate(wi(1:ndim, 1:state%space%max_dof), q(1:ndim, 1:state%space%max_Qdof))

   open(ilevel, file='data.vofx', status = 'UNKNOWN')
   write(ilevel, *)' ## Scalar levelse in format phiA phiB phiC phiAB phiBC phiCA '

   Qnum = 2
   L_rule =>  state%space%L_rule(Qnum)
   Qdof = L_rule%Qdof
   phi => L_rule%phi(1:state%space%max_dof, 1:Qdof)


   do i=1, grid%nelem
   !do i=1,1 !grid%nelem
      elem => grid%elem(i)

      dof = elem%dof

      do k=1,ndim
         ist = (k-1)*dof + 1
         wi(k, 1:dof) = elem%w(0,ist:ist+dof-1)

         ! evaluation of w in the Langrangian nodes
         do l=1, Qdof
            q(k, l) = dot_product( wi(k, 1:dof), phi(1:dof,l) )
         enddo

         !if(elem%i <= 20) then
         !write(*,'(a4,80es12.4)') 'wi.', q(1, 1:Qdof),wi(k,1:dof)
         !endif

      enddo

      ! output to file
      do k=1,ndim
         write(ilevel,'(6es16.8)') q(k,1),q (k,3), q(k,6), q(k,2), q(k,5), q(k,4)
      enddo
   enddo

   deallocate(wi, q)
   close(ilevel)

   !print*, 'stopped in  WriteOutput_LevelSet'
   !stop
 end subroutine WriteOutput_LevelSet




 !> compute the area given by \f$ \varphi \ge 0 \f$, \f$ \varphi \f$ is the Level Set function
 subroutine Comp_LevelSet_Area(area, length)
   real, intent(inout) :: area, length
   class(element), pointer :: elem     ! element where the basis is given
   type(Lagrang_rule), pointer :: L_rule
   real, dimension(:,:), pointer:: phi ! pointers to test functions
   real, dimension(:,:), allocatable :: wi, q
   integer :: i,  k,l, Qnum, Qdof, ist, dof, ipoc
   integer :: ilevel =12
   real :: area_part

   !print*,'start of Comp_LevelSet_Area( ) '
   Qnum = 2
   L_rule =>  state%space%L_rule(Qnum)
   Qdof = L_rule%Qdof
   phi => L_rule%phi(1:state%space%max_dof, 1:Qdof)

   allocate(wi(1:ndim, 1:state%space%max_dof), q(1:ndim, 1:state%space%max_Qdof))

   area = 0.
   length = 0.
   ipoc = 0

   do i=1, grid%nelem
   !do i=1,1 !grid%nelem
      elem => grid%elem(i)

      dof = elem%dof

      do k=1,ndim
         ist = (k-1)*dof + 1

         wi(k, 1:dof) = elem%w(0,ist:ist+dof-1)

         ! evaluation of w in the Langrangian nodes
         do l=1, Qdof
            q(k, l) = dot_product( wi(k, 1:dof), phi(1:dof,l) )
         enddo

         if(minval(q(k, 1:Qdof))  > 0.) then
            area = area + elem%area
            !write(141,*) elem%xc(1:2)
         elseif(maxval(q(k, 1:Qdof))  > 0.) then
            length = length + elem%diam / 2**(1.5)
            ipoc = ipoc + 1

            ! P1 approximation
            !call Eval_LevelSet_Area_part( q(k,1), q(k,3), q(k,6), area_part )
            !area = area + elem%area * area_part

            ! P2 approximation
            call Eval_LevelSet_Area_part( q(k,1), q(k,2), q(k,4), area_part )
            area = area + elem%area * area_part / 4

            call Eval_LevelSet_Area_part( q(k,2), q(k,3), q(k,5), area_part )
            area = area + elem%area * area_part / 4

            call Eval_LevelSet_Area_part( q(k,2), q(k,5), q(k,4), area_part )
            area = area + elem%area * area_part / 4

            call Eval_LevelSet_Area_part( q(k,4), q(k,5), q(k,6), area_part )
            area = area + elem%area * area_part / 4


            !write(*,'(a6,i5,16es12.4)') '#GQ#',i, area, area_part, minval(q(1, 1:Qdof)) , q(1, 1:Qdof)
            !write(142,*) elem%xc(1:2)
         endif

      enddo
      !if(minval(q(k, 1:Qdof))  > 0.) &
      !     write(*,'(a6,i5,16es12.4)') '#ED#',i, area, length, minval(q(1, 1:Qdof)) , q(1, 1:Qdof)
   enddo

   !write(*,'(a20,8es12.4)' ) ' Level Set Area = ', &
   !     area,  pi *0.25**2, abs(area -  pi *0.25**2), length, pi /2, 1.*ipoc, elem%diam
   !print*,'end of Comp_LevelSet_Area( ) '

 end subroutine Comp_LevelSet_Area

 !> evaluate of the appropriate part of the area of element given by the level set
 subroutine Eval_LevelSet_Area_part( q1, q2, q3, area_part )
   real, intent(in) :: q1, q2, q3
   real, intent(out) :: area_part
   real :: l1, l2


   if(q1 >= 0.) then
      if(q2 >= 0.) then
         if(q3 >= 0.) then
            ! q1 > 0, q2 > 0,  q3 > 0
            !print*, 'AREA_part:      ','EMPTY (1)'
            area_part = 1
         else
            ! q1 > 0, q2 > 0,  q3 < 0
            l1 = q3 / (q3 - q1)
            l2 = q3 / (q3 - q2)
            area_part = 1. - l1*l2
         endif

      else  ! q2 < 0
         if(q3 >= 0.) then
            ! q1 > 0, q2 < 0,  q3 > 0
            l1 = q2 / (q2 - q1)
            l2 = q2 / (q2 - q3)
            area_part = 1. - l1*l2

         else
            ! q1 > 0, q2 < 0,  q3 < 0
            l1 = q1 / (q1 - q2)
            l2 = q1 / (q1 - q3)
            area_part = l1*l2
         endif
      endif

   else ! q1 < 0
      if(q2 >= 0.) then
         if(q3 >= 0.) then
            ! q1 < 0, q2 > 0,  q3 > 0
            l1 = q1 / (q1 - q2)
            l2 = q1 / (q1 - q3)
            area_part = 1. - l1*l2

         else
            ! q1 < 0, q2 > 0,  q3 < 0
            l1 = q2 / (q2 - q1)
            l2 = q2 / (q2 - q3)
            area_part = l1*l2
         endif

      else  ! q2 < 0
         if(q3 >= 0.) then
            ! q1 < 0, q2 < 0,  q3 > 0
            l1 = q3 / (q3 - q1)
            l2 = q3 / (q3 - q2)
            area_part = l1*l2

         else
            ! q1 < 0, q2 < 0,  q3 < 0
            !print*, 'AREA_part:      ','EMPTY (2)'
            area_part = 0

         endif
      endif
   endif



 end subroutine Eval_LevelSet_Area_part


 !> setting of the constant initial conditions from the boundary conditions
 !> prescribed on inlet
 subroutine SetElementsIC()
   integer :: dof
   integer :: i, j, k, ibc, ipc, ibb
   real:: val, val1, xmin, xmax
   real, dimension(1:ndim) :: wi
   real, dimension(:,:), allocatable :: vv
   class(element), pointer :: elem
   real :: pL, pR, pp

   state%time%ctime = state%time%ttime

   if(.not. state%SP) then ! not saddle point problem

      if(state%modelName == 'scalar' .or. state%modelName == '2eqs') then     !!! SCALAR case
         if(state%type_IC .eq. 5) then
            ! LEVEL SET METHOD
            call SetElementsIC_LevelSet()

         else
            do i=1,grid%nelem
               call SetOneElementIC(grid%elem(i), .true. )
            enddo
         endif

      elseif(state%modelName == 'NSe' ) then !! Navier-Stokes or Euler equations

         if(state%type_IC .eq. 1) then  ! standard IC from BC
            ibc = -1
            do i=1, state%numBC
               if(state%BC(i)%ibc == 1) ibc = i
            enddo
            if(ibc <= 0) then
               print*,' Inlet BC in Set IC does not found !!!'
               stop
            endif

            do i=1,grid%nelem
               call SetOneElementConstantIC(grid%elem(i), state%BC(ibc)%ww(:) )
            enddo

         elseif( state%type_IC .eq. 3) then    ! 2D explosion

            do i=1,grid%nelem
               if(dot_product(grid%elem(i)%xc(:), grid%elem(i)%xc(:)) .le. 1) then
                  wi(1:4) = (/ 1.,  2.,  1.,  0.5 /)
               else
                  wi(1:4) = (/ 1.,  0.,  0.,  0.5 /)
               endif

               call SetOneElementConstantIC(grid%elem(i), wi(1:4) )
            enddo

         elseif( state%type_IC .eq. 10) then    ! SE1050
            allocate(vv (1:4, 1:ndim) )

            vv(1, 1:ndim) = state%BC(1)%ww(1:ndim)
            vv(2, 1:ndim) = state%BC(2)%ww(1:ndim)

            ipc = grid%iper(1,1)
            xmin = 1E+20
            xmax = -1E-20

            do i=1,grid%nelem
               elem => grid%elem(i)

               val = 1E+20

               do k=1,grid%nbelm ! seeking the orientation of the nearest periodic boundary

                  if(grid%b_edge(k)%ibc == ipc) then
                     val1 = abs(elem%xc(1) - grid%x(grid%b_edge(k)%lbn(1), 1) )
                     if(val1 < val) then
                        val = val1
                        ibb = k
                     endif
                     xmin = min(xmin,  grid%x(grid%b_edge(k)%lbn(1), 1) )
                     xmax = max(xmax,  grid%x(grid%b_edge(k)%lbn(1), 1) )
                  endif

               enddo

               ! tangential orientation of the nearest boundary edge
               vv(4, 1:2) = grid%x(grid%b_edge(ibb)%lbn(1), 1:2) - grid%x(grid%b_edge(ibb)%lbn(2), 1:2)
               val = VectorNorm(vv(4, 1:2) )
               vv(4, 1:2) = vv(4, 1:2) / val
               if(vv(4, 1) < 0) vv(4, 1:2) = - vv(4, 1:2)

               ! linear interpolation of the BC
               vv(3, 1:ndim) = vv(1, 1:ndim) &
                    + (vv(1, 1:ndim) - vv(2, 1:ndim) ) * (elem%xc(1)  - xmin)/(xmax - xmin)
               pL = pressure(4, vv(1,:))
               pR = pressure(4, vv(2,:))
               pp = pL + (pR - pL) * (elem%xc(1)  - xmin)/(xmax - xmin)
               vv(3, 4) = pp / state%model%kappa1 + 0.5*( vv(3, 2)**2 + vv(3, 3)**2) / vv(3, 1)

               ! setting the tangential orientation of the velocity
               val = VectorNorm(vv(3, 2:3) )
               vv(3, 2:3) = vv(4, 1:2) * val

               !if(elem%xc(1) < -0.055) then
               ! if(elem%xc(1) > 0.13) then
               !    write(*,'(a6,15es14.6)') 'vv1', vv(1,:), VectorNorm(vv(1, 2:3)),pressure(4, vv(1,:))
               !    write(*,'(a6,15es14.6)') 'vv2', vv(2,:), VectorNorm(vv(2, 2:3)),pressure(4, vv(2,:))
               !    write(*,'(a6,15es14.6)') 'vv3', vv(3,:), VectorNorm(vv(3, 2:3)),pressure(4, vv(3,:))
               !    write(*,'(a6,15es14.6)') 'tt', vv(4,1:2), VectorNorm(vv(4, 1:2))
               !    print*,'-------'
               ! endif

               !write(70,*) elem%xc(1:2)
               !write(70,*) elem%xc(1:2)+ vv(3, 2:3)*elem%diam/4
               !write(70,*) elem%xc(1:2)+ vv(4, 1:2)*0.01

               !write(70,*) grid%x(grid%b_edge(ibb)%lbn(1), 1:2) , ibb
               !write(70,*)

               call SetOneElementConstantIC(grid%elem(i), vv(3, 1:4) )

               !if(elem%xc(1) < -0.055) write(*,'(a6,12es12.4)') 'ic:',grid%elem(i)%w(0,:)
            enddo

            !stop 'EDHTSGWTSAYS'

         else ! other cases
            do i=1,grid%nelem
               call SetOneElementIC(grid%elem(i), .true. )
            enddo

         endif

      elseif(state%modelName == 'wet_steam' ) then !! Navier-Stokes or Euler equations

         if(state%type_IC .eq. 1) then  ! standard IC from BC
            ibc = -1
            do i=1, state%numBC
               if(state%BC(i)%ibc == 1) ibc = i
            enddo
            if(ibc <= 0) then
               print*,' Inlet BC in Set IC does not found !!!'
               stop
            endif

            do i=1,grid%nelem
               call SetOneElementConstantIC(grid%elem(i), state%BC(ibc)%ww(1:ndim) )
            enddo

         else
            print*,'Not yet implemented in set_sol.f90 WESQ wet steam'
            stop
         endif
      else
         print*,' set_sol.f90: (3) SetElementsIC, no IC for  ndim  =',ndim, state%modelName
         stop
      endif


   else ! SADDLE POINT problems

      if(state%modelName == 'incNS' ) then !! incompressible Navier-Stokes  equations

         if(state%type_IC .eq. 1) then  ! standard IC from BC
            ibc = -1
            do i=1, state%numBC
               if(state%BC(i)%ibc == 1) ibc = i
            enddo
            if(ibc <= 0) then
               print*,' Inlet BC in Set IC does not found !!!'
               stop
            endif

            do i=1,grid%nelem
               elem => grid%elem(i)
               dof = max(elem%dof, elem%dofP)

!!!allocate( elem%wSP(:,:,1:dof))
               elem%wSP(wV1:wP,0,1:dof) = 0.

               ! state%space%V_rule(state%space%Qdeg(0,1))%phi(1,1) is the value of the constant test functio
               elem%wSP(wV1, 0, 1) = state%BC(ibc)%ww(1) / (state%space%V_rule( state%space%Qdeg(0,1) )%phi(1,1) )
               elem%wSP(wV2, 0, 1) = state%BC(ibc)%ww(2) / (state%space%V_rule( state%space%Qdeg(0,1) )%phi(1,1) )
               elem%wSP(wP,  0, 1) = state%BC(ibc)%ww(3) / (state%space%V_rule( state%space%Qdeg(0,1) )%phi(1,1) )

               do k=1, state%time%deg
                  elem%wSP(wV1:wP, k, 1) =  elem%wSP(wV1:wP, 0, 1)
               enddo
               !write(*,'(a6,2i5,40es12.4)') 'EDW',elem%i,j, elem%w(0, (j-1)*dof + 1)
            enddo

         else
            print*,'Not yet implemented in set_sol.f90 WESQ'
            stop
         endif
      else
         print*,' set_sol.f90: SetElementsIC, no IC for  ndim  =',ndim, state%modelName
         stop
      endif
   end if

 end subroutine SetElementsIC

 !> reading the initial values from the file 'resultsx' obtained by ANGENER,
 !> \f$ P_0 \f$ approximation
 subroutine SetElementsIC_AMA()
   integer :: ifile=12
   integer :: i
   real, dimension(:), allocatable:: wi

   allocate(wi(1:ndim) )

    open(ifile, file='resultsx', status='OLD')

    do i=1,grid%nelem
       read(ifile,*) wi(1:ndim)

       call SetOneElementConstantIC(grid%elem(i), wi(1:ndim) )
    enddo

    deallocate(wi)

    close(ifile)
  end subroutine SetElementsIC_AMA

  !> reading the initial values from the file 'dgm.sol' containing solution
  !> in langangian nodes,
  !> \f$ P_k,\ k\ge 1 \f$ approximation, only for triangular grids !!!
  subroutine SetElementsIC_DGM()
    class(element), pointer :: elem        ! elem = element
    integer :: ifile=12
    integer :: i, k, l, deg, dof, nelem, Qnum, Qdof, iBC, iIN
    integer :: Ndeg, Ndof, RG_lev, ndim1
    real :: tau
    real, dimension(:,:), allocatable :: w, wi
    !real, dimension(:,:), allocatable :: psi, xi
    !real, dimension(:), allocatable :: w, wi, qi, f

    open(ifile, file='dgm.sol', status='OLD')


    read(ifile,*) nelem, ndim1, state%time%ttime,  tau, state%time%iter
    if(.not.  state%time%tau_fixed) then
       state%time%tau_new = tau
       state%time%tau_old = tau
    endif

    if(nelem .ne. grid%nelem .or. ndim1 .ne. ndim) then
       if(state%modelName == "wet_steam") then   ! for wet steam ndim == 8 and ndim1 == 4
         if(ndim /= 8 .or. ndim1 /= 4) then
           write(*,*) "Problem due to ndim /= 8 or ndim1 /= 4 in SetElementsIC_DGM() for wet steam"
           write(*,*) "ndim = ", ndim, "ndim1 = ", ndim1
           stop
         endif
       else
         print*,'Inconsistency between grid and grid.sol in  SetElementsIC_DGM()'
         print*, nelem, '=?',  grid%nelem,' | ', ndim1,'=?', ndim
         stop
       endif
    endif

    allocate( w(1:ndim, &
         1:(MaxDegreeImplemented + 1)*(MaxDegreeImplemented + 1)/2 ) )

    allocate(wi(1:ndim, &
         1:(MaxDegreeImplemented + 1)*(MaxDegreeImplemented + 1)/2 ) )

    if(state%modelName == 'wet_steam') then ! wet steam
      do iBC=1,state%numBC
        if(state%BC(iBC)%inout == 0) iIN = iBC  ! remember index for inlet
      end do
    endif

    do i=1,grid%nelem
       elem => grid%elem(i)

       do k=1, ndim1

          read(ifile,*) Ndeg, RG_lev, w(k, 1:(Ndeg+1)*(Ndeg+2)/2)

          if(k == 1)  Ndof = (Ndeg+1)*(Ndeg+2)/2
          if( Ndof /= (Ndeg+1)*(Ndeg+2)/2 ) then
             print*,' Different degree for different components in dgm.sol'
             print*,' elem%i = ',i, Ndeg, Ndof
             stop
          endif
       enddo

       if(state%modelName == 'wet_steam') then
         ! Initial values for wet steam equations are taken from Inlet
         w(5, :) = state%BC(iIN)%ww(5) ! BC
         w(6, :) = state%BC(iIN)%ww(6) ! BC
         w(7, :) = state%BC(iIN)%ww(7) ! BC
         w(8, :) = state%BC(iIN)%ww(8) ! BC
       endif

       write(22,*) i, w(:,1)

       dof =  elem%dof
       Qnum =  elem%Qnum

       if(Ndeg == elem%deg) then
          ! the same degree of polynomial approximation
          call Lagr2Basis(elem, Ndeg, w(1:ndim, 1:Ndof), Qnum, dof, wi(1:ndim, 1:dof) )
       else
          call Lagr2BasisDiff(elem, Ndeg,  w(1:ndim, 1:Ndof), Qnum, dof, wi(1:ndim, 1:dof) )
       endif

       do k=1,ndim
          elem%w(0, (k-1)*dof + 1: k*dof) = wi(k, 1:dof)
       enddo

    enddo
    close(ifile)
    deallocate (w, wi)

  end subroutine SetElementsIC_DGM

  !> setting of appropriate boundary conditions
  !> evaluation of coordinates integration nodes for computation of drag and lift
  subroutine SetConstBC(grid)
    use data_mod, grid_g => grid
    class(mesh), intent(inout), target :: grid
    class(element), pointer :: elem
    real, dimension(:,:), allocatable :: xi!, Fx
    real, dimension(1:nbDim) :: xc
    integer :: i, j, ib, k, Qdof, ie, je, jc
    real :: t

    do ib=1,grid%nbelm
       ie = grid%b_edge(ib)%itc
       je = grid%b_edge(ib)%jtc   ! je = 1, ..., flen

       elem => grid%elem(ie)


       jc = je                    ! jc = 1, ..., type
       if(elem%HGnode) jc = elem%HGface(1, je)

       !! grid%elem(ie)%iBC(:) allocated in mesh.f90, seek neighbours
       !!grid%elem(ie)%iBC(:): 0= fixed walls, >0 Inlet/Outlet

       grid%b_edge(ib)%BC = 0   ! Slip
       grid%b_edge(ib)%inout = -1   ! Slip

       do i=1,state%numBC
          if(state%BC(i)%ibc == grid%b_edge(ib)%ibc) then

             elem%iBC(je) = i

             grid%b_edge(ib)%BC = 1   ! Inlet Outlet
             grid%b_edge(ib)%inout = state%BC(i)%inout   ! Inlet Outlet

             ! inout = 0   ... inlet or Dirichlet BS as walls
             ! inout = 1   ... far-field outlet or Neumann BC
             ! inout = 2   ... channel outlet
             elem%tBC(je) = state%BC(i)%inout

             ! if(state%BC(i)%inout == 1 ) then
             !   elem%tBC(je) = 0             ! inlet, Dirichlet BC as walls
             ! elseif(state%BC(i)%inout == 1 .or. state%BC(i)%inout == 2) then
             !    elem%tBC(je) = 1             ! outlet Neumann BC
             ! endif

             goto 100
          endif
       enddo

       !if(ndim /= 4) then
       !   print*,'Slip BC for non Euler equations !!!!',' ndim =',ndim
       !   !stop
       !endif

100    continue


       ! evaluation of integ.  nodes for imprermeable walls for computation of Drag. or
       ! Dirichlet BC for scalar problem
       if(grid%b_edge(ib)%BC == 0 .or. state%modelName == 'scalar' .or.state%modelName == '2eqs' .or. &
            state%type_IC == 9 ) then

          ! barycentre of reference element
          if(elem%type == 3) xc(1:nbDim) = 1./3
          if(elem%type == 4) xc(1:nbDim) = 1./2

          Qdof = elem%face(fGdof,je)

          allocate(xi(1:Qdof, 1:nbDim) )

          allocate(grid%b_edge(ib)%x_div(1:Qdof,1:nbDim) )

          if(elem%type == 3) then ! triangle
             do k=1, Qdof
                ! Gauss integration node on (0,1)
                t = state%space%G_rule(elem%face(fGnum,je))%lambda(k)
                if(elem%HGnode) t = ResizeHG(t, elem%HGface(2, je))

                ! xi(1:nbDim) ... integration node on reference element
                if(jc == 1) xi(k,1:nbDim) = (/ t, 0./)
                if(jc == 2) xi(k,1:nbDim) = (/ 1-t, t/)
                if(jc == 3) xi(k,1:nbDim) = (/ 0., 1- t/)
             enddo
          elseif(elem%type == 4) then ! quadrilateral
             do k=1, Qdof
                t = state%space%G_rule(elem%face(fGnum,je))%lambda(k)
                if(elem%HGnode) t = ResizeHG(t, elem%HGface(2, je))

                ! xi(1:nbDim) ... barycentric coordinates of integration node
                if(jc == 1) xi(k,1:nbDim) = (/ t    , 0.   /)
                if(jc == 2) xi(k,1:nbDim) = (/ 1.   , t    /)
                if(jc == 3) xi(k,1:nbDim) = (/ 1.-t , 1.   /)
                if(jc == 4) xi(k,1:nbDim) = (/ 0.   , 1.-t  /)
             enddo

          else
             print*,'Only triang and quad in " SetConstBC()" implemented', elem%type
             stop
          endif

          call ComputeF(elem, Qdof, xi, grid%b_edge(ib)%x_div)

          deallocate(xi)
       endif

       !!write(*,'(a6,8i5)') 'BC:',ib, ie, je, elem%iBC(je), grid%b_edge(ib)%BC,  grid%b_edge(ib)%inout

    enddo
  end subroutine SetConstBC


  !> polynomial projection of a function given in integ nodes of elem
  !> function f given in integ nodes, we compute its polynomial projection and set their
  !> values in integ nodes again, i.e., if f is polynomial than Pf = f
  subroutine ElementPolynProj(elem, f, Pf)
   type(element), intent (inout) :: elem     ! element where the basis is given
   real, dimension(1:elem%Qdof, 1:nbDim, 1:ndim), intent(in) :: f
   real, dimension(1:elem%Qdof, 1:nbDim, 1:ndim), intent(out) :: Pf
   real, dimension(:,:), pointer :: phi
   real, dimension(:), allocatable :: wi, qi
   integer :: dof,  k, j, l, Qnum, Qdof, deg

   deg = elem%deg
   dof = elem%dof
   Qnum = elem%Qnum
   Qdof = state%space%V_rule(Qnum)%Qdof

   phi => state%space%V_rule(elem%Qnum)%phi(1:dof,1:Qdof)

   allocate( wi(1:dof) )
   allocate( qi(1:dof) )

   ! evaluating of the coefficients of basis expansion
   do k=1, ndim
      do j=1,nbDim
         qi(:) = 0.
         call IntegrateVectorB(elem, dof, f(1:Qdof, j, k), qi(1:dof) )

         ! wi are basis coefficients in integ nodes
         do l=1,dof
            wi(l) = dot_product(elem%MassInv%Mb(l,1:dof), qi(1:dof) )
         enddo

         do l=1,Qdof
            Pf(l, j, k) = dot_product(wi(1:dof), phi(1:dof, l) )
         enddo
      enddo
   enddo
   deallocate(wi, qi)

 end subroutine ElementPolynProj

 !> NSe: check the positive presure and the density
 subroutine CheckPhysicalProperties(elem)
   class(element), pointer,intent(inout) :: elem
   real, dimension(:,:), pointer :: phi
   real, dimension(:,:), allocatable :: wi, qi, w
   real, dimension(:,:), allocatable :: Fx
   integer :: Qnum, deg, dof, j, k, l, ist, Qdof, Qlen
   integer :: ind

   ind = 0

   deg = elem%deg
   dof = elem%dof
   Qnum = elem%deg

   if(elem%type == 4) stop "Only triangles impelemnted"
   !if(elem%type == 4) Qnum = Qnum + QnumOffset

   Qdof = state%space%L_rule(Qnum)%Qdof
   Qlen = Qdof
   !if(elem%type == 4) Qlen = Qlen/2  ! quandrilaterals forms two triangles

   allocate(wi(1:ndim, 1:Qdof), qi(1:ndim, 1:Qdof), w(1:4, 1:dof) )

   !allocate(Fx(1:Qdof, 1:2) )
   !call ComputeF(elem, Qdof,  state%space%L_rule(Qnum)%lambda(1:Qdof,1:nbDim), Fx(1:Qdof, 1:nbDim) )

   phi => state%space%L_rule(Qnum)%phi(1:dof, 1:Qdof)

   do k=1,ndim
      ist = (k-1)*dof
      !wi(k, 1:dof) = elem%w(0,ist+1 : ist+dof)

      ! evaluation of w in the Langrangian nodes
      do l=1, Qdof
         wi(k, l) = dot_product(  elem%w(0,ist+1 : ist+dof), phi(1:dof,l) )
      enddo

    enddo

    ! Check of the physical properties

    !do j = 0, elem%type - 3

    call Transform_W2Q_NSe(Qdof, wi(1:4, 1:Qdof), qi(1:4,1:Qdof) )

    do l=1,Qdof
       if(qi(1,l) <= 0. .or. qi(4,l) <= 0) then

          !do j=1,Qdof
          !   write(94, *) Fx(j, 1:2), qi(4, j), qi(1:4, j)
          !enddo

          ind = ind + 1
          !write(*,'(a4,2i5,2(4es12.4,a3))' ) 'CPP',elem%i, l, wi(1:4, l),'|',qi(1:4,l)
          write(*,'(a33,2es12.4,a2,2es12.4)') &
               '!! Negative presure or density: ',elem%xc(1:2),'|',qi(1,l), qi(4,l)
          !if(qi(1,l) <= -0.1) stop
          !if(qi(4,l) <= -0.1) stop

          !if(qi(1,l) <= 0.) qi(1,l) = 0.01
          !if(qi(4,l) <= 0.) qi(4,l) = 0.01

       endif
    enddo

    if(ind > 0) then
       !write(*,'(a4,i5,80es12.4)') 'elem%w:', elem%i, elem%w(0,:)
       !call WriteProgressOutput( 'ST', .false. )
       !state%isol = state%isol + 1

       ! recomputation back to the conservative variables in Lagr. integ nodes
       call Transform_Q2W_NSe(Qdof, qi(1:4, 1:Qdof), wi(1:4,1:Qdof) )

       ! transformation to the basis coeficients
       call Lagr2Basis(elem, elem%deg, wi(1:ndim, 1:Qdof), elem%Qnum, dof, w(1:ndim, 1:dof) )

       do k=1,ndim
          elem%w(0, (k-1)*dof + 1: k*dof) = w(k, 1:dof)
       enddo

       !write(*,'(a4,i5,80es12.4)') 'elem%w:', elem%i, elem%w(0,:)
       !print*

       !call WriteProgressOutput( 'ST', .false. )

       !stop

    endif

   !enddo

    deallocate(wi, qi, w)
    ! deallocate(Fx)

 end subroutine CheckPhysicalProperties

end module set_solution
