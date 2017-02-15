module space_mod
   use adaptation_mod
   use basis_mod
   use geometry
   use integration
   use paramets
   use estims_mod

   implicit none

   type, public :: Space_t

      integer :: deg
      !integer :: degP         ! not used - IVAN
      character(len=20) :: disc_space
      character(len=20) :: estim_space
      integer :: m_IPG                     ! type of IPG stabil: NIPG(1), SIPG(-1),IIPG(0)
      real :: sigma                        ! penalty parameter
      real :: pen_deg                      ! penalty parameter sigma/ h**pen_deg
      real :: C_W                          ! penalty parameter

      real, dimension(:,:,:), allocatable ::    GScoeff    ! GramSchmidt coeffs
      integer, dimension(:), allocatable :: ldeg  ! ldeg(i) = (i+1)*(i+2)/2 + 1
      integer, dimension(:,:), allocatable :: Qdeg  ! prescribed degree of quadrature for given deg

      integer :: max_Qdof                      ! maximal of elem(i)%Qdof
      integer :: max_dof                       ! maximal of elem(i)%dof


      type(volume_rule), dimension(:), allocatable :: V_rule   ! Dunavant integ rules
      type(Gauss_rule), dimension(:), allocatable :: G_rule    ! Gauss integ rules
      type(Lagrang_rule), dimension(:), allocatable :: L_rule ! Langrangian rule (for visual)

      logical:: gridS_allocated

      real :: h                                ! maximal diameter
      real :: diam                             ! diameter of the computational domain
      real :: domain_volume                    ! meassure of the computational domain

      logical :: deg_plus                      ! compute for deg or deg + plusDeg
      integer :: plusDeg                       ! how much to increase the polynomial degree (typically 1)

      !real :: tol_min, tol_max             ! tolerances for mesh adaptation

      logical :: adaptation, adaptationR     ! are we using adaptation -> should allocate adapt
      class(Adaptation_t), allocatable :: adapt
!      class(SpaceEstims_t), allocatable :: estim

   contains
      procedure :: allocAdaptation
!      procedure :: allocEstims
      procedure :: ComputeGSCoefficients
      procedure :: ComputeBasisFunctions
      procedure :: copyGeometry
      procedure :: createQuadRules
      procedure :: init => initSpace
   !   procedure :: initAdapt
      procedure :: initDGdata
      procedure :: test_Vrules
      procedure :: Reconstruction_test

   end type Space_t

   contains

   !used while reading the main data from .ini file
   subroutine initDGdata( this, disc_space, C_W, deg )
      class(Space_t), intent( inout ) :: this
      character(len=20), intent( in ) :: disc_space
      real, intent( in ) :: C_W
      integer, intent( in ) :: deg

      this%disc_space = disc_space
      this%C_W = C_W
      this%sigma = C_W
      this%deg = deg

      if(disc_space == 'NIPG') then
         this%m_IPG = 1
      elseif(disc_space == 'IIPG') then
        this%m_IPG = 0
      elseif(disc_space == 'SIPG') then
        this%m_IPG = -1
      elseif(disc_space == 'FEM') then
        this%m_IPG = 5
        print*,'FEM NOT TESTED !!!'
        stop
      else
        print*,'Bad space discretization, only SIPG, NIPG, IIPG are implemented'
        stop
      endif

      if( deg < 0) then
        print*,'Default approximation degree has to be nonnegative ',this%deg
        stop
      elseif( deg > MaxDegreeImplemented ) then
          write(*,*) 'The default this%deg=',deg,&
               ' is higher than MaxDegreeImplemented =', MaxDegreeImplemented
          stop

      else
        write(*,'(a4,a24, a4,a7,es12.4, a4,i1)') &
             '  # ','Space discretization by ', &
             disc_space, ', c_W =', C_W, ", P_", deg
      endif


   end subroutine initDGdata

   subroutine allocAdaptation( this , tol_max, tol_min, adapt_space, max_adapt_level, Lq )
      class( Space_t ), intent( inout ) :: this
      real, intent(in) :: tol_max, tol_min
      character(len=20) :: adapt_space
      integer :: max_adapt_level
      real, intent(in), optional :: Lq

      if ( adapt_space == 'none' .or. adapt_space == '-') then
         this%adaptation = .false.
         allocate( Adaptation_t :: this%adapt )
         call this%adapt%init( tol_max, tol_min, adapt_space, max_adapt_level )
      else
         this%adaptation = .true.
         allocate( Adaptation_t :: this%adapt )
         !print*, 'allocate Adapt type!'

         if ( present(Lq) ) then
            call this%adapt%init( tol_max, tol_min, adapt_space, max_adapt_level, Lq )
         else
            call this%adapt%init( tol_max, tol_min, adapt_space, max_adapt_level )
         endif
      endif

   end subroutine allocAdaptation

!   !> allocation of various types of this%estim, it depends on the estim_space technique
!   subroutine allocEstims( this, estim_space, Lq)
!      class(Space_t), intent(inout) :: this
!      character(len=20), intent(in) :: estim_space
!      real, intent(in) :: Lq
!
!      select case (estim_space)
!         case('DWR')
!            allocate( DWR_estim_t :: this%estim )
!            call this%estim%init( Lq )
!         case default
!         print*, 'Allocation of state%space%estim not done for other estim techniques than DWR.'
!      end select
!
!   end subroutine allocEstims

   subroutine createQuadRules( this )
      class (Space_t ), intent(inout) :: this
      integer :: i, j, j1

      !set in paramets
      !QnumOffset = maxVrule   !  for the old assessment

      !VOLUME RULE
      allocate ( this%V_rule( 1:maxVrule + maxGrule ) )
      this%V_rule(:)%def = .false.


      ! preparing of volume integration quadratures
      !print*,'ATTENTION in createQuadRules in o_space.f90, tested VERSION'
      ! SEE 3 lines bellow

      do i=1,maxVrule

      !do j = 0, MaxDegreeImplemented
      !    i = this%Qdeg(j, 1)

         if(.not. this%V_rule(i)%def) then
            call this%V_rule(i)%createVrule(i)
            this%V_rule(i)%Qnum = i
         endif

         !if(minval(this%V_rule(i)%weights(:) ) < 0. ) print *,'???? w',i
         !if(minval(this%V_rule(i)%lambda(:,:)) < 0. ) print *,'???? l',i
      enddo


      ! testing of the accuracy of quadrature rules on monomials
      ! call this%test_Vrules()
      ! stop 'Vrules tested! in o_space.f90'
      ! do j= 1,maxVrule + maxGrule
      !   print*, 'V-rule of deg', j, ' defined:', this%V_rule(j)%def
      ! enddos


   !GAUSS RULE
   ! preparing of edge integration quadratures
   allocate( Gauss_rule :: this%G_rule(1:maxGrule) )
   do i=1,maxGrule
      call this%G_rule(i)%createGaussRule(i)
      call this%G_rule(i)%initLegendrePolynomials( )
   end do



!     LEGENDRE polynomials - TEST
!     i = 12
!     do j=1, this%G_rule(i)%Qdof
!      print*, ':', this%G_rule(i)%lambda(j), this%G_rule(i)%Leg_phi(0:MaxDegreeImplemented, j)
!     enddo
!
!     do j=0,MaxDegreeImplemented
!        do j1=0,MaxDegreeImplemented
!           print*,'EDS',j,j1, dot_product(this%G_rule(i)%weights(:), &
!                this%G_rule(i)%Leg_phi(j, 1:) * this%G_rule(i)%Leg_phi(j1, 1:) )
!        enddo
!        print*, this%G_rule(i)%Leg_phi(j,0)
!     enddo
!     stop


    ! Bi-Gauss volume integration for quadrilaterals, only weights
    do i=1, maxGrule
       call this%V_rule(i+maxVrule)%create4V_rule(this%G_rule(i), i)
       this%V_rule(i+maxVrule)%Qnum = i + maxVrule
    enddo

    ! preparing of volume integration quadratures
    ! has to be < = QnumOffset, >QnumOffset are for quadrilateralls !!!

    allocate( Lagrang_rule :: this%L_rule( 0:maxLrule + maxVrule ) )

    if (maxLrule >= maxVrule) &
      stop 'Problem in createQuadRules - maxLrule >= maxVrule'

    do i=0,  maxLrule
       call this%L_rule(i)%createLagrangeRule(i)
       call this%L_rule(i + maxVrule )%createLagrangeRule(i + maxVrule)
    enddo

   end subroutine createQuadRules

   !> set space quadrature rules and basis functions
   subroutine initSpace( this )
      class (Space_t ), intent(inout) :: this
      integer :: i, j

      !FERROR TODO setting deg(1:ndim)

      ! array for the computation of the hierarchical projection error
      allocate ( this%ldeg(-1: MaxDegreeImplemented ) )

      do i=-1, MaxDegreeImplemented
       this%ldeg(i) =  DOFtriang(i)
      enddo

      !  this%ldeg( -1 : MaxDegreeImplemented ) = DOFtriang( (/ (i, i = -1, MaxDegreeImplemented) /) )

      ! for given degree of polynomial approximation gives the degree of Dunavant/Wandzura quadrature
      allocate(this%Qdeg(0:MaxDegreeImplemented, 1:2), source = 0 )


      if(MaxDegreeImplemented > 10) then
         print*,'Adding of Dunavant/Wandzura quadrature necessary in problem.f90'
         stop
      endif

    ! degree of approx: elem%deg: 0   1   2   3   4   5   6   7   8   9  10
      this%Qdeg(0:10, 1)   =   (/ 2,  5,  8, 12, 14, 17, 21, 22, 23, 23, 23 /)

      call this%createQuadRules()

!      this%max_HGlevel = 2

      if (this%adapt%adapt_type == 'HG') then
         this%adapt%max_HGlevel = 1
      else
         this%adapt%max_HGlevel = 0
      endif

      this%adapt%HG = 2**( this%adapt%max_HGlevel+1)-1

      call this%ComputeGSCoefficients()

      call this%ComputeBasisFunctions()

      !call this%V_rule(14)%printVrule( 'V_rule14' )

      this%gridS_allocated = .false.

      this%deg_plus = .false.
      ! how much to increase the polynomial degree
      this%plusDeg = 1
!          ! numerical Gram-Schmidt, general deg
!    call ComputeGramSchmidtCoefficients( )
!    call SetBasisFunctions(PHI_orthonormal)

   end subroutine initSpace



!> evaluation of the coefficients of the transformation to the orthonormal basis
  !> by GramSchmidt
  subroutine ComputeGSCoefficients( this )
    class( Space_t ), intent( inout ) :: this
    real, dimension(:,:), allocatable :: phi
    !type(volume_rule), pointer :: V_rule
    real, dimension(1:nbDim) :: xi, xc
    integer:: Qnum, Qdof, dof
    integer:: i, k, len, l

    dof = (MaxDegreeImplemented + 1) * (MaxDegreeImplemented + 2) / 2

    allocate( this%GScoeff(3:4, 1:dof, 1:dof) )

    do len = 3, 4  ! triangles or quadrilaterals
       if(len == 3) then
          Qnum = maxVrule
          xc(1:nbDim) = 1./3
       else
          Qnum = maxVrule + maxGrule
          xc(1:nbDim) = 1./2
       endif

       ! setting of test functions in integ nodes
       !V_rule => this%V_rule(Qnum)
       associate ( V_rule => this%V_rule(Qnum) )

       Qdof = V_rule%Qdof
       allocate(phi(1:dof,1:Qdof) )

       do l=1, Qdof ! values in l-th integ. node
          xi(1:nbDim) = V_rule%lambda(l,1:nbDim)
          call phi_vals(dof, xi-xc, phi(1:dof, l) )

          !write(90+i,'(25es12.4)')xi(1:nbDim)-xc(1:nbDim), phi(1:nbDim1, l)
       enddo

       ! Gram-Schmidt
       do k=1,dof
          do i=1, k-1
             this%GScoeff(len, i, k) = dot_product(V_rule%weights(1:Qdof),  &
                  phi(k, 1:Qdof) * phi(i, 1:Qdof) ) * 0.5
          enddo

          do i=1, k-1
             phi(k, 1:Qdof) = phi(k, 1:Qdof) - phi(i, 1:Qdof) * this%GScoeff(len, i, k)
          enddo

          this%GScoeff(len, k, k) = dot_product(V_rule%weights(1:Qdof),  &
               phi(k, 1:Qdof) * phi(k, 1:Qdof) ) * 0.5

          if( this%GScoeff(len, k, k) <= 0 ) then
             print*, 'Trouble in  ComputeGramSchmidtCoefficients',dof, Qnum, Qdof
             stop
          endif

          this%GScoeff(len, k, k) = ( this%GScoeff(len, k, k))**0.5

          phi(k, 1:Qdof) =  phi(k, 1:Qdof) / this%GScoeff(len, k, k)
       enddo

       !do l=1, Qdof ! values in l-th integ. node
       !   xi(1:nbDim) = V_rule%lambda(l,1:nbDim)
       !   write(90+len,'(25es12.4)') xi(1:nbDim), phi(1:nbDim1, l)
       !enddo

       deallocate(phi)

       end associate ! V_rule

    enddo
  end subroutine ComputeGSCoefficients



  !> setting of basis functions, evaluation of test functions and their derivatives
  !> on reference elements \f$\hat{K} \f$ in all integ node (volumes, edges)
  !>
  !> \f$ \hat{\varphi}_{d}(\hat{x}_{iq, l}),\
  !> \hat{\partial}_j \hat{\varphi}_{d}(\hat{x}_{iq, l} )\f$,
  !> \f$\hat{\varphi}_{d} \in \f$
  !> hierachical basis, \f$ \hat{x}_{iq, l}\f$: \f$l^{\rm th}\f$  node
  !> of \f$iq^{\rm th}\f$ quadrature rule,
  !> volume rules for triangles and quadrilateralls, edge rules.
  !> For hanging nodes (HG) each edge can be split in 2, 4, 8, ... equidistant pieces
  subroutine ComputeBasisFunctions(this)
    class(Space_t), intent (inout), target :: this
    type(volume_rule), pointer :: V_rule
    integer :: iq, iK, l, ie, len, i, j, ii
    integer :: max_deg, dof, Qdeg, Qdof
    real ::  t
    real, dimension (1:nbDim) :: x0, a0
    real, dimension(:,:), allocatable :: p0
    real, dimension(:,:,:), allocatable :: void
    logical :: V_rule_alloc



    dof = DOFtriang( MaxDegreeImplemented )

    do iq=1, maxVrule+maxGrule
       if( iq <= maxVrule) then
          ! triangular volume quadratures, only some of them are used
          len = 3
          V_rule_alloc = .false.
          do j=0, MaxDegreeImplemented
             if( this%Qdeg(j, 1) == iq) V_rule_alloc = .true.
          enddo

          !if(iq ==1) print*,'ATTENTION in createQuadRules in o_space.f90, tested VERSION (2)'
          V_rule_alloc = .true.

          if(.not. V_rule_alloc) goto 100

       elseif(iq >= maxVrule + 1) then
          ! quadrilateral  volume quadratures
          len = 4
       else
          ! NO quadratures
          goto 100
       endif

       V_rule => this%V_rule(iq)
       Qdeg = V_rule%Qdeg
       Qdof = V_rule%Qdof

       allocate(V_rule%phi(1:dof,1:Qdof) )
       allocate(V_rule%Dphi(1:dof,1:nbDim,1:Qdof) )

       allocate(p0(1:Qdof, 1:nbDim) )
       p0(1:Qdof, 1:nbDim) = V_rule%lambda(1:Qdof,1:nbDim)!values in integ nodes


       call PHIorthonormalNew(Qdof, p0, len, dof, this%GScoeff(len,1:dof,1:dof), &
       V_rule%phi(1:dof, 1:Qdof), V_rule%Dphi(1:dof,1:nbDim,1:Qdof) )

       ! print*,'plotting of the values of test functions in integ nodes'
       !if(iq == 15) then
       !   do i=1,dof
       !      do j=1,Qdof
       !         write(20+i, *) p0(j, 1:2), V_rule%phi(i, j)
       !      enddo
       !   enddo
       !endif

       deallocate(p0)

100    continue
    enddo ! iq

       !   print*, 'Max HGlevel: ',this%max_HGlevel, this%HG

    ! EDGE QUADRATURES
    do iq=1, maxGrule
       Qdeg = this%G_rule(iq)%Qdeg
       Qdof = this%G_rule(iq)%Qdof

       allocate(this%G_rule(iq)%phi(3:4,1:4,this%adapt%HG, 1:dof,1:Qdof) )
       allocate(this%G_rule(iq)%Dphi(3:4,1:4,this%adapt%HG, 1:dof, 1:nbDim,1:Qdof) )

       do iK = 3,4   ! triangle or quadrilateral
          do ie = 1, iK   ! loop through edges
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

             ii = 0

             do i=0,this%adapt%max_HGlevel
                do j=1, 2**i
                   ii = ii + 1

                   allocate(p0(1:Qdof, 1:nbDim) )
                   do l=1, Qdof  !loop though 1D Gauss quadrature rule
                      t = this%G_rule(iq)%lambda(l)

                      p0(l, 1:nbDim) = x0(1:nbDim) +  a0(1:nbDim) * (j-1 +t)/(2**i)
                   enddo ! l

                   call PHIorthonormalNew(Qdof, p0, iK, dof, this%GScoeff(iK,1:dof,1:dof), &
                        this%G_rule(iq)%phi( iK, ie, ii, 1:dof, 1:Qdof), &
                        this%G_rule(iq)%Dphi(iK, ie, ii, 1:dof, 1:nbDim, 1:Qdof) )

                   deallocate(p0)
                enddo
             enddo

          enddo ! ie
       enddo  ! iK
    enddo ! iq

!    do l=0, 3
!       write(*,'(10es12.4)')  this%G_rule(2)%phi(3,1,l+4,2:3,1)
!    enddo
!    write(*,'(a4,10e10.3)') '??>.', this%G_rule(2)%phi(3,2,1:10,2)


     !LAGRANGE  QUADRATURES
    do iq=0, maxLrule  + maxVrule
       if( iq <= maxLrule) then
          ! triangular volume quadratures
          len = 3
       elseif(iq >= maxVrule) then
          ! quadrilateral  volume quadratures
          len = 4
       else
          ! NO quadratures
          goto 200
       endif

       Qdeg = this%L_rule(iq)%Qdeg
       Qdof = this%L_rule(iq)%Qdof

       allocate(this%L_rule(iq)%phi(1:dof,1:Qdof) )
       allocate(this%L_rule(iq)%Dphi(1:dof,1:nbDim,1:Qdof) )

       allocate(p0(1:Qdof, 1:nbDim) )

       p0(1:Qdof, 1:nbDim) =  this%L_rule(iq)%lambda(1:Qdof,1:nbDim)  ! values in integ node

       call PHIorthonormalNew(Qdof, p0(1:Qdof, 1:nbDim), len, dof, this%GScoeff(len, 1:dof, 1:dof), &
            this%L_rule(iq)%phi(1:dof, 1:Qdof), this%L_rule(iq)%Dphi(1:dof,1:nbDim,1:Qdof) )


       deallocate(p0)
     !  deallocate(void)

200    continue


    enddo ! iq

  end subroutine ComputeBasisFunctions

  !> test the accuracy of quadrature rules on monomials
 subroutine test_Vrules( this )
   class( Space_t ), intent(in) :: this
   integer :: i, i1, l
   integer :: a,b, deg
   real :: area, coef, quad, x, y, val, exact, err
   logical :: line


   ! negative weights or nodes outsides of element

   ! All quadratures
   !do i = 1,maxVrule

   ! ONLY the USED ONES
   do i1 = 0, MaxDegreeImplemented
      i = this%Qdeg(i1, 1)

      do l = 1, this%V_rule(i)%Qdof

         if(minval(this%V_rule(i)%lambda(l,:)) < 0. ) &
              write(*,'(a6,2i5,3es12.4)') 'negQUA',i,l,this%V_rule(i)%lambda(l,:)
         if(maxval(this%V_rule(i)%lambda(l,:)) > 1. ) &
              write(*,'(a6,2i5,3es12.4)') 'negQUA',i,l,this%V_rule(i)%lambda(l,:)
         if(this%V_rule(i)%weights(l) < 0. ) &
              write(*,'(a6,2i5,3es12.4)') 'negQUA',i,l,this%V_rule(i)%weights(l)
      enddo
   enddo
   print*,'######## positivity of quadratures checked'


   area = 0.5D+00

   !do deg = 34, 36
   do deg = 0, 36
      do b = 0, deg

         a = deg - b
         write ( *, '(a6,i3,a1,i2,a1,i2,a28)') &
              'deg=', a+b,'(',a,'+',b,') --------------------------------'
         !
         !  Multiplying X**A * Y**B by COEF will give us an integrand
         !  whose integral is exactly 1.  This makes the error calculations easy.
         !
         coef = real ( a + b + 2, kind = 8 ) * real ( a + b + 1, kind = 8 )
         do i = 1, b
            coef = coef * real ( a + i, kind = 8 ) / real ( i, kind = 8 )
         end do

         line = .false.

         ! All quadratures
         do i = 1,maxVrule

         ! ONLY the USED ONES
         !do i1 = 1, MaxDegreeImplemented
         !   i = this%Qdeg(i1, 1)

            quad = 0.0D+00

            do l = 1,this%V_rule(i)%Qdof

               x = this%V_rule(i)%lambda(l,1)
               y = this%V_rule(i)%lambda(l,2)

               if ( a == 0 .and. b == 0 ) then
                  val = coef
               else if ( a == 0 .and. b /= 0 ) then
                  val = coef        * y**b
               else if ( a /= 0 .and. b == 0 ) then
                  val = coef * x**a
               else if ( a /= 0 .and. b /= 0 ) then
                  val = coef * x**a * y**b
               end if

               quad = quad + this%V_rule(i)%weights(l) * val

            end do

            quad = area * quad

            exact = 1.0D+00
            err = abs ( exact - quad )

            !if(i >= 20) then
               if(abs(err) < 1E-10) then
                  write ( *, '(a6,i3,a1,i2,a1,i2,2(a8,i8),2x,es18.10,2x,es14.6)' ) &
                       'deg=', a+b,'(',a,'+',b,') Qnum=', i, 'Qdof=', this%V_rule(i)%Qdof, quad, err
                  line = .true.
               else
                  write ( *, '(a6,i3,a1,i2,a1,i2,2(a8,i8),2x,es18.10,2x,es14.6, a6)' ) &
                       'deg=', a+b,'(',a,'+',b,') Qnum=', i,  'Qdof=', this%V_rule(i)%Qdof,quad, err,'!!!!!!'
               endif
            !endif



            !deallocate ( wtab )
            !deallocate ( xytab )

      end do
      !if(line) write ( *,*)
    end do

  end do


 end subroutine test_Vrules

 !> copies geometry datas from mesh to state%space
 subroutine copyGeometry( this, h, domain_volume )
   class( Space_t ) , intent( inout) :: this
   real, intent(in) :: h
   real, intent(in) :: domain_volume

   write(debug,*) 'Temporarily h, domain_volume are duplicated in mesh and %space. '

   this%h = h
   this%domain_volume = domain_volume

 end subroutine copyGeometry


 !> test for recontructions
 subroutine Reconstruction_test(this )
   class( Space_t ), intent( inout ), target :: this
   type(Gauss_rule), pointer :: G_rule, G_ruleE
   real, dimension(:,:), allocatable :: phiE, phi
   real, dimension(:,:), allocatable :: xi
   real, dimension(:,:,:), allocatable :: wi, wE, ww
   real, dimension(:,:), allocatable :: A, AA
   real, dimension(:), allocatable :: b, Ab, x

   real :: h, t,val
   integer :: N, i,j,k,l, p_K, Qnum, Qdof, QdofE, dofE, nr, nc, Rp_K

   h = 1.0  !0.1
   N = 3

   p_K = 2
   Qnum = 3

   Rp_K = 3



   G_rule  => this%G_rule(Qnum)
   G_ruleE => this%G_rule(10)

   Qdof = G_rule%Qdof
   QdofE = G_ruleE%Qdof
   dofE = G_ruleE%Qdeg

   ! basis functions of the given quadrature
   allocate( phi(0:MaxDegreeImplemented, 1:Qdof) )
   phi(0:MaxDegreeImplemented, 1:Qdof) = G_rule%Leg_phi(0:MaxDegreeImplemented, 1:Qdof)

   ! basis functions of the "exact" quadrature
   allocate( phiE(0:MaxDegreeImplemented, 1:QdofE) )
   phiE(0:MaxDegreeImplemented, 1:QdofE) = G_ruleE%Leg_phi(0:MaxDegreeImplemented, 1:QdofE)


   allocate( xi(1:N, 1:5) )

   do i=1,N
      xi(i, 1) = (i-1)*h
      xi(i, 2) = (i  )*h
      xi(i, 3) = xi(i,2) - xi(i,1)
      xi(i, 4) = (xi(i,2) - xi(i,1)) / 2
   enddo

   ! setting of the function in integ nodes
   allocate( wi(1:N, 1:Qnum, 1:5) )
   do i=1,N
      do j=1,Qnum
         t = xi(i,1) + h * G_rule%lambda(j)
         wi(i,j, 1) = F_proj(t)

         write(*,'(a6,2i5,8es12.4)') ' FUNC:',i,j,t,wi(i,j,1)
         write(97, *) t,wi(i,j,1)
      enddo
   enddo

   ! L2 projection
   ! setting of the function in integ nodes
   allocate( wE(1:N, 1:QdofE, 1:5) )
   allocate( ww(1:N, 0:p_K, 1:5) )
   print*
   do i=1,N
      do j=1,QdofE
         t = xi(i,1) + h * G_ruleE%lambda(j)
         wE(i,j, 1) = F_proj(t)

         !write(*,'(a6,2i5,8es12.4)') ' FUNC:',i,j,t,wE(i,j,1),phiE(0:4, j)
         write(99, *) t,wE(i,j,1),phiE(0:4, j)
      enddo
      print*

      ! projection in basis coeficients
      do j=0,p_K
         do l = 0, p_K
            val = dot_product(G_ruleE%weights(1:QdofE), phiE(j, 1:QdofE)*phiE(l,1:QdofE))
            print*,'Legendre',j,l,val
         enddo

         ! normalization
         val = dot_product(G_ruleE%weights(1:QdofE), phiE(j, 1:QdofE)*phiE(j,1:QdofE))

         ww(i,j,1) = dot_product(G_ruleE%weights(1:QdofE), phiE(j, 1:QdofE)*wE(i,1:QdofE, 1)) / val
      enddo

      write(*,'(a6,i5,8es12.4)') ' ProB:',i, ww(i,:,1)

      ! projection in integ nodes
      wi(i,1:Qdof,2) = matmul(ww(i,0:p_K, 1),  phi(0:p_K, 1:Qdof) )

      do j=1,Qdof
         t = xi(i,1) + h * G_rule%lambda(j)
         write(*,'(a6,2i5,8es12.4)') ' ProN:',i,j,t,wi(i,j,2)
         write(96, *) t,wi(i,j,2)
      enddo


      ! projection in many nodes
      do j = 1, QdofE
         t = xi(i,1) + h * G_ruleE%lambda(j)
         val =  dot_product(ww(i,0:p_K, 1),  phiE(0:p_K, j) )
         write(88,*) t, val
      enddo


   enddo

   ! reconstruction
   ! A number of nodes is the number of rows of A, number of DOF is the number of columns of A
   print*,'-------------------------------------------'
   nr = N*Qdof
   nc = Rp_K + 1

   print*,'nr = ', nr,',   nc = ', nc
   if(nr <  nc) stop 'nr (number of equations) is lower then nc (number of unknowns'
   allocate(A(1:nr, 1:nc) , AA(1:nc, 1:nc), b(1:nr), Ab(1:nc) , x(1:nc) )

   do i=1, N
      do j = 1, Qdof
         l = (i-1)*Qdof + j  ! index of the node
         t = xi(i,1) + h * G_rule%lambda(j)  ! x-variable of the node

         do k=1, nc   ! number of global basis over the whole stencil
            A(l, k) = t**(k-1)   ! canonical basis
         enddo

         b(l) = wi(i, j, 2) ! RHS, function, which has to be reconstructed, in integ nodes

         write(*,'(a6,3i5, 20es12.4)') 'A, b:',i,j,l, b(l), A(l,:)
      enddo
   enddo

   AA(1:nc, 1:nc) = matmul(transpose(A(1:nr, 1:nc)), A(1:nr, 1:nc) )
   Ab(1:nc) = matmul(transpose(A(1:nr, 1:nc)), b(1:nr) )

   print*,'                  TRANSPOSE  '
   do l=1,nc
      write(*,'(a6,i5, 20es12.4)') 'AA, b:',l, Ab(l), AA(l,:)
   enddo

   call MblockInverse(nc, AA)
   x(1:nc) = matmul(AA(1:nc, 1:nc), Ab(1:nc) )

   do i=1, N
      do j = 0, 10
         t = xi(i,1) + 1.*j/10 * (xi(i,2) - xi(i,1) )
         val = 0.
         do k=1,nc
            val = val + x(k) * t**(k-1)
         enddo
         write(98,*) t, val
      enddo
   enddo


   deallocate(xi,wi,wE, ww, phi, phiE, A, AA, b, Ab, x)

 end subroutine Reconstruction_test

 function F_proj(x)
   real :: F_proj
   real, intent(in) :: x

   !F_proj = x*x*x + 0.1*x*x - x
   F_proj = sin(1.5*x)

 end function F_proj



end module space_mod
