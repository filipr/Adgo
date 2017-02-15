!>subroutines for marking elements for hp adaptation
module marking

  use data_mod  ! contains type(mesh) :: grid for computation
  use element_mod
  use main_data
!  use euler_problem
!  use apost_estimation
!  use project_estimation

  implicit none

  public:: MarkTopElements
  public:: MarkElements
  public:: RedGreenMarking
  public:: HGMarking
contains

  !> mark elements for hp-adaptation based on the error estimates,
  !> subroutines marks only the elements with the highest error
  !> uses the arrays elem%estim_loc to refinement
  subroutine MarkTopElements( )
    class(element), pointer :: elem, elemD, elem1
    real :: estim_max, estim_min, tol_max, estim_tot
    real :: limit_min, limit_max
    integer :: i, j, lev
    real , dimension(:), allocatable :: est
    integer , dimension(:), allocatable :: iest
    logical :: only_h_adapt, only_p_adapt
    real :: ratio_min, ratio_max

    only_h_adapt = .false.
    only_p_adapt = .false.
    if(state%space%adapt%adapt_space == 'HGh') only_h_adapt = .true. ! only h-refinement
    if(state%space%adapt%adapt_space == 'HGp') only_p_adapt = .true. ! only p-refinement


    ! setting of the appropriate element error estimates
    do i=1,grid%nelem
       elem => grid%elem(i)
       elem%hsplit = 0
       elem%psplit = 0

       ! quantity for error control, depends on the method and the problem
       if( state%space%estim_space == 'pNeu')  then
          elem%estim_loc = elem%eta(resST, 1)

       elseif( state%space%estim_space == 'DWR')  then
!          elem%estim_loc = elem%eta(resST, 1)

       else
          if(state%time_dependent) then
             if( state%time%disc_time /= 'STDG') then
                elem%estim_loc = (elem%estim_loc + elem%jumpsJh)**0.5     !  [ VD, MATCOM ??]
             else
                elem%estim_loc = elem%eta(resST, 1)
             endif
          else
             if( state%time%disc_time /= 'STDG') then
                !elem%estim_loc = (elem%estim_loc + elem%jumpsJh)**0.5     !  [ VD, MATCOM ??]
                elem%estim_loc = (elem%estim_loc )**0.5     !  elem%jumpsJh is incorect, Maybe
             else
                elem%estim_loc = elem%eta(resS, 1) / state%time%tau(1)
             endif
          endif
       endif

    enddo


    ! setting of the global values
    if(state%time_dependent) then
       tol_max = state%space%adapt%tol_max * sqrt(state%time%tau(1)/ state%time%FinTime)
    else
       tol_max = state%space%adapt%tol_max
    endif

    estim_max = maxval(grid%elem(:)%estim_loc)
    estim_min = minval(grid%elem(:)%estim_loc)
    estim_tot = sqrt(sum(grid%elem(:)%estim_loc**2))

    write(*,'(a40, 6es12.4)')'estim_max, estim_min, estim_tot, tol_max:', &
         estim_max, estim_min, estim_tot, tol_max

    ! ordering of the elements according the values of estimator
    allocate(iest(1:grid%nelem), est(1:grid%nelem) )
    do i=1,grid%nelem
       iest(i) = i
       est(i) = grid%elem(i)%estim_loc
    enddo

    call order_estims(grid%nelem, iest, est)


    ! variant with elements with error at lest 0.1 % of the errors
    !limit_max = 0.5 * estim_max
    !limit_min =   5 * estim_min

    !ratio_max = 1.0
    !ratio_min = 0.0


    ! variant with 10% of elemets ordered with inceasing error
    limit_max = 0.
    !limit_min = 5 * estim_min
    limit_min = 5 * estim_max   ! NO COARSENING

!    ratio_max = 0.1    ! first 10% of elements refined

    ! FR CHANGED
!    print*, 'All 100% of the elements are refined in MarkTopElements'
!    ratio_max = 1.0   ! first 100% of elements refined
    ratio_max = 0.2   ! first 20% of elements refined
    ratio_min = 1 + 0.95  !0.9    ! last  10% of elements de-refined



    do j=1, grid%nelem !!  int(0.1 * grid%nelem)
       i = iest(j)
       elem => grid%elem(i)


       ! refinement
       if(elem%estim_loc >= limit_max .and. j <= int(ratio_max * grid%nelem) ) then
          if(only_h_adapt) then
             elem%hsplit = 4

          elseif(only_p_adapt) then
             elem%psplit = 1

          else
             if(elem%reg <= elem%regT0 .and. elem%deg < MaxDegreeImplemented-1 ) then

                elem%psplit = 1
                !write(100 + state%space%adapt%adapt_level + 1,*) elem%xc(:), grid%x(elem%face(idx, 1:3), 1:2)
                !write(*,'(a5,6es12.4)') 'p+', elem%xc(:), elem%estim_loc, elem%reg, elem%regT0

             else
                elem%hsplit = 4
                !write(200 + state%space%adapt%adapt_level + 1,*) elem%xc(:), grid%x(elem%face(idx, 1:3), 1:2)
                !write(*,'(a5,6es12.4)') 'h+', elem%xc(:), elem%estim_loc, elem%reg, elem%regT0

             endif
          endif

       ! de-refinement
       elseif(elem%estim_loc < limit_min .and. j >= int(ratio_min * grid%nelem) ) then

          if(only_h_adapt) then
             elem%hsplit = -4

          elseif(only_p_adapt) then
             elem%psplit = -1

          else
             if(elem%reg <= elem%regT0) then
                elem%hsplit = -1
                write(300 + state%space%adapt%adapt_level + 1,*) elem%xc(:), grid%x(elem%face(idx, 1:3), 1:2)

             else
                if(elem%deg > MinDegreeImplemented)  elem%hsplit = -1  ! elem%psplit = -1
                write(400 + state%space%adapt%adapt_level + 1,*) elem%xc(:), grid%x(elem%face(idx, 1:3), 1:2)
                !write(*,'(a5,6es12.4)') 'p-', elem%xc(:), elem%estim_loc, elem%reg, elem%regT0

             endif
          endif

       endif

       !write(100 + 10*state%space%adapt%adapt_level + 5 + elem%hsplit/4,*) elem%xc(:)
    enddo


    deallocate(iest, est)

    !WE ENFORCE THE H-DEREFINEMENT
    !if(state%space%adapt%adapt_level == 3) then
    !   grid%elem(:)%psplit = 0
    !   grid%elem(:)%hsplit = -1
    !endif


    ! checking of the marked elements for h-derefinement
    if(ratio_min < 1. .or. limit_min <  estim_max ) then
       do i = 1, grid%nelem
          elem => grid%elem(i)
          lev = elem%RGhistory

          ! only elements which were refined in previous computations and which are  direct
          ! daughters of their mother elements will be checked
          if(lev > 0 ) then
             if(elem%i == elem%RG(lev)%daughter(4)) then

                ! not necessary elemD%LAST is elem
                !if(elem%hsplit /= -1 ) goto 200 ! elem is not marked for DE-refinement

                ! we go though all daughter elements including the original one
                do j=1,elem%RG(lev)%subel
                   !write(*,'(a8,5i5)') 'elmD:',i, elem%i, lev, elem%RG(lev)%daughter(4)

                   elemD => grid%elem(elem%RG(lev)%daughter(j))
                   !write(*,'(a6, 14i5)') 'h-DREF',elem%i, elem%hsplit, &
                   !     j, elemD%i, elem%RGhistory, elemD%RGhistory, elemD%hsplit

                   if(elem%RGlevel /= elemD%RGlevel) goto 200   ! different level of RG refinement
                   if(elemD%hsplit /= -1 ) goto 200       ! this sub-element is not marked for refinement

                   if(elemD%per > 0) then
                      goto 200 ! DEREFINEMENT NOT IMPLEMENTED FOR PERIODIC
                      !elem1 => grid%elem(elemD%per)
                      !if(elem1%eta(resST, 1) > tol_min * scale  &
                      !  .or. elem%RGlevel /= elem1%RGlevel)
                   endif



                   ! NO derefinement
                enddo

                ! the error estimate of all daughter elements are smaller than tol_min,
                ! we mark them for derefinement

                grid%elem(elem%RG(lev)%daughter(1:elem%RG(lev)%subel))%hsplit = - 4

                write(500 + state%space%adapt%adapt_level + 1,*) elem%xc(:), grid%x(elem%face(idx, 1:3), 1:2)
                !write(*,'(a5,6es12.4)') 'h-', elem%xc(:), elem%estim_loc, elem%reg, elem%regT0

                !!goto 100
200             continue

                if (only_h_adapt) goto 100 ! NO h -> p substitution

                ! h -> p substitution
                ! not necessary elemD%LAST is elem
                !if( elem%hsplit /= 0 .and. elem%psplit /= 0  ) goto 100 ! NO h-p substitution
                !if( elem%reg > elem%regT0 ) goto 100 ! NO h-p substitution

                ! we go though all daughter elements including the original one
                do j=1,elem%RG(lev)%subel
                   elemD => grid%elem(elem%RG(lev)%daughter(j))
                   if(elemD%deg ==  MaxDegreeImplemented-1 .or. &
                        elemD%hsplit /= 0 .or. elemD%psplit /= 0 .or. &!B elemD%reg > 1. .or. &
                        elemD%reg > elemD%regT0 .or. elem%RGlevel /= elemD%RGlevel) goto 100
                   ! NO h -> p  -substitution
                enddo
                ! all daughter elements are not marker for p- nor h- adaptation,
                ! solution is regular, we derefine them and increase p

                print*,'?? h -> p ?',lev, elem%i, elemD%i

                grid%elem(elem%RG(lev)%daughter(1:elem%RG(lev)%subel))%hsplit = - 4
                grid%elem(elem%RG(lev)%daughter(1:elem%RG(lev)%subel))%psplit = + 1

                ! !write(600 +state%space%adapt%adapt_level + 1,*) elem%xc(:)

100             continue

             endif
          endif
       enddo
    end if !end of derefinement


    ! removing the marks of marked elements which can not be recoarsened
    do i = 1, grid%nelem
       elem => grid%elem(i)
       if(elem%hsplit == -1) elem%hsplit = 0
       if(elem%hsplit == -4) &
            write(700 + state%space%adapt%adapt_level + 1,*) elem%xc(:)
    enddo


149 continue


    if(state%space%adapt%adapt_type == 'HG') then
       call HGMarking()

    elseif(state%space%adapt%adapt_type == 'RG') then
       call RedGreenMarking( )

    else
       print*,'UNKNOWN marking strategy in marking.f90'
    endif

    j= 0
    do i = 1, grid%nelem
       elem => grid%elem(i)
       if(elem%hsplit == -4) then
          write(700 + state%space%adapt%adapt_level + 1,*) elem%xc(:)
          j = j+ 1
       endif

    enddo
    !if( j > 0) print*,'SOME ELEMENTS MARKED for DE-REFINEMENT !!!!!'

    print*,'$$$ - p-refinement', sum(grid%elem(:)%psplit, mask = grid%elem(:)%psplit == 1)
    print*,'$$$ - h-refinement', sum(grid%elem(:)%hsplit, mask = grid%elem(:)%hsplit == 4) / 4

    !stop '333e3e43e de43adapt ecde'

  end subroutine MarkTopElements

  !> order of elements according the values of estimators
  subroutine order_estims(nelem, iest, est)
    integer, intent(in) :: nelem
    integer, dimension(1:nelem), intent(inout) :: iest
    real, dimension(1:nelem), intent(inout) :: est
    integer :: ival,j,k, ipoc
    real :: val

    do j=1, nelem
       ipoc = 0
       do k=1, nelem-1
          if(est(k) < est(k+1) ) then
             val = est(k)
             est(k) = est(k+1)
             est(k+1) = val

             ival = iest(k)
             iest(k) = iest(k+1)
             iest(k+1) = ival

             ipoc = 1
          endif
       enddo
       if(ipoc == 0) goto 10
    enddo
    10 continue
    print*,'ordered after j = ', j

  end subroutine order_estims


  !> mark elements for hp-adaptation based on the error estimates
  !> uses the arrays elem%estim_loc to refine, but also other arrays for p-refinement
  subroutine MarkElements( )
    class(element), pointer :: elem, elemD, elem1
    integer :: lev
    integer :: i,j,k, ib, pK
    real :: tol_max, tol_min, dist_ref, dist_ref1
    real :: rez_estim, rez2
    real :: scale, s, ss, s1, s2, weig
    integer:: ityp, ipoc, ipc
    logical :: only_h_adapt, only_p_adapt


    only_h_adapt = .false.
    only_p_adapt = .false.
    if(state%space%adapt%adapt_space == 'HGh') only_h_adapt = .true. ! only h-refinement
    if(state%space%adapt%adapt_space == 'HGp') only_p_adapt = .true. ! only p-refinement

    if(state%space%adapt%adapt_space == 'RGh') only_h_adapt = .true. ! only h-refinement
    if(state%space%adapt%adapt_space == 'RGp') only_p_adapt = .true. ! only p-refinement

!    if (state%space%adapt%adapt_method == 'ALG' .or. state%space%adapt%adapt_method == 'ALG2' .or. state%type_IC == 8 ) then
!       print*,'CKECK in marking.f90'
!       only_h_adapt = .true.
!
!    endif

    !dist_ref = 1. * 1.5**(-state%space%adapt%adapt_level)
    !dist_ref = 1.2 * 2.**(-state%space%adapt%adapt_level)
    !dist_ref = 0.8 * 3.**(-state%space%adapt%adapt_level)  !A4
    dist_ref = 0.8 * 2.9**(-state%space%adapt%adapt_level)  !A5
    dist_ref = 1.0 * 3.0**(-state%space%adapt%adapt_level)
    !dist_ref = 0.75 * 4.**(-state%space%adapt%adapt_level)

    dist_ref1 = dist_ref * 2.2

    lev = 0

    !write(500,*) 'mark------------------------------------------'

    ! for singul, ityp=2 is much better, for linBL work both
    ! (ityp==1) scale = 1, (ityp==2) scale = 1./N, (ityp==3) scale = |K|/|Omega|
    !ityp = 1
    !ityp = 2
    ityp = 3
    !ityp = 4

    scale = 1.
    if(ityp == 2) scale = 1./grid%nelem**0.5

    if(ityp == 4) scale = 0.25

    if(state%space%estim_space == 'DWR') then
      print*, 'DWR TOL set as 1/#T_h in MarkElements, ityp = 2'
      ityp = 2
      scale = 1./grid%nelem
    endif

    !do i=1,nelem
    !   grid%elem(i)%estim = grid%elem(i)%errL2
    !   grid%elem(i)%estim1 = grid%elem(i)%errL2
    !enddo

    tol_max = state%space%adapt%tol_max
    tol_min = state%space%adapt%tol_min

!    print*, 'dsafd', tol_max, tol_min
    !write(*,'(a32,3es12.4)') ' # Tolerances for adaptation = ', tol_min,tol_max, scale*tol_max


    rez_estim = 0.

    do i = 1, grid%nelem
       elem => grid%elem(i)

       elem%hsplit = 0
       elem%psplit = 0
       pK = max(1,elem%deg)

       if (state%space%estim_space /= 'pNeu') then
          ! indicator of the local error
          if( state%time%disc_time /= 'STDG') then
             elem%estim_loc = (elem%estim_loc + elem%jumpsJh)**0.5
          elseif (state%space%estim_space == 'DWR') then
            ! no need of the sqrt in DWR
!            print*, i,'estim_loc:' , elem%estim_loc, 'rez_estim:',rez_estim

          else
             elem%estim_loc = elem%estim_loc**0.5
          endif

          ! indicator of the regularity
          elem%reg = elem%rezid / elem%area / elem%diam**(2*pK)* elem%diam**3
          !elem%reg = elem%rezid / ( elem%area * elem%diam**(2*pK+1) )

          ! seting tolerances for h-p decision
          elem%regT0 = elem%diam**(-2)    ! elem%reg <  elem%regT0  ==> p-refinement
          elem%regT1 = elem%diam**(-3)    ! elem%reg <  elem%regT1  ==> h->p substitution
          elem%regT2 = elem%diam**(-4)    ! elem%reg >  elem%regT2  ==> p-DErefinement

          !elem%regT0 = elem%diam**(-1)    ! elem%reg <  elem%regT0  ==> p-refinement
          !elem%regT1 = elem%diam**(-2)    ! elem%reg <  elem%regT1  ==> h->p substitution
          !elem%regT2 = elem%diam**(-3)    ! elem%reg >  elem%regT2  ==> p-DErefinement

          !elem%regT0 = elem%diam**( -0)    ! elem%reg <  elem%regT0  ==> p-refinement
          !elem%regT1 = elem%diam**( -1)    ! elem%reg <  elem%regT1  ==> h->p substitution
          !elem%regT2 = elem%diam**( -2)    ! elem%reg >  elem%regT2  ==> p-DErefinement

       elseif (state%space%estim_space == 'pNeu') then
          ! definitions of elem%eta(resST, 1), elem%reg, elem%regT0, ...
          ! MOVED directly to neumann.f90

          elem%estim_loc = elem%eta(resST, 1)
       endif

       if (state%space%estim_space == 'DWR') then
         rez_estim = rez_estim + elem%estim_loc
       else
         rez_estim = rez_estim + elem%estim_loc**2
       endif


       if(ityp == 3) scale = (elem%area/state%space%domain_volume)**0.5
       !write(100+state%space%adapt%adapt_level,*) elem%xc(:),elem%estim_loc,tol_max *scale, &
       !     elem%rezid,  elem%reg, elem%regT0, elem%regT1, elem%regT2, elem%diam, elem%area
       !!, elem%estim_loc**0.5, &
       !!elem%jumpsJh**0.5, (elem%estim_loc + elem%jumpsJh)**0.5, elem%reg, &
       !!     elem%reg, elem%reg/elem%diam, &
       !!     tol_max*(elem%area/state%space%domain_volume)**0.5, elem%errL2, &
       !!     (elem%errL2**2 + elem%errH1**2*state%model%Re1)**0.5, &
       !!     elem%rezid, elem%rezid*elem%limit_par, tol_max *scale , tol_min *scale




       ! IGNORING THE BOUNDARY CONDITIONS FOR SHOCK-VORTEX INTERACTION
       ! performed in subroutine ElementEdgeJumpProj in euler.f90
       !!if( state%type_IC == 8 )  elem%estim_loc = elem%estim_loc


    enddo

    if(ityp == 4) tol_max = maxval(grid%elem(:)%estim_loc)

    if (state%space%estim_space == 'DWR') then
      !do nothing
    else
      rez_estim = rez_estim**0.5
    endif

    ! for state%space%adapt%adapt_type == 'Ihp' the refinement is given in a different way directly
    ! in   subroutine IsotropicMetric( )
    if(state%space%adapt%adapt_type == 'Ihp') then
       print*,'Marking of elements for refinement done'

       RETURN

    endif


    ! marking element for hp-refinement
    do i = 1, grid%nelem
       elem => grid%elem(i)

       if(ityp == 3) scale = (elem%area/state%space%domain_volume)**0.5

       if(elem%estim_loc > tol_max *scale  ) then

!          if (state%space%estim_space == 'DWR') then
!            print*, 'Element' , i , 'is marked: ' , elem%estim_loc , '>' , tol_max *scale
!          endif

          if(elem%reg <= elem%regT0 .and. elem%deg < MaxDegreeImplemented-1 ) then
             ! solution is regular => p-refinement
             elem%psplit = 1

             !write(200+state%space%adapt%adapt_level,*) elem%xc(:),elem%estim_loc,elem%estim_loc,elem%reg

          else
             ! solution isn't regular => h-refinement
             elem%hsplit = 4

             if(elem%reg >= elem%regT2 .and. elem%deg > MinDegreeImplemented) &
                  elem%psplit = -1

             !write(300+state%space%adapt%adapt_level,*) elem%xc(:),elem%estim_loc,elem%estim_loc,elem%reg

          endif

       elseif(elem%estim_loc < tol_min*scale ) then

!          if (state%space%estim_space == 'DWR') then
!            print*, 'Element' , i , 'is NOT marked! '
!          endif
          ! NOT YET CHECKED
          if( elem%deg > MinDegreeImplemented)  elem%psplit = -1

          !write(400+state%space%adapt%adapt_level,*) elem%xc(:),elem%estim_loc,elem%estim_loc,elem%reg, elem%regT2
          !write(*,'(8es12.4)') elem%xc(:),elem%estim_loc,elem%estim_loc,elem%reg, elem%regT2
       endif

       !elem%hsplit = 4
       !elem%psplit = 0
       !elem%hsplit = 0
       !elem%psplit = 1

    enddo


    if (only_h_adapt) then  ! only h-refinement
       do i = 1, grid%nelem
          elem => grid%elem(i)
          if(elem%psplit == 1) then
             elem%psplit = 0
             elem%hsplit = 4
          elseif(elem%psplit == -1) then
              elem%psplit = 0
           endif
        enddo
     elseif (only_p_adapt) then  ! only p-refinement
       do i = 1, grid%nelem
          elem => grid%elem(i)
          if(elem%hsplit == 4) then
             elem%psplit = 1
             elem%hsplit = 0
          elseif(elem%hsplit == -4) then
              elem%hsplit = 0
           endif
        enddo
     endif



     ! no refinement on boundary outside of shock :-)))
     if( state%type_IC == 8 )  then
        do i = 1, grid%nelem
           elem => grid%elem(i)
           if( (abs(elem%xc(1) - 1.) > 0.02) .and.  &
                (elem%xc(2) > 1.75 .or. elem%xc(2) < 0.25) )elem%hsplit = 0
           !if( (abs(elem%xc(1) - 1.) > 0.05) .and.  &
           !     (elem%xc(2) > 1.8 .or. elem%xc(2) < 0.25) ) write(88,*) elem%xc
        enddo
     endif


    goto 149 ! WE SKIP DEREFINEMENT   !!!!!!
     if(only_p_adapt) goto 150
     if(state%space%adapt%adapt_type == 'RG') goto 150   ! NOT YET IMPLEMENTED

    ! marking element for h-derefinement
    do i = 1, grid%nelem
       elem => grid%elem(i)
       lev = elem%RGhistory

       ! only elements which were refined in previous computations and which are  direct
       ! daughters of their mother elements will be checked
       if(lev > 0 ) then
          if(elem%i == elem%RG(lev)%daughter(4)) then

             ! we go though all daughter elements including the original one
             do j=1,elem%RG(lev)%subel

                !write(*,'(a8,5i5)') 'elmD:',i, elem%i, lev, elem%RG(lev)%daughter(4)

                elemD => grid%elem(elem%RG(lev)%daughter(j))
                if(ityp == 3) scale = (elemD%area/state%space%domain_volume)**0.5

                if(elemD%estim_loc > tol_min * scale  &
                     .or. elem%RGlevel /= elemD%RGlevel) goto 200

                if(elemD%per > 0) then
                   goto 200 ! DEREFINEMENT NOT IMPLEMENTED FOR PERIODIC
                   !elem1 => grid%elem(elemD%per)
                   !if(elem1%estim_loc > tol_min * scale  &
                   !  .or. elem%RGlevel /= elem1%RGlevel)
                endif



                ! NO derefinement
             enddo

             ! the error estimate of all daughter elements are smaller than tol_min,
             ! we mark them for derefinement

             grid%elem(elem%RG(lev)%daughter(1:elem%RG(lev)%subel))%hsplit = - 4
             goto 100
200          continue

             if (only_h_adapt) goto 100 ! NO h -> p substitution

             ! h -> p substitution
             ! we go though all daughter elements including the original one
             do j=1,elem%RG(lev)%subel
                elemD => grid%elem(elem%RG(lev)%daughter(j))

                !!if(elemD%hsplit /= 0 .or.  elemD%reg > 1.&

                if(elemD%deg ==  MaxDegreeImplemented-1 .or. &
                     elemD%hsplit /= 0 .or. elemD%psplit /= 0 .or. &!B elemD%reg > 1. .or. &
                     elemD%reg > elemD%regT1 .or. elem%RGlevel /= elemD%RGlevel) goto 100
                ! NO adaptation
             enddo
             ! all daughter elements are not marker for p- nor h- adaptation,
             ! solution is regular, we derefine them and increase p

             !print*,'???',lev, elem%i, elemD%i

             grid%elem(elem%RG(lev)%daughter(1:elem%RG(lev)%subel))%hsplit = - 4
             grid%elem(elem%RG(lev)%daughter(1:elem%RG(lev)%subel))%psplit = + 1

             !write(600 +state%space%adapt%adapt_level,*) elem%xc(:)

100          continue

          endif
       endif
    enddo

149 continue

    !goto 150 ! for hp_steady, not for st_interpol !!!!!!!!!!!!!!

    ! if all neighbours of elements are marked for h-refinement, we marke also the
    ! original element
    do ib= 1, 3
       ipoc = 0
       do i = 1, grid%nelem
          elem => grid%elem(i)
          if(elem%hsplit == 0 ) then
             lev = 0
             do j=1,elem%flen
                if(elem%face(neigh, j) > 0 ) then
                   elemD => grid%elem(elem%face(neigh, j))
                   if(elemD%hsplit == 4) lev = lev + 1
                endif
             enddo
             !if(lev == elem%flen )then
             if(lev >= elem%flen -1 ) then
                !write(*,'(a6,4i5,3es12.4)') '>>!<<',ib,elem%i, lev, elem%flen, elem%xc(:)
                elem%hsplit = 4
                ipoc = ipoc + 1
             endif

             !if(lev == elem%flen ) write(99,*) elem%xc(:)
          endif

          ! if all neighbours has a larger level of refinent
          if(elem%hsplit == 0 ) then
             lev = 1
             do j=1,elem%flen
                if(elem%face(neigh, j) > 0 ) then
                   elemD => grid%elem(elem%face(neigh, j))
                   if(elemD%RGlevel <= elem%RGlevel) lev = 0
                endif
             enddo
             if(lev >= 1) then
                !write(*,'(a6,4i5,3es12.4)') '>><<',ib, elem%i, lev, elem%flen, elem%xc(:)
                elem%hsplit = 4
                ipc = ipoc + 1
             endif

          endif
       enddo
       if(ipoc == 0) goto 150
    enddo

    !if(state%space%adapt%adapt_level == 1) grid%elem(:)%psplit = 0
    !if(state%space%adapt%adapt_level == 1) grid%elem(:)%hsplit = -4

    !!!!call RedGreenMarking( )

    !!grid%elem(1:grid%nelem)%hsplit = 0
    !call PlotMesh(grid)  ! necessary !!!!!!!!!!!!
    !call PlotMarkedElements(1000+state%space%adapt%adapt_level)

150 continue

    !call PlotMarkedElements(1000+state%space%adapt%adapt_level)

    ! h-refinement for periodic elements, both periodic elements
    do i=1,grid%nelem
       elem => grid%elem(i)
       if(elem%per > 0) then
          elem1 => grid%elem(elem%per)

          if(elem%hsplit == 4) then
             elem1%hsplit = 4
             !write(*,*) elem%xc(:)
             !write(*,*)elem1%xc(:)
          endif

       endif
    enddo

    !      endif
    !   enddo
    !end if


    !do i=1,grid%nelem
    !   elem => grid%elem(i)
    !   if(elem%hsplit == 4)  write(800 +state%space%adapt%adapt_level,*) elem%xc(:)
    !   if(elem%hsplit == -4)  write(700 +state%space%adapt%adapt_level,*) elem%xc(:)

    !   !elem%psplit = 0
    !   !elem%hsplit = 0
    !   !if(state%space%adapt%adapt_level == 0 .and. elem%i == 7) elem%hsplit = 4
    !   !if(state%space%adapt%adapt_level == 1 .and. elem%i == 11) elem%hsplit = 4
    !enddo

    if(state%space%adapt%adapt_type == 'HG') then
       call HGMarking()
    elseif(state%space%adapt%adapt_type == 'RG') then
       call RedGreenMarking( )
    else
       print*,'UNKNOWN marking strategy in marking.f90'
    endif

    !call PlotMarkedElements(2000+state%space%adapt%adapt_level)

    !if(state%space%adapt%adapt_level >= 0) then
    !   grid%elem(:)%hsplit = 0
    !   grid%elem(:)%psplit = 0
    !endif


    state%space%adaptation = .false.   ! any adaptation
    state%space%adaptationR = .false.  ! only refinement
    if(count(grid%elem(1:grid%nelem)%hsplit == 4) > 0) state%space%adaptation = .true.
    if(count(grid%elem(1:grid%nelem)%hsplit == -4) >= 12) state%space%adaptation = .true.
    !if(count(grid%elem(1:grid%nelem)%hsplit == -4) >= 40) state%space%adaptation = .true.
    if(count(grid%elem(1:grid%nelem)%psplit /= 0) > 0) state%space%adaptation = .true.

    if(count(grid%elem(1:grid%nelem)%hsplit == 4) > 0) state%space%adaptationR = .true.
    if(count(grid%elem(1:grid%nelem)%psplit /= 0) > 0) state%space%adaptationR = .true.

    !print*,'##### ADAPT?', state%space%adaptation, count(grid%elem(1:grid%nelem)%hsplit == 4),  &
    !     count(grid%elem(1:grid%nelem)%hsplit == -4), &
    !     count(grid%elem(1:grid%nelem)%psplit /= 0)

    if(.not. state%space%adaptation) print*,'# NO (de)refinement detected'


  end subroutine MarkElements

  subroutine PlotMarkedElements(isol)
    integer :: isol
    character(len=15) :: solfile
    character(len=15) :: hpsplit
    integer :: i, chlen
    character(len=1) :: ch1
    character(len=2) :: ch2
    character(len=3) :: ch3
    character(len=4) :: ch4
    character(len=5) :: ch5

    solfile = 'rge'

    if(isol >= 0 .and. isol <= 9) then
       write( ch1, '(i1)' ) isol
       solfile(4:4) = ch1
       chlen = 4
    elseif(isol >= 10 .and. isol <= 99) then
       write( ch2, '(i2)' ) isol
       solfile(4:5) = ch2
       chlen = 5
    elseif(isol >= 100 .and. isol <= 999) then
       write( ch3, '(i3)' ) isol
       solfile(4:6) = ch3
       chlen = 6
    elseif(isol >= 1000 .and. isol <= 9999) then
       write( ch4, '(i4)' ) isol
       solfile(4:7) = ch4
       chlen = 7
    elseif(isol >= 10000 .and. isol <= 99999) then
       write( ch5, '(i5)' ) isol
       solfile(4:8) = ch5
       chlen = 8
    endif

    !!call PlotMesh(grid, 'mesh0')

    open(14,file=solfile, status="replace")
    open(15,file="marked.gn", status="replace")
    do i=1,grid%nelem
       hpsplit(1:7) = "..' at "
       if (grid%elem(i)%hsplit ==  4) hpsplit(1:1) = '+'
       if (grid%elem(i)%hsplit == -4) hpsplit(1:1) = '-'

       if (grid%elem(i)%psplit ==  1) hpsplit(2:2) = '+'
       if (grid%elem(i)%psplit == -1) hpsplit(2:2) = '-'


       !if (grid%elem(i)%i >= 0) then
       !if (grid%elem(i)%xc(2) < 0.2 .and. abs(grid%elem(i)%xc(1) -0.4) < 0.2) then
       write(14, *) grid%elem(i)%xc(:),i,grid%elem(i)%hsplit
       write(15,'(a13,i1,a1,i5,a8,es14.6,a3,es14.6, a25)') &
            "set label '(",grid%elem(i)%deg,")", &
            i,hpsplit(1:7),grid%elem(i)%xc(1)," , ",grid%elem(i)%xc(2), &
            ' center font "Symbol,6"'
    enddo
    close(14)
    !write(15,*) 'set terminal postscript landscape'
    write(15,*) 'set terminal postscript eps color'
    !write(15,*) 'set size 2., 2.'
    !write(15,*) 'set output "marked.ps" '
    write(15,*) 'set output "',solfile(1:chlen),'.eps" '
    !write(15,*) "plot 'mesh' w l"
    !write(15,*) "plot [0.2:0.65][0.:0.35]'mesh' w l"
    write(15,*) "plot [0.:1][0.:1]'mesh0' w l"
    !write(15,*) "plot [0.47:0.6][0.2:0.3]'mesh' w l"
    !write(15,*) "plot 'mesh' w l"
    !write(15,*) "plot [0.47:0.57] [0.22:0.33]  'mesh' w l"
    !write(15,*) "plot [0.35:0.55][0.1:0.25] 'mesh' w l"
    !write(15,*) "plot 'mesh' w l, '",solfile(1:chlen),"'"
    close(15)

    call system('gnuplot -persist "marked.gn"')

  end subroutine PlotMarkedElements

  !> marking the elements for the reg green refinement
  !> elem%hsplit = 4      +=> red refinement
  !> elem%hsplit = 1,2,3  ==>  1,2 or 3rd side will be refined green
  subroutine RedGreenMarking( )
    class(element), pointer :: elem, elem1
    integer :: i,j
    integer :: imark, imark2, isum, ifile
    integer :: isum1

    !ifile = 10
    ifile = 100+state%space%adapt%adapt_level*10




    !do i=1,grid%nelem
    !      elem => grid%elem(i)
    !      if(elem%hsplit == 4 )  write(ifile,*) elem%xc(:)
    !   end do

    !!call PlotMarkedElements(ifile)


    imark2 = 1
    do
       if (imark2 == 0) exit
       imark2 = 0

       do i=1,grid%nelem
          elem => grid%elem(i)

          if(elem%hsplit == 4) then
             isum = sum(grid%elem(elem%face(neigh,1:elem%flen))%psplit, &
                  mask = elem%face(neigh,:)> 0 )

             isum1 = sum(elem%face(neigh,:)/elem%face(neigh,:), mask = elem%face(neigh,:)> 0 )

             !print*,'%%%%%%%%%',elem%i, isum, isum1

             if(isum >= isum1 -1 ) then
                elem%hsplit = 0
                elem%psplit = 1
                imark = imark+1
                !!write(ifile,*) elem%xc(:)
             endif

          endif

          !!call PlotMarkedElements(ifile)

          !!print*,'in ###',ifile,imark2, imark


       enddo
    enddo



    ! seeking additional elements for red refinement
    imark2 = 1
    do
       if (imark2 == 0) exit
       imark2 = 0


       imark = 1
       do
          if (imark == 0) exit
          ifile = ifile + 1
          imark = 0
          do i=1,grid%nelem
             elem => grid%elem(i)

             isum = sum(grid%elem(elem%face(neigh,1:elem%flen))%hsplit, &
                  mask = elem%face(neigh,:)> 0 )

             if(isum > 4 .and. elem%hsplit < 4 ) then
                elem%hsplit = 4
                imark = imark+1
                !!write(ifile,*) elem%xc(:)
             endif

          enddo

          !!call PlotMarkedElements(ifile)

          !!print*,'in ###',ifile,imark2, imark


       enddo
       !!print*,'Inner loop exited',ifile

       !do i=1,grid%nelem
       !   elem => grid%elem(i)
       !   if(elem%hsplit == 4 )
       !end do

       ! marking the elements for green refinement
       ifile = ifile + 1


       do i=1,grid%nelem
          elem => grid%elem(i)
          if(elem%hsplit == 0 ) then
             do j=1,elem%flen
                if(elem%face(neigh,j)> 0 ) then
                   if(grid%elem(elem%face(neigh,j))%hsplit == 4) &
                        elem%hsplit = j
                endif
             enddo
             if(elem%RGtype == 'G' .and. elem%hsplit > 0) then
                elem%hsplit = 4
                imark2 = imark2 + 1
                !!write(ifile,*) elem%xc(:)
             endif

          endif
       enddo

       !print*,'out###',ifile,imark2, imark

    enddo  ! outer cycles imark2

    ! do i=1,grid%nelem
    !    elem => grid%elem(i)
    !    if(elem%hsplit == 4 ) &
    !         write(1000+10*state%space%adapt%adapt_level +1,*) elem%xc(:)
    !    if(elem%hsplit > 0 .and. elem%hsplit < 4  ) &
    !         write(1000+10*state%space%adapt%adapt_level +2,*) elem%xc(:)
    ! end do

  end subroutine RedGreenMarking

  !> marking the elements for the reg HG refinement,
  !> if elem%HGlevel == state%space%adapt%max_HGlevel_implemented and any neighbouring element
  !> is marked then we refine also the actual one element
  subroutine HGMarking( )
    class(element), pointer :: elem, elemD
    integer, dimension(:), allocatable :: RGlev
    integer :: i, ie, j, ip, id, j1
    integer :: imark, isum, ifile, HGlevel, RGhis
    integer :: minRGlev, maxRGlev
    logical :: derefine
    integer :: itest


    ! marking for p-refinement
    do i=1,grid%nelem
       elem => grid%elem(i)
       if(elem%psplit == 1 .and. elem%deg == MaxDegreeImplemented-1  .or. &
            elem%psplit == -1 .and. elem%deg == MinDegreeImplemented ) then
          elem%psplit = 0
       endif
    enddo

    ! marking for h-refinement
    allocate(RGlev(1:grid%nelem) )  ! future RG level

    do i=1,grid%nelem
       elem => grid%elem(i)

       if(elem%hsplit ==  4 .and. elem%RGlevel >= max_RGlevel) elem%hsplit = 0

       RGlev(i) = elem%RGlevel

       if(elem%hsplit ==  4 ) RGlev(i) = RGlev(i) + 1
       if(elem%hsplit == -4 ) RGlev(i) = RGlev(i) - 1

    enddo


    ! in order to avoid maximal number of HG nodes, we add some additional refinement
    imark = 1
    do
       if (imark == 0) exit
       imark = 0
       do i=1,grid%nelem
          elem => grid%elem(i)
          if(elem%hsplit <= 0 .and. elem%HGnode) then

             do j=1,elem%flen
                if(elem%face(neigh, j) > 0) then

                   if(RGlev(elem%face(neigh, j)) - RGlev(i)  > state%space%adapt%max_HGlevel) then
                      if(elem%hsplit == 0) RGlev(i) = RGlev(i) + 1
                      if(elem%hsplit <  0) RGlev(i) = RGlev(i) + 2
                      elem%hsplit = 4
                      imark = imark+1
                      !print*,'....',elem%i,elem%RGlevel,RGlev(i)
                      goto 100
                   endif
                endif
             enddo
100          continue
          endif
       enddo
       !print*,'!!!',imark
    enddo


    ! if one of the element originally marked for derefinement was in the
    ! previous cyclus marked for refinement
    ! then all sister elements will be unmarked for derefinement

    do i=1,grid%nelem
       elem => grid%elem(i)
       RGhis = elem%RGhistory

       derefine = .false.
       if(elem%hsplit < 0 .and. RGhis > 0) then
          ! we go though all daughter elements NOT including the original one
          do id=1,elem%RG(RGhis)%subel
             elemD => grid%elem(elem%RG(RGhis)%daughter(id))

             if(elem%RGlevel /= elemD%RGlevel) goto 15 ! NOT real sisters

             if(elem%hsplit /= elemD%hsplit) derefine = .true.
          enddo

          !if(derefine) then
          !   print*,'strange in estimation.f90'
          !   write(*,'(20i5)') elem%i, elem%hsplit,elem%RGlevel
          !   do id=1,elem%RG(RGhis)%subel
          !      elemD => grid%elem(elem%RG(RGhis)%daughter(id))
          !      write(*,'(20i5)') elemD%i, elemD%hsplit,elemD%RGlevel
          !   enddo
          !endif

          if(derefine) then
             do j=1,elem%RG(RGhis)%subel
                if(grid%elem(elem%RG(RGhis)%daughter(j) )%hsplit < 0) then
                   grid%elem(elem%RG(RGhis)%daughter(j) )%hsplit = 0

                   RGlev(elem%RG(RGhis)%daughter(j)) = &
                        RGlev(elem%RG(RGhis)%daughter(j)) + 1
                endif
             enddo
          endif
       endif
15     continue
    enddo

    ! removing of derefinement in order to avoid overflow of max_HGlevel
    imark = 1
    do
       if (imark == 0) exit
       imark = 0

       do i=1,grid%nelem
          elem => grid%elem(i)
          RGhis = elem%RGhistory

          if(elem%hsplit < 0 .and. RGhis > 0) then
             minRGlev = max_RGlevel + 1
             maxRGlev = 0

             ! we go though all daughter elements NOT including the original one
             do id=1,elem%RG(RGhis)%subel
                elemD => grid%elem(elem%RG(RGhis)%daughter(id))

                if(elem%RGlevel /= elemD%RGlevel) goto 20 ! NOT real sisters

                if(elemD%i /= elem%i) then
                   do j=1, elemD%flen
                      if(elemD%face(neigh,j) > 0 .and. elemD%face(neigh,j) /= i) then

                         minRGlev = min(minRGlev, RGlev(elemD%face(neigh,j)) )
                         maxRGlev = max(maxRGlev, RGlev(elemD%face(neigh,j)) )

                      endif
                   enddo
                endif
             enddo

             if(RGlev(i) < minRGlev  &
                  .or.  maxRGlev - RGlev(i)  > state%space%adapt%max_HGlevel) then

                grid%elem(elem%RG(RGhis)%daughter(1:elem%RG(RGhis)%subel) )%hsplit = 0

                RGlev(elem%RG(RGhis)%daughter(1:elem%RG(RGhis)%subel)) = &
                     RGlev(elem%RG(RGhis)%daughter(1:elem%RG(RGhis)%subel)) + 1
                imark = imark + 1
             endif
          end if
20        continue
       enddo

    enddo


    ! checking the maximal level of adaptation
    do i=1,grid%nelem
       elem => grid%elem(i)
       if(RGlev(i) > max_RGlevel) grid%elem(i)%hsplit = 0
    enddo


    deallocate(RGlev)
  end subroutine HGMarking


  !> marking the elements for the reg HG refinement,
  !> if elem%HGlevel == state%space%adapt%max_HGlevel and any neighbouring element
  !> is marked then we refine also the actual one element
  subroutine HGMarkingOLD( )
    class(element), pointer :: elem, elemD
    integer :: i, ie, j, ip, id, j1
    integer :: imark, isum, ifile, HGlevel, RGhis
    integer :: minRGlev, iRGlev
    logical :: refine
    integer :: itest

    itest = 2727

    ifile = 0

    i = itest
    if(i <= grid%nelem) &
         print*,'??? B0',i, grid%elem(i)%hsplit, grid%elem(i)%RGlevel, max_RGlevel

    ! checking the maximal level of adaptation
    do i=1,grid%nelem
       elem => grid%elem(i)
       if(elem%hsplit == 4 .and. &
            elem%RGlevel == max_RGlevel) elem%hsplit = 0
    enddo

    i = itest
    if(i <= grid%nelem) &
         print*,'??? B1',i, grid%elem(i)%hsplit, grid%elem(i)%RGlevel, max_RGlevel

    ! checking if all neighbours of unmarked element has higher RGlevel or are marked
    do i=1,grid%nelem
       elem => grid%elem(i)
       if(elem%hsplit == 0 ) then

          minRGlev = max_RGlevel + 1

          do j=1,elem%flen
             if(elem%face(neigh,j) > 0) then

                iRGlev = grid%elem(elem%face(neigh,j))%RGlevel
                if(grid%elem(elem%face(neigh,j))%hsplit == 4) iRGlev = iRGlev + 1

                minRGlev = min(minRGlev, iRGlev)
             endif
          enddo

          !if(elem%RGlevel < minRGlev .and. elem%RGlevel < max_RGlevel) then
          !   elem%hsplit = 4
          !   write(1000+state%space%adapt%adapt_level,'(2es12.4,a8,6i5)')&
          !        grid%elem(elem%face(neigh,j))%xc(:),&
          !        ' elem : ',elem%RGlevel, elem%hsplit
          !endif
       endif
    enddo

    !i = itest
    !if(i <= grid%nelem) &
    !    print*,'??? B2',i, grid%elem(i)%hsplit, grid%elem(i)%RGlevel, max_RGlevel

    !i = 2337
    !if(i <= grid%nelem) &
    !    print*,'??? B2',i, grid%elem(i)%hsplit, grid%elem(i)%RGlevel, max_RGlevel

    ! checking if all neighbours of elements marked for derefinement
    ! were refined or are marked
    do i=1,grid%nelem
       elem => grid%elem(i)
       RGhis = elem%RGhistory

       ! only elements which were refined in previous computations and
       ! which are  direct
       ! daughters of their mother elements will be checked
       if( elem%hsplit <= 0 .and. RGhis > 0) then
          minRGlev = max_RGlevel + 1

          ! we go though all daughter elements NOT including the original one
          do id=1,elem%RG(RGhis)%subel
             elemD => grid%elem(elem%RG(RGhis)%daughter(id))

             if(elem%RGlevel /= elemD%RGlevel) goto 20 ! NOT real sisters

             if(elemD%i /= elem%i) then
                do j=1, elemD%flen
                   if(elemD%face(neigh,j) > 0 .and. elemD%face(neigh,j) /= i) then

                      ! RG level of neighbouring element will be:
                      iRGlev = grid%elem(elemD%face(neigh,j))%RGlevel
                      if(grid%elem(elemD%face(neigh,j))%hsplit == 4) &
                           iRGlev = iRGlev + 1

                      minRGlev = min(minRGlev, iRGlev)

                      !write(*,'(a4,3i5,a1,6i5)') 'der : ',elem%i,elemD%i,j,'|', &
                      !     grid%elem(elemD%face(neigh,j))%RGlevel, elemD%RGlevel, &
                      !     grid%elem(elemD%face(neigh,j))%hsplit

                      !if(grid%elem(elemD%face(neigh,j))%RGlevel <= elemD%RGlevel .and. &
                      !     grid%elem(elemD%face(neigh,j))%hsplit <= 0) goto 20
                   endif
                enddo
             endif
          enddo

          if(minRGlev >= elem%RGlevel)  then ! no derefinement
             grid%elem(elem%RG(RGhis)%daughter(1:elem%RG(RGhis)%subel))%hsplit = 0
                          do id=1,elem%RG(RGhis)%subel
                elemD => grid%elem(elem%RG(RGhis)%daughter(id))

                !if(elemD%RGlevel == max_RGlevel) then
                   !print*,'Troubles in estimation (2)', id
                   write(2000+state%space%adapt%adapt_level,'(2es12.4,a8,6i5)') elemD%xc(:),&
                        ' elemD = ',elemD%i, elemD%RGlevel, elemD%hsplit
                   do j=1,elemD%flen
                      if(elemD%face(neigh,j) > 0) then
                         write(2000+state%space%adapt%adapt_level,'(2es12.4,a8,6i5)')&
                              grid%elem(elemD%face(neigh,j))%xc(:),&
                              ' elem1 : ',elemD%face(neigh,j), &
                              grid%elem(elemD%face(neigh,j))%RGlevel, &
                              grid%elem(elemD%face(neigh,j))%hsplit

                      endif

                   end do
                   write(2000+state%space%adapt%adapt_level,'(x)')
                !endif
             enddo

          endif


          if(minRGlev > elem%RGlevel .and. &
               elem%RGlevel < max_RGlevel ) then ! additional refinement
             grid%elem(elem%RG(RGhis)%daughter(1:elem%RG(RGhis)%subel))%hsplit = 4

             !if( elem%RG(RGhis)%daughter(1) == itest  .or. &
             !    elem%RG(RGhis)%daughter(2) == itest  .or. &
             !    elem%RG(RGhis)%daughter(3) == itest  .or. &
             !    elem%RG(RGhis)%daughter(4) == itest   ) then!

             !print*,'??? A3',i, elem%hsplit, elem%RGlevel, max_RGlevel
             !   do id=1,elem%RG(RGhis)%subel
             !      elemD => grid%elem(elem%RG(RGhis)%daughter(id))

             !      write(*,'(2es12.4,a8,6i5)') elemD%xc(:),&
             !           ' elemD = ',elemD%i, elemD%RGlevel, elemD%hsplit
             !      do j=1,elemD%flen
             !         if(elemD%face(neigh,j) > 0) then
             !            write(*,'(2es12.4,a8,6i5)')&
             !                 grid%elem(elemD%face(neigh,j))%xc(:),&
             !                 ' elem1 : ',elemD%face(neigh,j), &
             !                 grid%elem(elemD%face(neigh,j))%RGlevel, &
             !                 grid%elem(elemD%face(neigh,j))%hsplit
             !
             !         endif
             !      enddo
             !   enddo
             !endif

             !write(*,'(2es12.4,a6,3i5)') elem%xc(:),'%%%%%',elem%i,elem%RGlevel, minRGlev
          ! NO DErefinement

             do id=1,elem%RG(RGhis)%subel
                elemD => grid%elem(elem%RG(RGhis)%daughter(id))

                !if(elemD%RGlevel == max_RGlevel) then
                   !print*,'Troubles in estimation (2)', id
                   write(3000+state%space%adapt%adapt_level,'(2es12.4,a8,6i5)') elemD%xc(:),&
                        ' elemD = ',elemD%i, elemD%RGlevel, elemD%hsplit
                   do j=1,elemD%flen
                      if(elemD%face(neigh,j) > 0) then
                         write(3000+state%space%adapt%adapt_level,'(2es12.4,a8,6i5)')&
                              grid%elem(elemD%face(neigh,j))%xc(:),&
                              ' elem1 : ',elemD%face(neigh,j), &
                              grid%elem(elemD%face(neigh,j))%RGlevel, &
                              grid%elem(elemD%face(neigh,j))%hsplit

                      endif

                   end do
                   write(3000+state%space%adapt%adapt_level,'(x)')
                !endif
             enddo

          endif

          !print*,grid%elem(elem%RG(RGhis)%daughter(1:elem%RG(RGhis)%subel))%hsplit,'****'
20        continue
       endif
    enddo

    i = itest
    if(i <= grid%nelem) &
         print*,'??? B3',i, grid%elem(i)%hsplit, grid%elem(i)%RGlevel, max_RGlevel


    ! in order to avoid maximal number of HG nodes, we remove some derefinement
    imark = 1
    do
       if (imark == 0) exit
       !call PlotMarkedElements(ifile)
       !ifile = ifile + 1

       imark = 0
       do i=1,grid%nelem
          elem => grid%elem(i)
          RGhis = elem%RGhistory

          ! only elements which were refined in previous computations and
          ! which are  direct
          ! daughters of their mother elements will be checked
          if(elem%hsplit == -4 .and. RGhis > 0) then
             HGlevel = 0

             ! we go though all daughter elements including the original one
             do id=1,elem%RG(RGhis)%subel
                elemD => grid%elem(elem%RG(RGhis)%daughter(id))
                if(elemD%HGnode) then
                   do ie=1,elemD%type
                      do j=1 ,elemD%flen   !HERE may be accelerated !!!!!!!!!!!!!!!!!!!!

                         if(elemD%HGface(1,j) == ie) then
                            ip = log(1.0001*elemD%HGface(2, j) ) / log(2.)

                            ! if neighbour will be splited, HG level is increasing
                            if(elemD%face(neigh,j) > 0) then
                               if(grid%elem(elemD%face(neigh,j))%hsplit == 4) ip = ip + 1
                            endif

                            HGlevel = max(ip, HGlevel)
                         endif
                      enddo
                   enddo
                endif
             enddo

             !write(*,'(a6,20i5)') 'MarkHG',elem%RG(RGhis)%daughter(:)
             !write(*,'(a6,20i5)') 'MarkHG',HGlevel, state%space%adapt%max_HGlevel


             if(HGlevel >= state%space%adapt%max_HGlevel) then
                grid%elem(elem%RG(RGhis)%daughter(1:elem%RG(RGhis)%subel))%hsplit = 0

                do id=1,elem%RG(RGhis)%subel
                   elemD => grid%elem(elem%RG(RGhis)%daughter(id))

                !if(elemD%RGlevel == max_RGlevel) then
                   !print*,'Troubles in estimation (2)', id
                   write(4000+state%space%adapt%adapt_level,'(2es12.4,a8,6i5)') elemD%xc(:),&
                        ' elemD = ',elemD%i, elemD%RGlevel, elemD%hsplit
                   do j=1,elemD%flen
                      if(elemD%face(neigh,j) > 0) then
                         write(4000+state%space%adapt%adapt_level,'(2es12.4,a8,6i5)')&
                              grid%elem(elemD%face(neigh,j))%xc(:),&
                              ' elem1 : ',elemD%face(neigh,j), &
                              grid%elem(elemD%face(neigh,j))%RGlevel, &
                              grid%elem(elemD%face(neigh,j))%hsplit

                      endif

                   end do
                   write(4000+state%space%adapt%adapt_level,'(x)')
                !endif
             enddo

             endif
          endif
200       continue
       enddo
    enddo

    i = itest
    if(i <= grid%nelem) &
         print*,'??? B4',i, grid%elem(i)%hsplit, grid%elem(i)%RGlevel, max_RGlevel

    ! in order to avoid maximal number of HG nodes, we add some additional refinement
    imark = 1
    do
       if (imark == 0) exit
       !call PlotMarkedElements(ifile)
       !ifile = ifile + 1


       imark = 0
       do i=1,grid%nelem
          elem => grid%elem(i)

          if(elem%hsplit == 0 .and. elem%HGnode) then

             do ie=1,elem%type
                HGlevel = 0
                do j=1 ,elem%flen   !HERE may be accelerated !!!!!!!!!!!!!!!!!!!!
                   if(elem%HGface(1,j) == ie) then
                      ip = log(1.0001*elem%HGface(2, j) ) / log(2.)
                      HGlevel = max(ip, HGlevel)
                   endif
                enddo

                if(HGlevel >= state%space%adapt%max_HGlevel) then
                   do j=1 ,elem%flen   !HERE may be accelerated !!!!!!!!!!!!!!!!!!!!

                      if(elem%HGface(1,j) == ie .and. elem%face(neigh,j) >0) then
                         if(grid%elem(elem%face(neigh,j))%hsplit == 4 .and. &
                              elem%RGlevel < max_RGlevel) then

                            elem%hsplit = 4
                            imark = imark+1


                            write(5000+state%space%adapt%adapt_level,'(2es12.4,a8,6i5)')&
                                 grid%elem(elem%face(neigh,j))%xc(:),&
                                 ' elem : ',elem%RGlevel, elem%hsplit


                            goto 100
                         endif
                      endif
                   enddo
                endif
             enddo
          endif
100       continue
       enddo
    enddo


    i = itest
    if(i <= grid%nelem) &
         print*,'??? B5',i, grid%elem(i)%hsplit, grid%elem(i)%RGlevel, max_RGlevel

    ! CHECKING
    do i=1,grid%nelem
       elem => grid%elem(i)
       if(elem%RGlevel == max_RGlevel .and. elem%hsplit == 4) then
          print*,'# RGlevel will be exeeded'
          write(*,'(2es12.4,5i5)') elem%xc(:),elem%i, elem%RGlevel, elem%hsplit
          do j=1,elem%flen
             if(elem%face(neigh, j) > 0) then
                elemD => grid%elem(elem%face(neigh, j))
                write(*,'(2es12.4,5i5)') elemD%xc(:),elemD%i,elemD%RGlevel,elemD%hsplit,j
             endif
          enddo
          write(*,'(x)')
          stop
       endif
    enddo

  end subroutine HGMarkingOLD

end module marking



