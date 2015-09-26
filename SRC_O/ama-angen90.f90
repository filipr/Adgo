!> f90 version of the code ANGENER_77
module angener90
  use main_data  ! contains "type(mesh) ::  grid"   for computation
  use AMAdata
  use angener_sub
  use ama_hp_interpol

  implicit none
contains

  !> main subroutine
  subroutine ANGENER_90(metric, ndim)
    implicit none
    logical, intent(in) :: metric
    integer, intent(in) :: ndim
    integer :: i_rem, i_ins, i_swa, ichag, ichag1, nsqrt, ishiftmax
    integer :: noit, imin,  icha, nob, noflc, iter, iter2, iloc, ipocel, ipocel1
    integer :: ichagold, nob1, i, j, k, ipoc, iac_old, iac_new, ii, itet, is, imt
    real :: err, rminerr, rmaxrez, det
    real :: xrmin, xrmax, yrmin, yrmax, surface
    real :: x1, y1, x2, y2, x3, y3
    character(len=15) :: file1, file2
    character(len=5) :: ch5

    real, dimension(:, :), allocatable :: rgabc, rgabcOLD
    real, dimension(:), allocatable :: xi
    integer, dimension(:,:), pointer :: lnd, iae
    real, dimension(:), pointer :: x, y

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lnd => AMA%lnd(1:AMA%melem, 1:3)
    iae => AMA%iae(1:AMA%melem, 1:3)



    write(AMA%ifig1,*)'**************************************************'
    write(AMA%ifig1,*)'****                                          ****'
    write(AMA%ifig1,*)'****       A  N  G  E  N  E  R    3.2         ****'
    write(AMA%ifig1,*)'****                                          ****'
    write(AMA%ifig1,*)'**************************************************'

    noit = 0

    call ADJAC( )

    call COINS( )

    call SEEK_BOUN( )

    call CYKLE( )

    call CYKLE_BOUND( )

    call TEST( )

    if(metric) call METRIX(ndim,  surface)

    !     for back interpolation
    do j=1,ndim
       AMA%wpold(1:AMA%npoin,j) = AMA%wp(1:AMA%npoin,j+1)
       !write(21,*) x(i), y(i), wpold(i, ndim+1)
    enddo


    if(metric)  call ERROR1(ndim, surface )

    ! storing the original metric for refreshing
    allocate(rgabc(1:AMA%mpoin, 1:3), rgabcOLD(1:AMA%mpoin, 1:3))
    rgabcOLD(1:AMA%mpoin, 1) = AMA%rga(1:AMA%mpoin)
    rgabcOLD(1:AMA%mpoin, 2) = AMA%rgb(1:AMA%mpoin)
    rgabcOLD(1:AMA%mpoin, 3) = AMA%rgc(1:AMA%mpoin)


    ! drawing of ellipses
    ! file1 = 'metrixV00000'
    ! is = 0
    ! if(state%space%adapt%adapt_level > 0) is = int(log(1. * state%space%adapt%adapt_level)/log(10.))

    ! write( ch5, '(i5)' ) state%space%adapt%adapt_level  ! change the format if num_size /= 5 !!!
    ! file1(12-is: 12)  = ch5(5-is:5)

    ! imt = 118
    ! open(imt,file=file1,status='unknown')
    ! allocate(xi(1:2) )
    ! do  i=1,AMA%npoin
    !    xi(1) = x(i)
    !    xi(2) = y(i)
    !    call DrawEllips(imt, rgabcOLD(i,1:3), xi(1:2) )
    ! enddo
    ! deallocate(xi)
    ! close(imt)


    ! create a cartesian grid for an interpolation
    call SEEK_FRAME_INTERPOLATION(nsqrt, xrmin, xrmax, yrmin, yrmax, ishiftmax)

    !refreshing of the metric
    call INTERPOLATION(3, rgabcOLD(1:AMA%mpoin, 1:3), rgabc(1:AMA%mpoin, 1:3),  &
         nsqrt,  xrmin, xrmax, yrmin, yrmax, ishiftmax)

    !do i=1,npoin
    !   write(20+state%space%adapt%adapt_level,'(11es12.4)') &
    !        x(i),y(i), rgabcOLD(i, 1:3), rgabc(i, 1:3), rgabcOLD(i, 1:3)- rgabc(i, 1:3)
    !enddo

    !if(state%space%adapt%adapt_level == 2) stop

    !      AMA%ifig2 = 1000
    !      call WriteMetrix(AMA%mpoin, npoin, AMA%melem,  lnd,
    !     *     x, y)

    write(AMA%ifig1,*)'Iteration    quality          changes  ',  &
         ' GL^2:  min        max        nelem'
    write(AMA%ifig1,*)'-----------------------------------------',  &
         '---------------------------------'

    call QUALITY(ndim, err, rminerr,imin,rmaxrez)

    write(AMA%ifig1,99999) noit,AMA%errrez,'begin',0,AMA%glmin,AMA%glmax,AMA%nelem
    ichag = 0

    iter = 1

    call C_REP_BOUND(ndim, noit, iter,  ichag, icha)


    !!print*,'########## ANGENER SKIPPED in ama-angen.f90'
    !!goto 100

    !noflc = 30
    noflc = 60
    do 99 nob=1,noflc

       !refreshing of the metric
       !call INTERPOLATION(3, rgabcOLD(1:AMA%mpoin, 1:3), rgabc(1:AMA%mpoin, 1:3),  &
       !     nsqrt,  xrmin, xrmax, yrmin, yrmax, ishiftmax)

       !weits = 0.5
       !AMA%rga(1:AMA%mpoin) = AMA%rga(1:AMA%mpoin) * weits  + rgabc(1:AMA%mpoin, 1) * ( 1. - weits)
       !AMA%rgb(1:AMA%mpoin) = AMA%rgb(1:AMA%mpoin) * weits  + rgabc(1:AMA%mpoin, 2) * ( 1. - weits)
       !AMA%rgc(1:AMA%mpoin) = AMA%rgc(1:AMA%mpoin) * weits  + rgabc(1:AMA%mpoin, 3) * ( 1. - weits)


       i_rem = 0
       i_ins = 0
       i_swa = 0

       ichag = 0

       ichag1 = ichag
       iter = 5
       call C_REM_BOUND(ndim, noit, iter, ichag, icha)

       i_rem = i_rem + ichag - ichag1
       ichag1 = ichag

       iter = 10
       call C_MOVING(ndim, noit, iter)
       !         goto 100

       i_rem = i_rem + ichag - ichag1
       ichag1 = ichag

       !call TEST( )  ! SMAZ

       iloc = 0
       do iter2 = 1,15
          if(iloc == 0) then
             iter = 1
             call C_REMOVE(ndim, noit, iter, ichag, icha)
             !call TEST( )  ! SMAZ

             i_rem = i_rem + ichag - ichag1
             ichag1 = ichag

             if(icha == 0) iloc = -1
          endif
          !            goto 100

          if(iloc == 0) then
             iter = 5
             call C_SWAPPING(ndim, noit, iter,   ichag, icha)
             !call TEST( )  ! SMAZ

             i_swa = i_swa + ichag - ichag1
             ichag1 = ichag

             iter = 1
             call C_REP_BOUND(ndim, noit, iter,   ichag, icha)
             !call TEST( )  ! SMAZ

          endif
       enddo


       iter = 10
       call C_MOVING(ndim, noit, iter)


       iloc = 0
       do iter2 =1,15
          if(iloc == 0) then

             ichag1 = ichag

             iter = 1
             call C_INSERT(ndim,  noit, iter,   ichag, icha)
             !call TEST( )  ! SMAZ

             i_ins = i_ins + ichag - ichag1
             ichag1 = ichag

             !               if(iter2 == 2)   goto 100

             !     bad dimension
             if(icha == 0) iloc = -1
          endif



          if(iloc == 0) then
             iter = 5
             call C_SWAPPING(ndim, noit, iter,    ichag, icha)
             !call TEST( )  ! SMAZ


             i_swa = i_swa + ichag - ichag1
             ichag1 = ichag

             iter = 1
             call C_REP_BOUND(ndim, noit, iter,    ichag, icha)
             !call TEST( )  ! SMAZ
          endif
       enddo

       iter = 10
       call C_MOVING(ndim, noit, iter)


       ichag1 = ichag

       iter = 5
       call C_SWAPPING(ndim, noit, iter, ichag, icha)
       !call TEST( )  ! SMAZ

       i_swa = i_swa + ichag - ichag1
       ichag1 = ichag

       iter = 1
       call C_REP_BOUND(ndim, noit, iter,   ichag, icha)
       !call TEST( )  ! SMAZ

       !         if(nob == 1) goto 100

       !         call WriteMetrix(AMA%mpoin,  AMA%melem,  lnd,
       !     *        x, y, rga, rgb, rgc)


       write(AMA%ifig1, *)'Global AMA loop',nob, noflc,ichag, ichag1, ichagold
       !write(*, '(a7,i3, 5(a6, i5), a5,es10.2)')  'AMA ope:',nob,  &
       !     'swap:',i_swa, 'ins:', i_ins,   &
       !     'rem:', i_rem,'tot:',ichag, '#T_h:',AMA%nelem, 'qua:',AMA%errrez

       if( ichag == 0 .or. (ichag == ichagold .and. ichag .le. 20) )   &
            goto 101
       ichagold = ichag
99  enddo

101 continue
    !call TEST( )  ! SMAZ

    !     ... no DELAUNAY
    !goto 100
    !     ... to ensure the triangulation of Delaunay type
    do nob1 = 1,3
       ichag = 0

       iter =  30
       call C_MOVING(ndim, noit, iter)

       iter = 1
       if(nob1 == 1) call C_B_INSERT(ndim, noit, iter, ichag, icha)



       iter = 0
       call C_DELAUNAY(ndim, noit, iter, ichag, icha)
       if(ichag == 0) goto 100
    enddo

100 continue

    if(AMA%ifig .ge. 0 ) then
       print*, 'Number of plotted mesh sequence: ',AMA%ifig-1
    endif

    call PLOT( )

    call TEST( )


    !     ... ordering of elements
    call SEEK_NEIGH(AMA%npoin, AMA%nelem, AMA%lnd(1:AMA%nelem, 1:3), &
         AMA%iae(1:AMA%nelem, 1:3),&
         AMA%icyc(1:AMA%mpoin, 1:AMA%maxdeg) )

    !      igra = 11
    !      open (igra,STATUS='unknown', file='matrix_old')
    !      call DRAW_MATRIX(AMA%mbelm,igra)
    !      close(igra)

    !call COMPUTE_BAND(AMA%mbelm,iband)

    !!  minimization of the matrix band
    do i=1,AMA%nelem
       do j=1,3
          AMA%iae1(i,j) = 0
          AMA%lnd1(i,j) = 0
       enddo
       AMA%itrans(i) = 0
    enddo

    ipoc = 1
    iac_old = 1
    do j=1,3
       AMA%lnd1(ipoc,j) = AMA%lnd(iac_old,j)
    enddo
    AMA%iae1(iac_old,1) = -2
    iac_new = ipoc
    AMA%itrans(iac_new) = iac_old

    do i=2,AMA%nelem
       iac_old = AMA%itrans(iac_new)
       if(iac_old == 0) then
          print*,'iac_old is ZERO'
          stop
       endif
       do j=1,3
          ii = iae(iac_old,j)
          if(ii .gt. 0 ) then
             if(AMA%iae1(ii,1) == 0) then
                ipoc = ipoc + 1
                AMA%itrans(ipoc) = ii
                do k=1,3
                   AMA%lnd1(ipoc,k) = AMA%lnd(ii,k)
                enddo
                AMA%iae1(ii,1) = -2
             endif
          endif
       enddo
       iac_new = iac_new+1
    enddo

    call SEEK_NEIGH(AMA%npoin, AMA%nelem, AMA%lnd1(1:AMA%nelem, 1:3), &
         AMA%iae1(1:AMA%nelem, 1:3), &
         AMA%icyc(1:AMA%mpoin, 1:AMA%maxdeg) )

    !call SEEK_NEIGH(npoin,   &
    !     lnd1,iae1)

    !      igra = 11
    !      open (igra,STATUS='unknown', file='matrix_new')
    !      call DRAW_MATRIX(igra,AMA%lnd1,AMA%iae1)
    !      close(igra)

    !call COMPUTE_BAND(AMA%mbelm1,iband1)



    write(AMA%ifig1,*)'The final mesh:'
    write(AMA%ifig1,*)'number of elements:',AMA%nelem
    write(AMA%ifig1,*)'number of points  :',AMA%npoin
    write(AMA%ifig1,*)'number of boun. s.:',AMA%nbelm

    !print*,'AMA: interpolation'
    call INTERPOLATION(ndim+1, AMA%wpold, AMA%wp, nsqrt,  &
         xrmin, xrmax, yrmin, yrmax, ishiftmax)

    !print*,'AMA: recomputation'   !AW always within ADGFEM

    do i=1,AMA%nelem
       do  k=1,ndim + 1  ! including degree of polynomial approximation
          AMA%w(i,k) = sum( AMA%wp(AMA%lnd(i,1:3),k)  )/3
       enddo

       !write(300+state%space%adapt%adapt_level, *)  &
       !     sum( x(AMA%lnd(i,1:3))  )/3, sum( y(AMA%lnd(i,1:3))  )/3, &
       !     AMA%w(i,ndim + 1), int(AMA%w(i,ndim + 1) + 0.5)

    enddo


    !print*,'AMA: heck of "positivity"'
    ipocel = 0
    do 800 i=1,AMA%nelem

       det = POS_TEST(x(lnd(i,1)), y(lnd(i,1)), &
            x(lnd(i,2)), y(lnd(i,2)), x(lnd(i,3)), y(lnd(i,3)) )
       if (det .le. 1.)  then
          ipocel = ipocel + 1
       endif


800 enddo
    ipocel1 = ipocel
    !#print *,'Total number of violations of positivity:',ipocel


    !print*,'AMA: heck of "shaphness"'
    itet = 0
    ipocel = 0
    do  i=1,AMA%nelem
       x1 = x(lnd(i,1))
       y1 = y(lnd(i,1))
       x2 = x(lnd(i,2))
       y2 = y(lnd(i,2))
       x3 = x(lnd(i,3))
       y3 = y(lnd(i,3))
       call POS1TEST(x1,y1,x2,y2,x3,y3,itet)
       if(itet == 1) ipocel = ipocel + 1
       if(itet == 1) then
          write (41,*) x1,y1
          write (41,*) x2,y2
          write (41,*) x3,y3
          write (41,*) x1,y1
          write(41,'(x)')
       endif
    enddo
    print *,'# ANGENER finished, violations of positivity/sharpness: ',ipocel1,ipocel

    deallocate(rgabc, rgabcOLD)

!88888 format(a1,4x,2i5)
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
    return
  end subroutine ANGENER_90


  subroutine  COMPUTE_BAND(iband)
    implicit none
    integer, intent(inout) :: iband
    integer, dimension(:,:), pointer :: iae
    integer:: i,j,ii

    iae => AMA%iae(1:AMA%melem, 1:3)

    iband = 0
    do i=1,AMA%nelem
       do j=1,3
          ii = iae(i,j)
          if(ii .gt. 0) then
             iband = max(iband,abs(i-ii) )
          endif
       enddo
    enddo

    write(*,'(a6,i7,a14,i7,a8,f12.8)')   &
         'nelem=',AMA%nelem,', matrix band=',iband,  &
         ', ratio=',1.*iband/AMA%nelem
  end subroutine COMPUTE_BAND

  subroutine  DRAW_MATRIX(igra)
    implicit none
    integer, intent(inout) :: igra
    integer :: i,j
    integer, dimension(:,:), pointer :: iae

    iae => AMA%iae(1:AMA%melem, 1:3)

    do i=1,AMA%nelem
       write(igra,*) i, -i
       do j=1,3
          if(iae(i,j) .gt. 0) then
             write(igra,*) i, -iae(i,j)
          endif
       enddo
    enddo
  end subroutine DRAW_MATRIX



  subroutine SEEK_FRAME_INTERPOLATION(nsqrt, xrmin, xrmax, yrmin, yrmax, ishiftmax)
    implicit none
    real, dimension(:), pointer :: xold, yold
    integer, dimension(:, :), pointer ::  lndold
    !integer, dimension(1:AMA%melem, 1:4), intent(inout) ::  iaegr
    integer, intent(inout) ::  nsqrt, ishiftmax
    real, intent(inout) :: xrmin, xrmax, yrmin, yrmax
    integer :: i, ix, iy, iele, iel, ielem, iord, it, isum, inum, icom, ic, ip, il
    real :: xc, yc, sradiusmax, x1, y1, x2, y2, x3, y3, radius, sx, sy, rboxmax
    real :: tol_shift = 1E-05, rshift
    !print*,'AMA: coefficients for transformation'

    lndold =>  AMA%lndold(1:AMA%melem, 1:3)
    xold => AMA%xold( 1:AMA%mpoin )
    yold => AMA%yold( 1:AMA%mpoin )

    xrmax = -1E+30
    xrmin = -xrmax
    yrmax = xrmax
    yrmin = -yrmax
    do i=1,AMA%npoinold
       xc = xold(i)
       yc = yold(i)
       xrmax = max(xrmax,xc)
       xrmin = min(xrmin,xc)
       yrmax = max(yrmax,yc)
       yrmin = min(yrmin,yc)
    enddo

    ! slightly larger frame
    rshift = tol_shift * (xrmax - xrmin)
    xrmax = xrmax + rshift
    xrmin = xrmin - rshift

    rshift = tol_shift * (yrmax - yrmin)
    yrmax = yrmax + rshift
    yrmin = yrmin - rshift

    sradiusmax = 0.
    do i=1,AMA%nelemold
       AMA%iaegr(i,2) =  0
       AMA%iaegr(i,3) =  0
       x1 = xold(lndold(i,1))
       y1 = yold(lndold(i,1))
       x2 = xold(lndold(i,2))
       y2 = yold(lndold(i,2))
       x3 = xold(lndold(i,3))
       y3 = yold(lndold(i,3))
       radius = max(((x1-x2)**2 + (y1-y2)**2)**0.5,    &
            ((x2-x3)**2 + (y2-y3)**2)**0.5 ,  &
            ((x1-x3)**2 + (y1-y3)**2)**0.5, sradiusmax )
       sradiusmax = max(sradiusmax, radius)
    enddo
    sradiusmax = sradiusmax *2

    !write(*,'(a25,6es12.4)') 'AMA-solve: sradiusmax =', &
    !     sradiusmax, xrmin, xrmax, yrmin, yrmax

    nsqrt = AMA%nelemold**0.5
    do i=1,AMA%nelemold
       xc = (xold(lndold(i,1)) + xold(lndold(i,2)) +   &
            xold(lndold(i,3)) )/3
       yc = (yold(lndold(i,1)) + yold(lndold(i,2)) +   &
            yold(lndold(i,3)) )/3
       sx = (xc - xrmin)/(xrmax - xrmin)
       sy = (yc - yrmin)/(yrmax - yrmin)
       if(sx .ge. 1. .or. sy .ge. 1 .or.   &
            sx .le. 0 .or. sy .le. 0.) then
          print*, 'Trouble in INTERPOLATION'
          print *,xc,yc,sx,sy
          stop
       endif
       ix = int(sx*nsqrt) + 1
       iy = int(sy*nsqrt) + 1
       iele = (ix-1)*nsqrt + iy
       AMA%iaegr(i,1) = iele
       AMA%iaegr(iele,2) = AMA%iaegr(iele,2) + 1

       if(state%space%adapt%adapt_level  == -2) then
          write(*,'(a6,8i5)') 'FR:',i, iele, ix, iy, AMA%iaegr(iele,2)
       endif

    enddo
    rboxmax = max((xrmax - xrmin)/nsqrt,(yrmax - yrmin)/nsqrt )
    ishiftmax = sradiusmax/rboxmax + 1

    if(state%space%adapt%adapt_level  == -2) then
       print *,'BOXES:',nsqrt,nsqrt**2,AMA%nelemold
       print *,sradiusmax,rboxmax,sradiusmax/rboxmax,ishiftmax
       it = 14
       do i=1,AMA%nelemold
          if(AMA%iaegr(i,1) == it) then
             write(31,*)xold(lndold(i,1)), yold(lndold(i,1)),lndold(i,1),i
             write(31,*)xold(lndold(i,2)), yold(lndold(i,2)),lndold(i,2)
             write(31,*)xold(lndold(i,3)), yold(lndold(i,3)),lndold(i,3)
             write(31,*)xold(lndold(i,1)), yold(lndold(i,1))
             write(31,*)
          endif
       enddo
       print *,'pocet =',AMA%iaegr(it,2),it/nsqrt,mod(it,nsqrt)
    endif

    if(state%space%adapt%adapt_level  == -2) then
       write(*,'(a6,50i5)') 'iaegr2',AMA%iaegr(1:15,2)
    endif

    isum = 0
    do i=1,nsqrt**2
       inum = AMA%iaegr(i,2)
       AMA%iaegr(i,2) = isum + 1
       isum = isum + inum
    enddo

    if(state%space%adapt%adapt_level  == -2) then
       write(*,'(a6,50i5)') 'iaegr2',AMA%iaegr(1:15,2)
       write(*,'(a6,50i5)') 'iaegr1',AMA%iaegr(1:15,1)
    endif

    do i=1,AMA%nelemold
       icom = AMA%iaegr(i,1)
       iord = AMA%iaegr(icom,2) + AMA%iaegr(icom,3)
       AMA%iaegr(icom,3) = AMA%iaegr(icom,3) + 1
       AMA%iaegr(iord,4) = i
    enddo

    if(state%space%adapt%adapt_level  == -2) then
       write(*,'(a6,50i5)') 'iaegr3',AMA%iaegr(1:15,3)
       write(*,'(a6,50i5)') 'iaegr4',AMA%iaegr(1:15,4)
    endif

    do ic =1,nsqrt**2
       ip = AMA%iaegr(ic,2)
       if(ic .ne. nsqrt**2) then
          il = AMA%iaegr(ic+1,2) - 1
       else
          il = AMA%nelemold
       endif
       if(ip .gt. il .and. AMA%iaegr(ic,3) .ne. 0) then
          print *,'error for ic =',ic
          print *,'nozero'
          print *,ip,il,nsqrt**2
          stop
       endif
       do iel = ip,il
          ielem = AMA%iaegr(iel,4)
          if(AMA%iaegr(ielem,1) .ne. ic) then
             print *,'error for ic =',ic
             print *,'the elements don''t correspond'
             print *,ic,iel,ielem,AMA%iaegr(ielem,1)
             stop
          endif
       enddo
    enddo


    !if(state%space%adapt%adapt_level  == 6) then
    !   do i=1,AMA%nelem
    !      write(*,'(a6,50i5)') 'iaegr3',i, AMA%iaegr(i,1:4)
    !   enddo
    !endif


    ! frame found
  end subroutine SEEK_FRAME_INTERPOLATION

  !> interpolation of the solution from the old mesh to the new one
  subroutine INTERPOLATION(ndimL, wpold, wp, &
       nsqrt, xrmin, xrmax, yrmin, yrmax, ishiftmax)
    implicit none
    integer, intent(in):: ndimL
    real, dimension(:), pointer :: xold, yold
    integer, dimension(:, :), pointer ::   lndold
    real, dimension(1:AMA%mpoin, 1:ndimL ), intent(in) :: wpold
    !integer, dimension(1:AMA%melem, 1:4), intent(in) ::  iaegr
    real, dimension(1:AMA%mpoin, 1:ndimL ), intent(inout) :: wp
    integer, intent(in) ::  nsqrt,  ishiftmax
    real, intent(in) :: xrmin, xrmax, yrmin, yrmax
    integer, dimension(:), allocatable :: itli
    real, dimension(:,:), allocatable :: tlr
    integer :: ip, ishift, iposit, ipom, iiner, ix, iy, lx, ly, ic, il1, il2
    integer :: iel, i, i1, i2, i3, ie, imin, k, itest
    real :: x0, y0, xc, yc, det0, det1, det2, det3, sx, sy, rmin
    real :: xx0, yy0, xx1, yy1, xx2,yy2, xx3, yy3, epss, yyc, rlen
    integer, dimension(:,:), pointer :: lnd
    real, dimension(:), pointer :: x, y

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lnd => AMA%lnd(1:AMA%melem, 1:3)

    itest = -3

    lndold =>  AMA%lndold(1:AMA%melem, 1:3)
    xold => AMA%xold( 1:AMA%mpoin )
    yold => AMA%yold( 1:AMA%mpoin )

    allocate(itli(1:AMA%npoin), tlr(1:AMA%npoin, 1:3) )
    itli(:) = 0
    tlr(:,:) = 0.

    do ip=1,AMA%npoin
       !         if(mod(Ip,1000) == 1) print *,'***',ip,ishift
       ishift = 0
       x0 = x(ip)
       y0 = y(ip)
       iposit = 0

       !     some improvement for a unit circle with a crack
       if(AMA%icrack == 1) then
          if(y0 == 0) then
             do ipom = 1,AMA%nelem
                if(lnd(ipom,1) == ip .or. lnd(ipom,2) == ip  &
                     .or. lnd(ipom,3) == ip ) then
                   yc = (y(lnd(ipom,1))+y(lnd(ipom,2))+  &
                        y(lnd(ipom,3)))/3
                   if(yc .gt. 0) then
                      iposit = 1
                   else
                      iposit = -1
                   endif
                   goto 234
                endif
             enddo
234          continue
          endif
       endif


       sx = (x0 - xrmin)/(xrmax - xrmin)
       sy = (y0 - yrmin)/(yrmax - yrmin)
       iiner = 0
       if(sx .ge. 0. .and. sx .le. 1. .and.   &
            sy .ge. 0. .and. sy .le. 1.) then
          !     inner point
          iiner = 1
       endif

       if(itest == ip) then
          write(*,'(a8,3i5, 2es12.4)') 'itest', ishift, iiner, ip, sx, sy
       endif

       if(iiner /= 1) then
          print* ,'Trouble in INTERPOLATION ?????'
          print*,sx,sy, x0,y0
          stop

       else
          ix = int(sx*nsqrt) + 1
          iy = int(sy*nsqrt) + 1

          if(itest == ip) then
             write(32, *) x0, y0
             write(*,'(a8,3i5, 2es12.4)') 'itest', ishift, ix, iy, sx, sy
          endif

954       continue
          do lx = -ishift,ishift
             do ly = -ishift,ishift
                if( ix+lx .ge. 1 .and. ix+lx .le. nsqrt .and.  &
                     iy+ly .ge. 1 .and. iy+ly .le. nsqrt .and.  &
                     (abs(lx) == ishift .or. abs(ly) == ishift) )  &
                     then

                   ic = (ix-1+lx)*nsqrt + iy+ly
                   il1 = AMA%iaegr(ic,2)
                   if(ic .ne. nsqrt**2) then
                      il2 = AMA%iaegr(ic+1,2) - 1
                   else
                      il2 = AMA%nelemold
                   endif
                   !     we go through the elements belonging to ic

                   !print*,'$$$',il1,il2,AMA%nelemold, ic
                   do iel = il1,il2
                      i = AMA%iaegr(iel,4)
                      i1 = lndold(i,1)
                      i2 = lndold(i,2)
                      i3 = lndold(i,3)
                      xx0 = x0
                      yy0 = y0
                      xx1 = xold(i1)
                      yy1 = yold(i1)
                      xx2 = xold(i2)
                      yy2 = yold(i2)
                      xx3 = xold(i3)
                      yy3 = yold(i3)

                      if(iposit .ne. 0) then
                         yyc = (yy1 + yy2 + yy3)/3.
                         if(yyc*iposit .lt. 0) then
                            goto 235
                         endif
                      endif

                      det0 = xx3*(yy1-yy2)+xx1*(yy2-yy3)+xx2*(yy3-yy1)
                      det1 = xx0*(yy1-yy2)+xx1*(yy2-yy0)+xx2*(yy0-yy1)
                      det2 = xx0*(yy2-yy3)+xx2*(yy3-yy0)+xx3*(yy0-yy2)
                      det3 = xx0*(yy3-yy1)+xx3*(yy1-yy0)+xx1*(yy0-yy3)


                      if(det0 .le. 0.) then
                         write(112,'(2e14.6,2i7)')xx1,yy1,i1,i
                         write(112,'(2e14.6,i7)')xx2,yy2,i2
                         write(112,'(2e14.6,i7)')xx3,yy3,i3
                         write(112,'(3e14.6)')xold(i1),yold(i1),det0
                         write(112,'(x)')
                      endif
                      epss = 1E-03
                      if(abs(det1)/det0 .lt. epss) det1 = abs(det1)
                      if(abs(det2)/det0 .lt. epss) det2 = abs(det2)
                      if(abs(det3)/det0 .lt. epss) det3 = abs(det3)
                      if(det1 .ge. 0 .and. det2 .ge.0   &
                           .and. det3 .ge. 0 .and. det0 .gt. 0.) then
                         !     the node is in this triangle
                         itli(ip) = i
                         tlr(ip,1) = det2/det0
                         tlr(ip,2) = det3/det0
                         tlr(ip,3) = det1/det0
                         goto 209
                      endif
235                   continue
                   enddo
                endif
             enddo
          enddo

          ishift = ishift + 1
          if(ishift .le. ishiftmax) goto 954
          iiner = 0
209       continue
       endif

       if(iiner == 0) then
          !     outer point
          rmin = 100000.
          do ie =1,AMA%nelemold
             xc = (xold(lndold(ie,1)) + xold(lndold(ie,2)) +   &
                  xold(lndold(ie,3)))/3
             yc = (yold(lndold(ie,1)) + yold(lndold(ie,2)) +   &
                  yold(lndold(ie,3)))/3
             rlen = ((x0-xc)**2 + (y0-yc)**2 )**0.5
             if(iposit .ne. 0) then
                if(yc*iposit .lt. 0) then
                   goto 236
                endif
             endif

             if(rlen .lt. rmin) then
                rmin = rlen
                imin = ie
             endif
236          continue
          enddo
          i = imin
          i1 = lndold(i,1)
          i2 = lndold(i,2)
          i3 = lndold(i,3)
          xx0 = x0
          yy0 = y0
          xx1 = xold(i1)
          yy1 = yold(i1)
          xx2 = xold(i2)
          yy2 = yold(i2)
          xx3 = xold(i3)
          yy3 = yold(i3)

          det0 = xx3*(yy1-yy2)+xx1*(yy2-yy3)+xx2*(yy3-yy1)
          det1 = xx0*(yy1-yy2)+xx1*(yy2-yy0)+xx2*(yy0-yy1)
          det2 = xx0*(yy2-yy3)+xx2*(yy3-yy0)+xx3*(yy0-yy2)
          det3 = xx0*(yy3-yy1)+xx3*(yy1-yy0)+xx1*(yy0-yy3)
          !     the node is the nearest to this triangle
          itli(ip) = i
          tlr(ip,1) = det2/det0
          tlr(ip,2) = det3/det0
          tlr(ip,3) = det1/det0
          if(abs(tlr(ip,1)).gt. 1.5 .or. abs(tlr(ip,2)).gt. 1.5  &
               .or. abs(tlr(ip,3) ) .gt. 1.5 ) then
             tlr(ip,1) = 1./3
             tlr(ip,2) = 1./3
             tlr(ip,3) = 1./3
          endif
       endif
    enddo



    do i=1,AMA%npoin
       ie = itli(i)
       do k=1,ndimL   ! including degree of polynomial approximation
          wp(i,k) = tlr(i,1)*wpold(lndold(ie,1),k) +  &
               tlr(i,2)*wpold(lndold(ie,2),k) +  &
               tlr(i,3)*wpold(lndold(ie,3),k)
       enddo

       !         wp(i,1) = tlr(i,1)*wpold(lndold(ie,1),ndim+1) +
       !     *        tlr(i,2)*wpold(lndold(ie,2),ndim+1) +
       !     *        tlr(i,3)*wpold(lndold(ie,3),ndim+1)

    enddo

    !if(ndimL /= 3) then
    !   do i=1,AMA%npoinold
    !      ip = int(wpold(i, ndimL) + 0.5)
    !      write(2000+state%space%adapt%adapt_level*10 + ip, *) xold(i), yold(i), wp(i, ndimL)
    !   enddo
    !
    !   do i=1,AMA%npoin
    !      ip = int(wp(i, ndimL) + 0.5)
    !      write(3000+state%space%adapt%adapt_level*10 + ip, *) x(i), y(i), wp(i, ndimL)
    !   enddo
    !end if

    deallocate(itli, tlr )

  end subroutine INTERPOLATION


  subroutine C_DELAUNAY(ndim,  noit, iter,  ichag, icha)
    implicit none
    integer, intent(in) :: ndim, iter
    integer, intent(inout) :: noit, ichag, icha
    real:: err, err1, err1rez, rminerr, rmaxrez
    integer :: icy, imin

    err1 = 0.
    err1rez = 0.
    do 30 icy =1,iter
       noit = noit + 1
       call DELAUNAY(icha)

       call CYKLE_REP( )

       call CYKLE_BOUND()

       call QUALITY(ndim, err, rminerr, imin, rmaxrez)

       write(AMA%ifig1,99999)noit,AMA%errrez,'delaunay',icha, AMA%glmin, AMA%glmax, &
            AMA%nelem
       ichag = ichag + icha
       if(icha .le. 1 ) return
30  enddo
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
    return
  end subroutine C_DELAUNAY

  subroutine C_SWAPPING(ndim, noit, iter, ichag, icha)
    implicit none
    integer, intent(in) :: ndim, iter
    integer, intent(inout) :: noit, ichag, icha
    real:: err, err1, err1rez, rminerr, rmaxrez
    integer :: icy, imin
    err1 = 0.
    err1rez = 0.
    do  icy =1,iter
       noit = noit + 1

       call SWAPPING(icha, icy)

       call CYKLE_REP( )

       call CYKLE_BOUND()

       call QUALITY(ndim, err, rminerr, imin, rmaxrez )

       write(AMA%ifig1,99999)noit,AMA%errrez,'swapping',icha,AMA%glmin,AMA%glmax, &
            AMA%nelem
       ichag = ichag + icha
       if(icha .le. 1 ) return
    enddo
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
    return
  end subroutine C_SWAPPING

  subroutine C_MOVING(ndim, noit, iter)
    implicit none
    integer, intent(in) :: ndim, iter
    integer, intent(inout) :: noit
    real:: err, err1, err1rez, rminerr, rmaxrez
    integer :: icy, imin

    err1 = 0.
    err1rez = 1.E+35

    do 40 icy =1,iter
       noit = noit + 1
       call MOVING(icy )

       call QUALITY(ndim, err, rminerr, imin, rmaxrez)

       write(AMA%ifig1,99998) noit,AMA%errrez,'moving',AMA%glmin,AMA%glmax, AMA%nelem
       if( err1rez - AMA%errrez .lt. 1E-06 .or. icy == iter) then
          write(AMA%ifig1,99998) noit,AMA%errrez,'moving',AMA%glmin,AMA%glmax, &
               AMA%nelem
          return
       endif

       err1rez = AMA%errrez
40  enddo
99998 format (i6,2x,1(e16.8),2x,a9,5x,2x,2e12.4,i8)
    return
  end subroutine C_MOVING

  subroutine C_REM_BOUND(ndim, noit, iter, ichag, icha)
    implicit none
    integer, intent(in) :: ndim, iter
    integer, intent(inout) :: noit, ichag, icha
    real:: err, rminerr, rmaxrez
    integer :: it, imin

    do 20 it =1,iter
       noit = noit + 1
       call REM_BOUND(ndim, icha, it)

       call QUALITY(ndim, err, rminerr, imin, rmaxrez )

       write(AMA%ifig1,99999)noit,AMA%errrez,'rem_bou',icha,AMA%glmin,AMA%glmax, &
            AMA%nelem
       ichag = ichag + icha
       if(icha .le. 1) return
20  enddo
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
    return
  end subroutine C_REM_BOUND

  subroutine C_REMOVE(ndim, noit, iter, ichag, icha)
    implicit none
    integer, intent(in) :: ndim, iter
    integer, intent(inout) :: noit, ichag, icha
    real:: err, rminerr, rmaxrez
    integer :: icy, imin

    do 30 icy =1,iter
       noit = noit + 1
       call REMOVE(ndim, icha, icy)

       call QUALITY(ndim, err, rminerr, imin, rmaxrez )

       write(AMA%ifig1,99999)noit,AMA%errrez,'remove',icha,AMA%glmin,AMA%glmax,&
            AMA%nelem
       ichag = ichag + icha
       if(icha .le. 1) return
30  enddo
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
    return
  end subroutine C_REMOVE

  subroutine C_INSERT(ndim, noit, iter, ichag, icha)
    implicit none
    integer, intent(in) :: ndim, iter
    integer, intent(inout) :: noit, ichag, icha
    real:: err, rminerr, rmaxrez
    integer :: icy, imin

    do 30 icy =1,iter
       noit = noit + 1
       call INSERT(ndim, icha, icy)

       call CYKLE_REP( )

       call CYKLE_BOUND( )

       call QUALITY(ndim, err, rminerr, imin, rmaxrez )

       write(AMA%ifig1,99999)noit,AMA%errrez,'insert',icha,AMA%glmin,AMA%glmax,AMA%nelem
       ichag = ichag + icha
       if(icha .le. 1) return
30  enddo
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
    return
  end subroutine C_INSERT

  subroutine C_REP_BOUND(ndim, noit, iter, ichag, icha)
    implicit none
    integer, intent(in) :: ndim, iter
    integer, intent(inout) :: noit, ichag, icha
    real:: err, rminerr, rmaxrez
    integer :: icy, imin

    do 30 icy =1,iter
       noit = noit + 1

       call REP_BOUND( icha, icy )


       call QUALITY(ndim, err, rminerr, imin, rmaxrez )

       write(AMA%ifig1,99999)noit,AMA%errrez,'rep_boun',icha,AMA%glmin,AMA%glmax,&
            AMA%nelem
       ichag = ichag + icha
30  enddo
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
    return
  end subroutine C_REP_BOUND


  subroutine C_B_INSERT(ndim, noit, iter, ichag, icha)
    implicit none
    integer, intent(in) :: ndim, iter
    integer, intent(inout) :: noit, ichag, icha
    real:: err, rminerr, rmaxrez
    integer :: icy, imin

    do 30 icy =1,iter
       noit = noit + 1
       call INSERT_BOUND(ndim, icha, icy)

       call CYKLE_REP(  )

       call CYKLE_BOUND()

       call QUALITY(ndim, err, rminerr, imin, rmaxrez)

       write(AMA%ifig1,99999)noit,AMA%errrez,'B_insert',icha,AMA%glmin,AMA%glmax,&
            AMA%nelem
       ichag = ichag + icha
       if(icha .le. 1) return
30  enddo
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
    return
  end subroutine C_B_INSERT


  subroutine REP_BOUND( icha, icy )
    implicit none
    integer, intent(in) :: icy
    integer, intent(inout) :: icha
    integer :: iyi, i, j, j1, j2, kk, k0, k1, k2, irtk, il0, il1, il2,  irtkdel
    integer :: ill1, ill2, i1, i2, jj1, jj2, kl1, kl2, l, k, kk1, ll, ll1
    integer :: ib, ib1, ikk1, ipoc, itet
    integer lock(20)
    real :: x0, y0, x1, y1, x2, y2, det, rll, rl0, rlen0
    real :: rd, det123, xx1, yy1, xx2, yy2
    integer:: nelemold, nbelmold ! local variables
    real, dimension(:), pointer :: xb, yb
    integer, dimension(:,:), pointer :: ibb, ibp
    integer, dimension(:,:), pointer :: icyc
    integer, dimension(:,:), pointer :: lnd, iae, lbn
    integer, dimension(:), pointer :: ibc, itc
    real, dimension(:), pointer :: x, y
    integer :: ic_start, ic_end, ic_skip

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lbn => AMA%lbn(1:AMA%mbelm,1:2)
    ibc => AMA%ibc(1:AMA%mbelm)
    itc => AMA%itc(1:AMA%mbelm)

    lnd => AMA%lnd(1:AMA%melem, 1:3)
    iae => AMA%iae(1:AMA%melem, 1:3)

    icyc => AMA%icyc(1:AMA%mpoin, 1:AMA%maxdeg)
    ibp => AMA%ibp(1:AMA%mpoin, 1:2)
    ibb => AMA%ibb(1:AMA%mpoin, 1:3)

    xb => AMA%xb(1:AMA%ipoint)
    yb => AMA%yb(1:AMA%ipoint)

    icha = 0
    nelemold = AMA%nelem
    ipoc = 1

    if(mod(icy, 2) == 0) then
       ic_start = 1;         ic_end = nelemold;   ic_skip =  1
    else
       ic_start = nelemold;  ic_end = 1;          ic_skip = -1
    endif

!    do 10 iyi=1,nelemold
    do 10 iyi = ic_start, ic_end, ic_skip

       i = ipoc
       do 20 j=1,3
          if(iae(i,j) .lt. 0) then
             j1 = mod(j,3) + 1
             j2 = mod(j1,3) + 1
             k1 = lnd(i,j)
             if(ibb(k1,1) .gt. 0) then
                ikk1 = ibb(k1,2)
                if(ibb(k1,1) == AMA%ibpoin(ikk1-1)+1  .or.  &
                     ibb(k1,1) == AMA%ibpoin(ikk1)) then
                   !               NO moving of final or initial node of profiles
                   goto 21
                endif
             endif
             k2 = lnd(i,j1)
             k0 = lnd(i,j2)
             x1 = x(k1)
             y1 = y(k1)
             x2 = x(k2)
             y2 = y(k2)
             x0 = x(k0)
             y0 = y(k0)


             if(ibb(k1,1) .gt. 0 .and. ibb(k2,1) .gt. 0 )then

                if(ibb(k1,2)  .ne. ibb(k2,2) ) then
                   print *,'Troubles in insert.f'
                   print *,ibb(k1,2), ibb(k2,2)
                endif
                irtk = ibb(k1,2)

                if(k0 == -1289) then
                   print *,x0,y0,k0
                   print *,x1,y1,k1
                   print *,x2,y2,k2
                   print *,'@@@@@'
                endif

                !     we seek the point in the field [xb(i),yb(i)] i=1,ipoint
                !     ... rll length of edge k1,k2
                rll = ((x(k1) - x(k2))*(x(k1) - x(k2)) +  &
                     (y(k1) - y(k2))*(y(k1) - y(k2)) )**0.5
                rl0 = 1E+25*rll
                do 30 ll1=AMA%ibpoin(irtk-1)+1,AMA%ibpoin(irtk)
                   rlen0 = (x0 -xb(ll1))*(x0 -xb(ll1)) +   &
                        (y0 -yb(ll1))*(y0 -yb(ll1))
                   if(rlen0 .lt. rl0) then
                      rl0 = rlen0
                      il0 = ll1
                   endif
30              enddo
                il1 = ibb(k1,1)
                il2 = ibb(k2,1)


                if(il1 == il0 .or. il2 == il0) then
                   goto 20
                endif



                !     to verify: il0 must be between il1 and il2

                !                  if( ( il0 .gt. il1 .and. il0 .gt. il2) .or.
                !     *                 ( il0 .lt. il1 .and. il0 .lt. il2)) goto 20

                irtkdel = abs(AMA%ibpoin(irtk-1)+1 - AMA%ibpoin(irtk))
                ill1 = il1 - il0
                if(abs(ill1) .gt. irtkdel/2) then
                   if(il1 .gt. il0) then
                      ill1 = ill1 - irtkdel
                   else
                      ill1 = ill1 + irtkdel
                   endif
                   !                     print *,'@@@',il1,il0,il2,ill1,ill2
                endif
                ill2 = il2 - il0
                if(abs(ill2) .gt. irtkdel/2) then
                   if(il2 .gt. il0) then
                      ill2 = ill2 - irtkdel
                   else
                      ill2 = ill2 + irtkdel
                   endif
                   !                     print *,'@.@',il1,il0,il2,ill1,ill2
                endif

                if( ill1 * ill2  .gt. 0 ) goto 20
                !     il0 must be between il1 and il2

                !     ... computation of the length of il0 from the edge k1,k2
                x0 = xb(il0)
                y0 = yb(il0)
                det = x0*(y1-y2) + x1*(y2-y0) + x2*(y0-y1)
                rd = det/rll

                if(k0 == -1289) then
                   print *,'!!',rl0**0.5,rd,0.45*rd,icyc(k0,1)
                endif

                if(rl0**0.5 .lt. 0.95*rd .and. icyc(k0,1) .gt. 0) then
                   !                     print *,'!!!',k0,rl0**0.5,rd,rd1
                   !                     print *,x0,y0,k0,il0
                   !                     print *,x1,y1,k1,il1
                   !                     print *,x2,y2,k2,il2
                   !                     print *,xb(il0),yb(il0)
                   !                     print *

                   !     this is a candidate, we check the positivity
                   do 40 kk=1,icyc(k0,1)
                      if(icyc(k0,kk+1) .ne. k1) then
                         kk1 = mod(kk,icyc(k0,1)) + 1
                         xx1 = x(icyc(k0,kk+1))
                         yy1 = y(icyc(k0,kk+1))
                         xx2 = x(icyc(k0,kk1+1))
                         yy2 = y(icyc(k0,kk1+1))

                         det123 = POS_TEST(x0, y0, xx1, yy1, xx2, yy2)
                         !print *,'..RP .',det123
                         if( det123 .le. 1.) then
                            !     violation of positivity, go to next i
                            goto 20
                         endif

                         call POS2TEST(x0,y0,xx1,yy1,xx2,yy2,itet)
                         if(itet == 1) then
                            !     violation of positivity, go to next i
                            goto 20
                         endif

                      endif
40                 enddo

                   !     we remove this triangle
                   if(icyc(k1,1) .gt. 0 .or. icyc(k2,1) .gt. 0) then
                      print *,'very divny in REP_BOUN'
                      stop
                   endif

                   !                     print *,'@@@'
                   !                     print *, x(k0),y(k0)
                   !                     print *,xb(il0),yb(il0)


                   x(k0) = xb(il0)
                   y(k0) = yb(il0)
                   ibb(k0,1) = il0
                   do ll=1,AMA%nbp
                      if(il0 .ge. AMA%ibpoin(ll-1) +1 .and.  &
                           il0 .le. AMA%ibpoin(ll) ) then
                         ibb(k0,2) = ll
                         goto 109
                      endif
                   enddo
109                continue


                   i1 = iae(i,j1)
                   i2 = iae(i,j2)

                   if(i1 .lt. 0) then
                      print *,'ERROR1 in REP_BOUN'
                   else
                      jj1 = 0
                      do ll=1,3
                         if(iae(i1,ll) == i) jj1 = ll
                      enddo
                   endif
                   if(i2 .lt. 0) then
                      print *,'ERROR2 in REP_BOUN'
                   else
                      jj2 = 0
                      do ll=1,3
                         if(iae(i2,ll) == i) jj2 = ll
                      enddo
                   endif
                   if(jj1 == 0) then
                      print *,'ERROR3 in REP_BOUN'
                   else
                      iae(i1,jj1) = -2
                   endif
                   if(jj2 == 0) then
                      print *,'ERROR4 in REP_BOUN'
                   else
                      iae(i2,jj2) = -2
                   endif

                   do 101 l=1,AMA%nelem
                      do 111 k=1,3
                         if(iae(l,k) .gt. i) iae(l,k) = iae(l,k)-1
111                   enddo
101                enddo
                   do 100 l=i,AMA%nelem-1
                      do 110 k=1,3
                         lnd(l,k) = lnd(l+1,k)
                         iae(l,k) = iae(l+1,k)
110                   enddo
100                enddo

                   if(icyc(k1,2) == k2) then
                      do l=1,abs(icyc(k1,1))-1
                         icyc(k1,l+1) = icyc(k1,l+2)
                      enddo
                      icyc(k1,1) = icyc(k1,1) + 1
                   elseif(icyc(k1,abs(icyc(k1,1))+1) == k2)then
                      icyc(k1,1) = icyc(k1,1) + 1
                   else
                      print *,'ERROR k2 is not in the cykle of k1'
                      stop
                   endif
                   if(icyc(k2,2) == k1) then
                      do l=1,abs(icyc(k2,1))-1
                         icyc(k2,l+1) = icyc(k2,l+2)
                      enddo
                      icyc(k2,1) = icyc(k2,1) + 1
                   elseif(icyc(k2,abs(icyc(k2,1))+1) == k1)then
                      icyc(k2,1) = icyc(k2,1) + 1
                   else
                      print *,'ERROR k1 is not in the cykle of k2'
                      stop
                   endif

                   kl1 = 0
                   kl2 = 0
                   do l=1,icyc(k0,1)
                      lock(l) = icyc(k0,l+1)
                      if(lock(l) == k1) kl1 = l
                      if(lock(l) == k2) kl2 = l
                   enddo
                   if(kl1*kl2 == 0) then
                      print *,'k1 and k2 are not in cykle of k0'
                      stop
                   endif

                   if(abs(kl1-kl2) == 1) then
                      if(kl1 .lt. kl2) then
                         do l=1, icyc(k0,1) - kl1
                            icyc(k0,l+1) = lock(kl1+l)
                         enddo
                         do l=1,kl1
                            icyc(k0,icyc(k0,1)-kl1+l+1) = lock(l)
                         enddo
                      else
                         do l=1, icyc(k0,1) - kl2
                            icyc(k0,l+1) = lock(kl2+l)
                         enddo
                         do l=1,kl2
                            icyc(k0,icyc(k0,1)-kl2+l+1) = lock(l)
                         enddo
                      endif
                      icyc(k0,1) = -icyc(k0,1)
                   else
                      !     all is OK
                      icyc(k0,1) = -icyc(k0,1)
                   endif

                   !     repair of boundary files
                   do ib=1,AMA%nbelm
                      if(itc(ib) .gt. i) itc(ib) = itc(ib) - 1
                   enddo

                   nbelmold = AMA%nbelm
                   do 600 ib=1,nbelmold
                      if(lbn(ib,1) == k1   &
                           .and. lbn(ib,2) == k2) then
                         do ib1=0,AMA%nbelm-ib-1
                            lbn(AMA%nbelm+1-ib1,1) = lbn(AMA%nbelm-ib1,1)
                            lbn(AMA%nbelm+1-ib1,2) = lbn(AMA%nbelm-ib1,2)
                            ibc(AMA%nbelm+1-ib1) = ibc(AMA%nbelm-ib1)
                            itc(AMA%nbelm+1-ib1) = itc(AMA%nbelm-ib1)
                         enddo
                         lbn(ib,2) = k0
                         lbn(ib+1,1) = k0
                         lbn(ib+1,2) = k2
                         ibc(ib+1) = ibc(ib)
                         if(i2.gt.i) then
                            itc(ib) = i2-1
                         else
                            itc(ib) = i2
                         endif
                         if(i1.gt.i) then
                            itc(ib+1) = i1-1
                         else
                            itc(ib+1) = i1
                         endif
                         AMA%nbelm = AMA%nbelm + 1
                         goto 610
                      endif
600                enddo
                   print *,'boundary segment doesn''t found'
                   stop
610                continue
                   AMA%nelem = AMA%nelem - 1
                   icha = icha + 1

                   if(AMA%ifig .ge. 0 .and. mod(icha,5) == 0 ) then
                      call PLOT1()
                      AMA%ifig = AMA%ifig + 1
                   endif

                endif
             endif
          endif
21        continue
20     enddo
       if( i .gt. AMA%nelem) goto 2000
       ipoc = ipoc + 1
10  enddo
2000 continue

    return
  end subroutine REP_BOUND

  subroutine INSERT(ndim, icha, icy)
    implicit none
    integer, intent(in) :: ndim, icy
    integer, intent(inout) :: icha
    real xx(4),yy(4),rmax(3)
    integer jmax(3)
    integer:: nelemold ! local variable
    integer :: ice, itet, itest, ipoc, i, j, j1, j2, ii1, ii2, k, l, ll, ib, ib1
    integer :: ii, jj1, jj2, k1, k2, k3, k4, kk, kl1, ia1, ia2, il0, ipe
    integer :: il1, il2, il2new, ll1, ll11, ibper, je1, ke1, ke2, ke3
    integer :: imov, ja1, ja2, jmaxhelp, kl, jja1, jja2, jbak, je2, jel, jel1, jj, jjj
    integer :: idif, iel, iia1, iia2, nb
    real :: rlmax2, xi, yi, xi1, yi1, zi, zi1, a, b, c, det, rl0, rll, x0, y0
    real :: xperreal, yperreal, xe0, ye0, acc, rmaxhelp, xx0, yy0, rlen0
    real, dimension(:), pointer :: xb, yb
    integer, dimension(:,:), pointer :: ibb, ibp
    real, dimension(:), pointer :: rga, rgb, rgc
    real, dimension(:,:), pointer :: wp
    integer, dimension(:,:), pointer :: lnd, iae, lbn
    integer, dimension(:), pointer :: ibc, itc
    real, dimension(:), pointer :: x, y
    logical :: skip
    integer :: ic_start, ic_end, ic_skip

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lbn => AMA%lbn(1:AMA%mbelm,1:2)
    ibc => AMA%ibc(1:AMA%mbelm)
    itc => AMA%itc(1:AMA%mbelm)

    lnd => AMA%lnd(1:AMA%melem, 1:3)
    iae => AMA%iae(1:AMA%melem, 1:3)

    wp    => AMA%wp(   1:AMA%mpoin,1:ndim+1)

    rga => AMA%rga( 1:AMA%mpoin )
    rgb => AMA%rgb( 1:AMA%mpoin )
    rgc => AMA%rgc( 1:AMA%mpoin )

    ibp => AMA%ibp(1:AMA%mpoin, 1:2)
    ibb => AMA%ibb(1:AMA%mpoin, 1:3)

    xb => AMA%xb(1:AMA%ipoint)
    yb => AMA%yb(1:AMA%ipoint)

    icha = 0
    ice = 0

    itest = -3


    rlmax2 = 5.33

    do 5 i=1,AMA%melem
       AMA%nserr(i,1) = 0
5   enddo

    nelemold = AMA%nelem

    if(mod(icy, 2) == 0) then
       ic_start = 1;         ic_end = nelemold;   ic_skip =  1
    else
       ic_start = nelemold;  ic_end = 1;          ic_skip = -1
    endif

    !do 10 ipoc = 1,nelemold
    !do 10 ipoc = nelemold, 1, -1
    do 10 ipoc = ic_start, ic_end, ic_skip
       i = ipoc
       do 20 j=1,3
          j1 = mod(j,3) +1
          ii1 = lnd(i,j)
          ii2 = lnd(i,j1)
          xi = x(ii1)
          yi = y(ii1)
          xi1 = x(ii2)
          yi1 = y(ii2)
          zi = wp(ii1,1)
          zi1 = wp(ii2,1)
          a = (rga(ii1) + rga(ii2) )/2
          b = (rgb(ii1) + rgb(ii2) )/2
          c = (rgc(ii1) + rgc(ii2) )/2
          rmax(j) = ( a*(xi-xi1)*(xi-xi1) + c*(yi-yi1)*(yi-yi1)  &
               +2*b*(xi-xi1)*(yi-yi1))
          jmax(j) = j
20     enddo

       do 25 k=1,3
          do 26 l=1,2
             if(rmax(l) .lt. rmax(l+1) ) then
                rmaxhelp = rmax(l)
                rmax(l) = rmax(l+1)
                rmax(l+1) = rmaxhelp
                jmaxhelp = jmax(l)
                jmax(l) = jmax(l+1)
                jmax(l+1) = jmaxhelp
             endif
26        enddo
25     enddo

       !         acc = ACCUTE(x(lnd(i,1)),y(lnd(i,1)),x(lnd(i,2)),y(lnd(i,2)),
       !     *        x(lnd(i,3)),y(lnd(i,3)), 1.)
       !         if( (iae(i,1) .lt. 0 .or. iae(i,2) .lt. 0
       !     *        .or. iae(i,3) .lt. 0 ) .and. acc .gt. 1. ) then
       if(I == itest) then
          write(*,'(2e12.4,i5)') x(lnd(I,1)),y(lnd(I,1)),i
          write(*,'(3e12.4,i5)')   &
               x(lnd(I,2)),y(lnd(I,2)),rmax(1),jmax(1)
          write(*,'(3e12.4,i5)')   &
               x(lnd(I,3)),y(lnd(I,3)),rmax(2),jmax(2)
          write(*,'(3e12.4,i5)')   &
               x(lnd(I,1)),y(lnd(I,1)),rmax(3),jmax(3)
          print *
          print *,'#',jmax(1),jmax(2),jmax(3)
          print *,'#',rmax(1),rmax(2),rmax(3),AMA%nserr(i,1),acc
          print *,'#------------------------'
       endif


       do 28 l=1,3
          !     checking the dimension of arrays (even for periodical boundary)
          if(AMA%npoin .ge. AMA%mpoin-2 .or.AMA%nelem .ge. AMA%melem-4 .or.  &
               AMA%nbelm .ge. AMA%mbelm-2 ) then
             print *,'Dimension in insert full'
             print *,'nelem,AMA%melem=',AMA%nelem,AMA%melem
             print *,'npoin,mpoin=',AMA%npoin,AMA%mpoin
             print *,'nbelm,AMA%mbelm=',AMA%nbelm,AMA%mbelm
             return
          endif

          if(rmax(l) .ge. rlmax2) then
             !     ... NEW ACCUTE
             !            if(rmax(l) .ge. rlmax2 .or.
             !     *           (iae(i,jmax(l)) .lt. 0
             !     *           .and. acc*rmax(l)/2 .ge. rlmax2)  ) then

             j = jmax(l)
             !print*,'AAS###!@!',i,j,iae(i,j) ,AMA%nelem

             if(iae(i,j) .gt. 0) then
                !     for non boundary sides
                ii = iae(i,j)
                if(AMA%nserr(i,1) == 0 .and. AMA%nserr(ii,1) == 0) then
                   j1 = mod(j,3)+1
                   j2 = mod(j1,3)+1
                   do 30 jjj=1,3
                      if(iae(ii,jjj) == i) jj = jjj
30                 enddo
                   jj1 = mod(jj,3)+1
                   jj2 = mod(jj1,3)+1
                   k1 = lnd(i,j)
                   k2 = lnd(i,j1)
                   k3 = lnd(i,j2)
                   k4 = lnd(ii,jj2)
                   if(k2.ne.lnd(ii,jj) .or. k1.ne.lnd(ii,jj1))then
                      print *,'ERRROR in INSERT'
                      print *,i,k1,k2,k3,k4
                   endif

                   !     check if no validation of positivity
                   xx(1) = x(k1)
                   xx(2) = x(k4)
                   xx(3) = x(k2)
                   xx(4) = x(k3)
                   yy(1) = y(k1)
                   yy(2) = y(k4)
                   yy(3) = y(k2)
                   yy(4) = y(k3)
                   xx0 = (x(k1)+x(k2) )/2
                   yy0 = (y(k1)+y(k2) )/2

                   do 43 kl = 1,4
                      kl1 = mod(kl, 4) + 1

                      det = POS_TEST(xx0, yy0, xx(kl), yy(kl), xx(kl1), yy(kl1))
                      !print *,'..INS .',det
                      if( det .le. 1.) then
                         !     violation of positivity, go to next j
                         goto 28
                      endif
                      call POS1TEST(xx0,yy0,xx(kl),yy(kl),  &
                           xx(kl1),yy(kl1),itet)
                      if(itet == 1) then
                         !     violation of positivity, go to next i
                         goto 28
                      endif
43                 enddo

                   if(iae(i,j1) .gt. 0) then
                      ia1 = iae(i,j1)
                      do 40 kk =1,3
                         if(iae(ia1,kk) == i) ja1 = kk
40                    enddo
                   else
                      ia1 = -1
                   endif
                   if(iae(i,j2) .gt. 0) then
                      ia2 = iae(i,j2)
                      do 50 kk =1,3
                         if(iae(ia2,kk) == i) ja2 = kk
50                    enddo
                   else
                      ia2 = -1
                   endif
                   if(iae(ii,jj1) .gt. 0) then
                      iia1 = iae(ii,jj1)
                      do 60 kk =1,3
                         if(iae(iia1,kk) == ii) jja1 = kk
60                    enddo
                   else
                      iia1 = -1
                   endif
                   if(iae(ii,jj2) .gt. 0) then
                      iia2 = iae(ii,jj2)
                      do 70 kk =1,3
                         if(iae(iia2,kk) == ii) jja2 = kk
70                    enddo
                   else
                      iia2 = -2
                   endif

                   if(icha == -1) then
                      write(*,'(a2,8i5)') '@@',i,ii,k1,k2,k3,k4,  &
                           AMA%npoin,AMA%nelem
                   endif
                   x(AMA%npoin+1) = (x(k1)+x(k2))/2
                   y(AMA%npoin+1) = (y(k1)+y(k2))/2
                   ibb(AMA%npoin+1,1) = -1
                   ibb(AMA%npoin+1,2) = 0
                   ibb(AMA%npoin+1,3) = 0


                   lnd(i,j) = AMA%npoin+1
                   iae(i,j2) = AMA%nelem+1
                   iae(i,j) = AMA%nelem+2
                   lnd(ii,jj) = AMA%npoin+1
                   iae(ii,jj2) = AMA%nelem+2
                   iae(ii,jj) = AMA%nelem+1
                   lnd(AMA%nelem+1,1) = AMA%npoin+1
                   lnd(AMA%nelem+1,2) = k3
                   lnd(AMA%nelem+1,3) = k1
                   iae(AMA%nelem+1,1) = i
                   iae(AMA%nelem+1,2) = ia2
                   iae(AMA%nelem+1,3) = ii
                   if(ia2 .gt. 0) then
                      iae(ia2,ja2) = AMA%nelem+1
                   else
                      !     corection of adjacent triangle
                      do ib=1,AMA%nbelm
                         if(itc(ib) == i) itc(ib) = AMA%nelem+1
                      enddo
                   endif

                   lnd(AMA%nelem+2,1) = AMA%npoin+1
                   lnd(AMA%nelem+2,2) = k4
                   lnd(AMA%nelem+2,3) = k2
                   iae(AMA%nelem+2,1) = ii
                   iae(AMA%nelem+2,2) = iia2
                   iae(AMA%nelem+2,3) = i
                   if(iia2 .gt. 0) then
                      iae(iia2,jja2) = AMA%nelem+2
                   else
                      !     corection of adjacent triangle
                      do ib=1,AMA%nbelm
                         if(itc(ib) == ii) itc(ib) = AMA%nelem+2
                      enddo
                   endif

                   wp(AMA%npoin+1,1) = ( wp(k1,1) + wp(k2,1) )/2

                   rga(AMA%npoin+1) = (rga(k1) + rga(k2))/2
                   rgb(AMA%npoin+1) = (rgb(k1) + rgb(k2))/2
                   rgc(AMA%npoin+1) = (rgc(k1) + rgc(k2))/2
                   ibp(AMA%npoin+1,1) = 0
                   ibp(AMA%npoin+1,2) = 0

                   AMA%npoin = AMA%npoin + 1
                   AMA%nelem = AMA%nelem + 2
                   icha = icha + 1

                   if(AMA%ifig .ge. 0 .and. mod(icha,5) == 0 ) then
                      call PLOT1()
                      AMA%ifig = AMA%ifig + 1
                   endif

                   if(AMA%npoin .gt. AMA%mpoin .or. AMA%nelem .gt. AMA%melem) then
                      print *,'Error in dimension in insert'
                      print *,'nelem,AMA%melem=',AMA%nelem,AMA%melem
                      print *,'npoin,mpoin=',AMA%npoin,AMA%mpoin
                      stop
                   endif

                   !     ... in each triangle only one division
                   AMA%nserr(i,1) = -1
                   AMA%nserr(ii,1) = -1
                   AMA%nserr(AMA%nelem,1) = -1
                   AMA%nserr(AMA%nelem-1,1) = -1
                   goto 10
                endif
             else
                !     for boundary sides
                !if(I == itest)
                !print *,'$$',itest,AMA%nserr(i,1), lnd(I,1),lnd(I,2),lnd(i,3)
                if(AMA%nserr(i,1) == 0 )then
                   ipe = 0
999                j1 = mod(j,3)+1
                   j2 = mod(j1,3)+1
                   !print*,'###!@!',i,j
                   k1 = lnd(i,j)
                   k2 = lnd(i,j1)
                   k3 = lnd(i,j2)
                   ia1 = iae(i,j1)
                   ia2 = iae(i,j2)
                   il0 = -1
                   if(ipe == 1) goto 128
                   !     for periodic boundary, the second point

                   !     we seek the poin inthe field [xb(i),yb(i)] i=1,ipoint
                   rll = ((x(k1) - x(k2))*(x(k1) - x(k2)) +  &
                        (y(k1) - y(k2))*(y(k1) - y(k2)) )
                   x0 = (x(k1) + x(k2))/2
                   y0 = (y(k1) + y(k2))/2


                   !                     print*,'@@@@',x0,y0

                   if(ibb(k1,1) .gt. 0 .and. ibb(k2,1) .gt. 0   &
                        .and. ibb(k1,2) == ibb(k2,2) ) then
                      il1 = ibb(k1,1)
                      il2 = ibb(k2,1)
                      rl0 = 1E+25*(rll**0.5)

                      !     test if between il1 and il2 is a point
                      idif = AMA%ibpoin(ibb(k1,2))-AMA%ibpoin(ibb(k1,2)-1)-1
                      !                        print *,'@@@@@',il1,il2,idif
                      if(abs(il1-il2) == 1 ) then
                         !                        .or.
                         !     *                       abs(il1 - il2)  == idif ) then
                         print *,'We can not insert a new node,',  &
                              'there is few points on the profile'
                         print *,il1,il2,idif,AMA%ibpoin(ibb(k1,2)),  &
                              AMA%ibpoin(ibb(k1,2)-1)
                         print *,k1,x(k1),y(k1)
                         print *,k2,x(k2),y(k2)
                         goto 28
                      endif
                      if( il1 .gt. il2) then
                         !     0 node id between il1 and il2
                         il2new = il2 + AMA%ibpoin(ibb(k1,2))
                      else
                         il2new = il2
                      endif
                      do 213 ll1=il1,il2new
                         ll11 = ll1
                         if(ll11 .gt. AMA%ibpoin(ibb(k1,2)) )  &
                              ll11 = ll11 - AMA%ibpoin(ibb(k1,2))
                         rlen0 = (x0 -xb(ll11))*(x0 -xb(ll11)) +   &
                              (y0 -yb(ll11))*(y0 -yb(ll11))

                         if(rlen0 .lt. rl0) then
                            rl0 = rlen0
                            il0 = ll11
                         endif
213                   enddo
                      if(rl0 .gt. 0.3*rll) then
                         print *,'very divnyin INSERT.F',k1,k2
                         print *,x(k1),y(k1)
                         print *,x(k2),y(k2)
                         print *
                         print *,x0,y0
                         print *
                         print *,xb(il0),yb(il0)
                         print *,xb(il1),yb(il1)
                         print *,xb(il2),yb(il2)
                         stop
                      endif
                      x0 = xb(il0)
                      y0 = yb(il0)
                   endif

                   !if(i == itest) print *,'##',x(k1),y(k1),k1
                   !if(i == itest) print *,'##',x0,y0,il0
                   !if(i == itest) print *,'##',x(k3),y(k3),k3
                   !if(i == itest) print *,'##',x(k2),y(k2),k2


                   det = POS_TEST(x(k1), y(k1), x0, y0, x(k3), y(k3) )
                   !print *,'..INSa.',det

                   if( det .le. 1.) then
                      !     violation of positivity, go to next j
                      goto 28
                   endif

                   call POS1TEST(x0,y0,x(k1),y(k1),x(k3),y(k3),itet)
                   if(i == itest) print *,'? ?',itet
                   if(itet == 1) then
                      !     violation of positivity, go to next i
                      goto 28
                   endif


                   det = POS_TEST(x(k3), y(k3), x0, y0, x(k2), y(k2) )
                   !print *,'..INSa.',det

                   if( det .le. 1.) then
                      !     violation of positivity, go to next j
                      goto 28
                   endif

                   call POS1TEST(x0,y0,x(k2),y(k2),x(k3),y(k3),itet)
                   if(i == itest) print *,'? ?',itet,AMA%pos1
                   if(itet == 1) then
                      !     violation of positivity, go to next i
                      goto 28
                   endif

                   if(i == itest) write(*,'(a2,4e12.4)')  '##',x0,y0,det
                   if(i == itest) write(*,'(a2,3i5)')'@@', k1,ibp(k1,1),ibp(k1,2)
                   if(i == itest) write(*,'(a2,3i5)')'@@', k2,ibp(k2,1),ibp(k2,2)

                   skip = .false.
                   do nb =1,AMA%nbelm
                   !   print*,'#####',nb,lbn(nb,:), ibc(nb), itc(nb),AMA%iper(:,:)
                      if(itc(nb) == i) then
                         if(ibc(nb) == AMA%iper(1,1) .or. ibc(nb) == AMA%iper(1,2) .or. &
                              ibc(nb) == AMA%iper(2,1) .or. ibc(nb) == AMA%iper(2,2)) skip = .true.
                         goto 24
                      endif
                   enddo
                   print*,'The corresponding element does not found! FRRTY'
                   stop

24                 continue

                  !     periodic boundary
                   if(ibp(k1,1) .gt. 0 .and. ibp(k2,1) .gt. 0 .and. skip) then

                      if((ibp(k1,2) .ne. ibp(k2,2)) .and.  &
                           (ibp(k1,2) .ne. 3 .and. ibp(k2,2) .ne. 3))  &
                           goto 28

                      if(ibp(k1,2) == 3) then
                         ibper = ibp(k2,2)
                      else
                         ibper = ibp(k1,2)
                      endif

                      if(ibper == 1) then
                         xperreal = AMA%xper(1,1)
                         yperreal = AMA%xper(1,2)
                      else
                         xperreal = AMA%xper(2,1)
                         yperreal = AMA%xper(2,2)
                      endif


                      do 1000 iel =1,AMA%nelem
                         do 1001 jel=1,3
                            jel1 = mod(jel, 3) + 1
                            if(iae(iel,jel) .lt. 0) then
                               !write(*,'(A3,3i5,4(a2,2i5))')'#@@>', lnd(iel,1:3), &
                               !     '|',k1,k2,'|', ibp(k1,1), ibp(k2,1), &
                               !     '|', ibp(k1,2), ibp(k2,2), '|', jel,jel1
                               !write(*,*) x(lnd(iel,1)), y(lnd(iel,1))
                               !write(*,*) x(lnd(iel,2)), y(lnd(iel,2))
                               !write(*,*) x(lnd(iel,3)), y(lnd(iel,3))
                               !write(*,*) x(lnd(iel,1)), y(lnd(iel,1))
                               !write(*,*)
                               !write(*,*) x( ibp(k1,1)), y( ibp(k1,1))
                               !write(*,*) x( ibp(k2,1)), y( ibp(k2,1))

                               if(  (ibp(k1, 2) == 3 .and.  lnd(iel,jel) == ibp(k2,1) ).or. &
                                    (ibp(k2, 2) == 3 .and.  lnd(iel,jel1) == ibp(k1,1)) .or. &
                                    (ibp(k1, 2) /= 3 .and. ibp(k2,2) /= 3 .and. &
                                    lnd(iel,jel) == ibp(k2,1) .and. lnd(iel,jel1) == ibp(k1,1)) ) then

                                  !write(*,'(A3,16i5)')'#<>', lnd(iel,jel), ibp(k2,1), lnd(iel,jel1), ibp(k1,1)
                                  goto 1002
                               endif

                            endif

1001                     enddo
1000                  enddo
                      print*,'Corresponding element does not found in ama-anener90.f90'
                      !print *,x(k1),y(k1)
                      !print *,x(k2),y(k2),xperreal
                      !print *,x(k3),y(k3)
                      stop

1002                  continue

                      je1 = mod(jel,3)+1
                      je2 = mod(je1,3)+1
                      ke1 = lnd(iel,jel)
                      ke2 = lnd(iel,je1)
                      ke3 = lnd(iel,je2)

                      if(abs(x(k2) + xperreal - x(ke1) ).lt. 1E-05 .and.  &
                           abs(y(k2)+yperreal-y(ke1) ) .lt.  1E-05 ) then
                         imov = 1
                      elseif(abs(x(k2)-xperreal-x(ke1)) .lt. 1E-05 .and.  &
                           abs(y(k2)-yperreal-y(ke1)) .lt.   1E-05 ) then
                         imov = -1
                      else
                         print *,'BAD in insert in periodical points'
                         write(*,'(a6,6i5)') 'k1: ',k1, ibp(k1,1:2)
                         write(*,'(a6,6i5)') 'k2: ',k2, ibp(k2,1:2)
                         print *,i,k1,k2,k3
                         print *,iel,ke1,ke2,ke3
                         print *,x(k1),y(k1)
                         print *,x(k2),y(k2),xperreal
                         print *,x(k3),y(k3)
                         print*
                         print *,x(ke1),y(ke1),yperreal
                         print *,x(ke2),y(ke2)
                         print *,x(ke3),y(ke3)
                         print *,abs(x(k2) + xperreal - x(ke1) ),  &
                              abs(y(k2) + yperreal - y(ke1) ),  &
                              abs(x(k2) - xperreal - x(ke1) ),  &
                              abs(y(k2) - yperreal - y(ke1) )
                         stop
                      endif

                      xe0 = x0 + imov*xperreal
                      ye0 = y0 + imov*yperreal

                      det = POS_TEST(x(ke1), y(ke1), xe0, ye0, x(ke3), y(ke3) )
                      !print *,'..INSc.',det

                      if( det .le. 1.) then
                         !     violation of positivity, go to next j
                         goto 28
                      endif
                      call POS1TEST(xe0,ye0,x(ke1),y(ke1),  &
                           x(ke3),y(ke3),itet)
                      if(itet == 1) then
                         !     violation of positivity, go to next i
                         goto 28
                      endif

                      det = POS_TEST(x(ke3), y(ke3), xe0, ye0, x(ke2), y(ke2) )
                      !print *,'..INSd.',det

                      if( det .le. 1.) then
                         !     violation of positivity, go to next j
                         goto 28
                      endif
                      call POS1TEST(xe0,ye0,x(ke1),y(ke1),  &
                           x(ke2),y(ke2),itet)
                      if(itet == 1) then
                         !     violation of positivity, go to next i
                         goto 28
                      endif
                   endif

128                continue

                   if(ia1 .gt. 0) then
                      ja1 = 0
                      do 125 ll=1,3
                         if(iae(ia1,ll) == i) then
                            ja1 = ll
                         endif
125                   enddo
                      if(ja1 == 0) then
                         print *,'ERROR in INSERT_BOUNDARY-1'
                         stop
                      endif
                   endif

                   if(ipe == 0) then
                      x(AMA%npoin+1) = x0
                      y(AMA%npoin+1) = y0
                      ibb(AMA%npoin+1,1) = il0
                      ibb(AMA%npoin+1,3) = 0
                   elseif(ipe == 1) then
                      x(AMA%npoin+1) = xe0
                      y(AMA%npoin+1) = ye0
                      ibb(AMA%npoin+1,1) = il0
                      ibb(AMA%npoin+1,3) = 0
                   else
                      print *,'bad value of ipe =',ipe
                   endif

                   if(ibb(k1,2) == 0 .or. ibb(k2,2) == 0) then
                      ibb(AMA%npoin+1,2) = 0
                   elseif(ibb(k1,2) == ibb(k2,2)) then
                      ibb(AMA%npoin+1,2) = ibb(k1,2)
                   else
                      !                        print *,'error jkol1'
                   endif


                   wp(AMA%npoin+1,1) = (wp(k1,1) + wp(k2,1))/2
                   rga(AMA%npoin+1) = (rga(k1) + rga(k2))/2
                   rgb(AMA%npoin+1) = (rgb(k1) + rgb(k2))/2
                   rgc(AMA%npoin+1) = (rgc(k1) + rgc(k2))/2

                   lnd(i,j1) = AMA%npoin+1
                   iae(i,j1) = AMA%nelem+1

                   lnd(AMA%nelem+1,1) = AMA%npoin+1
                   lnd(AMA%nelem+1,2) = k2
                   lnd(AMA%nelem+1,3) = k3
                   iae(AMA%nelem+1,1) = -2
                   iae(AMA%nelem+1,2) = ia1
                   iae(AMA%nelem+1,3) = i

                   if(ia1 .gt. 0) iae(ia1,ja1) = AMA%nelem+1

                   !     the change in lbn, nbc,itc
                   ib1 = 0
                   do 265 ib=1,AMA%nbelm
                      if(lbn(ib,1) == k1 .and.lbn(ib,2) == k2)then
                         ib1 = ib
                         goto 266
                      endif
265                enddo
266                continue
                   if(ib1 == 0 ) then
                      print *,'ERROR in INSERT'
                      print *,'the boundary segment does not found'
                      stop
                   endif

                   do 276 ii=1,AMA%nbelm-ib1
                      ib = AMA%nbelm +2 -ii
                      lbn(ib,1) = lbn(ib-1,1)
                      lbn(ib,2) = lbn(ib-1,2)
                      ibc(ib) = ibc(ib-1)
                      itc(ib) = itc(ib-1)
276                enddo


                   lbn(ib1,2) = AMA%npoin + 1
                   lbn(ib1+1,1) = AMA%npoin + 1
                   lbn(ib1+1,2) = k2
                   ibc(ib1+1) = ibc(ib1)
                   itc(ib1+1) = AMA%nelem+1
                   if(ipe == 0) ibp(AMA%npoin+1,1) = 0
                   if(ipe == 0) ibp(AMA%npoin+1,2) = 0


                   AMA%npoin = AMA%npoin + 1
                   AMA%nelem = AMA%nelem + 1
                   AMA%nbelm = AMA%nbelm + 1

                   if(AMA%npoin .gt. AMA%mpoin .or.AMA%nelem .gt. AMA%melem .or.  &
                        AMA%nbelm .gt. AMA%mbelm ) then
                      print *,'ERROR in dimension in insert'
                      print *,'nelem,AMA%melem=',AMA%nelem,AMA%melem
                      print *,'npoin,mpoin=',AMA%npoin,AMA%mpoin
                      print *,'nbelm,AMA%mbelm=',AMA%nbelm,AMA%mbelm
                      stop
                   endif

                   AMA%nserr(i,1) = -1
                   AMA%nserr(AMA%nelem,1) = -1
                   if( ibp(k1,1) .gt. 0 .and. ibp(k2,1) .gt. 0 .and.   &
                        ipe == 0 .and.  skip) then
                      jbak = j
                      i = iel
                      j = jel
                      ibp(AMA%npoin,1) = AMA%npoin + 1
                      ibp(AMA%npoin+1,1) = AMA%npoin
                      ibp(AMA%npoin,2) = ibper
                      ibp(AMA%npoin+1,2) = ibper
                      ipe = 1
                      goto 999
                   endif
                   if(ipe == 1) then
                      ipe = 0
                      i = ipoc
                      j = jbak
                   endif
                   icha = icha + 1
                   goto 10
                endif
             endif
          endif
28     enddo
10  enddo
    !11 enddo

    return

  end subroutine INSERT



  subroutine INSERT_BOUND(ndim, icha, icy)
    implicit none
    integer, intent(in) :: ndim, icy
    integer, intent(inout) :: icha
    real :: rmax(3)
    integer jmax(3)
    integer:: nelemold ! local variable
    integer :: ice, itet, itest, ipoc, i, j, j1, j2, ii1, ii2, k, l, ll, ib, ib1
    integer :: ii, k1, k2, k3,  ia1, ia2, il0, ipe
    integer :: il1, il2, il2new, ll1, ll11, ibper, je1, ke1, ke2, ke3
    integer :: imov, ja1, jmaxhelp, jbak, je2, jel, jel1
    integer :: idif, iel, nb
    real :: rlmax2, xi, yi, xi1, yi1, zi, zi1, a, b, c, det, rl0, rll, x0, y0
    real :: xperreal, yperreal, xe0, ye0, acc, rmaxhelp, rlen0
    real, dimension(:), pointer :: xb, yb
    integer, dimension(:,:), pointer :: ibb, ibp
    real, dimension(:), pointer :: rga, rgb, rgc
    real, dimension(:,:), pointer :: wp
    integer, dimension(:,:), pointer :: lnd, iae, lbn
    integer, dimension(:), pointer :: ibc, itc
    real, dimension(:), pointer :: x, y
    logical :: skip
    integer :: ic_start, ic_end, ic_skip

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lbn => AMA%lbn(1:AMA%mbelm,1:2)
    ibc => AMA%ibc(1:AMA%mbelm)
    itc => AMA%itc(1:AMA%mbelm)


    lnd => AMA%lnd(1:AMA%melem, 1:3)
    iae => AMA%iae(1:AMA%melem, 1:3)

    wp    => AMA%wp(   1:AMA%mpoin,1:ndim+1)

    rga => AMA%rga( 1:AMA%mpoin )
    rgb => AMA%rgb( 1:AMA%mpoin )
    rgc => AMA%rgc( 1:AMA%mpoin )

    ibp => AMA%ibp(1:AMA%mpoin, 1:2)
    ibb => AMA%ibb(1:AMA%mpoin, 1:3)

    xb => AMA%xb(1:AMA%ipoint)
    yb => AMA%yb(1:AMA%ipoint)

    icha = 0
    ice = 0

    itest = -1183


    rlmax2 = 5.33

    do 5 i=1,AMA%melem
       AMA%nserr(i,1) = 0
5   enddo

    nelemold = AMA%nelem

    if(mod(icy, 2) == 0) then
       ic_start = 1;         ic_end = nelemold;   ic_skip =  1
    else
       ic_start = nelemold;  ic_end = 1;          ic_skip = -1
    endif

!    do 10 ipoc = 1,nelemold
    do 10 ipoc = ic_start, ic_end, ic_skip

       i = ipoc
       do 20 j=1,3
          j1 = mod(j,3) +1
          ii1 = lnd(i,j)
          ii2 = lnd(i,j1)
          xi = x(ii1)
          yi = y(ii1)
          xi1 = x(ii2)
          yi1 = y(ii2)
          zi = wp(ii1,1)
          zi1 = wp(ii2,1)
          a = (rga(ii1) + rga(ii2) )/2
          b = (rgb(ii1) + rgb(ii2) )/2
          c = (rgc(ii1) + rgc(ii2) )/2
          rmax(j) = ( a*(xi-xi1)*(xi-xi1) + c*(yi-yi1)*(yi-yi1)  &
               +2*b*(xi-xi1)*(yi-yi1))
          jmax(j) = j
20     enddo

       do 25 k=1,3
          do 26 l=1,2
             if(rmax(l) .lt. rmax(l+1) ) then
                rmaxhelp = rmax(l)
                rmax(l) = rmax(l+1)
                rmax(l+1) = rmaxhelp
                jmaxhelp = jmax(l)
                jmax(l) = jmax(l+1)
                jmax(l+1) = jmaxhelp
             endif
26        enddo
25     enddo

       !         if( (iae(i,1) .lt. 0 .or. iae(i,2) .lt. 0
       !     *        .or. iae(i,3) .lt. 0 ) .and. acc .gt. 1. ) then
       if(I == itest) then
          write(*,'(2e12.4,i5)') x(lnd(I,1)),y(lnd(I,1)),i
          write(*,'(3e12.4,i5)')   &
               x(lnd(I,2)),y(lnd(I,2)),rmax(1),jmax(1)
          write(*,'(3e12.4,i5)')   &
               x(lnd(I,3)),y(lnd(I,3)),rmax(2),jmax(2)
          write(*,'(3e12.4,i5)')   &
               x(lnd(I,1)),y(lnd(I,1)),rmax(3),jmax(3)
          print *
          print *,'#',jmax(1),jmax(2),jmax(3)
          print *,'#',rmax(1),rmax(2),rmax(3),AMA%nserr(i,1),acc
          print *,'#------------------------'
       endif


       do 28 l=1,3
          !     checking the dimension of arrays (even for periodical boundary)
          if(AMA%npoin .ge. AMA%mpoin-2 .or.AMA%nelem .ge. AMA%melem-4 .or.  &
               AMA%nbelm .ge. AMA%mbelm-2 ) then
             print *,'Dimension in insert_bound full'
             print *,'nelem,AMA%melem=',AMA%nelem,AMA%melem
             print *,'npoin,mpoin=',AMA%npoin,AMA%mpoin
             print *,'nbelm,AMA%mbelm=',AMA%nbelm,AMA%mbelm
             return
          endif

          acc = ACCUTE_I(x(lnd(i,1)),y(lnd(i,1)),  &
               x(lnd(i,2)),y(lnd(i,2)),x(lnd(i,3)),y(lnd(i,3)), l,1.)

          !            if(acc .gt. 1 .and. iae(i,l) .lt. 0) then
          !            if(iae(i,l) .lt. 0) then
          !               write(*,'(2e12.4,i5)') x(lnd(i,1)),y(lnd(i,1)),i
          !               write(*,'(2e12.4,i5)') x(lnd(i,1)),y(lnd(i,1)),iae(i,l)
          !               write(*,'(2e12.4)') x(lnd(i,1)),y(lnd(i,1))
          !               write(*,'(3e12.4)') x(lnd(i,1)),y(lnd(i,1)),acc
          !               write(*,'(x)')
          !            endif

          !            if(rmax(l) .ge. rlmax2) then
          !     ... NEW ACCUTE
          !     ... only non accute edges
          !write(91,'(5es12.4)') x(lnd(i,1)),y(lnd(i,1)),acc
          !write(91,'(5es12.4)') x(lnd(i,2)),y(lnd(i,2))
          !write(91,'(5es12.4)') x(lnd(i,3)),y(lnd(i,3))
          !write(91,'(5es12.4)') x(lnd(i,1)),y(lnd(i,1)),acc
          !write(91,'(x)')
          if(acc .gt. 10. ) then

            !               j = jmax(l)

             j = l
             if(iae(i,j) .gt. 0) then
                !     for non boundary sides
                !           no inserting

             else

                !write(92,'(3es12.4,i5)') x(lnd(i,1)),y(lnd(i,1)),acc, i
                !write(92,'(5es12.4)') x(lnd(i,2)),y(lnd(i,2))
                !write(92,'(5es12.4)') x(lnd(i,3)),y(lnd(i,3))
                !write(92,'(5es12.4)') x(lnd(i,1)),y(lnd(i,1)),acc
                !write(92,'(x)')

                !     for boundary sides
                if(I == itest) print *,'$$',itest,AMA%nserr(i,1),  &
                     lnd(I,1),lnd(I,2),lnd(i,3)
                if(AMA%nserr(i,1) == 0 )then
                   ipe = 0
999                j1 = mod(j,3)+1
                   j2 = mod(j1,3)+1
                   k1 = lnd(i,j)
                   k2 = lnd(i,j1)
                   k3 = lnd(i,j2)
                   ia1 = iae(i,j1)
                   ia2 = iae(i,j2)
                   il0 = -1
                   if(ipe == 1) goto 128
                   !     for periodic boundary, the second point

                   !     we seek the poin inthe field [xb(i),yb(i)] i=1,ipoint
                   rll = ((x(k1) - x(k2))*(x(k1) - x(k2)) +  &
                        (y(k1) - y(k2))*(y(k1) - y(k2)) )
                   x0 = (x(k1) + x(k2))/2
                   y0 = (y(k1) + y(k2))/2

                   if(ibb(k1,1) .gt. 0 .and. ibb(k2,1) .gt. 0   &
                        .and. ibb(k1,2) == ibb(k2,2) ) then
                      il1 = ibb(k1,1)
                      il2 = ibb(k2,1)
                      rl0 = 1E+25*(rll**0.5)

                      !     test if between il1 and il2 is a point
                      idif = AMA%ibpoin(ibb(k1,2))-AMA%ibpoin(ibb(k1,2)-1)-1
                      !                        print *,'@@@@@',il1,il2,idif
                      if(abs(il1-il2) == 1 ) then
                         !                        .or.
                         !     *                       abs(il1 - il2)  == idif ) then
                         print *,'We can not insert a new node,',  &
                              'there is few points on the profile'
                         print *,il1,il2,idif,AMA%ibpoin(ibb(k1,2)),  &
                              AMA%ibpoin(ibb(k1,2)-1)
                         print *,k1,x(k1),y(k1)
                         print *,k2,x(k2),y(k2)
                         goto 28
                      endif
                      if( il1 .gt. il2) then
                         !     0 node id between il1 and il2
                         il2new = il2 + AMA%ibpoin(ibb(k1,2))
                      else
                         il2new = il2
                      endif
                      do 213 ll1=il1,il2new
                         ll11 = ll1
                         if(ll11 .gt. AMA%ibpoin(ibb(k1,2)) )  &
                              ll11 = ll11 - AMA%ibpoin(ibb(k1,2))
                         rlen0 = (x0 -xb(ll11))*(x0 -xb(ll11)) +   &
                              (y0 -yb(ll11))*(y0 -yb(ll11))
                         if(rlen0 .lt. rl0) then
                            rl0 = rlen0
                            il0 = ll11
                         endif
213                   enddo
                      if(rl0 .gt. 0.3*rll) then
                         print *,'very divnyin INSERT.F',k1,k2
                         print *,x(k1),y(k1)
                         print *,x(k2),y(k2)
                         print *
                         print *,x0,y0
                         print *
                         print *,xb(il0),yb(il0)
                         print *,xb(il1),yb(il1)
                         print *,xb(il2),yb(il2)
                         stop
                      endif
                      x0 = xb(il0)
                      y0 = yb(il0)
                   endif

                   det = POS_TEST(x(k1), y(k1), x0, y0, x(k3), y(k3) )
                   !print *,'..IN B1.',det, x0, y0

                   if( det .le. 1.) then
                      !     violation of positivity, go to next j
                      goto 28
                   endif

                   call POS1TEST(x0,y0,x(k1),y(k1),x(k3),y(k3),itet)
                   if(i == itest) print *,'? ?',itet
                   if(itet == 1) then
                      !     violation of positivity, go to next i
                      goto 28
                   endif

                   det = POS_TEST(x(k3), y(k3), x0, y0, x(k2), y(k2) )
                   !print *,'..IN B2.',det, x0, y0

                   if( det .le. 1.) then
                      !     violation of positivity, go to next j
                      goto 28
                   endif
                   call POS1TEST(x0,y0,x(k2),y(k2),x(k3),y(k3),itet)
                   if(i == itest) print *,'? ?',itet,AMA%pos1
                   if(itet == 1) then
                      !     violation of positivity, go to next i
                      goto 28
                   endif

                   if(i == itest) then
                      write(*,'(a2,4e12.4)') '##',x0,y0,det
                      write(*,'(a2,8i5)')'@@',  k1,ibp(k1,1),ibp(k1,2),ibb(k1,:)
                      write(*,'(a2,8i5)')'@@',  k2,ibp(k2,1),ibp(k2,2),ibb(k2,:)
                      write(*,'(a2,8i5)')'!!',  i,j,iae(i,j),lbn(abs(iae(i,j)), 1:2)
                   endif

                   skip = .false.
                   do nb=1,AMA%nbelm
                      !print*,'#!!!#',nb,lbn(nb,:), ibc(nb), itc(nb),AMA%iper(:,:)
                      if(itc(nb) == i) then
                         if(ibc(nb) == AMA%iper(1,1) .or. ibc(nb) == AMA%iper(1,2) .or. &
                              ibc(nb) == AMA%iper(2,1) .or. ibc(nb) == AMA%iper(2,2)) skip = .true.
                         goto 24
                      endif
                   enddo
                   print*,'The corresponding element does not found! WERTY'
                   stop

24                 continue

                   !print*,'@@#  $$$   ', skip

                   !     periodic boundary
                   if(ibp(k1,1) .gt. 0 .and. ibp(k2,1) .gt. 0 .and. skip ) then
                      if((ibp(k1,2) .ne. ibp(k2,2)) .and.  &
                           (ibp(k1,2) .ne. 3 .and. ibp(k2,2) .ne. 3))  &
                           goto 28

                      if(ibp(k1,2) == 3) then
                         ibper = ibp(k2,2)
                      else
                         ibper = ibp(k1,2)
                      endif

                      if(ibper == 1) then
                         xperreal = AMA%xper(1,1)
                         yperreal = AMA%xper(1,2)
                      else
                         xperreal = AMA%xper(2,1)
                         yperreal = AMA%xper(2,2)
                      endif


                      do 1000 iel =1,AMA%nelem
                         do 1001 jel=1,3
                            jel1 = mod(jel, 3) + 1
                            if(iae(iel,jel) .lt. 0 ) then
                               !if( lnd(iel,jel) == ibp(k2,1) .and. lnd(iel,jel1) == ibp(k1,2)) then
                               if(  (ibp(k1, 2) == 3 .and.  lnd(iel,jel) == ibp(k2,1) ).or. &
                                    (ibp(k2, 2) == 3 .and.  lnd(iel,jel1) == ibp(k1,1)) .or. &
                                    (ibp(k1, 2) /= 3 .and. ibp(k2,2) /= 3 .and. &
                                    lnd(iel,jel) == ibp(k2,1) .and. lnd(iel,jel1) == ibp(k1,1)) ) then

                                  goto 1002
                               endif
                            endif

1001                     enddo
1000                  enddo
                      print*,'Corresponding element does not found in ama-anener90.f90 (2)'
                      stop

1002                  continue
                      je1 = mod(jel,3)+1
                      je2 = mod(je1,3)+1
                      ke1 = lnd(iel,jel)
                      ke2 = lnd(iel,je1)
                      ke3 = lnd(iel,je2)
                      if(abs(x(k2) + xperreal - x(ke1) ) < 1E-05 .and.  &
                           abs(y(k2)+yperreal-y(ke1) ) <  1E-05 ) then
                         imov = 1
                      elseif(abs(x(k2)-xperreal-x(ke1)) < 1E-05 .and.  &
                           abs(y(k2)-yperreal-y(ke1))  <   1E-05 ) then
                         imov = -1
                      else
                         print *,'BAD in insert in periodical points (2)'
                         !print *,i,k2,ke1
                         !print *,x(k2),y(k2),xperreal
                         !print *,x(ke1),y(ke1),yperreal
                         print *,abs(x(k2) + xperreal - x(ke1) ),  &
                              abs(y(k2) + yperreal - y(ke1) ),  &
                              abs(x(k2) - xperreal - x(ke1) ),  &
                              abs(y(k2) - yperreal - y(ke1) )
                         !stop

                      endif

                      xe0 = x0 + imov*xperreal
                      ye0 = y0 + imov*yperreal


                      det = POS_TEST(x(ke1), y(ke1), xe0, ye0, x(ke3), y(ke3) )
                      print *,'..IN B3.',det, xe0, ye0

                      if( det .le. 1.) then
                         !     violation of positivity, go to next j
                         goto 28
                      endif
                      call POS1TEST(xe0,ye0,x(ke1),y(ke1),  &
                           x(ke3),y(ke3),itet)
                      if(itet == 1) then
                         !     violation of positivity, go to next i
                         goto 28
                      endif

                      det = POS_TEST(x(ke3), y(ke3), xe0, ye0, x(ke2), y(ke2) )
                      print *,'..IN B4.',det, xe0, ye0

                      if( det .le. 1.) then
                         !     violation of positivity, go to next j
                         goto 28
                      endif
                      call POS1TEST(xe0,ye0,x(ke1),y(ke1),  &
                           x(ke2),y(ke2),itet)
                      if(itet == 1) then
                         !     violation of positivity, go to next i
                         goto 28
                      endif
                   endif

128                continue

                   if(ia1 .gt. 0) then
                      ja1 = 0
                      do 125 ll=1,3
                         if(iae(ia1,ll) == i) then
                            ja1 = ll
                         endif
125                   enddo
                      if(ja1 == 0) then
                         print *,'ERROR in INSERT_BOUNDARY-1'
                         stop
                      endif
                   endif

                   if(ipe == 0) then
                      x(AMA%npoin+1) = x0
                      y(AMA%npoin+1) = y0
                      ibb(AMA%npoin+1,1) = il0
                      ibb(AMA%npoin+1,3) = 0
                   elseif(ipe == 1) then
                      x(AMA%npoin+1) = xe0
                      y(AMA%npoin+1) = ye0
                      ibb(AMA%npoin+1,1) = il0
                      ibb(AMA%npoin+1,3) = 0
                   else
                      print *,'bad value of ipe =',ipe
                   endif

                   if(ibb(k1,2) == 0 .or. ibb(k2,2) == 0) then
                      ibb(AMA%npoin+1,2) = 0
                   elseif(ibb(k1,2) == ibb(k2,2)) then
                      ibb(AMA%npoin+1,2) = ibb(k1,2)
                   else
                      !                        print *,'error jkol1'
                   endif


                   wp(AMA%npoin+1,1) = (wp(k1,1) + wp(k2,1))/2
                   rga(AMA%npoin+1) = (rga(k1) + rga(k2))/2
                   rgb(AMA%npoin+1) = (rgb(k1) + rgb(k2))/2
                   rgc(AMA%npoin+1) = (rgc(k1) + rgc(k2))/2

                   lnd(i,j1) = AMA%npoin+1
                   iae(i,j1) = AMA%nelem+1

                   lnd(AMA%nelem+1,1) = AMA%npoin+1
                   lnd(AMA%nelem+1,2) = k2
                   lnd(AMA%nelem+1,3) = k3
                   iae(AMA%nelem+1,1) = -2
                   iae(AMA%nelem+1,2) = ia1
                   iae(AMA%nelem+1,3) = i

                   if(ia1 .gt. 0) iae(ia1,ja1) = AMA%nelem+1

                   !     the change in lbn, nbc,itc
                   ib1 = 0
                   do 265 ib=1,AMA%nbelm
                      if(lbn(ib,1) == k1 .and.lbn(ib,2) == k2)then
                         ib1 = ib
                         goto 266
                      endif
265                enddo
266                continue
                   if(ib1 == 0 ) then
                      print *,'ERROR in INSERT'
                      print *,'the boundary segment does not found'
                      stop
                   endif

                   do 276 ii=1,AMA%nbelm-ib1
                      ib = AMA%nbelm +2 -ii
                      lbn(ib,1) = lbn(ib-1,1)
                      lbn(ib,2) = lbn(ib-1,2)
                      ibc(ib) = ibc(ib-1)
                      itc(ib) = itc(ib-1)
276                enddo


                   lbn(ib1,2) = AMA%npoin + 1
                   lbn(ib1+1,1) = AMA%npoin + 1
                   lbn(ib1+1,2) = k2
                   ibc(ib1+1) = ibc(ib1)
                   itc(ib1+1) = AMA%nelem+1
                   if(ipe == 0) ibp(AMA%npoin+1,1) = 0
                   if(ipe == 0) ibp(AMA%npoin+1,2) = 0


                   AMA%npoin = AMA%npoin + 1
                   AMA%nelem = AMA%nelem + 1
                   AMA%nbelm = AMA%nbelm + 1

                   if(AMA%npoin .gt. AMA%mpoin .or.AMA%nelem .gt. AMA%melem .or.  &
                        AMA%nbelm .gt. AMA%mbelm ) then
                      print *,'ERROR in dimension in insert'
                      print *,'nelem,AMA%melem=',AMA%nelem,AMA%melem
                      print *,'npoin,mpoin=',AMA%npoin,AMA%mpoin
                      print *,'nbelm,AMA%mbelm=',AMA%nbelm,AMA%mbelm
                      stop
                   endif

                   AMA%nserr(i,1) = -1
                   AMA%nserr(AMA%nelem,1) = -1
                   if( ibp(k1,1) .gt. 0 .and. ibp(k2,1) .gt. 0 .and.   &
                       skip .and.  ipe == 0) then
                      jbak = j
                      i = iel
                      j = jel
                      ibp(AMA%npoin,1) = AMA%npoin + 1
                      ibp(AMA%npoin+1,1) = AMA%npoin
                      ibp(AMA%npoin,2) = ibper
                      ibp(AMA%npoin+1,2) = ibper
                      ipe = 1
                      goto 999
                   endif
                   if(ipe == 1) then
                      ipe = 0
                      i = ipoc
                      j = jbak
                   endif
                   icha = icha + 1
                   goto 10
                endif
             endif
          endif
28     enddo
10  enddo

    return
  end subroutine INSERT_BOUND

  subroutine METRIX(ndim,  surface)
    implicit none
    integer, intent(in) :: ndim
    real, intent(inout) :: surface
    integer :: i, j, k, ismoothing, is, l
    real :: xc, yc, rc, ro, u, v, p, xi, yi, ri, rkappa, rmeas, x3, y3
    real :: x1, y1, x2, y2
    real, dimension(:), allocatable :: tria, supp
    integer, dimension(:,:), pointer ::  ibp
    real, dimension(:,:), pointer :: wp, w
    integer, dimension(:,:), pointer :: lnd, iae, lbn
    integer, dimension(:), pointer :: ibc, itc
    real, dimension(:), pointer :: x, y

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lbn => AMA%lbn(1:AMA%mbelm,1:2)
    ibc => AMA%ibc(1:AMA%mbelm)
    itc => AMA%itc(1:AMA%mbelm)

    lnd => AMA%lnd(1:AMA%melem, 1:3)
    iae => AMA%iae(1:AMA%melem, 1:3)

    w     => AMA%w(   1:AMA%melem,1:ndim+1)
    wp    => AMA%wp(  1:AMA%mpoin,1:ndim+1)

    ibp => AMA%ibp(1:AMA%mpoin, 1:2)

    rkappa = 1.4
    allocate( tria(1:AMA%melem) , supp(1:AMA%mpoin) )

    surface = 0.
    do 2 i=1,AMA%nelem
       x1 = x(lnd(i,1))
       y1 = y(lnd(i,1))
       x2 = x(lnd(i,2))
       y2 = y(lnd(i,2))
       x3 = x(lnd(i,3))
       y3 = y(lnd(i,3))
       rmeas = (x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) )/2
       tria(i) = rmeas
       surface = surface + rmeas
       !write(21,*) (x1+x2+x3)/3, (y1+y2+y3)/3, w(i,1)

2   enddo

    if(AMA%ifv == 0) then
       do i = 1,AMA%npoin
          do k = 1,ndim
             wp(i,k+1) = w(i,k)
          enddo
       enddo
    else
       !     recomputation on the nodes and smoothing(for ismoothing > 1)
       if(AMA%ifv == 0) then
          ismoothing = 2
       else
          ismoothing = 1
       endif

       do 12 is=1,ismoothing
          do 10 i=1,AMA%npoin
             supp(i) = 0.
             do 15 j=1,ndim+1
                wp(i,j) = 0.
15           enddo
10        enddo
          do 20 i=1,AMA%nelem
             !xc = (x(lnd(i,1))+x(lnd(i,2))+x(lnd(i,3)))/3
             !yc = (y(lnd(i,1))+y(lnd(i,2))+y(lnd(i,3)))/3
             do 30 j=1,3
                !rl = ((xc-x(lnd(i,j)))**2 + (yc-y(lnd(i,j)))**2)**0.5
                do 40 k=1,ndim
                   wp(lnd(i,j),k+1) = wp(lnd(i,j),k+1) +w(i,k)*tria(i)
40              enddo
                supp(lnd(i,j)) = supp(lnd(i,j)) + tria(i)
30           enddo
20        enddo
          do 45 i=1,AMA%npoin
             do 60 k=2,ndim+1
                wp(i,k) = wp(i,k)/supp(i)
60           enddo
45        enddo



          do 48 i=1,AMA%nelem
             do 49 k=1,ndim
                w(i,k) = (wp(lnd(i,1),k+1) + wp(lnd(i,2),k+1) +   &
                     wp(lnd(i,3),k+1) )/3
49           enddo
48        enddo
12     enddo
    endif


    !     correction for the Navier - Stokes
    do i=1,AMA%npoin

       !write(22,*) x(i), y(i), wp(i,2)

       if(ibp(i,1) == -1 .and. ( AMA%ityp == 5)) then
          wp(i,3) = wp(i,3)*0.05
          wp(i,4) = wp(i,4)*0.05
       endif
       !     computation already done, we can forget this information
       !     for simplicity
       if(ibp(i,1) == -1) ibp(i,1) = 0
    enddo


    !     HERE IS POSSIBLE TO CHANGE THE USED QUANTITY FOR HESSIAN METRIXES
    !     for public
    !      do 50 i=1,AMA%npoin
    !         if(AMA%ityp == 0 ) then
    !     the uniform triangulation
    !            wp(i,1) = 1.0
    !         else
    !            wp(i,1) = wp(i,AMA%ityp+1)
    !         endif
    ! 50   enddo


    do 50 i=1,AMA%npoin
       !         write(99,*) x(i),y(i), wp(i,2),wp(i,3),wp(i,4), wp(i,5)


       if(AMA%ityp == 0 .or. AMA%ityp == 3) then
          !     the uniform triangulation
          wp(i,1) = 1.0
       elseif(AMA%ityp == -1) then
          !     the "exact mesh"
          xc = x(i)
          yc = y(i)
          !c2            wp(i,1) = (xc*xc + yc*yc)/2.
          !c3            wp(i,1) = (100*xc*xc + yc*yc)/2.
          !c4
          rc = (xc*xc+yc*yc)
          !c4
          wp(i,1) = 10*rc*exp(-10*(rc**0.5-1)**2)

       elseif(AMA%ityp == 1 .or. AMA%ityp == 4) then
          !     the testing quantity is the density
          wp(i,1) = wp(i,2)

       elseif(AMA%ityp == 2 ) then
          !     the testing quantity is the velocity
          wp(i,1) = ((wp(i,3)/wp(i,2))**2+(wp(i,4)/wp(i,2))**2 )**0.5
       elseif(AMA%ityp == 5 .or. AMA%ityp == 6) then
          !     the testing quantity is the Mach number
          ro = wp(i,2)
          if(ro .le. 0)then
             print *,'density zero on element ',i,'=',ro
             stop
          endif
          u  = wp(i,3)/ro
          v  = wp(i,4)/ro
          p  = (rkappa-1.)*(wp(i,5)-0.5*ro*(u*u+v*v))
          if( p .le. 0. ) then
             print *, ' Pressure <= 0.0 sur l''element ', i
             !               stop
             p = 0.001
          endif
          wp(i,1)=sqrt((u*u+v*v)*ro/(rkappa*p))

       elseif(AMA%ityp == -10)then
          if(x(i) < 0. ) then
             ri = (x(i)*x(i) + y(i)*y(i))**0.5
          elseif(x(i) >= 0. .and. x(i) .le. 1. ) then
             if(abs(y(i)) < 1.0) then
                ri = 1E+10
                do j=1,AMA%nbelm
                   do l=1,2
                      xi = x(lbn(j,l)) - x(i)
                      yi = y(lbn(j,l)) - y(i)
                      ri = min(ri, (xi*xi + yi*yi)**0.5 )
                   enddo
                enddo
             else
                ri = abs(y(i))
             endif
          else
             ri = (y(i)* y(i) + 0.025*x(i)*x(i))**0.5
          endif
          !ri = min( (x(i) -0.25)**2, (x(i) -0.5)**2, (x(i) -0.75)**2)  &
          !     + y(i)**2
          !wp(i,1) = 1.0 + 1E+8*exp(-2000*ri)
          !wp(i,1) = 100/(1+100 * ri**2)
          !wp(i,1) = 2.5*exp(-1 * max(0., ri-0.5 ) ) ! UA*
          wp(i,1) = 45*exp(-2.5 * max(0., ri-0.0 ) )   ! UB*
          if(x(i)*x(i) + y(i)*y(i) < 1E-2) wp(i,1) = wp(i,1)*5
          !write(99,*) x(i), y(i), ri
       else
          print *,'bad number of ityp, ityp = ',AMA%ityp
       endif

       !         write(21,'(i5,7e12.4)')
       !     *        i,x(i),y(i),wp(i,1),wp(i,2),wp(i,3),wp(i,4),wp(i,5)

50  enddo


    !     boundary layers
    if(AMA%ityp == 3) then
       !print*,'###',wp(lbn(100,1), 1),wp(lbn(100,2), 1)
       do i=1,AMA%nbelm
          if(ibc(i) == 3 .or. ibc(i) == 4.) then
             wp(lbn(i,1), 1) = wp(lbn(i,1), 1)*0.9
             !               wp(lbn(i,2), 1) = wp(lbn(i,2), 1)*0.8
             !               print*,'         ',i,ibc(i)
          endif
       enddo


       !print*,'###',wp(lbn(100,1), 1),wp(lbn(100,2), 1)

    endif


    !     periodic problems
    do 63 i=1,AMA%npoin
       do 64 k=1,ndim+1
          if(ibp(i,1) .gt. 0) then
             wp(i,k) = (wp(i,k) + wp(ibp(i,1),k))/2
             wp(ibp(i,1),k) = wp(i,k)
          endif
64     enddo
63  enddo

    !     now only a local array
    do 65 i=1,AMA%npoin
       supp(i) = 0.
65  enddo

    deallocate(tria, supp)
    return
  end subroutine METRIX



  subroutine ERROR1(ndim, surface)
    implicit none
    integer, intent(in) :: ndim
    real, intent(inout) :: surface
    real, dimension(:), allocatable :: dx, dy
    real, dimension(:), allocatable :: ra
    real, dimension(:), allocatable :: area, tria
    real :: dxi, dyi, rmeas, x1, x2, x3, y1, y2, y3, w1, w2, w3, areai, rl
    real :: a, b, c, disc, rlam1, rlam2, x11, x12, x21, x22, y12, y11, y21, y22
    real :: t11, t12, t21, t22, z11, z12, z21, z22, rdet, s, t
    real :: rkappa, err4, err5, err6, rmaxder
    real :: par, q, qbas, radius, rcharlen, Res,  Resminside, rmaxlambdas
    real :: rga1, rga2, rga3, rga4, rgai, rgb1, rgb2, rgb3, rgb4, rgbi
    real :: rgc1, rgc2, rgc3, rgc4, rgci, rgdi, rh1, rh2, rl1, rl2, rlampom
    real :: rlenrel, rmax, rmax1, rmax2, rmax3, rmaxderteor, rnorm, rminlenght
    real :: rn1, rn2, rradius, rwake, sstt, xc,  yc, xi, yi, xk, yk
    real :: xj, yj, rlen, xmax, xmin, ymax, ymin, xx, yy, ymax0, ymax1
    real :: eta, xlen, eps1, epsilon2, era1, era2, epsilon2p, dxfi, dyfi
    real :: derxxtrue, deryytrue, der, cor, beta, alpha, cc
    integer :: i, j, len, ii, ilkw, imax1, imax2, imax, imax3, imaxl, ipoc, imt
    integer :: ipr, is, itest, j1, je, je1, je2, k, k1, k2, l, n0, ie, icontrol
    integer :: i1, i2, ib, ismooth
    integer, dimension(:,:), pointer :: ibp
    real, dimension(:), pointer :: rga, rgb, rgc
    integer, dimension(:,:), pointer :: icyc
    real, dimension(:,:), pointer :: wp
    integer, dimension(:,:), pointer :: lnd, iae, lbn
    integer, dimension(:), pointer :: ibc, itc
    real, dimension(:), pointer :: x, y

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lbn => AMA%lbn(1:AMA%mbelm,1:2)
    ibc => AMA%ibc(1:AMA%mbelm)
    itc => AMA%itc(1:AMA%mbelm)


    lnd => AMA%lnd(1:AMA%melem, 1:3)
    iae => AMA%iae(1:AMA%melem, 1:3)

    wp    => AMA%wp(   1:AMA%mpoin,1:ndim+1)
    icyc => AMA%icyc(1:AMA%mpoin, 1:AMA%maxdeg)

    rga => AMA%rga( 1:AMA%mpoin )
    rgb => AMA%rgb( 1:AMA%mpoin )
    rgc => AMA%rgc( 1:AMA%mpoin )
    ibp => AMA%ibp(1:AMA%mpoin, 1:2)

    rkappa = 1.4
    allocate(dx(1:AMA%nelem),  dy(1:AMA%nelem) )
    allocate(ra(1:AMA%melem*3) )
    allocate(area(1:AMA%mpoin) )
    allocate(tria(1:AMA%melem) )

    !     dx the first derivative dw/dx
    !     dy                      dw/dy
    !     rga  the second der.    d2w/dx dx
    !     rgb                     d2w/dx dy
    !     rg!                     d2w/dy dy

    do 10 i=1,AMA%nelem
       x1 = x(lnd(i,1))
       y1 = y(lnd(i,1))
       w1 = wp(lnd(i,1),1)
       x2 = x(lnd(i,2))
       y2 = y(lnd(i,2))
       w2 = wp(lnd(i,2),1)
       x3 = x(lnd(i,3))
       y3 = y(lnd(i,3))
       w3 = wp(lnd(i,3),1)
       rmeas = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
       tria(i) = rmeas/2
       dx(i) = (y2-y1)*(w1+w2) + (y3-y2)*(w2+w3) + (y1-y3)*(w1+w3)
       dy(i) = (x1-x2)*(w1+w2) + (x2-x3)*(w2+w3) + (x3-x1)*(w1+w3)
       dx(i) = dx(i) /rmeas
       dy(i) = dy(i) /rmeas

       !         write(25,*) (x1+x2+x3)/3, (y1+y2+y3)/3, dx(i), dy(i)

       !         write(26,*) x(lnd(i,1)), y(lnd(i,1)), wp(lnd(i,1),1)
       !         write(26,*) x(lnd(i,2)), y(lnd(i,2)), wp(lnd(i,2),1)
       !         write(26,*) x(lnd(i,3)), y(lnd(i,3)), wp(lnd(i,3),1)
10  enddo


    rmaxder = 0.
    !     HERE INSERT deriv.f

    itest = -20

    err4 = 0.
    err5 = 0.
    err6 = 0.

    do 20 i=1,AMA%npoin
       areai = 0.
       rga(i) = 0.
       rgb(i) = 0.
       rgc(i) = 0.
       rgdi = 0.
       if(icyc(i,1) .gt. 0) then
          len = icyc(i,1)
       else
          len = -icyc(i,1)  - 1
       endif
       x3 = x(i)
       y3 = y(i)
       w3 = wp(i,1)
       do 30 j=1,len
          j1 = mod(j, abs(icyc(i,1))) + 1
          x1 = x(icyc(i,j+1))
          y1 = y(icyc(i,j+1))
          w1 = wp(icyc(i,j+1),1)
          x2 = x(icyc(i,j1+1))
          y2 = y(icyc(i,j1+1))
          w2 = wp(icyc(i,j1+1),1)
          rmeas = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
          areai = areai + rmeas/6
          dxi = (y2-y1)*(w1+w2) + (y3-y2)*(w2+w3) + (y1-y3)*(w1+w3)
          dyi = (x1-x2)*(w1+w2) + (x2-x3)*(w2+w3) + (x3-x1)*(w1+w3)
          dxfi = (y1-y2)/rmeas
          dyfi = (x2-x1)/rmeas

          if(i == itest) write(*,'(a2,i2,3e14.6,2e12.4)')  &
               '$$',j,w1,w2,w3,dxi/rmeas,dyi/rmeas
          if(i == itest) write(*,'(a2,4i5,4e12.4)')  &
               '  ',j,icyc(i,j+1),icyc(i,j1+1),i,x1,y1,x2,y2
          if(i == itest) write(*,'(a2,3e14.6)') '  ',x1,y1,w1
          if(i == itest) write(*,'(a2,3e14.6)') '  ',x2,y2,w2
          if(i == itest) write(*,'(a2,3e14.6)') '  ',x3,y3,w3

          rga(i) = rga(i) - dxi/2*dxfi
          rgb(i) = rgb(i) - dyi/2*dxfi
          rgc(i) = rgc(i) - dyi/2*dyfi
          rgdi = rgdi - dxi/2*dyfi

          !write(100 + AMA%adapt_level,*) &
          !            write(*,'(10es12.4)')  &
          !                 (x1+x2+x3)/3, (y1+y2+y3)/3,   dxi/rmeas, dyi/rmeas

          !            if(i == itest) print *,'**',j,dxi,dxfi,rga(i)

          if(icyc(i,1) .lt. 0 .and. ibp(i,1) == 0 ) then
             !     for the boundary points me must add the boundary sides
             if(j == 1) then
                rga(i) = rga(i) + dxi/rmeas*(y1-y3)/2
                rgb(i) = rgb(i) + dyi/rmeas*(y1-y3)/2
                rgc(i) = rgc(i) + dyi/rmeas*(x3-x1)/2
                rgdi = rgdi + dxi/rmeas*(x3-x1)/2
             endif
             if( j == len) then
                rga(i) = rga(i) + dxi/rmeas*(y3-y2)/2
                rgb(i) = rgb(i) + dyi/rmeas*(y3-y2)/2
                rgc(i) = rgc(i) + dyi/rmeas*(x2-x3)/2
                rgdi = rgdi + dxi/rmeas*(x2-x3)/2

             endif
          endif


          !            if(i == itest) print *,'**',j,dxi,dxfi,rga(i)

30     enddo
       area(i) = areai
       !         if(icyc(i,1) .gt. 0) then
       !         rradius = (x(i)*x(i) + y(i)*y(i))**0.5
       !         if(rradius .gt. 0) then
       !            derxxtrue = 3./16*rradius**(-1.75)
       !            deryytrue = derxxtrue
       !
       !            derxxtrue = -(1-y(i)**20)*90*x(I)**8
       !            deryytrue = -(1-x(i)**10)*380*y(i)**18
       !            derxytrue = 200*x(i)**9*y(i)**18
       !
       !            era1 = ( ( rga(i) ) /area(i)-derxxtrue)
       !            era2 = ( ( rgc(i) ) /area(i)-deryytrue)
       !            era3 = ( ( rgb(i) ) /area(i)-derxytrue)
       !            err4 = err4 + (era1*era1 + era2*era2 +
       !     *           2.*era3*era3)*area(i)
       !            err5 = err5 + (derxxtrue**2 + derxytrue**2 +
       !     *           2.*deryytrue**2)*area(i)
       !            err6 = err6 + ( rga(i)**2 + 2*rgb(i)**2 +
       !     *           rgc(i)**2)/area(i)
       !            write(98,'(i5,7e12.4)') i,area(i),
       !     *        rga(i)/area(i),derxxtrue,
       !     *        (era1*era1 + era2*era2)*area(i),
       !     *           (derxxtrue**2 + deryytrue**2)*area(i),
       !     *           err4,err5,err4/err5,err6
       !            endif
       !         endif

       if(AMA%ityp == -10) then
          rga(i) = wp(i,1)
          rgb(i) = 0.
          rgc(i) = rga(i)
       endif

       write(200 + AMA%adapt_level,*) x(i), y(i),   &
            rga(i)/area(i), rgb(i)/area(i), rgc(i)/area(i)

20  enddo

    !     for periodic boundary points
    if( AMA%xper(1,1) .gt. 0 .or. AMA%xper(1,2) .gt. 0) then
       do 18 i1=1,AMA%npoin
          i2 = ibp(i1,1)
          if( i2 .gt. 0 ) then
             if(rga(i1) .ne. rga(i2) .or. rgb(i1) .ne. rgb(i2) .or.   &
                  rgc(i1) .ne. rgc(i2) .or. area(i1) .ne. area(i2) )  &
                  then
                rga(i1) = rga(i1) + rga(i2)
                rgb(i1) = rgb(i1) + rgb(i2)
                rgc(i1) = rgc(i1) + rgc(i2)
                area(i1) = area(i1) + area(i2)
                rga(i2) = rga(i1)
                rgb(i2) = rgb(i1)
                rgc(i2) = rgc(i1)
                area(i2) = area(i1)
             endif
          endif
18     enddo
    endif

    !      print *,'Total error = ',err4**0.5,err5**0.5,
    !     *     err4**0.5/err5**0.5
    !      print *,'Total error = ',(err4/err6)**0.5

    do 15 i=1,AMA%npoin
       rga(i) = rga(i)/area(i)
       rgb(i) = rgb(i)/area(i)
       rgc(i) = rgc(i)/area(i)
       if(i == itest) print *,'..',area(i),rga(i)
       !         if(abs(rga(i)) + abs(rgc(i)) .gt. 5.)
       !         if(x(i) .gt. 0.95)
15  enddo


    !     ... comparing real second derivatives with their approximation
    do i=1,AMA%npoin
       rradius = (x(i)*x(i) + y(i)*y(i))**0.5
       if(rradius .gt. 0) then
          derxxtrue = 3./16*rradius**(-1.75)
       else
          derxxtrue = 1E+25
       endif
       deryytrue = derxxtrue

       !         derxxtrue = -(1-y(i)**20)*90*x(I)**8
       !         deryytrue = -(1-x(i)**10)*380*y(i)**18

       era1 = (abs(rga(i))-derxxtrue)
       era2 = (abs(rgc(i))-deryytrue)
       !         write(99,'(6e14.6)')
       !     *        x(i),y(i),rga(i),derxxtrue,rgc(i),deryytrue
       !     *        x(I),y(i),era1,era2,rga(i),rgc(i),derxxtrue
       !         rat = 250.*exp(-1*rradius)
       !         rga(i) = rat
       !         rgb(i) = rat
    enddo



    rmaxder = 0.
    !     improvement for boundary points
    goto 101

    do 100 i=1,AMA%npoin
       if(icyc(i,1) .lt. 0 .and. ibp(i,1) == 0) then
          len = abs(icyc(i,1))
          if( len .ge. 3) then
             ipoc = 0
             rga(i) = 0.
             rgb(i) = 0.
             rgc(i) = 0.
             rminlenght = 1E+38
             do 110 j=2,len-1
                if(icyc(icyc(i,j+1),1) .gt. 0) then
                   !     this point isn't boundary
                   xlen = ((x(i)-x(icyc(i,j+1)))**2 +  &
                        (y(i) -y(icyc(i,j+1)))**2 )**0.5
                   if(xlen .lt. rminlenght) then
                      rga(i) = rga(icyc(i,j+1))
                      rgb(i) = rgb(icyc(i,j+1))
                      rgc(i) = rgc(icyc(i,j+1))
                      rminlenght = xlen
                   endif
                   if( i == -170) then
                      print *,'**',j,icyc(i,j+1),rgc(icyc(i,j+1)),  &
                           xlen,rminlenght
                   endif
                endif
110          enddo
             !     we use the value from the nearest point
             if(ipoc .gt. 0 ) then
                rga(i) = rga(i)/ipoc
                rgb(i) = rgb(i)/ipoc
                rgc(i) = rgc(i)/ipoc
             else
                !     in this case, we let the second derivations = 0 !!!!

             endif
          elseif( len == 2) then
             i1 = icyc(i,2)
             i2 = icyc(i,3)
             do 120 ie =1,AMA%nelem
                do 130 je =1,3
                   je1 = mod(je,3)+1
                   if(lnd(ie,je) == i2 .and.   &
                        lnd(ie,je1) == i1) then
                      je2 = mod(je1,3) + 1
                      rga(i) = rga(lnd(ie,je2))
                      rgb(i) = rgb(lnd(ie,je2))
                      rgc(i) = rgc(lnd(ie,je2))
                      goto 140
                   endif
130             enddo
120          enddo
140          continue
          else
             print *,'ERROR len < 2 !!!!'
          endif
       endif
       rmaxder = max(rmaxder,abs(rga(i)),abs(rgb(i)),abs(rgc(i)) )
100 enddo
101 continue

    rnorm = 1.29903810567 *AMA%numel/surface
    rmaxderteor = AMA%epsilon1/3*rnorm
    der = rmaxderteor
    epsilon2p = AMA%epsilon1/AMA%p
    !      epsilon2p = p

    epsilon2 = max (1.,epsilon2p)

    write(AMA%ifig1,*)'   '
    write(AMA%ifig1,*)'Given data:'
    write(AMA%ifig1,*)'Used component for adaptation:',AMA%ityp,  &
         '  (0-uniform mesh)'
    write(AMA%ifig1,*)'Dimension of the solution:    ',ndim
    write(AMA%ifig1,*)'Positivity:                   ',AMA%pos
    write(AMA%ifig1,*)'Prescribed number of elements:',AMA%numel
    write(AMA%ifig1,*)'epsilon1:                     ',AMA%epsilon1
    write(AMA%ifig1,*)'p:                            ',AMA%p

    ipr = -1
    !     we compute the elements of matrix M from the second derivations

    rmax1 = 0.
    rmax2 = 0.
    rmax3 = 0.
    rmaxlambdas = 0.
    imaxl = 1

    !      itest = 50

    !      xmin1 = -0.15
    !      xmax1 = 0.05
    !      ymax1 = 0.01
    !      rmaa = 0.
    !      rmac = 0.
    !      rmac = 0.
    !      do i=1,AMA%npoin
    !         if(x(i) .gt. xmin1 .and. x(i).lt. xmax1 .and.
    !     *        y(i) .lt. ymax1) then
    !            rmaa = max(rmaa, rga(i) )
    !            rmab = max(rmaa, rgb(i) )
    !            rmac = max(rmaa, rgc(i) )
    !            write(49,*) x(i),y(i),rmaa,rmab,rmac
    !         endif
    !      enddo
    !      do i=1,AMA%npoin
    !         if(x(i) .gt. xmin1 .and. x(i).lt. xmax1 .and.
    !     *        y(i) .lt. ymax1) then
    !            rga(i) = 2*rmaa
    !            rgb(i) = 0.
    !            rgc(i) = 2*rmac
    !         endif
    !      enddo
    !      close(49)


    do 200 i=1,AMA%npoin

       if(rga(i) .gt. rmax1) then
          rmax1 = rga(i)
          imax1 = i
       endif
       if(rgb(i) .gt. rmax2) then
          rmax2 = rgb(i)
          imax2 = i
       endif
       if(rgc(i) .gt. rmax3) then
          rmax3 = rgc(i)
          imax3 = i
       endif

       !!write(*,'(a6,4es12.4)') '@@@@',rga(i), rgb(i), rgc(i)

       if(abs(rgb(i)) .lt. 1E-05) then
          !     the diagonal matrix

          eps1 = AMA%epsilon1/(epsilon2+max(rga(i),rgc(i)))

          rgb(i) = 0.
          rga(i) = (1. + eps1*abs(rga(i)))*rnorm
          rgc(i) = (1. + eps1*abs(rgc(i)))*rnorm
          if(rga(i) + rgc(i) .gt. rmaxlambdas)   &
               rmaxlambdas = rga(i)+rgc(i)
       else
          !     the eigenvalues rlam1 , rlam2
          if(i == itest) then
             print *,x(i),y(i)
             print *,'@#',rga(i),rgb(i),rgc(i)
          endif
          a = rga(i)
          b = rgb(i)
          c = rgc(i)

          disc = ((a - c)*(a - c) + 4*b*b)**0.5
          rlam1 = (a + c + disc)/2
          rlam2 = (a + c - disc)/2
          if(abs(rlam1) .lt. abs(rlam2) ) then
             rlampom = rlam1
             rlam1 = rlam2
             rlam2 = rlampom
          endif
          if(abs(rlam1) + abs(rlam2) .gt. rmaxlambdas) then
             rmaxlambdas = abs(rlam1) + abs(rlam2)
             imaxl = i
          endif
          x11 = b
          x21 = -(a-rlam1)
          rl1 = (x11*x11+x21*x21)**0.5
          x11 = x11/rl1
          x21 = x21/rl1

          x12 = b
          x22 = -(a-rlam2)
          rl2 = (x12*x12+x22*x22)**0.5
          x12 = x12/rl2
          x22 = x22/rl2

          rdet = x11*x22 - x21*x12

          !     inverse matrix
          y11 = x22/rdet
          y12 = -x12/rdet
          y21 = -x21/rdet
          y22 = x11/rdet

          z11 = abs(rlam1)*y11
          z12 = abs(rlam1)*y12
          z21 = abs(rlam2)*y21
          z22 = abs(rlam2)*y22

          t11 = x11*z11 + x12*z21
          t12 = x11*z12 + x12*z22
          t21 = x21*z11 + x22*z21
          t22 = x21*z12 + x22*z22

          eps1 = AMA%epsilon1/(epsilon2+max(t11,t22))

          if(i == itest) print *,' #',eps1,epsilon2,rnorm
          if(i == itest) print *,' #',eps1*t11+1.,rnorm

          rga(i) = (1. + eps1*t11)*rnorm
          rgb(i) = (eps1*(t12+t21)/2)*rnorm
          rgc(i) = (1. + eps1*t22)*rnorm

          if(i == itest) print *,'@#',rga(i),rgb(i),rgc(i)

          if(rga(i)*rgc(i) .le. rgb(i)*rgb(i) ) then
             print *,'ERROR in metric, nonpositive matrix'
             write(*,'(i5,3e14.6)') i,rga(i),rgb(i),rgc(i)
             print *,a,b,c
             print *,'----------------------'
             write(*,'(6e12.4)')  &
                  x11,x12,rlam1,0,y11,y12
             write(*,'(6e12.4)')  &
                  x21,x22,0,rlam2,y21,y22
             print *,'----------------------'
             write(*,'(4e12.4)')  &
                  z11,z12,t11,t12
             write(*,'(4e12.4)')  &
                  z21,z22,t21,t22
             print *,'----------------------'
             stop
          endif
          !write(23,*) x(i),y(i),rga(i),rgb(i),rgc(i)
       endif
200 enddo
    !      print *,'###############',rmaxlambdas,imaxl,
    !     *     x(imaxl),y(imaxl)

    !     CORRECTION of positivity a priori
    par = 1.0
    do i=1,AMA%npoin
       a = rga(i)
       b = rgb(i)
       c = rgc(i)
       disc = ((a - c)*(a - c) + 4*b*b)**0.5
       rlam1 = (a + c + disc)/2
       radius = 4*AMA%pos*AMA%pos*rlam1/par/par
       !         rh0 = 1./(rlam1)**0.5
       !         rh1 = rh0/3/AMA%pos*par
       !         radius = 1/rh1/rh1
       !     rgai,rgbi,rgci - cyrcle whose radius is the pos times the smallest
       !     side
       rgai = radius
       rgbi = 0.
       rgci = radius
       if(i == -20) then
          rgai = 4.
          rgbi = 2.
          rgci = 12.
          rga(i) = 3.
          rgb(i) = 1.5
          rgc(i) = 9.
          print *,x(i),y(i),rlam1,radius
          print *,rga(i),rgb(i),rgc(i)
          print *,rgai,rgbi,rgci
          call ELIPS(14,rga(i),rgb(i),rgc(i),x(i),y(i) )
          call ELIPS(14,rgai,rgbi,rgci,x(i),y(i) )
       endif
       if(i == itest) print *,'$$',rga(i),rgb(i),rgc(i)

       call INTERSECTION_METRIC(rga(i),rgb(i),rgc(i),  &
            rgai,rgbi,rgci)
       if(i == itest) print *,'$2',rga(i),rgb(i),rgc(i)
       !        if(i == 20) then
       !            call ELIPS(14,rga(i),rgb(i),rgc(i),x(i),y(i) )
       !            print *,rga(i),rgb(i),rgc(i)
       !            stop
       !         endif
    enddo
    !     end of CORRECTION of positivity a priori

    !     smoothing
    print*,'################### smoothing', 0
    ismooth = 0
    !ismooth = 6
    if(ismooth > 0) then
       do 500 ic = 1,ismooth
          do 400 i=1,AMA%nelem
             ra(i) = (rga(lnd(i,1)) + rga(lnd(i,2)) + rga(lnd(i,3)))/3
             ra(i+AMA%melem) =   &
                  (rgb(lnd(i,1)) + rgb(lnd(i,2)) + rgb(lnd(i,3)))/3
             ra(i+2*AMA%melem) =   &
                  (rgc(lnd(i,1)) + rgc(lnd(i,2)) + rgc(lnd(i,3)))/3
400       enddo

          do 410 i=1,AMA%npoin
             rga(i) = 0.
             rgb(i) = 0.
             rgc(i) = 0.
             area(i) = 0.
410       enddo

          do 420 i=1,AMA%nelem
             xc = (x(lnd(i,1))+x(lnd(i,2))+x(lnd(i,3)))/3
             yc = (y(lnd(i,1))+y(lnd(i,2))+y(lnd(i,3)))/3
             do 430 j=1,3
                k = lnd(i,j)
                !               rl = ((xc-x(k))**2 + (yc-y(k))**2)**0.5
                rga(k) = rga(k) + ra(i)*tria(i)
                rgb(k) = rgb(k) + ra(i+AMA%melem)*tria(i)
                rgc(k) = rgc(k) + ra(i+2*AMA%melem)*tria(i)
                area(k) = area(k) + tria(i)
430          enddo
420       enddo

          do 440 i=1,AMA%npoin
             rga(i) = rga(i)/area(i)
             rgb(i) = rgb(i)/area(i)
             rgc(i) = rgc(i)/area(i)
             if(i == itest) print *,'$3',rga(i),rgb(i),rgc(i)
440       enddo
500    enddo
    endif

    ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


    if((AMA%ityp == 3 .or. AMA%ityp == 4 .or. AMA%ityp == 5) .and.  &
         AMA%Re .gt. 0.) then


       !     part of courved boundary, where is fixed walls
       print*,'$$$$$$$$$$$$$  ATTENTION HERE'

       ilkw = 1

       !     generation the triangulation for the boundary layer for NS
       qbas = 1.25
       !     charakterictic length
       xmax = -100000.
       ymax0 = -100000.
       xmin = 100000.
       ymin = 100000.
       !         do i=1,ipoint
       !            if( i .gt. AMA%ibpoin(ilkw-1) .and.
       !     *           i .le. AMA%ibpoin(ilkw)) then
       !               xmax = max(xmax, xb(i) )
       !               ymax0 = max(ymax0, yb(i) )
       !               xmin = min(xmin, xb(i) )
       !               ymin = min(ymin, yb(i) )
       !            endif
       !         enddo

       !         print*, xmax, ymax0

       do i=1,AMA%nbelm
          do j=1,AMA%iwa
             if(ibc(i) == AMA%iwall(j) )then
                xmax = max(xmax, x(lbn(i,1)), x(lbn(i,2)) )
                ymax0 = max(ymax0, y(lbn(i,1)), y(lbn(i,2))  )
                xmin = min(xmin, x(lbn(i,1)), x(lbn(i,2))  )
                ymin = min(ymin, y(lbn(i,1)), y(lbn(i,2))  )
             endif
          enddo
       enddo

       !         print*, xmin, xmax, ymin, ymax0

       rcharlen = ((xmax -xmin)**2 + (ymax0 - ymin)**2)**0.5
       !     double precision
       x1 = AMA%xte(1,1)
       y1 = AMA%xte(1,2)
       x2 = AMA%xte(2,1)
       y2 = AMA%xte(2,2)

       !     blasius
       !     DMR
       !         rcharlen = 1.

       Resminside = rmaxlambdas/3
       Res = min(AMA%Re, 1000000*Resminside)
       !         Res = min(Re, 10000*Resminside)
       if(AMA%Re .gt. Res) print*,'Re  <  10000*Resminside !!!'
       Res = Res/rcharlen

       print*,'Reynolds =',Res,'  (Resminside, Re):',Resminside,AMA%Re

       !         ymax = ((3/rnorm)**0.5 - 1./sqrt(Res))/(qbas-1.)*rcharlen
       ymax1 = ((3/rnorm)**0.5 - 1./sqrt(Res))/(qbas-1.)*rcharlen
       ymax = 20./sqrt(Res)
       !         ymax = 60./sqrt(Res)
       n0 = 3

       imt = 18
       open(imt,file='mt0',status='unknown')

       !         print*,'$$$',rcharlen, xmin,xmax, ymax0, ymin,AMA%npoin

       do i=1,AMA%npoin
          xi = x(i)
          yi = y(i)
          rmax = 100000.

          do ib=1,AMA%nbelm
             do j=1, AMA%iwa
                if(ibc(ib) == AMA%iwall(j) )then

                   xj = (x(lbn(ib,1))+ x(lbn(ib,1)))/2
                   yj = (y(lbn(ib,1))+ y(lbn(ib,1)))/2

                   rlen = (xi-xj)**2 + (yi-yj)**2
                   if(rlen .lt. rmax) then
                      rmax = rlen
                      imax = ib
                   endif

                endif
             enddo
          enddo


          !           do k=AMA%ibpoin(ilkw-1)+1,AMA%ibpoin(ilkw)
          !               rlen = (xi-xb(k))**2 + (yi-yb(k))**2
          !               if(rlen .lt. rmax) then
          !                  rmax = rlen
          !                  imax = k
          !               endif
          !            enddo
          !     distance from profile
          rmax = rmax**0.5

          !     channel
          !            rmax = min ( (1- y(i)), y(i) )
          !            imax = -10
          !            n0 = 1
          !     blasius
          !            if(xi .ge. 0 ) then
          !               rmax = abs(yi)
          !            else
          !               rmax = (xi*xi + yi*yi/10.)**0.5
          !            endif

          x3 = xi
          y3 = yi

          sstt = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)
          s = ( (x1-x2)*(x1-x3) + (y1-y2)*(y1-y3) )/sstt
          t = ( x1*(y3-y2) + x2*(y1-y3) + x3*(y2-y1) )/sstt
          !     wake and boundary layer are investigated separetely,
          !     as some points need as to BL as to Wake


          !     BOUNDARY LAYERS
          !            print*,'$$$',i,s,rcharlen, rmax, ymax

          if(s .le. 0.01*rcharlen) then

             !               if(rmax .gt. ymax .and. rmax .lt. 5*ymax)
             !                   write(*,*) x(i),y(i),rmax,ymax
             rwake = 2*rmax + 1.
             if( rmax .lt. ymax) then
                !               if( rmax .lt. 5*ymax) then
                cor = 1./n0

                !                  if(imax .ge. 0) then
                !                     if(imax .gt. 1) then
                !                        ii1 = imax - 1
                !                     else
                !                        ii1 = ipoint - 1
                !                     endif
                !                     rn1 = -(yb(ii1) - yb(imax))
                !                     rn2 = xb(imax) - xb(ii1)
                !                  endif

                rn1 = y(lbn(imax,2)) - y(lbn(imax,1))
                rn2 = x(lbn(imax,1)) - x(lbn(imax,2))

                !                  write(*,'(i5,6e12.4)') i, xi,yi,
                !     *                 x(lbn(imax,1)), y(lbn(imax,1)),rn1,rn2

                !                  q = 5.*qbas
                q = (0.2+qbas)/1.2 -ymax/20.
                q = max(q, 0.)

                rlen = (rn1*rn1 + rn2*rn2)**0.5
                rn1 = rn1/rlen
                rn2 = rn2/rlen

                rh2 = cor/sqrt(Res) + rmax*(q - 1.)/rcharlen

                !               print*,'###',rh2,cor, rmax, q, rn1, rn2
                !                  rh2 = cor/sqrt(Res)+ (rcharlen/10. - cor/sqrt(Res))*
                !     *                 rmax/ymax
                rh2 = rh2*4.
                rh1 = rh2*(1+(1-24.*AMA%pos*AMA%pos)**0.5)/4/AMA%pos
                !     DMR
                !                  rh1 = rh2

                !                  write(*,*) x(i),y(i),q,rh1,rh2
                rgai = rn2*rn2/rh1/rh1 + rn1*rn1/rh2/rh2
                rgci = rn2*rn2/rh2/rh2 + rn1*rn1/rh1/rh1
                rgbi = rn1*rn2*(1./rh1/rh1 - 1./rh2/rh2)
                call INTERSECTION_METRIC(rga(i),rgb(i),rgc(i),  &
                     rgai,rgbi,rgci)

                !     for plotting of the elipses
                !                  if(xi .gt. -0.01 .and. xi .lt. 0.2 .and.
                !     *                 abs(y(i)) .lt. 0.02) then
                !                     call ELIPS(imt,rga(i),rgb(i),rgc(i),x(i),y(i) )
                call ELIPS(imt,rgai,rgbi,rgci,x(i),y(i) )
                !                     print*,'##',i,xi,yi,rmax, ymax
                !                  endif

             endif
          endif

          !     WAKES
          if(s .ge. 0.0) then
             xk = x1 + s*(x2-x1)
             yk = y1 + s*(y2-y1)
             rwake = ( ( xk - x3)**2 + (yk - y3)**2)**0.5
             if( rwake .lt. ymax) then
                !                  cor = (2.+ s*4.)*0.15/n0
                rlenrel = s*sqrt(sstt)/rcharlen
                !                  cor = (1.+rlenrel)/n0*2.
                cor = (1.+rlenrel)/n0*6.
                !                  write(46,*) xi,yi,cor
                rn1 = -(AMA%xte(1,2) - AMA%xte(2,2))
                rn2 = AMA%xte(2,1) - AMA%xte(1,1)
                rmax = rwake-0.02*rcharlen
                rmax = max (rmax,0.)
                !                  q = (0.2+qbas)/1.2 - ymax/20.
                !                  q = max(q, 0.)
                !                  q = qbas
                rlen = (rn1*rn1 + rn2*rn2)**0.5
                rn1 = rn1/rlen
                rn2 = rn2/rlen

                !                  rh2 = cor/sqrt(Res) + rmax*(q - 1.)/rcharlen
                rh2 = cor/sqrt(Res)+ (rcharlen/10. - cor/sqrt(Res))*  &
                     rmax/ymax
                rh1 = rh2*(1+(1-24.*AMA%pos*AMA%pos)**0.5)/4/AMA%pos
                !               write(33,*) x(i),y(i),q,rh1,rh2
                rgai = rn2*rn2/rh1/rh1 + rn1*rn1/rh2/rh2
                rgci = rn2*rn2/rh2/rh2 + rn1*rn1/rh1/rh1
                rgbi = rn1*rn2*(1./rh1/rh1 - 1./rh2/rh2)

                !               call ELIPS(64,rgai,rgbi,rgci,x(i),y(i) )
                call INTERSECTION_METRIC(rga(i),rgb(i),rgc(i),  &
                     rgai,rgbi,rgci)
             endif
          endif
       enddo

       close(imt)

    endif
    !      stop


    icontrol = -1
    if(icontrol == 1) then
       !     MESH GRADIENT CONTROL H-variation
       alpha = 2.0
       do is =1,4
          do k=1,3*AMA%nelem
             ra(k) = 1.
          enddo
          do i=1,AMA%nelem
             do j=1,3
                if(ra(3*(i-1)+j) .gt. 0) then
                   j1 = mod(j,3) + 1
                   k1 = lnd(i,j)
                   k2 = lnd(i,j1)
                   x1 = x(k1)
                   y1 = y(k1)
                   x2 = x(k2)
                   y2 = y(k2)
                   rga1 = rga(k1)
                   rgb1 = rgb(k1)
                   rgc1 = rgc(k1)
                   rga2 = rga(k2)
                   rgb2 = rgb(k2)
                   rgc2 = rgc(k2)
                   rl1 = (rga1*(x1-x2)**2 + 2*rgb1*(x1-x2)*(y1-y2) +  &
                        rgc1*(y1-y2)**2)**0.5
                   rl2 = (rga2*(x1-x2)**2 + 2*rgb2*(x1-x2)*(y1-y2) +  &
                        rgc2*(y1-y2)**2)**0.5
                   rga3 = rga2*(1+alpha*rl1)**(-2)
                   rgb3 = rgb2*(1+alpha*rl1)**(-2)
                   rgc3 = rgc2*(1+alpha*rl1)**(-2)
                   rga4 = rga1*(1+alpha*rl2)**(-2)
                   rgb4 = rgb1*(1+alpha*rl2)**(-2)
                   rgc4 = rgc1*(1+alpha*rl2)**(-2)
                   call INTERSECTION_METRIC(rga(k1),rgb(k1),rgc(k1),  &
                        rga3,rgb3,rgc3)
                   call INTERSECTION_METRIC(rga(k2),rgb(k2),rgc(k2),  &
                        rga4,rgb4,rgc4)

                   ra(3*(i-1)+j) = -1.
                   if(iae(i,j) .gt. 0) then
                      ii = iae(i,j)
                      do l=1,3
                         if(iae(ii,l) == i) ra(3*(ii-1)+l) = -1.
                      enddo
                   endif
                endif
             enddo
          enddo
       enddo
       !     end of MESH GRADIENT CONTROL

    elseif(icontrol == 2) then
       !     MESH GRADIENT CONTROL H-shock
       beta = 1.5
       do is =1,2
          do k=1,3*AMA%nelem
             ra(k) = 1.
          enddo
          do i=1,AMA%nelem
             do j=1,3
                if(ra(3*(i-1)+j) .gt. 0) then
                   j1 = mod(j,3) + 1
                   k1 = lnd(i,j)
                   k2 = lnd(i,j1)
                   x1 = x(k1)
                   y1 = y(k1)
                   x2 = x(k2)
                   y2 = y(k2)
                   rga1 = rga(k1)
                   rgb1 = rgb(k1)
                   rgc1 = rgc(k1)
                   rga2 = rga(k2)
                   rgb2 = rgb(k2)
                   rgc2 = rgc(k2)
                   rnorm = ((x2-x1)**2 + (y2-y1)**2)**0.5
                   xx = (x2 -x1)/rnorm
                   yy = (y2 -y1)/rnorm
                   rl1 = (rga1*xx**2 + 2*rgb1*xx*yy +  &
                        rgc1*yy**2)**0.5
                   rl2 = (rga2*xx**2 + 2*rgb2*xx*yy +  &
                        rgc2*yy**2)**0.5
                   rl = ((rga1+rga2)/2*(x1-x2)**2 +   &
                        (rgb1+rgb2)*(x1-x2)*(y1-y2) +  &
                        (rgc1+rgc2)/2*(y1-y2)**2)**0.5
                   if(rl2 .gt. rl1) then
                      cc = (rl2/rl1)**(1./rl)
                      if(cc .gt. beta) then
                         eta = (beta/cc)**rl
                         print *,'eta =',eta
                         rga(k2) = rga(k2)/eta/eta
                         rgb(k2) = rgb(k2)/eta/eta
                         rgc(k2) = rgc(k2)/eta/eta
                      endif
                   endif
                   ra(3*(i-1)+j) = -1.
                   if(iae(i,j) .gt. 0) then
                      ii = iae(i,j)
                      do l=1,3
                         if(iae(ii,l) == i) ra(3*(ii-1)+l) = -1.
                      enddo
                   endif
                endif
             enddo
          enddo
       enddo
       !     end of MESH GRADIENT CONTROL
    endif


    !     smoothing
    print*,'################### smoothing', 2
    do 1500 ic = 1,2
       !      do 1500 ic = 1,2
       do 1400 i=1,AMA%nelem
          ra(i) = (rga(lnd(i,1)) + rga(lnd(i,2)) + rga(lnd(i,3)))/3
          ra(i+AMA%melem) =   &
               (rgb(lnd(i,1)) + rgb(lnd(i,2)) + rgb(lnd(i,3)))/3
          ra(i+2*AMA%melem) =   &
               (rgc(lnd(i,1)) + rgc(lnd(i,2)) + rgc(lnd(i,3)))/3
1400   enddo

       do 1410 i=1,AMA%npoin
          rga(i) = 0.
          rgb(i) = 0.
          rgc(i) = 0.
          area(i) = 0.
1410   enddo

       do 1420 i=1,AMA%nelem
          xc = (x(lnd(i,1))+x(lnd(i,2))+x(lnd(i,3)))/3
          yc = (y(lnd(i,1))+y(lnd(i,2))+y(lnd(i,3)))/3
          do 1430 j=1,3
             k = lnd(i,j)
             !               rl = ((xc-x(k))**2 + (yc-y(k))**2)**0.5
             rga(k) = rga(k) + ra(i)*tria(i)
             rgb(k) = rgb(k) + ra(i+AMA%melem)*tria(i)
             rgc(k) = rgc(k) + ra(i+2*AMA%melem)*tria(i)
             area(k) = area(k) + tria(i)
1430      enddo
1420   enddo

       do 1440 i=1,AMA%npoin
          rga(i) = rga(i)/area(i)
          rgb(i) = rgb(i)/area(i)
          rgc(i) = rgc(i)/area(i)
          if(i == itest) print *,'$3',rga(i),rgb(i),rgc(i)
1440   enddo
1500 enddo


    !     ... no plotting (too big file)
    if(AMA%npoin .gt. 5000) return

    imt = 18
    open(imt,file='mt',status='unknown')
    do 63 i=1,AMA%npoin
       !     for plotting of the elipses
       if(ibp(i,1) .gt. 0) then
          rga(i) = (rga(i) + rga(ibp(i,1)))/2
          rga(ibp(i,1)) = rga(i)
          rgb(i) = (rgb(i) + rgb(ibp(i,1)))/2
          rgb(ibp(i,1)) = rgb(i)
          rgc(i) = (rgc(i) + rgc(ibp(i,1)))/2
          rgc(ibp(i,1)) = rgc(i)
       endif

       if(mod(i, 8) == 1)   &
            call ELIPS(imt,rga(i),rgb(i),rgc(i),x(i),y(i) )
       !         if(i == itest) print *,'$$',rga(i),rgb(i),rgc(i)
63  enddo
    !      stop
    ! 1313 continue

    deallocate(dx, dy, ra, area, tria)
    return
  end subroutine ERROR1



  subroutine SEEK_BOUN( )
    implicit none
    integer :: i, j, j1, ip, il0, k, l, kl0
    real :: yc, rlen0, rl0
    real :: xp, yp, rlen
    real, dimension(:), pointer :: xb, yb
    integer, dimension(:,:), pointer :: ibb
    integer, dimension(:,:), pointer :: lnd, iae, lbn
    integer, dimension(:), pointer :: ibc, itc
    real, dimension(:), pointer :: x, y

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lbn => AMA%lbn(1:AMA%mbelm,1:2)
    ibc => AMA%ibc(1:AMA%mbelm)
    itc => AMA%itc(1:AMA%mbelm)


    lnd => AMA%lnd(1:AMA%melem, 1:3)
    iae => AMA%iae(1:AMA%melem, 1:3)

    ibb => AMA%ibb(1:AMA%mpoin, 1:3)

    xb => AMA%xb(1:AMA%ipoint)
    yb => AMA%yb(1:AMA%ipoint)

    do i=1,AMA%npoin
       ibb(i,1) = -1
       ibb(i,2) = 0
    enddo


    !     no smooth boundary
    if(AMA%nbp == 0) return

    do i=1,AMA%nelem
       do j=1,3
          if(iae(i,j) .lt. 0) then
             j1 = mod(j,3) + 1
             ip = lnd(i,j)
             !               if(ip == 112) print *,'###'
             xp = x(ip)
             yp = y(ip)
             rlen0 = (xp - x(lnd(i,j1)))*(xp - x(lnd(i,j1))) +   &
                  (yp - y(lnd(i,j1)))*(yp - y(lnd(i,j1)))
             rl0 = 1.E+25*rlen0
             do k=1,AMA%nbp
                do l=AMA%ibpoin(k-1)+1,AMA%ibpoin(k)
                   rlen = (xp - xb(l))*(xp - xb(l)) +  &
                        (yp - yb(l))*(yp - yb(l))
                   if(rlen .lt. rl0) then
                      if(AMA%icrack .gt. 0 .and. xp .gt. -0.01) then
                         yc = (y(lnd(i,1))+y(lnd(i,2))+  &
                              y(lnd(i,3)))/3
                         if(( yc .gt. 0 .and. l .lt. AMA%ibpoin(k)/2)   &
                              .or. ( yc .lt. 0 .and.   &
                              l .gt. AMA%ibpoin(k)/2) )then
                            rl0 = rlen
                            il0 = l
                            kl0 = k
                         endif
                      else
                         rl0 = rlen
                         il0 = l
                         kl0 = k
                      endif
                   endif

                   !                     if(ip  == 2) then
                   !                        write(91,'(4i5,9e10.2)')
                   !     *                       ip,l,il0,kl0,x(ip),y(ip), xb(l),yb(l),
                   !     *                       rlen,rl0,rlen0,rl0/rlen0
                   !                     endif

                enddo
                !               if(ip == 2)
                !                  write(*,'(a3,3i5,6e12.4)')
                !     *                 '###',ip,il0,kl0,x(ip),y(ip),rl0,rlen0,rl0/rlen0
             enddo
             if(rl0/rlen0 .lt. 1E-03) then
                !     the node in 'profile' found
                ibb(ip,1) = il0
                ibb(ip,2) = kl0
             endif
          endif
       enddo
    enddo

    do i=1,AMA%npoin
       if(ibb(i,1) .gt.0 ) then

          !            write(21,'(2e12.4,3i8)')
          !     *        xb(ibb(i,1)),yb(ibb(i,1)),ibb(i,1),ibb(i,2),i

          !     corection of the boundary
          x(i) = xb(ibb(i,1))
          y(i) = yb(ibb(i,1))
       endif
    enddo

    return
  end subroutine SEEK_BOUN

  subroutine COINS(  )
    implicit none
    integer :: i, j, j1, j2, ii, jj, jjj, jj1, jj2, k1, k2, k3, k4, ia1, ia2
    integer :: ib, ja1, ja2, k
    integer, dimension(:,:), pointer :: lnd, iae, lbn
    integer, dimension(:), pointer :: ibc, itc

    lbn => AMA%lbn(1:AMA%mbelm,1:2)
    ibc => AMA%ibc(1:AMA%mbelm)
    itc => AMA%itc(1:AMA%mbelm)

    lnd => AMA%lnd(1:AMA%melem, 1:3)
    iae => AMA%iae(1:AMA%melem, 1:3)

    !     if some triangle has a two times iae(i,j) < 0 the me make a
    !     bascule such that this case doesn't appear

    do 10 i=1,AMA%nelem
       do 20 j=1,3
          j1 = mod(j,3) + 1
          j2 = mod(j1,3) + 1
          if( iae(i,j1) .lt. 0 .and. iae(i,j2) .lt. 0 ) then
             ii = iae(i,j)
             do 30 jjj=1,3
                if(iae(ii,jjj) == i) then
                   jj = jjj
                endif
30           enddo
             jj1 = mod(jj,3) + 1
             jj2 = mod(jj1,3) + 1
             k1 = lnd(i,j)
             k2 = lnd(i,j1)
             k3 = lnd(i,j2)
             k4 = lnd(ii,jj2)
             if(k1 .ne. lnd(ii,jj1) .or. k2 .ne. lnd(ii,jj) ) then
                print *,'ERROR in COINS 1'
             endif
             ia1 = iae(ii,jj1)
             if(ia1 .gt. 0) then
                do 40 k=1,3
                   if(iae(ia1,k) == ii) ja1 = k
40              enddo
             endif
             ia2 = iae(ii,jj2)
             if(ia2 .gt. 0) then
                do 50 k=1,3
                   if(iae(ia2,k) == ii) ja2 = k
50              enddo
             endif
             lnd(i,j1) = k4
             lnd(ii,jj1) = k3
             iae(i,j1) = ii
             iae(ii,jj) = -2
             iae(i,j) = ia1
             iae(ii,jj1) = i
             if(ia1 .gt. 0) iae(ia1,ja1) = i

             do ib=1,AMA%nbelm
                if(itc(ib) == i .and.   &
                     lbn(ib,1) == k2 .and.lbn(ib,2) ==k3)then
                   itc(ib) = ii
                endif
             enddo

          endif
20     enddo
10  enddo
    return
  end subroutine COINS


  subroutine ADJAC( )
    implicit none
    integer :: i, j, k, ii, jj, jj1, j1, il
    integer, dimension(:,:), pointer :: icyc
    integer, dimension(:,:), pointer :: lnd, iae

    lnd => AMA%lnd(1:AMA%melem, 1:3)
    iae => AMA%iae(1:AMA%melem, 1:3)
    icyc => AMA%icyc(1:AMA%mpoin, 1:AMA%maxdeg)

    !     the array icyc has here only the local meaning
    do i=1,AMA%npoin
       icyc(i,1) = 0
    enddo
    do i=1,AMA%nelem
       do j=1,3
          k = lnd(i,j)
          icyc(k,1) = icyc(k,1) + 1
          if(icyc(k,1) .gt. AMA%maxdeg) then
             print *,'Bad dimension in ADJAC'
             print *,'maxdeg < icyc(k,1)',AMA%maxdeg,icyc(k,1)
             print *,k
             stop
          endif
          icyc(k,icyc(k,1)+1) = i
       enddo
    enddo

    do 10 i=1,AMA%nelem
       iae(i,1) = -2
       iae(i,2) = -2
       iae(i,3) = -2
10  enddo
    do 20 i=1,AMA%nelem
       do 30 j=1,3
          if(iae(i,j) == -2)then
             j1 = mod(j,3) +1
             do 40 il=1,icyc(lnd(i,j),1)
                ii = icyc(lnd(i,j),il+1)
                if( ii .ne. i) then
                   do 50 jj=1,3
                      jj1 = mod(jj,3) +1
                      if( (lnd(i,j) == lnd(ii,jj1) ).and.  &
                           (lnd(i,j1) == lnd(ii,jj))) then
                         iae(i,j) = ii
                         iae(ii,jj) = i
                         goto 60
                      endif
50                 enddo
                endif
40           enddo
60           continue
          endif
30     enddo
20  enddo

    !      print *,'The seeking of neighbourhouds done'
    return
  end subroutine ADJAC

  subroutine CYKLE_BOUND( )
    implicit none
    integer :: i, j, j0, ip, l, itest, j1, j2, ii1, ii, ipoc, iost, ilen, imax, jmax
    real :: rmax, x1, y1, x2, y2, x3, y3, u1, v1, u2, v2, rcos
    integer, dimension(:,:), pointer :: icyc
    integer, dimension(:,:), pointer :: lnd, iae
    real, dimension(:), pointer :: x, y

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lnd => AMA%lnd(1:AMA%melem, 1:3)
    iae => AMA%iae(1:AMA%melem, 1:3)

    icyc => AMA%icyc(1:AMA%mpoin, 1:AMA%maxdeg)


    do 10 i=1,AMA%nelem
       do 20 j=1,3
          if(iae(i,j) .le. 0) then
             j0 = j
             j1 = mod(j0,3) + 1
             ip = lnd(i,j0)
             if(icyc(ip,1) .ne. -1) then
                print *,'ERROR in CYKLE_BOUND'
                print *,ip,icyc(ip,1),icyc(ip,2),icyc(ip,3)
                print *,x(ip),y(ip)
                return
             endif
             icyc(ip,2) = lnd(i,j1)
             ii = i
             !     seeking the following point
1000         ilen = abs(icyc(ip,1))
             itest = 0
             do 30 l=1,3
                if(lnd(ii,l) == ip) then
                   itest = 1
                   j0 = l
                   goto 40
                endif
30           enddo
40           continue
             j1 = mod(j0 ,3) + 1
             j2 = mod(j1,3) + 1
             icyc(ip,ilen+2) = lnd(ii,j2)
             icyc(ip,1) = icyc(ip,1) - 1
             if(abs(icyc(ip,1)) +1 .gt. AMA%maxdeg) then
                print *,'The lenght of cykles > maxdeg = ',AMA%maxdeg
                stop
             endif

             ii1 = iae(ii,j2)
             if(ii1 .gt. 0) then
                ii = ii1
                goto 1000
             endif
          endif
20     enddo
10  enddo


    !      test of angles
    ipoc = 0
    if(ipoc == 0) return

    rmax = 10.
    iost = 12134
    open(iost,file='ost',status='unknown')
    do 2000 i=1,AMA%nelem
       do 2010 j=1,3
          j1 = mod(j,3) + 1
          j2 = mod(j1,3) + 1
          x1 = x(lnd(i,j))
          y1 = y(lnd(i,j))
          x2 = x(lnd(i,j1))
          y2 = y(lnd(i,j1))
          x3 = x(lnd(i,j2))
          y3 = y(lnd(i,j2))
          u1 = x2 - x1
          v1 = y2 - y1
          u2 = x3 - x1
          v2 = y3 - y1
          rcos = (u1*u2+v1*v2)/((u1*u1+v1*v1)*(u2*u2+v2*v2))**0.5
          if (rcos .lt. -0.2) then
             write(iost,'(2e14.6)')x1,y1
             write(iost,'(2e14.6)')x2,y2
             write(iost,'(2e14.6)')x3,y3
             write(iost,'(2e14.6)')x1,y1
             write(iost,'(x)')

             ipoc = ipoc + 1
             if(rcos .lt. rmax) then
                rmax = rcos
                imax = i
                jmax = j
             endif
          endif
2010   enddo
2000 enddo
    write( *,'(a12,3i5,2e14.6)')   &
         'the number = ',ipoc,imax,jmax,rmax,acos(rmax)
    print *,'nelem = ',AMA%nelem
    return
  end subroutine CYKLE_BOUND

  subroutine CYKLE( )
    implicit none
    integer :: i, j, k, ip, ii, j1, k1, jj, jj1, jj2, jjj, ibeg1, ibeg2, num1, num2
    integer :: iel, k2
    integer, dimension(:,:), pointer :: ibp
    integer, dimension(:,:), pointer :: icyc
    integer, dimension(:,:), pointer :: lnd, iae, lbn
    integer, dimension(:), pointer :: ibc, itc
    real, dimension(:), pointer :: x, y

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lbn => AMA%lbn(1:AMA%mbelm,1:2)
    ibc => AMA%ibc(1:AMA%mbelm)
    itc => AMA%itc(1:AMA%mbelm)


    lnd => AMA%lnd(1:AMA%melem, 1:3)
    iae => AMA%iae(1:AMA%melem, 1:3)

    icyc => AMA%icyc(1:AMA%mpoin, 1:AMA%maxdeg)
    ibp => AMA%ibp(1:AMA%mpoin, 1:2)


    do 5 i=1,AMA%npoin
       icyc(i,1) = 0
5   enddo

    do 10 k=1,AMA%nbelm
       icyc(lbn(k,1),1) = -1
       icyc(lbn(k,2),1) = -1
10  enddo

    !  ... seeking the  cyclus around each point exept boundaries

    do 30 i=1,AMA%nelem
       do 40 j=1,3
          ip = lnd(i,j)
          if(icyc(ip,1) == 0) then
             !     initialisation
             j1 = mod(j,3) + 1
             k1 = lnd(i,j1)
             ii = i
             jj2 = mod(j1,3) + 1
             icyc(ip,1) = icyc(ip,1) + 1
             icyc(ip,icyc(ip,1) + 1) = k1
             !     repetition
45           icyc(ip,1) = icyc(ip,1) + 1
             icyc(ip,icyc(ip,1) + 1) = lnd(ii,jj2)
             ii = iae(ii,jj2)
             jj = 0
             do jjj=1,3
                if(lnd(ii,jjj) == ip ) jj = jjj
             enddo
             if(jj == 0) then
                print *,'Triangle for node',ip,'doesn''t found'
                print *,'@@',ip, ii,lnd(ii,1),lnd(ii,2),lnd(ii,3)
                print *,'!!',icyc(ip,1),icyc(ip,2),icyc(ip,3),  &
                     icyc(ip,4),icyc(ip,5),icyc(ip,6)
                print *,x(ip),y(ip)
                print *,x(icyc(ip,2)),y(icyc(ip,2))
                print *,x(icyc(ip,3)),y(icyc(ip,3))
                stop
             endif
             jj1 = mod(jj,3) + 1
             jj2 = mod(jj1,3) + 1
             if( lnd(ii,jj2) .ne. k1) goto 45
          endif
40     enddo
30  enddo

    !      print *,'The seeking of cykles done'

    if(AMA%xper(1,1) .gt. 0. .or. AMA%xper(1,2) .gt. 0.) then
       ibeg1 = 0
       ibeg2 = 0
       num1 = 0
       num2 = 0
       do i=1,AMA%nbelm
          if(ibeg1 == 0 .and. ibc(i) == AMA%iper(1,1)) ibeg1 = i
          if(ibeg2 == 0 .and. ibc(i) == AMA%iper(1,2)) ibeg2 = i
          if(ibc(i) == AMA%iper(1,1)) num1 = num1 + 1
          if(ibc(i) == AMA%iper(1,2)) num2 = num2 + 1
          if(ibc(i) == AMA%iper(1,1) .or. ibc(i) == AMA%iper(1,2)) then
             iel = itc(i)
             do j=1,3
                j1 = mod(j,3) +1
                if(lnd(iel,j) == lbn(i,1) .and.   &
                     lnd(iel,j1) == lbn(i,2) ) then
                   iae(iel,j) = -1
                   goto 2005
                endif
             enddo
2005         continue
          endif
       enddo
       if(num1 .ne. num2) then
          print *,'On the periodic boundary is not the same number ',  &
               'of points'
          stop
       endif
       if(num1 == 0) then
          print *,'On the periodic boundary is zero points'
          stop
       endif
       do i=1,num1
          k1 = lbn(ibeg1 -1 + i,1)
          k2 = lbn(ibeg2 + num1 - i,2)
          ibp(k1,1) = k2
          ibp(k2,1) = k1
          ibp(k1,2) = ibp(k1,2) + 1
          ibp(k2,2) = ibp(k2,2) + 1
          if( i == num1) then
             k1 = lbn(ibeg1 -1 + i,2)
             k2 = lbn(ibeg2 + num1 - i,1)
             ibp(k1,1) = k2
             ibp(k2,1) = k1
             ibp(k1,2) = ibp(k1,2) + 1
             ibp(k2,2) = ibp(k2,2) + 1
          endif
       enddo

       !         print *,'The periodic boundary done'
    endif

    !     the second pair of periodic boundaty component
    if(AMA%xper(2,1) .gt. 0. .or. AMA%xper(2,2) .gt. 0.) then
       ibeg1 = 0
       ibeg2 = 0
       num1 = 0
       num2 = 0
       do i=1,AMA%nbelm
          if(ibeg1 == 0 .and. ibc(i) == AMA%iper(2,1)) ibeg1 = i
          if(ibeg2 == 0 .and. ibc(i) == AMA%iper(2,2)) ibeg2 = i
          if(ibc(i) == AMA%iper(2,1)) num1 = num1 + 1
          if(ibc(i) == AMA%iper(2,2)) num2 = num2 + 1
          if(ibc(i) == AMA%iper(2,1) .or. ibc(i) == AMA%iper(2,2)) then
             iel = itc(i)
             do j=1,3
                j1 = mod(j,3) +1
                if(lnd(iel,j) == lbn(i,1) .and.   &
                     lnd(iel,j1) == lbn(i,2) ) then
                   iae(iel,j) = -1
                   goto 2006
                endif
             enddo
2006         continue
          endif
       enddo
       if(num1 .ne. num2) then
          print *,'On the periodic boundary is not the same number ',  &
               'of points'
          stop
       endif
       if(num1 == 0) then
          print *,'On the periodic boundary is zero points'
          stop
       endif
       !     ...ibp(k,2) = 3   ==> point is a periodic for both part of boundary
       do i=1,num1
          k1 = lbn(ibeg1 -1 + i,1)
          k2 = lbn(ibeg2 + num1 - i,2)
          ibp(k1,1) = k2
          ibp(k2,1) = k1
          ibp(k1,2) = ibp(k1,2) + 2
          ibp(k2,2) = ibp(k2,2) + 2
          !            print *,' .',k1,ibp(k1,1),k2,ibp(k2,1)
          if( i == num1) then
             k1 = lbn(ibeg1 -1 + i,2)
             k2 = lbn(ibeg2 + num1 - i,1)
             ibp(k1,1) = k2
             ibp(k2,1) = k1
             ibp(k1,2) = ibp(k1,2) + 2
             ibp(k2,2) = ibp(k2,2) + 2
          endif
       enddo

       !         print *,'The periodic boundary done'
    endif

    return
  end subroutine CYKLE

  subroutine TEST( )
    implicit none
    integer :: i, j, j1, ii, jj, jjj, jj1, k, it, k1, k2, ipp, len
    integer :: ic, ic1, ip1, ip2, itest, jpp, itt, iu, ipoc,idx
    real :: x1, y1, x2, y2, x3, y3, det, xl, yl
    integer, dimension(:,:), pointer :: icyc, ibp
    integer, dimension(:,:), pointer :: lnd, iae, lbn
    integer, dimension(:), pointer :: ibc, itc
    real, dimension(:), pointer :: x, y

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lbn => AMA%lbn(1:AMA%mbelm,1:2)
    ibc => AMA%ibc(1:AMA%mbelm)
    itc => AMA%itc(1:AMA%mbelm)


    lnd => AMA%lnd(1:AMA%melem, 1:3)
    iae => AMA%iae(1:AMA%melem, 1:3)

    icyc => AMA%icyc(1:AMA%mpoin, 1:AMA%maxdeg)
    ibp => AMA%ibp(1:AMA%mpoin, 1:2)


    write(AMA%ifig1,*)'Control of input/output files'

    !                 do i=1,AMA%nelem
    !                     write(*,'(7i5)') i,lnd(I,1),lnd(I,2),lnd(I,3),
    !     *                    iae(i,1),iae(i,2),iae(i,3)
    !                  enddo
    !                  print*,'********************'

    !     ... test of input datas
    do 600 i=1,AMA%nelem
       x1 = x(lnd(i,1))
       y1 = y(lnd(i,1))
       x2 = x(lnd(i,2))
       y2 = y(lnd(i,2))
       x3 = x(lnd(i,3))
       y3 = y(lnd(i,3))
       det = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
       if(det .le. 0. ) then
          print *,'Bad triangle in GT'
          print *,x1,y1,lnd(i,1),det
          print *,x2,y2,lnd(i,2)
          print *,x3,y3,lnd(i,3)
          stop
       endif
       do 610 j=1,3

          j1 = mod(j,3) + 1
          if(lnd(i,j) == lnd(i,j1) ) then
             print *,'ERROR in GT:'
             print *,i,'-th triangle contains the same nodes'
             write(*,'(3i5)') lnd(i,1),lnd(i,2),lnd(i,3)
             stop
          endif
          if(iae(i,j) == iae(i,j1) ) then
             print *,'ERROR in GT:'
             print *,'code found bad neighbours for', i,'-th triangle'
             write(*,'(3i5)') iae(i,1),iae(i,2),iae(i,3)
             stop
          endif
          if(iae(i,j) .gt. 0) then
             ii = iae(i,j)
             jj = 0
             do 620 jjj=1,3
                if(iae(ii,jjj) == i) jj = jjj
620          enddo
             if( jj == 0) then
                print *, 'ERROR in GT:'
                print *, 'the neighbourhouds are bad'
                write(*,'(7i5)') i,lnd(i,1),lnd(i,2),lnd(i,3),  &
                     iae(i,1),iae(i,2),iae(i,3)
                write(*,'(7i5)') ii,lnd(ii,1),lnd(ii,2),lnd(ii,3),  &
                     iae(ii,1),iae(ii,2),iae(ii,3)
                stop
             endif
             jj1 = mod(jj,3) + 1
             if(lnd(i,j) .ne. lnd(ii,jj1) .or.  &
                  lnd(i,j1) .ne. lnd(ii,jj) ) then
                print *,'ERROR in GT: '
                print *, 'the neighbourhouds points are bad'
                write(*,'(7i5)') i,lnd(i,1),lnd(i,2),lnd(i,3),  &
                     iae(i,1),iae(i,2),iae(i,3)
                write(*,'(7i5)') ii,lnd(ii,1),lnd(ii,2),lnd(ii,3),  &
                     iae(ii,1),iae(ii,2),iae(ii,3)
                stop
             endif
          endif
610    enddo
600 enddo

    !     checking of boundary
    do k=1,AMA%nbelm
       it = itc(k)
       k1 = lbn(k,1)
       k2 = lbn(k,2)

       if(it .gt. AMA%nelem) then
          print *,'ERROR in GT:'
          print *,'too high number itc'
          print *,k,it,k1,k2,AMA%nelem
          stop
       endif

       if((lnd(it,1) == k1 .and. lnd(it,2) == k2 .and.   &
            iae(it,1) .lt. 0) .or.  &
            (lnd(it,2) == k1 .and. lnd(it,3) == k2 .and.   &
            iae(it,2) .lt. 0) .or.  &
            (lnd(it,3) == k1 .and. lnd(it,1) == k2 .and.   &
            iae(it,3) .lt. 0))  then
          !     OK
       else
          print *,'ERROR in GT:'
          print *,k,'-th boundary segment is bad'
          print *,'lbn(k,1) and lbn(k,2) don''t coresponds with',  &
               'the field lnd(i,j)'
          print *,k,k1,k2
          stop
       endif
    enddo


    !     checking of the  cykles
    do 900 ipp =1,AMA%npoin-1
       if( icyc(ipp,1) .gt. 0 ) then
          len = icyc(ipp,1)
       else
          len = abs( icyc(ipp,1)) - 1
       endif
       do 910 ic = 1,len
          ic1 = mod(ic,abs(icyc(ipp,1)) ) + 1
          ip1 = icyc(ipp,ic+1)
          ip2 = icyc(ipp,ic1+1)
          itest = 0
          do 920 jpp=1,AMA%nelem
             if(( lnd(jpp,1) == ipp .and.  &
                  lnd(jpp,2) == ip1 .and.  &
                  lnd(jpp,3) == ip2 ) .or.  &
                  ( lnd(jpp,2) == ipp .and.  &
                  lnd(jpp,3) == ip1 .and.  &
                  lnd(jpp,1) == ip2 ) .or.  &
                  ( lnd(jpp,3) == ipp .and.  &
                  lnd(jpp,1) == ip1 .and.  &
                  lnd(jpp,2) == ip2 ) ) then
                itest = 1
                goto 922
             endif
920       enddo
922       continue
          if(itest == 0)then
             print *,'ERROR in CHECk in GT'
             print *,ipp,ip1,ip2,ic,ic1
             it = ipp
             print *,x(it),y(it)
             write(*,'(10i5)') it,icyc(it,1),  &
                  icyc(it,2),  &
                  icyc(it,3),icyc(it,4),icyc(it,5),  &
                  icyc(it,6),icyc(it,7),  &
                  icyc(it,8),icyc(it,9)
             do 912 itt=1,abs(icyc(it,1))
                print *,x(icyc(it,itt+1)),  &
                     y(icyc(it,itt+1)),icyc(it,itt+1)
912          enddo
             do iu=1,AMA%nelem
                if(lnd(iu,1) == it .or. lnd(iu,2) == it  &
                     .or. lnd(iu,3) == it) then
                   write(*,'(4i5)') iu,lnd(iu,1),lnd(iu,2),  &
                        lnd(iu,3)
                endif
             enddo
             stop
          endif
910    enddo
900 enddo

    if(AMA%xper(1,1) .gt. 0 .or. AMA%xper(1,2) .gt. 0) then
       ipoc = 0
       do i=1,AMA%npoin
          if(ibp(i,1) .gt. 0 .and. ibp(i,2) == 1) then
             ipoc = ipoc + 1
             if(i .ne. ibp(ibp(i,1),1) ) then
                print *,'ERROR in GT:'
                print *,'The periodic nodes in triang are ',  &
                     'not periodical-1:'
                print *,i,ibp(i,1)
                print *,ibp(i,1),ibp(ibp(i,1),1)
                print *,x(i),y(i)
                print *,x(ibp(i,1)),y(ibp(i,1))
                print *,x(ibp(ibp(i,1),1)),y(ibp(ibp(i,1),1))
                stop
             endif

             xl = ( abs(x(i) - x(ibp(i,1))) ) - AMA%xper(1,1)
             yl = ( abs(y(i) - y(ibp(i,1))) ) - AMA%xper(1,2)
             if( xl .gt. 1E-05 .or. yl .gt. 1E-05) then
                print *,'ERROR in GT, bad periodical boundary-1'
                print *,i,ibp(i,1),ibp(ibp(i,1),1)
                print *,x(i),y(i)
                print *,x(ibp(i,1)),y(ibp(i,1))
                print *,AMA%xper(1,1),AMA%xper(1,2),xl,yl
                stop
             endif
          endif
       enddo
       if(ipoc == 0) then
          print *,'ERROR in GT:'
          print *,'non periodic point!!!-1'
          stop
       endif
    endif

    if(AMA%xper(2,1) .gt. 0 .or. AMA%xper(2,2) .gt. 0) then
       ipoc = 0
       do i=1,AMA%npoin
          if(ibp(i,1) .gt. 0 .and. ibp(i,2) == 2) then
             ipoc = ipoc + 1
             if(i .ne. ibp(ibp(i,1),1) ) then
                print *,'ERROR in GT:'
                print *,'The periodic nodes in triang are ',  &
                     'not periodical-2:'
                print *,i,ibp(i,1)
                print *,ibp(i,1),ibp(ibp(i,1),1)
                print *,x(i),y(i)
                print *,x(ibp(i,1)),y(ibp(i,1))
                print *,x(ibp(ibp(i,1),1)),y(ibp(ibp(i,1),1))
                stop
             endif

             xl = ( abs(x(i) - x(ibp(i,1))) ) - AMA%xper(2,1)
             yl = ( abs(y(i) - y(ibp(i,1))) ) - AMA%xper(2,2)
             if( xl .gt. 1E-05 .or. yl .gt. 1E-05) then
                print *,'ERROR in GT, bad periodical boundary-2'
                print *,i,ibp(i,1),ibp(ibp(i,1),1)
                print *,x(i),y(i)
                print *,x(ibp(i,1)),y(ibp(i,1))
                print *,AMA%xper(2,1),AMA%xper(2,2),xl,yl
                stop
             endif
          endif
       enddo
       if(ipoc == 0) then
          print *,'ERROR in GT:'
          print *,'non periodic point!!!-2'
          stop
       endif

       ! corner points
       ipoc = 0
       do i=1,AMA%npoin
          if(ibp(i,1) .gt. 0 .and. ibp(i,2) == 3) then
             ipoc = ipoc + 1
             if(i .ne. ibp(ibp(i,1),1) ) then
                print *,'ERROR in GT:'
                print *,'The periodic nodes in triang are ',  &
                     'not periodical-3:'
                print *,i,ibp(i,1)
                print *,ibp(i,1),ibp(ibp(i,1),1)
                print *,x(i),y(i)
                print *,x(ibp(i,1)),y(ibp(i,1))
                print *,x(ibp(ibp(i,1),1)),y(ibp(ibp(i,1),1))
                stop
             endif

             xl = ( abs(x(i) - x(ibp(i,1))) ) - AMA%xper(2,1)
             yl = ( abs(y(i) - y(ibp(i,1))) ) - AMA%xper(2,2)
             if( xl .gt. 1E-05 .or. yl .gt. 1E-05) then
                print *,'ERROR in GT, bad periodical boundary-3'
                print *,i,ibp(i,1),ibp(ibp(i,1),1)
                print *,x(i),y(i)
                print *,x(ibp(i,1)),y(ibp(i,1))
                print *,AMA%xper(2,1),AMA%xper(2,2),xl,yl
                stop
             endif
          endif
       enddo
       if(ipoc == 0) then
          print *,'ERROR in GT:'
          print *,'non periodic point!!!-2'
          stop
       endif
    endif
    write(AMA%ifig1,*)'Input/output file is OK'
    return
  end subroutine TEST


  subroutine CYKLE_REP( )
    implicit none
    integer :: i, j, k, ip, j1, k1, ii, jj, jj1, jj2, jjj
    integer, dimension(:,:), pointer :: icyc
    integer, dimension(:,:), pointer :: lnd, iae, lbn
    integer, dimension(:), pointer :: ibc, itc

    lbn => AMA%lbn(1:AMA%mbelm,1:2)
    ibc => AMA%ibc(1:AMA%mbelm)
    itc => AMA%itc(1:AMA%mbelm)


    lnd => AMA%lnd(1:AMA%melem, 1:3)
    iae => AMA%iae(1:AMA%melem, 1:3)

    icyc => AMA%icyc(1:AMA%mpoin, 1:AMA%maxdeg)

    do 5 i=1,AMA%npoin
       icyc(i,1) = 0
5   enddo

    do 10 k=1,AMA%nbelm
       icyc(lbn(k,1),1) = -1
       icyc(lbn(k,2),1) = -1
10  enddo

    !  ... seeking the  cyclus around each point exept boundaries

    do 30 i=1,AMA%nelem
       do 40 j=1,3
          ip = lnd(i,j)
          if(icyc(ip,1) == 0) then
             !     initialisation
             j1 = mod(j,3) + 1
             k1 = lnd(i,j1)
             ii = i
             jj2 = mod(j1,3) + 1
             icyc(ip,1) = icyc(ip,1) + 1
             icyc(ip,icyc(ip,1) + 1) = k1
             !     repetition
45           icyc(ip,1) = icyc(ip,1) + 1
             if(icyc(ip,1) +1 .gt. AMA%maxdeg) then
                print *,'The lenght of cykles > maxdeg = ', AMA%maxdeg
                stop
             endif
             icyc(ip,icyc(ip,1) + 1) = lnd(ii,jj2)
             ii = iae(ii,jj2)
             jj = 0
             do jjj=1,3
                if(lnd(ii,jjj) == ip ) jj = jjj
             enddo
             if(jj == 0) then
                print *,'Triangle for node',ip,'doesn''t found'
                stop
             endif
             jj1 = mod(jj,3) + 1
             jj2 = mod(jj1,3) + 1
             if( lnd(ii,jj2) .ne. k1) goto 45
          endif
40     enddo
30  enddo
    return
  end subroutine CYKLE_REP


  subroutine QUALITY(ndim, err, rminerr, imin, rmaxrez)
    implicit none
    integer, intent(in) :: ndim
    integer, intent(inout) :: imin
    real, intent(inout) :: err, rminerr, rmaxrez
    integer :: ice, ipobound, ibo, i, i1, i2, i3, j
    real :: x1, x2, x3, y1, y2, y3, w1, w2, w3
    real :: ra1, rb1, rc1, ra2, rb2, rc2, ra3, rb3, rc3, err1, err2
    real, dimension(:), pointer :: rga, rgb, rgc
    real, dimension(:,:), pointer :: wp
    integer, dimension(:,:), pointer :: lnd, iae
    real, dimension(:), pointer :: x, y

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lnd => AMA%lnd(1:AMA%melem, 1:3)
    iae => AMA%iae(1:AMA%melem, 1:3)

    wp    => AMA%wp(   1:AMA%mpoin,1:ndim+1)

    rga => AMA%rga( 1:AMA%mpoin )
    rgb => AMA%rgb( 1:AMA%mpoin )
    rgc => AMA%rgc( 1:AMA%mpoin )


    ice = 0
    err = 0.
    rminerr = 1.
    AMA%errrez = 0.
    rmaxrez = 0.
    AMA%glmax = 0.
    AMA%glmin = 1000.
    ipobound = 0
    ibo = 0
    do 10 i=1,AMA%nelem
       i1 = lnd(i,1)
       i2 = lnd(i,2)
       i3 = lnd(i,3)
       x1 = x(i1)
       y1 = y(i1)
       x2 = x(i2)
       y2 = y(i2)
       x3 = x(i3)
       y3 = y(i3)
       w1 = wp(i1,1)
       w2 = wp(i2,1)
       w3 = wp(i3,1)
       ra1 = rga(i1)
       rb1 = rgb(i1)
       rc1 = rgc(i1)
       ra2 = rga(i2)
       rb2 = rgb(i2)
       rc2 = rgc(i2)
       ra3 = rga(i3)
       rb3 = rgb(i3)
       rc3 = rgc(i3)


       do j=1,3
          if( iae(i,j) .lt. -1) then
             ibo = j
             ipobound = ipobound + 1
          endif
       enddo

       call QUA2(x1, y1, w1, ra1, rb1, rc1, x2, y2, w2, ra2, rb2, rc2,  &
            x3, y3, w3, ra3, rb3, rc3, err1, err2, ice, ibo)

       ice = 0
       err = err + err1
       AMA%errrez = AMA%errrez + err2
       if(rminerr .gt. err1 ) then
          rminerr = err1
          imin = i
       endif
       rmaxrez = max(rmaxrez,err2)

       if( i == -1 )   &
            write(*,'(a12,i5,6es12.4)') 'ama-restA',i,AMA%errrez,err2,rmaxrez

10  enddo
    err = err/AMA%nelem
    AMA%errrez = AMA%errrez/(AMA%nelem*3 + ipobound)


    if(AMA%ifig .ge. 0) then
       call PLOT1()
       AMA%ifig = AMA%ifig + 1
    endif

    return
  end subroutine QUALITY

  subroutine PLOT1( )
    implicit none
    integer :: igra, i
    character*10 meshfile
    character*1 num
    character*2 nux
    character*3 nuy
    character*4 nuz
    character*5 nuw
    integer, dimension(:,:), pointer :: lnd
    real, dimension(:), pointer :: x, y

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lnd => AMA%lnd(1:AMA%melem, 1:3)

    meshfile(1:10) = 'mesh      '
    if( AMA%ifig .lt. 10 ) then
       write( num, '(i1)' ) AMA%ifig
       meshfile(5:5) =  num
    elseif( AMA%ifig .lt. 100 ) then
       write( nux, '(i2)' ) AMA%ifig
       meshfile(5:6) =  nux
    elseif( AMA%ifig .lt. 1000 ) then
       write( nuy, '(i3)' ) AMA%ifig
       meshfile(5:7) =  nuy
    elseif( AMA%ifig .lt. 10000 ) then
       write( nuz, '(i4)' ) AMA%ifig
       meshfile(5:8) =  nuz
    elseif( AMA%ifig .lt. 100000 ) then
       write( nuw, '(i5)' ) AMA%ifig
       meshfile(5:9) =  nuw
    else
       print *, ' Number of gnu-files has to be less then 99999'
       stop
    endif

    !      ... plot of a new mesh

    igra = 53
    open(igra, file=meshfile,status='UNKNOWN')

    do i=1,AMA%nelem
       write(igra,'(2e14.6)') x(lnd(i,1)),y(lnd(i,1))
       write(igra,'(2e14.6)') x(lnd(i,2)),y(lnd(i,2))
       write(igra,'(2e14.6)') x(lnd(i,3)),y(lnd(i,3))
       write(igra,'(2e14.6)') x(lnd(i,1)),y(lnd(i,1))
       write(igra,'(20x)')
    enddo
    close(igra)

  end subroutine PLOT1

  subroutine PLOT( )
    implicit none
    integer :: igra, i
    integer, dimension(:,:), pointer :: lnd
    real, dimension(:), pointer :: x, y

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lnd => AMA%lnd(1:AMA%melem, 1:3)

    !      ... plot of a new mesh
    igra = 53
    open(igra, file='mesh',status='UNKNOWN')
    !open(igra, file='mesh',status='replace')

    do i=1,AMA%nelem
       write(igra,'(2e14.6)') x(lnd(i,1)),y(lnd(i,1))
       write(igra,'(2e14.6)') x(lnd(i,2)),y(lnd(i,2))
       write(igra,'(2e14.6)') x(lnd(i,3)),y(lnd(i,3))
       write(igra,'(2e14.6)') x(lnd(i,1)),y(lnd(i,1))
       write(igra,'(20x)')
    enddo
    if(AMA%xper(1,1) .gt. 0 .or. AMA%xper(1,2) .gt. 0) then
       do  i=1,AMA%nelem
          write(igra,'(2e14.6)') x(lnd(i,1))+AMA%xper(1,1),y(lnd(i,1))+AMA%xper(1,2)
          write(igra,'(2e14.6)') x(lnd(i,2))+AMA%xper(1,1),y(lnd(i,2))+AMA%xper(1,2)
          write(igra,'(2e14.6)') x(lnd(i,3))+AMA%xper(1,1),y(lnd(i,3))+AMA%xper(1,2)
          write(igra,'(2e14.6)') x(lnd(i,1))+AMA%xper(1,1),y(lnd(i,1))+AMA%xper(1,2)
          write(igra,'(20x)')
       enddo
    endif

    if(AMA%xper(2,1) .gt. 0 .or. AMA%xper(2,2) .gt. 0) then
       do  i=1,AMA%nelem
          write(igra,'(2e14.6)') x(lnd(i,1))+AMA%xper(2,1),y(lnd(i,1))+AMA%xper(2,2)
          write(igra,'(2e14.6)') x(lnd(i,2))+AMA%xper(2,1),y(lnd(i,2))+AMA%xper(2,2)
          write(igra,'(2e14.6)') x(lnd(i,3))+AMA%xper(2,1),y(lnd(i,3))+AMA%xper(2,2)
          write(igra,'(2e14.6)') x(lnd(i,1))+AMA%xper(2,1),y(lnd(i,1))+AMA%xper(2,2)
          write(igra,'(20x)')
       enddo
    endif

    close(igra)

  end subroutine PLOT


  subroutine WriteMetrix( )
    implicit none
    integer :: i, j, ip
    real :: r_max, r_min
    real, dimension(:), pointer :: rga, rgb, rgc
    integer, dimension(:,:), pointer :: lnd
    real, dimension(:), pointer :: x, y

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lnd => AMA%lnd(1:AMA%melem, 1:3)

    rga => AMA%rga( 1:AMA%mpoin )
    rgb => AMA%rgb( 1:AMA%mpoin )
    rgc => AMA%rgc( 1:AMA%mpoin )

    print*,'####,  Writing file "fort.',AMA%ifig2
    !      do i=1,AMA%npoin
    !         r_max = max(rga(i), rgc(i) )
    !         r_min = min(rga(i), rgc(i) )
    !
    !         if(r_max > 1)
    !     *        write(AMA%ifig2, *) x(i), y(i), r_max, r_max/r_min
    !
    !      enddo


    do i=1,AMA%nelem

       do j=1,4
          ip = lnd(i, mod(j,3)+1)
          r_max = max(rga(ip), rgc(ip) )
          r_min = min(rga(ip), rgc(ip) )

          write(AMA%ifig2, *)   &
               x(ip), y(ip), r_max, r_max/r_min
       enddo
       write(AMA%ifig2, '(x)')
       write(AMA%ifig2, '(x)')
       write(AMA%ifig2, '(x)')


    enddo
    AMA%ifig2 = AMA%ifig2 + 1

  end subroutine WriteMetrix


  subroutine MOVING( icy  )
    implicit none
    integer, intent(in) ::  icy
    real, dimension(:), allocatable :: xloc, yloc
    integer, dimension(:), allocatable :: iloc

    integer :: i, j, j1, ii1, ii2, k, l, ilo, i0, i1, i2, numit, ierr, ice, len
    integer :: il0, il1, il2, ilk, ikk1, imov, iprint, itet, jbet, kk1, len1
    real :: det, quality1, qualityE, quality1rez,  qualityrez, xi, yi, xi1, yi1
    real :: a, b, c, d, x1, x2, y1, y2, xs, ys
    real :: rlen1, rlen2, ss, rcos, yq, xq, xp, yp, xnew, ynew
    real :: xold, yold, acc, rl1, rkg, x0, y0, xperreal, yperreal
    real, dimension(:), pointer :: xb, yb
    integer, dimension(:,:), pointer :: ibb, ibp
    real, dimension(:), pointer :: rga, rgb, rgc
    integer, dimension(:,:), pointer :: icyc
    integer, dimension(:,:), pointer :: lnd, iae, lbn
    integer, dimension(:), pointer :: ibc, itc
    real, dimension(:), pointer :: x, y
    integer :: ic_start, ic_end, ic_skip

    x => AMA%x(1:AMA%mpoin)
    y => AMA%y(1:AMA%mpoin)

    lbn => AMA%lbn(1:AMA%mbelm,1:2)
    ibc => AMA%ibc(1:AMA%mbelm)
    itc => AMA%itc(1:AMA%mbelm)


    lnd => AMA%lnd(1:AMA%melem, 1:3)
    iae => AMA%iae(1:AMA%melem, 1:3)

    icyc => AMA%icyc(1:AMA%mpoin, 1:AMA%maxdeg)

    rga => AMA%rga( 1:AMA%mpoin )
    rgb => AMA%rgb( 1:AMA%mpoin )
    rgc => AMA%rgc( 1:AMA%mpoin )

    ibp => AMA%ibp(1:AMA%mpoin, 1:2)
    ibb => AMA%ibb(1:AMA%mpoin, 1:3)

    xb => AMA%xb(1:AMA%ipoint)
    yb => AMA%yb(1:AMA%ipoint)


    allocate(xloc(1:25), yloc(1:25), iloc(1:25) )

    if(mod(icy, 2) == 0) then
       ic_start = 1;         ic_end = AMA%npoin;   ic_skip =  1
    else
       ic_start = AMA%npoin;  ic_end = 1;          ic_skip = -1
    endif

    !do 10 i=1,AMA%npoin
    do 10 i = ic_start, ic_end, ic_skip

       if(ibb(i,3) == -1) goto 10
       if(icyc(i,1) .gt. 0) then
          if(icyc(i,1) .gt. 25) then
             print *,'Error in dimension in MOVING',icyc(i,1)
             return
          endif
          xp = x(i)
          yp = y(i)
          xq = 0.
          yq = 0.
          qualityE = 0.
          do 20 j=1,icyc(i,1)
             j1 = mod(j,icyc(i,1) ) +1
             ii1 = icyc(i,j+1)
             ii2 = icyc(i,j1+1)
             xi = x(ii1)
             yi = y(ii1)
             xi1 = x(ii2)
             yi1 = y(ii2)
             a = (rga(ii1) + rga(ii2))/2
             b = (rgb(ii1) + rgb(ii2))/2
             c = (rgc(ii1) + rgc(ii2))/2
             d = ( (a*c-b*b)/3)**0.5
             if( d .le. 0) then
                print *,'ERROR in matrix in moving'
                print *,a,b,c
                stop
             endif

             xs = xi + ((d-b)*(xi1-xi) - c*(yi1-yi))/2/d
             ys = yi + ( a*(xi1-xi) + (d+b)*(yi1-yi))/2/d

             xloc(j) = xs
             yloc(j) = ys

             xq = xq + xloc(j)
             yq = yq + yloc(j)
20        enddo

          xq = xq/icyc(i,1)
          yq = yq/icyc(i,1)

          !     now we check the the triangles xi,xi1,xq are reals
          numit = 10
          k=numit
1240      continue
          ierr = 0
          xnew = xp + 1.*k/numit*(xq - xp)
          ynew = yp + 1.*k/numit*(yq - yp)
          do 1230 j=1,icyc(i,1)
             j1 = mod(j,icyc(i,1))+1
             ii1 = icyc(i,j+1)
             ii2 = icyc(i,j1+1)
             xi = x(ii1)
             yi = y(ii1)
             xi1 = x(ii2)
             yi1 = y(ii2)
             !     check the orientation and bad conditionality
             det = POS_TEST(xi, yi, xi1, yi1, xnew, ynew )
             !print *,'..MOVa.',det

             if( det .le. 1.) then
                ierr = 1
                goto 1231
             endif
             call POS1TEST(xi,yi,xi1,yi1,xnew,ynew,itet)
             if(itet == 1) then
                ierr = 1
                goto 1231
             endif

1230      enddo
1231      continue

          if(ierr == 1 .and. k .ge. 2) then
             k = k - 1
             goto 1240
          endif
          if( k == 1) then
             xq = xp
             yq = yp
          else
             xq = xnew
             yq = ynew
          endif

          qualityE = 0.
          qualityrez = 1E+35
          numit = 8
          ice = 0
          xold = xp
          yold = yp

          do 40 k=0,numit
             quality1 = 0.
             quality1rez = 0.
             xnew = xp + 1.*k/numit*(xq - xp)
             ynew = yp + 1.*k/numit*(yq - yp)
             do 30 j=1,icyc(i,1)
                j1 = mod(j,icyc(i,1))+1
                ii1 = icyc(i,j+1)
                ii2 = icyc(i,j1+1)
                xi = x(ii1)
                yi = y(ii1)
                xi1 = x(ii2)
                yi1 = y(ii2)

                !     check the orientation and bad conditionality
                det = POS_TEST(xi, yi, xi1, yi1, xnew, ynew )
                !print *,'..MOVb.',det

                if( det .le. 1.) then
                   ! print *,'VIOLATION of positivity - algorith ERROR'
                   goto 45
                endif
                call POS1TEST(xi,yi,xi1,yi1,xnew,ynew,itet)
                if(itet == 1) then
                   goto 45
                endif

                acc = ACCUTE(xi,yi,xi1,yi1,xnew,ynew, 1.)

                rl1 = (0.5*((rga(i)+rga(ii1))*(xi-xnew)*(xi-xnew)+  &
                     2*(rgb(i)+rgb(ii1))*(xi-xnew)*(yi-ynew)+  &
                     (rgc(i)+rgc(ii1))*(yi-ynew)*(yi-ynew)))**0.5

                quality1rez = quality1rez + (rl1 - 1.732050807)**2
                !     ... NEW ACCUTE
                !                  quality1rez = quality1rez
                !     *                 +(rl1 - 1.732050807)**2*acc

30           enddo

             if(quality1rez .ge. qualityrez ) goto 45

             xold = xnew
             yold = ynew
             qualityrez = quality1rez
40        enddo
45        continue

          x(i) = xold
          y(i) = yold
       endif
10  enddo

    !     for boundary points, the movement, only if the point is on the
    !     straight segment
    !     or on the boundary given by the field [xb(i),yb(i)] i=1,ipoint

    do 100 i=1,AMA%npoin

       if(ibb(i,3) == -1) goto 100
       if(ibb(i,1) .gt. 0) then
          ikk1 = ibb(i,2)
          if(ibb(i,1) == AMA%ibpoin(ikk1-1)+1  .or.  &
               ibb(i,1) == AMA%ibpoin(ikk1)) then
             !               NO moving of final or initial node of profiles
             goto 21
          endif
       endif

       if(icyc(i,1) .lt. 0 ) then
          len = abs(icyc(i,1))
          x1 = x(i)
          y1 = y(i)
          x2 = x(icyc(i,2))
          y2 = y(icyc(i,2))
          x0 = x(icyc(i,len+1))
          y0 = y(icyc(i,len+1))
          rlen1 = ((x2-x1)**2 + (y2-y1)**2)**0.5
          rlen2 = ((x0-x1)**2 + (y0-y1)**2)**0.5
          ss = ((x2-x1)*(x0-x1) + (y2-y1)*(y0-y1))
          rcos = ss/rlen1/rlen2
          !            if(i == 2) print *,'%%%%',rcos
            if( rcos .gt. -0.9999 .or. (ibb(i,1) .gt. 0 .and. &
                 !            if( rcos .gt. -0.9866 .or. (ibb(i,1) .gt. 0 .and.  &
                 (ibb(icyc(i,2),1) .gt. 0 .and.   &
                 ibb(icyc(i,len+1),1) .lt. 0 ) .or.   &
                 (ibb(icyc(i,2),1) .lt. 0 .and.   &
                 ibb(icyc(i,len+1),1) .gt. 0 ) )) then
!     too sharp angle or. begin or end of noclosed profile,
!     movement forbidden
               goto 100
            elseif ( ibb(i,1) .lt. 0  ) then
!     this points are on straight segments
!     the maximal value of numit can not be greater than 9 !!!!!!
               numit = 9
               ice = 0
               xold = x1
               yold = y1
               do 141 k=-numit,0
                  l = numit + k + 1
                  xloc(l) = x1 - 1.*k/(numit+1)*(x0 - x1)
                  yloc(l) = y1 - 1.*k/(numit+1)*(y0 - y1)
                  iloc(l) = -1
141            enddo
               do 146 k=1,numit
                  l = numit + k + 1
                  xloc(l) = x1 + 1.*k/(numit+1)*(x2 - x1)
                  yloc(l) = y1 + 1.*k/(numit+1)*(y2 - y1)
                  iloc(l) = -1
146            enddo
               ilo = 2*numit + 1
            else
               !     nodes are on curved part of boundary
               i2 = icyc(i,2)
               i0 = icyc(i,abs(icyc(i,1))+1)

               if(ibb(i,1) .lt. 0 .or. ibb(i0,1) .lt. 0   &
                    .or. ibb(i2,1) .lt.0) then
                  print *,'LOGICAL error in MOVING, ibb'
                  print *,i,i0,i2,ibb(i,1),ibb(i0,1),ibb(i2,1)
                  print *,x(i),y(i)  !,xb(ibb(i,1)),yb(ibb(i,1))
                  print *,x(i0),y(i0)  !,xb(ibb(i0,1)),yb(ibb(i0,1))
                  print *,x(i2),y(i2) !,xb(ibb(i2,1)),yb(ibb(i2,1))
                  stop
               endif

               il0 = ibb(i0,1)
               il1 = ibb(i,1)
               il2 = ibb(i2,1)
               ilk = ibb(i,2)

               if(ilk.ne. ibb(i0,2) .or. ilk .ne. ibb(i2,2) ) then
                  print *,'error jkol'
                  print *,x(i),y(i),i,ibb(i,1),ibb(i,2),ilk
                  print *,x(i0),y(i0),i0,ibb(i0,1),ibb(i0,2)
                  print *,x(i2),y(i2),i2,ibb(i2,1),ibb(i2,2)
                  stop
               endif

!     reordering in the case of end point
               if(il0 == AMA%ibpoin(ilk-1)+1 .and. il1 .gt. il2) then
                  il0 = AMA%ibpoin(ilk)
                  if(abs(xb(AMA%ibpoin(ilk-1)+1)-xb(AMA%ibpoin(ilk)) )   &
                       .gt. 1E-05 .or.  &
                       abs(yb(AMA%ibpoin(ilk-1)+1)-yb(AMA%ibpoin(ilk)) )   &
                       .gt. 1E-05 )then
                     print *,'error in profile, the end point is not'
                     print *,'the first, or bad reordering'
                     print *,ilk,(AMA%ibpoin(ilk-1)+1),AMA%ibpoin(ilk)
                     print *,'1',il0,il1,il2
                     stop
                  endif
               endif

               if(il0 == AMA%ibpoin(ilk) .and. il2 .gt. il1) then
                  il0 = AMA%ibpoin(ilk-1)+1
                  if(abs(xb(AMA%ibpoin(ilk-1)+1)-xb(AMA%ibpoin(ilk)) )   &
                       .gt. 1E-05 .or.  &
                       abs(yb(AMA%ibpoin(ilk-1)+1)-yb(AMA%ibpoin(ilk)) )   &
                       .gt. 1E-05 )then
                     print *,'error in profile, the end point is not'
                     print *,'the first, or bad reordering'
                     print *,'2',il0,il1,il2
                     stop
                  endif
               endif

               if(il2 == AMA%ibpoin(ilk-1)+1 .and. il1 .gt. il0) then
                  il2 = AMA%ibpoin(ilk)
                  if(abs(xb(AMA%ibpoin(ilk-1)+1)-xb(AMA%ibpoin(ilk)) )   &
                       .gt. 1E-05 .or.  &
                       abs(yb(AMA%ibpoin(ilk-1)+1)-yb(AMA%ibpoin(ilk)) )   &
                       .gt. 1E-05 )then
                     print *,'error in profile, the end point is not'
                     print *,'the first, or bad reordering'
                     print *,'3',il0,il1,il2
                     stop
                  endif
               endif

               if(il2 == AMA%ibpoin(ilk) .and. il0 .gt. il1) then
                  il2 = AMA%ibpoin(ilk-1)+1
                  if(abs(xb(AMA%ibpoin(ilk-1)+1)-xb(AMA%ibpoin(ilk)) )   &
                       .gt. 1E-05 .or.  &
                       abs(yb(AMA%ibpoin(ilk-1)+1)-yb(AMA%ibpoin(ilk)) )   &
                       .gt. 1E-05 )then
                     print *,'error in profile, the end point is not'
                     print *,'the first, or bad reordering'
                     print *,'4',il0,il1,il2
                     print *,x(i0),y(i0),i0
                     print *,x(i),y(i),i
                     print *,x(i2),y(i2),i2
                     stop
                  endif
               endif
               !     always il0 < il1 < il2 due to the orientation
               !     but may be il1 or il2 > AMA%ibpoin(ilk)
               if( il2 .gt. il0 .and. il2 - il0 .le. 20) then
                  ilo = il2 - il0 - 1
                  do 210 l=il0+1,il2-1
                     xloc(l-il0) = xb(l)
                     yloc(l-il0) = yb(l)
                     iloc(l-il0) = l
210               enddo
               endif
               if(il0 .gt. il2 .and. il2+AMA%ibpoin(ilk)-il0-2 .le. 20) then
                  ilo = il2 + AMA%ibpoin(ilk) - il0 -2
                  do 220 l=il0+1,AMA%ibpoin(ilk)-1
                     xloc(l-il0) = xb(l)
                     yloc(l-il0) = yb(l)
                     iloc(l-il0) = l
220               enddo
                  do l=1,il2-1
                     xloc(l+AMA%ibpoin(ilk) - il0 -1) = xb(l)
                     yloc(l+AMA%ibpoin(ilk) - il0 -1) = yb(l)
                     iloc(l+AMA%ibpoin(ilk) - il0 -1) = l
                  enddo
               endif

               if(il2 .gt. il0 .and. il2 - il0 .gt. 20) then
                  rkg = 1.*(il2-il0 -2) /19
                  do 230 k=0,19
                     l = il0 + 1 +int(k*rkg+0.5)
                     xloc(k+1) = xb(l)
                     yloc(k+1) = yb(l)
                     iloc(k+1) = l
230               enddo
                  ilo = 20
               endif
               if(il0.gt.il2 .and. il2+AMA%ibpoin(ilk)-il0-2 .gt. 20) then
                  rkg = 1.*(il2-il0+AMA%ibpoin(ilk) -2) /19
                  do 240 k=0,19
                     l = il2 + 1 +int(k*rkg+0.5)
                     if(l .ge. AMA%ibpoin(ilk)) l = l -AMA%ibpoin(ilk)+1
                     iprint = 0
                     if(iprint == 1) then
                        print *,'&&&',i,l,AMA%ibpoin(ilk)
                        print *,x(i),y(i),i
                        print *,x(i0),y(i0),i0
                        print *,x(i2),y(i2),i2
                        print *,'______________________'
                        do kk1 =1,abs(icyc(i,1))
                           print *,x(icyc(i,kk1+1)),y(icyc(i,kk1+1)),  &
                                icyc(i,kk1+1)
                        enddo
                        print *,'______________________'
                        print *,xb(il0),yb(il0),il0
                        print *,xb(il1),yb(il1),il1
                        print *,xb(il2),yb(il2),il2
                        print *,'Mischmatch detected in MOVING.f????'
                     endif
                     xloc(k+1) = xb(l)
                     yloc(k+1) = yb(l)
                     iloc(k+1) = l
240               enddo
                  ilo = 20
                  !                  stop
               endif
            endif

            qualityE = 0.
            qualityrez = 1E+35

            if(ibp(i,1) == 0 ) then
!     for nonperiodical points
               jbet = 0
               do 140 k=1,ilo
                  quality1 = 0.
                  quality1rez = 0.
                  xnew = xloc(k)
                  ynew = yloc(k)
                  do 130 j=1,len-1
                     j1 = mod(j,len )+1
                     ii1 = icyc(i,j+1)
                     ii2 = icyc(i,j1+1)
                     xi = x(icyc(i,j+1))
                     yi = y(icyc(i,j+1))
                     xi1 = x(icyc(i,j1+1))
                     yi1 = y(icyc(i,j1+1))
!     check the orientation and bad conditionality

                     det = POS_TEST(xi, yi, xi1, yi1, xnew, ynew )
                     !print *,'..MOVc.',det

                     if( det .le. 1.) goto 140

                     call POS1TEST(xi,yi,xi1,yi1,xnew,ynew,itet)
                     if(itet == 1) then
                        goto 140
                     endif

                     acc = ACCUTE(xi,yi,xi1,yi1,xnew,ynew,5.)

                     rl1 = (0.5*((rga(i)+rga(ii1))*(xi-xnew)*(xi-xnew)+  &
                          2*(rgb(i)+rgb(ii1))*(xi-xnew)*(yi-ynew)+  &
                          (rgc(i)+rgc(ii1))*(yi-ynew)*(yi-ynew)))**0.5

                     quality1rez = quality1rez + (rl1 - 1.732050807)**2
                     !     ... NEW ACCUTE
!                     quality1rez = quality1rez
!     *                    + (rl1 - 1.732050807)**2*acc

                     if(j == len -1) then
                        rl1 = (0.5*((rga(i)+rga(ii2))*  &
                             (xi1-xnew)*(xi1-xnew)+  &
                             2*(rgb(i)+rgb(ii2))*  &
                             (xi1-xnew)*(yi1-ynew)+  &
                             (rgc(i)+rgc(ii2))*  &
                             (yi1-ynew)*(yi1-ynew)))**0.5

                        quality1rez = quality1rez +   &
                             (rl1 - 1.732050807)**2
                     endif

                     ice = 0
130               enddo
                  quality1rez = quality1rez/len

                  if(quality1rez .lt. qualityrez ) then
                     qualityrez = quality1rez
                     jbet = k
                  endif
140            enddo
               if(jbet .gt. 0) then
                  x(i) = xloc(jbet)
                  y(i) = yloc(jbet)
                  ibb(i,1) = iloc(jbet)
               endif
            else
               !     for periodical points
               i1 = ibp(i,1)
               if(ibp(i,2) == 1) then
                  xperreal = AMA%xper(1,1)
                  yperreal = AMA%xper(1,2)
               else
                  xperreal = AMA%xper(2,1)
                  yperreal = AMA%xper(2,2)
               endif
               len1 = abs(icyc(i1,1))
               if(abs(x(i1) + xperreal - x(i) ) .lt. 1E-05 .and.  &
                    abs(y(i1) + yperreal - y(i) ) .lt. 1E-05 ) then
                  imov = 1
               elseif(abs(x(i1) - xperreal - x(i) ) .lt. 1E-05 .and.  &
                    abs(y(i1) - yperreal - y(i) ) .lt. 1E-05 ) then
                  imov = -1
               else
                  print *,'BAD in moving in periodical points'
                  print *,i,i1,imov
                  print *,x(i),y(i),xperreal
                  print *,x(i1),y(i1),yperreal
                  print *,abs(x(i1) + xperreal - x(i) ),  &
                       abs(y(i1) + yperreal - y(i) ),  &
                       abs(x(i1) - xperreal - x(i) ),  &
                       abs(y(i1) - yperreal - y(i) )
                  stop
               endif

               jbet = 0
               do 340 k=1,ilo
                  quality1 = 0.
                  quality1rez = 0.
                  xnew = xloc(k)
                  ynew = yloc(k)
                  do 330 j=1,len-1
                     j1 = mod(j,len )+1
                     ii1 = icyc(i,j+1)
                     ii2 = icyc(i,j1+1)
                     xi = x(icyc(i,j+1))
                     yi = y(icyc(i,j+1))
                     xi1 = x(icyc(i,j1+1))
                     yi1 = y(icyc(i,j1+1))

!     check the orientation and bad conditionality
                     det = POS_TEST(xi, yi, xi1, yi1, xnew, ynew )
                     !print *,'..MOVc.',det

                     if( det .le. 1.) goto 340

                     call POS1TEST(xi,yi,xi1,yi1,xnew,ynew,itet)
                     if(itet == 1) then
                        goto 340
                     endif

                     acc = ACCUTE(xi,yi,xi1,yi1,xnew,ynew, 1.)

                     rl1 = (0.5*((rga(i)+rga(ii1))*(xi-xnew)*(xi-xnew)+  &
                          2*(rgb(i)+rgb(ii1))*(xi-xnew)*(yi-ynew)+  &
                          (rgc(i)+rgc(ii1))*(yi-ynew)*(yi-ynew)))**0.5

                     quality1rez = quality1rez + (rl1 - 1.732050807)**2

!     ... NEW ACCUTE
!                     quality1rez = quality1rez
!     *                    +(rl1 - 1.732050807)**2*acc

                     if(j == len -1) then
                        rl1 = (0.5*((rga(i)+rga(ii2))*  &
                             (xi1-xnew)*(xi1-xnew)+  &
                             2*(rgb(i)+rgb(ii2))*  &
                             (xi1-xnew)*(yi1-ynew)+  &
                             (rgc(i)+rgc(ii2))*  &
                             (yi1-ynew)*(yi1-ynew)))**0.5

                        quality1rez = quality1rez +   &
                             (rl1 - 1.732050807)**2
                     endif

330               enddo
                  xnew = xloc(k) - imov*xperreal
                  ynew = yloc(k) - imov*yperreal
                  do 350 j=1,len1-1
                     j1 = mod(j,len1 )+1
                     ii1 = icyc(i1,j+1)
                     ii2 = icyc(i1,j1+1)
                     xi = x(icyc(i1,j+1))
                     yi = y(icyc(i1,j+1))
                     xi1 = x(icyc(i1,j1+1))
                     yi1 = y(icyc(i1,j1+1))

!     check the orientation and bad conditionality
                     det = POS_TEST(xi, yi, xi1, yi1, xnew, ynew )
                     !print *,'..MOVc.',det

                     if( det .le. 1.) goto 340

                     call POS1TEST(xi,yi,xi1,yi1,xnew,ynew,itet)
                     if(itet == 1) then
                        goto 340
                     endif

                     acc = ACCUTE(xi,yi,xi1,yi1,xnew,ynew, 1.)

                     rl1 = (0.5*((rga(i)+rga(ii1))*(xi-xnew)*(xi-xnew)+  &
                          2*(rgb(i)+rgb(ii1))*(xi-xnew)*(yi-ynew)+  &
                          (rgc(i)+rgc(ii1))*(yi-ynew)*(yi-ynew)))**0.5

                     quality1rez = quality1rez + (rl1 - 1.732050807)**2
!     ... NEW ACCUTE
!                     quality1rez = quality1rez
!     *                    +(rl1 - 1.732050807)**2*acc

                     if(j == len1 -1) then
                        rl1 = (0.5*((rga(i)+rga(ii2))*  &
                             (xi1-xnew)*(xi1-xnew)+  &
                             2*(rgb(i)+rgb(ii2))*  &
                             (xi1-xnew)*(yi1-ynew)+  &
                             (rgc(i)+rgc(ii2))*  &
                             (yi1-ynew)*(yi1-ynew)))**0.5

                        quality1rez = quality1rez +   &
                             (rl1 - 1.732050807)**2
                     endif
350               enddo
                  quality1rez = quality1rez/(len+len1+2)

                  if(quality1rez .lt. qualityrez ) then
                     qualityrez = quality1rez
                     jbet = k
                  endif
340            enddo
               if(jbet .gt. 0) then
                  x(i) = xloc(jbet)
                  y(i) = yloc(jbet)
                  x(i1) = xloc(jbet) - imov*xperreal
                  y(i1) = yloc(jbet) - imov*yperreal
               endif
            endif

         endif
21       continue
100   enddo

      deallocate(xloc, yloc, iloc )

      return
    end subroutine

    subroutine REM_BOUND(ndim, icha, icy)
      implicit none
      integer, intent(in) :: ndim, icy
      integer, intent(inout) :: icha
      integer, dimension(:), allocatable :: nsr, locyc
      integer :: ipoc, i, ikk1, ibo, ice, iyes, k, kik, k0, k1, k2, kk
      integer :: i0, i1, i2, j, j1, j2, kj1, kj2, kj3, kjj1, kjj2, kjj3
      integer :: ka1, ka2, ij, is, ib, ibstart1, ibstart2, itest, itet
      integer :: ka1j, ka2j, ib1, ie, iy
      real :: det, rep, err0, errrez0, err1, errrez1, err2, errrez2
      integer, dimension(:,:), pointer :: ibb, ibp
      real, dimension(:), pointer :: rga, rgb, rgc
      integer, dimension(:,:), pointer :: icyc
      real, dimension(:,:), pointer :: wp
      integer, dimension(:,:), pointer :: lnd, iae, lbn
      integer, dimension(:), pointer :: ibc, itc
      real, dimension(:), pointer :: x, y
      !integer :: ic_start, ic_end, ic_skip

      x => AMA%x(1:AMA%mpoin)
      y => AMA%y(1:AMA%mpoin)

      lbn => AMA%lbn(1:AMA%mbelm,1:2)
      ibc => AMA%ibc(1:AMA%mbelm)
      itc => AMA%itc(1:AMA%mbelm)


      lnd => AMA%lnd(1:AMA%melem, 1:3)
      iae => AMA%iae(1:AMA%melem, 1:3)

      wp    => AMA%wp(   1:AMA%mpoin,1:ndim+1)

      icyc => AMA%icyc(1:AMA%mpoin, 1:AMA%maxdeg)

      rga => AMA%rga( 1:AMA%mpoin )
      rgb => AMA%rgb( 1:AMA%mpoin )
      rgc => AMA%rgc( 1:AMA%mpoin )

      ibp => AMA%ibp(1:AMA%mpoin, 1:2)
      ibb => AMA%ibb(1:AMA%mpoin, 1:3)


      allocate(nsr(1:3), locyc(1:30) )

      icha = 0
      ipoc = 0
 10   ipoc = ipoc + 1
      if(ipoc .gt. AMA%npoin) return
      i = ipoc

      if(ibb(i,1) .gt.0 ) then
         ikk1 = ibb(i,2)
         if(ibb(i,1) == AMA%ibpoin(ikk1-1)+1  .or.  &
              ibb(i,1) == AMA%ibpoin(ikk1) .or. ibb(i,3) == -1) then
!     no removing, initial or final point of profiles
            goto 10
         endif
      endif
      if(icyc(i,1) == -3 ) then
         i0 = icyc(i,2)
         i1 = icyc(i,3)
         i2 = icyc(i,4)
         if(icyc(i1,1) .lt. 0) goto 101
         det = x(i0)*(y(i)-y(i2)) + x(i)*(y(i2)-y(i0)) +   &
              x(i2)*(y(i0)-y(i))

         rep =  ((x(i0)-x(i))**2 + (y(i0)-y(i))**2) +  &
              ((x(i)-x(i2))**2 + (y(i)-y(i2))**2) +  &
              ((x(i2)-x(i0))**2 + (y(i2)-y(i0))**2)

         if (abs(det)/rep .lt. 1e-02) then
!                  the points i0, i, i2 are in the straight,
!                  we can remove one triangle

            det = POS_TEST(x(i0), y(i0), x(i1), y(i1), x(i2), y(i2) )
            !print *,'..RBa.',det

            if( det .le. 1.) then
!               print *,'violation of positivity'
               goto 100
            endif
            call POS1TEST(x(i0),y(i0),x(i1),y(i1),x(i2),y(i2),itet)
            if(itet == 1) then
               goto 100
            endif

            ibo = 0
            call QUA2(x(i0),y(i0),wp(i0,1),rga(i0),rgb(i0),rgc(i0),  &
                 x(i1),y(i1),wp(i1,1),rga(i1),rgb(i1),rgc(i1),  &
                 x(i2),y(i2),wp(i2,1),rga(i2),rgb(i2),rgc(i2),  &
                 err0,errrez0,ice,ibo)
            call QUA2(x(i0),y(i0),wp(i0,1),rga(i0),rgb(i0),rgc(i0),  &
                 x(i1),y(i1),wp(i1,1),rga(i1),rgb(i1),rgc(i1),  &
                 x(i),y(i),wp(i,1),rga(i),rgb(i),rgc(i),  &
                 err2,errrez1,ice,ibo)
            call QUA2(x(i),y(i),wp(i,1),rga(i),rgb(i),rgc(i),  &
                 x(i1),y(i1),wp(i1,1),rga(i1),rgb(i1),rgc(i1),  &
                 x(i2),y(i2),wp(i2,1),rga(i2),rgb(i2),rgc(i2),  &
                 err1,errrez2,ice,ibo)

            iyes = 0
            if( 1.8*errrez0 .lt. errrez1 + errrez2) then
               iyes = 1
            endif

            if(ibp(i,1) .gt. 0 .and. iyes == 1) then
               k=ibp(i,1)
               kik = k
               if(icyc(k,1) .ne. -3) goto 100
               k0 = icyc(k,2)
               k1 = icyc(k,3)
               k2 = icyc(k,4)
               if(icyc(k1,1) .lt. 0) goto 100
               det = x(k0)*(y(k)-y(k2)) + x(k)*(y(k2)-y(k0)) +   &
                    x(k2)*(y(k0)-y(k))
               rep =  ((x(k0)-x(k))**2 + (y(k0)-y(k))**2) +  &
                    ((x(k)-x(k2))**2 + (y(k)-y(k2))**2) +  &
                    ((x(k2)-x(k0))**2 + (y(k2)-y(k0))**2)

               if (abs(det)/rep .lt. 1e-02) then
!     the points k0, k, k2 are in the straight,
!     we can remove one triangle
                  det = POS_TEST(x(k0), y(k0), x(k1), y(k1), x(k2), y(k2) )
                  !print *,'..RBb.',det

                  if( det .le. 1.) then
                     !     print *,'violation of positivity'
                     goto 100
                  endif
                  call POS1TEST(x(k0),y(k0),x(k1),y(k1),  &
                       x(k2),y(k2),itet)
                  if(itet == 1) then
                     goto 100
                  endif

                  ibo = 0
                  call QUA2(x(k0),y(k0),wp(k0,1),  &
                       rga(k0),rgb(k0),rgc(k0),  &
                       x(k1),y(k1),wp(k1,1),rga(k1),rgb(k1),rgc(k1),  &
                       x(k2),y(k2),wp(k2,1),rga(k2),rgb(k2),rgc(k2),  &
                       err0,errrez0,ice,ibo)
                  call QUA2(x(k0),y(k0),wp(k0,1),  &
                       rga(k0),rgb(k0),rgc(k0),  &
                       x(k1),y(k1),wp(k1,1),rga(k1),rgb(k1),rgc(k1),  &
                       x(k),y(k),wp(k,1),rga(k),rgb(k),rgc(k),  &
                       err2,errrez1,ice,ibo)
                  call QUA2(x(k),y(k),wp(k,1),rga(k),rgb(k),rgc(k),  &
                       x(k1),y(k1),wp(k1,1),rga(k1),rgb(k1),rgc(k1),  &
                       x(k2),y(k2),wp(k2,1),rga(k2),rgb(k2),rgc(k2),  &
                       err1,errrez2,ice,ibo)


                  if( 1.8*errrez0 .lt. errrez1 + errrez2 ) then
                     iyes = 2
                  else
                     iyes = 0
                  endif
               endif
            endif

            if( iyes .ge. 1) then
!     better quality, we  remove a boundary points
 999           i0 = icyc(i,2)
               i1 = icyc(i,3)
               i2 = icyc(i,4)
               itest = 0
               do 20 kk=1,AMA%nelem
                  do 30 j=1,3
                     j1 = mod(j,3) + 1
                     j2 = mod(j1,3) + 1
                     if(lnd(kk,j) == i0 .and. lnd(kk,j1) == i1 .and.  &
                          lnd(kk,j2) == i) then
                        k1 = kk
                        kj1 = j
                        kj2 = j1
                        kj3 = j2
                        itest = itest + 1
                     endif
!                     print *,'k1=',k1,kk,AMA%nelem,i0,i1,i
                     if(lnd(kk,j) == i2 .and. lnd(kk,j1) == i .and.  &
                          lnd(kk,j2) == i1) then
                        k2 = kk
                        kjj1 = j
                        kjj2 = j1
                        kjj3 = j2
                        itest = itest + 1
                     endif
                     if(itest == 2) goto 25
 30               enddo
 20            enddo
 25            continue

               ka1 = iae(k1,kj1)
               if(ka1 .gt. 0) then
                  do 40 j=1,3
                     if(iae(ka1,j) == k1) ka1j = j
 40               enddo
               endif

               ka2 = iae(k2,kjj3)
               if(ka2 .gt. 0) then
                  do 45 j=1,3
                     if(iae(ka2,j) == k2) ka2j = j
 45               enddo
               endif

               lnd(k1,kj3) = i2
               iae(k1,kj2) = ka2
               if(ka2 .gt. 0) iae(ka2,ka2j) = k1

               do ie = 1,AMA%nelem
                  do j=1,3
                     if(iae(ie,j) == AMA%nelem .and. k2 .ne. AMA%nelem)   &
                          iae(ie,j) = k2
                     if(lnd(ie,j) == AMA%npoin) lnd(ie,j) = i
                  enddo
               enddo

               do 70 ij =1,3
                  if(k2 .ne. AMA%nelem) then
                     lnd(k2,ij) = lnd(AMA%nelem,ij)
                     iae(k2,ij) = iae(AMA%nelem,ij)
                  endif
 70            enddo


               if (kik == AMA%npoin) kik = i

               x(i) = x(AMA%npoin)
               y(i) = y(AMA%npoin)
               ibb(i,1) = ibb(AMA%npoin,1)
               ibb(i,2) = ibb(AMA%npoin,2)
               ibb(i,3) = ibb(AMA%npoin,3)
               rga(i) = rga(AMA%npoin)
               rgb(i) = rgb(AMA%npoin)
               rgc(i) = rgc(AMA%npoin)
               wp(i,1) = wp(AMA%npoin,1)
               ibp(i,1) = ibp(AMA%npoin,1)
               ibp(i,2) = ibp(AMA%npoin,2)
               do iy=1,AMA%npoin
                  if(ibp(iy,1) == AMA%npoin) ibp(iy,1) = i
               enddo


               icyc(i0,abs(icyc(i0,1))+1) = i2
               icyc(i2,2) = i0
               do 200 k=1,icyc(i1,1)
                  if(icyc(i1,k+1) == i) then
                     kk = k
                     goto 210
                  endif
 200           enddo
 210           continue

               do 220 k=kk,icyc(i1,1)-1
                  icyc(i1,k+1) = icyc(i1,k+2)
 220           enddo
               icyc(i1,1) = icyc(i1,1) - 1

               do is = 1,AMA%npoin
                  do j=1,abs(icyc(is,1))
                     if(icyc(is,j+1) == AMA%npoin)icyc(is,j+1) =i
                  enddo
               enddo

               do 260 j=1,abs(icyc(AMA%npoin,1))+1
                  icyc(i,j)=icyc(AMA%npoin,j)
 260           enddo



!     reparation of lbn, ibc, itc
               do ib=1,AMA%nbelm
                  if(lbn(ib,2) == i) then
                     lbn(ib,2) = i0
                     itc(ib) = k1
                     ibstart1 = ib
                  elseif(lbn(ib,1) == i) then
                     ibstart2 = ib
                  endif
               enddo

!               print *,'ibstart =',ibstart1,ibstart2


               do 698 ib=1,AMA%nbelm
                  if(lbn(ib,1) == AMA%npoin) lbn(ib,1) = i
                  if(lbn(ib,2) == AMA%npoin) lbn(ib,2) = i
 698           enddo
               do 710 ib1 = ibstart2,AMA%nbelm-1
                  lbn(ib1,1) = lbn(ib1+1,1)
                  lbn(ib1,2) = lbn(ib1+1,2)
                  ibc(ib1) = ibc(ib1+1)
                  itc(ib1) = itc(ib1+1)
 710           enddo
               do ib1 = 1,AMA%nbelm-1
                  if(itc(ib1) == AMA%nelem) then
                     itc(ib1) = k2
                  endif
               enddo

               AMA%nbelm = AMA%nbelm - 1
               AMA%nelem = AMA%nelem - 1
               AMA%npoin = AMA%npoin - 1

               icha = icha + 1
               if(AMA%ifig .ge. 0 .and. mod(icha,5) == 0 ) then
                  call PLOT1()
                  AMA%ifig = AMA%ifig + 1
               endif


               if(iyes == 2) then
                  iyes = 3
                  i = kik
                  goto 999
               endif
            endif
 100        continue
         endif
101      continue
      endif
      if(ipoc .lt. AMA%npoin) goto 10

      deallocate( nsr, locyc )

    end subroutine REM_BOUND


    subroutine REMOVE(ndim, icha, icy)
      implicit none
      integer, intent(in) :: ndim, icy
      integer, intent(inout) :: icha
      integer, dimension(:), allocatable :: nsr, locyc
      integer :: i, ii, iia1, iia2, ikk1, ikk2, il, il0, ilen, ipoc, is
      integer :: j, j1, j2, jj2, itet, jja1, jja2, jjj, jminhelp
      integer :: k, k1, k2, k3, k4, k3len, kb0, kb1, kb2, ke1, ke2, ke3, kj
      integer :: kk, kk1, kk2, kl, kl1, kll, l, iemin, ig1, ig2, ikr1, ikr2
      integer :: ic1, icc, ie, iel, iemax, ib, ibbk1, ibbk11, ibbk2, ibbk111
      integer :: ilen1, ilen2, ilk, ilk1, ill, imov, ip, ip1, ip2, ipe
      integer :: ir1, ir2, irr, is1, is2, iskk2, it, itestit, itest2, iyi
      integer :: jb, jdo, ja1, ja2, ja21, je1, je2, jel, jj, jj1
      integer :: ia1, ia2
      real :: det, ss, ve1x, ve1y, ve2x, ve2y, rlmin2
      real :: rgai, rgbi, rgci, rgak, rgbk, rgck, rlen, rll, rminhelp, w0, we0
      real :: wpk1, x0, y0, x1, x2, y1, y2, xe0, ye0, xk1, yk1, xperreal, yperreal
      real :: xx0, xx1, xx2, yy0, yy1, yy2, a, b, c
      real rmin(3)
      integer jmin(3)
      integer:: nelemold ! local variable
      real, dimension(:), pointer :: xb, yb
      integer, dimension(:,:), pointer :: ibb, ibp
      real, dimension(:), pointer :: rga, rgb, rgc
      integer, dimension(:,:), pointer :: icyc
      real, dimension(:,:), pointer :: wp
      integer, dimension(:,:), pointer :: lnd, iae, lbn
      integer, dimension(:), pointer :: ibc, itc
      real, dimension(:), pointer :: x, y
      integer :: ic_start, ic_end, ic_skip

      x => AMA%x(1:AMA%mpoin)
      y => AMA%y(1:AMA%mpoin)

      lbn => AMA%lbn(1:AMA%mbelm,1:2)
      ibc => AMA%ibc(1:AMA%mbelm)
      itc => AMA%itc(1:AMA%mbelm)


      lnd => AMA%lnd(1:AMA%melem, 1:3)
      iae => AMA%iae(1:AMA%melem, 1:3)

      wp    => AMA%wp(   1:AMA%mpoin,1:ndim+1)
      icyc => AMA%icyc(1:AMA%mpoin, 1:AMA%maxdeg)

      rga => AMA%rga( 1:AMA%mpoin )
      rgb => AMA%rgb( 1:AMA%mpoin )
      rgc => AMA%rgc( 1:AMA%mpoin )

      ibp => AMA%ibp(1:AMA%mpoin, 1:2)
      ibb => AMA%ibb(1:AMA%mpoin, 1:3)

      xb => AMA%xb(1:AMA%ipoint)
      yb => AMA%yb(1:AMA%ipoint)

      allocate(nsr(1:3), locyc(1:50) )

      rlmin2 = 1.33

      icha = 0
      nelemold = AMA%nelem
      ipoc = 1


      if(mod(icy, 2) == 0) then
         ic_start = 1;         ic_end = nelemold;   ic_skip =  1
      else
         ic_start = nelemold;  ic_end = 1;          ic_skip = -1
      endif

      !do 20 iyi=1,nelemold
      do 20 iyi = ic_start, ic_end, ic_skip
         i = ipoc
         if( i .le. AMA%nelem) then
            do 27 jdo=1,3
               j = jdo
               j1 = mod(j,3) + 1
               k1 = lnd(i,j)
               k2 = lnd(i,j1)
               x1 = x(k1)
               y1 = y(k1)
               x2 = x(k2)
               y2 = y(k2)
               a = (rga(k1) + rga(k2) )/2
               b = (rgb(k1) + rgb(k2) )/2
               c = (rgc(k1) + rgc(k2) )/2
               rmin(j) = (a*(x1-x2)*(x1-x2) + c*(y1-y2)*(y1-y2)  &
                    + 2*b*(x1-x2)*(y1-y2))
               jmin(j) = j
 27         enddo

            do 25 k=1,3
               do 26 l=1,2
                  if(rmin(l) .gt. rmin(l+1) ) then
                     rminhelp = rmin(l)
                     rmin(l) = rmin(l+1)
                     rmin(l+1) = rminhelp
                     jminhelp = jmin(l)
                     jmin(l) = jmin(l+1)
                     jmin(l+1) = jminhelp
                  endif
 26            enddo
 25         enddo

            do 28 jdo=1,3
               j = jmin(jdo)
               j1 = mod(j,3) + 1
               k1 = lnd(i,j)
               k2 = lnd(i,j1)
               x1 = x(k1)
               y1 = y(k1)
               x2 = x(k2)
               y2 = y(k2)
               ig1 = 0
               ig2 = 0
               do 21 ib=1,AMA%nbelm
                  if(lbn(ib,1) == k1) ig1 = 1
                  if(lbn(ib,1) == k2) ig2 = 1
 21            enddo

               if(ig1 == 1 .and. icyc(k1,1) .ge. 0) then
                  print *,'ERROR in REMOVE - boundary point 1'
                  kl = k1
                  write(*,'(i5,2e12.4,2i5)') kl,x(kl),y(kl),  &
                       icyc(kl,1),ig1
                  do kl1 = 1,abs(icyc(kl,1))
                     kll = icyc(kl,kl1+1)
                     write(*,'(2e12.4,i5)') x(kll),y(kll),kll
                  enddo
                  stop
               endif
               if(ig2 == 1 .and. icyc(k2,1) .ge. 0) then
                  print *,'ERROR in REMOVE - boundary point 2'
                  kl = k2
                  write(*,'(i5,2e12.4,2i5)') kl,x(kl),y(kl),  &
                       icyc(kl,1),ig2
                  do kl1 = 1,abs(icyc(kl,1))
                     kll = icyc(kl,kl1+1)
                     write(*,'(2e12.4,i5)') x(kll),y(kll),kll
                  enddo

                  do ib=1,AMA%nbelm
                     if(lbn(ib,1) == kl) then
                        print *,'^^',lbn(ib,1),ib,x(lbn(ib,1)),  &
                             y(lbn(ib,1)),ibb(kl,1)
                     endif
                  enddo
                 stop
               endif

               if(iae(i,j) .gt. 0) then
!     for non boundary sides
                  if(ig1+ig2 .le. 1 .and. rmin(jdo) .le. rlmin2) then
!     we  remove this side (may be)

                     j2 = mod(j1,3)+1
                     ii = iae(i,j)
                     jj = 0
                     do 35 jjj=1,3
                        if(iae(ii,jjj) == i) jj = jjj
 35                  enddo
                     if(jj == 0) then
                        print *,'error 1342'
                        stop
                     endif
                     jj1 = mod(jj,3)+1
                     jj2 = mod(jj1,3)+1
                     k3 = lnd(i,j2)
                     k4 = lnd(ii,jj2)

                     if(k2.ne.lnd(ii,jj) .or. k1.ne.lnd(ii,jj1))then
                        print *,'ERRROR in REMOVE'
                        print *,i,k1,k2,k3,k4
                        print *,i,lnd(i,1),lnd(i,2),lnd(i,3)
                        print *,ii,lnd(ii,1),lnd(ii,2),lnd(ii,3)
                        print *,x(k1),y(k1)
                        print *,x(k2),y(k2)
                        print *,x(k3),y(k3)
                        print *,x(k4),y(k4)
                        stop
                     endif

                     if(iae(i,j1) .gt. 0) then
                        ia1 = iae(i,j1)
                        do 40 kk =1,3
                           if(iae(ia1,kk) == i) ja1 = kk
 40                     enddo
                     else
                        ia1 = -2
                     endif
                     if(iae(i,j2) .gt. 0) then
                        ia2 = iae(i,j2)
                        do 50 kk =1,3
                           if(iae(ia2,kk) == i) ja2 = kk
 50                     enddo
                     else
                        ia2 = -2
                     endif
                     if(iae(ii,jj1) .gt. 0) then
                        iia1 = iae(ii,jj1)
                        do 60 kk =1,3
                           if(iae(iia1,kk) == ii) jja1 = kk
 60                     enddo
                     else
                        iia1 = -2
                     endif
                     if(iae(ii,jj2) .gt. 0) then
                        iia2 = iae(ii,jj2)
                        do 70 kk =1,3
                           if(iae(iia2,kk) == ii) jja2 = kk
 70                     enddo
                     else
                        iia2 = -2
                     endif
!     here we check, so we will not obtain the triangle with two
!     boundary sides
                     if(ia2 .lt. 0)then
                        if(iae(ia1,1) .lt. 0 .or. iae(ia1,2) .lt. 0   &
                             .or.iae(ia1,3) .lt. 0 ) then
                           goto 28
                        endif
                     endif
                     if(ia1 .lt. 0)then
                        if(iae(ia2,1) .lt. 0 .or. iae(ia2,2) .lt. 0   &
                             .or.iae(ia2,3) .lt. 0 )then
                           goto 28
                        endif
                     endif

                     if(iia2 .lt. 0)then
                        if(iae(iia1,1) .lt. 0 .or. iae(iia1,2) .lt. 0   &
                             .or.iae(iia1,3) .lt. 0 ) then
                          goto 28
                       endif
                     endif
                     if(iia1 .lt. 0)then
                        if(iae(iia2,1) .lt. 0 .or. iae(iia2,2) .lt. 0   &
                             .or.iae(iia2,3) .lt. 0 ) then
                          goto 28
                       endif
                     endif

!     here we check that after emoving don't arise a corner triangle
                     if(iia1 == ia2 .and. ia2 .gt. 0 .and.  &
                          iia2 .lt. 0 .and. ia1 .lt. 0) then
                        goto 28
                     endif
                     if(ia1 == iia2 .and. iia2 .gt. 0 .and.  &
                          ia2 .lt. 0 .and. iia1 .lt. 0) then
                        goto 28
                     endif

                     if( ig1 == 1) then
!     k1 is on the boundary, no changes
                        xk1 = x(k1)
                        yk1 = y(k1)
                        wpk1 = wp(k1,1)
                        rgak = rga(k1)
                        rgbk = rgb(k1)
                        rgck = rgc(k1)
                     elseif( ig2 == 1) then
!     k2 on the boundary
                        xk1 = x(k2)
                        yk1 = y(k2)
                        wpk1 = wp(k2,1)
                        rgak = rga(k2)
                        rgbk = rgb(k2)
                        rgck = rgc(k2)
                     else
                        xk1 = (x(k1) + x(k2) )/2
                        yk1 = (y(k1) + y(k2) )/2
                        wpk1 = (wp(k1,1) + wp(k2,1) )/2
                        rgak = (rga(k1) + rga(k2) )/2
                        rgbk = (rgb(k1) + rgb(k2) )/2
                        rgck = (rgc(k1) + rgc(k2) )/2



                     endif

!     we check that if some new triangle will satisfy the positivity

                     xx0 = xk1
                     yy0 = yk1
!     for point k1
                     if(icyc(k1,1) .gt. 0) then
                        ilen = icyc(k1,1)
                     else
                        ilen = abs(icyc(k1,1)) - 1
                     endif
                     do 450 kk=1,ilen
                        kk1 = mod(kk,icyc(k1,1)) + 1
                        if( icyc(k1,kk+1) .ne. k2 .and.   &
                             icyc(k1,kk1+1) .ne. k2) then
                           xx1 = x(icyc(k1,kk+1))
                           yy1 = y(icyc(k1,kk+1))
                           xx2 = x(icyc(k1,kk1+1))
                           yy2 = y(icyc(k1,kk1+1))

                           det = POS_TEST(xx0, yy0, xx1, yy1, xx2, yy2 )
                           !print *,'..REMb.',det

                           if( det .le. 1.) then
!     violation of positivity, go to next j
                              goto 28
                           endif
                           call POS1TEST(xx0,yy0,xx1,yy1,xx2,yy2,itet)
                           if(itet == 1) then
                              goto 28
                           endif
                        endif
 450                 enddo
!     now for k2

                     if(icyc(k2,1) .gt. 0) then
                        ilen = icyc(k2,1)
                     else
                        ilen = abs(icyc(k2,1)) - 1
                     endif
                     do 460 kk=1,ilen
                        kk2 = mod(kk,icyc(k2,1)) + 1
                        if( icyc(k2,kk+1) .ne. k1 .and.   &
                             icyc(k2,kk2+1) .ne. k1) then
                           xx1 = x(icyc(k2,kk+1))
                           yy1 = y(icyc(k2,kk+1))
                           xx2 = x(icyc(k2,kk2+1))
                           yy2 = y(icyc(k2,kk2+1))

                           det = POS_TEST(xx0, yy0, xx1, yy1, xx2, yy2 )
                           !print *,'..REMc.',det

                           if( det .le. 1.) then

!     violation of positivity, go to next j
                              goto 28
                           endif
                           call POS1TEST(xx0,yy0,xx1,yy1,xx2,yy2,itet)
                           if(itet == 1) then
                              goto 28
                           endif


                        endif
 460                 enddo

                     kk1 = k1
                     kk2 = k2
                     if( ig2 == 1) then
                        kk1 = k2
                        kk2 = k1
                     endif

                     x(kk1) = xk1
                     y(kk1) = yk1

                     wp(kk1,1) = wpk1
                     rga(kk1) = rgak
                     rgb(kk1) = rgbk
                     rgc(kk1) = rgck

!     shiftting of all array after removing k2

                     do 6745 ip=1,AMA%npoin
                        if(ibp(ip,1) == kk2 .and.kk1 .ne. AMA%npoin)   &
                             ibp(ip,1) = kk1
                        if(ibp(ip,1) == AMA%npoin .and.kk2 .ne. AMA%npoin)   &
                             ibp(ip,1) = kk2
 6745                enddo


                     ibp(kk2,1) = ibp(AMA%npoin,1)
                     ibp(kk2,2) = ibp(AMA%npoin,2)
                     x(kk2) = x(AMA%npoin)
                     y(kk2) = y(AMA%npoin)
                     ibb(kk2,1) = ibb(AMA%npoin,1)
                     ibb(kk2,2) = ibb(AMA%npoin,2)
                     ibb(kk2,3) = ibb(AMA%npoin,3)
                     wp(kk2,1) = wp(AMA%npoin,1)
                     rga(kk2) = rga(AMA%npoin)
                     rgb(kk2) = rgb(AMA%npoin)
                     rgc(kk2) = rgc(AMA%npoin)

                     iemin = min (i,ii)
                     iemax = max (i,ii)

                     if(ia1 .gt. 0) iae(ia1,ja1) = ia2
                     if(ia2 .gt. 0) iae(ia2,ja2) = ia1
                     if(iia1 .gt. 0) iae(iia1,jja1) = iia2
                     if(iia2 .gt. 0) iae(iia2,jja2) = iia1

                     do 355 ie=1,AMA%nelem
                        do 360 kj=1,3
                           if(lnd(ie,kj) == kk2 .and. kk1 .ne. AMA%npoin)  &
                                lnd(ie,kj) = kk1
                           if(lnd(ie,kj) == AMA%npoin .and. kk2 .ne.AMA%npoin)  &
                                lnd(ie,kj) = kk2
 360                    enddo
                        if( iemax .lt. AMA%nelem-1 ) then
                           do kj=1,3
                              if(iae(ie,kj)==AMA%nelem-1) iae(ie,kj)=iemin
                              if(iae(ie,kj) == AMA%nelem) iae(ie,kj)=iemax
                           enddo
                        elseif(iemax == AMA%nelem-1 ) then
                           do kj=1,3
                              if(iae(ie,kj) == AMA%nelem) iae(ie,kj)=iemin
                           enddo
                        elseif(iemax == AMA%nelem .and.   &
                                iemin .lt. AMA%nelem-1) then
                           do kj=1,3
                              if(iae(ie,kj)==AMA%nelem-1) iae(ie,kj)=iemin
                           enddo
                        elseif(iemax == AMA%nelem .and.   &
                                iemin == AMA%nelem-1) then
!     no performance
                        else
                           print *,'LOGICAL error in REMOVE - 111&*('
                           stop
                        endif
 355                 enddo

                     if( iemax .lt. AMA%nelem-1 ) then
                        do kj=1,3
                           lnd(iemin,kj) = lnd(AMA%nelem-1,kj)
                           iae(iemin,kj) = iae(AMA%nelem-1,kj)
                           lnd(iemax,kj) = lnd(AMA%nelem,kj)
                           iae(iemax,kj) = iae(AMA%nelem,kj)
                        enddo
                     elseif(iemax == AMA%nelem-1 ) then
                        do kj=1,3
                           lnd(iemin,kj) = lnd(AMA%nelem,kj)
                           iae(iemin,kj) = iae(AMA%nelem,kj)
                        enddo
                     elseif(iemax == AMA%nelem .and.   &
                             iemin .lt. AMA%nelem-1) then
                        do kj=1,3
                           lnd(iemin,kj) = lnd(AMA%nelem-1,kj)
                           iae(iemin,kj) = iae(AMA%nelem-1,kj)
                        enddo
                     elseif(iemax == AMA%nelem .and.   &
                             iemin == AMA%nelem-1) then
!     no performance
                     else
                        print *,'LOGICAL error in REMOVE - &*('
                        stop
                     endif

                     do k=1,AMA%nbelm
                        if(itc(k) == i) then
                           if(ia1 .lt. 0 .and. ia2 .gt. 0) then
                              itc(k) = ia2
                           elseif(ia2 .lt. 0 .and. ia1 .gt.0) then
                              itc(k) = ia1
                           else
                              print *,'LOG. EROR in REM 145'
                              stop
                           endif
                        elseif(itc(k) == ii) then
                           if(iia1 .lt. 0 .and. iia2 .gt. 0) then
                              itc(k) = iia2
                           elseif(iia2 .lt. 0 .and. iia1 .gt. 0) then
                              itc(k) = iia1
                           else
                              print *,'LOG. EROR in REM 146'
                              stop
                           endif
                        endif
                     enddo

                     do 300 k=1,AMA%nbelm
                        do 310 l=1,2
                           if(lbn(k,l) == kk2 .and. kk1.ne.AMA%npoin)   &
                                lbn(k,l) = kk1
                           if(lbn(k,l) == AMA%npoin .and. kk2 .ne. AMA%npoin)  &
                                lbn(k,l) = kk2
 310                    enddo
                        if( iemax .lt. AMA%nelem-1 ) then
                           if(itc(k) == AMA%nelem-1)  itc(k)=iemin
                           if(itc(k) == AMA%nelem)  itc(k)=iemax
                        elseif(iemax == AMA%nelem-1 ) then
                           if(itc(k) == AMA%nelem)  itc(k)=iemin
                        elseif(iemax == AMA%nelem .and.   &
                                iemin .lt. AMA%nelem-1 ) then
                           if(itc(k) == AMA%nelem-1)  itc(k)=iemin
                        elseif(iemax == AMA%nelem .and.   &
                                iemin == AMA%nelem-1 ) then
!     no performance
                        else
                           print *,'ERRor in REMOVE 123'
                           print *,'LOGICAL error'
                           stop
                        endif
 300                 enddo
!     connection of two cykles
                     is = 0
                     do 700 ic=1,abs(icyc(kk2,1))
                        ic1 = mod(ic,abs(icyc(kk2,1)) ) +1
                        if( icyc(kk2,ic1+1) == kk1) is = ic
 700                 enddo
                     if(is == 0) then
                        print *,'ERROR in REMOVE in @#$%'
                        stop
                     endif
                     is1 = mod(is,abs(icyc(kk2,1)) ) +1
                     is2 = mod(is1,abs(icyc(kk2,1)) ) +1
                     ip = icyc(kk2,is+1)
                     ip1 = icyc(kk2,is1+1)
                     ip2 = icyc(kk2,is2+1)

                     if(ip1 .ne. kk1) then
                        print *,'ERROR in $%$%$%',kk1,kk2
                        print *,is,is1,is2
                        print *,ip,ip1,ip2
                        stop
                     endif

                     ir2 = 0
                     ir1 = 0
                     do 710 ic = 1,abs(icyc(kk1,1))
                        if(icyc(kk1,ic+1) == ip2) ir2 = ic
                        if(icyc(kk1,ic+1) == ip) ir1 = ic
 710                 enddo

                     if(ir1 == 0 .or. ir2 == 0) then
                        print *,'ERROR in REMOVE in $$$$'
                        stop
                     endif

                     locyc(1) = 0
                     irr = 0
                     do 720 ic =1,ir2
                        if(icyc(kk1,ic+1) .ne. kk2) then
                           irr = irr + 1
                           locyc(irr+1) = icyc(kk1,ic+1)
                           locyc(1) = locyc(1) + 1
                        endif
 720                 enddo

                     irr = locyc(1)
                     iskk2 = is2
 730                 ikk2 = mod(iskk2, abs(icyc(kk2,1) )) + 1
                     if( icyc(kk2,ikk2+1) .ne. ip) then
                        irr = irr + 1
                        locyc(1) = locyc(1) + 1
                        locyc(irr+1) = icyc(kk2,ikk2+1)
                        iskk2 = ikk2
                        goto 730
                     endif

                     if( ir1 .gt. 2) then
                        do 740 ic= ir1,abs(icyc(kk1,1) )
                           irr = irr + 1
                           locyc(irr+1) = icyc(kk1,ic+1)
                           locyc(1) = locyc(1) + 1
 740                    enddo
                     endif


                     if(locyc(1) + 1  .gt. AMA%maxdeg) then
                        print *,'Too long cykles in REMOVE'
                        print *,icyc(kk1,1),icyc(kk2,1),locyc(1), AMA%maxdeg
                        stop
                     endif
                     if(icyc(kk1,1) .gt. 0) then
                        icyc(kk1,1) = locyc(1)
                     else
                        icyc(kk1,1) = -locyc(1)
                     endif

                     do 750 ic=1,locyc(1)
                        icyc(kk1,ic+1) = locyc(ic+1)
 750                 enddo
!     end of connection of two cykles

                     do 500 ip =1,AMA%npoin
                        do 510 ic =1,abs(icyc(ip,1))
                           if(icyc(ip,ic+1) == kk2)icyc(ip,ic+1) = kk1
                           if(icyc(ip,ic+1) == AMA%npoin)   &
                                icyc(ip,ic+1) = kk2
 510                    enddo
 500                 enddo

                     do 530 ic=1,abs(icyc(AMA%npoin,1)) + 1
                        icyc(kk2,ic) = icyc(AMA%npoin,ic)
 530                 enddo

                     do 400 ip1 =1,2
                        if(ip1 == 1) then
                           ip = k3
                           if(k3 == AMA%npoin ) ip = kk2
                        endif
                        if(ip1 == 2) then
                           ip = k4
                           if(k4 == AMA%npoin ) ip = kk2
                        endif
                        do 410 ic =1,abs(icyc(ip,1))
                           ic1 = mod(ic,abs(icyc(ip,1)) ) + 1
                           if( icyc(ip,ic+1) == icyc(ip,ic1+1) ) then
                              do 420 icc= ic1,abs(icyc(ip,1))-1
                                 icyc(ip,icc+1) = icyc(ip,icc+2)
 420                          enddo
                              if(icyc(ip,1) .gt. 0) then
                                 icyc(ip,1) = icyc(ip,1) - 1
                              else
                                 icyc(ip,1) = icyc(ip,1) + 1
                              endif
                              goto 400
                           endif
 410                    enddo
 400                 enddo
!     end of shiftting

!     all is done, we can continue

                     AMA%npoin = AMA%npoin - 1
                     AMA%nelem = AMA%nelem - 2
                     icha = icha + 1

                     if(AMA%ifig .ge. 0 .and. mod(icha,5) == 0 ) then
                        call PLOT1()
                        AMA%ifig = AMA%ifig + 1
                     endif


                  endif
               else
!     for boundary sides

                  if( ig1 + ig2 .ne. 2) then
                     print *,'ERROR in REMOVE'
                     print *,'the boundary segment do not have ',  &
                          'the points on the boundary'
                     print *,i,k1,k2,ig1,ig2,icha
                     print *,lnd(i,1),lnd(i,2),lnd(i,3)
                     print *,AMA%nelem,AMA%npoin,AMA%nbelm
                     print *,AMA%melem,AMA%mpoin,AMA%mbelm
                     print *,x1,y1
                     print *,x2,y2
                     stop
                  endif

                  itest2 = 0
                  if(rmin(jdo) .le. rlmin2) then
!     we  remove this edge
                     itest2 = 1
                     ipe = 0
 999                 j1 = mod(j,3) +1
                     j2 = mod(j1,3)+1
                     k1 = lnd(i,j)
                     k2 = lnd(i,j1)
                     k3 = lnd(i,j2)
                     x1 = x(k1)
                     y1 = y(k1)
                     x2 = x(k2)
                     y2 = y(k2)
                     kb0 = 0
                     kb1 = 0
                     kb2 = 0


                     do 610 ib=1,AMA%nbelm
                        if(lbn(ib,1) == k1 .and.   &
                             lbn(ib,2) == k2) then
                           kb1 = ib
                        elseif(lbn(ib,2) == k1) then
                           kb0 = ib
                        elseif(lbn(ib,1) == k2) then
                           kb2 = ib
                        endif
 610                 enddo

                     if( kb0*kb1*kb2 == 0 ) then
                        print *,'The bounary segment does not found'
                        print *,i,k1,k2
                        print *,kb0,kb1,kb2
                        print *,x1,y1,x(k1),y(k1)
                        print *,x2,y2,x(k2),y(k2)
                        stop
                     endif

                     if(ipe == 1) goto 993

!     "sharp angle
                     ikr1 = 0
                     ikr2 = 0
                     ibbk1 = ibb(k1,1)
                     ibbk11 = ibb(k1,2)
                     ibbk2 = ibb(k2,1)
                     ilk = 0
                     ilk1 = 0

                     if(ibb(k1,1) .gt. 0) then
                        ikk1 = ibb(k1,2)
                        if(ibb(k1,1) == AMA%ibpoin(ikk1-1)+1  .or.  &
                             ibb(k1,1) == AMA%ibpoin(ikk1) .or.   &
                             ibb(k1,3) == -1) then
!               NO moving of final or initial node of profiles
                           ikr1 = 1
                           ilk = ibb(k1,2)
                           ilk1 = ibb(k1,3)
                        endif
                     endif

                     if(ibb(k2,1) .gt. 0) then
                        ikk1 = ibb(k2,2)
                        if(ibb(k2,1) == AMA%ibpoin(ikk1-1)+1  .or.  &
                             ibb(k2,1) == AMA%ibpoin(ikk1) .or.   &
                             ibb(k2,3) == -1) then
!               NO moving of final or initial node of profiles
                           ikr2 = 1
                           ilk = ibb(k2,2)
                           ilk1 = ibb(k2,3)
                        endif
                     endif

                     if(ibbk1 .gt. 0 .and.  ibbk2 .gt. 0) then
                        ilk = ibb(k1,2)
                        ilk1 = ibb(k1,3)
                        if(ilk .ne. ibb(k2,2) ) then
                           print *,'error jkoli3',ss
                           print *,x(k1),y(k1),ibb(k1,1),ibb(k1,2)
                           print *,x(k2),y(k2),ibb(k2,1),ibb(k2,2)
                           stop
                        endif
                     endif

                     if(ilk .gt. 0 ) then
                        if(  ibbk1 ==  AMA%ibpoin(ilk-1) + 1 .or.   &
                             ibbk1 == AMA%ibpoin(ilk)) then
!     the point k1 we can not remove or move
                           ikr1 = 1
                        endif
                        if(  ibbk2 ==  AMA%ibpoin(ilk-1) + 1 .or.   &
                             ibbk2 == AMA%ibpoin(ilk)) then
!     the point k2 we can not remove or move
                           ikr2 = 1
                        endif
                     endif


                     rlen = ((x2-x(lbn(kb0,1)))**2 +   &
                          (y2-y(lbn(kb0,1)))**2)
                     det = x(lbn(kb0,1))*(y1-y2) + x1*(y2-y(lbn(kb0,1)))  &
                          + x2*(y(lbn(kb0,1)) - y1)

                     ve1x = x(lbn(kb0,1)) - x1
                     ve1y = y(lbn(kb0,1)) - y1
                     ve2x = x2 - x1
                     ve2y = y2 - y1

                     ss = (ve1x*ve2x + ve1y*ve2y)/  &
                          (ve1x**2 + ve1y**2)**0.5/  &
                          (ve2x**2 + ve2y**2)**0.5

                     if(ss .gt. -0.96 ) then
!     the point k1 we can not remove or move
                        ikr1 = 1
                     endif
                     rlen = ((x1-x(lbn(kb2,2)))**2 +   &
                          (y1-y(lbn(kb2,2)))**2)
                     det = x(lbn(kb2,2))*(y1-y2) + x1*(y2-y(lbn(kb2,2)))  &
                          + x2*(y(lbn(kb2,2))-y1)
                     ve1x = x(lbn(kb2,2)) - x2
                     ve1y = y(lbn(kb2,2)) - y2
                     ve2x = x1 - x2
                     ve2y = y1 - y2
                     ss = (ve1x*ve2x + ve1y*ve2y)/  &
                          (ve1x**2 + ve1y**2)**0.5/  &
                          (ve2x**2 + ve2y**2)**0.5

                     if(ss .gt. -0.96 ) then
!     the point k2 we can not remove or move
                        ikr2 = 1
                     endif


                     if( ikr2 == 0 ) then
                        if(ikr1 == 0) then
                           x0 = (x(k1) + x(k2) )/2
                           y0 = (y(k1) + y(k2) )/2
                           w0 = (wp(k1,1) + wp(k2,1) )/2
                           rgai = (rga(k1) + rga(k2)) /2
                           rgbi = (rgb(k1) + rgb(k2)) /2
                           rgci = (rgc(k1) + rgc(k2)) /2
                           if(ibb(k1,1).gt.0 .and. ibb(k2,1).gt.0) then
                              if(ibb(k2,1) .gt. ibb(k1,1) ) then
                                 if(ibb(k2,1) - ibb(k1,1) .ge. 2) then
                                    il0 = int(ibb(k1,1) + ibb(k2,1))/2
                                    x0 = xb(il0)
                                    y0 = yb(il0)
                                    ibbk1 = il0
                                    ibbk111 = ilk1
                                 else
!     a few points
                                    x0 = x(k1)
                                    y0 = y(k1)
                                 endif
                              else

                                 if(abs(ibb(k1,1)-ibb(k2,1)-  &
                                      AMA%ibpoin(ilk)+AMA%ibpoin(ilk-1)+1)  &
                                      .ge. 2) then
                                    il0=int(ibb(k1,1)+ibb(k2,1)+  &
                                       AMA%ibpoin(ilk)-AMA%ibpoin(ilk-1)+1)/2
                                    if(il0 .ge. AMA%ibpoin(ilk))   &
                                         il0=il0-AMA%ibpoin(ilk)+  &
                                         AMA%ibpoin(ilk-1)+1
                                    x0 = xb(il0)
                                    y0 = yb(il0)
                                    ibbk1 = il0
                                    ibbk11 = ilk
                                    ibbk111 = ilk1
                                 else
!     a few points
                                    x0 = x(k1)
                                    y0 = y(k1)
                                 endif
                              endif
                           endif
                        else
                           x0 = x(k1)
                           y0 = y(k1)
                           w0 = wp(k1,1)
                           rgai = rga(k1)
                           rgbi = rgb(k1)
                           rgci = rgc(k1)
                        endif
                     else
                        if(ikr1 == 0) then
                           x0 = x(k2)
                           y0 = y(k2)
                           w0 = wp(k2,1)
                           rgai = rga(k2)
                           rgbi = rgb(k2)
                           rgci = rgc(k2)
                           ibbk1 = ibb(k2,1)
                           ibbk11 = ibb(k2,2)
                           ibbk111 = ibb(k2,3)
                        else
!     no removing
                           goto 28
                        endif
                     endif
                     rll = ((x(k1) - x(k2))*(x(k1) - x(k2)) +  &
                          (y(k1) - y(k2))*(y(k1) - y(k2)) )

!     test the violation of positivity

                     ilen1 = abs(icyc(k1,1)) - 1
                     do 801 ic = 2,ilen1
                        ic1 = ic + 1
                        x1 = x(icyc(k1,ic+1))
                        y1 = y(icyc(k1,ic+1))
                        x2 = x(icyc(k1,ic1+1))
                        y2 = y(icyc(k1,ic1+1))

                        det = POS_TEST(x1, y1, x2, y2, x0, y0 )
                        !print *,'..REMc.',det

                        if( det .le. 1.) then
                           goto 28
                        endif
                        call POS1TEST(x0,y0,x1,y1,x2,y2,itet)
                        if(itet == 1) then
                           goto 28
                        endif


 801                 enddo

                     ilen2 = abs(icyc(k2,1)) - 1
                     do 802 ic = 1,ilen2-1
                        ic1 = ic + 1
                        x1 = x(icyc(k2,ic+1))
                        y1 = y(icyc(k2,ic+1))
                        x2 = x(icyc(k2,ic1+1))
                        y2 = y(icyc(k2,ic1+1))

                        det = POS_TEST(x1, y1, x2, y2, x0, y0 )
                        !print *,'..REMd.',det

                        if( det .le. 1.) then
!     violation of positivity
                           goto 28
                        endif
                        call POS1TEST(x0,y0,x1,y1,x2,y2,itet)
                        if(itet == 1) then
                           goto 28
                        endif
 802                 enddo

!     periodic boundary
                     if(ibp(k1,1) .gt. 0 .and. ibp(k2,1) .gt. 0) then
                        if(ibp(k1,2) .ne. ibp(k2,2) ) goto 28
                        if(ibp(k2,2) == 1) then
                           xperreal = AMA%xper(1,1)
                           yperreal = AMA%xper(1,2)
                        else
                           xperreal = AMA%xper(2,1)
                           yperreal = AMA%xper(2,2)
                        endif

!                        print *,'@@',k1,ibp(k1,1),ibp(k1,2),x(k1),y(k1)
!                        print *,'@@',k2,ibp(k2,1),ibp(k2,2),x(k2),y(k2)

                        ipe = -1
                        do 1000 iel =1,AMA%nelem
                           do 1001 jel=1,3
                              if(iae(iel,jel) .lt. 0 .and.  &
                                   lnd(iel,jel) == ibp(k2,1))goto 1002
 1001                      enddo
 1000                   enddo
 1002                   continue
                        je1 = mod(jel,3)+1
                        je2 = mod(je1,3)+1
                        ke1 = lnd(iel,jel)
                        ke2 = lnd(iel,je1)
                        ke3 = lnd(iel,je2)
                        if(abs(x(k2) + xperreal - x(ke1) ) .lt.   &
                             1E-05 .and.  &
                             abs(y(k2)+yperreal-y(ke1) ) .lt.   &
                             1E-05 ) then
                           imov = 1
                        elseif(abs(x(k2)-xperreal-x(ke1)) .lt.   &
                                1E-05 .and.  &
                                abs(y(k2)-yperreal-y(ke1)).lt.   &
                                1E-05 ) then
                           imov = -1
                        else
                           print *,'BAD in remove in periodical points'
                           print *,AMA%xper(1,1),AMA%xper(1,2),AMA%xper(2,1),AMA%xper(2,2)
                           print *,i,k2,ke1,ibp(k2,1),ibp(k2,2)
                           print *,x(k2),y(k2),xperreal
                           print *,x(ke1),y(ke1),yperreal
                           print *,abs(x(k2) + xperreal - x(ke1) ),  &
                                abs(y(k2) + yperreal - y(ke1) ),  &
                                abs(x(k2) - xperreal - x(ke1) ),  &
                                abs(y(k2) - yperreal - y(ke1) )
                           stop
                        endif
                        xe0 = x0 + imov*xperreal
                        ye0 = y0 + imov*yperreal
                        we0 = w0


!     test the violation of positivity
                        ilen1 = abs(icyc(ke1,1)) - 1
                        do 811 ic = 2,ilen1
                           ic1 = ic + 1
                           x1 = x(icyc(ke1,ic+1))
                           y1 = y(icyc(ke1,ic+1))
                           x2 = x(icyc(ke1,ic1+1))
                           y2 = y(icyc(ke1,ic1+1))

                           det = POS_TEST(x1, y1, x2, y2, xe0, ye0 )
                           !print *,'..REMd.',det

                           if( det .le. 1.) then
                              goto 28
                           endif
                           call POS1TEST(xe0,ye0,x1,y1,x2,y2,itet)
                           if(itet == 1) then
                              goto 28
                           endif
 811                    enddo

                        ilen2 = abs(icyc(ke2,1)) - 1
                        do 812 ic = 1,ilen2-1
                           ic1 = ic + 1
                           x1 = x(icyc(ke2,ic+1))
                           y1 = y(icyc(ke2,ic+1))
                           x2 = x(icyc(ke2,ic1+1))
                           y2 = y(icyc(ke2,ic1+1))

                           det = POS_TEST(x1, y1, x2, y2, xe0, ye0 )
                           !print *,'..REMe.',det

                           if( det .le. 1.) then
                              goto 28
                           endif
                           call POS1TEST(xe0,ye0,x1,y1,x2,y2,itet)
                           if(itet == 1) then
                              goto 28
                           endif
 812                    enddo
                     endif
 993                 continue

                     ia1 = iae(i,j2)
                     if(ia1 .gt. 0) then
                        do 611 il =1,3
                           if(iae(ia1,il) == i) ja1 = il
 611                    enddo
                     endif
                     ia2 = iae(i,j1)
                     if(ia2 .gt. 0) then
                        do 612 il =1,3
                           if(iae(ia2,il) == i) ja2 = il
 612                    enddo
                     endif
                     if(ia1 .lt. 0 .or. ia2 .lt. 0 ) then
!                        print *,'PRoblem in remove'
!                        print *,icha
                        goto 28
                     endif

                     iae(ia1,ja1) = ia2

                     ja21 = mod(ja2,3) + 1
                     lnd(ia2,ja21) = k1
                     iae(ia2,ja2) = ia1

!     begin of shifting 2
                     ibb(k1,1) = ibbk1
                     ibb(k1,2) = ibbk11
                     ibb(k1,3) = ibbk111

                     if(ipe == 0 .or. ipe == -1) then
                        x(k1) = x0
                        y(k1) = y0
                        wp(k1,1) = w0
                     elseif(ipe == 1) then
                        x(k1) = xe0
                        y(k1) = ye0
                        wp(k1,1) = we0
                     else
                        print *,'bad number ipe=',ipe
                        stop
                     endif

                     rga(k1) = rgai
                     rgb(k1) = rgbi
                     rgc(k1) = rgci


                     if(ibp(k2,1) .gt. 0 .and. ibp(k1,1) .le. 0) then
!     in fact we remove k1 and k2 stay with new index k1
                        if(ikr1 ==1 .or. ikr2 == 0) then
                           print *,'MISHMATCH in REMOVE'
                           print *,ikr1,ikr2
                           stop
                        endif
                        ibp(k1,1) = ibp(k2,1)
                        ibp(k1,2) = ibp(k2,2)
                     endif
                     ibp(k2,1) = ibp(AMA%npoin,1)
                     ibp(k2,2) = ibp(AMA%npoin,2)
                     do 6645 ip=1,AMA%npoin
                        if(ibp(ip,1) == k2 .and.k1 .ne. AMA%npoin)then
                           ibp(ip,1) = k1
                        endif
                        if(ibp(ip,1) == AMA%npoin .and.k2 .ne. AMA%npoin)   &
                             ibp(ip,1) = k2
 6645                enddo


                     x(k2) = x(AMA%npoin)
                     y(k2) = y(AMA%npoin)
                     ibb(k2,1) = ibb(AMA%npoin,1)
                     ibb(k2,2) = ibb(AMA%npoin,2)
                     ibb(k2,3) = ibb(AMA%npoin,3)
                     wp(k2,1) = wp(AMA%npoin,1)
                     rga(k2) = rga(AMA%npoin)
                     rgb(k2) = rgb(AMA%npoin)
                     rgc(k2) = rgc(AMA%npoin)

                     do 620 ie =1,AMA%nelem
                        do 630 kj =1,3
                           if(lnd(ie,kj) == k2 .and. k1 .ne. AMA%npoin)  &
                                lnd(ie,kj) = k1
                           if(lnd(ie,kj) == AMA%npoin .and. k2 .ne.AMA%npoin)  &
                                lnd(ie,kj) = k2
                           if(iae(ie,kj) == AMA%nelem .and. i .ne. AMA%nelem)  &
                                iae(ie,kj)=i
 630                    enddo
 620                 enddo
                     do kj=1,3
                        lnd(i,kj) = lnd(AMA%nelem,kj)
                        iae(i,kj) = iae(AMA%nelem,kj)
                     enddo


!     connection of two cykles 2
                     locyc(1) = icyc(k2,1) + 1
                     ill = abs(locyc(1))
                     if(abs(locyc(1)) .gt. AMA%maxdeg) then
                        print *,'2 - ERROR too long icyc'
                        stop
                     endif
                     do 705 ic=1,ill
                        locyc(ic+1) = icyc(k2,ic+1)
 705                 enddo
                     do 715 ic=1,abs(icyc(k1,1)) - 2
                        locyc(ill+ic+1) = icyc(k1,ic+3)
 715                 enddo
                     locyc(1) = locyc(1) - (abs(icyc(k1,1)) - 2)

                     if(abs(locyc(1)) .gt. AMA%maxdeg) then
                        print *,'3 - ERROR too long icyc'
                        stop
                     endif
                     do 723 ic=1,abs(locyc(1)) + 1
                        icyc(k1,ic) = locyc(ic)
 723                 enddo

                     itestit = 0
                     k3len = abs(icyc(k3,1))

                     if(icyc(k3,k3len+1) == k2) then
                        icyc(k3,1) = icyc(k3,1) - 1
                        itestit = 1
                        goto 729
                     endif

                     do 725 ic=1,k3len - 1
                        if(icyc(k3,ic+1) == k2) then
                           do 727 ic1=ic,k3len - 1
                              icyc(k3,ic1+1) = icyc(k3,ic1+2)
 727                       enddo
                           if(icyc(k3,1) .gt. 0) then
                              icyc(k3,1) = icyc(k3,1) - 1
                           else
                              icyc(k3,1) = icyc(k3,1) + 1
                           endif
                           itestit = 1
                           goto 729
                        endif
 725                 enddo
 729                 continue

                     if(itestit == 0) then
                        print *,'ERROR k1 does not found in icyc of k3'
                        print *,icha,k1,k2,k3
                        print *,x(k1),y(k1)
                        print *,x(k2),y(k2)
                        print *,x(k3),y(k3)
                        it = k3
                        write(*,'(10i5)') it,icyc(it,1),icyc(it,2),  &
                             icyc(it,3),icyc(it,4),icyc(it,5),  &
                             icyc(it,6),icyc(it,7),  &
                             icyc(it,8),icyc(it,9)

                        stop
                     endif
!     end of connection of two cykles 2

                     do 735 ie = 1,AMA%npoin
                        do 737 ic=1,abs(icyc(ie,1))
                           if(icyc(ie,ic+1) == k2) icyc(ie,ic+1) = k1
                           if(icyc(ie,ic+1) == AMA%npoin)   &
                                icyc(ie,ic+1) = k2
 737                    enddo
 735                 enddo

                     do  ic=1,abs(icyc(AMA%npoin,1)) + 1
                        icyc(k2,ic) = icyc(AMA%npoin,ic)
                     enddo

                     if(k1 .ne. AMA%npoin) then
                        lbn(kb2,1) = k1
                     else
                        lbn(kb2,1) = k2
                        lbn(kb0,2) = k2
                     endif
                     do 810 ib=1,AMA%nbelm
                        do 815 jb=1,2
                           if(lbn(ib,jb) == AMA%npoin .and. k2.ne.AMA%npoin)  &
                                lbn(ib,jb) = k2
 815                    enddo
                        if(itc(ib) == AMA%nelem) itc(ib) = i
 810                 enddo
                     do 820 ib=kb1,AMA%nbelm-1
                        lbn(ib,1) = lbn(ib+1,1)
                        lbn(ib,2) = lbn(ib+1,2)
                        ibc(ib) = ibc(ib+1)
                        itc(ib) = itc(ib+1)
 820                 enddo

!     end of shifting 2
                     AMA%npoin = AMA%npoin - 1
                     AMA%nelem = AMA%nelem - 1
                     AMA%nbelm = AMA%nbelm - 1

                     if( ipe == -1) then
                        if(ke1 == AMA%npoin+1) ke1 = k2
                        if(ke2 == AMA%npoin+1) ke2 = k2
                        do 1110 iel =1,AMA%nelem
                           do 1111 jel=1,3
                              if(iae(iel,jel) .lt. 0 .and.  &
                                   lnd(iel,jel) == ke1)goto 1112
 1111                      enddo
 1110                   enddo
 1112                   continue
                        i = iel
                        j = jel
                        ipe = 1
                        goto 999
                     endif
                     if(ipe == 1) then
                        ipe = 0
                        i = ipoc
                     endif
                     icha = icha + 1

                     if(AMA%ifig .ge. 0 .and. mod(icha,5) == 0 ) then
                        call PLOT1( )
                        AMA%ifig = AMA%ifig + 1
                     endif

                     do i=1,AMA%npoin
                        if(ibb(i,1) .gt. 0) then
                           do j=i+1,AMA%npoin
                              if(ibb(i,1) == ibb(j,1)) then
                                 print *,'the same IBB,  icha =',icha
                                 print *,xb(ibb(i,1)),yb(ibb(i,1)),i,  &
                                      ibb(i,1)
                                 print *,xb(ibb(j,1)),yb(ibb(j,1)),j,  &
                                      ibb(j,1)
                              endif
                           enddo
                        endif
                     enddo

                  endif
               endif


               if( i .gt. AMA%nelem) goto 2000
28          enddo
            ipoc = ipoc + 1
         endif
20    enddo
2000  continue

      deallocate( nsr, locyc )

    end subroutine REMOVE


    subroutine DELAUNAY(icha)
      implicit none
      integer, intent(inout) ::  icha
      real :: x1, y1, x2, y2, x3, y3, x4, y4, det123, det134
      real :: detdel, z1, z2, z3, z4
      integer :: i, j, j1, j2, ii, jj, jjj, jj1, j0, k1, k2, k3, k4
      integer :: itci, itcii, iyii, iy, ib1, ice, itest, itet, iyi, iyiij, iyij
      integer, dimension(:,:), pointer :: iba
      integer, dimension(:,:), pointer :: lnd, iae, lbn
      integer, dimension(:), pointer :: ibc, itc
      real, dimension(:), pointer :: x, y

      x => AMA%x(1:AMA%mpoin)
      y => AMA%y(1:AMA%mpoin)

      lbn => AMA%lbn(1:AMA%mbelm,1:2)
      ibc => AMA%ibc(1:AMA%mbelm)
      itc => AMA%itc(1:AMA%mbelm)


      lnd => AMA%lnd(1:AMA%melem, 1:3)
      iae => AMA%iae(1:AMA%melem, 1:3)

      iba => AMA%iba(1:AMA%melem, 1:3)

      icha = 0
!      return

      do i=1,AMA%nelem
         do j=1,3
            iba(i,j) = 0
         enddo
      enddo

      itest = -1499

      ice = 0
      do 10 i=1,AMA%nelem
         do 20 j=1,3
            if(iae(i,j) .gt. 0 .and. iba(i,j) == 0) then
               j1 = mod(j,3) +1
               j2 = mod(j1,3) +1
               ii = iae(i,j)

               do 30 jjj=1,3
                  if(lnd(ii,jjj) == lnd(i,j1)) jj = jjj
 30            enddo
               iba(i,j) = 1
               iba(ii,jj) = 1
               jj1 = mod(jj,3) +1
               j0 = mod(jj1,3) +1
               k1 = lnd(i,j2)
               k2 = lnd(i,j)
               k3 = lnd(ii,j0)
               k4 = lnd(i,j1)
               if( k2 .ne. lnd(ii,jj1) .or. k4 .ne. lnd(ii,jj) ) then
                  print *,'ERROR !@#'
                  print *,k1,k2,k3,k4
                  return
               endif
               x1 = x(k1)
               y1 = y(k1)
               x2 = x(k2)
               y2 = y(k2)
               x3 = x(k3)
               y3 = y(k3)
               x4 = x(k4)
               y4 = y(k4)

               if( i == itest) then
!                  print *,x(lnd(i,j)),y(lnd(i,j)),j,ii,iba(i,j)
!                  print *,x(lnd(i,j1)),y(lnd(i,j1))
!                  print *,x(lnd(i,j2)),y(lnd(i,j2))
                  print *
                  print *,x1,y1
                  print *,x2,y2
                  print *,x3,y3
                  print *,x4,y4
               endif

!     we prohibid SWAPPING in case, where can appear alement with two
!     boundary segment
               if( (iae(i,j2) .lt. 0 .and. iae(ii,jj1) .lt. 0 ) .or.  &
                    (iae(i,j1) .lt. 0 .and. iae(ii,j0) .lt. 0) ) then
                  goto 20
               endif
!     we must still chech the orientation, i.e. the  cykles of points
!     k1,k2,k4 and k1,k4,k3 must have positive orientation

               !det123 = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
               !det134 = x1*(y3-y4) + x3*(y4-y1) + x4*(y1-y3)
               !
               !reps123 = AMA%pos*( ((x1-x3)**2 + (y1-y3)**2) +  &
               !     ((x2-x3)**2 + (y2-y3)**2) +  &
               !     ((x1-x2)**2 + (y1-y2)**2) )
               !reps134 = AMA%pos*( ((x1-x3)**2 + (y1-y3)**2) +  &
               !     ((x4-x3)**2 + (y4-y3)**2) +  &
               !     ((x1-x4)**2 + (y1-y4)**2) )
               !if( det123 .le. reps123 .or. det134 .le. reps134) then


               det123 = POS_TEST(x1, y1, x2, y2, x3, y3 )
               det134 = POS_TEST(x1, y1, x3, y3, x4, y4 )
               !print *,'..SWA.',det

               if( det123 .le. 1. .or. det134 .le. 1.) then

!                  print *,'violation of positivity-1'
                  goto 20
               endif
               call POS1TEST(x(k1),y(k1),x(k2),y(k2),x(k3),y(k3),itet)
               if(itet == 1) then
!                  print *,'violation of positivity-2'
                  goto 20
               endif
               call POS1TEST(x(k1),y(k1),x(k3),y(k3),x(k4),y(k4),itet)
               if(itet == 1) then
!                  print *,'violation of positivity-3'
                  goto 20
               endif


               z1 = x1*x1 + y1*y1
               z2 = x2*x2 + y2*y2
               z3 = x3*x3 + y3*y3
               z4 = x4*x4 + y4*y4

               detdel = (x1*(y2*z3-y3*z2) - x2*(y1*z3-y3*z1) +  &
                    x3*(y1*z2-y2*z1))   &
                    -  (x1*(y2*z4-y4*z2) - x2*(y1*z4-y4*z1) +  &
                    x4*(y1*z2-y2*z1))   &
                    +  (x1*(y3*z4-y4*z3) - x3*(y1*z4-y4*z1) +  &
                    x4*(y1*z3-y3*z1))   &
                    -  (x2*(y3*z4-y4*z3) - x3*(y2*z4-y4*z2) +  &
                    x4*(y2*z3-y3*z2))

               if( i == itest) print *, detdel

!               if(detdel .lt. -1E-05) then
               if(detdel .lt. -1.D-12) then


!                  print *,x1,y1,detdel
!                  print *,x2,y2
!                  print *,x3,y3
!                  print *,x4,y4
!                  call PLOT( )
!
!                  stop


!     we swap the diagonal
                  itci = 0
                  itcii = 0

                  iyii = iae(ii,jj1)
                  if(iyii .gt. 0) then
                     do 62 iy=1,3
                        if( iae(iyii,iy) == ii) then
                           iyiij = iy
                        endif
 62                  enddo
                  else
!     boundary segment, we seek itc
                     do ib1 =1,AMA%nbelm
                        if(itc(ib1) == ii) then
                           itci = ib1
                           goto 662
                        endif
                     enddo
                     print *,'boundary segment in DELANAY for itc not'
                     stop
                  endif
 662              continue

                  iyi = iae(i,j1)
                  if(iyi .gt. 0) then
                     do 63 iy=1,3
                        if( iae(iyi,iy) == i) then
                           iyij = iy
                        endif
 63                  enddo
                  else
                     do ib1 =1,AMA%nbelm
                        if(itc(ib1) == i) then
                           itcii = ib1
                           goto 663
                        endif
                     enddo
                     print *,'boundary segment in DELANAY for itc not2'
                     stop
                  endif
 663              continue
                  if(itci .gt. 0) itc(itci) = i
                  if(itcii .gt. 0) itc(itcii) = ii

                  lnd(i,j1) = k3
                  lnd(ii,jj1) = k1

                  iae(i,j) = iyii
                  iae(i,j1) = ii
                  if(iyii .gt.0 ) iae(iyii,iyiij) = i

                  iae(ii,jj) = iyi
                  iae(ii,jj1) = i
                  if(iyi .gt.0 ) iae(iyi,iyij) = ii

                  icha = icha+1
               endif
            endif
 20      enddo
 10   enddo
      return
    end subroutine DELAUNAY



    subroutine SWAPPING(icha, icy)
      implicit none
      integer, intent(in) ::  icy
      integer, intent(inout) ::  icha
      real :: x1, y1, x2, y2, x3, y3, x4, y4, det123, det134
      real :: acc_new, acc_old, epsround, rl13, rl24
      integer :: i, j, j1, j2, ii, jj, jjj, jj1, j0, k1, k2, k3, k4
      integer :: itci, itcii, iyii, iy, ib1, ice, itet, iyi, iyiij, iyij
      integer, dimension(:,:), pointer :: iba
      real, dimension(:), pointer :: rga, rgb, rgc
      integer, dimension(:,:), pointer :: icyc
      integer, dimension(:,:), pointer :: lnd, iae, lbn
      integer, dimension(:), pointer :: ibc, itc
      real, dimension(:), pointer :: x, y
      integer :: ic_start, ic_end, ic_skip

      x => AMA%x(1:AMA%mpoin)
      y => AMA%y(1:AMA%mpoin)

      lbn => AMA%lbn(1:AMA%mbelm,1:2)
      ibc => AMA%ibc(1:AMA%mbelm)
      itc => AMA%itc(1:AMA%mbelm)


      lnd => AMA%lnd(1:AMA%melem, 1:3)
      iae => AMA%iae(1:AMA%melem, 1:3)

      icyc => AMA%icyc(1:AMA%mpoin, 1:AMA%maxdeg)

      rga => AMA%rga( 1:AMA%mpoin )
      rgb => AMA%rgb( 1:AMA%mpoin )
      rgc => AMA%rgc( 1:AMA%mpoin )

      iba => AMA%iba(1:AMA%melem, 1:3)

      do i=1,AMA%nelem
         do j=1,3
            iba(i,j) = 0
         enddo
      enddo

      ice = 0
      icha = 0

      if(mod(icy, 2) == 0) then
         ic_start = 1;         ic_end = AMA%nelem;   ic_skip =  1
      else
         ic_start = AMA%nelem;  ic_end = 1;          ic_skip = -1
      endif

      !do 10 i=1,AMA%nelem
      do 10 i = ic_start, ic_end, ic_skip

         do 20 j=1,3

            if(iae(i,j) .gt. 0 .and. iba(i,j) == 0) then
               j1 = mod(j,3) +1
               j2 = mod(j1,3) +1
               ii = iae(i,j)
               do 30 jjj=1,3
                  if(lnd(ii,jjj) == lnd(i,j1)) jj = jjj
 30            enddo
               iba(i,j) = 1
               iba(ii,jj) = 1
               jj1 = mod(jj,3) +1
               j0 = mod(jj1,3) +1
               k1 = lnd(i,j2)
               k2 = lnd(i,j)
               k3 = lnd(ii,j0)
               k4 = lnd(i,j1)
               if( k2 .ne. lnd(ii,jj1) .or. k4 .ne. lnd(ii,jj) ) then
                  print *,'ERROR !@#'
                  print *,k1,k2,k3,k4
                  return
               endif
               x1 = x(k1)
               y1 = y(k1)
               x2 = x(k2)
               y2 = y(k2)
               x3 = x(k3)
               y3 = y(k3)
               x4 = x(k4)
               y4 = y(k4)

!     we prohibid SWAPPING in case, where can appear alement with two
!     boundary segment
               if( (iae(i,j2) .lt. 0 .and. iae(ii,jj1) .lt. 0 ) .or.  &
                    (iae(i,j1) .lt. 0 .and. iae(ii,j0) .lt. 0) ) then
                  goto 20
               endif
!     we must still chech the orientation, i.e. the  cykles of points
!     k1,k2,k4 and k1,k4,k3 must have positive orientation

               !det123 = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
               !det134 = x1*(y3-y4) + x3*(y4-y1) + x4*(y1-y3)
               !
               !reps123 = AMA%pos*( ((x1-x3)**2 + (y1-y3)**2) +  &
               !     ((x2-x3)**2 + (y2-y3)**2) +  &
               !     ((x1-x2)**2 + (y1-y2)**2) )
               !reps134 = AMA%pos*( ((x1-x3)**2 + (y1-y3)**2) +  &
               !     ((x4-x3)**2 + (y4-y3)**2) +  &
               !     ((x1-x4)**2 + (y1-y4)**2) )
               !if( det123 .le. reps123 .or. det134 .le. reps134) then


               det123 = POS_TEST(x1, y1, x2, y2, x3, y3 )
               det134 = POS_TEST(x1, y1, x3, y3, x4, y4 )
               !print *,'..SWA.',det

               if( det123 .le. 1. .or. det134 .le. 1.) then
!                  print *,'violation of positivity'
                  goto 20
               endif

               itet = 0
               call POS1TEST(x(k1),y(k1),x(k2),y(k2),x(k3),y(k3),itet)
               if(itet == 1) then
                  goto 20
               endif
               call POS1TEST(x(k1),y(k1),x(k3),y(k3),x(k4),y(k4),itet)
               if(itet == 1) then
                  goto 20
               endif

               acc_old = ACCUTE(x(k1),y(k1),x(k2),y(k2),x(k4),y(k4),1.)  &
                    *ACCUTE(x(k2),y(k2),x(k3),y(k3),x(k4),y(k4),1. )

               acc_new = ACCUTE(x(k1),y(k1),x(k2),y(k2),x(k3),y(k3),1. )  &
                    *ACCUTE(x(k1),y(k1),x(k3),y(k3),x(k4),y(k4),1. )
               ice = 0

               epsround = 5E-03

               rl13 = ( (rga(k1)+rga(k3))*(x1-x3)*(x1-x3) +  &
                    2*(rgb(k1)+rgb(k3))*(x1-x3)*(y1-y3)  +  &
                    (rgc(k1)+rgc(k3))*(y1-y3)*(y1-y3) )/2

               rl24 = ( (rga(k2)+rga(k4))*(x2-x4)*(x2-x4) +  &
                    2*(rgb(k2)+rgb(k4))*(x2-x4)*(y2-y4)  +  &
                    (rgc(k2)+rgc(k4))*(y2-y4)*(y2-y4) )/2


               if(abs(rl13 - 3.) .lt. 0.995 * abs(rl24 -3.) )then
!     ... NEW ACCUTE
!               if(abs(rl13 - 3.)*acc_new .lt.
!     *              0.995 * abs(rl24 -3.)*acc_old )then

!        ... we use a SWAPPING
                  itci = 0
                  itcii = 0

                  iyii = iae(ii,jj1)
                  if(iyii .gt. 0) then
                     do 62 iy=1,3
                        if( iae(iyii,iy) == ii) then
                           iyiij = iy
                        endif
 62                  enddo
                  else
!     boundary segment, we seek itc
                     do ib1 =1,AMA%nbelm
                        if(itc(ib1) == ii) then
                           itci = ib1
                           goto 662
                        endif
                     enddo
                     print *,'boundary segment in SWAPPING for itc not'
                     stop
                  endif
 662              continue

                  iyi = iae(i,j1)
                  if(iyi .gt. 0) then
                     do 63 iy=1,3
                        if( iae(iyi,iy) == i) then
                           iyij = iy
                        endif
 63                  enddo
                  else
                     do ib1 =1,AMA%nbelm
                        if(itc(ib1) == i) then
                           itcii = ib1
                           goto 663
                        endif
                     enddo
                     print *,'boundary segment in SWAPPING for itc not2'
                     stop
                  endif
 663              continue
                  if(itci .gt. 0) itc(itci) = i
                  if(itcii .gt. 0) itc(itcii) = ii

                  lnd(i,j1) = k3
                  lnd(ii,jj1) = k1

                  iae(i,j) = iyii
                  iae(i,j1) = ii
                  if(iyii .gt.0 ) iae(iyii,iyiij) = i

                  iae(ii,jj) = iyi
                  iae(ii,jj1) = i
                  if(iyi .gt.0 ) iae(iyi,iyij) = ii
                  icha = icha+1

                     if(AMA%ifig .ge. 0 .and. mod(icha,5) == 0 ) then
                        call PLOT1()
                        AMA%ifig = AMA%ifig + 1
                     endif

               endif
            endif
 20      enddo
 10   enddo

    end subroutine SWAPPING

    end module
