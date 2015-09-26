module angener_sub
  use main_data  ! contains "type(mesh) ::  grid"   for computation
  use AMAdata

  implicit none

contains

  subroutine SEEK_NEIGH(npoin, nelem, lnd, iae, icyc)
    integer, intent(in) :: npoin, nelem
    integer, dimension(1:nelem, 1:3), intent(in) :: lnd
    integer, dimension(1:nelem, 1:3), intent(out) :: iae
    integer, dimension(1:AMA%mpoin, 1:AMA%maxdeg), intent(inout) :: icyc
    integer :: i,j, ii, il, j1, jj, jj1, k
    !     ... seeking of neighbours

    do i=1,nelem
       do j=1,3
          iae(i,j) = -2
       enddo
    enddo

    do i=1,npoin
       icyc(i,1) = 0
    enddo

    !      print*,'23'
    do i=1,nelem
       do j=1,3
          k = lnd(i,j)
          icyc(k,1) = icyc(k,1) + 1
          if(icyc(k,1) .ge. AMA%maxdeg) then
             print *,'Bad dimension in ADJAC'
             print *,'maxdeg < icyc(k,1)',AMA%maxdeg,icyc(k,1)
             stop
          endif
          icyc(k,icyc(k,1)+1) = i
       enddo
    enddo
    !      print*,'24'

    do i=1,nelem
       do  j=1,3
          if(iae(i,j) == -2)then
             j1 = mod(j,3) +1
             do il=1,icyc(lnd(i,j),1)
                ii = icyc(lnd(i,j),il+1)
                if( ii .ne. i) then
                   do jj=1,3
                      jj1 = mod(jj,3) +1
                      if( (lnd(i,j) == lnd(ii,jj1) ).and.  &
                           (lnd(i,j1) == lnd(ii,jj))) then
                         !     print *,'**',i,j,lnd(i,j),lnd(i,j1),ii
                         iae(i,j) = ii
                         iae(ii,jj) = i
                         goto 60
                      endif
                   enddo
                endif
             enddo
60           continue
          endif
       enddo
    enddo

  end subroutine SEEK_NEIGH

  !>     ... compute the inner angles of a triangle [x1,y1], [x2,y2], [x3,y3]
  !>    ... and return \f$ 1 + \sum_{j=1}^3 ( (\cos a_i)^-)^2 \f$
  function ACCUTE_I(x1,y1,x2,y2,x3,y3,iedge,par)
    real :: ACCUTE_I
    real, intent(inout) :: x1,y1,x2,y2,x3,y3
    real, intent(in) ::  par
    integer, intent(inout) :: iedge
    real x(3),y(3),vax,vay,vbx,vby,ccos,ang
    integer :: j0, j1, j2

    ACCUTE_I = 1.
    x(1) = x1
    x(2) = x2
    x(3) = x3
    y(1) = y1
    y(2) = y2
    y(3) = y3


    j0 = mod(iedge,3) + 1
    j1 = mod(j0,3) + 1
    j2 = mod(j1,3) + 1
    vax = x(j0) - x(j1)
    vay = y(j0) - y(j1)
    vbx = x(j2) - x(j1)
    vby = y(j2) - y(j1)
    ccos = (vax*vbx + vay*vby)
    if(ccos .lt. 0.D+0) then
       ang = ccos*ccos/(vax*vax + vay*vay)/(vbx*vbx + vby*vby)
       ACCUTE_I = ACCUTE_I+12.*par*ang
    endif

  end function ACCUTE_I

  !>     ... compute the inner angles of a triangle [x1,y1], [x2,y2], [x3,y3]
  !>     ... and return  \f$ 1 + \sum_{j=1}^3 ( (\cos a_i)^-)^2  \f$
  function ACCUTE(x1,y1,x2,y2,x3,y3,par)
    real :: ACCUTE
    real, intent(inout) :: x1,y1,x2,y2,x3,y3
    real, intent(in) :: par
    real x(3),y(3),vax,vay,vbx,vby,ccos,ang
    integer :: j0, j1, j2

    ACCUTE = 1.
    x(1) = x1
    x(2) = x2
    x(3) = x3
    y(1) = y1
    y(2) = y2
    y(3) = y3


    do j0=1,3
       j1 = mod(j0,3) + 1
       j2 = mod(j1,3) + 1
       vax = x(j0) - x(j1)
       vay = y(j0) - y(j1)
       vbx = x(j2) - x(j1)
       vby = y(j2) - y(j1)
       ccos = (vax*vbx + vay*vby)
       if(ccos .lt. 0.D+0) then
          ang = ccos*ccos/(vax*vax + vay*vay)/(vbx*vbx + vby*vby)
          ACCUTE = ACCUTE+12.*par*ang
       endif
    enddo
  end function ACCUTE


  subroutine POS1TEST(x1,y1,x2,y2,x3,y3,itet)
    real, intent(inout) :: x1,y1,x2,y2,x3,y3
    integer, intent(inout) :: itet
    real :: rl1, rl2, rl3
    integer :: itet1

    if( AMA%pos1 .le. 1E-10) then
       itet = 0
       return
    endif

    if(itet == -1) then
       itet1 = -1
    else
       itet1 = 0
    endif
    itet1 = 0

    if(itet1 == -1) then
       print *,'***  ---'
       print *,x1,y1
       print *,x2,y2
       print *,x3,y3
    endif
    itet = 0
    rl1 = ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))**0.5
    rl2 = ((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3))**0.5
    rl3 = ((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1))**0.5
    if(rl1 + rl2 .lt. (1.+AMA%pos1)*rl3) itet = 1
    if(rl2 + rl3 .lt. (1.+AMA%pos1)*rl1) itet = 1
    if(rl3 + rl1 .lt. (1.+AMA%pos1)*rl2) itet = 1

    if(itet1 == -1) print *,'* *',itet

    if(itet /= 0) print*,'##POS1TEST ',AMA%pos1, itet
  end subroutine POS1TEST

  subroutine POS2TEST(x1,y1,x2,y2,x3,y3,itet)
    real, intent(inout) :: x1,y1,x2,y2,x3,y3
    integer, intent(inout) :: itet
    real :: rl1, rl2, rl3
    integer :: itet1

    if(itet == -1) then
       itet1 = -1
    else
       itet1 = 0
    endif
    if(itet1 == -1) then
       print *,'***'
       print *,x1,y1
       print *,x2,y2
       print *,x3,y3
    endif
    itet = 0
    rl1 = ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))**0.5
    rl2 = ((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3))**0.5
    rl3 = ((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1))**0.5
    if(rl1 + rl2 .lt. (1.+AMA%pos1)*rl3) itet = 1
    if(rl2 + rl3 .lt. (1.+AMA%pos1)*rl1) itet = 1
    if(rl3 + rl1 .lt. (1.+AMA%pos1)*rl2) itet = 1
    if(itet1 == -1) print *,'* *',itet
  end subroutine POS2TEST

  !> testing the positivity, if POST_TEST > 1 then traingle is OK else refuse adaptation
  function POS_TEST(x1, y1, x2, y2, x3, y3)
    real :: POS_TEST
    real, intent(inout) :: x1, y1, x2, y2, x3, y3
    real :: det, reps, a12, a23, a31, b12, b23, b31
    real :: l12, l23, l31, s1, s2, s3, g1, g2, g3
    real :: pi

    !pi = 2* acos(0.)

    !print *,'!!',x1, y1
    !print *,'!!',x2, y2
    !print *,'!!',x3, y3

    !x3 = 0.
    !y3 = 1.
    !x1=  0.
    !y1 = 0.
    !x2 = 1.
    !y2 = 0.

    a12 = x1 - x2
    a23 = x2 - x3
    a31 = x3 - x1

    b12 = y1 - y2
    b23 = y2 - y3
    b31 = y3 - y1

    l12 = a12 * a12 + b12 * b12
    l23 = a23 * a23 + b23 * b23
    l31 = a31 * a31 + b31 * b31

    ! scalar product of two vector sharing the vertex 1
    s1 = - a31 * a12 - b31 * b12
    s2 = - a12 * a23 - b12 * b23
    s3 = - a23 * a31 - b23 * b31

    det = x1 * b23 + x2* b31 + x3 * b12

    reps = l12 + l23 + l31

    ! variant of cosine theorem
    !g1 = - (l23 - l12 - l31  ) / ( 2 *(l12 * l31)**0.5 )

    ! variant of scalar product
    g1 = s1 / ((l31 * l12)**0.5)
    g2 = s2 / ((l12 * l23)**0.5)
    g3 = s3 / ((l23 * l31)**0.5)

    !print *,'!!',x3, y3
    !print *,'!!',x1, y1
    !print *,'!!',x2, y2

    !print*,'###',g1, g2, g3
    !print*,'###',acos(g1)/pi * 180 , acos(g2)/pi * 180, acos(g3)/pi * 180

    !stop
    if(min(g1, g2, g3) < 0.) then ! non-accute angle, cos(angle) < 0.
       POS_TEST = det / ( reps * AMA%posW)   ! stronger condition for non-accute angles
    else
       POS_TEST = det / ( reps * AMA%pos)
    endif

  end function POS_TEST


    subroutine ELIPS(ieli,a,b,c,xi,yi)
      integer, intent(in) :: ieli
      real, intent(in)  :: a,b,c,xi,yi
      integer :: num, i
      real :: t, x, y, rnorm, xs, ys

      num = 100
      do i=0,num
         t = 1.*i/num*6.283185307
         x = cos(t)
         y = sin(t)
         rnorm = (x*x*a + 2*x*y*b + y*y*c)**0.5
         xs = x/rnorm + xi
         ys = y/rnorm + yi
         write(ieli,*) xs,ys, i, xi,yi, ieli
         !write(*,*) xs,ys, i, xi,yi, ieli
      enddo
      write(ieli, *)'## '
      write(ieli, *)
      return
    end subroutine ELIPS




    subroutine QUA2(x1,y1,w1,ra1,rb1,rc1,x2,y2,w2,ra2,rb2,rc2,  &
         x3,y3,w3,ra3,rb3,rc3,err1,err2,ice,ibo )
      real, intent(in) :: x1,y1,w1,ra1,rb1,rc1,x2,y2,w2,ra2,rb2,rc2,  &
         x3,y3,w3,ra3,rb3,rc3
      real, intent(inout) :: err1, err2
      integer, intent(in) :: ice, ibo
      real:: sqrt3, rl1, rl2, rl3

      err1 = 1.
      sqrt3 = 1.7320508075

!      print*, 'ra1, ra2: ' , ra1 , ra2
!      print*, 'rb1, rb2: ' , rb1, rb2
!      print*, 'x1,x2,y1,y2: ', x1,x2,y1,y2

      rl1 = (((ra1+ra2)*(x1-x2)*(x1-x2) + 2*(rb1+rb2)*(x1-x2)*(y1-y2) +  &
           (rc1+rc2)*(y1-y2)*(y1-y2) )/2 )**0.5
      rl2 = (((ra2+ra3)*(x2-x3)*(x2-x3) + 2*(rb2+rb3)*(x2-x3)*(y2-y3) +  &
           (rc2+rc3)*(y2-y3)*(y2-y3) )/2 )**0.5
      rl3 = (((ra1+ra3)*(x1-x3)*(x1-x3) + 2*(rb1+rb3)*(x1-x3)*(y1-y3) +  &
           (rc1+rc3)*(y1-y3)*(y1-y3) )/2 )**0.5

      err2 = (sqrt3-rl1)**2 + (sqrt3-rl2)**2 + (sqrt3-rl3)**2

      if(ibo == 1) then
         err2 = err2 + (sqrt3-rl1)**2
      elseif(ibo == 2) then
         err2 = err2 + (sqrt3-rl2)**2
      elseif(ibo == 3) then
         err2 = err2 + (sqrt3-rl3)**2
      endif


      AMA%glmax = max(AMA%glmax, rl1**2,rl2**2,rl3**2)
      AMA%glmin = min(AMA%glmin, rl1**2,rl2**2,rl3**2)

      if(ice == 13) then
         write(*,'(4e14.6)') rl1,rl2,rl3,err2
         write(*,'(5e14.6)') x1,y1,ra1,rb1,rc1
         write(*,'(5e14.6)') x2,y2,ra2,rb2,rc2
         write(*,'(5e14.6)') x3,y3,ra3,rb3,rc3
      endif
      return
    end subroutine QUA2

end module angener_sub
