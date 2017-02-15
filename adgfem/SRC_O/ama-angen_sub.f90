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
    real, intent(in) :: x1,y1,x2,y2,x3,y3
    real, intent(in) :: par
    real x(3),y(3),vax,vay,vbx,vby,ccos,ang
    integer :: j0, j1, j2

    ACCUTE = 0.
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
          !ACCUTE = ACCUTE+ par*ang
          ACCUTE = acos(ccos / sqrt((vax*vax + vay*vay)* (vbx*vbx + vby*vby)) )
       endif

       !if( AMA%ifig == 9 )write(*,'(a10, i5, 8es12.4)') 'ACCUTE',AMA%ifig, ccos, ang , ACCUTE
    enddo
    !if( AMA%ifig == 9 )write(*,'(a20, 6es12.4)' )'././././././.,m,.', x1,y1,x2,y2,x3,y3
  end function ACCUTE

  !>     ... compute the inner angle of a triangle [x1,y1], [x2,y2], [x3,y3] at [x2, y2]
  function ANGLE(x1,y1,x2,y2,x3,y3)
    real :: ANGLE
    real, intent(in) :: x1,y1,x2,y2,x3,y3
    real x(3),y(3),vax,vay,vbx,vby,ccos,ang
    integer :: j0, j1, j2

    
    vax = x1 - x2
    vay = y1 - y2
    vbx = x3 - x2
    vby = y3 - y2
    
    if( (vax*vax + vay*vay) * (vbx*vbx + vby*vby) <= 0.) then
       ANGLE = 2* acos(0.)  ! patological case, we prohibit this adaptation operation
    else
       ccos = (vax*vbx + vay*vby)
       ANGLE = acos( max(-0.9999, ccos / sqrt( (vax*vax + vay*vay) * (vbx*vbx + vby*vby) ) ) )
    endif

  end function ANGLE


  !>     ... compute the inner angles of a triangle [x1,y1], [x2,y2], [x3,y3] at [x2, y2]
  !> first angle is oposite to the first edge, i.e, angle at [x3, y3]
  !> second angle is oposite to the 2nd edge, i.e, angle at [x3, y3]
  !> third angle is oposite to the 3rd edge, i.e, angle at [x3, y3]
  subroutine ALL_ANGLES(x1,y1,x2,y2,x3,y3, ang1, ang2, ang3)
    real, intent(in) :: x1,y1,x2,y2,x3,y3
    real, intent(inout) :: ang1, ang2, ang3

    ang2 = ANGLE(x3,y3,x1,y1,x2,y2)
    ang3 = ANGLE(x1,y1,x2,y2,x3,y3)
    ang1 = ANGLE(x2,y2,x3,y3,x1,y1)

  end subroutine ALL_ANGLES



  subroutine POS1TEST(x1,y1,x2,y2,x3,y3,itet, val)
    real, intent(inout) :: x1,y1,x2,y2,x3,y3
    integer, intent(inout) :: itet
    real, intent(inout), optional :: val
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

    if( present(val) ) then
       val = 0.
       val = min( (rl1 + rl2)/rl3, (rl2 + rl3)/rl2 , (rl3 + rl1)/rl2 )
       val = val - 1

    endif


    if(itet1 == -1) print *,'* *',itet

    !if(itet /= 0) print*,'##POS1TEST ',AMA%pos1, itet
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
    real :: pi, ccos

    !print*,'POS_TEST'
    pi = 2* acos(0.)

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

    !if(l12 <= 0. .or. l23 <= 0.) then
    !   write(*,'(a8, 50es12.4)') 'de3d',l12,l23
    !   print*, x1,y1
    !   print*, x2,y2
    !   print*, x3,y3
    !endif

    ! variant of scalar product
    g1 = s1 / max( sqrt(l31 * l12), 1E-15)
    g2 = s2 / max( sqrt(l12 * l23), 1E-15)
    g3 = s3 / max( sqrt(l23 * l31), 1E-15)

    !print *,'!!',x3, y3
    !print *,'!!',x1, y1
    !print *,'!!',x2, y2

    !print*,'###',g1, g2, g3
    !print*,'###',acos(g1)/pi * 180 , acos(g2)/pi * 180, acos(g3)/pi * 180

    !stop
    if(min(g1, g2, g3) < 0.) then ! non-accute angle, cos(angle) < 0.
       POS_TEST = det / ( reps * AMA%posW)   ! stronger condition for non-accute angles
       !print*,'POS_TEST  V_1',g1,g2,g3
    else
       POS_TEST = det / ( reps * AMA%pos)
       !print*,'POS_TEST  V_2',g1,g2,g3
    endif


    return
    
    ! NEW TEST
    !print*,'#E#SW@',g1,g2,g3, min(g1, g2, g3)

    ccos = acos( max(-0.9999,  min(g1, g2, g3)) ) ! -0.9999 in order to avoid patological case

    !print*, ccos, pi ! ccos/pi*180

    if(ccos > AMA%maximal_angle) then  ! too large angle
       POS_TEST = 0.1

       !print*,'----------------------------------'
       !print *,x3, y3, ccos/pi*180, AMA%maximal_angle/pi*180,'yu38dh3i'
       !print *,x1, y1
       !print *,x2, y2

    !else
    !   POS_TEST = 2.
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




    subroutine QUA2(x1, y1, w1, ra1, rb1, rc1, &
         x2, y2, w2, ra2, rb2, rc2,  &
         x3, y3, w3, ra3, rb3, rc3, err1, err2, ice, ibo )
      real, intent(in) :: x1,y1,w1,ra1,rb1,rc1,x2,y2,w2,ra2,rb2,rc2,  &
         x3,y3,w3,ra3,rb3,rc3
      real, intent(inout) :: err1, err2
      integer, intent(in) :: ice, ibo
      real:: sqrt3, rl1, rl2, rl3, obtuse, pi2
      real :: alpha

      pi2 = acos(-1.) /2.

      err1 = 1.
      sqrt3 = sqrt(3.) !1.7320508075

      !print*, 'ra1, ra2: ' , ra1 , ra2
      !print*, 'rb1, rb2: ' , rb1, rb2
      !print*, 'x1,x2,y1,y2: ', x1,x2,y1,y2

      ! if(  (ra1+ra2)*(rc1+rc2) - (rb1+rb2)**2 <= 0. .or.  &
      !      (ra2+ra3)*(rc2+rc3) - (rb2+rb3)**2  <= 0. .or.  &
      !      (ra3+ra1)*(rc3+rc1) - (rb3+rb1)**2 <= 0. ) then

      !    print*,'??>>?>,', ( (((ra1+ra2)*(x1-x2)*(x1-x2) + 2*(rb1+rb2)*(x1-x2)*(y1-y2) +  &
      !         (rc1+rc2)*(y1-y2)*(y1-y2) )/2 )), ( (((ra2+ra3)*(x2-x3)*(x2-x3) + 2*(rb2+rb3)*(x2-x3)*(y2-y3) +  &
      !         (rc2+rc3)*(y2-y3)*(y2-y3) )/2 )),  ( (((ra1+ra3)*(x1-x3)*(x1-x3) + 2*(rb1+rb3)*(x1-x3)*(y1-y3) +  &
      !         (rc1+rc3)*(y1-y3)*(y1-y3) )/2 ))
      !    print*,'??>>?>,', &
      !         (ra1+ra2)*(rc1+rc2) - (rb1+rb2)**2, &
      !         (ra2+ra3)*(rc2+rc3) - (rb2+rb3)**2, &
      !         (ra3+ra1)*(rc3+rc1) - (rb3+rb1)**2
      !    print*,'??>>?>,', ra1*rc1 - rb1**2,  ra2*rc2 - rb2**2,  ra3*rc3 - rb3**2
      !    print*
      ! endif

      !rl1 = sqrt(abs (((ra1+ra2)*(x1-x2)*(x1-x2) + 2*(rb1+rb2)*(x1-x2)*(y1-y2) +  &
      !     (rc1+rc2)*(y1-y2)*(y1-y2) )/2 ))
      !rl2 = sqrt(abs (((ra2+ra3)*(x2-x3)*(x2-x3) + 2*(rb2+rb3)*(x2-x3)*(y2-y3) +  &
      !     (rc2+rc3)*(y2-y3)*(y2-y3) )/2 ))
      !rl3 = sqrt(abs (((ra1+ra3)*(x1-x3)*(x1-x3) + 2*(rb1+rb3)*(x1-x3)*(y1-y3) +  &
      !     (rc1+rc3)*(y1-y3)*(y1-y3) )/2 ))

      rl1 = sqrt( EDGE_NORM_SQUARE (x1, y1, x2, y2, ra1, rb1, rc1, ra2, rb2, rc2) )

      rl2 = sqrt( EDGE_NORM_SQUARE (x2, y2, x3, y3, ra2, rb2, rc2, ra3, rb3, rc3) )

      rl3 = sqrt( EDGE_NORM_SQUARE (x3, y3, x1, y1, ra3, rb3, rc3, ra1, rb1, rc1) )

      ! TEST in l^1-norm
      !err2 = (sqrt3-rl1)**2 + (sqrt3-rl2)**2 + (sqrt3-rl3)**2

      alpha = 2.0
      err2 = (abs(sqrt3-rl1))**alpha + (abs(sqrt3-rl2))**alpha + (abs(sqrt3-rl3))**alpha


      if(ibo == 1) then
         err2 = err2 + (abs(sqrt3-rl1))**alpha
      elseif(ibo == 2) then
         err2 = err2 + (abs(sqrt3-rl2))**alpha
      elseif(ibo == 3) then
         err2 = err2 + (abs(sqrt3-rl3))**alpha
      endif


      AMA%glmax = max(AMA%glmax, rl1**2,rl2**2,rl3**2)
      AMA%glmin = min(AMA%glmin, rl1**2,rl2**2,rl3**2)

      if(ice == 13) then
         write(*,'(4e14.6)') rl1,rl2,rl3,err2
         write(*,'(5e14.6)') x1,y1,ra1,rb1,rc1
         write(*,'(5e14.6)') x2,y2,ra2,rb2,rc2
         write(*,'(5e14.6)') x3,y3,ra3,rb3,rc3
      endif

      obtuse = ACCUTE(x1, y1, x2, y2, x3, y3, 1.)


      !if(obtuse > -1.) write(*,*)'detde5xswiw63', obtuse, acos(-1.), obtuse /  acos(-1.) * 180
      if(obtuse > 1. ) then
         obtuse = (obtuse - pi2) / pi2   ! for obtuse angles, the quatity obtuse is between 0 and 1
         !!write(*,*)'detde5xswiw63', obtuse

         !err2 = err2 * (1 +  15 * obtuse**2)    ! quality reflets the obtuse angles
         !err2 = err2 * (1 +  100 * obtuse**3)    ! quality reflets the obtuse angles
      endif

      return
    end subroutine QUA2


    !> compute the norm of the edge [x1, y1] ---> [x2,y2] in the metric 
    !> (a1, b1, c1) in [x1, y1] and (a2,b2,c2) in [x2, y2] 
    function EDGE_NORM_SQUARE (x1, y1, x2, y2, a1, b1, c1, a2, b2, c2)
      real :: EDGE_NORM_SQUARE
      real, intent(in) :: x1, y1, x2, y2, a1, b1, c1, a2, b2, c2
      real :: l1, l2, a, b, c

      a = (a1+a2)/2
      b = (b1+b2)   ! / 2
      c = (c1+c2)/2

      l1 = x1 - x2
      l2 = y1 - y2

      EDGE_NORM_SQUARE =  a*l1*l1 + b*l1*l2 + c*l2*l2  
    end function EDGE_NORM_SQUARE

end module angener_sub
