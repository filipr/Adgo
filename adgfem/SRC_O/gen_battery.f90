  program gen_battery
    integer :: igrid, nn,ix, iy, mx, my, ipp, iee
    integer :: nelem, npoin, nbelm, ic
    real :: x(10), y(10), xx, yy
    real :: xs(4), ys(4), rgb(1:10, 1:3)
    integer :: nx(0:10), ny(0:10)
    real, dimension(:,:), allocatable :: xi
    integer, dimension(:,:), allocatable :: lnd, lbn
    character(len=4) ch1
    character(len=80) ch80

    COMMON/SIZES/DIM(4)
    DATA DIM/3, 9.,1.,1./   ! battery

    nx(:) = 0
    ny(:) = 0

    ch1 = '    '
    x(1) = 0.
    x(2) = 6.1
    x(3) = 6.5
    x(4) = 8.0
    x(5) = 8.4
    
    y(1) = 0.
    y(2) = 0.8
    y(3) = 1.6
    y(4) = 3.6
    y(5) = 18.8
    y(6) = 21.2
    y(7) = 23.2
    y(8) = 24.0

    igrid = 10

    ! file for the constrains
    open(igrid,file='battery_lines', status='UNKNOWN')
    write(igrid, *) 15
    write(igrid, *) 'NONE'
    write(igrid, '(i5, 4f10.2)')  1, x(1), y(2), x(2), y(2)
    write(igrid, '(i5, 4f10.2)')  2, x(2), y(2), x(3), y(2)
    write(igrid, '(i5, 4f10.2)')  3, x(3), y(2), x(4), y(2)
    write(igrid, '(i5, 4f10.2)')  4, x(4), y(2), x(4), y(7)
    write(igrid, '(i5, 4f10.2)')  5, x(1), y(7), x(4), y(7)

    write(igrid, '(i5, 4f10.2)')  6, x(1), y(3), x(2), y(3)
    write(igrid, '(i5, 4f10.2)')  7, x(1), y(4), x(2), y(4)
    write(igrid, '(i5, 4f10.2)')  8, x(1), y(5), x(2), y(5)
    write(igrid, '(i5, 4f10.2)')  9, x(1), y(6), x(2), y(6)
    write(igrid, '(i5, 4f10.2)') 10, x(2), y(6), x(3), y(6)

    write(igrid, '(i5, 4f10.2)') 11, x(3), y(2), x(3), y(6)
    write(igrid, '(i5, 4f10.2)') 12, x(2), y(2), x(2), y(3)
    write(igrid, '(i5, 4f10.2)') 13, x(2), y(3), x(2), y(4)
    write(igrid, '(i5, 4f10.2)') 14, x(2), y(4), x(2), y(5)
    write(igrid, '(i5, 4f10.2)') 15, x(2), y(5), x(2), y(6)

    close(igrid)

    open(igrid,file='battery', status='UNKNOWN')
    write(igrid,*)  x(1), y(1)
    write(igrid,*)  x(5), y(1)
    write(igrid,*)  x(5), y(8)
    write(igrid,*)  x(1), y(8)
    write(igrid,*)  x(1), y(1)
    write(igrid,'(x)')

    write(igrid,*)  x(1), y(2)
    write(igrid,*)  x(4), y(2)
    write(igrid,*)  x(4), y(7)
    write(igrid,*)  x(1), y(7)
    write(igrid,*)  x(1), y(2)
    write(igrid,'(x)')

    write(igrid,*)  x(3), y(2)
    write(igrid,*)  x(3), y(6)
    write(igrid,*)  x(1), y(6)
    write(igrid,'(x)')

    write(igrid,*)  x(2), y(2)
    write(igrid,*)  x(2), y(6)
    write(igrid,'(x)')

    write(igrid,*)  x(1), y(3)
    write(igrid,*)  x(2), y(3)
    write(igrid,'(x)')

    write(igrid,*)  x(1), y(4)
    write(igrid,*)  x(2), y(4)
    write(igrid,'(x)')

    write(igrid,*)  x(1), y(5)
    write(igrid,*)  x(2), y(5)
    write(igrid,'(x)')


    close(igrid)

    rgb(1, 1:3) = (/1., 0., 0. /)
    rgb(2, 1:3) = (/0., 1., 0. /)
    rgb(3, 1:3) = (/1., 0., 1. /)
    rgb(4, 1:3) = (/1., 1., 0. /)
    rgb(5, 1:3) = (/0., 0., 1. /)
    rgb(6, 1:3) = (/0., 1., 1. /)

    CALL BEGIN_GGG('Gfigure')
    CALL MICKEY_MOUSE(1, 1, x(1), x(5), y(1), y(8) )
    CALL FRAME

    xs(1) = x(1);   xs(2) = x(2);     xs(3) = xs(2);     xs(4) = xs(1); 
    ys(1) = y(2);   ys(3) = y(3);     ys(2) = ys(1);     ys(4) = ys(3); 
    ic = 5

    CALL RGB_CLOSE_FILL(xs(1:4), ys(1:4), 4,  rgb(ic,1), rgb(ic, 2), rgb(ic, 3) , .false.)
    xs(1) = x(1);   xs(2) = x(2);     xs(3) = xs(2);     xs(4) = xs(1); 
    ys(1) = y(3);   ys(3) = y(4);     ys(2) = ys(1);     ys(4) = ys(3); 
    ic = 2
    CALL RGB_CLOSE_FILL(xs(1:4), ys(1:4), 4,  rgb(ic,1), rgb(ic, 2), rgb(ic, 3) , .false.)

    xs(1) = x(1);   xs(2) = x(2);     xs(3) = xs(2);     xs(4) = xs(1); 
    ys(1) = y(4);   ys(3) = y(5);     ys(2) = ys(1);     ys(4) = ys(3); 
    ic = 3
    CALL RGB_CLOSE_FILL(xs(1:4), ys(1:4), 4,  rgb(ic,1), rgb(ic, 2), rgb(ic, 3) , .false.)

    xs(1) = x(1);   xs(2) = x(2);     xs(3) = xs(2);     xs(4) = xs(1); 
    ys(1) = y(5);   ys(3) = y(6);     ys(2) = ys(1);     ys(4) = ys(3); 
    ic = 2
    CALL RGB_CLOSE_FILL(xs(1:4), ys(1:4), 4,  rgb(ic,1), rgb(ic, 2), rgb(ic, 3) , .false.)

    xs(1) = x(2);   xs(2) = x(3);     xs(3) = xs(2);     xs(4) = xs(1); 
    ys(1) = y(2);   ys(3) = y(6);     ys(2) = ys(1);     ys(4) = ys(3); 
    ic = 4
    CALL RGB_CLOSE_FILL(xs(1:4), ys(1:4), 4,  rgb(ic,1), rgb(ic, 2), rgb(ic, 3) , .false.)

    xs(1) = x(1);   xs(2) = x(4);     xs(3) = xs(2);     xs(4) = xs(1); 
    ys(1) = y(6);   ys(3) = y(7);     ys(2) = ys(1);     ys(4) = ys(3); 
    ic = 5
    CALL RGB_CLOSE_FILL(xs(1:4), ys(1:4), 4,  rgb(ic,1), rgb(ic, 2), rgb(ic, 3) , .false.)

    xs(1) = x(3);   xs(2) = x(4);     xs(3) = xs(2);     xs(4) = xs(1); 
    ys(1) = y(2);   ys(3) = y(7);     ys(2) = ys(1);     ys(4) = ys(3); 
    ic = 5
    CALL RGB_CLOSE_FILL(xs(1:4), ys(1:4), 4,  rgb(ic,1), rgb(ic, 2), rgb(ic, 3) , .false.)

    xs(1) = x(1);   xs(2) = x(5);     xs(3) = xs(2);     xs(4) = xs(1); 
    ys(1) = y(1);   ys(3) = y(2);     ys(2) = ys(1);     ys(4) = ys(3); 
    ic = 1
    CALL RGB_CLOSE_FILL(xs(1:4), ys(1:4), 4,  rgb(ic,1), rgb(ic, 2), rgb(ic, 3) , .false.)

    xs(1) = x(1);   xs(2) = x(5);     xs(3) = xs(2);     xs(4) = xs(1); 
    ys(1) = y(7);   ys(3) = y(8);     ys(2) = ys(1);     ys(4) = ys(3); 
    ic = 1
    CALL RGB_CLOSE_FILL(xs(1:4), ys(1:4), 4,  rgb(ic,1), rgb(ic, 2), rgb(ic, 3) , .false.)

    xs(1) = x(4);   xs(2) = x(5);     xs(3) = xs(2);     xs(4) = xs(1); 
    ys(1) = y(1);   ys(3) = y(8);     ys(2) = ys(1);     ys(4) = ys(3); 
    ic = 1
    CALL RGB_CLOSE_FILL(xs(1:4), ys(1:4), 4,  rgb(ic,1), rgb(ic, 2), rgb(ic, 3) , .false.)


    CALL GREATER_BOUNDINGBOX(1.5,0.3,0.5,0.5)


    CALL END_GGG
    nn = 2
    
    ! horizontal
    nx(0) = 1
    nx(1) = 10 * nn
    nx(2) = 1* nn
    nx(3) = 2* nn
    nx(4) = 1* nn

    ! vertical
    ny(0) = 1
    ny(1) = 1 * nn
    ny(2) = 1 * nn
    ny(3) = 2 * nn
    ny(4) = 10 * nn
    ny(5) = 2 * nn
    ny(6) = 2 * nn
    ny(7) = 1 * nn

    ix = 4    ! number of panels
    iy = 7    ! number of panels

    mx = sum(nx(1:ix) ) 
    my = sum(ny(1:iy) ) 

    print*,'mx = ', mx
    print*,'my = ', my

    npoin = (mx + 1) * (my + 1)
    nelem = mx * my * 2
    nbelm = 2 *(mx + my)


    allocate(xi(1:npoin, 1:2) ) 
    allocate(lnd(1:nelem, 1:3) ) 

    ! nodes
    
    ipp = 0
    iee = 0 
    do i=0, ix
       do ii = 1, nx(i)
          if(i== 0) then
             xx = x(1)
          else
             xx = x(i) + 1. * ii / nx(i) *(x(i+1) - x(i) )
          endif
          
          do j=0, iy
             
             do jj = 1, ny(j)
                if(j== 0) then
                   yy = y(1)
                else
                   yy = y(j) + 1. * jj / ny(j) *(y(j+1) - y(j) )
                endif

                ipp = ipp + 1
                xi(ipp, 1) = xx
                xi(ipp, 2) = yy

                !print*,'???', i,j,ii,jj,ipp, xi(ipp, :)
                
                if(i > 0 .and. j > 0) then
                   iee = iee + 1
                   lnd(iee, 1) = ipp
                   lnd(iee, 2) = ipp - my - 1
                   lnd(iee, 3) = ipp - my - 2

                   iee = iee + 1
                   lnd(iee, 1) = ipp
                   lnd(iee, 2) = ipp - my - 2
                   lnd(iee, 3) = ipp - 1
                endif


             enddo
          enddo
       enddo
    enddo


    open(igrid, file = 'battery.grid', status = 'UNKNOWN')

    write(igrid, *)  npoin, nelem, nbelm, 4
    write(igrid, *)  0.,  0., 0, 0,  0., 0., 0, 0

    do ipp = 1, npoin
       write(igrid,*) xi(ipp, 1:2)
    enddo

    do ipp = 1, nelem
       write(igrid,*) 3, lnd(ipp, 1:3)
    enddo

    ! boundary
    ! left
    do ipp = my, 1, - 1
       write(igrid,*) ipp+1, ipp, 4
    enddo


    ! bottom
    do ipp = 1, mx
       write(igrid,*) 1 + (ipp-1)*(my+1),  1 + (ipp )*(my+1),  1
    enddo

    ! right
    do ipp =  my, 1, -1
       write(igrid,*) (mx+1)*(my+1) -ipp , (mx+1)*(my+1) - ipp +1 , 2
    enddo


    ! top
    do ipp = mx, 1, -1
       write(igrid,*) (ipp+1)*(my+1),  (ipp )*(my+1),  3
    enddo

  end program gen_battery
