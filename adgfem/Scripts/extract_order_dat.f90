  program extract_order_dat
    integer:: nelem, dof, degP, degT, degQ
    real :: h, tau 
    real, dimension(:), allocatable :: errL2, errH1

    allocate(errL2(1:3), errH1(1:3))

    idat = 11
    ires = 12
    ires2 = 13

    open(idat, file = 'order.bak', status='OLD')
    open(ires, file = 'orderL2.dat', status='UNKNOWN')
    open(ires2, file = 'orderH1.dat', status='UNKNOWN')

    do i=1,3
       read(idat, *) nelem, h, tau, dof, degP, degT, degQ, errL2(i), errH1(i)
    enddo
    
    write(ires, *) nelem, h, tau, dof, degP, degT, 0, errL2(1:3), errH1(1:3)
    write(ires2, *) nelem, h, tau, dof, degP, degT,1, errH1(1:3)

    deallocate(errL2, errH1)
    close(idat)
    close(ires)
    close(ires2)

  end program extract_order_dat
