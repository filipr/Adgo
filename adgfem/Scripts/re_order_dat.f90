  program repair_order_dat
    integer:: N_max, R_max
    real, dimension(:,:), allocatable :: x    
    integer, dimension(:,:), allocatable :: xi
    character*50 in_file, out_file

    integer:: i, ires, ires1, imesh

    if (command_argument_count() == 2) then
       call get_command_argument(1,in_file)
       
       call get_command_argument(2,out_file)
       
    else
       print*,'% Syntax:  re_order_dat  <input order.dat> <output orde.dat>'
       stop
    endif
    
    ires = 11
    open(ires, file=in_file,status='OLD')

    ires1 = 12
    open(ires1, file=out_file,status='UNKNOWN')

    
    N_max = 500
    R_max = 70

    allocate(x(1:N_max, 1:R_max), xi(1:N_max, 1:R_max)) 

    do k = 1,N_max
       read(ires, *, end=110) xi(k,1), x(k,2:3), xi(k,4:9), x(k,10:66)
    enddo
110 imesh = k-1


    do k = 1,imesh-1

       if( abs(x(k, 10) - x(k+1, 10)) > 1E-7) then
          write(ires1, 200 ) &
               xi(k,1), x(k,2:3), xi(k,4:9), x(k,10:67)

       endif
    enddo

200 format (i6,2es14.6,i9, 3i4, 2i5, es11.3, 60es16.6)


    close(ires)
    close(ires1)

  end program repair_order_dat
