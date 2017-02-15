program remake_ord

    integer:: N_max, R_max
    real, dimension(:,:), allocatable :: x    
    integer, dimension(:,:), allocatable :: xi
    character*1 iord
    character*50 in_file, out_file, remesh_file
    integer:: p_max, dof_tot, ncount, dof_min, dof_max, dof_old, i_remesh
    real :: Eip_aver, Eip_max, tau_aver, EL8L2, e1, e2, e3, a1, a2, a3, a0
    integer :: dof, it
    integer:: i, ires, ires1, ires2, imesh, il, k, itab, Miter, itab1


    itab = 13
    open(itab, file='tab.tex',status='UNKNOWN', position='append')

    itab1 = 14


    if (command_argument_count() == 4) then
       call get_command_argument(1,in_file)
       
       call get_command_argument(2,out_file)

       call get_command_argument(3,remesh_file)

       call get_command_argument(4,iord)

       if(iord == "1") then
          open(itab1, file='data.dat',status='OLD')
          read(itab1, *) e1, e2, e3, dof, it
          close(itab1)
       endif

    elseif (command_argument_count() == 1) then
       call get_command_argument(1,in_file)

       if(in_file == "tab_init") then
          write(itab,*)  '%{\footnotesize           %  vetsi '
          write(itab,*)  ' {\scriptsize             %  mensi ' 
          write(itab,*)  ' %{\tiny                   %  nejmensi '
          write(itab,*)  ' \begin{tabular}{cc|rcc|c|rrr|cc|ccc|r} '
          write(itab,*)  ' \hline '
          write(itab,*)  ' &  & & & & & \multicolumn{3}{|c|}{{$\Nhpm$}}  & '
          write(itab,*)  '\multicolumn{2}{|c|}{$\| \EIpm \|_{L^2(\Omega)}$} &  '
          write(itab,*)  '\multicolumn{3}{|c|}{$\|u -u_h\|$} & comput. \\'
          write(itab,*)  '  $ k $ & $\omega $ '
          !write(itab,*)  ' & $M_it$  '
          write(itab,*)  '& $\# \tau_m$ & $\# \Thpm$ & $ p_K^{\max}$'
          write(itab,*)  '  & $\overline \tau_m$  '
          write(itab,*)  '& min & aver & max'
          write(itab,*)  ' & aver  ' 
          write(itab,*)  ' & max   '
          write(itab,*)  ' & ${L^\infty(L^2)}$  '
          write(itab,*)  ' & ${L^2(L^2)}$       '                      
          write(itab,*)  ' & ${L^2(H^1)}$       '                      
          write(itab,*)  ' & time(s)                     '                      
          write(itab,*)  ' \\   \hline  '

       elseif(in_file == "hline") then
          write(itab,*)  '   \hline  '

       elseif(in_file == "tab_end") then
          write(itab,*)  '   \hline  '
          write(itab,*)  '    \end{tabular} '
          write(itab,*)  '    } %%% end of font size   '

       end if

       close(itab)
       return

    else
       print*,'% Syntax:  re_order_dat  <input order.dat> <output orde.dat> <remesh.dat>'
       print*, 'or'
       print*,'% Syntax:  re_order_dat  tab_init/ hline / tab_end '
       stop
    endif

    il = 0
    do k=1, 50
       if(in_file(k:k) == ' ') then
          il = k
          goto 13
       endif
    enddo

13  continue

    if(il == 0) then
       print*,'Touble in remake_ord.f90'
       stop
    endif

    !print*,il
    !print*,in_file
    !print*, in_file(il-18:il-11)
    
    ires = 11
    open(ires, file=in_file,status='OLD')

    ires1 = 12
    open(ires1, file=out_file,status='UNKNOWN')

    ires2 = 21
    open(ires2, file=remesh_file,status='UNKNOWN')

    
    N_max = 50000
    R_max = 70

    allocate(x(1:N_max, 1:R_max), xi(1:N_max, 1:R_max)) 
    
    ! reading of data
    do k = 1,N_max
       read(ires, *, end=110) xi(k,1), x(k,2:3), xi(k,4:9), x(k,10:66)
    enddo
110 imesh = k-1



    ! processing of data
    dof_old = 0
    i_remesh = 0

    p_max = 0
    dof_tot = 0
    dof_max = 0
    dof_min = 1000000
    ncount = 0
    Eip_aver = 0.
    Eip_max = 0.
    tau_aver = 0.
    EL8L2 = 0.

    do k = 1,imesh-1

       !if( x(k+1, 10) >  x(k, 10) ) then
       if( xi(k+1, 9) >  xi(k, 9) ) then

          ncount = ncount + 1
          p_max = max(p_max, xi(k,5) )
          dof_tot = dof_tot +  xi(k,4)
          dof_max = max(dof_max,  xi(k,4))
          dof_min = min(dof_min,  xi(k,4))

          tau_aver =  tau_aver + x(k,3)
          Eip_aver =  Eip_aver + x(k,23)
          Eip_max =  max(Eip_max,  x(k,23))

          EL8L2 = max (EL8L2, x(k, 11) )
          x(k, 31) = EL8L2

          ! mesh was changed
          if( xi(k,4) /= dof_old) then
             dof_old =  xi(k,4)
             i_remesh = i_remesh + 1
             write(ires2, 200 ) xi(k,1), x(k,2:3), xi(k,4:9), x(k,10:67)
          endif

          write(ires1, 200 ) xi(k,1), x(k,2:3), xi(k,4:9), x(k,10:67)


       endif
    enddo

    ! last step
    k = imesh
    EL8L2 = max (EL8L2, x(k, 11) )
    x(k, 31) = EL8L2

    write(ires1, 200 ) xi(k,1), x(k,2:3), xi(k,4:9), x(k,10:67)
    !print*,k, x(k,61)


200 format (i6,2es14.6,i9, 3i4, 2i5, es11.3, 60es16.6)


    close(ires)
    close(ires1)
    close(ires2)

    
    Miter = int( 1. * xi(k,8) /k +0.5)
    write(itab, 95) xi(k,6), &
         !in_file(il-18:il-11), & 
         '$10^{-',in_file(il-11:il-11),'}$', & 
         !Miter,  &
         ncount*Miter, &
         i_remesh, &
         p_max, &
         tau_aver/ncount, &
         dof_min, int(1.*dof_tot/ncount +0.5) ,dof_max,  &
         Eip_aver / ncount, Eip_max, x(k,31:33), x(k,61)

    if(iord == '1') then
       a0 = 1. /dof * int(1.*dof_tot/ncount +0.5)
       a1 = 2*log(e1 / x(k,31)) / log(a0)
       a2 = 2*log(e2/ x(k,32) ) / log(a0)
       a3 = 2*log(e3/ x(k,33) ) / log(a0)
       write(*,*), a0, a1,a2,a3
       write(itab, 96) a1,a2,a3
    endif

    
    close(itab)


    open(itab1, file='data.dat',status='UNKNOWN')
    write(itab1, *)x(k,31:33),int(1.*dof_tot/ncount +0.5), ncount
    close(itab1)

!95 format (i7, '&', a10,'&', 6(i7,' & '), 6(es10.2,' & '), f10.1, '\\') 
95 format (i7, '&', a6,a2,a2,'&', 3(i7,' & '), es10.2,' &' , 3(i7,' & '), 5(es10.2,' & '), f10.1, '\\') 
96 format( ' \multicolumn{2}{c|}{(EOC)}  & & & & & & & & & & ', 3( '{\bf (',f7.2, ')} & '), '\\ ')
end program remake_ord
