program make_gnu
  integer :: i,k
  integer :: Pmax, Mmax
  real, dimension(:,:,:), allocatable :: x
  character(len=40) :: ch, ch1, ch2, ch3, ch4, ch5, ch6
  character(len=240) :: chlong
  character(len=4) :: chmet
  character(len=2) :: chnorm, chm
  character(len=8) :: chnorm1
  character(len=12) :: chnorm2
  character(len=1) :: chP
  character(len=8) :: color
  logical:: lcolor

  if (command_argument_count() == 1) then
    call get_command_argument(1, color)
  else
     print*,'Syntax make_gnu <color>'
     stop
  endif

  lcolor = .false.
  if(color == "RGB" .or. color == "color") lcolor = .true.


  Qmax = 1 !4 ! polynomial degree 0 - Qmax
  Mmax = 4 !??
  allocate(x(1:Pmax, 1:Mmax, 0:6) )
  
  ida0 = 10
  ida1 = 11
  ida2 = 12
  ida3 = 13
  ida4 = 14
 ! ida5 = 15
 ! ida6 = 16

  open(ida1, file='orderP0.T0.dat', status='UNKNOWN')
  open(ida1, file='orderP0.T1.dat', status='UNKNOWN')
  open(ida2, file='orderP0.T2.dat', status='UNKNOWN')
  open(ida3, file='orderP0.T3.dat', status='UNKNOWN')
  open(ida4, file='orderP0.T4.dat', status='UNKNOWN')
 ! open(ida5, file='orderP0.T4.L2.dat', status='UNKNOWN')
 ! open(ida6, file='orderP6.T1.L2.dat', status='UNKNOWN')


  do k = 0, Qmax
     do i=1, Mmax
        read(10 + k, *) N, x(k, i, 0), xx, N1, N2, N3, N4, x(k,i, 1:6) 

     enddo
  enddo


  !do k = 1, Pmax
  !   do i=1,Mmax
  !      write(*, '(7es12.4)')  x(k, i, 0:6) 
  !   enddo
  !enddo

  ignu = 100
  open(ignu, file = 'conv_order.gnu', status='UNKNOWN')

  if(lcolor) then
     write(ignu,* ) 'set terminal postscript eps color enhanced'
     do k = 1, Pmax
        write(ignu,* ) 'set linetype ',k,'  lw 3.5'
     enddo

  write(ignu,* ) 'set logscale'
  write(ignu,* ) 'set size 0.45, 0.63' 

  else
     write(ignu,* ) 'set terminal postscript eps enhanced'
  write(ignu,* ) 'set logscale'
  write(ignu,* ) 'set size 0.55, 0.8' 
  endif


  
  do meth = 1,3
     if(meth == 1) chmet = 'SIPG'
     if(meth == 2) chmet = 'NIPG'
     if(meth == 3) chmet = 'IIPG'

     do norm = 1,2
        inum =  3* (norm -1) +  meth

        if(norm == 1) chnorm = 'L2'
        if(norm == 2) chnorm = 'H1'
        if(norm == 1) chnorm1 = 'L^2-norm'
        if(norm == 2) chnorm2 = 'H^1-seminorm'
        
        
        write(ignu,* ) 'set key right bottom'
        !write(ignu,* ) 'set xlabel "', chmet,' - ',chnorm1,'" '
        write(ignu,* ) "set output 'case-",chmet,'-',chnorm,".eps' "
        write(ignu,* ) 'unset label'

        y0old = 1000.
        ival = 1
        do k=1, Pmax
           x1 = x(k, Mmax, 0)
           x2 = x(k, Mmax-2, 0)
           
           y1 = x(k, Mmax, inum)
           y2 = x(k, Mmax-2, inum)
           
           order = log(y2/y1) / log(x2/x1)
           
           io = int (order + 0.25)
           io1 = int (order + 0.75)
           
           x0 = x1 *0.9 
           y0 = y1
           
           print*,'####', order, io, io1, y0old/y0
           
           if(y0old/y0 >  3.05) then ! in order to non-overwrite the labels
              if(io == io1) then
                 write(ignu, '(a11,i1,a6,i1,a7,es12.4,a1,es12.4,a6)') &
                      'set label ', k, ' "O(h^', io, ')" at ', x0, ',',y0," right"
              else
                 write(ignu, '(a11,i1,a7,i1,a10,es12.4,a1,es12.4,a6)') &
                      'set label ', k, ' "O(h^{', io+io1, '/2})" at ', x0, ',',y0," right"
                 
              endif
           y0old = y0
           endif
        enddo
        
        write(chm,'(i2)') 7 + inum 
        
        maxlen = 40
        ch(1:maxlen) = "'orderP1.T1.L2.dat' u 2: 8 t 'P_1' w lp,"

        do k=1,Pmax
           write(chP,'(i1)') k
           chlong( (k-1)*maxlen +1 : k * maxlen) = ch(1:maxlen) 
           chlong( (k-1)*maxlen +8 : (k-1)*maxlen +8) = chP
           chlong( (k-1)*maxlen +25 : (k-1)*maxlen +26) = chm
           chlong( (k-1)*maxlen +32 : (k-1)*maxlen +33) = chP
        enddo
        chlong( 240:240) = " "

        !print*,'@@@',chlong
        
        !ch1(1:40) = "'orderP1.T1.L2.dat' u 2:8 t 'P_1' w lp,  "
        !ch2(1:40) = "'orderP2.T1.L2.dat' u 2:8 t 'P_2' w lp,  "
        !ch3(1:40) = "'orderP3.T1.L2.dat' u 2:8 t 'P_3' w lp,  "
        !ch4(1:40) = "'orderP4.T1.L2.dat' u 2:8 t 'P_4' w lp,  "
        !ch5(1:40) = "'orderP5.T1.L2.dat' u 2:8 t 'P_5' w lp,  "
        !ch6(1:40) = "'orderP6.T1.L2.dat' u 2:8 t 'P_6' w lp  "
        !write(ignu,'(a2,6a40)' ) "p ", ch1, ch2, ch3, ch4, ch5, ch6

        if(lcolor) then
           write(ignu,'(a15,a240)' ) "p[0.008:0.25] ", chlong
        else
           write(ignu,'(a10,a240)' ) "p[ : 0.25] ", chlong
        endif

        write(ignu,* ) ' '
        

     enddo  ! do norm
  enddo ! do meth

end program make_gnu

  
