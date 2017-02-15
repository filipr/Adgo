program Barenblatt
  implicit none
  integer:: N, i
  real :: x, m, r2, rt, rtt, t, u, ev, a0, a1, ev1,ev0, rf, rft,rm
  character(len=20) :: typ

  
  if (command_argument_count() == 2) then
    call get_command_argument(1,typ)
    open(12,file='smaz', status='unknown')
    write(12,*) typ
    close(12)
    open(12,file='smaz', status='OLD')
    read(12,*) m
    close(12)

    call get_command_argument(2,typ)
    open(12,file='smaz', status='unknown')
    write(12,*) typ
    close(12)
    open(12,file='smaz', status='OLD')
    read(12,*) t
    close(12)
 else
    print *,'Syntax gen_exact_Barenblatx $par $time'
    stop
 endif

 open(30,file='BBexact', status='unknown')


 rm = (m-1.)/(4*m*m)
 rt = (t+1)**(-1./m)
 rtt = (-1./m) * (t+1)**(-1./m - 1)
    
    
 N = 200

 do i=0,N
    x = -6 + 12.*i / N
    
    r2 = 2*x**2
    ev = 1. - rm * r2 * rt
    
    if(ev < 0) then
       u = 0.
    else
       rf = 1./(t+1)
       rft = -1./(t+1)**2
       
       a0 = m/(m-1)
       a1 = a0 - 1.
       !a2 = a1 - 1.
       
       ev0 = ev**a0
       ev1 = a0 * ev**a1
       !ev2 = a0 *a1 * ev**a2
       
       u = ( rf * ev0)**(1./m)
    endif
    
    write(30,*) x, u
 end do

close(30)

end program Barenblatt
