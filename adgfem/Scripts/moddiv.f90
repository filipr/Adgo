program moddiv
  real :: a,b
  
  read(*,*) a, b
  
  !write(*,'(f10.5)') a/b
  !write(*,'(es12.4)') a/b
  write(*,'(i5)') mod(int(a), int(b))

end program Moddiv
