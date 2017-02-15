program DivideNumbers
  real :: a,b
  
  read(*,*) a, b
  
  !write(*,'(f10.5)') a/b
  !write(*,'(es12.4)') a/b
  write(*,'(es9.2)') min( a/b, 0.1)
  !write(*,'(es12.5)') a/b

end program DivideNumbers
