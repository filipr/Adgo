program Velocity_2
  real :: a
  real :: pi

  pi = 2 * asin(1.) 
  
  read(*,*) a
  
  write(*,'(es14.6)') sin(a/180 * pi)


end program Velocity_2
