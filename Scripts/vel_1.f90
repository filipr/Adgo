program Velocity_1
  real :: a
  real :: pi

  pi = 2 * asin(1.) 
  
  
  read(*,*) a
  
  write(*,'(es14.6)') cos(a/180* pi)


end program Velocity_1
