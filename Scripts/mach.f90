program Mach_number
  real :: M
  real :: kappa 

  kappa = 1.4
  
  
  read(*,*) M
  
  write(*,'(es14.6)') 1. / (kappa * M * M )


end program Mach_number

! M^2 = 1./(kappa * p)
! M^2 * p * kappa = 1
! p = 1./(kappa * M^2)
