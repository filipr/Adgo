program ExtrOrder
  integer:: maxmesh, ival, ivar     
  real, dimension (:,:), allocatable :: x, y, r, orders
  integer, dimension(:), allocatable :: ix, idof
  integer, dimension(:,:), allocatable :: ipar
  logical :: logEOC
  integer :: ife, ife_b, ife_e, ife_txt, iEOC, iline
  character*50 text(1:30)
  character*7  EOC
  character*10 typ

  ! ival = number of real values in one row in order.dat at the input 
  ! ncounts = number of real columns in tables at the output
  logEOC  = .false.

  iEOC = 3

  iline = 7  ! after "iline" lines the  \hline is plotted
  
  
  
  maxmesh = 500 
  ivar = 4
  ival = 30
  allocate( x(maxmesh,ival+ivar), y(maxmesh,ival+ivar),r(ival,3), ix(maxmesh) )
  allocate( idof(maxmesh), ipar(maxmesh,3), orders(0:maxmesh, ival) )
  
  
  imesh0 =1

  ! reading of the data
  ires = 11
  open(ires, file='order.dat',status='OLD')

  ife_txt = 23
  open (ife_txt,file='tab.txt', status='UNKNOWN')
  
  do k = 1,maxmesh
     read(ires, *, end=110) ix(k),x(k,1), x(k,3), idof(k), &
          ipar(k,1:3), (x(k,j),j=8, 27)
     
     
     write(ife_txt,*) ix(k),x(k,1), x(k,3), idof(k), &
          x(k,  8:11), &  ! 5..8   L2, H1, X, J at last time
          x(k, 12:13), &  ! 9..10  space time X, Y norms
          x(k, 14:17), &  ! 11..14 etas A, S, T, ST  for space time
          x(k, 18:21), &  ! 15..18 etas A, S, T, ST  at last time level
          x(k, 16)/ x(k, 15) ! 19  ratio eta^T/ eta^S

  enddo
110 imesh = k-1
  print*,'Number of records (lines) =',imesh
  close(ires)
  

  ! average step mesh h
!  x(1:imesh,2) = 12./((ix(1:imesh)/2)**0.5)   ! computationa domain size 12 x 12
!  x(1:imesh,2) = 2./((ix(1:imesh)/2)**0.5)   ! computationa domain size 2 x 2
!  x(1:imesh,2) = 1./((ix(1:imesh)/2)**0.5)   ! computationa domain size 1 x 1

  

end program ExtrOrder
