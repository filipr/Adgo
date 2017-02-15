  program sqrt3_order_dat

    use paramets

  character*150 input_file, output_file

  real, dimension (:,:), allocatable :: x, y,  z, ord
  integer, dimension(:), allocatable :: idof, ix
  integer, dimension(:,:), allocatable :: ipar
  integer:: maxmesh, ival, ivar

  if (command_argument_count() == 2) then
     call get_command_argument(1,input_file)
     call get_command_argument(2,output_file)

  else
     print*,'% Syntax:  sqrt3_order_dat "input_file" "output_file" '
     stop
  endif

  maxmesh = 500
  ivar = 4
  Nx = 10
  Ny0 =  max_errS + max_errSTnorm + max_eta  ! = 50?
  Ny = Ny0 + 2 !FR changed from 2
  Nz = 8

  !allocate( x(1:maxmesh,ival+ivar), y(maxmesh,ival+ivar),r(ival,3), ix(maxmesh) )
  print*,'::',maxmesh, Nx, Ny, Nz, Nx+Ny+Nz, max_errS , max_errSTnorm , max_eta
  allocate( x(1:maxmesh,1:Nx), y(1:maxmesh,1:Ny), ord(0:maxmesh, 1:Ny), z(1:maxmesh, 1:Nz) )
  allocate( idof(maxmesh), ipar(maxmesh, 1:5), ix(1:maxmesh) )

  ! reading of the data
  ires = 11
  open(ires, file=input_file,status='OLD')

  iout = 12
  open(iout, file=output_file,status='UNKNOWN')

  do k = 1,maxmesh
     read(ires, *, end=110) ix(k), x(k,2:3), idof(k), ipar(k, 1:5), x(k,10), y(k, 1:Ny0), z(k, 1:Nz)

     y(k,50) = 1000.
     if(y(k,2) > 0.) y(k, 50) = y(k, 36) / y(k, 2)
     
     write(iout, '(i6,2es14.6, es12.4, 3i4, 2i7, es11.3, 59es16.6, 2i7)') &
          ix(k), x(k,2:3), idof(k)**(1./3), ipar(k, 1:5), x(k,10), y(k, 1:Ny0), z(k, 1:Nz)

     !write(*,*) Nx, Ny, Nz, x(k, 1:4),y(k, 1:4)
  enddo

110 continue
  close(ires)
  close(iout)


  deallocate(x, y, ord, z, idof, ipar, ix)

end program sqrt3_order_dat
