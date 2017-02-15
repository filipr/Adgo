!     this program transforms seekd the maximal and minimal value
!     of the given physical quantity within the computational domain
program dgfem_range
  use plot_geom

  character*5 quantity
  character*20 ppmfile, file_type, type
  integer :: ifig
  character*20 figures(1:6)
  integer:: idat, idatx
  integer:: nelem, ndim, i, k, ideg, idof, ilev
  real :: qmin, qmax
  real :: w(30,4), q(30)
  integer:: irange, iframe

      
  idat = 10
  idatx = 9
  isol = 12
  
  
  open (idat, STATUS='old', file='movie.cfg')
  read(idat, *) quantity
  read(idat, *) irange
  read(idat, *) qmin
  read(idat, *) qmax
  read(idat, *) iframe
  read(idat, *) xmin,xmax
  read(idat, *) ymin, ymax
  read(idat, *) ppmfile
  read(idat, *) file_type
  read(idat, *) type
  read(idat,*) ifig
  do i=1, ifig
     read(idat,*) figures(i)
  enddo

  close(idat)

  !print*,'######',qmin, qmax

  open (idatx, STATUS='UNKNOWN', file='moviex.cfg')
  
  open (isol, STATUS='old', file='solx')
  
  read(isol, *) nelem, ndim

   do i=1,nelem
      do k=1,ndim
         read(isol,*) ideg, ilev,  w(1:(ideg+1)*(ideg+2)/2, k )
      enddo
      
      idof = (ideg+1)*(ideg+2)/2

      if(quantity == 'hp') then
         q(1:idof) = ideg
      else
         call ComputeQuantity(idof, ndim, w(1:idof,1:ndim), q(1:idof), quantity )
      endif

      do k=1,idof
         qmin = min(qmin, q(k) )
         qmax = max(qmax, q(k) )
      enddo
   enddo

  !print*,'######',qmin, qmax

  write(idatx, *) quantity
  write(idatx, *) irange
  write(idatx, *) qmin
  write(idatx, *) qmax
  write(idatx, *) iframe
  write(idatx, *) xmin,xmax
  write(idatx, *) ymin, ymax
  write(idatx, *) ppmfile
  write(idatx, *) file_type
  write(idatx, *) type
  write(idatx,*) ifig
  do i=1, ifig
     write(idatx,*) figures(i)
  enddo

  close(idatx)


end program dgfem_range


