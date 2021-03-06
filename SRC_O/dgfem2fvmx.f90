!     this program transforms the outputs of <Adgfem> (<tri>, <sol>) to the
!     file <fv_tri_res> needed for <elem_vis> visualization
!     if w is polynom of degree k on K => K s devdided on (k+1)^2 subtriangles
program dgfem2fvm
  use plot_geom

  character*5 quantity
  !character*16 time 
  real ::  time 

  integer, dimension(:, :), allocatable :: lnd
  real, dimension(:), allocatable :: x, y
  integer, dimension(:,:,:), allocatable :: subtri
  integer, dimension(:,:,:), allocatable :: subedge
  real, dimension(:,:,:), allocatable :: lambda

  integer, dimension(:), allocatable :: deg
  real, dimension(:,:), allocatable :: q

  real, dimension(:,:), allocatable :: rgb
  real, dimension(:,:), allocatable :: w

  integer:: itri, isol, ires, ipal, ipoc
  integer:: npoin, nelem, nnelem, ndim, i1,i2,i3,i4, k,j, ideg, idof
  integer :: max_dof, sub_tot, RG_level
  real :: xmin, xmax, ymin, ymax
  real :: r1,r2,r3,r4, lambda_c(1:3), qmin, qmax
  real ::  xl(3), yl(3), val
  integer :: il(3)
  
  !! max_deg  initiated and defined in plgeom.f90
  
  max_dof = (max_deg+1)*(max_deg+2)/2

      
  idat = 10
  itri = 11
  isol = 12
  ires = 13
  
  
  open (idat, STATUS='old', file='movie.cfg')
  read(idat, *) quantity
  close(idat)
  
  open (itri, STATUS='old', file='tri')
  open (isol, STATUS='old', file='solx')
  
  open (ires, STATUS='unknown', file='fv_tri_res')
  


   read(itri, *) npoin, nelem, i1,i2
   read(itri, *) r1,r2, i1,i2, r3,r4, i3,i4

   read(isol, *) nnelem, ndim, time
!   read(isol, *) nnelem, ndim
   if(nnelem .ne. nelem) then
      print*,' Different values of nelem in <tri> and <sol> files'
      print*,nelem,nnelem
      stop
   endif

  allocate(w(1:max_dof, 1:ndim) )


   allocate(subtri(0:max_deg, 1: max_deg**2, 1:3) )
   allocate(subedge(0:max_deg, 1:3, 1: max_deg + 1) )
   allocate(lambda(0:max_deg, 1: max_dof, 1:3) )
   !print*,'&',max_deg, max_dof
   call SetSubTri(max_deg, max_dof, subtri, lambda, subedge)
   !print*,'&'

   allocate(lnd(nelem,3) )
   allocate(x(npoin), y(npoin) )
   allocate(deg(nelem) )
   allocate(q(1:nelem, 1:max_dof) )

   read (itri, *) (x(k), y(k), k=1,npoin)    
   read (itri, *) ((lnd(k,j), j=1,3), k=1,nelem)  
   close(itri)

   xmin = 1E+20
   xmax = -1E+20
   ymin = 1E+20
   ymax = -1E+20

   do k=1,npoin
      xmin = min(xmin, x(k) )
      xmax = max(xmax, x(k) )
      ymin = min(ymin, y(k) )
      ymax = max(ymax, y(k) )
   enddo


   qmin = 1E+20
   qmax = -1E+20

   sub_tot = 0
   do i=1,nelem
      do k=1,ndim
         read(isol,*) ideg, RG_level, w(1:(ideg+1)*(ideg+2)/2, k )
         !read(isol,*) ideg, w(1:(ideg+1)*(ideg+2)/2, k )  ! OLD version

         !!!if(i < 3) write(*,*)'!!!', ideg, RG_level, w(1:(ideg+1)*(ideg+2)/2, k)
         if(ideg .gt. max_deg) then
            print*,'Degree ', ideg, ' is not impelemnted yet !'
            stop
         endif
      enddo
      
      if(quantity == 'hp') then
         deg(i) = 0
         idof = 1
         sub_tot = sub_tot + 1

         q(i,1) = ideg
      else
         deg(i) = ideg
         idof = (ideg+1)*(ideg+2)/2
         sub_tot = sub_tot + (max(1,ideg))**2
         
         call ComputeQuantity(idof, ndim, w(1:idof,1:ndim), q(i, 1:idof), quantity )

        !  if( i < 3) then
        !    do k=1,ndim
        !        write(*,'(a6,3i5,50es12.4)') '!@#$',i, k,idof, w(k,1:ndim)
        !     enddo
        !    write(*,'(a6,3i5,50es12.4)') '!@#$',i, 0,idof, q(i,1:idof)
        ! endif
      endif
      
      do k=1,idof
         qmin = min(qmin, q(i,k) )
         qmax = max(qmax, q(i,k) )
      enddo
   enddo

   !print*,'DGM2FVM', sub_tot, xmin, xmax, ymin, ymax, qmin, qmax, time

   write(ires,*) sub_tot, xmin, xmax, ymin, ymax, qmin, qmax, time
   
   ipoc = 0
   do i=1,nelem
      ideg = max(1, deg(i))
      idof = (ideg+1)*(ideg+2)/2

      do j=1,(ideg)**2
         il(1:3) = subtri(ideg, j, 1:3)
         
         lambda_c(1:3) = 0. ! barycentric coordinates of the centre of subtriangle
         
         do l=1,3
            xl(l) = dot_product(lambda(ideg, il(l), 1:3), x(lnd(i, 1:3)) )
            yl(l) = dot_product(lambda(ideg, il(l), 1:3), y(lnd(i, 1:3)) )
            !lambda_c(1:3) = lambda_c(1:3) + lambda(ideg+1, il(l), 1:3)
            !           write(20+j,*)xl(l),yl(l)
         enddo
         
        !        write(20+j,*)xl(1),yl(1)
        
        
        !lambda_c(1:3) = lambda_c(1:3) / 3
        
         if(deg(i) == 0) then
            val = q(i, 1) 
         else
            val = sum(q(i,il(1:3)) )/3
         endif

        !val = SetPolynomValue(ideg, idof, q(i,1:idof), lambda_c(1:3) )

         !if( i < 10) then
         !   print*,'!!!!dgfem2fvmx' ,i,val
         !endif

        ipoc = ipoc + 1
        write(ires,'(7es14.6,a5,2i6)') xl(1), yl(1), xl(2), yl(2), xl(3), yl(3), val, &
             '  ##=',ipoc,i


     enddo
!     stop

  enddo

  write(ires,*) nelem
  do i=1,nelem
     write(ires,'(6es14.6, i5)') &
          x(lnd(i,1)),y(lnd(i,1)),x(lnd(i,2)),y(lnd(i,2)),x(lnd(i,3)),y(lnd(i,3)),&
          deg(i)
  enddo
     

  close(ires)
  deallocate(lnd,x,deg,q)

end program dgfem2fvm




subroutine SetSubTri1(max_deg, max_dof, subtri, lambda)
  integer, intent(in) :: max_deg, max_dof
  integer,dimension(1:max_deg, 1: max_deg**2, 1:3), intent(inout) :: subtri
  real,dimension(1:max_deg, 1: max_dof, 1:3), intent(inout) :: lambda
  
  if(max_deg .ge. 1) then
     ! P_1 approximation
     subtri(1, 1, 1:3) = (/1, 2, 3 /)
     lambda(1, 1, 1:3) = (/  1.,  0.,  0. /)
     lambda(1, 2, 1:3) = (/  0.,  1.,  0. /)
     lambda(1, 3, 1:3) = (/  0.,  0.,  1. /)

     if(max_deg .ge. 2) then
        ! P_2 approximation
        subtri(2, 1, 1:3) = (/ 1, 4, 6 /)
        subtri(2, 2, 1:3) = (/ 2, 5, 4 /)
        subtri(2, 3, 1:3) = (/ 3, 6, 5 /)
        subtri(2, 4, 1:3) = (/ 4, 5, 6 /)

        lambda(2, 1, 1:3) = (/  1.,   0.,    0. /)
        lambda(2, 2, 1:3) = (/  0.,   1.,    0. /)
        lambda(2, 3, 1:3) = (/  0.,   0.,    1. /)
        lambda(2, 4, 1:3) = (/  1./2, 1./2,  0. /)
        lambda(2, 5, 1:3) = (/  0.,   1./2,  1./2 /)
        lambda(2, 6, 1:3) = (/  1./2, 0.,    1./2 /)

        if(max_deg .ge. 3) then
           ! P_3 approximation
           subtri(3, 1, 1:3) = (/ 1, 4, 9 /)
           subtri(3, 2, 1:3) = (/ 2, 6, 5 /)
           subtri(3, 3, 1:3) = (/ 3, 8, 7 /)
           subtri(3, 4, 1:3) = (/ 4, 5, 10 /)
           subtri(3, 5, 1:3) = (/ 6, 7, 10 /)
           subtri(3, 6, 1:3) = (/ 8, 9, 10 /)
           subtri(3, 7, 1:3) = (/ 4, 10, 9 /)
           subtri(3, 8, 1:3) = (/ 6, 10, 5 /)
           subtri(3, 9, 1:3) = (/ 8, 10, 7 /)

           lambda(3, 1, 1:3) = (/  1.,   0.,    0. /)
           lambda(3, 2, 1:3) = (/  0.,   1.,    0. /)
           lambda(3, 3, 1:3) = (/  0.,   0.,    1. /)
           lambda(3, 4, 1:3) = (/  2./3, 1./3,  0. /)
           lambda(3, 5, 1:3) = (/  1./3, 2./3,  0. /)
           lambda(3, 6, 1:3) = (/  0.,   2./3,  1./3 /)
           lambda(3, 7, 1:3) = (/  0.,   1./3,  2./3 /)
           lambda(3, 8, 1:3) = (/  1./3, 0. ,   2./3 /)
           lambda(3, 9, 1:3) = (/  2./3, 0. ,   1./3 /)
           lambda(3,10, 1:3) = (/  1./3, 1./3,    1./3 /)


           if(max_deg .ge. 4) then
              ! P_4 approximation
              subtri(4,  1, 1:3) = (/  1,  4, 12  /)
              subtri(4,  2, 1:3) = (/  4,  5, 13  /)
              subtri(4,  3, 1:3) = (/  5,  6, 14  /)
              subtri(4,  4, 1:3) = (/  6,  2,  7  /)
              subtri(4,  5, 1:3) = (/ 12,  4, 13  /)
              subtri(4,  6, 1:3) = (/ 13,  5, 14  /)
              subtri(4,  7, 1:3) = (/ 14,  6,  7  /)
              subtri(4,  8, 1:3) = (/ 12, 13, 11  /)
              subtri(4,  9, 1:3) = (/ 13, 14, 15  /)
              subtri(4, 10, 1:3) = (/  8, 14,  7  /)
              subtri(4, 11, 1:3) = (/ 11, 13, 15  /)
              subtri(4, 12, 1:3) = (/ 15, 14,  8  /)
              subtri(4, 13, 1:3) = (/ 10, 11, 15  /)
              subtri(4, 14, 1:3) = (/  9, 15,  8  /)
              subtri(4, 15, 1:3) = (/ 10, 15,  9  /)
              subtri(4, 16, 1:3) = (/  3, 10,  9  /)
              
              lambda(4,  1, 1:3) = (/  1.  , 0.    , 0.      /)
              lambda(4,  2, 1:3) = (/  0.  , 1.    , 0.      /)
              lambda(4,  3, 1:3) = (/  0.  , 0.    , 1.      /)
              lambda(4,  4, 1:3) = (/  3./4, 1./4  , 0.      /)
              lambda(4,  5, 1:3) = (/  1./2, 1./2  , 0.      /)
              lambda(4,  6, 1:3) = (/  1./4, 3./4  , 0.      /)
              lambda(4,  7, 1:3) = (/  0.  , 3./4  , 1./4    /)
              lambda(4,  8, 1:3) = (/  0.  , 1./2  , 1./2    /)
              lambda(4,  9, 1:3) = (/  0.  , 1./4  , 3./4    /)
              lambda(4, 10, 1:3) = (/  1./4, 0.    , 3./4    /)
              lambda(4, 11, 1:3) = (/  1./2, 0.    , 1./2    /)
              lambda(4, 12, 1:3) = (/  3./4, 0.    , 1./4    /)
              lambda(4, 13, 1:3) = (/  1./2, 1./4  , 1./4    /)
              lambda(4, 14, 1:3) = (/  1./4, 1./2  , 1./4    /)
              lambda(4, 15, 1:3) = (/  1./4, 1./4  , 1./2    /)


              if(max_deg .ge. 5) then
                 print*,'Not impelemnted'
                 stop
              endif
           endif
        endif
     endif
  endif

end subroutine SetSubTri1


function  SetPolynomValue(ideg, idof, q, lambda )
  real :: SetPolynomValue
  integer, intent(in) :: idof, ideg
  real, dimension(1:idof), intent(in) :: q
  real, dimension(1:3), intent(in) :: lambda

  if(ideg == 0) then
     SetPolynomValue = q(1)
  elseif(ideg == 1) then
     SetPolynomValue = dot_product(q(1:3), lambda(1:3))
  elseif(ideg == 2) then
     SetPolynomValue = &
          + q(1) * lambda(1) * (2*lambda(1) - 1.) &
          + q(2) * lambda(2) * (2*lambda(2) - 1.) &
          + q(3) * lambda(3) * (2*lambda(3) - 1.) &
          + q(4) * 4.* lambda(1) * lambda(2) &
          + q(5) * 4.* lambda(2) * lambda(3) &
          + q(6) * 4.* lambda(3) * lambda(1) 
  elseif(ideg == 3) then
     SetPolynomValue = &
          + q(1) * lambda(1) * (3*lambda(1) - 1.) * (3*lambda(1) - 2.) /2. &
          + q(2) * lambda(2) * (3*lambda(2) - 1.) * (3*lambda(2) - 2.) /2. &
          + q(3) * lambda(3) * (3*lambda(3) - 1.) * (3*lambda(3) - 2.) /2. &
          + q(4) * 4.5 * lambda(1) * lambda(2) * (3*lambda(1) -1)  &
          + q(5) * 4.5 * lambda(1) * lambda(2) * (3*lambda(2) -1)  &
          + q(6) * 4.5 * lambda(2) * lambda(3) * (3*lambda(2) -1)  &
          + q(7) * 4.5 * lambda(2) * lambda(3) * (3*lambda(3) -1)  &
          + q(8) * 4.5 * lambda(3) * lambda(1) * (3*lambda(3) -1)  &
          + q(9) * 4.5 * lambda(3) * lambda(1) * (3*lambda(1) -1)  &
          + q(10) * 27.* lambda(1) * lambda(2) * lambda(3) 
  else
     print*,'Not YET ImpleMented !'
     stop
  endif
end function SetPolynomValue
