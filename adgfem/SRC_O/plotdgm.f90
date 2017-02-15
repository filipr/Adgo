!     this program transforms the outputs of <Adgfem> (<tri>, <sol>) to the
!     file <fv_tri_res> needed for <elem_vis> visualization
!     if w is polynom of degree k on K => K s devdided on (k+1)^2 subtriangles
program plotdgm
  use plot_geom

  !> file containg grid in generalised ANGENER format
  character(len=50) :: trifile
  !> file containg output solution, i.e.. basis coefficients
  character(len=50) :: solfile
  !>
  character(len=50) :: gnufile

  character(len=150) :: chline

  character*5 quantity, quantity1
  character*6 ptype
  character*20 qisol
  !character*16 time 
  real ::  time 

  integer, dimension(:, :), allocatable :: lnd
  integer, dimension(:, :), allocatable :: lbn
  real, dimension(:,:), allocatable :: x
  integer, dimension(:,:,:), allocatable :: subtri
  integer, dimension(:,:,:), allocatable :: subedge
  real, dimension(:,:,:), allocatable :: lambda

  integer, dimension(:,:), allocatable :: deg
  integer, dimension(:), allocatable :: iwall
  real, dimension(:,:), allocatable :: q, q1
  real, dimension(:,:,:), allocatable :: w
  integer:: itri, isol, idat, ignu, iw, inum
  integer:: npoin, nelem, nnelem, ndim, i1,i2,i3,i4, k,j, ideg, idof, RG_lev, j1
  integer :: max_dof, sub_tot, nisol, irange
  real :: r1,r2,r3,r4,  qmin, qmax, rarea, rlen, rratio,a1,a2
  real :: xA(1:2), xB(1:2)
  integer :: ichar

  ndim = 4

  max_dof = (max_deg+1)*(max_deg+2)/2

  itri = 10
  isol = 11
  idat = 12
  ignu = 13

  if (command_argument_count() == 2) then
    call get_command_argument(1,trifile)
    open(itri,file=trifile,action='read')

    call get_command_argument(2,solfile)
    open(isol,file=solfile,action='read')
  else
     print*,'Syntax plotdgm <tri_file>  <sol_file>'
     stop
  endif
  

  open (idat, STATUS='old', file='plot.dgm')

  read(itri, *) npoin, nelem, nbelm,i2
  read(itri, *) r1,r2, i1,i2, r3,r4, i3,i4
  
  read(isol, *) nnelem, ndim, time
  !   read(isol, *) nnelem, ndim
  if(nnelem .ne. nelem) then
     print*,' Different values of nelem in <tri> and <sol> files'
     print*,nelem,nnelem
     stop
  endif
  

  allocate(subtri(0:max_deg, 1: max_deg**2, 1:3) )
  allocate(subedge(0:max_deg, 1:3, 1: max_deg+1) )
  allocate(lambda(0:max_deg, 1: max_dof, 1:3) )
  !print*,'&',max_deg, max_dof
  call SetSubTri(max_deg, max_dof, subtri, lambda, subedge)
  !print*,'&'


  allocate(lnd(nelem,3) )
  allocate(lbn(nbelm,5) )
  allocate(x(npoin, 1:2) )
  allocate(deg(nelem, 1:2) )
  allocate(q(nelem,1:max_dof) )

  allocate(w(1:nelem, 1:max_dof, 1:ndim) )
  


  read (itri, *) ((x(k, j), j=1,2), k=1,npoin)
  read (itri, *) ((lnd(k,j), j=1,3), k=1,nelem)  
  read (itri, *) ((lbn(k,j), j=1,3), k=1,nbelm)  
  close(itri)
  
!  xmin = 1E+20
!   xmax = -1E+20
!   ymin = 1E+20
!   ymax = -1E+20
!
!   do k=1,npoin
!      xmin = min(xmin, x(k) )
!      xmax = max(xmax, x(k) )
!      ymin = min(ymin, y(k) )
!      ymax = max(ymax, y(k) )
!   enddo


!   qmin = 1E+20
!   qmax = -1E+20

  call SeekBoundElem(nelem, nbelm, lnd, lbn)

   sub_tot = 0
   do i=1,nelem
      do k=1,ndim
         !print*,i,k, nelem, ndim
         read(isol,*) ideg, RG_lev, w(i, 1:(ideg+1)*(ideg+2)/2, k )

         deg(i,1) = ideg
         deg(i,2) = RG_lev

         if(ideg .gt. max_deg) then
            print*,'Degree ', ideg, ' is not impelemnted yet !'
            stop
         endif
      enddo
 
   enddo

   if(ndim == 4) read(isol, *) rho_infty, v_infty, p_infty, alpha_infty, theta_infty

   close(isol)

   inum = 0
   call SetGnuFilename(inum, gnufile)
   
   ! mesh writing
   open(ignu, file=gnufile, status='replace')

   open(99, file='streching', status='unknown', position='append')

   rratio = 1
   do i=1,nelem
      rmax = 0.
      do j=0,3
         j1 = mod(j,3)+1
         write(ignu, *) x(lnd(i,j1 ), 1:2), deg(i,1:2), i 
         if(j>=1) rmax = max(rmax, VectorNorm( x(lnd(i,j), 1:2) - x(lnd(i,j1), 1:2) ))
      enddo
      area = x(lnd(i,1 ), 1) * ( x(lnd(i,2 ), 2) - x(lnd(i,3 ), 2) ) &
           + x(lnd(i,2 ), 1) * ( x(lnd(i,3 ), 2) - x(lnd(i,1 ), 2) ) &
           + x(lnd(i,3 ), 1) * ( x(lnd(i,1 ), 2) - x(lnd(i,2 ), 2) ) 

      rratio = max(rratio, rmax**2 / area)

      !xc(1) = sum(x(lnd(i,1:3),1)) /3
      !xc(2) = sum(x(lnd(i,1:3),2)) /3
      !if(xc(1) > 0.9999 .and. xc(2) > 0.5002 .and. xc(2) < 0.5005) &
      !     print*,xc(:), i, '@@@'
      write(ignu,*)
      write(ignu,*)
   enddo
   close(ignu)
   write(*,'(2a6 )') gnufile(1:6),': mesh'

   write(99,*) nelem, rratio
   close(99)

   !call PlotHPColorMesh(npoin, nelem, x(1:npoin, 1:2), lnd(1:nelem, 1:3), deg(1:nelem, 1:2) )


   !reading of 'plot.dgm'
   read(idat, *) i1, a1, a2  ! not used !!!!!!!!!!!!!!!!!!!!!!!

   do i1=1, 1000
      ichar = 1
      read(idat, '(a150)') chline
      ! quantity for visualization
      call ExtractChar(chline, ichar, quantity)
      if(quantity == 'EOD') goto 110

      ! opening of the output file gnu.**
      inum = inum + 1
      call SetGnuFilename(inum, gnufile)
      open(ignu, file=gnufile, status='UNKNOWN')

      ! calculating of the quantity for visualization
      do i=1,nelem
         idof = (deg(i,1) + 1)*(deg(i,1) + 2) /2
         call ComputeQuantity(idof, ndim, w(i,1:idof,1:ndim), q(i, 1:idof), quantity )
      enddo

      ! type of figure: 3D, ISO, WALL, CUT
      call ExtractChar(chline, ichar, ptype)

      
      ! drawing ...
      if(ptype == '3D' ) then
         call Plot3D(ignu, nelem, npoin, max_deg, max_dof, &
              subtri, lambda, x, lnd, deg, q)


         write(*,'(a6,a2,a6, a5 )')gnufile(1:6),': ',quantity, ptype

      elseif(ptype == '3Dfr' ) then
         call Plot3DFrame(ignu, nelem, npoin, max_deg, max_dof, &
              subtri, lambda, x, lnd, deg, q)

         write(*,'(a6,a2,a6, a5 )')gnufile(1:6),': ',quantity, ptype

      elseif(ptype == 'VEC' ) then
         allocate(q1(nelem,1:max_dof) )

         ! second component of the vector
         call ExtractChar(chline, ichar, quantity1)
         do i=1,nelem
            idof = (deg(i,1) + 1)*(deg(i,1) + 2) /2
            call ComputeQuantity(idof, ndim, w(i,1:idof,1:ndim), q1(i, 1:idof), quantity1 )
         enddo


         call PlotVectors(ignu, nelem, npoin, max_deg, max_dof, &
              subtri, lambda, x, lnd, deg, q, q1)

         write(*,'(a6,a2,2a6, a5 )')gnufile(1:6),': ',quantity, quantity1, ptype
         deallocate(q1)

      elseif(ptype == 'ISO' ) then
         ! number of isolines
         call ExtractChar(chline, ichar, qisol)
         read(qisol,*) nisol

         ! range: automatic (irange = 0) or user given (irange /= 0)
         call ExtractChar(chline, ichar, qisol)
         read(qisol,*) irange

         if(irange /= 0) then
            ! range :: qmin
            call ExtractChar(chline, ichar, qisol)
            read(qisol,*) qmin
            
            ! range :: qmax
            call ExtractChar(chline, ichar, qisol)
            read(qisol,*) qmax
         endif

         call PlotIsolines(ignu, nelem, npoin, nbelm, max_deg, max_dof, &
              subtri, lambda, x, lnd, lbn, deg, q, nisol, irange, qmin, qmax, r1,r2,r3,r4)

         write(*,'(a6,a2,a6, a5, a9, i5,a13,es16.8,a1,es16.8,a2 )') &
              gnufile(1:6),': ',quantity, ptype,'num_iso =',nisol, &
              ',  range = ( ',qmin,",", qmax," )"

         open(22, file='plot_min_max', status='UNKNOWN', position='append')
         write(22,*), qmin, qmax, quantity
         close(22)

      elseif(ptype == 'CUT' ) then
         ! a cut along the line xA(1:2) <--> xB(1:2)
         call ExtractChar(chline, ichar, qisol)
         read(qisol,*) xA(1)
         call ExtractChar(chline, ichar, qisol)
         read(qisol,*) xA(2)

         call ExtractChar(chline, ichar, qisol)
         read(qisol,*) xB(1)
         call ExtractChar(chline, ichar, qisol)
         read(qisol,*) xB(2)
            
         call PlotCUT(ignu, nelem, npoin, max_deg, max_dof, &
              subtri, lambda, x, lnd, deg, q, xA, xB)

         write(*,'(a6,a2,a6, a5, a1, es12.4,a1,es12.4,a8,es12.4,a1,es12.4,a1 )')  &
              gnufile(1:6),': ',quantity, ptype, &
              "(",xA(1),",",xA(2),') <--> (',xB(1),",",xB(2),')'

      elseif(ptype == 'WALL' ) then
         call ExtractChar(chline, ichar, qisol)
         read(qisol,*) iw
         allocate(iwall(iw) )
         do i=1, iw
            call ExtractChar(chline, ichar, qisol)
            read(qisol,*) iwall(i)
         enddo

         call PlotWalls(ignu, nelem, npoin, nbelm, max_deg, max_dof, &
              lambda, subedge, x, lnd, lbn, deg, q, iw, iwall)
         
         write(*,'(a6,a2,a6, a5,a12,30 i5 )')gnufile(1:6),': ',quantity, ptype, &
              'bound. comp:',iwall(1:iw)
         
         deallocate(iwall)
      else
         print*,'Drawing method ',ptype,' is not implemented'
      endif
      close(ignu)

   enddo

110 continue


   deallocate(lnd,x,deg,q)
   
 end program plotdgm

 subroutine ExtractChar(chline, ichar, chname)
   character(len=150), intent(in) :: chline
   character(len=*), intent(inout) :: chname
   integer, intent(inout) :: ichar 

   integer :: i, j, is

   chname(1:len(chname)) = ' '

   is = ichar
   do i=is, len(chline)
      if(chline(i:i) /= ' ') goto 20
      ichar = ichar + 1
   enddo
   
20 continue
   

   is = ichar
   j = 0
   do i=is, len(chline)
      if(chline(i:i) /= ' ') then
         j = j + 1
         chname(j:j) = chline(i:i)

         !print*,'@@',i,j,chname(j:j), chline(i:i)
      else
         goto 10
      endif
      ichar = ichar + 1
   enddo
   
10 continue

   !print*,'next char =', ichar,'name =', chname(1:j)

 end subroutine ExtractChar


 subroutine SeekBoundElem(nelem, nbelm, lnd, lbn)
   integer, intent(in) :: nelem, nbelm
   integer, dimension(1:nelem, 1:3), intent(in)  :: lnd
   integer, dimension(1:nbelm, 1:5), intent(inout)  :: lbn

   integer :: i,ie, j, j1, k1, k2

   do i=1,nbelm
      k1 = lbn(i,1)
      k2 = lbn(i,2)

      do ie = 1,nelem
         do j=1,3
            j1 = mod(j, 3) + 1
            if(k1 == lnd(ie, j)  .and. k2 == lnd(ie, j1)) then
               lbn(i,4) = ie
               lbn(i,5) = j
               goto 10
            endif
         enddo
      enddo
10    continue

!      write(*,'(a2,6i5)') '!!!',i,lbn(i,1:5)

   enddo
 end subroutine SeekBoundElem
