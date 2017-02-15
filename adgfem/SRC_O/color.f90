
module color_figures
  use geometry
  use hp_adaptation

  real :: xmin, xmax, ymin, ymax
  integer:: num_color
  real, allocatable, dimension(:,:) :: rgb

  public:: InitColor
  public:: PlotMesh_hpColor
contains


  subroutine InitColor(grid)
    type(mesh), intent(in) :: grid
    integer :: ipal, i

    ! reading of colors
    ipal = 46
    open(ipal, file='paleta',status='unknown')
    read(ipal,*) num_color
    allocate(rgb(1:num_color, 1:3))

    print*,' # Number of colors = ',num_color
    do i=1,num_color
       read(ipal,*) rgb(i,1:3)
    enddo

    close(ipal)

    xmin = 1E+20
    xmax = -1E+20
    ymin = 1E+20
    ymax = -1E+20

    do k=1,grid%npoin
       xmin = min(xmin, grid%x(k,1) )
       xmax = max(xmax, grid%x(k,1) )
       ymin = min(ymin, grid%x(k,2) )
       ymax = max(ymax, grid%x(k,2) )
    enddo

    write(*,'(a19,es11.4,a1,es11.4,a5,es11.4,a1,es11.4,a1)') &
         ' # Maximal frame: [',xmin,',',xmax,'] x [',ymin,',',ymax,']'
    write(*,*)
  end subroutine InitColor


  !> plotting of hp mesh in color
  subroutine PlotMesh_hpColor(grid, isol)
    type(mesh), intent(in) :: grid
    class(element), pointer :: elem
    integer, intent(in) :: isol
    real, allocatable, dimension(:,:,:) :: xt
    integer, allocatable, dimension(:,:) :: flen

    !INTEGER :: PGOPEN, PGBEG
    !external :: PGENV, PGPAP, PGSCR, PGSCI, PGPOLY, PGPAGE, PGEND  ! from PGPLOT
    character(len=15) :: solfile
    character(len=1) :: ch1
    character(len=2) :: ch2
    character(len=3) :: ch3
    character(len=4) :: ch4
    character(len=5) :: ch5
    real :: rmin, rmax
    integer :: IER, i,j, ic

    solfile = 'hpmesh00000.ps'

    if(isol >= 0 .and. isol <= 9) then
       write( ch1, '(i1)' ) isol
       solfile(11:11) = ch1

    elseif(isol >= 10 .and. isol <= 99) then
       write( ch2, '(i2)' ) isol
       solfile(10:11) = ch2

    elseif(isol >= 100 .and. isol <= 999) then
       write( ch3, '(i3)' ) isol
       solfile(9:11) = ch3

    elseif(isol >= 1000 .and. isol <= 9999) then
       write( ch4, '(i4)' ) isol
       solfile(8:11) = ch4

    elseif(isol >= 10000 .and. isol <= 99999) then
       write( ch5, '(i5)' ) isol
       solfile(7:11) = ch5

    endif

    rmin = 0.
    rmin = 10.

    print*,'!!!  A', solfile

    allocate(xt(1:grid%nelem, maxval(grid%elem(:)%flen), 1:2))
    allocate(flen(1:grid%nelem, 1:2))

    do i=1,grid%nelem
       elem => grid%elem(i)
       flen(i,1) = elem%flen
       flen(i,2) = elem%deg
       do j=1,elem%flen
          xt(i,j,1:2) = grid%x(elem%face(idx, j),1:2)
       enddo
    enddo




    call Color_draw(solfile, num_color, rgb, xmin, xmax, ymin, ymax, &
         nelem, maxval(grid%elem(:)%flen), xt, flen)


  end subroutine PlotMesh_hpColor


end module color_figures
