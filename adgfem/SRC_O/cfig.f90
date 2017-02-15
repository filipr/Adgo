!> graphical library for plotting of color maps of the solution
!> has to be used in combination with script cfigx
program cfig
!  implicit none

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                  C
!C  DRIVER EXAMPLE FOR THE GGG PICTURE ENVIRONMENT  C
!C                                                  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  real, dimension(1:1000) :: XDES,YDES
!C ------- memory for lines and fillings -------

  COMMON/SIZES/DIM(4)
!  DATA DIM/4.,4.,1.,1./   ! scalar equation on square Bigger fonts
  DATA DIM/6.,6.,1.,1./   ! scalar equation on square
!  DATA DIM/9.,9.,1.,1./   ! scalar equation on square
!  DATA DIM/9.,3.,1.,1./   ! scalar equation on square
!  DATA DIM/12.,6.,1.,1./  ! NACA - Mital
  !C --------- Size of the pictures in cm and distances ------


  character*16 axes_type, file_type, type, time
  character*20 qval, quantity, quantity_name
  integer ival, in_fig, jn_fig, idig
  integer  iscale, nelem, ncolors
  real  xmin, xmax, ymin, ymax, rmin, rmax, rtime
  real  xmin1, xmax1, ymin1, ymax1
  real  xxmin, xxmax, yymin, yymax
  logical  with_deg, color
  integer ifig, ivalue, icfg, ifig2
  character*20 figures(1:6)
  real, dimension(:,:), allocatable :: wi
  real, dimension(:,:), allocatable :: rwi
  real, dimension(:,:), allocatable :: rgb
  integer :: rnelem, mancol, manbnd, ipal, i
  real :: pmin, pmax
  real :: colmax, colmin
  real :: XINCR, YINCR
  integer :: KX, KY, IDIGIT

!c     soubor s resenim v uzlovych bodech
      ivalue = 48
!c     ... konfiguracni soubor
      icfg = 49

!c     ... openings
      open( ivalue, file='fv_tri_res', status = 'OLD')
      open( icfg, file='movie.cfg', status='unknown')

!c     prectem prvni cast konfiguracniho s.
      read(icfg,*) quantity

      read(icfg,*) mancol
      if (mancol .eq. 1) then   !    print *, 'Manualni nastaveni barev'
         read(icfg,*) colmin
         read(icfg,*) colmax
      else                      !    print *, 'Automaticke nastaveni barev'
         read(icfg,*) colmin
         read(icfg,*) colmax
      endif
      read(icfg,*) manbnd
      if (manbnd .eq. 1) then   !print *, 'Manualni nastaveni vyrezu'
         read(icfg,*) xxmin, xxmax
         read(icfg,*) yymin, yymax
      else                      !print *, 'Automaticke nastaveni vyrezu: cela oblast'
         read(icfg,*) xxmin, xxmax
         read(icfg,*) yymin, yymax
      endif
      read(icfg,*) axes_type
      read(icfg,*) file_type  ! color/ BW
      read(icfg,*) type
      read(icfg,*) ifig
      do i=1,abs(ifig)
         read(icfg,*) figures(i)
      enddo
      close(icfg)

      !print*,'@@',xxmax, xxmin, yymax, yymin
      !print*,'#####################',DIM(:)
      val = DIM(1)
      DIM(1) = val* (xxmax - xxmin )/ (yymax - yymin)
      DIM(2) = val

      if(max(DIM(1) , DIM(2) ) > 10.) then
         DIM(1) = DIM(1) / 2
         DIM(2) = DIM(2) / 2
      endif

      if(type == 'battery') then
         DIM(1) = DIM(1) * 2.
         DIM(2) = DIM(2) * 2.
      endif

      !print*,'#####################',DIM(:)

      if(file_type == 'color' .or. file_type == 'COLOR' .or. &
           file_type == 'RGB' .or. file_type == 'rgb') then
         color = .true.
      else
         color = .false.
      endif

!c      print *,'Obrazek bude v souboru  ', file_name

!c     setting of the figures for visualization

!C     nacteni palety
      ipal = 46
      open(ipal, file='paleta',status='unknown')
      read(ipal,*) ncolors
      allocate(rgb(1:ncolors, 1:3) )

      if(ncolors .gt. 512) then
         print*,'Dimension of rgb is too small !'
         stop
      endif
      !print *,'Number of colors = ',ncolors
      do i=1,ncolors
         read(ipal,*) rgb(i,1),rgb(i,2),rgb(i,3)
      enddo
      close(ipal)

      !print *, 'Soubor paleta precten'

      read(ivalue,*) nelem, xmin, xmax, ymin, ymax, rmin, rmax, rtime

      ! reading of the data
      rmax = -10000.
      rmin = 1000.

      ival = 7
      allocate(wi(1:nelem, 1:ival) )
      do i=1,nelem
         read(ivalue,*) wi (i, 1:ival)
         rmin = min(rmin, wi(i,7) )
         rmax = max(rmax, wi(i,7) )
         !if( i < 10) then
         !print*,'@@',i,wi(i,7), rmin, rmax
         !endif

      enddo

      !print*,'#######',rmin, rmax
      if(quantity(1:3) == 'Alg') then
         rmin = 0.
         idig = int(- log(rmax) /log(10.) +1.) +1
         !print*, '@@@',idig, (rmax * 10**idig ), int(rmax * 10**idig + 1), &
         !     int(rmax * 10**idig + 1) * 10.** (-idig)

         rmax = int(rmax * 10**idig + 1)
         rmax = rmax * 10.** (-idig)

      endif

      if(abs(rmax - rmin) / (abs(rmax) + abs(rmin) + 0.1) < 1E-10) then
         if(rmax > 0 ) then
            rmax = rmax *1.1
            rmin = rmin *0.9
         else
            rmax = rmax *0.9
            rmin = rmin *1.1
         endif
      endif


      !print*,'#######',rmin, rmax

      read(ivalue,*) rnelem

      ! reading of the data of hp-mesh
      pmax = 0.
      pmin = 1000.
      allocate(rwi(1:rnelem, 1:ival) )
      do i=1,rnelem
         read(ivalue,*) rwi (i, 1:ival)
         !pmin = min(pmin, int(rwi(i,7)+0.2 ) )
         !pmax = max(pmax, int(rwi(i,7)+0.2 ) )
         pmin = min(pmin, rwi(i,7) )
         pmax = max(pmax, rwi(i,7) )
         !print*,'@@@@@', pmin,pmax,rwi(i,7)
      enddo
      if(pmin == pmax) then
         pmin = pmin - 1
         pmax = pmax + 1
      endif

      !pmin = 1
      !pmax = 8

      !pmin = max(0., pmin -1)
      !pmax = pmax + 1


!c     V pripade man. barev prenastavime rmin, rmax
      if (mancol .eq. 1) then
         rmin = colmin
         rmax = colmax

         pmin = colmin
         pmax = colmax
      endif

      !print*,'#######',rmin, rmax

!C     zvoleni pozadovaneho vyrezu
      if (manbnd .eq. 1) then
         xmax1 = xxmax
         xmin1 = xxmin
         ymax1 = yymax
         ymin1 = yymin
      else
         xmax1 = xmax
         xmin1 = xmin
         ymax1 = ymax
         ymin1 = ymin
      endif


!c      print *,'W - range:',rmin,rmax


      CALL BEGIN_GGG('Gfigure')
!C ------ INITIALISATION GRAPHICS --------
!C  The statement CALL BEGIN_GGG('name')  initiates
!C  the environment, which will produce a ps-file
!C  'name.ps' and a LaTeX input file 'name.inp'
!C  to be inserted as follows in the LaTeX file:
!C  \begin{figure}[h]
!C    \centering
!C     \input name.inp
!C    \caption{Test for GGGnugnat}\label{Gnugnat}
!C  \end{figure}

      XMIN=0.
      XMAX= 1.
      YMIN=0.
      YMAX=1.


      XMIN = xmin1
      XMAX = xmax1
      YMIN = ymin1
      YMAX = ymax1

      !write(*,*) '---------- Cfig library,  X, Y, R, P ranges: '
      write(*,200) 'Cfig: X,Y,R,P range:',xmin1,xmax1,ymin1,ymax1,rmin,rmax,pmin,pmax
200   format(a20, 4('(',es8.1,':',es8.1,'), ') )

!C ----- LOCAL COORDINATE LIMITS ----

!      CALL MICKEY_MOUSE(1,1,XMIN,XMAX,YMIN,YMAX)

!C ---- place picture on paper ------
!C  CALL MICKEY_MOUSE(I,J,XMIN,XMAX,YMIN,YMAX) places a picture
!C  of size SIZE(1),SIZE(2) cm on the paper like in a Mickey Mouse
!C  comics. Can be used several times. I,J are the positions
!C  like the indices of a matrix.
!C  XMIN,XMAX,YMIN,YMAX are the bounds of the local coordinates.
!C  SIZE(3),SIZE(4) are the distances between individual pictures.

      in_fig = 0
      jn_fig = 0

      do i=1,abs(ifig)
         ! one column
         in_fig = i
         jn_fig = 1

         ! one row
         !in_fig = 1
         !jn_fig = jn_fig + 1

         if(figures(i) == 'sol' .or. figures(i) == 'sol_VS' .or. &
              figures(i) == 'solM' .or. figures(i) == 'solM_VS' .or. &
              figures(i) == 'sol_HS' .or. figures(i) == 'solM_HS' ) then
            CALL MICKEY_MOUSE(in_fig, jn_fig, XMIN,XMAX,YMIN,YMAX)
            CALL FRAME

            ! grawing the results
            call PlotColorMaps(XMIN,XMAX,YMIN,YMAX, &
                 nelem, ival, wi, ncolors, rgb, rmin, rmax, quantity, .false., color )

            ! mesh in the same figure
            if(figures(i) == 'solM' .or. figures(i) == 'solM_VS' ) &
                 call PlotGrid(rnelem, ival, rwi )

            if(axes_type /= 'N' .and. axes_type /= 'n') &
                 call DrawAxes(XMIN, XMAX, YMIN, YMAX, axes_type )

            if(figures(i) == 'sol_VS' .or. figures(i) == 'solM_VS' ) then

               jn_fig = jn_fig + 1

               CALL MICKEY_MOUSE(in_fig,jn_fig,XMIN,XMAX,YMIN,YMAX)
               !      call frame

               call PlotVscale(XMIN,XMAX,YMIN,YMAX,ncolors, rgb, rmin, rmax, &
                    quantity, .false., color, rtime, type )
            endif

            ! horizontal scale
            if(figures(i) == 'sol_HS' .or. figures(i) == 'solM_HS' ) then

               in_fig = in_fig + 1

               CALL MICKEY_MOUSE(in_fig,jn_fig,XMIN,XMAX,YMIN,YMAX)
               !      call frame

               call PlotHscale(XMIN,XMAX,YMIN,YMAX,ncolors, rgb, rmin, rmax, &
                    quantity, .false., color, rtime, type )
            endif


         elseif(figures(i) == 'grid' .or. figures(i) == 'grid_VS' .or. &
              figures(i) == 'gridHP' .or. figures(i) == 'gridHP_VS' .or. &
              figures(i) == 'grid_HS' .or. figures(i) == 'gridHP_HS' ) then
            CALL MICKEY_MOUSE(in_fig, jn_fig, XMIN,XMAX,YMIN,YMAX)
            CALL FRAME


            if(figures(i) == 'grid' .or. figures(i) == 'grid_VS'.or. figures(i) == 'grid_HS') &
                 call PlotGrid(rnelem, ival, rwi )
            if(figures(i) == 'gridHP' .or. figures(i) == 'gridHP_VS'.or. figures(i) == 'gridHP_HS') &
                 call PlotColorMaps(XMIN,XMAX,YMIN,YMAX, &
                 rnelem,ival,rwi(1:rnelem, 1:7), &
                 ncolors, rgb, pmin, pmax, 'hp-mesh',.true., color)

            if(axes_type /= 'N' .and. axes_type /= 'n') &
                 call DrawAxes(XMIN, XMAX, YMIN, YMAX, axes_type)

            if(figures(i) == 'grid_VS' .or. figures(i) == 'gridHP_VS' ) then

               jn_fig = jn_fig + 1

               CALL MICKEY_MOUSE(in_fig,jn_fig,XMIN,XMAX,YMIN,YMAX)
               !      call frame

               call PlotVscale(XMIN,XMAX,YMIN,YMAX,ncolors, rgb, pmin, pmax, &
                    'hp-mesh', .true. , color, rtime, type )
            endif

            ! horizontal scale
            if(figures(i) == 'grid_HS' .or. figures(i) == 'gridHP_HS' ) then

               in_fig = in_fig + 1

               CALL MICKEY_MOUSE(in_fig,jn_fig,XMIN,XMAX,YMIN,YMAX)
               !      call frame

               call PlotHscale(XMIN,XMAX,YMIN,YMAX,ncolors, rgb, pmin, pmax, &
                    'hp-mesh', .true. , color, rtime, type )
            endif


         else
            print*,'unknown type in cfig.f90', figures(i)

        endif

      enddo







!      CALL MICKEY_MOUSE(2,1,XMIN,XMAX,YMIN,YMAX)
!      call frame

!!      call PlotColorMaps(rnelem, ival, rwi, ncolors, rgb, pmin, pmax, .false. )
!
!
!      CALL MICKEY_MOUSE(2,2,XMIN,XMAX,YMIN,YMAX)
!      call frame
!
!      call PlotVscale(XMIN,XMAX,YMIN,YMAX,ncolors, rgb, pmin, pmax, .true. )









!      CALL THICK_PIXEL(15)
!      CALL LINE(0.45,0.6,0.7,0.85)
!      CALL LINE(0.95,0.6,0.7,0.85)
!      CALL THICK_PIXEL(30)
!      CALL LINE(0.6,0.75,0.6,0.85)

!C ----- drawing thick lines ----
!C  CALL THICK_PIXEL(ITHICK) makes all subsequent lines
!C  ITHICK pixels thick.
!C  CALL LINE(XB,YB,XE,YE) draws a line connecting XB,YB with XE,YE.

!      NDES=101
!      DO IX=1,NDES
!         DESM1=NDES-1
!         X=0.1+0.7*(IX-1)/DESM1
!         Y=0.15+0.1*SIN(20.*X)
!         XDES(IX)=X
!         YDES(IX)=Y
!      END DO
!      CALL THICK_PIXEL(6)
!      CALL LING(XDES,YDES,NDES)

!C ----- drawing a curve -------
!C  CALL LING(XDES,YDES,NDES) connects NDES points with
!C  given coordinates XDES(I),YDES(I),I=1,NDES.


!      CALL THICK_PIXEL(4)
!      Y=0.15+0.1*SIN(20.*0.8)
!      CALL BUBBLE(0.8,Y,1.2,1.,'circle')

!C ----- draws a symbol -------
!C  CALL BUBBLE(X,Y,AMM,GREY,SYM): X,Y is the position,
!C  AMM size approximately in mm, GREY (between 0 and 1)
!C  brightness of filling, SYM form of symbol with many choices:
!C  'circle','square','losange','xstar','pstar','utri','dtri','umtri',
!C  'dmtri','hardy','laurel','mpent','mhex','pent','hex','mmhex'.


!      CALL RGB_BUBBLE(0.2,0.7,8.,1.,0.5,0.,'square')

!      CALL RGB_BUBBLE(0.3,0.2,5.,1.,0.5,0.,'square')

!C ---- a colored symbol ------
!C  the parameter GREY in BUBBLE is replaced by three
!C  "red,green,blue" brightness values.


!      CALL TEXT_LATEX_W(.7,.55,.5,.5,'$E=mc^2$')
!      CALL TEXT_LATEX(0.1,0.25,0.,0.,'important formula')

!C ------ write LaTeX formulas inside the picture ----
!C  CALL TEXT_TEX(X,Y,ALAM,AMU,CHAR) writes LaTeX formula CHAR
!C  at position X,Y. The parameters ALAM,AMU are the relative
!C  coordinates (between 0 and 1) of the position point X,Y
!C  inside the box and allow all types of centering.
!C  Values outside the interval [0,1] are also possible.
!C  CALL TEXT_TEX_W(X,Y,ALAM,AMU,CHAR) puts the formula
!C  inside a white box.


!      CALL THICK_PIXEL(6)
!      CALL ARROW(0.4,0.35,0.6,0.5,3.)

!C ----- draw an arrow -----
!C  CALL ARROW(XB,YB,XE,YE,AMM) has same arguments as LINE, with
!C  AMM indicating size of arrow tip.


!!      CALL GREATER_BOUNDINGBOX(0.75,0.5,0.6,0.4)
!      CALL GREATER_BOUNDINGBOX(1.0,-3.3,0.5,0.5)

!      CALL GREATER_BOUNDINGBOX(1.5,-3.3,0.5,0.5)
      CALL GREATER_BOUNDINGBOX(1.5,0.3,0.5,0.5)

!      CALL COLOR_BOUNDINGBOX('yellow')
!c      CALL RGB_COLOR_BOUNDINGBOX(1.,0.9,0.6)

!c ---- correct bounding box (LEFT,RIGHT,DOWN,UP; IN CM)----
!C  for fine adjustment of the picture inside the text


      CALL END_GGG

!C ---- and that's the end -----
!C  writes the necessary marcos and produces final files


      STOP
    END program cfig


!> ----- drawing the axes -------
!>  CALL AXEX(XORIG,YORIG,XMIN,XMAX,XINCR,KX) draws x-axix with
!>  XORIG,YORIG position of origin, XMIN,XMAX delimiters, XINCR
!>  distance of large graduations, KX number of subgraduations.
!>  CALL LABEL_X(XORIG,YORIG,XMIN,XMAX,XINCR,IDIGIT) writes labels,
!>  IDIGIT is number of digits,
    subroutine DrawAxes(XMIN, XMAX, YMIN, YMAX,axes_type )
      real, intent(in) :: XMIN, XMAX, YMIN, YMAX
      character (len = *), intent (in) :: axes_type
      integer :: KX, KY, IDIGIT
      real :: XINCR, YINCR

      XINCR = (XMAX - XMIN) / 2
      if(axes_type == 'F2') XINCR = (XMAX - XMIN) / 1
      !XINCR=0.5
      !print*,'XINCR=',XINCR, XMAX, XMIN
      KX=5
!E
      CALL AXEX(XMIN, YMIN, XMIN,XMAX,XINCR,KX)

      if(axes_type == 'F' .or. axes_type == 'F2' .or. axes_type == 'Fx') then
      !   IDIGIT=4
         IDIGIT=5
         if(axes_type == 'F2') IDIGIT=2
         CALL LABEL_X(XMIN, YMIN, XMIN,XMAX,XINCR,IDIGIT)

      else if(axes_type == 'E' .or. axes_type == 'Ex') then
         IDIGIT=6
         CALL LABEL_X_E(XMIN, YMIN, XMIN,XMAX,XINCR,IDIGIT)
      endif

      YINCR = (YMAX - YMIN) / 2
      if(axes_type == 'F2') YINCR = (YMAX - YMIN) / 1
      !YINCR=0.5
      KY=5

      CALL AXEY(XMIN, YMIN, YMIN,YMAX,YINCR,KY)

      if(axes_type == 'F' .or. axes_type == 'F2' .or. axes_type == 'Fy') then
         IDIGIT=6
         if(axes_type == 'F2') IDIGIT=2
         CALL LABEL_Y(XMIN, YMIN, YMIN,YMAX,YINCR,IDIGIT)

      elseif(axes_type == 'E' .or. axes_type == 'Ey') then
         IDIGIT=6
         CALL LABEL_Y_E(XMIN, YMIN, YMIN,YMAX,YINCR,IDIGIT)
      endif

      !return

      !F
      !KX=2
      !XINCR=0.5
      !CALL AXEX(0.,0.,XMIN,XMAX,XINCR,KX)
      !IDIGIT=5
      !CALL LABEL_X(0.,0.,XMIN,XMAX,XINCR,IDIGIT)
      !YINCR=0.5
      !KY=2
      !CALL AXEY(0.,0.,0.01,YMAX,YINCR,KY)
      !IDIGIT=5
      !CALL LABEL_Y(0.,0.,0.1,YMAX,YINCR,IDIGIT)
    end subroutine DrawAxes

    subroutine NodeInTriang(xp, xt, interior)
      real, dimension(1:2), intent(in) :: xp  ! node
      real, dimension(1:3, 1:2), intent(in) :: xt  ! vertices of triangle
      logical, intent(inout) :: interior
      real :: det
      integer :: j, j1


      !write(*,*) xt(1, 1:2)
      !write(*,*) xt(2, 1:2)
      !write(*,*) xt(3, 1:2)
      !write(*,*) xp(1:2)

      interior = .false.

      do j=1,3
         j1 = mod(j, 3) + 1
         det = xt(j, 1) * ( xt(j1, 2) - xp(2) ) + xt(j1, 1) * ( xp(2) - xt(j,2) )  + xp(1) * ( xt(j, 2) - xt(j1,2))
         if(det < 0.) return
      enddo
      interior = .true.

      !print*,'##############  inside',  NodeInTriang

    end subroutine NodeInTriang


    subroutine FindIntersectFrameTriang(XMIN, XMAX, YMIN, YMAX, x1, x2, x3, nod, xi )
      real, intent(in) :: XMIN,XMAX,YMIN,YMAX
      real, dimension(1:2), intent(in) :: x1,  x2,  x3   ! triange
      real, dimension(1:10, 1:2), intent(inout)  :: xi  ! polygon
      integer, intent(out) :: nod
      real, dimension(1:4, 1:2) :: xq
      real, dimension(1:3, 1:2) :: xt
      real, dimension(1:2) :: xp
      real :: shift
      integer:: i,j
      logical :: interior

      shift = 1E-6

      ! frame = quadrilaterall
      xq(1, 1) = XMIN
      xq(2, 1) = XMAX
      xq(3, 1) = XMAX
      xq(4, 1) = XMIN

      xq(1, 2) = YMIN
      xq(2, 2) = YMIN
      xq(3, 2) = YMAX
      xq(4, 2) = YMAX

      ! moving the vertices of triangle slihtly into interior
      xt(1, 1:2) = x1(1:2)*(1 - 2*shift) + x2(1:2) * shift + x3(1:2) * shift
      xt(2, 1:2) = x2(1:2)*(1 - 2*shift) + x3(1:2) * shift + x1(1:2) * shift
      xt(3, 1:2) = x3(1:2)*(1 - 2*shift) + x1(1:2) * shift + x2(1:2) * shift

      !write(40,*) xq(1, 1:2)
      !write(40,*) xq(2, 1:2)
      !write(40,*) xq(3, 1:2)
      !write(40,*) xq(4, 1:2)
      !write(40,*) xq(1, 1:2)

      !write(41,*) xt(1, 1:2)
      !write(41,*) xt(2, 1:2)
      !write(41,*) xt(3, 1:2)
      !write(41,*) xt(1, 1:2)

      nod = 0
      do i=1,3
         call InnerNode(XMIN,XMAX,YMIN,YMAX, xt(i,1), xt(i,2), interior )
         !print*,'@@@@@',i,interior, xt(i, 1:2)
         if(interior) then
            nod = nod + 1
            xi(nod, 1:2) = xt(i, 1:2)
         endif
      enddo

      !print*,'#######################################', nod
      if(nod < 3) then ! some inersection detected

         do i=1,4
!UNCOMMENT!
            call NodeInTriang(xq(i,1:2), xt(1:3,1:2), interior )
            !print*,'@@@@@',i,interior, xq(i, 1:2)
            if(interior) then
               nod = nod + 1
               xi(nod, 1:2) = xq(i, 1:2)
            endif
         enddo

         do i=1,4
            i1 = mod(i, 4) + 1
            do j=1,3
               j1 = mod(j, 3) + 1
               call SeekIntersectLine( xq(i, 1:2), xq(i1, 1:2), xt(j, 1:2), xt(j1,1:2), xp(1:2), interior)
               if(interior) then
                  nod = nod + 1
                  xi(nod, 1:2) = xp(1:2)
               endif
            enddo
         enddo

      end if

    end subroutine FindIntersectFrameTriang


    !> plot of the collor map, intersection of all elements with the frame
    !> XMIN,XMAX,YMIN,YMAX is sought
    subroutine PlotColorMaps( XMIN,XMAX,YMIN,YMAX, &
         nelem, ival, wi, ncolors, rgb, rmin, rmax, quantity,  wlines, color)
      COMMON/SIZES/DIM(4)
      real, intent(in) :: XMIN,XMAX,YMIN,YMAX
      integer, intent(in) :: nelem, ival, ncolors
      real, dimension(1:nelem, 1:ival), intent(in) :: wi
      real, dimension(1:ncolors, 1:3), intent(in) :: rgb
      real, intent(in) :: rmin, rmax
      character (len = *), intent(in) :: quantity
      logical, intent(in) :: wlines, color
      real, dimension(1:4, 1:2) :: xq
      real, dimension(1:10, 1:2) :: xi
      integer, dimension(1:10) :: ii, ibx
      real, dimension(1:3, 1:2) :: XDES
      logical, dimension(1:3) :: inner
      real, dimension(1:2) :: xc, xK
      logical :: inter, ident
      real :: ri, qc
      integer :: i,j, j1, j2, l, l1, ic, nnodes, ipoc, ip1, ip2, ip, iiner, ib, elem
      integer :: KX, KY, IDIGIT
      integer :: isub
      !character*20 quantity



      elem = -15174

      xq(1, 1) = XMIN
      xq(2, 1) = XMAX
      xq(3, 1) = XMAX
      xq(4, 1) = XMIN

      xq(1, 2) = YMIN
      xq(2, 2) = YMIN
      xq(3, 2) = YMAX
      xq(4, 2) = YMAX

      do i=1,nelem

         do j =1,3
            XDES(j, 1) = wi(i, 2*(j-1) + 1)
            XDES(j, 2) = wi(i, 2*(j-1) + 2)
         enddo

         !print*,'####'
         !XDES(1,1) = 0.56
         !XDES(1,2) = 0.01
         !XDES(2,1) = 0.58
         !XDES(2,2) = 0.05
         !XDES(3,1) = 0.42
         !XDES(3,2) = 0.02

         !do isub = 1,2

         call FindIntersectFrameTriang(XMIN, XMAX, YMIN, YMAX, XDES(1, :), XDES(2, :), XDES(3, :),&
              nnodes, xi(1:10, 1:2) )

            !print*,'NODES', nnodes

            if(nnodes > 0) then
               if(nnodes > 3) &
                    call CHECK_CONVEXITY(nnodes, xi(1:nnodes, 1:2) )

               ri = wi(i,7)
               ri = min( ri, rmax)
               ri = max( ri, rmin)

               if(rmax == rmin) then
                  ic = 1
                  qc = rmin
               else

                  ic = int((ri - rmin)*0.9999/(rmax-rmin)*(ncolors-1)) + 1
                  qc = (ri - rmin)*0.9999/(rmax-rmin)
                  qc = 1. - (qc*0.9 +0.05)
               endif

               ! color modifications
               if( quantity == 'hp-mesh') then
                  call Set_hp_color(ri, rmin, rmax-rmin, ic)
                  !if(rmax - rmin >= 5.5) then
                  !   if( abs(ri - 5.) < 0.1 ) then  ! DIFFERENCE in PlotVscale
                  !      icn = int(ic* 1.2)
                  !      !print*,'RGBC:', ri, ic, icn, qc, rmax, rmin
                  !      ic = icn
                  !   endif
                  !endif
               endif

               !if( i == elem) then
               !do j=1,nnodes
               !!   !write(*,'(2es12.4,i5,a8)') xi(j,1:2),j,'   !!'
               !   write(50,'(2es12.4,i5,a8)') xi(j,1:2),j,'   !!'
               !enddo
               !endif


               if(color) then
                  CALL RGB_CLOSE_FILL(xi(1:nnodes, 1), xi(1:nnodes, 2), nnodes, &
                       rgb(ic, 1), rgb(ic, 2), rgb(ic,3), wlines)
               else
                  CALL CLOSE_FILL(xi(1:nnodes, 1), xi(1:nnodes, 2), nnodes, &
                       qc, wlines)
               endif

            endif

         !enddo
            !stop
         enddo

      return

      ! OLD VARIANT, SHOULD BE REMOVED
      do i=1,nelem
         iiner = 0
         do j =1,3
            XDES(j, 1) = wi(i, 2*(j-1) + 1)
            XDES(j, 2) = wi(i, 2*(j-1) + 2)
            call InnerNode(XMIN,XMAX,YMIN,YMAX, XDES(j,1), XDES(j,2), inner(j) )
            if(inner(j) ) iiner = iiner + 1
         enddo

         xK(1) = sum(XDES(1:3, 1) ) / 3
         xK(2) = sum(XDES(1:3, 2) ) / 3

         nnodes = 0

         !if(xk(1) > 0.8 .and. xk(1) < 0.85 .and. xk(2) > -0.15 .and. xk(2)< -0.1) then
         !!!if(xk(1) > 0.92 .and. xk(1) < 1.99 .and. xk(2) > -0.25 .and. xk(2)< -0.01) then
         !   write(94, *) '###   elem=',i
         !   write(94, *) xdes(1, 1:2)
         !   write(94, *) xdes(2, 1:2)
         !   write(94, *) xdes(3, 1:2)
         !   write(94, *) xdes(1, 1:2)
         !endif

         !if(inner(1) .or. inner(2) .or. inner(3) ) then
         !   ! at least one inner node
         !   do j=1,3
         !      write(40,'(2es12.4,i5)') XDES(j,1:2),j
         !   enddo
         !   write(40,'(2es12.4,i5, 2es12.4)') XDES(1,1:2),0
         !   write(40,'(x)')
         !   write(45,'(2es12.4,i5)')  xK(1:2)
         !endif

         if( i == elem) then
            print*,'*********************  elem=', i, inner(:)
            write(*,'(a6,3es16.8)') 'xi:',XDES(1:3,1)
            write(*,'(a6,3es16.8)') 'yi:',XDES(1:3,2)
         endif

         ib = 0
         do j=1,3
            j1 = mod(j, 3) + 1
            j2 = mod(j1,3) + 1


            !print*,'@@@@',j,inner(j), inner(:), XDES(j, 1:2), XDES(j1, 1:2)

            if(inner(j) ) then
                  nnodes = nnodes + 1
                  xi(nnodes,1:2) = XDES(j, 1:2)
                  ii(nnodes) = 0
                  !print*,'::::',j,nnodes, xi(nnodes,:)
                  !if(i == elem) write(95,*) xi(nnodes, 1:2), ii(nnodes), 'A',j

                  if(.not. inner(j1) ) then
                     ipoc = 0
                     ident = .false.
                     do l=1,4
                        l1 = mod(l, 4) + 1
                        call SeekIntersectLine(xdes(j,:), xdes(j1,:), xq(l,:), xq(l1,:), xc, inter)

                        if(ipoc > 0) then
                           if(dot_product(xi(nnodes,:)- xc(:), xi(nnodes,:)- xc(:))/ &
                                dot_product(xq(1,:)-xq(3,:), xq(1,:)-xq(3,:)) < 1e-10) &
                                ident = .true.
                        endif

                        !if(i == elem) then
                        !   write(*,'(10es14.6)') xi(nnodes, :)-xc(:),xi(nnodes, :),xc(:)
                        !   write(*,'(10es14.6)') xq(1,:) - xq(3,:), xq(1,:), xq(3,:)
                        !   write(*,'(10es14.6)') &
                        !        dot_product(xi(nnodes,:)- xc(:), xi(nnodes,:)- xc(:)), &
                        !        dot_product(xq(1,:)-xq(3,:), xq(1,:)-xq(3,:)), &
                        !        dot_product(xi(nnodes,:)- xc(:), xi(nnodes,:)- xc(:)) /&
                        !        dot_product(xq(1,:)-xq(3,:), xq(1,:)-xq(3,:))
                        !endif


                        if(i == elem) write(*,*) '###',j,j1,j2,l,l1, ipoc, inter, ident

                        if(inter .and. .not. ident) then
                           !print*,'....',j,l, xc(:)

                           ib = ib + 1
                           ibx(ib) = l

                           ipoc = ipoc + 1
                           if(ipoc == 1) ip1 = l
                           if(ipoc == 2) ip2 = l

                           nnodes = nnodes + 1
                           xi(nnodes,1:2) = xc(1:2)
                           ii(nnodes) = l

                           if(i == elem) then
                           !!!   write(95,*) '## K12',xdes(j,:), xdes(j1,:)
                           !!!   write(95,*) '## D12',xq(l,:), xq(l1,:)
                           !   write(95,*) xi(nnodes, 1:2), ii(nnodes),'B',ib,ibx(ib),ipoc,j
                           endif

                        endif

                        if(ipoc == 2) then ! corner has to be included
                           !write(*,'(a1,3i5,6es12.4)' ) 'A', &
                           !     ip1, ip2, nnodes, xi(nnodes-1, 1:2), xi(nnodes, 1:2), &
                           !     xq(ip1, 1:2)

                           !nnodes = nnodes + 1
                           !xi(nnodes,1:2) = xi(nnodes-1,1:2)
                           !ip1 = mod(ip1, 4) + 1
                           !xi(nnodes-1, 1:2) = xq(ip1, 1:2)
                        endif


                     enddo
                     !!if(i == elem) write(95,*) '### ............'

                  endif
               else

                  ident = .false.
                  ipoc = 0
                  do l=1,4
                     l1 = mod(l, 4) + 1
                     call SeekIntersectLine(xdes(j,:), xdes(j1,:), xq(l,:), xq(l1,:), xc, inter)
                     if(ipoc > 0) then
                        if(dot_product(xi(nnodes,:)- xc(:), xi(nnodes,:)- xc(:))/ &
                             dot_product(xq(1,:)-xq(3,:), xq(1,:)-xq(3,:)) < 1e-10) &
                             ident = .true.
                     endif

                     if(i == elem) then
                        write(*,'(10es14.6)') xdes(j,:)
                        write(*,'(10es14.6)') xdes(j1,:)
                        write(*,'(10es14.6)') xi(nnodes, :)-xc(:),xi(nnodes, :),xc(:)
                        write(*,'(10es14.6)') xq(1,:) - xq(3,:), xq(1,:), xq(3,:)
                        write(*,'(10es14.6)') &
                             dot_product(xi(nnodes,:)- xc(:), xi(nnodes,:)- xc(:)), &
                             dot_product(xq(1,:)-xq(3,:), xq(1,:)-xq(3,:)), &
                             dot_product(xi(nnodes,:)- xc(:), xi(nnodes,:)- xc(:)) /&
                             dot_product(xq(1,:)-xq(3,:), xq(1,:)-xq(3,:))
                        write(*,*) l, l1, ipoc, inter, ident
                     endif

                     if(inter .and. .not. ident) then

                        !print*,',,,,',j,l, xc(:)

                        ib = ib + 1
                        ibx(ib) = l

                        ipoc = ipoc + 1
                        if(ipoc == 1) ip1 = l
                        if(ipoc == 2) ip2 = l

                        nnodes = nnodes + 1
                        xi(nnodes,1:2) = xc(1:2)
                        ii(nnodes) = l
                        !if(i == elem) write(95,*) xi(nnodes, 1:2), ii(nnodes),'C',ib,ibx(ib),ipoc,j

                     endif
                  enddo

                  if(ipoc == 2 .and. .not. inner(j2) ) then ! corner has to be included
                     !write(*,'(a1,3i5,6es12.4)' ) 'B', &
                     !     ip1, ip2, nnodes, xi(nnodes-1, 1:2), xi(nnodes, 1:2), &
                     !     xq(ip1, 1:2)
                     nnodes = nnodes + 1
                     xi(nnodes,1:2) = xi(nnodes-1,1:2)
                     ip = mod(min(ip1, ip2), 4) + 1     ! IT IS CORRECT???
                     xi(nnodes-1, 1:2) = xq(ip, 1:2)
                     ii(nnodes) = -ip
                     !if(i == elem) write(95,*) xi(nnodes-1, 1:2), ii(nnodes-1),'Da'
                     !if(i == elem) write(95,*) xi(nnodes, 1:2), ii(nnodes),'Db'
                  endif

               endif
            enddo

            if(ib == 2 .and. ibx(1) /= ibx(2) ) then
               ! intersection of two different edges, we include the corner
               ip = max(ibx(1), ibx(2))
               if( (ibx(1) == 1 .and. ibx(2) == 4)  &
                    .or. (ibx(2) == 1 .and. ibx(1) == 4) )then
                  ip = 1
               endif
               do l=1,nnodes
                  if( dot_product(xi(l, 1:2) - xq(ip, 1:2), xi(l, 1:2) - xq(ip, 1:2))/ &
                       dot_product(xq(1,:)-xq(3,:), xq(1,:)-xq(3,:)) < 1e-10) goto 20
               enddo
               nnodes = nnodes + 1
               xi(nnodes,1:2) = xq(ip, 1:2)
               !!if(i == elem) write(95,*) xi(nnodes, 1:2), ii(nnodes),'E'

               !print*,'////',ibx(1:2), xq(ip, :)

20             continue
            endif


            !do l=1,nnodes
            !   print*,'###',l,ii(l), xi(l,1:2)
            !enddo

         !!endif


         if(nnodes > 0) then
            if(nnodes > 3) &
                 call CHECK_CONVEXITY(nnodes, xi(1:nnodes, 1:2) )

            ri = wi(i,7)
            ri = min( ri, rmax)
            ri = max( ri, rmin)

            if(rmax == rmin) then
               ic = 1
               qc = rmin
            else

               ic = int((ri - rmin)*0.9999/(rmax-rmin)*(ncolors-1)) + 1
               qc = (ri - rmin)*0.9999/(rmax-rmin)
               qc = 1. - (qc*0.9 +0.05)
            endif

            !print*,'###',qc, rmax, rmin

!!            if(ri > 0.) print*,'!!!!',ri,rmax,rmin, ic

            if( i == elem) then
               do j=1,nnodes
                  write(*,'(2es12.4,i5,a8)') xi(j,1:2),j,'   !!'
               enddo
            endif

            !write(50,'(2es12.4,i5,2es12.4)') xi(1,1:2),0, xK(1:2)
            !write(50,'(x)')


            if(color) then
               CALL RGB_CLOSE_FILL(xi(1:nnodes, 1), xi(1:nnodes, 2), nnodes, &
                 rgb(ic, 1), rgb(ic, 2), rgb(ic,3), wlines)
            else
               CALL CLOSE_FILL(xi(1:nnodes, 1), xi(1:nnodes, 2), nnodes, &
                    qc, wlines)
            endif

         endif

      end do

    end subroutine PlotColorMaps

    !> check the convexity of the polygon
    subroutine CHECK_CONVEXITY(nnodes, xi )
      integer, intent(in) :: nnodes
      real, dimension(1:nnodes, 1:2), intent(inout) :: xi
      real, dimension(1:2) :: xp
      real :: r
      integer :: j, j1, j2, ipoc, i, ip

      ip = 0

      do i=1,nnodes
         ipoc = 0
         do j=1,nnodes
            j1 = mod(j, nnodes) + 1
            j2 = mod(j1,nnodes) + 1

            r = DET(xi(j, 1:2), xi(j1, 1:2), xi(j2, 1:2) )

            ! switching of the order
            if(r < 0) then
               xp(1:2) = xi(j2, 1:2)
               xi(j2, 1:2) = xi(j1, 1:2)
               xi(j1, 1:2) = xp(1:2)

               ipoc = ipoc + 1
               ip = ip + 1
            endif
         enddo
         if( ipoc == 0) goto 100 ! no change, all is OK
      enddo

100   continue

      !if(ip >0) then
      !   do j=1,nnodes
      !      write(96,*) xi(j, 1:2)
      !   enddo
      !   write(96,*) xi(1, 1:2)
      !   write(96,*)
      !endif

    end subroutine CHECK_CONVEXITY

    function DET (x1, x2, x3)
      real :: DET
      real, dimension (1:2), intent(in) :: x1, x2, x3

      DET = x1(1) * ( x2(2) - x3(2) ) + x2(1) * ( x3(2) - x1(2) )  + x3(1) * ( x1(2) - x2(2) )
    end function DET


    !> intersection of abscicas (xa, xb) and (xc, xd) is xi if inner = .true.
    subroutine SeekIntersectLine( xa, xb, xc, xd, xi, inter)
      real, dimension(1:2), intent(in) :: xa, xb, xc, xd
      real, dimension(1:2), intent(out) :: xi
      logical, intent(out) :: inter
      real, dimension (1:2, 1:2) :: A, A1
      real, dimension (1:2) :: b, t
      real :: det

      inter = .false.

      A(1:2,1) =   xb(1:2) - xa(1:2)
      A(1:2,2) = -(xd(1:2) - xc(1:2) )

      b(1:2) = xc(1:2) - xa(1:2)

      det = A(1,1) * A(2,2) - A(1,2) * A(2,1)
      if(det == 0) return  ! paralell lines

      A1(1,1) =  A(2,2) / det
      A1(1,2) = -A(1,2) / det
      A1(2,1) = -A(2,1) / det
      A1(2,2) =  A(1,1) / det

      t(1:2) = matmul(A1(1:2,1:2), b(1:2) )


      if(t(1) >= 0. .and. t(1) <= 1. .and. t(2) >= 0. .and. t(2) <= 1. ) then
         xi(1:2) = xa(1:2) + t(1) * (xb(1:2)  - xa(1:2) )
         inter = .true.

         !write(95,'(2es12.4,a5, 2es12.4)') &
         !     xi(1:2),'  |  ', t(1), t(2)
         !xi(1:2),'  |  ', xc(1:2) + t(2) * (xd(1:2)  - xc(1:2) )

         !print*,'%%%',matmul(A, A1)

         !write(100,*) xa(1:2)
         !write(100,*) xb(1:2)
         !write(100,*) '  '
         !write(100,*) xc(1:2)
         !write(100,*) xd(1:2)
         !write(100,*) '  '

         !write(101,*) xi(1:2)
      endif


      !write(*,'(2es12.4,a5, 2es12.4)') A(1,1:2),'  |  ',A1(1,1:2)
      !write(*,'(2es12.4,a5, 2es12.4)') A(2,1:2),'  |  ',A1(2,1:2)
      !print*,'----------------------------------'

    end subroutine SeekIntersectLine

    !> is (xi, yi) inner node of the fram XMIN,XMAX,YMIN,YMAX ? If YES then inner = .true.
    subroutine InnerNode(XMIN,XMAX,YMIN,YMAX, xi, yi, inner )
      real, intent(in) :: XMIN,XMAX,YMIN,YMAX, xi, yi
      logical, intent(out) :: inner

      inner = .false.
      if(xi <= XMAX .and. xi >= XMIN .and. yi <= YMAX .and. yi >= YMIN ) inner = .true.

    end subroutine InnerNode

    ! x1 is inner, x2 is outer
    subroutine SeekIntersect(XMIN,XMAX,YMIN,YMAX, xA ,xB, xc)
      real, intent(in) :: XMIN,XMAX,YMIN,YMAX
      real, dimension(1:2), intent(in) :: xA, xB
      real, dimension(1:2), intent(out) :: xc
      real :: r

      if(xB(1) > XMAX) then
         r = (XMAX - xA(1) ) / (xB(1) - xA(1) )
         xc(1) = XMAX
         xc(2) = xA(2) + r *( xB(2) - xA(2) )

      elseif(xB(1) < XMIN) then
         r = (xA(1) - XMIN ) / (xA(1) - xB(1) )
         xc(1) = XMIN
         xc(2) = xA(2) + r *( xB(2) - xA(2) )

      elseif(xB(2) > YMAX) then
         r = (YMAX - xA(2) ) / (xB(2) - xA(2) )
         xc(2) = YMAX
         xc(1) = xA(1) + r *( xB(1) - xA(1) )

      elseif(xB(2) < YMIN) then
         r = (xA(2) - YMIN ) / (xA(2) - xB(2) )
         xc(2) = YMIN
         xc(1) = xA(1) + r *( xB(1) - xA(1) )

      else
         print*,'Trouble in SeekIntersect'
         stop
      endif
    end subroutine SeekIntersect

    !> old subroutine, neighbouring elements are not solved
    subroutine PlotColorMapsOLD( XMIN,XMAX,YMIN,YMAX, &
         nelem, ival, wi, ncolors, rgb, rmin, rmax, quantity,  wlines)
      COMMON/SIZES/DIM(4)
      real, intent(in) :: XMIN,XMAX,YMIN,YMAX
      integer, intent(in) :: nelem, ival, ncolors
      real, dimension(1:nelem, 1:ival), intent(in) :: wi
      real, dimension(1:ncolors, 1:3), intent(in) :: rgb
      real, intent(in) :: rmin, rmax
      character (len = *), intent(in) :: quantity
      logical, intent(in) :: wlines
      real, dimension(1:3) :: XDES, YDES
      real :: xc, yc
      real :: ri
      integer :: i,j, ic
      integer :: KX, KY, IDIGIT

      do i=1,nelem
         do j =1,3
            XDES(j) = wi(i, 2*(j-1) + 1)
            YDES(j) = wi(i, 2*(j-1) + 2)
         enddo
         xc = sum(XDES(:) ) /3
         yc = sum(YDES(:) ) /3

         !print*,'*****',XMIN,XMAX,YMIN,YMAX
         !print*,'@@@@@', minval(XDES(:)), maxval(XDES(:)), minval(YDES(:)), maxval(YDES(:))


         !if(minval(XDES(:)) >= XMIN .and. maxval(XDES(:)) <= XMAX .and. &
         !     minval(YDES(:)) >= YMIN .and. maxval(YDES(:)) <= YMAX ) then

         if(xc >= XMIN .and. xc <= XMAX .and. yc >= YMIN .and. yc <= YMAX ) then

            ri = wi(i,7)
            ri = min( ri, rmax)
            ri = max( ri, rmin)

            ic = int((ri - rmin)*0.9999/(rmax-rmin)*(ncolors-1)) + 1

            !write(*,'(a4,12es10.2,i5)') 'rg?',xdes(1:3),ydes(1:3), wi(i,7),rmin,rmax, &
            !     rgb(ic,:) , ic
            !write(*,'(a4,2i5,12es10.2)') 'rg?',i,ic, wi(i,7),rmin,rmax, ri,rgb(ic,:)

            CALL RGB_CLOSE_FILL(XDES,YDES,3, rgb(ic, 1), rgb(ic, 2), rgb(ic,3), wlines)
         endif

      end do

    end subroutine PlotColorMapsOLD

    subroutine PlotGrid(nelem, ival, wi )
      COMMON/SIZES/DIM(4)
      integer, intent(in) :: nelem, ival
      real, dimension(1:nelem, 1:ival), intent(in) :: wi
      real, dimension(1:3) :: XDES, YDES
      real :: ri
      integer :: i,j, ic

      !CALL THICK_PIXEL(0.5)

      do i=1,nelem
         do j =1,3
            XDES(j) = wi(i, 2*(j-1) + 1)
            YDES(j) = wi(i, 2*(j-1) + 2)
         enddo

         CALL LINE(XDES(1), YDES(1), XDES(2), YDES(2) )
         CALL LINE(XDES(2), YDES(2), XDES(3), YDES(3) )
         CALL LINE(XDES(3), YDES(3), XDES(1), YDES(1) )

      end do

    end subroutine PlotGrid


    ! plotting of a vertical scale
    subroutine PlotVscale(XMIN,XMAX,YMIN,YMAX, ncolors, rgb, rmin, rmax, quantity, &
         discr, color, rtime, type )
      COMMON/SIZES/DIM(4)
      real, intent(in) :: XMIN,XMAX,YMIN,YMAX
      integer, intent(in) :: ncolors
      real, dimension(1:ncolors, 1:3), intent(in) :: rgb
      real, intent(in) :: rmin, rmax, rtime
      character (len = *), intent(in) :: quantity, type
      character (len = 180) :: quantity1
      logical, intent(in) :: discr, color
      real, dimension(1:4) :: XDES, YDES
      character(len = 20) :: text
      character(len = 1) :: ch1
      character(len = 5) :: ch5
      character(len = 10) :: ch10
      character(len = 12) :: ch12
      real :: ri, rx, ry, ry1, font_scale, rc
      integer :: i,j, ic, nc, tlen


      rx = 0.1
      ry = 0.8

      ry1 = 0.84
      ry2 = 0.95

      font_scale = 1.
      if(.not. discr) font_scale = 0.85

      if(type == 'battery') then
         font_scale = 1.25
         rx = 0.15
         ry = 0.9
      endif

      !font_scale = 1.

      call SCALE_CHAR(font_scale)

      !print*,'#####E##E$R$#', xmin

      XDES(1) = xmin
      XDES(2) = xmin + (xmax - xmin) * rx

      XDES(3) = XDES(2)
      XDES(4) = XDES(1)

      if(discr) then
         nc = rmax - rmin + 1
      else
         nc = ncolors
      endif

      do i=0,nc-1

         YDES(1) = ymin +  1.*i / nc * (ymax - ymin) * ry
         YDES(3) = ymin +  1.*(i+1) / nc * (ymax - ymin) * ry
         YDES(2) = YDES(1)
         YDES(4) = YDES(3)

         ic = int( 1.* i / (nc-1)  *(ncolors-1)) + 1
         rc = (1.* i / (nc-1))
         rc = 1. - (rc*0.9 +0.05)

         ! color modifications
         if( quantity == 'hp-mesh') then
            ri = rmin + 1.*i/(nc-1) * (rmax - rmin)
            call Set_hp_color(ri, rmin, rmax-rmin, ic)
            !print*,'###',i,ri, ic
         endif


         !write(*,'(a4,i5, 12es10.2,i5)') 'rg?',i,xdes(1:4),ydes(1:4),rmin,rmax, &
         !     rgb(ic,:)
         !write(*,'(a4,2i5,12es10.2)') 'rg?',i,ic, rmin,rmax, rgb(ic,:)

         if(color) then
            CALL RGB_CLOSE_FILL(XDES,YDES,4, rgb(ic, 1), rgb(ic, 2), rgb(ic,3), .false.)
         else
            CALL CLOSE_FILL(XDES,YDES,4, rc, .false.)
         endif

      end do


      ! lines of the scale
      XDES(3) =  XDES(2) + 0.5*(XDES(2) - XDES(1) )

      if(type == 'battery') XDES(3) =  XDES(2) + 0.25*(XDES(2) - XDES(1) )

      if(.not. discr) nc = 6
      if(discr) text(1:2) = 'P_'


      do i=0,nc

         YDES(1) = ymin +  1.*i / nc * (ymax - ymin) * ry
         YDES(2) = YDES(1)

         YDES(3) = ymin +  1.*(i-0.125) / nc * (ymax - ymin) * ry
         if(discr) YDES(3) = ymin +  1.*(i+0.25) / nc * (ymax - ymin) * ry

         !CALL THICK_PIXEL(1)
         CALL LINE(XDES(1), YDES(1), XDES(2), YDES(2) )


         if(discr) then
            write( ch1, '(i1)' ) int(rmin+i) ! change the format

            text(3:3) = ch1(1:1)
            tlen = 3
         else
            write( ch10, '(es10.2)') rmin + 1.*i/nc *(rmax-rmin) ! change the format
            if(quantity(1:3) == 'Alg') &
                 write( ch10, '(es10.1)') rmin + 1.*i/nc *(rmax-rmin)

            ! logarithmic scale
            if( quantity(1:1) == 'E' .and. len(quantity) > 0) &
                 write( ch10, '(es10.2)') 10**(rmin + 1.*i/nc *(rmax-rmin))

            text(1:10) = ch10(1:10)

            tlen = 10
            !print*,'####',ch10, tlen, rmin + 1.*i/nc * (rmax-rmin)
         endif

         if(.not. discr .or. i < nc) then
            call TEXT_W(XDES(3), YDES(3), text(1:tlen) )
         endif

      enddo
      call SCALE_CHAR(1./font_scale)


      YDES(1) = ymin +  (ymax - ymin) * ry1

      call Set_Quantity1( quantity, type, quantity1)

      font_scale = 7.
      !if(type == 'battery') font_scale = 0.1

      call SCALE_CHAR(font_scale)

      !call TEXT_W(XDES(1), YDES(1)-0.2, quantity1 )
      call TEXT_LATEX(XDES(1), YDES(1), 0.0, 0.0, quantity1 )
      call SCALE_CHAR(1./font_scale)

      YDES(1) = ymin +  (ymax - ymin) * ry2
      XDES(1) = XDES(1) - (XDES(2) - XDES(1)) / 1.5

      if(rtime < 1E+10) then
         font_scale = 1.
         call SCALE_CHAR(font_scale)
         write( ch10, '(es10.2)') min(rtime, 999999.)
         ch12(1:2) = 't='
         ch12(3:12) = ch10(1:10)
         !!!!call TEXT_LATEX(XDES(3), YDES(1), 0.0, 0.0, ch12 )
         call TEXT_LATEX(XDES(1), YDES(1), 0.0, 0.0, ch12 )
         call SCALE_CHAR(1./font_scale)
         !print*,'######',ch10, quantity1
      endif

    end subroutine PlotVscale



    ! plotting of a horozintal scale
    subroutine PlotHscale(XMIN,XMAX,YMIN,YMAX, ncolors, rgb, rmin, rmax, quantity, &
         discr, color, rtime, type )
      COMMON/SIZES/DIM(4)
      real, intent(in) :: XMIN,XMAX,YMIN,YMAX
      integer, intent(in) :: ncolors
      real, dimension(1:ncolors, 1:3), intent(in) :: rgb
      real, intent(in) :: rmin, rmax, rtime
      character (len = *), intent(in) :: quantity, type
      character (len = 180) :: quantity1
      logical, intent(in) :: discr, color
      real, dimension(1:4) :: XDES, YDES
      character(len = 20) :: text
      character(len = 1) :: ch1
      character(len = 5) :: ch5
      character(len = 10) :: ch10
      character(len = 12) :: ch12
      real :: ri, rx, rx1, rx2, font_scale, rc
      integer :: i,j, ic, nc, tlen


      rx = 0.8
      ry = 0.1

      rx1 = 0.85
      rx2 = 0.95

      font_scale = 0.8
      if(.not. discr) font_scale = 0.6

      !font_scale = 1.

      call SCALE_CHAR(font_scale)

      YDES(1) = ymax - 0.0 * (ymax - ymin) * ry
      YDES(2) = ymax - 1.5* (ymax - ymin) * ry

      YDES(3) = YDES(2)
      YDES(4) = YDES(1)

      if(discr) then
         nc = rmax - rmin + 1
      else
         nc = ncolors
      endif

      do i=0,nc-1

         xDES(1) = xmin +  1.*i / nc * (xmax - xmin) * rx
         xDES(3) = xmin +  1.*(i+1) / nc * (xmax - xmin) * rx
         xDES(2) = xDES(1)
         xDES(4) = xDES(3)

         ic = int( 1.* i / (nc-1)  *(ncolors-1)) + 1
         rc = (1.* i / (nc-1))
         rc = 1. - (rc*0.9 +0.05)

         !write(*,'(a4,i5, 12es10.2,i5)') 'rg?',i,xdes(1:4),ydes(1:4),rmin,rmax, &
         !     rgb(ic,:)
         !write(*,'(a4,2i5,12es10.2)') 'rg?',i,ic, rmin,rmax, rgb(ic,:)

         if(color) then
            CALL RGB_CLOSE_FILL(XDES,YDES,4, rgb(ic, 1), rgb(ic, 2), rgb(ic,3), .false.)
         else
            CALL CLOSE_FILL(XDES,YDES,4, rc, .false.)
         endif

      end do


      ! lines of the scale
      YDES(3) =  YDES(2) + 0.75*(YDES(2) - YDES(1) )

      if(.not. discr) nc = 6
      if(discr) text(1:2) = 'P_'


      do i=0,nc

         XDES(1) = xmin +  1.*i / nc * (xmax - xmin) * rx
         XDES(2) = XDES(1)

         XDES(3) = xmin +  1.*(i-0.25) / nc * (xmax - xmin) * rx
         if(discr) XDES(3) = xmin +  1.*(i+0.25) / nc * (xmax - xmin) * rx

         !CALL THICK_PIXEL(1)
         CALL LINE(XDES(1), YDES(1), XDES(2), YDES(2) )


         if(discr) then
            write( ch1, '(i1)' ) int(rmin+i) ! change the format

            text(3:3) = ch1(1:1)
            tlen = 3
         else
            write( ch10, '(es10.2)') rmin + 1.*i/nc *(rmax-rmin) ! change the format
            if(quantity(1:3) == 'Alg') &
                 write( ch10, '(es10.1)') rmin + 1.*i/nc *(rmax-rmin)

            ! logarithmic scale
            if( quantity(1:1) == 'E' .and. len(quantity) > 0) &
                 write( ch10, '(es10.2)') 10**(rmin + 1.*i/nc *(rmax-rmin))

            text(1:10) = ch10(1:10)

            tlen = 10
            !print*,'####',ch10, tlen, rmin + 1.*i/nc * (rmax-rmin)
         endif

         if(.not. discr .or. i < nc) then
            call TEXT_W(XDES(3), YDES(3), text(1:tlen) )
         endif

      enddo
      !call SCALE_CHAR(1./font_scale)
       font_scale = 4.
      call SCALE_CHAR(font_scale)


      XDES(1) = xmin +  (xmax - xmin) * rx1
      YDES(1) = ymax + 0. * (ymax - ymin) * ry


      call Set_Quantity1( quantity, type, quantity1)

      font_scale = 4.
      call SCALE_CHAR(font_scale)

      !call TEXT_W(XDES(1), YDES(1)-0.2, quantity1 )
      call TEXT_LATEX(XDES(1), YDES(1), 0.0, 0.0, quantity1 )
      call SCALE_CHAR(1./font_scale)

      XDES(1) = xmin +  (xmax - xmin) * rx1
      YDES(1) = ymax - 1.5 * (ymax - ymin) * ry

      if(rtime < 999999.) then
         font_scale = 10.
         call SCALE_CHAR(font_scale)
         write( ch10, '(es10.2)') min(rtime, 999999.)
         ch12(1:2) = 't='
         ch12(3:12) = ch10(1:10)
         !!!!call TEXT_LATEX(XDES(3), YDES(1), 0.0, 0.0, ch12 )
         call TEXT_LATEX(XDES(1), YDES(1), 0.0, 0.0, ch12 )
         call SCALE_CHAR(1./font_scale)
         !print*,'######',ch10, quantity1
      endif

    end subroutine PlotHscale


    subroutine Set_Quantity1( quantity, type, quantity1)
      character (len = *), intent(in) :: quantity, type
      character (len = *), intent(inout) :: quantity1

      quantity1(:) = quantity(:)
      if(quantity(:) == 'S1')   quantity1 = '{\Large $\eta_T^n $}'
      if(quantity(:) == 'S234') quantity1 = '{\Large $e_T^n$ }'
      if(quantity(:) == 'S5')   quantity1 = '{\Large $\eta_{{\rm sp},T}^n $}'
      if(quantity(:) == 'S6')   quantity1 = '{\Large $\eta_{{\rm tm},T}^n $}'

      if(quantity(:) == 'E6')   quantity1 = '{\Large $\| e_h\|_{L^\infty(K)} $}'
      if(quantity(:) == 'E7')   quantity1 = '{\Large $\| E^I\|_{L^\infty(K)} $}'
      if(quantity(:) == 'E8')   quantity1 = '{\Large $\| e_h\|_{L^2(K)} $}'
      !if(quantity(:) == 'E9')   quantity1 = '{\Large $\| E^I\|_{L^2(K)} $}'
      !if(quantity(:) == 'E9')   quantity1 = '{\Large $\| e^{\rm int}\|_{L^2(K)} $}'
      if(quantity(:) == 'E9')   quantity1 = '{\Large $g_T / G_T $}'
      if(quantity(:) == 'E6')   quantity1 = '{\Large $g_T / G_T^{T2} $}'

      if(quantity(:) == 'E3')   quantity1 = '{\Large $\eta_{T,{\rm tot}}/\omega_T $}'
      !if(quantity(:) == 'E6')   quantity1 = '{\Large $\eta_K^{{\rm tot},{p-1}} $}'
      !if(quantity(:) == 'E10')   quantity1 = '{\Large $\eta_K/h_K^{p_K}/|K|^{1/2} $}'
      if(quantity(:) == 'E16')   quantity1 = '{\Large $\eta_K/e_{h,K} $}'
      !if(quantity(:) == 'E12')   quantity1 = '{\Large $ e_{h,K}/h_K^{p_K}/|K|^{1/2} $}'

      if(quantity(:) == 'E10')  quantity1 = '{\Large $\eta_K $}'
      if(quantity(:) == 'E11')  quantity1 = '{\Large $\eta_K^{p+1} $}'
      !if(quantity(:) == 'S12')  quantity1 = '{\Large $(\eta_K -\eta_K^{p+1})/\eta_K^{p+1}$}'
      if(quantity(:) == 'S11')  quantity1 = '{\Large $\tilde{C}_K$}'
      if(quantity(:) == 'S12')  quantity1 = '{\Large $C_K$}'
      if(quantity(:) == 'S17')  quantity1 = '{\Large $ama_p$}'
      if(quantity(:) == 'S18')  quantity1 = '{\Large $s_K$}'
      if(quantity(:) == 'E14')  quantity1 = '{\LARGE $|e_h|_{H^1} $}'

      if(quantity(:) == 'hp-mesh')   quantity1 = '{\Large $hp$ mesh}'

      if(quantity(:) == 'Alg1')   quantity1 = '{\Large $e_T^{\rm tot} $}'
      if(quantity(:) == 'Alg2')   quantity1 = '{\Large $e_T^{\rm disc} $}'
      if(quantity(:) == 'Alg3')   quantity1 = '{\Large $e_T^{\rm alg} $}'
      if(quantity(:) == 'Alg4')   quantity1 = '{\Large $\eta_T^{\rm tot} $}'
      if(quantity(:) == 'Alg5')   quantity1 = '{\Large $\eta_T^{\rm disc} $}'
      if(quantity(:) == 'Alg6')   quantity1 = '{\Large $\eta_T^{\rm alg} $}'


      if(type(:) == 'pNeuA') then
          if(quantity(:) == 'S1') quantity1 = '{\LARGE $\nabla( u - u_h)$}'
          if(quantity(:) == 'S6') quantity1 = '{\LARGE $\eta_K^{\rm tot}$}'
          if(quantity(:) == 'S7') quantity1 = '{\LARGE $\|\nabla e_h\|_K$}'
       elseif(type(:) == 'pNeu') then
          if(quantity(:) == 'S3') quantity1 = '{\LARGE $\tilde{g}_K([u_h])$}'
          if(quantity(:) == 'S1') quantity1 = '{\LARGE $| u - u_h|$}'
          if(quantity(:) == 'S5') quantity1 = '{\LARGE $\eta_K^{\mathrm{tot}}$}'
          if(quantity(:) == 'S9') quantity1 = '{\LARGE $\frac{\|\Pi_{-1}\sigma\|_K}{\|\sigma\|_K}$}'
          if(quantity(:) == 'S10') quantity1 = '{\LARGE $\frac{\|\Pi_{-2}\sigma\|_K}{\|\sigma\|_K}$}'
          if(quantity(:) == 'S11') quantity1 = '{\LARGE $\frac{\|\Pi_{-3}\sigma\|_K}{\|\sigma\|_K}$}'
          if(quantity(:) == 'S12') quantity1 = '{\LARGE $\frac{\|\Pi_{-4}\sigma\|_K}{\|\sigma\|_K}$}'
          if(quantity(:) == 'S13') quantity1 = '{\LARGE $\frac{\|\Pi_{-1} s \|_K}{\| s\|_K}$}'
          if(quantity(:) == 'S14') quantity1 = '{\LARGE $\frac{\|\Pi_{-2} s \|_K}{\| s\|_K}$}'
          if(quantity(:) == 'S15') quantity1 = '{\LARGE $\frac{\|\Pi_{-3} s \|_K}{\| s\|_K}$}'
          if(quantity(:) == 'S16') quantity1 = '{\LARGE $\frac{\|\Pi_{-4} s \|_K}{\| s\|_K}$}'

          !if(quantity(:) == 'S10') quantity1 = '{\Large $\Pi_{-2}\sigma$}'
          !if(quantity(:) == 'S11') quantity1 = '{\Large $\Pi_{-3}\sigma$}'
          !if(quantity(:) == 'S12') quantity1 = '{\Large $\Pi_{-4}\sigma$}'
          !if(quantity(:) == 'S13') quantity1 = '{\Large $\Pi_{-1}s$}'
          !if(quantity(:) == 'S14') quantity1 = '{\Large $\Pi_{-2}s$}'
          !if(quantity(:) == 'S15') quantity1 = '{\Large $\Pi_{-3}s$}'
          !if(quantity(:) == 'S16') quantity1 = '{\Large $\Pi_{-4}s$}'
          if(quantity(:) == 'S19') quantity1 = '{\LARGE $h_K^{-1/2} \frac{\|\Pi_{-1}\sigma\|_K}{\|\Pi_{-2}\sigma\|_K}$}'
          if(quantity(:) == 'S20') quantity1 = '{\LARGE $\frac{1}{\sqrt{h_K}} \frac{\|\Pi_{-2}\sigma\|_K}{\|\Pi_{-3}\sigma\|_K}$}'
          if(quantity(:) == 'S21') quantity1 = '{\LARGE $\frac{1}{\sqrt{h_K}} \frac{\|\Pi_{-3}\sigma\|_K}{\|\Pi_{-4}\sigma\|_K}$}'
          if(quantity(:) == 'S22') quantity1 = '{\LARGE $\frac{1}{\sqrt{h_K}} \frac{\|\Pi_{-1} s\|_K }{\|\Pi_{-2} s \|_K}$}'
          if(quantity(:) == 'S23') quantity1 = '{\LARGE $\frac{1}{\sqrt{h_K}} \frac{\|\Pi_{-2} s\|_K }{\|\Pi_{-3} s \|_K}$}'
          if(quantity(:) == 'S24') quantity1 = '{\LARGE $\frac{1}{\sqrt{h_K}} \frac{\|\Pi_{-3} s\|_K }{\|\Pi_{-4} s \|_K}$}'

          if(quantity(:) == 'E7')  quantity1 = '{\Large $|e_h|_K $}'
          if(quantity(:) == 'E10')  quantity1 = '{\Large $\eta_K $}'
          if(quantity(:) == 'E11')  quantity1 = '{\Large $\eta_K /|e_h|_K $}'


       elseif(type(:) == 'HO_rec') then
          if(quantity(:) == 'E9') quantity1 = '{\LARGE $\eta_{L^2}(u_h) $}'
          if(quantity(:) == 'E10') quantity1 = '{\LARGE $\eta_{H^1}(u_h) $}'
          if(quantity(:) == 'E6') quantity1 = '{\LARGE $\eta_{H^1}(\Pi^{p-1}u_h) $}'
          !if(quantity(:) == 'E5') quantity1 = '{\LARGE $\eta_{L^2}^{p+1} $}'
          !if(quantity(:) == 'E6') quantity1 = '{\LARGE $\eta_{H^1}^{p+1} $}'

          if(quantity(:) == 'E13') quantity1 = '{\LARGE $\|e_h\|_{L^2} $}'
          if(quantity(:) == 'E14') quantity1 = '{\LARGE $|e_h|_{H^1} $}'
          if(quantity(:) == 'E15') quantity1 = '{\LARGE $\eta_{L^2}/\|e_h\|_{L^2} $}'
          if(quantity(:) == 'E16') quantity1 = '{\LARGE $\eta_{H^1}^{p+1}/|e_h|_{H^1} $}'
          !if(quantity(:) == 'E17') quantity1 = '{\LARGE $\eta_{H^1}^{p+2}/|e_h|_{H^1} $}'
          !if(quantity(:) == 'E17') quantity1 = '{\LARGE $\eta(u_h)/\eta(u_h^{p-1}) $}'
          if(quantity(:) == 'E17') quantity1 = '{\LARGE $\eta(u_h)/\omega_K $}'
          if(quantity(:) == 'E18') quantity1 = '{\LARGE $g_K/G_K^0 $}'
          if(quantity(:) == 'E19') quantity1 = '{\LARGE $g_K/G_K^2 $}'

      elseif( type(:) == 'DWR' ) then         ! FR for DWR
         if(quantity(:) == 'E2')   quantity1 = '{\Large $\eta_{S}^{p} $}'
         if(quantity(:) == 'E1')   quantity1 = '{\Large $\eta_{A}^{p} $}'
         if(quantity(:) == 'E4')   quantity1 = '{\Large $\eta_{S}^{d} $}'
         if(quantity(:) == 'E3')   quantity1 = '{\Large $\eta_{A}^{d} $}'
         if(quantity(:) == 'E5')   quantity1 = '{\Large $\eta_{S}^{av} $}'
         if(quantity(:) == 'E6')   quantity1 = '{\Large $\eta_{sign}^{p} $}'
         if(quantity(:) == 'E7')   quantity1 = '{\Large $\eta_{sign}^{d} $}'

      elseif( type(:) == 'RES' ) then         ! FR for DWR
         if(quantity(:) == 'E1')   quantity1 = '{\Large $\eta_{K,AS} $}'
         if(quantity(:) == 'E2')   quantity1 = '{\Large $\eta_{K,S} $}'
         if(quantity(:) == 'E3')   quantity1 = '{\Large $\eta_{K,T} $}'
         if(quantity(:) == 'E4')   quantity1 = '{\Large $\eta_{K,ST} $}'
         if(quantity(:) == 'E19')   quantity1 = '{\Large $\eta_T/\eta_{ST} $}'
         if(quantity(:) == 'E20')   quantity1 = '{\Large $\eta_T/\eta_{S} $}'

      endif

      if(type(:) == 'battery') then
         !quantity1 = ' ' !{\LARGE $u$}'
      endif


    end subroutine Set_Quantity1


    subroutine Set_hp_color(ri, rmin, range, ic)
      real, intent(in) :: ri   ! input polyn degree
      real, intent(in) :: rmin  ! minimal polyn degrees
      real, intent(in) :: range   ! range of polyn degrees
      integer, intent(inout) :: ic ! output rgb color between 1, 256
      integer :: pi, p_min

      pi = int (ri+0.5)
      p_min = int (rmin+0.5)

      ! at most five colors
      if(range <= 3.1) then
         if( pi == p_min + 0) ic = 1
         if( pi == p_min + 1) ic = 160
         if( pi == p_min + 2) ic = 200
         if( pi == p_min + 3) ic = 255

      elseif(range <= 4.1) then
         if( pi == p_min + 0) ic = 1
         if( pi == p_min + 1) ic = 90
         if( pi == p_min + 2) ic = 160
         if( pi == p_min + 3) ic = 200
         if( pi == p_min + 4) ic = 255

      else if(range <= 5.1) then
         if( pi == p_min + 0) ic = 1
         if( pi == p_min + 1) ic = 50
         if( pi == p_min + 2) ic = 90
         if( pi == p_min + 3) ic = 160
         if( pi == p_min + 4) ic = 200
         if( pi == p_min + 5) ic = 255

      else if(range <= 7.1) then
         if( pi == p_min + 0) ic = 1
         if( pi == p_min + 1) ic = 50
         if( pi == p_min + 2) ic = 100
         if( pi == p_min + 3) ic = 152  ! yellow
         if( pi == p_min + 4) ic = 168  ! orange
         if( pi == p_min + 5) ic = 210
         if( pi == p_min + 6) ic = 240
         if( pi == p_min + 7) ic = 255

      else if(range <= 8.1) then
         if( pi == p_min + 0) ic = 1
         if( pi == p_min + 1) ic = 50
         if( pi == p_min + 2) ic = 100
         if( pi == p_min + 3) ic = 152  ! yellow
         if( pi == p_min + 4) ic = 168  ! orange
         if( pi == p_min + 5) ic = 200
         if( pi == p_min + 6) ic = 225
         if( pi == p_min + 7) ic = 252
         if( pi == p_min + 8) ic = 255

      else
         if( pi == 0) ic = 1
         if( pi == 1) ic = 25
         if( pi == 2) ic = 50
         if( pi == 3) ic = 75
         if( pi == 4) ic = 100
         if( pi == 5) ic = 152
         if( pi == 6) ic = 168
         if( pi == 7) ic = 185
         if( pi == 8) ic = 220
         if( pi == 9) ic = 245
         if( pi ==10) ic = 255
      endif



    end subroutine Set_hp_color
