!> module for the space adaptation, refinement of triangulations
module adaptation_mod
  use paramets

  implicit none

  public :: SplitTriangToFourNew

  type :: Adaptation_t
     character(len=20) :: adapt_space
     character(len=8) :: adapt_method    !  AMA, RES
     character(len=8) :: adapt_type      !  Ahp, HG, RG, derived from adapt_space
     integer :: max_adapt_level           ! maximal level of mesh adaptation ??

     integer :: adapt_level               ! level of mesh adaptation
     logical :: adapt, adaptR             ! adapt the mesh ?
     real :: tol_min, tol_max             ! tolerances for mesh adaptation

     integer :: type_regularity
     real :: type_regularity_par

     !AMA
     integer, dimension(:,:), allocatable :: AMAhistory   !history of AMA computation
     integer :: stop_adaptation
     real :: Lq

     !HG
     integer :: max_HGlevel               ! maximal implemented HG level 2 for HG/ 0 for others
     integer :: HG                        ! maximal number of hanging nodes per edge

     !RG
     real, dimension(:,:,:), allocatable :: RGred ! red green refinement
     real, dimension(:,:,:,:), allocatable :: RGgreen ! red green refinement

   contains

     procedure :: init => initAdaptation
     procedure :: initRGadapt

  end type Adaptation_t

contains
  
  subroutine initAdaptation( this, tol_max, tol_min, adapt_space, max_adapt_level, Lq )
    class( Adaptation_t ), intent(inout) :: this
    real, intent(in) :: tol_max, tol_min
    character(len=20), intent(in) :: adapt_space
    integer, intent(in) :: max_adapt_level
    real, intent(in), optional :: Lq

    this%tol_max = tol_max
    this%tol_min = tol_min
    this%adapt_space = adapt_space
    this%max_adapt_level = max_adapt_level



    write(debug,*) ' initAdaptation -does connection to mesh work right'

    if (this%adapt_space == 'HGhp' .or. &
         this%adapt_space == 'HGh' .or. &
         this%adapt_space == 'HGp') then
       
       this%adapt_method = 'RES'
       this%adapt_type = 'HG'

       write(debug,*) 'FR: Dont we need initRGrefirement in HG too? '

    elseif (this%adapt_space == 'RGhp' .or. this%adapt_space == 'RGh' &
         .or. this%adapt_space == 'RGp') then
       this%adapt_method = 'RES'
       this%adapt_type = 'RG'

       call this%initRGadapt()

    elseif ( this%adapt_space == 'AMAhp' .or. this%adapt_space == 'AMAh') then
       this%adapt_method = 'Ahp'   !  this%adapt_method should be avoided in future
       this%adapt_type = 'Ahp'
       allocate ( this%AMAhistory(0:max(1, this%max_adapt_level), 1:5) )

    elseif(this%adapt_space == 'ANIhp' .or. this%adapt_space == 'ANIh') then
       this%adapt_method = 'ANI'  !  this%adapt_method should be avoided in future
       this%adapt_type = 'Ahp'
       allocate ( this%AMAhistory(0:max(1, this%max_adapt_level), 1:5) )

    elseif(this%adapt_space == 'IMAhp' .or. this%adapt_space == 'IMAh') then
       this%adapt_method = 'Ihp'   !  this%adapt_method should be avoided in future
       this%adapt_type = 'Ihp'
       allocate ( this%AMAhistory(0:max(1, this%max_adapt_level), 1:5) )

    elseif(this%adapt_space == 'none' .or. this%adapt_space == '-') then
       this%max_adapt_level = 0

    else
       print*,'Unknown mesh adaptation method, possibilities: '
       print*,'          HGh/HGp/HGhp / RGh/RGp/RGhp / AMAh/AMAhp  IMAhp/IMAh'
       stop
    endif
    !
    !  if( this%adapt_type == 'Ahp' .or.  this%adapt_type == 'Ihp')  allocate(AMA)

    if ( present(Lq) ) then
       this%Lq = Lq
    elseif ( this%adapt_type == 'Ahp' .or. this%adapt_type == 'Ihp' ) then
       stop 'CONTROL initAdaptation - Lq should be initialized for AMA,ANI,IMA'
    endif


    ! FERROR what is this for?
    write(debug,*) 'FR: What is tri_solA ???'
    !  this%tri_solA = .false.
    !  if(ch1 == 'Y' .or. ch1 == 'y' .or. ch1 == 'A' .or. ch1 == 'a') state%tri_solA = .true.


    !
    !  if( state%space%adapt%adapt_type == 'Ahp' .or.  state%space%adapt%adapt_type == 'Ihp')  allocate(AMA)
    !



  end subroutine initAdaptation

  ! copy of InitRGrefirement
  !> initialization of the red green refinement
  subroutine initRGadapt( this )
    class( Adaptation_t ), intent( inout ) :: this

    integer :: split = 4
    !real, dimension(1:3, 1:3) :: Rlambda

    if (allocated(this%RGred) ) then
       stop 'already allocated!'
    else
       allocate( this%RGred(0:split, 1:3, 1:3) )
    endif

    this%RGred(0, 1, 1:3) = (/0., 0., 1. /)
    this%RGred(0, 2, 1:3) = (/1., 0., 0. /)
    this%RGred(0, 3, 1:3) = (/0., 1., 0. /)

    call SplitTriangToFourNew(this%RGred(0, 1:3, 1:3), this%RGred( 1:split, 1:3, 1:3) )

    allocate(this%RGgreen(1:3, 1:nbDim, 1:3, 1:3) ) !!! ONLY for triangles

    !state%RG%green(0)%lambda(1, 1, 1:3) = (/0., 0., 1. /)
    !state%RG%green(0)%lambda(1, 2, 1:3) = (/1., 0., 0. /)
    !state%RG%green(0)%lambda(1, 3, 1:3) = (/0., 1., 0. /)

    ! split edge 1
    this%RGgreen(1, 1, 1, 1:3) = (/0.5, 0., 0.5 /)
    this%RGgreen(1, 1, 2, 1:3) = (/1., 0., 0. /)
    this%RGgreen(1, 1, 3, 1:3) = (/0., 1., 0. /)

    this%RGgreen(1, 2, 1, 1:3) = (/0., 0., 1. /)
    this%RGgreen(1, 2, 2, 1:3) = (/0.5, 0., 0.5 /)
    this%RGgreen(1, 2, 3, 1:3) = (/0., 1., 0. /)


    ! split edge 2
    this%RGgreen(2, 1, 1, 1:3) = (/0.5, 0.5, 0. /)
    this%RGgreen(2, 1, 2, 1:3) = (/0., 1., 0. /)
    this%RGgreen(2, 1, 3, 1:3) = (/0., 0., 1. /)

    this%RGgreen(2, 2, 1, 1:3) = (/1., 0., 0. /)
    this%RGgreen(2, 2, 2, 1:3) = (/0.5, 0.5, 0. /)
    this%RGgreen(2, 2, 3, 1:3) = (/0., 0., 1. /)

    ! split edge 3
    this%RGgreen(3, 1, 1, 1:3) = (/0., 0.5, 0.5 /)
    this%RGgreen(3, 1, 2, 1:3) = (/0., 0., 1. /)
    this%RGgreen(3, 1, 3, 1:3) = (/1., 0., 0. /)

    this%RGgreen(3, 2, 1, 1:3) = (/0., 1., 0. /)
    this%RGgreen(3, 2, 2, 1:3) = (/0., 0.5, 0.5 /)
    this%RGgreen(3, 2, 3, 1:3) = (/1., 0., 0. /)


  end subroutine initRGadapt

  !from problem
  !> split triangle with barycentric coordinates x to four triangles with zi)
  subroutine SplitTriangToFourNew(x, zi)
    real, dimension(1:3, 1:3), intent(in) :: x
    real, dimension(1:4, 1:3, 1:3), intent(out) :: zi


    zi(1, 1, 1:3) = x(1,1:3)
    zi(1, 2, 1:3) = (x(1,1:3) + x(2,1:3))/2
    zi(1, 3, 1:3) = (x(3,1:3) + x(1,1:3))/2

    zi(2, 1, 1:3) = x(2,1:3)
    zi(2, 2, 1:3) = (x(2,1:3) + x(3,1:3))/2
    zi(2, 3, 1:3) = (x(1,1:3) + x(2,1:3))/2

    zi(3, 1, 1:3) = x(3,1:3)
    zi(3, 2, 1:3) = (x(3,1:3) + x(1,1:3))/2
    zi(3, 3, 1:3) = (x(2,1:3) + x(3,1:3))/2

    zi(4, 1, 1:3) = (x(1,1:3) + x(2,1:3))/2
    zi(4, 2, 1:3) = (x(2,1:3) + x(3,1:3))/2
    zi(4, 3, 1:3) = (x(3,1:3) + x(1,1:3))/2

  end subroutine SplitTriangToFourNew

end module adaptation_mod
