!> basis geometrical entities which can be stored in parallel mode

module geometry
  use paramets

  implicit none
  !private
  public


  !> matrix block, to each element \f$ K\in{\cal T}_h \f$ corresponds several blocks
  type, public :: Mblock
     real, allocatable, dimension (:, :)  :: Mb  ! block of matrix

     contains

     procedure :: CopyMblock

  end type Mblock

  !> boundary edge, from input and the curved representation, if any

!  type, public :: bound_edge
!     integer :: lbn(4)        ! indexes of end nodes
!     integer :: ibc           ! index of boundary element
!     integer :: itc           ! index of adjacent element
!     integer :: jtc           ! inner index of adjacent element
!     integer :: BC            ! type of boundary condition
!     integer :: inout         ! 0=inlet , 1=outlet, -1=solid walls
!     integer:: icurv          ! icurv = 0, polygonal, icurv =k, P_{k+1} approx
!     real, pointer, dimension(:,:) :: x_inn ! coordinates of inner nodes
!     real, pointer, dimension(:,:) :: x_div ! coordinates of integ nodes for c_L,..
!  end type bound_edge
!
  !> curved parts of the boundary defined by the list of nodes
  type coord
     real, allocatable, dimension(:,:) :: pnt !2D points of a particular part
  end type coord

  type curved_boundary
     integer :: nbp   ! number of curved boundary parts  (=0 => polygonal boundary)
     type(coord), pointer, dimension(:) :: CBP
     ! points of curved boundary parts, CBP(i)%pnt(j,1/2) = coordinate of the j-th
     ! point of the i-th curved boundary part
     logical, allocatable, dimension(:) ::  closed  ! =1 closed boundary part, =0 nonclosed
  end type curved_boundary

  type, public :: MGlevel
    real, allocatable, dimension(:) :: rhs, sol
    integer :: nsize
  end type MGlevel


!  public:: InitElemDof
  public:: InitMblock
  public:: DeleteMblock
  public:: WriteMblock
  public:: WriteMblock_Screene
  public:: VertexIndex
  public:: EdgeIndex

  public :: GetFaceIndexes

  public :: WriteLocalMatrix

  public :: DOFtriang
  public :: DOFtriangR
  public :: DegFromDofTriang

  public:: operator(+), operator(-), operator(*)

  interface operator(+)
    module procedure addMblock
  end interface

  interface operator(-)
    module procedure subMblock
  end interface

  interface operator(*)
    module procedure mulMblock
  end interface

contains

  elemental function DOFtriang(deg)
    integer :: DOFtriang
    integer, intent(in) :: deg

    DOFtriang = (deg + 1) *(deg+2)/ 2
  end function DOFtriang

  elemental function DegFromDofTriang( dof )
    integer, intent(in) :: dof
    integer :: DegFromDofTriang

    DegFromDofTriang = int( ( -3. + sqrt( 1 + 8.*dof ) ) / 2. )
  end function

  function DOFtriangR(deg)
    real :: DOFtriangR
    real :: deg

    DOFtriangR = (deg + 1) *(deg+2)/ 2.
  end function DOFtriangR


  !> Initiate Mblock, allocation of the memory and setting to zero
  subroutine InitMblock(A, nr, nc)
    type(Mblock), intent(inout) :: A
    integer, intent(in) :: nr, nc

    allocate(A%Mb(nr, nc), source = 0.0)
!    A%Mb(:,:) = 0.
  end subroutine InitMblock

    !> Initiate Mblock, allocation of the memory and setting from old Mblock
  subroutine CopyMblock(A, B)
    class(Mblock), intent(inout) :: A
    type(Mblock), intent(in) :: B
    integer:: nr, nc

    if (.not. allocated(B%Mb) ) &
      stop 'CopyMblock - the source matrix is not Allocated'

    if ( allocated(A%Mb) ) &
      deallocate( A%Mb )

    nr = size(B%Mb(:,1))
    nc = size(B%Mb(1,:))

    allocate(A%Mb(nr, nc), source = B%Mb(1:nr,1:nc) )
  end subroutine CopyMblock


  !> Initiate Mblock, allocation of the memory and setting to zero
  subroutine DeleteMblock(A)
    type(Mblock), intent(inout) :: A

    deallocate(A%Mb)

  end subroutine DeleteMblock

  !> print Mblock to the screene
  subroutine WriteMblock( A)
    type(Mblock), intent(in) :: A
    integer:: i , ifile

    ifile = 10
    open(ifile, file='Mblock',status='UNKNOWN', position='append')

    write(ifile, '(a28,i5,a3,i5)')'------- Matrix block -------', &
         size(A%Mb,1),' x ',size(A%Mb,2)
    write(ifile, '(a2,20i11)') 'c:',1,2,3,4,5,6,7,8 !,9,10,11,12,13,14,15,16,17,18,19,20
!    write(ifile, '(a2,20i11)') 'c:',13,14,15,16,17,18,19,20

    do i=1, size(A%Mb,1)
    !do i=13, size(A%Mb,1)
    !do i=1, 3
       !write(ifile, '(40es9.1)') A%Mb(i,:)
       !write(ifile, '(i2,a2,40es11.3)') i,': ',A%Mb(i,13:)
       write(ifile, '(i2,a2,200es11.3)') i,': ',A%Mb(i,:)
    enddo
    write(ifile, *)' End of matrix block'
    close(ifile )

  end subroutine WriteMblock

  !> print Mblock to the screene
  subroutine WriteMblock_Screene( A)
    type(Mblock), intent(in) :: A
    integer:: i , ifile

    !ifile = 10
    !open(ifile, file='Mblock',status='UNKNOWN', position='append')

    write(*, '(a28,i5,a3,i5)')'------- Matrix block -------', &
         size(A%Mb,1),' x ',size(A%Mb,2)
    write(*, '(a2,20i11)') 'c:',1,2,3,4,5,6,7,8 !,9,10,11,12,13,14,15,16,17,18,19,20
!    write(*, '(a2,20i11)') 'c:',13,14,15,16,17,18,19,20

    do i=1, size(A%Mb,1)
    !do i=13, size(A%Mb,1)
    !do i=1, 3
       !write(*, '(40es9.1)') A%Mb(i,:)
       !write(*, '(i2,a2,40es11.3)') i,': ',A%Mb(i,13:)
       write(*, '(i2,a2,200es11.3)') i,': ',A%Mb(i,:)
    enddo
    write(*, *)' End of matrix block'

  end subroutine WriteMblock_Screene

  !> allocate block \f$ c\f$ and  \f$ c= a+b\f$, where \f$ a,b \f$ are blocks
  function addMblock(a,b) result(c)
    type(Mblock),intent(in):: a,b
    type(Mblock):: c

    allocate(c%Mb(size(a%Mb,1),size(a%Mb,2)))

    c%Mb = a%Mb + b%Mb

  end function addMblock

  !> allocate block \f$ c\f$ and  \f$ c= a-b\f$, where \f$ a,b \f$ are blocks
  function subMblock(a,b) result(c)
    type(Mblock),intent(in):: a,b
    type(Mblock):: c

    allocate(c%Mb(size(a%Mb,1),size(a%Mb,2)))

    c%Mb = a%Mb - b%Mb

  end function subMblock

  !> allocate block \f$ c\f$ and  \f$ c= a*b\f$, where \f$ a \f$ is real number
  !> and \f$ b\f$ is a block
  function mulMblock(a,b) result(c)
    real,intent(in):: a
    type(Mblock),intent(in):: b
    type(Mblock):: c

    allocate(c%Mb(size(b%Mb,1),size(b%Mb,2)))

    c%Mb = a * b%Mb

  end function

  !> for j-th vertex returns the index of Lagrangian nodes of degree Rdeg
  function VertexIndex(Rdeg, j)
    integer :: VertexIndex
    integer :: Rdeg, j

    select case(j)
    case(1)
       VertexIndex = 1
    case(2)
       VertexIndex = Rdeg + 1
    case(3)
       VertexIndex = (Rdeg+1)*(Rdeg+2)/2
    case(4:)
       print*,'Trouble in VertexIndex'
       stop
    end select

  end function VertexIndex

  !> for l-th node on the j-th edge returns the index of Lagrangian nodes of degree Rdeg
  function EdgeIndex(Rdeg, j, l)
    integer :: EdgeIndex
    integer :: Rdeg, j, l

    !HERE

    select case(j)
    case(1)
       EdgeIndex = l + 1
    case(2)
       EdgeIndex = (2*(Rdeg + 1) -l)*(l+1)/2
    case(3)
       EdgeIndex = (2*Rdeg +3 - (Rdeg - l) )* (Rdeg -l)/2   + 1
    case(4:)
       print*,'Trouble in EdgeIndex'
       stop
    end select

  end function EdgeIndex



 subroutine GetFaceIndexes(faces, nodes, j, flen)
    integer, dimension(:,:), intent(in) :: faces
    integer, dimension(:), intent(out) :: nodes
    integer, intent(in) :: flen
    integer, intent(in) :: j

    if (nbDim == 2) then
     nodes(1)= faces(idx,j)
     nodes(2)= faces(idx,mod(j,flen)+1)
    else
      if(flen .ne. 4) then !tetrahedra
        print*,'Only tetrahedra implmented in GetFaceIndexes'
        stop
      endif
        if (j == 1) then
          nodes(1)= faces(idx,2)
          nodes(2)= faces(idx,3)
          nodes(3)= faces(idx,4)
        elseif (j == 2) then
          nodes(1)= faces(idx,3)
          nodes(2)= faces(idx,1)
          nodes(3)= faces(idx,4)
        elseif (j == 3) then
          nodes(1)= faces(idx,1)
          nodes(2)= faces(idx,2)
          nodes(3)= faces(idx,4)
        else
          nodes(1)= faces(idx,3)
          nodes(2)= faces(idx,2)
          nodes(3)= faces(idx,1)
        endif
    endif

  end subroutine GetFaceIndexes


  subroutine WriteLocalMatrix( A, chfile)
    real, dimension(:,:), intent(in) :: A
    character(*) :: chfile
    integer:: i , ifile

    ifile = 10

    print*, ' WriteLocalMatrix ->', chfile

    open(ifile, file=chfile,status='UNKNOWN', position='append')

    write(ifile, '(a28,i5,a3,i5)')'------- Matrix block -------', &
         size(A,1),' x ',size(A,2)
    write(ifile, '(a2,20i11)') 'c:',1,2,3,4,5,6,7,8 !,9,10,11,12,13,14,15,16,17,18,19,20
!    write(ifile, '(a2,20i11)') 'c:',13,14,15,16,17,18,19,20

    do i=1, size(A,1)
       write(ifile, '(i2,a2,200es16.8)') i,': ',A(i,:)
    enddo
    write(ifile, *)' End of matrix block'
    close(ifile )

  end subroutine WriteLocalMatrix


end module geometry



