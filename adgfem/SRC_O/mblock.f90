module mblock_mod
   use paramets
   use geometry

implicit none

   type, public :: MblockST
      integer :: ndim ! ndim
!      integer(1:2, 1:2) :: dof  ! line dof - starting and ending index, column dof - starting and ending index
!      integer(1:2, 1:2) :: Tdof ! line and column Tdof
      integer, dimension(1:2) :: dof ! line and column dof
      integer, dimension(1:2) :: Tdof ! line and column dof
      real, allocatable, dimension (:, :)  :: Mb  ! block of matrix

   contains

      procedure :: init => InitMblockSquare
      procedure :: delete => DeleteMblockST
      procedure :: write => WriteMblockST

      procedure :: GetColumnIndex
      procedure :: GetRowIndex
      procedure :: CutFromMblockST
      !procedure :: copy => CopyMblockST
      !  procedure :: WriteMblock_Screene


   end type MblockST

   interface assignment(=)
      procedure :: CopyMblockST
   end interface

   public :: assignment(=)

!   interface operator(+)
!         procedure :: addMblockST
!   end interface
!   interface operator(-)
!         procedure :: subMblockST
!   end interface
!   interface operator(*)
!         procedure :: mulMblockST
!   end interface
!   public :: operator(+)
!   public :: operator(-)
!   public :: operator(*)


contains

!  !> allocate block \f$ c\f$ and  \f$ c= a+b\f$, where \f$ a,b \f$ are blocks
!  function addMblockST(a,b) result(c)
!    type(MblockST),intent(in):: a,b
!    type(MblockST):: c
!
!    allocate( MblockST :: c )
!    c%ndim = a%ndim
!    c%Tdof(1:2) = a%Tdof(1:2)
!    c%dof(1:2) = a%dof(1:2)
!    allocate(c%Mb(size(a%Mb,1),size(a%Mb,2)))
!
!    c%Mb = a%Mb + b%Mb
!
!  end function addMblockST


!> Initiate Mblock, allocation of the memory and setting to zero
!> SQUARE
   subroutine InitMblockSquare( this, nd, dof, Tdof)
      class(MblockST), intent(inout) :: this
      integer, intent(in) :: nd
      integer, intent(in) :: dof
      integer, intent(in) :: Tdof

      integer :: nsize

      nsize = nd*dof*Tdof
      this%ndim = nd
      this%dof(1:2) = dof
      this%Tdof(1:2) = Tdof

      allocate( this%Mb(nsize, nsize), source = 0.0 )

  end subroutine InitMblockSquare

!> Initiate Mblock, allocation of the memory and setting to zero
!> Rectangle
   subroutine InitMblockRectangle( this, nd, dof, Tdof)
      class(MblockST), intent(inout) :: this
      integer, intent(in) :: nd
      integer, dimension(1:2), intent(in) :: dof
      integer, dimension(1:2), intent(in) :: Tdof

      integer :: nr, nc

      nr = nd*dof(1)*Tdof(1)
      nc = nd*dof(2)*Tdof(2)

      this%ndim = nd
      this%dof(1:2) = dof(1:2)
      this%Tdof(1:2) = Tdof(1:2)

      allocate( this%Mb(nr, nc), source = 0.0 )

  end subroutine InitMblockRectangle

    !> Initiate Mblock, allocation of the memory and setting from old MblockST
  subroutine CopyMblockST( this, B)
    class(MblockST), intent(inout) :: this
    type(MblockST), intent(in) :: B
    integer :: nr, nc

    if (.not. allocated(B%Mb) ) &
      stop 'CopyMblock - the source matrix is not Allocated'

    if ( allocated(this%Mb) ) &
      deallocate( this%Mb )

    this%ndim = B%ndim
    this%dof(1:2) = B%dof(1:2)
    this%Tdof(1:2) = B%Tdof(1:2)


    nr = this%ndim * this%dof(1) * this%Tdof(1)
    nc = this%ndim * this%dof(2) * this%Tdof(2)

    allocate( this%Mb(nr, nc), source = B%Mb(1:nr,1:nc) )

  end subroutine CopyMblockST

   !> copy data from the older structure Mblock to MblockST
   subroutine CopyFromMblock( this, B, nd, dof, Tdof)
      class( MblockST ), intent(inout) :: this
      type( Mblock ), intent(in) :: B
      integer, intent(in) :: nd
      integer, dimension(1:2), intent(in) :: dof
      integer, dimension(1:2), intent(in) :: Tdof
      integer :: nr, nc

      if (.not. allocated(B%Mb) ) &
         stop 'CopyFromMblock - the source matrix is not Allocated'

      if ( allocated(this%Mb) ) &
         deallocate( this%Mb )

      this%ndim = nd
      this%dof(1:2) = dof(1:2)
      this%Tdof(1:2) = Tdof(1:2)

      nr = this%ndim * this%dof(1) * this%Tdof(1)
      nc = this%ndim * this%dof(2) * this%Tdof(2)

      if ( nr /= size( B%Mb(:,1) ) .or. nc /= size( B%Mb(1,:) ) ) &
         stop 'Wrong input arguments in CopyFromMblock!'

      allocate( this%Mb(1:nr,1:nc), source = B%Mb(1:nr,1:nc) )

   end subroutine CopyFromMblock

   !> copy data from the older structure Mblock to MblockST
   function CopyAllToMblock( BST) result( A )
      type( Mblock ) :: A
      type( MblockST ), intent(in) :: BST
      integer :: nr, nc

      if (.not. allocated(BST%Mb) ) &
         stop 'CopyToMblock - the source matrix is not Allocated'

      if ( allocated(A%Mb) ) &
         deallocate( A%Mb )
      nr = BST%ndim * BST%dof(1) * BST%Tdof(1)
      nc = BST%ndim * BST%dof(2) * BST%Tdof(2)

      allocate( A%Mb(1:nr,1:nc), source = BST%Mb(1:nr,1:nc) )

   end function CopyAllToMblock


    !> Cut part of MblockST
    !> nd = 1,2,3,4 or 0 - ALL
  subroutine CutFromMblockST( this, B, nd, dofFirst, dofLast, TdofFirst, TdofLast )
    class(MblockST), intent(inout) :: this
    type(MblockST), intent(in) :: B
    integer, intent(in) :: nd
    integer, dimension(1:2), intent(in) :: dofFirst
    integer, dimension(1:2), intent(in) :: dofLast
    integer, dimension(1:2), intent(in) :: TdofFirst
    integer, dimension(1:2), intent(in) :: TdofLast

    integer :: nr, nc, i, j , k
    integer :: A_i, A_j, B_i, B_j

    if (.not. allocated(B%Mb) ) &
      stop 'CutMblockST - the source matrix is not Allocated'

    if ( allocated(this%Mb) ) &
      deallocate( this%Mb )

    if (nd == 0)  then
      this%ndim = B%ndim
    else
      this%ndim = 1
    endif

    this%Tdof(1:2) = TdofLast(1:2) - TdofFirst(1:2) + 1
    this%dof(1:2) = dofLast(1:2) - dofFirst(1:2) + 1

    nr = this%ndim * this%Tdof(1)*this%dof(1)
    nc = this%ndim * this%Tdof(2)*this%dof(2)

    allocate( this%Mb(nr, nc) )

    ! all dimensions
    if (nd == 0) then

       do k = 1, this%ndim
         !rows
          do i = 1, this%Tdof(1)
            ! first row index of the current block
            B_i = B%GetRowIndex( k, dofFirst(1) , TdofFirst(1)+i-1 )
            A_i = this%GetRowIndex( k, 1, i )
            !columns
            do j = 1, this%Tdof(2)
               B_j = B%GetColumnIndex( k, dofFirst(2), TdofFirst(2)+j-1 )
               A_j = this%GetColumnIndex( k , 1, j)
               this%Mb( A_i : A_i+this%dof(1)-1 , A_j:A_j+ this%dof(2)-1 ) = &
                  B%Mb(  B_i: B_i + this%dof(1)-1, B_j:B_j+this%dof(2)-1 )
            end do !j
          end do !i
       end do !k

    ! dimension on index nd in matrix B
    else ! one dimension only

      !rows
       do i = 1,this%Tdof(1)
         ! first row index of the current block
         B_i = B%GetRowIndex( nd, dofFirst(1) , TdofFirst(1)+i-1 )
         A_i = this%GetRowIndex( 1, 1, i )
         !columns
         do j = 1, this%Tdof(2)
            B_j = B%GetColumnIndex( nd, dofFirst(2), TdofFirst(2)+j-1 )
            A_j = this%GetColumnIndex( 1 , 1, j)
            this%Mb( A_i : A_i+this%dof(1)-1 , A_j:A_j+ this%dof(2)-1 ) = &
               B%Mb(  B_i: B_i + this%dof(1)-1, B_j:B_j+this%dof(2)-1 )
         end do !j
       end do !i

    endif !nd


  end subroutine CutFromMblockST

   ! for given indices ndim, dof, tdof computes the corresponding row index in the matrix
  function GetRowIndex( this, nd, d , td ) result( row_index )
   class( MblockST ), intent(in) :: this
   integer, intent(in) :: nd, td, d

   integer :: row_index

   row_index = (nd-1) * this%dof(1) * this%Tdof(1) + (td-1) * this%dof(1) + d

  end function GetRowIndex

  ! for given indices ndim, dof, tdof computes the corresponding column index in the matrix
  function GetColumnIndex( this, nd, d , td ) result( col_index )
   class( MblockST ), intent(in) :: this
   integer, intent(in) :: nd, td, d

   integer :: col_index

   col_index = (nd-1) * this%dof(2) * this%Tdof(2) + (td-1) * this%dof(2) + d

  end function GetColumnIndex


  !> Initiate Mblock, allocation of the memory and setting to zero
  subroutine DeleteMblockST(A)
    class(MblockST), intent(inout) :: A

!    A%ndim = 0
!    A%dof = 0
!    A%Tdof = 0

    deallocate(A%Mb)

  end subroutine DeleteMblockST

   !> print Mblock to the screene
   subroutine WriteMblockST( A)
       class(MblockST), intent(in) :: A
       integer:: i , ifile

       ifile = 10
       open(ifile, file='MblockST.txt',status='UNKNOWN') ! , position='append')

       write(ifile, '(a28,i5,a3,i5)')'------- Matrix block -------', &
            size(A%Mb,1),' x ',size(A%Mb,2)

       write(ifile, '(i5,a3,i5,a3,i5,a3,i5,a3,i5)') &
            A%ndim, ' , ', A%dof(1), ' x ', A%dof(2), ' , ', A%Tdof(1), ' x ', A%Tdof(2)

       write(ifile, '(a2,20i11)') 'c:',1,2,3,4,5,6,7,8 !,9,10,11,12,13,14,15,16,17,18,19,20

       do i=1, size(A%Mb,1)
       !do i=13, size(A%Mb,1)
       !do i=1, 3
          !write(ifile, '(40es9.1)') A%Mb(i,:)
          !write(ifile, '(i2,a2,40es11.3)') i,': ',A%Mb(i,13:)
          write(ifile, '(i2,a2,200es9.1)') i,': ',A%Mb(i,:)
       enddo
       write(ifile, *)' End of matrix block'
       close(ifile )

       print*, 'MblockST written to the file MblockST.txt'

  end subroutine WriteMblockST


!> allocate block \f$ c\f$ and  \f$ c= a+b\f$, where \f$ a,b \f$ are blocks
!  function addMblockST(a,b) result(c)
!    type(MblockST),intent(in):: a,b
!    type(MblockST):: c
!
!    allocate(c%Mb(size(a%Mb,1),size(a%Mb,2)))
!
!    c%Mb = a%Mb + b%Mb
!
!  end function addMblockST
!
!  !> allocate block \f$ c\f$ and  \f$ c= a-b\f$, where \f$ a,b \f$ are blocks
!  function subMblockST(a,b) result(c)
!    type(MblockST),intent(in):: a,b
!    type(MblockST):: c
!
!    allocate(c%Mb(size(a%Mb,1),size(a%Mb,2)))
!
!    c%Mb = a%Mb - b%Mb
!
!  end function subMblockST
!
!  !> allocate block \f$ c\f$ and  \f$ c= a*b\f$, where \f$ a \f$ is real number
!  !> and \f$ b\f$ is a block
!  function mulMblockST(a,b) result(c)
!    real,intent(in):: a
!    type(MblockST),intent(in):: b
!    type(MblockST):: c
!
!    allocate(c%Mb(size(b%Mb,1),size(b%Mb,2)))
!
!    c%Mb = a * b%Mb
!
!  end function mulMblockST



end module mblock_mod


! TEST in o_main.f90

!   integer :: i,j, td, d, nd,jj
!   type(MblockST), allocatable :: AA, BB
!   type(Mblock), allocatable :: MM
!
!
!  !open(debug, file = debug_file, action="write", status="replace" )
!  open(debug, file = debug_file, action="write", status="UNKNOWN" )
!
!  call cpu_time(state%start_time)
!
!!  print*, 'InitElementDof - Tdof_plus = Tdof'
!
!  allocate( MblockST :: AA, BB )
!  nd =2
!  td = 3
!  d =2
!
!  call AA%init( nd, d, td )
!
!  do i = 1, nd
!   do j = 1,td
!      do jj =1,td
!         AA%Mb( (i-1)*d*td+ (j-1)*d + 1 : (i-1)*d*td+ j*d , &
!         (i-1)*d*td+ (jj-1)*d + 1 : (i-1)*d*td+ jj*d) = i*j*jj
!      print*, (i-1)*d*td+ (j-1)*d + 1 , (i-1)*d*td+ j*d
!      end do
!   end do
!  enddo
!
!  call BB%CutFromMblockST(AA, 1, (/1,1/),(/2,2/),(/1,1/),(/3,3/) )
!
!  allocate( MM )
!  MM = CopyAllToMblock(BB)
!
!  BB = AA
!
!  call WriteMblock(MM)
!
!
!  call BB%write()
!  deallocate( AA, BB, MM )
!  stop 'FILIP'
