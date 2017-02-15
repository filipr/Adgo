module mtxform

  use main_data
  use data_mod
  use paramets
  use geometry
  use element_mod

  use matrix_oper_int

IMPLICIT NONE

  ! fcn related to matrix
  public:: DG2HB			! DG into Harwel-Boeing format (compressed colum storage)
  private:: SORT1, PIKSRT, TWRIT 	! used by DG2HB
  public:: DG2CROWS 			! --- NOT TESTED YET --- !

  public:: CompareMatrix 		! compare element of DG and matrix stored in HB format
  public:: DajHBij			! get (i,j)-element
  public:: GetHBij                      ! get (i,j)-element

  public:: GetDGij

CONTAINS


  subroutine SORT1(A,N,iout)
  integer,dimension(:,:),intent(inout) ::    A     !Table to be sorted
  integer,intent(in) :: N,iout

  if (iout.eq.1) then
    print *,' '
    print *,'Table to be sorted:'
    call TWRIT(N,A)
  endif

  call PIKSRT(N,A)

  if (iout.eq.1) then
    print *,' '
    print *,'Sorted table (straight insertion):'
    call TWRIT(N,A)
  endif

  end subroutine SORT1

  ! Sorting an array with straight insertion method *
  subroutine PIKSRT(N,ARR)
  integer,intent(in) :: N
  integer,intent(inout) :: ARR(3,N)
  integer :: i,j,a(3)

  do j=2,N
    a = ARR(1:3,j)
    do i=j-1,1,-1
      if ( ARR(1,i)<=a(1)) goto 10
      ARR(1:3,i+1)=ARR(1:3,i)
    end do
    i=0
10  ARR(1:3,i+1 )= a
  end do
  return

  end subroutine PIKSRT


  !write table of size N to standard output
  subroutine TWRIT(N,ARR)
  integer::N
  integer :: ARR(3,N)
  integer :: i

  write(*,"(I6,I6,I6)") (ARR(1,i),ARR(2,i),ARR(3,i),i=1,N)
  return
  end subroutine TWRIT


  !> Store DG-matrix in Harwell-Boeing (compressed column storage) format
  !> with 1-based indexing
  subroutine DG2HB(Ap,Ai,Ax)
  integer,dimension(:),allocatable,intent(inout) :: Ap, Ai
  real,dimension(:),allocatable,intent(inout) :: Ax
  integer :: i,j,k,shiftidx,last
  integer :: allsize, facesize, ptrsize, ptridx
  class(element), pointer :: elem, elem1
  integer :: dof,ndof,ndof1
  integer,dimension(:,:),allocatable :: array
  real,dimension(:,:),allocatable :: accum
  integer :: iii


  allocate( accum(maxval(grid%elem%MGdof)*ndim, maxval(grid%elem%MGdof)*ndim) )
  allocate( array(3, 1+maxval(grid%elem%flen) ) )

  if( 0.eq.1) then
    print*,'#elem',grid%nelem
    do i=1,grid%nelem,1
      print*,'elem',i
      do j=1,grid%elem(i)%flen,1
        k = grid%elem(i)%face(neigh,j)
        print*,'sused',k
      enddo
    enddo
  end if

  ! evaluate sizes of arrays Ap,Ai,Ax
  allsize = 0
  ptrsize = 1	! first after, ie. we have 0 columns

  do i=1,grid%nelem,1
    elem => grid%elem(i)
    ptrsize = ptrsize + ndim*elem%MGdof
    facesize = elem%MGdof
    do j=1,elem%flen,1
      k = elem%face(neigh,j)
      if (k>0) then
        facesize = facesize + grid%elem(k)%MGdof
      endif
    enddo
    allsize = allsize + ndim**2*elem%MGdof*facesize
  enddo

  ! allocate arrays
  allocate(Ap(ptrsize), Ai(allsize), Ax(allsize))

  Ap=0
  Ai=0
  Ax=0.

  !print*,ptrsize-1,'?=',state%nsize

  ! print*,'fill arrays with data'
  shiftidx = 1
  ptridx = 0

  do i=1,grid%nelem,1
    elem => grid%elem(i)
    dof = elem%MGdof
    ndof = ndim*dof

    ! Array "array" is used to get information about column structure, ie.
    ! which block in column is above and/or bellow mtx diagonal
    ! array(1,:) stores information for faces
    array(1,:) = (/ i, elem%face(neigh, 1:elem%flen) /)
    ! array(2,:) stores information where corresponding flux matrix block
    ! is stored, ie. elem%block(*)
    array(2,:) = (/ (k, k=0,elem%flen,1) /)
    array(3,:) = (/ 0, elem%face(nei_i, 1:elem%flen) /)

    call SORT1(array(:,1:elem%flen+1),elem%flen+1,0)

    !print*,'set data for indices array Ai'
    last = 0;
    do j=1,elem%flen+1,1
      if ( array(1,j).gt.0 ) then
        iii = array(1,j)
        elem1 => grid%elem(iii)
        ndof1 = ndim*elem1%MGdof
        Ai(shiftidx+last : shiftidx+last+ndof1-1) = (/ (k, k=elem1%MGncv,elem1%MGncv+ndof1-1,1) /)
        last = last+ndof1
      endif
    enddo
    !print*,'copy Ai for all 2:ndof columns'
    do k=2,ndof,1
      Ai(shiftidx+(k-1)*last : shiftidx+k*last-1) = Ai(shiftidx : shiftidx+last-1)
    enddo

    ! evaluate diagonal matrix block
    accum(1:ndof,1:ndof) = elem%block(0)%Mb(1:ndof, 1:ndof)

    if ( eta.ne.0. ) then
      do k = 0,ndof-1,dof
        accum(k+1:k+dof,k+1:k+dof) = accum(k+1:k+dof,k+1:k+dof) + &
          eta * elem%Mass%Mb(1:dof, 1:dof)
      enddo
    endif

    !print*,'set data for values array Ax'
    do k=1,ndof,1
      ptridx = ptridx+1
      Ap(ptridx)=shiftidx

      do j=1,elem%flen+1,1

        if ( array(1,j).eq.i ) then ! equivalent (array(2,j)==0)
          !print*,'diagonal block'
          Ax(shiftidx:shiftidx+ndof-1) = accum(1:ndof,k)
          shiftidx = shiftidx+ndof
        elseif ( array(1,j).gt.0 ) then
          !print*,'offdiagonal block in column'
          elem1 => grid%elem( array(1,j) )
          ndof1 = ndim*elem1%MGdof
          Ax(shiftidx:shiftidx+ndof1-1) = grid%elem( array(1,j) )%block( array(3,j) )%Mb(1:ndof1,k)
          shiftidx = shiftidx+ndof1
        endif

      enddo !j=1,4,1
    enddo !k=1,ndof,1pMGLinSolverRecu4LINSOL
  enddo !i=1,grid%nelem,1

  Ap(ptrsize) = allsize+1

  end subroutine DG2HB


  !> Store DG-matrix in COMPRESSED ROW STORAGE format
  !> with 1-based indexing
  subroutine DG2CROWS(Ap,Ai,Ax)
  integer,dimension(:),allocatable,intent(inout) :: Ap, Ai
  real,dimension(:),allocatable,intent(inout) :: Ax
  integer :: i,j,k,shiftidx,last
  integer :: allsize, facesize, ptrsize, ptridx
  class(element), pointer :: elem, elem1
  integer :: dof,ndof,ndof1
  integer,dimension(:,:),allocatable :: array
  real,dimension(:,:),allocatable :: accum
  integer :: iii


  allocate( accum(maxval(grid%elem%MGdof)*ndim, maxval(grid%elem%MGdof)*ndim) )
  allocate( array(3, 1+maxval(grid%elem%flen) ) )

  if( 0.eq.1) then
    print*,'#elem',grid%nelem
    do i=1,grid%nelem,1
      print*,'elem',i
      do j=1,grid%elem(i)%flen,1
        k = grid%elem(i)%face(neigh,j)
        print*,'sused',k
      enddo
    enddo
  end if

  ! evaluate sizes of arrays Ap,Ai,Ax
  allsize = 0
  ptrsize = 1	! first after, ie. we have 0 columns

  do i=1,grid%nelem,1
    elem => grid%elem(i)
    ptrsize = ptrsize + ndim*elem%MGdof
    facesize = elem%MGdof
    do j=1,elem%flen,1
      k = elem%face(neigh,j)
      if (k>0) then
        facesize = facesize + grid%elem(k)%MGdof
      endif
    enddo
    allsize = allsize + ndim**2*elem%MGdof*facesize
  enddo

  ! allocate arrays
  allocate(Ap(ptrsize), Ai(allsize), Ax(allsize))

  Ap=0
  Ai=0
  Ax=0.

  !print*,ptrsize-1,'?=',state%nsize

  ! print*,'fill arrays with data'
  shiftidx = 1
  ptridx = 0

  do i=1,grid%nelem,1
    elem => grid%elem(i)
    dof = elem%MGdof
    ndof = ndim*dof

    ! Array "array" is used to get information about column structure, ie.
    ! which block in column is above and/or bellow mtx diagonal
    ! array(1,:) stores information for faces
    array(1,:) = (/ i, elem%face(neigh, 1:elem%flen) /)
    ! array(2,:) stores information where corresponding flux matrix block
    ! is stored, ie. elem%block(*)
    array(2,:) = (/ (k, k=0,elem%flen,1) /)
    array(3,:) = (/ 0, elem%face(nei_i, 1:elem%flen) /)

    call SORT1(array(:,1:elem%flen+1),elem%flen+1,0)

    !print*,'set data for indices array Ai'
    last = 0;
    do j=1,elem%flen+1,1
      if ( array(1,j).gt.0 ) then
        iii = array(1,j)
        elem1 => grid%elem(iii)
        ndof1 = ndim*elem1%MGdof
        Ai(shiftidx+last : shiftidx+last+ndof1-1) = (/ (k, k=elem1%MGncv,elem1%MGncv+ndof1-1,1) /)
        last = last+ndof1
      endif
    enddo
    !print*,'copy Ai for all 2:ndof columns'
    do k=2,ndof,1
      Ai(shiftidx+(k-1)*last : shiftidx+k*last-1) = Ai(shiftidx : shiftidx+last-1)
    enddo

    ! evaluate diagonal matrix block
    accum(1:ndof,1:ndof) = elem%block(0)%Mb(1:ndof, 1:ndof)

    if ( eta.ne.0. ) then
      do k = 0,ndof-1,dof
        accum(k+1:k+dof,k+1:k+dof) = accum(k+1:k+dof,k+1:k+dof) + &
          eta * elem%Mass%Mb(1:dof, 1:dof)
      enddo
    endif

    !print*,'set data for values array Ax'
    do k=1,ndof,1
      ptridx = ptridx+1
      Ap(ptridx)=shiftidx

      do j=1,elem%flen+1,1


        if ( array(1,j).eq.i ) then ! equivalent (array(2,j)==0)
          !print*,'diagonal block'
          Ax(shiftidx:shiftidx+ndof-1) = accum(k,1:ndof)
          shiftidx = shiftidx+ndof
        elseif ( array(1,j).gt.0 ) then
          !print*,'offdiagonal block in column'
          elem1 => grid%elem( array(1,j) )
          ndof1 = ndim*elem1%MGdof
          Ax(shiftidx:shiftidx+ndof1-1) = grid%elem( array(1,j) )%block( array(2,j) )%Mb(k,1:ndof1)
          shiftidx = shiftidx+ndof1
        endif

      enddo !j=1,4,1 pMGLinSolverRecu4LINSOL
    enddo !k=1,ndof,1
  enddo !i=1,grid%nelem,1

  Ap(ptrsize) = allsize+1

  end subroutine DG2CROWS






! Get element of DG matrix on position [i,j]
  !!!!!! wihtout Mass matrix !!!!!!
  function GetDGij(i,j)
    real :: GetDGij
    integer, intent(in) :: i,j
    class(element), pointer :: elemi,elemj
    integer :: cumsum,k

    if( (min(i,j)<1).or.(max(i,j)>state%nsize) ) then
      print*,'Wrong index - GetDJij'
      STOP
    end if

    elemi => null()
    elemj => null()

    ! find element for row "i"
    ! find element for column "j"
    cumsum = 0
    do k=1,grid%nelem,1

      cumsum = cumsum + ndim*grid%elem(k)%dof

      if( .not.associated(elemi) ) then
        if( cumsum.ge.i) then
          elemi => grid%elem(k)
          !print*,'I: elem,i,j',k,i,j
          if( associated(elemj) ) EXIT
        end if
      end if

      if( .not.associated(elemj) ) then
        if( cumsum.ge.j) then
          elemj => grid%elem(k)
          !print*,'J: elem,i,j',k,i,j
          if( associated(elemi) ) EXIT
        end if
      end if

    end do

    if( elemi%i.eq.elemj%i) then ! diagonal matrix block
        ! mass matrix

        ! flux matrix
        GetDGij = elemi%block(0)%Mb( 1+i-elemi%ncv , 1+j-elemj%ncv )
    elseif( any( elemi%face(neigh,1:elemi%flen).eq.elemj%i ) ) then ! nonzero off-diagonal matrix block
        do k=1,elemi%flen,1
          if( elemi%face(neigh,k).eq.elemj%i) then
            EXIT
          end if
        end do
        GetDGij = elemi%block(k)%Mb( 1+i-elemi%ncv , 1+j-elemj%ncv )
    else ! zero matrix block
      GetDGij = 0
    end if

  end function GetDGij



subroutine CompareMatrix
  integer :: i,j
  integer,allocatable,dimension(:) :: Ap
  integer,allocatable,dimension(:) :: Ai
  real,allocatable,dimension(:) :: Ax
  integer :: N


  grid%elem(:)%MGdeg = grid%elem(:)%deg
  grid%elem(:)%MGdof = grid%elem(:)%dof
  grid%elem(:)%MGncv = grid%elem(:)%ncv

  call DG2hb(Ap,Ai,Ax)
print*,'@@ eta',eta,state%time%tau(1)

print*,'Ax',Ax(1:3)
print*,'C=',grid%elem(1)%block(0)%Mb(1:3,1)
print*,'M=',grid%elem(1)%Mass%Mb(1:3,1), maxval(  grid%elem(1)%Mass%Mb )
pause
    do i=1,state%nsize,1
      do j=1,state%nsize,1
        print*,i,j,state%nsize
        if( GetDGij(i,j).eq.GetHBij(Ap,Ai,Ax,i,j) ) then

        else
          print*,'PROBLEM',i,j,GetDGij(i,j),GetHBij(Ap,Ai,Ax,i,j)
          PAUSE
        end if
      end do
    end do
  end subroutine CompareMatrix


  function DajHBij(i,j)
  real :: DajHBij
  integer,intent(in) :: i,j
  integer,allocatable,dimension(:) :: Ap
  integer,allocatable,dimension(:) :: Ai
  real,allocatable,dimension(:) :: Ax
  integer :: N

  grid%elem(:)%MGdeg = grid%elem(:)%deg
  grid%elem(:)%MGdof = grid%elem(:)%dof
  grid%elem(:)%MGncv = grid%elem(:)%ncv

  call DG2HB(Ap,Ai,Ax)

  DajHBij = GetHBij(Ap,Ai,Ax,i,j)

  end function

  ! Get element of HB matrix on position [i,j]
  ! HB = Harwell-Boeing matrix is described by arrays: Ap,Ai,Ax
  function GetHBij(Ap,Ai,Ax,i,j)
    real :: GetHBij
    integer,dimension(:),intent(in) :: Ap,Ai
    real,dimension(:),intent(in) :: Ax
    integer,intent(in) :: i,j
    integer :: k, nzidx

    nzidx = -1

!    print*,'i=',i
!    print*,Ai(Ap(j):Ap(j+1)-1)


    do k=Ap(j),Ap(j+1)-1,1
      if( Ai(k).eq.i) then
        nzidx = k
        EXIT
      end if
    end do

    if( nzidx.gt.-1) then
      GetHBij = Ax( nzidx )
    else
      GetHBij = 0
    end if

  end function GetHbij

end module mtxform
