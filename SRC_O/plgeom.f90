!> volume (triangular, quadrilateral) and edge quadrature rules
!>
!> contain also the values of the basis functions (phi) and their derivatives
!> (Dphi) in the corresponding integration nodes,

module plot_geom

  implicit none


  real :: rho_infty, v_infty, p_infty, alpha_infty, theta_infty

  public :: SetSubTri
  public :: ComputeQuantity

  public :: max_deg 
  public :: Plot3D
  public :: PlotIsolines
  public :: PlotHPColorMesh
  public :: SetHPmeshFileName

  integer, parameter ::  max_deg  = 10
contains
  
subroutine SetSubTri(max_deg, max_dof, subtri, lambda, subedge)
  integer, intent(in) :: max_deg, max_dof
  integer,dimension(0:max_deg, 1: max_deg**2, 1:3), intent(inout) :: subtri
  integer,dimension(0:max_deg, 1:3, 1: max_deg+1), intent(inout) :: subedge
  real,dimension(0:max_deg, 1: max_dof, 1:3), intent(inout) :: lambda
  integer :: i,j,k, ie, deg, ideg

  do ideg = 0, max_deg
     deg = max(1, ideg)
     ! barycentric coordinates of nodes
     k = 0
     do i=0, deg
        do j=0, deg - i
           k = k + 1
           lambda(ideg, k, 2:3) = (/ 1.*j/deg, 1.*i/deg /)
           lambda(ideg, k,   1) = 1. - sum(lambda(deg, k, 2:3) )
           !write(*,'(a5,4i5,3f8.4)') 'lam',ideg,k,i,j, lambda(ideg, k, 1:3)

           if(i == 0) subedge(ideg, 1, j+1) = k    ! subedges
           if(j == deg-i) subedge(ideg, 2, i+1) = k    ! subedges
           if(j == 0) subedge(ideg, 3, deg+1-i) = k    ! subedges
        enddo
     enddo

     !write(*,'(a5,2i5,a1,8i3)') 'edge',ideg,1,'|', subedge(ideg, 1, 1:deg+1)
     !write(*,'(a5,2i5,a1,8i3)') 'edge',ideg,2,'|', subedge(ideg, 2, 1:deg+1)
     !write(*,'(a5,2i5,a1,8i3)') 'edge',ideg,3,'|', subedge(ideg, 3, 1:deg+1)

     
     ! subtriangles
     ie = 0
     k = 0
     do i=0, deg
        do j=0, deg - i
           k = k + 1
           if( j < deg - i) then
              ie = ie + 1
              subtri(ideg, ie, 1) = k 
              subtri(ideg, ie, 2) = k + 1
              subtri(ideg, ie, 3) = k + deg - i + 1

              !write(*,'(a5,4i5,a1,5i3)') 'tri',ideg,k,i,j,'|', subtri(ideg, ie, 1:3)

              if(j < deg - i - 1 .and. i < deg - 1) then
                 ie = ie + 1
                 subtri(ideg, ie, 1) = k + 1
                 subtri(ideg, ie, 2) = k + deg - i + 2
                 subtri(ideg, ie, 3) = k + deg - i + 1
                 
                 !write(*,'(a5,4i5,a1,5i3)') '!!!',ideg,k,i,j,'|', subtri(ideg, ie, 1:3)
              endif
           endif

        enddo
     enddo
  enddo
end subroutine SetSubTri


subroutine ComputeQuantity(idof, ndim, w, q, quantity )
  integer, intent(in) :: idof, ndim
  character*5 quantity
  real, dimension(1:idof,1:ndim), intent(in) :: w
  real, dimension(1:idof), intent(out) :: q
  real :: rv2(1:100), p(1:100), kappa=1.4
  integer :: j,k, iqua

  if((quantity(1:1) == 'E' .or. quantity(1:1) == 'S' )  .and. len(quantity) > 0) then
     do j=2,len(quantity)
        if(quantity(j:j) /= ' ') k = j
     enddo
     if(k == 2) read(quantity(2:2),'(i1)') iqua
     if(k == 3) read(quantity(2:3),'(i2)') iqua

     !print*,'quantity = ',quantity, len(quantity), iqua
     !stop

     if(quantity(1:1) == 'S') then
        if(ndim < iqua) then
           q(1:idof) = 0.
        else
           q(1:idof) = w(1:idof,iqua)
        endif

     else !if(quantity(1:1) == 'E') then
        if(ndim < iqua) then
           q(1:idof) = 0.
        else
           q(1:idof) = log(max(1E-17, w(1:idof,iqua))) / log(10.)
        endif
     endif
     !print*,'quantity = ',quantity, len(quantity), iqua, q(1:3)

  ! density
  elseif(quantity .eq. 'RO' .or. quantity == 'S1'.or. quantity == 'Alg1') then
     q(1:idof) = w(1:idof,1)

  elseif(quantity .eq. 'E' .or. quantity == 'S4'.or. quantity == 'Alg4') then
     q(1:idof) = w(1:idof,4)

  elseif(quantity .eq. 'V') then
     q(1:idof) = (w(1:idof,2)**2 + w(1:idof,3)**2)**0.5/w(1:idof,1)

  elseif(quantity .eq. 'V1') then
     q(1:idof) = w(1:idof,2)/w(1:idof,1)

  elseif(quantity .eq. 'V2') then
     q(1:idof) = w(1:idof,3)/w(1:idof,1)

  elseif(quantity .eq. 'RV') then
     q(1:idof) = (w(1:idof,2)**2 + w(1:idof,3)**2)**0.5

  elseif(quantity .eq. 'RV1'  .or. quantity == 'S2'.or. quantity == 'Alg2') then
     q(1:idof) = w(1:idof,2)

  elseif(quantity .eq. 'RV2'  .or. quantity == 'S3'.or. quantity == 'Alg3') then
     q(1:idof) = w(1:idof,3)

  elseif(quantity .eq. 'RV2'  .or. quantity == 'S234') then
     q(1:idof) = w(1:idof,2) + w(1:idof,3) + w(1:idof,4)

  ! elseif(quantity == 'S5'.or. quantity == 'Alg5') then
  !    if(ndim < 5) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,5)
  !    endif

  ! elseif(quantity == 'S6' .or. quantity == 'Alg6' ) then
  !    if(ndim < 6) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,6)
  !    endif

  ! elseif(quantity == 'S7' .or. quantity == 'Alg7' ) then
  !    if(ndim < 7) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,7)
  !    endif

  ! elseif(quantity == 'E1' ) then
  !    if(ndim < 1) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = log( w(1:idof,1) ) / log(10.)
  !    endif

  ! elseif(quantity == 'E2' ) then
  !    if(ndim < 2) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = log( w(1:idof,2) ) / log(10.)
  !    endif

  ! elseif(quantity == 'E3' ) then
  !    if(ndim < 3) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = log( w(1:idof,3) ) / log(10.)
  !    endif

  ! elseif(quantity == 'E4' ) then
  !    if(ndim < 4) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = log( w(1:idof,4) ) / log(10.)
  !    endif

  ! elseif(quantity == 'E5' ) then
  !    if(ndim < 5) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = log( w(1:idof,5) ) / log(10.)
  !    endif

  ! elseif(quantity == 'E6' ) then
  !    if(ndim < 6) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = log( w(1:idof,6) ) / log(10.)
  !    endif

  ! elseif(quantity == 'E7' ) then
  !    if(ndim < 7) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = log( w(1:idof,7) ) / log(10.)
  !    endif

  ! elseif(quantity == 'E8' ) then
  !    if(ndim < 8) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = log( w(1:idof,8) ) / log(10.)
  !    endif

  ! elseif(quantity == 'E9' ) then
  !    if(ndim < 9) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = log( w(1:idof,9) ) / log(10.)
  !    endif

  ! elseif(quantity == 'E10' ) then
  !    if(ndim < 10) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = log( w(1:idof,10) ) / log(10.)
  !    endif

  ! elseif(quantity == 'E11' ) then
  !    if(ndim < 11) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = log( w(1:idof,11) ) / log(10.)
  !    endif

  ! elseif(quantity == 'E12' ) then
  !    if(ndim < 12) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = log( w(1:idof,12) ) / log(10.)
  !    endif

  ! elseif(quantity == 'S1'  ) then
  !    if(ndim < 1) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,1)
  !    endif

  ! elseif(quantity == 'S5'  ) then
  !    if(ndim < 5) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,5)
  !    endif

  ! elseif(quantity == 'S6'  ) then
  !    if(ndim < 6) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,5)
  !    endif

  ! elseif(quantity == 'S8'  ) then
  !    if(ndim < 8) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,8)
  !    endif

  ! elseif(quantity == 'S9'  ) then
  !    if(ndim < 9) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,9)
  !    endif

  ! elseif(quantity == 'S10'  ) then
  !    if(ndim < 10) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,10)
  !    endif
  ! elseif(quantity == 'S11'  ) then
  !    if(ndim < 11) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,11)
  !    endif
  ! elseif(quantity == 'S12'  ) then
  !    if(ndim < 12) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,12)
  !    endif

  ! elseif(quantity == 'S13'  ) then
  !    if(ndim < 13) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,13)
  !    endif

  ! elseif(quantity == 'S14'  ) then
  !    if(ndim < 14) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,14)
  !    endif

  ! elseif(quantity == 'S15'  ) then
  !    if(ndim < 15) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,15)
  !    endif

  ! elseif(quantity == 'S16'  ) then
  !    if(ndim < 16) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,16)
  !    endif

  ! elseif(quantity == 'S17'  ) then
  !    if(ndim < 17) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,17)
  !    endif

  ! elseif(quantity == 'S18'  ) then
  !    if(ndim < 18) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,18)
  !    endif

  ! elseif(quantity == 'S19'  ) then
  !    if(ndim < 19) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,19)
  !    endif

  ! elseif(quantity == 'S20'  ) then
  !    if(ndim < 20) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,20)
  !    endif

  ! elseif(quantity == 'S21'  ) then
  !    if(ndim < 21) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,21)
  !    endif

  ! elseif(quantity == 'S22'  ) then
  !    if(ndim < 22) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,22)
  !    endif

  ! elseif(quantity == 'S23'  ) then
  !    if(ndim < 23) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,23)
  !    endif

  ! elseif(quantity == 'S24'  ) then
  !    if(ndim < 24) then
  !       q(1:idof) = 0.
  !    else
  !       q(1:idof) = w(1:idof,24)
  !    endif

  elseif(quantity .eq. 'P') then
     rv2(1:idof) = (w(1:idof,2)**2 + w(1:idof,3)**2)/w(1:idof,1) 
     q(1:idof) = (kappa-1)* (w(1:idof,4) - rv2(1:idof)/2 )

  elseif(quantity .eq. 'PC') then
     rv2(1:idof) = (w(1:idof,2)**2 + w(1:idof,3)**2)/w(1:idof,1) 
     q(1:idof) = (kappa-1)* (w(1:idof,4) - rv2(1:idof)/2 )

     q(1:idof)  = (q(1:idof)  - p_infty)/(rho_infty * v_infty * v_infty/2)

  elseif(quantity .eq. 'M') then
     rv2(1:idof) = (w(1:idof,2)**2 + w(1:idof,3)**2)/w(1:idof,1)
     p(1:idof) = (kappa-1)* (w(1:idof,4) - rv2(1:idof)/2 )
     q(1:idof) = (rv2(1:idof) /p(1:idof) /kappa)**0.5

  elseif(quantity .eq. 'S') then
     rv2(1:idof) = (w(1:idof,2)**2 + w(1:idof,3)**2)/w(1:idof,1)
     p(1:idof) = (kappa-1)* (w(1:idof,4) - rv2(1:idof)/2 )
     q(1:idof) = log(p(1:idof) / w(1:idof,1)**kappa)
     !write(*,'(a4,10es12.4)') '????',q(1:idof)

     ! mass liquid fraction for wet_steam
  elseif(quantity .eq. 'w') then
     q(1:idof) = w(1:idof,5)/w(1:idof,1)

  else
     print*,'Unknown quantity ',quantity,' in ComputeQuantity'
     stop
  endif
end subroutine ComputeQuantity


  !> Setting of names for tri* sol* and commmand for plotting
  subroutine SetGnuFileName(inum, gnu_name)
    integer, intent (in) :: inum
    character(len=50),intent(inout) :: gnu_name
    character(len=2) :: ch2
    integer :: num_size, text_size,  file_size
    integer :: is
    
    text_size = 4
    num_size = 2
    file_size = text_size + num_size


    gnu_name = 'gnu.00'

    if(inum > 0) then
       is = int(log(1.*inum)/log(10.)) 
    else
       is = 0
    endif

    !print*,'!!!',inum,is, num_size+text_size-is, num_size+text_size, num_size-is, num_size

    write( ch2, '(i2)' ) inum  ! change the format if num_size /= 5 !!!
    gnu_name(num_size+text_size-is:num_size+text_size) = ch2(num_size-is: num_size)

    !print*,'######',gnu_name,'|',inum, ch2, is
    !print*,'######',sol_name,'|',inum, ch5, is
    
  end subroutine SetGnuFileName


  !> Setting of names for tri* sol* and commmand for plotting
  subroutine SetHPmeshFileName(inum, gnu_name)
    integer, intent (in) :: inum
    character(len=50),intent(inout) :: gnu_name
    character(len=2) :: ch2
    integer :: num_size, text_size, file_size
    integer :: is
    
    text_size = 4
    num_size = 2
    file_size = text_size + num_size


    gnu_name = 'hpP_00'

    if(inum > 0) then
       is = int(log(1.*inum)/log(10.)) 
    else
       is = 0
    endif

    !print*,'!!!',inum,is, num_size+text_size-is, num_size+text_size, num_size-is, num_size

    write( ch2, '(i2)' ) inum  ! change the format if num_size /= 5 !!!
    gnu_name(num_size+text_size-is:num_size+text_size) = ch2(num_size-is: num_size)

    !print*,'######',gnu_name,'|',inum, ch2, is
    !print*,'######',sol_name,'|',inum, ch5, is
    
  end subroutine SetHPmeshFileName

  subroutine Plot3D(ignu, nelem, npoin, max_deg, max_dof, &
              subtri, lambda, x, lnd, deg, q)
    integer, intent (in) :: ignu, nelem, npoin, max_deg, max_dof
    integer, dimension(0:max_deg, 1: max_deg**2, 1:3), intent(in) :: subtri
    real, dimension(0:max_deg, 1: max_dof, 1:3), intent(in) :: lambda
    real, dimension(1:npoin, 1:2), intent(in) :: x
    integer, dimension(1:nelem, 1:3), intent(in) :: lnd
    integer, dimension(1:nelem), intent(in) :: deg
    real, dimension(1:nelem, 1:max_dof), intent(in) :: q
    real :: xl(1:3, 1:2)
    integer :: i,k, l,  ideg, il(1:3)

    do i=1,nelem
       ideg = max(1, deg(i) )

       do k=1, ideg**2
          il(1:3) = subtri(ideg, k, 1:3)
          
          do l=1,3
             xl(l,1) = dot_product(lambda(ideg, il(l), 1:3), x(lnd(i, 1:3),1) )
             xl(l,2) = dot_product(lambda(ideg, il(l), 1:3), x(lnd(i, 1:3),2) )
          enddo
         
          do l=0,3
             if(deg(i) == 0) then
                write(ignu,'(3es14.6)') xl(mod(l,3) +1, 1:2), q(i,1)
             else
                write(ignu,'(3es14.6)') xl(mod(l,3) +1, 1:2), &
                     q(i,subtri(ideg,k,mod(l,3) + 1))
             endif
          enddo
          write(ignu,'(x)')
          write(ignu,'(x)')
       enddo
    enddo

  end subroutine Plot3D


  subroutine PlotVectors(ignu, nelem, npoin, max_deg, max_dof, &
              subtri, lambda, x, lnd, deg, qx, qy)
    integer, intent (in) :: ignu, nelem, npoin, max_deg, max_dof
    integer, dimension(0:max_deg, 1: max_deg**2, 1:3), intent(in) :: subtri
    real, dimension(0:max_deg, 1: max_dof, 1:3), intent(in) :: lambda
    real, dimension(1:npoin, 1:2), intent(in) :: x
    integer, dimension(1:nelem, 1:3), intent(in) :: lnd
    integer, dimension(1:nelem), intent(in) :: deg
    real, dimension(1:nelem, 1:max_dof), intent(in) :: qx, qy
    real, dimension(:,:), allocatable :: xl, R, R1
    real, dimension(:), allocatable :: v
    real :: ratio, alpha, ratio1
    integer :: i,k, l,  ideg, il(1:3)

    allocate(v(1:2),  xl(1:3, 1:2), R(1:2, 1:2), R1(1:2, 1:2) )

    alpha = 0.05
    R(1,1) = cos(alpha)
    R(2,1) = sin(alpha)
    R(1,2) =-sin(alpha)
    R(2,2) = cos(alpha)
    
    R1(1,1) = cos(alpha)
    R1(2,1) =-sin(alpha)
    R1(1,2) = sin(alpha)
    R1(2,2) = cos(alpha)
    
    do i=1,nelem
       ideg = max(1, deg(i) )

       do k=1, ideg**2
          il(1:3) = subtri(ideg, k, 1:3)
          
          do l=1,3
             xl(l,1) = dot_product(lambda(ideg, il(l), 1:3), x(lnd(i, 1:3),1) )
             xl(l,2) = dot_product(lambda(ideg, il(l), 1:3), x(lnd(i, 1:3),2) )
          enddo
         
          ratio = max(VectorNorm( xl(1,:) - xl(2,:)), VectorNorm( xl(1,:) - xl(3,:)), &
               VectorNorm( xl(2,:) - xl(3,:)) )
          ratio1 = ratio * 0.9

          do l=0,3
             if(deg(i) == 0) then
                v(1) = qx(i,1)
                v(2) = qy(i,1)
             else
                v(1) = qx(i,subtri(ideg,k,mod(l,3) + 1))
                v(2) = qy(i,subtri(ideg,k,mod(l,3) + 1))
             endif

             !if(xl(1,1) < 0.055) write(*,'(a6,6es12.4)') 'vvv:',v(1:2)
             
             write(ignu,'(3es14.6)') xl(mod(l,3) +1, 1:2)
             write(ignu,'(3es14.6)') xl(mod(l,3) +1, 1:2) + ratio * v(1:2) 
             write(ignu,'(3es14.6)') xl(mod(l,3) +1, 1:2) + ratio1 * matmul(R(1:2, 1:2), v(1:2) )
             write(ignu,'(3es14.6)') xl(mod(l,3) +1, 1:2) + ratio1 * matmul(R1(1:2, 1:2), v(1:2) )
             write(ignu,'(3es14.6)') xl(mod(l,3) +1, 1:2) + ratio * v(1:2) 
             write(ignu,'(x)')
             write(ignu,'(x)')

          enddo
          write(ignu,'(x)')
          write(ignu,'(x)')
       enddo
    enddo

    deallocate (v, xl, R, R1)
  end subroutine PlotVectors


  subroutine PlotWalls(ignu, nelem, npoin, nbelm, max_deg, max_dof, &
       lambda, subedge, x, lnd, lbn, deg, q, iw, iwall)
    integer, intent (in) :: ignu, nelem, npoin, max_deg, max_dof, nbelm, iw
    integer,dimension(0:max_deg, 1:3, 1: max_deg+1), intent(inout) :: subedge
    real, dimension(0:max_deg, 1: max_dof, 1:3), intent(in) :: lambda
    real, dimension(1:npoin, 1:2), intent(in) :: x
    integer, dimension(1:nelem, 1:3), intent(in) :: lnd
    integer, dimension(1:nbelm, 1:5), intent(in) :: lbn
    integer, dimension(1:nelem), intent(in) :: deg
    integer, dimension(1:iw), intent(in) :: iwall
    real, dimension(1:nelem, 1:max_dof), intent(in) :: q
    real :: xl(1:2)
    integer :: i, ib, k, l, j, ideg, il
    logical :: iplot 


    do ib=1,nbelm
       iplot = .false.
       
       do l=1, iw
          if(lbn(ib, 3) == iwall(l) ) iplot = .true.
       enddo

       if(iplot) then
          i = lbn(ib, 4)
          j = lbn(ib, 5)

          ideg = max(1, deg(i) )
          
          do k=1, ideg+1
             il = subedge(ideg, j, k)
             
             xl(1) = dot_product(lambda(ideg, il, 1:3), x(lnd(i, 1:3),1) )
             xl(2) = dot_product(lambda(ideg, il, 1:3), x(lnd(i, 1:3),2) )
             
             if(deg(i) == 0) then
                write(ignu,*) xl(1), q(i,1), xl(2),-1.4
             else
                write(ignu,*) xl(1), q(i,subedge(ideg,j, k)), xl(2),-1.4
             endif
          enddo

          write(ignu,'(x)')
          write(ignu,'(x)')
       end if
    enddo

  end subroutine PlotWalls

  subroutine PlotIsolines(ignu, nelem, npoin, nbelm, max_deg, max_dof, &
              subtri, lambda, x, lnd, lbn, deg, q, nisol, irange, qmin, qmax, &
              r1, r2, r3, r4)
    integer, intent (in) :: ignu, nelem, npoin, nbelm, max_deg, max_dof
    integer, dimension(0:max_deg, 1: max_deg**2, 1:3), intent(in) :: subtri
    real, dimension(0:max_deg, 1: max_dof, 1:3), intent(in) :: lambda
    real, dimension(1:npoin, 1:2), intent(in) :: x
    integer, dimension(1:nelem, 1:3), intent(in) :: lnd
    integer, dimension(1:nbelm, 1:5), intent(in) :: lbn
    integer, dimension(1:nelem), intent(in) :: deg
    real, dimension(1:nelem, 1:max_dof), intent(in) :: q
    integer, intent(in) :: nisol, irange
    real, intent(inout) :: qmin, qmax
    real, intent(inout) :: r1, r2, r3, r4  ! parameters of the peridocicity
    real, dimension(:), allocatable :: xshift

    real :: qval, ql(1:3)
    real :: xl(1:3, 1:2), xq(1:2)
    integer :: i,k, j, j1, l, iq, ideg, il(1:3), idof
    logical :: qfound
    integer :: ipoc, npoc

    allocate(xshift(1:2) )
    
    if(irange == 0) then
       qmin = 1E+20
       qmax = -1E+20
       do i=1,nelem
          idof = (deg(i) + 1)*(deg(i) + 2)/2
          qmax = max (qmax, maxval (q(i,1:idof) ) )
          qmin = min (qmin, minval (q(i,1:idof) ) )
       enddo
    endif

    ipoc = 0
    xshift(:) = 0.
    if(abs(r1) + abs(r2) > 1e-20) npoc = 2 


    10 continue
    do i=1,nelem
       ideg = max(1, deg(i) )

       do k=1, ideg**2
          il(1:3) = subtri(ideg, k, 1:3)
          
          do l=1,3
             xl(l,1) = dot_product(lambda(ideg, il(l), 1:3), x(lnd(i, 1:3),1) )
             xl(l,2) = dot_product(lambda(ideg, il(l), 1:3), x(lnd(i, 1:3),2) )
          enddo
         
          do l=1,3
             if(deg(i) == 0) then
                ql(l) = q(i,1)
             else
                ql(l) = q(i,subtri(ideg,k, l))
             endif
          enddo

          !do l=1,3
          !   write(*,'(a3,5es12.4)') 'w:',xl(l,1:2), ql(l), qmin, qmax
          !enddo

          do iq = 0, nisol
             qval = qmin + 1.*iq/nisol*(qmax - qmin)

             qfound = .false.

             do j=1,3
                j1 = mod(j,3) + 1

                if(abs(ql(j) - ql(j1) ) < 1E-6) then
                else


                   if(  (ql(j) <= qval .and. qval <= ql(j1)) .or. &
                        (ql(j) >= qval .and. qval >= ql(j1)) ) then
                      
                      !write(*,'(a3,3i5,3es12.4)') '?>?',iq,j,j1,ql(j),ql(j1),qval
                      

                      xq(1:2) = xl(j, 1:2) &
                           + (qval - ql(j))/(ql(j1) - ql(j) ) * (xl(j1,1:2) - xl(j,1:2) )

                      !write(ignu,'(3es16.8)') xq(1:2), qval
                      write(ignu,*) xq(1:2)+xshift(1:2), qval
                                               
                      
                      qfound = .true.
                   endif
                endif
             enddo
             if(qfound) then
                write(ignu,'(x)')
                write(ignu,'(x)')
             endif

          enddo
       enddo
    enddo

    ! border
    do i=1,nbelm
       if(npoc == 0 .or. (lbn(i,3) /= 5 .and. lbn(i,3) /= 6)  &
            .or. (lbn(i,3) == 6 .and. ipoc == 1) &
            .or. (lbn(i,3) == 5 .and. ipoc == 2) ) then
          write(ignu,'(3es14.6)') x(lbn(i,1), 1:2)+xshift(1:2), qmin
          write(ignu,'(3es14.6)') x(lbn(i,2), 1:2)+xshift(1:2), qmin
          write(ignu,'(x)')
          write(ignu,'(x)')
       endif
    enddo

    if(npoc > 0 .and. ipoc == 0 ) then
       xshift(1) =  r1
       xshift(2) =  r2
       ipoc = ipoc + 1
       goto 10
    elseif(npoc > 0 .and. ipoc == 1 ) then
       xshift(1) =  -r1
       xshift(2) =  -r2
       ipoc = ipoc + 1
       goto 10
    endif




    deallocate(xshift)

  end subroutine PlotIsolines


  subroutine PlotCUT(ignu, nelem, npoin, max_deg, max_dof, &
              subtri, lambda, x, lnd, deg, q, xA, xB)
    integer, intent (in) :: ignu, nelem, npoin, max_deg, max_dof
    integer, dimension(0:max_deg, 1: max_deg**2, 1:3), intent(in) :: subtri
    real, dimension(0:max_deg, 1: max_dof, 1:3), intent(in) :: lambda
    real, dimension(1:npoin, 1:2), intent(in) :: x
    real, dimension(1:2), intent(in) :: xA, xB
    integer, dimension(1:nelem, 1:3), intent(in) :: lnd
    integer, dimension(1:nelem), intent(in) :: deg
    real, dimension(1:nelem, 1:max_dof), intent(in) :: q
    real :: qval
    real :: xl(1:3, 1:2), ql(1:3)
    real :: x1(1:2), x2(1:2), xp(1:2)
    real :: wp(1:3, 1:3)
    integer :: i, j, j1, k, l,  ideg, il(1:3), ip
    logical :: inters, ifound

    do i=1,nelem

       ideg = max(1, deg(i) )

       do k=1, ideg**2
          il(1:3) = subtri(ideg, k, 1:3)
          
          do l=1,3
             xl(l,1) = dot_product(lambda(ideg, il(l), 1:3), x(lnd(i, 1:3),1) )
             xl(l,2) = dot_product(lambda(ideg, il(l), 1:3), x(lnd(i, 1:3),2) )
          enddo
         
          do l=1,3
             if(deg(i) == 0) then
                ql(l) = q(i,1)
             else
                ql(l) = q(i,subtri(ideg,k, l))
             endif
          enddo


          ifound = .false.
          ip = 0
          do j=1,3
             j1 = mod(j, 3) + 1

             x1(1:2) = xl(j, 1:2)
             x2(1:2) = xl(j1, 1:2)
          
             call IntersectOfLines(x1, x2, xA, xB, xp, inters)

             if(inters) then
                ip = ip + 1
                qval = ql(j) + Distance(xp, x1) / Distance(x2, x1) *(ql(j1) - ql(j) )

                wp(ip, 1:2) = xp(1:2)
                wp(ip, 3) = qval

                !write(ignu, *) xp(1),qval, xp(2)
                ifound = .true.
             endif
             
          enddo
          if(ifound) then
             if(ip == 1) then ! only touch of the elementwe do not draw the segment
             elseif(ip == 2) then
                write(ignu, *) wp(1,1), wp(1,3), wp(1,2)
                write(ignu, *) wp(2,1), wp(2,3), wp(2,2)

                write(ignu, *) 
                write(ignu, *) 

             elseif(ip == 3) then
                write(ignu, *) wp(1,1), wp(1,3), wp(1,2)
                write(ignu, *) wp(2,1), wp(2,3), wp(2,2)
                write(ignu, *) wp(3,1), wp(3,3), wp(3,2)

                write(ignu, *) 
                write(ignu, *) 
             endif

          endif

       enddo
    enddo

  end subroutine PlotCUT


  !> intersection of lines (a,b) and (c,d) is a node x
  subroutine IntersectOfLines(a, b, c, d, x, inters)
    real, dimension(1:2), intent(in):: a,b,c,d
    real, dimension(1:2), intent(out):: x
    logical, intent(out):: inters

    real ::st, s, t, rlen1, rlen2;

    st = (b(2)-a(2))*(d(1)-c(1))+(a(1)-b(1))*(d(2)-c(2));
    inters = .false.

    if(st /= 0.) then

       s = (b(1)*(c(2)-a(2)) + a(1)*(b(2)-c(2)) + c(1)*(a(2)-b(2)))/st
       t = (a(2)*(c(1)-d(1)) +c(2)*(d(1)-a(1)) + d(2)*(a(1)-c(1)))/st

       rlen1 = ((a(1)-b(1))*(a(1)-b(1)) + (a(2)-b(2))*(a(2)-b(2)))**0.5
       rlen2 = ((c(1)-d(1))*(c(1)-d(1)) + (c(2)-d(2))*(c(2)-d(2)))**0.5


       if( s >= 0. .and.  s <= 1. .and.  t >= 0. .and.  t <= 1. )then
          x(1) = a(1) + t*(b(1)-a(1))
          x(2) = a(2) + t*(b(2)-a(2))
          inters = .true.
       endif
    endif
  end subroutine IntersectOfLines
    


  subroutine  PlotHPColorMesh(npoin, nelem, x, lnd, deg)
    integer, intent(in) :: npoin, nelem
    real, dimension (1:npoin, 1:2), intent(in):: x
    integer, dimension(1:nelem, 1:3),intent(in) :: lnd
    integer, dimension (1:nelem, 1:2), intent(in) :: deg
    character(len=50) :: gnufile
    real, dimension(1:3,1:3) :: xi
    real :: q
    integer :: i,j, j1, j2, p, pmax, pmin, ifile

    ifile = 121
    
    q = 0.99

    pmin = 0
    pmax = 10
    do p=pmin, pmax
       call SetHPmeshFilename(p, gnufile)

       open(ifile+p, file=gnufile, status='replace')
       write(ifile+p,*) x(1,1:2)
    enddo

    do i=1, nelem
       do j=1,3
          j1 = mod(j, 3) + 1
          j2 = mod(j1, 3) + 1
          xi(j, 1:2) = q * x(lnd(i,j), 1:2)  &
               + (1-q)/2 * x(lnd(i,j1), 1:2)  &
               + (1-q)/2 * x(lnd(i,j2), 1:2)
       enddo

       write(ifile + deg(i,1),*)  xi(1, 1:2)
       write(ifile + deg(i,1),*)  xi(2, 1:2)
       write(ifile + deg(i,1),*)  xi(3, 1:2)
       write(ifile + deg(i,1),*)  xi(1, 1:2)
       write(ifile + deg(i,1),'(x)')
    enddo

       
    do p=pmin, pmax
       close(ifile+p)
    enddo

  end subroutine PlotHPColorMesh


  function VectorNorm(x)
    real :: VectorNorm
    real, dimension(:), intent(in) :: x
    
    VectorNorm = sqrt(dot_product(x,x))

  end function VectorNorm

  function Distance(p,q)
    real :: Distance
    real, dimension(:), intent(in) :: p,q
    if(size(p) /= size(q) ) then
       print*,'Different dimension of p and q in Distance(p,q)'
       Distance = 0.
    else
       Distance = sqrt(dot_product(q-p,q-p))
    endif
  end function Distance


end module plot_geom
