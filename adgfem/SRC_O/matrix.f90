!> parameters \f$ \eta = 1/\tau\f$ and \f$  \iota \f$ used in double problem
module matrix_oper_int
implicit none

  public:: eta, iota

  real:: eta = 0., iota(2) = (/0., 0./)

end module

!> matrix operations for sparse (block) matrices
module matrix_oper
  use lapack_oper
  use geometry
  use element_mod
  use define_state
  use main_data
  use data_mod

  implicit none

!!  public:: sMVprod            ! product of a sparse matrix and a vector OLD
  public:: MatrixNormFrobenius
  public:: null_precond              ! No preconditioner
  public:: BlockJacobi        ! block Jacobi preconditioner
  public:: bMassInvVprod      ! block inverse mass matrix - vector product
  public:: bMVprod            ! block matrix(M+tC)-vector product
  public:: bMVprodST	      ! block matrix(M+tC)-vector product in STDGM
  public:: bMVprodST_Dual	      ! block matrix(M+tC)-vector product in STDGM DUAL PROBLEM
  public:: MGbMVprod	      ! !! be sure grid%elem%MG*** is initialized !!
  public:: bMVprod2           ! double block matrix-vector product
  public:: bMVprodA           ! block matrix(A)-vector product
  public:: bMVprodOffC        ! block matrix-vector product only off diagonal terms
  public:: bMVmassprod        ! block mass matrix - vector product
  public:: EvalSSresid        ! eval steady-state residuum := C(w)w - q(w)
  !public:: EvalSSresidDirect  ! eval steady-state residuum := f(w) (= q(w) - C(w)w )
  public:: EvalSSresidExplicit! eval steady-state residuum := f(w) (= q(w) - C(w)w )
  public:: EvalWeightTDresid  ! eval weighted Time Dependent  residuum

  public:: VectorPrecondNorm  ! "preconditioned" norm of the vector
  public:: VectorScaleNorm    ! "scaled" norm of the vector
  public:: bMVnull            ! no preconditioner
  public:: bMVdiagprod        ! diag-block matrix-vector product
  public:: bMVdiagprod2       ! double diag-block matrix-vector product
  public:: bMVdiagprodST 		! diag-block matrix-vector product for STDGM
  public:: bMVdiagprodST_Dual ! diag-block matrix-vector product for STDGM DUAL PROBLEM
  public:: bMViLUprod         ! performs ILU preconditioning
  public:: bMViLUprod2        ! performs double ILU preconditioning
  public:: bMViLUprodST 		! performs ILU preconditioning for STDGM
  public:: bMViLUprodST_Dual  ! performs ILU preconditioning for STDGM DUAL PROBLEM



  public:: CopyBlocksSTtoBlockPlus ! copy bigger diagonal block to blockPlus, used for Ritz reconstruction
  public:: TaylorSolution     ! solution of (M+\tau C) w = b with Taylor serie
  public:: WriteMatrixA       ! write the matrix (M+\tau C)
  public:: WriteMatrixA_ST       ! write the matrix (1/tauM+C) for STDGM
  public:: WriteMatrixST_Blocks   ! write the nonzeroblocks of matrix (1/tau M+ C) for STDGM
  public:: WriteMatrixST       ! write the RefTimeMatrix for STDGM
  public:: test_bMVprodST     !testing bMVprodST
  public:: Write_rhsST        !write RHS for STDGM

  public:: WriteMatrixAMatlab ! write the matrix (M+\tau C) into a file (vector row indices, vector column indices, nonzero values), input for Matlab
  public:: WriteMatrixLU      ! write the matrix ILU

  public:: prodFEM            !  matrix-vector product for conforming FEM
  public:: diagFEM            !  diagonal preconditioner for conforming FEM

contains

  ! evaluate the Froenius norm of the block matrix elem%block(:)
  function MatrixNormFrobenius()
    real :: MatrixNormFrobenius
    class(element), pointer:: elem
    integer :: i,k,in,j, ndof, ndof1
    real :: val

    val = 0
    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof = elem%dof*ndim

       do j=1,ndof
          val = val + dot_product(elem%block(0)%Mb(j,1:ndof), elem%block(0)%Mb(j,1:ndof))
       enddo

       do k=1,elem%flen
          in = elem%face(neigh,k)
          if(in >0) then
             ndof1 = grid%elem(in)%dof * ndim

             do j=1,ndof
                val = val + dot_product(elem%block(k)%Mb(j,1:ndof1), &
                     elem%block(k)%Mb(j,1:ndof1) )
             enddo
          endif
       enddo

    enddo

    MatrixNormFrobenius = val**0.5

  end function MatrixNormFrobenius

  subroutine WriteMatrixA(eta)
    real, intent(in) :: eta
    class(mesh), pointer :: grid_print
    class(element), pointer:: elem,elem1 ! one element
    real, dimension(:),allocatable:: accum
    integer :: i,k,in,j,j1, k1, is, is1, dof, ndof, ndof1


    grid_print => grid
    if(state%local_problem) grid_print => gridL

    allocate(accum(1:state%nsize))

    print*,'--------------  matrix 1/tau*M + C  -----------------', state%nsize
    do i=1,grid_print%nelem
       elem => grid_print%elem(i)
       is = elem%ncv
       dof = elem%dof
       ndof = dof*ndim

       !do j=1,ndof
       do k1=1,ndim
          do j1= 1, dof

             j = (k1 -1)*dof + j1

             accum(:) = 0.

             if(eta /= 0.) then
                !print*,'^^^',j,dof, is, is + (j-1)/dof *dof, is + (j-1)/dof *dof + dof-1,'||', &
                !     size(elem%Mass%Mb, 1), size(elem%Mass%Mb, 2)

                accum(is + (j-1)/dof *dof: is + (j-1)/dof *dof + dof-1) &
                     = eta * elem%Mass%Mb(j1,1:dof)
             endif

             accum(is:is+ndof-1) =  accum(is:is+ndof-1) + elem%block(0)%Mb(j,1:ndof)


             do k=1,elem%flen
                in = elem%face(neigh,k)
                if(in >0) then
                   elem1 => grid_print%elem(in)
                   ndof1 = elem1%dof * ndim
                   is1 = elem1%ncv
                   accum(is1:is1+ndof1-1) = elem%block(k)%Mb(j,1:ndof1)
                endif
             enddo

             !  write(*,'(i5,a2,100es11.3)')is+j-1,': ',accum(:)
             write(*,'(i5,a2,500es9.1)')is+j-1,': ',accum(:) !
          enddo ! j1
       enddo ! k1
    enddo ! i
    print*,'-----------end of   matrix 1/tau*M + C  -----------------'

    deallocate(accum)
  end subroutine WriteMatrixA

  subroutine test_bMVprodST()
    real, dimension(:,:), allocatable :: b
    real, dimension(:), allocatable ::   x
    class(element), pointer :: elem, elem1
    integer :: i,j,k, m, ndof,ndof1,is,is1, face

    open (58, file="bMVprodST-ADGo", action ="write", status="replace")

    allocate( b( 1:state%nsize,1:state%nsize))
    allocate( x( 1:state%nsize))
    b(:,:) = 0

    do i = 1, state%nsize
       x(:) = 0
       x(i) = 1
       call bMVprodST(b(:,i),x,state%nsize)
    enddo !i

    ! 	do i = 1,state%nsize
    ! 		write(58,'(100es10.2)'), b(i,:)
    ! 	enddo

    do i = 1, grid%nelem
       elem => grid%elem(i)
       ndof = elem%Tdof * elem%dof * ndim
       is = elem%ncv
       write(58,*), 'Element ', i
       write(58,*), 'block(0):', is, is + ndof  - 1
       do j = 1, ndof
          write(58,'(200es12.4)'), b(is + j - 1 , is : is + ndof  - 1)
       enddo !j

       do j = 1,elem%flen
          face = elem%face(neigh,j)
          if (face > 0) then
             write(58,*), 'block(',j,'):'
             elem1 => grid%elem(face)
             ndof1 = elem1%dof * ndim * elem1%Tdof
             is1 = elem1%ncv

             do m = 1, ndof
            	write(58, '(200es12.4)' ) b(is + m - 1, is1 : is1 + ndof1 - 1)
             enddo !m

          endif
       enddo !j

       write(58,*), '--------------------------------'
    enddo !i

    close(58)
     print*, ''
     print*, 'End of test_bMVprodST - it shouldnt be used for large systems!!!!!!!'
    deallocate(x)
    deallocate(b)

  end subroutine test_bMVprodST

  subroutine Write_rhsST()
    real, dimension(:,:), allocatable :: b
    real, dimension(:), allocatable ::   x
    class(element), pointer :: elem, elem1
    integer :: ie, kvec,ivec,dof, k, l, ndof,Tdof

    !open (58, file="rhsST", action ="write", status="replace")
    open (58, file="bMVprodST", action ="write", status='UNKNOWN', position = 'append')

!  write(58, *), 'RHS ST'
!

!
!   do ie=1,grid%nelem
!      elem => grid%elem(ie)
!      dof = elem%dof_plus
!      Tdof = elem%Tdof
!      ndof = dof * Tdof * ndim

!      !b(ivec+1:ivec+ndof) = 0 !elem%rhsST()
!    !  kvec = ivec
!      write(58, *), 'elem=', ie, 'tdof=',Tdof, 'ndim=',ndim,'dof=', elem%dof, 'dof_plus=' , elem%dof_plus
!      do l = 1, Tdof
!         do k = 1, ndim
!            write(58, '(a5, i5, a5, i5, 100es12.4)') 'Tdof=',l,  'ndim=', k,  elem%rhsST(k,1:dof, l)
!          !  kvec = kvec + dof
!         enddo !k
!      enddo !l
!      write(58,*), '--------------------------------'
!     ! ivec = ivec + ndof
!   end do

  write(58,*), 'Newton b--------------------------------'
   ivec = 0

   do ie=1,grid%nelem
      elem => grid%elem(ie)
      dof = elem%dof
      ndof = dof * elem%Tdof * ndim


      kvec = ivec

      do l = 1, elem%Tdof
         do k = 1, ndim
             write(58, '(4i5, 200es12.4)') l , k ,kvec+1, kvec+dof, state%nlSolver%b(kvec+1:kvec + dof)
            kvec = kvec + dof
         enddo !k
      enddo !l
      write(58,*), '--------------------------------'
      ivec = ivec + ndof
   end do

   write(58,*), 'wSTfinAD--------------------------------'
   do ie= 1, grid%nelem
      elem => grid%elem(ie)

      do k = 1, ndim
         write(58, *) k , elem%wSTfinAD(k, 1:elem%dof)
      enddo
      write(58,*), '--------------------------------'
   enddo
!   write(58,*) 'wStfin_____________________'
!   do ie= 1, grid%nelem
!      elem => grid%elem(ie)
!
!      do k = 1, ndim
!         write(58, *) k , elem%wSTfin(k, 1:elem%Qdof)
!      enddo
!      write(58,*), '--------------------------------'
!   enddo

  close(58)
  end subroutine Write_rhsST



  subroutine WriteMatrixA_ST(eta)
    real, intent(in) :: eta
    class(element), pointer:: elem,elem1 ! one element
    real, dimension(:),allocatable:: accum
    integer :: i,j,k,l,m,n,p,mm,nn,kk,jj,pp,r
    integer :: is, is1, face, dof, ndof, ndof1,Tdof

    allocate(accum(1:state%nsize))

    associate( time => state%time )
    select type( time )
    type is ( TimeTDG_t )

       print*,'--------------  matrix 1/tau*M + C  -----------------', state%nsize
       do i=1,grid%nelem
          elem => grid%elem(i)

          is = elem%ncv

          dof = elem%dof
          Tdof = elem%Tdof
          ndof = dof * ndim
      !	print*, 'Velikost block0:', size(elem%block(0)%Mb(1,:))
      !	print*, 'Velikost blockST:', size(elem%blockST(0)%Mb(1,:))
         do m = 1, Tdof
            mm = (m-1)*ndof
            do j = 1, ndim
               do k = 1, dof
                  accum(:) = 0
                  if (eta > 0) then
                  do n = 1,Tdof
                     nn = (n-1) * ndof
                     do p = 1, ndim
                        pp = (p-1)*ndim
                        accum(is + nn + pp : is + nn + pp + dof - 1) = eta *  & !eta *
                             time%refTimeMatrix%Mb(m,n) ! * elem%Mass%Mb(k,1:dof)
                     enddo !p
                  enddo !n
                  endif

                  accum(is : is + Tdof * ndof) =  accum(is : is + Tdof * ndof) &
                      + elem%blockST(0)%Mb(mm + (j-1)*ndim + k, 1:Tdof*ndof )
                  ! write(*,'(i5,a2,100es11.4)')is+mm+(j-1)*ndim + k - 1,': ',accum(is : is + Tdof * ndof - 1)

                  do r=1,elem%flen
                     face = elem%face(neigh,r)
                     if(face > 0) then
                        elem1 => grid%elem(face)
                        ndof1 = elem1%dof * ndim * elem1%Tdof
                        is1 = elem1%ncv
                        accum(is1: is1 + ndof1 -1) = elem%blockST(r)%Mb(mm + (j-1)*ndim + k, 1:ndof1 )
                     endif
                  enddo !r

                  write(*,'(i5,a2,100es8.0)')is+mm+(j-1)*ndim + k - 1,': ',accum(:)

               enddo !k
            enddo !j
         enddo !m

      enddo !i
      class default
         stop 'Reference time matrix is allocated only for STDG'
      end select
      end associate

    deallocate(accum)
  end subroutine WriteMatrixA_ST

  !writes Matrix eta*MassST + blockST into file ../stdgm/Matrix
   subroutine WriteMatrixST_Blocks(eta)
    real, intent(in) :: eta
    class(element), pointer:: elem,elem1 ! one element
    real, dimension(:),allocatable:: accum
    integer :: i,j,k,l,m,n,p,mm,nn,kk,jj,pp,r
    integer :: is, is1, face, dof, ndof, ndof1,Tdof

    associate( time => state%time )
    select type( time )
    type is ( TimeTDG_t )

       !open (59, file="../Tests/Matrix", action ="write", status="replace")
       open (59, file="MatrixST-ADGo", action ="write", status="replace")
      !	if ( ierror /= 0 ) then
      !		print*, "Failed to open test.dat!"
      !		stop
      !	end if


       write(59,*) '--------------  matrix 1/tau*M + C  -----------------', state%nsize
       !only 1-st elem now

       do i=1, state%time%max_Tdof
          !write(59,'(i5,a2,100f8.2)') i,': ', state%time%refTimeMatrix%Mb(i,1:state%time%max_Tdof)
          write(59,'(i5,a2,100es12.4)') i,': ', time%refTimeMatrix%Mb(i,1:time%max_Tdof)
       enddo
       write(59,*) , 'eta=' , eta
       do i=1, grid%nelem
       !do i= 1721, 1726
         write(59,*) '--------------------'
         write(59,*) 'element(', i, ')'
          elem => grid%elem(i)

          is = elem%ncv

          dof = elem%dof
          Tdof = elem%Tdof
          ndof = dof * ndim
!          write(59,*) 'Mass:'
!          do m =1,dof
!             write(59,'(100es12.4)') elem%Mass%Mb(m,:)
!          enddo !m

          allocate(accum(1:ndof*Tdof))
          !	print*, 'Velikost block0:', size(elem%block(0)%Mb(1,:))
          !	print*, 'Velikost blockST:', size(elem%blockST(0)%Mb(1,:))
          write(59,*) 'massST:'
          do m = 1, Tdof
             mm = (m-1)*ndof
             do j = 1, ndim
                do k = 1, dof
                   accum(:) = 0
!                   if (eta > 0) then
                      do n = 1,Tdof
                         nn = (n-1) * ndof
                        ! do p = 1, ndim
                            ! <> 0 only when the same component meet
                            p = j
                            pp = (p-1)*dof
                            accum(nn + pp + 1 : nn + pp + dof) = eta   & !eta *
                                 * elem%Mass%Mb(k,1:dof) * time%refTimeMatrix%Mb(m,n)
                        ! enddo !p
                      enddo !n
!                   else
!
!                   print*, 'WriteMatrixST_Blocks not implemented for eta = 0'
!                   stop
!                   endif

                   !		accum(1: Tdof * ndof) =  accum(1 : Tdof * ndof) &
                   !				 + elem%blockST(0)%Mb(mm + (j-1)*ndim + k, 1:Tdof*ndof )

                   write(59,'(i5,a2,100f8.2)') mm + (j-1)*dof + k ,': ',accum(:)
                   !write(59,'(i5,a2,200es12.4)') mm + (j-1)*dof + k ,': ',accum(:)


                   !write(*,'(i5,a2,100es8.0)')is+mm+(j-1)*ndim + k - 1,': ',accum(:)

                enddo !k
             enddo !j
          enddo !m

          deallocate(accum)

          write(59,*) 'block(0):', size(elem%blockST(0)%Mb(:,1)), 'x' , size(elem%blockST(0)%Mb(1,:))
          do m = 1, Tdof*dof*ndim
             write(59,'(i5,a2,100f8.2)') m,': ', elem%blockST(0)%Mb(m,:)
             !write(59,'(i5,a2,200es12.4)') m,': ', elem%blockST(0)%Mb(m,:)
          enddo !m

   !offdiagonal block commented
         do r=1,elem%flen
             face = elem%face(neigh,r)
                		if(face > 0) then
                			write(59,*), 'block(',r,'):'
                 		   elem1 => grid%elem(face)
                   		ndof1 = elem1%dof * ndim * elem1%Tdof
                   		is1 = elem1%ncv
                   		allocate(accum(1:ndof1))
                   		do m = 1, ndof*Tdof
                   			accum(1: ndof1) = elem%blockST(r)%Mb( m, 1:ndof1 )
                   			write(59,'(i5,a2,100f8.2)'), m - 1,': ',accum(:)
!                   			write(59,'(i5,a2,100es10.2)'), m - 1,': ',accum(:)
                   		enddo !m
                   		deallocate(accum)
                		endif
         enddo !r

         enddo !i


      close(59)

   class default
      stop  'WriteMatrixST_Blocks for STDGM only'
   end select
   end associate

  end subroutine WriteMatrixST_Blocks

  subroutine WriteMatrixST()
    real, dimension(:),allocatable:: accum
    integer :: Tdof,i,j

    associate( time => state%time )
    select type( time )
    type is ( TimeTDG_t )

      Tdof = state%time%max_Tdof

      !allocate(accum(1:state%nsize))

      !print*, 'nelem=', grid%nelem

      print*,'--------------  Time part of Mass matrix -----------------', Tdof


         print*, 'Tdof = ', Tdof
         do i = 1,Tdof
            do j = 1,Tdof
               write(*,'(i5,a2,i5,a2,100es11.4)') i,',',j,': ', time%refTimeMatrix%Mb(i,j)
      !         write(*,'(a2, es5.5)') '  ',
      !
            enddo !j
         enddo !i =1,Tdof

    class default
      stop  'WriteMatrixST_Blocks for STDGM only'
    end select
    end associate


  end subroutine WriteMatrixST


    subroutine WriteMatrixAMatlab(eta)
    real, intent(in) :: eta
    class(element), pointer:: elem,elem1 ! one element
    real, dimension(:),allocatable:: accum
    integer :: i,k,in,j,is, is1, dof, ndof, ndof1
    character(len=7) :: MatrixA
    integer :: l

    allocate(accum(1:state%nsize))

    MatrixA = 'MatrixA'
    open(21, file=MatrixA, status='UNKNOWN', position = 'append')

    print*,'--------------  matrix 1/tau*M + C  -----------------', state%nsize
    do i=1,grid%nelem
       elem => grid%elem(i)
       is = elem%ncv
       dof = elem%dof
       ndof = dof*ndim

       do j=1,ndof
          accum(:) = 0.

          if(eta /= 0.) then
             print*,'^^^',j,dof, is + (j-1)/dof *dof, is + (j-1)/dof *dof + dof-1
             accum(is + (j-1)/dof *dof: is + (j-1)/dof *dof + dof-1) &
                  = eta * elem%Mass%Mb(j,1:dof)
          endif

          accum(is:is+ndof-1) =  accum(is:is+ndof-1) + elem%block(0)%Mb(j,1:ndof)


          do k=1,elem%flen
             in = elem%face(neigh,k)
             if(in >0) then
                elem1 => grid%elem(in)
                ndof1 = elem1%dof * ndim
                is1 = elem1%ncv
                accum(is1:is1+ndof1-1) = elem%block(k)%Mb(j,1:ndof1)
             endif
          enddo

          !write(*,'(i5,a2,100es11.3)')is+j-1,': ',accum(:)
          do l=1,state%nsize
            if ( accum(l)/= 0.  ) then
               write(21,'(i5,i5,5000es14.6)') is+j-1, l, accum(l)
            endif
          enddo !l

       enddo
    enddo
    print*,'-----------end of   matrix 1/tau*M + C  -----------------'

    close(21)

    deallocate(accum)
  end subroutine WriteMatrixAMatlab


  subroutine WriteMatrixLU()
    class(element), pointer:: elem,elem1 ! one element
    real, dimension(:),allocatable:: accum
    integer :: i,k,in,j,is, is1, dof, ndof, ndof1,Tdof,Tdof1

    allocate(accum(1:state%nsize))

    print*,'--------------  matrix ILU  -----------------'
    do i=1,grid%nelem
       elem => grid%elem(i)
       is = elem%ncv
       dof = elem%dof

       if (state%time%disc_time == 'STDG') then
       	Tdof = elem%Tdof
       else
       	Tdof = 1
       endif

       ndof = dof*ndim * Tdof

       do j=1,ndof
          accum(:) = 0.

          accum(is:is+ndof-1) =  accum(is:is+ndof-1) + elem%ILU(0)%Mb(j,1:ndof)


          do k=1,elem%flen
             in = elem%face(neigh,k)
             if(in >0) then
                elem1 => grid%elem(in)

                if (state%time%disc_time == 'STDG') then
                	Tdof1 = elem1%Tdof
                else
                	Tdof1 = 1
                endif

                ndof1 = elem1%dof * ndim * Tdof1
                is1 = elem1%ncv
                accum(is1:is1+ndof1-1) = elem%ILU(k)%Mb(j,1:ndof1)
             endif
          enddo

          write(*,'(i5,a2,100es12.4)')is+j-1,': ',accum(:)
       enddo
    enddo
    print*,'-----------end of   matrix ILU  -----------------'
    deallocate(accum)

  end subroutine WriteMatrixLU



  !> evaluation of the inverse (\f$ M^{-1} \f$ ) to Mblock Matrix \f$ M \f$,
  !> not further used
  subroutine  MblockLU(M, M1)
    type(Mblock), intent(in) :: M
    type(Mblock), intent(inout) :: M1
    real, dimension(:, :), allocatable :: L, U
    real, dimension(:), allocatable :: f, b, a
    integer :: dof, j, j1, k, k1

    if(size(M%Mb,1) /= size(M%Mb,2) .or. any(shape(M%Mb) /= shape(M1%Mb))) then
       print*,'Bad dimension in MblockLU in matrix.f90'
       print*,':',shape(M%Mb),shape(M%Mb)
       stop
    endif
    dof = size(M%Mb,1)

    ! TODO optimize

    allocate(L(1:dof, 1:dof), U(1:dof, 1:dof))
    allocate(a(1:dof), f(1:dof), b(1:dof))

    ! unit matrix
    L(1:dof, 1:dof) = 0.
    U(1:dof, 1:dof) = 0.
    do j=1, dof
       L(j,j) = 1.
    enddo

    ! let us compute the inverse elem%Mass by LU decomposition
    do k=1,dof
       do j=k,dof
          U(k,j) = M%Mb(k,j) - sum(L(k,1:k-1)*U(1:k-1,j) )
       enddo
       if(k /= dof) then
          do k1=k+1,dof
             L(k1,k) = (M%Mb(k1,k) - sum(L(k1,1:k-1)*U(1:k-1,k) )) /U(k,k)
          enddo
       endif
    enddo


    ! backward solution: LU a = e_k,  e_k is a canonical basis
    do k=1,dof
       f(1:dof) = 0.
       f(k) = 1.

       b(1) = f(1)/L(1,1)
       do j=2, dof
          b(j) = (f(j) - sum(b(1:j-1)*L(j,1:j-1) ) ) / L(j,j)
       enddo

       a(dof) = b(dof)/U(dof,dof)

       do j1 = 1, dof-1
          j = dof -j1
          a(j) = (b(j) - sum(a(j+1:dof) * U(j,j+1:dof) ) )/U(j,j)
       enddo

       M1%Mb(1:dof, k) = a(1:dof)
    enddo

    deallocate(L, U, f, b, a)

  end subroutine MblockLU

  !> block matrix - vector product: \f$ b = M^{-1}x \f$,  \f$ M^{-1}  \f$ is the
  !> inverse mass matrix in block form
  subroutine bMassInvVprod(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    class(element), pointer:: elem
    integer :: i, k,dof, ielem, is, ie

    do i=1,grid%nelem
       elem => grid%elem(i)
       dof = elem%dof

       ielem = elem%ncv
       do k=1,ndim
          is = ielem + (k-1)*dof
          ie = is + dof -1
          b(is:ie) =  matmul(elem%MassInv%Mb(1:dof,1:dof), x(is:ie) )
       enddo
    enddo

  end subroutine bMassInvVprod

  !> block matrix - vector product: \f$ b = (\eta M+ C_k) x \f$,
  !> \f$ M  \f$ is the block mass matrix grid.elem(*).Mass,
  !> \f$ C  \f$ is the block flux matrix grid.elem(*).block(*),
  !> \f$ \eta = 1/\tau_k\f$ is external
  subroutine bMVprod(b,x,nsize)
    use matrix_oper_int
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    class(element), pointer:: elem,elem1 ! one element
    integer :: i,j,k, ndof, ndof1, is, is1 !, l_accum
    real, dimension(:),allocatable :: accum
    !real, dimension(:),allocatable, save:: accum

    ! allocate accum once to accomodate for the largest dof.
    allocate(accum(maxval(grid%elem%dof) * ndim ) )
    !l_accum= maxval(grid%elem%dof) * ndim
    !if(size(accum) <= l_accum) then
    !   deallocate(accum)
    !   allocate(accum(1:l_accum))
    !endif

    !print*,'####!!!! eta=',eta
    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof1 = elem%dof
       ndof = elem%dof * ndim
       is = elem%ncv

       if (eta /= 0.) then
         do k = 0,ndof-1,ndof1
           accum(k+1:k+ndof1) = eta * matmul(elem%Mass%Mb(1:ndof1, 1:ndof1), &
                x(is+k: is+k+ndof1-1))
         enddo
         ! diagonal block
         accum(1:ndof) = accum(1:ndof) &
              + matmul(elem%block(0)%Mb(1:ndof, 1:ndof), x(is: is+ndof-1) )
       else
         ! diagonal block
         accum(1:ndof) = matmul(elem%block(0)%Mb(1:ndof, 1:ndof), x(is: is+ndof-1) )
       endif


       !! off-diagonal blocks
       do j=1,elem%flen
          k = grid%elem(i)%face(neigh,j)

          if(k > 0) then
             elem1 => grid%elem(k)
             ndof1 = elem1%dof * ndim
             is1 = elem1%ncv

             accum(1:ndof) = accum(1:ndof) &
                  + matmul(elem%block(j)%Mb(1:ndof, 1:ndof1), x(is1: is1+ndof1-1) )
          endif
       enddo
       b(is: is+ndof-1) = accum(1:ndof)

       ! if(i==1) then
       !    do k=1,ndof1
       !       write(*,'(a6,i5,300es12.4)') ' Mass:',k, elem%Mass%Mb(k, 1:ndof1)
       !    enddo
       ! write(*,'(a6,300es12.4)') ' x:',x(1:8)
       ! write(*,'(a6,300es12.4)') 'b=Ax:',b(1:8)

       ! endif


    enddo

    deallocate(accum)

  end subroutine bMVprod


  !> block matrix - vector product: \f$ b = (\eta M+ C_k) x \f$,
  !> \f$ M  \f$ is the block mass matrix grid.elem(*).Mass multiplied by state%time%refTimeMatrix,
  !> \f$ C  \f$ is the block flux matrix grid.elem(*).blockST(*),
  !> \f$ \eta = 1/\tau_k\f$ is external
  subroutine bMVprodST(b,x,nsize)
    use matrix_oper_int
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    class(element), pointer:: elem,elem1 ! one element
    integer :: i, j, k, l, m, n, mm, nn
    integer :: dof, dof1, ndof, ndof1, Tdof, is, is1
    real, dimension(:),allocatable :: accum

    associate( time => state%time )
    select type( time )
    type is ( TimeTDG_t )


       ! allocate accum once to accomodate for the largest dof.

       !allocate(accum(maxval(grid%elem%dof) * ndim ) )
       allocate( accum( state%space%max_dof * state%time%max_Tdof * ndim) )


       !print*,'####bMVprodST!!!! eta=',eta
       do i=1,grid%nelem
       ! 	print*, 'element', i
          accum(:) = 0.0
          elem => grid%elem(i)
          ndof1 = elem%dof
          Tdof = elem%Tdof
          ndof = elem%dof * ndim
          dof = ndof * Tdof

          is = elem%ncv

          if (eta /= 0.) then
             do m =1,Tdof
                mm = (m-1) * ndof
                do n = 1,Tdof
                   nn = (n-1) * ndof
                   do k = 0,ndof-1,ndof1
                      !write(*,'(2(a7,i5),a6,i5,a1,i5,3(a6,i5))') &
                      !     'elem',i,'ncv:', is, 'accum',mm + k + 1,':', mm+ k + ndof1, 'k =', k , 'm=', m, 'n=', n
                      !	print*,   (is+nn+k )
                      accum(mm + k + 1 : mm+ k + ndof1) = accum(mm + k + 1 : mm+ k + ndof1) + &
                           eta * time%refTimeMatrix%Mb(m,n) * &
                           matmul(elem%Mass%Mb(1:ndof1, 1:ndof1) , x(is+nn+k : is+nn+k+ndof1-1))

                   enddo !k

                enddo !n
             enddo !m

             ! diagonal block
             accum(1:dof) = accum(1:dof) &
                  + matmul(elem%blockST(0)%Mb(1:dof, 1:dof), x(is: is+dof-1) )

             !stop

          else
             ! diagonal block
             accum(1:dof) = matmul(elem%blockST(0)%Mb(1:dof, 1:dof), x(is: is+dof-1) )
          endif


         !! off-diagonal blocks
          do j=1,elem%flen
             k = grid%elem(i)%face(neigh,j)

             if(k > 0) then
                elem1 => grid%elem(k)
                dof1 = elem1%dof * ndim * elem1%Tdof
                is1 = elem1%ncv

                accum(1:dof) = accum(1:dof) &
                     + matmul(elem%blockST(j)%Mb(1:dof, 1:dof1), x(is1: is1+dof1-1) )
             endif
          enddo !j

          b(is: is+dof-1) = accum(1:dof)

          !if(dot_product(accum(1:dof), accum(1:dof)) > 1E-15 .and. state%time%recompute_back >= 2) &
          !     write(*,'(a8,2i5, 300es12.4)') 'Ax:',i, dof, accum(1:dof)
       enddo !i


       deallocate(accum)
    class default
      stop 'bMVprodST only for STDG method'

    end select
    end associate

    !print*,'#############E#E##E#E#', state%linSolver%iter
    !if(state%time%recompute_back == 2 .and. state%linSolver%iter >= 1 )stop "3d3ed388"

  end subroutine bMVprodST


  !> transposed block matrix - vector product: \f$ b = (\eta M+ C_k)^T x \f$,
  !> \f$ M  \f$ is the block mass matrix grid.elem(*).Mass multiplied by state%time%refTimeMatrix,
  !> \f$ C  \f$ is the block flux matrix grid.elem(*).blockST(*),
  !> \f$ \eta = 1/\tau_k\f$ is external
  subroutine bMVprodST_Dual(b,x,nsize)
    use matrix_oper_int
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
!    real, intent(in) :: eta ! 1 / \tau,  non/stationary problem
    class(element), pointer:: elem,elem1 ! one element
    integer :: i, j, k, l, m, n, mm, nn
    integer :: dof, dof1, ndof, ndof1, Tdof, is, is1, big_dof
    integer :: loc_neigh
    real, dimension(:),allocatable :: accum

    associate( time => state%time )
    select type( time )
    type is ( TimeTDG_t )

!       print*, 'ETA (bMVprodST_Dual)= ', eta

!      print*, 'nsize sa: ' , nsize, size(x(:))


       ! allocate accum once to accomodate for the largest dof.

       !allocate(accum(maxval(grid%elem%dof) * ndim ) )
       allocate( accum( state%space%max_dof * state%time%max_Tdof * ndim) )

       big_dof = 0
       !print*,'####!!!! eta=',eta
        ! stationary
       if (eta == 0.0) then
          do i=1,grid%nelem
            accum(:) = 0.0
            elem => grid%elem(i)
            ndof1 = elem%dof
            Tdof = elem%Tdof
            ndof = elem%dof * ndim
            dof = ndof * Tdof
            big_dof = big_dof + dof



            is = elem%ncv
            ! diagonal block - TRANSPOSED
            accum(1:dof) = matmul(  x(is: is+dof-1), elem%blockST(0)%Mb(1:dof, 1:dof) )

           !! off-diagonal blocks - we have to find them in elem(?)%blockST(??)
            do j=1,elem%flen
               k = elem%face(neigh,j)

               if(k > 0) then
                  elem1 => grid%elem(k) ! for Dual version purposes
                  dof1 = elem1%dof * ndim * elem1%Tdof ! same for dual and primal
                  is1 = elem1%ncv ! same for dual and primal
                  loc_neigh = elem%face(nei_i, j) ! local index of elem as an neighbor of elem1

!                  print*, 'DOF:' , dof, dof1
!                  print*, size( x(is1: is1+dof1-1) ) , size( elem1%blockST(loc_neigh)%Mb(1:dof1, 1) )
!
!                  print*, 'x: ' , x(is1: is1+dof1-1)
!
!                  print*, 'block:', elem1%blockST(loc_neigh)%Mb(1:dof1, 1:dof)

!                  print*, 'csdad' , is1, is1+dof1-1

                  accum(1:dof) = accum(1:dof) &
                       ! + matmul(elem%blockST(j)%Mb(1:dof, 1:dof1), x(is1: is1+dof1-1) ) - primal version
                      + matmul( x(is1: is1+dof1-1) , elem1%blockST(loc_neigh)%Mb(1:dof1, 1:dof) )! dual version
               endif
            enddo !j

            b(is: is+dof-1) = accum(1:dof)
         enddo !i

         if (big_dof /= nsize) then
            print*, 'wrong dimension in bMVprodST_dual' , big_dof, nsize
            stop
         endif
       else
         do i=1,grid%nelem
            accum(:) = 0.0
            elem => grid%elem(i)
            ndof1 = elem%dof
            Tdof = elem%Tdof
            ndof = elem%dof * ndim
            dof = ndof * Tdof

            is = elem%ncv
  !             do m =1,Tdof
  !                mm = (m-1) * ndof
  !                do n = 1,Tdof
  !                   nn = (n-1) * ndof
  !                   do k = 0,ndof-1,ndof1
  !                      !write(*,'(2(a7,i5),a6,i5,a1,i5,3(a6,i5))') &
  !                      !     'elem',i,'ncv:', is, 'accum',mm + k + 1,':', mm+ k + ndof1, 'k =', k , 'm=', m, 'n=', n
  !                      !	print*,   (is+nn+k )
  !                      accum(mm + k + 1 : mm+ k + ndof1) = accum(mm + k + 1 : mm+ k + ndof1) + &
  !                           eta * time%refTimeMatrix%Mb(m,n) * &
  !                           matmul(elem%Mass%Mb(1:ndof1, 1:ndof1) , x(is+nn+k : is+nn+k+ndof1-1))
  !
  !                   enddo !k
  !
  !                enddo !n
  !             enddo !m
  !
  !             ! diagonal block
  !             accum(1:dof) = accum(1:dof) &
  !                  + matmul(elem%blockST(0)%Mb(1:dof, 1:dof), x(is: is+dof-1) )
  !
           stop 'eta /= 0 not implemented in bMVprodST_Dual'

           !! off-diagonal blocks
           !F@R control offdiag blocks - structure
            do j=1,elem%flen
               k = grid%elem(i)%face(neigh,j)

               if(k > 0) then
                  elem1 => grid%elem(k)
                  dof1 = elem1%dof * ndim * elem1%Tdof ! same for dual and primal
                  is1 = elem1%ncv ! same for dual and primal
                  loc_neigh = elem%face(nei_i, j) ! local index of elem as an neighbor of elem1
                  accum(1:dof) = accum(1:dof) &
                       ! + matmul(elem%blockST(j)%Mb(1:dof, 1:dof1), x(is1: is1+dof1-1) ) - primal version
                      + matmul( x(is1: is1+dof1-1) , elem1%blockST(loc_neigh)%Mb(1:dof1, 1:dof) )! dual version
               endif
            enddo !j

            b(is: is+dof-1) = accum(1:dof)
         enddo !i

       endif ! eta

       deallocate(accum)
    class default
      stop 'bMVprodST_Dual only for STDG method'

    end select
    end associate

  end subroutine bMVprodST_Dual



  SUBROUTINE MGbMVprod(b,x,nsize)
  ! Block matrix-vector product in pMG cycle.
    use matrix_oper_int

    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    class(element), pointer:: elem,elem1 ! one element
!    integer :: i,j,k, ndof, ndof1, is, is1 !, l_accum
    !integer :: i,j,k, dof, ndof, is,is1, dof1,ndof1, is2
    integer:: i,j,k,ki,kj,dof,ndof,dof1,ndof1,MGdof,MGdof1,MGndof,MGndof1,MGis,MGis1
    real, dimension(:),allocatable :: accum
    !real,dimension(:,:),allocatable :: mtxdummy
    !real, dimension(:),allocatable, save:: accum

    ! allocate accum once to accomodate for the largest dof.
    allocate(accum(sum(grid%elem(:)%MGdof) * ndim ) )
    !l_accum= maxval(grid%elem%dof) * ndim
    !if(size(accum) <= l_accum) then
    !   deallocate(accum)
    !   allocate(accum(1:l_accum))
    !endif

    do  i=1,grid%nelem,1

        accum(:)=0.

        elem => grid%elem(i)

        MGdof = elem%MGdof  ! povodne ndof1
        MGndof = ndim*MGdof   != elem%MGdof * ndim
        MGis = elem%MGncv

        dof = elem%dof
        ndof = ndim*dof

        !! diagonal blocks
        if( eta /= 0.) then
            do  k = 0,ndof-1,dof
                accum(k+1:k+MGdof) = &
                    eta * matmul(elem%Mass%Mb(1:MGdof, 1:MGdof),x(MGis+k: Mgis+k+MGdof-1))
            end do
        end if

        !FIXME - nasobenie v rezime MG
        !accum(1:ndof) = accum(1:ndof) &
        !    + matmul(elem%block(0)%Mb(1:ndof, 1:ndof), x(is: is+ndof-1) )
        do  ki=0,ndim-1,1
            do  kj=0,ndim-1,1
                accum(ki*MGdof+1:(ki+1)*MGdof) = accum(ki*MGdof+1:(ki+1)*MGdof) &
                    + matmul( &
                    elem%block(0)%Mb(ki*dof+1:ki*dof+MGdof, kj*dof+1:kj*dof+MGdof), &
                    x(MGis+kj*MGdof: MGis+(kj+1)*MGdof-1) &
                    )
            end do
        end do


       !! off-diagonal blocks
       do j=1,elem%flen
          k = grid%elem(i)%face(neigh,j)

          if(k > 0) then
             elem1 => grid%elem(k)
             MGdof1 = elem1%MGdof
             MGndof1 = ndim*dof1
             MGis1 = elem1%MGncv
             dof1 = elem1%dof

             !FIXME - nasobenie v rezime MG
             !accum(1:ndof) = accum(1:ndof) &
             !     + matmul(elem%block(j)%Mb(1:ndof, 1:ndof1), x(is1: is1+ndof1-1) )
              do  ki=0,ndim-1,1
                  do  kj=0,ndim-1,1
                      accum(ki*MGdof+1:(ki+1)*MGdof) = accum(ki*MGdof+1:(ki+1)*MGdof) &
                          + matmul( &
                          elem%block(j)%Mb(ki*dof+1:ki*dof+MGdof, kj*dof1+1:kj*dof1+MGdof1), &
                          x(MGis1+kj*MGdof: MGis1+(kj+1)*MGdof1-1) &
                          )
                  end do
             end do
          endif
       enddo
       b(MGis: MGis+MGndof-1) = accum(1:MGndof)
    enddo

    deallocate(accum)

  END SUBROUTINE MGbMVprod

subroutine MGbMVprodDiag(b,x,nsize)
! Block matrix-vector product in pMG cycle.
    use matrix_oper_int
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    class(element), pointer:: elem
    integer :: i, ndof, is


    do i=1,grid%nelem
       elem => grid%elem(i)

       ndof = ndim * elem%MGdof  != elem%MGdof * ndim,  povodne ndof1
       is = elem%MGncv

       b(is: is+ndof-1) = matmul(elem%ILU(0)%Mb(1:ndof, 1:ndof), x(is: is+ndof-1) )
       !b(is: is+ndof-1) =  x(is: is+ndof-1)
    enddo


  end subroutine MGbMVprodDiag


  subroutine MGbMVprod2(b,x,nsize)
  ! Block matrix-vector product in pMG cycle.
    use matrix_oper_int
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    class(element), pointer:: elem,elem1 ! one element
!    integer :: i,j,k, ndof, ndof1, is, is1 !, l_accum
    integer :: i,j,k, dof, ndof, is,is1, dof1,ndof1, is2
   !integer:: i,j,k,dof,ndof,dof1,ndof1,MGdof,MGdof1,MGndof,MGndfo1,is1,is2
    real, dimension(:),allocatable :: accum
    !real, dimension(:),allocatable, save:: accum

    ! allocate accum once to accomodate for the largest dof.
    allocate(accum(maxval(grid%elem%dof) * ndim ) )
    !l_accum= maxval(grid%elem%dof) * ndim
    !if(size(accum) <= l_accum) then
    !   deallocate(accum)
    !   allocate(accum(1:l_accum))
    !endif

    accum(:)=0.

    do  i=1,grid%nelem,1
        elem => grid%elem(i)

        dof = elem%MGdof  ! povodne ndof1
        ndof = ndim*dof   != elem%MGdof * ndim
        is = elem%MGncv

        !! diagonal blocks
        if( eta /= 0.) then
            do  k = 0,ndof-1,dof
                accum(k+1:k+dof) = &
                    eta * matmul(elem%Mass%Mb(1:dof, 1:dof),x(is+k: is+k+dof-1))
            end do
        end if

        !FIXME - nasobeniev rezime MG
        !accum(1:ndof) = accum(1:ndof) &
        !    + matmul(elem%block(0)%Mb(1:ndof, 1:ndof), x(is: is+ndof-1) )
        do  k=0,ndof-1,dof
            accum(k+1:k+dof) = &
                matmul(elem%block(0)%Mb(1:dof, 1:dof), x(is+k: is+k+dof-1))
        end do


       !! off-diagonal blocks
       do j=1,elem%flen
          k = grid%elem(i)%face(neigh,j)

          if(k > 0) then
             elem1 => grid%elem(k)
             dof1 = elem1%MGdof
             ndof1 = dof1 * ndim
             is1 = elem1%MGncv

             !FIXME - nasobenie v rezime MG
             !accum(1:ndof) = accum(1:ndof) &
             !     + matmul(elem%block(j)%Mb(1:ndof, 1:ndof1), x(is1: is1+ndof1-1) )
              do  k=0,ndof-1,dof
                  accum(k+1:k+dof) = &
                      matmul(elem%block(j)%Mb(1:dof, 1:dof1), x(is1+k: is1+k+dof1-1))
             end do
          endif
       enddo
       b(is: is+ndof-1) = accum(1:ndof)
    enddo

    deallocate(accum)

  end subroutine MGbMVprod2


  !> double product by J. Hajek
  subroutine bMVprod2(b,x,nsize)
    use matrix_oper_int
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    real, dimension(1:nsize):: bm
    integer:: n

    n = nsize/2

    call bMVprod(b,x,n)
    call bMVprod(b(1+n),x(1+n),n)
    call bMVmassprod(bm,x,n)
    call bMVmassprod(bm(1+n),x(1+n),n)
    b(:n) = b(:n) + iota(1) * bm(n+1:)
    b(n+1:) = b(n+1:) + iota(2) * bm(:n)
  end subroutine bMVprod2

  !> diagonal block matrix - vector product: \f$ b = Mx \f$,  \f$ M  \f$ is the block
  !> mass matrix grid.elem(*).Mass
  subroutine bMVmassprod(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    class(element), pointer:: elem ! one element
    integer :: i, ndof, is

    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof = elem%dof * ndim
       is = elem%ncv

       ! FIXME
       ! diagonal block
       b(is: is+ndof-1) = matmul(elem%Mass%Mb, x(is: is+ndof-1) )
    enddo
  end subroutine bMVmassprod



  !> block matrix - vector product: \f$ b = Ax \f$,
  !>   \f$ A  \f$ is the block matrix grid.elem(*).block(*)
  subroutine bMVprodA(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    class(element), pointer:: elem,elem1 ! one element
    integer :: i,j,k,dof,  ndof, ndof1, is, is1

    real, dimension(:),allocatable:: accum

    ! allocate accum once to accomodate for the largest dof.
    allocate(accum(maxval(grid%elem%dof) * ndim))


    do i=1,grid%nelem
       elem => grid%elem(i)
       dof =  elem%dof
       ndof = dof * ndim
       is = elem%ncv

       ! diagonal block
       accum(1:ndof) = matmul(elem%block(0)%Mb(1:ndof, 1:ndof), x(is: is+ndof-1) )

       !! off-diagonal blocks
       do j=1,elem%flen
          k = grid%elem(i)%face(neigh,j)

          if(k > 0) then
             elem1 => grid%elem(k)
             ndof1 = elem1%dof * ndim
             is1 = elem1%ncv

             accum(1:ndof) = accum(1:ndof)  &
                  + matmul(elem%block(j)%Mb(1:ndof, 1:ndof1),  x(is1: is1+ndof1-1) )
          endif
       enddo
       b(is: is+ndof-1) = accum(1:ndof)
    enddo

    deallocate(accum)

  end subroutine bMVprodA

  !> block matrix - vector product: \f$ b = Ax \f$,
  !>   \f$ A  \f$ is the block matrix grid.elem(*).block(*)
  !> ONLY OFF DIAGONAL TERMS
  subroutine bMVprodOffC(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    class(element), pointer:: elem,elem1 ! one element
    integer :: i,j,k,dof,  ndof, ndof1, is, is1

    real, dimension(:),allocatable:: accum

    ! allocate accum once to accomodate for the largest dof.
    allocate(accum(maxval(grid%elem%dof) * ndim))


    !    print*,'B1'
    do i=1,grid%nelem
       elem => grid%elem(i)
       dof =  elem%dof
       ndof = dof * ndim
       is = elem%ncv

       accum(1:ndof) = 0.

       !! off-diagonal blocks
       do j=1,elem%flen
          k = grid%elem(i)%face(neigh,j)

          if(k > 0) then
             elem1 => grid%elem(k)
             ndof1 = elem1%dof * ndim
             is1 = elem1%ncv

             accum(1:ndof) = accum(1:ndof)  &
                  + matmul(elem%block(j)%Mb(1:ndof, 1:ndof1),  x(is1: is1+ndof1-1) )
          endif
       enddo
       b(is: is+ndof-1) = accum(1:ndof)
    enddo

    deallocate(accum)

  end subroutine bMVprodOffC


  !> evaluation of the Weighted Time Dependent residuum
  !> block matrix grid%elem(*)%block(*) + eta* grid%elem(*)%Mass(*)
  !> right-hand-side  b
  !> solution  x,
  !> res_new_old = (old residuum - new_residuum)
  subroutine EvalWeightTDresid(b, x, nsize, residuum, res_new_old, ires_no )
    use matrix_oper_int
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x, b
    real, intent(inout) :: residuum, res_new_old
    integer, intent(in) :: ires_no  ! = 1 => compute res_new_old
    class(element), pointer:: elem,elem1 ! one element
    integer :: i,j,k, ndof, ndof1, is, is1
    real, dimension(:),allocatable:: accum
    real :: max_res

    ! allocate accum once to accomodate for the largest dof.
    allocate(accum(maxval(grid%elem%dof) * ndim))

    max_res = 0.
    residuum = 0.
    res_new_old = 0.

    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof1 = elem%dof
       ndof = elem%dof * ndim
       is = elem%ncv

       ! diagonal block
       if (eta /= 0.) then
          do k = 0,ndof-1,ndof1
             accum(k+1:k+ndof1) = eta &
                  * matmul(elem%Mass%Mb(1:ndof1, 1:ndof1), x(is+k: is+k+ndof1-1))
          enddo
          accum(1:ndof) = accum(1:ndof) &
               + matmul(elem%block(0)%Mb(1:ndof, 1:ndof), x(is: is+ndof-1) )
       else
          accum(1:ndof) = matmul(elem%block(0)%Mb(1:ndof, 1:ndof), x(is: is+ndof-1) )
       endif


       !! off-diagonal blocks
       do j=1,elem%flen
          k = grid%elem(i)%face(neigh,j)

          if(k > 0) then
             elem1 => grid%elem(k)
             ndof1 = elem1%dof * ndim
             is1 = elem1%ncv

             accum(1:ndof) = accum(1:ndof) &
                  + matmul(elem%block(j)%Mb(1:ndof, 1:ndof1), x(is1: is1+ndof1-1) )
          endif
       enddo

       accum(1:ndof) = (accum(1:ndof) - b(is: is+ndof-1) )/ elem%area

       ! linear algebra residuum
       residuum = residuum + dot_product(accum(1:ndof), accum(1:ndof))

       ! difference od the old and new linear algebra residuum
       if(ires_no == 1) &
            res_new_old = res_new_old &
            + dot_product(accum(1:ndof) - elem%vec(res_vec,1:ndof), &
            accum(1:ndof) - elem%vec(res_vec,1:ndof))

       elem%vec(res_vec,1:ndof) = accum(1:ndof)

    enddo

    deallocate(accum)

    residuum = (residuum)**0.5 / state%space%domain_volume

    if(ires_no == 1) &
         res_new_old = (res_new_old)**0.5  / state%space%domain_volume

    !print*,'End of EvalWeightTDresid'

  end subroutine EvalWeightTDresid

  function EvalSSresid()
    real :: EvalSSresid
    class(element), pointer:: elem,elem1 ! one element
    integer :: i, j, k, ndof, ndof1

    real, dimension(:),allocatable:: accum

    ! allocate accum once to accomodate for the largest dof.
    allocate(accum(maxval(grid%elem%dof) * ndim))

    !    print*,'B1'
    EvalSSresid = 0.
    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof = elem%dof  * ndim

       ! diagonal block
       accum(1:ndof) = matmul(elem%block(0)%Mb(1:ndof, 1:ndof), elem%w(0,1:ndof) )

       !! off-diagonal blocks
       do j=1,elem%flen
          k = grid%elem(i)%face(neigh,j)

          if(k > 0) then
             elem1 => grid%elem(k)
             ndof1 = elem1%dof * ndim

             accum(1:ndof) = accum(1:ndof)  &
                  + matmul(elem%block(j)%Mb(1:ndof, 1:ndof1), elem1%w(0,1:ndof1) )
          endif
       enddo

       accum(1:ndof) = (accum(1:ndof) - elem%vec(rhs,1:ndof))/elem%area

       EvalSSresid = EvalSSresid + dot_product(accum(1:ndof), accum(1:ndof))


       !do j=1,ndim
       !write(*, '(a4,i5,20es12.4)') 'impl', i, accum(1:min(ndof,6) )

         !write(100+state%time%iter, '(2i5,20es12.4)') &
               !i, j, accum((j-1)*elem%dof+1 : j*elem%dof)
       !enddo
       !write(100+state%time%iter, '(i5,25es12.4)') i, elem%xc(1:nbDim), abs(accum(1:ndof))

    enddo


!    print*,'@@@',state%nsize, state%space%domain_volume, (EvalSSresid)**0.5 / state%space%domain_volume, &
!         (EvalSSresid/state%nsize)**0.5

    !EvalSSresid = (EvalSSresid)**0.5 / state%space%domain_volume
    EvalSSresid = (EvalSSresid/state%nsize)**0.5
    deallocate(accum)

  end function EvalSSresid

  function EvalSSresidExplicit()
    real :: EvalSSresidExplicit
    class(element), pointer:: elem ! one element
    integer :: i, ndof

    EvalSSresidExplicit = 0.

    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof = elem%dof  * ndim

       EvalSSresidExplicit = EvalSSresidExplicit &
            + dot_product(elem%vec(rhs,1:ndof), elem%vec(rhs,1:ndof))/elem%area**2

       !write(*, '(a4,i5,20es12.4)') 'expl',i, elem%vec(rhs,1:min(ndof, 6))/elem%area

    enddo

    !EvalSSresidExplicit = (EvalSSresidExplicit)**0.5 / state%space%domain_volume
    EvalSSresidExplicit = (EvalSSresidExplicit/state%nsize)**0.5

  end function EvalSSresidExplicit


  !> NO preconditioner: \f$ b = x \f$
  subroutine bMVnull(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b

    b(1:nsize) = x(1:nsize)
  end subroutine bMVnull

  !> diagonal block matrix - vector product: \f$ b = Ax \f$,  \f$ A  \f$ is the block
  !> matrix grid.elem(*).block(*)
  subroutine bMVdiagprod(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    class(element), pointer:: elem ! one element
    integer :: i, ndof, is

!    print*,'B1'
    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof = elem%dof * ndim
       is = elem%ncv

       ! diagonal block
       b(is: is+ndof-1) = matmul(elem%ILU(0)%Mb(1:ndof, 1:ndof), x(is: is+ndof-1) )

       !!call WriteMblock(elem%ILU(0) )
    enddo

  end subroutine bMVdiagprod

  !> diagonal block matrix - vector product for STDGM: \f$ b = Ax \f$,  \f$ A  \f$ is the block
  !> matrix grid.elem(*).block(*)
  subroutine bMVdiagprodST(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    class(element), pointer:: elem ! one element
    integer :: i, ndof, is

!    print*,'B1'
    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof = elem%dof * elem%Tdof * ndim
       is = elem%ncv

       ! diagonal block
       b(is: is+ndof-1) = matmul(elem%ILU(0)%Mb(1:ndof, 1:ndof), x(is: is+ndof-1) )

       !!call WriteMblock(elem%ILU(0) )
    enddo

  end subroutine bMVdiagprodST

    !> diagonal block matrix - vector product for Dual problem for STDGM: \f$ b = A^T * x \f$,  \f$ A  \f$ is the block
  !> matrix grid.elem(*).blockST(*)
  subroutine bMVdiagprodST_Dual(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    class(element), pointer:: elem ! one element
    integer :: i, ndof, iv

!    print*,'B1'
    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof = elem%dof * elem%Tdof * ndim
       iv = elem%ncv

       ! diagonal block - blockST is transposed for dual problem
       b(iv: iv+ndof-1) = matmul( x(iv: iv+ndof-1) , elem%ILU(0)%Mb(1:ndof, 1:ndof)  )

       !!call WriteMblock(elem%ILU(0) )
    enddo

  end subroutine bMVdiagprodST_Dual


  subroutine bMVdiagprod2(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    integer:: n
    n = nsize/2
    call bMVdiagprod(b,x,n)
    call bMVdiagprod(b(1+n),x(1+n),n)
  end subroutine bMVdiagprod2

  !> "preconditioned" norm of the vector \f$ \| x \|_P = \| Px\|_{\ell 2}\f$,
  !> where \f$ P \f$ is the actual matrix of ILU decomposition
  function VectorPrecondNorm(x)
    real :: VectorPrecondNorm
    real, dimension(:), intent(in) :: x
    real, dimension(:), allocatable :: y
    integer :: nsize

    nsize = size(x)
    allocate(y(1:nsize) )

!    if (state%dual) &
!      stop 'Vector precond norm not implemented for preconditioned Dual problem'

    if ( state%time%disc_time == 'STDG' ) then
       call bMViLUprodST( y, x, nsize )
    else
       call bMViLUprod(y, x ,nsize)
    endif
    VectorPrecondNorm = dot_product(y(:), y(:) )**0.5

    deallocate(y)

  end function VectorPrecondNorm

  !> "scaled" norm of the vector \f$ \| x \|_S = \| Sx\|_{\ell 2}\f$,
  !> where \f$ S \f$ is a scaling matrix
  function VectorScaleNorm(x)
    real :: VectorScaleNorm
    real, dimension(:), intent(in) :: x
    real, dimension(:), allocatable :: y
    class(element), pointer:: elem ! one element
    integer :: i, j, k, dof, is, is1
    real :: Re1

    Re1 = 0.
    if(state%model%Re > 0) Re1 = 1./state%model%Re
    allocate(y(1:state%space%max_dof) )

    VectorScaleNorm = 0.

    do i=1,grid%nelem
       elem => grid%elem(i)

       dof = elem%dof
       do k=1, ndim
          is  = elem%ncv + (k-1)*ndim
          is1 = is + dof - 1
          y(1:dof) = x(is:is1)  !/ (2*elem%area)
          do j=1,dof

             y(j) = y(j) / (elem%Mass%Mb(j,j) + elem%Stiff%Mb(j,j) * Re1 )**0.5
             !y(j) = y(j) / (elem%Mass%Mb(j,j) + elem%Stiff%Mb(j,j) )**0.5
          enddo
          VectorScaleNorm = VectorScaleNorm + dot_product( y(1:dof), y(1:dof) )

       enddo
    enddo

    !print*,'####',VectorScaleNorm, smaz, VectorScaleNorm / smaz

    VectorScaleNorm = VectorScaleNorm**0.5


    deallocate(y)

  end function VectorScaleNorm

  !> ILU preconditioning: \f$ b = (LU)^{-1}x \f$,  \f$ LU  \f$ is the incomplete LU block
  !> preconditioner having the same structure as matrix
  subroutine bMViLUprod(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    real, dimension(:), allocatable :: y
    class(element), pointer:: elem, elem1 ! one element
    type(Mblock) :: Loc

    integer :: i, ndof, is, j, i1, ndof1, is1, ii

    call InitMblock(Loc, grid%elem(1)%dof * ndim, grid%elem(1)%dof * ndim)

    allocate(y(1:nsize) )

    !! L solution
    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof = elem%dof * ndim
       is = elem%ncv

       y(is: is+ndof-1) = x(is: is+ndof-1)

       do j=1,elem%flen
          i1 = elem%face(neigh,j)
          if(i1 > 0 .and. i1 < i) then
             elem1 => grid%elem(i1)
             ndof1 = elem1%dof * ndim
             is1 = elem1%ncv

             y(is: is+ndof-1) = y(is: is+ndof-1) &
                  - matmul(elem%ILU(j)%Mb(1:ndof, 1:ndof1), y(is1: is1+ndof1-1) )

             if(is1 > is) print*,'%%%%%%%%%%%%%%%%%%%%%%%???',is,is1, i,i1

          endif
       enddo
    enddo

    !! U solution
    do ii=1,grid%nelem
       i = grid%nelem - ii + 1

       elem => grid%elem(i)
       ndof = elem%dof * ndim
       is = elem%ncv

       do j=1,elem%flen
          i1 = elem%face(neigh,j)

          if( i1 > i) then
             elem1 => grid%elem(i1)
             ndof1 = elem1%dof * ndim
             is1 = elem1%ncv

             y(is: is+ndof-1) = y(is: is+ndof-1) &
                  - matmul(elem%ILU(j)%Mb(1:ndof, 1:ndof1), b(is1: is1+ndof1-1) )

             if(is1 < is) print*,'UUU%%%%%%%%%%%%%%%%%%%%%%%???',is,is1, i,i1
          endif
       enddo

       if(ndof .ne. size(Loc%Mb,1)) then
          deallocate (Loc%Mb)
          call InitMblock(Loc, ndof, ndof)
       endif

       Loc%Mb(1:ndof,1:ndof) = grid%elem(i)%ILU(0)%Mb(1:ndof,1:ndof)
       call MblockInverse(ndof, Loc%Mb)

       b(is: is+ndof-1) = matmul(Loc%Mb(1:ndof,1:ndof), y(is: is+ndof-1) )

    enddo


    deallocate (Loc%Mb)
    deallocate(y)

  end subroutine bMViLUprod

  subroutine bMViLUprod2(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    integer:: n
    n = nsize/2
    call bMViLUprod(b,x,n)
    call bMViLUprod(b(1+n),x(1+n),n)
  end subroutine bMViLUprod2

  !> ILU preconditioning for STDGM: \f$ b = (LU)^{-1}x \f$,  \f$ LU  \f$ is the incomplete LU block
  !> preconditioner having the same structure as matrix
  subroutine bMViLUprodST(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    real, dimension(:), allocatable :: y
    class(element), pointer:: elem, elem1 ! one element
    type(Mblock) :: Loc
    integer :: i, ndof, is, j, i1, ndof1, is1, ii

    ndof = grid%elem(1)%dof * grid%elem(1)%Tdof * ndim
    call InitMblock(Loc, ndof, ndof)


    allocate( y(1:nsize) )

    !! L solution - L^{-1}: for elem (row) in elem%ILU(:) are saved the columns of this row
    do i=1,grid%nelem ! go through the rows of the matrix
       elem => grid%elem(i)
       ndof = elem%dof * elem%Tdof * ndim
       is = elem%ncv
       ! there is 1 on the diagonal of L^(-1)
       y(is: is+ndof-1) = x(is: is+ndof-1)

       ! go through the columns
       do j=1,elem%flen
          i1 = elem%face(neigh,j)
          if(i1 > 0 .and. i1 < i) then ! elem1 is before elem in the matrix
             elem1 => grid%elem(i1)
             ndof1 = elem1%dof * elem1%Tdof * ndim
             is1 = elem1%ncv

!             print*, 'ndof,ndof1, is1:' , ndof,ndof1, is1
!             print*, size( elem%ILU(j)%Mb(:,1) ) , size(elem%blockST(j)%Mb(:,1))
!             print*, elem%ILU(j)%Mb(1:ndof, 1:ndof1)
!             print*, 'dfs'
!             print*,y(is1: is1+ndof1-1)
!             print*, y(is: is+ndof-1)

             y(is: is+ndof-1) = y(is: is+ndof-1) &
                  - matmul(elem%ILU(j)%Mb(1:ndof, 1:ndof1), y(is1: is1+ndof1-1) )

             if(is1 > is) print*,'Problem in bMViLUProdST',is,is1, i,i1

          endif
       enddo
    enddo

    !! U solution
    do ii=1,grid%nelem
       i = grid%nelem - ii + 1

       elem => grid%elem(i)
       ndof = elem%dof * elem%Tdof * ndim
       is = elem%ncv

       do j=1,elem%flen
          i1 = elem%face(neigh,j)

          if( i1 > i) then
             elem1 => grid%elem(i1)
             ndof1 = elem1%dof * elem1%Tdof * ndim
             is1 = elem1%ncv

             y(is: is+ndof-1) = y(is: is+ndof-1) &
                  - matmul(elem%ILU(j)%Mb(1:ndof, 1:ndof1), b(is1: is1+ndof1-1) )

             if(is1 < is) print*,'Problem in bMViLUProdST',is,is1, i,i1
          endif
       enddo

       if(ndof .ne. size(Loc%Mb,1)) then
          deallocate (Loc%Mb)
          call InitMblock(Loc, ndof, ndof)
       endif

       Loc%Mb(1:ndof,1:ndof) = grid%elem(i)%ILU(0)%Mb(1:ndof,1:ndof)
       call MblockInverse(ndof, Loc%Mb)

       b(is: is+ndof-1) = matmul(Loc%Mb(1:ndof,1:ndof), y(is: is+ndof-1) )

    enddo


    deallocate (Loc%Mb)
    deallocate(y)

  end subroutine bMViLUprodST


  !> ILU preconditioning for STDGM DUAL PROBLEM: \f$ b = (LU)^{-T} x\f$,  \f$ LU  \f$ is the incomplete LU block
  !> we do not have L^{-1} neither , only implicitly as forward Gauss elimination)
  !> 1. solve: U^T z = x
  !> 2. multiply (G.e.): b = L^{-T} z
  !> the ordering of multiplication is reversed from the primal problem
  !> preconditioner having the same structure as matrix
  subroutine bMViLUprodST_Dual(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    real, dimension(:), allocatable :: z
    class(element), pointer:: elem, elem1 ! one element
    type(Mblock) :: Loc
    integer :: i, ii, ndof, dof, ncv, j, i1, ndof1, is, is1, loc_neigh

    ndof = grid%elem(1)%dof * grid%elem(1)%Tdof * ndim
    call InitMblock(Loc, ndof, ndof)

    allocate( z(1:nsize), source = 0.0 )


    !! U^Tz=x solution -> z
    do i=1,grid%nelem
       !i = ii ! grid%nelem - ii + 1
       elem => grid%elem(i)
       ndof = elem%dof * elem%Tdof * ndim
       is = elem%ncv
       z(is: is+ndof-1) = x(is:is+ndof-1)

       do j=1,elem%flen
          i1 = elem%face(neigh,j)

          if( i1>0 .and. i1 < i) then

             elem1 => grid%elem(i1)
             ndof1 = elem1%dof * elem1%Tdof * ndim
             is1 = elem1%ncv
             loc_neigh = elem%face(nei_i, j)
               ! zi = zi - u_ij^T zj
             z(is: is+ndof-1) = z(is: is+ndof-1) &
                - matmul( z(is1:is1+ndof1-1), elem1%ILU(loc_neigh)%Mb(1:ndof1, 1:ndof) )
!               - matmul(elem%ILU(j)%Mb(1:ndof, 1:ndof1), b(is1: is1+ndof1-1) )

             if(is1 > is) print*,'Problem in bMViLUProdST_Dual',is,is1, i,i1
          endif
       enddo

       if(ndof .ne. size(Loc%Mb,1)) then
          deallocate (Loc%Mb)
          call InitMblock(Loc, ndof, ndof)
       endif

       Loc%Mb(1:ndof,1:ndof) = grid%elem(i)%ILU(0)%Mb(1:ndof,1:ndof)
       call MblockInverse(ndof, Loc%Mb)
       ! z_i = z_i / u_ii
       z(is: is+ndof-1) = matmul( z(is: is+ndof-1), Loc%Mb(1:ndof,1:ndof) )

    enddo

    ! multiplication: b = L^{-T} z / we do not have L^{-1} => forward Gauss elimination of z
    do ii=1,grid%nelem ! go through the COLUMNS of the matrix (backwards)
       i = grid%nelem - ii + 1
       elem => grid%elem(i)
       ndof = elem%dof * elem%Tdof * ndim
       is = elem%ncv
       ! there is 1 on the diagonal of L^(-1)
       b(is: is+ndof-1) = z(is: is+ndof-1)

       ! go through the ROWS
       do j=1,elem%flen
          i1 = elem%face(neigh,j)
          if(i1 > i) then ! elem1 is below elem in the matrix L
             elem1 => grid%elem(i1)
             ndof1 = elem1%dof * elem1%Tdof * ndim
             is1 = elem1%ncv
             loc_neigh = elem%face(nei_i, j) ! local index of elem as an neighbor of elem1

             b(is: is+ndof-1) = b(is: is+ndof-1) &
               - matmul( b(is1: is1+ndof1-1), elem1%ILU(loc_neigh)%Mb(1:ndof1, 1:ndof) ) ! dual
!                  - matmul(elem%ILU(j)%Mb(1:ndof, 1:ndof1), y(is1: is1+ndof1-1) ) ! from primal

             if(is1 < is) print*,'Problem in bMViLUProdST',is,is1, i,i1

          endif
       enddo
    enddo


    deallocate (Loc%Mb)
    deallocate(z)

  end subroutine bMViLUprodST_Dual


  !>      null_precond  preconditioner
  subroutine null_precond(z,g,nsize,spars,irwst, idx)
    integer ::irwst(*), idx(*), nsize
    real z(*), g(*), spars(*)

    z(1:nsize)=g(1:nsize)
  end subroutine null_precond


  !>  Block Jacobi  preconditioner
  subroutine BlockJacobi(z,g,nsize,spars_prec,irwst, idx)
    integer :: i, j, iend, istart, irwst(*), idx(*), nsize
    real sum, z(*), g(*), spars_prec(*)


    !multiplication spars_prec*g

    do i=1,nsize
       istart=irwst(i)
       iend=irwst(i+1)-1
       !       print*,'###',i,istart,iend, j, idx(j),spars(j)
       sum=0.0
       do  j=istart,iend
          sum=sum+g(idx(j))*spars_prec(j)
          !print*,'###',i,istart,iend, j, idx(j),spars(j), sum
       enddo
       z(i)=sum
       !print*,'-----------------',i, sum
    enddo

  end subroutine BlockJacobi

  !> solution of \f$ (M+\tau C_k) x = b \f$ using Taylor serie,
  !> \f$ C_k = D_k + E_k \f$, $D_k is block diagonal part of \f$ C_k\f$
  subroutine TaylorSolution(nsize, x, b, tol, it, rezid, not_converge)
    !external prod
    integer, intent(in)  :: nsize
    integer, intent(out) :: it
    real, dimension(1:nsize), intent(inout) :: x, b
    real, intent(in) :: tol
    real, intent(out) :: rezid
    integer, intent(inout) :: not_converge   ! not converge = 1, converge = 0
    real, dimension(:,:), allocatable :: q
    integer :: itaylor

    allocate(q(1:nsize,1:nbDim) )

    not_converge = 0

!    call bMVprod(q(:,3), x(:), nsize)
!    q(:,3) = q(:,3) - b(:)
!    sum = dot_product(q(:,3) , q(:,3) )**0.5
!    print*,'reziduum',0,sum

    ! block multiplication $ q1 = (M+\tau D)^{-1} b $
    call bMVdiagprod(q(:,1), b(:), nsize)

    ! Taylor Series, zero-th step
    x(1:nsize) = q(:,1)

    ! iterative cycle
    do itaylor = 1, 500
       ! block matrix multiplication $ q1 = (M+\tau D)^{-1} E_k q1 $

       call bMVprodOffC(q(:,2), q(:,1), nsize)
       call bMVdiagprod(q(:,1), q(:,2), nsize)

       q(:,1) = - q(:,1)

       x(:) = x(:) + q(:,1)

       !call bMVprod(q(:,3), x(:), nsize)
       !q(:,3) = q(:,3) - b(:)
       !rezid = dot_product(q(:,3) , q(:,3) )**0.5

       rezid = (dot_product(q(:,1) , q(:,1) ) / dot_product(x(:) , x(:) ))**0.5
       !print*,'it rez sum',itaylor,rezid, sum

       if(rezid < tol) goto 100
       not_converge = 1
    enddo

100 continue

    it = itaylor

    deallocate( q )

  end subroutine TaylorSolution


  !> matrix - vector product for conforming FEM
  subroutine prodFEM(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    integer :: i, j

    b(:) = 0.
    do i=1,nsize
       do j= Mshape%irwst(i), Mshape%irwst(i+1)-1
          b(i) = b(i) + state%A(j) * x(Mshape%idx(j))
       enddo
    enddo

  end subroutine prodFEM

  !> diagonal preconditioner for conforming FEM
  subroutine diagFEM(b,x,nsize)
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b
    integer :: i, j

    b(:) = x(:)
    do i=1,nsize
       j= Mshape%irwst(i) ! first element is the diagonal one
       if(state%A(j) /= 0.)   b(i) = b(i) / state%A(j)

       if(Mshape%idx(j) /= i) then
          print*,'NOn-diagonal element in diagFEM !!!'
          stop
       end if
    enddo

  end subroutine diagFEM


  !> block matrix - vector product: \f$ b = (\eta M+ C_k) x \f$,
  !> \f$ M  \f$ is the block mass matrix grid.elem(*).Mass,
  !> \f$ C  \f$ is the block flux matrix grid.elem(*).block(*),
  !> \f$ \eta = 1/\tau_k\f$ is external
  subroutine bMVprod_SCHUR(b,x,nsize)
    use matrix_oper_int
    integer, intent (in):: nsize
    real, dimension(1:nsize), intent(in) :: x
    real, dimension(1:nsize), intent(inout) :: b

    b(1:nsize) = matmul(state%Schur_A(1:nsize, 1:nsize), x(1:nsize) )

  end subroutine bMVprod_SCHUR

  !> copy elem%blockST(0) to blockPlus - used for Ritz reconstruction
  subroutine CopyBlocksSTtoBlockPlus(  )
    class( element ), pointer :: elem
!    integer, intent(in) :: dof ! size of the block
    integer i,j, dof

    !print*, ' CopyBlocksSTtoBlockPlus called'

    do i = 1, grid%nelem
      elem => grid%elem(i)
      dof = size( elem%blockST(0)%Mb(:,1) )
      allocate( elem%blockPlus )
      call InitMblock( elem%blockPlus, dof, dof )
      elem%blockPlus%Mb(1:dof, 1:dof) = elem%blockST(0)%Mb(1:dof, 1:dof)
    end do !i

  end subroutine CopyBlocksSTtoBlockPlus

!  ! reduce the size of the blocks blockST to new size
!  subroutine reduceMatrixBlocks( elem, )

end module matrix_oper
