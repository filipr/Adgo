!> subroutines for linear algebra solvers
module lin_solvers
  use main_data
  use data_mod
  use matrix_oper
  use matrix_oper_int
  use gmres_solver
  use dir_solver
  use agmg_solver
  !  use pMultiGrid2

  implicit none

  public:: ComputeBlockDiagPrecond
  public:: ComputeBlockILUPrecond

  public:: ComputeBlockDiagPrecondST
  public:: ComputeBlockILUPrecondST

  public:: SolveBlockDiagonalSystem
  public:: SolveBlockLinearProblem
  public:: SolveBlockLinearDoubleProblem
  public:: SolveBlockLinearSTDGMProblem



contains

  !> evaluation of block diagonal preconditioner,
  !> LU decomposition of elem%block(0) using LAPACK
  subroutine ComputeBlockDiagPrecond()
    use matrix_oper_int
    integer :: i, k, dof, ndof

    do i=1,grid%nelem
       dof = grid%elem(i)%dof
       ndof = dof * ndim

       ! flux matrix
       grid%elem(i)%ILU(0)%Mb(1:ndof, 1:ndof) = grid%elem(i)%block(0)%Mb(1:ndof, 1:ndof)

       ! adding of mass matrix
       if (eta /= 0.) then
          do k = 0,ndof-1, dof
             grid%elem(i)%ILU(0)%Mb(k+1:k+dof, k+1:k+dof) = &
                  grid%elem(i)%ILU(0)%Mb(k+1:k+dof, k+1:k+dof) &
                  + eta * grid%elem(i)%Mass%Mb(1:dof,1:dof)
          enddo
       end if

       call MblockInverse(ndof, grid%elem(i)%ILU(0)%Mb(1:ndof, 1:ndof))

    enddo
  end subroutine ComputeBlockDiagPrecond

  !> evaluation of block diagonal preconditioner for STDG method,
  !> LU decomposition of elem%blockST(0) using LAPACK
  subroutine ComputeBlockDiagPrecondST()
    use matrix_oper_int
    type(Mblock) :: Temp
    class(element), pointer :: elem
    integer :: i, k, dof, Tdof, Sdof, ndof, m, n, im, in, ik

    associate( time => state%time )
      select type ( time)
         class is ( TimeTDG_t )

         dof = grid%elem(1)%dof * grid%elem(1)%Tdof
         call InitMblock(Temp, dof,dof)

         do i=1,grid%nelem
            elem => grid%elem(i)
            Sdof = elem%dof
            Tdof = elem%Tdof
            dof = Sdof * Tdof
            ndof = dof * ndim

            ! flux matrix
            elem%ILU(0)%Mb(1:ndof, 1:ndof) = elem%blockST(0)%Mb(1:ndof, 1:ndof)

            ! adding of mass matrix
            if (eta /= 0.) then

               ! bad variant
               !if(ndof .ne. size(Temp%Mb,1)) then
               !   deallocate (Temp%Mb)
               !   call InitMblock(Temp, dof, dof)
               !endif

               ! do m = 1,Tdof
               !    do n = 1,Tdof
               !       Temp%Mb((m-1)*Sdof + 1 : m*Sdof, (n-1)*Sdof + 1 : n*Sdof) = &
               !            eta * time%refTimeMatrix%Mb(m,n) * elem%Mass%Mb(1:Sdof,1:Sdof)
               !    enddo ! n
               ! enddo !m


               ! do k = 0,ndof-1, dof
               !    elem%ILU(0)%Mb(k+1:k+dof, k+1:k+dof) = &
               !         elem%ILU(0)%Mb(k+1:k+dof, k+1:k+dof) + Temp%Mb(1:dof,1:dof)
               ! enddo


               do m = 1,Tdof
                  do n = 1,Tdof
                     im = (m-1)* elem%dof * ndim
                     in = (n-1)* elem%dof * ndim
                     do k= 1, ndim
                        ik = (k-1)*Sdof

                        elem%ILU(0)%Mb(im + ik + 1: im + ik + Sdof, in + ik + 1: in + ik + Sdof) = &
                             elem%ILU(0)%Mb(im + ik + 1: im + ik + Sdof, in + ik + 1: in + ik + Sdof) &
                             + eta * time%refTimeMatrix%Mb(m,n) * elem%Mass%Mb(1:Sdof,1:Sdof)

                        !write(*,'(a6,6i5, a3,6i5)') '###WS',m,n,k,im,in, ik, '|', &
                        !     im + ik + 1, im + ik + Sdof, in + ik + 1, in + ik + Sdof
                     enddo

                  enddo ! n
               enddo !m
               !print*,'-------------------------------------'

            end if !eta

            !write(155,*), 'Element ', elem%i
            !do k = 1, ndof
            !   write(155,'(200es12.4)'), grid%elem(i)%ILU(0)%Mb(k, 1:ndof)
            !enddo !j

            call MblockInverse(ndof, grid%elem(i)%ILU(0)%Mb(1:ndof, 1:ndof))

         enddo

         !deallocate( Temp%Mb)

         !print *,'@@@@@@@@@@######'

         class default
         stop 'ComputeBlockDiagPrecondST this is for STDG only'


      end select
    end associate

  end subroutine ComputeBlockDiagPrecondST


  !> evaluation of block ILU preconditioner
  subroutine ComputeBlockILUPrecond()
    use matrix_oper_int
    type(Mblock) :: Loc
    !class(element), pointer:: elem  ! elem = element, elem1 = neighbour element
    integer :: i, in, j, k, dof, ndof

    dof = grid%elem(1)%dof * ndim
    call InitMblock(Loc, dof, dof)

    ! BLOCK ILU decomposition
    do i=1,grid%nelem
       dof = grid%elem(i)%dof
       ndof = dof * ndim

       ! diagonal term
       ! flux matrix
       grid%elem(i)%ILU(0)%Mb(1:ndof, 1:ndof) = grid%elem(i)%block(0)%Mb(1:ndof, 1:ndof)

       ! adding of mass matrix
       if (eta /= 0.) then
          do k = 0,ndof-1, dof
             grid%elem(i)%ILU(0)%Mb(k+1:k+dof, k+1:k+dof) = &
                  grid%elem(i)%ILU(0)%Mb(k+1:k+dof, k+1:k+dof) &
                  + eta * grid%elem(i)%Mass%Mb(1:dof,1:dof)
          enddo
       end if

       do k=1,grid%elem(i)%flen
          in = grid%elem(i)%face(neigh,k)

          if(in > 0 .and. in < i) then
             grid%elem(i)%ILU(0)%Mb(:,:) = grid%elem(i)%ILU(0)%Mb(:,:) &
                  - matmul(grid%elem(i)%ILU(k)%Mb(:,:), &
                  grid%elem(in)%ILU(grid%elem(i)%face(nei_i,k) )%Mb(:,:) )
          endif
       enddo

       ! off diagonal terms in the row
       do k=1,grid%elem(i)%flen
          in = grid%elem(i)%face(neigh,k)
          if(in > i) then
             grid%elem(i)%ILU(k)%Mb(:,:) = grid%elem(i)%block(k)%Mb(:,:)
          endif
       enddo


       ! off diagonal terms in the columns
       !dof = grid%elem(i)%dof * ndim   ERROR ????????????????
       ndof = grid%elem(i)%dof * ndim
       if(ndof .ne. size(Loc%Mb,1)) then
          deallocate (Loc%Mb)
          call InitMblock(Loc, ndof, ndof)
       endif

       !call MblockLU( grid%elem(i)%ILU(0), Loc )
       Loc%Mb(1:ndof,1:ndof) = grid%elem(i)%ILU(0)%Mb(1:ndof,1:ndof)

       call MblockInverse(ndof, Loc%Mb)

       do k=1,grid%elem(i)%flen
          in = grid%elem(i)%face(neigh,k)
          j = grid%elem(i)%face(nei_i,k)
          if(in > i) then

             !if(size( grid%elem(in)%ILU(j)%Mb, 1) /= size( grid%elem(in)%ILU(j)%Mb, 2))then
             !   write(*,'(a6,8i5)') '####', ndof, dof, &
             !        size( grid%elem(in)%ILU(j)%Mb, 1), &
             !        size( grid%elem(in)%ILU(j)%Mb, 2), &
             !        size( grid%elem(in)%block(j)%Mb, 1), &
             !        size( grid%elem(in)%block(j)%Mb, 2), &
             !        size(Loc%Mb, 1), size(Loc%Mb, 2)
             !endif


             grid%elem(in)%ILU(j)%Mb(:,:) &
                  = matmul(grid%elem(in)%block(j)%Mb(:,:), Loc%Mb(:,:) )


             !grid%elem(in)%ILU(j)%Mb(1:ndof,1:ndof) &
             !     = matmul(grid%elem(in)%block(j)%Mb(1:ndof,1:ndof), Loc%Mb(1:ndof,1:ndof) )
          endif
       enddo

       !
       !if(i .eq. 1891) then
       !if(state%time%iter >= 86 ) then
       !!   call WriteMblock(grid%elem(i)%block(0))
       !   call WriteMblock(grid%elem(i)%ILU(0))
       !   do k=1,grid%elem(i)%flen
       !      in = grid%elem(i)%face(neigh,k)
       !      if(in > 0 ) then
       !!       !  call WriteMblock(grid%elem(i)%block(k))
       !         call WriteMblock(grid%elem(i)%ILU(k))
       !         print*,'*********',i,k,in, grid%elem(in)%xc(:)
       !      endif
       !   enddo
       !   print*,'*******************',grid%elem(i)%i,'** ', grid%elem(i)%xc(:)
       !endif

    enddo

    deallocate (Loc%Mb)

    !print*,'ILU decomposition computed'
    !if(state%time%iter >= 86) stop

  end subroutine ComputeBlockILUPrecond


  !> evaluation of block ILU preconditioner
  subroutine ComputeBlockMGILUPrecond()
    use matrix_oper_int
    type(Mblock) :: Loc
    !class(element), pointer:: elem  ! elem = element, elem1 = neighbour element
    integer :: i, in, j, k, dof, ndof

    dof = grid%elem(1)%MGdof * ndim
    call InitMblock(Loc, dof, dof)

    ! BLOCK ILU decomposition
    do i=1,grid%nelem,1
       dof = grid%elem(i)%MGdof
       ndof = dof * ndim

       ! diagonal term
       ! flux matrix
       grid%elem(i)%XXX(0)%Mb(1:ndof, 1:ndof) = grid%elem(i)%block(0)%Mb(1:ndof, 1:ndof)

       ! adding of mass matrix
       if (eta /= 0.) then
          do k = 0,ndof-1, dof
             grid%elem(i)%XXX(0)%Mb(k+1:k+dof, k+1:k+dof) = &
                  grid%elem(i)%XXX(0)%Mb(k+1:k+dof, k+1:k+dof) &
                  + eta * grid%elem(i)%Mass%Mb(1:dof,1:dof)
          enddo
       end if

       do k=1,grid%elem(i)%flen
          in = grid%elem(i)%face(neigh,k)

          if(in > 0 .and. in < i) then
             grid%elem(i)%XXX(0)%Mb(:,:) = grid%elem(i)%XXX(0)%Mb(:,:) &
                  - matmul(grid%elem(i)%XXX(k)%Mb(:,:), &
                  grid%elem(in)%XXX(grid%elem(i)%face(nei_i,k) )%Mb(:,:) )
          endif
       enddo

       ! off diagonal terms in the row
       do k=1,grid%elem(i)%flen

          in = grid%elem(i)%face(neigh,k)
          if(in > i) then
             grid%elem(i)%XXX(k)%Mb(:,:) = grid%elem(i)%block(k)%Mb(:,:)
          endif
       enddo


       ! off diagonal terms in the columns
       !dof = grid%elem(i)%dof * ndim   ERROR ????????????????
       ndof = grid%elem(i)%MGdof * ndim
       if(ndof .ne. size(Loc%Mb,1)) then
          deallocate (Loc%Mb)
          call InitMblock(Loc, ndof, ndof)
       endif

       !call MblockLU( grid%elem(i)%ILU(0), Loc )
       Loc%Mb(1:ndof,1:ndof) = grid%elem(i)%XXX(0)%Mb(1:ndof,1:ndof)

       call MblockInverse(ndof, Loc%Mb)

       do k=1,grid%elem(i)%flen
          in = grid%elem(i)%face(neigh,k)
          j = grid%elem(i)%face(nei_i,k)
          if(in > i) then

             !if(size( grid%elem(in)%ILU(j)%Mb, 1) /= size( grid%elem(in)%ILU(j)%Mb, 2))then
             !   write(*,'(a6,8i5)') '####', ndof, dof, &
             !        size( grid%elem(in)%ILU(j)%Mb, 1), &
             !        size( grid%elem(in)%ILU(j)%Mb, 2), &
             !        size( grid%elem(in)%block(j)%Mb, 1), &
             !        size( grid%elem(in)%block(j)%Mb, 2), &
             !        size(Loc%Mb, 1), size(Loc%Mb, 2)
             !endif


             grid%elem(in)%XXX(j)%Mb(:,:) &
                  = matmul(grid%elem(in)%block(j)%Mb(:,:), Loc%Mb(:,:) )


             !grid%elem(in)%ILU(j)%Mb(1:ndof,1:ndof) &
             !     = matmul(grid%elem(in)%block(j)%Mb(1:ndof,1:ndof), Loc%Mb(1:ndof,1:ndof) )
          endif
       enddo

       !
       !if(i .eq. 1891) then
       !if(state%time%iter >= 86 ) then
       !!   call WriteMblock(grid%elem(i)%block(0))
       !   call WriteMblock(grid%elem(i)%ILU(0))
       !   do k=1,grid%elem(i)%flen
       !      in = grid%elem(i)%face(neigh,k)
       !      if(in > 0 ) then
       !!       !  call WriteMblock(grid%elem(i)%block(k))
       !         call WriteMblock(grid%elem(i)%ILU(k))
       !         print*,'*********',i,k,in, grid%elem(in)%xc(:)
       !      endif
       !   enddo
       !   print*,'*******************',grid%elem(i)%i,'** ', grid%elem(i)%xc(:)
       !endif

    enddo

    deallocate (Loc%Mb)

    !print*,'ILU decomposition computed'
    !if(state%time%iter >= 86) stop

  end subroutine ComputeBlockMGILUPrecond


  subroutine MatmulLU
    INTEGER :: N,i,j,k,j4i,k4i,j4k,krow,kcol
    REAL,dimension(:,:),allocatable :: Mtx
    INTEGER,dimension(:,:),allocatable  :: nnzirow,nnzjcol
    INTEGER:: nnzilen, nnzjlen

    allocate( Mtx(1,1) )
    allocate( nnzirow(2,4), nnzjcol(2,4) )
    call ComputeBlockILUPrecond
    call ComputeBlockXXXPrecond

    N = grid%nelem

    do i=1,N,1
       print*, ''
       ! Aii nonzero diagonal block

       ! finding nonzero off-diagonal blocks
       do j4i=1,grid%elem(i)%flen,1
          j=grid%elem(i)%face( neigh,j4i )
          if ( (0<j) ) then ! Aij nonzero off-diagonal block

             if ( (size(Mtx,1).NE.ndim*grid%elem(i)%dof).AND.(size(Mtx,2).NE.ndim*grid%elem(j)%dof) ) then
                deallocate(Mtx)
                allocate( Mtx(ndim*grid%elem(i)%dof,ndim*grid%elem(j)%dof) )
             end if
             Mtx=0.

             print*,'Aij',i,j
             !
             if ( size(nnzirow,1).LT.grid%elem(i)%flen+2 ) then
                deallocate(nnzirow)
                allocate( nnzirow(2,grid%elem(i)%flen+2) )
             end if
             !
             nnzirow(1,1:1)=(/i/)
             nnzirow(2,1:1)=(/0/)
             nnzilen = 1
             do k=1,grid%elem(i)%flen,1
                nnzilen = nnzilen+1
                nnzirow( 1,nnzilen ) = grid%elem(i)%face(neigh,k)
                nnzirow( 2,nnzilen ) = k
                if ( nnzirow(1,nnzilen)<0 ) then
                   nnzilen = nnzilen-1
                end if
             end do
             !
             if ( size(nnzjcol,1).LT.grid%elem(j)%flen+2 ) then
                deallocate(nnzjcol)
                allocate( nnzjcol(2,grid%elem(j)%flen+2) )
             end if
             !
             nnzjcol(1,1:1)=(/j/)
             nnzjcol(2,1:1)=(/0/)
             nnzjlen = 1
             do k=1,grid%elem(j)%flen,1
                nnzjlen = nnzjlen+1
                nnzjcol( 1,nnzjlen ) = grid%elem(j)%face(neigh,k)
                nnzjcol( 2,nnzjlen ) = k
                if ( nnzjcol(1,nnzjlen)<0 ) then
                   nnzjlen = nnzjlen-1
                end if
             end do
             !
             !
             do krow=1,nnzilen,1
                do kcol=1,nnzjlen,1
                   if (  nnzirow(1,krow).EQ.nnzjcol(1,kcol) ) then
                      k = nnzirow(1,krow)
                      if ( k.LE.min(i,j) ) then
                         if ( k.LT.i) then
                            print*,'+Lik*Ukj',i,k,j

                            if ( nnzjcol(2,kcol).NE.0 ) then
  Mtx = Mtx+matmul(grid%elem(i)%XXX( nnzirow(2,krow) )%Mb, &
       grid%elem( grid%elem(j)%face( neigh,nnzjcol(2,kcol) )  )%XXX( grid%elem(j)%face( nei_i,nnzjcol(2,kcol) ) )%Mb)
                            else
                               Mtx = Mtx+matmul(grid%elem(i)%XXX( nnzirow(2,krow) )%Mb, &
                                    grid%elem( j )%ILU( 0 )%Mb)
                            end if
                         elseif ( k.EQ.i )  then
                            print*,'+Uij',i,j
                            print*,size(Mtx,1),ndim*grid%elem(i)%dof,size(grid%elem(j)%XXX( nnzjcol(2,kcol) )%Mb,1)
                            print*,size(Mtx,2),ndim*grid%elem(j)%dof,size(grid%elem(j)%XXX( nnzjcol(2,kcol) )%Mb,2)
 Mtx = Mtx+grid%elem( grid%elem(j)%face( neigh,nnzjcol(2,kcol) ) )%XXX( grid%elem(j)%face( nei_i,nnzjcol(2,kcol) ) )%Mb
                         end if
                      end if
                   end if
                end do
             end do
             do krow=1,size(Mtx,1)
                do kcol=1,size(Mtx,2)
                   print*,Mtx(krow,kcol), grid%elem(i)%block(j4i)%Mb(krow,kcol)
                end do
             end do
          end if
       end do
    end do

  end subroutine MatmulLU

  subroutine CheckBlockILU
    INTEGER :: i,jid,j,k,N,nrow,ncol,dof
    REAL,dimension(:,:),allocatable :: Mtx,MtxA
    INTEGER,dimension(:),allocatable :: nnzi,nnzj

    call ComputeBlockILUPrecond
    call ComputeBlockXXXPrecond
    !
    N = grid%nelem
    !
    allocate( Mtx(1,1), MtxA(1,1) )
    !
    ! *** DIAGONAL *** blocks, i.e. Aii
    do i=1,N,1
       dof = grid%elem(i)%dof
       nrow = ndim*dof
       ncol = nrow
       if ( (size(Mtx,1).NE.nrow) .AND. (size(Mtx,2).NE.ncol) ) then
          deallocate( Mtx, MtxA )
          allocate( Mtx(nrow,ncol), MtxA(nrow,ncol) )
       end if
       !
       ! product for diagonal blocks
       ! Lii*Uii = Uii = Aii
       Mtx = grid%elem(i)%ILU(0)%Mb(:,:)
       !
       ! product for off-diagonal blocks
       do jid=1,grid%elem(i)%flen,1
          j = grid%elem(i)%face( neigh,jid )
          if ( (0<j).AND.(j<i) ) then
             print*,'j<i',j,i,jid,grid%elem(i)%face(nei_i,jid)
             Mtx = Mtx+matmul(grid%elem(i)%ILU(jid)%Mb(:,:), grid%elem(j)%ILU( grid%elem(i)%face(nei_i,jid) )%Mb(:,:))
          end if
       end do
       !
       ! Control ("original") matrix
       MtxA = grid%elem(i)%block(0)%Mb(:,:)
       if ( eta /= 0. ) then
          do k=0,ndim*dof-1,dof
             MtxA(k+1:k+dof, k+1:k+dof) = MtxA(k+1:k+dof, k+1:k+dof) &
                  + eta*grid%elem(i)%Mass%Mb(1:dof, 1:dof)
          end do
       end if
       !
       ! Comparision
       print*,'Nonzero block should be "identical"'
       do j=1,nrow,1
          do k=1,ncol,1
             print*,Mtx(j,k),'?=',MtxA(j,k) !, grid%elem(1)%ILU(0)%Mb(j,k), grid%elem(1)%XXX(0)%Mb(j,k),grid%elem(i)%block(0)%Mb(j,k)
          end do
       end do
    end do

    !  *** OFF-DIAGONAL *** blocks (on nonzero-block position)
    do i=1,N,1
       do jid=1,grid%elem(i)%flen,1
          j = grid%elem(i)%face( neigh,jid )
          if ( (0<j).AND.(j<i) ) then

          end if
       end do
    end do

    i=1


    ! off-diagonal blocks, i.e. Aij, i/=j
    do jid=1,grid%elem(i)%flen,1
       ! evaluate Aij= SUM(k=1:N) Lik*Ukj
       ! observation Lii*Uij=Uij
       !
       ! Dimension of block Aij is nrow-times-ncol
       j = grid%elem(i)%face( neigh,jid )
       !
       !print*,'i=',i,'j=',j
       if ( j.gt.0 ) then
          nrow = ndim*grid%elem(i)%dof
          ncol = ndim*grid%elem(j)%dof
          if ( (size(Mtx,1).NE.nrow) .AND. (size(Mtx,2).NE.ncol) ) then
             deallocate( Mtx )
             allocate( Mtx(nrow,ncol) )
             Mtx = 0.
          end if
          print*,size(grid%elem(i)%ILU(jid)%Mb(:,:),1), size(grid%elem(i)%ILU(jid)%Mb(:,:),2), &
               size(grid%elem(j)%ILU( grid%elem(i)%face(nei_i,jid) )%Mb(:,:),1), &
               size(grid%elem(j)%ILU( grid%elem(i)%face(nei_i,jid) )%Mb(:,:),2), &
               size(Mtx,1), size(Mtx,2)
          Mtx = Mtx+matmul(grid%elem(i)%ILU(jid)%Mb(:,:), grid%elem(j)%ILU( grid%elem(i)%face(nei_i,jid) )%Mb(:,:))
       end if
    end do
  end subroutine CheckBlockILU


  SUBROUTINE ComputeBlockXXXPrecond
    use matrix_oper_int
    type(Mblock) :: Mtx
    !class(element), pointer:: elem  ! elem = element, elem1 = neighbour element
    INTEGER :: i,j,k,kid,jid,dof,ndof
    INTEGER :: N, NNN
    !
    N=grid%nelem
    !
    ndof = ndim*grid%elem(1)%dof
    call InitMblock( Mtx, ndof, ndof )
    !
    ! Setup matrix blocks
    do i=1,N,1
       dof = grid%elem(i)%dof
       ndof = ndim*dof
       !
       ! diagonal blocks
       grid%elem(i)%XXX(0)%Mb(1:ndof,1:ndof) = grid%elem(i)%block(0)%Mb(1:ndof,1:ndof)
       !
       if ( eta/=0. ) then
          do k=0,ndof-1,dof
             grid%elem(i)%XXX(0)%Mb(k+1:k+dof, k+1:k+dof) &
                  = grid%elem(i)%XXX(0)%Mb(k+1:k+dof, k+1:k+dof) &
                  + eta*grid%elem(i)%Mass%Mb(1:dof, 1:dof)
          end do
       end if
       !
       ! off-diagonal blocks
       do k=1,grid%elem(i)%flen,1
          !kid = grid%elem(i)%face(nei_i,k)
          if ( grid%elem(i)%face(neigh,k)>0 ) then
             grid%elem(i)%XXX( k )%Mb(:,:) = grid%elem(i)%block( k )%Mb(:,:)
          end if
       end do
    end do

    return
    print*,'pokracujem'
    !
    ! ILU(0) algorithm done in place
    do i=2,N,1
       ! for k=1 to (i-1)
       do kid=1,grid%elem(i)%flen
          k = grid%elem(i)%face(neigh,kid)
          if ( 0<k .AND. k<i ) then
             ! Aik=Aik/Akk
             !
             NNN = ndim*grid%elem(k)%dof
             if ( NNN .NE. size(Mtx%Mb,1) ) then
                deallocate( Mtx%Mb )
                call InitMblock( Mtx, NNN, NNN)
             end if
             !
             Mtx%Mb(1:NNN, 1:NNN ) = grid%elem(k)%XXX(0)%Mb( 1:NNN, 1:NNN )
             call MblockInverse( ndim*grid%elem(k)%dof, Mtx%Mb )
             !
             grid%elem(i)%XXX( kid )%Mb(1:NNN,1:NNN) &
                  = matmul( grid%elem(i)%XXX( kid )%Mb(1:NNN,1:NNN) , Mtx%Mb(1:NNN, 1:NNN) )
             !
             ! for j=(k+1) to N
             do jid=1,grid%elem(i)%flen
                j = grid%elem(i)%face(neigh,jid)
                if ( k<j .AND. j<(N+1) ) then
                   ! Aij=Aij-Aik*Akj
                   grid%elem(i)%XXX( jid )%Mb(:,:) &
                        = grid%elem(i)%XXX( jid )%Mb(:,:) &
                        - matmul( &
                        grid%elem(i)%XXX( kid )%Mb(:,:) &
                        ,grid%elem(k)%XXX( grid%elem(i)%face(nei_i,kid) )%Mb(:,:) &
                        )
                end if
             end do
          end if
       end do
    end do
  END SUBROUTINE ComputeBlockXXXPrecond


  !> evaluation of block ILU preconditioner
  subroutine ComputeBlockILUPrecondST()
    use matrix_oper_int
    type(Mblock) :: Loc
    type(Mblock) :: Temp
    class(element), pointer :: elem
    integer :: i, in, j, k, dof, ndof , Tdof, Sdof, m, n, im, ik


    associate( time => state%time )
    select type ( time)
       class is ( TimeTDG_t )

       dof = grid%elem(1)%dof * grid%elem(1)%Tdof
       !call InitMblock( Temp, dof, dof )

       dof = dof * ndim
       call InitMblock(Loc, dof, dof)

       ! BLOCK ILU decomposition
       do i=1,grid%nelem
          elem => grid%elem(i)
          Sdof = elem%dof
          Tdof = elem%Tdof
          dof = Sdof * Tdof
          ndof = dof * ndim

          ! diagonal term
          ! flux matrix
          grid%elem(i)%ILU(0)%Mb(1:ndof, 1:ndof) = grid%elem(i)%blockST(0)%Mb(1:ndof, 1:ndof)

          ! adding of mass matrix
          if (eta /= 0.) then

             do m = 1,Tdof
                do n = 1,Tdof
                   im = (m-1)* elem%dof * ndim
                   in = (n-1)* elem%dof * ndim
                   do k= 1, ndim
                      ik = (k-1)*Sdof

                      elem%ILU(0)%Mb(im + ik + 1: im + ik + Sdof, in + ik + 1: in + ik + Sdof) = &
                           elem%ILU(0)%Mb(im + ik + 1: im + ik + Sdof, in + ik + 1: in + ik + Sdof) &
                           + eta * time%refTimeMatrix%Mb(m,n) * elem%Mass%Mb(1:Sdof,1:Sdof)

                      !write(*,'(a6,6i5, a3,6i5)') '###WS',m,n,k,im,in, ik, '|', &
                      !     im + ik + 1, im + ik + Sdof, in + ik + 1, in + ik + Sdof
                   enddo

                enddo ! n
             enddo !m

             ! if(ndof .ne. size(Temp%Mb,1)) then
             !    deallocate (Temp%Mb)
             !    call InitMblock(Temp, dof, dof)
             ! endif


             ! do m = 1,Tdof
             !    do n = 1,Tdof
             !       Temp%Mb((m-1)*Sdof + 1 : m*Sdof, (n-1)*Sdof + 1 : n*Sdof) = &
             !            eta * time%refTimeMatrix%Mb(m,n) * elem%Mass%Mb(1:Sdof,1:Sdof)
             !    enddo ! n
             ! enddo !m

             ! do k = 0,ndof-1, dof
             !    elem%ILU(0)%Mb(k+1:k+dof, k+1:k+dof) = &
             !         elem%ILU(0)%Mb(k+1:k+dof, k+1:k+dof) + Temp%Mb(1:dof,1:dof)
             ! enddo

          end if !eta


          do k=1,grid%elem(i)%flen
             in = grid%elem(i)%face(neigh,k)

             if(in > 0 .and. in < i) then
                grid%elem(i)%ILU(0)%Mb(:,:) = grid%elem(i)%ILU(0)%Mb(:,:) &
                     - matmul(grid%elem(i)%ILU(k)%Mb(:,:), &
                     grid%elem(in)%ILU(grid%elem(i)%face(nei_i,k) )%Mb(:,:) )

             endif
          enddo


          ! off diagonal terms in the row
          do k=1,grid%elem(i)%flen
             in = grid%elem(i)%face(neigh,k)
             if(in > i) then
                grid%elem(i)%ILU(k)%Mb(:,:) = grid%elem(i)%blockST(k)%Mb(:,:)
             endif
          enddo


          ! off diagonal terms in the columns
          !dof = grid%elem(i)%dof * ndim
          !ndof = grid%elem(i)%dof * ndim
          if(ndof .ne. size(Loc%Mb,1)) then
             deallocate (Loc%Mb)
             call InitMblock(Loc, ndof, ndof)
          endif

          !call MblockLU( grid%elem(i)%ILU(0), Loc )
          Loc%Mb(1:ndof,1:ndof) = elem%ILU(0)%Mb(1:ndof,1:ndof)

          call MblockInverse(ndof, Loc%Mb)

          do k = 1, elem%flen
             in = elem%face(neigh,k)
             j = elem%face(nei_i,k)

             if(in > i) then
                grid%elem(in)%ILU(j)%Mb(:,:) &
                     = matmul(grid%elem(in)%blockST(j)%Mb(:,:), Loc%Mb(:,:) )
             endif
          enddo !k

       enddo !i

       deallocate (Loc%Mb)
       !deallocate (Temp%Mb)

       !print*,'ILU decomposition computed'
       !if(time%iter >= 86) stop
       class default
       stop 'For STDGM only'
    end select
  end associate


end subroutine ComputeBlockILUPrecondST



!> solution of the block diagonal system
!> \f$ \frac{1}{\tau_k} {\bf M}{\bf w}^k = {\bf q}({\bf w}^{k-1}) \f$,
!> "steady-state residuum" norm_res
!> = \f$ \left(\frac{1}{dof} \sum_{i=1}^{dof} fluxes_i({\bf w})^2 \right)^{1/2} \f$
subroutine SolveBlockDiagonalSystem(eta, norm_res)
real, intent(in) :: eta
real, intent(out) :: norm_res
class(element), pointer:: elem ! one element
!real :: res2
integer  :: ie, k, dof,  i1, i2

norm_res = 0.

do ie=1,grid%nelem
   elem => grid%elem(ie)
   dof = elem%dof

   do k=1,ndim
      i1 = (k-1)*dof + 1
      i2 = k*dof

      elem%w(0,i1:i2) = matmul(elem%MassInv%Mb(1:dof,1:dof), &
           elem%vec(rhs, i1:i2)/eta + elem%vec(rhsM, i1:i2) )

      norm_res = norm_res + dot_product( elem%vec(rhs, i1:i2), elem%vec(rhs, i1:i2))

   end do
end do

norm_res = (norm_res / state%nsize)**0.5

end subroutine SolveBlockDiagonalSystem


! subroutine SolveBlockLinearProblemAAA(nsize, eta, b, x, rezid, tot_iter, &
!       precond_update, not_conv)

!   if(state%space%adapt%adapt_method /= 'ALG2') then
!      call SolveBlockLinearProblem(state%nsize, eta, Newton%b, Newton%x, &
!           state%linSolver%residuum, state%linSolver%iter, precond_update, state%linSolver%lin_solver_not_conv)
!   else

!      do i=1, max_i
!         state%first_GMRES_conv = ??

!          call SolveBlockLinearProblem(state%nsize, eta, Newton%b, Newton%x, &
!           state%linSolver%residuum, state%linSolver%iter, precond_update, state%linSolver%lin_solver_not_conv)

!          estim Alg error

!   endif

! end subroutine SolveBlockLinearProblemAAA

!> solution of large block linear algebraic problem \f$ (A + \eta\,M)x= b\f$
!>
!> \f$ A\f$ is a sparse block matrix given by elem%block(0:len),
!> \f$ x\f$ is vecor with initial guess of the solution (in) and the solution (out),
!> \f$ b\f$ is the righ-hand side
subroutine SolveBlockLinearProblem(nsize, eta, b, x, rezid, tot_iter, &
     precond_update, not_conv)
  use matrix_oper_int, mx_eta => eta
  !    use pMultiGrid, only: L2res, SolveViaUMFPACK
  integer :: nsize                               ! size of the algebraic problem
  real, intent(in):: eta
  real, dimension(1:nsize), intent(inout):: x    ! solution
  real, dimension(1:nsize), intent(inout):: b    ! RHS
  real, intent(inout) :: rezid                   ! reziduum
  integer, intent(inout) :: tot_iter             ! number of iterations
  logical, intent(in) :: precond_update          ! = .false. preconditioner is not update
  integer, intent(inout) :: not_conv             ! convergency
  character(len=1) :: precond  ! type of preconditioner: ' ', 'D', 'L', 'J'
  real :: rezid0
  integer:: iout, iter, i
  !real :: size(-1:3)
  integer:: restart = 30 !30    ! 50  ! GMRES restarted after 'restart' iterations  !45
  integer:: nloops = 50 !25 !30  !5     ! 10   ! maximal number of restarted cycles !100	  !40
  real :: t0, t1, t2

  if(state%space%adapt%adapt_method == 'ALG') then

     if (stop_crit == 'L' .or. stop_crit == 'G') then

        if (state%time%iter == 0) then
           restart = nu
           nloops = 1
           state%linSolver%tol = 1E-15 !1E-25
           print*, '# Initialization of "Restart", "nloops", and "state%linSolver%tol" for GMRES.'
        elseif (state%time%iter_loc == 1) then   ! (state%time%iter /= 1 .and. state%time%iter_loc == 1) then
           restart = nu
           nloops = 1
           state%linSolver%tol = 1E-15 !1E-25
           print*, '# Initialization of "Restart", "nloops", &
                and "state%linSolver%tol" for GMRES in the beginning of new adaptation.'
        endif

        if (state%time%iter_SC > -1) then  ! "state%time%iter_SC =" index of a time step when Stopping Criteria based on AEE are satisfied
           if (restart == nu) then
              print*, '# "Restart" and "nloops" for GMRES is set to a higher value due to a convergence purpose.'
              restart = nu*2
              nloops = 30
           endif
        endif

     elseif (stop_crit == 'N') then

        if (state%time%iter == 0) then
           restart = nu
           nloops = 1
           state%linSolver%tol = 1E-15 !1E-25
           print*, '# Initialization of "Restart", "nloops", and "state%linSolver%tol" for GMRES.'
        endif

        !restart = 50 !30
        !nloops = 20 !30
        !state%linSolver%tol = 1E-11

     else
        Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
     endif


     !print*, 'state%err(SSL8), state%conv_rez =', state%err(SSL8), state%conv_rez
  endif


  if (state%space%adapt%adapt_method == 'ALG2') then

     if (state%first_GMRES_conv) then
        !if (stop_crit /= 'N') then
        restart = 30 !nu
        nloops = 30 !1
        state%linSolver%tol = 1E-15
        !elseif (stop_crit == 'N') then
        !   if (state%time%iter == 0) then
        !     restart = nu
        !     nloops = 1
        !     state%linSolver%tol = 1E-15 !1E-25
        !     print*, '# Initialization of "Restart", "nloops", and "state%linSolver%tol" for GMRES.'
        !   endif
        !endif

        !if (state%time%iter_SC == 2) then
        !   restart = nu*2
        !   nloops = 30
        !endif

     else
        !if (stop_crit == 'L' .or. stop_crit == 'G') then

        if (state%time%iter == 0) then
           restart = nu
           nloops = 1
           state%linSolver%tol = 1E-15 !1E-25
           print*, '# Initialization of "Restart", "nloops", and "state%linSolver%tol" for GMRES.'
        elseif (state%time%iter_loc == 1) then   ! (state%time%iter /= 1 .and. state%time%iter_loc == 1) then
           restart = nu
           nloops = 1
           state%linSolver%tol = 1E-15 !1E-25
           print*, '# Initialization of "Restart", &
                "nloops", and "state%linSolver%tol" for GMRES in the beginning of new adaptation.'
        endif

        !elseif (stop_crit == 'N') then
        !   print*, 'Error! state%first_GMRES_conv == false together with stop_crit == N cannot happen!!!'
        !   stop

        !else
        !  print*, 'Error! Only three possibilities of stopping criterion: &
        !   & L (local), G (global), and N (classical stop. crit.).'

        !endif ! 'L' 'G'

     endif ! state%first_GMRES_conv

  endif  ! ALG2


  !state%linSolver%tol = 0.25D+00

  !!call EvalWeightTDresid(b, x, nsize, res1, res2, 1 )
  !if (.not. state%linSolver%tol_fixed) then
  !   if(state%time%iter == 0) then
  !      state%linSolver%tol = 0.1
  !   else
  !      !state%linSolver%tol = 0.1*(res1 - state%linSolver%residuum)/ state%linSolver%consistency
  !       !state%linSolver%tol = max(1E-4, min(res2/ state%linSolver%consistency, 0.75) )
  !
  !      state%linSolver%tol = max(1E-4, min(0.5*(res1/ state%linSolver%consistency)**1.5, 0.75) )
  !      !state%linSolver%tol = max(min(1D-04, state%err(SSL8)/1000.), 1D-15)
  !   endif
  !endif

  if(state%linSolver%tol <= 0.) then
     print*,'Zero tolerance state%linSolver%tol '
     stop
  endif
  !print*,'@@@@',state%linSolver%tol_fixed,state%linSolver%tol

  !!state%linSolver%consistency = res1

!!! choice of the method  Given by state%linSolver%name from .ini file
  !precond = ' '   ! GMRES + NO  preconditioning
  !precond = 'D'   ! GMRES + block diagonal preconditioning
  !precond = 'L'    ! GMRES + block ILU(0) preconditioning
  !precond = 'T'   ! Taylor expansis solution (does not work)


  if(state%space%adapt%adapt_method == 'ALG') precond = 'L' !' '

  iout = 0    ! no printing
  !iout = 1    ! GMRES printing

  mx_eta = eta
  !modified gmres ....
  !!if(precond .eq. ' ') then

  !print*,state%linSolver%name
  !PAUSE

  if(state%linSolver%name == "GMRES") then
     !restart = 50
     !nloops = 150
     !print*,'@@@@@ tol',state%linSolver%tol
     !call gmres(nsize, x, b, restart*nloops, state%linSolver%tol,  &
     call gmres(nsize, x, b, restart*nloops, 1.,  &
          bMVprod, bMVnull, restart,  state%linSolver%tol, iout, iter, rezid, &
          not_conv)
     !print*,'@@@@@ rez',rezid**0.5

     !elseif(precond .eq. 'D') then
  elseif(state%linSolver%name == "GMRES_D") then

     !do i=1,grid%nelem
     ! output for Zdenek Strakos, comparison of diagonal and off-digonal blocks
     !size(:) = 0.
     !do j=-1, grid%elem(i)%flen
     !   if(j == -1) then
     !      do k=1,grid%elem(i)%dof
     !         size(j) = size(j) +  &
     !              dot_product(grid%elem(i)%Mass%Mb(k,:), &
     !              grid%elem(i)%Mass%Mb(k,:) ) *ndim
     !      enddo
     !   elseif(j == 0) then
     !      do k=1,grid%elem(i)%dof*ndim
     !         size(j) = size(j) +  &
     !              dot_product(grid%elem(i)%block(0)%Mb(k,:), &
     !              grid%elem(i)%block(0)%Mb(k,:) )
     !      enddo
     !   elseif( grid%elem(i)%face(neigh,j) > 0 ) then
     !      do k=1,grid%elem(i)%dof*ndim
     !         size(j) = size(j) +  &
     !              dot_product(grid%elem(i)%block(j)%Mb(k,:), &
     !           grid%elem(i)%block(j)%Mb(k,:) )
     !      enddo
     !   endif
     !enddo
     !
     !size(-1:3) = size(-1:3)**0.5/ (grid%elem(i)%dof*ndim)
     !
     !write(100+state%time%iter,*) i, size(-1:3),grid%elem(i)%area
     !enddo

     if(precond_update) call ComputeBlockDiagPrecond( )

     call gmres(nsize, x, b, restart*nloops, state%linSolver%tol,  &
          bMVprod, bMVdiagprod, restart,  state%linSolver%tol, iout, iter, rezid, &
          not_conv)

     !elseif(precond .eq. 'L') then
  elseif(state%linSolver%name == "GMRES_ILU") then

     !call cpu_time(t0)
     !print*, 'HERE?'
     if(precond_update)  call ComputeBlockILUPrecond( )

     !call cpu_time(t1)

     !print*,'@@@@@ tol',state%linSolver%tol
     !call gmres(nsize, x, b, restart*nloops, 1.,  &   ! TEST FOR ALGEB
     !if (state%space%adapt%adapt_method == 'ALG2' .or. state%space%adapt%adapt_method == 'ALG') print*, 'restart, nloops = ', restart, nloops
     !print*, 'eta =', eta

     !print*,'gmres eta:   eta                      rezid                     state%MTol                         iter'
     !print*,'gmres eta:',0.0,rezid,state%MTol,iter,state%linSolver%tol_fixed

     call cpu_time(t1)
     call gmres(nsize, x, b, restart*nloops, state%linSolver%tol,  &
          bMVprod, bMViLUprod, restart,  state%linSolver%tol, iout, iter, rezid, &
          not_conv)
     call cpu_time(t2)

     !print*,'gmres eta:',t2-t1,rezid,state%MTol,iter,state%linSolver%tol_fixed

     !print*, 'rezid =', rezid
     !pause
     !print*,'@@@@@ rez',rezid**0.5
     !print*,'@@@@@ rez',rezid
     !print*,'@@@@@ state%linSolver%residuum',state%linSolver%residuum
     !call cpu_time(t2)

     !write(74,*) state%time%iter+1, t1-t0, t2-t1, t2-t0

     !!elseif(precond .eq. 'T') then
     !!
     !!call ComputeBlockDiagPrecond( )

     !!call TaylorSolution(nsize, x, b, state%linSolver%tol, iter, rezid, not_conv)

     !    elseif(state%linSolver%name == "ILU") then
     !
     !       call cpu_time(t0)
     !
     !       !if(precond_update)  call ComputeBlockILUPrecond( )
     !       if(precond_update)  call  ComputeBlockDiagPrecond()
     !       call cpu_time(t1)
     !
     !       !print*,'@@@@@ tol',state%linSolver%tol
     !       !call gmres(nsize, x, b, restart*nloops, 1.,  &   ! TEST FOR ALGEB
     !       !if (state%space%adapt%adapt_method == 'ALG2' .or. state%space%adapt%adapt_method == 'ALG') print*, 'restart, nloops = ', restart, nloops
     !       !print*, 'eta =', eta
     !
     !       !print*,'gmres eta:   eta                      rezid                     state%MTol                         iter'
     !       !print*,'gmres eta:',0.0,rezid,state%MTol,iter,state%linSolver%tol_fixed
     !
     !       call cpu_time(t1)
     !       !call gmres(nsize, x, b, restart*nloops, state%linSolver%tol,  &
     !       !     bMVprod, bMViLUprod, restart,  state%linSolver%tol, iout, iter, rezid, &
     !       !     not_conv)
     !       rezid0 = L2Res(x,b,nsize,bMVprod)**2
     !
     !       do i=1,1
     !          !call bMViLUprod(x,b,nsize)
     !          call bMVdiagprod(x,b,nsize)
     !       enddo
     !       rezid = L2Res(x,b,nsize,bMVprod)**2
     !
     !       call cpu_time(t2)
     !
     !       write(*,'(a6,l2,8es12.4)') '####', precond_update, t2-t1, t1 - t0, rezid0, rezid, rezid/ rezid0
     !
     !
     !    elseif ( state%linSolver%name == 'UMFPACK' ) then
     !        print*,'UMFPACK solver -- not fully tested yet'
     !        call SolveViaUMFPACK(nsize,x,b)
     !        iter = 1
     !        rezid = L2Res(x,b,nsize,bMVprod)**2
     !        !
     !    elseif ( state%linSolver%name == 'AGMG' ) then
     !        print*,'AGMG solver -- are parameters correct?'
     !        call SolveViaAGMG(nsize,x,b)
     !        iter = 1
     !        rezid = L2Res(x,b,nsize,bMVprod)**2
     !        !
     !   elseif ( state%linSolver%name == 'JACOBI' ) then
     !        !print*,'bJacobi iterations'
     !        rezid = 10*state%linSolver%tol**2
     !        iter = 0
     !        do  while ( rezid > state%MTol**2 )
     !            !call cpu_time(t1)
     !            call bJacobi(x,b,nsize)
     !            iter = iter+1
     !            rezid = L2Res(x,b,nsize,bMVprod)**2
     !            !call cpu_time(t2)
     !            !print*,'1x Jacobi iteration eta:',t2-t1
     !        end do
     !        !
     !    elseif ( state%linSolver%name == 'GS' ) then
     !        !print*,'bGS iterations'
     !        rezid = 10*state%linSolver%tol**2
     !        iter = 0
     !        do  while ( rezid > state%MTol**2 )
     !            !call cpu_time(t1)
     !            call bGS(x,b,nsize)
     !            iter = iter+1
     !            rezid = L2Res(x,b,nsize,bMVprod)**2
     !            !call cpu_time(t2)
     !            !print*,'1x G-S iteration eta:',t2-t1
     !        end do
     !
     !    elseif ( state%linSolver%name == 'MG1JACGS' ) then
     !
     !        MGRun%presmooth = 1
     !        MGRun%postsmooth = 1
     !
     !        !print*,'bGS iterations'
     !        rezid = 10*state%linSolver%tol**2
     !        iter = 0
     !        do  while ( rezid > state%MTol**2 )
     !            !call cpu_time(t1)
     !            call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'bGS',rezid,bJacobi)
     !            iter = iter+1
     !            rezid = L2Res(x,b,nsize,bMVprod)**2
     !            !call cpu_time(t2)
     !            !print*,'1x G-S iteration eta:',t2-t1
     !        end do
     !        !
     !    elseif ( state%linSolver%name == 'MG2JACGS' ) then
     !
     !        MGRun%presmooth = 2
     !        MGRun%postsmooth = 2
     !
     !        !print*,'bGS iterations'
     !        rezid = 10*state%linSolver%tol**2
     !        iter = 0
     !        do  while ( rezid > state%MTol**2 )
     !            !call cpu_time(t1)
     !            call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'bGS',rezid,bJacobi)
     !            iter = iter+1
     !            rezid = L2Res(x,b,nsize,bMVprod)**2
     !            !call cpu_time(t2)
     !            !print*,'1x G-S iteration eta:',t2-t1
     !        end do
     !        !
     !    elseif ( state%linSolver%name == 'MG3JACGS' ) then
     !
     !        MGRun%presmooth = 3
     !        MGRun%postsmooth = 3
     !
     !        !print*,'bGS iterations'
     !        rezid = 10*state%linSolver%tol**2
     !        iter = 0
     !        do  while ( rezid > state%MTol**2 )
     !            !call cpu_time(t1)
     !            call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'bGS',rezid,bJacobi)
     !            iter = iter+1
     !            rezid = L2Res(x,b,nsize,bMVprod)**2
     !            !call cpu_time(t2)
     !            !print*,'1x G-S iteration eta:',t2-t1
     !        end do
     !        !
     !    elseif ( state%linSolver%name == 'MG1JACJAC' ) then
     !
     !        MGRun%presmooth = 1
     !        MGRun%postsmooth = 1
     !
     !        !print*,'bGS iterations'
     !        rezid = 10*state%linSolver%tol**2
     !        iter = 0
     !        do  while ( rezid > state%MTol**2 )
     !            !call cpu_time(t1)
     !            call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'bJacobi',rezid,bJacobi)
     !            iter = iter+1
     !            rezid = L2Res(x,b,nsize,bMVprod)**2
     !            !call cpu_time(t2)
     !            !print*,'1x G-S iteration eta:',t2-t1
     !        end do
     !        !
     !    elseif ( state%linSolver%name == 'MG2JACJAC' ) then
     !
     !        MGRun%presmooth = 2
     !        MGRun%postsmooth = 2
     !
     !        !print*,'bGS iterations'
     !        rezid = 10*state%linSolver%tol**2
     !        iter = 0
     !        do  while ( rezid > state%MTol**2 )
     !            !call cpu_time(t1)
     !            call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'bJacobi',rezid,bJacobi)
     !            iter = iter+1
     !            rezid = L2Res(x,b,nsize,bMVprod)**2
     !            !call cpu_time(t2)
     !            !print*,'1x G-S iteration eta:',t2-t1
     !        end do
     !        !
     !    elseif ( state%linSolver%name == 'MG3JACJAC' ) then
     !
     !        MGRun%presmooth = 3
     !        MGRun%postsmooth = 3
     !
     !        !print*,'bGS iterations'
     !        rezid = 10*state%linSolver%tol**2
     !        iter = 0
     !        do  while ( rezid > state%MTol**2 )
     !
     !            !call cpu_time(t1)
     !            call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'bJacobi',rezid,bJacobi)
     !            iter = iter+1
     !            rezid = L2Res(x,b,nsize,bMVprod)**2
     !            !call cpu_time(t2)
     !            !print*,'1x G-S iteration eta:',t2-t1, rezid, state%MTol**2
     !        end do
     !        !
     !    elseif ( state%linSolver%name == 'MG1GSGS' ) then
     !
     !        MGRun%presmooth = 1
     !        MGRun%postsmooth = 1
     !
     !        !print*,'bGS iterations'
     !        rezid = 10*state%linSolver%tol**2
     !        iter = 0
     !        do  while ( rezid > state%MTol**2 )
     !            !call cpu_time(t1)
     !            call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'bGS',rezid,bGS)
     !            iter = iter+1
     !            rezid = L2Res(x,b,nsize,bMVprod)**2
     !            !call cpu_time(t2)
     !            !print*,'1x G-S iteration eta:',t2-t1
     !        end do
     !        !
     !    elseif ( state%linSolver%name == 'MG2GSGS' ) then
     !
     !        MGRun%presmooth = 2
     !        MGRun%postsmooth = 2
     !
     !        !print*,'bGS iterations'
     !        rezid = 10*state%linSolver%tol**2
     !        iter = 0
     !        do  while ( rezid > state%MTol**2 )
     !            !call cpu_time(t1)
     !            call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'bGS',rezid,bGS)
     !            iter = iter+1
     !            rezid = L2Res(x,b,nsize,bMVprod)**2
     !            !call cpu_time(t2)
     !            !print*,'1x G-S iteration eta:',t2-t1
     !        end do
     !        !
     !    elseif ( state%linSolver%name == 'MG3GSGS' ) then
     !
     !        MGRun%presmooth = 4
     !        MGRun%postsmooth = 4
     !
     !        !print*,'bGS iterations'
     !        rezid = 10*state%linSolver%tol**2
     !        rezid0 = L2Res(x,b,nsize,bMVprod)**2
     !        iter = 0
     !        do  while ( rezid/rezid0 > state%MTol**2 )
     !        ! VD
     !        !do  while ( rezid > state%MTol**2 )
     !           ! VD
     !           call cpu_time(t1)
     !           call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'bGS',rezid,bGS)
     !           iter = iter+1
     !           rezid = L2Res(x,b,nsize,bMVprod)**2
     !           call cpu_time(t2)
     !           write(*,'(a25,6es14.6)') '1x G-S iteration eta:',t2-t1, rezid, rezid0, rezid/rezid0, state%MTol**2
     !        end do
     !        !
     !    elseif( state%linSolver%name == 'MG_Jacobi1' ) then
     !        ! MG cycle with "1" block Jacobi step used as an exact solution procedure
     !        !
     !        rezid = 10*state%linSolver%tol**2
     !        iter = 0
     !        print*,'1x V-cyle eta:   eta                      rezid                     state%MTol                         iter'
     !        !
     !        do  while( rezid > state%MTol**2 )
     !            call cpu_time(t1)
     !            !
     !            !call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'bJacobi1',rezid)
     !            call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'GMRES',rezid,bJacobi)
     !            iter = iter+1
     !            rezid = L2Res(x,b,nsize,bMVprod)**2
     !            !
     !            call cpu_time(t2)
     !            print*,'1x V-cyle eta:',t2-t1,rezid,state%MTol,iter
     !        end do
     !        !
     !    elseif( state%linSolver%name == 'MG_Jacobi1' ) then
     !        ! MG cycle with "1" block Jacobi step used as an exact solution procedure
     !        !
     !        rezid = 10*state%linSolver%tol**2
     !        iter = 0
     !        print*,'1x V-cyle eta:   eta                      rezid                     state%MTol                         iter'
     !        !
     !        do  while( rezid > state%MTol**2 )
     !            call cpu_time(t1)
     !            !
     !            !call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'bJacobi1',rezid)
     !            call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'GMRES',rezid,bJacobi)
     !            iter = iter+1
     !            rezid = L2Res(x,b,nsize,bMVprod)**2
     !            !
     !            call cpu_time(t2)
     !            print*,'1x V-cyle eta:',t2-t1,rezid,state%MTol,iter
     !        end do
     !        !
     !    elseif( state%linSolver%name == 'MG_JacobiX' ) then
     !        ! MG cycle with "10" block Jacobi steps used as an exact solution procedure
     !        !
     !        rezid = 10*state%linSolver%tol**2
     !        iter = 0
     !        print*,'1x V-cyle eta:   eta                      rezid                     state%MTol                         iter'
     !        !
     !        do  while( rezid > state%MTol**2 )
     !            call cpu_time(t1)
     !            !
     !            call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'GMRES',rezid,bJacobi)
     !            iter = iter+1
     !            rezid = L2Res(x,b,nsize,bMVprod)**2
     !            !
     !            call cpu_time(t2)
     !            print*,'1x V-cyle eta:',t2-t1,rezid,state%MTol,iter
     !        end do
     !        !
     !    elseif( state%linSolver%name == 'MG_ILU_DO-NOT-USE' ) then
     !        ! DO NOTHING
     !    elseif( state%linSolver%name == 'MGxGMRES' ) then
     !        ! MG cycle with "10" block Jacobi steps used as an exact solution procedure
     !        !
     !        MGRun%presmooth = 2
     !        MGRun%postsmooth = 2
     !
     !
     !        rezid0 = L2Res(x,b,nsize,bMVprod)**2
     !        !rezid = 10*state%linSolver%tol**2
     !        rezid = rezid0 + 1
     !
     !        iter = 0
     !        print*,'@@@@ 1x V-cyle eta:   eta                      rezid                     state%MTol                         iter'
     !        print*,'@@@@',  MGRun%presmooth,  MGRun%postsmooth, rezid, rezid0
     !        !
     !        !do  while( rezid > state%MTol**2 )
     !        do  while ( rezid /rezid0 > state%MTol**2 )
     !           call cpu_time(t1)
     !           !
     !           !call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'bGS',rezid,bGS)
     !           call pMGLinSolverRecu4LINSOL(nsize,eta,b,x,maxval( grid%elem(:)%deg ),'GMRES',rezid,bGS)
     !            iter = iter+1
     !            rezid = L2Res(x,b,nsize,bMVprod)**2
     !            !
     !            call cpu_time(t2)
     !            write(*,'(a25,6es14.6)') '@@@@ 1x V-cyle eta:',t2-t1,rezid,rezid0, rezid/rezid0, state%MTol**2
     !            print*,'##############################################################'
     !        end do
     !        !
  else
     print*,'Unknown linear solver', state%linSolver%name,' in "lin_solvers.f90" -- STOP'
     STOP
  endif

  state%linSolver%iter_tot = state%linSolver%iter_tot + iter
  tot_iter = tot_iter + iter

  rezid = rezid**0.5

  !!call  EvalWeightTDresid(b, x, nsize, state%linSolver%residuum, res2, 0 )
  !call  EvalWeightTDresid(b, x, nsize, res1, res2, 0 )

  !!write(*,'(a6,i5,5es12.4)') 'GMRES:',state%linSolver%iter, state%linSolver%residuum, rezid, res1, &
  !!     res1 / VectorPrecondNorm(b(:) )

  !mx_eta = 0.

end subroutine SolveBlockLinearProblem

!> solution of large block linear algebraic problem
!> \f$ (A + \eta\,M)x_1 + \iota_2\,My = b_1, \iota_1\,Mx_1 + (A+\eta\,M)x_2 = b_2\f$
!>
!> \f$ A\f$ is a sparse block matrix given by elem%block(0:len),
!> \f$ x\f$ is vecor with initial guess of the solution (in) and the solution (out),
!> \f$ b\f$ is the righ-hand side
subroutine SolveBlockLinearDoubleProblem(nsize, eta, iota, b, x, rezid, tot_iter, not_conv)
  use matrix_oper_int, mx_eta => eta, mx_iota => iota
  integer :: nsize                               ! size of the algebraic problem
  real, intent(in):: eta, iota(2)
  real, dimension(1:nsize,1:nbDim), intent(inout):: x    ! solution
  real, dimension(1:nsize,1:nbDim), intent(in)   :: b    ! RHS
  real :: rezid                                  ! reziduum
  integer :: tot_iter                            ! number of iterations
  integer :: not_conv                            ! convergency
  character(len=1) :: precond  ! type of preconditioner: ' ', 'D', 'L', 'J'
  integer:: iout
  integer:: restart = 40    ! GMRES restarted after 'restart' iterations  !45
  integer:: nloops = 50    ! maximal number of restarted cycles !100	  !40
  !!real   :: tol  !=1D-05  !1D-03       ! tolerance for GMRES solver  tol=1D-08

  !state%linSolver%tol = max(min(1D-04, state%err(SSL8)/1000.), 1D-15)
  state%linSolver%tol = 1D-6
  if(state%modelName == 'scalar' .or.state%modelName == '2eqs') state%linSolver%tol = 1D-10

  if (.not. state%linSolver%tol_fixed) then
     state%linSolver%tol = max(min(1D-04, state%err(SSL8)/1000.), 1D-15)
  endif
  !precond = ' '   ! NO  preconditioning
  precond = 'D'   ! block diagonal preconditioning

  ! off-diagonal terms are not allocated !!!
  !precond = 'L'   ! block ILU(0) preconditioning

  iout=0    ! no printing

  mx_eta = eta
  mx_iota = iota
  !modified gmres ....
  if(precond .eq. ' ') then
     call gmres(2*nsize, x, b, restart*nloops, 2*state%linSolver%tol, &
          bMVprod2, bMVnull, restart,  2*state%linSolver%tol/100, iout, tot_iter, rezid, &
          not_conv)

  elseif(precond .eq. 'D') then
     call ComputeBlockDiagPrecond()

     call gmres(2*nsize, x, b, restart*nloops, 2*state%linSolver%tol, &
          bMVprod2, bMVdiagprod2, restart,  2*state%linSolver%tol/100, iout, tot_iter, rezid, &
          not_conv)

  elseif(precond .eq. 'L') then
     call ComputeBlockILUPrecond()

     call gmres(2*nsize, x, b, restart*nloops, 2*state%linSolver%tol, &
          bMVprod2, bMViLUprod2, restart,  2*state%linSolver%tol/100, iout, tot_iter, rezid, &
          not_conv)
  else
     print*,'Any other preconditioner not implemented in "euler.f90"'
     stop
  endif

  mx_eta = 0.
  mx_iota = 0.
end subroutine SolveBlockLinearDoubleProblem

!> solution of a linear problem of the form \f$ AX + MXK = B\f$, where
!>
!> \f$ A\f$ is a sparse block matrix given by elem%block(0:len),
!> \f$ M\f$ is the block diagonal mass matrix
!> \f$ K\f$ is a small dense matrix
!> \f$ X\f$ is a matrix with initial guess of the solution (in)
!> and the solution (out),
!> \f$ B\f$ is the right-hand side
subroutine SolveMixedSylvesterProblem(X,B,K,rezid,tot_iter,not_conv)
real,intent(in):: K(:,:), B(:,:)
real,intent(inout):: X(:,:)
real,intent(out):: rezid
integer,intent(out):: tot_iter,not_conv
real:: Z(size(K,1),size(K,1)),T(size(K,1),size(K,2))
real:: C(size(B,1),size(B,2)),d(size(B,1))
!real:: C1(size(B,1),size(B,2))
real:: rezi,rconde,rcondv
real,dimension(size(K,1)):: WR,WI
integer:: nk,info,i,toti, idummy
external:: dgeesx,xerbla ! sgeesx,
real:: work(size(K,1)*10)
integer:: iwork(1)
logical:: bwork(1)

nk = size(K,1)
if (nk /= size(K,2)) then
   stop 'dimension mismatch in SolveMixedSylvesterProblem'
end if

if (any(isnan(K))) then
   stop 'SolveMixedSylvesterProblem: NaNs in the K matrix'
endif
! form the Schur decomposition
T = K
!call la_geesx(T,WR,WI,Z,info = info)
if (kind(T) == 8) then
   call dgeesx('V','N',xerbla,'N',size(T,1),T,size(T,1),idummy,WR,WI,Z,size(Z,1),&
        rconde,rcondv,work,size(work),iwork,size(iwork),bwork,info)
   !else if (kind(T) == 4) then
   !   call sgeesx('V','N',xerbla,'N',size(T,1),T,size(T,1),idummy,WR,WI,Z,size(Z,1),&
   !  rconde,rcondv,work,size(work),iwork,size(iwork),bwork,info)
else
   stop 'unsupported precision for LAPACK'
endif

! transform to Schur basis
C = matmul(B,Z)
X = matmul(X,Z)

i = 1
rezid = 0
tot_iter = 0
main:do while (i <= nk)
   if (WI(i) == 0.) then
      ! real (single) subproblem
      if (i > 1) then
         call bMVmassprod(d,matmul(X(:,1:i-1),T(1:i-1,i)),size(d))
         C(:,i) = C(:,i) - d
      end if
      ! solve
      call SolveBlockLinearProblem(size(X,1),T(i,i),C(:,i),X(:,i),&
           rezi,toti, .true.,  not_conv)
      rezid = rezid + rezi
      tot_iter = tot_iter + toti
      if (not_conv /= 0) exit main
      i = i + 1
   else
      ! complex (double) subproblem
      if (i > 1) then
         call bMVmassprod(d,matmul(X(:,1:i-1),T(1:i-1,i)),size(d))
         C(:,i) = C(:,i) - d
         call bMVmassprod(d,matmul(X(:,1:i-1),T(1:i-1,i+1)),size(d))
         C(:,i+1) = C(:,i+1) - d
      end if
      ! solve
      call SolveBlockLinearDoubleProblem(size(X,1),T(i,i),[T(i+1,i),T(i,i+1)],&
           C(:,i:i+1),X(:,i:i+1),rezi,toti,not_conv)
      rezid = rezid + rezi
      tot_iter = tot_iter + toti
      if (not_conv /= 0) exit main
      i = i + 2
   end if
end do main

! transfer back to normal basis
X = matmul(X,transpose(Z))

! FIXME: WORKING (?)
! C = 0
! C1 = 0
! do i = 1,size(X,2)
!   call bMVprod(C(:,i),X(:,i),size(X,1))
!   call bMVmassprod(C1(:,i),X(:,i),size(X,1))
! end do
! C = C + matmul(C1,K) - B
! print *,'resid = ',sum(C**2) / sum(B**2)

end subroutine SolveMixedSylvesterProblem

!> form matrices of the problem
!>\f$\int_0^\tau (u(t)',v(t)) + (Au(t),v(t))\ dt + (u(0),v(0)) =
!> \int_0^\tau (b(t),v(t)) + (u_0,v(0))\f$
!> on input, B is b evaluated at quadrature points
!> on return, B is updated to be the rhs of Sylvester problem.
!> K is set so that the values at the quadrature points can be found by
!> solving A*X + M*X*K = B.
!> qp are quadrature points
!> qw are quadrature weights (sum up to 1)
!> tau is the interval length.
!> the quadrature rule should be exact for polynomials up to degree 2*(n-1)
subroutine PrepareTimeElementLinearProblem(tau,qw,qp,B,K,X0)
use helpers
real,intent(in):: tau,qw(:),qp(:),X0(:)
real,intent(inout):: B(:,:)
real,intent(out):: K(:,:)
real,dimension(size(qw),size(qw)):: DD,D
real:: MX0(size(X0))
real:: Y(size(qw))
integer:: n,i,j

n = size(qp)

! form a matrix of base polynomials in the divided differences form
DD = 0
forall(i=1:n) DD(i,i) = 1

do i=2,n
   do j=n,i,-1
      DD(:,j) = (DD(:,j) - DD(:,j-1)) / (qp(j) - qp(j+1-i))
   end do
end do

! form a matrix of derivatives - D(i,j) is phi_i'(x_j).
do j=1,n
   Y = DD(:,n)
   D(:,j) = 0
   do i = n-1,1,-1
      D(:,j) = Y + (qp(j) - qp(i)) * D(:,j)
      Y = DD(:,i) + (qp(j) - qp(i)) * Y
   end do
end do

! call saveoct('D',D,'N') (OK)

! calc values at zero
Y = DD(:,n)
do i = n-1,1,-1
   Y = DD(:,i) + (0 - qp(i)) * Y
end do

! call saveoct('Y',reshape(Y,[n,1]),'N') (OK)

! calculate matrix of the derivative form (u',v) + (u(0),v(0))
do i = 1,n
   K(:,i) = (D(:,i) + Y * Y(i) / qw(i)) / tau
end do

! call saveoct('K',K,'N') (OK)

! for M*X0
call bMVmassprod(MX0,X0,size(X0))

! update B
do i=1,n
   B(:,i) = B(:,i) + MX0 * Y(i) / (tau * qw(i))
end do

end subroutine PrepareTimeElementLinearProblem


!> solution of large block linear algebraic problem \f$ (A + \eta\,M)x= b\f$
!>
!> \f$ A\f$ is a sparse block matrix given by elem%block(0:len),
!> \f$ x\f$ is vecor with initial guess of the solution (in) and the solution (out),
!> \f$ b\f$ is the righ-hand side
subroutine SolveBlockLinearSTDGMProblem( nsize, eta, b, x, rezid, tot_iter, &
   precond_update, not_conv )
use matrix_oper_int, mx_eta => eta
integer :: nsize                               ! size of the algebraic problem
real, intent(in):: eta
real, dimension(1:nsize), intent(inout):: x    ! solution
real, dimension(1:nsize), intent(inout):: b    ! RHS
real, intent(inout) :: rezid                   ! reziduum
integer, intent(inout) :: tot_iter             ! number of iterations
logical, intent(in) :: precond_update          ! = .false. preconditioner is not update
integer, intent(inout) :: not_conv             ! convergency
character(len=20) :: precond  ! type of preconditioner: ' ', 'D', 'L', 'J'
integer:: iout, iter
real :: start, finish
integer:: restart = 30 !30    ! 50  ! GMRES restarted after 'restart' iterations  !45
integer:: nloops = 5  !5     ! 10   ! maximal number of restarted cycles !100	  !40

if(state%linSolver%tol <= 0.) then
   print*,'Zero tolerance state%linSolver%tol '
   stop
endif

!!state%linSolver%consistency = res1

!!! choice of the method
!precond = ' '   ! GMRES + NO  preconditioning
!precond = 'D'   ! GMRES + block diagonal preconditioning
! if (state%linSolver%name .ne. '') then

precond = state%linSolver%name

! else
! 	precond = 'GMRES_ILU'    ! GMRES + block ILU(0) preconditioning
! endif
!print*, precond

iout = 0    ! no printing
!iout = 1    ! GMRES printing

mx_eta = eta
!modified gmres ....
if(precond .eq. 'GMRES') then
   !restart = 50
   !nloops = 150
   !print*,'@@@@@ tol',state%linSolver%tol
   !call gmres(nsize, x, b, restart*nloops, state%linSolver%tol,  &
   !       write(debug,*) 'calling gmres with no preconditioning'
   !
   !       write(debug, *) 'nsize: ' , nsize
   !!       write(debug, *) 'X : ', x
   !!       write(debug,*) 'B : ' , b
   !       write(debug, *) 'TOL: ' , state%linSolver%tol


   !call cpu_time( start)
   call gmres(nsize, x, b, restart*nloops, 1.,  &
        bMVprodST, bMVnull, restart,  state%linSolver%tol, iout, iter, rezid, &
        not_conv)

   !      write(debug, *) '-------'
   !      !write(debug, *) 'X : ', x
   !      write(debug,*) 'restart, nloops: ', restart, nloops
   !      write(debug,*) 'state%linSolver%tol' , state%linSolver%tol
   !      write(debug,*) 'iout : ' , iout
   !      write(debug,*) 'iter, rezid:' , iter, rezid
   !      write(debug,*) 'not_conv : ' , not_conv
   !      write(debug,*) '||x|| = ' , VectorNorm(x)
   !      close(debug)
   !      stop 'STOPPING in lin_solvers'

   !call cpu_time( finish)
   !write(63,*) 'TIME FOR GMRES:', finish-start
   !write(63,*) '________________________________'

   !print*,'@@@@@ rez',rezid**0.5

elseif(precond .eq. 'GMRES_D') then

   if(precond_update) call ComputeBlockDiagPrecondST( )
   !write(63,*) 'calling GMRES with BD preconditioning'

   !call cpu_time( start)
   call gmres(nsize, x, b, restart*nloops, state%linSolver%tol,  &
        bMVprodST, bMVdiagprodST, restart,  state%linSolver%tol, iout, iter, rezid, &
        not_conv)

   !call cpu_time( finish)
   !write(63,*) 'TIME FOR GMRES:', finish-start
   !write(63,*) '________________________________'

elseif(precond .eq. 'GMRES_ILU') then
   !call cpu_time(t0)
!   write(*,*) 'calling gmres with ILU preconditioning'

   if(precond_update)  call ComputeBlockILUPrecondST( )

   !call cpu_time( start)
   !print*, 'Gmres!'
   call gmres(nsize, x, b, restart*nloops, state%linSolver%tol,  &
        bMVprodST, bMViLUprodST, restart,  state%linSolver%tol, iout, iter, rezid, &
        not_conv)
   !call cpu_time( finish)
   !print*, 'after GMRES'
   !write(63,*) 'TIME FOR GMRES:', finish-start
   !write(63,*) '________________________________'
   !print*, 'rezid =', rezid
   !pause
   !print*,'@@@@@ rez',rezid**0.5
   !print*,'@@@@@ rez',rezid
   !print*,'@@@@@ state%linSolver%residuum',state%linSolver%residuum
   !call cpu_time(t2)

   !write(74,*) state%time%iter+1, t1-t0, t2-t1, t2-t0

elseif(precond .eq. 'T') then
   print*, 'Taylor solution not implemented fo STDGM'
   stop


   call ComputeBlockDiagPrecond( )

   call TaylorSolution(nsize, x, b, state%linSolver%tol, iter, rezid, not_conv)

else
   print*,'Any other preconditioner not implemented in "euler.f90"', precond
   stop
endif

state%linSolver%iter_tot = state%linSolver%iter_tot + iter
tot_iter = tot_iter + iter

rezid = rezid**0.5

!    print*, 'END of SolveBlockLinearSTDGMProblem'

end subroutine SolveBlockLinearSTDGMProblem


function SteadyStateResidual(X,RHS) result(res)
  use matrix_oper_int
  real,intent(in):: X(:),RHS(:)
  real:: res,Y(size(X))
  call bMVprod(Y,X,size(X))
  res = VectorNorm(RHS - Y) / VectorNorm(RHS)
end function SteadyStateResidual



end module lin_solvers
