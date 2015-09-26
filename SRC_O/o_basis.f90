!> setting of basis functions
module basis_mod
   use paramets
!  use mesh_oper
!  use main_data
!  use blocks_integ

  implicit none

   public :: PHIorthonormalNew
   public :: phi_vals



!  public:: SetBasisFunctions
!  public:: Lagr2Basis
!  public:: Lagr2BasisDiff
!  public:: Lagr2BasisDiffOne
!  public:: SetTrulePhi
!  public:: EvalTrulePhi
!  public:: Write_TrulePhi
!
!  public:: SetEdgeNodes
!  public:: SetEdgeNodesHG


contains



  !> evalaution of test functions \f$ \varphi_i,\ i=1,\dots, dof \f$ in
  !> the node \f$ xi(1:nbDim) \f$
  subroutine phi_vals(dof, xi, phi, Dphi)
    integer, intent(in) :: dof    ! number of test functions
    real, dimension(1:nbDim), intent(in)  :: xi   ! coodrinates of an integ node
    real, dimension(1:dof), intent(out) :: phi        ! value of test function
    real, dimension(1:dof,1:nbDim), intent(out), optional :: Dphi   ! value of derivatives of test function
    real :: Fx, Fy, Dfx, Dfy
    integer :: i, k, l

    l = 0
    do k = 0, MaxDegreeImplemented
       do i= 0, k
          l = l + 1
          if(l <= dof) then
             Fx = xi(1)**(k - i)
             Fy = xi(2)**(i)
             phi(l) = Fx * Fy

             if(present(Dphi)) then
                if(k - i == 0) then
                   Dfx = 0
                else
                   DFx = (k -i) * xi(1)**(k - i - 1)
                endif
                if( i == 0) then
                   DFy = 0.
                else
                   DFy = i * xi(2)**(i-1)
                endif

                Dphi(l, 1) = DFx *  Fy
                Dphi(l, 2) =  Fx * DFy
             endif

          endif
       enddo
    enddo
  end subroutine phi_vals

   !> evaluation of orthonormal test functions in node \f$ (x, y) \f$
  subroutine PHIorthonormalNew(Qdof, xi, len, dof, GScoeff, phi, Dphi)
    integer, intent(in) :: Qdof           ! number of integ nodes
    real, dimension(1:Qdof,1:nbDim),intent(in) :: xi ! coordinates of integ. node
    integer, intent(in) :: len           ! triangle/square
    integer, intent(in) :: dof           ! number of test functions
    real, dimension(1:dof, 1:dof), intent(in) :: GScoeff ! precomputed coefficients saved in state%space
    real, dimension(1:dof,1:Qdof), intent(out) :: phi        ! value of test function
    real, dimension(1:dof,1:nbDim, 1:Qdof), optional :: Dphi   ! derivatives of test function
    real, dimension(1:nbDim) :: xc
    integer :: i, k, l

    if(dof > ( MaxDegreeImplemented+1)*( MaxDegreeImplemented+2)/2 ) then
       print*,'Implementaion in ORTH_PHI only forl deg= MaxDegreeImplemented'
       stop
    else
       if(len == 3) then ! triangle
          xc(1:nbDim) = 1./3
       elseif(len == 4) then !square
          xc(1:nbDim) = 1./2
       else
          print*,'Only triangle or square in ORTH_PHI implemented'
          stop
       endif

       if(present(Dphi)) then
          do l=1,Qdof
             call phi_vals(dof, xi(l, 1:nbDim) - xc(1:nbDim), phi(1:dof, l), &
                  Dphi(1:dof, 1:nbDim, l) )
          enddo
       else
          do l=1,Qdof
             call phi_vals(dof, xi(l, 1:nbDim) - xc(1:nbDim), phi(1:dof, l) )
          enddo
       endif

       ! Gram-Schmidt
       do k=1,dof
          do i=1, k-1
             phi(k, 1:Qdof) = phi(k, 1:Qdof) - phi(i, 1:Qdof) * GScoeff( i, k)
          enddo
          phi(k, 1:Qdof) =  phi(k, 1:Qdof) / GScoeff( k, k)
       enddo

       if(present(Dphi)) then
          do k=1,dof
             do i=1, k-1
                Dphi(k, 1:nbDim, 1:Qdof) = Dphi(k, 1:nbDim, 1:Qdof) &
                     - Dphi(i, 1:nbDim, 1:Qdof) * GScoeff( i, k)
             enddo
             Dphi(k, 1:nbDim, 1:Qdof) =  Dphi(k, 1:nbDim, 1:Qdof) / GScoeff( k, k)
          enddo
       endif
    end if

  end subroutine PHIorthonormalNew

end module basis_mod

