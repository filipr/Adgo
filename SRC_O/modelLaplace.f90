!> definition of Laplace problem
module modelLaplace

  use main_data
  use f_mapping
  use mesh_oper
  use blocks_integ

  implicit none

  public:: Set_f_s_Laplace
  public:: Set_A_s_Laplace
  public:: Set_Ppm_Laplace
  public:: Set_R_s_Laplace
  public:: Set_K_sk_Laplace

contains

  !> compute Vectors \f$ f_s = {\bf f_s}({\bf w}) ,\ s=1,2\f$
  !> for scalar equation
  subroutine Set_f_s_Laplace(ndimL, nbDim, Qdof, w, f_s, x )
      integer, intent(in) :: Qdof, ndimL, nbDim
      real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
      real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
      real, dimension(1:Qdof,1 :nbDim), intent(in) :: x

       f_s(:,:,:) = 0.

  end subroutine Set_f_s_Laplace



  !> compute matrices \f$ A_s = \frac{D{\bf f_s}({\bf w})}{D{\bf w}},\ s=1,2\f$
  !> for scalar equation
  subroutine Set_A_s_Laplace(ndimL, nbDim, Qdof, w, A_s ,xi)
      integer, intent(in) :: Qdof, nbdim, ndimL
      real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
      real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: A_s
         ! matrices A_s in  -- " --
      real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi

      A_s(:, : , :, :) = 0.

  end subroutine Set_A_s_Laplace

  !> compute matrices
  !> \f$ P^{\pm} = \left(\frac{D({\bf f_1}({\bf w})n_1+{\bf f_2}({\bf w})n_2}{D{\bf w}}
  !>  \right)^{\pm}\f$
  !> for scalar equation
  subroutine Set_Ppm_Laplace( ndimL, nbDim, Qdof, w, n, xi, Ppm, one_over_area, elem)
      integer, intent(in) :: Qdof, ndimL, nbDim
      real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
      real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: Ppm
                                            ! matrices Ppm in  -- " --
      real, dimension(1:Qdof, 1:nbDim), intent(in) :: n   ! outer normal
      real, dimension(1:Qdof, 1:nbDim),intent(in) ::  xi                    ! node on the edge?
      real, intent(in), optional :: one_over_area
      type(element), intent(inout), optional :: elem


      Ppm(1:Qdof, 1:nbDim, 1:ndimL, 1:ndimL) = 0.


  end subroutine Set_Ppm_Laplace


  !> compute matrices ndim x ndim K_sk, s,k=1,2 for scalar euation
  !> in integ nodes, LAPLACE case
  subroutine Set_K_sk_Laplace(ndimL, nbDim,Qdof, w, Dw, Re_1, K_sk)
    integer, intent(in) :: ndimL, nbDim, Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,ndimL), intent(inout) :: K_sk
    integer :: i,j


    !initialization
    K_sk(1:Qdof,1:nbDim,1:nbDim,1:ndimL,1:ndimL) = 0.

    ! diagonal structure
    do i=1,2
       do j=1,ndimL
          K_sk(1:Qdof,i, i, j,j) = Re_1(1:Qdof)
       enddo
    enddo
  end subroutine Set_K_sk_Laplace

  !> compute viscous fluxes R_s, s=1,2 for scalar euation
  !> in integ nodes
  subroutine Set_R_s_Laplace(ndimL, nbDim, Qdof, w, Dw, Re_1, R_s)
    integer, intent(in) :: ndimL, nbDim, Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    !real, intent(in) :: Re_1                     ! inverse of Reynolds number
    real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
    integer :: k

    do k=1,ndimL
       R_s(1:Qdof, 1, k) = Re_1(1:Qdof) * Dw(1:Qdof, k, 1)
       R_s(1:Qdof, 2, k) = Re_1(1:Qdof) * Dw(1:Qdof, k, 2)
    end do

  end subroutine Set_R_s_Laplace

end module modelLaplace
