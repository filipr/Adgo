      subroutine Set_Ppm( ndimL, nbDim, Qdof, w, n, xi, Ppm, one_over_area, elem)
         import :: element
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: Ppm
                                               ! matrices Ppm in  -- " --
         real, dimension(1:Qdof, 1:nbDim), intent(in) :: n   ! outer normal
         real, dimension(1:Qdof, 1:nbDim),intent(in) ::  xi                    ! node on the edge?
         real, intent(in), optional :: one_over_area
         type(element), intent(inout), optional :: elem
      end subroutine

      subroutine Set_R_s(ndimL, nbDim, Qdof, w, Dw, Re_1, R_s)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
      end subroutine Set_R_s


      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
      end subroutine Set_f_s

      subroutine Set_A_s(ndimL, nbDim, Qdof, w, A_s, xi)
         integer, intent(in) :: Qdof, nbdim, ndimL
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: A_s
         ! matrices A_s in  -- " --
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
      end subroutine


      subroutine Set_K_sk(ndimL, nbDim,Qdof, w, Dw, Re_1, K_sk)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,ndimL), intent(inout) :: K_sk
       end subroutine Set_K_sk

       subroutine Set_S(ndimL, nbDim, Qdof, xi, w, Dw, S)
          integer, intent(in) :: ndimL, nbDim, Qdof
          real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
          real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
          real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
          real, dimension(1:Qdof, 1:ndimL), intent(inout) :: S
       end subroutine

       subroutine Set_DS(ndimL, nbDim, Qdof, xi, w, Dw, DS)
          integer, intent(in) :: ndimL, nbDim, Qdof
          real, dimension(1:Qdof, 1:nbDim), intent(in) :: xi
          real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
          real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
          real, dimension(1:Qdof, 1:ndimL, 1:ndimL), intent(inout) :: DS
       end subroutine

       subroutine Set_NumFlux(ndimL, nbDim, Qdof, wi, wj, nc, xi, H, area_1, ie )
          integer, intent(in) :: Qdof, ndimL, nbDim, ie
          real, dimension(1:Qdof, 1:ndimL), intent(in):: wi, wj ! state  w in integ nodes
          real, dimension(1:Qdof, 1:nbDim), intent(in):: nc        ! outer normal in integ nodes
          real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
          real, dimension(1:Qdof,1:ndimL), intent(inout) :: H   ! numer flux H in  -- " --
          real, intent(in) :: area_1
      end subroutine


      subroutine BasisType(Qdof, xi, len, dof,  phi, Dphi)
         integer, intent(in) :: Qdof          ! number of integ nodes
         real, dimension (1:Qdof, :), intent(in) :: xi   ! coordinates of integ. node
         integer, intent(in) :: len           ! triangle/square
         integer, intent(in) :: dof           ! number of test functions
         real, dimension(1:dof,1:Qdof) :: phi        ! value of test function
         real, dimension(1:dof,:,1:Qdof), optional:: Dphi  ! value of derivatives of test function
      end subroutine BasisType
