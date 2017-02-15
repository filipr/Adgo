!> definition of 2D models which are simulated: Euler, Navier-Stokes
module modelTurb2e

  use main_data
  use f_mapping
  use mesh_oper
  use blocks_integ
  use model2DNS
!  use model3DNS

  implicit none

  public:: Set_f_s_Turb2e
  public:: Set_A_s_Turb2e
  public:: Set_Ppm_Turb2e
  public:: Set_R_s_Turb2e
  public:: Set_K_sk_Turb2e

contains

  !> compute the Euler fluxes \f$ f_s({\bf w}),\quad s=1,2\f$
  subroutine Set_f_s_Turb2e(ndimL, nbDim, Qdof, w, f_s, xi, ie)
    integer, intent(in) :: Qdof, nbDim, ndimL
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
                                               ! matrices A_s in  -- " --
    real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
    integer, intent(in) :: ie
    
    
    ! laminar terms
    if(nbDim == 2) then
       call Set_f_s_Euler(nbDim+2, nbDim, Qdof, w(1:Qdof,1:nbDim+2), &
            f_s(1:Qdof, 1:nbDim, 1:nbDim+2), xi, ie )
    elseif(nbDim == 3) then
      stop 'Stopped in Set_f_s_Turb2e - 3D not implemented'
!       call Set_f_s_Euler3D(nbDim+2, Qdof, w(1:Qdof,1:nbDim+2), &
!            f_s(1:Qdof, 1:nbDim, 1:nbDim+2) )
    endif

    !turbulence terms
    f_s(1:Qdof, 1:nbDim, nbDim+3:ndimL) = 0.

  end subroutine Set_f_s_Turb2e


  !> compute matrices \f$ A_s = \frac{D{\bf f_s}({\bf w})}{D{\bf w}},\quad s=1,2\f$
  !> for Euler fluxes
  subroutine Set_A_s_Turb2e(ndimL, nbDim, Qdof, w, A_s, xi, ie)
    integer, intent(in) :: ndimL, nbDim, Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: A_s
                                               ! matrices A_s in  -- " --
    real, dimension(1:Qdof,1 :nbDim), intent(in) :: xi
    integer, intent(in) :: ie

    ! laminar terms
    if(nbDim == 2) then
       call Set_A_s_Euler(nbDim+2, nbDim, Qdof, w(1:Qdof,1:nbDim+2), &
            A_s(1:Qdof, 1:nbDim, 1:nbDim+2, 1:nbDim+2), xi, ie )
    elseif(nbDim == 3) then
     stop 'Stopped in Set_A_s_Turb2e - 3D not implemented'
     !  call Set_A_s_Euler3D(nbDim+2, Qdof, w(1:Qdof,1:nbDim+2), &
     !       A_s(1:Qdof, 1:nbDim, 1:nbDim+2, 1:nbDim+2) )
    endif

    !turbulence terms
    A_s(1:Qdof, 1:nbDim, nbDim+3:ndimL, :) = 0.
    A_s(1:Qdof, 1:nbDim, :, nbDim+3:ndimL) = 0.

  end subroutine Set_A_s_Turb2e


  !> compute matrices 4x4 K_sk, s,k=1,2 for N.-S. equations
  !> in integ nodes
  subroutine Set_K_sk_Turb2e(ndimL, Qdof, w, Dw, Re_1, K_sk, xi)
    integer, intent(in) :: Qdof, ndimL
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes, not USED
    real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    !real, intent(in) :: Re_1                     ! inverse of Reynolds number
    real, dimension(1:Qdof,1:nbDim,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: K_sk
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
    real :: v1, v2, E, one_over_Rew1, gamPr
    integer :: i

    ! laminar terms
    if(nbDim == 2) then
       call Set_K_sk_NS(nbDim+2, nbDim, iRe, Qdof, w(1:Qdof,1:nbDim+2), Dw(1:Qdof,1:nbDim+2,1:nbDim), Re_1, &
            K_sk(1:Qdof, 1:nbDim, 1:nbDim, 1:nbDim+2, 1:nbDim+2),xi )
    elseif(nbDim == 3) then
       call Set_K_sk_NS(nbDim+2, nbDim, iRe, Qdof, w(1:Qdof,1:nbDim+2), Dw(1:Qdof,1:nbDim+2,1:nbDim), Re_1, &
            K_sk(1:Qdof, 1:nbDim, 1:nbDim, 1:nbDim+2, 1:nbDim+2), xi )
    endif

    !turbulence terms
    K_sk(1:Qdof, 1:nbDim, 1:nbDim, nbDim+3:ndimL, :) = 0.
    K_sk(1:Qdof, 1:nbDim, 1:nbDim, :, nbDim+3:ndimL) = 0.

  end subroutine Set_K_sk_Turb2e


  !> compute viscous fluxes R_s, s=1,2 for N.-S. equations
  !> in integ nodes
  subroutine Set_R_s_Turb2e(ndimL, nbDim, iRe, Qdof, w, Dw, Re_1, R_s, xi)
    integer, intent(in) :: ndimL, nbDim, iRe, Qdof
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
    real, dimension(1:iRe, 1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
    !real, intent(in) :: Re_1                     ! inverse of Reynolds number
    real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
    real, dimension(1:Qdof, 1:nbDim), intent(in):: xi ! physical coordinates
    real, dimension(:), allocatable :: u, v, oRe, e

    ! laminar terms
    if(nbDim == 2) then
       call Set_R_s_NS(nbDim+2, nbDim, iRe, Qdof, w(1:Qdof,1:nbDim+2), Dw(1:Qdof,1:nbDim+2, 1:nbDim), &
            Re_1, R_s(1:Qdof, 1:nbDim, 1:nbDim+2), xi )
    elseif(nbDim == 3) then
      stop 'Stopped in Set_Ppm_Turb2e - 3D not implemented'
      ! call Set_R_s_NS3D(Qdof, w(1:Qdof,1:nbDim+2), Dw(1:Qdof,1:nbDim+2, 1:nbDim), &
      !      Re_1, R_s(1:Qdof, 1:nbDim, 1:nbDim+2) )
    endif

    !turbulence terms
    R_s(1:Qdof, 1:nbDim, nbDim+3:ndimL) = 0.


  end subroutine Set_R_s_Turb2e


  !> compute matrices
  !> \f$ P^{\pm} = \left(\frac{D({\bf f_1}({\bf w})n_1+{\bf f_2}({\bf w})n_2}{D{\bf w}}
  !>  \right)^{\pm}\f$
  !> for the Euler equation
  subroutine Set_Ppm_Turb2e(grid,ndimL, Qdof, w, n, xi, Ppm, one_over_area, elem)
    class(mesh), intent(in) :: grid
    integer, intent(in) :: Qdof, ndimL
    type(element), intent(inout) :: elem
    real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
    real, intent(in) :: one_over_area
    real, dimension(1:Qdof,1:nbDim,1:ndimL,1:ndimL), intent(inout) :: Ppm
                                               ! matrices Ppm in  -- " --
    real, dimension(1:Qdof, 1:nbDim), intent(in) :: n   ! outer normal
    real, dimension(1:nbDim),intent(in) ::  xi                    ! node on the edge?

    ! laminar terms
    if(nbDim == 2) then
       call Set_Ppm_Euler(nbDim+2, nbDim, Qdof, w(1:Qdof,1:nbDim+2), n, xi,  &
            Ppm(1:Qdof, 1:nbDim, 1:nbDim+2, 1:nbDim+2), one_over_area, elem )
    elseif(nbDim == 3) then
         stop 'Stopped in Set_Ppm_Turb2e - 3D not implemented'
     !   call Set_Ppm_Euler3D(nbDim+2, Qdof, w(1:Qdof,1:nbDim+2), n, xi,  &
     !       Ppm(1:Qdof, 1:nbDim, 1:nbDim+2, 1:nbDim+2), one_over_area, elem%i )
    endif

    !turbulence terms
    Ppm(1:Qdof, 1:nbDim, nbDim+3:ndimL, :) = 0.
    Ppm(1:Qdof, 1:nbDim, :, nbDim+3:ndimL) = 0.

  end subroutine Set_Ppm_Turb2e

  subroutine Set_Ppm_Turb2e_Slip(ndimL, nbDim, e_Qdof, w_ein, n_e, Ppm )
    ! compute matrix Pp on a slip boundary
    integer, intent(in) :: e_Qdof, ndimL, nbDim
    real, dimension(1:e_Qdof, 1:ndimL), intent(in):: w_ein !state  w in #Qdof nodes
    real, dimension(1:e_Qdof,2:nbDim+1,1:ndimL), intent(inout) :: Ppm
                                               ! matrices Ppm in  -- " --
    real, dimension(1:e_Qdof,1:nbDim), intent(in) :: n_e      ! outer normal
    real :: kappa, kappa1
    real :: v(2), vv, nn(2)
    integer :: ie, i, j



    ! laminar terms
    if(nbDim == 2) then
       call Set_Ppm_Euler_Slip(nbDim+2, nbDim, e_Qdof, w_ein(1:e_Qdof,1:nbDim+2), n_e,  &
            Ppm(1:e_Qdof, 2:nbDim+1, 1:nbDim+2)  )
    elseif(nbDim == 3) then
       print*,'To be done in modelTurb2e.f90 !! '
    !   call Set_Ppm_Euler3D_Slip(nbDim+2, e_Qdof, w_ein(1:e_Qdof,1:nbDim+2), n_e,  &
    !        Ppm(1:e_Qdof, 2:nbDim+1, 1:nbDim+2)  )
    endif

    !turbulence terms
    Ppm(1:e_Qdof, 2:nbDim+1, nbDim+3:ndimL) = 0.

  end subroutine Set_Ppm_Turb2e_Slip




end module modelTurb2e
