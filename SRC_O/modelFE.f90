!> definition of 2D models for conforming FEM
module modelFE

  use main_data
  use f_mapping
  use mesh_oper
!  use blocks_integ

  implicit none

  public:: Scal_Lin
  public:: LinearConvection
  public:: LinearReaction
contains
    !FERROR Is this used? for FEM
    !> function for linear convection and reaction
  function Scal_Lin(imod, x, y, t)
    real :: Scal_Lin
    integer:: imod
    real ::  x,y,t

    if(imod == 0) Scal_Lin = 0.
    !if(imod == 1) Scal_Lin = -1.
    !if(imod == 2) Scal_Lin = -1.

    !if(imod == 1) Scal_Lin = 1.
    !if(imod == 2) Scal_Lin = 1.

    if(imod == 1) Scal_Lin = 1.
    if(imod == 2) Scal_Lin = 0.

    !Scal_Lin = 0.
  end function Scal_Lin


  !> linear convection (for implicit discretization)
  subroutine LinearConvection(Qdof, ndim, x, t, conv)
    integer, intent(in) :: Qdof, ndim
    real, dimension(1:Qdof, 1:nbDim), intent(in):: x        ! x cooedinates
    real, intent(in) :: t                               ! time
    real, dimension(1:Qdof, 1:nbDim, 1:ndim, 1:ndim ), intent(inout) :: conv
    integer :: l

    do l = 1, Qdof
       conv(l, 1, 1:ndim, 1:ndim) = Scal_Lin(1, x(l,1), x(l, 2), t)
       conv(l, 2, 1:ndim, 1:ndim) = Scal_Lin(2, x(l,1), x(l, 2), t)
    enddo

  end subroutine LinearConvection


  !> linear convection (for implicit discretization)
  subroutine LinearReaction(Qdof, ndim, x, t, reac)
    integer, intent(in) :: Qdof, ndim
    real, dimension(1:Qdof, 1:nbDim), intent(in):: x        ! x cooedinates
    real, intent(in) :: t                               ! time
    real, dimension(1:Qdof,1:ndim, 1:ndim ), intent(inout) :: reac
    integer :: l

    do l = 1,Qdof
       reac(l, 1:ndim, 1:ndim) = Scal_Lin(0, x(l,1), x(l, 2), t)
    enddo

  end subroutine LinearReaction

end module modelFE
