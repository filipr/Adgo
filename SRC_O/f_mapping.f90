!> operations with mappings \f$ F: \hat{K}\to K\f$ defined
!> in module geometry

module f_mapping

  use geometry
  use element_mod

  implicit none

  public:: SetF
  public:: ComputeF
  public:: ComputeDF

  public:: CheckElement
contains

  !> setting of mappings \f$ F: \hat{K}\to K\quad \forall K\in{\cal T}_h\f$, i.e.
  !> mapping of the reference element onto the physical one
  !>
  !> \f$ K= \f$ elem,  given by Langrangian nodes within  TRIANGLE or QUADRILATERALL,
  !> \f$ F \f$ can be linear, quadratic, cubic, ....,
  !> F(l,1:nbDim) means  the \f$(x,y)\f$ coordinates of the l\f$^{\rm th}\f$
  !> Langrangian nodes
  subroutine SetF(elem, nod, x)
    type(element), intent(inout):: elem   ! elem = element
    integer, intent(in) :: nod
    real, dimension(1:nod, 1:nbDim), intent(in) :: x
    integer :: deg
    real, dimension(:, :), pointer :: F

    F => elem%F%F(1:elem%F%dof, 1:nbDim)

    if(nod /= elem%F%dof) then
       print*,'Inconsistency in SetF:',nod, elem%F%dof
       stop
    endif

    elem%F%iFlin = .false.

    deg = elem%F%deg
    if(elem%type == 3 .and. nbDim == 2) then   !  TRIANGLES
       if(deg == 1) then
          F(1, 1:2) = x(1,1:2)
          F(2, 1:2) = x(2,1:2) - x(1,1:2)
          F(3, 1:2) = x(3,1:2) - x(1,1:2)

          elem%F%iFlin = .true.   ! linear mapping

       elseif(deg == 2) then
          F(1, 1:2) = x(1,1:2)
          F(2, 1:2) = -(x(2,1:2) -4*x(4,1:2) +3*F(1,1:2))
          F(3, 1:2) = -(x(3,1:2) -4*x(6,1:2) +3*F(1,1:2))
          F(4, 1:2) = 2*(x(2,1:2) + F(1,1:2) - 2*x(4,1:2))
          F(6, 1:2) = 2*(x(3,1:2) + F(1,1:2) - 2*x(6,1:2))
          F(5, 1:2) = 4*(x(5,1:2) - F(1,1:2) -(F(2,1:2)+F(3,1:2))/2  &
               - (F(4,1:2)+F(6,1:2))/4 )
       elseif(deg == 3) then
          F(1,1:2) = x(1,1:2)
          F(2,1:2) = -11./2.*x(1,1:2)+x(2,1:2)+9*x(4,1:2)-9./2.*x(5,1:2)
          F(3,1:2) = -11./2.*x(1,1:2)+x(3,1:2)-9./2.*x(8,1:2)+9*x(9,1:2)
          F(4,1:2) = 9*x(1,1:2)-9./2.*x(2,1:2)-45./2.*x(4,1:2)+18*x(5,1:2)
          F(5,1:2) = 18*x(1,1:2)-45./2.*x(4,1:2)+9./2.*x(5,1:2)-9./2.*x(6,1:2) &
               -9./2.*x(7,1:2)+9./2.*x(8,1:2)-45./2.*x(9,1:2)+27*x(10,1:2)
          F(6,1:2) = 9*x(1,1:2)-9./2.*x(3,1:2)+18*x(8,1:2)-45./2.*x(9,1:2)
          F(7,1:2) = -9./2.*x(1,1:2)+9./2.*x(2,1:2)+27./2.*x(4,1:2)-27./2.*x(5,1:2)
          F(8,1:2) = -27./2.*x(1,1:2)+27*x(4,1:2)-27./2.*x(5,1:2)+27./2.*x(6,1:2) &
               +27./2.*x(9,1:2)-27*x(10,1:2)
          F(9,1:2) = -27./2.*x(1,1:2)+27./2.*x(4,1:2)+27./2.*x(7,1:2)-27./2.*x(8,1:2)&
               +27*x(9,1:2)-27*x(10,1:2)
          F(10,1:2) = -9./2.*x(1,1:2)+9./2.*x(3,1:2)-27./2.*x(8,1:2)+27./2.*x(9,1:2)
       else
          print*,'Maximal mapping P_3 on triangles implemented'
          stop
       endif

    elseif(elem%type == 4 .and. nbDim == 3) then   !  TETRAHEDRA

       if(deg == 1) then
          F(1, 1:nbDim) = x(1,1:nbDim)
          F(2, 1:nbDim) = x(2,1:nbDim) - x(1,1:nbDim)
          F(3, 1:nbDim) = x(3,1:nbDim) - x(1,1:nbDim)
          F(4, 1:nbDim) = x(4,1:nbDim) - x(1,1:nbDim)

          elem%F%iFlin = .true.   ! linear mapping
       else
          print*,'Maximal mapping P_1 on tetrahedra implemented'
          stop
       endif

    elseif(elem%type == 4 .and. nbDim == 2 ) then   !  QUADRILATERALS

       if(deg == 1) then
          F(1, 1:2) = x(1,1:2)
          F(2, 1:2) = x(2,1:2) - x(1,1:2)
          F(3, 1:2) = x(3,1:2) - x(1,1:2)
          F(4, 1:2) = x(4,1:2) - x(2,1:2)  - x(3,1:2) + x(1,1:2)

          if(dot_product(F(4,1:2), F(4,1:2)) <= 1D-14) elem%F%iFlin = .true.   ! linear mapping
       elseif(deg == 2) then
          F(1, 1:2) = x(1,1:2)
          F(2, 1:2) = -3*x(1,1:2)+4*x(2,1:2)-x(3,1:2)
          F(3, 1:2) = -3*x(1,1:2)+4*x(4,1:2)-x(7,1:2)
          F(4, 1:2) = 9*x(1,1:2)-12*x(2,1:2)+3*x(3,1:2)-12*x(4,1:2) &
               +16*x(5,1:2)-4*x(6,1:2)+3*x(7,1:2)-4*x(8,1:2)+x(9,1:2)
          F(5, 1:2) = 2*x(1,1:2)-4*x(2,1:2)+2*x(3,1:2)
          F(6, 1:2) = -6*x(1,1:2)+12*x(2,1:2)-6*x(3,1:2)+8*x(4,1:2) &
               -16*x(5,1:2)+8*x(6,1:2)-2*x(7,1:2)+4*x(8,1:2)-2*x(9,1:2)
          F(7, 1:2) = 4*x(1,1:2)-8*x(2,1:2)+4*x(3,1:2)-8*x(4,1:2)+ &
               16*x(5,1:2)-8*x(6,1:2)+4*x(7,1:2)-8*x(8,1:2)+4*x(9,1:2)
          F(8, 1:2) = 2*x(1,1:2)-4*x(4,1:2)+2*x(7,1:2)
          F(9, 1:2) = -6*x(1,1:2)+8*x(2,1:2)-2*x(3,1:2)+12*x(4,1:2) &
               -16*x(5,1:2)+4*x(6,1:2)-6*x(7,1:2)+8*x(8,1:2)-2*x(9,1:2)

       elseif(deg == 3) then
          F(1,1:2) = x(1,1:2)
          F(2,1:2) = -5.5*x(1,1:2) +9*x(2,1:2) -4.5*x(3,1:2) +x(4,1:2)
          F(3,1:2) = -5.5*x(1,1:2) +9*x(5,1:2) -4.5*x(9,1:2) +x(13,1:2)
          F(4,1:2) = 121./4.*x(1,1:2) -49.5*x(2,1:2) +99./4.*x(3,1:2)-5.5*x(4,1:2) &
               -49.5*x(5,1:2) +81*x(6,1:2) -40.5*x(7,1:2) +9*x(8,1:2) &
               +99./4.*x(9,1:2) -40.5*x(10,1:2) +81./4.*x(11,1:2) -4.5*x(12,1:2) &
               -5.5*x(13,1:2) +9*x(14,1:2) -4.5*x(15,1:2) +x(16,1:2)
          F(5,1:2) = 9*x(1,1:2) -22.5*x(2,1:2) +18*x(3,1:2) -4.5*x(4,1:2)
          F(6,1:2) = -49.5*x(1,1:2) +495./4.*x(2,1:2) -99*x(3,1:2) +99./4.*x(4,1:2) &
               +81*x(5,1:2) -202.5*x(6,1:2) +162*x(7,1:2) -40.5*x(8,1:2)  &
               -40.5*x(9,1:2) +405./4.*x(10,1:2) -81*x(11,1:2) +81./4.*x(12,1:2) &
               +9*x(13,1:2) -22.5*x(14,1:2) +18*x(15,1:2) -4.5*x(16,1:2)
          F(7,1:2) = 81*x(1,1:2) -202.5*x(2,1:2) +162*x(3,1:2) -81./2.*x(4,1:2) &
               -405./2.*x(5,1:2) +2025./4.*x(6,1:2) -405*x(7,1:2) +405./4.*x(8,1:2) &
               +162*x(9,1:2) -405*x(10,1:2) +324*x(11,1:2) -81*x(12,1:2) &
               -81./2.*x(13,1:2) +405./4.*x(14,1:2) -81*x(15,1:2) +81./4.*x(16,1:2)
          F(8,1:2) = 9*x(1,1:2) -45./2.*x(5,1:2) +18*x(9,1:2) -4.5*x(13,1:2)
          F(9,1:2) = -99./2.*x(1,1:2) +81*x(2,1:2) -81./2.*x(3,1:2) +9*x(4,1:2) &
               +495./4.*x(5,1:2) -405./2.*x(6,1:2) +405./4.*x(7,1:2) -45./2.*x(8,1:2) &
               -99*x(9,1:2) +162*x(10,1:2) -81*x(11,1:2) +18*x(12,1:2) &
               +99./4.*x(13,1:2) -81./2.*x(14,1:2) +81./4.*x(15,1:2) -4.5*x(16,1:2)
          F(10,1:2) = -4.5*x(1,1:2) +27./2.*x(2,1:2) -27./2.*x(3,1:2) &
               +9./2.*x(4,1:2)
          F(11,1:2) = 99./4.*x(1,1:2) -297./4.*x(2,1:2) +297./4.*x(3,1:2) &
               -99./4.*x(4,1:2) -81./2.*x(5,1:2) +243./2.*x(6,1:2) -243./2.*x(7,1:2) &
               +81./2.*x(8,1:2) +81./4.*x(9,1:2) -243./4.*x(10,1:2) +243./4.*x(11,1:2) &
               -81./4.*x(12,1:2) -4.5*x(13,1:2) +27./2.*x(14,1:2) -27./2.*x(15,1:2) &
               +9./2.*x(16,1:2)
          F(12,1:2) = -81./2.*x(1,1:2) +243./2.*x(2,1:2) -243./2.*x(3,1:2) &
               +81./2*x(4,1:2) +405./4.*x(5,1:2) -1215./4.*x(6,1:2) +1215./4.*x(7,1:2) &
               -405./4.*x(8,1:2) -81*x(9,1:2) +243*x(10,1:2) -243*x(11,1:2) &
               +81*x(12,1:2) +81./4.*x(13,1:2) -243./4.*x(14,1:2) +243./4.*x(15,1:2) &
               -81./4.*x(16,1:2)
          F(13,1:2) = 81./4.*x(1,1:2) -243./4.*x(2,1:2) +243./4.*x(3,1:2) &
               -81./4.*x(4,1:2) -243./4.*x(5,1:2) +729./4.*x(6,1:2) -729./4.*x(7,1:2)&
               +243./4.*x(8,1:2)+ 243./4.*x(9,1:2)-729./4. *x(10,1:2) &
               +729./4.*x(11,1:2) -243./4.*x(12,1:2) -81./4.*x(13,1:2) &
               +243./4.*x(14,1:2) -243./4.*x(15,1:2) +81./4.*x(16,1:2)
          F(14,1:2) = -4.5*x(1,1:2) +27./2.*x(5,1:2) -27./2.*x(9,1:2)+ 9./2.*x(13,1:2)
          F(15,1:2) = 99./4.*x(1,1:2) -81./2.*x(2,1:2) +81./4.*x(3,1:2) &
               -4.5*x(4,1:2) -297./4.*x(5,1:2) +243./2.*x(6,1:2) -243./4.*x(7,1:2) &
               +27./2.*x(8,1:2)+297./4.*x(9,1:2) -243./2.*x(10,1:2) +243./4.*x(11,1:2)&
               -27./2.*x(12,1:2) -99./4.*x(13,1:2) +81./2.*x(14,1:2) -81./4.*x(15,1:2)&
               +9./2.*x(16,1:2)
          F(16,1:2) = -81./2.*x(1,1:2) +405./4.*x(2,1:2)- 81*x(3,1:2) +81./4.*x(4,1:2)&
               +243./2.*x(5,1:2) -1215./4.*x(6,1:2) +243*x(7,1:2) -243./4.*x(8,1:2)  &
               -243./2.*x(9,1:2)+ 1215./4.*x(10,1:2) -243*x(11,1:2) +243./4.*x(12,1:2) &
               +81./2.*x(13,1:2) -405./4.*x(14,1:2) +81*x(15,1:2) -81./4.*x(16,1:2)
       else
          print*,'Maximal mapping P_3 on quadrilateralls implemented'
          stop
       endif
    else
       print*,'Only triangles or quadrilaterals implemented in SetF'
       stop
    endif

  end subroutine SetF

  !> OLD, it is replaced by elem%ComputeF in o_element.f90 - only for 2D
  !> evaluation of \f$ F(x_i)\in R^2,\ x_i\in \hat{K},\ i=1..nod,  \f$
  !>
  !> xi= \f$x_i\in \hat{K}\f$ are arbitrary nodes within reference elements,
  !> \f$ F \f$ mapping of reference element \f$\hat{K}\f$ to the actual one \f$K\f$=elem,
  !> F can be linear, quadratic, cubic, ....
  subroutine ComputeF(elem, nod, xi, Fx)
    type(element), intent(in):: elem   ! elem = element
    integer, intent(in) :: nod            ! number of nodes
    real, dimension(1:nod, 1:nbDim), intent(in) :: xi
    real, dimension(1:nod, 1:nbDim), intent(out) :: Fx
    integer :: i, deg

    real, dimension(:, :), pointer :: F

    F => elem%F%F(1:elem%F%dof, 1:nbDim)
    deg = elem%F%deg

    if (nbDim == 2) then   ! 2D elements
       if(elem%type == 3) then  ! TRIANGLES
          do i=1, 2
             ! absolute term
             Fx(1:nod, i) = F(1, i)

             if(deg >= 1) then
                ! linear terms
                Fx(1:nod, i) = Fx(1:nod, i)  &
                     + F(2, i)*xi(1:nod, 1) + F(3, i)*xi(1:nod, 2)
                if(deg >= 2) then
                   ! quadratic terms
                   Fx(1:nod, i) = Fx(1:nod, i)  &
                        + F(4, i) * xi(1:nod, 1) * xi(1:nod, 1) &
                        + F(5, i) * xi(1:nod, 1) * xi(1:nod, 2) &
                        + F(6, i) * xi(1:nod, 2) * xi(1:nod, 2)
                   if(deg >= 3) then
                      ! cubic terms
                      Fx(1:nod, i) = Fx(1:nod, i)  &
                           + F(7, i) * xi(1:nod, 1) * xi(1:nod, 1) * xi(1:nod, 1) &
                           + F(8, i) * xi(1:nod, 1) * xi(1:nod, 1) * xi(1:nod, 2) &
                           + F(9, i) * xi(1:nod, 1) * xi(1:nod, 2) * xi(1:nod, 2) &
                           + F(10,i) * xi(1:nod, 2) * xi(1:nod, 2) * xi(1:nod, 2)
                      if(deg >= 4) then
                         ! 4-th order terms
                         Fx(1:nod, i) = Fx(1:nod, i)  &
                              + F(11,i)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1) &
                              + F(12,i)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2) &
                              + F(13,i)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2) &
                              + F(14,i)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,2) &
                              + F(15,i)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,2)
                         if(deg >=5) then
                            print*,'Only polynoms up to order 4 are implemented',&
                                 'in computeF for triangles'
                            stop
                         endif
                      endif
                   endif
                endif
             endif
          enddo

       elseif(elem%type == 4) then   ! QUADRILATERALS
          do i=1, 2
             ! absolute term
             Fx(1:nod, i) = F(1, i)

             if(deg >= 1) then
                ! bilinear terms
                Fx(1:nod, i) = Fx(1:nod, i) + F(2, i)*xi(1:nod, 1) &
                     + F(3, i)*xi(1:nod, 2) + F(4, i)*xi(1:nod, 1)*xi(1:nod, 2)
                if(deg >= 2) then
                   ! biquadratic terms
                   Fx(1:nod, i) = Fx(1:nod, i)  &
                        + F(5, i)*xi(1:nod, 1)*xi(1:nod, 1) &
                        + F(6, i)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 2) &
                        + F(7, i)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 2)*xi(1:nod, 2) &
                        + F(8, i)*xi(1:nod, 2)*xi(1:nod, 2) &
                        + F(9, i)*xi(1:nod, 1)*xi(1:nod, 2)*xi(1:nod, 2)
                   if(deg >= 3) then
                      ! bicubic terms
                      Fx(1:nod, i) = Fx(1:nod, i)  &
                           + F(10, i)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 1) &
                           + F(11, i)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 2) &
                           + F(12, i)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 2) &
                           *xi(1:nod, 2) &
                           + F(13, i)*xi(1:nod, 1)**3 *xi(1:nod, 2)**3 &
                           + F(14, i)*xi(1:nod, 2)*xi(1:nod, 2)*xi(1:nod, 2) &
                           + F(15, i)*xi(1:nod, 1)*xi(1:nod, 2)*xi(1:nod, 2)*xi(1:nod, 2) &
                           + F(16, i)*xi(1:nod, 1)*xi(1:nod, 1)*xi(1:nod, 2)*xi(1:nod, 2) &
                           *xi(1:nod, 2)
                      !print*,'# Nontested  Eval_QPolynoms for deg = 3 !!'
                      if(deg >=4) then
                         print*,'Only Qpolynoms up to order 3 are implemented in ComputeF_4'
                         stop
                      endif
                   endif
                endif
             endif
          enddo

       else
          print*,'Only triagles and quarilaterals are implemented in ComputeF'
          stop
       endif

    elseif(nbDim == 3) then ! 3D elements

       if(elem%type == 4) then  ! TRETRAHEDRA
          do i=1, nbDim
             ! absolute term
             Fx(1:nod, i) = F(1, i)

             if(deg >= 1) then
                ! linear terms
                Fx(1:nod, i) = Fx(1:nod, i)  &
                     + F(2, i)*xi(1:nod, 1) + F(3, i)*xi(1:nod, 2) + F(4, i)*xi(1:nod, 3)

                !    if(deg >= 2) then
                !     ! quadratic terms
                !     Fx(1:nod, i) = Fx(1:nod, i)  &
                !             + F(4, i) * xi(1:nod, 1) * xi(1:nod, 1) &
                !             + F(5, i) * xi(1:nod, 1) * xi(1:nod, 2) &
                !             + F(6, i) * xi(1:nod, 2) * xi(1:nod, 2)
                if(deg >=2) then
                   print*,'Only polynoms up to order 1 are implemented',&
                        'in computeF for tetrahedra'
                   stop
                endif
                !                    endif
             endif
          enddo
       else
          print*,'Only tetrahedra are implemented in ComputeF'
          stop
       endif
    endif

  end subroutine ComputeF

   !> OLD, it is replaced by elem%ComputeDF in o_element.f90 - only for 2D
  !> evaluation of \f$ \frac{D F(x_i)}{D \hat{x}} \in R^{2\times 2},
  !> \ x_i\in \hat{K},\ i=1..nod,  \f$
  !>
  !> xi= \f$x_i\in \hat{K}\f$ are arbitrary nodes within reference elements,
  !> \f$ F \f$ mapping of reference element \f$\hat{K}\f$ to the actual one \f$K\f$=elem,
  !> F can be linear, quadratic, cubic, ....
  subroutine ComputeDF(elem, nod, xi, DF  )
    type(element), intent(in):: elem   ! elem = element
    integer, intent(in) ::  nod
    real, dimension(1:nod, 1:nbDim), intent(in) :: xi
    real, dimension(1:nod, 1:nbDim, 1:nbDim ), intent(out) :: DF
    integer :: deg, i

    real, dimension(:, :), pointer :: F

    F => elem%F%F(1:elem%F%dof, 1:nbDim)

    deg = elem%F%deg

    !if(elem%i == 15 ) then
    !   do i=1,nod
    !      write(99,*) xi(i, 1:2)
    !   enddo
    !endif


    if (nbDim==2) then   ! 2D case

       if(elem%type == 3) then  ! TRIANGLES
          if(deg >= 1) then
             ! linnear mapping
             DF(1:nod,1,1) = F(2,1)
             DF(1:nod,1,2) = F(3,1)
             DF(1:nod,2,1) = F(2,2)
             DF(1:nod,2,2) = F(3,2)

             if(deg >= 2) then
                ! quadratic mapping
                DF(1:nod,1,1) = DF(1:nod,1,1) + 2*F(4,1)*xi(1:nod,1)+ F(5,1)*xi(1:nod,2)
                DF(1:nod,1,2) = DF(1:nod,1,2) + 2*F(6,1)*xi(1:nod,2)+ F(5,1)*xi(1:nod,1)
                DF(1:nod,2,1) = DF(1:nod,2,1) + 2*F(4,2)*xi(1:nod,1)+ F(5,2)*xi(1:nod,2)
                DF(1:nod,2,2) = DF(1:nod,2,2) + 2*F(6,2)*xi(1:nod,2)+ F(5,2)*xi(1:nod,1)

                if(deg >= 3) then
                   ! cubic mapping
                   DF(1:nod,1,1) = DF(1:nod,1,1) + 3*F(7,1)*xi(1:nod,1)*xi(1:nod,1) &
                        + 2*F(8,1)*xi(1:nod,1)*xi(1:nod,2) &
                        + F(9,1)*xi(1:nod,2)*xi(1:nod,2)
                   DF(1:nod,1,2) = DF(1:nod,1,2) + 3*F(10,1)*xi(1:nod,2)*xi(1:nod,2) &
                        + 2*F(9,1)*xi(1:nod,1)*xi(1:nod,2) &
                        + F(8,1)*xi(1:nod,1)*xi(1:nod,1)
                   DF(1:nod,2,1) = DF(1:nod,2,1) + 3*F(7,2)*xi(1:nod,1)*xi(1:nod,1) &
                        + 2*F(8,2)*xi(1:nod,1)*xi(1:nod,2) &
                        + F(9,2)*xi(1:nod,2)*xi(1:nod,2)
                   DF(1:nod,2,2) = DF(1:nod,2,2) + 3*F(10,2)*xi(1:nod,2)*xi(1:nod,2) &
                        + 2*F(9,2)*xi(1:nod,1)*xi(1:nod,2) &
                        + F(8,2)*xi(1:nod,1)*xi(1:nod,1)
                   !      print*,'P_3 approximation of boundary is implemented for', &
                   !                  ' Jacobian but not tested !!!'

                   if(deg >=4) then
                      print*,'ComputeDF_3 for deg > 3 is not implemented'
                      stop
                   endif
                endif
             endif
          endif

       elseif(elem%type == 4) then   ! QUADRILATERALS

          if(deg >= 1) then
             ! bi-linnear mapping
             DF(1:nod,1,1) = F(4,1)*xi(1:nod,2) + F(2,1)
             DF(1:nod,1,2) = F(4,1)*xi(1:nod,1) + F(3,1)

             DF(1:nod,2,1) = F(4,2)*xi(1:nod,2) + F(2,2)
             DF(1:nod,2,2) = F(4,2)*xi(1:nod,1) + F(3,2)

             if(deg >= 2) then
                ! bi-quadratic mapping
                DF(1:nod,1,1) = DF(1:nod,1,1) + &
                     2*F(5,1)*xi(1:nod,1)+ &
                     2*F(6,1)*xi(1:nod,1)*xi(1:nod,2) + &
                     2*F(7,1)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2) + &
                     F(9,1)*xi(1:nod,2)*xi(1:nod,2)

                DF(1:nod,1,2) = DF(1:nod,1,2) + &
                     F(6,1)*xi(1:nod,1)*xi(1:nod,1)+ &
                     2*F(7,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2) + &
                     2*F(8,1)*xi(1:nod,2) + &
                     2*F(9,1)*xi(1:nod,1)*xi(1:nod,2)

                DF(1:nod,2,1) = DF(1:nod,2,1) + &
                     2*F(5,2)*xi(1:nod,1)+ &
                     2*F(6,2)*xi(1:nod,1)*xi(1:nod,2) + &
                     2*F(7,2)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2) + &
                     F(9,2)*xi(1:nod,2)*xi(1:nod,2)

                DF(1:nod,2,2) = DF(1:nod,2,2) + &
                     F(6,2)*xi(1:nod,1)*xi(1:nod,1)+ &
                     2*F(7,2)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2) + &
                     2*F(8,2)*xi(1:nod,2) + &
                     2*F(9,2)*xi(1:nod,1)*xi(1:nod,2)

                if(deg >= 3) then
                   ! bi-cubic mapping
                   DF(1:nod,1,1) = DF(1:nod,1,1) + &
                        3*F(10,1)*xi(1:nod,1)*xi(1:nod,1)+ &
                        3*F(11,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2) + &
                        3*F(12,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2) + &
                        3*F(13,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2) &
                        *xi(1:nod,2)  + &
                        F(15,1)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,2) + &
                        2*F(16,1)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,2)


                   DF(1:nod,1,2) = DF(1:nod,1,2) + &
                        F(11,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)+ &
                        2*F(12,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2) + &
                        3*F(13,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2)&
                        *xi(1:nod,2)+ &
                        3*F(14,1)*xi(1:nod,2)*xi(1:nod,2) + &
                        3*F(15,1)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,1) + &
                        3*F(16,1)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,1)*xi(1:nod,1)

                   DF(1:nod,2,1) = DF(1:nod,2,1) + &
                        3*F(10,2)*xi(1:nod,1)*xi(1:nod,1)+ &
                        3*F(11,2)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2) + &
                        3*F(12,2)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2) + &
                        3*F(13,2)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2) &
                        *xi(1:nod,2)  + &
                        F(15,2)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,2) + &
                        2*F(16,2)*xi(1:nod,1)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,2)


                   DF(1:nod,2,2) = DF(1:nod,2,2) + &
                        F(11,2)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)+ &
                        2*F(12,2)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2) + &
                        3*F(13,2)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,1)*xi(1:nod,2)&
                        *xi(1:nod,2)+ &
                        3*F(14,2)*xi(1:nod,2)*xi(1:nod,2) + &
                        3*F(15,2)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,1) + &
                        3*F(16,2)*xi(1:nod,2)*xi(1:nod,2)*xi(1:nod,1)*xi(1:nod,1)

                   !print*,'#P_3 approximation of Q boundary is implemented for', &
                   !     ' Jacobian but not tested !!!'
                   if(deg >=4) then
                      print*,'Compute_DF4 for deg > 3 is not implemented'
                      stop
                   endif
                endif
             endif
          endif

       else
          print*,'Only triagles and quarilaterals are implemented in ComputeDF'
          stop

       endif

    elseif (nbDim==3) then  ! 3D case

       if(elem%type == 4) then  ! TETRAHEDRA
          if(deg >= 1) then
             ! linnear mapping
             do i=1,3
                DF(1:nod,i,1) = F(2,i)
                DF(1:nod,i,2) = F(3,i)
                DF(1:nod,i,3) = F(4,i)
             enddo

             !if(deg >= 2) then
             !   ! quadratic mapping
             !   DF(1:nod,1,1) = DF(1:nod,1,1) + 2*F(4,1)*xi(1:nod,1)+ F(5,1)*xi(1:nod,2)
             !   DF(1:nod,1,2) = DF(1:nod,1,2) + 2*F(6,1)*xi(1:nod,2)+ F(5,1)*xi(1:nod,1)
             !   DF(1:nod,2,1) = DF(1:nod,2,1) + 2*F(4,2)*xi(1:nod,1)+ F(5,2)*xi(1:nod,2)
             !   DF(1:nod,2,2) = DF(1:nod,2,2) + 2*F(6,2)*xi(1:nod,2)+ F(5,2)*xi(1:nod,1)
             if(deg >=2) then
                print*,'ComputeDF_3 for deg > 1 is not implemented'
                stop
             endif
             !endif
          endif
       else
          print*,'Only tetrahedra for 3D are implemented in ComputeDF'
          stop
       endif

    endif
  end subroutine ComputeDF

  !> draw a results of mapping \f$ F: \hat{K}\to K\f$ to file 'fort.(ifile)'
  !> visualizable by gnuplot, \f$ K= \f$ elem
  subroutine CheckElement(elem, ifile)
    type(element), intent(inout):: elem   ! elem = element
    integer, intent(in) :: ifile
    real, dimension(:,:), allocatable :: xi, Fxi
    integer :: i,j, itest, it

    itest = 20

    if(elem%type == 3) then
       allocate(xi(1:(itest+1)*(itest+2)/2,1:nbDim) )
       allocate(Fxi(1:(itest+1)*(itest+2)/2,1:nbDim) )
       it = 0
       do i=0,itest
          do j=0,itest-i
             it = it + 1
             xi(it,1) = 1.*i/itest
             xi(it,2) = 1.*j/itest
          enddo
       enddo

    elseif(elem%type == 4) then
       allocate(xi(1:(itest+1)**2,1:nbDim) )
       allocate(Fxi(1:(itest+1)**2,1:nbDim) )
       it = 0
       do i=0,itest
          do j=0,itest
             it = it + 1
             xi(it,1) = 1.*i/itest
             xi(it,2) = 1.*j/itest
          enddo
       enddo
    else
       ! no work
       allocate(xi(1:1,1:nbDim) )
       allocate(Fxi(1:1,1:nbDim) )
    endif

    call ComputeF(elem, it, xi(1:it,1:nbDim), Fxi(1:it,1:nbDim) )

    do i=1,it
       if(xi(i,1)*xi(1,2)*(1-xi(i,1) - xi(i,2)) < 1E-12) then
          write(ifile,*) Fxi(i,1:nbDim),xi(i,1:nbDim),i, xi(i,1)*xi(1,2)*(1-xi(i,1) - xi(i,2))
       endif
    enddo
    write(ifile,'(x)')

    deallocate(xi,Fxi)
  end subroutine CheckElement


end module f_mapping
