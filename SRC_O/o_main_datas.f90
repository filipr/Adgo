module data_mod
   use mesh_mod

   implicit none

   !> OOP global datas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> global file variables
  class(mesh), pointer :: grid, gridN, gridS, gridD, gridL


!  !> coordinates of nodes on curved boundary part
!  type(curved_boundary) :: curvbound

end module
