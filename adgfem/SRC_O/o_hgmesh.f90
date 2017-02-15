module hgmesh_mod
   use mesh_mod
   use element_mod

implicit none
   type, public, extends( mesh ) :: MeshHG_t
      integer :: HG                        ! maximal number of hanging nodes per edge

      contains

      procedure :: allocElem => allocElementsHG

   end type MeshHG_t

contains
!!!!!!!!!!!!!!!!!!  HG   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!§§
   subroutine allocElementsHG( this, nelem)
      class( MeshHG_t), intent(inout) :: this
      integer, intent(in) :: nelem

      print*, 'allocElementsHG'

         allocate( ElementHG_t :: this%elem(1:nelem) )


   end subroutine allocElementsHG


end module hgmesh_mod
