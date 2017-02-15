!> types and data for ANGENER
module AMAdata

  type AMAmesh
     integer :: ndim
     integer :: numel
     real :: p
     real :: epsilon1
     integer :: ityp
     real :: pos
     real :: pos1
     real :: posW
     integer :: adapt_level
     real, dimension (1:2,1:2) :: xper
     integer, dimension (1:2,1:2) :: iper
     integer :: ifv, icrack, iwa
     integer, dimension(:), allocatable :: iwall
     real :: glmin, glmax, errrez
     real :: Re
     real :: position_tolerance
     real, dimension (1:2,1:2) :: xte
     integer :: ifig, ifig1, ifig2
     integer :: melem, nelem, mpoin, npoin, maxdeg, ipoint
     integer :: nbelm, mbelm, nbc, nbp
     real, dimension(:), allocatable :: x, y, xold, yold, ra, supp
     real, dimension(:), allocatable :: rga, rgb, rgc, dx, dy, area, tria, xb, yb, f1, f2
     real, dimension(:,:), allocatable :: w, wp, wpold, tlr, bx, by
     integer, dimension(:), allocatable:: ima, iretk, ibc, itc, ibi, ibpoin, itrans, itli
     integer, dimension(:,:), allocatable:: lnd, lbn, icyc, lndold,iae, lnd1, iae1, nserr, ibb
     integer, dimension(:,:), allocatable:: ibp, iaegr, iba
     integer:: ipoin, inpoint,npoinold, nelemold
     character*1  ch
     integer:: nconstr
     real, dimension(:,:,:), allocatable:: constr
     integer, dimension(:), allocatable:: iconstr
     real :: maximal_angle 
     real :: maximal_angle_optimal 
  end type AMAmesh

  type(AMAmesh), allocatable, target :: AMA

end module AMAdata
