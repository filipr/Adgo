module elemental_estims_mod
   implicit none


   type, public :: Elemental_estims_t
     real :: estimA, estimS, estimT, estimST, estim_loc ! estimates of "residal errors"
     real :: jumpsJh                           ! \sum_{G \in \partial K} |G|^{-1} \int_G [u_h]^2 dS
     real :: rezid, rezid1                     ! reziduum for artificial viscosity
     real :: reg, regT0, regT1, regT2          ! regularity of the solution

     real :: errL2, errH1, errL8,interL8, interLq, interH1
     real, dimension(1:2, 1:3) :: errSnorm	  ! STDGM L2/H1 norms in three different nodes 1-endpoint, 2-node of the time quad, 3-1/2
     real, allocatable, dimension(:,:) :: eta    ! array of different error indicators

   end type Elemental_estims_t

!   type, public, extends ( Elemental_estims_t ) :: ResElemental_estims_t
!
!   end type ResElemental_estims_t

   contains

end module elemental_estims_mod
