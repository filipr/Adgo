module errp_estimates_mod
  use data_mod
  use dual_problem_mod
  use main_data  ! contains type(mesh) :: grid for computation
  use estimates


  implicit none

  public :: ComputeERRpEstimate

contains

  !> main routine for apost error estimate based on the Ritz reconstruction
  !> the reconstruction is the solution of local problems 
  !>  \f$ a(up,phi) = res(phi) \forall phi \in P^{p+1}(K) \f$ 
  !> the final reconstruction is \f$ u^+ = u_h + up\f$ , but we use only \f$ u^+ - u_h = up \f$ 
  !> in the estimates !
  subroutine ComputeERRpEstimate( grid, convfile )
    class( mesh ), intent(inout) :: grid
    character(len=*), intent(in) :: convfile
    class( element ), pointer :: elem
    logical :: flag, loc_implicitly
    integer :: i
    real, dimension(1:ndim) :: norm
    logical :: deg_plus

    call updateProblemPlus( grid, .true. )
    ! we need the larger matrix

    loc_implicitly = state%nlSolver%implicitly
    state%nlSolver%implicitly = .true.

    deg_plus = .false.
    call ComputeSTDGM_Terms( deg_plus )
    ! fill blockPlus

    call CopyBlocksSTtoBlockPlus()

    call updateProblemPlus( grid, .false. )
    ! smaller matrix again
    print*, '!!! ComputeSTDGM_Terms impl=True was computed only for the bigger system and not back for the smaller one!'
    !      call ComputeSTDGM_Terms( )

    ! fill the RHS
    state%nlSolver%implicitly = .false.
    !grid%elem(:)%deg_plus = .true.

    deg_plus = .true.
    call ComputeSTDGM_Terms(deg_plus )
    !grid%elem(:)%deg_plus = .false.
    state%nlSolver%implicitly = loc_implicitly

    !      print*, grid%elem(1)%rhsST(:,:,:)
    !      print*, grid%elem(1)%vec(1,:)

    do i = 1,grid%nelem
       elem => grid%elem(i)
       if ( associated(elem%wSTplus ) ) &
            stop 'wST is already allocated in ComputeERRpEstimate'
       allocate( elem%wSTplus( 1:ndim, 1:elem%dof_plus, 1:elem%Tdof ), source = 0.0 )
    enddo !i

    ! reconstruct the primal solution wST -> wSTplus
    ! fill the local blocks elem%blockPlus
    ! need to compute globally the larger matrix
    ! TODO:  really inefficient needs to be changed in future
    !solve the local problems
    call RitzReconstruction( grid, .true., flag )

    !global estimate
    state%estim(resS,1:ndim) = 0.0

    ! compute |u^+ - uh|_H^1(\om) = |up|_H^1(\om)
    do i = 1,grid%nelem
       elem => grid%elem(i)

       ! compute H1-seminorm
       norm(1:ndim) = evalL8H1STNorm_Elem( elem, ndim, elem%dof_plus, elem%Tdof_plus, &
            elem%wSTplus(1:ndim, 1:elem%dof_plus, 1:elem%Tdof_plus) )
       !         print*, 'semi-normH1 elem(',i,'):', norm(:)

       ! --->>> to elem%eta( ??? )
       state%estim( resS, 1:ndim ) = state%estim( resS, 1:ndim ) + norm(1:ndim)
       elem%eta( resS, 1:ndim) = sqrt(norm(1:ndim))
       elem%estim_loc = sqrt(norm(1))

    enddo !i

  end subroutine ComputeERRpEstimate



end module errp_estimates_mod
