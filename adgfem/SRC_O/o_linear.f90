module linear_mod

   implicit none

   type, public :: LinearSol_t
      character(len=20) :: name
      real :: residuum                   ! linear algebra residuum
      real :: consistency                ! linear algebra consistency
      integer :: iter                    ! number of iterations in LA solver
      integer :: iter_tot                ! total number of iterations in LA solver
      integer :: iter_tot_SC             ! total number of iterations in LA solver when Stopping Criteria based on AEE are satisfied
      real :: backward                   ! algebraic error using backward analysis
      real :: tol                         ! tolerance for linear algebra system
      logical :: tol_fixed                ! is tol fixed in *.ini
      integer::  lin_solver_not_conv       ! linear solver did not converges
      logical :: precond_update           ! whether the preconditioner should be updated before solving the linear alg. problem
      logical :: update_matrix            ! whether the matrix C should be updated before solving the linear alg. problem

      contains

      procedure :: init => initLinearSol

   end type LinearSol_t

   type, public, EXTENDS( LinearSol_t ) :: Gmres_t

   contains

      procedure :: init => initGmres

   end type
!FERROR problems with datas in extended data types
   type, public, EXTENDS( LinearSol_t ) :: MultiGrid_t
      logical :: MGLinSolver
      logical :: MGNlnSolver
      integer :: MaxMGlevel                ! maximal level for multigrid method
      integer :: CurMGlevel                ! current level for multigrid method
      integer, dimension(9) :: MGlvlSize
      integer, dimension(9) :: MGlvlDof
      integer, dimension(9) :: MGlvlDeg
      integer :: MGnsize                   ! vector size for MG

   contains

      procedure :: init => initMultigrid

   end type MultiGrid_t


   contains


   subroutine initLinearSol( this, name, tol , mg_type )
      class(LinearSol_t), intent ( inout ) :: this
      character(len=20), intent( in ) :: name
      real, intent( in ) :: tol
      character(len=10), intent( in ), optional :: mg_type

      print*, 'initLinSol should not be called - abstract only'
   end subroutine initLinearSol


   subroutine initGmres( this, name, tol, mg_type )
      class(Gmres_t), intent ( inout ) :: this
      character(len=20), intent( in ) :: name
      real, intent( in ) :: tol
      character(len=10), intent( in ), optional :: mg_type

      this%name = name
      this%tol = tol
      this%precond_update = .true.
      this%update_matrix = .true.

      this%tol_fixed = .true.
      if ( this%tol <= 0.D+00 ) this%tol_fixed = .false.


      if( this%name == "GMRES") then
        write(*,'(a54,es10.4)') &
             '  # GMRES linear iterative solver without prec, tol = ',this%tol

      elseif( this%name == "GMRES_D") then
        write(*,'(a62,es10.4)') &
             '  # GMRES linear iterative solver with block Diag prec, tol = ',this%tol

      elseif( this%name == "GMRES_ILU") then
        write(*,'(a58,es10.4)') &
             '  # GMRES linear iterative solver with ILU(0) prec, tol = ',this%tol

      elseif( this%name == "ILU") then
        write(*,'(a58,es10.4)') &
             '  # SMOOTHING ILU(0) prec, tol = ',this%tol
      else
        print*,'UNKNOWN linear iterative solver -- STOP'
        STOP
      endif

   end subroutine



   subroutine initMultigrid( this, name, tol, mg_type )
      class(MultiGrid_t), intent ( inout ) :: this
      character(len=20), intent( in ) :: name
      real, intent( in ) :: tol
      character(len=10), intent( in ), optional :: mg_type

      if ( .not. present( mg_type ) ) then
         stop 'PROBLEM: mg_type must be specified for MG methods'
      else

         this%tol = tol


         if( mg_type == 'UMFPACK' ) then
            write(*,*) ' # direct solver is used -- UMFPACK'
     !
         elseif( mg_type == 'AGMG' ) then
            write(*,*) ' # algebraic multigrid approach -- AGMG by Yvan Notay'
     !
         elseif( mg_type == 'JACOBI' ) then
            write(*,*), ' # iterative solver is used -- block Jacobi'
     !
         elseif( mg_type == 'GS' ) then
            write(*,*), ' # iterative solver is used -- block Gaus-Seidel'
     !
         elseif( mg_type == 'MG1JACJAC' ) then
            write(*,*), ' # linear multigrid approach -- MG with Jacobi as 1x smoother and exact solver'
     !
         elseif( mg_type == 'MG2JACJAC' ) then
            write(*,*), ' # linear multigrid approach -- MG with Jacobi as 2x smoother and exact solver'
     !
         elseif( mg_type == 'MG3JACJAC' ) then
            write(*,*), ' # linear multigrid approach -- MG with Jacobi as 3x smoother and exact solver'
     !
         elseif( mg_type == 'MG1GSGS' ) then
            write(*,*), ' # linear multigrid approach -- MG with G-S as 1x smoother and exact solver'
     !
         elseif( mg_type == 'MG2GSGS' ) then
            write(*,*), ' # linear multigrid approach -- MG with G-S as 2x smoother and exact solver'
     !
         elseif( mg_type == 'MG3GSGS' ) then
            write(*,*), ' # linear multigrid approach -- MG with G-S as 3x smoother and exact solver'
     !
         elseif( mg_type == 'MG1JACGS' ) then
            write(*,*), ' # linear multigrid approach -- with Jacobi as 1x smoother and G-S as exact solver'
     !
         elseif( mg_type == 'MG2JACGS' ) then
            write(*,*), ' # linear multigrid approach -- with Jacobi as 2x smoother and G-S as exact solver'
     !
         elseif( mg_type == 'MG3JACGS' ) then
            write(*,*), ' # linear multigrid approach -- with Jacobi as 3x smoother and G-S as exact solver'
     !
         elseif( mg_type == 'MG_Jacobi1' ) then
            print*,  'PROBLEM %lin_solver set from mg_type but it should be: MG_Jacobi1'
            write(*,*) ' # linear multigrid approach -- bJacobi smoother and 1xJacobi iter. as exact solution'
     !
         elseif( mg_type == 'MG_JacobiX' ) then
            print*,  'PROBLEM %lin_solver set from mg_type but it should be: MG_JacobiX'
            write(*,*) ' # linear multigrid approach -- bJacobi smoother and 10xJacobi iteration as exact solution'
     !
         elseif( mg_type == 'MGxGMRES' ) then
            print*,  'PROBLEM %lin_solver set from mg_type but it should be: MGxGMRES'
            write(*,*) ' # linear multigrid approach -- bJacobi smoother and GMRES as exact solution'
     !
         elseif( mg_type == 'MG_ILU') then
            print*,  'PROBLEM %lin_solver set from mg_type but it should be: MG_ILU'
            write(*,*) ' # linear multigrid approach -- ILU smoother'
     !
         else
            stop 'unknown multigrid method'
         endif

      endif

      this%tol_fixed = .true.
      if ( this%tol <= 0.D+00 ) this%tol_fixed = .false.

   end subroutine


end module linear_mod
