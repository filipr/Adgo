!> iinput / output subroutines
module io_sub
  use mesh_oper
  use main_data
  use eval_sol
  use model_oper
  use matrix_oper
  use stdgm_mod

  implicit none


  public:: WriteResults
  public:: WriteVecRHS
  public:: OutputDGFEMtri
  public:: OutputDGFEMsol
  public:: OutputDGFEMerrors
  public:: OutputDGFEMestims
  public:: SetCommandWithFileName1
  public:: SetFileNames
  public:: WriteProgressOutput

  public:: OutputElementSol
  public:: OutputElementErrors
  public:: OutputElementEstims
  public:: ElementReadResults

  public:: TEST_ILU
  public:: TEST_ILU1

  public :: WriteDiscribErrOutput
  public :: OutputDiscribErr
  public :: WriteDiscAlgErrOutput
  public :: WriteDiscribErrOutputMatlab
  public :: OutputDiscribErrMatlab

contains


 !> Output of computed results in file'rsolfile', coefficients of the basis
 !> expansion and also its \f$ \Pi_{h0} \f$ projection "results" for
 !> ANGENER code
 subroutine WriteResults(solfile)
   character(len=*), intent(in) :: solfile
   integer :: ifile=12
   integer :: i, k, ist, deg, dof
   real, dimension(1:ndim) :: wl

   open(ifile, file=solfile, status='UNKNOWN')
   write(ifile,*) grid%nelem, ndim, state%nsize, &
        state%time%iter, state%time%ttime, state%time%tau(1), state%err(L2)

   do i=1,grid%nelem
      deg = grid%elem(i)%deg
      dof = grid%elem(i)%dof
      do k=1,ndim
         ist = (k-1)*dof + 1
         write(ifile,*) deg, grid%elem(i)%w(0,ist: ist+dof-1)
      enddo
   enddo
   close(ifile)


   !print*, " output file 'results' for ANGENER"
   open(ifile, file='results', status='UNKNOWN')

   do i=1,grid%nelem
      call Eval_aver_w_Elem(grid%elem(i), wl)
      write(ifile,*) wl(1:ndim)
   enddo

   close(ifile)

 end subroutine WriteResults

 !> Output of elem%vec(rhs, :) into a file
 subroutine WriteVecRHS(solfile)
   character(len=*), intent(in) :: solfile
   integer :: ifile=12
   class(element), pointer :: elem
   integer :: i,k, ist,  dof

   print*,'-------writing file: ',solfile
   open(ifile, file=solfile, status='UNKNOWN')

   do i=1,grid%nelem
      elem => grid%elem(i)
      dof = elem%dof
      if(elem%deg_plus) dof = elem%dof_plus
      do k=1,ndim
         ist = (k-1)*dof + 1
         write(ifile,'(a5,i5,a3,i1,120es12.4)') 'elem=',i,' k=',k,&
              elem%vec(rhs,ist: ist+dof-1)
      enddo
   enddo
   close(ifile)


 end subroutine WriteVecRHS


  !> reading of results from solfile in ADGFEM format per each element
  subroutine ElementReadResults(solfile)
    character(len=*), intent(in) :: solfile
    class(element), pointer :: elem
    real, dimension(:), allocatable :: wi

    integer :: ifile=12
    integer :: i, k, l, deg, dof

    allocate(wi(1:100) )

    open(ifile, file=solfile, status='OLD')
    read(ifile,*) i, k, state%nsize, &
         state%time%iter, state%time%ttime, state%time%tau_old, state%err(L2)

    if(i /= grid%nelem) then
       print*,'Incompatible files "grid" and "sol", different nelem:',grid%nelem,i
       stop
    endif

    if(k /= ndim) then
       print*,'Incompatible files "grid" and "sol", different ndim:',ndim, k
       stop
    endif

    do i=1,grid%nelem
       elem => grid%elem(i)

       do k=1,ndim

          !read(ifile,*) deg, irt, wi(1:(deg+1)*(deg+2)/2)
          read(ifile,*) deg, wi(1:(deg+1)*(deg+2)/2)

          if(k == 1) then
             elem%deg = deg
             call elem%initElementDof( )
             dof = elem%dof

          else
             if(deg /=  elem%deg) then
                print*,'Error in ElementReadResults, different order of ', &
                     ' approximation for different components of w'
             endif
          endif

          if(k == 1) then
             allocate(elem%w(0:state%time%deg+1,1:dof*ndim) )
          endif
          do l=0, state%time%deg+1
             elem%w(l, (k-1)*dof + 1 : k * dof) = wi(1:dof)
          enddo
       enddo
    enddo
    close(ifile)
    deallocate(wi)

  end subroutine ElementReadResults


  !> output of triangulation for DGFEM for visualization into file "tri" in ANGENER format
  !>
  !> triangles are simply copied, quadrilaterals are divided onto 2 elements
  subroutine OutputDGFEMtri(gridfile)
    character(len=*), intent(in) :: gridfile
    integer, parameter :: ifile = 11
    integer:: i,j, il(3), dof, deg

    !print*,'Saving :',gridfile

    ! output  of a triangulation
    open(ifile, file = gridfile, status = 'unknown')

 !   write(ifile, *) grid%npoin, nelem, grid%nbelm, grid%nbc
    write(ifile, *) grid%npoin, grid%n3elem, grid%nbelm, grid%nbc
    write(ifile, '(3x,2(2es14.6,2i5))') grid%xper(1,1:nbDim),grid%iper(1,1:nbDim),&
          grid%xper(2, 1:nbDim),grid%iper(2,1:nbDim)

    do i=1,grid%npoin
       write(ifile,'(3es14.6)') grid%x(i,1:nbDim)
    enddo

    do i=1,grid%nelem
       deg = grid%elem(i)%deg
       dof = grid%elem(i)%dof
       do j=1,grid%elem(i)%type -2
          if(grid%elem(i)%type == 3) then
             il(1:3) = (/1, 2, 3/)

             if(grid%elem(i)%HGnode)&
                  il(1:3) = grid%elem(i)%HGvertex(1:3)

          else if(grid%elem(i)%type == 4) then
             if(j == 1) then
                !!il(1:3) = (/1, 2, 3/)
                !il(1:3) = (/2, 4, 1/)
                il(1:3) = (/1, 2, 4/)

                if(grid%elem(i)%HGnode)then
                   il(1:nbDim) = grid%elem(i)%HGvertex(1:nbDim)
                   il(3) = grid%elem(i)%HGvertex(4)
                endif

             else
                !il(1:3) = (/1, 3, 4/)
                il(1:3) = (/2, 3, 4/)

                if(grid%elem(i)%HGnode)&
                     il(1:3) = grid%elem(i)%HGvertex(2:4)

             endif
          else
             print*,'Elem%type > 4 is not yet implemented in OutputDGFEM'
             stop
          endif

          write(ifile,*) grid%elem(i)%face(idx,il(1:3))
       enddo
    enddo


    do i=1,grid%nbelm
       write(ifile,*)  grid%b_edge(i)%lbn(1:2),grid%b_edge(i)%ibc
    enddo

    close(ifile)

  end subroutine OutputDGFEMtri

  !> draw error per element
  subroutine OutputElementErrors(elem, ifile1, ifile2)
    type(element), intent (inout) :: elem
    integer, intent(in) :: ifile1, ifile2
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:), allocatable :: wi, q
    real, dimension(:,:), allocatable :: xi, Fx, qExact
    integer :: Qnum, deg, dof, k, j, l, ist, Qdof, Qlen,RGlevel

    deg = elem%deg
    dof = elem%dof
    Qnum = max(1, elem%deg)
    RGlevel = 0

    if(elem%type == 4) Qnum = Qnum + QnumOffset

    Qdof = state%space%L_rule(Qnum)%Qdof
    Qlen = Qdof
    if(elem%type == 4) Qlen = Qlen/2  ! quandrilaterals forms two triangles

    allocate(wi(1:ndim, 1:dof), q(1:Qdof, 1:ndim))
    allocate(xi(1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim), qExact(1:Qdof, 1:ndim) )

    phi => state%space%L_rule(Qnum)%phi(1:dof, 1:Qdof)

    ! evaluation of the exact solution in Lag. integ. nodes
    xi(1:Qdof, 1:nbDim) = state%space%L_rule(Qnum)%lambda(1:Qdof,1:nbDim)   !
    call ComputeF(elem, Qdof, xi(1:Qdof, 1:nbDim), Fx(1:Qdof, 1:nbDim) )

    call Exact_Sol(Qdof, Fx(1:Qdof, 1:nbDim), qExact(1:Qdof, 1:ndim), state%time%ctime)
    !do j=1, Qdof
    !   call Exact_Scalar(Fx(j,1:nbDim), qExact(1:ndim, j), state%time%ctime)
    !   !call Der_Exact_Scalar(Fx(j,1:nbDim), DwExact(j, 1:ndim, 1:nbDim), state%time%ctime)
    !enddo

    do k=1,ndim
       ist = (k-1)*dof + 1
       wi(k, 1:dof) = elem%w(0,ist:ist+dof-1)

       ! evaluation of w - wExact in the Langrangian nodes
       do j=1, Qdof
          q(j, k) = qExact(j,k) - dot_product( wi(k, 1:dof), phi(1:dof,j) )
       enddo
    enddo

    ! output to file
    do l = 0, elem%type - 3
       do j=2, 2  !!!!!2   ! j = 1 error,  j==2 abs(error)

          do k=1,ndim
             ! correction
             !       if (deg > 3) then
             !          write(ifile,'(i4,100es14.6)') 3, q(k,l*Qlen + 1: (l+1)*Qlen)
             !       else

             if( j== 1) then
                write(ifile1,'(2i5,100es14.6)')deg, RGlevel, qExact(l*Qlen + 1: (l+1)*Qlen, k)
                write(ifile2,'(2i5,100es14.6)')deg, RGlevel, q(l*Qlen + 1: (l+1)*Qlen, k)
             else
                write(ifile1,'(2i5,100es14.6)')deg, RGlevel, abs(qExact(l*Qlen + 1: (l+1)*Qlen, k))

                write(ifile2,'(2i5,100es14.6)')deg, RGlevel, abs(q(l*Qlen + 1: (l+1)*Qlen, k))
             endif
!       endif
!!!!!!!!!!!!!!!!!!!!!!!
          enddo
       enddo
    enddo

    deallocate(wi, q)

  end subroutine OutputElementErrors




  !> called by subroutine "OutputDGFEMsol" per each element \f$K\in {\cal T}_h \f$
  subroutine OutputElementSol(elem, ifile)
    type(element), intent (inout) :: elem
    integer, intent(in) :: ifile
    real, dimension(:,:), pointer :: phi
    real, dimension(:,:), allocatable :: wi, q
    integer :: Qnum, deg, dof, k, l, ist, Qdof, Qlen, RGlevel

    deg = elem%deg
    dof = elem%dof
    Qnum = max(1, elem%deg)
    Qnum = elem%deg
    RGlevel = 0

! correction
!    Qnum = min(3, elem%deg)
!!!!!!!!!!!!!!!!!!!!!!!!


    if(elem%type == 4) Qnum = Qnum + QnumOffset

    Qdof = state%space%L_rule(Qnum)%Qdof
    Qlen = Qdof
    if(elem%type == 4) Qlen = Qlen/2  ! quandrilaterals forms two triangles

    allocate(wi(1:ndim, 1:dof), q(1:ndim, 1:Qdof) )

    phi => state%space%L_rule(Qnum)%phi(1:dof, 1:Qdof)

    do k=1,ndim
       ist = (k-1)*dof + 1
       wi(k, 1:dof) = elem%w(0,ist:ist+dof-1)

       ! evaluation of w in the Langrangian nodes
       do l=1, Qdof
          q(k, l) = dot_product( wi(k, 1:dof), phi(1:dof,l) )
       enddo

       !if(elem%i <= 20) then
       !   write(*,'(a4,80es12.4)') 'wi.', q(1, 1:Qdof),wi(k,1:dof)
       !endif

    enddo

    ! output to file
    do l = 0, elem%type - 3
       do k=1,ndim
          !write(ifile,'(i,100es15.6') max(1,deg), q(k,l*Qlen + 1: (l+1)*Qlen)
!          write(ifile,'(i4,100es15.6)') deg, q(k,l*Qlen + 1: (l+1)*Qlen)
!          write(ifile,'(i4,100es15.6)') Qnum, q(k,l*Qlen + 1: (l+1)*Qlen)
! correction
!       if (deg > 3) then
!          write(ifile,'(i4,100es15.6)') 3, q(k,l*Qlen + 1: (l+1)*Qlen)
!       else
!          write(ifile,'(2i5,100es16.8)') deg,  elem%RGlevel, q(k,l*Qlen + 1: (l+1)*Qlen)
          write(ifile,*) deg,  RGlevel, q(k,l*Qlen + 1: (l+1)*Qlen), elem%i, k
!       endif
!!!!!!!!!!!!!!!!!!!!!!!
       enddo
    enddo

    deallocate(wi, q)

  end subroutine OutputElementSol


   !> draw estimator per element
  subroutine OutputElementEstims(ndimL, elem, ifile)
    integer, intent(in) :: ndimL   ! number of quantities at output
    type(element), intent (inout) :: elem
    integer, intent(in) :: ifile
    real, dimension(:), allocatable :: q
    integer :: deg, l, k, pK, is, iss, RGlevel
    real :: loc_tol

    RGlevel = 0

    loc_tol = state%space%adapt%tol_min /(1.05 * grid%nelem)**0.5
    allocate (q(1:ndimL) )
    q(:) = 0.

    if (state%space%adapt%adapt_method /= 'ALG' .and. state%space%adapt%adapt_method /= 'ALG2') then

       if( (state%modelName == 'scalar' .or.state%modelName == '2eqs')  &
            .and. state%space%estim_space == 'pNeu' ) then
          is =  P_F_p1 - 1  ! = 9
          iss = P_potP      ! = 18

          q(1:iss) = elem%eta(1:iss, 1)


          do l=1,3
             if(q(is+l) > 0. .and. q(is +1 +l) > 0.) then
                !q(iss +l ) = log(q(11+l) / q(is+l)) / log (2.)
                q(iss +l ) = q(is+l) / q(is+1+l) / elem%diam**0.5
             else
                q(iss+l) = 0.
             endif

             if(q(is++l) > 0. .and. q(is+5+l) > 0.) then
                !q(iss+3+l ) = log(q(14+l) / q(13+l)) / log (2.)
                q(iss+3+l ) = q(is+4+l) / q(is+5+l) / elem%diam**0.5
             else
                q(iss+3+l) = 0.
             endif
             !if(elem%i == 1) print*,'#ESW@#$',iss +l, is+l, is+1+l,'|',iss+3+l, is+4+l,  is+5+l
          enddo
          pK = max(1,elem%deg)
          !q(3) = elem%rezid / elem%area / elem%diam**(2*pK)* elem%diam**5
          !q(3) = elem%rezid / elem%area / elem%diam**(2*pK-1)
          q(7) = elem%errH1
          !q(8) = elem%eta(P_pot,1) / elem%eta(P_potP,1)
          q(3) = elem%eta(P_tot,1) / loc_tol  !state%space%adapt%tol_min * grid%nelem**0.5!  / ( state%space%adapt%tol_min * (elem%area/state%space%domain_volume)**0.5 )
          !q(6) = elem%eta(P_potP,1)/loc_tol ! (state%space%adapt%tol_min /(1.05 * grid%nelem)**0.5 ! state%space%adapt%tol_min * grid%nelem**0.5


          q(8:9) = 2.
          !if(elem%deg >= 2) q(8)  = elem%eta(P_F_p1, 1) / elem%eta(P_F_p2, 1) ! /  elem%diam**0.5
          q(8) = elem%eta(P_tot,1) ! * (elem%area/state%space%domain_volume)**0.5 )
          q(6) = elem%reg / elem%regT2
          q(9) = elem%reg / elem%regT0
          !q(10) = elem%eta(P_tot,1) / elem%diam**elem%deg / elem%area**0.5
          q(11) = elem%eta(P_tot,1) / elem%errH1
          !q(12) = elem%errH1 / elem%diam**elem%deg / elem%area**0.5
          q(10) = elem%eta(P_tot,1)
          q(12) = elem%eta(P_potP,1)

          ! for Honza Papez
          !if(elem%i ==1) write(88,*) ' numb     error        estim          souradnice teziste'
          !write(88, '(i5, 4es16.8)' ) elem%i, elem%errH1, q(8), elem%xc(:)

          !!q(3) = elem%rezid / elem%area / elem%diam**(2*pK-3) / elem%diam**2

          !write(*,'(a6,i5,40es10.2)') '###',elem%i, q(3),  elem%rezid, elem%area, elem%diam**(2*pK-5), elem%diam, elem%reg, elem%reg1

       elseif( (state%modelName == 'scalar' .or.state%modelName == '2eqs') &
            .and. state%space%estim_space == 'HO_rec') then ! scalar equation
          iss = HO_trunc_H1_p2 ! = 12
          q(1:iss) = (abs( elem%eta(1:iss, 1)) )**0.5

          q(3) = elem%estimA
          q(11) = elem%estim_loc
          q(12) = elem%psplit
          q(13) =  elem%errL2
          q(14) =  elem%errH1
          q(15) = (abs(elem%eta(HO_estim_L2_p2, 1))) **0.5 / elem%errL2
          q(16) = (abs(elem%eta(HO_estim_H1_p2, 1)))**0.5 / elem%errH1
          !q(17) = elem%eta(HO_estim_H1_p1, 1)**0.5 / elem%eta(HO_estim_H1_p2, 1)**0.5
          q(17) = elem%estimST/ (state%space%adapt%tol_min /(1.05 * grid%nelem)**0.5)
          !                           !state%tol_min*grid%nelem**0.5
          !q(18) = elem%eta(HO_estim_H1_p1, 1) / elem%eta(HO_estim_H1_p2, 1) / (0.75 * elem%diam**0.5 )
          !q(18) = elem%reg /elem%regT0
          !q(19) = elem%reg /elem%regT2
          k = HO_estim_H1_p2
          l = HO_estim_H1_p0
          if(elem%eta(l, 1) > 0 .and. elem%eta(k, 1)  > 0. ) &
               q(18) = (elem%eta(l, 1) / elem%eta(k, 1) )**(1.0/(elem%deg - 1)) * DOFtriang(elem%deg -1)

          l = HO_estim_H1_p1
          if(elem%eta(l, 1) > 0 .and. elem%eta(k, 1)  > 0. ) &
               q(19) = (elem%eta(l, 1) / elem%eta(k, 1) )**(1.0/(elem%deg - 0)) * DOFtriang(elem%deg -0)

          l = HO_estim_H1_p2
          if(elem%eta(l, 1) > 0 .and. elem%eta(k, 1)  > 0. ) &
               q(20) = (elem%eta(l, 1) / elem%eta(k, 1) )**(1.0/(elem%deg + 1)) * DOFtriang(elem%deg +1)


       elseif( (state%modelName == 'scalar' .or.state%modelName == '2eqs') &
            .and. state%space%estim_space /= 'RES' .and. state%space%adapt%adapt_method /= 'ANI' &
            .and. state%space%adapt%adapt_method /= 'Ahp' ) then ! scalar equation
          q(1) = elem%estimST            ! value from estimDual.f90 ( already square root )
          q(2) = elem%eta(eN1,1)**0.5   ! values from errorFlux.f90
          q(3) = elem%eta(eN2,1)**0.5
          q(4) = elem%eta(NC1n,1)**0.5
          q(5) = elem%eta(DFnS,1)**0.5
          q(6) = elem%eta(DFnT,1)**0.5
       else
          q(1) = elem%estimA
          q(2) = elem%estimS
          q(3) = elem%estimT
          q(4) = elem%estimST
          q(5) = elem%estim_loc   !elem%estimS / elem%estimT
          q(6) = elem%errL8
          q(7) = elem%interL8
          q(8) = elem%errL2
          q(9) = abs(elem%interLq)**0.5
       endif
    else
       q(1) = elem%errTot
       q(2) = elem%errDisc
       q(3) = elem%errAlg
       q(4) = elem%estTot
       q(5) = elem%estDisc
       q(6) = elem%estAlg

    endif



    !print*,'###',elem%i, q(:)
    ! one value per element, hence the minimal aceptable degree is 1
    deg = 1

    ! output to file
    do l = 0, elem%type - 3
       do k=1,ndimL
          write(ifile,'(2i5,100es14.6)')deg, RGlevel, q(k), q(k), q(k)
       enddo
    enddo

    deallocate( q)

  end subroutine OutputElementEstims



  !> output of results for DGFEM visualization in file "sol(num)", (num)=isol
  !>
  !> format of file = format of "dgfem" package,
  !> solution of \f$ {\bf w}\in{\bf S}_{hp} \f$ in Lagrangian nodes
  subroutine OutputDGFEMsol(solfile)
    character(len=*), intent(in) :: solfile
    integer  :: ifile1 = 12
    integer:: i, nelem

    nelem = 0
    do i=1,grid%nelem
       nelem = nelem + grid%elem(i)%type -2
    enddo

    ! output of results
    open(ifile1, file = solfile, status = 'unknown')
    !write(ifile1, *) nelem, ndim
    write(ifile1, *) nelem, ndim, state%time%ttime, state%time%tau(1), state%time%iter ! new

    do i=1,grid%nelem
       call OutputElementSol(grid%elem(i), ifile1)
    enddo
!FERROR commented? division by 0
!
!    write(ifile1,'(5es12.4, a30)') state%rho_infty,  state%v_infty, state%p_infty, &
!           state%alpha_infty, state%theta_infty,'  (rho_8, v_8, p_8, a_8, T_8);'
!    write(ifile1,'(2es12.4,a15)') state%BC(1:2)%press_extrap, ' (p_extrap);'
!    write(ifile1,'(2es12.4,a15,2es12.4, a15, es12.4, a15)') &
!         state%BC(1:2)%rho_extrap, ' (rho_extrap);', &
!         state%BC(1:2)%press_extrap/state%p_infty, ' (p_extrap/p8);', &
!         state%BC(2)%press_extrap/state%BC(1)%press_extrap, ' (p_out/p_in);'
!    write(ifile1, '(4es12.4, a8)') state%BC(1)%w(1:4), '  BC_in'
!    write(ifile1, '(4es12.4, a8)') state%BC(2)%w(1:4), ' BC_out'
!    write(ifile1, 99 ) state%BC(1)%w(1:4), state%BC(2)%w(1:4), &
!         state%BC(1:2)%press_extrap, &
!         state%BC(1:2)%rho_extrap, state%BC(1:2)%press_extrap/state%p_infty, &
!         state%BC(2)%press_extrap/state%BC(1)%press_extrap
!    close(ifile1)
!
!99  format ( 4es12.4, ' |', 4es12.4 ,' ||', 2es12.4, ' |', 2es12.4, ' .', 2es12.4, ' ||', es12.4)


    if(state%modelName == 'NSe') then
       if(state%numBC > 1) then
          write(ifile1,'(5es14.6, a30)') state%rho_infty,  state%v_infty, state%p_infty, &
               state%alpha_infty, state%theta_infty,'  (rho_8, v_8, p_8, a_8, T_8);'
          write(ifile1,'(2es12.4,a15)') state%BC(1:2)%press_extrap, ' (p_extrap);'
          write(ifile1,'(2es12.4,a15,2es12.4, a15, es12.4, a15)') &
               state%BC(1:2)%rho_extrap, ' (rho_extrap);', &
               state%BC(1:2)%press_extrap/state%p_infty, ' (p_extrap/p8);', &
               state%BC(2)%press_extrap/state%BC(1)%press_extrap, ' (p_out/p_in);'
          write(ifile1, '(4es12.4, a8)') state%BC(1)%w(1:4), '  BC_in'
          write(ifile1, '(4es12.4, a8)') state%BC(2)%w(1:4), ' BC_out'
          write(ifile1, 99 ) state%BC(1)%w(1:4), state%BC(2)%w(1:4), &
               state%BC(1:2)%press_extrap, &
               state%BC(1:2)%rho_extrap, state%BC(1:2)%press_extrap/state%p_infty, &
               state%BC(2)%press_extrap/state%BC(1)%press_extrap
       else
          write(ifile1,'(5es14.6, a30)') state%rho_infty,  state%v_infty, state%p_infty, &
               state%alpha_infty, state%theta_infty,'  (rho_8, v_8, p_8, a_8, T_8);'
          write(ifile1, '(4es14.6, a8)') state%BC(1)%w(1:4), '  BC_FF'
       endif

    endif

    close(ifile1)

99  format ( 4es12.4, ' |', 4es12.4 ,' ||', 2es12.4, ' |', 2es12.4, ' .', 2es12.4, ' ||', es12.4)

  end subroutine OutputDGFEMsol



  !> output of error for DGFEM visualization in file "sol(num)", (num)=isol
  !>
  !> format of file = format of "dgfem" package,
  !> exact solution of \f$ {\bf w}_h\in{\bf S}_{hp} \f$ in Lagrangian nodes
  !> computational error \f$ {\bf w} - {\bf w}_h \f$ in Lagrangian nodes
  subroutine OutputDGFEMerrors(exact_file, error_file)
    character(len=*), intent(in) :: exact_file, error_file
    integer  :: ifile1 = 12, ifile2 = 13
    integer:: i, nelem

    !print*,'Saving :',exact_file
    !print*,'Saving :',error_file

    nelem = 0
    do i=1,grid%nelem
       nelem = nelem + grid%elem(i)%type -2
    enddo

    ! output of results
    open(ifile1, file = exact_file, status = 'unknown')
    open(ifile2, file = error_file, status = 'unknown')

    !write(ifile1, *) nelem, ndim
    write(ifile1, *) nelem, ndim, state%time%ttime, state%time%tau(1), state%time%iter ! new
    write(ifile2, *) nelem, ndim, state%time%ttime, state%time%tau(1), state%time%iter ! new

    do i=1,grid%nelem
       call OutputElementErrors(grid%elem(i), ifile1, ifile2)
    enddo

    write(ifile1,*) state%rho_infty,  state%v_infty, state%p_infty, &
           state%alpha_infty, state%theta_infty
    write(ifile2,*) state%rho_infty,  state%v_infty, state%p_infty, &
           state%alpha_infty, state%theta_infty

    close(ifile1)
    close(ifile2)

  end subroutine OutputDGFEMerrors


  !> output of estimates and errors for DGFEM visualization in file "sol(num)", (num)=isol
  !> format of file = format of "dgfem" package,
  !> values computed in estimDual.f90 and  errorFlux.f90
  subroutine OutputDGFEMestims(estim_file)
    character(len=*), intent(in) :: estim_file
    integer  :: ifile1 = 12
    integer:: i, nelem, ndimL

    if(state%space%estim_space == 'pNeu' ) then
       ndimL =  24  !max_eta + 8 ?

    elseif(state%space%estim_space == 'HO_rec' ) then
       ndimL =  20  !max_eta

    elseif (state%space%adapt%adapt_method == 'ALG' .or. state%space%adapt%adapt_method == 'ALG2') then
       ndimL = 6  !we write at output ndimL quantities, defined in OutputElementEstims
    else
       ndimL = 9  !we write at output ndimL quantities, defined in OutputElementEstims
    endif

    !print*,'Saving :',estim_file

    nelem = 0
    do i=1,grid%nelem
       nelem = nelem + grid%elem(i)%type -2
    enddo

    ! output of results
    open(ifile1, file = estim_file, status = 'unknown')

    !write(ifile1, *) nelem, ndim
    write(ifile1, *) nelem, ndimL, state%time%ttime, state%time%tau(1), state%time%iter ! new

    do i=1,grid%nelem
       call OutputElementEstims(ndimL, grid%elem(i), ifile1)
    enddo

    write(ifile1,*) state%rho_infty,  state%v_infty, state%p_infty, &
           state%alpha_infty, state%theta_infty

    close(ifile1)

  end subroutine OutputDGFEMestims


  !> Open file (mesh, tri, sol) with a succesive name
  subroutine SetCommandWithFileName1(command_name, T, inumber) !, ifile)
    character(len=30),intent(inout) :: command_name
    character(len=1), intent(in) :: T
    integer, intent(in) :: inumber!, ifile
    character(len=5) :: ch5
    integer :: num_size, text_size, inum
    integer :: is

    text_size = 12
    num_size = 5

    inum = inumber
    if(inumber < 0) inum = 10**num_size + inumber

    !print*,'????',inum, inumber

    if(T == 'T') then
       command_name = 'mv tri0 triA00000    '
    elseif(T == 'S') then
       command_name = 'mv sol0 solA00000    '
    else
       print*,'Bad input value "T" in SetCommandWithFileName in problem.f90'
       stop
    endif

    if(inum > 0) then
       is = int(log(1.*inum)/log(10.))
    else
       is = 0
    endif

    !print*,'!!!',inum,is, num_size+text_size-is, num_size+text_size, num_size-is, num_size

    write( ch5, '(i5)' ) inum  ! change the format if num_size /= 5 !!!
    command_name(num_size+text_size-is:num_size+text_size) = ch5(num_size-is: num_size)

    !print*,'######',command_name(1:num_size+text_size),'|',inum, ch5, is

  end subroutine SetCommandWithFileName1

  !> Setting of names for tri* sol*
  !> if command = .true. then setting of the commmand for plotting
  subroutine SetFileNames(command, command_name, tri_name, sol_name, exa_name, err_name, &
       est_name, Aname)
    logical, intent(in) :: command, Aname
    character(len=50),intent(inout) :: command_name, tri_name, sol_name, exa_name, err_name
    character(len=50),intent(inout) :: est_name
    character(len=5) :: ch5
    integer :: num_size, text_size, long_text_size, file_size, inum
    integer :: is

    long_text_size = 10
    text_size = 4
    num_size = 5
    file_size = text_size + num_size

    !print*,'###',state%space%adapt%max_adapt_level,  state%time%OutTime

    !if(state%space%adapt%max_adapt_level == 0 .or. state%time%OutTime > 0.) then
    if(.not. Aname) then
       inum = state%isol
       tri_name = 'tri-00000    '
       sol_name = 'sol-00000    '
       exa_name = 'exa-00000    '
       err_name = 'err-00000    '
       est_name = 'est-00000    '
    else
       inum = state%space%adapt%adapt_level
       if(state%space%adapt%adapt_level < 0) inum = 10**num_size + state%space%adapt%adapt_level

       tri_name = 'triA00000    '
       sol_name = 'solA00000    '
       exa_name = 'exaA00000    '
       err_name = 'errA00000    '
       est_name = 'estA00000    '
    endif

    if(inum > 0) then
       is = int(log(1.*inum)/log(10.))
    else
       is = 0
    endif

    !print*,'!!!',inum,is, num_size+text_size-is, num_size+text_size, num_size-is, num_size

    write( ch5, '(i5)' ) inum  ! change the format if num_size /= 5 !!!
    tri_name(num_size+text_size-is:num_size+text_size) = ch5(num_size-is: num_size)
    sol_name(num_size+text_size-is:num_size+text_size) = ch5(num_size-is: num_size)
    exa_name(num_size+text_size-is:num_size+text_size) = ch5(num_size-is: num_size)
    err_name(num_size+text_size-is:num_size+text_size) = ch5(num_size-is: num_size)
    est_name(num_size+text_size-is:num_size+text_size) = ch5(num_size-is: num_size)

    if(command) then
       command_name(1:long_text_size+2+2*file_size) &
            = '../Plotdgm                                    '

       command_name(long_text_size+2:long_text_size+file_size+1) = &
            tri_name(1:file_size)

       command_name(long_text_size+file_size + 3:long_text_size+2*file_size+2) = &
            sol_name(1:file_size)

       command_name(long_text_size+2*file_size+3:long_text_size+2*file_size+10) = &
            ' > smaz'
    endif
  end subroutine SetFileNames

  !> write the files 'triAxxxxx', 'solAxxxxx',
  !> interpolation at t = state%time%ttime-tau(1) + ratio * tau(1)
  !> tri-file : output_type = T, TS, ST
  !> sol-file : output type = S, TS, ST
  subroutine WriteProgressOutput_time( output_type, ratio )
    character(len=*), intent(in) :: output_type
    real, intent(in):: ratio
    class(element), pointer :: elem
    real, dimension(:), allocatable :: Tphi
    integer :: i, k, ndof, Tdof, dof

    ! storing of array elem%w(0,:)
    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof = elem%dof * ndim
       allocate (elem%wS(0:0, 1:ndof) )
       elem%wS(0,1:ndof)  = elem%w(0,1:ndof)
    enddo

    ! interpolation of the solution at time
    if (state%time%disc_time /= 'STDG') then
       ! BDF
       do i=1,grid%nelem
          elem => grid%elem(i)
          ndof = elem%dof * ndim
          elem%w(0,1:ndof) =  elem%w(1,1:ndof)*(1-ratio) + elem%w(0,1:ndof)*ratio
       enddo
    else
       ! STDGM
       do i=1,grid%nelem
          elem => grid%elem(i)
          dof = elem%dof
          Tdof = elem%Tdof
          allocate(Tphi(1:Tdof))

           !TODO FR: Trule basis is evaluated by EvalLegendrePolynomials now
          !call EvalTrulePhi(ratio, Tdof, Tphi(1:Tdof))


          !!TEST
          !do k=0, 10
          !   call EvalTrulePhi(1.*k/10, Tdof, Tphi(1:Tdof))
          !   !write(299,'(120es10.2)') 1.*k,Tphi(1:Tdof)
          !enddo
          !stop

          do k = 1, ndim
             !write(*,'(120es10.2)') elem%wST(k,1:dof,1:Tdof)
             !write(*,'(120es10.2)') Tphi(1:Tdof)

             elem%w(0,dof*(k-1) + 1 :dof*k) = matmul(elem%wST(k,1:dof,1:Tdof), Tphi(1:Tdof) )
          enddo

          deallocate(Tphi)
       enddo
    endif


    ! writting of the output
    call WriteProgressOutput( output_type, .false. )


    ! restoring of array elem%w(0,:)
    do i=1,grid%nelem
       elem => grid%elem(i)
       ndof = elem%dof * ndim
       elem%w(0,1:ndof)  = elem%wS(0,1:ndof)
       deallocate (elem%wS)
    enddo

  end subroutine WriteProgressOutput_time



  !> write the files 'triAxxxxx', 'solAxxxxx'
  !> tri-file : output_type = T, TS, ST
  !> sol-file : output type = S, TS, ST
  subroutine WriteProgressOutput( output_type, Aname )
    character(len=*), intent(in) :: output_type
    logical, intent(in) ::  Aname
    character(len=50) :: command_name, tri_name, sol_name, exa_name, err_name, est_name
    character(len=120) :: output_line
    character(len=1) :: ch1
    integer:: i, lol

    call SetFileNames(.false., command_name, tri_name, sol_name, exa_name, err_name, &
         est_name, Aname)

    do lol=1,len(output_line)
       output_line(lol:lol) = ' '
    enddo

    lol = 10
    output_line(1:lol) = ' Saving: '


    if (state%time%disc_time == 'STDG') then
       do i = 1, grid%nelem
          call Transfer_wST_to_w_Elem(grid%elem(i), 0, grid%elem(i)%TQnum)
       enddo
    endif

    do i=1, len_trim(output_type)
       ch1(1:1) = output_type(i:i)

       if(ch1 == 'T') call OutputDGFEMtri(tri_name)
       if(ch1 == 'S') call OutputDGFEMsol(sol_name)
       if(ch1 == 'E') call OutputDGFEMerrors(exa_name, err_name)
       if(ch1 == 'A') call OutputDGFEMestims(est_name)

       if(ch1 == 'T') output_line(lol+1: lol+9) = tri_name
       if(ch1 == 'S') output_line(lol+1: lol+9) = sol_name
       if(ch1 == 'E') output_line(lol+1: lol+9) = exa_name
       if(ch1 == 'A') output_line(lol+1: lol+9) = est_name

       lol = lol + 11

       if(ch1 == 'E') then
          output_line(lol+1: lol+9) = err_name
          lol = lol+11
       endif

    enddo

    write(*,*) output_line(1:lol)

    ! if(output_type == 'T' .or. output_type == 'TS' .or. output_type == 'ST' .or. &
    !      output_type == 'TE' .or. output_type == 'TSE' .or. output_type == 'STE' ) &
    !      then
    !    call OutputDGFEMtri(tri_name)

    !    open(44, file='commands', status='UNKNOWN', position = 'append')
    !    write(44,'(a4,a10, a6)') 'cp ',tri_name,' tri'
    !    write(44,'(a4,a10, a6)') 'cp ',sol_name,' sol'
    !    close(44)
    ! endif

    ! if(output_type == 'S' .or. output_type == 'TS' .or. output_type == 'ST' .or. &
    !      output_type == 'SE' .or. output_type == 'TSE' .or. output_type == 'STE') &
    !      call OutputDGFEMsol(sol_name)

    ! if(output_type == 'E' .or. output_type == 'TE' .or. output_type == 'TSE' .or. &
    !      output_type == 'STE'.or. output_type == 'SE' ) &
    !      call OutputDGFEMerrors(exa_name, err_name)

    ! if(output_type == 'A') call OutputDGFEMestims(est_name)

  end subroutine WriteProgressOutput



  !> write the files 'Tot-xxxxxx.bb', 'Disc-xxxxxx.bb', 'Alg-xxxxxx.bb',
  !> 'EstTot-xxxxxx.bb', 'EstDisc-xxxxxx.bb', and 'EstAlg-xxxxxx.bb' : input for Medit software
  !> Tot-file : output_type = T
  !> Disc-file : output type = D
  !> Alg-file : output type = A
  !> EstTot-file : output_type = EstT
  !> EstDisc-file : output_type = EstD
  !> EstAlg-file : output_type = EstA
  subroutine WriteDiscribErrOutput( output_type )
    character(len=*), intent(in) :: output_type
    character(len=50) :: err_name
    character(len=3) :: ch3
    integer :: text_size, is, num_size


    if(state%space%adapt%adapt_level > 0) then
        is = int(log(1.*state%space%adapt%adapt_level)/log(10.))
    else
        is = 0
    endif

    num_size = 3
    write( ch3, '(i3)' ) state%space%adapt%adapt_level

    if(output_type == 'T') then
         err_name = 'Tot-000    '
         text_size = 4
         err_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         err_name(num_size+text_size+2:num_size+text_size+4) = '.bb'
    elseif (output_type == 'D') then
         err_name = 'Disc-000    '
         text_size = 5
         err_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         err_name(num_size+text_size+2:num_size+text_size+4) = '.bb'
    elseif (output_type == 'A') then
         err_name = 'Alg-000    '
         text_size = 4
         err_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         err_name(num_size+text_size+2:num_size+text_size+4) = '.bb'
    elseif (output_type == 'EstT') then
         err_name = 'EstTot-000    '
         text_size = 7
         err_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         err_name(num_size+text_size+2:num_size+text_size+4) = '.bb'
    elseif (output_type == 'EstD') then
         err_name = 'EstDisc-000    '
         text_size = 8
         err_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         err_name(num_size+text_size+2:num_size+text_size+4) = '.bb'
    elseif (output_type == 'EstA') then
         err_name = 'EstAlg-000    '
         text_size = 7
         err_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         err_name(num_size+text_size+2:num_size+text_size+4) = '.bb'
    else
      Print*, 'Unknown output_type. Possible choices: T, D, A, EstT, EstD, EstA'

    endif



    call OutputDiscribErr(err_name, output_type)


  end subroutine WriteDiscribErrOutput



  subroutine OutputDiscribErr(err_name, output_type)
    character(len=50) :: err_name
    character(len=*), intent(in) :: output_type
    integer :: i


    open(20+state%time%iter, file=err_name, status='UNKNOWN', position = 'append')
    write(20+state%time%iter,'(a3, i6, a2)') '2 1', grid%nelem, '1'
    do i=1, grid%nelem
       if(output_type == 'T') then
         write(20+state%time%iter,'(es14.6)') grid%elem(i)%errTot
       elseif (output_type == 'D') then
         write(20+state%time%iter,'(es14.6)') grid%elem(i)%errDisc
       elseif (output_type == 'A') then
         write(20+state%time%iter,'(es14.6)') grid%elem(i)%errAlg
       elseif (output_type == 'EstT') then
         write(20+state%time%iter,'(es14.6)') grid%elem(i)%estTot
       elseif (output_type == 'EstD') then
         write(20+state%time%iter,'(es14.6)') grid%elem(i)%estDisc
       elseif (output_type == 'EstA') then
         write(20+state%time%iter,'(es14.6)') grid%elem(i)%estAlg
       endif
    enddo
    close(20+state%time%iter)

  end subroutine OutputDiscribErr



  subroutine WriteDiscribErrOutputMatlab( output_type )
    character(len=*), intent(in) :: output_type
    character(len=50) :: err_name
    character(len=3) :: ch3
    integer :: text_size, is, num_size


    if(state%space%adapt%adapt_level > 0) then
        is = int(log(1.*state%space%adapt%adapt_level)/log(10.))
    else
        is = 0
    endif

    num_size = 3
    write( ch3, '(i3)' ) state%space%adapt%adapt_level

    if(output_type == 'T') then
         err_name = 'Tot-000    '
         text_size = 4
         err_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         err_name(num_size+text_size+2:num_size+text_size+5) = '.txt'
    elseif (output_type == 'D') then
         err_name = 'Disc-000    '
         text_size = 5
         err_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         err_name(num_size+text_size+2:num_size+text_size+5) = '.txt'
    elseif (output_type == 'A') then
         err_name = 'Alg-000    '
         text_size = 4
         err_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         err_name(num_size+text_size+2:num_size+text_size+5) = '.txt'
    elseif (output_type == 'EstT') then
         err_name = 'EstTot-000    '
         text_size = 7
         err_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         err_name(num_size+text_size+2:num_size+text_size+5) = '.txt'
    elseif (output_type == 'EstD') then
         err_name = 'EstDisc-000    '
         text_size = 8
         err_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         err_name(num_size+text_size+2:num_size+text_size+5) = '.txt'
    elseif (output_type == 'EstA') then
         err_name = 'EstAlg-000    '
         text_size = 7
         err_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
         if (stop_crit == 'L') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'L'
         elseif (stop_crit == 'G') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'G'
         elseif (stop_crit == 'N') then
            err_name(num_size+text_size+1:num_size+text_size+2) = 'N'
         endif
         err_name(num_size+text_size+2:num_size+text_size+5) = '.txt'
    else
      Print*, 'Unknown output_type. Possible choices: T, D, A, EstT, EstD, EstA'

    endif



    call OutputDiscribErrMatlab(err_name, output_type)



  end subroutine WriteDiscribErrOutputMatlab



  subroutine OutputDiscribErrMatlab(err_name, output_type)
    character(len=50) :: err_name
    character(len=*), intent(in) :: output_type
    integer :: i


    open(20+state%time%iter, file=err_name, status='UNKNOWN', position = 'append')

    do i=1, grid%nelem
       if(output_type == 'T') then
         write(20+state%time%iter,'(es14.6)') grid%elem(i)%errTot
       elseif (output_type == 'D') then
         write(20+state%time%iter,'(es14.6)') grid%elem(i)%errDisc
       elseif (output_type == 'A') then
         write(20+state%time%iter,'(es14.6)') grid%elem(i)%errAlg
       elseif (output_type == 'EstT') then
         write(20+state%time%iter,'(es14.6)') grid%elem(i)%estTot
       elseif (output_type == 'EstD') then
         write(20+state%time%iter,'(es14.6)') grid%elem(i)%estDisc
       elseif (output_type == 'EstA') then
         write(20+state%time%iter,'(es14.6)') grid%elem(i)%estAlg
       endif
    enddo
    close(20+state%time%iter)



  end subroutine OutputDiscribErrMatlab




  subroutine WriteDiscAlgErrOutput()
    real, dimension(:,:), allocatable :: sum_err
    integer :: i, k
    integer :: text_size, is, num_size
    character(len=50) :: graf_name
    character(len=3) :: ch3

    allocate(sum_err(1:2, 1:ndim) )
    sum_err(:,:) = 0.

    do i= 1, grid%nelem
       do k=1, ndim
          sum_err(1, k) = sum_err(1, k) + grid%elem(i)%errDisc**2
          sum_err(2, k) = sum_err(2, k) + grid%elem(i)%errAlg**2
       enddo

    enddo

    do k=1, ndim
       sum_err(1, k) = sum_err(1, k)**0.5
       sum_err(2, k) = sum_err(2, k)**0.5
    enddo

    if(state%space%adapt%adapt_level > 0) then
      is = int(log(1.*state%space%adapt%adapt_level)/log(10.))
    else
      is = 0
    endif

    num_size = 3
    write( ch3, '(i3)' ) state%space%adapt%max_adapt_level

    graf_name = 'graf_DiscErr-A000    '
    text_size = 14
    graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
    if (stop_crit == 'L') then
       graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
    elseif (stop_crit == 'G') then
       graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
    elseif (stop_crit == 'N') then
       graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
    endif
    open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
    write(20+state%time%iter,'(i7, es14.6)') state%linSolver%iter_tot_SC, sum_err(1, 1)       ! the estimate according to which step is compared with the actual error
    if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '

    graf_name = 'graf_DiscErr-A000    '
    text_size = 14
    graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
    if (stop_crit == 'L') then
       graf_name(num_size+text_size+1:num_size+text_size+9) = 'L-adapit'
    elseif (stop_crit == 'G') then
       graf_name(num_size+text_size+1:num_size+text_size+9) = 'G-adapit'
    elseif (stop_crit == 'N') then
       graf_name(num_size+text_size+1:num_size+text_size+9) = 'N-adapit'
    endif
    open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
    write(20+state%time%iter,'(i7, es14.6)') state%space%adapt%adapt_level, sum_err(1, 1)       ! the estimate according to which step is compared with the actual error
    if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '


    graf_name = 'graf_AlgErr-A000    '
    text_size = 13
    graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
    if (stop_crit == 'L') then
       graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
    elseif (stop_crit == 'G') then
       graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
    elseif (stop_crit == 'N') then
       graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
    endif
    open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
    write(20+state%time%iter,'(i7, es14.6)') state%linSolver%iter_tot_SC, sum_err(2, 1)       ! the estimate according to which step is compared with the actual error
    if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '

    graf_name = 'graf_AlgErr-A000    '
    text_size = 13
    graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
    if (stop_crit == 'L') then
       graf_name(num_size+text_size+1:num_size+text_size+9) = 'L-adapit'
    elseif (stop_crit == 'G') then
       graf_name(num_size+text_size+1:num_size+text_size+9) = 'G-adapit'
    elseif (stop_crit == 'N') then
       graf_name(num_size+text_size+1:num_size+text_size+9) = 'N-adapit'
    endif
    open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
    write(20+state%time%iter,'(i7, es14.6)') state%space%adapt%adapt_level, sum_err(2, 1)       ! the estimate according to which step is compared with the actual error
    if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '

    close(20+state%time%iter)

    deallocate(sum_err)

  end subroutine WriteDiscAlgErrOutput



  subroutine TEST_ILU()
  type(Mblock), dimension(:, :), allocatable :: A, ALU, L,U
  type(Mblock) :: G, GLU, H, Loc
  real, dimension(:), allocatable :: p, b, f
  integer :: nb = 3
  integer :: nsize = 3
  integer :: i,j,k,k1, k2

  print*,'test ILU'

  allocate( A(1:nb, 1:nb), ALU(1:nb, 1:nb), L(1:nb, 1:nb), U(1:nb, 1:nb) )

  do i=1,nb
     do j=1,nb
        call InitMblock(A(i,j), nsize, nsize)
        call InitMblock(ALU(i,j), nsize, nsize)
        call InitMblock(L(i,j), nsize, nsize)
        call InitMblock(U(i,j), nsize, nsize)

        do k=1,nsize
           do k1=1,nsize
              A(i,j)%Mb(k,k1) = k*k1 +k + 1./(abs(i-j) +1)
              if(i==j)  A(i,j)%Mb(k,k1) =  4/(abs(k-k1) +1)

              if(i==j .and. k==k1)  L(i,j)%Mb(k,k1) = 1.
           enddo
        enddo
        !!!call WriteMblock(A(i,j))
     enddo
  enddo

  call InitMblock(G, nsize*nb, nsize*nb)
  call InitMblock(GLU, nsize*nb, nsize*nb)
  call InitMblock(H, nsize*nb, nsize*nb)



  do i=1,nb
     do j=1,nb

        G%Mb( (i-1)*nsize+1:i*nsize,(j-1)*nsize+1:j*nsize) = A(i,j)%Mb(1:nsize,1:nsize)

     enddo
  enddo

  print*,'Global matrix'
  call WriteMblock(G)


  call MblockLU( G, GLU )
  print*,'Global matrix inverse'
  call WriteMblock(GLU)


  print*,'test block ILU'

  ! local matrix for inversion
  call InitMblock(Loc, nsize, nsize)

  do k=1,nb
     do j=k,nb
        ALU(k,j)%Mb(:,:) = A(k,j)%Mb(:,:)
        do k1=1,k-1
           ALU(k,j)%Mb(:,:) = ALU(k,j)%Mb(:,:) &
                - matmul(ALU(k,k1)%Mb(:,:), ALU(k1,j)%Mb(:,:) )
        enddo
     enddo

     do i=k+1,nb
        if(k .ne. nb) then
           ALU(i,k)%Mb(:,:) = A(i,k)%Mb(:,:)

           do k1=1,k-1
              ALU(i,k)%Mb(:,:) = ALU(i,k)%Mb(:,:) &
                   - matmul(ALU(i,k1)%Mb(:,:), ALU(k1,k)%Mb(:,:) )

           enddo

           call MblockLU( ALU(k,k), Loc )

           ALU(i,k)%Mb(:,:) = matmul(ALU(i,k)%Mb(:,:), Loc%Mb(:,:) )

        endif
     enddo
  enddo

  do i=1,nb
     do j=1,nb
        if(i<=j) U(i,j)%Mb(:,:) = ALU(i,j)%Mb(:,:)
        if(i> j) L(i,j)%Mb(:,:) = ALU(i,j)%Mb(:,:)

        GLU%Mb( (i-1)*nsize+1:i*nsize,(j-1)*nsize+1:j*nsize) = ALU(i,j)%Mb(1:nsize,1:nsize)
     enddo
  enddo

  print*,'Global LU decomposition'
  call WriteMblock(GLU)


  do i=1,nb
     do j=1,nb
        Loc%Mb(:,:) = 0.
        do k=1,nb
           Loc%Mb(:,:) = Loc%Mb(:,:) + matmul(L(i,k)%Mb(:,:), U(k,j)%Mb(:,:) )
        enddo

        H%Mb((i-1)*nsize+1:i*nsize,(j-1)*nsize+1:j*nsize) = Loc%Mb(1:nsize, 1:nsize)
     enddo
  enddo
  print*,'Multiplication of LU'
  call WriteMblock(H)



  allocate( p(nsize*nb), b(nsize*nb), f(nsize*nb) )

  ! backward solution per BLOCKS: LU a = e_k,  e_k is a canonical basis
  do k=1, nb
     do i=1,nsize
        f(1:nb*nsize) = 0.
        f((k-1)*nsize + i) = 1.

        !print*,k,i,'***********************'
        !print*,f(:)

        b(1:nsize) = f(1:nsize) !!!/GLU%Mb(1,1)== 1 for matrix L
        do j=2, nb
           b((j-1)*nsize+1:j*nsize) = f((j-1)*nsize+1:j*nsize)
           do k1 = 1, j-1
              b((j-1)*nsize+1:j*nsize) = b((j-1)*nsize+1:j*nsize) &
                   - matmul(L(j,k1)%Mb(1:nsize,1:nsize), b((k1-1)*nsize+1:k1*nsize) )
           enddo
        enddo

        !print*,b(:)

        call MblockLU( U(nb,nb), Loc )
        p((nb-1)*nsize +1: nb*nsize) = matmul(Loc%Mb(:,:), b((nb-1)*nsize +1: nb*nsize) )

        do k1 = 1, nb -1
           j = nb -k1
           p((j-1)*nsize +1: j*nsize) = b((j-1)*nsize +1: j*nsize)

           !print*,'****',p(:)
           do k2 = j+1, nb
              p((j-1)*nsize +1: j*nsize) =  p((j-1)*nsize +1: j*nsize) &
                   - matmul(U(j,k2)%Mb(1:nsize,1:nsize), p((k2-1)*nsize +1: k2*nsize) )
           !print*,'****',p(:), U(j,k2)%Mb(1:nsize,1:nsize), p((k2-1)*nsize +1: k2*nsize), &
           !     matmul(U(j,k2)%Mb(1:nsize,1:nsize), p((k2-1)*nsize +1: k2*nsize) )

           enddo

           call MblockLU( U(j,j), Loc )
           p((j-1)*nsize +1: j*nsize) = matmul(Loc%Mb(:,:), p((j-1)*nsize +1: j*nsize) )

           !print*,'****',p(:)

        enddo
        !print*,p(:)

        H%Mb(1:nb*nsize, (k-1)*nsize+i) = p(1:nb*nsize)
     enddo
  enddo

  print*,'Inverse of G using LU decomposition'
  call WriteMblock(H)

  do i=1,nb*nsize
     write(*,'(20es10.3)') matmul(H%Mb(i,1:nb*nsize), G%Mb(1:nb*nsize,1:nb*nsize) )
  enddo

  deallocate( p, b, f )
  deallocate( A, ALU, L, U )

end subroutine TEST_ILU

subroutine TEST_ILU1()
  type(Mblock), dimension(:, :), allocatable :: A, ALU, L,U
  type(Mblock) :: G, GLU, H, Loc
  integer, dimension(:,:), allocatable :: isparse, iLsparse
  real, dimension(:), allocatable :: p, b, f
  integer :: nb
  integer :: nsize
  integer :: i,j,k,k1, in, k2, dof

  print*,'test ILU1'

  dof = grid%elem(1)%dof
  nb = grid%nelem
  nsize = ndim * dof

  allocate( A(1:nb, 1:nb), ALU(1:nb, 1:nb), L(1:nb, 1:nb), U(1:nb, 1:nb) )
  allocate(isparse(1:nb,1:nb), iLsparse(1:nb*nsize,1:nb*nsize) )

  do i=1,nb
     do j=1,nb
        call InitMblock(A(i,j), nsize, nsize)
        call InitMblock(ALU(i,j), nsize, nsize)
        call InitMblock(L(i,j), nsize, nsize)
        call InitMblock(U(i,j), nsize, nsize)

        do k=1,nsize
           do k1=1,nsize
              if(i==j .and. k==k1)  L(i,j)%Mb(k,k1) = 1.
           enddo
        enddo
        !!!call WriteMblock(A(i,j))
     enddo
  enddo

  isparse(:,:) = 0
  iLsparse(:,:) = 0

  do i=1,grid%nelem
     ! diagonal blocks
     do k = 1,ndim
        A(i,i)%Mb( (k-1)*dof +1:k*dof, (k-1)*dof +1:k*dof) = grid%elem(i)%Mass%Mb(1:dof,1:dof)*16
        grid%elem(i)%block(0)%Mb( (k-1)*dof +1:k*dof, (k-1)*dof +1:k*dof) = grid%elem(i)%Mass%Mb(1:dof,1:dof)*16
        isparse(i,i) = 1
     enddo
     ! off-diagonal blocks
     do k=1,grid%elem(i)%flen
        j = grid%elem(i)%face(neigh,k)
        if(j > 0) then
           do k1=1,nsize
              A(i,j)%Mb( k1,k1) =  grid%elem(i)%area*grid%elem(j)%area * 16
              grid%elem(i)%block(k)%Mb( k1,k1) =  grid%elem(i)%area*grid%elem(j)%area * 16
              isparse(i,j) = 1
           enddo
        endif
     enddo
  enddo


  call InitMblock(G, nsize*nb, nsize*nb)
  call InitMblock(GLU, nsize*nb, nsize*nb)
  call InitMblock(H, nsize*nb, nsize*nb)


  do i=1,nb
     do j=1,nb

        G%Mb( (i-1)*nsize+1:i*nsize,(j-1)*nsize+1:j*nsize) = A(i,j)%Mb(1:nsize,1:nsize)
        if(isparse(i,j) == 1) iLsparse( (i-1)*nsize+1:i*nsize,(j-1)*nsize+1:j*nsize) = 1
     enddo
  enddo

  print*,'Global matrix'
  call WriteMblock(G)

  do i=1,nb*nsize
     write(*,'(50i3)')  iLsparse(i,:)
  enddo

!  call MblockLU( G, GLU )
!  print*,'Global matrix inverse'
!  call WriteMblock(GLU)


  print*,'test block ILU'

  ! local matrix for inversion
  call InitMblock(Loc, nsize, nsize)

  do k=1,nb
     do j=k,nb
        ALU(k,j)%Mb(:,:) = A(k,j)%Mb(:,:)
        do k1=1,k-1
           ALU(k,j)%Mb(:,:) = ALU(k,j)%Mb(:,:) &
                - matmul(ALU(k,k1)%Mb(:,:), ALU(k1,j)%Mb(:,:) )
        enddo
     enddo

     do i=k+1,nb
        if(k .ne. nb) then
           ALU(i,k)%Mb(:,:) = A(i,k)%Mb(:,:)

           do k1=1,k-1
              ALU(i,k)%Mb(:,:) = ALU(i,k)%Mb(:,:) &
                   - matmul(ALU(i,k1)%Mb(:,:), ALU(k1,k)%Mb(:,:) )

           enddo

           call MblockLU( ALU(k,k), Loc )

           ALU(i,k)%Mb(:,:) = matmul(ALU(i,k)%Mb(:,:), Loc%Mb(:,:) )

        endif
     enddo
  enddo

  do i=1,nb
     do j=1,nb
        if(i<=j) U(i,j)%Mb(:,:) = ALU(i,j)%Mb(:,:)
        if(i> j) L(i,j)%Mb(:,:) = ALU(i,j)%Mb(:,:)

        GLU%Mb( (i-1)*nsize+1:i*nsize,(j-1)*nsize+1:j*nsize) = ALU(i,j)%Mb(1:nsize,1:nsize)
     enddo
  enddo

  print*,'Global LU decomposition'
  call WriteMblock(GLU)


  do i=1,nb
     do j=1,nb
        Loc%Mb(:,:) = 0.
        do k=1,nb
           Loc%Mb(:,:) = Loc%Mb(:,:) + matmul(L(i,k)%Mb(:,:), U(k,j)%Mb(:,:) )
        enddo

        H%Mb((i-1)*nsize+1:i*nsize,(j-1)*nsize+1:j*nsize) = Loc%Mb(1:nsize, 1:nsize)
     enddo
  enddo
  print*,'Multiplication of LU'
  call WriteMblock(H)

  print*
  print*,'._._._._._._._'


  do i=1,nb
     do j=1,nb
        ALU(i,j)%Mb(:,:) = 0.
     enddo
  enddo

 ! BLOCK ILU decomposition
  do i=1,grid%nelem
     ! diagonal term
     grid%elem(i)%ILU(0)%Mb(:,:) = grid%elem(i)%block(0)%Mb(:,:)

     do k=1,grid%elem(i)%flen
        in = grid%elem(i)%face(neigh,k)

        if(in > 0 .and. in < i) then
           grid%elem(i)%ILU(0)%Mb(:,:) = grid%elem(i)%ILU(0)%Mb(:,:) &
                - matmul(grid%elem(i)%ILU(k)%Mb(:,:), &
                grid%elem(in)%ILU(grid%elem(i)%face(nei_i,k) )%Mb(:,:) )
        endif
     enddo

     ! off diagonal terms in the row
     do k=1,grid%elem(i)%flen
        in = grid%elem(i)%face(neigh,k)
        if(in > i) then
           grid%elem(i)%ILU(k)%Mb(:,:) = grid%elem(i)%block(k)%Mb(:,:)
        endif
     enddo


     ! off diagonal terms in the columns
     call MblockLU( grid%elem(i)%ILU(0), Loc )

     do k=1,grid%elem(i)%flen
        in = grid%elem(i)%face(neigh,k)
        j = grid%elem(i)%face(nei_i,k)
        if(in > i) then
           grid%elem(in)%ILU(j)%Mb(:,:) &
                = matmul(grid%elem(in)%block(j)%Mb(:,:), Loc%Mb(:,:) )
        endif
     enddo
  enddo

  GLU%Mb(:,:) = 0.
  do i=1,grid%nelem
     GLU%Mb( (i-1)*nsize+1:i*nsize,(i-1)*nsize+1:i*nsize) = grid%elem(i)%ILU(0)%Mb(1:nsize,1:nsize)

     do j=1,grid%elem(i)%flen
        in = grid%elem(i)%face(neigh,j)
        if(in > 0) then
           GLU%Mb( (i-1)*nsize+1:i*nsize,(in-1)*nsize+1:in*nsize) = grid%elem(i)%ILU(j)%Mb(1:nsize,1:nsize)

        endif
     enddo
  enddo

  print*,'Global Block ILU decomposition'
  call WriteMblock(GLU)

  do i=1,nb
     do j=1,nb
        if(i<=j) U(i,j)%Mb(:,:) = GLU%Mb( (i-1)*nsize+1:i*nsize,(j-1)*nsize+1:j*nsize)
        if(i> j) L(i,j)%Mb(:,:) = GLU%Mb( (i-1)*nsize+1:i*nsize,(j-1)*nsize+1:j*nsize)

     enddo
  enddo



  do i=1,nb
     do j=1,nb
        Loc%Mb(:,:) = 0.
        do k=1,nb
           Loc%Mb(:,:) = Loc%Mb(:,:) + matmul(L(i,k)%Mb(:,:), U(k,j)%Mb(:,:) )
        enddo

        H%Mb((i-1)*nsize+1:i*nsize,(j-1)*nsize+1:j*nsize) = Loc%Mb(1:nsize, 1:nsize)
     enddo
  enddo
  print*,'Multiplication of LU'
  call WriteMblock(H)

  print*,'**************** difference LU and iLU multiplications'
  do i=1,nb*nsize
     !write(*,'(20es9.2)') G%Mb(i,:) - H%Mb(i,:)
     do j=1,nb*nsize
        if(iLsparse(i,j) <= 1 .and. abs(G%Mb(i,j) - H%Mb(i,j)) > 0.) &
             write(*,'(3i5,es14.7)') iLsparse(i,j), i,j,G%Mb(i,j) - H%Mb(i,j)
     enddo

  enddo
  print*,'**************** difference LU and iLU multiplications'

  stop


  allocate( p(nsize*nb), b(nsize*nb), f(nsize*nb) )

  ! backward solution per BLOCKS: LU a = e_k,  e_k is a canonical basis
  do k=1, nb
     do i=1,nsize
        f(1:nb*nsize) = 0.
        f((k-1)*nsize + i) = 1.

        !print*,k,i,'***********************'
        !print*,f(:)

        b(1:nsize) = f(1:nsize) !!!/GLU%Mb(1,1)== 1 for matrix L
        do j=2, nb
           b((j-1)*nsize+1:j*nsize) = f((j-1)*nsize+1:j*nsize)
           do k1 = 1, j-1
              b((j-1)*nsize+1:j*nsize) = b((j-1)*nsize+1:j*nsize) &
                   - matmul(L(j,k1)%Mb(1:nsize,1:nsize), b((k1-1)*nsize+1:k1*nsize) )
           enddo
        enddo

        !print*,b(:)

        call MblockLU( U(nb,nb), Loc )
        p((nb-1)*nsize +1: nb*nsize) = matmul(Loc%Mb(:,:), b((nb-1)*nsize +1: nb*nsize) )

        do k1 = 1, nb -1
           j = nb -k1
           p((j-1)*nsize +1: j*nsize) = b((j-1)*nsize +1: j*nsize)

           !print*,'****',p(:)
           do k2 = j+1, nb
              p((j-1)*nsize +1: j*nsize) =  p((j-1)*nsize +1: j*nsize) &
                   - matmul(U(j,k2)%Mb(1:nsize,1:nsize), p((k2-1)*nsize +1: k2*nsize) )
           !print*,'****',p(:), U(j,k2)%Mb(1:nsize,1:nsize), p((k2-1)*nsize +1: k2*nsize), &
           !     matmul(U(j,k2)%Mb(1:nsize,1:nsize), p((k2-1)*nsize +1: k2*nsize) )

           enddo

           call MblockLU( U(j,j), Loc )
           p((j-1)*nsize +1: j*nsize) = matmul(Loc%Mb(:,:), p((j-1)*nsize +1: j*nsize) )

           !print*,'****',p(:)

        enddo
        !print*,p(:)

        H%Mb(1:nb*nsize, (k-1)*nsize+i) = p(1:nb*nsize)
     enddo
  enddo

  print*,'Inverse of G using LU decomposition'
  call WriteMblock(H)

  do i=1,nb*nsize
     write(*,'(20es10.3)') matmul(H%Mb(i,1:nb*nsize), G%Mb(1:nb*nsize,1:nb*nsize) )
  enddo

  deallocate (A, ALU, L,U)
  deallocate ( p, b, f)

end subroutine TEST_ILU1




end module io_sub
