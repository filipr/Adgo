!> aposteriori error estimation subroutines for Vohralik approach for linear problem
!> estimates in the energy norm using RTN flux reconstruction including algebraic error
module alg_estim
  use rav_tho_nedelec
  use eval_sol
  use apost_estimation
  use euler_problem
  use loc_rav_tho_ned
  use eval_rav_tho_ned
  use dual_estim
  use error_subs

  implicit none

  public:: ComputeAlgEstim         !
  public:: ComputeElemAlgEstim     !
  !public:: ConstructFluxCD
  public:: ResidualElemAlgEstimCDSS      ! SS = Steady State
  public:: DataOscillationAlgEstimSS     ! SS = Steady State
  public:: FluxElemDiscEstimCDSS         ! SS = Steady State
  public:: FluxElemAlgEstimCDSS         ! SS = Steady State
  public:: FluxElemTotEstimCDSS         ! SS = Steady State
  !public:: ComputeLocRTNMomentumsElem
  !public:: FluxMM2RTN
  public:: RemainderElemAlgEstimSS       ! SS = Steady State
  public:: PotentNCElemEstimSS       ! SS = Steady State
  public:: FluxNCElemDiscEstimSS       ! SS = Steady State
  public:: FluxNCEstimSS       ! SS = Steady State
  public:: FluxNCElemEstimSSnew    ! SS = Steady State
  public:: ElemFaceFNC
  public:: ComputeFaceFNCestim
  public:: BoundTraceIneq
  public:: ComputeElemRTNMomentumCD
  public:: InitElemRTNMomentumCD

  !public:: CheckRTNMomentums   !! should be created based on ConstructFluxCD

contains
  !> estimates the computational error using the approach based on the flux reconstruction
  subroutine ComputeAlgEstim( )
    class(element), pointer :: elem
    type(basis_rtn_fe), pointer:: loc_RTN
    !real, dimension(:,:), allocatable :: MMinvRE ! inverse of Momentum matrix on elem
    integer :: dof, Fdeg, Fdof, iter_est_loc
    integer :: i
    real, dimension(:,:), allocatable :: estim, sum_estim

    real, dimension(:,:), pointer :: phi, phiT   !for TEST
    real, dimension(:,:), allocatable :: wi, q !for TEST
    real :: val
    integer :: k, Qdof, ist, l !for TEST
    character(len=50) :: gridfile
    integer :: text_size, is, num_size
    character(len=50) :: graf_name
    character(len=3) :: ch3
    character(len=12) :: ch12
    real :: DiscErr, AlgErr
    !real :: pok_Alg

    !pok_Alg = 0

    !print*,'##############  alg_estim.f90 state%time%iter ',state%time%iter
    !print*,'##############  alg_estim.f90 state%time%iter_loc',state%time%iter_loc

    !if (state%time%iter /= 1 .and. state%time%iter_loc == 1) state%linSolver%iter_tot = state%linSolver%iter_tot_SC

    !print*,'state%space%adapt%adapt_level, state%time%iter, state%time%iter_SC, state%linSolver%iter_tot, state%linSolver%iter_tot_SC:',  &
    !                                  state%space%adapt%adapt_level, state%time%iter, state%time%iter_SC, state%linSolver%iter_tot, state%linSolver%iter_tot_SC

    !write(*,'(a8, i2)') 'state%time%iter =', state%time%iter


    !Fdeg = state%space%deg-1
    Fdeg = state%space%deg


    !Fdof = state%RTN(Fdeg)%dof
    Fdof = SetRTNdof(Fdeg)

    loc_RTN => state%loc_RTN(Fdeg)
    if(.not. loc_RTN%defined ) call Init_Loc_RTN(state%loc_RTN(Fdeg), Fdeg )

    !do i = 1, grid%nelem
    !   elem => grid%elem(i)
    !   print*, 'Vypis matice predpodmineni:'
    !   do l=1, 6
    !     write(*, '(6es16.8)') elem%ILU(0)%Mb(l,:)
    !   enddo

       !write(*,*) 'elem%i, elem%HGnode, elem%xc(1:2) =', elem%i, elem%HGnode, elem%xc(1:2)
    !enddo
    !pause

    !if (state%time%iter == 4 .or. state%time%iter == 5 .or. state%time%iter == 6) then
    !       write(*,*) 'state%stop_rem', state%time%iter, state%stop_rem
    !       write(*,*) 'state%stop_alg', state%time%iter, state%stop_alg
    !       pause
    !endif

    do i = 1, grid%nelem
       elem => grid%elem(i)
       call ComputeElemRTNMomentumCD(elem, Fdeg, Fdof) ! for one elem RTN flux momentums of other elements are needed in FluxNCElemDiscEstimSS
       !if (elem%i == 1) write(*, '(a38, 15es14.6)') 'elem%RTNflux(0, 1:ndim, 1:Fdof)', elem%RTNflux(0, 1:ndim, 1:Fdof)
       !pause
    enddo  ! end of i=1,grid%nelem


    state%stop_rem = .true. !initialization after storing of the older values of elem%RTNflux for ALL ELEMENTS
    state%stop_rem_glob = .true. !initialization after storing of the older values of elem%RTNflux for ALL ELEMENTS

    allocate(estim(1:max_eta, 1:ndim), sum_estim(1:max_eta, 1:ndim) )
    sum_estim(:,:) = 0.
    estim(:,:) = 0.

    !if(state%time%iter == 1) state%err(IC_L2) = 0.

    !if(state%time%iter > 0) then           OLD VERSION 2
    !  ! computation of flux nonconformity estimator (discretization and algebraic components) for all elements at once
    !  call FluxNCEstimSS(Fdeg, Fdof, .true.)  ! discretization FNC
    !  call FluxNCEstimSS(Fdeg, Fdof, .false.)  ! algebraic FNC
    !endif

    do i = 1, grid%nelem
       elem => grid%elem(i)
       !allocate (elem%res_func(1:ndim, 1:elem%Qdof) )    DONE in sub. PrepareOneElement in problem.f90
       call Eval_res_func_Elem(elem) ! residual function r_T^i at volume integ nodes
       call ComputeElemAlgEstim(elem, Fdeg, Fdof, estim(1:max_eta, 1:ndim) )


       if(state%time%iter > 0) then
          if (estim(Rem, 1) > gamma_rem * (estim(Disc, 1) + estim(Alg, 1)) ) then
             !print*, '# Local quasi-equilibration condition is not satisfied.'
             state%stop_rem = .false.

             !return                          ! \nu := \nu + \nu (quasi-equilibration condition is not satisfied at least for one mesh element)
          endif

          if (estim(Alg, 1) > gamma_alg * estim(Disc, 1) ) then
             !print*, '# Local algebraic stopping criterion is not satisfied.'
             state%stop_alg = .false.
          endif

          !if (state%time%iter == 1 .and. i < 11 ) then
          !write(*,'(i3, a8, es16.8)')  i, 'etares =', estim(Resid, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etaFD =', estim(FD, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etarem =', estim(Rem, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etaPNC =', estim(PNC, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etaFA =', estim(FA, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etaFNCD =', estim(FNCD, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etaFNCA =', estim(FNCA, 1:ndim)
          !write(*,'(i3, a24, es16.8)')  state%time%iter, 'estim(Disc, 1:ndim) =', estim(Disc, 1:ndim)
          !write(*,'(i3, a24, es16.8)')  state%time%iter, 'estim(Alg, 1:ndim) =', estim(Alg, 1:ndim)
          !write(*,'(i3, a24, es16.8)')  state%time%iter, 'estim(Rem, 1:ndim) =', estim(Rem, 1:ndim)
          !pause
          !endif

          !if (state%time%iter == 57) then
          !    gridfile = 'square01.strX.HG.grid'
          !    call WriteMesh(gridfile, grid)
          !endif

          !if (state%time%iter < 0 ) then
          !write(*,'(i3, a8, es16.8)')  i, 'etares =', estim(Resid, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etaFD =', estim(FD, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etarem =', estim(Rem, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etaPNC =', estim(PNC, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etaFA =', estim(FA, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etaFNCD =', estim(FNCD, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etaFNCA =', estim(FNCA, 1:ndim)
          !write(*,'(i3, a24, es16.8)')  state%time%iter, 'estim(Disc, 1:ndim) =', estim(Disc, 1:ndim)
          !write(*,'(i3, a24, es16.8)')  state%time%iter, 'estim(Alg, 1:ndim) =', estim(Alg, 1:ndim)
          !write(*,'(i3, a24, es16.8)')  state%time%iter, 'estim(Rem, 1:ndim) =', estim(Rem, 1:ndim)
          !pause
          !endif

          !if (state%time%iter_loc == 4 .and. elem%i < 11 ) then
          !write(*,*) 'elem%HGnode =', elem%HGnode
          !write(*,'(i3, a8, es16.8)')  i, 'etares =', estim(Resid, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etaFD =', estim(FD, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etarem =', estim(Rem, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etaPNC =', estim(PNC, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etaFA =', estim(FA, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etaFNCD =', estim(FNCD, 1:ndim)
          !write(*,'(i3, a8, es16.8)')  i, 'etaFNCA =', estim(FNCA, 1:ndim)
          !write(*,'(i3, a24, es16.8)')  state%time%iter, 'estim(Disc, 1:ndim) =', estim(Disc, 1:ndim)
          !write(*,'(i3, a24, es16.8)')  state%time%iter, 'estim(Alg, 1:ndim) =', estim(Alg, 1:ndim)
          !write(*,'(i3, a24, es16.8)')  state%time%iter, 'estim(Rem, 1:ndim) =', estim(Rem, 1:ndim)
          !pause
          !endif

          !if (i == 120) then
          !write(*,'(i3, a24, es16.8)')  state%time%iter, 'estim(Disc, 1:ndim) =', estim(Disc, 1:ndim)
          !write(*,'(i3, a24, es16.8)')  state%time%iter, 'estim(Alg, 1:ndim) =', estim(Alg, 1:ndim)
          !write(*,'(i3, a24, es16.8)')  state%time%iter, 'estim(Rem, 1:ndim) =', estim(Rem, 1:ndim)
          !endif
          !if (i == 5) stop
       endif

       if(state%time%iter > 0) then
          elem%estimST = abs(estim(Disc, 1))
          sum_estim(Disc, 1:ndim) = sum_estim(Disc, 1:ndim) + (estim(Disc, 1:ndim)**2 )
          sum_estim(Alg, 1:ndim) = sum_estim(Alg, 1:ndim) + (estim(Alg, 1:ndim)**2 )
          sum_estim(Rem, 1:ndim) = sum_estim(Rem, 1:ndim) + (estim(Rem, 1:ndim)**2 )

          sum_estim(Resid, 1:ndim) = sum_estim(Resid, 1:ndim) + (estim(Resid, 1:ndim)**2 )
          sum_estim(FD, 1:ndim) = sum_estim(FD, 1:ndim) + (estim(FD, 1:ndim)**2 )
          sum_estim(PNC, 1:ndim) = sum_estim(PNC, 1:ndim) + (estim(PNC, 1:ndim)**2 )
          sum_estim(FA, 1:ndim) = sum_estim(FA, 1:ndim) + (estim(FA, 1:ndim)**2 )
          sum_estim(FNCD, 1:ndim) = sum_estim(FNCD, 1:ndim) + (estim(FNCD, 1:ndim)**2 )
          sum_estim(FNCA, 1:ndim) = sum_estim(FNCA, 1:ndim) + (estim(FNCA, 1:ndim)**2 )
          sum_estim(FT, 1:ndim) = sum_estim(FT, 1:ndim) + (  ( estim(FT, 1:ndim)+estim(Resid, 1:ndim) + estim(FNCT, 1:ndim) )**2 )
                                                           ! sum of total flux, residual and total flux nonconformity estimators!!!!
          !if (estim(Resid, 1) > 0.1) then
          !write(*,*) 'elem%i = ', elem%i
          !write(*,'(a16, 30es12.4)') 'estim(Resid, 1:ndim)**2:', estim(Resid, 1:ndim)**2
          !write(*,'(a16, 30es12.4)') 'sum_estim(Resid, 1:ndim):', sum_estim(Resid, 1:ndim)
          !pause
          !endif

          !pok_Alg =  pok_Alg + (estim(FA, 1) + estim(FNCA, 1) )**2

          !time dep: sum_estim(IC, 1:ndim) = sum_estim(IC, 1:ndim) + estim(IC, 1:ndim)
       endif

    enddo  ! end of i=1,grid%nelem
    !pause

    state%estim(Disc:Rem, 1:ndim) = sum_estim(Disc:Rem, 1:ndim)

    state%estim(Tot, 1:ndim) = ( (sum_estim(FT, 1:ndim)**0.5 + state%estim(Rem, 1:ndim)**0.5 )**2 + sum_estim(PNC, 1:ndim) )**0.5
                           ! guaranteed upper bound   !in sum_estim(FT, 1:ndim) also Residual estimator and total flux nonconformity estimator included!!!

    !state%estim(Tot, 1:ndim) = (2*state%estim(Disc, 1:ndim) )**0.5 &         ! guaranteed upper bound
    !     + state%estim(Alg, 1:ndim)**0.5 + state%estim(Rem, 1:ndim)**0.5

    if(state%time%iter > 0) then
          if (state%estim(Rem, 1)**0.5 > gamma_rem * (state%estim(Disc, 1)**0.5 + state%estim(Alg, 1)**0.5) ) then
             !print*, '# Global quasi-equilibration condition is not satisfied.'
             state%stop_rem_glob = .false.
             !return                          ! \nu := \nu + \nu (quasi-equilibration condition is not satisfied at least for one mesh element)
          endif

          if (state%estim(Alg, 1)**0.5 > gamma_alg * state%estim(Disc, 1)**0.5 ) then
             !print*, '# Global algebraic stopping criterion is not satisfied.'
             !write(*,*) 'ALG ERROR:', state%estim(Alg, 1)**0.5
             !write(*,*) 'DISC ERROR:', state%estim(Disc, 1)**0.5

             state%stop_alg_glob = .false.
          endif
    endif


    !!!TEST
    if (state%time%iter < -8) then

        if(state%space%adapt%adapt_level > 0) then
             is = int(log(1.*state%space%adapt%adapt_level)/log(10.))
        else
             is = 0
        endif

        num_size = 3
        write( ch3, '(i3)' ) state%space%adapt%adapt_level

        graf_name = 'graf_Refin-000    '
             text_size = 11
             graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
             if (stop_crit == 'L') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
             elseif (stop_crit == 'G') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
             elseif (stop_crit == 'N') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
             endif
             open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')

             do i=1,grid%nelem
               elem => grid%elem(i)
               write(20+state%time%iter,'(2i7, es14.6)') state%linSolver%iter_tot - nu, elem%i, elem%estimST
             enddo
             write(20+state%time%iter,'(a11)') ' '
             write(20+state%time%iter,'(a11)') '***********************************************************'
             write(20+state%time%iter,'(a11)') ' '
             close(20+state%time%iter)

    endif
    !!!TEST


    if (state%time%iter > 0) then
       !if (state%time%iter == 1 .or. state%time%iter_loc /= 1) then
          if (state%time%iter_SC == -1) then  ! computational solution has not been obtained yet

         if(state%space%adapt%adapt_level > 0) then
             is = int(log(1.*state%space%adapt%adapt_level)/log(10.))
             !if ( stop_crit /= 'N' ) then
             !   iter_est_loc = (state%time%iter_loc - 2)*nu  ! in which step the estimate is computed
             !else
             !   iter_est_loc = state%linSolver%iter_tot ! no stop criteria applied
             !endif
          else
             is = 0
             !if ( stop_crit /= 'N' ) then
             !   iter_est_loc = (state%time%iter_loc - 1)*nu  ! in which step the estimate is computed
             !else
             !   iter_est_loc = state%linSolver%iter_tot ! no stop criteria applied
             !endif
         endif

         if ( stop_crit /= 'N' ) then
            !iter_est_loc = (state%time%iter_loc - 1)*nu  ! in which step the estimate is computed
            iter_est_loc = state%linSolver%iter_tot - nu ! convergence of the estimators in one graph
         else
            !iter_est_loc = (state%time%iter_loc - 1)*nu  !(state%linSolver%iter - nu) should be the same ! classical stop criteria applied
            iter_est_loc = state%linSolver%iter_tot - nu ! convergence of the estimators in one graph
         endif

         !if (state%space%adapt%adapt_level == 2 .or. state%space%adapt%adapt_level == 3) then
         !   if (iter_est_loc > 200) then
         !      write(*,'(i3, l3)')  iter_est_loc, state%stop_rem
         !      pause
         !   endif
         !endif
         !write(*,'(i3)')  state%time%iter
         !write(*,'(a24, es16.8)')  'sum_estim(Alg, 1:ndim) =', sum_estim(Alg, 1:ndim)**0.5
         !write(*,'(a24, es16.8)')  'pok_Alg**0.5 =', pok_Alg**0.5
         !pause

         num_size = 3
         write( ch3, '(i3)' ) state%space%adapt%adapt_level
         !write(*, *) 'ch3', ch3

         call ComputeDiscAlgError(DiscErr, AlgErr)

         state%err(algeb) = AlgErr     ! algebraic error at the current iteration

         do l=FNCD, Tot
          !write(*,*) 'l:', l
          if (l == Resid) then
             graf_name = 'graf_Resid-000    '
             text_size = 11
             graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
             if (stop_crit == 'L') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
             elseif (stop_crit == 'G') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
             elseif (stop_crit == 'N') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
             endif
             open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')


             if (stop_crit == 'L') then
                if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Resid, 1)**0.5
                if(state%stop_rem .and. state%stop_alg) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'G') then
                if((state%stop_rem_glob) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Resid, 1)**0.5
                if(state%stop_rem_glob .and. state%stop_alg_glob) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'N') then
                !if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Resid, 1)**0.5
                if ((state%err(SSL8) < state%conv_rez) .and. (state%err(algeb) < state%conv_rez)) write(20+state%time%iter,*) ' '
             else
                    Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
             endif


          elseif (l == FD) then
             graf_name = 'graf_FD-000    '
             text_size = 8
             graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
             if (stop_crit == 'L') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
             elseif (stop_crit == 'G') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
             elseif (stop_crit == 'N') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
             endif
             open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')

             if (stop_crit == 'L') then
                if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FD, 1)**0.5
                if(state%stop_rem .and. state%stop_alg) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'G') then
                if((state%stop_rem_glob) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FD, 1)**0.5
                if(state%stop_rem_glob .and. state%stop_alg_glob) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'N') then
                !if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FD, 1)**0.5
                if ((state%err(SSL8) < state%conv_rez) .and. (state%err(algeb) < state%conv_rez)) write(20+state%time%iter,*) ' '
             else
                    Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
             endif

          elseif (l == Rem) then
             graf_name = 'graf_Rem-000    '
             text_size = 9
             graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
             if (stop_crit == 'L') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
             elseif (stop_crit == 'G') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
             elseif (stop_crit == 'N') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
             endif
             open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')

             if (stop_crit == 'L') then
                if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Rem, 1)**0.5
                if(state%stop_rem .and. state%stop_alg) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'G') then
                if((state%stop_rem_glob) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Rem, 1)**0.5
                if(state%stop_rem_glob .and. state%stop_alg_glob) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'N') then
                !if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Rem, 1)**0.5
                if ((state%err(SSL8) < state%conv_rez) .and. (state%err(algeb) < state%conv_rez)) write(20+state%time%iter,*) ' '
             else
                    Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
             endif


          elseif (l == PNC) then
             graf_name = 'graf_PNC-000    '
             text_size = 9
             graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
             if (stop_crit == 'L') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
             elseif (stop_crit == 'G') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
             elseif (stop_crit == 'N') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
             endif
             open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')

             if (stop_crit == 'L') then
                if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(PNC, 1)**0.5
                if(state%stop_rem .and. state%stop_alg) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'G') then
                if((state%stop_rem_glob) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(PNC, 1)**0.5
                if(state%stop_rem_glob .and. state%stop_alg_glob) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'N') then
                !if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(PNC, 1)**0.5
                if ((state%err(SSL8) < state%conv_rez) .and. (state%err(algeb) < state%conv_rez)) write(20+state%time%iter,*) ' '
             else
                    Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
             endif

          elseif (l == FA) then
             graf_name = 'graf_FA-000    '
             text_size = 8
             graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
             if (stop_crit == 'L') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
             elseif (stop_crit == 'G') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
             elseif (stop_crit == 'N') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
             endif
             open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')

             if (stop_crit == 'L') then
                if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FA, 1)**0.5
                if(state%stop_rem .and. state%stop_alg) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'G') then
                if((state%stop_rem_glob) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FA, 1)**0.5
                if(state%stop_rem_glob .and. state%stop_alg_glob) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'N') then
                !if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FA, 1)**0.5
                if ((state%err(SSL8) < state%conv_rez) .and. (state%err(algeb) < state%conv_rez)) write(20+state%time%iter,*) ' '
             else
                    Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
             endif

          elseif (l == FNCD) then
             graf_name = 'graf_FNCD-000    '
             text_size = 10
             graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
             if (stop_crit == 'L') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
             elseif (stop_crit == 'G') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
             elseif (stop_crit == 'N') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
             endif
             open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')

             if (stop_crit == 'L') then
                if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FNCD, 1)**0.5
                if(state%stop_rem .and. state%stop_alg) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'G') then
                if((state%stop_rem_glob) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FNCD, 1)**0.5
                if(state%stop_rem_glob .and. state%stop_alg_glob) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'N') then
                !if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FNCD, 1)**0.5
                if ((state%err(SSL8) < state%conv_rez) .and. (state%err(algeb) < state%conv_rez)) write(20+state%time%iter,*) ' '
             else
                    Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
             endif

          elseif (l == FNCA) then
             graf_name = 'graf_FNCA-000    '
             text_size = 10
             graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
             if (stop_crit == 'L') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
             elseif (stop_crit == 'G') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
             elseif (stop_crit == 'N') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
             endif
             open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')

             if (stop_crit == 'L') then
                if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FNCA, 1)**0.5
                if(state%stop_rem .and. state%stop_alg) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'G') then
                if((state%stop_rem_glob) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FNCA, 1)**0.5
                if(state%stop_rem_glob .and. state%stop_alg_glob) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'N') then
                !if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FNCA, 1)**0.5
                if ((state%err(SSL8) < state%conv_rez) .and. (state%err(algeb) < state%conv_rez)) write(20+state%time%iter,*) ' '
             else
                    Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
             endif


          elseif (l == Disc) then
             graf_name = 'graf_Disc-000    '
             text_size = 10
             graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
             if (stop_crit == 'L') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
             elseif (stop_crit == 'G') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
             elseif (stop_crit == 'N') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
             endif
             open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')

             if (stop_crit == 'L') then
                if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Disc, 1)**0.5
                if(state%stop_rem .and. state%stop_alg) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'G') then
                if((state%stop_rem_glob) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Disc, 1)**0.5
                if(state%stop_rem_glob .and. state%stop_alg_glob) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'N') then
                !if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Disc, 1)**0.5
                if ((state%err(SSL8) < state%conv_rez) .and. (state%err(algeb) < state%conv_rez)) write(20+state%time%iter,*) ' '
             else
                    Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
             endif


          elseif (l == Alg) then
             graf_name = 'graf_Alg-000    '
             text_size = 9
             graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
             if (stop_crit == 'L') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
             elseif (stop_crit == 'G') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
             elseif (stop_crit == 'N') then
                graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
             endif
             open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')

             if (stop_crit == 'L') then
                if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Alg, 1)**0.5
                if(state%stop_rem .and. state%stop_alg) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'G') then
                if((state%stop_rem_glob) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Alg, 1)**0.5
                if(state%stop_rem_glob .and. state%stop_alg_glob) write(20+state%time%iter,*) ' '
             elseif (stop_crit == 'N') then
                !if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Alg, 1)**0.5
                if ((state%err(SSL8) < state%conv_rez) .and. (state%err(algeb) < state%conv_rez)) write(20+state%time%iter,*) ' '
             else
                    Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
             endif

           elseif (l == Tot) then
              graf_name = 'graf_Est-000    '
              text_size = 9
              graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
              if (stop_crit == 'L') then
                 graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
              elseif (stop_crit == 'G') then
                 graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
              elseif (stop_crit == 'N') then
                 graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
              endif
              open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')

              if (stop_crit == 'L') then
                     if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                        write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, state%estim(Tot, 1:ndim)
                     if(state%stop_rem .and. state%stop_alg) write(20+state%time%iter,*) ' '
              elseif (stop_crit == 'G') then
                     if((state%stop_rem_glob) .or. (state%time%iter_loc == 1))&
                        write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, state%estim(Tot, 1:ndim)
                     if(state%stop_rem_glob .and. state%stop_alg_glob) write(20+state%time%iter,*) ' '
              elseif (stop_crit == 'N') then
                     !if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                        write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, state%estim(Tot, 1:ndim)
                     if ((state%err(SSL8) < state%conv_rez) .and. &
                        (state%err(algeb) < state%conv_rez)) write(20+state%time%iter,*) ' '
              else
                     Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
              endif

          endif


         enddo

         graf_name = 'graf_ErrXnorm-000    '
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

         if (stop_crit == 'L') then
            if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                 write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, state%errSTnorm(L2H1)**0.5
            if(state%stop_rem .and. state%stop_alg) write(20+state%time%iter,*) ' '
         elseif (stop_crit == 'G') then
                if((state%stop_rem_glob) .or. (state%time%iter_loc == 1))&
                     write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, state%errSTnorm(L2H1)**0.5
                if(state%stop_rem_glob .and. state%stop_alg_glob) write(20+state%time%iter,*) ' '
         elseif (stop_crit == 'N') then
                !if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                      write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, state%errSTnorm(L2H1)**0.5
                if ((state%err(SSL8) < state%conv_rez) .and. (state%err(algeb) < state%conv_rez)) write(20+state%time%iter,*) ' '
         else
                Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
         endif


         if (state%space%adapt%adapt_method == 'ALG2') then

            graf_name = 'graf_DiscErr-000    '
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

            if (stop_crit == 'L') then
               if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                      write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, DiscErr
               if(state%stop_rem .and. state%stop_alg) write(20+state%time%iter,*) ' '
            elseif (stop_crit == 'G') then
                   if((state%stop_rem_glob) .or. (state%time%iter_loc == 1))&
                      write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, DiscErr
                   if(state%stop_rem_glob .and. state%stop_alg_glob) write(20+state%time%iter,*) ' '
            elseif (stop_crit == 'N') then
                   !if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                      write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, DiscErr
                   if ((state%err(SSL8) < state%conv_rez) .and. (state%err(algeb) < state%conv_rez)) write(20+state%time%iter,*) ' '
            else
                   Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
            endif



            graf_name = 'graf_AlgErr-000    '
            text_size = 12
            graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
            if (stop_crit == 'L') then
               graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
            elseif (stop_crit == 'G') then
               graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
            elseif (stop_crit == 'N') then
               graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
            endif
            open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')

            if (stop_crit == 'L') then
               if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                   write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, AlgErr
               if(state%stop_rem .and. state%stop_alg) write(20+state%time%iter,*) ' '
            elseif (stop_crit == 'G') then
                   if((state%stop_rem_glob) .or. (state%time%iter_loc == 1))&
                      write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, AlgErr
                   if(state%stop_rem_glob .and. state%stop_alg_glob) write(20+state%time%iter,*) ' '
            elseif (stop_crit == 'N') then
                   !if((state%stop_rem) .or. (state%time%iter_loc == 1))&
                      write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, AlgErr
                   if ((state%err(SSL8) < state%conv_rez) .and. (state%err(algeb) < state%conv_rez)) write(20+state%time%iter,*) ' '
            else
                   Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
            endif

         endif

         close(20+state%time%iter)
         !write(*, *) 'int(log(1.*state%isol)/log(10.))', int(log(1.*state%isol)/log(10.))
         !write(*, *) 'log(1.*state%isol)/log(10.)', log(1.*state%isol)/log(10.)
         !write(*, *) 'state%isol', state%isol
         !pause
        endif
      !endif
    endif

    !if (state%time%iter > 0) then
    !   if (.not. state%stop_rem) then
    !      print*, '# Local quasi-equilibration condition is not satisfied.'
    !   endif
    !   if (.not. state%stop_alg) then
    !      print*, '# Local algebraic stopping criterion is not satisfied.'
    !   endif
    !
    !   if (.not. state%stop_rem_glob) then
    !      print*, '# Global quasi-equilibration condition is not satisfied.'
    !   endif
    !   if (.not. state%stop_alg_glob) then
    !      print*, '# Global algebraic stopping criterion is not satisfied.'
    !   endif
    !endif

    ! NEW print OUTPUT:
    ch12(1:12) = '............'
    if (state%time%iter > 0) then
       if (.not. state%stop_rem) ch12(2:3) = 'Le'
       !print*, '# Local quasi-equilibration condition is not satisfied.'
       if (.not. state%stop_alg) ch12(5:6) = 'La'
       !print*, '# Local algebraic stopping criterion is not satisfied.'

       if (.not. state%stop_rem_glob) ch12(8:9) = 'Ge'
       !print*, '# Global quasi-equilibration condition is not satisfied.'
       if (.not. state%stop_alg_glob) ch12(11:12) = 'Ga'
       !print*, '# Global algebraic stopping criterion is not satisfied.'
    endif

    write(*,'(a40, 2i5, a16)') &
         '## alg_estim.f90: iter, iter_loc, criter ', &
         state%time%iter,state%time%iter_loc,ch12(1:12)  !, &
         !estim(Rem, 1) / (estim(Disc, 1) + estim(Alg, 1)), &
         !estim(Alg, 1) / estim(Disc, 1)


    !write(*,'(i3, a24, es16.8)')  state%time%iter, 'sum_estim(Disc, 1:ndim) =', sum_estim(Disc, 1:ndim)
    !write(*,'(i3, a24, es16.8)')  state%time%iter, 'sum_estim(Alg, 1:ndim) =', sum_estim(Alg, 1:ndim)
    !write(*,'(i3, a24, es16.8)')  state%time%iter, 'sum_estim(Rem, 1:ndim) =', sum_estim(Rem, 1:ndim)

    !if (state%time%iter == 1) then
    !  write(*,'(a27)') '!!!STOP TEST resid estim'
    !  stop
    !endif

    !if(state%time%iter == 1) state%err(IC_L2) = state%err(IC_L2)**0.5



    !if (state%time%iter == 10) then
    !   write(*,'(i3, a24, es16.8)')  state%time%iter, 'sum_estim(Disc, 1:ndim) =', sum_estim(Disc, 1:ndim)
    !   write(*,'(i3, a24, es16.8)')  state%time%iter, 'sum_estim(Alg, 1:ndim) =', sum_estim(Alg, 1:ndim)
    !   write(*,'(i3, a24, es16.8)')  state%time%iter, 'sum_estim(Rem, 1:ndim) =', sum_estim(Rem, 1:ndim)
    !   pause
    !endif

    !write(*,'(i3, a29, 3es16.8)')  state%time%iter, 'state%estim(Disc:Rem, 1:ndim) =', state%estim(Disc:Rem, 1:ndim)
    !write(*,'(a33, 3es16.8)')  'state%estim(Disc:Rem, 1:ndim)**0.5 =', state%estim(Disc:Rem, 1:ndim)**0.5

    !!!TEST ALG vs. ALG2

    if (state%space%adapt%adapt_level > 0 .and. state%time%iter_loc == 1) then
       write(*,'(a29, 3es16.8)') 'state%estim(Resid, 1:ndim) =', sum_estim(Resid, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FD, 1:ndim) =', sum_estim(FD, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(PNC, 1:ndim) =', sum_estim(PNC, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FA, 1:ndim) =', sum_estim(FA, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FNCD, 1:ndim) =', sum_estim(FNCD, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FNCA, 1:ndim) =', sum_estim(FNCA, 1:ndim)**0.5

       write(*,'(a29, 3es16.8)') 'state%estim(Disc, 1:ndim) =', state%estim(Disc, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(Alg, 1:ndim) =', state%estim(Alg, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(Rem, 1:ndim) =', state%estim(Rem, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(Tot, 1:ndim) =', state%estim(Tot, 1:ndim)
       !pause
    endif

    if (stop_crit == 'L') then

    if(state%stop_rem .and. state%stop_alg) then
       write(*,'(a29, 3es16.8)') 'state%estim(Resid, 1:ndim) =', sum_estim(Resid, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FD, 1:ndim) =', sum_estim(FD, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(PNC, 1:ndim) =', sum_estim(PNC, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FA, 1:ndim) =', sum_estim(FA, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FNCD, 1:ndim) =', sum_estim(FNCD, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FNCA, 1:ndim) =', sum_estim(FNCA, 1:ndim)**0.5

       write(*,'(a29, 3es16.8)') 'state%estim(Disc, 1:ndim) =', state%estim(Disc, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(Alg, 1:ndim) =', state%estim(Alg, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(Rem, 1:ndim) =', state%estim(Rem, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(Tot, 1:ndim) =', state%estim(Tot, 1:ndim)
    endif

    elseif (stop_crit == 'G') then

    if(state%stop_rem_glob .and. state%stop_alg_glob) then
       write(*,'(a29, 3es16.8)') 'state%estim(Resid, 1:ndim) =', sum_estim(Resid, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FD, 1:ndim) =', sum_estim(FD, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(PNC, 1:ndim) =', sum_estim(PNC, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FA, 1:ndim) =', sum_estim(FA, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FNCD, 1:ndim) =', sum_estim(FNCD, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FNCA, 1:ndim) =', sum_estim(FNCA, 1:ndim)**0.5

       write(*,'(a29, 3es16.8)') 'state%estim(Disc, 1:ndim) =', state%estim(Disc, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(Alg, 1:ndim) =', state%estim(Alg, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(Rem, 1:ndim) =', state%estim(Rem, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(Tot, 1:ndim) =', state%estim(Tot, 1:ndim)
    endif

    elseif (stop_crit == 'N') then

    if (state%stop_rem_glob .and. state%stop_alg_glob) then           ! state%linSolver%lin_solver_not_conv == 0
       write(*,'(a29, 3es16.8)') 'state%estim(Resid, 1:ndim) =', sum_estim(Resid, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FD, 1:ndim) =', sum_estim(FD, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(PNC, 1:ndim) =', sum_estim(PNC, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FA, 1:ndim) =', sum_estim(FA, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FNCD, 1:ndim) =', sum_estim(FNCD, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(FNCA, 1:ndim) =', sum_estim(FNCA, 1:ndim)**0.5

       write(*,'(a29, 3es16.8)') 'state%estim(Disc, 1:ndim) =', state%estim(Disc, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(Alg, 1:ndim) =', state%estim(Alg, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(Rem, 1:ndim) =', state%estim(Rem, 1:ndim)**0.5
       write(*,'(a29, 3es16.8)') 'state%estim(Tot, 1:ndim) =', state%estim(Tot, 1:ndim)
    endif

    else
      Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
    endif



    !write(100,'(a6,4es12.4)') 'sum?',state%estim(Disc:Rem, 1)**0.5, state%estim(Tot, 1:ndim)

    if (stop_crit == 'L') then

    if (state%time%iter > 0) then
      !if (state%time%iter == 1 .or. state%time%iter_loc /= 1) then
         if (state%stop_rem) then
            if (state%stop_alg) then
               if (state%time%iter_SC == -1) then   ! "state%time%iter_SC =" index of a time step when Stopping Criteria based on AEE are satisfied
                  print*, '# Computational solution based on AEE ALG satisfying local stop criteria obtained.'
                  ! no additional GMRES step!!! guaranteed by ln 862 in compute.f90

                  state%time%iter_SC = state%time%iter - 1  ! state%time%iter corresponds to "i + nu" NOT "i"!!
                                                  ! state%time%iter_SC initialized in sub. InitProblem in problem.f90
                                                  ! and reinit. for adaptation in sub. ReprepareProblem in adaptation.f90
                  state%linSolver%iter_tot_SC = state%linSolver%iter_tot - nu   ! total number of GMRES iterations corresponding to "i"
               endif
            endif
         endif
      !endif
    endif

    elseif (stop_crit == 'G') then

    if (state%time%iter > 0) then
      !if (state%time%iter == 1 .or. state%time%iter_loc /= 1) then
         if (state%stop_rem_glob) then
            if (state%stop_alg_glob) then
               if (state%time%iter_SC == -1) then   ! "state%time%iter_SC =" index of a time step when Stopping Criteria based on AEE are satisfied
                  print*, '# Computational solution based on AEE ALG satisfying global stop criteria obtained.'
                  ! no additional GMRES step!!! guaranteed by ln 862 in compute.f90

                  state%time%iter_SC = state%time%iter - 1  ! state%time%iter corresponds to "i + nu" NOT "i"!!
                                                  ! state%time%iter_SC initialized in sub. InitProblem in problem.f90
                                                  ! and reinit. for adaptation in sub. ReprepareProblem in adaptation.f90
                  state%linSolver%iter_tot_SC = state%linSolver%iter_tot - nu   ! total number of GMRES iterations corresponding to "i"
               endif
            endif
         endif
      !endif
    endif

    elseif (stop_crit == 'N') then  ! "state%linSolver%lin_solver_not_conv == 0" is never satisfied because "nu" iterations are not
                                    ! sufficient for GMRES convergence!!!
    if (state%time%iter > 0) then
      !if (state%time%iter == 1 .or. state%time%iter_loc /= 1) then
        if ((state%err(SSL8) < state%conv_rez) .and. (state%err(algeb) < state%conv_rez)) then !(state%err(SSL8) < state%conv_rez .and. state%linSolver%lin_solver_not_conv == 0) then
            print*, '# Computational solution obtained. Classical stop criteria have been applied.'
            state%time%iter_SC = state%time%iter - 1

            state%linSolver%iter_tot_SC = state%linSolver%iter_tot - nu  ! total number of GMRES iterations corresponding to "i"

        endif
      !endif
    endif


    else
      Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (classical stop. crit.).'
    endif



    if (state%time%iter == state%time%iter_SC + 1 .and. state%time%iter > 0) then

       do i = 1, grid%nelem
          elem => grid%elem(i)
          !allocate(elem%wc(0:1, 1:ndim*elem%Qdof) )    DONE in sub. PrepareOneElement in problem.f90

          if (state%space%adapt%adapt_method == 'ALG') then
             ! storing of the values of the computational solution u_h^i, computation will continue to obtain u_h
             if (stop_crit == 'N') then
                elem%wc(0, :) = elem%w(1, :)  ! CHANGED!!! elem%w(1, :) is considered as a comput sol
                elem%wc(1, :) = elem%w(2, :)  ! to be able to compute its estimators using information in "+ nu" iteration

             else
                elem%wc(0, :) = elem%w(1, :)  ! in elem%w(0, :) is stored u_h^{i+nu}!!!
                elem%wc(1, :) = elem%w(2, :)
             endif

          endif

          ! storing of the estimate of total error u - u_h^i, discretization error u - u_h, and algebraic error u_h - u_h^i
          elem%estTot = elem%eta(Tot, 1)
          elem%estDisc = elem%eta(Disc, 1)
          elem%estAlg = elem%eta(Alg, 1)
          elem%estRem = elem%eta(Rem, 1)
       enddo

       !!state%errSTnorm(L2X_sc) = state%errSTnorm(L2H1)   ! SQUARED!!!   state%errXnorm = state%err(H1_old)**2 in sub. PassToNewTimeStep
       state%estimTot_SC = state%estim(Tot, 1)

       if (state%space%adapt%max_adapt_level > 0) then  ! saving values of errors, estimates, effectivity indices, ... at the end of one adapt cycle

           if(state%space%adapt%adapt_level > 0) then
             is = int(log(1.*state%space%adapt%adapt_level)/log(10.))
           else
             is = 0
           endif

           if ( stop_crit /= 'N' ) then
              iter_est_loc = state%linSolver%iter_tot - nu  ! in which step of the adapt process the estimate is computed,
                                                    ! i.e. the number of GMRES iterations is cumulated
           else
              iter_est_loc = state%linSolver%iter_tot - nu ! classical stop criteria applied
           endif

           if (state%model%icase == 35) then
              num_size = 3
              write( ch3, '(i3)' ) state%space%adapt%adapt_level

              graf_name = 'Res-000    '
              text_size = 4
              graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
              if (stop_crit == 'L') then
                 graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
              elseif (stop_crit == 'G') then
                 graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
              elseif (stop_crit == 'N') then
                 graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
              endif
              graf_name(num_size+text_size+2:num_size+text_size+5) = '.txt'
              open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
              do i=1, grid%nelem
                 write(20+state%time%iter,'(es14.6)') grid%elem(i)%eta(Resid, 1)
              enddo
              close(20+state%time%iter)
           endif


           num_size = 3
           write( ch3, '(i3)' ) state%space%adapt%max_adapt_level

           do l=FNCD, Rem
            !write(*,*) 'l:', l
            if (l == Resid) then
               graf_name = 'graf_Resid-A000    '
               text_size = 12
               graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
               if (stop_crit == 'L') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
               elseif (stop_crit == 'G') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
               elseif (stop_crit == 'N') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
               endif
               open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
               write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Resid, 1)**0.5
               if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '


            elseif (l == FD) then
               graf_name = 'graf_FD-A000    '
               text_size = 9
               graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
               if (stop_crit == 'L') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
               elseif (stop_crit == 'G') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
               elseif (stop_crit == 'N') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
               endif
               open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
               write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FD, 1)**0.5
               if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '

            elseif (l == Rem) then
               graf_name = 'graf_Rem-A000    '
               text_size = 10
               graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
               if (stop_crit == 'L') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
               elseif (stop_crit == 'G') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
               elseif (stop_crit == 'N') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
               endif
               open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
               write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Rem, 1)**0.5
               if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '


            elseif (l == PNC) then
               graf_name = 'graf_PNC-A000    '
               text_size = 10
               graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
               if (stop_crit == 'L') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
               elseif (stop_crit == 'G') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
               elseif (stop_crit == 'N') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
               endif
               open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
               write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(PNC, 1)**0.5
               if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '

            elseif (l == FA) then
               graf_name = 'graf_FA-A000    '
               text_size = 9
               graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
               if (stop_crit == 'L') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
               elseif (stop_crit == 'G') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
               elseif (stop_crit == 'N') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
               endif
               open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
               write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FA, 1)**0.5
               if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '

            elseif (l == FNCD) then
               graf_name = 'graf_FNCD-A000    '
               text_size = 11
               graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
               if (stop_crit == 'L') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
               elseif (stop_crit == 'G') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
               elseif (stop_crit == 'N') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
               endif
               open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
               write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FNCD, 1)**0.5
               if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '

            elseif (l == FNCA) then
               graf_name = 'graf_FNCA-A000    '
               text_size = 11
               graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
               if (stop_crit == 'L') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
               elseif (stop_crit == 'G') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
               elseif (stop_crit == 'N') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
               endif
               open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
               write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(FNCA, 1)**0.5
               if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '


            elseif (l == Disc) then
               graf_name = 'graf_Disc-A000    '
               text_size = 11
               graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
               if (stop_crit == 'L') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
               elseif (stop_crit == 'G') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
               elseif (stop_crit == 'N') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
               endif
               open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
               write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Disc, 1)**0.5
               if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '

               graf_name = 'graf_Disc-A000    '
               text_size = 11
               graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
               if (stop_crit == 'L') then
                  graf_name(num_size+text_size+1:num_size+text_size+9) = 'L-adapit'
               elseif (stop_crit == 'G') then
                  graf_name(num_size+text_size+1:num_size+text_size+9) = 'G-adapit'
               elseif (stop_crit == 'N') then
                  graf_name(num_size+text_size+1:num_size+text_size+9) = 'N-adapit'
               endif
               open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
               write(20+state%time%iter,'(i7, es14.6)') state%space%adapt%adapt_level, sum_estim(Disc, 1)**0.5
               if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '

            elseif (l == Alg) then
               graf_name = 'graf_Alg-A000    '
               text_size = 10
               graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
               if (stop_crit == 'L') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
               elseif (stop_crit == 'G') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
               elseif (stop_crit == 'N') then
                  graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
               endif
               open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
               write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, sum_estim(Alg, 1)**0.5
               if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '

               graf_name = 'graf_Alg-A000    '
               text_size = 10
               graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
               if (stop_crit == 'L') then
                  graf_name(num_size+text_size+1:num_size+text_size+9) = 'L-adapit'
               elseif (stop_crit == 'G') then
                  graf_name(num_size+text_size+1:num_size+text_size+9) = 'G-adapit'
               elseif (stop_crit == 'N') then
                  graf_name(num_size+text_size+1:num_size+text_size+9) = 'N-adapit'
               endif
               open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
               write(20+state%time%iter,'(i7, es14.6)') state%space%adapt%adapt_level, sum_estim(Alg, 1)**0.5
               if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '

            endif


           enddo

           graf_name = 'graf_ErrXnorm-A000    '
           text_size = 15
           graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
           if (stop_crit == 'L') then
              graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
           elseif (stop_crit == 'G') then
              graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
           elseif (stop_crit == 'N') then
              graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
           endif
           open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
           write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, state%errSTnorm(L2H1)**0.5
           if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '

           graf_name = 'graf_ErrXnorm-A000    '
           text_size = 15
           graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
           if (stop_crit == 'L') then
              graf_name(num_size+text_size+1:num_size+text_size+9) = 'L-adapit'
           elseif (stop_crit == 'G') then
              graf_name(num_size+text_size+1:num_size+text_size+9) = 'G-adapit'
           elseif (stop_crit == 'N') then
              graf_name(num_size+text_size+1:num_size+text_size+9) = 'N-adapit'
           endif
           open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
           write(20+state%time%iter,'(i7, es14.6)') state%space%adapt%adapt_level, state%errSTnorm(L2H1)**0.5
           if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '

           graf_name = 'graf_Est-A000    '
           text_size = 10
           graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
           if (stop_crit == 'L') then
              graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
           elseif (stop_crit == 'G') then
              graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
           elseif (stop_crit == 'N') then
              graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
           endif
           open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
           write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, state%estim(Tot, 1:ndim)
           if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '

           graf_name = 'graf_Est-A000    '
           text_size = 10
           graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
           if (stop_crit == 'L') then
              graf_name(num_size+text_size+1:num_size+text_size+9) = 'L-adapit'
           elseif (stop_crit == 'G') then
              graf_name(num_size+text_size+1:num_size+text_size+9) = 'G-adapit'
           elseif (stop_crit == 'N') then
              graf_name(num_size+text_size+1:num_size+text_size+9) = 'N-adapit'
           endif
           open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
           write(20+state%time%iter,'(i7, es14.6)') state%space%adapt%adapt_level, state%estim(Tot, 1:ndim)
           if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '


           graf_name = 'graf_IE-A000    '
           text_size = 9
           graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
           if (stop_crit == 'L') then
              graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
           elseif (stop_crit == 'G') then
              graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
           elseif (stop_crit == 'N') then
              graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
           endif
           open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
           write(20+state%time%iter,'(i7, es14.6)') iter_est_loc, state%estim(Tot, 1:ndim)/(state%errSTnorm(L2H1)**0.5)
           if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '


           graf_name = 'graf_ErrXnormdof-A000    '
           text_size = 18
           graf_name(num_size+text_size-is:num_size+text_size) = ch3(num_size-is: num_size)
           if (stop_crit == 'L') then
              graf_name(num_size+text_size+1:num_size+text_size+2) = 'L'
           elseif (stop_crit == 'G') then
              graf_name(num_size+text_size+1:num_size+text_size+2) = 'G'
           elseif (stop_crit == 'N') then
              graf_name(num_size+text_size+1:num_size+text_size+2) = 'N'
           endif
           open(20+state%time%iter, file=graf_name, status='UNKNOWN', position = 'append')
           write(20+state%time%iter,'(e14.3, es14.6)') (ndim*sum(grid%elem(:)%dof))**(-0.5), state%errSTnorm(L2H1)**0.5
           if (state%space%adapt%adapt_level == state%space%adapt%max_adapt_level) write(20+state%time%iter,*) ' '


           close(20+state%time%iter)


       endif

    endif




    deallocate(estim, sum_estim)

  end subroutine ComputeAlgEstim

  !> precomputation of RTN flux momentums on elem in the first iteration step of an adaptive cycle
  subroutine InitElemRTNMomentumCD(elem, Fdeg, Fdof)
    type(element), intent(inout) :: elem
    real, dimension(:,:), allocatable :: MMinvRE ! inverse of Momentum matrix on elem
    integer, intent(in) :: Fdeg, Fdof


    if ( state%modelName == 'scalar' .or.state%modelName == '2eqs') then

       allocate(MMinvRE(1:Fdof, 1:Fdof) )
       call ComputeLocRTNMomentumsElem(elem, Fdeg, Fdof, MMinvRE)  ! RE = real element
       call ConstructFluxCD(Set_f_s_scalar, Set_Ppm_scalar, Set_R_s_scalar,Set_K_sk_scalar,&
                   elem, Fdeg, Fdof, MMinvRE, elem%RTNflux(0, 1:ndim, 1:Fdof) )

       deallocate(MMinvRE)

    else
       print*,'Case nbDim ==',nbDim,' & ndim ==',ndim,' not implemented in alg_estim.f90'
       stop
    endif


  end subroutine InitElemRTNMomentumCD


  !> precomputation of RTN flux momentums on elem
  subroutine ComputeElemRTNMomentumCD(elem, Fdeg, Fdof)
    type(element), intent(inout) :: elem
    real, dimension(:,:), allocatable :: MMinvRE ! inverse of Momentum matrix on elem
    integer, intent(in) :: Fdeg, Fdof
    integer :: ib


    allocate(MMinvRE(1:Fdof, 1:Fdof) )

    ! evaluate the momentums of the local RTN basis on K
    !call ComputeRTNMomentumsRealElem(elem, Fdeg, Fdof, MMinvRE)  ! RE = real element
    call ComputeLocRTNMomentumsElem(elem, Fdeg, Fdof, MMinvRE)  ! RE = real element

       !if (elem%i > 5 .and. elem%i < 9) then
       !   write(*,*) 'elem%i, elem%HGnode =', elem%i, elem%HGnode
       !   print*, 'v------------After ComputeLocRTNMomentumsElem sub.----------------v'
       !   write(*,'(a19)') 'MMinvRE(ib, it-vse):'
       !   do ib= 1, Fdof
       !     write(*,'(8es14.6)') MMinvRE(ib, 1:Fdof)
       !   enddo
       !   pause
       !endif


    ! flux at k-th time level
    if( state%modelName == 'scalar' .or.state%modelName == '2eqs') then

       ! allocation of arays for storing of RTN fluxes, RTNflux(l,*,*),at t_{k-l} level
       if(state%time%iter == 0) then

          !allocate(elem%RTNflux(0:1, 1:ndim, 1:Fdof) )  DONE in sub. PrepareOneElement in problem.f90

          state%time%ctime = 0.
          call ConstructFluxCD(Set_f_s_scalar, Set_Ppm_scalar, Set_R_s_scalar,Set_K_sk_scalar,&
               elem, Fdeg, Fdof, MMinvRE, elem%RTNflux(0, 1:ndim, 1:Fdof) )

          deallocate(MMinvRE)
          return

          !else if (state%time%iter > 1 .and. state%time%iter_loc == 1) then !first iteration step in an adaptation

          !   state%time%ctime = state%time%ttime

             !if (elem%i == 1 .and. state%time%iter_loc == 1) then
             !   print*,'B state%space%m_IPG ==',state%space%m_IPG
             !   pause
             !endif

          !   call ConstructFluxCD(Set_f_s_scalar, Set_Ppm_scalar, Set_R_s_scalar,Set_K_sk_scalar,&
          !        elem, Fdeg, Fdof, MMinvRE, elem%RTNflux(0, 1:ndim, 1:Fdof) )

          !   deallocate(MMinvRE)
          !   return

       else
             !write(*, '(a27, 10es14.6)') 'elem%RTNflux(1, 1, 1:Fdof)', elem%RTNflux(1, 1, 1:Fdof)
             !write(*, '(a27, 10es14.6)') 'elem%RTNflux(0, 1, 1:Fdof)', elem%RTNflux(0, 1, 1:Fdof)

             !if (elem%i == 1 .and. state%time%iter_loc == 1) then
             !    print*,'C state%space%m_IPG ==',state%space%m_IPG
             !    pause
             !endif


             !if (state%time%iter == 1) write(*,*) state%stop_rem

             ! storing of the older values
          if (stop_crit == 'L') then
             if (state%stop_rem) then      !else during increasing \nu previous information remains unchanged
                elem%RTNflux(1, 1:ndim, 1:Fdof) = elem%RTNflux(0, 1:ndim, 1:Fdof)
             endif
          elseif (stop_crit == 'G') then
             if (state%stop_rem_glob) then  !else during increasing \nu previous information remains unchanged
                elem%RTNflux(1, 1:ndim, 1:Fdof) = elem%RTNflux(0, 1:ndim, 1:Fdof)
             endif
          elseif (stop_crit == 'N') then
             if (state%stop_rem) then
                elem%RTNflux(1, 1:ndim, 1:Fdof) = elem%RTNflux(0, 1:ndim, 1:Fdof)
             endif
          else
             Print*, 'Only three possibilities of stopping criterion: L (local), G (global), and N (no stop. crit.).'
          endif


             !!!state%stop_rem = .true. !initialization CANNOT BE HERE - older values only for ONE element has been stored!

             state%time%ctime = state%time%ttime
             ! old version
             !call ConstructFluxCD(Set_f_s_scalar, Set_Ppm_scalar, Set_R_s_scalar,Set_K_sk_scalar,&
             !     elem, Fdeg, Fdof, MMinvRE, fluxi(1:ndim, 1:Fdof))

             ! new version
             call ConstructFluxCD(Set_f_s_scalar, Set_Ppm_scalar, Set_R_s_scalar,Set_K_sk_scalar,&
                  elem, Fdeg, Fdof, MMinvRE, elem%RTNflux(0, 1:ndim, 1:Fdof) )


             deallocate(MMinvRE)

          endif

       else
          print*,'Case nbDim ==',nbDim,' & ndim ==',ndim,' not implemented in alg_estim.f90'
          stop
       endif

       !if(elem%i <= 1) then
       !   print*,state%stop_rem
       !   !write(*,'(a4, i5,20es12.4)')'w0',elem%i, elem%w(0,: )
       !   !write(*,'(a4, i5,20es12.4)')'w1',elem%i, elem%w(1,: )
       !   !write(*,'(a4, i5,20es12.4)')'w01',elem%i, elem%w(0,: )- elem%w(1,: )
       !   !write(199,'(a4, i5,20es12.4)')'w01',elem%i, elem%w(0,: )- elem%w(1,: )
       !   !print*,'------------------------'
       !   write(*,'(a4, i5,20es12.4)')'RT1',elem%i, elem%RTNflux(1, 1,:)
       !   write(*,'(a4, i5,20es12.4)')'RT0',elem%i, elem%RTNflux(0, 1,:)
       !   print*
       !endif

  end subroutine ComputeElemRTNMomentumCD

  !> compute ALGebraic error included estimates in energy norm using RTN flux reconstruction
  !> for one element
  subroutine ComputeElemAlgEstim(elem, Fdeg, Fdof, estim )
    type(element), intent(inout) :: elem
    real, dimension(1:max_eta, 1:ndim), intent(inout) :: estim
    real, dimension(:,:), allocatable :: fluxi   ! basis coefficients of flux
    real, dimension(:), allocatable :: fluxmomentum  !for testing
    !real, dimension(:,:,:), allocatable :: CDfluxes   ! physical fluxes
    integer, intent(in) :: Fdeg, Fdof

    ! scaling factor computed in errorFlux.f90
    !elem%CTn = state%time%FinTime/state%time%tau(1)**2 + (elem%CK + elem%Cb)/elem%diam**2


    allocate(fluxi(1:ndim, 1:Fdof) )
    ! allocate(CDfluxes(1:elem%Qdof, 1:nbDim, 1:ndim)  )

    ! flux at k-th time level
    if(state%modelName == 'scalar' .or.state%modelName == '2eqs') then

       ! allocation of arays for storing of RTN fluxes, RTNflux(l,*,*),at t_{k-l} level
       if(state%time%iter > 0) then
         !if (state%time%iter == 1 .or. state%time%iter_loc /= 1) then ! first iteration step or any iteration step in an adaptation except for the first one

          state%time%ctime = state%time%ttime


          call PotentNCElemEstimSS(Set_R_s_scalar, elem, estim(PNC, 1:ndim) )

          if (state%model%icase == 35 .or. state%model%icase == 43) then
             call DataOscillationAlgEstimSS(elem, Fdeg, Fdof, estim(Resid, 1:ndim))
                                        else
             call ResidualElemAlgEstimCDSS(elem, Fdeg, Fdof, fluxi, estim(Resid, 1:ndim))
          endif

          elem%eta(Resid, 1:ndim) = estim(Resid, 1:ndim)    !saving of values of elem residual

          call FluxElemDiscEstimCDSS(Set_f_s_scalar, Set_R_s_scalar, &
               elem, Fdeg, Fdof, fluxi,  estim(FD, 1:ndim))

          !if(elem%i<=5)write(*,'(a2, i5,20es12.4)') &
          !     'Be',elem%i,estim(FA,1:ndim), estim(:,1)
          call FluxElemAlgEstimCDSS(elem, Fdeg, Fdof, estim(FA, 1:ndim))

          if(elem%i<= -2) then
             !write(*,'(a2, i5,20es12.4)') &
             !  'RT',elem%i, elem%RTNflux(1, 1,:) - elem%RTNflux(0, 1,:)
             write(*,'(a2, i5,20es12.4)') &
               'Af',elem%i,estim(FA,1:ndim), estim(:,1)
          endif

          call FluxElemTotEstimCDSS(Set_f_s_scalar, Set_R_s_scalar, &
               elem, Fdeg, Fdof, estim(FT, 1:ndim))

          !if (state%space%adapt%adapt_method == 'ALG' .and. stop_crit == 'N') then
          !   estim(FA, 1:ndim) = 0
          !endif

          !call FluxNCElemDiscEstimSS(elem, Fdeg, Fdof, .true., estim(FNCD, 1:ndim))  OLD VERSION 1

          !call FluxNCElemDiscEstimSS(elem, Fdeg, Fdof, .false., estim(FNCA, 1:ndim))  OLD VERSION 1

          !estim(FNCD, 1:ndim) = elem%estimFNCD(1:ndim) OLD VERSION 2

          !estim(FNCA, 1:ndim) = elem%estimFNCA(1:ndim) OLD VERSION 2

          call FluxNCElemEstimSSnew(elem, Fdeg, Fdof, 'dis')

          call FluxNCElemEstimSSnew(elem, Fdeg, Fdof, 'alg')

          call FluxNCElemEstimSSnew(elem, Fdeg, Fdof, 'tot')

          estim(FNCD, 1:ndim) = elem%estimFNCD(1:ndim)

          estim(FNCA, 1:ndim) = elem%estimFNCA(1:ndim)

          estim(FNCT, 1:ndim) = elem%estimFNCT(1:ndim)

          !if (state%space%adapt%adapt_method == 'ALG' .and. stop_crit == 'N') then
          !   estim(FNCA, 1:ndim) = 0
          !endif


          estim(Disc, 1:ndim) = estim(PNC, 1:ndim) + estim(Resid, 1:ndim) &
                             + estim(FNCD, 1:ndim) + estim(FD, 1:ndim)

          estim(Alg, 1:ndim) = estim(FA, 1:ndim) + estim(FNCA, 1:ndim)


          call RemainderElemAlgEstimSS(elem, estim(Rem, 1:ndim))


          elem%eta(Disc, 1:ndim) = estim(Disc, 1:ndim)

          elem%eta(Alg, 1:ndim) = estim(Alg, 1:ndim)

          elem%eta(Rem, 1:ndim) = estim(Rem, 1:ndim)

          elem%eta(Tot, 1:ndim) = estim(Disc, 1:ndim) + estim(Alg, 1:ndim) + estim(Rem, 1:ndim)

         !endif
       endif

    else
       print*,'Case nbDim ==',nbDim,' & ndim ==',ndim,' not implemented in alg_estim.f90'
       stop
    endif


    ! computed in errorFlux.f90
    !call NonConfEstimCD(elem, estim(NC1n,1:ndim) )
    !estim(NC1n,1:ndim) = elem%eta(NC1n, 1:ndim)
    !estim(NC1n,1:ndim) = 0.   ! completly done in errorF.f90

    ! we sum the same values for each time step
    !time dep: estim(IC,1:ndim) = elem%eta(IC,1)**2 / (elem%CTn*state%time%tau(1) )

    ! \eta_space include {\not nonconformity term} and redidual term!!!
    !estim(DFnS,1:ndim) = ( estim(DFnS,1:ndim)**0.5 + estim(Rn, 1:ndim)**0.5 &
    !     + elem%eta(NC1n, 1:ndim)**0.5 )**2
    !estim(DFnS,1:ndim) = ( estim(DFnS,1:ndim)**0.5 + estim(Rn, 1:ndim)**0.5 )**2

    ! for adaptation
    !elem%eta(DFn, 1:ndim)  = estim(DFn, 1:ndim)
    !elem%eta(DFnS, 1:ndim)  = estim(DFnS, 1:ndim)
    !elem%eta(DFnT, 1:ndim)  = estim(DFnT, 1:ndim)
    !elem%eta(Rn, 1:ndim)  = estim(Rn, 1:ndim)

    !elem%estimST = elem%eta(DFn, 1)**0.5 + elem%eta(Rn, 1)**0.5  + elem%eta(NC1n, 1)**0.5


    !write(*,'(a6,i5,5es12.4)') &
    !     '######',elem%i, elem%estim, elem%eta(DFn, 1), elem%eta(Rn, 1), elem%eta(NC1n, 1)

    !write(1000+state%time%iter,*) elem%xc(:),estim(DFn, 1:ndim)**0.5, estim(Rn, 1:ndim)**0.5
    !write(*,'(a12,i5,6es12.4)') 'eta DF, R :', elem%i, &
    !     estim(DFn, 1:ndim)**0.5, estim(Rn, 1:ndim)**0.5, estim(NC1n, 1:ndim)**0.5

    deallocate(fluxi) !,CDfluxes )

  end subroutine ComputeElemAlgEstim








  !> evaluate resiual estimator
  !> \f$ \eta_{R,T}^i := c_P,T h_T \|f -\nabla\cdot {\bf t}^{i}_{h} - r_h^{i+nu}\|_{T} \f$
  !Time dep: \f$ (c_P h_T m_T \|f -\partial_t u_{h\tau} -\nabla\cdot {\bf t}_{h\tau}\|_{T\times I_m} )^2\f$
  subroutine ResidualElemAlgEstimCDSS(elem, Fdeg, Fdof, fluxi, etaRnT)
    type(element), intent(in) :: elem
    integer, intent(in) :: Fdeg, Fdof
    real, dimension(1:ndim, 1:Fdof), intent(in) :: fluxi
    real, dimension(1:ndim), intent(inout) :: etaRnT
    type(basis_rtn_fe), pointer :: loc_RTN
    type(Gauss_rule), pointer :: G_rule
    real, dimension(:,:), allocatable :: DivPsi ! divergence of local RTN basis functions
    real, dimension(:,:), allocatable :: wt     ! time derivative of the approximate sol
    real, dimension(:,:), allocatable :: RTNflux_loc  ! flux reconstruction pw linear in time
    real, dimension(:,:), allocatable :: temp, temp1, x, Fx, f
    real, dimension(:,:), pointer :: phi
    real, dimension(1:ndim) :: val, val2
    real :: t, mT
    integer :: k, l, it, dof, Qdof, Qnum, Gnum, Gdof

    dof = elem%dof
    Qdof = elem%Qdof
    Qnum = elem%Qnum

    loc_RTN => state%loc_RTN(Fdeg)

    ! righ hand side (source terms)
    allocate(x(1:Qdof, 1:nbDim))
    x(1:Qdof, 1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)

    allocate(Fx(1:Qdof, 1:nbDim))
    allocate(f(1:ndim, 1:Qdof) )
    allocate(RTNflux_loc(1:ndim, 1:Fdof) )

    ! Fx contains physical coordinates of integ nodes
    call ComputeF(elem, Qdof, x, Fx)

    ! time derivative of the approximate solution
    allocate( wt(1:Qdof, 1:ndim) )
    call Eval_w_t_Elem(elem, wt)


    ! divergence of the RTN basis functions on elem in integ nodes
    allocate(DivPsi(1:Fdof, 1:Qdof) )
    call Eval_Loc_RTN_Div_Elem(elem, loc_RTN, DivPsi )

    allocate( temp(1:ndim, 1:Qdof), temp1(1:ndim, 1:Qdof) )
    ! time independent
!    do k=1, ndim
!       temp(k, 1:Qdof) = - matmul(fluxi(k, 1:Fdof), DivPsi(1:Fdof, 1:Qdof) )
!
!       temp(k, 1:Qdof) = temp(k, 1:Qdof)  - wt(1:Qdof, k)
!
!       !if(elem%i == 2) &
!       !call PlotElemFunctionQ(740+state%time%iter, elem, 'V', Qnum, Qdof,  temp(k,1:Qdof) )
!
!    enddo


    Gnum = 1   ! code for time dependent problem adapted for stationary one
    !Time dep: Gnum = 2   ! 15 is the maximal one !!!!!!!!!
    G_rule => state%space%G_rule(Gnum)
    Gdof = G_rule%Qdof

    etaRnT(1:ndim) = 0.
    ! integration over Gauss integ nodes in time
    do it=1, Gdof
       !Time dep: t = G_rule%lambda(it)
       t = 0. !SS
       state%time%ctime = state%time%ttime - state%time%tau(1) * t

       ! time derivative and the divergence of the flux
       !if(elem%i < 2) then
       !   write(*,'(a8,i5,30es12.4)') '..old',elem%i, fluxi(1, 1:Fdof)
       !   write(*,'(a8,i5,30es12.4)') '..new',elem%i, elem%RTNflux(0, 1, 1:Fdof)
       !   write(*,'(a8,i5,30es12.4)') '..new',elem%i, elem%RTNflux(1, 1, 1:Fdof)
       !   write(*,'(a8,i5,30es12.4)') 'diff ',elem%i, fluxi(1, 1:Fdof)-elem%RTNflux(0, 1, 1:Fdof)
        !  write(*,'(a8,i5,30es12.4)') 'diff2 ',elem%i,elem%RTNflux(1, 1, 1:Fdof) -elem%RTNflux(0, 1, 1:Fdof)
        !  write(*,*)'---------------------------------------------------'
       !endif

       !write(*,'(a2, es12.4)') 't=', t

       ! linear interpolation
       !! t = 0.  ! PIECEWISE CONSTANT RECONSTRUCTION
       RTNflux_loc(1:ndim, 1:Fdof) = t * elem%RTNflux(1, 1:ndim, 1:Fdof) &
            + ( 1. - t)  *  elem%RTNflux(0, 1:ndim, 1:Fdof)

       do k=1, ndim
          temp(k, 1:Qdof) = - matmul(RTNflux_loc(k, 1:Fdof), DivPsi(1:Fdof, 1:Qdof) )

          !Time dep: temp(k, 1:Qdof) = temp(k, 1:Qdof)  - wt(1:Qdof, k)
       enddo

       ! source (RHS) term
       do l=1, Qdof
          call RHS_Scalar(Fx(l, 1:nbDim), f(1:ndim, l), state%time%ctime)
       enddo


       temp1(1:ndim, 1:Qdof) = temp(1:ndim, 1:Qdof) + f(1:ndim, 1:Qdof) - elem%res_func(1:ndim, 1:Qdof)

       !SS
       !if (elem%i < 21) then
          !write(*,*) 'elem%xc = ', elem%xc
          !write(*,*) 'elem%Qnum = ', elem%Qnum
          !write(*,*) 'elem%Qdof = ', elem%Qdof
          !do l=1, Qdof
             !write(*,'(a2, i2, a1, es12.4)') 'f(', l, ')', f(1:ndim, l)
             !write(*,'(a8, i2, a1, es12.4)') '-divthi(', l, ')', temp(1:ndim, l)
             !write(*,'(a10, i2, a1, es12.4)') 'f-divthi(', l, ')', f(1:ndim, l) + temp(1:ndim, l)
            !write(*,'(a20, 17es12.4)') 'elem%res_func( 1:Qdof )', elem%res_func(1:ndim, 1:Qdof)
             !write(*,'(a20, 17es12.4)') 'f-divthi-elem%res_func(1:Qdof)', temp1(1:ndim, 1:Qdof)
             !print*, ' '
          !enddo
          !pause
       !endif

       !temp1(1:ndim, 1:Qdof) = 0.
       !do l=1,nbDim
       !   do k =1, Qdof
       !     temp1(1:ndim, k) = temp1(1:ndim, k) + (temp(1:ndim, k)*x(k, l) )**2
       !   enddo !k
       !enddo !l

       !call IntegrateVectorFunction(elem, temp1(1:ndim, 1:Qdof), val(1:ndim))
       !write(*,'(a18, es12.4)') 'aver -divthi * x', val(1:ndim)

       !temp1(1:ndim, 1:Qdof) = 0.
       !do l=1,nbDim
       !   do k =1, Qdof
       !     temp1(1:ndim, k) = temp1(1:ndim, k) + (f(1:ndim, k)*x(k, l) )**2
       !   enddo !k
       !enddo !l

       !call IntegrateVectorFunction(elem, temp1(1:ndim, 1:Qdof), val2(1:ndim))
       !write(*,'(a18, es12.4)') 'aver f * x', val2(1:ndim)
       !pause

       !temp1(1:ndim, 1:Qdof) = 0.
       !do l=1,nbDim
       !   do k =1, Qdof
       !     temp1(1:ndim, k) = temp1(1:ndim, k) + (temp(1:ndim, k)*(x(k, l)**2) )**2
       !   enddo !k
       !enddo !l

       !call IntegrateVectorFunction(elem, temp1(1:ndim, 1:Qdof), val(1:ndim))
       !write(*,'(a18, es12.4)') 'aver -divthi * x^2', val(1:ndim)

       !temp1(1:ndim, 1:Qdof) = 0.
       !do l=1,nbDim
       !   do k =1, Qdof
       !     temp1(1:ndim, k) = temp1(1:ndim, k) + (f(1:ndim, k)*(x(k, l)**2) )**2
       !   enddo !k
       !enddo !l

       !call IntegrateVectorFunction(elem, temp1(1:ndim, 1:Qdof), val2(1:ndim))
       !write(*,'(a18, es12.4)') 'aver f * x^2', val2(1:ndim)
       !pause


       !call IntegrateVectorFunction(elem, temp(1:ndim, 1:Qdof), val(1:ndim))
       !write(*,'(a15, es12.4)') 'aver -divthi', val(1:ndim)
       !call IntegrateVectorFunction(elem, f(1:ndim, 1:Qdof), val2(1:ndim))
       !write(*,'(a10, es12.4)') 'aver f', val2(1:ndim)
       !write(*,'(a26, es12.4)') 'aver -divthi + aver f', val(1:ndim) + val2(1:ndim)
       !pause



       !call PlotElemFunctionQ(840+state%time%iter, elem, 'V', Qnum, Qdof,  temp(1,1:Qdof) )


       call IntegrateSquareVectorFunction(elem, temp1(1:ndim, 1:Qdof), val(1:ndim) )

       !Time dep: etaRnT(1:ndim) = etaRnT(1:ndim) + G_rule%weights(it) * val(1:ndim)    ! summing eta_disc over T not eta_R,T


       !if(elem%i == grid%nelem) write(100,'(a8,i5, 6es12.4)') 'val :',it, val(:)
    enddo ! it =1,Gdof ! time integration


    !Time dep: etaRnT(1:ndim) = etaRnT(1:ndim) * state%time%tau(1)    ! length of the time interval

    ! old version May 2011
    !mT = state%time%tau(1)/(state%time%FinTime**0.5 * elem%diam)
    !if(elem%CK >0) mT = min(mT, 1./elem%CK**0.5)
    !if(elem%Cb >0) mT = min(mT, 1./elem%Cb**0.5)
    !mT = (mT* elem%diam/ Pi)**2
    !etaRnT(1:ndim) = etaRnT(1:ndim) * mT


    ! new version June 2011
    !etaRnT(1:ndim) = etaRnT(1:ndim) / (Pi * Pi * elem%CTn)  ! etaRnt stores \eta_R^2 !!!
    ! new version October 2011 Poincare -> Friedrichs
    !Time dep: etaRnT(1:ndim) = etaRnT(1:ndim) / (2. * elem%CTn)  ! etaRnt stores \eta_R^2 !!!

    etaRnT(1:ndim) = elem%diam * elem%CP * (val(1:ndim)**0.5)   !elem%CP defined in problem.f90

    !!write(*,*) 'elem%i = ', elem%i
    !!write(*,*) 'elem%xc = ', elem%xc
    !!write(*,'(a8, es12.4)') 'etaRnT =', etaRnT(1:ndim)
    !!pause

    !write(*,'(a6,i5,12es12.4)') '#Rn#',elem%i, &
    !     1./(Pi * elem%CTn**0.5), mT,1./(Pi * elem%CTn**0.5)-  mT, &
    !     state%time%tau(1)/(state%time%FinTime**0.5 * elem%diam), 1./elem%CK**0.5,1./elem%Cb**0.5

    !if(elem%i == grid%nelem) write(100,'(a8,i5, 6es12.4)') 'etaRn:',it, etaRnT(:), &
    !     state%time%tau(1), mT, (state%time%tau(1)/Pi)**2,elem%CK, elem%Cb


    deallocate(DivPsi, temp, temp1, x, Fx, f, RTNflux_loc)

    !write(*,'(a6,i5,12es14.6)') 'etaRn:',elem%i, etaRnT(1:ndim), mT

  end subroutine ResidualElemAlgEstimCDSS


  subroutine ElementPolynProjstate(elem, f, Pf)
   type(element), intent (inout) :: elem     ! element where the basis is given
   real, dimension(1:elem%Qdof, 1:nbDim, 1:ndim), intent(in) :: f
   real, dimension(1:elem%Qdof, 1:nbDim, 1:ndim), intent(out) :: Pf
   real, dimension(:,:), pointer :: phi
   real, dimension(:), allocatable :: wi, qi
   integer :: dof,  k, j, l, Qnum, Qdof, deg

   deg = state%space%deg
   dof = (deg+1)*(deg+2)/2
   Qnum = elem%Qnum
   Qdof = state%space%V_rule(Qnum)%Qdof

   phi => state%space%V_rule(Qnum)%phi(1:dof,1:Qdof)

   allocate( wi(1:dof) )
   allocate( qi(1:dof) )

   ! evaluating of the coefficients of basis expansion
   do k=1, ndim
      do j=1,nbDim
         qi(:) = 0.
         call IntegrateVectorB(elem, dof, f(1:Qdof, j, k), qi(1:dof) )

         ! wi are basis coefficients in integ nodes
         do l=1,dof
            wi(l) = dot_product(elem%MassInv%Mb(l,1:dof), qi(1:dof) )
         enddo

         do l=1,Qdof
            Pf(l, j, k) = dot_product(wi(1:dof), phi(1:dof, l) )
         enddo
      enddo
   enddo
   deallocate(wi, qi)

 end subroutine ElementPolynProjstate

  !> evaluate Data oscillation
  !> \f$ \eta_{R,T}^i := c_P,T h_T \|f - \Pi_f\|_{T} \f$
  subroutine DataOscillationAlgEstimSS(elem, Fdeg, Fdof, etaRnT)
    type(element), intent(inout) :: elem
    integer, intent(in) :: Fdeg, Fdof
    real, dimension(1:ndim), intent(inout) :: etaRnT
    type(Gauss_rule), pointer :: G_rule
    real, dimension(:,:), allocatable :: temp, temp1, x, Fx, f
    real, dimension(:,:,:), allocatable :: f_vec, Pf_vec
    real, dimension(:,:), pointer :: phi
    real, dimension(1:ndim) :: val, val2
    real :: t, mT
    integer :: k, l, it, dof, Qdof, Qnum, Gnum, Gdof

    real, dimension(:), allocatable :: qi

    !dof = elem%dof
    Qdof = elem%Qdof
    Qnum = elem%Qnum


    ! righ hand side (source terms)
    allocate(x(1:Qdof, 1:nbDim))
    x(1:Qdof, 1:nbDim) = state%space%V_rule(Qnum)%lambda(1:Qdof,1:nbDim)

    allocate(Fx(1:Qdof, 1:nbDim))
    allocate(f(1:ndim, 1:Qdof) )
    allocate(f_vec(1:Qdof, 1:nbDim, 1:ndim) )
    allocate(Pf_vec(1:Qdof, 1:nbDim, 1:ndim) )

    allocate( qi(1:dof) )

    ! Fx contains physical coordinates of integ nodes
    call ComputeF(elem, Qdof, x, Fx)



    allocate( temp1(1:ndim, 1:Qdof) )
    ! time independent
!    do k=1, ndim
!       temp(k, 1:Qdof) = - matmul(fluxi(k, 1:Fdof), DivPsi(1:Fdof, 1:Qdof) )
!
!       temp(k, 1:Qdof) = temp(k, 1:Qdof)  - wt(1:Qdof, k)
!
!       !if(elem%i == 2) &
!       !call PlotElemFunctionQ(740+state%time%iter, elem, 'V', Qnum, Qdof,  temp(k,1:Qdof) )
!
!    enddo


    Gnum = 1   ! code for time dependent problem adapted for stationary one
    !Time dep: Gnum = 2   ! 15 is the maximal one !!!!!!!!!
    G_rule => state%space%G_rule(Gnum)
    Gdof = G_rule%Qdof

    etaRnT(1:ndim) = 0.
    ! integration over Gauss integ nodes in time
    do it=1, Gdof
       !Time dep: t = G_rule%lambda(it)
       t = 0. !SS
       state%time%ctime = state%time%ttime - state%time%tau(1) * t

       ! time derivative and the divergence of the flux
       !if(elem%i < 2) then
       !   write(*,'(a8,i5,30es12.4)') '..old',elem%i, fluxi(1, 1:Fdof)
       !   write(*,'(a8,i5,30es12.4)') '..new',elem%i, elem%RTNflux(0, 1, 1:Fdof)
       !   write(*,'(a8,i5,30es12.4)') '..new',elem%i, elem%RTNflux(1, 1, 1:Fdof)
       !   write(*,'(a8,i5,30es12.4)') 'diff ',elem%i, fluxi(1, 1:Fdof)-elem%RTNflux(0, 1, 1:Fdof)
        !  write(*,'(a8,i5,30es12.4)') 'diff2 ',elem%i,elem%RTNflux(1, 1, 1:Fdof) -elem%RTNflux(0, 1, 1:Fdof)
        !  write(*,*)'---------------------------------------------------'
       !endif

       !write(*,'(a2, es12.4)') 't=', t



       ! source (RHS) term
       do l=1, Qdof
          call RHS_Scalar(Fx(l, 1:nbDim), f(1:ndim, l), state%time%ctime)
       enddo

       do k=1, ndim
         f_vec(1:Qdof, 1, k) = f(k, 1:Qdof)
       enddo

       do k=2, nbDim
          f_vec(1:Qdof, k, 1:ndim) = 0.
       enddo

       if (state%space%deg /= elem%deg) write(*,*) 'Problem: state%space%deg /= elem%deg. ElementPolynProj cannot be used!!'

       call ElementPolynProjstate(elem, f_vec, Pf_vec)


       do k=1, ndim
         temp1(k, 1:Qdof) = f(k, 1:Qdof) - Pf_vec(1:Qdof, 1, k)
       enddo


       !TEST
       !qi(:) = 0.
       !call IntegrateVectorB(elem, dof, f(1, 1:Qdof), qi(1:dof) )
       !write(*,*) 'f * phi = ', qi(1:dof)
       !qi(:) = 0.
       !call IntegrateVectorB(elem, dof, Pf_vec(1:Qdof, 1, 1), qi(1:dof) )
       !write(*,*) 'Pf * phi = ', qi(1:dof)
       !pause



       !temp1(1:ndim, 1:Qdof) = 0.
       !do l=1,nbDim
       !   do k =1, Qdof
       !     temp1(1:ndim, k) = temp1(1:ndim, k) + (temp(1:ndim, k)*x(k, l) )**2
       !   enddo !k
       !enddo !l

       !call IntegrateVectorFunction(elem, temp1(1:ndim, 1:Qdof), val(1:ndim))
       !write(*,'(a18, es12.4)') 'aver -divthi * x', val(1:ndim)

       !temp1(1:ndim, 1:Qdof) = 0.
       !do l=1,nbDim
       !   do k =1, Qdof
       !    temp1(1:ndim, k) = temp1(1:ndim, k) + (f(1:ndim, k)*x(k, l) )**2
       !   enddo !k
       !enddo !l

       !call IntegrateVectorFunction(elem, temp1(1:ndim, 1:Qdof), val2(1:ndim))
       !write(*,'(a18, es12.4)') 'aver f * x', val2(1:ndim)
       !pause

       !temp1(1:ndim, 1:Qdof) = 0.
       !do l=1,nbDim
       !   do k =1, Qdof
       !     temp1(1:ndim, k) = temp1(1:ndim, k) + (temp(1:ndim, k)*(x(k, l)**2) )**2
       !   enddo !k
       !enddo !l

       !call IntegrateVectorFunction(elem, temp1(1:ndim, 1:Qdof), val(1:ndim))
       !write(*,'(a18, es12.4)') 'aver -divthi * x^2', val(1:ndim)

       !temp1(1:ndim, 1:Qdof) = 0.
       !do l=1,nbDim
       !   do k =1, Qdof
       !     temp1(1:ndim, k) = temp1(1:ndim, k) + (f(1:ndim, k)*(x(k, l)**2) )**2
       !   enddo !k
       !enddo !l

       !call IntegrateVectorFunction(elem, temp1(1:ndim, 1:Qdof), val2(1:ndim))
       !write(*,'(a18, es12.4)') 'aver f * x^2', val2(1:ndim)
       !pause


       !call IntegrateVectorFunction(elem, temp(1:ndim, 1:Qdof), val(1:ndim))
       !write(*,'(a15, es12.4)') 'aver -divthi', val(1:ndim)
       !call IntegrateVectorFunction(elem, f(1:ndim, 1:Qdof), val2(1:ndim))
       !write(*,'(a10, es12.4)') 'aver f', val2(1:ndim)
       !write(*,'(a26, es12.4)') 'aver -divthi + aver f', val(1:ndim) + val2(1:ndim)
       !pause



       !call PlotElemFunctionQ(840+state%time%iter, elem, 'V', Qnum, Qdof,  temp(1,1:Qdof) )


       call IntegrateSquareVectorFunction(elem, temp1(1:ndim, 1:Qdof), val(1:ndim) )

       !Time dep: etaRnT(1:ndim) = etaRnT(1:ndim) + G_rule%weights(it) * val(1:ndim)    ! summing eta_disc over T not eta_R,T


       !if(elem%i == grid%nelem) write(100,'(a8,i5, 6es12.4)') 'val :',it, val(:)
    enddo ! it =1,Gdof ! time integration


    !Time dep: etaRnT(1:ndim) = etaRnT(1:ndim) * state%time%tau(1)    ! length of the time interval



    etaRnT(1:ndim) = elem%diam * elem%CP * (val(1:ndim)**0.5)   !elem%CP defined in problem.f90

    !SS
    if (etaRnT(1) < 0.) then
      write(*,*) 'elem%i = ', elem%i
      write(*,*) 'elem%deg = ', elem%deg
       write(*,*) 'elem%xc = ', elem%xc
       write(*,*) 'elem%Qnum = ', elem%Qnum
       write(*,*) 'elem%Qdof = ', elem%Qdof
       write(*,'(a8, es12.4)') 'etaRnT =', etaRnT(1:ndim)
       print*, ' '
       !do l=1, Qdof
          !write(*,'(a2, i2, a1, es12.4)') 'f(', l, ')', f(1:ndim, l)
          !write(*,'(a8, i2, a1, es12.4)') '-divthi(', l, ')', temp(1:ndim, l)
          !write(*,'(a10, i2, a1, es12.4)') 'f-divthi(', l, ')', f(1:ndim, l) + temp(1:ndim, l)
          write(*,'(a20, 175es12.4)') 'f-\Pi_f', temp1(1:ndim, 1:Qdof)
          print*, ' '
       !enddo
       pause
    endif


    deallocate(temp1, x, Fx, f, f_vec, Pf_vec)

    deallocate( qi )

  end subroutine DataOscillationAlgEstimSS



  ! evaluate flux convective diffusive estimator
  !> evaluation of total flux (diffusive) estimator
  !> \f$ \|\nabla u_h^i +  {\bf t}_{h}^i\|_{T}\f$
  subroutine FluxElemTotEstimCDSS(Set_f_s, Set_R_s, elem, Fdeg, Fdof, etaDFn)
    interface
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
      end subroutine Set_f_s

      subroutine Set_R_s(ndimL, nbDim, Qdof, w, Dw, Re_1, R_s)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         !real, intent(in) :: Re_1                     ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
      end subroutine Set_R_s


    end interface
    type(element), intent(in) :: elem
    integer, intent(in) :: Fdeg, Fdof
    !real, dimension(1:ndim, 1:Fdof), intent(in) :: fluxi
    !!real, dimension(1:elem%Qdof, 1:nbDim, 1:ndim), intent(in) :: CDfluxes
    real, dimension(:, :, :), allocatable :: fluxes, fluxesF
    !Time dep: real, dimension(1:3, 1:ndim), intent(inout) :: etaDFn ! 1 = total, 2 = space, 3 = time
    real, dimension(1:ndim), intent(inout) :: etaDFn !SS
    type(basis_rtn_fe), pointer :: loc_RTN
    type(Gauss_rule), pointer :: G_rule
    real, dimension(:,:,:), allocatable :: psi  ! local RTN basis functions
    real, dimension(:,:), allocatable :: wt     ! time derivative of the approximate sol
    real, dimension(:,:), allocatable :: RTNflux_loc  ! flux reconstruction pw linear in time
    real, dimension(:,:), allocatable :: x, Fx, f
    real, dimension(:,:,:), allocatable :: f_s, R_s, temp
    real, dimension(1:ndim) :: val
    real :: t, mT
    integer :: it, k, l, dof, Qdof, Qnum, Gnum, Gdof

    dof = elem%dof
    Qdof = elem%Qdof
    Qnum = elem%Qnum

    loc_RTN => state%loc_RTN(Fdeg)

    ! divergence of the RTN basis functions on elem in integ nodes
    allocate(psi(1:Fdof, 1:nbDim, 1:Qdof) )
    call Eval_Loc_RTN_Elem(elem, loc_RTN, psi )

    allocate( temp(1:3, 1:ndim, 1:Qdof) )  ! first index: 1 = total, 2 = space, 3 = time
    allocate(fluxes(1:Qdof, 1:nbDim, 1:ndim), fluxesF(1:Qdof, 1:nbDim, 1:ndim) )
    allocate(R_s(1:Qdof, 1:nbDim, 1:ndim) )
    allocate(f_s(1:Qdof, 1:nbDim, 1:ndim) )

    allocate(RTNflux_loc(1:ndim, 1:Fdof) )

    ! fluxes at time level t_k
    state%time%ctime = state%time%ttime
    call Eval_R_s_Elem(Set_R_s, elem, R_s(1:Qdof, 1:nbDim, 1:ndim) )

    call Eval_f_s_Elem(Set_f_s, elem, f_s(1:Qdof, 1:nbDim, 1:ndim) )

    fluxesF(1:Qdof, 1:nbDim, 1:ndim) = f_s(1:Qdof, 1:nbDim, 1:ndim) &
            - R_s(1:Qdof, 1:nbDim, 1:ndim)

    ! setting of integ  rule
    Gnum = 1   ! code for time dependent problem adapted for stationary one
    !Time dep: Gnum = 2   ! 15 is the maximal one !!!!!!!!!
    G_rule => state%space%G_rule(Gnum)
    Gdof = G_rule%Qdof

    !Time dep: etaDFn(1:3,1:ndim) = 0.
    etaDFn(1:ndim) = 0.  !SS

    ! integration over Gauss integ nodes in time
    do it=1, Gdof
       !Time dep: t = G_rule%lambda(it)
       t = 0. !SS   evaluation of entities at the present iteration step
       state%time%ctime = state%time%ttime - state%time%tau(1) * t

       call Eval_R_s_Elem_at_time(Set_R_s, elem, R_s(1:Qdof, 1:nbDim, 1:ndim), 0. )   ! t=0 => t_{k-1} - corresponds to \nabla u_h^i
       !write(*,'(a6,2i5,20es14.6)') 'Rs=',elem%i,1,R_s(1:Qdof, 1, 1:ndim)
       !write(*,'(a6,2i5,20es14.6)') 'Rs=',elem%i,2,R_s(1:Qdof, 2, 1:ndim)

       call Eval_f_s_Elem_at_time(Set_f_s, elem, f_s(1:Qdof, 1:nbDim, 1:ndim), 0. )
       !write(*,'(a6,2i5,20es14.6)') 'fs=',elem%i,1,f_s(1:Qdof, 1, 1:ndim)
       !write(*,'(a6,2i5,20es14.6)') 'fs=',elem%i,2,f_s(1:Qdof, 2, 1:ndim)


       ! used in flux estimator
       fluxes(1:Qdof, 1:nbDim, 1:ndim) = f_s(1:Qdof, 1:nbDim, 1:ndim) &
            - R_s(1:Qdof, 1:nbDim, 1:ndim)

       ! linear interpolation
       !! t = 0.  ! PIECEWISE CONSTANT RECONSTRUCTION
       RTNflux_loc(1:ndim, 1:Fdof) = t * elem%RTNflux(1, 1:ndim, 1:Fdof) &
            + ( 1. - t)  *  elem%RTNflux(0, 1:ndim, 1:Fdof)

       do k=1, ndim
          temp(1:3, k, 1:Qdof) = 0.
          do l=1,nbDim
             ! total estimate
             temp(1, k, 1:Qdof) = temp(1, k, 1:Qdof) + &
                  (matmul(RTNflux_loc(k, 1:Fdof), psi(1:Fdof, l, 1:Qdof) ) &
                  - fluxes(1:Qdof, l, k) )**2

             ! space estimate (can be outside of the cycles)
             temp(2, k, 1:Qdof) = temp(1, k, 1:Qdof) + &
                  (matmul(RTNflux_loc(k, 1:Fdof), psi(1:Fdof, l, 1:Qdof) ) &
                  - fluxesF(1:Qdof, l, k) )**2

             ! time estimate
             temp(3, k, 1:Qdof) = temp(1, k, 1:Qdof) + &
                  (fluxesF(1:Qdof, l, k) - fluxes(1:Qdof, l, k) )**2

             !write(*,'(a8,2i5,20es12.4)') 'theta',elem%i,l, &
             !     matmul(fluxi(k, 1:Fdof), psi(1:Fdof, l, 1:Qdof) )
             !
             !write(*,'(a8,2i5,20es12.4)') 'flux   ',elem%i,l, fluxes(1:Qdof, l, k)
             !write(*,'(a8,2i5,20es12.4)') 'fluxCD ',elem%i,l, CDfluxes(1:Qdof, l, k)


          enddo
          !call PlotElemFunctionQ(840+state%time%iter, elem, 'V', Qnum, Qdof,  temp(k,1:Qdof) )

       enddo


       !Time dep: do l=1,3
          !Time dep: call IntegrateVectorFunction(elem, temp(l, 1:ndim, 1:Qdof), val(1:ndim))
          call IntegrateVectorFunction(elem, temp(1, 1:ndim, 1:Qdof), val(1:ndim))
          etaDFn(1:ndim) = val(1:ndim)**0.5 !SS
          !Time dep: etaDFn(l, 1:ndim) = etaDFn(l, 1:ndim) + G_rule%weights(it) * val(1:ndim)
       !Time dep: enddo

    enddo

    !Time dep: etaDFn(1:3, 1:ndim) = etaDFn(1:3, 1:ndim) * state%time%tau(1)    ! length of the time interval

    ! etaDFn stores \eta_DF^2 !!
    !Time dep: etaDFn(1:3, 1:ndim) = etaDFn(1:3, 1:ndim) / (elem%diam**2 * elem%CTn)

    !write(*,'(a6,i5,12es12.4)') '#DFn',elem%i, &
    !     1./(elem%diam * elem%CTn**0.5), mT,1./(elem%diam * elem%CTn**0.5)-  mT

    deallocate(psi, temp, R_s, f_s, fluxes, fluxesF, RTNflux_loc )

    !write(*,'(a6,i5,12es14.6)') 'etaDF:',elem%i, etaDFn(1:ndim)

  end subroutine FluxElemTotEstimCDSS


  ! evaluate flux convective diffusive estimator
  !> evaluation of discretization flux (diffusive) estimator
  !> \f$ \|\nabla u_h^i +  {\bf d}_{h}^i\|_{T}\f$
  subroutine FluxElemDiscEstimCDSS(Set_f_s, Set_R_s, elem, Fdeg, Fdof, fluxi, &
       !!CDfluxes,
       etaDFn)
    interface
      subroutine Set_f_s(ndimL, nbDim, Qdof, w, f_s, x )
         integer, intent(in) :: Qdof, ndimL, nbDim
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof,1:nbDim,1:ndimL), intent(inout) :: f_s
         real, dimension(1:Qdof,1 :nbDim), intent(in) :: x
      end subroutine Set_f_s
      subroutine Set_R_s(ndimL, nbDim, Qdof, w, Dw, Re_1, R_s)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
      end subroutine Set_R_s
    end interface

    type(element), intent(in) :: elem
    integer, intent(in) :: Fdeg, Fdof
    real, dimension(1:ndim, 1:Fdof), intent(in) :: fluxi
    !!real, dimension(1:elem%Qdof, 1:nbDim, 1:ndim), intent(in) :: CDfluxes
    real, dimension(:, :, :), allocatable :: fluxes, fluxesF
    !Time dep: real, dimension(1:3, 1:ndim), intent(inout) :: etaDFn ! 1 = total, 2 = space, 3 = time
    real, dimension(1:ndim), intent(inout) :: etaDFn !SS
    type(basis_rtn_fe), pointer :: loc_RTN
    type(Gauss_rule), pointer :: G_rule
    real, dimension(:,:,:), allocatable :: psi  ! local RTN basis functions
    real, dimension(:,:), allocatable :: wt     ! time derivative of the approximate sol
    real, dimension(:,:), allocatable :: RTNflux_loc  ! flux reconstruction pw linear in time
    real, dimension(:,:), allocatable :: x, Fx, f
    real, dimension(:,:,:), allocatable :: f_s, R_s, temp
    real, dimension(1:ndim) :: val
    real :: t, mT
    integer :: it, k, l, dof, Qdof, Qnum, Gnum, Gdof

    dof = elem%dof
    Qdof = elem%Qdof
    Qnum = elem%Qnum

    loc_RTN => state%loc_RTN(Fdeg)

    ! divergence of the RTN basis functions on elem in integ nodes
    allocate(psi(1:Fdof, 1:nbDim, 1:Qdof) )
    call Eval_Loc_RTN_Elem(elem, loc_RTN, psi )

    allocate( temp(1:3, 1:ndim, 1:Qdof) )  ! first index: 1 = total, 2 = space, 3 = time
    allocate(fluxes(1:Qdof, 1:nbDim, 1:ndim), fluxesF(1:Qdof, 1:nbDim, 1:ndim) )
    allocate(R_s(1:Qdof, 1:nbDim, 1:ndim) )
    allocate(f_s(1:Qdof, 1:nbDim, 1:ndim) )

    allocate(RTNflux_loc(1:ndim, 1:Fdof) )

    ! fluxes at time level t_k
    state%time%ctime = state%time%ttime
    call Eval_R_s_Elem(Set_R_s, elem, R_s(1:Qdof, 1:nbDim, 1:ndim) )

    call Eval_f_s_Elem(Set_f_s, elem, f_s(1:Qdof, 1:nbDim, 1:ndim) )

    fluxesF(1:Qdof, 1:nbDim, 1:ndim) = f_s(1:Qdof, 1:nbDim, 1:ndim) &
            - R_s(1:Qdof, 1:nbDim, 1:ndim)

    ! setting of integ  rule
    Gnum = 1   ! code for time dependent problem adapted for stationary one
    !Time dep: Gnum = 2   ! 15 is the maximal one !!!!!!!!!
    G_rule => state%space%G_rule(Gnum)
    Gdof = G_rule%Qdof

    !Time dep: etaDFn(1:3,1:ndim) = 0.
    etaDFn(1:ndim) = 0.  !SS

    ! integration over Gauss integ nodes in time
    do it=1, Gdof
       !Time dep: t = G_rule%lambda(it)
       t = 1. !SS   evaluation of entities at the previous iteration step
       state%time%ctime = state%time%ttime - state%time%tau(1) * t

       call Eval_R_s_Elem_at_time(Set_R_s, elem, R_s(1:Qdof, 1:nbDim, 1:ndim), 1.-t )
       !write(*,'(a6,2i5,20es14.6)') 'Rs=',elem%i,1,R_s(1:Qdof, 1, 1:ndim)
       !write(*,'(a6,2i5,20es14.6)') 'Rs=',elem%i,2,R_s(1:Qdof, 2, 1:ndim)

       call Eval_f_s_Elem_at_time(Set_f_s, elem, f_s(1:Qdof, 1:nbDim, 1:ndim), 1.-t )
       !write(*,'(a6,2i5,20es14.6)') 'fs=',elem%i,1,f_s(1:Qdof, 1, 1:ndim)
       !write(*,'(a6,2i5,20es14.6)') 'fs=',elem%i,2,f_s(1:Qdof, 2, 1:ndim)


       ! used in flux estimator
       fluxes(1:Qdof, 1:nbDim, 1:ndim) = f_s(1:Qdof, 1:nbDim, 1:ndim) &
            - R_s(1:Qdof, 1:nbDim, 1:ndim)

       ! linear interpolation
       !! t = 0.  ! PIECEWISE CONSTANT RECONSTRUCTION
       RTNflux_loc(1:ndim, 1:Fdof) = t * elem%RTNflux(1, 1:ndim, 1:Fdof) &
            + ( 1. - t)  *  elem%RTNflux(0, 1:ndim, 1:Fdof)

       do k=1, ndim
          temp(1:3, k, 1:Qdof) = 0.
          do l=1,nbDim
             ! total estimate
             temp(1, k, 1:Qdof) = temp(1, k, 1:Qdof) + &
                  (matmul(RTNflux_loc(k, 1:Fdof), psi(1:Fdof, l, 1:Qdof) ) &
                  - fluxes(1:Qdof, l, k) )**2

             ! space estimate (can be outside of the cycles)
             temp(2, k, 1:Qdof) = temp(1, k, 1:Qdof) + &
                  (matmul(RTNflux_loc(k, 1:Fdof), psi(1:Fdof, l, 1:Qdof) ) &
                  - fluxesF(1:Qdof, l, k) )**2

             ! time estimate
             temp(3, k, 1:Qdof) = temp(1, k, 1:Qdof) + &
                  (fluxesF(1:Qdof, l, k) - fluxes(1:Qdof, l, k) )**2

             !write(*,'(a8,2i5,20es12.4)') 'theta',elem%i,l, &
             !     matmul(fluxi(k, 1:Fdof), psi(1:Fdof, l, 1:Qdof) )
             !
             !write(*,'(a8,2i5,20es12.4)') 'flux   ',elem%i,l, fluxes(1:Qdof, l, k)
             !write(*,'(a8,2i5,20es12.4)') 'fluxCD ',elem%i,l, CDfluxes(1:Qdof, l, k)


          enddo
          !call PlotElemFunctionQ(840+state%time%iter, elem, 'V', Qnum, Qdof,  temp(k,1:Qdof) )

       enddo


       !Time dep: do l=1,3
          !Time dep: call IntegrateVectorFunction(elem, temp(l, 1:ndim, 1:Qdof), val(1:ndim))
          call IntegrateVectorFunction(elem, temp(1, 1:ndim, 1:Qdof), val(1:ndim))
          etaDFn(1:ndim) = val(1:ndim)**0.5 !SS
          !Time dep: etaDFn(l, 1:ndim) = etaDFn(l, 1:ndim) + G_rule%weights(it) * val(1:ndim)
       !Time dep: enddo

    enddo

    !Time dep: etaDFn(1:3, 1:ndim) = etaDFn(1:3, 1:ndim) * state%time%tau(1)    ! length of the time interval

    ! etaDFn stores \eta_DF^2 !!
    !Time dep: etaDFn(1:3, 1:ndim) = etaDFn(1:3, 1:ndim) / (elem%diam**2 * elem%CTn)

    !write(*,'(a6,i5,12es12.4)') '#DFn',elem%i, &
    !     1./(elem%diam * elem%CTn**0.5), mT,1./(elem%diam * elem%CTn**0.5)-  mT

    deallocate(psi, temp, R_s, f_s, fluxes, fluxesF, RTNflux_loc )

    !write(*,'(a6,i5,12es14.6)') 'etaDF:',elem%i, etaDFn(1:ndim)

  end subroutine FluxElemDiscEstimCDSS


  ! evaluate flux convective diffusive estimator
  !> evaluation of algebraic flux (diffusive) estimator
  !> \f$ \|{\bf a}_{h}^i\|_{T}\f$
  subroutine FluxElemAlgEstimCDSS(elem, Fdeg, Fdof, etaDFn)
    type(element), intent(in) :: elem
    integer, intent(in) :: Fdeg, Fdof
    !real, dimension(1:ndim, 1:Fdof), intent(in) :: fluxi
    !!real, dimension(1:elem%Qdof, 1:nbDim, 1:ndim), intent(in) :: CDfluxes
    real, dimension(:, :, :), allocatable :: fluxes, fluxesF
    !Time dep: real, dimension(1:3, 1:ndim), intent(inout) :: etaDFn ! 1 = total, 2 = space, 3 = time
    real, dimension(1:ndim), intent(inout) :: etaDFn !SS
    type(basis_rtn_fe), pointer :: loc_RTN
    type(Gauss_rule), pointer :: G_rule
    real, dimension(:,:,:), allocatable :: psi  ! local RTN basis functions
    real, dimension(:,:), allocatable :: wt     ! time derivative of the approximate sol
    real, dimension(:,:), allocatable :: RTNflux_loc  ! flux reconstruction pw linear in time
    real, dimension(:,:), allocatable :: x, Fx, f
    real, dimension(:,:,:), allocatable :: f_s, R_s, temp
    real, dimension(1:ndim) :: val
    real :: t, mT
    integer :: it, k, l, dof, Qdof, Qnum, Gnum, Gdof

    dof = elem%dof
    Qdof = elem%Qdof
    Qnum = elem%Qnum

    loc_RTN => state%loc_RTN(Fdeg)

    ! divergence of the RTN basis functions on elem in integ nodes
    allocate(psi(1:Fdof, 1:nbDim, 1:Qdof) )
    call Eval_Loc_RTN_Elem(elem, loc_RTN, psi )

    allocate( temp(1:3, 1:ndim, 1:Qdof) )  ! first index: 1 = total, 2 = space, 3 = time
    allocate(fluxes(1:Qdof, 1:nbDim, 1:ndim), fluxesF(1:Qdof, 1:nbDim, 1:ndim) )
    allocate(R_s(1:Qdof, 1:nbDim, 1:ndim) )
    allocate(f_s(1:Qdof, 1:nbDim, 1:ndim) )

    allocate(RTNflux_loc(1:ndim, 1:Fdof) )

    ! fluxes at time level t_k
    state%time%ctime = state%time%ttime
    !call Eval_R_s_Elem(Set_R_s, elem, R_s(1:Qdof, 1:nbDim, 1:ndim) )

    !call Eval_f_s_Elem(Set_f_s, elem, f_s(1:Qdof, 1:nbDim, 1:ndim) )

    fluxesF(1:Qdof, 1:nbDim, 1:ndim) = f_s(1:Qdof, 1:nbDim, 1:ndim) &
            - R_s(1:Qdof, 1:nbDim, 1:ndim)


    ! setting of integ  rule
    Gnum = 1   ! code for time dependent problem adapted for stationary one
    !Time dep: Gnum = 2   ! 15 is the maximal one !!!!!!!!!
    G_rule => state%space%G_rule(Gnum)
    Gdof = G_rule%Qdof

    !Time dep: etaDFn(1:3,1:ndim) = 0.
    etaDFn(1:ndim) = 0.  !SS

    ! integration over Gauss integ nodes in time
    do it=1, Gdof
       !Time dep: t = G_rule%lambda(it)
       t = 1. !SS   evaluation of entities at the previous iteration step
       state%time%ctime = state%time%ttime - state%time%tau(1) * t

       !call Eval_R_s_Elem_at_time(Set_R_s, elem, R_s(1:Qdof, 1:nbDim, 1:ndim), 1.-t )
       !write(*,'(a6,2i5,20es14.6)') 'Rs=',elem%i,1,R_s(1:Qdof, 1, 1:ndim)
       !write(*,'(a6,2i5,20es14.6)') 'Rs=',elem%i,2,R_s(1:Qdof, 2, 1:ndim)

       !call Eval_f_s_Elem_at_time(Set_f_s, elem, f_s(1:Qdof, 1:nbDim, 1:ndim), 1.-t )
       !write(*,'(a6,2i5,20es14.6)') 'fs=',elem%i,1,f_s(1:Qdof, 1, 1:ndim)
       !write(*,'(a6,2i5,20es14.6)') 'fs=',elem%i,2,f_s(1:Qdof, 2, 1:ndim)


       ! used in flux estimator
       fluxes(1:Qdof, 1:nbDim, 1:ndim) = f_s(1:Qdof, 1:nbDim, 1:ndim) &
            - R_s(1:Qdof, 1:nbDim, 1:ndim)

       ! linear interpolation
       !! t = 0.  ! PIECEWISE CONSTANT RECONSTRUCTION
       RTNflux_loc(1:ndim, 1:Fdof) = t * elem%RTNflux(1, 1:ndim, 1:Fdof) &
            + ( 1. - t)  *  elem%RTNflux(0, 1:ndim, 1:Fdof)

       do k=1, ndim
          temp(1:3, k, 1:Qdof) = 0.
          do l=1,nbDim
             ! total estimate
             temp(1, k, 1:Qdof) = temp(1, k, 1:Qdof) + &
                  (matmul(elem%RTNflux(0, k, 1:Fdof), psi(1:Fdof, l, 1:Qdof) ) &
                  -  matmul(elem%RTNflux(1, k, 1:Fdof), psi(1:Fdof, l, 1:Qdof) ))**2

             !write(*,*) 'elem%i =', elem%i
             !write(*,'(a20, 16es14.6)') 'd_h^nu - d_h 1 K', matmul(elem%RTNflux(0, k, 1:Fdof), psi(1:Fdof, 1, 1:Qdof) ) &
             !     -  matmul(elem%RTNflux(1, k, 1:Fdof), psi(1:Fdof, 1, 1:Qdof) )
             !write(*,'(a20, 16es14.6)') 'd_h^nu - d_h 2 K', matmul(elem%RTNflux(0, k, 1:Fdof), psi(1:Fdof, 2, 1:Qdof) ) &
             !     -  matmul(elem%RTNflux(1, k, 1:Fdof), psi(1:Fdof, 2, 1:Qdof) )
             !pause


             ! space estimate (can be outside of the cycles)
             temp(2, k, 1:Qdof) = temp(1, k, 1:Qdof) + &
                  (matmul(RTNflux_loc(k, 1:Fdof), psi(1:Fdof, l, 1:Qdof) ) &
                  - fluxesF(1:Qdof, l, k) )**2

             ! time estimate
             temp(3, k, 1:Qdof) = temp(1, k, 1:Qdof) + &
                  (fluxesF(1:Qdof, l, k) - fluxes(1:Qdof, l, k) )**2

             !write(*,'(a8,2i5,20es12.4)') 'theta',elem%i,l, &
             !     matmul(fluxi(k, 1:Fdof), psi(1:Fdof, l, 1:Qdof) )
             !
             !write(*,'(a8,2i5,20es12.4)') 'flux   ',elem%i,l, fluxes(1:Qdof, l, k)
             !write(*,'(a8,2i5,20es12.4)') 'fluxCD ',elem%i,l, CDfluxes(1:Qdof, l, k)


          enddo
          !call PlotElemFunctionQ(840+state%time%iter, elem, 'V', Qnum, Qdof,  temp(k,1:Qdof) )

          !write(*,*) 'elem%i =', elem%i
          !write(*,'(a20, 16es14.6)') 'd_h^nu - d_h K **2', temp(1, k, 1:Qdof)
          !pause

       enddo

       !Time dep: do l=1,3
          !Time dep: call IntegrateVectorFunction(elem, temp(l, 1:ndim, 1:Qdof), val(1:ndim))
          call IntegrateVectorFunction(elem, temp(1, 1:ndim, 1:Qdof), val(1:ndim))
          etaDFn(1:ndim) = val(1:ndim)**0.5 !SS
          !Time dep: etaDFn(l, 1:ndim) = etaDFn(l, 1:ndim) + G_rule%weights(it) * val(1:ndim)
       !Time dep: enddo

    enddo

    !Time dep: etaDFn(1:3, 1:ndim) = etaDFn(1:3, 1:ndim) * state%time%tau(1)    ! length of the time interval

    ! etaDFn stores \eta_DF^2 !!
    !Time dep: etaDFn(1:3, 1:ndim) = etaDFn(1:3, 1:ndim) / (elem%diam**2 * elem%CTn)

    !write(*,'(a6,i5,12es12.4)') '#DFn',elem%i, &
    !     1./(elem%diam * elem%CTn**0.5), mT,1./(elem%diam * elem%CTn**0.5)-  mT

    deallocate(psi, temp, R_s, f_s, fluxes, fluxesF, RTNflux_loc )

    !write(*,'(a6,i5,12es14.6)') 'etaDF:',elem%i, etaDFn(1:ndim)

  end subroutine FluxElemAlgEstimCDSS


  !> whether or not we will compute flux nonconformity estimator \f$ \eta_Tface \f$,
  !> Is Tface a "face of maximal length" in its direction?
  subroutine ElemFaceFNC(elem, Tface, faceFNC, whichcase)
    type(element), intent(in) :: elem
    class(element), pointer :: elem1
    integer, intent(in) :: Tface
    logical, intent(inout) :: faceFNC
    integer, intent(out) :: whichcase ! 0= no case, 1= Tface is a side of both elem and elem1, 2= Tface is a subedge of the side of elem

    integer :: ii, ie1

    !real, dimension(1:elem%face(fGdof, Tface))

    faceFNC = .false.  ! init: FNC estimator will not be computed for Tface
    whichcase = 0      ! init: no case, FNC estimator will not be computed for Tface

    if (.not. elem%HGnode) then ! Tface without HG nodes

       ii = elem%face(neigh, Tface)
       if( ii > 0) then  !! inner face
          elem1 => grid%elem(ii)
          ie1 = elem%face(nei_i,Tface)

          if (.not. elem1%HGnode) then
             faceFNC = .true.            ! Tface is a side of both elem and elem1 (i.e. Tface is not a subedge.)
             whichcase = 1
          endif ! not elem1%HGnode

          if (elem1%HGnode) then
              if (elem1%HGface(2,ie1) == 1) then ! Tface without HG nodes
                 faceFNC = .true.            ! Tface is a side of both elem and elem1 (i.e. Tface is not a subedge.)
                 whichcase = 1
              endif
          endif  ! elem1%HGnode

       endif ! ii


    endif  ! not elem%HGnode


    if (elem%HGnode) then
        if (elem%HGface(2,Tface) == 1) then ! Tface without HG nodes

           ii = elem%face(neigh, Tface)
           if( ii > 0) then  !! inner face
              elem1 => grid%elem(ii)
              ie1 = elem%face(nei_i,Tface)

              if (.not. elem1%HGnode) then
                 faceFNC = .true.            ! Tface is a side of both elem and elem1 (i.e. Tface is not a subedge.)
                 whichcase = 1
              endif ! not elem1%HGnode

              if (elem1%HGnode) then
                  if (elem1%HGface(2,ie1) == 1) then ! Tface without HG nodes
                     faceFNC = .true.            ! Tface is a side of both elem and elem1 (i.e. Tface is not a subedge.)
                     whichcase = 1
                  endif
              endif  ! elem1%HGnode

           endif ! ii

        endif

    endif  ! elem%HGnode

    if (elem%HGnode) then
       if (elem%HGface(2,Tface) > 1) then ! Tface possess at least one HG node
          faceFNC = .true.            ! Tface is a subedge of the side ie of elem
          whichcase = 2
       endif
    endif ! elem%HGnode


  end subroutine ElemFaceFNC


  subroutine BoundTraceIneqnew(elem, H_Tface, C_Gamma_K, neigh, facesharing)
    type(element), intent(in) :: elem
    real, intent(in) :: H_Tface  ! diameter of the side of the simplex elem
    logical, intent(in) :: neigh  ! trace constant for elem or its neigbor
    integer, intent(in) :: facesharing
    real, intent(out) :: C_Gamma_K
    real :: rohat ! the diameter of the inscribed ball of reference simplex
    real :: C_sd, alpha, C_HT

    if (nbDim == 2) then

      alpha = 0.730276      !! see p. 1486 in "A POSTERIORI ERROR ESTIMATIONS OF SOME CELL-CENTERED FINITE VOLUME METHODS"

      rohat = 2./(2 + (2**0.5))

      !write(*, '(a8, es12.5)') 'rohat = ', rohat

    else
      alpha = (2 * (11 + 4*(6**0.5) ) )**0.25

      !write(*, '(a8, es12.5)') 'alpha = ', alpha

      alpha = alpha / ((Pi * tanh(Pi) )**0.5 )

      !write(*, '(a8, es12.5)') 'alpha = ', alpha

      rohat = 2. * 1. * (6**0.5)/12.

      !write(*, '(a8, es12.5)') 'rohat = ', rohat

    endif


    C_sd = (alpha * (1./rohat) * (1./nbDim**0.5) )**2

    !write(*, '(a8, es12.5)') 'C_sd = ', C_sd

    if (.not. neigh) then
       C_Gamma_K = H_Tface * (elem%diam**2)
       C_Gamma_K = C_Gamma_K / (elem%area *  H_Tface )
       C_Gamma_K = (C_Gamma_K * C_sd)**0.5
                     else    ! neighbor of elem
       if (facesharing == 2) then
          C_Gamma_K = H_Tface * (elem%diam**2)
          C_Gamma_K = C_Gamma_K / (elem%area *  H_Tface )
          C_Gamma_K = (C_Gamma_K * C_sd)**0.5
       elseif (facesharing > 2) then
          C_Gamma_K = H_Tface * ( (elem%diam * (2**elem%RGlevel) )**2)
          C_Gamma_K = C_Gamma_K / ( (elem%area * (4**elem%RGlevel) ) *  H_Tface )
          C_Gamma_K = (C_Gamma_K * C_sd)**0.5
       else
          Print*, '## TROUBLE with facesharing in sub. BoundTraceIneqnew in alg_estim.f90'
       endif

    endif ! not neigh




    !write(*, '(i2, a14, es12.5, a14, es12.5)') elem%i, 'elem%diam = ', elem%diam, 'elem%area = ',elem%area

  end subroutine BoundTraceIneqnew


  subroutine BoundTraceIneqnew2(elem, H_Tface, C_Gamma_K)
    type(element), intent(in) :: elem
    real, intent(in) :: H_Tface  ! diameter of the side of the simplex elem
    real, intent(out) :: C_Gamma_K
    real :: rohat ! the diameter of the inscribed ball of reference simplex
    real :: C_sd, alpha

    if (nbDim == 2) then

      alpha = 0.730276      !! see p. 1486 in "A POSTERIORI ERROR ESTIMATIONS OF SOME CELL-CENTERED FINITE VOLUME METHODS"

      rohat = 2./(2 + (2**0.5))

      !write(*, '(a8, es12.5)') 'rohat = ', rohat

    else
      alpha = (2 * (11 + 4*(6**0.5) ) )**0.25

      !write(*, '(a8, es12.5)') 'alpha = ', alpha

      alpha = alpha / ((Pi * tanh(Pi) )**0.5 )

      !write(*, '(a8, es12.5)') 'alpha = ', alpha

      rohat = 2. * 1. * (6**0.5)/12.

      !write(*, '(a8, es12.5)') 'rohat = ', rohat

    endif


    C_sd = (alpha * (1./rohat) * (1./nbDim**0.5) )**2

    !write(*, '(a8, es12.5)') 'C_sd = ', C_sd

    C_Gamma_K = H_Tface * (elem%diam**2)
    C_Gamma_K = C_Gamma_K / (elem%area *  H_Tface )
    C_Gamma_K = (C_Gamma_K * C_sd)**0.5

    !write(*, '(i2, a14, es12.5, a14, es12.5)') elem%i, 'elem%diam = ', elem%diam, 'elem%area = ',elem%area

  end subroutine BoundTraceIneqnew2


  subroutine ComputeFaceFNCestim(elem, Tface1, Tface2, ie, Fdeg, Fdof, disc, facesharing)
    type(element), intent(inout) :: elem
    integer, intent(in) :: Tface1, Tface2, ie, Fdeg, Fdof, facesharing
    logical, intent(in) :: disc   ! discretization or algebraic flux nonconformity

    real, dimension(1:ndim) :: val, etaFNCie
    integer, dimension(:), allocatable :: smallElemsidx   ! indices of "small" elements sharing ie with elem (for whichcase = 2)
                                                          ! or index of one element sharing ie with elem (for whichcase = 1)
    type(basis_rtn_fe), pointer :: loc_RTN
    real, dimension(:,:,:), allocatable :: temp, psi, psi1, Rpsi1
    real, dimension(:,:), allocatable :: RTNfluxalg, RTNfluxalg1
    class(element), pointer ::   elem1  ! elem1 = neigh element
    type(Gauss_rule), pointer :: G_rule
    integer :: Gdof, k, l, ii, Gdof1, Mdof, ie1, i
    integer :: HGidx, Tface, idx
    real :: omega_edge
    real :: C_Gamma_K, C_Gamma_K1 ! constants from the trace inequality for elem K and K1
    real :: H_Tface   ! diameter of the side of elem containing Tface

    etaFNCie(:) = 0.

    H_Tface = 0.

    omega_edge = 0.5

    idx = 0  ! for indexing smallElemsidx

    allocate(smallElemsidx(1: facesharing - 1) )

    loc_RTN => state%loc_RTN(Fdeg)

    allocate( psi(1:Fdof, 1:nbDim, 1: maxval(elem%face(fGdof,:) ))  )

    allocate( temp(1:ndim, 1:nbDim, 1:maxval(elem%face(fGdof,:) )) )


    psi(:, :, :) = 0.

    do Tface = Tface1, Tface2

       if (facesharing == 2) then ! whichcase = 1, no HG node is present on ie
          HGidx = 1
          ! local basis RTN functions at edge integration nodes, no HG node is present on ie
          call Eval_Loc_RTN_Edge(elem, loc_RTN, ie, psi(1:Fdof, 1:nbDim, 1:elem%face(fGdof, ie) ) )
       elseif (facesharing > 2) then ! whichcase = 2, at least one HG node is present on ie
          HGidx = elem%HGface(2,Tface)
          call Eval_Loc_RTN_Edge_HG(elem, loc_RTN, Tface, HGidx, psi(1:Fdof, 1:nbDim, 1:elem%face(fGdof, Tface) ), .false. )
       else
          print*, '## TROUBLE with facesharing in sub. ComputeFaceFNCestim in alg_estim.f90'
       endif

       G_rule => state%space%G_rule(elem%face(fGnum, Tface)) !choice of edge quadrature
       Gdof = state%space%G_rule(elem%face(fGnum, Tface))%Qdof !number of edge integration nodes

       ii = elem%face(neigh, Tface)
       if( ii < 1) then
         Print*, '## Trouble in sub. ComputeFaceFNCestim. Tface cannot be a boundary face!'
         stop
       endif
       elem1 => grid%elem(ii)
       ie1 = elem%face(nei_i,Tface)
       Gdof1 = elem1%face(fGdof,ie1)

       if(Gdof /= Gdof1) print*,'## Trouble in ComputeFaceFNCestim with number of edge integ nodes.'

       allocate( psi1(1:Fdof, 1:nbDim, 1:Gdof1) )
       psi1(:, :, :) = 0.

       if (elem1%HGnode) then
          if (elem1%HGface(2, ie1) > 1) then
             Print*, '## Trouble in sub. ComputeFaceFNCestim. ie1 cannot possess any HG nodes!'
             stop
          endif
       endif

       allocate( Rpsi1(1:Fdof, 1:nbDim, 1: Gdof1) )
       Rpsi1(:, :, :) = 0.
       call Eval_Loc_RTN_Edge(elem1, loc_RTN, ie1, Rpsi1(1:Fdof, 1:nbDim, 1:Gdof1 ) )  !no possibility true/false in Eval_Loc_RTN_Edge

       ! reordering
       do l=1, Gdof1
          psi1(1:Fdof, 1:nbDim, l) = Rpsi1(1:Fdof, 1:nbDim, Gdof1 -l + 1)
       enddo
       deallocate(Rpsi1)

       Mdof = Gdof

       do k=1, ndim
          temp(k, 1:nbDim, 1:maxval(elem%face(fGdof,:) )) = 0.  !temp allocated for 1:maxval(elem%face(fGdof,:) )
          do l=1,nbDim
             if (disc) then  ! discretization flux nonconformity
                temp(k, l, 1:Mdof) = matmul(elem%RTNflux(1, k, 1:Fdof), psi(1:Fdof, l, 1:Mdof) ) &
                                     -  matmul(elem1%RTNflux(1, k, 1:Fdof), psi1(1:Fdof, l, 1:Mdof) ) !elem1 RTNflux at previous lvl saved in elem1%RTNflux(1,..)

                !print*,''
                !write(*, '(a27, 15es14.6)') 'elem1%RTNflux(1, 1, 1:Fdof)', elem1%RTNflux(1, 1, 1:Fdof)
                !print*,''

             else ! algebraic flux nonconformity
                allocate( RTNfluxalg(1:ndim, 1:Fdof), RTNfluxalg1(1:ndim, 1:Fdof) )
                RTNfluxalg(k, 1:Fdof) = elem%RTNflux(0, k, 1:Fdof) - elem%RTNflux(1, k, 1:Fdof)
                RTNfluxalg1(k, 1:Fdof) = elem1%RTNflux(0, k, 1:Fdof) - elem1%RTNflux(1, k, 1:Fdof)

                if (elem%i == 1) then
                   !write(*, '(a24, 15es14.6)') 'RTNfluxalg(k, 1:Fdof) =', RTNfluxalg(k, 1:Fdof)
                   !print*,''
                   !write(*, '(a24, 15es14.6)')  'RTNfluxalg1(k, 1:Fdof) =', RTNfluxalg1(k, 1:Fdof)
                   !print*,''
                endif

                temp(k, l, 1:Mdof) = matmul(RTNfluxalg(k, 1:Fdof), psi(1:Fdof, l, 1:Mdof) ) &
                                     -  matmul(RTNfluxalg1(k, 1:Fdof), psi1(1:Fdof, l, 1:Mdof) )
                deallocate(RTNfluxalg, RTNfluxalg1)
             endif

          enddo ! l
       enddo ! k

       deallocate( psi1 )

       idx = idx + 1

       smallElemsidx(idx) = ii  ! saving of index of an element sharing ie with elem

       do k=1, ndim
          call IntegrateSquareFunctionNormalEdge(elem, Tface, temp(k, 1:nbDim, 1:elem%face(fGdof,Tface)), val(k))
       enddo ! k

       etaFNCie(1:ndim) = etaFNCie(1:ndim) + val(1:ndim)**0.5  ! sum over subedges of the jump of the normal component of d_h^i/a_h^i

       H_Tface = H_Tface + elem%dn(Tface)  ! sum over diameters of subedges of side of elem

    enddo ! Tface

    etaFNCie(1:ndim) = etaFNCie(1:ndim)**2 ! norm of the jump of the normal component of d_h^i/a_h^i on ie is squared

    call BoundTraceIneqnew(elem, H_Tface, C_Gamma_K, .false., facesharing)

    elem1 => grid%elem(smallElemsidx(idx))

    call BoundTraceIneqnew(elem1, H_Tface, C_Gamma_K1, .true., facesharing)

    ! How can I compute C_Gamma_K for a coarse element sharing "the whole" ie with elem if it is not in the actual grid any more??

    etaFNCie(1:ndim) = etaFNCie(1:ndim) * H_Tface * (omega_edge**2) * ((C_Gamma_K**2) + (C_Gamma_K1**2) )   ! eta_\Gamma**2 from the article (for discret. or algebraic component)

    if (idx /= facesharing - 1) then
       Print*, '## Trouble in sub. ComputeFaceFNCestim with the number of elements sharing a side with elem.'
       stop
    endif

    ! distribute the contribution of eta_\Gamma**2 among elements sharing ie with elem
    do i=1, idx  ! loop through indices of elements sharing ie with elem (except for elem)
       elem1 => grid%elem(smallElemsidx(i))
       if (disc) then  ! discretization flux nonconformity
          elem1%estimFNCD(1:ndim) = elem1%estimFNCD(1:ndim) + (etaFNCie(1:ndim) / facesharing)
                 else  ! algebraic flux nonconformity
          elem1%estimFNCA(1:ndim) = elem1%estimFNCA(1:ndim) + (etaFNCie(1:ndim) / facesharing)
       endif
    enddo !i

    if (disc) then  ! discretization flux nonconformity
       elem%estimFNCD(1:ndim) = elem%estimFNCD(1:ndim) + (etaFNCie(1:ndim) / facesharing)
              else  ! algebraic flux nonconformity
       elem%estimFNCA(1:ndim) = elem%estimFNCA(1:ndim) + (etaFNCie(1:ndim) / facesharing)
    endif



    deallocate(smallElemsidx)
    deallocate(psi)
    deallocate(temp)

  end subroutine ComputeFaceFNCestim


  subroutine ComputeFaceFNCestimnew(elem, Tface1, Tface2, ie, Fdeg, Fdof, sort, facesharing)
    type(element), intent(inout) :: elem
    integer, intent(in) :: Tface1, Tface2, ie, Fdeg, Fdof, facesharing
    character(len=3) :: sort   ! discretization, algebraic, or total flux nonconformity

    real, dimension(1:ndim) :: val, etaFNCie
    type(basis_rtn_fe), pointer :: loc_RTN
    real, dimension(:,:,:), allocatable :: temp, psi, psi1, Rpsi1
    real, dimension(:,:), allocatable :: RTNfluxalg, RTNfluxalg1
    class(element), pointer ::   elem1  ! elem1 = neigh element
    type(Gauss_rule), pointer :: G_rule
    integer :: Gdof, k, l, ii, Gdof1, Mdof, ie1, i
    integer :: HGidx, Tface
    real :: omega_edge
    real :: C_Gamma_K ! constants from the trace inequality for elem K
    real :: H_Tface   ! diameter of the side of elem containing Tface

    etaFNCie(:) = 0.

    H_Tface = 0.

    omega_edge = 0.5

    loc_RTN => state%loc_RTN(Fdeg)

    allocate( psi(1:Fdof, 1:nbDim, 1: maxval(elem%face(fGdof,:) ))  )

    allocate( temp(1:ndim, 1:nbDim, 1:maxval(elem%face(fGdof,:) )) )


    psi(:, :, :) = 0.

    do Tface = Tface1, Tface2

       if (facesharing == 2) then ! whichcase = 1, no HG node is present on ie
          HGidx = 1
          ! local basis RTN functions at edge integration nodes, no HG node is present on ie
          call Eval_Loc_RTN_Edge(elem, loc_RTN, ie, psi(1:Fdof, 1:nbDim, 1:elem%face(fGdof, ie) ) )
       elseif (facesharing > 2) then ! whichcase = 2, at least one HG node is present on ie
          HGidx = elem%HGface(2,Tface)
          call Eval_Loc_RTN_Edge_HG(elem, loc_RTN, Tface, HGidx, psi(1:Fdof, 1:nbDim, 1:elem%face(fGdof, Tface) ), .false. )
       else
          print*, '## TROUBLE with facesharing in sub. ComputeFaceFNCestimnew in alg_estim.f90'
       endif

       G_rule => state%space%G_rule(elem%face(fGnum, Tface)) !choice of edge quadrature
       Gdof = state%space%G_rule(elem%face(fGnum, Tface))%Qdof !number of edge integration nodes

       ii = elem%face(neigh, Tface)
       if( ii < 1) then
         Print*, '## Trouble in sub. ComputeFaceFNCestimnew. Tface cannot be a boundary face!'
         stop
       endif
       elem1 => grid%elem(ii)
       ie1 = elem%face(nei_i,Tface)
       Gdof1 = elem1%face(fGdof,ie1)

       if(Gdof /= Gdof1) print*,'## Trouble in ComputeFaceFNCestimnew with number of edge integ nodes.'

       allocate( psi1(1:Fdof, 1:nbDim, 1:Gdof1) )
       psi1(:, :, :) = 0.

       if (elem1%HGnode) then
          if (elem1%HGface(2, ie1) > 1) then
             Print*, '## Trouble in sub. ComputeFaceFNCestimnew. ie1 cannot possess any HG nodes!'
             stop
          endif
       endif

       allocate( Rpsi1(1:Fdof, 1:nbDim, 1: Gdof1) )
       Rpsi1(:, :, :) = 0.
       call Eval_Loc_RTN_Edge(elem1, loc_RTN, ie1, Rpsi1(1:Fdof, 1:nbDim, 1:Gdof1 ) )  !no possibility true/false in Eval_Loc_RTN_Edge

       ! reordering
       do l=1, Gdof1
          psi1(1:Fdof, 1:nbDim, l) = Rpsi1(1:Fdof, 1:nbDim, Gdof1 -l + 1)
       enddo
       deallocate(Rpsi1)

       Mdof = Gdof

       do k=1, ndim
          temp(k, 1:nbDim, 1:maxval(elem%face(fGdof,:) )) = 0.  !temp allocated for 1:maxval(elem%face(fGdof,:) )
          do l=1,nbDim
             if (sort == 'dis') then  ! discretization flux nonconformity
                temp(k, l, 1:Mdof) = matmul(elem%RTNflux(1, k, 1:Fdof), psi(1:Fdof, l, 1:Mdof) ) &
                                     -  matmul(elem1%RTNflux(1, k, 1:Fdof), psi1(1:Fdof, l, 1:Mdof) ) !elem1 RTNflux at previous lvl saved in elem1%RTNflux(1,..)

                !write(*, '(a30, i2, a5, 6es14.6)') 'temp(1, 1:elem%face(fGdof,', Tface, ')) =', temp(1, l, 1:Mdof)
                !print*,''
                !write(*, '(a27, 15es14.6)') 'elem1%RTNflux(1, 1, 1:Fdof)', elem1%RTNflux(1, 1, 1:Fdof)
                !print*,''


             else if (sort == 'alg') then ! algebraic flux nonconformity
                allocate( RTNfluxalg(1:ndim, 1:Fdof), RTNfluxalg1(1:ndim, 1:Fdof) )
                RTNfluxalg(k, 1:Fdof) = elem%RTNflux(0, k, 1:Fdof) - elem%RTNflux(1, k, 1:Fdof)
                RTNfluxalg1(k, 1:Fdof) = elem1%RTNflux(0, k, 1:Fdof) - elem1%RTNflux(1, k, 1:Fdof)

                if (elem%i == 1) then
                   !write(*, '(a24, 15es14.6)') 'RTNfluxalg(k, 1:Fdof) =', RTNfluxalg(k, 1:Fdof)
                   !print*,''
                   !write(*, '(a24, 15es14.6)')  'RTNfluxalg1(k, 1:Fdof) =', RTNfluxalg1(k, 1:Fdof)
                   !print*,''
                endif

                temp(k, l, 1:Mdof) = matmul(RTNfluxalg(k, 1:Fdof), psi(1:Fdof, l, 1:Mdof) ) &
                                     -  matmul(RTNfluxalg1(k, 1:Fdof), psi1(1:Fdof, l, 1:Mdof) )
                deallocate(RTNfluxalg, RTNfluxalg1)

             else if (sort == 'tot') then ! total flux nonconformity
                temp(k, l, 1:Mdof) = matmul(elem%RTNflux(0, k, 1:Fdof), psi(1:Fdof, l, 1:Mdof) ) &
                                     -  matmul(elem1%RTNflux(0, k, 1:Fdof), psi1(1:Fdof, l, 1:Mdof) ) ! t_h^i = d_h^i + d_h^{i + nu} - d_h^i = d_h^{i + nu}
             endif

          enddo ! l
       enddo ! k

       if (state%space%adapt%adapt_level == -1) then
          if (sort == 'disc') then
          if (elem%i == 6) then
             if (elem%face(neigh, Tface) == 10) then
                temp(1, 1, 1:Mdof) = matmul(elem%RTNflux(1, 1, 1:Fdof), psi(1:Fdof, 1, 1:Mdof) )
                temp(1, 2, 1:Mdof) = matmul(elem%RTNflux(1, 1, 1:Fdof), psi(1:Fdof, 2, 1:Mdof) )
                write(*,*) 'disc', disc
                write(*,*) 'elem%i', elem%i
                write(*,*) 'elem%n(Tface, 1:2)', elem%n(Tface, 1:2)
                write(*,*) '(elem%n(Tface, 1)/elem%dn(Tface))', (elem%n(Tface, 1)/elem%dn(Tface))
                write(*, '(a30, i2, a5, 6es14.6)') 'temp(1, 1:elem%face(fGdof,', Tface, ')) =', temp(1, 1, 1:Mdof)
                write(*, '(a30, i2, a5, 6es14.6)') 'temp(1, 1:elem%face(fGdof,', Tface, ')) =', temp(1, 2, 1:Mdof)
                write(*, '(a30, 6es14.6)') 'dh1  n1 R', temp(1, 1, 1:Mdof)*(elem%n(Tface, 1)/elem%dn(Tface))
                write(*, '(a30, 6es14.6)') 'dh2  n2 R',  temp(1, 2, 1:Mdof)*(elem%n(Tface, 2)/elem%dn(Tface))
                pause
          endif
             endif

          if (elem%i == 6) then
             if (elem%face(neigh, Tface) == 10) then
                temp(1, 1, 1:Mdof) = matmul(elem1%RTNflux(1, 1, 1:Fdof), psi1(1:Fdof, 1, 1:Mdof) )
                temp(1, 2, 1:Mdof) = matmul(elem1%RTNflux(1, 1, 1:Fdof), psi1(1:Fdof, 2, 1:Mdof) )
                write(*,*) 'disc', disc
                write(*,*) 'elem%i', elem%i
                write(*,*) 'elem%n(Tface, 1:2)', elem%n(Tface, 1:2)
                write(*,*) '(elem%n(Tface, 1)/elem%dn(Tface))', (elem%n(Tface, 1)/elem%dn(Tface))
                write(*, '(a30, i2, a5, 6es14.6)') 'temp(1, 1:elem%face(fGdof,', Tface, ')) =', temp(1, 1, 1:Mdof)
                write(*, '(a30, i2, a5, 6es14.6)') 'temp(1, 1:elem%face(fGdof,', Tface, ')) =', temp(1, 2, 1:Mdof)
                write(*, '(a30, 6es14.6)') 'dh1  n1 L', temp(1, 1, 1:Mdof)*(elem%n(Tface, 1)/elem%dn(Tface))
                write(*, '(a30, 6es14.6)') 'dh2  n2 L',  temp(1, 2, 1:Mdof)*(elem%n(Tface, 2)/elem%dn(Tface))
                pause
          endif
             endif


          if (elem%i == 6) then
                if (elem%face(neigh, Tface) == 11) then
                   temp(1, 1, 1:Mdof) = matmul(elem%RTNflux(1, 1, 1:Fdof), psi(1:Fdof, 1, 1:Mdof) )
                   temp(1, 2, 1:Mdof) = matmul(elem%RTNflux(1, 1, 1:Fdof), psi(1:Fdof, 2, 1:Mdof) )
                   write(*,*) 'disc', disc
                   write(*,*) 'elem%i', elem%i
                   write(*,*) 'elem%n(Tface, 1:2)', elem%n(Tface, 1:2)
                   write(*,*) '(elem%n(Tface, 1)/elem%dn(Tface))', (elem%n(Tface, 1)/elem%dn(Tface))
                   write(*, '(a30, i2, a5, 6es14.6)') 'temp(1, 1:elem%face(fGdof,', Tface, ')) =', temp(1, 1, 1:Mdof)
                   write(*, '(a30, i2, a5, 6es14.6)') 'temp(1, 1:elem%face(fGdof,', Tface, ')) =', temp(1, 2, 1:Mdof)
                   write(*, '(a30, 6es14.6)') 'dh1  n1 R', temp(1, 1, 1:Mdof)*(elem%n(Tface, 1)/elem%dn(Tface))
                   write(*, '(a30, 6es14.6)') 'dh2  n2 R',  temp(1, 2, 1:Mdof)*(elem%n(Tface, 2)/elem%dn(Tface))
                   pause
          endif
             endif

          if (elem%i == 6) then
                if (elem%face(neigh, Tface) == 11) then
                   temp(1, 1, 1:Mdof) = matmul(elem1%RTNflux(1, 1, 1:Fdof), psi1(1:Fdof, 1, 1:Mdof) )
                   temp(1, 2, 1:Mdof) = matmul(elem1%RTNflux(1, 1, 1:Fdof), psi1(1:Fdof, 2, 1:Mdof) )
                   write(*,*) 'disc', disc
                   write(*,*) 'elem%i', elem%i
                   write(*,*) 'elem%n(Tface, 1:2)', elem%n(Tface, 1:2)
                   write(*,*) '(elem%n(Tface, 1)/elem%dn(Tface))', (elem%n(Tface, 1)/elem%dn(Tface))
                   write(*, '(a30, i2, a5, 6es14.6)') 'temp(1, 1:elem%face(fGdof,', Tface, ')) =', temp(1, 1, 1:Mdof)
                   write(*, '(a30, i2, a5, 6es14.6)') 'temp(1, 1:elem%face(fGdof,', Tface, ')) =', temp(1, 2, 1:Mdof)
                   write(*, '(a30, 6es14.6)') 'dh1  n1 L', temp(1, 1, 1:Mdof)*(elem%n(Tface, 1)/elem%dn(Tface))
                   write(*, '(a30, 6es14.6)') 'dh2  n2 L',  temp(1, 2, 1:Mdof)*(elem%n(Tface, 2)/elem%dn(Tface))
                   pause
          endif
             endif

              endif
              endif

       deallocate( psi1 )


       do k=1, ndim
          call IntegrateSquareFunctionNormalEdge(elem, Tface, temp(k, 1:nbDim, 1:elem%face(fGdof,Tface)), val(k))
       enddo ! k

       etaFNCie(1:ndim) = etaFNCie(1:ndim) + val(1:ndim)**0.5  ! sum over subedges of the jump of the normal component of d_h^i/a_h^i

       H_Tface = H_Tface + elem%dn(Tface)  ! sum over diameters of subedges of side of elem

    enddo ! Tface


    call BoundTraceIneqnew2(elem, H_Tface, C_Gamma_K)

    !if (facesharing == 2) then
    !  write(*,'(a30, i5, es12.4)') 'elem%i, etaFNCie(1:ndim):', elem%i, etaFNCie(1:ndim)
    !  pause
    !endif

    etaFNCie(1:ndim) = etaFNCie(1:ndim) * ((H_Tface)**0.5) * C_Gamma_K

    if (facesharing == 2) then ! whichcase = 1, no HG node is present on ie
       etaFNCie(1:ndim) = etaFNCie(1:ndim) * omega_edge
    endif


    if (sort == 'dis') then  ! discretization flux nonconformity
       elem%estimFNCD(1:ndim) = elem%estimFNCD(1:ndim) + etaFNCie(1:ndim)
              else if (sort == 'alg') then ! algebraic flux nonconformity
       elem%estimFNCA(1:ndim) = elem%estimFNCA(1:ndim) + etaFNCie(1:ndim)
              else if (sort == 'tot') then ! total flux nonconformity
       elem%estimFNCT(1:ndim) = elem%estimFNCT(1:ndim) + etaFNCie(1:ndim)
    endif

    deallocate(psi)
    deallocate(temp)

  end subroutine ComputeFaceFNCestimnew


  subroutine FluxNCEstimSS(Fdeg, Fdof, disc)
    class(element), pointer :: elem
    integer, intent(in) :: Fdeg, Fdof
    logical, intent(in) :: disc   ! discretization or algebraic flux nonconformity

    integer :: HGidx, Tface, Tface1, Tface2, HGidx1, ie, i
    integer, dimension(:), allocatable :: HGvertex
    logical :: faceFNC
    integer :: whichcase, facesharing

    do i = 1, grid%nelem
       elem => grid%elem(i)

       if (disc) then
          elem%estimFNCD(1:ndim) = 0.      ! init
       else
          elem%estimFNCA(1:ndim) = 0.
       endif
    enddo ! i

    do i = 1, grid%nelem
       elem => grid%elem(i)

       if(elem%HGnode) then
          if(elem%type /= 3) print*,' TROUBLES in alg_estim.f90 with HG'

          allocate( HGvertex(4) )
          HGvertex(1:3)   = elem%HGvertex(1:3 )
          HGvertex(4) = elem%flen + 1

          !write(*,'(i5,a2,10i5)') elem%i,'|',elem%HGvertex(:)
          !write(*,'(i5,a2,10i5)') elem%i,'|',HGvertex(:)
          !write(*,'(i5,a2,10i5)') elem%i,'|',elem%HGface(1,:)
          !write(*,'(i5,a2,10i5)') elem%i,'|',elem%HGface(2,:)
          !print*,'-------------------'
          !write(*,'(a5,i2,a15)') 'Elem', elem%i, 'has HG nodes.'
       endif ! elem%HGnode

       do ie = 1, 3!loop through edges of triangle  (meshes with HG nodes  included)

          if(.not. elem%HGnode) then
             Tface1 = ie   ! only one element edge on the triangle edge
             Tface2 = ie
             HGidx = 1
             ! local basis RTN functions at edge integration nodes, no HG node is present on ie
             !!!call Eval_Loc_RTN_Edge(elem, loc_RTN, ie, psi(1:Fdof, 1:nbDim, 1:elem%face(fGdof, ie) ) )
          else
             Tface1 = HGvertex(ie)
             Tface2 = HGvertex(ie+1) - 1
          endif

          do Tface = Tface1, Tface2

             call ElemFaceFNC(elem, Tface, faceFNC, whichcase)
             if (faceFNC) then
                if (whichcase == 1) then
                   facesharing = 2    ! two elements share Tface
                endif
                if (whichcase == 2) then
                   facesharing = HGvertex(ie+1) - HGvertex(ie) + 1   ! Tface is a subedge of the side ie of elem, ie is shared by
                endif                                                ! HGvertex(ie+1) - 1 - [HGvertex(ie) - 1] + 1 elements
                call ComputeFaceFNCestim(elem, Tface1, Tface2, ie, Fdeg, Fdof, .true., facesharing)    ! discretization flux nonconformity estimator
                call ComputeFaceFNCestim(elem, Tface1, Tface2, ie, Fdeg, Fdof, .false., facesharing)   ! algebraic flux nonconformity estimator
                goto 100   ! do not go through the rest of subedges of ie side, FNC estimator for ie has been already computed
             endif ! faceFNC

          enddo ! Tface

100       continue

       enddo ! ie  (side of elem)

       if(elem%HGnode) deallocate(HGvertex)

    enddo !i

    ! square root of contributions to elem
    do i = 1, grid%nelem
       elem => grid%elem(i)

       if (disc) then  ! discretization flux nonconformity
          elem%estimFNCD(1:ndim) = elem%estimFNCD(1:ndim)**0.5
                 else  ! algebraic flux nonconformity
          elem%estimFNCA(1:ndim) = elem%estimFNCA(1:ndim)**0.5
       endif

    enddo ! i



  end subroutine FluxNCEstimSS


  subroutine FluxNCElemEstimSSnew(elem, Fdeg, Fdof, sort)
    type(element), intent(inout) :: elem
    integer, intent(in) :: Fdeg, Fdof
    character(len=3) :: sort   ! discretization, algebraic, or total flux nonconformity

    integer :: HGidx, Tface, Tface1, Tface2, HGidx1, ie, i
    integer, dimension(:), allocatable :: HGvertex
    logical :: faceFNC
    integer :: whichcase, facesharing


       if (sort == 'dis') then
          elem%estimFNCD(1:ndim) = 0.      ! init
       else if (sort == 'alg') then
          elem%estimFNCA(1:ndim) = 0.
       else if (sort == 'tot') then
          elem%estimFNCT(1:ndim) = 0.
       endif

       if(elem%HGnode) then
          if(elem%type /= 3) print*,' TROUBLES in alg_estim.f90 with HG'

          allocate( HGvertex(4) )
          HGvertex(1:3)   = elem%HGvertex(1:3 )
          HGvertex(4) = elem%flen + 1

          !write(*,'(i5,a2,10i5)') elem%i,'|',elem%HGvertex(:)
          !write(*,'(i5,a2,10i5)') elem%i,'|',HGvertex(:)
          !write(*,'(i5,a2,10i5)') elem%i,'|',elem%HGface(1,:)
          !write(*,'(i5,a2,10i5)') elem%i,'|',elem%HGface(2,:)
          !print*,'-------------------'
          !write(*,'(a5,i2,a15)') 'Elem', elem%i, 'has HG nodes.'
       endif ! elem%HGnode

       do ie = 1, 3!loop through edges of triangle  (meshes with HG nodes  included)

          if(.not. elem%HGnode) then
             Tface1 = ie   ! only one element edge on the triangle edge
             Tface2 = ie
             HGidx = 1
             ! local basis RTN functions at edge integration nodes, no HG node is present on ie
             !!!call Eval_Loc_RTN_Edge(elem, loc_RTN, ie, psi(1:Fdof, 1:nbDim, 1:elem%face(fGdof, ie) ) )
          else
             Tface1 = HGvertex(ie)
             Tface2 = HGvertex(ie+1) - 1
          endif

          do Tface = Tface1, Tface2

             call ElemFaceFNC(elem, Tface, faceFNC, whichcase)
             if (faceFNC) then
                if (whichcase == 1) then
                   facesharing = 2    ! two elements share Tface
                endif
                if (whichcase == 2) then
                   facesharing = HGvertex(ie+1) - HGvertex(ie) + 1   ! Tface is a subedge of the side ie of elem, ie is shared by
                endif                                                ! HGvertex(ie+1) - 1 - [HGvertex(ie) - 1] + 1 elements
                if (sort == 'dis') then
                call ComputeFaceFNCestimnew(elem, Tface1, Tface2, ie, Fdeg, Fdof, 'dis', facesharing)    ! discretization flux nonconformity estimator
                else if (sort == 'alg') then
                call ComputeFaceFNCestimnew(elem, Tface1, Tface2, ie, Fdeg, Fdof, 'alg', facesharing)   ! algebraic flux nonconformity estimator
                else if (sort == 'tot') then
                call ComputeFaceFNCestimnew(elem, Tface1, Tface2, ie, Fdeg, Fdof, 'tot', facesharing)   ! total flux nonconformity estimator
                endif
                goto 100   ! do not go through the rest of subedges of ie side, FNC estimator for ie has been already computed
             endif ! faceFNC

          enddo ! Tface

100       continue


       enddo ! ie  (side of elem)

       if(elem%HGnode) deallocate(HGvertex)


  end subroutine FluxNCElemEstimSSnew


  !> evaluation of discretization or algebraic flux nonconformity (diffusive) estimator, i.e.
  !> \f$ \eta_{FNCD,T}^i := \sum_{\gamma \in \mathcal(F)_T} C_{\gamma, T, b} w_{\gamma} h_{\gamma}^{1/2} \|[{\bf d}_{h}^i {\cdot} {\bf n}] \|_{\gamma}\f$ or
  !> \f$ \eta_{FNCA,T}^i := \sum_{\gamma \in \mathcal(F)_T} C_{\gamma, T, b} w_{\gamma} h_{\gamma}^{1/2} \|[({\bf d}_{h}^{i+\nu} - {\bf d}_{h}^i) {\cdot} {\bf n}] \|_{\gamma}\f$
  subroutine FluxNCElemDiscEstimSS(elem, Fdeg, Fdof, disc, etaFNCD)
    type(element), intent(in) :: elem
    integer, intent(in) :: Fdeg, Fdof
    logical, intent(in) :: disc   ! discretization or algebraic flux nonconformity
    type(basis_rtn_fe), pointer :: loc_RTN
    real, dimension(1:ndim), intent(inout) :: etaFNCD
    real, dimension(1:ndim) :: val, etaFNCDie
    class(element), pointer ::   elem1  ! elem1 = neigh element
    real, dimension(:,:,:), allocatable :: temp, extract, psi, psi1, Rpsi1
    real, dimension(:,:), allocatable :: RTNfluxalg, RTNfluxalg1
    real :: C_edgTb, omega_edge, C_HT, max_C_HT
    integer :: ie, ie1

    type(Gauss_rule), pointer :: G_rule
    integer :: Gdof, k, l, ii, Gdof1, Mdof
    integer :: HGidx, Tface, Tface1, Tface2, HGidx1
    integer, dimension(:), allocatable :: HGvertex


    etaFNCD(:) = 0.

    C_HT = 0.
    max_C_HT = 0.

    loc_RTN => state%loc_RTN(Fdeg)

    allocate( psi(1:Fdof, 1:nbDim, 1: maxval(elem%face(fGdof,:) ))  )

    allocate( temp(1:ndim, 1:nbDim, 1:maxval(elem%face(fGdof,:) )) )

    allocate( extract(1:ndim, 1:nbDim, 1:maxval(elem%face(fGdof,:) )) )

    if(elem%HGnode) then
       if(elem%type /= 3) print*,' TROUBLES in alg_estim.f90 with HG'

       allocate( HGvertex(4) )
       HGvertex(1:3)   = elem%HGvertex(1:3 )
       HGvertex(4) = elem%flen + 1

       !write(*,'(i5,a2,10i5)') elem%i,'|',elem%HGvertex(:)
       !write(*,'(i5,a2,10i5)') elem%i,'|',HGvertex(:)
       !write(*,'(i5,a2,10i5)') elem%i,'|',elem%HGface(1,:)
       !write(*,'(i5,a2,10i5)') elem%i,'|',elem%HGface(2,:)
       !print*,'-------------------'
       !write(*,'(a5,i2,a15)') 'Elem', elem%i, 'has HG nodes.'
    endif

    !write(*, '(a18, i2)') 'state%time%iter = ', state%time%iter

    !TEST
    !etaFNCDie(:) = 0.
    !val(:)= 0.
    !if (elem%i == 1) then
    !do ie = 1, 3
    !temp(1, 1, 1:elem%face(fGdof, ie)) = 1.
    !call Eval_Loc_RTN_Edge(elem, loc_RTN, ie, psi(1:Fdof, 1:nbDim, 1:elem%face(fGdof, ie) ) )
    !call IntegrateFunctionNormalEdge(elem, ie, psi(15, 1:nbDim, 1:elem%face(fGdof,ie)), temp(1, 1, 1:elem%face(fGdof, ie)), val(1) )
    !etaFNCDie(1) = etaFNCDie(1) + val(1)
    !enddo
    !write(*,*) 'etaFNCDie(1)', etaFNCDie(1)
    !pause
    !endif
    !EST

    !write(*,'(a10,i3)') '*****Elem = ', elem%i

    do ie = 1, 3!loop through edges of triangle  (meshes with HG nodes  included)

       etaFNCDie(:) = 0.
       psi(:, :, :) = 0.



       if(.not. elem%HGnode) then
          Tface1 = ie   ! only one element edge on the triangle edge
          Tface2 = ie
          HGidx = 1
          ! local basis RTN functions at edge integration nodes, no HG node is present on ie
          call Eval_Loc_RTN_Edge(elem, loc_RTN, ie, psi(1:Fdof, 1:nbDim, 1:elem%face(fGdof, ie) ) )
       else
          Tface1 = HGvertex(ie)
          Tface2 = HGvertex(ie+1) - 1
       endif

       do Tface = Tface1, Tface2
          if( elem%HGnode) then
             HGidx = elem%HGface(2,Tface)
             if (elem%HGface(2,Tface) == 1) then !no HG node is present on Tface
                call Eval_Loc_RTN_Edge(elem, loc_RTN, ie, psi(1:Fdof, 1:nbDim, 1:elem%face(fGdof, ie) ) )
             else !HG node is present on ie (i.e. Tface is a subedge)
                call Eval_Loc_RTN_Edge_HG(elem, loc_RTN, Tface, HGidx, psi(1:Fdof, 1:nbDim, 1:elem%face(fGdof, Tface) ), .false. )
             endif
          endif


          !write(*,'(a4,6i5)') 'HG:',elem%i, ie, Tface, Tface1, Tface2, HGidx

          G_rule => state%space%G_rule(elem%face(fGnum, Tface)) !choice of edge quadrature
          Gdof = state%space%G_rule(elem%face(fGnum, Tface))%Qdof !number of edge integration nodes


          !if (elem%HGface(1, ie) == 1) write(*, '(a20, i2)') 'elem%HGface(1, ie)', elem%HGface(1, ie)



          ii = elem%face(neigh, Tface)
          if( ii > 0) then  !! inner face
             elem1 => grid%elem(ii)
             ie1 = elem%face(nei_i,Tface)
             Gdof1 = elem1%face(fGdof,ie1)

             if(Gdof /= Gdof1) print*,'## Trouble in FluxNCElemDiscEstimSS.'


             if (HGidx > 1) then
                !if (grid%elem( elem%face(neigh, Tface1) )%face(fGdof, elem%face(nei_i,Tface1) )&
                !    + grid%elem( elem%face(neigh, Tface2) )%face(fGdof, elem%face(nei_i,Tface2) )&
                !     /= elem%face(fGdof,ie)  ) print*,'## Trouble in FluxNCElemDiscEstimSS with Gdof.'
                !write(*,'(a100,i5)') 'grid%elem( elem%face(neigh, Tface1) )%face(fGdof, elem%face(nei_i,Tface1) )', &
                !                    grid%elem( elem%face(neigh, Tface1) )%face(fGdof, elem%face(nei_i,Tface1) )
                !write(*,'(a100,i5)')  'grid%elem( elem%face(neigh, Tface2) )%face(fGdof, elem%face(nei_i,Tface2)', &
                !                    grid%elem( elem%face(neigh, Tface2) )%face(fGdof, elem%face(nei_i,Tface2) )
                !write(*,'(a10,i5)') 'elem%face(fGdof,ie)', elem%face(fGdof,ie)
                !pause

             endif

             C_HT = elem%diam / elem1%diam

             if (C_HT > max_C_HT) max_C_HT = C_HT


             allocate( psi1(1:Fdof, 1:nbDim, 1:Gdof1) )
             psi1(:, :, :) = 0.

             if( elem1%HGnode) then
               HGidx1 = elem1%HGface(2,ie1)

               if (elem1%HGface(2, ie1) > 1) then
                  call Eval_Loc_RTN_Edge_HG(elem1, loc_RTN, ie1, HGidx1, psi1(1:Fdof, 1:nbDim, 1:Gdof1 ), .true. )
               else
                  allocate( Rpsi1(1:Fdof, 1:nbDim, 1: Gdof1) )
                  Rpsi1(:, :, :) = 0.
                  call Eval_Loc_RTN_Edge(elem1, loc_RTN, ie1, Rpsi1(1:Fdof, 1:nbDim, 1:Gdof1 ) )  !no possibility true/false in Eval_Loc_RTN_Edge

                  ! reordering
                  do l=1, Gdof1
                     psi1(1:Fdof, 1:nbDim, l) = Rpsi1(1:Fdof, 1:nbDim, Gdof1 -l + 1)
                  enddo
                  deallocate(Rpsi1)
               endif

                              else
               allocate( Rpsi1(1:Fdof, 1:nbDim, 1: Gdof1) )
               Rpsi1(:, :, :) = 0.
               call Eval_Loc_RTN_Edge(elem1, loc_RTN, ie1, Rpsi1(1:Fdof, 1:nbDim, 1:Gdof1 ) )  !no possibility true/false in Eval_Loc_RTN_Edge

               ! reordering
               do l=1, Gdof1
                  psi1(1:Fdof, 1:nbDim, l) = Rpsi1(1:Fdof, 1:nbDim, Gdof1 -l + 1)
               enddo
               deallocate(Rpsi1)
             endif


             Mdof = min( Gdof1, Gdof )

             if (elem%i < 4 .and. (state%time%iter == 5 .or. state%time%iter == 6) ) then
             !do l=1, Fdof
             !   write(*, '(a5, i2, a20, 12es14.6)') 'psi1(', l, ' 1, 1:Mdof)', psi1(l, 1, 1:Mdof)
             !   write(*, '(a5, i2, a20, 12es14.6)') 'psi1(', l, ' 2, 1:Mdof)', psi1(l, 2, 1:Mdof)
             !   print*,''
             !   write(*, '(a5, i2, a20, 12es14.6)') 'psi(', l, ' 1, 1:Mdof)', psi(l, 1, 1:Mdof)
             !   write(*, '(a5, i2, a20, 12es14.6)') 'psi(', l, ' 2, 1:Mdof)', psi(l, 2, 1:Mdof)
             !   pause
             !enddo
             !write(*,*) elem%HGnode
             !if (elem%i == 1) then
             !    temp(1, 1, 1:Mdof) =  elem%RTNflux(1, 1, 4)*psi(4, 1, 1:Mdof)
             !    temp(1, 2, 1:Mdof) =  elem1%RTNflux(1, 1, 7)*psi1(7, 1, 1:Mdof)
             !    write(*, '(a10, 15es14.6)'), 'rozdil odpov 1 mom', temp(1, 1, 1:Mdof) - temp(1, 2, 1:Mdof)
             !endif
             !pause
             !!print*,'********************'
             !!write(*,'(a8, i2, a4, i2, a14, 2es10.2)') 'elem%i=', elem%i, 'ie=', ie, 'elem%xc(1:2)', elem%xc(1:2)
             !!write(*,'(a8, i2, a4, i2, a14, 2es10.2)') 'elem1%i=', elem1%i, 'ie1=', ie1, 'elem1%xc(1:2)', elem1%xc(1:2)
             !!write(*, '(a40, 15es14.6)'), 'elem%RTNflux(1, 1, 1:Fdof)', elem%RTNflux(1, 1, 1:Fdof)
             !!write(*, '(a40, 15es14.6)'), 'elem1%RTNflux(1, 1, 1:Fdof)', elem1%RTNflux(1, 1, 1:Fdof)
             !do l=1, Fdof
             !   write(*, '(a4, 15es14.6)'), 'psi', psi(l, 2, 1:Mdof)
             !   write(*, '(a4, 15es14.6)'), 'psi1', psi1(l, 2, 1:Mdof)
             !enddo
             if (disc) then
                !print*, ' '
                !print*, ' '
                !print*, 'elem%i=', elem%i
                !print*, 'elem1%i=', elem1%i
                !temp(1, 1, 1:Mdof) = matmul(elem%RTNflux(1, 1, 1:Fdof), psi(1:Fdof, 1, 1:Mdof) )
                !temp(1, 2, 1:Mdof) = matmul(elem%RTNflux(1, 1, 1:Fdof), psi(1:Fdof, 2, 1:Mdof) )
                !temp(1, 1, 1:Mdof) = temp(1, 1, 1:Mdof) * elem%n(ie, 1) + temp(1, 2, 1:Mdof) * elem%n(ie, 2)
                !write(*, '(a40, 12es14.6)') 'd_h * n L', temp(1, 1, 1:Mdof)
                !print*, ' '
                !!!write(*, '(a28, 15es14.6)') 'elem1%RTNflux(1, 1, 1:Fdof)', elem1%RTNflux(1, 1, 1:Fdof)
                !temp(1, 1, 1:Mdof) = matmul(elem1%RTNflux(1, 1, 1:Fdof), psi1(1:Fdof, 1, 1:Mdof) )
                !temp(1, 2, 1:Mdof) = matmul(elem1%RTNflux(1, 1, 1:Fdof), psi1(1:Fdof, 2, 1:Mdof) )
                !temp(1, 1, 1:Mdof) = temp(1, 1, 1:Mdof) * elem%n(ie, 1) + temp(1, 2, 1:Mdof) * elem%n(ie, 2)
                !write(*, '(a40, 12es14.6)') 'd_h * n R', temp(1, 1, 1:Mdof)
                !print*, ' '
                !pause
                        else
                !print*, ' '
                !print*, ' '
                !print*, 'elem%i=', elem%i
                !print*, 'elem1%i=', elem1%i
                !temp(1, 1, 1:Mdof) = matmul(elem%RTNflux(0, 1, 1:Fdof) - elem%RTNflux(1, 1, 1:Fdof), psi(1:Fdof, 1, 1:Mdof) )
                !temp(1, 2, 1:Mdof) = matmul(elem%RTNflux(0, 1, 1:Fdof) - elem%RTNflux(1, 1, 1:Fdof), psi(1:Fdof, 2, 1:Mdof) )
                !temp(1, 1, 1:Mdof) = temp(1, 1, 1:Mdof) * elem%n(ie, 1) + temp(1, 2, 1:Mdof) * elem%n(ie, 2)
                !write(*, '(a40, 12es14.6)') '(d_h^nu - d_h) * n L', temp(1, 1, 1:Mdof)
                !print*, ' '
                !temp(1, 1, 1:Mdof) = matmul(elem1%RTNflux(0, 1, 1:Fdof) - elem1%RTNflux(1, 1, 1:Fdof), psi1(1:Fdof, 1, 1:Mdof) )
                !temp(1, 2, 1:Mdof) = matmul(elem1%RTNflux(0, 1, 1:Fdof) - elem1%RTNflux(1, 1, 1:Fdof), &
                !                                                                          psi1(1:Fdof, 2, 1:Mdof) )
                !temp(1, 1, 1:Mdof) = temp(1, 1, 1:Mdof) * elem%n(ie, 1) + temp(1, 2, 1:Mdof) * elem%n(ie, 2)
                !write(*, '(a40, 12es14.6)') '(d_h^nu - d_h) * n R', temp(1, 1, 1:Mdof)
                !print*, ' '
                !pause
               endif
             endif

             do k=1, ndim
                temp(k, 1:nbDim, 1:maxval(elem%face(fGdof,:) )) = 0.  !temp allocated for 1:maxval(elem%face(fGdof,:) )
                do l=1,nbDim
                   if (disc) then  ! discretization flux nonconformity
                      temp(k, l, 1:Mdof) = matmul(elem%RTNflux(1, k, 1:Fdof), psi(1:Fdof, l, 1:Mdof) ) &
                                           -  matmul(elem1%RTNflux(1, k, 1:Fdof), psi1(1:Fdof, l, 1:Mdof) ) !elem1 RTNflux at previous lvl saved in elem1%RTNflux(1,..)

                      !write(*, '(a30, i2, a5, 6es14.6)') 'temp(1, 1:elem%face(fGdof,', Tface, ')) =', temp(1, l, 1:Mdof)
                      !print*,''
                      !write(*, '(a27, 15es14.6)') 'elem1%RTNflux(1, 1, 1:Fdof)', elem1%RTNflux(1, 1, 1:Fdof)
                      !print*,''

                   else ! algebraic flux nonconformity
                      allocate( RTNfluxalg(1:ndim, 1:Fdof), RTNfluxalg1(1:ndim, 1:Fdof) )
                      RTNfluxalg(k, 1:Fdof) = elem%RTNflux(0, k, 1:Fdof) - elem%RTNflux(1, k, 1:Fdof)
                      RTNfluxalg1(k, 1:Fdof) = elem1%RTNflux(0, k, 1:Fdof) - elem1%RTNflux(1, k, 1:Fdof)

                      if (elem%i == 1) then
                         !write(*, '(a24, 15es14.6)') 'RTNfluxalg(k, 1:Fdof) =', RTNfluxalg(k, 1:Fdof)
                         !print*,''
                         !write(*, '(a24, 15es14.6)')  'RTNfluxalg1(k, 1:Fdof) =', RTNfluxalg1(k, 1:Fdof)
                         !print*,''
                      endif

                      temp(k, l, 1:Mdof) = matmul(RTNfluxalg(k, 1:Fdof), psi(1:Fdof, l, 1:Mdof) ) &
                                           -  matmul(RTNfluxalg1(k, 1:Fdof), psi1(1:Fdof, l, 1:Mdof) )
                      deallocate(RTNfluxalg, RTNfluxalg1)
                   endif

                enddo ! l


             enddo ! k



             !if (elem%i < 3) then
             !if (Tface < 2) then
                !write(*, '(a9, i2)') 'ie = ', ie
                !write(*, '(a9, i2)') 'Tface = ', Tface
                !write(*, '(a27, 10es14.6)') 'temp(1, 1, 1:Mdof)', temp(1, 1, 1:Mdof)
                !write(*, '(a27, 10es14.6)') 'temp(1, 2, 1:Mdof)', temp(1, 2, 1:Mdof)
                !write(*, '(a27, 10es14.6)') 'elem%RTNflux(1, 1, 1:Fdof)', elem%RTNflux(1, 1, 1:Fdof)
                !if (elem%RTNfluxeval) then
                !   write(*, '(a28, 10es14.6)') 'elem1%RTNflux(1, 1, 1:Fdof)', elem1%RTNflux(1, 1, 1:Fdof)
                !else
                !   write(*, '(a28, 10es14.6)') 'elem1%RTNflux(0, 1, 1:Fdof)', elem1%RTNflux(0, 1, 1:Fdof)
                !endif
                !write(*, '(a9, i2)') 'mdof = ', mdof
                !do l=1, Fdof
                !   write(*, '(a5, i2, a20, 12es14.6)') 'psi(', l, ' 1, 1:Mdof), psi1', psi(l, 1, 1:Mdof), psi1(l, 1, 1:Mdof)
                !   write(*, '(a5, i2, a20, 12es14.6)') 'psi(', l, ' 2, 1:Mdof), psi1', psi(l, 2, 1:Mdof), psi1(l, 2, 1:Mdof)
                !   print*,'                                   '
                !enddo
                !print*,'                                   '
                !extract(1, 1, 1:Mdof) = matmul(elem%RTNflux(1, 1, 1:Fdof), psi(1:Fdof, 1, 1:Mdof) )
                !write(*, '(a30, 12es14.6)') 'theta(1,:)psi(:, 1:Mdof) =', extract(1, 1, 1:Mdof)
                !extract(1, 1, 1:Mdof) = matmul(elem1%RTNflux(1, 1, 1:Fdof), psi1(1:Fdof, 1, 1:Mdof) )
                !write(*, '(a30, 12es14.6)') 'theta1(1,:)psi1(:, 1:Mdof) =', extract(1, 1, 1:Mdof)
                !print*,''
                !extract(1, 2, 1:Mdof) = matmul(elem%RTNflux(1, 1, 1:Fdof), psi(1:Fdof, 2, 1:Mdof) )
                !write(*, '(a30, 12es14.6)') 'theta(2,:)psi(:, 1:Mdof) =', extract(1, 2, 1:Mdof)
                !extract(1, 2, 1:Mdof) = matmul(elem1%RTNflux(1, 1, 1:Fdof), psi1(1:Fdof, 2, 1:Mdof) )
                !write(*, '(a30, 12es14.6)') 'theta1(2,:)psi1(:, 1:Mdof) =', extract(1, 2, 1:Mdof)
                !pause
                !print*,'                                   '
                !temp(1, 2, 1:Mdof) = matmul(elem%RTNflux(1, 1, 1:Fdof), psi(1:Fdof, 2, 1:Mdof) )
                !write(*, '(a30, i2, a5, 12es14.6)') 'temp(1, 1, 1:elem%face(fGdof,', 2, ')) =', temp(1, 2, 1:Mdof)
                !temp(1, 2, 1:Mdof) = matmul(elem1%RTNflux(1, 1, 1:Fdof), psi1(1:Fdof, 2, 1:Mdof) )
                !write(*, '(a30, i2, a5, 12es14.6)') 'temp(1, 2, 1:elem%face(fGdof,', 2, ')) =', temp(1, 2, 1:Mdof)
             !endif
             !endif

             deallocate( psi1 )

             omega_edge = 0.5

          else  !! boundary face

             Mdof = Gdof

             do k=1, ndim
                temp(k, 1:nbDim, 1:maxval(elem%face(fGdof,:) )) = 0.  !temp allocated for 1:maxval(elem%face(fGdof,:) )
             enddo ! k

             omega_edge = 1.

          endif

          !write(*, '(a30, i2, a5, i2, a10, i2)') 'elem%face(fGdof,', Tface, ') =', elem%face(fGdof,Tface), 'Gdof = ', Gdof

          if (elem%face(fGdof,Tface) /= Mdof ) print*,'## Trouble in FluxNCElemDiscEstimSS with Mdof.'

          !write(*, '(a30, i2, a5, 8es14.6)') 'temp(1, 1, 1:elem%face(fGdof,', Tface, ')) =', temp(1, 1, 1:elem%face(fGdof,Tface))
          !write(*, '(a30, i2, a5, 8es14.6)') 'temp(1, 2, 1:elem%face(fGdof,', Tface, ')) =', temp(1, 2, 1:elem%face(fGdof,Tface))
          !write(*,'(a5, i2)') 'nbDim', nbDim

          do k=1, ndim
             call IntegrateSquareFunctionNormalEdge(elem, Tface, temp(k, 1:nbDim, 1:elem%face(fGdof,Tface)), val(k))
          enddo ! k


          if( elem%HGnode) then
             if (elem%HGface(2,Tface) == 1) then
                 !write(*, '(a9, i2)') 'Tface = ', Tface

                 !write(*, '(a30, i2, a5, es14.6)') '\|[{\bf d}_{h}^i {\cdot} {\bf n}] \|_{', Tface, '} = ', val(1:ndim)**0.5

             endif
          endif

          if(.not. elem%HGnode) then

            !write(*, '(a9, i2)') 'Tface = ', Tface

            !write(*, '(a30, i2, a5, es14.6)') '\|[{\bf d}_{h}^i {\cdot} {\bf n}] \|_{', Tface, '} = ', val(1:ndim)**0.5

          endif

          etaFNCDie = etaFNCDie + val(1:ndim)**0.5

          !write(*, '(a10, i2, a3, es14.6)') 'etaFNCD', ie, ' = ', etaFNCDie

          !if(elem%i == 2) &
          !     write(*,'(a6,2i5,20es14.6)') &
          !     'Rflux=',elem%i,elem%face(neigh,Tface),-Rflux(1:Gdof, 1:ndim)



          !if(elem%i == 2) &
          !     write(*,'(a6,2i5,20es14.6)') 'Cflux=', &
          !     elem%i,elem%face(neigh,Tface),Cflux(1:Gdof, 1:ndim)



          !if(elem%i == 2) &
          !     write(*,'(a6,2i5,20es14.6)') 'jump=',elem%i,elem%face(neigh,Tface), &
          !     jump_wi(1:Gdof, 1:ndim), state%space%sigma * state%model%Re1


       enddo !Tface

       C_HT = max_C_HT

       call BoundTraceIneq(elem, ie, HGidx, C_edgTb)

       !write(*, '(a10, es14.6)') 'C_edgTb', C_edgTb

       !write(*, '(a10, F3.1)') 'omega_edge', omega_edge

       !write(*, '(a18, es14.6)') 'elem%dn(ie)**0.5', elem%dn(ie)**0.5

       !!!etaFNCDie = C_edgTb * omega_edge * (elem%dn(ie)**0.5) * etaFNCDie

      ! write(*, '(a10, i3)') 'elem =', elem%i

       !write(*, '(a10, es14.6)') 'C_edgTb =', C_edgTb

       !write(*, '(a10, i2, a11, es14.6)') 'etaFNCDie', ie, '* const = ', etaFNCDie

       etaFNCD = etaFNCD + etaFNCDie

    enddo  ! ie


    deallocate(temp, extract, psi)

  end subroutine FluxNCElemDiscEstimSS


  !> evaluation of algebraic remainder estimator
  !> \f$ \eta_{rem,T}^i := c_{F,Omega} \|r_h^{i+nu}\|_{T} \f$
  subroutine RemainderElemAlgEstimSS(elem, etaRem)
    type(element), intent(in) :: elem
    real, dimension(1:ndim), intent(inout) :: etaRem
    real, dimension(1:ndim) :: val
    integer :: Qdof
    real :: C_FOm

    val(:) = 0.
    Qdof = elem%Qdof
    call IntegrateSquareVectorFunction(elem, elem%res_func(1:ndim, 1:Qdof), val(1:ndim) )

    !Friedrichs' constant
    !C_FOm = (1./Pi) * &
    !     ( 1./(maxval(grid%x(:,1)) - minval(grid%x(:,1) ) )&
    !     + 1./(maxval(grid%x(:,2)) - minval(grid%x(:,2) ) ) )**0.5

    C_FOm = (1./Pi) * (1./nbDim)**0.5  ! Friedrichs' constant, 1 = length of domain edges,  FOR DIFFERENT DOMAIN HAS TO BE MODIFIED!!!
    !write(*,'(a7, es12.4)') 'C_FOm =', C_FOm
    !write(*,'(a9, es12.4)') 'state%space%h =', state%space%h

    !write(*,'(a20, es12.4)') '\|r_h^{i+nu}\|_{T} =', val(1:ndim)**0.5

    etaRem(1:ndim) = C_FOm * (val(1:ndim)**0.5)

  end subroutine RemainderElemAlgEstimSS


  !> evaluation of potential nonconformity estimator
  !> \f$ \eta_{PNC,T}^i := \|\nabla (u_h^{i} - I_{AV}(u_h^{i}) )\|_{T} \f$
  subroutine PotentNCElemEstimSS(Set_R_s, elem, etaPNC)
    interface
      subroutine Set_R_s(ndimL, nbDim, Qdof, w, Dw, Re_1, R_s)
         integer, intent(in) :: ndimL, nbDim, Qdof
         real, dimension(1:Qdof, 1:ndimL), intent(in):: w !state  w in #Qdof nodes
         real, dimension(1:Qdof, 1:ndimL, 1:nbDim), intent(in):: Dw !state  Dw in #Qdof nodes
         real, dimension(1:Qdof), intent(in) :: Re_1        ! inverse of Reynolds number
         real, dimension(1:Qdof, 1:nbDim, 1:ndimL), intent(inout) :: R_s
      end subroutine Set_R_s

    end interface

    type(element), intent(in) :: elem
    real, dimension(1:ndim), intent(inout) :: etaPNC
    real, dimension(:,:,:,:), allocatable :: Oswald  ! Oswald inter in Lagr. nodes
    real, dimension(:,:), allocatable :: potential  ! potential reconstr. in basis
    real, dimension(:,:,:), allocatable :: R_s
    real, dimension(:,:), allocatable :: temp
    real, dimension(:,:,:), allocatable :: Dpotent
    integer ::  Rdeg, Rdof, Qnum, Qdof, k, l
    real, dimension(1:ndim) :: val


    ! degree of polynomial reconstruction
    Rdeg = maxval(grid%elem(:)%deg)
    Rdof = (Rdeg+1)*(Rdeg+2)/2

    !Qnum = 2*Rdeg
    Qnum = state%space%Qdeg(Rdeg,1)

    Qdof = elem%Qdof

    ! construction of the Oswald, contains values in Lagr. nodes
    allocate(Oswald(0:1,1:grid%nelem,1:ndim, 1:Rdof) )

    ! Oswald at i-th time level
    call ConstructOswald(1, Rdeg, Rdof, Oswald(1, 1:grid%nelem, 1:ndim, 1:Rdof))

    allocate(potential(1:ndim, 1:Rdof) )

    ! recomputation of the Oswald interpolation from Lagr. nodes into basis
    call Lagr2BasisDiff(elem, Rdeg,  Oswald(1, elem%i, 1:ndim, 1:Rdof), Qnum, Rdof, &
         potential(1:ndim, 1:Rdof) )

    allocate(R_s(1:Qdof, 1:nbDim, 1:ndim) )

    call Eval_R_s_Elem(Set_R_s, elem, R_s(1:Qdof, 1:nbDim, 1:ndim) )

    allocate( temp(1:ndim, 1:Qdof) )

    allocate(Dpotent(1:ndim, 1:Qdof, 1:nbDim ) )

    ! evaluation of the gradient of potential in integ. nodes
    call Eval_DVec_Elem(elem, Rdof, potential(1:ndim, 1:Rdof), &
         Dpotent(1:ndim, 1:Qdof, 1:nbDim) )

    do k=1, ndim
          temp(k, 1:Qdof) = 0.
          do l=1,nbDim
             temp(k, 1:Qdof) = temp(k, 1:Qdof) + &
                  (R_s(1:Qdof, l, k) - Dpotent(k, 1:Qdof, l ) )**2
          enddo !l

    enddo !k


    call IntegrateVectorFunction(elem, temp(1:ndim, 1:Qdof), val(1:ndim))
    etaPNC(1:ndim) = val(1:ndim)**0.5 !SS

    !write(*,'(a50, es12.4)') '\|\nabla (u_h^{i} - I_{AV}(u_h^{i}) )\|_{T} =', val(1:ndim)**0.5

    deallocate(Oswald, R_s, potential, temp, Dpotent)


  end subroutine PotentNCElemEstimSS


  subroutine BoundTraceIneq(elem, ie, HGidx, C_edgTb)
    type(element), intent(in) :: elem
    integer, intent(in) :: ie  ! index of the edge of the simplex elem
    integer, intent(in) :: HGidx ! whether or not HG are present at ie
    real, intent(out) :: C_edgTb
    real :: rohat ! the diameter of the inscribed ball of reference simplex
    real :: C_sd, alpha, C_HT

    if (nbDim == 2) then

      alpha = 0.730276      !! see p. 1486 in "A POSTERIORI ERROR ESTIMATIONS OF SOME CELL-CENTERED FINITE VOLUME METHODS"

      rohat = 2./(2 + (2**0.5))

      !write(*, '(a8, es12.5)') 'rohat = ', rohat

    else
      alpha = (2 * (11 + 4*(6**0.5) ) )**0.25

      !write(*, '(a8, es12.5)') 'alpha = ', alpha

      alpha = alpha / ((Pi * tanh(Pi) )**0.5 )

      !write(*, '(a8, es12.5)') 'alpha = ', alpha

      rohat = 2. * 1. * (6**0.5)/12.

      !write(*, '(a8, es12.5)') 'rohat = ', rohat

    endif


    C_sd = (alpha * (1./rohat) * (1./nbDim**0.5) )**2

    !write(*, '(a8, es12.5)') 'C_sd = ', C_sd

    C_edgTb = elem%dn(ie) * (elem%diam**2)
    C_edgTb = C_edgTb / (elem%area *  elem%dn(ie) )
    C_edgTb = (C_edgTb * C_sd)**0.5

    !write(*, '(i2, a14, es12.5, a14, es12.5)') elem%i, 'elem%diam = ', elem%diam, 'elem%area = ',elem%area

    if (HGidx == 1) then !! ie is the only one element edge on the triangle edge

        !write(*, '(a10, i2, a12, es12.5)') 'HGidx = ', HGidx, 'C_edgTb = ', C_edgTb
    else


        C_edgTb = C_edgTb * (C_HT**((nbDim + 1)/2.) ) ! * maxval(...) is missing !!, C_HT is not evaluated yet

        !write(*, '(a10, i2, a12, es12.5)') 'HGidx = ', HGidx, 'C_edgTb = ', C_edgTb

        !print*, '***********Sub. BoundTraceIneq************'
        !write(*, '(a10, i3)') 'elem =', elem%i
        !write(*, '(a10, es12.5)') 'C_HT =', C_HT
        !write(*, '(a30, es12.5)') 'C_HT**((nbDim + 1)/2.) =', C_HT**((nbDim + 1)/2.)
        !write(*, '(a10, es12.5)') 'C_edgTb =', C_edgTb
        !pause

    endif


  end subroutine BoundTraceIneq





end module alg_estim
