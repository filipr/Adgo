!> set the regularity parameters for adaptation
module regularity
  use mesh_oper
  use main_data
  use data_mod
  use integration
  use blocks_integ
  use model_oper
  use basis
  use f_mapping
  use eval_sol

  implicit none

  public:: Set_Mesh_regularity
  public:: Set_Elem_regularity
  public:: Regularity_hp_derefinement

contains
  !>  set the regularity parameters for adaptation for all mesh
  subroutine Set_Mesh_regularity( )
    class(element), pointer:: elem      ! elem = element
    integer :: i

    do i=1, grid%nelem
       elem => grid%elem(i)
       call Set_Elem_regularity(elem)
    enddo
  end subroutine Set_Mesh_regularity

       
  !>  set the regularity parameters for adaptation on element
  subroutine Set_Elem_regularity(elem)
    type(element), intent(inout):: elem      ! elem = element
    integer ::  l1, l2, rn 
    real :: param, val, val1, qq, estPP
    real :: err_h_deref, err_p_deref, slope_h, slope_p
    logical :: singularity
    real, dimension(:,:), allocatable :: xi, yi
    real :: sumx, sumy, sumxy, sumxx, sumyy, an, bn

    ! indicator of the regularity
    elem%reg =  1E-15
    elem%reg1 = 1E-15
    elem%reg2 = 1E-15
    
    param = state%space%adapt%type_regularity_par

    if( state%space%estim_space == 'pNeu') then
       ! variant used in the 1st revision of SISC
       !if(elem%deg >= 1) elem%reg  = elem%eta(P_tot, 1) / elem%eta(P_potP, 1) !/  elem%diam**0.5

       ! testing variant for the 1st revision of SISC
       if(elem%deg >= 1) then
          call Regularity_only_p_derefinement(elem, err_p_deref)
          elem%reg  = elem%eta(P_tot, 1) / err_p_deref
          elem%eta(P_potP, 1) = err_p_deref
          !write(*,'(a10, i5, 30es12.4)') &
          !     '#####DE#',elem%i,  &
          !     elem%eta(P_tot, 1) / elem%eta(P_potP, 1), &
          !     elem%eta(P_tot, 1) / err_p_deref, &
          !     elem%eta(P_potP, 1), err_p_deref  ,  elem%eta(P_potP, 1) / err_p_deref 
       endif


    elseif( state%space%estim_space == 'RES' .or.  state%space%estim_space == 'inter' ) then
       elem%reg = 1.

    else
       stop 'non0,plimented estim_space  in regularity.f90 8343yfdh843id'
    endif


    !Setting of elem%regT0:    elem%reg <  elem%regT0  ==> p-refinement
    
    select case(state%space%adapt%type_regularity )

    case(-1) 
       if(elem%i == 1) print*,'! ideal hp-refinement based on the known analyticity '
       call Detect_apriori_known_singularity(elem, singularity) 
     
       if(singularity) then
          elem%regT0 = elem%reg / 2.  ! h-refinement
       else
          elem%regT0 = elem%reg * 2.  ! p-refinement
       endif

    case(0)  ! DEFAULT
       elem%regT0 = 0.5 

    case(1)  
       if(elem%i == 1) print*,'! Babuska indicator, par =', param
       elem%regT0 = param   ! BABUSKA
       !elem%regT0 = 0.3   ! BABUSKA

    case(2) 
       if(elem%i == 1) print*,'! DEV - Dolejsi, Ern, Vohralik, par =', param

       elem%regT0 = param* elem%diam**0.5  ! param=0.75 

       !elem%regT0 = (elem%deg -1. )/ elem%deg * elem%diam**0.5  ! CAN WORKS nicely, VERIFY
       !if(elem%deg == 1) elem%regT0 = elem%reg * 2  ! p -refinement 

    case(3)
       if(elem%i == 1) print*,'!PRIOR2P indicator, see [Mitchel, McClain]'
       if(elem%deg > 1) then
          val1 = 1. &
               - log(elem%eta(P_tot, 1) / elem%eta(P_potP, 1)) / log(1.*elem%deg / (elem%deg-1))

          if(val1 < elem%deg + 1) then
              elem%regT0 = elem%reg / 2.  ! h-refinement
           else
              elem%regT0 = elem%reg * 2.  ! p-refinement
           endif
          !write(*,'(a20, 4i5, 10es12.4)') &
          !     'elem_regu:',elem%i, elem%deg, int(val+0.5), int(val1+0.5), val, val1, &
          !     sqrt(dot_product(elem%xc(:), elem%xc(:) ))
        else
          elem%regT0 = elem%reg * 2.  ! p-refinement
       endif

    case(4)

       if(elem%i == 1) print*,'! ANALYTICITY test Houston'
       if(elem%deg == 1) then
           elem%regT0 = elem%reg * 2.  ! p-refinement
        else
           l1 = elem%dof - elem%deg 
           l2 = elem%dof
           val = dot_product(elem%w(0, l1:l2), elem%w(0, l1:l2) ) / elem%area**2
           val1 = log((2.*elem%deg + 1)/(2*val)) / 2 / log(1.*elem%deg) 

           if(elem%deg <= val1 -1) then
              elem%regT0 = elem%reg * 2.  ! p-refinement
           else
              elem%regT0 = elem%reg / 2.  ! h-refinement
           endif
           !write(*,'(a8, 6es12.4)') 'regul:', val, val1, elem%xc(:)
        endif

     case(5)  ! approach close to [Dolejsi, Solin AMC 2016??]
        if(elem%i == 1) print*,'! Approach minimal DoF'
        qq =  (elem%eta(P_tot, 1) / elem%eta(P_potP, 1)  )**(2./(elem%deg -1) )
        val = (elem%deg + 2.)/(elem%deg ) * qq
        
        !write(*,'(a7, i5, 6es12.4)') 'regul:',elem%i, elem%eta(P_tot, 1) / elem%eta(P_potP, 1), qq, val

        !if(val <= 1.) then  
        !   elem%regT0 = elem%reg * 2.  ! p-refinement
        !   !write(*,'(a8, 6es12.4)') 'regul:', val, elem%xc(:)
        !else
        !   elem%regT0 = elem%reg / 2.  ! h-refinement
        !   !write(*,'(a8, 6es12.4)') 'regul:', val, elem%xc(:), sqrt(dot_product(elem%xc(:), elem%xc(:)))
        !endif

        elem%regT0 = 0.95 * (1. * elem%deg /(elem%deg +2)) **((elem%deg -1) / 2.)

        !if(elem%reg >= elem%regT0) &
        !     write(*,'(a7, 2i5, 3(3es12.4, a2 ))') 'regul??',elem%i,  elem%deg, &
        !     elem%eta(P_tot, 1) / elem%eta(P_potP, 1), qq, val, "|", &
        !     elem%reg, elem%regT0, elem%reg / elem%regT0,"|", sqrt(dot_product(elem%xc(:), elem%xc(:) ))


    case(6)
       if(elem%i == 1) print*,'!modified PRIOR2P regularity indicator'
       if(elem%deg > 1) then
          val = 1. &
               - log(elem%eta(P_tot, 1) / elem%eta(P_potP, 1)/elem%diam)/ log(1.*elem%deg /(elem%deg-1))

          val1 = 1. &
               - log(elem%eta(P_tot, 1) / elem%eta(P_potP, 1) ) / log(1.*elem%deg /(elem%deg-1))

          if(val1 <  elem%deg + 1) then
              elem%regT0 = elem%reg / 2.  ! h-refinement
           else
              elem%regT0 = elem%reg * 2.  ! p-refinement
           endif

           !write(*,'(a20, 4i5, 3es12.4,a3,10es12.4)') &
           !    'elem_regu:',elem%i, elem%deg, int(val+0.5), int(val1+0.5), val, val1, &
           !    sqrt(dot_product(elem%xc(:), elem%xc(:) )),'|||', &
           !    log(elem%eta(P_tot, 1) / elem%eta(P_potP, 1) ) , &
           !    log(elem%eta(P_tot, 1) / elem%eta(P_potP, 1)/elem%diam), log(1.*elem%deg /(elem%deg-1))
        else
          elem%regT0 = elem%reg * 2.  ! p-refinement
       endif

    case(7)
       if(elem%i == 1) print*,'! minimal DoF extrapolated to p+1'
       if(elem%deg > 1) then

          val =  log(elem%eta(P_tot, 1) / elem%eta(P_potP, 1) ) &
               * log(1.*elem%deg /(elem%deg+1)) / log(1.*(elem%deg-1) /elem%deg  )

          estPP = elem%eta(P_tot, 1) * exp(val)

                  
        !write(*,'(a7, i5, 6es12.4)') 'regul:',elem%i, elem%eta(P_tot, 1) / elem%eta(P_potP, 1), qq, val

        !if(val <= 1.) then  
        !   elem%regT0 = elem%reg * 2.  ! p-refinement
        !   !write(*,'(a8, 6es12.4)') 'regul:', val, elem%xc(:)
        !else
        !   elem%regT0 = elem%reg / 2.  ! h-refinement
        !   !write(*,'(a8, 6es12.4)') 'regul:', val, elem%xc(:), sqrt(dot_product(elem%xc(:), elem%xc(:)))
        !endif

        elem%regT0 =  (1. * (elem%deg+1) /(elem%deg +3)) **((elem%deg ) / 2.)

        elem%reg  = estPP/ elem%eta(P_tot, 1)
        
        !write(*,'(a20, 4i5, 3es12.4,a3,10es12.4)') &
        !     'elem_regu:',elem%i, elem%deg, int(val+0.5), int(val1+0.5), &
        !     elem%eta(P_potP, 1), elem%eta(P_tot, 1) , estPP , &
        !     '|||', elem%reg, elem%regT0, estPP/ elem%eta(P_tot, 1)

        else
          elem%regT0 = elem%reg * 2.  ! p-refinement
       endif

    case(8)
       if(elem%i == 1) print*,'! minimal DoF extrapolated to p+1 from 3-steps'

       if(elem%deg > 2) then
          ! linear regreation of the error estims, we assume that e = c * p^s
          rn = 3
          allocate(xi(1:rn, 1:2), yi(1:rn, 1:2) )
          
          xi(1, 1) = elem%deg - 2
          xi(2, 1) = elem%deg - 1
          xi(3, 1) = elem%deg  
          xi(1, 2) =  elem%eta(P_potPP, 1)
          xi(2, 2) =  elem%eta(P_potP, 1)
          xi(3, 2) =  elem%eta(P_tot, 1)

          yi(1:rn, 1:2)  = log( xi(1:rn, 1:2)  )

          sumx = sum(yi(1:rn, 1))
          sumy = sum(yi(1:rn, 2))

          sumxx = sum(yi(1:rn, 1)*yi(1:rn, 1))
          sumxy = sum(yi(1:rn, 1)*yi(1:rn, 2))
          !sumyy = sum(yi(1:rn, 2)*yi(1:rn, 2))

          an = (sumxy - sumx * sumy / rn) /(sumxx - sumx*sumx/ rn )
          bn = (sumy - an * sumx) / rn
          

          estPP = exp(bn + an * log(elem%deg + 1.) )

          ! write(*,'(a7, i5,es12.4, a3,  6es12.4)') 'regul:',elem%i, estPP, '|', &
          !      elem%eta(P_tot, 1), elem%eta(P_potP, 1),  elem%eta(P_potPP, 1)

          ! if( sqrt(dot_product(elem%xc, elem%xc)) < 10.25 ) then
          !    l1 = 40+state%space%adapt%adapt_level
          !    l2 = 30+state%space%adapt%adapt_level
          !    write(l2,*)  elem%xc,estPP
          !    write(l1,*)  elem%xc,estPP
          !    write(l1,*)  elem%xc, elem%eta(P_tot, 1)
          !    write(l1,*)  elem%xc, elem%eta(P_potP, 1) 
          !    write(l1,*)  elem%xc, elem%eta(P_potPP, 1), &
          !         elem%eta(P_potP, 1)/ elem%eta(P_tot, 1),  elem%eta(P_potPP, 1) / elem%eta(P_potP, 1)
          !    write(l1, '(x)') 
          !    write(l1, '(x)') 
          !    write(l1, '(x)') 
          ! endif

          !stop "3e3d"
          deallocate(xi,yi)

          ! setting of the type of the refinement based on the minimal DOF
          elem%regT0 =  (1. * (elem%deg+1) /(elem%deg +3)) **((elem%deg ) / 2.)

          elem%reg  = estPP/ elem%eta(P_tot, 1)

       !endif


       elseif(elem%deg == 2) then

          val =  log(elem%eta(P_tot, 1) / elem%eta(P_potP, 1) ) &
               * log(1.*elem%deg /(elem%deg+1)) / log(1.*(elem%deg-1) /elem%deg  )

          estPP = elem%eta(P_tot, 1) * exp(val)

                  
          !write(*,'(a7, i5, 6es12.4)') 'regul:',elem%i, elem%eta(P_tot, 1) / elem%eta(P_potP, 1), qq, val
          
          !if(val <= 1.) then  
          !   elem%regT0 = elem%reg * 2.  ! p-refinement
          !   !write(*,'(a8, 6es12.4)') 'regul:', val, elem%xc(:)
          !else
          !   elem%regT0 = elem%reg / 2.  ! h-refinement
          !   !write(*,'(a8, 6es12.4)') 'regul:', val, elem%xc(:), sqrt(dot_product(elem%xc(:), elem%xc(:)))
          !endif
          
          elem%regT0 =  (1. * (elem%deg+1) /(elem%deg +3)) **((elem%deg ) / 2.)
          
          elem%reg  = estPP/ elem%eta(P_tot, 1)
          
          !write(*,'(a20, 4i5, 3es12.4,a3,10es12.4)') &
          !     'elem_regu:',elem%i, elem%deg, int(val+0.5), int(val1+0.5), &
          !     elem%eta(P_potP, 1), elem%eta(P_tot, 1) , estPP , &
          !     '|||', elem%reg, elem%regT0, estPP/ elem%eta(P_tot, 1)
          
       else
          elem%regT0 = elem%reg * 2.  ! p-refinement
       endif

    case(9)
       if(elem%i == 1) print*,'! maximal exponential degree  of convergence'
       if(elem%deg > 1) then
          
          rn = 3
          allocate(xi(1:rn, 1:2) )

          ! p-1, h
          xi(1, 1) = elem%eta(P_potP, 1)                    ! error
          xi(1, 2) = (elem%deg + 0) * (elem%deg + 1) / 2    ! DOF

          ! p-1, h/2
          xi(2, 1) = elem%eta(P_potP, 1) / 2**(elem%deg-1)  ! error
          xi(2, 2) = (elem%deg + 0) * (elem%deg + 1) * 2    ! DOF

          ! p, h
          xi(3, 1) = elem%eta(P_tot, 1)                     ! error
          xi(3, 2) = (elem%deg + 1) * (elem%deg + 2) / 2    ! DOF

          ! EOC for h-refinement
          val = log(xi(1, 1) / xi(2, 1)) /( xi(2, 2)**(1./3) - xi(1, 2)**(1./3) )

          ! EOC for p-refinement
          val1 = log(xi(1, 1) / xi(3, 1)) /( xi(3, 2)**(1./3) - xi(1, 2)**(1./3) )

          !val =  log(elem%eta(P_tot, 1) / elem%eta(P_potP, 1) ) &
          !     * log(1.*elem%deg /(elem%deg+1)) / log(1.*(elem%deg-1) /elem%deg  )

                  
          !write(*,'(a7, i5, 6es12.4)') 'regul:',elem%i, elem%eta(P_tot, 1), elem%eta(P_potP, 1), val, val1

          if(val1 >= val) then  
             elem%regT0 = elem%reg * 2.  ! p-refinement
             !   !write(*,'(a8, 6es12.4)') 'regul:', val, elem%xc(:)
          else
             elem%regT0 = elem%reg / 2.  ! h-refinement
             !   !write(*,'(a8, 6es12.4)') 'regul:', val, elem%xc(:), sqrt(dot_product(elem%xc(:), elem%xc(:)))
          endif

        else
          elem%regT0 = elem%reg * 2.  ! p-refinement
       endif

    case(10)
       if(elem%i == 1) print*,'! minimal DoF extrapolated to p+1 from potential'
       if(elem%deg > 1) then

          val =  log(elem%eta(P_pot, 1) / elem%eta(P_potP, 1) ) &
               * log(1.*elem%deg /(elem%deg+1)) / log(1.*(elem%deg-1) /elem%deg  )
          
          estPP = elem%eta(P_pot, 1) * exp(val)
          
          
          elem%regT0 =  (1. * (elem%deg+1) /(elem%deg +3)) **((elem%deg ) / 2.)
          
          elem%reg  = estPP/ elem%eta(P_pot, 1)
          
          !write(*,'(a10, 3i5, 4es12.4,a3,10es12.4)') &
          !     'elem_regu:',elem%i, elem%deg, int(val+0.5), &
          !     elem%eta(P_potP, 1), elem%eta(P_pot, 1) , estPP ,  elem%eta(P_tot, 1),  &
          !     '|||', elem%reg, elem%regT0, estPP/ elem%eta(P_pot, 1)

          
          ! VERSION ONLY TO p
          !!qq =  (elem%eta(P_pot, 1) / elem%eta(P_potP, 1)  )**(2./(elem%deg -1) )
          !!val = (elem%deg + 2.)/(elem%deg ) * qq
          
          elem%reg  = elem%eta(P_pot, 1) / elem%eta(P_potP, 1)
          elem%regT0 = 0.95 * (1. * elem%deg /(elem%deg +2)) **((elem%deg -1) / 2.)

          
        else
          elem%regT0 = elem%reg * 2.  ! p-refinement
       endif

    case(11)
       if(elem%i == 1) print*,'! hp-derefinement indication'

       call Regularity_hp_derefinement(elem, err_h_deref, err_p_deref)

       !slope_h = (err_h_deref - elem%eta(P_tot, 1) ) / (1.5*(elem%deg +1)*(elem%deg + 2) )
       !slope_p = (err_p_deref - elem%eta(P_tot, 1) ) / ( elem%deg +1 )


       rn = 3
       allocate(xi(1:rn, 1:2) )
       
       ! p-1, h
       xi(1, 1) = err_p_deref                            ! error
       xi(1, 2) = (elem%deg + 0) * (elem%deg + 1) / 2    ! DOF
       
       ! p, 2 * h 
       xi(2, 1) = err_h_deref                            ! error
       xi(2, 2) = (elem%deg + 1) * (elem%deg + 2) / 8    ! DOF

       ! p, h
       xi(3, 1) = elem%eta(P_tot, 1)                     ! error
       xi(3, 2) = (elem%deg + 1) * (elem%deg + 2) / 2    ! DOF

       ! EOC for h-refinement
       slope_h = log(xi(2, 1) / xi(3, 1)) /( xi(3, 2)**(1./3) - xi(2, 2)**(1./3) )

       ! EOC for p-refinement
       slope_p = log(xi(1, 1) / xi(3, 1)) /( xi(3, 2)**(1./3) - xi(1, 2)**(1./3) )

       !write(*,'(a7, i5, 6es12.4)') 'regul:',elem%i, elem%eta(P_tot, 1), err_h_deref, err_p_deref,slope_h, slope_p

       if(slope_p >= slope_h) then  
          elem%regT0 = elem%reg * 2.  ! p-refinement
          !   !write(*,'(a8, 6es12.4)') 'regul:', val, elem%xc(:)
       else
          elem%regT0 = elem%reg / 2.  ! h-refinement
          !   !write(*,'(a8, 6es12.4)') 'regul:', val, elem%xc(:), sqrt(dot_product(elem%xc(:), elem%xc(:)))
       endif
       
       deallocate(xi)

     case  default
       print*,'EEE', state%space%adapt%type_regularity
       stop "trouble in  Set_Elem_regularity"

    end select

    elem%regT2 = elem%diam**(-0.5)    ! elem%reg >  elem%regT2  ==> p-DErefinement

       ! 1st variant
       !if(elem%deg >= 2) elem%reg  = elem%eta(P_F_p1, 1) / elem%eta(P_F_p2, 1) !/  elem%diam**0.5
       !if(elem%deg >= 3) elem%reg1 = elem%eta(P_F_p2, 1) / elem%eta(P_F_p3, 1) !/  elem%diam**0.5

       ! Babuska variant
       !if(elem%deg >= 2) elem%reg1 = elem%eta(P_F_p1, 1) / elem%eta(P_F_p2, 1) /  elem%diam**0.5
       !if(elem%deg >= 3) elem%reg2 = elem%eta(P_F_p2, 1) / elem%eta(P_F_p3, 1) /  elem%diam**0.5

       ! the following already NECESSARY
       !elem%regT0 = 0.75 * elem%diam**0.5       ! elem%reg <  elem%regT0  ==> p-refinement
       !elem%regT0 = 0.75 * elem%diam**0.5       ! elem%reg <  elem%regT0  ==> p-refinement SISC_1

       
       !elem%regT0 = 0.6       ! elem%reg <  elem%regT0  ==> p-refinement  BABUSKA

       !elem%regT0 = 0.75 * elem%diam**0.5       ! elem%reg <  elem%regT0  ==> p-refinement  SISC_2 ?
       !elem%regT0 = 1.00 * elem%diam**0.5       ! elem%reg <  elem%regT0  ==> p-refinement  SISC_2 ?

       !!!elem%regT0 = elem%deg /(elem%deg +1.) * elem%diam**0.5  ! CAN WORKS nicely, VERIFY
       ! test the following
       !elem%regT0 = (elem%deg -1. )/ elem%deg * elem%diam**0.5  ! CAN WORKS nicely, VERIFY
       !if(elem%deg == 1) elem%regT0 = elem%reg * 2  ! p -refinement 

       !elem%regT0 = 1.0 * (elem%area / elem%diam)**0.5       ! elem%reg <  elem%regT0  ==> p-refinement  SISC_2 ?
       !elem%regT0 = 1.25 *sqrt(4./sqrt(3.)*elem%area) ! elem%reg <  elem%regT0  ==> p-refinement  SISC_2 ?

       !if ( state%space%adapt%adapt_type == 'HG') &  ! elem%reg <  elem%regT0  ==> p-ref
       !     elem%regT0 = 1.5 * elem%diam**0.5    ! starting with p=2, En projection
       !     !elem%regT0 = 1.25 * elem%diam**0.5    ! starting with p=1, L^2-proj
       !     !elem%regT0 = 1.5 * elem%diam**0.5    ! starting with p=2, L^2-proj


       !!!elem%regT0 = 0.95 * elem%diam**0.5 ! case J      ! elem%reg <  elem%regT0  ==> p-refinement
       !elem%regT0 = 1.5*elem%diam**0.5       ! elem%reg <  elem%regT0  ==> p-refinement
       !elem%regT1 = 1.                   ! elem%reg <  elem%regT1  ==> h->p substitution

       !elem%regT0 = elem%diam**1.0                ! elem%reg <  elem%regT0  ==> p-refinement
       !elem%regT1 = elem%diam**(0.5)    ! elem%reg <  elem%regT1  ==> h->p substitution
       !elem%regT2 = 1.  !0.9  !elem%diam**(0.0)    ! elem%reg >  elem%regT2  ==> p-DErefinement









  end subroutine Set_Elem_regularity


  subroutine Regularity_Smoothing( )
    class(element), pointer :: elem, elem1
    real, dimension(:,:), allocatable :: regul !list of corresponding elements
    integer :: ismoothing, is, i, j, l1, l2
    real :: normF, normS, weightE, weightR

    ! smoothing of the regul parameter
    allocate(regul(1:grid%nelem, 1:12) )

    ismoothing = 1 !1 ! number of smoothing cycles
    weightE = 0.0  !0.05  ! weight for error estim of the neighbouring elements
    !weightR = 0.25 !0.75  ! weight for regularity of the neighbouring elements
    weightR = 0.00 !0.75  ! weight for regularity of the neighbouring elements

    if ( state%space%adapt%adapt_type == 'HG' .or. state%space%adapt%adapt_type == 'RG') then
       weightE = 0.0 
       weightR = 0.0 
    endif

    do is = 1, ismoothing
       regul(1:grid%nelem, 1) = grid%elem(1:grid%nelem)%reg
       regul(1:grid%nelem, 2) = grid%elem(1:grid%nelem)%reg1
       regul(1:grid%nelem, 3) = grid%elem(1:grid%nelem)%reg2

       regul(1:grid%nelem, 6) = grid%elem(1:grid%nelem)%estim_loc

       regul(1:grid%nelem, 11) = 1.
       regul(1:grid%nelem, 12) = 1.

       do i = 1, grid%nelem
          elem => grid%elem(i)
          do j=1,elem%flen
             if(elem%face(neigh,j) > 0) then
                elem1 => grid%elem(elem%face(neigh,j))

                regul(i,1) = regul(i,1) + weightR * elem1%reg
                regul(i,2) = regul(i,2) + weightR * elem1%reg1
                regul(i,3) = regul(i,3) + weightR * elem1%reg2
                regul(i,11) = regul(i,11) + weightR

                regul(i,6) = regul(i,6) + weightE * elem1%eta(resST, 1)
                regul(i,12) = regul(i,12) + weightE
             endif
          enddo
       enddo

       do i = 1, grid%nelem
          elem => grid%elem(i)

          !write(*,'(a5,i5, 3es12.4)') 'DESW',i, elem%reg, elem%reg1
          elem%reg     = regul(i,1) / regul(i,11)
          !elem%reg = max ( elem%reg , regul(i,1) / regul(i,11))

          elem%reg1    = regul(i,2) / regul(i,11)
          elem%reg2    = regul(i,3) / regul(i,11)

          elem%eta(resST, 1) = regul(i,6) / regul(i,12)
          !write(*,'(a5,i5, 3es12.4)') 'DESW',i, elem%reg, elem%reg1, elem%regT0
          !print*

       enddo
    enddo ! do is

    deallocate(regul)

    !print*,'Stopped in Neumann'
    !stop
  end subroutine Regularity_Smoothing


  !> indication of the hp-adaptation based on the formalhp-derefinement
  !> we evaluate the differences::  
  !> \f$ \| w_h - \Pi^{p-1}_h \w_h\|\f$  and \f$ \| w_h - \Pi^{p}_{2h} \w_h\|\f$ 
  !> the temporary arrays contains the following values
  !> elem%wSS(1, 0, :) = actual solution
  !> elem%wSS(1, 1, :) = p-1 projection
  !> elem%wSS(1, 2, :) = p-2 projection
  subroutine Regularity_hp_derefinement(elem, err_h_deref, err_p_deref)
    class(element), intent(inout) :: elem
    real, intent(inout) :: err_h_deref, err_p_deref
    class(element), pointer :: elem1
    type(volume_rule), pointer :: V_rule
    real, dimension(:,:), pointer :: phi0 ! the test functions
    real, dimension(:,:,:), allocatable :: Dphi0 ! the test functions
    real, dimension(:,:), allocatable :: phi ! the test functions
    real, dimension(:,:), allocatable :: phiW ! the test functions multiplied by weights
    real, dimension(:,:,:), allocatable :: Dphi ! Der of test functions
    real, dimension(:,:,:), allocatable :: DphiW ! Der of test functions multiplied by weights
    real, dimension(:,:), allocatable :: xi ! barycentric coordinates
    real, dimension(:,:), allocatable :: Fx ! real physical coordinates
    real, dimension(:,:), allocatable :: Mass ! mass matrix
    real, dimension(:,:), allocatable :: Stiff ! mass matrix
    real, dimension(:,:), allocatable :: rhs ! right-hand-side
    real, dimension(:,:), allocatable :: rhsSff ! right-hand-side
    real, dimension(:), allocatable :: weights, wi
    real, dimension(:,:), allocatable :: Dwi
    real, dimension(:,:), allocatable :: wR  ! reconstructed solution
    real, dimension(:,:), allocatable :: DwiR  ! reconstructed solution
    integer :: ndimL, deg, dof, Qdof, Qdof1, dof1, N_elem
    integer :: i, i1, j, k, l, l1, l2, n, n1
    integer:: p_min, p_max, p_diff, itype
    real :: weight, L2weight, err, errK, area

    ndimL = ndim
    deg = elem%deg
    dof = elem%dof

    V_rule => state%space%V_rule(elem%Qnum)
    Qdof = V_rule%Qdof

    ! Mass(i,j) = \int_{K and its neighbours} \phi_i \phi_j dx
    ! rhs(i) = \int_{K and its neighbours} \phi_i w_h dx
    allocate(Mass(1:dof, 1:dof), Stiff(1:dof, 1:dof))
    allocate(rhs(1:dof, 1:ndimL), rhsSff(1:dof, 1:ndimL) )

    ! reconstructed solution in basis coefficients
    allocate(wR(1:ndimL, 1:dof) )

    err = 0.   ! difference bewteen w_h ands its interpolant


    do itype = 1, 2 ! itype == 1 ==> we compute the projection; itype == 2 ==> we compute the error

       if(itype == 1) then
          Mass(:,:) = 0.
          Stiff(:,:) = 0.
          rhs(:,:) = 0.
          rhsSff(:,:) = 0.
          
          !integerals over elem
          Mass(1:dof,1:dof)  = elem%Mass%Mb(1:dof,1:dof)
          Stiff(1:dof,1:dof) = elem%Stiff%Mb(1:dof,1:dof)
          
          area = elem%area
       endif

       !integerals over elem
       phi0 => V_rule%phi(1:dof, 1:Qdof)

       allocate( Dphi0(1:dof, 1:nbDim, 1:Qdof) )
       call Eval_Dphi(elem, dof, Dphi0(1:dof, 1:nbDim, 1:Qdof) )


       ! integration of (ww, \phi)_K, and (D\w, D\phi)_K
       allocate( wi(1:Qdof) , Dwi(1:Qdof, 1:nbDim) )
       allocate( DwiR(1:Qdof, 1:nbDim) )    ! reconstructed solution in integ nodes

       do k=1, ndimL
          ! product (wi, phi0)
          wi(1:Qdof) = matmul(elem%wSS(1, 0, (k-1)*dof+1 : k*dof), phi0(1:dof, 1:Qdof) )
          
          if(itype == 1) &
               call IntegrateVectorB(elem, dof, wi(1:Qdof), rhs(1:dof, k) )
          
          ! product (Dwi, Dphi0)
          do n=1,nbDim
             Dwi(1:Qdof, n) = matmul(elem%wSS(1, 0, (k-1)*dof+1 : k*dof), Dphi0(1:dof, n, 1:Qdof))
          enddo
          
          if(itype == 1) &
               call IntegrateVectorD(elem, dof, Dwi(1:Qdof, 1:nbDim), rhsSff(1:dof, k))

          if(itype == 2) then
             allocate (weights(1:Qdof)  )
             call Eval_V_Weights(elem, weights(1:Qdof))

             do n=1,nbDim
                DwiR(1:Qdof, n) = matmul(wR(k, 1:dof), Dphi0(1:dof, n, 1:Qdof))
             enddo
             
             err = err + dot_product(weights(1:Qdof) , &
                  ( Dwi(1:Qdof, 1) - DwiR(1:Qdof, 1) ) **2  + ( Dwi(1:Qdof, 2) - DwiR(1:Qdof, 2) ) **2 )

             errK = err

             !write(*,'(a7, 30es12.4)') 'Dwi, 1',Dwi(1:Qdof,1) 
             !write(*,'(a7, 30es12.4)') 'DwiR, 1',DwiR(1:Qdof,1) 
             !write(*,'(a7, 30es12.4)') 'Dwi, 2',Dwi(1:Qdof,2) 
             !write(*,'(a7, 30es12.4)') 'DwiR, 2',DwiR(1:Qdof,2) 

             !print*,'###', elem%i, 0, sqrt(err)


             ! p-derefinement, projection into p-1
             do n=1,nbDim
                DwiR(1:Qdof, n) = matmul(elem%wSS(1, 1, (k-1)*dof+1 : k*dof), Dphi0(1:dof, n, 1:Qdof))
             enddo
             err_p_deref =  dot_product(weights(1:Qdof) , &
                  ( Dwi(1:Qdof, 1) - DwiR(1:Qdof, 1) ) **2  + ( Dwi(1:Qdof, 2) - DwiR(1:Qdof, 2) ) **2 )

             err_p_deref = sqrt(err_p_deref)

             deallocate(weights)
          endif
       enddo
  
          

       deallocate(wi, Dwi, Dphi0, DwiR )

       ! !!goto 111

       ! adding of integrals over neighbouring elements
       ! phi0  basis functions from elem1 in elem1's integ nodes
       ! phi   basis functions from elem  in elem1's integ nodes, i.e., outside of elem
       
       
       ! we go over neighbours of elem
       N_elem = elem%flen
       do j=1, N_elem
          i1 = elem%face(neigh, j)
          
          ! N_elem = elem%isupp
          ! do j=1, N_elem
          !    i1 = elem%supp(j, 1)
          !    if(i1 > 0) then
          !       p_min = min(p_min, grid%elem(i1)%deg)
          !       p_max = max(p_max, grid%elem(i1)%deg)
          !    endif
          ! enddo
          ! p_diff = p_max - p_min
          
          
          ! N_elem = elem%isupp   ! elements sharing a vertex
          
          ! do j=1, N_elem
          !i1 = elem%supp(j, 1)
          
          if(i1 > 0) then
             weight = 1.
             !if(elem%supp(j, 2) == 1) weight = 1.
             !if(elem%supp(j, 2) == 2) weight = 0.02  !0.5  !0.5   !0.05  ! 0.005
             !if(elem%supp(j, 2) == 2) weight = 0.02 * (p_diff -1)  ! AMC 2015, HO_reconstruction
             
             if(elem%deg == 0) weight = 1.

             if(weight > 0.) then
                !if(elem%supp(j,2)==2 .and. p_diff > 0) &
                !     print*,'Weight-',elem%i,j,i1, elem%supp(j,2) ,p_diff , weight
                
                elem1 => grid%elem(i1)
                dof1  = elem1%dof
                
                !V_rule => state%space%V_rule(elem%Qnum)  ! quadrature on the neighbour
                V_rule => state%space%V_rule(elem%Qnum)   ! currect quandature
                Qdof1 = V_rule%Qdof
                
                
                ! Fx integ nodes on elem1 - real physical cordinates
                ! xi barycentric coordinates of Fx (integ nodes on elem1) with respect elem
                allocate(Fx(1:Qdof1, 1:2), xi(1:Qdof1, 1:2) )
                call ComputeF(elem1, Qdof1, V_rule%lambda(1:Qdof1, 1:2), Fx(1:Qdof1, 1:2) )
                call BarycCoord(elem, Qdof1, Fx(1:Qdof1, 1:2),  xi(1:Qdof1, 1:2) )
                
                allocate(phi(1:dof, 1:Qdof1) ) ! value of the test function in integ nodes
                allocate(Dphi(1:dof, 1:nbDim, 1:Qdof1) ) ! value of Deriv test functions
                
                call Eval_phi_Qnode(elem, dof, Qdof1, xi(1:Qdof1, 1:nbDim), &
                     phi(1:dof, 1:Qdof1), Dphi(1:dof, 1:nbDim, 1:Qdof1) )
                

                allocate(wi(1:Qdof1),  weights(1:Qdof1) )
                allocate(Dwi(1:Qdof1, 1:nbDim)  )
                allocate(DwiR(1:Qdof1, 1:nbDim) )    ! reconstructed solution in integ nodes

                call Eval_V_Weights_plus(elem1, V_rule, weights(1:Qdof1))
                
                allocate(phiW(1:dof, 1:Qdof1) )           ! phi multiplied by weights
                allocate(DphiW(1:dof, 1:nbDim, 1:Qdof1) ) ! Dphi multiplied by weights
                
                do l=1, dof
                   phiW(l, 1:Qdof1) = phi(l, 1:Qdof1) * weights(1:Qdof1)
                   
                   do n=1,nbDim
                      DphiW(l, n, 1:Qdof1) = Dphi(l, n, 1:Qdof1) * weights(1:Qdof1)
                   enddo
                enddo
                

                ! adding of Mass and Stiff  matrices
                if(itype == 1) then
                   do l1 = 1, dof
                      do l2 = 1, dof
                         Mass(l1, l2) = Mass(l1, l2) &
                              + weight * dot_product(phiW(l1, 1:Qdof1), phi(l2, 1:Qdof1) )

                         do n=1,nbDim
                            Stiff(l1, l2) = Stiff(l1, l2) &
                                 + weight * dot_product(DphiW(l1, n, 1:Qdof1), Dphi(l2, n, 1:Qdof1) )
                         enddo
                      enddo
                   enddo

                   !area of the patch
                   area = area + elem1%area

                endif

                ! adding of rhs
                phi0 => V_rule%phi(1:dof1, 1:Qdof1)

                allocate( Dphi0(1:dof1, 1:nbDim, 1:Qdof1) )
                call Eval_Dphi_plus(elem1, V_rule, dof1, Dphi0(1:dof1, 1:nbDim, 1:Qdof1) )


                do k=1, ndimL
                   ! wi = values of w at integ nodes of elem1
                   wi(1:Qdof1) = matmul(elem1%wSS(1, 0, (k-1)*dof1+1 : k*dof1),  phi0(1:dof1, 1:Qdof1) )
                   !write(*,'(a8, 2i5, 20es12.4)') 'de3:',elem%i,k, wi(1), (elem1%xc(1)/40)**2, elem1%xc
                   ! product (wi, phi0)

                    if(itype == 1) then
                       do l=1, dof
                          rhs(l, k) = rhs(l, k)  &
                               + weight * dot_product(wi(1:Qdof1), phiW(l, 1:Qdof1) )
                       enddo
                    endif
                   
                   ! Dwi in integ nodes
                   do n=1,nbDim
                      Dwi(1:Qdof1, n) = matmul(elem1%wSS(1, 0, (k-1)*dof1+1 : k*dof1), &
                           Dphi0(1:dof1, n, 1:Qdof1))
                   enddo


                   if(itype == 1) then
                      ! product (Dwi, Dphi0)
                      do l=1, dof
                         do n=1,nbDim
                            rhsSff(l, k) = rhsSff(l, k)  &
                                 + weight * dot_product(Dwi(1:Qdof1, n), DphiW(l, n, 1:Qdof1) )
                         enddo
                      enddo

                   else !!!   (itype == 2) 
                      do n=1,nbDim
                         ! bellow Dphi is correct, in Dwi is Dphi0 !!!
                         DwiR(1:Qdof1, n) = matmul(wR(k, 1:dof), Dphi(1:dof, n, 1:Qdof1)) 
                      enddo
             
                      err = err + dot_product(weights(1:Qdof1) , &
                           ( Dwi(1:Qdof1, 1) - DwiR(1:Qdof1, 1) ) **2  &
                           + ( Dwi(1:Qdof1, 2) - DwiR(1:Qdof1, 2) ) **2 )

                      !print*,'###', elem%i, j, sqrt(err)

                   endif

                enddo

                deallocate(Fx, xi, phi, phiW, wi, weights)
                deallocate(Dphi, DphiW, Dwi, Dphi0, DwiR)

             endif   ! if (weight > 0.) then

          endif  !(i1 > 0)

       enddo ! j=1, elem%flen ( N_elem)

!!!111 continue
       
       if(itype == 1) then
          ! the H1-norm (comment for the L^2-norm)
          L2weight = 1.0
          Mass(1:dof, 1:dof) = Mass(1:dof, 1:dof) +  L2weight * Stiff(1:dof, 1:dof)
          rhs(1:dof, 1:ndimL) = rhs(1:dof, 1:ndimL) +  L2weight * rhsSff(1:dof, 1:ndimL)
          
          call SolveLocalMatrixProblem(dof, Mass(1:dof, 1:dof), ndimL, rhs(1:dof, 1:ndimL) )
          wR(1:ndimL, 1:dof) = transpose( rhs(1:dof, 1:ndimL))
          

          !call PlotElemFunction3D(10, elem, elem%dof, wR(1, 1:elem%dof))    
       endif

       !print*,'###', elem%i, j, sqrt(err)

    end do  ! do itype = 1, 2
    
    !err_h_deref = sqrt(err * elem%area/ area)
    err_h_deref = sqrt(errK)

    !write(*,'(a10, i5, 30es12.4)') 'Der: h, p:',elem%i, sqrt(errK), &
    !     err_h_deref, err_p_deref, elem%eta(P_tot, 1)

    !stop "d383d38yd38kdjw"


    deallocate( Mass, Stiff, rhs, rhsSff, wR )




  end subroutine Regularity_hp_derefinement




  !> indication of the hp-adaptation based on the formalhp-derefinement
  !> we evaluate the differences::  
  !> \f$ \| w_h - \Pi^{p-1}_h \w_h\|\f$, \f$ \| w_h - \Pi^{p}_{2h} \w_h\|\f$ 
  !> and \f$ \| w_h - \Pi^{p-1}_{2h} \w_h\|\f$
  !> \f$ 2 h \f$ corresponds to the patches from elements sharing only faces or at least nodes
  !> the temporary arrays contains the following values
  !> elem%wSS(1, 0, :) = actual solution (input)
  !> elem%wSS(1, 1, :) = p-1 projection  (input)
  subroutine Regularity_full_hp_derefinement(elem, slope_h, slope_p) 
    class(element), intent(inout) :: elem
    real, intent(inout) :: slope_h, slope_p
    class(element), pointer :: elem1
    type(volume_rule), pointer :: V_rule
    real, dimension(:,:), pointer :: phi0 ! the test functions
    real, dimension(:,:,:), allocatable :: Dphi0 ! the test functions
    real, dimension(:,:), allocatable :: phi ! the test functions
    real, dimension(:,:), allocatable :: phiW ! the test functions multiplied by weights
    real, dimension(:,:,:), allocatable :: Dphi ! Der of test functions
    real, dimension(:,:,:), allocatable :: DphiW ! Der of test functions multiplied by weights
    real, dimension(:,:), allocatable :: xi ! barycentric coordinates
    real, dimension(:,:), allocatable :: Fx ! real physical coordinates
    real, dimension(:,:), allocatable :: MM ! final matrix matrix
    real, dimension(:,:), allocatable :: Mass ! mass matrix
    real, dimension(:,:), allocatable :: Stiff ! stiff matrix
    real, dimension(:,:), allocatable :: bb! right-hand-side
    real, dimension(:,:), allocatable :: rhs ! right-hand-side
    real, dimension(:,:), allocatable :: rhsSff ! right-hand-side
    real, dimension(:), allocatable :: weights
    real, dimension(:,:), allocatable :: Dwi
    real, dimension(:,:, :), allocatable :: wR  ! reconstructed solution
    real, dimension(:,:), allocatable :: DwiR  ! reconstructed solution
    integer :: ndimL, deg, dof, dofM,  Qdof, Qdof1, dof1, N_elem
    integer :: i, i1, j, k, l, l1, l2, n, n1
    integer:: p_min, p_max, p_diff, itype, rn
    !real :: err_h_deref, err_p_deref, err_hp_deref
    real :: eK_h, eK_p, eK_hp
    real :: weight, L2weight, area
    real :: err_2h_p, err_2h_p1, err_h_p1, errK, errKp
    real :: dof_2h_p, dof_2h_p1, dof_h_p1
    logical :: iprint

    iprint = .false.
    if( elem%i == -1 .or. elem%i == -38) then
    !if( dot_product( grid%x(elem%face(idx, 1), 1:2 ),  grid%x(elem%face(idx, 1), 1:2 ) )<1E-8 .or. &
    !     dot_product( grid%x(elem%face(idx, 2), 1:2 ),  grid%x(elem%face(idx, 2), 1:2 ) )< 1E-8 .or. &
    !     dot_product( grid%x(elem%face(idx, 3), 1:2 ),  grid%x(elem%face(idx, 3), 1:2 ) )< 1E-8 ) then
       iprint = .true.
    endif


    ndimL = ndim
    deg = elem%deg
    dof = elem%dof
    dofM = dof - deg - 1  ! projection of degree deg-1

    V_rule => state%space%V_rule(elem%Qnum)
    Qdof = V_rule%Qdof

    ! Mass(i,j) = \int_{K and its neighbours} \phi_i \phi_j dx
    ! rhs(i) = \int_{K and its neighbours} \phi_i w_h dx
    allocate(Mass(1:dof, 1:dof), Stiff(1:dof, 1:dof), MM(1:dof, 1:dof))
    allocate(rhs(1:dof, 1:ndimL), rhsSff(1:dof, 1:ndimL), bb(1:dof, 1:ndimL) )

    ! reconstructed solution in basis coefficients
    allocate(wR(1:2, 1:ndimL, 1:dof) )   ! wR(1, :, :) = reconstruction of degree deg
    !                                    ! wR(2, :, :) = reconstruction of degree deg-1
    wR = 0.


    err_2h_p  = 0.   ! difference bewteen w_h ands its interpolant of degree p   on the patch
    err_2h_p1 = 0.   ! difference bewteen w_h ands its interpolant of degree p-1 on the patch
    err_h_p1  = 0.   ! difference bewteen w_h ands its interpolant of degree p-1 on each element

    dof_2h_p  = 0.
    dof_2h_p1 = 0.
    dof_h_p1  = 0.

    do itype = 1, 2 ! itype == 1 ==> we compute the projection; itype == 2 ==> we compute the error

       if(itype == 1) then
          Mass(:,:) = 0.
          Stiff(:,:) = 0.
          rhs(:,:) = 0.
          rhsSff(:,:) = 0.
          
          !integerals over elem
          Mass(1:dof,1:dof)  = elem%Mass%Mb(1:dof,1:dof)
          Stiff(1:dof,1:dof) = elem%Stiff%Mb(1:dof,1:dof)
          
          area = elem%area
       endif

       !integerals over elem
       phi0 => V_rule%phi(1:dof, 1:Qdof)

       allocate( Dphi0(1:dof, 1:nbDim, 1:Qdof) )
       call Eval_Dphi(elem, dof, Dphi0(1:dof, 1:nbDim, 1:Qdof) )


       ! integration of (ww, \phi)_K, and (D\w, D\phi)_K
       allocate( Dwi(1:Qdof, 0:nbDim) )   ! Dwi(:, 0) = wi, 
       allocate( DwiR(1:Qdof, 1:nbDim) )    ! reconstructed solution in integ nodes

       do k=1, ndimL
          ! product (wi, phi0)
          Dwi(1:Qdof, 0) = matmul(elem%wSS(1, 0, (k-1)*dof+1 : k*dof), phi0(1:dof, 1:Qdof) )
          
          ! product (Dwi, Dphi0)
          do n=1,nbDim
             Dwi(1:Qdof, n) = matmul(elem%wSS(1, 0, (k-1)*dof+1 : k*dof), Dphi0(1:dof, n, 1:Qdof))
          enddo

          if(itype == 1) then
             call IntegrateVectorB(elem, dof, Dwi(1:Qdof, 0), rhs(1:dof, k) )
          
             call IntegrateVectorD(elem, dof, Dwi(1:Qdof, 1:nbDim), rhsSff(1:dof, k))

          endif


          if(itype == 2) then
             allocate (weights(1:Qdof)  )
             call Eval_V_Weights(elem, weights(1:Qdof))

             ! h-derefinement
             do n=1,nbDim
                DwiR(1:Qdof, n) = matmul(wR(1, k, 1:dof), Dphi0(1:dof, n, 1:Qdof))
             enddo
             
             err_2h_p = err_2h_p + dot_product(weights(1:Qdof) , &
                  ( Dwi(1:Qdof, 1) - DwiR(1:Qdof, 1) ) **2  + ( Dwi(1:Qdof, 2) - DwiR(1:Qdof, 2) ) **2 )

             eK_h =  dot_product(weights(1:Qdof) , &
                  ( Dwi(1:Qdof, 1) - DwiR(1:Qdof, 1) ) **2  + ( Dwi(1:Qdof, 2) - DwiR(1:Qdof, 2) ) **2 )

             dof_2h_p = dof_2h_p + dof


             ! hp-derefinement
             do n=1,nbDim
                DwiR(1:Qdof, n) = matmul(wR(2, k, 1:dof), Dphi0(1:dof, n, 1:Qdof))
             enddo
             
             err_2h_p1 = err_2h_p1 + dot_product(weights(1:Qdof) , &
                  ( Dwi(1:Qdof, 1) - DwiR(1:Qdof, 1) ) **2  + ( Dwi(1:Qdof, 2) - DwiR(1:Qdof, 2) ) **2 )

             eK_hp = dot_product(weights(1:Qdof) , &
                  ( Dwi(1:Qdof, 1) - DwiR(1:Qdof, 1) ) **2  + ( Dwi(1:Qdof, 2) - DwiR(1:Qdof, 2) ) **2 )

             dof_2h_p1 = dof_2h_p1 + dof - deg - 1


             ! p-derefinement, projection into p-1
             do n=1,nbDim
                DwiR(1:Qdof, n) = matmul(elem%wSS(1, 1, (k-1)*dof+1 : k*dof) &
                     - elem%wSS(1, 0, (k-1)*dof+1 : k*dof), Dphi0(1:dof, n, 1:Qdof))
             enddo
             
             err_h_p1 = err_h_p1 + dot_product(weights(1:Qdof) , &
                  ( DwiR(1:Qdof, 1) ) **2  + (DwiR(1:Qdof, 2) ) **2 )

             eK_p = dot_product(weights(1:Qdof) , &
                  ( DwiR(1:Qdof, 1) ) **2  + (DwiR(1:Qdof, 2) ) **2 )

             dof_h_p1 = dof_h_p1 + dof - deg - 1


             if(iprint) then
                allocate(Fx(1:Qdof, 1:2) )
                call ComputeF(elem, Qdof, V_rule%lambda(1:Qdof, 1:2), Fx(1:Qdof, 1:2) )

                write(*,'(a4, i5, 30es12.4)') 'wSS0:',dof, elem%wSS(1, 0,1:dof)
                write(*,'(a4, i5, 30es12.4)') 'wSS1:',dof, elem%wSS(1, 1,1:dof)
                
                write(*,'(a4, i5, 30es12.4)') 'wR1:',dof, wR(1, 1, 1:dof)
                write(*,'(a4, i5, 30es12.4)') 'wR2:',dof, wR(2, 1, 1:dof)

                do n=1,Qdof
                   write(71, *) Fx(n, 1:2), &
                        dot_product(elem%wSS(1, 0, 1:dof), phi0(1:dof,  n )) &
                        - dot_product(elem%wSS(1, 1, 1:dof), phi0(1:dof,  n ))
                   
                   write(81, *) Fx(n, 1:2), &
                        dot_product(elem%wSS(1, 0, 1:dof), phi0(1:dof,  n )) , &
                        sqrt( dot_product(elem%wSS(1, 0, 1:dof), Dphi0(1:dof, 1,  n ))**2 &
                        + dot_product(elem%wSS(1, 0, 1:dof), Dphi0(1:dof, 2,  n ))**2) , &
                        dot_product(elem%wSS(1, 0, 1:dof), Dphi0(1:dof, 1,  n )), &
                        dot_product(elem%wSS(1, 0, 1:dof), Dphi0(1:dof, 2,  n ))
                   
                   write(82, *) Fx(n, 1:2), &
                        dot_product(elem%wSS(1, 1, 1:dof), phi0(1:dof,  n )) , &
                        sqrt( dot_product(elem%wSS(1, 1, 1:dof), Dphi0(1:dof, 1,  n ))**2 &
                        + dot_product(elem%wSS(1, 1, 1:dof), Dphi0(1:dof, 2,  n ))**2) , &
                        dot_product(elem%wSS(1, 1, 1:dof), Dphi0(1:dof, 1,  n )), &
                        dot_product(elem%wSS(1, 1, 1:dof), Dphi0(1:dof, 2,  n ))
                   
                   write(83, *) Fx(n, 1:2), dot_product(wR(1, 1, 1:dof), phi0(1:dof,  n )) , &
                        sqrt( dot_product(wR(1, 1, 1:dof), Dphi0(1:dof, 1,  n ))**2 &
                        + dot_product(wR(1, 1, 1:dof), Dphi0(1:dof, 2,  n ))**2) , &
                        dot_product(wR(1, 1, 1:dof), Dphi0(1:dof, 1,  n )), &
                        dot_product(wR(1, 1, 1:dof), Dphi0(1:dof, 2,  n ))
                   write(84, *) Fx(n, 1:2), dot_product(wR(2, 1, 1:dof), phi0(1:dof,  n )) , &
                        sqrt( dot_product(wR(2, 1, 1:dof), Dphi0(1:dof, 1,  n ))**2 &
                        + dot_product(wR(2, 1, 1:dof), Dphi0(1:dof, 2,  n ))**2) , &
                        dot_product(wR(2, 1, 1:dof), Dphi0(1:dof, 1,  n )), &
                        dot_product(wR(2, 1, 1:dof), Dphi0(1:dof, 2,  n ))
                enddo
                deallocate(Fx)
             endif

             if(iprint ) &
                  write(*,'(a6,2i5,es12.4,a1, 3es12.4,a2,3es12.4)') &
                  'hp;',elem%i, elem%deg, sqrt(dot_product(elem%xc, elem%xc)), '|', &
                  sqrt(eK_h), sqrt(eK_p), sqrt(eK_hp),'||', &
                  sqrt(eK_hp)/ sqrt(eK_h), sqrt(eK_hp)/ sqrt(eK_p)

             deallocate(weights)
          endif
       enddo
  
       deallocate(Dwi, Dphi0, DwiR )


 
       ! adding of integrals over neighbouring elements
       ! phi0  basis functions from elem1 in elem1's integ nodes
       ! phi   basis functions from elem  in elem1's integ nodes, i.e., outside of elem
       
       
       ! we go over neighbours of elem
       !N_elem = elem%flen
       !do j=1, N_elem
       !   i1 = elem%face(neigh, j)
          
          ! N_elem = elem%isupp
          ! do j=1, N_elem
          !    i1 = elem%supp(j, 1)
          !    if(i1 > 0) then
          !       p_min = min(p_min, grid%elem(i1)%deg)
          !       p_max = max(p_max, grid%elem(i1)%deg)
          !    endif
          ! enddo
          ! p_diff = p_max - p_min
          
          
       N_elem = elem%isupp   ! elements sharing a vertex
       do j=1, N_elem
          i1 = elem%supp(j, 1)
          
          if(i1 > 0) then
             weight = 1.
             !if(elem%supp(j, 2) == 1) weight = 1.
             !if(elem%supp(j, 2) == 2) weight = 0.02  !0.5  !0.5   !0.05  ! 0.005
             !if(elem%supp(j, 2) == 2) weight = 0.02 * (p_diff -1)  ! AMC 2015, HO_reconstruction
             
             if(elem%deg == 0) weight = 1.

             if(weight > 0.) then
                !if(elem%supp(j,2)==2 .and. p_diff > 0) &
                !     print*,'Weight-',elem%i,j,i1, elem%supp(j,2) ,p_diff , weight
                
                elem1 => grid%elem(i1)
                dof1  = elem1%dof
                
                V_rule => state%space%V_rule(elem%Qnum)   ! correct quandature
                Qdof1 = V_rule%Qdof
                
                
                ! Fx integ nodes on elem1 - real physical cordinates
                ! xi barycentric coordinates of Fx (integ nodes on elem1) with respect elem
                allocate(Fx(1:Qdof1, 1:2), xi(1:Qdof1, 1:2) )
                call ComputeF(elem1, Qdof1, V_rule%lambda(1:Qdof1, 1:2), Fx(1:Qdof1, 1:2) )
                call BarycCoord(elem, Qdof1, Fx(1:Qdof1, 1:2),  xi(1:Qdof1, 1:2) )
                
                allocate(phi(1:dof, 1:Qdof1) ) ! value of the test function in integ nodes
                allocate(Dphi(1:dof, 1:nbDim, 1:Qdof1) ) ! value of Deriv test functions
                
                call Eval_phi_Qnode(elem, dof, Qdof1, xi(1:Qdof1, 1:nbDim), &
                     phi(1:dof, 1:Qdof1), Dphi(1:dof, 1:nbDim, 1:Qdof1) )
                

                allocate(weights(1:Qdof1) )
                allocate(Dwi( 1:Qdof1, 0:nbDim)  )
                allocate(DwiR(1:Qdof1, 1:nbDim) )    ! reconstructed solution in integ nodes

                call Eval_V_Weights_plus(elem1, V_rule, weights(1:Qdof1))
                
                allocate(phiW(1:dof, 1:Qdof1) )           ! phi multiplied by weights
                allocate(DphiW(1:dof, 1:nbDim, 1:Qdof1) ) ! Dphi multiplied by weights
                
                do l=1, dof
                   phiW(l, 1:Qdof1) = phi(l, 1:Qdof1) * weights(1:Qdof1)
                   
                   do n=1,nbDim
                      DphiW(l, n, 1:Qdof1) = Dphi(l, n, 1:Qdof1) * weights(1:Qdof1)
                   enddo
                enddo
                

                ! adding of Mass and Stiff  matrices
                if(itype == 1) then
                   do l1 = 1, dof
                      do l2 = 1, dof
                         Mass(l1, l2) = Mass(l1, l2) &
                              + weight * dot_product(phiW(l1, 1:Qdof1), phi(l2, 1:Qdof1) )

                         do n=1,nbDim
                            Stiff(l1, l2) = Stiff(l1, l2) &
                                 + weight * dot_product(DphiW(l1, n, 1:Qdof1), Dphi(l2, n, 1:Qdof1) )
                         enddo
                      enddo
                   enddo

                   !area of the patch
                   area = area + elem1%area

                endif

                phi0 => V_rule%phi(1:dof1, 1:Qdof1)

                allocate( Dphi0(1:dof1, 1:nbDim, 1:Qdof1) )
                call Eval_Dphi_plus(elem1, V_rule, dof1, Dphi0(1:dof1, 1:nbDim, 1:Qdof1) )


                do k=1, ndimL
                   ! wi = values of w at integ nodes of elem1
                   Dwi(1:Qdof1, 0) = matmul(elem1%wSS(1, 0, (k-1)*dof1+1: k*dof1),phi0(1:dof1, 1:Qdof1))
                   ! Dwi in integ nodes
                   do n=1,nbDim
                      Dwi(1:Qdof1, n) = matmul(elem1%wSS(1, 0, (k-1)*dof1+1 : k*dof1), &
                           Dphi0(1:dof1, n, 1:Qdof1))
                   enddo

                   ! adding of rhs
                   if(itype == 1) then
                      ! product (wi, phi0)
                      do l=1, dof
                         rhs(l, k) = rhs(l, k)  &
                              + weight * dot_product(Dwi(1:Qdof1, 0), phiW(l, 1:Qdof1) )

                         ! product (Dwi, Dphi0)
                         do n=1,nbDim
                            rhsSff(l, k) = rhsSff(l, k)  &
                                 + weight * dot_product(Dwi(1:Qdof1, n), DphiW(l, n, 1:Qdof1) )
                         enddo
                      enddo
                   endif

                   if(itype == -2) then

                      ! h-derefinement
                      do n=1,nbDim
                         ! bellow Dphi is correct, in Dwi is Dphi0 !!!
                         DwiR(1:Qdof1, n) = matmul(wR(1, k, 1:dof), Dphi(1:dof, n, 1:Qdof1)) 
                      enddo
             
                      err_2h_p = err_2h_p + dot_product(weights(1:Qdof1) , &
                           ( Dwi(1:Qdof1, 1) - DwiR(1:Qdof1, 1) ) **2  &
                           + ( Dwi(1:Qdof1, 2) - DwiR(1:Qdof1, 2) ) **2 )

                      eK_h = dot_product(weights(1:Qdof1) , &
                           ( Dwi(1:Qdof1, 1) - DwiR(1:Qdof1, 1) ) **2  &
                           + ( Dwi(1:Qdof1, 2) - DwiR(1:Qdof1, 2) ) **2 )

                      ! hp-derefinement
                      do n=1,nbDim
                         ! bellow Dphi is correct, in Dwi is Dphi0 !!!
                         DwiR(1:Qdof1, n) = matmul(wR(2, k, 1:dof), Dphi(1:dof, n, 1:Qdof1)) 
                      enddo
             
                      err_2h_p1 = err_2h_p1 + dot_product(weights(1:Qdof1) , &
                           ( Dwi(1:Qdof1, 1) - DwiR(1:Qdof1, 1) ) **2  &
                           + ( Dwi(1:Qdof1, 2) - DwiR(1:Qdof1, 2) ) **2 )

                      eK_hp = dot_product(weights(1:Qdof1) , &
                           ( Dwi(1:Qdof1, 1) - DwiR(1:Qdof1, 1) ) **2  &
                           + ( Dwi(1:Qdof1, 2) - DwiR(1:Qdof1, 2) ) **2 )

                      ! p-derefinement, projection into p-1
                      do n=1,nbDim
                         DwiR(1:Qdof1, n) = matmul(elem1%wSS(1, 1, (k-1)*dof1+1 : k*dof1) &
                              - elem1%wSS(1, 0, (k-1)*dof1+1 : k*dof1), Dphi0(1:dof1, n, 1:Qdof1))
                      enddo
             
                      err_h_p1 = err_h_p1 + dot_product(weights(1:Qdof1) , &
                           ( DwiR(1:Qdof1, 1) ) **2  + (DwiR(1:Qdof1, 2) ) **2 )

                      eK_p =  dot_product(weights(1:Qdof1) , &
                           ( DwiR(1:Qdof1, 1) ) **2  + (DwiR(1:Qdof1, 2) ) **2 )


                      if(iprint) then
                         if(j == 1) then
                            write(*,'(a4, i5, 30es12.4)') 'wSS0:',dof1, elem1%wSS(1, 0,1:dof1)
                            write(*,'(a4, i5, 30es12.4)') 'wSS1:',dof1, elem1%wSS(1, 1,1:dof1)
                            
                            write(*,'(a4, i5, 30es12.4)') 'wR1:',dof, wR(1, 1, 1:dof)
                            write(*,'(a4, i5, 30es12.4)') 'wR2:',dof, wR(2, 1, 1:dof)
                         endif

                         do n=1,Qdof1
                            write(71, *) Fx(n, 1:2), &
                                 dot_product(elem1%wSS(1, 0, 1:dof1), phi0(1:dof1,  n )) &
                                 - dot_product(elem1%wSS(1, 1, 1:dof1), phi0(1:dof1,  n ))

                            write(81, *) Fx(n, 1:2), &
                                 dot_product(elem1%wSS(1, 0, 1:dof1), phi0(1:dof1,  n )) , &
                                 sqrt( dot_product(elem1%wSS(1, 0, 1:dof1), Dphi0(1:dof1, 1,  n ))**2 &
                                 + dot_product(elem1%wSS(1, 0, 1:dof1), Dphi0(1:dof1, 2,  n ))**2) , &
                                 dot_product(elem1%wSS(1, 0, 1:dof1), Dphi0(1:dof1, 1,  n )), &
                                 dot_product(elem1%wSS(1, 0, 1:dof1), Dphi0(1:dof1, 2,  n ))

                            write(82, *) Fx(n, 1:2), &
                                 dot_product(elem1%wSS(1, 1, 1:dof1), phi0(1:dof1,  n )) , &
                                 sqrt( dot_product(elem1%wSS(1, 1, 1:dof1), Dphi0(1:dof1, 1,  n ))**2 &
                                 + dot_product(elem1%wSS(1, 1, 1:dof1), Dphi0(1:dof1, 2,  n ))**2) , &
                                 dot_product(elem1%wSS(1, 1, 1:dof1), Dphi0(1:dof1, 1,  n )), &
                                 dot_product(elem1%wSS(1, 1, 1:dof1), Dphi0(1:dof1, 2,  n ))

                            write(83, *) Fx(n, 1:2), dot_product(wR(1, 1, 1:dof), phi(1:dof,  n )) , &
                                 sqrt( dot_product(wR(1, 1, 1:dof), Dphi(1:dof, 1,  n ))**2 &
                                 + dot_product(wR(1, 1, 1:dof), Dphi(1:dof, 2,  n ))**2) , &
                                 dot_product(wR(1, 1, 1:dof), Dphi(1:dof, 1,  n )), &
                                 dot_product(wR(1, 1, 1:dof), Dphi(1:dof, 2,  n ))
                            write(84, *) Fx(n, 1:2), dot_product(wR(2, 1, 1:dof), phi(1:dof,  n )) , &
                                 sqrt( dot_product(wR(2, 1, 1:dof), Dphi(1:dof, 1,  n ))**2 &
                                 + dot_product(wR(2, 1, 1:dof), Dphi(1:dof, 2,  n ))**2) , &
                                 dot_product(wR(2, 1, 1:dof), Dphi(1:dof, 1,  n )), &
                                 dot_product(wR(2, 1, 1:dof), Dphi(1:dof, 2,  n ))
                         enddo
                      endif

                      !if(iprint ) &
                      !     write(*,'(a6,2i5,es12.4,a1, 3es12.4,a2,3es12.4)') &
                      !     'hp:',elem1%i, elem1%deg, sqrt(dot_product(elem1%xc, elem1%xc)), '|', &
                      !     sqrt(eK_h), sqrt(eK_p), sqrt(eK_hp), '||', &
                      !     sqrt(eK_hp)/ sqrt(eK_h), sqrt(eK_hp)/ sqrt(eK_p)

                   endif

                   ! adding of DOF for p/refinements
                   if(itype == 2) dof_h_p1 = dof_h_p1 + dof1 - elem1%deg - 1
                      

                enddo

                deallocate(Fx, xi, phi, phiW, weights)
                deallocate(Dphi, DphiW, Dwi, Dphi0, DwiR)

             endif   ! if (weight > 0.) then

          endif  !(i1 > 0)

       enddo ! j=1, elem%flen ( N_elem)

!!!111 continue

       if(itype == 1) then
          ! the H1-norm (comment for the L^2-norm)
          L2weight = 1.0
          Mass(1:dof, 1:dof) = Mass(1:dof, 1:dof) +  L2weight * Stiff(1:dof, 1:dof)
          rhs(1:dof, 1:ndimL) = rhs(1:dof, 1:ndimL) +  L2weight * rhsSff(1:dof, 1:ndimL)

          ! deg (p) projection
          MM(1:dof, 1:dof) = Mass(1:dof, 1:dof)
          bb(1:dof, 1:ndimL) = rhs(1:dof, 1:ndimL)

          call SolveLocalMatrixProblem(dof, MM(1:dof, 1:dof), ndimL, bb(1:dof, 1:ndimL) )
          wR(1, 1:ndimL, 1:dof) = transpose( bb(1:dof, 1:ndimL))

          if(dofM > 0) then
             ! deg (p-) projection
             !MM(1:dofM, 1:dofM) = Mass(1:dofM, 1:dofM)
             !bb(1:dofM, 1:ndimL) = rhs(1:dofM, 1:ndimL)

             !call SolveLocalMatrixProblem(dofM, MM(1:dofM, 1:dofM), ndimL, bb(1:dofM, 1:ndimL) )
             !wR(2, 1:ndimL, 1:dof) = transpose( bb(1:dofM, 1:ndimL))

             wR(2, 1:ndimL, 1:dofM) =  wR(1, 1:ndimL, 1:dofM)
          endif

          if(iprint) then
             call PlotElemFunction3D(91, elem, elem%dof, elem%wSS(1, 0, 1:elem%dof) )
             call PlotElemFunction3D(92, elem, elem%dof, elem%wSS(1, 1, 1:elem%dof) )

             call PlotElemFunction3D(93, elem, elem%dof, wR(1, 1, 1:elem%dof))    
             call PlotElemFunction3D(94, elem, elem%dof, wR(2, 1, 1:elem%dof))    


             call PlotElemFunction3D(95, elem, elem%dof, elem%wSS(1, 0, 1:elem%dof)-elem%wSS(1, 1, 1:elem%dof))
             call PlotElemFunction3D(96, elem, elem%dof, elem%wSS(1, 0, 1:elem%dof)-wR(1, 1, 1:elem%dof))
             call PlotElemFunction3D(97, elem, elem%dof, elem%wSS(1, 0, 1:elem%dof)-wR(2, 1, 1:elem%dof))
             call PlotElemFunction3D(98, elem, elem%dof, wR(1, 1, 1:elem%dof) -wR(2, 1, 1:elem%dof))    

          endif

       endif

       !print*,'###', elem%i, j, sqrt(err_2h_p)

    end do  ! do itype = 1, 2


    err_2h_p = sqrt(err_2h_p)
    err_h_p1 = sqrt(err_h_p1)
    err_2h_p1 = sqrt(err_2h_p1)

    !err_2h_p = sqrt(err_2h_p * elem%area/ area)
    !err_2h_p1 = sqrt(err_2h_p1 * elem%area/ area)
    !err_2h_p = sqrt(errK)
    !err_2h_p1 = sqrt(errKp)

    if(iprint ) then
       write(*,'(a6,2i5,es12.4,a1, 3es12.4,a2,3es12.4)') &
            'hp:',elem%i, elem%deg, sqrt(dot_product(elem%xc, elem%xc)), '|', &
            err_2h_p, err_h_p1, err_2h_p1, '||', &
            err_2h_p1/ err_2h_p, err_2h_p1/ err_h_p1
       !print*,'--------------------------------'
    endif



    deallocate( Mass, Stiff, rhs, rhsSff, MM, bb, wR )


    ! computing of slopes
    rn = 3
    allocate(xi(1:rn, 1:2) )

    ! p-1, h
    xi(1, 1) = err_h_p1                            ! error
    xi(1, 2) = dof_h_p1
    !xi(1, 2) = 1.*(elem%deg + 0) * (elem%deg + 1) / 2    ! DOF

    ! p, 2 * h 
    xi(2, 1) = err_2h_p                            ! error
    xi(2, 2) = dof_2h_p
    !xi(2, 2) = 1.*(elem%deg + 1) * (elem%deg + 2) / 8    ! DOF

    ! p-1, 2 * h
    xi(3, 1) = err_2h_p1                              ! error
    xi(3, 2) = dof_2h_p1
    !xi(3, 2) = 1.*(elem%deg + 0) * (elem%deg + 1) / 8    ! DOF

    !write(*,'(a8, 2i5,3es12.4,3(a2, 3es12.4))') &
    !     'decayL',elem%i, elem%dof, err_2h_p, err_h_p1, err_2h_p1, &
    !     '|', xi(1:3, 2),'|', dof_h_p1, dof_2h_p, dof_2h_p1


    ! EOC for h-refinement
    slope_h = log(xi(1, 1) / xi(3, 1)) /( xi(3, 2)**(1./3) - xi(1, 2)**(1./3) )

    ! EOC for p-refinement
    slope_p = log(xi(2, 1) / xi(3, 1)) /( xi(3, 2)**(1./3) - xi(2, 2)**(1./3) )


     ! write(100 + state%space%adapt%adapt_level, *)  xi(3, 2), err_2h_p1
     ! write(100 + state%space%adapt%adapt_level, *)  xi(1, 2), err_h_p1
     ! write(100 + state%space%adapt%adapt_level, *) '   '

     ! write(200 + state%space%adapt%adapt_level, *)  xi(3,2), err_2h_p1
     ! write(200 + state%space%adapt%adapt_level, *)  xi(2,2), err_2h_p
     ! write(200 + state%space%adapt%adapt_level, *) '   '


    if(iprint ) then
       write(*,'(a6,2i5,es12.4,a1, 3es12.4,a2,3es12.4)') &
            'slopes',elem%i, elem%deg, sqrt(dot_product(elem%xc, elem%xc)), '|', &
            err_2h_p, err_h_p1, err_2h_p1, '||', &
            slope_h, slope_p
       print*,'--------------------------------'
    endif

    !if(elem%i == grid%nelem ) stop "d383d38yd38kdjw"

    deallocate(xi)



  end subroutine Regularity_full_hp_derefinement




  !> indication of the hp-adaptation based on the formal p-derefinement
  !> we evaluate the differences::  
  !> \f$ \| w_h - \Pi^{p-1}_h \w_h\|\f$  
  !> the temporary arrays contains the following values
  !> elem%wSS(1, 0, :) = actual solution
  !> elem%wSS(1, 1, :) = p-1 projection
  subroutine Regularity_only_p_derefinement(elem, err_p_deref)
    class(element), intent(inout) :: elem
    real, intent(inout) :: err_p_deref
    class(element), pointer :: elem1
    type(volume_rule), pointer :: V_rule
    real, dimension(:,:), pointer :: phi0 ! the test functions
    real, dimension(:,:,:), allocatable :: Dphi0 ! the test functions
    real, dimension(:), allocatable :: weights
    real, dimension(:,:), allocatable :: Dwi
    real, dimension(:,:), allocatable :: DwiR  ! reconstructed solution
    integer :: ndimL, deg, dof, Qdof
    integer :: k, n
    real :: err_p

    ndimL = ndim
    deg = elem%deg
    dof = elem%dof

    V_rule => state%space%V_rule(elem%Qnum)
    Qdof = V_rule%Qdof

    allocate (weights(1:Qdof)  )
    call Eval_V_Weights(elem, weights(1:Qdof))
          
          
    !integerals over elem
    phi0 => V_rule%phi(1:dof, 1:Qdof)

    allocate( Dphi0(1:dof, 1:nbDim, 1:Qdof) )
    call Eval_Dphi(elem, dof, Dphi0(1:dof, 1:nbDim, 1:Qdof) )


    ! integration of (ww, \phi)_K, and (D\w, D\phi)_K
    allocate( Dwi(1:Qdof, 1:nbDim) )
    allocate( DwiR(1:Qdof, 1:nbDim) )    ! reconstructed solution in integ nodes

    err_p_deref = 0.

    do k=1, ndimL
          
       ! product (Dwi, Dphi0)
       do n=1,nbDim
          Dwi(1:Qdof, n) = matmul(elem%wSS(1, 0, (k-1)*dof+1 : k*dof), Dphi0(1:dof, n, 1:Qdof))
       enddo
          

       ! p-derefinement, projection into p-1
       do n=1,nbDim
          DwiR(1:Qdof, n) = matmul(elem%wSS(1, 1, (k-1)*dof+1 : k*dof), Dphi0(1:dof, n, 1:Qdof))
       enddo
       err_p =  dot_product(weights(1:Qdof) , &
            ( Dwi(1:Qdof, 1) - DwiR(1:Qdof, 1) ) **2  + ( Dwi(1:Qdof, 2) - DwiR(1:Qdof, 2) ) **2 )


       err_p_deref = err_p_deref + err_p
    enddo
  

    err_p_deref = sqrt(err_p_deref)

    deallocate(weights, Dwi, DwiR)


  end subroutine Regularity_only_p_derefinement


end module regularity
