program SetOrder
  integer:: maxmesh, ival, ivar     
  real, dimension (:,:), allocatable :: x, y, r, orders
  integer, dimension(:), allocatable :: ix, idof
  integer, dimension(:,:), allocatable :: ipar
  logical :: logEOC
  integer :: ife, ife_b, ife_e, ife_txt, iEOC, iline
  character*50 text(1:30)
  character*7  EOC
  character*15 typ
  real :: told
  
  ! ival = number of real values in one row in order.dat at the input 
  ! ncounts = number of real columns in tables at the output
  logEOC  = .false.

  iEOC = 3

  !iline = 30  ! after "iline" lines the  \hline is plotted
  !iline = 21  ! after "iline" lines the  \hline is plotted
  iline = 16  ! after "iline" lines the  \hline is plotted
  !iline = 5  ! after "iline" lines the  \hline is plotted

  if (command_argument_count() == 2) then
    call get_command_argument(1,typ)
    open(12,file='smaz', status='unknown')
    write(12,*) typ
    close(12)
    open(12,file='smaz', status='OLD')
    read(12,*) indexi
    close(12)

    call get_command_argument(2,typ)

    !print*,'%Insert dependence variable: '
    !print*,'%   1 (h_max), 2 (h_aver), 3 (tau),  4 (DOF), 5(ST DOF), 6(p)'
    !!read(*,*) indexi
    !write(*,*) ' Dependence variable = ', indexi
  else
    print*,'% Syntax:  Setorderx  <indexi> <type_of_table>'
    print*,'% indexi:  1 (h_max), 2 (h_aver), 3 (tau),  4 (DOF), 5(ST DOF), 6(p)'
    print*,'% type of table: RES, RES-ST-book, AMA, AMApap, AMAst.... (see .f90)'
    stop
 endif

 !typ = 'STDGM'     ! new ouput for ST DGM 
 !typ = 'RES'
 !typ = 'RES-ST-book'   ! space - time error estimates vortex for BOOK
 !typ = 'RESoT'   ! space - time error estimates
 !typ = 'RESoTpap' ! space - time error estimates, paper BDF-DGFEM IJNMF - time convergence
 !typ = 'RESoSpap' ! space - time error estimates, paper BDF-DGFEM IJNMF  - space convrgence
 !typ = 'RESoSS'   ! space - time error estimates
 !typ = 'RESo'    ! space RES estimation -- L2, H1, X, J, \Pi J 
 !typ = 'RESo2'    ! space RES estimation -- L2, H1, X, J,
 !typ = 'RESo3'    ! space RES estimation -- L2, H1, X, J for the paper hp-steady3
 !typ = 'RESo4'    ! space RES estimation -- L2, H1, X, J for the paper hp-steady4
 !typ = 'AMA' 
 !!typ = 'AMApap' 
 !!typ = 'AMAst' 
 !!typ = 'CFD-AMA' 
 !typ = 'book-hp'
 !typ = 'book-IPG'
 !typ = 'book-IPG2'
 !typ = 'cDLM'
 !typ = 'RESa'   ! with algebraic error
 !typ = 'RESo1'    ! space RES estimation -- L2, X
 !typ = 'RESbeam'
 !typ = 'RTN'
 !typ = 'HEL'
 !typ = 'HELo'    ! paper ENUMATH2011
 !typ = 'DUA'
 !typ = 'DUAo1'
 !typ = 'DUAo'
 !typ = 'DUAp'
 !typ = 'DUAtot'
 !typ = 'DUAtotO'
 !typ = 'DUApapO'   ! paper D. Ern Vohralik SINUM
 !typ = 'DUApapO1'   ! paper D. Ern Vohralik SINUM - revisison

 
  if(typ == 'RTN') then
     ncounts = 14 + 0    
     ival = 14
  endif
  
  if(typ == 'HEL') then
     ncounts = 13 + 0    
     ival = 13
  endif

  if(typ == 'HELo') then ! number of columns in order.dat
     ncounts = 6    
     ival = 13
  endif
  
  if(typ == 'DUA' .or. typ == 'DUAo1' ) then
     ncounts = 12 + 2 + 0    
     ival = 14
     if(typ == 'DUAo1' ) then
        !logEOC  = .true.
        iEOC = 12
     endif

  endif
  
  if(typ == 'DUAo') then
     ncounts = 8 + 2
     ival = 14
     logEOC  = .true.
     iEOC = 8
  endif

  if(typ == 'DUAtot') then
     ncounts = 13 + 2 + 0    
     ival = 18
  endif
  
  if(typ == 'DUAtotO') then
     ncounts = 13 + 2 + 0    
     ival = 18
     !logEOC  = .true.
     iEOC = 13
  endif
  
  if(typ == 'DUApapO') then
     ncounts = 6 + 1 + 0        ! number of columns in tab.tex 
     ival = 21   ! number of columns in order.dat minus 7 
     !logEOC  = .true.
     iEOC = 13
  endif
  
  if(typ == 'DUApapO1') then
     !ncounts = 12 + 1 + 0        ! number of columns in tab.tex 
     ncounts = 9        ! number of columns in tab.tex 
     ival = 23   ! number of columns in order.dat minus 7 
     !logEOC  = .true.
     iEOC = 23
  endif
  
 if(typ == 'DUAp') then
     ncounts = 7
     ival = 14
  endif
  

  if(typ == 'RESoTpap') then
     ncounts = 7 + 1 !+ 2       ! number of columns in tab.tex 
     ival = 20    ! number of columns in order.dat
  endif

  if(typ == 'RESoSpap') then
     ncounts = 7 + 1 + 1       ! number of columns in tab.tex 
     ival = 20    ! number of columns in order.dat
  endif

  if(typ == 'STDGM') then
     ncounts = 7       ! number of columns in tab.tex 
     ival = 7    ! number of columns in order.dat
     iEOC = 7
     logEOC  = .true.
  endif

  if(typ == 'RESo'  ) then
     ncounts = 10
     ival = 19 + 1 ! ratio J_p / J_{p-1}
     logEOC  = .true.
     iEOC = 6
  endif

  if(typ == 'RESo2'  ) then
     ncounts = 10   ! number of columns in tab.tex 
     ival = 19 + 2 ! number of columns in order.dat
     logEOC  = .true.
     iEOC = 7
  endif

  if(typ == 'RESo3'  ) then
     ncounts = 5   ! number of columns in tab.tex 
     ival = 20 + 2 ! number of columns in order.dat
     !ival = 27 + 2 ! number of columns in order.dat
     logEOC  = .true.
     iEOC = 3
  endif

  if(typ == 'RES-ST-book') then
     ncounts = 4  + 0     ! number of columns in tab.tex 
     ival = 20    ! number of columns in order.dat
     logEOC  = .false.
     iEOC = 4

  endif

  if(typ == 'RES' .or. typ == 'RESoT' .or. typ == 'RESoSS') then
     ncounts = 16       ! number of columns in tab.tex 
     ival = 20    ! number of columns in order.dat
     logEOC  = .false.
     !!if(typ == 'RESoT' ) logEOC = .true.  !!! CORRECT !!!
     iEOC = 11

  endif

  if(typ == 'pNeu0') then
     ncounts = 10       ! number of columns in tab.tex 
     ival = 26    ! number of columns in order.dat
     logEOC  = .false.
     !!if(typ == 'RESoT' ) logEOC = .true.  !!! CORRECT !!!
     iEOC = 10
  endif

  if(typ == 'pNeu') then
     ncounts = 18       ! number of columns in tab.tex 
     ival = 26    ! number of columns in order.dat
     logEOC  = .false.
     !!if(typ == 'RESoT' ) logEOC = .true.  !!! CORRECT !!!
     iEOC = 10
  endif

  if(typ == 'RESo4'  ) then
     ncounts = 6   ! number of columns in tab.tex 
     ival = 20 + 2 ! number of columns in order.dat
     !ival = 27 + 2 ! number of columns in order.dat
     logEOC  = .true.
     iEOC = 3
  endif

  if(typ == 'AMA'  ) then
     ncounts = 11   ! number of columns in tab.tex 
     ival = 28 ! number of columns in order.dat
     logEOC  = .true.
     iEOC = 7
  endif


  if(typ == 'AMAst'  ) then
     ncounts = 8   ! number of columns in tab.tex 
     ival = 39 ! number of columns in order.dat
     logEOC  = .false.
     iEOC = 7
     iline = 1000000
  endif

  if(typ == 'AMApap'  ) then
     ncounts = 4   ! number of columns in tab.tex 
     ival = 29 ! number of columns in order.dat
     logEOC  = .true.
     iEOC = 4
  endif

  if(typ == 'CFD-AMA'  ) then
     ncounts = 6   ! number of columns in tab.tex 
     ival = 36 ! number of columns in order.dat
     logEOC  = .false.
     iEOC = 0
  endif

  if(typ == 'book-hp'  ) then
     ncounts = 4   ! number of columns in tab.tex 
     ival = 29 ! number of columns in order.dat
     logEOC  = .true.
     iEOC = 3
  endif

  if(typ == 'book-IPG'  ) then
     ncounts = 4   ! number of columns in tab.tex 
     ival = 3 ! number of columns in order.dat
     logEOC  = .true.
     iEOC = 3
  endif
  if(typ == 'book-IPG2'  ) then
     ncounts = 6   ! number of columns in tab.tex 
     ival = 6 ! number of columns in order.dat
     logEOC  = .true.
     iEOC = 6
  endif

  if(typ == 'cDLM'  ) then
     ncounts = 7   ! number of columns in tab.tex 
     ival = 36 ! number of columns in order.dat
     logEOC  = .false.
     iEOC = 4
  endif

  if(typ == 'RESo1' ) then
     ncounts = 4    
     ival = 19 
     logEOC  = .true.
     iEOC = 2
  endif

  if(typ == 'RESa'  ) then
     ncounts = 7    
     ival = 16 
     !logEOC  = .true.
     logEOC  = .false.
     iEOC = 4
  endif

  if(typ == 'RESbeam') then
     ncounts = 5    
     ival = 13 

     logEOC  = .false.
  endif


  !print*,'####',ival, typ
  !stop
  
  maxmesh = 500 
  ivar = 4
  allocate( x(maxmesh,ival+ivar), y(maxmesh,ival+ivar),r(ival,3), ix(maxmesh) )
  allocate( idof(maxmesh), ipar(maxmesh,3), orders(0:maxmesh, ival) )
  
  if(typ == 'RTN') then
     text(1) = "$\|e_h(T)\|_0$"
     text(2) = "$|e_h(T)|_1$"
     text(3) = "$\|e_h\|_X$"
     text(4) = "$\|\partial_t e_h\|_{X'}$"
     text(5) = "$\|e_h\|_Y$"
     text(6) = "$\eta_{DF} $"
     text(7) = "$\eta_{R} $"
     text(8) = "$\eta_{DFR} $"
     text(9) = "$\eta_{NC1} $"
     text(10) = "$\eta_{NC2} $"
     text(11) = "$\eta_{Osc} $"
     text(12) = "$\eta_{IC} $"
     text(13) = "$\eta_{tot} $"
     text(14) = "$i_{\rm eff} $"
     text(15) = "CPU(s)"
     !text(15) = "$\eta_{tot}^{old}$)"
  end if

  if(typ == 'HEL') then
     text(1) = "$\|e_h(T)\|_0$"
     text(2) = "$|e_h(T)|_1$"
     text(3) = "$\|e_h\|_X$"
     text(4) = "$\|e_h\|_Y$"
     !text(4) = "$(\|e_h(T)\|_0^2 + \|e_h\|_X^2)^{1/2}$"
     text(5) = "$\eta_{rez} $"
     text(6) = "$\eta_{h[u_h]^2}=\eta_2 $"
     text(7) = "$\eta_{[\nabla u_h\cdot n]} $"
     text(8) = "$\eta_{[u_h]^2/h} $"
     text(9) = "$\eta_{1} $"
     text(10) = "$\eta_{1+2} $"
     text(11) = "$\eta_{IC} $"
     text(12) = "$\eta_{tot} $"
     text(13) = "$i_{\rm eff} $"
     text(14) = "CPU(s)"  ! not used
  end if

  if(typ == 'HELo') then
     text(1) = "$\|e_h\|_Y$"
     text(2) = "$\eta_{1} $"
     text(3) = "$\eta_{2} $"
     text(4) = "$\eta_{\rm IC} $"
     text(5) = "$\eta_{\rm tot} $"
     text(6) = "$i_{\rm eff} $"
  end if
  
  if(typ == 'DUA' .or. typ == 'DUAo1') then
     text(1) = "$\|e_h(T)\|_0$"
     text(2) = "$|e_h(T)|_1$"
     text(3) = "$\|e_h\|_X$"
     !text(3) = "$\|[e_h]\|_0/h_\Gamma$"
     text(4) = "$\|e_h\|_Y$"
     text(5) = "$\|e_h\|_F$"
     text(6) = "$\|e_h\|_{NC}$"
     text(7) = "$\eta_{DF} $"
     text(8) = "$\eta_{R} $"
     text(9) = "$\eta_{DFR} $"
     text(10) = "$\eta_{NC} $"
     text(11) = "$\eta_{IC} $"
     text(12) = "$\eta_{sum} $"
     !text(13) = "$\eta_{ic} $"
     text(13) = "$i_{{\rm eff}_Y} $"
     text(14) = "$i_{{\rm eff}_F} $"
     !text(13) = "CPU(s)"  ! not used
  end if

  if(typ == 'DUAtot' .or. typ == 'DUAtotO') then
     text(1) = "norm1"
     text(2) = "norm2"
     text(3) = "norm'1"
     text(4) = "norm'2"
     text(5) = "norm'3"
     text(6) = "$\|e_h\|_Y$"
     text(7) = "$\|e_h\|_F$"
     text(8) = "$\|e_h\|_{NC}=\eta_{NC}$"
     text(9) = "$\eta_{DF} $"
     text(10) = "$\eta_{R} $"
     text(11) = "$\eta_{DFR} $"
     !text(12) = "$\eta_{NC} $"
     text(12) = "$\eta_{IC} $"
     text(13) = "$\eta_{sum} $"
     text(14) = "$i_{{\rm eff}_Y} $"
     text(15) = "$i_{{\rm eff}_F} $"
     !text(13) = "CPU(s)"  ! not used
  end if

  if(typ == 'DUApapO') then
     !text(1) = "$\|e_h\|_X$"
     text(1) = "$\|e\|_{\rm FR}$"
     text(2) = "$\|e\|_{\rm NC}=\eta_{NC}$"
     text(3) = "$\eta_{\rm F} $"
     text(4) = "$\eta_{\rm R} $"
     !text(6) = "$\eta_{DFR} $"
     text(5) = "$\eta_{\rm IC} $"
     text(6) = "$\eta $"
     text(7) = "$i_{{\rm eff}} $"
     !text(13) = "CPU(s)"  ! not used
  end if

  if(typ == 'DUApapO1') then
     text(1) = "$J_{u,{\rm FR}}(u_{h\tau}) $"
     !text(2) = "$J_{u,{\rm FR}}(u_{h\tau})^\prime $"
     !text(1) = "$\|e\|_{\rm FR}$"
     text(2) = "$\eta_{\rm F} $"
     text(3) = "$\eta_{\rm R} $"
     text(4) = "$\eta_{\rm NC}$"
     !text(6) = "$\eta_{DFR} $"
     text(5) = "$\eta_{\rm IC} $"
     text(6) = "$\eta_{\rm qd} $"
     text(7) = "$\eta $"
     text(8) = "$i_{{\rm eff}} $"
     !text(10) = "$i_{{\rm eff}}\prime $"
     !text(9) = "$C_{Tn} $"
     text(9) = "$i_{{\rm eff,FR}} $"
     !text(12) = "$N_{\rm splits} $"
     !text(13) = "$\|e_h\|_X $"
     !text(13) = "CPU(s)"  ! not used
  end if
  

  if(typ == 'DUAp') then
     text(1) = "$\|e_h\|_F$"
     text(2) = "$\|e_h\|_{NC}$"
     text(3) = "$\eta_{DF} $"
     text(4) = "$\eta_{R} $"
     text(5) = "$\eta_{sum} $"
     !text(6) = "$\eta_{IC} $"
     text(6) = "$\eta_{tot} $"
     text(7) = "$i_{{\rm eff}_F} $"
     !text(13) = "CPU(s)"  ! not used
  end if

  if(typ == 'DUAo') then
     text(1) = "$\|e_h\|_X$"
     text(2) = "$\|e_h\|_F$"
     text(3) = "$\|e_h\|_{NC}$"
     text(4) = "$\eta_{DF} $"
     text(5) = "$\eta_{R} $"
     text(6) = "$\eta_{DFR} $"
     text(7) = "$\eta_{IC} $"
     text(8) = "$\eta_{sum} $"
     text(9) = "$i_{{\rm eff}_F} $"
     text(10) = "$i_{{\rm eff}_X} $"
     !text(13) = "CPU(s)"  ! not used
  end if
  
  if(typ == 'pNeu' .or. typ == 'pNeu0' ) then
     text(1) = "$\|e_h\|_0$"
     text(2) = "$|e_h|_1$"
     text(3) = "$\|\nabla u_h + \sigma\|$"
     text(4) = "$\frac{h_K}{\pi}\|f-\nabla\cdot\sigma\|$"
     !text(5) = "$\eta_{FR}$"
     text(5) = "$\|\nabla(u_h-s_h)\|$"
     text(6) = "$\eta_{tot}$"
     text(7) = "$i_{\mathrm{eff}}$"
     text(8) = "$\|\nabla s + \sigma\|$"
     text(9) = "$\|\nabla (s -u) \|$"
     text(10) = "$\|\nabla u + \sigma\|$"
  endif
  if(typ == 'pNeu' ) then
     text(11) = "$\| \Pi_{-1}\sigma\|$"
     text(12) = "$\| \Pi_{-2}\sigma\|$"
     text(13) = "$\| \Pi_{-3}\sigma\|$"
     text(14) = "$\| \Pi_{-4}\sigma\|$"
     text(15) = "$\| \Pi_{-1}\nabla s\|$"
     text(16) = "$\| \Pi_{-2}\nabla s\|$"
     text(17) = "$\| \Pi_{-3}\nabla s\|$"
     text(18) = "$\| \Pi_{-4}\nabla s\|$"

  endif

  if(typ == 'RES' .or. typ == 'RESoT' .or. typ == 'RESoSS') then
     text(1) = "$\|e_h\|_0$"
     text(2) = "$\|e_h\|_1$"
     text(3) = "$\|e_h(t)\|_{\varepsilon}$"
     text(4) = "$\sum_{\Gamma}\|u_h\|^2_\Gamma$"   !! "$|e_h|_1$"
     text(5) = "$\|e_h\|_X$"
     text(6) = "$\|e_h\|_Y$"
     text(7) = "$\eta_{\rm A}$"
     text(8) = "$\eta_{\rm S}$"
     text(9) = "$\eta_{\rm T}$"
     text(10) = "$\eta_{\rm ST}$"
     if (typ == 'RESoSS') then
        text(7) = "$\eta_{\rm A}^{\rm ss}$"
        text(8) = "$\eta_{\rm S}^{\rm ss}$"
        text(9) = "$\eta_{\rm T}^{\rm ss}$"
        text(10) = "$\eta_{\rm ST}^{\rm ss}$"
     endif

     text(11) = "$\eta_{\rm time}$"
     text(12) = "$i_{ST/X}$"
     text(13) = "$i_{ST/Y}$"
     text(14) = "$\eta_T/\eta_{time}$"
     text(15) = "$\eta_T/\eta_S$"
     text(16) = "CPU(s)"
  endif


  if(typ == 'RES-ST-book') then
     text(1) = "$\|e_h(T)\|_{L^2(\Omega)}$"
     text(2) = "$|e_h(T)|_{H^1(\Omega)}$"
     text(3) = "$\|e_{h\tau}\|_{L^2(Q_T)}$"
     text(4) = "$\|e_{h\tau}\|_{L^2(0,T;\,H^1(\Omega))}$"
  endif

  if(typ == 'RESoTpap' ) then
     text(1) = "$\errX$"
     text(2) = "$\errY$"
     text(3) = "$\etaA$"
     text(4) = "$\etaS$"
     text(5) = "$\etaT$"
     text(6) = "$\etaST$"
     text(7) = "$\etaT/\etaS$"
     !text(8) = "$i_X$"
     !text(9) = "$i_Y$"
     text(8) = "CPU(s)"
  endif

  if(typ == 'STDGM' ) then
     text(1) = "$\|e_h\|_{L^\infty(L^2)}$"
     text(2) = "$\|e_h\|_{L^2(L^2)}$"
     text(3) = "$\|e_h\|_{L^2(H^1)}$"
     text(4) = "$\|e_h\|_{L^2(H^1_0)}$"
     text(5) = "$\|e_h(t_m)\|_{L^2(\Omega)}$"
     text(6) = "$\|e_h(t_r)\|_{L^2(\Omega)}$"
     text(7) = "$\|e_h(t_i)\|_{L^2(\Omega)}$"
  endif

  if(typ == 'RESoSpap' ) then
     text(1) = "$\errX$"
     text(2) = "$\errY$"
     text(3) = "$\etaA$"
     text(4) = "$\etaS$"
     text(5) = "$\etaT$"
     text(6) = "$\etaST$"
     text(7) = "$\etaT/\etaS$"
     text(8) = "$i_X$"
     !text(9) = "$i_Y$"
     text(9) = "CPU(s)"
  endif

  if(typ == 'RESo' .or. typ == 'RESbeam') then
     text(1) = "$\|e_h\|_{L^2(\Omega)}$"
     text(2) = "$|e_h|_{H^1(\Omega)}$"
     text(3) = "$\|e_h\|_X$"
     text(4) = "$\|e_h\|_J$"
     text(5) = "$\|\Pi e_h\|_J$"
     text(6) = "$\|\Pi e_h\|_J / \| e_h\|_J$ "
     text(7) = "$h^{\mu-1}/p^{s-3/2} $ "
     text(8) = "$\eta(u_h)$"
     text(9) = "$i_{\rm eff} $"
     text(10) = "CPU(s)"
  endif

  if(typ == 'RESo2') then
     text(1) = "$\|e_h\|_{L^2(\Omega)}$"
     text(2) = "$|e_h|_{H^1(\Omega)}$"
     text(3) = "$\|e_h\|_X$"
     text(4) = "$\|e_h\|_J$"
     text(5) = "$\eta(u_h)$"
     text(6) = "$\|e_h\|_{X+J}$"
     text(7) = "$\eta(u_h)_{X+J}$"
     text(8) = "$i_{\rm eff} $"
     text(9) = "$i_{\rm eff}^\prime $"
     text(10) = "CPU(s)"
  endif

  if(typ == 'RESo3') then
     text(1) = "$\ehX$"
     text(2) = "$\ehJ$"
     text(3) = "$\etaX$"
     text(4) = "$i_{\rm eff} $"
     text(5) = "CPU(s)"
  endif
  if(typ == 'RESo4') then
     text(1) = "$\ehX$"
     text(2) = "$\ehJ$"
     text(3) = "$\etaX$"
     text(4) = "$i_{\rm eff}^\rho $"
     text(5) = "$i_{\rm eff}^\eta $"
     text(6) = "CPU(s)"
  endif

  if(typ == 'AMA') then
     text(1) = "$\|e_h\|_\infty$"
     text(2) = "$|E_I(u_h)|_\infty$"
     text(3) = "$\|e_h\|_{L^2(\Omega)}$"
     text(4) = "$|E_I(u_h)|_{L^2}$"
     text(5) = "$\|e_h\|_{H^1(\Omega)}$"
     text(6) = "$\ehX$"
     !!!!!text(6) = "$\ehJ$"
     text(7) = "$\etaX$"
     text(8) = "$i_{\rm eff}^{\infty} $"
     text(9) = "$i_{\rm eff}^\rho $"
     text(10) = "$i_{\rm eff}^\eta $"
     text(11) = "CPU(s)"
  endif

  if(typ == 'AMApap') then
     !text(1) = "$|E_I(u_h)|_{L^2}$"
     text(1) = "$\|\ehp\|_{L^{\infty}}$" !(\Omega)}$"
     text(2) = "$\|\ehp\|_{L^2}$" !(\Omega)}$"
     text(3) = "$\|\ehp\|_{H^1}$" !(\Omega)}$"
     text(4) = "estim "  !!!$|E_I(u_h)|_\infty$"
     !text(6) = "$\ehX$"
     !!!!!text(6) = "$\ehJ$"
     !text(7) = "$\etaX$"
     !text(8) = "$i_{\rm eff}^{\infty} $"
     !text(9) = "$i_{\rm eff}^\rho $"
     !text(10) = "$i_{\rm eff}^\eta $"
     !text(5) = "CPU(s)"
  endif

  if(typ == 'AMAst') then
     !text(1) = "$|E_I(u_h)|_{L^2}$"
     text(1) = "$\|\ehp\|_{L^{\infty}}$" !(\Omega)}$"
     text(2) = "$\|\ehp\|_{L^2}$" !(\Omega)}$"
     text(3) = "$\|\ehp\|_{H^1}$" !(\Omega)}$"
     text(4) = "$\|\ehp\|_{L^2(L^2)}$" !(\Omega)}$"
     text(5) = "$\|\ehp\|_{L^2(H^1)}$" !(\Omega)}$"
     text(6) = "$\|E_I\|_{L^2}$"  !!!$|E_I(u_h)|_\infty$"
     text(7) = "$\|E_I\|_{L^\infty}$"  !!!$|E_I(u_h)|_\infty$"
     text(8) = "CPU(s)"
  endif

  if(typ == 'CFD-AMA') then
     text(1) = "$c_D$"
     text(2) = "$c_L$"
     text(3) = "$c_M$"
     text(4) = "$\|E_I(u_h)\|_{L^\infty}$"
     text(5) = "$\|E_I(u_h)\|_{L^q}$"
     text(6) = "CPU(s)"
  endif

  if(typ == 'book-hp') then
     !text(1) = "$|E_I(u_h)|_{L^2}$"
     text(1) = "$\|e_h\|_{L^{\infty}(\Omega)}$"
     text(2) = "$\|e_h\|_{L^2(\Omega)}$"
     text(3) = "$\|e_h\|_{H^1(\Omega)}$"
     !text(4) = "estim "  !!!$|E_I(u_h)|_\infty$"
     !text(6) = "$\ehX$"
     !!!!!text(6) = "$\ehJ$"
     !text(7) = "$\etaX$"
     !text(8) = "$i_{\rm eff}^{\infty} $"
     !text(9) = "$i_{\rm eff}^\rho $"
     !text(10) = "$i_{\rm eff}^\eta $"
     text(4) = "CPU(s)"
  endif

  if(typ == 'book-IPG') then
     text(1) = "$\|e_h\|_{L^2(\Omega)}$"
     text(2) = "$|e_h|_{H^1(\Omega,\cT_h)}$"
     text(3) = "$|e_h|_{H^1(\Omega)}$"
  endif
  if(typ == 'book-IPG2') then
     text(1) = "$\|e_h\|_{L^2}$"
     text(2) = "$\|e_h\|_{L^2}$"
     text(3) = "$\|e_h\|_{L^2}$"
     text(4) = "$|e_h|_{H^1}$"
     text(5) = "$|e_h|_{H^1}$"
     text(6) = "$|e_h|_{H^1}$"
  endif


  if(typ == 'cDLM') then
     !text(1) = "$|E_I(u_h)|_{L^2}$"
     text(1) = "$c_D$" 
     text(2) = "$\Delta c_D$" 
     text(3) = "$c_L$" 
     text(4) = "$\Delta c_L$" 
     text(5) = "$c_M$" 
     text(6) = "$\Delta c_M$" 
     !text(3) = "$\|\ehp\|_{H^1}$" !(\Omega)}$"
     !text(4) = "estim "  !!!$|E_I(u_h)|_\infty$"
     !text(6) = "$\ehX$"
     !!!!!text(6) = "$\ehJ$"
     !text(7) = "$\etaX$"
     !text(8) = "$i_{\rm eff}^{\infty} $"
     !text(9) = "$i_{\rm eff}^\rho $"
     !text(10) = "$i_{\rm eff}^\eta $"
     text(7) = "CPU(s)"
  endif

  if(typ == 'RESo1' ) then
     text(1) = "$\|e_h\|_{L^2(\Omega)}$"
     !text(2) = "$|e_h|_{H^1(\Omega)}$"
     text(2) = "$\|e_h\|_X$"
     text(3) = "$\eta(u_h)$"
     text(4) = "$i_{\rm eff} $"
  endif

  if(typ == 'RESa' ) then
     text(1) = "$\|e_h\|_{L^2(\Omega)}$"
     text(2) = "$|e_h|_{H^1(\Omega)}$"
     text(3) = "$|e_h|_{Algeb}$"
     text(4) = "LA res"
     text(5) = "$\|e_h\|_X$"
     text(6) = "$\eta(u_h)$"
     text(7) = "$i_{\rm eff} $"
  endif

  
  EOC = 'EOC'
  
  imesh0 =1

  ! reading of the data
  ires = 11
  open(ires, file='order.dat',status='OLD')
  
  do k = 1,maxmesh
     read(ires, *, end=110) ix(k),x(k,1), x(k,3), idof(k), ipar(k,1:3), &
          (x(k,j),j=ivar+1,ival+ivar)

     !write(*,*) ix(k),x(k,1), x(k,3), idof(k), ipar(k,1:3), (x(k,j)
     !,j=ivar+1,ival+ivar)

     if(typ == 'DUApapO1') then
        read(ires,*) x(k, 22 ), ipar(k, 3) ! we rewrite some quanity
        !  !!!!
        !read(ires,*) x(k, 23 ), ipar(k, 3) ! we rewrite some quanity
        ! !!!!
     endif

     if(typ == 'RESo' ) x(k,ival + ivar) = x(k,22) / x(k,8)

     if(typ == 'RESo2'.or. typ == 'RESo3'.or. typ == 'RESo4' ) then
        x(k,ival + ivar - 1) = (x(k, 7)**2 + x(k,8)**2)**0.5
        x(k,ival + ivar - 0) = (x(k,16)**2 + x(k,8)**2)**0.5
     endif
     !write(*, *) ix(k),x(k,1), x(k,3), idof(k), ipar(k,1:3), (x(k,j)
     !,j=ivar+1,ival+ivar)


  enddo
110 imesh = k-1
  print*,'Number of records (lines) =',imesh
  close(ires)
  

  ! average step mesh h
  !   size 12 x 12
  !x(1:imesh,2) = 12./((ix(1:imesh)/2)**0.5)   ! computational domain

  !   size 2 x 2
  !x(1:imesh,2) = 2./((ix(1:imesh)/2)**0.5)   ! computational domain

  !  size 10 x 10
  x(1:imesh,2) = 10./((ix(1:imesh)/2)**0.5)   ! computational domain

  !   size 1 x 1
  !x(1:imesh,2) = 1./((ix(1:imesh)/2)**0.5)   ! computational domain

!    print*,'###',ix(1:imesh)
!  print*,'###',x(1:imesh,2)

  if(indexi == 4) & ! 1/sqrt(DOF)
       x(1:imesh,4) = 1./((idof(1:imesh))**0.5)
  
  if(indexi == 5) then ! 1/(ST DOF)
     x(1:imesh,4) = (x(1:imesh, 3) / idof(1:imesh) ) **(1./3)
     !write(*,'(12es12.4)') x(1:imesh, 3) 
     !write(*,'(12es12.4)') 1.*idof(1:imesh)
     !write(*,'(12es12.4)') (x(1:imesh, 3) / idof(1:imesh))**(-1)
     !write(*,'(12es12.4)') x(1:imesh, 4) 
     indexi = 4
  endif
  if(indexi == 6) then ! p - degree of polynomial approximation
     x(1:imesh,4) = ipar(1:imesh,1)
     indexi = 4
  endif

  ! computing of the logarithm
  y( 1:imesh, 1:ival+ivar ) = log( x( 1:imesh, 1:ival+ivar ) )
  
  
  
  ! computation of the regresion
  do l=1,ival
     ! global
     call Regres(imesh-imesh0 +1, y(imesh0:imesh,indexi),&
          & y(imesh0:imesh, ivar+l), orders(0, l) )
     
     ! locals
     do k=imesh0,imesh-1
        call Regres(2, y(k:k+1,indexi), y(k:k+1, ivar+l), orders(k+1,&
             & l) )
     enddo
  end do
  
  !x(imesh0:imesh, 2) = ix(imesh0:imesh)
  
  ! SQUARE ROOT of number of elements is N 
  !ix(imesh0:imesh) = (ix(imesh0:imesh)/2)**0.5
  !idof(imesh0:imesh) = (idof(imesh0:imesh))**0.5
  
  ife_b = 25
  open (ife_b,file='tab_beg.tex', status='UNKNOWN')
  
  ife_e = 26
  open (ife_e,file='tab_end.tex', status='UNKNOWN')
  
  ife = 22
  open (ife,file='tab.tex', status='UNKNOWN')
  
  ife_txt = 23
  open (ife_txt,file='tab.txt', status='UNKNOWN')
  
  
  if(typ == 'AMA' .or. typ == 'AMApap' .or. typ == 'book-hp' .or. typ == 'pNeu0'&
       .or. typ == 'STDGM' .or. typ == 'AMAst' ) then
     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"
  elseif( typ == 'CFD-AMA') then
     write(ife_b, *) "{\footnotesize           %  vetsi"
     write(ife_b, *) "%{\scriptsize             %  mensi"
     write(ife_b, *) "%{\tiny                   %  nejmensi"
  else
     write(ife_b, *) "%{\footnotesize           %  vetsi"
     write(ife_b, *) "%{\scriptsize             %  mensi"
     write(ife_b, *) "{\tiny                   %  nejmensi"
  endif

  if(typ == 'RTN') then
     !write(ife_b, *) "\begin{tabular}{rc|c|ccccc|ccccccc|c|c}"
     !write(ife_b, *) "\begin{tabular}{rc|c|ccccc|cccccccc|c|c}"
     write(ife_b, *) "\begin{tabular}{rc|c|ccccc|cccccccc|c}"
     
     write(ife_b, *) "\hline"
     write(ife_b, *) "$N$  & $\tau$& $P_k$  "
     
  endif

  if(typ == 'HEL') then
     write(ife_b, *) "\begin{tabular}{rcc|c|cccc|cccccccc|c}"
     write(ife_b, *) "\hline"
     write(ife_b, *) "$N$  & $h$ & $\tau$& $P_k$  "
  endif

  if(typ == 'HELo') then
     write(ife_b, *) "\begin{tabular}{ccc|ccccc|c}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " $P_k$  & $h_m$ & $\tau_m$  "
  endif

  if(typ == 'DUA' .or. typ == 'DUAo1') then
     write(ife_b, *) "\begin{tabular}{rcc|cccccc|cccccc|cc}"
     
     write(ife_b, *) "\hline"
     !write(ife_b, *) "$N$  & $h$ & $\tau$& $P_k$  "
     write(ife_b, *) "$N$  & $\tau$& $P_k$  "
     
  endif

  if(typ == 'DUAtot') then
     write(ife_b, *) "\begin{tabular}{rcc|cccccccc|ccccc|cc}"
     
     write(ife_b, *) "\hline"
     !write(ife_b, *) "$N$  & $h$ & $\tau$& $P_k$  "
     write(ife_b, *) "$N$  & $\tau$& $P_k$  "
     
  endif

  if(typ == 'DUAtotO') then
     write(ife_b, *) "\begin{tabular}{rcc|cccccccc|ccccc|cc}"
     
     write(ife_b, *) "\hline"
     !write(ife_b, *) "$N$  & $h$ & $\tau$& $P_k$  "
     write(ife_b, *) "$N$  & $\tau$& $P_k$  "
     
  endif

  if(typ == 'DUApapO') then
     write(ife_b, *) "\begin{tabular}{ccc|cc|cccc|c}"
     
     write(ife_b, *) "\hline"
     !write(ife_b, *) "$N$  & $h$ & $\tau$& $P_k$  "
     write(ife_b, *) "$h$  & $\tau$& $\mathbb{P}_k$  "
     !write(ife_b, *) "$N$  & $\tau$& $P_k$  "
     
  endif

  if(typ == 'DUApapO1') then
     !!write(ife_b, *) "\begin{tabular}{ccc|cc|ccccc|cc|ccc|c}"
     write(ife_b, *) "\begin{tabular}{cc|c|cccccc|cc}"
     
     write(ife_b, *) "\hline"
     !!write(ife_b, *) "$h$  & $\tau$& $\mathbb{P}_k$  "
     write(ife_b, *) " $m$  & $\mathbb{P}_p$  "
     
  endif

  if(typ == 'DUAo') then
     write(ife_b, *) "&
          &\begin{tabular}{rcc|cc|cc|cc||cc|cc|cc|cc|cc||cc}"
     
     write(ife_b, *) "\hline"
     !write(ife_b, *) "$N$  & $h$ & $\tau$& $P_k$  "
     write(ife_b, *) "$N$  & $\tau$& $P_k$  "
     
  endif

  if(typ == 'DUAp') then
     write(ife_b, *) "\begin{tabular}{rcc|cc|cccc|c}"
     
     write(ife_b, *) "\hline"
     !write(ife_b, *) "$N$  & $h$ & $\tau$& $P_k$  "
     write(ife_b, *) "$N$  & $\tau$& $P_k$  "
     
  endif

  if(typ == 'pNeu0' ) then
     write(ife_b, *) "\begin{tabular}{|ccc|cc|cccc|r|ccc|}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " $h$ &  $P_k$  & dof"
  endif
  if(typ == 'pNeu' ) then
     write(ife_b, *) "\begin{tabular}{|ccc|cc|cccc|r|ccc|cccc|cccc|}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " $h$ &  $P_k$  & dof "
  endif

  if(typ == 'RES'  .or. typ == 'RESoT'.or. typ == 'RESoSS' ) then
     !write(ife_b, *) "\begin{tabular}{cr|cccc|cccc|cccc|c}"
     write(ife_b, *) "\begin{tabular}{cccc|cccc|cc|cccc|c|cccc|r}"
     write(ife_b, *) "\hline"
     !write(ife_b, *) "$N$  & $\tau $ & ${DOF} $   & $P_k$ & $T_m$ "
     write(ife_b, *) " $h$ & $\tau $ &  $P_k$ & $T_m$ "
     !write(ife_b, *) "$N$  & ${DOF}$ & $\tau$& $P_k$ & $T_m$ "
     !write(ife_b, *) "$N$  & $\tau$& $P_k$ & $T_m$ "
  endif

  if(typ == 'RES-ST-book') then 
     write(ife_b, *) "\begin{tabular}{cccc|cc|cc}"
     write(ife_b, *) "\hline"
     !write(ife_b, *) "$N$  & $\tau $ & ${DOF} $   & $P_k$ & $T_m$ "
     write(ife_b, *) " $h$ & $\tau $ &  $P_k$ & $T_m$ "
     !write(ife_b, *) "$N$  & ${DOF}$ & $\tau$& $P_k$ & $T_m$ "
     !write(ife_b, *) "$N$  & $\tau$& $P_k$ & $T_m$ "
  endif

  if(typ == 'STDGM') then 
     write(ife_b, *) "\begin{tabular}{cccc|cc|cc|cc|cc|cc|cc|cc}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " $h$ & $\tau $ &  $P_k$ & $T_m$ "
  endif

  if(typ == 'RESoTpap') then
     write(ife_b, *) "\begin{tabular}{ccc|cc|cccc|cr}"
     write(ife_b, *) "\hline"
     write(ife_b, *) "  $\tau $ &  $P_p$ & $ m $ "
  endif

  if(typ == 'RESoSpap') then
     write(ife_b, *) "\begin{tabular}{ccc|cc|cccc|ccr}"
     write(ife_b, *) "\hline"
     write(ife_b, *) "  $ h $ &  $P_p$ & $ m$ "
  endif

  if(typ == 'RESo') then
     write(ife_b, *) "\begin{tabular}{ccrr|cc|cc|cc|cc|cc|cc|c|cc|c}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " {\rm lev} & $p$ & $\#\Tr$  & $\dof_h $"  !!! &
!!!  $P_k$ & $T_m$ "
     !write(ife_b, *) "$N$  & ${DOF} $"  !!! & $P_k$ & $T_m$ "
     !write(ife_b, *) "$N$  & ${DOF}$ & $\tau$& $P_k$ & $T_m$ "
     !write(ife_b, *) "$N$  & $\tau$& $P_k$ & $T_m$ "
  endif

  if(typ == 'RESo2') then
     write(ife_b, *) "&
          &\begin{tabular}{ccrr|cc|cc||cc|cc|cc||cc|cc||cc|r}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " {\rm lev} & $p$ & $\#\Tr$  & $\dof_h $"  !!! &
!!!  $P_k$ & $T_m$ "
  endif

  if(typ == 'RESo3') then
     write(ife_b, *) "\begin{tabular}{crr||cc|cc|cc|c|r}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " {\rm lev} &  $\#\Tr$  & $\dof_h $"  !!! &
!!!  $P_k$ & $T_m$ "
  endif
  if(typ == 'RESo4') then
     write(ife_b, *) "\begin{tabular}{crr||cc|cc|cc|cc|r}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " {\rm lev} &  $\#\Tr$  & $\dof_h $"  !!! &
!!!  $P_k$ & $T_m$ "
  endif

  if(typ == 'AMA') then
     write(ife_b, *) "&
          &\begin{tabular}{crr||cc|cc|cc|cc|cc|cc|cc|ccc|r}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " {\rm lev} &  $\#\Tr$  & $\dof_h $"  !!! &
!!!  $P_k$ & $T_m$ "
  endif

  if(typ == 'AMApap') then
     write(ife_b, *) "\begin{tabular}{crr||cc|cc|cc|cc}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " $\ell$ &  $N_h$  & $\Nhp $"  !!! &
     !write(ife_b, *) " {\rm lev} &  $\#\Tr$  & $\Nhp $"  !!! &
!!!  $P_k$ & $T_m$ "
  endif

  if(typ == 'AMAst') then
     write(ife_b, *) "\begin{tabular}{crrr||ccc|cc|cc|r}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " $\ell$ &  $N_h$  & $\Nhp $ & $t$"  !!! &
     !write(ife_b, *) " {\rm lev} &  $\#\Tr$  & $\Nhp $"  !!! &
!!!  $P_k$ & $T_m$ "
  endif

  if(typ == 'CFD-AMA') then
     write(ife_b, *) "&
          &\begin{tabular}{crr|ccc|cc|r}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " {\rm lev} &  $\#\Tr$  & $\dof_h $"  !!! &
!!!  $P_k$ & $T_m$ "
  endif


  if(typ == 'book-hp') then
     write(ife_b, *) "\begin{tabular}{ccrr||cc|cc|cc|r}"
     write(ife_b, *) "\hline"
     !write(ife_b, *) " $\ell$ &  $p$ & $N_h$  & $\Nhp $"  !!! &
     write(ife_b, *) " level & \ $p$ \ & $\#\Tr$  &\quad $\mbox{DOF} $\ "  !!! &
  endif

  if(typ == 'book-IPG') then
     write(ife_b, *) "\begin{tabular}{cc|cc|cc|cc}"
     write(ife_b, *) "\hline"
     write(ife_b, *) "  & & \multicolumn{2}{c|}{SIPG} & \multicolumn{2}{c|}{NIPG} & \multicolumn{2}{c}{IIPG} \\"
     write(ife_b, *) "\hline"
     write(ife_b, *) " \ $p$ \ & $h$  &  "  !!! &
  endif

  if(typ == 'book-IPG2') then
     write(ife_b, *) "\begin{tabular}{cc||cc|cc|cc||cc|cc|cc}"
     write(ife_b, *) "\hline"
     write(ife_b, *) "  & & \multicolumn{6}{c||}{$L^2(\Omega)$-norm} & \multicolumn{6}{c|}{$H^1(\Omega,\cT_h)$-seminorm} \\"
     write(ife_b, *) "\hline"
     write(ife_b, *) "  & & \multicolumn{2}{c|}{SIPG} & \multicolumn{2}{c|}{NIPG} & \multicolumn{2}{c||}{IIPG} "
     write(ife_b, *) "  & \multicolumn{2}{c|}{SIPG} & \multicolumn{2}{c|}{NIPG} & \multicolumn{2}{c}{IIPG} \\"
     write(ife_b, *) "\hline"
     write(ife_b, *) " \ $p$ \ & $h$   "  !!! &
     !write(ife_b, *) " \ $p$ \ & $\ell$   "  !!! &
  endif

  if(typ == 'cDLM') then
     write(ife_b, *) "\begin{tabular}{crr||cc|cc|cc|r}"
     write(ife_b, *) "\hline"
     write(ife_b, *) "  $N_h$  & $P_k$ & $\Nhp $"  !!! &
     !write(ife_b, *) " {\rm lev} &  $\#\Tr$  & $\Nhp $"  !!! &
!!!  $P_k$ & $T_m$ "
  endif

  if(typ == 'RESo1') then
     write(ife_b, *) "\begin{tabular}{crr|cc|cc|cc}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " {\rm lev} & $\#\Tr$  & $\dof_h $"  !!! & $P_k$
!!!  & $T_m$ "
     !write(ife_b, *) "$N$  & ${DOF} $"  !!! & $P_k$ & $T_m$ "
     !write(ife_b, *) "$N$  & ${DOF}$ & $\tau$& $P_k$ & $T_m$ "
     !write(ife_b, *) "$N$  & $\tau$& $P_k$ & $T_m$ "
  endif

  if(typ == 'RESa') then
     if(logEOC)       write(ife_b, *) "&
          &\begin{tabular}{crrrr|cc|cc|cc|cc|cc}"
     if(.not. logEOC) write(ife_b, *) "&
          &\begin{tabular}{crrrr|c|ccc|c|cc}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " {\rm lev} & $\#\Tr$  & $\dof_h $ & Mtol & LA&
          & it" !!! & $P_k$ & $T_m$ "
     !write(ife_b, *) "$N$  & ${DOF} $"  !!! & $P_k$ & $T_m$ "
     !write(ife_b, *) "$N$  & ${DOF}$ & $\tau$& $P_k$ & $T_m$ "
     !write(ife_b, *) "$N$  & $\tau$& $P_k$ & $T_m$ "

  endif


  if(typ == 'RESbeam') then
     write(ife_b, *) "\begin{tabular}{crr|cc|cc|c}"
     write(ife_b, *) "\hline"
     write(ife_b, *) " {\rm lev} & $\#\ \Tr$  & $ \dof_h $"  !!! &
!!!  $P_k$ & $T_m$ "
     !write(ife_b, *) "$N$  & ${DOF} $"  !!! & $P_k$ & $T_m$ "
     !write(ife_b, *) "$N$  & ${DOF}$ & $\tau$& $P_k$ & $T_m$ "
     !write(ife_b, *) "$N$  & $\tau$& $P_k$ & $T_m$ "
  endif

  if(typ == 'book-IPG') then
     if(ipar(1, 3) == 0 ) then
        write(ife_b, *)  text(1) , '& EOC &',text(1) , '& EOC &',text(1) , '& EOC'
     else if(ipar(1, 3) == 1 ) then
        write(ife_b, *)  text(2) , '& EOC &',text(2) , '& EOC &',text(2) , '& EOC'
     endif
  else
     do i = 1, ncounts
        !if(logEOC )       write(ife_b, *) "& ", text(i), "& ", EOC
        if(logEOC .and. i <=iEOC ) then
           write(ife_b, *) "& ", text(i), "& ", EOC
        else
           write(ife_b, *) "& ", text(i) 
        endif
     enddo
  end if

  write(ife_b, *) "\\   \hline"
  if(typ == 'book-IPG' .or.typ == 'book-IPG2' ) write(ife_b, *) "  \hline"

  if(typ == 'RESoTpap' .and. ipar(imesh0, 2) == 1)   write(ife, *) " &
       &\hline  % FRD"
  if(typ == 'RESoSpap' .and. ipar(imesh0, 2) == 1)   write(ife, *) " &
       &\hline  % FRD"
  
  i = imesh0
  
88 format(es9.1 )
98 format(es10.2 )
99 format(es10.2, ' & 'es10.2 )
101 format(i8, '\hspace{2mm} & 'es12.3, '\hspace{2mm}' )
103 format(i8, ' & ' i5 )
111 format(i8, ' & 'es11.3 )
102 format(i8, ' & ', i8, ' & ',es9.1 )
151 format(i8, 2(' & 'es10.2) )
100 format(2(' & ', i3) )
80  format( i3 )
90  format(1(' & ', i3) )
95  format(1(' & ', i8) )
  
201 format(  1(' & ', es10.2) )
202 format(  2(' & ', es10.2) )
203 format(  3(' & ', es10.2) )
204 format(  4(' & ', es10.2) )
205 format(  5(' & ', es10.2) )
206 format(  6(' & ', es10.2) )
208 format(  8(' & ', es10.2) )
210 format( 10(' & ', es10.2) )
211 format( 11(' & ', es10.2) )
212 format( 12(' & ', es10.2) )
213 format( 13(' & ', es10.2) )
282 format(  2(' & ', es11.3) )
292 format(  2(' & ', es12.4) )

284 format(  4(' & ', es12.4) )

1282 format(  2(' & {\bf ', es11.3,'}') )

1201 format(  1(' & {\bf ', es10.2,'}') )
1202 format(  2(' & {\bf ', es10.2,'}') )
1203 format(  3(' & {\bf ', es10.2,'}') )
1204 format(  4(' & {\bf ', es10.2,'}') )

251 format(  1(' & ', es9.1) )
252 format(  2(' & ', es9.1) )
253 format(  3(' & ', es9.1) )
254 format(  4(' & ', es9.1) )
1251 format(  1(' & {\bf ', es9.1,'}') )
1252 format(  2(' & {\bf ', es9.1,'}') )
1254 format(  4(' & {\bf ', es9.1,'}') )

266 format( 3(' &', es12.4, '&',es9.1) )

1602 format(  2(' & {\bf ', f6.2,'}') )


301 format( '&', f8.1, '\\' )
309 format( '&', i8, '\\' )
300 format('\\')
  
501 format( 1(' & ', f5.2) )
504 format( 4(' & ', f5.1) )
  
600 format( 1(' & ',f6.2) )
601 format( 1(' & (', f6.2,')') )
602 format( 2(' & (', f6.2,')') )
603 format( 3(' & (', f6.2,')') )
604 format( 4(' & (', f6.2,')') )
605 format( 5(' & (', f6.2,')') )
606 format( 6(' & (', f6.2,')') )
607 format( 7(' & (', f6.2,')') )
608 format( 8(' & (', f6.2,')') )
610 format(10(' & (', f6.2,')') )
612 format(12(' & (', f6.2,')') )

622 format( 2(' & {\tiny(', f7.2,')}') )
623 format( 3(' & {\tiny(', f7.2,')}') )
624 format( 4(' & {\tiny(', f7.2,')}') )
625 format( 5(' & {\tiny(', f7.2,')}') )
627 format( 7(' & {\tiny(', f7.2,')}') )

701 format( 1(' & ', f8.2) )
704 format( 4(' & ', f8.2) )

721 format( 1(' & ', f10.4) )
741 format( 1(' & ', f12.6) )

752 format( ' &\hspace{2mm} ', es12.4, ' & \hspace{2mm}', f8.2,'\hspace{2mm}' )
792 format( ' &\hspace{2mm} ', es12.4, ' & \hspace{2mm}', ' -- \hspace{2mm}' )

742 format( ' & ', es11.1, ' & ', f8.2 )
782 format( ' & ', es11.1, ' & ', ' -- ' )

  
802 format(i8, ' & ', i8 )
803 format(i3, '&' , i8, ' & ', i8 )
804 format(i3, '&' , i8, ' & ', i8, '\ &', i8,'\quad ' )
8041 format(i3, '&' , '--  & ', i8, '\ &', i8,'\quad ')

813 format(i8, '&' , i8, ' & ', i8 )

902 format( ' & ', es10.2, ' & ', f6.2 )
952 format( ' & ', es10.2, ' & ', f7.1 )
992 format( ' & ', es10.2, ' & ', ' -- ' )

999 format( '& & \hspace{0mm}', f8.2, '& & \hspace{0mm}', f8.2, '& & \hspace{2mm}', f8.2, '\hspace{2mm} \\ ' )

  

  do i=imesh0,imesh
     if(typ == 'book-IPG' ) then
        if( i == 1 ) then
           write(ife, *) ' \vspace{-2mm} & & & & & & & \\'     
        endif
     endif

     ! leading data: N, \tau, DOF ,...
     if(typ == 'RTN') then
        write(ife, 101) ix(i), x(i,3)

        ! P_k, T_m, ..
        !write(ife, 100) ipar(i, 1:2)
        write(ife, 90) ipar(i, 1:1)
     elseif(typ == 'HEL') then


        !write(ife, 101) ix(i), x(i,3)
        write(ife, 151) ix(i), x(i,2),x(i,3)

        ! P_k, T_m, ..
        !write(ife, 100) ipar(i, 1:2)
        write(ife, 90) ipar(i, 1:1)

     elseif(typ =='HELo') then

        ! P_k, T_m, ..
        !write(ife, 100) ipar(i, 1:2)
        write(ife, 80) ipar(i, 1:1)

        !write(ife, 101) ix(i), x(i,3)
        !write(ife, 151) ix(i), x(i,2),x(i,3)
        write(ife, 202)  x(i,2),x(i,3)

     elseif(typ == 'DUA' .or. typ == 'DUAo'.or. typ == 'DUAo1' .or.&
          & typ == 'DUAtot' .or. typ == 'DUAtotO') then
        !write(ife, 101) ix(i), x(i,3)
        write(ife, 99) x(i,2), x(i,3)
        !write(ife, 151) ix(i), x(i,1),x(i,3)

        ! P_k, T_m, ..
        !write(ife, 100) ipar(i, 1:2)
        write(ife, 90) ipar(i, 1:1)

     elseif( typ == 'DUAp' .or. typ == 'DUApapO') then
        !write(ife, 101) ix(i), x(i,3)
        write(ife, 99) x(i,1)/2**0.5, x(i,3)
        !write(ife, 100) ipar(i, 1:2)         ! P_k, T_m, ..
        write(ife, 90) ipar(i, 1:1)           ! P_k, ..

     elseif( typ == 'DUApapO1') then
        !!write(ife, 99) x(i,1)/2**0.5, x(i,3)
        write(ife, 80) i 
        write(ife, 90) ipar(i, 1:1)           ! P_k, ..

     elseif( typ == 'pNeu' .or. typ=='pNeu0') then
        write(ife, 88) x(i,1)/2**0.5   ! h
        write(ife, 90) ipar(i, 1:1)           ! P_k, ..
        write(ife, 95) idof(i)           ! dof

     elseif(typ == 'RESo') then
        write(ife, 804) i-1, ipar(i,1), ix(i), idof(i)
     elseif(typ == 'RESo2') then
        write(ife, 804) i-1, ipar(i,1), ix(i), idof(i)
     elseif(typ == 'RESo3' .or. typ == 'RESo4') then
        !ITAB
        write(ife, 803) i-1,  ix(i), idof(i)
        !write(ife, 803)  ipar(i,1), ix(i), idof(i)
     elseif(typ == 'AMA' ) then
        !ITAB
        write(ife, 803) i-1,  ix(i), idof(i)
        !write(ife, 803)  ipar(i,1), ix(i), idof(i)
     elseif(typ == 'AMApap' ) then
        !ITAB
        write(ife, 803) i-1,  ix(i), idof(i)

     elseif(typ == 'AMAst' ) then
        !ITAB
        if( abs(told - x(i,34))> 1E-5 )  write(ife,*) "\hline"
        write(ife, 803) i-1,  ix(i), idof(i)
        write(ife, '(a1, f10.2)') '&', x(i,34)
        told = x(i,34)

     elseif(typ == 'CFD-AMA' ) then
        !ITAB
        write(ife, 803) i-1,  ix(i), idof(i)
        !write(ife, 803)  ipar(i,1), ix(i), idof(i)

     elseif(typ == 'book-IPG' ) then
        !ITAB
        write(ife, 101) ipar(i, 1), x(i,1)

     elseif(typ == 'book-IPG2' ) then
        !ITAB
        !write(ife, 103) ipar(i, 1), i
        write(ife, 101) ipar(i, 1), x(i,1)

     elseif(typ == 'book-hp' ) then
        !ITAB
        if(imesh < 10) then
           write(ife, 804) i-1,  ipar(i, 1:1),  ix(i), idof(i)
        else
           write(ife, 8041) i-1,   ix(i), idof(i)
        endif

     elseif(typ == 'cDLM') then
       write(ife, 813)  ix(i), ipar(i,1), idof(i)

     elseif(typ == 'RESo1' ) then
        write(ife, 803) i-1, ix(i), idof(i)
     elseif(typ == 'RESoTpap') then
        !write(ife, 98) x(i,2)          ! h
        write(ife, 98) x(i,3)          ! tau
        write(ife, 100) ipar(i, 1:2)   ! P_k, T_m, ..
        !write(ife, '(es10.2,2(i3))') x(i,3),ipar(i,1:2)
     elseif(typ == 'RESoSpap') then
        write(ife, 98) x(i,2)          ! h
        !write(ife, 98) x(i,3)          ! tau
        write(ife, 100) ipar(i, 1:2)   ! P_k, T_m, ..
        !write(ife, '(es10.2,2(i3))') x(i,3),ipar(i,1:2)
     elseif(typ == 'RES-ST-book') then
        write(ife, 99) x(i,2), x(i,3)  ! h_aver, tau
        write(ife, 100) ipar(i, 1:2)   ! P_k, T_m, ..
     elseif(typ == 'STDGM') then
        write(ife, 99) x(i,2), x(i,3)  ! h_aver, tau
        write(ife, 100) ipar(i, 1:2)   ! P_k, T_m, ..
     else
        !write(ife, 101) ix(i), x(i,3)
        !write(ife, 803) ix(i),  idof(i)
        !write(ife, 803) i-1, ix(i), idof(i)
        !write(ife, 102) ix(i), idof(i), x(i,3)
        write(ife, 99) x(i,1), x(i,3)  ! h, tau
        !write(ife, 98) x(i,3)  ! tau

        ! P_k, T_m, ..
        write(ife, 100) ipar(i, 1:2)
        !write(ife, 90) ipar(i, 1)

        if(typ == 'RESa') write(ife, 201) x(i,15) ! GMRES tol
        if(typ == 'RESa') write(ife, 95) int(x(i,20) ) ! GMRES iter

     endif

     
     ! particular data RTN
     if(typ == 'RTN') then
        write(ife, 213) x(i,5:17) 
        !write(ife, 501) x(i,17) / x(i,9)
        write(ife, 701) x(i,17) / x(i,9)
        
        ! CPU + end of the line
        !write(ife, 201) x(i,23) 
        write(ife, '(a3, es10.2)') '%& ', x(i,23) 
        write(ife, 300)
        !write(ife, 301) x(i,23) 
     endif
     
     ! particular data HEL
     if(typ == 'HEL') then
        write(ife, 212) x(i,5:16) 
        !write(ife, 501) x(i,16) / x(i,8)
        write(ife, 721) x(i,16) / x(i,8)
        
        ! CPU + end of the line
        !write(ife, 201) x(i,15) 
        write(ife, '(a3, es10.2)') '%& ', x(i,15) 
        write(ife, 300)
        !write(ife, 301) x(i,23) 
     endif

     if(typ == 'HELo') then
        write(ife, 205) x(i,8), x(i,13), x(i,10), x(i,15),x(i,16) 
        write(ife, 721) x(i,16) / x(i,8)
        
        ! CPU + end of the line
        !write(ife, 201) x(i,15) 
        write(ife, '(a3, es10.2)') '%& ', x(i,15) 
        write(ife, 300)
        !write(ife, 301) x(i,23) 
        if(i > 1 ) then
           write(ife,*) '   & \multicolumn{2}{c|}{(EOC)} '
           write(ife, 605 )orders(i,4),orders(i,9),orders(i,6),&
                & orders(i,11),orders(i,12)
           write(ife, 300)
        endif
     endif

     ! particular data DUA
     if(typ == 'DUA' .or. typ == 'DUAo1') then
        write(ife, 212) x(i,5:16) 
        write(ife, 721) x(i,16) / (x(i,8) + x(i,10) )
        write(ife, 721) x(i,16) / (x(i,9) + x(i,10) )

        ! CPU + end of the line
        !write(ife, 201) x(i,15) 
        write(ife, '(a3, es10.2)') '%& ', x(i,23) 
        write(ife, 300)
        !write(ife, 301) x(i,23) 

        if( typ == 'DUAo1' .and. i > 1) then
           write(ife,*) ' &  & (EOC) '
           write(ife, 612) orders(i,1:12)
           write(ife, 300)
        endif
     endif

     ! particular data DUAtot
     if(typ == 'DUAtot' .or. typ == 'DUAtotO' ) then
        write(ife, 208) x(i,17:21), x(i,8:10)   ! eN1, eN2, eNp1,
        !  eNp2, eNp3, Ynorm, Fnorm, NC
        write(ife, 205) x(i,11:13),x(i,15:16)   ! \eta_DF, \eta_R,
        ! \eta_DFR, \eta_IC, \eta_tot
        write(ife, 721) x(i,16) / (x(i,8) + x(i,10) )
        write(ife, 721) x(i,16) / (x(i,9) + x(i,10) )

        ! CPU + end of the line
        !write(ife, 201) x(i,15) 
        write(ife, '(a3, es10.2)') '%& ', x(i,23) 
        write(ife, 300)
        !write(ife, 301) x(i,23) 

        if(typ == 'DUAtotO' .and. i > 1 ) then
           write(ife,*) ' &  & (EOC) '
           write(ife, 608) orders(i,13:17), orders(i,4:6) 
           write(ife, 605) orders(i,7:9), orders(i,11:12) 
           write(ife, 300)
        endif
     endif


     ! particular data DUAtot
     if(typ == 'DUApap' .or. typ == 'DUApapO' ) then
        write(ife, 202)  x(i,9:10)   ! Fnorm, NC
        write(ife, 204) x(i,11:12),x(i,15:16)  ! \eta_DF, \eta_R,
        ! \eta_IC, \eta_tot
        !!!write(ife, 721) x(i,16) / (x(i,8) + x(i,10) )
        write(ife, 721) x(i,16) / (x(i,9) + x(i,10) )

        ! CPU + end of the line
        !write(ife, 201) x(i,15) 
        write(ife, '(a3, es10.2)') '%& ', x(i,23) 
        write(ife, 300)
        !write(ife, 301) x(i,23) 

        if(typ == 'DUApapO' .and. i > 1 ) then
           write(ife,*) ' &  \multicolumn{2}{c|}{(EOC)} '
           write(ife, 602)  orders(i,5:6) 
           write(ife, 604) orders(i,7:8), orders(i,11:12) 
           write(ife, 300)
        endif
     endif

     ! particular data DUApapO1
     if(typ == 'DUApapO1' ) then

        write(ife, 201)  x(i,22)   ! dual norm J_u(u_h\tau)
        !write(ife, 201)  x(i,23)   ! dual norm J_u(u_h\tau)
        !!write(ife, 202)  x(i,9:10)   ! Fnorm, NC
        write(ife, 204) x(i,11:12), x(i, 10), x(i,15)  ! \eta_DF,
        ! \eta_R, \eta_NC, \eta_IC
        write(ife, 201)  x(i,25)   ! quadrature
        write(ife, 201) x(i,16) !!! \eta_tot
        write(ife, 721) x(i,16) / (x(i,22) + x(i,10) ) ! effectivity
        !  index NEW
        !write(ife, 721) x(i,16) / (x(i,23) + x(i,10) ) ! effectivity
        ! index NEW
        write(ife, 721) x(i,16) / (x(i,9) + x(i,10) )  ! effectivity
        !  index OLD
        !write(ife, 201)  x(i,26)   ! CTn
        !write(ife, 90) ipar(i, 3)  ! number of space-time
        ! interelement splittings
        !write(ife, 201)  x(i,7)   ! X-norm = (L^2(0,T; H1)

        ! CPU + end of the line
        !write(ife, 201) x(i,15) 
        write(ife, '(a3, es10.2)') '%& ', x(i,23) 
        write(ife, 300)
        !write(ife, 301) x(i,23) 

        if( i > 1 ) then
           !write(ife,*) ' &  \multicolumn{2}{c|}{(EOC)} '
           write(ife,*) ' &  \multicolumn{1}{c|}{(EOC)} '
           write(ife, 601)  orders(i,18) 
           !write(ife, 601)  orders(i,19) 
           !write(ife, 602)  orders(i,5:6) 
           write(ife, 603) orders(i,7:8), orders(i,6)
           write(ife, 601)  orders(i,11) 
           write(ife, 601)  orders(i,21) 
           write(ife, 601)  orders(i,12) 
           !write(ife, *) ' &  & ' 
           !write(ife, *) ' &  & '
           !write(ife, 601)  orders(i,3) 
           write(ife, 300)
        endif
     endif


     ! particular data DUA
     if(typ == 'DUAp') then
        write(ife, 206) x(i,9:12) , x(i, 15), x(i, 16)
        !write(ife, 721) x(i,16) / x(i,8)
        write(ife, 721) x(i,16) / (x(i,9) + x(i,10) )
        
        ! CPU + end of the line
        !write(ife, 201) x(i,15) 
        write(ife, '(a3, es10.2)') '%& ', x(i,23) 
        write(ife, 300)
        !write(ife, 301) x(i,23) 
     endif
     
     ! particular data  RES
     if(typ == 'RES' .or. typ == 'RESoT' .or. typ == 'RESoSS')  then

        !write(ife, 206) x(i,5:10)         ! errors
        write(ife, 204) x(i,5:8)         ! errors
        write(ife, 292) x(i,9:10)         ! errors

        if( typ /= 'RESoSS') then
           write(ife, 204) x(i,11:14)          ! estimates
        else
           write(ife, 204) x(i,15:18)          ! estimates
        endif
        !write(ife, 704) x(i,9) / x(i,5), x(i,10) / (x(i,5)**2 + x(i
        !,6)**2)**0.5, x(i,11) / x(i,7), x(i,12) / x(i,8)         !
        ! effectivity indexes

        write(ife, 201) x(i,23)         ! time error
        write(ife, 701) x(i,14) / x(i,9 )  ! index ST / X
        write(ife, 701) x(i,14) / x(i,10)  ! index ST / Y
        write(ife, 201) x(i,13) / x(i,23)  ! index T  / time error
        write(ife, 201) x(i,13) / x(i,12)  ! index T  / S
        
        ! end with CPU time
        write(ife, 301) x(i,24) 
        ! end without CPU time
        !write(ife, 300)

        if(typ == 'RESoT' .and. i > 1)  then
           write(ife,*) ' &  & \multicolumn{2}{c|}{(EOC)} '
           !write(ife,*) ' & &  (EOC) '
           write(ife, 610) orders(i,1:10)
           write(ife, 601) orders(i,19)
           write(ife,*) ' & &  & & &'
           write(ife, 300)
        endif

        if(typ == 'RESoSS' .and. i > 1)  then
           write(ife,*) ' &  \multicolumn{2}{c|}{(EOC)} '
           !write(ife,*) ' & &  (EOC) '
           write(ife, 610) orders(i,1:6), orders(i,11:14)
           write(ife, 300)
        endif

     endif



     ! particular data  RES
     if(typ == 'pNeu' .or. typ == 'pNeu0')  then

        write(ife, 202) x(i,5:6)         ! errors
        !write(ife, 205) x(i,15:19)        ! estims
        write(ife, 202) x(i,15:16)        ! estims
        write(ife, 202) x(i,18:19)        ! estims

        ! effectivity indexes

        write(ife, 701) x(i,19) / x(i,6 )  ! index 

        write(ife, 203) x(i,20:22)        ! estims
        if(typ == 'pNeu' ) then
           write(ife, 204) x(i,23:26)        ! estims
           write(ife, 204) x(i,27:30)        ! estims
        endif

        write(ife, 300)

        if( i > 1)  then
           write(ife,*) '  \multicolumn{3}{|c|}{{\scriptsize (EOC)}} '
           !write(ife,*) ' & &  (EOC) '
           !write(ife, 627) orders(i,1:2), orders(i,11:15)
           write(ife, 622) orders(i,1:2)
           write(ife, 622) orders(i,11:12)
           write(ife, 622) orders(i,14:15)
           write(ife,*) ' & '
           write(ife, 623) orders(i,16:18)
           if(typ == 'pNeu' ) then
              write(ife, 624) orders(i,19:22)
              write(ife, 624) orders(i,23:26)
           endif

           write(ife, 300)
        endif

     endif


     ! particular data  RES
     if(typ == 'RES-ST-book') then 

        !if( i == 1)  then
        !   write(ife, 992) x(i,5)
        !   write(ife, 992) x(i,6)
        !   write(ife, 992) x(i,10)
        !   write(ife, 992) x(i, 9)
        !else
        !   write(ife, 902) x(i,5), orders(i,1)
        !   write(ife, 902) x(i,6), orders(i,2)
        !   write(ife, 902) x(i,10), orders(i,4)
        !   write(ife, 902) x(i, 9), orders(i,3)
        !endif

        !write(ife, 206) x(i,5:10)         ! errors
        write(ife, 202) x(i,5:6)         ! errors
        write(ife, 202) x(i,10), x(i,9)         ! errors

        
        ! end with CPU time
        !write(ife, 301) x(i,24) 
        ! end without CPU time
        write(ife, 300)

        if( i > 1)  then
           write(ife,*) ' &  & \multicolumn{2}{c|}{(EOC)} '
           !write(ife,*) ' & &  (EOC) '
           write(ife, 604) orders(i,1:2), orders(i,6), orders(i,5)
           write(ife, 300)
        endif

     endif


     ! particular data  STDGM
     if(typ == 'STDGM')  then
        if(i > 1) then
           write(ife, 952) x(i,5), orders(i,1)  ! L^\infty(L^2)
           write(ife, 952) x(i,6), orders(i,2)  ! L^2(L^2)
           write(ife, 952) x(i,7), orders(i,3)  ! L^2(H^1)
           write(ife, 952) x(i,8), orders(i,4)  ! L^2(H^1_0) 
           write(ife, 952) x(i,9), orders(i,5)  ! 
           write(ife, 952) x(i,10), orders(i,6)  !
           write(ife, 952) x(i,11), orders(i,7)  ! L^2(H^1_0) 

        else
           write(ife, 992) x(i,5)
           write(ife, 992) x(i,6)
           write(ife, 992) x(i,7)
           write(ife, 992) x(i,8)
           write(ife, 992) x(i,9) 
           write(ife, 992) x(i,10)
           write(ife, 992) x(i,11)
        endif

        ! end with CPU time
        !write(ife, 309) int(x(i,24) )
        ! end without CPU time
        write(ife, 300)

     endif

     if(typ == 'RESoTpap' .or. typ == 'RESoSpap')  then
        !write(ife, 292) x(i,9:10)         ! errors
        if(x(i,13) / x(i,12) < -1E-02) then
        !if(x(i,13) / x(i,12) < 1E-02) then  !! BOLTED FONT
           !write(ife, 1282) x(i,9:10)         ! errors
           write(ife, 1202) x(i,9:10)         ! errors
           write(ife, 1204) x(i,11:14)          ! estimates
           write(ife, 1201) x(i,13) / x(i,12)  ! index T  / S
           !write(ife, 1602) x(i,14)/x(i,9), x(i,14)/x(i,10)    ! eff
           ! index 
        else
           !write(ife, 282) x(i,9:10)         ! errors
           write(ife, 202) x(i,9:10)         ! errors
           write(ife, 204) x(i,11:14)          ! estimates
           write(ife, 201) x(i,13) / x(i,12)  ! index T  / S
           !write(ife, 602) x(i,14)/x(i,9), x(i,14)/x(i,10)    ! eff
           ! index 
        endif

        !write(ife, 204) x(i,11:14)          ! estimates
        !write(ife, 201) x(i,23)         ! time error
        if(typ == 'RESoSpap') then
           write(ife, 701) x(i,14) / x(i,9 )  ! index ST / X
           !!write(ife, 701) x(i,14) / x(i,10)  ! index ST / Y
        endif

        !write(ife, 201) x(i,13) / x(i,23)  ! index T  / time error
        !write(ife, 201) x(i,13) / x(i,12)  ! index T  / S
        
        ! end with CPU time
        write(ife, 309) int(x(i,24) )
        ! end without CPU time
        !write(ife, 300)

        !if(typ == 'RESoT' .and. i > 1)  then
        !   write(ife,*) ' &  & \multicolumn{2}{c|}{(EOC)} '
        !   !write(ife,*) ' & &  (EOC) '
        !   write(ife, 610) orders(i,1:10)
        !   write(ife, 601) orders(i,19)
        !   write(ife,*) ' & &  & & &'
        !   write(ife, 300)
        !endif

     endif

     if( typ == 'DUAo') then
        !if(logEOC) then
           ! errors
           if(i > 1) then
              write(ife, 902) x(i,7), orders(i,3)  ! X
              write(ife, 902) x(i,9), orders(i,5)  ! F
              write(ife, 902) x(i,10), orders(i,6)  ! NC
              write(ife, 902) x(i,11), orders(i,7)  ! DF
              write(ife, 902) x(i,12), orders(i,8)  ! R
              write(ife, 902) x(i,13), orders(i,9)  ! DFR
              write(ife, 902) x(i,15), orders(i,11)  ! IC
              write(ife, 902) x(i,16), orders(i,12)  ! tot
           else
              write(ife, 992) x(i,7)
              write(ife, 992) x(i,9)
              write(ife, 992) x(i,10)
              write(ife, 992) x(i,11)
              write(ife, 992) x(i,12)
              write(ife, 992) x(i,13)
              write(ife, 992) x(i,15)
              write(ife, 992) x(i,16)
           endif
           write(ife, 721) x(i,16) / (x(i,9) + x(i,10) )
           write(ife, 721) x(i,16) /  x(i,7) 
           
           write(ife, '(a3, es10.2)') '%& ', x(i,23) 
           write(ife, 300)

        !endif
     endif

     if( typ == 'RESa') then
        if(logEOC) then
           ! errors
           if(i > 1) then
              write(ife, 902) x(i,5), orders(i,1)  ! L2
              write(ife, 902) x(i,6), orders(i,2)  ! H1
              write(ife, 902) x(i,16),orders(i,12)  ! algeb
              write(ife, 902) x(i,7), orders(i,3)  ! X
           else
              write(ife, 992) x(i,5)
              write(ife, 992) x(i,6)
              write(ife, 992) x(i,16) !algeb
              write(ife, 992) x(i,7)
           endif
        else

           write(ife, 205) x(i,5), x(i,6), x(i,16), x(i,17), x(i,7)
           
        endif

        write(ife, 902) x(i,11), x(i,11) / x(i,7)  ! \eta, i_{eff}
        
        write(ife, 300) 

     endif


     if(typ == 'RESo' .or. typ == 'RESo1') then
        ! errors
        if(i > 1) then
           !if(typ == 'RESo') 
           write(ife, 902) x(i,5), orders(i,1)  ! L2
           if(typ == 'RESo' ) write(ife, 902) x(i,6), orders(i,2)  !
           !  H1
           write(ife, 902) x(i,7), orders(i,3)  ! X
           if(typ == 'RESo' ) write(ife, 902) x(i,8), orders(i,4)   !J
           if(typ == 'RESo' ) write(ife, 902) x(i,22), orders(i,18)&
                & !Pi J
           if(typ == 'RESo' ) write(ife, 902) x(i,24), orders(i,20)&
           & !J/ Pi J
        else
           !if(typ == 'RESo') 
           write(ife, 992) x(i,5)
           if(typ == 'RESo' ) write(ife, 992) x(i,6)
           write(ife, 992) x(i,7)
           if(typ == 'RESo' ) write(ife, 992) x(i,8)
           if(typ == 'RESo' ) write(ife, 992) x(i,22)
           if(typ == 'RESo' ) write(ife, 902) x(i,24) ! J/ Pi J
        endif

        if(typ == 'RESo' ) write(ife, 201) (x(i,1)**(min (ipar(i,1)&
             &+1., 5.25)-1)/(ipar(i,1)**(5.25-3./2))) / x(i,8)

        !print*,'$$$$',x(i,1),ipar(i,1), min (ipar(i,1)+1., 5.25)-1,
        ! x(i,1)**(min (ipar(i,1)+1., 5.25)-1), (x(i,1)**(min (ipar(i
        ! ,1)+1., 5.25)-1)/(ipar(i,1)**(5.25-3./2))) 

        write(ife, 902) x(i,16), x(i,16) / x(i,7)  ! \eta, i_{eff}


        if(typ == 'RESo' ) then
           write(ife, 301) x(i,23) 
        else
           write(ife, 300) 
        endif

     endif

     if(typ == 'RESo2' .or. typ == 'RESo3'.or. typ == 'RESo4' ) then
        ! errors
        if(i > 1) then
           if (typ == 'RESo2')  write(ife, 902) x(i,5), orders(i,1) &
                & ! L2
           if (typ == 'RESo2')  write(ife, 902) x(i,6), orders(i,2) &
           & ! H1
           write(ife, 902) x(i,7), orders(i,3)  ! X
           write(ife, 902) x(i,8), orders(i,4)   !J
           !write(ife, 902) x(i,22), orders(i,18) !Pi J
           !write(ife, 902) x(i,24), orders(i,20) !J/ Pi J
           
            write(ife, 902) x(i,16), orders(i,12)   !eta

            if (typ == 'RESo2')  write(ife, 902) x(i,24), orders(i&
                 &,20)  ! ||e_h||_X  +  ||e_h||_J
            if (typ == 'RESo2')  write(ife, 902) x(i,25), orders(i&
                 &,21)  ! eta   +  ||e_h||_J

        else
           if (typ == 'RESo2')  write(ife, 992) x(i,5)
           if (typ == 'RESo2')  write(ife, 992) x(i,6)
           write(ife, 992) x(i,7)
           write(ife, 992) x(i,8)
           !write(ife, 992) x(i,22)
           write(ife, 992) x(i,16) !eta

           if (typ == 'RESo2')  write(ife, 992) x(i,ival + ivar - 1)&
                & ! ||e_h||_X  +  ||e_h||_J
           if (typ == 'RESo2')  write(ife, 992) x(i,ival + ivar - 0)&
           & ! eta   +  ||e_h||_J

        endif


        if (typ == 'RESo2')  write(ife, 601) x(i,16) / x(i,7)  ! \eta
        ! , i_{eff}


        if (typ == 'RESo4') write(ife, 600) x(i, 16) / x(i,7) !
        !  i_{eff}

        write(ife, 600) (x(i, 16)**2 + x(i,8)**2)**0.5 / (x(i,7)**2 +&
             & x(i,8)**2)**0.5 ! i_{eff}

        write(ife, 301) x(i,24) 

        !write(*,'(i5,5es12.4)') i, x(i,22:26)
        !write(*,'(i5,10es12.4)') i, x(i,7),x(i,8),x(i,16), (x(k, 16)
        !**2 + x(k,8)**2)**0.5 , (x(k,7)**2 + x(k,8)**2)**0.5, (x(k,
        ! 16)**2 + x(k,8)**2)**0.5 / (x(k,7)**2 + x(k,8)**2)**0.5
     endif


     if(typ == 'AMA' ) then
        ! errors
        if(i > 1) then
           write(ife, 902) x(i,27), orders(i,23)  ! L^\infty
           write(ife, 902) x(i,31), orders(i,27)  ! E_I L^\infty
           write(ife, 902) x(i,5), orders(i,1)  ! L^2
           write(ife, 902) x(i,32), orders(i,28)  ! E_I L^2
           write(ife, 902) x(i,6), orders(i,2)   !J
           write(ife, 902) x(i,7), orders(i,3)  ! X
           !write(ife, 902) x(i,8), orders(i,4)   !J
           !write(ife, 902) x(i,22), orders(i,18) !Pi J
           !write(ife, 902) x(i,24), orders(i,20) !J/ Pi J

           write(ife, 902) x(i,16), orders(i,12)   !eta


        else
           write(ife, 992) x(i,27)
           write(ife, 992) x(i,31)

           write(ife, 992) x(i,5)
           write(ife, 992) x(i,32)

           write(ife, 992) x(i,6)
           write(ife, 992) x(i,7)
           !write(ife, 992) x(i,8)
           !write(ife, 992) x(i,22)
           write(ife, 992) x(i,16) !eta
        endif


        write(ife, 600) x(i, 31) / x(i,27) ! i_{eff}


        write(ife, 600) x(i, 16) / x(i,7) ! i_{eff}

        write(ife, 600) (x(i, 16)**2 + x(i,8)**2)**0.5 / (x(i,7)**2 + x(i,8)**2)**0.5 ! i_{eff}

        write(ife, 301) x(i,24) 

        !write(*,'(i5,5es12.4)') i, x(i,22:26)
        !write(*,'(i5,10es12.4)') i, x(i,7),x(i,8),x(i,16), (x(k, 16)**2 + x(k,8)**2)**0.5 , &
        !      (x(k,7)**2 + x(k,8)**2)**0.5, &
        !       (x(k, 16)**2 + x(k,8)**2)**0.5 / (x(k,7)**2 + x(k,8)**2)**0.5
     endif

     if(typ == 'AMApap' ) then
        ! errors
        if(i > 1) then
           !write(ife, 902) x(i,32), orders(i,28)  ! E_I L^2
           !write(ife, 902) x(i,31), orders(i,27)  ! E_I L^\infty
           write(ife, 952) x(i,27), orders(i,23)  ! L^\infty
           write(ife, 952) x(i,5), orders(i,1)  ! L^2
           write(ife, 952) x(i,6), orders(i,2)   !H^1
           write(ife, 952) x(i,33), orders(i,29)  ! E_I L^q or eta
           !write(ife, 902) x(i,16), orders(i,12)   !eta

        else
           !write(ife, 992) x(i,32)
           !write(ife, 992) x(i,31)
           write(ife, 992) x(i,27)
           write(ife, 992) x(i,5)
           write(ife, 992) x(i,6)
           write(ife, 992) x(i,33)
           !write(ife, 992) x(i,16) !eta
        endif

        !write(ife, 600) x(i, 31) / x(i,27) ! i_{eff}
        !write(ife, 600) x(i, 16) / x(i,7) ! i_{eff}
        !write(ife, 600) (x(i, 16)**2 + x(i,8)**2)**0.5 / (x(i,7)**2 + x(i,8)**2)**0.5 ! i_{eff}

        !write(ife, 301) x(i,24)  ! CPU time
        write(ife, 300) 

     endif

     if(typ == 'AMAst' ) then
        ! errors
        !if(i > 1) then
        !   write(ife, 952) x(i,27), orders(i,23)  ! L^\infty
        !   write(ife, 952) x(i,5), orders(i,1)  ! L^2
        !   write(ife, 952) x(i,6), orders(i,2)   !H^1
        !   write(ife, 952) x(i,9), orders(i,2)   !L2(L2)
        !   write(ife, 952) x(i,10), orders(i,2)   ! L2(H^1)
        !   write(ife, 952) x(i,32), orders(i,29)  ! E_I L^q or eta
        !   write(ife, 952) x(i,31), orders(i,29)  ! E_I L^q or eta

        !else
           write(ife, 201) x(i,27)
           write(ife, 201) x(i,5)
           write(ife, 201) x(i,6)
           write(ife, 201) x(i,9)
           write(ife, 201) x(i,10)
           write(ife, 201) x(i,32)
           write(ife, 201) x(i,31)
        !endif

        !write(ife, 600) x(i, 31) / x(i,27) ! i_{eff}
        !write(ife, 600) x(i, 16) / x(i,7) ! i_{eff}
        !write(ife, 600) (x(i, 16)**2 + x(i,8)**2)**0.5 / (x(i,7)**2 + x(i,8)**2)**0.5 ! i_{eff}

        write(ife, 301) x(i,24)  ! CPU time
        !write(ife, 300) 

     endif


     if(typ == 'CFD-AMA' ) then
        ! errors
        write(ife, 741) x(i,35)
        write(ife, 741) x(i,36)
        write(ife, 741) x(i,37)
        write(ife, 201) x(i,31)
        write(ife, 201) x(i,32)
        !write(ife, 992) x(i,33)
        !write(ife, 992) x(i,16) !eta
        
        write(ife, 301) x(i,24)  ! CPU time
        !!write(ife, 300) 
        
     endif


     if(typ == 'book-hp' ) then
        ! errors
        if(i > 1) then
           !write(ife, 902) x(i,32), orders(i,28)  ! E_I L^2
           !write(ife, 902) x(i,31), orders(i,27)  ! E_I L^\infty
           write(ife, 952) x(i,27), orders(i,23)  ! L^\infty
           write(ife, 952) x(i,5), orders(i,1)  ! L^2
           write(ife, 952) x(i,6), orders(i,2)   !H^1
           !write(ife, 952) x(i,33), orders(i,29)  ! E_I L^q or eta
           !write(ife, 902) x(i,16), orders(i,12)   !eta

        else
           !write(ife, 992) x(i,32)
           !write(ife, 992) x(i,31)
           write(ife, 992) x(i,27)
           write(ife, 992) x(i,5)
           write(ife, 992) x(i,6)
           !write(ife, 992) x(i,33)
           !write(ife, 992) x(i,16) !eta
        endif

        !write(ife, 600) x(i, 31) / x(i,27) ! i_{eff}
        !write(ife, 600) x(i, 16) / x(i,7) ! i_{eff}
        !write(ife, 600) (x(i, 16)**2 + x(i,8)**2)**0.5 / (x(i,7)**2 + x(i,8)**2)**0.5 ! i_{eff}

        write(ife, 301) x(i,24)  ! CPU time
        !write(ife, 300) 

     endif


     if(typ == 'book-IPG' ) then
        ! errors
        if(i > 1) then
           write(ife, 752) x(i,5), orders(i,1)  ! SIPG
           write(ife, 752) x(i,6), orders(i,2)  ! NIPG
           write(ife, 752) x(i,7), orders(i,3)  ! IIPG
        else
           write(ife, 792) x(i,5)  ! SIPG
           write(ife, 792) x(i,6)  ! NIPG
           write(ife, 792) x(i,7)  ! IIPG
        endif

        write(ife, 300) 

     endif

     if(typ == 'book-IPG2' ) then
        ! errors
        if(i > 1) then
           write(ife, 742) x(i,5), orders(i,1)  ! SIPG
           write(ife, 742) x(i,6), orders(i,2)  ! NIPG
           write(ife, 742) x(i,7), orders(i,3)  ! IIPG
           write(ife, 742) x(i,8), orders(i,4)  ! SIPG
           write(ife, 742) x(i,9), orders(i,5)  ! NIPG
           write(ife, 742) x(i,10), orders(i,6)  ! IIPG
        else
           write(ife, 782) x(i,5)  ! SIPG
           write(ife, 782) x(i,6)  ! NIPG
           write(ife, 782) x(i,7)  ! IIPG
           write(ife, 782) x(i,8)  ! SIPG
           write(ife, 782) x(i,9)  ! NIPG
           write(ife, 782) x(i,10)  ! IIPG
        endif

        write(ife, 300) 

     endif


     if(typ == 'cDLM' ) then
        ! errors
        write(ife, 266) x(i,35),x(i,38),x(i,36),x(i,39),x(i,37),x(i,40)

        write(ife, 301) x(i,24)  ! CPU time

     endif

     if(typ == 'RESbeam') then
        ! errors
        write(ife, 203) x(i,5),x(i,6),x(i,7)
        write(ife, 902) x(i,11), x(i,11) / x(i,7)  ! \eta, i_{eff}

        write(ife, 300) 
     endif

     write(ife_txt, *) x(i,indexi), x(i,4), 0., x(i,5), 0.  !!,x(i,8), 0.
     
     if(typ == 'RESo' .or. typ == 'RESo1'.or. typ == 'RESo2' ) then
        if( i == imesh )  write(ife, *) "\hline"
     elseif(typ == 'RESa' ) then
        if(i==1 .or. mod(i,iline) == 0 .or. i == imesh  ) write(ife, *) "\hline"

     else if(typ == 'book-IPG' ) then
         if( i == imesh ) then
            write(ife, *) ' \vspace{-3mm} & & & & & & & \\'
            !write(ife, *) '\bigskip'
            write(ife,* ) ' \multicolumn{2}{c|}{GEOC} '
            write(ife, 999)  orders(0,1:3) !, orders(0,5)
            write(ife,*) '\hline'
         endif

     else
        if(mod(i,iline) == 0 .or. i == imesh )  write(ife, *) "\hline"

        !if(mod(i,iline) == 0  .or. ipar(i,2) == 3)  write(ife, *) "\hline"
        !if(i == imesh )  write(ife_e, *) "\hline"
     endif
  enddo

  write(ife_txt, *) i, orders(0,1), orders(0,2) !, orders(0,5)


  write(ife_e, *) "\end{tabular}"
  write(ife_e, *) " } %%% end of font size  "
  
end program SetOrder

    subroutine Regres(n, x, y, order )
      integer, intent(in) :: n     ! number of datas
      real, dimension(1:n), intent(in) :: x, y   ! arrays of data
      real, intent(out) :: order 
      real :: sumx, sumy, sumx2, sumy2, sumxy
      real :: rb, ra, rr
      integer :: i

      sumx = sum( x(1:n) )
      sumy = sum( y(1:n) )
      sumx2 = dot_product( x(1:n), x(1:n) )
      sumy2 = dot_product( y(1:n), y(1:n) )
      sumxy = dot_product( x(1:n), y(1:n) )

      rb = (n*sumxy -sumx*sumy)/(n*sumx2-sumx*sumx)
      ra = (sumy - rb*sumx)/n
      rr = (n*sumxy-sumx*sumy)/sqrt((n*sumx2-sumx*sumx)*(n*sumy2-sumy*sumy))
      ra = exp(ra)
            
      order = rb

    end subroutine Regres
