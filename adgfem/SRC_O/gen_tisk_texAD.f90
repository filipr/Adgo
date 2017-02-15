      program gen_tisk

      character*27 Z_incl
      character*3 Z_end
      character*5 Z_ps, adapt_method, estim_method
      character*9 Z_plot, Z_plot1, IPG, chA
      character*20 type,model, method
      character*20 empty, type_tab
      character*6 time_method
      logical vertical_tab
      character*50 grid_file,ini_file, prof_file

      idat = 11
      itex = 12
      itex1 = 18
      igrid = 13
      iini = 14

      Z_incl = '   \includegraphics[width='
      Z_end = 'cm]'
      Z_ps = '.eps}'

      open(idat,file='num',status = 'OLD')
      read(idat,*) n_start, n_end
      read(idat,*) type
      read(idat,*) method
      read(idat,*) ini_file
      read(idat,*) type_tab

      num = n_end - n_start + 1
      close(idat)

      ! for long tables, we plot it vertically
      vertical_tab = .false.
      if(type_tab == 'AMA')  vertical_tab = .true.

      if(type_tab == 'pNeu')  vertical_tab = .true.
      if(type_tab == 'RES_STDGM')  vertical_tab = .true.

      print*,'##  gen_tisk_texAD.f90 settings: ',type, type_tab, vertical_tab, ini_file


      open(iini,file=ini_file, status = 'OLD')
      read(iini,*) model, Re, icase, par
      read(iini,*) Ttime
      read(iini,*) empty
      read(iini,*) nbDim, grid_file
      read(iini,*) iQ_k, prof_file
      read(iini,*) empty
      read(iini,*) empty
      read(iini,*) IPG, C_W, iP_k
      read(iini,*) time_method, iBDF
      read(iini,*) chA, BDF_tol
      read(iini,*) estim_method, tol_max, tol_min, rLq
      read(iini,*) adapt_method, adapt_iter
      read(iini,*) maxiter
      read(iini,*) empty           ! newton
      read(iini,*) empty           ! GMRES
      read(iini,*) out_time
      read(iini,*) ipoc

      if(model == 'NSe') then
         do i=1, ipoc
            read(iini,*) iw, iw1, ro, u, v, p
         enddo
      else
         do i=1, ipoc
            read(iini,*) iw, iw1, ro
         enddo
      endif

      read(iini,*) stab1, stab2, stab3, stab4

      read(iini,*) tol1, tol2, tol3, nu


      close(iini)

      open(igrid,file=grid_file, status = 'OLD')
      read(igrid,*) npoin, nelem, nbelm, nbc
      close(igrid)


      open(itex,file='tisk.tex',status = 'UNKNOWN')
      open(itex1,file='des.tex',status = 'UNKNOWN')

      write(itex,*) '\documentclass[10pt]{article}'
      write(itex,*) '\usepackage{epsfig}'
      write(itex,*) '\usepackage{rotating}'
      write(itex,*) '\usepackage{amssymb}'
      write(itex,*) '\usepackage{amsmath}'
      write(itex,*) '\usepackage{GGGraphics}'
      write(itex,*) '\hoffset=-41mm'
      write(itex,*) '\voffset=-8mm'
      write(itex,*) '\topmargin -2cm'
      write(itex,*) '%\leftmargin -25mm'
      write(itex,*) '\textheight 265mm'
      write(itex,*) '\textwidth 195mm'
      write(itex,*) '\def\Tr{{T_h}}'
      write(itex,*) '\def\Th{{T_h}}'
      write(itex,*) '\def\ehX{{\|e_h\|_X}}'
      write(itex,*) '\def\ehJ{{\|e_h\|_J}}'
      write(itex,*) '\def\hp{{hp}}'
      write(itex,*) '\def\etaX{{\eta}}'
      write(itex,*) '\def\dof{{\sf dof}}'
      write(itex,*) '\def\Nhp{{N_{hp}}}'
      write(itex,*) '\def\ehp{{e_{hp}}}'
      write(itex,*) '\def\ehp{{e_{hp}}}'

      write(itex,*) '\def\etaA{{\eta_\mathrm{A}}}'
      write(itex,*) '\def\etaS{{\eta_\mathrm{S}}}'
      write(itex,*) '\def\etaT{{\eta_\mathrm{T}}}'
      write(itex,*) '\def\etaST{{\eta_\mathrm{ST}}}'

      write(itex,*) '\def\Eh{{\mathcal{E}_h}}'
      write(itex,*) '\def\ieff{{i_{\mathrm{eff}}}}'

      write(itex,*) '\def\errX{{\|e_h\|_X}}'
      write(itex,*) '\def\errY{{\|e_h\|_Y}}'
      write(itex,*) '\def\Thp {{\cal T}_{hp}}'

      write(itex,*) '\begin{document}'
      write(itex,*) '\thispagestyle{empty}'
      write(itex,*) ' '
      write(itex,*) '%\leftmargin=+0.5in'
      write(itex,*) '%two figures one next to one'
      write(itex,*) '%\hspace*{56 mm} %insert the text'

!      write(itex,*) '\input {description.tex}'
      write(itex,*) '\input{des.tex}'
      write(itex,'(x)')

      write(itex1,*) '%\hspace{50mm}'
      if(type .eq. 'step') then
         write(itex1,*) '{\Large \bf Forward facing step -- ', &
              'Woodward-Collella  }'
      elseif(type .eq. 'gamm') then
         write(itex1,*) '{\Large \bf GAMM channel }'
      elseif(type .eq. 'battery') then
         write(itex1,*) '{\Large \bf battery }'
      elseif(type .eq. 'battery_sim') then
         write(itex1,*) '{\Large \bf simplified battery }'
      elseif(type .eq. 'naca') then
         write(itex1,*) '{\Large \bf NACA0012 profile }'
      elseif(type .eq. 'vortex') then
         write(itex1,*) '{\Large \bf Inviscid isentropic vortex propagation  }'
      elseif(type .eq. 'se1050') then
         write(itex1,*) '{\Large \bf SE 1050 profile cascade  }'
      elseif(type .eq. 'sod') then
         write(itex1,*) '{\Large \bf Sod''s tube }'
      elseif(type .eq. 'scalar') then
         write(itex1,*) '{\Large \bf scalar  }'
      elseif(type .eq. 'scalarST') then
         write(itex1,*) '{\Large \bf scalarST  }'
      elseif(type .eq. 'scalarL') then
         write(itex1,*) '{\Large \bf scalarL  }'
      elseif(type .eq. 'ARC') then
         write(itex1,*) '{\Large \bf ARC  }'
      elseif(type .eq. 'ARC2') then
         write(itex1,*) '{\Large \bf ARC2  }'
      elseif(type .eq. 'scalar11') then
         write(itex1,*) '{\Large \bf scalar $(-1,1)^2$ }'
      elseif(type .eq. 'singul') then
         write(itex1,*) '{\Large \bf singul  }'
      elseif(type .eq. 'ALG2') then
         write(itex1,*) '{\Large \bf Algebraic error estimates  }'
      elseif(type .eq. 'pulse') then
         write(itex1,*) '{\Large \bf Gaussian pulse  }'
      elseif(type .eq. 'scalar66') then
         write(itex1,*) '{\Large \bf  Barenblatt, porus media flow  }'
      elseif(type .eq. 'pNeu') then
         write(itex1,*) '{\Large \bf $p$-robust  }'
      elseif(type .eq. 'STDGM') then
         write(itex1,*) '{\Large \bf STDGM  }'
      else
         write(itex1,*) '{\Large \bf unknown problem }'
      endif
      write(itex1,*) ''
      write(itex1,*) '%\hspace{150mm}'
       write(itex1,*) '\today'

      write(itex1,*) ''
      write(itex1,*) '%\hspace{50mm}'

      if (time_method .eq. 'STDG') then
         write(itex1,'(a8,i6,a5,i1,a17,i2,a20)') &
              'nelem = ',nelem,', $P_',iP_k,'$ approximation, ', &
              iBDF,'- STDGM'

      else

         write(itex1,'(a8,i6,a5,i1,a17,i2,a20)') &
              'nelem = ',nelem,', $P_',iP_k,'$ approximation, ',&
              iBDF,'--step BDF Formula'
      endif

      if(out_time > 0) then
         write(itex1,'(a35,f7.4 )')', visualization after $\Delta t =$',&
              out_time
      endif


      write(itex1,*) ''

      write(itex1,*) '%\hspace{50mm}'
      if(model == 'NSe') then
         rmach = ((u*u + v*v)/(1.4*p/ro))**0.5
         alpha = asin( v/(v*v + u*u)**0.5)/asin(1.)*90

         write(itex1,'(a25,f10.2,a5)')  &
             '$M_{\mathrm{in}} = ', rmach,'$, '
         write(itex1,'(a25,f10.2,a5)')  &
             '$\alpha  = ', alpha,'$, '
         if(Re == 0.) then
            write(itex1, *) 'inviscid \\'
         else
            write(itex1,'(a25,f10.0,a5)')  &
                '$\mathrm{Re} = $', Re,'\\ '
         endif

      elseif(model == 'sca' .or. model=='scalar' ) then
         write(itex1,'(a24, i3, a2)')  &
             'scalar problem, icase = ',icase,': '

         if(icase == 2) then
            write(itex1,*) 'moving front, ', &
                '$u=\frac{1}{1+\mbox{exp}[(x_1+x_2-1-t)/\varepsilon]}$'
            write(itex1,'(a24,f10.6,a4)') &
                'with $\varepsilon=', 1./Re,'$,'
     		elseif(icase == 1) then
            write(itex1,*)'$-\nabla\cdot(K(u)\nabla u)+u \nabla u=f$,'
            write(itex1,*)  &
                'linear equation (only diffusion) $u = x(1-x)y(1-y)$'
                !'linear parabolic equation $u=sin(2 \Pi t x y)$'
         elseif(icase == 4) then
            write(itex1,*)'$-\nabla\cdot(K(u)\nabla u)+u \nabla u=f$,'
            write(itex1,*) 'singular corner $u=r^\alpha xy(1-x)(1-y)$'
            write(itex1,'(a14,f6.2,a4)') 'with $\alpha=', par,'$,'
         elseif(icase == 3) then
            write(itex1,*) 'linear boundary layer $\varepsilon=$'
            write(itex1,'(es10.2,a4)') 1./Re,',   '
         elseif(icase == 7) then
            write(itex1,*) '$\partial_t u -\Delta u = 0$, '
            write(itex1,*) '$\Omega = (0,1)^2$ '
            write(itex1,*) '$u=\exp(2t + x_1 + x_2)$, '
         elseif(icase == 14) then
            write(itex1,*) 'nonlinear difusion, regular solution '
            write(itex1,*) '[Houston,Sulli, Robson 07],'
         elseif(icase == 15) then
            write(itex1,*) 'singular corner L-shape domain'
            write(itex1,*) 'nonlinear difusion'
            write(itex1,*) '[Houston,Sulli, Robson 07]'
            write(itex1, *) '$u=r^{2/3} \sin(2\phi/3)$,'
         elseif(icase == 16) then
            write(itex1,*) 'John, Knobloch'
            write(itex1,*) 'linear difusion convection-diffusion, '
            write(itex1,'(a12es10.2,a4)') '$\varepsilon=',1./Re,'$,'
         elseif(icase == 24) then
            write(itex1,*) 'Barenblatt, porus media'
            write(itex1,*) '$K(u) = \kappa |u|^{\kappa-1}I$'
            write(itex1,'(a12es10.2,a4)') '$\kappa=',par,'$'
         elseif(icase == 28) then
            write(itex1,*) 'Two boundary layers:'
            write(itex1,*) '$-\Delta u = f$,'
            write(itex1,*) '$u=(1 - x_1)^{10}(1 - x_2)^{20}$,'
            write(itex1,*) 'linear difusion convection-diffusion, ARC '
            write(itex1,'(a12es10.2,a4)') '$\varepsilon=',1./Re,'$,'
         elseif(icase == 29) then
            write(itex1,*) 'Knopp, Lube,  Rapin CMAME 02'
            write(itex1,*) 'linear difusion convection-diffusion, ARC '
            write(itex1,'(a12es10.2,a4)') '$\varepsilon=',1./Re,'$,'
         elseif(icase == 30) then
            write(itex1,*) 'Generalized Knopp, Lube,  Rapin CMAME 02'
            write(itex1,*) 'linear difusion convection-diffusion, ARC2 '
            write(itex1,'(a14,es10.2,a4)') '$\varepsilon=',1./Re,'$,'
         elseif(icase == 35) then
            write(itex1,*) '$-\Delta u = f$,'
            write(itex1,*) 'singular corner $u=r^\alpha xy(1-x)(1-y)$'
            write(itex1,'(a14,f6.2,a4)') 'with $\alpha=', par,'$,'
         elseif(icase == 36) then
            write(itex1,*) '$-\Delta u = 0$, '
            write(itex1,*) 'L-shaped domain, '
            write(itex1,*) '$u=r^{3/2}\sin(2\varphi/3)$, '
         elseif(icase == 37) then
            write(itex1,*) '$-\Delta u = 0$, '
            write(itex1,*) '$\Omega=(0,1)^2$, '
            write(itex1,*) '$u=\sin(2\pi x_1) \sin(2\pi x_2)$'
         elseif(icase == 38) then
            write(itex1,*) 'Two boundary layers:'
            write(itex1,*) '$-\Delta u = f$,'
            write(itex1,*) '$u=(1 - x_1)^{10}(1 - x_2)^{10}$,'
            write(itex1,*) 'linear difusion convection-diffusion, ARC '
            write(itex1,'(a12es10.2,a4)') '$\varepsilon=',1./Re,'$,'
         elseif(icase == 39) then
            write(itex1,*) 'exponential boundary layers, '
            write(itex1,'(x)')
            write(itex1,*) '$-\Delta u = f$,'
            write(itex1,*) '$u =[c_1+c_2(1-x_1)+ \exp(-\alpha x_1)], ',  &
                '[c_1+c_2 (1 - x_2)  + \exp(-\alpha x_2)]$, '
            write(itex1,'(x)')
            write(itex1,*) '$c_1 = -\exp(-\alpha)$', &
                '$c_2=-1-c_1$, $\alpha='
            write(itex1,'(f8.2,a4)') par,'$,  '
         elseif(icase == 40) then
            write(itex1,*) 'LL-shaped domain, '
            write(itex1,*) '$-\Delta u = 0$, '
            write(itex1,*) '$u= r^{2/3} \sin(2\varphi/3)$, '
         elseif(icase == 42) then
            write(itex1,*) 'time dependent convection-diffusion, '
            write(itex1,*) '$u= 16 \frac{e^{mt}-1.}{e^m-1.}x_1 x_2(1-x_1) (1-x_2)$, '
            write(itex1,*) 'where m = 10'
         elseif(icase == 10) then
            write(itex1,*) 'singul, '
            write(itex1,*) '$-\Delta u = f$, '
            write(itex1,*) '$u= 2 r^\alpha x_1 x_2(1-x_1) (1-x_2)$, '
            write(itex1,'(a10,f8.2, a6)' ) '$\alpha = ' , par,' $, '
            write(itex1,'(a10,f7.1, a24)')'$u\in H^{', par+3,  &
                '+\epsilon}(\Omega) $'

         elseif(icase == 44) then
            write(itex1,*) 'reaction - convection equation'
            write(itex1,*) '$b\cdot \nabla u + \sigma u = f$, '
            write(itex1,'(a30,f8.1,a4)') &
                '$b=(1,1),\ \sigma = $ ',par,',  '

            !write(itex1,*)'$u=\exp(-100[(x_1-t-1/4)^2 +(x_2-t-1/4)^2])$'
            write(itex1,*) &
                '$u={0.1+\exp(10 t)} xy(1-x)(1-y)$'
!     *           '$u=\frac{0.1+\exp(10 t)}{0.1 + \exp(10)}xy(1-x)(1-y)$'

         elseif(icase == 45) then
            write(itex1,*) 'convection(-reaction) equation'
            write(itex1,*) '$\partial_t  u+b\cdot \nabla u', &
                '+\sigma u = f$, '
            write(itex1,'(a30,f8.1,a4)') &
                '$b=(1,1),\ \sigma = $ ',par,',  '

            write(itex1,*)'$u=\exp(-100[(x_1-t-1/4)^2 +(x_2-t-1/4)^2])$'
        elseif(icase == 48) then
            write(itex1,*) 'convection diffusion equation'
            write(itex1,*) '$\partial_t  u+ \partial_{x_1} u + \partial_{x_2} u ', &
                '-\varepsilon \Delta u = f$, '
            write(itex1,'(a30,f8.1,a4)') &
                '$\varepsilon = $ ', 1./Re,',  '

            write(itex1,*)'$u=(1+t) x_1(x_1 - 1) x_2 (x_2 -1) $'
        elseif(icase == 49) then
           write(itex1,*) 'convection diffusion equation'
           write(itex1,*) '$\partial_t  u+ \partial_{x_1} u + \partial_{x_2} u ', &
                '-\varepsilon \Delta u = f$, '
           write(itex1,'(a30,f8.1,a4)') &
                '$\varepsilon = $ ', 1./Re,',  '

           write(itex1,*)'$u=(0.1+\exp(10 t)) x_1(x_1 - 1) x_2 (x_2 -1) $'

        elseif(icase == 50) then
           write(itex1,*) '$-\Delta u = 0$, '
           write(itex1,*) '$\Omega=(0,1)^2$, '
           write(itex1,*) '$u=\sin( \pi x_1) \sin( \pi x_2)$'

        elseif(icase == 53) then
            write(itex1,*) '$-\Delta u = f$, '
            write(itex1,*) '$u= \exp(x_1 x_2)$ '
        elseif(icase == 64) then
            write(itex1,*)'$-\nabla\cdot(\nabla u)=f$,'
            write(itex1,*)  &
                'point value test case in a=(0.25,0.25) $u = -1/2 log( |x-a|^2 + \varepsilon^2 ) $'
        elseif(icase == 69) then
            write(itex1,*)'$\partial_t \rho- \nabla\cdot( K(|\nabla\rho|)\nabla \rho)= f$,'
            write(itex1,*)'$K(\xi) = \frac{2}{1+\sqrt{1+4\xi}}$,'
            write(itex1,'(x)')
            write(itex1,*)'$\rho(x,t) =\exp(-2t)x_1(1-x_1)x_2(1-x_2)$,'
         else
            write(itex1,*) '!!unknown problem '
            write(itex1,*) '(add in gen tisk texAD.f in SRC)'
        endif
        write(itex1,'(x)')
      endif

      if(adapt_iter >0) then
         write(itex1,*) 'adaptation method:  ', adapt_method,', '
         write(itex1,'(a10,es12.4)') ' tol = ', tol_max
      endif

      write(itex1,*) '\vskip 2mm'

      if(type /= 'step' .and. type /= 'pulse'  &
          .and. type /= 'ALG2'  &
          .and. type /= 'shock' .and. type /= 'gamm' ) then

         if(vertical_tab) write(itex,*) '\begin{sidewaystable}[h]'
         write(itex,*) '%\input{des.tex}'
         write(itex,*) '\input{tabA.tex}'
         if(vertical_tab) write(itex,*) '\end{sidewaystable}'
         write(itex,*) '\vskip 5mm'
      endif

      num_fig_b = 1
      num_fig_e = 5

      if(type .eq. 'step') then
         !num_fig_b = 1
         !num_fig_e = 2
         num_fig_b = 3
         num_fig_e = 3
      endif

      if(type .eq. 'pulse') then
         !num_fig_b = 1
         !num_fig_e = 2
         num_fig_b = 3
         num_fig_e = 4
      endif

      if(type .eq. 'gamm') then
         num_fig_b = 1
         num_fig_e = 3
      endif

      if(type .eq. 'sod') then
         num_fig_b = 1
         num_fig_e = 3
      endif

      if(type .eq. 'naca') then
         num_fig_b = -1
         num_fig_e = 6
      endif

      if(type .eq. 'battery' ) then
         num_fig_b = 1
         num_fig_e = 4
      endif

      if(type .eq. 'battery_sim') then
         num_fig_b = 1
         num_fig_e = 3 ! 4
      endif

      if(type .eq. 'pNeu') then
         num_fig_b = 0
         num_fig_e = 7
      endif

      if(type .eq. 'scalar' .or. type .eq. 'singul' &
         .or.  type .eq. 'scalarL'.or.  type .eq. 'BL'  &
          .or.  type .eq. 'scalarBL'.or.  type .eq. 'ARC' &
          .or.  type .eq. 'ARC' .or. type .eq. 'ARC2' ) then
         num_fig_b = 0
         num_fig_e = 7
      endif

      if(type .eq. 'JK') then
         num_fig_b = 0
         num_fig_e = 5
      endif

      if(type .eq. 'scalar11' ) then
         num_fig_b = 0
         num_fig_e = 6
      endif

      if(type .eq. 'scalarST' ) then
      	write(itex1,*) 'scalarST in \verb|gen_tisk|'
         num_fig_b = 0
         num_fig_e = 3
      endif

      if(type .eq. 'scalar66' ) then
         num_fig_b = 0
         num_fig_e = 3
      endif


      if(type .eq. 'shock' ) then
         num_fig_b = 1
         num_fig_e = 4
      endif

      if(type .eq. 'vortex' ) then
         num_fig_b = 0
         num_fig_e = 2
      endif

      if(type .eq. 'se1050' ) then
         num_fig_b = 0
         num_fig_e = 3
      endif

      if(method .eq. 'meshes' ) then
         num_fig_b = 1
         num_fig_e = 5
      endif

      !different type of drawing
      if(type .eq. 'ALG2') then

         write(itex,*) "\includegraphics[width= 8.00cm] {estims-L.eps}"
         write(itex,*) "\includegraphics[width= 8.00cm] {errs-L.eps}"
         write(itex,'(x)')
         write(itex,*) '\vspace{10mm}'

         num_fig_b = 1
         num_fig_e = 0

         width = 8.

         do k=  n_start, n_end
            write(itex,'(x)')
            write(itex,*) '{\Large \bf Mesh adaptation level = ', k,'}'
            write(itex,'(x)')

            do i=1, 3
               if(i==1) then
                  Z_plot ='{err_Tot-'
                  Z_plot1='{est_Tot-'
               elseif(i==2) then
                  Z_plot ='{err_Dis-'
                  Z_plot1='{est_Dis-'
               elseif(i==3) then
                  Z_plot ='{err_Alg-'
                  Z_plot1='{est_Alg-'
               endif

               call WRITE_LINE(itex,Z_incl,Z_end,Z_plot, Z_ps,k, width)
               write(itex,*) '\hspace{5mm}'
               call WRITE_LINE(itex,Z_incl,Z_end,Z_plot1,Z_ps,k, width)
               write(itex,'(x)')
            enddo

            write(itex,*) '\newpage'
            write(itex,'(x)')
         enddo


      endif

      do ifig = num_fig_b, num_fig_e

         if(method == 'meshes') then
            ncol = 3
            width = 6.
            if(ifig .eq. 1) then
               Z_plot='  {meshA-'

            elseif(ifig .eq. 2) then
               Z_plot='  {meshB-'

            elseif(ifig .eq. 3) then
               Z_plot='  {meshC-'

            elseif(ifig .eq. 4) then
               Z_plot='  {meshD-'

            elseif(ifig .eq. 5) then
               Z_plot='  {meshE-'
            endif

         else

            if(type .eq. 'scalar' .or. type .eq. 'scalarL'  &
                 .or. type .eq. 'singul'.or. type .eq. 'BL' &
                 .or. type .eq. 'scalarBL' .or. type .eq. 'ARC'  &
                 .or. type .eq. 'ARC2' .or. type .eq. 'pNeu' ) then
               if(ifig .eq. 0) then
                  ncol = 3
                  width = 6.
                  !               Z_plot='   {mesh-'
                  Z_plot=' {hpmesh-'
               elseif(ifig .eq. 1) then
                  ncol = 3
                  width = 6.
                  !               Z_plot='   {mesh-'
                  Z_plot='{hpmeshD-'
               elseif(ifig .eq. 2) then
                  ncol = 3
                  width = 6.
                  Z_plot='{hpmeshE-'
               elseif(ifig .eq. 3) then
                  ncol = 3
                  width = 6.
                  Z_plot='{hpmeshF-'
               elseif(ifig .eq. 4) then
                  ncol = 3
                  width = 6.
                  Z_plot='{hpmeshG-'
               elseif(ifig .eq. 5) then
                  ncol = 3
                  width = 6.
                  !Z_plot='  {U_ISO-'
                  Z_plot='{hpmeshH-'
               elseif(ifig .eq. 6) then
                  ncol = 3
                  width = 6.
                  !Z_plot='  {U_CUT-'
                  Z_plot='{hpmeshI-'
               elseif(ifig .eq. 7) then
                  ncol = 3
                  width = 6.
                  !Z_plot=' {U_CUTa-'
                  Z_plot='{hpmeshJ-'
               endif

               if(type .eq. 'ARC2') then
                  ncol = 2
                  width = 9.
               endif

            elseif( type .eq.'JK' ) then
               if(ifig .eq. 0) then
                  ncol = 3
                  width = 6.
                  !               Z_plot='   {mesh-'
                  Z_plot=' {hpmesh-'
               elseif(ifig .eq. 1) then
                  ncol = 3
                  width = 6.
                  !               Z_plot='   {mesh-'
                  Z_plot='{hpmeshD-'
               elseif(ifig .eq. 2) then
                  ncol = 3
                  width = 6.
                  !               Z_plot='   {mesh-'
                  Z_plot='{hpmeshE-'
               elseif(ifig .eq. 3) then
                  ncol = 3
                  width = 6.
                  Z_plot='  {U_ISO-'
               elseif(ifig .eq. 4) then
                  ncol = 3
                  width = 6.
                  Z_plot='  {U_CUT-'
               elseif(ifig .eq. 5) then
                  ncol = 3
                  width = 6.
                  Z_plot=' {U_CUTa-'
               endif


            elseif(type .eq. 'scalar11' ) then
               if(ifig .eq. 0) then
                  ncol = 3
                  width = 6.
                  !               Z_plot='   {mesh-'
                  Z_plot=' {hpmesh-'
               elseif(ifig .eq. 1) then
                  ncol = 3
                  width = 6.
                  !               Z_plot='   {mesh-'
                  Z_plot='{hpmeshD-'

               elseif(ifig .eq. 2) then
                  ncol = 3
                  width = 6.
                  Z_plot='{hpmeshS-'

               elseif(ifig .eq. 3) then
                  ncol = 3
                  width = 6.
                  Z_plot='{hpmeshE-'

               elseif(ifig .eq. 4) then
                  ncol = 3
                  width = 6.
                  Z_plot='{hpmeshF-'

               elseif(ifig .eq. 5) then
                  ncol = 3
                  width = 6.
                  Z_plot='{hpmeshG-'

               elseif(ifig .eq. 6) then
                  ncol = 3
                  width = 6.
                  Z_plot='{hpmeshH-'

               elseif(ifig .eq. 7) then
                  ncol = 3
                  width = 6.
                  Z_plot='{hpmeshI-'

               elseif(ifig .eq. 8) then
                  ncol = 3
                  width = 6.
                  Z_plot='{hpmeshJ-'

               !elseif(ifig .eq. 4) then
               !   ncol = 3
               !   width = 6.
               !   Z_plot='  {U_ISO-'
               !elseif(ifig .eq. 5) then
               !   ncol = 3
               !   width = 6.
               !   Z_plot='  {U_CUT-'
                  !elseif(ifig .eq. 4) then
                  !   ncol = 3
                  !   width = 6.
                  !   Z_plot=' {U_CUTa-'
               endif

            elseif(type .eq. 'scalar66' ) then
               if(ifig .eq. 0) then
                  ncol = 3
                  width = 6.
                  !               Z_plot='   {mesh-'
                  Z_plot=' {hpmesh-'
               elseif(ifig .eq. 1) then
                  ncol = 3
                  width = 6.
                  !               Z_plot='   {mesh-'
                  Z_plot='{hpmeshD-'

               ! elseif(ifig .eq. 2) then
               !    ncol = 3
               !    width = 6.
               !    Z_plot='{hpmeshS-'

               ! elseif(ifig .eq. 3) then
               !    ncol = 3
               !    width = 6.
               !    Z_plot='{hpmeshE-'

               ! elseif(ifig .eq. 4) then
               !    ncol = 3
               !    width = 6.
               !    Z_plot='{hpmeshF-'

               ! elseif(ifig .eq. 5) then
               !    ncol = 3
               !    width = 6.
               !    Z_plot='{hpmeshG-'

               ! elseif(ifig .eq. 6) then
               !    ncol = 3
               !    width = 6.
               !    Z_plot='{hpmeshH-'

               ! elseif(ifig .eq. 7) then
               !    ncol = 3
               !    width = 6.
               !    Z_plot='{hpmeshI-'

               ! elseif(ifig .eq. 8) then
               !    ncol = 3
               !    width = 6.
               !    Z_plot='{hpmeshJ-'

                  elseif(ifig .eq. 2) then
                  ncol = 3
                  width = 6.
                  Z_plot='  {U_ISO-'
               elseif(ifig .eq. 3) then
                  ncol = 3
                  width = 6.
                  Z_plot='  {U_CUT-'
                  !elseif(ifig .eq. 4) then
                  !   ncol = 3
                  !   width = 6.
                  !   Z_plot=' {U_CUTa-'
               endif

            elseif(type .eq. 'scalarST') then
               if(ifig .eq. 0) then
                  ncol = 2
                  width = 9.
                  Z_plot='   {mesh-'
                  !               Z_plot='{hpmesh-'
               elseif(ifig .eq. 1) then
                  ncol = 2
                  width = 9.
                  Z_plot='  {U_ISO-'
                  !               Z_plot='{hpmeshD-'
               elseif(ifig .eq. 2) then
                  ncol = 3
                  width = 9.
                  Z_plot='   {U_3D-'
               elseif(ifig .eq. 3) then
                  ncol = 3
                  width = 9.
                  Z_plot='  {U_CUT-'
               endif


            elseif(type .eq. 'shock' ) then
               if(ifig .eq. 1) then
                  ncol = 2
                  width = 9.
                  Z_plot='  {P_ISO-'
               elseif(ifig .eq. 2) then
                  ncol = 3
                  width = 6.
                  Z_plot='  {P_CUT-'
               elseif(ifig .eq. 3) then
                  ncol = 2
                  width = 9.
                  !               Z_plot='   {mesh-'
                  Z_plot=' {hpmesh-'
               elseif(ifig .eq. 4) then
                  ncol = 2
                  width = 9.
                  !               Z_plot='   {mesh-'
                  Z_plot='{hpmeshS-'
               endif

            elseif(type .eq. 'step') then
               if(ifig .eq. 1) then
                  ncol = 1
                  width = 18.
                  !               Z_plot='   {mesh-'
                  Z_plot=' {hpmesh-'
               elseif(ifig .eq. 2) then
                  ncol = 1
                  width = 18.
                  !               Z_plot='   {mesh-'
                  Z_plot='{hpmeshD-'
               elseif(ifig .eq. 3) then
                  ncol = 2
                  width = 9.
                  !               Z_plot='   {mesh-'
                  Z_plot=' {RO_ISO-'
               endif

            elseif(type .eq. 'pulse') then
               if(ifig .eq. 1) then
                  ncol = 2
                  width = 9.
                  !               Z_plot='   {mesh-'
                  Z_plot=' {hpmesh-'
               elseif(ifig .eq. 2) then
                  ncol = 2
                  width = 9.
                  !               Z_plot='   {mesh-'
                  Z_plot='{hpmeshD-'
               elseif(ifig .eq. 3) then
                  ncol = 3
                  width = 6.
                  !               Z_plot='   {mesh-'
                  Z_plot='  {P_ISO-'
               elseif(ifig .eq. 4) then
                  ncol = 3
                  width = 6.
                  !               Z_plot='   {mesh-'
                  Z_plot='  {M_ISO-'
               endif

            elseif(type .eq. 'vortex') then
               if(ifig .eq. 0) then
                  ncol = 3
                  width = 6.
                  Z_plot='   {mesh-'
               elseif(ifig .eq. 1) then
                  ncol = 3
                  width = 6.
                  Z_plot='  {M_ISO-'
               elseif(ifig .eq. 2) then
                  ncol = 3
                  width = 6.
                  Z_plot='  {M_CUT-'
               endif

            elseif(type .eq. 'se1050') then
               if(ifig .eq. 0) then
                  ncol = 3
                  width = 6.
                  Z_plot='   {mesh-'
               elseif(ifig .eq. 1) then
                  ncol = 3
                  width = 6.
                  Z_plot='  {M_ISO-'
               elseif(ifig .eq. 2) then
                  ncol = 3
                  width = 6.
                  Z_plot=' {M_ISOa-'
               elseif(ifig .eq. 3) then
                  ncol = 3
                  width = 6.
                  Z_plot='  {P_CUT-'
               endif

            elseif(type .eq. 'naca') then
               if(ifig .eq. -1) then
                  ncol = 2
                  width = 9.
                  Z_plot='{hpmeshT-'
               elseif(ifig .eq. 0) then
                  ncol = 2
                  width = 9.
                  Z_plot=' {hpmesh-'
               elseif(ifig .eq. 1) then
                  ncol = 2
                  width = 9.
                  Z_plot='{hpmeshD-'
               elseif(ifig .eq. 2) then
                  ncol = 2
                  width = 9.
                  Z_plot='{hpmeshE-'
               elseif(ifig .eq. 3) then
                  ncol = 2
                  width = 9.
                  Z_plot='  {M_ISO-'
               elseif(ifig .eq. 4) then
                  ncol = 3
                  width = 6.
                  Z_plot='{P_WALLS-'
               elseif(ifig .eq. 5) then
                  ncol = 2
                  width = 9.
                  Z_plot='{hpmeshF-'
               elseif(ifig .eq. 6) then
                  ncol = 2
                  width = 9.
                  Z_plot='{hpmeshG-'
               endif

            elseif(type .eq. 'battery') then
               if(ifig .eq. 1) then
                  ncol = 4
                  width = 4.
                  Z_plot=' {hpmesh-'
               elseif(ifig .eq. 2) then
                  ncol = 4
                  width = 4.
                  Z_plot='  {U_sol-'
               elseif(ifig .eq. 3) then
                  ncol = 4
                  width = 4.
                  Z_plot='  {U_ISO-'
               elseif(ifig .eq. 4) then
                  ncol = 4
                  width = 4.
                  Z_plot='{hpmeshD-'
               endif

            elseif(type .eq. 'battery_sim') then
               if(ifig .eq. 1) then
                  ncol = 2
                  width = 9.
                  Z_plot=' {hpmesh-'
               elseif(ifig .eq. 2) then
                  ncol = 2
                  width = 9.
                  Z_plot='{hpmeshE-'
               !elseif(ifig .eq. 3) then
               !   ncol = 4
               !   width = 4.
               !   Z_plot='  {U_ISO-'
               elseif(ifig .eq. 3) then
                  ncol = 2
                  width = 9.
                  Z_plot='{hpmeshF-'

               elseif(ifig .eq. 4) then
                  ncol = 3
                  width = 6.
                  Z_plot='{hpmeshG-'
               endif

            elseif(type .eq. 'gamm') then
               if(ifig .eq. 1) then
                  ncol = 2
                  width = 9.
                  Z_plot=' {hpmesh-'
                  !elseif(ifig .eq. 2) then
                  !   ncol = 2
                  !   width = 9.
                  !   Z_plot='{hpmeshE-'
                  !elseif(ifig .eq. 3) then
                  !   ncol = 1
                  !   width = 14.
                  !   Z_plot='  {M_ISO-'
               elseif(ifig .eq. 2) then
                  ncol = 2
                  width = 9.
                  Z_plot='  {M_ISO-'
               elseif(ifig .eq. 3) then
                  ncol = 3
                  width = 6.
                  Z_plot='{P_WALLS-'
               endif

            else
               if(ifig .eq. 1) then
                  ncol = 6
                  width = 3.
                  Z_plot='{M_WALLS-'
                  if(type .eq. 'sod') Z_plot='{R_WALLS-'
               elseif(ifig .eq. 2) then
                  ncol = 6
                  width = 3.
                  Z_plot='{P_WALLS-'
               elseif(ifig .eq. 3) then
                  ncol = 3
                  width = 6.
                  Z_plot='  {M_ISO-'
                  if(type .eq. 'step') then
                     ncol = 2
                     width = 9.
                  endif
                  if(type .eq. 'sod') then
                     ncol = 6
                     width = 3.
                     Z_plot='{V_WALLS-'
                  endif
               elseif(ifig .eq. 4) then
                  if(type .eq. 'step') then
                     ncol = 2
                     width = 9.
                     Z_plot='   {M_3D-'
                  else
                     goto 100
                  endif

               elseif(ifig .eq. 5) then
                  if(type .ne. 'step') then
                     ncol = 3
                     width = 6.
                     Z_plot='  {P_ISO-'
                  else
                     ncol =2
                     width = 9.
                     Z_plot=' {RO_ISO-'
                  endif

               endif
            endif

         endif

         write(itex,'(x)')
         write(itex,*) '\begin{center}'
         write(itex,*) '{\Large \bf \verb|', Z_plot,'x}, x=', &
             n_start,',...,',n_end,' |}'
         write(itex,'(x)')
         write(itex,*) '\vspace{2mm}'
         write(itex,'(x)')
         write(itex,*) '\end{center}'
         write(itex,'(x)')


!         print*,'####',ifig, num, ncol, num/ncol,
!     *        n_end -  (n_end - mod( num, ncol)+1 )

         do i = 1, num/ncol
            write(itex,*) '\begin{center}'
            write(itex,*) '{\tiny ',(i-1)*ncol + n_start,'\dots', &
                i*ncol - 1 + n_start,' }'
            do j = 1, ncol
               k = (i-1)*ncol + j - 1 + n_start
               call WRITE_LINE(itex,Z_incl,Z_end,Z_plot,Z_ps,k, width)
            enddo
            write(itex,*) '\end{center}'
         enddo

         if(n_end -  (n_end - mod( num, ncol)+1) .ge. 0) then
            write(itex,*) '\begin{center}'
            write(itex,*) '{\tiny ',n_end - mod( num, ncol)+1,'\dots', &
                n_end,' }'
            do k = n_end - mod( num, ncol)+1, n_end
               call WRITE_LINE(itex,Z_incl,Z_end,Z_plot,Z_ps,k, width)
            enddo
            write(itex,*) '\end{center}'
         endif

         if(ifig .lt. num_fig_e)  write(itex,*) '\newpage'
 100     continue
      enddo

      if(method /= 'meshes' ) then
         if(type .ne. 'shock' .and. type .ne. 'gamm'  &
              .and. type .ne. 'ALG2' .and. type .ne. 'vortex'.and. type .ne. 'se1050') then
            write(itex,*) "\includegraphics[width= 14.00cm] {errs1.eps}"
            write(itex,'(x)')

            if(type .eq. 'scalar' .or. type .eq. 'scalarL'  &
                 .or. type .eq. 'scalar11' .or. type .eq. 'BL'  &
                 .or. type .eq. 'scalar66'  &
                 .or. type .eq. 'singul' .or. type .eq. 'JK'  &
                 .or. type .eq. 'scalarBL'.or. type .eq. 'battery' &
                 .or. type .eq. 'battery_sim' .or. type .eq. 'ARC' ) then

               write(itex,*) "\includegraphics[width= 14.00cm] {errs2.eps}"
               write(itex,'(x)')
               write(itex,*) "\includegraphics[width= 14.00cm] {errs3.eps}"

            elseif(type .eq. 'pNeu' ) then
               !write(itex,*) "\includegraphics[width= 14.00cm] {errs2.eps}"

               !write(itex,'(x)')

            endif
         endif
      endif

      write(itex,*) '\end{document}'
      close(itex)

    end program gen_tisk

      subroutine WRITE_LINE(itex,Z_incl,Z_end,Z_plot,Z_ps,k, width)
      character*27 Z_incl
      character*3 Z_end
      character*5 Z_ps
      character*5 Z_num, ch5
      character*9 Z_plot

      real width
      integer itex, k
      integer text_size, num_size, inum, is

      text_size = 0
      num_size = 5

      inum = k
      if(k < 0) inum = 10**num_size + k

      Z_num = '00000'


      if(inum > 0) then
         is = int(log(1.*inum)/log(10.))
      else
         is = 0
      endif


      write( ch5, '(i5)' ) inum ! change the format if num_size /= 5 !!!
      Z_num(num_size+text_size-is:num_size+text_size)  &
          = ch5(num_size-is: num_size)



      write(itex,'(a27,f8.2,a3,a9,a5,a5)')  &
             Z_incl,width,Z_end,Z_plot,Z_num,Z_ps

!      if(k<= 9) then
!         write(itex,'(a27,f8.2,a3,a9,i1,a5)')
!     *        Z_incl,width,Z_end,Z_plot,k,Z_ps
!
!      elseif(k<= 99) then
!         write(itex,'(a27,f8.2,a3,a9,i2,a5)')
!     *        Z_incl,width,Z_end,Z_plot,k,Z_ps
!
!      elseif(k<= 999) then
!         write(itex,'(a27,f8.2,a3,a9,i3,a5)')
!     *        Z_incl,width,Z_end,Z_plot,k,Z_ps

!      else
!         print*,'$#$%##$$  ',k
!         stop
!      endif

    end subroutine WRITE_LINE
