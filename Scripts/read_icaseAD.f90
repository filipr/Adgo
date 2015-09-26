      program gen_tisk

      character*27 Z_incl
      character*3 Z_end
      character*5 Z_ps
      character*9 Z_plot
      character*10 type
      character*20 empty
      character*1 time_method

      character*50 grid_file,ini_file, prof_file, chtype
      character*50 chA, chB, chC, model, IPG
      integer ::  itype
      iini = 14

      !if (command_argument_count()== 2) then
      !   call get_command_argument(1, ini_file)
      !   call get_command_argument(2, chtype)
      !   input = 10
      !   open(iini,file=ini_file, status = 'OLD')
      !else
      !   print*,' Syntax error:  iread_icase  <file.ini>  < itype> ' 
      !   stop
      !endif
      !if(chtype == "1") itype = 1
      !if(chtype == "2") itype = 2
      !if(chtype == "3") itype = 3

      !read(*,*) ini_file
      !open(iini,file=ini_file, status = 'OLD')

      read(*,*) itype


      open(iini,file='smaz.ini', status = 'OLD')

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
      read(iini,*) empty
      read(iini,*) empty
      read(iini,*) maxiter         
      read(iini,*) empty           ! newton
      read(iini,*) empty           ! GMRES
      read(iini,*) out_time
      read(iini,*) nBC
      do i=1,nBC
         read(iini,*) empty
      enddo
      read(iini,*) empty
      read(iini,*) tol1, tol2, tol3, nu

      close(iini)

      if(itype == 1) write(*,*) icase
      if(itype == 2) write(*,*) par
      if(itype == 3) write(*,*) Re
      if(itype == 4) write(*,'(f10.2)') tol1
      if(itype == 5) write(*,'(f10.2)') tol2
      if(itype == 6) write(*,'(f10.2)') tol3
      if(itype == 7) write(*,'(i5)') nu
      
    end program gen_tisk
