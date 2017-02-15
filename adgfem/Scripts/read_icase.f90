      program gen_tisk

      character*27 Z_incl
      character*3 Z_end
      character*5 Z_ps
      character*9 Z_plot
      character*10 type
      character*20 empty
      character*1 time_method

      character*50 grid_file,ini_file, chtype
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

      open(iini,file='smaz.ini', status = 'OLD')
      read(*,*) itype

      read(iini,*) nbDim
      read(iini,*) grid_file
      read(iini,*) iQ_k
      read(iini,*) empty
      read(iini,*) empty
      read(iini,*) ndim, icase, par
      read(iini,*) iP_k
      read(iini,*) time_method, iBDF
      read(iini,*) empty
      read(iini,*) tol1, tol2, tol3, nu
      read(iini,*) empty
      read(iini,*) init
      read(iini,*) empty
      read(iini,*) empty
      read(iini,*) maxiter
      read(iini,*) empty
      read(iini,*) tol
      read(iini,*) empty
      read(iini,*) empty
      read(iini,*) empty
      read(iini,*) Ttime
      read(iini,*) empty
      read(iini,*) empty
      read(iini,*) out_time
      read(iini,*) empty
      read(iini,*) CFL
      read(iini,*) BDF_tol
      read(iini,*) Re
      close(iini)

      if(itype == 1) write(*,*) icase
      if(itype == 2) write(*,*) par
      if(itype == 3) write(*,*) Re
      if(itype == 4) write(*,'(f10.2)') tol1
      if(itype == 5) write(*,'(f10.2)') tol2
      if(itype == 6) write(*,'(f10.2)') tol3
      if(itype == 7) write(*,'(i5)') nu
      
    end program gen_tisk
