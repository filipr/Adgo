      program gen_tisk

      character*27 Z_incl
      character*3 Z_end
      character*5 Z_ps
      character*9 Z_plot
      character*10 type
      character*20 empty
      character*1 time_method

      character*50 grid_file,ini_file
      
      iini = 14

      open(iini,file=ini_file, status = 'OLD')
      read(iini,*) nbDim
      read(iini,*) grid_file
      read(iini,*) iQ_k
      read(iini,*) empty
      read(iini,*) empty
      read(iini,*) ndim, icase, par
      read(iini,*) iP_k
      read(iini,*) time_method, iBDF
      read(iini,*) empty
      read(iini,*) empty
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

      write(icase)
      
    end program gen_tisk
