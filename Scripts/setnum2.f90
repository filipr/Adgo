program SetNumber 
  integer:: inum
  character(len=1) T
  character(len=30) :: file_name
  character(len=5) :: ch5
  integer :: num_size, text_size
  integer :: is
  
  !read(*,*) inumber, T
  read(*,*) inumber
  
  !text_size = 4
  text_size = 0
  num_size = 2
  
  inum = inumber
  if(inumber < 0) inum = 10**num_size + inumber
  
  !print*,'????',inum, inumber
  
  !if(T == 'T') then
  !   file_name = 'triA00000    '
  !elseif(T == 'S') then
  !   file_name = 'solA00000    '
  !else
  !   print*,'Bad input value "T" in SetCommandWithFileName in problem.f90'
  !   stop
  !endif

  file_name = '00       '

  
  if(inum > 0) then
     is = int(log(1.*inum)/log(10.)) 
  else
     is = 0
  endif
  
  !print*,'!!!',inum,is, num_size+text_size-is, num_size+text_size, num_size-is, num_size
  
  write( ch5, '(i2)' ) inum  ! change the format if num_size /= 5 !!!
  file_name(num_size+text_size-is:num_size+text_size) = ch5(num_size-is: num_size)
  
  
  print*,file_name

end program SetNumber
