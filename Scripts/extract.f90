!> main program
program extract

  character(len=50) :: command_name, tri_name, sol_name, exa_name, err_name, est_name
  character(len=1) :: ch1
  character(len=3) :: ch3
!  character(len=5) :: ch5
  character(len=10) :: ch10
  character(len=20) :: chA
  integer :: i, i0, isca
  !> real time instants (CPU)
  real, dimension(:,:), allocatable :: xi
  integer:: input, output, n
  character(len=128):: input_file, output_file
  logical :: inside

  if (command_argument_count() == 6) then
    call get_command_argument(1,input_file)
    input = 10
    open(input,file=input_file,action='read')

    call get_command_argument(2,output_file)
    output = 11
    open(output,file=output_file, status='replace')

    call get_command_argument(3,chA)
    read(chA, *) x1

    call get_command_argument(4,chA)
    read(chA, *) x2

    call get_command_argument(5,chA)
    read(chA, *) y1

    call get_command_argument(6,chA)
    read(chA, *) y2


  else
     print*,'Syntax extract <in_file>   <out_file>   <x1> <x2> < y1> <y2>'
  endif

  close(output)

  n = 100
  allocate(xi(0:n, 1:2) )

10 continue

  inside = .true.
  do i=0,n
     read(input, *, end=110) xi(i,1:2)
     
     if( xi(i, 1) < x1 .or. xi(i, 1) > x2 .or.  xi(i, 2) < y1 .or. xi(i, 2) > y2 ) inside = .false.
  enddo
  read(input,*) chA


  if(inside) then
     open(output,file=output_file, status='OLD', position = 'append')

     
     do i=0,n
        write(output, *) xi(i,1:2)
     enddo
     write(output,*) "##"
     write(output,*) "   "

     close(output)

  endif


  !print*,'##', chA
  goto 10

110 continue


end program


