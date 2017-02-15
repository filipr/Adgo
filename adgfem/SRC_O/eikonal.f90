!> modul contains setings for pedestrians model important fo BR algorithmh
module PedestrianSettings
  ! real :: po=1 ! is not used in BR algorithm
  ! real :: gamma=2 ! is not used BR algorithm
  real, parameter,private :: vmax=2 ! maximal value in v=v(rho)
  real, parameter, private :: alpha=7.5 
  real, parameter, private :: rhomax=9
  ! real,parameter, private :: relax_time=0.61 ! is not used in BR algorithm

  public :: VRhoF,CRhoF ! public functions

contains
  !>Dependency of the magnitude of the velocity on the density ---------------------------------------------
  function VRhoF(rho)
    real, intent(in) :: rho
    VRhoF=vmax*exp(-alpha*(rho/rhomax)**2) !v(rho)=vmax*e^(-alpha*(rho/rhomax)^2)) 


  end function VRhoF
  !>Cost function, dependency on the density---------------------------------------------------------------
  function CRhoF(rho)
    real, intent(in) :: rho
    CRhoF=1/VRhoF(rho)  ! c(rho)=1/V(rho)
  end function CRhoF
end module PedestrianSettings


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> FIFO implementation for vertices  of the FVM mesh
module vertqueue
  integer, dimension(:), allocatable,private :: vertexq !array for storing vertices in FIFO maner
  logical, dimension(:), allocatable,private :: isinq !array for information about vertices in the queue
  integer,private :: queuesize,first,last,qcount
contains

  !>create queue, alocate memory for queue
  subroutine QueueCreate(sizeofqueue)
    integer,intent(in) :: sizeofqueue
    queuesize=sizeofqueue	! size of the queue
    !write (*,*) 'Queue alocated'
    allocate(vertexq(1:queuesize), isinq(1:queuesize) ) ! alocate memory
    call QueueInitialize() ! prepare queue
  end subroutine QueueCreate
  
  subroutine QueueDeallocate()
    
    deallocate(vertexq, isinq)

  end subroutine QueueDeallocate

  !------------------------------------------------------------------------------------------------
  !> reinitialize created queue
  subroutine QueueInitialize()
    isinq=.FALSE. ! set all vertices to be out
    first=1
    last=1
    qcount=0 ! set count vertices in empty queue to 0
    !		do item=1,queuesize
    !			print *,item,isinq(item)
    !		end do
  end subroutine QueueInitialize
  !------------------------------------------------------------------------------------------------
  !>adds index of vertex into the queue
  subroutine QueueEnque(item)
    integer,intent(in) :: item

    !print*,'e3e3 ISINQ',item, size(isinq)
    if( .NOT.(isinq(item))) then ! check if vertex is not in the queue
       isinq(item)=.TRUE. ! mark vertex to be in the queue
       vertexq(last)=item ! add vertex to the array
       qcount=qcount+1
       !			print *, 'Add item into queue    ', item,'pos', last,'first',first,qcount
       last=last+1
       if(qcount>queuesize) then ! TODO comment
          print *,'ERROR queue is full' !queue is full for testing TODO comment!!
       end if! TODO comment
       if(last >queuesize) then ! check if end of the array is reached
          last=1
       end if

       !		else
       !			print *, 'Error: Vertex is in queue!',item
    endif
  end subroutine QueueEnque
  !------------------------------------------------------------------------------------------------		


  subroutine QueueDeque(item)
    integer::item
    if(qcount>0) then !in queue are vertices
       item=vertexq(first)
       !			print *, 'Remove item from queue ', item,'pos', last,'first',first,qcount
       isinq(item)=.FALSE.
       first=first+1
       qcount=qcount-1
       if(first >queuesize) then ! check if end of the array is reached
          first=1
       end if
    else
       print *, 'Error: Queue is empty',item
    endif
    return 
  end subroutine QueueDeque
  !------------------------------------------------------------------------------------------------
  integer function  QueueCount()
    !		print *,'Qcount:',qcount,'first:',first,'last:',last-1
    QueueCount=qcount		
    return 
  end function QueueCount
  !------------------------------------------------------------------------------------------------			
  subroutine QueueTestInfo()
    print *,'QUEUE TEST:'
    print *,'queue count', qcount
    do item=1,queuesize
       print *,item,isinq(item),vertexq(item)
    end do

  end subroutine QueueTestInfo
end module vertqueue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!>  Module for Bornemann Rash algorithm
module BRAlgorithm
  type Edge ! edge
     integer :: s,e !starting anf ending point of edge
  end type Edge

  type KSetInfo
     type(Edge), dimension(:),allocatable :: edges ! set of edge arround each vertex P_i
     real, dimension(:),allocatable :: calpha ! array of cosine alpha
     real, dimension(:),allocatable :: cbeta ! array of cosine beta
     logical :: isoutflow=.FALSE.
     integer :: count=0
  end type KSetInfo

  integer, dimension(:), allocatable,private :: sigmao ! set of internal vertices having neighbour on outflow, allocated in PrepareSigmao
  integer, private :: sigmaosize=0 ! size of sigmao set
  real, dimension(:),allocatable,private	:: phi ! array for storing potential
  real, dimension(:),allocatable,private	:: crho ! array for storing values of C(rho) on vertices of the mesh

  real, parameter, private :: BRTol=1e-12; ! tolerance of BR algorithm
  real, parameter, private :: BRInfinity=1e18 !this number is considered as infinity
  type(KSetInfo), dimension(:),allocatable,private :: bkpi ! set of edges forming boundary of K_pi


  public :: PrepareSigmao,PrepareBKPi,InitializeBRAlgorithm,PrepareVertexInfo,BRSolve,ComputeVelocity,Export ! public subroutines
  private :: AddEdgeToBKPi  ! private

contains
  !> prepare sigmo set ---------------------------------------------------------------------------------------
  subroutine PrepareSigmao(npoin,max_nbp, loc_ibp, ibp)
    implicit none
    integer, intent(in) :: npoin ! number of the grid nodes
    integer, intent(in) :: max_nbp ! maximal number of neighbouring nodes
    integer, dimension(1:npoin, 0:max_nbp), intent(in) ::loc_ibp !loc_ibp(i,0) = number of neigh nodes
    !loc_ibp(i,1:) = indexes of ^^^^^^^
    integer, dimension(1:npoin, 1), intent(in) ::ibp !loc_ibp(i) = -1 ==> internal node
    !loc_ibp(i) =  0 ==> node on boundary NOT OUTFLOW
    !loc_ibp(i) =  1 ==> node on boundary  OUTFLOW
    integer :: tmp,volume,vert,neibindex,neibvert

    allocate(sigmao(1:npoin)) ! alocate memory for sigmao set
    sigmaosize = 0
    !print*,'AAA:prepare sigmao set, sigmaosize is upper index ',sigmaosize

    do vert=1,npoin
       if(ibp(vert,1)==-1) then !if vertex is internal
          !			print *, 'Check vertex',vert
          do neibindex=1,loc_ibp(vert,0) ! check neighb of the vertex
             neibvert=loc_ibp(vert,neibindex) ! get index
             if(ibp(neibvert,1)==1) then ! neighb vertex is on the outflow
                sigmaosize=sigmaosize+1
                !print *, 'write',vert,'to',sigmaosize
                sigmao(sigmaosize)=vert 
                exit ! break inner loop, vertex is marked to be in sigmao			
             endif
          enddo
       endif
    enddo
    !print *, 'SigmaO set prepared: Count of vertices in this set ',sigmaosize
    
  end subroutine PrepareSigmao

  !> prepare set of edges forming boundary for each vertex : partial KP_i------------------------------------
  subroutine PrepareBKPi(npoin, nelem, max_nbp, x,lnd,ibp)
    implicit none
    integer, intent(in) :: npoin ! number of the grid nodes
    integer, intent(in) :: nelem ! number of the grid triangles
    integer, intent(in) :: max_nbp ! maximal number of neighbouring nodes
    real, dimension(1:npoin, 1:2), intent(in) :: x  ! coordinates of the nodes
    integer, dimension(1:nelem, 1:3), intent(in) :: lnd ! indexes of nodes forming each triangle
    integer, dimension(1:npoin, 1), intent(in) ::ibp !loc_ibp(i) = -1 ==> internal node
    !loc_ibp(i) =  0 ==> node on boundary NOT OUTFLOW
    !loc_ibp(i) =  1 ==> node on boundary  OUTFLOW

    integer :: volume,pi,pl,pk,tmp
    allocate(bkpi(1:npoin))	!alocate memory for information related to each vertex

    !print*,'...eder4..mh allocate(bkpi(1:npoin))?', allocated(bkpi)
    do volume=1,nelem	!foreach volume in set of volumes
       pi=lnd(volume,1)
       pk=lnd(volume,2)
       pl=lnd(volume,3)

       call AddEdgeToBKPi(npoin,pi,pk,pl,max_nbp,x(1:npoin, 1:2)) ! add edge <pk,pl> to pi
       call AddEdgeToBKPi(npoin,pk,pi,pl,max_nbp,x(1:npoin, 1:2)) ! add edge <pi,pl> to pk
       call AddEdgeToBKPi(npoin,pl,pi,pk,max_nbp,x(1:npoin, 1:2)) ! add edge <pi,pk> to pl
    enddo

    !storing info about outflow vertices
    do tmp=1,npoin ! foreach vertex of the mesh set if vertex is on outflow
       if(ibp(tmp,1)==1) then ! vertex is on the outflow
          bkpi(tmp)%isoutflow=.TRUE. !store this information
          !			print *,tmp,x(tmp,1),x(tmp,2)
       end if
    enddo
    !for testing
    !		do tmp=1,npoin
    !			print *,"vertex:",tmp
    !			do pi=1,bkpi(tmp)%count
    !				print *,bkpi(tmp)%edges(pi)%s,bkpi(tmp)%edges(pi)%e
    !			enddo
    !		enddo	
    ! TESTING		
    !print *, 'BKPi INFOS  for each vertex prepared'

    !print*,'...eder4..mh allocate(bkpi(1:npoin))?', allocated(bkpi),bkpi(1)%isoutflow

  end subroutine PrepareBKPi



  !> prepare potetntial and crho  array, alocate memory for potential and for crho--------------------------------------------
  subroutine PrepareVertexInfo(npoin)
    implicit none
    integer, intent(in) :: npoin ! number of the grid nodes
    allocate(phi(1:npoin),crho(1:npoin))
    ! TESTING			print *, 'Allocated memory for potential phi and for C(rho)'
  end subroutine PrepareVertexInfo

  !> initialize BR algorithm for iterate--------------------------------------------	
  !In original BR algorithm is potential at the vertices on the outflow set to zero, elsewhere is set to infinity.
  !In case of iterative usage for pedestrian flow we set last solution of the eikonnal equation as initial condition
  !for tho current step.
  subroutine InitializeBRAlgorithm(npoin,density,init)
    use vertqueue ! for queue of vertices
    use pedestriansettings ! parametres of the model
    implicit none
    integer, intent(in) :: npoin ! number of the grid nodes
    real, dimension(1:npoin), intent(in) :: density  ! density at vertexes - INPUT
    logical,intent(in) :: init ! .TRUE. use default initial condition,else use solution from last step
    integer :: vertex ! temporary variable for index of vertices
    integer :: pom 
    !set initial condition


    if(init) then !if is set, then use default initial condition, else last solution is used as initial condition
       ! TESTING			print *, 'Use default initialization'
       do vertex=1,npoin
          !!print*,'###edr4e323e3',vertex,npoin, allocated(bkpi)
          !!print*,'###edr4e323e3',vertex,npoin, allocated(bkpi),bkpi(vertex)%isoutflow
          if(bkpi(vertex)%isoutflow) then !in case vertex is on outflow
             phi(vertex)=0 
             !					print *,vertex,'set phi to ', phi(vertex) ! for testing
          else
             phi(vertex)=BRInfinity ! REMARK aproximation of infinity
          endif
       end do
    endif
    call QueueInitialize() ! initialize queue, important 
    do vertex=1,sigmaosize
       call QueueEnque(sigmao(vertex)) ! insert vertices from sigmao to queue
       !print *,'add vertex into queue',sigmao(vertex)
    end do
    ! TESTING	
    !print *,'Initial count vertices in the queue', QueueCount()

    do vertex=1,npoin ! fill array C(rho) on the vertices
       crho(vertex)= crhof(density(vertex)) !
    end do
    ! TESTING		print *, 'BR-algorithm initialized'

  end subroutine InitializeBRAlgorithm
  !>Solve eikonal equation via BR algorithm--------------------------------------------------------------------
  subroutine BRSolve(npoin,max_nbp,x,loc_ibp)
    use vertqueue ! for queue of vertices
    use pedestriansettings ! parametres of the model
    implicit none	
    integer, intent(in) :: npoin ! number of the grid nodes
    integer, intent(in) :: max_nbp ! maximal number of neighbouring nodes
    real, dimension(1:npoin, 1:2), intent(in) :: x  ! coordinates of the nodes
    integer, dimension(1:npoin, 0:max_nbp), intent(in) ::loc_ibp !loc_ibp(i,0) = number of neigh nodes
    !loc_ibp(i,1:) = indexes of ^^^^^^^

    integer :: pi !actual vertex
    integer :: pk !first vertex forming oposite edges to pi
    integer :: pl !second vertex forming oposite edges to pi
    integer :: neibindex,neibvert ! for storing vertex and index
    real :: iknorm,lknorm,ilnorm,delta,ksi,cl,cb,phij,minphij
    integer :: tmp ! temporary variables for edges forming PKP_i

    do while(QueueCount()>0)
       call QueueDeque(pi) !remove vertex from the queue

       minphij=BRInfinity ! REMARK approximation of infinity
       do tmp=1,bkpi(pi)%count !foreach edges in oposite to pi-minimize over this edge
          pk=bkpi(pi)%edges(tmp)%s
          pl=bkpi(pi)%edges(tmp)%e
          !print *, 'edge: ',pk,pl
          iknorm=crho(pi)*sqrt((x(pi,1)-x(pk,1))**2+(x(pi,2)-x(pk,2))**2) ! ||pi-pk||
          lknorm=crho(pi)*sqrt((x(pl,1)-x(pk,1))**2+(x(pl,2)-x(pk,2))**2) ! ||pl-pk||
          ilnorm=crho(pi)*sqrt((x(pi,1)-x(pl,1))**2+(x(pi,2)-x(pl,2))**2) ! ||pi-pl||
          delta=(phi(pl)-phi(pk))/lknorm
          cl=bkpi(pi)%calpha(tmp)
          cb=bkpi(pi)%cbeta(tmp) 
          !!!!ksi=cl*delta+sqrt((1-cl**2)*(1-delta**2))
          if (cl <= delta) then
             phij = phi (pk) + iknorm;
          else
             if (delta > -cb) then
	        ksi=cl*delta+sqrt((1-cl**2)*(1-delta**2))
                phij = phi (pk) + ksi * iknorm
             else
                phij = phi (pl) + ilnorm
             endif
          endif

          if (phij < minphij) then
             minphij = phij
          endif
       enddo
       if((phi(pi)-minphij)>BRTol) then
          do neibindex=1,loc_ibp(pi,0) ! check neighb of the vertex
             neibvert=loc_ibp(pi,neibindex) ! get index
             if(bkpi(neibvert)%isoutflow .EQV. .FALSE.) then
                call QueueEnque(neibvert) !put vertex into queue

             endif

          end do
       endif
       phi(pi)=minphij
    enddo

    !call	QueueTestInfo()
  end subroutine BRSolve
  !-------------------------------------------------------------------------------------------------------------------------------------
  !>Compute normalized gradient on each triange from piecevise linear solution. This gradient is normalized and then the velocity is computed.
  !Meaning of the computed velocity is V(rho)*mu
  subroutine ComputeVelocity(npoin, nelem, max_nbp, x,lnd,density,velocity)
    use pedestriansettings ! parametres of the model
    implicit none
    integer, intent(in) :: npoin ! number of the grid nodes
    integer, intent(in) :: nelem ! number of the grid triangles
    integer, intent(in) :: max_nbp ! maximal number of neighbouring nodes
    real, dimension(1:npoin, 1:2), intent(in) :: x  ! coordinates of the nodes
    integer, dimension(1:nelem, 1:3), intent(in) :: lnd ! indexes of nodes forming each triangle

    real, dimension(1:npoin), intent(in) :: density  ! density at vertexes - INPUT
    real, dimension(1:nelem, 1:2), intent(out) :: velocity  ! velocity on elements - OUTPUT
    real :: pix,piy,pjx,pjy,pkx,pky ! x,y coords of vertices of triangle
    real :: difphiij,difphiik !diference phi(pi)-phi(pj), phi(pi)-phi(pk)
    real :: detA1,detA2 !determinants for Cramers rule detA-neni potreba
    real :: norm,Velrho
    integer :: vol,pi,pj,pk
    do vol=1,nelem ! for each triangle
       pi=lnd(vol,1)
       pj=lnd(vol,2)
       pk=lnd(vol,3)
       !obtain position of vertices
       pix=x(pi,1)
       piy=x(pi,2)

       pjx=x(pj,1)
       pjy=x(pj,2)

       pkx=x(pk,1)
       pky=x(pk,2)
       !diferences of potential between vertices of triangle
       difphiij=phi(pi)-phi(pj)
       difphiik=phi(pi)-phi(pk)
       !determinants for Cramers rule			
       !  detA=(pix-pjx)*(piy-pky)-(piy-pjy)*(pix-pkx) !neni potreba, pri normalizaci se to vykrati
       detA1=difphiij*(piy-pky)-(piy-pjy)*difphiik !
       detA2=(pix-pjx)*difphiik-difphiij*(pix-pkx) ! 

       velocity(vol,1)=detA1!first component of gradient
       velocity(vol,2)=detA2!second component of gradient

       !normalizing gradient direction 
       norm=sqrt(velocity(vol,1)**2+velocity(vol,2)**2)

       ! REMARK direction of velocity mu is in OPOSITE direction of normalized gradient
       velocity(vol,1)=-velocity(vol,1)/norm
       velocity(vol,2)=-velocity(vol,2)/norm			

       ! HERE vector mu is computed. Next part is computing size i.e. mu*V(rho)
       !HACK-In FVM contex we need density over each volume, HERE we use average density avgdens=(rho(pi)+rho(pj)+rho(pk))/3
       !NEXT line compute average and then dependency of the velocity on the rho is evaluated		
       Velrho=VrhoF((density(pi)+density(pj)+density(pk))/3) !compute velocity based on averaged density
       velocity(vol,1)=velocity(vol,1)*Velrho
       velocity(vol,2)=velocity(vol,2)*Velrho
    end do

  end subroutine ComputeVelocity
!----------------------------------------------------------------------------------------------------

!> Deallocate all array used by vertex queue: 
  subroutine BRDeallocate()
    use  vertqueue
    integer :: test, vertex
    
    deallocate(sigmao,phi,crho, STAT=test) ! deallocate both arrays and test 
    if(test .NE. 0) then
       stop 'Deallocation BR arrays problem:sigmao,phi,crho'
    endif
    do vertex=1,size(bkpi) ! for each vertex in vertex inforamtion set
       deallocate(bkpi(vertex)%calpha,bkpi(vertex)%cbeta,bkpi(vertex)%edges,STAT=test)
       if(test .NE. 0) then
          stop 'Deallocation BR arrays problem-bkpi structure: calpha,cbeta,edges'
       endif
    end do
    deallocate(bkpi, STAT=test)
    if(test .NE. 0) then
       stop 'Deallocation BR arrays problem-bkpi structure'
    endif
  end subroutine BRDeallocate


  !> Export potential to gnuplot format---------------------------------------
  !Subroutine to export piecewise linear data on the mesh for gnuplot visualisation
  subroutine Export(npoin, nelem, max_nbp, x,lnd,velocity)
    implicit none
    integer, intent(in) :: npoin ! number of the grid nodes
    integer, intent(in) :: nelem ! number of the grid triangles
    integer, intent(in) :: max_nbp ! maximal number of neighbouring nodes
    real, dimension(1:npoin, 1:2), intent(in) :: x  ! coordinates of the nodes
    integer, dimension(1:nelem, 1:3), intent(in) :: lnd ! indexes of nodes forming each triangle
    real, dimension(1:nelem,1:2), intent(in) :: velocity  ! density at vertexes - INPUT

    integer :: pi,pj,pk,i;
    OPEN(UNIT=12, FILE="phi.dta", ACTION="write", STATUS="replace")
    DO i=1,nelem
       pi=lnd(i,1)
       pj=lnd(i,2)
       pk=lnd(i,3)
       WRITE(12,*) x(pi,1),x(pi,2),phi(pi)
       WRITE(12,*) x(pj,1),x(pj,2),phi(pj)
       WRITE(12,*) 
       WRITE(12,*) x(pk,1),x(pk,2),phi(pk)
       WRITE(12,*) x(pk,1),x(pk,2),phi(pk)
       WRITE(12,*) 
       WRITE(12,*) 
    END DO
    close(12)

    OPEN(UNIT=12, FILE="vx.dta", ACTION="write", STATUS="replace")
    DO i=1,nelem
       pi=lnd(i,1)
       pj=lnd(i,2)
       pk=lnd(i,3)
       WRITE(12,*) x(pi,1),x(pi,2),velocity(i,1)
       WRITE(12,*) x(pj,1),x(pj,2),velocity(i,1)
       WRITE(12,*) 
       WRITE(12,*) x(pk,1),x(pk,2),velocity(i,1)
       WRITE(12,*) x(pk,1),x(pk,2),velocity(i,1)
       WRITE(12,*) 
       WRITE(12,*) 
    END DO
    close(12)
    OPEN(UNIT=12, FILE="vy.dta", ACTION="write", STATUS="replace")
    DO i=1,nelem
       pi=lnd(i,1)
       pj=lnd(i,2)
       pk=lnd(i,3)
       WRITE(12,*) x(pi,1),x(pi,2),velocity(i,2)
       WRITE(12,*) x(pj,1),x(pj,2),velocity(i,2)
       WRITE(12,*) 
       WRITE(12,*) x(pk,1),x(pk,2),velocity(i,2)
       WRITE(12,*) x(pk,1),x(pk,2),velocity(i,2)
       WRITE(12,*) 
       WRITE(12,*) 
    END DO
    close(12)

  end subroutine Export


  !***********************----------- PRIVATE SUBROUTINES -------------*********************
  subroutine AddEdgeToBKPi(npoin,pi,pk,pl,max_nbp,x)
    implicit none
    integer, intent(in) :: npoin ! number of the grid nodes
    integer,intent(in) :: pi
    integer,intent(in) :: pk
    integer,intent(in) :: pl
    integer, intent(in) :: max_nbp ! maximal number of neighbouring nodes
    real, dimension(1:npoin, 1:2), intent(in) :: x  ! coordinates of the nodes
    real:: cl,cb ! cos alpha, cos beta
    real :: iknorm,lknorm,ilnorm
    integer ::tmp
    !------prepare inner boundary over vertex 
    if(bkpi(pi)%count == 0) then !no edge added to information
       allocate(bkpi(pi)%edges(1:max_nbp),bkpi(pi)%calpha(1:max_nbp),bkpi(pi)%cbeta(1:max_nbp) )! allocate memory for edges,calpha,cbeta		
    endif
    !topological information
    bkpi(pi)%count=bkpi(pi)%count+1
    tmp=bkpi(pi)%count
    bkpi(pi)%edges(tmp)%s=pk
    bkpi(pi)%edges(tmp)%e=pl
    !computation of geometrical infox--cos alpha, cos beta
    iknorm=sqrt((x(pi,1)-x(pk,1))**2+(x(pi,2)-x(pk,2))**2) ! ||pi-pk||
    lknorm=sqrt((x(pl,1)-x(pk,1))**2+(x(pl,2)-x(pk,2))**2) ! ||pl-pk||
    ilnorm=sqrt((x(pi,1)-x(pl,1))**2+(x(pi,2)-x(pl,2))**2) ! ||pi-pl||


    cl=((x(pi,1)-x(pk,1))*(x(pl,1)-x(pk,1))+(x(pi,2)-x(pk,2))*(x(pl,2)-x(pk,2)))/(iknorm*lknorm); ! compute (Pi-Pk)*(Pl-Pk)/(||pi-pk||*||pl-pk||)
    cb=((x(pi,1)-x(pl,1))*(x(pk,1)-x(pl,1))+(x(pi,2)-x(pl,2))*(x(pk,2)-x(pl,2)))/(ilnorm*lknorm); ! compute (Pi-Pl)*(Pk-Pl)/(||pi-pl||*||pl-pk||)
    bkpi(pi)%calpha(tmp)=cl
    bkpi(pi)%cbeta(tmp)=cb

  end subroutine AddEdgeToBKPi


end module BRAlgorithm
