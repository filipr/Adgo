!> subroutines for interpolation in space-time discontinuous Galerkin method
module st_interpol
  use main_data
!  use problem_oper
  use geometry
  use mesh_oper

  implicit none

  ! public :: EvalTimeJumps
  !public :: SeekIntersecLine
  public :: IntersectGridTri
  ! public :: MaxDiam
  public :: IntersectTris
  public :: PosNodes
  public :: Copy
  public :: AreaTri
  public :: GetBarCoord
  type :: intersect
     !already tested; have intersect
     logical,allocatable, dimension(:)	:: done, had
     !count of points which intersects K
     integer,allocatable, dimension(:)	:: NumPnt
     !points which intersects K in order which makes polygonal
     real,allocatable, dimension(:,:,:) :: IntPnt
  end type intersect

  type :: linter
     ! count of points in IntPnt, number of elements i and K
     integer			:: NumPnt,i,K
     ! intersection points
     real,dimension(1:6,1:3)	:: IntPnt
     !pointer to next linter in the linked list
     type(linter), pointer	:: next
  end type linter

  !linked list of triangles which will be counted
  type :: link
     integer 		 :: i
     type(link), pointer :: next
  end type link

  !remember which vertex is of this type and how many vertices is of this type
  !in the triangle
  !types for inside, edge, vertex, outsided will be used
  type :: pos
     integer 			:: num		!number of this type
     logical, dimension (1:3)	:: is		!is this edge this type
  end type pos



contains

  ! artificialy ceated new grid, readed from a file
  subroutine SetNewGrid(gridN, gridfile )
    type(mesh)		:: gridN
    character(len=*)	:: gridfile

    !    allocate( gridN)
    call gridN%read( gridfile )
    !    call ReadMesh(gridfile, gridN)

    call gridN%seekNeighbours()
    !print*,'Neigbours found'

    if (nbDim==2) then
       call SeekCurvedBoundary(gridN)
       ! call SplineCurvedBoundary( )
       !print*,'Curved boundary found'
    endif

    ! not necessary - is called in prepare problem!
    !call gridN%computeGeometry()
    !call ComputeGeometry(gridN)
    !print*,'Geometry computed'

   stop 'FR PrepareProblem(gridN, rsolfileAB) COMMENTED'
   ! call PrepareProblem(gridN, 'rsolfileAB')



  end subroutine SetNewGrid





  !find all triangles with non-empty intersect with triangle with number K in gridK
  ! all the coordinates of the intersections are savedin inter, numbers of triangles
  !which has non-empty intrsection with K are saved in array triangles, NumTri is
  ! number of these triangles
  subroutine IntersectGridTri(grid,gridK,K,inter,triangles,NumTri)
    type(intersect),intent(inout)		:: inter	!triangles & intersections
    type(mesh), intent(in)		:: grid, gridK
    integer,intent(in)			:: K
    integer,dimension(:),intent(out)	:: triangles	!triangles with non-empty 								intersection with K
    !type(link),pointer,intent(out)		:: f
    integer,intent(out)			:: NumTri	!number of triangles
    integer				:: i,j,l,m
    real					:: radius	!radius of circle, which 								decide if the bar. is close
    real, dimension (1:2)			:: bar,bari	!bar. of K and elem. i
    type(link),pointer			:: current
    type(link),pointer			:: first,help
    real, dimension (1:3,1:2)		:: nodesK,nodesi
    logical				:: fin		!the end of whole subroutine
    logical				:: yet		!help for K is included in i
    real					:: max_diam, min_diam
    logical				:: sam
    !logical,dimension(:),intent(inout) :: u
    logical				:: nothing

    !initiate

    !allocate(triangles(1:grid%nelem))
    !allocate(inter%NumPnt(1:grid%nelem))
    !allocate(inter%done(1:grid%nelem))
    !allocate(inter%had(1:grid%nelem))
    !allocate(inter%IntPnt(1:grid%nelem,1:6,1:2))

    !print*,"starting of initiation", state%space%adapt%adapt_level

    sam=.false.
    nothing=.true.
    ! print*,"noth"
    triangles(1:grid%nelem) = 0
    ! print*,"tri"
    inter%NumPnt(1:grid%nelem) = 0
    inter%done(1:grid%nelem)=.false.
    inter%had(1:grid%nelem)=.false.
    inter%IntPnt(1:grid%nelem,1:6,1:2) = 0
    fin=.false.
    NumTri=0
    ! print*,"heeere"
    nullify(first)
    ! print*,"nulify"
    nullify(current)
    ! print*,"initiate"
    !nullify(c)
    !nullify(f)

    max_diam = MaxDiam(grid)
    min_diam=MinDiam(grid)
    if (K==1) then
       print*,"max_diam", max_diam
       print*,"min_diam",min_diam
    end if
    radius=(max_diam)
    !print*,"radius", radius
    i=1
    !barycentre of K element
    bar(1:2) = gridK%elem(K)%xc(1:2)
    !print*,"barycentre of K is", bar(1:2)




    !open(11,file='barycentre')
    do
       !print*,"i",i
       !end or if nothing was found, then try with twice bigger radius once more
       if ( (i==grid%nelem+1) .and. (NumTri<=0) )then
          if (nothing .eqv. .true.) then
             i=1
             radius=2*max_diam
             nothing=.false.
          else
             print*, "Find no triangle with non-empty intersection with K"
             !u(k)=.false.
             return
          end if
       else if (i==grid%nelem+1) then
          !	         print*, "finished"
          return
       end if

       !barycentre of i element
       bari(1:2) = grid%elem(i)%xc(1:2)
       !write(11,*) "barycentre of",i,"is", bari(1:2)
       !write(11,*) "dist is", Dist(bar,bari)
       ! if distance is less than or equal to radius then try to find intersection with K
       if (Dist(bar,bari) <= radius) then
          nothing=.false.
          !if the intersection of K and i element is non-empty
          !then it's accomplished & add i element to triangles
          !		  print*,"IntersectTris on element",i,"from grid and ",K,"from gridK"
          !		  print*,"inter%NumPnt(i)",inter%NumPnt(i)
          call IntersectTris(i, K, grid, gridK, inter, fin)
          !		  print*,"inter%had(i) 1st",inter%had(i)
          !		  PRINT*,"inter%NumPnt(i) 1st",inter%NumPnt(i)
          !	  	  print*,"NumTri",NumTri
          if (inter%had(i) .eqv. .true.) then
             NumTri=NumTri+1
             !			print*,"NumTri",NumTri
             !			print*,"triangles",triangles(NumTri)
             !call Put(i,f,c)
             triangles(NumTri)=i
             !			print*,"triangles(NumTri)",triangles(NumTri)
             ! case of 3 intersectionpoints checking f we're done or not
             if (inter%NumPnt(i) == 3) then
                nodesK(1:3,1:2)=gridK%x(gridK%elem(K)%face(idx,1:3),1:2)
                !				print*,"nodesK",nodesK
                !				print*,"IntPnt",inter%IntPnt(i,1:3,1:2)
                !if Same is true then all nodesK have been added
                !so we`re done
                call Same(nodesK,inter%IntPnt(i,1:3,1:2),sam)
                !				print*,"sam",sam
                if (sam .eqv. .true.) then
                   !					print*,"ended for same"
                   return
                end if
             end if
             !fin is alternative of same, fin is fill inside the subroutin
             if (fin .eqv. .true.) then
                !				print*,"Ended for fin"
                return
             end if
             !			print*,"exit"
             exit
          end if
       end if
       i=i+1
    end do
    !  print*,"sousedi",grid%elem(i)%face(neigh,1:3)
    !adding neighbours to the list
    add1:do j=1,3
       if (grid%elem(i)%face(neigh,j) > 0) then
          allocate(current)
          call Put(grid%elem(i)%face(neigh,j),first,current)
          !		print*, "Neighbour number",grid%elem(i)%face(neigh,j),"was added"
       end if
    end do add1
    ! if first triangle has no neighbours then we`re finished
    if ((grid%elem(i)%face(neigh,1) < 0) .and. (grid%elem(i)%face(neigh,2) < 0) &
         .and. (grid%elem(i)%face(neigh,3) < 0) ) then
       !print*, "No more triangles to be computed"
       return
    end if

    do
       !print*,"inter%done(first%i)",inter%done(first%i)
       !if i-element hasn't been counted then count him
       if (associated(first) .eqv. .true.) then
          if (inter%done(first%i) .eqv. .false.) then
             i=first%i
             !		print*,"IntersectTris called on element",i
             !count intersection of two triangles (i&K) from grid and gridK
             call IntersectTris(i, K, grid, gridK, inter, fin)
             if (inter%had(i) .eqv. .true.) then
                !			print*,i,"element has intersection"
                NumTri=NumTri+1
                triangles(NumTri)=i
                !call Put(i,f,c)
                !			print*,"inter%NumPnt(i)",inter%NumPnt(i)
                ! case of 3 intersectionpoints checking f we're done or not
                if (inter%NumPnt(i) == 3) then
                   nodesK(1:3,1:2)=gridK%x(gridK%elem(K)%face(idx,1:3),1:2)
                   call Same(nodesK,inter%IntPnt(i,1:3,1:2),sam)
                   if (sam .eqv. .true.) then
                      !					!print*,"ended for same,second case"
                      return
                   end if
                end if
                if (fin .eqv. .true.) then
                   !				!print*,"Ended for fin, second case"
                   return
                end if
                !exit
             end if
             !if i-element has non-empty intersection with K,
             !then add his neighbours to the list
             if (inter%had(i) .eqv. .true.) then
                add2:do j=1,3
                   if (grid%elem(i)%face(neigh,j) > 0) then
                      allocate(current)
                      call Put(grid%elem(i)%face(neigh,j),first,current)
                      !print*, "Neighbour number",grid%elem(i)%face(neigh,j),"was added"
                   end if
                   if ((grid%elem(i)%face(neigh,1) < 0) .and. (grid%elem(i)%face(neigh,2) < 0) &
                        .and. (grid%elem(i)%face(neigh,3) < 0) .and. (associated(first%next) .eqv. .false.)) then
                      !print*, "No more triangles to be computed"
                      return
                   end if
                end do add2
             end if
  	  end if
       end if
       !if there's another element in the list, then continue
       !otherwise we're done
       if ((associated(first%next) .eqv. .true.).and.(inter%done(first%i) .eqv. .true.)) then
          help => first%next
          deallocate(first)
          allocate(first)
          first => help
       else
          if (associated(first%next) .eqv. .false.) then
             exit
          end if
       end if
    end do
  end subroutine IntersectGridTri



  !count concrete intersection of two triangles, fulfil type intersect
  subroutine IntersectTris(i,K,grid,gridK,inter,fin)
    real, dimension(1:3,1:2) 			:: nodesK, nodesi
    integer,intent(in)		 		:: i,K		!number of element i,K
    type(mesh), intent(in)			:: grid, gridK  !grids
    type(intersect), intent(inout)		:: inter	!table to fulfil
    logical, intent(inout)			:: fin		!variable to end subroutine
    character(len=2), dimension(1:3)		:: post		!position of vertices
    type(pos)					:: ins,edg,ver,ous!numbers of inside, edge,
    !vertex, outside vertices
    integer					:: j,l,m,n,cnt,cnt2,cnt3,edge,edge2,vrt,ed,ou,vr,vr1,vr2,vr3,vr4,lin
    logical, dimension(1:3)			:: inner,inner2
    real, dimension(1:3,1:2)			:: xi,xi2
    real, dimension(1:2)				:: x,x2,finpoint
    logical					:: inn,inn2,difer
    integer, dimension(1:2)			:: edgi, e
    logical, dimension(1:6)			:: line,dif
    character(len=2)				:: posi
    real, dimension(1:3)				:: edg1,edg2,edg3
    logical, dimension(1:3,1:3)			:: inner3
    real, dimension(1:3,1:3,1:2)			:: xi3
    real, parameter				:: TOL=0.000005


    !initiate

    fin=.false.

    if (inter%done(i) .eqv. .true.) then
       !print*,"return for done called"
       return
    end if

    inter%done(i)=.true.
    !vertices of K and i triangle
    nodesi(1:3,1:2)=grid%x(grid%elem(i)%face(idx,1:3),1:2)
    nodesK(1:3,1:2)=gridK%x(gridK%elem(K)%face(idx,1:3),1:2)

    !open(13,file="cur-ele")
    !write(13,*),nodesi(1,1:2),i
    !write(13,*),nodesi(2,1:2),i
    !  write(13,*),nodesi(3,1:2),i
    ! write(13,*),nodesi(4,1:2),i
    !close(13)
   !open(21,file="Pos")
   !i-element and K-element are the same then we fill nodesK and we are done

    if ((abs(nodesK(1,1)-nodesi(1,1))<= TOL) .and. (abs(nodesK(1,2)-nodesi(1,2))<= TOL) .and. &
         (abs(nodesK(2,1)-nodesi(2,1))<= TOL) .and. (abs(nodesK(2,2)-nodesi(2,2))<= TOL) .and. &
         (abs(nodesK(3,1)-nodesi(3,1))<= TOL) .and. (abs(nodesK(3,2)-nodesi(3,2))<= TOL))  &
         then
       inter%had(i)=.true.
       inter%NumPnt(i)=3
       inter%IntPnt(i,1:3,1:2)=nodesK(1:3,1:2)
       fin=.true.
       !print*, i, "and", K, "are exactly the same triangles"
       !print*,nodesK(1,1)-nodesi(1,1),nodesK(1,2)-nodesi(1,2),nodesK(2,1)-nodesi(2,1),nodesK(2,2)-nodesi(2,2)
       !print*,nodesK(3,1)-nodesi(3,1),nodesK(3,2)-nodesi(3,2)
       return
    end if
    !print*,"on",i,"element IntersectTris has been called"


   !if i-element and K-element aren't the same...
    !we shall find out where are they
    !subroutine finds out the position of nodesi in nodesK, fill post with position
    !of nodesi(just like "VA","AB","OU","IN") and ins,edg,ver,ous are also filled
    call PosNodes(nodesi(1:3,1:2),nodesK,post,ins,edg,ver,ous)
    if (k==36) then
       open(14,file="Pos")
       write(14,*),"nodesi,i",i
       write(14,*),nodesi(1,1:2)
       write(14,*),nodesi(2,1:2)
       write(14,*),nodesi(3,1:2)
       write(14,*),"post", post
    end if

    !MAIN

    !selection by number of inside vertices
    select case(ins%num)
       ! i is included in K
       ! case of ins%num
    case(3)
       inter%had(i) = .true.
       !print*,"inter%NumPnt(i)=3 3rd"
       inter%NumPnt(i) = 3
       inter%IntPnt(i,1:3,1:2) = nodesi(1:3,1:2)
       !fin=.true.
       !print*, "3 inside vertices, we're done for",i, "and", K
       return
       ! case of ins%num
    case(2)
       !print*,"2 inside"
       !if the 3rd vertex is on edge or it's on the vertex then
       ! intersection are all nodesi and we are done
       if ( (edg%num == 1) .or. (ver%num == 1) ) then
          inter%had(i) = .true.
          inter%NumPnt(i) = 3
          !print*,"inter%NumPnt(i)=3,,4"
          inter%IntPnt(i,1:3,1:2) = nodesi(1:3,1:2)
          !print*, "3rd on edge or vertex"
          return
       end if

       !the 3rd vertex is outside ---> 4 intersection points
       inter%had(i) = .true.
       inter%NumPnt(i) = 4
       !print*,"inter%NumPnt(i)=4"
       !which vertex is the one outside
       j=1
       do
          if (ous%is(j) .eqv. .true.) then
             exit
          end if
          j=j+1
       end do
       !if j>3 then some error must have occured, cause we have just 3 vertices
       if (j>3) then
          print*, "ERROR in IntersectTris"
       end if

       !if j isn't one then we should swap nodes as if the one outside was node number one
       if ( j/=1 ) then
          call SwapNodes(nodesi,j)
       end if
       !1
       call IntersectLineTri1(nodesi(1,1:2),nodesi(2,1:2),nodesK,inn,x)
       call IntersectLineTri1(nodesi(1,1:2),nodesi(3,1:2),nodesK,inn2,x2)
       !add the inside vertices and one intersection point
       if ( inn .eqv. .true.) then
          inter%IntPnt(i,1,1:2) = x(1:2)
          !print*,"adding inter with 1&2",x(1:2)
          inter%IntPnt(i,2:3,1:2) = nodesi(2:3,1:2)
          !print*,"adding nodesi",nodesi(2,1:2),nodesi(3,1:2)
       else
          print*, "error case ins%num=2"
       end if
       !add the last intersection point
       if ( inn2 .eqv. .true.) then
          inter%IntPnt(i,4,1:2) = x2(1:2)
          !print*,"adding inter with 1&3",x2(1:2)
       else
          print*,"Errorcase ins%num=2"
       end if
       ! knowing if some nodesK is in the nodesi, if yes and if it differs
       ! from already added then add him too
       do l=1,3
          call Inside2(nodesK(l,1:2),nodesi,line(l))
          if (line(l).eqv..true.) then
             call Differ(nodesK(l,1:2),inter%IntPnt(i,1,1:2),dif(1))
             call Differ(nodesK(l,1:2),inter%IntPnt(i,2,1:2),dif(2))
             call Differ(nodesK(l,1:2),inter%IntPnt(i,3,1:2),dif(3))
             call Differ(nodesK(l,1:2),inter%IntPnt(i,4,1:2),dif(4))
             if ((dif(1).eqv..false.).and.(dif(2).eqv..false.).and.(dif(3).eqv..false.)&
                  .and.(dif(4).eqv..false.))then
                inter%NumPnt(i)=inter%NumPnt(i)+1
                inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(l,1:2)
                !print*,"adding NodesK",l,nodesK(l,1:2)
             end if
          end if
      end do

       ! case of ins%num
       !one vertex inside
    case(1)
       !print*, "1 inside vertex"
       !number of vertices on the edge & vertex
       select case(edg%num+ver%num)
          !number of vertices on the edge & vertex
          !1 inside, 2 on edge or on the vertex
       case(2)
          !print*, " 2 on edge or vertex"
          inter%had(i) = .true.
          inter%NumPnt(i) = 3
          !!print*,"inter%NumPnt(i)=3,,5"
          inter%IntPnt(i,1:3,1:2) = nodesi(1:3,1:2)
          !all nodesi have been added
          return
          !number of vertices on the edge & vertex
          ! 1 inside, 1 outside, 1 on the edge or on the vertex
       case(1)
          !print*,"1 inside, 1 outside, 1 on the edge or on the vertex"
          !which vertex is the one outside
          j=1
          do
             if (ous%is(j) .eqv. .true.) then
                exit
             end if
             j=j+1
          end do
          !if j>3 then some error must have occured, cause we have just 3 vertices
          if (j>3) then
             print*, "ERROR in IntersectTris"
          end if
          !print*,"j",j
          !if j isn't one then we should swap nodes as if the one outside was node number one
          if ( j/=1 ) then
             call SwapNodes(nodesi,j)
          end if
          !add the vertices inside and on the edge
          inter%had(i) = .true.
          inter%NumPnt(i) = 2
          inter%IntPnt(i,1:2,1:2) = nodesi(2:3,1:2)
          !print*,"adding inside and edge",nodesi(2:3,1:2)
          !if 1 & 3 has intersect with K, different from 3 then add it
          call IntersectLineTri3(nodesi(1,1:2),nodesi(3,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
          if  (inn .eqv. .true.) then
             inter%NumPnt(i) = inter%NumPnt(i)+1
             inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
          end if
          !if 1 & 2 has intersect with K, different from 2 then add it
          call IntersectLineTri3(nodesi(1,1:2),nodesi(2,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
          if  (inn .eqv. .true.) then
             inter%NumPnt(i) = inter%NumPnt(i)+1
             inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
          end if
          if (inter%NumPnt(i) >=3) then
             return
          else
             print*, "Error in 1 inside 1 on the edge/vertex, 1 outside"
          end if


          !number of vertices on the edge & vertex
          !1 inside, 2 outside
       case(0)
          !print*," 1 inside, 2 outside"
          !which one is the one inside
          j=1
          do
             if (ins%is(j) .eqv. .true.) then
                exit
             end if
             j=j+1
          end do
          !if j>3 then some error must have occured, cause we have just 3 vertices
          if (j>3) then
             print*, "ERROR in IntersectTris"
          end if
          !if j isn't one then we should swap nodes as if the one outside was node number one
          if ( j/=1 ) then
             call SwapNodes(nodesi,j)
             !print*,"SwapNodes called on nodesi and", j
          end if
          if ((k==36).and.(i==695)) then
             open(20,file="ins")
             write(20,*),"K-oo inside, 2outside",K,i

          end if
          if (j==2) then
             call SwapNodesKJ(nodesi,2,3)
             if ((k==36).and.(i==695))then
                write(20,*),"swap called"
             end if
          end if

          !the one inside is for sure the one for the intersect
          inter%NumPnt(i) = 1
          inter%had(i) = .true.
          inter%IntPnt(i,1,1:2) = nodesi(1,1:2)

          !print*,"adding of nodesi(1)",nodesi(1,1:2)

          call IntersectLineTri2(nodesi(2,1:2),nodesi(3,1:2),nodesK,inner,xi,cnt)
          if ((k==36).and.(i==695)) then
             !write(20,*),nodesi(1,1:2)
             !write(20,*),nodesi(2,1:2)
             !write(20,*),nodesi(3,1:2)
             write(20,*),"cnt",cnt
             write(20,*),"nodesi",1,nodesi(1,1:2)
          end if
          !print*,"cnt",cnt

          !selecting according to value of cnt = "number of intersection of nodesi 2 & 3 and K"

          select case(cnt)
             !B
             !number of intersection of nodesi 2 & 3 and K
          case(0)
             !have to find intersects with 1 & 2 and 1 & 3
             !finding and adding intersect 1 & 2 to the intersectionts points
             edge=0
             call IntersectLineTri3b(nodesi(1,1:2),nodesi(2,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2),edge)
             edge2=0
             !finding and adding intersect 1 & 3 to the intersectionts points
             call IntersectLineTri3b(nodesi(1,1:2),nodesi(3,1:2),nodesK,inn2,x2,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2),edge2)
             !print*,"inn",inn
            if ((k==36).and.(i==695)) then
               write(20,*),"inn,inn2",inn,inn2
               write(20,*),"edge,edge2",edge,edge2
            end if
             if ( inn .eqv. .true. ) then
                inter%NumPnt(i) = inter%NumPnt(i)+1
                inter%IntPnt(i, inter%NumPnt(i),1:2) = x(1:2)
                !if ((K==200).and.(i==171)) then
                !	write(20,*),"NumPnt",inter%NumPnt(i)
                !	write(20,*),"x(1:2)",x(1:2)
                !end if

                !print*,"adding point-intersect with 1&2",x(1:2)

             else
                print*, "Error in case B"
                if (inn2 .eqv. .true.) then
                   inter%NumPnt(i) = inter%NumPnt(i)+1
                   inter%IntPnt(i, inter%NumPnt(i),1:2) = x2(1:2)
                end if
                do l=1,3
                   call Inside(nodesK(l,1:2),nodesi,line(l))
                   !if ((K==200).and.(i==171)) then
                   !write(20,*),"line(l)",line(l)
                   !end if
                   if (line(l).eqv..true.) then

                      call Differ(nodesK(l,1:2),inter%IntPnt(i,inter%NumPnt(i),1:2),dif(1))
                      !	if ((K==200).and.(i==171)) then
                      !	write(20,*),"dif",dif(1)

                      !	end if
                      if (dif(1).eqv..false.) then
                         inter%NumPnt(i)=inter%NumPnt(i)+1
                         inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(l,1:2)
                         !print*,"add nodesK-dole",l
                      end if
                   end if
                end do
                return
             end if

             !print*,"inn2",inn
             !inn .eqv. .true. and has been already added
             if ( inn2 .eqv. .true. ) then
                if (edge == edge2) then
                   inter%NumPnt(i) = inter%NumPnt(i)+1
                   inter%IntPnt(i, inter%NumPnt(i),1:2) = x2(1:2)
                   !if ((K==200).and.(i==171)) then
                   !write(20,*),"NumPnt",inter%NumPnt(i)
                   !write(20,*),"x2(1:2)",x2(1:2)
                   !end if

                   !print*,"adding point-intersect with 1&3",x2(1:2)

                   return
                else
                   !inter%NumPnt(i) = inter%NumPnt(i)+1
                   !if ((edge<=3).and.(edge2<=3)) then
                   !print*,"edge+edge2",edge+edge2
                   !print*,"edge",edge
                   !gets vertices of edge
                   call GetVertex(edge,vr1,vr2)
                   !get third vertex from two
                   call GetThirdVertex(vr1,vr2,vr)
                   !if nodesK(vr1) then vr2 then the 3rd
                   !are inside nodesi then add him
                   call Inside(nodesK(vr1,1:2),nodesi,line(1))
                   if (line(1).eqv..true.) then
                      inter%NumPnt(i)=inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(vr1,1:2)
                      !print*,"adding nodesK-nahore vr1",vr1
                   end if
                   call Inside(nodesK(vr2,1:2),nodesi,line(2))
                   if (line(2).eqv..true.) then
                      inter%NumPnt(i)=inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(vr2,1:2)
                      !print*,"adding nodesK-nahore vr2",vr2
                   end if

                   call Inside(nodesK(vr,1:2),nodesi,line(2))

                   if (line(2).eqv..true.) then
                      inter%NumPnt(i)=inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(vr,1:2)
                      !print*,"adding nodesK-nahore vr",vr
                   end if

                   inter%NumPnt(i) = inter%NumPnt(i)+1
                   inter%IntPnt(i, inter%NumPnt(i),1:2) = x2(1:2)
                   !print*,"adding point-intersect with 1&3",x2(1:2)
                   return
                end if
             else
                ! if some nodesK is inside then add him

                do l=1,3
                   call Inside(nodesK(l,1:2),nodesi,line(l))
                   if (line(l).eqv..true.) then
                      call Differ(nodesK(l,1:2),x(1:2),dif(1))
                      if (dif(1).eqv..false.) then
                         inter%NumPnt(i)=inter%NumPnt(i)+1
                         inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(l,1:2)
                      end if
                      !print*,"add nodesK-dole",l
                   end if
                end do
                print*, "Error in case B"
             end if
             !number of intersection of nodesi 2 & 3 and K
          case(1)
             !have to find intersects with 1 & 2 and 1 & 3
             !finding and adding intersect 1 & 2 to the intersectionts points
             call IntersectLineTri3(nodesi(1,1:2),nodesi(2,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
             if ( inn .eqv. .true. ) then
                inter%NumPnt(i) = inter%NumPnt(i)+1
                inter%IntPnt(i, inter%NumPnt(i),1:2) = x(1:2)
                !adding the intersect of nodesi 2,3
                do l=1,3
                   if (inner(l) .eqv. .true.) then
                      inter%NumPnt(i) = inter%NumPnt(i)+1
                      inter%IntPnt(i, inter%NumPnt(i),1:2) = xi(l,1:2)
                   end if
                end do
             else
                print*, "Error in case C"
             end if
             !finding and adding intersect 1 & 3 to the intersectionts points
             call IntersectLineTri3(nodesi(1,1:2),nodesi(3,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
             if ( inn .eqv. .true. ) then
                inter%NumPnt(i) = inter%NumPnt(i)+1
                inter%IntPnt(i, inter%NumPnt(i),1:2) = x(1:2)
             else
                print*, "Error in case C"
             end if
             !number of intersection of nodesi 2 & 3 and K
          case(2)
             !first and second point
             call IntersectLineTri3b(nodesi(1,1:2),nodesi(2,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2),edge)
             if (inn.eqv..true.)then
                inter%NumPnt(i) = inter%NumPnt(i)+1
                 inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                if((k==36).and.(i==695)) then
                   write(20,*),"intersection with 1&2", inter%NumPnt(i),x(1:2)
                   write(20,*),"edge",edge
                end if
             end if
             !print*,"inn",inn
             !print*,"adding point",x(1:2)

             !print*,"edge",edge
             !intersection of nodesi 1 & 2 is the edge number
             select case(edge)
             case(1)
                !3rd point is on the K 1 & 2
                call SeekIntersect(nodesi(2,1:2),nodesi(3,1:2),nodesK(1,1:2),nodesK(2,1:2),x,inn)

             case(2)
                !3rd point is on the K 2 & 3
                call SeekIntersect(nodesi(2,1:2),nodesi(3,1:2),nodesK(2,1:2),nodesK(3,1:2),x,inn)
             case(3)
                !3rd point is on the K 1 & 3
                call SeekIntersect(nodesi(2,1:2),nodesi(3,1:2),nodesK(1,1:2),nodesK(3,1:2),x,inn)
             case default
                print*, "Error in A"
             end select
             line(2)=.true.
             !3rd point
             if (inn .eqv..true.) then
                inter%NumPnt(i) = inter%NumPnt(i)+1

                inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                line(2)=.false.
                if((k==36).and.(i==695)) then
                   write(20,*),"intersection with edge", inter%NumPnt(i),x(1:2)
                   write(20,*),"edge",edge
                end if
             end if
             if (edge==3) then
                edge=1
             else
                edge=edge+1
             end if
             edgi(1)=edge
             call Inside2(nodesK(edge,1:2),nodesi,line(1))
             !line(4)=.false.
             if((k==36).and.(i==695)) then
                write(20,*),"inside edge",edge,line(1)
             end if
             !if ((inn.eqv..false.).and.(line(1).eqv..false.)) then
             !	if((k==36).and.(i==695)) then
             !	write(20,*),"other edge"
             !	end if
             !	line(4)=.true.
             !	if(edge/=3) then
             !		edge=edge+1
             !	else
             !		edge=1
             !	end if
             !call Inside2(nodesK(edge,1:2),nodesi,line(1))
             !if((k==36).and.(i==695)) then
             !write(20,*),"inside edge",edge,line(1)
             !end if
             !end if
             if (line(1).eqv..true.) then
                call Differ(nodesK(edge,1:2),inter%IntPnt(i,1,1:2),dif(1))
                call Differ(nodesK(edge,1:2),inter%IntPnt(i,2,1:2),dif(2))
                call Differ(nodesK(edge,1:2),inter%IntPnt(i,3,1:2),dif(3))
                if ((dif(1).eqv..false.).and.(dif(2).eqv..false.).and. &
                     (dif(3).eqv..false.)) then
                   inter%NumPnt(i)=inter%NumPnt(i)+1
                   inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(edge,1:2)
                   if((k==36).and.(i==695)) then
                      write(20,*),"numpnt",inter%NumPnt(i)
                      write(20,*),"adding nodesK",edge, nodesK(edge,1:2)
                   end if
                end if
             end if
             !print*,"adding point",x(1:2)


             call IntersectLineTri3b(nodesi(1,1:2),nodesi(3,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2), m)
             finpoint(1:2)=x(1:2)

             !print*,"inn",inn

             !print*,"adding point",x(1:2)

             !print*,"edge",edge

             !intersection of nodesi 1 & 2 is the edge number

             select case(edge)
             case(1)
                !4th point is on the K 1 & 2
                call SeekIntersect(nodesi(2,1:2),nodesi(3,1:2),nodesK(1,1:2),nodesK(2,1:2),x,inn)

             case(2)
                !4th point is on the K 2 & 3
                call SeekIntersect(nodesi(2,1:2),nodesi(3,1:2),nodesK(2,1:2),nodesK(3,1:2),x,inn)
             case(3)
                !4th point is on the K 1 & 3
                call SeekIntersect(nodesi(2,1:2),nodesi(3,1:2),nodesK(1,1:2),nodesK(3,1:2),x,inn)
             case default
                print*, "Error in A"
             end select
             !4th point
             !call Differ(x(1:2),inter%IntPnt(i,5,1:2),dif(1))
             !if (dif(1).eqv..false.) then
             !	inter%NumPnt(i) = inter%NumPnt(i)+1

             !	inter%IntPnt(i,inter%NumPnt(i)-1,1:2) = x(1:2)
             !	if((K==200).and.(i==171)) then

             !		write(20,*),"intersection with 2&3",l, inter%NumPnt(i)-1,x(1:2)
             !	end if
             !else
             difer=.false.
             n=0
             do l=1,3
                if (inner(l).eqv..true.) then
                   dif=.false.
                   do j=1,3
                      call Differ(xi(l,1:2),inter%IntPnt(i,j,1:2),dif(j))
                   end do
                   if ((dif(1).eqv..false.).and.(dif(2).eqv..false.).and.&
                        (dif(3).eqv..false.))then
                      !if (line(4).eqv..false.) then
                      !		if (difer.eqv..false.) then
                      n=n+1
                      xi2(n,1:2)=xi(l,1:2)
                      e(n)=l

                      !inter%NumPnt(i) = inter%NumPnt(i)+1
                      !inter%IntPnt(i,inter%NumPnt(i),1:2) = xi(l,1:2)
                      !inter%IntPnt(i,inter%NumPnt(i)-1,1:2) = xi(l,1:2)
                      if((k==36).and.(i==695)) then

                         write(20,*),"xi(l",l, inter%NumPnt(i),xi(l,1:2)
                      end if
                      !difer=.true.
                      !else
                      !inter%NumPnt(i) = inter%NumPnt(i)+1
                      !x(1:2)=inter%IntPnt(i,inter%NumPnt(i)-1,1:2)
                      !inter%IntPnt(i,inter%NumPnt(i)-1,1:2) = xi(l,1:2)
                      !inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                      !if((k==36).and.(i==695)) then

                      !		write(20,*),"xi(l",l, inter%NumPnt(i),xi(l,1:2)
                      !		end if
                      !		end if
                      !		else
                      !		inter%NumPnt(i) = inter%NumPnt(i)+1
                      !		inter%IntPnt(i,inter%NumPnt(i),1:2) = xi(l,1:2)
                      !end if


                   end if
                end if
             end do
             if (n==1) then
                inter%NumPnt(i) = inter%NumPnt(i)+1
                inter%IntPnt(i,inter%NumPnt(i),1:2) = xi2(1,1:2)
             else
                if (n==2) then
                   call FillInOrder(e(1:2),inter%NumPnt(i),inter%IntPnt(i,1:6,1:2),xi2(1,1:2),xi2(2,1:2))
                end if
             end if
             if ((k==36).and.(i==695)) then
                do l=1,inter%NumPnt(i)
                   write(20,*),l,inter%IntPnt(i,l,1:2)
                end do
             end if
             if (edge==3) then
                edge=1
             else
                edge=edge+1
             end if
             edgi(2)=edge
             !print*,"adding point",x(1:2)
             call Inside2(nodesK(edge,1:2),nodesi,line(1))
             !write(20,*),"inside edge",edge,line(1)
             if (line(1).eqv..true.) then
                dif=.false.
                do l=1,inter%NumPnt(i)
                   call Differ(nodesK(edge,1:2),inter%IntPnt(i,l,1:2),dif(l))
                end do
                call Differ(nodesK(edge,1:2),finpoint(1:2),dif(6))
                if ((dif(1).eqv..false.).and.(dif(2).eqv..false.).and. &
                     (dif(3).eqv..false.).and.(dif(5).eqv..false.) .and. &
                     (dif(6).eqv..false.)) then
                   inter%NumPnt(i)=inter%NumPnt(i)+1
                   inter%IntPnt(i,inter%NumPnt(i),1:2) = nodesK(edge,1:2)
                   !x(1:2)=inter%IntPnt(i,inter%NumPnt(i)-1,1:2)
                   !inter%IntPnt(i,inter%NumPnt(i)-1,1:2)=nodesK(edge,1:2)
                   !inter%IntPnt(i,inter%NumPnt(i),1:2)=x(1:2)
                   if((k==36).and.(i==695)) then
                      write(20,*),"addin nodesK",edge,nodesk(edge,1:2)
                   end if
                end if
             end if
             call GetThirdVertex2(edgi(1),edgi(2),vr,vr2,inn)
             if((k==36).and.(i==695)) then
                write(20,*),"inn,vr,vr2",vr,vr2,inn
             end if
             call Inside2(nodesK(vr,1:2),nodesi,line(1))
             if (line(1).eqv..true.) then
                dif=.false.
                do l=1,inter%NumPnt(i)
                   call Differ(nodesK(vr,1:2),inter%IntPnt(i,l,1:2),dif(l))
                end do
                call Differ(nodesK(edge,1:2),finpoint(1:2),dif(6))
                if ((dif(1).eqv..false.).and.(dif(2).eqv..false.).and. &
                     (dif(3).eqv..false.).and.(dif(5).eqv..false.) .and. &
                     (dif(6).eqv..false.)) then
                   !if (inter%NumPnt(i)==5) then
                   inter%NumPnt(i)=inter%NumPnt(i)+1
                   !x(1:2)=inter%IntPnt(i,inter%NumPnt(i)-1,1:2)
                   !inter%IntPnt(i,inter%NumPnt(i)-1,1:2)=nodesK(vr,1:2)
                   !inter%IntPnt(i,inter%NumPnt(i),1:2)=x(1:2)
                   inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(vr,1:2)
                   if((k==36).and.(i==695)) then

                      write(20,*),"nodesK(vr)",vr,inter%NumPnt(i),nodesK(vr,1:2)
                   end if
                   !else
                   !if (inter%NumPnt(i)<6) then
                   !inter%NumPnt(i)=inter%NumPnt(i)+1
                   !inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(vr,1:2)
                   !if((k==36).and.(i==695)) then

                   !write(20,*),"nodesK(vr)",vr,inter%NumPnt(i),nodesK(vr,1:2)
                   !end if
                   !else
                   !return
                   !end if
                   if((k==36).and.(i==695)) then

                      write(20,*),"nodesK(vr)",vr,inter%NumPnt(i),nodesK(vr,1:2)
                   end if
                   !write(20,*),"addin nodesK",edge,nodesk(edge,1:2)
                   !end if
                end if
             end if
             if (inn.eqv..true.) then
                call Inside2(nodesK(vr2,1:2),nodesi,line(1))
                if ((line(1).eqv..true.).and.(inter%NumPnt(i)<6)) then
                   dif=.false.
                   do l=1,inter%NumPnt(i)
                      call Differ(nodesK(vr2,1:2),inter%IntPnt(i,l,1:2),dif(l))
                   end do
                   call Differ(nodesK(edge,1:2),finpoint(1:2),dif(6))
                   if ((dif(1).eqv..false.).and.(dif(2).eqv..false.).and. &
                        (dif(3).eqv..false.).and.(dif(5).eqv..false.) .and. &
                        (dif(6).eqv..false.)) then
                      !if (inter%NumPnt(i)==5) then
                      !inter%NumPnt(i)=inter%NumPnt(i)+1
                      !x(1:2)=inter%IntPnt(i,inter%NumPnt(i)-1,1:2)
                      !inter%IntPnt(i,inter%NumPnt(i)-1,1:2)=nodesK(vr2,1:2)
                      !inter%IntPnt(i,inter%NumPnt(i),1:2)=x(1:2)
                      !if((k==36).and.(i==695)) then
                      !write(20,*),"nodesK(vr)",vr,inter%NumPnt(i),nodesK(vr2,1:2)
                      !end if
                      !else

                      if (inter%NumPnt(i)<6) then
                         inter%NumPnt(i)=inter%NumPnt(i)+1
                         inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(vr2,1:2)
                         if((k==36).and.(i==695)) then
                            write(20,*),"nodesK(vr)",vr,inter%NumPnt(i),nodesK(vr2,1:2)
                         end if
                         !else
                         !return
                      end if
                      !if((k==36).and.(i==695)) then

                      !write(20,*),"nodesK(vr)",vr,inter%NumPnt(i),nodesK(vr2,1:2)
                      !end if
                      !write(20,*),"addin nodesK",edge,nodesk(edge,1:2)
                      !end if
                   end if
                end if
             end if
             !5th point

             call IntersectLineTri3b(nodesi(1,1:2),nodesi(3,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2), m)
             !if (line(2).eqv..false.) then
             inter%NumPnt(i) = inter%NumPnt(i)+1

             inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
             if((k==36).and.(i==695)) then
                write(20,*),"numpnt",inter%NumPnt(i)
             end if
             !else
             !	inter%NumPnt(i) = inter%NumPnt(i)+1
             !
             !				inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
             !				if((k==36).and.(i==695)) then
             !				write(20,*),"numpnt",inter%NumPnt(i)
             !				end if
             !			end if
             if((k==36).and.(i==695)) then
                write(20,*),"inn",inn
                write(20,*),"intersection with 1&3",line(4), x(1:2)
                write(20,*),"edge",edge
             end if

          end select

          !number of intersection of nodesi 2 & 3 and K

       case default
          print*, "ERROR in IntersectTris - number of intersection of nodesi 2 & 3 is out of range"
       end select
       ! case of ins%num

    case(0)
       !print*,"edg%num + ver%num",edg%num + ver%num
       select case(edg%num + ver%num)
          ! case of edg%num + ver%num
       case(3)
          !print*,"edg%num+ver%num=3"
          inter%had(i) = .true.
          inter%NumPnt(i) = 3
          !print*,"inter%NumPnt(i)=3,,6"
          inter%IntPnt(i,1:3,1:2) = nodesi(1:3,1:2)
          !print*,"intpnt",inter%IntPnt(i,1:3,1:2)
          !print*,"return called"
          return
          ! case of edg%num + ver%num
          ! 2 on edge or vetrex
       case(2)
          select case(edg%num)
          case(2)
             !finding of the two nodes on the edge
             m=1
             do l=1,3
                if ( edg%is(l) .eqv. .true. ) then
                   edgi(m) = l
                   m = m+1
                end if
                if (ous%is(l) .eqv. .true.) then
                   ou=l
                end if
             end do
             !print*,"post",post(edgi(1)),post(edgi(2))
             !these two on the same edge?
             !YES---> just these two are the intersection
             if ( post(edgi(1)) == post(edgi(2)) ) then
                inter%had(i) = .true.
                inter%NumPnt(i) = 2
                inter%IntPnt(i,1,1:2) = nodesi(edgi(1),1:2)
                inter%IntPnt(i,2,1:2) = nodesi(edgi(2),1:2)
                !print*,"adding 2 edges ones",nodesi(edgi(1),1:2),nodesi(edgi(2),1:2)
                call IntersectLineTri3(nodesi(edgi(2),1:2),nodesi(ou,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
                call IntersectLineTri3(nodesi(edgi(1),1:2),nodesi(ou,1:2),nodesK,inn2,x2,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
                if (inn .eqv. .true.) then
                   inter%NumPnt(i)=inter%NumPnt(i)+1
                   inter%IntPnt(i,inter%NumPnt(i),1:2)=x(1:2)
                   !print*,"adding inter point of ",edgi(2),ou,x(1:2)
                end if
                if (inn2 .eqv. .true.) then
                   inter%NumPnt(i)=inter%NumPnt(i)+1
                   inter%IntPnt(i,inter%NumPnt(i),1:2)=x2(1:2)
                   !print*,"adding inter point of ",edgi(1),ou,x2(1:2)
                end if
                return
                ! NO--not on the same edge
             else
                !the outside vertex is number one
                !print*,"OneIsOne calling on ous and nodesi"
                call OneIsOne(ous,nodesi)
                inter%had(i) = .true.
                inter%NumPnt(i) = 2
                inter%IntPnt(i,1:2,1:2) = nodesi(2:3,1:2)
                !print*,"adding 2 edges ones",nodesi(2,1:2),nodesi(3,1:2)
                call IntersectLineTri3(nodesi(1,1:2),nodesi(3,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
                call IntersectLineTri3(nodesi(1,1:2),nodesi(2,1:2),nodesK,inn2,x2,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))

                !print*,"inn,inn2",inn,inn2
                if (inn .eqv. .true.) then
                   !2 intersection points, case U
                   if (inn2 .eqv. .true.) then
                      inter%NumPnt(i) = inter%NumPnt(i) + 1
                      inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                      inter%NumPnt(i) = inter%NumPnt(i) + 1
                      inter%IntPnt(i,inter%NumPnt(i),1:2) = x2(1:2)
                      !print*,"adding 2",x(1:2),x2(1:2)
                      return
                   else
                      !1 intersection point, case W1,X1
                      inter%NumPnt(i) = inter%NumPnt(i) + 1
                      inter%IntPnt(i,3,1:2) = x(1:2)
                      !print*,"adding inter",x(1:2)
                      call GetSharedVertex(post(edgi(1)),post(edgi(2)),vrt)
                      call Inside(nodesK(vrt,1:2),nodesi,line(1))
                      if (line(1) .eqv. .true.) then
                         inter%NumPnt(i)=inter%NumPnt(i)+1
                         inter%IntPnt(i,inter%NumPnt(i),1:2) = nodesK(vrt,1:2)
                         !print*,"adding nodesK",nodesK(vrt,1:2)
                      end if
                      !print*,"adding nodesK",nodesK(vrt,1:2)
                      return
                   end if
                else
                   !inner(1) = .false.
                   !W2,X2
                   if (inn2 .eqv. .true.) then
                      !inter%NumPnt(i)=4
                      call GetSharedVertex(post(edgi(1)),post(edgi(2)),vrt)
                      call Inside(nodesK(vrt,1:2),nodesi,line(1))
                      if (line(1) .eqv. .true.) then
                         inter%NumPnt(i)=inter%NumPnt(i)+1
                         inter%IntPnt(i,inter%NumPnt(i),1:2) = nodesK(vrt,1:2)
                         !print*,"adding nodesK",nodesK(vrt,1:2)
                      end if
                      inter%NumPnt(i)=inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2) = x2(1:2)
                      !print*,"adding inter",x2(1:2)
                      return
                   else
                      !V1 i V2
                      call GetSharedVertex(post(edgi(1)),post(edgi(2)),vrt)
                      !print*,"inter%NumPnt(i)=3"
                      inter%NumPnt(i)=3
                      inter%IntPnt(i,3,1:2)=nodesK(vrt,1:2)
                      !print*,"adding 1,shared vertex",nodesK(vrt,1:2)
                      return
                   end if

                end if

             end if

            !case of edg%num
          case(1)
             do j=1,3

               if (ous%is(j) .eqv. .true.) then
                   ou=j
                end if
                if (edg%is(j) .eqv. .true.) then
                   ed=j
                end if
                if (ver%is(j) .eqv. .true.) then
                   vr=j
                end if
             end do
             !add the first two
             inter%had(i)=.true.
             inter%NumPnt(i)=2
             inter%IntPnt(i,1,1:2)=nodesi(ed,1:2)
             inter%IntPnt(i,2,1:2)=nodesi(vr,1:2)
             !print*,"adding edge and vertex"

             call IntersectLineTri3(nodesi(ou,1:2),nodesi(ed,1:2),nodesK,inner(1),xi(1,1:2),inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
             call IntersectLineTri3(nodesi(ou,1:2),nodesi(vr,1:2),nodesK,inner(2),xi(2,1:2),inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
             !((2))
             if (inner(1) .eqv. .true.) then
                if (inner(2) .eqv. .false.) then
                   inter%NumPnt(i)=inter%NumPnt(i)+1
                   inter%IntPnt(i,inter%NumPnt(i),1:2)=xi(1,1:2)
                   return
                else
                  !print*, "Error, in ((2))"
                   inter%NumPnt(i)=4
                   inter%IntPnt(i,3,1:2)=xi(2,1:2)
                   inter%IntPnt(i,4,1:2)=xi(1,1:2)
                   return
                end if
             else

                !((3))
                if (inner(2) .eqv. .true.) then
                   inter%NumPnt(i)=inter%NumPnt(i)+1
                   inter%IntPnt(i,inter%NumPnt(i),1:2)=xi(2,1:2)
                   return
                else
                   do l=1,3
                      call Inside2(nodesK(l,1:2),nodesi,line(l))
                      if (line(l) .eqv. .true.) then
                         call Differ(inter%IntPnt(i,1,1:2),nodesK(l,1:2),dif(1))
                         call Differ(inter%IntPnt(i,2,1:2),nodesK(l,1:2),dif(2))
                         if ((dif(1).eqv..false.).and. (dif(2).eqv..false.))then
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(l,1:2)
                            return
                         end if
                      end if
                   end do
                end if
             end if

             !case of edg%num
             !2 vertices, 1 outside
          case(0)
             !add 2 vertices
             inter%had(i)=.true.
             !print*, "2vertices,1 outside"
             !print*,ver%is(1:3)
             !print*,inter%NumPnt(i)
             do l=1,3
                if ( ver%is(l) .eqv. .true. ) then
                   !print*,"adding of",l,"vertex", nodesi(l,1:2)
                   inter%NumPnt(i) = inter%NumPnt(i)+1
                   inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesi(l,1:2)
                end if
             end do
             !print*,"2 vertices added"
             do l=1,3
                if (ous%is(l) .eqv. .true.) then
                   call SwapNodes(nodesi,l)
                   exit
                end if
             end do
             !print*,l,"the one outside"
             !has it intersection or not?
             !2 or 3 points
             call IntersectLineTri3(nodesi(1,1:2),nodesi(2,1:2),nodesK,inner(1),xi(1,1:2),inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
             if (inner(1) .eqv. .true.) then
                inter%NumPnt(i) = inter%NumPnt(i)+1
                inter%IntPnt(i,inter%NumPnt(i),1:2)=xi(1,1:2)
                !print*,"inter%NumPnt(i)",inter%NumPnt(i)
                !print*,"inter%IntPnt(i,inter%NumPnt(i),1:2)",inter%IntPnt(i,inter%NumPnt(i),1:2)
                !print*,"third intersection",xi(1,1:2)
                return
             end if
             call IntersectLineTri3(nodesi(1,1:2),nodesi(3,1:2),nodesK,inner(1),xi(1,1:2),inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
             if (inner(1) .eqv. .true.) then
                inter%NumPnt(i) = inter%NumPnt(i)+1
                inter%IntPnt(i,inter%NumPnt(i),1:2)=xi(1,1:2)
                !print*,"third intersection-2nd case",xi(1,1:2)
                return
             end if

             do l=1,3
                call Inside(nodesK(l,1:2),nodesi,line(1))
                !print*,"Inside on",l,line(1)
                if (line(1) .eqv. .true.) then
                   call Differ(nodesK(l,1:2),inter%IntPnt(i,1,1:2),dif(1))
                   call Differ(nodesK(l,1:2),inter%IntPnt(i,2,1:2),dif(2))
                   if ((dif(1).eqv..false.).and.(dif(2).eqv..false.))then
                      inter%NumPnt(i)=inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(l,1:2)
                      !print*,"nodesK",l,"added"
                      !print*,nodesK(l,1:2)
                      return
                   end if
                end if
             end do
             return
          end select
          ! case of edg%num + ver%num
          !1 edge or vertex, 2 outside
       case(1)
          !print*,"1 edge or vertex, 2 outside"
          !print*,"edg%num",edg%num
          !(edg%num == 1)
          if (edg%num == 1) then
             j=1
             do
                if (edg%is(j) .eqv. .true.) then
                   posi=post(j)
                   exit
                end if
                j=j+1
             end do
             !if j>3 then some error must have occured, cause we have just 3 vertices
             if (j>3) then
                print*, "ERROR in IntersectTris"
             end if
             !if j isn't one then we should swap nodes as if the one outside was node number one
             if ( j/=1 ) then
                call SwapNodes(nodesi,j)
                !print*,"SwapNodes called on nodesi and",j
             end if
             !adding the one on the edge
             inter%had(i)=.true.
             inter%NumPnt(i)=inter%NumPnt(i)+1
             inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesi(1,1:2)
             !print*,"adding the one on the edge",nodesi(1,1:2)

             !intersection of outside vertices and K
             call IntersectLineTri2(nodesi(2,1:2),nodesi(3,1:2),nodesK,inner,xi,cnt)
             !print*,"number of intersection2&3 and nodesK",cnt
             !select for number of intersection 2&3 and K
             select case(cnt)
                !number of intersection 2&3 and K
             case(2)
                !2
                call IntersectLineTri3b(nodesi(1,1:2),nodesi(2,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2),edge)
                !print*,"edge",edge
                !((2a))
                if (inn .eqv. .true.) then
                   inter%NumPnt(i) = inter%NumPnt(i)+1
                   inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                   !print*,"adding intersection of 1&2",x(1:2)
                   !intersection of nodesi 1 & 2 is the edge number
                   !print*,"edge",edge
                   select case(edge)
                   case(1)
                      !3rd point is on the K 1 & 2
                      call SeekIntersect(nodesi(2,1:2),nodesi(3,1:2),nodesK(1,1:2),nodesK(2,1:2),x,inn)
                   case(2)
                      !3rd point is on the K 2 & 3
                      call SeekIntersect(nodesi(2,1:2),nodesi(3,1:2),nodesK(2,1:2),nodesK(3,1:2),x,inn)
                   case(3)
                      !3rd point is on the K 1 & 3
                      call SeekIntersect(nodesi(2,1:2),nodesi(3,1:2),nodesK(1,1:2),nodesK(3,1:2),x,inn)
                   case default
                      print*, "Error in A"
                   end select
                   !3
                   inter%NumPnt(i) = inter%NumPnt(i)+1
                   inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                   !print*,"adding intersection with edge",x(1:2)


                   call IntersectLineTri3(nodesi(2,1:2),nodesi(3,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
                   if (inn .eqv. .true.) then
                      inter%NumPnt(i) = inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                      !print*,"adding intersection with 2&3",x(1:2)
                   end if
                   !5
                   call IntersectLineTri3(nodesi(1,1:2),nodesi(3,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
                   if (inn .eqv. .true.) then
                      inter%NumPnt(i) = inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                      !print*,"adding intersection with 1&3",x(1:2)
                   end if
                   return

                   !25/4/((2b))
                else
                   call IntersectLineTri3(nodesi(1,1:2),nodesi(3,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
                   inter%NumPnt(i) = inter%NumPnt(i)+1
                   inter%IntPnt(i,inter%NumPnt(i)+1,1:2) = x(1:2)
                  m=1
                  do j=1,3
                      if (inner(j) .eqv. .true.) then
                        edgi(m)=j
                        xi2(m,1:2)=xi(j,1:2)
                        m=m+1
                        !inter%NumPnt(i) = inter%NumPnt(i)+1
                         !inter%IntPnt(i,inter%NumPnt(i),1:2) = xi(j,1:2)
                      end if
                   end do
                  call FillInOrder(edgi(1:2),inter%NumPnt(i),inter%IntPnt(i,1:6,1:2),xi2(1,1:2),xi2(2,1:2))
                  return
                end if

               !number of intersection 2&3 and K
             case(1)
                call IntersectLineTri3(nodesi(1,1:2),nodesi(2,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
                call IntersectLineTri3(nodesi(1,1:2),nodesi(3,1:2),nodesK,inn2,x2,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
                if (inn .eqv. .true.) then
                   !add intersection with 1&2
                   inter%NumPnt(i) = inter%NumPnt(i)+1
                   inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                   !add intersection with 2&3
                   do j=1,3
                      if (inner(j) .eqv. .true.) then
                         inter%NumPnt(i) = inter%NumPnt(i)+1
                         inter%IntPnt(i,inter%NumPnt(i),1:2) = xi(j,1:2)
                      end if
                   end do
                   !add intersection with 1&3
                   if (inn2 .eqv. .true.) then
                      inter%NumPnt(i) = inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2) = x2(1:2)
                      return
                   end if
                else
                   !should be needless,there should allways be such intersection
                   if (inn2 .eqv. .true.) then
                      inter%NumPnt(i) = inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2) = x2(1:2)
                      do j=1,3
                         if (inner(j) .eqv. .true.) then
                            inter%NumPnt(i) = inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i),1:2) = xi(j,1:2)
                         end if
                      end do
                   end if
                end if

                !number of intersection 2&3 and K
             case(0)
                call IntersectLineTri3b(nodesi(1,1:2),nodesi(2,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2),edgi(1))
                !print*,"Intersect of nodesi 1&2 and nodesK",inn

                call IntersectLineTri3b(nodesi(1,1:2),nodesi(3,1:2),nodesK,inn2,x2,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2),edgi(2))
                !if (K==287) then
                !open(17,file="edg")
                !		write(17,*),"K",K
                !		write(17,*),"edgi",edgi(1:2)
                !		write(17,*),"inn",inn,inn2
                !		write(17,*),"x",x(1:2)
                !		write(17,*),"x2",x2(1:2)
                !	close(17)
                !	end if
                !print*,"Intersect of nodesi 1&3 and nodesK",inn2

                !print*,"edgi1,2",edgi(1),edgi(2)

                if (inn .eqv. .true.) then
                   inter%NumPnt(i) = inter%NumPnt(i)+1
                   inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                   !print*,"adding of",x(1:2)
                   if (inn2 .eqv. .true.) then
                      if (edgi(1)/=edgi(2)) then
                         !((0.a))
                         call GetSharedVertex2(edgi(1),edgi(2),vrt)
                         !print*,"vrt",vrt
                         inter%NumPnt(i) = inter%NumPnt(i)+1
                         inter%IntPnt(i,inter%NumPnt(i),1:2) = nodesK(vrt,1:2)
                         inter%NumPnt(i) = inter%NumPnt(i)+1
                         inter%IntPnt(i,inter%NumPnt(i),1:2) = x2(1:2)
                         return
                         !((0.b))
                      else
                         inter%NumPnt(i) = inter%NumPnt(i)+1
                         inter%IntPnt(i,inter%NumPnt(i),1:2) = x2(1:2)
                        !write(17,*),"tady"
                         return
                      end if
                   else
                      !((0.c))

                      !open(18,file="edgi")
                      !write(18,*),"K",K
                      !write(18,*),"edgi",edgi(1:2)
                      !write(18,*),"inn",inn,inn2
                      !write(18,*),"x",x(1:2)
                      !write(18,*),"x2",x2(1:2)
                      !adding these nodesK which are inside triangle

                      do j=1,3

                         call Inside2(nodesK(j,1:2),nodesi,line(j))
                         !write(17,*),"line(j)",j,line(j)
                         !nodesK(j,1:2) is inside nodesi

                         if (line(j) .eqv. .true.) then
                            !write(17,*),j,"nodesk is inside",nodesK(j,1:2)
                            dif(1:3) = .false.

                            !is he also different from already added ones?

                            do l=1,inter%NumPnt(i)

                               call Differ(nodesK(j,1:2),inter%IntPnt(i,l,1:2),dif(l))

                            end do
                            !write(18,*),"Differ",dif(1:3)
                            if ((dif(1) .eqv. .false.) .and. (dif(2) .eqv. .false.) .and. (dif(3) .eqv. .false.)) then

                               inter%NumPnt(i) = inter%NumPnt(i)+1

                               inter%IntPnt(i,inter%NumPnt(i),1:2) = nodesK(j,1:2)
                               !write(17,*),"adding nodesK",j,nodesK(j,1:2)
                               !print*,"adding nodesK",j,nodesK(j,1:2)

                            end if

                         end if

                      end do
                      !close(17)
                      !close(18)
                      return

                   end if

                else

                   !((0.d))

                   if (inn2 .eqv. .true.) then
                      !add intersect 1&3 and K
                      inter%NumPnt(i) = inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2) = x2(1:2)
                      !print*,"adding intersection with 1&3",x2(1:2)

                      do j=1,3
                         call Inside2(nodesK(j,1:2),nodesi,line(j))
                         !nodesK(j,1:2) is inside nodesi
                         if (line(j) .eqv. .true.) then
                            dif(1:3) = .false.
                            !is he also different from already added ones?
                            do l=1,inter%NumPnt(i)
                               call Differ(nodesK(j,1:2),inter%IntPnt(i,l,1:2),dif(l))
                            end do
                            if ((dif(1) .eqv. .false.) .and. (dif(2) .eqv. .false.) .and. (dif(3) .eqv. .false.)) then
                               inter%NumPnt(i) = inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i),1:2) = nodesK(j,1:2)
                               !print*,"adding nodesK",j,nodesK(j,1:2)
                            end if
                         end if
                      end do
                      return
                   else
                      !getting know if we shall add something or not(case ((0.e)) or ((0.f)))
                      l=1
                      do
                         !print*,"call InLine"
                         call InLine(nodesK(l,1:2),nodesi(1,1:2),nodesi(2,1:2),line(1))
                         !print*,"line",line(1)
                         if (line(1) .eqv. .true.) then
                            !print*,"Inline nodesi 1&2 and nodesK",l
                            call Differ(nodesK(l,1:2),inter%IntPnt(i,1,1:2),dif(1))
                            if (dif(1) .eqv. .false.) then
                               inter%NumPnt(i) = inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i),1:2) = nodesK(l,1:2)
                            end if
                            return
                         end if
                         call InLine(nodesK(l,1:2),nodesi(1,1:2),nodesi(3,1:2),line(2))
                         !print*,"line",line(1)
                         if (line(2) .eqv. .true.) then
                            !print*,"Inline nodesi 1&3 and nodesK",l
                            call Differ(nodesK(l,1:2),inter%IntPnt(i,1,1:2),dif(1))
                            if (dif(1) .eqv. .false.) then
                               inter%NumPnt(i) = inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i),1:2) = nodesK(l,1:2)
                            end if
                            return
                         end if
                         if (l==3) then
                            exit
                            return
                         end if
                         l=l+1
                      end do
                   end if
                end if
             end select

             !!(edg%num /= 1)-->ver%num=1
             !! 23.5.
          else
             j=1
             do
                if (ver%is(j) .eqv. .true.) then
                   !print*,j,"th vertex is on vertex"
                   posi=post(j)
                   exit
                end if
                j=j+1
             end do
             !if j>3 then some error must have occured, cause we have just 3 vertices
             if (j>3) then
                print*, "ERROR in IntersectTris"
             end if
             !if j isn't one then we should swap nodes as if the one outside was node number one
             if ( j/=1 ) then
                !print*,"SwapNodes called on",j,"vertex"
                call SwapNodes(nodesi,j)
             end if
             !adding the one on the vertex
             inter%had(i)=.true.
             inter%NumPnt(i)=inter%NumPnt(i)+1
             inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesi(1,1:2)
             !print*,"adding nodesi(1,1:2)",nodesi(1,1:2)
             !if ((K==287) .and.(i==342)) then
             !	open(31,file="spec")
             !	write(31,*),"nodesi",inter%NumPnt(i),nodesi(1,1:2)
             !end if


             !intersection of outside vertices and K

             call IntersectLineTri2(nodesi(2,1:2),nodesi(3,1:2),nodesK,inner,xi,cnt)

             !print*,"count of intersection of nodesi 2&3 and nodesK",cnt

             !!print*,"nodesi(2,1:2)",nodesi(2,1:2)

             !!print*,"nodesi(3,1:2)",nodesi(3,1:2)

             !select for number of intersection 2&3 and K
             !if ((K==287) .and.(i==342)) then
             !open(31,file="spec")
             !	write(31,*),"cnt",cnt
             !	do j=1,3
             !		if (inner(j).eqv..true.) then
             !			write(31,*),xi(j,1:2)
             !		end if
             !	end do
             !end if
             select case(cnt)

                !number of intersection 2&3 and K

             case(2)

                !2

                call IntersectLineTri3b(nodesi(1,1:2),nodesi(2,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2),edge)
                !if ((K==287) .and.(i==342)) then
                !write(31,*),"inn",inn
                !write(31,*),"edge",edge
                !write(31,*),"x",x(1:2)
                !end if
                !((2a))

                if (inn .eqv. .true.) then

                   inter%NumPnt(i) = inter%NumPnt(i)+1

                   inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                   !if ((K==287) .and.(i==342)) then
                   !	write(31,*),"x",inter%NumPnt(i),x(1:2)
                   !end if
                   !intersection of nodesi 1 & 2 is the edge number

                   select case(edge)

                   case(1)

                      !3rd point is on the K 1 & 2

                      call SeekIntersect(nodesi(2,1:2),nodesi(3,1:2),nodesK(1,1:2),nodesK(2,1:2),x,inn)
                      !if ((K==287) .and.(i==342)) then
                      !write(31,*),"inn",inn
                      !write(31,*),"x",x(1:2)
                      !end if


                   case(2)
                      !3rd point is on the K 2 & 3
                      call SeekIntersect(nodesi(2,1:2),nodesi(3,1:2),nodesK(2,1:2),nodesK(3,1:2),x,inn)

                   case(3)
                      !3rd point is on the K 1 & 3
                      call SeekIntersect(nodesi(2,1:2),nodesi(3,1:2),nodesK(1,1:2),nodesK(3,1:2),x,inn)

                   case default
                      print*, "Error in A"
                   end select
                   if (inn .eqv. .true.) then
                      inter%NumPnt(i) = inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                      !if ((K==287) .and.(i==342)) then
                      !write(31,*),"inn",inn
                      !write(31,*),"x",inter%NumPnt(i),x(1:2)
                      !end if
                   end if

                   !5
                   call IntersectLineTri3(nodesi(1,1:2),nodesi(3,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
                   !write(31,*),"1&3",K
                   if (inn.eqv..true.) then
                      inter%NumPnt(i) = inter%NumPnt(i)+1

                      inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                      !	if ((K==287) .and.(i==342)) then
                      !write(31,*),"inn",inn
                      !	write(31,*),"x",inter%NumPnt(i),x(1:2)
                      !	end if
                   end if
                   if (inn.eqv..false.) then
                      !call IntersectLineTri2(nodesi(2,1:2),nodesi(3,1:2),nodesK,inner,xi,cnt)
                      do l=1,3
                         if (inner(l).eqv..true.) then
                            call Differ(xi(l,1:2),inter%IntPnt(i,1,1:2),dif(1))
                            call Differ(xi(l,1:2),inter%IntPnt(i,2,1:2),dif(2))
                            call Differ(xi(l,1:2),inter%IntPnt(i,3,1:2),dif(3))
                            !write(31,*),"l,dif",l,dif(1:3)
                            if ((dif(1).eqv..false.).and.(dif(2).eqv..false.)&
                                 .and.(dif(3).eqv..false.)) then
                               inter%NumPnt(i) = inter%NumPnt(i)+1

                               inter%IntPnt(i,inter%NumPnt(i),1:2) = xi(l,1:2)
                               !	if ((K==287) .and.(i==342)) then
                               !write(31,*),"inn",inn
                               !	write(31,*),"x",inter%NumPnt(i),xi(l,1:2)
                               !	end if
                               return
                            end if
                         end if

                      end do
                   end if

                   !((2b))

                else
                   call IntersectLineTri3b(nodesi(1,1:2),nodesi(3,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2),edge)
                   inter%NumPnt(i) = inter%NumPnt(i)+1
                   inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                   select case(edge)
                   case(1)
                      !3rd point is on the K 1 & 2
                      call SeekIntersect(nodesi(2,1:2),nodesi(3,1:2),nodesK(1,1:2),nodesK(2,1:2),x,inn)

                   case(2)
                      !3rd point is on the K 2 & 3
                      call SeekIntersect(nodesi(2,1:2),nodesi(3,1:2),nodesK(2,1:2),nodesK(3,1:2),x,inn)

                   case(3)
                      !3rd point is on the K 1 & 3
                      call SeekIntersect(nodesi(2,1:2),nodesi(3,1:2),nodesK(1,1:2),nodesK(3,1:2),x,inn)

                   case default
                      print*, "Error in A"
                   end select
                   if (inn .eqv. .true.) then
                      inter%NumPnt(i) = inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                   end if
                   call IntersectLineTri3(nodesi(1,1:2),nodesi(3,1:2),nodesK,inn,x,inter%NumPnt(i),inter%IntPnt(i,1:6,1:2))
                   inter%NumPnt(i) = inter%NumPnt(i)+1
                   inter%IntPnt(i,inter%NumPnt(i),1:2) = x(1:2)
                end if

                !number of intersection 2&3 and K
             case(1)
                do l=1,3
                   call Inside2(nodesK(l,1:2),nodesi,inn)
                   !print*,"nodesK",l,"inside",inn
                   call Differ(nodesK(l,1:2),nodesi(1,1:2),dif(l))
                   !print*,"nodesK",l,"differ",dif(l)
                   if ((dif(l) .eqv. .false.) .and. (inn .eqv. .true.))then
                      !print*,"adding of nodesK(l,1:2)",nodesK(l,1:2)
                      inter%NumPnt(i) = inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2) = nodesK(l,1:2)
                   end if
                end do
                !!23/5/1.s
                ! if all the nodesK aren't inside then we should find third intersection
                if (inter%NumPnt(i)<3) then
                   call IntersectLineTri2(nodesi(1,1:2),nodesi(2,1:2),nodesK,inner,xi,cnt)
                   !print*,"cnt of nodesi 1&2",cnt
                   if (cnt>1) then
                      !print*,"cnt>1"
                      do l=1,3
                         !print*,"inner",l,"::",inner(l)
                         if (inner(l) .eqv. .true.) then
                            call Differ(xi(l,1:2),nodesi(1,1:2),difer)
                            !print*,"difer",difer
                            if (difer .eqv. .false.) then
                               !print*,"adding of",xi(l,1:2)
                               inter%NumPnt(i) = inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i),1:2) = xi(l,1:2)
                               return
                            end if
                         end if
                      end do
                   end if
                   call IntersectLineTri2(nodesi(1,1:2),nodesi(3,1:2),nodesK,inner,xi,cnt)
                   !print*,"cnt of nodesi 1&2",cnt
                   if (cnt>1) then
                      !print*,"cnt2>1"
                      do l=1,3
                         if (inner(l) .eqv. .true.) then
                            !print*,"inner",l,"::",inner(l)
                            call Differ(xi(l,1:2),nodesi(1,1:2),difer)
                            !print*,"difer",difer
                            if (difer .eqv. .false.) then
                               !print*,"adding of-2nd",xi(l,1:2)
                               inter%NumPnt(i) = inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i),1:2) = xi(l,1:2)
                               return
                            end if
                         end if
                      end do
                   end if
                end if

                !number of intersection 2&3 and K
             case(0)
                call Inside2(nodesK(1,1:2),nodesi(1:3,1:2),inner(1))
                !print*,"nodesK(1,1:2) is inside nodesi",inner(1)
                call Inside2(nodesK(2,1:2),nodesi(1:3,1:2),inner(2))
                !print*,"nodesK(2,1:2) is inside nodesi",inner(2)
                call Inside2(nodesK(3,1:2),nodesi(1:3,1:2),inner(3))
                !print*,"nodesK(3,1:2) is inside nodesi",inner(3)
                do l=1,3
                   call Differ(nodesK(l,1:2),nodesi(1,1:2),dif(l))
                   ! if the points inside and it's different from the one which we have
                   ! already added then add him
                   if ((inner(l) .eqv. .true.) .and. (dif(l) .eqv. .false.)) then
                      inter%NumPnt(i) = inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2) = nodesK(l,1:2)
                   end if
                end do
                if (posi=="VB") then
                   call SwapNodes(nodesK,2)
                   !print*,"SwapNodes called on 2 nodesK"
                end if
                if (posi=="VC") then
                   call SwapNodes(nodesK,3)
                   !print*,"SwapNodes called on 3 nodesK"
                end if
                call IntersectLineTri2(nodesK(2,1:2),nodesK(3,1:2),nodesi,inner,xi,cnt)
                !print*,"cnt",cnt
                do l=1,3
                   if (inner(l) .eqv. .true.) then
                      !print*,"xi",l,":",xi(l,1:2)
                   end if
                end do
                do l=1,3
                   if (inner(l) .eqv. .true.) then
                      dif=.false.
                      do j=1,inter%NumPnt(i)
                         call Differ(xi(l,1:2),inter%IntPnt(i,j,1:2),dif(j))
                      end do
                      if ((dif(1) .eqv. .false.).and.(dif(2).eqv..false.).and.(dif(3).eqv..false.)&
                           .and.(dif(4).eqv..false.)) then
                         inter%NumPnt(i) = inter%NumPnt(i)+1
                         inter%IntPnt(i,inter%NumPnt(i),1:2) = xi(l,1:2)
                      end if
                   end if
                end do
                return

             end select
             !!(edg%num == 1)
          end if
          ! case of edg%num + ver%num
       case(0)
          call IntersectLineTri2(nodesi(1,1:2),nodesi(2,1:2),nodesK,inner3(1,1:3),xi3(1,1:3,1:2),cnt)
          call IntersectLineTri2(nodesi(2,1:2),nodesi(3,1:2),nodesK,inner3(2,1:3),xi3(2,1:3,1:2),cnt2)
          call IntersectLineTri2(nodesi(1,1:2),nodesi(3,1:2),nodesK,inner3(3,1:3),xi3(3,1:3,1:2),cnt3)
          !if ((K==287).and.(i==343)) then
          !open(19,file="ous")
          !write(19,*),"K",K
          !write(19,*),"cnt,cnt2,cnt3",cnt,cnt2,cnt3
          !end if
          !print*,"inner",inner3(1,1:3),"xi",xi3(1,1:3,1:2)

          !print*,"cnt",cnt

          !print*,"inner",inner3(2,1:3),"xi",xi3(2,1:3,1:2)

          !print*,"cnt2",cnt2

          !print*,"inner",inner3(3,1:3),"xi",xi3(3,1:3,1:2)

          !print*,"cnt3",cnt3

          !select according to number of intersections points in these two triangles

          select case(cnt+cnt2+cnt3)

             !cnt+cnt2+cnt3
          case(6)
             inter%had(i)=.true.
             !adding six intersection in right order
             m=1
             do j=1,3
                if (inner3(1,j).eqv..true.) then
                   edgi(m)=j
                   m=m+1
                   inner3(1,j)=.false.
                end if
             end do
             !print*,"edgi",edgi(1:2)
             if((edgi(1)==2).or.(edgi(2)==2))then
                if (edgi(1)<edgi(2)) then
                   inter%NumPnt(i)=2
                   inter%IntPnt(i,1,1:2)=xi3(1,edgi(1),1:2)
                   inter%IntPnt(i,2,1:2)=xi3(1,edgi(2),1:2)
                   ed=edgi(2)
                else
                   inter%NumPnt(i)=2
                   inter%IntPnt(i,1,1:2)=xi3(1,edgi(2),1:2)
                   inter%IntPnt(i,2,1:2)=xi3(1,edgi(1),1:2)
                   ed=edgi(1)
                end if
             else
                if (edgi(1)==3) then
                   inter%NumPnt(i)=2
                   inter%IntPnt(i,1,1:2)=xi3(1,edgi(1),1:2)
                   inter%IntPnt(i,2,1:2)=xi3(1,edgi(2),1:2)
                   ed=edgi(2)
                else
                   inter%NumPnt(i)=2
                   inter%IntPnt(i,1,1:2)=xi3(1,edgi(2),1:2)
                   inter%IntPnt(i,2,1:2)=xi3(1,edgi(1),1:2)
                   ed=edgi(1)
                end if
             end if

             do l=2,3
                do j=1,3
                   if (inner3(l,j).eqv..true.) then
                      !print*,"l,j",l,j
                      if ((inter%NumPnt(i)==2).or.(inter%NumPnt(i)==4)) then
                         if (ed==j) then
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(l,j,1:2)
                            inner3(l,j)=.false.
                            inn=.false.
                            !print*,"added"
                         else
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i)+1,1:2)=xi3(l,j,1:2)
                            inner3(l,j)=.false.
                            inn=.true.
                            ed=j
                            !print*,"added +1"
                         end if
                      else
                         if (inn.eqv..false.) then
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(l,j,1:2)
                            inner3(l,j)=.false.
                            ed=j
                            !print*,"added"
                         else
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i)-1,1:2)=xi3(l,j,1:2)
                            inner3(l,j)=.false.
                            !print*,"added -1"
                         end if
                      end if
                   end if
                end do
             end do
          case(5)
             !adding 5 intersection points in right order
             !firstly add the intersection of the edge which has just one intersection
             !then add the others edges in right order
             inter%had(i)=.true.
             if (cnt==1) then
                j=1
                do
                   if (inner3(1,j).eqv..true.) then

                     inter%NumPnt(i)=inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(1,j,1:2)
                      inner3(1,j)=.false.
                      l=1
                      do
                         call Differ(nodesK(l,1:2),xi3(1,j,1:2),dif(1))
                         if (dif(1).eqv..true.) then
                            exit
                         end if
                         l=l+1
                      end do
                      exit
                   end if
                   j=j+1
                end do
                !l is next edge then the one last added
                inn=.false.

                do j=1,3
                   if (inner3(2,j).eqv..true.) then
                      if (inter%NumPnt(i)==1) then
                         if (l==j) then
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(2,j,1:2)
                            inner3(2,j)=.false.
                            inn=.false.
                         else
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i)+1,1:2)=xi3(2,j,1:2)
                            inner3(2,j)=.false.
                            inn=.true.
                            ed=j
                         end if
                      else
                         if (inn.eqv..false.) then
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(2,j,1:2)
                            inner3(2,j)=.false.
                            ed=j
                         else
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i)-1,1:2)=xi3(2,j,1:2)
                            inner3(2,j)=.false.
                         end if
                      end if
                   end if
                end do

                do j=1,3
                   if (inner3(3,j).eqv..true.) then
                      if (inter%NumPnt(i)==3) then
                         if (ed==j) then
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(3,j,1:2)
                            inner3(3,j)=.false.
                            inn=.false.
                         else
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i)+1,1:2)=xi3(3,j,1:2)
                            inner3(3,j)=.false.
                            inn=.true.
                         end if
                      else
                         if (inn.eqv..false.) then
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(3,j,1:2)

                           inner3(3,j)=.false.
                         else
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i)-1,1:2)=xi3(3,j,1:2)
                            inner3(3,j)=.false.
                         end if
                      end if
                   end if
                end do
             else
                if (cnt2==1) then
                   j=1
                   do
                      if (inner3(2,j).eqv..true.) then
                         inter%NumPnt(i)=inter%NumPnt(i)+1
                         inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(2,j,1:2)
                         inner3(2,j)=.false.
                         l=1
                         do
                            call Differ(nodesK(l,1:2),xi3(2,j,1:2),dif(1))
                            !print*,"dif(1),l",dif(1),l
                            if (dif(1).eqv..true.) then
                               exit
                            end if
                            l=l+1
                         end do
                         exit
                      end if
                      j=j+1
                   end do
                   !l is next edge then the one last added
                   inn=.false.
                   !print*,"l",l
                   do j=1,3
                      if (inner3(3,j).eqv..true.) then
                         !print*,"3 and",j
                         if (inter%NumPnt(i)==1) then
                            if (l==j) then
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(3,j,1:2)
                               inner3(3,j)=.false.
                               inn=.false.
                               !print*,"added",inter%NumPnt(i)
                            else
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i)+1,1:2)=xi3(3,j,1:2)
                               inner3(3,j)=.false.
                               inn=.true.
                               ed=j
                               !print*,"added to +1,ed",ed,inter%NumPnt(i)

                            end if
                         else
                            if (inn.eqv..false.) then
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(3,j,1:2)
                               inner3(3,j)=.false.
                               ed=j
                               !print*,"added",inter%NumPnt(i),"ed",ed
                            else
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i)-1,1:2)=xi3(3,j,1:2)
                               inner3(3,j)=.false.
                               !print*,"added -1",inter%NumPnt(i)
                            end if
                         end if
                      end if
                   end do
                   !print*,"ed",ed
                   do j=1,3
                      if (inner3(1,j).eqv..true.) then
                         !print*,"1 and",j
                         if (inter%NumPnt(i)==3) then
                            if (ed==j) then
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(1,j,1:2)
                               inner3(1,j)=.false.
                               inn=.false.
                               !print*,"added",inter%NumPnt(i)
                            else
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i)+1,1:2)=xi3(1,j,1:2)
                               inner3(1,j)=.false.
                               inn=.true.
                               !print*,"added +1",inter%NumPnt(i)
                            end if
                         else
                            if (inn.eqv..false.) then
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(1,j,1:2)
                               inner3(1,j)=.false.
                               !print*,"added",inter%NumPnt(i)
                            else
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i)-1,1:2)=xi3(1,j,1:2)
                               inner3(1,j)=.false.
                               !print*,"added -1",inter%NumPnt(i)
                            end if
                         end if
                      end if
                   end do
                else
                   j=1
                   do
                      if (inner3(3,j).eqv..true.) then
                         inter%NumPnt(i)=inter%NumPnt(i)+1
                         inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(3,j,1:2)
                         inner3(3,j)=.false.
                         l=1
                         do
                            call Differ(nodesK(l,1:2),xi3(3,j,1:2),dif(1))
                            if (dif(1).eqv..true.) then
                               exit
                            end if
                            l=l+1
                         end do
                         exit
                      end if
                      j=j+1
                   end do
                   !l is next edge then the one last added
                   inn=.false.

                   do j=1,3
                      if (inner3(1,j).eqv..true.) then
                         if (inter%NumPnt(i)==1) then
                            if (l==j) then
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(1,j,1:2)
                               inner3(1,j)=.false.
                               inn=.false.
                            else
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i)+1,1:2)=xi3(1,j,1:2)
                               inner3(1,j)=.false.
                               inn=.true.
                               ed=j
                            end if
                         else
                            if (inn.eqv..false.) then
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(1,j,1:2)
                               inner3(1,j)=.false.
                               ed=j
                            else
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i)-1,1:2)=xi3(1,j,1:2)
                               inner3(1,j)=.false.
                            end if
                         end if
                      end if
                   end do

                   do j=1,3
                      if (inner3(2,j).eqv..true.) then
                         if (inter%NumPnt(i)==3) then
                            if (ed==j) then
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(2,j,1:2)
                               inner3(2,j)=.false.
                               inn=.false.
                            else
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i)+1,1:2)=xi3(2,j,1:2)
                               inner3(2,j)=.false.
                               inn=.true.
                            end if
                         else
                            if (inn.eqv..false.) then
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(2,j,1:2)
                               inner3(2,j)=.false.
                            else
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%IntPnt(i,inter%NumPnt(i)-1,1:2)=xi3(2,j,1:2)
                               inner3(2,j)=.false.
                            end if
                         end if
                      end if
                   end do
                end if
             end if

             !cnt+cnt2+cnt3
          case(4)
             !!print*,"cnt+cnt2+cnt3",cnt+cnt2+cnt3
             !print*,"cnt",cnt,"cnt2",cnt2,"cnt3",cnt3
             inter%had(i)=.true.
             !are the intersect doubled on the edge or not?
             edg1(1:3)=0
             do l=1,3
                do j=1,3
                   if (inner3(l,j) .eqv. .true.) then
                      edg1(j)=edg1(j)+1
                   end if
                end do
             end do
             !print*,"edg1",edg1(1:3)
             !((26/4/4.4 ))
             !count of interpoints on edges are 2,2,1

            if ((cnt/=0) .and. (cnt2 /=0) .and. (cnt3/=0) ) then
                inter%NumPnt(i)=0
                if (cnt==2) then
                   !adding the ones which are one 1&2
                   !edg(1),edgi(2) are numbers of edges of K
                   do j=1,3
                      if (inner3(1,j) .eqv. .true.) then
                         inter%NumPnt(i)=inter%NumPnt(i)+1
                         inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(1,j,1:2)
                         !print*,"adding xi3 1",j
                         edgi(1)=j
                         !print*,"edgi(1)",edgi(1)
                         inner3(1,j)=.false.
                      end if
                   end do
                   !!print*,"edgi(1)",edgi(1)
                   !adding intersects with 2 and with edge edgi(1)
                   j=2
                   do
                      if (inner3(j,edgi(1)) .eqv. .true.) then
                         inter%NumPnt(i)=inter%NumPnt(i)+1
                         inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(j,edgi(1),1:2)
                         !print*,"adding xi3",j,edgi(1)
                         inner3(j,edgi(1))=.false.
                         exit
                      end if
                      if (j==3) then
                         !print*,"calling GetVertex",edgi(1)
                         call GetVertex(edgi(1),vr1,vr2)
                         call Inside2(nodesK(vr1,1:2),nodesi,inn)
                         call Inside2(nodesK(vr2,1:2),nodesi,inn2)
                         !print*,"Inside2 on",vr1,inn
                         !print*,"Inside2 on",vr2,inn2
                         if (inn .eqv. .true.) then
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(vr1,1:2)
                         end if
                         if (inn2 .eqv. .true.) then
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(vr2,1:2)
                         end if
                         if ((inn .eqv. .true.).and. (inn2 .eqv. .true.))then
                            return
                         end if
                         exit
                      end if
                      j=j+1
                   end do
                   do l=1,3
                      do j=1,3
                         if (inner3(l,j) .eqv. .true.) then
                            call Differ(inter%intPnt(i,inter%NumPnt(i),1:2),xi3(l,j,1:2),dif(1))
                            if (dif(1) .eqv. .false.) then
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(l,j,1:2)
                               return
                            end if
                         end if
                      end do
                   end do
                end if

                if (cnt2==2) then
                   !adding the ones which are on 2&3
                   !edg(1),edgi(2) are numbers of edges of K
                   do j=1,3
                      if (inner3(2,j) .eqv. .true.) then
                         inter%NumPnt(i)=inter%NumPnt(i)+1
                         inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(2,j,1:2)
                         !print*,"adding xi3 2",j
                         edgi(1)=j
                         !print*,"edgi(1)",edgi(1)
                         inner3(2,j)=.false.
                      end if
                   end do
                   !print*,"edgi(1)",edgi(1)
                   !adding else
                   !print*,"inner3",inner3(1,1:3),inner3(2,1:3),inner3(3,1:3)
                   j=1
                   do
                      if (inner3(j,edgi(1)) .eqv. .true.) then
                         inter%NumPnt(i)=inter%NumPnt(i)+1
                         inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(j,edgi(1),1:2)
                         !print*,"adding xi3",j,edgi(1)
                         inner3(j,edgi(1))=.false.
                         exit
                      end if
                      if (j==3) then
                         call GetVertex(edgi(1),vr1,vr2)
                         call Inside2(nodesK(vr1,1:2),nodesi,inn)
                         call Inside2(nodesK(vr2,1:2),nodesi,inn2)
                         if (inn .eqv. .true.) then
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(vr1,1:2)
                         end if
                         if (inn2 .eqv. .true.) then
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(vr2,1:2)
                         end if
                         if ((inn .eqv. .true.).and. (inn2 .eqv. .true.))then
                            return
                         end if

                         exit
                      end if
                      j=j+2
                   end do

                   do l=1,3
                      do j=1,3
                         if (inner3(l,j) .eqv. .true.) then
                            call Differ(inter%intPnt(i,inter%NumPnt(i),1:2),xi3(l,j,1:2),dif(1))
                            if (dif(1) .eqv. .false.) then
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(l,j,1:2)
                            end if
                            return
                         end if
                      end do
                   end do
                end if
                !cnt3=2
                if (cnt3==2) then
                   !adding the ones which are on 1&3
                   !edg(1),edgi(2) are numbers of edges of K
                   do j=1,3
                      if (inner3(3,j) .eqv. .true.) then
                         inter%NumPnt(i)=inter%NumPnt(i)+1
                         inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(3,j,1:2)
                         !print*,"adding xi3 3",j
                         edgi(1)=j
                         !print*,"edgi(1)",edgi(1)
                         inner3(3,j)=.false.
                      end if
                   end do
                   !print*,"edgi(1)",edgi(1)
                   !adding else
                   j=1
                   do
                      if (inner3(j,edgi(1)) .eqv. .true.) then

                          inter%NumPnt(i)=inter%NumPnt(i)+1
                          inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(j,edgi(1),1:2)

                        !print*,"adding xi3",j,edgi(1)
                         inner3(j,edgi(1))=.false.
                         exit
                      end if
                      if (j==2) then
                         call GetVertex(edgi(1),vr1,vr2)
                         call Inside2(nodesK(vr1,1:2),nodesi,inn)
                         call Inside2(nodesK(vr2,1:2),nodesi,inn2)
                         if (inn .eqv. .true.) then
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(vr1,1:2)
                         end if
                         if (inn2 .eqv. .true.) then
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(vr2,1:2)
                         end if
                         if ((inn .eqv. .true.).and. (inn2 .eqv. .true.))then
                            return
                         end if
                         exit
                      end if
                      j=j+1
                   end do

                   do l=1,3
                      do j=1,3
                         if (inner3(l,j) .eqv. .true.) then
                            call Differ(inter%intPnt(i,inter%NumPnt(i),1:2),xi3(l,j,1:2),dif(1))
                            if (dif(1) .eqv. .false.) then
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(l,j,1:2)
                               return
                            end if
                         end if
                      end do
                   end do
                end if
                !end if

             else
                !((26/4/4.5))
                !deciding if there's inside nodesK
                inn2=.false.
                do j=1,3
                   inn=.false.
                   call Inside(nodesK(j,1:2),nodesi,inn)
                   !print*,"inn",inn
                   if (inn .eqv. .true.) then
                      !print*,"Inside nodesK",vr
                      vr=j
                      inn2=.true.
                   end if
                end do
                if (inn2 .eqv..false.) then
                   m=1
                   j=1
                   do
                      if (m==4) then
                         exit
                      end if

                     do
                         if (inner3(m,j) .eqv. .true.) then
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(m,j,1:2)
                            inner3(m,j)=.false.
                            !2 consecutive
                            if (j==3) then
                               ed=1
                            else
                               ed=j+1
                            end if
                            if (inner3(m,ed) .eqv. .true.) then
                               inter%NumPnt(i)=inter%NumPnt(i)+1
                               inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(m,ed,1:2)
                               inner3(m,ed)=.false.
                               !print*,"ed",ed
                               inn=.false.
                               !rest
                               do n=1,3
                                  do l=1,3
                                     if ((inner3(n,l) .eqv. .true.)) then
                                        !print*,"n,l",n,l
                                        if (inn.eqv..false.) then
                                           !print*,"n,l",n,l
                                           if ((ed==l).or.(inter%NumPnt(i)>=3) )then
                                              inter%NumPnt(i)=inter%NumPnt(i)+1
                                              inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(n,l,1:2)
                                              inner3(n,l)=.false.
                                           else
                                              inter%NumPnt(i)=inter%NumPnt(i)+1
                                              inter%intPnt(i,inter%NumPnt(i)+1,1:2)=xi3(n,l,1:2)
                                              inner3(n,l)=.false.
                                              inn=.true.
                                           end if
                                        else
                                           inter%NumPnt(i)=inter%NumPnt(i)+1
                                           inter%intPnt(i,inter%NumPnt(i)-1,1:2)=xi3(n,l,1:2)
                                           inner3(n,l)=.false.
                                           return
                                        end if
                                     end if

                                  end do
                               end do
                               if (inter%NumPnt(i)>=4) then
                                  return
                               end if
                               !not consecutive, add 2,3 and then find the fourth one
                            else
                               inn2=.false.
                               do l=1,3
                                  if (inner3(m+1,l) .eqv. .true.) then
                                     inter%NumPnt(i)=inter%NumPnt(i)+1
                                     inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(m+1,l,1:2)
                                     inner3(m+1,l)=.false.
                                     !print*,"adding m+1",m+1,l
                                     inn2=.true.
                                  end if
                               end do
                               if (inn2.eqv..false.) then
                                  !print*,"m",m
                                  if (m==1) then
                                     m=3
                                  else
                                     if (m==2) then
                                        m=1
                                     else
                                        m=2
                                     end if
                                  end if
                                  do l=1,3
                                     if (inner3(m,l) .eqv. .true.) then
                                        inter%NumPnt(i)=inter%NumPnt(i)+1
                                        inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(m,l,1:2)
                                        inner3(m,l)=.false.
                                          !print*,"adding m+2",m,l
                                     end if
                                  end do
                               end if

                               do n=1,3
                                  do l=1,3
                                     if (inner3(n,l) .eqv. .true. )then
                                        inter%NumPnt(i)=inter%NumPnt(i)+1
                                        inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(n,l,1:2)
                                        inner3(n,l)=.false.
                                        !print*,"adding",n,l
                                     end if
                                  end do
                               end do
                               if (inter%NumPnt(i)>=4)then
                                  return
                               end if
                            end if
                         end if
                         if (j==3) then
                            j=1
                            m=m+1
                            exit
                         end if
                         j=j+1


                      end do
                   end do
                   !some vertex is inside
                else
                   !add the vertex inside
                   inter%NumPnt(i)=1
                   inter%IntPnt(i,1,1:2)=nodesK(vr,1:2)
                end if
                !add the intersection with the edge of elem. K number ed
                j=1
                do
                   if(inner3(j,vr).eqv..true.)then
                      inter%NumPnt(i)=2
                      inter%IntPnt(i,2,1:2)=xi3(j,vr,1:2)
                      inner3(j,vr)=.false.
                      exit
                   end if
                   j=j+1
                end do
                !add the intersection with edge of elem. i number j
                l=1
                do
                   if (inner3(j,l).eqv..true.) then
                      inter%NumPnt(i)=3
                      inter%IntPnt(i,3,1:2)=xi3(j,l,1:2)
                      inner3(j,l)=.false.
                      exit
                   end if
                   l=l+1
                end do
                !intersection with edge of elem K number 1
                j=1
                do
                   if (inner3(j,l).eqv..true.)then
                      inter%NumPnt(i)=4
                      inter%IntPnt(i,4,1:2)=xi3(j,l,1:2)
                      inner3(j,l)=.false.
                      exit
                   end if
                   j=j+1
                end do
                !add the intersection with edge of elem. i number j
                l=1
                do
                   if (inner3(j,l).eqv..true.)then
                      inter%NumPnt(i)=5
                      inter%IntPnt(i,5,1:2)=xi3(j,l,1:2)
                      inner3(j,l)=.false.
                      exit
                   end if
                   l=l+1
                end do


             end if


             !cnt+cnt2+cnt3
          case(3)
             !print*,"cnt+cnt2+cnt3",cnt+cnt2+cnt3
             inter%had(i)=.true.
             inter%NumPnt(i)=0
             call Inside2(nodesK(1,1:2),nodesi,inn)
             !print*,"inside for nodesk 1",inn
             if (inn .eqv. .true.) then
                inter%NumPnt(i)=inter%NumPnt(i)+1
                inter%intPnt(i,inter%NumPnt(i),1:2)=nodesK(1,1:2)
             end if
             !do l=1,3
             !adding intersections points
             m=1
             do j=1,3
                if (inner3(j,1) .eqv. .true.) then
                   dif=.false.
                   do l=1,inter%NumPnt(i)
                      call Differ(xi3(j,1,1:2),inter%IntPnt(i,l,1:2),dif(l))
                   end do
                   call Differ(xi3(j,1,1:2),nodesK(2,1:2),line(1))
                   call Differ(xi3(j,1,1:2),nodesK(3,1:2),line(2))
                   if ((dif(1).eqv..false.).and.(dif(2).eqv..false.).and.(dif(3).eqv..false.)&
                        .and.(dif(4).eqv..false.).and.(line(1).eqv..false.).and.(line(2).eqv..false.)) then
                      !inter%NumPnt(i)=inter%NumPnt(i)+1
                      !inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(j,1,1:2)
                      edgi(m)=j
                      m=m+1
                      !print*,"adding",j,"1 intersection"
                   end if
                end if
             end do
             !print*,"m",m
             if (m/=1) then
                if (m==3) then
                   !print*,"edgi",edgi(1:2)
                   call FillInOrder(edgi(1:2),inter%NumPnt(i),inter%IntPnt(i,1:6,1:2),xi3(edgi(1),1,1:2),xi3(edgi(2),1,1:2))
                else
                   inter%NumPnt(i)=inter%NumPnt(i)+1
                   inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(edgi(1),1,1:2)
                end if
             end if
             call Inside2(nodesK(2,1:2),nodesi,inn)
             !print*,"inside for nodesk 2",inn
             if (inn .eqv. .true.) then
                dif=.false.
                do l=1,inter%NumPnt(i)
                   call Differ(nodesK(2,1:2),inter%IntPnt(i,l,1:2),dif(l))
                end do

                if ((dif(1).eqv..false.).and.(dif(2).eqv..false.).and.(dif(3).eqv..false.)&
                     .and.(dif(4).eqv..false.)) then
                   inter%NumPnt(i)=inter%NumPnt(i)+1
                   inter%intPnt(i,inter%NumPnt(i),1:2)=nodesK(2,1:2)
                   !print*,"adding nodesK 2"
                end if
             end if
             m=1
             do j=1,3
                if (inner3(j,2) .eqv. .true.) then
                   dif=.false.
                   do l=1,inter%NumPnt(i)
                      call Differ(xi3(j,2,1:2),inter%IntPnt(i,l,1:2),dif(l))
                   end do
                   call Differ(xi3(j,1,1:2),nodesK(3,1:2),line(1))
                   if ((dif(1).eqv..false.).and.(dif(2).eqv..false.).and.(dif(3).eqv..false.)&
                        .and.(dif(4).eqv..false.).and.(line(1).eqv..false.)) then
                      edgi(m)=j
                      m=m+1
                      !print*,"adding",j,"2 intersection"
                   end if
                end if
             end do
             !print*,"m",m
             if (m/=1) then
                if (m==3) then
                   !print*,"edgi",edgi(1:2)
                   call FillInOrder(edgi(1:2),inter%NumPnt(i),inter%IntPnt(i,1:6,1:2),xi3(edgi(1),2,1:2),xi3(edgi(2),2,1:2))
                else
                   inter%NumPnt(i)=inter%NumPnt(i)+1
                   inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(edgi(1),2,1:2)
                end if
             end if

             call Inside2(nodesK(3,1:2),nodesi,inn)
             !print*,"inside for nodesk 3",inn
             if (inn .eqv. .true.) then
                dif=.false.
                do l=1,inter%NumPnt(i)
                   call Differ(nodesK(3,1:2),inter%IntPnt(i,l,1:2),dif(l))
                end do

                if ((dif(1).eqv..false.).and.(dif(2).eqv..false.).and.(dif(3).eqv..false.)&
                     .and.(dif(4).eqv..false.)) then
                   inter%NumPnt(i)=inter%NumPnt(i)+1
                   inter%intPnt(i,inter%NumPnt(i),1:2)=nodesK(3,1:2)
                   !print*,"adding nodesK 3"
                end if
             end if

             m=1
             do j=1,3
                if (inner3(j,3) .eqv. .true.) then
                   dif=.false.
                   do l=1,inter%NumPnt(i)
                      call Differ(xi3(j,3,1:2),inter%IntPnt(i,l,1:2),dif(l))
                   end do

                   if ((dif(1).eqv..false.).and.(dif(2).eqv..false.).and.(dif(3).eqv..false.)&
                        .and.(dif(4).eqv..false.)) then
                      edgi(m)=j
                      m=m+1
                      !print*,"adding",j,"3 intersection"
                   end if
                end if
             end do
             !print*,"m",m
             if (m/=1) then
                if (m==3) then
                   !print*,"edgi",edgi(1:2)
                   call FillInOrder(edgi(1:2),inter%NumPnt(i),inter%IntPnt(i,1:6,1:2),xi3(edgi(1),3,1:2),xi3(edgi(2),3,1:2))
                else
                   inter%NumPnt(i)=inter%NumPnt(i)+1
                   inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(edgi(1),3,1:2)
                end if
             end if
             return
             !cnt+cnt2+cnt3
          case(2)
             inter%had(i)=.true.
             !print*,"cnt+cnt2+cnt3",cnt+cnt2+cnt3
             !print*,"cnt",cnt
             !print*,"cnt2",cnt2
             !print*,"cnt3",cnt3
             !which one is the one with no intersection point
             do l=1,3
                do j=1,3
                   if (inner3(l,j) .eqv. .true.) then
                      inter%NumPnt(i)=inter%NumPnt(i)+1
                      inter%IntPnt(i,inter%NumPnt(i),1:2)=xi3(l,j,1:2)
                      !if ((K==287).and.(i==343))then
                      !write(19,*),"adding point",xi3(l,j,1:2)
                      !end if
                      !edge of K
                      edge=j
                   end if
                end do
             end do
             ! getting know how many of nodesK are inside of i element
             lin=0
             do l=1,3
                call Inside2(nodesK(l,1:2),nodesi,line(l))
                !if ((K==287).and.(i==343))then
                !	write(19,*),"line",line(l)
                !	end if
                if (line(l) .eqv. .true.) then

                   lin=lin+1

                end if
             end do
             !print*,"lin",lin
             select case(lin)
             case(3)
                inter%NumPnt(i)=3
                inter%IntPnt(i,1:3,1:2)=nodesK(1:3,1:2)
                fin=.true.
                return
             case(2)
                call GetVertex(edge,vr1,vr2)
                do l=1,3
                   dif=.false.
                   ! if l-th of nodesK is inside and equal to vr1 or vr2
                   ! and hasn't been already added then add him
                   if ((line(l) .eqv. .true.) .and. ((l==vr1) .or. (l==vr2)) ) then
                      do j=1,2
                         call Differ(inter%IntPnt(i,j,1:2),nodesK(l,1:2),dif(j))
                      end do
                      if ((dif(1) .eqv. .false.).and.(dif(2) .eqv. .false.)) then
                         inter%NumPnt(i)=inter%NumPnt(i)+1
                         inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(l,1:2)
                         line(l)=.false.
                      end if
                   end if
                end do
                ! if there aren't four added points then find the 4th inside
                if (inter%NumPnt(i)<4) then
                   do l=1,3
                      dif=.false.
                      if (line(l) .eqv. .true.) then
                         do j=1,3
                            call Differ(inter%IntPnt(i,j,1:2),nodesK(l,1:2),dif(j))
                         end do
                         if ((dif(1) .eqv. .false.).and.(dif(2) .eqv. .false.)) then
                            inter%NumPnt(i)=inter%NumPnt(i)+1
                            inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(l,1:2)
                            line(l)=.false.
                         end if
                      end if
                   end do
                end if
                return
             case(1)
                inter%NumPnt(i)=inter%NumPnt(i)+1
                do l=1,3
                   if (line(l) .eqv. .true.) then
                      inter%IntPnt(i,inter%NumPnt(i),1:2)=nodesK(l,1:2)
                      !write(19,*),"addong nodesK",l,nodesK(l,1:2)
                      return
                   end if
                end do
             case default
                !print*,"mistake in 3outside,2 intersections"
             end select

             !cnt+cnt2+cnt3
          case(1)
             !print*,"cnt+cnt2+cnt3",cnt+cnt2+cnt3
             do l=1,3
                call Inside2(nodesK(l,1:2),nodesi,line(l))
                !print*,"line(l)",line(l)
             end do
             if ( (line(1) .eqv. .true.) .and. (line(2) .eqv. .true.) .and. (line(3) .eqv. .true.) ) then
                inter%had(i)=.true.
                inter%NumPnt(i)=3
                !print*,"adding nodesK(1:3)"
                inter%IntPnt(i,1:3,1:2)=nodesK(1:3,1:2)
                fin=.true.
             else
                do l=1,3
                   do j=1,3
                      if (inner3(l,j) .eqv. .true.) then
                         inter%had(i)=.true.
                         inter%NumPnt(i)=inter%NumPnt(i)+1
                         inter%intPnt(i,inter%NumPnt(i),1:2)=xi3(l,j,1:2)
                         !print*,"adding the only one adding point"
                      end if
                   end do
                end do
             end if
             !cnt+cnt2+cnt3
          case(0)
             !print*,"cnt+cnt2+cnt3",cnt+cnt2+cnt3
             do l=1,3
                call Inside2(nodesK(l,1:2),nodesi,line(l))
               !if ((k==287).and.(i==343)) then
               !write(19,*),"line(l)",line(l)
               !end if

             end do
             !close(19)
             if ( (line(1) .eqv. .true.) .and. (line(2) .eqv. .true.) .and. (line(3) .eqv. .true.) ) then
                inter%had(i)=.true.
                inter%NumPnt(i)=3
                !print*,"inter%NumPnt(i)=3, second time"
                inter%IntPnt(i,1:3,1:2)=nodesK(1:3,1:2)
                fin=.true.
             else
                !write(19,*),"false"
                inter%had(i)=.false.
             end if
             return
             !cnt+cnt2+cnt3
          end select

          ! case of edg%num + ver%num
       case default
          print*, "Error - edg%num out of range"
       end select
    case default
       print*, "ERROR in IntersectTris - number of inside vertices is out of range"
    end select


  end subroutine IntersectTris


  !compare two triangles if they are same
  subroutine Same(tri1,tri2,sam)
    real,dimension(1:3,1:2),intent(in)	:: tri1,tri2
    logical,intent(out)			:: sam
    real, parameter			:: TOL=0.000005
    integer				:: l,k,j
    logical,dimension(1:9)		:: dif

    sam=.false.
    j=1
    do l=1,3
       do k=1,3
          call Differ(tri1(l,1:2),tri2(k,1:2),dif(j))
          j=j+1
       end do
    end do
    if(( (dif(1).eqv. .true.) .or. (dif(2).eqv. .true.) .or. (dif(3).eqv. .true.) )  .and. &
         ( (dif(4).eqv. .true.) .or. (dif(5).eqv. .true.) .or. (dif(6).eqv. .true.) )  .and. &
         ( (dif(7).eqv. .true.) .or. (dif(8).eqv. .true.) .or. (dif(9).eqv. .true.) )  ) then
       sam=.true.
    else
       sam=.false.
    end if
  end subroutine Same


  !subroutine to find intersection of line and triangle in case there's just one intersection point
  subroutine IntersectLineTri1(xa,xb,nodesK,inner,xi)
    real, dimension(1:2), intent(in)			:: xa,xb
    real, dimension(1:3,1:2), intent(in)			:: nodesK
    logical, intent(out)					:: inner
    real, dimension(1:2), intent(out)			:: xi


    inner = .false.

    !xa, xb intersects nodesK 1 & 2
    call SeekIntersect(xa,xb,nodesK(1,1:2),nodesK(2,1:2),xi,inner)
    if (inner .eqv. .true.) then
       return
    end if
    !xa, xb intersects nodesK 2 & 3 unless we have already found intersect with 1 & 2
    call SeekIntersect(xa,xb,nodesK(2,1:2),nodesK(3,1:2),xi,inner)
    if (inner .eqv. .true.) then
       return
    end if

    !xa,xb intersects nodesK 1 & 3
    call SeekIntersect(xa,xb,nodesK(1,1:2),nodesK(3,1:2),xi,inner)
    if (inner .eqv. .true.) then
       return
    end if

  end subroutine IntersectLineTri1




  !better version, which also passes number of edge with witch the intersect is
  subroutine IntersectLineTri1b(xa,xb,nodesK,inner,xi,edge)
    real, dimension(1:2), intent(in)			:: xa,xb
    real, dimension(1:3,1:2), intent(in)			:: nodesK
    logical, intent(out)					:: inner
    real, dimension(1:2), intent(out)			:: xi
    integer, intent(out)					:: edge


    inner = .false.
    edge = 0


    !xa, xb intersects nodesK 1 & 2
    call SeekIntersect(xa,xb,nodesK(1,1:2),nodesK(2,1:2),xi,inner)
    if (inner .eqv. .true.) then
       edge = 1
       return
    end if
    !xa, xb intersects nodesK 2 & 3 unless we have already found intersect with 1 & 2
    call SeekIntersect(xa,xb,nodesK(2,1:2),nodesK(3,1:2),xi,inner)
    if (inner .eqv. .true.) then
       edge = 2
       return
    end if

    !xa,xb intersects nodesK 1 & 3
    call SeekIntersect(xa,xb,nodesK(1,1:2),nodesK(3,1:2),xi,inner)
    if (inner .eqv. .true.) then
       edge = 3
       return
    end if

  end subroutine IntersectLineTri1b

  !subroutine to find intersection of line and triangle in case there's more than one intersection point
  !cnt is number of intersection points with number of edge on which the intersect is
  subroutine IntersectLineTri2b(xa,xb,nodesK,inner,xi,cnt,edge)
    real, dimension(1:2), intent(in)			:: xa,xb
    real, dimension(1:3,1:2), intent(in)			:: nodesK
    logical,dimension(1:3), intent(out)			:: inner
    real, dimension(1:3,1:2), intent(out)			:: xi
    integer, intent(inout)				:: cnt
    integer,dimension(1:3), intent(out)			:: edge
    logical						:: dif,dif2

    inner(1:3) = .false.
    edge=0
    cnt=0

    !xa, xb intersects nodesK 1 & 2
    call SeekIntersect(xa,xb,nodesK(1,1:2),nodesK(2,1:2),xi(1,1:2),inner(1))
    if (inner(1) .eqv. .true.) then
       cnt = cnt+1
       edge(1)=1
    end if
    !xa, xb intersects nodesK 2 & 3
    call SeekIntersect(xa,xb,nodesK(2,1:2),nodesK(3,1:2),xi(2,1:2),inner(2))
    if (inner(2) .eqv. .true.) then
       call Differ(xi(1,1:2),xi(2,1:2),dif)
       if (dif .eqv. .false.) then
          cnt = cnt+1
          edge(2)=2
       end if
    end if

    !xa,xb intersects nodesK 1 & 3
    call SeekIntersect(xa,xb,nodesK(1,1:2),nodesK(3,1:2),xi(3,1:2),inner(3))
    if (inner(3) .eqv. .true.) then
       call Differ(xi(1,1:2),xi(3,1:2),dif)
       call Differ(xi(2,1:2),xi(3,1:2),dif2)
       if ((dif .eqv. .false.) .and. (dif2 .eqv. .false.)) then
          cnt = cnt+1
          edge(3)=3
       end if
    end if

  end subroutine IntersectLineTri2b

  !subroutine to find intersection of line and triangle in case there's more than one intersection point
  !cnt is number of intersection points
  subroutine IntersectLineTri2(xa,xb,nodesK,inner,xi,cnt)
    real, dimension(1:2), intent(in)			:: xa,xb
    real, dimension(1:3,1:2), intent(in)			:: nodesK
    logical,dimension(1:3), intent(out)			:: inner
    real, dimension(1:3,1:2), intent(out)			:: xi
    integer, intent(inout)				:: cnt
    logical						:: dif,dif2


    inner(1:3) = .false.
    cnt=0

    !xa, xb intersects nodesK 1 & 2
    call SeekIntersect(xa,xb,nodesK(1,1:2),nodesK(2,1:2),xi(1,1:2),inner(1))
    !!print*,"inner(1)",inner(1)
    if (inner(1) .eqv. .true.) then
       !!print*,"intersect with nodesK 1&2)"
       cnt = cnt+1
    end if

    !xa,xb intersects nodesK 2 & 3
    call SeekIntersect(xa,xb,nodesK(2,1:2),nodesK(3,1:2),xi(2,1:2),inner(2))
    !!print*,"inner(2)",inner(2)
    if (inner(2) .eqv. .true.) then
       ! !print*,"intersect with nodesK 2&3)"
       if (inner(1) .eqv. .true.) then
          call Differ(xi(1,1:2),xi(2,1:2),dif)
          ! !print*,"dif",dif
          if (dif .eqv. .false.) then
             cnt = cnt+1
          else
             inner(2)=.false.
          end if
       else
          cnt=cnt+1
       end if
    end if

    !xa, xb intersects nodesK 1 & 3
    call SeekIntersect(xa,xb,nodesK(1,1:2),nodesK(3,1:2),xi(3,1:2),inner(3))
    !!print*,"inner(3)",inner(3)
    if (inner(3) .eqv. .true.) then
       dif=.false.
       dif2=.false.
       !print*,"intersect with nodesK 1&3)"
       if (inner(1) .eqv. .true.) then
          call Differ(xi(1,1:2),xi(3,1:2),dif)
          !!print*,"dif 1&3",dif
       end if
       if (inner(2) .eqv. .true.) then
          call Differ(xi(2,1:2),xi(3,1:2),dif2)
          !!print*,"dif2 1&3",dif2
       end if
       if ((dif .eqv. .false.) .and. (dif2 .eqv. .false.)) then
          cnt = cnt+1
       else
          inner(3)=.false.
       end if
    end if

  end subroutine IntersectLineTri2

  !subroutine to find intersection of line and triangle in case there's just one intersection point
  !without the points which are already added into inter
  subroutine IntersectLineTri3(xa,xb,nodesK,inner,xi,NumPnt,IntPnt)
    real, dimension(1:2), intent(in)			:: xa,xb
    real, dimension(1:3,1:2), intent(in)			:: nodesK
    logical, intent(out)					:: inner
    real, dimension(1:2), intent(out)			:: xi
    integer,intent(in)					:: NumPnt
    real, dimension(1:6,1:2),intent(in)			:: IntPnt
    logical, dimension(1:6)				:: difer
    logical						:: dif
    integer						:: l

    inner = .false.

    !xa, xb intersects nodesK 1 & 2
    call SeekIntersect(xa,xb,nodesK(1,1:2),nodesK(2,1:2),xi,inner)
    if (inner .eqv. .true.) then
       !write(31,*),"intersection found xi,1&2",xi
       difer=.false.
       if (NumPnt>0) then
          do l=1,NumPnt
             call Differ(xi,IntPnt(l,1:2),difer(l))
         end do
         !write(31,*),"differ",difer(1:NumPnt)
         if ((difer(1) .eqv. .false.) .and.(difer(2) .eqv. .false.) .and. (difer(3) .eqv. .false.).and. &
             (difer(4) .eqv. .false.) .and. (difer(5) .eqv. .false.) .and. (difer(6) .eqv. .false.))then
              return
          else
             inner=.false.
             !write(31,*), "fallse"
          end if
       end if

    end if
    !xa, xb intersects nodesK 2 & 3 unless we have already found intersect with 1 & 2
    call SeekIntersect(xa,xb,nodesK(2,1:2),nodesK(3,1:2),xi,inner)
    if (inner .eqv. .true.) then
       !write(31,*),"intersection found xi,2&3",xi
       difer=.false.
       if (NumPnt>0) then
          do l=1,NumPnt
             call Differ(xi,IntPnt(l,1:2),difer(l))
          end do
          !write(31,*),"differ",difer(1:NumPnt)
          if ((difer(1) .eqv. .false.) .and.(difer(2) .eqv. .false.) .and. (difer(3) .eqv. .false.).and. &
               (difer(4) .eqv. .false.) .and. (difer(5) .eqv. .false.) .and. (difer(6) .eqv. .false.))then
             return
          else
             inner=.false.
             !write(31,*), "fallse"
          end if

       end if

    end if
    !xa,xb intersects nodesK 1 & 3
    call SeekIntersect(xa,xb,nodesK(1,1:2),nodesK(3,1:2),xi,inner)
    if (inner .eqv. .true.) then
       !write(31,*),"intersection found xi,1&3",xi
       difer=.false.
       if (NumPnt>0) then
          do l=1,NumPnt
             call Differ(xi,IntPnt(l,1:2),difer(l))
          end do
          !write(31,*),"differ",difer(1:NumPnt)
          if ((difer(1) .eqv. .false.) .and.(difer(2) .eqv. .false.) .and. (difer(3) .eqv. .false.).and. &
               (difer(4) .eqv. .false.) .and. (difer(5) .eqv. .false.) .and. (difer(6) .eqv. .false.))then
             return
          else
             !write(31,*),"false"
             inner=.false.
          end if

       end if

    end if
    !print*, "No intersection in Intersectline3 found"


  end subroutine IntersectLineTri3



  !subroutine to find intersection of line and triangle in case there's just one intersection point
  !without the points which are already added into inter
  !plus edge is number of edge with wich thwe intersect is
  subroutine IntersectLineTri3b(xa,xb,nodesK,inner,xi,NumPnt,IntPnt,edge)
    real, dimension(1:2), intent(in)			:: xa,xb
    real, dimension(1:3,1:2), intent(in)			:: nodesK
    logical, intent(out)					:: inner
    real, dimension(1:2), intent(out)			:: xi
    integer,intent(in)					:: NumPnt
    real, dimension(1:6,1:2),intent(in)			:: IntPnt
    integer, intent(out)					:: edge
    logical, dimension(1:6)				:: difer
    logical						:: dif
    integer						:: l

    inner = .false.

    !xa, xb intersects nodesK 1 & 2
    call SeekIntersect(xa,xb,nodesK(1,1:2),nodesK(2,1:2),xi,inner)
    if (inner .eqv. .true.) then
       !!print*,"intersection found xi,1&2",xi
       dif=.false.
       do l=1,NumPnt
          call Differ(xi,IntPnt(l,1:2),difer(l))
          if (difer(l) .eqv. .true.) then
             dif=.true.
             inner = .false.
          end if
       end do
       if (dif .eqv. .false.) then
          inner=.true.
          !!print*,"IntersectLineTri3-intersetion with nodesK 1&2"
          edge=1
          return
       end if
    end if
    !xa, xb intersects nodesK 2 & 3 unless we have already found intersect with 1 & 2
    call SeekIntersect(xa,xb,nodesK(2,1:2),nodesK(3,1:2),xi,inner)
    if (inner .eqv. .true.) then
       !!print*,"intersection found xi 1&3",xi
       dif=.false.
       do l=1,NumPnt
          call Differ(xi,IntPnt(l,1:2),difer(l))
          if (difer(l) .eqv. .true.) then
             dif=.true.
             inner=.false.
          end if
       end do
       if (dif .eqv. .false.) then
          inner=.true.
          !!print*,"IntersectLineTri3-intersetion with nodesK 2&3"
          edge =2
          return
       end if
    end if

    !xa,xb intersects nodesK 1 & 3
    call SeekIntersect(xa,xb,nodesK(1,1:2),nodesK(3,1:2),xi,inner)
    if (inner .eqv. .true.) then
       dif=.false.
       !!print*,"intersection found xi,1&3",xi
       do l=1,NumPnt
          call Differ(xi,IntPnt(l,1:2),difer(l))
          if (difer(l) .eqv. .true.) then
             dif=.true.
             inner=.false.
          end if
       end do
       if (dif .eqv. .false.) then
          inner=.true.
          !!print*,"IntersectLineTri3-intersetion with nodesK 2&3"
          edge = 3
          return
       end if
    end if
    edge=4
    !print*, "No intersection in Intersectline3 found"

  end subroutine IntersectLineTri3b


  !subroutine which says if the points are the same or not
  subroutine Differ(xa,xb,dif)
    real, dimension(1:2),intent(in)	:: xa,xb
    logical,intent(out)			:: dif
    real, parameter			:: TOL=0.000001

    if ( (abs(xa(1)-xb(1))<TOL) .and. (abs(xa(2)-xb(2))<TOL) ) then
       dif=.true.
    else
       dif=.false.
    end if
  end subroutine Differ

  !subroutine, which finds out, where is the node D in view of K (inside,outside,on edge, in vertex)
  subroutine PosNodes(nodesi,nodesK,post,ins,edg,ver,ous)
    type(pos),intent(inout)			:: ins, edg, ver, ous
    real, dimension(1:3,1:2),intent(in)		:: nodesK,nodesi
    character(len=2),dimension(1:3),intent(out)	:: post		!position of vertex D
    real, dimension(1:2)				:: t, xi	!parameter d, intersection point
    logical					:: inner
    real, parameter				:: TOL=0.000005
    integer					:: j


    !initiate ins,edg,ver,ots
    ins%num=0
    ins%is(1:3)=.false.
    edg%num=0
    edg%is(1:3)=.false.
    ver%num=0
    ver%is(1:3)=.false.
    ous%num=0
    ous%is(1:3)=.false.
    do j=1,3
       !initiate
       !!print*,j
       post(j) = "00"
       t(1) = -5
       t(2) = -5

       !if C is not the same as D then call seeking of intersection
       if ((abs(nodesK(3,1)-nodesi(j,1))<= TOL) .and. (abs(nodesK(3,2)-nodesi(j,2))<= TOL)) then
          ! D is vertex C
          post(j)="VC"
          !		!print*,"VC"
          ver%num = ver%num+1
          ver%is(j) = .true.
          cycle
       else
          call SeekIntersectLine(nodesK(1,1:2), nodesK(2,1:2), nodesK(3,1:2), nodesi(j,1:2), xi, t, inner)
          !print*,"j,t,inner",j,t, inner
       end if
       !according to values of parameter vD, we know, were the vertex D is in view of triangle K
       !D is inside K
       if ( (TOL < t(1)) .and. (t(1) < 1-TOL) .and. (t(2) > 1+TOL) ) then
          post(j)="IN"
          !	!print*,"IN"
          ins%num = ins%num+1
          ins%is(j) = .true.
          cycle
       end if
       !D is outside K
       if ( (t(1) > 1 + TOL ) .or. (t(1) < - TOL) .or. ((TOL < t(2)).and. (t(2) < 1-TOL)) .or.    &
            (t(2) < -TOL).or.  (inner .eqv. .false.) ) then
          post(j)="OU"
          !	!print*,"OU"
          ous%num = ous%num+1
          ous%is(j) = .true.
          cycle
       end if
       !D is in AB, especially D=B
       if (abs(t(2)-1) <= TOL) then
          if ( abs(t(1)-1) <= TOL ) then
             post(j)="VB"
             !		!print*,"VB"
             ver%num = ver%num+1
             ver%is(j) = .true.
             cycle
             return
          else
             if (abs(t(1)) <= TOL) then
                post(j)="VA"
                !				!print*,"VA"
                ver%num = ver%num+1
                ver%is(j) = .true.
                cycle
             else
                post(j)="AB"
                !				!print*,"AB"
                edg%num = edg%num+1
                edg%is(j) = .true.
                cycle
             end if
             return
          end if
       end if

       !D is in BC, intersection point is B
       !if t(1)=1 and t(2)/=1
       if ( (abs(t(1)-1) <= TOl) .and. (abs(t(2)-1) >= TOL) )then
          post(j)="BC"
          !	!print*,"BC"
          !	!print*,"edg%num",edg%num
          edg%num = edg%num+1
          edg%is(j) = .true.
          cycle
       end if
       !D is in CA, intersection point is A, especially D is A
       !if t(1)=0;t(2)=1
       if (abs(t(1)) <= TOL) then
          if (abs(t(2)-1) <= TOL) then
             post(j)="VA"
             !		!print*,"VA"
             ver%num = ver%num+1
             ver%is(j) = .true.
             cycle
          else
             post(j)="CA"
             !		!print*,"CA"
             edg%num = edg%num+1
             edg%is(j) = .true.
             cycle
          end if
       end if
    end do
  end subroutine PosNodes

  !> intersection of lines (xa, xb) and (xc, xd) is xi if inner = .true.
  !inter is false, when the lines are parallel
  subroutine SeekIntersectLine(xa, xb, xc, xd, xi,t, inter)
    real, dimension(1:2), intent(in) :: xa, xb, xc, xd
    real, dimension(1:2), intent(out) :: xi,t
    logical, intent(out) :: inter
    real, dimension (1:2, 1:2) :: A, A1
    real, dimension (1:2) :: b
    !real, dimension (1:2), intent(out) ::t
    real :: det

    inter = .false.

    A(1:2,1) =   xb(1:2) - xa(1:2)
    A(1:2,2) = -(xd(1:2) - xc(1:2) )

    b(1:2) = xc(1:2) - xa(1:2)

    det = A(1,1) * A(2,2) - A(1,2) * A(2,1)

    if(det == 0) return  ! paralell lines


    A1(1,1) =  A(2,2) / det
    A1(1,2) = -A(1,2) / det
    A1(2,1) = -A(2,1) / det
    A1(2,2) =  A(1,1) / det

    t(1:2) = matmul(A1(1:2,1:2), b(1:2) )
    xi(1:2) = xa(1:2) + t(1) * (xb(1:2)  - xa(1:2) )
    inter = .true.

  end subroutine SeekIntersectLine


  ! intersection of lines, return true if it's on the abcissa xa,xb and xc,xd
  subroutine SeekIntersect(xa, xb, xc, xd, xi, inter)
    real, dimension(1:2), intent(in) :: xa, xb, xc, xd
    real, dimension(1:2), intent(out) :: xi
    logical, intent(inout) :: inter
    real, dimension (1:2, 1:2) :: A, A1
    real, dimension (1:2) :: b, t
    real :: det
    real, parameter			:: TOL=0.000005

    inter = .false.

    A(1:2,1) =   xb(1:2) - xa(1:2)
    A(1:2,2) = -(xd(1:2) - xc(1:2) )

    b(1:2) = xc(1:2) - xa(1:2)

    det = A(1,1) * A(2,2) - A(1,2) * A(2,1)

    if(det == 0) return  ! paralell lines


    A1(1,1) =  A(2,2) / det
    A1(1,2) = -A(1,2) / det
    A1(2,1) = -A(2,1) / det
    A1(2,2) =  A(1,1) / det

    t(1:2) = matmul(A1(1:2,1:2), b(1:2) )
    !!print*,"t",t(1:2)

    if(t(1) >= -TOL .and. t(1) <= (1+TOL) .and. t(2) >= -TOL .and. t(2) <= (1+TOL) )then
       xi(1:2) = xa(1:2) + t(1) * (xb(1:2)  - xa(1:2) )
       inter = .true.
    end if
  end subroutine SeekIntersect


  !seek intersect without the edge points
  subroutine SeekIntersectOutOf(xa, xb, xc, xd, xi, inter)
    real, dimension(1:2), intent(in) :: xa, xb, xc, xd
    real, dimension(1:2), intent(out) :: xi
    logical, intent(inout) :: inter
    real, dimension (1:2, 1:2) :: A, A1
    real, dimension (1:2) :: b, t
    real :: det
    real, parameter	:: TOL=0.000005
    inter = .false.

    A(1:2,1) =   xb(1:2) - xa(1:2)
    A(1:2,2) = -(xd(1:2) - xc(1:2) )

    b(1:2) = xc(1:2) - xa(1:2)

    det = A(1,1) * A(2,2) - A(1,2) * A(2,1)

    if(det == 0) then
       !print*,"SeekIntersectOutOf det == 0"
       return  ! paralell lines
    end if

    A1(1,1) =  A(2,2) / det
    A1(1,2) = -A(1,2) / det
    A1(2,1) = -A(2,1) / det
    A1(2,2) =  A(1,1) / det

    t(1:2) = matmul(A1(1:2,1:2), b(1:2) )
    !    !print*,"t",t(1:2)
    if((t(1) > TOL ).and. (t(1) < 1-TOL) .and. (t(2) > TOL) .and. (t(2) < 1-TOL )) then
       xi(1:2) = xa(1:2) + t(1) * (xb(1:2)  - xa(1:2) )
       inter = .true.
    end if
  end subroutine SeekIntersectOutOf

  !Get vertices of edge
  subroutine GetVertex(edge,vr1,vr2)
    integer,intent(in)		:: edge
    integer, intent(out)		:: vr1,vr2

    select case(edge)
    case(1)
       vr1=1
       vr2=2
    case(2)
       vr1=2
       vr2=3
    case(3)
       vr1=3
       vr2=1
    end select
  end subroutine GetVertex


  ! from the numbers of edges returns number of the third edge and her vertices
  subroutine GetOtherEdge(num1,num2,num3,vr1,vr2)
    integer, intent(in)		:: num1, num2
    integer, intent(out)		:: num3
    integer, intent(out)		:: vr1,vr2

    select case(num1+num2)
    case(3)
       num3=3
       vr1=2
       vr2=3
    case(4)
       num3=2
       vr1=1
       vr2=3
    case(5)
       num3=1
       vr1=1
       vr2=2
    end select
  end subroutine GetOtherEdge

  !post of edge to number of edge
  subroutine GetNumber(post,num)
    character(len=2), intent(in)		:: post
    integer,intent(out)			:: num
    if (post=="AB") then
       num=1
    else
       if (post=="CA") then
          num=2
       else
          num=3
       end if
    end if
  end subroutine GetNumber

  !subroutine which returns ins true, when x is inside the triangle
  subroutine Inside(x,triang,ins)
    real, dimension(1:2),intent(in)			:: x
    real, dimension(1:3,1:2),intent(in)			:: triang
    logical, intent(out)					:: ins
    real, dimension(1:3)					:: lambda
    real							:: detT
    real, parameter					:: TOL=0.000005
    ins=.false.
    detT = ( ( (triang(1,1)-triang(3,1))*(triang(2,2)-triang(3,2)) ) - ( (triang(1,2)-triang(3,2))*(triang(2,1)-triang(3,1)) ) )
    lambda(1)=((triang(2,2)-triang(3,2))*(x(1)-triang(3,1)) + (triang(3,1)-triang(2,1))*(x(2)-triang(3,2)) )/detT
    lambda(2)=((triang(3,2)-triang(1,2))*(x(1)-triang(3,1)) + (triang(1,1)-triang(3,1))*(x(2)-triang(3,2)) )/detT
    lambda(3)=1-lambda(1)-lambda(2)
    !print*,"lambda",lambda(1),lambda(2),lambda(3)
    if ( (TOL < lambda(1)) .and. (lambda(1) < 1-TOL) .and. (TOL < lambda(2)) .and. (lambda(2) < 1-TOL) &
         .and.(TOL < lambda(3)) .and. (lambda(3) < 1-TOL) )then
       ins = .true.
    end if

  end subroutine Inside

  !subroutine which returns ins true, when x is inside the triangle
  ! with the ending points, like edges and vertices
  subroutine Inside2(x,triang,ins)
    real, dimension(1:2),intent(in)			:: x
    real, dimension(1:3,1:2),intent(in)			:: triang
    logical, intent(out)					:: ins
    real, dimension(1:3)					:: lambda
    real							:: detT
    real, parameter					:: TOL=0.000005

    ins=.false.
    detT = ( ( (triang(1,1)-triang(3,1))*(triang(2,2)-triang(3,2)) ) - ( (triang(1,2)-triang(3,2))*(triang(2,1)-triang(3,1)) ) )
    !  !print*,"detT",detT
    lambda(1)=((triang(2,2)-triang(3,2))*(x(1)-triang(3,1)) + (triang(3,1)-triang(2,1))*(x(2)-triang(3,2)) )/detT
    lambda(2)=((triang(3,2)-triang(1,2))*(x(1)-triang(3,1)) + (triang(1,1)-triang(3,1))*(x(2)-triang(3,2)) )/detT
    lambda(3)=1-lambda(1)-lambda(2)
    !print*,"lambda",lambda(1:3)
    ! write(20,*),"lambda",lambda(1:3)
    if ( (-TOL <= lambda(1)) .and. (lambda(1) <= (1+TOL)) .and. (-TOL <= lambda(2)) .and. (lambda(2) <= (1+TOL))  &
         .and.(-TOL <= lambda(3)) .and. (lambda(3) <= (1+TOL)) ) then
       ins = .true.
    end if

  end subroutine Inside2


  ! from numbers of vertices returns third vertex
  subroutine GetOtherVer(V,vr1,vr2)
    character(len=2),intent(in)		:: V
    integer,intent(out)			:: vr1,vr2

    if (V=="VA") then
       vr1=2
       vr2=3
    else
       if (V=="VB") then
          vr1=1
          vr2=3
       else
          vr1=1
          vr2=2
       end if
    end if

  end subroutine GetOtherVer


  !compare vertex and edge and returns inn, when it's true, they're on the same edge
  subroutine CompareEdge(V,E,inn)
    character(len=2),intent(in)		:: V,E
    logical, intent(out)			:: inn

    inn=.false.
    if (V=="VA") then
       if ( (E=="AB") .or. (E=="CA") ) then
          inn=.true.
          return
       end if
    end if
    if (V=="VB") then
       if ( (E=="AB") .or. (E=="BC") ) then
          inn=.true.
          return
       end if
    end if
    if (V=="VC") then
       if ( (E=="BC") .or. (E=="CA") ) then
          inn=.true.
          return
       end if
    end if

  end subroutine CompareEdge

  !to particular vertex get his edges
  subroutine GetEdges(D,V,W)
    character(len=2),intent(in)		:: D
    character(len=2),intent(out)		:: V,W

    if (D=="VA") then
       V="AB"
       W="CA"
    end if
    if (D=="VB") then
       V="AB"
       W="BC"
    end if
    if (D=="VC") then
       V="BC"
       W="CA"
    end if

  end subroutine GetEdges


  ! from two edges returns shared vertex
  ! made from letters
  subroutine GetSharedVertex(D,E,vrt)
    character(len=2),intent(in)		:: D,E
    integer, intent(out)			:: vrt

    if (D==E) then
       print*, "Error in GetSharedVertex, two same vertices"
       return
    end if
    if (D=="AB") then
       if (E=="BC") then
          vrt=2
       else
          !E="CA"
          vrt=1
       end if
       return
    end if
    if (D=="BC") then
       if (E=="AB") then
          vrt=2
       else
          !E=CA
          vrt=3
       end if
       return
    end if
    if (D=="CA") then
       if (E=="AB") then
          vrt=1
       else
          !E=BC
          vrt=3
       end if
       return
    end if
    !print*,"Bad task in GEtSharedVertex, unknown vertices"

  end subroutine GetSharedVertex


  ! from two edges returns shared vertex
  ! made from numbers
  subroutine GetSharedVertex2(D,E,vrt)
    integer,intent(in)			:: D,E
    integer, intent(out)			:: vrt


    if (D==E) then
       print*, "Error in GetSharedVertex, two same vertices"
       return
    end if

    if (D==1) then
       if (E==3) then
          vrt=1
       else
          !E="CA"
          vrt=2
       end if
       return
    end if
    if (D==3) then
       if (E==1) then
          vrt=1
       else
          !E=CA
          vrt=3
       end if
       return
    end if
    if (D==2) then
       if (E==1) then
          vrt=2
       else
          !E=BC
          vrt=3
       end if
       return
    end if
    !print*,"Bad task in GEtSharedVertex, unknown vertices"
  end subroutine GetSharedVertex2

  subroutine GetThirdVertex(vr1,vr2,vr3)
    integer,intent(in)	:: vr1,vr2
    integer,intent(out)	:: vr3

    select case(vr1+vr2)
    case(5)
       vr3=1
    case(4)
       vr3=2
    case(3)
       vr3=3
    case default
       print*,"Error in GetThirdVertex"
    end select
  end subroutine GetThirdVertex

  subroutine GetThirdVertex2(vr1,vr2,vr3,vr4,two)

    integer,intent(in)	:: vr1,vr2

    integer,intent(out)	:: vr3, vr4

    logical, intent(out)	:: two

    two=.false.
    if (vr1==vr2)then
       two=.true.
       select case(vr1)
       case(1)
          vr3=2
          vr4=3
       case(2)
          vr3=1
          vr4=3

       case(3)
          vr3=1
          vr4=2
       end select
       return
    end if


    select case(vr1+vr2)

    case(5)

       vr3=1

    case(4)

       vr3=2

    case(3)

       vr3=3

    case default

       print*,"Error in GetThirdVertex"

    end select

  end subroutine GetThirdVertex2



  !subroutine, which makes the special node node number one
  subroutine OneIsOne(ins,nodesi)
    type(pos),intent(in)				:: ins
    real, dimension(1:3,1:2),intent(inout)	:: nodesi
    integer					:: j


    !which one is the special one
    j=1
    do
       if (ins%is(j) .eqv. .true.) then
          exit
       end if
       j=j+1
    end do
    !if j>3 then some error must have occured, cause we have just 3 vertices
    if (j>3) then
       print*, "ERROR in OneIsOne"
    end if
    !if j isn't one then we should swap nodes as if the special one was node number one
    if ( j/=1 ) then
       call SwapNodes(nodesi,j)
    end if

  end subroutine OneIsOne



  !swap j nodes with the first nodes

  subroutine SwapNodes(nodesi,j)
    real, dimension(1:3,1:2),intent(inout)	:: nodesi
    real, dimension(1:2)				:: nodes
    integer, intent(in)				:: j


    nodes(1:2) = nodesi(1,1:2)
    nodesi(1,1:2) = nodesi(j,1:2)
    nodesi(j,1:2) = nodes(1:2)
    !  !print*,"nodesi(1,1:2)",nodesi(1,1:2)
    !  !print*,"nodesi(j,1:2)  j",j,nodesi(j,1:2)

  end subroutine SwapNodes



  !swap j nodes with k nodes

  subroutine SwapNodesKJ(nodesi,j,k)
    real, dimension(1:3,1:2),intent(inout)	:: nodesi
    real, dimension(1:2)				:: nodes
    integer, intent(in)				:: j,k

    nodes(1:2) = nodesi(j,1:2)
    nodesi(j,1:2) = nodesi(k,1:2)
    nodesi(k,1:2) = nodes(1:2)

  end subroutine SwapNodesKJ



  !subroutine which finds out if x is in line xa,xb

  subroutine InLine(x,xa,xb,line)
    real, dimension(1:2),intent(in)	:: x,xa,xb
    logical,intent(out)			:: line
    real, parameter			:: TOL=0.000005
    real, dimension(1:2)			:: s

    line=.false.
    !xa=xb
    if ( (abs(xa(1)-xb(1)) <= TOL) .and. (abs(xa(2)-xb(2)) <= TOL) )then
       !print*, "xa,xb are the same points, it's not a line"
       return
    end if
    !x=xa
    if ( (abs(xa(1)-x(1)) <= TOL) .and. (abs(xa(2)-x(2)) <= TOL) )then
       !print*,"xa=x"
       line = .true.
       return
    end if
    !x=xb
    if ( (abs(xb(1)-x(1)) <= TOL) .and. (abs(xb(2)-x(2)) <= TOL) ) then
       line = .true.
       !print*,"xb=x"
       return
    end if
    !xa and xb have the same x-coordinate
    if (abs(xa(1)-xb(1)) <= TOL) then
       !x have the same x-coord
       if  (abs(xa(1)-x(1)) <= TOL) then
          s(2) = ((xa(2)-x(2))/(xb(2) - xa(2)))
          !0<=s<=1 ---> true, else false, anyway it's done
          if ( (TOL<s(2)).and. (s(2)<1-TOL) ) then
             !print*,"s(2)",s(2)
             line=.true.
             return
          end if
       end if
       return
    end if
    !xa and xb have the same y-coordinate
    if (abs(xa(2)-xb(2)) <= TOL) then
       !x have the same x-coord
       if  (abs(xa(2)-x(2)) <= TOL) then
          s(1) = ((xa(1)-x(1))/(xb(1) - xa(1)))
          !0<=s<=1 ---> true, else false, anyway it's done
          if ( (TOL<s(1)).and. (s(1)<1-TOL) ) then
             !print*,"s(1)",s(1)
             line=.true.
          end if
       end if
       return
    end if
    !NO same coordinates of any point
    s(1) = (xa(1)-x(1))/(xb(1) - xa(1))
    s(2) = (xa(2)-x(2))/(xb(2) - xa(2))
    !!print*,"xa(2)-x(2)",xa(2)-x(2),"xb(2) - xa(2)",xb(2) - xa(2),"s(2)",s(2)
    !!print*,"xa",xa
    !!print*,"xb",xb
    !!print*,"x",x
    !!print*,"s",s
    if (( abs(s(1)-s(2))<=TOL ).and. (TOL<s(2)).and. (s(2)<1-TOL) .and. (TOL<s(1)).and. (s(1)<1-TOL) ) then
       line=.true.
    end if

  end subroutine InLine



  ! distance of x and y

  real function Dist(x,y)
    real, dimension(1:2)	:: x,y

    dist = sqrt(  ( x(1)-y(1) )**2+( x(2)-y(2) )**2 )

  end function Dist



  !find maximum of element`s diameteres

  real function MaxDiam(grid)
    real 			:: max_diam
    integer 		:: i
    type(mesh)		:: grid
    class(element), pointer:: elem

    max_diam=0
    do i=1,grid%nelem
       !!print*,"Element",i,"has diam", grid%elem(i)%diam
       if ((grid%elem(i)%diam)>(max_diam)) then
          max_diam = grid%elem(i)%diam
       end if
    end do
    !!print*, "maximum diameter is", max_diam
    MaxDiam=max_diam

  end function MaxDiam

  real function MinDiam(grid)
    real 			:: min_diam
    integer 		:: i
    type(mesh)		:: grid
    class(element), pointer:: elem

    min_diam=10
    do i=1,grid%nelem
       !!print*,"Element",i,"has diam", grid%elem(i)%diam
       if ((grid%elem(i)%diam)<(min_diam)) then
          min_diam = grid%elem(i)%diam
       end if
    end do
    !!print*, "maximum diameter is", max_diam
    MinDiam=min_diam
  end function MinDiam


  !> subroutine which has two triangles coordinates and barycentric coord. of
  !> a point in tri1 and returns bar. coord of the point in tri2
  subroutine GetBarCoord(tri1, tri2, lambda, lambda2)
    real, dimension(1:3,1:2),intent(in)	:: tri1, tri2
    real, dimension(1:2),intent(in)	:: lambda
    real, dimension(1:2),intent(out)	:: lambda2
    real, dimension(1:3)			:: lambda1
    real, dimension(1:2)			:: y
    real					:: detT, Dx, Dy
    integer 				:: i
    real,parameter			:: TOL=0.000001

    !print*,"GetBarCoord called"
    lambda1(1:2)=lambda(1:2)
    lambda1(3)=1-lambda(1)-lambda(2)
    do i=1,3
       if ( (lambda1(i)>1+TOL) .or. (lambda1(i) < -TOL) ) then
          print*, "Error in GetBarCoord, lambda out of range, point out of triangle"
          return
       end if
    end do

    !determinant of the transformation
    detT = ( ( (tri2(1,1)-tri2(3,1))*(tri2(2,2)-tri2(3,2)) ) - ( (tri2(1,2)-tri2(3,2))*(tri2(2,1)-tri2(3,1)) ) )

    !print*,"detT",detT
    !physical coordinates
    y(1) = lambda1(1)*tri1(1,1) + lambda1(2)*tri1(2,1) + lambda1(3)*tri1(3,1)
    y(2) = lambda1(1)*tri1(1,2) + lambda1(2)*tri1(2,2) + lambda1(3)*tri1(3,2)

    !KS  lambda2(1)=((tri2(2,2)-tri2(3,2))*(y(1)-tri2(3,1)) + (tri2(3,1)-tri2(2,1))*(y(2)-tri2(3,2)) )/detT
    !KS  lambda2(2)=((tri2(3,2)-tri2(1,2))*(y(1)-tri2(3,1)) + (tri2(1,1)-tri2(3,1))*(y(2)-tri2(3,2)) )/detT
    !KS  !lambda2(3)=1-lambda(1)-lambda(2)


    Dx = ((tri2(2,2)-tri2(3,2))*(y(1)-tri2(3,1)) + (tri2(3,1)-tri2(2,1))*(y(2)-tri2(3,2)) )/detT
    Dy = ((tri2(3,2)-tri2(1,2))*(y(1)-tri2(3,1)) + (tri2(1,1)-tri2(3,1))*(y(2)-tri2(3,2)) )/detT
    !KSlambda2(1)= Dx
    !KSlambda2(2)= Dy
    lambda2(1)= Dy
    lambda2(2)= 1. - Dx - Dy
    !lambda2(3)=1-lambda(1)-lambda(2)


    ! USE subrotine NodeInTriangle in errorFlux.f90

  end subroutine GetBarCoord



  !subroutine which puts i to current and curren to the top of linked list

  subroutine Put(i,first,current)
    integer,intent(in)		:: i
    type(link),pointer,intent(out)	:: current
    type(link),pointer,intent(inout)	:: first

   !allocate(current)
    current%i=i
    current%next => first
    first => current
  end subroutine Put



  !copy inter to linked list

  subroutine Copy(K,inter,triangles,NumTri,firs,curren)
    type(intersect),intent(in)		:: inter
    integer, intent(in)			:: K
    !type(link),pointer,intent(inout)	:: f
    integer,dimension(:),intent(in)	:: triangles
    integer,intent(in)			:: NumTri
    type(linter),pointer,intent(inout)	:: firs,curren
    integer				:: j
    type(link),pointer			:: h

    do j=1,NumTri
       !if (associated(f).eqv..true.)then
       allocate(curren)
       curren%K = K
       curren%i = triangles(j)
       !print*,"curren%i",curren%i
       curren%NumPnt = inter%NumPnt(triangles(j))
       !print*,"curren%NumPnt",curren%NumPnt
       curren%IntPnt(1:inter%NumPnt(triangles(j)),1:2) = inter%IntPnt(triangles(j),1:inter%NumPnt(triangles(j)),1:2)
       !print*,"current%IntPnt",curren%IntPnt(1:inter%NumPnt(triangles(j)),1:2)
       curren%next => firs
       !!print*,"next",curren%next%i
       firs => curren
       ! else
       !exit
       ! end if
       !	if (associated(f%next).eqv..true.) then
       !		h=>f%next
       !		deallocate(f)
       !		allocate(f)
       !		f=>h
       !	else
       !		exit
       !	end if
    end do

  end subroutine Copy



  !subroutine which fills IntPnt and NumPnt with points x1,x2, according to edgi,

  ! which is number of edge on which the intersections lies

  subroutine FillInOrder(edgi,NumPnt,IntPnt,x1,x2)
    integer,dimension(1:2),intent(in)	:: edgi
    integer,intent(inout)			:: NumPnt
    real,dimension(1:6,1:2),intent(inout) :: IntPnt
    real,dimension(1:2),intent(in)	:: x1,x2
    !print*,"edgi",edgi(1:2)

    if((edgi(1)==2).or.(edgi(2)==2))then
       if (edgi(1)<edgi(2)) then
          NumPnt=NumPnt+2
          IntPnt(NumPnt-1,1:2)=x1(1:2)
          IntPnt(NumPnt,1:2)=x2(1:2)
       else
          NumPnt=NumPnt+2
          IntPnt(NumPnt-1,1:2)=x2(1:2)
          IntPnt(NumPnt,1:2)=x1(1:2)
       end if
    else
       if (edgi(1)==3) then
          NumPnt=NumPnt+2
          IntPnt(NumPnt-1,1:2)=x1(1:2)
          IntPnt(NumPnt,1:2)=x2(1:2)
       else
          !print*,"x2,pak x1",x2,x1
          NumPnt=NumPnt+2
          IntPnt(NumPnt-1,1:2)=x2(1:2)
          IntPnt(NumPnt,1:2)=x1(1:2)
          !print*,"IntPnt",IntPnt(NumPnt-1:NumPnt,1:2)
       end if
    end if

  end subroutine FillInOrder


  !compute area of triangle tri
  subroutine AreaTri(tri,area)
    real,dimension(1:3,1:2),intent(in)	:: tri
    real,intent(out)			:: area
    real					:: s,a,b,c
    a = Dist(tri(1,1:2),tri(2,1:2))
    b = Dist(tri(2,1:2),tri(3,1:2))
    c = Dist(tri(1,1:2),tri(3,1:2))
    s = ((a+b+c)/2)
    area = Sqrt(s*(s-a)*(s-b)*(s-c))
  end subroutine AreaTri

end module st_interpol

