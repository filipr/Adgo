!> space high-order  iterpolation on unstructured grids

module ama_interpol

  use main_data
  use problem_oper
  use geometry
  use mesh_oper
  use st_interpol

  implicit none

  public :: InterpolDGsolution

contains


  !> evaluate \f$ L^2\f$-projection of piecewise polynomial function from gridO to gridN,
  !> new grid: pointer gridN,
  !> old grid: pointer gridO
   subroutine InterpolDGsolution(gridN, gridO )
    type(mesh), intent(inout)	:: gridN, gridO
    class(element), pointer :: elem, elemN
    integer                :: i, j, j1, K, k1, l, m, NumTri, NumSubTri, NumPnt, o, p, q, qN
    integer                :: Qdof, Qnum, dof, dofN, max_dof, max_Qdof
    real, dimension(:,:), allocatable :: xii,xiK  ! barycentric coords in elem of gridN and grid
    real, dimension(:,:), allocatable :: phiLi,phiLK  ! values of test functions in xii, xiK
    real, dimension(:), allocatable   :: weights, w, wO, wN
    real, dimension(:,:), pointer 	:: phi
    type(intersect)			:: inter
    integer, dimension(:), allocatable	:: triangles! numbers of triangles in gridN which have
						! non-empty intersection with i element of grid
    real, dimension(1:4,1:3,1:2)	:: tris	! polygonal intersection divided into triangles
    real, dimension(1:3,1:2)		:: nodesK, nodesi ! vertices of current triangles
    real				:: val, area, sum_area
    real, dimension(1:2)		:: lambda,lambda2
    !type(link),pointer			:: f,h



    !print*,'Nelems: ',gridN%nelem, gridO%nelem

    !call PlotMesh(gridO,  'mesh-old')
    !call PlotMesh(gridN, 'mesh-new')

    allocate(triangles(1:gridN%nelem))
    allocate(inter%NumPnt(1:gridN%nelem))
    allocate(inter%done(1:gridN%nelem))
    allocate(inter%had(1:gridN%nelem))
    allocate(inter%IntPnt(1:gridN%nelem,1:6,1:2))

    do i=1, gridN%nelem
       elemN => gridN%elem(i)
       elemN%vec(rhs,:) = 0.
    enddo


    !allocating for the biggest possible size
    !max_dof=0
    !max_Qdof=0
    !do K=1,gridO%nelem
    !    elem => gridO%elem(K)
    !	Qnum = elem%Qnum
    !	print*,"max_dof",elem%dof
    !	if (elem%dof>max_dof) then
    !	    max_dof = elem%dof
    !	end if
    !	if (state%space%V_rule(Qnum)%Qdof > max_Qdof) then
    !		max_Qdof = state%space%V_rule(Qnum)%Qdof
    !	end if
    !end do
    !print*,"Qdof",Qdof
    !print*,"nbDim",nbDim

    max_dof  = max ( maxval(gridO%elem(:)%dof), maxval(gridN%elem(:)%dof) )
    max_Qdof = max( maxval(gridO%elem(:)%Qdof), maxval(gridN%elem(:)%Qdof) )

    allocate(xii(1:max_Qdof, 1:nbDim))   ! array of Qdof barycenntric coordinates
    allocate(xiK(1:max_Qdof, 1:nbDim))

    allocate( phiLi(1:max_dof, 1:max_Qdof) )     ! values of test functions in bar. coor.
    allocate( phiLK(1:max_dof, 1:max_Qdof) )
    allocate( weights(1:max_Qdof) )

    sum_area = 0.
    do K=1,gridO%nelem
       nodesK(1:3,1:2)=gridO%x(gridO%elem(K)%face(idx,1:3),1:2)
       elem => gridO%elem(K)

       ! subroutine which finds intersections of gridN and i-element of grid
       ! result are saved into inter and NumTri is number of triangles with
       ! non-empty intersection and triangles are the numbers of these triangles
       ! in the gridN
       call IntersectGridTri(gridN,gridO,K,inter,triangles,NumTri)

       do m=1,NumTri
          elemN => gridN%elem(triangles(m))
          dofN = elemN%dof

          nodesi(1:3,1:2)=gridN%x(gridN%elem(triangles(m))%face(idx,1:3),1:2)
          !select case according to number of intersection points
          ! ---> according to type of polygonal
          !making triangles from polygons
          !print*,",inter%NumPnt(triangles(m))",inter%NumPnt(triangles(m))
          select case (inter%NumPnt(triangles(m)))
          case(6)
             !4 triangles
             do o=1,4
                tris(o,1,1:2) = inter%IntPnt(triangles(m),1,1:2)
             end do
             tris(1:4,2,1:2) = inter%IntPnt(triangles(m),2:5,1:2)
             tris(1:4,3,1:2) = inter%IntPnt(triangles(m),3:6,1:2)
          case(5)
             ! 3 triangles
             do o=1,3
                tris(o,1,1:2) = inter%IntPnt(triangles(m),1,1:2)
             end do
             tris(1:3,2,1:2) = inter%IntPnt(triangles(m),2:4,1:2)
             tris(1:3,3,1:2) = inter%IntPnt(triangles(m),3:5,1:2)
          case(4)
             ! 2 triangles
             do o=1,2
                tris(o,1,1:2) = inter%IntPnt(triangles(m),1,1:2)
             end do
             tris(1:2,2,1:2) = inter%IntPnt(triangles(m),2:3,1:2)
             tris(1:2,3,1:2) = inter%IntPnt(triangles(m),3:4,1:2)
          case(3)
             ! 1 triangle
             tris(1,1:3,1:2) = inter%IntPnt(triangles(m),1:3,1:2)
          case(2)
             ! line
             cycle
          case(1)
             ! just a point, measure is 0
             cycle
          case default
             cycle
          end select

          !print*,"case ended"

          dof = elem%dof   ! number of test functions
          Qnum = elem%Qnum ! quadrature
          Qdof = elem%Qdof ! number of integ nodes

          ! Qdof = state%space%V_rule(Qnum)%Qdof
          ! !!  Qdof=max_Qdof  -- nemusi byt dobre !!
          !print*,"dof,Qdof",dof,Qdof
          phi => state%space%V_rule(Qnum)%phi(1:dof, 1:Qdof)

          !do for every triangle of the polygon
          NumSubTri = inter%NumPnt(triangles(m))-2 ! number of sub-triangles of the polygon
          if( NumSubTri >=1 ) then
             do o=1,NumSubTri
                !write(100+K, *) tris(o,1,1:2)
                !write(100+K, *) tris(o,2,1:2)
                !write(100+K, *) tris(o,3,1:2)
                !write(100+K, *) tris(o,1,1:2)
                !write(100+k, '(x)')
                !print*,'***********************************'
		!counts
		do p=1,Qdof  ! new barycentric coordinates, in K (old grid) and i (new grid)
                   lambda(1:2)=state%space%V_rule(Qnum)%lambda(p,1:2)
                   call GetBarCoord(tris(o,1:3,1:2),nodesK,lambda,lambda2(1:2))
                   !print*,"zde lambda2,p",lambda2(1:2),p
                   xiK(p,1:2)=lambda2(1:2)
                   !print*,"az tady"
                   call GetBarCoord(tris(o,1:3,1:2),nodesi,lambda,lambda2(1:2))
                   xii(p,1:2)=lambda2(1:2)

                   !write(1000+K, *) xiK(p,1:2), xii(p,1:2)

                end do

                call PHI_orthonormal(Qdof, nbDim, xii(1:Qdof, 1:nbDim), 3, dofN, phiLi(1:dofN, 1:Qdof) )
                call PHI_orthonormal(Qdof, nbDim, xiK(1:Qdof, 1:nbDim), 3, dof,  phiLK(1:dof, 1:Qdof) )


                !  if(i == 1) then
                !     write(*,'(10(a8,i5))') 'elem=',i,', Qnum =',Qnum,', Qdof=',Qdof
                !     write(*,'(20es12.4)') state%space%V_rule(Qnum)%weights(1:Qdof)
                !    write(*,'(20es12.4)') state%space%V_rule(Qnum)%lambda(1:Qdof,1)
                !     write(*,'(20es12.4)') state%space%V_rule(Qnum)%lambda(1:Qdof,2)
                !     write(*,'(20es12.4)') state%space%V_rule(Qnum)%lambda(1:Qdof,3)
                !     do j=1,dof
                !        print*,'----------------'
                !        write(*,'(20es12.4)') phi(j, 1:Qdof)
                !        write(*,'(20es12.4)') phiLi(j, 1:Qdof)
                !	     write(*,'(20es12.4)') phiLK(j, 1:Qdof)
                !          enddo
                !       endif

                ! call Eval_V_Weights(elem, weights(1:Qdof) )
                call AreaTri(tris(o,1:3,1:2),area)
                sum_area = sum_area + area
                !write(*,'(a6,i5,5es12.4)') 'area:',o, area, sum_area, elem%area, &
                !     sum_area - elem%area
                weights(1:Qdof)=state%space%V_rule(elem%Qnum)%weights(1:elem%Qdof)*area

                do qN = 1,dofN
                   do q = 1, dof
                      val = dot_product(weights(1:Qdof), phiLi(qN, 1:Qdof) * phiLK(q, 1:Qdof))

                      do k1=1,ndim
                         elemN%vec(rhs, (k1-1)*ndim + qN)  = elemN%vec(rhs, (k1-1)*ndim + qN)  &
                              + val * elem%w(0, (k1-1)*ndim + q )

                      end do
                   end do
                end do

             end do !end do for o to number of tris
          end if

       end do !end do for m to NumTri
       !stop
    enddo  ! end for K=1, gridO%nelem

    deallocate(xii, xiK, phiLi, phiLK, weights )

    ! setting of the new solution
    do i=1,gridN%nelem
       elem => gridN%elem(i)
       dof = elem%dof

       do k1 = 1,ndim
          elem%w(0, (k1-1)*dof +1 : k1*dof) &
               = matmul( elem%MassInv%Mb(1:dof, 1:dof),elem%vec(rhs, (k1-1)*dof+1 : k1*dof ) )
       enddo
    enddo

    ! test on mass conservation
    allocate(w(1:ndim), wO(1:ndim), wN(1:ndim) )
    wO(:) = 0.
    wN(:) = 0.

    ! old mesh
    do i=1,gridO%nelem
       elem => gridO%elem(i)
       call  Eval_aver_w_Elem(elem, w)
       wO(1:ndim) = wO(1:ndim) + w(1:ndim)*elem%area
       !write(*,'(a4,i5,3es12.4)') '@@@',i,w(1), elem%area, wO(1)
    enddo

    ! new mesh
    do i=1,gridN%nelem
       elem => gridN%elem(i)
       call  Eval_aver_w_Elem(elem, w)
       wN(1:ndim) = wN(1:ndim) + w(1:ndim)*elem%area
       !write(*,'(a4,i5,3es12.4)') '@@@',i,w(1), elem%area, wN(1)
    enddo

    print*,' # end subroutine InterpolDGsolution -- conservation:'
    write(*,'(a9, 10es14.6)') ' # w_old', wO(1:ndim)
    write(*,'(a9, 10es14.6)') ' # w_new', wN(1:ndim)
    write(*,'(a9, 10es14.6)') ' # diff:', wN(1:ndim)-wO(1:ndim)

    !print*,' # subroutine InterpolDGsolution -- finished'
    !stop

  end subroutine InterpolDGsolution

end module ama_interpol

