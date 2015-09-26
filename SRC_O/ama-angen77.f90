!> f90 version of the code ANGENER_77
module angener77

  use AMAdata 

contains
  
  subroutine ANGENER_77(metric, ifv,  &
       ndim, melem,nelem,mpoin,npoin,maxdeg, npoinold, nelemold,  &
       nbelm,mbelm, nbc, iorref,xold,yold,lndold,itli,tlr,iaegr,  &
       x,y,lnd,ima,iretk, lbn,ibc,itc,bx,by,ibi,icyc,  &
       xnew,ynew,wpold,ra,w,wp,supp,iae,rga,rgb,rgc,  &
       nserr,dx,dy,area,tria,ipoint,nbp,xb,yb,ibpoin,ibp,f1,f2,   &
       iflag,iba,ibb,icrack, lnd1, iae1, itrans, iwall ) 

    dimension x(mpoin),y(mpoin),lnd(melem,3),ima(melem),iretk(melem),  &
         lbn(mbelm,2),ibc(mbelm),itc(mbelm),bx(mbelm,2),by(mbelm,2),  &
         icyc(mpoin,maxdeg),xold(mpoin),yold(mpoin),lndold(melem,3),  &
         xnew(mpoin),ynew(mpoin), ibi(mbelm),iae(melem,3),  &
         ra(melem*3),wp(mpoin,ndim+1),supp(mpoin),w(1:melem,1:ndim+1),  &
         wpold(mpoin,ndim+1), lnd1(melem,3), iae1(melem,3),  &
         rga(mpoin),rgb(mpoin),rgc(mpoin),nserr(melem*3,2),  &
         dx(melem),dy(melem),area(mpoin),tria(melem),ibb(mpoin,3),  &
         xb(ipoint),yb(ipoint),ibpoin(0:nbp),  &
         ibp(mpoin,2),f1(melem),f2(melem), itrans(melem),  &
         itli(mpoin),tlr(mpoin,3),iaegr(melem,4), iba(melem, 3),   &
         iwall(10)
    logical :: metric
    integer :: i_rem, i_ins, i_swa, ichag, ichag1
    !!integer:: npoinold, nelemold
    
    write(AMA%ifig1,*)'**************************************************'
    write(AMA%ifig1,*)'****                                          ****'
    write(AMA%ifig1,*)'****       A  N  G  E  N  E  R    3.2         ****'
    write(AMA%ifig1,*)'****                                          ****'
    write(AMA%ifig1,*)'**************************************************'
    
    noit = 0
    
    call ADJAC_77(melem,nelem,mpoin,npoin,lnd,iae,maxdeg,icyc)
    
    call COINS_77(melem,nelem,mpoin,npoin,lnd,iae,nbelm,mbelm,lbn,itc)
    
    call SEEK_BOUN_77(melem,nelem,mpoin,npoin,ipoint,nbp,xb,yb,lnd,iae,  &
         ibb,ibpoin,x,y, icrack)
    
    call CYKLE_77(melem,nelem,mpoin,npoin,maxdeg,  &
         x,y,lnd,iae,icyc,  &
         lbn,nbelm,mbelm,ibp,itc,ibc)
    
    call CYKLE_BOUND_77(melem,nelem,mpoin,npoin,maxdeg,  &
         x,y,lnd,iae,icyc)
    
    call TEST_77(melem,nelem,mpoin,npoin,maxdeg,x,y,  &
         lnd,iae,icyc,lbn,nbelm,mbelm,ibp,itc)
    
    if(metric)   &
         call METRIX_77(ndim, melem,nelem,mpoin,npoin,nbelm,mbelm,  &
         maxdeg, x,y,lnd,iae,w,wp,supp,icyc,ra,surface,ibp,ifv,  &
         ipoint,nbp, ibb,ibpoin,xb,yb,tria, lbn, ibc)
    
    
    
    !     for back interpolation
    do i=1,npoin
       do j=1,ndim
          wpold(i,j) = wp(i,j+1)
       enddo
       !write(21,*) x(i), y(i), wpold(i, ndim+1)
    enddo
    

    if(metric)   &
         call ERROR1_77(ndim,melem,nelem,mpoin,npoin,maxdeg,  &
         dx(1:melem),dy(1:melem),area(1:mpoin),tria,noit,  &
         x(1:mpoin),y(1:mpoin),lnd(1:melem, 1:3),iae(1:melem,1:3),  &
         wp(1:mpoin, 1:ndim+1),icyc(1:mpoin,1:maxdeg),ra,  &
         rga(1:mpoin),rgb(1:mpoin),rgc(1:mpoin),surface,ibp,  &
         ipoint,nbp,ibb,ibpoin,xb,yb,ifv, iwa, iwall,  &
         nbelm,mbelm,lbn, ibc)
    
    !      AMA%ifig2 = 1000
    !      call WriteMetrix_77(mpoin, npoin, melem, nelem, lnd, 
    !     *     x, y, rga, rgb, rgc)
    
    write(AMA%ifig1,*)'Iteration    quality          changes  ',  &
         ' GL^2:  min        max        nelem'
    write(AMA%ifig1,*)'-----------------------------------------',  &
         '---------------------------------'
    
    call QUALITY_77(ndim,melem,nelem,mpoin,npoin,x,y,lnd,iae,noit,err,  &
         rminerr,imin,rmaxrez,ra,wp,rga,rgb,rgc)
    
    write(AMA%ifig1,99999) noit,AMA%errrez,'begin',0,AMA%glmin,AMA%glmax,nelem
    ichag = 0
    
    iter = 1
    
    call C_REP_BOUND_77(ndim,melem,nelem,mpoin,npoin,maxdeg,nbelm,mbelm,  &
         noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
         rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
         ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp)
    


    !      goto 100
    !noflc = 30
    noflc = 60
    do 99 nob=1,noflc
       i_rem = 0 
       i_ins = 0
       i_swa = 0
       
       ichag = 0
       
       ichag1 = ichag
       iter = 5
       call C_REM_BOUND_77(ndim, melem,nelem,mpoin,npoin,  &
            maxdeg,nbelm,mbelm,  &
            noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
            rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
            ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp)
       
       i_rem = i_rem + ichag - ichag1
       ichag1 = ichag
       
       iter = 10
       call C_MOVING_77(ndim, melem,nelem,mpoin,npoin,  &
            maxdeg,nbelm,mbelm,  &
            noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
            rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
            ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp)
       !         goto 100
       
       i_rem = i_rem + ichag - ichag1
       ichag1 = ichag
       
       
       iloc = 0
       do iter2 = 1,15
          if(iloc == 0) then
             iter = 1
             call C_REMOVE_77(ndim, melem,nelem,mpoin,npoin,  &
                  maxdeg,nbelm,mbelm,  &
                  noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
                  rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,  &
                  area,tria,surface,  &
                  ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp,iba,ier)
             
             i_rem = i_rem + ichag - ichag1
             ichag1 = ichag
             
             if(icha == 0) iloc = -1
          endif
          !            goto 100
          
          if(iloc == 0) then
             iter = 5
             call C_SWAPPING_77(ndim, melem,nelem,mpoin,npoin,  &
                  maxdeg,nbelm,mbelm,  &
                  noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
                  rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,  &
                  area,tria,surface,  &
                  ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp,iba,ier)
             
             i_swa = i_swa + ichag - ichag1
             ichag1 = ichag
             
             iter = 1
             call C_REP_BOUND_77(ndim, melem,nelem,mpoin,npoin,  &
                  maxdeg,nbelm,mbelm,  &
                  noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
                  rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
                  ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp)
             
          endif
       enddo
       
       
       iter = 10
       call C_MOVING_77(ndim, melem,nelem,mpoin,npoin,maxdeg,nbelm,mbelm,  &
            noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
            rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
            ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp)
       

       iloc = 0
       do iter2 =1,15
          if(iloc == 0) then
             
             ichag1 = ichag
             
             iter = 1
             call C_INSERT_77(ndim, melem,nelem,mpoin,npoin,  &
                  maxdeg,nbelm,mbelm,  &
                  noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
                  rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,  &
                  area,tria,surface,  &
                  ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp,iba,ier)
             
             i_ins = i_ins + ichag - ichag1
             ichag1 = ichag
             
             !               if(iter2 == 2)   goto 100
             
             !     bad dimension
             if(icha == 0) iloc = -1
          endif
          
          
          
          if(iloc == 0) then
             iter = 5
             call C_SWAPPING_77(ndim, melem,nelem,mpoin,npoin,  &
                  maxdeg,nbelm,mbelm,  &
                  noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
                  rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,  &
                  area,tria,surface,  &
                  ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp,iba,ier)
             
             
             i_swa = i_swa + ichag - ichag1
             ichag1 = ichag
             
             iter = 1
             call C_REP_BOUND_77(ndim, melem,nelem,mpoin,npoin,  &
                  maxdeg,nbelm,mbelm,  &
                  noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
                  rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,  &
                  surface,  &
                  ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp)
             
             
             !               if(iter2 == 2)   goto 100
             
          endif
       enddo
       
       !         goto 100
       
       iter = 10
       call C_MOVING_77(ndim, melem,nelem,mpoin,npoin,maxdeg,nbelm,mbelm,  &
            noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
            rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
            ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp)
       
       
       ichag1 = ichag
       
       iter = 5
       call C_SWAPPING_77(ndim, melem,nelem,mpoin,npoin,  &
            maxdeg,nbelm,mbelm,  &
            noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
            rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
            ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp,iba,ier)
       
       i_swa = i_swa + ichag - ichag1
       ichag1 = ichag
       
       iter = 1
       call C_REP_BOUND_77(ndim, melem,nelem,mpoin,npoin,  &
            maxdeg,nbelm,mbelm,  &
            noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
            rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
            ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp)
       
       !         if(nob == 1) goto 100
       
       !         call WriteMetrix_77(mpoin, npoin, melem, nelem, lnd, 
       !     *        x, y, rga, rgb, rgc)
       

       write(AMA%ifig1, *)'Stop AMA',nob, noflc,ichag, ichagold
       write(*, '(5(a7, i5), a8,es12.4)')  'AMA ope:',nob,  &
            'swap:',i_swa, 'ins:', i_ins,   &
            'rem:', i_rem,'tot:',ichag, 'qua:',AMA%errrez

       if( ichag == 0 .or. (ichag == ichagold .and.  &
            ichag .le. 20) )   &
            goto 101
       ichagold = ichag
99  enddo
       
101 continue
    
    !     ... no DELAUNAY
    goto 100
    !     ... to ensure the triangulation of Delaunay type
    do nob1 = 1,1
       ichag = 0
       
       iter =  30
       call C_MOVING_77(ndim, melem,nelem,mpoin,npoin,maxdeg,nbelm,mbelm,  &
            noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
            rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
            ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp)
       
       iter = 1
       if(nob1 == 1)  &
            call C_B_INSERT_77(ndim, melem,nelem,mpoin,npoin,  &
            maxdeg,nbelm,mbelm,  &
            noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
            rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,  &
            area,tria,surface,  &
            ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp,iba,ier)
       
       
       iter = 0
       call C_DELAUNAY_77(ndim, melem,nelem,mpoin,npoin,  &
            maxdeg,nbelm,mbelm,  &
            noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
            rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
            ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp,iba,ier)
       if(ichag == 0) goto 100
    enddo
    
100 continue
    
    if(AMA%ifig .ge. 0 ) then
       print*, 'Number of plotted mesh sequence: ',AMA%ifig-1
    endif
    
    call PLOT_77(melem,nelem,mpoin,npoin,x,y,lnd)
    
    call TEST_77(melem,nelem,mpoin,npoin,maxdeg,x,y,  &
         lnd,iae,icyc,lbn,nbelm,mbelm,ibp,itc)
    
    
    !     ... ordering of elements
    call SEEK_NEIGH_77(melem,mpoin,mbelm,npoin,nelem,maxdeg,lnd,iae,icyc)
    !      print*,'25'
    
    !      igra = 11
    !      open (igra,STATUS='unknown', file='matrix_old')
    !      call DRAW_MATRIX_77(melem,mpoin,mbelm,igra,nelem,lnd,iae)
    !      close(igra)
    
    !      print*,'26'
    
    call COMPUTE_BAND_77(melem,mpoin,mbelm,nelem,lnd,iae,iband)
    
    !!  simple algorithmus
    do i=1,nelem
       do j=1,3
          iae1(i,j) = 0
          lnd1(i,j) = 0
       enddo
       itrans(i) = 0
    enddo
    
    ipoc = 1
    iac_old = 1
    do j=1,3
       lnd1(ipoc,j) = lnd(iac_old,j)
    enddo
    iae1(iac_old,1) = -2
    iac_new = ipoc
    itrans(iac_new) = iac_old
    
       
    do i=2,nelem
       iac_old = itrans(iac_new)
       if(iac_old == 0) then
          print*,'iac_old is ZERO'
          stop
       endif
       do j=1,3
          ii = iae(iac_old,j)
          if(ii .gt. 0 ) then
             if(iae1(ii,1) == 0) then
                ipoc = ipoc + 1
                itrans(ipoc) = ii
                do k=1,3
                   lnd1(ipoc,k) = lnd(ii,k)
                enddo
                iae1(ii,1) = -2
             endif
          endif
       enddo
       iac_new = iac_new+1
    enddo
    
       !      do i=1,nelem
       !         write(*,'(i5,2(a1,3i5))')i,'|',lnd(i,1),lnd(i,2),lnd(i,3),
       !     *        '|',iae(i,1),iae(i,2),iae(i,3)
       !      enddo
       !      print*,'*****************************'
       !      do i=1,nelem
       !         write(*,'(i5,2(a1,3i5),a1,i5)')
       !     *        i,'|',lnd1(i,1),lnd1(i,2),lnd1(i,3),
       !     *        '|',iae1(i,1),iae1(i,2),iae1(i,3),'|',itrans(i)
       !      enddo

    call SEEK_NEIGH_77(melem,mpoin,mbelm,npoin,   &
         nelem,maxdeg,lnd1,iae1,icyc)
       !      print*,'25'
       
       !      igra = 11
       !      open (igra,STATUS='unknown', file='matrix_new')
       !      call DRAW_MATRIX_77(melem,mpoin,mbelm,igra,nelem,lnd1,iae1)
       !      close(igra)
       
       !      print*,'26'
       
    call COMPUTE_BAND_77(melem,mpoin,mbelm,nelem,lnd,iae1,iband1)


       !      do i=1,npoin
       !         write(itrix,*) x(i),y(i)
       !      enddo
       !
       !      do i=1,nelem
       !         write(itrix,*) lnd1(i,1),lnd1(i,2),lnd1(i,3)
       !      enddo
       !      
       !
       !      do  k=1,nbelm
       !         write(itrix,* ) lbn(k,1),lbn(k,2),ibc(k)
       !      enddo
       !      close(itrix)
       
       
       
    write(AMA%ifig1,*)'The final mesh:'
    write(AMA%ifig1,*)'number of elements:',nelem
    write(AMA%ifig1,*)'number of points  :',npoin
    write(AMA%ifig1,*)'number of boun. s.:',nbelm
    
       !      write (ivltx,*) npoin, nelem, nbelm, nbc
       !      write (ivltx,'(2e14.7,2i3,2e14.7,2i3)') 
       !     *     AMA%xper(1,1),AMA%xper(1,2), AMA%iper(1,1),AMA%iper(1,2),AMA%xper(2,1),AMA%xper(2,2),AMA%iper(2,1),AMA%iper(2,2)
       !      do k=1,npoin
       !         write (ivltx,*) x(k),y(k)
       !      enddo
       !      do k=1,nelem
       !cc         write (ivltx,*) lnd(k,1),lnd(k,2),lnd(k,3)
       !!     ... ordering
       !         write (ivltx,*) lnd1(k,1),lnd1(k,2),lnd1(k,3)
       !      enddo
       

       !ivltx = 131
       !do  k=1,nbelm
       !   write (ivltx,*) lbn(k,1),lbn(k,2),ibc(k)
       !enddo
       !      write (ivltx,'(a4,i6,3e12.4)') '****',AMA%numel,epsilon1,p,AMA%pos
       !      close(ivltx)
       
    !print*,'AMA: interpolation'
    call INTERPOLATION_77(ndim,mpoin,melem,nelem,npoin,  &
         nelemold,npoinold,  &
         x,y,lnd,  &
         xold,yold,lndold,wpold,itli,tlr,wp,iaegr,icrack)
    
    

       !print*,'AMA: recomputation'   !AW always within ADGFEM
       !AW if(AMA%ityp .gt. 0) then
       !AW if( ifv == 1) then
    do i=1,nelem
       do  k=1,ndim + 1  ! including degree of polynomial approximation
          w(i,k) = (wp(lnd(i,1),k)  &
               + wp(lnd(i,2),k) + wp(lnd(i,3),k))/3
       enddo
    enddo
    
    !            do ii=1,nelem
       !              i = ii
       !!     ... ordering
       !               i = itrans(ii)
       !               write (iresx,'(4e16.8)' )   ( w(i,k),k=1,ndim )
       !            enddo
       !         else
       !            do  i=1,npoin
       !               write(iresx,'(4e16.8)') ( wp(i,k+1),k=1,ndim)
       !            enddo
       !AW endif
       !     close(iresx)
       !AW endif
       
       
    !print*,'AMA: heck of "positivity"'
    ipocel = 0
    do 800 i=1,nelem
       x1 = x(lnd(i,1))
       y1 = y(lnd(i,1))
       x2 = x(lnd(i,2))
       y2 = y(lnd(i,2))
       x3 = x(lnd(i,3))
       y3 = y(lnd(i,3))
       det = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
       reps = AMA%pos*( ((x1-x3)**2 + (y1-y3)**2) +  &
            ((x2-x3)**2 + (y2-y3)**2) +  &
            ((x1-x2)**2 + (y1-y2)**2) )
       if (det .le. reps) then
          ipocel = ipocel + 1
       endif
800 enddo
    print *,'Total number of violations of positivity:',ipocel
       
       
    !print*,'AMA: heck of "shaphness"'
    itet = 0
    ipocel = 0
    do  i=1,nelem
       x1 = x(lnd(i,1))
       y1 = y(lnd(i,1))
       x2 = x(lnd(i,2))
       y2 = y(lnd(i,2))
       x3 = x(lnd(i,3))
       y3 = y(lnd(i,3))
       call POS1TEST_77(x1,y1,x2,y2,x3,y3,itet)
       if(itet == 1) ipocel = ipocel + 1
       if(itet == 1) then
          write (41,*) x1,y1
          write (41,*) x2,y2
          write (41,*) x3,y3
          write (41,*) x1,y1
          write(41,'(x)')
       endif
    enddo
    print *,'Total number of violations of sharpness: ',ipocel
       
       
88888 format(a1,4x,2i5)                                                 
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
    return
  end subroutine ANGENER_77


  subroutine  COMPUTE_BAND_77(melem,mpoin,mbelm,nelem,lnd,iae,iband)
    dimension lnd(melem,3),iae(melem,3)
      
    iband = 0
    do i=1,nelem
       do j=1,3
          ii = iae(i,j)
          if(ii .gt. 0) then
             iband = max(iband,abs(i-ii) )
          endif
       enddo
    enddo
    
    write(*,'(a6,i7,a14,i7,a8,f12.8)')   &
         'nelem=',nelem,', matrix band=',iband,  &
         ', ratio=',1.*iband/nelem
  end subroutine COMPUTE_BAND_77
  
  subroutine  DRAW_MATRIX_77(melem,mpoin,mbelm,igra,nelem,lnd,iae)
    dimension lnd(melem,3),iae(melem,3)
    
    do i=1,nelem
       write(igra,*) i, -i
       do j=1,3
          if(iae(i,j) .gt. 0) then
             write(igra,*) i, -iae(i,j)
          endif
       enddo
    enddo
  end subroutine DRAW_MATRIX_77

  subroutine SEEK_NEIGH_77(melem,mpoin,mbelm,  &
       npoin,nelem,maxdeg,lnd,iae,icyc)
    dimension lnd(melem,3),iae(melem,3),icyc(mpoin,maxdeg)
    !     ... seeking of neighbours
    do i=1,nelem
       do j=1,3
          iae(i,j) = -2
       enddo
    enddo
    
    do i=1,npoin
       icyc(i,1) = 0
    enddo
    
    !      print*,'23'
    do i=1,nelem
       do j=1,3
          k = lnd(i,j)
          icyc(k,1) = icyc(k,1) + 1
          if(icyc(k,1) .ge. maxdeg) then
             print *,'Bad dimension in ADJAC'
             print *,'maxdeg < icyc(k,1)',maxdeg,icyc(k,1)
             stop
          endif
          icyc(k,icyc(k,1)+1) = i
       enddo
    enddo
    !      print*,'24'
      
    do i=1,nelem
       do  j=1,3
          if(iae(i,j) == -2)then
             j1 = mod(j,3) +1
             do il=1,icyc(lnd(i,j),1)
                ii = icyc(lnd(i,j),il+1)
                if( ii .ne. i) then
                   do jj=1,3
                      jj1 = mod(jj,3) +1
                      if( (lnd(i,j) == lnd(ii,jj1) ).and.  &
                           (lnd(i,j1) == lnd(ii,jj))) then
                         !     print *,'**',i,j,lnd(i,j),lnd(i,j1),ii
                         iae(i,j) = ii
                         iae(ii,jj) = i
                         goto 60
                      endif
                   enddo
                endif
             enddo
60           continue
          endif
       enddo
    enddo

  end subroutine SEEK_NEIGH_77

  subroutine POS1TEST_77(x1,y1,x2,y2,x3,y3,itet)

    if(itet == -1) then
       itet1 = -1
    else
       itet1 = 0
    endif
    itet1 = 0
    
    if(itet1 == -1) then
       print *,'***  ---'
       print *,x1,y1
       print *,x2,y2
       print *,x3,y3
    endif
    itet = 0
    rl1 = ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))**0.5
    rl2 = ((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3))**0.5
    rl3 = ((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1))**0.5
    if(rl1 + rl2 .lt. (1.+AMA%pos1)*rl3) itet = 1
    if(rl2 + rl3 .lt. (1.+AMA%pos1)*rl1) itet = 1
    if(rl3 + rl1 .lt. (1.+AMA%pos1)*rl2) itet = 1
    
    if(itet1 == -1) print *,'* *',itet
  end subroutine POS1TEST_77

  subroutine POS2TEST_77(x1,y1,x2,y2,x3,y3,itet)
    real *8 x1,y1,x2,y2,x3,y3
    
    if(itet == -1) then
       itet1 = -1
    else
       itet1 = 0
    endif
    if(itet1 == -1) then
       print *,'***'
       print *,x1,y1
       print *,x2,y2
       print *,x3,y3
    endif
    itet = 0
    rl1 = ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))**0.5
    rl2 = ((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3))**0.5
    rl3 = ((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1))**0.5
    if(rl1 + rl2 .lt. (1.+AMA%pos1)*rl3) itet = 1
    if(rl2 + rl3 .lt. (1.+AMA%pos1)*rl1) itet = 1
    if(rl3 + rl1 .lt. (1.+AMA%pos1)*rl2) itet = 1
    if(itet1 == -1) print *,'* *',itet
  end subroutine POS2TEST_77


  subroutine INTERPOLATION_77(ndim,mpoin,melem,  &
       nelem,npoin,nelemold,npoinold,x,y,lnd,  &
       xold,yold,lndold,wpold,itli,tlr,wp,iaegr,icrack)
    
    dimension  x(mpoin),y(mpoin),lnd(melem,3),  &
         xold(mpoin),yold(mpoin),lndold(melem,3),  &
         itli(mpoin),tlr(mpoin,3),wpold(mpoin,ndim+1),  &
         iaegr(melem,4),wp(mpoin,ndim+1)
    
    
    !print*,'AMA: coefficients for transformation'
    !write(*,'(30i8)') ndim,mpoin,melem,nelem,npoin,nelemold,npoinold

    !     from global mesh to the local one'
    do i=1,npoin
       itli(i) = 0
       tlr(i,1) = 0.
       tlr(i,2) = 0.
       tlr(i,3) = 0.
    enddo
    
    xrmax = -10000.
    xrmin = 100000.
    yrmax = -10000.
    yrmin = 100000.
    do i=1,npoinold
       xc = xold(i)
       yc = yold(i)
       xrmax = max(xrmax,xc)
       xrmin = min(xrmin,xc)
       yrmax = max(yrmax,yc)
       yrmin = min(yrmin,yc)
    enddo
    
    sradiusmax = 0.
    do i=1,nelemold
       iaegr(i,2) =  0
       iaegr(i,3) =  0
       x1 = xold(lndold(i,1))
       y1 = yold(lndold(i,1))
       x2 = xold(lndold(i,2))
       y2 = yold(lndold(i,2))
       x3 = xold(lndold(i,3))
       y3 = yold(lndold(i,3))
       radius = max(((x1-x2)**2 + (y1-y2)**2)**0.5,    &
            ((x2-x3)**2 + (y2-y3)**2)**0.5 ,  &
            ((x1-x3)**2 + (y1-y3)**2)**0.5, sradiusmax )
       sradiusmax = max(sradiusmax, radius)
    enddo
    sradiusmax = sradiusmax *2
    
    !print*,'AMA-solve: sradiusmax =', sradiusmax
    nsqrt = nelemold**0.5 
    do i=1,nelemold
       xc = (xold(lndold(i,1)) + xold(lndold(i,2)) +   &
            xold(lndold(i,3)) )/3
       yc = (yold(lndold(i,1)) + yold(lndold(i,2)) +   &
            yold(lndold(i,3)) )/3
       sx = (xc - xrmin)/(xrmax - xrmin)
       sy = (yc - yrmin)/(yrmax - yrmin)
       if(sx .ge. 1. .or. sy .ge. 1 .or.   &
            sx .le. 0 .or. sy .le. 0.) then
          print *,xc,yc,sx,sy
          stop
       endif
       ix = sx*nsqrt + 1
       iy = sy*nsqrt + 1
       iele = (ix-1)*nsqrt + iy
       iaegr(i,1) = iele
       iaegr(iele,2) = iaegr(iele,2) + 1
    enddo
    rboxmax = max((xrmax - xrmin)/nsqrt,(yrmax - yrmin)/nsqrt )
    ishiftmax = sradiusmax/rboxmax + 1
    
    !print *,'BOXES:',nsqrt,nsqrt**2,nelemold
    !print *,sradiusmax,rboxmax,sradiusmax/rboxmax,ishiftmax
    it = -10
    do i=1,nelemold
       if(iaegr(i,1) == it) then
          print *,xold(lndold(i,1)), yold(lndold(i,1)),lndold(i,1),i
          print *,xold(lndold(i,2)), yold(lndold(i,2)),lndold(i,2)
          print *,xold(lndold(i,3)), yold(lndold(i,3)),lndold(i,3)
          print *,xold(lndold(i,1)), yold(lndold(i,1))
          print *
       endif
    enddo
    !print *,'pocet =',iaegr(it,2),it/nsqrt,mod(it,nsqrt)
    isum = 0
    do i=1,nsqrt**2
       inum = iaegr(i,2)
       iaegr(i,2) = isum + 1
       isum = isum + inum
    enddo
    
    do i=1,nelemold
       icom = iaegr(i,1)
       iord = iaegr(icom,2) + iaegr(icom,3)
       iaegr(icom,3) = iaegr(icom,3) + 1
       iaegr(iord,4) = i
    enddo
    
    do ic =1,nsqrt**2
       ip = iaegr(ic,2)
       if(ic .ne. nsqrt**2) then
          il = iaegr(ic+1,2) - 1
       else
          il = nelemold
       endif
       if(ip .gt. il .and. iaegr(ic,3) .ne. 0) then
          print *,'error for ic =',ic
          print *,'nozero'
          print *,ip,il,nsqrt**2
          stop
       endif
       do iel = ip,il
          ielem = iaegr(iel,4)
          if(iaegr(ielem,1) .ne. ic) then
             print *,'error for ic =',ic
             print *,'the elements don''t correspond'
             print *,ic,iel,ielem,iaegr(ielem,1)
             stop
          endif
       enddo
    enddo
    
    do ip=1,npoin
       !         if(mod(Ip,1000) == 1) print *,'***',ip,ishift
       ishift = 0
       x0 = x(ip)
       y0 = y(ip)
       iposit = 0
       
       !     some improvement for a unit circle with a crack
       if(icrack == 1) then
          if(y0 == 0) then
             do ipom = 1,nelem
                if(lnd(ipom,1) == ip .or. lnd(ipom,2) == ip  &
                     .or. lnd(ipom,3) == ip ) then
                   yc = (y(lnd(ipom,1))+y(lnd(ipom,2))+  &
                        y(lnd(ipom,3)))/3
                   if(yc .gt. 0) then
                      iposit = 1
                   else
                      iposit = -1
                   endif
                   goto 234
                endif
             enddo
234          continue
          endif
       endif

       
       sx = (x0 - xrmin)/(xrmax - xrmin)
       sy = (y0 - yrmin)/(yrmax - yrmin)
       iiner = 0
       if(sx .gt. 0. .and. sx .lt. 1. .and.   &
            sy .gt. 0. .and. sy .lt. 1.) then
          !     inner point
          iiner = 1
       endif
       if(iiner == 1) then
          ix = sx*nsqrt + 1
          iy = sy*nsqrt + 1
954       continue
          do lx = -ishift,ishift
             do ly = -ishift,ishift
                if( ix+lx .ge. 1 .and. ix+lx .le. nsqrt .and.  &
                     iy+ly .ge. 1 .and. iy+ly .le. nsqrt .and.  &
                     (abs(lx) == ishift .or. abs(ly) == ishift) )  &
                     then
                   
                   ic = (ix-1+lx)*nsqrt + iy+ly
                   il1 = iaegr(ic,2)
                   if(ic .ne. nsqrt**2) then
                      il2 = iaegr(ic+1,2) - 1
                   else
                      il2 = nelemold
                   endif
                   !     we go through the elements belonging to ic
                   do iel = il1,il2
                      i = iaegr(iel,4)
                      i1 = lndold(i,1) 
                      i2 = lndold(i,2) 
                      i3 = lndold(i,3) 
                      xx0 = x0
                      yy0 = y0
                      xx1 = xold(i1)
                      yy1 = yold(i1)
                      xx2 = xold(i2)
                      yy2 = yold(i2)
                      xx3 = xold(i3)
                      yy3 = yold(i3)
                      
                      if(iposit .ne. 0) then
                         yyc = (yy1 + yy2 + yy3)/3.
                         if(yyc*iposit .lt. 0) then
                            goto 235
                         endif
                      endif
                      
                      det0 = xx3*(yy1-yy2)+xx1*(yy2-yy3)+xx2*(yy3-yy1)
                      det1 = xx0*(yy1-yy2)+xx1*(yy2-yy0)+xx2*(yy0-yy1)
                      det2 = xx0*(yy2-yy3)+xx2*(yy3-yy0)+xx3*(yy0-yy2)
                      det3 = xx0*(yy3-yy1)+xx3*(yy1-yy0)+xx1*(yy0-yy3)
                      
                      
                      if(det0 .le. 0.) then
                         write(112,'(2e14.6,2i7)')xx1,yy1,i1,i
                         write(112,'(2e14.6,i7)')xx2,yy2,i2
                         write(112,'(2e14.6,i7)')xx3,yy3,i3
                         write(112,'(3e14.6)')xold(i1),yold(i1),det0
                         write(112,'(x)')
                      endif
                      epss = 1E-03
                      if(abs(det1)/det0 .lt. epss) det1 = abs(det1)
                      if(abs(det2)/det0 .lt. epss) det2 = abs(det2)
                      if(abs(det3)/det0 .lt. epss) det3 = abs(det3)
                      if(det1 .ge. 0 .and. det2 .ge.0   &
                           .and. det3 .ge. 0 .and. det0 .gt. 0.) then
                         !     the node is in this triangle
                         itli(ip) = i
                         tlr(ip,1) = det2/det0
                         tlr(ip,2) = det3/det0
                         tlr(ip,3) = det1/det0
                         goto 209
                      endif
235                   continue
                   enddo
                endif
             enddo
          enddo
          
          ishift = ishift + 1
          if(ishift .le. ishiftmax) goto 954
          iiner = 0
209       continue
       endif

       if(iiner == 0) then
          !     outer point
          rmin = 100000.
          do ie =1,nelemold
             xc = (xold(lndold(ie,1)) + xold(lndold(ie,2)) +   &
                  xold(lndold(ie,3)))/3
             yc = (yold(lndold(ie,1)) + yold(lndold(ie,2)) +   &
                  yold(lndold(ie,3)))/3
             rlen = ((x0-xc)**2 + (y0-yc)**2 )**0.5
             if(iposit .ne. 0) then
                if(yc*iposit .lt. 0) then
                   goto 236
                endif
             endif
             
             if(rlen .lt. rmin) then
                rmin = rlen
                imin = ie
             endif
236          continue
          enddo
          i = imin
          i1 = lndold(i,1) 
          i2 = lndold(i,2) 
          i3 = lndold(i,3) 
          xx0 = x0
          yy0 = y0
          xx1 = xold(i1)
          yy1 = yold(i1)
          xx2 = xold(i2)
          yy2 = yold(i2)
          xx3 = xold(i3)
          yy3 = yold(i3)
          
          det0 = xx3*(yy1-yy2)+xx1*(yy2-yy3)+xx2*(yy3-yy1)
          det1 = xx0*(yy1-yy2)+xx1*(yy2-yy0)+xx2*(yy0-yy1)
          det2 = xx0*(yy2-yy3)+xx2*(yy3-yy0)+xx3*(yy0-yy2)
          det3 = xx0*(yy3-yy1)+xx3*(yy1-yy0)+xx1*(yy0-yy3)
          !     the node is the nearest to this triangle
          itli(ip) = i
          tlr(ip,1) = det2/det0
          tlr(ip,2) = det3/det0
          tlr(ip,3) = det1/det0
          if(abs(tlr(ip,1)).gt. 1.5 .or. abs(tlr(ip,2)).gt. 1.5  &
               .or. abs(tlr(ip,3) ) .gt. 1.5 ) then
             tlr(ip,1) = 1./3
             tlr(ip,2) = 1./3
             tlr(ip,3) = 1./3
          endif
       endif
    enddo
    
    do i=1,npoin
       ie = itli(i)
       do k=1,ndim+1   ! including degree of polynomial approximation
          wp(i,k) = tlr(i,1)*wpold(lndold(ie,1),k) +  &
               tlr(i,2)*wpold(lndold(ie,2),k) +  &
               tlr(i,3)*wpold(lndold(ie,3),k) 
       enddo

!         wp(i,1) = tlr(i,1)*wpold(lndold(ie,1),ndim+1) +
!     *        tlr(i,2)*wpold(lndold(ie,2),ndim+1) +
!     *        tlr(i,3)*wpold(lndold(ie,3),ndim+1) 
       
    enddo

  end subroutine INTERPOLATION_77


    subroutine C_DELAUNAY_77(ndim, melem,nelem,mpoin,npoin,  &
         maxdeg,nbelm,mbelm,  &
         noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
         rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
         ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp,iba,ier)
      
      dimension x(mpoin),y(mpoin),lnd(melem,3),iae(melem,3),  &
           lbn(mbelm,2),ibc(mbelm),itc(mbelm),  &
           icyc(mpoin,maxdeg),  &
           ra(melem*3),wp(mpoin,ndim+1),  &
           rga(mpoin),rgb(mpoin),rgc(mpoin),nserr(melem*3,2),  &
           dx(melem),dy(melem),area(mpoin),tria(melem),  &
           xb(ipoint),yb(ipoint),ibp(mpoin,2),ibb(mpoin,3),ibpoin(0:nbp)
      integer iba(melem,3)
      
      err1 = 0.
      err1rez = 0.
      do 30 icy =1,iter 
         noit = noit + 1
         call DELAUNAY_77(ndim, melem,nelem,mpoin,npoin,x,y,lnd,iae,icha,  &
              ra,wp,rga,rgb,rgc,iba,nbelm,mbelm,lbn,ibc,itc,  &
              maxdeg,icyc,ibp )
         call CYKLE_REP_77(melem,nelem,mpoin,npoin,maxdeg,  &
              x,y,lnd,iae,icyc,lbn,nbelm,mbelm,ibp)
         call CYKLE_BOUND_77(melem,nelem,mpoin,npoin,maxdeg,  &
              x,y,lnd,iae,icyc)
         call QUALITY_77(ndim, melem,nelem,mpoin,npoin,x,y,lnd,iae,noit,  &
              err,rminerr,imin,rmaxrez,ra,wp,  &
              rga,rgb,rgc)
         write(AMA%ifig1,99999)noit,AMA%errrez,'delaunay',icha,AMA%glmin,AMA%glmax,nelem
         ichag = ichag + icha
         if(icha .le. 1 ) return
30    enddo
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
      return
    end subroutine C_DELAUNAY_77

    subroutine C_SWAPPING_77(ndim, melem,nelem,mpoin,npoin,  &
         maxdeg,nbelm,mbelm,  &
         noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
         rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
         ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp,iba,ier)
      
      dimension x(mpoin),y(mpoin),lnd(melem,3),iae(melem,3),  &
           lbn(mbelm,2),ibc(mbelm),itc(mbelm),  &
           icyc(mpoin,maxdeg),  &
           ra(melem*3),wp(mpoin,ndim+1),  &
           rga(mpoin),rgb(mpoin),rgc(mpoin),nserr(melem*3,2),  &
           dx(melem),dy(melem),area(mpoin),tria(melem),  &
           xb(ipoint),yb(ipoint),ibp(mpoin,2),ibb(mpoin,3),ibpoin(0:nbp)
      integer iba(melem,3)

      err1 = 0.
      err1rez = 0.
      do  icy =1,iter 
         noit = noit + 1

         call SWAPPING_77(ndim, melem,nelem,mpoin,npoin,x,y,lnd,iae,icha,  &
              ra,wp,rga,rgb,rgc,iba,nbelm,mbelm,lbn,ibc,itc,  &
              maxdeg,icyc,ibp )

         call CYKLE_REP_77(melem,nelem,mpoin,npoin,maxdeg,  &
              x,y,lnd,iae,icyc,lbn,nbelm,mbelm,ibp)

         call CYKLE_BOUND_77(melem,nelem,mpoin,npoin,maxdeg,  &
              x,y,lnd,iae,icyc)
         call QUALITY_77(ndim, melem,nelem,mpoin,npoin,x,y,lnd,iae,noit,  &
              err,rminerr,imin,rmaxrez,ra,wp,  &
              rga,rgb,rgc)
         write(AMA%ifig1,99999)noit,AMA%errrez,'swapping',icha,AMA%glmin,AMA%glmax,nelem
         ichag = ichag + icha
         if(icha .le. 1 ) return
30    enddo
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
      return
    end subroutine C_SWAPPING_77

    subroutine C_MOVING_77(ndim, melem,nelem,mpoin,npoin,  &
         maxdeg,nbelm,mbelm,  &
         noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
         rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
         ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp)
      
      dimension x(mpoin),y(mpoin),lnd(melem,3),iae(melem,3),  &
           lbn(mbelm,2),ibc(mbelm),itc(mbelm),  &
           icyc(mpoin,maxdeg),  &
           ra(melem*3),wp(mpoin,ndim+1),  &
           rga(mpoin),rgb(mpoin),rgc(mpoin),nserr(melem*3,2),  &
           dx(melem),dy(melem),area(mpoin),tria(melem),  &
           xb(ipoint),yb(ipoint),ibp(mpoin,2),ibb(mpoin,3),ibpoin(0:nbp)

      err1 = 0.
      err1rez = 1.E+35

      do 40 icy =1,iter
         noit = noit + 1
         call MOVING_77(ndim, melem,nelem,mpoin,npoin,maxdeg,x,y,lnd,iae,  &
              ra,icyc,wp,rga,rgb,rgc,noit,ipoint,nbp,xb,yb,ibb,  &
              ibpoin,ichag, mbelm,nbelm,ibp)
         call QUALITY_77(ndim, melem,nelem,mpoin,npoin,x,y,lnd,iae,noit,  &
              err,rminerr,imin,rmaxrez,ra,wp,  &
              rga,rgb,rgc)
         write(AMA%ifig1,99998) noit,AMA%errrez,'moving',AMA%glmin,AMA%glmax,nelem
         if( err1rez - AMA%errrez .lt. 1E-06 .or. icy == iter) then
            write(AMA%ifig1,99998) noit,AMA%errrez,'moving',AMA%glmin,AMA%glmax,nelem
            return
         endif

         err1rez = AMA%errrez
40    enddo
99998 format (i6,2x,1(e16.8),2x,a9,5x,2x,2e12.4,i8)
      return
    end subroutine C_MOVING_77

    subroutine C_REM_BOUND_77(ndim, melem,nelem,mpoin,npoin,  &
         maxdeg,nbelm,mbelm,  &
         noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
         rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
         ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp)
      
      dimension x(mpoin),y(mpoin),lnd(melem,3),iae(melem,3),  &
           lbn(mbelm,2),ibc(mbelm),itc(mbelm),  &
           icyc(mpoin,maxdeg),  &
           ra(melem*3),wp(mpoin,ndim+1),  &
           rga(mpoin),rgb(mpoin),rgc(mpoin),nserr(melem*3,2),  &
           dx(melem),dy(melem),area(mpoin),tria(melem),  &
           xb(ipoint),yb(ipoint),ibp(mpoin,2),ibb(mpoin,3),ibpoin(0:nbp)

      do 20 it =1,iter
         noit = noit + 1
         call REM_BOUND_77(ndim, melem,nelem,mpoin,npoin,  &
              maxdeg,x,y,lnd,iae,  &
              ra,icyc,rminerr,imin,rmaxrez,icha,nserr,  &
              wp,lbn,itc,ibc,nbelm,mbelm,rga,rgb,rgc,ibp,ibb,ibpoin,  &
              ipoint,nbp)
         call QUALITY_77(ndim,melem,nelem,mpoin,npoin,x,y,lnd,iae,noit,  &
              err,rminerr,imin,rmaxrez,ra,wp,  &
              rga,rgb,rgc)
         write(AMA%ifig1,99999)noit,AMA%errrez,'rem_bou',icha,AMA%glmin,AMA%glmax,nelem
         ichag = ichag + icha
         if(icha .le. 1) return
20    enddo
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
      return
    end subroutine C_REM_BOUND_77

    subroutine C_REMOVE_77(ndim, melem,nelem,mpoin,npoin,  &
         maxdeg,nbelm,mbelm,  &
         noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
         rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
         ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp,iba,ier)
      
      dimension x(mpoin),y(mpoin),lnd(melem,3),iae(melem,3),  &
           lbn(mbelm,2),ibc(mbelm),itc(mbelm),  &
           icyc(mpoin,maxdeg),  &
           ra(melem*3),wp(mpoin,ndim+1),  &
           rga(mpoin),rgb(mpoin),rgc(mpoin),nserr(melem*3,2),  &
           dx(melem),dy(melem),area(mpoin),tria(melem),  &
           xb(ipoint),yb(ipoint),ibp(mpoin,2),ibb(mpoin,3),ibpoin(0:nbp)
      integer iba(melem,3)

      do 30 icy =1,iter
         noit = noit + 1
         call REMOVE_77(ndim,melem,nelem,mpoin,npoin,maxdeg,x,y,lnd,iae,  &
              ra,icyc,rminerr,imin,rmaxrez,icha,nserr,  &
              wp,lbn,ibc,itc,nbelm,mbelm,rga,rgb,rgc,ipoint,nbp,xb,yb,  &
              ibb,ibpoin,ibp)
         call QUALITY_77(ndim,melem,nelem,mpoin,npoin,x,y,lnd,iae,noit,  &
              err,rminerr,imin,rmaxrez,ra,wp,  &
              rga,rgb,rgc)
         write(AMA%ifig1,99999)noit,AMA%errrez,'remove',icha,AMA%glmin,AMA%glmax,nelem
         ichag = ichag + icha
         if(icha .le. 1) return
30    enddo
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
      return
    end subroutine C_REMOVE_77

    subroutine C_INSERT_77(ndim,melem,nelem,mpoin,npoin,  &
         maxdeg,nbelm,mbelm,  &
         noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
         rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
         ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp,iba,ier)
      
      dimension x(mpoin),y(mpoin),lnd(melem,3),iae(melem,3),  &
           lbn(mbelm,2),ibc(mbelm),itc(mbelm),  &
           icyc(mpoin,maxdeg),  &
           ra(melem*3),wp(mpoin,ndim+1),  &
           rga(mpoin),rgb(mpoin),rgc(mpoin),nserr(melem*3,2),  &
           dx(melem),dy(melem),area(mpoin),tria(melem),  &
           xb(ipoint),yb(ipoint),ibp(mpoin,2),ibb(mpoin,3),ibpoin(0:nbp)
      integer iba(melem,3)

      do 30 icy =1,iter
         noit = noit + 1
         call INSERT_77(ndim,melem,nelem,mpoin,npoin,maxdeg,x,y,lnd,iae,  &
              ra,icyc,rminerr,imin,rmaxrez,icha,nserr,  &
              wp,rga,rgb,rgc,lbn,itc,ibc,nbelm,mbelm,  &
              ipoint,nbp,xb,yb,ibb,ibpoin,  &
              ibp,ier)
         call CYKLE_REP_77(melem,nelem,mpoin,npoin,maxdeg,  &
              x,y,lnd,iae,icyc,lbn,nbelm,mbelm,ibp)
         call CYKLE_BOUND_77(melem,nelem,mpoin,npoin,maxdeg,  &
              x,y,lnd,iae,icyc)
         call QUALITY_77(ndim, melem,nelem,mpoin,npoin,x,y,lnd,iae,noit,  &
              err,rminerr,imin,rmaxrez,ra,wp,  &
              rga,rgb,rgc)
         write(AMA%ifig1,99999)noit,AMA%errrez,'insert',icha,AMA%glmin,AMA%glmax,nelem
         ichag = ichag + icha
         if(icha .le. 1) return
30    enddo
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
      return
    end subroutine C_INSERT_77

    subroutine C_REP_BOUND_77(ndim,melem,nelem,mpoin,npoin,  &
         maxdeg,nbelm,mbelm,  &
         noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
         rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
         ipoint,nbp,xb,yb,ibb,ibpoin, ichag, icha,ibp)

      dimension x(mpoin),y(mpoin),lnd(melem,3),iae(melem,3),  &
           lbn(mbelm,2),ibc(mbelm),itc(mbelm),  &
           icyc(mpoin,maxdeg),  &
           ra(melem*3),wp(mpoin,ndim+1),  &
           rga(mpoin),rgb(mpoin),rgc(mpoin),nserr(melem*3,2),  &
           dx(melem),dy(melem),area(mpoin),tria(melem),  &
           xb(ipoint),yb(ipoint),ibp(mpoin,2),ibb(mpoin,3),ibpoin(0:nbp)

      do 30 icy =1,iter
         noit = noit + 1
         call REP_BOUND_77(ndim,melem,nelem,mpoin,npoin,  &
              maxdeg,x,y,lnd,iae,  &
              ra,icyc,rminerr,imin,rmaxrez,icha,nserr,  &
              wp,rga,rgb,rgc,lbn,itc,ibc,nbelm,mbelm,  &
              ipoint,nbp,xb,yb,ibb,ibpoin,  &
              ibp)
         !     the subroutine REP_BOUND is in insert.f
         call QUALITY_77(ndim,melem,nelem,mpoin,npoin,x,y,lnd,iae,noit,  &
              err,rminerr,imin,rmaxrez,ra,wp,  &
              rga,rgb,rgc)
         write(AMA%ifig1,99999)noit,AMA%errrez,'rep_boun',icha,AMA%glmin,AMA%glmax,nelem
         ichag = ichag + icha
30    enddo
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
      return
    end subroutine C_REP_BOUND_77
    
    subroutine C_B_INSERT_77(ndim,melem,nelem,mpoin,npoin,  &
         maxdeg,nbelm,mbelm,  &
         noit, iter, x,y,lnd,iae,wp,icyc,ra,  &
         rga,rgb,rgc,lbn,ibc,itc,nserr,dx,dy,area,tria,surface,  &
         ipoint,nbp,xb,yb,ibb,ibpoin,ichag, icha,ibp,iba,ier)
      
      dimension x(mpoin),y(mpoin),lnd(melem,3),iae(melem,3),  &
           lbn(mbelm,2),ibc(mbelm),itc(mbelm),  &
           icyc(mpoin,maxdeg),  &
           ra(melem*3),wp(mpoin,ndim+1),  &
           rga(mpoin),rgb(mpoin),rgc(mpoin),nserr(melem*3,2),  &
           dx(melem),dy(melem),area(mpoin),tria(melem),  &
           xb(ipoint),yb(ipoint),ibp(mpoin,2),ibb(mpoin,3),ibpoin(0:nbp)
      integer iba(melem,3)
      
      do 30 icy =1,iter
         noit = noit + 1
         call INSERT_BOUND_77(ndim,melem,nelem,mpoin,npoin,maxdeg,  &
              x,y,lnd,iae,  &
              ra,icyc,rminerr,imin,rmaxrez,icha,nserr,  &
              wp,rga,rgb,rgc,lbn,itc,ibc,nbelm,mbelm,  &
              ipoint,nbp,xb,yb,ibb,ibpoin,  &
              ibp,ier)
         call CYKLE_REP_77(melem,nelem,mpoin,npoin,maxdeg,  &
              x,y,lnd,iae,icyc,lbn,nbelm,mbelm,ibp)
         call CYKLE_BOUND_77(melem,nelem,mpoin,npoin,maxdeg,  &
              x,y,lnd,iae,icyc)
         call QUALITY_77(ndim, melem,nelem,mpoin,npoin,x,y,lnd,iae,noit,  &
              err,rminerr,imin,rmaxrez,ra,wp,  &
              rga,rgb,rgc)
         write(AMA%ifig1,99999)noit,AMA%errrez,'B_insert',icha,AMA%glmin,AMA%glmax,nelem
         ichag = ichag + icha
         if(icha .le. 1) return
30    enddo
99999 format (i6,2x,1(e16.8),2x,a9,i5,2x,2e12.4,i8)
      return
    end subroutine C_B_INSERT_77


    subroutine REP_BOUND_77(ndim,melem,nelem,mpoin,npoin,  &
         maxdeg,x,y,lnd,iae,  &
         ra,icyc,rminerr,imin,rmaxrez,  &
         icha,nserr,wp,rga,rgb,rgc,lbn,itc,ibc,nbelm,mbelm,  &
         ipoint,nbp,xb,yb,ibb,ibpoin,ibp)
      dimension lnd(melem,3),iae(melem,3),x(mpoin),y(mpoin),  &
           ra(mpoin),wp(mpoin,ndim+1),lbn(mbelm,2),  &
           itc(mbelm),ibc(mbelm),  &
           icyc(mpoin,maxdeg),  &
           nserr(melem*3,2),rga(mpoin),rgb(mpoin),rgc(mpoin),  &
           xb(ipoint),yb(ipoint),ibp(mpoin,2),ibb(mpoin,3),ibpoin(0:nbp)
      integer lock(20)
      real xx(4),yy(4),rmax(3),jmax(3)
      real*8 x0,y0,x1,y1,x2,y2,det,rll,rl0,rl1,rl2,rlen0,  &
           rlen1,rlen2,rd,reps,det123,xx1,yy1,xx2,yy2
      icha = 0
      nelemold = nelem
      ipoc = 1
      do 10 iyi=1,nelemold
         i = ipoc
         do 20 j=1,3
            if(iae(i,j) .lt. 0) then
               j1 = mod(j,3) + 1
               j2 = mod(j1,3) + 1
               k1 = lnd(i,j)
               if(ibb(k1,1) .gt. 0) then 
                  ikk1 = ibb(k1,2)
                  if(ibb(k1,1) == ibpoin(ikk1-1)+1  .or.  &
                        ibb(k1,1) == ibpoin(ikk1)) then
!               NO moving of final or initial node of profiles
                     goto 21
                  endif
               endif
               k2 = lnd(i,j1)
               k0 = lnd(i,j2)
               x1 = x(k1)
               y1 = y(k1)
               x2 = x(k2)
               y2 = y(k2)
               x0 = x(k0)
               y0 = y(k0)


               if(ibb(k1,1) .gt. 0 .and. ibb(k2,1) .gt. 0 )then

                  if(ibb(k1,2)  .ne. ibb(k2,2) ) then
                     print *,'Troubles in insert.f'
                     print *,ibb(k1,2), ibb(k2,2) 
                  endif
                  irtk = ibb(k1,2)
                  
                  if(k0 == -1289) then
                     print *,x0,y0,k0
                     print *,x1,y1,k1
                     print *,x2,y2,k2
                     print *,'@@@@@'
                  endif

!     we seek the point in the field [xb(i),yb(i)] i=1,ipoint
!     ... rll length of edge k1,k2
                  rll = ((x(k1) - x(k2))*(x(k1) - x(k2)) +  &
                       (y(k1) - y(k2))*(y(k1) - y(k2)) )**0.5
                  rl0 = 1E+25*rll
                  do 30 ll1=ibpoin(irtk-1)+1,ibpoin(irtk)
                     rlen0 = (x0 -xb(ll1))*(x0 -xb(ll1)) +   &
                          (y0 -yb(ll1))*(y0 -yb(ll1))
                     if(rlen0 .lt. rl0) then
                        rl0 = rlen0
                        il0 = ll1
                     endif
30                enddo
                  il1 = ibb(k1,1)
                  il2 = ibb(k2,1)


                  if(il1 == il0 .or. il2 == il0) then
                     goto 20
                  endif


                  
!     to verify: il0 must be between il1 and il2

!                  if( ( il0 .gt. il1 .and. il0 .gt. il2) .or.
!     *                 ( il0 .lt. il1 .and. il0 .lt. il2)) goto 20

                  irtkdel = abs(ibpoin(irtk-1)+1 - ibpoin(irtk))
                  ill1 = il1 - il0
                  if(abs(ill1) .gt. irtkdel/2) then
                     if(il1 .gt. il0) then
                        ill1 = ill1 - irtkdel
                     else
                        ill1 = ill1 + irtkdel
                     endif
!                     print *,'@@@',il1,il0,il2,ill1,ill2
                  endif
                  ill2 = il2 - il0
                  if(abs(ill2) .gt. irtkdel/2) then
                     if(il2 .gt. il0) then
                        ill2 = ill2 - irtkdel
                     else
                        ill2 = ill2 + irtkdel
                     endif
!                     print *,'@.@',il1,il0,il2,ill1,ill2
                  endif

                  if( ill1 * ill2  .gt. 0 ) goto 20
!     il0 must be between il1 and il2

!     ... computation of the length of il0 from the edge k1,k2
                  x0 = xb(il0)
                  y0 = yb(il0)
                  det = x0*(y1-y2) + x1*(y2-y0) + x2*(y0-y1)
                  rd = det/rll

                  if(k0 == -1289) then
                     print *,'!!',rl0**0.5,rd,0.45*rd,icyc(k0,1)
                  endif

                  if(rl0**0.5 .lt. 0.95*rd .and. icyc(k0,1) .gt. 0) then
!                     print *,'!!!',k0,rl0**0.5,rd,rd1
!                     print *,x0,y0,k0,il0
!                     print *,x1,y1,k1,il1
!                     print *,x2,y2,k2,il2
!                     print *,xb(il0),yb(il0)
!                     print *

!     this is a candidate, we check the positivity
                     do 40 kk=1,icyc(k0,1)
                        if(icyc(k0,kk+1) .ne. k1) then
                           kk1 = mod(kk,icyc(k0,1)) + 1
                           xx1 = x(icyc(k0,kk+1))
                           yy1 = y(icyc(k0,kk+1))
                           xx2 = x(icyc(k0,kk1+1))
                           yy2 = y(icyc(k0,kk1+1))
                           reps = AMA%pos*( ((x0-xx2)**2 + (y0-yy2)**2) +  &
                                ((xx1-xx2)**2 + (yy1-yy2)**2) +  &
                                ((x0-xx1)**2 + (y0-yy1)**2) )
                           det123 = x0*(yy1-yy2) + xx1*(yy2-y0) +   &
                                xx2*(y0-yy1) 
!                           print *,'...',det123,reps
                           if( det123 .le. reps) then
!     violation of positivity, go to next i
                              goto 20
                           endif
                           call POS2TEST_77(x0,y0,xx1,yy1,xx2,yy2,itet)
                           if(itet == 1) then
!     violation of positivity, go to next i
                              goto 20
                           endif

                        endif
40                   enddo

!     we remove this triangle
                     if(icyc(k1,1) .gt. 0 .or. icyc(k2,1) .gt. 0) then
                        print *,'very divny in REP_BOUN'
                        stop
                     endif

!                     print *,'@@@'
!                     print *, x(k0),y(k0)
!                     print *,xb(il0),yb(il0)


                     x(k0) = xb(il0)
                     y(k0) = yb(il0)
                     ibb(k0,1) = il0
                     do ll=1,nbp
                        if(il0 .ge. ibpoin(ll-1) +1 .and.  &
                             il0 .le. ibpoin(ll) ) then
                           ibb(k0,2) = ll
                           goto 109
                        endif
                     enddo
 109                 continue
                     

                     i1 = iae(i,j1)
                     i2 = iae(i,j2)

                     if(i1 .lt. 0) then
                        print *,'ERROR1 in REP_BOUN'
                     else
                        jj1 = 0
                        do ll=1,3
                           if(iae(i1,ll) == i) jj1 = ll
                        enddo
                     endif
                     if(i2 .lt. 0) then
                        print *,'ERROR2 in REP_BOUN'
                     else
                        jj2 = 0
                        do ll=1,3
                           if(iae(i2,ll) == i) jj2 = ll
                        enddo
                     endif
                     if(jj1 == 0) then
                        print *,'ERROR3 in REP_BOUN'
                     else
                        iae(i1,jj1) = -2
                     endif
                     if(jj2 == 0) then
                        print *,'ERROR4 in REP_BOUN'
                     else
                        iae(i2,jj2) = -2
                     endif

                     do 101 l=1,nelem
                        do 111 k=1,3
                           if(iae(l,k) .gt. i) iae(l,k) = iae(l,k)-1
111                     enddo
101                  enddo
                     do 100 l=i,nelem-1
                        do 110 k=1,3
                           lnd(l,k) = lnd(l+1,k)
                           iae(l,k) = iae(l+1,k)
110                     enddo
100                  enddo

                     if(icyc(k1,2) == k2) then
                        do l=1,abs(icyc(k1,1))-1
                           icyc(k1,l+1) = icyc(k1,l+2)
                        enddo
                        icyc(k1,1) = icyc(k1,1) + 1
                     elseif(icyc(k1,abs(icyc(k1,1))+1) == k2)then
                        icyc(k1,1) = icyc(k1,1) + 1
                     else
                        print *,'ERROR k2 is not in the cykle of k1'
                        stop
                     endif
                     if(icyc(k2,2) == k1) then
                        do l=1,abs(icyc(k2,1))-1
                           icyc(k2,l+1) = icyc(k2,l+2)
                        enddo
                        icyc(k2,1) = icyc(k2,1) + 1
                     elseif(icyc(k2,abs(icyc(k2,1))+1) == k1)then
                        icyc(k2,1) = icyc(k2,1) + 1
                     else
                        print *,'ERROR k1 is not in the cykle of k2'
                        stop
                     endif

                     kl1 = 0
                     kl2 = 0
                     do l=1,icyc(k0,1)
                        lock(l) = icyc(k0,l+1)
                        if(lock(l) == k1) kl1 = l
                        if(lock(l) == k2) kl2 = l
                     enddo
                     if(kl1*kl2 == 0) then
                        print *,'k1 and k2 are not in cykle of k0'
                        stop
                     endif

                     if(abs(kl1-kl2) == 1) then
                        if(kl1 .lt. kl2) then
                           do l=1, icyc(k0,1) - kl1
                              icyc(k0,l+1) = lock(kl1+l)
                           enddo
                           do l=1,kl1
                              icyc(k0,icyc(k0,1)-kl1+l+1) = lock(l)
                           enddo
                        else
                           do l=1, icyc(k0,1) - kl2
                              icyc(k0,l+1) = lock(kl2+l)
                           enddo
                           do l=1,kl2
                              icyc(k0,icyc(k0,1)-kl2+l+1) = lock(l)
                           enddo
                        endif
                        icyc(k0,1) = -icyc(k0,1)
                     else
!     all is OK
                        icyc(k0,1) = -icyc(k0,1)
                     endif
                     
!     repair of boundary files
                     do ib=1,nbelm
                        if(itc(ib) .gt. i) itc(ib) = itc(ib) - 1
                     enddo

                     nbelmold = nbelm
                     do 600 ib=1,nbelmold
                        if(lbn(ib,1) == k1   &
                             .and. lbn(ib,2) == k2) then
                           do ib1=0,nbelm-ib-1
                              lbn(nbelm+1-ib1,1) = lbn(nbelm-ib1,1)
                              lbn(nbelm+1-ib1,2) = lbn(nbelm-ib1,2)
                              ibc(nbelm+1-ib1) = ibc(nbelm-ib1)
                              itc(nbelm+1-ib1) = itc(nbelm-ib1)
                           enddo
                           lbn(ib,2) = k0
                           lbn(ib+1,1) = k0
                           lbn(ib+1,2) = k2
                           ibc(ib+1) = ibc(ib)
                           if(i2.gt.i) then
                              itc(ib) = i2-1
                           else
                              itc(ib) = i2
                           endif
                           if(i1.gt.i) then
                              itc(ib+1) = i1-1
                           else
                              itc(ib+1) = i1
                           endif
                           nbelm = nbelm + 1
                           goto 610
                        endif
600                  enddo
                     print *,'boundary segment doesn''t found'
                     stop
610                  continue
                     nelem = nelem - 1
                     icha = icha + 1
                     
                     if(AMA%ifig .ge. 0 .and. mod(icha,5) == 0 ) then
                        call PLOT1_77(melem,nelem,mpoin,npoin,x,y,lnd)
                        AMA%ifig = AMA%ifig + 1
                     endif
                     
                  endif
               endif
            endif
21          continue
20       enddo
         if( i .gt. nelem) goto 2000
         ipoc = ipoc + 1
10    enddo
2000  continue

      return
    end subroutine REP_BOUND_77

    subroutine INSERT_77(ndim,melem,nelem,mpoin,npoin,  &
         maxdeg,x,y,lnd,iae,  &
         ra,icyc,rminerr,imin,rmaxrez,  &
         icha,nserr,wp,rga,rgb,rgc,lbn,itc,ibc,nbelm,mbelm,  &
         ipoint,nbp,xb,yb,ibb,ibpoin,ibp,ier)
      dimension lnd(melem,3),iae(melem,3),x(mpoin),y(mpoin),  &
           ra(mpoin),wp(mpoin,ndim+1),  &
           lbn(mbelm,2),itc(mbelm),ibc(mbelm),  &
           icyc(mpoin,maxdeg),  &
           nserr(melem*3,2),rga(mpoin),rgb(mpoin),rgc(mpoin),  &
           xb(ipoint),yb(ipoint),ibp(mpoin,2),ibb(mpoin,3),ibpoin(0:nbp)
      real xx(4),yy(4),rmax(3)
      integer jmax(3)
      
      ier = 0
      
      icha = 0
      ice = 0
      
      itest = -3
      
      
      rlmax2 = 5.33
      
      do 5 i=1,melem
         nserr(i,1) = 0
5     enddo
      
      nelemold = nelem
      do 10 ipoc = 1,nelemold
         i = ipoc
         do 20 j=1,3
            j1 = mod(j,3) +1
            ii1 = lnd(i,j)
            ii2 = lnd(i,j1)
            xi = x(ii1)
            yi = y(ii1)
            xi1 = x(ii2)
            yi1 = y(ii2)
            zi = wp(ii1,1)
            zi1 = wp(ii2,1)
            a = (rga(ii1) + rga(ii2) )/2
            b = (rgb(ii1) + rgb(ii2) )/2
            c = (rgc(ii1) + rgc(ii2) )/2
            rmax(j) = ( a*(xi-xi1)*(xi-xi1) + c*(yi-yi1)*(yi-yi1)  &
                 +2*b*(xi-xi1)*(yi-yi1))
            jmax(j) = j
20       enddo

         do 25 k=1,3
            do 26 l=1,2
               if(rmax(l) .lt. rmax(l+1) ) then
                  rmaxhelp = rmax(l)
                  rmax(l) = rmax(l+1)
                  rmax(l+1) = rmaxhelp
                  jmaxhelp = jmax(l)
                  jmax(l) = jmax(l+1)
                  jmax(l+1) = jmaxhelp
               endif
26          enddo
25       enddo
         
         !         acc = ACCUTE_77(x(lnd(i,1)),y(lnd(i,1)),x(lnd(i,2)),y(lnd(i,2)),
         !     *        x(lnd(i,3)),y(lnd(i,3)), 1.) 
         !         if( (iae(i,1) .lt. 0 .or. iae(i,2) .lt. 0 
         !     *        .or. iae(i,3) .lt. 0 ) .and. acc .gt. 1. ) then
         if(I == itest) then
            write(*,'(2e12.4,i5)') x(lnd(I,1)),y(lnd(I,1)),i
            write(*,'(3e12.4,i5)')   &
                 x(lnd(I,2)),y(lnd(I,2)),rmax(1),jmax(1)
            write(*,'(3e12.4,i5)')   &
                 x(lnd(I,3)),y(lnd(I,3)),rmax(2),jmax(2)
            write(*,'(3e12.4,i5)')   &
                 x(lnd(I,1)),y(lnd(I,1)),rmax(3),jmax(3)
            print *
            print *,'#',jmax(1),jmax(2),jmax(3)
            print *,'#',rmax(1),rmax(2),rmax(3),nserr(i,1),acc
            print *,'#------------------------'
         endif


         do 28 l=1,3
            !     checking the dimension of arrays (even for periodical boundary)
            if(npoin .ge. mpoin-2 .or.nelem .ge. melem-4 .or.  &
                 nbelm .ge. mbelm-2 ) then
               print *,'Dimension in insert full'
               print *,'nelem,melem=',nelem,melem
               print *,'npoin,mpoin=',npoin,mpoin
               print *,'nbelm,mbelm=',nbelm,mbelm
               ier = -4
               return
            endif
            
            if(rmax(l) .ge. rlmax2) then
               !     ... NEW ACCUTE
               !            if(rmax(l) .ge. rlmax2 .or. 
               !     *           (iae(i,jmax(l)) .lt. 0 
               !     *           .and. acc*rmax(l)/2 .ge. rlmax2)  ) then

               j = jmax(l)
               if(iae(i,j) .gt. 0) then
                  !     for non boundary sides
                  ii = iae(i,j)
                  if(nserr(i,1) == 0 .and. nserr(ii,1) == 0) then
                     j1 = mod(j,3)+1
                     j2 = mod(j1,3)+1
                     do 30 jjj=1,3
                        if(iae(ii,jjj) == i) jj = jjj
30                   enddo
                     jj1 = mod(jj,3)+1
                     jj2 = mod(jj1,3)+1
                     k1 = lnd(i,j)
                     k2 = lnd(i,j1)
                     k3 = lnd(i,j2)
                     k4 = lnd(ii,jj2)
                     if(k2.ne.lnd(ii,jj) .or. k1.ne.lnd(ii,jj1))then
                        print *,'ERRROR in INSERT'
                        print *,i,k1,k2,k3,k4
                     endif
                     
                     !     check if no validation of positivity
                     xx(1) = x(k1)
                     xx(2) = x(k4)
                     xx(3) = x(k2)
                     xx(4) = x(k3)
                     yy(1) = y(k1)
                     yy(2) = y(k4)
                     yy(3) = y(k2)
                     yy(4) = y(k3)
                     xx0 = (x(k1)+x(k2) )/2
                     yy0 = (y(k1)+y(k2) )/2

                     do 43 kl = 1,4
                        kl1 = mod(kl, 4) + 1
                        reps = AMA%pos*( (xx0-xx(kl))**2+(yy0-yy(kl))**2 +  &
                             (xx(kl)-xx(kl1))**2+(yy(kl)-yy(kl1))**2 +  &
                             (xx0-xx(kl1))**2 + (yy0-yy(kl1))**2 )
                        det = xx0*(yy(kl)-yy(kl1)) +   &
                             xx(kl)*(yy(kl1)-yy0) + xx(kl1)*(yy0-yy(kl))
                        if( det .le. reps) then
!     violation of positivity, go to next j
                           goto 28
                        endif
                        call POS1TEST_77(xx0,yy0,xx(kl),yy(kl),  &
                             xx(kl1),yy(kl1),itet)
                        if(itet == 1) then
!     violation of positivity, go to next i
                           goto 28
                        endif
43                   enddo
                     
                     if(iae(i,j1) .gt. 0) then
                        ia1 = iae(i,j1)
                        do 40 kk =1,3
                           if(iae(ia1,kk) == i) ja1 = kk
40                      enddo
                     else
                        ia1 = -1
                     endif
                     if(iae(i,j2) .gt. 0) then
                        ia2 = iae(i,j2)
                        do 50 kk =1,3
                           if(iae(ia2,kk) == i) ja2 = kk
50                      enddo
                     else
                        ia2 = -1
                     endif
                     if(iae(ii,jj1) .gt. 0) then
                        iia1 = iae(ii,jj1)
                        do 60 kk =1,3
                           if(iae(iia1,kk) == ii) jja1 = kk
60                      enddo
                     else
                        iia1 = -1
                     endif
                     if(iae(ii,jj2) .gt. 0) then
                        iia2 = iae(ii,jj2)
                        do 70 kk =1,3
                           if(iae(iia2,kk) == ii) jja2 = kk
70                      enddo
                     else
                        iia2 = -2
                     endif
                     
                     if(icha == -1) then
                        write(*,'(a2,8i5)') '@@',i,ii,k1,k2,k3,k4,  &
                             npoin,nelem
                     endif
                     x(npoin+1) = (x(k1)+x(k2))/2
                     y(npoin+1) = (y(k1)+y(k2))/2
                     ibb(npoin+1,1) = -1
                     ibb(npoin+1,2) = 0
                     ibb(npoin+1,3) = 0

                     
                     lnd(i,j) = npoin+1
                     iae(i,j2) = nelem+1
                     iae(i,j) = nelem+2
                     lnd(ii,jj) = npoin+1
                     iae(ii,jj2) = nelem+2
                     iae(ii,jj) = nelem+1
                     lnd(nelem+1,1) = npoin+1
                     lnd(nelem+1,2) = k3
                     lnd(nelem+1,3) = k1
                     iae(nelem+1,1) = i
                     iae(nelem+1,2) = ia2
                     iae(nelem+1,3) = ii
                     if(ia2 .gt. 0) then
                        iae(ia2,ja2) = nelem+1
                     else
                        !     corection of adjacent triangle
                        do ib=1,nbelm
                           if(itc(ib) == i) itc(ib) = nelem+1
                        enddo
                     endif
                     
                     lnd(nelem+2,1) = npoin+1
                     lnd(nelem+2,2) = k4
                     lnd(nelem+2,3) = k2
                     iae(nelem+2,1) = ii
                     iae(nelem+2,2) = iia2
                     iae(nelem+2,3) = i
                     if(iia2 .gt. 0) then
                        iae(iia2,jja2) = nelem+2
                     else
                        !     corection of adjacent triangle
                        do ib=1,nbelm
                           if(itc(ib) == ii) itc(ib) = nelem+2
                        enddo
                     endif
                     
                     wp(npoin+1,1) = ( wp(k1,1) + wp(k2,1) )/2
                     
                     rga(npoin+1) = (rga(k1) + rga(k2))/2
                     rgb(npoin+1) = (rgb(k1) + rgb(k2))/2
                     rgc(npoin+1) = (rgc(k1) + rgc(k2))/2
                     ibp(npoin+1,1) = 0
                     ibp(npoin+1,2) = 0
                     
                     npoin = npoin + 1 
                     nelem = nelem + 2
                     icha = icha + 1
                     
                     if(AMA%ifig .ge. 0 .and. mod(icha,5) == 0 ) then
                        call PLOT1_77(melem,nelem,mpoin,npoin,x,y,lnd)
                        AMA%ifig = AMA%ifig + 1
                     endif
                     
                     if(npoin .gt. mpoin .or. nelem .gt. melem) then
                        print *,'Error in dimension in insert'
                        print *,'nelem,melem=',nelem,melem
                        print *,'npoin,mpoin=',npoin,mpoin
                        stop
                     endif
                     
                     !     ... in each triangle only one division
                     nserr(i,1) = -1
                     nserr(ii,1) = -1
                     nserr(nelem,1) = -1
                     nserr(nelem-1,1) = -1
                     goto 10
                  endif
               else
                  !     for boundary sides
                  if(I == itest) print *,'$$',itest,nserr(i,1),  &
                       lnd(I,1),lnd(I,2),lnd(i,3)
                  if(nserr(i,1) == 0 )then
                     ipe = 0
999                  j1 = mod(j,3)+1
                     j2 = mod(j1,3)+1
                     k1 = lnd(i,j)
                     k2 = lnd(i,j1)
                     k3 = lnd(i,j2)
                     ia1 = iae(i,j1)
                     ia2 = iae(i,j2)
                     il0 = -1
                     if(ipe == 1) goto 128
!     for periodic boundary, the second point

!     we seek the poin inthe field [xb(i),yb(i)] i=1,ipoint
                     rll = ((x(k1) - x(k2))*(x(k1) - x(k2)) +  &
                          (y(k1) - y(k2))*(y(k1) - y(k2)) )
                     x0 = (x(k1) + x(k2))/2 
                     y0 = (y(k1) + y(k2))/2 


!                     print*,'@@@@',x0,y0
                     
                     if(ibb(k1,1) .gt. 0 .and. ibb(k2,1) .gt. 0   &
                          .and. ibb(k1,2) == ibb(k2,2) ) then
                        il1 = ibb(k1,1)
                        il2 = ibb(k2,1)
                        rl0 = 1E+25*(rll**0.5)
                        
                        !     test if between il1 and il2 is a point
                        idif = ibpoin(ibb(k1,2))-ibpoin(ibb(k1,2)-1)-1
                        !                        print *,'@@@@@',il1,il2,idif
                        if(abs(il1-il2) == 1 ) then
                           !                        .or. 
                           !     *                       abs(il1 - il2)  == idif ) then
                           print *,'We can not insert a new node,',  &
                                'there is few points on the profile'
                           print *,il1,il2,idif,ibpoin(ibb(k1,2)),  &
                                ibpoin(ibb(k1,2)-1)
                           print *,k1,x(k1),y(k1)
                           print *,k2,x(k2),y(k2)
                           goto 28
                        endif
                        if( il1 .gt. il2) then
!     0 node id between il1 and il2
                           il2new = il2 + ibpoin(ibb(k1,2))
                        else
                           il2new = il2
                        endif
                        do 213 ll1=il1,il2new
                           ll11 = ll1
                           if(ll11 .gt. ibpoin(ibb(k1,2)) )  &
                                ll11 = ll11 - ibpoin(ibb(k1,2))
                           rlen0 = (x0 -xb(ll11))*(x0 -xb(ll11)) +   &
                                (y0 -yb(ll11))*(y0 -yb(ll11))

                           if(rlen0 .lt. rl0) then
                              rl0 = rlen0
                              il0 = ll11
                           endif
 213                    enddo
                        if(rl0 .gt. 0.3*rll) then
                           print *,'very divnyin INSERT.F',k1,k2
                           print *,x(k1),y(k1)
                           print *,x(k2),y(k2)
                           print *
                           print *,x0,y0
                           print *
                           print *,xb(il0),yb(il0)
                           print *,xb(il1),yb(il1)
                           print *,xb(il2),yb(il2)
                           stop
                        endif
                        x0 = xb(il0)
                        y0 = yb(il0)
                     endif

                     if(i == itest) print *,'##',x0,y0,il0
                     if(i == itest) print *,'##',x(k1),y(k1),k1
                     if(i == itest) print *,'##',x(k2),y(k2),k2
                     if(i == itest) print *,'##',x(k3),y(k3),k3

                     reps = AMA%pos*( (x0-x(k1))**2+(y0-y(k1))**2 +  &
                          (x(k1)-x(k3))**2+(y(k1)-y(k3))**2 +  &
                          (x0-x(k3))**2 + (y0-y(k3))**2 )
                     det = x(k1)*(y0-y(k3)) + x0*(y(k3)-y(k1)) +   &
                          x(k3)*(y(k1)-y0)

!                     if(i .ne. itest) print *,'???',det,reps

                     if( det .le. reps ) then
!     violation of positivity, go to next j
                        goto 28
                     endif
                     call POS1TEST_77(x0,y0,x(k1),y(k1),x(k3),y(k3),itet)
                     if(i == itest) print *,'? ?',itet
                     if(itet == 1) then
!     violation of positivity, go to next i
                        goto 28
                     endif


                     reps = AMA%pos*( (x0-x(k3))**2+(y0-y(k3))**2 +  &
                          (x(k3)-x(k2))**2+(y(k3)-y(k2))**2 +  &
                          (x0-x(k2))**2 + (y0-y(k2))**2 )
                     det = x(k3)*(y0-y(k2)) + x0*(y(k2)-y(k3)) +   &
                          x(k2)*(y(k3)-y0)
!                     if(i == itest) print *,'???',det,reps
                     if( det .le. reps ) then
!     violation of positivity, go to next j
                        goto 28
                     endif
                     call POS1TEST_77(x0,y0,x(k2),y(k2),x(k3),y(k3),itet)
                     if(i == itest) print *,'? ?',itet,AMA%pos1
                     if(itet == 1) then
!     violation of positivity, go to next i
                        goto 28
                     endif

                     if(i == itest) write(*,'(a2,4e12.4)')  &
                             '##',x0,y0,det,reps
                     if(i == itest) write(*,'(a2,3i5)')'@@',  &
                          k1,ibp(k1,1),ibp(k1,2)
                     if(i == itest) write(*,'(a2,3i5)')'@@',  &
                          k2,ibp(k2,1),ibp(k2,2)
!     periodic boundary
                     if(ibp(k1,1) .gt. 0 .and. ibp(k2,1) .gt. 0 ) then
                        if((ibp(k1,2) .ne. ibp(k2,2)) .and.  &
                             (ibp(k1,2) .ne. 3 .and. ibp(k2,2) .ne. 3))  &
                             goto 28

                        if(ibp(k1,2) == 3) then
                           ibper = ibp(k2,2)
                        else
                           ibper = ibp(k1,2)
                        endif

                        if(ibper == 1) then
                           xperreal = AMA%xper(1,1)
                           yperreal = AMA%xper(1,2)
                        else
                           xperreal = AMA%xper(2,1)
                           yperreal = AMA%xper(2,2)
                        endif                  


                        do 1000 iel =1,nelem
                           do 1001 jel=1,3
                              if(iae(iel,jel) .lt. 0 .and.  &
                                   lnd(iel,jel) == ibp(k2,1))goto 1002
 1001                      enddo
 1000                   enddo
 1002                   continue
                        je1 = mod(jel,3)+1
                        je2 = mod(je1,3)+1
                        ke1 = lnd(iel,jel)
                        ke2 = lnd(iel,je1)
                        ke3 = lnd(iel,je2)
                        if(abs(x(k2) + xperreal - x(ke1) ) .lt.   &
                             1E-05 .and.  &
                             abs(y(k2)+yperreal-y(ke1) ) .lt.   &
                             1E-05 ) then
                           imov = 1
                        elseif(abs(x(k2)-xperreal-x(ke1)) .lt.   &
                                1E-05 .and.  &
                                abs(y(k2)-yperreal-y(ke1)).lt.   &
                                1E-05 ) then
                           imov = -1
                        else
                           print *,'BAD in insert in periodical points'
                           print *,i,k2,ke1
                           print *,x(k2),y(k2),xperreal
                           print *,x(ke1),y(ke1),yperreal
                           print *,abs(x(k2) + xperreal - x(ke1) ),  &
                                abs(y(k2) + yperreal - y(ke1) ),  &
                                abs(x(k2) - xperreal - x(ke1) ),  &
                                abs(y(k2) - yperreal - y(ke1) )
                           stop
                        endif
                        
                        xe0 = x0 + imov*xperreal
                        ye0 = y0 + imov*yperreal


                        reps = AMA%pos*( (xe0-x(ke1))**2+(ye0-y(ke1))**2 +  &
                             (x(ke1)-x(ke3))**2+(y(ke1)-y(ke3))**2 +  &
                             (xe0-x(ke3))**2 + (ye0-y(ke3))**2 )
                        det = x(ke1)*(ye0-y(ke3)) + xe0*(y(ke3)-y(ke1))  &
                             + x(ke3)*(y(ke1)-ye0)

                     if(i == itest) write(*,'(a2,4e12.4)')  &
                             '# ',xe0,ye0,det,reps

                        if( det .le. reps ) then
!     violation of positivity, go to next j
                           goto 28
                        endif
                        call POS1TEST_77(xe0,ye0,x(ke1),y(ke1),  &
                             x(ke3),y(ke3),itet)
                        if(itet == 1) then
                           !     violation of positivity, go to next i
                           goto 28
                        endif
                        
                        reps = AMA%pos*( (xe0-x(ke3))**2+(ye0-y(ke3))**2 +  &
                             (x(ke3)-x(ke2))**2+(y(ke3)-y(ke2))**2 +  &
                             (xe0-x(ke2))**2 + (ye0-y(ke2))**2 )
                        det = x(ke3)*(ye0-y(ke2)) + xe0*(y(ke2)-y(ke3))  &
                             + x(ke2)*(y(ke3)-ye0)
                        if( det .le. reps ) then
                           !     violation of positivity, go to next j
                           goto 28
                        endif
                        call POS1TEST_77(xe0,ye0,x(ke1),y(ke1),  &
                             x(ke2),y(ke2),itet)
                        if(itet == 1) then
                           !     violation of positivity, go to next i
                           goto 28
                        endif
                     endif
                     
128                  continue
                     
                     if(ia1 .gt. 0) then
                        ja1 = 0
                        do 125 ll=1,3
                           if(iae(ia1,ll) == i) then
                              ja1 = ll
                           endif
125                     enddo
                        if(ja1 == 0) then 
                           print *,'ERROR in INSERT_BOUNDARY-1'
                           stop
                        endif
                     endif
                     
                     if(ipe == 0) then
                        x(npoin+1) = x0
                        y(npoin+1) = y0
                        ibb(npoin+1,1) = il0
                        ibb(npoin+1,3) = 0
                     elseif(ipe == 1) then
                        x(npoin+1) = xe0
                        y(npoin+1) = ye0
                        ibb(npoin+1,1) = il0
                        ibb(npoin+1,3) = 0
                     else
                        print *,'bad value of ipe =',ipe
                     endif
                     
                     if(ibb(k1,2) == 0 .or. ibb(k2,2) == 0) then
                        ibb(npoin+1,2) = 0
                     elseif(ibb(k1,2) == ibb(k2,2)) then
                        ibb(npoin+1,2) = ibb(k1,2)
                     else
                        !                        print *,'error jkol1'
                     endif


                     wp(npoin+1,1) = (wp(k1,1) + wp(k2,1))/2
                     rga(npoin+1) = (rga(k1) + rga(k2))/2
                     rgb(npoin+1) = (rgb(k1) + rgb(k2))/2
                     rgc(npoin+1) = (rgc(k1) + rgc(k2))/2
                     
                     lnd(i,j1) = npoin+1
                     iae(i,j1) = nelem+1
                     
                     lnd(nelem+1,1) = npoin+1
                     lnd(nelem+1,2) = k2
                     lnd(nelem+1,3) = k3
                     iae(nelem+1,1) = -2
                     iae(nelem+1,2) = ia1
                     iae(nelem+1,3) = i

                     if(ia1 .gt. 0) iae(ia1,ja1) = nelem+1

!     the change in lbn, nbc,itc
                     ib1 = 0
                     do 265 ib=1,nbelm
                        if(lbn(ib,1) == k1 .and.lbn(ib,2) == k2)then
                           ib1 = ib
                           goto 266
                        endif
265                  enddo
266                  continue
                     if(ib1 == 0 ) then
                        print *,'ERROR in INSERT'
                        print *,'the boundary segment does not found'
                        stop
                     endif
                     
                     do 276 ii=1,nbelm-ib1
                        ib = nbelm +2 -ii
                        lbn(ib,1) = lbn(ib-1,1)
                        lbn(ib,2) = lbn(ib-1,2)
                        ibc(ib) = ibc(ib-1)
                        itc(ib) = itc(ib-1)
276                  enddo
                     

                     lbn(ib1,2) = npoin + 1
                     lbn(ib1+1,1) = npoin + 1
                     lbn(ib1+1,2) = k2
                     ibc(ib1+1) = ibc(ib1)
                     itc(ib1+1) = nelem+1
                     if(ipe == 0) ibp(npoin+1,1) = 0
                     if(ipe == 0) ibp(npoin+1,2) = 0


                     npoin = npoin + 1
                     nelem = nelem + 1
                     nbelm = nbelm + 1

                     if(npoin .gt. mpoin .or.nelem .gt. melem .or.  &
                          nbelm .gt. mbelm ) then
                        print *,'ERROR in dimension in insert'
                        print *,'nelem,melem=',nelem,melem
                        print *,'npoin,mpoin=',npoin,mpoin
                        print *,'nbelm,mbelm=',nbelm,mbelm
                        stop
                     endif
                     
                     nserr(i,1) = -1
                     nserr(nelem,1) = -1
                     if( ibp(k1,1) .gt. 0 .and. ibp(k2,1) .gt. 0 .and.   &
                          ipe == 0) then
                        jbak = j
                        i = iel
                        j = jel
                        ibp(npoin,1) = npoin + 1
                        ibp(npoin+1,1) = npoin
                        ibp(npoin,2) = ibper
                        ibp(npoin+1,2) = ibper
                        ipe = 1
                        goto 999
                     endif
                     if(ipe == 1) then
                        ipe = 0
                        i = ipoc
                        j = jbak
                     endif
                     icha = icha + 1
                     goto 10
                  endif
               endif
            endif
28       enddo
10    enddo
      !11 enddo

      return
      
    end subroutine INSERT_77



    subroutine INSERT_BOUND_77(ndim,melem,nelem,mpoin,npoin,  &
         maxdeg,x,y,lnd,iae,  &
         ra,icyc,rminerr,imin,rmaxrez,  &
         icha,nserr,wp,rga,rgb,rgc,lbn,itc,ibc,nbelm,mbelm,  &
         ipoint,nbp,xb,yb,ibb,ibpoin,ibp,ier)
      dimension lnd(melem,3),iae(melem,3),x(mpoin),y(mpoin),  &
           ra(mpoin),wp(mpoin,ndim+1),  &
           lbn(mbelm,2),itc(mbelm),ibc(mbelm),  &
           icyc(mpoin,maxdeg),  &
           nserr(melem*3,2),rga(mpoin),rgb(mpoin),rgc(mpoin),  &
           xb(ipoint),yb(ipoint),ibp(mpoin,2),ibb(mpoin,3),ibpoin(0:nbp)
      real xx(4),yy(4),rmax(3)
      integer jmax(3)

      ier = 0

      icha = 0
      ice = 0

      itest = -3


      rlmax2 = 5.33

      do 5 i=1,melem
         nserr(i,1) = 0
 5    enddo

      nelemold = nelem
      do 10 ipoc = 1,nelemold
         i = ipoc
         do 20 j=1,3
            j1 = mod(j,3) +1
            ii1 = lnd(i,j)
            ii2 = lnd(i,j1)
            xi = x(ii1)
            yi = y(ii1)
            xi1 = x(ii2)
            yi1 = y(ii2)
            zi = wp(ii1,1)
            zi1 = wp(ii2,1)
            a = (rga(ii1) + rga(ii2) )/2
            b = (rgb(ii1) + rgb(ii2) )/2
            c = (rgc(ii1) + rgc(ii2) )/2
            rmax(j) = ( a*(xi-xi1)*(xi-xi1) + c*(yi-yi1)*(yi-yi1)  &
                 +2*b*(xi-xi1)*(yi-yi1))
            jmax(j) = j
 20      enddo

         do 25 k=1,3
            do 26 l=1,2
               if(rmax(l) .lt. rmax(l+1) ) then
                  rmaxhelp = rmax(l)
                  rmax(l) = rmax(l+1)
                  rmax(l+1) = rmaxhelp
                  jmaxhelp = jmax(l)
                  jmax(l) = jmax(l+1)
                  jmax(l+1) = jmaxhelp
               endif
 26         enddo
 25      enddo

!         if( (iae(i,1) .lt. 0 .or. iae(i,2) .lt. 0 
!     *        .or. iae(i,3) .lt. 0 ) .and. acc .gt. 1. ) then
         if(I == itest) then
            write(*,'(2e12.4,i5)') x(lnd(I,1)),y(lnd(I,1)),i
            write(*,'(3e12.4,i5)')   &
                 x(lnd(I,2)),y(lnd(I,2)),rmax(1),jmax(1)
            write(*,'(3e12.4,i5)')   &
                 x(lnd(I,3)),y(lnd(I,3)),rmax(2),jmax(2)
            write(*,'(3e12.4,i5)')   &
                 x(lnd(I,1)),y(lnd(I,1)),rmax(3),jmax(3)
            print *
            print *,'#',jmax(1),jmax(2),jmax(3)
            print *,'#',rmax(1),rmax(2),rmax(3),nserr(i,1),acc
            print *,'#------------------------'
         endif


         do 28 l=1,3
!     checking the dimension of arrays (even for periodical boundary)
            if(npoin .ge. mpoin-2 .or.nelem .ge. melem-4 .or.  &
                 nbelm .ge. mbelm-2 ) then
               print *,'Dimension in insert_bound full'
               print *,'nelem,melem=',nelem,melem
               print *,'npoin,mpoin=',npoin,mpoin
               print *,'nbelm,mbelm=',nbelm,mbelm
               ier = -4
               return
            endif

            acc = ACCUTE_I_77(x(lnd(i,1)),y(lnd(i,1)),  &
                 x(lnd(i,2)),y(lnd(i,2)),x(lnd(i,3)),y(lnd(i,3)), l,1.) 

!            if(acc .gt. 1 .and. iae(i,l) .lt. 0) then
!            if(iae(i,l) .lt. 0) then
!               write(*,'(2e12.4,i5)') x(lnd(i,1)),y(lnd(i,1)),i
!               write(*,'(2e12.4,i5)') x(lnd(i,1)),y(lnd(i,1)),iae(i,l)
!               write(*,'(2e12.4)') x(lnd(i,1)),y(lnd(i,1))
!               write(*,'(3e12.4)') x(lnd(i,1)),y(lnd(i,1)),acc
!               write(*,'(x)')
!            endif
               
!            if(rmax(l) .ge. rlmax2) then
!     ... NEW ACCUTE
!     ... only non accute edges
            if(acc .gt. 10. ) then 

!               j = jmax(l)

               j = l
               if(iae(i,j) .gt. 0) then
!     for non boundary sides
!           no inserting

               else          
!     for boundary sides
                  if(I == itest) print *,'$$',itest,nserr(i,1),  &
                       lnd(I,1),lnd(I,2),lnd(i,3)
                  if(nserr(i,1) == 0 )then
                     ipe = 0
 999                 j1 = mod(j,3)+1
                     j2 = mod(j1,3)+1
                     k1 = lnd(i,j)
                     k2 = lnd(i,j1)
                     k3 = lnd(i,j2)
                     ia1 = iae(i,j1)
                     ia2 = iae(i,j2)
                     il0 = -1
                     if(ipe == 1) goto 128
!     for periodic boundary, the second point

!     we seek the poin inthe field [xb(i),yb(i)] i=1,ipoint
                     rll = ((x(k1) - x(k2))*(x(k1) - x(k2)) +  &
                          (y(k1) - y(k2))*(y(k1) - y(k2)) )
                     x0 = (x(k1) + x(k2))/2 
                     y0 = (y(k1) + y(k2))/2 

                     if(ibb(k1,1) .gt. 0 .and. ibb(k2,1) .gt. 0   &
                          .and. ibb(k1,2) == ibb(k2,2) ) then
                        il1 = ibb(k1,1)
                        il2 = ibb(k2,1)
                        rl0 = 1E+25*(rll**0.5)

!     test if between il1 and il2 is a point
                        idif = ibpoin(ibb(k1,2))-ibpoin(ibb(k1,2)-1)-1
!                        print *,'@@@@@',il1,il2,idif
                        if(abs(il1-il2) == 1 ) then
!                        .or. 
!     *                       abs(il1 - il2)  == idif ) then
                           print *,'We can not insert a new node,',  &
                                'there is few points on the profile'
                           print *,il1,il2,idif,ibpoin(ibb(k1,2)),  &
                                ibpoin(ibb(k1,2)-1)
                           print *,k1,x(k1),y(k1)
                           print *,k2,x(k2),y(k2)
                           goto 28
                        endif
                        if( il1 .gt. il2) then
!     0 node id between il1 and il2
                           il2new = il2 + ibpoin(ibb(k1,2))
                        else
                           il2new = il2
                        endif
                        do 213 ll1=il1,il2new
                           ll11 = ll1
                           if(ll11 .gt. ibpoin(ibb(k1,2)) )  &
                                ll11 = ll11 - ibpoin(ibb(k1,2))
                           rlen0 = (x0 -xb(ll11))*(x0 -xb(ll11)) +   &
                                (y0 -yb(ll11))*(y0 -yb(ll11))
                           if(rlen0 .lt. rl0) then
                              rl0 = rlen0
                              il0 = ll11
                           endif
 213                    enddo
                        if(rl0 .gt. 0.3*rll) then
                           print *,'very divnyin INSERT.F',k1,k2
                           print *,x(k1),y(k1)
                           print *,x(k2),y(k2)
                           print *
                           print *,x0,y0
                           print *
                           print *,xb(il0),yb(il0)
                           print *,xb(il1),yb(il1)
                           print *,xb(il2),yb(il2)
                           stop
                        endif
                        x0 = xb(il0)
                        y0 = yb(il0)
                     endif

!                     if(i .ne. itest) print *,'##',x0,y0,il0
!                     if(i .ne. itest) print *,'##',x(k1),y(k1),k1
!                     if(i .ne. itest) print *,'##',x(k2),y(k2),k2
!                     if(i .ne. itest) print *,'##',x(k3),y(k3),k3

                     reps = AMA%pos*( (x0-x(k1))**2+(y0-y(k1))**2 +  &
                          (x(k1)-x(k3))**2+(y(k1)-y(k3))**2 +  &
                          (x0-x(k3))**2 + (y0-y(k3))**2 )
                     det = x(k1)*(y0-y(k3)) + x0*(y(k3)-y(k1)) +   &
                          x(k3)*(y(k1)-y0)

                     if(i == itest) print *,'???',det,reps

                     if( det .le. reps ) then
!     violation of positivity, go to next j
                        goto 28
                     endif
                     call POS1TEST_77(x0,y0,x(k1),y(k1),x(k3),y(k3),itet)
                     if(i == itest) print *,'? ?',itet
                     if(itet == 1) then
!     violation of positivity, go to next i
                        goto 28
                     endif


                     reps = AMA%pos*( (x0-x(k3))**2+(y0-y(k3))**2 +  &
                          (x(k3)-x(k2))**2+(y(k3)-y(k2))**2 +  &
                          (x0-x(k2))**2 + (y0-y(k2))**2 )
                     det = x(k3)*(y0-y(k2)) + x0*(y(k2)-y(k3)) +   &
                          x(k2)*(y(k3)-y0)
                     if(i == itest) print *,'???',det,reps
                     if( det .le. reps ) then
!     violation of positivity, go to next j
                        goto 28
                     endif
                     call POS1TEST_77(x0,y0,x(k2),y(k2),x(k3),y(k3),itet)
                     if(i == itest) print *,'? ?',itet,AMA%pos1
                     if(itet == 1) then
!     violation of positivity, go to next i
                        goto 28
                     endif

                     if(i == itest) write(*,'(a2,4e12.4)')  &
                             '##',x0,y0,det,reps
                     if(i == itest) write(*,'(a2,3i5)')'@@',  &
                          k1,ibp(k1,1),ibp(k1,2)
                     if(i == itest) write(*,'(a2,3i5)')'@@',  &
                          k2,ibp(k2,1),ibp(k2,2)
!     periodic boundary
                     if(ibp(k1,1) .gt. 0 .and. ibp(k2,1) .gt. 0 ) then
                        if((ibp(k1,2) .ne. ibp(k2,2)) .and.  &
                             (ibp(k1,2) .ne. 3 .and. ibp(k2,2) .ne. 3))  &
                             goto 28

                        if(ibp(k1,2) == 3) then
                           ibper = ibp(k2,2)
                        else
                           ibper = ibp(k1,2)
                        endif

                        if(ibper == 1) then
                           xperreal = AMA%xper(1,1)
                           yperreal = AMA%xper(1,2)
                        else
                           xperreal = AMA%xper(2,1)
                           yperreal = AMA%xper(2,2)
                        endif                  


                        do 1000 iel =1,nelem
                           do 1001 jel=1,3
                              if(iae(iel,jel) .lt. 0 .and.  &
                                   lnd(iel,jel) == ibp(k2,1))goto 1002
 1001                      enddo
 1000                   enddo
 1002                   continue
                        je1 = mod(jel,3)+1
                        je2 = mod(je1,3)+1
                        ke1 = lnd(iel,jel)
                        ke2 = lnd(iel,je1)
                        ke3 = lnd(iel,je2)
                        if(abs(x(k2) + xperreal - x(ke1) ) .lt.   &
                             1E-05 .and.  &
                             abs(y(k2)+yperreal-y(ke1) ) .lt.   &
                             1E-05 ) then
                           imov = 1
                        elseif(abs(x(k2)-xperreal-x(ke1)) .lt.   &
                                1E-05 .and.  &
                                abs(y(k2)-yperreal-y(ke1)).lt.   &
                                1E-05 ) then
                           imov = -1
                        else
                           print *,'BAD in insert in periodical points'
                           print *,i,k2,ke1
                           print *,x(k2),y(k2),xperreal
                           print *,x(ke1),y(ke1),yperreal
                           print *,abs(x(k2) + xperreal - x(ke1) ),  &
                                abs(y(k2) + yperreal - y(ke1) ),  &
                                abs(x(k2) - xperreal - x(ke1) ),  &
                                abs(y(k2) - yperreal - y(ke1) )
                           stop
                        endif
                        
                        xe0 = x0 + imov*xperreal
                        ye0 = y0 + imov*yperreal


                        reps = AMA%pos*( (xe0-x(ke1))**2+(ye0-y(ke1))**2 +  &
                             (x(ke1)-x(ke3))**2+(y(ke1)-y(ke3))**2 +  &
                             (xe0-x(ke3))**2 + (ye0-y(ke3))**2 )
                        det = x(ke1)*(ye0-y(ke3)) + xe0*(y(ke3)-y(ke1))  &
                             + x(ke3)*(y(ke1)-ye0)

                     if(i == itest) write(*,'(a2,4e12.4)')  &
                             '# ',xe0,ye0,det,reps

                        if( det .le. reps ) then
!     violation of positivity, go to next j
                           goto 28
                        endif
                        call POS1TEST_77(xe0,ye0,x(ke1),y(ke1),  &
                             x(ke3),y(ke3),itet)
                        if(itet == 1) then
!     violation of positivity, go to next i
                           goto 28
                        endif

                        reps = AMA%pos*( (xe0-x(ke3))**2+(ye0-y(ke3))**2 +  &
                             (x(ke3)-x(ke2))**2+(y(ke3)-y(ke2))**2 +  &
                             (xe0-x(ke2))**2 + (ye0-y(ke2))**2 )
                        det = x(ke3)*(ye0-y(ke2)) + xe0*(y(ke2)-y(ke3))  &
                             + x(ke2)*(y(ke3)-ye0)
                        if( det .le. reps ) then
!     violation of positivity, go to next j
                           goto 28
                        endif
                        call POS1TEST_77(xe0,ye0,x(ke1),y(ke1),  &
                             x(ke2),y(ke2),itet)
                        if(itet == 1) then
!     violation of positivity, go to next i
                           goto 28
                        endif
                     endif

 128                 continue

                     if(ia1 .gt. 0) then
                        ja1 = 0
                        do 125 ll=1,3
                           if(iae(ia1,ll) == i) then
                              ja1 = ll
                           endif
 125                    enddo
                        if(ja1 == 0) then 
                           print *,'ERROR in INSERT_BOUNDARY-1'
                           stop
                        endif
                     endif

                     if(ipe == 0) then
                        x(npoin+1) = x0
                        y(npoin+1) = y0
                        ibb(npoin+1,1) = il0
                        ibb(npoin+1,3) = 0
                     elseif(ipe == 1) then
                        x(npoin+1) = xe0
                        y(npoin+1) = ye0
                        ibb(npoin+1,1) = il0
                        ibb(npoin+1,3) = 0
                     else
                        print *,'bad value of ipe =',ipe
                     endif

                     if(ibb(k1,2) == 0 .or. ibb(k2,2) == 0) then
                        ibb(npoin+1,2) = 0
                     elseif(ibb(k1,2) == ibb(k2,2)) then
                        ibb(npoin+1,2) = ibb(k1,2)
                     else
!                        print *,'error jkol1'
                     endif


                     wp(npoin+1,1) = (wp(k1,1) + wp(k2,1))/2
                     rga(npoin+1) = (rga(k1) + rga(k2))/2
                     rgb(npoin+1) = (rgb(k1) + rgb(k2))/2
                     rgc(npoin+1) = (rgc(k1) + rgc(k2))/2
                     
                     lnd(i,j1) = npoin+1
                     iae(i,j1) = nelem+1

                     lnd(nelem+1,1) = npoin+1
                     lnd(nelem+1,2) = k2
                     lnd(nelem+1,3) = k3
                     iae(nelem+1,1) = -2
                     iae(nelem+1,2) = ia1
                     iae(nelem+1,3) = i

                     if(ia1 .gt. 0) iae(ia1,ja1) = nelem+1

!     the change in lbn, nbc,itc
                     ib1 = 0
                     do 265 ib=1,nbelm
                        if(lbn(ib,1) == k1 .and.lbn(ib,2) == k2)then
                           ib1 = ib
                           goto 266
                        endif
 265                 enddo
 266                 continue
                     if(ib1 == 0 ) then
                        print *,'ERROR in INSERT'
                        print *,'the boundary segment does not found'
                        stop
                     endif
                     
                     do 276 ii=1,nbelm-ib1
                        ib = nbelm +2 -ii
                        lbn(ib,1) = lbn(ib-1,1)
                        lbn(ib,2) = lbn(ib-1,2)
                        ibc(ib) = ibc(ib-1)
                        itc(ib) = itc(ib-1)
 276                 enddo
                     

                     lbn(ib1,2) = npoin + 1
                     lbn(ib1+1,1) = npoin + 1
                     lbn(ib1+1,2) = k2
                     ibc(ib1+1) = ibc(ib1)
                     itc(ib1+1) = nelem+1
                     if(ipe == 0) ibp(npoin+1,1) = 0
                     if(ipe == 0) ibp(npoin+1,2) = 0


                     npoin = npoin + 1
                     nelem = nelem + 1
                     nbelm = nbelm + 1

                     if(npoin .gt. mpoin .or.nelem .gt. melem .or.  &
                          nbelm .gt. mbelm ) then
                        print *,'ERROR in dimension in insert'
                        print *,'nelem,melem=',nelem,melem
                        print *,'npoin,mpoin=',npoin,mpoin
                        print *,'nbelm,mbelm=',nbelm,mbelm
                        stop
                     endif
                     
                     nserr(i,1) = -1
                     nserr(nelem,1) = -1
                     if( ibp(k1,1) .gt. 0 .and. ibp(k2,1) .gt. 0 .and.   &
                          ipe == 0) then
                        jbak = j
                        i = iel
                        j = jel
                        ibp(npoin,1) = npoin + 1
                        ibp(npoin+1,1) = npoin
                        ibp(npoin,2) = ibper
                        ibp(npoin+1,2) = ibper
                        ipe = 1
                        goto 999
                     endif 
                     if(ipe == 1) then
                        ipe = 0
                        i = ipoc
                        j = jbak
                     endif
                     icha = icha + 1
                     goto 10
                  endif
               endif
            endif
28       enddo
10    enddo

      return
    end subroutine INSERT_BOUND_77

      subroutine METRIX_77(ndim,melem,nelem,mpoin,npoin,nbelm,mbelm,  &
           maxdeg,x,y,lnd,iae,w,wp,supp,icyc,ra,surface,ibp,ifv,  &
           ipoint,nbp,ibb,ibpoin,xb,yb,tria, lbn, ibc)
      dimension lnd(melem,3),iae(melem,3),x(mpoin),y(mpoin),  &
           ra(melem*3),wp(mpoin,ndim+1),supp(mpoin),w(melem,ndim+1),  &
           icyc(mpoin,maxdeg),ibp(mpoin,2),ibpoin(0:nbp),  &
           ibb(mpoin,3),xb(ipoint),yb(ipoint),tria(melem),  &
           lbn(mbelm,2),ibc(mbelm)

      real*8 x0,y0,z0,x1,y1,z1,x2,y2,z2,deta,detb,detc,det0

      rkappa = 1.4

      surface = 0.
      do 2 i=1,nelem
         x1 = x(lnd(i,1))
         y1 = y(lnd(i,1))
         x2 = x(lnd(i,2))
         y2 = y(lnd(i,2))
         x3 = x(lnd(i,3))
         y3 = y(lnd(i,3))
         rmeas = (x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) )/2
         tria(i) = rmeas
         surface = surface + rmeas
         !write(21,*) (x1+x2+x3)/3, (y1+y2+y3)/3, w(i,1)
 2    enddo

      if(ifv == 0) then
         do i = 1,npoin
            do k = 1,ndim
               wp(i,k+1) = w(i,k)
            enddo
         enddo
      else
!     recomputation on the nodes and smoothing(for ismoothing > 1)
         if(ifv == 0) then
            ismoothing = 2
         else
            ismoothing = 1
         endif

         do 12 is=1,ismoothing
            do 10 i=1,npoin
               supp(i) = 0.
               do 15 j=1,ndim+1
                  wp(i,j) = 0.
 15            enddo
 10         enddo
            do 20 i=1,nelem
               !xc = (x(lnd(i,1))+x(lnd(i,2))+x(lnd(i,3)))/3
               !yc = (y(lnd(i,1))+y(lnd(i,2))+y(lnd(i,3)))/3
               do 30 j=1,3
                  !rl = ((xc-x(lnd(i,j)))**2 + (yc-y(lnd(i,j)))**2)**0.5
                  do 40 k=1,ndim
                     wp(lnd(i,j),k+1) = wp(lnd(i,j),k+1) +w(i,k)*tria(i)
 40               enddo
                  supp(lnd(i,j)) = supp(lnd(i,j)) + tria(i)
 30            enddo
 20         enddo
            do 45 i=1,npoin
               do 60 k=2,ndim+1
                  wp(i,k) = wp(i,k)/supp(i)
 60            enddo
 45         enddo



            do 48 i=1,nelem
               do 49 k=1,ndim
                  w(i,k) = (wp(lnd(i,1),k+1) + wp(lnd(i,2),k+1) +   &
                       wp(lnd(i,3),k+1) )/3
 49            enddo
 48         enddo
 12      enddo
      endif


!     correction for the Navier - Stokes
      do i=1,npoin

         !write(22,*) x(i), y(i), wp(i,2)

         if(ibp(i,1) == -1 .and. ( AMA%ityp == 5)) then
           wp(i,3) = wp(i,3)*0.05
           wp(i,4) = wp(i,4)*0.05
         endif
!     computation already done, we can forget this information
!     for simplicity
         if(ibp(i,1) == -1) ibp(i,1) = 0
      enddo


!     HERE IS POSSIBLE TO CHANGE THE USED QUANTITY FOR HESSIAN METRIXES
!     for public
!      do 50 i=1,npoin
!         if(AMA%ityp == 0 ) then
!     the uniform triangulation
!            wp(i,1) = 1.0
!         else
!            wp(i,1) = wp(i,AMA%ityp+1)
!         endif
! 50   enddo

      
      do 50 i=1,npoin
!         write(99,*) x(i),y(i), wp(i,2),wp(i,3),wp(i,4), wp(i,5)


         if(AMA%ityp == 0 .or. AMA%ityp == 3) then
!     the uniform triangulation
            wp(i,1) = 1.0
         elseif(AMA%ityp == -1) then 
!     the "exact mesh"
            xc = x(i)
            yc = y(i)
!c2            wp(i,1) = (xc*xc + yc*yc)/2.
!c3            wp(i,1) = (100*xc*xc + yc*yc)/2.
!c4
            rc = (xc*xc+yc*yc)
!c4
            wp(i,1) = 10*rc*exp(-10*(rc**0.5-1)**2)

         elseif(AMA%ityp == 1 .or. AMA%ityp == 4) then
!     the testing quantity is the density
            wp(i,1) = wp(i,2)
         elseif(AMA%ityp == 2 ) then
!     the testing quantity is the velocity
            wp(i,1) = ((wp(i,3)/wp(i,2))**2+(wp(i,4)/wp(i,2))**2 )**0.5
         elseif(AMA%ityp == 5 .or. AMA%ityp == 6) then
!     the testing quantity is the Mach number
            ro = wp(i,2)
            if(ro .le. 0)then
               print *,'density zero on element ',i,'=',ro
               stop
            endif
            u  = wp(i,3)/ro                                            
            v  = wp(i,4)/ro                                             
            p  = (rkappa-1.)*(wp(i,5)-0.5*ro*(u*u+v*v))                 
            if( p .le. 0. ) then                                      
               print *, ' Pressure <= 0.0 sur l''element ', i       
!               stop
               p = 0.001
            endif                                                      
            wp(i,1)=sqrt((u*u+v*v)*ro/(rkappa*p))               

         elseif(AMA%ityp == -10)then
            if(x(i) < 0. ) then
               ri = (x(i)*x(i) + y(i)*y(i))**0.5
            elseif(x(i) >= 0. .and. x(i) .le. 1. ) then
               if(abs(y(i)) < 1.0) then
                  ri = 1E+10
                  do j=1,nbelm
                     do l=1,2
                        xi = x(lbn(j,l)) - x(i)
                        yi = y(lbn(j,l)) - y(i)
                        ri = min(ri, (xi*xi + yi*yi)**0.5 )
                     enddo
                  enddo
               else
                  ri = abs(y(i))
               endif
            else
               ri = (y(i)* y(i) + 0.025*x(i)*x(i))**0.5
            endif
            !ri = min( (x(i) -0.25)**2, (x(i) -0.5)**2, (x(i) -0.75)**2)  &
            !     + y(i)**2 
            !wp(i,1) = 1.0 + 1E+8*exp(-2000*ri)
            !wp(i,1) = 100/(1+100 * ri**2)
            !wp(i,1) = 2.5*exp(-1 * max(0., ri-0.5 ) ) ! UA*
            wp(i,1) = 45*exp(-2.5 * max(0., ri-0.0 ) )   ! UB*
            if(x(i)*x(i) + y(i)*y(i) < 1E-2) wp(i,1) = wp(i,1)*5
            !write(99,*) x(i), y(i), ri
         else
            print *,'bad number of ityp, ityp = ',AMA%ityp
         endif

!         write(21,'(i5,7e12.4)') 
!     *        i,x(i),y(i),wp(i,1),wp(i,2),wp(i,3),wp(i,4),wp(i,5)
         
 50   enddo


!     boundary layers
      if(AMA%ityp == 3) then
         !print*,'###',wp(lbn(100,1), 1),wp(lbn(100,2), 1)
         do i=1,nbelm
            if(ibc(i) == 3 .or. ibc(i) == 4.) then
               wp(lbn(i,1), 1) = wp(lbn(i,1), 1)*0.9
!               wp(lbn(i,2), 1) = wp(lbn(i,2), 1)*0.8
!               print*,'         ',i,ibc(i)
            endif
         enddo


         !print*,'###',wp(lbn(100,1), 1),wp(lbn(100,2), 1)

      endif


!     periodic problems
      do 63 i=1,npoin
         do 64 k=1,ndim+1
            if(ibp(i,1) .gt. 0) then
               wp(i,k) = (wp(i,k) + wp(ibp(i,1),k))/2
               wp(ibp(i,1),k) = wp(i,k)
            endif
 64      enddo
 63   enddo

!     now only a local array
      do 65 i=1,npoin
         ra(i) = 0.
         supp(i) = 0.
 65   enddo


      return
    end subroutine METRIX_77

    subroutine ERROR_77(ndim,melem,nelem,mpoin,npoin,  &
         maxdeg,dx,dy,area,  &
         tria,noit,f1,f2,  &
         x,y,lnd,iae,w,wp,icyc,ra,  &
         rga,rgb,rgc,surface, ibp)

      dimension lnd(melem,3),iae(melem,3),x(mpoin),y(mpoin),  &
           ra(melem*3),wp(mpoin,5),dx(melem),dy(melem),  &
           icyc(mpoin,maxdeg),rga(mpoin),rgb(mpoin),rgc(mpoin),  &
           area(mpoin),tria(melem),  &
           ibp(mpoin,1),w(melem,ndim)

      dimension  f1(melem,4),f2(melem,4),  &
           nasobx(3),nasoby(3),sum(3),rmaxi(4),  &
           aa(100,3),er(3),rlong(3)
!     *     aa(ndim,3),er(3),rlong(3)

!      real*8  x1,y1,x2,y2,x3,y3,d1,d2,d3,e1,e2,e3
      real*8 a,b,c,disc,rlam1,rlam2,xm11,xm12,xm21,xm22,  &
           ym12,ym11,ym21,ym22,  &
           t11,t12,t21,t22,z11,z12,z21,z22,rdet

      rkappa = 1.4

      ielk = -3638

      ihel = 231
      open (ihel,file='hel',status='unknown')
      
      do 5 k=1,nelem
         pom=(rkappa-1)*(w(k,4)-0.5/w(k,1)*(w(k,2)*w(k,2)+  &
              w(k,3)*w(k,3)))
         pom1=(rkappa-1)*0.5*(w(k,2)*w(k,2)+w(k,3)*w(k,3))/w(k,1)
         f1(k,1)=w(k,2)
         f1(k,2)=w(k,2)*w(k,2)/w(k,1)+pom
         f1(k,3)=w(k,2)*w(k,3)/w(k,1)
         f1(k,4)=w(k,2)/w(k,1)*(rkappa*w(k,4)-pom1)
       
         f2(k,1)=w(k,3)
         f2(k,2)=f1(k,3)
         f2(k,3)=w(k,3)*w(k,3)/w(k,1)+pom
         f2(k,4)=w(k,3)/w(k,1)*(rkappa*w(k,4)-pom1)
 5    enddo

      rmin = 10.
      rmax = 0.

!     the area of ideal triangle
      rnorm = 1.29903810567 *AMA%numel/surface


! vypocet souradnic vrcholu trojuhelnika a stedu jeho hran
      do 10 k=1,nelem
         x1=x(lnd(k,1)) 
         y1=y(lnd(k,1)) 
         x2=x(lnd(k,2)) 
         y2=y(lnd(k,2)) 
         x3=x(lnd(k,3)) 
         y3=y(lnd(k,3)) 
         rlong(1) = ((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) )**0.5
         rlong(2) = ((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3) )**0.5
         rlong(3) = ((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) )**0.5

  
         x12=0.5*(x1+x2) 
         y12=0.5*(y1+y2) 
         x23=0.5*(x3+x2) 
         y23=0.5*(y3+y2) 
         x13=0.5*(x1+x3) 
         y13=0.5*(y1+y3) 
  
! poradova cisla sousedu
         l1=iae(k,1)
         l2=iae(k,2)
         l3=iae(k,3)

!      ...    turbine cascade
!  vypocet teziste trojuhelnika a zjisteni,zda neni na 
! periodicke hranici (je-li, pak nasob=1 nebo -1 )
         rperiod = 0.05511679
!       ***
         xik = 1./3*( x1 + x2 + x3 )  
         yik = 1./3*( y1 + y2 + y3 )  

         do 15 j=1,3
            nasobx(j)=0
            nasoby(j)=0
            kj = iae(k,j)
            if(kj.gt.0)then         
               xckj = 1./3*( x(lnd(kj,1)) + x(lnd(kj,2))   &
                    + x(lnd(kj,3)) ) 
               yckj = 1./3*( y(lnd(kj,1)) + y(lnd(kj,2))   &
                    + y(lnd(kj,3)) )

               if((yckj-yik) .gt. AMA%xper(1,2)/4)nasoby(j)=-1           
               if((yik-yckj) .gt. AMA%xper(1,2)/4)nasoby(j)=1 
               if((xckj-xik) .gt. AMA%xper(1,1)/4)nasobx(j)=-1           
               if((xik-xckj) .gt. AMA%xper(1,1)/4)nasobx(j)=1 
            endif
 15      enddo

!         print *, 'za vypoctem tezist  '     
! postupne u kazdeho souseda vytvoreni sestiuhelniku a ten se pak 
! posle do sixtriang, kde se vypocte norma fi

         if(l1.gt.0)then 
            do 2 j=1,3
               if(iae(l1,j)==k) lj=j-1
 2          enddo
            if(lj == 0) lj=3        
            xx12=0.5*(x1+x(lnd(l1,lj))+nasobx(1)*AMA%xper(1,1))
            yy12=0.5*(y1+y(lnd(l1,lj))+nasoby(1)*AMA%xper(1,2))
            xxx12=0.5*(x2+x(lnd(l1,lj))+nasobx(1)*AMA%xper(1,1))
            yyy12=0.5*(y2+y(lnd(l1,lj))+nasoby(1)*AMA%xper(1,2))

            call SIXTRIANG_77(x12,y12,x1,y1,xx12,yy12,xxx12,yyy12,x2,y2,    &
                 x23,y23,x13,y13,sum1) 
         else
            sum1=1.
         endif 
         
         if(l2.gt.0)then 
            do 3 j=1,3
               if(iae(l2,j)==k) lj=j-1
 3          enddo
            if(lj == 0) lj=3        
            xx23=0.5*(x2+x(lnd(l2,lj))+nasobx(2)*AMA%xper(1,1))
            yy23=0.5*(y2+y(lnd(l2,lj))+nasoby(2)*AMA%xper(1,2))
            xxx23=0.5*(x3+x(lnd(l2,lj))+nasobx(2)*AMA%xper(1,1))
            yyy23=0.5*(y3+y(lnd(l2,lj))+nasoby(2)*AMA%xper(1,2))

            call SIXTRIANG_77(x23,y23,x12,y12,x2,y2,xx23,yy23,xxx23,  &
                 yyy23,x3,y3,x13,y13,sum2) 
         else
            sum2=1.
         endif
         
         if(l3.gt.0)then 
            do 4 j=1,3
               if(iae(l3,j)==k) lj=j-1
 4          enddo
            if(lj == 0) lj=3        
            xx13=0.5*(x3+x(lnd(l3,lj))+nasobx(3)*AMA%xper(1,1))
            yy13=0.5*(y3+y(lnd(l3,lj))+nasoby(3)*AMA%xper(1,2))
            xxx13=0.5*(x1+x(lnd(l3,lj))+nasoby(3)*AMA%xper(1,1))
            yyy13=0.5*(y1+y(lnd(l3,lj))+nasoby(3)*AMA%xper(1,2))

            call SIXTRIANG_77(x13,y13,x1,y1,x12,y12,x23,y23,x3,y3,    &
                 xx13,yy13,xxx13,yyy13,sum3) 
         else
            sum3=1.
         endif

!         open(42,status='unknown' , file='zkousky')
!         if(sum1.gt.0) write(42,'(i5,3e16.5)') k,sum1,sum2,sum3

         sum(1)=sqrt(sum1)
         sum(2)=sqrt(sum2)
         sum(3)=sqrt(sum3)

! vypocet forem a , jejich abs hodnota a vydeleni normou fi pro tri
! sestiuheniky vztahujici se k trojuhelniku k 
         do 20 l=1,4
            do 30 j=1,3
               if(iae(k,j).gt.0)then
                  j1 = mod(j,3) + 1
                  rnu1kj = y(lnd(k,j1)) - y(lnd(k,j))
                  rnu2kj = x(lnd(k,j)) - x(lnd(k,j1))
                  kj=iae(k,j)
                  aa(l,j)=0.5*(rnu1kj*(f1(k,l)-f1(kj,l))+  &
                       rnu2kj*(f2(k,l)-f2(kj,l)))
                  aa(l,j)=abs(aa(l,j))
                  aa(l,j)=aa(l,j)/sum(j)
               else
                  aa(l,j)=0.
               endif 
               if(k == ielk) then
                  write(*,'(4i5)') k,l,j,iae(k,j)
                  write(*,'(4e12.4)')   &
                       aa(l,j),sum(j),0.5*(rnu1kj*(f1(k,l)-f1(kj,l))),  &
                       0.5*(rnu2kj*(f2(k,l)-f2(kj,l))),  &
                       0.5*(rnu1kj*(f1(k,l)-f1(kj,l))+  &
                       rnu2kj*(f2(k,l)-f2(kj,l)))
                  print *,'_________________' 
               endif
 30         enddo

            do 35 j=1,3
               if( aa(l,j) == 0) then
                  aa(l,j) = (aa(l,1) + aa(l,2) + aa(l,3))/2
               endif
               rmin = min (rmin,aa(1,j))
               rmax = max (rmax,aa(1,j))
 35         enddo
!            if( l== 1) 
!     *           write(*,'(2i5,3e12.4)') k,j,aa(1,1),aa(1,2),aa(1,3)


!            maximum ze tri vyrazu prave vypocteni coz je eta(l)
!            if( aa(l,1) .gt. aa(l,2))then
!               rmaxi(l)=aa(l,1)
!            else
!               rmaxi(l)=aa(l,2)
!            endif
!            
!            if(aa(l,3) .gt. rmaxi(l)) rmaxi(l)= aa(l,3)
 20      enddo   

!         eh1 = (aa(1,1)**2+aa(2,1)**2+aa(3,1)**2+aa(4,1)**2)**0.5
!         eh2 = (aa(1,2)**2+aa(2,2)**2+aa(3,2)**2+aa(4,2)**2)**0.5
!         eh3 = (aa(1,3)**2+aa(2,3)**2+aa(3,3)**2+aa(4,3)**2)**0.5

         eh1 = aa(1,1)
         eh2 = aa(1,2)
         eh3 = aa(1,3)

!         if(max(eh1,eh2,eh3) .gt. 0.002) then
         if(k == ielk) then
            write(ihel,'(2e14.6)') x1,y1
            write(ihel,'(2e14.6)') x2,y2
            write(ihel,'(2e14.6)') x3,y3
            write(ihel,'(2e14.6)') x1,y1
            write(ihel,'(x)')
         endif

         do 50 j=1,3
!            er(j) = + 1.E+03*
!     *           (aa(1,j)**2+aa(2,j)**2+aa(3,j)**2+aa(4,j)**2)**0.5
            er(j) = + 1.E+03 * aa(1,j)
 50      enddo

         e1 = er(1)
         e2 = er(2)
         e3 = er(3)


         d1 = (y2*e1*y1-e1*y3*y1-e1*y2*y3+e1*y3**2+y1**2*e2-e2*y3*y1-  &
              y1*e2*y2+e2*y2*y3-y2*y1*e3+e3*y3*y1+e3*y2**2-e3*y2*y3)/  &
              (y2**2*x1**2-2*y2*y3*x1**2+y3**2*x1**2-2*y2**2*x1*x3+  &
              2*y1*y2*x1*x3-2*x2*y2*y1*x1+2*x1*x2*y2*y3+2*y2*y3*x1*x3-  &
              2*x1*x2*y3**2+2*x1*x2*y3*y1-2*x1*y1*x3*y3+y2**2*x3**2+  &
              2*y1*y2*x2*x3-2*x2*y2*y3*x3-2*y1*y2*x3**2-2*x2**2*y3*y1-  &
              2*x2*x3*y1**2+x3**2*y1**2+x2**2*y1**2+x2**2*y3**2+  &
              2*x2*x3*y3*y1)

         d2 = -(e1*x1*y2-e1*x1*y3-e1*x3*y2+2*e1*x3*y3+e1*x2*y1-e1*y1*x3  &
              -e1*x2*y3-e2*x1*y2+2*e2*x1*y1-e2*x1*y3+e2*x3*y2-e2*x2*y1+  &
              e2*x2*y3-e2*y1*x3-e3*x1*y2+e3*x1*y3+2*e3*x2*y2-e3*x3*y2-  &
              e3*x2*y1-e3*x2*y3+e3*y1*x3)/  &
              (y2**2*x1**2-2*y2*y3*x1**2+y3**2*x1**2-2*y2**2*x1*x3+  &
              2*y1*y2*x1*x3-2*x2*y2*y1*x1+2*x1*x2*y2*y3+2*y2*y3*x1*x3-  &
              2*x1*x2*y3**2+2*x1*x2*y3*y1-2*x1*y1*x3*y3+y2**2*x3**2+  &
              2*y1*y2*x2*x3-2*x2*y2*y3*x3-2*y1*y2*x3**2-2*x2**2*y3*y1-  &
              2*x2*x3*y1**2+x3**2*y1**2+x2**2*y1**2+  &
              x2**2*y3**2+2*x2*x3*y3*y1)/2

         d3 = (e1*x2*x1-e1*x1*x3-e1*x2*x3+e1*x3**2+e2*x1**2-e2*x1*x2-  &
              e2*x1*x3+e2*x2*x3-x1*x2*e3+e3*x1*x3+e3*x2**2-e3*x2*x3)/  &
              (y2**2*x1**2-2*y2*y3*x1**2+y3**2*x1**2-2*y2**2*x1*x3+  &
              2*y1*y2*x1*x3-2*x2*y2*y1*x1+2*x1*x2*y2*y3+2*y2*y3*x1*x3-  &
              2*x1*x2*y3**2+2*x1*x2*y3*y1-2*x1*y1*x3*y3+y2**2*x3**2+  &
              2*y1*y2*x2*x3-2*x2*y2*y3*x3-2*y1*y2*x3**2-2*x2**2*y3*y1-  &
              2*x2*x3*y1**2+x3**2*y1**2+x2**2*y1**2+x2**2*y3**2+  &
              2*x2*x3*y3*y1)

!      print *,k,d1,d2,d3,d1*d3 - d2*d2

!         if(max(eh1,eh2,eh3) .gt. 0.003) then
!            print *,'BAD metric in ERROR'
         if(k == ielk) then
            write(*,'(2i5,4e12.4)') k,j,d1,d2,d3,d1*d3-d2*d2
            print *,x1,y1,rlong(1)
            print *,x2,y2,rlong(2)
            print *,x3,y3,rlong(3)
            print *,eh1,eh2,eh3
            print *,e1,e2,e3
            print *,d1,d2,d3
            print *,d1*(x2-x1)**2+ 2*d2*(x2-x1)*(y2-y1)+ d3*(y2-y1)**2
            print *,d1*(x2-x3)**2+ 2*d2*(x2-x3)*(y2-y3)+ d3*(y2-y3)**2
            print *,d1*(x3-x1)**2+ 2*d2*(x3-x1)*(y3-y1)+ d3*(y3-y1)**2
            stop
         endif

         ra(k) = d1
         ra(k+melem) = d2
         ra(k+2*melem) = d3


!     etae(l) tvori vektor, z nich pak norma je etae(k) 
!         s=0
!         do 40 l=1,4
!            s=s+rmaxi(l)*rmaxi(l)
! 40      enddo
         
! nasobim konstantou, aby se magicke cislo dalo ve file.dt4 
! zapsat 
!       etae(k)=1000000000*sqrt(s)
!         etae(k)=10000*sqrt(s)
 10   enddo
      
      do 210 i=1,npoin
         rga(i) = 0.
         rgb(i) = 0.
         rgc(i) = 0.
         area(i) = 0.
 210  enddo
      
      do 220 i=1,nelem
         xc = (x(lnd(i,1))+x(lnd(i,2))+x(lnd(i,3)))/3
         yc = (y(lnd(i,1))+y(lnd(i,2))+y(lnd(i,3)))/3
         do 230 j=1,3
            k = lnd(i,j)
            rl = ((xc-x(k))**2 + (yc-y(k))**2)**0.5
            rga(k) = rga(k) + ra(i)/rl
            rgb(k) = rgb(k) + ra(i+melem)/rl
            rgc(k) = rgc(k) + ra(i+2*melem)/rl
            area(k) = area(k) + 1./rl
 230     enddo
 220  enddo
         
      do 240 i=1,npoin
         rga(i) = rga(i)/area(i)
         rgb(i) = rgb(i)/area(i)
         rgc(i) = rgc(i)/area(i)
         if( i .lt. -20) then
            write(*,'(i5,4e14.6)') i,rga(i),rgb(i),rgc(i),wp(i,1)
         endif
 240  enddo
      
      rkon = 750.
      do 260 i=1,npoin
         if(ibp(i,1) .gt. 0) then
            rga(i) = (rga(i) + rga(ibp(i,1)))/2
            rga(ibp(i,1)) = rga(i)
            rgb(i) = (rgb(i) + rgb(ibp(i,1)))/2
            rgb(ibp(i,1)) = rgb(i)
            rgc(i) = (rgc(i) + rgc(ibp(i,1)))/2
            rgc(ibp(i,1)) = rgc(i)
         endif
         if(abs(rga(i)) .gt. rkon .or. abs(rgc(i)) .gt. rkon) then
            write(ihel,'(2e14.6)') x(i),y(i)
         endif
            
 260  enddo
      

      ipr = -1
      rmax = 0.

      do 200 i=1,npoin
!      do 200 i=1,1
         if(abs(rgb(i)) .lt. 1E-05) then
!     the diagonal metrix
            rgb(i) = 0.
            rga(i) = (1. + AMA%epsilon1*abs(rga(i)))*rnorm
            rgc(i) = (1. + AMA%epsilon1*abs(rgc(i)))*rnorm
         else
!     the eigenvalues rlam1 , rlam2
            a = rga(i)
            b = rgb(i)
            c = rgc(i)

            disc = ((a - c)*(a - c) + 4*b*b)**0.5
            rlam1 = (a + c + disc)/2
            rlam2 = (a + c - disc)/2
            if(abs(rlam1) .lt. abs(rlam2) ) then
               rlampom = rlam1
               rlam1 = rlam2
               rlam2 = rlampom
            endif

            xm11 = b
            xm21 = -(a-rlam1)
            rl1 = (xm11*xm11+xm21*xm21)**0.5
            xm11 = xm11/rl1
            xm21 = xm21/rl1

            xm12 = b
            xm22 = -(a-rlam2)
            rl2 = (xm12*xm12+xm22*xm22)**0.5
            xm12 = xm12/rl2
            xm22 = xm22/rl2

            rdet = xm11*xm22 - xm21*xm12

!     inverse matrix
            ym11 = xm22/rdet
            ym12 = -xm12/rdet
            ym21 = -xm21/rdet
            ym22 = xm11/rdet

            z11 = abs(rlam1)*ym11
            z12 = abs(rlam1)*ym12
            z21 = abs(rlam2)*ym21
            z22 = abs(rlam2)*ym22
            
            t11 = xm11*z11 + xm12*z21
            t12 = xm11*z12 + xm12*z22
            t21 = xm21*z11 + xm22*z21
            t22 = xm21*z12 + xm22*z22

            if(i == ipr) then
               print *,'***********************',i,'    *************'
               print *,a,b,c,rlam1,rlam2
               print *,a-rlam1,b
               print *,b,c-rlam1
               print *,'            '
               print *,'----------------------'
               write(*,'(6e12.4)')  &
                    xm11,xm12,rlam1,0,ym11,ym12
               write(*,'(6e12.4)')  &
                    xm21,xm22,0,rlam2,ym21,ym22
               print *,'----------------------'
               write(*,'(4e12.4)')  &
                    z11,z12,t11,t12
               write(*,'(4e12.4)')  &
                    z21,z22,t21,t22
               print *,'----------------------'
            endif

            if(t11 .gt. rmax .or. t22 .gt. rmax) then
               rmax = max(t11,t22)
               imax = i
            endif

!            rkon = 500.
!            if(abs(t11) .gt. rkon .or. abs(t22) .gt. rkon) then
!               write(ihel,'(2e14.6)') x(i),y(i)
!            endif

!     Gamm channel
!            eps1 = 7500./(7500.+max(t11,t22))

!     profile
!            eps1 = 250./(10000.+max(t11,t22))

            eps1 = AMA%epsilon1

            rga(i) = (1. + eps1*t11)*rnorm 
            rgb(i) = (eps1*(t12+t21)/2)*rnorm
            rgc(i) = (1. + eps1*t22)*rnorm

            if(rga(i)*rgc(i) .le. rgb(i)*rgb(i) ) then
               print *,'ERROR in metric, nonpositive matrix'
               write(*,'(i5,3e14.6)') i,rga(i),rgb(i),rgc(i)
               print *,a,b,c
               print *,'----------------------'
               write(*,'(6e12.4)')  &
                    xm11,xm12,rlam1,0,ym11,ym12
               write(*,'(6e12.4)')  &
                    xm21,xm22,0,rlam2,ym21,ym22
               print *,'----------------------'
               write(*,'(4e12.4)')  &
                    z11,z12,t11,t12
               write(*,'(4e12.4)')  &
                    z21,z22,t21,t22
               print *,'----------------------'
               stop
            endif

         endif

!     for NACA we improve
         rimp = 1.

!         rrrl = ( x(i)*x(i) + y(i)*y(i))**0.5 + 
!     *        ( (x(i)-1.)*(x(i)-1.) + y(i)*y(i))**0.5 
!         if(rrrl .lt. 1.) then
!            print *,'??',x(i),y(i),rrrl
!            rrrl = 1.0
!         endif
!         rimp = ((rrrl*rrrl )/(rrrl*rrrl-0.75))*2

         rga(i) = rga(i)*rimp
         rgb(i) = rgb(i)*rimp
         rgc(i) = rgc(i)*rimp

!         write(*,'(i5,7e12.4)') 
!     *        i,x(i),y(i),rrrl,rimp,rga(i),rgb(i),rgc(i)
            
!         if( i .gt. 80 .and. i .lt. 6130) then
         if( i .lt. -20) then
            write(*,'(i5,4e14.6)') i,rga(i),rgb(i),rgc(i),wp(i,1)
         endif
 200  enddo

!      print *,'MAX =',imax,rmax

!     smoothing
      do 500 ic = 1,6
         do 400 i=1,nelem
            ra(i) = (rga(lnd(i,1)) + rga(lnd(i,2)) + rga(lnd(i,3)))/3
            ra(i+melem) =   &
                 (rgb(lnd(i,1)) + rgb(lnd(i,2)) + rgb(lnd(i,3)))/3
            ra(i+2*melem) =   &
                 (rgc(lnd(i,1)) + rgc(lnd(i,2)) + rgc(lnd(i,3)))/3
 400     enddo
         
         do 410 i=1,npoin
            rga(i) = 0.
            rgb(i) = 0.
            rgc(i) = 0.
            area(i) = 0.
 410     enddo
      
         do 420 i=1,nelem
            xc = (x(lnd(i,1))+x(lnd(i,2))+x(lnd(i,3)))/3
            yc = (y(lnd(i,1))+y(lnd(i,2))+y(lnd(i,3)))/3
            do 430 j=1,3
               k = lnd(i,j)
               rl = ((xc-x(k))**2 + (yc-y(k))**2)**0.5
               rga(k) = rga(k) + ra(i)/rl
               rgb(k) = rgb(k) + ra(i+melem)/rl
               rgc(k) = rgc(k) + ra(i+2*melem)/rl
               area(k) = area(k) + 1./rl
 430        enddo
 420     enddo
         
         do 440 i=1,npoin
            rga(i) = rga(i)/area(i)
            rgb(i) = rgb(i)/area(i)
            rgc(i) = rgc(i)/area(i)
            if( i .lt. -20) then
               write(*,'(i5,4e14.6)') i,rga(i),rgb(i),rgc(i),wp(i,1)
            endif
 440     enddo
 500  enddo
     
      do 63 i=1,npoin
         if(ibp(i,1) .gt. 0) then
            rga(i) = (rga(i) + rga(ibp(i,1)))/2
            rga(ibp(i,1)) = rga(i)
            rgb(i) = (rgb(i) + rgb(ibp(i,1)))/2
            rgb(ibp(i,1)) = rgb(i)
            rgc(i) = (rgc(i) + rgc(ibp(i,1)))/2
            rgc(ibp(i,1)) = rgc(i)
         endif
 63   enddo

!      if(noit == -1) then
!         imet = 44
!         open(imet, file='metric',status='unknown')
!         write(imet,*) 'The metric in each point:'
!         do 350 i=1,npoin
!            write(imet,'(i5,4e14.6)') i,rga(i),rgb(i),rgc(i),wp(i,1)
! 350     enddo
!      endif

      return
    end subroutine ERROR_77

    subroutine SOLVE6_1_77(rmat,bmat,xmat)
      real*8 rmat(6,6),bmat(6),xmat(6),xmatold(6),sum,err,errold
      
      tau = 0.0005
      do i=1,6
         xmatold(i) = 0.D+00
      enddo

      errold = 100000.
      it = 1
10    it = it + 1
      err = 0.D+00
      do i=1,6
         sum = 0.D+00
         do j=1,6
            sum = sum + rmat(i,j)*xmatold(j)
         enddo
         xmat(i) = xmatold(i) - tau*(sum - bmat(i))
         err = err + (xmat(i)-xmatold(i))**2
         xmatold(i) = xmat(i)
      enddo
      !         do j=1,6
      !            err = err + (xmat(j)-xmatold(j))**2
      !            xmatold(j) = xmat(j)
      !         enddo
      if(mod(it,5000) == 1)   &
           write(*,'(i10,7e11.3)') it,(xmat(ij), ij=1,6),err
      if(err .lt. 1E-25 .and. it .gt. 1) goto 20
      errold = err
      goto 10
20    continue
      return
    end subroutine SOLVE6_1_77

    subroutine SOLVE6_77(rmat,bmat,xmat)
      real*8 rmat(6,6),bmat(6),xmat(6),xmatold(7),quot
      
      do i=1,6
         quot = rmat(i,i)
         if(quot == 0. ) then
            ip = i - 1
11          ip = ip + 1 
            if(quot == 0. .and.  ip .lt. 6) then
               !     reordering
               do l=1,6
                  xmatold(l) = rmat(i,l)
               enddo
               xmatold(7) = bmat(i)
               do k=i,5
                  do l=1,6
                     rmat(k,l) = rmat(k+1,l)
                  enddo
                  bmat(k) = bmat(k+1)
               enddo
               do l=1,6
                  rmat(6,l) = xmatold(l) 
               enddo
               bmat(6) = xmatold(7) 
               quot = rmat(i,i)
               goto 11
            elseif(quot == 0 .and. ip == 6) then
               print *,'Singular matrix'
               stop
            endif
         endif
         
         bmat(i) = bmat(i)/quot
         do j = i,6
            rmat(i,j) = rmat(i,j)/quot
         enddo
         do j= 1,6
            if(i .ne. j) then
               quot = rmat(j,i)
               do l =1,6
                  rmat(j,l) = rmat(j,l) - quot*rmat(i,l)
               enddo
               bmat(j) = bmat(j) - quot*bmat(i)
            endif
         enddo
         !         do ii = 1,6
         !            write(*,'(7f9.4)') (rmat(ii,ij), ij=1,6),bmat(ii)
         !         enddo
         !         print *,'----------------------------------------'
      enddo
      do j=1,6
         xmat(j) = bmat(j)
      enddo
      return
    end subroutine SOLVE6_77



    subroutine ERROR1_77(ndim, melem,nelem,mpoin,npoin,  &
         maxdeg,dx,dy,area,tria,noit,  &
         x,y,lnd,iae,wp,icyc,ra,rga,rgb,rgc,surface,  &
         ibp,ipoint,nbp,ibb,ibpoin,xb,yb,ifv, iwa, iwall,  &
         nbelm,mbelm,lbn, ibc )
      dimension lnd(melem,3),iae(melem,3),x(mpoin),y(mpoin),  &
           ra(melem*3),wp(mpoin,ndim+1),dx(melem),dy(melem),  &
           icyc(mpoin,maxdeg),rga(mpoin),rgb(mpoin),rgc(mpoin),  &
           area(mpoin),tria(melem),ibp(mpoin,2),ibb(mpoin,3),  &
           xb(ipoint),yb(ipoint),ibpoin(0:nbp), iwall(1:10),  &
           lbn(mbelm,2),ibc(mbelm)
      real*8 dxi,dyi,rmeas,x1,x2,x3,y1,y2,y3,w1,w2,w3,areai,rl
      real*8 a,b,c,disc,rlam1,rlam2,x11,x12,x21,x22,y12,y11,y21,y22,  &
           t11,t12,t21,t22,z11,z12,z21,z22,rdet,s,t
      real*8 rmat(6,6),bmat(6),xmat(6)
      real*8 sumx, sumy, sumx2, sumxy, sumy2, sumx3, sumx2y, sumxy2,  &
           sumy3, sumx4, sumx3y, sumx2y2, sumxy3, sumy4, sumu, sumux,  &
           sumuy, sumux2, sumuxy, sumuy2

      rkappa = 1.4

!     dx the first derivative dw/dx
!     dy                      dw/dy
!     rga  the second der.    d2w/dx dx
!     rgb                     d2w/dx dy
!     rg!                     d2w/dy dy       

      do 10 i=1,nelem
         x1 = x(lnd(i,1))
         y1 = y(lnd(i,1))
         w1 = wp(lnd(i,1),1)
         x2 = x(lnd(i,2))
         y2 = y(lnd(i,2))
         w2 = wp(lnd(i,2),1)
         x3 = x(lnd(i,3))
         y3 = y(lnd(i,3))
         w3 = wp(lnd(i,3),1)
         rmeas = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
         tria(i) = rmeas/2
         dx(i) = (y2-y1)*(w1+w2) + (y3-y2)*(w2+w3) + (y1-y3)*(w1+w3) 
         dy(i) = (x1-x2)*(w1+w2) + (x2-x3)*(w2+w3) + (x3-x1)*(w1+w3) 
         dx(i) = dx(i) /rmeas
         dy(i) = dy(i) /rmeas

!         write(25,*) (x1+x2+x3)/3, (y1+y2+y3)/3, dx(i), dy(i)

!         write(26,*) x(lnd(i,1)), y(lnd(i,1)), wp(lnd(i,1),1)
!         write(26,*) x(lnd(i,2)), y(lnd(i,2)), wp(lnd(i,2),1)
!         write(26,*) x(lnd(i,3)), y(lnd(i,3)), wp(lnd(i,3),1)
 10   enddo


      rmaxder = 0.
!     HERE INSERT deriv.f

      itest = -20

      err4 = 0.
      err5 = 0.
      err6 = 0.

      do 20 i=1,npoin
         areai = 0.
         rga(i) = 0.
         rgb(i) = 0.
         rgc(i) = 0.
         rgdi = 0.
         if(icyc(i,1) .gt. 0) then
            len = icyc(i,1) 
         else
            len = -icyc(i,1)  - 1
         endif
         x3 = x(i)
         y3 = y(i)
         w3 = wp(i,1)
         do 30 j=1,len
            j1 = mod(j, abs(icyc(i,1))) + 1
            x1 = x(icyc(i,j+1))
            y1 = y(icyc(i,j+1))
            w1 = wp(icyc(i,j+1),1)
            x2 = x(icyc(i,j1+1))
            y2 = y(icyc(i,j1+1))
            w2 = wp(icyc(i,j1+1),1)
            rmeas = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
            areai = areai + rmeas/6
            dxi = (y2-y1)*(w1+w2) + (y3-y2)*(w2+w3) + (y1-y3)*(w1+w3) 
            dyi = (x1-x2)*(w1+w2) + (x2-x3)*(w2+w3) + (x3-x1)*(w1+w3) 
            dxfi = (y1-y2)/rmeas
            dyfi = (x2-x1)/rmeas

            if(i == itest) write(*,'(a2,i2,3e14.6,2e12.4)')  &
                 '$$',j,w1,w2,w3,dxi/rmeas,dyi/rmeas
            if(i == itest) write(*,'(a2,4i5,4e12.4)')  &
                 '  ',j,icyc(i,j+1),icyc(i,j1+1),i,x1,y1,x2,y2
            if(i == itest) write(*,'(a2,3e14.6)') '  ',x1,y1,w1
            if(i == itest) write(*,'(a2,3e14.6)') '  ',x2,y2,w2
            if(i == itest) write(*,'(a2,3e14.6)') '  ',x3,y3,w3

            rga(i) = rga(i) - dxi/2*dxfi
            rgb(i) = rgb(i) - dyi/2*dxfi
            rgc(i) = rgc(i) - dyi/2*dyfi
            rgdi = rgdi - dxi/2*dyfi

!            write(100 + AMA%adapt_level,*) (x1+x2+x3)/3, (y1+y2+y3)/3, 
!     *           dxi/rmeas, dyi/rmeas
            
!            if(i == itest) print *,'**',j,dxi,dxfi,rga(i)

            if(icyc(i,1) .lt. 0 .and. ibp(i,1) == 0 ) then
!     for the boundary points me must add the boundary sides
               if(j == 1) then
                  rga(i) = rga(i) + dxi/rmeas*(y1-y3)/2
                  rgb(i) = rgb(i) + dyi/rmeas*(y1-y3)/2
                  rgc(i) = rgc(i) + dyi/rmeas*(x3-x1)/2
                  rgdi = rgdi + dxi/rmeas*(x3-x1)/2
               endif
               if( j == len) then
                  rga(i) = rga(i) + dxi/rmeas*(y3-y2)/2
                  rgb(i) = rgb(i) + dyi/rmeas*(y3-y2)/2
                  rgc(i) = rgc(i) + dyi/rmeas*(x2-x3)/2
                  rgdi = rgdi + dxi/rmeas*(x2-x3)/2

               endif
            endif


!            if(i == itest) print *,'**',j,dxi,dxfi,rga(i)

 30      enddo
         area(i) = areai
!         if(icyc(i,1) .gt. 0) then
!         rradius = (x(i)*x(i) + y(i)*y(i))**0.5
!         if(rradius .gt. 0) then
!            derxxtrue = 3./16*rradius**(-1.75)
!            deryytrue = derxxtrue
!            
!            derxxtrue = -(1-y(i)**20)*90*x(I)**8
!            deryytrue = -(1-x(i)**10)*380*y(i)**18
!            derxytrue = 200*x(i)**9*y(i)**18
!            
!            era1 = ( ( rga(i) ) /area(i)-derxxtrue)
!            era2 = ( ( rgc(i) ) /area(i)-deryytrue)
!            era3 = ( ( rgb(i) ) /area(i)-derxytrue)
!            err4 = err4 + (era1*era1 + era2*era2 + 
!     *           2.*era3*era3)*area(i)
!            err5 = err5 + (derxxtrue**2 + derxytrue**2 +
!     *           2.*deryytrue**2)*area(i)
!            err6 = err6 + ( rga(i)**2 + 2*rgb(i)**2 +
!     *           rgc(i)**2)/area(i)
!            write(98,'(i5,7e12.4)') i,area(i),
!     *        rga(i)/area(i),derxxtrue,
!     *        (era1*era1 + era2*era2)*area(i),
!     *           (derxxtrue**2 + deryytrue**2)*area(i),
!     *           err4,err5,err4/err5,err6
!            endif
!         endif

         if(AMA%ityp == -10) then
            rga(i) = wp(i,1)
            rgb(i) = 0.
            rgc(i) = rga(i)
         endif

         write(200 + AMA%adapt_level,*) x(i), y(i),   &
              rga(i)/area(i), rgb(i)/area(i), rgc(i)/area(i)

 20   enddo

!     for periodic boundary points
      if( AMA%xper(1,1) .gt. 0 .or. AMA%xper(1,2) .gt. 0) then
         do 18 i1=1,npoin
            i2 = ibp(i1,1)
            if( i2 .gt. 0 ) then
               if(rga(i1) .ne. rga(i2) .or. rgb(i1) .ne. rgb(i2) .or.   &
                    rgc(i1) .ne. rgc(i2) .or. area(i1) .ne. area(i2) )  &
                    then
                  rga(i1) = rga(i1) + rga(i2)
                  rgb(i1) = rgb(i1) + rgb(i2)
                  rgc(i1) = rgc(i1) + rgc(i2)
                  area(i1) = area(i1) + area(i2)
                  rga(i2) = rga(i1) 
                  rgb(i2) = rgb(i1) 
                  rgc(i2) = rgc(i1) 
                  area(i2) = area(i1)
               endif
            endif
 18      enddo
      endif

!      print *,'Total error = ',err4**0.5,err5**0.5,
!     *     err4**0.5/err5**0.5
!      print *,'Total error = ',(err4/err6)**0.5

      do 15 i=1,npoin
         rga(i) = rga(i)/area(i)
         rgb(i) = rgb(i)/area(i)
         rgc(i) = rgc(i)/area(i)
         if(i == itest) print *,'..',area(i),rga(i)
!         if(abs(rga(i)) + abs(rgc(i)) .gt. 5.) 
!         if(x(i) .gt. 0.95) 
 15   enddo


!     ... comparing real second derivatives with their approximation
      do i=1,npoin
         rradius = (x(i)*x(i) + y(i)*y(i))**0.5
         if(rradius .gt. 0) then
            derxxtrue = 3./16*rradius**(-1.75)
         else
            derxxtrue = 1E+25
         endif
         deryytrue = derxxtrue

!         derxxtrue = -(1-y(i)**20)*90*x(I)**8
!         deryytrue = -(1-x(i)**10)*380*y(i)**18

         era1 = (abs(rga(i))-derxxtrue)
         era2 = (abs(rgc(i))-deryytrue)
!         write(99,'(6e14.6)')
!     *        x(i),y(i),rga(i),derxxtrue,rgc(i),deryytrue
!     *        x(I),y(i),era1,era2,rga(i),rgc(i),derxxtrue
!         rat = 250.*exp(-1*rradius)
!         rga(i) = rat
!         rgb(i) = rat
      enddo



      rmaxder = 0.
!     improvement for boundary points
      goto 101

      do 100 i=1,npoin
         if(icyc(i,1) .lt. 0 .and. ibp(i,1) == 0) then
            len = abs(icyc(i,1))
            if( len .ge. 3) then
               ipoc = 0
               rga(i) = 0.
               rgb(i) = 0.
               rgc(i) = 0.
               rminlenght = 1E+38
               do 110 j=2,len-1
                  if(icyc(icyc(i,j+1),1) .gt. 0) then
!     this point isn't boundary
                     xlen = ((x(i)-x(icyc(i,j+1)))**2 +  &
                          (y(i) -y(icyc(i,j+1)))**2 )**0.5
                     if(xlen .lt. rminlenght) then
                        rga(i) = rga(icyc(i,j+1))
                        rgb(i) = rgb(icyc(i,j+1))
                        rgc(i) = rgc(icyc(i,j+1))
                        rminlenght = xlen
                     endif
                     if( i == -170) then
                        print *,'**',j,icyc(i,j+1),rgc(icyc(i,j+1)),  &
                             xlen,rminlenght 
                     endif
                  endif
 110           enddo
!     we use the value from the nearest point
               if(ipoc .gt. 0 ) then
                  rga(i) = rga(i)/ipoc
                  rgb(i) = rgb(i)/ipoc
                  rgc(i) = rgc(i)/ipoc
               else
!     in this case, we let the second derivations = 0 !!!!

               endif
            elseif( len == 2) then
               i1 = icyc(i,2)
               i2 = icyc(i,3)
               do 120 ie =1,nelem
                  do 130 je =1,3
                     je1 = mod(je,3)+1
                     if(lnd(ie,je) == i2 .and.   &
                          lnd(ie,je1) == i1) then
                        je2 = mod(je1,3) + 1
                        rga(i) = rga(lnd(ie,je2))
                        rgb(i) = rgb(lnd(ie,je2))
                        rgc(i) = rgc(lnd(ie,je2))
                        goto 140
                     endif
 130              enddo
 120           enddo
 140           continue
            else
               print *,'ERROR len < 2 !!!!'
            endif
         endif
         rmaxder = max(rmaxder,abs(rga(i)),abs(rgb(i)),abs(rgc(i)) )
 100  enddo
 101  continue

      rnorm = 1.29903810567 *AMA%numel/surface
      rmaxderteor = AMA%epsilon1/3*rnorm
      der = rmaxderteor
      epsilon2p = AMA%epsilon1/AMA%p
!      epsilon2p = p

      epsilon2 = max (1.,epsilon2p)

      write(AMA%ifig1,*)'   '
      write(AMA%ifig1,*)'Given data:'
      write(AMA%ifig1,*)'Used component for adaptation:',AMA%ityp,  &
           '  (0-uniform mesh)'
      write(AMA%ifig1,*)'Dimension of the solution:    ',ndim
      write(AMA%ifig1,*)'Solution associated:          ',ifv,  &
           '  (0-cell vertex, 1-cell centered)'
      write(AMA%ifig1,*)'Positivity:                   ',AMA%pos
      write(AMA%ifig1,*)'Prescribed number of elements:',AMA%numel
      write(AMA%ifig1,*)'epsilon1:                     ',AMA%epsilon1
      write(AMA%ifig1,*)'p:                            ',p

      ipr = -1
!     we compute the elements of matrix M from the second derivations 

      rmax1 = 0.
      rmax2 = 0.
      rmax3 = 0.
      rmaxlambdas = 0.
      imaxl = 1

!      itest = 50

!      xmin1 = -0.15
!      xmax1 = 0.05
!      ymax1 = 0.01
!      rmaa = 0.
!      rmac = 0.
!      rmac = 0.
!      do i=1,npoin
!         if(x(i) .gt. xmin1 .and. x(i).lt. xmax1 .and.
!     *        y(i) .lt. ymax1) then
!            rmaa = max(rmaa, rga(i) )
!            rmab = max(rmaa, rgb(i) )
!            rmac = max(rmaa, rgc(i) )
!            write(49,*) x(i),y(i),rmaa,rmab,rmac
!         endif
!      enddo
!      do i=1,npoin
!         if(x(i) .gt. xmin1 .and. x(i).lt. xmax1 .and.
!     *        y(i) .lt. ymax1) then
!            rga(i) = 2*rmaa
!            rgb(i) = 0.
!            rgc(i) = 2*rmac
!         endif
!      enddo
!      close(49)


      do 200 i=1,npoin

         if(rga(i) .gt. rmax1) then
            rmax1 = rga(i)
            imax1 = i
         endif
         if(rgb(i) .gt. rmax2) then
            rmax2 = rgb(i)
            imax2 = i
         endif
         if(rgc(i) .gt. rmax3) then
            rmax3 = rgc(i)
            imax3 = i
         endif

         if(abs(rgb(i)) .lt. 1E-05) then
!     the diagonal matrix

            eps1 = AMA%epsilon1/(epsilon2+max(rga(i),rgc(i)))

            rgb(i) = 0.
            rga(i) = (1. + eps1*abs(rga(i)))*rnorm
            rgc(i) = (1. + eps1*abs(rgc(i)))*rnorm
            if(rga(i) + rgc(i) .gt. rmaxlambdas)   &
                 rmaxlambdas = rga(i)+rgc(i)
         else
!     the eigenvalues rlam1 , rlam2
            if(i == itest) then
               print *,x(i),y(i)
               print *,'@#',rga(i),rgb(i),rgc(i)
            endif
            a = rga(i)
            b = rgb(i)
            c = rgc(i)

            disc = ((a - c)*(a - c) + 4*b*b)**0.5
            rlam1 = (a + c + disc)/2
            rlam2 = (a + c - disc)/2
            if(abs(rlam1) .lt. abs(rlam2) ) then
               rlampom = rlam1
               rlam1 = rlam2
               rlam2 = rlampom
            endif
            if(abs(rlam1) + abs(rlam2) .gt. rmaxlambdas) then
               rmaxlambdas = abs(rlam1) + abs(rlam2)
               imaxl = i
            endif
            x11 = b
            x21 = -(a-rlam1)
            rl1 = (x11*x11+x21*x21)**0.5
            x11 = x11/rl1
            x21 = x21/rl1

            x12 = b
            x22 = -(a-rlam2)
            rl2 = (x12*x12+x22*x22)**0.5
            x12 = x12/rl2
            x22 = x22/rl2

            rdet = x11*x22 - x21*x12

!     inverse matrix
            y11 = x22/rdet
            y12 = -x12/rdet
            y21 = -x21/rdet
            y22 = x11/rdet

            z11 = abs(rlam1)*y11
            z12 = abs(rlam1)*y12
            z21 = abs(rlam2)*y21
            z22 = abs(rlam2)*y22
            
            t11 = x11*z11 + x12*z21
            t12 = x11*z12 + x12*z22
            t21 = x21*z11 + x22*z21
            t22 = x21*z12 + x22*z22

            eps1 = AMA%epsilon1/(epsilon2+max(t11,t22))

            if(i == itest) print *,' #',eps1,epsilon2,rnorm
            if(i == itest) print *,' #',eps1*t11+1.,rnorm

            rga(i) = (1. + eps1*t11)*rnorm 
            rgb(i) = (eps1*(t12+t21)/2)*rnorm
            rgc(i) = (1. + eps1*t22)*rnorm

            if(i == itest) print *,'@#',rga(i),rgb(i),rgc(i)

            if(rga(i)*rgc(i) .le. rgb(i)*rgb(i) ) then
               print *,'ERROR in metric, nonpositive matrix'
               write(*,'(i5,3e14.6)') i,rga(i),rgb(i),rgc(i)
               print *,a,b,c
               print *,'----------------------'
               write(*,'(6e12.4)')  &
                    x11,x12,rlam1,0,y11,y12
               write(*,'(6e12.4)')  &
                    x21,x22,0,rlam2,y21,y22
               print *,'----------------------'
               write(*,'(4e12.4)')  &
                    z11,z12,t11,t12
               write(*,'(4e12.4)')  &
                    z21,z22,t21,t22
               print *,'----------------------'
               stop
            endif
            !write(23,*) x(i),y(i),rga(i),rgb(i),rgc(i)
         endif
 200  enddo
!      print *,'###############',rmaxlambdas,imaxl,
!     *     x(imaxl),y(imaxl)

!     CORRECTION of positivity a priori
      par = 1.0
      do i=1,npoin
         a = rga(i)
         b = rgb(i)
         c = rgc(i)
         disc = ((a - c)*(a - c) + 4*b*b)**0.5
         rlam1 = (a + c + disc)/2
         radius = 4*AMA%pos*AMA%pos*rlam1/par/par
!         rh0 = 1./(rlam1)**0.5
!         rh1 = rh0/3/AMA%pos*par
!         radius = 1/rh1/rh1
!     rgai,rgbi,rgci - cyrcle whose radius is the pos times the smallest
!     side
         rgai = radius
         rgbi = 0.
         rgci = radius
         if(i == -20) then
            rgai = 4.
            rgbi = 2.
            rgci = 12.
            rga(i) = 3. 
            rgb(i) = 1.5
            rgc(i) = 9.
            print *,x(i),y(i),rlam1,radius
            print *,rga(i),rgb(i),rgc(i)
            print *,rgai,rgbi,rgci
            call ELIPS_77(14,rga(i),rgb(i),rgc(i),x(i),y(i) )
            call ELIPS_77(14,rgai,rgbi,rgci,x(i),y(i) )
          endif
          if(i == itest) print *,'$$',rga(i),rgb(i),rgc(i)
          
          call INTERSECTION_METRIC_77(rga(i),rgb(i),rgc(i),  &
               rgai,rgbi,rgci)
          if(i == itest) print *,'$2',rga(i),rgb(i),rgc(i)
!        if(i == 20) then
!            call ELIPS_77(14,rga(i),rgb(i),rgc(i),x(i),y(i) )
!            print *,rga(i),rgb(i),rgc(i)
!            stop
!         endif
      enddo
!     end of CORRECTION of positivity a priori

!     smoothing
      print*,'################### smoothing', 0
      do 500 ic = 1,0
!!      do 500 ic = 1,6
         do 400 i=1,nelem
            ra(i) = (rga(lnd(i,1)) + rga(lnd(i,2)) + rga(lnd(i,3)))/3
            ra(i+melem) =   &
                 (rgb(lnd(i,1)) + rgb(lnd(i,2)) + rgb(lnd(i,3)))/3
            ra(i+2*melem) =   &
                 (rgc(lnd(i,1)) + rgc(lnd(i,2)) + rgc(lnd(i,3)))/3
 400     enddo
         
         do 410 i=1,npoin
            rga(i) = 0.
            rgb(i) = 0.
            rgc(i) = 0.
            area(i) = 0.
 410     enddo
      
         do 420 i=1,nelem
            xc = (x(lnd(i,1))+x(lnd(i,2))+x(lnd(i,3)))/3
            yc = (y(lnd(i,1))+y(lnd(i,2))+y(lnd(i,3)))/3
            do 430 j=1,3
               k = lnd(i,j)
!               rl = ((xc-x(k))**2 + (yc-y(k))**2)**0.5
               rga(k) = rga(k) + ra(i)*tria(i)
               rgb(k) = rgb(k) + ra(i+melem)*tria(i)
               rgc(k) = rgc(k) + ra(i+2*melem)*tria(i)
               area(k) = area(k) + tria(i)
 430        enddo
 420     enddo
         
         do 440 i=1,npoin
            rga(i) = rga(i)/area(i)
            rgb(i) = rgb(i)/area(i)
            rgc(i) = rgc(i)/area(i)
         if(i == itest) print *,'$3',rga(i),rgb(i),rgc(i)
 440     enddo
 500  enddo

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


      if((AMA%ityp == 3 .or. AMA%ityp == 4 .or. AMA%ityp == 5) .and.  &
           AMA%Re .gt. 0.) then


!     part of courved boundary, where is fixed walls
         print*,'$$$$$$$$$$$$$  ATTENTION HERE'

         ilkw = 1

!     generation the triangulation for the boundary layer for NS
         qbas = 1.25
!     charakterictic length
         xmax = -100000.
         ymax0 = -100000.
         xmin = 100000.
         ymin = 100000.
!         do i=1,ipoint
!            if( i .gt. ibpoin(ilkw-1) .and. 
!     *           i .le. ibpoin(ilkw)) then
!               xmax = max(xmax, xb(i) )
!               ymax0 = max(ymax0, yb(i) )
!               xmin = min(xmin, xb(i) )
!               ymin = min(ymin, yb(i) )
!            endif
!         enddo
         
!         print*, xmax, ymax0

         do i=1,nbelm
            do j=1,iwa
               if(ibc(i) == iwall(j) )then
                  xmax = max(xmax, x(lbn(i,1)), x(lbn(i,2)) )
                  ymax0 = max(ymax0, y(lbn(i,1)), y(lbn(i,2))  )
                  xmin = min(xmin, x(lbn(i,1)), x(lbn(i,2))  )
                  ymin = min(ymin, y(lbn(i,1)), y(lbn(i,2))  )
               endif
            enddo 
         enddo

!         print*, xmin, xmax, ymin, ymax0

         rcharlen = ((xmax -xmin)**2 + (ymax0 - ymin)**2)**0.5
!     double precision
         x1 = AMA%xte(1,1)
         y1 = AMA%xte(1,2)
         x2 = AMA%xte(2,1)
         y2 = AMA%xte(2,2)

!     blasius
!     DMR
!         rcharlen = 1.
         
         Resminside = rmaxlambdas/3
         Res = min(AMA%Re, 1000000*Resminside)
!         Res = min(Re, 10000*Resminside)
         if(AMA%Re .gt. Res) print*,'Re  <  10000*Resminside !!!'
         Res = Res/rcharlen

         print*,'Reynolds =',Res,'  (Resminside, Re):',Resminside,AMA%Re

!         ymax = ((3/rnorm)**0.5 - 1./sqrt(Res))/(qbas-1.)*rcharlen
         ymax1 = ((3/rnorm)**0.5 - 1./sqrt(Res))/(qbas-1.)*rcharlen
         ymax = 20./sqrt(Res)
!         ymax = 60./sqrt(Res)
         n0 = 3

         imt = 18
         open(imt,file='mt0',status='unknown')

!         print*,'$$$',rcharlen, xmin,xmax, ymax0, ymin,npoin

         do i=1,npoin
            xi = x(i)
            yi = y(i)
            rmax = 100000.

            do ib=1,nbelm
               do j=1,iwa
                  if(ibc(ib) == iwall(j) )then

                     xj = (x(lbn(ib,1))+ x(lbn(ib,1)))/2
                     yj = (y(lbn(ib,1))+ y(lbn(ib,1)))/2

                     rlen = (xi-xj)**2 + (yi-yj)**2
                     if(rlen .lt. rmax) then
                        rmax = rlen
                        imax = ib
                     endif

                  endif
               enddo 
            enddo


!           do k=ibpoin(ilkw-1)+1,ibpoin(ilkw)
!               rlen = (xi-xb(k))**2 + (yi-yb(k))**2
!               if(rlen .lt. rmax) then
!                  rmax = rlen
!                  imax = k
!               endif
!            enddo
!     distance from profile
            rmax = rmax**0.5

!     channel
!            rmax = min ( (1- y(i)), y(i) )
!            imax = -10
!            n0 = 1
!     blasius
!            if(xi .ge. 0 ) then
!               rmax = abs(yi)
!            else
!               rmax = (xi*xi + yi*yi/10.)**0.5
!            endif

            x3 = xi
            y3 = yi

            sstt = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)
            s = ( (x1-x2)*(x1-x3) + (y1-y2)*(y1-y3) )/sstt
            t = ( x1*(y3-y2) + x2*(y1-y3) + x3*(y2-y1) )/sstt
!     wake and boundary layer are investigated separetely,
!     as some points need as to BL as to Wake


!     BOUNDARY LAYERS            
!            print*,'$$$',i,s,rcharlen, rmax, ymax

            if(s .le. 0.01*rcharlen) then

!               if(rmax .gt. ymax .and. rmax .lt. 5*ymax)
!                   write(*,*) x(i),y(i),rmax,ymax
               rwake = 2*rmax + 1.
               if( rmax .lt. ymax) then
!               if( rmax .lt. 5*ymax) then
                  cor = 1./n0

!                  if(imax .ge. 0) then
!                     if(imax .gt. 1) then
!                        ii1 = imax - 1
!                     else
!                        ii1 = ipoint - 1
!                     endif
!                     rn1 = -(yb(ii1) - yb(imax))
!                     rn2 = xb(imax) - xb(ii1) 
!                  endif

                  rn1 = y(lbn(imax,2)) - y(lbn(imax,1))
                  rn2 = x(lbn(imax,1)) - x(lbn(imax,2))

!                  write(*,'(i5,6e12.4)') i, xi,yi, 
!     *                 x(lbn(imax,1)), y(lbn(imax,1)),rn1,rn2

!                  q = 5.*qbas
                  q = (0.2+qbas)/1.2 -ymax/20.
                  q = max(q, 0.)

                  rlen = (rn1*rn1 + rn2*rn2)**0.5
                  rn1 = rn1/rlen
                  rn2 = rn2/rlen

               rh2 = cor/sqrt(Res) + rmax*(q - 1.)/rcharlen

!               print*,'###',rh2,cor, rmax, q, rn1, rn2
!                  rh2 = cor/sqrt(Res)+ (rcharlen/10. - cor/sqrt(Res))*
!     *                 rmax/ymax
                  rh2 = rh2*4.
                  rh1 = rh2*(1+(1-24.*AMA%pos*AMA%pos)**0.5)/4/AMA%pos
!     DMR
!                  rh1 = rh2

!                  write(*,*) x(i),y(i),q,rh1,rh2
                  rgai = rn2*rn2/rh1/rh1 + rn1*rn1/rh2/rh2
                  rgci = rn2*rn2/rh2/rh2 + rn1*rn1/rh1/rh1
                  rgbi = rn1*rn2*(1./rh1/rh1 - 1./rh2/rh2)
                  call INTERSECTION_METRIC_77(rga(i),rgb(i),rgc(i),  &
                       rgai,rgbi,rgci)

!     for plotting of the elipses
!                  if(xi .gt. -0.01 .and. xi .lt. 0.2 .and. 
!     *                 abs(y(i)) .lt. 0.02) then
!                     call ELIPS_77(imt,rga(i),rgb(i),rgc(i),x(i),y(i) )
                     call ELIPS_77(imt,rgai,rgbi,rgci,x(i),y(i) )
!                     print*,'##',i,xi,yi,rmax, ymax
!                  endif

               endif
            endif

!     WAKES
            if(s .ge. 0.0) then
               xk = x1 + s*(x2-x1)
               yk = y1 + s*(y2-y1)
               rwake = ( ( xk - x3)**2 + (yk - y3)**2)**0.5
               if( rwake .lt. ymax) then
!                  cor = (2.+ s*4.)*0.15/n0
                  rlenrel = s*sqrt(sstt)/rcharlen
!                  cor = (1.+rlenrel)/n0*2.
                  cor = (1.+rlenrel)/n0*6.
!                  write(46,*) xi,yi,cor
                  rn1 = -(AMA%xte(1,2) - AMA%xte(2,2))
                  rn2 = AMA%xte(2,1) - AMA%xte(1,1)
                  rmax = rwake-0.02*rcharlen
                  rmax = max (rmax,0.)
!                  q = (0.2+qbas)/1.2 - ymax/20.
!                  q = max(q, 0.)
!                  q = qbas
                  rlen = (rn1*rn1 + rn2*rn2)**0.5
                  rn1 = rn1/rlen
                  rn2 = rn2/rlen

!                  rh2 = cor/sqrt(Res) + rmax*(q - 1.)/rcharlen
                  rh2 = cor/sqrt(Res)+ (rcharlen/10. - cor/sqrt(Res))*  &
                       rmax/ymax
                  rh1 = rh2*(1+(1-24.*AMA%pos*AMA%pos)**0.5)/4/AMA%pos
!               write(33,*) x(i),y(i),q,rh1,rh2
                  rgai = rn2*rn2/rh1/rh1 + rn1*rn1/rh2/rh2
                  rgci = rn2*rn2/rh2/rh2 + rn1*rn1/rh1/rh1
                  rgbi = rn1*rn2*(1./rh1/rh1 - 1./rh2/rh2)

!               call ELIPS_77(64,rgai,rgbi,rgci,x(i),y(i) )
                  call INTERSECTION_METRIC_77(rga(i),rgb(i),rgc(i),  &
                       rgai,rgbi,rgci)
               endif
            endif
         enddo

      close(imt)

      endif
!      stop


      icontrol = -1
      if(icontrol == 1) then
!     MESH GRADIENT CONTROL H-variation
         alpha = 2.0
         do is =1,4
            do k=1,3*nelem
               ra(k) = 1.
            enddo
            do i=1,nelem
               do j=1,3
                  if(ra(3*(i-1)+j) .gt. 0) then
                     j1 = mod(j,3) + 1
                     k1 = lnd(i,j)
                     k2 = lnd(i,j1)
                     x1 = x(k1)
                     y1 = y(k1)
                     x2 = x(k2)
                     y2 = y(k2)
                     rga1 = rga(k1)
                     rgb1 = rgb(k1)
                     rgc1 = rgc(k1)
                     rga2 = rga(k2)
                     rgb2 = rgb(k2)
                     rgc2 = rgc(k2)
                     rl1 = (rga1*(x1-x2)**2 + 2*rgb1*(x1-x2)*(y1-y2) +  &
                          rgc1*(y1-y2)**2)**0.5
                     rl2 = (rga2*(x1-x2)**2 + 2*rgb2*(x1-x2)*(y1-y2) +  &
                          rgc2*(y1-y2)**2)**0.5
                     rga3 = rga2*(1+alpha*rl1)**(-2)
                     rgb3 = rgb2*(1+alpha*rl1)**(-2)
                     rgc3 = rgc2*(1+alpha*rl1)**(-2)
                     rga4 = rga1*(1+alpha*rl2)**(-2)
                     rgb4 = rgb1*(1+alpha*rl2)**(-2)
                     rgc4 = rgc1*(1+alpha*rl2)**(-2)
                     call INTERSECTION_METRIC_77(rga(k1),rgb(k1),rgc(k1),  &
                          rga3,rgb3,rgc3)
                     call INTERSECTION_METRIC_77(rga(k2),rgb(k2),rgc(k2),  &
                          rga4,rgb4,rgc4)
                     
                     ra(3*(i-1)+j) = -1.
                     if(iae(i,j) .gt. 0) then
                        ii = iae(i,j)
                        do l=1,3
                           if(iae(ii,l) == i) ra(3*(ii-1)+l) = -1.
                        enddo
                     endif
                  endif
               enddo
            enddo
         enddo
!     end of MESH GRADIENT CONTROL

      elseif(icontrol == 2) then
!     MESH GRADIENT CONTROL H-shock
         beta = 1.5
         do is =1,2
            do k=1,3*nelem
               ra(k) = 1.
            enddo
            do i=1,nelem
               do j=1,3
                  if(ra(3*(i-1)+j) .gt. 0) then
                     j1 = mod(j,3) + 1
                     k1 = lnd(i,j)
                     k2 = lnd(i,j1)
                     x1 = x(k1)
                     y1 = y(k1)
                     x2 = x(k2)
                     y2 = y(k2)
                     rga1 = rga(k1)
                     rgb1 = rgb(k1)
                     rgc1 = rgc(k1)
                     rga2 = rga(k2)
                     rgb2 = rgb(k2)
                     rgc2 = rgc(k2)
                     rnoem = ((x2-x1)**2 + (y2-y1)**2)**0.5
                     xx = (x2 -x1)/rnorm
                     yy = (y2 -y1)/rnorm
                     rl1 = (rga1*xx**2 + 2*rgb1*xx*yy +  &
                          rgc1*yy**2)**0.5
                     rl2 = (rga2*xx**2 + 2*rgb2*xx*yy +  &
                          rgc2*yy**2)**0.5
                     rl = ((rga1+rga2)/2*(x1-x2)**2 +   &
                          (rgb1+rgb2)*(x1-x2)*(y1-y2) +  &
                          (rgc1+rgc2)/2*(y1-y2)**2)**0.5
                     if(rl2 .gt. rl1) then
                        cc = (rl2/rl1)**(1./rl)
                        if(cc .gt. beta) then
                           eta = (beta/cc)**rl
                           print *,'eta =',eta
                           rga(k2) = rga(k2)/eta/eta
                           rgb(k2) = rgb(k2)/eta/eta
                           rgc(k2) = rgc(k2)/eta/eta
                        endif
                     endif
                     ra(3*(i-1)+j) = -1.
                     if(iae(i,j) .gt. 0) then
                        ii = iae(i,j)
                        do l=1,3
                           if(iae(ii,l) == i) ra(3*(ii-1)+l) = -1.
                        enddo
                     endif
                  endif
               enddo
            enddo
         enddo
!     end of MESH GRADIENT CONTROL
      endif


!     smoothing
      print*,'################### smoothing', 2
      do 1500 ic = 1,2
!      do 1500 ic = 1,2
         do 1400 i=1,nelem
            ra(i) = (rga(lnd(i,1)) + rga(lnd(i,2)) + rga(lnd(i,3)))/3
            ra(i+melem) =   &
                 (rgb(lnd(i,1)) + rgb(lnd(i,2)) + rgb(lnd(i,3)))/3
            ra(i+2*melem) =   &
                 (rgc(lnd(i,1)) + rgc(lnd(i,2)) + rgc(lnd(i,3)))/3
 1400     enddo
         
         do 1410 i=1,npoin
            rga(i) = 0.
            rgb(i) = 0.
            rgc(i) = 0.
            area(i) = 0.
 1410     enddo
      
         do 1420 i=1,nelem
            xc = (x(lnd(i,1))+x(lnd(i,2))+x(lnd(i,3)))/3
            yc = (y(lnd(i,1))+y(lnd(i,2))+y(lnd(i,3)))/3
            do 1430 j=1,3
               k = lnd(i,j)
!               rl = ((xc-x(k))**2 + (yc-y(k))**2)**0.5
               rga(k) = rga(k) + ra(i)*tria(i)
               rgb(k) = rgb(k) + ra(i+melem)*tria(i)
               rgc(k) = rgc(k) + ra(i+2*melem)*tria(i)
               area(k) = area(k) + tria(i)
 1430        enddo
 1420     enddo
         
         do 1440 i=1,npoin
            rga(i) = rga(i)/area(i)
            rgb(i) = rgb(i)/area(i)
            rgc(i) = rgc(i)/area(i)
         if(i == itest) print *,'$3',rga(i),rgb(i),rgc(i)
 1440     enddo
 1500  enddo


!     ... no plotting (too big file)
       if(npoin .gt. 5000) return

      imt = 18
      open(imt,file='mt',status='unknown')
      do 63 i=1,npoin
!     for plotting of the elipses
         if(ibp(i,1) .gt. 0) then
            rga(i) = (rga(i) + rga(ibp(i,1)))/2
            rga(ibp(i,1)) = rga(i)
            rgb(i) = (rgb(i) + rgb(ibp(i,1)))/2
            rgb(ibp(i,1)) = rgb(i)
            rgc(i) = (rgc(i) + rgc(ibp(i,1)))/2
            rgc(ibp(i,1)) = rgc(i)
         endif
         
         if(mod(i, 8) == 1)   &
              call ELIPS_77(imt,rga(i),rgb(i),rgc(i),x(i),y(i) )
!         if(i == itest) print *,'$$',rga(i),rgb(i),rgc(i)
 63   enddo
!      stop
! 1313 continue


      return
    end subroutine ERROR1_77

    subroutine ELIPS_77(ieli,a,b,c,xi,yi)
      !print*,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
      num = 100
      do i=0,num
         t = 1.*i/num*6.283185307
         x = cos(t)
         y = sin(t)
         rnorm = (x*x*a + 2*x*y*b + y*y*c)**0.5
         xs = x/rnorm + xi
         ys = y/rnorm + yi
         write(ieli,*) xs,ys, i, xi,yi, ieli
         !write(*,*) xs,ys, i, xi,yi, ieli
      enddo
      write(ieli, *)'## '
      write(ieli, *)
      return
    end subroutine ELIPS_77
         


    subroutine INTERSECTION_METRIC_77(ra1,rb1,rc1,ra2,rb2,rc2)

      real*8 a1,a2,b1,b2,c1,c2,p1,p2,p3,p4,d,rmat1(2),rmat2(2),  &
           rlam1,rlam2,pi1,pi2,pi3,pi4
      
      a1 = ra1
      b1 = rb1
      c1 = rc1
      a2 = ra2
      b2 = rb2
      c2 = rc2
      
      p3 = 1.D+00
      p4 = 1.D+00

      d = (a1*c2 - a2*c1)**2 -4.D+00*(a1*b2-b1*a2)*(b1*c2-c1*b2)
      if((a1*b2-b1*a2) .ne. 0 .and. d .ge. 0.) then
         p2 = (a2*c1-a1*c2 - d**0.5)/2.D+00/(a1*b2-b1*a2)
      elseif(a2*c1-a1*c2 .ne. 0. .and. d .ge. 0. ) then
         p2 = (b1*c2-c1*b2)/(a2*c1-a1*c2)
      else
!         print *,'zero polynom'
         if(a1 .le. a2) then
            ra1 = ra2
            rb1 = rb2
            rc1 = rc2
         endif
         return
      endif

      if(p2*a2 + b2 .ne. 0.) then
         p1 = -(b2*p2 + c2)/(p2*a2 + b2)
      else
         if(a1 .le. a2) then
            ra1 = ra2
            rb1 = rb2
            rc1 = rc2
         endif
         return
      endif
         
      rmat1(1) = p1**2*a1+2*p1*p3*b1+p3**2*c1
      rmat112 = p2*p1*a1+p2*p3*b1+p4*p1*b1+p4*p3*c1
      rmat121 = p2*p1*a1+p2*p3*b1+p4*p1*b1+p4*p3*c1
      rmat1(2) = p2**2*a1+2*p2*p4*b1+p4**2*c1

      rmat2(1) = p1**2*a2+2*p1*p3*b2+p3**2*c2
      rmat212 = p2*p1*a2+p2*p3*b2+p4*p1*b2+p4*p3*c2
      rmat221 = p2*p1*a2+p2*p3*b2+p4*p1*b2+p4*p3*c2
      rmat2(2) = p2**2*a2+2*p2*p4*b2+p4**2*c2

      rlam1 = max(rmat1(1),rmat2(1) )
      rlam2 = max(rmat1(2),rmat2(2) )

      if(rlam1 .le. 0. .or. rlam2 .le. 0) then
         print *,'nonpositive eigenvalues'
         stop
      endif

      d = p1*p4 - p2*p3
      if(d == 0.) then
         print *,'zero deter.'
         stop
      endif
      pi1 = p4/d
      pi2 = -p2/d
      pi3 = -p3/d
      pi4 = p1/d

      ra1 = pi1**2*rlam1+pi3**2*rlam2
      rb1 = pi1*rlam1*pi2+pi3*rlam2*pi4
      rc1 = pi2**2*rlam1+pi4**2*rlam2

      return
    end subroutine INTERSECTION_METRIC_77


    subroutine ERROR2_77(melem,nelem,mpoin,npoin,maxdeg,dx,dy,area,  &
         tria,noit,  &
         x,y,lnd,iae,wp,icyc,ra,rga,rgb,rgc,surface,  &
         ibp)
      dimension lnd(melem,3),iae(melem,3),x(mpoin),y(mpoin),  &
           ra(melem*3),wp(mpoin,5),dx(melem),dy(melem),  &
           icyc(mpoin,maxdeg),rga(mpoin),rgb(mpoin),rgc(mpoin),  &
           area(mpoin),tria(melem),ibp(mpoin,2)
      real*8 dxi,dyi,rmeas,x1,x2,x3,y1,y2,y3,w1,w2,w3,areai,rl
      real*8 a,b,c,disc,rlam1,rlam2,x11,x12,x21,x22,y12,y11,y21,y22,  &
           t11,t12,t21,t22,z11,z12,z21,z22,rdet

      rkappa = 1.4



!     dx the first derivative dw/dx
!     dy                      dw/dy
!     rga  the second der.    d2w/dx dx
!     rgb                     d2w/dx dy
!     rg!                     d2w/dy dy       

      do 1600 i1=1,nelem
         do 1610 ji=1,3
            j1 = mod(ji,3) + 1
            if(iae(i1,ji) .lt. 0 .and. iae(i1,j1) .lt. 0) then
               print *,'$$$',icha,i1,ji,j1
               print *,lnd(i1,1),lnd(i1,2),lnd(i1,3)
               print *,x(lnd(i1,1)),y(lnd(i1,1))
               print *,x(lnd(i1,2)),y(lnd(i1,2))
               print *,x(lnd(i1,3)),y(lnd(i1,3))
            endif
 1610    enddo
 1600 enddo



      do 10 i=1,nelem
         x1 = x(lnd(i,1))
         y1 = y(lnd(i,1))
         w1 = wp(lnd(i,1),1)
         x2 = x(lnd(i,2))
         y2 = y(lnd(i,2))
         w2 = wp(lnd(i,2),1)
         x3 = x(lnd(i,3))
         y3 = y(lnd(i,3))
         w3 = wp(lnd(i,3),1)
         rmeas = x1*y2 + x2*y3 + x3*y1 - x1*y3  -x2*y1 - x3*y2
         tria(i) = rmeas/2
         if( rmeas .le. 0) then
            print *,'bad triangle in ERROR'
            print *,i,lnd(i,1),lnd(i,2),lnd(i,3)
            print *,x1,y1,rmeas
            print *,x2,y2
            print *,x3,y3
            stop
         endif
         dx(i) = (y2-y1)*(w1+w2) + (y3-y2)*(w2+w3) + (y1-y3)*(w1+w3) 
         dy(i) = (x1-x2)*(w1+w2) + (x2-x3)*(w2+w3) + (x3-x1)*(w1+w3) 
         dx(i) = dx(i) /rmeas
         dy(i) = dy(i) /rmeas

         if( i == -10) then
            write(*,'(4i5 )')i,lnd(i,1),lnd(i,2),lnd(i,3)
            write(*,'(4e14.6 )')x1,y1,w1,rmeas
            write(*,'(3e14.6 )')x2,y2,w2
            write(*,'(3e14.6 )')x3,y3,w3
            write(*,'(3e14.6 )')(x1-x2),(x2-x3),(x3-x1)
            write(*,'(3e14.6 )')(w1+w2),(w2+w3),(w1+w3)
            write(*,'(3e14.6 )')(y2-y1),(y3-y2),(y1-y3)
            write(*,'(2e14.6 )')dx(i),dy(i)

         endif

 10   enddo

      do 20 i=1,npoin
         areai = 0.
         rga(i) = 0.
         rgb(i) = 0.
         rgc(i) = 0.
         rgdi = 0.
         if(icyc(i,1) .gt. 0) then
            len = icyc(i,1) 
         else
            len = -icyc(i,1)  - 1
         endif
!         if(i == -236)    print *,x(i),y(i),wp(i,1),i,len
         x3 = x(i)
         y3 = y(i)
         w3 = wp(i,1)
         do 30 j=1,len
            j1 = mod(j, abs(icyc(i,1))) + 1
            x1 = x(icyc(i,j+1))
            y1 = y(icyc(i,j+1))
            w1 = wp(icyc(i,j+1),1)
            x2 = x(icyc(i,j1+1))
            y2 = y(icyc(i,j1+1))
            w2 = wp(icyc(i,j1+1),1)
            rmeas = x1*y2 + x2*y3 + x3*y1 - x1*y3  -x2*y1 - x3*y2
            areai = areai + rmeas/6
            dxi = (y2-y1)*(w1+w2) + (y3-y2)*(w2+w3) + (y1-y3)*(w1+w3) 
            dyi = (x1-x2)*(w1+w2) + (x2-x3)*(w2+w3) + (x3-x1)*(w1+w3) 
            dxfi = (y1-y2)/rmeas
            dyfi = (x2-x1)/rmeas

            rga(i) = rga(i) - dxi/2*dxfi
            rgb(i) = rgb(i) - dyi/2*dxfi
            rgc(i) = rgc(i) - dyi/2*dyfi
            rgdi = rgdi - dxi/2*dyfi
            
!            itest = 13
!            if(i == itest ) then
!               write(*,'(a1,2i5,4e14.6)')
!     *              ' ',i,icyc(i,j+1), rga(i),rgb(i),rgdi,rgc(i)
!               write(*,'(3i5,4e14.6)') i,icyc(i,j+1),icyc(i,j1+1),
!     *              w3,w1,w2,x1
!               write(*,'(5e14.6)') x1,y1,w1,dxi,dyi
!               write(*,'(5e14.6)') x2,y2,w2,dxfi,dyfi
!               write(*,'(4e14.6)') x3,y3,w3
!               write(*,'(3e14.6)') x(i),y(i),wp(i,1)
!            endif

            if(icyc(i,1) .lt. 0 .and. ibp(i,1) == 0 ) then
!     for the boundary points me must add the boundary sides
               if(j == 1) then
                  rga(i) = rga(i) + dxi/rmeas*(y1-y3)/2
                  rgb(i) = rgb(i) + dyi/rmeas*(y1-y3)/2
                  rgc(i) = rgc(i) + dyi/rmeas*(x3-x1)/2
                  rgdi = rgdi + dxi/rmeas*(x3-x1)/2
!                  if(i == itest) write(*,'(a1,2i5,4e14.6)')
!     *                 'f',i,icyc(i,j+1), rga(i),rgc(i),dxi,dyi
               endif
               if( j == len) then
                  rga(i) = rga(i) + dxi/rmeas*(y3-y2)/2
                  rgb(i) = rgb(i) + dyi/rmeas*(y3-y2)/2
                  rgc(i) = rgc(i) + dyi/rmeas*(x2-x3)/2
                  rgdi = rgdi + dxi/rmeas*(x2-x3)/2

!                  if(i == itest) write(*,'(a1,2i5,4e14.6)')
!     *                 'l',i,icyc(i,j+1), rga(i),rgc(i),dxi/rmeas,
!     *                 dyi/rmeas
               endif
            endif
 30      enddo
         area(i) = areai
         
 20   enddo

!     for periodic boundary points
      if( AMA%xper(1,1) .gt. 0 .or. AMA%xper(1,2) .gt. 0) then
         do 18 i1=1,npoin
            i2 = ibp(i1,1)
            if( i2 .gt. 0 ) then
               if(rga(i1) .ne. rga(i2) .or. rgb(i1) .ne. rgb(i2) .or.   &
                    rgc(i1) .ne. rgc(i2) .or. area(i1) .ne. area(i2) )  &
                    then
                  rga(i1) = rga(i1) + rga(i2)
                  rgb(i1) = rgb(i1) + rgb(i2)
                  rgc(i1) = rgc(i1) + rgc(i2)
                  area(i1) = area(i1) + area(i2)
                  rga(i2) = rga(i1) 
                  rgb(i2) = rgb(i1) 
                  rgc(i2) = rgc(i1) 
                  area(i2) = area(i1)
               endif
            endif
 18      enddo
      endif

      do 15 i=1,npoin
         rga(i) = rga(i)/area(i)
         rgb(i) = rgb(i)/area(i)
         rgc(i) = rgc(i)/area(i)

!cx         u = x(i)
!cx         v = y(i)

!cx         rga(i) = 4*(-1+3*u*u)/(1+u*u)**3
!cx         rgb(i) = 0.
!cx         rgc(i) = 0.
!         write(10,'(5e14.6)') x(i),y(i),rga(i),rgb(i),rgc(i)
 15   enddo

      rmaxder = 0.
!     improvement for boundary points

      do 100 i=1,npoin
         if(icyc(i,1) .lt. 0 .and. ibp(i,1) == 0) then
            len = abs(icyc(i,1))
            if( len .ge. 3) then
               ipoc = 0
               rga(i) = 0.
               rgb(i) = 0.
               rgc(i) = 0.
               rminlenght = 1E+38
               do 110 j=2,len-1
                  if(icyc(icyc(i,j+1),1) .gt. 0) then
!     this point isn't boundary
!                     rga(i) = rga(i) + rga(icyc(i,j+1))
!                     rgb(i) = rgb(i) + rgb(icyc(i,j+1))
!                     rgc(i) = rgc(i) + rgc(icyc(i,j+1))
!                     ipoc = ipoc + 1
                     xlen = ((x(i)-x(icyc(i,j+1)))**2 +  &
                          (y(i) -y(icyc(i,j+1)))**2 )**0.5
                     if(xlen .lt. rminlenght) then
                        rga(i) = rga(icyc(i,j+1))
                        rgb(i) = rgb(icyc(i,j+1))
                        rgc(i) = rgc(icyc(i,j+1))
                        rminlenght = xlen
                     endif
                     if( i == -170) then
                        print *,'**',j,icyc(i,j+1),rgc(icyc(i,j+1)),  &
                             xlen,rminlenght 
                     endif
                  endif
 110           enddo
!     we use the value from the nearest point
               if(ipoc .gt. 0 ) then
                  rga(i) = rga(i)/ipoc
                  rgb(i) = rgb(i)/ipoc
                  rgc(i) = rgc(i)/ipoc
               else
!     in this case, we let the second derivations = 0 !!!!

               endif
            elseif( len == 2) then
               i1 = icyc(i,2)
               i2 = icyc(i,3)
               do 120 ie =1,nelem
                  do 130 je =1,3
                     je1 = mod(je,3)+1
                     if(lnd(ie,je) == i2 .and.   &
                          lnd(ie,je1) == i1) then
                        je2 = mod(je1,3) + 1
                        rga(i) = rga(lnd(ie,je2))
                        rgb(i) = rgb(lnd(ie,je2))
                        rgc(i) = rgc(lnd(ie,je2))
                        goto 140
                     endif
 130              enddo
 120           enddo
 140           continue
            else
               print *,'ERROR len < 2 !!!!'
            endif
         endif
!         rga(i) = 0.
!         rgb(i) = 0.
!         rgc(i) = 0.
!         if(icyc(i,1) .lt. 0 .and. ibp(i,1) == 0) rgc(i) = 100.
!         write(11,'(5e14.6)') x(i),y(i),rga(i),rgb(i),rgc(i)
         rmaxder = max(rmaxder,abs(rga(i)),abs(rgb(i)),abs(rgc(i)) )
 100  enddo
!      stop

!      return
!     the area of ideal triangle
      rnorm = 1.29903810567 *AMA%numel/surface
!     Gamm channel
!      epsilon1 = 7500.
!      epsilon2 = 7500.

!     NACA
!      epsilon1 = 10000.
!      epsilon2 = 250.

!     profile
!      epsilon1 = 250.
!      epsilon2 = 10000.

!     automatic generation
!      p = 0.95
!      epsilon1 = rmaxder*3/rnorm*10.
!      epsilon1 = min(10000.,epsilon1)
!      epsilon1 = 10000.
!      epsilon1 = max (50.,epsilon1)

      rmaxderteor = AMA%epsilon1/3*rnorm
      der = rmaxderteor
!      der = min(rmaxderteor,rmaxder)

      epsilon2p = der/10*(AMA%epsilon1 + 1.)*(1.-AMA%p)/  &
           (AMA%p*(AMA%epsilon1 + 1.)-1.)

      epsilon2 = max (1.,epsilon2p)

      write(AMA%ifig1,*)'                    p =',p,'    ityp =',AMA%ityp
      write(AMA%ifig1,*)'   number of elements =',AMA%numel
      write(AMA%ifig1,*)'rmaxder,rmaxderteor,c =',rmaxder,rmaxderteor,rnorm
      write(AMA%ifig1,*)'epsilon1 and epsilon2 = ',AMA%epsilon1,epsilon2,'(',  &
           epsilon2p,')'

      ipr = -1
!     we compute the elements of matrix M from the second derivations 

      rmax1 = 0.
      rmax2 = 0.
      rmax3 = 0.
      
      do 200 i=1,npoin
!      do 200 i=1,1

         if(rga(i) .gt. rmax1) then
            rmax1 = rga(i)
            imax1 = i
         endif
         if(rgb(i) .gt. rmax2) then
            rmax2 = rgb(i)
            imax2 = i
         endif
         if(rgc(i) .gt. rmax3) then
            rmax3 = rgc(i)
            imax3 = i
         endif

         if(abs(rgb(i)) .lt. 1E-05) then
!     the diagonal metrix

            eps1 = AMA%epsilon1/(epsilon2+max(rga(i),rgc(i)))

            rgb(i) = 0.
            rga(i) = (1. + eps1*abs(rga(i)))*rnorm
            rgc(i) = (1. + eps1*abs(rgc(i)))*rnorm
         else
!     the eigenvalues rlam1 , rlam2
            a = rga(i)
            b = rgb(i)
            c = rgc(i)

            disc = ((a - c)*(a - c) + 4*b*b)**0.5
            rlam1 = (a + c + disc)/2
            rlam2 = (a + c - disc)/2
            if(abs(rlam1) .lt. abs(rlam2) ) then
               rlampom = rlam1
               rlam1 = rlam2
               rlam2 = rlampom
            endif

            x11 = b
            x21 = -(a-rlam1)
            rl1 = (x11*x11+x21*x21)**0.5
            x11 = x11/rl1
            x21 = x21/rl1

            x12 = b
            x22 = -(a-rlam2)
            rl2 = (x12*x12+x22*x22)**0.5
            x12 = x12/rl2
            x22 = x22/rl2

            rdet = x11*x22 - x21*x12

!     inverse matrix
            y11 = x22/rdet
            y12 = -x12/rdet
            y21 = -x21/rdet
            y22 = x11/rdet

            z11 = abs(rlam1)*y11
            z12 = abs(rlam1)*y12
            z21 = abs(rlam2)*y21
            z22 = abs(rlam2)*y22
            
            t11 = x11*z11 + x12*z21
            t12 = x11*z12 + x12*z22
            t21 = x21*z11 + x22*z21
            t22 = x21*z12 + x22*z22

            if(i == ipr) then
               print *,'***********************',i,'    *************'
               print *,a,b,c,rlam1,rlam2
               print *,a-rlam1,b
               print *,b,c-rlam1
               print *,'            '
               print *,'----------------------'
               write(*,'(6e12.4)')  &
                    x11,x12,rlam1,0,y11,y12
               write(*,'(6e12.4)')  &
                    x21,x22,0,rlam2,y21,y22
               print *,'----------------------'
               write(*,'(4e12.4)')  &
                    z11,z12,t11,t12
               write(*,'(4e12.4)')  &
                    z21,z22,t21,t22
               print *,'----------------------'
            endif

            eps1 = AMA%epsilon1/(epsilon2+max(t11,t22))

            rga(i) = (1. + eps1*t11)*rnorm 
            rgb(i) = (eps1*(t12+t21)/2)*rnorm
            rgc(i) = (1. + eps1*t22)*rnorm

            if(rga(i)*rgc(i) .le. rgb(i)*rgb(i) ) then
               print *,'ERROR in metric, nonpositive matrix'
               write(*,'(i5,3e14.6)') i,rga(i),rgb(i),rgc(i)
               print *,a,b,c
               print *,'----------------------'
               write(*,'(6e12.4)')  &
                    x11,x12,rlam1,0,y11,y12
               write(*,'(6e12.4)')  &
                    x21,x22,0,rlam2,y21,y22
               print *,'----------------------'
               write(*,'(4e12.4)')  &
                    z11,z12,t11,t12
               write(*,'(4e12.4)')  &
                    z21,z22,t21,t22
               print *,'----------------------'
               stop
            endif

         endif

 200  enddo

!      print *,'*************************'
!      print *,imax1,rmax1,x(imax1),y(imax1)
!      print *,imax2,rmax2,x(imax2),y(imax2)
!      print *,imax3,rmax3,x(imax3),y(imax3)


!     smoothing
      do 500 ic = 1,2
         do 400 i=1,nelem
            ra(i) = (rga(lnd(i,1)) + rga(lnd(i,2)) + rga(lnd(i,3)))/3
            ra(i+melem) =   &
                 (rgb(lnd(i,1)) + rgb(lnd(i,2)) + rgb(lnd(i,3)))/3
            ra(i+2*melem) =   &
                 (rgc(lnd(i,1)) + rgc(lnd(i,2)) + rgc(lnd(i,3)))/3
 400     enddo
         
         do 410 i=1,npoin
            rga(i) = 0.
            rgb(i) = 0.
            rgc(i) = 0.
            area(i) = 0.
 410     enddo
      
         do 420 i=1,nelem
!            xc = (x(lnd(i,1))+x(lnd(i,2))+x(lnd(i,3)))/3
!            yc = (y(lnd(i,1))+y(lnd(i,2))+y(lnd(i,3)))/3
            do 430 j=1,3
               k = lnd(i,j)
!               rl = ((xc-x(k))**2 + (yc-y(k))**2)**0.5
               rga(k) = rga(k) + ra(i)*tria(i)
               rgb(k) = rgb(k) + ra(i+melem)*tria(i)
               rgc(k) = rgc(k) + ra(i+2*melem)*tria(i)
               area(k) = area(k) + tria(i)
 430        enddo
 420     enddo
         
         do 440 i=1,npoin
            rga(i) = rga(i)/area(i)
            rgb(i) = rgb(i)/area(i)
            rgc(i) = rgc(i)/area(i)
            if( i .lt. -20) then
               write(*,'(i5,4e14.6)') i,rga(i),rgb(i),rgc(i),wp(i,1)
            endif
 440     enddo
 500  enddo
     
      do 63 i=1,npoin
         if(ibp(i,1) .gt. 0) then
            rga(i) = (rga(i) + rga(ibp(i,1)))/2
            rga(ibp(i,1)) = rga(i)
            rgb(i) = (rgb(i) + rgb(ibp(i,1)))/2
            rgb(ibp(i,1)) = rgb(i)
            rgc(i) = (rgc(i) + rgc(ibp(i,1)))/2
            rgc(ibp(i,1)) = rgc(i)
         endif
!         write(12,'(5e14.6)') x(i),y(i),rga(i),rgb(i),rgc(i)
 63   enddo
!      stop
      return
    end subroutine ERROR2_77


    subroutine SIXTRIANG_77(xp,yp,x1,y1,x2,y2,x3,y3,x4,  &
         y4,x5,y5,x6,y6,sum)
      !         print *, 'uvnitr sixtriang '     
! vypocet normy fi v prostoru H(0,1), fi ma hodnotu 1 v xp,yp
! a 0 v xi,yi, je linearni
      
      call TRIANG_77(xp,yp,x1,y1,x2,y2,s1)
      call TRIANG_77(xp,yp,x2,y2,x3,y3,s2)
      call TRIANG_77(xp,yp,x3,y3,x4,y4,s3)
      call TRIANG_77(xp,yp,x4,y4,x5,y5,s4)
      call TRIANG_77(xp,yp,x5,y5,x6,y6,s5)
      call TRIANG_77(xp,yp,x6,y6,x1,y1,s6)

      sum=s1+s2+s3+s4+s5+s6
      return
    end subroutine SIXTRIANG_77


    subroutine  TRIANG_77(xp,yp,xa,ya,xb,yb,s)
      
      !         print *, 'uvnitr triang  '     
      ! vypocet integralu pro jednotlivy trojuhelnicek, ktery je 
      ! casti sestiuhelniku
      pom=(xp-xb)*(ya-yb)-(xa-xb)*(yp-yb)
      T=(ya-yb)/pom
      R=(xb-xa)/pom
      
      a=sqrt((xa-xp)*(xa-xp)+(ya-yp)*(ya-yp))
      b=sqrt((xb-xp)*(xb-xp)+(yb-yp)*(yb-yp))
      c=sqrt((xa-xb)*(xa-xb)+(ya-yb)*(ya-yb))

      cb=(c*c+b*b-a*a)/(2*c)
      v=sqrt(b*b-cb*cb)
      s=0.
!       s=v*c/12.
!        we have the space H(0,1), we use a seminorm, 
!             which is equivalent
      s=s+v*c*0.5*(T*T+R*R)
      return
    end subroutine TRIANG_77


    subroutine SEEK_BOUN_77(melem,nelem,mpoin,npoin,ipoint,nbp,xb,yb,  &
         lnd,iae,ibb,ibpoin,x,y, icrack)
      dimension lnd(melem,3),iae(melem,3),xb(ipoint),yb(ipoint),  &
           ibb(mpoin,3),x(mpoin),y(mpoin),ibpoin(0:nbp)
      
      do i=1,npoin
         ibb(i,1) = -1
         ibb(i,2) = 0
      enddo
      

!     no smooth boundary
      if(nbp == 0) return
     
      do i=1,nelem
         do j=1,3
            if(iae(i,j) .lt. 0) then
               j1 = mod(j,3) + 1
               ip = lnd(i,j)
!               if(ip == 112) print *,'###'
               xp = x(ip)
               yp = y(ip)
               rlen0 = (xp - x(lnd(i,j1)))*(xp - x(lnd(i,j1))) +   &
                    (yp - y(lnd(i,j1)))*(yp - y(lnd(i,j1)))
               rl0 = 1.E+25*rlen0
               do k=1,nbp
                  do l=ibpoin(k-1)+1,ibpoin(k)
                     rlen = (xp - xb(l))*(xp - xb(l)) +  &
                          (yp - yb(l))*(yp - yb(l))
                     if(rlen .lt. rl0) then
                        if(icrack .gt. 0 .and. xp .gt. -0.01) then
                           yc = (y(lnd(i,1))+y(lnd(i,2))+  &
                                y(lnd(i,3)))/3
                           if(( yc .gt. 0 .and. l .lt. ibpoin(k)/2)   &
                                .or. ( yc .lt. 0 .and.   &
                                l .gt. ibpoin(k)/2) )then
                              rl0 = rlen
                              il0 = l
                              kl0 = k
                           endif
                        else
                           rl0 = rlen
                           il0 = l
                           kl0 = k
                        endif
                     endif
                     
!                     if(ip  == 2) then
!                        write(91,'(4i5,9e10.2)')
!     *                       ip,l,il0,kl0,x(ip),y(ip), xb(l),yb(l), 
!     *                       rlen,rl0,rlen0,rl0/rlen0
!                     endif

                  enddo
!               if(ip == 2) 
!                  write(*,'(a3,3i5,6e12.4)')
!     *                 '###',ip,il0,kl0,x(ip),y(ip),rl0,rlen0,rl0/rlen0
               enddo
               if(rl0/rlen0 .lt. 1E-03) then
!     the node in 'profile' found
                  ibb(ip,1) = il0
                  ibb(ip,2) = kl0
               endif
            endif
         enddo
      enddo

       do i=1,npoin
         if(ibb(i,1) .gt.0 ) then

!            write(21,'(2e12.4,3i8)') 
!     *        xb(ibb(i,1)),yb(ibb(i,1)),ibb(i,1),ibb(i,2),i

!     corection of the boundary            
            x(i) = xb(ibb(i,1))
            y(i) = yb(ibb(i,1))
         endif       
      enddo

      return
    end subroutine SEEK_BOUN_77

    subroutine COINS_77(melem,nelem,mpoin,npoin,lnd,iae,  &
         nbelm,mbelm,lbn,itc)
      dimension lnd(melem,3),iae(melem,3),lbn(mbelm,2),itc(mbelm)
      
      !     if some triangle has a two times iae(i,j) < 0 the me make a 
      !     bascule such that this case doesn't appear

      do 10 i=1,nelem
         do 20 j=1,3
            j1 = mod(j,3) + 1
            j2 = mod(j1,3) + 1
            if( iae(i,j1) .lt. 0 .and. iae(i,j2) .lt. 0 ) then
               ii = iae(i,j)
               do 30 jjj=1,3
                  if(iae(ii,jjj) == i) then
                     jj = jjj
                  endif
 30            enddo
               jj1 = mod(jj,3) + 1
               jj2 = mod(jj1,3) + 1
               k1 = lnd(i,j)
               k2 = lnd(i,j1)
               k3 = lnd(i,j2)
               k4 = lnd(ii,jj2)
               if(k1 .ne. lnd(ii,jj1) .or. k2 .ne. lnd(ii,jj) ) then
                  print *,'ERROR in COINS 1'
               endif
               ia1 = iae(ii,jj1)
               if(ia1 .gt. 0) then
                  do 40 k=1,3
                     if(iae(ia1,k) == ii) ja1 = k
 40               enddo
               endif
               ia2 = iae(ii,jj2)
               if(ia2 .gt. 0) then
                  do 50 k=1,3
                     if(iae(ia2,k) == ii) ja2 = k
 50               enddo
               endif
               lnd(i,j1) = k4
               lnd(ii,jj1) = k3
               iae(i,j1) = ii
               iae(ii,jj) = -2
               iae(i,j) = ia1
               iae(ii,jj1) = i
               if(ia1 .gt. 0) iae(ia1,ja1) = i

               do ib=1,nbelm
                  if(itc(ib) == i .and.   &
                       lbn(ib,1) == k2 .and.lbn(ib,2) ==k3)then
                     itc(ib) = ii
                  endif
               enddo

            endif
 20      enddo
 10   enddo
      return
    end subroutine COINS_77

    
    subroutine ADJAC_77(melem,nelem,mpoin,npoin,lnd,iae,maxdeg,icyc)
      dimension lnd(melem,3),iae(melem,3),icyc(mpoin,maxdeg)
      
      !     the array icyc has here only the local meaning
      do i=1,npoin
         icyc(i,1) = 0
      enddo
      do i=1,nelem
         do j=1,3
            k = lnd(i,j)
            icyc(k,1) = icyc(k,1) + 1
            if(icyc(k,1) .gt. maxdeg) then
               print *,'Bad dimension in ADJAC'
               print *,'maxdeg < icyc(k,1)',maxdeg,icyc(k,1)
               print *,k
               stop
            endif
            icyc(k,icyc(k,1)+1) = i
         enddo
      enddo
      
      do 10 i=1,nelem
         iae(i,1) = -2
         iae(i,2) = -2
         iae(i,3) = -2
 10   enddo
      do 20 i=1,nelem
         do 30 j=1,3
            if(iae(i,j) == -2)then
               j1 = mod(j,3) +1
               do 40 il=1,icyc(lnd(i,j),1)
                  ii = icyc(lnd(i,j),il+1)
                  if( ii .ne. i) then
                     do 50 jj=1,3
                        jj1 = mod(jj,3) +1
                        if( (lnd(i,j) == lnd(ii,jj1) ).and.  &
                             (lnd(i,j1) == lnd(ii,jj))) then
                           iae(i,j) = ii
                           iae(ii,jj) = i
                           goto 60
                        endif
 50                  enddo
                  endif
 40            enddo
 60            continue
            endif
 30      enddo
 20   enddo

!      print *,'The seeking of neighbourhouds done'
      return
    end subroutine ADJAC_77

    subroutine CYKLE_BOUND_77(melem,nelem,mpoin,npoin,maxdeg,x,y,  &
         lnd,iae,icyc)
      dimension x(mpoin),y(mpoin),lnd(melem,3),  &
           icyc(mpoin,maxdeg),iae(melem,3)
      

      do 10 i=1,nelem
         do 20 j=1,3
            if(iae(i,j) .le. 0) then
               j0 = j
               j1 = mod(j0,3) + 1
               ip = lnd(i,j0)
               if(icyc(ip,1) .ne. -1) then
                  print *,'ERROR in CYKLE_BOUND'
                  print *,ip,icyc(ip,1),icyc(ip,2),icyc(ip,3)
                  print *,x(ip),y(ip)
                  return
               endif
               icyc(ip,2) = lnd(i,j1)
               ii = i
!     seeking the following point
 1000          ilen = abs(icyc(ip,1))
               itest = 0
               do 30 l=1,3
                  if(lnd(ii,l) == ip) then
                     itest = 1
                     j0 = l
                     goto 40
                  endif
 30            enddo
 40            continue
               j1 = mod(j0 ,3) + 1
               j2 = mod(j1,3) + 1
               icyc(ip,ilen+2) = lnd(ii,j2)
               icyc(ip,1) = icyc(ip,1) - 1
               if(abs(icyc(ip,1)) +1 .gt. maxdeg) then
                  print *,'The lenght of cykles > maxdeg = ',maxdeg
                  stop
               endif

               ii1 = iae(ii,j2)
               if(ii1 .gt. 0) then
                  ii = ii1
                  goto 1000
               endif
            endif
 20      enddo
 10   enddo

      
!      test of angles
      ipoc = 0
      if(ipoc == 0) return

      rmax = 10.
      iost = 12134
      open(iost,file='ost',status='unknown')
      do 2000 i=1,nelem
         do 2010 j=1,3
            j1 = mod(j,3) + 1
            j2 = mod(j1,3) + 1
            x1 = x(lnd(i,j))
            y1 = y(lnd(i,j))
            x2 = x(lnd(i,j1))
            y2 = y(lnd(i,j1))
            x3 = x(lnd(i,j2))
            y3 = y(lnd(i,j2))
            u1 = x2 - x1
            v1 = y2 - y1
            u2 = x3 - x1
            v2 = y3 - y1
            rcos = (u1*u2+v1*v2)/((u1*u1+v1*v1)*(u2*u2+v2*v2))**0.5
            if (rcos .lt. -0.2) then
               write(iost,'(2e14.6)')x1,y1
               write(iost,'(2e14.6)')x2,y2
               write(iost,'(2e14.6)')x3,y3
               write(iost,'(2e14.6)')x1,y1
               write(iost,'(x)')

               ipoc = ipoc + 1
               if(rcos .lt. rmax) then
                  rmax = rcos
                  imax = i
                  jmax = j
               endif
            endif
 2010    enddo
 2000 enddo
      write( *,'(a12,3i5,2e14.6)')   &
           'the number = ',ipoc,imax,jmax,rmax,acos(rmax)
      print *,'nelem = ',nelem
      return
    end subroutine CYKLE_BOUND_77

      subroutine CYKLE_77(melem,nelem,mpoin,npoin,maxdeg,x,y,  &
           lnd,iae,icyc,lbn,nbelm,mbelm,ibp,itc,ibc)
      dimension x(mpoin),y(mpoin),lnd(melem,3),lbn(mbelm,2),  &
           icyc(mpoin,maxdeg),  &
           iae(melem,3),ibp(mpoin,2),ibc(mbelm),itc(mbelm)

      

      do 5 i=1,npoin
         icyc(i,1) = 0
 5    enddo

      do 10 k=1,nbelm
         icyc(lbn(k,1),1) = -1
         icyc(lbn(k,2),1) = -1
 10   enddo

!  ... seeking the  cyclus around each point exept boundaries

      do 30 i=1,nelem
         do 40 j=1,3
            ip = lnd(i,j)
            if(icyc(ip,1) == 0) then
!     initialisation
               j1 = mod(j,3) + 1
               k1 = lnd(i,j1)
               ii = i
               jj2 = mod(j1,3) + 1
               icyc(ip,1) = icyc(ip,1) + 1
               icyc(ip,icyc(ip,1) + 1) = k1
!     repetition
 45            icyc(ip,1) = icyc(ip,1) + 1
               icyc(ip,icyc(ip,1) + 1) = lnd(ii,jj2)
               ii = iae(ii,jj2)
               jj = 0
               do jjj=1,3
                  if(lnd(ii,jjj) == ip ) jj = jjj
               enddo
               if(jj == 0) then
                  print *,'Triangle for node',ip,'doesn''t found'
                  print *,'@@',ip, ii,lnd(ii,1),lnd(ii,2),lnd(ii,3)
                  print *,'!!',icyc(ip,1),icyc(ip,2),icyc(ip,3),  &
                       icyc(ip,4),icyc(ip,5),icyc(ip,6)
                  print *,x(ip),y(ip)
                  print *,x(icyc(ip,2)),y(icyc(ip,2))
                  print *,x(icyc(ip,3)),y(icyc(ip,3))
                  stop
               endif
               jj1 = mod(jj,3) + 1
               jj2 = mod(jj1,3) + 1
               if( lnd(ii,jj2) .ne. k1) goto 45
            endif
 40      enddo
 30   enddo

!      print *,'The seeking of cykles done'

      if(AMA%xper(1,1) .gt. 0. .or. AMA%xper(1,2) .gt. 0.) then
         ibeg1 = 0
         ibeg2 = 0
         num1 = 0
         num2 = 0
         do i=1,nbelm
            if(ibeg1 == 0 .and. ibc(i) == AMA%iper(1,1)) ibeg1 = i
            if(ibeg2 == 0 .and. ibc(i) == AMA%iper(1,2)) ibeg2 = i
            if(ibc(i) == AMA%iper(1,1)) num1 = num1 + 1
            if(ibc(i) == AMA%iper(1,2)) num2 = num2 + 1
            if(ibc(i) == AMA%iper(1,1) .or. ibc(i) == AMA%iper(1,2)) then
               iel = itc(i)
               do j=1,3
                  j1 = mod(j,3) +1
                  if(lnd(iel,j) == lbn(i,1) .and.   &
                       lnd(iel,j1) == lbn(i,2) ) then
                     iae(iel,j) = -1
                     goto 2005
                  endif
               enddo
 2005          continue
            endif
         enddo
         if(num1 .ne. num2) then
            print *,'On the periodic boundary is not the same number ',  &
                 'of points'
            stop
         endif
         if(num1 == 0) then
            print *,'On the periodic boundary is zero points'
            stop
         endif
         do i=1,num1
            k1 = lbn(ibeg1 -1 + i,1)
            k2 = lbn(ibeg2 + num1 - i,2)
            ibp(k1,1) = k2
            ibp(k2,1) = k1
            ibp(k1,2) = ibp(k1,2) + 1
            ibp(k2,2) = ibp(k2,2) + 1
            if( i == num1) then
               k1 = lbn(ibeg1 -1 + i,2)
               k2 = lbn(ibeg2 + num1 - i,1)
               ibp(k1,1) = k2
               ibp(k2,1) = k1
               ibp(k1,2) = ibp(k1,2) + 1
               ibp(k2,2) = ibp(k2,2) + 1
            endif
         enddo

!         print *,'The periodic boundary done'
      endif

!     the second pair of periodic boundaty component
      if(AMA%xper(2,1) .gt. 0. .or. AMA%xper(2,2) .gt. 0.) then
         ibeg1 = 0
         ibeg2 = 0
         num1 = 0
         num2 = 0
         do i=1,nbelm
            if(ibeg1 == 0 .and. ibc(i) == AMA%iper(2,1)) ibeg1 = i
            if(ibeg2 == 0 .and. ibc(i) == AMA%iper(2,2)) ibeg2 = i
            if(ibc(i) == AMA%iper(2,1)) num1 = num1 + 1
            if(ibc(i) == AMA%iper(2,2)) num2 = num2 + 1
            if(ibc(i) == AMA%iper(2,1) .or. ibc(i) == AMA%iper(2,2)) then
               iel = itc(i)
               do j=1,3
                  j1 = mod(j,3) +1
                  if(lnd(iel,j) == lbn(i,1) .and.   &
                       lnd(iel,j1) == lbn(i,2) ) then
                     iae(iel,j) = -1
                     goto 2006
                  endif
               enddo
 2006          continue
            endif
         enddo
         if(num1 .ne. num2) then
            print *,'On the periodic boundary is not the same number ',  &
                 'of points'
            stop
         endif
         if(num1 == 0) then
            print *,'On the periodic boundary is zero points'
            stop
         endif
!     ...ibp(k,2) = 3   ==> point is a periodic for both part of boundary
         do i=1,num1
            k1 = lbn(ibeg1 -1 + i,1)
            k2 = lbn(ibeg2 + num1 - i,2)
            ibp(k1,1) = k2
            ibp(k2,1) = k1
            ibp(k1,2) = ibp(k1,2) + 2
            ibp(k2,2) = ibp(k2,2) + 2
!            print *,' .',k1,ibp(k1,1),k2,ibp(k2,1)
            if( i == num1) then
               k1 = lbn(ibeg1 -1 + i,2)
               k2 = lbn(ibeg2 + num1 - i,1)
               ibp(k1,1) = k2
               ibp(k2,1) = k1
               ibp(k1,2) = ibp(k1,2) + 2
               ibp(k2,2) = ibp(k2,2) + 2
            endif
         enddo

!         print *,'The periodic boundary done'
      endif

      return
    end subroutine CYKLE_77

      subroutine TEST_77(melem,nelem,mpoin,npoin,maxdeg,x,y,  &
           lnd,iae,icyc,lbn,nbelm,mbelm,ibp,itc)
      dimension x(mpoin),y(mpoin),lnd(melem,3),lbn(mbelm,2),  &
           icyc(mpoin,maxdeg),itc(mbelm),  &
           iae(melem,3),ibp(mpoin,2)

      
      write(AMA%ifig1,*)'Control of input/output files'

!                 do i=1,nelem
!                     write(*,'(7i5)') i,lnd(I,1),lnd(I,2),lnd(I,3),
!     *                    iae(i,1),iae(i,2),iae(i,3)
!                  enddo
!                  print*,'********************'
 
!     ... test of input datas
      do 600 i=1,nelem
         x1 = x(lnd(i,1))
         y1 = y(lnd(i,1))
         x2 = x(lnd(i,2))
         y2 = y(lnd(i,2))
         x3 = x(lnd(i,3))
         y3 = y(lnd(i,3))
         det = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
         if(det .le. 0. ) then
            print *,'Bad triangle in GT'
            print *,x1,y1,lnd(i,1),det
            print *,x2,y2,lnd(i,2)
            print *,x3,y3,lnd(i,3)
            stop
         endif
         do 610 j=1,3

            j1 = mod(j,3) + 1
            if(lnd(i,j) == lnd(i,j1) ) then
               print *,'ERROR in GT:'
               print *,i,'-th triangle contains the same nodes'
               write(*,'(3i5)') lnd(i,1),lnd(i,2),lnd(i,3)
               stop
            endif
            if(iae(i,j) == iae(i,j1) ) then
               print *,'ERROR in GT:'
               print *,'code found bad neighbours for', i,'-th triangle'
               write(*,'(3i5)') iae(i,1),iae(i,2),iae(i,3)
               stop
            endif
            if(iae(i,j) .gt. 0) then
               ii = iae(i,j)
               jj = 0
               do 620 jjj=1,3
                  if(iae(ii,jjj) == i) jj = jjj
 620           enddo
               if( jj == 0) then
                  print *, 'ERROR in GT:'
                  print *, 'the neighbourhouds are bad'
                  write(*,'(7i5)') i,lnd(i,1),lnd(i,2),lnd(i,3),  &
                       iae(i,1),iae(i,2),iae(i,3)
                  write(*,'(7i5)') ii,lnd(ii,1),lnd(ii,2),lnd(ii,3),  &
                       iae(ii,1),iae(ii,2),iae(ii,3)
                  stop
               endif
               jj1 = mod(jj,3) + 1
               if(lnd(i,j) .ne. lnd(ii,jj1) .or.  &
                    lnd(i,j1) .ne. lnd(ii,jj) ) then
                  print *,'ERROR in GT: '
                  print *, 'the neighbourhouds points are bad'
                  write(*,'(7i5)') i,lnd(i,1),lnd(i,2),lnd(i,3),  &
                       iae(i,1),iae(i,2),iae(i,3)
                  write(*,'(7i5)') ii,lnd(ii,1),lnd(ii,2),lnd(ii,3),  &
                       iae(ii,1),iae(ii,2),iae(ii,3)
                  stop
               endif
            endif
 610     enddo
 600  enddo
      
!     checking of boundary
      do k=1,nbelm
         it = itc(k)
         k1 = lbn(k,1)
         k2 = lbn(k,2)

         if(it .gt. nelem) then
            print *,'ERROR in GT:'
            print *,'too high number itc'
            print *,k,it,k1,k2,nelem
            stop
         endif

         if((lnd(it,1) == k1 .and. lnd(it,2) == k2 .and.   &
              iae(it,1) .lt. 0) .or.  &
              (lnd(it,2) == k1 .and. lnd(it,3) == k2 .and.   &
              iae(it,2) .lt. 0) .or.  &
              (lnd(it,3) == k1 .and. lnd(it,1) == k2 .and.   &
              iae(it,3) .lt. 0))  then
!     OK
         else
            print *,'ERROR in GT:'
            print *,k,'-th boundary segment is bad'
            print *,'lbn(k,1) and lbn(k,2) don''t coresponds with',  &
                 'the field lnd(i,j)'
            print *,k,k1,k2
            stop
         endif
      enddo


!     checking of the  cykles
      do 900 ipp =1,npoin-1
         if( icyc(ipp,1) .gt. 0 ) then
            len = icyc(ipp,1)
         else
            len = abs( icyc(ipp,1)) - 1
         endif
         do 910 ic = 1,len
            ic1 = mod(ic,abs(icyc(ipp,1)) ) + 1
            ip1 = icyc(ipp,ic+1)
            ip2 = icyc(ipp,ic1+1)
            itest = 0
            do 920 jpp=1,nelem
               if(( lnd(jpp,1) == ipp .and.  &
                    lnd(jpp,2) == ip1 .and.  &
                    lnd(jpp,3) == ip2 ) .or.  &
                    ( lnd(jpp,2) == ipp .and.  &
                    lnd(jpp,3) == ip1 .and.  &
                    lnd(jpp,1) == ip2 ) .or.  &
                    ( lnd(jpp,3) == ipp .and.  &
                    lnd(jpp,1) == ip1 .and.  &
                    lnd(jpp,2) == ip2 ) ) then
                  itest = 1
                  goto 922
               endif
 920        enddo
 922        continue
            if(itest == 0)then
               print *,'ERROR in CHECk in GT'
               print *,ipp,ip1,ip2,ic,ic1
               it = ipp
               print *,x(it),y(it)
               write(*,'(10i5)') it,icyc(it,1),  &
                    icyc(it,2),  &
                    icyc(it,3),icyc(it,4),icyc(it,5),  &
                    icyc(it,6),icyc(it,7),  &
                    icyc(it,8),icyc(it,9)
               do 912 itt=1,abs(icyc(it,1))
                  print *,x(icyc(it,itt+1)),  &
                       y(icyc(it,itt+1)),icyc(it,itt+1)
 912           enddo
               do iu=1,nelem
                  if(lnd(iu,1) == it .or. lnd(iu,2) == it  &
                       .or. lnd(iu,3) == it) then
                     write(*,'(4i5)') iu,lnd(iu,1),lnd(iu,2),  &
                          lnd(iu,3)
                  endif
               enddo
               stop
            endif
 910     enddo
 900  enddo

      if(AMA%xper(1,1) .gt. 0 .or. AMA%xper(1,2) .gt. 0) then
         ipoc = 0
         do i=1,npoin
            if(ibp(i,1) .gt. 0 .and. ibp(i,2) == 1) then
               ipoc = ipoc + 1
               if(i .ne. ibp(ibp(i,1),1) ) then
                  print *,'ERROR in GT:'
                  print *,'The periodic nodes in triang are ',  &
                       'not periodical-1:'
                  print *,i,ibp(i,1)
                  print *,ibp(i,1),ibp(ibp(i,1),1)
                  print *,x(i),y(i)
                  print *,x(ibp(i,1)),y(ibp(i,1))
                  print *,x(ibp(ibp(i,1),1)),y(ibp(ibp(i,1),1))
                  stop
               endif

               xl = ( abs(x(i) - x(ibp(i,1))) ) - AMA%xper(1,1)
               yl = ( abs(y(i) - y(ibp(i,1))) ) - AMA%xper(1,2)
               if( xl .gt. 1E-05 .or. yl .gt. 1E-05) then
                  print *,'ERROR in GT, bad periodical boundary-1'
                  print *,i,ibp(i,1),ibp(ibp(i,1),1)
                  print *,x(i),y(i)
                  print *,x(ibp(i,1)),y(ibp(i,1))
                  print *,AMA%xper(1,1),AMA%xper(1,2),xl,yl
                  stop
               endif
            endif
         enddo
         if(ipoc == 0) then
            print *,'ERROR in GT:'
            print *,'non periodic point!!!-1'
            stop
         endif
      endif

      if(AMA%xper(2,1) .gt. 0 .or. AMA%xper(2,2) .gt. 0) then
         ipoc = 0
         do i=1,npoin
            if(ibp(i,1) .gt. 0 .and. ibp(i,2) == 2) then
               ipoc = ipoc + 1
               if(i .ne. ibp(ibp(i,1),1) ) then
                  print *,'ERROR in GT:'
                  print *,'The periodic nodes in triang are ',  &
                       'not periodical-2:'
                  print *,i,ibp(i,1)
                  print *,ibp(i,1),ibp(ibp(i,1),1)
                  print *,x(i),y(i)
                  print *,x(ibp(i,1)),y(ibp(i,1))
                  print *,x(ibp(ibp(i,1),1)),y(ibp(ibp(i,1),1))
                  stop
               endif

               xl = ( abs(x(i) - x(ibp(i,1))) ) - AMA%xper(2,1)
               yl = ( abs(y(i) - y(ibp(i,1))) ) - AMA%xper(2,2)
               if( xl .gt. 1E-05 .or. yl .gt. 1E-05) then
                  print *,'ERROR in GT, bad periodical boundary-2'
                  print *,i,ibp(i,1),ibp(ibp(i,1),1)
                  print *,x(i),y(i)
                  print *,x(ibp(i,1)),y(ibp(i,1))
                  print *,AMA%xper(2,1),AMA%xper(2,2),xl,yl
                  stop
               endif
            endif
         enddo
         if(ipoc == 0) then
            print *,'ERROR in GT:'
            print *,'non periodic point!!!-2'
            stop
         endif
      endif
      write(AMA%ifig1,*)'Input/output file is OK'
      return
    end subroutine TEST_77


    subroutine CYKLE_REP_77(melem,nelem,mpoin,npoin,maxdeg,x,y,  &
         lnd,iae,icyc,lbn,nbelm,mbelm,ibp)
      dimension x(mpoin),y(mpoin),lnd(melem,3),lbn(mbelm,2),  &
           icyc(mpoin,maxdeg),iae(melem,3),ibp(mpoin,2)

      do 5 i=1,npoin
         icyc(i,1) = 0
 5    enddo

      do 10 k=1,nbelm
         icyc(lbn(k,1),1) = -1
         icyc(lbn(k,2),1) = -1
 10   enddo

!  ... seeking the  cyclus around each point exept boundaries

      do 30 i=1,nelem
         do 40 j=1,3
            ip = lnd(i,j)
            if(icyc(ip,1) == 0) then
!     initialisation
               j1 = mod(j,3) + 1
               k1 = lnd(i,j1)
               ii = i
               jj2 = mod(j1,3) + 1
               icyc(ip,1) = icyc(ip,1) + 1
               icyc(ip,icyc(ip,1) + 1) = k1
!     repetition
 45            icyc(ip,1) = icyc(ip,1) + 1
               if(icyc(ip,1) +1 .gt. maxdeg) then
                  print *,'The lenght of cykles > maxdeg = ', maxdeg
                  stop
               endif
               icyc(ip,icyc(ip,1) + 1) = lnd(ii,jj2)
               ii = iae(ii,jj2)
               jj = 0
               do jjj=1,3
                  if(lnd(ii,jjj) == ip ) jj = jjj
               enddo
               if(jj == 0) then
                  print *,'Triangle for node',ip,'doesn''t found'
                  stop
               endif
               jj1 = mod(jj,3) + 1
               jj2 = mod(jj1,3) + 1
               if( lnd(ii,jj2) .ne. k1) goto 45
            endif
 40      enddo
 30   enddo
      return
    end subroutine CYKLE_REP_77

    subroutine QUA2_77(x1,y1,w1,ra1,rb1,rc1,x2,y2,w2,ra2,rb2,rc2,  &
         x3,y3,w3,ra3,rb3,rc3,err1,err2,ice,ibo )
      dimension rra(3),rrb(3),rrc(3)
      
      err1 = 1.
      sqrt3 = 1.7320508075

      rl1 = (((ra1+ra2)*(x1-x2)*(x1-x2) + 2*(rb1+rb2)*(x1-x2)*(y1-y2) +  &
           (rc1+rc2)*(y1-y2)*(y1-y2) )/2 )**0.5
      rl2 = (((ra2+ra3)*(x2-x3)*(x2-x3) + 2*(rb2+rb3)*(x2-x3)*(y2-y3) +  &
           (rc2+rc3)*(y2-y3)*(y2-y3) )/2 )**0.5
      rl3 = (((ra1+ra3)*(x1-x3)*(x1-x3) + 2*(rb1+rb3)*(x1-x3)*(y1-y3) +  &
           (rc1+rc3)*(y1-y3)*(y1-y3) )/2 )**0.5

      err2 = (sqrt3-rl1)**2 + (sqrt3-rl2)**2 + (sqrt3-rl3)**2

      if(ibo == 1) then
         err2 = err2 + (sqrt3-rl1)**2
      elseif(ibo == 2) then
         err2 = err2 + (sqrt3-rl2)**2 
      elseif(ibo == 3) then
         err2 = err2 + (sqrt3-rl3)**2 
      endif
      

      AMA%glmax = max(AMA%glmax, rl1**2,rl2**2,rl3**2)
      AMA%glmin = min(AMA%glmin, rl1**2,rl2**2,rl3**2)

      if(ice == 13) then
         write(*,'(4e14.6)') rl1,rl2,rl3,err2
         write(*,'(5e14.6)') x1,y1,ra1,rb1,rc1
         write(*,'(5e14.6)') x2,y2,ra2,rb2,rc2
         write(*,'(5e14.6)') x3,y3,ra3,rb3,rc3
      endif
      return
    end subroutine QUA2_77

    subroutine QUALITY_77(ndim,melem,nelem,mpoin,npoin,x,y,lnd,iae,iter,  &
         err,rminerr,imin, rmaxrez,ra,wp,rga,rgb,rgc)
      dimension x(mpoin),y(mpoin),lnd(melem,3),iae(nelem,3),  &
           ra(melem*3),wp(mpoin,ndim+1),  &
           rga(mpoin),rgb(mpoin),rgc(mpoin)
      
      ice = 0
      err = 0.
      rminerr = 1.
      AMA%errrez = 0.
      rmaxrez = 0.
      AMA%glmax = 0.
      AMA%glmin = 1000.
      ipobound = 0
      ibo = 0
      do 10 i=1,nelem
         i1 = lnd(i,1)
         i2 = lnd(i,2)
         i3 = lnd(i,3)
         x1 = x(i1)
         y1 = y(i1)
         x2 = x(i2)
         y2 = y(i2)
         x3 = x(i3)
         y3 = y(i3)
         w1 = wp(i1,1)
         w2 = wp(i2,1)
         w3 = wp(i3,1)
         ra1 = rga(i1)
         rb1 = rgb(i1)
         rc1 = rgc(i1)
         ra2 = rga(i2)
         rb2 = rgb(i2)
         rc2 = rgc(i2)
         ra3 = rga(i3)
         rb3 = rgb(i3)
         rc3 = rgc(i3)
         

         do j=1,3
            if( iae(i,j) .lt. -1) then
               ibo = j
               ipobound = ipobound + 1 
            endif
         enddo

         call QUA2_77(x1,y1,w1,ra1,rb1,rc1,x2,y2,w2,ra2,rb2,rc2,  &
              x3,y3,w3,ra3,rb3,rc3,err1,err2,ice,ibo)

         ice = 0
         err = err + err1
         AMA%errrez = AMA%errrez + err2
         if(rminerr .gt. err1 ) then
            rminerr = err1
            imin = i
         endif
         rmaxrez = max(rmaxrez,err2)

         if( i == -1 )   &
              write(*,'(a12,i5,6es12.4)') 'ama-restA',i,AMA%errrez,err2,rmaxrez

 10   enddo
      err = err/nelem
      AMA%errrez = AMA%errrez/(nelem*3 + ipobound)


      if(AMA%ifig .ge. 0) then
         call PLOT1_77(melem,nelem,mpoin,npoin,x,y,lnd)
         AMA%ifig = AMA%ifig + 1
      endif

      return
    end subroutine QUALITY_77

    subroutine PLOT1_77(melem,nelem,mpoin,npoin,x,y,lnd)
      dimension x(mpoin),y(mpoin),lnd(melem,3)

      
      character*10 meshfile
      character*1 num
      character*2 nux
      character*3 nuy
      character*4 nuz
      character*5 nuw

      meshfile(1:10) = 'mesh      '
      if( AMA%ifig .lt. 10 ) then
         write( num, '(i1)' ) AMA%ifig
         meshfile(5:5) =  num
      elseif( AMA%ifig .lt. 100 ) then
         write( nux, '(i2)' ) AMA%ifig
         meshfile(5:6) =  nux
      elseif( AMA%ifig .lt. 1000 ) then
         write( nuy, '(i3)' ) AMA%ifig
         meshfile(5:7) =  nuy
      elseif( AMA%ifig .lt. 10000 ) then
         write( nuz, '(i4)' ) AMA%ifig
         meshfile(5:8) =  nuz
      elseif( AMA%ifig .lt. 100000 ) then
         write( nuw, '(i5)' ) AMA%ifig
         meshfile(5:9) =  nuw
      else
         print *, ' Number of gnu-files has to be less then 99999'
         stop
      endif

!      ... plot of a new mesh

      igra = 53
      open(igra, file=meshfile,status='UNKNOWN')
      
      do i=1,nelem
         write(igra,'(2e14.6)') x(lnd(i,1)),y(lnd(i,1))
         write(igra,'(2e14.6)') x(lnd(i,2)),y(lnd(i,2))
         write(igra,'(2e14.6)') x(lnd(i,3)),y(lnd(i,3))
         write(igra,'(2e14.6)') x(lnd(i,1)),y(lnd(i,1))
         write(igra,'(20x)')
      enddo
      close(igra)

    end subroutine PLOT1_77

    subroutine PLOT_77(melem,nelem,mpoin,npoin,x,y,lnd)
      dimension x(mpoin),y(mpoin),lnd(melem,3)

      

!      ... plot of a new mesh

      igra = 53
      open(igra, file='mesh',status='UNKNOWN')
      
      do i=1,nelem
         write(igra,'(2e14.6)') x(lnd(i,1)),y(lnd(i,1))
         write(igra,'(2e14.6)') x(lnd(i,2)),y(lnd(i,2))
         write(igra,'(2e14.6)') x(lnd(i,3)),y(lnd(i,3))
         write(igra,'(2e14.6)') x(lnd(i,1)),y(lnd(i,1))
         write(igra,'(20x)')
      enddo
      if(AMA%xper(1,1) .gt. 0 .or. AMA%xper(1,2) .gt. 0) then
         do  i=1,nelem
            write(igra,'(2e14.6)') x(lnd(i,1))+AMA%xper(1,1),y(lnd(i,1))+AMA%xper(1,2)
            write(igra,'(2e14.6)') x(lnd(i,2))+AMA%xper(1,1),y(lnd(i,2))+AMA%xper(1,2)
            write(igra,'(2e14.6)') x(lnd(i,3))+AMA%xper(1,1),y(lnd(i,3))+AMA%xper(1,2)
            write(igra,'(2e14.6)') x(lnd(i,1))+AMA%xper(1,1),y(lnd(i,1))+AMA%xper(1,2)
            write(igra,'(20x)')
         enddo
      endif

      if(AMA%xper(2,1) .gt. 0 .or. AMA%xper(2,2) .gt. 0) then
         do  i=1,nelem
            write(igra,'(2e14.6)') x(lnd(i,1))+AMA%xper(2,1),y(lnd(i,1))+AMA%xper(2,2)
            write(igra,'(2e14.6)') x(lnd(i,2))+AMA%xper(2,1),y(lnd(i,2))+AMA%xper(2,2)
            write(igra,'(2e14.6)') x(lnd(i,3))+AMA%xper(2,1),y(lnd(i,3))+AMA%xper(2,2)
            write(igra,'(2e14.6)') x(lnd(i,1))+AMA%xper(2,1),y(lnd(i,1))+AMA%xper(2,2)
            write(igra,'(20x)')
         enddo
      endif
      return
    end subroutine PLOT_77


    subroutine WriteMetrix_77(mpoin, npoin, melem, nelem, lnd,   &
           x, y, rga, rgb, rgc)
      integer :: mpoin, npoin, melem, nelem
      dimension x(mpoin),y(mpoin), rga(mpoin),rgb(mpoin),rgc(mpoin),  &
           lnd(melem, 3)

      print*,'####,  Writing file "fort.',AMA%ifig2
!      do i=1,npoin
!         r_max = max(rga(i), rgc(i) )
!         r_min = min(rga(i), rgc(i) )
!
!         if(r_max > 1) 
!     *        write(AMA%ifig2, *) x(i), y(i), r_max, r_max/r_min
!
!      enddo


      do i=1,nelem

         do j=1,4
            ip = lnd(i, mod(j,3)+1)
            r_max = max(rga(ip), rgc(ip) )
            r_min = min(rga(ip), rgc(ip) )

            write(AMA%ifig2, *)   &
                 x(ip), y(ip), r_max, r_max/r_min
         enddo
         write(AMA%ifig2, '(x)') 
         write(AMA%ifig2, '(x)') 
         write(AMA%ifig2, '(x)') 
         
         
      enddo
      AMA%ifig2 = AMA%ifig2 + 1
      
    end subroutine WriteMetrix_77


    subroutine MOVING_77(ndim, melem,nelem,mpoin,npoin,  &
         maxdeg,x,y,lnd,iae,  &
         ra,icyc,wp,rga,rgb,rgc,noit,ipoint,nbp,xb,yb,ibb,ibpoin,  &
         ichag,mbelm,nbelm,ibp)
      dimension lnd(melem,3),iae(melem,3),x(mpoin),y(mpoin),  &
           ra(mpoin),wp(mpoin,ndim),  &
           rga(mpoin),rgb(mpoin),rgc(mpoin),  &
           icyc(mpoin,maxdeg),  &
           xb(ipoint),yb(ipoint),ibp(mpoin,2),ibb(mpoin,3),ibpoin(0:nbp)
      real xloc(25),yloc(25)
      integer iloc(25)
      real*8 a,b,c,d,det0,det1,det2,det3,xs,ys,yq,xq,  &
           rlen0,rlen1,rlen2,ss,rcos


      do 10 i=1,npoin
         if(ibb(i,3) == -1) goto 10
         if(icyc(i,1) .gt. 0) then
            if(icyc(i,1) .gt. 25) then
               print *,'Error in dimension in MOVING',icyc(i,1)
               return
            endif
            xp = x(i)
            yp = y(i)
            xq = 0.
            yq = 0.
            qualityE = 0.
            do 20 j=1,icyc(i,1)
               j1 = mod(j,icyc(i,1) ) +1
               ii1 = icyc(i,j+1)
               ii2 = icyc(i,j1+1)
               xi = x(ii1)
               yi = y(ii1)
               xi1 = x(ii2)
               yi1 = y(ii2)
               a = (rga(ii1) + rga(ii2))/2
               b = (rgb(ii1) + rgb(ii2))/2
               c = (rgc(ii1) + rgc(ii2))/2
               d = ( (a*c-b*b)/3)**0.5
               if( d .le. 0) then
                  print *,'ERROR in matrix in moving'
                  print *,a,b,c
                  stop
               endif

               xs = xi + ((d-b)*(xi1-xi) - c*(yi1-yi))/2/d
               ys = yi + ( a*(xi1-xi) + (d+b)*(yi1-yi))/2/d

               xloc(j) = xs
               yloc(j) = ys

               xq = xq + xloc(j)
               yq = yq + yloc(j)
 20         enddo

            xq = xq/icyc(i,1)
            yq = yq/icyc(i,1)

!     now we check the the triangles xi,xi1,xq are reals
            numit = 10
            k=numit
 1240       continue
            ierr = 0
            xnew = xp + 1.*k/numit*(xq - xp) 
            ynew = yp + 1.*k/numit*(yq - yp) 
            do 1230 j=1,icyc(i,1)
               j1 = mod(j,icyc(i,1))+1
               ii1 = icyc(i,j+1)
               ii2 = icyc(i,j1+1)
               xi = x(ii1)
               yi = y(ii1)
               xi1 = x(ii2)
               yi1 = y(ii2)
!     check the orientation and bad conditionality
               detiii = xi*(yi1-ynew) + xi1*(ynew-yi) + xnew*(yi-yi1)
               reps = AMA%pos*( ((xi1-xi)**2 + (yi1-yi)**2) +  &
                    ((xi-xnew)**2 + (yi-ynew)**2) +  &
                    ((xi1-xnew)**2 + (yi1-ynew)**2) )
               if (detiii .le. reps) then 
                  ierr = 1
                  goto 1231
               endif
               call POS1TEST_77(xi,yi,xi1,yi1,xnew,ynew,itet)
               if(itet == 1) then
                  ierr = 1
                  goto 1231
               endif

1230        enddo
1231        continue

            if(ierr == 1 .and. k .ge. 2) then
               k = k - 1
               goto 1240
            endif
            if( k == 1) then
               xq = xp
               yq = yp
            else
               xq = xnew
               yq = ynew
            endif

            qualityE = 0.
            qualityrez = 1E+35
            numit = 8
            ice = 0
            xold = xp
            yold = yp

            do 40 k=0,numit
               quality1 = 0.
               quality1rez = 0.
               xnew = xp + 1.*k/numit*(xq - xp) 
               ynew = yp + 1.*k/numit*(yq - yp) 
               do 30 j=1,icyc(i,1)
                  j1 = mod(j,icyc(i,1))+1
                  ii1 = icyc(i,j+1)
                  ii2 = icyc(i,j1+1)
                  xi = x(ii1)
                  yi = y(ii1)
                  xi1 = x(ii2)
                  yi1 = y(ii2)
!     check the orientation and bad conditionality
                  detiii = xi*(yi1-ynew) + xi1*(ynew-yi) + xnew*(yi-yi1)
                  reps = AMA%pos*( ((xi1-xi)**2 + (yi1-yi)**2) +  &
                    ((xi-xnew)**2 + (yi-ynew)**2) +  &
                    ((xi1-xnew)**2 + (yi1-ynew)**2) )
                  if (detiii .le. reps) then 
!                     print *,'VIOLATION of positivity - algorith ERROR'
                     goto 45
                  endif
                  call POS1TEST_77(xi,yi,xi1,yi1,xnew,ynew,itet)
                  if(itet == 1) then
                     goto 45
                  endif

                  acc = ACCUTE_77(xi,yi,xi1,yi1,xnew,ynew, 1.)

                  rl1 = (0.5*((rga(i)+rga(ii1))*(xi-xnew)*(xi-xnew)+  &
                       2*(rgb(i)+rgb(ii1))*(xi-xnew)*(yi-ynew)+  &
                       (rgc(i)+rgc(ii1))*(yi-ynew)*(yi-ynew)))**0.5

                  quality1rez = quality1rez + (rl1 - 1.732050807)**2
!     ... NEW ACCUTE
!                  quality1rez = quality1rez 
!     *                 +(rl1 - 1.732050807)**2*acc

30             enddo

               if(quality1rez .ge. qualityrez ) goto 45
               
               xold = xnew
               yold = ynew
               qualityrez = quality1rez
40          enddo
45          continue

            x(i) = xold
            y(i) = yold
         endif
10    enddo
         
      !     for boundary points, the movement, only if the point is on the
      !     straight segment
      !     or on the boundary given by the field [xb(i),yb(i)] i=1,ipoint

      do 100 i=1,npoin
         
         if(ibb(i,3) == -1) goto 100
         if(ibb(i,1) .gt. 0) then 
            ikk1 = ibb(i,2)
            if(ibb(i,1) == ibpoin(ikk1-1)+1  .or.  &
                 ibb(i,1) == ibpoin(ikk1)) then
               !               NO moving of final or initial node of profiles
               goto 21
            endif
         endif
         
         if(icyc(i,1) .lt. 0 ) then
            len = abs(icyc(i,1))
            x1 = x(i)
            y1 = y(i)
            x2 = x(icyc(i,2))
            y2 = y(icyc(i,2))
            x0 = x(icyc(i,len+1))
            y0 = y(icyc(i,len+1))
            rlen1 = ((x2-x1)**2 + (y2-y1)**2)**0.5
            rlen2 = ((x0-x1)**2 + (y0-y1)**2)**0.5
            ss = ((x2-x1)*(x0-x1) + (y2-y1)*(y0-y1))
            rcos = ss/rlen1/rlen2
            !            if(i == 2) print *,'%%%%',rcos
            if( rcos .gt. -0.9999 .or. (ibb(i,1) .gt. 0 .and. &
                 !            if( rcos .gt. -0.9866 .or. (ibb(i,1) .gt. 0 .and.  &
                 (ibb(icyc(i,2),1) .gt. 0 .and.   &
                 ibb(icyc(i,len+1),1) .lt. 0 ) .or.   &
                 (ibb(icyc(i,2),1) .lt. 0 .and.   &
                 ibb(icyc(i,len+1),1) .gt. 0 ) )) then
!     too sharp angle or. begin or end of noclosed profile,
!     movement forbidden
               goto 100
            elseif ( ibb(i,1) .lt. 0  ) then
!     this points are on straight segments
!     the maximal value of numit can not be greater than 9 !!!!!!
               numit = 9
               ice = 0
               xold = x1
               yold = y1
               do 141 k=-numit,0
                  l = numit + k + 1
                  xloc(l) = x1 - 1.*k/(numit+1)*(x0 - x1) 
                  yloc(l) = y1 - 1.*k/(numit+1)*(y0 - y1) 
                  iloc(l) = -1
141            enddo
               do 146 k=1,numit
                  l = numit + k + 1
                  xloc(l) = x1 + 1.*k/(numit+1)*(x2 - x1) 
                  yloc(l) = y1 + 1.*k/(numit+1)*(y2 - y1) 
                  iloc(l) = -1
146            enddo
               ilo = 2*numit + 1
            else
               !     nodes are on curved part of boundary
               i2 = icyc(i,2)
               i0 = icyc(i,abs(icyc(i,1))+1)
               
               if(ibb(i,1) .lt. 0 .or. ibb(i0,1) .lt. 0   &
                    .or. ibb(i2,1) .lt.0) then
                  print *,'LOGICAL error in MOVING, ibb'
                  print *,i,i0,i2,ibb(i,1),ibb(i0,1),ibb(i2,1)
                  print *,x(i),y(i)  !,xb(ibb(i,1)),yb(ibb(i,1))
                  print *,x(i0),y(i0)  !,xb(ibb(i0,1)),yb(ibb(i0,1))
                  print *,x(i2),y(i2) !,xb(ibb(i2,1)),yb(ibb(i2,1))
                  stop
               endif
               
               il0 = ibb(i0,1)
               il1 = ibb(i,1)
               il2 = ibb(i2,1)
               ilk = ibb(i,2)
               
               if(ilk.ne. ibb(i0,2) .or. ilk .ne. ibb(i2,2) ) then
                  print *,'error jkol'
                  print *,x(i),y(i),i,ibb(i,1),ibb(i,2),ilk
                  print *,x(i0),y(i0),i0,ibb(i0,1),ibb(i0,2)
                  print *,x(i2),y(i2),i2,ibb(i2,1),ibb(i2,2)
                  stop
               endif
         
!     reordering in the case of end point
               if(il0 == ibpoin(ilk-1)+1 .and. il1 .gt. il2) then
                  il0 = ibpoin(ilk)
                  if(abs(xb(ibpoin(ilk-1)+1)-xb(ibpoin(ilk)) )   &
                       .gt. 1E-05 .or.  &
                       abs(yb(ibpoin(ilk-1)+1)-yb(ibpoin(ilk)) )   &
                       .gt. 1E-05 )then
                     print *,'error in profile, the end point is not'
                     print *,'the first, or bad reordering'
                     print *,ilk,(ibpoin(ilk-1)+1),ibpoin(ilk)
                     print *,'1',il0,il1,il2
                     stop
                  endif
               endif

               if(il0 == ibpoin(ilk) .and. il2 .gt. il1) then
                  il0 = ibpoin(ilk-1)+1
                  if(abs(xb(ibpoin(ilk-1)+1)-xb(ibpoin(ilk)) )   &
                       .gt. 1E-05 .or.  &
                       abs(yb(ibpoin(ilk-1)+1)-yb(ibpoin(ilk)) )   &
                       .gt. 1E-05 )then
                     print *,'error in profile, the end point is not'
                     print *,'the first, or bad reordering'
                     print *,'2',il0,il1,il2
                     stop
                  endif
               endif
               
               if(il2 == ibpoin(ilk-1)+1 .and. il1 .gt. il0) then
                  il2 = ibpoin(ilk)
                  if(abs(xb(ibpoin(ilk-1)+1)-xb(ibpoin(ilk)) )   &
                       .gt. 1E-05 .or.  &
                       abs(yb(ibpoin(ilk-1)+1)-yb(ibpoin(ilk)) )   &
                       .gt. 1E-05 )then
                     print *,'error in profile, the end point is not'
                     print *,'the first, or bad reordering'
                     print *,'3',il0,il1,il2
                     stop
                  endif
               endif

               if(il2 == ibpoin(ilk) .and. il0 .gt. il1) then
                  il2 = ibpoin(ilk-1)+1
                  if(abs(xb(ibpoin(ilk-1)+1)-xb(ibpoin(ilk)) )   &
                       .gt. 1E-05 .or.  &
                       abs(yb(ibpoin(ilk-1)+1)-yb(ibpoin(ilk)) )   &
                       .gt. 1E-05 )then
                     print *,'error in profile, the end point is not'
                     print *,'the first, or bad reordering'
                     print *,'4',il0,il1,il2
                     print *,x(i0),y(i0),i0
                     print *,x(i),y(i),i
                     print *,x(i2),y(i2),i2
                     stop
                  endif
               endif
               !     always il0 < il1 < il2 due to the orientation
               !     but may be il1 or il2 > ibpoin(ilk)
               if( il2 .gt. il0 .and. il2 - il0 .le. 20) then
                  ilo = il2 - il0 - 1
                  do 210 l=il0+1,il2-1
                     xloc(l-il0) = xb(l)
                     yloc(l-il0) = yb(l)
                     iloc(l-il0) = l
210               enddo
               endif
               if(il0 .gt. il2 .and. il2+ibpoin(ilk)-il0-2 .le. 20) then
                  ilo = il2 + ibpoin(ilk) - il0 -2
                  do 220 l=il0+1,ibpoin(ilk)-1
                     xloc(l-il0) = xb(l)
                     yloc(l-il0) = yb(l)
                     iloc(l-il0) = l
220               enddo
                  do l=1,il2-1
                     xloc(l+ibpoin(ilk) - il0 -1) = xb(l)
                     yloc(l+ibpoin(ilk) - il0 -1) = yb(l)
                     iloc(l+ibpoin(ilk) - il0 -1) = l
                  enddo
               endif
               
               if(il2 .gt. il0 .and. il2 - il0 .gt. 20) then
                  rkg = 1.*(il2-il0 -2) /19
                  do 230 k=0,19
                     l = il0 + 1 +int(k*rkg+0.5)
                     xloc(k+1) = xb(l)
                     yloc(k+1) = yb(l)
                     iloc(k+1) = l
230               enddo
                  ilo = 20
               endif
               if(il0.gt.il2 .and. il2+ibpoin(ilk)-il0-2 .gt. 20) then
                  rkg = 1.*(il2-il0+ibpoin(ilk) -2) /19
                  do 240 k=0,19
                     l = il2 + 1 +int(k*rkg+0.5)
                     if(l .ge. ibpoin(ilk)) l = l -ibpoin(ilk)+1
                     iprint = 0
                     if(iprint == 1) then
                        print *,'&&&',i,l,ibpoin(ilk)
                        print *,x(i),y(i),i
                        print *,x(i0),y(i0),i0
                        print *,x(i2),y(i2),i2
                        print *,'______________________'
                        do kk1 =1,abs(icyc(i,1))
                           print *,x(icyc(i,kk1+1)),y(icyc(i,kk1+1)),  &
                                icyc(i,kk1+1)
                        enddo
                        print *,'______________________'
                        print *,xb(il0),yb(il0),il0
                        print *,xb(il1),yb(il1),il1
                        print *,xb(il2),yb(il2),il2
                        print *,'Mischmatch detected in MOVING.f????' 
                     endif
                     xloc(k+1) = xb(l)
                     yloc(k+1) = yb(l)
                     iloc(k+1) = l
240               enddo
                  ilo = 20
                  !                  stop
               endif
            endif

            qualityE = 0.
            qualityrez = 1E+35

            if(ibp(i,1) == 0 ) then
!     for nonperiodical points
               jbet = 0
               do 140 k=1,ilo
                  quality1 = 0.
                  quality1rez = 0.
                  xnew = xloc(k)
                  ynew = yloc(k)
                  do 130 j=1,len-1
                     j1 = mod(j,len )+1
                     ii1 = icyc(i,j+1)
                     ii2 = icyc(i,j1+1)
                     xi = x(icyc(i,j+1))
                     yi = y(icyc(i,j+1))
                     xi1 = x(icyc(i,j1+1))
                     yi1 = y(icyc(i,j1+1))
!     check the orientation and bad conditionality
                     detiii = xi*(yi1-ynew)+xi1*(ynew-yi)+xnew*(yi-yi1)
                     reps = AMA%pos*( ((xi1-xi)**2 + (yi1-yi)**2) +  &
                          ((xi-xnew)**2 + (yi-ynew)**2) +  &
                          ((xi1-xnew)**2 + (yi1-ynew)**2) )
                     
                     if (detiii .le. reps) goto 140
                     call POS1TEST_77(xi,yi,xi1,yi1,xnew,ynew,itet)
                     if(itet == 1) then
                        goto 140
                     endif

                     acc = ACCUTE_77(xi,yi,xi1,yi1,xnew,ynew,5.)
                     
                     rl1 = (0.5*((rga(i)+rga(ii1))*(xi-xnew)*(xi-xnew)+  &
                          2*(rgb(i)+rgb(ii1))*(xi-xnew)*(yi-ynew)+  &
                          (rgc(i)+rgc(ii1))*(yi-ynew)*(yi-ynew)))**0.5
                     
                     quality1rez = quality1rez + (rl1 - 1.732050807)**2
                     !     ... NEW ACCUTE
!                     quality1rez = quality1rez 
!     *                    + (rl1 - 1.732050807)**2*acc

                     if(j == len -1) then
                        rl1 = (0.5*((rga(i)+rga(ii2))*  &
                             (xi1-xnew)*(xi1-xnew)+  &
                             2*(rgb(i)+rgb(ii2))*  &
                             (xi1-xnew)*(yi1-ynew)+  &
                             (rgc(i)+rgc(ii2))*  &
                             (yi1-ynew)*(yi1-ynew)))**0.5
                        
                        quality1rez = quality1rez +   &
                             (rl1 - 1.732050807)**2
                     endif
                     
                     ice = 0
130               enddo
                  quality1rez = quality1rez/len
                  
                  if(quality1rez .lt. qualityrez ) then
                     qualityrez = quality1rez
                     jbet = k
                  endif
140            enddo
               if(jbet .gt. 0) then
                  x(i) = xloc(jbet)
                  y(i) = yloc(jbet)
                  ibb(i,1) = iloc(jbet)
               endif
            else
               !     for periodical points
               i1 = ibp(i,1)
               if(ibp(i,2) == 1) then
                  xperreal = AMA%xper(1,1)
                  yperreal = AMA%xper(1,2)
               else
                  xperreal = AMA%xper(2,1)
                  yperreal = AMA%xper(2,2)
               endif                  
               len1 = abs(icyc(i1,1))
               if(abs(x(i1) + xperreal - x(i) ) .lt. 1E-05 .and.  &
                    abs(y(i1) + yperreal - y(i) ) .lt. 1E-05 ) then
                  imov = 1
               elseif(abs(x(i1) - xperreal - x(i) ) .lt. 1E-05 .and.  &
                    abs(y(i1) - yperreal - y(i) ) .lt. 1E-05 ) then
                  imov = -1
               else
                  print *,'BAD in moving in periodical points'
                  print *,i,i1,imov
                  print *,x(i),y(i),xperreal
                  print *,x(i1),y(i1),yperreal
                  print *,abs(x(i1) + xperreal - x(i) ),  &
                       abs(y(i1) + yperreal - y(i) ),  &
                       abs(x(i1) - xperreal - x(i) ),  &
                       abs(y(i1) - yperreal - y(i) )
                  stop
               endif
               
               jbet = 0
               do 340 k=1,ilo
                  quality1 = 0.
                  quality1rez = 0.
                  xnew = xloc(k)
                  ynew = yloc(k)
                  do 330 j=1,len-1
                     j1 = mod(j,len )+1
                     ii1 = icyc(i,j+1)
                     ii2 = icyc(i,j1+1)
                     xi = x(icyc(i,j+1))
                     yi = y(icyc(i,j+1))
                     xi1 = x(icyc(i,j1+1))
                     yi1 = y(icyc(i,j1+1))
!     check the orientation and bad conditionality
                     detiii = xi*(yi1-ynew)+xi1*(ynew-yi)+xnew*(yi-yi1)
                     reps = AMA%pos*( ((xi1-xi)**2 + (yi1-yi)**2) +  &
                          ((xi-xnew)**2 + (yi-ynew)**2) +  &
                          ((xi1-xnew)**2 + (yi1-ynew)**2) )
                     
                     if (detiii .le. reps) goto 340
                     call POS1TEST_77(xi,yi,xi1,yi1,xnew,ynew,itet)
                     if(itet == 1) then
                        goto 340
                     endif
                     
                     acc = ACCUTE_77(xi,yi,xi1,yi1,xnew,ynew, 1.)

                     rl1 = (0.5*((rga(i)+rga(ii1))*(xi-xnew)*(xi-xnew)+  &
                          2*(rgb(i)+rgb(ii1))*(xi-xnew)*(yi-ynew)+  &
                          (rgc(i)+rgc(ii1))*(yi-ynew)*(yi-ynew)))**0.5
                     
                     quality1rez = quality1rez + (rl1 - 1.732050807)**2

!     ... NEW ACCUTE
!                     quality1rez = quality1rez 
!     *                    +(rl1 - 1.732050807)**2*acc

                     if(j == len -1) then
                        rl1 = (0.5*((rga(i)+rga(ii2))*  &
                             (xi1-xnew)*(xi1-xnew)+  &
                             2*(rgb(i)+rgb(ii2))*  &
                             (xi1-xnew)*(yi1-ynew)+  &
                             (rgc(i)+rgc(ii2))*  &
                             (yi1-ynew)*(yi1-ynew)))**0.5
                        
                        quality1rez = quality1rez +   &
                             (rl1 - 1.732050807)**2
                     endif
                     
330               enddo
                  xnew = xloc(k) - imov*xperreal
                  ynew = yloc(k) - imov*yperreal
                  do 350 j=1,len1-1
                     j1 = mod(j,len1 )+1
                     ii1 = icyc(i1,j+1)
                     ii2 = icyc(i1,j1+1)
                     xi = x(icyc(i1,j+1))
                     yi = y(icyc(i1,j+1))
                     xi1 = x(icyc(i1,j1+1))
                     yi1 = y(icyc(i1,j1+1))
!     check the orientation and bad conditionality
                     detiii = xi*(yi1-ynew)+xi1*(ynew-yi)+xnew*(yi-yi1)
                     reps = AMA%pos*( ((xi1-xi)**2 + (yi1-yi)**2) +  &
                          ((xi-xnew)**2 + (yi-ynew)**2) +  &
                          ((xi1-xnew)**2 + (yi1-ynew)**2) )
                     
                     if (detiii .le. reps) goto 340
                     call POS1TEST_77(xi,yi,xi1,yi1,xnew,ynew,itet)
                     if(itet == 1) then
                        goto 340
                     endif

                     acc = ACCUTE_77(xi,yi,xi1,yi1,xnew,ynew, 1.)
                     
                     rl1 = (0.5*((rga(i)+rga(ii1))*(xi-xnew)*(xi-xnew)+  &
                          2*(rgb(i)+rgb(ii1))*(xi-xnew)*(yi-ynew)+  &
                          (rgc(i)+rgc(ii1))*(yi-ynew)*(yi-ynew)))**0.5
                     
                     quality1rez = quality1rez + (rl1 - 1.732050807)**2
!     ... NEW ACCUTE
!                     quality1rez = quality1rez 
!     *                    +(rl1 - 1.732050807)**2*acc

                     if(j == len1 -1) then
                        rl1 = (0.5*((rga(i)+rga(ii2))*  &
                             (xi1-xnew)*(xi1-xnew)+  &
                             2*(rgb(i)+rgb(ii2))*  &
                             (xi1-xnew)*(yi1-ynew)+  &
                             (rgc(i)+rgc(ii2))*  &
                             (yi1-ynew)*(yi1-ynew)))**0.5
                        
                        quality1rez = quality1rez +   &
                             (rl1 - 1.732050807)**2
                     endif
350               enddo
                  quality1rez = quality1rez/(len+len1+2)
                  
                  if(quality1rez .lt. qualityrez ) then
                     qualityrez = quality1rez
                     jbet = k
                  endif
340            enddo
               if(jbet .gt. 0) then
                  x(i) = xloc(jbet)
                  y(i) = yloc(jbet)
                  x(i1) = xloc(jbet) - imov*xperreal
                  y(i1) = yloc(jbet) - imov*yperreal
               endif
            endif

         endif
21       continue
100   enddo
      return
    end subroutine

    subroutine REM_BOUND_77(ndim,melem,nelem,mpoin,npoin,  &
         maxdeg,x,y,lnd,iae,  &
         ra,icyc,rminerr,imin,rmaxrez,  &
         icha,nserr,wp,lbn,itc,ibc,nbelm,mbelm,  &
         rga,rgb,rgc,ibp,ibb,ibpoin,ipoint,nbp)
      dimension lnd(melem,3),iae(melem,3),x(mpoin),y(mpoin),  &
           ra(mpoin),wp(mpoin,ndim+1),rga(mpoin),rgb(mpoin),rgc(mpoin),  &
           icyc(mpoin,maxdeg),  &
           nserr(melem*3,2),lbn(mbelm,2),itc(mbelm),  &
           ibc(mbelm),ibp(mpoin,2),ibb(mpoin,3),ibpoin(0:nbp)
      integer nsr(3),locyc(30)

      icha = 0
      ipoc = 0
 10   ipoc = ipoc + 1
      if(ipoc .gt. npoin) return
      i = ipoc

      if(ibb(i,1) .gt.0 ) then
         ikk1 = ibb(i,2)
         if(ibb(i,1) == ibpoin(ikk1-1)+1  .or.  &
              ibb(i,1) == ibpoin(ikk1) .or. ibb(i,3) == -1) then
!     no removing, initial or final point of profiles
            goto 10
         endif
      endif       
      if(icyc(i,1) == -3 ) then
         i0 = icyc(i,2)
         i1 = icyc(i,3)
         i2 = icyc(i,4)
         if(icyc(i1,1) .lt. 0) goto 101
         det = x(i0)*(y(i)-y(i2)) + x(i)*(y(i2)-y(i0)) +   &
              x(i2)*(y(i0)-y(i))

         rep =  ((x(i0)-x(i))**2 + (y(i0)-y(i))**2) +  &
              ((x(i)-x(i2))**2 + (y(i)-y(i2))**2) +  &
              ((x(i2)-x(i0))**2 + (y(i2)-y(i0))**2) 

         if (abs(det)/rep .lt. 1e-02) then
!                  the points i0, i, i2 are in the straight, 
!                  we can remove one triangle

            det = x(i0)*(y(i1)-y(i2)) + x(i1)*(y(i2)-y(i0)) +   &
                 x(i2)*(y(i0)-y(i1))
            reps = AMA%pos*( ((x(i0)-x(i1))**2 + (y(i0)-y(i1))**2) +  &
                    ((x(i1)-x(i2))**2 + (y(i1)-y(i2))**2) +  &
                    ((x(i2)-x(i0))**2 + (y(i2)-y(i0))**2) )


            if( det .le. reps ) then
!               print *,'violation of positivity'
               goto 100
            endif
            call POS1TEST_77(x(i0),y(i0),x(i1),y(i1),x(i2),y(i2),itet)
            if(itet == 1) then
               goto 100
            endif

            ibo = 0
            call QUA2_77(x(i0),y(i0),wp(i0,1),rga(i0),rgb(i0),rgc(i0),  &
                 x(i1),y(i1),wp(i1,1),rga(i1),rgb(i1),rgc(i1),  &
                 x(i2),y(i2),wp(i2,1),rga(i2),rgb(i2),rgc(i2),  &
                 err0,errrez0,ice,ibo)
            call QUA2_77(x(i0),y(i0),wp(i0,1),rga(i0),rgb(i0),rgc(i0),  &
                 x(i1),y(i1),wp(i1,1),rga(i1),rgb(i1),rgc(i1),  &
                 x(i),y(i),wp(i,1),rga(i),rgb(i),rgc(i),  &
                 err2,errrez1,ice,ibo)
            call QUA2_77(x(i),y(i),wp(i,1),rga(i),rgb(i),rgc(i),  &
                 x(i1),y(i1),wp(i1,1),rga(i1),rgb(i1),rgc(i1),  &
                 x(i2),y(i2),wp(i2,1),rga(i2),rgb(i2),rgc(i2),  &
                 err1,errrez2,ice,ibo)
            
            iyes = 0
            if( 1.8*errrez0 .lt. errrez1 + errrez2) then
               iyes = 1
            endif

            if(ibp(i,1) .gt. 0 .and. iyes == 1) then
               k=ibp(i,1)
               kik = k
               if(icyc(k,1) .ne. -3) goto 100
               k0 = icyc(k,2)
               k1 = icyc(k,3)
               k2 = icyc(k,4)
               if(icyc(k1,1) .lt. 0) goto 100
               det = x(k0)*(y(k)-y(k2)) + x(k)*(y(k2)-y(k0)) +   &
                    x(k2)*(y(k0)-y(k))
               rep =  ((x(k0)-x(k))**2 + (y(k0)-y(k))**2) +  &
                    ((x(k)-x(k2))**2 + (y(k)-y(k2))**2) +  &
                    ((x(k2)-x(k0))**2 + (y(k2)-y(k0))**2) 

               if (abs(det)/rep .lt. 1e-02) then
!     the points k0, k, k2 are in the straight, 
!     we can remove one triangle

                  det = x(k0)*(y(k1)-y(k2)) + x(k1)*(y(k2)-y(k0)) +   &
                       x(k2)*(y(k0)-y(k1))
                  reps = AMA%pos*( ((x(k0)-x(k1))**2 + (y(k0)-y(k1))**2) +  &
                       ((x(k1)-x(k2))**2 + (y(k1)-y(k2))**2) +  &
                       ((x(k2)-x(k0))**2 + (y(k2)-y(k0))**2) )
                  if( det .le. reps ) then
!     print *,'violation of positivity'
                     goto 100
                  endif
                  call POS1TEST_77(x(k0),y(k0),x(k1),y(k1),  &
                       x(k2),y(k2),itet)
                  if(itet == 1) then
                     goto 100
                  endif

                  ibo = 0
                  call QUA2_77(x(k0),y(k0),wp(k0,1),  &
                       rga(k0),rgb(k0),rgc(k0),  &
                       x(k1),y(k1),wp(k1,1),rga(k1),rgb(k1),rgc(k1),  &
                       x(k2),y(k2),wp(k2,1),rga(k2),rgb(k2),rgc(k2),  &
                       err0,errrez0,ice,ibo)
                  call QUA2_77(x(k0),y(k0),wp(k0,1),  &
                       rga(k0),rgb(k0),rgc(k0),  &
                       x(k1),y(k1),wp(k1,1),rga(k1),rgb(k1),rgc(k1),  &
                       x(k),y(k),wp(k,1),rga(k),rgb(k),rgc(k),  &
                       err2,errrez1,ice,ibo)
                  call QUA2_77(x(k),y(k),wp(k,1),rga(k),rgb(k),rgc(k),  &
                       x(k1),y(k1),wp(k1,1),rga(k1),rgb(k1),rgc(k1),  &
                       x(k2),y(k2),wp(k2,1),rga(k2),rgb(k2),rgc(k2),  &
                       err1,errrez2,ice,ibo)
            

                  if( 1.8*errrez0 .lt. errrez1 + errrez2 ) then
                     iyes = 2
                  else
                     iyes = 0
                  endif
               endif  
            endif

            if( iyes .ge. 1) then
!     better quality, we  remove a boundary points
 999           i0 = icyc(i,2)
               i1 = icyc(i,3)
               i2 = icyc(i,4)
               itest = 0
               do 20 kk=1,nelem
                  do 30 j=1,3
                     j1 = mod(j,3) + 1
                     j2 = mod(j1,3) + 1
                     if(lnd(kk,j) == i0 .and. lnd(kk,j1) == i1 .and.  &
                          lnd(kk,j2) == i) then
                        k1 = kk
                        kj1 = j
                        kj2 = j1
                        kj3 = j2
                        itest = itest + 1
                     endif
!                     print *,'k1=',k1,kk,nelem,i0,i1,i
                     if(lnd(kk,j) == i2 .and. lnd(kk,j1) == i .and.  &
                          lnd(kk,j2) == i1) then
                        k2 = kk
                        kjj1 = j
                        kjj2 = j1
                        kjj3 = j2
                        itest = itest + 1
                     endif
                     if(itest == 2) goto 25
 30               enddo
 20            enddo
 25            continue

               ka1 = iae(k1,kj1)
               if(ka1 .gt. 0) then
                  do 40 j=1,3
                     if(iae(ka1,j) == k1) ka1j = j
 40               enddo
               endif

               ka2 = iae(k2,kjj3)
               if(ka2 .gt. 0) then
                  do 45 j=1,3
                     if(iae(ka2,j) == k2) ka2j = j
 45               enddo
               endif

               lnd(k1,kj3) = i2
               iae(k1,kj2) = ka2
               if(ka2 .gt. 0) iae(ka2,ka2j) = k1

               do ie = 1,nelem
                  do j=1,3
                     if(iae(ie,j) == nelem .and. k2 .ne. nelem)   &
                          iae(ie,j) = k2
                     if(lnd(ie,j) == npoin) lnd(ie,j) = i
                  enddo
               enddo

               do 70 ij =1,3
                  if(k2 .ne. nelem) then
                     lnd(k2,ij) = lnd(nelem,ij)
                     iae(k2,ij) = iae(nelem,ij)
                  endif
 70            enddo

               
               if (kik == npoin) kik = i
               
               x(i) = x(npoin)
               y(i) = y(npoin)
               ibb(i,1) = ibb(npoin,1)
               ibb(i,2) = ibb(npoin,2)
               ibb(i,3) = ibb(npoin,3)
               rga(i) = rga(npoin)
               rgb(i) = rgb(npoin)
               rgc(i) = rgc(npoin)
               wp(i,1) = wp(npoin,1)
               ibp(i,1) = ibp(npoin,1)
               ibp(i,2) = ibp(npoin,2)
               do iy=1,npoin
                  if(ibp(iy,1) == npoin) ibp(iy,1) = i
               enddo
                  

               icyc(i0,abs(icyc(i0,1))+1) = i2
               icyc(i2,2) = i0
               do 200 k=1,icyc(i1,1)
                  if(icyc(i1,k+1) == i) then
                     kk = k
                     goto 210
                  endif
 200           enddo
 210           continue

               do 220 k=kk,icyc(i1,1)-1
                  icyc(i1,k+1) = icyc(i1,k+2)
 220           enddo
               icyc(i1,1) = icyc(i1,1) - 1

               do is = 1,npoin
                  do j=1,abs(icyc(is,1))
                     if(icyc(is,j+1) == npoin)icyc(is,j+1) =i
                  enddo
               enddo  
               
               do 260 j=1,abs(icyc(npoin,1))+1
                  icyc(i,j)=icyc(npoin,j)
 260           enddo



!     reparation of lbn, ibc, itc
               do ib=1,nbelm
                  if(lbn(ib,2) == i) then
                     lbn(ib,2) = i0
                     itc(ib) = k1
                     ibstart1 = ib
                  elseif(lbn(ib,1) == i) then
                     ibstart2 = ib
                  endif
               enddo

!               print *,'ibstart =',ibstart1,ibstart2


               do 698 ib=1,nbelm
                  if(lbn(ib,1) == npoin) lbn(ib,1) = i
                  if(lbn(ib,2) == npoin) lbn(ib,2) = i
 698           enddo
               do 710 ib1 = ibstart2,nbelm-1
                  lbn(ib1,1) = lbn(ib1+1,1)
                  lbn(ib1,2) = lbn(ib1+1,2)
                  ibc(ib1) = ibc(ib1+1)
                  itc(ib1) = itc(ib1+1)
 710           enddo
               do ib1 = 1,nbelm-1
                  if(itc(ib1) == nelem) then
                     itc(ib1) = k2
                  endif
               enddo

               nbelm = nbelm - 1
               nelem = nelem - 1
               npoin = npoin - 1

               icha = icha + 1
                     if(AMA%ifig .ge. 0 .and. mod(icha,5) == 0 ) then
                        call PLOT1_77(melem,nelem,mpoin,npoin,x,y,lnd)
                        AMA%ifig = AMA%ifig + 1
                     endif


               if(iyes == 2) then
                  iyes = 3
                  i = kik
                  goto 999
               endif
            endif
 100        continue
         endif
101      continue
      endif
      if(ipoc .lt. npoin) goto 10
      return
    end subroutine REM_BOUND_77


    subroutine REMOVE_77(ndim,melem,nelem,mpoin,npoin,  &
         maxdeg,x,y,lnd,iae,  &
         ra,icyc,rminerr,imin,rmaxrez,  &
         icha,nserr,wp,lbn,ibc,itc,nbelm,mbelm,  &
         rga,rgb,rgc,ipoint,nbp,xb,yb,ibb,ibpoin,ibp)
      dimension lnd(melem,3),iae(melem,3),x(mpoin),y(mpoin),  &
           ra(mpoin),wp(mpoin,ndim+1),rga(mpoin),rgb(mpoin),rgc(mpoin),  &
           icyc(mpoin,maxdeg),  &
           nserr(melem*3,2),lbn(mbelm,2),itc(mbelm),ibc(mbelm),  &
           xb(ipoint),yb(ipoint),ibp(mpoin,2),ibb(mpoin,3),ibpoin(0:nbp)
      integer nsr(3),locyc(50)
      real*8 rlen0,rlen1,rlen2,rcos,ss,ve1x,ve1y,ve2x,ve2y
      real rmin(3)
      integer jmin(3)
      
      rlmin2 = 1.33
      
      icha = 0
      nelemold = nelem
      ipoc = 1
      do 20 iyi=1,nelemold
         i = ipoc
         if( i .le. nelem) then
            do 27 jdo=1,3
               j = jdo
               j1 = mod(j,3) + 1
               k1 = lnd(i,j)
               k2 = lnd(i,j1)
               x1 = x(k1)
               y1 = y(k1)
               x2 = x(k2)
               y2 = y(k2)
               a = (rga(k1) + rga(k2) )/2
               b = (rgb(k1) + rgb(k2) )/2
               c = (rgc(k1) + rgc(k2) )/2
               rmin(j) = (a*(x1-x2)*(x1-x2) + c*(y1-y2)*(y1-y2)  &
                    + 2*b*(x1-x2)*(y1-y2))
               jmin(j) = j
 27         enddo

            do 25 k=1,3
               do 26 l=1,2
                  if(rmin(l) .gt. rmin(l+1) ) then
                     rminhelp = rmin(l)
                     rmin(l) = rmin(l+1)
                     rmin(l+1) = rminhelp
                     jminhelp = jmin(l)
                     jmin(l) = jmin(l+1)
                     jmin(l+1) = jminhelp
                  endif
 26            enddo
 25         enddo
            
            do 28 jdo=1,3
               j = jmin(jdo)
               j1 = mod(j,3) + 1
               k1 = lnd(i,j)
               k2 = lnd(i,j1)
               x1 = x(k1)
               y1 = y(k1)
               x2 = x(k2)
               y2 = y(k2)
               ig1 = 0
               ig2 = 0
               do 21 ib=1,nbelm
                  if(lbn(ib,1) == k1) ig1 = 1
                  if(lbn(ib,1) == k2) ig2 = 1
 21            enddo
               
               if(ig1 == 1 .and. icyc(k1,1) .ge. 0) then
                  print *,'ERROR in REMOVE - boundary point 1'
                  kl = k1
                  write(*,'(i5,2e12.4,2i5)') kl,x(kl),y(kl),  &
                       icyc(kl,1),ig1
                  do kl1 = 1,abs(icyc(kl,1))
                     kll = icyc(kl,kl1+1)
                     write(*,'(2e12.4,i5)') x(kll),y(kll),kll
                  enddo
                  stop
               endif
               if(ig2 == 1 .and. icyc(k2,1) .ge. 0) then
                  print *,'ERROR in REMOVE - boundary point 2'
                  kl = k2
                  write(*,'(i5,2e12.4,2i5)') kl,x(kl),y(kl),  &
                       icyc(kl,1),ig2
                  do kl1 = 1,abs(icyc(kl,1))
                     kll = icyc(kl,kl1+1)
                     write(*,'(2e12.4,i5)') x(kll),y(kll),kll
                  enddo
                  
                  do ib=1,nbelm
                     if(lbn(ib,1) == kl) then
                        print *,'^^',lbn(ib,1),ib,x(lbn(ib,1)),  &
                             y(lbn(ib,1)),ibb(kl,1)
                     endif
                  enddo
                 stop
               endif
               
               if(iae(i,j) .gt. 0) then
!     for non boundary sides
                  if(ig1+ig2 .le. 1 .and. rmin(jdo) .le. rlmin2) then
!     we  remove this side (may be)

                     j2 = mod(j1,3)+1
                     ii = iae(i,j)
                     jj = 0
                     do 35 jjj=1,3
                        if(iae(ii,jjj) == i) jj = jjj
 35                  enddo
                     if(jj == 0) then
                        print *,'error 1342'
                        stop
                     endif
                     jj1 = mod(jj,3)+1
                     jj2 = mod(jj1,3)+1
                     k3 = lnd(i,j2)
                     k4 = lnd(ii,jj2)

                     if(k2.ne.lnd(ii,jj) .or. k1.ne.lnd(ii,jj1))then
                        print *,'ERRROR in REMOVE'
                        print *,i,k1,k2,k3,k4
                        print *,i,lnd(i,1),lnd(i,2),lnd(i,3)
                        print *,ii,lnd(ii,1),lnd(ii,2),lnd(ii,3)
                        print *,x(k1),y(k1)
                        print *,x(k2),y(k2)
                        print *,x(k3),y(k3)
                        print *,x(k4),y(k4)
                        stop
                     endif
                     
                     if(iae(i,j1) .gt. 0) then
                        ia1 = iae(i,j1)
                        do 40 kk =1,3
                           if(iae(ia1,kk) == i) ja1 = kk
 40                     enddo
                     else
                        ia1 = -2
                     endif
                     if(iae(i,j2) .gt. 0) then
                        ia2 = iae(i,j2)
                        do 50 kk =1,3
                           if(iae(ia2,kk) == i) ja2 = kk
 50                     enddo
                     else
                        ia2 = -2
                     endif
                     if(iae(ii,jj1) .gt. 0) then
                        iia1 = iae(ii,jj1)
                        do 60 kk =1,3
                           if(iae(iia1,kk) == ii) jja1 = kk
 60                     enddo
                     else
                        iia1 = -2
                     endif
                     if(iae(ii,jj2) .gt. 0) then
                        iia2 = iae(ii,jj2)
                        do 70 kk =1,3
                           if(iae(iia2,kk) == ii) jja2 = kk
 70                     enddo
                     else
                        iia2 = -2
                     endif
!     here we check, so we will not obtain the triangle with two 
!     boundary sides
                     if(ia2 .lt. 0)then
                        if(iae(ia1,1) .lt. 0 .or. iae(ia1,2) .lt. 0   &
                             .or.iae(ia1,3) .lt. 0 ) then
                           goto 28
                        endif
                     endif
                     if(ia1 .lt. 0)then
                        if(iae(ia2,1) .lt. 0 .or. iae(ia2,2) .lt. 0   &
                             .or.iae(ia2,3) .lt. 0 )then
                           goto 28
                        endif
                     endif

                     if(iia2 .lt. 0)then
                        if(iae(iia1,1) .lt. 0 .or. iae(iia1,2) .lt. 0   &
                             .or.iae(iia1,3) .lt. 0 ) then
                          goto 28
                       endif
                     endif
                     if(iia1 .lt. 0)then
                        if(iae(iia2,1) .lt. 0 .or. iae(iia2,2) .lt. 0   &
                             .or.iae(iia2,3) .lt. 0 ) then
                          goto 28
                       endif
                     endif

!     here we check that after emoving don't arise a corner triangle
                     if(iia1 == ia2 .and. ia2 .gt. 0 .and.  &
                          iia2 .lt. 0 .and. ia1 .lt. 0) then
                        goto 28
                     endif
                     if(ia1 == iia2 .and. iia2 .gt. 0 .and.  &
                          ia2 .lt. 0 .and. iia1 .lt. 0) then
                        goto 28
                     endif

                     if( ig1 == 1) then
!     k1 is on the boundary, no changes
                        xk1 = x(k1)
                        yk1 = y(k1)
                        wpk1 = wp(k1,1)
                        rgak = rga(k1)
                        rgbk = rgb(k1)
                        rgck = rgc(k1)
                     elseif( ig2 == 1) then
!     k2 on the boundary
                        xk1 = x(k2)
                        yk1 = y(k2)
                        wpk1 = wp(k2,1)
                        rgak = rga(k2)
                        rgbk = rgb(k2)
                        rgck = rgc(k2)
                     else
                        xk1 = (x(k1) + x(k2) )/2
                        yk1 = (y(k1) + y(k2) )/2
                        wpk1 = (wp(k1,1) + wp(k2,1) )/2
                        rgak = (rga(k1) + rga(k2) )/2
                        rgbk = (rgb(k1) + rgb(k2) )/2
                        rgck = (rgc(k1) + rgc(k2) )/2


                  
                     endif
                  
!     we check that if some new triangle will satisfy the positivity

                     xx0 = xk1
                     yy0 = yk1
!     for point k1
                     if(icyc(k1,1) .gt. 0) then
                        ilen = icyc(k1,1)
                     else
                        ilen = abs(icyc(k1,1)) - 1
                     endif
                     do 450 kk=1,ilen
                        kk1 = mod(kk,icyc(k1,1)) + 1
                        if( icyc(k1,kk+1) .ne. k2 .and.   &
                             icyc(k1,kk1+1) .ne. k2) then
                           xx1 = x(icyc(k1,kk+1))
                           yy1 = y(icyc(k1,kk+1))
                           xx2 = x(icyc(k1,kk1+1))
                           yy2 = y(icyc(k1,kk1+1))
                           reps = AMA%pos*( ((xx0-xx2)**2 + (yy0-yy2)**2) +  &
                                ((xx1-xx2)**2 + (yy1-yy2)**2) +  &
                                ((xx0-xx1)**2 + (yy0-yy1)**2) )
                           det123 = xx0*(yy1-yy2) + xx1*(yy2-yy0) +   &
                                xx2*(yy0-yy1) 
                           
                           if( det123 .le. reps) then
!     violation of positivity, go to next j
                              goto 28
                           endif
                           call POS1TEST_77(xx0,yy0,xx1,yy1,xx2,yy2,itet)
                           if(itet == 1) then
                              goto 28
                           endif
                        endif
 450                 enddo
!     now for k2 

                     if(icyc(k2,1) .gt. 0) then
                        ilen = icyc(k2,1)
                     else
                        ilen = abs(icyc(k2,1)) - 1
                     endif
                     do 460 kk=1,ilen
                        kk2 = mod(kk,icyc(k2,1)) + 1
                        if( icyc(k2,kk+1) .ne. k1 .and.   &
                             icyc(k2,kk2+1) .ne. k1) then
                           xx1 = x(icyc(k2,kk+1))
                           yy1 = y(icyc(k2,kk+1))
                           xx2 = x(icyc(k2,kk2+1))
                           yy2 = y(icyc(k2,kk2+1))
                           reps = AMA%pos*( ((xx0-xx2)**2 + (yy0-yy2)**2) +  &
                                ((xx1-xx2)**2 + (yy1-yy2)**2) +  &
                                ((xx0-xx1)**2 + (yy0-yy1)**2) )
                           det123 = xx0*(yy1-yy2) + xx1*(yy2-yy0) +   &
                                xx2*(yy0-yy1) 
                           
                           if( det123 .le. reps) then
!     violation of positivity, go to next j
                              goto 28
                           endif
                           call POS1TEST_77(xx0,yy0,xx1,yy1,xx2,yy2,itet)
                           if(itet == 1) then
                              goto 28
                           endif


                        endif
 460                 enddo

                     kk1 = k1
                     kk2 = k2
                     if( ig2 == 1) then
                        kk1 = k2
                        kk2 = k1
                     endif

                     x(kk1) = xk1
                     y(kk1) = yk1

                     wp(kk1,1) = wpk1
                     rga(kk1) = rgak
                     rgb(kk1) = rgbk
                     rgc(kk1) = rgck

!     shiftting of all array after removing k2 

                     do 6745 ip=1,npoin
                        if(ibp(ip,1) == kk2 .and.kk1 .ne. npoin)   &
                             ibp(ip,1) = kk1
                        if(ibp(ip,1) == npoin .and.kk2 .ne. npoin)   &
                             ibp(ip,1) = kk2
 6745                enddo
                     
                     
                     ibp(kk2,1) = ibp(npoin,1) 
                     ibp(kk2,2) = ibp(npoin,2) 
                     x(kk2) = x(npoin)
                     y(kk2) = y(npoin)
                     ibb(kk2,1) = ibb(npoin,1)
                     ibb(kk2,2) = ibb(npoin,2)
                     ibb(kk2,3) = ibb(npoin,3)
                     wp(kk2,1) = wp(npoin,1)
                     rga(kk2) = rga(npoin)
                     rgb(kk2) = rgb(npoin)
                     rgc(kk2) = rgc(npoin)

                     iemin = min (i,ii)
                     iemax = max (i,ii)

                     if(ia1 .gt. 0) iae(ia1,ja1) = ia2
                     if(ia2 .gt. 0) iae(ia2,ja2) = ia1
                     if(iia1 .gt. 0) iae(iia1,jja1) = iia2
                     if(iia2 .gt. 0) iae(iia2,jja2) = iia1
                  
                     do 355 ie=1,nelem
                        do 360 kj=1,3
                           if(lnd(ie,kj) == kk2 .and. kk1 .ne. npoin)  &
                                lnd(ie,kj) = kk1
                           if(lnd(ie,kj) == npoin .and. kk2 .ne.npoin)  &
                                lnd(ie,kj) = kk2
 360                    enddo
                        if( iemax .lt. nelem-1 ) then
                           do kj=1,3
                              if(iae(ie,kj)==nelem-1) iae(ie,kj)=iemin
                              if(iae(ie,kj) == nelem) iae(ie,kj)=iemax
                           enddo
                        elseif(iemax == nelem-1 ) then
                           do kj=1,3
                              if(iae(ie,kj) == nelem) iae(ie,kj)=iemin
                           enddo
                        elseif(iemax == nelem .and.   &
                                iemin .lt. nelem-1) then
                           do kj=1,3
                              if(iae(ie,kj)==nelem-1) iae(ie,kj)=iemin
                           enddo
                        elseif(iemax == nelem .and.   &
                                iemin == nelem-1) then
!     no performance
                        else
                           print *,'LOGICAL error in REMOVE - 111&*('
                           stop
                        endif
 355                 enddo

                     if( iemax .lt. nelem-1 ) then
                        do kj=1,3
                           lnd(iemin,kj) = lnd(nelem-1,kj)
                           iae(iemin,kj) = iae(nelem-1,kj)
                           lnd(iemax,kj) = lnd(nelem,kj)
                           iae(iemax,kj) = iae(nelem,kj)
                        enddo
                     elseif(iemax == nelem-1 ) then
                        do kj=1,3
                           lnd(iemin,kj) = lnd(nelem,kj)
                           iae(iemin,kj) = iae(nelem,kj)
                        enddo
                     elseif(iemax == nelem .and.   &
                             iemin .lt. nelem-1) then
                        do kj=1,3
                           lnd(iemin,kj) = lnd(nelem-1,kj)
                           iae(iemin,kj) = iae(nelem-1,kj)
                        enddo
                     elseif(iemax == nelem .and.   &
                             iemin == nelem-1) then
!     no performance
                     else
                        print *,'LOGICAL error in REMOVE - &*('
                        stop
                     endif

                     do k=1,nbelm
                        if(itc(k) == i) then
                           if(ia1 .lt. 0 .and. ia2 .gt. 0) then
                              itc(k) = ia2
                           elseif(ia2 .lt. 0 .and. ia1 .gt.0) then
                              itc(k) = ia1
                           else
                              print *,'LOG. EROR in REM 145'
                              stop
                           endif
                        elseif(itc(k) == ii) then
                           if(iia1 .lt. 0 .and. iia2 .gt. 0) then
                              itc(k) = iia2
                           elseif(iia2 .lt. 0 .and. iia1 .gt. 0) then
                              itc(k) = iia1
                           else
                              print *,'LOG. EROR in REM 146'
                              stop
                           endif
                        endif
                     enddo

                     do 300 k=1,nbelm
                        do 310 l=1,2
                           if(lbn(k,l) == kk2 .and. kk1.ne.npoin)   &
                                lbn(k,l) = kk1
                           if(lbn(k,l) == npoin .and. kk2 .ne. npoin)  &
                                lbn(k,l) = kk2
 310                    enddo
                        if( iemax .lt. nelem-1 ) then
                           if(itc(k) == nelem-1)  itc(k)=iemin
                           if(itc(k) == nelem)  itc(k)=iemax
                        elseif(iemax == nelem-1 ) then
                           if(itc(k) == nelem)  itc(k)=iemin
                        elseif(iemax == nelem .and.   &
                                iemin .lt. nelem-1 ) then
                           if(itc(k) == nelem-1)  itc(k)=iemin
                        elseif(iemax == nelem .and.   &
                                iemin == nelem-1 ) then
!     no performance
                        else
                           print *,'ERRor in REMOVE 123'
                           print *,'LOGICAL error'
                           stop
                        endif
 300                 enddo
!     connection of two cykles
                     is = 0
                     do 700 ic=1,abs(icyc(kk2,1))
                        ic1 = mod(ic,abs(icyc(kk2,1)) ) +1
                        if( icyc(kk2,ic1+1) == kk1) is = ic
 700                 enddo
                     if(is == 0) then
                        print *,'ERROR in REMOVE in @#$%'
                        stop
                     endif
                     is1 = mod(is,abs(icyc(kk2,1)) ) +1
                     is2 = mod(is1,abs(icyc(kk2,1)) ) +1
                     ip = icyc(kk2,is+1)
                     ip1 = icyc(kk2,is1+1)
                     ip2 = icyc(kk2,is2+1)
                     
                     if(ip1 .ne. kk1) then
                        print *,'ERROR in $%$%$%',kk1,kk2
                        print *,is,is1,is2
                        print *,ip,ip1,ip2
                        stop
                     endif
                     
                     ir2 = 0
                     ir1 = 0
                     do 710 ic = 1,abs(icyc(kk1,1))
                        if(icyc(kk1,ic+1) == ip2) ir2 = ic
                        if(icyc(kk1,ic+1) == ip) ir1 = ic
 710                 enddo
                     
                     if(ir1 == 0 .or. ir2 == 0) then
                        print *,'ERROR in REMOVE in $$$$'
                        stop
                     endif

                     locyc(1) = 0
                     irr = 0
                     do 720 ic =1,ir2
                        if(icyc(kk1,ic+1) .ne. kk2) then
                           irr = irr + 1
                           locyc(irr+1) = icyc(kk1,ic+1)
                           locyc(1) = locyc(1) + 1
                        endif
 720                 enddo
                     
                     irr = locyc(1)
                     iskk2 = is2
 730                 ikk2 = mod(iskk2, abs(icyc(kk2,1) )) + 1
                     if( icyc(kk2,ikk2+1) .ne. ip) then
                        irr = irr + 1
                        locyc(1) = locyc(1) + 1
                        locyc(irr+1) = icyc(kk2,ikk2+1)
                        iskk2 = ikk2
                        goto 730
                     endif
                     
                     if( ir1 .gt. 2) then
                        do 740 ic= ir1,abs(icyc(kk1,1) )
                           irr = irr + 1
                           locyc(irr+1) = icyc(kk1,ic+1)
                           locyc(1) = locyc(1) + 1
 740                    enddo
                     endif


                     if(locyc(1) + 1  .gt. maxdeg) then
                        print *,'Too long cykles in REMOVE'
                        print *,icyc(kk1,1),icyc(kk2,1),locyc(1),maxdeg
                        stop
                     endif
                     if(icyc(kk1,1) .gt. 0) then
                        icyc(kk1,1) = locyc(1) 
                     else
                        icyc(kk1,1) = -locyc(1) 
                     endif
                     
                     do 750 ic=1,locyc(1)
                        icyc(kk1,ic+1) = locyc(ic+1)
 750                 enddo
!     end of connection of two cykles

                     do 500 ip =1,npoin
                        do 510 ic =1,abs(icyc(ip,1))
                           if(icyc(ip,ic+1) == kk2)icyc(ip,ic+1) = kk1
                           if(icyc(ip,ic+1) == npoin)   &
                                icyc(ip,ic+1) = kk2
 510                    enddo
 500                 enddo
                     
                     do 530 ic=1,abs(icyc(npoin,1)) + 1
                        icyc(kk2,ic) = icyc(npoin,ic)
 530                 enddo
                     
                     do 400 ip1 =1,2
                        if(ip1 == 1) then
                           ip = k3
                           if(k3 == npoin ) ip = kk2
                        endif
                        if(ip1 == 2) then
                           ip = k4
                           if(k4 == npoin ) ip = kk2
                        endif
                        do 410 ic =1,abs(icyc(ip,1))
                           ic1 = mod(ic,abs(icyc(ip,1)) ) + 1
                           if( icyc(ip,ic+1) == icyc(ip,ic1+1) ) then
                              do 420 icc= ic1,abs(icyc(ip,1))-1
                                 icyc(ip,icc+1) = icyc(ip,icc+2)
 420                          enddo
                              if(icyc(ip,1) .gt. 0) then
                                 icyc(ip,1) = icyc(ip,1) - 1
                              else
                                 icyc(ip,1) = icyc(ip,1) + 1
                              endif
                              goto 400
                           endif
 410                    enddo
 400                 enddo
!     end of shiftting

!     all is done, we can continue

                     npoin = npoin - 1
                     nelem = nelem - 2
                     icha = icha + 1

                     if(AMA%ifig .ge. 0 .and. mod(icha,5) == 0 ) then
                        call PLOT1_77(melem,nelem,mpoin,npoin,x,y,lnd)
                        AMA%ifig = AMA%ifig + 1
                     endif

                     
                  endif
               else
!     for boundary sides

                  if( ig1 + ig2 .ne. 2) then
                     print *,'ERROR in REMOVE'
                     print *,'the boundary segment do not have ',  &
                          'the points on the boundary'
                     print *,i,k1,k2,ig1,ig2,icha
                     print *,lnd(i,1),lnd(i,2),lnd(i,3)
                     print *,nelem,npoin,nbelm
                     print *,melem,mpoin,mbelm
                     print *,x1,y1
                     print *,x2,y2
                     stop
                  endif

                  itest2 = 0
                  if(rmin(jdo) .le. rlmin2) then
!     we  remove this edge
                     itest2 = 1
                     ipe = 0
 999                 j1 = mod(j,3) +1
                     j2 = mod(j1,3)+1
                     k1 = lnd(i,j)
                     k2 = lnd(i,j1)
                     k3 = lnd(i,j2)
                     x1 = x(k1)
                     y1 = y(k1)
                     x2 = x(k2)
                     y2 = y(k2)
                     kb0 = 0
                     kb1 = 0
                     kb2 = 0


                     do 610 ib=1,nbelm
                        if(lbn(ib,1) == k1 .and.   &
                             lbn(ib,2) == k2) then
                           kb1 = ib
                        elseif(lbn(ib,2) == k1) then
                           kb0 = ib 
                        elseif(lbn(ib,1) == k2) then
                           kb2 = ib
                        endif
 610                 enddo

                     if( kb0*kb1*kb2 == 0 ) then
                        print *,'The bounary segment does not found'
                        print *,i,k1,k2
                        print *,kb0,kb1,kb2
                        print *,x1,y1,x(k1),y(k1)
                        print *,x2,y2,x(k2),y(k2)
                        stop
                     endif

                     if(ipe == 1) goto 993

!     "sharp angle
                     ikr1 = 0
                     ikr2 = 0
                     ibbk1 = ibb(k1,1)
                     ibbk11 = ibb(k1,2)
                     ibbk2 = ibb(k2,1)
                     ilk = 0
                     ilk1 = 0

                     if(ibb(k1,1) .gt. 0) then 
                        ikk1 = ibb(k1,2)
                        if(ibb(k1,1) == ibpoin(ikk1-1)+1  .or.  &
                             ibb(k1,1) == ibpoin(ikk1) .or.   &
                             ibb(k1,3) == -1) then
!               NO moving of final or initial node of profiles
                           ikr1 = 1
                           ilk = ibb(k1,2)
                           ilk1 = ibb(k1,3)
                        endif
                     endif

                     if(ibb(k2,1) .gt. 0) then 
                        ikk1 = ibb(k2,2)
                        if(ibb(k2,1) == ibpoin(ikk1-1)+1  .or.  &
                             ibb(k2,1) == ibpoin(ikk1) .or.   &
                             ibb(k2,3) == -1) then
!               NO moving of final or initial node of profiles
                           ikr2 = 1
                           ilk = ibb(k2,2)
                           ilk1 = ibb(k2,3)
                        endif
                     endif

                     if(ibbk1 .gt. 0 .and.  ibbk2 .gt. 0) then
                        ilk = ibb(k1,2)
                        ilk1 = ibb(k1,3)
                        if(ilk .ne. ibb(k2,2) ) then
                           print *,'error jkoli3',ss
                           print *,x(k1),y(k1),ibb(k1,1),ibb(k1,2)
                           print *,x(k2),y(k2),ibb(k2,1),ibb(k2,2)
                           stop
                        endif
                     endif

                     if(ilk .gt. 0 ) then
                        if(  ibbk1 ==  ibpoin(ilk-1) + 1 .or.   &
                             ibbk1 == ibpoin(ilk)) then
!     the point k1 we can not remove or move
                           ikr1 = 1
                        endif
                        if(  ibbk2 ==  ibpoin(ilk-1) + 1 .or.   &
                             ibbk2 == ibpoin(ilk)) then
!     the point k2 we can not remove or move
                           ikr2 = 1
                        endif
                     endif


                     rlen = ((x2-x(lbn(kb0,1)))**2 +   &
                          (y2-y(lbn(kb0,1)))**2)
                     det = x(lbn(kb0,1))*(y1-y2) + x1*(y2-y(lbn(kb0,1)))  &
                          + x2*(y(lbn(kb0,1)) - y1)

                     ve1x = x(lbn(kb0,1)) - x1
                     ve1y = y(lbn(kb0,1)) - y1
                     ve2x = x2 - x1
                     ve2y = y2 - y1

                     ss = (ve1x*ve2x + ve1y*ve2y)/  &
                          (ve1x**2 + ve1y**2)**0.5/  &
                          (ve2x**2 + ve2y**2)**0.5

                     if(ss .gt. -0.96 ) then
!     the point k1 we can not remove or move
                        ikr1 = 1
                     endif
                     rlen = ((x1-x(lbn(kb2,2)))**2 +   &
                          (y1-y(lbn(kb2,2)))**2)
                     det = x(lbn(kb2,2))*(y1-y2) + x1*(y2-y(lbn(kb2,2)))  &
                          + x2*(y(lbn(kb2,2))-y1)
                     ve1x = x(lbn(kb2,2)) - x2
                     ve1y = y(lbn(kb2,2)) - y2
                     ve2x = x1 - x2
                     ve2y = y1 - y2
                     ss = (ve1x*ve2x + ve1y*ve2y)/  &
                          (ve1x**2 + ve1y**2)**0.5/  &
                          (ve2x**2 + ve2y**2)**0.5

                     if(ss .gt. -0.96 ) then
!     the point k2 we can not remove or move
                        ikr2 = 1
                     endif


                     if( ikr2 == 0 ) then
                        if(ikr1 == 0) then
                           x0 = (x(k1) + x(k2) )/2
                           y0 = (y(k1) + y(k2) )/2
                           w0 = (wp(k1,1) + wp(k2,1) )/2
                           rgai = (rga(k1) + rga(k2)) /2
                           rgbi = (rgb(k1) + rgb(k2)) /2
                           rgci = (rgc(k1) + rgc(k2)) /2
                           if(ibb(k1,1).gt.0 .and. ibb(k2,1).gt.0) then
                              if(ibb(k2,1) .gt. ibb(k1,1) ) then
                                 if(ibb(k2,1) - ibb(k1,1) .ge. 2) then
                                    il0 = int(ibb(k1,1) + ibb(k2,1))/2
                                    x0 = xb(il0)
                                    y0 = yb(il0)
                                    ibbk1 = il0
                                    ibbk111 = ilk1
                                 else
!     a few points
                                    x0 = x(k1)
                                    y0 = y(k1)
                                 endif
                              else
                           
                                 if(abs(ibb(k1,1)-ibb(k2,1)-  &
                                      ibpoin(ilk)+ibpoin(ilk-1)+1)  &
                                      .ge. 2) then
                                    il0=int(ibb(k1,1)+ibb(k2,1)+  &
                                       ibpoin(ilk)-ibpoin(ilk-1)+1)/2
                                    if(il0 .ge. ibpoin(ilk))   &
                                         il0=il0-ibpoin(ilk)+  &
                                         ibpoin(ilk-1)+1
                                    x0 = xb(il0)
                                    y0 = yb(il0)
                                    ibbk1 = il0
                                    ibbk11 = ilk
                                    ibbk111 = ilk1
                                 else
!     a few points
                                    x0 = x(k1)
                                    y0 = y(k1)
                                 endif
                              endif
                           endif
                        else
                           x0 = x(k1)
                           y0 = y(k1)
                           w0 = wp(k1,1)
                           rgai = rga(k1) 
                           rgbi = rgb(k1) 
                           rgci = rgc(k1) 
                        endif
                     else
                        if(ikr1 == 0) then
                           x0 = x(k2)
                           y0 = y(k2)
                           w0 = wp(k2,1)
                           rgai = rga(k2) 
                           rgbi = rgb(k2) 
                           rgci = rgc(k2) 
                           ibbk1 = ibb(k2,1)
                           ibbk11 = ibb(k2,2)
                           ibbk111 = ibb(k2,3)
                        else
!     no removing
                           goto 28
                        endif
                     endif
                     rll = ((x(k1) - x(k2))*(x(k1) - x(k2)) +  &
                          (y(k1) - y(k2))*(y(k1) - y(k2)) )

!     test the violation of positivity

                     ilen1 = abs(icyc(k1,1)) - 1
                     do 801 ic = 2,ilen1 
                        ic1 = ic + 1
                        x1 = x(icyc(k1,ic+1))
                        y1 = y(icyc(k1,ic+1))
                        x2 = x(icyc(k1,ic1+1))
                        y2 = y(icyc(k1,ic1+1))
                        reps = AMA%pos*( (x0-x1)**2+(y0-y1)**2 +  &
                             (x1-x2)**2+(y1-y2)**2 +  &
                             (x0-x2)**2 + (y0-y2)**2 )
                        det = x1*(y2-y0) + x0*(y1-y2) + x2*(y0-y1)
                        if( det .le. reps) then
                           goto 28
                        endif
                        call POS1TEST_77(x0,y0,x1,y1,x2,y2,itet)
                        if(itet == 1) then
                           goto 28
                        endif


 801                 enddo

                     ilen2 = abs(icyc(k2,1)) - 1
                     do 802 ic = 1,ilen2-1
                        ic1 = ic + 1
                        x1 = x(icyc(k2,ic+1))
                        y1 = y(icyc(k2,ic+1))
                        x2 = x(icyc(k2,ic1+1))
                        y2 = y(icyc(k2,ic1+1))
                        reps = AMA%pos*( (x0-x1)**2+(y0-y1)**2 +  &
                             (x1-x2)**2+(y1-y2)**2 +  &
                             (x0-x2)**2 + (y0-y2)**2 )
                        det = x1*(y2-y0) + x0*(y1-y2) + x2*(y0-y1)
                        if( det .le. reps) then
!     violation of positivity
                           goto 28
                        endif
                        call POS1TEST_77(x0,y0,x1,y1,x2,y2,itet)
                        if(itet == 1) then
                           goto 28
                        endif
 802                 enddo

!     periodic boundary
                     if(ibp(k1,1) .gt. 0 .and. ibp(k2,1) .gt. 0) then
                        if(ibp(k1,2) .ne. ibp(k2,2) ) goto 28
                        if(ibp(k2,2) == 1) then
                           xperreal = AMA%xper(1,1)
                           yperreal = AMA%xper(1,2)
                        else
                           xperreal = AMA%xper(2,1)
                           yperreal = AMA%xper(2,2)
                        endif                  

!                        print *,'@@',k1,ibp(k1,1),ibp(k1,2),x(k1),y(k1)
!                        print *,'@@',k2,ibp(k2,1),ibp(k2,2),x(k2),y(k2)

                        ipe = -1
                        do 1000 iel =1,nelem
                           do 1001 jel=1,3
                              if(iae(iel,jel) .lt. 0 .and.  &
                                   lnd(iel,jel) == ibp(k2,1))goto 1002
 1001                      enddo
 1000                   enddo
 1002                   continue
                        je1 = mod(jel,3)+1
                        je2 = mod(je1,3)+1
                        ke1 = lnd(iel,jel)
                        ke2 = lnd(iel,je1)
                        ke3 = lnd(iel,je2)
                        if(abs(x(k2) + xperreal - x(ke1) ) .lt.   &
                             1E-05 .and.  &
                             abs(y(k2)+yperreal-y(ke1) ) .lt.   &
                             1E-05 ) then
                           imov = 1
                        elseif(abs(x(k2)-xperreal-x(ke1)) .lt.   &
                                1E-05 .and.  &
                                abs(y(k2)-yperreal-y(ke1)).lt.   &
                                1E-05 ) then
                           imov = -1
                        else
                           print *,'BAD in remove in periodical points'
                           print *,AMA%xper(1,1),AMA%xper(1,2),AMA%xper(2,1),AMA%xper(2,2)
                           print *,i,k2,ke1,ibp(k2,1),ibp(k2,2)
                           print *,x(k2),y(k2),xperreal
                           print *,x(ke1),y(ke1),yperreal
                           print *,abs(x(k2) + xperreal - x(ke1) ),  &
                                abs(y(k2) + yperreal - y(ke1) ),  &
                                abs(x(k2) - xperreal - x(ke1) ),  &
                                abs(y(k2) - yperreal - y(ke1) )
                           stop
                        endif
                        xe0 = x0 + imov*xperreal
                        ye0 = y0 + imov*yperreal
                        we0 = w0
                        

!     test the violation of positivity
                        ilen1 = abs(icyc(ke1,1)) - 1
                        do 811 ic = 2,ilen1 
                           ic1 = ic + 1
                           x1 = x(icyc(ke1,ic+1))
                           y1 = y(icyc(ke1,ic+1))
                           x2 = x(icyc(ke1,ic1+1))
                           y2 = y(icyc(ke1,ic1+1))
                           reps = AMA%pos*( (xe0-x1)**2+(ye0-y1)**2 +  &
                                (x1-x2)**2+(y1-y2)**2 +  &
                                (xe0-x2)**2 + (ye0-y2)**2 )
                           det = x1*(y2-ye0) + xe0*(y1-y2) + x2*(ye0-y1)
                           if( det .le. reps) then
                              goto 28
                           endif
                           call POS1TEST_77(xe0,ye0,x1,y1,x2,y2,itet)
                           if(itet == 1) then
                              goto 28
                           endif
 811                    enddo

                        ilen2 = abs(icyc(ke2,1)) - 1
                        do 812 ic = 1,ilen2-1
                           ic1 = ic + 1
                           x1 = x(icyc(ke2,ic+1))
                           y1 = y(icyc(ke2,ic+1))
                           x2 = x(icyc(ke2,ic1+1))
                           y2 = y(icyc(ke2,ic1+1))
                           reps = AMA%pos*( (xe0-x1)**2+(ye0-y1)**2 +  &
                                (x1-x2)**2+(y1-y2)**2 +  &
                                (xe0-x2)**2 + (ye0-y2)**2 )
                           det = x1*(y2-ye0) + xe0*(y1-y2) + x2*(ye0-y1)
                           if( det .le. reps) then
                              goto 28
                           endif
                           call POS1TEST_77(xe0,ye0,x1,y1,x2,y2,itet)
                           if(itet == 1) then
                              goto 28
                           endif
 812                    enddo
                     endif
 993                 continue

                     ia1 = iae(i,j2)
                     if(ia1 .gt. 0) then
                        do 611 il =1,3
                           if(iae(ia1,il) == i) ja1 = il
 611                    enddo
                     endif
                     ia2 = iae(i,j1)
                     if(ia2 .gt. 0) then
                        do 612 il =1,3
                           if(iae(ia2,il) == i) ja2 = il
 612                    enddo
                     endif
                     if(ia1 .lt. 0 .or. ia2 .lt. 0 ) then
!                        print *,'PRoblem in remove'
!                        print *,icha
                        goto 28
                     endif

                     iae(ia1,ja1) = ia2
                     
                     ja21 = mod(ja2,3) + 1
                     lnd(ia2,ja21) = k1
                     iae(ia2,ja2) = ia1

!     begin of shifting 2
                     ibb(k1,1) = ibbk1
                     ibb(k1,2) = ibbk11
                     ibb(k1,3) = ibbk111

                     if(ipe == 0 .or. ipe == -1) then
                        x(k1) = x0
                        y(k1) = y0
                        wp(k1,1) = w0
                     elseif(ipe == 1) then
                        x(k1) = xe0
                        y(k1) = ye0
                        wp(k1,1) = we0
                     else
                        print *,'bad number ipe=',ipe
                        stop
                     endif

                     rga(k1) = rgai
                     rgb(k1) = rgbi
                     rgc(k1) = rgci


                     if(ibp(k2,1) .gt. 0 .and. ibp(k1,1) .le. 0) then
!     in fact we remove k1 and k2 stay with new index k1
                        if(ikr1 ==1 .or. ikr2 == 0) then
                           print *,'MISHMATCH in REMOVE'
                           print *,ikr1,ikr2
                           stop
                        endif
                        ibp(k1,1) = ibp(k2,1)
                        ibp(k1,2) = ibp(k2,2)
                     endif
                     ibp(k2,1) = ibp(npoin,1)
                     ibp(k2,2) = ibp(npoin,2)
                     do 6645 ip=1,npoin
                        if(ibp(ip,1) == k2 .and.k1 .ne. npoin)then 
                           ibp(ip,1) = k1
                        endif
                        if(ibp(ip,1) == npoin .and.k2 .ne. npoin)   &
                             ibp(ip,1) = k2
 6645                enddo


                     x(k2) = x(npoin)
                     y(k2) = y(npoin)
                     ibb(k2,1) = ibb(npoin,1)
                     ibb(k2,2) = ibb(npoin,2)
                     ibb(k2,3) = ibb(npoin,3)
                     wp(k2,1) = wp(npoin,1)
                     rga(k2) = rga(npoin)
                     rgb(k2) = rgb(npoin)
                     rgc(k2) = rgc(npoin)

                     do 620 ie =1,nelem
                        do 630 kj =1,3
                           if(lnd(ie,kj) == k2 .and. k1 .ne. npoin)  &
                                lnd(ie,kj) = k1
                           if(lnd(ie,kj) == npoin .and. k2 .ne.npoin)  &
                                lnd(ie,kj) = k2
                           if(iae(ie,kj) == nelem .and. i .ne. nelem)  &
                                iae(ie,kj)=i
 630                    enddo
 620                 enddo
                     do kj=1,3
                        lnd(i,kj) = lnd(nelem,kj)
                        iae(i,kj) = iae(nelem,kj)
                     enddo


!     connection of two cykles 2
                     locyc(1) = icyc(k2,1) + 1
                     ill = abs(locyc(1))
                     if(abs(locyc(1)) .gt. maxdeg) then
                        print *,'2 - ERROR too long icyc'
                        stop
                     endif
                     do 705 ic=1,ill
                        locyc(ic+1) = icyc(k2,ic+1)
 705                 enddo
                     do 715 ic=1,abs(icyc(k1,1)) - 2
                        locyc(ill+ic+1) = icyc(k1,ic+3)
 715                 enddo
                     locyc(1) = locyc(1) - (abs(icyc(k1,1)) - 2)

                     if(abs(locyc(1)) .gt. maxdeg) then
                        print *,'3 - ERROR too long icyc'
                        stop
                     endif
                     do 723 ic=1,abs(locyc(1)) + 1
                        icyc(k1,ic) = locyc(ic)
 723                 enddo

                     itestit = 0
                     k3len = abs(icyc(k3,1))

                     if(icyc(k3,k3len+1) == k2) then
                        icyc(k3,1) = icyc(k3,1) - 1
                        itestit = 1
                        goto 729
                     endif

                     do 725 ic=1,k3len - 1
                        if(icyc(k3,ic+1) == k2) then
                           do 727 ic1=ic,k3len - 1
                              icyc(k3,ic1+1) = icyc(k3,ic1+2)
 727                       enddo
                           if(icyc(k3,1) .gt. 0) then
                              icyc(k3,1) = icyc(k3,1) - 1
                           else
                              icyc(k3,1) = icyc(k3,1) + 1
                           endif
                           itestit = 1
                           goto 729
                        endif
 725                 enddo
 729                 continue

                     if(itestit == 0) then
                        print *,'ERROR k1 does not found in icyc of k3'
                        print *,icha,k1,k2,k3
                        print *,x(k1),y(k1)
                        print *,x(k2),y(k2)
                        print *,x(k3),y(k3)
                        it = k3
                        write(*,'(10i5)') it,icyc(it,1),icyc(it,2),  &
                             icyc(it,3),icyc(it,4),icyc(it,5),  &
                             icyc(it,6),icyc(it,7),  &
                             icyc(it,8),icyc(it,9)

                        stop
                     endif
!     end of connection of two cykles 2

                     do 735 ie = 1,npoin
                        do 737 ic=1,abs(icyc(ie,1))
                           if(icyc(ie,ic+1) == k2) icyc(ie,ic+1) = k1
                           if(icyc(ie,ic+1) == npoin)   &
                                icyc(ie,ic+1) = k2
 737                    enddo
 735                 enddo

                     do  ic=1,abs(icyc(npoin,1)) + 1
                        icyc(k2,ic) = icyc(npoin,ic)
                     enddo

                     if(k1 .ne. npoin) then
                        lbn(kb2,1) = k1
                     else
                        lbn(kb2,1) = k2
                        lbn(kb0,2) = k2
                     endif
                     do 810 ib=1,nbelm
                        do 815 jb=1,2
                           if(lbn(ib,jb) == npoin .and. k2.ne.npoin)  &
                                lbn(ib,jb) = k2
 815                    enddo
                        if(itc(ib) == nelem) itc(ib) = i
 810                 enddo
                     do 820 ib=kb1,nbelm-1
                        lbn(ib,1) = lbn(ib+1,1)
                        lbn(ib,2) = lbn(ib+1,2)
                        ibc(ib) = ibc(ib+1)
                        itc(ib) = itc(ib+1)
 820                 enddo

!     end of shifting 2
                     npoin = npoin - 1
                     nelem = nelem - 1
                     nbelm = nbelm - 1

                     if( ipe == -1) then
                        if(ke1 == npoin+1) ke1 = k2
                        if(ke2 == npoin+1) ke2 = k2 
                        do 1110 iel =1,nelem
                           do 1111 jel=1,3
                              if(iae(iel,jel) .lt. 0 .and.  &
                                   lnd(iel,jel) == ke1)goto 1112
 1111                      enddo
 1110                   enddo
 1112                   continue
                        i = iel
                        j = jel
                        ipe = 1
                        goto 999
                     endif 
                     if(ipe == 1) then
                        ipe = 0
                        i = ipoc
                     endif
                     icha = icha + 1

                     if(AMA%ifig .ge. 0 .and. mod(icha,5) == 0 ) then
                        call PLOT1_77(melem,nelem,mpoin,npoin,x,y,lnd)
                        AMA%ifig = AMA%ifig + 1
                     endif

                     do i=1,npoin
                        if(ibb(i,1) .gt. 0) then
                           do j=i+1,npoin
                              if(ibb(i,1) == ibb(j,1)) then
                                 print *,'the same IBB,  icha =',icha
                                 print *,xb(ibb(i,1)),yb(ibb(i,1)),i,  &
                                      ibb(i,1)
                                 print *,xb(ibb(j,1)),yb(ibb(j,1)),j,  &
                                      ibb(j,1)
                              endif
                           enddo
                        endif
                     enddo

                  endif
               endif
               

               if( i .gt. nelem) goto 2000
28          enddo
            ipoc = ipoc + 1
         endif
20    enddo
2000  continue


      return
    end subroutine REMOVE_77
    

    function ACCUTE_I_77(x1,y1,x2,y2,x3,y3,iedge,par)
      !     ... compute the inner angles of a triangle [x1,y1], [x2,y2], [x3,y3]
      !     ... and return 1 + \sum_{j=1}^3 ( (\cos a_i)^-)^2
      real*8 x(3),y(3),vax,vay,vbx,vby,ccos,ang
      
      ACCUTE_I = 1.
      x(1) = x1
      x(2) = x2
      x(3) = x3
      y(1) = y1
      y(2) = y2
      y(3) = y3

!      print*,x(1),y(1)
!      print*,x(2),y(2)
!      print*,x(3),y(3)

      j0 = mod(iedge,3) + 1
      j1 = mod(j0,3) + 1
      j2 = mod(j1,3) + 1
      vax = x(j0) - x(j1)
      vay = y(j0) - y(j1)
      vbx = x(j2) - x(j1)
      vby = y(j2) - y(j1)
      ccos = (vax*vbx + vay*vby)
      if(ccos .lt. 0.D+0) then
         ang = ccos*ccos/(vax*vax + vay*vay)/(vbx*vbx + vby*vby)
         ACCUTE_I = ACCUTE_I+12.*par*ang
      endif
!     print*
!     print*,vax,vay
!     print*,vbx,vby
!     print*, ccos,ang
!     print*,'--------------------'

!      print*,'ACCUTE_I = ',ACCUTE_I
!      print*,'***************************************'
      return 
    end function ACCUTE_I_77

    function ACCUTE_77(x1,y1,x2,y2,x3,y3,par)
!     ... compute the inner angles of a triangle [x1,y1], [x2,y2], [x3,y3]
!     ... and return 1 + \sum_{j=1}^3 ( (\cos a_i)^-)^2
      real*8 x(3),y(3),vax,vay,vbx,vby,ccos,ang

      ACCUTE = 1.
      x(1) = x1
      x(2) = x2
      x(3) = x3
      y(1) = y1
      y(2) = y2
      y(3) = y3

!      print*,x(1),y(1)
!      print*,x(2),y(2)
!      print*,x(3),y(3)

      do j0=1,3
         j1 = mod(j0,3) + 1
         j2 = mod(j1,3) + 1
         vax = x(j0) - x(j1)
         vay = y(j0) - y(j1)
         vbx = x(j2) - x(j1)
         vby = y(j2) - y(j1)
         ccos = (vax*vbx + vay*vby)
         if(ccos .lt. 0.D+0) then
            ang = ccos*ccos/(vax*vax + vay*vay)/(vbx*vbx + vby*vby)
            ACCUTE = ACCUTE+12.*par*ang
         endif
!         print*
!         print*,vax,vay
!         print*,vbx,vby
!         print*, ccos,ang
!         print*,'--------------------'
      enddo
!      print*,'ACCUTE = ',ACCUTE
!      print*,'***************************************'
      return 
    end function ACCUTE_77

    subroutine DELAUNAY_77(ndim, melem,nelem,mpoin,npoin,  &
         x,y,lnd,iae,icha,  &
         ra,wp,rga,rgb,rgc,iba,nbelm,mbelm,lbn,ibc,itc,  &
         maxdeg,icyc,ibp)
      dimension lnd(melem,3),iae(melem,3),x(mpoin),y(mpoin),  &
           lbn(mbelm,2),ibc(mbelm),itc(mbelm),  &
           ra(melem*3),wp(mpoin,ndim+1),   &
           icyc(mpoin,maxdeg),ibp(mpoin,2),  &
           rga(mpoin),rgb(mpoin),rgc(mpoin)
      integer iba(melem,3)
      real *8 x1,y1,x2,y2,x3,y3,x4,y4,reps123,reps134,det123,det134,  &
           detdel, z1, z2, z3, z4


      icha = 0
!      return

      do i=1,nelem
         do j=1,3
            iba(i,j) = 0
         enddo
      enddo

      itest = -1499

      ice = 0
      do 10 i=1,nelem
         do 20 j=1,3
            if(iae(i,j) .gt. 0 .and. iba(i,j) == 0) then
               j1 = mod(j,3) +1
               j2 = mod(j1,3) +1
               ii = iae(i,j)

               do 30 jjj=1,3
                  if(lnd(ii,jjj) == lnd(i,j1)) jj = jjj
 30            enddo
               iba(i,j) = 1
               iba(ii,jj) = 1
               jj1 = mod(jj,3) +1
               j0 = mod(jj1,3) +1
               k1 = lnd(i,j2)
               k2 = lnd(i,j)
               k3 = lnd(ii,j0)
               k4 = lnd(i,j1)
               if( k2 .ne. lnd(ii,jj1) .or. k4 .ne. lnd(ii,jj) ) then
                  print *,'ERROR !@#'
                  print *,k1,k2,k3,k4
                  return
               endif
               x1 = x(k1)
               y1 = y(k1)
               x2 = x(k2)
               y2 = y(k2)
               x3 = x(k3)
               y3 = y(k3)
               x4 = x(k4)
               y4 = y(k4)

               if( i == itest) then
!                  print *,x(lnd(i,j)),y(lnd(i,j)),j,ii,iba(i,j)
!                  print *,x(lnd(i,j1)),y(lnd(i,j1))
!                  print *,x(lnd(i,j2)),y(lnd(i,j2))
                  print *
                  print *,x1,y1
                  print *,x2,y2
                  print *,x3,y3
                  print *,x4,y4
               endif

!     we prohibid SWAPPING in case, where can appear alement with two
!     boundary segment
               if( (iae(i,j2) .lt. 0 .and. iae(ii,jj1) .lt. 0 ) .or.  &
                    (iae(i,j1) .lt. 0 .and. iae(ii,j0) .lt. 0) ) then
                  goto 20
               endif
!     we must still chech the orientation, i.e. the  cykles of points
!     k1,k2,k4 and k1,k4,k3 must have positive orientation

               det123 = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
               det134 = x1*(y3-y4) + x3*(y4-y1) + x4*(y1-y3) 

               reps123 = AMA%pos*( ((x1-x3)**2 + (y1-y3)**2) +  &
                    ((x2-x3)**2 + (y2-y3)**2) +  &
                    ((x1-x2)**2 + (y1-y2)**2) )
               reps134 = AMA%pos*( ((x1-x3)**2 + (y1-y3)**2) +  &
                    ((x4-x3)**2 + (y4-y3)**2) +  &
                    ((x1-x4)**2 + (y1-y4)**2) )
               if( det123 .le. reps123 .or. det134 .le. reps134) then
!                  print *,'violation of positivity-1'
                  goto 20
               endif
               call POS1TEST_77(x(k1),y(k1),x(k2),y(k2),x(k3),y(k3),itet)
               if(itet == 1) then
!                  print *,'violation of positivity-2'
                  goto 20
               endif
               call POS1TEST_77(x(k1),y(k1),x(k3),y(k3),x(k4),y(k4),itet)
               if(itet == 1) then
!                  print *,'violation of positivity-3'
                  goto 20
               endif


               z1 = x1*x1 + y1*y1
               z2 = x2*x2 + y2*y2
               z3 = x3*x3 + y3*y3
               z4 = x4*x4 + y4*y4

               detdel = (x1*(y2*z3-y3*z2) - x2*(y1*z3-y3*z1) +  &
                    x3*(y1*z2-y2*z1))   &
                    -  (x1*(y2*z4-y4*z2) - x2*(y1*z4-y4*z1) +  &
                    x4*(y1*z2-y2*z1))   &
                    +  (x1*(y3*z4-y4*z3) - x3*(y1*z4-y4*z1) +  &
                    x4*(y1*z3-y3*z1))   &
                    -  (x2*(y3*z4-y4*z3) - x3*(y2*z4-y4*z2) +  &
                    x4*(y2*z3-y3*z2)) 

               if( i == itest) print *, detdel

!               if(detdel .lt. -1E-05) then
               if(detdel .lt. -1.D-12) then


!                  print *,x1,y1,detdel
!                  print *,x2,y2
!                  print *,x3,y3
!                  print *,x4,y4
!                  call PLOT_77(melem,nelem,mpoin,npoin,x,y,lnd)
!
!                  stop


!     we swap the diagonal
                  itci = 0
                  itcii = 0

                  iyii = iae(ii,jj1)
                  if(iyii .gt. 0) then
                     do 62 iy=1,3
                        if( iae(iyii,iy) == ii) then
                           iyiij = iy
                        endif
 62                  enddo
                  else
!     boundary segment, we seek itc
                     do ib1 =1,nbelm
                        if(itc(ib1) == ii) then
                           itci = ib1
                           goto 662
                        endif
                     enddo
                     print *,'boundary segment in DELANAY for itc not'
                     stop
                  endif
 662              continue

                  iyi = iae(i,j1)
                  if(iyi .gt. 0) then
                     do 63 iy=1,3
                        if( iae(iyi,iy) == i) then
                           iyij = iy
                        endif
 63                  enddo
                  else
                     do ib1 =1,nbelm
                        if(itc(ib1) == i) then
                           itcii = ib1
                           goto 663
                        endif
                     enddo
                     print *,'boundary segment in DELANAY for itc not2'
                     stop
                  endif
 663              continue
                  if(itci .gt. 0) itc(itci) = i
                  if(itcii .gt. 0) itc(itcii) = ii
                  
                  lnd(i,j1) = k3
                  lnd(ii,jj1) = k1
                  
                  iae(i,j) = iyii
                  iae(i,j1) = ii
                  if(iyii .gt.0 ) iae(iyii,iyiij) = i

                  iae(ii,jj) = iyi
                  iae(ii,jj1) = i
                  if(iyi .gt.0 ) iae(iyi,iyij) = ii

                  icha = icha+1
               endif
            endif
 20      enddo
 10   enddo
      return
    end subroutine DELAUNAY_77



    subroutine SWAPPING_77(ndim, melem,nelem,mpoin,npoin,  &
         x,y,lnd,iae,icha,  &
         ra,wp,rga,rgb,rgc,iba,nbelm,mbelm,lbn,ibc,itc,  &
         maxdeg,icyc,ibp)
      dimension lnd(melem,3),iae(melem,3),x(mpoin),y(mpoin),  &
           lbn(mbelm,2),ibc(mbelm),itc(mbelm),  &
           ra(melem*3),wp(mpoin,ndim+1),   &
           icyc(mpoin,maxdeg),ibp(mpoin,2),  &
           rga(mpoin),rgb(mpoin),rgc(mpoin)
      integer iba(melem,3)
      real *8 x1,y1,x2,y2,x3,y3,x4,y4,reps123,reps134,det123,det134

      do i=1,nelem
         do j=1,3
            iba(i,j) = 0
         enddo
      enddo

      ice = 0
      icha = 0
      do 10 i=1,nelem
         do 20 j=1,3

            if(iae(i,j) .gt. 0 .and. iba(i,j) == 0) then
               j1 = mod(j,3) +1
               j2 = mod(j1,3) +1
               ii = iae(i,j)
               do 30 jjj=1,3
                  if(lnd(ii,jjj) == lnd(i,j1)) jj = jjj
 30            enddo
               iba(i,j) = 1
               iba(ii,jj) = 1
               jj1 = mod(jj,3) +1
               j0 = mod(jj1,3) +1
               k1 = lnd(i,j2)
               k2 = lnd(i,j)
               k3 = lnd(ii,j0)
               k4 = lnd(i,j1)
               if( k2 .ne. lnd(ii,jj1) .or. k4 .ne. lnd(ii,jj) ) then
                  print *,'ERROR !@#'
                  print *,k1,k2,k3,k4
                  return
               endif
               x1 = x(k1)
               y1 = y(k1)
               x2 = x(k2)
               y2 = y(k2)
               x3 = x(k3)
               y3 = y(k3)
               x4 = x(k4)
               y4 = y(k4)
               w1 = wp(k1,1)
               w2 = wp(k2,1)
               w3 = wp(k3,1)
               w4 = wp(k4,1)

!     we prohibid SWAPPING in case, where can appear alement with two
!     boundary segment
               if( (iae(i,j2) .lt. 0 .and. iae(ii,jj1) .lt. 0 ) .or.  &
                    (iae(i,j1) .lt. 0 .and. iae(ii,j0) .lt. 0) ) then
                  goto 20
               endif
!     we must still chech the orientation, i.e. the  cykles of points
!     k1,k2,k4 and k1,k4,k3 must have positive orientation

               det123 = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
               det134 = x1*(y3-y4) + x3*(y4-y1) + x4*(y1-y3) 

               reps123 = AMA%pos*( ((x1-x3)**2 + (y1-y3)**2) +  &
                    ((x2-x3)**2 + (y2-y3)**2) +  &
                    ((x1-x2)**2 + (y1-y2)**2) )
               reps134 = AMA%pos*( ((x1-x3)**2 + (y1-y3)**2) +  &
                    ((x4-x3)**2 + (y4-y3)**2) +  &
                    ((x1-x4)**2 + (y1-y4)**2) )
               if( det123 .le. reps123 .or. det134 .le. reps134) then
!                  print *,'violation of positivity'
                  goto 20
               endif

               itet = 0
               call POS1TEST_77(x(k1),y(k1),x(k2),y(k2),x(k3),y(k3),itet)
               if(itet == 1) then
                  goto 20
               endif
               call POS1TEST_77(x(k1),y(k1),x(k3),y(k3),x(k4),y(k4),itet)
               if(itet == 1) then
                  goto 20
               endif

               acc_old = ACCUTE_77(x(k1),y(k1),x(k2),y(k2),x(k4),y(k4),1.)  &
                    *ACCUTE_77(x(k2),y(k2),x(k3),y(k3),x(k4),y(k4),1. )

               acc_new = ACCUTE_77(x(k1),y(k1),x(k2),y(k2),x(k3),y(k3),1. )  &
                    *ACCUTE_77(x(k1),y(k1),x(k3),y(k3),x(k4),y(k4),1. )
               ice = 0

               epsround = 5E-03

               rl13 = ( (rga(k1)+rga(k3))*(x1-x3)*(x1-x3) +  &
                    2*(rgb(k1)+rgb(k3))*(x1-x3)*(y1-y3)  +  &
                    (rgc(k1)+rgc(k3))*(y1-y3)*(y1-y3) )/2

               rl24 = ( (rga(k2)+rga(k4))*(x2-x4)*(x2-x4) +  &
                    2*(rgb(k2)+rgb(k4))*(x2-x4)*(y2-y4)  +  &
                    (rgc(k2)+rgc(k4))*(y2-y4)*(y2-y4) )/2


               if(abs(rl13 - 3.) .lt. 0.995 * abs(rl24 -3.) )then
!     ... NEW ACCUTE
!               if(abs(rl13 - 3.)*acc_new .lt. 
!     *              0.995 * abs(rl24 -3.)*acc_old )then

!        ... we use a SWAPPING
                  itci = 0
                  itcii = 0

                  iyii = iae(ii,jj1)
                  if(iyii .gt. 0) then
                     do 62 iy=1,3
                        if( iae(iyii,iy) == ii) then
                           iyiij = iy
                        endif
 62                  enddo
                  else
!     boundary segment, we seek itc
                     do ib1 =1,nbelm
                        if(itc(ib1) == ii) then
                           itci = ib1
                           goto 662
                        endif
                     enddo
                     print *,'boundary segment in SWAPPING for itc not'
                     stop
                  endif
 662              continue

                  iyi = iae(i,j1)
                  if(iyi .gt. 0) then
                     do 63 iy=1,3
                        if( iae(iyi,iy) == i) then
                           iyij = iy
                        endif
 63                  enddo
                  else
                     do ib1 =1,nbelm
                        if(itc(ib1) == i) then
                           itcii = ib1
                           goto 663
                        endif
                     enddo
                     print *,'boundary segment in SWAPPING for itc not2'
                     stop
                  endif
 663              continue
                  if(itci .gt. 0) itc(itci) = i
                  if(itcii .gt. 0) itc(itcii) = ii
                  
                  lnd(i,j1) = k3
                  lnd(ii,jj1) = k1
                  
                  iae(i,j) = iyii
                  iae(i,j1) = ii
                  if(iyii .gt.0 ) iae(iyii,iyiij) = i

                  iae(ii,jj) = iyi
                  iae(ii,jj1) = i
                  if(iyi .gt.0 ) iae(iyi,iyij) = ii
                  icha = icha+1

                     if(AMA%ifig .ge. 0 .and. mod(icha,5) == 0 ) then
                        call PLOT1_77(melem,nelem,mpoin,npoin,x,y,lnd)
                        AMA%ifig = AMA%ifig + 1
                     endif

               endif
            endif
 20      enddo
 10   enddo

      return
    end subroutine SWAPPING_77

    end module 
