! -*- f90 -*-
!
      module shoreline

      parameter(tol=0.000000001)

      real*8,allocatable,dimension(:) :: slx,sly
      real,allocatable,dimension(:,:) :: slbnd
      integer,allocatable,dimension(:) :: sli,slt,polyi,polyn
      integer nt
      real mapbndx(4),mapbndy(4),maxbndx,maxbndy,minbndx,minbndy

      contains
      subroutine read_shoreline(file)

      integer pass,nn,ns
      character(len=120) file
      character(len=12) h1,h2

      minbndx=-180.
      maxbndx=360.
      minbndy=-90.
      maxbndy=90.

      if (allocated(slx)) then
        deallocate(slx,sly,sli,slt)
        deallocate(slbnd,polyi,polyn)
      endif

      open(10,file=file,status='old')
      do pass=1,2
        nn=0
        nt=0
        do while (.true.)
          read(10,*,end=100) h1,h2,ns
          if (h1(1:3).eq.'Map') then
            do i=1,4
              read(10,*)mapbndx(i),mapbndy(i)
            enddo
            minbndx=minval(mapbndx)
            maxbndx=maxval(mapbndx)
            minbndy=minval(mapbndy)
            maxbndy=maxval(mapbndy)
            cycle
          endif
          nt=nt+1
          if (pass.eq.2) then
            polyi(nt)=nn+1
            polyn(nt)=ns
          endif
          do i=1,ns
            nn=nn+1
            read(10,*)xdum,ydum
            if (pass.eq.2) then
              slx(nn)=xdum
              sly(nn)=ydum
              sli(nn)=nt
              slt(nn)=1
              slbnd(nt,1)=min(slbnd(nt,1),xdum)-tol
              slbnd(nt,2)=max(slbnd(nt,2),xdum)+tol
              slbnd(nt,3)=min(slbnd(nt,3),ydum)-tol
              slbnd(nt,4)=max(slbnd(nt,4),ydum)+tol
            endif
          enddo
        enddo
100     continue
        if (pass.eq.1) then
          allocate (slx(nn),sly(nn),sli(nn),slt(nn))
          allocate (slbnd(nt,4),polyi(nt),polyn(nt))
          do ibnd=1,nt
            slbnd(ibnd,:)=(/360.,0.,90.,-90./)
          enddo
          rewind(10)
        endif
      enddo
      close(10)
      end subroutine


      subroutine gridinpoly(mask,xx,yy,nx,ny)
        integer nx,ny,ipoly,iseg
        real mask(nx,ny)
        real*8 xx(nx),yy(ny),x1,x2,y1,y2,x,y
        logical test

        mask=1.0
        do ix=1,nx
          do iy=1,ny
            test=.false.
            do ipoly=1,nt
              if (test) cycle
              if (xx(ix).lt.slbnd(ipoly,1)) cycle
              if (xx(ix).gt.slbnd(ipoly,2)) cycle
              if (yy(iy).lt.slbnd(ipoly,3)) cycle
              if (yy(iy).gt.slbnd(ipoly,4)) cycle
              do iseg=polyi(ipoly),polyi(ipoly)+polyn(ipoly)-2
                x1=slx(iseg)
                y1=sly(iseg)
                x2=slx(iseg+1)
                y2=sly(iseg+1)
                x=xx(ix)
                y=yy(iy)
                if ((y2.gt.y).neqv.(y1.gt.y)) then
                  if (x.lt.((x1-x2)*(y-y2)/(y1-y2)+x2)) test=.not.test
                endif
              enddo
              if (test) then
                mask(ix,iy)=0.0
                cycle
              endif
            enddo
          enddo
        enddo
      end subroutine


      subroutine sl_intersect_wrap(pos,post,poso,psc,refloat,np)

      real*8 :: pos(np,2),post(np,2),poso(np,2),pxt(np),pyt(np)
      integer :: psc(np)

      pxt=post(:,1)
      pyt=post(:,2)
      call sl_intersect(pos(:,1),pos(:,2),pxt,pyt,psc,refloat,np)
      poso(:,1)=pxt
      poso(:,2)=pyt

      end subroutine


      subroutine sl_intersect(px,py,pxt,pyt,psc,refloat,np)

      parameter(tol=0.000000001)
      integer np,iseg,ip
      real*8 :: px(np),py(np),pxt(np),pyt(np)
      integer :: psc(np)
      real*8 a1,b1,c1,a2,b2,c2,det,x1,x2,y1,y2,xi,yi,dsx,dsx0
      real refloat
      ip=1

      parloop: do while (ip.le.np)
! Check for out of bounds
        if (psc(ip).lt.-1) then
            pxt(ip)=px(ip)
            pyt(ip)=py(ip)
            ip=ip+1
            cycle parloop
        elseif(pxt(ip).lt.minbndx.or.pxt(ip).gt.maxbndx.or.pyt(ip).lt.minbndy.or.pyt(ip).gt.maxbndy) then
            print*, 'aqui1 - stickou a shoreline = -9'
            px(ip)=pxt(ip)
            py(ip)=pyt(ip)
            psc(ip)=-9
            cycle parloop
        endif
! Check for refloat
        if (psc(ip).gt.1) then
          if (refloat.gt.0) then
            psc(ip)=0
          else
            psc(ip)=psc(ip)+1
            pxt(ip)=px(ip)
            pyt(ip)=py(ip)
            ip=ip+1
            cycle parloop
          endif
        endif
        dsx=1.e20
        do ipoly=1,nt
          if (pxt(ip).lt.slbnd(ipoly,1)) cycle
          if (pxt(ip).gt.slbnd(ipoly,2)) cycle
          if (pyt(ip).lt.slbnd(ipoly,3)) cycle
          if (pyt(ip).gt.slbnd(ipoly,4)) cycle
          do iseg=polyi(ipoly),polyi(ipoly)+polyn(ipoly)-2
            x1=slx(iseg)
            y1=sly(iseg)
            x2=slx(iseg+1)
            y2=sly(iseg+1)
            if (min(px(ip),pxt(ip)).gt.max(x1,x2)) cycle
            if (max(px(ip),pxt(ip)).lt.min(x1,x2)) cycle
            if (min(py(ip),pyt(ip)).gt.max(y1,y2)) cycle
            if (max(py(ip),pyt(ip)).lt.min(y1,y2)) cycle
            a1=x1-x2
            b1=y1-y2
            a2=px(ip)-pxt(ip)
            b2=py(ip)-pyt(ip)
            det=a1*b2-a2*b1
            if(det.ne.0) then
              c1=x1*y2-y1*x2
              c2=px(ip)*pyt(ip)-pxt(ip)*py(ip)
              xi=(c1*a2-a1*c2)/det
              yi=(c1*b2-b1*c2)/det
!              write(*,'(4(f12.7,f12.7,a),f12.7,f12.7)') x1,y1,'->',x2,y2,'  ',px(ip),py(ip),'->',pxt(ip),pyt(ip),':',xi,yi
!              print*,ip,iseg,a1,b1,a2,b2,det
!              print*,(x1-xi),(xi-x2),(px(ip)-xi),(xi-pxt(ip)),tol
              if ((x1-xi)*(x2-xi).lt.tol.and.(px(ip)-xi)*(pxt(ip)-xi).lt.tol) then
                if ((y1-yi)*(y2-yi).lt.tol.and.(py(ip)-yi)*(pyt(ip)-yi).lt.tol) then
                  if (psc(ip).eq.0.and.det.gt.0.) then
                    psc(ip)=1
                    cycle parloop
                  else
                    dsx0=sqrt((xi-px(ip))**2+(yi-py(ip))**2)
                    if (dsx0.lt.dsx) then
                       psc(ip)=2
                       pxt(ip)=xi
                       pyt(ip)=yi
                       dsx=dsx0
                    endif
                  endif
                endif
              endif
!           else
!              print*,'Determinant 0 ',ip
            endif
          enddo
        enddo
        ip=ip+1
      enddo parloop
      end subroutine



      end module shoreline
