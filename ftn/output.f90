! -*- f90 -*-
      module outgrid
      
      integer ognx,ogny
      real*8 ogdx,ogdy
      
      real*8,allocatable,dimension(:) :: xlon,ylat,gmfx
      real,allocatable,dimension(:,:) :: cconc,mask
      
      parameter(hfra=0.1,gmfy=111133.)
    
      contains
      
      subroutine init_outgrid(x1,x2,y1,y2,nx,ny)
        use shoreline
        real*8 x1,x2,y1,y2
        integer nx,ny  
        ognx=nx
        ogny=ny
        allocate(xlon(nx),ylat(ny),cconc(nx,ny),mask(nx,ny),gmfx(ny))
        cconc=0.
        ogdx=(x2-x1)/ognx
        ogdy=(y2-y1)/ogny
        do i=1,nx
          xlon(i)=ogdx*(i-1)+x1
        enddo
        do i=1,ny
          ylat(i)=ogdy*(i-1)+y1
          gmfx(i)=gmfy*cos(0.01745329*ylat(i))
        enddo
        call gridinpoly(mask,xlon,ylat,ognx,ogny)
      end subroutine
      
      subroutine calc_dens(px,py,mass,np)
        use m_rnkpar
        integer np,ix,iy,inh,ind(np),ndx,ndy
        real*8 px(np),py(np),dxp(np),dyp(np),r(np),mass(np)
        real*8 dx,dy,xmax,xmin,ymax,ymin,lr,dxr,dyr,df
        real sum1,sum2

        cconc=0.
        xmin=minval(px(1:np))
        xmax=maxval(px(1:np))
        ymin=minval(py(1:np))
        ymax=maxval(py(1:np))
        dx=xmax-xmin
        dy=ymax-ymin
        xmin=xmin-dx
        xmax=xmax+dx
        ymin=ymin-dy
        ymax=ymax+dy
        inh=int(hfra*np)+1
        do ix=1,ognx
          do iy=1,ogny
!            print*,ix,iy
            if (mask(ix,iy).eq.0.or.xlon(ix).lt.xmin.or.xlon(ix).gt.xmax.or.ylat(iy).lt.ymin.or.ylat(iy).gt.ymax) cycle
            dxp=abs(px(1:np)-xlon(ix))*gmfx(iy)
            dyp=abs(py(1:np)-ylat(iy))*gmfy
            r=sqrt(dxp*dxp+dyp*dyp)
            call D_rnkpar(r, ind, inh)
            lr=r(ind(inh))
!            lr=min(lr,2*std(r(ind(1:inh)),inh))
            do ip=1,inh
              if (r(ind(ip)).gt.lr) cycle
                dx=dxp(ind(ip))/lr
                dy=dyp(ind(ip))/lr    
                cconc(ix,iy)=cconc(ix,iy)+0.5625*mass(ind(ip))*(1-dx*dx)*(1-dy*dy)
            enddo
            cconc(ix,iy)=cconc(ix,iy)/(lr*lr)
            dx=gmfx(iy)*ogdx
            dy=gmfy*ogdy
            ndx=int(dx/lr)+1
            ndy=int(dy/lr)+1
            sum1=0.
            sum2=0.
            do i=max(1,ix-ndx),min(ognx,ix+ndx)
              do j=max(1,iy-ndy),min(ogny,iy+ndy)
                dxr=(dx*(i-ix))**2
                dyr=(dy*(j-iy))**2
                if (sqrt(dxr+dyr).gt.lr) cycle
                  df=(1.-dxr)*(1.-dyr)
                  sum1=sum1+0.5625*df
                  sum2=sum2+0.5625*df*mask(i,j)
              enddo
            enddo
            cconc(ix,iy)=sum1*cconc(ix,iy)/sum2
          enddo
        enddo
        print*, 'calc_dens cconc(1,1) = ', cconc(1,1)
      end subroutine
      
      
      end module outgrid
