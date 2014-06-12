! -*- f90 -*-
      subroutine interph(var,px,py,x0,y0,idx,idy,vout,np,nx0,ny0)
      integer np,ip,ii,jj,nx0,ny0
      real var(ny0,nx0),vout(np)
      real px(np),py(np),idx,idy
      real sdx,sdy,sdx1,sdy1,ppx,ppy,fac,tmpint,tmp

      do ip=1,np
        ppx=idx*(px(ip)-x0)+1.
        ppy=idy*(py(ip)-y0)+1.
!        print*,px(ip),py(ip),ppx,ppy
        tmpint=0.
        fac=0.
        sdx=mod(ppx,1.0)
        sdy=mod(ppy,1.0)
        sdx1=(1.-sdx)
        sdy1=(1.-sdy)
        ii=min(max(1,int(ppx)),nx0)
        jj=min(max(1,int(ppy)),ny0)
        if (abs(var(jj,ii)).lt.9999) then
          tmp=sdx1*sdy1
          tmpint=tmpint+tmp*var(jj,ii)
          fac=fac+tmp
        endif
        if (jj.lt.ny0.and.abs(var(jj+1,ii)).lt.9999) then
          tmp=sdx1*sdy
          tmpint=tmpint+tmp*var(jj+1,ii)
          fac=fac+tmp
        endif
        if (ii.lt.nx0.and.abs(var(jj,ii+1)).lt.9999) then
          tmp=sdx*sdy1
          tmpint=tmpint+tmp*var(jj,ii+1)
          fac=fac+tmp
        endif
        if (ii.lt.nx0.and.jj.lt.ny0.and.abs(var(jj+1,ii+1)).lt.9999) then
          tmp=sdx*sdy
          tmpint=tmpint+tmp*var(jj+1,ii+1)
          fac=fac+tmp
        endif
        if (fac.gt.0) then
          tmpint=tmpint/fac
        else
          tmpint=0.
        endif
        vout(ip)=tmpint
      enddo
      return
      end subroutine
      
      subroutine interp3D(var,px,py,pz,x0,y0,idx,idy,z,vout,np,nx0,ny0,nz0)
      integer np,ip,iz,ii,jj,kk,nx0,ny0,nz0
      real var(nz0,ny0,nx0),z(nz0),vout(np),idx,idy
      real*8 px(np),py(np),pz(np)
      real vfac,fac,tmpint
      logical zup
      
      zup=(nz0.gt.1.and.z(2).gt.z(1))
      do ip=1,np
        iz=0
        if (zup) then
          do while (iz.lt.nz0.and.pz(ip).gt.z(iz+1))
            iz=iz+1
          enddo
        else
          do while (iz.lt.nz0.and.pz(ip).lt.z(iz+1))
            iz=iz+1
          enddo
        endif
        ppx=idx*(px(ip)-x0)+1.
        ppy=idy*(py(ip)-y0)+1.
        vup=0.
        vdown=0.
        sdx=mod(ppx,1.0)
        sdy=mod(ppy,1.0)
        sdx1=(1.-sdx)
        sdy1=(1.-sdy)
        ii=min(max(1,int(ppx)),nx0)
        jj=min(max(1,int(ppy)),ny0)
        do kk=max(iz,1),min(iz+1,nz0)
          fac=0.
          tmpint=0.
          if (abs(var(kk,jj,ii)).lt.9999) then
            tmp=sdx1*sdy1
            tmpint=tmpint+tmp*var(kk,jj,ii)
            fac=fac+tmp
          endif
          if (jj.lt.ny0.and.abs(var(kk,jj+1,ii)).lt.9999) then
            tmp=sdx1*sdy
            tmpint=tmpint+tmp*var(kk,jj+1,ii)
            fac=fac+tmp
          endif
          if (ii.lt.nx0.and.abs(var(kk,jj,ii+1)).lt.9999) then
            tmp=sdx*sdy1
            tmpint=tmpint+tmp*var(kk,jj,ii+1)
            fac=fac+tmp
          endif
          if (ii.lt.nx0.and.jj.lt.ny0.and.abs(var(kk,jj+1,ii+1)).lt.9999) then
            tmp=sdx*sdy
            tmpint=tmpint+tmp*var(kk,jj+1,ii+1)
            fac=fac+tmp
          endif
          if (fac.gt.0) then
            tmpint=tmpint/fac
          else
            tmpint=1e20
          endif
          vdown=vup
          vup=tmpint
        enddo
        if (iz.eq.0.or.iz.eq.nz0.or.vdown.gt.1e10) then
          vout(ip)=vup
        else if (vup.gt.1e10) then
          vout(ip)=vdown
        else
          vfac=(pz(ip)-z(iz))/(z(iz+1)-z(iz))
          vout(ip)=vfac*vup+(1-vfac)*vdown
        endif
        if (vout(ip).gt.1e10) vout(ip)=0.
      enddo
      return
      end subroutine
      
      
      subroutine interpz(var,pz,z,vout,np,nz0)
      integer np,ip,iz,kk,nz0
      real var(nz0,np),z(nz0),vout(np)
      real*8 pz(np)
      real vfac
      logical zup
      
      zup=(nz0.gt.1.and.z(2).gt.z(1))
      do ip=1,np
        iz=0
        if (zup) then
          do while (iz.lt.nz0.and.pz(ip).gt.z(iz+1))
            iz=iz+1
          enddo
        else
          do while (iz.lt.nz0.and.pz(ip).lt.z(iz+1))
            iz=iz+1
          enddo
        endif
        vup=0.
        vdown=0.
        do kk=max(iz,1),min(iz+1,nz0)
          vdown=vup
          vup=var(kk,ip)
        enddo
        if (iz.eq.0.or.iz.eq.nz0.or.vdown.gt.1e10) then
          vout(ip)=vup
        else if (vup.gt.1e10) then
          vout(ip)=vdown
        else
          vfac=(pz(ip)-z(iz))/(z(iz+1)-z(iz))
          vout(ip)=vfac*vup+(1-vfac)*vdown
        endif
        if (vout(ip).gt.1e10) vout(ip)=0.
      enddo
      return
      end subroutine
      
      
      