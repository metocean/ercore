      module slipvel
      
      parameter(g=9.81,mu20=0.001002)
      
      contains
      
      function visco_water(temp)
      
      real temp,visco_water,a,b
      b=20-temp
      a=b*(1.2378-0.001303*b+3.06e-6*b**2+2.55e-8*b**3)/(temp+96)
      visco_water=0.001002*exp(a)
      
      end function
      
! Helper function calculate bubble slip in ellipsoid regime
! P=g*ddens/ift
! Q=(visco/visco_water(temp))**-0.14
! R=(g*visco**4*dd/dens/ift**3)**-0.149

      function ellip_slip(de,visco,dens,P,Q,R)
      
      real ellip_slip,de,visco,dens,P,Q,R
      real Eo,J,H
      
      Eo=P*de**2
      H=1.3333*Eo*R*Q**(-0.14)
      if (H.gt.59.3) then
        J=3.42*H**0.441
      else
        J=0.94*H**0.757
      endif
      ellip_slip=visco*R*(J-0.857)/de/dens
      
      end function
      
! Function to calculate bubble slip        
      subroutine bubble_slip(db,dens,ddens,temp,visco,ift,slip,np)
      
      integer ip,np
      real, dimension(np) :: db,dens,ddens,temp,slip
      real ift,visco
      real R,Nd,W,x1,y1,y2,a2,b1,b2,dd,v,Q,P,dc,dq
      
      do ip=1,np
      
      if (visco.le.0) then
        v=visco_water(temp(ip))
      else
        v=visco
      endif
      
      dd=abs(ddens(ip))/dens(ip)
      if (db(ip).le.0.001) then
        Nd=4*g*abs(ddens(ip))*dens(ip)*db(ip)**3/(3*v*v)
        if (Nd.le.73) then
          R=Nd/24-1.7569e-4*Nd**2+6.9252e-7*Nd**3-2.3027e-10*Nd**4
        else
          W=log10(Nd)
          if (Nd.le.580) then
            R=10**(-1.7095+1.33438*W-0.11592*W*W)
          else
            R=10**(-1.81391+1.34671*W-0.12427*W*W+0.006344*W**3)
          endif
        endif
        slip(ip)=sign(R*v/dens(ip)/db(ip),ddens(ip))
        !print*,Nd,db(ip),R,slip(ip)
      else
        if (visco.gt.0) then
          Q=(v/visco_water(temp(ip)))**(-0.14)
        else
          Q=1
        endif
        R=(g*v**4*dd/dens(ip)/ift**3)**(-0.149)
        P=g*abs(ddens(ip))/ift
        dq=sqrt(0.75*59.3/P/Q/R)
        y1=log10(ellip_slip(dq,v,dens(ip),P,Q,R))
        y1=v*R*19.8125/dens(ip)/dq
        y2=log10(ellip_slip(0.015,v,dens(ip),P,Q,R))
        b1=log10(0.711*(g*dd)*0.5)
        x1=log10(dq)
        a2=(y2-y1)/(log10(0.015)-x1)
        b2=y1-a2*x1
        dc=10**((b2-b1)/(0.5-a2))
        !print*,dq,dc,g
        if (db(ip).gt.dc) then
          slip(ip)=sign(0.711*(g*db(ip)*dd)**0.5,ddens(ip))
        else
          slip(ip)=sign(ellip_slip(db(ip),v,dens(ip),P,Q,R),ddens(ip))
        endif
      endif
      
      enddo
      
      return
      end subroutine
      
      end module slipvel
      
      
      
