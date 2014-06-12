! -*- f90 -*-
module plume

real,parameter :: PI=4.*atan(1.0)
real,parameter :: PI2=2*PI
real :: ES=2.0
real :: EF=1.0
real :: alpha1=0.055,alpha2=0.6,alpha3=0.055,alhpa4=0.5

contains
! Forced entrainment as per Yapa et al.
subroutine entrainmentF(u,ua,b,h,ddens,np,Q,np1) 

integer np
real u(np1,3),ua(np,3),b(np1),h(np1),ddens(np1)
real Q(np1)
real vmod,va,modv,ubar,pibdb,tbds,pib2,alpha,beta
real sin_phi,cos_phi,sin_theta,cos_theta
real Qs,Qx,Qy,Qz

beta=1./ES**2
do i=np+1,1,-1
  vmod=sqrt(u(i,1)**2+u(i,2)**2+u(i,3)**2)
  va=(u(i,1)*ua(i,1)+u(i,2)*ua(i,2)**2+u(i,3)*ua(i,3))/vmod
  modv=sqrt((u(i,1)-ua(i,1))**2+(u(i,2)-ua(i,2))**2+(u(i,3)-ua(i,3))**2)
  ubar=sqrt(u(i,1)**2+u(i,2)**2)+1e-10
  sin_phi1=sin_phi
  cos_phi1=cos_phi
  sin_theta1=sin_theta
  cos_theta1=cos_theta
  sin_phi=u(i,3)/vmod
  cos_phi=ubar/vmod
  sin_theta=min(u(i,2)/ubar,1.)
  cos_theta=min(u(i,1)/ubar,1.)
  if (i.gt.np) cycle
  pibdb=PI*b(i)*(b(i)-b(i+1))
  tbds=2*b(i)*h(i)
  pib2=0.5*PI*b(i)**2
  alpha=(0.08+7.67*sin_phi*ddens(i)*b(i)/modv**2)/(0.707+3.536*va/modv)
  Qs=PI2*alpha*b(i)*h(i)*modv
  Qx=abs(ua(i,1))*(pibdb*abs(cos_phi*cos_theta)+tbds*sqrt(1-cos_theta**2*cos_phi**2)+pib2*(cos_phi*cos_theta-cos_phi1*cos_theta1))
  Qy=abs(ua(i,2))*(pibdb*abs(cos_phi*sin_theta)+tbds*sqrt(1-sin_theta**2*cos_phi**2)+pib2*(cos_phi*sin_theta-cos_phi1*sin_theta1))
  Qz=abs(ua(i,3))*(pibdb*abs(sin_phi)+tbds*cos_phi+pib2*(sin_phi-sin_phi))
  Q(i)=Qs+EF*(Qx+Qy+Qz)
enddo

end subroutine

!Jirka(2004) formulation 
subroutine entrainmentJ(u,ua,b,h,ddens,np,Q,np1) 

integer np,i
real u(np1,3),ua(np,3),b(np1),h(np1),ddens(np1)
real Q(np1)
real vmod,va,modv,ubar,pib2,Frl2,Qp,Qw,Qt,u1,u2,u3
real sin_sigma,cos_sigma,sin_theta,cos_theta,cos_thetasig

do i=1,np
  vmod=sqrt(u(i,1)**2+u(i,2)**2+u(i,3)**2)
  va=sqrt(ua(i,1)*ua(i,1)+ua(i,2)*ua(i,2)**2+ua(i,3)*ua(i,3))
  u1=u(i,1)-ua(i,1)
  u2=u(i,3)-ua(i,2)
  u3=u(i,2)-ua(i,3)
  modv=sqrt(u1**2+u2**2+u3**2)
  ubar=sqrt(u1**2+u2**2)+1e-10
  
  sin_theta=u3/modv
  cos_theta=ubar/modv
  sin_sigma=min(u2/ubar,1.)
  cos_sigma=min(u1/ubar,1.)
  
!  print*,sin_theta,cos_theta,sin_sigma,cos_sigma
  pib2=2*PI*b(i)
  
  if (ddens(i).gt.0) then
    Frl2=abs(modv*modv/ddens(i)/b(i))
    Qp=alpha2*sin_theta/Frl2
  else
    Qp=0
  endif
  cos_thetasig=abs(cos_theta*cos_sigma)
  Qw=alpha3*va*cos_thetasig/(modv+va)
  Qt=alpha4*sqrt(1-cos_thetasig**2)*cos_thetasig
  Q(i)=2*PI*b(i)*h(i)*(modv*(alpha1+Qp+Qw)+va*Qt)
enddo

end subroutine

end module