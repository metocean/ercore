      program ERcore
      
      use shoreline
      use griddata
      use release
      use outgrid
      use ncout
      
      parameter(d16=1/6.,pi=3.141592653589793,pi2=2*pi)
      parameter(r2d=180./pi,d2r=1./r2d,arad=r2d/6367456.)
      
      real,allocatable :: kx(:,:),ky(:,:),srand(:)
      real*8,allocatable :: pxt(:),pyt(:)
      character(len=120) infile,windfile,curfile,wavefile,sstfile,depfile,bndfile,outfile,ncoutfile
      character(len=19) tstart,tend
      integer difftype,itout,incout,np,dims(2)
      real diffval,diffac,dt
      real wndx,wndy,curx,cury
      real wspunc,wdirunc,dcurunc,ccurunc,diffunc
      real*8 mfy,diffacx,diffacy
      real*8 ts,te,tt,xlims(2),ylims(2),dts,dt2
      logical bnd

      wavefile=''
      sstfile='' 
      windfile=''
      curfile=''
      depfile=''

      namelist/files/windfile,curfile,wavefile,sstfile,bndfile,outfile
      namelist/time/tstart,tend,dt,itout
      namelist/movers/windspd,winddir,curspd,curdir,difftype,diffval,hs,tp,sst
      namelist/uncertainty/iunc,wspunc,wdirunc,dcurunc,ccurunc,diffunc
      namelist/ncout/ncoutfile,incout,xlims,ylims,dims
      
      windspd=0.
      curspd=0.
      winddir=0.
      curdir=0.
      hs=0.
      tp=10.
      sst=20.
      iunc=.false.
      wspdunc=0.
      wdirunc=0
      ccurunc=0.
      dcurunc=0.
      diffunc=1.
      difftype=1     
 
      if (iargc().gt.0) then
        call getarg(1,infile)
      else
        infile='ercore.nml'
      endif
      
      write(*,*),'Starting ERcore'
      write(*,'(a,a)') 'Using namelist file: ',infile

      open(11,file=infile)
        read(11,nml=files)
        read(11,nml=time)
        read(11,nml=movers)
        call init_release(11)
        read(11,nml=ncout)
      close(11)
      
      winddir=d2r*(270.-winddir)      
      curx=curspd*(cos(d2r*(270.-curdir)))
      cury=curspd*(sin(d2r*(270.-curdir)))
      diffunc=sqrt(diffunc)      
      wdirunc=pi*wdirunc/180.

      icur=(curfile.ne.'')
      iwind=(windfile.ne.'')
      ihs=(wavefile.ne.'')
      isst=(sstfile.ne.'')
      idep=(depfile.ne.'')
      
      allocate(dg(5))

      if (bndfile.eq.'') then
        write(*,'(/,a)')'No boundary file '
        bnd=.false.
      else
        bnd=.true.
        write(*,'(/,a,a)')'Reading boundary file: ',bndfile
        call read_shoreline(bndfile)
        write(*,'(i,a,i,a)') size(slx),' boundary points in ',maxval(sli),' sections'
        write(*,'(a,4(x,g))'),'Boundary limits: ',minval(slx),maxval(slx),minval(sly),maxval(sly)
      endif
      
      if (iwind) then
        write(*,'(/,a,a)')'Reading wind grid file: ',windfile
        call init_grid(windfile,1)
        write(*,'(a,2i)') 'Grid dimensions :',dg(1)%nx,dg(1)%ny
        write(*,'(a,4(x,g))'),'Boundary limits: ',minval(dg(1)%xx),maxval(dg(1)%xx),minval(dg(1)%yy),maxval(dg(1)%yy)
      endif
      
      if (icur) then
        write(*,'(/,a,a)')'Reading current grid file: ',curfile
        call init_grid(curfile,2)
        write(*,'(a,2i)') 'Grid dimensions :',dg(2)%nx,dg(2)%ny
        write(*,'(a,4(x,g))'),'Boundary limits: ',minval(dg(2)%xx),maxval(dg(2)%xx),minval(dg(2)%yy),maxval(dg(2)%yy)
      endif

      if (ihs) then
        write(*,'(/,a,a)')'Reading wave grid file: ',wavefile
        call init_grid(wavefile,3)
        write(*,'(a,2i)') 'Grid dimensions :',dg(3)%nx,dg(3)%ny
        write(*,'(a,4(x,g))'),'Boundary limits:',minval(dg(3)%xx),maxval(dg(3)%xx),minval(dg(3)%yy),maxval(dg(3)%yy)
      endif

      if (isst) then
        write(*,'(/,a,a)')'Reading sst grid file: ',sstfile
        call init_grid(curfile,4)
        write(*,'(a,2i)') 'Grid dimensions :',dg(4)%nx,dg(4)%ny
        write(*,'(a,4(x,g))'),'Boundary limits:',minval(dg(4)%xx),maxval(dg(4)%xx),minval(dg(4)%yy),maxval(dg(4)%yy)
      endif

      if (idep) then
        write(*,'(/,a,a)')'Reading depth grid file: ',depfile
        call init_grid(curfile,5)
        write(*,'(a,2i)') 'Grid dimensions :',dg(5)%nx,dg(5)%ny
        write(*,'(a,4(x,g))'),'Boundary limits:',minval(dg(5)%xx),maxval(dg(5)%xx),minval(dg(5)%yy),maxval(dg(5)%yy)
      endif

      
      if (.not.iwind) write(*,'(a,g,a,g)'),'Constant wind: ',windspd,' m/s From[deg]: ',winddir
      if (.not.icur) write(*,'(a,g,a,g)'),'Constant current: ',curspd,' m/s From[deg]: ',curdir
      if (.not.iwind) write(*,'(a,g,a,g)'),'Constant wave Height[m]: ',hs,' m Period[s]: ',tp
      if (.not.icur) write(*,'(a,g,a)'),'Constant SST: ',sst,' C'
     

      if (.not.iwind.and..not.icur.and..not.ihs.and..not.isst)  timebase=udtime("0001-01-01 00:00:00")
 
      if (difftype.eq.1) then
        diffacx=2*sqrt(6*diffval*dt)
        diffacy=2*sqrt(6*diffval*dt)
      else if (difftype.eq.2) then
        diffacx=sqrt(2*diffval*dt)
        diffacy=sqrt(2*diffval*dt)
      endif
!      print*,diffacx,diffacy
      
      do ir=1,nr
        call dump_release(ir)
      enddo
      
       allocate(kx(4,npmax),ky(4,npmax),pxt(npmax),pyt(npmax))
      write(*,'(/,a,f,a)')'Diffusion ',diffval,' m^2/s'      
      
      if (outfile.ne.'')then
        open(20,file=outfile,status='replace')
        write(20,*)nr
        do ir=1,nr
          write(20,'(I8,2F12.5,2F12.1)')rel(ir)%np,rel(ir)%lon,rel(ir)%lat,rel(ir)%rstart,rel(ir)%rend
        enddo
      endif
      
      write(*,'(/,a,a)')'Initialialising output grid file: ',ncoutfile
      call init_outgrid(xlims(1),xlims(2),ylims(1),ylims(2),dims(1),dims(2))
      call init_ncout(ncoutfile,nr)
      write(*,'(a,2i)') 'Grid dimensions :',dims
      write(*,'(a,4(x,f9.4))'),'Boundary limits: ',xlims,ylims
      
      write(*,'(/,a,a)')'Starting simulation at: ',tstart
      ts=udtime(tstart)-timebase
      te=udtime(tend)-timebase
!      print*,ts,te

      it=0
      tt=ts
      dts=dt/86400.
      dt2=0.5*dt
      if (iwind) call read_slab(tt,1)
      if (icur) call read_slab(tt,2)
      mfy=arad

      do while (tt.lt.te)
        write(*,'(a,f12.5)') 'Timestep :',tt-ts
        call add_release(it,dt)
      
        if(mod(it,itout).eq.0)then
          call poutput(tt-ts,dts)
        endif
        if (incout.gt.0) then
        if(mod(it,incout).eq.0)then
          do ir=1,nr
            call calc_dens(rel(ir)%px,rel(ir)%py,rel(ir)%mass,rel(ir)%npr)
            call output_dens(tt,ir)
          enddo
        endif
        endif
        
        tt=tt+dts
        it=it+1
        
        if (iwind) call read_slab(tt,1)
        if (icur) call read_slab(tt,2)
        if (ihs) call read_slab(tt,3)
        if (isst) call read_slab(tt,4)
  
        do ir=1,nr
        
        np=rel(ir)%npr
        rel(ir)%age(1:np)=rel(ir)%age(1:np)+dt
        px=>rel(ir)%px
        py=>rel(ir)%py
        psc=>rel(ir)%psc
        
        allocate(srand(np),mfx(np))        
        mfx=arad/cos(d2r*py(1:np))
        
        if (curx.ne.0) then
          pxt(1:np)=px(1:np)+dt*curx*mfx
          pyt(1:np)=py(1:np)+dt*cury*mfy
        elseif (icur) then
          call setgrid(2)
          kx(1,1:np)=dt*interph(uuo,px,py,np)
          ky(1,1:np)=dt*interph(vvo,px,py,np)
          pxt(1:np)=px(1:np)+0.5*kx(1,1:np)*mfx
          pyt(1:np)=py(1:np)+0.5*ky(1,1:np)*mfy
          kx(2,1:np)=dt2*(interph(uuo,pxt,pyt,np)+interph(uu,pxt,pyt,np))
          ky(2,1:np)=dt2*(interph(vvo,pxt,pyt,np)+interph(vv,pxt,pyt,np))
          pxt(1:np)=px(1:np)+0.5*kx(2,1:np)*mfx
          pyt(1:np)=py(1:np)+0.5*ky(2,1:np)*mfy
          kx(3,1:np)=dt2*(interph(uuo,pxt,pyt,np)+interph(uu,pxt,pyt,np))
          ky(3,1:np)=dt2*(interph(vvo,pxt,pyt,np)+interph(vv,pxt,pyt,np))
          pxt(1:np)=px(1:np)+kx(3,1:np)*mfx
          pyt(1:np)=py(1:np)+ky(3,1:np)*mfy
          kx(4,1:np)=dt*interph(uu,pxt,pyt,np)
          ky(4,1:np)=dt*interph(vv,pxt,pyt,np)
          kx(1,1:np)=d16*(kx(1,1:np)+2*(kx(2,1:np)+kx(3,1:np))+kx(4,1:np))
          ky(1,1:np)=d16*(ky(1,1:np)+2*(ky(2,1:np)+ky(3,1:np))+ky(4,1:np))
          if (iunc) then
                call random_number(srand)
                kx(3,1:np)=1.+dcurunc*2*(srand-0.5)
                call random_number(srand)
                kx(4,1:np)=ccurunc*2*(srand-0.5)
                kx(2,1:np)=kx(1,1:np)*kx(3,1:np)+kx(4,1:np)*ky(1,1:np)
                ky(2,1:np)=ky(1,1:np)*kx(3,1:np)-kx(4,1:np)*kx(1,1:np)
                pxt(1:np)=px(1:np)+kx(2,1:np)*mfx
                pyt(1:np)=py(1:np)+ky(2,1:np)*mfy
          else
                pxt(1:np)=px(1:np)+kx(1,1:np)*mfx
                pyt(1:np)=py(1:np)+ky(1,1:np)*mfy
          endif
        else
          pxt(1:np)=px(1:np)
          pyt(1:np)=py(1:np)
        endif
        
        if (iwind.or.windspd.ne.0.) then
          if (iwind) then
            call setgrid(1)
            kx(1,1:np)=0.5*(interph(uuo,px,py,np)+interph(uu,pxt,pyt,np))
            ky(1,1:np)=0.5*(interph(vvo,px,py,np)+interph(vv,pxt,pyt,np))
            kx(2,1:np)=sqrt(kx(1,1:np)**2+ky(1,1:np)**2)
            kx(3,1:np)=atan2(ky(1,1:np),kx(1,1:np))
          else
            kx(2,1:np)=windspd
            kx(3,1:np)=winddir
          endif
          call random_number(srand)
          kx(4,1:np)=2*(rel(ir)%cw_max-rel(ir)%cw_min)*(srand-0.5)
          where (kx(4,1:np).lt.0)
            kx(4,1:np)=kx(4,1:np)-rel(ir)%cw_min
          elsewhere
            kx(4,1:np)=kx(4,1:np)+rel(ir)%cw_min
          end where
          kx(3,1:np)=kx(3,1:np)+asin(kx(4,1:np))
          call random_number(srand)
          kx(2,1:np)=kx(2,1:np)*(rel(ir)%dw_min+(rel(ir)%dw_max-rel(ir)%dw_min)*srand)
          if (iunc) then
                call random_number(srand)
                kx(2,1:np)=kx(2,1:np)*(1.+2*wspunc*(0.5-srand))
                call random_gauss(srand,np)
                kx(3,1:np)=kx(3,1:np)+wdirunc*srand
          endif
          pxt(1:np)=pxt(1:np)+dt*kx(2,1:np)*cos(kx(3,1:np))*mfx
          pyt(1:np)=pyt(1:np)+dt*kx(2,1:np)*sin(kx(3,1:np))*mfy
        endif
        
        if (difftype.eq.1) then
          call random_number(srand)
          kx(1,1:np)=diffacx*(srand-0.5)
          call random_number(srand)
          ky(1,1:np)=diffacy*(srand-0.5)
        else if (difftype.eq.2) then
          call random_gauss(srand,np)
          kx(1,1:np)=diffacx*srand
          call random_gauss(srand,np)
          ky(1,1:np)=diffacy*srand
        endif
        if (difftype.gt.0) then
          if (iunc) then
                kx(1,1:np)=diffunc*kx(1,1:np)
                ky(1,1:np)=diffunc*ky(1,1:np)
          endif
          pxt(1:np)=pxt(1:np)+kx(1,1:np)*mfx
          pyt(1:np)=pyt(1:np)+ky(1,1:np)*mfy
        endif
        
        call sl_intersect(px,py,pxt(1:np),pyt(1:np),psc,rel(ir)%refloat,np)

        rel(ir)%px(1:np)=pxt
        rel(ir)%py(1:np)=pyt
 
        deallocate(srand)
        
        enddo
        
      enddo
      
      write(*,'(/,a,a)')'Ended simulation at: ',tend
      
      close(20)
      call close_ncout()      
      
      end


      subroutine random_gauss(rrand,np)
        parameter(pi=3.141592653589793,pi2=2*pi)
        real rrand(np),rand0(np),rand1(np)

        call random_number(rand0)
        call random_number(rand1)

        rrand=sqrt(-2*log(rand0))*cos(pi2*rand1)
      end subroutine
