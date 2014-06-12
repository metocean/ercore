      module griddata
      
      type gridstr
        real,pointer,dimension(:) :: xx,yyz
        real,pointer,dimension(:,:) :: uu0,vv0,uu1,vv1
        real,pointer,dimension(:,:) :: uu,vv,uuo,vvo
        real*8,pointer,dimension(:) :: ttime
        real x0,y0,dx,dy,idx,idy
        real*8 time0,time1
        integer ncfile,iu,iv
        integer nt,nx,ny
        integer timeind
      end type gridstr
      
      type(gridstr), target, allocatable :: dg(:)
      
      logical icur,iwind
      real*8 timebase,timebase0
      integer nx0,ny0
      real x0,y0,idx,idy
      real,pointer:: uu(:,:),vv(:,:),uuo(:,:),vvo(:,:)
      real,pointer :: uu0(:,:),vv0(:,:),uu1(:,:),vv1(:,:)
      
      contains
      
      
      subroutine init_grid(file,ig)
        use netcdf
        character(len=120) file
        character(len=20) varstr
        character(len=30) timebasestr
        integer ierr,ig,ilon,ilat,itime
        real*8 time1 
        integer,pointer:: nx,ny
        nx=>dg(ig)%nx
        ny=>dg(ig)%ny
        
        ierr=nf90_open(path=trim(file),mode=nf90_nowrite,ncid=dg(ig)%ncfile)
        if (ierr.ne.nf90_noerr) write(*,'(4a)') 'Cannot open ',trim(file),' > ',nf90_strerror(ierr)
        
        ierr=nf90_inquire_dimension(dg(ig)%ncfile,1,varstr,dg(ig)%nt)
        if (varstr(1:4).ne.'time') write(*,'(a)') 'Time dimension not found'
        ierr=nf90_inquire_dimension(dg(ig)%ncfile,2,varstr,dg(ig)%ny)
        if (varstr(1:3).ne.'lat') write(*,'(a)') 'Latitude dimension not found'
        ierr=nf90_inquire_dimension(dg(ig)%ncfile,3,varstr,dg(ig)%nx)
        if (varstr(1:3).ne.'lon') write(*,'(a)') 'Longitude dimension not found'
        
        ierr=nf90_inq_varid(dg(ig)%ncfile,'time',itime)
        if (ierr.ne.nf90_noerr) write(*,'(2a)') 'Cannot find time variable [time] > ',nf90_strerror(ierr)
        ierr=nf90_inq_varid(dg(ig)%ncfile,'lat',ilat)
        if (ierr.ne.nf90_noerr) write(*,'(2a)') 'Cannot find latitude variable [lat] > ',nf90_strerror(ierr)
        ierr=nf90_inq_varid(dg(ig)%ncfile,'lon',ilon)
        if (ierr.ne.nf90_noerr) write(*,'(2a)') 'Cannot find longitude variable [lon]> ',nf90_strerror(ierr)
        
        if (itime.eq.0.or.ilon.eq.0.or.ilat.eq.0) then
          stop 'Netcdf coordinate error'
        endif
        
        if (ig.eq.1) then
          ierr=nf90_inq_varid(dg(ig)%ncfile,'air_u',dg(ig)%iu)
          if (ierr.ne.nf90_noerr) then
            ierr=nf90_inq_varid(dg(ig)%ncfile,'ugrd10m',dg(ig)%iu)
            if (ierr.ne.nf90_noerr) write(*,'(2a)') 'Cannot find east wind variable [air_u | ugrd10m] > ',nf90_strerror(ierr)
          endif
          ierr=nf90_inq_varid(dg(ig)%ncfile,'air_v',dg(ig)%iv)
          if (ierr.ne.nf90_noerr) then
            ierr=nf90_inq_varid(dg(ig)%ncfile,'vgrd10m',dg(ig)%iv)
            if (ierr.ne.nf90_noerr) write(*,'(2a)') 'Cannot find north wind variable [air_v | vgrd10m] > ',nf90_strerror(ierr)
          endif
        else if (ig.eq.2) then
          ierr=nf90_inq_varid(dg(ig)%ncfile,'water_u',dg(ig)%iu)
          if (ierr.ne.nf90_noerr) then
            ierr=nf90_inq_varid(dg(ig)%ncfile,'us',dg(ig)%iu)
            if (ierr.ne.nf90_noerr) then
              ierr=nf90_inq_varid(dg(ig)%ncfile,'um',dg(ig)%iu)
              if (ierr.ne.nf90_noerr) write(*,'(2a)') 'Cannot find east wind variable [water_u | us | um] > ',nf90_strerror(ierr)
            endif
          endif
          ierr=nf90_inq_varid(dg(ig)%ncfile,'water_v',dg(ig)%iv)
          if (ierr.ne.nf90_noerr) then
            ierr=nf90_inq_varid(dg(ig)%ncfile,'vs',dg(ig)%iv)
            if (ierr.ne.nf90_noerr) then
              ierr=nf90_inq_varid(dg(ig)%ncfile,'vm',dg(ig)%iv)
              if (ierr.ne.nf90_noerr) write(*,'(2a)') 'Cannot find north wind variable [water_v | vs | vm] > ',nf90_strerror(ierr)
            endif
          endif
        else if (ig.eq.3) then 
          ierr=nf90_inq_varid(dg(ig)%ncfile,'hs',dg(ig)%iu)
          if (ierr.ne.nf90_noerr) write(*,'(2a)') 'Cannot find wave height [hs] > ',nf90_strerror(ierr)
          ierr=nf90_inq_varid(dg(ig)%ncfile,'tp',dg(ig)%iv)
          if (ierr.ne.nf90_noerr) write(*,'(2a)') 'Cannot find wave period [tp] > ',nf90_strerror(ierr)
        elseif (ig.eq.4) then
          ierr=nf90_inq_varid(dg(ig)%ncfile,'sst',dg(ig)%iu)
          if (ierr.ne.nf90_noerr) write(*,'(2a)') 'Cannot find SST [sst] > ',nf90_strerror(ierr)
          dg(ig)%iv=-1
        elseif (ig.eq.5) then
          ierr=nf90_inq_varid(dg(ig)%ncfile,'dep',dg(ig)%iu)
          if (ierr.ne.nf90_noerr) write(*,'(2a)') 'Cannot find depth [dep] > ',nf90_strerror(ierr)
          dg(ig)%iv=-1
        endif

        if (dg(ig)%iu.eq.0.or.dg(ig)%iv.eq.0) then
          stop 'Netcdf variable error'
        endif
        
        timebase0=timebase
        ierr=nf90_get_att(dg(ig)%ncfile,itime,'units',timebasestr)
        if (ierr.ne.nf90_noerr) write(*,'(2a)') 'Time units not specified > ',nf90_strerror(ierr)
        if (timebasestr(1:10).ne.'days since'.or.len(trim(timebasestr)).ne.30) then
          stop 'Time units error - must be specified as "days since"'
        else
          timebase=udtime(timebasestr(12:30))
          if (timebase.lt.udtime('1582-10-15 00:00:00')) timebase=timebase-2
        endif
        if (timebase0.gt.0.and.timebase.ne.timebase0) then
          print*, timebase,timebase0
          stop 'All input files must share the same time units'
        endif
       
        allocate(dg(ig)%xx(nx),dg(ig)%yy(ny))
        allocate(dg(ig)%uu(nx,ny),dg(ig)%vv(nx,ny),dg(ig)%uuo(nx,ny),dg(ig)%vvo(nx,ny))
        allocate(dg(ig)%uu0(nx,ny),dg(ig)%vv0(nx,ny),dg(ig)%uu1(nx,ny),dg(ig)%vv1(nx,ny))
        allocate(dg(ig)%ttime(dg(ig)%nt))
        
        ierr=nf90_get_var(dg(ig)%ncfile,ilon,dg(ig)%xx)
        ierr=nf90_get_var(dg(ig)%ncfile,ilat,dg(ig)%yy)
        ierr=nf90_get_var(dg(ig)%ncfile,itime,dg(ig)%ttime)
        
        dg(ig)%dx=dg(ig)%xx(2)-dg(ig)%xx(1)
        dg(ig)%dy=dg(ig)%yy(2)-dg(ig)%yy(1)
        dg(ig)%x0=dg(ig)%xx(1)
        dg(ig)%y0=dg(ig)%yy(1)
        dg(ig)%idx=1./dg(ig)%dx
        dg(ig)%idy=1./dg(ig)%dy
                      
        dg(ig)%timeind=0
        dg(ig)%time1=dg(ig)%ttime(1)
        write(*,'(A,F15.6)') 'First time in file: ',dg(ig)%time1
        
        time1=dg(ig)%time1
        call read_slab(time1,ig)
        dg(ig)%uu0=dg(ig)%uu1
        dg(ig)%vv0=dg(ig)%vv1
        dg(ig)%time0=dg(ig)%time1

      end subroutine
      
      subroutine read_slab(ctime,igrid)
      
      use netcdf
      
      real*8 ctime,dtime
      integer ierr,igrid
      real*8,pointer:: time0,time1
      
      call setgrid(igrid)
      time0=>dg(igrid)%time0
      time1=>dg(igrid)%time1
      
      do while (time1.le.ctime.and.dg(igrid)%timeind.lt.dg(igrid)%nt)
        dg(igrid)%timeind=dg(igrid)%timeind+1
        dg(igrid)%time0=dg(igrid)%time1
        dg(igrid)%time1=dg(igrid)%ttime(dg(igrid)%timeind)
        uu0=uu1
        vv0=vv1
        ierr=nf90_get_var(dg(igrid)%ncfile,dg(igrid)%iu,uu1,(/1,1,dg(igrid)%timeind/),(/nx0,ny0,1/))
        ierr=nf90_get_var(dg(igrid)%ncfile,dg(igrid)%iv,vv1,(/1,1,dg(igrid)%timeind/),(/nx0,ny0,1/))
        if (ierr.ne.nf90_noerr) write(*,'(2a)') 'Variable read error ',nf90_strerror(ierr)
      enddo
      
      uuo=uu
      vvo=vv
      if (ctime.le.time0) then
        uu=uu0
        vv=vv0
      elseif (ctime.ge.time1) then
        uu=uu1
        vv=vv1
      else
        dtime=(ctime-time0)/(time1-time0)
        uu=(1.0-dtime)*uu0+dtime*uu1
        vv=(1.0-dtime)*vv0+dtime*vv1
      endif
      
      end subroutine
      
      
      subroutine setgrid(igrid)
      integer igrid
      nx0=dg(igrid)%nx
      ny0=dg(igrid)%ny
      x0=dg(igrid)%x0
      y0=dg(igrid)%y0
      idx=dg(igrid)%idx
      idy=dg(igrid)%idy
      uu=>dg(igrid)%uu
      vv=>dg(igrid)%vv
      uuo=>dg(igrid)%uuo
      vvo=>dg(igrid)%vvo
      uu0=>dg(igrid)%uu0
      vv0=>dg(igrid)%vv0
      uu1=>dg(igrid)%uu1
      vv1=>dg(igrid)%vv1
      end subroutine
      
      end module griddata
      
      
      
      function jday(j,m,d)
	integer j,m,d,c,jday

	if (m.lt.3) then
		m=m + 12
		j = j - 1
	endif
	c=2-floor(j/100.) + floor(j/400.);
	jday= floor(1461*(j + 4716)/4.)+floor(153*(m + 1)/5.)+d+c-1525; 
	return
	end function
      
	function udtime(strtime)

	real*8 udtime
	character*19 strtime
	integer year,month,day,hour,min,sec

	read(strtime(1:4),*) year
	read(strtime(6:7),*) month
	read(strtime(9:10),*) day
	read(strtime(12:13),*) hour
	read(strtime(15:16),*) min
	read(strtime(18:19),*) sec
	 
	udtime=dble(jday(year,month,day))+hour/24.0+min/1440.0+sec/86400.0
	return
	end function
      
      
