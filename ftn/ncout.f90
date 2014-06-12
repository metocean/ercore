    module ncout
    use netcdf
    
    integer ncid,ierr,itind
    integer dtimeid,dlonid,dlatid,densid
    integer vtimeid,vlonid,vlatid,vensid,vconcid
    
    contains
      subroutine close_ncout()
        ierr=nf90_close(ncid)
      end subroutine
      
      subroutine init_ncout(filename,nr)
        character*120 filename
        integer nr,nrel(nr),maskid
        
        if (filename.eq.'')then
          ncid=0
          return
        endif
        itind=0
        
        do i=1,nr
          nrel(i)=i
        enddo
        
        ierr=nf90_create(filename,nf90_clobber,ncid)
        ierr=nf90_def_dim(ncid,'time',nf90_unlimited,dtimeid)
        ierr=nf90_def_dim(ncid,'lon',ognx,dlonid)
        ierr=nf90_def_dim(ncid,'lat',ogny,dlatid)
        ierr=nf90_def_dim(ncid,'ens',nr,densid)
        
        call def_ncvar('time',1,(/dtimeid/),vtimeid,'time','days since 1-1-1',nf90_double)
        call def_ncvar('lon',1,(/dlonid/),vlonid,'longitude','degrees_east',nf90_double)
        call def_ncvar('lat',1,(/dlatid/),vlatid,'latitude','degrees_north',nf90_double)
        call def_ncvar('ens',1,(/densid/),vensid,'release','release_number',nf90_int)
        call def_ncvar('mask',2,(/dlonid,dlatid/),maskid,'wetmask','fraction',nf90_int)
        call def_ncvar('dens',4,(/dlonid,dlatid,dtimeid,densid/),vconcid,'particle density','',nf90_float)
        
        ierr=nf90_enddef(ncid)
        
        ierr=nf90_put_var(ncid,vlonid,xlon)
        ierr=nf90_put_var(ncid,vlatid,ylat)
        ierr=nf90_put_var(ncid,vensid,nrel)
        ierr=nf90_put_var(ncid,maskid,mask)
      end subroutine
      
      
      subroutine def_ncvar(name,nvdims,vdims,varid,long_name,units,nf_type)
        integer vdims(nvdims),ierr
        integer nvdims,varid,nf_type
        character*(*) name,long_name,units
        ierr=nf90_def_var(ncid,name,nf_type,vdims,varid)
        ierr=nf90_put_att(ncid,varid,'long_name',long_name)
        ierr=nf90_put_att(ncid,varid,'units',units)
        return
      end subroutine
      
      subroutine output_dens(ttime,ir)
        if (ir.eq.1) itind=itind+1
          ierr=nf90_put_var(ncid,vtimeid,ttime,(/itind/))
          ierr=nf90_put_var(ncid,vconcid,cconc,(/1,1,itind,ir/),(/ognx,ogny,1,1/))
        endif
      end subroutine
      
    end module ncout