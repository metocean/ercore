      module release
      
      integer nr,npmax
      
      type relstr
        integer :: np,sticky,npsplit,npr
        real :: lon,lat,rstart,rend,dw_min,dw_max,cw_min,cw_max,refloat,halflife
        real*8, pointer, dimension(:) :: px,py,mass
        real, pointer, dimension(:) :: srand,age
        integer, pointer, dimension(:) :: psc
      end type relstr
      
      type(relstr), target, allocatable :: rel(:)
      
      real*8,pointer,dimension(:) :: px,py,mfx
      integer,pointer,dimension(:) :: psc
      
      contains
      
      
      subroutine init_release(fn)
      integer fn
      integer, dimension(10) :: np,sticky
      real, dimension(10) :: lon,lat,rstart,rend,dw_min,dw_max,cw_min,cw_max,refloat,halflife
          
      namelist/releases/nr,np,lon,lat,rstart,rend,dw_min,dw_max,cw_min,cw_max,sticky,refloat,halflife
      read(fn,nml=releases)
      
      allocate(rel(nr))
      npmax=maxval(np)
           
      do ir=1,nr
        rel(ir)%np=np(ir)
        rel(ir)%sticky=sticky(ir)
        rel(ir)%lon=dble(lon(ir))
        rel(ir)%lat=dble(lat(ir))
        rel(ir)%rstart=rstart(ir)
        rel(ir)%rend=rend(ir)
        rel(ir)%dw_min=dw_min(ir)
        rel(ir)%dw_max=dw_max(ir)
        rel(ir)%cw_min=cw_min(ir)
        rel(ir)%cw_max=cw_max(ir)
        rel(ir)%refloat=refloat(ir)
        rel(ir)%halflife=halflife(ir)
        allocate(rel(ir)%px(np(ir)),rel(ir)%py(np(ir)),rel(ir)%psc(np(ir)))
        allocate(rel(ir)%mass(np(ir)),rel(ir)%age(np(ir)))
        rel(ir)%px=dble(lon(ir))
        rel(ir)%py=dble(lat(ir))
        rel(ir)%psc=0
        rel(ir)%mass=1.0
        rel(ir)%age=0.
        rel(ir)%npr=0
      enddo
      
      end subroutine
      
      
      
      subroutine dump_release(ir)
        integer ir
    
        write(*,'(/,a,i0)') 'Release ',ir
        write(*,'(a,i0,a)') 'with ',rel(ir)%np,' particles'
        write(*,'(a,f12.4,a,f12.4)') 'Release location: ',rel(ir)%lon,',',rel(ir)%lat
        write(*,'(a,f20.4,a,f20.4)') 'Time: ',rel(ir)%rstart,' s after start until ',rel(ir)%rend,' s after start'
        write(*,'(a,f12.4,a,f12.4)') 'Downwind factor min: ',rel(ir)%dw_min,' max: ',rel(ir)%dw_max
        write(*,'(a,f12.4,a,f12.4)') 'Crosswind factor min: ',rel(ir)%cw_min,' max: ',rel(ir)%cw_max
        write(*,'(a,f12.4,a)') 'Refloat time: ',rel(ir)%refloat,' s '
        write(*,'(a,f12.4,a)') 'Halflife: ',rel(ir)%halflife,' s '
      
      end subroutine
      
      
      subroutine add_release(it,dt)
        integer it,nrel,npr
        real dt,t1,t2
        
        do ir=1,nr
          if (rel(ir)%npr.eq.rel(ir)%np) cycle
          t1=(it*dt-rel(ir)%rstart)
          t2=(rel(ir)%rend-rel(ir)%rstart)
          if (t1.ge.0.) then
            npr=rel(ir)%npr
            if (t2.le.0) then
              nrel=rel(ir)%np
            else
              nrel=floor(min(t1/t2,1.)*rel(ir)%np)-npr
            endif
            rel(ir)%npr=npr+nrel
          endif
        enddo
        
      end subroutine
      

      
      subroutine poutput(tt,dt)
      real*8 tt,dt
      
        do ir=1,nr
          write(20,'(I8,X,F12.5,X,I8)')ir,tt,rel(ir)%npr
          do ip=1,rel(ir)%npr
            write(20,'(3F12.5,I8,F12.5)')rel(ir)%px(ip),rel(ir)%py(ip),rel(ir)%mass(ip),dt*rel(ir)%psc(ip),rel(ir)%age(ip)/86400.
          enddo
        enddo
      
      end subroutine
      
      
      
      end module release

