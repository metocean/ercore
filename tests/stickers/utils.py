def make_dep (x,y,dep,filedep='dep.nc'):
    import netCDF4
    dep = dep.astype('float32')
    nc = netCDF4.Dataset(filedep,'w')
    nc.createDimension('lon', len(x))
    nc.createDimension('lat', len(y))
    nc.createVariable('lon','float32',('lon'))
    nc.createVariable('lat','float32',('lat'))
    nc.createVariable('dep','float32',('lat','lon'))
    nc.variables['lon'][:] = x
    nc.variables['lat'][:] = y
    nc.variables['dep'][:] = dep
    nc.close()


def make_curr2d (x,y,u,v,time,units = 'days since 1970-01-01 00:00:00:', filecurr='cur.nc'):
    import netCDF4
    u = u.astype('float32')
    v = v.astype('float32')
    nc = netCDF4.Dataset(filecurr,'w')
    nc.createDimension('lon', len(x))
    nc.createDimension('lat', len(y))
    nc.createDimension('time', len(time))
    nc.createVariable('lon','float32',('lon'))
    nc.createVariable('lat','float32',('lat'))
    t = nc.createVariable('time','float32',('time'))
    t.units = units
    nc.createVariable('u','float32',('time', 'lat','lon'))
    nc.createVariable('v','float32',('time', 'lat','lon'))
    nc.variables['lon'][:] = x
    nc.variables['lat'][:] = y
    nc.variables['time'][:] = time
    nc.variables['u'][:] = u
    nc.variables['v'][:] = v
    nc.close()