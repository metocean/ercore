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