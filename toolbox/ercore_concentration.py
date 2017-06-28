#!/usr/bin/env python
# import requried libraries
# from read_result import read_ercore
import sys
import numpy
import pyproj
from ercore.shoreline import Shoreline,Boundary
from _flib_ercore import shoreline
import pdb
from scipy.io import loadmat
from scipy.spatial import ConvexHull
from matplotlib.path import Path
import netCDF4
import datetime


def read_ercore_outputs(fname):
    """ >> allow to use wildcard in fname to combine several files e.g. monthly _1,_2 etc...
    # NEW FORMAT : Time    id  x   y   z   state   age mass    zbottom elev
    # OLD FORMAT : Time id  x   y   z   state   age mass  
    """  
    data=numpy.loadtxt(fname,delimiter = '\t',dtype=None,skiprows = 1)
    data_susp = convert_ercore_coords(data[data[:,5] == 1,:],epsg_code_cartesian)
    data_depos = convert_ercore_coords(data[data[:,5] == -1,:]epsg_code_cartesian)
    time = numpy.unique(data[:,0])
    return time,data_susp,data_depos

def convert_ercore_coords(data):
    """data is data matrix from ercore [Mx 8 or 10 columns]
    # epsg_code = projection used to convert colums 3 and 4 [x,y]
    # wgs84 epsg_code=2193 , nztm epsg_code=2193
    # 
    # ercore outputs are in geograpgical coordinates wgs84
    # convert to cartesian coordinates using epsg code
    # default is nzdg2000=nztm / see spatialreference.org for others
    """
    p_cart = pyproj.Proj(init=epsg_cart_str)  # nzgd 2000 = NZTM
    x,y = p_cart(data[:,2],data[:,3]) # convert to cartesian coordinates
    data[:,2] = x[:]
    data[:,3] = y[:]
    return data

def convert_xy_coords(x,y,epsg_code_from = 'epsg:2193',epsg_code_to = 'epsg:4326'):
    """data is data matrix from ercore [Mx 8 or 10 columns]
    # epsg_code = projection used to convert colums 3 and 4 [x,y]
    # wgs84 epsg_code=2193 , nztm epsg_code=2193
    # 
    # ercore outputs are in geograpgical coordinates wgs84
    # convert to cartesian coordinates using epsg code
    # default is nzdg2000=nztm / see spatialreference.org for others
    """
    # epsg_str_from=('espg:%i' % (epsg_code_from))
    # epsg_str_to=('espg:%i' % (epsg_code_to))

    proj_from = pyproj.Proj(init = epsg_code_from)  # nzgd 2000 = NZTM
    proj_to = pyproj.Proj(init = epsg_code_to)  # nzgd 2000 = NZTM
    xconv,yconv = pyproj.transform(proj_from, proj_to, x, y)
    return xconv,yconv

def read_shoreline_ercore(bnd_file,x,y):
    """read GNOME shoreline as used in ercore, and flag points on land"""
    shoreline.read_shoreline(bnd_file)
    slx = shoreline.slx
    sly = shoreline.sly
    polyi = shoreline.polyi # id of shoreline segment
    polyn = shoreline.polyn # nb node of each segement
    points = numpy.ones([x.shape[0],2]) # allocation for correct format in Path
    on_land = numpy.zeros(x.shape[0],dtype=bool)
    points[:,0] = x[:]
    points[:,1] = y[:]
    # check if points are within shoreline polys or not
    # loop through all the different shoreline polygon
    for i in range(0,len(polyi)):
        polygon = numpy.ones([polyn[i],2])
        polygon[:,0] = slx[range(polyi[i]-1,polyn[i]+polyi[i]-1)]
        polygon[:,1] = sly[range(polyi[i]-1,polyn[i]+polyi[i]-1)]
        path = Path(polygon)  # matplotlib.path.Path
        on_land_tmp = path.contains_points(points)# path = matplotlib.path.Path(polygon)
        on_land = on_land | on_land_tmp # 
    # import pdb;pdb.set_trace()
    # numpy.savetxt('test',points[on_land])
    return on_land,shoreline

def conc_calc(pmatrix,x_recep,y_recep,nh = 1./20.,shoreline_file = '',epsg_code = 2193):
""" compute particle concentration
# 
# pmatrix = matrix with particle data [time id x y z etc..]
# x_recep = x-coordinates of receptor grid - array
# y_recep  = y-coordinates of receptor grid - array
# nh=1./20. = fraction used to define the number of particles to include in conc computation, will use nh-th closest
# shoreline='' = shoreline file, optional
"""
    
    #load shoreline, if applicable
    if shoreline_file :
        epsg_str=('espg:%i' % (epsg_code))
        shoreline.read_shoreline(shoreline_file) 
        x_shore,y_shore = convert_xy_coords(shoreline.slx,shoreline.sly,epsg_code_from='epsg:4326',epsg_code_to=epsg_cart_str)
        # unit matrix for shore proximity computation
        xunit = numpy.linspace(-1., 1., (1.0/0.05)+1)
        yunit = numpy.linspace(-1., 1., (1.0/0.05)+1)
        xxu, yyu = numpy.meshgrid(xunit, yunit)

    # load particles position
    xp = pmatrix[:,2]
    yp = pmatrix[:,3]
    np = len(xp)
    load = numpy.ones(np)
    # allocate concentration,nb_part arrays
    conc = numpy.zeros(x_recep.shape)
    nb_part = numpy.zeros(x_recep.shape)

    #if np == 0 or less than 3 , return null concentration (not enough to computea convex hull)
    if np <3: return conc,nb_part,np
    # define convex hull of particle cloud 
    # and flag receptors outside of it to skip them
    conv_hull = ConvexHull(numpy.vstack([xp,yp]).T)

    hull = Path(numpy.vstack([xp[conv_hull.vertices],yp[conv_hull.vertices]]).T)  # matplotlib.path.Path
    in_hull =hull.contains_points(numpy.vstack([xr,yr]).T)# path = matplotlib.path.Path(polygon)

    # loop through receptors
    for r in range(0,len(x_recep)):

        # skip if receptor is on land
        if on_land[r]: conc[r]=0.; continue
       # skip if receptor is outside of particle cloud convex hull
        if not in_hull[r]: conc[r]=0.; continue
    
        #compute distance of every particle from receptor
        dx = numpy.abs(xp-xr[r])
        dy =  numpy.abs(yp-yr[r]);
        dist = numpy.sqrt(dx**2+dy**2)
        #sort distance
        dist_sort = numpy.sort(dist)
        # take the nh-th closest particles - the subset will be used for concentration computation
        # **numpy.argsort  returns id that would sort `dist`
        part2use = numpy.argsort(dist)[0:int(numpy.round(np*nh))]
    
        # define kernel bandwidths as per RL3 method Vitali et al. 2006:
        # for each direction x and y, the bandwidth is defined as the
        # minimum value between i) the maximum projected distance of the
        # particles within the neighbourhood and ii) twice the standard
        # deviation of the projected distances within the neighbourhood.    
      
        lx = numpy.max(dx[part2use])
        ly = numpy.max(dy[part2use])
        lx = numpy.minimum(lx,2*numpy.std(dx[part2use]))
        ly = numpy.minimum(ly,2*numpy.std(dy[part2use]))
        # The aspect ratio (e.g. lx/ly) of the bandwidths are limited to be no
        # greater than 5:1 to prevent unrealistically elongated kernels, with
        # the smaller value increased.
        lx = numpy.maximum(lx,0.2*ly)
        ly = numpy.maximum(ly,0.2*lx)
      
        # find particles within kernel bandwidths.
        in_kernel = (dx<=lx) & (dy<=ly)

        if not in_kernel.any() or in_kernel.sum()<5:
            # if 0 or less than 5 particles in neighborhood - skip to next receptor
            conc[r]=0.
            continue
        else:
            #kernel function
            Kx = 0.75*(1-(dx[in_kernel]/lx)**2);
            Ky = 0.75*(1-(dy[in_kernel]/ly)**2);
            #concentration
            conc[r] = numpy.sum( (load[in_kernel]*Kx*Ky)/(lx*ly) )   
            nb_part[r] = in_kernel.sum()
            
            # correct for land proximity - to test
            if shoreline_file :
                Ldx = numpy.abs( x_shore -xr[r] )
                Ldy = numpy.abs( y_shore -yr[r] )
                inear = (Ldx<=lx) & (Ldy <=ly)
                
                if inear.any() & (conc[r]>0.) : 
                    #grid bandwidth
                    x1 = xr[r]+lx*numpy.hstack(xxu)
                    y1 = yr[r]+ly*numpy.hstack(yyu)
                    # check if points are within shoreline poly
                    # convert grids to lon/lat before doing the shoreline check
                    xx1,yy1 = convert_xy_coords(x1,y1,epsg_code_from=epsg_cart_str,epsg_code_to='epsg:4326')
                    inland,shore = read_shoreline_ercore(shoreline_file,xx1,yy1) # flag points on land
                    Ldx = numpy.abs(x1[inland]-xr[r])
                    Ldy = numpy.abs(y1[inland]-yr[r])
                    LKx = 0.75*(1-(Ldx/lx)**2)
                    LKy = 0.75*(1-(Ldy/ly)**2)
                    A = (1-0.0025*numpy.sum(numpy.sum(LKx*LKy)))
                    # apply correctionto concentration
                    conc[r]=conc[r]/A
    
    return conc,nb_part,np

def netcdf_file_unstruct_init(fname,lon,lat,tri,lev=None):
    """netcdf-writing subfunction for unstructured grid"""
    # Dimensions
    dataset = netCDF4.Dataset(fname, 'w',  format ='NETCDF4_CLASSIC') 
    dim_node = dataset.createDimension('node', len(lon))
    # dim_lat = dataset.createDimension('lat', len(lon))
    # dim_lon = dataset.createDimension('lon', len(lat))
    dim_nb_vertex = dataset.createDimension('nb_vertex', 3) # 3 = triangle
    dim_nb_elem = dataset.createDimension('nb_elem',tri.shape[0]) 
    dim_time = dataset.createDimension('time', None) # use unlimited time dimension

    if lev is not None:    
        level = dataset.createDimension('level', len(lev))
        var_lev = dataset.createVariable('lev', numpy.int32, ('level',)) 
    # Variables
    # Dataset.createVariable(<var_id>, <type>, <dimensions>)
    var_time = dataset.createVariable('time', numpy.float64, ('time',))
    var_lat = dataset.createVariable('lat', numpy.float64,   ('node',))
    var_lon = dataset.createVariable('lon', numpy.float64, ('node',))
    var_x = dataset.createVariable('easting', numpy.float64, ('node',))
    var_y = dataset.createVariable('northing', numpy.float64, ('node',))
    var_dep = dataset.createVariable('dep', numpy.float64, ('node',))
    var_elem= dataset.createVariable('elem', numpy.float64, ('nb_elem','nb_vertex'))

    if lev is not None :
    #var_conc = dataset.createVariable('np_total_depos', numpy.float64,('time','node'))            
    #var_conc = dataset.createVariable('np_total_depos', numpy.float64,('time','level','node'))   
    #var_conc = dataset.createVariable('np_total_depos', numpy.float64,('time','level','lat','lon'))    
        var_conc_susp = dataset.createVariable('conc_susp', numpy.float64,('time','level','node'))
        var_conc_depos = dataset.createVariable('conc_depos', numpy.float64,('time','level','node'))
        var_np_used_susp = dataset.createVariable('np_used_susp', numpy.float64,('time','level','node'))
        var_np_used_depos = dataset.createVariable('np_used_depos', numpy.float64,('time','level','node'))
        var_np_total_susp = dataset.createVariable('np_total_susp', numpy.float64,('time','level','node'))
        var_np_total_depos  = dataset.createVariable('np_total_depos', numpy.float64,('time','level','node'))
    else:
        var_conc_susp = dataset.createVariable('conc_susp', numpy.float64,('time','node'))
        var_conc_depos = dataset.createVariable('conc_depos', numpy.float64,('time','node'))
        var_np_used_susp = dataset.createVariable('np_used_susp', numpy.float64,('time','node'))
        var_np_used_depos = dataset.createVariable('np_used_depos', numpy.float64,('time','node'))
        var_np_total_susp = dataset.createVariable('np_total_susp', numpy.float64,('time','node'))
        var_np_total_depos = dataset.createVariable('np_total_depos', numpy.float64,('time','node'))  
    
    # Fill content of basic variables
    x,y=convert_xy_coords(lon,lat,epsg_code_from='epsg:4326',epsg_code_to=epsg_cart_str) # from WGS84 to epsg code input
    var_lon[:] = lon
    var_lat[:] = lat
    var_x[:] = x
    var_y[:] = y
    var_elem[:] = tri

    # Global Attributes 
    dataset.description = 'particle concentration processing'  
    dataset.history = '' 
    dataset.source = '' 

    # Variable Attributes  
    var_time.units = 'days since 0001-01-01 00:00:00'
    var_time.standard_name = 'time'
    var_lat.units = 'degree_north'
    var_lat.standard_name = 'latitude'
    var_lon.units = 'degree_east'
    var_lon.standard_name = 'longitude'  
    var_x.units = 'meters'
    var_x.standard_name = 'eastings'  
    var_y.units = 'meters'
    var_y.standard_name = 'northings' 
    var_dep.units = 'meters'
    var_dep.standard_name = 'water_depth' 
    var_dep.convention= 'negative_down'
    var_elem.units = ''
    var_elem.standard_name = 'triangular_element_matrix' 

    var_conc_susp.units = 'load/m3'
    var_conc_susp.standard_name = 'suspended_particle_concentration' 
    var_conc_susp.description= 'each particle is assumed to have load=1' 

    var_conc_depos.units = 'load/m3'
    var_conc_depos.standard_name = 'deposited_particle_concentration' 
    var_conc_depos.description= 'each particle is assumed to have load=1' 

    var_np_used_susp.units = ''
    var_np_used_susp.standard_name = 'number_particle_used_to_compute_conc_susp'
    var_np_used_susp.description= 'number of particles used to compute conc_susp at each receptors'  

    var_np_used_depos.units = ''
    var_np_used_depos.standard_name = 'number_particle_used_to_compute_conc_depos'
    var_np_used_depos.description= 'number of particles used to compute conc_depos at each receptor, each time step'  

    var_np_total_susp.units = ''
    var_np_total_susp.standard_name = 'total_suspended_particles'
    var_np_total_susp.description= 'total number of suspended particles'

    var_np_total_depos.units = ''
    var_np_total_depos.standard_name = 'total_deposited_particles'
    var_np_total_depos.description= 'total number of deposited particles'

    # w_nc_var.setncatts({'long_name': u"mean Daily Air temperature departure",\
    #                 'units': u"degK", 'level_desc': u'Surface',\
    #                 'var_desc': u"Air temperature departure",\
    #                 'statistic': u'Mean\nM'})

    return dataset 

def netcdf_file_regular_init(fname,lon,lat,lev=None):
    """netcdf-writing subfunction for structured/regular grid"""
    pass


# MAIN---------------------------------------------------------------------------------------------------------

if __name__  ==  '__main__':

    """The function computes the suspended and deposited concentration fields 
    for a given ERcore outputs file with particles position, for a given grid of receptors
    usage : ercore_concentration.py input_ercore_filename output_netcdf_name receptor_grid_filename epsg_code
    
    input_ercore_filename : name of ERcore ascii output file 
    output_netcdf_name : name of output netcdf file with concentration fields
    receptor_grid_filename : name of file containing the receptor grid info - <.mat> file for now
    epsg_code_cartesian : epsg code for conversion of lon/lat particles position to cartesian coordinates
    see  http://www.epsg-registry.org/ for code
    NZGD2000 / New Zealand Transverse Mercator 2000 is 2193
    WGS84 is 4326
    
    """
    # INPUTS
    fname_in = '/home_old/simon/0201_lyt_ercore_out/site5/2003/ercore.dispo_site5_surface4m_class1_2003_1.out'   
    fname_out = '/home/simon/ercore_tests/concentration_computation/dispo_site5_surface4m_class1_2003_1.nc'
    receptor_grid = 'home/simon/ercore_tests/concentration_computation/grid3_net.mat'
    # fname_in = int(sys.argv[0])
    # fname_in = int(sys.argv[1])
    # receptor_grid = int(sys.argv[2])

    epsg_code_cartesian=2193
    epsg_cart_str=('espg:%i' % (epsg_code_cartesian))
    global  epsg_code_cartesian
    global  epsg_cart_str

    # load ERCORE outputs
    time,susp,depos = read_ercore_outputs(fname_in)
    
    # load Receptor grids
    # - can be same as a netcdf grid
    # - or extents + resolution
    # - or input x/y
    # - matlab .mat file with receptors and triangle matrix ([elem x 3])

    #receptors=loadmat('lyt_channel_new_net.mat')
    receptors = loadmat(receptor_grid)

    # reformat information
    lonr = receptors['recep']['lon'][0][0][:].squeeze()
    latr = receptors['recep']['lat'][0][0][:].squeeze()
    xr = receptors['recep']['x'][0][0][:].squeeze()
    yr = receptors['recep']['y'][0][0][:].squeeze()
    zr = receptors['recep']['z'][0][0][:].squeeze()
    tri = receptors['tri'][:].squeeze()
    # the matlab file must have fields recep tri 
    # recep has several sub-fields recep.lon recep.lat recep.x recep.y recep.z
    # access subfields like this : receptors['recep']['lon']

    # load shoreline file, if applicable and flag points on land
    bnd_file = 'southnz_shoreline.bnd' # this will be in lat/lon
    on_land,shore = read_shoreline_ercore(bnd_file,lonr,latr)

    nc_dataset = netcdf_file_unstruct_init(fname_out,lonr,latr,tri,lev=None)
    # add some info to netcdf
    nc_dataset.variables['easting'][:] = xr
    nc_dataset.variables['northing'][:] = yr
    nc_dataset.variables['dep'][:] = zr
    
    for t_cnt,t in enumerate(time[0:20]):
        id_susp=susp[:,0] == t#numpy.where(susp[:,0] == t)
        id_depos=depos[:,0] == t#numpy.where(depos[:,0] == t)  

        # compute concentration, no land proximity correction for now
        conc_susp,npart_susp,np_susp_total = conc_calc(susp[id_susp,:],xr,yr,nh=1./20.,shoreline_file=''bnd_file'')
        # integrate concetrtaion of suspedned particle over water depth        
        conc_susp=conc_susp/zr

        conc_dep,npart_depos,np_depos_total = conc_calc(depos[id_depos,:],xr,yr,nh=1./20.)

        tt=netCDF4.num2date(t,units = nc_dataset.variables['time'].units)
        print tt.isoformat(' ')
        #import pdb;pdb.set_trace()

        # write to netcdf file
        nc_dataset.variables['conc_susp'][t_cnt,:] = conc_susp
        nc_dataset.variables['conc_depos'][t_cnt,:] = conc_dep

        nc_dataset.variables['np_used_susp'][t_cnt,:] = npart_susp
        nc_dataset.variables['np_used_depos'][t_cnt,:] = npart_depos

        nc_dataset.variables['np_total_susp'][t_cnt,:] = np_susp_total
        nc_dataset.variables['np_total_depos'][t_cnt,:]= np_depos_total

        nc_dataset.variables['time'][t_cnt] = t # the time saved in ercore output files is alread CF-compliant days since 1-1-1
    
    #close netcdf file
    nc_dataset.close()

        # TO DO 
        # >>>>> include cases for when there are no particle / no
        # >>>>> write results to a CF- netcdf file for further processing
        # >>>>> save concentration susp and depos, number of particles for each etc.., depth, 
        #  >>> include option for concetration computation at different level in water colums
