#!/usr/bin/env python
# import requried libraries
# from read_result import read_ercore
import sys
import os
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
    """ TO DO : >> allow to use wildcard in fname to combine several files e.g. monthly _1,_2 etc...
    # NEW FORMAT : Time    id  x   y   z   state   age mass    zbottom elev
    # OLD FORMAT : Time id  x   y   z   state   age mass  
    """  
    data=numpy.loadtxt(fname,delimiter = '\t',dtype=None,skiprows = 1)
    data_susp = convert_ercore_coords(data[data[:,5] == 1,:])
    data_depos = convert_ercore_coords(data[data[:,5] == -1,:])
    time = numpy.unique(data[:,0])
    return time,data_susp,data_depos

def convert_ercore_coords(data):
    """data is data matrix from ercore [Mx 8 or 10 columns]
    # epsg_cart_str = projection used to convert colums 3 and 4 [x,y]
    # wgs84 epsg_code=2193 , nztm epsg_code=2193
    # 
    # ercore outputs are in geograpgical coordinates wgs84
    # convert to cartesian coordinates using epsg code
    # default is nzdg2000=nztm / see spatialreference.org for others
    """
    
    p_cart = pyproj.Proj(init=epsg_cart_str)  # 
    x,y = p_cart(data[:,2],data[:,3]) # convert particles lon/lat to cartesian coordinates
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
    # shoreline = shoreline file, optional
    """

    #load shoreline, if applicable
    if shoreline_file :
        shoreline.read_shoreline(shoreline_file)
        # get cartesian shoreline coordinates 
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
    in_hull =hull.contains_points(numpy.vstack([x_recep,y_recep]).T)# path = matplotlib.path.Path(polygon)

    # loop through receptors
    for r in range(0,len(x_recep)):

        # skip if receptor is on land
        if on_land[r]: conc[r]=0.; continue
       # skip if receptor is outside of particle cloud convex hull
        if not in_hull[r]: conc[r]=0.; continue
    
        #compute distance of every particle from receptor
        dx = numpy.abs(xp-x_recep[r])
        dy =  numpy.abs(yp-y_recep[r]);
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
                Ldx = numpy.abs( x_shore -x_recep[r] )
                Ldy = numpy.abs( y_shore -y_recep[r] )
                inear = (Ldx<=lx) & (Ldy <=ly)

                if inear.any() & (conc[r]>0.) : 
                    #create local grid bandwidth
                    x1 = x_recep[r]+lx*numpy.hstack(xxu)
                    y1 = y_recep[r]+ly*numpy.hstack(yyu)
                    # check if points are within shoreline poly
                    # convert grids coords to lon/lat before doing the shoreline check
                    xx1,yy1 = convert_xy_coords(x1,y1,epsg_code_from=epsg_cart_str,epsg_code_to='epsg:4326')
                    inland,shore = read_shoreline_ercore(shoreline_file,xx1,yy1) # flag points on land
                    Ldx = numpy.abs(x1[inland]-x_recep[r])
                    Ldy = numpy.abs(y1[inland]-y_recep[r])
                    LKx = 0.75*(1-(Ldx/lx)**2)
                    LKy = 0.75*(1-(Ldy/ly)**2)
                    A = (1-0.0025*numpy.sum(numpy.sum(LKx*LKy)))
                    # apply correctionto concentration
                    conc[r]=conc[r]/A
    
    return conc,nb_part,np

def netcdf_file_unstruct_init(fname,recep_grid,conc_level=['dav']):
    """netcdf-writing subfunction for unstructured grid
       
       When there are more than one conc_level to store
       variables will have different levels in dimensions

    """
    # check that conc_level is a list
    if  not isinstance(conc_level,list) : conc_level= [conc_level]

    # Dimensions
    dataset = netCDF4.Dataset(fname, 'w',  format ='NETCDF4_CLASSIC')
    # dataset = netCDF4.Dataset(fname, 'w') 
    dim_node = dataset.createDimension('node', len(recep_grid['lonr']))
    # dim_lat = dataset.createDimension('lat', len(lon))
    # dim_lon = dataset.createDimension('lon', len(lat))
    dim_nb_vertex = dataset.createDimension('nb_vertex', 3) # 3 = triangle
    dim_nb_elem = dataset.createDimension('nb_elem',recep_grid['tri'].shape[0]) 
    dim_time = dataset.createDimension('time', None) # use unlimited time dimension
    dim_lev = dataset.createDimension('level', len(conc_level)) # 
    dim_nchar = dataset.createDimension('nchar', 4)

    # Variables
    # Dataset.createVariable(<var_id>, <type>, <dimensions>)
    var_time = dataset.createVariable('time', numpy.float64, ('time',))
    var_lat = dataset.createVariable('lat', numpy.float64,   ('node',))
    var_lon = dataset.createVariable('lon', numpy.float64, ('node',))
    var_x = dataset.createVariable('easting', numpy.float64, ('node',))
    var_y = dataset.createVariable('northing', numpy.float64, ('node',))
    var_dep = dataset.createVariable('dep', numpy.float64, ('node',))
    var_elem = dataset.createVariable('elem', numpy.float64, ('nb_elem','nb_vertex'))
    # var_lev = dataset.createVariable('lev', str, ('level'))
    var_lev = dataset.createVariable('lev', 'S1', ('level','nchar'))
    # netCDF4.stringtochar()
  

    var_conc_susp = dataset.createVariable('conc_susp', numpy.float64,('time','level','node')) 
    var_np_used_susp = dataset.createVariable('np_used_susp', numpy.int32,('time','level','node'))
    var_np_total_susp = dataset.createVariable('np_total_susp', numpy.int32,('time','level','node'))

    var_conc_depos = dataset.createVariable('conc_depos', numpy.float64,('time','node'))
    var_np_used_depos = dataset.createVariable('np_used_depos', numpy.int32,('time','node'))
    var_np_total_depos  = dataset.createVariable('np_total_depos', numpy.int32,('time','node'))
    
    # Fill content of basic variables
    x,y=convert_xy_coords(recep_grid['lonr'],recep_grid['latr'],epsg_code_from='epsg:4326',epsg_code_to=epsg_cart_str) # from WGS84 to epsg code input
    var_lon[:] = recep_grid['lonr']
    var_lat[:] = recep_grid['latr']
    var_x[:] = recep_grid['xr']
    var_y[:] = recep_grid['yr']
    var_elem[:] = recep_grid['tri'] # expect tri = [NB_ELEM x 3]
    var_dep[:] = recep_grid['zr']
    for lev_cnt,lev_name in enumerate(conc_level):
        # str_out = netCDF4.stringtochar(np.array(['test'], 'S4'))
        str_out = netCDF4.stringtochar(numpy.array(conc_level[lev_cnt], 'S4'))
        # var_lev[lev_cnt,:] = conc_level[lev_cnt]
        var_lev[lev_cnt,:] = str_out

    # Global Attributes 
    dataset.Convention = 'CF'  
    dataset.title = '' 
    dataset.institution = 'MetOcean'
    dataset.source = ''
    dataset.history = ''
    dataset.references = ''
    dataset.comment = ''

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

    var_np_total_susp.units = ''
    var_np_total_susp.standard_name = 'total_suspended_particles'
    var_np_total_susp.description= 'total number of suspended particles'

    var_np_used_susp.units = ''
    var_np_used_susp.standard_name = 'number_particle_used_to_compute_conc_susp'
    var_np_used_susp.description= 'number of particles used to compute conc_susp at each receptors'

    var_conc_depos.units = 'load/m3'
    var_conc_depos.standard_name = 'deposited_particle_concentration' 
    var_conc_depos.description= 'each particle is assumed to have load=1' 

    var_np_used_depos.units = ''
    var_np_used_depos.standard_name = 'number_particle_used_to_compute_conc_depos'
    var_np_used_depos.description= 'number of particles used to compute conc_depos at each receptor, each time step'  

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

def load_recep_mat(fname):
    """function to load a matlab-based receptor grid
    the matlab file must have fields recep tri 
    recep has several sub-fields recep.lon recep.lat recep.x recep.y recep.z tri
    access subfields like this : receptors['recep']['lon']
    """

    receptors = loadmat(fname) # matlab file with receptor grid
    # reformat information
    lonr = receptors['recep']['lon'][0][0][:].squeeze()
    latr = receptors['recep']['lat'][0][0][:].squeeze()
    # xr = receptors['recep']['x'][0][0][:].squeeze()
    # yr = receptors['recep']['y'][0][0][:].squeeze()
    xr,yr = convert_xy_coords(lonr,latr,epsg_code_from = 'epsg:4326',epsg_code_to = epsg_str_cart) 
    zr = receptors['recep']['z'][0][0][:].squeeze()
    tri = receptors['tri'][:].squeeze()

    # Package to a dictionnary
    recep = dict([('lonr', lonr),('latr', latr),('xr', xr),('yr', yr),('zr', zr),('tri', tri)])

    # return lonr,latr,xr,yr,zr,tri
    return recep

def load_recep_nc(fname):
    """function to load a netcdf-based receptor grid"""

    dataset = netCDF4.Dataset(fname, 'r')
    if 'node' in dataset.dimensions.keys(): # unstructured grid
        lonr = dataset.variables['lon'][:]
        latr = dataset.variables['lat'][:]
        # convert from wgs84 to cartesian coordinates
        xr,yr = convert_xy_coords(lonr,latr,epsg_code_from = 'epsg:4326',epsg_code_to = epsg_cart_str)

        if 'nv' in dataset.variables.keys(): 
            tri = dataset.variables['nv'][:,:]
            tri = tri.T.astype(int) # format matrxi to [N_ELEM x 3]
        elif 'elem' in dataset.variables.keys():
            tri = dataset.variables['elem'][:,:]
            tri = tri.T.astype(int) # format matrxi to [N_ELEM x 3]
        # add depth variable
        if 'dep' in dataset.variables.keys():           
            zr =  dataset.variables['dep'][:]
            # if  'positive' in dataset.variables['dep'].ncattrs():
            #     if dataset.variables['dep'].positive == 'down': 
            #         zr = -zr # use postive down convention as in ERcore
    else : # regular grid
        lonr = dataset.variables['lon'][:]
        latr = dataset.variables['lat'][:]
        # convert from wgs84 to cartesian coordinates
        xr,yr = convert_xy_coords(lonr,latr,epsg_code_from = 'epsg:4326',epsg_code_to = epsg_cart_str)
        lonr,latr = numpy.meshgrid(lonr,latr)
        xr,yr = numpy.meshgrid(lonr,latr)
        tri = []
        # add depth variable
        if 'dep' in dataset.variables.keys():           
            zr = dataset.variables['dep'][:,:]
        
        #  May need to set as 1 column matrix, and reshape later ? 
        #  to decide when adding content to netcdf_file_regular_init 
    
    # Package to a dictionnary
    recep = dict([('lonr', lonr),('latr', latr),('xr', xr),('yr', yr),('zr', zr),('tri', tri)])

    # return lonr,latr,xr,yr,zr,tri
    return recep

def load_recep_d3d_nc(fname):
    """function to load a netcdf-based receptor grid
       generated using the Delft3d GUI 
       The function assumes 
       - a flexible mesh netcdf file ending with <*_net.nc> as output by the GUI
       - geographic coordinates (WGS84)

       Use some function to convert the network matrix to triangle
       https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/OpenEarthTools/openearthtools/io/dflowfm/
    
       *** functional, but for now the delft3d gui doesnt allow interpolating depths to the grid, so cannot 
       directly output zr in this function. 
       A temporary workaround is to input a 2-entry list as fname e.g ['d3d_grid_net.nc','selfe.nc']
       - the first for the nc-delft3d file
       - the second for another file with lon,lat,dep infos - can be 
                - a netcdf file with lon,lat,dep infos e.g. SELFE grid, or ROMS, or 
                - a text file with 3 columns [lon,lat,dep]
       ***  

    """
    import dflowfm_io
    dataset = netCDF4.Dataset(fname[0], 'r')
    
    lonr = dataset.variables['NetNode_x'][:]
    latr = dataset.variables['NetNode_y'][:]
    elem_network = dataset.variables['NetElemNode'][:,:]
    # shift node id to python convention i.e. first node is lonr[0] and NOT lonr[1]
    
    elem_network = elem_network -1 # Not needed  > handled within  dflowfm_io.patch2tri

    elem_network = numpy.where(elem_network[:]>=0,elem_network,-1)
    # pad to width 6  i.e.  [nb_elem x 6] to comply with dflowfm_io.patch2tri
    elem_network6 = -1 * numpy.ones([elem_network.shape[0], 6],dtype = numpy.int32)
    elem_network6[:,:4] = elem_network[:,:]

    # convert from wgs84 to cartesian coordinates
    xr,yr = convert_xy_coords(lonr,latr,epsg_code_from = 'epsg:4326',epsg_code_to = epsg_cart_str)

    # delft3d Flexible Mesh can use triangale, quad, hexagons etc.. so the network matrix can be up to 
    #  [nb_elem x 6]
    # Here we convert to [nb_elem x 3] to be consistent with triangular mesh, and easier plotting later on
    # use some tools from Deltares 

    tri,map3,ntyp = dflowfm_io.patch2tri(lonr,latr,elem_network6)
    tri  = tri + 1 # shift back to 1-based indexing for netcdf writing 
    
    # Depth interpolation, from second file---temporary hack
    import scipy.interpolate
    if fname[1].endswith('nc'): #  depth info in netcdf file
        dset_dep = netCDF4.Dataset(fname[1], 'r')
        lond = dset_dep.variables['lon'][:]
        latd = dset_dep.variables['lat'][:]
        zr = dset_dep.variables['dep'][:]
        # interpolate to receptor grid 
        #  create interpolator
        points = numpy.vstack([lond , latd]).T
        # scipy.interpolate.griddata(points, values, xi, method='linear', fill_value=nan, rescale=False)
        zr = scipy.interpolate.griddata(numpy.vstack([lond , latd]).T , zr , numpy.vstack([lonr , latr]).T, 
                                       method = 'linear',
                                       fill_value = 0.)
    else: #  depth info in  ascii-file (<.txt>,<.xyz>, etc..)
        # NOT TESTED YET
        dset_dep = numpy.loadtxt(fname[1])
        zr = scipy.interpolate.griddata(dset_dep[:,:2] , dset_dep[:,:2]  , numpy.vstack([lonr , latr]).T, 
                                       method = 'linear',
                                       fill_value = 0.)
    
    # Package to a dictionnary
    recep = dict([('lonr', lonr),('latr', latr),('xr', xr),('yr', yr),('zr', zr),('tri', tri)])

    # return lonr,latr,xr,yr,zr,tri
    return recep

def subset_particles(pmatrix,recep_grid,timestep,level2process=''):
    """ subset particles matrix based on input timestep and level2process
    # 
    # pmatrix = matrix with particle data [time id x y z etc..]
    # recep_grid = receptor grid - dictionary
    # level2process can be 
    # 'dav'  : depth-averaged 
    # 'surf' : surface layer
    # 'mid'  : mid-depth layer
    # 'bot'  : bottom layer
    #
    #  default layer thickness is 3m
    """
    layer_thick=3
    id_time = pmatrix[:,0] == t
    pmatrix_subset_time = pmatrix[id_time,:]


    if level2process == '': return pmatrix_subset_time,recep_grid['zr']# include all particles
    if level2process == 'dav': return pmatrix_subset_time,recep_grid['zr']# include all particles 
    if not id_time.any(): return pmatrix_subset_time,recep_grid['zr']#

    else: # subset based on water depth at particle location
        import scipy.interpolate
        depth_at_part = scipy.interpolate.griddata(numpy.vstack([recep_grid['xr'] , recep_grid['yr']]).T ,
                                                   recep_grid['zr'] ,
                                                   numpy.vstack([pmatrix[id_time,2] , pmatrix[id_time,3]]).T, 
                                                   method = 'linear',
                                                   fill_value = 0.)

        # depth_at_part is expected to be positive down, while ERcore output depth is negative down
        
        if level2process == 'surf': # surface layer
            id_dep = (pmatrix_subset_time[:,4] >= -layer_thick) # & (pmatrix_subset[:,5] <= ) - second bit is probably not necessary - taken care of by ERcore itsell
        elif level2process == 'mid': # mid-depth layer
            id_dep = ( (pmatrix_subset_time[:,4] >= -depth_at_part-layer_thick/2) & (pmatrix_subset_time[:,5] <= -depth_at_part+layer_thick/2) )
        elif level2process == 'bot': # bottom layer
            id_dep = (pmatrix_subset_time[:,4] <= -depth_at_part+layer_thick)

        return pmatrix_subset_time[id_dep,:],layer_thick


# MAIN---------------------------------------------------------------------------------------------------------

if __name__  ==  '__main__':

    """The function computes the suspended and deposited concentration fields 
    for a given ERcore outputs file with particles position, for a given grid of receptors
    usage : ercore_concentration.py input_ercore_filename output_netcdf_name receptor_grid_filename epsg_code_cartesian shoreline_fname
    
    -input_ercore_filename : name of ERcore ascii output file 

    -output_netcdf_name : name of output netcdf file with concentration fields

    -receptor_grid_filename : name of file containing the receptor grid
        <.mat> file : generated from a delft3d flexible mesh using ercore_make_recep_grid.m
        <.nc> file :  use the grid from that netcdf file e.g. SELFE grid
        <._net.nc> :  netcdf file with receptor grid generated by the Delft FM RGFGRID GUI
                      > the GUI doesnt allow depth interp for now 
                      so needs to add a file with depth infos , netcdf or ascii
                      see infos load_recep_d3d_nc function above - temporary hack

    -epsg_code_cartesian : epsg code for conversion of lon/lat particles position to cartesian coordinates
        see  http://www.epsg-registry.org/ for code
        NZGD2000 / New Zealand Transverse Mercator 2000 is 2193
        WGS84 is 4326

    -shoreline_fname : shoreline file <.bnd> GNOME-formatted as used in ERcore
    
    """
    # INPUTS

    # TEST USING DISPOSAL OUTPUTS, MATLAB FILE
    fname_in = '/home_old/simon/0201_lyt_ercore_out/site5/2003/ercore.dispo_site5_surface4m_class1_2003_1.out'   
    fname_out = '/home/simon/ercore_tests/concentration_computation/dispo_site5_surface4m_class1_2003_1.nc'
    receptor_grid = '/home/simon/ercore_tests/concentration_computation/grid3_net.mat'
    bnd_file = '/home/simon/ercore_tests/concentration_computation/southnz_shoreline.bnd' # this will be in lat/lon

    
    # test using dredging outputs + selfe grid from netcdf
    receptor_grid = '/home/simon/0201_ercore_lyttelton/NEW_RUNS_DISPO/channel_runs/lyt_Original_cons.nc'
    fname_in = '/home/simon/ercore_tests/concentration_computation/ercore.dredging_site5_overflow_watercolum_new_class1.out'   
    fname_out = '/home/simon/ercore_tests/concentration_computation/ercore.dredging_site5_overflow_watercolum_new_class1.nc'
    
    receptor_grid = '/home/simon/ercore_tests/concentration_computation/d3d_grid_lyt_net.nc'
    fname_out = '/home/simon/ercore_tests/concentration_computation/dispo_site5_surface4m_class1_2003_1_D3D.nc'
    receptor_grid = ['/home/simon/ercore_tests/concentration_computation/d3d_grid_lyt_net.nc',
                    '/home/simon/0201_ercore_lyttelton/NEW_RUNS_DISPO/channel_runs/lyt_Original_cons.nc']

    # fname_in = int(sys.argv[0])
    # fname_in = int(sys.argv[1])
    # receptor_grid = int(sys.argv[2])
    # bnd_file = int(sys.argv[3])

    level2process = None
    level2process = ['dav','surf','mid','bot']
    # can be 
    # 'dav'  : depth-averaged 
    # 'surf' : surface layer
    # 'mid'  : mid-depth layer
    # 'bot'  : bottom layer
    # or a combination of options ['dav','surf','bot']
    # default layer thickness is 3m



    global  epsg_code_cartesian
    global  epsg_cart_str
    epsg_code_cartesian=2193
    epsg_cart_str=('epsg:%i' % (epsg_code_cartesian))

    # load ERCORE outputs
    time,susp,depos = read_ercore_outputs(fname_in)

    # load receptor grid : recep_grid is a dictionnary
    if receptor_grid[0].endswith('mat'):  # matlab file with receptor grid info
        recep_grid = load_recep_mat(receptor_grid)

    elif receptor_grid[0].endswith('nc'): 
        path, filename = os.path.split(receptor_grid[0])
        if filename[-6:-3]=='net' : # flexible mesh generated from delft3d GUI
            recep_grid = load_recep_d3d_nc(receptor_grid) 
        else : # use the same grid as a netcdf file e.g. SELFE/SCHISM
            recep_grid = load_recep_nc(receptor_grid[0])  

    # load shoreline file, if applicable and flag points on land
    on_land,shore = read_shoreline_ercore(bnd_file,recep_grid['lonr'],recep_grid['latr'])
    
    # initialise netcdf file
    nc_dataset = netcdf_file_unstruct_init(fname_out,recep_grid,conc_level=level2process)
    
    for t_cnt,t in enumerate(time[0:5]): # time loop

        tt=netCDF4.num2date(t,units = nc_dataset.variables['time'].units)
        print tt.isoformat(' ')
        
        # Suspended particle concentration
        # 
        for lev_cnt,lev in enumerate(level2process): # level loop
            # subset particles at time t, and level lev
            susp_subset,layer_thickness = subset_particles(susp,recep_grid,t,lev)
            # compute concentration
            conc_susp,npart_susp,np_susp_total = conc_calc(susp_subset,recep_grid['xr'],recep_grid['yr'],nh=1./20.,shoreline_file=bnd_file)       
            conc_susp[recep_grid['zr']<=0] = 0. # assume zr is positive down 
            # integrate concentration of suspended particles over water depth 
            conc_susp[:] = conc_susp[:] / numpy.abs(layer_thickness) # layer_thickness is total depth if 'dav', 3m otherwise
            # write to netcdf file
            nc_dataset.variables['conc_susp'][t_cnt,lev_cnt,:] = conc_susp
            nc_dataset.variables['np_used_susp'][t_cnt,lev_cnt,:] = npart_susp
            nc_dataset.variables['np_total_susp'][t_cnt,lev_cnt,:] = np_susp_total

        # Deposited particle concentration
        # 
        # subset particles at time t, and level lev
        depos_subset,layer_thickness = subset_particles(depos,recep_grid,t) # if not lev input > subset only in time       
        # subset particles at time t, and level lev
        conc_dep,npart_depos,np_depos_total = conc_calc(depos_subset,recep_grid['xr'],recep_grid['yr'],nh=1./20.)
        conc_dep[recep_grid['zr']<=0] = 0. # assume zr is positive down
        # write to netcdf file
        nc_dataset.variables['conc_depos'][t_cnt,:] = conc_dep
        nc_dataset.variables['np_used_depos'][t_cnt,:] = npart_depos
        nc_dataset.variables['np_total_depos'][t_cnt,:]= np_depos_total          

        # write time 
        nc_dataset.variables['time'][t_cnt] = t # the time saved in ercore output files is alread CF-compliant days since 1-1-1
    
    #close netcdf file
    nc_dataset.close()

        # TO DO 
        #  >>> include option for concetration computation at different level in water colums
