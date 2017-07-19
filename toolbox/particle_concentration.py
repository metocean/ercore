#!/usr/bin/env python
# import requried libraries
# from read_result import read_ercore
import sys
import os
import numpy
import pyproj
from ercore.shoreline import Shoreline,Boundary
from _flib_ercore import shoreline
from scipy.io import loadmat
from scipy.spatial import ConvexHull
from matplotlib.path import Path
import netCDF4
import datetime

def read_and_convert_ascii_outputs(fname,epsg_code_to):
    """ TO DO : >> allow to use wildcard in fname to combine several files e.g. monthly _1,_2 etc...
    # NEW FORMAT : Time    id  x   y   z   state   age mass    zbottom elev
    # OLD FORMAT : Time id  x   y   z   state   age mass  
    """
    # epsg_code_to_str = ('epsg:%i' % (epsg_code_to))
    # epsg_code_from_str = ('epsg:4326') 
    data = numpy.loadtxt(fname,delimiter = '\t',dtype=None,skiprows = 1)
    # convert particles coordinates from geographic wgs84 to cartesian
    x,y = convert_xy_coords(data[:,2],data[:,3],epsg_code_from = 4326,epsg_code_to = epsg_code_to) # ercore outputs will always be in lon/lat WGS84 coord
    data[:,2] = x
    data[:,3] = y
    data_susp= data[data[:,5] == 1,:]
    data_depos = data[data[:,5] == -1,:]
    time = numpy.unique(data[:,0])
    return time,data_susp,data_depos

def convert_xy_coords(x,y,epsg_code_from = 2193,epsg_code_to = 4326):
    """x,y are [Nx1] arrays
    # epsg_code_from : epsg code of input x,y coordinates 
    # epsg_code_to : epsg code of coordinates system x,y will be converted to 
    # 
    # ercore outputs are in geograpgical coordinates wgs84
    # convert to cartesian coordinates using epsg code
    # default is nzdg2000=nztm / see spatialreference.org for others
    """
    # convert to string code use in pyproj
    epsg_code_from_str=('epsg:%i' % (epsg_code_from))
    epsg_code_to_str=('epsg:%i' % (epsg_code_to))

    proj_from = pyproj.Proj(init = epsg_code_from_str)  # nzgd 2000 = NZTM
    proj_to = pyproj.Proj(init = epsg_code_to_str)  # nzgd 2000 = NZTM
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

def load_recep_mat(fname,epsg_code_from = 4326,epsg_code_to = 2193):
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
    xr,yr = convert_xy_coords(lonr,latr,epsg_code_from = epsg_code_from ,epsg_code_to = epsg_code_to) 
    zr = receptors['recep']['z'][0][0][:].squeeze()
    tri = receptors['tri'][:].squeeze()

    # Package to a dictionnary
    recep = dict([('lonr', lonr),('latr', latr),('xr', xr),('yr', yr),('zr', zr),('tri', tri)])

    # return lonr,latr,xr,yr,zr,tri
    return recep

def load_recep_nc(fname,epsg_code_from = 4326,epsg_code_to = 2193):
    """function to load a netcdf-based receptor grid"""

    dataset = netCDF4.Dataset(fname, 'r')
    if 'node' in dataset.dimensions.keys(): # unstructured grid
        lonr = dataset.variables['lon'][:]
        latr = dataset.variables['lat'][:]
        # convert from wgs84 to cartesian coordinates
        xr,yr = convert_xy_coords(lonr,latr,epsg_code_from = epsg_code_from ,epsg_code_to = epsg_code_to) 

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
        xr,yr = convert_xy_coords(lonr,latr,epsg_code_from = epsg_code_from ,epsg_code_to = epsg_code_to) 
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

def load_recep_d3d_nc(fname,epsg_code_from = 4326,epsg_code_to = 2193):
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
    xr,yr = convert_xy_coords(lonr,latr,epsg_code_from = epsg_code_from ,epsg_code_to = epsg_code_to) 

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

class ComputeConcentration(object):
    """The function computes time-varying suspended and deposited concentration fields 
    for a given ERcore outputs file with particles position, for a given grid of receptors
    usage : ercore_concentration.py input_ercore_filename output_netcdf_name receptor_grid_filename epsg_code_cartesian shoreline_fname
    
    - input_ercore_filename : name of ERcore ascii output file 

    - output_netcdf_name : name of output netcdf file with concentration fields

    - receptor_grid_filename : name of file containing the receptor grid
        <.mat> file : generated from a delft3d flexible mesh using ercore_make_recep_grid.m
        <.nc> file :  use the grid from that netcdf file e.g. SELFE grid
        <._net.nc> :  netcdf file with receptor grid generated by the Delft FM RGFGRID GUI
                      > the GUI doesnt allow depth interp for now 
                      so needs to add a file with depth infos , can be netcdf or ascii
                      see infos load_recep_d3d_nc function above - temporary hack

    - levels_to_process : can be 
        'dav'  : depth-averaged 
        'surf' : surface layer
        'mid'  : mid-depth layer
        'bot'  : bottom layer
        or combinations ['dav','surf','mid','bot']
        default is 'dav'

    - epsg_code_cartesian : epsg code for conversion of lon/lat particles position to cartesian coordinates
        see  http://www.epsg-registry.org/ for code
        NZGD2000 / New Zealand Transverse Mercator 2000 is 2193
        WGS84 is 4326

    - shoreline_fname : shoreline file <.bnd> GNOME-formatted as used in ERcore

    - neighborhood_ratio : ratio of particle to keep to define neighborhood of a given receptor
                           default is 1/20

    - layer_thickness :  in meters for the levels surf,mid,bot - default is 3.0 meters

    ** could become a base class for future processing options e.g. density, connectivity
    ** class ComputeDensity(ComputeConcentration): 
    ** then change or keep subfunctions to load receptors, subset particles etc.. as needed
    ** may be relevant to rename that initial class to a more general name: e.g. ProcessingTool
    **
    ** class ComputeConcentration(ProcessingTool): 
    ** class ComputeDensity(ProcessingTool): 
    
    """  
    def __init__(self,ercore_output_fname = None,
                    output_netcdf_fname = None,
                    receptor_grid_fname = None,
                    epsg_code_cartesian = 2193,
                    levels_to_process = ['dav'],
                    shoreline_fname = None,
                    neighborhood_ratio = 1./20,
                    layer_thickness = 3.0):

        # super(ComputeConcentration, self).__init__()
        self.ercore_output_fname = ercore_output_fname
        self.output_netcdf_fname = output_netcdf_fname
        self.receptor_grid_fname = receptor_grid_fname
        self.epsg_code_native = 4326
        self.epsg_code_native_str = ('epsg:%i' % (self.epsg_code_native)) 
        self.epsg_code_cartesian = epsg_code_cartesian
        self.epsg_cart_str = ('epsg:%i' % (epsg_code_cartesian)) 

        self.levels_to_process = levels_to_process
        self.shoreline_fname = shoreline_fname
        self.neighborhood_ratio = neighborhood_ratio
        self.layer_thickness=layer_thickness
        

        # load ERCORE outputs
        self.time,self.susp,self.depos = read_and_convert_ascii_outputs(self.ercore_output_fname,self.epsg_code_cartesian)      
        # load receptor grid : recep_grid is a dictionnary
        # 
        # Matlab file with receptor grid info
        if self.receptor_grid_fname[0].endswith('mat'):  
            self.receptor_grid = load_recep_mat(self.receptor_grid,epsg_code_from = self.epsg_code_native,epsg_code_to = self.epsg_code_cartesian)
        
        if isinstance(self.receptor_grid_fname,list):
            self.path, self.filename = os.path.split(self.receptor_grid_fname[0])
        else:
            self.path, self.filename = os.path.split(self.receptor_grid_fname)
        
        # Netcdf file : delft FM grid, or use grid read from a standard nc file
        if (self.filename.endswith('nc'))  & ( self.filename[-6:-3]=='net')  : 
            # flexible mesh generated from delft3d GUI
            # temporary hack to be able to have depth for the delft FM grid - need to input two files
            # [grid,file_with_depth] see load_recep_d3d_nc
            self.receptor_grid = load_recep_d3d_nc(self.receptor_grid_fname,epsg_code_from = self.epsg_code_native,epsg_code_to = self.epsg_code_cartesian)
        else:   # use the same grid as a netcdf file e.g. SELFE/SCHISM
            self.receptor_grid = load_recep_nc(self.receptor_grid_fname , epsg_code_from = self.epsg_code_native,epsg_code_to = self.epsg_code_cartesian)  

        
        # load shoreline file, if applicable and flag points on land
        # it is expected that the shoreline point coordinates is wgs84
        self.on_land,self.shore = read_shoreline_ercore(self.shoreline_fname,self.receptor_grid['lonr'],self.receptor_grid['latr'])        
        # initialise netcdf file
        self.nc_dataset = self.netcdf_file_unstruct_init()
        # subfunction looping through time, computing concentration, and wirting to file
        self.compute_and_write_to_ncfile()
        #close netcdf file
        self.nc_dataset.close()


    def conc_calc(self,pmatrix):
        """ compute particle concentration
        # now using data from the class itself i.e. : self.* 
        # 
        # pmatrix = matrix with particle data [time id x y z etc..]
        # x_recep = x-coordinates of receptor grid - array
        # y_recep  = y-coordinates of receptor grid - array
        # neighborhood_ratio=1./20. = fraction used to define the number of particles to include in conc computation, will use nh-th closest
        # shoreline = shoreline file, optional
        """

        #load shoreline, if applicable
        if self.shoreline_fname :
            shoreline.read_shoreline(self.shoreline_fname)
            # get cartesian shoreline coordinates 
            x_shore,y_shore = convert_xy_coords(shoreline.slx,shoreline.sly,epsg_code_from=4326,epsg_code_to=self.epsg_code_cartesian)
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
        conc = numpy.zeros(self.receptor_grid['xr'].shape)
        nb_part = numpy.zeros(self.receptor_grid['xr'].shape)

        #if np == 0 or less than 3 , return null concentration (not enough to computea convex hull)
        if np <3: return conc,nb_part,np
        # define convex hull of particle cloud 
        # and flag receptors outside of it to skip them
        conv_hull = ConvexHull(numpy.vstack([xp,yp]).T)

        hull = Path(numpy.vstack([xp[conv_hull.vertices],yp[conv_hull.vertices]]).T)  # matplotlib.path.Path
        in_hull =hull.contains_points(numpy.vstack([self.receptor_grid['xr'],self.receptor_grid['yr']]).T)# path = matplotlib.path.Path(polygon)

        # loop through receptors
        for r in range(0,len(self.receptor_grid['xr'])):

            # skip if receptor is on land
            if self.on_land[r]: conc[r]=0.; continue
           # skip if receptor is outside of particle cloud convex hull
            if not in_hull[r]: conc[r]=0.; continue
        
            #compute distance of every particle from receptor
            dx = numpy.abs(xp-self.receptor_grid['xr'][r])
            dy =  numpy.abs(yp-self.receptor_grid['yr'][r]);
            dist = numpy.sqrt(dx**2+dy**2)
            #sort distance
            dist_sort = numpy.sort(dist)
            # take the neighborhood_ratio-th closest particles - the subset will be used for concentration computation
            # **numpy.argsort  returns id that would sort `dist`
            part2use = numpy.argsort(dist)[0:int(numpy.round(np*self.neighborhood_ratio))]
        
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
                if self.shoreline_fname :
                    Ldx = numpy.abs( x_shore -self.receptor_grid['xr'][r] )
                    Ldy = numpy.abs( y_shore -self.receptor_grid['yr'][r] )
                    inear = (Ldx<=lx) & (Ldy <=ly)

                    if inear.any() & (conc[r]>0.) : 
                        #create local grid bandwidth
                        x1 = self.receptor_grid['xr'][r]+lx*numpy.hstack(xxu)
                        y1 = self.receptor_grid['yr'][r]+ly*numpy.hstack(yyu)
                        # check if points are within shoreline poly
                        # convert grids coords to lon/lat before doing the shoreline check
                        xx1,yy1 = convert_xy_coords(x1,y1,epsg_code_from=self.epsg_code_cartesian,epsg_code_to=4326)
                        inland,shore = read_shoreline_ercore(self.shoreline_fname,xx1,yy1) # flag points on land
                        Ldx = numpy.abs(x1[inland]-self.receptor_grid['xr'][r])
                        Ldy = numpy.abs(y1[inland]-self.receptor_grid['yr'][r])
                        LKx = 0.75*(1-(Ldx/lx)**2)
                        LKy = 0.75*(1-(Ldy/ly)**2)
                        A = (1-0.0025*numpy.sum(numpy.sum(LKx*LKy)))
                        # apply correctionto concentration
                        conc[r]=conc[r]/A
        
        return conc,nb_part,np

    def netcdf_file_unstruct_init(self):
        """netcdf-writing subfunction for unstructured grid
           
           When there are more than one conc_level to store
           variables will have different levels in dimensions

        """
        # check that conc_level is a list
        if  not isinstance(self.levels_to_process,list) : self.levels_to_process = [self.levels_to_process]
    
        # Dimensions
        dataset = netCDF4.Dataset(self.output_netcdf_fname, 'w',  format ='NETCDF4_CLASSIC')
        # dataset = netCDF4.Dataset(fname, 'w') 
        dim_node = dataset.createDimension('node', len(self.receptor_grid['lonr']))
        # dim_lat = dataset.createDimension('lat', len(lon))
        # dim_lon = dataset.createDimension('lon', len(lat))
        dim_nb_vertex = dataset.createDimension('nb_vertex', 3) # 3 = triangle
        dim_nb_elem = dataset.createDimension('nb_elem',self.receptor_grid['tri'].shape[0]) 
        dim_time = dataset.createDimension('time', None) # use unlimited time dimension
        dim_lev = dataset.createDimension('level', len(self.levels_to_process)) # 
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
        var_np_total_susp = dataset.createVariable('np_total_susp', numpy.int32,('time','level'))

        var_conc_depos = dataset.createVariable('conc_depos', numpy.float64,('time','node'))
        var_np_used_depos = dataset.createVariable('np_used_depos', numpy.int32,('time','node'))
        var_np_total_depos  = dataset.createVariable('np_total_depos', numpy.int32,('time'))
        
        # Fill content of basic variables
        var_lon[:] = self.receptor_grid['lonr']
        var_lat[:] = self.receptor_grid['latr']
        var_x[:] = self.receptor_grid['xr']
        var_y[:] = self.receptor_grid['yr']
        var_elem[:] = self.receptor_grid['tri'] # expect tri = [NB_ELEM x 3]
        var_dep[:] = self.receptor_grid['zr']
        for lev_cnt,lev_name in enumerate(self.levels_to_process):
            # str_out = netCDF4.stringtochar(np.array(['test'], 'S4'))
            str_out = netCDF4.stringtochar(numpy.array(self.levels_to_process[lev_cnt], 'S4'))
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

    def netcdf_file_regular_init(self):
        """netcdf-writing subfunction for structured/regular grid"""
        pass

    def subset_particles(self,pmatrix,timestep,level2process=''):
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
        #  returns :
        # - subset of the ERcore particle matrix for time = timestep and level = level2process
        # - the corresponding layer thickness
        # 
        """
        id_time = pmatrix[:,0] == timestep
        pmatrix_subset_time = pmatrix[id_time,:]

        if level2process == '': return pmatrix_subset_time # include all particles
        if level2process == 'dav': return pmatrix_subset_time # include all particles 
        
        if not id_time.any(): 
            return pmatrix_subset_time,self.receptor_grid['zr']#
        else: # subset based on water depth at particle location
            import scipy.interpolate
            depth_at_part = scipy.interpolate.griddata(numpy.vstack([self.receptor_grid['xr'] , self.receptor_grid['yr']]).T ,
                                                       self.receptor_grid['zr'] ,
                                                       numpy.vstack([pmatrix[id_time,2] , pmatrix[id_time,3]]).T, 
                                                       method = 'linear',
                                                       fill_value = 0.)

            # depth_at_part is expected to be positive down, while ERcore output depth is negative down
            
            if level2process == 'surf': # surface layer
                id_dep = (pmatrix_subset_time[:,4] >= -self.layer_thickness) # & (pmatrix_subset[:,5] <= ) - second bit is probably not necessary - taken care of by ERcore itsell
            elif level2process == 'mid': # mid-depth layer
                id_dep = ( (pmatrix_subset_time[:,4] >= -depth_at_part-self.layer_thickness/2) & (pmatrix_subset_time[:,5] <= -depth_at_part+self.layer_thickness/2) )
            elif level2process == 'bot': # bottom layer
                id_dep = (pmatrix_subset_time[:,4] <= -depth_at_part+self.layer_thickness)

            return pmatrix_subset_time[id_dep,:]
    
    def compute_and_write_to_ncfile(self):

        for t_cnt,t in enumerate(self.time[0:5]): # time loop

            tt=netCDF4.num2date(t,units = self.nc_dataset.variables['time'].units)
            print tt.isoformat(' ')
            # Suspended particle concentration
            # 
            for lev_cnt,lev in enumerate(self.levels_to_process): # level loop
                # subset particles at time t, and level lev
                susp_subset = self.subset_particles(self.susp,t,lev)
                # compute concentration
                conc_susp,npart_susp,np_susp_total = self.conc_calc(susp_subset)       
                conc_susp[self.receptor_grid['zr']<=0] = 0. # assume zr is positive down 
                # integrate concentration of suspended particles over water depth
                if lev == 'dav': 
                    conc_susp[:] = conc_susp[:] / numpy.abs(self.receptor_grid['zr']) 
                else:
                    conc_susp[:] = conc_susp[:] / numpy.abs(self.layer_thickness)
                # write to netcdf file
                self.nc_dataset.variables['conc_susp'][t_cnt,lev_cnt,:] = conc_susp
                self.nc_dataset.variables['np_used_susp'][t_cnt,lev_cnt,:] = npart_susp
                self.nc_dataset.variables['np_total_susp'][t_cnt,lev_cnt] = np_susp_total

            # Deposited particle concentration
            # 
            # subset particles at time t, and level lev
            depos_subset = self.subset_particles(self.depos,t) # if not lev input > subset only in time       
            # subset particles at time t, and level lev
            conc_dep,npart_depos,np_depos_total = self.conc_calc(depos_subset)
            conc_dep[self.receptor_grid['zr']<=0] = 0. # assume zr is positive down
            # write to netcdf file
            self.nc_dataset.variables['conc_depos'][t_cnt,:] = conc_dep
            self.nc_dataset.variables['np_used_depos'][t_cnt,:] = npart_depos
            self.nc_dataset.variables['np_total_depos'][t_cnt]= np_depos_total          
            # write time 
            self.nc_dataset.variables['time'][t_cnt] = t # the time saved in ercore output files is alread CF-compliant days since 1-1-1
            # Write all buffered data to disk
            self.nc_dataset.sync()        


class ComputeProbabilisticConcentration(ComputeConcentration):
    """ComputeProbabilisticConcentration is based on the ComputeConcentration class
    This class allows computing a probabilistic concentration field based on a large number of successive 
    particle clouds
    The main difference with the ComputeConcentration class is at file/particles position loading stage
    Here we loop through a user-define set of files (using wildcards in the file name) to form a large particle cloud
    combing outputs from many files
    """

    def __init__(self,ercore_output_fname = None,
                    output_netcdf_fname = None,
                    receptor_grid_fname = None,
                    epsg_code_cartesian = 2193,
                    levels_to_process = ['dav'],
                    shoreline_fname = None,
                    neighborhood_ratio = 1./20,
                    layer_thickness = 3.0):

        # super(ComputeConcentration, self).__init__()
        self.ercore_output_fname = ercore_output_fname
        self.output_netcdf_fname = output_netcdf_fname
        self.receptor_grid_fname = receptor_grid_fname
        self.epsg_code_native = 4326
        self.epsg_code_native_str = ('epsg:%i' % (self.epsg_code_native)) 
        self.epsg_code_cartesian = epsg_code_cartesian
        self.epsg_cart_str = ('epsg:%i' % (epsg_code_cartesian)) 

        self.levels_to_process = levels_to_process
        self.shoreline_fname = shoreline_fname
        self.neighborhood_ratio = neighborhood_ratio
        self.layer_thickness=layer_thickness
        
        # get a list of all files matching the wildcard used in the input ercore_output_fname
        import glob
        filelist = glob.glob(self.ercore_output_fname)
        self.susp = []
        self.depos = []
        self.time = []

        for f in filelist:
            # load all ERCORE outputs, and appends matrices 
            t,susp,depos = read_and_convert_ascii_outputs(f,self.epsg_code_cartesian)
            self.time.append(time)
            self.susp.append(susp)
            self.depos.append(depos)
        # then continue as in ComputeConcentration.__init__

        
        # load receptor grid : recep_grid is a dictionnary
        # 
        # Matlab file with receptor grid info
        if self.receptor_grid_fname[0].endswith('mat'):  
            self.receptor_grid = load_recep_mat(self.receptor_grid,epsg_code_from = self.epsg_code_native,epsg_code_to = self.epsg_code_cartesian)
        
        if isinstance(self.receptor_grid_fname,list):
            self.path, self.filename = os.path.split(self.receptor_grid_fname[0])
        else:
            self.path, self.filename = os.path.split(self.receptor_grid_fname)
        
        # Netcdf file : delft FM grid, or use grid read from a standard nc file
        if (self.filename.endswith('nc'))  & ( self.filename[-6:-3]=='net')  : 
            # flexible mesh generated from delft3d GUI
            # temporary hack to be able to have depth for the delft FM grid - need to input two files
            # [grid,file_with_depth] see load_recep_d3d_nc
            self.receptor_grid = load_recep_d3d_nc(self.receptor_grid_fname,epsg_code_from = self.epsg_code_native,epsg_code_to = self.epsg_code_cartesian)
        else:   # use the same grid as a netcdf file e.g. SELFE/SCHISM
            self.receptor_grid = load_recep_nc(self.receptor_grid_fname , epsg_code_from = self.epsg_code_native,epsg_code_to = self.epsg_code_cartesian)  

        
        # load shoreline file, if applicable and flag points on land
        # it is expected that the shoreline point coordinates is wgs84
        self.on_land,self.shore = read_shoreline_ercore(self.shoreline_fname,self.receptor_grid['lonr'],self.receptor_grid['latr'])        
        # initialise netcdf file
        self.nc_dataset = self.netcdf_file_unstruct_init()
        # subfunction looping through time, computing concentration, and wirting to file
        self.compute_and_write_to_ncfile()
        #close netcdf file
        self.nc_dataset.close()