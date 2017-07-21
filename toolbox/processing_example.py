#  Template file showing how to use the different post-processing classes

# from toolbox import particle_concentration
import particle_concentration

test_case = 2

if test_case ==1:
# TIME VARYING CONCENTRATION FIELDS--------------------------------------------------------------------------------------

        # TEST USING DISPOSAL OUTPUTS, MATLAB FILE
    fname_in = '/home_old/simon/0201_lyt_ercore_out/site5/2003/ercore.dispo_site5_surface4m_class1_2003_1.out'   
    fname_out = '/home/simon/ercore_tests/concentration_computation/dispo_site5_surface4m_class1_2003_1.nc'
    receptor_grid = '/home/simon/ercore_tests/concentration_computation/grid3_net.mat'
    bnd_file = '/home/simon/ercore_tests/concentration_computation/southnz_shoreline.bnd' # this will be in lat/lon


    # test using dredging outputs + selfe grid from netcdf
    receptor_grid = '/home/simon/0201_ercore_lyttelton/NEW_RUNS_DISPO/channel_runs/lyt_Original_cons.nc'
    fname_in = '/home/simon/ercore_tests/concentration_computation/ercore.dredging_site5_overflow_watercolum_new_class1.out'   
    fname_out = '/home/simon/ercore_tests/concentration_computation/ercore.dredging_site5_overflow_watercolum_new_class1.nc'



    receptor_grid = ['/home/simon/ercore_tests/concentration_computation/d3d_grid_lyt_net.nc',
                    '/home/simon/0201_ercore_lyttelton/NEW_RUNS_DISPO/channel_runs/lyt_Original_cons.nc']
    fname_out = '/home/simon/ercore_tests/concentration_computation/dispo_site5_surface4m_class1_2003_1_D3D_LEVELS.nc'

    level2process = None
    level2process = ['dav','surf','mid','bot']



    particle_concentration.ComputeConcentration(ercore_output_fname = fname_in,
                                                output_netcdf_fname = fname_out,
                                                receptor_grid_fname = receptor_grid,
                                                epsg_code_cartesian = 2193,
                                                levels_to_process = level2process,
                                                shoreline_fname = bnd_file,
                                                neighborhood_ratio = 1./20,
                                                layer_thickness = 3.0)

elif test_case == 2 :
# PROBABILISTIC/COMBINED CONCENTRATION FIELDS--------------------------------------------------------------------------------------


    pth = '/home/simon/ercore_tests/concentration_computation/'
    fname_in = pth + 'ercore.dredging_site5_overflow_watercolum_new_class1_*.out' # note the start symbol to indicate all months shoul be loaded
    receptor_grid = ['/home/simon/ercore_tests/concentration_computation/d3d_grid_lyt_net.nc',
                    '/home/simon/0201_ercore_lyttelton/NEW_RUNS_DISPO/channel_runs/lyt_Original_cons.nc']
    fname_out = '/home/simon/ercore_tests/concentration_computation/dispo_site5_surface4m_ALL_CLASSES_PROBABILISTIC.nc'
    bnd_file = '/home/simon/ercore_tests/concentration_computation/southnz_shoreline.bnd' # this will be in lat/lon

    level2process = None
    level2process = ['dav','surf','mid','bot']


    particle_concentration.ComputeProbabilisticConcentration(ercore_output_fname = fname_in,
                                                output_netcdf_fname = fname_out,
                                                receptor_grid_fname = receptor_grid,
                                                epsg_code_cartesian = 2193,
                                                levels_to_process = level2process,
                                                shoreline_fname = bnd_file,
                                                neighborhood_ratio = 1./20,
                                                layer_thickness = 3.0)