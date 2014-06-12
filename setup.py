#!/usr/bin/env python

def ext_configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('',parent_package,top_path)
    config.add_extension('ercore._flib_ercore',
                         sources = ['ftn/flib.pyf','ftn/shoreline.f90','ftn/plume.f90','ftn/slipvel.f90','ftn/rnkpar.f90','ftn/interp.f90','ftn/output.f90'],
                        )
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    k={'name':'ERcore',
      'packages':['ercore','ercore.materials','ercore.lib'],
      'version':'1.0',
      'author':'MetOcean Solutions Ltd.',
      'description':'Lagrangian particle simulation model'}
    k.update(ext_configuration(top_path='').todict())
    setup(**k)
