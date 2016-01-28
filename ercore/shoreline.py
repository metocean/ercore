#!/usr/bin/env python
import numpy,os
from ercore import ERCoreException,ERConfigException
from _flib_ercore import shoreline


class ShorelineException(ERCoreException):
  pass

class Shoreline:
  """Class for shoreline data"""
  is3d=False
  refloat=0.
  def __init__(self,id,file,**options):
    """Constructor for Shoreline class
    Arguments:
      id: Object id
      file: Filename of shoreline file (must be in bnd format)
    """
    if not os.path.isfile(file):
      raise ShorelineException('Shoreline file %s does not exist' % (file))
    shoreline.read_shoreline(file)
    self.id=id
    self.__dict__.update(options)
    self.slx=shoreline.slx
    self.sly=shoreline.sly
    self.polyi=shoreline.polyi
    self.polyn=shoreline.polyn

  def intersect(self,pos,post,stat):
    """Check for intersection of particles with shoreline
    Arguments
      pos: Position matrix at time t
      post: Position matrix at time t+1
      stat: Matrix of particle statuses
    Returns:
      Matrix of intersection positions
    """
    poso=numpy.hstack((shoreline.sl_intersect_wrap(pos[:,:2],post[:,:2],stat,self.refloat),post[:,2:3]))
    return poso

class Boundary:
  """Class to define simulation boundary"""
  def __init__(self,id,bnd,**options):
    """Constructor for Boundary class
    Arguments:
      id: Object id
      bnd: Boundary array as [minx,maxx,miny,maxy]
    """
    self.id = id
    self.x1=bnd[0]
    self.x2=bnd[1]
    self.y1=bnd[2]
    self.y2=bnd[3]

  def intersect(self,pos,post,stat):
    """Check for intersection of particles with boundary
    Arguments
      pos: Position matrix at time t
      post: Position matrix at time t+1
      stat: Matrix of particle statuses
    Returns:
      Matrix of intersection positions
    """
    outp=post[:,:3]
    xhigh=post[:,0]>self.x2
    xlow=post[:,0]<self.x1
    yhigh=post[:,1]>self.y2
    ylow=post[:,1]<self.y1
    stat[xhigh | xlow | yhigh | ylow]=-9
    outp[xhigh,0]=self.x2
    outp[xlow,0]=self.x1
    outp[yhigh,1]=self.y2
    outp[ylow,1]=self.y1
    return outp
