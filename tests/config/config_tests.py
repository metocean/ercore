#!/usr/bin/env python
import sys,os
sys.path.insert(1,os.path.join(os.path.dirname(__file__),'..','..'))
from ercore import ERcore
from ercore.materials import *
from ercore.materials.plumes import *
from ercore.fields import *


ercore=ERcore()
# ercore.readYAML('config1.yaml',globals())
ercore.readYAML('../plume/plume_test_1.yaml',globals()) 
# this takes care of the model initialization 
# you can check materials by doing :
# ercore.materials
#

import pdb;pdb.set_trace()
print ercore

