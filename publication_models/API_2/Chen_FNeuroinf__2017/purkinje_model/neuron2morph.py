#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

import steps.interface

from steps.geom import *

import sys

_, HOC_FILE, MORPH_FILE = sys.argv

morph = Morph.LoadHOC(HOC_FILE)
morph.Save(MORPH_FILE)
