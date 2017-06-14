#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################
import sys
import steps.utilities.morph_support as morph_support
import cPickle

HOC_FILE = sys.argv[1]
MORPH_FILE = sys.argv[2]

moprhdata = morph_support.hoc2morph(HOC_FILE)
morph_file = open(MORPH_FILE, 'w')
cPickle.dump(moprhdata, morph_file)
morph_file.close()
