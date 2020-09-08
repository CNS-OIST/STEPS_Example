import steps.interface

#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

from __future__ import print_function
import steps
# WARNING: Using a variable name that is reserved (['time']).
import time
from steps.rng import *
from steps.sim import *

from extra.constants import *
from extra import data_presets
import sys
import os
try:
    import cPickle as pickle
except:
    import pickle

if len(sys.argv) == 2:
    RESULT_DIR = sys.argv[1]
else:
    RESULT_DIR = "result_branch_background"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
MESH_FILE = "meshes/branch.inp"

SIM_TIME = 30.0e-3

########################### GET BIOCHEMICAL MODEL ###############################
import CaBurst_model

mdl = CaBurst_model.getModel()

########################### MESH & BRANCH MAPPING ###########################

import CaBurst_geom
mesh = CaBurst_geom.getGeom(MESH_FILE)

########################### Recording ###########################

try: os.mkdir(RESULT_DIR)
except: pass

########################### SIMULATION INITIALIZATION ###########################

# WARNING: Using a variable name that is reserved (['r']).
r = RNG('mt19937', 512, int(time.time()))

# WARNING: Using a variable name that is reserved (['r']).
sim = Simulation('Tetexact', mdl, mesh, r)

sim.cyto.Mg.Conc = Mg_conc

surfarea = sim.memb.Area
pumpnbs = 6.022141e12*surfarea

sim.memb.Pump.Count = round(pumpnbs)
sim.memb.CaPump.Count = 0

sim.cyto.iCBsf.Conc = iCBsf_conc
sim.cyto.iCBCaf.Conc = iCBCaf_conc
sim.cyto.iCBsCa.Conc = iCBsCa_conc
sim.cyto.iCBCaCa.Conc = iCBCaCa_conc

sim.cyto.CBsf.Conc = CBsf_conc
sim.cyto.CBCaf.Conc = CBCaf_conc
sim.cyto.CBsCa.Conc = CBsCa_conc
sim.cyto.CBCaCa.Conc = CBCaCa_conc

sim.cyto.PV.Conc = PV_conc
sim.cyto.PVCa.Conc = PVCa_conc
sim.cyto.PVMg.Conc = PVMg_conc

############################################################################

print("Simulating model, it will take a while...")

# WARNING: Using a variable name that is reserved (['time', 'time']).
start_time = time.time()
# WARNING: Using a variable name that is reserved (['run']).
sim.run(SIM_TIME)
# WARNING: Using a variable name that is reserved (['time', 'time']).
time_cost = (time.time()  - start_time)

performance_file = open(RESULT_DIR + '/ssa_performance.csv', 'w')
performance_file.write("Time Cost,%f" % (time_cost))
performance_file.write("\n")
performance_file.close()

