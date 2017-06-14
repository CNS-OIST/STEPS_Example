#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

import steps
import time
import steps.rng as srng
import steps.solver as ssolver

from extra.constants import *
from extra import data_presets
import sys
import os
import cPickle

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

r = srng.create_mt19937(512)
r.initialize(int(time.time()))

sim = ssolver.Tetexact(mdl, mesh, r)

sim.setCompConc('cyto', 'Mg', Mg_conc)

surfarea = sim.getPatchArea('memb')
pumpnbs = 6.022141e12*surfarea

sim.setPatchCount('memb', 'Pump', round(pumpnbs))
sim.setPatchCount('memb', 'CaPump', 0)

sim.setCompConc('cyto', 'iCBsf', iCBsf_conc)
sim.setCompConc('cyto', 'iCBCaf', iCBCaf_conc)
sim.setCompConc('cyto', 'iCBsCa', iCBsCa_conc)
sim.setCompConc('cyto', 'iCBCaCa', iCBCaCa_conc)

sim.setCompConc('cyto', 'CBsf', CBsf_conc)
sim.setCompConc('cyto', 'CBCaf', CBCaf_conc)
sim.setCompConc('cyto', 'CBsCa', CBsCa_conc)
sim.setCompConc('cyto', 'CBCaCa', CBCaCa_conc)

sim.setCompConc('cyto', 'PV', PV_conc)
sim.setCompConc('cyto', 'PVCa', PVCa_conc)
sim.setCompConc('cyto', 'PVMg', PVMg_conc)

############################################################################

print "Simulating model, it will take a while..."

start_time = time.time()
sim.run(SIM_TIME)
time_cost = (time.time()  - start_time)

performance_file = open(RESULT_DIR + '/ssa_performance.csv', 'w')
performance_file.write("Time Cost,%f" % (time_cost))
performance_file.write("\n")
performance_file.close()

