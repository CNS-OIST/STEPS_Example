#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

import steps.interface

import time

from steps.rng import *
from steps.sim import *
from steps.geom import *

from extra.constants import *
import sys
import os

try:
    _, RESULT_DIR = sys.argv
except:
    RESULT_DIR = "result_fullcell_background"

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
MESH_FILE = "meshes/fullcell.inp"

SIM_TIME = 30.0e-5

########################### GET BIOCHEMICAL MODEL ###############################
import CaBurst_model

mdl = CaBurst_model.getModel()

########################### MESH & BRANCH MAPPING ###########################

import CaBurst_geom
mesh = CaBurst_geom.getGeom(MESH_FILE)

########################### Recording ###########################

os.makedirs(RESULT_DIR, exist_ok=True)

########################### PARTITIONING ###########################
partition_file = 'meshes/partition/fullcell.metis.epart.' + str(MPI.nhosts)
partition = MetisPartition(mesh, partition_file, default_tris=mesh.surface)

########################### CREATE SOLVER ###########################

rng = RNG('mt19937', 512, int(time.time() * MPI.rank))

sim = Simulation('TetOpSplit', mdl, mesh, rng, partition)

sim.newRun()

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

if MPI.rank == 0:
    print("Simulating model, it will take a while if running with small amount of processes...")

start_time = time.time()
sim.run(SIM_TIME)
time_cost = (time.time()  - start_time)

with open(RESULT_DIR + '/proc_%i.csv' % MPI.rank, 'w', 1) as proc_file:
    proc_file.write("SimTime,CompTime,SyncTime,IdleTime,nIteration\n")
    proc_file.write("%f,%f,%f,%f,%i\n" % (time_cost, sim.solver.getCompTime(), sim.solver.getSyncTime(), sim.solver.getIdleTime(), sim.solver.getNIteration()))

if MPI.rank == 0:
    with open(RESULT_DIR + '/performance_%iprocs.csv' % MPI.nhosts, 'w') as performance_file:
        performance_file.write("Time Cost,%f" % (time_cost))
        performance_file.write("\n")

