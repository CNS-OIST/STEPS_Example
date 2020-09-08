import steps.interface

#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

                                                                                                                                                  
from __future__ import print_function
from steps.rng import *

# WARNING: Using a variable name that is reserved (['time']).
import time
from extra.constants import *
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

SIM_TIME = 30.0e-5

########################### GET BIOCHEMICAL MODEL ###############################
import CaBurst_model

mdl = CaBurst_model.getModel()

########################### MESH & BRANCH MAPPING ###########################

import CaBurst_geom
mesh = CaBurst_geom.getGeom(MESH_FILE)

########################### Recording ###########################
if MPI.rank == 0:
    try: os.mkdir(RESULT_DIR)
    except: pass

########################### PARTITIONING ###########################
partition_file = 'meshes/partition/branch.metis.epart.' + str(MPI.nhosts)
mpi_tet_partitions = metis_support.readPartition(partition_file)
# WARNING: partitionTris was incorporated into LinearMeshPartition or MetisPartition.
...

########################### CREATE SOLVER ###########################

# WARNING: Using a variable name that is reserved (['r']).
r = RNG('mt19937', 512, int(time.time()) * MPI.rank)

# WARNING: Using a variable name that is reserved (['r']).
sim = Simulation('TetOpSplit', mdl, mesh, r, MeshPartition(tet_hosts=mpi_tet_partitions, tri_hosts=mpi_tri_partitions, wm_hosts=[]))

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

# WARNING: Using a variable name that is reserved (['time', 'time']).
start_time = time.time()
# WARNING: Using a variable name that is reserved (['run']).
sim.run(SIM_TIME)
# WARNING: Using a variable name that is reserved (['time', 'time']).
time_cost = (time.time()  - start_time)

proc_file = open(RESULT_DIR + '/proc_%i.csv' % MPI.rank, 'w', 1)
proc_file.write("SimTime,CompTime,SyncTime,IdleTime,nIteration\n")
proc_file.write("%f,%f,%f,%f,%i\n" % (time_cost, sim.getCompTime(), sim.getSyncTime(), sim.getIdleTime(), sim.getNIteration()))
proc_file.close()

if MPI.rank == 0:
    performance_file = open(RESULT_DIR + '/performance_%iprocs.csv' % MPI.nhosts, 'w')
    performance_file.write("Time Cost,%f" % (time_cost))
    performance_file.write("\n")
    performance_file.close()

