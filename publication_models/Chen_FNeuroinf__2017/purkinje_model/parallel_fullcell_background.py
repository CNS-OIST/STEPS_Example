#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

                                                                                                                                                  
from __future__ import print_function
import steps.mpi
import steps.utilities.geom_decompose as gd
import steps.rng as srng
import steps.mpi.solver as mpisolver

import time
from extra.constants import *
from steps.utilities import metis_support
import sys
import os
import cPickle

if len(sys.argv) == 2:
    RESULT_DIR = sys.argv[1]
else:
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
if steps.mpi.rank == 0:
    try: os.mkdir(RESULT_DIR)
    except: pass

########################### PARTITIONING ###########################
partition_file = 'meshes/partition/fullcell.metis.epart.' + str(steps.mpi.nhosts)
mpi_tet_partitions = metis_support.readPartition(partition_file)
mpi_tri_partitions = gd.partitionTris(mesh, mpi_tet_partitions, mesh.getSurfTris())

########################### CREATE SOLVER ###########################

r = srng.create_mt19937(512)
r.initialize(int(time.time()) * steps.mpi.rank)

sim = mpisolver.TetOpSplit(mdl, mesh, r, tet_hosts = mpi_tet_partitions, tri_hosts = mpi_tri_partitions)

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
if steps.mpi.rank == 0:
    print("Simulating model, it will take a while if running with small amount of processes...")

start_time = time.time()
sim.run(SIM_TIME)
time_cost = (time.time()  - start_time)

proc_file = open(RESULT_DIR + '/proc_%i.csv' % (steps.mpi.rank), 'w', 0)
proc_file.write("SimTime,CompTime,SyncTime,IdleTime,nIteration\n")
proc_file.write("%f,%f,%f,%f,%i\n" % (time_cost, sim.getCompTime(), sim.getSyncTime(), sim.getIdleTime(), sim.getNIteration()))
proc_file.close()

if steps.mpi.rank == 0:
    performance_file = open(RESULT_DIR + '/performance_%iprocs.csv' % (steps.mpi.nhosts), 'w')
    performance_file.write("Time Cost,%f" % (time_cost))
    performance_file.write("\n")
    performance_file.close()

