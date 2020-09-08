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
    RESULT_DIR = "result_branch_influx"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

MESH_FILE = "meshes/branch.inp"
MORPH_FILE = "meshes/branch.morph"
CA_CONC_PRESET = "extra/CaConc_PRESET1"
CA_P_CURR_DATA_FILE = "extra/CaPCurr_DATASET1"

SIM_TIME = 30.0e-3
PRESET_DATA_START_TIME = 50.0e-3
INFLUX_UPDATE_INTERVAL= 1.0e-3
DATA_RECORD_INTERVAL = 0.02e-3

########################### GET BIOCHEMICAL MODEL ###############################
import CaBurst_model

mdl = CaBurst_model.getModel()

########################### MESH & BRANCH MAPPING ###########################

import CaBurst_geom
mesh, rois, roi_areas, roi_vols = CaBurst_geom.getGeom(MESH_FILE, MORPH_FILE)

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

########################### LOAD PRESET CALCIUM INFLUX DATA ###########################

# convert P type calcium current data to calcium influx profile
ca_curr_data = data_presets.readData(CA_P_CURR_DATA_FILE)
ca_influx_profile = data_presets.genCaInfluxProfile(ca_curr_data, roi_areas, roi_vols, PRESET_DATA_START_TIME, PRESET_DATA_START_TIME + SIM_TIME + INFLUX_UPDATE_INTERVAL, INFLUX_UPDATE_INTERVAL)

# load preset background calcium concerntrations
ca_conc_preset_file = open(CA_CONC_PRESET, 'r')
ca_conc_preset = pickle.load(ca_conc_preset_file)
ca_conc_preset_file.close()

if MPI.rank == 0:
    conc_result = open(RESULT_DIR + '/calcium_conc_%iprocs.dat' % MPI.nhosts, 'w', 1)
    conc_result.write('#Entries: Time ')
    for roi in rois:
        conc_result.write('%s ' % (roi))
    conc_result.write('\n')

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


for roi in rois:
    getattr(sim, roi).Ca.Conc = ca_conc_preset[roi]

n_tpns = int(SIM_TIME / DATA_RECORD_INTERVAL) + 1
update_tpns = int(INFLUX_UPDATE_INTERVAL / DATA_RECORD_INTERVAL)

influx_change_tpns = []
# for each series of the profile, the first data point is the influx start time, following by the influx rate of each roi for the specific time window
for series in ca_influx_profile['Data']:
    influx_change_tpns.append(series[0] - PRESET_DATA_START_TIME)

n_influx_changes = len(influx_change_tpns)
next_influx_change_tpn = 0

############################################################################
if MPI.rank == 0:
    print("Simulating model, it will take a while if running with small amount of processes...")
    print("Note: You can monitor the current simulating time point and data in your result directory, but please do not write anything to the files yourself.")

# WARNING: Using a variable name that is reserved (['time', 'time']).
start_time = time.time()

for l in range(n_tpns):
    sim_endtime = DATA_RECORD_INTERVAL * l
    # WARNING: Using a variable name that is reserved (['run']).
    sim.run(sim_endtime)
    if MPI.rank == 0:
        conc_result.write('%.6g' %(sim_endtime) + ' ')
    
    for roi in rois:
        conc = getattr(sim, roi).Ca.Conc
        if MPI.rank == 0:
            conc_result.write('%.6e' %(conc) + ' ')

    if MPI.rank == 0:
        conc_result.write('\n')

    # need to change influx rate
    if l % update_tpns == 0:
        for r in range(len(rois)):
            # WARNING: Using a variable name that is reserved (['r', 'r']).
            getattr(sim, rois[r]).CaInflux['fwd'].K = ca_influx_profile['Data'][next_influx_change_tpn][r + 1]
        next_influx_change_tpn += 1
# WARNING: Using a variable name that is reserved (['time', 'time']).
time_cost = (time.time()  - start_time)
conc_result.close()

proc_file = open(RESULT_DIR + '/proc_%i.csv' % MPI.rank, 'w', 1)
proc_file.write("SimTime,CompTime,SyncTime,IdleTime,nIteration\n")
proc_file.write("%f,%f,%f,%f,%i\n" % (time_cost, sim.getCompTime(), sim.getSyncTime(), sim.getIdleTime(), sim.getNIteration()))
proc_file.close()

if MPI.rank == 0:
    performance_file = open(RESULT_DIR + '/performance_%iprocs.csv' % MPI.nhosts, 'w')
    performance_file.write("Time Cost,%f" % (time_cost))
    performance_file.write("\n")
    performance_file.close()

