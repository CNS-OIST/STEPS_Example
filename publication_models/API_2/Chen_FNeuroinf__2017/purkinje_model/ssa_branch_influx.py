#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

import steps.interface

import time

from steps.rng import *
from steps.geom import *
from steps.sim import *
from steps.saving import *

from extra.constants import *
from extra import data_presets
import sys
import os
import pickle

try:
    _, RESULT_DIR = sys.argv
except:
    RESULT_DIR = "result_branch_influx"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
MESH_FILE = "meshes/branch.inp"
MORPH_FILE = "meshes/branch.morph"
CA_CONC_PRESET = "extra/CaConc_PRESET1"
CA_P_CURR_DATA_FILE = "extra/CaPCurr_DATASET1"

SIM_TIME = 30.0e-6

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

os.makedirs(RESULT_DIR, exist_ok=True)

########################### LOAD PRESET CALCIUM INFLUX DATA ###########################

# convert P type calcium current data to calcium influx profile

ca_curr_data = data_presets.readData(CA_P_CURR_DATA_FILE)
ca_influx_profile = data_presets.genCaInfluxProfile(ca_curr_data, roi_areas, roi_vols, PRESET_DATA_START_TIME, PRESET_DATA_START_TIME + SIM_TIME, INFLUX_UPDATE_INTERVAL)

# load preset background calcium concerntrations
with open(CA_CONC_PRESET, 'rb') as ca_conc_preset_file:
    ca_conc_preset = pickle.load(ca_conc_preset_file)

########################### SIMULATION INITIALIZATION ###########################

rng = RNG('mt19937', 512, int(time.time()))

sim = Simulation('Tetexact', mdl, mesh, rng)

rs = ResultSelector(sim)

rs1 = rs.LIST(*rois).Ca.Conc
rs1.toFile(os.path.join(RESULT_DIR, 'ssa_calcium_conc.dat.bin'))

sim.toSave(rs1, dt=DATA_RECORD_INTERVAL)

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

for roi in rois:
    sim.LIST(roi).Ca.Conc = ca_conc_preset[roi]

n_tpns = int(SIM_TIME / DATA_RECORD_INTERVAL) + 1
update_tpns = int(INFLUX_UPDATE_INTERVAL / DATA_RECORD_INTERVAL)

influx_change_tpns = []
# for each series of the profile, the first data point is the influx start time, following by the influx rate of each roi for the specific time window
for series in ca_influx_profile['Data']:
    influx_change_tpns.append(series[0] - PRESET_DATA_START_TIME)

n_influx_changes = len(influx_change_tpns)
next_influx_change_tpn = 0

############################################################################

print("Simulating model, it will take a while...")

start_time = time.time()

for l in range(n_tpns):
    sim_endtime = DATA_RECORD_INTERVAL * l
    sim.run(sim_endtime)
    
    # need to change influx rate
    if l % update_tpns == 0:
        print("update influx rates")
        for r in range(len(rois)):
            sim.LIST(rois[r]).CaInflux['fwd'].K = ca_influx_profile['Data'][next_influx_change_tpn][r + 1]
        next_influx_change_tpn += 1

time_cost = (time.time()  - start_time)

with open(RESULT_DIR + '/ssa_performance.csv', 'w') as performance_file:
    performance_file.write("Time Cost,%f" % (time_cost) + "\n")

############################################################################
# Needed to maintain backward compatibility with file formats

with open(os.path.join(RESULT_DIR, 'ssa_calcium_conc.dat'), 'w') as f:
    f.write('#Entries: Time ')
    for roi in rois:
        f.write('%s ' % (roi))
    f.write('\n')

    for t, row in zip(rs1.time[0], rs1.data[0]):
        f.write('%.6g' %(sim_endtime) + ' ')
        for val in row:
            f.write('%.6e' %(val) + ' ')
        f.write('\n')
