#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

import datetime
import steps
import steps.mpi
import steps.model as smod
import steps.utilities.geom_decompose as gd
import steps.mpi.solver as parallel_solver
import math
import steps.utilities.meshio as meshio
import steps.geom as stetmesh
import steps.rng as srng
import time
import numpy
import os
import sys

if len(sys.argv) == 4:
    MOLECULE_RATIO = float(sys.argv[1])
    MESHFILE = sys.argv[2]
    RESULT_DIR = sys.argv[3]
else:
    MOLECULE_RATIO = 1.0
    MESHFILE = '10x10x100_3363tets.inp'
    RESULT_DIR = "result"

# The initial molecule counts
N0A = 1000
N0B = 2000
N0C = 3000
N0D = 4000
N0E = 5000
N0F = 6000
N0G = 7000
N0H = 8000
N0I = 9000
N0J = 10000

ENDTIME = 20.0


########################################################################
# Biochemical Model

def gen_model():
    
    mdl = smod.Model()
    
    # The chemical species
    A = smod.Spec('A', mdl)
    B = smod.Spec('B', mdl)
    C = smod.Spec('C', mdl)
    D = smod.Spec('D', mdl)
    E = smod.Spec('E', mdl)
    F = smod.Spec('F', mdl)
    G = smod.Spec('G', mdl)
    H = smod.Spec('H', mdl)
    I = smod.Spec('I', mdl)
    J = smod.Spec('J', mdl)

    volsys = smod.Volsys('vsys',mdl)


    R1 = smod.Reac('R1', volsys, lhs = [A, B], rhs = [C],  kcst = 1000.0e6)
    R2 = smod.Reac('R2', volsys, lhs = [C],  rhs = [A,B], kcst = 100)
    R3 = smod.Reac('R3', volsys, lhs = [C, D], rhs = [E], kcst = 100e6)
    R4 = smod.Reac('R4', volsys, lhs = [E], rhs = [C,D], kcst = 10)

    R5 = smod.Reac('R5', volsys, lhs = [F, G], rhs = [H], kcst = 10e6)
    R6 = smod.Reac('R6', volsys, lhs = [H], rhs = [F,G], kcst = 1)
    R7 = smod.Reac('R7', volsys, lhs = [H, I], rhs = [J],  kcst = 1e6)
    R8 = smod.Reac('R8', volsys, lhs = [J],  rhs = [H,I], kcst = 0.1*10)


    # The diffusion rules
    D1 = smod.Diff('D1', volsys, A,  100e-12)
    D2 = smod.Diff('D2', volsys, B,  90e-12)
    D3 = smod.Diff('D3', volsys, C, 80e-12)
    D4 = smod.Diff('D4', volsys, D, 70e-12)
    D5 = smod.Diff('D5', volsys, E, 60e-12)
    D6 = smod.Diff('D6', volsys, F,  50e-12)
    D7 = smod.Diff('D7', volsys, G,  40e-12)
    D8 = smod.Diff('D8', volsys, H,  30e-12)
    D9 = smod.Diff('D9', volsys, I,  20e-12)
    D10 = smod.Diff('D10', volsys, J, 10e-12)
    
    return mdl

########################################################################
# Geometry
def gen_geom():
    mesh = meshio.importAbaqus("meshes/" + MESHFILE, 1e-6)[0]
    ntets = mesh.countTets()
    comp = stetmesh.TmComp('comp', mesh, range(ntets))
    comp.addVolsys('vsys')
    
    return mesh

########################################################################

m = gen_model()
g = gen_geom()

########################################################################
rng = srng.create('mt19937', 512)
rng.initialize(int(time.time()%4294967295)) # The max unsigned long

# Partitioning
tet_hosts = gd.linearPartition(g, [1,5,steps.mpi.nhosts/5])

sim = parallel_solver.TetOpSplit(m, g, rng, False, tet_hosts)

########################################################################
# recording
sim_result_dir = RESULT_DIR + "/%e_%s_%ihosts" % (MOLECULE_RATIO, MESHFILE, steps.mpi.nhosts)

if steps.mpi.rank == 0:
    try: os.mkdir(RESULT_DIR)
    except: pass
    sim_result_dir = RESULT_DIR + "/%e_%s_%ihosts" % (MOLECULE_RATIO, MESHFILE, steps.mpi.nhosts)
    try: os.mkdir(sim_result_dir)
    except: pass
    summary_file = open(sim_result_dir + "/result.csv", 'w', 1)
    summary_file.write("Simulation Time,")

########################################################################

# Set initial conditions
sim.setCompCount('comp', 'A', N0A * MOLECULE_RATIO)
sim.setCompCount('comp', 'B', N0B * MOLECULE_RATIO)
sim.setCompCount('comp', 'C', N0C * MOLECULE_RATIO)
sim.setCompCount('comp', 'D', N0D * MOLECULE_RATIO)
sim.setCompCount('comp', 'E', N0E * MOLECULE_RATIO)
sim.setCompCount('comp', 'F', N0F * MOLECULE_RATIO)
sim.setCompCount('comp', 'G', N0G * MOLECULE_RATIO)
sim.setCompCount('comp', 'H', N0H * MOLECULE_RATIO)
sim.setCompCount('comp', 'I', N0I * MOLECULE_RATIO)
sim.setCompCount('comp', 'J', N0J * MOLECULE_RATIO)
start_time = time.time()
sim.run(ENDTIME)
end_time = (time.time()  - start_time)
proc_file = open(sim_result_dir + "/proc_%i.csv" % (steps.mpi.rank), 'w', 1)
proc_file.write("SimTime,CompTime,SyncTime,IdleTime,nIteration\n")
proc_file.write("%f,%f,%f,%f,%i\n" % (end_time, sim.getCompTime(), sim.getSyncTime(), sim.getIdleTime(), sim.getNIteration()))
proc_file.close()
if steps.mpi.rank == 0:
    summary_file.write("%f\n" % (end_time))
    summary_file.close()

