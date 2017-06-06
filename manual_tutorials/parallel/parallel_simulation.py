# This script is provided as example for STEPS user manual.
# License: GPL2.0
# Contact: Dr. Weiliang Chen, w.chen@oist.jp

import steps
import steps.model as smod
import steps.utilities.meshio as meshio
import steps.geom as stetmesh
import steps.rng as srng
import time
import os

import steps.mpi
import steps.mpi.solver as parallel_solver
import steps.utilities.geom_decompose as gd

MESHFILE = '10x10x100_3363tets.inp'
RESULT_DIR = "parallel_result"

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
RECORDING_INTERVAL = 1.0

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
    mesh = meshio.importAbaqus(MESHFILE, 1e-6)[0]
    ntets = mesh.countTets()
    comp = stetmesh.TmComp('comp', mesh, range(ntets))
    comp.addVolsys('vsys')
    
    return mesh

########################################################################

m = gen_model()
g = gen_geom()

###################### Geometry Partitioning ###########################

tet_hosts = gd.linearPartition(g, [1, 1, steps.mpi.nhosts])
if steps.mpi.rank == 0:
    gd.validatePartition(g, tet_hosts)
    gd.printPartitionStat(tet_hosts)

########################################################################
# recording
if steps.mpi.rank == 0:
    try: os.mkdir(RESULT_DIR)
    except: pass

    summary_file = open(RESULT_DIR + "/result.csv", 'w', 0)
    summary_file.write("Simulation Time,A,B,C,D,E,F,G,H,I,J\n")

########################################################################

rng = srng.create('mt19937', 512)
rng.initialize(int(time.time()%4294967295))

sim = parallel_solver.TetOpSplit(m, g, rng, parallel_solver.EF_NONE, tet_hosts)

# Set initial conditions
sim.setCompCount('comp', 'A', N0A)
sim.setCompCount('comp', 'B', N0B)
sim.setCompCount('comp', 'C', N0C)
sim.setCompCount('comp', 'D', N0D)
sim.setCompCount('comp', 'E', N0E)
sim.setCompCount('comp', 'F', N0F)
sim.setCompCount('comp', 'G', N0G)
sim.setCompCount('comp', 'H', N0H)
sim.setCompCount('comp', 'I', N0I)
sim.setCompCount('comp', 'J', N0J)

n_tpns = int(ENDTIME / RECORDING_INTERVAL) + 1

for l in range(n_tpns):
    sim_endtime = RECORDING_INTERVAL * l
    sim.run(sim_endtime)
    current_simtime = sim.getTime()
    A_count = sim.getCompCount('comp', 'A')
    B_count = sim.getCompCount('comp', 'B')
    C_count = sim.getCompCount('comp', 'C')
    D_count = sim.getCompCount('comp', 'D')
    E_count = sim.getCompCount('comp', 'E')
    F_count = sim.getCompCount('comp', 'F')
    G_count = sim.getCompCount('comp', 'G')
    H_count = sim.getCompCount('comp', 'H')
    I_count = sim.getCompCount('comp', 'I')
    J_count = sim.getCompCount('comp', 'J')

    if steps.mpi.rank == 0:
        summary_file.write("%f,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i\n" % \
        (current_simtime, A_count, B_count, C_count, D_count, \
        E_count, F_count, G_count, H_count, I_count, J_count))

if steps.mpi.rank == 0:
    summary_file.close()

