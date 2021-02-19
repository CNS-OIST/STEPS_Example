#  This script is provided as example for STEPS user manual.
#  License: GPL2.0
#  Contact: Dr. Weiliang Chen, w.chen@oist.jp

import steps.interface

from steps.model import *
from steps.sim import *
from steps.saving import *
from steps.geom import *
from steps.rng import *
import time
import os

dirPath = os.path.dirname(os.path.abspath(__file__))
MESHFILE = os.path.join(dirPath, '../../meshes/parallel/10x10x100_3363tets.inp')
RESULT_DIR = os.path.join(dirPath, "serial_result")

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
    
    mdl = Model()
    r = ReactionManager()
    with mdl:
        # The chemical species
        SA, SB, SC, SD, SE, SF, SG, SH, SI, SJ = Species.Create()

        volsys = VolumeSystem.Create()

        with volsys:

            SA + SB <r[1]> SC
            SC + SD <r[2]> SE
            SF + SG <r[3]> SH
            SH + SI <r[4]> SJ
            r[1].K = 1000000000.0, 100.0
            r[2].K = 100000000.0, 10.0
            r[3].K = 10000000.0, 1.0
            r[4].K = 1000000.0, 0.1*10

            # The diffusion rules
            D1 =  Diffusion.Create(SA, 1e-10)
            D2 =  Diffusion.Create(SB, 9e-11)
            D3 =  Diffusion.Create(SC, 8e-11)
            D4 =  Diffusion.Create(SD, 7e-11)
            D5 =  Diffusion.Create(SE, 6e-11)
            D6 =  Diffusion.Create(SF, 5e-11)
            D7 =  Diffusion.Create(SG, 4e-11)
            D8 =  Diffusion.Create(SH, 3e-11)
            D9 =  Diffusion.Create(SI, 2e-11)
            D10 = Diffusion.Create(SJ, 1e-11)

    return mdl

########################################################################
# Geometry
def gen_geom():
    mesh = TetMesh.LoadAbaqus(MESHFILE, scale=1e-06)
    with mesh:
        comp = Compartment.Create(mesh.tets, 'volsys')
    
    return mesh

########################################################################

m = gen_model()
g = gen_geom()

########################################################################
rng = RNG('mt19937', 512, int(time.time()%4294967295))
sim = Simulation('Tetexact', m, g, rng)

########################################################################
# recording

os.makedirs(RESULT_DIR, exist_ok=True)

rs = ResultSelector(sim)

rs1 = rs.comp.LIST('SA', 'SB', 'SC', 'SD', 'SE', 'SF', 'SG', 'SH', 'SI', 'SJ').Count

rs1.toFile(os.path.join(RESULT_DIR, "result.dat"))

sim.toSave(rs1, dt=RECORDING_INTERVAL)

########################################################################

sim.newRun()

# Set initial conditions
sim.comp.SA.Count = N0A
sim.comp.SB.Count = N0B
sim.comp.SC.Count = N0C
sim.comp.SD.Count = N0D
sim.comp.SE.Count = N0E
sim.comp.SF.Count = N0F
sim.comp.SG.Count = N0G
sim.comp.SH.Count = N0H
sim.comp.SI.Count = N0I
sim.comp.SJ.Count = N0J

sim.run(ENDTIME)

# This last part is only present for backwards compatibility with the scripts created with API_1.
# We need to save to text files, like in the original script.
with open(os.path.join(RESULT_DIR, "result.csv"), 'w', 1) as f:
    f.write("Simulation Time,A,B,C,D,E,F,G,H,I,J\n")
    for t, row in zip(rs1.time[0], rs1.data[0]):
        f.write(f'{t},')
        f.write(','.join(map(str, map(int, row))) + '\n')

