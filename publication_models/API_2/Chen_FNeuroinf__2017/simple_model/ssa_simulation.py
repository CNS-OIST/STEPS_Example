#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

import steps.interface

from steps.model import *
from steps.sim import *
from steps.geom import *
from steps.rng import *
import math
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
            r[1].setRates(1000000000.0, 100.0)
            r[2].setRates(100000000.0, 10.0)
            r[3].setRates(10000000.0, 1.0)
            r[4].setRates(1000000.0, 0.1*10)

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
    mesh = TetMesh.LoadAbaqus(os.path.join("meshes", MESHFILE), scale=1e-06)
    with mesh:
        comp = TetComp.Create(mesh.tets, 'volsys')
    
    return mesh

########################################################################

m = gen_model()
g = gen_geom()

########################################################################
rng = RNG('mt19937', 512, int(time.time()%4294967295))
sim = Simulation('Tetexact', m, g, rng)

########################################################################
# recording
sim_result_dir = os.path.join(RESULT_DIR,  f"ssa_{MOLECULE_RATIO}_{MESHFILE}")

os.makedirs(sim_result_dir, exist_ok=True)

summary_file = open(os.path.join(sim_result_dir, "result.csv"), 'w', 1)
summary_file.write("Simulation Time,")

########################################################################

sim.newRun()

# Set initial conditions
sim.comp.SA.Count = N0A * MOLECULE_RATIO
sim.comp.SB.Count = N0B * MOLECULE_RATIO
sim.comp.SC.Count = N0C * MOLECULE_RATIO
sim.comp.SD.Count = N0D * MOLECULE_RATIO
sim.comp.SE.Count = N0E * MOLECULE_RATIO
sim.comp.SF.Count = N0F * MOLECULE_RATIO
sim.comp.SG.Count = N0G * MOLECULE_RATIO
sim.comp.SH.Count = N0H * MOLECULE_RATIO
sim.comp.SI.Count = N0I * MOLECULE_RATIO
sim.comp.SJ.Count = N0J * MOLECULE_RATIO

start_time = time.time()

sim.run(ENDTIME)

end_time = (time.time()  - start_time)

summary_file.write("%f\n" % (end_time))
summary_file.close()

