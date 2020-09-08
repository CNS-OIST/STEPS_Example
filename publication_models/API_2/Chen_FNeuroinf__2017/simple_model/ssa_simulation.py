import steps.interface

#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

import datetime
import steps
from steps.model import *
from steps.sim import *
import math
from steps.geom import *
from steps.rng_Modif import *
# WARNING: Using a variable name that is reserved (['time']).
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
    # WARNING: Using a variable name that is reserved (['r']).
    r = ReactionManager()
    with mdl:
        
        # The chemical species
        
        # The chemical species
        # WARNING: Using a variable name that is reserved (['A']).
        A_Modif, B, C_Modif, D_Modif, E, F, G, H_Modif, I_Modif, J = Species.Create()

        volsys = VolumeSystem.Create()
    with volsys, mdl:


        # WARNING: Using a variable name that is reserved (['A', 'C']).
        A_Modif + B <r['R1']> C_Modif ; r['R1'].setRates(1000000000.0, 100.0)
        # WARNING: Using a variable name that is reserved (['C', 'D']).
        C_Modif + D_Modif <r['R3']> E ; r['R3'].setRates(100000000.0, 10.0)

        # WARNING: Using a variable name that is reserved (['H']).
        F + G <r['R5']> H_Modif ; r['R5'].setRates(10000000.0, 1.0)
        # WARNING: Using a variable name that is reserved (['H', 'I']).
        H_Modif + I_Modif <r['R7']> J ; r['R7'].setRates(1000000.0, 0.1*10)


        # The diffusion rules
        # WARNING: Using a variable name that is reserved (['A']).
        D1 =         Diffusion.Create(A_Modif, 1e-10)

        D2 =         Diffusion.Create(B, 9e-11)

        # WARNING: Using a variable name that is reserved (['C']).
        D3 =         Diffusion.Create(C_Modif, 8e-11)

        # WARNING: Using a variable name that is reserved (['D']).
        D4 =         Diffusion.Create(D_Modif, 7e-11)

        D5 =         Diffusion.Create(E, 6e-11)

        D6 =         Diffusion.Create(F, 5e-11)

        D7 =         Diffusion.Create(G, 4e-11)

        # WARNING: Using a variable name that is reserved (['H']).
        D8 =         Diffusion.Create(H_Modif, 3e-11)

        # WARNING: Using a variable name that is reserved (['I']).
        D9 =         Diffusion.Create(I_Modif, 2e-11)

        D10 =         Diffusion.Create(J, 1e-11)

    
    return mdl

########################################################################
# Geometry
def gen_geom():
    mesh = TetMesh.LoadAbaqus("meshes/" + MESHFILE, scale=1e-06)
    ntets = len(mesh.tets)
    with mesh:
        comp = TetComp.Create(range(ntets), 'volsys')
    
    return mesh

########################################################################

m = gen_model()
g = gen_geom()

########################################################################
# WARNING: Using a variable name that is reserved (['rng']).
rng_Modif = RNG('mt19937', 512, int(time.time()%4294967295))
# WARNING: Using a variable name that is reserved (['rng']).
sim = Simulation('Tetexact', m, g, rng_Modif)

########################################################################
# recording
sim_result_dir = RESULT_DIR + "/ssa_%e_%s" % (MOLECULE_RATIO, MESHFILE)

try: os.mkdir(RESULT_DIR)
except: pass
try: os.mkdir(sim_result_dir)
except: pass
summary_file = open(sim_result_dir + "/result.csv", 'w', 1)
summary_file.write("Simulation Time,")

########################################################################

# Set initial conditions
sim.comp.A_Modif.Count = N0A * MOLECULE_RATIO
sim.comp.B.Count = N0B * MOLECULE_RATIO
sim.comp.C_Modif.Count = N0C * MOLECULE_RATIO
sim.comp.D_Modif.Count = N0D * MOLECULE_RATIO
sim.comp.E.Count = N0E * MOLECULE_RATIO
sim.comp.F.Count = N0F * MOLECULE_RATIO
sim.comp.G.Count = N0G * MOLECULE_RATIO
sim.comp.H_Modif.Count = N0H * MOLECULE_RATIO
sim.comp.I_Modif.Count = N0I * MOLECULE_RATIO
sim.comp.J.Count = N0J * MOLECULE_RATIO
# WARNING: Using a variable name that is reserved (['time', 'time']).
start_time = time.time()
# WARNING: Using a variable name that is reserved (['run']).
sim.run(ENDTIME)
# WARNING: Using a variable name that is reserved (['time', 'time']).
end_time = (time.time()  - start_time)
summary_file.write("%f\n" % (end_time))
summary_file.close()

