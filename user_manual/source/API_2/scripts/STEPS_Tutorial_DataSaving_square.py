# Example: Data recording and analysis / 6.2. Spatial Brusselator model
# http://steps.sourceforge.net/manual/API_2/STEPS_Tutorial_DataSaving.html#Spatial-Brusselator-model

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

TF = 10
CF = 5e4

#########################
# Model, geom, Simulation
#########################

DCST = 1e-13
A = 0.3
B = 3

mdl2 = Model()
r = ReactionManager()
with mdl2:
    X, Y = Species.Create()
    vsys = VolumeSystem.Create()
    with vsys:
        None <r[1]> X >r[2]> Y
        r[1].K = A * TF / CF, TF
        r[2].K = B * TF

        2*X + Y >r[3]> 3*X
        r[3].K = TF * (CF ** 2)

        Diffusion(X, DCST)
        Diffusion(Y, DCST)

ENDT = 1.5
SIM_DT = 0.01

mesh = TetMesh.LoadGmsh('../meshes/square_9k.msh', scale=1e-6)
with mesh:
    comp = Compartment.Create(mesh.tets, vsys)

rng = RNG('mt19937', 512, 1234)

sim = Simulation('Tetexact', mdl2, mesh, rng)

rs = ResultSelector(sim)

concs = rs.TETS().LIST(X, Y).Conc

sim.toSave(concs, dt=SIM_DT)

#########################
# Run simulation
#########################

with XDMFHandler('Brusselator_square') as hdf:
    sim.toDB(hdf, f'square_A{A}_B{B}_D{DCST}', A=A, B=B, DCST=DCST)

    sim.newRun()

    sim.comp.X.Conc = 5e-6
    sim.comp.Y.Conc = 1.4e-4

    sim.run(ENDT)
