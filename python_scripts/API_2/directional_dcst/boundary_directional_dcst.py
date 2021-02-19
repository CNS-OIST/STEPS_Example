#  ___license_placeholder___

import steps.interface

""" Example of directional dcst."""

from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.sim import *
from steps.utils import *

import os
import time

DCST = 0.2e-9

model = Model()
r = ReactionManager()
with model:
    SA = Species.Create()
    vsys = VolumeSystem.Create()
    with vsys:
        D_a = Diffusion.Create(SA, DCST)

dirPath = os.path.dirname(os.path.abspath(__file__))
tetFile = os.path.join(dirPath, "../../meshes/directional_dcst/mesh_tet.inp")
triFile = os.path.join(dirPath, "../../meshes/directional_dcst/mesh_tri.inp")
confFile = os.path.join(dirPath, "../../meshes/directional_dcst/mesh_conf")

mesh = TetMesh.LoadAbaqus((tetFile, triFile), scale=1e-6, shadow_mesh=confFile)

with mesh:
    comp1 = Compartment.Create(mesh.v1_tets.tets, vsys)
    comp2 = Compartment.Create(mesh.v2_tets.tets, vsys)
    
    db = DiffBoundary.Create(mesh.boundary.tris)

rng = RNG('mt19937', 512, int(time.time() % 4294967295))

sim = Simulation('Tetexact', model, mesh, rng)

print("Set directonal dcst from comp1 to comp2, and from comp2 to comp1 to 0...")

sim.newRun()
sim.comp1.SA.Count = 100
sim.comp2.SA.Count = 20
sim.db.SA.Dcst = 0.0
print("V1 Count: ", sim.comp1.SA.Count)
print("V2 Count: ", sim.comp2.SA.Count)
sim.run(1)
print("V1 Count: ", sim.comp1.SA.Count)
print("V2 Count: ", sim.comp2.SA.Count)

print("Set directonal dcst from comp1 to comp2 to 1/10 of DCST, and 0 from comp2 to comp1...")

sim.newRun()
sim.comp1.SA.Count = 100
sim.db(direc=comp2).SA.Dcst = DCST / 10
sim.db(direc=comp1).SA.Dcst = 0
print("V1 Count: ", sim.comp1.SA.Count)
print("V2 Count: ", sim.comp2.SA.Count)
sim.run(1)
print("V1 Count: ", sim.comp1.SA.Count)
print("V2 Count: ", sim.comp2.SA.Count)
