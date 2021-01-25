#  ___license_placeholder___

import steps.interface

""" Example of tetrahedron directional dcst."""

from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.sim import *
from steps.utils import *

import os
import time

DCST = 0.2e-9

#  define model
model = Model()
r = ReactionManager()
with model:
    SA = Species.Create()
    vsys = VolumeSystem.Create()
    with vsys:
        D_a = Diffusion.Create(SA, DCST)

#  setup geometry
dirPath = os.path.dirname(os.path.abspath(__file__))
tetFile = os.path.join(dirPath, "mesh_tet.inp")
triFile = os.path.join(dirPath, "mesh_tri.inp")
confFile = os.path.join(dirPath, "mesh_conf")

mesh = TetMesh.LoadAbaqus((tetFile, triFile), scale=1e-6, shadow_mesh=confFile)

with mesh:
    comp = Compartment.Create(mesh.tets, vsys)

    pairing = {}
    for tri in mesh.boundary.tris:
        if tri.tetNeighbs[0] in mesh.boundary_tets_1.tets:
            pairing[tri] = tuple(tri.tetNeighbs)
        else:
            pairing[tri] = tuple(reversed(tri.tetNeighbs))
    
rng = RNG('mt19937', 512, int(time.time() % 4294967295))

sim = Simulation('Tetexact', model, mesh, rng)

print("Set dcst from v1 to v2 to 0...")
sim.newRun()
for tri in pairing.keys():
    tet1, tet2 = pairing[tri]
    sim.TET(tet1).D_a(direc=tet2).D = 0
    # use this to get directional dcst
    # print(sim.TET(tet1).D_a(direc=tet2).D)
    sim.TET(tet1).SA.Count = 10
print("V1 Count: ", sim.v1_tets.SA.Count)
print("V2 Count: ", sim.v2_tets.SA.Count)
sim.run(1)
print("V1 Count: ", sim.v1_tets.SA.Count)
print("V2 Count: ", sim.v2_tets.SA.Count)

print("Set dcst from v1 to v2 to 1/10 of DCST...")
sim.newRun()
for tri in pairing.keys():
    tet1, tet2 = pairing[tri]
    sim.TET(tet1).D_a(direc=tet2).D = DCST / 10
    # use this to get directional dcst
    # print(sim.TET(tet1).D_a(direc=tet2).D)
    sim.TET(pairing[tri][0]).SA.Count = 10
print("V1 Count: ", sim.v1_tets.SA.Count)
print("V2 Count: ", sim.v2_tets.SA.Count)
sim.run(1)
print("V1 Count: ", sim.v1_tets.SA.Count)
print("V2 Count: ", sim.v2_tets.SA.Count)

print("Set nondirectional dcst...")
sim.newRun()
for tri in pairing.keys():
    tet1, tet2 = pairing[tri]
    sim.TET(tet1).D_a.D = DCST
    # use this to get nondirectional dcst
    # print(sim.TET(tet1).D_a.D)
    sim.TET(tet1).SA.Count = 10
print("V1 Count: ", sim.v1_tets.SA.Count)
print("V2 Count: ", sim.v2_tets.SA.Count)
sim.run(1)
print("V1 Count: ", sim.v1_tets.SA.Count)
print("V2 Count: ", sim.v2_tets.SA.Count)

