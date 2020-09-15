#  ___license_placeholder___

import steps.interface

""" Example of triangle directional dcst."""

import time
from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.sim import *

DCST = 0.2e-9

#  define model
model = Model()
r = ReactionManager()
with model:
    SA = Species.Create()
    ssys = SurfaceSystem.Create()
    with ssys:
        D_a = Diffusion.Create(SA, DCST)

DCST = 0.2e-9

mesh = TetMesh.LoadAbaqus(("mesh_tet.inp", "mesh_tri.inp"), scale=1e-6, shadow_mesh="mesh_conf")

with mesh:
    comp1 = TetComp.Create(mesh.v1_tets.tets)
    patch1 = TetPatch.Create(mesh.boundary.tris, comp1, None, ssys)

neigh_tris = mesh.neigh_tri.tris
focus_tri = mesh.focus_tri.tris[0]

rng = RNG('mt19937', 512, int(time.time() % 4294967295))

sim = Simulation('Tetexact', model, mesh, rng)

print("Set dcst from focus_tri to all neighbor tris to 0...")
sim.newRun()
for tri in neigh_tris:
    sim.TRI(focus_tri).D_a(direc=tri).D = 0
    print(sim.TRI(focus_tri).D_a(direc=tri).D)
sim.TRI(focus_tri).SA.Count = 10
print("Patch Count: ", sim.patch1.SA.Count)
print("tri Count: ", sim.TRI(focus_tri).SA.Count)
sim.run(1)
print("Patch Count: ", sim.patch1.SA.Count)
print("tri Count: ", sim.TRI(focus_tri).SA.Count)

print("Set dcst from focus_tri to all neighbor tris to 1/10 of DCST...")
sim.newRun()
for tri in neigh_tris:
    sim.TRI(focus_tri).D_a(direc=tri).D = DCST / 10
    print(sim.TRI(focus_tri).D_a(direc=tri).D)
sim.TRI(focus_tri).SA.Count = 10
print("Patch Count: ", sim.patch1.SA.Count)
print("tri Count: ", sim.TRI(focus_tri).SA.Count)
sim.run(1)
print("Patch Count: ", sim.patch1.SA.Count)
print("tri Count: ", sim.TRI(focus_tri).SA.Count)

print("Set nondirectional dcst...")
sim.newRun()
for tri in neigh_tris:
    sim.TRI(focus_tri).D_a.D = DCST
    print(sim.TRI(focus_tri).D_a.D)
sim.TRI(focus_tri).SA.Count = 10
print("Patch Count: ", sim.patch1.SA.Count)
print("tri Count: ", sim.TRI(focus_tri).SA.Count)
sim.run(1)
print("Patch Count: ", sim.patch1.SA.Count)
print("tri Count: ", sim.TRI(focus_tri).SA.Count)
