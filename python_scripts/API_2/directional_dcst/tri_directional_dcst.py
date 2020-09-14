import steps.interface

#  ___license_placeholder___

""" Example of triangle directional dcst."""

from __future__ import print_function

# WARNING: Using a variable name that is reserved (['time']).
import time
from steps.geom import *
from steps.model_Modif import *
from steps.rng_Modif import *
from steps.sim import *

DCST = 0.2e-9

#  define model
# WARNING: Using a variable name that is reserved (['model']).
model_Modif = Model()
# WARNING: Using a variable name that is reserved (['r']).
r = ReactionManager()
with model_Modif:
    # WARNING: Using a variable name that is reserved (['A', 'model']).
    A_Modif = Species.Create()
    # WARNING: Using a variable name that is reserved (['model']).
    surfsys = SurfaceSystem.Create()
with surfsys, model_Modif:
    # WARNING: Using a variable name that is reserved (['A']).
    D_a =     Diffusion.Create(A_Modif, DCST)

DCST = 0.2e-9

mesh = meshio.importAbaqus2("mesh_tet.inp", "mesh_tri.inp", 1e-6, "mesh_conf")

boundary_tris = mesh.getROIData("boundary")
v1_tets = mesh.getROIData("v1_tets")
with mesh:
    
    comp1 = TetComp.Create(v1_tets)
    
    patch1 = TetPatch.Create(boundary_tris, comp1, None, 'ssys')

neigh_tris = mesh.getROIData("neigh_tri")
focus_tri = mesh.getROIData("focus_tri")[0]

# WARNING: Using a variable name that is reserved (['rng']).
rng_Modif = RNG('mt19937', 512, int(time.time() % 4294967295))

# WARNING: Using a variable name that is reserved (['solver', 'model', 'rng']).
solver_Modif = Simulation('Tetexact', model_Modif, mesh, rng_Modif)

print("Set dcst from focus_tri to all neighbor tris to 0...")
for tri in neigh_tris:
    # WARNING: Using a variable name that is reserved (['solver']).
    solver_Modif.setTriSDiffD(focus_tri, "D_a", 0, tri)
    # WARNING: Using a variable name that is reserved (['solver']).
    print(solver_Modif.getTriSDiffD(focus_tri, "D_a", tri))
# WARNING: Using a variable name that is reserved (['solver']).
solver_Modif.TRI(focus_tri).A_Modif.Count = 10
# WARNING: Using a variable name that is reserved (['solver']).
print("Patch Count: ", solver_Modif.patch1.A_Modif.Count)
# WARNING: Using a variable name that is reserved (['solver']).
print("tri Count: ", solver_Modif.TRI(focus_tri).A_Modif.Count)
# WARNING: Using a variable name that is reserved (['solver', 'run']).
solver_Modif.run(1)
# WARNING: Using a variable name that is reserved (['solver']).
print("Patch Count: ", solver_Modif.patch1.A_Modif.Count)
# WARNING: Using a variable name that is reserved (['solver']).
print("tri Count: ", solver_Modif.TRI(focus_tri).A_Modif.Count)

print("Set dcst from focus_tri to all neighbor tris to 1/10 of DCST...")
# WARNING: Using a variable name that is reserved (['solver']).
solver_Modif.newRun()
for tri in neigh_tris:
    # WARNING: Using a variable name that is reserved (['solver']).
    solver_Modif.setTriSDiffD(focus_tri, "D_a", DCST / 10, tri)
    # WARNING: Using a variable name that is reserved (['solver']).
    print(solver_Modif.getTriSDiffD(focus_tri, "D_a", tri))
# WARNING: Using a variable name that is reserved (['solver']).
solver_Modif.TRI(focus_tri).A_Modif.Count = 10
# WARNING: Using a variable name that is reserved (['solver']).
print("Patch Count: ", solver_Modif.patch1.A_Modif.Count)
# WARNING: Using a variable name that is reserved (['solver']).
print("tri Count: ", solver_Modif.TRI(focus_tri).A_Modif.Count)
# WARNING: Using a variable name that is reserved (['solver', 'run']).
solver_Modif.run(1)
# WARNING: Using a variable name that is reserved (['solver']).
print("Patch Count: ", solver_Modif.patch1.A_Modif.Count)
# WARNING: Using a variable name that is reserved (['solver']).
print("tri Count: ", solver_Modif.TRI(focus_tri).A_Modif.Count)

print("Set nondirectional dcst...")
# WARNING: Using a variable name that is reserved (['solver']).
solver_Modif.newRun()
for tri in neigh_tris:
    # WARNING: Using a variable name that is reserved (['solver']).
    solver_Modif.setTriSDiffD(focus_tri, "D_a", DCST)
    # WARNING: Using a variable name that is reserved (['solver']).
    print(solver_Modif.getTriSDiffD(focus_tri, "D_a", tri))
# WARNING: Using a variable name that is reserved (['solver']).
solver_Modif.TRI(focus_tri).A_Modif.Count = 10
# WARNING: Using a variable name that is reserved (['solver']).
print("Patch Count: ", solver_Modif.patch1.A_Modif.Count)
# WARNING: Using a variable name that is reserved (['solver']).
print("tri Count: ", solver_Modif.TRI(focus_tri).A_Modif.Count)
# WARNING: Using a variable name that is reserved (['solver', 'run']).
solver_Modif.run(1)
# WARNING: Using a variable name that is reserved (['solver']).
print("Patch Count: ", solver_Modif.patch1.A_Modif.Count)
# WARNING: Using a variable name that is reserved (['solver']).
print("tri Count: ", solver_Modif.TRI(focus_tri).A_Modif.Count)
