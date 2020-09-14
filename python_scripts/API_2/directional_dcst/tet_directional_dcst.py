#  ___license_placeholder___

import steps.interface

""" Example of tetrahedron directional dcst."""

import time
from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.sim import *
from steps.utils import *

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
mesh = TetMesh.LoadAbaqus(("mesh_tet.inp", "mesh_tri.inp"), scale=1e-6, shadow_mesh="mesh_conf")

with mesh:
    comp = TetComp.Create(mesh.tets, vsys)

    pairing = {}
    for tri in mesh.boundary.tris:
        if tri.tetNeighbs[0] in mesh.boundary_tets_1.tets:
            pairing[tri] = tuple(tri.tetNeighbs)
        else:
            pairing[tri] = tuple(reversed(tri.tetNeighbs))
    
rng = RNG('mt19937', 512, int(time.time() % 4294967295))

sim = Simulation('Tetexact', model, mesh, rng)


print("Set dcst from v1 to v2 to 0...")
for tri in pairing.keys():
    sim.TET(pairing[tri][0]).D_a.D = Params(0.0, direction_tet=pairing[tri][1])



    # use this to get directional dcst
    # print(solver.getTetDiffD(pairing[tri][0], "D_a", pairing[tri][1]))
    sim.TET(pairing[tri][0]).SA.Count = 10
print("V1 Count: ", sim.v1_tets.SA.Count)
print("V2 Count: ", sim.v2_tets.SA.Count)
sim.run(1)
print("V1 Count: ", sim.v1_tets.SA.Count)
print("V2 Count: ", sim.v2_tets.SA.Count)

print("Set dcst from v1 to v2 to 1/10 of DCST...")
sim.newRun()
for tri in pairing.keys():
    sim.TET(pairing[tri][0]).D_a.D = Params(dcst=DCST / 10, direction_tet=pairing[tri][1])
    # use this to get directional dcst
    # print(solver.getTetDiffD(pairing[tri][0], "D_a", pairing[tri][1]))
    sim.TET(pairing[tri][0]).SA.Count = 10
print("V1 Count: ", sim.v1_tets.SA.Count)
print("V2 Count: ", sim.v2_tets.SA.Count)
sim.run(1)
print("V1 Count: ", sim.v1_tets.SA.Count)
print("V2 Count: ", sim.v2_tets.SA.Count)

print("Set nondirectional dcst...")
sim.newRun()
for tri in pairing.keys():
    sim.TET(pairing[tri][0]).D_a.D = DCST
    # use this to get nondirectional dcst
    # print(solver.getTetDiffD(pairing[tri][0], "D_a"))
    sim.TET(pairing[tri][0]).SA.Count = 10
print("V1 Count: ", sim.v1_tets.SA.Count)
print("V2 Count: ", sim.v2_tets.SA.Count)
sim.run(1)
print("V1 Count: ", sim.v1_tets.SA.Count)
print("V2 Count: ", sim.v2_tets.SA.Count)




# import steps.interface

# #  ___license_placeholder___

# """ Example of tetrahedron directional dcst."""

# from __future__ import print_function

# # WARNING: Using a variable name that is reserved (['time']).
# import time
# from steps.geom import *
# from steps.model_Modif import *
# from steps.rng_Modif import *
# from steps.sim import *

# DCST = 0.2e-9

# #  define model
# # WARNING: Using a variable name that is reserved (['model']).
# model_Modif = Model()
# # WARNING: Using a variable name that is reserved (['r']).
# r = ReactionManager()
# with model_Modif:
    # # WARNING: Using a variable name that is reserved (['A', 'model']).
    # A_Modif = Species.Create()
    # # WARNING: Using a variable name that is reserved (['model']).
    # volsys = VolumeSystem.Create()
# with volsys, model_Modif:
    # # WARNING: Using a variable name that is reserved (['A']).
    # D_a =     Diffusion.Create(A_Modif, DCST)


# #  setup geometry
# mesh = meshio.importAbaqus2("mesh_tet.inp", "mesh_tri.inp", 1e-6, "mesh_conf")
# with mesh:
    # comp = TetComp.Create(range(len(mesh.tets)), 'vsys')

# #  boundary triangles splitting v1 and v2
# boundary_tris = mesh.getROIData("boundary")
# #  tetrahedrons in v1 and adjancent to the boundary
# boundary_tets1 = mesh.getROIData("boundary_tets_1")
# #  tetrahedrons in v2 and adjancent to the boundary
# boundary_tets2 = mesh.getROIData("boundary_tets_2")

# #  pairing their indices
# pairing = {}
# for tri in boundary_tris:
    # neigh_tets = mesh.tris[tri].tetNeighbs.indices
    # if neigh_tets[0] in boundary_tets1:
        # pairing[tri] = (neigh_tets[0], neigh_tets[1])
    # else:
        # pairing[tri] = (neigh_tets[1], neigh_tets[0])

# # WARNING: Using a variable name that is reserved (['rng']).
# rng_Modif = RNG('mt19937', 512, int(time.time() % 4294967295))

# # WARNING: Using a variable name that is reserved (['solver', 'model', 'rng']).
# solver_Modif = Simulation('Tetexact', model_Modif, mesh, rng_Modif)

# print("Set dcst from v1 to v2 to 0...")
# for tri in pairing.keys():
    # # WARNING: Using a variable name that is reserved (['solver']).
    # solver_Modif.TET(pairing[tri][0]).D = Params(diff=D_a, dcst=0.0, direction_tet=pairing[tri][1])
    # # use this to get directional dcst
    # # print(solver.getTetDiffD(pairing[tri][0], "D_a", pairing[tri][1]))
    # # WARNING: Using a variable name that is reserved (['solver']).
    # solver_Modif.TET(pairing[tri][0]).A_Modif.Count = 10
# # WARNING: Using a variable name that is reserved (['solver']).
# print("V1 Count: ", solver_Modif.v1_tets.A_Modif.Count)
# # WARNING: Using a variable name that is reserved (['solver']).
# print("V2 Count: ", solver_Modif.v2_tets.A_Modif.Count)
# # WARNING: Using a variable name that is reserved (['solver', 'run']).
# solver_Modif.run(1)
# # WARNING: Using a variable name that is reserved (['solver']).
# print("V1 Count: ", solver_Modif.v1_tets.A_Modif.Count)
# # WARNING: Using a variable name that is reserved (['solver']).
# print("V2 Count: ", solver_Modif.v2_tets.A_Modif.Count)

# print("Set dcst from v1 to v2 to 1/10 of DCST...")
# # WARNING: Using a variable name that is reserved (['solver']).
# solver_Modif.newRun()
# for tri in pairing.keys():
    # # WARNING: Using a variable name that is reserved (['solver']).
    # solver_Modif.TET(pairing[tri][0]).D = Params(diff=D_a, dcst=DCST / 10, direction_tet=pairing[tri][1])
    # # use this to get directional dcst
    # # print(solver.getTetDiffD(pairing[tri][0], "D_a", pairing[tri][1]))
    # # WARNING: Using a variable name that is reserved (['solver']).
    # solver_Modif.TET(pairing[tri][0]).A_Modif.Count = 10
# # WARNING: Using a variable name that is reserved (['solver']).
# print("V1 Count: ", solver_Modif.v1_tets.A_Modif.Count)
# # WARNING: Using a variable name that is reserved (['solver']).
# print("V2 Count: ", solver_Modif.v2_tets.A_Modif.Count)
# # WARNING: Using a variable name that is reserved (['solver', 'run']).
# solver_Modif.run(1)
# # WARNING: Using a variable name that is reserved (['solver']).
# print("V1 Count: ", solver_Modif.v1_tets.A_Modif.Count)
# # WARNING: Using a variable name that is reserved (['solver']).
# print("V2 Count: ", solver_Modif.v2_tets.A_Modif.Count)

# print("Set nondirectional dcst...")
# # WARNING: Using a variable name that is reserved (['solver']).
# solver_Modif.newRun()
# for tri in pairing.keys():
    # # WARNING: Using a variable name that is reserved (['solver']).
    # solver_Modif.TET(pairing[tri][0]).D = Params(diff=D_a, dcst=DCST)
    # # use this to get nondirectional dcst
    # # print(solver.getTetDiffD(pairing[tri][0], "D_a"))
    # # WARNING: Using a variable name that is reserved (['solver']).
    # solver_Modif.TET(pairing[tri][0]).A_Modif.Count = 10
# # WARNING: Using a variable name that is reserved (['solver']).
# print("V1 Count: ", solver_Modif.v1_tets.A_Modif.Count)
# # WARNING: Using a variable name that is reserved (['solver']).
# print("V2 Count: ", solver_Modif.v2_tets.A_Modif.Count)
# # WARNING: Using a variable name that is reserved (['solver', 'run']).
# solver_Modif.run(1)
# # WARNING: Using a variable name that is reserved (['solver']).
# print("V1 Count: ", solver_Modif.v1_tets.A_Modif.Count)
# # WARNING: Using a variable name that is reserved (['solver']).
# print("V2 Count: ", solver_Modif.v2_tets.A_Modif.Count)
