###___license_placeholder___###

""" Example of tetrahedron directional dcst."""

import random
import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as solv
from steps.utilities import meshio
import time

DCST = 0.2e-9

# define model
model = smodel.Model()
A = smodel.Spec('A', model)
volsys = smodel.Volsys('vsys', model)
D_a = smodel.Diff('D_a', volsys, A)
D_a.setDcst(DCST)

# setup geometry
mesh = meshio.importAbaqus2("mesh_tet.inp", "mesh_tri.inp", 1e-6, "mesh_conf")[0]
comp = sgeom.TmComp("comp", mesh, range(mesh.ntets))
comp.addVolsys("vsys")

# boundary triangles splitting v1 and v2
boundary_tris = mesh.getROIData("boundary")
# tetrahedrons in v1 and adjancent to the boundary
boundary_tets1 = mesh.getROIData("boundary_tets_1")
# tetrahedrons in v2 and adjancent to the boundary
boundary_tets2 = mesh.getROIData("boundary_tets_2")

# pairing their indices
pairing = {}
for tri in boundary_tris:
    neigh_tets = mesh.getTriTetNeighb(tri)
    if neigh_tets[0] in boundary_tets1:
        pairing[tri] = (neigh_tets[0], neigh_tets[1])
    else:
        pairing[tri] = (neigh_tets[1], neigh_tets[0])

rng = srng.create('mt19937', 512)
rng.initialize(int(time.time()%4294967295))

solver = solv.Tetexact(model, mesh, rng)

print "Set dcst from v1 to v2 to 0..."
for tri in pairing.keys():
    solver.setTetDiffD(pairing[tri][0], "D_a", 0, pairing[tri][1])
    # use this to get directional dcst
    #print solver.getTetDiffD(pairing[tri][0], "D_a", pairing[tri][1])
    solver.setTetCount(pairing[tri][0], "A", 10)
print "V1 Count: ", solver.getROICount("v1_tets", "A")
print "V2 Count: ", solver.getROICount("v2_tets", "A")
solver.run(1)
print "V1 Count: ", solver.getROICount("v1_tets", "A")
print "V2 Count: ", solver.getROICount("v2_tets", "A")

print "Set dcst from v1 to v2 to 1/10 of DCST..."
solver.reset()
for tri in pairing.keys():
    solver.setTetDiffD(pairing[tri][0], "D_a", DCST / 10, pairing[tri][1])
    # use this to get directional dcst
    #print solver.getTetDiffD(pairing[tri][0], "D_a", pairing[tri][1])
    solver.setTetCount(pairing[tri][0], "A", 10)
print "V1 Count: ", solver.getROICount("v1_tets", "A")
print "V2 Count: ", solver.getROICount("v2_tets", "A")
solver.run(1)
print "V1 Count: ", solver.getROICount("v1_tets", "A")
print "V2 Count: ", solver.getROICount("v2_tets", "A")

print "Set nondirectional dcst..."
solver.reset()
for tri in pairing.keys():
    solver.setTetDiffD(pairing[tri][0], "D_a", DCST)
    # use this to get nondirectional dcst
    #print solver.getTetDiffD(pairing[tri][0], "D_a")
    solver.setTetCount(pairing[tri][0], "A", 10)
print "V1 Count: ", solver.getROICount("v1_tets", "A")
print "V2 Count: ", solver.getROICount("v2_tets", "A")
solver.run(1)
print "V1 Count: ", solver.getROICount("v1_tets", "A")
print "V2 Count: ", solver.getROICount("v2_tets", "A")

