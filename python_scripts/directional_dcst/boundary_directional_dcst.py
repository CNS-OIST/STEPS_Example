###___license_placeholder___###

""" Example of directional dcst."""

from __future__ import print_function
import random
import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as solv
from steps.utilities import meshio
import time

DCST = 0.2e-9

model = smodel.Model()
A = smodel.Spec('A', model)
volsys = smodel.Volsys('vsys', model)
D_a = smodel.Diff('D_a', volsys, A)
D_a.setDcst(DCST)
    
mesh = meshio.importAbaqus2("mesh_tet.inp", "mesh_tri.inp", 1e-6, "mesh_conf")[0]
    
boundary_tris = mesh.getROIData("boundary")
v1_tets = mesh.getROIData("v1_tets")
v2_tets = mesh.getROIData("v2_tets")
    
comp1 = sgeom.TmComp("comp1", mesh, v1_tets)
comp2 = sgeom.TmComp("comp2", mesh, v2_tets)
    
comp1.addVolsys("vsys")
comp2.addVolsys("vsys")
    
db = sgeom.DiffBoundary("boundary", mesh, boundary_tris)
rng = srng.create('mt19937', 512)
rng.initialize(int(time.time()%4294967295))
    
solver = solv.Tetexact(model, mesh, rng)

print("Set directonal dcst from comp1 to comp2, and from comp2 to comp1 to 0...")
solver.setCompCount("comp1", "A", 100)
solver.setCompCount("comp2", "A", 20)
solver.setDiffBoundaryDcst("boundary", "A", 0)
print("V1 Count: ", solver.getCompCount("comp1", "A"))
print("V2 Count: ", solver.getCompCount("comp2", "A"))
solver.run(1)
print("V1 Count: ", solver.getCompCount("comp1", "A"))
print("V2 Count: ", solver.getCompCount("comp2", "A"))

print("Set directonal dcst from comp1 to comp2 to 1/10 of DCST, and 0 from comp2 to comp1...")
solver.reset()
solver.setCompCount("comp1", "A", 100)
solver.setDiffBoundaryDcst("boundary", "A", DCST / 10, "comp2")
solver.setDiffBoundaryDcst("boundary", "A", 0.0, "comp1")
print("V1 Count: ", solver.getCompCount("comp1", "A"))
print("V2 Count: ", solver.getCompCount("comp2", "A"))
solver.run(1)
print("V1 Count: ", solver.getCompCount("comp1", "A"))
print("V2 Count: ", solver.getCompCount("comp2", "A"))
