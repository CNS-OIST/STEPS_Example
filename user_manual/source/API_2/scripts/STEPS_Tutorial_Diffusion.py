# Example: Diffusion in volumes and on surfaces
# http://steps.sourceforge.net/manual/API_2/STEPS_Tutorial_Diffusion.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

import numpy as np

# Diffusion constants (m^2/s)
DCST  = 20e-12
SDCST = 15e-12

# Numbers of molecules injected in center
NINJECT  = 10000
SNINJECT = 10000

# SImulation end time
ENDT = 0.1

# Number of runs
NBRUNS = 100

# Number of elements to sample data from
NSAMPLES = 2000

#########################
# Model setup
#########################

mdl = Model()
r = ReactionManager()

with mdl:
    X, Y = Species.Create()
    vsys = VolumeSystem.Create()
    ssys = SurfaceSystem.Create()
    with vsys:
        Diffusion(X, DCST)

    with ssys:
        Diffusion(Y, SDCST)

#########################
# Geom setup
#########################

mesh = TetMesh.LoadAbaqus('../meshes/cylinder_86k.inp', scale=1e-6)

with mesh:
    topTets = TetList(tet for tet in mesh.tets if tet.center.z > 0)

    comp1 = Compartment.Create(topTets, vsys)
    comp2 = Compartment.Create(mesh.tets - comp1.tets, vsys)

    diffb = DiffBoundary.Create(comp1.surface & comp2.surface)

    patchSurf = TriList(tri for tri in mesh.surface if tri.center.z == mesh.bbox.max.z)

    leftTris = TriList(tri for tri in patchSurf if tri.center.x > 0)

    patch1 = Patch.Create(leftTris, comp1, comp2, ssys)
    patch2 = Patch.Create(patchSurf - leftTris, comp1, comp2, ssys)

    sdiffb = SDiffBoundary.Create(patch1.edges & patch2.edges, patch1, patch2)

centerTet = mesh.tets[0, 0, 0]
topTri = min(patchSurf, key=lambda tri: np.linalg.norm(tri.center))

randTets = centerTet.neighbs
while len(randTets) < NSAMPLES:
    pos = np.random.random(3) * (mesh.bbox.max - mesh.bbox.min) + mesh.bbox.min
    try:
        tet = mesh.tets[pos]
        if tet not in randTets:
            randTets.append(tet)
    except KeyError:
        continue

leftTets = TetList([tet for tet in mesh.tets if tet.center.x > 0])
diffSurf = leftTets.surface & (mesh.tets - leftTets).surface

diffSlice = TetList(mesh=mesh)
triInds = []
for tri in diffSurf:
    for tet in tri.tetNeighbs:
        diffSlice.append(tet)
        triInds.append(tri.idx)

#########################
# Geom visualization
#########################

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def plotTriangles(ax, tris, color):
    ax.add_collection(Poly3DCollection(
        [tri.verts for tri in tris],
        facecolor=color,
        edgecolors='black',
        linewidth=0.1
    ))

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection='3d')

plotTriangles(ax, mesh.surface, (0.5, 0.5, 0.5, 0.1))
plotTriangles(ax, patch1.tris,  (0.5, 1.0, 0.5, 0.5))
plotTriangles(ax, patch2.tris,  (0.5, 0.5, 1.0, 0.5))
plotTriangles(ax, [tri for tri in diffb.tris if tri.center.x > 0],  (1.0, 0.5, 0.5, 0.5))
plotTriangles(ax, [tri for tri in diffb.tris if tri.center.x <= 0],  (1.0, 0.5, 0.5, 0.5))
plotTriangles(ax, [tri for tri in diffSurf if tri.center.z > 0], (1.0, 1.0, 0.5, 0.5))
plotTriangles(ax, [tri for tri in diffSurf if tri.center.z <= 0], (1.0, 1.0, 0.5, 0.5))

ax.set_xlim(mesh.bbox.min.x, mesh.bbox.max.x)
ax.set_ylim(mesh.bbox.min.y, mesh.bbox.max.y)
ax.set_zlim(mesh.bbox.min.z, mesh.bbox.max.z)
ax.set_xlabel('x position [m]')
ax.set_ylabel('y position [m]')
ax.set_zlabel('z position [m]')
plt.show()

#########################
# Simulation setup
#########################

rng = RNG('mt19937', 512, 2903)

sim = Simulation('Tetexact', mdl, mesh, rng)

rs = ResultSelector(sim)

XConcRand   = rs.TETS(randTets).X.Conc
YCountPatch = rs.TRIS(patchSurf).Y.Count
XCountSlice = rs.TETS(diffSlice).X.Count

XConcRand.metaData['distToCenter'] = [np.linalg.norm(tet.center) for tet in randTets]
YCountPatch.metaData['distToCenter'] = [
    np.linalg.norm(tri.center[:2]) for tri in patchSurf
]
YCountPatch.metaData['triArea'] = [tri.Area for tri in patchSurf]
XCountSlice.metaData['surfTri'] = triInds
XCountSlice.metaData['tetVol'] = [tet.Vol for tet in diffSlice]

sim.toSave(XConcRand, timePoints=[ENDT])
sim.toSave(XCountSlice, YCountPatch, dt=0.01)

XCountSlice.toFile('XCountSlice.dat')
XConcRand.toFile('XConcRand.dat')
YCountPatch.toFile('YCountPatch.dat')

#########################
# Run simulation
#########################

for r in range(NBRUNS):
    sim.newRun()
    sim.TET(centerTet).X.Count = NINJECT
    sim.TRI(topTri).Y.Count = SNINJECT
    sim.diffb.X.DiffusionActive = True
    sim.sdiffb.Y.DiffusionActive = True
    sim.run(ENDT)

for r in range(NBRUNS):
    sim.newRun()
    sim.TET(centerTet).X.Count = NINJECT
    sim.TRI(topTri).Y.Count = SNINJECT
    sim.run(ENDT)

for r in range(NBRUNS):
    sim.newRun()
    sim.TET(centerTet).X.Count = NINJECT
    sim.TRI(topTri).Y.Count = SNINJECT
    sim.diffb.X.Dcst = DCST / 10
    sim.sdiffb.Y.Dcst = SDCST / 10
    sim.run(ENDT)
