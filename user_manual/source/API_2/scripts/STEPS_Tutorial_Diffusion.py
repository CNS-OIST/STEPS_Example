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
NBRUNS = 10

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

#########################
# Plotting results
#########################

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib import colormaps
import matplotlib.cm as cm
import numpy as np

plt.figure(figsize=(10, 7))

distances = XConcRand.metaData['distToCenter']
concentrs = np.mean(XConcRand.data[0:NBRUNS, -1, :], axis=0) * 1e6

plt.scatter(distances * 1e6, concentrs, s=10, label='Simulations')

def volDiffFunc(r, t):
    # Analytical value of C(r, t)
    Crt = NINJECT * np.exp(-(r**2) / (4 * DCST * t)) / (8 * (np.pi * DCST * t)**1.5)
    # Convert from number of molecules per m^3 to uM
    return Crt * (1e3 / 6.022e23)

rvals = np.linspace(0, max(distances), 100)
plt.plot(rvals * 1e6, [volDiffFunc(r, ENDT) for r in rvals], 'r', label='Theory')

plt.xlabel('Distance [μm]')
plt.ylabel('Concentration [μM]')
plt.title(f'Unbounded volume diffusion, time = {ENDT}s')
plt.legend()
plt.show()

# # # # # # # # # # # # # # # # # #

plt.figure(figsize=(10, 7))

distances  = YCountPatch.metaData['distToCenter']
counts     = np.mean(YCountPatch.data[0:NBRUNS, -1, :], axis=0)
countsPerA = counts / YCountPatch.metaData['triArea'] * 1e-12

plt.scatter(distances * 1e6, countsPerA, s=10, label='Simulations')

def surfDiffFunc(r, t):
    # Analytical value of C(r, t)
    Crt = SNINJECT * np.exp(-(r**2) / (4 * SDCST * t)) / (4 * np.pi * SDCST * t)
    # Convert from number per m^2 to number per um^2
    return Crt * 1e-12

rvals = np.linspace(0, max(distances), 100)
plt.plot(rvals * 1e6, [surfDiffFunc(r, ENDT) for r in rvals], 'r', label='Theory')

plt.xlabel('Distance [μm]')
plt.ylabel('Concentration [1/μm²]')
plt.title(f'Unbounded surface diffusion, time = {ENDT}s')
plt.legend()
plt.show()

# # # # # # # # # # # # # # # # # #

def PlotProjectedSurface(fig, ax, triangles, values, proj, xlabel, ylabel, clabel):
    cmap = colormaps['viridis']
    maxVal = max(values)

    projTris = [[v @ proj.T for v in tri.verts] for tri in triangles]
    ax.add_collection(PolyCollection(
        projTris,
        facecolor=[cmap(v / maxVal) for v in values],
        edgecolors='black',
        linewidth=0.1
    ))

    fig.colorbar(
        cm.ScalarMappable(norm=Normalize(vmin=0, vmax=maxVal),cmap=cmap),
        ax=ax, label=clabel
    )

    minPoint = np.min(projTris, axis=(0, 1))
    maxPoint = np.max(projTris, axis=(0, 1))
    ax.set_xlim(minPoint[0], maxPoint[0])
    ax.set_ylim(minPoint[1], maxPoint[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

# # # # # # # # # # # # # # # # # #

def PlotSliceConc(results, runStart, runEnd, mesh, tinds, proj, xlabel, ylabel, clabel):
    fig = plt.figure(figsize=(6, 15))

    for i, tind in enumerate(tinds):
        ax = fig.add_subplot(len(tinds), 1, i+1)
        ax.set_aspect("equal")

        triangles = TriList(results.metaData['surfTri'][::2], mesh=mesh)

        mean = np.mean(results.data[runStart:runEnd, tind, :], axis=0)
        tetCounts = np.reshape(mean, (len(mean) // 2, 2))
        tetVols = np.reshape(results.metaData['tetVol'], (len(mean) // 2, 2))
        values = np.sum(tetCounts, axis=1) / 6.02214076e23 / np.sum(tetVols, axis=1) * 1e3

        PlotProjectedSurface(fig, ax, triangles, values, proj, xlabel, ylabel, clabel)

        ax.set_title(f'time = {results.time[0][tind]}s')
        if i < len(tinds) - 1:
            ax.set_xlabel(None)
            ax.set_xticklabels([])

    plt.show()

# # # # # # # # # # # # # # # # # #

tinds = np.linspace(1, len(XCountSlice.time[0]) - 1, 3).astype(int)
proj = np.array([
    [0, 1e6, 0],
    [0, 0, 1e6]
])
xlbl = 'y position [μm]'
ylbl = 'z position [μm]'
clbl = 'X Concentration [μM]'

PlotSliceConc(XCountSlice, 0, NBRUNS, mesh, tinds, proj, xlbl, ylbl, clbl)
PlotSliceConc(XCountSlice, NBRUNS, 2*NBRUNS, mesh, tinds, proj, xlbl, ylbl, clbl)
PlotSliceConc(XCountSlice, 2*NBRUNS, 3*NBRUNS, mesh, tinds, proj, xlbl, ylbl, clbl)

# # # # # # # # # # # # # # # # # #

def PlotPatchConc(results, runStart, runEnd, mesh, tinds, proj, xlabel, ylabel, clabel):
    fig = plt.figure(figsize=(6, 15))

    for i, tind in enumerate(tinds):
        ax = fig.add_subplot(len(tinds), 1, i+1)
        ax.set_aspect("equal")

        triangles = TriList([ind for ind in results.metaData['loc_id']], mesh=mesh)
        values = np.mean(results.data[runStart:runEnd, tind, :], axis=0)
        values /= results.metaData['triArea'] * 1e-12

        PlotProjectedSurface(fig, ax, triangles, values, proj, xlabel, ylabel, clabel)

        ax.set_title(f'time = {results.time[0][tind]}s')
        if i < len(tinds) - 1:
            ax.set_xlabel(None)
            ax.set_xticklabels([])

    plt.show()

# # # # # # # # # # # # # # # # # #

proj = np.array([
    [1e6, 0, 0],
    [0, 1e6, 0]
])
xlbl = 'x position [μm]'
ylbl = 'y position [μm]'
clbl = 'Y Concentration [/μm²]'

PlotPatchConc(YCountPatch, 0, NBRUNS, mesh, tinds, proj, xlbl, ylbl, clbl)
PlotPatchConc(YCountPatch, NBRUNS, 2*NBRUNS, mesh, tinds, proj, xlbl, ylbl, clbl)
PlotPatchConc(YCountPatch, 2*NBRUNS, 3*NBRUNS, mesh, tinds, proj, xlbl, ylbl, clbl)
