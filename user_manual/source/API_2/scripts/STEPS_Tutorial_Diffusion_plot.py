# Example: Diffusion in volumes and on surfaces / 5.6. Plotting the results
# http://steps.sourceforge.net/manual/API_2/STEPS_Tutorial_Diffusion.html#Plotting-the-results

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.saving import *
from steps.geom import *

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import numpy as np

DCST= 20e-12
SDCST = 15e-12

NINJECT = 10000
SNINJECT = 10000
NBRUNS = 100

XConcRand = ResultSelector.FromFile('XConcRand.dat')

TEND = XConcRand.time[0, -1]

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
plt.plot(rvals * 1e6, [volDiffFunc(r, TEND) for r in rvals], 'r', label='Theory')

plt.xlabel('Distance [μm]')
plt.ylabel('Concentration [μM]')
plt.title(f'Unbounded volume diffusion, time = {TEND}s')
plt.legend()
plt.show()

# # # # # # # # # # # # # # # # # #

YCountPatch = ResultSelector.FromFile('YCountPatch.dat')

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
plt.plot(rvals * 1e6, [surfDiffFunc(r, TEND) for r in rvals], 'r', label='Theory')

plt.xlabel('Distance [μm]')
plt.ylabel('Concentration [1/μm²]')
plt.title(f'Unbounded surface diffusion, time = {TEND}s')
plt.legend()
plt.show()

# # # # # # # # # # # # # # # # # #

def PlotProjectedSurface(fig, ax, triangles, values, proj, xlabel, ylabel, clabel):
    cmap = cm.get_cmap('viridis')
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

mesh = TetMesh.LoadAbaqus('meshes/cylinder_86k.inp', 1e-6)

XCountSlice = ResultSelector.FromFile('XCountSlice.dat')

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

YCountPatch = ResultSelector.FromFile('YCountPatch.dat')

PlotPatchConc(YCountPatch, 0, NBRUNS, mesh, tinds, proj, xlbl, ylbl, clbl)
PlotPatchConc(YCountPatch, NBRUNS, 2*NBRUNS, mesh, tinds, proj, xlbl, ylbl, clbl)
PlotPatchConc(YCountPatch, 2*NBRUNS, 3*NBRUNS, mesh, tinds, proj, xlbl, ylbl, clbl)
