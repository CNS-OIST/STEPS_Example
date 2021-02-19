# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Example: Surface diffusion boundary (with Tetexact)
# http://steps.sourceforge.net/manual/surface_diffusion.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.saving import *
from steps.utils import *
from steps.rng import *

import pylab
import numpy as np
import os
import sys

dirPath = os.path.dirname(os.path.abspath(__file__))
meshPath = os.path.join(dirPath, '../../meshes/surface_diffusion_boundary/coin_10r_1h_13861')

# Number of iterations; plotting dt; sim endtime:
NITER = 10

# The data collection time increment (s)
DT = 1.0

# The simulation endtime (s)
INT = 21.0

# Number of molecules injected in centre
NINJECT = 1000	

# The diffusion constant for our diffusing species (m^2/s)
DCST = 0.08e-12


########################################################################

# Create the biochemical model
mdl = Model()
with mdl:
    X = Species.Create()
    ssys = SurfaceSystem.Create()
    with ssys:
        Diffusion(X, DCST)

########################################################################

mesh = TetMesh.Load(meshPath)

with mesh:
    # Sort patch triangles as those of positive z: A +ve x, B -ve x
    Atris, Btris = TriList(), TriList()
    for tri in mesh.surface:
        if all(v.z > 0 for v in tri.verts):
            if tri.center.x > 0:
                Atris.append(tri)
            else:
                Btris.append(tri)

    cyto = Compartment.Create(mesh.tets)
    # Create the patches
    patchA = Patch.Create(Atris, cyto, None, ssys)
    patchB = Patch.Create(Btris, cyto, None, ssys)

    # Find the set of bars that connect the two patches as the intersecting bars
    sdiffb = SDiffBoundary.Create(Atris.edges & Btris.edges, patchA, patchB)

# Find the central tri
ctri = None
allTris = Atris | Btris
for tri in mesh.tets[0, 0, 0.5e-6].faces:
    if tri in allTris:
        ctri = tri

trirads_A = [np.linalg.norm(ctri.center - tri.center)*1e6 for tri in Atris]
trirads_B = [-np.linalg.norm(ctri.center - tri.center)*1e6 for tri in Btris]
triareas_A = [tri.Area*1e12 for tri in Atris]
triareas_B = [tri.Area*1e12 for tri in Btris]

########################################################################

# Create random number generator object
rng = RNG('mt19937', 512, 1234)

# Create solver object
sim = Simulation('Tetexact', mdl, mesh, rng)

rs = ResultSelector(sim)

CountA = rs.TRIS(Atris).X.Count
CountB = rs.TRIS(Btris).X.Count

sim.toSave(CountA, CountB, dt=DT)

# Run NITER number of iterations:
for j in range(NITER):
    print("Running iteration", j)
    sim.newRun()

    sim.TRI(ctri).X.Count = NINJECT

    sim.sdiffb.X.DiffusionActive = True
    sim.sdiffb(direc=patchB).X.Dcst = 0.008e-12

    sim.run(INT)

########################################################################

def plotres(tidx, CountA, CountB):
    if (tidx >= len(CountA.time[0])):
        print("Time index is out of range.")
        return

    res_A_mean = np.mean(CountA.data[:, tidx], axis=0)
    res_B_mean = np.mean(CountB.data[:, tidx], axis=0)
    
    pylab.plot(trirads_A, res_A_mean / triareas_A, 'bo', label='patchA')
    pylab.plot(trirads_B, res_B_mean / triareas_B, 'ro',  label='patchB')

    pylab.xlabel('Radial distance ($\mu$m)')
    pylab.ylabel('Concentration (/$\mu$m$^2$)')
    pylab.xlim(-10,10)
    pylab.ylim(0)
    pylab.legend()
    pylab.show()

########################################################################

plotres(20, CountA, CountB)
