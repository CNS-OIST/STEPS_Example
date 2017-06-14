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

import steps.model as smodel
import steps.geom as stetmesh
import steps.utilities.meshio as smeshio
import steps.rng as srng
import steps.solver as solvmod

import pylab
import math


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

def gen_model():
    
    mdl = smodel.Model()
    X = smodel.Spec('X', mdl)
    ssys = smodel.Surfsys('ssys', mdl)
    diff_X = smodel.Diff('diffX', ssys, X,  DCST)
    
    return mdl

########################################################################

mesh = smeshio.loadMesh('meshes/coin_10r_1h_13861')[0]

ntets = mesh.countTets()
comp = stetmesh.TmComp('cyto', mesh, range(ntets))


alltris = mesh.getSurfTris()

# Sort patch triangles as those of positive z: A +ve x, B -ve x
patchA_tris = []
patchB_tris = []
patchA_bars = set()
patchB_bars = set()

for t in alltris:
    vert0, vert1, vert2 = mesh.getTri(t)
    if (mesh.getVertex(vert0)[2] > 0.0 \
        and mesh.getVertex(vert1)[2] > 0.0 \
        and mesh.getVertex(vert2)[2] > 0.0):
        if mesh.getTriBarycenter(t)[0] > 0.0:
            patchA_tris.append(t)
            bar = mesh.getTriBars(t)
            patchA_bars.add(bar[0])
            patchA_bars.add(bar[1])
            patchA_bars.add(bar[2])
        else:
            patchB_tris.append(t)
            bar = mesh.getTriBars(t)
            patchB_bars.add(bar[0])
            patchB_bars.add(bar[1])
            patchB_bars.add(bar[2])

# Create the patch
patchA = stetmesh.TmPatch('patchA', mesh, patchA_tris, icomp = comp)
patchA.addSurfsys('ssys')
patchB = stetmesh.TmPatch('patchB', mesh, patchB_tris, icomp = comp)
patchB.addSurfsys('ssys')

# Find the set of bars that connect the two patches as the intersecting bars
barsDB = patchA_bars.intersection(patchB_bars)
barsDB=list(barsDB)

# Create the surface diffusion boundary
diffb = stetmesh.SDiffBoundary('sdiffb', mesh, barsDB, (patchA, patchB))

# Find the central tri
ctetidx = mesh.findTetByPoint([0.0, 0.0, 0.5e-6])
ctet_trineighbs = mesh.getTetTriNeighb(ctetidx)
ctri_idx=-1
for t in ctet_trineighbs: 
    if t in patchA_tris+patchB_tris:
        ctri_idx = t
cbaryc = mesh.getTriBarycenter(ctri_idx)

# Record the tri radii from centre and areas for patchA and patchB
trirads_A = pylab.zeros(len(patchA_tris))
trirads_B = pylab.zeros(len(patchB_tris))
triareas_A = pylab.zeros(len(patchA_tris))
triareas_B = pylab.zeros(len(patchB_tris))

for i in range(len(patchA_tris)):
    baryc = mesh.getTriBarycenter(patchA_tris[i])
    r2 = math.pow((baryc[0]-cbaryc[0]),2) + \
            math.pow((baryc[1]-cbaryc[1]),2) + \
                math.pow((baryc[2]-cbaryc[2]),2)
    r = math.sqrt(r2)
    # Convert to microns
    trirads_A[i] = r*1.0e6
    triareas_A[i] = mesh.getTriArea(patchA_tris[i])*1.0e12

for i in range(len(patchB_tris)):
    baryc = mesh.getTriBarycenter(patchB_tris[i])
    r2 = math.pow((baryc[0]-cbaryc[0]),2) + \
            math.pow((baryc[1]-cbaryc[1]),2) + \
                math.pow((baryc[2]-cbaryc[2]),2)
    r = math.sqrt(r2)
    # Convert to microns
    trirads_B[i] = -r*1.0e6
    triareas_B[i] = mesh.getTriArea(patchB_tris[i])*1.0e12

########################################################################

# Create the biochemical model
model = gen_model()

# Create rnadom number generator object
rng = srng.create('mt19937', 512)
rng.initialize(234)

# Create solver object
sim = solvmod.Tetexact(model, mesh, rng)

# Create the simulation data structures
tpnts = pylab.arange(0.0, INT, DT)
ntpnts = tpnts.shape[0]
res_A = pylab.zeros((NITER, ntpnts, len(patchA_tris)))
res_B = pylab.zeros((NITER, ntpnts, len(patchB_tris)))

# Run NITER number of iterations:
for j in range(NITER):
    print "Running iteration", j
    sim.reset()
    sim.setTriCount(ctri_idx, 'X', NINJECT)
    sim.setSDiffBoundaryDiffusionActive('sdiffb', 'X', True)
    sim.setSDiffBoundaryDcst('sdiffb', 'X', 0.008e-12 , 'patchB')
    
    for i in range(ntpnts):
        sim.run(tpnts[i])
        for k in range(len(patchA_tris)):
            res_A[j, i, k] = sim.getTriCount(patchA_tris[k], 'X')/ \
                            triareas_A[k]
        for k in range(len(patchB_tris)):
            res_B[j, i, k] = sim.getTriCount(patchB_tris[k], 'X')/ \
                            triareas_B[k]

res_A_mean = pylab.mean(res_A, axis = 0)
res_B_mean = pylab.mean(res_B, axis = 0)

########################################################################

def plotres(tidx):
    if (tidx >= INT/DT):
        print "Time index is out of range."
        return
    
    pylab.plot(trirads_A, res_A_mean[tidx], 'bo', label='patchA')
    pylab.plot(trirads_B, res_B_mean[tidx], 'ro',  label='patchB')

    pylab.xlabel('Radial distance ($\mu$m)')
    pylab.ylabel('Concentration (/$\mu$m$^2$)')
    t = tpnts[tidx]
    pylab.xlim(-10,10)
    pylab.ylim(0)
    pylab.legend()
    pylab.show()

########################################################################

plotres(20)
