####################################################################################
#
#     STEPS - STochastic Engine for Pathway Simulation
#     Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#     Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
#     See the file AUTHORS for details.
#     This file is part of STEPS.
#
#     STEPS is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License version 2,
#     as published by the Free Software Foundation.
#
#     STEPS is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################
###

#  Example: Surface diffusion (with Tetexact)
#  http://steps.sourceforge.net/manual/surface_diffusion.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from __future__ import print_function

import math

import pylab
import steps.geom as stetmesh
import steps.model as smodel
import steps.rng as srng
import steps.solver as solvmod
import steps.utilities.meshio as smeshio

#  Number of iterations; plotting dt; sim endtime:
NITER = 100

#  The data collection time increment (s)
DT = 1.0

#  The simulation endtime (s)
INT = 21.0

#  Number of molecules injected in centre
NINJECT = 1000

#  The diffusion constant for our diffusing species (m^2/s)
DCST = 0.08e-12


########################################################################


def gen_model():

    mdl = smodel.Model()
    A = smodel.Spec('A', mdl)
    ssys = smodel.Surfsys('ssys', mdl)
    diff_A = smodel.Diff('diffA', ssys, A, DCST)

    return mdl


########################################################################


def gen_geom():
    mesh = smeshio.loadMesh('meshes/coin_10r_1h_13861')[0]

    ntets = mesh.countTets()
    comp = stetmesh.TmComp('cyto', mesh, range(ntets))

    alltris = mesh.getSurfTris()

    # Sort patch triangles as those of positive z
    patch_tris = []
    for t in alltris:
        vert0, vert1, vert2 = mesh.getTri(t)
        if (
            mesh.getVertex(vert0)[2] > 0.0
            and mesh.getVertex(vert1)[2] > 0.0
            and mesh.getVertex(vert2)[2] > 0.0
        ):
            patch_tris.append(t)

    # Create the patch
    patch = stetmesh.TmPatch('patch', mesh, patch_tris, icomp=comp)
    patch.addSurfsys('ssys')

    patch_tris_n = len(patch_tris)
    trirads = pylab.zeros(patch_tris_n)
    triareas = pylab.zeros(patch_tris_n)

    # Find the central tri
    ctetidx = mesh.findTetByPoint([0.0, 0.0, 0.5e-6])
    ctet_trineighbs = mesh.getTetTriNeighb(ctetidx)
    ctri_idx = -1
    for t in ctet_trineighbs:
        if t in patch_tris:
            ctri_idx = t

    # Now find the distance of the centre of each tri to the central tri
    cbaryc = mesh.getTriBarycenter(ctri_idx)
    for i in range(patch_tris_n):
        baryc = mesh.getTriBarycenter(patch_tris[i])
        r2 = (
            math.pow((baryc[0] - cbaryc[0]), 2)
            + math.pow((baryc[1] - cbaryc[1]), 2)
            + math.pow((baryc[2] - cbaryc[2]), 2)
        )
        r = math.sqrt(r2)
        # Convert to microns
        trirads[i] = r * 1.0e6
        triareas[i] = mesh.getTriArea(patch_tris[i]) * 1.0e12

    return mesh, patch_tris, patch_tris_n, ctri_idx, trirads, triareas


########################################################################

model = gen_model()
tmgeom, patch_tris, patch_tris_n, ctri_idx, trirads, triareas = gen_geom()

rng = srng.create('mt19937', 512)
rng.initialize(234)

#  Create solver object
sim = solvmod.Tetexact(model, tmgeom, rng)

tpnts = pylab.arange(0.0, INT, DT)
ntpnts = tpnts.shape[0]

#  Create the  data structure: iterations x time points x tri samples
res = pylab.zeros((NITER, ntpnts, patch_tris_n))

#  Run NITER number of iterations:
for j in range(NITER):
    if not j % 10:
        print("Running iteration", j)
    sim.reset()
    sim.setTriCount(ctri_idx, 'A', NINJECT)
    for i in range(ntpnts):
        sim.run(tpnts[i])
        for k in range(patch_tris_n):
            res[j, i, k] = sim.getTriCount(patch_tris[k], 'A') / triareas[k]

res_mean = pylab.mean(res, axis=0)

########################################################################


def plotres(tidx):
    if tidx >= INT / DT:
        print("Time index is out of range.")
        return

    pylab.scatter(trirads, res_mean[tidx], s=2)
    pylab.xlabel('Radial distance ($\mu$m)')
    pylab.ylabel('Concentration (/$\mu$m$^2$)')
    t = tpnts[tidx]
    pylab.title('Unbounded surface diffusion. Time: ' + str(t) + 's')
    plotanlyt(t)
    pylab.xlim(0, 10)
    pylab.ylim(0)
    pylab.show()


########################################################################


def plotanlyt(t):
    segs = 100
    anlytconc = pylab.zeros((segs))
    radialds = pylab.zeros((segs))
    maxrad = 0.0
    for i in trirads:
        if i > maxrad:
            maxrad = i
    maxrad *= 1e-6
    intervals = maxrad / segs
    rad = 0.0
    for i in range((segs)):
        anlytconc[i] = (
            (NINJECT / (4 * math.pi * DCST * t))
            * (math.exp((-1.0 * (rad * rad)) / (4 * DCST * t)))
            * 1e-12
        )
        radialds[i] = rad * 1e6
        rad += intervals
    pylab.plot(radialds, anlytconc, color='red')


########################################################################

plotres(20)
