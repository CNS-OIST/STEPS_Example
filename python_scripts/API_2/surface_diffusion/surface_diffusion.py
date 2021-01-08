####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   
###

# Example: Surface diffusion (with Tetexact)
# http://steps.sourceforge.net/manual/surface_diffusion.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.saving import *
from steps.rng import *

import pylab
import math
import numpy as np
import os
import sys

dirPath = os.path.dirname(os.path.abspath(__file__))
meshPath = os.path.join(dirPath, 'meshes/coin_10r_1h_13861')

# Number of iterations; plotting dt; sim endtime:
NITER = 100

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
    mdl = Model()
    with mdl:
        SA = Species.Create()
        ssys = SurfaceSystem.Create()
        with ssys:
            Diffusion(SA, DCST)
    return mdl

########################################################################

def gen_geom(mdl):
    mesh = TetMesh.Load(meshPath)

    patch_tris = [tri for tri in mesh.surface if all(v.z > 0 for v in tri.verts)]

    with mesh:
        # Create a compartment containing all tetrahedron
        cyto = TetComp.Create(mesh.tets)

        # Create the patch
        patch = TetPatch.Create(patch_tris, cyto, None, mdl.ssys)

    # Find the central tri
    ctet = mesh.tets[0, 0, 0.5e-6]
    ctri = None
    for face in ctet.faces:
        if face in patch_tris:
            ctri = face

    # Now find the distance of the centre of each tri to the central tri
    trirads = []
    triareas = []
    center = ctri.center
    for tri in patch_tris:
        trirads.append(np.linalg.norm(center - tri.center)*1.0e6)
        triareas.append(tri.Area*1.0e12)
    
    return mesh, patch_tris, ctri, trirads, triareas

########################################################################

model = gen_model()
tmgeom, patch_tris, ctri, trirads, triareas = gen_geom(model)

rng = RNG('mt19937', 512, 234) 

sim = Simulation('Tetexact', model, tmgeom, rng)

rs = ResultSelector(sim)

ACount = rs.TRIS(patch_tris).SA.Count

sim.toSave(ACount, dt=DT)

# Run NITER number of iterations:
for j in range(NITER):
    if not j%10: print("Running iteration", j)
    sim.newRun()
    sim.TRI(ctri).SA.Count = NINJECT
    sim.run(INT)

tpnts = ACount.time[0]

res_mean = np.mean(ACount.data, 0)

########################################################################

def plotres(tidx):
    if (tidx >= INT/DT):
        print("Time index is out of range.")
        return

    pylab.scatter(trirads, res_mean[tidx,:] / triareas, s=2)
    pylab.xlabel('Radial distance ($\mu$m)')            
    pylab.ylabel('Concentration (/$\mu$m$^2$)')
    t = tpnts[tidx]
    pylab.title('Unbounded surface diffusion. Time: ' + str(t) + 's')
    plotanlyt(t)
    pylab.xlim(0,10)
    pylab.ylim(0)
    pylab.show()

########################################################################

def plotanlyt(t):     
    segs = 100     
    anlytconc = pylab.zeros((segs))     
    radialds = pylab.zeros((segs))     
    maxrad = 0.0     
    for i in trirads:         
        if (i > maxrad): maxrad = i     
    maxrad *= 1e-6     
    intervals = maxrad/segs     
    rad = 0.0     
    for i in range((segs)):         
        anlytconc[i]=(NINJECT/(4*math.pi*DCST*t))* \
            (math.exp((-1.0*(rad*rad))/(4*DCST*t)))*1e-12    
        radialds[i] = rad*1e6         
        rad += intervals     
    pylab.plot(radialds, anlytconc, color = 'red')

########################################################################

plotres(20)

