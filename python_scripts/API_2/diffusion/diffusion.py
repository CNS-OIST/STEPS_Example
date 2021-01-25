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

#  Example: Unbounded diffusion
#  http://steps.sourceforge.net/manual/diffusion.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.saving import *
from steps.rng import *

import math
import numpy
import pylab
import time
import os
import sys

########################################################################

dirPath = os.path.dirname(os.path.abspath(__file__))
meshPath = os.path.join(dirPath, 'meshes/sphere_rad10_11Ktets')

# The number of iterations to run 
NITER = 10
# The data collection time increment (s)
DT = 0.001
# The simulation endtime (s)
INT = 0.101

# The number of molecules to be injected into the centre
NINJECT = 10000    

# The number of tetrahedral elements to sample data from. 
SAMPLE = 2000    

# The diffusion constant for our diffusing species (m^2/s)
DCST= 20.0e-12

########################################################################

def gen_model():         
    mdl = Model()
    with mdl:
        SA = Species.Create()
        cytosolv = VolumeSystem.Create()
        with cytosolv:
            Diffusion(SA, DCST)
    return mdl

########################################################################

def gen_geom(mdl):

    print("Loading mesh...")
    mesh = TetMesh.Load(meshPath)
    print("Mesh Loaded")

    with mesh:
        # Create a compartment containing all tetrahedron
        cyto = Compartment.Create(mesh.tets, mdl.cytosolv)

        print("Finding tetrahedron samples...")
        # List to hold tetrahedrons
        tets = TetList()

        # Fetch the central tetrahedron index and store:
        ctet = mesh.tets[0, 0, 0]
        tets.append(ctet)
        tets += ctet.neighbs

        # Find the maximum and minimum coordinates of the mesh
        bmin = mesh.bbox.min
        bmax = mesh.bbox.max

        # Run a loop until we have stored all tet indices we require
        while len(tets) < SAMPLE:
            # Pick a random position within the bounding box
            pos = bmin + (bmax - bmin) * numpy.random.random(3)
            if pos in mesh.tets:
                tets.append(mesh.tets[pos])

        # Find the radial distance of the tetrahedrons to mesh center:
        tetrads = numpy.array([numpy.linalg.norm(tet.center)*1e6 for tet in tets])

        print("Tetrahedron samples found")

    return mesh, tets, tetrads

########################################################################

model = gen_model()
mesh, tets, tetrads = gen_geom(model)

rng = RNG('mt19937', 512, 2903)

sim = Simulation('Tetexact', model, mesh, rng)

rs = ResultSelector(sim)

AConc = rs.TETS(tets).SA.Conc

sim.toSave(AConc, dt=DT)

for i in range(NITER):
    sim.newRun()
    print("Running iteration", i)
    # Inject all molecules into the central tet:
    sim.TET(0, 0, 0).SA.Count = NINJECT
    sim.run(INT)

########################################################################

def plotres(tidx, tpnts, res_mean, tetrads):
    if (tidx >= INT/DT):
        print("Time index is out of range.")
        return
    
    pylab.scatter(tetrads, res_mean[tidx,:], s=2)
    pylab.xlabel('Radial distance of tetrahedron ($\mu$m)')            
    pylab.ylabel('Concentration in tetrahedron ($\mu$M)')
    t = tpnts[tidx]
    pylab.title('Unbounded diffusion. Time: ' + str(t) + 's')
    plotanlyt(t)
    pylab.xlim(0.0, 10.0)
    pylab.ylim(0.0)
    pylab.show()

########################################################################

def plotanlyt(t):     
    segs = 100     
    anlytconc = numpy.zeros((segs))     
    radialds = numpy.zeros((segs))     
    maxrad = 0.0     
    for i in tetrads:         
        if (i > maxrad): maxrad = i     
    maxrad *= 1e-6     
    intervals = maxrad/segs     
    rad = 0.0     
    for i in range((segs)):         
        # Find the conc from analytical solution, and convert to mol/L         
        anlytconc[i]=1.0e3*(1/6.022e23)* \
                ((NINJECT/(math.pow((4*math.pi*DCST*t),1.5)))* \
                (math.exp((-1.0*(rad*rad))/(4*DCST*t))))         
        radialds[i] = rad*1e6         
        rad += intervals     
    pylab.plot(radialds, anlytconc, color = 'red')

########################################################################

res_mean = numpy.mean(AConc.data, axis=0)*1e6

plotres(100, AConc.time[0], res_mean, tetrads)
