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

# Example: Diffusion boundary
# http://steps.sourceforge.net/manual/diffusion_boundary.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.saving import *
from steps.utils import *
from steps.rng import *

import numpy
import os
import sys
from pylab import *

########################################################################

dirPath = os.path.dirname(os.path.abspath(__file__))
meshPath = os.path.join(dirPath, '../../meshes/diffusion_boundary/cyl_len10_diam1')

DT = 0.001

########################################################################

def gen_model():
    # Create the model container object
    mdl = Model()
    with mdl:
        # Create the chemical species
        X, Y = Species.Create()
        # Create separate volume systems for compartments A and B
        vsysA, vsysB = VolumeSystem.Create()
        # Describe diffusion of molecules in compartments A and B
        for vsys in [vsysA, vsysB]:
            with vsys:
                Diffusion(X, 0.1e-9)
                Diffusion(Y, 0.1e-9)

    # Return the container object
    return mdl

########################################################################

def gen_geom(mdl):
    mesh = TetMesh.Load(meshPath)

    with mesh:
        tetsA, tetsB = TetList(), TetList()
        z_min, z_max = mesh.bbox.min.z, mesh.bbox.max.z
        z_mid = z_min+(z_max-z_min)/2.0

        tetsA = TetList(tet for tet in mesh.tets if tet.center.z < z_mid)
        tetsB = mesh.tets - tetsA

        # Create the mesh compartments
        compA, compB = Compartment.Create(
            Params(tetsA, mdl.vsysA), 
            Params(tetsB, mdl.vsysB)
        )

        diffb = DiffBoundary.Create(tetsA.surface & tetsB.surface)

    return mesh

########################################################################

mdl = gen_model()
mesh = gen_geom(mdl)

rng = RNG('mt19937', 256, 432) 

sim = Simulation('Tetexact', mdl, mesh, rng)

rs = ResultSelector(sim)

XCount = rs.TETS().X.Count
YCount = rs.TETS().Y.Count

sim.toSave(XCount, YCount, dt=DT)

sim.newRun()

sim.TET(0, 0, -4.99e-6).X.Count = 1000
sim.TET(0, 0, 4.99e-6).Y.Count = 500

sim.diffb.Y.DiffusionActive = True

sim.run(0.1)

########################################################################

def plot_binned(mesh, resX, resY, bin_n = 50):
    zbound_min = mesh.bbox.min.z
    z_tets = np.array([(tet.center.z - zbound_min)*1e6 for tet in mesh.tets])
    vol_tets = np.array([tet.Vol * 1e18 for tet in mesh.tets])

    bins = numpy.histogram_bin_edges(z_tets, bin_n)
    tet_bins = numpy.digitize(z_tets, bins)

    # Ignore empty bins
    with numpy.errstate(invalid='ignore'):
        # Compute average position for each bin
        bins_pos = numpy.bincount(tet_bins, weights=z_tets) / numpy.bincount(tet_bins)
        # Compute total volume of each bin
        bin_vols = numpy.bincount(tet_bins, weights=vol_tets)

        Xbinned  = numpy.bincount(tet_bins, weights=resX) / bin_vols
        Ybinned  = numpy.bincount(tet_bins, weights=resY) / bin_vols

    scatter(bins_pos, Xbinned, label = 'X', color = 'blue')
    scatter(bins_pos, Ybinned, label = 'Y', color = 'red')

    xlabel('Z axis (microns)', fontsize=16)
    ylabel('Bin concentration (N/um^3)', fontsize=16)
    ylim(0)    
    xlim(0, 10)
    legend(numpoints=1)
    show()

########################################################################

plot_binned(mesh, XCount.data[0, -1], YCount.data[0, -1], 50)

