########################################################################
#                                                                      #
#                            Anomalous Diffusion                       #
#                                                                      #
########################################################################

import steps.interface

########################################################################
# Create Model

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
from steps.visual import *

import time

mdl = Model()
r = ReactionManager()

with mdl:
    X = Species.Create()
    vsys = VolumeSystem.Create()
    with vsys:
        dif_X = Diffusion.Create(X, 2e-09)

########################################################################
# Create Gemoetry

tetmesh = TetMesh.LoadAbaqus('2_20_0.7.inp', scale=1e-06, ebs=None, shadow_mesh="2_20_0.7_conf")

########################################################################
# Create Random number generator

rng = RNG('mt19937', 512, int(time.time()%4294967295))

########################################################################
# Initialize simulation

sim = Simulation('Tetexact', mdl, tetmesh, rng)

sim.injection.X.Count = 2000

########################################################################
# Visualization

rs = ResultSelector(sim)

# Create control
sc = SimControl(end_time = 1.0, upd_interval = 0.00001)

with sc:
    with SimDisplay('Show Spine Species'):
        # Static mesh element
        ElementDisplay(rs.dend, color=[0, 0, 1, 0.2])

        # Dynamic element
        ElementDisplay(rs.LIST('dend', 'shaft').X, color=[1.0, 0.0, 0.0, 1.0], spec_size = 0.1)

    with SimDisplay('Hide Spine Species'):
        ElementDisplay(rs.dend, color=[0, 0, 1, 0.2])
        ElementDisplay(rs.shaft.X, color=[1.0, 0.0, 0.0, 1.0], spec_size = 0.1)

    with PlotDisplay('Plots'):
        SpatialPlot(rs.TETS(tetmesh.shaft.tets).X.Count, axis=[0, 0, 1], nbins=100)
        
# Enter visualization loop
sc.run()

