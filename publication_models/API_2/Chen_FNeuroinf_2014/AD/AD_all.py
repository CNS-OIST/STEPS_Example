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


meshes = []
rngs = []
for d in [0, 2, 4, 8]:
    meshes.append(TetMesh.LoadAbaqus(f'{d}_20_0.7.inp', scale=1e-6, ebs=None, shadow_mesh=f'{d}_20_0.7_conf'))
    rngs.append(RNG('mt19937', 512, int(time.time())))

sims = []
for mesh, rng in zip(meshes, rngs):
    sim = Simulation('Tetexact', mdl, mesh, rng)
    sim.injection.X.Count = 2000
    sims.append(sim)

########################################################################
# Visualization

# Create control
sc = SimControl(end_time = 1.0, upd_interval = 0.00001)
with sc:
    display_names = ["Density %i" % (d) for d in [0,2,4,8]]

    distPlots = PlotDisplay('Plots')

    for sim, dname in zip(sims, display_names):
        rs = ResultSelector(sim)

        with SimDisplay(dname):
            ElementDisplay(rs.dend, color=[0, 0, 1, 0.2])
            ElementDisplay(rs.shaft.X, color=[1, 0, 0, 1], spec_size = 0.1)

        with distPlots:
            mesh = sim.geom
            SpatialPlot(rs.TETS(mesh.shaft.tets).X.Count, title=f"{dname}_shaft", x_label=('Z Coordinate', 'm'), axis=[0, 0, 1], nbins=40)
            SpatialPlot(rs.TETS(mesh.dend.tets).X.Count, title=f"{dname}_dend", x_label=('Z Coordinate', 'm'), axis=[0, 0, 1], nbins=40)
            NewRow()
        
# Enter visualization loop
sc.run()
