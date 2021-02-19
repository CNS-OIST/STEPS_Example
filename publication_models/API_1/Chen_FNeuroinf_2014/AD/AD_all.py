########################################################################
#                                                                      #
#                            Anomalous Diffusion                       #
#                                                                      #
########################################################################


########################################################################
# Create Model

import steps.model as smodel

mdl = smodel.Model()
X = smodel.Spec('X', mdl)
dend_vsys = smodel.Volsys('vsys', mdl)
dif_X = smodel.Diff('diffX', dend_vsys, X, 2e-9)

########################################################################
# Create Gemoetry

import steps.utilities.meshio as meshio

meshes = []
for d in [0, 2, 4, 8]:
    meshes.append(meshio.importAbaqus('%i_20_0.7.inp' % (d), 1e-6, None, "%i_20_0.7_conf" % (d))[0])

########################################################################
# Create Random number generator

import steps.rng as srng
import time

rngs = []
for d in [0, 2, 4, 8]:
    rng = srng.create('mt19937', 512)
    rng.initialize(int(time.time()))
    rngs.append(rng)

########################################################################
# Initialize simulation

import steps.solver as ssolver
sims = []
for d in range(4):
    sim = ssolver.Tetexact(mdl, meshes[d], rngs[d])
    sim.setROICount("injection", "X", 2000)
    sims.append(sim)

########################################################################
# Visualization

import pyqtgraph as pg
import steps.visual as visual

app = pg.mkQApp()

# Create plot displays
display_names = ["Density %i" % (d) for d in [0,2,4,8]]
plots = visual.PlotDisplay("Molecule Distribution", size = (400, 800))
displays = []

labelStyle = {'color': '#ffffff', 'font-size': '16px'}

for d in range(4):
    # Create simulation display
    display = visual.SimDisplay(display_names[d], w = 600, h = 200)

    # Create static mesh component
    comp_view = visual.VisualCompMesh("dend_mesh", display, meshes[d], "dend", color = [0.0, 0.0, 1.0, 0.2])
    
    # Create dynamic species component
    spec_view = visual.VisualROITetsSpec("X_shaft", display, meshes[d], sims[d], "shaft", "X", [1.0, 0.0, 0.0, 1.0], spec_size = 0.1)
    
    # Position adjustment
    display.rotateItems(90, 1,0,0)
    display.widget.orbit(-45, 0)
    
    # Add to display list
    displays.append(display)
    
    # Creat plots
    plot = plots.addROISpecDist("<span style='font-size: 16pt'>" + display_names[d] + "_shaft", meshes[d], sims[d], "shaft", "X", axis = "z", y_range = (0, 500), nbins = 40,brush=(50,50,200,100))
    plot.showGrid(x=True, y=True)
    plot.setLabel('bottom', 'Z Coordinate', 'm', **labelStyle)
    plot = plots.addCompSpecDist("<span style='font-size: 16pt'>" + display_names[d] + "_dend", meshes[d], sims[d], "dend", "X", axis = "z", y_range = (0, 500), nbins = 40,brush=(50,50,200,100))
    plot.showGrid(x=True, y=True)
    plot.setLabel('bottom', 'Z Coordinate', 'm', **labelStyle)
    if d != 3:
        plots.nextRow()

# Create control
x= visual.SimControl(sims, displays, [plots], end_time = 1.0, upd_interval = 0.00001)

# Enter visualization loop
app.exec_()
