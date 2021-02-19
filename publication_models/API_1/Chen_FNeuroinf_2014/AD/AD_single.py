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

tetmesh = meshio.importAbaqus('2_20_0.7.inp', 1e-6, None, "2_20_0.7_conf")[0]

########################################################################
# Create Random number generator

import steps.rng as srng
import time

rng = srng.create('mt19937', 512)
rng.initialize(int(time.time()%4294967295)) # The max unsigned long

########################################################################
# Initialize simulation

import steps.solver as ssolver
sim = ssolver.Tetexact(mdl, tetmesh, rng)

sim.setROICount("injection", "X", 2000)

########################################################################
# Visualization

import pyqtgraph as pg
import steps.visual as visual

app = pg.mkQApp()

# Create simulation displays
display1 = visual.SimDisplay("Show Spine Species")
display2 = visual.SimDisplay("Hide Spine Species")

# Create plot display
plots = visual.PlotDisplay("Plots", size = (1024, 800))

# Create plot
plots.addROISpecDist("Shaft_dist", tetmesh, sim, "shaft", "X", axis = "z", y_range = (0, 300), nbins = 100,brush=(50,50,200,100))

# Create static mesh component
comp_view = visual.VisualCompMesh("dend_mesh", display1, tetmesh, "dend", color = [0.0, 0.0, 1.0, 0.2])
display2.addItem(comp_view)
comp_view.rotate(90, 1,0,0)

# Create dynamic species components
comp_spec_view = visual.VisualCompSpec("X_dend", display1, tetmesh, sim, "dend", "X", [1.0, 0.0, 0.0, 1.0], spec_size = 0.1)
comp_spec_view.rotate(90, 1,0,0)
roi_spec_view = visual.VisualROITetsSpec("X_shaft", display2, tetmesh, sim, "shaft", "X", [1.0, 0.0, 0.0, 1.0], spec_size = 0.1)
roi_spec_view.rotate(90, 1,0,0)

# Visualization position adjustment
display1.widget.orbit(-45, 0)
display2.widget.orbit(-45, 0)

# Create control
x= visual.SimControl([sim], [display1, display2], [plots], end_time = 1.0, upd_interval = 0.00001)

# Enter visualization loop
app.exec_()
