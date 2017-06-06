# IP3 receptor mesh simulation

import steps.model as smodel
import steps.geom as swm
import steps.rng as srng
import steps.solver as ssolver

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# DIFFUSION

# Source:
#   Allbritton, N.L., Meyer, T., and Stryer, L. (1992). 
#   Range of messenger action of calcium ion and inositol 
#   1,4,5-triphosphate. Science 258, 1812-1815.
DCST_Ca = 0.065e-9
DCST_IP3 = 0.283e-9

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

import ip3r_model 

# Import model
mdl = ip3r_model.getModel()

volsys = smodel.Volsys('vsys', mdl)

# Fetch reference to Calcium and IP3 Spec objects
Ca = mdl.getSpec('Ca')
IP3 = mdl.getSpec('IP3')

# Create diffusion rules
Ca_diff = smodel.Diff('Ca_diff', volsys, Ca, DCST_Ca)
IP3_diff = smodel.Diff('IP3_diff', volsys, IP3, DCST_IP3)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Import mesh
import steps.utilities.meshio as meshio

mesh = meshio.loadMesh("ip3r_mesh")[0]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Create random number generator
r = srng.create('mt19937', 512)
r.initialize(456)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Create reaction-diffusion solver object
sim = ssolver.Tetexact(mdl, mesh, r)

# Setup initial condition
sim.setCompConc('cyt', 'Ca', 3.30657e-8)
sim.setCompConc('cyt', 'IP3', 2.5e-6)
sim.setCompConc('ER', 'Ca', 150e-6)
sim.setPatchCount('memb', 'R', 16)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Visualization
import pyqtgraph as pg
import steps.visual as visual

# Visualization initialization
app = pg.mkQApp()

# Create plot display
plots = visual.PlotDisplay("IP3 Receptor Model", size = (600, 400))

# Create Plots
pen = pg.mkPen(color=(255,255,255), width=2)
p = plots.addCompSpecPlot("<span style='font-size: 16pt'>Ca_cyt", sim, "cyt", "Ca", data_size = 1000,y_range= [0, 1e-5], measure = "conc", pen=(255, 0.647 * 255, 0))
p.getAxis('left').setPen(pen)
p.getAxis('bottom').setPen(pen)
p.showGrid(x=True, y=True)
labelStyle = {'color': '#ffffff', 'font-size': '16px'}
p.setLabel('bottom', 'Time', 's', **labelStyle)

plots.nextRow()

p = plots.addPatchSpecPlot("<span style='font-size: 16pt'>Ropen_memb", sim, "memb", "Ropen", data_size = 1000,y_range= [0, 10], pen=(255, 0, 255))
p.getAxis('left').setPen(pen)
p.getAxis('bottom').setPen(pen)
p.showGrid(x=True, y=True)
p.setLabel('bottom', 'Time', 's', **labelStyle)

# Create simulation displays
ER_display = visual.SimDisplay("ER", w = 600, h = 400)
cytIP3_display = visual.SimDisplay("Cyt IP3", w = 600, h = 400)
cytCa_display = visual.SimDisplay("Cyt Calcium", w = 600, h = 400)
memb_display = visual.SimDisplay("memb", w = 600, h = 400)
full_display = visual.SimDisplay("Full View", w = 600, h = 400)

# Create static mesh components
ER_view = visual.VisualCompMesh("ER", full_display, mesh, "ER", color = [0.678, 1.000, 0.184, 0.05])
cyt_view = visual.VisualCompMesh("cyt", full_display, mesh, "cyt", color = [0.941, 1.000, 0.941, 0.05])
memb_view = visual.VisualPatchMesh("memb", full_display, mesh, "memb", color = [1.000, 0.973, 0.863, 0.05])

# Create dynamic species components
Ca_ER = visual.VisualCompSpec("Ca_ER", full_display, mesh, sim, "ER", "Ca", [1.000, 0.647, 0.000, 1.0], spec_size = 0.005)
IP3_cyt = visual.VisualCompSpec("IP3_cyt", full_display, mesh, sim, "cyt", "IP3", [1.0, 0.0, 0.0, 1.0], spec_size = 0.005)
Ca_cyt = visual.VisualCompSpec("Ca_cyt", full_display, mesh, sim, "cyt", "Ca", [1.000, 0.647, 0.000, 1.0], spec_size = 0.005)
IP3R_MEMB = visual.VisualPatchChannel("IP3R_memb", full_display, mesh, sim, "memb", {"R" : [0.0, 0.0, 1.0, 1.0], "RIP3" : [1.0, 0.0, 1.0, 0.2], "Ropen" : [1.0, 0.0, 1.0, 1.0], "RCa" : [0.0, 0.0, 1.0, 0.8], "R2Ca" : [0.0, 0.0, 1.0, 0.6], "R3Ca" : [0.0, 0.0, 1.0, 0.4], "R4Ca" : [0.0, 0.0, 1.0, 0.2]}, spec_size = 0.01)

# Add associated components to individual displays
ER_display.addItem(ER_view)
ER_display.addItem(Ca_ER)

cytCa_display.addItem(cyt_view)
cytCa_display.addItem(Ca_cyt)

cytIP3_display.addItem(cyt_view)
cytIP3_display.addItem(IP3_cyt)

memb_display.addItem(memb_view)
memb_display.addItem(IP3R_MEMB)

# Add simulation and displays to control
x = visual.SimControl([sim], [ER_display, cytIP3_display, cytCa_display, memb_display, full_display],[plots], end_time= 1.0, upd_interval = 0.0001)

# Enter visualization loop
app.exec_()

