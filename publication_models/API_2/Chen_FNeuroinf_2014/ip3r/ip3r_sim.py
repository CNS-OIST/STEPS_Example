import steps.interface

# IP3 receptor mesh simulation

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *

###############################################################################

def getModel():
	mdl = Model()
	# WARNING: Using a variable name that is reserved (['r']).
	r = ReactionManager()
	with mdl:
		
		# chemical species objects
		# Calcium
		# IP3
		
		# receptor state objects
		# IP3 receptor in 'naive' state
		# bound IP3 
		# bound IP3 and Ca (open)
		# 1 bound Ca to inactivation site
		# 2 bound Ca to inactivation sites
		# 3 bound Ca to inactivation sites
		
		# chemical species objects
		Ca, IP3, R, RIP3, Ropen, RCa, R2Ca, R3Ca, R4Ca = Species.Create()				# Calcium
		
		surfsys = SurfaceSystem.Create()
	with surfsys, mdl:
		
		# The 'forward' binding reactions: 
			
		# The 'backward' binding reactions:
		R.s + IP3.o <r['R_bind_IP3_f']> RIP3.s
		RIP3.s + Ca.o <r['RIP3_bind_Ca_f']> Ropen.s
		R.s + Ca.o <r['R_bind_Ca_f']> RCa.s
		RCa.s + Ca.o <r['RCa_bind_Ca_f']> R2Ca.s
		R2Ca.s + Ca.o <r['R2Ca_bind_Ca_f']> R3Ca.s
		R3Ca.s + Ca.o <r['R3Ca_bind_Ca_f']> R4Ca.s
		
		# Ca ions passing through open IP3R channel
		Ca.i + Ropen.s <r['R_Ca_channel_f']> Ropen.s + Ca.o
	
	# The reaction constants
	r['R_bind_IP3_f'].setRates(1000000000.0, 25800.0)
	r['RIP3_bind_Ca_f'].setRates(8000000000.0, 2000.0)
	r['R_bind_Ca_f'].setRates(8889000.0, 5.0)
	r['RCa_bind_Ca_f'].setRates(20000000.0, 10.0)
	r['R2Ca_bind_Ca_f'].setRates(40000000.0, 15.0)
	r['R3Ca_bind_Ca_f'].setRates(60000000.0, 20.0)
	
	# Corresponds to Ca input ~ 20000/ms for open receptor
	r['R_Ca_channel_f'].setRates(8000000.0, 8000000.0)          
	
	return mdl

###############################################################################

def getGeom():
	mesh  = sgeom.Geom()
	cyt = sgeom.Comp('cyt', mesh)
	cyt.vol = 0.1e-18
	
	ER = sgeom.Comp('ER', mesh)					
	ER.setVol(0.02e-18)
	
	ERmemb = sgeom.Patch('ERmemb', mesh, ER, cyt)
	ERmemb.addSurfsys('surfsys')
	
	return mesh

###############################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# DIFFUSION

# Source:
#   Allbritton, N.L., Meyer, T., and Stryer, L. (1992). 
#   Range of messenger action of calcium ion and inositol 
#   1,4,5-triphosphate. Science 258, 1812-1815.
DCST_Ca = 0.065e-9
DCST_IP3 = 0.283e-9

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# import ip3r_model 

# Import model
mdl = getModel()
with mdl:
	
	volsys = VolumeSystem.Create()

# Fetch reference to Calcium and IP3 Spec objects
Ca = mdl.getSpec('Ca')
IP3 = mdl.getSpec('IP3')
with volsys, mdl:
	
	# Create diffusion rules
	Ca_diff = 	Diffusion.Create(Ca, DCST_Ca)

	IP3_diff = 	Diffusion.Create(IP3, DCST_IP3)


mesh = TetMesh.Load('ip3r_mesh')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Create random number generator
# WARNING: Using a variable name that is reserved (['r']).
r = RNG('mt19937', 512, 456)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Create reaction-diffusion solver object
# WARNING: Using a variable name that is reserved (['r']).
sim = Simulation('Tetexact', mdl, mesh, r)

# Setup initial condition
sim.cyt.Ca.Count = 1
sim.cyt.IP3.Conc = 2.5e-06
sim.ER.Ca.Conc = 0.00015
sim.ERmemb.R.Count = 16

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

