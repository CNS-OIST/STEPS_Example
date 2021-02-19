# IP3 receptor mesh simulation

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
from steps.visual import *

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
with mdl:

    vsys = VolumeSystem.Create()
    with vsys:

        # Create diffusion rules (fetch reference to Ca and IP3 from mdl)
        Ca_diff =  Diffusion.Create(mdl.Ca,  DCST_Ca)
        IP3_diff = Diffusion.Create(mdl.IP3, DCST_IP3)

mesh = TetMesh.Load('ip3r_mesh')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Create random number generator
rng = RNG('mt19937', 512, 456)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Create reaction-diffusion solver object
sim = Simulation('Tetexact', mdl, mesh, rng)

# Setup initial condition
sim.cyt.Ca.Count = 1
sim.cyt.IP3.Conc = 2.5e-06
sim.ER.Ca.Conc = 0.00015
sim.memb.R.Count = 16

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Visualization

rs = ResultSelector(sim)

# Create control
sc = SimControl(end_time = 1.0, upd_interval = 0.0001)

with sc:
    # Plots
    with PlotDisplay('IP3 Receptor Model'):
        TimePlot(rs.cyt.Ca.Conc, title='Ca_cyt', pen=(255, 0.647 * 255, 0), data_size=1000, y_range=[0, 1e-5], y_label=('Concentration', 'M'))
        NewRow()
        TimePlot(rs.memb.Ropen.Count, title='Ropen_memb', pen=(255, 0, 255), data_size=1000, y_range=[0, 10], y_label=('Concentration', 'M'))

    # 3D displays
    ER_d, CytIP3_d, CytCa_d, memb_d, full_d = SimDisplay.Create('ER', 'Cyt IP3', 'Cyt Calcium', 'memb', 'Full view')

    with ER_d:
        # Static mesh element
        ElementDisplay(rs.ER, color=[0.678, 1, 0.184, 0.05])
        # Dynamic element
        ElementDisplay(rs.ER.Ca, color=[1, 0.647, 0, 1], spec_size = 0.005)

    with CytIP3_d:
        ElementDisplay(rs.cyt, color=[0.941, 1, 0.941, 0.05])
        ElementDisplay(rs.cyt.IP3, color=[1, 0, 0, 1], spec_size = 0.005)

    with CytCa_d:
        ElementDisplay(rs.cyt, color=[0.941, 1, 0.941, 0.05])
        ElementDisplay(rs.cyt.Ca, color=[1, 0.647, 0, 1], spec_size = 0.005)

    with memb_d:
        ElementDisplay(rs.memb, color=[1, 0.973, 0.863, 0.05])

        # Different colors depending on receptor state
        spec2Col = {"R" : [0.0, 0.0, 1.0, 1.0], "RIP3" : [1.0, 0.0, 1.0, 0.2], 
            "Ropen" : [1.0, 0.0, 1.0, 1.0], "RCa" : [0.0, 0.0, 1.0, 0.8], 
            "R2Ca" : [0.0, 0.0, 1.0, 0.6], "R3Ca" : [0.0, 0.0, 1.0, 0.4], 
            "R4Ca" : [0.0, 0.0, 1.0, 0.2]
        }
        ElementDisplay(rs.memb.MATCH('R.*'), spec_size = 0.01, color=lambda spec: spec2Col[spec.name])

    # Merge all sub-displays to the full one
    full_d.merge(ER_d, CytIP3_d, CytCa_d, memb_d)
        
# Enter visualization loop
sc.run()
