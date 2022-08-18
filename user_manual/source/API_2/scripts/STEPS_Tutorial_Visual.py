# Example: Simulation visualization
# http://steps.sourceforge.net/manual/API_2/STEPS_Tutorial_Visual.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
from steps.visual import *

# Source:
#   Allbritton, N.L., Meyer, T., and Stryer, L. (1992).
#   Range of messenger action of calcium ion and inositol
#   1,4,5-triphosphate. Science 258, 1812-1815.
DCST_Ca = 0.065e-9
DCST_IP3 = 0.283e-9

#########################
# Model, geom, Simulation
#########################

mdl = Model()
r = ReactionManager()

with mdl:
    # chemical species objects
    Ca, IP3, R, RIP3, Ropen, RCa, R2Ca, R3Ca, R4Ca = Species.Create()

    ssys = SurfaceSystem.Create()
    with ssys:
        # IP3 and activating Ca binding
        R.s    + IP3.o <r[1]> RIP3.s
        RIP3.s + Ca.o  <r[2]> Ropen.s
        r[1].K = 1000000000.0, 25800.0
        r[2].K = 8000000000.0, 2000.0

        # Inactivating Ca binding
        R.s    + Ca.o <r[3]> RCa.s
        RCa.s  + Ca.o <r[4]> R2Ca.s
        R2Ca.s + Ca.o <r[5]> R3Ca.s
        R3Ca.s + Ca.o <r[6]> R4Ca.s
        r[3].K = 8889000.0, 5.0
        r[4].K = 20000000.0, 10.0
        r[5].K = 40000000.0, 15.0
        r[6].K = 60000000.0, 20.0

        # Ca ions passing through open IP3R channel
        Ca.i + Ropen.s <r[1]> Ropen.s + Ca.o
        # Corresponds to Ca input ~ 20000/ms for open receptor
        r[1].K = 8000000.0, 8000000.0

    vsys = VolumeSystem.Create()
    with vsys:
        # Create diffusion rules
        Ca_diff =  Diffusion.Create(Ca,  DCST_Ca)
        IP3_diff = Diffusion.Create(IP3, DCST_IP3)

mesh = TetMesh.Load('../meshes/ip3r_mesh')

# Create random number generator
rng = RNG('mt19937', 512, 456)

# Create reaction-diffusion solver object
sim = Simulation('Tetexact', mdl, mesh, rng)

sim.newRun()

# Setup initial condition
sim.cyt.Ca.Count = 1
sim.cyt.IP3.Conc = 2.5e-06
sim.ER.Ca.Conc = 0.00015
sim.memb.R.Count = 16

rs = ResultSelector(sim)

#########################
# Visualization
#########################

# Create control
sc = SimControl(end_time = 1.0, upd_interval = 0.0001)

with sc:
    # Plots
    plots_d = PlotDisplay('IP3 Receptor Model')
    with plots_d:

        TimePlot(rs.cyt.Ca.Conc,
            title='Ca_cyt',
            pen=(255, 0.647 * 255, 0),
            data_size=1000,
            y_range=[0, 15e-6],
            y_label=('Concentration', 'M')
        )

        TimePlot(rs.memb.Ropen.Count,
            title='Ropen_memb',
            pen=(255, 0, 255),
            data_size=1000,
            y_range=[0, 10],
            y_label=('Count', '#')
        )

        NewRow()

        SpatialPlot(rs.TETS(mesh.cyt.tets).Ca.Count,
            title='Ca_cyt distribution anlong y-axis',
            axis=[0, 1, 0],
            y_range=[0, 30]
        )

        TimePlot(rs.memb.MATCH('R.*').Count,
            title='Receptor states',
            data_size=1000,
            y_range=[0, 17],
            y_label=('Count', '#')
        )

    # 3D displays
    ER_d, CytIP3_d, CytCa_d, memb_d, full_d = SimDisplay.Create(
        'ER', 'Cyt IP3', 'Cyt Calcium', 'memb', 'Full view'
    )

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
        ElementDisplay(
            rs.memb.MATCH('R.*'), spec_size = 0.01,
            color=lambda spec: spec2Col[spec.name]
        )

    # Merge all sub-displays to the full one
    full_d.merge(ER_d, CytIP3_d, CytCa_d, memb_d)

# Enter visualization loop
sc.run()
