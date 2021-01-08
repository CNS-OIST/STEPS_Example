#################################################################################
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#  Example: Hodgkin-Huxley Action Potential propagation model
#  Author Iain Hepburn
#  http://steps.sourceforge.net/manual/memb_pot.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # IMPORTS # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.saving import *
from steps.rng import *

import numpy as np
import os
import math

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# # # # # # # # # # # # # # # # # # PARAMETERS  # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # CHANNELS  # # # # # # # # # # # # # # # # # #

# Potassium conductance = 0.036 S/cm2
# Sodium conductance = 0.120 S/cm2

# Potassium single-channel conductance
K_G = 20.0e-12 # Siemens

# Potassium channel density
K_ro = 18.0e12 # per square meter

# Potassium reversal potential
K_rev = -77e-3 # volts

# Sodium single-channel conductance
Na_G = 20.0e-12 # Siemens

# Sodium channel density
Na_ro = 60.0e12 # per square meter

# Sodium reversal potential
Na_rev = 50e-3 # volts

# Leak single-channel conductance
L_G = 0.3e-12 # Siemens

# Leak density
L_ro = 10.0e12 # per square meter

# Leak reveral potential
leak_rev = -54.4e-3 # volts


# A table of potassium channel population factors: 
# n0, n1, n2, n3, n4
K_facs = [ 0.21768, 0.40513, 0.28093, 0.08647, 0.00979 ]

# A table of sodium channel population factors
# m0h0, m1h0, m2h0, m3h0, m0h1, m1h1, m2h1, m3h1:
Na_facs = [[0.34412, 0.05733, 0.00327, 6.0e-05],
           [0.50558, 0.08504, 0.00449, 0.00010]]

# # # # # # # # # # # # # # # # # # MESH  # # # # # # # # # # # # # # # # # # # # 

dirPath = os.path.dirname(os.path.abspath(__file__))
meshPath = os.path.join(dirPath, 'meshes/axon_cube_L1000um_D443nm_equiv0.5_19087tets.inp')

# # # # # # # # # # # # # # # SIMULATION CONTROLS # # # # # # # # # # # # # # # #

# Temperature for gating kinetics
celsius = 20.0

# Current injection
Iclamp = 50.0e-12 #	amps

# Voltage range for gating kinetics in Volts
Vrange = [-100.0e-3, 50e-3, 1e-4]

# The simulation dt
DT_sim = 1.0e-4 # seconds

# The time until which the simulation should be run
ENDT = 4.0e-3

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # BIOCHEMICAL MODEL # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

model = Model()

r = ReactionManager()

with model:
    ssys = SurfaceSystem.Create()

    #  Potassium channel
    Ko, Kc = SubUnitState.Create()
    KSU = SubUnit.Create([Ko, Kc])
    VGKC = Channel.Create([KSU]*4)

    # Sodium channel
    Na_mo, Na_mc, Na_hi, Na_ha = SubUnitState.Create()
    NamSU, NahSU = SubUnit.Create(
        [Na_mo, Na_mc],
        [Na_hi, Na_ha]
    )
    VGNaC = Channel.Create([NamSU, NamSU, NamSU, NahSU])

    # Leak channel
    lsus = SubUnitState.Create()
    Leak = Channel.Create([lsus])

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    # Hodgkin-Huxley gating kinetics

    # Temperature dependence
    thi = math.pow(3.0, ((celsius-6.3)/10.0)) * 1.0e3

    _a_n = VDepRate(lambda V: thi*((0.01*(10-(V*1e3+65))/(math.exp((10-(V*1e3+65))/10)-1))), vrange=Vrange)
    _b_n = VDepRate(lambda V: thi*((0.125*math.exp(-(V*1e3+65)/80))), vrange=Vrange)
    
    _a_m = VDepRate(lambda V: thi*((0.1*(25-(V*1e3+65))/(math.exp((25-(V*1e3+65))/10)-1))), vrange=Vrange)
    _b_m = VDepRate(lambda V: thi*((4*math.exp(-(V*1e3+65)/18))), vrange=Vrange)

    _a_h = VDepRate(lambda V: thi*((0.07*math.exp(-(V*1e3+65)/20))), vrange=Vrange)
    _b_h = VDepRate(lambda V: thi*((1/(math.exp((30-(V*1e3+65))/10)+1))), vrange=Vrange)

    with ssys:

        with VGKC[...]:
            Kc.s <r[1]> Ko.s
            r[1].setRates(_a_n, _b_n)

        with VGNaC[...]:
            Na_hi.s <r[1]> Na_ha.s
            r[1].setRates(_a_h, _b_h)
            
            Na_mc.s <r[1]> Na_mo.s
            r[1].setRates(_a_m, _b_m)

        # Create ohmic current objects
        VGKC_I = OhmicCurr.Create(VGKC[Ko, Ko, Ko, Ko], K_G, K_rev)
        VGNaC_I = OhmicCurr.Create(VGNaC[Na_mo, Na_mo, Na_mo, Na_ha], Na_G, Na_rev)
        Leak_I = OhmicCurr.Create(Leak[lsus], L_G, leak_rev)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # TETRAHEDRAL MESH  # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

mesh = TetMesh.LoadAbaqus(meshPath, scale=1e-6)

# # # # # # # # # # # # # # # MESH MANIPULATION # # # # # # # # # # # # # # # # #

facetris = TriList(tri for tri in mesh.tris if tri.center.z == mesh.bbox.min.z)
injverts = facetris.verts

if MPI.rank == 0:
    print("Found ", len(injverts), "I_inject vertices")
    print("Found ", len(facetris), "triangles on bottom face")

memb_tris = mesh.surface - facetris

# The points along (z) axis at which to record potential
pot_pos = np.arange(mesh.bbox.min.z, mesh.bbox.max.z, 10e-6)
pot_tet = [mesh.tets[(0, 0, z)] for z in pot_pos]

# # # # # # # # # # # # # # # GEOMETRY OBJECTS  # # # # # # # # # # # # # # # # #

with mesh:
    # Create cytosol compartment
    cyto = TetComp.Create(mesh.tets)

    # Create the patch and associate with surface system ssys
    patch = TetPatch.Create(memb_tris, cyto, None, ssys)

    # Create the membrane across which the potential will be solved
    membrane = Membrane.Create([patch], opt_method = 1)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rng = RNG('mt19937', 512, 1234)

partition = LinearMeshPartition(mesh, 1, 1, MPI.nhosts)

sim = Simulation('TetOpSplit', model, mesh, rng, True, partition)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rs = ResultSelector(sim)

CellPot = rs.TETS(pot_tet).V

CellPot.metaData['tetzpos'] = pot_pos

CellPot.toFile('CellPot_only.dat')

sim.toSave(CellPot, dt=DT_sim)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

sim.newRun()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Inject channels
surfarea = sim.patch.Area

for state in VGNaC:
    prop = Na_facs[state.Count(Na_ha)][state.Count(Na_mo)]
    sim.patch.VGNaC[state].Count = Na_ro * surfarea * prop

for state in VGKC:
    prop = K_facs[state.Count(Ko)]
    sim.patch.VGKC[state].Count = K_ro * surfarea * prop

sim.patch.Leak[lsus].Count = L_ro * surfarea

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Set some simulation variables:

# Set dt for membrane potential calculation to 0.01ms
sim.EfieldDT = 1.0e-5

# Initialize potential to -65mV
sim.membrane.Potential = -65e-3

# Set capacitance of the membrane to 1 uF/cm^2 = 0.01 F/m^2
sim.membrane.Capac = 1.0e-2

# Set resistivity of the conduction volume to 100 ohm.cm = 1 ohm.meter
sim.membrane.VolRes = 1.0

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Set the current clamp

sim.VERTS(injverts).IClamp = Iclamp/len(injverts)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Run the simulation
sim.run(ENDT)

