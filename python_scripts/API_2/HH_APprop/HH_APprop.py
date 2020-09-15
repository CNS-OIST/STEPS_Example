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

meshPath = 'meshes/axon_cube_L1000um_D443nm_equiv0.5_19087tets.inp'

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
    n0, n1, n2, n3, n4 = SubUnitState.Create()
    nKSU = SubUnit.Create([n0, n1, n2, n3, n4])
    VGKC = Channel.Create([nKSU])

    # Sodium channel
    h0, h1, m0, m1, m2, m3 = SubUnitState.Create()
    mNaSU, hNaSU = SubUnit.Create(
        [m0, m1, m2, m3],
        [h0, h1],
    )
    VGNaC = Channel.Create([mNaSU, hNaSU])

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
            n0.s <r[1]> n1.s <r[2]> n2.s <r[3]> n3.s <r[4]> n4.s
            r[1].setRates(4*_a_n,   _b_n)
            r[2].setRates(3*_a_n, 2*_b_n)
            r[3].setRates(2*_a_n, 3*_b_n)
            r[4].setRates(  _a_n, 4*_b_n)

        with VGNaC[...]:
            h0.s <r[1]> h1.s
            r[1].setRates(_a_h, _b_h)

            m0.s <r[1]> m1.s <r[2]> m2.s <r[3]> m3.s
            r[1].setRates(3*_a_m,   _b_m)
            r[2].setRates(2*_a_m, 2*_b_m)
            r[3].setRates(  _a_m, 3*_b_m)

        # Create ohmic current objects
        VGKC_I = OhmicCurr.Create(VGKC[n4], K_G, K_rev)
        VGNaC_I = OhmicCurr.Create(VGNaC[m3, h1], Na_G, Na_rev)
        Leak_I = OhmicCurr.Create(Leak[lsus], L_G, leak_rev)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # TETRAHEDRAL MESH  # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

mesh = TetMesh.LoadAbaqus(meshPath, scale=1e-6)

# # # # # # # # # # # # # # # MESH MANIPULATION # # # # # # # # # # # # # # # # #

facetris = TriList(tri for tri in mesh.tris if tri.center.z == mesh.bbox.min.z)
injverts = facetris.verts

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

sim = Simulation('Tetexact', model, mesh, rng, True)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rs = ResultSelector(sim)

NaCurrs = rs.TRIS(memb_tris).VGNaC_I.I
KCurrs = rs.TRIS(memb_tris).VGKC_I.I
CellPot = rs.TETS(pot_tet).V

NaCurrs.metaData['trizpos'] = [tri.center.z for tri in memb_tris]
KCurrs.metaData['trizpos'] = [tri.center.z for tri in memb_tris]
NaCurrs.metaData['triarea'] = [tri.Area for tri in memb_tris]
KCurrs.metaData['triarea'] = [tri.Area for tri in memb_tris]
CellPot.metaData['tetzpos'] = pot_pos

NaCurrs.toFile('NaCurrs.dat')
KCurrs.toFile('KCurrs.dat')
CellPot.toFile('CellPot.dat')

sim.toSave(NaCurrs, KCurrs, CellPot, dt=DT_sim)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

sim.newRun()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Inject channels
surfarea = sim.patch.Area

for i, m in enumerate(mNaSU):
    for j, h in enumerate(hNaSU):
        sim.patch.VGNaC[m, h].Count = Na_ro * surfarea * Na_facs[j][i]

for i, n in enumerate(nKSU):
    sim.patch.VGKC[n].Count = K_ro * surfarea * K_facs[i]

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

