# Example: Simulating membrane potential
# http://steps.sourceforge.net/manual/API_2/STEPS_Tutorial_Efield.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.saving import *
from steps.rng import *

import numpy as np
import math

# Potassium conductance = 0.036 S/cm2

# Potassium single-channel conductance
K_G = 20.0e-12 # Siemens

# Potassium channel density
K_ro = 18.0e12 # per square meter

# Potassium reversal potential
K_rev = -77e-3 # volts

# Sodium conductance = 0.120 S/cm2

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

# Temperature for gating kinetics
celsius = 20.0

# Current injection
Iclamp = 50.0e-12 # amps

# Voltage range for gating kinetics in Volts
Vrange = [-100.0e-3, 50e-3, 1e-4]

def HHRateFunction(A, B, C, D, F, H, V):
    num = A + B * V
    denom = C + H * math.exp((V + D) / F)
    if num == denom == 0:
        return F * B / (H * math.exp((V + D) / F))
    else:
        return num / denom

# The simulation dt
DT_sim = 1.0e-4 # seconds

# The time until which the simulation should be run
ENDT = 4.0e-3

#########################
# Model setup
#########################

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

    thi = math.pow(3.0, ((celsius-6.3)/10.0))

    _a_n = VDepRate(
        lambda V: thi * 1e3 * HHRateFunction(-0.55, -0.01, -1, 55, -10, 1, V*1e3),
        vrange=Vrange
    )
    _b_n = VDepRate(
        lambda V: thi * 1e3 * HHRateFunction(1, 0, 0, 65, 80, 8, V*1e3),
        vrange=Vrange
    )

    _a_m = VDepRate(
        lambda V: thi * 1e3 * HHRateFunction(-4, -0.1, -1, 40, -10, 1, V*1e3),
        vrange=Vrange
    )
    _b_m = VDepRate(
        lambda V: thi * 1e3 * HHRateFunction(1, 0, 0, 65, 18, 0.25, V*1e3),
        vrange=Vrange
    )

    _a_h = VDepRate(
        lambda V: thi * 1e3 * HHRateFunction(1, 0, 0, 65, 20, 1 / 0.07, V*1e3),
        vrange=Vrange
    )
    _b_h = VDepRate(
        lambda V: thi * 1e3 * HHRateFunction(1, 0, 1, 35, -10, 1, V*1e3),
        vrange=Vrange
    )

    with ssys:
        with VGKC[...]:
            Kc.s <r[1]> Ko.s
            r[1].K = _a_n, _b_n

        with VGNaC[...]:
            Na_hi.s <r[1]> Na_ha.s
            r[1].K = _a_h, _b_h

            Na_mc.s <r[1]> Na_mo.s
            r[1].K = _a_m, _b_m
        
        VGKC_I = OhmicCurr.Create(VGKC[Ko, Ko, Ko, Ko], K_G, K_rev)
        VGNaC_I = OhmicCurr.Create(VGNaC[Na_mo, Na_mo, Na_mo, Na_ha], Na_G, Na_rev)
        Leak_I = OhmicCurr.Create(Leak[lsus], L_G, leak_rev)

#########################
# Geom setup
#########################

mesh = TetMesh.LoadAbaqus('../meshes/axon.inp', scale=1e-6)

with mesh:
    facetris = TriList([tri for tri in mesh.tris if tri.center.z == mesh.bbox.min.z])
    injverts = facetris.verts

    memb_tris = mesh.surface - facetris

    # The points along (z) axis at which to record potential
    pot_pos = np.arange(mesh.bbox.min.z, mesh.bbox.max.z, 10e-6)
    pot_tet = TetList(mesh.tets[0, 0, z] for z in pot_pos)

    cyto = Compartment.Create(mesh.tets)
    patch = Patch.Create(memb_tris, cyto, None, ssys)

    # Create the membrane across which the potential will be solved
    membrane = Membrane.Create([patch])

#########################
# Simulation setup
#########################

rng = RNG('mt19937', 512, 1234)

sim = Simulation('Tetexact', model, mesh, rng, True)

rs = ResultSelector(sim)

NaCurrs = rs.TRIS(memb_tris).VGNaC_I.I
KCurrs = rs.TRIS(memb_tris).VGKC_I.I
CellPot = rs.TETS(pot_tet).V

NaCurrs.metaData['trizpos'] = [tri.center.z for tri in memb_tris]
KCurrs.metaData['trizpos'] = [tri.center.z for tri in memb_tris]

NaCurrs.metaData['triarea'] = [tri.Area for tri in memb_tris]
KCurrs.metaData['triarea'] = [tri.Area for tri in memb_tris]

CellPot.metaData['tetzpos'] = pot_pos

sim.toSave(NaCurrs, KCurrs, CellPot, dt=DT_sim)

#########################
# Run simulation
#########################

with HDF5Handler('Efield') as hdf:
    sim.toDB(hdf, 'TetexactSim')

    sim.newRun()
    
    # Inject channels
    surfarea = sim.patch.Area
    
    for state in VGNaC:
        prop = Na_facs[state.Count(Na_ha)][state.Count(Na_mo)]
        sim.patch.VGNaC[state].Count = Na_ro * surfarea * prop
    
    for state in VGKC:
        prop = K_facs[state.Count(Ko)]
        sim.patch.VGKC[state].Count = K_ro * surfarea * prop
    
    sim.patch.Leak[lsus].Count = L_ro * surfarea
    
    # Set dt for membrane potential calculation to 0.01ms
    sim.EfieldDT = 1.0e-5
    
    # Initialize potential to -65mV
    sim.membrane.Potential = -65e-3
    
    # Set capacitance of the membrane to 1 uF/cm^2 = 0.01 F/m^2
    sim.membrane.Capac = 1.0e-2
    
    # Set resistivity of the conduction volume to 100 ohm.cm = 1 ohm.meter
    sim.membrane.VolRes = 1.0
    
    # Set the current clamp
    sim.VERTS(injverts).IClamp = Iclamp/len(injverts)
    
    # Run the simulation
    sim.run(ENDT)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#########################
# TetODE Simulation
#########################

sim = Simulation('TetODE', model, mesh, calcMembPot=True)

rs = ResultSelector(sim)

CellPot = rs.TETS(pot_tet).V

CellPot.metaData['tetzpos'] = pot_pos

sim.toSave(CellPot)

with HDF5Handler('Efield') as hdf:
    sim.toDB(hdf, 'TetODESim')

    sim.newRun()
    
    # Inject channels
    surfarea = sim.patch.Area
    
    for state in VGNaC:
        prop = Na_facs[state.Count(Na_ha)][state.Count(Na_mo)]
        sim.patch.VGNaC[state].Count = Na_ro * surfarea * prop
    
    for state in VGKC:
        prop = K_facs[state.Count(Ko)]
        sim.patch.VGKC[state].Count = K_ro * surfarea * prop
    
    sim.patch.Leak[lsus].Count = L_ro * surfarea
    
    # Initialize potential to -65mV
    sim.membrane.Potential = -65e-3
    
    # Set capacitance of the membrane to 1 uF/cm^2 = 0.01 F/m^2
    sim.membrane.Capac = 1.0e-2
    
    # Set resistivity of the conduction volume to 100 ohm.cm = 1 ohm.meter
    sim.membrane.VolRes = 1.0
    
    # Set the current clamp
    sim.VERTS(injverts).IClamp = Iclamp/len(injverts)
    
    sim.setTolerances(1e-3, 1e-4)
    
    # Run the simulation
    EFDt = 1e-5
    for i in range(int(ENDT // EFDt) + 1):
        sim.run(i * EFDt)
        if i % 10 == 0:
            CellPot.save()
