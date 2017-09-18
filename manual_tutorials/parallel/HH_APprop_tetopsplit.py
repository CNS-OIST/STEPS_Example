# This script is provided as example for STEPS user manual.
# License: GPL2.0
# Parallel modification for user manual
# Contact: Dr. Weiliang Chen, w.chen@oist.jp

# Original
# Example: Hodgkin-Huxley Action Potential propagation model
# Author Iain Hepburn
# http://steps.sourceforge.net/manual/memb_pot.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# # # # # # # # # # # # # # # # # # # IMPORTS # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from __future__ import print_function
import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.utilities.meshio as meshio

# MPI stuff
import steps.mpi
import steps.mpi.solver as mpi_solver
import steps.utilities.geom_decompose as gd


from pylab import *
import numpy
import math
import time
from random import *

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
Na_facs = [ 0.34412, 0.05733, 0.00327, 6.0e-05, \
                0.50558, 0.08504, 0.00449, 0.00010 ]

# # # # # # # # # # # # # # # # # # MESH  # # # # # # # # # # # # # # # # # # # # 

meshfile_ab = 'axon_cube_L1000um_D443nm_equiv0.5_19087tets.inp'

# # # # # # # # # # # # # # # SIMULATION CONTROLS # # # # # # # # # # # # # # # #

# Temperature for gating kinetics
celsius = 20.0		

# Current injection
Iclamp = 50.0e-12 #	amps

# Voltage range for gating kinetics in Volts
Vrange = [-100.0e-3, 50e-3, 1e-4]

# The number of simulation time-points
N_timepoints = 41

# The simulation dt
DT_sim = 1.0e-4 # seconds

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # BIOCHEMICAL MODEL # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

mdl = smodel.Model()
ssys = smodel.Surfsys('ssys', mdl)

# Potassium channel
K = smodel.Chan('K', mdl)
K_n0 = smodel.ChanState('K_n0', mdl, K)		
K_n1 = smodel.ChanState('K_n1', mdl, K)
K_n2 = smodel.ChanState('K_n2', mdl, K)
K_n3 = smodel.ChanState('K_n3', mdl, K)
K_n4 = smodel.ChanState('K_n4', mdl, K)

# Sodium channel
Na = smodel.Chan('Na', mdl)
Na_m0h0 = smodel.ChanState('Na_m0h0', mdl, Na)
Na_m1h0 = smodel.ChanState('Na_m1h0', mdl, Na)
Na_m2h0 = smodel.ChanState('Na_m2h0', mdl, Na)
Na_m3h0 = smodel.ChanState('Na_m3h0', mdl, Na)
Na_m0h1 = smodel.ChanState('Na_m0h1', mdl, Na)
Na_m1h1 = smodel.ChanState('Na_m1h1', mdl, Na)
Na_m2h1 = smodel.ChanState('Na_m2h1', mdl, Na)
Na_m3h1 = smodel.ChanState('Na_m3h1', mdl, Na)

# Leak channel
L = smodel.Chan('L', mdl)
Leak = smodel.ChanState('Leak', mdl, L)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Hodgkin-Huxley gating kinetics

# Temperature dependence
thi = math.pow(3.0, ((celsius-6.3)/10.0))

_a_n = lambda mV: thi*((0.01*(10-(mV+65.))/(math.exp((10-(mV+65.))/10.)-1)))

_b_n = lambda mV: thi*((0.125*math.exp(-(mV+65.)/80.)))

_a_m = lambda mV: thi*((0.1*(25-(mV+65.))/(math.exp((25-(mV+65.))/10.)-1)))

_b_m = lambda mV: thi*((4.*math.exp(-(mV+65.)/18.)))


_a_h = lambda mV: thi*((0.07*math.exp(-(mV+65.)/20.)))

_b_h = lambda mV: thi*((1./(math.exp((30-(mV+65.))/10.)+1)))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

Kn0n1 = smodel.VDepSReac('Kn0n1', ssys, slhs = [K_n0], srhs = [K_n1], \
                            k=lambda V: 1.0e3 *4.*_a_n(V*1.0e3), vrange = Vrange)
Kn1n2 = smodel.VDepSReac('Kn1n2', ssys, slhs = [K_n1], srhs = [K_n2], \
                            k=lambda V: 1.0e3 *3.*_a_n(V*1.0e3), vrange = Vrange)
Kn2n3 = smodel.VDepSReac('Kn2n3', ssys, slhs = [K_n2], srhs = [K_n3], \
                            k=lambda V: 1.0e3 *2.*_a_n(V*1.0e3), vrange = Vrange)
Kn3n4 = smodel.VDepSReac('Kn3n4', ssys, slhs = [K_n3], srhs = [K_n4], \
                            k=lambda V: 1.0e3 *1.*_a_n(V*1.0e3), vrange = Vrange)

Kn4n3 = smodel.VDepSReac('Kn4n3', ssys, slhs = [K_n4], srhs = [K_n3], \
                            k=lambda V: 1.0e3 *4.*_b_n(V*1.0e3), vrange = Vrange)
Kn3n2 = smodel.VDepSReac('Kn3n2', ssys, slhs = [K_n3], srhs = [K_n2], \
                            k=lambda V: 1.0e3 *3.*_b_n(V*1.0e3), vrange = Vrange)
Kn2n1 = smodel.VDepSReac('Kn2n1', ssys, slhs = [K_n2], srhs = [K_n1], \
                            k=lambda V: 1.0e3 *2.*_b_n(V*1.0e3), vrange = Vrange)
Kn1n0 = smodel.VDepSReac('Kn1n0', ssys, slhs = [K_n1], srhs = [K_n0], \
                            k=lambda V: 1.0e3 *1.*_b_n(V*1.0e3), vrange = Vrange)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

Na_m0h1_m1h1 = smodel.VDepSReac('Na_m0h1_m1h1', ssys, \
                                slhs=[Na_m0h1], srhs=[Na_m1h1], \
                                k=lambda V:1.0e3*3.*_a_m(V*1.0e3), vrange=Vrange)
Na_m1h1_m2h1 = smodel.VDepSReac('Na_m1h1_m2h1', ssys, \
                                slhs=[Na_m1h1], srhs=[Na_m2h1], \
                                k=lambda V:1.0e3*2.*_a_m(V*1.0e3), vrange=Vrange)
Na_m2h1_m3h1 = smodel.VDepSReac('Na_m2h1_m3h1', ssys, \
                                slhs=[Na_m2h1], srhs=[Na_m3h1], \
                                k=lambda V:1.0e3*1.*_a_m(V*1.0e3), vrange=Vrange)

Na_m3h1_m2h1 = smodel.VDepSReac('Na_m3h1_m2h1', ssys, \
                                slhs=[Na_m3h1], srhs=[Na_m2h1], \
                                k=lambda V:1.0e3*3.*_b_m(V*1.0e3), vrange=Vrange)
Na_m2h1_m1h1 = smodel.VDepSReac('Na_m2h1_m1h1', ssys, \
                                slhs=[Na_m2h1], srhs=[Na_m1h1], \
                                k=lambda V:1.0e3*2.*_b_m(V*1.0e3), vrange=Vrange)
Na_m1h1_m0h1 = smodel.VDepSReac('Na_m1h1_m0h1', ssys, \
                                slhs=[Na_m1h1], srhs=[Na_m0h1], \
                                k=lambda V:1.0e3*1.*_b_m(V*1.0e3), vrange=Vrange)

Na_m0h0_m1h0 = smodel.VDepSReac('Na_m0h0_m1h0', ssys, \
                                slhs=[Na_m0h0], srhs=[Na_m1h0], \
                                k=lambda V:1.0e3*3.*_a_m(V*1.0e3), vrange=Vrange)
Na_m1h0_m2h0 = smodel.VDepSReac('Na_m1h0_m2h0', ssys, \
                                slhs=[Na_m1h0], srhs=[Na_m2h0], \
                                k=lambda V:1.0e3*2.*_a_m(V*1.0e3), vrange=Vrange)
Na_m2h0_m3h0 = smodel.VDepSReac('Na_m2h0_m3h0', ssys, \
                                slhs=[Na_m2h0], srhs=[Na_m3h0], \
                                k=lambda V:1.0e3*1.*_a_m(V*1.0e3), vrange=Vrange)

Na_m3h0_m2h0 = smodel.VDepSReac('Na_m3h0_m2h0', ssys, \
                                slhs=[Na_m3h0], srhs=[Na_m2h0], \
                                k=lambda V:1.0e3*3.*_b_m(V*1.0e3), vrange=Vrange) 
Na_m2h0_m1h0 = smodel.VDepSReac('Na_m2h0_m1h0', ssys, \
                                slhs=[Na_m2h0], srhs=[Na_m1h0], \
                                k=lambda V:1.0e3*2.*_b_m(V*1.0e3), vrange=Vrange)
Na_m1h0_m0h0 = smodel.VDepSReac('Na_m1h0_m0h0', ssys, \
                                slhs=[Na_m1h0], srhs=[Na_m0h0], \
                                k=lambda V:1.0e3*1.*_b_m(V*1.0e3), vrange=Vrange)

Na_m0h0_m0h1 = smodel.VDepSReac('Na_m0h0_m0h1', ssys, \
                                slhs=[Na_m0h0], srhs=[Na_m0h1], \
                                k=lambda V:1.0e3*_a_h(V*1.0e3), vrange=Vrange)
Na_m1h0_m1h1 = smodel.VDepSReac('Na_m1h0_m1h1', ssys, \
                                slhs=[Na_m1h0], srhs=[Na_m1h1], \
                                k=lambda V:1.0e3*_a_h(V*1.0e3), vrange=Vrange)
Na_m2h0_m2h1 = smodel.VDepSReac('Na_m2h0_m2h1', ssys, \
                                slhs=[Na_m2h0], srhs=[Na_m2h1], \
                                k=lambda V:1.0e3*_a_h(V*1.0e3), vrange=Vrange)
Na_m3h0_m3h1 = smodel.VDepSReac('Na_m3h0_m3h1', ssys, \
                                slhs=[Na_m3h0], srhs=[Na_m3h1], \
                                k=lambda V:1.0e3*_a_h(V*1.0e3), vrange=Vrange)

Na_m0h1_m0h0 = smodel.VDepSReac('Na_m0h1_m0h0', ssys, \
                                slhs=[Na_m0h1], srhs=[Na_m0h0], \
                                k=lambda V:1.0e3*_b_h(V*1.0e3), vrange=Vrange)
Na_m1h1_m1h0 = smodel.VDepSReac('Na_m1h1_m1h0', ssys, \
                                slhs=[Na_m1h1], srhs=[Na_m1h0], \
                                k=lambda V:1.0e3*_b_h(V*1.0e3), vrange=Vrange)
Na_m2h1_m2h0 = smodel.VDepSReac('Na_m2h1_m2h0', ssys, \
                                slhs=[Na_m2h1], srhs=[Na_m2h0], \
                                k=lambda V:1.0e3*_b_h(V*1.0e3), vrange=Vrange)
Na_m3h1_m3h0 = smodel.VDepSReac('Na_m3h1_m3h0', ssys, \
                                slhs=[Na_m3h1], srhs=[Na_m3h0], \
                                k=lambda V:1.0e3*_b_h(V*1.0e3), vrange=Vrange)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Create ohmic current objects

OC_K = smodel.OhmicCurr('OC_K', ssys, chanstate=K_n4, g=K_G, erev=K_rev)	
OC_Na = smodel.OhmicCurr('OC_Na', ssys, chanstate=Na_m3h1, g=Na_G, erev=Na_rev)
OC_L = smodel.OhmicCurr('OC_L', ssys, chanstate=Leak, g=L_G, erev=leak_rev) 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # TETRAHEDRAL MESH  # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

mesh = meshio.importAbaqus(meshfile_ab, 1e-6)[0]
tet_hosts = gd.binTetsByAxis(mesh, steps.mpi.nhosts)
tri_hosts = gd.partitionTris(mesh, tet_hosts, mesh.getSurfTris())

# # # # # # # # # # # # # # # MESH MANIPULATION # # # # # # # # # # # # # # # # #

# Find the vertices for the current clamp and store in a list
injverts = []
for i in range(mesh.nverts):
	if ((mesh.getVertex(i)[2] < (mesh.getBoundMin()[2]+0.1e-6))):
		injverts.append(i)
if steps.mpi.rank ==0: print("Found ", injverts.__len__(), "I_inject vertices")

facetris = []
for i in range(mesh.ntris):
	tri = mesh.getTri(i) 
	if ((tri[0] in injverts) and (tri[1] in injverts) and (tri[2] in injverts)):
		facetris.append(i)
if steps.mpi.rank ==0: print("Found ", facetris.__len__(), "triangles on bottom face")

memb_tris = list(mesh.getSurfTris()) 

# Remove triangles on bottom face from membrane triangles
for t in facetris: memb_tris.remove(t)

# The points along (z) axis at which to record potential
pot_pos = numpy.arange(mesh.getBoundMin()[2], mesh.getBoundMax()[2], 10e-6)
pot_n = len(pot_pos)

pot_tet = numpy.zeros(pot_n, dtype = 'uint')

i=0
for p in pot_pos:
    # Axis is aligned with z-axis
    pot_tet[i] = mesh.findTetByPoint([0.0, 0.0, pot_pos[i]])
    i=i+1

# # # # # # # # # # # # # # # GEOMETRY OBJECTS  # # # # # # # # # # # # # # # # #

# Create cytosol compartment
cyto = sgeom.TmComp('cyto', mesh, range(mesh.ntets))

# Create the patch and associate with surface system 'ssys'
patch = sgeom.TmPatch('patch', mesh, memb_tris, cyto)
patch.addSurfsys('ssys')

# Create the membrane across which the potential will be solved
membrane = sgeom.Memb('membrane', mesh, [patch], opt_method = 1)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Create the random number generator
r = srng.create('mt19937',512)
r.initialize(int(time.time()%10000))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Create solver object using SuperLU EField solver
sim = mpi_solver.TetOpSplit(mdl, mesh, r, mpi_solver.EF_DV_SLUSYS, tet_hosts, tri_hosts)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Inject channels

surfarea = sim.getPatchArea('patch')

sim.setPatchCount('patch', 'Na_m0h0', Na_ro*surfarea*Na_facs[0])	 
sim.setPatchCount('patch', 'Na_m1h0', Na_ro*surfarea*Na_facs[1])	 
sim.setPatchCount('patch', 'Na_m2h0', Na_ro*surfarea*Na_facs[2])	 
sim.setPatchCount('patch', 'Na_m3h0', Na_ro*surfarea*Na_facs[3])	 
sim.setPatchCount('patch', 'Na_m0h1', Na_ro*surfarea*Na_facs[4])	 
sim.setPatchCount('patch', 'Na_m1h1', Na_ro*surfarea*Na_facs[5])	 
sim.setPatchCount('patch', 'Na_m2h1', Na_ro*surfarea*Na_facs[6])	 
sim.setPatchCount('patch', 'Na_m3h1', Na_ro*surfarea*Na_facs[7])

sim.setPatchCount('patch', 'K_n0', K_ro*surfarea*K_facs[0])
sim.setPatchCount('patch', 'K_n1', K_ro*surfarea*K_facs[1])			
sim.setPatchCount('patch', 'K_n2', K_ro*surfarea*K_facs[2])			
sim.setPatchCount('patch', 'K_n3', K_ro*surfarea*K_facs[3])
sim.setPatchCount('patch', 'K_n4', K_ro*surfarea*K_facs[4])

sim.setPatchCount('patch', 'Leak', L_ro * surfarea)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Set some simulation variables:

# Set dt for membrane potential calculation to 0.01ms
sim.setEfieldDT(1.0e-5)

# Initialize potential to -65mV
sim.setMembPotential('membrane', -65e-3)

# Set capacitance of the membrane to 1 uF/cm^2 = 0.01 F/m^2
sim.setMembCapac('membrane', 1.0e-2)

# Set resistivity of the conduction volume to 100 ohm.cm = 1 ohm.meter
sim.setMembVolRes('membrane', 1.0)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Set the current clamp
niverts = injverts.__len__()
for t in injverts:
    sim.setVertIClamp(t, Iclamp/niverts) 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Create result structures
if steps.mpi.rank ==0:
    res = numpy.zeros((N_timepoints, pot_n))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Run the simulation
for l in range(N_timepoints):
    if steps.mpi.rank ==0: print("Tpnt: ", l,"/", N_timepoints)

    sim.run(DT_sim*l)
    
    if steps.mpi.rank ==0:
        for p in range(pot_n):
            res[l,p] = sim.getTetV(int(pot_tet[p]))*1.0e3

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if steps.mpi.rank ==0:
    results = (res,  pot_pos)
    tpnt = arange(0.0, N_timepoints*DT_sim, DT_sim)

    for tidx in (0, 10,20,30,40):
        plot(results[1]*1e6, results[0][tidx], \
            label=str(1e3*tidx*DT_sim)+'ms', linewidth=3)
    legend(numpoints=1)
    xlim(0, 1000)
    ylim(-80,40)
    xlabel('Z-axis (um)')
    ylabel('Membrane potential (mV)')
    show()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

