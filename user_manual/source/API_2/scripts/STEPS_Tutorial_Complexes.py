# Example: Multi-state complexes / 7.3. IP3 receptor model
# http://steps.sourceforge.net/manual/API_2/STEPS_Tutorial_Complexes.html#IP3-receptor-model

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.saving import *
from steps.rng import *

nAvog = 6.02214076e23

nbIP3R = 5
nbPumps = 5

c0 = 2e-6
c1 = 0.185

cytVol = 1.6572e-19
ERVol = cytVol * c1

a1 = 400e6
a2 = 0.2e6
a3 = 400e6
a4 = 0.2e6
a5 = 20e6

b1 = 0.13e-6 * a1
b2 = 1.049e-6 * a2
b3 = 943.4e-9 * a3
b4 = 144.5e-9 * a4
b5 = 82.34e-9 * a5

v1 = 6
v2 = 0.11

v3 = 0.9e-6
k3 = 0.1e-6

rp = v3 * 1e3 * cytVol * nAvog / nbPumps / 2
rb = 10 * rp
rf = (rb + rp) / (k3 ** 2)

kip3 = 1e3 * nAvog * ERVol * v1 / nbIP3R

#########################
# Model setup
#########################

mdl = Model()
r = ReactionManager()

with mdl:
    Ca, IP3, ERPump, ERPump2Ca = Species.Create()

    R000, R100, R010, R001, R110, R101, R111, R011 = SubUnitState.Create()

    IP3RSU = SubUnit.Create([R000, R100, R010, R001, R110, R101, R111, R011])

    IP3R = Complex.Create([IP3RSU, IP3RSU, IP3RSU, IP3RSU], statesAsSpecies=True)

    ssys = SurfaceSystem.Create()

    with ssys:
        # Ca2+ passing through open IP3R channel
        IP3R_1 = IP3R.get()
        IP3R_1[R110, R110, R110, :].s + Ca.i <r['caflx']> IP3R_1[R110, R110, R110, :].s + Ca.o
        r['caflx'].K = kip3, kip3

        # IP3R subunits reaction network
        with IP3R[...]:
            R000.s + IP3.o <r[1]> R100.s
            R000.s + Ca.o <r[2]> R010.s
            R000.s + Ca.o <r[3]> R001.s

            R100.s + Ca.o <r[4]> R110.s
            R100.s + Ca.o <r[5]> R101.s

            R010.s + IP3.o <r[6]> R110.s
            R010.s + Ca.o <r[7]> R011.s

            R001.s + IP3.o <r[8]> R101.s
            R001.s + Ca.o <r[9]> R011.s

            R110.s + Ca.o <r[10]> R111.s
            R101.s + Ca.o <r[11]> R111.s

            R011.s + IP3.o <r[12]> R111.s

            r[1].K = a1, b1
            r[2].K = a5, b5
            r[3].K = a4, b4
            r[4].K = a5, b5
            r[5].K = a2, b2
            r[6].K = a1, b1
            r[7].K = a4, b4
            r[8].K = a3, b3
            r[9].K = a5, b5
            r[10].K = a2, b2
            r[11].K = a5, b5
            r[12].K = a3, b3

        # Ca2+ leak
        Ca.i <r[1]> Ca.o
        r[1].K = v2, c1 * v2

        2*Ca.o + ERPump.s <r[1]> ERPump2Ca.s >r[2]> 2*Ca.i + ERPump.s
        r[1].K = rf, rb
        r[2].K = rp

#########################
# Geom setup
#########################

geom = Geometry()

with geom:
    cyt, ER = Compartment.Create()
    cyt.Vol = cytVol
    ER.Vol = ERVol

    memb = Patch.Create(ER, cyt, ssys)
    memb.Area = 0.4143e-12

#########################
# Simulation setup
#########################

rng = RNG('mt19937', 512, 7233)

sim = Simulation('Wmdirect', mdl, geom, rng)

rs = ResultSelector(sim)

cytCa = rs.cyt.Ca.Conc

caFlux = rs.SUM(rs.memb.caflx['fwd'].Extent) << rs.SUM(rs.memb.caflx['bkw'].Extent)

IP3RStates =   rs.memb.IP3R[~R110, ~R110, ~R110, ~R110].Count
IP3RStates <<= rs.memb.IP3R[ R110, ~R110, ~R110, ~R110].Count
IP3RStates <<= rs.memb.IP3R[ R110,  R110, ~R110, ~R110].Count
IP3RStates <<= rs.memb.IP3R[ R110,  R110,  R110, ~R110].Count
IP3RStates <<= rs.memb.IP3R[ R110,  R110,  R110,  R110].Count

sim.toSave(cytCa, caFlux, IP3RStates, dt=0.05)

#########################
# Run simulation
#########################

ENDT = 10.0

sim.newRun()

# Initial conditions
sim.cyt.Ca.Conc = 3.30657e-8
sim.cyt.IP3.Conc = 0.2e-6
sim.ER.Ca.Conc = c0/c1
sim.memb.ERPump.Count = nbPumps
sim.memb.IP3R[R000, R000, R000, R000].Count = nbIP3R

sim.run(ENDT)

#########################
# Plotting results
#########################

from matplotlib import pyplot as plt
import numpy as np

plt.figure(figsize=(10, 7))
plt.plot(cytCa.time[0], cytCa.data[0]*1e6)
plt.legend(cytCa.labels)
plt.xlabel('Time [s]')
plt.ylabel('Concentration [μM]')
plt.show()

plt.figure(figsize=(10, 7))
plt.plot(caFlux.time[0], caFlux.data[0])
plt.legend(caFlux.labels)
plt.xlabel('Time [s]')
plt.ylabel('Reaction extent')
plt.show()

# # # # # # # # # # # # # # # # # #

n = 20

plt.figure(figsize=(10, 7))
for i in range(IP3RStates.data[0].shape[1]):
    sig = IP3RStates.data[0, :, i]
    avg = np.convolve(sig, np.ones(n) / n, 'valid')
    tme = IP3RStates.time[0, n//2:-n//2+1]

    plt.plot(tme, avg, color=f'C{i}', label=IP3RStates.labels[i])
    plt.plot(IP3RStates.time[0], sig, '--', linewidth=1, color=f'C{i}', alpha=0.4)

plt.legend(loc=1)
plt.xlabel('Time [s]')
plt.ylabel('Count')
plt.show()
