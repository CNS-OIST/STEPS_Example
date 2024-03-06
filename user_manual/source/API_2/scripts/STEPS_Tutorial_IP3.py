# Example: Surface-Volume Reactions
# http://steps.sourceforge.net/manual/API_2/STEPS_Tutorial_IP3.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.saving import *
from steps.rng import *

r = ReactionManager()

#########################
# Model setup
#########################

mdl = Model()
with mdl:
    Ca, IP3, R, RIP3, Ropen, RCa, R2Ca, R3Ca, R4Ca = Species.Create()
    surfsys = SurfaceSystem.Create()

    with surfsys:
        # IP3 and activating Ca binding
        R.s    + IP3.o <r['r1']> RIP3.s
        RIP3.s + Ca.o  <r['r2']> Ropen.s
        r['r1'].K = 1000e6, 25800
        r['r2'].K = 8000e6, 2000

        # Inactivating Ca binding
        R.s    + Ca.o <r['r3']> RCa.s
        RCa.s  + Ca.o <r['r4']> R2Ca.s
        R2Ca.s + Ca.o <r['r5']> R3Ca.s
        R3Ca.s + Ca.o <r['r6']> R4Ca.s
        r['r3'].K = 8.889e6, 5
        r['r4'].K = 20e6, 10
        r['r5'].K = 40e6, 15
        r['r6'].K = 60e6, 20

        # Ca ions passing through open IP3R channel
        Ca.i + Ropen.s >r[1]> Ropen.s + Ca.o
        r[1].K = 2e8

#########################
# Geom setup
#########################

geom = Geometry()
with geom:
    # Create the cytosol and Endoplasmic Reticulum compartments
    cyt, ER = Compartment.Create()
    cyt.Vol = 1.6572e-19
    ER.Vol = 1.968e-20

    # ER is the 'inner' compartment, cyt is the 'outer' compartment
    memb = Patch.Create(ER, cyt, surfsys)
    memb.Area = 0.4143e-12

#########################
# Simulation setup
#########################

rng = RNG('mt19937', 512, 7233)

sim = Simulation('Wmdirect', mdl, geom, rng)

rs = ResultSelector(sim)

Rstates = rs.memb.MATCH('R.*').Count

Reacs = rs.memb.MATCH('r[1-6]')['fwd'].Extent + rs.memb.MATCH('r[1-6]')['bkw'].Extent

Rstates.labels = [l.split('.')[1] for l in Rstates.labels]

sim.toSave(Rstates, Reacs, dt=0.001)

#########################
# Run simulation
#########################

NITER = 100
ENDT = 0.201

for i in range (0, NITER):
    sim.newRun()

    sim.cyt.Ca.Conc = 3.30657e-8
    sim.cyt.IP3.Count = 6
    sim.ER.Ca.Conc = 150e-6
    sim.ER.Ca.Clamped = True
    sim.memb.R.Count = 160

    sim.run(ENDT)

#########################
# Plotting results
#########################

from matplotlib import pyplot as plt
import numpy as np

plt.figure(figsize=(10, 7))

RopenInd = Rstates.labels.index('Ropen')
RopenData = Rstates.data[:, :, RopenInd]

time = Rstates.time[0]
mean = np.mean(RopenData, axis=0)
std = np.std(RopenData, axis=0)

plt.plot(time, mean, linewidth=2, label='Average')
plt.fill_between(time, mean - std, mean + std, alpha=0.2, label='Std. Dev.')

for t, d in zip(Rstates.time, RopenData):
    plt.plot(t, d, color='grey', linewidth=0.1, zorder=-1)

plt.ylim(0)
plt.margins(0, 0.05)
plt.xlabel('Time [s]')
plt.ylabel('Number of open IP3R')
plt.legend()
plt.show()

plt.figure(figsize=(10, 7))

time = Rstates.time[0]
mean = np.mean(Rstates.data, axis=0)
std = np.std(Rstates.data, axis=0)

plt.plot(time, mean, linewidth=2)
for m, s in zip(mean.T, std.T):
    plt.fill_between(time, m - s, m + s, alpha=0.2)

plt.legend(Rstates.labels)
plt.xlabel('Time [s]')
plt.ylabel('Number of receptors')
plt.ylim(0)
plt.margins(0, 0.05)
plt.show()

plt.figure(figsize=(10, 7))

time = Reacs.time[0]
dt = time[1] - time[0]
meanDeriv = np.mean(np.gradient(Reacs.data, dt, axis=1), axis=0)

plt.stackplot(time, meanDeriv.T)

plt.legend([f'd{l} / dt' for l in Reacs.labels])
plt.margins(0, 0.05)
plt.xlabel('Time [s]')
plt.ylabel('Total reaction rate [1/s]')
plt.show()
