# Example: Well-mixed reaction system
# http://steps.sourceforge.net/manual/API_2/STEPS_Tutorial_wm.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.saving import *
from steps.rng import *

from matplotlib import pyplot as plt
import numpy as np

mdl = Model()

r = ReactionManager()

with mdl:
    molA, molB, molC = Species.Create()

    vsys = VolumeSystem.Create()

    with vsys:
        molA + molB <r['r1']> molC
        r['r1'].K = 0.3e6, 0.7

geom = Geometry()

with geom:
    comp = Compartment.Create(vsys, 1.6667e-21)

rng = RNG('mt19937', 256, 1234)

sim = Simulation('Wmdirect', mdl, geom, rng)

# Manual data saving

sim.newRun()

sim.comp.molA.Conc = 31.4e-6
sim.comp.molB.Conc = 22.3e-6

tpnts = np.arange(0.0, 2.001, 0.001)
values = []
for t in tpnts:
    sim.run(t)
    values.append(sim.comp.LIST(molA, molB, molC).Count)

plt.figure(figsize=(12, 8))
plt.plot(tpnts, values)
plt.legend(['molA', 'molB', 'molC'])
plt.xlabel('Time [s]')
plt.ylabel('#molecules')
plt.show()

# Automatic data saving

rs = ResultSelector(sim)

saver = rs.comp.LIST(molA, molB, molC).Count

sim.toSave(saver, dt=0.001)

NITER = 100

for i in range(NITER):
    sim.newRun()

    sim.comp.molA.Conc = 31.4e-6
    sim.comp.molB.Conc = 22.3e-6

    sim.run(2.0)

plt.plot(saver.time[0], saver.data[0])
plt.legend(saver.labels)
plt.xlabel('Time [s]')
plt.show()

plt.plot(saver.time[0], np.mean(saver.data, axis=0))
plt.legend(saver.labels)
plt.xlabel('Time [s]')
plt.show()

#####

for i in range(NITER):
    sim.newRun()

    sim.comp.molA.Conc = 31.4e-6
    sim.comp.molB.Conc = 22.3e-6

    sim.run(1.0)

    # Add 10 molecules of species A
    sim.comp.molA.Count += 10

    sim.run(2.0)

plt.plot(saver.time[0], np.mean(saver.data[NITER:], 0))
plt.legend(saver.labels)
plt.xlabel('Time [s]')
plt.show()

#####

sim.newRun()
sim.comp.molA.Conc = 31.4e-6
sim.comp.molB.Conc = 22.3e-6

sim.run(0.1)

sim.comp.molA.Clamped = True

sim.run(0.6)

sim.comp.molA.Clamped = False

sim.run(2.0)
    
plt.plot(saver.time[0], saver.data[-1])
plt.legend(saver.labels)
plt.xlabel('Time [s]')
plt.show()

#####

for i in range(NITER):
    sim.newRun()
    sim.comp.molA.Conc = 31.4e-6
    sim.comp.molB.Conc = 22.3e-6

    sim.run(2.0)

    sim.comp.r1['fwd'].Active = False
    sim.run(4.0)
    sim.comp.r1['fwd'].Active = True
    sim.run(6.0)
    sim.comp.r1['bkw'].Active = False
    sim.run(8.0)
    sim.comp.r1['bkw'].Active = True
    sim.run(10.0)
    sim.comp.r1['fwd'].Active = False
    sim.comp.r1['bkw'].Active = False
    sim.run(12.0)

plt.plot(saver.time[-1], np.mean(saver.data[-NITER:], 0))
plt.legend(saver.labels)
plt.xlabel('Time [s]')
plt.show()
