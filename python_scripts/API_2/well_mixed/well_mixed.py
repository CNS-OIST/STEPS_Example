####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   
###

# Example: Well-mixed reaction system
# http://steps.sourceforge.net/manual/well_mixed.html

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

rng = RNG('mt19937', 256, seed=23412)

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

plt.plot(saver.time[0], np.mean(saver.data, 0))
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
