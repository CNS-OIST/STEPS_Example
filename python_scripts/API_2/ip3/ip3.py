import steps.interface

from steps.model import *
from steps.geom import *
from steps.sim import *
from steps.saving import *
from steps.rng import *

import sys
import pylab
import numpy

r = ReactionManager()

#########################
# Model setup
#########################
mdl = Model()
with mdl:
    # Species 
    Ca, IP3, R, RIP3, Ropen, RCa, R2Ca, R3Ca, R4Ca = Species.Create()
    # Surface system
    surfsys = SurfaceSystem.Create()
    # Reactions
    with surfsys:
        # IP3 binding
        (IP3.o + R.s <r[1]> RIP3.s) + Ca.o <r[2]> Ropen.s
        r[1].K = 1000e6, 25800
        r[2].K = 8000e6, 2000

        # Ca binding
        (((R.s + Ca.o <r[7]> RCa.s) + Ca.o <r[8]> R2Ca.s) + Ca.o <r[3]> R3Ca.s) + Ca.o <r[4]> R4Ca.s
        r[7].K = 8.889e6, 5
        r[8].K = 20e6, 10
        r[3].K = 40e6, 15
        r[4].K = 60e6, 20

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

    print('Inner compartment to memb is', memb.innerComp.name)
    print('Outer compartment to memb is', memb.outerComp.name)

#########################
# RNG setup
#########################
rng = RNG('mt19937', 512, 7233)

#########################
# Solver setup
#########################
sim = Simulation('Wmdirect', mdl, geom, rng)
        
NITER = 100
ENDT = 0.201

rs = ResultSelector(sim)

RopenDat = rs.memb.Ropen.Count

CaConc = rs.cyt.Ca.Conc

sim.toSave(RopenDat, CaConc, dt=0.001)

#########################
# Run simulation
#########################

for i in range (0, NITER):
    sim.newRun()

    sim.cyt.Ca.Conc = 3.30657e-8
    sim.cyt.IP3.Count = 6
    sim.ER.Ca.Conc = 150e-6
    sim.ER.Ca.Clamped = True
    sim.memb.R.Count = 160

    sim.run(ENDT)

    pylab.plot(RopenDat.time[-1], RopenDat.data[-1], color = 'blue', linewidth = 0.1)

datTime = RopenDat.time[0]

res_mean = numpy.mean(RopenDat.data, 0)
res_std = numpy.std(RopenDat.data, 0)
res_std1 = res_mean + res_std
res_std2 = res_mean - res_std

pylab.plot(datTime, res_mean, color = 'black', linewidth = 2.0, label = 'mean')
pylab.plot(datTime, res_std1, color = 'gray', linewidth = 1.0, label='std')
pylab.plot(datTime, res_std2,color = 'gray', linewidth = 1.0)

pylab.xlabel('Time (sec)')
pylab.ylabel('# IP3 receptors in open state')
pylab.title('IP3 receptor model: %d iterations with Wmdirect'%NITER)
pylab.ylim(0)
pylab.legend()

pylab.show()
        
