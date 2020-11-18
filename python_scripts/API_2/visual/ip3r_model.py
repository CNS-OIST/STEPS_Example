# IP3 receptor model

import steps.interface

from steps.model import *
from steps.geom import *

###############################################################################

def getModel():
    mdl = Model()
    r = ReactionManager()

    with mdl:
        # chemical species objects
        Ca, IP3, R, RIP3, Ropen, RCa, R2Ca, R3Ca, R4Ca = Species.Create()
        
        ssys = SurfaceSystem.Create()

        with ssys:
            
            # Binding reactions
            (R.s + IP3.o <r[1]> RIP3.s) + Ca.o <r[2]> Ropen.s
            r[1].setRates(1000000000.0, 25800.0)
            r[2].setRates(8000000000.0, 2000.0)

            (((R.s + Ca.o <r[3]> RCa.s) + Ca.o <r[4]> R2Ca.s) + Ca.o <r[5]> R3Ca.s) + Ca.o <r[6]> R4Ca.s
            r[3].setRates(8889000.0, 5.0)
            r[4].setRates(20000000.0, 10.0)
            r[5].setRates(40000000.0, 15.0)
            r[6].setRates(60000000.0, 20.0)
            
            # Ca ions passing through open IP3R channel
            Ca.i + Ropen.s <r[1]> Ropen.s + Ca.o
            # Corresponds to Ca input ~ 20000/ms for open receptor
            r[1].setRates(8000000.0, 8000000.0)          
    
    return mdl

###############################################################################

def getGeom():
    mesh  = Geometry()

    with mesh:
        cyt = Compartment.Create(vol=0.1e-18)
        ER = Compartment.Create(vol=0.02e-18)
        
        memb = Patch.Create(ER, cyt, 'ssys')
    
    return mesh

###############################################################################
