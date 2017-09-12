# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Okinawa Institute of Science and Technology, Japan.
#
# This script runs on STEPS 2.x http://steps.sourceforge.net
#
# H Anwar, I Hepburn, H Nedelescu, W Chen and E De Schutter
# Stochastic calcium mechanisms cause dendritic calcium spike variability
# J Neuroscience 2013
#
# discrete.py : Simply injects channels into membrane discretely 
#
# Script author: Iain Hepburn
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from __future__ import print_function
from random import *
seed(17)

def inj_discrete(sim, memb_tris, spec_string, ninject, patch_string=''):
    print("")
    if patch_string: print("want to inject", ninject, spec_string, "into", patch_string)
    
    while (ninject > len(memb_tris)):
        for t in memb_tris: sim.setTriCount(t, spec_string, sim.getTriCount(t, spec_string)+1)
        ninject-= len(memb_tris)
    
    chosen_tris=[]
    ninjected = 0
    while (ninjected != ninject):
        tri = choice(memb_tris)
        if not tri in chosen_tris:
            sim.setTriCount(tri, spec_string, sim.getTriCount(tri, spec_string)+1)
            chosen_tris.append(tri)
            ninjected+=1
            print("Injected", sim.getTriCount(tri, spec_string), spec_string, "into tri", tri)
    print("Injected", sim.getPatchCount(patch_string, spec_string))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
