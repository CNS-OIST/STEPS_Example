import steps.interface

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
# *StochasticHH.py : The stochastic Hodgkin-Huxley model, used in the 
# above study. 
#
# Script authors: Haroon Anwar and Iain Hepburn
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# USAGE
#
# $ python StochasticHH.py *mesh* *root* *iter_n* 
#  
#  *mesh* is the tetrahedral mesh (10um to 160um cylinder)
#  *root* is the path to the location for data storage 
#  *iter_n* (is intened to be an integer) is an identifier number for each 
#     simulation iteration.
#
# E.g: 
# $ python StochasticHH.py Cylinder2_dia2um_L10um_outer0_3um_0.3shell_0.3size_19156tets_adaptive.inp ~/stochHHsims/ 1
#
#
# OUTPUT 
#
# In (root)/data/StochasticHH/(mesh)/(iter_n+time) directory 
# 2 data files will be recorded. Each file contains one row for every 
# time-point at which data is recorded, organised into the following columns:
# 
# currents.dat
# Time (ms), Na current, K current, leak current
# (current units are Amps/m^2)
#
# voltage.dat
# Time (ms), voltage at mesh centre (mV)
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


from __future__ import print_function
import math
# WARNING: Using a variable name that is reserved (['time']).
import time
from random import *
from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *

import meshes.gettets as gettets
from extra.constants_hh import *

import sys
import os

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

meshfile_ab, root, iter_n = sys.argv[1], sys.argv[2], sys.argv[3]

if meshfile_ab == 'Cylinder2_dia2um_L160um_outer0_0.3shell_0.3size_279152tets_adaptive.inp':
    cyl160=True
else:
    cyl160=False

########################### BIOCHEMICAL MODEL ###############################

mdl = Model()
# WARNING: Using a variable name that is reserved (['r']).
r = ReactionManager()
with mdl:
    
    #surface systems (No need for a vol sys)                                                                                                                                                                                                                                                            
    ssys = SurfaceSystem.Create()
    K_n0, K_n1, K_n2, K_n3, K_n4 = SubUnitState.Create()
    
    # Potassium channel
    Kchan = Channel.Create([K_n0, K_n1, K_n2, K_n3, K_n4])
with ssys, mdl, Kchan[...]:
    
    
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    K_n0.s <r['Kn0n1']> K_n1.s ; r['Kn0n1'].setRates(VDepRate(lambda V: 1.0e3 *4.*a_n(V*1.0e3)* Qt, vrange=Vrange), VDepRate(lambda V: 1.0e3 *1.*b_n(V*1.0e3)* Qt, vrange=Vrange))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    K_n1.s <r['Kn1n2']> K_n2.s ; r['Kn1n2'].setRates(VDepRate(lambda V: 1.0e3 *3.*a_n(V*1.0e3)* Qt, vrange=Vrange), VDepRate(lambda V: 1.0e3 *2.*b_n(V*1.0e3)* Qt, vrange=Vrange))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    K_n2.s <r['Kn2n3']> K_n3.s ; r['Kn2n3'].setRates(VDepRate(lambda V: 1.0e3 *2.*a_n(V*1.0e3)* Qt, vrange=Vrange), VDepRate(lambda V: 1.0e3 *3.*b_n(V*1.0e3)* Qt, vrange=Vrange))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    
    K_n3.s <r['Kn3n4']> K_n4.s ; r['Kn3n4'].setRates(VDepRate(lambda V: 1.0e3 *1.*a_n(V*1.0e3)* Qt, vrange=Vrange), VDepRate(lambda V: 1.0e3 *4.*b_n(V*1.0e3)* Qt, vrange=Vrange))
with ssys:
    
    OC_K = OhmicCurr.Create(Kchan[K_n4], K_G, K_rev)
with mdl:
    Na_m0h0, Na_m1h0, Na_m2h0, Na_m3h0, Na_m0h1, Na_m1h1, Na_m2h1, Na_m3h1 = SubUnitState.Create()
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #                
    
    # Sodium channel
    Nachan = Channel.Create([Na_m0h0, Na_m1h0, Na_m2h0, Na_m3h0, Na_m0h1, Na_m1h1, Na_m2h1, Na_m3h1])
with ssys, mdl, Nachan[...]:
    
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    Na_m0h1.s <r['Na_m0h1_m1h1']> Na_m1h1.s ; r['Na_m0h1_m1h1'].setRates(VDepRate(lambda V:1.0e3*3.*a_m(V*1.0e3)* Qt, vrange=Vrange), VDepRate(lambda V:1.0e3*1.*b_m(V*1.0e3)* Qt, vrange=Vrange))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    Na_m1h1.s <r['Na_m1h1_m2h1']> Na_m2h1.s ; r['Na_m1h1_m2h1'].setRates(VDepRate(lambda V:1.0e3*2.*a_m(V*1.0e3)* Qt, vrange=Vrange), VDepRate(lambda V:1.0e3*2.*b_m(V*1.0e3)* Qt, vrange=Vrange))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    
    Na_m2h1.s <r['Na_m2h1_m3h1']> Na_m3h1.s ; r['Na_m2h1_m3h1'].setRates(VDepRate(lambda V:1.0e3*1.*a_m(V*1.0e3)* Qt, vrange=Vrange), VDepRate(lambda V:1.0e3*3.*b_m(V*1.0e3)* Qt, vrange=Vrange))
    
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    Na_m0h0.s <r['Na_m0h0_m1h0']> Na_m1h0.s ; r['Na_m0h0_m1h0'].setRates(VDepRate(lambda V:1.0e3*3.*a_m(V*1.0e3)* Qt, vrange=Vrange), VDepRate(lambda V:1.0e3*1.*b_m(V*1.0e3)* Qt, vrange=Vrange))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    Na_m1h0.s <r['Na_m1h0_m2h0']> Na_m2h0.s ; r['Na_m1h0_m2h0'].setRates(VDepRate(lambda V:1.0e3*2.*a_m(V*1.0e3)* Qt, vrange=Vrange), VDepRate(lambda V:1.0e3*2.*b_m(V*1.0e3)* Qt, vrange=Vrange))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    
    Na_m2h0.s <r['Na_m2h0_m3h0']> Na_m3h0.s ; r['Na_m2h0_m3h0'].setRates(VDepRate(lambda V:1.0e3*1.*a_m(V*1.0e3)* Qt, vrange=Vrange), VDepRate(lambda V:1.0e3*3.*b_m(V*1.0e3)* Qt, vrange=Vrange))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    
    Na_m0h0.s <r['Na_m0h0_m0h1']> Na_m0h1.s ; r['Na_m0h0_m0h1'].setRates(VDepRate(lambda V:1.0e3*a_h(V*1.0e3)* Qt, vrange=Vrange), VDepRate(lambda V:1.0e3*b_h(V*1.0e3)* Qt, vrange=Vrange))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    Na_m1h0.s <r['Na_m1h0_m1h1']> Na_m1h1.s ; r['Na_m1h0_m1h1'].setRates(VDepRate(lambda V:1.0e3*a_h(V*1.0e3)* Qt, vrange=Vrange), VDepRate(lambda V:1.0e3*b_h(V*1.0e3)* Qt, vrange=Vrange))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    Na_m2h0.s <r['Na_m2h0_m2h1']> Na_m2h1.s ; r['Na_m2h0_m2h1'].setRates(VDepRate(lambda V:1.0e3*a_h(V*1.0e3)* Qt, vrange=Vrange), VDepRate(lambda V:1.0e3*b_h(V*1.0e3)* Qt, vrange=Vrange))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    Na_m3h0.s <r['Na_m3h0_m3h1']> Na_m3h1.s ; r['Na_m3h0_m3h1'].setRates(VDepRate(lambda V:1.0e3*a_h(V*1.0e3)* Qt, vrange=Vrange), VDepRate(lambda V:1.0e3*b_h(V*1.0e3)* Qt, vrange=Vrange))
with ssys:
    
    OC_Na = OhmicCurr.Create(Nachan[Na_m3h1], Na_G, Na_rev)
with mdl:
    Leak = SubUnitState.Create()
    
    # Leak channel
    L = Channel.Create([Leak])
with ssys:
    
    OC_L = OhmicCurr.Create(L[Leak], L_G, leak_rev)


##################################

########### MESH & COMPARTMENTALIZATION #################

##########Import Mesh

mesh = TetMesh.Load('./meshes/'+meshfile_ab)

outer_tets = range(len(mesh.tets))

###USE OF gettets
#getcyl(tetmesh, rad,  zmin, zmax, binnum=120, x = 0.0, y = 0.0):
inner_tets = gettets.getcyl(mesh, 1e-6, -200e-6, 200e-6)[0]

for i in inner_tets: outer_tets.remove(i)
assert(outer_tets.__len__() + inner_tets.__len__() == len(mesh.tets))

print(outer_tets.__len__(), " tets in outer compartment")
print(inner_tets.__len__(), " tets in inner compartment")

# Record voltage from the central tetrahedron
# WARNING: findTetByPoint was replaced by the new TetList syntax and trying to access a point outside the mesh will now raise an exception.
cent_tet = mesh.tets[0.0, 0.0, 0.0].idx
with mesh:
    
    ########## Create an intracellular compartment i.e. cytosolic compartment
    
    cyto = TetComp.Create(inner_tets)
    
    if cyl160:
        # Ensure that we use points a small distance inside the boundary:
        LENGTH = mesh.bbox.max[2] - mesh.bbox.min[2]
        boundminz = mesh.bbox.min[2] + LENGTH/len(mesh.tets)
        boundmaxz = mesh.bbox.max[2] - LENGTH/len(mesh.tets)

        memb_tris = list(mesh.surface.indices)
        minztris = []
        maxztris = []
        for tri in memb_tris:
            zminboundtri = True
            zmaxboundtri = True
            tritemp = mesh.tris[tri].verts.indices
            trizs = [0.0, 0.0, 0.0]
            trizs[0] = mesh.verts[tritemp[0]][2]
            trizs[1] = mesh.verts[tritemp[1]][2]
            trizs[2] = mesh.verts[tritemp[2]][2]
            for j in range(3):
                if (trizs[j]>boundminz):
                    zminboundtri = False
            if (zminboundtri):
                minztris.append(tri)
                continue
            for j in range(3):
                if (trizs[j]< boundmaxz):
                    zmaxboundtri = False
            if (zmaxboundtri):
                maxztris.append(tri)

        for t in minztris: memb_tris.remove(t)
        for t in maxztris: memb_tris.remove(t)
        
    else:
        print('Finding connecting triangles...')
        out_tris = set()
        for i in outer_tets:
                tritemp = mesh.tets[i].faces.indices
                for j in range(4): out_tris.add(tritemp[j])

        in_tris = set()
        for i in inner_tets:
                tritemp = mesh.tets[i].faces.indices
                for j in range(4): in_tris.add(tritemp[j])

        memb_tris = out_tris.intersection(in_tris)
        memb_tris = list(memb_tris)
    
    print(len(memb_tris), " surface triangles.")
    
    ########## Create a membrane as a surface mesh
    memb = TetPatch.Create(memb_tris, cyto, None, 'ssys')
    
    print("Area: ", memb.getArea())
    
    print("Creating membrane..")
    membrane = Membrane.Create([memb])
print("Membrane created.")

# # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

# WARNING: Using a variable name that is reserved (['r']).
r = RNG('mt19937', 512, 7)

# WARNING: Using a variable name that is reserved (['r']).
sim = Simulation('Tetexact', mdl, mesh, r, calcMembPot=True)

print("Resetting simulation object..")
sim.newRun()

print("Injecting molecules..")

sim.Temp = TEMPERATURE+273.15

surfarea = sim.memb.Area

sim.memb.L[Leak].Count = round(L_ro * surfarea)

sim.memb.Nachan[Na_m0h0].Count = round(Na_ro*surfarea*Na_facs[0])
sim.memb.Nachan[Na_m1h0].Count = round(Na_ro*surfarea*Na_facs[1])
sim.memb.Nachan[Na_m2h0].Count = round(Na_ro*surfarea*Na_facs[2])
sim.memb.Nachan[Na_m3h0].Count = round(Na_ro*surfarea*Na_facs[3])
sim.memb.Nachan[Na_m0h1].Count = round(Na_ro*surfarea*Na_facs[4])
sim.memb.Nachan[Na_m1h1].Count = round(Na_ro*surfarea*Na_facs[5])
sim.memb.Nachan[Na_m2h1].Count = round(Na_ro*surfarea*Na_facs[6])
sim.memb.Nachan[Na_m3h1].Count = round(Na_ro*surfarea*Na_facs[7])

sim.memb.Kchan[K_n0].Count = round(K_ro*surfarea*K_facs[0])
sim.memb.Kchan[K_n1].Count = round(K_ro*surfarea*K_facs[1])
sim.memb.Kchan[K_n2].Count = round(K_ro*surfarea*K_facs[2])
sim.memb.Kchan[K_n3].Count = round(K_ro*surfarea*K_facs[3])
sim.memb.Kchan[K_n4].Count = round(K_ro*surfarea*K_facs[4])

print('Leak', round(L_ro * surfarea))

print('Na_m0h0', round(Na_ro*surfarea*Na_facs[0]))
print('Na_m1h0', round(Na_ro*surfarea*Na_facs[1]))
print('Na_m2h0', round(Na_ro*surfarea*Na_facs[2]))
print('Na_m3h0', round(Na_ro*surfarea*Na_facs[3]))
print('Na_m0h1', round(Na_ro*surfarea*Na_facs[4]))
print('Na_m1h1', round(Na_ro*surfarea*Na_facs[5]))
print('Na_m2h1', round(Na_ro*surfarea*Na_facs[6]))
print('Na_m3h1', round(Na_ro*surfarea*Na_facs[7]))

print('K_n0', round(K_ro*surfarea*K_facs[0]))
print('K_n1', round(K_ro*surfarea*K_facs[1]))
print('K_n2', round(K_ro*surfarea*K_facs[2]))
print('K_n3', round(K_ro*surfarea*K_facs[3]))
print('K_n4', round(K_ro*surfarea*K_facs[4]))

print("Targeted Injection: ", round(Na_ro*surfarea), "Na channels")

print("Targeted Injection: ", round(K_ro*surfarea), "K channels")

print("Targeted Injection: ", round(L_ro*surfarea), "Leak channels")

sim.EfieldDT = EF_DT

sim.membrane.Potential = init_pot

sim.membrane.VolRes = Ra

sim.membrane.Capac = memb_capac


#### Recording #####

# WARNING: Using a variable name that is reserved (['time']).
c=time.ctime()

dc = c.split()[1]+c.split()[2]+'_'+c.split()[3]+'_'+c.split()[4]
dc= dc.replace(':', '_')

try: os.mkdir(root+'data')
except: pass
try: os.mkdir(root+'data/' +  'StochasticHH')
except: pass
try: os.mkdir(root+'data/' +  'StochasticHH/'+meshfile_ab)
except: pass 

os.mkdir(root+'data/' +  'StochasticHH/'+meshfile_ab+'/'+iter_n+'__'+dc )


datfile =  open(root+'data/' +  'StochasticHH/'+meshfile_ab+'/'+iter_n+'__'+dc + '/currents.dat', 'w')
datfile2 = open(root+'data/' +  'StochasticHH/'+meshfile_ab+'/'+iter_n+'__'+dc + '/voltage.dat', 'w')

# WARNING: Using a variable name that is reserved (['r']).
r.initialize(100*int(iter_n))

for l in range(NTIMEPOINTS):
    print("Tpnt: ", l)

    # WARNING: Using a variable name that is reserved (['run']).
    sim.run(TIMECONVERTER*l)

    tcur_Na = 0.0
    tcur_K = 0.0
    tcur_L = 0.0

    for m in memb_tris:
        tcur_Na = tcur_Na + sim.TRI(m).OC_Na.I 
        tcur_K = tcur_K + sim.TRI(m).OC_K.I
        tcur_L = tcur_L + sim.TRI(m).OC_L.I
    
    datfile.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile.write('%.6g' %((tcur_Na*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_K*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_L*1.0e-1)/surfarea) + ' ')
    datfile.write('\n')
    
    datfile2.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile2.write('%.6g' %(sim.TET(cent_tet).V*1.0e3) + ' ')
    datfile2.write('\n')

## END
