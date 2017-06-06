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


import math
import time
from random import *
import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.utilities.meshio as meshio
import steps.solver as ssolver

import meshes.gettets as gettets
from extra.constants_hh import *

import sys
import os

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

meshfile_ab, root, iter_n = sys.argv[1], sys.argv[2], sys.argv[3]

if meshfile_ab == 'Cylinder2_dia2um_L160um_outer0_0.3shell_0.3size_279152tets_adaptive.inp': cyl160=True
else: cyl160=False

########################### BIOCHEMICAL MODEL ###############################

mdl = smodel.Model()

#surface systems (No need for a vol sys)                                                                                                                                                                                                                                                            
ssys = smodel.Surfsys('ssys', mdl)

# Potassium channel
Kchan = smodel.Chan('Kchan', mdl)
K_n0 = smodel.ChanState('K_n0', mdl, Kchan)
K_n1 = smodel.ChanState('K_n1', mdl, Kchan)
K_n2 = smodel.ChanState('K_n2', mdl, Kchan)
K_n3 = smodel.ChanState('K_n3', mdl, Kchan)
K_n4 = smodel.ChanState('K_n4', mdl, Kchan)


Kn0n1 = smodel.VDepSReac('Kn0n1', ssys, slhs = [K_n0], srhs = [K_n1], \
                            k=lambda V: 1.0e3 *4.*a_n(V*1.0e3)* Qt, vrange = Vrange)
Kn1n2 = smodel.VDepSReac('Kn1n2', ssys, slhs = [K_n1], srhs = [K_n2], \
                            k=lambda V: 1.0e3 *3.*a_n(V*1.0e3)* Qt, vrange = Vrange)
Kn2n3 = smodel.VDepSReac('Kn2n3', ssys, slhs = [K_n2], srhs = [K_n3], \
                            k=lambda V: 1.0e3 *2.*a_n(V*1.0e3)* Qt, vrange = Vrange)
Kn3n4 = smodel.VDepSReac('Kn3n4', ssys, slhs = [K_n3], srhs = [K_n4], \
                            k=lambda V: 1.0e3 *1.*a_n(V*1.0e3)* Qt, vrange = Vrange)

Kn4n3 = smodel.VDepSReac('Kn4n3', ssys, slhs = [K_n4], srhs = [K_n3], \
                            k=lambda V: 1.0e3 *4.*b_n(V*1.0e3)* Qt, vrange = Vrange)
Kn3n2 = smodel.VDepSReac('Kn3n2', ssys, slhs = [K_n3], srhs = [K_n2], \
                            k=lambda V: 1.0e3 *3.*b_n(V*1.0e3)* Qt, vrange = Vrange)
Kn2n1 = smodel.VDepSReac('Kn2n1', ssys, slhs = [K_n2], srhs = [K_n1], \
                            k=lambda V: 1.0e3 *2.*b_n(V*1.0e3)* Qt, vrange = Vrange)
Kn1n0 = smodel.VDepSReac('Kn1n0', ssys, slhs = [K_n1], srhs = [K_n0], \
                            k=lambda V: 1.0e3 *1.*b_n(V*1.0e3)* Qt, vrange = Vrange)

OC_K = smodel.OhmicCurr('OC_K', ssys, chanstate=K_n4, g=K_G, erev=K_rev)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #                

# Sodium channel
Nachan = smodel.Chan('Nachan', mdl)
Na_m0h0 = smodel.ChanState('Na_m0h0', mdl, Nachan)
Na_m1h0 = smodel.ChanState('Na_m1h0', mdl, Nachan)
Na_m2h0 = smodel.ChanState('Na_m2h0', mdl, Nachan)
Na_m3h0 = smodel.ChanState('Na_m3h0', mdl, Nachan)
Na_m0h1 = smodel.ChanState('Na_m0h1', mdl, Nachan)
Na_m1h1 = smodel.ChanState('Na_m1h1', mdl, Nachan)
Na_m2h1 = smodel.ChanState('Na_m2h1', mdl, Nachan)
Na_m3h1 = smodel.ChanState('Na_m3h1', mdl, Nachan)

Na_m0h1_m1h1 = smodel.VDepSReac('Na_m0h1_m1h1', ssys, \
                                slhs=[Na_m0h1], srhs=[Na_m1h1], \
                                k=lambda V:1.0e3*3.*a_m(V*1.0e3)* Qt, vrange=Vrange)
Na_m1h1_m2h1 = smodel.VDepSReac('Na_m1h1_m2h1', ssys, \
                                slhs=[Na_m1h1], srhs=[Na_m2h1], \
                                k=lambda V:1.0e3*2.*a_m(V*1.0e3)* Qt, vrange=Vrange)
Na_m2h1_m3h1 = smodel.VDepSReac('Na_m2h1_m3h1', ssys, \
                                slhs=[Na_m2h1], srhs=[Na_m3h1], \
                                k=lambda V:1.0e3*1.*a_m(V*1.0e3)* Qt, vrange=Vrange)

Na_m3h1_m2h1 = smodel.VDepSReac('Na_m3h1_m2h1', ssys, \
                                slhs=[Na_m3h1], srhs=[Na_m2h1], \
                                k=lambda V:1.0e3*3.*b_m(V*1.0e3)* Qt, vrange=Vrange)
Na_m2h1_m1h1 = smodel.VDepSReac('Na_m2h1_m1h1', ssys, \
                                slhs=[Na_m2h1], srhs=[Na_m1h1], \
                                k=lambda V:1.0e3*2.*b_m(V*1.0e3)* Qt, vrange=Vrange)
Na_m1h1_m0h1 = smodel.VDepSReac('Na_m1h1_m0h1', ssys, \
                                slhs=[Na_m1h1], srhs=[Na_m0h1], \
                                k=lambda V:1.0e3*1.*b_m(V*1.0e3)* Qt, vrange=Vrange)

Na_m0h0_m1h0 = smodel.VDepSReac('Na_m0h0_m1h0', ssys, \
                                slhs=[Na_m0h0], srhs=[Na_m1h0], \
                                k=lambda V:1.0e3*3.*a_m(V*1.0e3)* Qt, vrange=Vrange)
Na_m1h0_m2h0 = smodel.VDepSReac('Na_m1h0_m2h0', ssys, \
                                slhs=[Na_m1h0], srhs=[Na_m2h0], \
                                k=lambda V:1.0e3*2.*a_m(V*1.0e3)* Qt, vrange=Vrange)
Na_m2h0_m3h0 = smodel.VDepSReac('Na_m2h0_m3h0', ssys, \
                                slhs=[Na_m2h0], srhs=[Na_m3h0], \
                                k=lambda V:1.0e3*1.*a_m(V*1.0e3)* Qt, vrange=Vrange)

Na_m3h0_m2h0 = smodel.VDepSReac('Na_m3h0_m2h0', ssys, \
                                slhs=[Na_m3h0], srhs=[Na_m2h0], \
                                k=lambda V:1.0e3*3.*b_m(V*1.0e3)* Qt, vrange=Vrange)
Na_m2h0_m1h0 = smodel.VDepSReac('Na_m2h0_m1h0', ssys, \
                                slhs=[Na_m2h0], srhs=[Na_m1h0], \
                                k=lambda V:1.0e3*2.*b_m(V*1.0e3)* Qt, vrange=Vrange)
Na_m1h0_m0h0 = smodel.VDepSReac('Na_m1h0_m0h0', ssys, \
                                slhs=[Na_m1h0], srhs=[Na_m0h0], \
                                k=lambda V:1.0e3*1.*b_m(V*1.0e3)* Qt, vrange=Vrange)
Na_m0h0_m0h1 = smodel.VDepSReac('Na_m0h0_m0h1', ssys, \
                                slhs=[Na_m0h0], srhs=[Na_m0h1], \
				k=lambda V:1.0e3*a_h(V*1.0e3)* Qt, vrange=Vrange)
Na_m1h0_m1h1 = smodel.VDepSReac('Na_m1h0_m1h1', ssys, \
                                slhs=[Na_m1h0], srhs=[Na_m1h1], \
                                k=lambda V:1.0e3*a_h(V*1.0e3)* Qt, vrange=Vrange)
Na_m2h0_m2h1 = smodel.VDepSReac('Na_m2h0_m2h1', ssys, \
                                slhs=[Na_m2h0], srhs=[Na_m2h1], \
                                k=lambda V:1.0e3*a_h(V*1.0e3)* Qt, vrange=Vrange)
Na_m3h0_m3h1 = smodel.VDepSReac('Na_m3h0_m3h1', ssys, \
                                slhs=[Na_m3h0], srhs=[Na_m3h1], \
                                k=lambda V:1.0e3*a_h(V*1.0e3)* Qt, vrange=Vrange)

Na_m0h1_m0h0 = smodel.VDepSReac('Na_m0h1_m0h0', ssys, \
                                slhs=[Na_m0h1], srhs=[Na_m0h0], \
                                k=lambda V:1.0e3*b_h(V*1.0e3)* Qt, vrange=Vrange)
Na_m1h1_m1h0 = smodel.VDepSReac('Na_m1h1_m1h0', ssys, \
                                slhs=[Na_m1h1], srhs=[Na_m1h0], \
                                k=lambda V:1.0e3*b_h(V*1.0e3)* Qt, vrange=Vrange)
Na_m2h1_m2h0 = smodel.VDepSReac('Na_m2h1_m2h0', ssys, \
                                slhs=[Na_m2h1], srhs=[Na_m2h0], \
                                k=lambda V:1.0e3*b_h(V*1.0e3)* Qt, vrange=Vrange)
Na_m3h1_m3h0 = smodel.VDepSReac('Na_m3h1_m3h0', ssys, \
                                slhs=[Na_m3h1], srhs=[Na_m3h0], \
                                k=lambda V:1.0e3*b_h(V*1.0e3)* Qt, vrange=Vrange)

OC_Na = smodel.OhmicCurr('OC_Na', ssys, chanstate=Na_m3h1, g=Na_G, erev=Na_rev)

# Leak channel
L = smodel.Chan('L', mdl)
Leak = smodel.ChanState('Leak', mdl, L)

OC_L = smodel.OhmicCurr('OC_L', ssys, chanstate=Leak, g=L_G, erev=leak_rev)


##################################

########### MESH & COMPARTMENTALIZATION #################

##########Import Mesh

mesh = meshio.loadMesh('./meshes/'+meshfile_ab)[0]

outer_tets = range(mesh.ntets)

###USE OF gettets
#getcyl(tetmesh, rad,  zmin, zmax, binnum=120, x = 0.0, y = 0.0):
inner_tets = gettets.getcyl(mesh, 1e-6, -200e-6, 200e-6)[0]

for i in inner_tets: outer_tets.remove(i)
assert(outer_tets.__len__() + inner_tets.__len__() == mesh.ntets)

print outer_tets.__len__(), " tets in outer compartment"
print inner_tets.__len__(), " tets in inner compartment"

# Record voltage from the central tetrahedron
cent_tet = mesh.findTetByPoint([0.0,0.0,0.0])

########## Create an intracellular compartment i.e. cytosolic compartment

cyto = sgeom.TmComp('cyto', mesh, inner_tets)

if cyl160:
    # Ensure that we use points a small distance inside the boundary:
    LENGTH = mesh.getBoundMax()[2] - mesh.getBoundMin()[2]
    boundminz = mesh.getBoundMin()[2] + LENGTH/mesh.ntets
    boundmaxz = mesh.getBoundMax()[2] - LENGTH/mesh.ntets

    memb_tris = list(mesh.getSurfTris())
    minztris = []
    maxztris = []
    for tri in memb_tris:
        zminboundtri = True
        zmaxboundtri = True
        tritemp = mesh.getTri(tri)
        trizs = [0.0, 0.0, 0.0]
        trizs[0] = mesh.getVertex(tritemp[0])[2]
        trizs[1] = mesh.getVertex(tritemp[1])[2]
        trizs[2] = mesh.getVertex(tritemp[2])[2]
        for j in range(3):
            if (trizs[j]>boundminz): zminboundtri = False
        if (zminboundtri):
            minztris.append(tri)
            continue
        for j in range(3):
            if (trizs[j]< boundmaxz): zmaxboundtri = False
        if (zmaxboundtri):
            maxztris.append(tri)

    for t in minztris: memb_tris.remove(t)
    for t in maxztris: memb_tris.remove(t)
    
else:
    print 'Finding connecting triangles...'
    out_tris = set()
    for i in outer_tets:
            tritemp = mesh.getTetTriNeighb(i)
            for j in range(4): out_tris.add(tritemp[j])

    in_tris = set()
    for i in inner_tets:
            tritemp = mesh.getTetTriNeighb(i)
            for j in range(4): in_tris.add(tritemp[j])

    memb_tris = out_tris.intersection(in_tris)
    memb_tris = list(memb_tris)

print len(memb_tris), " surface triangles."

########## Create a membrane as a surface mesh
memb = sgeom.TmPatch('memb', mesh, memb_tris, cyto)
memb.addSurfsys('ssys')

print "Area: ", memb.getArea()

print "Creating membrane.."
membrane = sgeom.Memb('membrane', mesh, [memb])
print "Membrane created."

# # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

r = srng.create_mt19937(512)
r.initialize(7)

sim = ssolver.Tetexact(mdl, mesh, r, True)

print "Resetting simulation object.."
sim.reset()

print "Injecting molecules.."

sim.setTemp(TEMPERATURE+273.15)

surfarea = sim.getPatchArea('memb')

sim.setPatchCount('memb', 'Leak', round(L_ro * surfarea))

sim.setPatchCount('memb', 'Na_m0h0', round(Na_ro*surfarea*Na_facs[0]))
sim.setPatchCount('memb', 'Na_m1h0', round(Na_ro*surfarea*Na_facs[1]))
sim.setPatchCount('memb', 'Na_m2h0', round(Na_ro*surfarea*Na_facs[2]))
sim.setPatchCount('memb', 'Na_m3h0', round(Na_ro*surfarea*Na_facs[3]))
sim.setPatchCount('memb', 'Na_m0h1', round(Na_ro*surfarea*Na_facs[4]))
sim.setPatchCount('memb', 'Na_m1h1', round(Na_ro*surfarea*Na_facs[5]))
sim.setPatchCount('memb', 'Na_m2h1', round(Na_ro*surfarea*Na_facs[6]))
sim.setPatchCount('memb', 'Na_m3h1', round(Na_ro*surfarea*Na_facs[7]))

sim.setPatchCount('memb', 'K_n0', round(K_ro*surfarea*K_facs[0]))
sim.setPatchCount('memb', 'K_n1', round(K_ro*surfarea*K_facs[1]))
sim.setPatchCount('memb', 'K_n2', round(K_ro*surfarea*K_facs[2]))
sim.setPatchCount('memb', 'K_n3', round(K_ro*surfarea*K_facs[3]))
sim.setPatchCount('memb', 'K_n4', round(K_ro*surfarea*K_facs[4]))

print 'Leak', round(L_ro * surfarea)

print 'Na_m0h0', round(Na_ro*surfarea*Na_facs[0])
print 'Na_m1h0', round(Na_ro*surfarea*Na_facs[1])
print 'Na_m2h0', round(Na_ro*surfarea*Na_facs[2])
print 'Na_m3h0', round(Na_ro*surfarea*Na_facs[3])
print 'Na_m0h1', round(Na_ro*surfarea*Na_facs[4])
print 'Na_m1h1', round(Na_ro*surfarea*Na_facs[5])
print 'Na_m2h1', round(Na_ro*surfarea*Na_facs[6])
print 'Na_m3h1', round(Na_ro*surfarea*Na_facs[7])

print 'K_n0', round(K_ro*surfarea*K_facs[0])
print 'K_n1', round(K_ro*surfarea*K_facs[1])
print 'K_n2', round(K_ro*surfarea*K_facs[2])
print 'K_n3', round(K_ro*surfarea*K_facs[3])
print 'K_n4', round(K_ro*surfarea*K_facs[4])

print "Targeted Injection: ", round(Na_ro*surfarea), "Na channels"

print "Targeted Injection: ", round(K_ro*surfarea), "K channels"

print "Targeted Injection: ", round(L_ro*surfarea), "Leak channels"

sim.setEfieldDT(EF_DT)

sim.setMembPotential('membrane', init_pot)

sim.setMembVolRes('membrane', Ra)

sim.setMembCapac('membrane',memb_capac)


#### Recording #####

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

r.initialize(100*int(iter_n))

for l in range(NTIMEPOINTS):
    print "Tpnt: ", l

    sim.run(TIMECONVERTER*l)

    tcur_Na = 0.0
    tcur_K = 0.0
    tcur_L = 0.0

    for m in memb_tris:
        tcur_Na = tcur_Na + sim.getTriOhmicI(m,'OC_Na') 
        tcur_K = tcur_K + sim.getTriOhmicI(m,'OC_K')
        tcur_L = tcur_L + sim.getTriOhmicI(m,'OC_L')
    
    datfile.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile.write('%.6g' %((tcur_Na*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_K*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_L*1.0e-1)/surfarea) + ' ')
    datfile.write('\n')
    
    datfile2.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile2.write('%.6g' %(sim.getTetV(cent_tet)*1.0e3) + ' ')
    datfile2.write('\n')

## END
