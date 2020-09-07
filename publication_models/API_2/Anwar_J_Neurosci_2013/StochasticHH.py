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

import steps.interface

import math
import time
from random import *
from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

from extra.constants_hh import *

import sys
import os

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

_, meshfile_ab, root, iter_n = sys.argv

if meshfile_ab == 'Cylinder2_dia2um_L160um_outer0_0.3shell_0.3size_279152tets_adaptive.inp':
    cyl160=True
else:
    cyl160=False

########################### BIOCHEMICAL MODEL ###############################

mdl = Model()
r = ReactionManager()

with mdl:
    
    ssys = SurfaceSystem.Create()

    # Potassium channel
    n0, n1, n2, n3, n4 = SubUnitState.Create()
    Kn = SubUnit.Create([n0, n1, n2, n3, n4])
    Kchan = Channel.Create([Kn])

    _a_n = VDepRate(lambda V: 1.0e3 * a_n(V*1.0e3)* Qt, vrange=Vrange)
    _b_n = VDepRate(lambda V: 1.0e3 * b_n(V*1.0e3)* Qt, vrange=Vrange)

    # Sodium channel
    m0, m1, m2, m3, h0, h1 = SubUnitState.Create()
    Nam, Nah = SubUnit.Create([m0, m1, m2, m3], [h0, h1])
    Nachan = Channel.Create([Nam, Nah])

    _a_m = VDepRate(lambda V:1.0e3*a_m(V*1.0e3)* Qt, vrange=Vrange)
    _b_m = VDepRate(lambda V:1.0e3*b_m(V*1.0e3)* Qt, vrange=Vrange)
    _a_h = VDepRate(lambda V:1.0e3*a_h(V*1.0e3)* Qt, vrange=Vrange)
    _b_h = VDepRate(lambda V:1.0e3*b_h(V*1.0e3)* Qt, vrange=Vrange)

    # Leak channel
    Leak = SubUnitState.Create()
    L = Channel.Create([Leak])

    with ssys:
        with Kchan[...]:
            n0.s <r[1]> n1.s <r[2]> n2.s <r[3]> n3.s <r[4]> n4.s
            r[1].setRates(4 * _a_n, 1 * _b_n)
            r[2].setRates(3 * _a_n, 2 * _b_n)
            r[3].setRates(2 * _a_n, 3 * _b_n)
            r[4].setRates(1 * _a_n, 4 * _b_n)

        with Nachan[...]:
            h0.s <r[1]> h1.s
            r[1].setRates(_a_h, _b_h)

            m0.s <r[1]> m1.s <r[2]> m2.s <r[3]> m3.s
            r[1].setRates(3*_a_m,   _b_m)
            r[2].setRates(2*_a_m, 2*_b_m)
            r[3].setRates(  _a_m, 3*_b_m)

        OC_K = OhmicCurr.Create(Kchan[n4], K_G, K_rev)
        OC_Na = OhmicCurr.Create(Nachan[m3, h1], Na_G, Na_rev)
        OC_L = OhmicCurr.Create(L[Leak], L_G, leak_rev)


##################################

########### MESH & COMPARTMENTALIZATION #################

##########Import Mesh
mesh = TetMesh.Load('./meshes/'+meshfile_ab)

with mesh:
    rad, zmin, zmax = 1e-6, -200e-6, 200e-6
    inner_tets, outer_tets = TetList(), TetList()
    for t in mesh.tets:
        c = t.center
        if zmin <= c.z <= zmax and c.x**2 + c.y**2 <= rad**2:
            inner_tets.append(t)
        else:
            outer_tets.append(t)

    print(len(outer_tets), " tets in outer compartment")
    print(len(inner_tets), " tets in inner compartment")

    # Record voltage from the central tetrahedron
    cent_tet = mesh.tets[0.0, 0.0, 0.0]

    ########## Create an intracellular compartment i.e. cytosolic compartment
    cyto = TetComp.Create(inner_tets)

    if cyl160:
        # Ensure that we use points a small distance inside the boundary:
        minz, maxz = mesh.bbox.min.z, mesh.bbox.max.z
        memb_tris = TriList(tri for tri in mesh_stock.surface if minz < tri.center.z < maxz)
    else:
        print('Finding connecting triangles...')
        memb_tris = inner_tets.surface & outer_tets.surface

    ########## Create a membrane as a surface mesh
    memb = TetPatch.Create(memb_tris, cyto, None, ssys)

    # For EField calculation
    print("Creating membrane..")
    membrane = Membrane.Create([memb])
    print("Membrane created.")

###### TRANSLATION TOKEN

# # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

rng = RNG('mt19937', 512, 7)

sim = Simulation('Tetexact', mdl, mesh, rng, calcMembPot=True)

#### Recording #####

dc = time.strftime('%b%d_%H_%M_%S_%Y')

runPath = os.path.join(root, 'data/StochasticHH/', meshfile_ab, f'{iter_n}__{dc}')
os.makedirs(runPath, exist_ok=True)

rs = ResultSelector(sim)

rs1 = rs.SUM(rs.TRIS(memb_tris).OC_Na.I) <<\
      rs.SUM(rs.TRIS(memb_tris).OC_K.I) <<\
      rs.SUM(rs.TRIS(memb_tris).OC_L.I)

rs2 = rs.TET(cent_tet).V

rs1.toFile(os.path.join(runPath, 'currents.dat.bin'))
rs2.toFile(os.path.join(runPath, 'voltage.dat.bin'))

sim.toSave(rs1, rs2, dt=TIMECONVERTER)

print("Resetting simulation object..")
sim.newRun()

print("Injecting molecules..")

sim.Temp = TEMPERATURE+273.15

surfarea = sim.memb.Area

sim.memb.L[Leak].Count = round(L_ro * surfarea)

for h, hsu in enumerate(Nah):
    for m, msu in enumerate(Nam):
        sim.memb.Nachan[msu, hsu].Count = round(Na_ro*surfarea*Na_facs[h*4 + m])

for n, ksu in enumerate(Kn):
    sim.memb.Kchan[ksu].Count = round(K_ro*surfarea*K_facs[n])

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

rng.initialize(100*int(iter_n))

for l in range(NTIMEPOINTS):
    print("Tpnt: ", l)

    sim.run(TIMECONVERTER*l)

# This last part is only present for backwards compatibility with the scripts created with API_1.
# We need to save to text files, like in the original script.

with open(os.path.join(runPath, 'currents.dat'), 'w') as f:
    for t, row in zip(rs1.time[0], rs1.data[0]):
        f.write('%.6g' % (t * 1e3) + ' ')
        for val in row:
            f.write('%.6g' % (val * 0.1 / surfarea) + ' ')
        f.write('\n')

with open(os.path.join(runPath, 'voltage.dat'), 'w') as f:
    for t, row in zip(rs2.time[0], rs2.data[0]):
        f.write('%.6g' % (t * 1e3) + ' ')
        for val in row:
            f.write('%.6g' % (val * 1e3) + ' ')
        f.write('\n')
