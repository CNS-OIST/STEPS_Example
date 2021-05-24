# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Rallpack1 model
# Author Iain Hepburn

# Updated for STEPS 3.6.0
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from __future__ import print_function
import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.utilities.meshio as meshio
import steps.solver as ssolver

import math
from random import *
from pylab import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Print the voltage RMS of correct and STEPS solutions
def RMS(correct, steps):
    rms = 0    
    nts = len(correct)
    
    for t in range(nts):
        rms+= math.pow((correct[t]-steps[t]), 2)
    
    rms/=nts
    rms=math.sqrt(rms)
    
    print(rms)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

meshfile = 'axon_cyl_L1000um_D1000nm_38819tets'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Leak conductance, Siemens/m^2
L_G = 0.25

# Diameter
diam = 1.0

# The current injection in amps
Iinj = 1.0*0.1e-9

# Ohm.m
Ra = 1.00

# Total leak conductance for ideal cylinder:
surfarea_cyl = diam*math.pi*1000*1e-12

# Total leak conductance
L_G_tot = L_G*surfarea_cyl

# Leak reveral potential
leak_rev = -65.0e-3

# # # # # # # # # # # # # # # # SIMULATION CONTROLS # # # # # # # # # # # # # # 

# The EField calculation time step (seconds)
EF_DT = 1.0e-5	

# The simulation dt (seconds); must be larger than EField dt
SIM_DT = 5.0e-5

# Sim end time (seconds)
SIM_END = 0.25

# The number of sim 'time points'; * SIM_DT = sim end time
SIM_NTPNTS = int(SIM_END/SIM_DT)+1

# # # # # # # # # # # # # DATA COLLECTION # # # # # # # # # # # # # # # # # # 

# Length of the mesh, in m
LENGTH = 1000.0e-6

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

mdl = smodel.Model()
ssys = smodel.Surfsys('ssys', mdl)

L = smodel.Chan('L', mdl)
Leak = smodel.ChanState('Leak', mdl, L)

mesh = meshio.loadMesh('meshes/'+meshfile)[0]

cyto = sgeom.TmComp('cyto', mesh, range(mesh.ntets))


# Find the tets connected to the bottom face
# First find all the tets with ONE face on a boundary
boundtets = []
#store the 0to3 index of the surface triangle for each of these boundary tets
bt_srftriidx = []

for i in range(mesh.ntets):
	tettemp = mesh.getTetTetNeighb(i)
	if (tettemp[0] == sgeom.UNKNOWN_TET or tettemp[1] == sgeom.UNKNOWN_TET or tettemp[2] == sgeom.UNKNOWN_TET or tettemp[3] == sgeom.UNKNOWN_TET): 
		boundtets.append(i)
		templist = []
		if (tettemp[0] == sgeom.UNKNOWN_TET): 
			templist.append(0)
		if (tettemp[1] == sgeom.UNKNOWN_TET): 
			templist.append(1)
		if (tettemp[2] == sgeom.UNKNOWN_TET): 
			templist.append(2)
		if (tettemp[3] == sgeom.UNKNOWN_TET): 
			templist.append(3)
		bt_srftriidx.append(templist)

assert (boundtets.__len__() == bt_srftriidx.__len__())

# Find the tets on the z=0 and z=1000um boundaries, and the triangles
# Note: Storing the minZ tets now, but not yet using them 
minztets = []
minztris = []
maxztris = []
minzverts=set([])
maxzverts=set([])

boundminz = mesh.getBoundMin()[2] + LENGTH/mesh.ntets
boundmaxz = mesh.getBoundMax()[2] - LENGTH/mesh.ntets

for i in range(boundtets.__len__()):
    # get the boundary triangle
    for btriidx in bt_srftriidx[i]:
        zminboundtri = True
        tribidx = mesh.getTetTriNeighb(boundtets[i])[btriidx]
        tritemp = mesh.getTri(tribidx)
        trizs = [0.0, 0.0, 0.0]
        trizs[0] = mesh.getVertex(tritemp[0])[2]
        trizs[1] = mesh.getVertex(tritemp[1])[2]
        trizs[2] = mesh.getVertex(tritemp[2])[2]
        for j in range(3):
            if (trizs[j]>boundminz): zminboundtri = False
        if (zminboundtri): 
            minztets.append(boundtets[i])
            minztris.append(tribidx)    
            minzverts.add(tritemp[0])
            minzverts.add(tritemp[1])
            minzverts.add(tritemp[2])            
            continue
        
        zmaxboundtri = True
        for j in range(3):
            if (trizs[j]< boundmaxz): zmaxboundtri = False
        if (zmaxboundtri): 
            maxztris.append(tribidx)   
            maxzverts.add(tritemp[0])
            maxzverts.add(tritemp[1])
            maxzverts.add(tritemp[2])  

n_minztris = len(minztris)
assert(n_minztris > 0)

minzverts = list(minzverts)
maxzverts = list(maxzverts)

n_minzverts = len(minzverts)
assert(n_minzverts > 0)

min0dist=1000
min0distvert = -1
max0dist=-1
max0distvert = -1

for v in minzverts: 
    vert = mesh.getVertex(v)
    dist = math.sqrt((vert[0]*vert[0]) + (vert[1]*vert[1]))
    if dist < min0dist: 
        min0dist = dist
        min0distvert = v
    if dist > max0dist:
        max0dist = dist
        max0distvert = v

min1000dist=1000
min1000distvert = -1
max1000dist=-1
max1000distvert = -1

for v in maxzverts: 
    vert = mesh.getVertex(v)
    dist = math.sqrt((vert[0]*vert[0]) + (vert[1]*vert[1]))
    if dist < min1000dist: 
        min1000dist = dist
        min1000distvert = v
    if dist > max1000dist:
        max1000dist = dist
        max1000distvert = v

memb_tris = list(mesh.getSurfTris())

# Doing this now, so will inject into first little z section
for t in minztris: memb_tris.remove(t)
for t in maxztris: memb_tris.remove(t)


# Create the membrane with the tris removed at faces
memb = sgeom.TmPatch('memb', mesh, memb_tris, cyto)
memb.addSurfsys('ssys')

corr_fac_area = memb.getArea()/surfarea_cyl

membrane = sgeom.Memb('membrane', mesh, [memb] )

# Set the single-channel conductance:
g_leak_sc = L_G_tot/len(memb_tris)
OC_L = smodel.OhmicCurr('OC_L', ssys, chanstate = Leak, erev = leak_rev, g = g_leak_sc) 

# Create random number generator
r = srng.create('mt19937',512)
r.initialize(7)

# Create solver object
sim = ssolver.Tetexact(mdl, mesh, r, True)


surfarea_mesh = sim.getPatchArea('memb')
surfarea_cyl = 1.0*math.pi*1000*1e-12


vol_cyl = math.pi*(diam/2.0)*(diam/2.0)*1000*1e-18
vol_mesh = sim.getCompVol('cyto')

corr_fac_vol = vol_mesh/vol_cyl

RES_POT = zeros(( SIM_NTPNTS, 4))

print("Running simulation...")
sim.reset()

print("Injecting molecules..")
for t in memb_tris: sim.setTriCount(t, 'Leak', 1)

    
sim.setEfieldDT(EF_DT)
sim.setMembPotential('membrane', -65e-3)
sim.setMembVolRes('membrane', Ra)
sim.setMembCapac('membrane', 0.01/corr_fac_area)

for v in minzverts: sim.setVertIClamp(v, Iinj/n_minzverts)

for l in range(SIM_NTPNTS):
    if not l%100: print("Sim time (ms): ", SIM_DT*l*1.0e3)
        
    sim.run(SIM_DT*l)
    
    # Write at the beginning and end with vertex voltages
    RES_POT[l,0] = sim.getVertV(min0distvert)*1.0e3
    RES_POT[l,1] = sim.getVertV(max0distvert)*1.0e3
    RES_POT[l,2] = sim.getVertV(min1000distvert)*1.0e3
    RES_POT[l,3] = sim.getVertV(max1000distvert)*1.0e3

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Load the correct solutions

# At 0um- the end of the mesh
ifile_benchmark_x0 = open('data/rallpack1_correct/v0', 'r')
# At 1000um- the end of the mesh
ifile_benchmark_x1000 = open('data/rallpack1_correct/vx', 'r')

tpnt_benchmark = []

v_benchmark_x0 = []
v_benchmark_x1000 = []

lines_benchmark_x0 = ifile_benchmark_x0.readlines()

# Read in mv and ms
for line_benchmark_x0 in lines_benchmark_x0:
    nums = line_benchmark_x0.split()
    tpnt_benchmark.append(float(nums[0])*1e3)
    v_benchmark_x0.append(float(nums[1])*1e3)
    
lines_benchmark_x1000 = ifile_benchmark_x1000.readlines()

for line_benchmark_x1000 in lines_benchmark_x1000:
    nums = line_benchmark_x1000.split()
    v_benchmark_x1000.append(float(nums[1])*1e3)
    

print("Voltage RMS at 0um (mV):",)
RMS(v_benchmark_x0, mean( (RES_POT[:,0], RES_POT[:,1]), axis=0))
print("Voltage RMS at 1000um (mV):",)
RMS(v_benchmark_x1000, mean( (RES_POT[:,3], RES_POT[:,2]), axis=0))

TPNTS = arange(0.0, SIM_NTPNTS*SIM_DT*1.0e3, SIM_DT*1.0e3)
TPNTS.resize(SIM_NTPNTS)

subplot(211)
plot(tpnt_benchmark, v_benchmark_x0, 'k-' ,label = 'Correct, 0um', linewidth=3)
plot(TPNTS, RES_POT[:,0],'r--', label = 'STEPS, 0um', linewidth=3)
legend(loc='best')
ylabel('Potential (mV)')
subplot(212)
plot(tpnt_benchmark, v_benchmark_x1000, 'k-' ,label = 'Correct, 1000um', linewidth=3)
plot(TPNTS, RES_POT[:,2],'r--', label = 'STEPS, 1000um', linewidth=3)
legend(loc='best')
ylabel('Potential (mV)')
xlabel('Time (ms)')
show()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
