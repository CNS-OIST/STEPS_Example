# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Rallpack3 model
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
import time
from random import *
from pylab import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def stats(bench, cdata):
    peaktimes_bench=[]
    peaktimes_cdata=[]
    
    prev_v = -1000
    climbing=True
    
    nts = len(bench)
    for t in range(nts):
        if climbing == True: 
            if bench[t] < prev_v and bench[t] > 0:
                peaktimes_bench.append((t-1)*0.005)
                climbing = False
        else: 
            if bench[t] > prev_v:
                climbing = True
        prev_v = bench[t]
        
    prev_v = -1000
    climbing=True
        
    nts = len(cdata)
    for t in range(nts):
        if climbing == True: 
            if cdata[t] < prev_v and cdata[t] > 0:
                peaktimes_cdata.append((t-1)*0.005)
                climbing = False
        else: 
            if cdata[t] > prev_v:
                climbing = True   
        prev_v = cdata[t]    

    time_diff = 0
    
    nps = min([len(peaktimes_bench), len(peaktimes_cdata)])
    
    for p in range(nps):
        time_diff+=abs(peaktimes_bench[p]-peaktimes_cdata[p])
    
    time_diff/=nps    

    print("Number of peaks", nps)
    print("Mean absolute peak time difference:", time_diff, 'ms')
    
    
    rms = 0
        
    if len(bench) != len(cdata):
        print("Warning: data different length tpnts", len(bench), len(cdata))

    nts = len(bench)
    for t in range(nts):
        rms+= math.pow((bench[t]-cdata[t]), 2)
    
    rms/=nts
    rms=math.sqrt(rms)
    
    print("Full root mean square:", rms, 'mV')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

meshfile ='axon_cube_L1000um_D866m_1135tets'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Potassium conductance, Siemens/m^2
K_G = 360
# Sodium conductance, Siemens/m^2
Na_G = 1200
# Leak conductance, Siemens/m^2
L_G = 0.25

# Potassium reversal potential, V
K_rev = -77e-3
# Sodium reversal potential, V
Na_rev = 50e-3
# Leak reveral potential, V
leak_rev = -65.0e-3

# Potassium channel density
K_ro = 18.0e12
# Sodium channel density
Na_ro = 60.0e12


# Total leak conductance for ideal cylinder:
surfarea_cyl = 1.0*math.pi*1000*1e-12
L_G_tot = L_G*surfarea_cyl


# A table of potassium density factors at -65mV, found in getpops. n0, n1, n2, n3, n4
K_FACS = [0.216750577045, 0.40366011853, 0.281904943772, \
            0.0874997924409, 0.0101845682113 ]

# A table of sodium density factors. m0h1, m1h1, m2h1, m3h1, m0h0, m1h0, m2h0, m3h0
NA_FACS = [0.343079175644, 0.0575250437508, 0.00321512825945, 5.98988373918e-05, \
            0.506380603793, 0.0849062503811, 0.00474548939393, 8.84099403236e-05]


# Ohm.m
Ra = 1.0

# # # # # # # # # # # # # # # # SIMULATION CONTROLS # # # # # # # # # # # # # # 

# The simulation dt (seconds); for TetODE this is equivalent to EField dt
SIM_DT = 5.0e-6

# Sim end time (seconds)
SIM_END = 0.25

# The number of sim 'time points'; * SIM_DT = sim end time
SIM_NTPNTS = int(SIM_END/SIM_DT)+1

# The current injection in amps
Iinj = 0.1e-9

# # # # # # # # # # # # # DATA COLLECTION # # # # # # # # # # # # # # # # # # 

# record potential at the two extremes along (z) axis 
POT_POS = array([ 0.0, 1.0e-03])

POT_N = len(POT_POS)

# Length of the mesh, in m
LENGTH = 1000.0e-6

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

mdl = smodel.Model()
ssys = smodel.Surfsys('ssys', mdl)

# K channel
K = smodel.Chan('K', mdl)
K_n0 = smodel.ChanState('K_n0', mdl, K)		
K_n1 = smodel.ChanState('K_n1', mdl, K)
K_n2 = smodel.ChanState('K_n2', mdl, K)
K_n3 = smodel.ChanState('K_n3', mdl, K)
K_n4 = smodel.ChanState('K_n4', mdl, K)

# Na channel
Na = smodel.Chan('Na', mdl)
Na_m0h1 = smodel.ChanState('Na_m0h1', mdl, Na)
Na_m1h1 = smodel.ChanState('Na_m1h1', mdl, Na)
Na_m2h1 = smodel.ChanState('Na_m2h1', mdl, Na)
Na_m3h1 = smodel.ChanState('Na_m3h1', mdl, Na)
Na_m0h0 = smodel.ChanState('Na_m0h0', mdl, Na)
Na_m1h0 = smodel.ChanState('Na_m1h0', mdl, Na)
Na_m2h0 = smodel.ChanState('Na_m2h0', mdl, Na)
Na_m3h0 = smodel.ChanState('Na_m3h0', mdl, Na)

# Leak
L = smodel.Chan('L', mdl)
Leak = smodel.ChanState('Leak', mdl, L)

# Gating kinetics 
_a_m = lambda mV: ((((0.1 * (25 -(mV + 65.)) / (math.exp((25 - (mV + 65.)) / 10.) - 1)))))
_b_m = lambda mV: ((((4. * math.exp(-((mV + 65.) / 18.))))))
_a_h = lambda mV: ((((0.07 * math.exp((-(mV + 65.) / 20.))))))
_b_h = lambda mV: ((((1. / (math.exp((30 - (mV + 65.)) / 10.) + 1)))))
_a_n = lambda mV: ((((0.01 * (10 -(mV + 65.)) / (math.exp((10 - (mV + 65.)) / 10.) - 1)))))
_b_n = lambda mV: ((((0.125 * math.exp(-(mV + 65.) / 80.)))))

Kn0n1  = smodel.VDepSReac('Kn0n1', ssys, slhs = [K_n0], srhs = [K_n1], k = lambda V: 1.0e3 *4.*_a_n(V*1.0e3))  
Kn1n2  = smodel.VDepSReac('Kn1n2', ssys, slhs = [K_n1], srhs = [K_n2], k = lambda V: 1.0e3 *3.*_a_n(V*1.0e3))  
Kn2n3  = smodel.VDepSReac('Kn2n3', ssys, slhs = [K_n2], srhs = [K_n3], k = lambda V: 1.0e3 *2.*_a_n(V*1.0e3)) 
Kn3n4  = smodel.VDepSReac('Kn3n4', ssys, slhs = [K_n3], srhs = [K_n4], k = lambda V: 1.0e3 *1.*_a_n(V*1.0e3)) 

Kn4n3  = smodel.VDepSReac('Kn4n3', ssys, slhs = [K_n4], srhs = [K_n3], k = lambda V: 1.0e3 *4.*_b_n(V*1.0e3)) 
Kn3n2  = smodel.VDepSReac('Kn3n2', ssys, slhs = [K_n3], srhs = [K_n2], k = lambda V: 1.0e3 *3.*_b_n(V*1.0e3)) 
Kn2n1  = smodel.VDepSReac('Kn2n1', ssys, slhs = [K_n2], srhs = [K_n1], k = lambda V: 1.0e3 *2.*_b_n(V*1.0e3)) 
Kn1n0  = smodel.VDepSReac('Kn1n0', ssys, slhs = [K_n1], srhs = [K_n0], k = lambda V: 1.0e3 *1.*_b_n(V*1.0e3))  

Na_m0h1_m1h1 = smodel.VDepSReac('Na_m0h1_m1h1', ssys, slhs=[Na_m0h1], srhs=[Na_m1h1], k = lambda V: 1.0e3 * 3. * _a_m(V * 1.0e3))
Na_m1h1_m2h1 = smodel.VDepSReac('Na_m1h1_m2h1', ssys, slhs=[Na_m1h1], srhs=[Na_m2h1], k = lambda V: 1.0e3 * 2. * _a_m(V * 1.0e3))
Na_m2h1_m3h1 = smodel.VDepSReac('Na_m2h1_m3h1', ssys, slhs=[Na_m2h1], srhs=[Na_m3h1], k = lambda V: 1.0e3 * 1. * _a_m(V * 1.0e3))

Na_m3h1_m2h1 = smodel.VDepSReac('Na_m3h1_m2h1', ssys, slhs=[Na_m3h1], srhs=[Na_m2h1], k = lambda V: 1.0e3 * 3. * _b_m(V * 1.0e3))
Na_m2h1_m1h1 = smodel.VDepSReac('Na_m2h1_m1h1', ssys, slhs=[Na_m2h1], srhs=[Na_m1h1], k = lambda V: 1.0e3 * 2. * _b_m(V * 1.0e3))
Na_m1h1_m0h1 = smodel.VDepSReac('Na_m1h1_m0h1', ssys, slhs=[Na_m1h1], srhs=[Na_m0h1], k = lambda V: 1.0e3 * 1. * _b_m(V * 1.0e3))

Na_m0h0_m1h0 = smodel.VDepSReac('Na_m0h0_m1h0', ssys, slhs=[Na_m0h0], srhs=[Na_m1h0], k = lambda V: 1.0e3 * 3. * _a_m(V * 1.0e3))
Na_m1h0_m2h0 = smodel.VDepSReac('Na_m1h0_m2h0', ssys, slhs=[Na_m1h0], srhs=[Na_m2h0], k = lambda V: 1.0e3 * 2. * _a_m(V * 1.0e3))
Na_m2h0_m3h0 = smodel.VDepSReac('Na_m2h0_m3h0', ssys, slhs=[Na_m2h0], srhs=[Na_m3h0], k = lambda V: 1.0e3 * 1. * _a_m(V * 1.0e3))
    
Na_m3h0_m2h0 = smodel.VDepSReac('Na_m3h0_m2h0', ssys, slhs=[Na_m3h0], srhs=[Na_m2h0], k = lambda V: 1.0e3 * 3. * _b_m(V * 1.0e3))
Na_m2h0_m1h0 = smodel.VDepSReac('Na_m2h0_m1h0', ssys, slhs=[Na_m2h0], srhs=[Na_m1h0], k = lambda V: 1.0e3 * 2. * _b_m(V * 1.0e3))
Na_m1h0_m0h0 = smodel.VDepSReac('Na_m1h0_m0h0', ssys, slhs=[Na_m1h0], srhs=[Na_m0h0], k = lambda V: 1.0e3 * 1. * _b_m(V * 1.0e3))

Na_m0h1_m0h0 = smodel.VDepSReac('Na_m0h1_m0h0', ssys, slhs=[Na_m0h1], srhs=[Na_m0h0], k = lambda V: 1.0e3 *_a_h(V * 1.0e3))
Na_m1h1_m1h0 = smodel.VDepSReac('Na_m1h1_m1h0', ssys, slhs=[Na_m1h1], srhs=[Na_m1h0], k = lambda V: 1.0e3 *_a_h(V * 1.0e3))
Na_m2h1_m2h0 = smodel.VDepSReac('Na_m2h1_m2h0', ssys, slhs=[Na_m2h1], srhs=[Na_m2h0], k = lambda V: 1.0e3 *_a_h(V * 1.0e3))
Na_m3h1_m3h0 = smodel.VDepSReac('Na_m3h1_m3h0', ssys, slhs=[Na_m3h1], srhs=[Na_m3h0], k = lambda V: 1.0e3 *_a_h(V * 1.0e3))
    
Na_m0h0_m0h1 = smodel.VDepSReac('Na_m0h0_m0h1', ssys, slhs=[Na_m0h0], srhs=[Na_m0h1], k = lambda V: 1.0e3 *_b_h(V * 1.0e3))
Na_m1h0_m1h1 = smodel.VDepSReac('Na_m1h0_m1h1', ssys, slhs=[Na_m1h0], srhs=[Na_m1h1], k = lambda V: 1.0e3 *_b_h(V * 1.0e3))
Na_m2h0_m2h1 = smodel.VDepSReac('Na_m2h0_m2h1', ssys, slhs=[Na_m2h0], srhs=[Na_m2h1], k = lambda V: 1.0e3 *_b_h(V * 1.0e3))
Na_m3h0_m3h1 = smodel.VDepSReac('Na_m3h0_m3h1', ssys, slhs=[Na_m3h0], srhs=[Na_m3h1], k = lambda V: 1.0e3 *_b_h(V * 1.0e3))


OC_K = smodel.OhmicCurr('OC_K', ssys, chanstate = K_n4, erev = K_rev, g = K_G/K_ro ) 
OC_Na = smodel.OhmicCurr('OC_Na', ssys, chanstate = Na_m3h0, erev = Na_rev, g = Na_G/Na_ro ) 

# Mesh geometry
mesh = meshio.loadMesh('./meshes/'+meshfile)[0]

cyto = sgeom.TmComp('cyto', mesh, range(mesh.ntets))

# The tetrahedrons from which to record potential
POT_TET = zeros(POT_N, dtype = 'uint')

i=0
for p in POT_POS:
    # Assuming axiz aligned with z-axis
    POT_TET[i] = mesh.findTetByPoint([0.0, 0.0, POT_POS[i]])
    i=i+1

# Find the tets connected to the bottom face
# First find all the tets with ONE face on a boundary
boundtets = []
#store the 0to3 index of the surface triangle for each of these boundary tets
bt_srftriidx = []

for i in range(mesh.ntets):
	tettemp = mesh.getTetTetNeighb(i)
	if (tettemp[0] ==sgeom.UNKNOWN_TET or tettemp[1] == sgeom.UNKNOWN_TET or tettemp[2] == sgeom.UNKNOWN_TET or tettemp[3] == sgeom.UNKNOWN_TET):
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
minztets = []
minztris = []
maxztris = []
minzverts=set([])

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

n_minztris = len(minztris)
assert(n_minztris > 0)
minzverts = list(minzverts)
n_minzverts = len(minzverts)
assert(n_minzverts > 0)

memb_tris = list(mesh.getSurfTris())

# Doing this now, so will inject into first little z section
for t in minztris: memb_tris.remove(t)
for t in maxztris: memb_tris.remove(t)

# Create the membrane with the tris removed at faces
memb = sgeom.TmPatch('memb', mesh, memb_tris, cyto)
memb.addSurfsys('ssys')

membrane = sgeom.Memb('membrane', mesh, [memb], opt_method=2, search_percent = 100.0)

# Set the single-channel conductance:
g_leak_sc = L_G_tot/len(memb_tris)
OC_L = smodel.OhmicCurr('OC_L', ssys, chanstate = Leak, erev = leak_rev, g = g_leak_sc) 


# Create the solver objects
sim = ssolver.TetODE(mdl, mesh, calcMembPot= True)  
sim.setTolerances(1.0e-6, 1e-6)


surfarea_mesh = sim.getPatchArea('memb')
surfarea_cyl = 1.0*math.pi*1000*1e-12
corr_fac_area = surfarea_mesh/surfarea_cyl

vol_cyl = math.pi*0.5*0.5*1000*1e-18
vol_mesh = sim.getCompVol('cyto')
corr_fac_vol = vol_mesh/vol_cyl

RES_POT = zeros(( SIM_NTPNTS, POT_N))

print("\nRunning simulation")

for t in memb_tris: sim.setTriCount(t, 'Leak', 1)

sim.setPatchCount('memb', 'Na_m0h1', (Na_ro*surfarea_cyl*NA_FACS[0]))
sim.setPatchCount('memb', 'Na_m1h1', (Na_ro*surfarea_cyl*NA_FACS[1]))
sim.setPatchCount('memb', 'Na_m2h1', (Na_ro*surfarea_cyl*NA_FACS[2]))
sim.setPatchCount('memb', 'Na_m3h1', (Na_ro*surfarea_cyl*NA_FACS[3]))
sim.setPatchCount('memb', 'Na_m0h0', (Na_ro*surfarea_cyl*NA_FACS[4]))
sim.setPatchCount('memb', 'Na_m1h0', (Na_ro*surfarea_cyl*NA_FACS[5]))
sim.setPatchCount('memb', 'Na_m2h0', (Na_ro*surfarea_cyl*NA_FACS[6]))
sim.setPatchCount('memb', 'Na_m3h0', (Na_ro*surfarea_cyl*NA_FACS[7]))
sim.setPatchCount('memb', 'K_n0', (K_ro*surfarea_cyl*K_FACS[0]))
sim.setPatchCount('memb', 'K_n1', (K_ro*surfarea_cyl*K_FACS[1]))			
sim.setPatchCount('memb', 'K_n2', (K_ro*surfarea_cyl*K_FACS[2]))			
sim.setPatchCount('memb', 'K_n3', (K_ro*surfarea_cyl*K_FACS[3]))
sim.setPatchCount('memb', 'K_n4', (K_ro*surfarea_cyl*K_FACS[4]))

sim.setMembPotential('membrane', -65e-3)
sim.setMembVolRes('membrane', Ra*corr_fac_vol)
sim.setMembCapac('membrane', 0.01/corr_fac_area)

for v in minzverts: sim.setVertIClamp(v, Iinj/n_minzverts)

for l in range(SIM_NTPNTS):
    if not l%200: print("Sim time (ms): ", SIM_DT*l*1.0e3)
    
    sim.run(SIM_DT*l)
    
    for p in range(POT_N):
        RES_POT[l,p] = sim.getTetV(int(POT_TET[p]))*1.0e3

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Benchmark
# At 0um- the end of the mesh
ifile_benchmark_x0 = open('data/rallpack3_benchmark/rallpack3_0_0.001dt_1000seg', 'r')

# At 1000um- the end of the mesh
ifile_benchmark_x1000 = open('data/rallpack3_benchmark/rallpack3_1000_0.001dt_1000seg', 'r')

tpnt_benchmark = []
v_benchmark_x0 = []
v_benchmark_x1000 = []

lines_benchmark_x0 = ifile_benchmark_x0.readlines()[2:]

# Read in mv and ms
for line_benchmark_x0 in lines_benchmark_x0:
    nums = line_benchmark_x0.split()
    tpnt_benchmark.append(float(nums[0]))
    v_benchmark_x0.append(float(nums[1]))
    
lines_benchmark_x1000 = ifile_benchmark_x1000.readlines()[2:]

for line_benchmark_x1000 in lines_benchmark_x1000:
    nums = line_benchmark_x1000.split()
    v_benchmark_x1000.append(float(nums[1]))

# Get rid of the last point which seems to be missing from STEPS
v_benchmark_x0 = v_benchmark_x0[:-1]
tpnt_benchmark = tpnt_benchmark[:-1]
v_benchmark_x1000 = v_benchmark_x1000[:-1]

print("\nComparison at 0um (ms):",)
stats(v_benchmark_x0, RES_POT[:,0])
print("\nComparison at 1000um (ms):",)
stats(v_benchmark_x1000, RES_POT[:,1])

TPNTS = arange(0.0, SIM_NTPNTS*SIM_DT*1.0e3, SIM_DT*1.0e3)
TPNTS.resize(SIM_NTPNTS)


subplot(211)
plot(tpnt_benchmark, v_benchmark_x0, 'k-' ,label = 'Benchmark, 0um', linewidth=3)
plot(TPNTS, RES_POT[:,0],'r--', label = 'STEPS, 0um', linewidth=3)
legend(loc='best')
ylabel('Potential (mV)')
subplot(212)
plot(tpnt_benchmark, v_benchmark_x1000, 'k-' ,label = 'Benchmark, 1000um', linewidth=3)
plot(TPNTS, RES_POT[:,1],'r--', label = 'STEPS, 1000um', linewidth=3)
legend(loc='best')
ylabel('Potential (mV)')
xlabel('Time (ms)')
show()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
