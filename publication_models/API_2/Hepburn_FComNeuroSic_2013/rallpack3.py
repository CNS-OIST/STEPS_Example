# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Rallpack3 model
# Author Iain Hepburn

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

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

# Length of the mesh, in m
LENGTH = 1000.0e-6

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

mdl = Model()
r = ReactionManager()
with mdl:
    ssys = SurfaceSystem.Create()

    # K channel
    n0, n1, n2, n3, n4 = SubUnitState.Create()
    VGKC = Channel.Create([n0, n1, n2, n3, n4])

    # Na channel
    h0, h1, m0, m1, m2, m3 = SubUnitState.Create()
    mNaSU, hNaSU = SubUnit.Create(
        [m0, m1, m2, m3],
        [h0, h1],
    )
    VGNaC = Channel.Create([mNaSU, hNaSU])

    # Leak
    leaksus = SubUnitState.Create()
    Leak = Channel.Create([leaksus])

    # Gating kinetics 
    _a_n = VDepRate(lambda V: 1e3*((0.01*(10-(V*1e3+65))/(math.exp((10-(V*1e3+65))/10)-1))))
    _b_n = VDepRate(lambda V: 1e3*((0.125*math.exp(-(V*1e3+65)/80))))
    _a_m = VDepRate(lambda V: 1e3*((0.1*(25-(V*1e3+65))/(math.exp((25-(V*1e3+65))/10)-1))))
    _b_m = VDepRate(lambda V: 1e3*((4*math.exp(-(V*1e3+65)/18))))
    _a_h = VDepRate(lambda V: 1e3*((0.07*math.exp(-(V*1e3+65)/20))))
    _b_h = VDepRate(lambda V: 1e3*((1/(math.exp((30-(V*1e3+65))/10)+1))))

    with ssys:
        with VGKC[...]:
            n0.s <r[1]> n1.s <r[2]> n2.s <r[3]> n3.s <r[4]> n4.s
            r[1].setRates(4*_a_n,   _b_n)
            r[2].setRates(3*_a_n, 2*_b_n)
            r[3].setRates(2*_a_n, 3*_b_n)
            r[4].setRates(  _a_n, 4*_b_n)

        with VGNaC[...]:
            h1.s <r[1]> h0.s
            r[1].setRates(_a_h, _b_h)

            m0.s <r[1]> m1.s <r[2]> m2.s <r[3]> m3.s
            r[1].setRates(3*_a_m,   _b_m)
            r[2].setRates(2*_a_m, 2*_b_m)
            r[3].setRates(  _a_m, 3*_b_m)

        VGKC_I = OhmicCurr.Create(VGKC[n4], K_G / K_ro, K_rev)
        VGNaC_I = OhmicCurr.Create(VGNaC[m3, h0], Na_G / Na_ro, Na_rev)

# Mesh geometry
mesh = TetMesh.Load('./meshes/'+meshfile)

with mesh:
    
    cyto = TetComp.Create(range(len(mesh.tets)))
    
    # The tetrahedrons from which to record potential
    POT_TET = TetList(mesh.tets[0, 0, z] for z in POT_POS)
    
    minz, maxz = mesh.bbox.min.z, mesh.bbox.max.z
    memb_tris = TriList(tri for tri in mesh.surface if minz < tri.center.z < maxz)

    sides = mesh.surface - memb_tris
    minzverts = VertList(vert for vert in sides.verts if vert.z <= minz)
    maxzverts = VertList(vert for vert in sides.verts if vert.z >= maxz)
    
    # Create the membrane with the tris removed at faces
    memb = TetPatch.Create(memb_tris, cyto, None, ssys)
    
    membrane = Membrane.Create([memb], opt_method=2, search_percent=100.0)

# Set the single-channel conductance:
g_leak_sc = L_G_tot/len(memb_tris)
with ssys:
    OC_L = OhmicCurr.Create(Leak[leaksus], g_leak_sc, leak_rev) 


# Create the solver objects
sim = Simulation('TetODE', mdl, mesh, None, calcMembPot=True)
sim.setTolerances(1.0e-6, 1e-6)


surfarea_mesh = sim.memb.Area
surfarea_cyl = 1.0*math.pi*1000*1e-12
corr_fac_area = surfarea_mesh/surfarea_cyl

vol_cyl = math.pi*0.5*0.5*1000*1e-18
vol_mesh = sim.cyto.Vol
corr_fac_vol = vol_mesh/vol_cyl

# Data saving

rs = ResultSelector(sim)

Vrs = rs.TETS(POT_TET).V

sim.toSave(Vrs, dt=SIM_DT)

sim.newRun()

print("\nRunning simulation")

sim.TRIS(memb_tris).Leak[leaksus].Count = 1

sim.memb.VGNaC[m0, h1].Count = (Na_ro*surfarea_cyl*NA_FACS[0])
sim.memb.VGNaC[m1, h1].Count = (Na_ro*surfarea_cyl*NA_FACS[1])
sim.memb.VGNaC[m2, h1].Count = (Na_ro*surfarea_cyl*NA_FACS[2])
sim.memb.VGNaC[m3, h1].Count = (Na_ro*surfarea_cyl*NA_FACS[3])
sim.memb.VGNaC[m0, h0].Count = (Na_ro*surfarea_cyl*NA_FACS[4])
sim.memb.VGNaC[m1, h0].Count = (Na_ro*surfarea_cyl*NA_FACS[5])
sim.memb.VGNaC[m2, h0].Count = (Na_ro*surfarea_cyl*NA_FACS[6])
sim.memb.VGNaC[m3, h0].Count = (Na_ro*surfarea_cyl*NA_FACS[7])
sim.memb.VGKC[n0].Count = (K_ro*surfarea_cyl*K_FACS[0])
sim.memb.VGKC[n1].Count = (K_ro*surfarea_cyl*K_FACS[1])
sim.memb.VGKC[n2].Count = (K_ro*surfarea_cyl*K_FACS[2])
sim.memb.VGKC[n3].Count = (K_ro*surfarea_cyl*K_FACS[3])
sim.memb.VGKC[n4].Count = (K_ro*surfarea_cyl*K_FACS[4])

sim.membrane.Potential = -65e-3
sim.membrane.VolRes = Ra*corr_fac_vol
sim.membrane.Capac = 0.01/corr_fac_area

sim.VERTS(minzverts).IClamp = Iinj / len(minzverts)

for l in range(SIM_NTPNTS):
    if not l%200:
        print("Sim time (ms): ", SIM_DT*l*1.0e3)
    
    sim.run(SIM_DT*l)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

RES_POT = Vrs.data[0] * 1e3

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
