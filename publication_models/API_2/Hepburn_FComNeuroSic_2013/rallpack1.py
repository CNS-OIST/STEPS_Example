# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Rallpack1 model
# Author Iain Hepburn

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

import math
from random import *
from pylab import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Print the voltage RMS of correct and STEPS solutions
def RMS(correct, steps):
    rms = 0    
    nts = len(correct)
    
    for t in range(nts):
        rms += math.pow((correct[t]-steps[t]), 2)
    
    rms /= nts
    rms = math.sqrt(rms)
    
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

mdl = Model()
r = ReactionManager()
with mdl:
    ssys = SurfaceSystem.Create()

mesh = TetMesh.Load('meshes/'+meshfile)
with mesh:
    
    cyto = TetComp.Create(mesh.tets)
    
    minz, maxz = mesh.bbox.min.z, mesh.bbox.max.z
    memb_tris = TriList(tri for tri in mesh.surface if minz < tri.center.z < maxz)

    sides = mesh.surface - memb_tris
    minzverts = VertList(vert for vert in sides.verts if vert.z <= minz)
    maxzverts = VertList(vert for vert in sides.verts if vert.z >= maxz)

    min0distvert = min(minzverts, key=lambda vert: np.linalg.norm(vert[:2]))
    max0distvert = max(minzverts, key=lambda vert: np.linalg.norm(vert[:2]))

    min1000distvert = min(maxzverts, key=lambda vert: np.linalg.norm(vert[:2]))
    max1000distvert = max(maxzverts, key=lambda vert: np.linalg.norm(vert[:2]))

    # Create the membrane with the tris removed at faces
    memb = TetPatch.Create(memb_tris, cyto, None, ssys)
    
    corr_fac_area = memb.Area/surfarea_cyl
    
    membrane = Membrane.Create([memb])

# Set the single-channel conductance:
with mdl:
    Leak = SubUnitState.Create()
    L = Channel.Create([Leak])

    g_leak_sc = L_G_tot/len(memb_tris)
    with ssys:
        OC_L = OhmicCurr.Create(L[Leak], g_leak_sc, leak_rev) 

# Create random number generator
rng = RNG('mt19937', 512, 7)

# Create solver object
sim = Simulation('Tetexact', mdl, mesh, rng, calcMembPot=True)


surfarea_mesh = sim.memb.Area
surfarea_cyl = 1.0*math.pi*1000*1e-12


vol_cyl = math.pi*(diam/2.0)*(diam/2.0)*1000*1e-18
vol_mesh = sim.cyto.Vol

corr_fac_vol = vol_mesh/vol_cyl

# Data saving 

rs = ResultSelector(sim)

recordVerts = [min0distvert, max0distvert, min1000distvert, max1000distvert]
Vrs = rs.VERTS(recordVerts).V

sim.toSave(Vrs, dt=SIM_DT)

print("Running simulation...")
sim.newRun()

print("Injecting molecules..")
sim.TRIS(memb_tris).L[Leak].Count = 1

    
sim.EfieldDT = EF_DT
sim.membrane.Potential = -65e-3
sim.membrane.VolRes = Ra
sim.membrane.Capac = 0.01/corr_fac_area

sim.VERTS(minzverts).IClamp = Iinj / len(minzverts)

for l in range(SIM_NTPNTS):
    if l % 100 == 0:
        print("Sim time (ms): ", SIM_DT*l*1.0e3)
        
    sim.run(SIM_DT*l)

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


RES_POT = Vrs.data[0] * 1e3

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
