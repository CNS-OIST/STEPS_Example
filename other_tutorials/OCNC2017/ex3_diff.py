# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Example 3: Unbounded diffusion 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

from __future__ import print_function
import math
import time

from pylab import *

import steps.model as smodel
import steps.solver as solvmod
import steps.utilities.meshio as meshio 
import steps.geom as stetmesh
import steps.rng as srng

########################################################################

# Number of iterations; plotting dt; sim endtime:
NITER = 10
DT = 0.001
INT = 0.2

# The Abaqus file to import
MESHFILE = 'sp_10r_3875'

# The diffusion constant
DCST = 0.02e-9

# Number of molecules injected in centre 
NINJECT = 1000	

# The number of tetrahedrons to sample, selected at random
SAMPLE = 500	

# Create the array of tet indices to be found at random
tetidxs = zeros(SAMPLE, dtype = 'int')

# Further create the array of tet barycentre distance to centre
tetrads = zeros(SAMPLE)

# Make random number generator avaialble for initial conditions
rng = srng.create('mt19937', 256) 
rng.initialize(123) 

########################################################################

beg_time = time.time()

# Function to plot simulation run time
def printtime(end_time):

	totsecs = int(end_time-beg_time)
	sec = totsecs%60
	totmin = totsecs/60
	min = totmin%60
	hours = totmin/60
	
	print('Simulation time: %d h, %d min, %d sec' %(hours, min, sec))
	
########################################################################

def gen_model():
   
    # Create model container object 
    mdl = smodel.Model()
    
    # Create species 'X'
    X = smodel.Spec('X', mdl)
    
    # Create volume system               
    cytosolv = smodel.Volsys('cytosolv', mdl)
    
    # Create diffusion rule
    dif_X = smodel.Diff('diffX', cytosolv, X,  DCST)

    return mdl
	
########################################################################

def gen_geom(meshfile):
	print("Loading mesh...")
	mesh = meshio.importAbaqus('meshes/'+meshfile, 1.0e-6)[0]
	print("Mesh loaded.")
	
	ntets = mesh.countTets()
    
    # Create cytosol compartment
	cyto = stetmesh.TmComp('cyto', mesh, range(ntets))
	cyto.addVolsys('cytosolv')
	
	# Now fill the array holding the tet indices to sample at random
	if (SAMPLE > ntets):
		print("SAMPLE larger than total number of tetrahedrons.")
		return
	if (SAMPLE < 5):
		print("SAMPLE must be greater than 4.")
	
	# First add the centre tet and its 4 neighbours
	ctetidx = mesh.findTetByPoint([0.0, 0.0, 0.0])
	tetidxs[0] = ctetidx
	neighbs = mesh.getTetTetNeighb(ctetidx)
	tetidxs[1] = neighbs[0]
	tetidxs[2] = neighbs[1]
	tetidxs[3] = neighbs[2]
	tetidxs[4] = neighbs[3]
	
	filled = 5
	
    # Find the rest, considering spherical geoemetry
	while (filled < SAMPLE):
		max = mesh.getBoundMax()
		min = mesh.getBoundMin()
		
		maxX2 = 0.0
		maxY2 = 0.0
		maxZ2 = 0.0
		
		if (max[0] > -min[0]) : maxX2 = abs(math.pow(max[0], 1.0/2))
		else : maxX2 = abs(math.pow(abs(min[0]), 1.0/2))
		if (max[1] > -min[1]) : maxY2 = abs(math.pow(max[1], 1.0/2))
		else : maxY2 = abs(math.pow(abs(min[1]), 1.0/2))
		if (max[2] > -min[2]) : maxZ2 = abs(math.pow(max[2], 1.0/2))
		else : maxZ2 = abs(math.pow(abs(min[2]), 1.0/2))
		
		rnx = rng.getUnfII()
		rny = rng.getUnfII()
		rnz = rng.getUnfII()
		
		signx = rng.getUnfII()
		signy = rng.getUnfII()
		signz = rng.getUnfII()
		
		if (signx >= 0.5) : xpnt = math.pow((maxX2*rnx), 2)
		else : xpnt = -(math.pow((maxX2*rnx), 2))
		
		if (signy >= 0.5) : ypnt = math.pow((maxY2*rny), 2)
		else : ypnt = -(math.pow((maxY2*rny), 2))
		
		if (signz >= 0.5) : zpnt = math.pow((maxZ2*rnz), 2)
		else : zpnt = -(math.pow((maxZ2*rnz), 2))
		
		idx = mesh.findTetByPoint([xpnt, ypnt, zpnt])
		
		if (idx == -1): continue
		if (idx not in tetidxs): 
			tetidxs[filled] = idx
			filled += 1
	
	tetidxs.sort()
	
	# Now find the distances to the centre of the centre tet (at 0,0,0)
	cbaryc = mesh.getTetBarycenter(ctetidx)
	for i in range(SAMPLE):
		baryc = mesh.getTetBarycenter(int(tetidxs[i]))
		r2 = math.pow((baryc[0]-cbaryc[0]),2) + \
			math.pow((baryc[1]-cbaryc[1]),2) + \
			math.pow((baryc[2]-cbaryc[2]),2)
		r = math.sqrt(r2)
		# Convert to microns
		tetrads[i] = r*1.0e6
	
	return mesh

########################################################################

m = gen_model()
g = gen_geom(MESHFILE)

# Fetch the index of the centre tet
ctetidx = g.findTetByPoint([0.0, 0.0, 0.0])

# And fetch the total number of tets to make the data structures
ntets = g.countTets()

# Create solver object
sim = solvmod.Tetexact(m, g, rng)

tpnts = arange(0.0, INT, DT)
ntpnts = tpnts.shape[0]

# Create the data structure to record concentrations
res = zeros((NITER, ntpnts, SAMPLE))

# Run the simulation and record data
for j in range(NITER):
    sim.reset()
    
    # Inject molecules into central tetrahedron
    sim.setTetCount(ctetidx, 'X', NINJECT)
    for i in range(ntpnts):
        sim.run(tpnts[i])
        
        for k in range(SAMPLE):
            # Record concentration from individual tetrahedrons
            res[j, i, k] = sim.getTetConc(int(tetidxs[k]), 'X')*1.0e6
    print('%d / %d' % (j + 1, NITER))
    printtime(time.time())

itermeans = mean(res, axis = 0)

########################################################################

# Plotting function for recorded concentrations
def plotconc(tidx):
	if (tidx >= INT/DT):
		print("Time index is out of range.")
		return
	
	scatter(tetrads, itermeans[tidx], s=2)
	xlim=(0.0, 10.0)
	ylim=(0.0)
	xlabel('Radial distance ($\mu$m)')
	ylabel('Mean concentration ($\mu$M)')
	t = tpnts[tidx]
	title ('Time: '+str(t)+'s')
	_plotdetc(t)
	show()

########################################################################

# Plotting function for analytical diffusion for comparison
def _plotdetc(timepnt):
	npnts = 200
	detconc = zeros((npnts))
	radialds = zeros((npnts))
	maxrad = 0.0
	for i in tetrads:
		if (i > maxrad): maxrad = i
	maxrad *= 1e-6
	intervals = maxrad/npnts
	rad = 0.0
	for i in range((npnts)):
		# Find the analytical concentration and convert to mol/L
		detconc[i] = 1.0e3*(1/6.022e23) * \
			((NINJECT/(math.pow((4*math.pi*DCST*timepnt),1.5)))*\
			(math.exp((-1.0*(rad*rad))/(4*DCST*timepnt))))
		radialds[i] = rad*1e6
		rad += intervals
	plot(radialds, detconc, color = 'red')
	
########################################################################

print("Number of time points: ", int(INT/DT))

########################################################################
