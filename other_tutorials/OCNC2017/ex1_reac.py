# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Example 1: Second-order reaction, well-mixed simulation 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Import biochemical model module
import steps.model as smod

# Create model container
mdl = smod.Model()

# Create chemical species
A = smod.Spec('A', mdl)				
B = smod.Spec('B', mdl)
C = smod.Spec('C', mdl)

# Create reaction set container
vsys = smod.Volsys('vsys', mdl)

# Create reaction
# A + B - > C with rate 200 /uM.s
reac_f = smod.Reac('reac_f', vsys, lhs=[A,B], rhs = [C])
reac_f.setKcst(200e6)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Import geometry module
import steps.geom as sgeom

# Create well-mixed geometry container
wmgeom = sgeom.Geom()

# Create cytosol compartment
cyt = sgeom.Comp('cyt', wmgeom)

# Give volume to cyt (m^3)
cyt.setVol(1.0e-18)

# Assign reaction set to compartment
cyt.addVolsys('vsys')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Import random number generator module
import steps.rng as srng

# Create random number generator, with buffer size as 256
r = srng.create('mt19937', 256)

# Initialise with some seed
r.initialize(899)

# Could use time to get random seed
#import time
#r.initialize(int(time.time()))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Import solver module
import steps.solver as ssolv

# Create well-mixed stochastic solver object
sim_direct = ssolv.Wmdirect(mdl, wmgeom, r)

# Import numpy 
import numpy as np

# Create time point Numpy array
tpnt = np.arange(0.0, 0.501, 0.001)

# Create data array, initialized with zeros
res_direct = np.zeros([int(0.501/0.001), 3])

# Reset the simulation (actually only important for multiple iterations)
sim_direct.reset()

# Initialize the number of 'A' molecules to 10
sim_direct.setCompCount('cyt', 'A', 10)       

# Or you can set the concentration (M), as for 'B'
sim_direct.setCompConc('cyt', 'B', 0.0332e-6)  

# Run simulation and record data
for t in range(0,int(0.501/0.001)):
	sim_direct.run(tpnt[t])
	res_direct[t,0] = sim_direct.getCompCount('cyt', 'A')
	res_direct[t,1] = sim_direct.getCompCount('cyt', 'B')		
	res_direct[t,2] = sim_direct.getCompCount('cyt', 'C')		

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Create well-mixed deterministic solver
sim_rk4 = ssolv.Wmrk4(mdl, wmgeom, r)

# Set the integration time-step (s)
sim_rk4.setRk4DT(0.00001)

# Repeat the simulation process for the deterministic solver
res_rk4 = np.zeros([int(0.501/0.001), 3])
sim_rk4.reset()
sim_rk4.setCompCount('cyt', 'A', 10)               
sim_rk4.setCompConc('cyt', 'B', 0.0332e-6)  
for t in range(0,int(0.501/0.001)):
	sim_rk4.run(tpnt[t])
	res_rk4[t,0] = sim_rk4.getCompCount('cyt', 'A')
	res_rk4[t,1] = sim_rk4.getCompCount('cyt', 'B')		
	res_rk4[t,2] = sim_rk4.getCompCount('cyt', 'C')
		
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Import matplotlib via pylab
from pylab import *

# Use Matlpotlib functions to plot data from both simulations
plot(tpnt, res_direct[:,0], label = 'A')
plot(tpnt, res_direct[:,1], label = 'B')
plot(tpnt, res_direct[:,2], label = 'C')

plot(tpnt, res_rk4[:,0], color='black')
plot(tpnt, res_rk4[:,1], color='black')
plot(tpnt, res_rk4[:,2], color='black')

ylabel('Number of molecules')
xlabel('Time (sec)')
legend()
show()

