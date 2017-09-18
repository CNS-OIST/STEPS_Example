# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Example 2: Surface reactions: IP3 receptor model 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

from __future__ import print_function
import steps.rng as srng
import steps.solver as ssolver

from pylab import *

import ex2_ip3model as ip3r_model 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Simulation control variables
NITER = 100
T_END = 0.201
DT = 0.001
POINTS = int(T_END/DT)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Import model
mdl = ip3r_model.getModel()

# Import geometry 
geom = ip3r_model.getGeom()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Create random number generator
r = srng.create('mt19937', 512)
r.initialize(654)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Create well-mixed solver object
sim = ssolver.Wmdirect(mdl, geom, r)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Create Numpy array data structures
tpnt = arange(0.0, T_END, DT)
res_m = zeros([NITER, POINTS, 2])
res_std1 = zeros([POINTS, 2])
res_std2 = zeros([POINTS, 2])

# Run the simulation
for i in range (0, NITER):
    
    # Reset the simulation object
    sim.reset()
    
    # Set initial conditions
    sim.setCompConc('cyt', 'Ca', 3.30657e-8)               
    sim.setCompConc('cyt', 'IP3', 2.5e-6)   
    sim.setCompConc('ER', 'Ca', 150e-6)
    sim.setPatchCount('memb', 'R', 16)
    
    # Clamp ER calcium to fixed concentration
    sim.setCompClamped('ER', 'Ca', True)
        
    for t in range(0,POINTS):
        sim.run(tpnt[t])            
        
        # Record data
        res_m[i,t,0] = sim.getCompConc('cyt', 'Ca')*1e6		
        res_m[i,t,1] = sim.getPatchCount('memb', 'Ropen')

print("Ran ", NITER, "sim iterations ")

# Numpy array manipulation
mean_res = mean(res_m, 0)
std_res = std(res_m, 0)

def plotCa():
    # Matplotlib plotting functions
    plot(tpnt, mean_res[:,0], color = 'black', \
        linewidth = 1.0, label = 'mean')
    res_std1 = mean_res[:,0] + std_res[:,0]
    res_std2 = mean_res[:,0]- std_res[:,0]

    plot(tpnt, res_std1, color = 'gray', linewidth = 0.5, \
        label = 'standard deviation')
    plot(tpnt, res_std2,color = 'gray', linewidth = 0.5)

    xlabel('Time (sec)')
    ylabel('cytosolic Ca concentration ($\mu$M)')
    title('%d iterations' %NITER)
    ylim(0)
    legend()
    show()

def plotIP3R():
    # Matplotlib plotting functions
    plot(tpnt, mean_res[:,1], color = 'black', \
        linewidth = 1.0, label = 'mean')
    res_std1 = mean_res[:,1] + std_res[:,1]
    res_std2 = mean_res[:,1]- std_res[:,1]

    plot(tpnt, res_std1, color = 'gray', linewidth = 0.5, \
        label = 'standard deviation')
    plot(tpnt, res_std2,color = 'gray', linewidth = 0.5)

    xlabel('Time (sec)')
    ylabel('# open receptors')
    title('%d iterations' %NITER)
    ylim(0)
    legend()
    show()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
