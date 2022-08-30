# Example: Stochastic Calcium Burst model with GHK currents / 11.6. Plotting the results
# http://steps.sourceforge.net/manual/API_2/STEPS_Tutorial_CaBurst.html#Plotting-the-results

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.saving import *

from matplotlib import pyplot as plt
import numpy as np

with HDF5Handler('Caburst') as hdf:
    Currents, Pot, CaConcs, BKstates =  hdf['CaBurstSim'].results

    # Membrane potential
    plt.figure(figsize=(10, 7))
    for r in range(len(Pot.time)):
        plt.plot(Pot.time[r], 1e3 * Pot.data[r,:,0], label=f'Run {r}')
    plt.legend()
    plt.xlabel('Time [s]')
    plt.ylabel('Membrane Potential [mV]')
    plt.show()
    
    # Currents
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    for i in range(len(Currents.data[0,0,:])):
        ax = axs[i//2][i%2]
        for r in range(len(Currents.time)):
            ax.plot(Currents.time[r], 1e12 * Currents.data[r,:,i], label=f'Run {r}')
            ax.set_title(Currents.labels[i].split('.')[-2])
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Current [nA]')
    plt.legend()
    
    # Calcium
    fig, axs = plt.subplots(1, 2, squeeze=True, sharey=True, figsize=(10, 4))
    for i in range(len(CaConcs.data[0,0,:])):
        for r in range(len(CaConcs.time)):
            axs[i].plot(CaConcs.time[r], 1e6 * CaConcs.data[r,:,i], label=f'Run {r}')
        axs[i].set_title(['Cytoplasm', 'Submembrane tetrahedrons'][i])
        axs[i].set_xlabel('Time [s]')
        if i == 0:
            axs[i].set_ylabel('Calcium concentration [uM]')
    plt.legend()
    plt.show()
    
    # BK channel states for run 1
    rind = 1
    plt.figure(figsize=(10, 7))
    totCount = np.sum(BKstates.data[rind,:,:], axis=1)
    for c in range(len(BKstates.data[rind, 0, :])):
        plt.plot(BKstates.time[rind], 100 * BKstates.data[rind,:,c] / totCount)
    plt.legend(BKstates.labels)
    plt.xlabel('Time [s]')
    plt.ylabel('BK channels states distribution [%]')
    plt.show()
