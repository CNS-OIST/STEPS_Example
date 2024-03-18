# Example: Parallel simulations
# http://steps.sourceforge.net/manual/API_2/STEPS_Tutorial_MPI.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.saving import *

from matplotlib import pyplot as plt
import numpy as np

def plotPotential(CellPot, tidx):
    plt.plot(
        np.array(CellPot.metaData['tetzpos']) * 1e6,
        CellPot.data[0, tidx, :] * 1e3,
        label=f'{CellPot.time[0, tidx]*1e3} ms'
    )

def plotCurrents(NaCurrs, KCurrs, tidx, nbins=100):
    for results, currName in zip([NaCurrs, KCurrs], ['Na', 'K']):
        data = results.data[0, tidx, :] * 1e12
        pos = np.array(results.metaData['trizpos']) * 1e6
        areas = np.array(results.metaData['triarea']) * 1e12
        bins = np.histogram_bin_edges(pos, nbins)
        dig = np.digitize(pos, bins)
        # Ignore empty bins
        with np.errstate(invalid='ignore'):
            meanData = np.bincount(dig, weights=data) / np.bincount(dig, weights=areas)
            meanPos  = np.bincount(dig, weights=pos) / np.bincount(dig)
        plt.plot(meanPos, meanData, label=f'{currName} {results.time[0, tidx]*1e3} ms')

# # # # # # # # # # # # # # # # # #

with HDF5Handler('Efield_MPI') as hdf:
    NaCurrs, KCurrs, CellPot = hdf['TetOpSplitSim'].results

    plt.figure(figsize=(10, 7))
    plotPotential(CellPot, 10)
    plotPotential(CellPot, 20)
    plotPotential(CellPot, 30)
    plt.xlabel('Z-axis (um)')
    plt.ylabel('Membrane potential (mV)')
    plt.legend()
    plt.show()
    
    # # # # # # # # # # # # # # # # # #
    
    plt.figure(figsize=(10, 7))
    plotCurrents(NaCurrs, KCurrs, 10)
    plotCurrents(NaCurrs, KCurrs, 20)
    plotCurrents(NaCurrs, KCurrs, 30)
    plt.xlabel('Z-axis (um)')
    plt.ylabel('Current  (pA/um^2)')
    plt.legend()
    plt.show()
