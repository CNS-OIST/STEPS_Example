#################################################################################
#
#     STEPS - STochastic Engine for Pathway Simulation
#     Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#     Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
#     See the file AUTHORS for details.
#     This file is part of STEPS.
#
#     STEPS is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License version 2,
#     as published by the Free Software Foundation.
#
#     STEPS is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#  Example: Hodgkin-Huxley Action Potential propagation model
#  Author Iain Hepburn
#  http://steps.sourceforge.net/manual/memb_pot.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.saving import *

from matplotlib import pyplot as plt
import numpy as np

# Load results from files

NaCurrs = ResultSelector.FromFile('NaCurrs.dat')
KCurrs = ResultSelector.FromFile('KCurrs.dat')
CellPot = ResultSelector.FromFile('CellPot.dat')

NBINS = 100

def getTIdx(saver, t):
    return min(enumerate(saver.time[0]), key=lambda v: abs(v[1] - t))[0]

def plotPotential(t):
    tidx = getTIdx(CellPot, t)
    plt.plot(
        CellPot.metaData['tetzpos'] * 1e6, 
        CellPot.data[0, tidx, :] * 1e3, 
        label=f'{CellPot.time[0, tidx]*1e3} ms'
    )

def plotCurrents(t):
    tidx = getTIdx(NaCurrs, t)
    for results, currName in zip([NaCurrs, KCurrs], ['Na', 'K']):
        data = results.data[0, tidx, :] * 1e12
        pos = results.metaData['trizpos'] * 1e6
        areas = results.metaData['triarea'] * 1e12
        bins = np.histogram_bin_edges(pos, NBINS)
        dig = np.digitize(pos, bins)
        # Ignore empty bins
        with np.errstate(invalid='ignore'):
            meanData = np.bincount(dig, weights=data) / np.bincount(dig, weights=areas)
            meanPos  = np.bincount(dig, weights=pos) / np.bincount(dig)
        plt.plot(meanPos, meanData, label=f'{currName} {results.time[0, tidx]*1e3} ms')
    
plotPotential(1e-3)
plotPotential(2e-3)
plotPotential(3e-3)
plt.xlabel('Z-axis (um)')
plt.ylabel('Membrane potential (mV)')
plt.legend()
plt.show()
    
plotCurrents(1e-3)
plotCurrents(2e-3)
plotCurrents(3e-3)
plt.xlabel('Z-axis (um)')
plt.ylabel('Current  (pA/um^2)')
plt.legend()
plt.show()

