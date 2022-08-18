# Example: Data recording and analysis
# http://steps.sourceforge.net/manual/API_2/STEPS_Tutorial_DataSaving.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *

import numpy as np

TF = 10
CF = 5e4
A = 1
B = 1

#########################
# Model, geom, Simulation
#########################

mdl = Model()
r = ReactionManager()
with mdl:
    X, Y = Species.Create()
    vsys = VolumeSystem.Create()
    with vsys:
        None <r['r1']> X >r['r2']> Y
        r['r1'].K = A * TF / CF, TF
        r['r2'].K = B * TF

        2*X + Y >r['r3']> 3*X
        r['r3'].K = TF * (CF ** 2)

geom = Geometry()
with geom:
    comp = Compartment.Create(vsys, vol=1e-19)

ENDT = 20
SIM_DT = 0.05

rng = RNG('mt19937', 512, 1234)

sim = Simulation('Wmrk4', mdl, geom, rng)
sim.setDT(5e-5)

rs = ResultSelector(sim)

concs = rs.comp.LIST(X, Y).Conc

sim.toSave(concs, dt=SIM_DT)

#########################
# Run simulation
#########################

# Simulation parameters
AVals = np.linspace(0, 2, 21).round(2)
BVals = np.linspace(0, 5, 21).round(2)

with HDF5Handler('Brusselator_wm') as hdf:
    for A in AVals:
        for B in BVals:
            sim.toDB(hdf, f'WM_A{A}_B{B}', A=A, B=B)

            sim.newRun()

            sim.comp.r1.K = (A * TF / CF), TF
            sim.comp.r2.K = B * TF
            sim.comp.r3.K = TF * (CF ** 2)

            sim.run(ENDT)

#########################
# Read data
#########################

from matplotlib import pyplot as plt

with HDF5Handler('Brusselator_wm') as hdf:
    group = hdf.get(A=0.3, B=3)

    concs, = group.results

    fig = plt.figure(figsize=(10,7))

    # Plot data for run 0
    plt.plot(concs.time[0], 1e3 * concs.data[0,:,:])

    plt.xlabel('Time [s]')
    plt.ylabel('Concentration [mM]')
    plt.legend(concs.labels)
    plt.show()

# # # # # # # # # # # # # # # # # #

pltAVals = [0.3, 0.5, 1, 1.5]
pltBVals = [4, 3, 2, 1]

with HDF5Handler('Brusselator_wm') as hdf:
    fig, axes = plt.subplots(4, 4, figsize=(10,10), sharey=True, sharex=True)

    # Plot data
    for i, B in enumerate(pltBVals):
        for j, A in enumerate(pltAVals):
            concs, = hdf.get(A=A, B=B).results
            axes[i][j].plot(concs.time[0], concs.data[0,:,:] * 1e3)

    # Add legend and labels
    axes[0][-1].legend(concs.labels)
    for k in range(4):
        axes[0][k].set_title(f'A = {pltAVals[k]}')
        axes[-1][k].set_xlabel('Time [s]', size='large')
        axes[k][0].set_ylabel(f'B = {pltBVals[k]}\nConcentr. [mM]', size='large')
    plt.show()

# # # # # # # # # # # # # # # # # #

tSkip = int(5 // SIM_DT)

with HDF5Handler('Brusselator_wm') as hdf:
    fig, axes = plt.subplots(1, len(pltAVals), figsize=(10, 2.5), sharey=True)

    for A, ax in zip(pltAVals, axes):
        BVals = []
        minMaxVals = []
        for group in hdf.filter(A=A):
            BVals.append(group.B)
            YConc = 1e3 * group.results[0].data[0,tSkip:,1]
            minMaxVals.append([np.max(YConc), np.min(YConc)])

        Bthresh = 1 + A ** 2
        ax.plot(BVals, minMaxVals, 'o')
        ax.axvline(Bthresh, color='r')
        ax.set_xlabel('B Parameter [AU]')
        ax.set_title(f'A = {A}')
    axes[0].set_ylabel('Y Concentration [mM]')
    plt.legend(['max Y', 'min Y', '$B_{thresh}$'], loc=2)

plt.show()

# # # # # # # # # # # # # # # # # #

from scipy.signal import find_peaks

with HDF5Handler('Brusselator_wm') as hdf:
    # Get all possible values of parameters A and B
    Avals = sorted(hdf.parameters['A'])
    Bvals = sorted(hdf.parameters['B'])

    # Compute amplitudes and frequencies of X peaks
    amplitudes = np.zeros((len(Avals), len(Bvals)))
    frequencies = np.zeros((len(Avals), len(Bvals)))
    for i, B in enumerate(Bvals):
        for j, A in enumerate(Avals):
            concs, = hdf.get(A=A, B=B).results
            time = concs.time[0,tSkip:]
            XConc = concs.data[0,tSkip:,0] * 1e3

            # Compute amplitude
            amplitudes[i,j] = np.max(XConc) - np.min(XConc)

            # Compute frequency
            peaks, _ = find_peaks(XConc, prominence=amplitudes[i,j] * 0.7)
            if len(peaks) > 1:
                frequencies[i,j] = 1 / np.mean(np.diff(time[peaks]))

fig = plt.figure(figsize=(10, 8))
plt.pcolormesh(Avals, Bvals, amplitudes)
plt.plot(Avals, [(1 + A ** 2) for A in Avals], '-r')
plt.plot(*np.meshgrid(pltAVals, pltBVals), '*', color='orange')
plt.xlabel('A Parameter [AU]')
plt.ylabel('B Parameter [AU]')
plt.title('Amplitude')
plt.colorbar(label='Amplitude [uM]')

fig = plt.figure(figsize=(10, 8))
plt.pcolormesh(Avals, Bvals, frequencies, cmap=plt.get_cmap('cividis'))
plt.plot(Avals, [(1 + A ** 2) for A in Avals], '-r')
plt.plot(*np.meshgrid(pltAVals, pltBVals), '*', color='orange')
plt.xlabel('A Parameter [AU]')
plt.ylabel('B Parameter [AU]')
plt.colorbar(label='Frequency [Hz]')
plt.title('Frequency')

plt.show()
