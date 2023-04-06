# Example: Surface-Volume Reactions / 4.6. Loading saved data
# http://steps.sourceforge.net/manual/API_2/STEPS_Tutorial_IP3.html#Loading-saved-data

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.saving import *
from matplotlib import pyplot as plt
import numpy as np

ldRstates = ResultSelector.FromFile('Rstates.dat')
ldReacs   = ResultSelector.FromFile('Reacs.dat')

plt.figure(figsize=(10, 7))

RopenInd = ldRstates.labels.index('Ropen')
RopenData = ldRstates.data[:, :, RopenInd]

time = ldRstates.time[0]
mean = np.mean(RopenData, axis=0)
std = np.std(RopenData, axis=0)

plt.plot(time, mean, linewidth=2, label='Average')
plt.fill_between(time, mean - std, mean + std, alpha=0.2, label='Std. Dev.')

for t, d in zip(ldRstates.time, RopenData):
    plt.plot(t, d, color='grey', linewidth=0.1, zorder=-1)

plt.ylim(0)
plt.margins(0, 0.05)
plt.xlabel('Time [s]')
plt.ylabel('Number of open IP3R')
plt.legend()
plt.show()

plt.figure(figsize=(10, 7))

time = ldRstates.time[0]
mean = np.mean(ldRstates.data, axis=0)
std = np.std(ldRstates.data, axis=0)

plt.plot(time, mean, linewidth=2)
for m, s in zip(mean.T, std.T):
    plt.fill_between(time, m - s, m + s, alpha=0.2)

plt.legend(ldRstates.labels)
plt.xlabel('Time [s]')
plt.ylabel('Number of receptors')
plt.ylim(0)
plt.margins(0, 0.05)
plt.show()

plt.figure(figsize=(10, 7))

time = ldReacs.time[0]
dt = time[1] - time[0]
meanDeriv = np.mean(np.gradient(ldReacs.data, dt, axis=1), axis=0)

plt.stackplot(time, meanDeriv.T)

plt.legend([f'd{l} / dt' for l in ldReacs.labels])
plt.margins(0, 0.05)
plt.xlabel('Time [s]')
plt.ylabel('Total reaction rate [1/s]')
plt.show()
