####################################################################################
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
###

#  Example: SBML import
#  http://steps.sourceforge.net/manual/sbml_importer.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.sim import *
from steps.rng import *
from steps.saving import *

from matplotlib import pyplot as plt
import os

def runSBMLmod(solver):
    rng = RNG('mt19937', 256, 7233)

    dirPath = os.path.dirname(os.path.abspath(__file__))
    modelFile = os.path.join(dirPath, 'biomodel/BIOMD0000000098.xml')

    sim = SBMLSimulation(solver, modelFile, rng, 0.001, volunits_def=1e-3, volume_def=1e-18)

    rs = ResultSelector(sim)
    allComps = rs.ALL().ALL(Species).Conc

    sim.toSave(allComps, dt=0.01)

    sim.newRun()
    sim.run(10)

    plt.figure(figsize=(15, 10))
    plt.plot(allComps.time[0], allComps.data[0])
    plt.legend(allComps.labels)
    plt.xlabel('Time (s)')
    plt.ylabel('Conc (M)')
    plt.show()

runSBMLmod('Wmrk4')
runSBMLmod('Wmdirect')

