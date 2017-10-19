####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   
###

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# Example: Hodgkin-Huxley Action Potential propagation model
# Author Iain Hepburn
# http://steps.sourceforge.net/manual/memb_pot.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

from __future__ import print_function
from pylab import *

from HH_APprop_tetode import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

tpnt = arange(0.0, N_timepoints*DT_sim, DT_sim)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

def plotVz(tidx):
    if (tidx >= tpnt.size): 
        print('Time index out of range')
        return
    plot(results[1]*1e6, results[0][tidx], \
         label=str(1e3*tidx*DT_sim)+'ms', linewidth=3)
    legend(numpoints=1)
    xlim(0, 1000)
    ylim(-80,40)
    xlabel('Z-axis (um)')
    ylabel('Membrane potential (mV)')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# END
