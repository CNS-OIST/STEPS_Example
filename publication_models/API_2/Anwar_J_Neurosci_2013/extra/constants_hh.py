import steps.interface

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Okinawa Institute of Science and Technology, Japan.
#
# This script runs on STEPS 2.x http://steps.sourceforge.net
#
# H Anwar, I Hepburn, H Nedelescu, W Chen and E De Schutter
# Stochastic calcium mechanisms cause dendritic calcium spike variability
# J Neuroscience 2013
#
# constants_hh.py : provides a set of parameters and other constants for the
# Hodgkin-Huxley model in the above study. 
# It is intended that this file is not altered. 
#
# Script authors: Haroon Anwar and Iain Hepburn
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import math

# # # # # # # # # # # # # # # # SIMULATION CONTROLS # # # # # # # # # # # # #

EF_DT = 1.0e-5          # The EField dt
NTIMEPOINTS =  5000 

TIMECONVERTER =  1.0e-5

NITER = 1

############################ PARAMETERS ################################

init_pot = -65e-3

TEMPERATURE = 20.0

Q10 = 3

Qt = math.pow(Q10, ((TEMPERATURE-6.3)/10))

########## BULK RESISTIVITY ##########

Ra = 1.0

########## MEMBRANE CAPACITANCE ##########

memb_capac = 1.0e-2


# # # # # # # # # # # # # # # # # # CHANNELS  # # # # # # # # # # # # # # # #

# Voltage range for gating kinetics in Volts
Vrange = [-100.0e-3, 50e-3, 1e-4]

# Hodgkin-Huxley gating kinetics

def a_n(V):
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    return  ((0.01*(10-(V+65.))/(math.exp((10-(V+65.))/10.)-1)))

def b_n(V):
    # WARNING: Using a variable name that is reserved (['V']).
    return ((0.125*math.exp(-(V+65.)/80.)))

def a_m(V):
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    return ((0.1*(25-(V+65.))/(math.exp((25-(V+65.))/10.)-1)))

def b_m(V):
    # WARNING: Using a variable name that is reserved (['V']).
    return ((4.*math.exp(-(V+65.)/18.)))

def a_h(V):
    # WARNING: Using a variable name that is reserved (['V']).
    return ((0.07*math.exp(-(V+65.)/20.)))

def b_h(V):
    # WARNING: Using a variable name that is reserved (['V']).
    return ((1./(math.exp((30-(V+65.))/10.)+1)))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Potassium conductance = 0.036 S/cm2
# Sodium conductance = 0.120 S/cm2

# Potassium single-channel conductance
K_G = 20.0e-12 # Siemens

# Potassium channel density
K_ro = 18.0e12 # per square meter

# Potassium reversal potential
K_rev = -77e-3 # volts

# Sodium single-channel conductance
Na_G = 20.0e-12 # Siemens

# Sodium channel density
Na_ro = 60.0e12 # per square meter

# Sodium reversal potential
Na_rev = 50e-3 # volts

# Leak single-channel conductance
L_G = 1.0e-12 # Siemens

# Leak density
L_ro = 10.0e12 # per square meter

# Leak reveral potential
leak_rev = -50.0e-3 # volts


# A table of potassium channel initial population factors:
# n0, n1, n2, n3, n4
K_facs = [ 0.21768, 0.40513, 0.28093, 0.08647, 0.00979 ]

# A table of sodium channel initial population factors
# m0h0, m1h0, m2h0, m3h0, m0h1, m1h1, m2h1, m3h1:
Na_facs = [ 0.34412, 0.05733, 0.00327, 6.0e-05, \
                0.50558, 0.08504, 0.00449, 0.00010 ]
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


