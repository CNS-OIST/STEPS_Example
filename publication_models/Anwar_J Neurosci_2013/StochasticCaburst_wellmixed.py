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
# *StochasticCaburst_wellmixed.py : The stochastic calcium burst model, with 
# well-mixed cytosolic compartments.
#
# Script authors: Haroon Anwar and Iain Hepburn
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# USAGE
#
# $ python StochasticCaburst_wellmixed.py *mesh* *root* *iter_n* 
#  
#  *mesh* is the tetrahedral mesh (10um to 160um cylinder)
#  *root* is the path to the location for data storage 
#  *iter_n* (is intened to be an integer) is an identifier number for each 
#     simulation iteration. 
#
# E.g: 
# $ python StochasticCaburst_wellmixed.py Cylinder2_dia2um_L10um_outer0_3um_0.3shell_0.3size_19156tets_adaptive.inp ~/stochcasims/ 1
#
#
# OUTPUT 
#
# In (root)/data/StochasticCaburst_wellmixed/(mesh)/(iter_n+time) directory 
# 3 data files will be recorded. Each file contains one row for every 
# time-point at which data is recorded, organised into the following columns:
# 
# currents.dat
# Time (ms), P-type current, T-type current, BK current, SK current
# (current units are Amps/m^2)
#
# voltage.dat
# Time (ms), voltage at mesh centre (mV)
#
# calcium.dat
# Time (ms), number of calcium ions in submembrane shell,
# calcium concentration in outer submembrane (micromolar). 
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import math
import time
from random import *
import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.utilities.meshio as meshio
import steps.solver as ssolver
import os

import meshes.gettets as gettets
from extra.constants import *
import extra.curr_funcs as cf

import sys

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

meshfile_ab, root, iter_n = sys.argv[1], sys.argv[2], sys.argv[3]

if meshfile_ab == 'Cylinder2_dia2um_L160um_outer0_0.3shell_0.3size_279152tets_adaptive.inp': cyl160=True
else: cyl160=False

########################### BIOCHEMICAL MODEL ###############################

# Two models required: Stochastic and deterministic
mdl_stoch = smodel.Model()
mdl_WM = smodel.Model()

# Calcium
Ca = smodel.Spec('Ca', mdl_WM)
Ca.setValence(2)


# Pump
Pump = smodel.Spec('Pump', mdl_WM)
# CaPump
CaPump = smodel.Spec('CaPump', mdl_WM)

# iCBsf
iCBsf = smodel.Spec('iCBsf', mdl_WM)
# iCBsCa
iCBsCa = smodel.Spec('iCBsCa', mdl_WM)
# iCBCaf
iCBCaf = smodel.Spec('iCBCaf', mdl_WM)
# iCBCaCa
iCBCaCa = smodel.Spec('iCBCaCa', mdl_WM)

# CBsf
CBsf = smodel.Spec('CBsf', mdl_WM)
# CBsCa
CBsCa = smodel.Spec('CBsCa', mdl_WM)
# CBCaf
CBCaf = smodel.Spec('CBCaf', mdl_WM)
# CBCaCa
CBCaCa = smodel.Spec('CBCaCa', mdl_WM)

# PV
PV = smodel.Spec('PV', mdl_WM)
# PVMg
PVMg = smodel.Spec('PVMg', mdl_WM)
# PVCa
PVCa = smodel.Spec('PVCa', mdl_WM)
# Mg
Mg = smodel.Spec('Mg', mdl_WM)

ssys_diff = smodel.Surfsys('ssys_diff', mdl_WM)

ssys_chans = smodel.Surfsys('ssys_chans', mdl_WM)

vsys_buffs = smodel.Volsys('vsys_buffs', mdl_WM)

#for EField
ssys_stoch = smodel.Surfsys('ssys_stoch', mdl_stoch)


# Diffusions 
diff_Ca_inward = smodel.SReac('diff_Ca_inward', ssys_diff, olhs=[Ca], irhs=[Ca], kcst=0)
diff_Ca_outward = smodel.SReac('diff_Ca_outward', ssys_diff, ilhs=[Ca], orhs=[Ca], kcst=0)

diff_CBsf_inward = smodel.SReac('diff_CBsf_inward', ssys_diff, olhs=[CBsf], irhs=[CBsf], kcst=0)
diff_CBsf_outward = smodel.SReac('diff_CBsf_outward', ssys_diff, ilhs=[CBsf], orhs=[CBsf], kcst=0)

diff_CBsCa_inward = smodel.SReac('diff_CBsCa_inward', ssys_diff, olhs=[CBsCa], irhs=[CBsCa], kcst=0)
diff_CBsCa_outward = smodel.SReac('diff_CBsCa_outward', ssys_diff, ilhs=[CBsCa], orhs=[CBsCa], kcst=0)

diff_CBCaf_inward = smodel.SReac('diff_CBCaf_inward', ssys_diff, olhs=[CBCaf], irhs=[CBCaf], kcst=0)
diff_CBCaF_outward = smodel.SReac('diff_CBCaf_outward', ssys_diff, ilhs=[CBCaf], orhs=[CBCaf], kcst=0)

diff_CBCaCa_inward = smodel.SReac('diff_CBCaCa_inward', ssys_diff, olhs=[CBCaCa], irhs=[CBCaCa], kcst=0)
diff_CBCaCa_outward = smodel.SReac('diff_CBCaCa_outward', ssys_diff, ilhs=[CBCaCa], orhs=[CBCaCa], kcst=0)

diff_PV_inward = smodel.SReac('diff_PV_inward', ssys_diff, olhs=[PV], irhs=[PV], kcst=0)
diff_PV_outward = smodel.SReac('diff_PV_outward', ssys_diff, ilhs=[PV], orhs=[PV], kcst=0)

diff_PVCa_inward = smodel.SReac('diff_PVCa_inward', ssys_diff, olhs=[PVCa], irhs=[PVCa], kcst=0)
diff_PVCa_outward = smodel.SReac('diff_PVCa_outward', ssys_diff, ilhs=[PVCa], orhs=[PVCa], kcst=0)

diff_PVMg_inward = smodel.SReac('diff_PVMg_inward', ssys_diff, olhs=[PVMg], irhs=[PVMg], kcst=0)
diff_PVMg_outward = smodel.SReac('diff_PVMg_outward', ssys_diff, ilhs=[PVMg], orhs=[PVMg], kcst=0)


#Pump
PumpD_f = smodel.SReac('PumpD_f', ssys_chans, ilhs=[Ca], slhs=[Pump], srhs=[CaPump])
PumpD_f.setKcst(P_f_kcst)

PumpD_b = smodel.SReac('PumpD_b', ssys_chans, slhs=[CaPump], irhs=[Ca], srhs=[Pump])
PumpD_b.setKcst(P_b_kcst)

PumpD_k = smodel.SReac('PumpD_k', ssys_chans, slhs=[CaPump], srhs=[Pump])
PumpD_k.setKcst(P_k_kcst)

#iCBsf-fast
iCBsf1_f = smodel.Reac('iCBsf1_f', vsys_buffs, lhs=[Ca,iCBsf], rhs=[iCBsCa], kcst = iCBsf1_f_kcst)
iCBsf1_b = smodel.Reac('iCBsf1_b', vsys_buffs, lhs=[iCBsCa], rhs=[Ca, iCBsf], kcst = iCBsf1_b_kcst)

#iCBsCa
iCBsCa_f = smodel.Reac('iCBsCa_f', vsys_buffs, lhs=[Ca,iCBsCa], rhs=[iCBCaCa], kcst = iCBsCa_f_kcst)
iCBsCa_b = smodel.Reac('iCBsCa_b', vsys_buffs, lhs=[iCBCaCa], rhs=[Ca,iCBsCa], kcst = iCBsCa_b_kcst)

#iCBsf_slow
iCBsf2_f = smodel.Reac('iCBsf2_f', vsys_buffs, lhs=[Ca,iCBsf], rhs=[iCBCaf], kcst = iCBsf2_f_kcst)
iCBsf2_b = smodel.Reac('iCBsf2_b', vsys_buffs, lhs=[iCBCaf], rhs=[Ca,iCBsf], kcst = iCBsf2_b_kcst)

#iCBCaf
iCBCaf_f = smodel.Reac('iCBCaf_f', vsys_buffs, lhs=[Ca,iCBCaf], rhs=[iCBCaCa], kcst = iCBCaf_f_kcst)
iCBCaf_b = smodel.Reac('iCBCaf_b', vsys_buffs, lhs=[iCBCaCa], rhs=[Ca,iCBCaf], kcst = iCBCaf_b_kcst)

#CBsf-fast
CBsf1_f = smodel.Reac('CBsf1_f', vsys_buffs, lhs=[Ca,CBsf], rhs=[CBsCa], kcst = CBsf1_f_kcst)
CBsf1_b = smodel.Reac('CBsf1_b', vsys_buffs, lhs=[CBsCa], rhs=[Ca,CBsf], kcst = CBsf1_b_kcst)

#CBsCa
CBsCa_f = smodel.Reac('CBsCa_f', vsys_buffs, lhs=[Ca,CBsCa], rhs=[CBCaCa], kcst = CBsCa_f_kcst)
CBsCa_b = smodel.Reac('CBsCa_b', vsys_buffs, lhs=[CBCaCa], rhs=[Ca,CBsCa], kcst = CBsCa_b_kcst)

#CBsf_slow
CBsf2_f = smodel.Reac('CBsf2_f', vsys_buffs, lhs=[Ca,CBsf], rhs=[CBCaf], kcst = CBsf2_f_kcst)
CBsf2_b = smodel.Reac('CBsf2_b', vsys_buffs, lhs=[CBCaf], rhs=[Ca,CBsf], kcst = CBsf2_b_kcst)

#CBCaf
CBCaf_f = smodel.Reac('CBCaf_f', vsys_buffs, lhs=[Ca,CBCaf], rhs=[CBCaCa], kcst = CBCaf_f_kcst)
CBCaf_b = smodel.Reac('CBCaf_b', vsys_buffs, lhs=[CBCaCa], rhs=[Ca,CBCaf], kcst = CBCaf_b_kcst)

#PVca
PVca_f = smodel.Reac('PVca_f', vsys_buffs, lhs=[Ca,PV], rhs=[PVCa], kcst = PVca_f_kcst)
PVca_b = smodel.Reac('PVca_b', vsys_buffs, lhs=[PVCa], rhs=[Ca,PV], kcst = PVca_b_kcst)

#PVmg
PVmg_f = smodel.Reac('PVmg_f', vsys_buffs, lhs=[Mg,PV], rhs=[PVMg], kcst = PVmg_f_kcst)
PVmg_b = smodel.Reac('PVmg_b', vsys_buffs, lhs=[PVMg], rhs=[Mg,PV], kcst = PVmg_b_kcst)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # CHANNELS  # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###### CaP channel ############## 

CaP_m0 = smodel.Spec('CaP_m0', mdl_WM)
CaP_m1 = smodel.Spec('CaP_m1', mdl_WM)
CaP_m2 = smodel.Spec('CaP_m2', mdl_WM)
CaP_m3 = smodel.Spec('CaP_m3', mdl_WM)


CaPm0m1 = smodel.SReac('CaPm0m1', ssys_chans, slhs = [CaP_m0], srhs = [CaP_m1], kcst= 0.0)
CaPm1m2 = smodel.SReac('CaPm1m2', ssys_chans, slhs = [CaP_m1], srhs = [CaP_m2], kcst= 0.0)
CaPm2m3 = smodel.SReac('CaPm2m3', ssys_chans, slhs = [CaP_m2], srhs = [CaP_m3], kcst= 0.0)

CaPm3m2 = smodel.SReac('CaPm3m2', ssys_chans, slhs = [CaP_m3], srhs = [CaP_m2], kcst= 0.0)
CaPm2m1 = smodel.SReac('CaPm2m1', ssys_chans, slhs = [CaP_m2], srhs = [CaP_m1], kcst= 0.0)
CaPm1m0 = smodel.SReac('CaPm1m0', ssys_chans, slhs = [CaP_m1], srhs = [CaP_m0], kcst= 0.0)


######## CaT channel ##########  

CaT_m0h0 = smodel.Spec('CaT_m0h0', mdl_WM)
CaT_m0h1 = smodel.Spec('CaT_m0h1', mdl_WM)
CaT_m1h0 = smodel.Spec('CaT_m1h0', mdl_WM)
CaT_m1h1 = smodel.Spec('CaT_m1h1', mdl_WM)
CaT_m2h0 = smodel.Spec('CaT_m2h0', mdl_WM)
CaT_m2h1 = smodel.Spec('CaT_m2h1', mdl_WM)


CaTm0h0_m1h0 = smodel.SReac('CaTm0h0_m1h0', ssys_chans, slhs = [CaT_m0h0], srhs = [CaT_m1h0], kcst=0.0)
CaTm1h0_m2h0 = smodel.SReac('CaTm1h0_m2h0', ssys_chans, slhs = [CaT_m1h0], srhs = [CaT_m2h0], kcst=0.0)

CaTm2h0_m1h0 = smodel.SReac('CaTm2h0_m1h0', ssys_chans, slhs = [CaT_m2h0], srhs = [CaT_m1h0], kcst=0.0)
CaTm1h0_m0h0 = smodel.SReac('CaTm1h0_m0h0', ssys_chans, slhs = [CaT_m1h0], srhs = [CaT_m0h0], kcst=0.0)

CaTm0h1_m1h1 = smodel.SReac('CaTm0h1_m1h1', ssys_chans, slhs = [CaT_m0h1], srhs = [CaT_m1h1], kcst=0.0)
CaTm1h1_m2h1 = smodel.SReac('CaTm1h1_m2h1', ssys_chans, slhs = [CaT_m1h1], srhs = [CaT_m2h1], kcst=0.0)

CaTm2h1_m1h1 = smodel.SReac('CaTm2h1_m1h1', ssys_chans, slhs = [CaT_m2h1], srhs = [CaT_m1h1], kcst=0.0)
CaTm1h1_m0h1 = smodel.SReac('CaTm1h1_m0h1', ssys_chans, slhs = [CaT_m1h1], srhs = [CaT_m0h1], kcst=0.0)


CaTm0h0_m0h1 = smodel.SReac('CaTm0h0_m0h1', ssys_chans, slhs = [CaT_m0h0], srhs = [CaT_m0h1], kcst=0.0)
CaTm1h0_m1h1 = smodel.SReac('CaTm1h0_m1h1', ssys_chans, slhs = [CaT_m1h0], srhs = [CaT_m1h1], kcst=0.0)
CaTm2h0_m2h1 = smodel.SReac('CaTm2h0_m2h1', ssys_chans, slhs = [CaT_m2h0], srhs = [CaT_m2h1], kcst=0.0)

CaTm2h1_m2h0 = smodel.SReac('CaTm2h1_m2h0', ssys_chans, slhs = [CaT_m2h1], srhs = [CaT_m2h0], kcst=0.0)
CaTm1h1_m1h0 = smodel.SReac('CaTm1h1_m1h0', ssys_chans, slhs = [CaT_m1h1], srhs = [CaT_m1h0], kcst=0.0)
CaTm0h1_m0h0 = smodel.SReac('CaTm0h1_m0h0', ssys_chans, slhs = [CaT_m0h1], srhs = [CaT_m0h0], kcst=0.0)

##### BK channel ####################

BK_C0 = smodel.Spec('BK_C0', mdl_WM)
BK_C1 = smodel.Spec('BK_C1', mdl_WM)
BK_C2 = smodel.Spec('BK_C2', mdl_WM)
BK_C3 = smodel.Spec('BK_C3', mdl_WM)
BK_C4 = smodel.Spec('BK_C4', mdl_WM)
BK_O0 = smodel.Spec('BK_O0', mdl_WM)
BK_O1 = smodel.Spec('BK_O1', mdl_WM)
BK_O2 = smodel.Spec('BK_O2', mdl_WM)
BK_O3 = smodel.Spec('BK_O3', mdl_WM)
BK_O4 = smodel.Spec('BK_O4', mdl_WM)


BKCAC0 = smodel.SReac('BKCAC0', ssys_chans, slhs = [BK_C0], ilhs = [Ca], srhs = [BK_C1], kcst = c_01)
BKCAC1 = smodel.SReac('BKCAC1', ssys_chans, slhs = [BK_C1], ilhs = [Ca], srhs = [BK_C2], kcst = c_12)
BKCAC2 = smodel.SReac('BKCAC2', ssys_chans, slhs = [BK_C2], ilhs = [Ca], srhs = [BK_C3], kcst = c_23)
BKCAC3 = smodel.SReac('BKCAC3', ssys_chans, slhs = [BK_C3], ilhs = [Ca], srhs = [BK_C4], kcst = c_34)

BKC0 = smodel.SReac('BKC0', ssys_chans, slhs = [BK_C1], srhs = [BK_C0], irhs=[Ca], kcst = c_10)
BKC1 = smodel.SReac('BKC1', ssys_chans, slhs = [BK_C2], srhs = [BK_C1], irhs=[Ca], kcst = c_21)
BKC2 = smodel.SReac('BKC2', ssys_chans, slhs = [BK_C3], srhs = [BK_C2], irhs=[Ca], kcst = c_32)
BKC3 = smodel.SReac('BKC3', ssys_chans, slhs = [BK_C4], srhs = [BK_C3], irhs=[Ca], kcst = c_43)

BKCAO0 = smodel.SReac('BKCAO0', ssys_chans, slhs = [BK_O0], ilhs = [Ca], srhs = [BK_O1], kcst = o_01)
BKCAO1 = smodel.SReac('BKCAO1', ssys_chans, slhs = [BK_O1], ilhs = [Ca], srhs = [BK_O2], kcst = o_12)
BKCAO2 = smodel.SReac('BKCAO2', ssys_chans, slhs = [BK_O2], ilhs = [Ca], srhs = [BK_O3], kcst = o_23)
BKCAO3 = smodel.SReac('BKCAO3', ssys_chans, slhs = [BK_O3], ilhs = [Ca], srhs = [BK_O4], kcst = o_34)

BKO0 = smodel.SReac('BKO0', ssys_chans, slhs = [BK_O1], srhs = [BK_O0], irhs=[Ca], kcst = o_10)
BKO1 = smodel.SReac('BKO1', ssys_chans, slhs = [BK_O2], srhs = [BK_O1], irhs=[Ca], kcst = o_21)
BKO2 = smodel.SReac('BKO2', ssys_chans, slhs = [BK_O3], srhs = [BK_O2], irhs=[Ca], kcst = o_32)
BKO3 = smodel.SReac('BKO3', ssys_chans, slhs = [BK_O4], srhs = [BK_O3], irhs=[Ca], kcst = o_43)

BKC0O0 = smodel.SReac('BKC0O0', ssys_chans, slhs = [BK_C0], srhs = [BK_O0], kcst=0.0)
BKC1O1 = smodel.SReac('BKC1O1', ssys_chans, slhs = [BK_C1], srhs = [BK_O1], kcst=0.0)
BKC2O2 = smodel.SReac('BKC2O2', ssys_chans, slhs = [BK_C2], srhs = [BK_O2], kcst=0.0)
BKC3O3 = smodel.SReac('BKC3O3', ssys_chans, slhs = [BK_C3], srhs = [BK_O3], kcst=0.0)
BKC4O4 = smodel.SReac('BKC4O4', ssys_chans, slhs = [BK_C4], srhs = [BK_O4], kcst=0.0)

BKO0C0 = smodel.SReac('BKO0C0', ssys_chans, slhs = [BK_O0], srhs = [BK_C0], kcst=0.0)
BKO1C1 = smodel.SReac('BKO1C1', ssys_chans, slhs = [BK_O1], srhs = [BK_C1], kcst=0.0)
BKO2C2 = smodel.SReac('BKO2C2', ssys_chans, slhs = [BK_O2], srhs = [BK_C2], kcst=0.0)
BKO3C3 = smodel.SReac('BKO3C3', ssys_chans, slhs = [BK_O3], srhs = [BK_C3], kcst=0.0)
BKO4C4 = smodel.SReac('BKO4C4', ssys_chans, slhs = [BK_O4], srhs = [BK_C4], kcst=0.0)


###### SK channel ################## DETERMINISTIC

SK_C1 = smodel.Spec('SK_C1', mdl_WM)
SK_C2 = smodel.Spec('SK_C2', mdl_WM)
SK_C3 = smodel.Spec('SK_C3', mdl_WM)
SK_C4 = smodel.Spec('SK_C4', mdl_WM)
SK_O1 = smodel.Spec('SK_O1', mdl_WM)
SK_O2 = smodel.Spec('SK_O2', mdl_WM)


SKCAC1 = smodel.SReac('SKCAC1', ssys_chans, slhs = [SK_C1], ilhs = [Ca], srhs = [SK_C2], kcst = dirc2_t)
SKCAC2 = smodel.SReac('SKCAC2', ssys_chans, slhs = [SK_C2], ilhs = [Ca], srhs = [SK_C3], kcst = dirc3_t)
SKCAC3 = smodel.SReac('SKCAC3', ssys_chans, slhs = [SK_C3], ilhs = [Ca], srhs = [SK_C4], kcst = dirc4_t)

SKC1 = smodel.SReac('SKC1', ssys_chans, slhs = [SK_C2], srhs = [SK_C1], irhs=[Ca], kcst = invc1_t)
SKC2 = smodel.SReac('SKC2', ssys_chans, slhs = [SK_C3], srhs = [SK_C2], irhs=[Ca], kcst = invc2_t)
SKC3 = smodel.SReac('SKC3', ssys_chans, slhs = [SK_C4], srhs = [SK_C3], irhs=[Ca], kcst = invc3_t)

SKC3O1 = smodel.SReac('SKC3O1', ssys_chans, slhs = [SK_C3], srhs = [SK_O1], kcst = diro1_t)
SKC4O2 = smodel.SReac('SKC4O2', ssys_chans, slhs = [SK_C4], srhs = [SK_O2], kcst = diro2_t)

SKO1C3 = smodel.SReac('SKO1C3', ssys_chans, slhs = [SK_O1], srhs = [SK_C3], kcst = invo1_t)
SKO2C4 = smodel.SReac('SKO2C4', ssys_chans, slhs = [SK_O2], srhs = [SK_C4], kcst = invo2_t)


###### Leak current channel #####

L = smodel.Chan('L', mdl_stoch)
Leak = smodel.ChanState('Leak', mdl_stoch, L)

OC_L = smodel.OhmicCurr('OC_L', ssys_stoch, chanstate = Leak, erev = L_rev, g = L_G)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

########### MESH & COMPARTMENTALIZATION #################

##########Import Mesh for EField calculation

# For stochastic sim:
mesh_stoch = meshio.loadMesh('./meshes/'+meshfile_ab)[0]

outer_tets = range(mesh_stoch.ntets)
inner_tets = gettets.getcyl(mesh_stoch, 1e-6, -200e-6, 200e-6)[0]

for i in inner_tets: outer_tets.remove(i)

print outer_tets.__len__(), " tets in outer compartment"
print inner_tets.__len__(), " tets in inner compartment"

# Record voltage from the central tetrahedron
cent_tet = mesh_stoch.findTetByPoint([0.0,0.0,0.0])

#Define geometrical constants for a compartment with concentric shells

Length = 80e-6

shell0_rmax = 1.0e-6
shell1_rmax = 0.9e-6
shell2_rmax = 0.7e-6
shell3_rmax = 0.5e-6
shell4_rmax = 0.3e-6
shell5_rmax = 0.1e-6

#Geometry container object:

wmgeom = sgeom.Geom()

#Set up compartments

shell0 = sgeom.Comp('shell0', wmgeom, vol=math.pi*((shell0_rmax**2) - (shell1_rmax**2))*Length)
shell1 = sgeom.Comp('shell1', wmgeom, vol=math.pi*((shell1_rmax**2) - (shell2_rmax**2))*Length)
shell2 = sgeom.Comp('shell2', wmgeom, vol=math.pi*((shell2_rmax**2) - (shell3_rmax**2))*Length)
shell3 = sgeom.Comp('shell3', wmgeom, vol=math.pi*((shell3_rmax**2) - (shell4_rmax**2))*Length)
shell4 = sgeom.Comp('shell4', wmgeom, vol=math.pi*((shell4_rmax**2) - (shell5_rmax**2))*Length)
shell5 = sgeom.Comp('shell5', wmgeom, vol=math.pi*((shell5_rmax**2))*Length)

#Set up patches: lower number compartment is 'outer', higher number is 'inner'
ring0 = sgeom.Patch('ring0', wmgeom, icomp=shell0, area=2*math.pi*shell0_rmax*Length)

ring1 = sgeom.Patch('ring1', wmgeom, ocomp=shell0, icomp=shell1, area=2*math.pi*shell1_rmax*Length)
ring2 = sgeom.Patch('ring2', wmgeom, ocomp=shell1, icomp=shell2, area=2*math.pi*shell2_rmax*Length)
ring3 = sgeom.Patch('ring3', wmgeom, ocomp=shell2, icomp=shell3, area=2*math.pi*shell3_rmax*Length)
ring4 = sgeom.Patch('ring4', wmgeom, ocomp=shell3, icomp=shell4, area=2*math.pi*shell4_rmax*Length)
ring5 = sgeom.Patch('ring5', wmgeom, ocomp=shell4, icomp=shell5, area=2*math.pi*shell5_rmax*Length)

#For all patches, add the surface system

ring0.addSurfsys('ssys_chans')

ring1.addSurfsys('ssys_diff')
ring2.addSurfsys('ssys_diff')
ring3.addSurfsys('ssys_diff')
ring4.addSurfsys('ssys_diff')
ring5.addSurfsys('ssys_diff')

shell0.addVolsys('vsys_buffs')
shell1.addVolsys('vsys_buffs')
shell2.addVolsys('vsys_buffs')
shell3.addVolsys('vsys_buffs')
shell4.addVolsys('vsys_buffs')
shell5.addVolsys('vsys_buffs')

########## Create an intracellular compartment i.e. cytosolic compartment

cyto_stoch = sgeom.TmComp('cyto_stoch', mesh_stoch, inner_tets)


if cyl160:
    # Ensure that we use points a small distance inside the boundary:
    LENGTH = mesh_stoch.getBoundMax()[2] - mesh_stoch.getBoundMin()[2]
    boundminz = mesh_stoch.getBoundMin()[2] + LENGTH/mesh_stoch.ntets
    boundmaxz = mesh_stoch.getBoundMax()[2] - LENGTH/mesh_stoch.ntets

    memb_tris = list(mesh_stoch.getSurfTris())
    minztris = []
    maxztris = []
    for tri in memb_tris:
        zminboundtri = True
        zmaxboundtri = True
        tritemp = mesh_stoch.getTri(tri)
        trizs = [0.0, 0.0, 0.0]
        trizs[0] = mesh_stoch.getVertex(tritemp[0])[2]
        trizs[1] = mesh_stoch.getVertex(tritemp[1])[2]
        trizs[2] = mesh_stoch.getVertex(tritemp[2])[2]
        for j in range(3):
            if (trizs[j]>boundminz): zminboundtri = False
        if (zminboundtri):
            minztris.append(tri)
            continue
        for j in range(3):
            if (trizs[j]< boundmaxz): zmaxboundtri = False
        if (zmaxboundtri):
            maxztris.append(tri)

    for t in minztris: memb_tris.remove(t)
    for t in maxztris: memb_tris.remove(t)
    
else:
    print 'Finding connecting triangles...'
    out_tris = set()
    for i in outer_tets:
            tritemp = mesh_stoch.getTetTriNeighb(i)
            for j in range(4): out_tris.add(tritemp[j])

    in_tris = set()
    for i in inner_tets:
            tritemp = mesh_stoch.getTetTriNeighb(i)
            for j in range(4): in_tris.add(tritemp[j])

    memb_tris = out_tris.intersection(in_tris)
    memb_tris = list(memb_tris)

border_tets = []
border_tets_vols = 0.0

for i in inner_tets:
    tritemp = mesh_stoch.getTetTriNeighb(i)
    for t in tritemp:
        if t in memb_tris:
            border_tets.append(i)
            border_tets_vols+=mesh_stoch.getTetVol(i)
            break

print "Border tet vols:", border_tets_vols

########## Create a membrane as a surface mesh

# Stochastic sim:
memb_stoch = sgeom.TmPatch('memb_stoch', mesh_stoch, memb_tris, cyto_stoch)
memb_stoch.addSurfsys('ssys_stoch')

# For EField calculation
print "Creating membrane.."
membrane = sgeom.Memb('membrane', mesh_stoch, [memb_stoch])
print "Membrane created."

# # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

r = srng.create_mt19937(512)
r.initialize(int(time.time()%1000))

# Random-number generator 
r_dummy = srng.create_mt19937(512)
r_dummy.initialize(int(time.time()%1000))

print "Creating tetexact solver..."
sim_stoch = ssolver.Tetexact(mdl_stoch, mesh_stoch, r, True)

print "Creating WM solver"
sim_WM = ssolver.Wmdirect(mdl_WM, wmgeom, r_dummy)


print "Resetting simulation objects.."
sim_stoch.reset()
sim_WM.reset()


print "Injecting molecules.."

sim_stoch.setTemp(TEMPERATURE+273.15)

sim_WM.setCompConc('shell0', 'Ca', Ca_iconc)
sim_WM.setCompConc('shell1', 'Ca', Ca_iconc)
sim_WM.setCompConc('shell2', 'Ca', Ca_iconc)
sim_WM.setCompConc('shell3', 'Ca', Ca_iconc)
sim_WM.setCompConc('shell4', 'Ca', Ca_iconc)
sim_WM.setCompConc('shell5', 'Ca', Ca_iconc)

sim_WM.setCompConc('shell0', 'Mg', Mg_conc)
sim_WM.setCompConc('shell1', 'Mg', Mg_conc)
sim_WM.setCompConc('shell2', 'Mg', Mg_conc)
sim_WM.setCompConc('shell3', 'Mg', Mg_conc)
sim_WM.setCompConc('shell4', 'Mg', Mg_conc)
sim_WM.setCompConc('shell5', 'Mg', Mg_conc)

surfarea = sim_WM.getPatchArea('ring0')


pumpnbs = 6.022141e12*surfarea

sim_WM.setPatchCount('ring0', 'Pump', pumpnbs)
sim_WM.setPatchCount('ring0', 'CaPump', 0)

print "Injected ", sim_WM.getPatchCount('ring0', 'Pump'), "pumps"

sim_WM.setCompConc('shell0', 'iCBsf',    iCBsf_conc)
sim_WM.setCompConc('shell0', 'iCBsCa',   iCBsCa_conc)
sim_WM.setCompConc('shell0', 'iCBCaf',   iCBCaf_conc)
sim_WM.setCompConc('shell0', 'iCBCaCa',  iCBCaCa_conc)

sim_WM.setCompConc('shell0', 'CBsf',     CBsf_conc)
sim_WM.setCompConc('shell0', 'CBsCa',    CBsCa_conc)
sim_WM.setCompConc('shell0', 'CBCaf',    CBCaf_conc)
sim_WM.setCompConc('shell0', 'CBCaCa',   CBCaCa_conc)

sim_WM.setCompConc('shell0', 'PV',       PV_conc)
sim_WM.setCompConc('shell0', 'PVCa',     PVCa_conc)
sim_WM.setCompConc('shell0', 'PVMg',     PVMg_conc)

sim_WM.setCompConc('shell1', 'iCBsf',    iCBsf_conc)
sim_WM.setCompConc('shell1', 'iCBsCa',   iCBsCa_conc)
sim_WM.setCompConc('shell1', 'iCBCaf',   iCBCaf_conc)
sim_WM.setCompConc('shell1', 'iCBCaCa',  iCBCaCa_conc)

sim_WM.setCompConc('shell1', 'CBsf',     CBsf_conc)
sim_WM.setCompConc('shell1', 'CBsCa',    CBsCa_conc)
sim_WM.setCompConc('shell1', 'CBCaf',    CBCaf_conc)
sim_WM.setCompConc('shell1', 'CBCaCa',   CBCaCa_conc)

sim_WM.setCompConc('shell1', 'PV',       PV_conc)
sim_WM.setCompConc('shell1', 'PVCa',     PVCa_conc)
sim_WM.setCompConc('shell1', 'PVMg',     PVMg_conc)

sim_WM.setCompConc('shell2', 'iCBsf',    iCBsf_conc)
sim_WM.setCompConc('shell2', 'iCBsCa',   iCBsCa_conc)
sim_WM.setCompConc('shell2', 'iCBCaf',   iCBCaf_conc)
sim_WM.setCompConc('shell2', 'iCBCaCa',  iCBCaCa_conc)

sim_WM.setCompConc('shell2', 'CBsf',     CBsf_conc)
sim_WM.setCompConc('shell2', 'CBsCa',    CBsCa_conc)
sim_WM.setCompConc('shell2', 'CBCaf',    CBCaf_conc)
sim_WM.setCompConc('shell2', 'CBCaCa',   CBCaCa_conc)

sim_WM.setCompConc('shell2', 'PV',       PV_conc)
sim_WM.setCompConc('shell2', 'PVCa',     PVCa_conc)
sim_WM.setCompConc('shell2', 'PVMg',     PVMg_conc)

sim_WM.setCompConc('shell3', 'iCBsf',    iCBsf_conc)
sim_WM.setCompConc('shell3', 'iCBsCa',   iCBsCa_conc)
sim_WM.setCompConc('shell3', 'iCBCaf',   iCBCaf_conc)
sim_WM.setCompConc('shell3', 'iCBCaCa',  iCBCaCa_conc)

sim_WM.setCompConc('shell3', 'CBsf',     CBsf_conc)
sim_WM.setCompConc('shell3', 'CBsCa',    CBsCa_conc)
sim_WM.setCompConc('shell3', 'CBCaf',    CBCaf_conc)
sim_WM.setCompConc('shell3', 'CBCaCa',   CBCaCa_conc)

sim_WM.setCompConc('shell3', 'PV',       PV_conc)
sim_WM.setCompConc('shell3', 'PVCa',     PVCa_conc)
sim_WM.setCompConc('shell3', 'PVMg',     PVMg_conc)

sim_WM.setCompConc('shell4', 'iCBsf',    iCBsf_conc)
sim_WM.setCompConc('shell4', 'iCBsCa',   iCBsCa_conc)
sim_WM.setCompConc('shell4', 'iCBCaf',   iCBCaf_conc)
sim_WM.setCompConc('shell4', 'iCBCaCa',  iCBCaCa_conc)

sim_WM.setCompConc('shell4', 'CBsf',     CBsf_conc)
sim_WM.setCompConc('shell4', 'CBsCa',    CBsCa_conc)
sim_WM.setCompConc('shell4', 'CBCaf',    CBCaf_conc)
sim_WM.setCompConc('shell4', 'CBCaCa',   CBCaCa_conc)

sim_WM.setCompConc('shell4', 'PV',       PV_conc)
sim_WM.setCompConc('shell4', 'PVCa',     PVCa_conc)
sim_WM.setCompConc('shell4', 'PVMg',     PVMg_conc)

sim_WM.setCompConc('shell5', 'iCBsf',    iCBsf_conc)
sim_WM.setCompConc('shell5', 'iCBsCa',   iCBsCa_conc)
sim_WM.setCompConc('shell5', 'iCBCaf',   iCBCaf_conc)
sim_WM.setCompConc('shell5', 'iCBCaCa',  iCBCaCa_conc)

sim_WM.setCompConc('shell5', 'CBsf',     CBsf_conc)
sim_WM.setCompConc('shell5', 'CBsCa',    CBsCa_conc)
sim_WM.setCompConc('shell5', 'CBCaf',    CBCaf_conc)
sim_WM.setCompConc('shell5', 'CBCaCa',   CBCaCa_conc)

sim_WM.setCompConc('shell5', 'PV',       PV_conc)
sim_WM.setCompConc('shell5', 'PVCa',     PVCa_conc)
sim_WM.setCompConc('shell5', 'PVMg',     PVMg_conc)


# CaP
sim_WM.setPatchCount('ring0', 'CaP_m0' , round(CaP_ro*surfarea*CaP_m0_p))
sim_WM.setPatchCount('ring0', 'CaP_m1' , round(CaP_ro*surfarea*CaP_m1_p))
sim_WM.setPatchCount('ring0', 'CaP_m2' , round(CaP_ro*surfarea*CaP_m2_p))
sim_WM.setPatchCount('ring0', 'CaP_m3' , round(CaP_ro*surfarea*CaP_m3_p))

print "Injected  ", CaP_ro*surfarea, "CaP channels"

# CaT

# From cstate: CaT_m2h0 conducting
sim_WM.setPatchCount('ring0', 'CaT_m0h0' , round(CaT_ro*surfarea*CaT_m0h0_p))
sim_WM.setPatchCount('ring0', 'CaT_m1h0' , round(CaT_ro*surfarea*CaT_m1h0_p))
sim_WM.setPatchCount('ring0', 'CaT_m2h0' , round(CaT_ro*surfarea*CaT_m2h0_p))
sim_WM.setPatchCount('ring0', 'CaT_m0h1' , round(CaT_ro*surfarea*CaT_m0h1_p))
sim_WM.setPatchCount('ring0', 'CaT_m1h1' , round(CaT_ro*surfarea*CaT_m1h1_p))
sim_WM.setPatchCount('ring0', 'CaT_m2h1' , round(CaT_ro*surfarea*CaT_m2h1_p))

print "Injected  ", CaT_ro*surfarea, "CaT channels"

# BK
sim_WM.setPatchCount('ring0', 'BK_C0' , round(BK_ro*surfarea*BK_C0_p))
sim_WM.setPatchCount('ring0', 'BK_C1' , round(BK_ro*surfarea*BK_C1_p))
sim_WM.setPatchCount('ring0', 'BK_C2' , round(BK_ro*surfarea*BK_C2_p))
sim_WM.setPatchCount('ring0', 'BK_C3' , round(BK_ro*surfarea*BK_C3_p))
sim_WM.setPatchCount('ring0', 'BK_C4' , round(BK_ro*surfarea*BK_C4_p))

sim_WM.setPatchCount('ring0', 'BK_O0' , round(BK_ro*surfarea*BK_O0_p))
sim_WM.setPatchCount('ring0', 'BK_O1' , round(BK_ro*surfarea*BK_O1_p))
sim_WM.setPatchCount('ring0', 'BK_O2' , round(BK_ro*surfarea*BK_O2_p))
sim_WM.setPatchCount('ring0', 'BK_O3' , round(BK_ro*surfarea*BK_O3_p))
sim_WM.setPatchCount('ring0', 'BK_O4' , round(BK_ro*surfarea*BK_O4_p))


print "Injected  ", BK_ro*surfarea, "BK channels"

# SK
sim_WM.setPatchCount('ring0', 'SK_C1' , round(SK_ro*surfarea*SK_C1_p))
sim_WM.setPatchCount('ring0', 'SK_C2' , round(SK_ro*surfarea*SK_C2_p))
sim_WM.setPatchCount('ring0', 'SK_C3' , round(SK_ro*surfarea*SK_C3_p))
sim_WM.setPatchCount('ring0', 'SK_C4' , round(SK_ro*surfarea*SK_C4_p))

sim_WM.setPatchCount('ring0', 'SK_O1' , round(SK_ro*surfarea*SK_O1_p))
sim_WM.setPatchCount('ring0', 'SK_O2' , round(SK_ro*surfarea*SK_O2_p))

#set the rate constants for diffusion (Diffusion is modeled as surface reaction here)

#for Ca diffusion

sim_WM.setPatchSReacK('ring1', 'diff_Ca_inward' ,(DCST*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_Ca_inward' ,(DCST*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_Ca_inward' ,(DCST*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_Ca_inward' ,(DCST*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_Ca_inward' ,(DCST*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))

sim_WM.setPatchSReacK('ring1', 'diff_Ca_outward' ,(DCST*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_Ca_outward' ,(DCST*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_Ca_outward' ,(DCST*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_Ca_outward' ,(DCST*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_Ca_outward' ,(DCST*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6)))

#for CBsf diffusin
sim_WM.setPatchSReacK('ring1', 'diff_CBsf_inward' ,(DCB*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_CBsf_inward' ,(DCB*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_CBsf_inward' ,(DCB*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_CBsf_inward' ,(DCB*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_CBsf_inward' ,(DCB*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))

sim_WM.setPatchSReacK('ring1', 'diff_CBsf_outward' ,(DCB*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_CBsf_outward' ,(DCB*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_CBsf_outward' ,(DCB*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_CBsf_outward' ,(DCB*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_CBsf_outward' ,(DCB*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6)))

#for CBsCa diffusion

sim_WM.setPatchSReacK('ring1', 'diff_CBsCa_inward' ,(DCB*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_CBsCa_inward' ,(DCB*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_CBsCa_inward' ,(DCB*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_CBsCa_inward' ,(DCB*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_CBsCa_inward' ,(DCB*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))

sim_WM.setPatchSReacK('ring1', 'diff_CBsCa_outward' ,(DCB*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_CBsCa_outward' ,(DCB*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_CBsCa_outward' ,(DCB*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_CBsCa_outward' ,(DCB*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_CBsCa_outward' ,(DCB*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6)))

#for CBCaf diffusion

sim_WM.setPatchSReacK('ring1', 'diff_CBCaf_inward' ,(DCB*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_CBCaf_inward' ,(DCB*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_CBCaf_inward' ,(DCB*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_CBCaf_inward' ,(DCB*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_CBCaf_inward' ,(DCB*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))

sim_WM.setPatchSReacK('ring1', 'diff_CBCaf_outward' ,(DCB*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_CBCaf_outward' ,(DCB*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_CBCaf_outward' ,(DCB*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_CBCaf_outward' ,(DCB*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_CBCaf_outward' ,(DCB*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6)))

#for CBCaCa diffusion

sim_WM.setPatchSReacK('ring1', 'diff_CBCaCa_inward' ,(DCB*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_CBCaCa_inward' ,(DCB*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_CBCaCa_inward' ,(DCB*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_CBCaCa_inward' ,(DCB*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_CBCaCa_inward' ,(DCB*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))

sim_WM.setPatchSReacK('ring1', 'diff_CBCaCa_outward' ,(DCB*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_CBCaCa_outward' ,(DCB*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_CBCaCa_outward' ,(DCB*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_CBCaCa_outward' ,(DCB*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_CBCaCa_outward' ,(DCB*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6)))

#for PV diffusion

sim_WM.setPatchSReacK('ring1', 'diff_PV_inward' ,(DPV*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_PV_inward' ,(DPV*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_PV_inward' ,(DPV*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_PV_inward' ,(DPV*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_PV_inward' ,(DPV*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))

sim_WM.setPatchSReacK('ring1', 'diff_PV_outward' ,(DPV*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_PV_outward' ,(DPV*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_PV_outward' ,(DPV*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_PV_outward' ,(DPV*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_PV_outward' ,(DPV*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6)))

#for PVCa diffusion

sim_WM.setPatchSReacK('ring1', 'diff_PVCa_inward' ,(DPV*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_PVCa_inward' ,(DPV*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_PVCa_inward' ,(DPV*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_PVCa_inward' ,(DPV*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_PVCa_inward' ,(DPV*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))

sim_WM.setPatchSReacK('ring1', 'diff_PVCa_outward' ,(DPV*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_PVCa_outward' ,(DPV*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_PVCa_outward' ,(DPV*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_PVCa_outward' ,(DPV*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_PVCa_outward' ,(DPV*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6)))

#for PVMg diffusion

sim_WM.setPatchSReacK('ring1', 'diff_PVMg_inward' ,(DPV*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_PVMg_inward' ,(DPV*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_PVMg_inward' ,(DPV*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_PVMg_inward' ,(DPV*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_PVMg_inward' ,(DPV*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))

sim_WM.setPatchSReacK('ring1', 'diff_PVMg_outward' ,(DPV*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring2', 'diff_PVMg_outward' ,(DPV*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring3', 'diff_PVMg_outward' ,(DPV*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring4', 'diff_PVMg_outward' ,(DPV*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6)))
sim_WM.setPatchSReacK('ring5', 'diff_PVMg_outward' ,(DPV*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6)))



print "Injected ", SK_ro*surfarea, "SK channels"

sim_stoch.setPatchCount('memb_stoch', 'Leak', round(L_ro * surfarea))
print "Injected  ", (L_ro * sim_stoch.getPatchArea('memb_stoch')), "Leak channels"


sim_stoch.setEfieldDT(EF_DT)

sim_stoch.setMembPotential('membrane', init_pot)

sim_stoch.setMembVolRes('membrane', Ra)

#cm = 1.5uF/cm2 -> 1.5e-6F/1e-4m2 ->1.5e-2 F/m2
sim_stoch.setMembCapac('membrane',memb_capac)


#### Recording #####

c=time.ctime()

dc = c.split()[1]+c.split()[2]+'_'+c.split()[3]+'_'+c.split()[4]
dc= dc.replace(':', '_')

try: os.mkdir(root+'data')
except: pass
try: os.mkdir(root+'data/' +  'StochasticCaburst_wellmixed')
except: pass
try: os.mkdir(root+'data/' +  'StochasticCaburst_wellmixed/'+meshfile_ab)
except: pass 

os.mkdir(root+'data/' +  'StochasticCaburst_wellmixed/'+meshfile_ab+'/'+iter_n+'__'+dc )


datfile =  open(root+'data/' +  'StochasticCaburst_wellmixed/'+meshfile_ab+'/'+iter_n+'__'+dc + '/currents.dat', 'w')
datfile2 = open(root+'data/' +  'StochasticCaburst_wellmixed/'+meshfile_ab+'/'+iter_n+'__'+dc + '/voltage.dat', 'w')
datfile3 = open(root+'data/' +  'StochasticCaburst_wellmixed/'+meshfile_ab+'/'+iter_n+'__'+dc + '/calcium.dat', 'w')


for l in range(NTIMEPOINTS):
    print "Tpnt: ", l
    
    # IN this sim V should be constant everywhere
    V = sim_stoch.getTriV(memb_tris[0])
    
    sim_WM.setPatchSReacK('ring0','CaPm0m1', 1.0e3 *3.* alpha_cap(V*1.0e3)*Qt)
    sim_WM.setPatchSReacK('ring0','CaPm1m2', 1.0e3 *2.* alpha_cap(V*1.0e3)*Qt)
    sim_WM.setPatchSReacK('ring0','CaPm2m3', 1.0e3 *1.* alpha_cap(V*1.0e3)*Qt)
    sim_WM.setPatchSReacK('ring0','CaPm3m2', 1.0e3 *3.* beta_cap(V*1.0e3)*Qt)
    sim_WM.setPatchSReacK('ring0','CaPm2m1', 1.0e3 *2.* beta_cap(V*1.0e3)*Qt)
    sim_WM.setPatchSReacK('ring0','CaPm1m0', 1.0e3 *1.* beta_cap(V*1.0e3)*Qt)
    
    sim_WM.setPatchSReacK('ring0', 'CaTm0h0_m1h0', 1.0e3 *2.* alpham_cat(V*1.0e3))
    sim_WM.setPatchSReacK('ring0', 'CaTm1h0_m2h0', 1.0e3 *1.* alpham_cat(V*1.0e3))
    
    sim_WM.setPatchSReacK('ring0', 'CaTm2h0_m1h0', 1.0e3 *2.* betam_cat(V*1.0e3))   
    sim_WM.setPatchSReacK('ring0', 'CaTm1h0_m0h0', 1.0e3 *1.* betam_cat(V*1.0e3)) 
    
    sim_WM.setPatchSReacK('ring0', 'CaTm0h0_m0h1', 1.0e3 *1.* alphah_cat(V*1.0e3))
    sim_WM.setPatchSReacK('ring0', 'CaTm1h0_m1h1', 1.0e3 *1.* alphah_cat(V*1.0e3))
    sim_WM.setPatchSReacK('ring0', 'CaTm2h0_m2h1', 1.0e3 *1.* alphah_cat(V*1.0e3))
    
    sim_WM.setPatchSReacK('ring0', 'CaTm2h1_m2h0', 1.0e3 *1.* betah_cat(V*1.0e3)) 
    sim_WM.setPatchSReacK('ring0', 'CaTm1h1_m1h0', 1.0e3 *1.* betah_cat(V*1.0e3))   
    sim_WM.setPatchSReacK('ring0', 'CaTm0h1_m0h0', 1.0e3 *1.* betah_cat(V*1.0e3)) 
    
    sim_WM.setPatchSReacK('ring0', 'CaTm0h1_m1h1', 1.0e3 *2.* alpham_cat(V*1.0e3))
    sim_WM.setPatchSReacK('ring0', 'CaTm1h1_m2h1', 1.0e3 *1.* alpham_cat(V*1.0e3))
    
    sim_WM.setPatchSReacK('ring0', 'CaTm2h1_m1h1', 1.0e3 *2.* betam_cat(V*1.0e3))   
    sim_WM.setPatchSReacK('ring0', 'CaTm1h1_m0h1', 1.0e3 *1.* betam_cat(V*1.0e3)) 
    
    
    sim_WM.setPatchSReacK('ring0', 'BKC0O0', f_0(V))
    sim_WM.setPatchSReacK('ring0', 'BKC1O1', f_1(V))
    sim_WM.setPatchSReacK('ring0', 'BKC2O2', f_2(V))
    sim_WM.setPatchSReacK('ring0', 'BKC3O3', f_3(V))
    sim_WM.setPatchSReacK('ring0', 'BKC4O4', f_4(V))
    sim_WM.setPatchSReacK('ring0', 'BKO0C0', b_0(V))
    sim_WM.setPatchSReacK('ring0', 'BKO1C1', b_1(V))
    sim_WM.setPatchSReacK('ring0', 'BKO2C2', b_2(V))
    sim_WM.setPatchSReacK('ring0', 'BKO3C3', b_3(V))
    sim_WM.setPatchSReacK('ring0', 'BKO4C4', b_4(V))
    
    sim_WM.run(TIMECONVERTER*l)

    # Now do the communication between the sims
    
    Si = sim_WM.getCompConc('shell0', 'Ca')

    So = Ca_oconc
    
    # Get the single-channel currents first
    tcur_CaP_sc = cf.getGHKI(CaP_P, V, 2, TEMPERATURE+273.15, Si*1.0e3, So*1.0e3)
    tcur_CaT_sc = cf.getGHKI(CaT_P, V, 2, TEMPERATURE+273.15, Si*1.0e3, So*1.0e3)
    tcur_BK_sc = cf.getOhmI(V, BK_rev, BK_G)
    
    tcur_SK_sc = cf.getOhmI(V, SK_rev, SK_G)
    
    tcur_CaP = tcur_CaP_sc*(sim_WM.getPatchCount('ring0', 'CaP_m3'))
    # alpha is to h1
    tcur_CaT = tcur_CaT_sc*(sim_WM.getPatchCount('ring0', 'CaT_m2h1'))
    tcur_BK = tcur_BK_sc*((sim_WM.getPatchCount('ring0', 'BK_O0')+\
                             sim_WM.getPatchCount('ring0', 'BK_O1')+\
                                sim_WM.getPatchCount('ring0', 'BK_O2')+\
                                  sim_WM.getPatchCount('ring0', 'BK_O3')+\
                                    sim_WM.getPatchCount('ring0', 'BK_O4')))
    tcur_SK = tcur_SK_sc*((sim_WM.getPatchCount('ring0', 'SK_O1')+\
                              sim_WM.getPatchCount('ring0', 'SK_O2')))

    nmtris = len(memb_tris)
    
    for t in memb_tris: sim_stoch.setTriIClamp(t, ((tcur_CaP+tcur_CaT+tcur_BK+tcur_SK)/nmtris))
    
    sim_stoch.run(TIMECONVERTER*l)
    
    ca_count_inj =  -1.0*((tcur_CaP+tcur_CaT)*TIMECONVERTER)/(2*E_CHARGE)
    actual_inj= 0.0
    
    t_count = sim_WM.getCompCount('shell0', 'Ca')
    sim_WM.setCompCount('shell0', 'Ca', t_count + ca_count_inj)
    
    ca_shell = sim_WM.getCompCount('shell0', 'Ca')
    
    datfile.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile.write('%.6g' %((tcur_CaP*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_CaT*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_BK*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_SK*1.0e-1)/surfarea) + ' ')  
    datfile.write('\n')

    datfile2.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile2.write('%.6g' %(sim_stoch.getTetV(cent_tet)*1.0e3) + ' ')
    datfile2.write('\n')
    
    datfile3.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile3.write('%.6g' %sim_WM.getCompCount('shell0', 'Ca') +' ')
    datfile3.write('%.6g' %(Si*1.0e6)+ ' ')
    datfile3.write('\n')

datfile.close()
datfile2.close()
datfile3.close()


