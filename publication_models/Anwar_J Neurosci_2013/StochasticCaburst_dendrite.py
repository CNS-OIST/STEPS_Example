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
# *StochasticCaburst_dendrite.py : The stochastic calcium burst 
# model, on the realistic dendrite morphology. Each dendritic branch is
# modeled with well-mixed cytosolic compartments.
#
# Script authors: Haroon Anwar and Iain Hepburn
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# USAGE
#
# $ python StochasticCaburst_dendrite.py *root* *iter_n* 
#
#  *root* is the path to the location for data storage 
#  *iter_n* (is intened to be an integer) is an identifier number for each 
#     simulation iteration
#
# E.g: python StochasticCaburst_dendrite.py ~/stochcasims/ 1
#
#
# OUTPUT 
#
# In (root)/data/StochasticCaburst_dendrite/chop_mesh.inp/(iter_n+time) directory 
# 3 data files will be recorded. Each file contains one row for every 
# time-point at which data is recorded, organised into the following columns:
# 
# currents.dat
# Time (ms), (for every compartment): P-type current, T-type current, BK current, SK current
# (current units are Amps/m^2)
#
# voltage.dat
# Time (ms), voltage at mesh centre (mV)
#
# calcium.dat
# Time (ms), (for every compartment): number of calcium ions in submembrane shell, 
#  calcium concentration in submembrane shell (micromolar)
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

root, iter_n =  sys.argv[1], sys.argv[2]

cp_times = [0.0, 0.1, 0.2, 0.3, 0.4]

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
# # # # # # # /Users/anwar/Copy/ModelDB_scripts/extra/geom_info_STEPS.dat# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

########### MESH & COMPARTMENTALIZATION #################

##########Import Mesh for EField calculation

# For stochastic sim:

meshfile_ab = 'chop_mesh.inp'
mesh_stoch, nodep, tetp, trip = meshio.importAbaqus('./meshes/'+meshfile_ab,1e-6)
tetgroups = tetp.blocksToGroups()

inner_tets = range(mesh_stoch.ntets)

print inner_tets.__len__(), " tets in inner compartment"


#Define geometrical constants for all (100) compartments with concentric shells

geom_info_file = file("./extra/geom_info_STEPS.dat","r")

geom_info = geom_info_file.read()

LengthNdiam = geom_info.split()

CUBIT_IDs = []
Length = []
Diam = []

count = 0

for i in LengthNdiam:
    if count==0:
        CUBIT_IDs.append(int(i))
        count = count+1
    elif count==1:
        Length.append(float(i)*1e-6)
        count=count+1
    elif count==2:
        Diam.append(float(i)*1e-6)
        count = 0

Comp_tetIDs = [None]*len(Length)

for i in range(len(Length)):
    Comp_tetIDs[i] = tetgroups['EB'+str(int(CUBIT_IDs[i]))]
    print len(Comp_tetIDs[i])


Nannulis = []
Shells_radii = [None]*len(Length)
Shells_vols = [None]*len(Length)
Rings_areas = [None]*len(Length)
DScales_ShellIn = [None]*len(Length)
DScales_ShellOut = [None]*len(Length)

for i in range(len(Length)):
    Shells_radii[i] = []
    Shells_vols[i] = []
    Rings_areas[i] = []
    DScales_ShellIn[i] = []
    DScales_ShellOut[i] = []
    Nannulis.append(round((Diam[i]/0.4e-6)+1))
    if Nannulis[i]-((Diam[i]/0.4e-6)+1)==0.5:
        Nannulis[i]=Nannulis[i]-1
    for j in range(int(Nannulis[i])):
        if j==0:
            Shells_radii[i].append(Diam[i]/2)
        elif j==1:
            Shells_radii[i].append((Diam[i]/2)-0.1e-6)
        else:
            Shells_radii[i].append(Shells_radii[i][j-1]-0.2e-6)
        Rings_areas[i].append(2*math.pi*Shells_radii[i][j]*Length[i])

    for j in range(int(Nannulis[i])):
        if j==int(Nannulis[i])-1:
            Shells_vols[i].append(math.pi*((Shells_radii[i][j]**2))*Length[i])
        else:
            Shells_vols[i].append(math.pi*((Shells_radii[i][j]**2) - (Shells_radii[i][j+1]**2))*Length[i])

    for j in range(int(Nannulis[i])-1):
        if j==int(Nannulis[i]-2):
            DScales_ShellOut[i].append(2*Shells_radii[i][j+1]/((Shells_radii[i][j+1]**2)*(0.1e-6+Shells_radii[i][j+1])))
            DScales_ShellIn[i].append(2*Shells_radii[i][j+1]/((Shells_radii[i][j]**2 - Shells_radii[i][j+1]**2)*(0.1e-6+Shells_radii[i][j+1])))
        else:
            DScales_ShellOut[i].append(2*Shells_radii[i][j+1]/((Shells_radii[i][j+1]**2 - Shells_radii[i][j+2]**2)*0.2e-6))
            DScales_ShellIn[i].append(2*Shells_radii[i][j+1]/((Shells_radii[i][j]**2 - Shells_radii[i][j+1]**2)*0.2e-6))

shells = [None]*len(Length)
rings = [None]*len(Length)

#Geometry container object:

wmgeom = sgeom.Geom()

for i in range(len(Length)):
    shells[i]=[None]*int(Nannulis[i])
    rings[i]=[None]*int(Nannulis[i])
    for j in range(int(Nannulis[i])):
        shells[i][j] = sgeom.Comp('shells'+str(i)+str(j), wmgeom, vol=Shells_vols[i][j])
        if j==0:
            rings[i][j] = sgeom.Patch('rings'+str(i)+str(j), wmgeom, icomp=shells[i][j], area=Rings_areas[i][j])
        else:
            rings[i][j] = sgeom.Patch('rings'+str(i)+str(j), wmgeom, ocomp=shells[i][j-1], icomp=shells[i][j], area=Rings_areas[i][j])


#For all patches, add the surface system

for i in range(len(Length)):
    for j in range(int(Nannulis[i])):
        if j==0:
            rings[i][j].addSurfsys('ssys_chans')
        else:
            rings[i][j].addSurfsys('ssys_diff')

        shells[i][j].addVolsys('vsys_buffs')


########## Create an intracellular compartment i.e. cytosolic compartment

cyto_stoch = sgeom.TmComp('cyto_stoch', mesh_stoch, inner_tets)


########## Finding the triangles comprising the memberane

memb_tris = list(mesh_stoch.getSurfTris())

Compborder_triIDs = [None]*len(Length)
Compborder_tri_areas = [0.0]*len(Length)

for i in range(len(Length)):
    Compborder_triIDs[i]=[]
    for j in range(len(Comp_tetIDs[i])):
        tritemp = mesh_stoch.getTetTriNeighb(Comp_tetIDs[i][j])
        for tri in tritemp:
            if tri in memb_tris:
                Compborder_triIDs[i].append(tri)
                Compborder_tri_areas[i] = Compborder_tri_areas[i] + mesh_stoch.getTriArea(tri)
                break

########## Create a membrane as a surface mesh

# Stochastic sim:
memb_stoch = sgeom.TmPatch('memb_stoch', mesh_stoch, memb_tris, cyto_stoch)
memb_stoch.addSurfsys('ssys_stoch')

# For EField calculation
print "Creating membrane.."
membrane = sgeom.Memb('membrane', mesh_stoch, [memb_stoch], opt_file_name = './meshes/'+meshfile_ab+"_optimalidx")
print "Membrane created."

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

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

for i in range(len(Length)):
    for j in range(int(Nannulis[i])):
        sim_WM.setCompConc('shells'+str(i)+str(j), 'Ca', Ca_iconc)
        sim_WM.setCompConc('shells'+str(i)+str(j), 'Mg', Mg_conc)

        sim_WM.setCompConc('shells'+str(i)+str(j), 'iCBsf',    iCBsf_conc)
        sim_WM.setCompConc('shells'+str(i)+str(j), 'iCBsCa',   iCBsCa_conc)
        sim_WM.setCompConc('shells'+str(i)+str(j), 'iCBCaf',   iCBCaf_conc)
        sim_WM.setCompConc('shells'+str(i)+str(j), 'iCBCaCa',  iCBCaCa_conc)

        sim_WM.setCompConc('shells'+str(i)+str(j), 'CBsf',     CBsf_conc)
        sim_WM.setCompConc('shells'+str(i)+str(j), 'CBsCa',    CBsCa_conc)
        sim_WM.setCompConc('shells'+str(i)+str(j), 'CBCaf',    CBCaf_conc)
        sim_WM.setCompConc('shells'+str(i)+str(j), 'CBCaCa',   CBCaCa_conc)

        sim_WM.setCompConc('shells'+str(i)+str(j), 'PV',       PV_conc)
        sim_WM.setCompConc('shells'+str(i)+str(j), 'PVCa',     PVCa_conc)
        sim_WM.setCompConc('shells'+str(i)+str(j), 'PVMg',     PVMg_conc)
    
        if j==0:
            surfarea = sim_WM.getPatchArea('rings'+str(i)+str(j))
            pumpnbs = 6.022141e12*surfarea

            sim_WM.setPatchCount('rings'+str(i)+str(j), 'Pump', pumpnbs)
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'CaPump', 0)

            # CaP
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'CaP_m0' , round(CaP_ro*surfarea*CaP_m0_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'CaP_m1' , round(CaP_ro*surfarea*CaP_m1_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'CaP_m2' , round(CaP_ro*surfarea*CaP_m2_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'CaP_m3' , round(CaP_ro*surfarea*CaP_m3_p))

            print "Injected  ", CaP_ro*surfarea, "CaP channels"

            # CaT

            # From cstate: CaT_m2h0 conducting
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'CaT_m0h0' , round(CaT_ro*surfarea*CaT_m0h0_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'CaT_m1h0' , round(CaT_ro*surfarea*CaT_m1h0_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'CaT_m2h0' , round(CaT_ro*surfarea*CaT_m2h0_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'CaT_m0h1' , round(CaT_ro*surfarea*CaT_m0h1_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'CaT_m1h1' , round(CaT_ro*surfarea*CaT_m1h1_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'CaT_m2h1' , round(CaT_ro*surfarea*CaT_m2h1_p))

            print "Injected  ", CaT_ro*surfarea, "CaT channels"

            # BK
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'BK_C0' , round(BK_ro*surfarea*BK_C0_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'BK_C1' , round(BK_ro*surfarea*BK_C1_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'BK_C2' , round(BK_ro*surfarea*BK_C2_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'BK_C3' , round(BK_ro*surfarea*BK_C3_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'BK_C4' , round(BK_ro*surfarea*BK_C4_p))

            sim_WM.setPatchCount('rings'+str(i)+str(j), 'BK_O0' , round(BK_ro*surfarea*BK_O0_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'BK_O1' , round(BK_ro*surfarea*BK_O1_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'BK_O2' , round(BK_ro*surfarea*BK_O2_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'BK_O3' , round(BK_ro*surfarea*BK_O3_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'BK_O4' , round(BK_ro*surfarea*BK_O4_p))


            print "Injected  ", BK_ro*surfarea, "BK channels"

            # SK
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'SK_C1' , round(SK_ro*surfarea*SK_C1_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'SK_C2' , round(SK_ro*surfarea*SK_C2_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'SK_C3' , round(SK_ro*surfarea*SK_C3_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'SK_C4' , round(SK_ro*surfarea*SK_C4_p))

            sim_WM.setPatchCount('rings'+str(i)+str(j), 'SK_O1' , round(SK_ro*surfarea*SK_O1_p))
            sim_WM.setPatchCount('rings'+str(i)+str(j), 'SK_O2' , round(SK_ro*surfarea*SK_O2_p))
            
            print "Injected ", SK_ro*surfarea, "SK channels"
        
        else:
            #set the rate constants for diffusion (Diffusion is modeled as surface reaction here)

            #for Ca diffusion

            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_Ca_inward' ,(DCST*DScales_ShellIn[i][j-1]))
            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_Ca_outward' ,(DCST*DScales_ShellOut[i][j-1]))

            #for CBsf diffusin
            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_CBsf_inward' ,(DCB*DScales_ShellIn[i][j-1]))
            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_CBsf_outward' ,(DCB*DScales_ShellOut[i][j-1]))

            #for CBsCa diffusion

            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_CBsCa_inward' ,(DCB*DScales_ShellIn[i][j-1]))
            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_CBsCa_outward' ,(DCB*DScales_ShellOut[i][j-1]))

            #for CBCaf diffusion

            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_CBCaf_inward' ,(DCB*DScales_ShellIn[i][j-1]))
            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_CBCaf_outward' ,(DCB*DScales_ShellOut[i][j-1]))

            #for CBCaCa diffusion

            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_CBCaCa_inward' ,(DCB*DScales_ShellIn[i][j-1]))
            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_CBCaCa_outward' ,(DCB*DScales_ShellOut[i][j-1]))

            #for PV diffusion

            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_PV_inward' ,(DPV*DScales_ShellIn[i][j-1]))
            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_PV_outward' ,(DPV*DScales_ShellOut[i][j-1]))

            #for PVCa diffusion

            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_PVCa_inward' ,(DPV*DScales_ShellIn[i][j-1]))
            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_PVCa_outward' ,(DPV*DScales_ShellOut[i][j-1]))

            #for PVMg diffusion

            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_PVMg_inward' ,(DPV*DScales_ShellIn[i][j-1]))
            sim_WM.setPatchSReacK('rings'+str(i)+str(j), 'diff_PVMg_outward' ,(DPV*DScales_ShellOut[i][j-1]))


surfarea = sim_stoch.getPatchArea('memb_stoch')
#Compute the surface area for full mesh rather than individual compartments
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
try: os.mkdir(root+'data/' +  'StochasticCaburst_dendrite')
except: pass
try: os.mkdir(root+'data/' +  'StochasticCaburst_dendrite/'+meshfile_ab)
except: pass 

os.mkdir(root+'data/' +  'StochasticCaburst_dendrite/'+meshfile_ab+'/'+iter_n+'__'+dc )


datfile =  open(root+'data/' +  'StochasticCaburst_dendrite/'+meshfile_ab+'/'+iter_n+'__'+dc + '/currents.dat', 'w')
datfile2 = open(root+'data/' +  'StochasticCaburst_dendrite/'+meshfile_ab+'/'+iter_n+'__'+dc + '/voltage.dat', 'w')
datfile3 = open(root+'data/' +  'StochasticCaburst_dendrite/'+meshfile_ab+'/'+iter_n+'__'+dc + '/calcium.dat', 'w')

btime = time.time()
for l in range(NTIMEPOINTS):
    print "Tpnt: ", l

    if TIMECONVERTER*l in cp_times:
        sim_stoch.checkpoint(root+'data/' +  'StochasticCaburst_dendrite/'+meshfile_ab+'/'+iter_n+'__'+dc + '/checkpoint__sim_stoch__'+  str(TIMECONVERTER*l))
        sim_WM.checkpoint(root+'data/' +  'StochasticCaburst_dendrite/'+meshfile_ab+'/'+iter_n+'__'+dc + '/checkpoint__sim_WM__'+  str(TIMECONVERTER*l))
        datfile.flush()
        datfile2.flush()
        datfile3.flush()
    
    datfile.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile2.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile3.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    # IN this sim V should be constant everywhere in the compartment
    for i in range(len(Length)):
        V = sim_stoch.getTriV(Compborder_triIDs[i][0])
    
        sim_WM.setPatchSReacK('rings'+str(i)+str(0),'CaPm0m1', 1.0e3 *3.* alpha_cap(V*1.0e3)*Qt)
        sim_WM.setPatchSReacK('rings'+str(i)+str(0),'CaPm1m2', 1.0e3 *2.* alpha_cap(V*1.0e3)*Qt)
        sim_WM.setPatchSReacK('rings'+str(i)+str(0),'CaPm2m3', 1.0e3 *1.* alpha_cap(V*1.0e3)*Qt)
        sim_WM.setPatchSReacK('rings'+str(i)+str(0),'CaPm3m2', 1.0e3 *3.* beta_cap(V*1.0e3)*Qt)
        sim_WM.setPatchSReacK('rings'+str(i)+str(0),'CaPm2m1', 1.0e3 *2.* beta_cap(V*1.0e3)*Qt)
        sim_WM.setPatchSReacK('rings'+str(i)+str(0),'CaPm1m0', 1.0e3 *1.* beta_cap(V*1.0e3)*Qt)
    
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'CaTm0h0_m1h0', 1.0e3 *2.* alpham_cat(V*1.0e3))
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'CaTm1h0_m2h0', 1.0e3 *1.* alpham_cat(V*1.0e3))
    
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'CaTm2h0_m1h0', 1.0e3 *2.* betam_cat(V*1.0e3))   
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'CaTm1h0_m0h0', 1.0e3 *1.* betam_cat(V*1.0e3)) 
    
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'CaTm0h0_m0h1', 1.0e3 *1.* alphah_cat(V*1.0e3))
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'CaTm1h0_m1h1', 1.0e3 *1.* alphah_cat(V*1.0e3))
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'CaTm2h0_m2h1', 1.0e3 *1.* alphah_cat(V*1.0e3))
    
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'CaTm2h1_m2h0', 1.0e3 *1.* betah_cat(V*1.0e3)) 
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'CaTm1h1_m1h0', 1.0e3 *1.* betah_cat(V*1.0e3))   
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'CaTm0h1_m0h0', 1.0e3 *1.* betah_cat(V*1.0e3)) 
    
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'CaTm0h1_m1h1', 1.0e3 *2.* alpham_cat(V*1.0e3))
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'CaTm1h1_m2h1', 1.0e3 *1.* alpham_cat(V*1.0e3))
    
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'CaTm2h1_m1h1', 1.0e3 *2.* betam_cat(V*1.0e3))   
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'CaTm1h1_m0h1', 1.0e3 *1.* betam_cat(V*1.0e3)) 
    
    
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'BKC0O0', f_0(V))
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'BKC1O1', f_1(V))
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'BKC2O2', f_2(V))
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'BKC3O3', f_3(V))
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'BKC4O4', f_4(V))
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'BKO0C0', b_0(V))
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'BKO1C1', b_1(V))
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'BKO2C2', b_2(V))
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'BKO3C3', b_3(V))
        sim_WM.setPatchSReacK('rings'+str(i)+str(0), 'BKO4C4', b_4(V))
    
    sim_WM.run(TIMECONVERTER*l)

    # Now do the communication between the sims
    for i in range(len(Length)):
        V = sim_stoch.getTriV(Compborder_triIDs[i][0])
        Si = sim_WM.getCompConc('shells'+str(i)+str(0), 'Ca')

        So = Ca_oconc
    
        # Get the single-channel currents first
        tcur_CaP_sc = cf.getGHKI(CaP_P, V, 2, TEMPERATURE+273.15, Si*1.0e3, So*1.0e3)
        tcur_CaT_sc = cf.getGHKI(CaT_P, V, 2, TEMPERATURE+273.15, Si*1.0e3, So*1.0e3)
        tcur_BK_sc = cf.getOhmI(V, BK_rev, BK_G)
    
        tcur_SK_sc = cf.getOhmI(V, SK_rev, SK_G)
    
        tcur_CaP = tcur_CaP_sc*(sim_WM.getPatchCount('rings'+str(i)+str(0), 'CaP_m3'))
        # alpha is to h1
        tcur_CaT = tcur_CaT_sc*(sim_WM.getPatchCount('rings'+str(i)+str(0), 'CaT_m2h1'))
        tcur_BK = tcur_BK_sc*((sim_WM.getPatchCount('rings'+str(i)+str(0), 'BK_O0')+\
                                  sim_WM.getPatchCount('rings'+str(i)+str(0), 'BK_O1')+\
                                     sim_WM.getPatchCount('rings'+str(i)+str(0), 'BK_O2')+\
                                        sim_WM.getPatchCount('rings'+str(i)+str(0), 'BK_O3')+\
                                           sim_WM.getPatchCount('rings'+str(i)+str(0), 'BK_O4')))
        tcur_SK = tcur_SK_sc*((sim_WM.getPatchCount('rings'+str(i)+str(0), 'SK_O1')+\
                                  sim_WM.getPatchCount('rings'+str(i)+str(0), 'SK_O2')))

        AreaRatio = Compborder_tri_areas[i]/Rings_areas[i][0]
        nmtris = len(Compborder_triIDs[i])
    
        for tri in Compborder_triIDs[i]: sim_stoch.setTriIClamp(tri, (((tcur_CaP+tcur_CaT+tcur_BK+tcur_SK)*AreaRatio)/nmtris))

        ca_count_inj =  -1.0*((tcur_CaP+tcur_CaT)*TIMECONVERTER)/(2*E_CHARGE)

        t_count = sim_WM.getCompCount('shells'+str(i)+str(0), 'Ca')
        sim_WM.setCompCount('shells'+str(i)+str(0), 'Ca', t_count + ca_count_inj)

        datfile.write('%.6g' %((tcur_CaP*1.0e-1)/Compborder_tri_areas[i]) + ' ')
        datfile.write('%.6g' %((tcur_CaT*1.0e-1)/Compborder_tri_areas[i]) + ' ')
        datfile.write('%.6g' %((tcur_BK*1.0e-1)/Compborder_tri_areas[i]) + ' ')
        datfile.write('%.6g' %((tcur_SK*1.0e-1)/Compborder_tri_areas[i]) + ' ')

        datfile3.write('%.6g' %(sim_WM.getCompCount('shells'+str(i)+str(0), 'Ca')) +' ')
        datfile3.write('%.6g' %(sim_WM.getCompConc('shells'+str(i)+str(0), 'Ca')*1.0e6)+ ' ')

    datfile.write('\n')
    datfile3.write('\n')
    
    sim_stoch.run(TIMECONVERTER*l)
    
    for i in range(len(Length)):
        V = sim_stoch.getTriV(Compborder_triIDs[i][0])
        datfile2.write('%.6g' %(V*1.0e3) + ' ')
    datfile2.write('\n')


datfile.close()
datfile2.close()
datfile3.close()

