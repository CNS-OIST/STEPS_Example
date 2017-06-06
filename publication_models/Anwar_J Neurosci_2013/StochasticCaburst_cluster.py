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
# *StochasticCaburst_cluster.py : The stochastic calcium burst model with 
# P-type calcium channels clustered around BK channels
#
# Script authors: Haroon Anwar and Iain Hepburn
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# USAGE
#
# $ python StochasticCaburst_cluster.py *mesh* *root* *iter_n* *clusterSize*
#  
#  *mesh* is the tetrahedral mesh (10um to 160um cylinder)
#  *root* is the path to the location for data storage 
#  *iter_n* (is intened to be an integer) is an identifier number for each 
#     simulation iteration. 
#  *clusterSize* is the size of the P-type channel clusters
#
# E.g: 
# $ python StochasticCaburst_cluster.py Cylinder2_dia2um_L10um_outer0_3um_0.3shell_0.3size_19156tets_adaptive.inp ~/stochcasims/ 1 4
#
#
# OUTPUT 
#
# In (root)/data/StochasticCaburst_cluster/(mesh)/(iter_n+time) directory 
# 5 data files will be recorded. Each file contains one row for every 
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
# Time (ms), calcium concentration in submembrane (micromolar),
# number of calcium ions in submembrane. 
#
# OpenBKandCa.dat
# Time (ms), (for every BK triangle):  BK_O0, BK_O1, BK_O2, BK_O3, BK_O4, number
# of calcium ions in inner tetrahedron, calcium concentration in inner tetrahedron.   
#
# ChannelsDistribution.dat 
# The channel distributions. Please see script for details. 
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import math
import time
import random
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

meshfile_ab, root, iter_n, clusterSize = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

if meshfile_ab == 'Cylinder2_dia2um_L160um_outer0_0.3shell_0.3size_279152tets_adaptive.inp': cyl160=True
else: cyl160=False

clusterSize = int(clusterSize)

########################### BIOCHEMICAL MODEL ###############################

# Two models required: Stochastic and deterministic
mdl = smodel.Model()

# Calcium
Ca = smodel.Spec('Ca', mdl)
Ca.setValence(2)


# Pump
Pump = smodel.Spec('Pump', mdl)
# CaPump
CaPump = smodel.Spec('CaPump', mdl)

# iCBsf
iCBsf = smodel.Spec('iCBsf', mdl)
# iCBsCa
iCBsCa = smodel.Spec('iCBsCa', mdl)
# iCBCaf
iCBCaf = smodel.Spec('iCBCaf', mdl)
# iCBCaCa
iCBCaCa = smodel.Spec('iCBCaCa', mdl)

# CBsf
CBsf = smodel.Spec('CBsf', mdl)
# CBsCa
CBsCa = smodel.Spec('CBsCa', mdl)
# CBCaf
CBCaf = smodel.Spec('CBCaf', mdl)
# CBCaCa
CBCaCa = smodel.Spec('CBCaCa', mdl)

# PV
PV = smodel.Spec('PV', mdl)
# PVMg
PVMg = smodel.Spec('PVMg', mdl)
# PVCa
PVCa = smodel.Spec('PVCa', mdl)
# Mg
Mg = smodel.Spec('Mg', mdl)

# Vol/surface systems
vsys = smodel.Volsys('vsys', mdl)
ssys = smodel.Surfsys('ssys', mdl)


# Diffusions 
diff_Ca = smodel.Diff('diff_Ca', vsys, Ca)
diff_Ca.setDcst(DCST)
diff_CBsf = smodel.Diff('diff_CBsf', vsys, CBsf)
diff_CBsf.setDcst(DCB)
diff_CBsCa = smodel.Diff('diff_CBsCa', vsys, CBsCa)
diff_CBsCa.setDcst(DCB)
diff_CBCaf = smodel.Diff('diff_CBCaf', vsys, CBCaf)
diff_CBCaf.setDcst(DCB)
diff_CBCaCa = smodel.Diff('diff_CBCaCa', vsys, CBCaCa)
diff_CBCaCa.setDcst(DCB)
diff_PV = smodel.Diff('diff_PV', vsys, PV)
diff_PV.setDcst(DPV)
diff_PVCa = smodel.Diff('diff_PVCa', vsys, PVCa)
diff_PVCa.setDcst(DPV)
diff_PVMg = smodel.Diff('diff_PVMg', vsys, PVMg)
diff_PVMg.setDcst(DPV)


#Pump
PumpD_f = smodel.SReac('PumpD_f', ssys, ilhs=[Ca], slhs=[Pump], srhs=[CaPump])
PumpD_f.setKcst(P_f_kcst)

PumpD_b = smodel.SReac('PumpD_b', ssys, slhs=[CaPump], irhs=[Ca], srhs=[Pump])
PumpD_b.setKcst(P_b_kcst)

PumpD_k = smodel.SReac('PumpD_k', ssys, slhs=[CaPump], srhs=[Pump])
PumpD_k.setKcst(P_k_kcst)

#iCBsf-fast
iCBsf1_f = smodel.Reac('iCBsf1_f', vsys, lhs=[Ca,iCBsf], rhs=[iCBsCa], kcst = iCBsf1_f_kcst)
iCBsf1_b = smodel.Reac('iCBsf1_b', vsys, lhs=[iCBsCa], rhs=[Ca, iCBsf], kcst = iCBsf1_b_kcst)

#iCBsCa
iCBsCa_f = smodel.Reac('iCBsCa_f', vsys, lhs=[Ca,iCBsCa], rhs=[iCBCaCa], kcst = iCBsCa_f_kcst)
iCBsCa_b = smodel.Reac('iCBsCa_b', vsys, lhs=[iCBCaCa], rhs=[Ca,iCBsCa], kcst = iCBsCa_b_kcst)

#iCBsf_slow
iCBsf2_f = smodel.Reac('iCBsf2_f', vsys, lhs=[Ca,iCBsf], rhs=[iCBCaf], kcst = iCBsf2_f_kcst)
iCBsf2_b = smodel.Reac('iCBsf2_b', vsys, lhs=[iCBCaf], rhs=[Ca,iCBsf], kcst = iCBsf2_b_kcst)

#iCBCaf
iCBCaf_f = smodel.Reac('iCBCaf_f', vsys, lhs=[Ca,iCBCaf], rhs=[iCBCaCa], kcst = iCBCaf_f_kcst)
iCBCaf_b = smodel.Reac('iCBCaf_b', vsys, lhs=[iCBCaCa], rhs=[Ca,iCBCaf], kcst = iCBCaf_b_kcst)

#CBsf-fast
CBsf1_f = smodel.Reac('CBsf1_f', vsys, lhs=[Ca,CBsf], rhs=[CBsCa], kcst = CBsf1_f_kcst)
CBsf1_b = smodel.Reac('CBsf1_b', vsys, lhs=[CBsCa], rhs=[Ca,CBsf], kcst = CBsf1_b_kcst)

#CBsCa
CBsCa_f = smodel.Reac('CBsCa_f', vsys, lhs=[Ca,CBsCa], rhs=[CBCaCa], kcst = CBsCa_f_kcst)
CBsCa_b = smodel.Reac('CBsCa_b', vsys, lhs=[CBCaCa], rhs=[Ca,CBsCa], kcst = CBsCa_b_kcst)

#CBsf_slow
CBsf2_f = smodel.Reac('CBsf2_f', vsys, lhs=[Ca,CBsf], rhs=[CBCaf], kcst = CBsf2_f_kcst)
CBsf2_b = smodel.Reac('CBsf2_b', vsys, lhs=[CBCaf], rhs=[Ca,CBsf], kcst = CBsf2_b_kcst)

#CBCaf
CBCaf_f = smodel.Reac('CBCaf_f', vsys, lhs=[Ca,CBCaf], rhs=[CBCaCa], kcst = CBCaf_f_kcst)
CBCaf_b = smodel.Reac('CBCaf_b', vsys, lhs=[CBCaCa], rhs=[Ca,CBCaf], kcst = CBCaf_b_kcst)

#PVca
PVca_f = smodel.Reac('PVca_f', vsys, lhs=[Ca,PV], rhs=[PVCa], kcst = PVca_f_kcst)
PVca_b = smodel.Reac('PVca_b', vsys, lhs=[PVCa], rhs=[Ca,PV], kcst = PVca_b_kcst)

#PVmg
PVmg_f = smodel.Reac('PVmg_f', vsys, lhs=[Mg,PV], rhs=[PVMg], kcst = PVmg_f_kcst)
PVmg_b = smodel.Reac('PVmg_b', vsys, lhs=[PVMg], rhs=[Mg,PV], kcst = PVmg_b_kcst)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # CHANNELS  # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###### CaP channel ############## 

CaPchan = smodel.Chan('CaPchan', mdl)

CaP_m0 = smodel.ChanState('CaP_m0', mdl, CaPchan)
CaP_m1 = smodel.ChanState('CaP_m1', mdl, CaPchan)
CaP_m2 = smodel.ChanState('CaP_m2', mdl, CaPchan)
CaP_m3 = smodel.ChanState('CaP_m3', mdl, CaPchan)


CaPm0m1 = smodel.VDepSReac('CaPm0m1', ssys, slhs = [CaP_m0], srhs = [CaP_m1], k= lambda V: 1.0e3 *3.* alpha_cap(V*1.0e3)* Qt)
CaPm1m2 = smodel.VDepSReac('CaPm1m2', ssys, slhs = [CaP_m1], srhs = [CaP_m2], k= lambda V: 1.0e3 *2.* alpha_cap(V*1.0e3)* Qt)
CaPm2m3 = smodel.VDepSReac('CaPm2m3', ssys, slhs = [CaP_m2], srhs = [CaP_m3], k= lambda V: 1.0e3 *1.* alpha_cap(V*1.0e3)* Qt)

CaPm3m2 = smodel.VDepSReac('CaPm3m2', ssys, slhs = [CaP_m3], srhs = [CaP_m2], k= lambda V: 1.0e3 *3.* beta_cap(V*1.0e3)* Qt)
CaPm2m1 = smodel.VDepSReac('CaPm2m1', ssys, slhs = [CaP_m2], srhs = [CaP_m1], k= lambda V: 1.0e3 *2.* beta_cap(V*1.0e3)* Qt)
CaPm1m0 = smodel.VDepSReac('CaPm1m0', ssys, slhs = [CaP_m1], srhs = [CaP_m0], k= lambda V: 1.0e3 *1.* beta_cap(V*1.0e3)* Qt)

if cyl160:
    OC_CaP = smodel.GHKcurr('OC_CaP', ssys, CaP_m3, Ca, virtual_oconc = Ca_oconc, computeflux = True)
else:
    OC_CaP = smodel.GHKcurr('OC_CaP', ssys, CaP_m3, Ca, computeflux = True)

OC_CaP.setP(CaP_P)

######## CaT channel ##########  

CaTchan = smodel.Chan('CaTchan', mdl)

CaT_m0h0 = smodel.ChanState('CaT_m0h0', mdl, CaTchan)
CaT_m0h1 = smodel.ChanState('CaT_m0h1', mdl, CaTchan)
CaT_m1h0 = smodel.ChanState('CaT_m1h0', mdl, CaTchan)
CaT_m1h1 = smodel.ChanState('CaT_m1h1', mdl, CaTchan)
CaT_m2h0 = smodel.ChanState('CaT_m2h0', mdl, CaTchan)
CaT_m2h1 = smodel.ChanState('CaT_m2h1', mdl, CaTchan)

CaTm0h0_m1h0 = smodel.VDepSReac('CaTm0h0_m1h0', ssys, slhs = [CaT_m0h0], srhs = [CaT_m1h0], k= lambda V: 1.0e3 *2.* alpham_cat(V*1.0e3))
CaTm1h0_m2h0 = smodel.VDepSReac('CaTm1h0_m2h0', ssys, slhs = [CaT_m1h0], srhs = [CaT_m2h0], k= lambda V: 1.0e3 *1.* alpham_cat(V*1.0e3))

CaTm2h0_m1h0 = smodel.VDepSReac('CaTm2h0_m1h0', ssys, slhs = [CaT_m2h0], srhs = [CaT_m1h0], k= lambda V: 1.0e3 *2.* betam_cat(V*1.0e3))
CaTm1h0_m0h0 = smodel.VDepSReac('CaTm1h0_m0h0', ssys, slhs = [CaT_m1h0], srhs = [CaT_m0h0], k= lambda V: 1.0e3 *1.* betam_cat(V*1.0e3))

CaTm0h1_m1h1 = smodel.VDepSReac('CaTm0h1_m1h1', ssys, slhs = [CaT_m0h1], srhs = [CaT_m1h1], k= lambda V: 1.0e3 *2.* alpham_cat(V*1.0e3))
CaTm1h1_m2h1 = smodel.VDepSReac('CaTm1h1_m2h1', ssys, slhs = [CaT_m1h1], srhs = [CaT_m2h1], k= lambda V: 1.0e3 *1.* alpham_cat(V*1.0e3))

CaTm2h1_m1h1 = smodel.VDepSReac('CaTm2h1_m1h1', ssys, slhs = [CaT_m2h1], srhs = [CaT_m1h1], k= lambda V: 1.0e3 *2.* betam_cat(V*1.0e3))
CaTm1h1_m0h1 = smodel.VDepSReac('CaTm1h1_m0h1', ssys, slhs = [CaT_m1h1], srhs = [CaT_m0h1], k= lambda V: 1.0e3 *1.* betam_cat(V*1.0e3))

CaTm0h0_m0h1 = smodel.VDepSReac('CaTm0h0_m0h1', ssys, slhs = [CaT_m0h0], srhs = [CaT_m0h1], k= lambda V: 1.0e3 *1.* alphah_cat(V*1.0e3))
CaTm1h0_m1h1 = smodel.VDepSReac('CaTm1h0_m1h1', ssys, slhs = [CaT_m1h0], srhs = [CaT_m1h1], k= lambda V: 1.0e3 *1.* alphah_cat(V*1.0e3))
CaTm2h0_m2h1 = smodel.VDepSReac('CaTm2h0_m2h1', ssys, slhs = [CaT_m2h0], srhs = [CaT_m2h1], k= lambda V: 1.0e3 *1.* alphah_cat(V*1.0e3))

CaTm2h1_m2h0 = smodel.VDepSReac('CaTm2h1_m2h0', ssys, slhs = [CaT_m2h1], srhs = [CaT_m2h0], k= lambda V: 1.0e3 *1.* betah_cat(V*1.0e3))
CaTm1h1_m1h0 = smodel.VDepSReac('CaTm1h1_m1h0', ssys, slhs = [CaT_m1h1], srhs = [CaT_m1h0], k= lambda V: 1.0e3 *1.* betah_cat(V*1.0e3))
CaTm0h1_m0h0 = smodel.VDepSReac('CaTm0h1_m0h0', ssys, slhs = [CaT_m0h1], srhs = [CaT_m0h0], k= lambda V: 1.0e3 *1.* betah_cat(V*1.0e3))

if cyl160:
    OC_CaT = smodel.GHKcurr('OC_CaT', ssys, CaT_m2h1, Ca, virtual_oconc = Ca_oconc, computeflux = True)
else:
    OC_CaT = smodel.GHKcurr('OC_CaT', ssys, CaT_m2h1, Ca, computeflux = True)

OC_CaT.setP(CaT_P)

##### BK channel ####################
BKchan = smodel.Chan('BKchan', mdl)

BK_C0 = smodel.ChanState('BK_C0', mdl, BKchan)
BK_C1 = smodel.ChanState('BK_C1', mdl, BKchan)
BK_C2 = smodel.ChanState('BK_C2', mdl, BKchan)
BK_C3 = smodel.ChanState('BK_C3', mdl, BKchan)
BK_C4 = smodel.ChanState('BK_C4', mdl, BKchan)
BK_O0 = smodel.ChanState('BK_O0', mdl, BKchan)
BK_O1 = smodel.ChanState('BK_O1', mdl, BKchan)
BK_O2 = smodel.ChanState('BK_O2', mdl, BKchan)
BK_O3 = smodel.ChanState('BK_O3', mdl, BKchan)
BK_O4 = smodel.ChanState('BK_O4', mdl, BKchan)

BKCAC0 = smodel.SReac('BKCAC0', ssys, slhs = [BK_C0], ilhs = [Ca], srhs = [BK_C1], kcst = c_01)
BKCAC1 = smodel.SReac('BKCAC1', ssys, slhs = [BK_C1], ilhs = [Ca], srhs = [BK_C2], kcst = c_12)
BKCAC2 = smodel.SReac('BKCAC2', ssys, slhs = [BK_C2], ilhs = [Ca], srhs = [BK_C3], kcst = c_23)
BKCAC3 = smodel.SReac('BKCAC3', ssys, slhs = [BK_C3], ilhs = [Ca], srhs = [BK_C4], kcst = c_34)

BKC0 = smodel.SReac('BKC0', ssys, slhs = [BK_C1], srhs = [BK_C0], irhs = [Ca], kcst = c_10)
BKC1 = smodel.SReac('BKC1', ssys, slhs = [BK_C2], srhs = [BK_C1], irhs = [Ca], kcst = c_21)
BKC2 = smodel.SReac('BKC2', ssys, slhs = [BK_C3], srhs = [BK_C2], irhs = [Ca], kcst = c_32)
BKC3 = smodel.SReac('BKC3', ssys, slhs = [BK_C4], srhs = [BK_C3], irhs = [Ca], kcst = c_43)

BKCAO0 = smodel.SReac('BKCAO0', ssys, slhs = [BK_O0], ilhs = [Ca], srhs = [BK_O1], kcst = o_01)
BKCAO1 = smodel.SReac('BKCAO1', ssys, slhs = [BK_O1], ilhs = [Ca], srhs = [BK_O2], kcst = o_12)
BKCAO2 = smodel.SReac('BKCAO2', ssys, slhs = [BK_O2], ilhs = [Ca], srhs = [BK_O3], kcst = o_23)
BKCAO3 = smodel.SReac('BKCAO3', ssys, slhs = [BK_O3], ilhs = [Ca], srhs = [BK_O4], kcst = o_34)

BKO0 = smodel.SReac('BKO0', ssys, slhs = [BK_O1], srhs = [BK_O0], irhs = [Ca], kcst = o_10)
BKO1 = smodel.SReac('BKO1', ssys, slhs = [BK_O2], srhs = [BK_O1], irhs = [Ca], kcst = o_21)
BKO2 = smodel.SReac('BKO2', ssys, slhs = [BK_O3], srhs = [BK_O2], irhs = [Ca], kcst = o_32)
BKO3 = smodel.SReac('BKO3', ssys, slhs = [BK_O4], srhs = [BK_O3], irhs = [Ca], kcst = o_43)

BKC0O0 = smodel.VDepSReac('BKC0O0', ssys, slhs = [BK_C0], srhs = [BK_O0], k=lambda V: f_0(V))
BKC1O1 = smodel.VDepSReac('BKC1O1', ssys, slhs = [BK_C1], srhs = [BK_O1], k=lambda V: f_1(V))
BKC2O2 = smodel.VDepSReac('BKC2O2', ssys, slhs = [BK_C2], srhs = [BK_O2], k=lambda V: f_2(V))
BKC3O3 = smodel.VDepSReac('BKC3O3', ssys, slhs = [BK_C3], srhs = [BK_O3], k=lambda V: f_3(V))
BKC4O4 = smodel.VDepSReac('BKC4O4', ssys, slhs = [BK_C4], srhs = [BK_O4], k=lambda V: f_4(V))

BKO0C0 = smodel.VDepSReac('BKO0C0', ssys, slhs = [BK_O0], srhs = [BK_C0], k=lambda V: b_0(V))
BKO1C1 = smodel.VDepSReac('BKO1C1', ssys, slhs = [BK_O1], srhs = [BK_C1], k=lambda V: b_1(V))
BKO2C2 = smodel.VDepSReac('BKO2C2', ssys, slhs = [BK_O2], srhs = [BK_C2], k=lambda V: b_2(V))
BKO3C3 = smodel.VDepSReac('BKO3C3', ssys, slhs = [BK_O3], srhs = [BK_C3], k=lambda V: b_3(V))
BKO4C4 = smodel.VDepSReac('BKO4C4', ssys, slhs = [BK_O4], srhs = [BK_C4], k=lambda V: b_4(V))

OC_BK0 = smodel.OhmicCurr('OC_BK0', ssys, chanstate = BK_O0, erev = BK_rev, g = BK_G )
OC_BK1 = smodel.OhmicCurr('OC_BK1', ssys, chanstate = BK_O1, erev = BK_rev, g = BK_G )
OC_BK2 = smodel.OhmicCurr('OC_BK2', ssys, chanstate = BK_O2, erev = BK_rev, g = BK_G )
OC_BK3 = smodel.OhmicCurr('OC_BK3', ssys, chanstate = BK_O3, erev = BK_rev, g = BK_G )
OC_BK4 = smodel.OhmicCurr('OC_BK4', ssys, chanstate = BK_O4, erev = BK_rev, g = BK_G )


###### SK channel ################## DETERMINISTIC
SKchan = smodel.Chan('SKchan', mdl)

SK_C1 = smodel.ChanState('SK_C1', mdl, SKchan)
SK_C2 = smodel.ChanState('SK_C2', mdl, SKchan)
SK_C3 = smodel.ChanState('SK_C3', mdl, SKchan)
SK_C4 = smodel.ChanState('SK_C4', mdl, SKchan)
SK_O1 = smodel.ChanState('SK_O1', mdl, SKchan)
SK_O2 = smodel.ChanState('SK_O2', mdl, SKchan)


SKCAC1 = smodel.SReac('SKCAC1', ssys, slhs = [SK_C1], ilhs = [Ca], srhs = [SK_C2], kcst = dirc2_t)
SKCAC2 = smodel.SReac('SKCAC2', ssys, slhs = [SK_C2], ilhs = [Ca], srhs = [SK_C3], kcst = dirc3_t)
SKCAC3 = smodel.SReac('SKCAC3', ssys, slhs = [SK_C3], ilhs = [Ca], srhs = [SK_C4], kcst = dirc4_t)

SKC1 = smodel.SReac('SKC1', ssys, slhs = [SK_C2], srhs = [SK_C1], irhs = [Ca], kcst = invc1_t)
SKC2 = smodel.SReac('SKC2', ssys, slhs = [SK_C3], srhs = [SK_C2], irhs = [Ca], kcst = invc2_t)
SKC3 = smodel.SReac('SKC3', ssys, slhs = [SK_C4], srhs = [SK_C3], irhs = [Ca], kcst = invc3_t)

SKC3O1 = smodel.SReac('SKC3O1', ssys, slhs = [SK_C3], srhs = [SK_O1], kcst = diro1_t)
SKC4O2 = smodel.SReac('SKC4O2', ssys, slhs = [SK_C4], srhs = [SK_O2], kcst = diro2_t)

SKO1C3 = smodel.SReac('SKO1C3', ssys, slhs = [SK_O1], srhs = [SK_C3], kcst = invo1_t)
SKO2C4 = smodel.SReac('SKO2C4', ssys, slhs = [SK_O2], srhs = [SK_C4], kcst = invo2_t)

OC1_SK = smodel.OhmicCurr('OC1_SK', ssys, chanstate = SK_O1, erev = SK_rev, g = SK_G )
OC2_SK = smodel.OhmicCurr('OC2_SK', ssys, chanstate = SK_O2, erev = SK_rev, g = SK_G )

###### Leak current channel #####

L = smodel.Chan('L', mdl)
Leak = smodel.ChanState('Leak', mdl, L)

OC_L = smodel.OhmicCurr('OC_L', ssys, chanstate = Leak, erev = L_rev, g = L_G)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

########### MESH & COMPARTMENTALIZATION #################

##########Import Mesh

mesh = meshio.loadMesh('./meshes/'+meshfile_ab)[0]

outer_tets = range(mesh.ntets)
inner_tets = gettets.getcyl(mesh, 1e-6, -200e-6, 200e-6)[0]

for i in inner_tets: outer_tets.remove(i)
assert(outer_tets.__len__() + inner_tets.__len__() == mesh.ntets)

print outer_tets.__len__(), " tets in outer compartment"
print inner_tets.__len__(), " tets in inner compartment"

# Record voltage from the central tetrahedron
cent_tet = mesh.findTetByPoint([0.0,0.0,0.0])

########## Create an intracellular compartment i.e. cytosolic compartment

cyto = sgeom.TmComp('cyto', mesh, inner_tets)
cyto.addVolsys('vsys')

if not cyl160: outer = sgeom.TmComp('outer', mesh, outer_tets)

if cyl160:
    # Ensure that we use points a small distance inside the boundary:
    LENGTH = mesh.getBoundMax()[2] - mesh.getBoundMin()[2]
    boundminz = mesh.getBoundMin()[2] + LENGTH/mesh.ntets
    boundmaxz = mesh.getBoundMax()[2] - LENGTH/mesh.ntets

    memb_tris = list(mesh.getSurfTris())
    minztris = []
    maxztris = []
    for tri in memb_tris:
        zminboundtri = True
        zmaxboundtri = True
        tritemp = mesh.getTri(tri)
        trizs = [0.0, 0.0, 0.0]
        trizs[0] = mesh.getVertex(tritemp[0])[2]
        trizs[1] = mesh.getVertex(tritemp[1])[2]
        trizs[2] = mesh.getVertex(tritemp[2])[2]
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
            tritemp = mesh.getTetTriNeighb(i)
            for j in range(4): out_tris.add(tritemp[j])

    in_tris = set()
    for i in inner_tets:
            tritemp = mesh.getTetTriNeighb(i)
            for j in range(4): in_tris.add(tritemp[j])

    memb_tris = out_tris.intersection(in_tris)
    memb_tris = list(memb_tris)

memb_tris_bk = []
memb_tris_sk = []
memb_tris_cat = []
memb_tris_cap = []

surfarea = 0.0

for i in memb_tris:
    surfarea = surfarea + mesh.getTriArea(i)

random.seed(7)

while (len(memb_tris_bk)<round(BK_ro*surfarea)):
    ctriID = random.choice(memb_tris)
    if ctriID not in memb_tris_bk:
        memb_tris_bk.append(ctriID)

while (len(memb_tris_sk)<round(SK_ro*surfarea)):
    ctriID = random.choice(memb_tris)
    if ctriID not in memb_tris_sk:
        memb_tris_sk.append(ctriID)

while (len(memb_tris_cat)<round(CaT_ro*surfarea)):
    ctriID = random.choice(memb_tris)
    if ctriID not in memb_tris_cat:
        memb_tris_cat.append(ctriID)
                    
while (len(memb_tris_cap)<(round(CaP_ro*surfarea)-(clusterSize*round(BK_ro*surfarea)))):
    ctriID = random.choice(memb_tris)
    if ctriID not in memb_tris_bk:
        memb_tris_cap.append(ctriID)




########## Find the submembrane tets

memb_tet_neighb = []
for i in memb_tris:
    tettemp = mesh.getTriTetNeighb(i)
    for j in tettemp:
        memb_tet_neighb.append(j)

submemb_tets = []
for i in memb_tet_neighb:
    if i in inner_tets:
        submemb_tets.append(i)

print len(submemb_tets)

vol = 0.0

for i in submemb_tets:
    vol = vol + mesh.getTetVol(i)

print 'Volume of submembrane region is', vol

submemb_tets_surftris = dict()

for m in submemb_tets:
    tris = mesh.getTetTriNeighb(m)
    for t in tris:
        if t in memb_tris:
            submemb_tets_surftris[m] = t
            break

assert(len(submemb_tets_surftris.values()) == len(submemb_tets))

for i in range(len(memb_tris)):
    ctri = memb_tris[i]
    ctet = submemb_tets[i]
    tettemp = mesh.getTriTetNeighb(ctri)
    if not ctet in tettemp:
        print 'Tri and Tet do not correspond to each other'


border_tets = []
border_tets_vols = 0.0

for i in inner_tets:
    tritemp = mesh.getTetTriNeighb(i)
    for t in tritemp:
        if t in memb_tris:
            border_tets.append(i)
            border_tets_vols+=mesh.getTetVol(i)
            break

print "Border tet vols:", border_tets_vols


########## Create a membrane as a surface mesh
if cyl160: 
    memb = sgeom.TmPatch('memb', mesh, memb_tris, cyto)
else:
    memb = sgeom.TmPatch('memb', mesh, memb_tris, cyto, outer)
memb.addSurfsys('ssys')

# For EField calculation
print "Creating membrane.."
membrane = sgeom.Memb('membrane', mesh, [memb])
print "Membrane created."

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

r = srng.create_mt19937(512)
r.initialize(7)

print "Creating tetexact solver..."
sim = ssolver.Tetexact(mdl, mesh, r, True)

print "Resetting simulation objects.."
sim.reset()

print "Injecting molecules.."

sim.setTemp(TEMPERATURE+273.15)

if not cyl160: 
    sim.setCompConc('outer', 'Ca', Ca_oconc)
    sim.setCompClamped('outer', 'Ca', True)
    
sim.setCompConc('cyto', 'Ca', Ca_iconc)

print "Calcium concentration is: ", sim.getCompConc('cyto', 'Ca')
print "No. of Ca molecules is: ", sim.getCompCount('cyto', 'Ca')

sim.setCompConc('cyto', 'Mg', Mg_conc)


surfarea = sim.getPatchArea('memb')

pumpnbs = 6.022141e12*surfarea

sim.setPatchCount('memb', 'Pump', pumpnbs)
sim.setPatchCount('memb', 'CaPump', 0)

print "Injected ", sim.getPatchCount('memb', 'Pump'), "pumps"

sim.setCompConc('cyto', 'iCBsf',    iCBsf_conc)
sim.setCompConc('cyto', 'iCBsCa',   iCBsCa_conc)
sim.setCompConc('cyto', 'iCBCaf',   iCBCaf_conc)
sim.setCompConc('cyto', 'iCBCaCa',  iCBCaCa_conc)

sim.setCompConc('cyto', 'CBsf',     CBsf_conc)
sim.setCompConc('cyto', 'CBsCa',    CBsCa_conc)
sim.setCompConc('cyto', 'CBCaf',    CBCaf_conc)
sim.setCompConc('cyto', 'CBCaCa',   CBCaCa_conc)

sim.setCompConc('cyto', 'PV',       PV_conc)
sim.setCompConc('cyto', 'PVCa',     PVCa_conc)
sim.setCompConc('cyto', 'PVMg',     PVMg_conc)

bk_c0_count = 0
bk_c1_count = 0
bk_c2_count = 0
bk_c3_count = 0
bk_c4_count = 0

bk_o0_count = 0
bk_o1_count = 0
bk_o2_count = 0
bk_o3_count = 0
bk_o4_count = 0

for i in memb_tris_bk:
    if bk_c0_count<round(BK_ro*surfarea*BK_C0_p):
        sim.setTriCount(i, 'BK_C0', sim.getTriCount(i, 'BK_C0') + 1)
        bk_c0_count = bk_c0_count + 1
    elif bk_c1_count<round(BK_ro*surfarea*BK_C1_p):
        sim.setTriCount(i, 'BK_C1', sim.getTriCount(i, 'BK_C1') + 1)
        bk_c1_count = bk_c1_count + 1
    elif bk_c2_count<round(BK_ro*surfarea*BK_C2_p):
        sim.setTriCount(i, 'BK_C2', sim.getTriCount(i, 'BK_C2') + 1)
        bk_c2_count = bk_c2_count + 1
    elif bk_c3_count<round(BK_ro*surfarea*BK_C3_p):
        sim.setTriCount(i, 'BK_C3', sim.getTriCount(i, 'BK_C3') + 1)
        bk_c3_count = bk_c3_count + 1
    elif bk_c4_count<round(BK_ro*surfarea*BK_C4_p):
        sim.setTriCount(i, 'BK_C4', sim.getTriCount(i, 'BK_C4') + 1)
        bk_c4_count = bk_c4_count + 1
    elif bk_o0_count<round(BK_ro*surfarea*BK_O0_p):
        sim.setTriCount(i, 'BK_O0', sim.getTriCount(i, 'BK_O0') + 1)
        bk_o0_count = bk_o0_count + 1
    elif bk_o1_count<round(BK_ro*surfarea*BK_O1_p):
        sim.setTriCount(i, 'BK_O1', sim.getTriCount(i, 'BK_O1') + 1)
        bk_o1_count = bk_o1_count + 1
    elif bk_o2_count<round(BK_ro*surfarea*BK_O2_p):
        sim.setTriCount(i, 'BK_O2', sim.getTriCount(i, 'BK_O2') + 1)
        bk_o2_count = bk_o2_count + 1
    elif bk_o3_count<round(BK_ro*surfarea*BK_O3_p):
        sim.setTriCount(i, 'BK_O3', sim.getTriCount(i, 'BK_O3') + 1)
        bk_o3_count = bk_o3_count + 1
    elif bk_o4_count<round(BK_ro*surfarea*BK_O4_p):
        sim.setTriCount(i, 'BK_O4', sim.getTriCount(i, 'BK_O4') + 1)
        bk_o4_count = bk_o4_count + 1
    else:
        print 'More tris picked up by algorithm than the number of BK channels'

sk_c1_count = 0
sk_c2_count = 0
sk_c3_count = 0
sk_c4_count = 0

sk_o1_count = 0
sk_o2_count = 0

for i in memb_tris_sk:
    if sk_c1_count<round(SK_ro*surfarea*SK_C1_p):
        sim.setTriCount(i, 'SK_C1', sim.getTriCount(i, 'SK_C1') + 1)
        sk_c1_count = sk_c1_count + 1
    elif sk_c2_count<round(SK_ro*surfarea*SK_C2_p):
        sim.setTriCount(i, 'SK_C2', sim.getTriCount(i, 'SK_C2') + 1)
        sk_c2_count = sk_c2_count + 1
    elif sk_c3_count<round(SK_ro*surfarea*SK_C3_p):
        sim.setTriCount(i, 'SK_C3', sim.getTriCount(i, 'SK_C3') + 1)
        sk_c3_count = sk_c3_count + 1
    elif sk_c4_count<round(SK_ro*surfarea*SK_C4_p):
        sim.setTriCount(i, 'SK_C4', sim.getTriCount(i, 'SK_C4') + 1)
        sk_c4_count = sk_c4_count + 1
    elif sk_o1_count<round(SK_ro*surfarea*SK_O1_p):
        sim.setTriCount(i, 'SK_O1', sim.getTriCount(i, 'SK_O1') + 1)
        sk_o1_count = sk_o1_count + 1
    elif sk_o2_count<round(SK_ro*surfarea*SK_O2_p):
        sim.setTriCount(i, 'SK_O2', sim.getTriCount(i, 'SK_O2') + 1)
        sk_o2_count = sk_o2_count + 1
    else:
        print 'More tris picked up by algorithm than the number of SK channels'


cat_m0h0_count = 0
cat_m1h0_count = 0
cat_m2h0_count = 0
cat_m0h1_count = 0
cat_m1h1_count = 0
cat_m2h1_count = 0

for i in memb_tris_cat:
    if cat_m0h0_count<round(CaT_ro*surfarea*CaT_m0h0_p):
        sim.setTriCount(i, 'CaT_m0h0', sim.getTriCount(i, 'CaT_m0h0') + 1)
        cat_m0h0_count = cat_m0h0_count + 1
    elif cat_m1h0_count<round(CaT_ro*surfarea*CaT_m1h0_p):
        sim.setTriCount(i, 'CaT_m1h0', sim.getTriCount(i, 'CaT_m1h0') + 1)
        cat_m1h0_count = cat_m1h0_count + 1
    elif cat_m2h0_count<round(CaT_ro*surfarea*CaT_m2h0_p):
        sim.setTriCount(i, 'CaT_m2h0', sim.getTriCount(i, 'CaT_m2h0') + 1)
        cat_m2h0_count = cat_m2h0_count + 1
    elif cat_m0h1_count<round(CaT_ro*surfarea*CaT_m0h1_p):
        sim.setTriCount(i, 'CaT_m0h1', sim.getTriCount(i, 'CaT_m0h1') + 1)
        cat_m0h1_count = cat_m0h1_count + 1
    elif cat_m1h1_count<round(CaT_ro*surfarea*CaT_m1h1_p):
        sim.setTriCount(i, 'CaT_m1h1', sim.getTriCount(i, 'CaT_m1h1') + 1)
        cat_m1h1_count = cat_m1h1_count + 1
    elif cat_m2h1_count<round(CaT_ro*surfarea*CaT_m2h1_p):
        sim.setTriCount(i, 'CaT_m2h1', sim.getTriCount(i, 'CaT_m2h1') + 1)
        cat_m2h1_count = cat_m2h1_count + 1
    else:
        print 'More tris picked up by algorithm than the number of CaT channels'

cap_m0_count = 0
cap_m1_count = 0
cap_m2_count = 0
cap_m3_count = 0

if clusterSize>0:
    for i in memb_tris_bk:
        count = 0
        while count<clusterSize:
            if cap_m0_count<round(CaP_ro*surfarea*CaP_m0_p):
                sim.setTriCount(i, 'CaP_m0', sim.getTriCount(i, 'CaP_m0') + 1)
                cap_m0_count = cap_m0_count + 1
		count = count +1
            elif cap_m1_count<round(CaP_ro*surfarea*CaP_m1_p):
                sim.setTriCount(i, 'CaP_m1', sim.getTriCount(i, 'CaP_m1') + 1)
                cap_m1_count = cap_m1_count + 1
		count = count +1
            elif cap_m2_count<round(CaP_ro*surfarea*CaP_m2_p):
                sim.setTriCount(i, 'CaP_m2', sim.getTriCount(i, 'CaP_m2') + 1)
                cap_m2_count = cap_m2_count + 1
		count = count +1
            elif cap_m3_count<round(CaP_ro*surfarea*CaP_m3_p):
                sim.setTriCount(i, 'CaP_m3', sim.getTriCount(i, 'CaP_m3') + 1)
                cap_m3_count = cap_m3_count + 1
		count = count +1
            else:
                print 'Cluster size is larger than the number of CaP channels available'

for i in memb_tris_cap:
    if cap_m0_count<round(CaP_ro*surfarea*CaP_m0_p):
        sim.setTriCount(i, 'CaP_m0', sim.getTriCount(i, 'CaP_m0') + 1)
        cap_m0_count = cap_m0_count + 1
    elif cap_m1_count<round(CaP_ro*surfarea*CaP_m1_p):
        sim.setTriCount(i, 'CaP_m1', sim.getTriCount(i, 'CaP_m1') + 1)
        cap_m1_count = cap_m1_count + 1
    elif cap_m2_count<round(CaP_ro*surfarea*CaP_m2_p):
        sim.setTriCount(i, 'CaP_m2', sim.getTriCount(i, 'CaP_m2') + 1)
        cap_m2_count = cap_m2_count + 1
    elif cap_m3_count<round(CaP_ro*surfarea*CaP_m3_p):
        sim.setTriCount(i, 'CaP_m3', sim.getTriCount(i, 'CaP_m3') + 1)
        cap_m3_count = cap_m3_count + 1
    else:
        print 'More tris picked up by the algorithm than the number of CaP channels available'
                                                                                                                                                                            

sim.setPatchCount('memb', 'Leak', int(L_ro * surfarea))
print "Injected  ", int(L_ro * sim.getPatchArea('memb')), "Leak channels"

memb_triID_withBK=[]
memb_countBK_pertriID=[]

memb_tetID_withBK=[]

count = 0
for m in memb_tris:
    BKchans=sim.getTriCount(m,'BK_C0')+sim.getTriCount(m,'BK_C1')+sim.getTriCount(m,'BK_C2')+sim.getTriCount(m,'BK_C3')+sim.getTriCount(m,'BK_C4')+sim.getTriCount(m,'BK_O0')+sim.getTriCount(m,'BK_O1')+sim.getTriCount(m,'BK_O2')+sim.getTriCount(m,'BK_O3')+sim.getTriCount(m,'BK_O4')
    if (BKchans>0):
        memb_triID_withBK.append(m)
        memb_countBK_pertriID.append(BKchans)
        memb_tetID_withBK.append(submemb_tets[count])
    count = count+1
                                    


sim.setEfieldDT(EF_DT)
sim.setMembPotential('membrane', init_pot)
sim.setMembVolRes('membrane', Ra)
#cm = 1.5uF/cm2 -> 1.5e-6F/1e-4m2 ->1.5e-2 F/m2
sim.setMembCapac('membrane',memb_capac)


#### Recording #####

c=time.ctime()

dc = c.split()[1]+c.split()[2]+'_'+c.split()[3]+'_'+c.split()[4]
dc= dc.replace(':', '_')

try: os.mkdir(root+'data')
except: pass
try: os.mkdir(root+'data/' +  'StochasticCaburst_cluster')
except: pass
try: os.mkdir(root+'data/' +  'StochasticCaburst_cluster/'+meshfile_ab)
except: pass 

os.mkdir(root+'data/' +  'StochasticCaburst_cluster/'+meshfile_ab+'/'+iter_n+'__'+dc )


datfile =  open(root+'data/' +  'StochasticCaburst_cluster/'+meshfile_ab+'/'+iter_n+'__'+dc + '/currents.dat', 'w')
datfile2 = open(root+'data/' +  'StochasticCaburst_cluster/'+meshfile_ab+'/'+iter_n+'__'+dc + '/voltage.dat', 'w')
datfile3 = open(root+'data/' +  'StochasticCaburst_cluster/'+meshfile_ab+'/'+iter_n+'__'+dc + '/calcium.dat', 'w')
datfile4 = open(root+'data/' +  'StochasticCaburst_cluster/'+meshfile_ab+'/'+iter_n+'__'+dc + '/OpenBKandCa.dat', 'w')
datfile5 = open(root+'data/' +  'StochasticCaburst_cluster/'+meshfile_ab+'/'+iter_n+'__'+dc + '/ChannelsDistribution.dat', 'w')

for m in memb_tris:
    tri_center = mesh.getTriBarycenter(m)
    cap_m0_chans = sim.getTriCount(m,'CaP_m0')
    cap_m1_chans = sim.getTriCount(m,'CaP_m1')
    cap_m2_chans = sim.getTriCount(m,'CaP_m2')
    cap_m3_chans = sim.getTriCount(m,'CaP_m3')
    cat_m0h0_chans = sim.getTriCount(m,'CaT_m0h0')
    cat_m1h0_chans = sim.getTriCount(m,'CaT_m1h0')
    cat_m2h0_chans = sim.getTriCount(m,'CaT_m2h0')
    cat_m0h1_chans = sim.getTriCount(m,'CaT_m0h1')
    cat_m1h1_chans = sim.getTriCount(m,'CaT_m1h1')
    cat_m2h1_chans = sim.getTriCount(m,'CaT_m2h1')
    bk_c0_chans = sim.getTriCount(m,'BK_C0')
    bk_c1_chans = sim.getTriCount(m,'BK_C1')
    bk_c2_chans = sim.getTriCount(m,'BK_C2')
    bk_c3_chans = sim.getTriCount(m,'BK_C3')
    bk_c4_chans = sim.getTriCount(m,'BK_C4')
    bk_o0_chans = sim.getTriCount(m,'BK_O0')
    bk_o1_chans = sim.getTriCount(m,'BK_O1')
    bk_o2_chans = sim.getTriCount(m,'BK_O2')
    bk_o3_chans = sim.getTriCount(m,'BK_O3')
    bk_o4_chans = sim.getTriCount(m,'BK_O4')
    sk_c1_chans = sim.getTriCount(m,'SK_C1')
    sk_c2_chans = sim.getTriCount(m,'SK_C2')
    sk_c3_chans = sim.getTriCount(m,'SK_C3')
    sk_c4_chans = sim.getTriCount(m,'SK_C4')
    sk_o1_chans = sim.getTriCount(m,'SK_O1')
    sk_o2_chans = sim.getTriCount(m,'SK_O2')
    datfile5.write('%.6g' %(tri_center[0]) + ' ')
    datfile5.write('%.6g' %(tri_center[1]) + ' ')
    datfile5.write('%.6g' %(tri_center[2]) + ' ')
    datfile5.write('%.6g' %(cap_m0_chans) + ' ')
    datfile5.write('%.6g' %(cap_m1_chans) + ' ')
    datfile5.write('%.6g' %(cap_m2_chans) + ' ')
    datfile5.write('%.6g' %(cap_m3_chans) + ' ')
    datfile5.write('%.6g' %(cat_m0h0_chans) + ' ')
    datfile5.write('%.6g' %(cat_m1h0_chans) + ' ')
    datfile5.write('%.6g' %(cat_m2h0_chans) + ' ')
    datfile5.write('%.6g' %(cat_m0h1_chans) + ' ')
    datfile5.write('%.6g' %(cat_m1h1_chans) + ' ')
    datfile5.write('%.6g' %(cat_m2h1_chans) + ' ')
    datfile5.write('%.6g' %(bk_c0_chans) + ' ')
    datfile5.write('%.6g' %(bk_c1_chans) + ' ')
    datfile5.write('%.6g' %(bk_c2_chans) + ' ')
    datfile5.write('%.6g' %(bk_c3_chans) + ' ')
    datfile5.write('%.6g' %(bk_c4_chans) + ' ')
    datfile5.write('%.6g' %(bk_o0_chans) + ' ')
    datfile5.write('%.6g' %(bk_o1_chans) + ' ')
    datfile5.write('%.6g' %(bk_o2_chans) + ' ')
    datfile5.write('%.6g' %(bk_o3_chans) + ' ')
    datfile5.write('%.6g' %(bk_o4_chans) + ' ')
    datfile5.write('%.6g' %(sk_c1_chans) + ' ')
    datfile5.write('%.6g' %(sk_c2_chans) + ' ')
    datfile5.write('%.6g' %(sk_c3_chans) + ' ')
    datfile5.write('%.6g' %(sk_c3_chans) + ' ')
    datfile5.write('%.6g' %(sk_o1_chans) + ' ')
    datfile5.write('%.6g' %(sk_o2_chans) + ' ')
    datfile5.write('\n')
															

r.initialize(int(time.time()%1000))

for l in range(NTIMEPOINTS):
    print "Tpnt: ", l

    sim.run(TIMECONVERTER*l)
    
    tcur_CaP = 0.0
    tcur_CaT = 0.0
    tcur_BK = 0.0
    tcur_SK = 0.0
    tca_count = 0.0

    So = Ca_oconc
    
    for m in submemb_tets:
        ctriID = submemb_tets_surftris[m]
        tcur_CaP = tcur_CaP + sim.getTriGHKI(ctriID,'OC_CaP')
        tcur_CaT = tcur_CaT + sim.getTriGHKI(ctriID,'OC_CaT')
        tcur_BK = tcur_BK + sim.getTriOhmicI(ctriID,'OC_BK0') \
	    + sim.getTriOhmicI(ctriID,'OC_BK1') \
	    + sim.getTriOhmicI(ctriID,'OC_BK2') \
	    + sim.getTriOhmicI(ctriID,'OC_BK3') \
	    + sim.getTriOhmicI(ctriID,'OC_BK4')
        tcur_SK = tcur_SK + sim.getTriOhmicI(ctriID,'OC1_SK') + sim.getTriOhmicI(ctriID,'OC2_SK')
        tca_count = tca_count + sim.getTetCount(m,'Ca')
    
    datfile.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile.write('%.6g' %((tcur_CaP*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_CaT*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_BK*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_SK*1.0e-1)/surfarea) + ' ')  
    datfile.write('\n')

    datfile2.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile2.write('%.6g' %(sim.getTetV(cent_tet)*1.0e3) + ' ')
    datfile2.write('\n')
    
    datfile3.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile3.write('%.6g' %(((tca_count/AVOGADRO)/(border_tets_vols*1.0e3))*1.0e6) +' ')
    datfile3.write('%.6g' %(tca_count)+ ' ')
    datfile3.write('\n')

    datfile4.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    for i in range(len(memb_triID_withBK)):
        datfile4.write('%.6g' %sim.getTriCount(memb_triID_withBK[i], 'BK_O0') + ' ')
        datfile4.write('%.6g' %sim.getTriCount(memb_triID_withBK[i], 'BK_O1') + ' ')
        datfile4.write('%.6g' %sim.getTriCount(memb_triID_withBK[i], 'BK_O2') + ' ')
        datfile4.write('%.6g' %sim.getTriCount(memb_triID_withBK[i], 'BK_O3') + ' ')
        datfile4.write('%.6g' %sim.getTriCount(memb_triID_withBK[i], 'BK_O4') + ' ')
        datfile4.write('%.6g' %sim.getTetCount(memb_tetID_withBK[i], 'Ca') + ' ')
        datfile4.write('%.6g' %sim.getTetConc(memb_tetID_withBK[i], 'Ca') +' ')
    datfile4.write('\n')
            

datfile.close()
datfile2.close()
datfile3.close()
datfile4.close()
datfile5.close()


