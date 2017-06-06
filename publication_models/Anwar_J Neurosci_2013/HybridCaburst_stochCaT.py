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
# *HybridCaburst_stochCaT.py : A hybrid calcium burst model with stochastic T-type 
# calcium channels and everything else deterministic.
#
# Script authors: Haroon Anwar and Iain Hepburn
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# USAGE
#
# $ python HybridCaburst_stochCaT.py *mesh* *root* *iter_n* 
#  
#  *mesh* is the tetrahedral mesh (10um to 160um cylinder)
#  *root* is the path to the location for data storage 
#  *iter_n* (is intened to be an integer) is an identifier number for each 
#     simulation iteration.
#
# E.g: python HybridCaburst_stochCaT.py Cylinder2_dia2um_L10um_outer0_3um_0.3shell_0.3size_19156tets_adaptive.inp ~/stochcasims/ 1
#
#
# OUTPUT 
#
# In (root)/data/HybridCaburst_stochCaT/(mesh)/(iter_n+time) directory 
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
# Time (ms), determinstic calcium concentration in submembrane (micromolar),
# stochastic calcium concentration in submembrane (micromolar),
# number of calcium ions in submembrane in deterministic solver, 
# number of calcium ions in submembrane in stochastic solver. 
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
mdl_det = smodel.Model()

# Calcium
Ca_det = smodel.Spec('Ca_det', mdl_det)
Ca_det.setValence(2)

Ca_stoch = smodel.Spec('Ca_stoch', mdl_stoch)
Ca_stoch.setValence(2)

# Pump
Pump = smodel.Spec('Pump', mdl_det)
# CaPump
CaPump = smodel.Spec('CaPump', mdl_det)

# iCBsf
iCBsf = smodel.Spec('iCBsf', mdl_det)
# iCBsCa
iCBsCa = smodel.Spec('iCBsCa', mdl_det)
# iCBCaf
iCBCaf = smodel.Spec('iCBCaf', mdl_det)
# iCBCaCa
iCBCaCa = smodel.Spec('iCBCaCa', mdl_det)

# CBsf
CBsf = smodel.Spec('CBsf', mdl_det)
# CBsCa
CBsCa = smodel.Spec('CBsCa', mdl_det)
# CBCaf
CBCaf = smodel.Spec('CBCaf', mdl_det)
# CBCaCa
CBCaCa = smodel.Spec('CBCaCa', mdl_det)

# PV
PV = smodel.Spec('PV', mdl_det)
# PVMg
PVMg = smodel.Spec('PVMg', mdl_det)
# PVCa
PVCa = smodel.Spec('PVCa', mdl_det)
# Mg
Mg = smodel.Spec('Mg', mdl_det)

# Vol/surface systems

vsys_stoch = smodel.Volsys('vsys_stoch', mdl_stoch)
ssys_stoch = smodel.Surfsys('ssys_stoch', mdl_stoch)

vsys_det = smodel.Volsys('vsys_det', mdl_det)
ssys_det = smodel.Surfsys('ssys_det', mdl_det)


# Diffusions 
diff_Ca = smodel.Diff('diff_Ca', vsys_det, Ca_det)
diff_Ca.setDcst(DCST)
diff_CBsf = smodel.Diff('diff_CBsf', vsys_det, CBsf)
diff_CBsf.setDcst(DCB)
diff_CBsCa = smodel.Diff('diff_CBsCa', vsys_det, CBsCa)
diff_CBsCa.setDcst(DCB)
diff_CBCaf = smodel.Diff('diff_CBCaf', vsys_det, CBCaf)
diff_CBCaf.setDcst(DCB)
diff_CBCaCa = smodel.Diff('diff_CBCaCa', vsys_det, CBCaCa)
diff_CBCaCa.setDcst(DCB)
diff_PV = smodel.Diff('diff_PV', vsys_det, PV)
diff_PV.setDcst(DPV)
diff_PVCa = smodel.Diff('diff_PVCa', vsys_det, PVCa)
diff_PVCa.setDcst(DPV)
diff_PVMg = smodel.Diff('diff_PVMg', vsys_det, PVMg)
diff_PVMg.setDcst(DPV)


#Pump
PumpD_f = smodel.SReac('PumpD_f', ssys_det, ilhs=[Ca_det], slhs=[Pump], srhs=[CaPump])
PumpD_f.setKcst(P_f_kcst)

PumpD_b = smodel.SReac('PumpD_b', ssys_det, slhs=[CaPump], irhs=[Ca_det], srhs=[Pump])
PumpD_b.setKcst(P_b_kcst)

PumpD_k = smodel.SReac('PumpD_k', ssys_det, slhs=[CaPump], srhs=[Pump])
PumpD_k.setKcst(P_k_kcst)

#iCBsf-fast
iCBsf1_f = smodel.Reac('iCBsf1_f', vsys_det, lhs=[Ca_det,iCBsf], rhs=[iCBsCa], kcst = iCBsf1_f_kcst)
iCBsf1_b = smodel.Reac('iCBsf1_b', vsys_det, lhs=[iCBsCa], rhs=[Ca_det, iCBsf], kcst = iCBsf1_b_kcst)

#iCBsCa
iCBsCa_f = smodel.Reac('iCBsCa_f', vsys_det, lhs=[Ca_det,iCBsCa], rhs=[iCBCaCa], kcst = iCBsCa_f_kcst)
iCBsCa_b = smodel.Reac('iCBsCa_b', vsys_det, lhs=[iCBCaCa], rhs=[Ca_det,iCBsCa], kcst = iCBsCa_b_kcst)

#iCBsf_slow
iCBsf2_f = smodel.Reac('iCBsf2_f', vsys_det, lhs=[Ca_det,iCBsf], rhs=[iCBCaf], kcst = iCBsf2_f_kcst)
iCBsf2_b = smodel.Reac('iCBsf2_b', vsys_det, lhs=[iCBCaf], rhs=[Ca_det,iCBsf], kcst = iCBsf2_b_kcst)

#iCBCaf
iCBCaf_f = smodel.Reac('iCBCaf_f', vsys_det, lhs=[Ca_det,iCBCaf], rhs=[iCBCaCa], kcst = iCBCaf_f_kcst)
iCBCaf_b = smodel.Reac('iCBCaf_b', vsys_det, lhs=[iCBCaCa], rhs=[Ca_det,iCBCaf], kcst = iCBCaf_b_kcst)

#CBsf-fast
CBsf1_f = smodel.Reac('CBsf1_f', vsys_det, lhs=[Ca_det,CBsf], rhs=[CBsCa], kcst = CBsf1_f_kcst)
CBsf1_b = smodel.Reac('CBsf1_b', vsys_det, lhs=[CBsCa], rhs=[Ca_det,CBsf], kcst = CBsf1_b_kcst)

#CBsCa
CBsCa_f = smodel.Reac('CBsCa_f', vsys_det, lhs=[Ca_det,CBsCa], rhs=[CBCaCa], kcst = CBsCa_f_kcst)
CBsCa_b = smodel.Reac('CBsCa_b', vsys_det, lhs=[CBCaCa], rhs=[Ca_det,CBsCa], kcst = CBsCa_b_kcst)

#CBsf_slow
CBsf2_f = smodel.Reac('CBsf2_f', vsys_det, lhs=[Ca_det,CBsf], rhs=[CBCaf], kcst = CBsf2_f_kcst)
CBsf2_b = smodel.Reac('CBsf2_b', vsys_det, lhs=[CBCaf], rhs=[Ca_det,CBsf], kcst = CBsf2_b_kcst)

#CBCaf
CBCaf_f = smodel.Reac('CBCaf_f', vsys_det, lhs=[Ca_det,CBCaf], rhs=[CBCaCa], kcst = CBCaf_f_kcst)
CBCaf_b = smodel.Reac('CBCaf_b', vsys_det, lhs=[CBCaCa], rhs=[Ca_det,CBCaf], kcst = CBCaf_b_kcst)

#PVca
PVca_f = smodel.Reac('PVca_f', vsys_det, lhs=[Ca_det,PV], rhs=[PVCa], kcst = PVca_f_kcst)
PVca_b = smodel.Reac('PVca_b', vsys_det, lhs=[PVCa], rhs=[Ca_det,PV], kcst = PVca_b_kcst)

#PVmg
PVmg_f = smodel.Reac('PVmg_f', vsys_det, lhs=[Mg,PV], rhs=[PVMg], kcst = PVmg_f_kcst)
PVmg_b = smodel.Reac('PVmg_b', vsys_det, lhs=[PVMg], rhs=[Mg,PV], kcst = PVmg_b_kcst)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # CHANNELS  # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###### CaP channel ############## 

CaP_m0 = smodel.Spec('CaP_m0', mdl_det)
CaP_m1 = smodel.Spec('CaP_m1', mdl_det)
CaP_m2 = smodel.Spec('CaP_m2', mdl_det)
CaP_m3 = smodel.Spec('CaP_m3', mdl_det)


CaPm0m1 = smodel.SReac('CaPm0m1', ssys_det, slhs = [CaP_m0], srhs = [CaP_m1], kcst= 0.0)
CaPm1m2 = smodel.SReac('CaPm1m2', ssys_det, slhs = [CaP_m1], srhs = [CaP_m2], kcst= 0.0)
CaPm2m3 = smodel.SReac('CaPm2m3', ssys_det, slhs = [CaP_m2], srhs = [CaP_m3], kcst= 0.0)

CaPm3m2 = smodel.SReac('CaPm3m2', ssys_det, slhs = [CaP_m3], srhs = [CaP_m2], kcst= 0.0)
CaPm2m1 = smodel.SReac('CaPm2m1', ssys_det, slhs = [CaP_m2], srhs = [CaP_m1], kcst= 0.0)
CaPm1m0 = smodel.SReac('CaPm1m0', ssys_det, slhs = [CaP_m1], srhs = [CaP_m0], kcst= 0.0)


######## CaT channel ##########  

CaTchan = smodel.Chan('CaTchan', mdl_stoch)

CaT_m0h0 = smodel.ChanState('CaT_m0h0', mdl_stoch, CaTchan)
CaT_m0h1 = smodel.ChanState('CaT_m0h1', mdl_stoch, CaTchan)
CaT_m1h0 = smodel.ChanState('CaT_m1h0', mdl_stoch, CaTchan)
CaT_m1h1 = smodel.ChanState('CaT_m1h1', mdl_stoch, CaTchan)
CaT_m2h0 = smodel.ChanState('CaT_m2h0', mdl_stoch, CaTchan)
CaT_m2h1 = smodel.ChanState('CaT_m2h1', mdl_stoch, CaTchan)

CaTm0h0_m1h0 = smodel.VDepSReac('CaTm0h0_m1h0', ssys_stoch, slhs = [CaT_m0h0], srhs = [CaT_m1h0], k= lambda V: 1.0e3 *2.* alpham_cat(V*1.0e3))
CaTm1h0_m2h0 = smodel.VDepSReac('CaTm1h0_m2h0', ssys_stoch, slhs = [CaT_m1h0], srhs = [CaT_m2h0], k= lambda V: 1.0e3 *1.* alpham_cat(V*1.0e3))

CaTm2h0_m1h0 = smodel.VDepSReac('CaTm2h0_m1h0', ssys_stoch, slhs = [CaT_m2h0], srhs = [CaT_m1h0], k= lambda V: 1.0e3 *2.* betam_cat(V*1.0e3))
CaTm1h0_m0h0 = smodel.VDepSReac('CaTm1h0_m0h0', ssys_stoch, slhs = [CaT_m1h0], srhs = [CaT_m0h0], k= lambda V: 1.0e3 *1.* betam_cat(V*1.0e3))

CaTm0h1_m1h1 = smodel.VDepSReac('CaTm0h1_m1h1', ssys_stoch, slhs = [CaT_m0h1], srhs = [CaT_m1h1], k= lambda V: 1.0e3 *2.* alpham_cat(V*1.0e3))
CaTm1h1_m2h1 = smodel.VDepSReac('CaTm1h1_m2h1', ssys_stoch, slhs = [CaT_m1h1], srhs = [CaT_m2h1], k= lambda V: 1.0e3 *1.* alpham_cat(V*1.0e3))

CaTm2h1_m1h1 = smodel.VDepSReac('CaTm2h1_m1h1', ssys_stoch, slhs = [CaT_m2h1], srhs = [CaT_m1h1], k= lambda V: 1.0e3 *2.* betam_cat(V*1.0e3))
CaTm1h1_m0h1 = smodel.VDepSReac('CaTm1h1_m0h1', ssys_stoch, slhs = [CaT_m1h1], srhs = [CaT_m0h1], k= lambda V: 1.0e3 *1.* betam_cat(V*1.0e3))


CaTm0h0_m0h1 = smodel.VDepSReac('CaTm0h0_m0h1', ssys_stoch, slhs = [CaT_m0h0], srhs = [CaT_m0h1], k= lambda V: 1.0e3 *1.* alphah_cat(V*1.0e3))
CaTm1h0_m1h1 = smodel.VDepSReac('CaTm1h0_m1h1', ssys_stoch, slhs = [CaT_m1h0], srhs = [CaT_m1h1], k= lambda V: 1.0e3 *1.* alphah_cat(V*1.0e3))
CaTm2h0_m2h1 = smodel.VDepSReac('CaTm2h0_m2h1', ssys_stoch, slhs = [CaT_m2h0], srhs = [CaT_m2h1], k= lambda V: 1.0e3 *1.* alphah_cat(V*1.0e3))

CaTm2h1_m2h0 = smodel.VDepSReac('CaTm2h1_m2h0', ssys_stoch, slhs = [CaT_m2h1], srhs = [CaT_m2h0], k= lambda V: 1.0e3 *1.* betah_cat(V*1.0e3))
CaTm1h1_m1h0 = smodel.VDepSReac('CaTm1h1_m1h0', ssys_stoch, slhs = [CaT_m1h1], srhs = [CaT_m1h0], k= lambda V: 1.0e3 *1.* betah_cat(V*1.0e3))
CaTm0h1_m0h0 = smodel.VDepSReac('CaTm0h1_m0h0', ssys_stoch, slhs = [CaT_m0h1], srhs = [CaT_m0h0], k= lambda V: 1.0e3 *1.* betah_cat(V*1.0e3))

if cyl160:
    OC_CaT = smodel.GHKcurr('OC_CaT', ssys_stoch, CaT_m2h1, Ca_stoch, virtual_oconc = Ca_oconc, computeflux = True)
else:
    OC_CaT = smodel.GHKcurr('OC_CaT', ssys_stoch, CaT_m2h1, Ca_stoch, computeflux = True)
OC_CaT.setP(CaT_P)


##### BK channel ####################

BK_C0 = smodel.Spec('BK_C0', mdl_det)
BK_C1 = smodel.Spec('BK_C1', mdl_det)
BK_C2 = smodel.Spec('BK_C2', mdl_det)
BK_C3 = smodel.Spec('BK_C3', mdl_det)
BK_C4 = smodel.Spec('BK_C4', mdl_det)
BK_O0 = smodel.Spec('BK_O0', mdl_det)
BK_O1 = smodel.Spec('BK_O1', mdl_det)
BK_O2 = smodel.Spec('BK_O2', mdl_det)
BK_O3 = smodel.Spec('BK_O3', mdl_det)
BK_O4 = smodel.Spec('BK_O4', mdl_det)


BKCAC0 = smodel.SReac('BKCAC0', ssys_det, slhs = [BK_C0], ilhs = [Ca_det], srhs = [BK_C1], kcst = c_01)
BKCAC1 = smodel.SReac('BKCAC1', ssys_det, slhs = [BK_C1], ilhs = [Ca_det], srhs = [BK_C2], kcst = c_12)
BKCAC2 = smodel.SReac('BKCAC2', ssys_det, slhs = [BK_C2], ilhs = [Ca_det], srhs = [BK_C3], kcst = c_23)
BKCAC3 = smodel.SReac('BKCAC3', ssys_det, slhs = [BK_C3], ilhs = [Ca_det], srhs = [BK_C4], kcst = c_34)

BKC0 = smodel.SReac('BKC0', ssys_det, slhs = [BK_C1], srhs = [BK_C0], irhs=[Ca_det], kcst = c_10)
BKC1 = smodel.SReac('BKC1', ssys_det, slhs = [BK_C2], srhs = [BK_C1], irhs=[Ca_det], kcst = c_21)
BKC2 = smodel.SReac('BKC2', ssys_det, slhs = [BK_C3], srhs = [BK_C2], irhs=[Ca_det], kcst = c_32)
BKC3 = smodel.SReac('BKC3', ssys_det, slhs = [BK_C4], srhs = [BK_C3], irhs=[Ca_det], kcst = c_43)

BKCAO0 = smodel.SReac('BKCAO0', ssys_det, slhs = [BK_O0], ilhs = [Ca_det], srhs = [BK_O1], kcst = o_01)
BKCAO1 = smodel.SReac('BKCAO1', ssys_det, slhs = [BK_O1], ilhs = [Ca_det], srhs = [BK_O2], kcst = o_12)
BKCAO2 = smodel.SReac('BKCAO2', ssys_det, slhs = [BK_O2], ilhs = [Ca_det], srhs = [BK_O3], kcst = o_23)
BKCAO3 = smodel.SReac('BKCAO3', ssys_det, slhs = [BK_O3], ilhs = [Ca_det], srhs = [BK_O4], kcst = o_34)

BKO0 = smodel.SReac('BKO0', ssys_det, slhs = [BK_O1], srhs = [BK_O0], irhs=[Ca_det], kcst = o_10)
BKO1 = smodel.SReac('BKO1', ssys_det, slhs = [BK_O2], srhs = [BK_O1], irhs=[Ca_det], kcst = o_21)
BKO2 = smodel.SReac('BKO2', ssys_det, slhs = [BK_O3], srhs = [BK_O2], irhs=[Ca_det], kcst = o_32)
BKO3 = smodel.SReac('BKO3', ssys_det, slhs = [BK_O4], srhs = [BK_O3], irhs=[Ca_det], kcst = o_43)

BKC0O0 = smodel.SReac('BKC0O0', ssys_det, slhs = [BK_C0], srhs = [BK_O0], kcst=0.0)
BKC1O1 = smodel.SReac('BKC1O1', ssys_det, slhs = [BK_C1], srhs = [BK_O1], kcst=0.0)
BKC2O2 = smodel.SReac('BKC2O2', ssys_det, slhs = [BK_C2], srhs = [BK_O2], kcst=0.0)
BKC3O3 = smodel.SReac('BKC3O3', ssys_det, slhs = [BK_C3], srhs = [BK_O3], kcst=0.0)
BKC4O4 = smodel.SReac('BKC4O4', ssys_det, slhs = [BK_C4], srhs = [BK_O4], kcst=0.0)

BKO0C0 = smodel.SReac('BKO0C0', ssys_det, slhs = [BK_O0], srhs = [BK_C0], kcst=0.0)
BKO1C1 = smodel.SReac('BKO1C1', ssys_det, slhs = [BK_O1], srhs = [BK_C1], kcst=0.0)
BKO2C2 = smodel.SReac('BKO2C2', ssys_det, slhs = [BK_O2], srhs = [BK_C2], kcst=0.0)
BKO3C3 = smodel.SReac('BKO3C3', ssys_det, slhs = [BK_O3], srhs = [BK_C3], kcst=0.0)
BKO4C4 = smodel.SReac('BKO4C4', ssys_det, slhs = [BK_O4], srhs = [BK_C4], kcst=0.0)


###### SK channel ################## DETERMINISTIC

SK_C1 = smodel.Spec('SK_C1', mdl_det)
SK_C2 = smodel.Spec('SK_C2', mdl_det)
SK_C3 = smodel.Spec('SK_C3', mdl_det)
SK_C4 = smodel.Spec('SK_C4', mdl_det)
SK_O1 = smodel.Spec('SK_O1', mdl_det)
SK_O2 = smodel.Spec('SK_O2', mdl_det)


SKCAC1 = smodel.SReac('SKCAC1', ssys_det, slhs = [SK_C1], ilhs = [Ca_det], srhs = [SK_C2], kcst = dirc2_t)
SKCAC2 = smodel.SReac('SKCAC2', ssys_det, slhs = [SK_C2], ilhs = [Ca_det], srhs = [SK_C3], kcst = dirc3_t)
SKCAC3 = smodel.SReac('SKCAC3', ssys_det, slhs = [SK_C3], ilhs = [Ca_det], srhs = [SK_C4], kcst = dirc4_t)

SKC1 = smodel.SReac('SKC1', ssys_det, slhs = [SK_C2], srhs = [SK_C1], irhs=[Ca_det], kcst = invc1_t)
SKC2 = smodel.SReac('SKC2', ssys_det, slhs = [SK_C3], srhs = [SK_C2], irhs=[Ca_det], kcst = invc2_t)
SKC3 = smodel.SReac('SKC3', ssys_det, slhs = [SK_C4], srhs = [SK_C3], irhs=[Ca_det], kcst = invc3_t)

SKC3O1 = smodel.SReac('SKC3O1', ssys_det, slhs = [SK_C3], srhs = [SK_O1], kcst = diro1_t)
SKC4O2 = smodel.SReac('SKC4O2', ssys_det, slhs = [SK_C4], srhs = [SK_O2], kcst = diro2_t)

SKO1C3 = smodel.SReac('SKO1C3', ssys_det, slhs = [SK_O1], srhs = [SK_C3], kcst = invo1_t)
SKO2C4 = smodel.SReac('SKO2C4', ssys_det, slhs = [SK_O2], srhs = [SK_C4], kcst = invo2_t)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

########### MESH & COMPARTMENTALIZATION #################

##########Import Mesh

# For stochastic sim:
mesh_stoch = meshio.loadMesh('./meshes/'+meshfile_ab)[0]


# For determinstic sim:
mesh_det = meshio.loadMesh('./meshes/'+meshfile_ab)[0]


outer_tets = range(mesh_stoch.ntets)
inner_tets = gettets.getcyl(mesh_stoch, 1e-6, -200e-6, 200e-6)[0]

for i in inner_tets: outer_tets.remove(i)

print outer_tets.__len__(), " tets in outer compartment"
print inner_tets.__len__(), " tets in inner compartment"

# Record voltage from the central tetrahedron
cent_tet = mesh_stoch.findTetByPoint([0.0,0.0,0.0])

########## Create an intracellular compartment i.e. cytosolic compartment

cyto_stoch = sgeom.TmComp('cyto_stoch', mesh_stoch, inner_tets)
cyto_stoch.addVolsys('vsys_stoch')

if not cyl160: outer_stoch = sgeom.TmComp('outer_stoch', mesh_stoch, outer_tets)

cyto_det = sgeom.TmComp('cyto_det', mesh_det, inner_tets)
cyto_det.addVolsys('vsys_det')

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

########## Find the submembrane tets

memb_tet_neighb = []
for i in memb_tris:
    tettemp = mesh_stoch.getTriTetNeighb(i)
    for j in tettemp:
        memb_tet_neighb.append(j)

submemb_tets = []
for i in memb_tet_neighb:
    if i in inner_tets:
        submemb_tets.append(i)

print len(submemb_tets)

vol = 0.0

for i in submemb_tets:
    vol = vol + mesh_stoch.getTetVol(i)

print 'Volume of submembrane region is', vol

submemb_tets_surftris = dict()

for m in submemb_tets:
    tris = mesh_stoch.getTetTriNeighb(m)
    for t in tris:
        if t in memb_tris:
            submemb_tets_surftris[m] = t
            break

assert(len(submemb_tets_surftris.values()) == len(submemb_tets))
					    

########## Create a membrane as a surface mesh

# Stochastic sim:
if cyl160: 
    memb_stoch = sgeom.TmPatch('memb_stoch', mesh_stoch, memb_tris, cyto_stoch)
else:
    memb_stoch = sgeom.TmPatch('memb_stoch', mesh_stoch, memb_tris, cyto_stoch, outer_stoch)
memb_stoch.addSurfsys('ssys_stoch')

# Determinsitic sim:

memb_det = sgeom.TmPatch('memb_det', mesh_det, memb_tris, cyto_det)
memb_det.addSurfsys('ssys_det')


# For EField calculation
print "Creating membrane.."
membrane = sgeom.Memb('membrane', mesh_stoch, [memb_stoch])
print "Membrane created."

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

r = srng.create_mt19937(512)
r.initialize(7)

r_dummy = srng.create_mt19937(512)
r_dummy.initialize(7)

print "Creating Tet exact solver"

#Creating two solvers
sim_stoch = ssolver.Tetexact(mdl_stoch, mesh_stoch, r, True)

print "Creating Tet ODE solver"

sim_det = ssolver.TetODE(mdl_det, mesh_det, r_dummy)

sim_det.setTolerances(1.0e-3, 1.0e-3)


print "Resetting simulation objects.."
sim_stoch.reset()

print "Injecting molecules.."

sim_stoch.setTemp(TEMPERATURE+273.15)

if not cyl160:
    sim_stoch.setCompConc('outer_stoch', 'Ca_stoch', Ca_oconc)
    sim_stoch.setCompClamped('outer_stoch', 'Ca_stoch', True)

sim_stoch.setCompConc('cyto_stoch', 'Ca_stoch', Ca_iconc)

print "Calcium concentration in stochastic simulation is: ", sim_stoch.getCompConc('cyto_stoch', 'Ca_stoch')
print "No. of Ca molecules in stochastic simulation is: ", sim_stoch.getCompCount('cyto_stoch', 'Ca_stoch')

sim_det.setCompConc('cyto_det', 'Ca_det', Ca_iconc)

print "Calcium concentration in deterministic simulation is: ", sim_det.getCompConc('cyto_det', 'Ca_det')
print "No. of Ca molecules in deterministic simulation is: ", sim_det.getCompCount('cyto_det', 'Ca_det')

sim_det.setCompConc('cyto_det', 'Mg', Mg_conc)


surfarea = sim_stoch.getPatchArea('memb_stoch')

pumpnbs = 6.022141e12*surfarea

sim_det.setPatchCount('memb_det', 'Pump', round(pumpnbs))
sim_det.setPatchCount('memb_det', 'CaPump', 0)

print "Injected ", sim_det.getPatchCount('memb_det', 'Pump'), "pumps"

sim_det.setCompConc('cyto_det', 'iCBsf',    iCBsf_conc)
sim_det.setCompConc('cyto_det', 'iCBsCa',   iCBsCa_conc)
sim_det.setCompConc('cyto_det', 'iCBCaf',   iCBCaf_conc)
sim_det.setCompConc('cyto_det', 'iCBCaCa',  iCBCaCa_conc)

sim_det.setCompConc('cyto_det', 'CBsf',     CBsf_conc)
sim_det.setCompConc('cyto_det', 'CBsCa',    CBsCa_conc)
sim_det.setCompConc('cyto_det', 'CBCaf',    CBCaf_conc)
sim_det.setCompConc('cyto_det', 'CBCaCa',   CBCaCa_conc)

sim_det.setCompConc('cyto_det', 'PV',       PV_conc)
sim_det.setCompConc('cyto_det', 'PVCa',     PVCa_conc)
sim_det.setCompConc('cyto_det', 'PVMg',     PVMg_conc)


# CaP

sim_det.setPatchCount('memb_det', 'CaP_m0' , round(CaP_ro*surfarea*CaP_m0_p))
sim_det.setPatchCount('memb_det', 'CaP_m1' , round(CaP_ro*surfarea*CaP_m1_p))
sim_det.setPatchCount('memb_det', 'CaP_m2' , round(CaP_ro*surfarea*CaP_m2_p))
sim_det.setPatchCount('memb_det', 'CaP_m3' , round(CaP_ro*surfarea*CaP_m3_p))


print "CaP_m0 ", round(CaP_ro*surfarea*CaP_m0_p)
print "CaP_m1 ", round(CaP_ro*surfarea*CaP_m1_p)
print "CaP_m2 ", round(CaP_ro*surfarea*CaP_m2_p)
print "CaP_m3 ", round(CaP_ro*surfarea*CaP_m3_p)

print "Targeted Injection: ", round(CaP_ro*surfarea), "CaP channels"

# CaT

# From cstate: CaT_m2h0 conducting
sim_stoch.setPatchCount('memb_stoch', 'CaT_m0h0' , round(CaT_ro*surfarea*CaT_m0h0_p))
sim_stoch.setPatchCount('memb_stoch', 'CaT_m1h0' , round(CaT_ro*surfarea*CaT_m1h0_p))
sim_stoch.setPatchCount('memb_stoch', 'CaT_m2h0' , round(CaT_ro*surfarea*CaT_m2h0_p))
sim_stoch.setPatchCount('memb_stoch', 'CaT_m0h1' , round(CaT_ro*surfarea*CaT_m0h1_p))
sim_stoch.setPatchCount('memb_stoch', 'CaT_m1h1' , round(CaT_ro*surfarea*CaT_m1h1_p))
sim_stoch.setPatchCount('memb_stoch', 'CaT_m2h1' , round(CaT_ro*surfarea*CaT_m2h1_p))

print "CaT_m0h0 ", round(CaT_ro*surfarea*CaT_m0h0_p)
print "CaT_m1h0 ", round(CaT_ro*surfarea*CaT_m1h0_p)
print "CaT_m2h0 ", round(CaT_ro*surfarea*CaT_m2h0_p)

print "CaT_m0h1 ", round(CaT_ro*surfarea*CaT_m0h1_p)
print "CaT_m1h1 ", round(CaT_ro*surfarea*CaT_m1h1_p)
print "CaT_m2h1 ", round(CaT_ro*surfarea*CaT_m2h1_p)

print " Target Injection  ", CaT_ro*surfarea, "CaT channels"

# BK
sim_det.setPatchCount('memb_det', 'BK_C0' , round(BK_ro*surfarea*BK_C0_p))
sim_det.setPatchCount('memb_det', 'BK_C1' , round(BK_ro*surfarea*BK_C1_p))
sim_det.setPatchCount('memb_det', 'BK_C2' , round(BK_ro*surfarea*BK_C2_p))
sim_det.setPatchCount('memb_det', 'BK_C3' , round(BK_ro*surfarea*BK_C3_p))
sim_det.setPatchCount('memb_det', 'BK_C4' , round(BK_ro*surfarea*BK_C4_p))

sim_det.setPatchCount('memb_det', 'BK_O0' , round(BK_ro*surfarea*BK_O0_p))
sim_det.setPatchCount('memb_det', 'BK_O1' , round(BK_ro*surfarea*BK_O1_p))
sim_det.setPatchCount('memb_det', 'BK_O2' , round(BK_ro*surfarea*BK_O2_p))
sim_det.setPatchCount('memb_det', 'BK_O3' , round(BK_ro*surfarea*BK_O3_p))
sim_det.setPatchCount('memb_det', 'BK_O4' , round(BK_ro*surfarea*BK_O4_p))


print "Injected  ", BK_ro*surfarea, "BK channels"

# SK
sim_det.setPatchCount('memb_det', 'SK_C1' , round(SK_ro*surfarea*SK_C1_p))
sim_det.setPatchCount('memb_det', 'SK_C2' , round(SK_ro*surfarea*SK_C2_p))
sim_det.setPatchCount('memb_det', 'SK_C3' , round(SK_ro*surfarea*SK_C3_p))
sim_det.setPatchCount('memb_det', 'SK_C4' , round(SK_ro*surfarea*SK_C4_p))

sim_det.setPatchCount('memb_det', 'SK_O1' , round(SK_ro*surfarea*SK_O1_p))
sim_det.setPatchCount('memb_det', 'SK_O2' , round(SK_ro*surfarea*SK_O2_p))


print "Injected ", SK_ro*surfarea, "SK channels"


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
try: os.mkdir(root+'data/' +  'HybridCaburst_stochCaT')
except: pass
try: os.mkdir(root+'data/' +  'HybridCaburst_stochCaT/'+meshfile_ab)
except: pass 

os.mkdir(root+'data/' +  'HybridCaburst_stochCaT/'+meshfile_ab+'/'+iter_n+'__'+dc )


datfile =  open(root+'data/' +  'HybridCaburst_stochCaT/'+meshfile_ab+'/'+iter_n+'__'+dc + '/currents.dat', 'w')
datfile2 = open(root+'data/' +  'HybridCaburst_stochCaT/'+meshfile_ab+'/'+iter_n+'__'+dc + '/voltage.dat', 'w')
datfile3 = open(root+'data/' +  'HybridCaburst_stochCaT/'+meshfile_ab+'/'+iter_n+'__'+dc + '/calcium.dat', 'w')


r.initialize(int(time.time()%1000))


btime = time.time()
for l in range(NTIMEPOINTS):
    print "Tpnt: ", l

    #1) RUN STOCHASTIC SIMULATION i.e. compute currents and update stochastic calcium concentration
    sim_stoch.run(TIMECONVERTER*l)

    #2) READ STOCHASTIC CA and 3) SET DETERMINISTIC CA AND RATE CONSTANTS FOR DETERMINISTIC CHANNELS
    for m in submemb_tets:
        Si = sim_stoch.getTetConc(m,'Ca_stoch')
        sim_det.setTetConc(m,'Ca_det',Si)

    #Assuming this sim V is not constant everwhere
    for m in submemb_tets:
        ctriID = submemb_tets_surftris[m]

        V = sim_stoch.getTriV(ctriID)

        #3) Set the rate constants and RUN THE DETERMINISTIC SIMULATION

        sim_det.setTriSReacK(ctriID,'CaPm0m1', 1.0e3 *3.* alpha_cap(V*1.0e3)*Qt)
        sim_det.setTriSReacK(ctriID,'CaPm1m2', 1.0e3 *2.* alpha_cap(V*1.0e3)*Qt)
        sim_det.setTriSReacK(ctriID,'CaPm2m3', 1.0e3 *1.* alpha_cap(V*1.0e3)*Qt)

        sim_det.setTriSReacK(ctriID,'CaPm3m2', 1.0e3 *3.* beta_cap(V*1.0e3)*Qt)
        sim_det.setTriSReacK(ctriID,'CaPm2m1', 1.0e3 *2.* beta_cap(V*1.0e3)*Qt)
        sim_det.setTriSReacK(ctriID,'CaPm1m0', 1.0e3 *1.* beta_cap(V*1.0e3)*Qt)

        sim_det.setTriSReacK(ctriID, 'BKC0O0', f_0(V))
        sim_det.setTriSReacK(ctriID, 'BKC1O1', f_1(V))
        sim_det.setTriSReacK(ctriID, 'BKC2O2', f_2(V))
        sim_det.setTriSReacK(ctriID, 'BKC3O3', f_3(V))
        sim_det.setTriSReacK(ctriID, 'BKC4O4', f_4(V))
        sim_det.setTriSReacK(ctriID, 'BKO0C0', b_0(V))
        sim_det.setTriSReacK(ctriID, 'BKO1C1', b_1(V))
        sim_det.setTriSReacK(ctriID, 'BKO2C2', b_2(V))
        sim_det.setTriSReacK(ctriID, 'BKO3C3', b_3(V))
        sim_det.setTriSReacK(ctriID, 'BKO4C4', b_4(V))

    #4) RUN DETERMINISTIC SIMULATION
    sim_det.run(TIMECONVERTER*l)

    # Now do the communication between the sims

    #5)READ DETERMINISTIC CHANNELS & THEN COMPUTE CURRENT USING DETERMINISTIC GHK (could be stochastic)
    So = Ca_oconc
    # i) For each tet in submembrane, find the corresponding triID
    # ii) For each tri, compute GHK current for each channel
    # iii) Count the channel states / Spec in open states for each of the triID and compute the total current of that channel

    tcur_CaP = 0.0
    tcur_CaT = 0.0
    tcur_BK = 0.0
    tcur_SK = 0.0

    tca_count_det = 0.0
    tca_count_stoch = 0.0
    
    for m in submemb_tets:
        ctriID = submemb_tets_surftris[m]
        V = sim_stoch.getTriV(ctriID)
        Si = sim_det.getTetConc(m,'Ca_det')
        cur_CaP_sc = cf.getGHKI(CaP_P, V, 2, TEMPERATURE+273.15, Si*1.0e3, So*1.0e3)
        cur_BK_sc = cf.getOhmI(V, BK_rev, BK_G)
        cur_SK_sc = cf.getOhmI(V, SK_rev, SK_G)
        cur_L_sc = cf.getOhmI(V, L_rev, L_G)

        cur_CaP = cur_CaP_sc*(sim_det.getTriCount(ctriID, 'CaP_m3'))
        cur_BK =  cur_BK_sc*(sim_det.getTriCount(ctriID, 'BK_O0') + sim_det.getTriCount(ctriID, 'BK_O1') + sim_det.getTriCount(ctriID, 'BK_O2') + sim_det.getTriCount(ctriID, 'BK_O3') + sim_det.getTriCount(ctriID, 'BK_O4'))
        cur_SK = cur_SK_sc*(sim_det.getTriCount(ctriID, 'SK_O1') + sim_det.getTriCount(ctriID, 'SK_O2'))
        #cur_L corresponding to each surftri has been corrected in the following script line
        cur_L = cur_L_sc*(round(L_ro * sim_det.getPatchArea('memb_det')))*(sim_stoch.getTriArea(ctriID)/sim_det.getPatchArea('memb_det'))

        ca_count_inj = -1.0*((cur_CaP)*TIMECONVERTER)/(2*E_CHARGE)
        sim_stoch.setTetCount(m, 'Ca_stoch', ca_count_inj+sim_det.getTetCount(m,'Ca_det'))
        sim_stoch.setTriIClamp(ctriID, cur_CaP+cur_BK+cur_SK+cur_L)
	
        tcur_CaT = tcur_CaT +  sim_stoch.getTriGHKI(ctriID,'OC_CaT')
        tcur_CaP = tcur_CaP + cur_CaP
        tcur_BK = tcur_BK + cur_BK
        tcur_SK = tcur_SK + cur_SK
        tca_count_det = tca_count_det + sim_det.getTetCount(m,'Ca_det')
        tca_count_stoch = tca_count_stoch + sim_stoch.getTetCount(m,'Ca_stoch')
        
    
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
    datfile3.write('%.6g' %(((tca_count_det/AVOGADRO)/(vol*1.0e3))*1.0e6) + ' ')
    datfile3.write('%.6g' %(((tca_count_stoch/AVOGADRO)/(vol*1.0e3))*1.0e6) + ' ')
    datfile3.write('%.6g' %tca_count_det + ' ')
    datfile3.write('%.6g' %tca_count_stoch + ' ')
    datfile3.write('\n')
    

datfile.close()
datfile2.close()
datfile3.close()


