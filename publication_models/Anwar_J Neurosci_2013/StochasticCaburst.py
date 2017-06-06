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
# *StochasticCaburst.py : The spatial stochastic calcium burst model, used in the 
# above study. 
#
# Script authors: Haroon Anwar and Iain Hepburn
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# USAGE
#
# $ python StochasticCaburst.py *mesh* *root* *iter_n* 
#  
#  *mesh* is the tetrahedral mesh (10um to 160um cylinder)
#  *root* is the path to the location for data storage 
#  *iter_n* (is intened to be an integer) is an identifier number for each 
#     simulation iteration. iter_n is also used to initialize the random
#     number generator.
#
# E.g: 
# $ python StochasticCaburst.py Cylinder2_dia2um_L10um_outer0_3um_0.3shell_0.3size_19156tets_adaptive.inp ~/stochcasims/ 1
#
#
# OUTPUT 
#
# In (root)/data/StochasticCaburst/(mesh)/(iter_n+time) directory 
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
# Time (ms), calcium concentration in submembrane (micromolar),
# number of calcium ions in submembrane. 
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

import meshes.gettets as gettets
from extra.constants import *

import sys
import os

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

meshfile_ab, root, iter_n = sys.argv[1], sys.argv[2], sys.argv[3]

if meshfile_ab == 'Cylinder2_dia2um_L160um_outer0_0.3shell_0.3size_279152tets_adaptive.inp': cyl160=True
else: cyl160=False

########################### BIOCHEMICAL MODEL ###############################

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
iCBsf1_b = smodel.Reac('iCBsf1_b', vsys, lhs=[iCBsCa], rhs=[Ca,iCBsf], kcst = iCBsf1_b_kcst)

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


##### BK channel ################################


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


###### SK channel ##################


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

##################################

########### MESH & COMPARTMENTALIZATION #################

##########Import Mesh

mesh = meshio.loadMesh('./meshes/'+meshfile_ab)[0]

outer_tets = range(mesh.ntets)

###USE OF gettets
#getcyl(tetmesh, rad,  zmin, zmax, binnum=120, x = 0.0, y = 0.0):

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

print len(memb_tris), " surface triangles."

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


########## Create a membrane as a surface mesh
if cyl160: 
    memb = sgeom.TmPatch('memb', mesh, memb_tris, cyto)
else:
    memb = sgeom.TmPatch('memb', mesh, memb_tris, cyto, outer)

memb.addSurfsys('ssys')

print "Area: ", memb.getArea()

print "Creating membrane.."
membrane = sgeom.Memb('membrane', mesh, [memb], opt_file_name = './meshes/'+meshfile_ab+"_optimalidx")
print "Membrane created."

# # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

r = srng.create_mt19937(512)
r.initialize(7)

sim = ssolver.Tetexact(mdl, mesh, r, True)

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

#Total pump is 1e-15 mol/cm2 ---> 1e-11 mol/m2
#pumpnbs per unit area (im m2) is Total pump times AVOGADRO's NUMBER (1e-11 mol/m2 * 6.022e23 /mol )
pumpnbs = 6.022141e12*surfarea
print "Number of pump molecules: ", pumpnbs

sim.setPatchCount('memb', 'Pump', round(pumpnbs))
sim.setPatchCount('memb', 'CaPump', 0)

sim.setCompConc('cyto', 'iCBsf', iCBsf_conc)
sim.setCompConc('cyto', 'iCBCaf', iCBCaf_conc)
sim.setCompConc('cyto', 'iCBsCa', iCBsCa_conc)
sim.setCompConc('cyto', 'iCBCaCa', iCBCaCa_conc)

sim.setCompConc('cyto', 'CBsf', CBsf_conc)
sim.setCompConc('cyto', 'CBCaf', CBCaf_conc)
sim.setCompConc('cyto', 'CBsCa', CBsCa_conc)
sim.setCompConc('cyto', 'CBCaCa', CBCaCa_conc)

sim.setCompConc('cyto', 'PV', PV_conc)
sim.setCompConc('cyto', 'PVCa', PVCa_conc)
sim.setCompConc('cyto', 'PVMg', PVMg_conc)

sim.setPatchCount('memb', 'CaP_m0' , round(CaP_ro*surfarea*CaP_m0_p))
sim.setPatchCount('memb', 'CaP_m1' , round(CaP_ro*surfarea*CaP_m1_p))
sim.setPatchCount('memb', 'CaP_m2' , round(CaP_ro*surfarea*CaP_m2_p))
sim.setPatchCount('memb', 'CaP_m3' , round(CaP_ro*surfarea*CaP_m3_p))


print "CaP_m0 ", round(CaP_ro*surfarea*CaP_m0_p)
print "CaP_m1 ", round(CaP_ro*surfarea*CaP_m1_p)
print "CaP_m2 ", round(CaP_ro*surfarea*CaP_m2_p)
print "CaP_m3 ", round(CaP_ro*surfarea*CaP_m3_p)

print "Targeted Injection: ", round(CaP_ro*surfarea), "CaP channels"

sim.setPatchCount('memb', 'CaT_m0h0' , round(CaT_ro*surfarea*CaT_m0h0_p))
sim.setPatchCount('memb', 'CaT_m1h0' , round(CaT_ro*surfarea*CaT_m1h0_p))
sim.setPatchCount('memb', 'CaT_m2h0' , round(CaT_ro*surfarea*CaT_m2h0_p))
sim.setPatchCount('memb', 'CaT_m0h1' , round(CaT_ro*surfarea*CaT_m0h1_p))
sim.setPatchCount('memb', 'CaT_m1h1' , round(CaT_ro*surfarea*CaT_m1h1_p))
sim.setPatchCount('memb', 'CaT_m2h1' , round(CaT_ro*surfarea*CaT_m2h1_p))

print "m0h0", round(CaT_ro*surfarea*CaT_m0h0_p)
print "m1h0", round(CaT_ro*surfarea*CaT_m1h0_p)
print "m2h0", round(CaT_ro*surfarea*CaT_m2h0_p)
print "m0h1", round(CaT_ro*surfarea*CaT_m0h1_p)
print "m1h1", round(CaT_ro*surfarea*CaT_m1h1_p)
print "m2h1", round(CaT_ro*surfarea*CaT_m2h1_p)

print "Targeted Injection: ", round(CaT_ro*surfarea), "CaT channels"

sim.setPatchCount('memb', 'BK_C0' , round(BK_ro*surfarea*BK_C0_p))
sim.setPatchCount('memb', 'BK_C1' , round(BK_ro*surfarea*BK_C1_p))
sim.setPatchCount('memb', 'BK_C2' , round(BK_ro*surfarea*BK_C2_p))
sim.setPatchCount('memb', 'BK_C3' , round(BK_ro*surfarea*BK_C3_p))
sim.setPatchCount('memb', 'BK_C4' , round(BK_ro*surfarea*BK_C4_p))

sim.setPatchCount('memb', 'BK_O0' , round(BK_ro*surfarea*BK_O0_p))
sim.setPatchCount('memb', 'BK_O1' , round(BK_ro*surfarea*BK_O1_p))
sim.setPatchCount('memb', 'BK_O2' , round(BK_ro*surfarea*BK_O2_p))
sim.setPatchCount('memb', 'BK_O3' , round(BK_ro*surfarea*BK_O3_p))
sim.setPatchCount('memb', 'BK_O4' , round(BK_ro*surfarea*BK_O4_p))

print "BK_C0 ", round(BK_ro*surfarea*BK_C0_p)
print "BK_C1 ", round(BK_ro*surfarea*BK_C1_p)
print "BK_C2 ", round(BK_ro*surfarea*BK_C2_p)
print "BK_C3 ", round(BK_ro*surfarea*BK_C3_p)
print "BK_C4 ", round(BK_ro*surfarea*BK_C4_p)

print "BK_O0 ", round(BK_ro*surfarea*BK_O0_p)
print "BK_O1 ", round(BK_ro*surfarea*BK_O1_p)
print "BK_O2 ", round(BK_ro*surfarea*BK_O2_p)
print "BK_O3 ", round(BK_ro*surfarea*BK_O3_p)
print "BK_O4 ", round(BK_ro*surfarea*BK_O4_p)

print "Targeted Injection: ", round(BK_ro*surfarea), "BK channels"

sim.setPatchCount('memb', 'SK_C1' , round(SK_ro*surfarea*SK_C1_p))
sim.setPatchCount('memb', 'SK_C2' , round(SK_ro*surfarea*SK_C2_p))
sim.setPatchCount('memb', 'SK_C3' , round(SK_ro*surfarea*SK_C3_p))
sim.setPatchCount('memb', 'SK_C4' , round(SK_ro*surfarea*SK_C4_p))

sim.setPatchCount('memb', 'SK_O1' , round(SK_ro*surfarea*SK_O1_p))
sim.setPatchCount('memb', 'SK_O2' , round(SK_ro*surfarea*SK_O2_p))

print "SK_C1 ", round(SK_ro*surfarea*SK_C1_p)
print "SK_C2 ", round(SK_ro*surfarea*SK_C2_p)
print "SK_C3 ", round(SK_ro*surfarea*SK_C3_p)
print "SK_C4 ", round(SK_ro*surfarea*SK_C4_p)

print "SK_O1 ", round(SK_ro*surfarea*SK_O1_p)
print "SK_O2 ", round(SK_ro*surfarea*SK_O2_p)

print "Targeted Injection ", round(SK_ro*surfarea), "SK channels"

sim.setPatchCount('memb', 'Leak', int(L_ro * surfarea))
print "Injected ", int(L_ro * sim.getPatchArea('memb')), "Leak channels"


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
try: os.mkdir(root+'data/' +  'StochasticCaburst')
except: pass
try: os.mkdir(root+'data/' +  'StochasticCaburst/'+meshfile_ab)
except: pass 

os.mkdir(root+'data/' +  'StochasticCaburst/'+meshfile_ab+'/'+iter_n+'__'+dc )


datfile =  open(root+'data/' +  'StochasticCaburst/'+meshfile_ab+'/'+iter_n+'__'+dc + '/currents.dat', 'w')
datfile2 = open(root+'data/' +  'StochasticCaburst/'+meshfile_ab+'/'+iter_n+'__'+dc + '/voltage.dat', 'w')
datfile3 = open(root+'data/' +  'StochasticCaburst/'+meshfile_ab+'/'+iter_n+'__'+dc + '/calcium.dat', 'w')


r.initialize(100*int(iter_n))

for l in range(NTIMEPOINTS):
    print "Tpnt: ", l
    
    sim.run(TIMECONVERTER*l)

    tcur_CaP = 0.0
    tcur_CaT = 0.0
    tcur_BK = 0.0
    tcur_SK = 0.0
    tca_count = 0.0
    for m in submemb_tets:
        tca_count = tca_count + sim.getTetCount(m,'Ca')
    for m in memb_tris:
        tcur_CaP = tcur_CaP + sim.getTriGHKI(m,'OC_CaP') 
        tcur_CaT = tcur_CaT + sim.getTriGHKI(m,'OC_CaT')
        tcur_BK = tcur_BK + sim.getTriOhmicI(m,'OC_BK0') \
            + sim.getTriOhmicI(m,'OC_BK1') \
            + sim.getTriOhmicI(m,'OC_BK2') \
            + sim.getTriOhmicI(m,'OC_BK3') \
            + sim.getTriOhmicI(m,'OC_BK4')
        tcur_SK = tcur_SK + sim.getTriOhmicI(m,'OC1_SK') + sim.getTriOhmicI(m,'OC2_SK')
    
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
    datfile3.write('%.6g' %(sim.getCompConc('cyto', 'Ca')*1.0e6) + ' ')
    datfile3.write('%.6g' %tca_count + ' ')
    datfile3.write('\n')
    

datfile.close()
datfile2.close()
datfile3.close()
## END
