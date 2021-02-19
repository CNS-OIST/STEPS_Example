#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################


import steps.model as smodel
from extra.constants import *

########################### BIOCHEMICAL MODEL ###############################
def getModel():
    mdl = smodel.Model()

    # Calcium
    Ca = smodel.Spec('Ca', mdl)

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

    # Ca Influx converted from P Type current
    CaInflux = smodel.Reac('CaInflux', vsys, lhs = [], rhs = [Ca], kcst = 0.0)

    return mdl

