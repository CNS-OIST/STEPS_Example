#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

import steps.interface

from steps.model import *
from extra.constants import *

########################### BIOCHEMICAL MODEL ###############################
def getModel():
    mdl = Model()
    with mdl:

        # Species
        Ca, Pump, CaPump, iCBsf, iCBsCa, iCBCaf, iCBCaCa, CBsf, CBsCa, CBCaf, CBCaCa, PV, PVMg, PVCa, Mg = Species.Create()

        # Vol/surface systems
        vsys = VolumeSystem.Create()
        ssys = SurfaceSystem.Create()

        with vsys:

            diff_Ca =     Diffusion.Create(Ca, DCST)
            diff_CBsf =   Diffusion.Create(CBsf, DCB)
            diff_CBsCa =  Diffusion.Create(CBsCa, DCB)
            diff_CBCaf =  Diffusion.Create(CBCaf, DCB)
            diff_CBCaCa = Diffusion.Create(CBCaCa, DCB)
            diff_PV =     Diffusion.Create(PV, DPV)
            diff_PVCa =   Diffusion.Create(PVCa, DPV)
            diff_PVMg =   Diffusion.Create(PVMg, DPV)

        with ssys:

            #Pump
            Ca.i + Pump.s <r['PumpD_f']> CaPump.s
            r['PumpD_f'].setRates(P_f_kcst, P_b_kcst)

            CaPump.s >r['PumpD_k']> Pump.s
            r['PumpD_k'].setRates(P_k_kcst)

        with vsys, mdl:

            #iCBsf-fast and iCBsCa
            (Ca + iCBsf <r['iCBsf1_f']> iCBsCa) + Ca <r['iCBsCa_f']> iCBCaCa
            r['iCBsf1_f'].setRates(iCBsf1_f_kcst, iCBsf1_b_kcst)
            r['iCBsCa_f'].setRates(iCBsCa_f_kcst, iCBsCa_b_kcst)

            #iCBsf_slow and iCBCaf
            (Ca + iCBsf <r['iCBsf2_f']> iCBCaf) + Ca <r['iCBCaf_f']> iCBCaCa
            r['iCBsf2_f'].setRates(iCBsf2_f_kcst, iCBsf2_b_kcst)
            r['iCBCaf_f'].setRates(iCBCaf_f_kcst, iCBCaf_b_kcst)

            #CBsf-fast and CBsCa
            (CBsf + Ca <r['CBsf1_f']> CBsCa) + Ca <r['CBsCa_f']> CBCaCa
            r['CBsf1_f'].setRates(CBsf1_f_kcst, CBsf1_b_kcst)
            r['CBsCa_f'].setRates(CBsCa_f_kcst, CBsCa_b_kcst)

            #CBsf_slow and CBCaf
            (CBsf + Ca <r['CBsf2_f']> CBCaf) + Ca <r['CBCaf_f']> CBCaCa
            r['CBsf2_f'].setRates(CBsf2_f_kcst, CBsf2_b_kcst)
            r['CBCaf_f'].setRates(CBCaf_f_kcst, CBCaf_b_kcst)

            #PVca
            Ca + PV <r['PVca_f']> PVCa
            r['PVca_f'].setRates(PVca_f_kcst, PVca_b_kcst)

            #PVmg
            Mg + PV <r['PVmg_f']> PVMg
            r['PVmg_f'].setRates(PVmg_f_kcst, PVmg_b_kcst)

            # Ca Influx converted from P Type current
            None >r['CaInflux']> Ca
            r['CaInflux'].setRates(0.0)

    return mdl

