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
# *HybridCaburst_stochBK.py : A hybrid calcium burst model with stochastic BK 
# channels and everything else deterministic.
#
# Script authors: Haroon Anwar and Iain Hepburn
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# USAGE
#
# $ python HybridCaburst_stochBK.py *mesh* *root* *iter_n* 
#  
#  *mesh* is the tetrahedral mesh (10um to 160um cylinder)
#  *root* is the path to the location for data storage 
#  *iter_n* (is intened to be an integer) is an identifier number for each 
#     simulation iteration. 
#
# E.g: python HybridCaburst_stochBK.py Cylinder2_dia2um_L10um_outer0_3um_0.3shell_0.3size_19156tets_adaptive.inp ~/stochcasims/ 1
#
#
# OUTPUT 
#
# In (root)/data/HybridCaburst_stochBK/(mesh)/(iter_n+time) directory 
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

import steps.interface

import math
import time
from random import *
from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
import os

import extra.curr_funcs as cf
from extra.constants import *

import sys
import numpy as np

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

_, meshfile_ab, root, iter_n = sys.argv

if meshfile_ab == 'Cylinder2_dia2um_L160um_outer0_0.3shell_0.3size_279152tets_adaptive.inp':
    cyl160=True
else:
    cyl160=False

########################### BIOCHEMICAL MODEL ###############################

r = ReactionManager()

mdl_stoch = Model()
mdl_det = Model()

with mdl_stoch:
    Ca_stoch = Species.Create(valence=2)

    vsys_stoch = VolumeSystem.Create()
    ssys_stoch = SurfaceSystem.Create()

    with ssys_stoch:
        ##### BK channel ####################
        BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O0, BK_O1, BK_O2, BK_O3, BK_O4 = SubUnitState.Create()
        BKchan = Channel.Create([BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O0, BK_O1, BK_O2, BK_O3, BK_O4])

        with BKchan[...]:
            (((BK_C0.s + Ca_stoch.i <r['BKCAC0']> BK_C1.s)\
                       + Ca_stoch.i <r['BKCAC1']> BK_C2.s)\
                       + Ca_stoch.i <r['BKCAC2']> BK_C3.s)\
                       + Ca_stoch.i <r['BKCAC3']> BK_C4.s
            r['BKCAC0'].setRates(c_01, c_10)
            r['BKCAC1'].setRates(c_12, c_21)
            r['BKCAC2'].setRates(c_23, c_32)
            r['BKCAC3'].setRates(c_34, c_43)

            (((BK_O0.s + Ca_stoch.i <r['BKCAO0']> BK_O1.s)\
                       + Ca_stoch.i <r['BKCAO1']> BK_O2.s)\
                       + Ca_stoch.i <r['BKCAO2']> BK_O3.s)\
                       + Ca_stoch.i <r['BKCAO3']> BK_O4.s
            r['BKCAO0'].setRates(o_01, o_10)
            r['BKCAO1'].setRates(o_12, o_21)
            r['BKCAO2'].setRates(o_23, o_32)
            r['BKCAO3'].setRates(o_34, o_43)
            
            BK_C0.s <r['BKC0O0']> BK_O0.s
            BK_C1.s <r['BKC1O1']> BK_O1.s
            BK_C2.s <r['BKC2O2']> BK_O2.s
            BK_C3.s <r['BKC3O3']> BK_O3.s
            BK_C4.s <r['BKC4O4']> BK_O4.s
            r['BKC0O0'].setRates(VDepRate(f_0), VDepRate(b_0))
            r['BKC1O1'].setRates(VDepRate(f_1), VDepRate(b_1))
            r['BKC2O2'].setRates(VDepRate(f_2), VDepRate(b_2))
            r['BKC3O3'].setRates(VDepRate(f_3), VDepRate(b_3))
            r['BKC4O4'].setRates(VDepRate(f_4), VDepRate(b_4))
        

        openStates = BK_O0 | BK_O1 | BK_O2 | BK_O3 | BK_O4
        OC_BK = OhmicCurr.Create(BKchan[openStates], BK_G, BK_rev)
        # OC_BK0 = OhmicCurr.Create(BKchan[BK_O0], BK_G, BK_rev)
        # OC_BK1 = OhmicCurr.Create(BKchan[BK_O1], BK_G, BK_rev)
        # OC_BK2 = OhmicCurr.Create(BKchan[BK_O2], BK_G, BK_rev)
        # OC_BK3 = OhmicCurr.Create(BKchan[BK_O3], BK_G, BK_rev)
        # OC_BK4 = OhmicCurr.Create(BKchan[BK_O4], BK_G, BK_rev)


with mdl_det:
    Ca_det = Species.Create(valence=2)
    iCBsf, iCBsCa, iCBCaf, iCBCaCa, CBsf, CBsCa, CBCaf, CBCaCa, PV, PVMg, PVCa, Mg, Pump, CaPump,\
        CaP_m0, CaP_m1, CaP_m2, CaP_m3, CaT_m0h0, CaT_m0h1, CaT_m1h0, CaT_m1h1, CaT_m2h0, \
        CaT_m2h1, SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2 = Species.Create()
    
    # Vol/surface systems
    vsys_det = VolumeSystem.Create()
    ssys_det = SurfaceSystem.Create()

    with vsys_det:
        diff_Ca     = Diffusion.Create(Ca_det, DCST)
        diff_CBsf   = Diffusion.Create(CBsf, DCB)
        diff_CBsCa  = Diffusion.Create(CBsCa, DCB)
        diff_CBCaf  = Diffusion.Create(CBCaf, DCB)
        diff_CBCaCa = Diffusion.Create(CBCaCa, DCB)
        diff_PV     = Diffusion.Create(PV, DPV)
        diff_PVCa   = Diffusion.Create(PVCa, DPV)
        diff_PVMg   = Diffusion.Create(PVMg, DPV)

        (Ca_det + iCBsf <r['iCBsf1_f']> iCBsCa) + Ca_det <r['iCBsCa_f']> iCBCaCa
        (Ca_det + iCBsf <r['iCBsf2_f']> iCBCaf) + Ca_det <r['iCBCaf_f']> iCBCaCa
        r['iCBsf1_f'].setRates(iCBsf1_f_kcst, iCBsf1_b_kcst)
        r['iCBsCa_f'].setRates(iCBsCa_f_kcst, iCBsCa_b_kcst)
        r['iCBsf2_f'].setRates(iCBsf2_f_kcst, iCBsf2_b_kcst)
        r['iCBCaf_f'].setRates(iCBCaf_f_kcst, iCBCaf_b_kcst)

        (CBsf + Ca_det <r['CBsf1_f']> CBsCa) + Ca_det <r['CBsCa_f']> CBCaCa
        (CBsf + Ca_det <r['CBsf2_f']> CBCaf) + Ca_det <r['CBCaf_f']> CBCaCa
        r['CBsf1_f'].setRates(CBsf1_f_kcst, CBsf1_b_kcst)
        r['CBsCa_f'].setRates(CBsCa_f_kcst, CBsCa_b_kcst)
        r['CBsf2_f'].setRates(CBsf2_f_kcst, CBsf2_b_kcst)
        r['CBCaf_f'].setRates(CBCaf_f_kcst, CBCaf_b_kcst)

        Ca_det + PV <r['PVca_f']> PVCa
        Mg + PV <r['PVmg_f']> PVMg
        r['PVca_f'].setRates(PVca_f_kcst, PVca_b_kcst)
        r['PVmg_f'].setRates(PVmg_f_kcst, PVmg_b_kcst)

    with ssys_det:
        #Pump
        Ca_det.i + Pump.s <r['PumpD_f']> CaPump.s >r['PumpD_k']> Pump.s
        r['PumpD_f'].setRates(P_f_kcst, P_b_kcst)
        r['PumpD_k'].setRates(P_k_kcst)

        
        CaP_m0.s <r['CaPm0m1']> CaP_m1.s <r['CaPm1m2']> CaP_m2.s <r['CaPm2m3']> CaP_m3.s
        r['CaPm0m1'].setRates(0.0, 0.0)
        r['CaPm1m2'].setRates(0.0, 0.0)
        r['CaPm2m3'].setRates(0.0, 0.0)
    
        CaT_m0h0.s <r['CaTm0h0_m1h0']> CaT_m1h0.s <r['CaTm1h0_m2h0']> CaT_m2h0.s <r['CaTm2h0_m2h1']> CaT_m2h1.s
        r['CaTm0h0_m1h0'].setRates(0.0, 0.0)
        r['CaTm1h0_m2h0'].setRates(0.0, 0.0)
        r['CaTm2h0_m2h1'].setRates(0.0, 0.0)

        CaT_m1h0.s <r['CaTm1h0_m1h1']> CaT_m1h1.s
        r['CaTm1h0_m1h1'].setRates(0.0, 0.0)
        
        CaT_m0h0.s <r['CaTm0h0_m0h1']> CaT_m0h1.s <r['CaTm0h1_m1h1']> CaT_m1h1.s <r['CaTm1h1_m2h1']> CaT_m2h1.s
        r['CaTm0h0_m0h1'].setRates(0.0, 0.0)
        r['CaTm1h1_m2h1'].setRates(0.0, 0.0)
        r['CaTm0h1_m1h1'].setRates(0.0, 0.0)
    

        
        ((SK_C1.s + Ca_det.i <r['SKCAC1']> SK_C2.s)\
                  + Ca_det.i <r['SKCAC2']> SK_C3.s)\
                  + Ca_det.i <r['SKCAC3']> SK_C4.s
        r['SKCAC1'].setRates(dirc2_t, invc1_t)
        r['SKCAC2'].setRates(dirc3_t, invc2_t)
        r['SKCAC3'].setRates(dirc4_t, invc3_t)
        
        SK_C3.s <r['SKC3O1']> SK_O1.s
        SK_C4.s <r['SKC4O2']> SK_O2.s
        r['SKC3O1'].setRates(diro1_t, invo1_t)
        r['SKC4O2'].setRates(diro2_t, invo2_t)

##################################

########### MESH & COMPARTMENTALIZATION #################

##########Import Mesh

mesh_stoch = TetMesh.Load('./meshes/'+meshfile_ab)
mesh_det = TetMesh.Load('./meshes/'+meshfile_ab)

with mesh_stoch as mesh:
    # Use mesh_stoch for geometrical operations

    rad, zmin, zmax = 1e-6, -200e-6, 200e-6
    inner_tets, outer_tets = TetList(), TetList()
    for t in mesh.tets:
        c = t.center
        if zmin <= c.z <= zmax and c.x**2 + c.y**2 <= rad**2:
            inner_tets.append(t)
        else:
            outer_tets.append(t)

    print(len(outer_tets), " tets in outer compartment")
    print(len(inner_tets), " tets in inner compartment")

    # Record voltage from the central tetrahedron
    cent_tet = mesh.tets[0.0, 0.0, 0.0]

    if cyl160:
        # Ensure that we use points a small distance inside the boundary:
        minz, maxz = mesh.bbox.min.z, mesh.bbox.max.z
        memb_tris = TriList(tri for tri in mesh_stock.surface if minz < tri.center.z < maxz)
    else:
        print('Finding connecting triangles...')
        memb_tris = inner_tets.surface & outer_tets.surface

    submemb_tets = TetList()
    for tri in memb_tris:
        submemb_tets += tri.tetNeighbs
    submemb_tets = submemb_tets & inner_tets

    print(len(submemb_tets))

    vol = sum(tet.Vol for tet in submemb_tets)
    print('Volume of submembrane region is', vol)

    submemb_tris = TriList()
    for tet in submemb_tets:
        for tri in tet.faces:
            if tri in memb_tris:
                submemb_tris.append(tri)
                break

    assert(len(submemb_tris) == len(submemb_tets))

    
    ########## Create an intracellular compartment i.e. cytosolic compartment
    cyto_stoch = TetComp.Create(inner_tets.indices, vsys_stoch)
    ########## Create a membrane as a surface mesh
    memb_stoch = TetPatch.Create(memb_tris, cyto_stoch, None, 'ssys_stoch')

    # For EField calculation
    print("Creating membrane..")
    membrane = Membrane.Create([memb_stoch])
    print("Membrane created.")
    print("Area: ", memb_stoch.Area)

with mesh_det:
    cyto_det = TetComp.Create(inner_tets.indices, vsys_det)
    memb_det = TetPatch.Create(memb_tris.indices, cyto_det, None, 'ssys_det')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

rng1 = RNG('mt19937', 512, 7)
rng2 = RNG('mt19937', 512, 7)

#Creating two solvers
sim_stoch = Simulation('Tetexact', mdl_stoch, mesh_stoch, rng1, calcMembPot=True)

sim_det = Simulation('TetODE', mdl_det, mesh_det, rng2)
sim_det.setTolerances(1.0e-3, 1.0e-3)

print("Resetting simulation objects..")
sim_stoch.newRun()
sim_det.newRun()

print("Injecting molecules..")

sim_stoch.Temp = TEMPERATURE+273.15

sim_det.cyto_det.Ca_det.Conc = Ca_iconc
sim_stoch.cyto_stoch.Ca_stoch.Conc = Ca_iconc

print("Calcium concentration in stochastic simulation is: ", sim_stoch.cyto_stoch.Ca_stoch.Conc)
print("No. of Ca molecules in stochastic simulation is: ", sim_stoch.cyto_stoch.Ca_stoch.Count)

print("Calcium concentration in deterministic simulation is: ", sim_det.cyto_det.Ca_det.Conc)
print("No. of Ca molecules in deterministic simulation is: ", sim_det.cyto_det.Ca_det.Count)

sim_det.cyto_det.Mg.Conc = Mg_conc

surfarea = sim_stoch.memb_stoch.Area

pumpnbs = 6.022141e12*surfarea

sim_det.memb_det.Pump.Count = round(pumpnbs)
sim_det.memb_det.CaPump.Count = 0

sim_det.cyto_det.iCBsf.Conc = iCBsf_conc
sim_det.cyto_det.iCBCaf.Conc = iCBCaf_conc
sim_det.cyto_det.iCBsCa.Conc = iCBsCa_conc
sim_det.cyto_det.iCBCaCa.Conc = iCBCaCa_conc

sim_det.cyto_det.CBsf.Conc = CBsf_conc
sim_det.cyto_det.CBCaf.Conc = CBCaf_conc
sim_det.cyto_det.CBsCa.Conc = CBsCa_conc
sim_det.cyto_det.CBCaCa.Conc = CBCaCa_conc

sim_det.cyto_det.PV.Conc = PV_conc
sim_det.cyto_det.PVCa.Conc = PVCa_conc
sim_det.cyto_det.PVMg.Conc = PVMg_conc

# CaP

sim_det.memb_det.CaP_m0.Count = round(CaP_ro*surfarea*CaP_m0_p)
sim_det.memb_det.CaP_m1.Count = round(CaP_ro*surfarea*CaP_m1_p)
sim_det.memb_det.CaP_m2.Count = round(CaP_ro*surfarea*CaP_m2_p)
sim_det.memb_det.CaP_m3.Count = round(CaP_ro*surfarea*CaP_m3_p)


print("CaP_m0 ", round(CaP_ro*surfarea*CaP_m0_p))
print("CaP_m1 ", round(CaP_ro*surfarea*CaP_m1_p))
print("CaP_m2 ", round(CaP_ro*surfarea*CaP_m2_p))
print("CaP_m3 ", round(CaP_ro*surfarea*CaP_m3_p))

print("Targeted Injection: ", round(CaP_ro*surfarea), "CaP channels")

# CaT

# From cstate: CaT_m2h0 conducting
sim_det.memb_det.CaT_m0h0.Count = round(CaT_ro*surfarea*CaT_m0h0_p)
sim_det.memb_det.CaT_m1h0.Count = round(CaT_ro*surfarea*CaT_m1h0_p)
sim_det.memb_det.CaT_m2h0.Count = round(CaT_ro*surfarea*CaT_m2h0_p)
sim_det.memb_det.CaT_m0h1.Count = round(CaT_ro*surfarea*CaT_m0h1_p)
sim_det.memb_det.CaT_m1h1.Count = round(CaT_ro*surfarea*CaT_m1h1_p)
sim_det.memb_det.CaT_m2h1.Count = round(CaT_ro*surfarea*CaT_m2h1_p)

print("Injected  ", CaT_ro*surfarea, "CaT channels")

# BK
sim_stoch.memb_stoch.BKchan[BK_C0].Count = round(BK_ro*surfarea*BK_C0_p)
sim_stoch.memb_stoch.BKchan[BK_C1].Count = round(BK_ro*surfarea*BK_C1_p)
sim_stoch.memb_stoch.BKchan[BK_C2].Count = round(BK_ro*surfarea*BK_C2_p)
sim_stoch.memb_stoch.BKchan[BK_C3].Count = round(BK_ro*surfarea*BK_C3_p)
sim_stoch.memb_stoch.BKchan[BK_C4].Count = round(BK_ro*surfarea*BK_C4_p)

sim_stoch.memb_stoch.BKchan[BK_O0].Count = round(BK_ro*surfarea*BK_O0_p)
sim_stoch.memb_stoch.BKchan[BK_O1].Count = round(BK_ro*surfarea*BK_O1_p)
sim_stoch.memb_stoch.BKchan[BK_O2].Count = round(BK_ro*surfarea*BK_O2_p)
sim_stoch.memb_stoch.BKchan[BK_O3].Count = round(BK_ro*surfarea*BK_O3_p)
sim_stoch.memb_stoch.BKchan[BK_O4].Count = round(BK_ro*surfarea*BK_O4_p)


print("Injected  ", BK_ro*surfarea, "BK channels")

# SK
sim_det.memb_det.SK_C1.Count = round(SK_ro*surfarea*SK_C1_p)
sim_det.memb_det.SK_C2.Count = round(SK_ro*surfarea*SK_C2_p)
sim_det.memb_det.SK_C3.Count = round(SK_ro*surfarea*SK_C3_p)
sim_det.memb_det.SK_C4.Count = round(SK_ro*surfarea*SK_C4_p)

sim_det.memb_det.SK_O1.Count = round(SK_ro*surfarea*SK_O1_p)
sim_det.memb_det.SK_O2.Count = round(SK_ro*surfarea*SK_O2_p)


print("Injected ", SK_ro*surfarea, "SK channels")

sim_stoch.EfieldDT = EF_DT
sim_stoch.membrane.Potential = init_pot
sim_stoch.membrane.VolRes = Ra
#cm = 1.5uF/cm2 -> 1.5e-6F/1e-4m2 ->1.5e-2 F/m2
sim_stoch.membrane.Capac = memb_capac


#### Recording #####

dc = time.strftime('%b%d_%H_%M_%S_%Y')

runPath = os.path.join(root, 'data/HybridCaburst_stochBK/', meshfile_ab, f'{iter_n}__{dc}')
os.makedirs(runPath, exist_ok=True)

rng1.initialize(int(time.time()%1000))

# TMP NOT TRANSLATED YET
datfile  = open(os.path.join(runPath, 'currents.dat'), 'w')
datfile2 = open(os.path.join(runPath, 'voltage.dat'), 'w')
datfile3 = open(os.path.join(runPath, 'calcium.dat'), 'w')
# END TMP

stets = submemb_tets.indices
stris = submemb_tris.indices

for l in range(NTIMEPOINTS):
    print("Tpnt: ", l)

    #1) RUN STOCHASTIC SIMULATION i.e. compute currents and update stochastic calcium concentration
    sim_stoch.run(TIMECONVERTER*l)

    #2) READ STOCHASTIC CA and 3) SET DETERMINISTIC CA AND RATE CONSTANTS FOR DETERMINISTIC CHANNELS
    sim_det.TETS(stets).Ca_det.Conc = sim_stoch.TETS(stets).Ca_stoch.Conc

    #Assuming this sim V is not constant everwhere
    allPots = sim_stoch.TRIS(stris).V

    #3) Set the rate constants and RUN THE DETERMINISTIC SIMULATION
    sim_det.TRIS(stris).CaPm0m1['fwd'].K = [1.0e3 *3.* alpha_cap(V*1.0e3)*Qt for V in allPots]
    sim_det.TRIS(stris).CaPm1m2['fwd'].K = [1.0e3 *2.* alpha_cap(V*1.0e3)*Qt for V in allPots]
    sim_det.TRIS(stris).CaPm2m3['fwd'].K = [1.0e3 *1.* alpha_cap(V*1.0e3)*Qt for V in allPots]
    
    sim_det.TRIS(stris).CaPm2m3['bkw'].K = [1.0e3 *3.* beta_cap(V*1.0e3)*Qt for V in allPots]
    sim_det.TRIS(stris).CaPm1m2['bkw'].K = [1.0e3 *2.* beta_cap(V*1.0e3)*Qt for V in allPots]
    sim_det.TRIS(stris).CaPm0m1['bkw'].K = [1.0e3 *1.* beta_cap(V*1.0e3)*Qt for V in allPots]
    
    sim_det.TRIS(stris).CaTm0h0_m1h0['fwd'].K = [1.0e3 *2.* alpham_cat(V*1.0e3) for V in allPots]
    sim_det.TRIS(stris).CaTm1h0_m2h0['fwd'].K = [1.0e3 *1.* alpham_cat(V*1.0e3) for V in allPots]
    
    sim_det.TRIS(stris).CaTm1h0_m2h0['bkw'].K = [1.0e3 *2.* betam_cat(V*1.0e3) for V in allPots]
    sim_det.TRIS(stris).CaTm0h0_m1h0['bkw'].K = [1.0e3 *1.* betam_cat(V*1.0e3) for V in allPots]
    
    sim_det.TRIS(stris).CaTm0h0_m0h1['fwd'].K = [1.0e3 *1.* alphah_cat(V*1.0e3) for V in allPots]
    sim_det.TRIS(stris).CaTm1h0_m1h1['fwd'].K = [1.0e3 *1.* alphah_cat(V*1.0e3) for V in allPots]
    sim_det.TRIS(stris).CaTm2h0_m2h1['fwd'].K = [1.0e3 *1.* alphah_cat(V*1.0e3) for V in allPots]
    
    sim_det.TRIS(stris).CaTm2h0_m2h1['bkw'].K = [1.0e3 *1.* betah_cat(V*1.0e3) for V in allPots]
    sim_det.TRIS(stris).CaTm1h0_m1h1['bkw'].K = [1.0e3 *1.* betah_cat(V*1.0e3) for V in allPots]
    sim_det.TRIS(stris).CaTm0h0_m0h1['bkw'].K = [1.0e3 *1.* betah_cat(V*1.0e3) for V in allPots]

    sim_det.TRIS(stris).CaTm0h1_m1h1['fwd'].K = [1.0e3 *2.* alpham_cat(V*1.0e3) for V in allPots]
    sim_det.TRIS(stris).CaTm1h1_m2h1['fwd'].K = [1.0e3 *1.* alpham_cat(V*1.0e3) for V in allPots]

    sim_det.TRIS(stris).CaTm1h1_m2h1['bkw'].K = [1.0e3 *2.* betam_cat(V*1.0e3) for V in allPots]
    sim_det.TRIS(stris).CaTm0h1_m1h1['bkw'].K = [1.0e3 *1.* betam_cat(V*1.0e3) for V in allPots]

    #4) RUN DETERMINISTIC SIMULATION
    sim_det.run(TIMECONVERTER*l)

    # Now do the communication between the sims

    #5)READ DETERMINISTIC CHANNELS & THEN COMPUTE CURRENT USING DETERMINISTIC GHK (could be stochastic)
    So = Ca_oconc
    # i) For each tet in submembrane, find the corresponding triID
    # ii) For each tri, compute GHK current for each channel
    # iii) Count the channel states / Spec in open states for each of the triID and compute the total current of that channel 

    allCa = sim_det.TETS(stets).Ca_det.Conc

    currs_CaP = np.array([
        nb * cf.getGHKI(CaP_P, V, 2, TEMPERATURE+273.15, Si*1.0e3, So*1.0e3) 
        for V, Si, nb in zip(allPots, allCa, sim_det.TRIS(stris).CaP_m3.Count)
    ])

    currs_CaT = np.array([
        nb * cf.getGHKI(CaT_P, V, 2, TEMPERATURE+273.15, Si*1.0e3, So*1.0e3) 
        for V, Si, nb in zip(allPots, allCa, sim_det.TRIS(stris).CaT_m2h1.Count)
    ])

    currs_BK = np.array(sim_stoch.TRIS(stris).OC_BK.I)

    allSK = np.array(sim_det.TRIS(stris).LIST(SK_O1, SK_O2).Count).reshape(len(stris), 2)
    currs_SK = np.array([
        nb * cf.getOhmI(V, SK_rev, SK_G) 
        for V, nb in zip(allPots, np.sum(allSK, axis=1))
    ])

    membArea = sim_det.memb_det.Area
    currs_L = np.array([
        cf.getOhmI(V, L_rev, L_G) * round(L_ro * membArea) * (area / membArea)
        for V, area in zip(allPots, sim_stoch.TRIS(stris).Area)
    ])

    tcur_CaP = sum(currs_CaP)
    tcur_CaT = sum(currs_CaT)
    tcur_BK = sum(currs_BK)
    tcur_SK = sum(currs_SK)
    tca_count_det = sum(sim_det.TETS(stets).Ca_det.Count)
    tca_count_stoch = sum(sim_stoch.TETS(stets).Ca_stoch.Count)

    # Update sim stoch
    sim_stoch.TETS(stets).Ca_stoch.Count = sim_det.TETS(stets).Ca_det.Count - (currs_CaP + currs_CaT) * TIMECONVERTER / (2 * E_CHARGE)
    sim_stoch.TRIS(stris).IClamp = currs_CaP + currs_CaT + currs_SK + currs_L



##### TRANSLATION PROGRESSION TOKEN

    datfile.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile.write('%.6g' %((tcur_CaP*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_CaT*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_BK*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_SK*1.0e-1)/surfarea) + ' ') 
    datfile.write('\n')

    datfile2.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile2.write('%.6g' %(sim_stoch.TET(cent_tet).V*1.0e3) + ' ')
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

## END
