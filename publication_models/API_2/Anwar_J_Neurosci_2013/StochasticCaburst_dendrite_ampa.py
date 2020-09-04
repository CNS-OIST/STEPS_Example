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
# *StochasticCaburst_dendrite_ampa.py : The stochastic calcium burst 
# model, on the realistic dendrite morphology. Each dendritic branch is 
# modeled with well-mixed cytosolic compartments. Burst is synaptically-evoked
# with an ampa current.
#
# Script authors: Haroon Anwar and Iain Hepburn
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# USAGE
#
# $ python StochasticCaburst_dendrite_ampa.py *root* *iter_n* 
#
#  *root* is the path to the location for data storage 
#  *iter_n* (is intened to be an integer) is an identifier number for each 
#     simulation iteration
#
# E.g: python StochasticCaburst_dendrite_ampa.py ~/stochcasims/ 1
#
#
# OUTPUT 
#
# In (root)/data/StochasticCaburst_dendrite_ampa/chop_mesh.inp/(iter_n+time) directory 
# 3 data files will be recorded. Each file contains one row for every 
# time-point at which data is recorded, organised into the following columns:
# 
# currents.dat
# Time (ms), (for every compartment): P-type current, T-type current, BK current, SK current,
# AMPA current 
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

import steps.interface

import math
import time
from random import *
from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
import os
import numpy as np

import extra.curr_funcs as cf
from extra.constants_withampa import *

import sys

######Glutamate transient#######
# Reference (Rudolph et al. 2011)
#Units (mM)
with open("./extra/Glut_Pulse_MVR.dat","r") as f:
    Glut = [0.0]*25001
    count = 0
    for i in f.read().split():
        Glut[count] = float(i)
        count = count+1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

_, root, iter_n = sys.argv

########################### BIOCHEMICAL MODEL ###############################

# Two models required: Stochastic and deterministic
mdl_stoch = Model()
mdl_WM = Model()

r = ReactionManager()

with mdl_WM:
    
    # Calcium
    Ca = Species.Create(valence=2)

    Pump, CaPump, iCBsf, iCBsCa, iCBCaf, iCBCaCa, CBsf, CBsCa, CBCaf, CBCaCa, PV, PVMg, PVCa, Mg,\
        CaP_m0, CaP_m1, CaP_m2, CaP_m3, CaT_m0h0, CaT_m0h1, CaT_m1h0, CaT_m1h1, CaT_m2h0, CaT_m2h1,\
        BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O0, BK_O1, BK_O2, BK_O3, BK_O4, SK_C1, SK_C2, SK_C3, \
        SK_C4, SK_O1, SK_O2, AMPA_C, AMPA_C1, AMPA_C2, AMPA_D1, AMPA_D2, AMPA_O = Species.Create()
    
    
    ssys_diff, ssys_chans = SurfaceSystem.Create()
    vsys_buffs = VolumeSystem.Create()

    with ssys_diff:
        # Diffusions 
        Ca.o <r['diff_Ca_inward']> Ca.i
        CBsf.o <r['diff_CBsf_inward']> CBsf.i
        CBsCa.o <r['diff_CBsCa_inward']> CBsCa.i
        CBCaf.o <r['diff_CBCaf_inward']> CBCaf.i
        CBCaCa.o <r['diff_CBCaCa_inward']> CBCaCa.i
        PV.o <r['diff_PV_inward']> PV.i
        PVCa.o <r['diff_PVCa_inward']> PVCa.i
        PVMg.o <r['diff_PVMg_inward']> PVMg.i

        r['diff_Ca_inward'].setRates(0.0, 0.0)
        r['diff_CBsf_inward'].setRates(0.0, 0.0)
        r['diff_CBsCa_inward'].setRates(0.0, 0.0)
        r['diff_CBCaf_inward'].setRates(0.0, 0.0)
        r['diff_CBCaCa_inward'].setRates(0.0, 0.0)
        r['diff_PV_inward'].setRates(0.0, 0.0)
        r['diff_PVCa_inward'].setRates(0.0, 0.0)
        r['diff_PVMg_inward'].setRates(0.0, 0.0)

    with ssys_chans:
        #Pump
        Ca.i + Pump.s <r['PumpD_f']> CaPump.s
        r['PumpD_f'].setRates(P_f_kcst, P_b_kcst)
        
        CaPump.s >r['PumpD_k']> Pump.s
        r['PumpD_k'].setRates(P_k_kcst)

        ###### CaP channel ##############
    
        CaP_m0.s <r['CaPm0m1']> CaP_m1.s <r['CaPm1m2']> CaP_m2.s <r['CaPm2m3']> CaP_m3.s
        r['CaPm0m1'].setRates(0.0, 0.0)
        r['CaPm1m2'].setRates(0.0, 0.0)
        r['CaPm2m3'].setRates(0.0, 0.0)

        ######## CaT channel ##########

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

        ##### BK channel ####################

        (((BK_C0.s + Ca.i <r['BKCAC0']> BK_C1.s)\
                   + Ca.i <r['BKCAC1']> BK_C2.s)\
                   + Ca.i <r['BKCAC2']> BK_C3.s)\
                   + Ca.i <r['BKCAC3']> BK_C4.s
        r['BKCAC0'].setRates(c_01, c_10)
        r['BKCAC1'].setRates(c_12, c_21)
        r['BKCAC2'].setRates(c_23, c_32)
        r['BKCAC3'].setRates(c_34, c_43)

        (((BK_O0.s + Ca.i <r['BKCAO0']> BK_O1.s)\
                   + Ca.i <r['BKCAO1']> BK_O2.s)\
                   + Ca.i <r['BKCAO2']> BK_O3.s)\
                   + Ca.i <r['BKCAO3']> BK_O4.s
        r['BKCAO0'].setRates(o_01, o_10)
        r['BKCAO1'].setRates(o_12, o_21)
        r['BKCAO2'].setRates(o_23, o_32)
        r['BKCAO3'].setRates(o_34, o_43)
        
        BK_C0.s <r['BKC0O0']> BK_O0.s
        BK_C1.s <r['BKC1O1']> BK_O1.s
        BK_C2.s <r['BKC2O2']> BK_O2.s
        BK_C3.s <r['BKC3O3']> BK_O3.s
        BK_C4.s <r['BKC4O4']> BK_O4.s
        r['BKC0O0'].setRates(0, 0)
        r['BKC1O1'].setRates(0, 0)
        r['BKC2O2'].setRates(0, 0)
        r['BKC3O3'].setRates(0, 0)
        r['BKC4O4'].setRates(0, 0)

        ###### SK channel ################## DETERMINISTIC

        ((SK_C1.s + Ca.i <r['SKCAC1']> SK_C2.s)\
                  + Ca.i <r['SKCAC2']> SK_C3.s)\
                  + Ca.i <r['SKCAC3']> SK_C4.s
        r['SKCAC1'].setRates(dirc2_t, invc1_t)
        r['SKCAC2'].setRates(dirc3_t, invc2_t)
        r['SKCAC3'].setRates(dirc4_t, invc3_t)
        
        SK_C3.s <r['SKC3O1']> SK_O1.s
        SK_C4.s <r['SKC4O2']> SK_O2.s
        r['SKC3O1'].setRates(diro1_t, invo1_t)
        r['SKC4O2'].setRates(diro2_t, invo2_t)

        ###### AMPA channel ###########

        AMPA_C.s <r['AMPACC1']> AMPA_C1.s <r['AMPAC1C2']> AMPA_C2.s <r['AMPAC2O']> AMPA_O.s
        r['AMPACC1'].setRates(0.0, ru1)
        r['AMPAC1C2'].setRates(0.0, ru2)
        r['AMPAC2O'].setRates(ro, rc)
        
        AMPA_C1.s <r['AMPAD1']> AMPA_D1.s
        AMPA_C2.s <r['AMPAD2']> AMPA_D2.s
        r['AMPAD1'].setRates(rd, rr)
        r['AMPAD2'].setRates(rd, rr)

    with vsys_buffs:
        
        #iCBsf-fast
        Ca + iCBsf <r['iCBsf1_f']> iCBsCa
        r['iCBsf1_f'].setRates(iCBsf1_f_kcst, iCBsf1_b_kcst)
        
        #iCBsCa
        Ca + iCBsCa <r['iCBsCa_f']> iCBCaCa
        r['iCBsCa_f'].setRates(iCBsCa_f_kcst, iCBsCa_b_kcst)
        
        #iCBsf_slow
        Ca + iCBsf <r['iCBsf2_f']> iCBCaf
        r['iCBsf2_f'].setRates(iCBsf2_f_kcst, iCBsf2_b_kcst)
        
        #iCBCaf
        Ca + iCBCaf <r['iCBCaf_f']> iCBCaCa
        r['iCBCaf_f'].setRates(iCBCaf_f_kcst, iCBCaf_b_kcst)
        
        #CBsf-fast
        CBsf + Ca <r['CBsf1_f']> CBsCa
        r['CBsf1_f'].setRates(CBsf1_f_kcst, CBsf1_b_kcst)
        
        #CBsCa
        CBsCa + Ca <r['CBsCa_f']> CBCaCa
        r['CBsCa_f'].setRates(CBsCa_f_kcst, CBsCa_b_kcst)
        
        #CBsf_slow
        CBsf + Ca <r['CBsf2_f']> CBCaf
        r['CBsf2_f'].setRates(CBsf2_f_kcst, CBsf2_b_kcst)
        
        #CBCaf
        CBCaf + Ca <r['CBCaf_f']> CBCaCa
        r['CBCaf_f'].setRates(CBCaf_f_kcst, CBCaf_b_kcst)
        
        #PVca
        Ca + PV <r['PVca_f']> PVCa
        r['PVca_f'].setRates(PVca_f_kcst, PVca_b_kcst)
        
        #PVmg
        Mg + PV <r['PVmg_f']> PVMg
        r['PVmg_f'].setRates(PVmg_f_kcst, PVmg_b_kcst)

with mdl_stoch:
    #for EField
    ssys_stoch = SurfaceSystem.Create()

    ###### Leak current channel #####
    Leak = SubUnitState.Create()
    L = Channel.Create([Leak])

    with ssys_stoch:
        OC_L = OhmicCurr.Create(L[Leak], L_G, L_rev)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

########### MESH & COMPARTMENTALIZATION #################

##########Import Mesh for EField calculation

# For stochastic sim:

meshfile_ab = 'chop_mesh.inp'
mesh_stoch = TetMesh.LoadAbaqus('./meshes/'+meshfile_ab, scale=1e-06)

with mesh_stoch:
    ########## Create an intracellular compartment i.e. cytosolic compartment
    cyto_stoch = TetComp.Create(mesh_stoch.tets)
    print(f'{len(cyto_stoch.tets)} tets in inner compartment')

    #Define geometrical constants for all (100) compartments with concentric shells
    with open('./extra/geom_info_STEPS.dat','r') as f:
        data = np.array(list(map(float, f.read().split())))

        CUBIT_IDs, Length, Diam = data.reshape(len(data) // 3, 3).T
        Length *= 1e-6
        Diam *= 1e-6

    Comp_tetIDs = [mesh_stoch.tetGroups[f'EB{int(i)}'] for i in CUBIT_IDs]

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

    ########## Finding the triangles comprising the memberane
    Compborder_triIDs = []
    Compborder_tri_areas = []
    for i in range(len(Length)):
        Compborder_triIDs.append(TriList())
        Compborder_tri_areas.append(0)
        for tet in Comp_tetIDs[i]:
            for tri in tet.faces:
                if tri in mesh_stoch.surface:
                    Compborder_triIDs[-1].append(tri)
                    Compborder_tri_areas[-1] += tri.Area
                    break

    ########## Create a membrane as a surface mesh

    memb_stoch = TetPatch.Create(mesh_stoch.surface, cyto_stoch, None, ssys_stoch)
    
    # For EField calculation
    print("Creating membrane..")
    membrane = Membrane.Create([memb_stoch], opt_file_name='./meshes/'+meshfile_ab+"_optimalidx")
    print("Membrane created.")


#Geometry container object:

wmgeom = Geometry()
with wmgeom:
    shells = [None]*len(Length)
    rings = [None]*len(Length)
    for i in range(len(Length)):
        shells[i]=[None]*int(Nannulis[i])
        rings[i]=[None]*int(Nannulis[i])
        for j in range(int(Nannulis[i])):
            shells[i][j] = Compartment(vsys_buffs, Shells_vols[i][j], name=f'shells{i}{j}')
            if j==0:
                rings[i][j] = Patch(shells[i][j], None, ssys_chans, Rings_areas[i][j], name=f'rings{i}{j}')
            else:
                rings[i][j] = Patch(shells[i][j], shells[i][j-1], ssys_diff, Rings_areas[i][j], name=f'rings{i}{j}')




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

rng = RNG('mt19937', 512, int(time.time()%1000))

# Random-number generator 
r_dummy = RNG('mt19937', 512, int(time.time()%1000))

print("Creating tetexact solver...")
sim_stoch = Simulation('Tetexact', mdl_stoch, mesh_stoch, rng, calcMembPot=True)

print("Creating WM solver")
sim_WM = Simulation('Wmdirect', mdl_WM, wmgeom, r_dummy)

dc = time.strftime('%b%d_%H_%M_%S_%Y')
runPath = os.path.join(root, 'data/StochasticCaburst_dendrite_ampa/', meshfile_ab, f'{iter_n}__{dc}')
os.makedirs(runPath, exist_ok=True)

sim_stoch.autoCheckpoint(0.1, os.path.join(runPath, 'checkpoint__sim_stoch__'))
sim_WM.autoCheckpoint(0.1, os.path.join(runPath, 'checkpoint__sim_WM__'))




print("Resetting simulation objects..")
sim_stoch.newRun()
sim_WM.newRun()

    
print("Injecting molecules..")

sim_stoch.Temp = TEMPERATURE+273.15

allShells = sim_WM.MATCH('shells\d+')

allShells.Ca.Conc = Ca_iconc
allShells.Mg.Conc = Mg_conc

allShells.iCBsf.Conc = iCBsf_conc
allShells.iCBsCa.Conc = iCBsCa_conc
allShells.iCBCaf.Conc = iCBCaf_conc
allShells.iCBCaCa.Conc = iCBCaCa_conc

allShells.CBsf.Conc = CBsf_conc
allShells.CBsCa.Conc = CBsCa_conc
allShells.CBCaf.Conc = CBCaf_conc
allShells.CBCaCa.Conc = CBCaCa_conc

allShells.PV.Conc = PV_conc
allShells.PVCa.Conc = PVCa_conc
allShells.PVMg.Conc = PVMg_conc

for i in range(len(Length)):

    ## For j == 0 ##
    outerRing = sim_WM.LIST(rings[i][0])
    # TODO replace long root by this root

    surfarea = rings[i][0].Area
    pumpnbs = 6.022141e12 * surfarea

    outerRing.Pump.Count = pumpnbs
    outerRing.CaPump.Count = 0

    # CaP
    outerRing.CaP_m0.Count = round(CaP_ro*surfarea*CaP_m0_p)
    outerRing.CaP_m1.Count = round(CaP_ro*surfarea*CaP_m1_p)
    outerRing.CaP_m2.Count = round(CaP_ro*surfarea*CaP_m2_p)
    outerRing.CaP_m3.Count = round(CaP_ro*surfarea*CaP_m3_p)

    print("Injected  ", CaP_ro*surfarea, "CaP channels")

    # CaT only in spiny dendrites
    if i>1:
        # From cstate: CaT_m2h0 conducting
        outerRing.CaT_m0h0.Count = round(CaT_ro*surfarea*CaT_m0h0_p)
        outerRing.CaT_m1h0.Count = round(CaT_ro*surfarea*CaT_m1h0_p)
        outerRing.CaT_m2h0.Count = round(CaT_ro*surfarea*CaT_m2h0_p)
        outerRing.CaT_m0h1.Count = round(CaT_ro*surfarea*CaT_m0h1_p)
        outerRing.CaT_m1h1.Count = round(CaT_ro*surfarea*CaT_m1h1_p)
        outerRing.CaT_m2h1.Count = round(CaT_ro*surfarea*CaT_m2h1_p)
        print("Injected  ", CaT_ro*surfarea, "CaT channels")
    else:
        outerRing.AMPA_C.Count = round(AMPA_receptors)
        outerRing.AMPA_C1.Count = 0
        outerRing.AMPA_C2.Count = 0
        outerRing.AMPA_O.Count = 0
        outerRing.AMPA_D1.Count = 0
        outerRing.AMPA_D2.Count = 0
        print("Injected  ", AMPA_receptors, "AMPA receptors")

    # BK
    outerRing.BK_C0.Count = round(BK_ro*surfarea*BK_C0_p)
    outerRing.BK_C1.Count = round(BK_ro*surfarea*BK_C1_p)
    outerRing.BK_C2.Count = round(BK_ro*surfarea*BK_C2_p)
    outerRing.BK_C3.Count = round(BK_ro*surfarea*BK_C3_p)
    outerRing.BK_C4.Count = round(BK_ro*surfarea*BK_C4_p)

    outerRing.BK_O0.Count = round(BK_ro*surfarea*BK_O0_p)
    outerRing.BK_O1.Count = round(BK_ro*surfarea*BK_O1_p)
    outerRing.BK_O2.Count = round(BK_ro*surfarea*BK_O2_p)
    outerRing.BK_O3.Count = round(BK_ro*surfarea*BK_O3_p)
    outerRing.BK_O4.Count = round(BK_ro*surfarea*BK_O4_p)
    print("Injected  ", BK_ro*surfarea, "BK channels")

    # SK
    outerRing.SK_C1.Count = round(SK_ro*surfarea*SK_C1_p)
    outerRing.SK_C2.Count = round(SK_ro*surfarea*SK_C2_p)
    outerRing.SK_C3.Count = round(SK_ro*surfarea*SK_C3_p)
    outerRing.SK_C4.Count = round(SK_ro*surfarea*SK_C4_p)

    outerRing.SK_O1.Count = round(SK_ro*surfarea*SK_O1_p)
    outerRing.SK_O2.Count = round(SK_ro*surfarea*SK_O2_p)
    
    print("Injected ", SK_ro*surfarea, "SK channels")

    for j in range(1, int(Nannulis[i])):
        #set the rate constants for diffusion (Diffusion is modeled as surface reaction here)
        ring = sim_WM.LIST(rings[i][j])

        #for Ca diffusion
        ring.diff_Ca_inward['fwd'].K = (DCST*DScales_ShellIn[i][j-1])
        ring.diff_Ca_inward['bkw'].K = (DCST*DScales_ShellOut[i][j-1])

        #for CBsf diffusin
        ring.diff_CBsf_inward['fwd'].K = (DCB*DScales_ShellIn[i][j-1])
        ring.diff_CBsf_inward['bkw'].K = (DCB*DScales_ShellOut[i][j-1])

        #for CBsCa diffusion
        ring.diff_CBsCa_inward['fwd'].K = (DCB*DScales_ShellIn[i][j-1])
        ring.diff_CBsCa_inward['bkw'].K = (DCB*DScales_ShellOut[i][j-1])

        #for CBCaf diffusion
        ring.diff_CBCaf_inward['fwd'].K = (DCB*DScales_ShellIn[i][j-1])
        ring.diff_CBCaf_inward['bkw'].K = (DCB*DScales_ShellOut[i][j-1])

        #for CBCaCa diffusion
        ring.diff_CBCaCa_inward['fwd'].K = (DCB*DScales_ShellIn[i][j-1])
        ring.diff_CBCaCa_inward['bkw'].K = (DCB*DScales_ShellOut[i][j-1])

        #for PV diffusion
        ring.diff_PV_inward['fwd'].K = (DPV*DScales_ShellIn[i][j-1])
        ring.diff_PV_inward['bkw'].K = (DPV*DScales_ShellOut[i][j-1])

        #for PVCa diffusion
        ring.diff_PVCa_inward['fwd'].K = (DPV*DScales_ShellIn[i][j-1])
        ring.diff_PVCa_inward['bkw'].K = (DPV*DScales_ShellOut[i][j-1])

        #for PVMg diffusion
        ring.diff_PVMg_inward['fwd'].K = (DPV*DScales_ShellIn[i][j-1])
        ring.diff_PVMg_inward['bkw'].K = (DPV*DScales_ShellOut[i][j-1])

surfarea = sim_stoch.memb_stoch.Area
#Compute the surface area for full mesh rather than individual compartments
sim_stoch.memb_stoch.L[Leak].Count = round(L_ro * surfarea)
print("Injected  ", (L_ro * sim_stoch.memb_stoch.Area), "Leak channels")

sim_stoch.EfieldDT = EF_DT

sim_stoch.membrane.Potential = init_pot

sim_stoch.membrane.VolRes = Ra

#cm = 1.5uF/cm2 -> 1.5e-6F/1e-4m2 ->1.5e-2 F/m2
sim_stoch.membrane.Capac = memb_capac

#### Recording #####


datfile =  open(os.path.join(runPath, 'currents.dat'), 'w')
datfile2 = open(os.path.join(runPath, 'voltage.dat'), 'w')
datfile3 = open(os.path.join(runPath, 'calcium.dat'), 'w')

for l in range(NTIMEPOINTS):
    print("Tpnt: ", l)

    datfile.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile2.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile3.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')

    # IN this sim V should be constant everywhere in the compartment
    for i in range(len(Length)):
        V = sim_stoch.TRI(Compborder_triIDs[i][0]).V

        outerRing = sim_WM.LIST(rings[i][0])
    
        outerRing.CaPm0m1['fwd'].K = 1.0e3 *3.* alpha_cap(V*1.0e3)*Qt
        outerRing.CaPm1m2['fwd'].K = 1.0e3 *2.* alpha_cap(V*1.0e3)*Qt
        outerRing.CaPm2m3['fwd'].K = 1.0e3 *1.* alpha_cap(V*1.0e3)*Qt
        outerRing.CaPm2m3['bkw'].K = 1.0e3 *3.* beta_cap(V*1.0e3)*Qt
        outerRing.CaPm1m2['bkw'].K = 1.0e3 *2.* beta_cap(V*1.0e3)*Qt
        outerRing.CaPm0m1['bkw'].K = 1.0e3 *1.* beta_cap(V*1.0e3)*Qt

        if i>1:    
            outerRing.CaTm0h0_m1h0['fwd'].K = 1.0e3 *2.* alpham_cat(V*1.0e3)
            outerRing.CaTm1h0_m2h0['fwd'].K = 1.0e3 *1.* alpham_cat(V*1.0e3)
    
            outerRing.CaTm1h0_m2h0['bkw'].K = 1.0e3 *2.* betam_cat(V*1.0e3)
            outerRing.CaTm0h0_m1h0['bkw'].K = 1.0e3 *1.* betam_cat(V*1.0e3)
    
            outerRing.CaTm0h0_m0h1['fwd'].K = 1.0e3 *1.* alphah_cat(V*1.0e3)
            outerRing.CaTm1h0_m1h1['fwd'].K = 1.0e3 *1.* alphah_cat(V*1.0e3)
            outerRing.CaTm2h0_m2h1['fwd'].K = 1.0e3 *1.* alphah_cat(V*1.0e3)
    
            outerRing.CaTm2h0_m2h1['bkw'].K = 1.0e3 *1.* betah_cat(V*1.0e3)
            outerRing.CaTm1h0_m1h1['bkw'].K = 1.0e3 *1.* betah_cat(V*1.0e3)
            outerRing.CaTm0h0_m0h1['bkw'].K = 1.0e3 *1.* betah_cat(V*1.0e3)
    
            outerRing.CaTm0h1_m1h1['fwd'].K = 1.0e3 *2.* alpham_cat(V*1.0e3)
            outerRing.CaTm1h1_m2h1['fwd'].K = 1.0e3 *1.* alpham_cat(V*1.0e3)
    
            outerRing.CaTm1h1_m2h1['bkw'].K = 1.0e3 *2.* betam_cat(V*1.0e3)
            outerRing.CaTm0h1_m1h1['bkw'].K = 1.0e3 *1.* betam_cat(V*1.0e3)

        else:
            outerRing.AMPACC1['fwd'].K = 1.0e-3 *rb*Glut[l]
            outerRing.AMPAC1C2['fwd'].K = 1.0e-3 *rb*Glut[l]

        outerRing.BKC0O0['fwd'].K = f_0(V)
        outerRing.BKC1O1['fwd'].K = f_1(V)
        outerRing.BKC2O2['fwd'].K = f_2(V)
        outerRing.BKC3O3['fwd'].K = f_3(V)
        outerRing.BKC4O4['fwd'].K = f_4(V)
        outerRing.BKC0O0['bkw'].K = b_0(V)
        outerRing.BKC1O1['bkw'].K = b_1(V)
        outerRing.BKC2O2['bkw'].K = b_2(V)
        outerRing.BKC3O3['bkw'].K = b_3(V)
        outerRing.BKC4O4['bkw'].K = b_4(V)
    
    sim_WM.run(TIMECONVERTER*l)

    # Now do the communication between the sims
    for i in range(len(Length)):
        V = sim_stoch.TRI(Compborder_triIDs[i][0]).V

        outerRing = sim_WM.LIST(rings[i][0])
        outerShell = sim_WM.LIST(shells[i][0])

        Si = outerShell.Ca.Conc

        So = Ca_oconc
    
        # Get the single-channel currents first
        tcur_CaP_sc = cf.getGHKI(CaP_P, V, 2, TEMPERATURE+273.15, Si*1.0e3, So*1.0e3)
        tcur_BK_sc = cf.getOhmI(V, BK_rev, BK_G)
        tcur_SK_sc = cf.getOhmI(V, SK_rev, SK_G)

        if i>1:
            tcur_CaT_sc = cf.getGHKI(CaT_P, V, 2, TEMPERATURE+273.15, Si*1.0e3, So*1.0e3)
            tcur_CaT = tcur_CaT_sc * outerRing.CaT_m2h1.Count
            tcur_AMPA = 0.0
        else:
            tcur_CaT = 0.0
            tcur_AMPA_sc = cf.getOhmI(V, AMPA_rev, AMPA_G)
            tcur_AMPA = tcur_AMPA_sc * outerRing.AMPA_O.Count
    
        tcur_CaP = tcur_CaP_sc * outerRing.CaP_m3.Count
        # alpha is to h1
        tcur_BK = tcur_BK_sc * sum(outerRing.MATCH('^BK_O[0-4]$').Count)
        tcur_SK = tcur_SK_sc * sum(outerRing.MATCH('^SK_O[1-2]$').Count)

        AreaRatio = Compborder_tri_areas[i]/Rings_areas[i][0]
        nmtris = len(Compborder_triIDs[i])
    
        for tri in Compborder_triIDs[i]: 
            sim_stoch.TRI(tri).IClamp = (tcur_CaP+tcur_CaT+tcur_BK+tcur_SK+tcur_AMPA)*AreaRatio/nmtris

        ca_count_inj =  -1.0*((tcur_CaP+tcur_CaT)*TIMECONVERTER)/(2*E_CHARGE)

        t_count = outerShell.Ca.Count
        outerShell.Ca.Count = t_count + ca_count_inj

        datfile.write('%.6g' %((tcur_CaP*1.0e-1)/Compborder_tri_areas[i]) + ' ')
        datfile.write('%.6g' %((tcur_CaT*1.0e-1)/Compborder_tri_areas[i]) + ' ')
        datfile.write('%.6g' %((tcur_BK*1.0e-1)/Compborder_tri_areas[i]) + ' ')
        datfile.write('%.6g' %((tcur_SK*1.0e-1)/Compborder_tri_areas[i]) + ' ')
        datfile.write('%.6g' %((tcur_AMPA*1.0e-1)/Compborder_tri_areas[i]) + ' ')
        
        datfile3.write('%.6g' %outerShell.Ca.Count +' ')
        datfile3.write('%.6g' %(outerShell.Ca.Conc*1.0e6)+ ' ')

    datfile.write('\n')
    datfile3.write('\n')
    
    sim_stoch.run(TIMECONVERTER*l)
    
    for i in range(len(Length)):
        V = sim_stoch.TRI(Compborder_triIDs[i][0]).V
        datfile2.write('%.6g' %(V*1.0e3) + ' ')
    datfile2.write('\n')


datfile.close()
datfile2.close()
datfile3.close()
