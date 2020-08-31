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

from __future__ import print_function
import math
# WARNING: Using a variable name that is reserved (['time']).
import time
from random import *
from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
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
mdl_stoch = Model()
# WARNING: Using a variable name that is reserved (['r']).
r = ReactionManager()
mdl_WM = Model()
with mdl_WM:
    
    # Calcium
    Ca = Species.Create(valence=2)
    
    
    # Pump
    # CaPump
    
    # iCBsf
    # iCBsCa
    # iCBCaf
    # iCBCaCa
    
    # CBsf
    # CBsCa
    # CBCaf
    # CBCaCa
    
    # PV
    # PVMg
    # PVCa
    # Mg
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # CHANNELS  # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    
    ###### CaP channel ############## 
    
    
    
    ######## CaT channel ##########  
    
    
    ##### BK channel ####################
    
    
    
    ###### SK channel ################## DETERMINISTIC
    
    
    
    # Pump
    Pump, CaPump, iCBsf, iCBsCa, iCBCaf, iCBCaCa, CBsf, CBsCa, CBCaf, CBCaCa, PV, PVMg, PVCa, Mg, CaP_m0, CaP_m1, CaP_m2, CaP_m3, CaT_m0h0, CaT_m0h1, CaT_m1h0, CaT_m1h1, CaT_m2h0, CaT_m2h1, BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O0, BK_O1, BK_O2, BK_O3, BK_O4, SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2 = Species.Create()
    
    ssys_diff = SurfaceSystem.Create()
    
    ssys_chans = SurfaceSystem.Create()
    
    vsys_buffs = VolumeSystem.Create()
with mdl_stoch:
    
    #for EField
    ssys_stoch = SurfaceSystem.Create()
with ssys_diff, mdl_stoch:
    
    
    # Diffusions 
    Ca.o <r['diff_Ca_inward']> Ca.i ; r['diff_Ca_inward'].setRates(0.0, 0.0)
    
    CBsf.o <r['diff_CBsf_inward']> CBsf.i ; r['diff_CBsf_inward'].setRates(0.0, 0.0)
    
    CBsCa.o <r['diff_CBsCa_inward']> CBsCa.i ; r['diff_CBsCa_inward'].setRates(0.0, 0.0)
    
    CBCaf.o <r['diff_CBCaf_inward']> CBCaf.i ; r['diff_CBCaf_inward'].setRates(0.0, 0.0)
    
    CBCaCa.o <r['diff_CBCaCa_inward']> CBCaCa.i ; r['diff_CBCaCa_inward'].setRates(0.0, 0.0)
    
    PV.o <r['diff_PV_inward']> PV.i ; r['diff_PV_inward'].setRates(0.0, 0.0)
    
    PVCa.o <r['diff_PVCa_inward']> PVCa.i ; r['diff_PVCa_inward'].setRates(0.0, 0.0)
    
    PVMg.o <r['diff_PVMg_inward']> PVMg.i ; r['diff_PVMg_inward'].setRates(0.0, 0.0)
with ssys_chans, mdl_stoch:
    
    
    #Pump
    
    Ca.i + Pump.s <r['PumpD_f']> CaPump.s
    r['PumpD_f'].setRates(P_f_kcst, P_b_kcst)
    
    CaPump.s >r['PumpD_k']> Pump.s
r['PumpD_k'].setRates(P_k_kcst)
with vsys_buffs, mdl_stoch:
    
    #iCBsf-fast
    Ca + iCBsf <r['iCBsf1_f']> iCBsCa ; r['iCBsf1_f'].setRates(iCBsf1_f_kcst, iCBsf1_b_kcst)
    
    #iCBsCa
    Ca + iCBsCa <r['iCBsCa_f']> iCBCaCa ; r['iCBsCa_f'].setRates(iCBsCa_f_kcst, iCBsCa_b_kcst)
    
    #iCBsf_slow
    Ca + iCBsf <r['iCBsf2_f']> iCBCaf ; r['iCBsf2_f'].setRates(iCBsf2_f_kcst, iCBsf2_b_kcst)
    
    #iCBCaf
    Ca + iCBCaf <r['iCBCaf_f']> iCBCaCa ; r['iCBCaf_f'].setRates(iCBCaf_f_kcst, iCBCaf_b_kcst)
    
    #CBsf-fast
    CBsf + Ca <r['CBsf1_f']> CBsCa ; r['CBsf1_f'].setRates(CBsf1_f_kcst, CBsf1_b_kcst)
    
    #CBsCa
    CBsCa + Ca <r['CBsCa_f']> CBCaCa ; r['CBsCa_f'].setRates(CBsCa_f_kcst, CBsCa_b_kcst)
    
    #CBsf_slow
    CBsf + Ca <r['CBsf2_f']> CBCaf ; r['CBsf2_f'].setRates(CBsf2_f_kcst, CBsf2_b_kcst)
    
    #CBCaf
    CBCaf + Ca <r['CBCaf_f']> CBCaCa ; r['CBCaf_f'].setRates(CBCaf_f_kcst, CBCaf_b_kcst)
    
    #PVca
    Ca + PV <r['PVca_f']> PVCa ; r['PVca_f'].setRates(PVca_f_kcst, PVca_b_kcst)
    
    #PVmg
    Mg + PV <r['PVmg_f']> PVMg ; r['PVmg_f'].setRates(PVmg_f_kcst, PVmg_b_kcst)
with mdl_WM:
    pass
with ssys_chans, mdl_stoch:
    
    
    CaP_m0.s <r['CaPm0m1']> CaP_m1.s ; r['CaPm0m1'].setRates(0.0, 0.0)
    CaP_m1.s <r['CaPm1m2']> CaP_m2.s ; r['CaPm1m2'].setRates(0.0, 0.0)
    
    CaP_m2.s <r['CaPm2m3']> CaP_m3.s ; r['CaPm2m3'].setRates(0.0, 0.0)
with mdl_WM:
    pass
with ssys_chans, mdl_stoch:
    
    
    CaT_m0h0.s <r['CaTm0h0_m1h0']> CaT_m1h0.s ; r['CaTm0h0_m1h0'].setRates(0.0, 0.0)
    
    CaT_m1h0.s <r['CaTm1h0_m2h0']> CaT_m2h0.s ; r['CaTm1h0_m2h0'].setRates(0.0, 0.0)
    
    CaT_m0h1.s <r['CaTm0h1_m1h1']> CaT_m1h1.s ; r['CaTm0h1_m1h1'].setRates(0.0, 0.0)
    
    CaT_m1h1.s <r['CaTm1h1_m2h1']> CaT_m2h1.s ; r['CaTm1h1_m2h1'].setRates(0.0, 0.0)
    
    
    CaT_m0h0.s <r['CaTm0h0_m0h1']> CaT_m0h1.s ; r['CaTm0h0_m0h1'].setRates(0.0, 0.0)
    CaT_m1h0.s <r['CaTm1h0_m1h1']> CaT_m1h1.s ; r['CaTm1h0_m1h1'].setRates(0.0, 0.0)
    
    CaT_m2h0.s <r['CaTm2h0_m2h1']> CaT_m2h1.s ; r['CaTm2h0_m2h1'].setRates(0.0, 0.0)
with mdl_WM:
    pass
with ssys_chans, mdl_stoch:
    
    
    
    Ca.i + BK_C0.s <r['BKCAC0']> BK_C1.s ; r['BKCAC0'].setRates(c_01, c_10)
    Ca.i + BK_C1.s <r['BKCAC1']> BK_C2.s ; r['BKCAC1'].setRates(c_12, c_21)
    Ca.i + BK_C2.s <r['BKCAC2']> BK_C3.s ; r['BKCAC2'].setRates(c_23, c_32)
    Ca.i + BK_C3.s <r['BKCAC3']> BK_C4.s ; r['BKCAC3'].setRates(c_34, c_43)
    
    
    Ca.i + BK_O0.s <r['BKCAO0']> BK_O1.s ; r['BKCAO0'].setRates(o_01, o_10)
    Ca.i + BK_O1.s <r['BKCAO1']> BK_O2.s ; r['BKCAO1'].setRates(o_12, o_21)
    Ca.i + BK_O2.s <r['BKCAO2']> BK_O3.s ; r['BKCAO2'].setRates(o_23, o_32)
    Ca.i + BK_O3.s <r['BKCAO3']> BK_O4.s ; r['BKCAO3'].setRates(o_34, o_43)
    
    
    BK_C0.s <r['BKC0O0']> BK_O0.s ; r['BKC0O0'].setRates(0.0, 0.0)
    BK_C1.s <r['BKC1O1']> BK_O1.s ; r['BKC1O1'].setRates(0.0, 0.0)
    BK_C2.s <r['BKC2O2']> BK_O2.s ; r['BKC2O2'].setRates(0.0, 0.0)
    BK_C3.s <r['BKC3O3']> BK_O3.s ; r['BKC3O3'].setRates(0.0, 0.0)
    BK_C4.s <r['BKC4O4']> BK_O4.s ; r['BKC4O4'].setRates(0.0, 0.0)
with mdl_WM:
    pass
with ssys_chans, mdl_stoch:
    
    
    
    Ca.i + SK_C1.s <r['SKCAC1']> SK_C2.s ; r['SKCAC1'].setRates(dirc2_t, invc1_t)
    Ca.i + SK_C2.s <r['SKCAC2']> SK_C3.s ; r['SKCAC2'].setRates(dirc3_t, invc2_t)
    Ca.i + SK_C3.s <r['SKCAC3']> SK_C4.s ; r['SKCAC3'].setRates(dirc4_t, invc3_t)
    
    
    SK_C3.s <r['SKC3O1']> SK_O1.s ; r['SKC3O1'].setRates(diro1_t, invo1_t)
    SK_C4.s <r['SKC4O2']> SK_O2.s ; r['SKC4O2'].setRates(diro2_t, invo2_t)
with mdl_stoch:
    Leak = SubUnitState.Create()
    
    
    ###### Leak current channel #####
    
    L = Channel.Create([Leak])
with ssys_stoch:
    
    OC_L = OhmicCurr.Create(L[Leak], L_G, L_rev)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # /Users/anwar/Copy/ModelDB_scripts/extra/geom_info_STEPS.dat# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

########### MESH & COMPARTMENTALIZATION #################

##########Import Mesh for EField calculation

# For stochastic sim:

meshfile_ab = 'chop_mesh.inp'
trip = TetMesh.LoadAbaqus('./meshes/'+meshfile_ab, scale=1e-06)
tetgroups = tetp.blocksToGroups()

inner_tets = range(len(mesh_stoch.tets))

print(inner_tets.__len__(), " tets in inner compartment")


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
    print(len(Comp_tetIDs[i]))


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

wmgeom = Geometry()
with ssys_stoch:
    
    OC_L = OhmicCurr.Create(L[Leak], L_G, L_rev)

for i in range(len(Length)):
    shells[i]=[None]*int(Nannulis[i])
    rings[i]=[None]*int(Nannulis[i])
    for j in range(int(Nannulis[i])):
        shells[i][j] = Compartment.Create(None, Shells_vols[i][j])
        if j==0:
            rings[i][j] = Patch.Create(shells[i][j], None, None, Rings_areas[i][j])
        else:
            rings[i][j] = Patch.Create(shells[i][j], shells[i][j-1], None, Rings_areas[i][j])


#For all patches, add the surface system

for i in range(len(Length)):
    for j in range(int(Nannulis[i])):
        if j==0:
            rings[i][j].addSurfsys('ssys_chans')
        else:
            rings[i][j].addSurfsys('ssys_diff')

        shells[i][j].addVolsys('vsys_buffs')
with mesh_stoch:


########## Create an intracellular compartment i.e. cytosolic compartment

    cyto_stoch = TetComp.Create(inner_tets)
    
    
    ########## Finding the triangles comprising the memberane
    
    memb_tris = list(mesh_stoch.surface.indices)
    
    Compborder_triIDs = [None]*len(Length)
    Compborder_tri_areas = [0.0]*len(Length)
    
    for i in range(len(Length)):
        Compborder_triIDs[i]=[]
        for j in range(len(Comp_tetIDs[i])):
            tritemp = mesh_stoch.tets[Comp_tetIDs[i][j]].faces.indices
            for tri in tritemp:
                if tri in memb_tris:
                    Compborder_triIDs[i].append(tri)
                    Compborder_tri_areas[i] = Compborder_tri_areas[i] + mesh_stoch.tris[tri].Area
                    break

########## Create a membrane as a surface mesh

# Stochastic sim:
    memb_stoch = TetPatch.Create(memb_tris, cyto_stoch, None, 'ssys_stoch')
    
    # For EField calculation
    print("Creating membrane..")
    membrane = Membrane.Create([memb_stoch], opt_file_name='./meshes/'+meshfile_ab+"_optimalidx")
print("Membrane created.")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

# WARNING: Using a variable name that is reserved (['r']).
r = RNG('mt19937', 512, int(time.time()%1000))

# Random-number generator 
r_dummy = RNG('mt19937', 512, int(time.time()%1000))

print("Creating tetexact solver...")
# WARNING: Using a variable name that is reserved (['r']).
sim_stoch = Simulation('Tetexact', mdl_stoch, mesh_stoch, r, calcMembPot=True)

print("Creating WM solver")
sim_WM = Simulation('Wmdirect', mdl_WM, wmgeom, r_dummy)

print("Resetting simulation objects..")
sim_stoch.newRun()
sim_WM.newRun()

print("Injecting molecules..")

sim_stoch.Temp = TEMPERATURE+273.15

for i in range(len(Length)):
    for j in range(int(Nannulis[i])):
        getattr(sim_WM, 'shells'+str(i)+str(j)).Ca.Conc = Ca_iconc
        getattr(sim_WM, 'shells'+str(i)+str(j)).Mg.Conc = Mg_conc

        getattr(sim_WM, 'shells'+str(i)+str(j)).iCBsf.Conc = iCBsf_conc
        getattr(sim_WM, 'shells'+str(i)+str(j)).iCBsCa.Conc = iCBsCa_conc
        getattr(sim_WM, 'shells'+str(i)+str(j)).iCBCaf.Conc = iCBCaf_conc
        getattr(sim_WM, 'shells'+str(i)+str(j)).iCBCaCa.Conc = iCBCaCa_conc

        getattr(sim_WM, 'shells'+str(i)+str(j)).CBsf.Conc = CBsf_conc
        getattr(sim_WM, 'shells'+str(i)+str(j)).CBsCa.Conc = CBsCa_conc
        getattr(sim_WM, 'shells'+str(i)+str(j)).CBCaf.Conc = CBCaf_conc
        getattr(sim_WM, 'shells'+str(i)+str(j)).CBCaCa.Conc = CBCaCa_conc

        getattr(sim_WM, 'shells'+str(i)+str(j)).PV.Conc = PV_conc
        getattr(sim_WM, 'shells'+str(i)+str(j)).PVCa.Conc = PVCa_conc
        getattr(sim_WM, 'shells'+str(i)+str(j)).PVMg.Conc = PVMg_conc
    
        if j==0:
            surfarea = getattr(sim_WM, 'rings'+str(i)+str(j)).Area
            pumpnbs = 6.022141e12*surfarea

            getattr(sim_WM, 'rings'+str(i)+str(j)).Pump.Count = pumpnbs
            getattr(sim_WM, 'rings'+str(i)+str(j)).CaPump.Count = 0

            # CaP
            getattr(sim_WM, 'rings'+str(i)+str(j)).CaP_m0.Count = round(CaP_ro*surfarea*CaP_m0_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).CaP_m1.Count = round(CaP_ro*surfarea*CaP_m1_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).CaP_m2.Count = round(CaP_ro*surfarea*CaP_m2_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).CaP_m3.Count = round(CaP_ro*surfarea*CaP_m3_p)

            print("Injected  ", CaP_ro*surfarea, "CaP channels")

            # CaT

            # From cstate: CaT_m2h0 conducting
            getattr(sim_WM, 'rings'+str(i)+str(j)).CaT_m0h0.Count = round(CaT_ro*surfarea*CaT_m0h0_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).CaT_m1h0.Count = round(CaT_ro*surfarea*CaT_m1h0_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).CaT_m2h0.Count = round(CaT_ro*surfarea*CaT_m2h0_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).CaT_m0h1.Count = round(CaT_ro*surfarea*CaT_m0h1_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).CaT_m1h1.Count = round(CaT_ro*surfarea*CaT_m1h1_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).CaT_m2h1.Count = round(CaT_ro*surfarea*CaT_m2h1_p)

            print("Injected  ", CaT_ro*surfarea, "CaT channels")

            # BK
            getattr(sim_WM, 'rings'+str(i)+str(j)).BK_C0.Count = round(BK_ro*surfarea*BK_C0_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).BK_C1.Count = round(BK_ro*surfarea*BK_C1_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).BK_C2.Count = round(BK_ro*surfarea*BK_C2_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).BK_C3.Count = round(BK_ro*surfarea*BK_C3_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).BK_C4.Count = round(BK_ro*surfarea*BK_C4_p)

            getattr(sim_WM, 'rings'+str(i)+str(j)).BK_O0.Count = round(BK_ro*surfarea*BK_O0_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).BK_O1.Count = round(BK_ro*surfarea*BK_O1_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).BK_O2.Count = round(BK_ro*surfarea*BK_O2_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).BK_O3.Count = round(BK_ro*surfarea*BK_O3_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).BK_O4.Count = round(BK_ro*surfarea*BK_O4_p)


            print("Injected  ", BK_ro*surfarea, "BK channels")

            # SK
            getattr(sim_WM, 'rings'+str(i)+str(j)).SK_C1.Count = round(SK_ro*surfarea*SK_C1_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).SK_C2.Count = round(SK_ro*surfarea*SK_C2_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).SK_C3.Count = round(SK_ro*surfarea*SK_C3_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).SK_C4.Count = round(SK_ro*surfarea*SK_C4_p)

            getattr(sim_WM, 'rings'+str(i)+str(j)).SK_O1.Count = round(SK_ro*surfarea*SK_O1_p)
            getattr(sim_WM, 'rings'+str(i)+str(j)).SK_O2.Count = round(SK_ro*surfarea*SK_O2_p)
            
            print("Injected ", SK_ro*surfarea, "SK channels")
        
        else:
            #set the rate constants for diffusion (Diffusion is modeled as surface reaction here)

            #for Ca diffusion

            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_Ca_inward['fwd'].K = (DCST*DScales_ShellIn[i][j-1])
            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_Ca_inward['bkw'].K = (DCST*DScales_ShellOut[i][j-1])

            #for CBsf diffusin
            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_CBsf_inward['fwd'].K = (DCB*DScales_ShellIn[i][j-1])
            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_CBsf_inward['bkw'].K = (DCB*DScales_ShellOut[i][j-1])

            #for CBsCa diffusion

            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_CBsCa_inward['fwd'].K = (DCB*DScales_ShellIn[i][j-1])
            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_CBsCa_inward['bkw'].K = (DCB*DScales_ShellOut[i][j-1])

            #for CBCaf diffusion

            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_CBCaf_inward['fwd'].K = (DCB*DScales_ShellIn[i][j-1])
            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_CBCaf_inward['bkw'].K = (DCB*DScales_ShellOut[i][j-1])

            #for CBCaCa diffusion

            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_CBCaCa_inward['fwd'].K = (DCB*DScales_ShellIn[i][j-1])
            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_CBCaCa_inward['bkw'].K = (DCB*DScales_ShellOut[i][j-1])

            #for PV diffusion

            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_PV_inward['fwd'].K = (DPV*DScales_ShellIn[i][j-1])
            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_PV_inward['bkw'].K = (DPV*DScales_ShellOut[i][j-1])

            #for PVCa diffusion

            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_PVCa_inward['fwd'].K = (DPV*DScales_ShellIn[i][j-1])
            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_PVCa_inward['bkw'].K = (DPV*DScales_ShellOut[i][j-1])

            #for PVMg diffusion

            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_PVMg_inward['fwd'].K = (DPV*DScales_ShellIn[i][j-1])
            getattr(sim_WM, 'rings'+str(i)+str(j)).diff_PVMg_inward['bkw'].K = (DPV*DScales_ShellOut[i][j-1])


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

# WARNING: Using a variable name that is reserved (['time']).
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

# WARNING: Using a variable name that is reserved (['time', 'time']).
btime = time.time()
for l in range(NTIMEPOINTS):
    print("Tpnt: ", l)

    if TIMECONVERTER*l in cp_times:
        # WARNING: Using a variable name that is reserved (['checkpoint']).
        sim_stoch.checkpoint(root+'data/' +  'StochasticCaburst_dendrite/'+meshfile_ab+'/'+iter_n+'__'+dc + '/checkpoint__sim_stoch__'+  str(TIMECONVERTER*l))
        # WARNING: Using a variable name that is reserved (['checkpoint']).
        sim_WM.checkpoint(root+'data/' +  'StochasticCaburst_dendrite/'+meshfile_ab+'/'+iter_n+'__'+dc + '/checkpoint__sim_WM__'+  str(TIMECONVERTER*l))
        datfile.flush()
        datfile2.flush()
        datfile3.flush()
    
    datfile.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile2.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile3.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    # IN this sim V should be constant everywhere in the compartment
    for i in range(len(Length)):
        # WARNING: Using a variable name that is reserved (['V']).
        V = sim_stoch.TRI(Compborder_triIDs[i][0]).V
    
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaPm0m1['fwd'].K = 1.0e3 *3.* alpha_cap(V*1.0e3)*Qt
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaPm1m2['fwd'].K = 1.0e3 *2.* alpha_cap(V*1.0e3)*Qt
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaPm2m3['fwd'].K = 1.0e3 *1.* alpha_cap(V*1.0e3)*Qt
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaPm2m3['bkw'].K = 1.0e3 *3.* beta_cap(V*1.0e3)*Qt
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaPm1m2['bkw'].K = 1.0e3 *2.* beta_cap(V*1.0e3)*Qt
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaPm0m1['bkw'].K = 1.0e3 *1.* beta_cap(V*1.0e3)*Qt
    
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaTm0h0_m1h0['fwd'].K = 1.0e3 *2.* alpham_cat(V*1.0e3)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaTm1h0_m2h0['fwd'].K = 1.0e3 *1.* alpham_cat(V*1.0e3)
    
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaTm1h0_m2h0['bkw'].K = 1.0e3 *2.* betam_cat(V*1.0e3)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaTm0h0_m1h0['bkw'].K = 1.0e3 *1.* betam_cat(V*1.0e3)
    
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaTm0h0_m0h1['fwd'].K = 1.0e3 *1.* alphah_cat(V*1.0e3)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaTm1h0_m1h1['fwd'].K = 1.0e3 *1.* alphah_cat(V*1.0e3)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaTm2h0_m2h1['fwd'].K = 1.0e3 *1.* alphah_cat(V*1.0e3)
    
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaTm2h0_m2h1['bkw'].K = 1.0e3 *1.* betah_cat(V*1.0e3)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaTm1h0_m1h1['bkw'].K = 1.0e3 *1.* betah_cat(V*1.0e3)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaTm0h0_m0h1['bkw'].K = 1.0e3 *1.* betah_cat(V*1.0e3)
    
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaTm0h1_m1h1['fwd'].K = 1.0e3 *2.* alpham_cat(V*1.0e3)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaTm1h1_m2h1['fwd'].K = 1.0e3 *1.* alpham_cat(V*1.0e3)
    
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaTm1h1_m2h1['bkw'].K = 1.0e3 *2.* betam_cat(V*1.0e3)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).CaTm0h1_m1h1['bkw'].K = 1.0e3 *1.* betam_cat(V*1.0e3)
    
    
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).BKC0O0['fwd'].K = f_0(V)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).BKC1O1['fwd'].K = f_1(V)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).BKC2O2['fwd'].K = f_2(V)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).BKC3O3['fwd'].K = f_3(V)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).BKC4O4['fwd'].K = f_4(V)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).BKC0O0['bkw'].K = b_0(V)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).BKC1O1['bkw'].K = b_1(V)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).BKC2O2['bkw'].K = b_2(V)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).BKC3O3['bkw'].K = b_3(V)
        # WARNING: Using a variable name that is reserved (['V']).
        getattr(sim_WM, 'rings'+str(i)+str(0)).BKC4O4['bkw'].K = b_4(V)
    
    # WARNING: Using a variable name that is reserved (['run']).
    sim_WM.run(TIMECONVERTER*l)

    # Now do the communication between the sims
    for i in range(len(Length)):
        # WARNING: Using a variable name that is reserved (['V']).
        V = sim_stoch.TRI(Compborder_triIDs[i][0]).V
        Si = getattr(sim_WM, 'shells'+str(i)+str(0)).Ca.Conc

        So = Ca_oconc
    
        # Get the single-channel currents first
        # WARNING: Using a variable name that is reserved (['V']).
        tcur_CaP_sc = cf.getGHKI(CaP_P, V, 2, TEMPERATURE+273.15, Si*1.0e3, So*1.0e3)
        # WARNING: Using a variable name that is reserved (['V']).
        tcur_CaT_sc = cf.getGHKI(CaT_P, V, 2, TEMPERATURE+273.15, Si*1.0e3, So*1.0e3)
        # WARNING: Using a variable name that is reserved (['V']).
        tcur_BK_sc = cf.getOhmI(V, BK_rev, BK_G)
    
        # WARNING: Using a variable name that is reserved (['V']).
        tcur_SK_sc = cf.getOhmI(V, SK_rev, SK_G)
    
        tcur_CaP = tcur_CaP_sc*getattr(sim_WM, 'rings'+str(i)+str(0)).CaP_m3.Count
        # alpha is to h1
        tcur_CaT = tcur_CaT_sc*getattr(sim_WM, 'rings'+str(i)+str(0)).CaT_m2h1.Count
        tcur_BK = tcur_BK_sc*((getattr(sim_WM, 'rings'+str(i)+str(0)).BK_O0.Count+\
                                  getattr(sim_WM, 'rings'+str(i)+str(0)).BK_O1.Count+\
                                     getattr(sim_WM, 'rings'+str(i)+str(0)).BK_O2.Count+\
                                        getattr(sim_WM, 'rings'+str(i)+str(0)).BK_O3.Count+\
                                           getattr(sim_WM, 'rings'+str(i)+str(0)).BK_O4.Count))
        tcur_SK = tcur_SK_sc*((getattr(sim_WM, 'rings'+str(i)+str(0)).SK_O1.Count+\
                                  getattr(sim_WM, 'rings'+str(i)+str(0)).SK_O2.Count))

        AreaRatio = Compborder_tri_areas[i]/Rings_areas[i][0]
        nmtris = len(Compborder_triIDs[i])
    
        for tri in Compborder_triIDs[i]: sim_stoch.setTriIClamp(tri, (((tcur_CaP+tcur_CaT+tcur_BK+tcur_SK)*AreaRatio)/nmtris))

        ca_count_inj =  -1.0*((tcur_CaP+tcur_CaT)*TIMECONVERTER)/(2*E_CHARGE)

        t_count = getattr(sim_WM, 'shells'+str(i)+str(0)).Ca.Count
        getattr(sim_WM, 'shells'+str(i)+str(0)).Ca.Count = t_count + ca_count_inj

        datfile.write('%.6g' %((tcur_CaP*1.0e-1)/Compborder_tri_areas[i]) + ' ')
        datfile.write('%.6g' %((tcur_CaT*1.0e-1)/Compborder_tri_areas[i]) + ' ')
        datfile.write('%.6g' %((tcur_BK*1.0e-1)/Compborder_tri_areas[i]) + ' ')
        datfile.write('%.6g' %((tcur_SK*1.0e-1)/Compborder_tri_areas[i]) + ' ')

        datfile3.write('%.6g' %getattr(sim_WM, 'shells'+str(i)+str(0)).Ca.Count +' ')
        datfile3.write('%.6g' %(getattr(sim_WM, 'shells'+str(i)+str(0)).Ca.Conc*1.0e6)+ ' ')

    datfile.write('\n')
    datfile3.write('\n')
    
    # WARNING: Using a variable name that is reserved (['run']).
    sim_stoch.run(TIMECONVERTER*l)
    
    for i in range(len(Length)):
        # WARNING: Using a variable name that is reserved (['V']).
        V = sim_stoch.TRI(Compborder_triIDs[i][0]).V
        # WARNING: Using a variable name that is reserved (['V']).
        datfile2.write('%.6g' %(V*1.0e3) + ' ')
    datfile2.write('\n')


datfile.close()
datfile2.close()
datfile3.close()

