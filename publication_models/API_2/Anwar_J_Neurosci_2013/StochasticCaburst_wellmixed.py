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

import steps.interface

import math
import time
from random import *
from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
import os

from extra.constants import *
import extra.curr_funcs as cf

import sys

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

_, meshfile_ab, root, iter_n = sys.argv

if meshfile_ab == 'Cylinder2_dia2um_L160um_outer0_0.3shell_0.3size_279152tets_adaptive.inp':
    cyl160=True
else:
    cyl160=False

########################### BIOCHEMICAL MODEL ###############################

# Two models required: Stochastic and deterministic
mdl_stoch = Model()
mdl_WM = Model()

r = ReactionManager()

with mdl_WM:
    
    # Calcium
    Ca = Species.Create(valence=2)
    
    # Species
    Pump, CaPump, iCBsf, iCBsCa, iCBCaf, iCBCaCa, CBsf, CBsCa, CBCaf, CBCaCa, PV, PVMg, PVCa, Mg,\
        CaP_m0, CaP_m1, CaP_m2, CaP_m3, CaT_m0h0, CaT_m0h1, CaT_m1h0, CaT_m1h1, CaT_m2h0, CaT_m2h1,\
        BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O0, BK_O1, BK_O2, BK_O3, BK_O4, SK_C1, SK_C2, SK_C3, \
        SK_C4, SK_O1, SK_O2 = Species.Create()
    
    # Vol/surface systems
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

        r['diff_Ca_inward'].K = 0.0, 0.0
        r['diff_CBsf_inward'].K = 0.0, 0.0
        r['diff_CBsCa_inward'].K = 0.0, 0.0
        r['diff_CBCaf_inward'].K = 0.0, 0.0
        r['diff_CBCaCa_inward'].K = 0.0, 0.0
        r['diff_PV_inward'].K = 0.0, 0.0
        r['diff_PVCa_inward'].K = 0.0, 0.0
        r['diff_PVMg_inward'].K = 0.0, 0.0

    with ssys_chans:
        #Pump
        Ca.i + Pump.s <r['PumpD_f']> CaPump.s
        r['PumpD_f'].K = P_f_kcst, P_b_kcst
        
        CaPump.s >r['PumpD_k']> Pump.s
        r['PumpD_k'].K = P_k_kcst

        ###### CaP channel ##############
    
        CaP_m0.s <r['CaPm0m1']> CaP_m1.s <r['CaPm1m2']> CaP_m2.s <r['CaPm2m3']> CaP_m3.s
        r['CaPm0m1'].K = 0.0, 0.0
        r['CaPm1m2'].K = 0.0, 0.0
        r['CaPm2m3'].K = 0.0, 0.0

        ######## CaT channel ##########

        CaT_m0h0.s <r['CaTm0h0_m1h0']> CaT_m1h0.s <r['CaTm1h0_m2h0']> CaT_m2h0.s <r['CaTm2h0_m2h1']> CaT_m2h1.s
        r['CaTm0h0_m1h0'].K = 0.0, 0.0
        r['CaTm1h0_m2h0'].K = 0.0, 0.0
        r['CaTm2h0_m2h1'].K = 0.0, 0.0

        CaT_m1h0.s <r['CaTm1h0_m1h1']> CaT_m1h1.s
        r['CaTm1h0_m1h1'].K = 0.0, 0.0
        
        CaT_m0h0.s <r['CaTm0h0_m0h1']> CaT_m0h1.s <r['CaTm0h1_m1h1']> CaT_m1h1.s <r['CaTm1h1_m2h1']> CaT_m2h1.s
        r['CaTm0h0_m0h1'].K = 0.0, 0.0
        r['CaTm1h1_m2h1'].K = 0.0, 0.0
        r['CaTm0h1_m1h1'].K = 0.0, 0.0

        ##### BK channel ####################

        (((BK_C0.s + Ca.i <r['BKCAC0']> BK_C1.s)\
                   + Ca.i <r['BKCAC1']> BK_C2.s)\
                   + Ca.i <r['BKCAC2']> BK_C3.s)\
                   + Ca.i <r['BKCAC3']> BK_C4.s
        r['BKCAC0'].K = c_01, c_10
        r['BKCAC1'].K = c_12, c_21
        r['BKCAC2'].K = c_23, c_32
        r['BKCAC3'].K = c_34, c_43

        (((BK_O0.s + Ca.i <r['BKCAO0']> BK_O1.s)\
                   + Ca.i <r['BKCAO1']> BK_O2.s)\
                   + Ca.i <r['BKCAO2']> BK_O3.s)\
                   + Ca.i <r['BKCAO3']> BK_O4.s
        r['BKCAO0'].K = o_01, o_10
        r['BKCAO1'].K = o_12, o_21
        r['BKCAO2'].K = o_23, o_32
        r['BKCAO3'].K = o_34, o_43
        
        BK_C0.s <r['BKC0O0']> BK_O0.s
        BK_C1.s <r['BKC1O1']> BK_O1.s
        BK_C2.s <r['BKC2O2']> BK_O2.s
        BK_C3.s <r['BKC3O3']> BK_O3.s
        BK_C4.s <r['BKC4O4']> BK_O4.s
        r['BKC0O0'].K = 0, 0
        r['BKC1O1'].K = 0, 0
        r['BKC2O2'].K = 0, 0
        r['BKC3O3'].K = 0, 0
        r['BKC4O4'].K = 0, 0

        ###### SK channel ################## DETERMINISTIC

        ((SK_C1.s + Ca.i <r['SKCAC1']> SK_C2.s)\
                  + Ca.i <r['SKCAC2']> SK_C3.s)\
                  + Ca.i <r['SKCAC3']> SK_C4.s
        r['SKCAC1'].K = dirc2_t, invc1_t
        r['SKCAC2'].K = dirc3_t, invc2_t
        r['SKCAC3'].K = dirc4_t, invc3_t
        
        SK_C3.s <r['SKC3O1']> SK_O1.s
        SK_C4.s <r['SKC4O2']> SK_O2.s
        r['SKC3O1'].K = diro1_t, invo1_t
        r['SKC4O2'].K = diro2_t, invo2_t

    with vsys_buffs:
        
        #iCBsf-fast
        Ca + iCBsf <r['iCBsf1_f']> iCBsCa
        r['iCBsf1_f'].K = iCBsf1_f_kcst, iCBsf1_b_kcst
        
        #iCBsCa
        Ca + iCBsCa <r['iCBsCa_f']> iCBCaCa
        r['iCBsCa_f'].K = iCBsCa_f_kcst, iCBsCa_b_kcst
        
        #iCBsf_slow
        Ca + iCBsf <r['iCBsf2_f']> iCBCaf
        r['iCBsf2_f'].K = iCBsf2_f_kcst, iCBsf2_b_kcst
        
        #iCBCaf
        Ca + iCBCaf <r['iCBCaf_f']> iCBCaCa
        r['iCBCaf_f'].K = iCBCaf_f_kcst, iCBCaf_b_kcst
        
        #CBsf-fast
        CBsf + Ca <r['CBsf1_f']> CBsCa
        r['CBsf1_f'].K = CBsf1_f_kcst, CBsf1_b_kcst
        
        #CBsCa
        CBsCa + Ca <r['CBsCa_f']> CBCaCa
        r['CBsCa_f'].K = CBsCa_f_kcst, CBsCa_b_kcst
        
        #CBsf_slow
        CBsf + Ca <r['CBsf2_f']> CBCaf
        r['CBsf2_f'].K = CBsf2_f_kcst, CBsf2_b_kcst
        
        #CBCaf
        CBCaf + Ca <r['CBCaf_f']> CBCaCa
        r['CBCaf_f'].K = CBCaf_f_kcst, CBCaf_b_kcst
        
        #PVca
        Ca + PV <r['PVca_f']> PVCa
        r['PVca_f'].K = PVca_f_kcst, PVca_b_kcst
        
        #PVmg
        Mg + PV <r['PVmg_f']> PVMg
        r['PVmg_f'].K = PVmg_f_kcst, PVmg_b_kcst

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

##########Import Mesh

mesh_stoch = TetMesh.Load('./meshes/'+meshfile_ab)

with mesh_stoch:
    rad, zmin, zmax = 1e-6, -200e-6, 200e-6
    inner_tets, outer_tets = TetList(), TetList()
    for t in mesh_stoch.tets:
        c = t.center
        if zmin <= c.z <= zmax and c.x**2 + c.y**2 <= rad**2:
            inner_tets.append(t)
        else:
            outer_tets.append(t)

    print(len(outer_tets), " tets in outer compartment")
    print(len(inner_tets), " tets in inner compartment")

    # Record voltage from the central tetrahedron
    cent_tet = mesh_stoch.tets[0.0, 0.0, 0.0]

    ########## Create an intracellular compartment i.e. cytosolic compartment
    cyto_stoch = Compartment.Create(inner_tets)

    if cyl160:
        # Ensure that we use points a small distance inside the boundary:
        minz, maxz = mesh_stoch.bbox.min.z, mesh_stoch.bbox.max.z
        memb_tris = TriList(tri for tri in mesh_stock.surface if minz < tri.center.z < maxz)
    else:
        print('Finding connecting triangles...')
        memb_tris = inner_tets.surface & outer_tets.surface

    ########## Create a membrane as a surface mesh
    memb_stoch = Patch.Create(memb_tris, cyto_stoch, None, ssys_stoch)

    # For EField calculation
    print("Creating membrane..")
    membrane = Membrane.Create([memb_stoch])
    print("Membrane created.")

#Define geometrical constants for a compartment with concentric shells
Length = 80e-6

shell0_rmax = 1.0e-6
shell1_rmax = 0.9e-6
shell2_rmax = 0.7e-6
shell3_rmax = 0.5e-6
shell4_rmax = 0.3e-6
shell5_rmax = 0.1e-6

#Geometry container object:

wmgeom = Geometry()
with wmgeom:
    
    #Set up compartments
    
    shell0 = Compartment.Create('vsys_buffs', math.pi*((shell0_rmax**2) - (shell1_rmax**2))*Length)
    shell1 = Compartment.Create('vsys_buffs', math.pi*((shell1_rmax**2) - (shell2_rmax**2))*Length)
    shell2 = Compartment.Create('vsys_buffs', math.pi*((shell2_rmax**2) - (shell3_rmax**2))*Length)
    shell3 = Compartment.Create('vsys_buffs', math.pi*((shell3_rmax**2) - (shell4_rmax**2))*Length)
    shell4 = Compartment.Create('vsys_buffs', math.pi*((shell4_rmax**2) - (shell5_rmax**2))*Length)
    shell5 = Compartment.Create('vsys_buffs', math.pi*((shell5_rmax**2))*Length)
    
    #Set up patches: lower number compartment is 'outer', higher number is 'inner'
    ring0 = Patch.Create(shell0, None, 'ssys_chans', 2*math.pi*shell0_rmax*Length)
    
    ring1 = Patch.Create(shell1, shell0, 'ssys_diff', 2*math.pi*shell1_rmax*Length)
    ring2 = Patch.Create(shell2, shell1, 'ssys_diff', 2*math.pi*shell2_rmax*Length)
    ring3 = Patch.Create(shell3, shell2, 'ssys_diff', 2*math.pi*shell3_rmax*Length)
    ring4 = Patch.Create(shell4, shell3, 'ssys_diff', 2*math.pi*shell4_rmax*Length)
    ring5 = Patch.Create(shell5, shell4, 'ssys_diff', 2*math.pi*shell5_rmax*Length)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

rng = RNG('mt19937', 512, seed=int(time.time()%1000))

# Random-number generator 
r_dummy = RNG('mt19937', 512, seed=int(time.time()%1000))

print("Creating tetexact solver...")
sim_stoch = Simulation('Tetexact', mdl_stoch, mesh_stoch, rng, calcMembPot=True)

print("Creating WM solver")
sim_WM = Simulation('Wmdirect', mdl_WM, wmgeom, r_dummy)

#### Recording #####

dc = time.strftime('%b%d_%H_%M_%S_%Y')

runPath = os.path.join(root, 'data/StochasticCaburst_wellmixed/', meshfile_ab, f'{iter_n}__{dc}')
os.makedirs(runPath, exist_ok=True)

datfile =  open(os.path.join(runPath, 'currents.dat'), 'w')
datfile2 = open(os.path.join(runPath, 'voltage.dat'), 'w')
datfile3 = open(os.path.join(runPath, 'calcium.dat'), 'w')

print("Resetting simulation objects..")
sim_stoch.newRun()
sim_WM.newRun()

print("Injecting molecules..")

sim_stoch.Temp = TEMPERATURE + 273.15

sim_WM.shell0.Ca.Conc = Ca_iconc
sim_WM.shell1.Ca.Conc = Ca_iconc
sim_WM.shell2.Ca.Conc = Ca_iconc
sim_WM.shell3.Ca.Conc = Ca_iconc
sim_WM.shell4.Ca.Conc = Ca_iconc
sim_WM.shell5.Ca.Conc = Ca_iconc

sim_WM.shell0.Mg.Conc = Mg_conc
sim_WM.shell1.Mg.Conc = Mg_conc
sim_WM.shell2.Mg.Conc = Mg_conc
sim_WM.shell3.Mg.Conc = Mg_conc
sim_WM.shell4.Mg.Conc = Mg_conc
sim_WM.shell5.Mg.Conc = Mg_conc

surfarea = sim_WM.ring0.Area


pumpnbs = 6.022141e12*surfarea

sim_WM.ring0.Pump.Count = pumpnbs
sim_WM.ring0.CaPump.Count = 0

print("Injected ", sim_WM.ring0.Pump.Count, "pumps")

sim_WM.shell0.iCBsf.Conc = iCBsf_conc
sim_WM.shell0.iCBsCa.Conc = iCBsCa_conc
sim_WM.shell0.iCBCaf.Conc = iCBCaf_conc
sim_WM.shell0.iCBCaCa.Conc = iCBCaCa_conc

sim_WM.shell0.CBsf.Conc = CBsf_conc
sim_WM.shell0.CBsCa.Conc = CBsCa_conc
sim_WM.shell0.CBCaf.Conc = CBCaf_conc
sim_WM.shell0.CBCaCa.Conc = CBCaCa_conc

sim_WM.shell0.PV.Conc = PV_conc
sim_WM.shell0.PVCa.Conc = PVCa_conc
sim_WM.shell0.PVMg.Conc = PVMg_conc

sim_WM.shell1.iCBsf.Conc = iCBsf_conc
sim_WM.shell1.iCBsCa.Conc = iCBsCa_conc
sim_WM.shell1.iCBCaf.Conc = iCBCaf_conc
sim_WM.shell1.iCBCaCa.Conc = iCBCaCa_conc

sim_WM.shell1.CBsf.Conc = CBsf_conc
sim_WM.shell1.CBsCa.Conc = CBsCa_conc
sim_WM.shell1.CBCaf.Conc = CBCaf_conc
sim_WM.shell1.CBCaCa.Conc = CBCaCa_conc

sim_WM.shell1.PV.Conc = PV_conc
sim_WM.shell1.PVCa.Conc = PVCa_conc
sim_WM.shell1.PVMg.Conc = PVMg_conc

sim_WM.shell2.iCBsf.Conc = iCBsf_conc
sim_WM.shell2.iCBsCa.Conc = iCBsCa_conc
sim_WM.shell2.iCBCaf.Conc = iCBCaf_conc
sim_WM.shell2.iCBCaCa.Conc = iCBCaCa_conc

sim_WM.shell2.CBsf.Conc = CBsf_conc
sim_WM.shell2.CBsCa.Conc = CBsCa_conc
sim_WM.shell2.CBCaf.Conc = CBCaf_conc
sim_WM.shell2.CBCaCa.Conc = CBCaCa_conc

sim_WM.shell2.PV.Conc = PV_conc
sim_WM.shell2.PVCa.Conc = PVCa_conc
sim_WM.shell2.PVMg.Conc = PVMg_conc

sim_WM.shell3.iCBsf.Conc = iCBsf_conc
sim_WM.shell3.iCBsCa.Conc = iCBsCa_conc
sim_WM.shell3.iCBCaf.Conc = iCBCaf_conc
sim_WM.shell3.iCBCaCa.Conc = iCBCaCa_conc

sim_WM.shell3.CBsf.Conc = CBsf_conc
sim_WM.shell3.CBsCa.Conc = CBsCa_conc
sim_WM.shell3.CBCaf.Conc = CBCaf_conc
sim_WM.shell3.CBCaCa.Conc = CBCaCa_conc

sim_WM.shell3.PV.Conc = PV_conc
sim_WM.shell3.PVCa.Conc = PVCa_conc
sim_WM.shell3.PVMg.Conc = PVMg_conc

sim_WM.shell4.iCBsf.Conc = iCBsf_conc
sim_WM.shell4.iCBsCa.Conc = iCBsCa_conc
sim_WM.shell4.iCBCaf.Conc = iCBCaf_conc
sim_WM.shell4.iCBCaCa.Conc = iCBCaCa_conc

sim_WM.shell4.CBsf.Conc = CBsf_conc
sim_WM.shell4.CBsCa.Conc = CBsCa_conc
sim_WM.shell4.CBCaf.Conc = CBCaf_conc
sim_WM.shell4.CBCaCa.Conc = CBCaCa_conc

sim_WM.shell4.PV.Conc = PV_conc
sim_WM.shell4.PVCa.Conc = PVCa_conc
sim_WM.shell4.PVMg.Conc = PVMg_conc

sim_WM.shell5.iCBsf.Conc = iCBsf_conc
sim_WM.shell5.iCBsCa.Conc = iCBsCa_conc
sim_WM.shell5.iCBCaf.Conc = iCBCaf_conc
sim_WM.shell5.iCBCaCa.Conc = iCBCaCa_conc

sim_WM.shell5.CBsf.Conc = CBsf_conc
sim_WM.shell5.CBsCa.Conc = CBsCa_conc
sim_WM.shell5.CBCaf.Conc = CBCaf_conc
sim_WM.shell5.CBCaCa.Conc = CBCaCa_conc

sim_WM.shell5.PV.Conc = PV_conc
sim_WM.shell5.PVCa.Conc = PVCa_conc
sim_WM.shell5.PVMg.Conc = PVMg_conc


# CaP
sim_WM.ring0.CaP_m0.Count = round(CaP_ro*surfarea*CaP_m0_p)
sim_WM.ring0.CaP_m1.Count = round(CaP_ro*surfarea*CaP_m1_p)
sim_WM.ring0.CaP_m2.Count = round(CaP_ro*surfarea*CaP_m2_p)
sim_WM.ring0.CaP_m3.Count = round(CaP_ro*surfarea*CaP_m3_p)

print("Injected  ", CaP_ro*surfarea, "CaP channels")

# CaT

# From cstate: CaT_m2h0 conducting
sim_WM.ring0.CaT_m0h0.Count = round(CaT_ro*surfarea*CaT_m0h0_p)
sim_WM.ring0.CaT_m1h0.Count = round(CaT_ro*surfarea*CaT_m1h0_p)
sim_WM.ring0.CaT_m2h0.Count = round(CaT_ro*surfarea*CaT_m2h0_p)
sim_WM.ring0.CaT_m0h1.Count = round(CaT_ro*surfarea*CaT_m0h1_p)
sim_WM.ring0.CaT_m1h1.Count = round(CaT_ro*surfarea*CaT_m1h1_p)
sim_WM.ring0.CaT_m2h1.Count = round(CaT_ro*surfarea*CaT_m2h1_p)

print("Injected  ", CaT_ro*surfarea, "CaT channels")

# BK
sim_WM.ring0.BK_C0.Count = round(BK_ro*surfarea*BK_C0_p)
sim_WM.ring0.BK_C1.Count = round(BK_ro*surfarea*BK_C1_p)
sim_WM.ring0.BK_C2.Count = round(BK_ro*surfarea*BK_C2_p)
sim_WM.ring0.BK_C3.Count = round(BK_ro*surfarea*BK_C3_p)
sim_WM.ring0.BK_C4.Count = round(BK_ro*surfarea*BK_C4_p)

sim_WM.ring0.BK_O0.Count = round(BK_ro*surfarea*BK_O0_p)
sim_WM.ring0.BK_O1.Count = round(BK_ro*surfarea*BK_O1_p)
sim_WM.ring0.BK_O2.Count = round(BK_ro*surfarea*BK_O2_p)
sim_WM.ring0.BK_O3.Count = round(BK_ro*surfarea*BK_O3_p)
sim_WM.ring0.BK_O4.Count = round(BK_ro*surfarea*BK_O4_p)


print("Injected  ", BK_ro*surfarea, "BK channels")

# SK
sim_WM.ring0.SK_C1.Count = round(SK_ro*surfarea*SK_C1_p)
sim_WM.ring0.SK_C2.Count = round(SK_ro*surfarea*SK_C2_p)
sim_WM.ring0.SK_C3.Count = round(SK_ro*surfarea*SK_C3_p)
sim_WM.ring0.SK_C4.Count = round(SK_ro*surfarea*SK_C4_p)

sim_WM.ring0.SK_O1.Count = round(SK_ro*surfarea*SK_O1_p)
sim_WM.ring0.SK_O2.Count = round(SK_ro*surfarea*SK_O2_p)

#set the rate constants for diffusion (Diffusion is modeled as surface reaction here)

#for Ca diffusion

sim_WM.ring1.diff_Ca_inward['fwd'].K = (DCST*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_Ca_inward['fwd'].K = (DCST*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_Ca_inward['fwd'].K = (DCST*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_Ca_inward['fwd'].K = (DCST*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_Ca_inward['fwd'].K = (DCST*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))

sim_WM.ring1.diff_Ca_inward['bkw'].K = (DCST*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_Ca_inward['bkw'].K = (DCST*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_Ca_inward['bkw'].K = (DCST*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_Ca_inward['bkw'].K = (DCST*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_Ca_inward['bkw'].K = (DCST*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6))

#for CBsf diffusin
sim_WM.ring1.diff_CBsf_inward['fwd'].K = (DCB*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_CBsf_inward['fwd'].K = (DCB*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_CBsf_inward['fwd'].K = (DCB*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_CBsf_inward['fwd'].K = (DCB*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_CBsf_inward['fwd'].K = (DCB*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))

sim_WM.ring1.diff_CBsf_inward['bkw'].K = (DCB*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_CBsf_inward['bkw'].K = (DCB*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_CBsf_inward['bkw'].K = (DCB*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_CBsf_inward['bkw'].K = (DCB*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_CBsf_inward['bkw'].K = (DCB*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6))

#for CBsCa diffusion

sim_WM.ring1.diff_CBsCa_inward['fwd'].K = (DCB*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_CBsCa_inward['fwd'].K = (DCB*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_CBsCa_inward['fwd'].K = (DCB*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_CBsCa_inward['fwd'].K = (DCB*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_CBsCa_inward['fwd'].K = (DCB*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))

sim_WM.ring1.diff_CBsCa_inward['bkw'].K = (DCB*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_CBsCa_inward['bkw'].K = (DCB*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_CBsCa_inward['bkw'].K = (DCB*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_CBsCa_inward['bkw'].K = (DCB*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_CBsCa_inward['bkw'].K = (DCB*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6))

#for CBCaf diffusion

sim_WM.ring1.diff_CBCaf_inward['fwd'].K = (DCB*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_CBCaf_inward['fwd'].K = (DCB*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_CBCaf_inward['fwd'].K = (DCB*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_CBCaf_inward['fwd'].K = (DCB*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_CBCaf_inward['fwd'].K = (DCB*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))

sim_WM.ring1.diff_CBCaf_inward['bkw'].K = (DCB*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_CBCaf_inward['bkw'].K = (DCB*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_CBCaf_inward['bkw'].K = (DCB*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_CBCaf_inward['bkw'].K = (DCB*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_CBCaf_inward['bkw'].K = (DCB*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6))

#for CBCaCa diffusion

sim_WM.ring1.diff_CBCaCa_inward['fwd'].K = (DCB*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_CBCaCa_inward['fwd'].K = (DCB*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_CBCaCa_inward['fwd'].K = (DCB*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_CBCaCa_inward['fwd'].K = (DCB*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_CBCaCa_inward['fwd'].K = (DCB*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))

sim_WM.ring1.diff_CBCaCa_inward['bkw'].K = (DCB*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_CBCaCa_inward['bkw'].K = (DCB*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_CBCaCa_inward['bkw'].K = (DCB*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_CBCaCa_inward['bkw'].K = (DCB*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_CBCaCa_inward['bkw'].K = (DCB*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6))

#for PV diffusion

sim_WM.ring1.diff_PV_inward['fwd'].K = (DPV*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_PV_inward['fwd'].K = (DPV*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_PV_inward['fwd'].K = (DPV*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_PV_inward['fwd'].K = (DPV*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_PV_inward['fwd'].K = (DPV*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))

sim_WM.ring1.diff_PV_inward['bkw'].K = (DPV*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_PV_inward['bkw'].K = (DPV*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_PV_inward['bkw'].K = (DPV*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_PV_inward['bkw'].K = (DPV*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_PV_inward['bkw'].K = (DPV*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6))

#for PVCa diffusion

sim_WM.ring1.diff_PVCa_inward['fwd'].K = (DPV*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_PVCa_inward['fwd'].K = (DPV*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_PVCa_inward['fwd'].K = (DPV*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_PVCa_inward['fwd'].K = (DPV*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_PVCa_inward['fwd'].K = (DPV*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))

sim_WM.ring1.diff_PVCa_inward['bkw'].K = (DPV*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_PVCa_inward['bkw'].K = (DPV*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_PVCa_inward['bkw'].K = (DPV*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_PVCa_inward['bkw'].K = (DPV*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_PVCa_inward['bkw'].K = (DPV*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6))

#for PVMg diffusion

sim_WM.ring1.diff_PVMg_inward['fwd'].K = (DPV*2*shell1_rmax)/((shell0_rmax**2 - shell1_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_PVMg_inward['fwd'].K = (DPV*2*shell2_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_PVMg_inward['fwd'].K = (DPV*2*shell3_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_PVMg_inward['fwd'].K = (DPV*2*shell4_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_PVMg_inward['fwd'].K = (DPV*2*shell5_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))

sim_WM.ring1.diff_PVMg_inward['bkw'].K = (DPV*2*shell1_rmax)/((shell1_rmax**2 - shell2_rmax**2)*(0.2e-6))
sim_WM.ring2.diff_PVMg_inward['bkw'].K = (DPV*2*shell2_rmax)/((shell2_rmax**2 - shell3_rmax**2)*(0.2e-6))
sim_WM.ring3.diff_PVMg_inward['bkw'].K = (DPV*2*shell3_rmax)/((shell3_rmax**2 - shell4_rmax**2)*(0.2e-6))
sim_WM.ring4.diff_PVMg_inward['bkw'].K = (DPV*2*shell4_rmax)/((shell4_rmax**2 - shell5_rmax**2)*(0.2e-6))
sim_WM.ring5.diff_PVMg_inward['bkw'].K = (DPV*2*shell5_rmax)/((shell5_rmax**2)*(0.2e-6))

print("Injected ", SK_ro*surfarea, "SK channels")

sim_stoch.memb_stoch.L[Leak].Count = round(L_ro * surfarea)
print("Injected  ", (L_ro * sim_stoch.memb_stoch.Area), "Leak channels")


sim_stoch.EfieldDT = EF_DT
sim_stoch.membrane.Potential = init_pot
sim_stoch.membrane.VolRes = Ra
#cm = 1.5uF/cm2 -> 1.5e-6F/1e-4m2 ->1.5e-2 F/m2
sim_stoch.membrane.Capac = memb_capac

for l in range(NTIMEPOINTS):
    print("Tpnt: ", l)
    # IN this sim V should be constant everywhere
    V = sim_stoch.TRI(memb_tris[0]).V
    
    sim_WM.ring0.CaPm0m1['fwd'].K = 1.0e3 *3.* alpha_cap(V*1.0e3)*Qt
    sim_WM.ring0.CaPm1m2['fwd'].K = 1.0e3 *2.* alpha_cap(V*1.0e3)*Qt
    sim_WM.ring0.CaPm2m3['fwd'].K = 1.0e3 *1.* alpha_cap(V*1.0e3)*Qt
    sim_WM.ring0.CaPm2m3['bkw'].K = 1.0e3 *3.* beta_cap(V*1.0e3)*Qt
    sim_WM.ring0.CaPm1m2['bkw'].K = 1.0e3 *2.* beta_cap(V*1.0e3)*Qt
    sim_WM.ring0.CaPm0m1['bkw'].K = 1.0e3 *1.* beta_cap(V*1.0e3)*Qt
    
    sim_WM.ring0.CaTm0h0_m1h0['fwd'].K = 1.0e3 *2.* alpham_cat(V*1.0e3)
    sim_WM.ring0.CaTm1h0_m2h0['fwd'].K = 1.0e3 *1.* alpham_cat(V*1.0e3)
    
    sim_WM.ring0.CaTm1h0_m2h0['bkw'].K = 1.0e3 *2.* betam_cat(V*1.0e3)
    sim_WM.ring0.CaTm0h0_m1h0['bkw'].K = 1.0e3 *1.* betam_cat(V*1.0e3)
    
    sim_WM.ring0.CaTm0h0_m0h1['fwd'].K = 1.0e3 *1.* alphah_cat(V*1.0e3)
    sim_WM.ring0.CaTm1h0_m1h1['fwd'].K = 1.0e3 *1.* alphah_cat(V*1.0e3)
    sim_WM.ring0.CaTm2h0_m2h1['fwd'].K = 1.0e3 *1.* alphah_cat(V*1.0e3)
    
    sim_WM.ring0.CaTm2h0_m2h1['bkw'].K = 1.0e3 *1.* betah_cat(V*1.0e3)
    sim_WM.ring0.CaTm1h0_m1h1['bkw'].K = 1.0e3 *1.* betah_cat(V*1.0e3)
    sim_WM.ring0.CaTm0h0_m0h1['bkw'].K = 1.0e3 *1.* betah_cat(V*1.0e3)
    
    sim_WM.ring0.CaTm0h1_m1h1['fwd'].K = 1.0e3 *2.* alpham_cat(V*1.0e3)
    sim_WM.ring0.CaTm1h1_m2h1['fwd'].K = 1.0e3 *1.* alpham_cat(V*1.0e3)
    
    sim_WM.ring0.CaTm1h1_m2h1['bkw'].K = 1.0e3 *2.* betam_cat(V*1.0e3)
    sim_WM.ring0.CaTm0h1_m1h1['bkw'].K = 1.0e3 *1.* betam_cat(V*1.0e3)
    
    
    sim_WM.ring0.BKC0O0['fwd'].K = f_0(V)
    sim_WM.ring0.BKC1O1['fwd'].K = f_1(V)
    sim_WM.ring0.BKC2O2['fwd'].K = f_2(V)
    sim_WM.ring0.BKC3O3['fwd'].K = f_3(V)
    sim_WM.ring0.BKC4O4['fwd'].K = f_4(V)
    sim_WM.ring0.BKC0O0['bkw'].K = b_0(V)
    sim_WM.ring0.BKC1O1['bkw'].K = b_1(V)
    sim_WM.ring0.BKC2O2['bkw'].K = b_2(V)
    sim_WM.ring0.BKC3O3['bkw'].K = b_3(V)
    sim_WM.ring0.BKC4O4['bkw'].K = b_4(V)
    
    sim_WM.run(TIMECONVERTER*l)

    # Now do the communication between the sims
    
    Si = sim_WM.shell0.Ca.Conc

    So = Ca_oconc
    
    # Get the single-channel currents first
    tcur_CaP_sc = cf.getGHKI(CaP_P, V, 2, TEMPERATURE+273.15, Si*1.0e3, So*1.0e3)
    tcur_CaT_sc = cf.getGHKI(CaT_P, V, 2, TEMPERATURE+273.15, Si*1.0e3, So*1.0e3)
    tcur_BK_sc = cf.getOhmI(V, BK_rev, BK_G)
    tcur_SK_sc = cf.getOhmI(V, SK_rev, SK_G)
    
    tcur_CaP = tcur_CaP_sc*sim_WM.ring0.CaP_m3.Count
    # alpha is to h1
    tcur_CaT = tcur_CaT_sc*sim_WM.ring0.CaT_m2h1.Count
    tcur_BK = tcur_BK_sc * sum(sim_WM.ring0.MATCH('^BK_O[0-4]$').Count)
    tcur_SK = tcur_SK_sc * sum(sim_WM.ring0.MATCH('^SK_O[1-2]$').Count)

    nmtris = len(memb_tris)
    
    for tri in memb_tris: 
        sim_stoch.TRI(tri).IClamp = (tcur_CaP+tcur_CaT+tcur_BK+tcur_SK)/nmtris
    
    sim_stoch.run(TIMECONVERTER*l)
    
    ca_count_inj =  -1.0*((tcur_CaP+tcur_CaT)*TIMECONVERTER)/(2*E_CHARGE)
    actual_inj= 0.0
    
    sim_WM.shell0.Ca.Count += ca_count_inj
    
    ca_shell = sim_WM.shell0.Ca.Count
    
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
    datfile3.write('%.6g' %sim_WM.shell0.Ca.Count +' ')
    datfile3.write('%.6g' %(Si*1.0e6)+ ' ')
    datfile3.write('\n')

datfile.close()
datfile2.close()
datfile3.close()

