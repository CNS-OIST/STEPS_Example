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

import sys

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

_, meshfile_ab, root, iter_n = sys.argv

if meshfile_ab == 'Cylinder2_dia2um_L160um_outer0_0.3shell_0.3size_279152tets_adaptive.inp':
    cyl160=True
else:
    cyl160=False

########################### BIOCHEMICAL MODEL ###############################

mdl = Model()
r = ReactionManager()

with mdl:
    
    # Calcium
    Ca = Species.Create(valence=2)
    
    # Species
    Pump, CaPump, iCBsf, iCBsCa, iCBCaf, iCBCaCa, CBsf, CBsCa, CBCaf, CBCaCa, PV, PVMg, PVCa, Mg = Species.Create()
    
    # Vol/surface systems
    vsys = VolumeSystem.Create()
    ssys = SurfaceSystem.Create()

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # CHANNELS  # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    ###### CaP channel ############## 
    CaP_m0, CaP_m1, CaP_m2, CaP_m3 = SubUnitState.Create()
    CaPchan = Channel.Create([CaP_m0, CaP_m1, CaP_m2, CaP_m3])

    ######## CaT channel ##########  
    CaT_m0h0, CaT_m0h1, CaT_m1h0, CaT_m1h1, CaT_m2h0, CaT_m2h1 = SubUnitState.Create()
    CaTchan = Channel.Create([CaT_m0h0, CaT_m0h1, CaT_m1h0, CaT_m1h1, CaT_m2h0, CaT_m2h1])

    ##### BK channel ####################
    BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O0, BK_O1, BK_O2, BK_O3, BK_O4 = SubUnitState.Create()
    BKchan = Channel.Create([BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O0, BK_O1, BK_O2, BK_O3, BK_O4])

    ###### SK channel ##################
    SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2 = SubUnitState.Create()
    SKchan = Channel.Create([SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2])

    ###### Leak current channel #####
    Leak = SubUnitState.Create()
    L = Channel.Create([Leak])

    with vsys:
        
        # Diffusions 
        diff_Ca =     Diffusion.Create(Ca, DCST)
        diff_CBsf =   Diffusion.Create(CBsf, DCB)
        diff_CBsCa =  Diffusion.Create(CBsCa, DCB)
        diff_CBCaf =  Diffusion.Create(CBCaf, DCB)
        diff_CBCaCa = Diffusion.Create(CBCaCa, DCB)
        diff_PV =     Diffusion.Create(PV, DPV)
        diff_PVCa =   Diffusion.Create(PVCa, DPV)
        diff_PVMg =   Diffusion.Create(PVMg, DPV)

        (Ca + iCBsf <r['iCBsf1_f']> iCBsCa) + Ca <r['iCBsCa_f']> iCBCaCa
        (Ca + iCBsf <r['iCBsf2_f']> iCBCaf) + Ca <r['iCBCaf_f']> iCBCaCa
        r['iCBsf1_f'].K = iCBsf1_f_kcst, iCBsf1_b_kcst
        r['iCBsCa_f'].K = iCBsCa_f_kcst, iCBsCa_b_kcst
        r['iCBsf2_f'].K = iCBsf2_f_kcst, iCBsf2_b_kcst
        r['iCBCaf_f'].K = iCBCaf_f_kcst, iCBCaf_b_kcst

        (CBsf + Ca <r['CBsf1_f']> CBsCa) + Ca <r['CBsCa_f']> CBCaCa
        (CBsf + Ca <r['CBsf2_f']> CBCaf) + Ca <r['CBCaf_f']> CBCaCa
        r['CBsf1_f'].K = CBsf1_f_kcst, CBsf1_b_kcst
        r['CBsCa_f'].K = CBsCa_f_kcst, CBsCa_b_kcst
        r['CBsf2_f'].K = CBsf2_f_kcst, CBsf2_b_kcst
        r['CBCaf_f'].K = CBCaf_f_kcst, CBCaf_b_kcst

        Ca + PV <r['PVca_f']> PVCa
        Mg + PV <r['PVmg_f']> PVMg
        r['PVca_f'].K = PVca_f_kcst, PVca_b_kcst
        r['PVmg_f'].K = PVmg_f_kcst, PVmg_b_kcst

    with ssys:
    
        #Pump
        Ca.i + Pump.s <r['PumpD_f']> CaPump.s >r['PumpD_k']> Pump.s
        r['PumpD_f'].K = P_f_kcst, P_b_kcst
        r['PumpD_k'].K = P_k_kcst

        with CaPchan[...]:
            CaP_m0.s <r['CaPm0m1']> CaP_m1.s <r['CaPm1m2']> CaP_m2.s <r['CaPm2m3']> CaP_m3.s
            r['CaPm0m1'].K = VDepRate(lambda V: 1.0e3 *3.* alpha_cap(V*1.0e3)* Qt), VDepRate(lambda V: 1.0e3 *1.* beta_cap(V*1.0e3)* Qt)
            r['CaPm1m2'].K = VDepRate(lambda V: 1.0e3 *2.* alpha_cap(V*1.0e3)* Qt), VDepRate(lambda V: 1.0e3 *2.* beta_cap(V*1.0e3)* Qt)
            r['CaPm2m3'].K = VDepRate(lambda V: 1.0e3 *1.* alpha_cap(V*1.0e3)* Qt), VDepRate(lambda V: 1.0e3 *3.* beta_cap(V*1.0e3)* Qt)

        with CaTchan[...]:
            CaT_m0h0.s <r['CaTm0h0_m1h0']> CaT_m1h0.s <r['CaTm1h0_m2h0']> CaT_m2h0.s <r['CaTm2h0_m2h1']> CaT_m2h1.s
            r['CaTm0h0_m1h0'].K = VDepRate(lambda V: 1.0e3 *2.* alpham_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *1.* betam_cat(V*1.0e3))
            r['CaTm1h0_m2h0'].K = VDepRate(lambda V: 1.0e3 *1.* alpham_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *2.* betam_cat(V*1.0e3))
            r['CaTm2h0_m2h1'].K = VDepRate(lambda V: 1.0e3 *1.* alphah_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *1.* betah_cat(V*1.0e3))

            CaT_m1h0.s <r['CaTm1h0_m1h1']> CaT_m1h1.s
            r['CaTm1h0_m1h1'].K = VDepRate(lambda V: 1.0e3 *1.* alphah_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *1.* betah_cat(V*1.0e3))
            
            CaT_m0h0.s <r['CaTm0h0_m0h1']> CaT_m0h1.s <r['CaTm0h1_m1h1']> CaT_m1h1.s <r['CaTm1h1_m2h1']> CaT_m2h1.s
            r['CaTm0h0_m0h1'].K = VDepRate(lambda V: 1.0e3 *1.* alphah_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *1.* betah_cat(V*1.0e3))
            r['CaTm1h1_m2h1'].K = VDepRate(lambda V: 1.0e3 *1.* alpham_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *2.* betam_cat(V*1.0e3))
            r['CaTm0h1_m1h1'].K = VDepRate(lambda V: 1.0e3 *2.* alpham_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *1.* betam_cat(V*1.0e3))

        with BKchan[...]:
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
            r['BKC0O0'].K = VDepRate(f_0), VDepRate(b_0)
            r['BKC1O1'].K = VDepRate(f_1), VDepRate(b_1)
            r['BKC2O2'].K = VDepRate(f_2), VDepRate(b_2)
            r['BKC3O3'].K = VDepRate(f_3), VDepRate(b_3)
            r['BKC4O4'].K = VDepRate(f_4), VDepRate(b_4)

        with SKchan[...]:
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
        
        if cyl160:
            OC_CaP = GHKCurr.Create(CaPchan[CaP_m3], Ca, CaP_P, virtual_oconc=Ca_oconc, computeflux=True)
            OC_CaT = GHKCurr.Create(CaTchan[CaT_m2h1], Ca, CaT_P, virtual_oconc=Ca_oconc, computeflux=True)
        else:
            OC_CaP = GHKCurr.Create(CaPchan[CaP_m3], Ca, CaP_P, computeflux=True)
            OC_CaT = GHKCurr.Create(CaTchan[CaT_m2h1], Ca, CaT_P, computeflux=True)

        OC_BK = OhmicCurr.Create(BKchan[BK_O0|BK_O1|BK_O2|BK_O3|BK_O4], BK_G, BK_rev)
        OC_SK = OhmicCurr.Create(SKchan[SK_O1|SK_O2], SK_G, SK_rev)

        OC_L = OhmicCurr.Create(L[Leak], L_G, L_rev)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

########### MESH & COMPARTMENTALIZATION #################

##########Import Mesh

mesh = TetMesh.Load('./meshes/'+meshfile_ab)

with mesh:
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

    ########## Create an intracellular compartment i.e. cytosolic compartment
    cyto = Compartment.Create(inner_tets, vsys)

    if not cyl160:
        outer = Compartment.Create(outer_tets)

    if cyl160:
        # Ensure that we use points a small distance inside the boundary:
        minz, maxz = mesh.bbox.min.z, mesh.bbox.max.z
        memb_tris = TriList(tri for tri in mesh_stock.surface if minz < tri.center.z < maxz)
    else:
        print('Finding connecting triangles...')
        memb_tris = inner_tets.surface & outer_tets.surface

    print(len(memb_tris), " surface triangles.")

    ########## Find the submembrane tets
    submemb_tets = TetList()
    for tri in memb_tris:
        submemb_tets += tri.tetNeighbs
    submemb_tets = submemb_tets & inner_tets

    print(len(submemb_tets))

    vol = sum(tet.Vol for tet in submemb_tets)
    print('Volume of submembrane region is', vol)

    ########## Create a membrane as a surface mesh
    if cyl160: 
        memb = Patch.Create(memb_tris, cyto, None, ssys)
    else:
        memb = Patch.Create(memb_tris, cyto, outer, ssys)

    print("Area: ", memb.Area)

    # For EField calculation
    print("Creating membrane..")
    membrane = Membrane.Create([memb], opt_file_name = './meshes/'+meshfile_ab+"_optimalidx")
    print("Membrane created.")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

rng = RNG('mt19937', 512, seed=int(time.time()%1000))

print("Creating tetexact solver...")
sim = Simulation('Tetexact', mdl, mesh, rng, calcMembPot=True)

#### Recording #####

dc = time.strftime('%b%d_%H_%M_%S_%Y')

runPath = os.path.join(root, 'data/StochasticCaburst/', meshfile_ab, f'{iter_n}__{dc}')
os.makedirs(runPath, exist_ok=True)

rs = ResultSelector(sim)

rs1 = rs.SUM(rs.TRIS(memb_tris).OC_CaP.I) << \
      rs.SUM(rs.TRIS(memb_tris).OC_CaT.I) << \
      rs.SUM(rs.TRIS(memb_tris).OC_BK.I) << \
      rs.SUM(rs.TRIS(memb_tris).OC_SK.I)

rs2 = rs.TET(cent_tet).V

rs3 = rs.cyto.Ca.Conc << rs.SUM(rs.TETS(submemb_tets).Ca.Count)

rs1.toFile(os.path.join(runPath, 'currents.dat.bin'))
rs2.toFile(os.path.join(runPath, 'voltage.dat.bin'))
rs3.toFile(os.path.join(runPath, 'calcium.dat.bin'))

sim.toSave(rs1, rs2, rs3, dt=TIMECONVERTER)

print("Resetting simulation objects..")
sim.newRun()

print("Injecting molecules..")

sim.Temp = TEMPERATURE + 273.15

if not cyl160: 
    sim.outer.Ca.Conc = Ca_oconc
    sim.outer.Ca.Clamped = True
    
sim.cyto.Ca.Conc = Ca_iconc

print("Calcium concentration is: ", sim.cyto.Ca.Conc)
print("No. of Ca molecules is: ", sim.cyto.Ca.Count)

sim.cyto.Mg.Conc = Mg_conc


surfarea = sim.memb.Area

pumpnbs = 6.022141e12*surfarea

sim.memb.Pump.Count = round(pumpnbs)
sim.memb.CaPump.Count = 0

sim.cyto.iCBsf.Conc = iCBsf_conc
sim.cyto.iCBsCa.Conc = iCBsCa_conc
sim.cyto.iCBCaf.Conc = iCBCaf_conc
sim.cyto.iCBCaCa.Conc = iCBCaCa_conc

sim.cyto.CBsf.Conc = CBsf_conc
sim.cyto.CBsCa.Conc = CBsCa_conc
sim.cyto.CBCaf.Conc = CBCaf_conc
sim.cyto.CBCaCa.Conc = CBCaCa_conc

sim.cyto.PV.Conc = PV_conc
sim.cyto.PVCa.Conc = PVCa_conc
sim.cyto.PVMg.Conc = PVMg_conc

sim.memb.CaPchan[CaP_m0].Count = round(CaP_ro*surfarea*CaP_m0_p)
sim.memb.CaPchan[CaP_m1].Count = round(CaP_ro*surfarea*CaP_m1_p)
sim.memb.CaPchan[CaP_m2].Count = round(CaP_ro*surfarea*CaP_m2_p)
sim.memb.CaPchan[CaP_m3].Count = round(CaP_ro*surfarea*CaP_m3_p)


print("CaP_m0 ", round(CaP_ro*surfarea*CaP_m0_p))
print("CaP_m1 ", round(CaP_ro*surfarea*CaP_m1_p))
print("CaP_m2 ", round(CaP_ro*surfarea*CaP_m2_p))
print("CaP_m3 ", round(CaP_ro*surfarea*CaP_m3_p))

print("Targeted Injection: ", round(CaP_ro*surfarea), "CaP channels")

sim.memb.CaTchan[CaT_m0h0].Count = round(CaT_ro*surfarea*CaT_m0h0_p)
sim.memb.CaTchan[CaT_m1h0].Count = round(CaT_ro*surfarea*CaT_m1h0_p)
sim.memb.CaTchan[CaT_m2h0].Count = round(CaT_ro*surfarea*CaT_m2h0_p)
sim.memb.CaTchan[CaT_m0h1].Count = round(CaT_ro*surfarea*CaT_m0h1_p)
sim.memb.CaTchan[CaT_m1h1].Count = round(CaT_ro*surfarea*CaT_m1h1_p)
sim.memb.CaTchan[CaT_m2h1].Count = round(CaT_ro*surfarea*CaT_m2h1_p)

print("m0h0", round(CaT_ro*surfarea*CaT_m0h0_p))
print("m1h0", round(CaT_ro*surfarea*CaT_m1h0_p))
print("m2h0", round(CaT_ro*surfarea*CaT_m2h0_p))
print("m0h1", round(CaT_ro*surfarea*CaT_m0h1_p))
print("m1h1", round(CaT_ro*surfarea*CaT_m1h1_p))
print("m2h1", round(CaT_ro*surfarea*CaT_m2h1_p))

print("Targeted Injection: ", round(CaT_ro*surfarea), "CaT channels")

sim.memb.BKchan[BK_C0].Count = round(BK_ro*surfarea*BK_C0_p)
sim.memb.BKchan[BK_C1].Count = round(BK_ro*surfarea*BK_C1_p)
sim.memb.BKchan[BK_C2].Count = round(BK_ro*surfarea*BK_C2_p)
sim.memb.BKchan[BK_C3].Count = round(BK_ro*surfarea*BK_C3_p)
sim.memb.BKchan[BK_C4].Count = round(BK_ro*surfarea*BK_C4_p)

sim.memb.BKchan[BK_O0].Count = round(BK_ro*surfarea*BK_O0_p)
sim.memb.BKchan[BK_O1].Count = round(BK_ro*surfarea*BK_O1_p)
sim.memb.BKchan[BK_O2].Count = round(BK_ro*surfarea*BK_O2_p)
sim.memb.BKchan[BK_O3].Count = round(BK_ro*surfarea*BK_O3_p)
sim.memb.BKchan[BK_O4].Count = round(BK_ro*surfarea*BK_O4_p)

print("BK_C0 ", round(BK_ro*surfarea*BK_C0_p))
print("BK_C1 ", round(BK_ro*surfarea*BK_C1_p))
print("BK_C2 ", round(BK_ro*surfarea*BK_C2_p))
print("BK_C3 ", round(BK_ro*surfarea*BK_C3_p))
print("BK_C4 ", round(BK_ro*surfarea*BK_C4_p))

print("BK_O0 ", round(BK_ro*surfarea*BK_O0_p))
print("BK_O1 ", round(BK_ro*surfarea*BK_O1_p))
print("BK_O2 ", round(BK_ro*surfarea*BK_O2_p))
print("BK_O3 ", round(BK_ro*surfarea*BK_O3_p))
print("BK_O4 ", round(BK_ro*surfarea*BK_O4_p))

print("Targeted Injection: ", round(BK_ro*surfarea), "BK channels")

sim.memb.SKchan[SK_C1].Count = round(SK_ro*surfarea*SK_C1_p)
sim.memb.SKchan[SK_C2].Count = round(SK_ro*surfarea*SK_C2_p)
sim.memb.SKchan[SK_C3].Count = round(SK_ro*surfarea*SK_C3_p)
sim.memb.SKchan[SK_C4].Count = round(SK_ro*surfarea*SK_C4_p)

sim.memb.SKchan[SK_O1].Count = round(SK_ro*surfarea*SK_O1_p)
sim.memb.SKchan[SK_O2].Count = round(SK_ro*surfarea*SK_O2_p)

print("SK_C1 ", round(SK_ro*surfarea*SK_C1_p))
print("SK_C2 ", round(SK_ro*surfarea*SK_C2_p))
print("SK_C3 ", round(SK_ro*surfarea*SK_C3_p))
print("SK_C4 ", round(SK_ro*surfarea*SK_C4_p))

print("SK_O1 ", round(SK_ro*surfarea*SK_O1_p))
print("SK_O2 ", round(SK_ro*surfarea*SK_O2_p))

print("Targeted Injection ", round(SK_ro*surfarea), "SK channels")

sim.memb.L[Leak].Count = int(L_ro * surfarea)
print("Injected ", int(L_ro * sim.memb.Area), "Leak channels")

sim.EfieldDT = EF_DT
sim.membrane.Potential = init_pot
sim.membrane.VolRes = Ra
#cm = 1.5uF/cm2 -> 1.5e-6F/1e-4m2 ->1.5e-2 F/m2
sim.membrane.Capac = memb_capac

for l in range(NTIMEPOINTS):
    print("Tpnt: ", l)
    sim.run(TIMECONVERTER*l)

# This last part is only present for backwards compatibility with the scripts created with API_1.
# We need to save to text files, like in the original script.

with open(os.path.join(runPath, 'currents.dat'), 'w') as f:
    for t, row in zip(rs1.time[0], rs1.data[0]):
        f.write('%.6g' % (t * 1e3) + ' ')
        for val in row:
            f.write('%.6g' % (val * 0.1 / surfarea) + ' ')
        f.write('\n')

with open(os.path.join(runPath, 'voltage.dat'), 'w') as f:
    for t, row in zip(rs2.time[0], rs2.data[0]):
        f.write('%.6g' % (t * 1e3) + ' ')
        for val in row:
            f.write('%.6g' % (val * 1e3) + ' ')
        f.write('\n')

with open(os.path.join(runPath, 'calcium.dat'), 'w') as f:
    for t, row in zip(rs3.time[0], rs3.data[0]):
        f.write('%.6g' % (t * 1e3) + ' ')
        conc, val = row
        f.write('%.6g' % (conc * 1e6) + ' ')
        f.write('%.6g' % val + ' ')
        f.write('\n')

