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

import steps.interface

import math
import time
import random
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

_, meshfile_ab, root, iter_n, clusterSize = sys.argv

if meshfile_ab == 'Cylinder2_dia2um_L160um_outer0_0.3shell_0.3size_279152tets_adaptive.inp':
    cyl160=True
else:
    cyl160=False

clusterSize = int(clusterSize)

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
        r['iCBsf1_f'].setRates(iCBsf1_f_kcst, iCBsf1_b_kcst)
        r['iCBsCa_f'].setRates(iCBsCa_f_kcst, iCBsCa_b_kcst)
        r['iCBsf2_f'].setRates(iCBsf2_f_kcst, iCBsf2_b_kcst)
        r['iCBCaf_f'].setRates(iCBCaf_f_kcst, iCBCaf_b_kcst)

        (CBsf + Ca <r['CBsf1_f']> CBsCa) + Ca <r['CBsCa_f']> CBCaCa
        (CBsf + Ca <r['CBsf2_f']> CBCaf) + Ca <r['CBCaf_f']> CBCaCa
        r['CBsf1_f'].setRates(CBsf1_f_kcst, CBsf1_b_kcst)
        r['CBsCa_f'].setRates(CBsCa_f_kcst, CBsCa_b_kcst)
        r['CBsf2_f'].setRates(CBsf2_f_kcst, CBsf2_b_kcst)
        r['CBCaf_f'].setRates(CBCaf_f_kcst, CBCaf_b_kcst)

        Ca + PV <r['PVca_f']> PVCa
        Mg + PV <r['PVmg_f']> PVMg
        r['PVca_f'].setRates(PVca_f_kcst, PVca_b_kcst)
        r['PVmg_f'].setRates(PVmg_f_kcst, PVmg_b_kcst)

    with ssys:
    
        #Pump
        Ca.i + Pump.s <r['PumpD_f']> CaPump.s >r['PumpD_k']> Pump.s
        r['PumpD_f'].setRates(P_f_kcst, P_b_kcst)
        r['PumpD_k'].setRates(P_k_kcst)

        with CaPchan[...]:
            CaP_m0.s <r['CaPm0m1']> CaP_m1.s <r['CaPm1m2']> CaP_m2.s <r['CaPm2m3']> CaP_m3.s
            r['CaPm0m1'].setRates(VDepRate(lambda V: 1.0e3 *3.* alpha_cap(V*1.0e3)* Qt), VDepRate(lambda V: 1.0e3 *1.* beta_cap(V*1.0e3)* Qt))
            r['CaPm1m2'].setRates(VDepRate(lambda V: 1.0e3 *2.* alpha_cap(V*1.0e3)* Qt), VDepRate(lambda V: 1.0e3 *2.* beta_cap(V*1.0e3)* Qt))
            r['CaPm2m3'].setRates(VDepRate(lambda V: 1.0e3 *1.* alpha_cap(V*1.0e3)* Qt), VDepRate(lambda V: 1.0e3 *3.* beta_cap(V*1.0e3)* Qt))

        with CaTchan[...]:
            CaT_m0h0.s <r['CaTm0h0_m1h0']> CaT_m1h0.s <r['CaTm1h0_m2h0']> CaT_m2h0.s <r['CaTm2h0_m2h1']> CaT_m2h1.s
            r['CaTm0h0_m1h0'].setRates(VDepRate(lambda V: 1.0e3 *2.* alpham_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *1.* betam_cat(V*1.0e3)))
            r['CaTm1h0_m2h0'].setRates(VDepRate(lambda V: 1.0e3 *1.* alpham_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *2.* betam_cat(V*1.0e3)))
            r['CaTm2h0_m2h1'].setRates(VDepRate(lambda V: 1.0e3 *1.* alphah_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *1.* betah_cat(V*1.0e3)))

            CaT_m1h0.s <r['CaTm1h0_m1h1']> CaT_m1h1.s
            r['CaTm1h0_m1h1'].setRates(VDepRate(lambda V: 1.0e3 *1.* alphah_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *1.* betah_cat(V*1.0e3)))
            
            CaT_m0h0.s <r['CaTm0h0_m0h1']> CaT_m0h1.s <r['CaTm0h1_m1h1']> CaT_m1h1.s <r['CaTm1h1_m2h1']> CaT_m2h1.s
            r['CaTm0h0_m0h1'].setRates(VDepRate(lambda V: 1.0e3 *1.* alphah_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *1.* betah_cat(V*1.0e3)))
            r['CaTm1h1_m2h1'].setRates(VDepRate(lambda V: 1.0e3 *1.* alpham_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *2.* betam_cat(V*1.0e3)))
            r['CaTm0h1_m1h1'].setRates(VDepRate(lambda V: 1.0e3 *2.* alpham_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *1.* betam_cat(V*1.0e3)))

        with BKchan[...]:
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
            r['BKC0O0'].setRates(VDepRate(f_0), VDepRate(b_0))
            r['BKC1O1'].setRates(VDepRate(f_1), VDepRate(b_1))
            r['BKC2O2'].setRates(VDepRate(f_2), VDepRate(b_2))
            r['BKC3O3'].setRates(VDepRate(f_3), VDepRate(b_3))
            r['BKC4O4'].setRates(VDepRate(f_4), VDepRate(b_4))

        with SKchan[...]:
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

    surfarea = sum(tri.Area for tri in memb_tris)
    random.seed(7)

    memb_tris_bk = TriList(random.sample(memb_tris.indices, round(BK_ro*surfarea)))
    memb_tris_sk = TriList(random.sample(memb_tris.indices, round(SK_ro*surfarea)))
    memb_tris_cat = TriList(random.sample(memb_tris.indices, round(CaT_ro*surfarea)))
    memb_tris_cap = TriList(random.choices((memb_tris - memb_tris_bk).indices, k=round(CaP_ro*surfarea)-(clusterSize*round(BK_ro*surfarea))))

    ########## Find the submembrane tets
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

    border_tets_vols = sum(tet.Vol for tet in set(submemb_tets))
    print("Border tet vols:", border_tets_vols)

    ########## Create a membrane as a surface mesh
    if cyl160: 
        memb = Patch.Create(memb_tris, cyto, None, ssys)
    else:
        memb = Patch.Create(memb_tris, cyto, outer, ssys)

    # For EField calculation
    print("Creating membrane..")
    membrane = Membrane.Create([memb])
    print("Membrane created.")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

seed = 7
#TODO
# seed = int(time.time()%1000)
rng = RNG('mt19937', 512, seed)

#_________ mdl.stepsModel __________

print("Creating tetexact solver...")
sim = Simulation('Tetexact', mdl, mesh, rng, calcMembPot=True)

#### Recording #####

dc = time.strftime('%b%d_%H_%M_%S_%Y')

runPath = os.path.join(root, 'data/StochasticCaburst_cluster/', meshfile_ab, f'{iter_n}__{dc}')
os.makedirs(runPath, exist_ok=True)

rs = ResultSelector(sim)

rs1 = rs.SUM(rs.TRIS(memb_tris).OC_CaP.I) << \
      rs.SUM(rs.TRIS(memb_tris).OC_CaT.I) << \
      rs.SUM(rs.TRIS(memb_tris).OC_BK.I) << \
      rs.SUM(rs.TRIS(memb_tris).OC_SK.I)

rs2 = rs.TET(cent_tet).V

rs3 = rs.SUM(rs.TETS(submemb_tets).Ca.Count)

rs4 = rs.TRIS(memb_tris).LIST(*BKchan[BK_O0|BK_O1|BK_O2|BK_O3|BK_O4]).Count << rs.TETS(submemb_tets).Ca.Count

rs5 = rs.TRIS(memb_tris).LIST(*CaPchan, *CaTchan, *BKchan, *SKchan).Count

rs1.toFile(os.path.join(runPath, 'currents.dat.bin'))
rs2.toFile(os.path.join(runPath, 'voltage.dat.bin'))
rs3.toFile(os.path.join(runPath, 'calcium.dat.bin'))
rs4.toFile(os.path.join(runPath, 'OpenBKandCa.dat.bin'))
rs5.toFile(os.path.join(runPath, 'ChannelsDistribution.dat.bin'))

sim.toSave(rs1, rs2, rs3, rs4, dt=TIMECONVERTER)
sim.toSave(rs5)

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

sim.memb.Pump.Count = pumpnbs
sim.memb.CaPump.Count = 0

print("Injected ", sim.memb.Pump.Count, "pumps")

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

def initChannels(sim, tris, fracs, subs, ro, chan, cs):
    cnts = [sim.memb.LIST(chan)[sb].Count for sb in subs]
    for tri in tris:
        for c in range(cs):
            added = False
            for i in range(len(cnts)):
                if cnts[i] < round(ro * surfarea * fracs[i]):
                    sim.TRI(tri).LIST(chan)[subs[i]].Count += 1
                    cnts[i] += 1
                    added = True
                    break
            if not added:
                print(f'More tris picked up by algorithm than the number of {chan}')

initChannels(
    sim, memb_tris_bk, 
    [BK_C0_p, BK_C1_p, BK_C2_p, BK_C3_p, BK_C4_p, BK_O1_p, BK_O2_p, BK_O3_p, BK_O4_p], 
    [BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O1, BK_O2, BK_O3, BK_O4],
    BK_ro, BKchan, 1
)

initChannels(
    sim, memb_tris_sk, 
    [SK_C1_p, SK_C2_p, SK_C3_p, SK_C4_p, SK_O1_p, SK_O2_p], 
    [SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2], 
    SK_ro, SKchan, 1
)

initChannels(
    sim, memb_tris_cat, 
    [CaT_m0h0_p, CaT_m1h0_p, CaT_m2h0_p, CaT_m0h1_p, CaT_m1h1_p, CaT_m2h1_p], 
    [CaT_m0h0, CaT_m1h0, CaT_m2h0, CaT_m0h1, CaT_m1h1, CaT_m2h1], 
    CaT_ro, CaTchan, 1
)

initChannels(
    sim, memb_tris_bk, 
    [CaP_m0_p, CaP_m1_p, CaP_m2_p, CaP_m3_p], 
    [CaP_m0, CaP_m1, CaP_m2, CaP_m3], 
    CaP_ro, CaPchan, clusterSize
)

initChannels(
    sim, memb_tris_cap, 
    [CaP_m0_p, CaP_m1_p, CaP_m2_p, CaP_m3_p], 
    [CaP_m0, CaP_m1, CaP_m2, CaP_m3], 
    CaP_ro, CaPchan, 1
)

sim.memb.L[Leak].Count = int(L_ro * surfarea)
print("Injected  ", int(L_ro * sim.memb.Area), "Leak channels")

memb_countBK_pertriID = sim.TRIS(memb_tris).BKchan[...].Count

sim.EfieldDT = EF_DT
sim.membrane.Potential = init_pot
sim.membrane.VolRes = Ra
#cm = 1.5uF/cm2 -> 1.5e-6F/1e-4m2 ->1.5e-2 F/m2
sim.membrane.Capac = memb_capac

# Save the rs5 information only once
rs5.save()

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
        for val in row:
            f.write('%.6g' % (val / AVOGADRO / (border_tets_vols*1.0e3) * 1e6) + ' ')
            f.write('%.6g' % val + ' ')
        f.write('\n')

with open(os.path.join(runPath, 'ChannelsDistribution.dat'), 'w') as f:
    for t, row in zip(rs4.time[0], rs4.data[0]):
        f.write('%.6g' % (t * 1e3) + ' ')
        n = len(memb_tris)
        for i, tri in enumerate(memb_tris):
            if memb_countBK_pertriID[i] > 0:
                for j in range(5):
                    f.write('%.6g' % row[i*5 + j] + ' ')
                f.write('%.6g' % row[n * 5 + i] + ' ')
                f.write('%.6g' % (row[n * 5 + i]/AVOGADRO/submemb_tets[i].Vol) + ' ')
        f.write('\n')

n = len(list(CaPchan)) + len(list(CaTchan)) + len(list(BKchan)) + len(list(SKchan))
with open(os.path.join(runPath, 'OpenBKandCa.dat'), 'w') as f:
    for t, row in zip(rs5.time[0], rs5.data[0]):
        for i, val in enumerate(row):
            if i % n == 0:
                for pos in memb_tris[i // n].center:
                    f.write('%.6g' % (pos) + ' ')
            f.write('%.6g' % val + ' ')
        f.write('\n')

