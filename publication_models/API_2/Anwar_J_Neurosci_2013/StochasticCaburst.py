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

from __future__ import print_function
import math
# WARNING: Using a variable name that is reserved (['time']).
import time
from random import *
from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *

import meshes.gettets as gettets
from extra.constants import *

import sys
import os

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

meshfile_ab, root, iter_n = sys.argv[1], sys.argv[2], sys.argv[3]

if meshfile_ab == 'Cylinder2_dia2um_L160um_outer0_0.3shell_0.3size_279152tets_adaptive.inp':
    cyl160=True
else:
    cyl160=False

########################### BIOCHEMICAL MODEL ###############################

mdl = Model()
# WARNING: Using a variable name that is reserved (['r']).
r = ReactionManager()
with mdl:
    
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
    
    # Pump
    Pump, CaPump, iCBsf, iCBsCa, iCBCaf, iCBCaCa, CBsf, CBsCa, CBCaf, CBCaCa, PV, PVMg, PVCa, Mg = Species.Create()
    
    # Vol/surface systems
    vsys = VolumeSystem.Create()
    ssys = SurfaceSystem.Create()
with vsys, mdl:
    
    diff_Ca =     Diffusion.Create(Ca, DCST)

    diff_CBsf =     Diffusion.Create(CBsf, DCB)

    diff_CBsCa =     Diffusion.Create(CBsCa, DCB)

    diff_CBCaf =     Diffusion.Create(CBCaf, DCB)

    diff_CBCaCa =     Diffusion.Create(CBCaCa, DCB)

    diff_PV =     Diffusion.Create(PV, DPV)

    diff_PVCa =     Diffusion.Create(PVCa, DPV)

    diff_PVMg =     Diffusion.Create(PVMg, DPV)

with ssys, mdl:
    
    #Pump
    
    Ca.i + Pump.s <r['PumpD_f']> CaPump.s
    r['PumpD_f'].setRates(P_f_kcst, P_b_kcst)
    
    CaPump.s >r['PumpD_k']> Pump.s
r['PumpD_k'].setRates(P_k_kcst)
with vsys, mdl:
    
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
with mdl:
    
    CaP_m0, CaP_m1, CaP_m2, CaP_m3 = SubUnitState.Create()
    
    
    ###### CaP channel ##############
    
    CaPchan = Channel.Create([CaP_m0, CaP_m1, CaP_m2, CaP_m3])
with ssys, mdl, CaPchan[...]:
    
    
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    CaP_m0.s <r['CaPm0m1']> CaP_m1.s ; r['CaPm0m1'].setRates(VDepRate(lambda V: 1.0e3 *3.* alpha_cap(V*1.0e3)* Qt), VDepRate(lambda V: 1.0e3 *1.* beta_cap(V*1.0e3)* Qt))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    CaP_m1.s <r['CaPm1m2']> CaP_m2.s ; r['CaPm1m2'].setRates(VDepRate(lambda V: 1.0e3 *2.* alpha_cap(V*1.0e3)* Qt), VDepRate(lambda V: 1.0e3 *2.* beta_cap(V*1.0e3)* Qt))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    
    CaP_m2.s <r['CaPm2m3']> CaP_m3.s ; r['CaPm2m3'].setRates(VDepRate(lambda V: 1.0e3 *1.* alpha_cap(V*1.0e3)* Qt), VDepRate(lambda V: 1.0e3 *3.* beta_cap(V*1.0e3)* Qt))

if cyl160:
    OC_CaP = smodel.GHKcurr('OC_CaP', ssys, CaP_m3, Ca, virtual_oconc = Ca_oconc, computeflux = True)
else:
    OC_CaP = GHKCurr.Create(CaPchan[CaP_m3], Ca, CaP_P, computeflux=True)
with mdl:
    
    CaT_m0h0, CaT_m0h1, CaT_m1h0, CaT_m1h1, CaT_m2h0, CaT_m2h1 = SubUnitState.Create()
    
    ######## CaT channel ##########
    
    CaTchan = Channel.Create([CaT_m0h0, CaT_m0h1, CaT_m1h0, CaT_m1h1, CaT_m2h0, CaT_m2h1])
with ssys, mdl, CaTchan[...]:
    
    
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    CaT_m0h0.s <r['CaTm0h0_m1h0']> CaT_m1h0.s ; r['CaTm0h0_m1h0'].setRates(VDepRate(lambda V: 1.0e3 *2.* alpham_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *1.* betam_cat(V*1.0e3)))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    
    CaT_m1h0.s <r['CaTm1h0_m2h0']> CaT_m2h0.s ; r['CaTm1h0_m2h0'].setRates(VDepRate(lambda V: 1.0e3 *1.* alpham_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *2.* betam_cat(V*1.0e3)))
    
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    CaT_m0h1.s <r['CaTm0h1_m1h1']> CaT_m1h1.s ; r['CaTm0h1_m1h1'].setRates(VDepRate(lambda V: 1.0e3 *2.* alpham_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *1.* betam_cat(V*1.0e3)))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    
    CaT_m1h1.s <r['CaTm1h1_m2h1']> CaT_m2h1.s ; r['CaTm1h1_m2h1'].setRates(VDepRate(lambda V: 1.0e3 *1.* alpham_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *2.* betam_cat(V*1.0e3)))
    
    
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    CaT_m0h0.s <r['CaTm0h0_m0h1']> CaT_m0h1.s ; r['CaTm0h0_m0h1'].setRates(VDepRate(lambda V: 1.0e3 *1.* alphah_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *1.* betah_cat(V*1.0e3)))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    CaT_m1h0.s <r['CaTm1h0_m1h1']> CaT_m1h1.s ; r['CaTm1h0_m1h1'].setRates(VDepRate(lambda V: 1.0e3 *1.* alphah_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *1.* betah_cat(V*1.0e3)))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    
    CaT_m2h0.s <r['CaTm2h0_m2h1']> CaT_m2h1.s ; r['CaTm2h0_m2h1'].setRates(VDepRate(lambda V: 1.0e3 *1.* alphah_cat(V*1.0e3)), VDepRate(lambda V: 1.0e3 *1.* betah_cat(V*1.0e3)))

if cyl160:
    OC_CaT = smodel.GHKcurr('OC_CaT', ssys, CaT_m2h1, Ca, virtual_oconc = Ca_oconc, computeflux = True)
else:
    OC_CaT = GHKCurr.Create(CaTchan[CaT_m2h1], Ca, CaT_P, computeflux=True)
with mdl:
    
    BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O0, BK_O1, BK_O2, BK_O3, BK_O4 = SubUnitState.Create()
    
    
    ##### BK channel ################################
    
    
    BKchan = Channel.Create([BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O0, BK_O1, BK_O2, BK_O3, BK_O4])
with ssys, mdl, BKchan[...]:
    
    
    
    Ca.i + BK_C0.s <r['BKCAC0']> BK_C1.s ; r['BKCAC0'].setRates(c_01, c_10)
    Ca.i + BK_C1.s <r['BKCAC1']> BK_C2.s ; r['BKCAC1'].setRates(c_12, c_21)
    Ca.i + BK_C2.s <r['BKCAC2']> BK_C3.s ; r['BKCAC2'].setRates(c_23, c_32)
    Ca.i + BK_C3.s <r['BKCAC3']> BK_C4.s ; r['BKCAC3'].setRates(c_34, c_43)
    
    
    Ca.i + BK_O0.s <r['BKCAO0']> BK_O1.s ; r['BKCAO0'].setRates(o_01, o_10)
    Ca.i + BK_O1.s <r['BKCAO1']> BK_O2.s ; r['BKCAO1'].setRates(o_12, o_21)
    Ca.i + BK_O2.s <r['BKCAO2']> BK_O3.s ; r['BKCAO2'].setRates(o_23, o_32)
    Ca.i + BK_O3.s <r['BKCAO3']> BK_O4.s ; r['BKCAO3'].setRates(o_34, o_43)
    
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    
    BK_C0.s <r['BKC0O0']> BK_O0.s ; r['BKC0O0'].setRates(VDepRate(lambda V: f_0(V)), VDepRate(lambda V: b_0(V)))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    BK_C1.s <r['BKC1O1']> BK_O1.s ; r['BKC1O1'].setRates(VDepRate(lambda V: f_1(V)), VDepRate(lambda V: b_1(V)))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    BK_C2.s <r['BKC2O2']> BK_O2.s ; r['BKC2O2'].setRates(VDepRate(lambda V: f_2(V)), VDepRate(lambda V: b_2(V)))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    BK_C3.s <r['BKC3O3']> BK_O3.s ; r['BKC3O3'].setRates(VDepRate(lambda V: f_3(V)), VDepRate(lambda V: b_3(V)))
    # WARNING: Using a variable name that is reserved (['V', 'V']).
    BK_C4.s <r['BKC4O4']> BK_O4.s ; r['BKC4O4'].setRates(VDepRate(lambda V: f_4(V)), VDepRate(lambda V: b_4(V)))
with ssys:
    
    OC_BK0 = OhmicCurr.Create(BKchan[BK_O0], BK_G, BK_rev)
    OC_BK1 = OhmicCurr.Create(BKchan[BK_O1], BK_G, BK_rev)
    OC_BK2 = OhmicCurr.Create(BKchan[BK_O2], BK_G, BK_rev)
    OC_BK3 = OhmicCurr.Create(BKchan[BK_O3], BK_G, BK_rev)
    OC_BK4 = OhmicCurr.Create(BKchan[BK_O4], BK_G, BK_rev)
with mdl:
    
    SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2 = SubUnitState.Create()
    
    
    ###### SK channel ##################
    
    
    SKchan = Channel.Create([SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2])
with ssys, mdl, SKchan[...]:
    
    
    
    Ca.i + SK_C1.s <r['SKCAC1']> SK_C2.s ; r['SKCAC1'].setRates(dirc2_t, invc1_t)
    Ca.i + SK_C2.s <r['SKCAC2']> SK_C3.s ; r['SKCAC2'].setRates(dirc3_t, invc2_t)
    Ca.i + SK_C3.s <r['SKCAC3']> SK_C4.s ; r['SKCAC3'].setRates(dirc4_t, invc3_t)
    
    
    SK_C3.s <r['SKC3O1']> SK_O1.s ; r['SKC3O1'].setRates(diro1_t, invo1_t)
    SK_C4.s <r['SKC4O2']> SK_O2.s ; r['SKC4O2'].setRates(diro2_t, invo2_t)
with ssys:
    
    OC1_SK = OhmicCurr.Create(SKchan[SK_O1], SK_G, SK_rev)
    OC2_SK = OhmicCurr.Create(SKchan[SK_O2], SK_G, SK_rev)
with mdl:
    Leak = SubUnitState.Create()
    
    
    ###### Leak current channel #####
    
    L = Channel.Create([Leak])
with ssys:
    
    OC_L = OhmicCurr.Create(L[Leak], L_G, L_rev)

##################################

########### MESH & COMPARTMENTALIZATION #################

##########Import Mesh

mesh = TetMesh.Load('./meshes/'+meshfile_ab)

outer_tets = range(len(mesh.tets))

###USE OF gettets
#getcyl(tetmesh, rad,  zmin, zmax, binnum=120, x = 0.0, y = 0.0):

inner_tets = gettets.getcyl(mesh, 1e-6, -200e-6, 200e-6)[0]

for i in inner_tets: outer_tets.remove(i)
assert(outer_tets.__len__() + inner_tets.__len__() == len(mesh.tets))

print(outer_tets.__len__(), " tets in outer compartment")
print(inner_tets.__len__(), " tets in inner compartment")


# Record voltage from the central tetrahedron
# WARNING: findTetByPoint was replaced by the new TetList syntax and trying to access a point outside the mesh will now raise an exception.
cent_tet = mesh.tets[0.0, 0.0, 0.0].idx
with mesh:
    
    ########## Create an intracellular compartment i.e. cytosolic compartment
    
    cyto = TetComp.Create(inner_tets, 'vsys')
    
    if not cyl160:
        outer = TetComp.Create(outer_tets)
    
    if cyl160:
        # Ensure that we use points a small distance inside the boundary:
        LENGTH = mesh.bbox.max[2] - mesh.bbox.min[2]
        boundminz = mesh.bbox.min[2] + LENGTH/len(mesh.tets)
        boundmaxz = mesh.bbox.max[2] - LENGTH/len(mesh.tets)

        memb_tris = list(mesh.surface.indices)
        minztris = []
        maxztris = []
        for tri in memb_tris:
            zminboundtri = True
            zmaxboundtri = True
            tritemp = mesh.tris[tri].verts.indices
            trizs = [0.0, 0.0, 0.0]
            trizs[0] = mesh.verts[tritemp[0]][2]
            trizs[1] = mesh.verts[tritemp[1]][2]
            trizs[2] = mesh.verts[tritemp[2]][2]
            for j in range(3):
                if (trizs[j]>boundminz):
                    zminboundtri = False
            if (zminboundtri):
                minztris.append(tri)
                continue
            for j in range(3):
                if (trizs[j]< boundmaxz):
                    zmaxboundtri = False
            if (zmaxboundtri):
                maxztris.append(tri)

        for t in minztris: memb_tris.remove(t)
        for t in maxztris: memb_tris.remove(t)
    
    else:
        print('Finding connecting triangles...')
        out_tris = set()
        for i in outer_tets:
                tritemp = mesh.tets[i].faces.indices
                for j in range(4): out_tris.add(tritemp[j])

        in_tris = set()
        for i in inner_tets:
                tritemp = mesh.tets[i].faces.indices
                for j in range(4): in_tris.add(tritemp[j])

        memb_tris = out_tris.intersection(in_tris)
        memb_tris = list(memb_tris)
    
    print(len(memb_tris), " surface triangles.")
    
    ########## Find the submembrane tets
    
    memb_tet_neighb = []
    for i in memb_tris:
        tettemp = mesh.tris[i].tetNeighbs.indices
        for j in tettemp:
            memb_tet_neighb.append(j)

    submemb_tets = []
    for i in memb_tet_neighb:
        if i in inner_tets:
            submemb_tets.append(i)

    print(len(submemb_tets))
    
    vol = 0.0
    
    for i in submemb_tets:
        vol = vol + mesh.tets[i].Vol
    
    print('Volume of submembrane region is', vol)
    
    
    ########## Create a membrane as a surface mesh
    if cyl160: 
        memb = sgeom.TmPatch('memb', mesh, memb_tris, cyto)
    else:
        memb = TetPatch.Create(memb_tris, cyto, outer, 'ssys')
    
    print("Area: ", memb.getArea())
    
    print("Creating membrane..")
    membrane = Membrane.Create([memb], opt_file_name='./meshes/'+meshfile_ab+"_optimalidx")
print("Membrane created.")

# # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

# WARNING: Using a variable name that is reserved (['r']).
r = RNG('mt19937', 512, 7)

# WARNING: Using a variable name that is reserved (['r']).
sim = Simulation('Tetexact', mdl, mesh, r, calcMembPot=True)

sim.newRun()

print("Injecting molecules..")

sim.Temp = TEMPERATURE+273.15

if not cyl160: 
    sim.outer.Ca.Conc = Ca_oconc
    sim.outer.Ca.Clamped = True
    
sim.cyto.Ca.Conc = Ca_iconc

print("Calcium concentration is: ", sim.cyto.Ca.Conc)
print("No. of Ca molecules is: ", sim.cyto.Ca.Count)
sim.cyto.Mg.Conc = Mg_conc

surfarea = sim.memb.Area

#Total pump is 1e-15 mol/cm2 ---> 1e-11 mol/m2
#pumpnbs per unit area (im m2) is Total pump times AVOGADRO's NUMBER (1e-11 mol/m2 * 6.022e23 /mol )
pumpnbs = 6.022141e12*surfarea
print("Number of pump molecules: ", pumpnbs)

sim.memb.Pump.Count = round(pumpnbs)
sim.memb.CaPump.Count = 0

sim.cyto.iCBsf.Conc = iCBsf_conc
sim.cyto.iCBCaf.Conc = iCBCaf_conc
sim.cyto.iCBsCa.Conc = iCBsCa_conc
sim.cyto.iCBCaCa.Conc = iCBCaCa_conc

sim.cyto.CBsf.Conc = CBsf_conc
sim.cyto.CBCaf.Conc = CBCaf_conc
sim.cyto.CBsCa.Conc = CBsCa_conc
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


#### Recording #####

# WARNING: Using a variable name that is reserved (['time']).
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


# WARNING: Using a variable name that is reserved (['r']).
r.initialize(100*int(iter_n))

for l in range(NTIMEPOINTS):
    print("Tpnt: ", l)
    
    # WARNING: Using a variable name that is reserved (['run']).
    sim.run(TIMECONVERTER*l)

    tcur_CaP = 0.0
    tcur_CaT = 0.0
    tcur_BK = 0.0
    tcur_SK = 0.0
    tca_count = 0.0
    for m in submemb_tets:
        tca_count = tca_count + sim.TET(m).Ca.Count
    for m in memb_tris:
        tcur_CaP = tcur_CaP + sim.TRI(m).OC_CaP.I 
        tcur_CaT = tcur_CaT + sim.TRI(m).OC_CaT.I
        tcur_BK = tcur_BK + sim.TRI(m).OC_BK0.I \
            + sim.TRI(m).OC_BK1.I \
            + sim.TRI(m).OC_BK2.I \
            + sim.TRI(m).OC_BK3.I \
            + sim.TRI(m).OC_BK4.I
        tcur_SK = tcur_SK + sim.TRI(m).OC1_SK.I + sim.TRI(m).OC2_SK.I
    
    datfile.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile.write('%.6g' %((tcur_CaP*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_CaT*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_BK*1.0e-1)/surfarea) + ' ')
    datfile.write('%.6g' %((tcur_SK*1.0e-1)/surfarea) + ' ') 
    datfile.write('\n')

    datfile2.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile2.write('%.6g' %(sim.TET(cent_tet).V*1.0e3) + ' ')
    datfile2.write('\n')
         
    datfile3.write('%.6g' %(1.0e3*TIMECONVERTER*l) + ' ')
    datfile3.write('%.6g' %(sim.cyto.Ca.Conc*1.0e6) + ' ')
    datfile3.write('%.6g' %tca_count + ' ')
    datfile3.write('\n')
    

datfile.close()
datfile2.close()
datfile3.close()
## END
