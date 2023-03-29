import steps.interface

from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *
from steps.utils import *

import numpy as np

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###########################################################
# Simulation Parameters
###########################################################

ENDT = 0.06

DT =  2.0e-5

EF_DT = 5.0e-6

SEED = 1234

###########################################################
# Model Parameters
###########################################################

TEMPERATURE = 34.0 + 273.15
Q10 = 3

FARADAY = 96485.3365     # C/mol
R = 8.3144621            # J/mol K
AVOGADRO = 6.02214129e23 # /mol

Qt = Q10 ** ((TEMPERATURE - (23 + 273.15)) / 10)
Qt_mslo = Q10 ** ((TEMPERATURE - (25 + 273.15))/10)

#######################################
# Membrane Parameters
#######################################

init_pot = Parameter(-60, 'mV', Description='Initial membrane potential')
Ra = Parameter(235.7*1.0e-2, 'ohm m', Description='Bulk resistivity')
Cm = 0.64e-2
memb_capac_proximal = Parameter(1.2*Cm, 'F m^-2', Description='Smooth membrane capacitance')
memb_capac_spiny = Parameter(5.3*Cm, 'F m^-2', Description='Spiny membrane capacitance')

#######################################
# CaP channels parameters
#######################################

CaP_P = Parameter(2.5e-2, 'um^3 s^-1', Description='CaP single channel permeability')
CaP_ro = Parameter(38, 'um^-2', Description='CaP channels density')

# Reaction rates

vhalfm = -29.458 # mV
cvm = 8.429      # mV

def minf_cap(mV):
    vhalfm = -29.458
    cvm = 8.429
    return 1 / (1 + np.exp(-(mV - vhalfm) / cvm))

def tau_cap(mV):
    if mV >= -40:
        return 0.2702 + 1.1622 * np.exp(-(mV + 26.798) ** 2 / 164.19)
    else:
        return 0.6923 * np.exp(mV / 1089.372)

alpha_cap = VDepRate.Create(
    lambda V: (minf_cap(V * 1e3) / tau_cap(V * 1e3)) * Qt * 1e3
)
beta_cap = VDepRate.Create(
    lambda V: (1 - minf_cap(V * 1e3)) / tau_cap(V * 1e3) * Qt * 1e3
)

# Initial conditions
CaP_p = [0.92402, 0.073988, 0.0019748, 1.7569e-05]

#######################################
# CaT channels parameters
#######################################

CaT_P = Parameter(1.65e-2, 'um^3 s^-1', Description='CaT single channel permeability')
CaT_ro = Parameter(1.9636, 'um^-2', Description='CaT channels density')

# Reaction rates

def minf_cat(mV):
    vhalfm = -52
    cvm = -5
    return 1 / (1 + np.exp((mV - vhalfm) / cvm))

def taum_cat(mV):
    if mV > -90:
        return 1 + 1 / (np.exp((mV + 40) / 9) + np.exp(-(mV + 102) / 18))
    else:
        return 1

def hinf_cat(mV):
    vhalfh = -72
    cvh = 7
    return 1 / (1 + np.exp((mV - vhalfh) / cvh))

def tauh_cat(mV):
    return (15 + 1 / (np.exp((mV + 32) / 7)))

alpham_cat = VDepRate.Create(lambda V: minf_cat(V * 1e3) / taum_cat(V * 1e3) * 1e3)
betam_cat = VDepRate.Create(lambda V: (1 - minf_cat(V * 1e3)) / taum_cat(V * 1e3) * 1e3)

alphah_cat = VDepRate.Create(lambda V: hinf_cat(V * 1e3) / tauh_cat(V * 1e3) * 1e3)
betah_cat = VDepRate.Create(lambda V: (1 - hinf_cat(V * 1e3)) / tauh_cat(V * 1e3) * 1e3)

# Initial conditions
CaT_p = [
    [0.58661, 0.23687, 0.023912],   # h0
    [0.10564, 0.042658, 0.0043063], # h1
]

#######################################
# BK channels parameters
#######################################

BK_G = Parameter(210, 'pS', Description='BK single channel conductance')
BK_ro = Parameter(1.5 * 2.0238, 'um^-2', Description='BK channels density')
BK_rev = Parameter(-77, 'mV', Description='BK channel reversal potential')

# Reaction rates

#Units (1)
Qo = 0.73
Qc = -0.67

#Units (/s)
pf0 = 2.39
pf1 = 5.4918
pf2 = 24.6205
pf3 = 142.4546
pf4 = 211.0220

pb0 = 3936
pb1 = 687.3251
pb2 = 234.5875
pb3 = 103.2204
pb4 = 11.6581

#Units(/M)
k1 = 1.0e6

#Units(/s)
onoffrate = 1.0e3

L0 = 1806

#Units (M)
Kc = 8.63e-6
Ko = 0.6563e-6

BK_f = k1*onoffrate*Qt_mslo
BKo_b = Ko*k1*onoffrate*Qt_mslo
BKc_b = Kc*k1*onoffrate*Qt_mslo

BK_f0 = VDepRate.Create(
    lambda V: pf0 * Qt_mslo * (np.exp((Qo * FARADAY * V) / (R * TEMPERATURE)))
)
BK_f1 = VDepRate.Create(
    lambda V: pf1 * Qt_mslo * (np.exp((Qo * FARADAY * V) / (R * TEMPERATURE)))
)
BK_f2 = VDepRate.Create(
    lambda V: pf2 * Qt_mslo * (np.exp((Qo * FARADAY * V) / (R * TEMPERATURE)))
)
BK_f3 = VDepRate.Create(
    lambda V: pf3 * Qt_mslo * (np.exp((Qo * FARADAY * V) / (R * TEMPERATURE)))
)
BK_f4 = VDepRate.Create(
    lambda V: pf4 * Qt_mslo * (np.exp((Qo * FARADAY * V) / (R * TEMPERATURE)))
)
BK_oc_f = [BK_f0, BK_f1, BK_f2, BK_f3, BK_f4]

BK_b0 = VDepRate.Create(
    lambda V: pb0 * Qt_mslo * (np.exp((Qc * FARADAY * V) / (R * TEMPERATURE)))
)
BK_b1 = VDepRate.Create(
    lambda V: pb1 * Qt_mslo * (np.exp((Qc * FARADAY * V) / (R * TEMPERATURE)))
)
BK_b2 = VDepRate.Create(
    lambda V: pb2 * Qt_mslo * (np.exp((Qc * FARADAY * V) / (R * TEMPERATURE)))
)
BK_b3 = VDepRate.Create(
    lambda V: pb3 * Qt_mslo * (np.exp((Qc * FARADAY * V) / (R * TEMPERATURE)))
)
BK_b4 = VDepRate.Create(
    lambda V: pb4 * Qt_mslo * (np.exp((Qc * FARADAY * V) / (R * TEMPERATURE)))
)
BK_oc_b = [BK_b0, BK_b1, BK_b2, BK_b3, BK_b4]

# Initial conditions
BK_p = [
    [0.99997,    4.3619e-07, 4.1713e-09, 4.4449e-11, 6.3132e-14],
    [2.5202e-05, 1.1765e-06, 6.6148e-08, 2.4392e-09, 4.0981e-11],
]

#######################################
# SK channels parameters
#######################################

SK_G = Parameter(10, 'pS', Description='SK single channel conductance')
SK_ro = Parameter(0.31, 'um^-2', Description='SK channels density')
SK_rev = Parameter(-77, 'mV', Description='SK channel reversal potential')

# Reaction rates

#Units (/s)
invc1 = 80
invc2 = 80
invc3 = 200

invo1 = 1000
invo2 = 100

diro1 = 160
diro2 = 1200

#Units ( /s M)

dirc2 = 200e6
dirc3 = 160e6
dirc4 = 80e6

invc1_t = invc1*Qt
invc2_t = invc2*Qt
invc3_t = invc3*Qt

invo1_t = invo1*Qt
invo2_t = invo2*Qt

diro1_t = diro1*Qt
diro2_t = diro2*Qt

dirc2_t = dirc2*Qt/3.0
dirc3_t = dirc3*Qt/3.0
dirc4_t = dirc4*Qt/3.0

# Intital conditions
SK_C1_p= 0.96256
SK_C2_p= 0.036096
SK_C3_p= 0.0010829
SK_C4_p= 6.4973e-06

SK_O1_p= 0.00017326
SK_O2_p= 7.7967e-05

#######################################
# AMPA channels parameters
#######################################

AMPA_G = Parameter(7, 'pS', Description="AMPA channel conductance")
AMPA_Ro = Parameter(20, 'um^-2', Description="AMPA channek density")
AMPA_rev = Parameter(0, 'mV', Description="AMPA channel reversal potential")

#Units (/s M)

rb = 13e6

#Units (/s)

ru1 = 0.0059e3
ru2 = 86e3
ro = 2.7e3
rc = 0.2e3
rd = 0.9e3
rr = 0.064e3

# Glutamate transient

# Units (s)

glut_start = 5e-3

# Units (M), Reference (Rudolph et al. 2011)
def Glut(t):
    if t >= 0:
        return 12e-3 * (0.86 * np.exp(-t / 0.4e-3) + 0.14 * np.exp(-t / 4.2e-3))
    else:
        return 0

#######################################
# Leak channels parameters
#######################################

L_G = Parameter(0.4, 'pS', Description='Leak single channel conductance')
L_ro_proximal = Parameter(0.25, 'um^-2', Description='Proximal leak channels density')
L_ro_spiny = Parameter(1, 'um^-2', Description='Proximal leak channels density')
L_rev = Parameter(-61, 'mV', Description='Leak channel reversal potential')

#######################################
# Ca pump channels parameters
#######################################

P_ro = Parameter(6.022141, 'um^-2', Description='Ca2+ pump density')

# Reaction rates

P_f = 3e9
P_b = 1.75e4
P_k = 7.255e4

#######################################
# Calcium buffering parameters
#######################################

# Ca concentrations

Ca_oconc = 2e-3
Ca_iconc = 45e-9

# Mg concentrations

Mg_conc = 590e-6

# Buffer concentrations

iCBsf_conc = 27.704e-6
iCBCaf_conc = 2.6372e-6
iCBsCa_conc= 1.5148e-6
iCBCaCa_conc= 0.14420e-6

CBsf_conc= 110.82e-6
CBCaf_conc= 10.549e-6
CBsCa_conc= 6.0595e-6
CBCaCa_conc= 0.57682e-6

PV_conc= 3.2066e-6
PVCa_conc= 16.252e-6
PVMg_conc= 60.541e-6

# Diffusion constants

DCST = 0.223e-9 # Ca
DCB = 0.028e-9  # Calbindin (CB)
DPV = 0.043e-9  # Parvalbumin (PV)

# Reaction rates

CBf_f_kcst = 4.35e7
CBf_b_kcst = 35.8

CBs_f_kcst = 0.55e7
CBs_b_kcst = 2.6

PVca_f = 10.7e7
PVca_b = 0.95

PVmg_f_kcst = 0.8e6
PVmg_b_kcst = 25

###########################################################
# Mesh Parameters
###########################################################

mesh_file = '../meshes/caburst_dist/split_256/CNG_segmented_2_split_256'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###########################################################
# Biochemical model
###########################################################

mdl = Model()
with mdl:
    # Species
    Pump, CaPump, PV, PVMg, PVCa, Mg = Species.Create()
    Ca = Species.Create(valence=2)

    # Calbindin
    CBs, CBf, CBsCa, CBfCa, CBmob, CBimmob = SubUnitState.Create()
    CBsSU, CBfSU, CBmobSU = SubUnit.Create(
        [CBs, CBsCa], [CBf, CBfCa], [CBmob, CBimmob]
    )
    CB = Complex.Create([CBsSU, CBfSU, CBmobSU])

    # Channels
    CaPc, CaPo = SubUnitState.Create()
    CaP_SU = SubUnit.Create([CaPc, CaPo])
    CaPchan = Channel.Create([CaP_SU, CaP_SU, CaP_SU])

    BK, BKCa, BKopen, BKclose = SubUnitState.Create()
    BKCaSU = SubUnit.Create([BK, BKCa])
    BKocSU = SubUnit.Create([BKopen, BKclose])
    BKchan = Channel.Create([BKCaSU, BKCaSU, BKCaSU, BKCaSU, BKocSU])

    SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2 = SubUnitState.Create()
    SKchan = Channel.Create([SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2])

    AMPA_C, AMPA_C1, AMPA_C2, AMPA_D1, AMPA_D2, AMPA_O = SubUnitState.Create()
    AMPA = Channel.Create([AMPA_C, AMPA_C1, AMPA_C2, AMPA_D1, AMPA_D2, AMPA_O])

    Leak = SubUnitState.Create()
    L = Channel.Create([Leak])
    
    r = ReactionManager()

    vsys = VolumeSystem.Create()
    with vsys:
        # PVCa
        PV + Ca <r[1]> PVCa
        r[1].K = PVca_f, PVca_b

        # PVMg
        PV + Mg <r[1]> PVMg
        r[1].K = PVmg_f_kcst, PVmg_b_kcst

        with CB[...]:
            # Fast binding
            CBf + Ca <r[1]> CBfCa
            r[1].K = CBf_f_kcst, CBf_b_kcst

            # Slow binding
            CBs + Ca <r[1]> CBsCa
            r[1].K = CBs_f_kcst, CBs_b_kcst

        diff_Ca   = Diffusion.Create(Ca, DCST)
        diff_CB   = Diffusion.Create(CB[:, :, CBmob], DCB)
        diff_PV   = Diffusion.Create(PV, DPV)
        diff_PVCa = Diffusion.Create(PVCa, DPV)
        diff_PVMg = Diffusion.Create(PVMg, DPV)

    ssys = SurfaceSystem.Create()
    with ssys:
        # Ca Pump
        Pump.s + Ca.i <r[1]> CaPump.s >r[2]> Pump.s
        r[1].K = P_f, P_b
        r[2].K = P_k
    
        # CaP channel
        with CaPchan[...]:
            CaPc.s <r[1]> CaPo.s
            r[1].K = VDepRate(alpha_cap), VDepRate(beta_cap)
        OC_CaP = GHKCurr.Create(
            CaPchan[CaPo, CaPo, CaPo], Ca, CaP_P,
            computeflux=True, virtual_oconc=Ca_oconc
        )
    
        # BK channel
        with BKchan[..., BKclose]:
            BK.s + Ca.i <r[1]> BKCa.s
            r[1].K = BK_f, BKc_b

        with BKchan[..., BKopen]:
            BK.s + Ca.i <r[1]> BKCa.s
            r[1].K = BK_f, BKo_b

        with BKchan[...]:
            BKclose.s <r[1]> BKopen.s
            r[1].K = (
                CompDepRate(lambda s: BK_oc_f[s.Count(BKCa)], [BKchan]),
                CompDepRate(lambda s: BK_oc_b[s.Count(BKCa)], [BKchan])
            )
        OC_BK = OhmicCurr.Create(BKchan[..., BKopen], BK_G, BK_rev)

        # SK channel
        with SKchan[...]:
            ((SK_C1.s + Ca.i <r[1]> SK_C2.s)\
                      + Ca.i <r[2]> SK_C3.s)\
                      + Ca.i <r[3]> SK_C4.s
            r[1].K = dirc2_t, invc1_t
            r[2].K = dirc3_t, invc2_t
            r[3].K = dirc4_t, invc3_t
            
            SK_C3.s <r[1]> SK_O1.s
            SK_C4.s <r[2]> SK_O2.s
            r[1].K = diro1_t, invo1_t
            r[2].K = diro2_t, invo2_t
        OC_SK = OhmicCurr.Create(SKchan[SK_O1|SK_O2], SK_G, SK_rev)

        # AMPA channel
        with AMPA[...]:
            AMPA_C.s <r['A1']> AMPA_C1.s <r['A2']> AMPA_C2.s <r[3]> AMPA_O.s
            r['A1'].K = 0, ru1
            r['A2'].K = 0, ru2
            r[3].K = ro, rc

            AMPA_C1.s <r[1]> AMPA_D1.s
            AMPA_C2.s <r[2]> AMPA_D2.s
            r[1].K = rd, rr
            r[2].K = rd, rr
        OC_AMPAR1 = OhmicCurr.Create(AMPA[AMPA_O], AMPA_G, AMPA_rev)

        # Leak current channel
        OC_L = OhmicCurr.Create(L[Leak], L_G, L_rev)
    
###########################################################
# Mesh and compartmentalization
###########################################################

mesh = DistMesh(mesh_file, 1e-6)

with mesh.asLocal(owned=False):
    cyto = Compartment.Create(mesh.tets, vsys)
    cyto.Conductivity = 1 / Ra

    smoothTets = mesh.tetGroups['smooth']
    spinyTets = mesh.tetGroups['spiny']

    smooth = Patch.Create(smoothTets.surface & mesh.surface, cyto, None, ssys)
    spiny = Patch.Create(spinyTets.surface & mesh.surface, cyto, None, ssys)

    memb_smooth = Membrane.Create([smooth], capacitance=memb_capac_proximal)
    memb_spiny = Membrane.Create([spiny], capacitance=memb_capac_spiny)

###########################################################
# Simulation
###########################################################

rng = RNG('mt19937', 512, SEED)

sim = Simulation(
    'DistTetOpSplit', mdl, mesh, rng,
    searchMethod=NextEventSearchMethod.GIBSON_BRUCK
)

rs = ResultSelector(sim)

with mesh.asLocal():
    Pots = rs.VERTS(mesh.surface.verts).V
    CaConc = rs.TETS().Ca.Conc
    CaConc.metaData['Vols'] = [tet.Vol for tet in mesh.tets]
    Currents = rs.TRIS(mesh.surface).OC_CaP.I

sim.toSave(Pots, CaConc, Currents, dt=DT)

with XDMFHandler('DistCaburst') as hdf:
    sim.toDB(hdf, 'CaBurstSim')

    sim.newRun()

    # Set temperature for ghk reactions
    sim.Temp = TEMPERATURE + 273.15

    # Setting initial conditions
    simPatches = [sim.smooth, sim.spiny]
    smoothArea = Parameter(smooth.Area, 'm^2')
    spinyArea = Parameter(spiny.Area, 'm^2')
    for patch, area in zip(simPatches, [smoothArea, spinyArea]):
        patch.Pump.Count = round(P_ro * area)

        for s in CaPchan[...]:
            patch.LIST(s).Count = round(CaP_ro*area*CaP_p[s.Count(CaPo)])

        for s in BKchan[...]:
            isOpen, nbCa = s.Count(BKopen), s.Count(BKCa)
            patch.LIST(s).Count = round(BK_ro*area*BK_p[isOpen][nbCa])

        patch.SKchan[SK_C1].Count = round(SK_ro*area*SK_C1_p)
        patch.SKchan[SK_C2].Count = round(SK_ro*area*SK_C2_p)
        patch.SKchan[SK_C3].Count = round(SK_ro*area*SK_C3_p)
        patch.SKchan[SK_C4].Count = round(SK_ro*area*SK_C4_p)
        patch.SKchan[SK_O1].Count = round(SK_ro*area*SK_O1_p)
        patch.SKchan[SK_O2].Count = round(SK_ro*area*SK_O2_p)

    sim.spiny.L[Leak].Count = round(L_ro_spiny * spinyArea)
    sim.smooth.L[Leak].Count = round(L_ro_proximal * smoothArea)

    sim.smooth.AMPA[AMPA_C].Count = round(AMPA_Ro * smoothArea)

    sim.cyto.Ca.Conc = Ca_iconc
    sim.cyto.Mg.Conc = Mg_conc

    sim.cyto.CB[CBs,   CBf,   CBimmob].Conc = iCBsf_conc
    sim.cyto.CB[CBsCa, CBf,   CBimmob].Conc = iCBCaf_conc
    sim.cyto.CB[CBs,   CBfCa, CBimmob].Conc = iCBsCa_conc
    sim.cyto.CB[CBsCa, CBfCa, CBimmob].Conc = iCBCaCa_conc

    sim.cyto.CB[CBs,   CBf,   CBmob].Conc = CBsf_conc
    sim.cyto.CB[CBsCa, CBf,   CBmob].Conc = CBCaf_conc
    sim.cyto.CB[CBs,   CBfCa, CBmob].Conc = CBsCa_conc
    sim.cyto.CB[CBsCa, CBfCa, CBmob].Conc = CBCaCa_conc

    sim.cyto.PV.Conc = PV_conc
    sim.cyto.PVCa.Conc = PVCa_conc
    sim.cyto.PVMg.Conc = PVMg_conc

    sim.EfieldDT = EF_DT

    sim.ALL(Membrane).Potential = init_pot

    for t in np.arange(0, ENDT, DT):
        if MPI.rank == 0:
            print(f'Run timestep {t}')
        # Update AMPA receptors reaction rate
        sim.smooth.A1['fwd'].K = rb * Glut(t - glut_start)
        sim.smooth.A2['fwd'].K = rb * Glut(t - glut_start)

        sim.run(t)

