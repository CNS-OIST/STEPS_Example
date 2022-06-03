.. _stoch_spikes_API2:

*************************************************
Stochastic Calcium Burst model with GHK currents
*************************************************

*This chapter is based on the publication "Anwar H, Hepburn I, Nedelescu H, Chen W, De Schutter E (2013) Stochastic 
Calcium Mechanisms Cause Dendritic Calcium Spike Variability. The Journal of Neuroscience, 33(40): 15848-15867". The model 
file is StochasticCaburst.py which is available at http://senselab.med.yale.edu/modeldb/ShowModel.asp?model=150635 and should 
be run within the complete modelling package because it depends on other files such as /extra/constants.py.*


The simulation script described in this chapter is also available at `STEPS_Example repository <https://github.com/CNS-OIST/STEPS_Example/tree/master/publication_models/API_2/Anwar_J%20Neurosci_2013>`_.

This chapter builds on previous chapters by simulating a model that includes both a reaction-diffusion component as well as  
electrical excitability. As described in [#f1]_, the two are closely coupled and this model contains ion channels where
activation is both voltage-dependent and calcium-dependent. In addition, the calcium ions form an important part of the 
current across the membrane, further coupling the reaction-diffusion component with the electrical excitability. This 
chapter introduces an important new object: the GHK Current object, which is described in some detail in section :ref:`ptype`.

As in previous chapters we will go through the script, looking in some depth at new concepts, but only brief explanations
will be offered of things that have been described in previous chapters. 


Modelling solution
==================

At the top of the script, as usual, we import some modules, including STEPS modules and some self-written helper modules. This
import includes the extra/constants.py file, which includes all important parameters for the module such as physical 
constants, membrane properties, kinetic properties of the channels, initial conditions and so on::

    import steps.interface

    from steps.model import *
    from steps.geom import *
    from steps.rng import *
    from steps.sim import *
    from steps.saving import *

    from extra.constants import *

    import math
    import os
    import sys
    import time
    
.. _constants_API2:

extra/constants.py file
-----------------------

Although not part of this model script, for easy reference to many important parameters used in the model this is the contents of the imported 
extra.constants module. What these parameters are and how they are used will become apparent as we go through the model. All unit are S.I. units, 
with the exception of concentration where units are molar::

    import math

    # # # # # # # # # # # # # # # # SIMULATION CONTROLS # # # # # # # # # # # # #

    EF_DT = 2.0e-5          # The EField dt
    NTIMEPOINTS =  25001 
     
    TIMECONVERTER =  2.0e-5

    NITER = 1

    ############################ PARAMETERS ################################

    init_pot = -60e-3

    TEMPERATURE = 34.0

    Q10 = 3

    # Faraday constant: unit of FARADAY is C/mol 
    # Source: http://physics.nist.gov/cgi-bin/cuu/Value?f 24/2/2012
    FARADAY = 96485.3365

    # Molar Gas Constant: unit of R is J/mol K
    # Source: http://physics.nist.gov/cgi-bin/cuu/Value?r 24/2/2012
    R = 8.3144621

    # Avogadro constant: unit of AVOGADRO is /mol
    # Source: http://physics.nist.gov/cgi-bin/cuu/Value?na 24/2/2012
    AVOGADRO = 6.02214129e23

    # Elementary charge: unit of E_CHARGE is C
    # Source: http://physics.nist.gov/cgi-bin/cuu/Value?e 24/2/2012
    E_CHARGE = 1.602176565e-19


    #FOR MSLO, THERE IS A NEW VALUE FOR Qt wrt to 25 degC
    Qt = math.pow(Q10, ((TEMPERATURE-23)/10))
    Qt_mslo = math.pow(Q10, ((TEMPERATURE-25)/10))

    ########## BULK RESISTIVITY ##########

    Ra = 235.7*1.0e-2

    ########## MEMBRANE CAPACITANCE ##########

    memb_capac = 1.5e-2


    ########## CaP channels density & permiability per channel ##########

    # CaP_P is permiability per channel (m3/s)
    # CaP_ro is channel/surface area (/m2)
    # P in Ca Dynamics model is 0.95e-4 cm/s --> 0.95e-6 m/s

    CaP_P = 2.5e-20 
    CaP_ro = 3.8e13

    ##########CaP channel parameters ####################

    #Units (mV)
    vhalfm = -29.458
    cvm = 8.429

    def minf_cap(V):
        #Units (mV)
        vhalfm = -29.458
        cvm = 8.429
        vshift = 0.0
        
        return (1.0/(1.0 + math.exp(-(V-vhalfm-vshift)/cvm)))

    def tau_cap(V):
        vshift = 0.0
        if (V-vshift) >= -40:
            return (0.2702 + 1.1622 * math.exp(-(V+26.798-vshift)*(V+26.798-vshift)/164.19))
        else:
            return (0.6923 * math.exp((V-vshift)/1089.372))

    def alpha_cap(V):
        return (minf_cap(1e3 * V)/tau_cap(1e3 * V)) * Qt * 1e3

    def beta_cap(V):
        return ((1.0-minf_cap(V * 1e3))/tau_cap(V * 1e3)) * Qt * 1e3


    ## Intitial conditions

    CaP_m0_p = 0.92402
    CaP_m1_p = 0.073988
    CaP_m2_p = 0.0019748
    CaP_m3_p = 1.7569e-05


    ########## CaT channels density & permiability per channel ##########

    # CaT_P is permiability per channel (m3/s)
    # CaT_ro is channel/surface area (/m2)
    # P in Ca Dynamics model is 6.2e-6 cm/s -->6.2e-8 m/s

    CaT_P = 1.65e-20
    CaT_ro = 3.7576e12


    def minf_cat(V):
        #Units (mV)
        vhalfm = -52.0
        cvm = -5.0
        vshift = 0.0
        
        return (1.0/(1.0 + math.exp((V-vhalfm-vshift)/cvm)))

    def taum_cat(V):
        vshift = 0.0
        if V > -90.0:
            return (1.0 + 1.0 / (math.exp((V+40.0-vshift)/9.0) + math.exp(-(V+102.0-vshift)/18.0)))
        else:
            return 1.0

    def hinf_cat(V):
        vhalfh = -72.0
        cvh = 7.0
        vshift = 0.0
        return (1.0/(1.0 + math.exp((V-vhalfh-vshift)/cvh)))

    def tauh_cat(V):
        vshift = 0.0
        return (15.0 + 1.0 / (math.exp((V+32.0-vshift)/7.0)))

    def alpham_cat(V):
        return (minf_cat(1e3 * V)/taum_cat(1e3 * V)) * 1e3

    def betam_cat(V):
        return ((1-minf_cat(1e3 * V))/taum_cat(1e3 * V)) * 1e3

    def alphah_cat(V):
        return (hinf_cat(1e3 * V)/tauh_cat(1e3 * V)) * 1e3

    def betah_cat(V):
        return ((1-hinf_cat(1e3 * V))/tauh_cat(1e3 * V)) * 1e3

    ## Initial conditions

    CaT_m0h0_p = 0.58661
    CaT_m1h0_p = 0.23687
    CaT_m2h0_p = 0.023912
    CaT_m0h1_p = 0.10564
    CaT_m1h1_p = 0.042658
    CaT_m2h1_p = 0.0043063

    ########## BK channels density & conductance per channel ##########

    # Total conductance = BK_G (conductance/channel) * BK_ro (channel/surface area)
    # BK in Ca Dynamics model is 4.25e-2 S/cm2 --> 4.25e2 S/m2


    BK_G = 2.1e-10
    BK_ro = 2.0238e12
    BK_rev = -77e-3

    ######### BK channel parameters ######################

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


    c_01 = 4.*k1*onoffrate*Qt_mslo
    c_12 = 3.*k1*onoffrate*Qt_mslo
    c_23 = 2.*k1*onoffrate*Qt_mslo
    c_34 = 1.*k1*onoffrate*Qt_mslo
    o_01 = 4.*k1*onoffrate*Qt_mslo
    o_12 = 3.*k1*onoffrate*Qt_mslo
    o_23 = 2.*k1*onoffrate*Qt_mslo
    o_34 = 1.*k1*onoffrate*Qt_mslo

    c_10 = 1.*Kc*k1*onoffrate*Qt_mslo
    c_21 = 2.*Kc*k1*onoffrate*Qt_mslo
    c_32 = 3.*Kc*k1*onoffrate*Qt_mslo
    c_43 = 4.*Kc*k1*onoffrate*Qt_mslo
    o_10 = 1.*Ko*k1*onoffrate*Qt_mslo
    o_21 = 2.*Ko*k1*onoffrate*Qt_mslo
    o_32 = 3.*Ko*k1*onoffrate*Qt_mslo
    o_43 = 4.*Ko*k1*onoffrate*Qt_mslo


    f_0 = lambda mV: pf0*Qt_mslo*(math.exp((Qo* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
    f_1 = lambda mV: pf1*Qt_mslo*(math.exp((Qo* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
    f_2 = lambda mV: pf2*Qt_mslo*(math.exp((Qo* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
    f_3 = lambda mV: pf3*Qt_mslo*(math.exp((Qo* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
    f_4 = lambda mV: pf4*Qt_mslo*(math.exp((Qo* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))

    b_0 = lambda mV: pb0*Qt_mslo*(math.exp((Qc* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
    b_1 = lambda mV: pb1*Qt_mslo*(math.exp((Qc* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
    b_2 = lambda mV: pb2*Qt_mslo*(math.exp((Qc* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
    b_3 = lambda mV: pb3*Qt_mslo*(math.exp((Qc* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
    b_4 = lambda mV: pb4*Qt_mslo*(math.exp((Qc* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))


    # Initial conditions
    BK_C0_p= 0.99997
    BK_C1_p= 4.3619e-07
    BK_C2_p= 4.1713e-09
    BK_C3_p= 4.4449e-11
    BK_C4_p= 6.3132e-14

    BK_O0_p= 2.5202e-05
    BK_O1_p= 1.1765e-06
    BK_O2_p= 6.6148e-08
    BK_O3_p= 2.4392e-09
    BK_O4_p= 4.0981e-11

    ########## SK channel density & conductance per channel #############

    # Total conductance = SK_G (conductance/channel) * SK_ro (channel/surface area)
    # SK in Ca Dynamics model is 3.1e-4 S/cm2 --> 3.1 S/m2


    SK_G = 1.0e-11
    SK_ro = 31.0e10

    SK_rev = -77e-3

    ######### SK channel parameters ###################

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


    ######### leak current channel density & conductance per channel ########
    # Total conductance = 1e-6 S/cm2 --> 1e-2 S/m2

    L_G = 4.0e-14
    L_ro = 25.0e10

    L_rev = -61e-3


    ######### Pump parameters ###################

    P_f_kcst = 3e9
    P_b_kcst = 1.75e4
    P_k_kcst = 7.255e4


    ############################CALCIUM BUFFERING MODEL################################

    ########## Ca concentrations #########

    Ca_oconc = 2e-3
    Ca_iconc = 45e-9

    ########## Mg concentrations #########

    Mg_conc = 590e-6

    ########## Buffer concentrations #############

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

    # Diffusion constant of Calcium
    DCST = 0.223e-9
    # Diffusion constant of Calbindin (CB)
    DCB = 0.028e-9
    # Diffusion constant of Parvalbumin (PV)
    DPV = 0.043e-9

    #iCBsf-fast
    iCBsf1_f_kcst = 4.35e7
    iCBsf1_b_kcst = 35.8

    #iCBsCa
    iCBsCa_f_kcst = 0.55e7
    iCBsCa_b_kcst = 2.6

    #iCBsf_slow
    iCBsf2_f_kcst = 0.55e7
    iCBsf2_b_kcst = 2.6

    #iCBCaf
    iCBCaf_f_kcst = 4.35e7
    iCBCaf_b_kcst = 35.8

    #CBsf-fast
    CBsf1_f_kcst = 4.35e7
    CBsf1_b_kcst = 35.8

    #CBsCa
    CBsCa_f_kcst = 0.55e7
    CBsCa_b_kcst = 2.6

    #CBsf_slow
    CBsf2_f_kcst = 0.55e7
    CBsf2_b_kcst = 2.6

    #CBCaf
    CBCaf_f_kcst = 4.35e7
    CBCaf_b_kcst = 35.8

    #PVca
    PVca_f_kcst = 10.7e7
    PVca_b_kcst = 0.95

    #PVmg
    PVmg_f_kcst = 0.8e6
    PVmg_b_kcst = 25

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


Command line execution
----------------------


Back to the StochasticCaburst.py script, and next we take a slightly new approach to previous models, in which we utilise command line arguments (``sys.argv``).
The 'zeroth' argument (``sys.argv[0]``) is the script pathname, then other arguments are command line arguments which are read as 
strings. In this sense it is intended that the StochasticCaburst.py script is run from the command line with a statement such as::

    $ python StochasticCaburst.py Cylinder2_dia2um_L10um_outer0_3um_0.3shell_0.3size_19156tets_adaptive.inp ~/stochcasims/ 1

and, in the present form, this script can not be run interactively as previous examples can. In the script, we read command line arguments thusly::

    _, meshfile_ab, root, iter_n = sys.argv

Looking at the above example, sys.argv[0] contains the name of the script and is discarded (by attributing it to `_`), sys.argv[1] would be the string 'Cylinder2_dia2um_L10um_outer0_3um_0.3shell_0.3size_19156tets_adaptive.inp' 
and will be stored as variable ``meshfile_ab``, 
sys.argv[2] will be the string '~/stochcasims/' and stored as variable ``root`` which, as we will see, defines where to store
simulation output, and sys.argv[3] will be the string '1' and stored as variable ``iter_n`` which is involved in random 
number initialisation (after conversion to an integer) and data storage. We will look at these variables later in the script as they are used. 

The last thing to be done before moving onto the biochemical model description is to set a flag to detect whether we are 
dealing with the 160um mesh (available meshes that are intended for this simulation are cylinders all of diameter 2um and lengths of 10um, 
20um, 40um, 80um and 160um). This flag is necessary because the 160um mesh has slightly different properties 
than the other meshes in that it contains no extracellular tetrahedral compartment, and so there will be some different behaviours throughout the script 
depending on whether we are using the 160um mesh or not::

    cyl160 = meshfile_ab == 'Cylinder2_dia2um_L160um_outer0_0.3shell_0.3size_279152tets_adaptive.inp'

Model specification
===================

Since this is a relatively large model we will split its description up into two sections, after first creating the parent 
model container object, one volume system and one surface system::

    mdl = Model()
    r = ReactionManager()
    with mdl:
        
        # Vol/surface systems
        vsys = VolumeSystem.Create()
        ssys = SurfaceSystem.Create()

.. _calc_dyn_API2:

Calcium dynamics
--------------------

The lines of code we will highlight in this chapter do not necessarily appear in the same order in the full script, but changing the order, while not having any impact on the actual model, will help us understand how the code works. We will add a `# ...` line when lines from the original scripts are ommited for clarity.

The following lines of code describe the calcium and calcium buffer reactions and diffusion. Since these are 'ordinary' 
dynamics with no voltage-dependence we will not look look at this part in detail. A more detailed explanation is offered
in [#f1]_ and [#f2]_. Most parameters come from the :ref:`constants`::

    # Calcium
    Ca = Species.Create(valence=2)
    
    # Species
    Pump, CaPump, iCBsf, iCBsCa, iCBCaf, iCBCaCa, CBsf, CBsCa, CBCaf, CBCaCa, PV, PVMg, PVCa, Mg = Species.Create()

    # ...

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

.. _ptype_API2:

P-type Calcium channel
----------------------

The P-type calcium channel is a different type of ion channel to those we have seen before. In previous chapters we saw 
Hodgkin-Huxley sodium and potassium channels that conducted an Ohmic current. The sodium and potassium ions in that situation 
were not explicitly simulated, which was reasonable because those ions were not involved in other processes we were 
interested in, and we could assume their concentrations inside and outside the cell were not altered significantly during their 
conduction. However, with calcium we need a different approach. Here calcium is involved in intracellular processes such as 
potassium channel-activation (as we will see), buffering and diffusion, and so we must simulate the influx of calcium through these
P-type channels. Furthermore, the Ohmic approximation is no longer sufficient for our purposes. The large differences between 
intracellular and extracellular concentration along with large changes in intracellular concentration mean that, in effect, channel 
conductance has some voltage and concentration dependence and is described much better by the GHK flux equation. The GHK flux equation itself 
is derived under certain simplifying assumptions that are good approximations for many ion channels, specifically 
those where channel occupancy and competition are negligible. Please see [#f3]_ for further discussion on the use of the GHK flux 
equation and the behaviour of the GHK current object in STEPS. It is worth noting that use of the GHK flux equation means that 
(instead of conductance) we must specify the channel's permeability, which can be more difficult to parameterize. 

The P-type calcium channel kinetics are described in detail in [#f1]_. To create the channel we first describe the channel states::
 
    ###### CaP channel ############## 
    CaP_m0, CaP_m1, CaP_m2, CaP_m3 = SubUnitState.Create()
    CaPchan = Channel.Create([CaP_m0, CaP_m1, CaP_m2, CaP_m3])
    
and the voltage-dependent kinetics. Remember for each of these discrete channels this voltage will be read from the local voltage
across the membrane triangle where the channel resides::

    with ssys:

        # ...

        with CaPchan[...]:
            CaP_m0.s <r['CaPm0m1']> CaP_m1.s <r['CaPm1m2']> CaP_m2.s <r['CaPm2m3']> CaP_m3.s
            r['CaPm0m1'].K = 3 * VDepRate(alpha_cap), 1 * VDepRate(beta_cap)
            r['CaPm1m2'].K = 2 * VDepRate(alpha_cap), 2 * VDepRate(beta_cap)
            r['CaPm2m3'].K = 1 * VDepRate(alpha_cap), 3 * VDepRate(beta_cap)

We come to creating our GHK current object (:class:`steps.API_1.model.GHKcurr`). This object will calculate single-channel current for a given
channel state by the GHK flux equation:

.. math::
     I_{s}=P_{s}z_{s}^{2}\frac{V_{m}F^{2}}{RT}\frac{[S]_{i}-[S]_{o}exp(-z_{s}V_{m}F/RT)}{1-exp(-z_{s}V_{m}F/RT)}
    :label: 9.1.2

where :math:`I_{s}` is the single-channel current (amps) of ion S, :math:`P_{s}` is the single-channel permeability of ion S (:math:`m^{3}.s^{-1}`), :math:`z_{s}` is the valence of ion S, :math:`V_{m}` is the membrane voltage (volts), F is the Faraday constant, R is the gas constant, T is temperature (Kelvin), :math:`[S]_{i}` is the intracellular concentration of ion S (:math:`mol.m^{-3}`) and :math:`[S]_{o}` is the extracellular concentration of ion S (:math:`mol.m^{-3}`).

When a GHK current is applied in STEPS it (optionally) results in movement of ions between the 'outer' and 'inner' compartments, the direction of which will depend 
on the sign of the current and the valence of the ions. 

Many of the values required for calculating a GHK current are simulation variables, such as concentrations and voltage, simulation constants such as 
temperature, or fixed constants such as the Faraday constant and the gas constant. Such values are either known or can be found by STEPS during runtime and so are not part of 
object construction, with the exception of single-channel permeability which we will come to later. Like we saw in :doc:`Interface_Tutorial_5_Efield`, we then create the `GHKCurr` object. There are also optional keyword arguments ('virtual_oconc' and 'computeflux') and we'll 
see that which of these optional arguments are used depends on whether the mesh has an extracellular 'outer' compartment available (e.g. the 10um, 20um, 40um and 
80um meshes) or not (e.g. the 160um mesh)::

    with ssys:

        # ...

        if cyl160:
            OC_CaP = GHKCurr.Create(CaPchan[CaP_m3], Ca, CaP_P, virtual_oconc=Ca_oconc, computeflux=True)
        else:
            OC_CaP = GHKCurr.Create(CaPchan[CaP_m3], Ca, CaP_P, computeflux=True)

First let's look at the 'virtual_oconc' argument. This option allows us to not explicitly model the extracellular ('outer') concentration of the ion, useful because
often the extracellular compartment is not modelled. This option, rather, allows a fixed 'outer' concentration for the ion to be 
specified and that number will be used in the GHK flux calculations. The value of the parameter ``Ca_oconc`` in the extra.constants module is 2mM, 
so when the 160um mesh is used (when the ``cyl160`` flag is True) where there is no extracellular compartment, the extracellular concentration of Ca2+ in 
all GHK flux calculations will be 2mM. 

The second optional argument is 'computeflux'. This flag (which defaults to True) tells STEPS whether to model this GHK current process as ion transport 
or not. If 'computeflux' is True, then the calculated GHK current will result in transport of ions between the 'outer' and 'inner' compartments. 
For example, if over some 0.01ms time step, somewhere on the membrane a mean current of approximately 1.6pA is calculated through a membrane channel to which a GHK current is applied, 
then for an ion of valence 2+ this means that 50 ions moved from one compartment to the other. The direction of movement depends on the signs of the current
and the ion valence. The movement only occurs between surface tetrahedrons surrounding the membrane triangles in which the channels reside and so, for ions 
where this kind of process occurs, for accuracy it is necessary to model diffusion of these ions at least within the inner compartment 
and often within both compartments. This can be an expensive computation, particularly where concentrations are in the millimolar range, which shows the value of the 'computeflux'
flag- if the GHK flux is applied to an ion which does not have any other particularly important effects in the model other than its effect on membrane 
excitability (a possible example is potassium) then it may be a good labour-saver to clamp 'inner' and 'outer' concentrations of the ion and turn off the transport 
of ions as an approximation. However, in this model if we set 'computeflux' to False then the result would be no intracellular calcium, which is 
obviously not desirable, and so the 'computeflux' flag is set to True, as it usually will be for most ions in most models.  

For calcium (and only for calcium) we used `valence=2` to specify a valence of 2. 'Valence' can be an ambiguous term, but 
here it means the net elementary electrical charge per ion, which in this example for Ca2+ is +2. Negative valences can of course be specified by 
using a negative number. It is essential that this function is called to set a valence for any ion that will be used for a GHK current in the simulation- 
if no valence is specified the result will be an error. 

The third argument corresponds to single-channel permeability. Because conductance is not constant for a GHK current (apart from under certain unusual 
conditions) one value for a conductance parameter does not suffice. However, since single-channel permeability is often rather a difficult parameter
to define, STEPS does provide functionality for estimating the permeability. So we have two options for setting single-channel permeability: 
giving it directly to the GHKCurr constructor; or giving the result of a call to :func:`steps.API_2.model.GHKcurr.PInfo` to estimate the permeability from data. The first is straightforward and simply means providing single-channel 
permeability in S.I. units of cubic metres / second. In this model the parameter can be found in the :ref:`constants` and takes the value 
2.5e-20 cubic metres / second [#f4]_. If we did not provide the permeability during construction, we could also write::

    OC_CaP.P = CaP_P

The second option, the :func:`steps.API_2.model.GHKcurr.PInfo` function, requires some explanation. In effect, the conductance of a channel that is modelled 
by the GHK flux equation varies with 
voltage (:ref:`Figure 10.1 <figure_10_1>`) with a dependence on the 'outer' and 'inner' concentrations of the ion (in fact conductance is only constant with voltage 
when these concentrations are equal), as well as weakly on temperature. 

.. _figure_10_1_API2:

.. figure:: images/GHK_K.png
   :height: 6.0in
   :width: 8.0in

   `Figure 10.1: A single-channel GHK flux in the physiological range for a typical monovalent cation compared to an Ohmic approximation. The GHK flux is calculated with single-channel permeability of 9e-20 cubic metres / second, fixed extracellular concentration of 4mM, fixed intracellular concentration of 155mM and temperature of 20 Celsius. The single-channel Ohmic conductance is 20pS with reversal potential -77mV.`  


STEPS is able to estimate single-channel permeability from single-channel conductance, but for STEPS to do so the user must supply 
information about the conditions under which the conductance was measured, and in theory this should be enough to find the single-channel permeability since it is 
assumed constant (although there are occasions when permeability too can have some weak voltage dependence [#f3]_, 
which is, however, currently not possible to model with STEPS). Specifically, the :func:`steps.API_2.model.GHKcurr.PInfo` function requires arguments of:
estimated single-channel conductance [#f5]_ (units: Siemens), one voltage within the range at which conductance was measured (Volts), temperature (Kelvin), 'outer' concentration 
of the ion (molar), and 'inner' concentration of the ion (molar). Since the valence of the ion is known it is not necessary to supply that information to 
the :func:`steps.API_2.model.GHKcurr.PInfo` function. So, for example, for some GHKcurrent object called ``K_GHK``, if we measured single-channel conductance 
as 20pS in a small voltage range around -22mV at 20 degrees Celsius (293.15 Kelvin) with an estimated extracellular ion concentration of 4mM and 
intracellular concentration of 155mM, then we could create the current with::

    K_Pinfo = GHKCurr.PInfo(g = 20e-12, V = -22e-3, T = 293.15, oconc = 4e-3, iconc = 155e-3)
    K_GHK = GHKCurr.Create(Kchan[Ko], K, K_Pinfo)

and the single-channel permeability would be set to approximately 9e-20 cubic metres / second. The behaviour of such a channel is shown in :ref:`Figure 10.1 <figure_10_1>`.

We are now familiar, through aspects discussed so far in this chapter and other chapters, with most of the concepts applied for this model, so 
a very detailed description is not necessary for most remaining parts of the model. We move on to our other three ion channels in the model.

T-type Calcium channel
----------------------

Like the P-type Calcium channel, transitions between channel states of the T-type Calcium channel are voltage-dependent and we model the calcium current as a GHK current::

    ######## CaT channel ##########  
    CaT_m0h0, CaT_m0h1, CaT_m1h0, CaT_m1h1, CaT_m2h0, CaT_m2h1 = SubUnitState.Create()
    CaTchan = Channel.Create([CaT_m0h0, CaT_m0h1, CaT_m1h0, CaT_m1h1, CaT_m2h0, CaT_m2h1])

    # ...

    with ssys:

        # ...

        with CaTchan[...]:
            CaT_m0h0.s <r['CaTm0h0_m1h0']> CaT_m1h0.s <r['CaTm1h0_m2h0']> CaT_m2h0.s <r['CaTm2h0_m2h1']> CaT_m2h1.s
            r['CaTm0h0_m1h0'].K = 2 * VDepRate(alpham_cat), 1 * VDepRate(betam_cat)
            r['CaTm1h0_m2h0'].K = 1 * VDepRate(alpham_cat), 2 * VDepRate(betam_cat)
            r['CaTm2h0_m2h1'].K = 1 * VDepRate(alphah_cat), 1 * VDepRate(betah_cat)

            CaT_m1h0.s <r['CaTm1h0_m1h1']> CaT_m1h1.s
            r['CaTm1h0_m1h1'].K = VDepRate(alphah_cat), VDepRate(betah_cat)
            
            CaT_m0h0.s <r['CaTm0h0_m0h1']> CaT_m0h1.s <r['CaTm0h1_m1h1']> CaT_m1h1.s <r['CaTm1h1_m2h1']> CaT_m2h1.s
            r['CaTm0h0_m0h1'].K = 1 * VDepRate(alphah_cat), 1 * VDepRate(betah_cat)
            r['CaTm1h1_m2h1'].K = 1 * VDepRate(alpham_cat), 2 * VDepRate(betam_cat)
            r['CaTm0h1_m1h1'].K = 2 * VDepRate(alpham_cat), 1 * VDepRate(betam_cat)

        # ...

        if cyl160:
            OC_CaT = GHKCurr.Create(CaTchan[CaT_m2h1], Ca, CaT_P, virtual_oconc=Ca_oconc, computeflux=True)
        else:
            OC_CaT = GHKCurr.Create(CaTchan[CaT_m2h1], Ca, CaT_P, computeflux=True)


BK-type Calcium-activated Potassium channel
-------------------------------------------

The BK channel in the model undergoes both voltage-dependent and non-voltage dependent processes. This is an example of Channel States interacting with Species through surface reactions, 
 Here we will notice that Channel States (e.g. ``BK_C0``) appear alongside Species (``Ca``) in reactions:: 

    ##### BK channel ####################
    BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O0, BK_O1, BK_O2, BK_O3, BK_O4 = SubUnitState.Create()
    BKchan = Channel.Create([BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O0, BK_O1, BK_O2, BK_O3, BK_O4])

    # ...

    with ssys:

        # ...

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

        # ...

        OC_BK = OhmicCurr.Create(BKchan[BK_O0|BK_O1|BK_O2|BK_O3|BK_O4], BK_G, BK_rev)


:class:`steps.API_2.model.OhmicCurr` objects are applied to 5 
different channel states, demonstrating the support for multiple conducting/permeable states for a channel.


SK-type Calcium-activated Potassium channel
-------------------------------------------

The SK channel does not have any voltage dependence, and contains two conducting states::

    ###### SK channel ##################
    SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2 = SubUnitState.Create()
    SKchan = Channel.Create([SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2])

    # ...

    with ssys:

        # ...

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

        # ...

        OC_SK = OhmicCurr.Create(SKchan[SK_O1|SK_O2], SK_G, SK_rev)


Leak channel
------------

Although another option for setting the leak would have been to (later) set `sim.membrane.Res`, the leak conductance is described as a leak channel,  ::

    ###### Leak current channel #####
    Leak = SubUnitState.Create()
    L = Channel.Create([Leak])

    # ...

    with ssys:

        # ...

        OC_L = OhmicCurr.Create(L[Leak], L_G, L_rev)


Geometry specification
======================

This model is set up for relatively simple geometry- cylinders of diameter 2um and varying lengths from 10um to 160um. As discussed previously, 
the 160um cylinder does not include an extracellular compartment within the mesh whereas the 10um, 20um, 40um and 80um cylinders do, so the initialisation is slightly different for the 160um mesh compared 
to the others. 

First we separate the tetrahedrons into the 'inner' tetrahedrons, which will form the cytosolic compartment, and the 'outer' tetrahedrons, which will form the 
extracellular compartment. We do that by finding the tetrahedrons within a 1um radius along the cylinder axis (which is along the z axis) to form the inner 
compartment, and exclude those tetrahedrons from the complete set to find the outer compartment::

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

        # Record voltage from the central tetrahedron
        cent_tet = mesh.tets[0.0, 0.0, 0.0]

        ########## Create an intracellular compartment i.e. cytosolic compartment
        cyto = Compartment.Create(inner_tets, vsys)

        if not cyl160:
            outer = Compartment.Create(outer_tets)


We also find the tetrahedron at the centre of the mesh `cent_tet`, where we will record voltage from during simulation [#f6]_.
Now we move on to find the triangles that form the cell membrane between the intracellular and extracellular compartments. The way to do that is different
depending on whether we are using the 160um cylinder or not. The circular faces at each end of the cylinder are excluded from the membrane::

    with mesh:

        # ...

        if cyl160:
            # Ensure that we use points a small distance inside the boundary:
            minz, maxz = mesh.bbox.min.z, mesh.bbox.max.z
            memb_tris = TriList(tri for tri in mesh_stock.surface if minz < tri.center.z < maxz)
        else:
            print('Finding connecting triangles...')
            memb_tris = inner_tets.surface & outer_tets.surface


We also find the submembrane tetrahedrons, that is all tetrahedrons connected to a membrane triangle from the intracellular side::

    with mesh:

        # ...

        ########## Find the submembrane tets
        submemb_tets = TetList()
        for tri in memb_tris:
            submemb_tets += tri.tetNeighbs
        submemb_tets = submemb_tets & inner_tets

        print(len(submemb_tets))
        print('Volume of submembrane region is', submemb_tets.Vol)

    ########## Find the submembrane tets


And we are ready to create our membrane. That is a :class:`steps.API_2.geom.Membrane` object for which we are able to model electrical excitability by adding 
ion channels and solving potential across the membrane and within the intracellular conduction volume. For details of the method see [#f3]_.  
First we need to create a patch, named a little confusingly ``memb``::

    ########## Create a membrane
    if cyl160: 
        memb = Patch.Create(memb_tris, cyto, None, ssys)
    else:
        memb = Patch.Create(memb_tris, cyto, outer, ssys)


And then we create the membrane. Here we take advantage of previously found and stored connectivity optimisation (by function :func:`steps.API_1.solver.Tetexact.saveMembOpt`), 
in files such as /meshes/Cylinder2_dia2um_L10um_outer0_3um_0.3shell_0.3size_19156tets_adaptive.inp_optimalidx. Connectivity optimisation is discussed in [#f3]_ and :doc:`/memb_pot`::

    # For EField calculation
    membrane = Membrane.Create([memb], opt_file_name = './meshes/'+meshfile_ab+"_optimalidx")


Simulation with Tetexact
========================

Initialization
--------------

Before simulation we create the random number generator, but this time initialise with a fixed number. The reason for doing this is to ensure that initial conditions 
(placement of ion channels etc) is the same for each simulation iteration, so that the stochastic effects observed are purely from stochastic kinetics, and not due to 
different arrangements of channels::

    rng = RNG('mt19937', 512, seed=int(time.time()%1000))


We create the spatial stochastic solver object, turning the voltage calculation on by setting the 'calcMembPot' flag to True. Recall that we need one (and only one) 
:class:`steps.API_2.geom.Membrane` object to exist in the geometry for this to work::

    sim = Simulation('Tetexact', mdl, mesh, rng, calcMembPot=True)

We then create relevant folders and declare the data to be saved, as seen in previous chapters::

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


We can now start a new run::

    sim.newRun()


Next we see a new simulation property, :func:`steps.API_2.sim.Simulation.Temp`, which holds the simulation temperature. Currently, this will only influence any GHK flux rates, and 
will have no influence on any other kinetics. The value for ``TEMPERATURE`` is 34 degrees Celsius and we need to set temperature in Kelvin, following the usual S.I. rule 
in STEPS::

    sim.Temp = TEMPERATURE + 273.15


Next we inject ions, buffers and channels. Most values appear in the :ref:`constants`::

    if not cyl160: 
        sim.outer.Ca.Conc = Ca_oconc
        sim.outer.Ca.Clamped = True
        
    sim.cyto.Ca.Conc = Ca_iconc

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

    sim.memb.CaTchan[CaT_m0h0].Count = round(CaT_ro*surfarea*CaT_m0h0_p)
    sim.memb.CaTchan[CaT_m1h0].Count = round(CaT_ro*surfarea*CaT_m1h0_p)
    sim.memb.CaTchan[CaT_m2h0].Count = round(CaT_ro*surfarea*CaT_m2h0_p)
    sim.memb.CaTchan[CaT_m0h1].Count = round(CaT_ro*surfarea*CaT_m0h1_p)
    sim.memb.CaTchan[CaT_m1h1].Count = round(CaT_ro*surfarea*CaT_m1h1_p)
    sim.memb.CaTchan[CaT_m2h1].Count = round(CaT_ro*surfarea*CaT_m2h1_p)

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

    sim.memb.SKchan[SK_C1].Count = round(SK_ro*surfarea*SK_C1_p)
    sim.memb.SKchan[SK_C2].Count = round(SK_ro*surfarea*SK_C2_p)
    sim.memb.SKchan[SK_C3].Count = round(SK_ro*surfarea*SK_C3_p)
    sim.memb.SKchan[SK_C4].Count = round(SK_ro*surfarea*SK_C4_p)

    sim.memb.SKchan[SK_O1].Count = round(SK_ro*surfarea*SK_O1_p)
    sim.memb.SKchan[SK_O2].Count = round(SK_ro*surfarea*SK_O2_p)

    sim.memb.L[Leak].Count = int(L_ro * surfarea)


 
And finally set some parameters for the 'E-Field' voltage calculation using solver functions, in order: `sim.EfieldDT` to set the communication 
time-step between the voltage calculation and the reaction-diffusion 
simulation (from :ref:`constants`: ``EF_DT`` = 0.02e-3 seconds), `sim.membrane.Potential` to set the initial membrane potential (``init_pot`` = -60e-3 V), `sim.membrane.VolRes` to set resistivity of the conduction volume enclosed by the 
membrane (``Ra`` = 2.357 ohm.m) and `sim.membrane.Capac` to set membrane capacitance (``memb_capac`` = 1.5e-2 F/m\ :sup:`2`\)::

    sim.EfieldDT = EF_DT
    sim.membrane.Potential = init_pot
    sim.membrane.VolRes = Ra
    sim.membrane.Capac = memb_capac


Running the simulation
----------------------

At last we are ready to run the simulation, which is achieved simply by calls to function :func:`steps.API_2.sim.Simulation.run` within a simulation loop.::     

    for l in range(NTIMEPOINTS):
        print("Tpnt: ", l)
        sim.run(TIMECONVERTER*l)


:ref:`Figure 10.2 <figure_10_2>` shows five example runs of this script using the 80um mesh with different random number seeds. This model is analysed in much more depth in [#f1]_.

.. _figure_10_2_API2:

.. figure:: images/stochstoch80V.png
   :height: 6.0in
   :width: 8.0in

   `Figure 10.2: Voltage traces from five iterations on the 80um mesh.` 



.. rubric:: Footnotes
.. [#f1] Anwar H, Hepburn I, Nedelescu H, Chen W, De Schutter E (2013) Stochastic Calcium Mechanisms Cause Dendritic Calcium Spike Variability. The Journal of Neuroscience, 33(40): 15848-15867, doi: 10.1523/​JNEUROSCI.1722-13.2013. 
.. [#f2] Anwar H, Hong S, De Schutter E (2012) Controlling Ca2+-activated K+ channels with models of Ca2+ buffering in Purkinje cells. Cerebellum, 11(3):681-93, doi: 10.1007/s12311-010-0224-3. 
.. [#f3] Hepburn I, Cannon R and De Schutter E (2013) Efficient calculation of the quasi-static electrical potential on a tetrahedral mesh and its implementation in STEPS. Frontiers in Computational Neuroscience: 7:129, doi: 10.3389/fncom.2013.00129
.. [#f4] The same considerations for converting membrane permeability to single-channel permeability apply as for conductance discussed in :doc:`/memb_pot`, requiring some estimate of the channel density.
.. [#f5] Since it is assumed that conductance is measured by estimating the slope of an I-V curve over some small voltage range, the conductance will be treated as a slope conductance for the purposes of single-channel permeability estimation.
.. [#f6] We may record voltage from anywhere on the membrane surface or within the 'conduction volume' (here and in most models the conduction volume is the cytosolic compartment). 


